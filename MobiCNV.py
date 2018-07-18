##!/usr/local/bin/ python3.5
# coding: utf-8

# ==============================================================================

### IMPORT LIBRARIES
import sys		# system command
import csv		# read, write csv
import re		# regex
import argparse	# for options
import os		# for options
import pprint   # print data structure
import numpy as np
import math
import xlsxwriter
import vcf
#from operator import itemgetter

# ==============================================================================





#funcion to build differents dicts from coverage files
def build_dict(sample, line, region_number, prm, psm, key, chr_semaph):
	if(sample not in psm):
		psm[sample] = {"rawDocSum": float(line[4])}
	else:
		psm[sample]["rawDocSum"] += float(line[4])
	#deal with chr1 or 1
	chrom = line[0]
	chr_reg = r'^chr'
	if not re.match(chr_reg, chrom):
		chr_semaph = True
	if chr_semaph:
		chrom = "chr" + chrom
	#create dictionnary with rawDoc values per sample per region
	coordinate = (key,chrom,line[1],line[2],line[3])
	if coordinate not in prm :
		region_number += 1
		prm[coordinate] = {sample : {"rawDoc" : float(line[4])}}
	else :
		prm[coordinate][sample] = {"rawDoc" : float(line[4])}	
	return(region_number, prm, psm, chr_semaph)

#############

#function to populate per_region_metrics dictionnary with mean coverage per exons
#two metrics:
#1- full mean including all samples per region just to be printed in final file for informative purpose
#2- mean excluding current sample for normalisation purpose
def exon_mean(prm):
	for coordinate in prm :
		for sample_name in prm[coordinate]:
			total = full_total = sample_number = 0
			for sample_other in prm[coordinate]:
				if sample_other != sample_name:
					total += prm[coordinate][sample_other]['rawDoc']
					sample_number +=1
				full_total += prm[coordinate][sample_other]['rawDoc']
			#regionMeanOtherSamples = total / sample_number
			prm[coordinate][sample_name]["regionMeanDoc"] = int(full_total / (sample_number+1))
			prm[coordinate][sample_name]["regionMeanOtherSamples"] = round(float(total / sample_number), 3)
	return(prm)

#############


def compute_ratio(psm, prm, region_number, VcfDir, variants):
	#loop to compute per region the mean coverage of all regions except the ROI
	for sample_name in psm:
		psm[sample_name]["regionMeanOtherSamplesSum"] = 0
		for coordinate in prm :
			#compute mean sum per exon
			psm[sample_name]["regionMeanOtherSamplesSum"] += prm[coordinate][sample_name]["regionMeanOtherSamples"]
		#mean for all samples except current sample
		psm[sample_name]["totalMeanOtherSample"] = round(psm[sample_name]["regionMeanOtherSamplesSum"] / region_number, 3)
	#loop to compute normalised values per region (vakue for the ROI and for all others excluding ROI)
	for sample_name in psm:	
		for coordinate in prm :
			#normalisation per exon
			prm[coordinate][sample_name]["normalisedMeanOtherSamples"] = round(prm[coordinate][sample_name]["regionMeanOtherSamples"] / psm[sample_name]["totalMeanOtherSample"], 3)
			prm[coordinate][sample_name]["normalisedRegion"] = round(float(prm[coordinate][sample_name]['rawDoc'] / psm[sample_name]["meanRawDoc"]), 3)
	#computes final ratio
	for coordinate in prm :
		for sample_name in prm[coordinate]:			
			try :
				normalised_ratio = float(prm[coordinate][sample_name]["normalisedRegion"]) / float(prm[coordinate][sample_name]["normalisedMeanOtherSamples"])
			except ZeroDivisionError :
				normalised_ratio = float(0) 
			prm[coordinate][sample_name]["normalisedRatio"] = round(float(normalised_ratio), 3)
	#computes stdev for final ratio
	for coordinate in prm :
	  	for sample_name in prm[coordinate]:
	  		std_dev=[]
	  		for sample_name2 in prm[coordinate]:
	  			if(sample_name != sample_name2):
	  				std_dev.append(prm[coordinate][sample_name2]["normalisedRatio"])
	  		prm[coordinate][sample_name]["ratioStdev"]=round(float(np.std(std_dev)), 3)
	#interpretation
	#<0.3 => hom del
	#>1.7 => hom dup
	#betwwen 0.8 and 1.2 => normal
	#between 0.3 and 0.8 =>supect het del then
	#	between 1-2.5sigma and 1 => normal
	#	<1-2.5sigma => het del
	#between 1.2 and 1.7 =>supect het del then
	#	between 1 and 1+2.5sigma => normal
	#	>1+2.5sigma => het dup
	het_high = 1.3
	het_low = 0.7
	hom_high = 1.7
	hom_low = 0.3
	xfactor = 2
	#dev purpose
	r=0
	s=0
	for coordinate in prm :
		#print(coordinate)
		for sample_name in prm[coordinate]:
			ratio = prm[coordinate][sample_name]["normalisedRatio"]
			if ratio < hom_low:
				prm[coordinate][sample_name]["MobiAdvice"] = "HomDel"
			elif ratio > hom_high:
				prm[coordinate][sample_name]["MobiAdvice"] = "HomDup"
			elif ratio > het_low and ratio < het_high:
				prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
			elif ratio >= hom_low and ratio <= het_low:
				dynamic_threshold = 1 - (xfactor * prm[coordinate][sample_name]["ratioStdev"])
				if ratio < dynamic_threshold:
					vcf_semaph = 0
					#check against VCF info (het variants) if relevant
					if VcfDir != False:
						for var_pos in variants:
							#chr must match
							if var_pos[0] == coordinate[1]:
								#if var_pos[0] == 'chrX' and "gender" in psm[sample_name] and psm[sample_name]["gender"] == 'female':
								#	break
								if int(var_pos[1]) >= int(coordinate[2]) and int(var_pos[1]) <= int(coordinate[3]):
									prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
									vcf_semaph = 1
									r+=1
									#print(var_pos[0], coordinate[1], var_pos[1], coordinate[2], coordinate[3])
									break
					if vcf_semaph == 0:
						prm[coordinate][sample_name]["MobiAdvice"] = "HetDel"
						s+=1
				else:
					prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
			elif ratio >= het_high and ratio <= hom_high:
				dynamic_threshold = 1 + (xfactor * prm[coordinate][sample_name]["ratioStdev"])
				if ratio > dynamic_threshold:
					prm[coordinate][sample_name]["MobiAdvice"] = "HetDup"
				else:
					prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
	print("HetDels: " + str(s) + " vcfed: " + str(r))
	return (psm, prm)

#############

#function to basically format new Excel sheets
def format_sheet(sheet, last_col, styles):
	sheet.freeze_panes(1, 6)
	sheet.set_row(0, 20, styles[4])#style5
	sheet.set_column('A:F', 15, styles[4])
	sheet.set_column('E:E', 25)
	sheet.set_column('F:F', 15, styles[5])
	sheet.set_column(last_col+2, last_col+2, 15)
	sheet.set_column(last_col+3, last_col+3, 15)
	sheet.write(3, last_col+2, 'Legend:', styles[4])
	sheet.write(4, last_col+2, '', styles[0])
	sheet.write(5, last_col+2, '', styles[1])
	sheet.write(6, last_col+2, '', styles[2])
	sheet.write(7, last_col+2, '', styles[3])
	sheet.write(4, last_col+3, 'Homozygous deletion', styles[4])
	sheet.write(5, last_col+3, 'Heterozygous deletion', styles[4])
	sheet.write(6, last_col+3, 'Heterozygous duplication', styles[4])
	sheet.write(7, last_col+3, 'Homozygous duplication', styles[4])

#############

#function to add conditional formatting for mean global DoC of the regions (3_traffic_lights)
def add_conditionnal_format(worksheet, threshold, start, end):
	#add a conditionnal format to mean DoC
	if (start == 0):
		start = 2
	cell_range = "F" + str(start) + ":F" + str(end)
	#http://xlsxwriter.readthedocs.io/working_with_conditional_formats.html
	worksheet.conditional_format(
		cell_range,
		{'type': 'icon_set',
		 'icon_style': '3_traffic_lights',
		 'icons': [{'criteria': '>=', 'type': 'number', 'value': 100},
		{'criteria': '<=',  'type': 'number', 'value': threshold}]}
	)
	#xlsxwriter < 1.0.0 => cell formatting
# 	worksheet.conditional_format(cell_range, {'type': 'cell',
#                                           'criteria': '>=',
#                                           'value': threshold,
#                                           'format': format2})
# 	worksheet.conditional_format(cell_range, {'type': 'cell',
#                                           'criteria': '<',
#                                           'value': threshold,
#                                           'format': format1})
#############

#function to actually create result worksheets
def print_worksheet(name, last_col, last_col_2_hide, workbook, prm, quality, reduced_regions, panel_regions, low_cov_regions, Panel, panel_list, psmX, styles, number_of_file):
	#sheet creation
	worksheet = workbook.add_worksheet(str(name))
	format_sheet(worksheet, last_col, styles)
	if quality == "summary":
		worksheet.activate()
	#i=first row
	i=j=0
	#dict iterations
	for region in sorted(prm):
		#j=col
		#j=0
		if i == 0:
			#sheet headers
			headers = ("Region Number", "Chromosome", "Start", "End", "Annotation", "Mean DoC")
			for text in headers:
				worksheet.write(i, j, text, styles[4])
				j+=1
			for sample in prm[region]:
				worksheet.write(i, j, sample + "_rawDoc", styles[4])
				worksheet.write(i, j+number_of_file, sample + "_regionMeanOtherSamples", styles[4])
				worksheet.write(i, j+(2*number_of_file), sample + "_normalisedRegion", styles[4])
				worksheet.write(i, j+(3*number_of_file), sample + "_normalisedMeanOtherSamples", styles[4])
				worksheet.write(i, j+(4*number_of_file), sample + "_ratioStdev", styles[4])
				worksheet.write(i, j+(5*number_of_file), sample + "_normalisedRatio", styles[4])
				j+=1
			i+=1
		#else:
		j=0
		for row_header in region:
			worksheet.write(i, j, row_header)
			if j == 4 and Panel:
				for gene in panel_list:
					if re.compile(r'.*' + gene + '.*').search(row_header) :
						panel_regions[region] = prm[region]
			j+=1
		m=10
		for sample in prm[region]:
			j+=1
			if i == 1 and sample in psmX:
				style = styles[7]
				if psmX[sample]["gender"] == 'female':
					style = styles[6]
				worksheet.write(m, last_col+2, sample, style)
				worksheet.write(m, last_col+3, psmX[sample]["gender"], style)
				worksheet.write(m, last_col+4, round(psmX[sample]["xratio"], 2), style)
				m+=1	
			worksheet.write(i, j, prm[region][sample]["rawDoc"])
			worksheet.write(i, j+number_of_file, prm[region][sample]["regionMeanOtherSamples"])
			worksheet.write(i, j+(2*number_of_file), prm[region][sample]["normalisedRegion"])
			worksheet.write(i, j+(3*number_of_file), prm[region][sample]["normalisedMeanOtherSamples"])
			worksheet.write(i, j+(4*number_of_file), prm[region][sample]["ratioStdev"])
			#define cell style - default style5: just bold
			cell_style = styles[4]
			if quality == "global" and prm[region][sample]["MobiAdvice"] != 'Normal':
				if region not in reduced_regions:
					#here we optionally check whether the sample's vcf contains heterozygous (PASS) variants in the region of interest
					#here we check if number_of_sample > 4 and meanDoC all samples == 0 => we don't put that region in the summary sheet
					if number_of_file <= 3:
						reduced_regions[region] = prm[region]
					elif prm[region][sample]["regionMeanDoc"] > 100:
						reduced_regions[region] = prm[region]
					elif prm[region][sample]["regionMeanDoc"] > 0:
						low_cov_regions[region] = prm[region]
			if prm[region][sample]["MobiAdvice"] == "HomDel":
				cell_style = styles[0]
			elif prm[region][sample]["MobiAdvice"] == "HetDel":
				cell_style = styles[1]
			elif prm[region][sample]["MobiAdvice"] == "HetDup":
				cell_style = styles[2]
			elif prm[region][sample]["MobiAdvice"] == "HomDup":
				cell_style = styles[3]
			worksheet.write(i, j+(5*number_of_file), prm[region][sample]["normalisedRatio"], cell_style)
			last_sample = sample
		#mean global doc for region
		worksheet.write(i, 5, prm[region][last_sample]["regionMeanDoc"])
		i+=1

	worksheet.set_column(6, last_col_2_hide, None, None, {'level': 1, 'hidden': True})
	add_conditionnal_format(worksheet, 50, 2, len(list(prm))+1)
	worksheet.protect()

	if quality == "global":
		return reduced_regions, low_cov_regions, panel_regions
	
#############

def main():
	
	#for developping purpose
	#pp = pprint.PrettyPrinter(indent=4, depth=6)
	
	#added david 27/03/2018 - deal args
	parser = argparse.ArgumentParser(usage='python MobiCNV.py [-i PATH/TO/DoC_FILE_DIR/ -t (tsv|csv) -p PATH/TO/GENE/TXT/FILE -o PATH/TO/OUTPUT/FILE -v PATH/TO/VCFs_DIR/]', description='Compute CNV predictions from Illumina Depth of Coverage files (generated by MSR or LRM in DNA enrichment mode. Can also deal with samtools bedcov slightly modified files. Generates an Excel File with all computed values in Autosomes, X and Y sheets if relevant, in addition to a summary sheet which displays aberrant regions +- 1. Optionally, an additional panel sheet focuses on a list of genes passed as argument as a txt file (one gene per line).')
	parser.add_argument('-i', '--input', default='.', help='Path to the directory containing the coverage files. The sample name will be deduced from the file names which should be sample_coverage.[ct]sv.')
	parser.add_argument('-p', '--panel', default='', help='If a text file with gene names is provided, and if the coverage files include gene annotation, a supplementary sheet will be created focusing on the genes of interest.')
	parser.add_argument('-t', '--type', default='csv', help='Can be csv or tsv.')
	parser.add_argument('-o', '--output', default='mobicnv.xlsx')
	parser.add_argument('-v', '--vcf', default='', help='Path to the directory containing optional VCF files. Can be the same as coverage directory, but must be specified for the VCF search to be activated. Data can be centralised in a single VCF, or in mutliple VCFs. MobiCNV will search for samples in any present VCF based on the name of the coverage files.')
	args = parser.parse_args()
	
	#Variables declaration
	#print(args.path)
	#sys.exit()
	Path = args.input
	Panel = args.panel
	Type = args.type
	OutFile = args.output
	VcfDir = args.vcf
	ext = "csv"
	delim = ","
	if Type == "tsv":
		ext = "tsv"
		delim = "\t"
	if Panel == '':
		Panel = False
	if VcfDir == '':
		VcfDir = False
	else:
		vcf_list = os.listdir(VcfDir)
	file_list = os.listdir(Path)
	number_of_file = 0
	region_number = 0
	per_region_metrics = {}
	per_sample_metrics = {}
	total_mean = 0
	per_region_metrics_ChrX = {}
	per_sample_metrics_ChrX = {}
	total_mean_ChrX = 0
	region_number_ChrX = 0
	per_region_metrics_ChrY = {}
	per_sample_metrics_ChrY = {}
	total_mean_ChrY = 0
	region_number_ChrY = 0
	variants = {}
	chr_semaph = False	
	
	print("\nFiles that will be considered:\n")
	
	for i in file_list:
		#### à fonctionnariser
		regex = re.compile(r'^([^\.].+)[\._]coverage\.%s$'%ext)
		match_file = regex.search(os.path.basename(i))
		if match_file:
			#for each file
			number_of_file += 1
			sample = match_file.group(1)
			print("Sample: " + sample)
			print("Metrics file: " + i)
			key = 1
			#we read and each line is converted into a dictionnary
			#2 dicts, one per sample for global metrics and one more detailed for per region metrics
			#Each dict type (region, sample) is in fact composed of up to 3 dicts, one for autosomes, one for chrX, and one for chrY
			#coz more simple to handle when writing in separate tabs in excel
			with open(Path + i, 'r') as csvfile:
				csvreader = csv.reader(csvfile, delimiter =delim)
				for line in csvreader:
					line = [w.replace(',', '.') for w in line]
					expression = r'^.*#.*$' # remove header
					#autosomes
					if not re.match(expression, line[0]) and line[0] != "chrX" and line[0] != 'X' and line[0] != "chrY" and line[0] != "Y":
						(region_number, per_region_metrics, per_sample_metrics, chr_semaph) = build_dict(sample, line, region_number, per_region_metrics, per_sample_metrics, key, chr_semaph)
						key+=1
					#chrX
					elif not re.match(expression, line[0]) and (line[0] == "chrX" or line[0] == "X"):
						(region_number_ChrX, per_region_metrics_ChrX, per_sample_metrics_ChrX, chr_semaph) = build_dict(sample, line, region_number_ChrX, per_region_metrics_ChrX, per_sample_metrics_ChrX, key, chr_semaph)
						key+=1
					#chrY
					elif not re.match(expression, line[0]) and (line[0] == "chrY" or line[0] == "Y"):
						(region_number_ChrY, per_region_metrics_ChrY, per_sample_metrics_ChrY, chr_semaph) = build_dict(sample, line, region_number_ChrY, per_region_metrics_ChrY, per_sample_metrics_ChrY, key, chr_semaph)
						key+=1
	
			#we build vcf dict if possible chr-pos => {sample, status (1)}
			if VcfDir != False:
				#vcf_regexp = re.compile(r''+sample+'\..*\.?vcf\.?g?z?$')
				chr_reg = r'^chr'
				vcf_regexp = re.compile(r'.+\.vcf\.?g?z?$')
				for vcf_file in vcf_list:
					found_sample = False
					match_vcf = vcf_regexp.search(os.path.basename(vcf_file))
					if match_vcf:
						#print("Associated VCF: " + vcf_file)
						vcf_reader = vcf.Reader(open(VcfDir + vcf_file, 'r'))
						test_record = next(vcf_reader)
						#if test_record.genotype(sample):
						for vcf_calls in test_record.samples:
							if vcf_calls.sample == sample:
						#if sample in test_record.samples.sample:
								#ok our sample is here  we can read the entire vcf
								print("Associated VCF: " + vcf_file)
								vcf_chr_semaph = False
								for record in vcf_reader:
									#we store only heterozygous calls or all cals on chrX to treat male dels on X chr
									sample_call = record.genotype(sample)
									#if sample_call.is_het == True or (record.CHROM == 'chrX' or record.CHROM == 'X'):
									if sample_call.is_het == True:
										chrom = record.CHROM
										if not re.match(chr_reg, chrom):
											vcf_chr_semaph = True
										if vcf_chr_semaph:
											chrom = "chr" + chrom
										position = (chrom, record.POS)
										variants[position] = {sample: 1}
										#print(record.CHROM + "-" + str(record.POS) + "-" + str(sample) + "-" + str(record.heterozygosity))
								found_sample = True
								break
						if found_sample:
							break
			print("")
	#############
	
	#############
	#pp.pprint(per_sample_metrics)
	#sys.exit()
	
	
	
	#############
	# we iterate sample dict to compute the global Doc
	for sample_name in per_sample_metrics:
			per_sample_metrics[sample_name]["meanRawDoc"] = per_sample_metrics[sample_name]["rawDocSum"]/ region_number
	# FOR X and in addition we compute on the fly the x ratio which determinates the gender
	#then the ratio and predicted gender are added to the per_sample_metrics_chrX dict ONLY
	# FOR Y we look for Y regions in female samples and print a warning if found unconsistant
	if region_number_ChrX > 0:
		for sample_name in per_sample_metrics_ChrX:
				per_sample_metrics_ChrX[sample_name]["meanRawDoc"] = round(per_sample_metrics_ChrX[sample_name]["rawDocSum"]/ region_number_ChrX, 3)
				per_sample_metrics_ChrX[sample_name]["xratio"] = round(float(per_sample_metrics_ChrX[sample_name]["meanRawDoc"] / per_sample_metrics[sample_name]["meanRawDoc"]), 3)
				per_sample_metrics_ChrX[sample_name]["gender"] = "male"
				if (per_sample_metrics_ChrX[sample_name]["xratio"] > 0.65):
					per_sample_metrics_ChrX[sample_name]["gender"] = "female"
	if region_number_ChrY > 0:
		for sample_name in per_sample_metrics_ChrY:
				per_sample_metrics_ChrY[sample_name]["meanRawDoc"] = per_sample_metrics_ChrY[sample_name]["rawDocSum"]/ region_number_ChrY
				if (per_sample_metrics_ChrY[sample_name]["meanRawDoc"] > 1 and per_sample_metrics_ChrX[sample_name]["gender"] != "male"):
					print("\nWARNING Gender inconsistancy for " + sample_name + " reads on Y chr with X ratio > 0.65\n")
					per_sample_metrics_ChrX[sample_name]["gender"] = "male/female"
	#############
	
	#############
	# we iterate sample dict to compute mean DoCs per exons and normalised Doc
	
	per_region_metrics = exon_mean(per_region_metrics)
	if region_number_ChrX > 0:
		per_region_metrics_ChrX = exon_mean(per_region_metrics_ChrX)
	if region_number_ChrY > 0:
		per_region_metrics_ChrY = exon_mean(per_region_metrics_ChrY)
	#############
	
	#############
	# we iterate sample dict to compute mean sum per sample
	# except for the pateint's current region - added to patient dict
	
	
	(per_sample_metrics, per_region_metrics) = compute_ratio(per_sample_metrics, per_region_metrics, region_number, VcfDir, variants)
	if region_number_ChrX > 0:
		(per_sample_metrics_ChrX, per_region_metrics_ChrX) = compute_ratio(per_sample_metrics_ChrX, per_region_metrics_ChrX, region_number_ChrX, VcfDir, variants)
	if region_number_ChrY > 0:
		(per_sample_metrics_ChrY, per_region_metrics_ChrY) = compute_ratio(per_sample_metrics_ChrY, per_region_metrics_ChrY, region_number_ChrY, VcfDir, variants)
	
	# pp.pprint(per_region_metrics)
	#sys.exit()
	#############
	
	#############
	#We build a small list of genes of interest
	panel_list = []
	if Panel:
		panel = open(Panel, 'r')
		# panelreader = csv.DictReader(panel, delimiter='\t')
		for gene in panel :
			panel_list.append(gene.rstrip())
	###########
	
	###########
	### Excel conversion
	
	
	workbook = xlsxwriter.Workbook(OutFile)
	#Define style of cells
	style1 = workbook.add_format({'bold': True, 'bg_color': '#FF3333', 'locked' : True})
	style2 = workbook.add_format({'bold': True, 'bg_color': '#FFC25E', 'locked' : True})
	style3 = workbook.add_format({'bold': True, 'bg_color': '#5EBBFF', 'locked' : True})
	style4 = workbook.add_format({'bold': True, 'bg_color': '#8F5EFF', 'locked' : True})
	style5 = workbook.add_format({'bold': True, 'locked' : True})
	style6 = workbook.add_format({'bold': True, 'num_format': 1, 'locked' : True})
	style7 = workbook.add_format({'bold': True, 'color': '#F791E7', 'locked' : True})
	style8 = workbook.add_format({'bold': True, 'color': '#5EBBFF', 'locked' : True})
	styles = (style1, style2, style3, style4, style5, style6, style7, style8)
	#http://xlsxwriter.readthedocs.io/example_conditional_format.html#ex-cond-format
	# Add a format. Light red fill with dark red text.
	format1 = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
	# Add a format. Green fill with dark green text.
	format2 = workbook.add_format({'bg_color': '#C6EFCE', 'font_color': '#006100'})
	
	last_col = number_of_file*6 + 5
	#sample number x number of data blocks + 6 first cols to show - 1 coz numbering begins at 0 (A = 0)
	last_col_2_hide = number_of_file*5 + 5
	#sample number x number of data blocks to hide + 6 first cols to show - 1 coz numbering begins at 0 (A = 0)
	
	
	
	print("\nBuilding Excel File:")
	print("Autosomes worksheet...")
	summary_regions = {}
	low_cov_summary_regions = {}
	panel_regions = {}
	
	(summary_regions, low_cov_summary_regions, panel_regions) = print_worksheet('Autosomes', last_col, last_col_2_hide, workbook, per_region_metrics, 'global', summary_regions, panel_regions, low_cov_summary_regions, Panel, panel_list, per_sample_metrics_ChrX, styles, number_of_file)
	
	if region_number_ChrX > 0:
		print("ChrX worksheet...")
		(summary_regions, low_cov_summary_regions, panel_regions) = print_worksheet('Chromosome X', last_col, last_col_2_hide, workbook, per_region_metrics_ChrX, 'global', summary_regions, panel_regions, low_cov_summary_regions, Panel, panel_list, per_sample_metrics_ChrX, styles, number_of_file)
	if region_number_ChrY > 0:
		print("ChrY worksheet...")
		(summary_regions, low_cov_summary_regions, panel_regions) = print_worksheet('Chromosome Y', last_col, last_col_2_hide, workbook, per_region_metrics_ChrY, 'global', summary_regions, panel_regions, low_cov_summary_regions, Panel, panel_list, per_sample_metrics_ChrX, styles, number_of_file)
	
	#pp.pprint(summary_regions)
	
	print("Summary worksheet...")
	print_worksheet('Summary', last_col, last_col_2_hide, workbook, summary_regions, 'summary', '', '', '', '', '', per_sample_metrics_ChrX, styles, number_of_file)
	if len(low_cov_summary_regions) > 0:
		print("Low Coverage Summary Worksheet...")
		print_worksheet('LowCovSummary', last_col, last_col_2_hide, workbook, low_cov_summary_regions, 'low_cov_summary', '', '', '', '', '', per_sample_metrics_ChrX, styles, number_of_file)
	if Panel:
		print("Panel worksheet...")
		print_worksheet('Panel', last_col, last_col_2_hide, workbook, panel_regions, 'panel', '', '', '', '', '', per_sample_metrics_ChrX, styles, number_of_file)
	#pp.pprint(per_sample_metrics_ChrX)
	
	workbook.close()
	
	
	sys.exit()

if __name__ == "__main__":
    main()

#
#def write_small_worksheets(selected, start, first_row, small_worksheet, col_list, last_col, regex_r, threshold_del_hmz, threshold_del_htz, threshold_dup_htz, threshold_dup_hmz):
#	#called inside writing_total to wirte summary and panel sheets
#	i = 0
#	uniq_selected = list(set(selected))
#	for column in col_list:
#		j = first_row
#		if(start > 0):
#			j = start
#		for item in range(len(column)):
#			if (item in uniq_selected):
#				if regex_r.search(column[0]):
#					if(item > 0):
#						if(float(column[item]) <= threshold_del_hmz):
#							small_worksheet.write(j, i , column[item],style1)
#						elif(float(column[item]) <= threshold_del_htz):
#							small_worksheet.write(j, i, column[item],style2)
#						elif(float(column[item]) <= threshold_dup_htz):
#							small_worksheet.write(j, i, column[item], style5)
#						elif(float(column[item]) <= threshold_dup_hmz):
#							small_worksheet.write(j, i, column[item],style3)
#						else :
#							small_worksheet.write(j, i, column[item],style4)
#					else:
#						small_worksheet.write(j, i, column[item], style5)
#				else:
#					try:
#						small_worksheet.write(j,i,int(column[item]))
#					except ValueError:
#						small_worksheet.write(j,i,column[item])
#				j+=1
#		i+=1
#	small_worksheet.write(9, last_col+2, "Sample ID", style5)
#	small_worksheet.write(9, last_col+3, "Predicted Gender", style5)
#	small_worksheet.write(9, last_col+4, "X ratio", style5)
#	if region_number_ChrX > 0 or region_number_ChrY > 0:
#		m = 10
#		for sample_name in dict_gender:
#			style = style8
#			if dict_gender[sample_name]["gender"] == 'female':
#				style = style7
#			small_worksheet.write(m, last_col+2, sample_name, style)
#			small_worksheet.write(m, last_col+3, dict_gender[sample_name]["gender"], style)
#			small_worksheet.write(m, last_col+4, round(dict_gender[sample_name]["xratio"], 2), style)
#			m += 1
#	return (uniq_selected, j)
#
#
#
#def writing_total(worksheet, txt_file, threshold_del_hmz, threshold_del_htz, threshold_dup_htz, threshold_dup_hmz, last_col_2_hide, last_col, start1=0, start2=0):
#	# TODO: change it to parameter and save worksheet
#	# if panel:
#
#	worksheet = workbook.add_worksheet(str(worksheet))
#	format_sheet(worksheet, last_col)
#	# worksheet.freeze_panes(1, 5)
#	# worksheet.set_row(0, 20, style5)
#	# worksheet.set_column('A:E', 15, style5)
#	# worksheet.set_column('D:D', 25)
#	#structure data from txt
#	f = open(str(txt_file), 'r+')
#	row_list = []
#	for row in f:
#		row_list.append(row.split('\t'))
#	# pp.pprint(row_list)
#
#	column_list = zip(*row_list)
#	column_list2 = zip(*row_list)
#	column_list3 = zip(*row_list)
#	i = 0
#	l = 0
#	interesting = []
#	gene4interest = []
#	regex_ratio = re.compile(r'(.*)_ratio$')
#	regex_noCNV = re.compile(r'no CNV')
#	regex_region = re.compile(r'RegionID')
#	for column in column_list:
#		for item in range(len(column)):
#			if regex_region.search(column[0]):
#				if (Panel != False):
#					for gene in liste_panel:
#						if re.compile(r'.*' + gene + '.*').search(column[item]) :
#							# print (column[item])
#							gene4interest.append(item)
#				else:
#					gene4interest.append(item)
#			if regex_ratio.search(column[0]):
#				if (item > 0) :
#					if(float(column[item]) <= threshold_del_hmz):
#						worksheet.write(item, i , column[item],style1)
#						for k in range(item-1,item+2):
#							if (k > 0):
#								interesting.append(k)
#					elif(float(column[item]) <= threshold_del_htz):
#						worksheet.write(item, i, column[item],style2)
#						for k in range(item-1,item+2):
#							if (k > 0):
#								interesting.append(k)
#					elif(float(column[item]) <= threshold_dup_htz):
#						worksheet.write(item, i, column[item], style5)
#					elif(float(column[item]) <= threshold_dup_hmz):
#						worksheet.write(item, i, column[item],style3)
#						for k in range(item-1,item+2):
#							if (k > 0):
#								interesting.append(k)
#					else :
#						worksheet.write(item, i, column[item],style4)
#						for k in range(item-1,item+2):
#							if (k > 0):
#								interesting.append(k)
#				else:
#					worksheet.write(item, i, column[item], style5)
#			else:
#				try:
#					worksheet.write(item,i,int(column[item]))
#				except ValueError:
#					worksheet.write(item,i,column[item])
#			if (item == 0):
#				summary.write(item,i,column[item], style5)
#				# if panel:
#				if (Panel != False):
#					worksheet2.write(item,i,column[item], style5)
#		i+=1
#	#i = 0
#	(uniq_interesting, j) = write_small_worksheets(interesting, start1, 0, summary, column_list2, last_col, regex_ratio, threshold_del_hmz, threshold_del_htz, threshold_dup_htz, threshold_dup_hmz)
#	if (Panel != False):
#		(uniq_interesting_panel, l) = write_small_worksheets(gene4interest, start2, 1, worksheet2, column_list3, last_col,  regex_ratio, threshold_del_hmz, threshold_del_htz, threshold_dup_htz, threshold_dup_hmz)
#		worksheet2.set_column(5,last_col_2_hide, None, None, {'level': 1, 'hidden': True})
#		add_conditionnal_format(worksheet2, 50, start2, start2 + len(gene4interest))
#	# uniq_interesting = list(set(interesting))
#	# for column in column_list2:
#	# 	j = 0
#	# 	if(start1 > 0):
#	# 		j = start1
#	# 	for item in range(len(column)):
#	# 		if (item in uniq_interesting):
#	# 			if regex_ratio.search(column[0]):
#	# 				if(item > 0):
#	# 					if(float(column[item]) <= threshold_del_hmz):
#	# 						summary.write(j, i , column[item],style1)
#	# 					elif(float(column[item]) <= threshold_del_htz):
#	# 						summary.write(j, i, column[item],style2)
#	# 					elif(float(column[item]) <= threshold_dup_htz):
#	# 						summary.write(j, i, column[item], style5)
#	# 					elif(float(column[item]) <= threshold_dup_hmz):
#	# 						summary.write(j, i, column[item],style3)
#	# 					else :
#	# 						summary.write(j, i, column[item],style4)
#	# 				else:
#	# 					summary.write(j, i, column[item], style5)
#	# 			else:
#	# 				summary.write(j,i,column[item])
#	# 			j+=1
#	# 	i+=1
#	
#
#	# i = 0
#	# # part to dev
#	# if (Panel != False):
#	# 	#l=0
#	# 	uniq_interesting_panel = list(set(gene4interest))
#	# 	for column in column_list3:
#	# 		l = 1
#	# 		if (start2 > 0 ):
#	# 			l = start2
#	# 		for item in range(len(column)):
#	# 			if (item in uniq_interesting_panel):
#	# 				if regex_ratio.search(column[0]):
#	# 					if(item > 0):
#	# 						if(float(column[item]) <= threshold_del_hmz):
#	# 							worksheet2.write(l, i , column[item],style1)
#	# 						elif(float(column[item]) <= threshold_del_htz):
#	# 							worksheet2.write(l, i, column[item],style2)
#	# 						elif(float(column[item]) <= threshold_dup_htz):
#	# 							worksheet2.write(l, i, column[item], style5)
#	# 						elif(float(column[item]) <= threshold_dup_hmz):
#	# 							worksheet2.write(l, i, column[item],style3)
#	# 						else :
#	# 							worksheet2.write(l, i, column[item],style4)
#	# 					else:
#	# 						worksheet2.write(l, i, column[item], style5)
#	# 				else:
#	# 					worksheet2.write(l,i,column[item])
#	# 				l+=1
#	# 		i+=1
#	# 	i = 0
#	# 	worksheet2.set_column('F:BM', None, None, {'level': 1, 'hidden': True})
#	# 	add_conditionnal_format(worksheet2, 50, start2, start2 + len(gene4interest))
#	worksheet.set_column('F:BM', None, None, {'level': 1, 'hidden': True})
#	worksheet.set_column(5, last_col_2_hide, None, None, {'level': 1, 'hidden': True})
#	add_conditionnal_format(worksheet, 50, 2, len(row_list))
#	worksheet.protect()
#	summary.set_column(5, last_col_2_hide, None, None, {'level': 1, 'hidden': True})
#	add_conditionnal_format(summary, 50, start1, start1 + len(uniq_interesting))
#	return (j,l)
#
##worksheet for autosomes
#(start1, start2) = writing_total('Autosomes','cnv_analysis_sorted.txt', 0.3, 0.7, 1.3, 1.7, last_col_2_hide, last_col)
## print (start1, start2)
##worksheet for ChrX
#start3 = 0
#start4 = 0
#if region_number_ChrX > 0:
#	(start3, start4) = writing_total('Chromosome_X','cnv_analysis_ChrX_sorted.txt', 0.3, 0.7, 1.3, 1.7, last_col_2_hide, last_col, start1, start2)
#if region_number_ChrY > 0:
#	if start3 > 0 and start4 > 0:
#		writing_total('Chromosome_Y','cnv_analysis_ChrY_sorted.txt', 0.3, 0.7, 1.3, 1.7, last_col_2_hide, last_col, start3, start4)
#	else:
#		writing_total('Chromosome_Y','cnv_analysis_ChrY_sorted.txt', 0.3, 0.7, 1.3, 1.7, last_col_2_hide, last_col, start1, start2)
#
#if (Panel != False):
#	worksheet2.protect()
#summary.protect()
#workbook.close()
##remove temporary files
#os.system("rm cnv_analysis.txt cnv_analysis_sorted.txt cnv_analysis_ChrX.txt cnv_analysis_ChrX_sorted.txt cnv_analysis_ChrY.txt cnv_analysis_ChrY_sorted.txt")
#print("\nDone!!!\n")
