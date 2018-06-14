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
#from operator import itemgetter

# ==============================================================================

pp = pprint.PrettyPrinter(indent=4, depth=6)
#added david 27/03/2018 - deal with path passed as arg
parser = argparse.ArgumentParser(usage='python MobiCNV.py [-i PATH/TO/DoC_FILE/DIR/ -t (tsv|csv) -p PATH/TO/GENE/TXT/FILE -o PATH/TO/OUTPUT/FILE]', description='Compute CNV predictions from Illumina Depth of Coverage files (generated by MSR or LRM in DNA enrichment mode. Can also deal with samtools bedcov slightly modified files. Generates an Excel File with all computed values in Autosomes, X and Y sheets if relevant, in addition to a summary sheet which displays aberrant regions +- 1. Optionally, an additional panel sheet focuses on a list of genes passed as argument as a txt file (one gene per line).')
parser.add_argument('-i', '--input', default='.')
parser.add_argument('-p', '--panel', default='')
parser.add_argument('-t', '--type', default='csv')
parser.add_argument('-o', '--output', default='mobicnv.xlsx')
args = parser.parse_args()
#print(args.path)
#sys.exit()
#Path = "/Users/david/Downloads/henri_cnv/03_2018/170925H27JNL/"
Path = args.input
Panel = args.panel
Type = args.type
OutFile = args.output
ext = "csv"
delim = ","
if Type == "tsv":
	ext = "tsv"
	delim = "\t"
if Panel == '':
	Panel = False

filelist = os.listdir(Path)
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

print("\nFiles that will be considered:\n")


def build_dict(baseName, line, region_number, prm, psm, key):
	if(baseName not in psm):
		psm[baseName] = {"rawDocSum": float(line[4])}
	else:
		psm[baseName]["rawDocSum"] += float(line[4])
	# création dictionnaire avec les valeurs rawDocs
	coordinate = (key,line[0],line[1],line[2],line[3])
	if coordinate not in prm :
		region_number += 1
		prm[coordinate] = {baseName : {"rawDoc" : float(line[4])}}
	else :
		prm[coordinate][baseName] = {"rawDoc" : float(line[4])}	
	return(region_number, prm, psm)

#############
for i in filelist:
	#### à fonctionnariser
	regex = re.compile(r'^([^\.].+)[\._]coverage\.%s$'%ext)
	matchObj = regex.search(os.path.basename(i))
	if(matchObj):
		#for each file
		print(i)
		number_of_file += 1
		baseName = matchObj.group(1)
		keyAuto = keyX = keyY = 1
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
				if not re.match(expression, line[0]) and line[0] != "chrX" and line[0] != "chrY":
					(region_number, per_region_metrics, per_sample_metrics) = build_dict(baseName, line, region_number, per_region_metrics, per_sample_metrics, keyAuto)
					keyAuto+=1
				#chrX
				elif not re.match(expression, line[0]) and (line[0] == "chrX" or line[0] == "X"):
					(region_number_ChrX, per_region_metrics_ChrX, per_sample_metrics_ChrX) = build_dict(baseName, line, region_number_ChrX, per_region_metrics_ChrX, per_sample_metrics_ChrX, keyX)
					keyX+=1
				#chrY
				elif not re.match(expression, line[0]) and (line[0] == "chrY" or line[0] == "Y"):
					(region_number_ChrY, per_region_metrics_ChrY, per_sample_metrics_ChrY) = build_dict(baseName, line, region_number_ChrY, per_region_metrics_ChrY, per_sample_metrics_ChrY, keyY)
					keyY+=1

#############

#############
#pp.pprint(per_sample_metrics)
#pp.pprint(per_region_metrics_ChrX)
#sys.exit()


# Itération sur le dictionnaire patient pour calculer la moyenne par  patient
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
				print("\n\nWARNING Gender inconsistancy for " + sample_name + "reads on Y chr with X ratio > 0.65\n\n")
				per_sample_metrics_ChrX[sample_name]["gender"] = "male/female"
#############
# Itération sur le dictionnaire patient pour calculer la moyenne par exon et l'exon normalisé

#pp.pprint(per_sample_metrics_ChrX)
#sys.exit()
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

per_region_metrics = exon_mean(per_region_metrics)
if region_number_ChrX > 0:
	per_region_metrics_ChrX = exon_mean(per_region_metrics_ChrX)
if region_number_ChrY > 0:
	per_region_metrics_ChrY = exon_mean(per_region_metrics_ChrY)


#############
# Itération sur le dictionnaire régions pour calculer la somme moyenne par patient
# sauf pour les exons du patient - Somme ajoutée dans le dict patients


def compute_ratio(psm, prm, region_number):
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
	xfactor = 2.5
	for coordinate in prm :
	  	for sample_name in prm[coordinate]:
	  		ratio = prm[coordinate][sample_name]["normalisedRatio"]
	  		if ratio < 0.3:
	  			prm[coordinate][sample_name]["MobiAdvice"] = "HomDel"
	  		elif ratio > 1.7:
	  			prm[coordinate][sample_name]["MobiAdvice"] = "HomDup"
	  		elif ratio > 0.8 and ratio < 1.2:
	  			prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
	  		elif ratio >= 0.3 and ratio <= 0.8:
	  			dynamic_threshold = 1 - (xfactor * prm[coordinate][sample_name]["ratioStdev"])
	  			if ratio < dynamic_threshold:
	  				prm[coordinate][sample_name]["MobiAdvice"] = "HetDel"
	  			else:
	  				prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
	  		elif ratio >= 1.2 and ratio <= 1.7:
	  			dynamic_threshold = 1 + (xfactor * prm[coordinate][sample_name]["ratioStdev"])
	  			if ratio > dynamic_threshold:
	  				prm[coordinate][sample_name]["MobiAdvice"] = "HetDup"
	  			else:
	  				prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
	return (psm, prm)

(per_sample_metrics, per_region_metrics) = compute_ratio(per_sample_metrics, per_region_metrics, region_number)
if region_number_ChrX > 0:
	(per_sample_metrics_ChrX, per_region_metrics_ChrX) = compute_ratio(per_sample_metrics_ChrX, per_region_metrics_ChrX, region_number_ChrX)
if region_number_ChrY > 0:
	(per_sample_metrics_ChrY, per_region_metrics_ChrY) = compute_ratio(per_sample_metrics_ChrY, per_region_metrics_ChrY, region_number_ChrY)

#pp.pprint(per_region_metrics)
#pp.pprint(per_region_metrics_ChrX)
#sys.exit()




#def write_csv_file(file_in, file_out, dict_r, dict_m):
#	with open(file_in, 'w') as csv_file:
#		#csv.writer(csv_file)
#		header = "Chr\tStart\tEnd\tRegionID\tMean DoC\t"
#		for coordinate in dict_r :
#			for sample_name in dict_r[coordinate]:
#				header += str(sample_name) + "_rawDoc" + "\t"
#			for sample_name in dict_r[coordinate]:
#			 	header += str(sample_name) + "_moy_exon" + "\t"
#			for sample_name in dict_r[coordinate]:
#			 	header += str(sample_name) + "_patient_normalise" + "\t"
#			for sample_name in dict_r[coordinate]:
#			 	header += str(sample_name) + "_exon_normalise" + "\t"
#			for sample_name in dict_r[coordinate]:
#			 	header += str(sample_name) + "_stdev" + "\t"
#			for sample_name in dict_r[coordinate]:
#			 	header += str(sample_name) + "_ratio" + "\t"
#			break
#		csv_file.write(header + "\n")
#		for coordinate in dict_r :
#			csv_file.write(
#				str(coordinate[0]) + "\t" +
#				str(coordinate[1]) + "\t" +
#				str(coordinate[2]) + "\t" +
#				str(coordinate[3]) + "\t" +
#				str((dict_m[coordinate]["exonMeanDoc"])) + "\t"
#				)
#			for sample_name in dict_r[coordinate] :
#				csv_file.write(str(dict_r[coordinate][sample_name]["rawDoc"]) + "\t")
#			for sample_name in dict_r[coordinate] :
#				csv_file.write(str(round(dict_r[coordinate][sample_name]["exonMeanOtherSamples"], 2)) + "\t")
#			for sample_name in dict_r[coordinate] :
#				csv_file.write(str(round(dict_r[coordinate][sample_name]["patient_normalise"],2)) + "\t")
#			for sample_name in dict_r[coordinate] :
#				csv_file.write(str(round(dict_r[coordinate][sample_name]["exon_normalise"],2)) + "\t")
#			for sample_name in dict_r[coordinate] :
#				csv_file.write(str(round(dict_r[coordinate][sample_name]["stdev"],3)) + "\t")
#			for sample_name in dict_r[coordinate] :
#				csv_file.write(str(round(dict_r[coordinate][sample_name]["ratio_normalise"],3)) + "\t")
#			csv_file.write("\n")
#	os.system("sort -k1.4n -k2,2n -k3,3n " + file_in + " > " + file_out)
#
#write_csv_file('cnv_analysis.txt', 'cnv_analysis_sorted.txt', per_region_metrics, dict_mean)
#if region_number_ChrX > 0:
#	write_csv_file('cnv_analysis_ChrX.txt', 'cnv_analysis_ChrX_sorted.txt', per_region_metrics_ChrX, dict_mean)
#if region_number_ChrY > 0:
#	write_csv_file('cnv_analysis_ChrY.txt', 'cnv_analysis_ChrY_sorted.txt', per_region_metrics_ChrY, dict_mean)

# with open('cnv_analysis.txt', 'w') as csv_file:
# 	writer = csv.writer(csv_file)
# 	header = "Chr\tStart\tEnd\tRegionID\tMean DoC\t"
# 	for coordinate in per_region_metrics :
# 		for sample_name in per_region_metrics[coordinate]:
# 			header += str(sample_name) + "_rawDoc" + "\t"
# 		for sample_name in per_region_metrics[coordinate]:
# 			header += str(sample_name) + "_moy_exon" + "\t"
# 		for sample_name in per_region_metrics[coordinate]:
# 			header += str(sample_name) + "_patient_normalise" + "\t"
# 		for sample_name in per_region_metrics[coordinate]:
# 			header += str(sample_name) + "_exon_normalise" + "\t"
# 		for sample_name in per_region_metrics[coordinate]:
# 			header += str(sample_name) + "_ratio" + "\t"
# 		break
# 	csv_file.write(header + "\n")
# 	for coordinate in per_region_metrics :
# 		csv_file.write(
# 			str(coordinate[0]) + "\t" +
# 			str(coordinate[1]) + "\t" +
# 			str(coordinate[2]) + "\t" +
# 			str(coordinate[3]) + "\t" +
# 			str(round(dict_mean[coordinate]["exonMeanDoc"], 2)) + "\t"
# 			)
# 		for sample_name in per_region_metrics[coordinate] :
# 			csv_file.write(str(per_region_metrics[coordinate][sample_name]["rawDoc"]) + "\t")
# 		for sample_name in per_region_metrics[coordinate] :
# 			csv_file.write(str(round(per_region_metrics[coordinate][sample_name]["exonMeanOtherSamples"], 2)) + "\t")
# 		for sample_name in per_region_metrics[coordinate] :
# 			csv_file.write(str(round(per_region_metrics[coordinate][sample_name]["patient_normalise"],2)) + "\t")
# 		for sample_name in per_region_metrics[coordinate] :
# 			csv_file.write(str(round(per_region_metrics[coordinate][sample_name]["exon_normalise"],2)) + "\t")
# 		for sample_name in per_region_metrics[coordinate] :
# 			csv_file.write(str(round(per_region_metrics[coordinate][sample_name]["ratio_normalise"],3)) + "\t")
# 		csv_file.write("\n")

# with open('cnv_analysis_ChrX.txt', 'w') as csv_file:
# 	writer = csv.writer(csv_file)
# 	header = "Chr\tStart\tend\tRegionID\tMean DoC\t"
# 	for coordinate in per_region_metrics_ChrX :
# 		for sample_name in per_region_metrics_ChrX[coordinate]:
# 			header += str(sample_name) + "_rawDoc" + "\t"
# 		for sample_name in per_region_metrics_ChrX[coordinate]:
# 			header += str(sample_name) + "_moy_exon" + "\t"
# 		for sample_name in per_region_metrics_ChrX[coordinate]:
# 			header += str(sample_name) + "_patient_normalise" + "\t"
# 		for sample_name in per_region_metrics_ChrX[coordinate]:
# 			header += str(sample_name) + "_exon_normalise" + "\t"
# 		for sample_name in per_region_metrics_ChrX[coordinate]:
# 			header += str(sample_name) + "_ratio" + "\t"
# 		break
# 	csv_file.write(header + "\n")
# 	for coordinate in per_region_metrics_ChrX :
# 		csv_file.write(
# 			str(coordinate[0]) + "\t" +
# 			str(coordinate[1]) + "\t" +
# 			str(coordinate[2]) + "\t" +
# 			str(coordinate[3]) + "\t" +
# 			str(round(dict_mean[coordinate]["exonMeanDoc"], 2)) + "\t"
# 			)
# 		for sample_name in per_region_metrics_ChrX[coordinate] :
# 			csv_file.write(str(per_region_metrics_ChrX[coordinate][sample_name]["rawDoc"]) + "\t")
# 		for sample_name in per_region_metrics_ChrX[coordinate] :
# 			csv_file.write(str(round(per_region_metrics_ChrX[coordinate][sample_name]["exonMeanOtherSamples"], 2)) + "\t")
# 		for sample_name in per_region_metrics_ChrX[coordinate] :
# 			csv_file.write(str(round(per_region_metrics_ChrX[coordinate][sample_name]["patient_normalise"],2)) + "\t")
# 		for sample_name in per_region_metrics_ChrX[coordinate] :
# 			csv_file.write(str(round(per_region_metrics_ChrX[coordinate][sample_name]["exon_normalise"],2)) + "\t")
# 		for sample_name in per_region_metrics_ChrX[coordinate] :
# 			csv_file.write(str(round(per_region_metrics_ChrX[coordinate][sample_name]["ratio_normalise"],3)) + "\t")
# 		csv_file.write("\n")

# Sorting files, works on Mac only, if you are poor (like Charles) with other OS, go fuck yourself
#os.system("sort -k1.4n -k2,2n -k3,3n cnv_analysis.txt > cnv_analysis_sorted.txt")
#os.system("sort -k1.4n -k2,2n -k3,3n cnv_analysis_ChrX.txt > cnv_analysis_ChrX_sorted.txt")
#


#We build a small list of genes of interest
###########
if (Panel != False):
	panel = open(Panel, 'r')
	# panelreader = csv.DictReader(panel, delimiter='\t')
	liste_panel = []
	for gene in panel :
		liste_panel.append(gene.rstrip())
###########


### Excel conversion


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
#http://xlsxwriter.readthedocs.io/example_conditional_format.html#ex-cond-format
# Add a format. Light red fill with dark red text.
format1 = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
# Add a format. Green fill with dark green text.
format2 = workbook.add_format({'bg_color': '#C6EFCE', 'font_color': '#006100'})

last_col = number_of_file*6 + 5
#sample number x number of data blocks + 6 first cols to show - 1 coz numbering begins at 0 (A = 0)
last_col_2_hide = number_of_file*5 + 5
#sample number x number of data blocks to hide + 6 first cols to show - 1 coz numbering begins at 0 (A = 0)

def format_sheet(sheet, last_col):
	sheet.freeze_panes(1, 6)
	sheet.set_row(0, 20, style5)
	sheet.set_column('A:F', 15, style5)
	sheet.set_column('E:E', 25)
	sheet.set_column('F:F', 15, style6)
	sheet.set_column(last_col+2, last_col+2, 15)
	sheet.set_column(last_col+3, last_col+3, 15)
	sheet.write(3, last_col+2, 'Legend:', style5)
	sheet.write(4, last_col+2, '', style1)
	sheet.write(5, last_col+2, '', style2)
	sheet.write(6, last_col+2, '', style3)
	sheet.write(7, last_col+2, '', style4)
	sheet.write(4, last_col+3, 'Homozygous deletion', style5)
	sheet.write(5, last_col+3, 'Heterozygous deletion', style5)
	sheet.write(6, last_col+3, 'Heterozygous duplication', style5)
	sheet.write(7, last_col+3, 'Homozygous duplication', style5)


#worksheet summary to get only interesting stuff
#summary = workbook.add_worksheet(str("Summary"))
#format_sheet(summary, last_col)

#if (Panel != False):
#	worksheet2 = workbook.add_worksheet(str("Panel"))
#	format_sheet(worksheet2, last_col)
	# worksheet2.freeze_panes(1, 5)
	# worksheet2.set_row(0, 20, style5)
	# worksheet2.set_column('A:E', 15, style5)
	# worksheet2.set_column('D:D', 25)

def add_conditionnal_format(worksheet, threshold, start, end):
	#add a conditionnal format to raw DoC
	if (start == 0):
		start = 2
	cell_range = "F" + str(start) + ":F" + str(end)
	#http://xlsxwriter.readthedocs.io/working_with_conditional_formats.html
	worksheet.conditional_format(
		cell_range,
		{'type': 'icon_set',
		 'icon_style': '3_traffic_lights',
		 'icons': [{'criteria': '>=', 'type': 'number',     'value': 100},
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



def print_worksheet(name, last_col, last_col_2_hide, workbook, prm, quality, reduced_regions):
	#sheet creation
	worksheet = workbook.add_worksheet(str(name))
	format_sheet(worksheet, last_col)
	#i=first row
	i=0
	#dict iterations
	for region in sorted(prm):
		#j=col
		j=0
		if i == 0:
			#sheet headers
			headers = ("Region Number", "Chromosome", "Start", "End", "Annotation", "Mean DoC")
			for text in headers:
				worksheet.write(i, j, text, style5)
				j+=1
			for sample in prm[region]:
				worksheet.write(i, j, sample + "_rawDoc", style5)
				worksheet.write(i, j+number_of_file, sample + "_regionMeanOtherSamples", style5)
				worksheet.write(i, j+(2*number_of_file), sample + "_normalisedRegion", style5)
				worksheet.write(i, j+(3*number_of_file), sample + "_normalisedMeanOtherSamples", style5)
				worksheet.write(i, j+(4*number_of_file), sample + "_ratioStdev", style5)
				worksheet.write(i, j+(5*number_of_file), sample + "_normalisedRatio", style5)
				j+=1
		else:
			for row_header in region:
				worksheet.write(i, j, row_header)
				j+=1
			for sample in prm[region]:
				j+=1
				worksheet.write(i, j, prm[region][sample]["rawDoc"])
				worksheet.write(i, j+number_of_file, prm[region][sample]["regionMeanOtherSamples"])
				worksheet.write(i, j+(2*number_of_file), prm[region][sample]["normalisedRegion"])
				worksheet.write(i, j+(3*number_of_file), prm[region][sample]["normalisedMeanOtherSamples"])
				worksheet.write(i, j+(4*number_of_file), prm[region][sample]["ratioStdev"])
				#define cell style - default style5: just bold
				cell_style = style5
				if quality == "global":
					reduced_regions[region] = {sample: prm[region][sample]}

				if prm[region][sample]["MobiAdvice"] == "HomDel":
					cell_style = style1
				elif prm[region][sample]["MobiAdvice"] == "HetDel":
					cell_style = style2
				elif prm[region][sample]["MobiAdvice"] == "HetDup":
					cell_style = style3
				elif prm[region][sample]["MobiAdvice"] == "HomDup":
					cell_style = style4
				worksheet.write(i, j+(5*number_of_file), prm[region][sample]["normalisedRatio"], cell_style)
				last_sample = sample
			#mean global doc for region
			worksheet.write(i, 5, prm[region][last_sample]["regionMeanDoc"])
		i+=1

	worksheet.set_column(6, last_col_2_hide, None, None, {'level': 1, 'hidden': True})
	add_conditionnal_format(worksheet, 50, 2, len(list(prm)))
	worksheet.protect()

	if quality == "global":
		return reduced_regions

print("\nBuilding Excel File:")
print("Autosomes worksheet...")
summary_regions = {}
summary_regions = print_worksheet('Autosomes', last_col, last_col_2_hide, workbook, per_region_metrics, 'global', summary_regions)

if region_number_ChrX > 0:
	print("ChrX worksheet...")
	summary_regions = print_worksheet('Chromosome X', last_col, last_col_2_hide, workbook, per_region_metrics_ChrX, 'global', summary_regions)
if region_number_ChrY > 0:
	print("ChrY worksheet...")
	summary_regions = print_worksheet('Chromosome Y', last_col, last_col_2_hide, workbook, per_region_metrics_ChrY, 'global', summary_regions)

print("Summary worksheet...")
print_worksheet('Summary', last_col, last_col_2_hide, workbook, summary_regions, 'reduced', '')


workbook.close()


sys.exit()







def write_small_worksheets(selected, start, first_row, small_worksheet, col_list, last_col, regex_r, threshold_del_hmz, threshold_del_htz, threshold_dup_htz, threshold_dup_hmz):
	#called inside writing_total to wirte summary and panel sheets
	i = 0
	uniq_selected = list(set(selected))
	for column in col_list:
		j = first_row
		if(start > 0):
			j = start
		for item in range(len(column)):
			if (item in uniq_selected):
				if regex_r.search(column[0]):
					if(item > 0):
						if(float(column[item]) <= threshold_del_hmz):
							small_worksheet.write(j, i , column[item],style1)
						elif(float(column[item]) <= threshold_del_htz):
							small_worksheet.write(j, i, column[item],style2)
						elif(float(column[item]) <= threshold_dup_htz):
							small_worksheet.write(j, i, column[item], style5)
						elif(float(column[item]) <= threshold_dup_hmz):
							small_worksheet.write(j, i, column[item],style3)
						else :
							small_worksheet.write(j, i, column[item],style4)
					else:
						small_worksheet.write(j, i, column[item], style5)
				else:
					try:
						small_worksheet.write(j,i,int(column[item]))
					except ValueError:
						small_worksheet.write(j,i,column[item])
				j+=1
		i+=1
	small_worksheet.write(9, last_col+2, "Sample ID", style5)
	small_worksheet.write(9, last_col+3, "Predicted Gender", style5)
	small_worksheet.write(9, last_col+4, "X ratio", style5)
	if region_number_ChrX > 0 or region_number_ChrY > 0:
		m = 10
		for sample_name in dict_gender:
			style = style8
			if dict_gender[sample_name]["gender"] == 'female':
				style = style7
			small_worksheet.write(m, last_col+2, sample_name, style)
			small_worksheet.write(m, last_col+3, dict_gender[sample_name]["gender"], style)
			small_worksheet.write(m, last_col+4, round(dict_gender[sample_name]["xratio"], 2), style)
			m += 1
	return (uniq_selected, j)



def writing_total(worksheet, txt_file, threshold_del_hmz, threshold_del_htz, threshold_dup_htz, threshold_dup_hmz, last_col_2_hide, last_col, start1=0, start2=0):
	# TODO: change it to parameter and save worksheet
	# if panel:

	worksheet = workbook.add_worksheet(str(worksheet))
	format_sheet(worksheet, last_col)
	# worksheet.freeze_panes(1, 5)
	# worksheet.set_row(0, 20, style5)
	# worksheet.set_column('A:E', 15, style5)
	# worksheet.set_column('D:D', 25)
	#structure data from txt
	f = open(str(txt_file), 'r+')
	row_list = []
	for row in f:
		row_list.append(row.split('\t'))
	# pp.pprint(row_list)

	column_list = zip(*row_list)
	column_list2 = zip(*row_list)
	column_list3 = zip(*row_list)
	i = 0
	l = 0
	interesting = []
	gene4interest = []
	regex_ratio = re.compile(r'(.*)_ratio$')
	regex_noCNV = re.compile(r'no CNV')
	regex_region = re.compile(r'RegionID')
	for column in column_list:
		for item in range(len(column)):
			if regex_region.search(column[0]):
				if (Panel != False):
					for gene in liste_panel:
						if re.compile(r'.*' + gene + '.*').search(column[item]) :
							# print (column[item])
							gene4interest.append(item)
				else:
					gene4interest.append(item)
			if regex_ratio.search(column[0]):
				if (item > 0) :
					if(float(column[item]) <= threshold_del_hmz):
						worksheet.write(item, i , column[item],style1)
						for k in range(item-1,item+2):
							if (k > 0):
								interesting.append(k)
					elif(float(column[item]) <= threshold_del_htz):
						worksheet.write(item, i, column[item],style2)
						for k in range(item-1,item+2):
							if (k > 0):
								interesting.append(k)
					elif(float(column[item]) <= threshold_dup_htz):
						worksheet.write(item, i, column[item], style5)
					elif(float(column[item]) <= threshold_dup_hmz):
						worksheet.write(item, i, column[item],style3)
						for k in range(item-1,item+2):
							if (k > 0):
								interesting.append(k)
					else :
						worksheet.write(item, i, column[item],style4)
						for k in range(item-1,item+2):
							if (k > 0):
								interesting.append(k)
				else:
					worksheet.write(item, i, column[item], style5)
			else:
				try:
					worksheet.write(item,i,int(column[item]))
				except ValueError:
					worksheet.write(item,i,column[item])
			if (item == 0):
				summary.write(item,i,column[item], style5)
				# if panel:
				if (Panel != False):
					worksheet2.write(item,i,column[item], style5)
		i+=1
	#i = 0
	(uniq_interesting, j) = write_small_worksheets(interesting, start1, 0, summary, column_list2, last_col, regex_ratio, threshold_del_hmz, threshold_del_htz, threshold_dup_htz, threshold_dup_hmz)
	if (Panel != False):
		(uniq_interesting_panel, l) = write_small_worksheets(gene4interest, start2, 1, worksheet2, column_list3, last_col,  regex_ratio, threshold_del_hmz, threshold_del_htz, threshold_dup_htz, threshold_dup_hmz)
		worksheet2.set_column(5,last_col_2_hide, None, None, {'level': 1, 'hidden': True})
		add_conditionnal_format(worksheet2, 50, start2, start2 + len(gene4interest))
	# uniq_interesting = list(set(interesting))
	# for column in column_list2:
	# 	j = 0
	# 	if(start1 > 0):
	# 		j = start1
	# 	for item in range(len(column)):
	# 		if (item in uniq_interesting):
	# 			if regex_ratio.search(column[0]):
	# 				if(item > 0):
	# 					if(float(column[item]) <= threshold_del_hmz):
	# 						summary.write(j, i , column[item],style1)
	# 					elif(float(column[item]) <= threshold_del_htz):
	# 						summary.write(j, i, column[item],style2)
	# 					elif(float(column[item]) <= threshold_dup_htz):
	# 						summary.write(j, i, column[item], style5)
	# 					elif(float(column[item]) <= threshold_dup_hmz):
	# 						summary.write(j, i, column[item],style3)
	# 					else :
	# 						summary.write(j, i, column[item],style4)
	# 				else:
	# 					summary.write(j, i, column[item], style5)
	# 			else:
	# 				summary.write(j,i,column[item])
	# 			j+=1
	# 	i+=1
	

	# i = 0
	# # part to dev
	# if (Panel != False):
	# 	#l=0
	# 	uniq_interesting_panel = list(set(gene4interest))
	# 	for column in column_list3:
	# 		l = 1
	# 		if (start2 > 0 ):
	# 			l = start2
	# 		for item in range(len(column)):
	# 			if (item in uniq_interesting_panel):
	# 				if regex_ratio.search(column[0]):
	# 					if(item > 0):
	# 						if(float(column[item]) <= threshold_del_hmz):
	# 							worksheet2.write(l, i , column[item],style1)
	# 						elif(float(column[item]) <= threshold_del_htz):
	# 							worksheet2.write(l, i, column[item],style2)
	# 						elif(float(column[item]) <= threshold_dup_htz):
	# 							worksheet2.write(l, i, column[item], style5)
	# 						elif(float(column[item]) <= threshold_dup_hmz):
	# 							worksheet2.write(l, i, column[item],style3)
	# 						else :
	# 							worksheet2.write(l, i, column[item],style4)
	# 					else:
	# 						worksheet2.write(l, i, column[item], style5)
	# 				else:
	# 					worksheet2.write(l,i,column[item])
	# 				l+=1
	# 		i+=1
	# 	i = 0
	# 	worksheet2.set_column('F:BM', None, None, {'level': 1, 'hidden': True})
	# 	add_conditionnal_format(worksheet2, 50, start2, start2 + len(gene4interest))
	worksheet.set_column('F:BM', None, None, {'level': 1, 'hidden': True})
	worksheet.set_column(5, last_col_2_hide, None, None, {'level': 1, 'hidden': True})
	add_conditionnal_format(worksheet, 50, 2, len(row_list))
	worksheet.protect()
	summary.set_column(5, last_col_2_hide, None, None, {'level': 1, 'hidden': True})
	add_conditionnal_format(summary, 50, start1, start1 + len(uniq_interesting))
	return (j,l)

#worksheet for autosomes
(start1, start2) = writing_total('Autosomes','cnv_analysis_sorted.txt', 0.3, 0.7, 1.3, 1.7, last_col_2_hide, last_col)
# print (start1, start2)
#worksheet for ChrX
start3 = 0
start4 = 0
if region_number_ChrX > 0:
	(start3, start4) = writing_total('Chromosome_X','cnv_analysis_ChrX_sorted.txt', 0.3, 0.7, 1.3, 1.7, last_col_2_hide, last_col, start1, start2)
if region_number_ChrY > 0:
	if start3 > 0 and start4 > 0:
		writing_total('Chromosome_Y','cnv_analysis_ChrY_sorted.txt', 0.3, 0.7, 1.3, 1.7, last_col_2_hide, last_col, start3, start4)
	else:
		writing_total('Chromosome_Y','cnv_analysis_ChrY_sorted.txt', 0.3, 0.7, 1.3, 1.7, last_col_2_hide, last_col, start1, start2)

if (Panel != False):
	worksheet2.protect()
summary.protect()
workbook.close()
#remove temporary files
os.system("rm cnv_analysis.txt cnv_analysis_sorted.txt cnv_analysis_ChrX.txt cnv_analysis_ChrX_sorted.txt cnv_analysis_ChrY.txt cnv_analysis_ChrY_sorted.txt")
print("\nDone!!!\n")
