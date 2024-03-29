##!/usr/local/bin/ python3.5
# coding: utf-8

# ==============================================================================

# IMPORT LIBRARIES
import sys        # system command
import csv        # read, write csv
import re        # regex
import argparse    # for options
import os        # for options
# import pprint   # print data structure
import numpy as np
# import math
import xlsxwriter
import vcf
# from operator import itemgetter

# ==============================================================================




# funcion to build differents dicts from coverage files


def build_dict(sample, line, region_number, prm, psm, key, chr_semaph):
    if(sample not in psm):
        psm[sample] = {"rawDocSum": float(line[4])}
    else:
        psm[sample]["rawDocSum"] += float(line[4])
    # deal with chr1 or 1
    chrom = line[0]
    chr_reg = r'^chr'
    if not re.match(chr_reg, chrom):
        chr_semaph = True
    if chr_semaph:
        chrom = "chr" + chrom
    # create dictionnary with rawDoc values per sample per region
    coordinate = (key, chrom, line[1], line[2], line[3])
    if coordinate not in prm:
        region_number += 1
        prm[coordinate] = {sample: {"rawDoc": float(line[4])}}
    else:
        prm[coordinate][sample] = {"rawDoc": float(line[4])}
    return(region_number, prm, psm, chr_semaph)

#############

# function to populate per_region_metrics dictionnary with mean coverage per exons
# two metrics:
# 1- full mean including all samples per region just to be printed in final file for informative purpose
# 2- mean excluding current sample for normalisation purpose


def exon_mean(prm):
    for coordinate in prm:
        for sample_name in prm[coordinate]:
            total = full_total = sample_number = 0
            for sample_other in prm[coordinate]:
                if sample_other != sample_name and \
                        float(prm[coordinate][sample_other]['rawDoc']) > 0:
                    total += prm[coordinate][sample_other]['rawDoc']
                    sample_number += 1
                full_total += prm[coordinate][sample_other]['rawDoc']
            # regionMeanOtherSamples = total / sample_number
            prm[coordinate][sample_name]["regionMeanDoc"] = int(full_total / (sample_number+1))
            if sample_number > 0:
                prm[coordinate][sample_name]["regionMeanOtherSamples"] = round(float(total / sample_number), 3)
            else:
                prm[coordinate][sample_name]["regionMeanOtherSamples"] = 0
    return(prm)

#############


def compute_ratio(psm, prm, region_number, VcfDir, variants, chr_type, het_high, het_low, hom_high, hom_low):
    # loop to compute per region the mean coverage of all regions except the ROI
    for sample_name in psm:
        psm[sample_name]["regionMeanOtherSamplesSum"] = 0
        for coordinate in prm:
            # compute mean sum per exon
            psm[sample_name]["regionMeanOtherSamplesSum"] += prm[coordinate][sample_name]["regionMeanOtherSamples"]
        # mean for all samples except current sample
        psm[sample_name]["totalMeanOtherSample"] = round(psm[sample_name]["regionMeanOtherSamplesSum"] / region_number, 3)
    # loop to compute normalised values per region (vakue for the ROI and for all others excluding ROI)
    for sample_name in psm:
        for coordinate in prm:
            # normalisation per exon
            try:
                prm[coordinate][sample_name]["normalisedMeanOtherSamples"] = round(prm[coordinate][sample_name]["regionMeanOtherSamples"] / psm[sample_name]["totalMeanOtherSample"], 3)
            except ZeroDivisionError:
                prm[coordinate][sample_name]["normalisedMeanOtherSamples"] = float(0)
            try:
                prm[coordinate][sample_name]["normalisedRegion"] = round(float(prm[coordinate][sample_name]['rawDoc'] / psm[sample_name]["meanRawDoc"]), 3)
            except ZeroDivisionError:
                prm[coordinate][sample_name]["normalisedRegion"] = float(0)
    # computes final ratio
    for coordinate in prm:
        for sample_name in prm[coordinate]:
            try:
                normalised_ratio = float(prm[coordinate][sample_name]["normalisedRegion"]) / float(prm[coordinate][sample_name]["normalisedMeanOtherSamples"])
            except ZeroDivisionError:
                normalised_ratio = float(0)
            prm[coordinate][sample_name]["normalisedRatio"] = round(float(normalised_ratio), 3)
    # computes stdev for final ratio
    for coordinate in prm:
        for sample_name in prm[coordinate]:
            std_dev = []
            for sample_name2 in prm[coordinate]:
                if(sample_name != sample_name2) and \
                        float(prm[coordinate][sample_name2]["normalisedRatio"]) > 0:
                    std_dev.append(prm[coordinate][sample_name2]["normalisedRatio"])
            if std_dev == []:
                std_dev.append(0)
            prm[coordinate][sample_name]["ratioStdev"] = round(float(np.std(std_dev)), 3)
    # interpretation
    # <0.3 => hom del
    # >1.7 => hom dup
    # betwwen 0.8 and 1.2 => normal
    # between 0.3 and 0.8 =>supect het del then
    #    between 1-2.5sigma and 1 => normal
    #    <1-2.5sigma => het del
    # between 1.2 and 1.7 =>supect het del then
    #    between 1 and 1+2.5sigma => normal
    #    >1+2.5sigma => het dup
    # het_high = 1.3
    # het_low = 0.7
    # hom_high = 1.7
    # hom_low = 0.3
    xfactor = 2
    # dev purpose
    r = 0
    s = 0
    t = 0
    u = 0
    v = 0
    w = 0
    for coordinate in prm:
        # print(coordinate)
        for sample_name in prm[coordinate]:
            ratio = prm[coordinate][sample_name]["normalisedRatio"]
            if ratio < hom_low:
                prm[coordinate][sample_name]["MobiAdvice"] = "HomDel"
                v += 1
            elif ratio > hom_high:
                prm[coordinate][sample_name]["MobiAdvice"] = "HomDup"
                w += 1
            elif ratio > het_low and ratio < het_high:
                prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
            elif ratio >= hom_low and ratio <= het_low:
                dynamic_threshold = 1 - (xfactor * prm[coordinate][sample_name]["ratioStdev"])
                if ratio < dynamic_threshold:
                    vcf_semaph = 0
                    # check against VCF info (het variants) if relevant
                    if VcfDir is not False:
                        for var_pos in variants:
                            # chr must match
                            if var_pos[0] == coordinate[1] and int(var_pos[1]) >= int(coordinate[2]) and int(var_pos[1]) <= int(coordinate[3]) and sample_name in variants[var_pos]:
                                prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
                                vcf_semaph = 1
                                r += 1
                                # print(var_pos[0], coordinate[1], var_pos[1], coordinate[2], coordinate[3])
                                break
                    if vcf_semaph == 0:
                        prm[coordinate][sample_name]["MobiAdvice"] = "HetDel"
                        s += 1
                else:
                    prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
            elif ratio >= het_high and ratio <= hom_high:
                dynamic_threshold = 1 + (xfactor * prm[coordinate][sample_name]["ratioStdev"])
                if ratio > dynamic_threshold:
                    vcf_semaph = 0
                    # check against VCF info (het variants) if relevant
                    if VcfDir is not False:
                        for var_pos in variants:
                            # chr must match
                            if var_pos[0] == coordinate[1] and int(var_pos[1]) >= int(coordinate[2]) and int(var_pos[1]) <= int(coordinate[3]) and sample_name in variants[var_pos]:
                                if variants[var_pos][sample_name]['AB'] > 0.4 and variants[var_pos][sample_name]['AB'] < 0.6:
                                    # if variants[var_pos][key] > 0.4 and variants[var_pos][key] < 0.6:
                                    prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
                                    vcf_semaph = 1
                                    t += 1
                                    # print(var_pos[0] + "-" + str(var_pos[1]) + "-" + sample_name + "-" + str(variants[var_pos][sample_name]['AB']));
                                    break
                    if vcf_semaph == 0:
                        prm[coordinate][sample_name]["MobiAdvice"] = "HetDup"
                        u += 1
                    # prm[coordinate][sample_name]["MobiAdvice"] = "HetDup"
                else:
                    prm[coordinate][sample_name]["MobiAdvice"] = "Normal"
    print(chr_type + " - HetDels: " + str(s) + " vcfed: " + str(r))
    print(chr_type + " - HomDels: " + str(v))
    print(chr_type + " - HetDups: " + str(u) + " vcfed: " + str(t))
    print(chr_type + " - HomDups: " + str(w))
    return (psm, prm)

#############

# function to basically format new Excel sheets


def format_sheet(sheet, last_col, styles):
    sheet.freeze_panes(1, 6)
    sheet.set_row(0, 20, styles[4]) # style5
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

# function to add conditional formatting for mean global DoC of the regions (3_traffic_lights)


def add_conditionnal_format(worksheet, threshold, start, end):
    # add a conditionnal format to mean DoC
    if (start == 0):
        start = 2
    cell_range = "F" + str(start) + ":F" + str(end)
    # http://xlsxwriter.readthedocs.io/working_with_conditional_formats.html
    worksheet.conditional_format(
        cell_range,
        {'type': 'icon_set',
         'icon_style': '3_traffic_lights',
         'icons': [{'criteria': '>=', 'type': 'number', 'value': 100},
        {'criteria': '<=',  'type': 'number', 'value': threshold}]}
    )
    # xlsxwriter < 1.0.0 => cell formatting
#     worksheet.conditional_format(cell_range, {'type': 'cell',
#                                           'criteria': '>=',
#                                           'value': threshold,
#                                           'format': format2})
#     worksheet.conditional_format(cell_range, {'type': 'cell',
#                                           'criteria': '<',
#                                           'value': threshold,
#                                           'format': format1})
#############

# function to actually create result worksheets


def print_worksheet(name, last_col, last_col_2_hide, workbook, prm, quality, reduced_regions, panel_regions, low_cov_regions, Panel, panel_list, psmX, styles, number_of_file):
    # sheet creation
    worksheet = workbook.add_worksheet(str(name))
    format_sheet(worksheet, last_col, styles)
    if quality == "summary":
        worksheet.activate()
    # i=first row
    i = j = 0
    # dict iterations
    for region in sorted(prm):
        # j=col
        # j=0
        if i == 0:
            # sheet headers
            headers = ("Region Number", "Chromosome", "Start", "End", "Annotation", "Mean DoC")
            for text in headers:
                worksheet.write(i, j, text, styles[4])
                j += 1
            for sample in prm[region]:
                worksheet.write(i, j, sample + "_rawDoc", styles[4])
                worksheet.write(i, j+number_of_file, sample + "_regionMeanOtherSamples", styles[4])
                worksheet.write(i, j+(2*number_of_file), sample + "_normalisedRegion", styles[4])
                worksheet.write(i, j+(3*number_of_file), sample + "_normalisedMeanOtherSamples", styles[4])
                worksheet.write(i, j+(4*number_of_file), sample + "_ratioStdev", styles[4])
                worksheet.write(i, j+(5*number_of_file), sample + "_normalisedRatio", styles[4])
                j += 1
            i += 1
        # else:
        j = 0
        for row_header in region:
            worksheet.write(i, j, row_header)
            if j == 4 and Panel:
                for gene in panel_list:
                    if re.compile(r'.*' + gene + '.*').search(row_header):
                        panel_regions[region] = prm[region]
            j += 1
        m = 10
        for sample in prm[region]:
            j += 1
            if i == 1 and sample in psmX:
                style = styles[7]
                if psmX[sample]["gender"] == 'female':
                    style = styles[6]
                worksheet.write(m, last_col+2, sample, style)
                worksheet.write(m, last_col+3, psmX[sample]["gender"], style)
                worksheet.write(m, last_col+4, round(psmX[sample]["xratio"], 2), style)
                m += 1
            worksheet.write(i, j, prm[region][sample]["rawDoc"])
            worksheet.write(i, j+number_of_file, prm[region][sample]["regionMeanOtherSamples"])
            worksheet.write(i, j+(2*number_of_file), prm[region][sample]["normalisedRegion"])
            worksheet.write(i, j+(3*number_of_file), prm[region][sample]["normalisedMeanOtherSamples"])
            worksheet.write(i, j+(4*number_of_file), prm[region][sample]["ratioStdev"])
            # define cell style - default style5: just bold
            # and prm[region][sample]["MobiAdvice"] != 'Normal'
            cell_style = styles[4]
            if quality == "global":
                if region not in reduced_regions:
                    # here we optionally check whether the sample's vcf contains heterozygous (PASS) variants in the region of interest
                    # here we check if number_of_sample > 4 and meanDoC all samples == 0 => we don't put that region in the summary sheet
                    if number_of_file <= 3 and prm[region][sample]["MobiAdvice"] != 'Normal':
                        reduced_regions[region] = prm[region]
                    elif prm[region][sample]["regionMeanDoc"] >= 100 and prm[region][sample]["MobiAdvice"] != 'Normal':
                        reduced_regions[region] = prm[region]
                    elif prm[region][sample]["regionMeanDoc"] < 100:
                        low_cov_regions[region] = prm[region]
            # LowCovSummary version
            # if quality == "global" and prm[region][sample]["MobiAdvice"] != 'Normal':
            #     if region not in reduced_regions:
            #         #here we optionally check whether the sample's vcf contains heterozygous (PASS) variants in the region of interest
            #         #here we check if number_of_sample > 4 and meanDoC all samples == 0 => we don't put that region in the summary sheet
            #         if number_of_file <= 3:
            #             reduced_regions[region] = prm[region]
            #         elif prm[region][sample]["regionMeanDoc"] >= 100:
            #             reduced_regions[region] = prm[region]
            #         elif prm[region][sample]["regionMeanDoc"] >= 0:
            #             low_cov_regions[region] = prm[region]
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
        # mean global doc for region
        worksheet.write(i, 5, prm[region][last_sample]["regionMeanDoc"])
        i += 1

    worksheet.set_column(6, last_col_2_hide, None, None, {'level': 1, 'hidden': True})
    add_conditionnal_format(worksheet, 50, 2, len(list(prm))+1)
    worksheet.protect()

    if quality == "global":
        return reduced_regions, low_cov_regions, panel_regions

#############


def main():

    # for developping purpose
    # pp = pprint.PrettyPrinter(indent=4, depth=6)

    # added david 27/03/2018 - deal args
    parser = argparse.ArgumentParser(usage='python MobiCNV.py [-i PATH/TO/DoC_FILE_DIR/ -t (tsv|csv) -p PATH/TO/GENE/TXT/FILE -o PATH/TO/OUTPUT/FILE -v PATH/TO/VCFs_DIR/]', description='Compute CNV predictions from Illumina Depth of Coverage files (generated by MSR or LRM in DNA enrichment mode. Can also deal with samtools bedcov slightly modified files. Generates an Excel File with all computed values in Autosomes, X and Y sheets if relevant, in addition to a summary sheet which displays aberrant regions +- 1. Optionally, an additional panel sheet focuses on a list of genes passed as argument as a txt file (one gene per line).')
    parser.add_argument('-i', '--input', default='.', required=True, help='Path to the directory containing the coverage files. The sample name will be deduced from the file names which should be sample_coverage.[ct]sv.')
    parser.add_argument('-p', '--panel', default='', help='If a text file with gene names is provided, and if the coverage files include gene annotation, a supplementary sheet will be created focusing on the genes of interest.')
    parser.add_argument('-t', '--type', default='csv', help='Can be csv or tsv.')
    parser.add_argument('-o', '--output', default='mobicnv.xlsx')
    parser.add_argument('-v', '--vcf', default='', help='Path to the directory containing optional VCF files. Can be the same as coverage directory, but must be specified for the VCF search to be activated. Data can be centralised in a single VCF, or in mutliple VCFs. MobiCNV will search for samples in any present VCF based on the name of the coverage files.')
    parser.add_argument('-hedu', '--het-dup', type=float, default='1.3', help='Threshold for heterozygous duplication ratio. Default 1.3.')
    parser.add_argument('-hede', '--het-del', type=float, default='0.7', help='Threshold for heterozygous deletion ratio. Default 0.7.')
    parser.add_argument('-hodu', '--hom-dup', type=float, default='1.7', help='Threshold for homozygous duplication ratio. Default 1.7.')
    parser.add_argument('-hode', '--hom-del', type=float, default='0.3', help='Threshold for homozygous deletion ratio. Default 0.3.')
    args = parser.parse_args()

    # Variables declaration
    # print(args.path)
    # sys.exit()
    # if !args.input:
    #    sys.exit('You must provide an input path, try python MobiCNV.py -h to get some help')
    if os.path.isdir(args.input):
        Path = args.input
    else:
        sys.exit('invalid input path, please check your command')
    Panel = args.panel
    Type = args.type
    OutFile = args.output
    VcfDir = args.vcf
    het_high = args.het_dup
    het_low = args.het_del
    hom_high = args.hom_dup
    hom_low = args.hom_del
    ext = "csv"
    delim = ","
    if Type == "tsv":
        ext = "tsv"
        delim = "\t"
    if Panel == '':
        Panel = False
    # if VcfDir == '':
    #    VcfDir = False
    if os.path.isdir(VcfDir):
        vcf_list = os.listdir(VcfDir)
    else:
        print
        VcfDir = False
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
        regex = re.compile(r'^([^\.].+)[\._]coverage\.%s$'%ext)
        match_file = regex.search(os.path.basename(i))
        if match_file:
            # for each file
            number_of_file += 1
            sample = match_file.group(1)
            # filter for Illumina extensions (_S1, _S2...)
            illumina_ext = re.compile(r'^(\w+)_S\d{1,2}$')
            illumina_match = illumina_ext.search(sample)
            if illumina_match:
                sample = illumina_match.group(1)
            print("Sample: " + sample)
            print("Metrics file: " + i)
            key = 1
            # we read and each line is converted into a dictionnary
            # 2 dicts, one per sample for global metrics and one more detailed for per region metrics
            # Each dict type (region, sample) is in fact composed of up to 3 dicts, one for autosomes, one for chrX, and one for chrY
            # coz more simple to handle when writing in separate tabs in excel
            with open(Path + i, 'r') as csvfile:
                csvreader = csv.reader(csvfile, delimiter =delim)
                for line in csvreader:
                    line = [w.replace(',', '.') for w in line]
                    expression = r'^.*#.*$'  # remove header
                    # autosomes
                    if not re.match(expression, line[0]) and line[0] != "chrX" and line[0] != 'X' and line[0] != "chrY" and line[0] != "Y":
                        (region_number, per_region_metrics, per_sample_metrics, chr_semaph) = build_dict(sample, line, region_number, per_region_metrics, per_sample_metrics, key, chr_semaph)
                        key += 1
                    # chrX
                    elif not re.match(expression, line[0]) and (line[0] == "chrX" or line[0] == "X"):
                        (region_number_ChrX, per_region_metrics_ChrX, per_sample_metrics_ChrX, chr_semaph) = build_dict(sample, line, region_number_ChrX, per_region_metrics_ChrX, per_sample_metrics_ChrX, key, chr_semaph)
                        key += 1
                    # chrY
                    elif not re.match(expression, line[0]) and (line[0] == "chrY" or line[0] == "Y"):
                        (region_number_ChrY, per_region_metrics_ChrY, per_sample_metrics_ChrY, chr_semaph) = build_dict(sample, line, region_number_ChrY, per_region_metrics_ChrY, per_sample_metrics_ChrY, key, chr_semaph)
                        key += 1

            # we build vcf dict if possible chr-pos => {sample, status (1)}
            if VcfDir is not False:
                chr_reg = r'^chr'
                # vcf_regexp = re.compile(r'.+\.vcf\.?g?z?$')
                # vcf_regexp = re.compile(r'.+\.vcf$')
                # vcf_gz_regexp = re.compile(r'.+\.vcf\.gz$')
                for vcf_file in vcf_list:
                    found_sample = False
                    # match_vcf = vcf_regexp.search(os.path.basename(vcf_file))
                    # if re.match(vcf_regexp, os.path.basename(vcf_file)):
                    if vcf_file.lower().endswith('.vcf'):
                        # if match_vcf:
                        # print("Associated VCF: " + vcf_file)
                        # vcf_reader = vcf.Reader(open(VcfDir + vcf_file, 'rb'))
                        # b opens in binary mode - to be modified
                        vcf_reader = vcf.Reader(open(VcfDir + vcf_file, 'r'))
                    elif vcf_file.lower().endswith('.vcf.gz'):
                        # elif re.match(vcf_gz_regexp, os.path.basename(vcf_file)):
                        vcf_reader = vcf.Reader(open(VcfDir + vcf_file, 'rb'))
                    else:
                        continue
                    test_record = next(vcf_reader)
                    # if test_record.genotype(sample):
                    for vcf_calls in test_record.samples:
                        # print(vcf_calls.sample)
                        if vcf_calls.sample == sample:
                            # if sample in test_record.samples.sample:
                            # ok our sample is here  we can read the entire vcf
                            print("Associated VCF: " + vcf_file)
                            vcf_chr_semaph = False
                            for record in vcf_reader:
                                # we store only heterozygous calls
                                sample_call = record.genotype(sample)
                                # if sample_call.is_het == True or (record.CHROM == 'chrX' or record.CHROM == 'X'): and record.FILTER == 'PASS'
                                # FILTER PASS returns empty record.FILTER
                                if sample_call.is_het is True and not record.FILTER:
                                    chrom = record.CHROM
                                    if not re.match(chr_reg, chrom):
                                        vcf_chr_semaph = True
                                    if vcf_chr_semaph:
                                        chrom = "chr" + chrom
                                    position = (chrom, record.POS)
                                    if position not in variants:
                                        # variants[position] = {sample: 1}
                                        variants[position] = {sample: {"AB": 1}}
                                    else:
                                        variants[position].update({sample: {"AB": 1}})
                                    if hasattr(sample_call, 'AD') and sample_call['AD'][0] > 0 and sample_call['AD'][1] > 0:
                                        variant_ab = round(float(sample_call['AD'][1]) / (float(sample_call['AD'][0]) + float(sample_call['AD'][1])), 3)
                                        variants[position].update({sample: {"AB": variant_ab}})
                                    # print(record.CHROM + "-" + str(record.POS) + "-" + str(sample) + "-" + str(record.heterozygosity) + "-" + str(sample_call['AD'][0]) + "-" + str(sample_call['AD']))
                            found_sample = True
                            break
                    if found_sample:
                        break
            print("")
    #############

    #############
    # pprint.pprint(variants)
    # sys.exit()



    #############
    # we iterate sample dict to compute the global Doc
    for sample_name in per_sample_metrics:
            per_sample_metrics[sample_name]["meanRawDoc"] = per_sample_metrics[sample_name]["rawDocSum"]/ region_number
    # FOR X and in addition we compute on the fly the x ratio which determinates the gender
    # then the ratio and predicted gender are added to the per_sample_metrics_chrX dict ONLY
    # FOR Y we look for Y regions in female samples and print a warning if found unconsistant
    if region_number_ChrX > 0:
        for sample_name in per_sample_metrics_ChrX:
            per_sample_metrics_ChrX[sample_name]["meanRawDoc"] = round(per_sample_metrics_ChrX[sample_name]["rawDocSum"]/ region_number_ChrX, 3)
            print(sample_name)
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

    (per_sample_metrics, per_region_metrics) = compute_ratio(per_sample_metrics, per_region_metrics, region_number, VcfDir, variants, 'Autosomes', het_high, het_low, hom_high, hom_low)
    if region_number_ChrX > 0:
        (per_sample_metrics_ChrX, per_region_metrics_ChrX) = compute_ratio(per_sample_metrics_ChrX, per_region_metrics_ChrX, region_number_ChrX, VcfDir, variants, 'ChrX', het_high, het_low, hom_high, hom_low)
    if region_number_ChrY > 0:
        (per_sample_metrics_ChrY, per_region_metrics_ChrY) = compute_ratio(per_sample_metrics_ChrY, per_region_metrics_ChrY, region_number_ChrY, VcfDir, variants, 'ChrY', het_high, het_low, hom_high, hom_low)

    # pp.pprint(per_region_metrics)
    # sys.exit()
    #############

    #############
    # We build a small list of genes of interest
    panel_list = []
    if Panel:
        panel = open(Panel, 'r')
        # panelreader = csv.DictReader(panel, delimiter='\t')
        for gene in panel:
            panel_list.append(gene.rstrip())
    ###########

    ###########
    # Excel conversion

    workbook = xlsxwriter.Workbook(OutFile)
    # Define style of cells
    style1 = workbook.add_format({'bold': True, 'bg_color': '#FF3333', 'locked': True})
    style2 = workbook.add_format({'bold': True, 'bg_color': '#FFC25E', 'locked': True})
    style3 = workbook.add_format({'bold': True, 'bg_color': '#5EBBFF', 'locked': True})
    style4 = workbook.add_format({'bold': True, 'bg_color': '#8F5EFF', 'locked': True})
    style5 = workbook.add_format({'bold': True, 'locked': True})
    style6 = workbook.add_format({'bold': True, 'num_format': 1, 'locked': True})
    style7 = workbook.add_format({'bold': True, 'color': '#F791E7', 'locked': True})
    style8 = workbook.add_format({'bold': True, 'color': '#5EBBFF', 'locked': True})
    styles = (style1, style2, style3, style4, style5, style6, style7, style8)
    # http://xlsxwriter.readthedocs.io/example_conditional_format.html#ex-cond-format
    # Add a format. Light red fill with dark red text.
    format1 = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
    # Add a format. Green fill with dark green text.
    format2 = workbook.add_format({'bg_color': '#C6EFCE', 'font_color': '#006100'})

    last_col = number_of_file*6 + 5
    # sample number x number of data blocks + 6 first cols to show - 1 coz numbering begins at 0 (A = 0)
    last_col_2_hide = number_of_file*5 + 5
    # sample number x number of data blocks to hide + 6 first cols to show - 1 coz numbering begins at 0 (A = 0)

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

    # pp.pprint(summary_regions)

    print("Summary worksheet...")
    print_worksheet('Summary', last_col, last_col_2_hide, workbook, summary_regions, 'summary', '', '', '', '', '', per_sample_metrics_ChrX, styles, number_of_file)
    if len(low_cov_summary_regions) > 0:
        print("Low Coverage Summary Worksheet...")
        print_worksheet('LowCovRegions', last_col, last_col_2_hide, workbook, low_cov_summary_regions, 'low_cov_summary', '', '', '', '', '', per_sample_metrics_ChrX, styles, number_of_file)
    if Panel:
        print("Panel worksheet...")
        print_worksheet('Panel', last_col, last_col_2_hide, workbook, panel_regions, 'panel', '', '', '', '', '', per_sample_metrics_ChrX, styles, number_of_file)
    # pp.pprint(per_sample_metrics_ChrX)

    workbook.close()
    sys.exit()


if __name__ == "__main__":
    main()
