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
#import xlwt
import xlsxwriter
#import xlrd

# ==============================================================================

print("\n","||||||||||||||||||||||||||||||||||||||||||||||||||||")
print("\n","||||||||||||||||||||||||||||||||||||||||||||||||||||", "\n")
pp = pprint.PrettyPrinter(indent=4, depth=6)

#added david 27/03/2018 - deal with path passed as arg
parser = argparse.ArgumentParser(description='deal with path.')
parser.add_argument('-p', '--path', default='.')
args = parser.parse_args()
#print(args.path)
#sys.exit()

#Path = "/Users/david/Downloads/henri_cnv/03_2018/170925H27JNL/"
Path = args.path


filelist = os.listdir(Path)
number_of_file = 0
sample_name = []
counter = 0
dict_regions = {}
dict_patients = {}
moyenne_totale = 0
dict_regions_ChrX = {}
dict_patients_ChrX = {}
moyenne_totale_ChrX = 0
counter_ChrX = 0

#############
for i in filelist:
	if i.endswith("coverage.tsv"):  # You could also add "and i.startswith('f')
		number_of_file += 1
		regex = re.compile(r'(.*)[\._]coverage\.tsv')
		matchObj = regex.search(os.path.basename(i))
		if(matchObj):
			baseName = matchObj.group(1)
			sample_name.append(baseName)
		else:
			sys.exit("Error [1] : file \'" + i + "\' does not respect format (sample.coverage.csv).")
		with open(Path + i, 'r') as csvfile:
			csvreader = csv.reader(csvfile, delimiter ='\t')
			for line in csvreader:
				line = [w.replace(',', '.') for w in line]
				expression = r'^.*#.*$' # remove header
				if not re.match(expression, line[0]) and line[0] != "chrX" and line[0] != "chrY":
					#création moyenne par patient
					if(baseName not in dict_patients):
						dict_patients[baseName] = {"sum_patient": float(line[4])}
					else:
						dict_patients[baseName]["sum_patient"] += float(line[4])
					# création dictionnaire avec les valeurs brutes
					coordinate = (line[0],line[1],line[2], line[3])
					if coordinate not in dict_regions :
						counter += 1
						dict_regions[coordinate] = {baseName : {"Brute" : float(line[4])}}
					else :
						dict_regions[coordinate][baseName] = {"Brute" : float(line[4])}
				else :
					if not re.match(expression, line[0]) and line[0] == "chrX":
						if(baseName not in dict_patients_ChrX):
							dict_patients_ChrX[baseName] = {"sum_patient": float(line[4])}
						else:
							dict_patients_ChrX[baseName]["sum_patient"] += float(line[4])
						# création dictionnaire avec les valeurs brutes
						coordinate = (line[0],line[1],line[2],line[3])
						if coordinate not in dict_regions_ChrX :
							counter_ChrX += 1
							dict_regions_ChrX[coordinate] = {baseName : {"Brute" : float(line[4])}}
						else :
							dict_regions_ChrX[coordinate][baseName] = {"Brute" : float(line[4])}

#############

# Itération sur le dictionnaire régions pour calculer la moyenne totale
sum_of_brute = 0
for coordinate in dict_regions :
	#print (coordinate)
	for sample_other in dict_regions[coordinate]:
		sum_of_brute += dict_regions[coordinate][sample_other]['Brute']

moyenne_totale = sum_of_brute/(counter*number_of_file)

# FOR X

sum_of_brute_X = 0
for coordinate in dict_regions_ChrX :
	for sample_other in dict_regions_ChrX[coordinate]:
		sum_of_brute += dict_regions_ChrX[coordinate][sample_other]['Brute']

moyenne_totale_ChrX = sum_of_brute/(counter_ChrX*number_of_file)

#############

# Itération sur le dictionnaire patient pour calculer la moyenne par par patient
for sample_name in dict_patients:
		dict_patients[sample_name]["moyenne_patient"] = dict_patients[sample_name]["sum_patient"]/ counter
# FOR X
for sample_name in dict_patients_ChrX:
		dict_patients_ChrX[sample_name]["moyenne_patient"] = dict_patients_ChrX[sample_name]["sum_patient"]/ counter_ChrX

#############
# Itération sur le dictionnaire patient pour calculer la moyenne par exon et l'exon normalisé
for coordinate in dict_regions :
	for sample_name in dict_regions[coordinate]:
		total = 0
		sample_number = 0
		for sample_other in dict_regions[coordinate]:
			if sample_other != sample_name:
				total += dict_regions[coordinate][sample_other]['Brute']
				sample_number +=1
				moyenne_exon = total / sample_number
				dict_regions[coordinate][sample_name]["moyenne_exon"] = float(moyenne_exon)

# FOR X
for coordinate in dict_regions_ChrX :
	for sample_name in dict_regions_ChrX[coordinate]:
		total = 0
		sample_number = 0
		for sample_other in dict_regions_ChrX[coordinate]:
			if sample_other != sample_name:
				total += dict_regions_ChrX[coordinate][sample_other]['Brute']
				sample_number +=1
				moyenne_exon = total / sample_number
				dict_regions_ChrX[coordinate][sample_name]["moyenne_exon"] = float(moyenne_exon)

#############
# Itération sur le dictionnaire régions pour calculer la somme moyenne par patient
# sauf pour les exons du patient - Somme ajoutée dans le dict patients
for sample_name in dict_patients:
	dict_patients [sample_name]["somme_moyenne_exon"] = 0
	for coordinate in dict_regions :
		dict_patients [sample_name]["somme_moyenne_exon"] += dict_regions[coordinate][sample_name]["moyenne_exon"]

for sample_name in dict_patients:
	dict_patients [sample_name]["moyenne_totale_sans_le_patient"] = dict_patients [sample_name]["somme_moyenne_exon"] / counter

# FOR X
for sample_name in dict_patients_ChrX:
	dict_patients_ChrX [sample_name]["somme_moyenne_exon"] = 0
	for coordinate in dict_regions_ChrX :
		dict_patients_ChrX [sample_name]["somme_moyenne_exon"] += dict_regions_ChrX[coordinate][sample_name]["moyenne_exon"]

for sample_name in dict_patients_ChrX:
	dict_patients_ChrX [sample_name]["moyenne_totale_sans_le_patient"] = dict_patients_ChrX [sample_name]["somme_moyenne_exon"] / counter_ChrX

#############
# Itération sur le dictionnaire régions pour ajouter la valeur "exon_normalise"
# Moyenne exon / Moyenne tota
#@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@
for sample_name in dict_patients:
	for coordinate in dict_regions :
		dict_regions [coordinate][sample_name]["exon_normalise"] = dict_regions[coordinate][sample_name]["moyenne_exon"] / dict_patients[sample_name]["moyenne_totale_sans_le_patient"]

# FOR X

for sample_name in dict_patients_ChrX:
	for coordinate in dict_regions_ChrX :
		dict_regions_ChrX [coordinate][sample_name]["exon_normalise"] = dict_regions_ChrX[coordinate][sample_name]["moyenne_exon"] / dict_patients_ChrX[sample_name]["moyenne_totale_sans_le_patient"]



#############
# Normalisation par patient : valeur couverture exon / moyenne patient

for coordinate in dict_regions :
	for sample_name in dict_regions[coordinate]:
		patient_normalise = dict_regions[coordinate][sample_name]['Brute'] / dict_patients[sample_name]["moyenne_patient"]
		dict_regions[coordinate][sample_name]["patient_normalise"] = float(patient_normalise)

# FOR X

for coordinate in dict_regions_ChrX :
	for sample_name in dict_regions_ChrX[coordinate]:
		patient_normalise = dict_regions_ChrX[coordinate][sample_name]['Brute'] / dict_patients_ChrX[sample_name]["moyenne_patient"]
		dict_regions_ChrX[coordinate][sample_name]["patient_normalise"] = float(patient_normalise)

#############

# Calcul ratio normalise

for coordinate in dict_regions :
	for sample_name in dict_regions[coordinate]:
		try :
			ratio_normalise = float(dict_regions[coordinate][sample_name]["patient_normalise"]) / float(dict_regions[coordinate][sample_name]["exon_normalise"])
		except ZeroDivisionError :
			ratio_normalise = float(0)
		dict_regions[coordinate][sample_name]["ratio_normalise"] = float(ratio_normalise)
# FOR X
for coordinate in dict_regions_ChrX :
	for sample_name in dict_regions_ChrX[coordinate]:
		try :
			ratio_normalise = float(dict_regions_ChrX[coordinate][sample_name]["patient_normalise"]) / float(dict_regions_ChrX[coordinate][sample_name]["exon_normalise"])
		except ZeroDivisionError :
			ratio_normalise = float(0)
		dict_regions_ChrX[coordinate][sample_name]["ratio_normalise"] = float(ratio_normalise)

pp.pprint(dict_patients)
#############
# Rules (à vraiment commencer et finir)
# for coordinate in dict_regions :
# 	for sample_name in dict_regions[coordinate]:
# 		dict_regions[coordinate][sample_name]["interpretation"] = "no CNV"
# 		if  0.3 <= dict_regions[coordinate][sample_name]["ratio_normalise"] <= 0.7 and dict_regions[coordinate][sample_name]["Brute"] > 300:
# 			dict_regions[coordinate][sample_name]["interpretation"] = "del htz"
# 		if  1.3 <= dict_regions[coordinate][sample_name]["ratio_normalise"] <= 1.7 and dict_regions[coordinate][sample_name]["Brute"] > 300:
# 			dict_regions[coordinate][sample_name]["interpretation"] = "dup htz"
# 		if  1.8 <= dict_regions[coordinate][sample_name]["ratio_normalise"] <= 2.2 and dict_regions[coordinate][sample_name]["Brute"] > 300:
# 			dict_regions[coordinate][sample_name]["interpretation"] = "del htz"



# pp.pprint(dict_regions)

# Writing output file

with open('cnv_analysis.txt', 'w') as csv_file:
	writer = csv.writer(csv_file)
	header = "Chr\tStart\tEnd\tRegionID\t"
	for coordinate in dict_regions :
		for sample_name in dict_regions[coordinate]:
			header += str(sample_name) + "_brute" + "\t"
		for sample_name in dict_regions[coordinate]:
			header += str(sample_name) + "_moy_exon" + "\t"
		for sample_name in dict_regions[coordinate]:
			header += str(sample_name) + "_patient_normalise" + "\t"
		for sample_name in dict_regions[coordinate]:
			header += str(sample_name) + "_exon_normalise" + "\t"
		for sample_name in dict_regions[coordinate]:
			header += str(sample_name) + "_ratio" + "\t"
		break
	csv_file.write(header + "\n")
	for coordinate in dict_regions :
		csv_file.write(
			str(coordinate[0]) + "\t" +
			str(coordinate[1]) + "\t" +
			str(coordinate[2]) + "\t" +
			str(coordinate[3]) + "\t"
			)
		for sample_name in dict_regions[coordinate] :
			csv_file.write(str(dict_regions[coordinate][sample_name]["Brute"]) + "\t")
		for sample_name in dict_regions[coordinate] :
			csv_file.write(str(round(dict_regions[coordinate][sample_name]["moyenne_exon"], 2)) + "\t")
		for sample_name in dict_regions[coordinate] :
			csv_file.write(str(round(dict_regions[coordinate][sample_name]["patient_normalise"],2)) + "\t")
		for sample_name in dict_regions[coordinate] :
			csv_file.write(str(round(dict_regions[coordinate][sample_name]["exon_normalise"],2)) + "\t")
		for sample_name in dict_regions[coordinate] :
			csv_file.write(str(round(dict_regions[coordinate][sample_name]["ratio_normalise"],3)) + "\t")
		csv_file.write("\n")

with open('cnv_analysis_ChrX.txt', 'w') as csv_file:
	writer = csv.writer(csv_file)
	header = "Chr\tStart\tend\tRegionID\t"
	for coordinate in dict_regions_ChrX :
		for sample_name in dict_regions_ChrX[coordinate]:
			header += str(sample_name) + "_brute" + "\t"
		for sample_name in dict_regions_ChrX[coordinate]:
			header += str(sample_name) + "_moy_exon" + "\t"
		for sample_name in dict_regions_ChrX[coordinate]:
			header += str(sample_name) + "_patient_normalise" + "\t"
		for sample_name in dict_regions_ChrX[coordinate]:
			header += str(sample_name) + "_exon_normalise" + "\t"
		for sample_name in dict_regions_ChrX[coordinate]:
			header += str(sample_name) + "_ratio" + "\t"
		break
	csv_file.write(header + "\n")
	for coordinate in dict_regions_ChrX :
		csv_file.write(
			str(coordinate[0]) + "\t" +
			str(coordinate[1]) + "\t" +
			str(coordinate[2]) + "\t" +
			str(coordinate[3]) + "\t"
			)
		for sample_name in dict_regions_ChrX[coordinate] :
			csv_file.write(str(dict_regions_ChrX[coordinate][sample_name]["Brute"]) + "\t")
		for sample_name in dict_regions_ChrX[coordinate] :
			csv_file.write(str(round(dict_regions_ChrX[coordinate][sample_name]["moyenne_exon"], 2)) + "\t")
		for sample_name in dict_regions_ChrX[coordinate] :
			csv_file.write(str(round(dict_regions_ChrX[coordinate][sample_name]["patient_normalise"],2)) + "\t")
		for sample_name in dict_regions_ChrX[coordinate] :
			csv_file.write(str(round(dict_regions_ChrX[coordinate][sample_name]["exon_normalise"],2)) + "\t")
		for sample_name in dict_regions_ChrX[coordinate] :
			csv_file.write(str(round(dict_regions_ChrX[coordinate][sample_name]["ratio_normalise"],3)) + "\t")
		csv_file.write("\n")

# Sorting files, works on Mac only, if you are poor (like Charles) with other OS, go fuck yourself
os.system("sort -k1.4n -k2,2n -k3,3n cnv_analysis.txt > cnv_analysis_sorted.txt")
os.system("sort -k1.4n -k2,2n -k3,3n cnv_analysis_ChrX.txt > cnv_analysis_ChrX_sorted.txt")
#
### Excel convertion


workbook = xlsxwriter.Workbook('cnv_analysis_sorted.xlsx')
# Define style of cells
style1 = workbook.add_format({'bold': True, 'bg_color': 'coral'})
style2 = workbook.add_format({'bold': True, 'bg_color': 'red'})
style3 = workbook.add_format({'bold': True, 'bg_color': 'blue'})
style4 = workbook.add_format({'bold': True, 'bg_color': 'purple'})
style5 = workbook.add_format({'bold': True})
#structure data from txt
f = open('cnv_analysis_sorted.txt', 'r+')
row_list = []
for row in f:
    row_list.append(row.split('\t'))
column_list = zip(*row_list)
column_list2 = zip(*row_list)
#worksheet for autosomes
worksheet = workbook.add_worksheet('Autosomes')
worksheet.freeze_panes(1, 4)
worksheet.set_row(0, 20, style5)
worksheet.set_column('A:D', 20, style5)
#structure data from txt for X chr
g = open('cnv_analysis_ChrX_sorted.txt', 'r+')
row_list_x = []
for row in g:
    row_list_x.append(row.split('\t'))
column_list_x = zip(*row_list_x)
column_list_x2 = zip(*row_list_x)
#worksheet for X chr
worksheet_x = workbook.add_worksheet('Chromosome X')
worksheet_x.freeze_panes(1, 4)
worksheet_x.set_row(0, 20, style5)
worksheet_x.set_column('A:D', 20, style5)
#worksheet summary to get only interesting stuff
summary = workbook.add_worksheet('Summary')
summary.freeze_panes(1, 4)
summary.set_row(0, 20, style5)
summary.set_column('A:D', 20, style5)

i = 0
interesting = []
regex_ratio = re.compile(r'(.*)_ratio$')
regex_noCNV = re.compile(r'no CNV')
for column in column_list:
	for item in range(len(column)):
		if regex_ratio.search(column[0]):
			if(item > 0):#to exclude header
				if(float(column[item]) <= 0.3):
					worksheet.write(item, i , column[item],style1)
					for k in range(item-1,item+2):
						if (k > 0):
							interesting.append(k)
				elif(float(column[item]) <= 0.75):
					worksheet.write(item, i, column[item],style2)
					for k in range(item-1,item+2):
						if (k > 0):
							interesting.append(k)
				elif(float(column[item]) <= 1.3):
					worksheet.write(item, i, column[item], style5)
				elif(float(column[item]) <= 1.7):
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
			worksheet.write(item,i,column[item])
		if (item == 0):
			summary.write(item,i,column[item], style5)
	i+=1
	
#2nd iteration on the same data to print only interesting regions surrounding dels/dup
i = 0
j = 1
#remove double values
uniq_interesting = list(set(interesting))
for column in column_list2:
	j = 1
	for item in range(len(column)):
		if (item in uniq_interesting):
			if regex_ratio.search(column[0]):
				if(item > 0):
					if(float(column[item]) <= 0.3):
						summary.write(j, i , column[item],style1)
					elif(float(column[item]) <= 0.75):
						summary.write(j, i, column[item],style2)
					elif(float(column[item]) <= 1.3):
						summary.write(j, i, column[item], style5)
					elif(float(column[item]) <= 1.7):
						summary.write(j, i, column[item],style3)
					else :
						summary.write(j, i, column[item],style4)
				else:
					summary.write(j, i, column[item], style5)
			else:
				summary.write(j,i,column[item])
			j+=1
	i+=1
	
#same for X chr
interesting = []
i = 0
for column in column_list_x:
	for item in range(len(column)):
		if regex_ratio.search(column[0]):
			if(item > 0):
				if(float(column[item]) <= 0.3):
					worksheet_x.write(item, i, column[item],style1)
					for k in range(item-1,item+2):
						if (k > 0):
							interesting.append(k)
				elif(float(column[item]) <= 0.75):
					worksheet_x.write(item, i, column[item],style2)
					for k in range(item-1,item+2):
						if (k > 0):
							interesting.append(k)
				elif(float(column[item]) <= 1.3):
					worksheet_x.write(item, i, column[item], style5)
				elif(float(column[item]) <= 1.7):
					worksheet_x.write(item, i, column[item],style3)
					for k in range(item-1,item+2):
						if (k > 0):
							interesting.append(k)
				else :
					worksheet_x.write(item, i, column[item],style4)
					for k in range(item-1,item+2):
						if (k > 0):
							interesting.append(k)
			else:
				worksheet_x.write(item, i, column[item], style5)
		else:
			worksheet_x.write(item, i, column[item])
	i+=1
	
uniq_interesting = list(set(interesting))
i = 0
l = j + 1  #keep current line number in memory with one empty line
for column in column_list_x2:
	l = j + 1
	for item in range(len(column)):
		if (item in uniq_interesting):
			if regex_ratio.search(column[0]):
				if(item > 0):
					if(float(column[item]) <= 0.3):
						summary.write(l, i , column[item],style1)
					elif(float(column[item]) <= 0.75):
						summary.write(l, i, column[item],style2)
					elif(float(column[item]) <= 1.3):
						summary.write(l, i, column[item], style5)
					elif(float(column[item]) <= 1.7):
						summary.write(l, i, column[item],style3)
					else :
						summary.write(l, i, column[item],style4)
				else:
					summary.write(l, i, column[item], style5)
			else:
				summary.write(l,i,column[item])
			l+=1
	i+=1
summary.set_column('E:BL', None, None, {'level': 1, 'hidden': True})
worksheet.set_column('E:BL', None, None, {'level': 1, 'hidden': True})
worksheet_x.set_column('E:BL', None, None, {'level': 1, 'hidden': True})
workbook.close()
