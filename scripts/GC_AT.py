#!/usr/bin/python
# script that output the preferred codons for 4-fould aa
# between xyA versus xyT and
# between xyC versus xyG

from Bio.SeqIO import parse
from scipy import stats
from scipy import nan
from scipy import mean
from scipy import nanmean
from scipy import nanstd
import sys


seuil_pvalue = 0.05


codonTable = {}
codonTable["GC"] = {} 
codonTable["CG"] = {} 
codonTable["GG"] = {} 
codonTable["CT"] = {} 
codonTable["CC"] = {} 
codonTable["TC"] = {} 
codonTable["AC"] = {} 
codonTable["GT"] = {} 

bases = ["A", "T", "G", "C"]

for i in codonTable:
	for j in bases:
		codonTable[i][j] = {}
		codonTable[i][j]["expression"] = []
		codonTable[i][j]["occurence"] = []
		codonTable[i][j]["pos"] = -1
		

input = open("output_summarized.txt", "r")
header = input.readline().strip().split("\t")


for i in codonTable:
	for j in bases:
		codonTable[i][j]["pos"] = int(header.index(i+j))


expression = header.index("GCcorrected_expression")
longueur = header.index("Ltot")


for loci in input:
	loci = loci.strip().split("\t")
	for i in codonTable:
		for j in bases:
			codonTable[i][j]["expression"].append(float(loci[expression]))
			pos = codonTable[i][j]["pos"]
			codonTable[i][j]["occurence"].append(float(loci[pos])/(float(loci[longueur]*3)))


res = {}
for i in codonTable:
	res[i] = {}
	for j in bases:
		res[i][j] = {}
		r, p = stats.pearsonr(codonTable[i][j]["expression"], codonTable[i][j]["occurence"])
		res[i][j]["R"] = round(r, 5)
		res[i][j]["pvalue"] = p


print("preferred\tpearsonR\tpvalue\tunpreferred\tpearsonR\tpalue")
for i in res:
		if res[i]["A"]["R"] > res[i]["T"]["R"]:
			j = "A"
			k = "T"
		else:
			j = "T"
			k = "A"
		if res[i][j]["pvalue"] < seuil_pvalue and res[i][j]["R"]>0:
			print("{0}{1}\t{2}\t{3}\t{4}{5}\t{6}\t{7}".format(i, j, res[i][j]["R"], res[i][j]["pvalue"], i, k, res[i][k]["R"], res[i][k]["pvalue"]))
		if res[i]["G"]["R"] > res[i]["C"]["R"]:
			j = "G"
			k = "C"
		else:
			j = "C"
			k = "G"
		if res[i][j]["pvalue"] < seuil_pvalue and res[i][j]["R"]>0:
			print("{0}{1}\t{2}\t{3}\t{4}{5}\t{6}\t{7}".format(i, j, res[i][j]["R"], res[i][j]["pvalue"], i, k, res[i][k]["R"], res[i][k]["pvalue"]))


#print("\n")
#for i in res:
#	for j in bases:
#		print("{0}{1}\t{2}\t{3}".format(i, j, res[i][j]["R"], res[i][j]["pvalue"]))

