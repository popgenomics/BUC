#!/usr/bin/python
#./prefered_codons.py Caenorhabditis_brenneri#solo200.cds
import os
import sys
from numpy.random import choice
from Bio.SeqIO import parse
import numpy as np
import matplotlib.pyplot as plt
from pylab import savefig

cdsFile = sys.argv[1]
#cdsFile = "Caenorhabditis_brenneri#solo200.cds"

def codonAnalyse(x, nInd, threshold): # x = dataset[loci_i], nInd = twice the number of sampled individuals, threshold = value in [0 - nInd[
	# function that loop over codons of one alignment of sequences
	# watch for positions with polymorphism of two synonymous codons
	# object to return:
	res = {}
	res['nCodons'] = 0 # number of polymorphic codon in the studied locus
	res['pos'] = [] # vector with positions of the polymorphic codons
	res['pref'] = [] # vector with codes of the preferred codons
	res['unpref'] = [] # vector with codes of the unpreferred codons
	res['nPref'] = [] # vector with the number of individuals with the preferred variant at this position
	res['nUnpref'] = [] # vector with the number of individuals with the unPreferred variant at this position
	res['ENcP'] = [] # vector with the ENc' value of the locus
	res['expression'] = [] # vector with the expression value of the locus
	# retain nInd - threshold sequences
	nPos = len(x[0])
	toRemove = nPos%3
	if toRemove==0:
		toRemove = 3 # remove the last codon
	cnt = -1
	for i in range(0, nPos - toRemove, 3):
		cnt += 1
		codon = []
		for j in range(x['nSeq'] + 1):
			codonTMP = str(x[j][i:(i+3)])
			if "N" not in codonTMP:
				codon.append(codonTMP)
		if len(codon) > (nInd - threshold):
			sample = choice(len(codon), nInd-threshold, replace=False)
			codon2 = []
			for j in sample:
				codon2.append(codon[j])
			content = list(set(codon2))
			if(len(content) != 2):
				# reject everything which is not bi-allelic polymorphism
				continue
			if(len(content) == 2):
				# check that the polymorphism is synonymous
				tri1 = content[0]
				tri2 = content[1]
				if tri1 not in translaTable or tri2 not in translaTable:
					continue
				aa1 = translaTable[tri1]
				aa2 = translaTable[tri2]
				if(aa1 != aa2):
					# start the next step of the loop if the polymorphism is not synonymous
					continue
				else: # if synonymous polymorphism
					if aa1 not in keptAA: # if the codon is not in the list of retained aa
						continue
					else:
						prefCodon = codons[aa1]['prefCodon']
						unprefCodon = codons[aa1]['unprefCodon']
						nPref = codon2.count(prefCodon)
						nUnpref = codon2.count(unprefCodon)
						res['nCodons'] += 1
						res['pos'].append(cnt)
						res['pref'].append(prefCodon)
						res['unpref'].append(unprefCodon)
						res['nPref'].append(nPref) 
						res['nUnpref'].append(nUnpref)
						res['ENcP'].append(x['Ncp'])
						res['expression'].append(x['expression'])
	return(res)


# codon table: provides synonymous codon for an amino acyl 'i'
codonTable = {}
codonTable["A"] = ("GCA", "GCT", "GCG", "GCC")
codonTable["R1"] = ("CGA","CGT","CGG","CGC")
codonTable["R2"] = ("AGA", "AGG")
codonTable["N"] = ("AAT", "AAC")
codonTable["D"] = ("GAT", "GAC")
codonTable["C"] = ("TGT", "TGC")
codonTable["Q"] = ("CAA", "CAG")
codonTable["E"] = ("GAA", "GAG")
codonTable["G"] = ("GGA", "GGT", "GGG", "GGC")
codonTable["H"] = ("CAT", "CAC")
codonTable["I"] = ("ATT", "ATC", "ATA")
codonTable["L1"] = ("CTA", "CTT", "CTG", "CTC")
codonTable["L2"] = ("TTA", "TTG")
codonTable["K"] = ("AAA", "AAG")
codonTable["F"] = ("TTT", "TTC")
codonTable["P"] = ("CCA", "CCT", "CCG", "CCC")
codonTable["S1"] = ("TCA", "TCT", "TCG", "TCC")
codonTable["S2"] = ("AGT", "AGC")
codonTable["T"] = ("ACA", "ACT", "ACG", "ACC")
codonTable["Y"] = ("TAT", "TAC")
codonTable["V"] = ("GTA", "GTT", "GTG", "GTC")


# translaTable: provides amino acyl for a codon 'i'
translaTable = {}
for i in codonTable:
	for j in codonTable[i]:
		translaTable[j] = i


# detect the prefered synonymous codons
command = "../scripts/prefered_codons.R"
testCommand = os.system(command) # launch the R code "prefered_codons.R" producing "output_prefered_codons.txt"


codons = {} # contains preferred and unpreferred codons for aa with found preferred codons


# record the R output into python
infile = "output_prefered_codons.txt"
input = open(infile, "r")
for i in input:
	if "bestCodon" not in i:
		i = i.strip().split(" ")
		aa = i[0]
		prefCodon = i[1]
		unprefCodon = i[2]
		coefficient = float(i[3])
		pvalue = float(i[4])
		if coefficient > 0 and pvalue < 0.001:
			codons[aa] = {}
			codons[aa]["prefCodon"] = prefCodon
			codons[aa]["unprefCodon"] = unprefCodon
			codons[aa]["coefficient"] = coefficient
			codons[aa]["pvalue"] = pvalue
input.close()


# keptAA: list of retained amino acyls: with only 2 synonymous codons AND IF there is a significant preferred codon (based on correlation)
keptAA = []
for i in codons:
	if codons[i]['prefCodon'] != "NNN" and codons[i]['unprefCodon'] != "NNN":
		keptAA.append(i)


# parse the cds file:
input = parse(cdsFile, "fasta")


dataset = {}


for i in input:
	tmp = i.id.split("|")
	contig = tmp[0]
	species = tmp[1]
	individual = tmp[2]
	allele = tmp[3]
	if contig not in dataset:
		dataset[contig] = {}
		dataset[contig]["nSeq"] = -1
	dataset[contig]["nSeq"] += 1
	nSeq = dataset[contig]["nSeq"]
	dataset[contig][nSeq] = i.seq
input.close()


nInd = []
for i in dataset:
	nInd.append(dataset[i]["nSeq"])
nInd = max(nInd) + 1

if nInd < 6:
	threshold = 0
if nInd == 6:
	threshold = 1
if nInd == 7:
	threshold = 2
if nInd == 8:
	threshold = 3
if nInd == 9:
	threshold = 4
if nInd == 10:
	threshold = 5
if nInd >= 11:
	threshold = 6


# get the estimated expression
infile = "output_summarized.txt"
input = open(infile, "r")
i = input.readline()
i = i.strip().split("\t")
pos_Expression = i.index("GCcorrected_expression")
pos_Ncp = i.index("Ncp")
pos_loci = i.index("loci")


for i in input:
	i = i.strip().split("\t")
	loci = i[pos_loci]
	expression = float(i[pos_Expression])
	encp = float(i[pos_Ncp])
	if loci in dataset:
		dataset[loci]["expression"] = expression
		dataset[loci]["Ncp"] = encp
input.close()


# loop to analyse spectrum over the whole dataset
res = {}
res['loci'] = []
res['pos'] = []
res['pref'] = []
res['unpref'] = []
res['nPref'] = []
res['nUnpref'] = []
res['ENcP'] = []
res['expression'] = []


nPolymorphicPosition = 0
for i in dataset:
	if 'Ncp' not in dataset[i] or 'expression' not in dataset[i]:
		continue
	a = codonAnalyse(dataset[i], nInd, threshold)
	if a['nCodons'] > 0:
		for j in range(a['nCodons']):
			nPolymorphicPosition += 1
			res['loci'].append(i)
			res['pos'].append(a['pos'][j])
			res['pref'].append(a['pref'][j])
			res['unpref'].append(a['unpref'][j])
			res['nPref'].append(a['nPref'][j])
			res['nUnpref'].append(a['nUnpref'][j])
			res['ENcP'].append(a['ENcP'][j])
			res['expression'].append(a['expression'][j])


# output:
output = "loci\tpos\texpression\tEncp\tpref\tunpref\tnPref\tnUnpref\n"
for i in range(nPolymorphicPosition):
	output += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(res['loci'][i], res['pos'][i], res['expression'][i], res['ENcP'][i], res['pref'][i], res['unpref'][i], res['nPref'][i], res['nUnpref'][i])

outfile = open("output_pref_unpref_polymorphism.txt", "w")
outfile.write(output)
outfile.close()

if len(res['loci']) < 100:
	sys.exit(0)

# plot
N = max(res['nPref'] + res['nUnpref'])
ind = np.arange(N)  # the x locations for the groups
width = 0.25       # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)

table_nPref = []
table_nUnpref = []
for i in ind+1:
	table_nPref.append(res['nPref'].count(i))
	table_nUnpref.append(res['nUnpref'].count(i))

rects1 = ax.bar(ind, table_nUnpref, width, color='r')
rects2 = ax.bar(ind+width, table_nPref, width, color='g')

ax.set_xlim(left = -0.5)
ax.set_xlim(right = N)
ax.set_xticks(ind+width)
ax.set_xticklabels(ind+1)

savefig("codon_frequency_spectrum.pdf", bbox_inches='tight')

