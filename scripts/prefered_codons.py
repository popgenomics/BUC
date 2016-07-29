#!/usr/bin/python
import os
import sys
from Bio.SeqIO import parse


#cdsFile = sys.argv[1]
cdsFile = "Caenorhabditis_brenneri#solo200.cds"


# codon table
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


# translaTable
translaTable = {}
for i in codonTable:
	for j in codonTable[i]:
		translaTable[j] = i


# detect the prefered synonymous codons
command = "../scripts/prefered_codons.R"
testCommand = os.system(command)


codons = {} # contains preferred and unpreferred codons for aa with found preferred codons


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


# keptAA: list of retained amino acyls: with only 2 synonymous codons + significant preferred codon (based on correlation)
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


def codonAnalyse(x): # x = dataset[loci_i]
	nPos = len(x[0])
	toRemove = nPos%3
	if toRemove==0:
		toRemove = 3 # remove the last codon
	for i in range(0, nPos - toRemove, 3):
		codon = []
		for j in range(x['nSeq']+1):
			codonTMP = str(x[j][i:(i+3)])
			if "N" not in codonTMP:
				codon.append(codonTMP)
	return(codon)


