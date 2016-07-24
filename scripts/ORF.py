#!/usr/bin/python
# ./UTR.py Trachemys_scripta#solo200.utr3 Trachemys_scripta#solo200.utr5 
from Bio.SeqIO import parse
from numpy import nan
import sys


def countGC(x):
	res = {}
	nA = x.upper().count("A")
	nT = x.upper().count("T")
	nC = x.upper().count("C")
	nG = x.upper().count("G")
	Ltot = nA + nT + nG + nC
	if Ltot > 20:
		GCpercentage = (nC + nG) / (Ltot*1.0)
		res["GC"] = GCpercentage
		res["L"] = Ltot
	else:
		res["GC"] = nan
		res["L"] = nan
	return(res)


orfFile = sys.argv[1]
covFile = sys.argv[2]
csvFile = sys.argv[3]

#orfFile = "Tripylus_abatoides#solo200.orf.new"  
#covFile = "Tripylus_abatoides#solo200.cov" 
#csvFile = "Tripylus_abatoides#solo200.csv"

# coverage
cov = {}
input = open(covFile, "r")
for i in input:
	i = i.strip().split("\t")
	cov[i[0]] = i[1]
input.close()

# csv: buc
csv = {}
input = open(csvFile, "r")
tmp = input.readline()
for i in input:
	i = i.strip().split("\t")
	gene = i[0][1::]
	csv[gene] = {}
	csv[gene]["GC"] = i[2]
	csv[gene]["GC1"] = i[3]
	csv[gene]["GC2"] = i[4]
	csv[gene]["GC3"] = i[5]
	csv[gene]["Nc"] = i[6]
	csv[gene]["NcP"] = i[7]
	csv[gene]["NcFDS"] = i[8]
	csv[gene]["NcPFDS"] = i[9]
	if i[10] == "****":
		i[10] = "nan"
	csv[gene]["Um"] = i[10]

orf = {}
# deals with ORF
input = parse(orfFile, "fasta")
for i in input:
	gene = i.id.split("|")[0]
	if gene not in orf:
		orf[gene] = ""
	orf[gene] += str(i.seq)
input.close()


res = "loci\tLtot\tGCtot\tL12\tGC12\tL3\tGC3\tnReads\tGC_cubt\tGC1_cubt\tGC2_cubt\tGC3_cubt\tENc_cubt\tENcP_cubt\tENcFDS_cubt\tENcPFDS_cubt\tUm_cubt\n"
for i in orf:
	tmp = countGC(orf[i])
	tmp12 = countGC(orf[i][0::3] + orf[i][1::3])
	tmp3 = countGC(orf[i][2::3])
	res = res + "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t".format(i, tmp["L"], tmp["GC"], tmp12["L"], tmp12["GC"], tmp3["L"], tmp3["GC"], cov[i])
	res = res + "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(csv[i]["GC"], csv[i]["GC1"], csv[i]["GC2"], csv[i]["GC3"], csv[i]["Nc"], csv[i]["NcP"], csv[i]["NcFDS"], csv[i]["NcPFDS"], csv[i]["Um"])

outfile = "output_GC_ORF.txt"
output = open(outfile, "w")
output.write(res)
output.close()

