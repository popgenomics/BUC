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

utr3File = sys.argv[1]
utr5File = sys.argv[2]

#utr3File = "Trachemys_scripta#solo200.utr3"
#utr5File = "Trachemys_scripta#solo200.utr5"

utr = {}
utr3 = {}
utr5 = {}

# deals with utr3
input = parse(utr3File, "fasta")
for i in input:
	gene = i.id.split("|")[0]
	if gene not in utr:
		utr[gene] = ""
		utr3[gene] = ""
		utr5[gene] = ""
	utr[gene] += str(i.seq)
	utr3[gene] += str(i.seq)
input.close()

# deals with utr5
input = parse(utr5File, "fasta")
for i in input:
	gene = i.id.split("|")[0]
	if gene not in utr:
		utr[gene] = ""
		utr3[gene] = ""
		utr5[gene] = ""
	utr[gene] += str(i.seq)
	utr5[gene] += str(i.seq)
input.close()

res = "loci\tLutr\tGCutr\tLutr3\tGCutr3\tLutr5\tGCutr5\n"
for i in utr:
	tmp = countGC(utr[i])
	tmp3 = countGC(utr3[i])
	tmp5 = countGC(utr5[i])
	res = res + "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(i, tmp["L"], tmp["GC"], tmp3["L"], tmp3["GC"], tmp5["L"], tmp5["GC"])

outfile = "output_GC_UTR.txt"
output = open(outfile, "w")
output.write(res)
output.close()

