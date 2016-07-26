#!/usr/bin/python
from scipy import stats
from scipy import nan
from scipy import mean
from scipy import nanmean
from scipy import nanstd
import sys

dataset = sys.argv[1]
NCBI = sys.argv[2]
NCBI = " ".join(NCBI.split("_"))

def correlation(x,y):
	a, b = [], []
	for i in range(len(x)):
		if x[i]!=nan and y[i]!=nan:
			a.append(x[i])
			b.append(y[i])
	r, p = stats.pearsonr(a, b)
	rho, p = stats.spearmanr(a, b)
	tau, p = stats.kendalltau(a, b)
	return([round(r, 5), round(rho, 5), round(tau, 5)])
	
def stat(x):
	return([round(nanmean(x), 5), round(nanstd(x), 5)])


orfFile = "output_summarized.txt"

orf = {}
input = open(orfFile, "r")
i = input.readline()
header = i.strip().split()

for i in header:
	orf[i] = []

for i in input:
	i = i.strip().split("\t")
	for j in range(len(header)):
		if j==0:
			orf[header[j]].append(i[j])
		if j in range(1, len(header)+1, 1):
			if i[j] == "nan" or "*" in i[j]:
				orf[header[j]].append(nan)
			else:
				orf[header[j]].append(float(i[j]))
input.close()

nLoci = len(orf["loci"])
Ltot = stat(orf["Ltot"])
GC = stat(orf["GCtot"])
GC3 = stat(orf["GC3"])
GC12 = stat(orf["GC12"])
GCutr = stat(orf["GCutr"])


res = "dataset\tNCBIname\tnLoci\t"
res2 = dataset + "\t" + NCBI + "\t" + str(nLoci) + "\t"

for i in ["Ltot", "GCtot", "GC12", "GC3", "GCutr"]:
	res2 += "\t".join([str(k) for k in stat(orf[i])]) + "\t"
	for j in ["avg", "std"]:
		res = res + "{0}\t".format(i + "_" + j)


for i in ["GC12", "GC3", "GCutr", "Nc_JN", "Ncp_JN", "ENc_cubt", "ENcP_cubt", "ENcFDS_cubt", "ENcPFDS_cubt", "Um_cubt", "Ltot", "nReads"]:
	for j in ["r", "rho", "tau"]:
		res += "{0}_CorrectedExpression_{1}\t".format(j, "".join(i.split("_")))
	res2 += "\t".join([ str(k) for k in correlation(orf["GCcorrected_expression"], orf[i])]) + "\t"


for i in ["GC12", "GC3", "GCutr", "Nc_JN", "Ncp_JN", "ENc_cubt", "ENcP_cubt", "ENcFDS_cubt", "ENcPFDS_cubt", "Um_cubt", "Ltot"]:
	for j in ["r", "rho", "tau"]:
		res += "{0}_nCountedReads_{1}\t".format(j, "".join(i.split("_")))
	res2 += "\t".join([ str(k) for k in correlation(orf["nReads"], orf[i])]) + "\t"

res += "r_GC3_GCutr\trho_GC3_GCutr\ttau_GC3_GCutr\n"
res2 += "\t".join([ str(k) for k in correlation(orf["GC3"], orf["GCutr"])]) + "\n"

outfile = open("output_final_1.txt", "w")
outfile.write(res + res2)
outfile.close()

