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
	return([r, rho, tau])
	
def stat(x):
	return([nanmean(x), nanstd(x)])


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
GCutr = stat(orf["GCutr"])

#nReads_GC12 = correlation(orf["nReads"], orf["GC12"])
#GCcorrectedExpression_GC12 = correlation(orf["GCcorrected_expression"], orf["GC12"])
#GCcorrectedExpression_GC3 = correlation(orf["GCcorrected_expression"], orf["GC3"])
#GCcorrectedExpression_GCutr = correlation(orf["GCcorrected_expression"], orf["GCutr"])
#GCcorrectedExpression_ENc = correlation(orf["GCcorrected_expression"], orf["ENc_cubt"])
#GCcorrectedExpression_ENcP = correlation(orf["GCcorrected_expression"], orf["ENcP_cubt"])
#GCcorrectedExpression_ENcFDS = correlation(orf["GCcorrected_expression"], orf["ENcFDS_cubt"])
#GCcorrectedExpression_ENcPFDS = correlation(orf["GCcorrected_expression"], orf["ENcPFDS_cubt"])
#GCcorrectedExpression_Um = correlation(orf["GCcorrected_expression"], orf["Um_cubt"])
#GCcorrectedExpression_L = correlation(orf["GCcorrected_expression"], orf["Ltot"])

res = "dataset\tNCBIname\tnLoci\t"
res2 = dataset + "\t" + NCBI + "\t" + str(nLoci) + "\t"

for i in ["Ltot", "GCtot", "GC3", "GCutr"]:
	res2 += "\t".join([str(k) for k in stat(orf[i])]) + "\t"
	for j in ["avg", "std"]:
		res = res + "{0}\t".format(i + "_" + j)


for i in ["GC12", "GC3", "GCutr", "ENc_cubt", "ENc_cubt", "ENcP_cubt", "ENcFDS_cubt", "ENcPFDS_cubt", "Um_cubt", "Ltot", "nReads"]:
	for j in ["r", "rho", "tau"]:
		res += "{0}.expression_{1}\t".format(j, i.split("_")[0])
	res2 += "\t".join([ str(k) for k in correlation(orf["GCcorrected_expression"], orf[i])]) + "\t"


for i in ["GC12", "GC3", "GCutr", "ENc_cubt", "ENc_cubt", "ENcP_cubt", "ENcFDS_cubt", "ENcPFDS_cubt", "Um_cubt", "Ltot"]:
	for j in ["r", "rho", "tau"]:
		res += "{0}.nReads_{1}\t".format(j, i.split("_")[0])
	res2 += "\t".join([ str(k) for k in correlation(orf["nReads"], orf[i])]) + "\t"

res += "r.GC3_GCutr\trho.GC3_GCutr\ttau.GC3_GCutr\n"
res2 += "\t".join([ str(k) for k in correlation(orf["GC3"], orf["GCutr"])]) + "\n"

outfile = open("output_final_2.txt", "w")
outfile.write(res + res2)
outfile.close()

