#!/usr/bin/python

infile = "output_summarized.txt"
input = open(infile, "r")


header = input.readline().strip().split("\t")


res = {}
nTot = 0
for i in range(36, 100, 1):
	res[header[i]] = {}
	res[header[i]]["count"] = 0


nLoci = 0
for i in input:
	nLoci += 1
	i = i.strip().split("\t")
	for j in range(36, 100, 1):
		res[header[j]]["count"] += int(i[j])
		nTot += int(i[j])
input.close()


for i in res:
	res[i]["freq"] = round(1000 * res[i]["count"]/(nTot*1.0), 1)


GC1 = 0
GC2 = 0
GC3 = 0

count1 = 0
count2 = 0
bases = ['T', 'C', 'A', 'G']

message = "{0} CDS's ({1} codons)\nfields: [triplet] [frequency: per thousand] ([number])\n\n".format(nLoci, nTot)

for i in bases:
	for j in bases:
		for k in bases:
			count1 += 1
			count2 += 1
			tmp = i + k + j
			message += "{0}\t{1} ({2})".format(tmp, res[tmp]["freq"], res[tmp]["count"])
			if count1 % 4 == 0:
				message += "\n"
			else:
				message += "\t"
			if count2 % 16 == 0:
				message += "\n"
			if i == "G" or i == "C":
				GC1 += res[tmp]["freq"] 
			if k == "G" or k == "C":
				GC2 += res[tmp]["freq"] 
			if j == "G" or j == "C":
				GC3 += res[tmp]["freq"] 

message += "\nCoding GC {0:.2f}% 1st letter GC {1:.2f}% 2nd letter GC {2:.2f}% 3rd letter GC {3:.2f}%\n".format((GC1+GC2+GC3)/30.0, GC1/10.0, GC2/10.0, GC3/10.0)

print(message)

