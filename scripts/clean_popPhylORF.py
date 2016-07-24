#! /usr/bin/python

import sys
file = sys.argv[1]

input = open(file, "r")
res = ""
for i in input:
	if i[0] == ">":
		i = i.split(" ")[0] + "\n"
	res = res + i
input.close()

output = open(file + ".new", "w")
output.write(res)
output.close()

