a = "ATGGGGCCCAAATTT"

def countCodon(x):
	codons = {}
	for i in range(0, len(x)-3, 3): 
		tmp = a[i:(i+3)]
		if tmp not in codons:
			codons[tmp] = 0
		codons[tmp] += 1
	return(codons)

