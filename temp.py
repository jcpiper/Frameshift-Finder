

	
def genSeqs():
	base = ['A', 'C', 'T', 'G']
	seqs = []
	for x in range(0,4):
		a = base[x] + base[x] + base[x]
		for i in range(0,4):
			b = a + base[i] + base[i] + base[i]
			for n in range(0,4):
				c = b + base[n]
				seqs.append(c)
		
	return seqs
	
	

seqs = genSeqs()
dict = open('sequenceDictionary.txt', 'w')
for i in seqs:
	dict.write(i + '\n')
	
dict.close()