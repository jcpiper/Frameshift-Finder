## desktop version of algorithm

import os

# find any slippery sequences of form XXXYYYZ where Z can be anything (including X or Y) and X and Y must be unique
# Ex. sequences: GGGAAAG, GGGAAAA, GGGAAAT	
def genSeqs():
	base = ['A', 'C', 'T', 'G']
	seqs = []
	for x in range(0,4):
		#Generate XXX portion of pattern
		a = base[x] + base[x] + base[x]
		for i in range(0,4):
			#Generate YYY portion of pattern, i.e str is now XXXYYY
			b = a + base[i] + base[i] + base[i]
			for n in range(0,4):
				#Generate Z portion of pattern, end up w/ full XXXYYYZ sequence
				c = b + base[n]
				seqs.append(c)
		
	return seqs
# searches orf for each slippery sequence, returns any matches
def findSlipSeq(orf, seqs):
	for x in seqs:
		result = orf.find(x)
		if result >= 0:
			match = [x, result, result + 7, 'Have to figure out how to determine whether +1 or -1']
			return match
	return 'No match found'
	
	
## Read in fasta file
dna = open('test.txt', 'r')

### File to be sent to glimmer ###
orfs = open('orfs.csv', 'w+')


## read entire sequence into a string
sequence = dna.read()
dna.close()
sequence = sequence.upper()

seqs = genSeqs()

## find first start codon
while(len(sequence) > 3):
	print('starting iteration')
	print(len(sequence))
	start = 0
	atg = sequence.find('ATG')
	if (atg >= 0):
		start = atg
	gtg = sequence.find('GTG')

	if (gtg >= 0 and gtg < start):
		start = gtg

	# make sure start isn't mistakenly left as 0
	if (atg < 0 and gtg >= 0):
		start = gtg

	ttg = sequence.find('TTG')

	if (ttg >= 0 and ttg < start):
		start = ttg

	# make sure start isn't mistakenly left as 0
	if(start == 0 and gtg < 0 and atg < 0 and ttg >= 0):
		start = ttg
	elif (gtg < 0 and atg < 0 and ttg < 0):
		print('breaking out of loop?')
		print('atg: ' + str(atg))
		break #break out of loop if there are no more stop codons
		
	# Loop thru sequence until stop codon encountered
	curr = start

	codon = sequence[curr:curr+3]
	rframe = []
	otherFrame = False

	while(codon != 'TGA' and codon !='TAA' and codon !='TAG'):
		print(codon)
		codon = sequence[curr:curr+3]
		rframe.append(codon)
		sequence = sequence[curr+3:]
		#check for another start codon
		if (codon == 'ATG' or codon == 'GTG' or codon == 'TTG'):
			otherFrame = True
	
	# rframe.append(codon)
	if (otherFrame):
		orfs.write(','.join(rframe))
		orfs.write(',')
		print(findSlipSeq(''.join(rframe), seqs))
		
	print(rframe)
	print('remaining sequence: ' + sequence)
	

#Remove ',' from end of file
orfs.seek(-1, os.SEEK_END)
orfs.truncate()
orfs.close()

