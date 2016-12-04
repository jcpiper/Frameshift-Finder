## desktop version of algorithm

import os


###################################### FUNCTIONS TO FIND SLIPPERY SEQUENCE MATCHES IN ORFS ######################################
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
def findSlipSeq(orf, seqs, shift):
	for x in seqs:
		result = orf.find(x)
		if result >= 0:
			if shift:
				match = [x, result, result + 7, '-1 Frameshift']
			else:
				match = [x, result, result + 7, '+1 Frameshift']
			return match
	return 'No match found'
	
###################################### ###################################### ######################################	

###################################### START OF SCRIPT	######################################
## Read in fasta file
dna = open('etude.fasta', 'r')

### File to be sent to glimmer ###
orfs = open('orfs.csv', 'w+')


## read entire sequence into a string
sequence = dna.read()
dna.close()
sequence = sequence.upper()

seqs = genSeqs()

## Loop thru entire genome and find all orfs w/ potential frameshifts
## Keep in mind that bacteriophage genomes are circular, therefore eof != end of genome, must "circle back" to beginning of file and proceed to first in frame stop codon
## ^ how to do this? 
while(len(sequence) > 3):
	print('starting iteration')
	print(len(sequence))
	## find next start codon
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
		break #break out of loop if there are no more start codons
		
	# Loop thru sequence until stop codon encountered
	curr = start

	codon = sequence[curr:curr+3]
	rframe = []
	otherFrame = False

	while(codon != 'TGA' and codon !='TAA' and codon !='TAG' and len(sequence) > 3):
		print(codon)
		codon = sequence[curr:curr+3]
		rframe.append(codon)
		
		#check for another start codon in frame
		# if (codon == 'ATG' or codon == 'GTG' or codon == 'TTG'):
			# otherFrame = True
			
		#check for another start codon in +1 and -1 frames
		# shift is holds a boolean value --> if True: -1 shift, if False: +1 shift
		backFrame = sequence[curr-1:curr+2]
		if (backFrame == 'ATG' or backFrame == 'GTG' or backFrame == 'TTG'):
			otherFrame = True
			shift = True
		forwardFrame = sequence[curr+1:curr+4]
		if (forwardFrame == 'ATG' or forwardFrame == 'GTG' or forwardFrame == 'TTG'):
			otherFrame = True
			shift = False
		#trim string to remove this codon
		sequence = sequence[curr+3:]
		
	# rframe.append(codon)
	if (otherFrame):
		orfs.write(','.join(rframe))
		orfs.write(',')
		print(findSlipSeq(''.join(rframe), seqs, shift))
	
	# debugging output:
	if (len(rframe) > 1):
		print(rframe)
		print('remaining sequence: ' + sequence)
	else:
		break

#Remove ',' from end of file
orfs.seek(-1, os.SEEK_END)
orfs.truncate()
orfs.close()

