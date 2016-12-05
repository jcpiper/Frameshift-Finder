
## desktop version of algorithm

import os


###################################### FUNCTIONS TO FIND SLIPPERY SEQUENCE MATCHES IN ORFS ######################################
# generate all possible slippery sequences of form XXXYYYZ where Z can be anything (including X or Y) and X and Y must be unique
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

## Ignore fasta header
dna.readline()
## read entire sequence into a string
sequence = dna.read()
# strip string of whitespace
sequence.strip()
genomeSize = len(sequence)

readCount = 1 #counter for # of times file is read (should max out at 2)
sequence = sequence.upper()

seqs = genSeqs()

wrappedAround = False

#counter to keep track of bp coordinates
bp = 0

## Loop thru entire genome and find all orfs w/ potential frameshifts
## Keep in mind that bacteriophage genomes are circular, therefore eof != end of genome, must "circle back" to beginning of file and proceed to first in frame stop codon
## ^ how to do this? 

#queue of potential orfs to process
start = []
# list of coordinates of all start codons already investigated
visited = []
curr = 0
while(len(sequence) > 3):
	## debugging output
	# print('starting iteration')
	# print(len(sequence))
	# print('current location: ' + str(genomeSize - len(sequence)))
	## find next start codon
	if not start:
		sequence = sequence[curr:]
		del visited[:]
		## debugging output
		# print('trimmed string!')
		
		# find first ATG start codon
		atg = sequence.find('ATG')
		if (atg >= 0):
			start.append(atg)
		
		# find first GTG start codon
		gtg = sequence.find('GTG')
		if (gtg >= 0):
			start.append(gtg)

		# find first TTG start codon
		ttg = sequence.find('TTG')
		if (ttg >= 0):
			start.append(ttg)
		
		# start processing at first start codon
		start.sort()
		curr = start.pop(0)
		# clear list
		del start[:]
		## debugging output
		# print('Start location: ' + str(curr))
		# print('start contents: ' + str(start))
		# print('Visited? ' + str(visited))
	else:
		curr = start.pop(0)
		## debugging output
		# print('Start Location: ' + str(curr))
		# print('start contents: ' + str(start))
		# print('visited? ' + str(visited))
	
	# Loop thru sequence until stop codon encountered
	codon = sequence[curr:curr+3]
	# bp += 3
	
	# array to hold current reading frame
	rframe = []
	# flag to indicate presence of an out-of-frame start codon within this orf
	otherFrame = False
	
	# flags to indicate whether the interior orfs within this orf have been added to the queue already
	minusFrameVisited = False
	plusFrameVisited = False
	while(codon != 'TGA' and codon !='TAA' and codon !='TAG'):
		## debugging output
		# print(codon + ': ' + str(curr))
		######## HANDLING CIRCULAR NATURE OF GENOME HERE #########
		#circle back thru genome if end of file is not stop codon
		if (len(sequence) - curr < 3 and readCount < 2):
			print('going back to start of genome')
			length = 3 - (len(sequence) - curr)
			print('taking first ' + str(length) + ' characters from start of file')
			print('Remaining sequence: ' + sequence[curr:])
			codon = sequence[curr:]
			# return to beginning of fasta file
			dna.seek(0)
			dna.readline()
			sequence = dna.read()
			readCount += 1
			dna.close()
			codon += sequence[:length]
			print('codon branching end and beginning: ' + codon)
			print('initial sequence = ' + sequence[:length])
			
			#Flag to exit after next stop codon
			wrappedAround = True
			#reset bp counter to current location at beginning of file
			bp = length
			# handle looking for out of frame start codons when circling back
			last = rframe[len(rframe)-1]
			backframe = last[len(last)-1:] + codon[0:2]
			if (backFrame == 'ATG' or backFrame == 'GTG' or backFrame == 'TTG'):
				otherFrame = True
				shift = True
			forwardFrame = sequence[0:3]
			if (forwardFrame == 'ATG' or forwardFrame == 'GTG' or forwardFrame == 'TTG'):
				otherFrame = True
				shift = False
			sequence = sequence[length:]
			curr = length
			rframe.append(codon)
			# print('orf after looping back to beginning of file: ' + str(rframe))
		######### ######### ######### #########
		else:			
			codon = sequence[curr:curr+3]
			rframe.append(codon)
			
			#check for another start codon in +1 and -1 frames
			# shift holds a boolean value --> if True: -1 shift, if False: +1 shift
			# check whether this is the first occurrence of a start in this frame (minusFrameVisited or plusFrameVisited)
			backFrame = sequence[curr-1:curr+2]
			if (backFrame == 'ATG' or backFrame == 'GTG' or backFrame == 'TTG'):
				otherFrame = True
				shift = True
				if not minusFrameVisited:
					minusFrameVisited = True
					if visited.count(curr-1) is 0:
						start.append(curr-1)
						visited.append(curr-1)
				
			forwardFrame = sequence[curr+1:curr+4]
			if (forwardFrame == 'ATG' or forwardFrame == 'GTG' or forwardFrame == 'TTG'):
				otherFrame = True
				shift = False
				if not plusFrameVisited:
					plusFrameVisited = True
					if visited.count(curr+1) is 0:
						start.append(curr+1)
						visited.append(curr+1)

			# move forward in frame
			curr += 3
		
	# rframe.append(codon)
	if (otherFrame):
		orfs.write(','.join(rframe))
		orfs.write(',')
		print(findSlipSeq(''.join(rframe), seqs, shift))
	
	# debugging output:
	# if (len(rframe) > 1):
		# print('ORF: ' + str(rframe))
		# print('remaining sequence: ' + sequence)
	# else:
		# break
	
	#Break out of loop if we already wrapped around to the beginning of the file
	if wrappedAround:
		break
#Remove ',' from end of file
orfs.seek(-1, os.SEEK_END)
orfs.truncate()
orfs.close()
