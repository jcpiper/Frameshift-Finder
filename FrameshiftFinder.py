## desktop version of algorithm

## Read in fasta file
dna = open('test.txt', 'r')

### File to be sent to glimmer ###
orfs = open('orfs.csv', 'w+')


## read entire sequence into a string
sequence = dna.read()
dna.close()
sequence.upper()

## find first start codon
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
	
# Loop thru sequence until stop codon encountered
curr = start

codon = sequence[curr:curr+3]
rframe = []
otherFrame = False

while(codon != 'TGA' and codon !='TAA' and codon !='TAG'):
	print(codon)
	rframe.append(codon)
	curr += 3
	codon = sequence[curr:curr+3]
	#check for another start codon
	if (codon == 'ATG' or codon == 'GTG' or codon == 'TTG'):
		otherFrame = True
	

print(rframe)
orfs.close()
