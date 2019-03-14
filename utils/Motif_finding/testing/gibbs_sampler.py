import random
from fractions import Fraction

DNA = ['AGCTACGACGAGCTCG','GCTACTTCGGCTAGCT','GTCCATACGATGAAGT','ATTCGACGGCTATAGC','GGTAGGACGAAGCTT']

def rand_motif(Dna_col,k):
	motifs = []
	for seq in Dna_col:
		start = random.choice([x for x in range(0,len(seq)-k+1)])
		motif = seq[start:start+k]
		motifs.append(motif)
	return motifs


def findcounts(motifs):
	countsdict = {'A':[1]*len(motifs[0]),'C':[1]*len(motifs[0]),'G':[1]*len(motifs[0]),'T':[1]*len(motifs[0])}

	for motif in motifs:
		i = 0
		for base in motif:
			countsdict[base][i] += 1
			i += 1

	return countsdict


def findprofile(motifs):
	counts = findcounts(motifs)
	totals = 0
	for base in counts:
		totals += counts[base][0]
	profile = {'A':[Fraction(x,totals) for x in counts['A']],
	'C':[Fraction(x,totals) for x in counts['C']],
	'G':[Fraction(x,totals) for x in counts['G']],
	'T':[Fraction(x,totals) for x in counts['T']]}

	return profile


def Gibbs(Dna_col, k):
	init_motifs = rand_motif(Dna_col,k)

	while True:
		indices = [x for x in range(0,len(Dna_col))]
		removal = random.choice(indices)
		"""temp_col = Dna_col[0:removal] + Dna_col[removal+1:len(Dna_col)]"""
		temp_motifs = init_motifs[0:removal] + init_motifs[removal+1:len(init_motifs)]
		
	return




motifs = rand_motif(DNA,6)
for motif in motifs:
	print(motif)
profile = findprofile(motifs)
for letter in profile:
	print(letter + ': ' + str(profile[letter]))