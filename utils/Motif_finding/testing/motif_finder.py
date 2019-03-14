import random

DNA = ['AGCTACGACGAGCTCG','GCTACTTCGGCTAGCT','GTCCATACGATGAAGT','ATTCGACGGCTATAGC','GGTAGGACGAAGCTT']


def randommotifs(DNAstringlist,k):
	motifs = []
	for string in DNAstringlist:
		start = random.choice([x for x in range(0,len(string)-k+1)])
		motif = string[start:start+k]
		motifs.append(motif)
	return motifs


def genprofile(motiflist):
	countsdict = {'A':[1/(len(motiflist)+4)]*len(motiflist[0]),
	'C':[1/(len(motiflist)+4)]*len(motiflist[0]),
	'G':[1/(len(motiflist)+4)]*len(motiflist[0]),
	'T':[1/(len(motiflist)+4)]*len(motiflist[0])}

	for motif in motiflist:
		i = 0
		for base in motif:
			countsdict[base][i] += 1/(len(motiflist)+4)
			i += 1
	return countsdict


def genmotifs(profile,DNAstringlist):
	k = len(profile['A'])
	motifs = []
	for seq in DNAstringlist:
		mostlikelykmer = ''
		maxlikelyhood = 0
		for i in range(0,len(seq)-k+1):
			kmer = seq[i:i+k]
			likelyhood = 1
			i = 0
			for base in kmer:
				likelyhood = likelyhood * profile[base][i]
				i += 1
			if likelyhood > maxlikelyhood:
				mostlikelykmer = kmer
				maxlikelyhood = likelyhood
		motifs.append(mostlikelykmer)
	return motifs


def findconsensus(motifs):
	profile = genprofile(motifs)
	consensus = ''
	for i in range(0,len(profile['A'])):
		maxc = 0
		consensusbase = None
		for base in profile:
			if profile[base][i] >= maxc:
				maxc = profile[base][i]
				consensusbase = base
		consensus += consensusbase
	return consensus


def getmotifscore(motifs,consensus):
	score = 0
	for motif in motifs:
		i = 0
		for base in motif:
			if base != consensus[i]:
				score += 1
	return score

def findbestmotif(DNAstringlist,k):
	motifs = randommotifs(DNA,k)

	i = 0
	while i < 1000:
		prof = genprofile(motifs)
		potentialmotifs = genmotifs(prof,DNAstringlist)
		motifs = potentialmotifs
		consensus = findconsensus(potentialmotifs)
		i += 1

	return findconsensus(motifs)
		

print(findbestmotif(DNA,6))

