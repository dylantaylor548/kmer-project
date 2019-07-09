import argparse
import copy
from itertools import groupby


def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


def fasta2dict(file_path):
	seqdict = {}

	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seqdict[header] = sequence
	return seqdict


def fasta2seq(file_path):
	seq = ''

	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seq = sequence
	return seq


def get_hamming_dist(seq_1,seq_2):
	dist = 0
	compare = ''

	seq1 = copy.deepcopy(seq_1)
	seq2 = copy.deepcopy(seq_2)

	dif = len(seq2) - len(seq1)
	if dif > 0:
		seq1 = seq1 + ('-')*(dif)
	elif dif < 0:
		seq2 = seq2 + ('-')*((-1)*(dif))
	
	for i in range(0, len(seq1)):
		if seq1[i] != seq2[i]:
			dist += 1
			compare += ' '
		else:
			compare += '|'

	return dist, compare


def get_prot_seq(dna_seq):
	prot_seq = ''

	codons = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	} 

	for i in range(0,len(dna_seq),3):
		codon = dna_seq[i:i+3]
		if codon in codons:
			prot_seq += codons[codon]
		else:
			prot_seq += '*'

	return prot_seq


################################################################################################################################################

def main():
	parser = argparse.ArgumentParser(description='Compares a given sequence against an entire FASTA file to determine minimum Hamming Distance')
	parser.add_argument('-s','--sequence',help='A query sequence, either DNA or Protein',required=True)
	parser.add_argument('-t','--query_type',help='Insert "dna" if a DNA sequence or "prot" if a protein sequence',required=True)
	parser.add_argument('-f','--fasta',help='A FASTA file containing the sequences you wish to test against',required=True)
	args = parser.parse_args()

	test_seqs = fasta2dict(args.fasta)

	query = fasta2seq(args.sequence)
	
	min_hamming = None
	min_align = None
	min_seq_id = None
	min_seq = None

	if args.query_type == 'dna':
		for seq_id in test_seqs:
			hamming, align = get_hamming_dist(query,test_seqs[seq_id])
			if min_hamming == None or hamming < min_hamming:
				min_hamming = hamming
				min_align = align
				min_seq_id = seq_id

		min_seq = test_seqs[min_seq_id]
	
	elif args.query_type == 'prot':
		min_hamming = None
		min_align = None
		min_seq_id = None

		for seq_id in test_seqs:
			hamming, align = get_hamming_dist(query,get_prot_seq(test_seqs[seq_id]))
			if min_hamming == None or hamming < min_hamming:
				min_hamming = hamming
				min_align = align
				min_seq_id = seq_id

		min_seq = get_prot_seq(test_seqs[min_seq_id])

	max_len = max(len(query),len(min_seq))

	lines = 0
	if max_len % 100 == 0:
		lines = int(max_len/100)
	else:
		lines = int(max_len/100) + 1

	print(' ')
	for i in range(0,lines-1):
		print('	' + query[(100*i):(100*(i+1))])
		print('	' + min_align[(100*i):(100*(i+1))])
		print('	' + min_seq[(100*i):(100*(i+1))])
		print('')
	print('	' + query[(100*(lines-1)):])
	print('	' + min_align[(100*(lines-1)):])
	print('	' + min_seq[(100*(lines-1)):])

	print('\nThe minimum Hamming distance was: ' + str(min_hamming) + '\n')

################################################################################################################################################

main()