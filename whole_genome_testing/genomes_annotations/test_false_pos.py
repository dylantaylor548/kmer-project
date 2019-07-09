from itertools import groupby, izip
import argparse
import string


def fna_iter(fna_name):
	fh = open(fna_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


# Reads in a fasta file as a dictionary with keys as sequence IDs and values as the sequence associated with that ID
def fna_read(file_path):
	seqdict = {}
	fiter = fna_iter(file_path)
	for ff in fiter:
		seqdict[ff[0].split(' ')[0]] = ff[1]

	return seqdict


tab = string.maketrans("ACTGN","TGACN")

def get_rev_complement(sequence):
	return sequence.translate(tab)[::-1]


def hamming_dist(seq1,seq2):
	min_len = min(len(seq1),len(seq2))
	seq1 = seq1[0:min_len]
	seq2 = seq2[0:min_len]
	return sum(c1 != c2 for c1,c2 in izip(seq1,seq2))


def min_hamming_dist(sequence,seqdict):
	min_dist = None
	min_seq_id = ''
	direction = 'forward'
	for seq_id in seqdict:
		forward_dist = hamming_dist(sequence,seqdict[seq_id])
		reverse_dist = hamming_dist(get_rev_complement(sequence),seqdict[seq_id])
		if min_dist == None or min(forward_dist, reverse_dist) < min_dist:
			min_dist = min(forward_dist,reverse_dist)
			if reverse_dist < forward_dist:
				direction = 'reverse'
			if forward_dist < reverse_dist:
				direction = 'forward'
			min_seq_id = seq_id
			if min_dist == 0:
				return min_dist
	return min_dist, min_seq_id, direction


def print_align(test_seq, query_seq, line_length=100):
	lines = int(max(len(test_seq),len(query_seq))/100)
	for i in range(0,lines+1):
		test_sub = test_seq[i*line_length:(i+1)*line_length]
		query_sub = query_seq[i*line_length:(i+1)*line_length]
		align_sub = ''
		for j in range(0,max(len(test_sub),len(query_sub))):
			if (j >= len(test_sub)) or (j >= len(query_sub)) or (test_sub[j] != query_sub[j]):
				align_sub += ' '
			elif test_sub[j] == query_sub[j]:
				align_sub += '|'
		print('test_seq         ' + test_sub + '\n                 ' + align_sub + '\nquery_max_hit    ' + query_sub + '\n')

#####################################################################################################################


def main():
	parser = argparse.ArgumentParser(description='Determines most alike hit from a fna gene in a set of query sequences')
	parser.add_argument('-w','--whole_genome',help='A file containing the whole genome with the potential false positive',required=True)
	parser.add_argument('-q','--query_seqs',help='A .fasta (or similar) file containing the query sequences to test the potential false positive against',required=True)
	parser.add_argument('-s','--start',help='The start site of the potential hit in the whole genome',required=True)
	parser.add_argument('-e','--end',help='The end site of the potential hit in the whole genome',required=True)
	args = parser.parse_args()

	query_seqs = fna_read(args.query_seqs)
	whole_genome = fna_read(args.whole_genome)

	for genome_id in whole_genome:
		if len(whole_genome[genome_id]) >= int(args.end):
			hit_seq = whole_genome[genome_id][int(args.start)-1:int(args.end)]
			min_hamming, min_seq, direction = min_hamming_dist(hit_seq,query_seqs)
			if direction == 'reverse':
				hit_seq = get_rev_complement(hit_seq)
			print('\n')
			print_align(query_seqs[min_seq],hit_seq)
			print("Minimum Hamming Distance was " + str(min_hamming) + " in sequence " + min_seq)
#####################################################################################################################

main()