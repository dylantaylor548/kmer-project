import argparse
from itertools import groupby
import string
tab = string.maketrans("ACTGN","TGACN")


def get_align_score(seq1,seq2,match_score,mismatch_score,gap_score):

	score_matrix = [[[0,0,0] for j in range(0,len(seq2)+1)] for i in range(0,len(seq1)+1)]

	counter = 0
	for j in range(0,len(seq2)+1):
		score_matrix[0][j] = [counter*gap_score]*3
		counter += 1

	counter = 0
	for i in range(0,len(seq1)+1):
		score_matrix[i][0] = [counter*gap_score]*3
		counter += 1


	for i in range(1, len(seq1)+1):
		for j in range(1, len(seq2)+1):
			if seq1[i-1] == seq2[j-1]:
				score_matrix[i][j][0] = max(score_matrix[i-1][j-1]) + match_score
			else:
				score_matrix[i][j][0] = max(score_matrix[i-1][j-1]) + mismatch_score
			score_matrix[i][j][1] = max(score_matrix[i][j-1]) + gap_score
			score_matrix[i][j][2] = max(score_matrix[i-1][j]) + gap_score

	return max(score_matrix[len(seq1)][len(seq2)])


def fna_iter(fna_name):
	fh = open(fna_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


def fna_read(file_path):
	seqdict = {}
	fiter = fna_iter(file_path)
	for ff in fiter:
		seqdict[ff[0].split(' ')[0]] = ff[1]

	return seqdict


def get_rev_complement(sequence):
	return sequence.translate(tab)[::-1]


def avg(elements):
	return float(sum(elements))/len(elements)


######################################################################################################################################################################

def main():
	parser = argparse.ArgumentParser(description='Returns a maximum alignment score for two given sequences \n(NOTE: the range of alignment scores will be dependent on the length of the input sequences and assigned alignment values')
	parser.add_argument('-m','--match_score',help='The score assigned when two matching bases align (should be positive)',required=False,default='3')
	parser.add_argument('-n','--mismatch_score',help='The score assigned when two non-matching bases are aligned (should be negative)',required=False,default='-1')
	parser.add_argument('-g','--gap_score',help='The score assigned when there is a gap in one of the two sequences (should be negative)',required=False,default='-5')
	parser.add_argument('-f','--fna_file',help='A file containing the sequences you would like to compare (typically FASTA or FNA)',required=True)
	args = parser.parse_args()

	sequences = fna_read(args.fna_file)
	for seq_id in sequences:
		sequences[seq_id] = [sequences[seq_id],get_rev_complement(sequences[seq_id])]
	alignments = {}
	
	m_score = float(args.match_score)
	n_score = float(args.mismatch_score)
	g_score = float(args.gap_score)
	
	for seq_id1 in sequences:
		alignments[seq_id1] = []
		for seq_id2 in sequences:
			if seq_id1 != seq_id2:
				alignments[seq_id1].append(min(get_align_score(sequences[seq_id1][0],sequences[seq_id2][0],m_score,n_score,g_score),get_align_score(sequences[seq_id1][0],sequences[seq_id2][1],m_score,n_score,g_score)))


	for seq_id in alignments:
		print(seq_id + ' alignment score: ' + str(avg(alignments[seq_id])) + '\n')

######################################################################################################################################################################

main()