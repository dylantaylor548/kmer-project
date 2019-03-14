
from itertools import permutations
from copy import deepcopy
import argparse


def gen_kmerlist(file_path):
	oligos = []

	with open(file_path) as opened_file:
		for line in opened_file:
			line = line.strip()
			oligos.append(line)

	return oligos


def gen_overlap_dict(kmerlist):
	perms = permutations(kmerlist,2)

	overlap_dict = {}

	for permutation in perms:
		kmer1 = permutation[0]
		kmer2 = permutation[1]

		overlap_count = 0 
		
		i = 0
		while i <= min([len(kmer1),len(kmer2)]):
			suffix = kmer1[-i:]
			prefix = kmer2[:i]
			if suffix == prefix:
				overlap_count = i
			i += 1

		overlap_dict[kmer1 + ',' + kmer2] = overlap_count

	return overlap_dict


def gen_synth_seq(overlap_dict,synth_size):

	synth_seq_long = ''

	max_perm = ''
	max_overlap = 0

	for oligo_perm in overlap_dict:
		if overlap_dict[oligo_perm] >= max_overlap:
			max_overlap = overlap_dict[oligo_perm]
			max_perm = oligo_perm

	max_perm_list = max_perm.split(',')
	
	overall_left_oligo = max_perm_list[0]
	overall_right_oligo = max_perm_list[1]

	synth_seq_long = overall_left_oligo + overall_right_oligo[max_overlap:]
	print("incorporated 2 kmers")

	for oligo_combo in deepcopy(overlap_dict):
		left_oligo = oligo_combo.split(',')[0]
		right_oligo = oligo_combo.split(',')[1]
		if (left_oligo == overall_left_oligo) or (right_oligo == overall_right_oligo) or (left_oligo == overall_right_oligo and right_oligo == overall_left_oligo):
			del overlap_dict[oligo_combo]


	i = 2
	while overlap_dict != {}:

		left_oligo_combo = ''
		right_oligo_combo = ''
		left_oligo_overlap = 0
		right_oligo_overlap = 0

		for oligo_combo in overlap_dict:
			left_oligo = oligo_combo.split(',')[0]
			right_oligo = oligo_combo.split(',')[1]
			if left_oligo == overall_right_oligo:
				if overlap_dict[oligo_combo] >= right_oligo_overlap:
					right_oligo_combo = oligo_combo
					right_oligo_overlap = overlap_dict[oligo_combo]
			if right_oligo == overall_left_oligo:
				if overlap_dict[oligo_combo] >= left_oligo_overlap:
					left_oligo_combo = oligo_combo
					left_oligo_overlap = overlap_dict[oligo_combo]

		new_left_oligo = left_oligo_combo.split(',')[0]
		new_right_oligo = right_oligo_combo.split(',')[1]

		if left_oligo_overlap >= right_oligo_overlap:
			if left_oligo_overlap > 0:
				synth_seq_long = new_left_oligo[:-left_oligo_overlap] + synth_seq_long
			elif left_oligo_overlap == 0:
				synth_seq_long = new_left_oligo + synth_seq_long

			for oligo_combo in deepcopy(overlap_dict):
				if oligo_combo.startswith(new_left_oligo) or (overall_left_oligo in oligo_combo) or (oligo_combo.startswith(overall_right_oligo) and oligo_combo.endswith(new_left_oligo)):
					del overlap_dict[oligo_combo]
			i += 1
			print("incorporated " + str(i) + " kmers")

			overall_left_oligo = new_left_oligo

		elif right_oligo_overlap > left_oligo_overlap:
			synth_seq_long = synth_seq_long + new_right_oligo[right_oligo_overlap:]

			for oligo_combo in deepcopy(overlap_dict):
				if oligo_combo.endswith(new_right_oligo) or (overall_right_oligo in oligo_combo) or (oligo_combo.startswith(new_right_oligo) and oligo_combo.endswith(overall_left_oligo)):
					del overlap_dict[oligo_combo]
			i += 1
			print("incorporated " + str(i) + " kmers")

			overall_right_oligo = new_right_oligo

		
	return synth_seq_long

################################################################################

def main():
	parser = argparse.ArgumentParser(description='Finds a short set of sequences that collectively comprise all kmers within an input kmer list. These synthetic sequences are used to determine relative abundancies when not all sequences are present in a sample.')
	parser.add_argument("-kl","--kmer_list", help="A .csv file containing the kmer list you are interested in compressing", required=True)
	parser.add_argument("-n","--seq_len", help="Desired maximum length of synthetic sequences", required=False, default=None)
	parser.add_argument("-o","--out_file", help="A .csv file to output the synthetic sequences", required=False)
	args = parser.parse_args()

	kmerlist = gen_kmerlist(args.kmer_list)
	print(str(len(kmerlist)) + " kmers in list")
	overlap_dict = gen_overlap_dict(kmerlist)
	synth_seq = gen_synth_seq(overlap_dict,args.seq_len)
	print("Your synthetic sequence is:\n" + synth_seq)

################################################################################

main()