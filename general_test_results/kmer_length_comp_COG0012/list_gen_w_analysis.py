import os
import copy
from fractions import Fraction
import time
import argparse
from itertools import groupby


# Generates the set of all kmers of a given size from an input string (sequence)
def kmerize_seq(seq,k_size):
	kmers = set()

	if (k_size <= len(seq)) and (k_size >= 1):
		for start_loc in range(0, len(seq)-k_size+1):
			kmers.add(seq[start_loc:start_loc+k_size])
		return kmers
	else:
		print('Selected kmer length was either less than 1 or longer than the length of the input sequence')


# Iterates through a fasta file (or file in a similar format) yielding sequence name and actual sequence
# for each item (i.e. gene sequence) in the file
def fasta_iter(file_path):
	fh = open(file_path)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == '>'))
	for header in faiter:
		seq_name = header.next()[1:].strip()
		sequence = ''.join(s.strip() for s in faiter.next())
		yield seq_name, sequence


# Generates dictionaries for (1) the kmers in each sequence and (2) the sequences that contain
# each kmer
def gen_kmer_dicts(file_path, k_size):
	kmers_by_seq = {}
	seqs_by_kmer = {}

	fiter = fasta_iter(file_path)
	for ff in fiter:
		seq_name = ff[0]
		sequence = ff[1]
		kmers = kmerize_seq(sequence, k_size)

		kmers_by_seq[seq_name] = kmers
		
		for kmer in kmers:
			seqs_by_kmer.setdefault(kmer, set())
			seqs_by_kmer[kmer].add(seq_name)

	return kmers_by_seq, seqs_by_kmer


# Generates a dictionary of default coverage (i.e. 0) for each sequence in the kmers_by_seq dictionary
def gen_seq_cov(kmers_by_seq):
	seq_counts = {}

	for seq_name in kmers_by_seq:
		seq_counts[seq_name] = 0
	
	return seq_counts


# Generates a dictionary of default worth (i.e. the number of sequences covered) for 
# each kmers in the seqs_by_kmer dictionary
def gen_worth_dict(seqs_by_kmer):
	worth_dict = {}

	for kmer in seqs_by_kmer:
		worth_dict[kmer] = len(seqs_by_kmer[kmer])

	return worth_dict


# Determines if current coverage of all sequences containing a given kmer are
# above the cutoff coverage (True) or if any are at or below coverage (False)
def check_coverage(kmer, seqs_by_kmer, seq_cov, cutoff):
	for seq_name in seqs_by_kmer[kmer]:
		if seq_cov[seq_name] <= cutoff:
			return False
	return True


# Determines the G/C content of an input DNA sequence (as a fraction)
def getGC(seq):
	GC = 0

	for base in seq:
		if base == 'G' or base == 'C':
			GC += 1

	return Fraction(GC,len(seq))


# Determines the length of the longest run of a single character (base)
# in a string (sequence)
def get_longest_run(seq):
	current_run = 1
	long_run = 0

	for i in range(1, len(seq)):
		if seq[i] == seq[i-1]:
			current_run += 1
		else:
			if current_run > long_run:
				long_run = current_run
			current_run = 1
	if current_run > long_run:
		return current_run

	return long_run


# Determines the kmer with the highest worth, largest number of covered sequences, GC content closest to 50%,
# shortest run of one base in that order (then smallest lexicographically)
def find_kmax(seqs_by_kmer, worth_dict, k_size):
	kmax = ''
	worthcopy = copy.deepcopy(worth_dict)
	
	for kmer in worth_dict:
		if worth_dict[kmer] <= 0:
			del worthcopy[kmer]
		else:
			if kmax == '':
				kmax = kmer
			elif worth_dict[kmer] > worth_dict[kmax]:
				kmax = kmer
			elif worth_dict[kmer] == worth_dict[kmax]:
				kmer_hits = len(seqs_by_kmer[kmer])
				kmax_hits = len(seqs_by_kmer[kmax])
				if kmer_hits > kmax_hits:
					kmax = kmer
				elif kmer_hits == kmax_hits:
					kmer_GC_off50 = abs(getGC(kmer) - Fraction(1,2))
					kmax_GC_off50 = abs(getGC(kmax) - Fraction(1,2))
					if kmer_GC_off50 < kmax_GC_off50:
						kmax = kmer
					elif kmer_GC_off50 == kmax_GC_off50:
						kmer_long_run = get_longest_run(kmer)
						kmax_long_run = get_longest_run(kmax)
						if kmer_long_run < kmax_long_run:
							kmax = kmer
						elif kmer_long_run == kmax_long_run:
							kmax = min(kmer,kmax)

	return kmax, worthcopy


def gen_rep_kmer_list(file_path, k_size, cutoff):
	rep_list = {}

	print("Generating kmers...")
	kmers_by_seq, seqs_by_kmer = gen_kmer_dicts(file_path, k_size)
	print("Done.\n")
	seq_cov = gen_seq_cov(kmers_by_seq)
	worth_dict = gen_worth_dict(seqs_by_kmer)

	cov_counter = len(seq_cov)
	while cov_counter > 0:
		if seqs_by_kmer == {}:
			print("Your value of c is too high, not enough unique kmers were found.")
			return
		else:
			start_time = time.time()
			kmax, worth_dict = find_kmax(seqs_by_kmer, worth_dict, k_size)
			rep_list[kmax] = len(seqs_by_kmer[kmax])
			run_time = time.time() - start_time
			print(str(run_time) + " seconds : Highest coverage kmer is: " + kmax + ", with worth: " + str(worth_dict[kmax]) + ", coverage: " + str(len(seqs_by_kmer[kmax])) + ", and GC content: " + str(getGC(kmax)))
			for seq_name in seqs_by_kmer[kmax]:
				for kmer in kmers_by_seq[seq_name]:
					if seq_cov[seq_name] < cutoff:
						worth_dict[kmer] -= Fraction(1,cutoff)
				seq_cov[seq_name] += 1
				if seq_cov[seq_name] == cutoff:
					cov_counter -= 1
			del worth_dict[kmax]

	print("\nChecking output for redundant kmers...\n")
	redundant_kmers = []
	for kmer in rep_list:
		if check_coverage(kmer, seqs_by_kmer, seq_cov, cutoff):
			redundant_kmers.append(kmer)

	while redundant_kmers != []:
		least_coverage = None
		
		for kmer in redundant_kmers:
			for seq_name in seqs_by_kmer[kmer]:
				if least_coverage == None or seq_cov[seq_name] < least_coverage:
					least_coverage = seq_cov[seq_name]

		most_redundant_kmer = None
		least_cov_counts = None

		for kmer in redundant_kmers:
			count = 0
			for seq_name in seqs_by_kmer[kmer]:
				if seq_cov[seq_name] == least_coverage:
					count += 1
			if least_cov_counts == None or count < least_cov_counts:
				most_redundant_kmer = kmer
				least_cov_counts = count

		del rep_list[most_redundant_kmer]
		redundant_kmers.remove(most_redundant_kmer)
		print(most_redundant_kmer + ' is redundant...')

		for seq_name in seqs_by_kmer[most_redundant_kmer]:
			seq_cov[seq_name] -= 1

		new_redun_kmers = []
		for kmer in redundant_kmers:
			if check_coverage(kmer, seqs_by_kmer, seq_cov, cutoff):
				new_redun_kmers.append(kmer)
		redundant_kmers = new_redun_kmers

	return rep_list


################################################################ vvv This is the part that does the code vvv


def main():
	parser = argparse.ArgumentParser(description="Finds representative kmers for a set of sequences")
	parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required=True)
	parser.add_argument("-o","--out_file", help="Output file with representative kmers",required=True)
	parser.add_argument("-c","--cutoff_value", help="cutoff value - it is the number of times each sequence needs to be covered",required=False,default='1')
	parser.add_argument("-k","--kmer_len", help="The length of the oligos you would like to print on the chip",required=False,default='21')
	args = parser.parse_args()

	start_time = time.time()

	rep_list = gen_rep_kmer_list(args.fasta_file, int(args.kmer_len), int(args.cutoff_value))		

	print("There are " + str(len(rep_list)) + " " + str(args.kmer_len) + "-mers in the representative list.")

	end_time = time.time()
	runtime = end_time - start_time
	seconds = int(runtime % 60)
	minutes = int(runtime/60)

	kmers_by_seq, seqs_by_kmer = gen_kmer_dicts(args.fasta_file, int(args.kmer_len))
	print("It took " + str(minutes) + " minutes and " + str(seconds) + " seconds to generate a list from " + str(len(kmers_by_seq)) + " sequences.")

	fw = open(args.out_file, 'w')
	fw.write('runtime (secs),' + str(runtime) + '\n\n')

	if rep_list:
		for kmer in rep_list:
			fw.write(kmer + ',' + str(rep_list[kmer]) + '\n') 

################################################################ ^^^ This is the part that does the code ^^^



if __name__ == '__main__':
	main()