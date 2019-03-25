
import os
import copy
from fractions import Fraction
import time
import argparse
from itertools import groupby

# Takes a single dna sequence and returns all kmers as a dictionary,
# with each kmer's value being a list of it's locations within the sequence
def kmerize(dna,kmer_size):
	kmers = {}
	dna = dna.upper()
	if (kmer_size <= len(dna) and kmer_size >= 1):
		for start in range(0,len(dna)-kmer_size+1,1):
			kmer = dna[start:start+kmer_size]
			if kmer in kmers:
				kmers[kmer].append(start)
			elif kmer not in kmers:
				kmers[kmer] = [start]
		return kmers


#Reads a fasta file
def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


# For a fasta file containing all known diversity of a gene of interest, creates a dictionary
# where the keys are the sequence id's of each sequence in the file and the values are themselves dictionaries
# with keys as the kmers that appear in that sequence and values as a list of locations within that sequence the
# kmer appears
def kmerize_fasta(file_path,kmer_size):
	kmerdict = {}
	print("Generating " + str(kmer_size)+ "-mers\n...")

	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		kmers_in_dna = kmerize(sequence,kmer_size)
		kmerdict[header] = kmers_in_dna
		# kmerdict[header] = kmers_in_dna
	print("\nAll " + str(kmer_size)+ "-mers generated!\n")
	return kmerdict


# Reverses the kmerdict such that the keys are all the kmers in the fasta file
# and the values are themselves dictionaries with keys as the sequences that contain
# that kmer and the values being a list of the (starting) locations of that kmer 
# within that sequence
def reverse_dict(kmerdict):
	rev_dict = {}
	for seq in kmerdict:
		for kmer in kmerdict[seq]:
			if kmer in rev_dict:
				if seq in rev_dict[kmer]:
					rev_dict[kmer][seq].append(kmerdict[seq][kmer])
				else:
					rev_dict[kmer][seq] = kmerdict[seq][kmer]
			else:
				rev_dict[kmer] = {}
				rev_dict[kmer][seq] = kmerdict[seq][kmer]
	return rev_dict


# Creates a dictionary with keys as input sequences (by id's) and values default to 0
def generate_seq_counts(kmerdict):
	seq_counts = {}
	for key in kmerdict:
		seq_counts.setdefault(key,0)
	return seq_counts


# This effectively checks whether if the input kmer were removed from the final
# list of kmers, would the coverage of any of the genes drop below the required
# coverage level
def checkcoverage(kmer,seq_kmers_dict,coverage_dict,cutoff):
	for seq in seq_kmers_dict:
		if kmer in seq_kmers_dict[seq]:
			if coverage_dict[seq] - 1 < cutoff:
				return False
	return True


# For a given kmer, checks whether that kmer could possibly end up on the same fragment
# as another kmer that has already been chosen
def test_kmer_on_frag(kmer,rev_dict,chosen_locs,frag_size,k_len):

	for seq in rev_dict[kmer]:
		for kmer_loc in rev_dict[kmer][seq]:
			for loc in chosen_locs[seq]:
				if abs(kmer_loc - loc) <= (frag_size - k_len):
					return False

	return True


# Determines the GC content (as a decimal) of a string
def getGC(dna):
	GC = 0
	dna = dna.upper()
	i = 0
	while i < len(dna):
		base = dna[i]
		if base == 'G' or base == 'C':
			GC = GC + 1
		i += 1
	content = Fraction(GC, len(dna))
	return content


# Finds the longest run of a single letter (base) within a string (DNA sequence)
def getLongRun(dna):
	current_run = 1
	long_run = 0
	i = 1
	while i < len(dna):
		if dna[i] == dna[i-1]:
			current_run += 1
		else:
			if current_run > long_run:
				long_run = current_run
			current_run = 1
		i += 1
	if current_run > long_run:
		long_run = current_run
	return long_run


# In a reference kmer worth dictionary and reverse dictionary, returns the kmer with the highest worth.
# In the case of ties, returns the kmer that covers the most sequences. Then the kmer with GC content
# closest to 50%, then the kmer with the shortest run of one base. In the case of further ties,
# returns the kmer with the lowest lexographic value
def findkmax(rev_dict,worth_dict,chosen_kmer_locs,frag_size,k_len):
	kmax = ''
	worthcopy = copy.deepcopy(worth_dict)
	for kmer in worth_dict:
		if worth_dict[kmer] <= 0:
			del worthcopy[kmer]
		else:
			test_kmer = test_kmer_on_frag(kmer,rev_dict,chosen_kmer_locs,frag_size,k_len)
			if test_kmer == True:
				if kmax == '':
					kmax = kmer
				elif worth_dict[kmer] > worth_dict[kmax]:
					kmax = kmer
				elif worth_dict[kmer] == worth_dict[kmax]:
					if len(rev_dict[kmer]) > len(rev_dict[kmax]):
						kmax = kmer
					elif len(rev_dict[kmer]) == len(rev_dict[kmax]):
						if abs(getGC(kmer) - Fraction(1,2)) < abs(getGC(kmax) - Fraction(1,2)):
							kmax = kmer
						elif abs(getGC(kmer) - Fraction(1,2)) == abs(getGC(kmax) - Fraction(1,2)):
							if getLongRun(kmer) < getLongRun(kmax):
								kmax = kmer
							elif getLongRun(kmer) == getLongRun(kmax):
								kmax = min(kmer,kmax)

	for seq in rev_dict[kmax]:
		for loc in rev_dict[kmax][seq]:
			chosen_kmer_locs[seq].append(loc)

	return kmax, worthcopy, chosen_kmer_locs


# Given a kmer dictionary and a cutoff, returns a list of kmers such that each
# sequence in the dictionary is represented in the list by at least the cutoff
# (This has been tested and works)
def rep_kmers_indict(kmerdict,cutoff,frag_size,k_len):
	rev_dict = reverse_dict(kmerdict)
	seq_counts = generate_seq_counts(kmerdict)
	rep_kmer_list = []
	
	worth_dict = {}
	for kmer in rev_dict:
		worth_dict.setdefault(kmer,len(rev_dict[kmer]))

	chosen_kmer_locs = {}
	for seq in kmerdict:
		chosen_kmer_locs[seq] = []

	cov_counter = len(seq_counts)
	while cov_counter > 0:
		if rev_dict == {}:
			print("Your value of c is too high, not enough unique kmers were found.")
			return
		else:
			kmer_max,worth_dict,chosen_kmer_locs = findkmax(rev_dict,worth_dict,chosen_kmer_locs,frag_size,k_len)
			rep_kmer_list.append(kmer_max)
			print(time.strftime("%c") + ": Highest coverage kmer is: " + kmer_max + ", with worth: " + str(worth_dict[kmer_max]) + ", coverage: " + str(len(rev_dict[kmer_max])) + ", and GC content: " + str(getGC(kmer_max)))
			for seq in rev_dict[kmer_max]:
				for kmer in kmerdict[seq]:
					if seq_counts[seq] < cutoff:
						worth_dict[kmer] -= Fraction(1,cutoff)
				seq_counts[seq] += 1
				if seq_counts[seq] == cutoff:
					cov_counter -= 1
			del worth_dict[kmer_max]

	print("\nChecking output for redundant kmers...\n")
	redundant_kmers = []
	for kmer in rep_kmer_list:
		if checkcoverage(kmer,kmerdict,seq_counts,cutoff):
			redundant_kmers.append(kmer)

	while redundant_kmers != []:
		
		least_coverage = None
		for kmer in redundant_kmers:
			for seq in rev_dict[kmer]:
				if least_coverage == None:
					least_coverage = seq_counts[seq]
				elif seq_counts[seq] < least_coverage:
					least_coverage = seq_counts[seq]
		
		most_redundant_kmer = ''
		least_cov_counts = None

		for kmer in redundant_kmers:
			count = 0
			for seq in rev_dict[kmer]:
				if seq_counts[seq] == least_coverage:
					count += 1
			if least_cov_counts == None:
				most_redundant_kmer = kmer
				least_cov_counts = count
			elif count < least_cov_counts:
				most_redundant_kmer = kmer
				least_cov_counts = count

		rep_kmer_list.remove(most_redundant_kmer)
		redundant_kmers.remove(most_redundant_kmer)
		print(most_redundant_kmer + " is redundant...")

		for seq in rev_dict[most_redundant_kmer]:
			seq_counts[seq] -= 1
		
		new_redun_kmers = []
		for kmer in redundant_kmers:
			if checkcoverage(kmer,kmerdict,seq_counts,cutoff):
				new_redun_kmers.append(kmer)
		redundant_kmers = new_redun_kmers

	return rep_kmer_list



################################################################ vvv This is the part that does the code vvv


def main():
	parser = argparse.ArgumentParser(description="Finds representative kmers for a set of sequences")
	parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required=True)
	parser.add_argument("-o","--out_file", help="Output file with representative kmers",required=True)
	parser.add_argument("-c","--cutoff_value", help="cutoff value - it is the number of times each sequence needs to be covered",required=False,default='1')
	parser.add_argument("-k","--kmer_len", help="The length of the oligos you would like to print on the chip",required=False,default='21')
	parser.add_argument("-s","--frag_size", help="The expected size of the fragments you will be washing your chip with",required=False,default='100')
	args = parser.parse_args()

	start_time = time.time()

	kmerized_fasta = {}
	kmerized_fasta = kmerize_fasta(args.fasta_file,int(args.kmer_len))

	rep_list = rep_kmers_indict(kmerized_fasta, int(args.cutoff_value),int(args.frag_size), int(args.kmer_len))
	fw = open(args.out_file, 'w')
	if rep_list:
		for element in rep_list:
			fw.write(element+'\n')
			"""line = element
			for seq in reverse_dict(kmerized_fasta)[element]:
				line += ',' + seq
			line += '\n'
			fw.write(line)"""


	print("There are " + str(len(rep_list)) + " " + str(args.kmer_len) + "-mers in the representative list.")
	
	end_time = time.time()
	runtime = end_time - start_time
	seconds = int(runtime % 60)
	minutes = int(runtime/60)
	print("It took " + str(minutes) + " minutes and " + str(seconds) + " seconds to generate a list from " + str(len(kmerized_fasta)) + " sequences.")

################################################################ ^^^ This is the part that does the code ^^^



if __name__ == '__main__':
	main()
