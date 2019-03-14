
import os
import copy
from fractions import Fraction
import time
import argparse
from itertools import groupby

# Takes a single dna sequence and returns all kmers as a list,
# with each kmer appearing only once (This has been tested and works)
def kmerize(dna,kmer_size):
    kmers = set()
    kmers_filtered = []
    dna = dna.upper()
    if (kmer_size <= len(dna) and kmer_size >= 1):
        for start in range(0,len(dna)-kmer_size+1,1):
            kmer = dna[start:start+kmer_size]
            kmers.add(kmer)
        kmers_filtered = list(kmers)
        return kmers_filtered


#Reads a fasta file
def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


# For a directory containing files each of one DNA sequence, creates a dictionary
# with keywords: file_name and values: list of all kmers in that file
# (This has been tested and works)
def kmerize_directory(file_path,kmer_size):
    kmerdict = {}
    print("Generating " + str(kmer_size)+ "-mers\n...")

    fiter = fasta_iter(file_path)
    for ff in fiter:
        header = ff[0]
        sequence = ff[1]
        kmers_in_dna = kmerize(sequence,kmer_size)
        kmerdict[header] = dict((x, 1) for x in kmers_in_dna)
        # kmerdict[header] = kmers_in_dna
    print("\nAll " + str(kmer_size)+ "-mers generated!\n")
    return kmerdict


# Removes kmers from a reverse dictionary who's sequences are a subset of another kmers
# While this is functional, it does not appear to help runtime at all
def rem_redundant(seq_kmers_dict,note=False):
    dict_copy = copy.deepcopy(seq_kmers_dict)
    i = 0
    for kmer in dict_copy:
        if note == True:
            print(str(i) + "/" + str(len(dict_copy)) + "  " + str((100*i)/len(dict_copy)) + "% done")
        for kmer2 in seq_kmers_dict:
            redunlist =[]
            if kmer != kmer2:
                if set(dict_copy[kmer]) < set(dict_copy[kmer2]):
                    redunlist.append(kmer)
                elif set(dict_copy[kmer2]) < set(dict_copy[kmer]):
                    redunlist.append(kmer2)
        for kmers in redunlist:
            if kmers in seq_kmers_dict:
                del seq_kmers_dict[kmers]
        i += 1
    return seq_kmers_dict


# Reverses a dictionary such that the elements of the value lists become keys
# and the keys become value lists of the elements
# (This has been tested and works)
def reverse_dict(kmerdict):
    rev_dict = {}
    for seqs in kmerdict:
        for kmers in kmerdict[seqs]:
            if kmers in rev_dict:
                rev_dict[kmers][seqs] = 1
            else:
                rev_dict[kmers] = {}
                rev_dict[kmers][seqs] = 1
    # for seq, kmer_list in kmerdict.items():
    #     for kmer in kmer_list:
    #         rev_dict.setdefault(kmer,[]).append(seq)
    return rev_dict


# Creates a dictionary with keys as input sequences and values default to 0
# (This has been tested and works)
def generate_seq_counts(kmerdict):
    seq_counts = {}
    for key in kmerdict:
        seq_counts.setdefault(key,0)
    return seq_counts


# This effectively checks whether if the input kmer were removed from the final
# list of kmers, would the coverage of any of the genes drop below the required
# coverage level (This has been tested and works for what it does)
def checkcoverage(kmer,seq_kmers_dict,coverage_dict,cutoff):
    for seq in seq_kmers_dict:
        if kmer in seq_kmers_dict[seq]:
            if coverage_dict[seq] - 1 < cutoff:
                return False
    return True


# In a reference kmer worth dictionary and reverse dictionary, returns the kmer with the highest worth.
# In the case of ties, returns the kmer that covers the most sequences. In the case of further ties,
# returns the kmer with the lowest lexographic value
# (This has been tested and works)
def findkmax(rev_dict,seq_w_kmer):
    kmax = ''
    worthcopy = copy.deepcopy(seq_w_kmer)
    if worthcopy == {}:
    	print("we've run out of kmers")
    for kmer in seq_w_kmer:
        if seq_w_kmer[kmer] == 0:
            del worthcopy[kmer]
        else:
            if kmax == '':
                kmax = kmer
            elif seq_w_kmer[kmer] > seq_w_kmer[kmax]:
                kmax = kmer
            elif seq_w_kmer[kmer] == seq_w_kmer[kmax]:
                if len(rev_dict[kmer]) > len(rev_dict[kmax]):
                    kmax = kmer
                elif len(rev_dict[kmer]) == len(rev_dict[kmax]):
                    kmax = min(kmer,kmax)
    return kmax, worthcopy


# Returns the reverse complement of a DNA sequence
def find_rev_complement(sequence):
    complement = ''
    for base in sequence:
        if base == 'A':
            complement += 'T'
        if base == 'C':
            complement += 'G'
        if base == 'G':
            complement += 'C'
        if base == 'T':
            complement += 'A'

    rev_complement = ''
    for base in complement:
        rev_complement = base + rev_complement
    return rev_complement


# Checks a kmer (as a string) to see if it meets the requirements of oligo design
# Currently, these requirements are limited to: the kmer has at most 3 of same base in a row,
# the kmer has at most 3 of the same dinucleotide in a row,
# the kmer is not self_complementary for four or more bases
def check_kmer_limitations(kmer):
    
    check = True

    for i in range(0,len(kmer)-4):
        substring = kmer[i:i+4]
        if substring == 4 * substring[0]:
            return False

    for i in range(0,len(kmer)-8):
        substring = kmer[i:i+8]
        if substring == 4 * substring[:2]:
            return False

    # Currently checks to see if the kmer is self-complementary over six bases
    for i in range(0,len(kmer)-6):
        substring = kmer[i:i+6]
        rev_complement = find_rev_complement(substring)
        if rev_complement in kmer:
            return False

    return check


# Given a kmer dictionary and a cutoff, returns a list of kmers such that each
# sequence in the dictionary is represented in the list by at least the cutoff
# (This has been tested and works)
def rep_kmers_indict(kmerdict,cutoff):
    rev_dict = reverse_dict(kmerdict)

    print(time.strftime("%c") + ": Removing kmers that do not meet oligo standards...")
    for kmer in rev_dict.keys():
        if not check_kmer_limitations(kmer):
            del rev_dict[kmer]
    print(time.strftime("%c") + ": Done. Kmers cleansed...")

    seq_counts = generate_seq_counts(kmerdict)
    rep_kmer_list = []
    seq_w_kmer = {}
    for kmer in rev_dict:
        seq_w_kmer.setdefault(kmer,len(rev_dict[kmer]))
    cov_counter = len(seq_counts)
    while cov_counter > 0:
        if rev_dict == {}:
            print("Your value of c is too high, not enough unique kmers were found.")
            return
        else:
            kmer_max,seq_w_kmer = findkmax(rev_dict,seq_w_kmer)

            if kmer_max == "":
            	print("\nRan out of kmers that meet the requirements...\nThe sequences that were not covered are:\n")
            	for sequence in seq_counts:
            		if seq_Counts[sequence] == 0:
            			print sequence
            	return rep_kmer_list


            rep_kmer_list.append(kmer_max)
            print(time.strftime("%c") + ": Highest coverage kmer is: " + kmer_max + ", with worth: " + str(seq_w_kmer[kmer_max]) + " and coverage: " + str(len(rev_dict[kmer_max])))
            for seq in rev_dict[kmer_max]:
                for kmer in kmerdict[seq]:
                    if seq_counts[seq] < cutoff:
                    	if kmer in seq_w_kmer:
                        	seq_w_kmer[kmer] -= Fraction(1,cutoff)
                seq_counts[seq] += 1
                if seq_counts[seq] == cutoff:
                    cov_counter -= 1
            del seq_w_kmer[kmer_max]

    print("\nChecking output for redundant kmers...\n")
    rep_list_copy = copy.deepcopy(rep_kmer_list)
    for kmer in rep_list_copy:
        if checkcoverage(kmer,kmerdict,seq_counts,cutoff):
            print(kmer + " is redundant.")
            rep_kmer_list.remove(kmer)
            for seq in kmerdict:
                if kmer in kmerdict[seq]:
                    seq_counts[seq] -= 1
                    
    return rep_kmer_list



################################################################ vvv This is the part that does the code vvv


def main():
    parser = argparse.ArgumentParser(description="Finds representative kmers for a set of sequences")
    parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required=True)
    parser.add_argument("-o","--out_file", help="Output file with representative kmers",required=True)
    parser.add_argument("-c","--cutoff_value", help="cutoff value - it is the number of times each sequence needs to be covered",required=False,default='1')
    parser.add_argument("-k","--kmer_len", help="kmer length",required=False,default='21')
    args = parser.parse_args()

    start_time = time.time()

    kmerized_dir = {}
    kmerized_dir = kmerize_directory(args.fasta_file,int(args.kmer_len))

    rep_list = rep_kmers_indict(kmerized_dir, int(args.cutoff_value))
    fw = open(args.out_file, 'w')
    if rep_list:
        for element in rep_list:
            fw.write(element+'\n')
            # print(element)


    print("There are " + str(len(rep_list)) + " " + str(args.kmer_len) + "-mers in the representative list.")
    
    end_time = time.time()
    runtime = end_time - start_time
    seconds = int(runtime % 60)
    minutes = int(runtime/60)
    print("It took " + str(minutes) + " minutes and " + str(seconds) + " seconds to generate a list from " + str(len(kmerized_dir)) + " sequences.")

################################################################ ^^^ This is the part that does the code ^^^



if __name__ == '__main__':
    main()