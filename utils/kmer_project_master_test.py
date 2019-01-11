
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
        kmerdict[header] = kmers_in_dna
    print("\nAll " + str(kmer_size)+ "-mers generated!\n")
    return kmerdict


# Removes kmers from a reverse dictionary who's sequences are a subset of another kmers
def rem_redundant(seq_kmers_dict):
    dict_copy = copy.deepcopy(seq_kmers_dict)
    for kmer in dict_copy:
        for kmer2 in dict_copy:
            if kmer != kmer2:
                if set(dict_copy[kmer]) < set(dict_copy[kmer2]):
                    if kmer in seq_kmers_dict:
                        del seq_kmers_dict[kmer]
                elif set(dict_copy[kmer2]) < set(dict_copy[kmer]):
                    if kmer2 in seq_kmers_dict:
                        del seq_kmers_dict[kmer2]
    return seq_kmers_dict


# Reverses a dictionary such that the elements of the value lists become keys
# and the keys become value lists of the elements
# (This has been tested and works)
def reverse_dict(kmerdict):
    rev_dict = {}
    for seq, kmer_list in kmerdict.items():
        for kmer in kmer_list:
            rev_dict.setdefault(kmer,[]).append(seq)
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
    for kmer in seq_w_kmer:
        if kmax == '':
            kmax = kmer
        elif seq_w_kmer[kmer] > seq_w_kmer[kmax]:
            kmax = kmer
        elif seq_w_kmer[kmer] == seq_w_kmer[kmax]:
            if len(rev_dict[kmer]) > len(rev_dict[kmax]):
                kmax = kmer
            elif len(rev_dict[kmer]) == len(rev_dict[kmax]):
                kmax = min(kmer,kmax)
    return kmax


# Given a kmer dictionary and a cutoff, returns a list of kmers such that each
# sequence in the dictionary is represented in the list by at least the cutoff
# (This has been tested and works)
def rep_kmers_indict(kmerdict,cutoff):
    rev_dict = reverse_dict(kmerdict)
    if cutoff == 1:
	print("Removing subsets...")
        rev_dict = rem_redundant(rev_dict)
	print("Done!")
    seq_counts = generate_seq_counts(kmerdict)
    rep_kmer_list = []
    seq_w_kmer = {}
    for kmer in rev_dict:
        seq_w_kmer.setdefault(kmer,len(rev_dict[kmer]))
    while not all(count >= cutoff for count in seq_counts.values()):
        if rev_dict == {}:
            print("Your value of c is too high, not enough unique kmers were found.")
            return
        else:
            kmer_max = findkmax(rev_dict,seq_w_kmer)
            rep_kmer_list.append(kmer_max)
            print("Highest coverage kmer is: " + kmer_max + ", with worth: " + str(seq_w_kmer[kmer_max]) + " and coverage: " + str(len(rev_dict[kmer_max])))
            for seq in seq_counts:
                if seq in rev_dict[kmer_max]:
                    for kmer in seq_w_kmer:
                        if (seq in rev_dict[kmer]) and (seq_counts[seq] < cutoff):
                            seq_w_kmer[kmer] = seq_w_kmer[kmer] - Fraction(1,cutoff)
                    seq_counts[seq] += 1

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

    kmerized_dir = {}
    kmerized_dir = kmerize_directory(args.fasta_file,int(args.kmer_len))

    rep_list = rep_kmers_indict(kmerized_dir, int(args.cutoff_value))
    fw = open(args.out_file, 'w')
    if rep_list:
        for element in rep_list:
            fw.write(element+'\n')
            # print(element)


    print("There are " + str(len(rep_list)) + " " + str(args.kmer_len) + "-mers in the representative list.")

################################################################ ^^^ This is the part that does the code ^^^



if __name__ == '__main__':
    main()

