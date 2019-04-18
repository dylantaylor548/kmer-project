import numpy as np
import matplotlib.pyplot as plt
import os
import copy
from fractions import Fraction
import time
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
        header = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
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
            kmer_max = findkmax(rev_dict,seq_w_kmer)
            rep_kmer_list.append(kmer_max)
            print("Highest coverage kmer is: " + kmer_max + ", with worth: " + str(seq_w_kmer[kmer_max]) + " and coverage: " + str(len(rev_dict[kmer_max])))
            for seq in rev_dict[kmer_max]:
                for kmer in kmerdict[seq]:
                    if seq_counts[seq] < cutoff:
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


############################################################################################################### vvv This is the part that does the code vvv

directory = "C:/Users/Dylan/Desktop/Pop_Lab/kmer-project/data/COG0088/COG0088_10000.fasta"
directory_out = "C:/Users/Dylan/Desktop"
coverage_min = int(input("Input minimum coverage: "))
coverage_max = int(input("Input maximum coverage: "))
interval = int(input("What interval of sequences would you like to calculate runtime for? "))

num_sequences = 8000

rep_list_lens = []
rel_rep_list_lens = []

kstart_time = time.time()
kmerdir = kmerize_directory(directory,21)
kend_time = time.time()
kmerize_time = kend_time - kstart_time

fig1, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1)

for c in range(coverage_min,coverage_max+1):
    run_times = []
    fasta_lens = []
    rel_run_times = []
    for sequences in range(0,num_sequences+1,interval):
        print("Analyzing " + str(sequences) + " sequences at desired coverage c=" + str(c))
        kmerized_dir = {k:v for k,v in list(kmerdir.items())[:sequences]}
        start_time = time.time()
        rep_list = rep_kmers_indict(kmerized_dir,c)
        end_time = time.time()
        run_time = end_time - start_time
        run_time += ((kmerize_time * sequences) / num_sequences)
        fasta_lens.append(sequences)
        run_times.append(run_time)
        rel_run_times.append(run_time/c)
    ax1.plot(fasta_lens,run_times,"o",label='c='+str(c))
    ax2.plot(fasta_lens,rel_run_times,"o",label='c='+str(c))

ax1.set_title("Runtime as a Function of fasta length and c at k = 21")
ax2.set_title("Relative Runtime")
ax1.set_xlabel("# of sequences")
ax1.set_ylabel("Run time (seconds)")
ax1.legend()
ax2.legend()
ax2.set_xlabel("# of sequences")
ax2.set_ylabel("Relative run time (seconds/c)")
plt.tight_layout()
plt.savefig(directory_out+"/Runtime_Analysis"+".png")
plt.show()


############################################################################################################### ^^^ This is the part that does the code ^^^
