
import os
import copy
"""from Bio import SeqIO"""
from fractions import Fraction
import time

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


# For a directory containing files each of one DNA sequence, creates a dictionary
# with keywords: file_name and values: list of all kmers in that file
# (This has been tested and works)
"""def kmerize_directory(directory,kmer_size):
    kmerdict = {}
    print("Generating " + str(kmer_size)+ "-mers\n...")
    for file_name in os.listdir(directory):
        if file_name.endswith('.fasta'):
            file_path = (directory + '/' + file_name)
            for seq_record in SeqIO.parse(file_path,'fasta'):
                seq_name = str(seq_record.id).rstrip('\n')
                sequence = str(seq_record.seq).rstrip('\n')
                kmers_in_dna = kmerize(sequence,kmer_size)
                kmerdict[seq_name] = kmers_in_dna
    print("\nAll " + str(kmer_size)+ "-mers generated!\n")
    return kmerdict"""
def kmerize_directory(directory,kmer_size):
    kmerdict = {}
    print("Generating " + str(kmer_size)+ "-mers\n...")
    for file_name in os.listdir(directory):
        if file_name.endswith('.fasta'):
            file_path = (directory + '/' + file_name)
            with open(file_path) as opened_file:
                for line in opened_file:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        sequence_name = line[1:]
                        continue
                    sequence = line
                    kmers_in_dna= kmerize(sequence,kmer_size)
                    kmerdict[sequence_name] = kmers_in_dna
    print("\nAll " + str(kmer_size)+ "-mers generated!\n")
    return kmerdict

# Takes a dictionary whose values are lists of kmers and removes the kmers present
# in all of the lists from each list
# (This has been tested and works, but it may be commented out in the actual code
# if repeating kmers are allowed)
def rem_redun_kmer(kmerdict):
    all_kmers = []
    filtered_dict = copy.deepcopy(kmerdict)
    for kmer_list in filtered_dict.values():
        for kmer in kmer_list:
            all_kmers.append(kmer)
    for kmer in all_kmers:
        if all_kmers.count(kmer) == len(filtered_dict):
            for key in filtered_dict:
                if kmer in filtered_dict[key]:
                    filtered_dict[key].remove(kmer)
    return filtered_dict


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
                            
    """print("\nChecking output for redundant kmers...\n")"""
    rep_list_copy = copy.deepcopy(rep_kmer_list)
    for kmer in rep_list_copy:
        if checkcoverage(kmer,kmerdict,seq_counts,cutoff):
            print(kmer + " is redundant.")
            rep_kmer_list.remove(kmer)
                    
    return rep_kmer_list



############################################################################################################### vvv This is the part that does the code vvv

kmerized_dir = {}

directory = input("Please select the directory containing your .dna files ")

cutoffn = int(input("Please select the minimum kmer coverage for each of your sequences: "))

while True:
    kmer_len = int(input("Please select your desired kmer length: "))
    print("\nLet's see if we can generate " + str(kmer_len) + "-mers for all of your sequences...\n")
    kmerized_dir = kmerize_directory(directory,kmer_len)

    if all(value == None for value in kmerized_dir.values()):
        print("The kmer size you have chosen is either greater than the length of one of your sequences, or less than one.")
        continue
    else:
        break


# This will be commented out if we can choose non-unique kmers
"""kmerized_dir = rem_redun_kmer(kmerized_dir)
print("Done removing non-unique kmers!")"""

rep_list = rep_kmers_indict(kmerized_dir,cutoffn)

if rep_list:
    for element in rep_list:
        print(element)

print("There are " + str(len(rep_list)) + " " + str(kmer_len) + "-mers in the representative list.")

############################################################################################################### ^^^ This is the part that does the code ^^^

