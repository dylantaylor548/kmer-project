import os
import copy
from fractions import Fraction
import sys
# on 500

# For a directory containing files each of one DNA sequence, creates a dictionary
# with keywords: file_name and values: list of all kmers in that file
# (This has been tested and works)

# run this on a file with a specific kmer size for loop 15-31 every 4th one
def kmerize(dna,kmer_size):
    kmers = set()
    kmers_filtered = [] # list
    dna = dna.upper()
    if (kmer_size <= len(dna) and kmer_size >= 1):
        for start in range(0,len(dna)-kmer_size+1,1):
            kmer = dna[start:start+kmer_size]
            kmers.add(kmer)
        kmers_filtered = list(kmers)
        return kmers_filtered

def kmerize_file(file_name,kmer_size):
    kmerdict = {}
    print("Generating " + str(kmer_size)+ "-mers\n...")
    file_path = (file_name)
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
    return kmerdict # puts output list in dictionary w key as name, value as list of kmers for that sequence

# reverses dictionary takes a dictionary that has keys as sequence names and values as list of kmers and
# returns a dictionary that keys are values and values as keys
def reverse_dict(kmerdict):
    rev_dict = {}
    for seq, kmer_list in kmerdict.items():
        for kmer in kmer_list:
            rev_dict.setdefault(kmer,[]).append(seq)
    return rev_dict

# takes a non reversed dictionary and makes values as 0
def generate_seq_counts(kmerdict):
    seq_counts = {}
    for key in kmerdict:
        seq_counts.setdefault(key,0)
    return seq_counts

# if it goes through all of the sequences and are above cutoff
# no sequences at cutoff
def checkcoverage(kmer,seq_kmers_dict,coverage_dict,cutoff):
    for seq in seq_kmers_dict:
        if kmer in seq_kmers_dict[seq]:
            if coverage_dict[seq] - 1 < cutoff:
                return False
    return True

# outputs kmer with highest worth
def findkmax(rev_dict,seq_w_kmer):
    kmax = ''
    for kmer in seq_w_kmer:
        if kmax == '':
            kmax = kmer
        elif seq_w_kmer[kmer] > seq_w_kmer[kmax]:
            kmax = kmer
        # case of a tie, picks one with higher coverage
        elif seq_w_kmer[kmer] == seq_w_kmer[kmax]:
            if len(rev_dict[kmer]) > len(rev_dict[kmax]):
                kmax = kmer
            elif len(rev_dict[kmer]) == len(rev_dict[kmax]):
                kmax = min(kmer,kmax)
    return kmax


def rep_kmers_indict(kmerdict,cutoff):
    rev_dict = reverse_dict(kmerdict)
    seq_counts = generate_seq_counts(kmerdict)
    rep_kmer_dict = {}
    seq_w_kmer = {}
    # makes worth dictionary
    for kmer in rev_dict:
        seq_w_kmer.setdefault(kmer,len(rev_dict[kmer]))
    # keep doing this while the coverage does not meet the cutoff
    while not all(count >= cutoff for count in seq_counts.values()):
        kmer_max = findkmax(rev_dict,seq_w_kmer)
        # growing list of representative kmers (we pick it)
        rep_kmer_dict[kmer_max] = seq_w_kmer[kmer_max]
        for seq in rev_dict[kmer_max]:
            # for all of the sequences in the coverage dict, if that sequence is
            # represented by that chosen kmer
                # for each kmer in the worth dictionary, if that seq is below
                # cutoff coverage, it will drop the worth of that kmer by 1/c
            for kmer in seq_w_kmer:
                if (seq in rev_dict[kmer]) and (seq_counts[seq] < cutoff):
                    seq_w_kmer[kmer] = seq_w_kmer[kmer] - Fraction(1,cutoff)
            seq_counts[seq] += 1
            # save the seq_w_kmer[kmer_max] somehow
        del seq_w_kmer[kmer_max]

    """print("\nChecking output for redundant kmers...\n")"""
    rep_list_copy = copy.deepcopy(rep_kmer_dict)
    for kmer in rep_list_copy.keys():
        # for each seq in the forward dict
        if checkcoverage(kmer,kmerdict,seq_counts,cutoff):
            del rep_kmer_dict[kmer]

    return rep_kmer_dict,rev_dict

# call kmerize file on a file with a for loop for kmer len
#rep_kmers_indict call with cutoff 1
# k len, c (1), worth, coverage
k = 15
c = 1
f = open('kmerdata.csv','w')
f.write('kmerlen'+','+'c'+','+'worth'+','+'coverage'+'\n')

while(k <= 31):
    dict = kmerize_file("C:/Users/Dylan/Desktop/Pop_Lab/kmer-project/TIGR02012_500.fasta", k)
    rep_dict,rev_dict = rep_kmers_indict(dict,c)
    for kmer in rep_dict:
        f.write(str(k)+','+str(c)+','+str(rep_dict[kmer])+','+str(len(rev_dict[kmer]))+'\n')
    k = k + 4

f.close()
