import random
import copy
from fractions import Fraction

##
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

##
def reverse_dict(kmerdict):
    rev_dict = {}
    for seq, kmer_list in kmerdict.items():
        for kmer in kmer_list:
            rev_dict.setdefault(kmer,[]).append(seq)
    return rev_dict

##
def generate_seq_counts(kmerdict):
    seq_counts = {}
    for key in kmerdict:
        seq_counts.setdefault(key,0)
    return seq_counts

##
def checkcoverage(kmer,seq_kmers_dict,coverage_dict,cutoff):
    for seq in seq_kmers_dict:
        if kmer in seq_kmers_dict[seq]:
            if coverage_dict[seq] - 1 < cutoff:
                return False
    return True

##
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

##
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
            for seq in seq_counts:
                if seq in rev_dict[kmer_max]:
                    for kmer in seq_w_kmer:
                        if (seq in rev_dict[kmer]) and (seq_counts[seq] < cutoff):
                            seq_w_kmer[kmer] = seq_w_kmer[kmer] - Fraction(1,cutoff)
                    seq_counts[seq] += 1

            del seq_w_kmer[kmer_max]
                            
    rep_list_copy = copy.deepcopy(rep_kmer_list)
    for kmer in rep_list_copy:
        if checkcoverage(kmer,kmerdict,seq_counts,cutoff):
            rep_kmer_list.remove(kmer)
                    
    return rep_kmer_list

##
def partition_fasta(file_path,partition_n,destination_folder,k=21,c=1):
    file = open(file_path)
    file_dict = {}
    
    for line in file:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            sequence_name = line
            continue
        sequence = line
        file_dict[sequence_name] = sequence
        
    sequences = [key for key in file_dict.keys()]
    random.shuffle(sequences)

    partition_length = int(len(sequences)/partition_n)

    for x in range(0,partition_n):
        partition_destination_path = destination_folder + '/partition_' + str(x+1) + '.fasta'
        replist_path = destination_folder + '/partition_' + str(x+1) + '_replist.txt'
        f_partition = open(partition_destination_path,'w')
        f_replist = open(replist_path,'w')
        
        partition_x = sequences[(x*partition_length):((x*partition_length)+partition_length)]
        partition_dict = {}

        print("Generating kmer list for partition " + str(x+1) + "...")
        
        for sequence in partition_x:
            f_partition.write(sequence+'\n'+file_dict[sequence]+'\n')
            partition_dict[sequence] = kmerize(file_dict[sequence],k)
        f_partition.close()

        partition_list = rep_kmers_indict(partition_dict,c)

        for kmer in partition_list:
            f_replist.write(kmer + '\n')
        f_replist.close()

        print("List Generated! There were " + str(len(partition_list)) + " kmers in the representative list.")

        


#######################################################################################3

chosen_file = input('Choose a file to sample: ')
partition_num = int(input('Choose the number of partitions you would like to make: '))
destination = input('Choose a directory to place the partitions in: ')

partition_fasta(chosen_file,partition_num,destination)
