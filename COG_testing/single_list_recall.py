def fasta_to_dict(file_path):
    fasta_dict = {}
    with open(file_path) as opened_file:
        for line in opened_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                sequence_name = line[1:]
                continue
            sequence = line
            fasta_dict[sequence_name] = sequence
    return fasta_dict

def kmerfile_to_list(file_path):
    kmerlist = []
    with open(file_path) as opened_file:
        for line in opened_file:
            if not line:
                continue
            kmer = line.strip()
            kmerlist.append(kmer)
    return kmerlist

##############################################################################

partition_file = "C:/Users/dylta/Desktop/COG0090.fna"
replist_file = "C:/Users/dylta/Desktop/test_list.txt"

partition_dict = fasta_to_dict(partition_file)
partition_kmerlist = kmerfile_to_list(replist_file)

print("The length of the representative list is: " + str(len(partition_kmerlist)))

i = 0
for sequence_name in partition_dict:
    if all(kmer not in partition_dict[sequence_name] for kmer in partition_kmerlist):
        i+=1

false_positive = len(partition_dict.keys())-i

print("There were " + str(false_positive) + " false positives in " + str(len(partition_dict.keys())) + " sequences.")
            
