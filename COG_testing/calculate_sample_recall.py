import os

def directory_to_dictlist(directory):
    dict_list = []
    for file_name in os.listdir(directory):
        partition_dict = {}
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
                    partition_dict[sequence_name] = sequence
                dict_list.append(partition_dict)
    return dict_list

def directory_to_replistdict(directory):
    replistdict = {}
    for file_name in os.listdir(directory):
        kmer_list = []
        if file_name.endswith('.txt'):
            file_path = file_path = (directory + '/' + file_name)
            with open(file_path) as opened_file:
                for line in opened_file:
                    if not line:
                        continue
                    kmer = line.strip()
                    kmer_list.append(kmer)
                replistdict[file_name.rstrip('.txt')] = kmer_list
    return replistdict

#############################################

partition_directory = "C:/Users/dylta/Desktop/Pop lab/kmer-project/COG_Recall_Testing/COG0088/COG0088-partitions_c2"
replist_directory = "C:/Users/dylta/Desktop/Pop lab/kmer-project/COG_Recall_Testing/COG0088/COG0088-replists_c2"
output_location = "C:/Users/dylta/Desktop/COG0088_Recall.csv"

partition_dictlist = directory_to_dictlist(partition_directory)
replistdict = directory_to_replistdict(replist_directory)

num_partitions = len(partition_dictlist)

f = open(output_location,'w')

first_line = ''
for i in range(1,num_partitions+1):
    first_line += ',Partition_' + str(i)
first_line += ',Recall\n'

f.write(first_line)


for replist in replistdict.keys():
    file_line = replist
    total_coverage = 0
    max_coverage = 0
    seq_counted = 0
    for partition_dict in partition_dictlist:
        coverage_count = 0
        for sequence in partition_dict:
            seq_counted += 1
            for kmer in replistdict[replist]:
                if kmer in partition_dict[sequence]:
                    coverage_count += 1
                    break
        file_line += ',' + str(coverage_count)
        total_coverage += coverage_count
        if coverage_count > max_coverage:
            max_coverage = coverage_count
    recall = ((total_coverage - max_coverage)*100)/(seq_counted - max_coverage)
    file_line += ',' + str(recall) + '\n'
    f.write(file_line)

f.close()
    
