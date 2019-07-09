import os
import sys
import argparse
import copy


def main():
    parser = argparse.ArgumentParser(description="Generates kmer lists from a directory of training fasta partitions")
    parser.add_argument("-rd","--replist_directory",help="A directory containing representative kmer lists",required=True)
    parser.add_argument("-pd","--partition_directory",help="A directory containing partitioned fasta files",required=True)
    parser.add_argument("-s","--stats_suf",help="Output statistics file suffix",required=False,default='stats')
    parser.add_argument("-m","--sum_suf",help="Output summary file suffix",required=False,default='summary')
    args = parser.parse_args()
    
    dircopy = copy.deepcopy(os.listdir(args.replist_directory))
    for file_name in dircopy:
        if file_name.endswith('replist.csv'):
			file_path = args.replist_directory + "/" + file_name
			prefix = file_name[:-11]
			cmd = 'python kmerlist_recall.py -f ' + args.partition_directory + '/' + prefix + 'test.fasta -r ' + file_path + ' -s ' + args.replist_directory + '/' + prefix + args.stats_suf + ' -m ' + args.replist_directory + '/' + prefix + args.sum_suf
			os.system(cmd)

if __name__ == '__main__':
    main()
