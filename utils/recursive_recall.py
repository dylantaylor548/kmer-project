import os
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Generates kmer lists from a directory of training fasta partitions")
    parser.add_argument("-pd","--partition_directory",help="A directory containing test and trial partitions of a fasta file",required=True)
    parser.add_argument("-o","--out_list_dir",help="Output directory name for the generated training lists",required=False,default="recall_stats")
    parser.add_argument("-stp","--stats_pre",help="Output statistics file suffix",required=False,default='stats')
    parser.add_argument("-sump","--sum_pre",help="Output summary file suffix",required=False,default='summary')
    args = parser.parse_args()

    output_path = args.partition_directory + "/" + args.out_list_dir

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    for file_name in os.listdir(args.partition_directory):
        if file_name.endswith('replist.csv'):
			file_path = args.partition_directory + "/" + file_name
			prefix = file_path[:-11]
			cmd = 'python kmerlist_recall.py -f ' + prefix + 'test.fasta -r ' + prefix + 'replist.csv -s ' + output_path + '/' + file_name[:-11] + args.stats_pre + ' -m ' + output_path + '/' + file_name[:-11] + args.sum_pre
			os.system(cmd)

if __name__ == '__main__':
    main()
