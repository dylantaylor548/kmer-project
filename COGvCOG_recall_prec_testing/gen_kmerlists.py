import os, errno
import sys
import argparse

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def main():
    parser = argparse.ArgumentParser(description="Generates kmer lists from a directory of training fasta partitions")
    parser.add_argument("-pd","--partition_directory",help="A directory containing test and trial partitions of a fasta file",required=True)
    parser.add_argument("-o","--out_list_dir",help="Output directory for the generated training lists",required=False,default="training_lists")
    parser.add_argument("-c","--cutoff_value", help="cutoff value - it is the number of times each sequence needs to be covered",required=False,default='1')
    parser.add_argument("-k","--kmer_len", help="kmer length",required=False,default='21')
    args = parser.parse_args()
    if not os.path.exists(args.out_list_dir):
        os.makedirs(args.out_list_dir)

    script_path = get_script_path()
    part_dir = script_path + "/" + str(args.partition_directory)
    
    for file_name in os.listdir(part_dir):
        if file_name.endswith('train.fasta'):
            file_path = part_dir + "/" + file_name
            cmd = "python ../src/kmer_project_master.py -f " + str(file_path) + " -o " + args.out_list_dir + "/" + file_name[:-11] + "replist.csv -c " + args.cutoff_value + " -k " + args.kmer_len
            os.system(cmd)

if __name__ == '__main__':
    main()
