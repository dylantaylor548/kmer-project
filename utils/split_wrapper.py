import os, errno
import sys
import argparse

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def main():
    parser = argparse.ArgumentParser(description="Cross validation setup")
    parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required = True)
    parser.add_argument("-r","--train_test_ratio", help="Training to test ratio (between 0 to 1, default=0.8)",required = False, default = 0.8)
    parser.add_argument("-n ","--num_trials", help="Number of random trials (default = 10)",required = False, default = 10)
    parser.add_argument("-o","--out_dir", help="Output directory", required = False, default = 'temp')
    parser.add_argument("-pre","--prefix", help="Output file prefix", required = False, default ='trial')
    args = parser.parse_args()
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    script_path = get_script_path()
    for t in range(0, int(args.num_trials)):
        op_file_path = str(args.out_dir) + '/' + str(args.prefix) + '_' + str(t)
        cmd = str(script_path) +'/split.sh ' + args.fasta_file + ' ' + str(args.train_test_ratio) + ' ' + op_file_path
        os.system(cmd)

if __name__ == '__main__':
    main()