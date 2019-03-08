from itertools import groupby
import argparse
from math import floor
import copy
from re import finditer

def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


def gen_seqdict(file_path):
	seqdict = {}

	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seqdict[header] = sequence

	return seqdict


def gen_kmerlist(file_path):
	oligos = []

	with open(file_path) as opened_file:
		for line in opened_file:
			line = line.strip()
			oligos.append(line)

	return oligos


def gen_bins(seqdict,binsize,klen):
	seq_lens = [len(x)-klen for x in seqdict.values()]
	max_length = max(seq_lens)

	bin_dict = {}

	for i in range(0,max_length,binsize):
		bin_name = str(i) + "to" + str(i+binsize-1)
		bin_dict[bin_name] = 0

	return bin_dict


def get_bin(value,bindict):
	target = []
	for bin_name in bindict:
		ranges = [int(x) for x in bin_name.split('to')]
		if value in range(ranges[0],ranges[1]+1):
			target = bin_name
			break

	return target

##########################################################################################

def main():
	parser = argparse.ArgumentParser(description="Outputs a graph of oligo locations within input sequences")
	parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required=True)
	parser.add_argument("-l","--oligo_list", help="A list of oligos (in csv format) from the sequences in the fasta file",required=True)
	parser.add_argument("-o","--output_csv", help="An output containing the frequency of kmers at each binned location",required=True)
	parser.add_argument("-b","--bin_size", help="The size of the bins (in bp) for location within the genome",required=False,default='10')
	args = parser.parse_args()

	seqdict = gen_seqdict(args.fasta_file)
	oligo_list = gen_kmerlist(args.oligo_list)

	klen = len(oligo_list[0])

	bins = gen_bins(seqdict, int(args.bin_size), klen)

	count = 0
	for sequence in seqdict:
		count += 1
		if count % 100 == 0:
			print(str(count) + " sequences checked...")
		for oligo in oligo_list:
			if oligo in seqdict[sequence]:
				locs = [m.start() for m in finditer(oligo,seqdict[sequence])]
				for loc in locs:
					bin_name = get_bin(loc,bins)
					bins[bin_name] += 1
	
	bin_copy = copy.deepcopy(bins)

	f = open(args.output_csv,'w')

	print("compiling data...")

	while bin_copy != {}:
		smallest_bin = 0
		small_bin_name = None
		for bin_name in bin_copy:
			if (int(bin_name.split('to')[0]) <= smallest_bin) or (small_bin_name == None):
				smallest_bin = int(bin_name.split('to')[0])
				small_bin_name = bin_name
		line = small_bin_name + "," + str(smallest_bin) + "," + str(bin_copy[small_bin_name]) + "\n"
		f.write(line)
		del bin_copy[small_bin_name]

	f.close()
				
##########################################################################################

main()