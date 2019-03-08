from itertools import groupby
import argparse
from math import floor


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
		print(ranges)
		if value in range(ranges[0],ranges[1]+1):
			target = bin_name
			break

	return target

##########################################################################################

def main():
	parser = argparse.ArgumentParser(description="Outputs a graph of oligo locations within input sequences")
	parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required=True)
	parser.add_argument("-oll","--oligo_list", help="A list of oligos (in csv format) from the sequences in the fasta file",required=True)
	parser.add_argument("-b","--bin_size", help="The size of the bins (in bp) for location within the genome",required=False,default='10')
	args = parser.parse_args()

	seqdict = gen_seqdict(args.fasta_file)
	oligo_list = gen_kmerlist(args.oligo_list)

	bins = gen_bins(seqdict,args.bin_size, len(oligo_list[0]))

	for sequence in seqdict:
		for oligo in oligo_list:
			if oligo in sequence:
				location = sequence.index(oligo)
				bin_name = get_bin(location,bindict)
				bins[bin_name] += 1
	
	for bin_name in bins:
		print(bin_name + " : " + str(bins[bin_name]))
				
##########################################################################################

main()