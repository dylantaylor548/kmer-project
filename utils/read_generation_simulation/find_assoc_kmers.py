import argparse
from itertools import  groupby


def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


def fasta2dict(file_path):
	seqdict = {}
	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seqdict[header] = sequence
	return seqdict


def oligos2list(file_path):
	oligos = []
	fo = open(file_path)
	for line in fo:
		oligos.append(line.strip())
	return oligos


def reverse_dict(dictionary):
	rev_dict = {}
	for key in dictionary:
		for value in dictionary[key]:
			rev_dict.setdefault(value,[])
			rev_dict[value].append(key)
	return rev_dict

def main():
	parser = argparse.ArgumentParser(description='Determines the presence of a select list of kmers within a select set of sequences')
	parser.add_argument('-f','--fasta',help='A file containing the sequences of interest',required=True)
	parser.add_argument('-l','--kmer_list',help='A list of kmers (in .csv format) to be found within your selected sequences',required=True)
	parser.add_argument('-o','--output',help='Output containing report of kmers associated with sequences of interest',required=True)
	args=parser.parse_args()


	sequences = fasta2dict(args.fasta)
	oligos = oligos2list(args.kmer_list)

	fw = open(args.output,'w')
	fw.write('seq_id,kmers,shared_seqs')

	kmers_in_seqs = {}
	for seq in sequences:
		kmers_in_seqs[seq] = []
		for kmer in oligos:
			if kmer in sequences[seq]:
				kmers_in_seqs[seq].append(kmer)

	seqs_w_kmer = reverse_dict(kmers_in_seqs)

	for seq_id in kmers_in_seqs:
		kmer_count = str(len(kmers_in_seqs[seq_id]))
		shared_seqs = set()
		shared_seqs.add(seq_id)
		for kmer in kmers_in_seqs[seq_id]:
			for seq in seqs_w_kmer[kmer]:
				shared_seqs.add(seq)
		shared_seqs_count = str(len(shared_seqs) - 1)
		print(seq_id + ',' + kmer_count + ',' + str(shared_seqs))
		fw.write('\n' + seq_id + ',' + kmer_count + ',' + shared_seqs_count)
	fw.close()


main()