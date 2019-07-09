import argparse
from itertools import groupby


#Reads a fasta file
def fna_iter(fna_name):
	fh = open(fna_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


# Reads in a fasta file as a dictionary with keys as sequence IDs and values as the sequence associated with that ID
def fna_read(file_path):
	seq_id = ''
	genome_seq = ''
	fiter = fna_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seq_id = header
		genome_seq = sequence

	return seq_id, genome_seq


##############################################################################################################################################

def main():
	parser = argparse.ArgumentParser(description='From a representative oligo list, determines the location of these kmers in a whole genome')
	parser.add_argument('-f','--whole_fna',help='A .fna or .fasta (or similar file) containing the sequence of a whole genome',required=True)
	args = parser.parse_args()

	seq_id, genome_seq = fna_read(args.whole_fna)

	alphabet = set()
	for base in genome_seq:
		alphabet.add(base)

	print("Bases are: " + str(alphabet))

##############################################################################################################################################


main()