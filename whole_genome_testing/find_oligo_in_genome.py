import argparse
from itertools import groupby


# Returns a list of oligos from a file containing those oligos
def oligos2list(replist_file):
	oligos = []
	fr = open(replist_file)
	for line in fr:
		oligo = line.strip()
		oligos.append(oligo)
	fr.close()
	return oligos


def get_rev_complement(sequence):
	rev_complement = ''
	alphabet = ['A','C','T','G']
	for base in sequence:
		rev_complement = alphabet[(alphabet.index(base)-2)%4] + rev_complement
	return rev_complement


#Reads a fasta file
def fna_iter(fna_name):
	fh = open(fna_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


# Reads in a fasta file as a dictionary with keys as sequence IDs and values as the sequence associated with that ID
def fna2dict(file_path):
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
	parser.add_argument('-l','--oligo_list',help='A file containing representative oligos',required=True)
	parser.add_argument('-f','--whole_fna',help='A .fna or .fasta (or similar file) containing the sequence of a whole genome',required=True)
	parser.add_argument('-g','--whole_gff',help='A .gff file containing annotations of the input whole genome file',required=False)
	args = parser.parse_args()

	seq_id, genome_seq = fna2dict(args.whole_fna)
	oligos = oligos2list(args.oligo_list)
	oligo_len = len(oligos[0])


	locations = {}
	i = 0
	while i <= len(genome_seq)-oligo_len:
		frame = genome_seq[i:i+oligo_len]
		for oligo in oligos:
			if oligo == frame:
				locations.setdefault(oligo,[])
				locations[oligo].append(i)
			rev_comp = get_rev_complement(oligo)
			if rev_comp == frame:
				locations.setdefault('*' + rev_comp,[])
				locations['*' + rev_comp].append(i)
		i += 1

	print(locations)

##############################################################################################################################################

main()