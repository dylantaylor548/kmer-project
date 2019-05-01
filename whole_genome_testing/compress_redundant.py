import argparse
from itertools import groupby


#Reads a fasta file
def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


# Reads in a fasta file as a dictionary with keys as sequence IDs and values as the sequence associated with that ID
def fasta2dict_redundant(file_path):
	seqdict = {}
	sequences = set()
	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		if sequence in sequences:
			print('The sequence associated with ' + header + " is already represented in this file")
		else:
			seqdict[header] = sequence
			sequences.add(sequence)
	return seqdict


######################################################################################################################################

def main():
	parser = argparse.ArgumentParser(description='Removes identical sequences from a file in .fasta format or a similar format')
	parser.add_argument('-f','--input_file',help='A .fasta file (or similar file) with potentially identical sequences',required=True)
	parser.add_argument('-o','--output_file',help='The file destination of the compressed input file',required=True)
	args = parser.parse_args()

	reduced_fasta = fasta2dict_redundant(args.input_file)

	fw = open(args.output_file,'w')
	for seq_id in reduced_fasta:
		fw.write('>' + seq_id + '\n' + reduced_fasta[seq_id] + '\n')
	fw.close()

######################################################################################################################################

main()