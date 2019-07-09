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
def fasta2dict_redundant(file_path,genome_cutoff):
	seqdict = {}
	sequences = set()
	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		if sequence in sequences:
			print('The sequence associated with ' + header + " is already represented in this file")
		elif genome_cutoff != None and len(sequence) > int(genome_cutoff):
			print('The sequence associated with ' + header + ' is ' + str(len(sequence)) + ' bases long and has been designated a whole genome')
		else:
			seqdict[header] = sequence
			sequences.add(sequence)
	return seqdict


######################################################################################################################################

def main():
	parser = argparse.ArgumentParser(description='Removes identical sequences from a file in .fasta format or a similar format and optionally removes whole genomes')
	parser.add_argument('-f','--input_file',help='A .fasta file (or similar file) with potentially identical sequences',required=True)
	parser.add_argument('-o','--output_file',help='The file destination of the compressed input file',required=True)
	parser.add_argument('-w','--whole_genome_cutoff',help='Max length a sequence can be before it is designated a whole genome',required=False,default=None)
	args = parser.parse_args()

	reduced_fasta = fasta2dict_redundant(args.input_file,args.whole_genome_cutoff)

	fw = open(args.output_file,'w')
	for seq_id in reduced_fasta:
		fw.write('>' + seq_id + '\n' + reduced_fasta[seq_id] + '\n')
	fw.close()

######################################################################################################################################

main()
