import argparse

def store_kmerlist(csv_filepath):
	oligo_list = []
	with open(csv_filepath) as opened_file:
		for line in opened_file:
			oligo_list.append(line)
	return oligo_list

def main():
	parser = argparse.ArgumentParser(description="Calculates recall and precision statistics for an oligo list against a metagenomic sample")
	parser.add_argument("-f","--sample_fasta",help='A .fasta file containing the sequences in your metagenomic "sample"',required=True)
	parser.add_argument("--id","--identifier",help='The identifier for the correct sequence within the fasta file (i.e. "COG0088")',required=True)
	parser.add_argument("-l","--oligo_list",help='A .csv file containing the representative list from your training data set\n(i.e. only the sequences you are interested in capturing)',required=True)
	parser.add_argument("-o","--output_file",help='An output .csv file to store data about the recall and precision statistics',required=False)
	args = parser.parse_args()

	oligo_list = store_kmerlist(args.oligo_list)



main()