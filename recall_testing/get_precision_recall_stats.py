import argparse
from itertools import groupby
from fractions import Fraction


def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


def store_seqs_w_id(file_path,identifier):
	interest_dict = {}
	noise_dict = {}

	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		if identifier in header:
			interest_dict[header] = sequence
		else:
			noise_dict[header] = sequence
	return interest_dict, noise_dict


def store_kmerlist(oligo_csv_path):
	oligo_list = []
	with open(oligo_csv_path) as opened_file:
		for line in opened_file:
			line = line.strip()
			oligo_list.append(line)
	return oligo_list


def find_rev_complement(sequence):
    complement = ''
    for base in sequence:
        if base == 'A':
            complement += 'T'
        if base == 'C':
            complement += 'G'
        if base == 'G':
            complement += 'C'
        if base == 'T':
            complement += 'A'

    rev_complement = ''
    for base in complement:
        rev_complement = base + rev_complement
    return rev_complement


##################################################################################

def main():
	parser = argparse.ArgumentParser(description="Calculates recall and precision statistics for an oligo list against a metagenomic sample")
	parser.add_argument("-f","--sample_fasta",help='A .fasta file containing the sequences in your metagenomic "sample"',required=True)
	parser.add_argument("-id","--identifier",help='The identifier for the correct sequence within the fasta file (i.e. "COG0088")',required=True)
	parser.add_argument("-l","--oligo_list",help='A .csv file containing the representative list from your training data set\n(i.e. only the sequences you are interested in capturing)',required=True)
	parser.add_argument("-o","--output_file",help='An output .csv file to store data about the recall and precision statistics',required=False,default=None)
	args = parser.parse_args()

	oligo_list = store_kmerlist(args.oligo_list)

	interest_dict, noise_dict = store_seqs_w_id(args.sample_fasta,args.identifier)

	if args.output_file != None:
		f = open(args.output_file,'w')
		line1 = "false_pos_seq_name,COG_name,kmers_hit\n"
		f.write(line1)

	true_pos = 0
	true_neg = 0
	false_pos = 0
	false_neg = 0

	for seq_name in interest_dict:
		covered = False
		for oligo in oligo_list:
			complement = find_rev_complement(oligo)
			if (oligo in interest_dict[seq_name]) or (complement in interest_dict[seq_name]):
				covered = True
				true_pos += 1
				break
		if not covered:
			false_neg += 1

	false_pos_dict = {}

	for seq_name in noise_dict:
		covered = False
		for oligo in oligo_list:
			complement = find_rev_complement(oligo)
			if args.output_file == None:
				if (oligo in noise_dict[seq_name]) or (complement in noise_dict[seq_name]):
					covered = True
					false_pos += 1
					break
			else:
				if (oligo in noise_dict[seq_name]) or (complement in noise_dict[seq_name]):
					covered = True
					if seq_name not in false_pos_dict:
						false_pos += 1
						false_pos_dict[seq_name] = [oligo]
					else:
						false_pos_dict[seq_name].append(oligo)
		if not covered:
			true_neg += 1

	recall = None
	precision = None
	if true_pos != 0:
		recall = float(Fraction(true_pos,(true_pos + false_neg))*100)
		precision = float(Fraction(true_pos,(true_pos + false_pos))*100)
	else:
		recall = 0
		precision = 0

	print("Recall: " + str(recall) + "%")
	print("Precision: " + str(precision) + "%")
	print("Of the " + str(len(noise_dict)) + " sequences that were not " + args.identifier + ", we identified " + str(true_neg) + " as true negatives")

	if args.output_file != None:
		for seq_name in false_pos_dict:
			COG_name = seq_name[-7:]
			kmers = ','.join(false_pos_dict[seq_name])
			line = seq_name + ',' + COG_name + ',' + kmers + '\n'
			f.write(line)
		f.close() 
##################################################################################

main()