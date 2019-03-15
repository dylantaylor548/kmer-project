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
	parser.add_argument("-o","--outfile",help='An output .csv file to store sequences that were retrieved and the kmers that retrieved them',required=False,default=None)
	args = parser.parse_args()

	oligo_list = store_kmerlist(args.oligo_list)

	interest_dict, noise_dict = store_seqs_w_id(args.sample_fasta,args.identifier)

	true_pos_seqs = []
	true_neg_seqs = []
	false_pos_seqs = []
	false_neg_seqs = []

	hit_dict = {}

	for seq_name in interest_dict:
		covered = False
		for oligo in oligo_list:
			if (oligo in interest_dict[seq_name]):
				covered = True
				hit_dict.setdefault(seq_name,[])
				hit_dict[seq_name].append(oligo)
		if covered == True:
			true_pos_seqs.append(seq_name)
		if covered == False:
			false_neg_seqs.append(seq_name)

	for seq_name in noise_dict:
		covered = False
		for oligo in oligo_list:
			complement = find_rev_complement(oligo)
			if oligo in noise_dict[seq_name]:
				covered = True
				hit_dict.set_default(seq_name,[])
				hit_dict[seq_name].append(oligo)
			if complement in noise_dict[seq_name]:
				covered = True
				hit_dict.setdefault(seq_name,[])
				hit_dict[seq_name].append("*" + complement)
		if covered == True:
			false_pos_seqs.append(seq_name)
		if covered == False:
			true_neg_seqs.append(seq_name)

	if args.outfile != None:
		f = open(args.outfile,'w')
		first_line = "COG_id,sequence_id,kmers_in_seq\n"
		f.write(first_line)

		for seq_name in hit_dict:
			COG_id = seq_name[-7:]
			line = COG_id + ',' + seq_name + ',' + (',').join(hit_dict[seq_name]) + '\n'
			f.write(line)

		f.close()

	recall_frac = Fraction(len(true_pos_seqs),(len(true_pos_seqs) + len(false_neg_seqs)))
	precision_frac = Fraction(len(true_pos_seqs),(len(true_pos_seqs) + len(false_pos_seqs)))


	recall_str = str(100 * float(recall_frac))[:5]
	precision_str = str(100 * float(precision_frac))[:5]

	print("\nRecall is " + recall_str + "%")
	print("Of " + str(len(interest_dict)) + " " + args.identifier + " sequences in the sample, " + str(len(true_pos_seqs)) + " were captured.\n")
	print("Precision is " + precision_str + "%")
	print("Of " + str(len(noise_dict)) + " non-" + args.identifier + " sequences in the sample, " + str(len(false_pos_seqs)) + " were captured as false positives")
	if args.outfile != None:
		print("These sequences will be included in the output .csv file")

##################################################################################

main()