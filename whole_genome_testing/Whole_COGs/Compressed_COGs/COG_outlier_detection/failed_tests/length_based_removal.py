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
def fasta2dict(file_path):
	seqdict = {}
	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seqdict[ff[0]] = ff[1]
	return seqdict


def mean(items):
	return float(sum(items))/len(items)


def std_dev(items):
	return sum([(item - mean(items))**2 for item in items])/len(items)


def median(items):
	if len(items)%2 == 0:
		return (float(items[len(items)/2]) + float(items[(len(items)/2)-1]))/2
	else:
		return items[int(len(items)/2)]


def get_range(items,k_mult):
	items.sort()
	middle = median(items)
	bottom_half = []
	top_half = []
	for item in items:
		if item <= middle:
			bottom_half.append(item)
		if item >= middle:
			top_half.append(item)
	q1 = median(bottom_half)
	q3 = median(top_half)
	iqr = q3 - q1

	return q1-(k_mult*iqr) , q3+(k_mult*iqr)

######################################################################################################################################


def main():
	parser = argparse.ArgumentParser(description='Removes sequences from a fasta file whose length is dissimilar to the length of the other sequences in the file')
	parser.add_argument('-f','--input_file',help='A .fasta file (or similar file) with sequences potentially incorrectly identified',required=True)
	parser.add_argument('-o','--output_file',help='The file destination of the compressed input file',required=True)
	parser.add_argument('-k','--outlier_constant',help='The constant used to determine outlier lengths in the data',required=False,default='1.5')
	args = parser.parse_args()

	uncorrected_fasta = fasta2dict(args.input_file)

	seq_lens = []
	for seq_id in uncorrected_fasta:
		seq_lens.append(len(uncorrected_fasta[seq_id]))


	acceptable_min , acceptable_max = get_range(seq_lens,float(args.outlier_constant))
	print('Mean sequence length was ' + str(mean(seq_lens)) + ' bases.')
	print('Finding genes with lengths between ' + str(acceptable_min) + ' bases and ' + str(acceptable_max) + ' bases...')

	fw = open(args.output_file,'w')

	count = 0
	for seq_id in uncorrected_fasta:
		sequence = uncorrected_fasta[seq_id]
		if len(sequence) >= acceptable_min and len(sequence) <= acceptable_max:
			fw.write('>' + seq_id + '\n' + sequence + '\n')
		else:
			count += 1
	print('Done. Removed ' + str(count) + ' sequences of ' + str(len(uncorrected_fasta)) + ' initial sequences.')

	fw.close()		





######################################################################################################################################

main()
