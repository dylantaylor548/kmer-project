import os
from itertools import groupby
from fractions import Fraction


def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


def read_fasta(file_path):

	seqdict = {}

	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seqdict[header] = sequence
	return seqdict


#################################################################################################################

cutoff = 10000

cwd = os.getcwd()

"""whole_genomes = set()

for file_name in os.listdir(cwd):
	if file_name.endswith('sample_test.fasta'):
		file_path = cwd + '/' + file_name
		seqdict = read_fasta(file_path)
		for seq in seqdict:
			if len(seqdict[seq]) > cutoff:
				whole_genomes.add(seq)
		print('done...')

print('\n\nHere are the whole genomes:\n')
whole_genomes = list(whole_genomes)
for seq in whole_genomes:
	print(seq)

print(' ')"""
whole_genomes = ['s28485_COG0124','s29885_COG0052','s27023_COG0495','s9848_COG0081','s7501_COG0081','s4786_COG0552','s33011_COG0085','s30601_COG0552','s40863_COG0201','s28476_COG0080','s37557_COG0049','s34414_COG0552','s32698_COG0049','s32743_COG0087']

for file_name in os.listdir(cwd):
	if file_name.endswith('.report'):
		print('Analyzing ' + file_name)
		COG_count = 0
		file_path = cwd + '/' + file_name
		with open(file_path) as opened_file:
			for line in opened_file:
				for seq_id in whole_genomes:
					if seq_id in line:
						COG_count += 1
		if COG_count != 0:
			outfile_path = cwd + '/' + file_name[:-7] + '.out'
			file_stuff = []
			with open(outfile_path) as opened_file:
				for line in opened_file:
					file_stuff.append(line)
			print('The length of file_stuff is ' + str(len(file_stuff)))
			if len(file_stuff) == 7:
				raw_precision_stats = file_stuff[5].split(' ')
				raw_recall_stats = file_stuff[2].split(' ')
				precision_percent = file_stuff[4].split(' ')

				neg_seqs = int(raw_precision_stats[1]) - COG_count
				true_pos = int(raw_recall_stats[7]) - COG_count
				false_pos = int(raw_precision_stats[7]) - COG_count
				precision = Fraction(true_pos,(true_pos + false_pos))
				precision_str = str(100 * float(precision))[:5]

				raw_precision_stats[1] = str(neg_seqs)
				raw_recall_stats[7] = str(true_pos)
				raw_precision_stats[7] = str(false_pos)
				precision_percent[2] = precision_str + '%\n'

				file_stuff[5] = ' '.join(raw_precision_stats)
				file_stuff[2] = ' '.join(raw_recall_stats)
				file_stuff[4] = ' '.join(precision_percent)

				new_file = outfile_path + '.new'
				fw = open(new_file,'w')
				for line in file_stuff:
					fw.write(line)
				fw.close()