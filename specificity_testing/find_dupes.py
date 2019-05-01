import os
from itertools import groupby


def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


def test_fasta(file_path):
	sequences = {}
	count = 0

	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		sequences.setdefault(sequence,0)
		sequences[sequence] += 1
		count += 1

	if len(sequences) < count:
		print(str(len(sequences)) + " sequences are represented " + str(count) + " times in this file...\n")
	else:
		print('clear.\n')
######################################################

cwd = os.getcwd()

print("\n\n")
for file_name in os.listdir(cwd):
	if file_name.endswith("_train_seqs.fasta"):
		print("Analzing " + file_name[:7] + " for duplicate sequences")
		file_path = cwd +"/" + file_name
		test_fasta(file_path)