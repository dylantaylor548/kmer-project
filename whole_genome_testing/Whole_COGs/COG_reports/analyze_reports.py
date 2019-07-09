import argparse
import copy
from itertools import groupby, izip
import string
tab = string.maketrans("ACTGN","TGACN")


def avg(listonums):
	total = 0.0 
	for num in listonums:
		total += num
	return round(total/len(listonums),2)


def print_gene_counts(counts):
	counts_copy = copy.deepcopy(counts)
	print('')
	while counts_copy != {}:
		max_appearance = None
		for gene_name in counts_copy:
			if max_appearance == None or counts[gene_name][0] > counts[max_appearance][0]:
				max_appearance = gene_name
		print(max_appearance + ' (' + counts[max_appearance][2] + ') was hit ' + str(counts[max_appearance][0]) + ' times, with ' + str(avg(counts[max_appearance][1])) + ' oligos on average')
		del counts_copy[max_appearance]
	print('')
	return


def print_iter_counts(iter_counts):
	i = 0
	max_gene_name_len = 0
	max_description_len = 0
	max_hit_len = 0
	max_oligos = 0
	for gene_stuff in iter_counts:
		if len(gene_stuff[0][0]) > max_gene_name_len:
			max_gene_name_len = len(gene_stuff[0][0])
		if len(gene_stuff[3]) > max_description_len:
			max_description_len = len(gene_stuff[3])
		if len(str(gene_stuff)[1]) > max_hit_len:
			max_hit_len = len(str(gene_stuff[1]))
		if len("%.2f"%avg(gene_stuff[2])) > max_oligos:
			max_oligos = len("%.2f"%avg(gene_stuff[2]))
	print('')
	for gene_stuff in iter_counts:
		printline = '	[' + str(i) + ']' + (' ')*(len(str(len(iter_counts)-1)) - len(str(i)) + 2) + gene_stuff[0][0] + (' ')*(max_gene_name_len - len(gene_stuff[0][0]) + 2) + '(' + gene_stuff[3] + ')' + (' ')*(max_description_len - len(gene_stuff[3]) + 2) + str(gene_stuff[1]) + ' hit(s)' + (' ')*(max_hit_len - len(str(gene_stuff[1])) + 2) + "%.2f"%avg(gene_stuff[2]) + ' oligo(s) on average' + (' ')*(max_oligos - len("%.2f"%avg(gene_stuff[2])) + 2) + "%.3f"%avg(gene_stuff[4]) + ' average hamming distance'
		print(printline)
		i += 1
	print('')
	return


def gene_counts2list(counts):
	gene_counts_list = []
	counts_copy = copy.deepcopy(counts)
	while counts_copy != {}:
		max_appearance = None
		for gene_name in counts_copy:
			if max_appearance == None or counts[gene_name][0] > counts[max_appearance][0]:
				max_appearance = gene_name
		gene_counts_list.append([[max_appearance], counts[max_appearance][0], counts[max_appearance][1], counts[max_appearance][2], counts[max_appearance][3]])
		del counts_copy[max_appearance]
	return gene_counts_list


def output_counts(countsaslist,outfile):
	fw = open(outfile,'w')
	fw.write('Gene name(s),Description,Hits,Average Oligo presence,Average Hamming Distance\n')
	for gene_stuff in countsaslist:
		fw.write(';'.join(gene_stuff[0]) + ',' + gene_stuff[3] + ',' + str(gene_stuff[1]) + ',' + str(avg(gene_stuff[2])) + ',' + str(avg(gene_stuff[4])) + '\n')
	fw.close()
	return


def hamming_dist(seq1,seq2):
	min_len = min(len(seq1),len(seq2))
	seq1 = seq1 + ' '*(len(seq2) - min_len)
	seq2 = seq2 + ' '*(len(seq1) - min_len)
	return sum(c1 != c2 for c1,c2 in izip(seq1,seq2))



def get_rev_complement(sequence):
	return sequence.translate(tab)[::-1]


def fna_iter(fna_name):
	fh = open(fna_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


# Reads in a fasta file as a dictionary with keys as sequence IDs and values as the sequence associated with that ID
def fna_read(file_path):
	seqdict = {}
	fiter = fna_iter(file_path)
	for ff in fiter:
		seqdict[ff[0].split(' ')[0]] = ff[1]

	return seqdict


def min_hamming_dist(sequence,seqdict):
	min_dist = None
	for seq_id in seqdict:
		forward_dist = hamming_dist(sequence,seqdict[seq_id])
		reverse_dist = hamming_dist(get_rev_complement(sequence),seqdict[seq_id])
		if min_dist == None or min(forward_dist, reverse_dist) < min_dist:
			min_dist = min(forward_dist,reverse_dist)
			if min_dist == 0:
				return min_dist
	return min_dist

########################################################################################################################################################

def main():
	parser = argparse.ArgumentParser(description='From a report of hits from a replist on an annoted whole genome, determines the presence of each hit')
	parser.add_argument('-r','--report',help='A file containing info about where whole genomes where hit by oligos from a replit',required=True)
	parser.add_argument('-o','--output',help='A .csv file containing analysis of the input report',required=True)
	parser.add_argument('-f','--fasta',help='A .fasta (or similar filetype) containing the diversity of the gene of interest the report is based on',required=True)
	args = parser.parse_args()

	query_seqs = fna_read(args.fasta)

	genes_by_id = {}
	oligo_hits = {}

	current_org = ''
	for line in open(args.report):
		if line.startswith('#!genome-build-accession'):
			current_org = line.split(':')[1].strip()
		if line.startswith('	'):
			line_elems = line.strip().split('	')
			if line_elems != ['None']:
				if line_elems[0] == 'gene':
					genes_by_id.setdefault(current_org + line_elems[1],[])
					genes_by_id[current_org + line_elems[1]] = [line_elems[2]] + genes_by_id[current_org + line_elems[1]]
					oligo_hits[current_org + line_elems[1]] = line_elems[3]
				if line_elems[0] == 'CDS':
					genes_by_id.setdefault(current_org + line_elems[1],[])
					genes_by_id[current_org + line_elems[1]].append(line_elems[2])
					genes_by_id[current_org + line_elems[1]].append(line_elems[3])
					genes_by_id[current_org + line_elems[1]].append(line_elems[4])
				if line_elems[0] == 'ncRNA' or line_elems[0] == 'non-coding RNA' or line_elems[0] == 'non-coding':
					genes_by_id.setdefault(current_org + 'gene' + line_elems[1][3:],[])
					genes_by_id[current_org + 'gene' + line_elems[1][3:]].append('non-coding')
					genes_by_id[current_org + 'gene' + line_elems[1][3:]].append(line_elems[2])
					genes_by_id[current_org + 'gene' + line_elems[1][3:]].append(line_elems[3])

	for gene_id in genes_by_id:
		if len(genes_by_id[gene_id]) == 1:
			genes_by_id[gene_id].append('non-coding')
			genes_by_id[gene_id].append(oligo_hits[gene_id])
			genes_by_id[gene_id].append('')


	gene_counts = {}
	possible_genes = {}
	assigned_names = {}



	test_counter = 1
	min_hamming = None
	for gene_id in genes_by_id:
		gene_name = genes_by_id[gene_id][0]
		gene_description = genes_by_id[gene_id][1]
		gene_hits = genes_by_id[gene_id][2]
		if len(genes_by_id[gene_id]) < 4:
			print(gene_id + ' : ' + str(genes_by_id[gene_id]))
		min_hamming = min_hamming_dist(genes_by_id[gene_id][3], query_seqs)
		print("completed gene #" + str(test_counter) + " (" + gene_id + "    " + gene_name + "). Minimum hamming distance was " + str(min_hamming))
		test_counter += 1
		if gene_name in gene_counts:
			gene_counts[gene_name][0] += 1
			gene_counts[gene_name][1].append(int(gene_hits[0]))
			gene_counts[gene_name][3].append(min_hamming)
		elif gene_name not in gene_counts:
			gene_counts[gene_name] = [1,[int(gene_hits[0])],gene_description,[min_hamming]]

	gene_counts_iter = gene_counts2list(gene_counts)

	finished = False
	while not finished:
		print_iter_counts(gene_counts_iter)
		second_check = False
		while not second_check:
			combine = raw_input('If you would like to combine item 3 into item 5 for example, type "3 5" without quotes (using a space).\nIf you are done compressing, type "Done" without quotes (case-sensitive).\n')
			if combine == 'Done':
				finished = True
				break
			elif len(combine.split(' ')) == 2 and all([int(x) in range(0,len(gene_counts_iter)) for x in combine.split(' ')]):
				selections = [int(x) for x in combine.split(' ')]
				second_check = True
				gene_counts_iter[selections[1]][0] += (gene_counts_iter[selections[0]][0])
				gene_counts_iter[selections[1]][1] += gene_counts_iter[selections[0]][1]
				gene_counts_iter[selections[1]][2] += gene_counts_iter[selections[0]][2]
				gene_counts_iter[selections[1]][4] += gene_counts_iter[selections[0]][4]
				del gene_counts_iter[selections[0]]
			else:
				print("Make sure your input is valid.\n")
	output_counts(gene_counts_iter,args.output)

########################################################################################################################################################

main()