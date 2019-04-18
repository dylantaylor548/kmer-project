from itertools import groupby
import random
from numpy.random import choice
from fractions import Fraction
from math import floor
import argparse


def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq


def fasta2dict(file_path):
	seqdict = {}
	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seqdict[header] = sequence
	return seqdict

def chip2list(file_path):
	chip_list = []
	ch = open(file_path)
	for line in ch:
		line = line.split(',')
		copies = [line[0]] * int(line[1])
		chip_list += copies
	return chip_list

########################################################################################

def main():
	parser = argparse.ArgumentParser(description='Simulates washing a metagenomic sample over our oligo chip, and outputs a fastq file of reads')
	parser.add_argument('-f','--fasta_file',help='The .fasta file of sequence diversity used to generate your oligo list',required=True)
	parser.add_argument('-c','--chip_file',help='The .csv file representing the abundance of oligos on your chip',required=True)
	parser.add_argument('-s','--spikes_file',help='The .fasta file containing the spike sequence(s) used to validate abundances',required=True)
	parser.add_argument('-o','--out_file',help='An output file containing your reads in .fastq format',required=True)
	parser.add_argument('-l','--frag_size',help='The estimated size of fragments generated for sequencing (i.e. read length)',required=False,default='200')
	args = parser.parse_args()

	fasta = args.fasta_file
	chip_file = args.chip_file
	spikes = args.spikes_file
	outfile = args.out_file
	frag_size = int(args.frag_size)
	
	read_choices = 50

	seqs = fasta2dict(fasta)
	spike_seqs = fasta2dict(spikes)
	for spike in spike_seqs:
		seqs[spike] = spike_seqs[spike]

	chip = chip2list(chip_file)

	simulated_presence = {}
	for seq in seqs:
		presence = 10000 * random.choice([x for x in range(50,read_choices+1)])  # This has been edited such that all sequences will be present 500000 times
		if seq in spike_seqs:
			presence = 10000 * random.choice([x for x in range(50,read_choices+1)])  # This has been edited such that the spike sequnece will be present 500000 times
		simulated_presence[seq] = presence

	output = open('C:/Users/Dylan/Desktop/tmp_sim_pres.csv','w')
	for seq in simulated_presence:
		output.write(seq + ',' + str(simulated_presence[seq]) + '\n')
	output.close()

	print('Generating fragments from .fasta and sequence presence data...')
	possible_frags = {}
	for seq in simulated_presence:
		for i in range(0,len(seqs[seq])-frag_size+1):
			frag = seqs[seq][i:i+frag_size]
			for kmer in list(set(chip)):
				if kmer in frag:
					possible_frags.setdefault(frag + "." + seq,0)
					possible_frags[frag + "." + seq] += simulated_presence[seq]
					break

	print('Useful fragments generated...')

	random.shuffle(chip)

	fw = open(outfile,'w')
	i = 1

	frag_list = list(possible_frags.keys())
	frag_presence = list(possible_frags.values())
	total_presence = sum(frag_presence)
	frag_presence = [Fraction(x,total_presence) for x in frag_presence]
	while chip != []:
		chosen_oligo = random.choice(chip)
		chosen_frag = choice(frag_list,p=frag_presence)
		read = chosen_frag.split(".")[0]
		if chosen_oligo in read:
			seq = chosen_frag.split(".")[1]
			chip.remove(chosen_oligo)
			print(str(len(chip)) + " oligos left on chip")
			possible_frags[chosen_frag] -= 1	
			fw.write('@NR_00000000'[:12-len(str(i))] + str(i) + '.' + seq + '\n')
			print('@NR_00000000'[:12-len(str(i))] + str(i) + '.' + seq + '\n')
			fw.write(read + '\n')
			fw.write('+\n')
			fw.write(('!' * len(read) + '\n'))
			i += 1
			frag_presence = list(possible_frags.values())
			total_presence = sum(frag_presence)
			frag_presence = [Fraction(x,total_presence) for x in frag_presence]
	fw.close()

########################################################################################

main()