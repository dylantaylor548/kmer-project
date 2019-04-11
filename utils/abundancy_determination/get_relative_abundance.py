import argparse
from itertools import groupby
from math import floor


def FASTQ_iter(fastq_name):
	fq = open(fastq_name)
	fqiter = (x[1] for x in groupby(fq, lambda line: line[0] == "@"))
	for header in fqiter:
		header = header.__next__()[1:].strip()
		stuff = (''.join(s.strip() for s in fqiter.__next__())).split('+')
		read = stuff[0]
		quality = stuff[1]
		yield header, read, quality


def fastq2dict(file_path):
	readsdict = {}
	read_qual = {}
	fiter = FASTQ_iter(file_path)
	for ff in fiter:
		header = ff[0]
		read = ff[1]
		quality = ff[2]
		readsdict[header] = read
		read_qual[header] = quality
	return readsdict, read_qual


def fasta_iter(fasta_name):
	fh = open(fasta_name)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		header = header.__next__()[1:].strip()
		seq = "".join(s.strip() for s in faiter.__next__())
		yield header, seq


def fasta2dict(file_path):
	seqdict = {}
	fiter = fasta_iter(file_path)
	for ff in fiter:
		header = ff[0]
		sequence = ff[1]
		seqdict[header] = sequence
	return seqdict


def chip2oligos(file_path):
	chipcov = {}
	chip = open(file_path)
	for line in chip:
		line = line.strip().split(',')
		oligo = line[0]
		chip_coverage = int(line[1])
		chipcov[oligo] = chip_coverage
	return chipcov


def gen_kmerdict(seqdict,oligos):
	kmerdict = {}
	for oligo in oligos:
		seqs = []
		for seq in seqdict:
			sequence = seqdict[seq]
			if oligo in sequence:
				seqs.append(seq)
		kmerdict[oligo] = seqs
	return kmerdict


def condense_read_info(reads_in_seqs_dict):
	seq_presence = {}
	for oligo in reads_in_seqs_dict:
		seq_presence[oligo] = {}
		for seq in reads_in_seqs_dict[oligo]:
			total_pres = 0
			for read in reads_in_seqs_dict[oligo][seq]:
				total_pres += reads_in_seqs_dict[oligo][seq][read]
			seq_presence[oligo][seq] = total_pres
	return seq_presence


##############################################################################################################################################################################################

def main():
	"""parser = argparse.ArgumentParser(description='Given a set of reads from an oligo chip, returns a report regarding the relative abundance of strains in the sample that yielded the reads')
	parser.add_argument('-q','--fastq_reads',help='A FASTQ file containing the isolated reads',required=True)
	parser.add_argument('-f','--fasta_file',help='A .fasta file containing all known diversity of your marker gene of interest',required=True)
	parser.add_argument('-s','--spike_seqs',help='A .fasta file containing the sequences used to spike your chip assay',required=True)
	parser.add_argument('-c','--chip',help='A .csv file containing your generated oligos AND their presence on the chip',required=True)
	args = parser.parse_args()"""

	fasta_file = 'C:/Users/Dylan/Desktop/Pop_Lab/kmer_project/data/TIGR02012/TIGR02012_100.fasta'
	fastq_reads = 'C:/Users/Dylan/Desktop/TIGR_sample.fq'
	chip = 'C:/Users/Dylan/Desktop/tmp_chip.csv'
	"""spike_seqs = 'C:/Users/Dylan/Desktop/Pop_Lab/kmer_project/data/COG0088/COG0088_10000.fasta'"""

	seqdict = fasta2dict(fasta_file)
	readsdict, read_qual = fastq2dict(fastq_reads)
	"""spikes = fasta2dict(spike_seqs)"""
	
	known_seqs = {}
	for seq in seqdict:
		known_seqs[seq] = seqdict[seq]
	"""for seq in spikes:
		known_seqs[seq] = spikes[seq]"""

	chipcovdict = chip2oligos(chip)

	kmerdict = gen_kmerdict(known_seqs,chipcovdict)

	seqs_bound_oligos = {}
	for oligo in chipcovdict:
		seqs_bound_oligos[oligo] = {}

	overassigned_reads = {}
	for oligo in chipcovdict:
		overassigned_reads[oligo] = {}

	for read in readsdict:
		read_seq = readsdict[read]
		oligo_matches = []
		for oligo in chipcovdict:
			if oligo in read_seq:
				oligo_matches.append(oligo)
		if len(oligo_matches) == 1:
			oligo_match = oligo_matches[0]
			seq_matches = []
			for seq in kmerdict[oligo_match]:
				if read_seq in known_seqs[seq]:
					seq_matches.append(seq)
			if len(seq_matches) == 1:
				seq_match = seq_matches[0]
				seqs_bound_oligos[oligo_match].setdefault(seq_match,{})
				seqs_bound_oligos[oligo_match][seq_match].setdefault(read_seq,0)
				seqs_bound_oligos[oligo_match][seq_match][read_seq] += 1
			elif len(seq_matches) > 1:
				overassigned_reads[oligo_match].setdefault(read_seq,0)
				overassigned_reads[oligo_match][read_seq] += 1

	avg_read_cov = {}
	for oligo in chipcovdict:
		avg_read_cov[oligo] = {}
		for seq in seqs_bound_oligos[oligo]:
			tot_reads = 0
			counter = 0
			for read in seqs_bound_oligos[oligo][seq]:
				tot_reads += seqs_bound_oligos[oligo][seq][read]
				counter += 1
			avg_read_cov[oligo][seq] = float(tot_reads/counter)

	for oligo in overassigned_reads:
		for read in overassigned_reads[oligo]:
			avail_reads = overassigned_reads[oligo][read]
			needed_reads = 0
			for seq in kmerdict[oligo]:
				if read in known_seqs[seq]:
					if seq in avg_read_cov[oligo]:
						needed_reads += avg_read_cov[oligo][seq]
			if needed_reads <= avail_reads:
				for seq in kmerdict[oligo]:
					if read in known_seqs[seq]:
						if seq in avg_read_cov[oligo]:
							seqs_bound_oligos[oligo][seq][read] = floor(avg_read_cov[oligo][seq])
							overassigned_reads[oligo][read] -= floor(avg_read_cov[oligo][seq])
			elif needed_reads > avail_reads:
				for seq in kmerdict[oligo]:
					if read in known_seqs[seq]:
						if seq in avg_read_cov[oligo]:
							seqs_bound_oligos[oligo][seq][read] = floor(avg_read_cov[oligo][seq]*(avail_reads/needed_reads))
							overassigned_reads[oligo][read] -= floor(avg_read_cov[oligo][seq]*(avail_reads/needed_reads))
	# Once this step is done, we should be left with a dictionary of oligos and reads containing those oligos that have not been attributed to any sequences.
	# Now what we want to do is to assign these reads evenly to those sequences that have not had any reads assigned to them so far

	for oligo in overassigned_reads:
		for read in overassigned_reads[oligo]:
			avail_reads = overassigned_reads[oligo][read]
			seqs_left = []
			for seq in kmerdict[oligo]:
				if read in known_seqs[seq]:
					if seq not in avg_read_cov[oligo]:
						seqs_left.append(seq)
			if len(seqs_left) != 0:
				quant = floor(avail_reads/len(seqs_left))
				for seq in seqs_left:
					seqs_bound_oligos[oligo].setdefault(seq,{})
					seqs_bound_oligos[oligo][seq].setdefault(read,0)
					seqs_bound_oligos[oligo][seq][read] += quant
					overassigned_reads[oligo][read] -= quant

	seq_presence = condense_read_info(seqs_bound_oligos)

	print(seq_presence)
						



##############################################################################################################################################################################################

if __name__ == '__main__':
	main()