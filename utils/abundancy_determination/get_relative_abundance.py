import argparse
from itertools import groupby


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

	for read in readsdict:
		read_seq = readsdict[read]
		oligo_matches = []
		for oligo in chipcovdict:
			if oligo in read_seq:
				oligo_matches.append(oligo)
		if len(oligo_matches) == 1:
			seq_matches = []
			for seq in kmerdict[oligo_matches[0]]:
				if read_seq in known_seqs[seq]:
					seq_matches.append(seq)
			if len(seq_matches) == 1:
				seqs_bound_oligos[oligo_matches[0]].setdefault(seq_matches[0],0)
				seqs_bound_oligos[oligo_matches[0]][seq_matches[0]] += 1
			elif len(seq_matches) > 1:
				overassigned_reads.setdefault(read_seq,0)
				overassigned_reads[read_seq] += 1
				for seq in seq_matches:
					seqs_bound_oligos[oligo_matches[0]].setdefault(seq,0)

	for read_seq, read_count in overassigned_reads.items():
		for oligo in chipcovdict:
			if oligo in read_seq:
				seq_matches = []
				for seq_id  in seqs_bound_oligos[oligo]:
					sequence = known_seqs[seq_id]
					if read_seq in sequence:
						seq_matches.append(seq_id)
				for seq_id in seq_matches:
					quant = (read_count/len(seq_matches))
					seqs_bound_oligos[oligo][seq_id] += quant
				break


	for oligo in seqs_bound_oligos:
		print(oligo)
		print(seqs_bound_oligos[oligo])
		print('')

##############################################################################################################################################################################################

if __name__ == '__main__':
	main()