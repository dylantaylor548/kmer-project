import argparse
import random


def gen_kmerlist(file_path):
	oligos = []

	with open(file_path) as opened_file:
		for line in opened_file:
			line = line.strip()
			oligos.append(line)

	return oligos


def gen_spike_tag(oligos):
	tag = ''
	i = 0
	while i < len(oligos[0]):
		base_presence = {'A':0,'C':0,'G':0,'T':0}
		for oligo in oligos:
			base = oligo[i]
			base_presence[base] += 1
		min_pres = None
		min_base = None
		for base in base_presence:
			if min_pres == None:
				min_base = base
				min_pres = base_presence[base]
			elif base_presence[base] < min_pres:
				min_base = base
				min_pres = base_presence[base]
			elif base_presence[base] == min_pres:
				if random.choice([True,False]):
					min_base = base
					min_pres = base_presence[base]
		tag = tag + min_base
		i += 1
	if tag not in oligos:
		return tag
	else:
		while True:
			print('tag in oligo list, finding new tag...')
			alphabet = ['A','C','G','T']
			new_tag = ''
			for i in range(0,len(oligos[0])):
				new_tag = new_tag + random.choice(alphabet)
			if new_tag not in oligos:
				return new_tag


# Generates a list of spike sequences that do not exceed a given length and that have exactly one oligo in any fragment of a given length.
# Each spike sequence contains the spike tag oligo (same in each spike sequence but unique to the spike sequences) surrounded on each side
# by 5 random bases allowing distinction between spike sequences on fragments containing the spike tag oligo.
def gen_spikes(oligos, frag_len, seq_len, spike_tag):
	alphabet = ['A','T','C','G']
	spikes = []
	ids = []
	
	i = 0
	running_len = frag_len
	while True:
		spike_id = ''
		while True:
			spike_id = random.choice(alphabet) + random.choice(alphabet) + random.choice(alphabet) + random.choice(alphabet) + random.choice(alphabet) + spike_tag + random.choice(alphabet) + random.choice(alphabet) + random.choice(alphabet) + random.choice(alphabet) + random.choice(alphabet)
			if spike_id not in ids:
				ids += spike_id
				break
		first_base = spike_id[0]
		spike = ''
		if first_base == 'A' or first_base == 'T':
			spike = (frag_len-len(oligos[i])-5)*'G'
		elif first_base == 'G' or first_base == 'C':
			spike = (frag_len-len(oligos[i])-5)*'A'
		spike = spike + spike_id
		current_end = spike[-1:]
		next_start = oligos[i][0]
		useable_bases = [x for x in alphabet if x != current_end and x != next_start]
		spike = spike + (frag_len-len(oligos[i])-5)*random.choice(useable_bases)

		while len(spike) + frag_len <= seq_len:
			oligo = oligos[i]
			spike = spike + oligo
			last_base = oligo[-1:]
			next_base = ''
			if i+1 < len(oligos):
				next_base = oligos[i+1][0]
			new_useable_bases = [x for x in alphabet if x != last_base and x != next_base]
			spike = spike + (frag_len-len(oligos[i]))*random.choice(new_useable_bases)
			i += 1
			if i == len(oligos):
				spikes.append(spike)
				return spikes
		spikes.append(spike)


################################################################################

def main():
	parser = argparse.ArgumentParser(description='Finds a short set of sequences that collectively comprise all kmers within an input kmer list. These synthetic sequences are used to determine relative abundancies when not all sequences are present in a sample.')
	parser.add_argument("-kl","--kmer_list", help="A .csv file containing the kmer list you are interested in compressing", required=True)
	parser.add_argument("-n","--seq_len", help="Desired maximum length of synthetic sequences", required=True)
	parser.add_argument("-s","--frag_size",help="Estimated size of reads produced by fragmenting DNA sample",required=True)
	parser.add_argument("-o","--out_file", help="A .fasta file to output the synthetic sequences", required=True)
	parser.add_argument("-i","--fasta_id",help='An identification marker for your output sequences',required=True)
	parser.add_argument("-t","--spike_tag",help='A unique oligo sequence added to each spike sequence',required=False,default=None)
	args = parser.parse_args()

	kmerlist = gen_kmerlist(args.kmer_list)
	print('\n\n' + str(len(kmerlist)) + " kmers in list")
	print("\nGenerating spike sequences...")

	spike_tag = None
	if args.spike_tag == None:
		spike_tag = gen_spike_tag(kmerlist)
	else:
		spike_tag = args.spike_tag

	spikes = gen_spikes(kmerlist, int(args.frag_size), int(args.seq_len), spike_tag)
	print("Done.\n\n" + str(len(spikes)) + " spike sequences generated. Spike tag is '" + spike_tag + "'")

	fw = open(args.out_file,'w')
	i = 1
	for spike in spikes:
		seq_id = ('0000'+str(i))[-4:] + "_" + args.fasta_id
		fw.write('>' + seq_id + '\n' + spike + '\n')
		i += 1

	fw.close()


################################################################################

main()