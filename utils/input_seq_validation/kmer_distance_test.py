from itertools import groupby
import argparse
import copy


def kmerize_seq(seq,k_size):
	kmers = set()

	for start_loc in range(0, len(seq)-k_size+1):
		kmers.add(seq[start_loc:start_loc+k_size])
	
	return kmers


def fasta_iter(file_path):
	fh = open(file_path)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == '>'))
	for header in faiter:
		seq_name = header.next()[1:].strip()
		sequence = ''.join(s.strip() for s in faiter.next())
		yield seq_name, sequence


def gen_kmer_dict(file_path, k_size):
	kmers_by_seq = {}
	seqs_by_seq = {}

	fiter = fasta_iter(file_path)
	for ff in fiter:
		kmers_by_seq[ff[0]] = kmerize_seq(ff[1], k_size)
		seqs_by_seq[ff[0]] = ff[1]

	return kmers_by_seq, seqs_by_seq


def mean(items):
	return float(sum(items))/len(items)


def std_dev(items):
	return sum([(item - mean(items))**2 for item in items])/len(items)


def get_median(dist_dict):
	total = 0
	for distance in dist_dict:
		total+= dist_dict[distance]

	if total % 2 == 0:
		counter = 0
		dist_dict_copy = copy.deepcopy(dist_dict)
		while True:
			min_dist = None
			for distance in dist_dict_copy:
				if min_dist == None or float(distance) < float(min_dist):
					min_dist = distance
			counter += dist_dict_copy[min_dist]
			if counter == int(total/2):
				lower_value = float(min_dist)
				del dist_dict_copy[min_dist]
				upper_value = None
				for distance in dist_dict_copy:
					if upper_value == None or float(distance) < float(upper_value):
						upper_value = distance
				upper_value = float(upper_value)
				return float(upper_value + lower_value)/2
			elif counter > int(total/2):
				return float(min_dist)
			else:
				del dist_dict_copy[min_dist]

	elif total % 2 != 0:
		counter = 0
		dist_dict_copy = copy.deepcopy(dist_dict)
		while True:
			min_dist = None
			for distance in dist_dict_copy:
				if min_dist == None or float(distance) < float(min_dist):
					min_dist = distance
			counter += dist_dict_copy[min_dist]
			if counter >= int(total/2) + 1:
				return float(min_dist)
			else:
				del dist_dict_copy[min_dist]


def get_range(dist_dict, k_mult):
	middle = get_median(dist_dict)
	print('Median is %s' %(str(middle)))
	bottom_half = {}
	top_half = {}

	for distance in dist_dict:
		if float(distance) >= middle:
			top_half[distance] = dist_dict[distance]
		if float(distance) <= middle:
			bottom_half[distance] = dist_dict[distance]

	q1 = get_median(bottom_half)
	q3 = get_median(top_half)
	iqr = q3 - q1

	print('Minimum value is %s' %(str(q1-(iqr*k_mult))))
	print('Maximum value is %s' %(str(q3+(iqr*k_mult))))

	return q1-(iqr*k_mult), q3+(iqr*k_mult)

#############################################################################

def main():
	parser = argparse.ArgumentParser(description="Determines kmer distance between sequences in input file")
	parser.add_argument("-f","--fna_file", help="A FNA file of sequences",required=True)
	parser.add_argument("-p","--out_prefix", help="Prefix for the output files. A 'prefix'_true.fna file and a 'prefix'_false.fna file will be created",required=True)
	parser.add_argument("-k","--kmer_len", help="The integer length of the kmers to be assessed",required=False,default='9')
	args = parser.parse_args()

	seqs_as_kmers, seqs = gen_kmer_dict(args.fna_file, int(args.kmer_len))

	seq_distances = {}
	all_distances = {}

	for seq_id1 in seqs_as_kmers:
		print('Analyzing sequence: ' + seq_id1)
		seq_distances.setdefault(seq_id1,0)
		for seq_id2 in seqs_as_kmers:
			if seq_id2 != seq_id1:
				seq1_kmers = seqs_as_kmers[seq_id1]
				seq2_kmers = seqs_as_kmers[seq_id2]
				distance = float(len(seq1_kmers.intersection(seq2_kmers)))
				seq_distances[seq_id1] += distance
				all_distances.setdefault(str(distance),0)
				all_distances[str(distance)] += 1
				
	min_dist, max_dist = get_range(all_distances, 0.1)

	fto = open(args.out_prefix + '_true.fna','w')
	ffo = open(args.out_prefix + '_false.fna','w')

	removal_counts = 0

	for seq_id in seq_distances:
		mean_dist = seq_distances[seq_id]/(len(seq_distances)-1)
		if mean_dist >= min_dist:
			"""print('Sequence: ' + seq_id + ' had average identity of ' + str(mean_dist) + ' and was identified as a true variant of this sequence')"""
			fto.write('>' + seq_id + '\n' + seqs[seq_id] + '\n')
		elif mean_dist < min_dist:
			print('Sequence: ' + seq_id + ' had average identity of ' + str(mean_dist) + ' and was identified as a false variant of this sequence')
			ffo.write('>' + seq_id + '\n' + seqs[seq_id] + '\n')
			removal_counts += 1

	fto.close()
	ffo.close()

	print(str(removal_counts) + ' sequences were removed from the input file.')

	
#############################################################################

main()