import argparse
from itertools import groupby


# Returns a list of oligos from a file containing those oligos
def oligos2list(replist_file):
	oligos = []
	fr = open(replist_file)
	for line in fr:
		oligo = line.strip()
		oligos.append(oligo)
	fr.close()
	return oligos


def get_rev_complement(sequence):
	rev_complement = ''
	alphabet = ['A','C','T','G']
	for base in sequence:
		if base in alphabet:
			rev_complement = alphabet[(alphabet.index(base)-2)%4] + rev_complement
		elif base == 'N':
			rev_complement = 'N' + rev_complement
	return rev_complement


#Reads a fasta file
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


##############################################################################################################################################

def main():
	parser = argparse.ArgumentParser(description='From a representative oligo list, determines the location of these kmers in a whole genome')
	parser.add_argument('-l','--oligo_list',help='A file containing representative oligos',required=True)
	parser.add_argument('-f','--whole_fna',help='A .fna or .fasta (or similar file) containing the sequence of a whole genome',required=True)
	parser.add_argument('-g','--whole_gff',help='A .gff file containing annotations of the input whole genome file',required=False)
	parser.add_argument('-o','--output',help='An output file containing information about which elements of the whole genome were hit by the kmers in the list')
	args = parser.parse_args()

	genome_seqs = fna_read(args.whole_fna)
	oligos = oligos2list(args.oligo_list)
	oligo_len = len(oligos[0])
	oligos = set(oligos)

	locations = {}
	for seq_id in genome_seqs:
		locations.setdefault(seq_id,{})
		genome_seq = genome_seqs[seq_id]
		for i in range(0,len(genome_seq)-oligo_len+1):
			frame = genome_seq[i:i+oligo_len]
			if frame in oligos:
				locations[seq_id].setdefault(frame,[])
				locations[seq_id][frame].append(i)
				print("Oligo found at location: " + str(i) + " in sequence: " + seq_id)
			rev_comp =  get_rev_complement(frame)
			if rev_comp in oligos:
				locations[seq_id].setdefault('*' + rev_comp,[])
				locations[seq_id]['*' + rev_comp].append(i)
				print("Reverse complement found at location: " + str(i) + " in sequence: " + seq_id)



	matches = {}
	for seq_id in locations:
		matches.setdefault(seq_id,{})
		for match in locations[seq_id]:
			for loc in locations[seq_id][match]:
				check = False
				for line in open(args.whole_gff):
					if not line.startswith('#'):
						line_elems = line.split('	')
						if line_elems[0] == seq_id:
							start_loc = int(line_elems[3])
							end_loc = int(line_elems[4])
							if (loc + 1) in range(start_loc,end_loc+1):
								feature = line_elems[2]
								if feature != 'region':
									feature_seq = genome_seqs[seq_id][start_loc-1:end_loc]
									check = True
									matches[seq_id].setdefault(feature,{})
									feat_props = line_elems[8].split(';')
									feat_name = ''
									if feature == 'gene':
										for prop in feat_props:
											if prop.startswith('ID='):
												feat_name += prop[3:] + '	'
											elif prop.startswith('Name='):
												feat_name += prop[5:]
												break
									elif feature == 'CDS':
										for prop in feat_props:
											if prop.startswith('Parent='):
												feat_name += prop[7:] + '	'
											elif prop.startswith('product='):
												feat_name += prop[8:]
												break
									else:
										for prop in feat_props:
											if prop.startswith('ID='):
												feat_name = prop[3:]
												break

									matches[seq_id][feature].setdefault(feat_name,[[],feature_seq])
									matches[seq_id][feature][feat_name][0].append(match)
							
				if check == False:
					print('Found one in intergenic space')


	seq_identities = {}
	for line in open(args.whole_gff):
		if not line.startswith('#'):
			line_elems = line.split('	')
			if line_elems[2] == 'region':
				feat_props = line_elems[8].split(';')
				for prop in feat_props:
					if prop.startswith('genome='):
						seq_identities[line_elems[0]] = prop[7:]


	fo = open(args.output,'w')

	for line in open(args.whole_gff):
		if line.startswith('#!ge'):
			fo.write(line)
	fo.write('##query-gene-list ' + args.oligo_list[:-5] + '\n')

	
	for seq_id in matches:
		fo.write(seq_id + '	' + seq_identities[seq_id] + '\n')
		if matches[seq_id] != {}:
			for feature in matches[seq_id]:
				for feat_name in matches[seq_id][feature]:
					fo.write('	' + feature + '	' + feat_name + '	' + str(len(matches[seq_id][feature][feat_name][0])) + ' oligo matches' + '	' + matches[seq_id][feature][feat_name][1] + '\n')
		elif matches[seq_id] == {}:
			fo.write('	None\n')

	fo.close()



##############################################################################################################################################

main()