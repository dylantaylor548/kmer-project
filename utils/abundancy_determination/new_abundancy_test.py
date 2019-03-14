from math import floor
import random
import copy


# Generates a chip when given a dictionary containing sequence matches for each oligo, as well as a max size (# of oligos) for the output chip
def make_chip(chip_dict,chip_size,setup='proportional'):
	chip = []

	if setup == 'proportional':
		chip_intermediate = []
		for oligo in chip_dict:
			for i in range(0,len(chip_dict[oligo])):
				chip_intermediate.append(oligo)
		chip_mult = floor(chip_size/len(chip_intermediate))
		chip = chip_intermediate * chip_mult
		print("Using " + str(len(chip)) + " spots on the chip")

	elif setup == 'even':
		possible_seqs = set()
		for oligo in chip_dict:
			for seq in chip_dict[oligo]:
				possible_seqs.add(seq)
		possible_seqs = list(possible_seqs)

		chip_mult = floor(chip_size/len(possible_seqs))
		chip = possible_seqs * chip_mult
		print("Using " + str(len(chip)) + " spots on the chip")

	else:
		print("'" + str(setup) + "' is not a valid input parameter for 'setup'")
		return	

	return chip


# Simulates washing an input sample (in the form of a list containing all the sequences in the sample) over an input chip
def simulate_chip_washing(chip_dict,chip,input_seqs):
	seq_abundances = {}
	for oligo in chip_dict:
		for seq in chip_dict[oligo]:
			seq_abundances[seq] = 0

	while chip != []:
		pick_oligo = random.choice(chip)
		pick_seq =  random.choice(input_seqs)

		if pick_seq in chip_dict[pick_oligo]:
			seq_abundances[pick_seq] += 1
			chip.remove(pick_oligo)
			input_seqs.remove(pick_seq)
		else:
			continue

	return seq_abundances


############################################################################################################

chip_dict = {
	'kmer1' : ['seq1','seq2'],
	'kmer2' : ['seq2','seq3','seq4']
}
seq_sample = (['seq1'] * 100000) + (['seq2'] * 100000) + (['seq3'] * 200000) + (['seq4'] * 1000000)
chip_size = 1000
trials = 25

############################################################################################################

f = open('C:/Users/Dylan/Desktop/abundance_test.csv','w')
possible_seqs = list(set(seq_sample))
firstline = ','+ ','.join(possible_seqs) + '\n'
f.write(firstline)

for i in range(1,trials+1):
	print("Doing Trial " + str(i))
	seq_sample_copy = copy.deepcopy(seq_sample)
	chip = make_chip(chip_dict,chip_size)
	abundances = simulate_chip_washing(chip_dict,chip,seq_sample_copy)
	line = 'Trial ' + str(i)
	for seq in possible_seqs:
		line += ',' + str(abundances[seq])
	f.write(line + '\n')

f.close()
