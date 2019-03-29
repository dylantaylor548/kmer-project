from __future__ import print_function
import argparse
from fractions import Fraction
from math import floor

def main():
	parser = argparse.ArgumentParser(description='Creates a chip (a .csv of the oligos in the abundances they will appear on the chip) from an oligo list, and a list of input sequences')
	parser.add_argument("-l","--oligo_list",help='A .csv file containing the oligos generated from your input sequences',required=True)
	parser.add_argument("-f","--fasta_file",help="The .fasta file used to generate the oligo list",required=False,default=None)
	parser.add_argument("-s","--spike_seq",help="A .fasta file of synthesized spike sequences for your sample",required=False,default=None)
	parser.add_argument("-c","--chip_size",help="The number of oligos you will be able to print onto your chip",required=True)
	parser.add_argument("-o","--output_file",help="An output .csv file that will store the abundances of each oligo on your chip")
	args = parser.parse_args()

	print("\n\nCreating chip...",end='')

	oligos = {}
	with open(args.oligo_list,'r') as oligo_file:
		for line in oligo_file:
			line = line.strip()
			oligos[line] = 0
	oligo_file.close()

	if args.fasta_file != None:
		with open(args.fasta_file,'r') as fasta_file:
			for line in fasta_file:
				if not line.startswith('>'):
					for oligo in oligos:
						if oligo in line:
							oligos[oligo] += 1
		fasta_file.close()

	if args.spike_seq != None:
		with open(args.spike_seq,'r') as spike_seq:
			for line in spike_seq:
				if not line.startswith('>'):
					for oligo in oligos:
						if oligo in line:
							oligos[oligo] += 1
		spike_seq.close()

	total_presence = 0
	for oligo in oligos:
		total_presence += oligos[oligo]

	for oligo in oligos:
		oligos[oligo] = int(floor(int(args.chip_size) * Fraction(oligos[oligo],total_presence)))

	total_oligos = 0
	fw = open(args.output_file,'w')
	for oligo in oligos:
		line = oligo + ',' + str(oligos[oligo]) + '\n'
		fw.write(line)
		total_oligos += oligos[oligo]

	print(" done.\n" + str(total_oligos) + " spaces used on chip")

if __name__ == '__main__':
	main()