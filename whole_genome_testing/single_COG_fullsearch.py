import argparse
import os


def main():
	parser = argparse.ArgumentParser(description='From a representative oligo list, determines the location of these kmers in a set of whole genomes')
	parser.add_argument('-l','--oligo_list',help='A file containing representative oligos',required=True)
	parser.add_argument('-d','--whole_genomes',help='A directory containing the genomes you are interested in (.fna) and their annotations (.gff)',required=True)
	parser.add_argument('-o','--output',help='Output report with data about all the matches for input gene list',required=True)
	args = parser.parse_args()

	fw = open(args.output,'w')

	for file_name in os.listdir(args.whole_genomes):
		if file_name.endswith('.fna'):
			prefix = file_name[:-4]
			gff_file = prefix + '.gff'

			os.system('python find_oligo_in_genome.py -l ' + args.oligo_list + ' -f ' + args.whole_genomes + '/' + file_name + ' -g ' + args.whole_genomes + '/' + gff_file + ' -o ' + prefix + '.report')

			print('\n' + prefix + '.report')

			for line in open(prefix + '.report'):
				fw.write(line)
			fw.write('\n')
			os.system('rm ' + prefix + '.report')
	fw.close()

main()
