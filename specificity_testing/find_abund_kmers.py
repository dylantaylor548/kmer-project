import argparse

def main():
	parser = argparse.ArgumentParser(description='Finds the kmers that hit the most sequences (true and false positives) in a sample')
	parser.add_argument('-r','--report_file',help='A report file containing data about the kmers that hit each sequence',required=True)
	parser.add_argument('-t','--top_results',help='The number of most represented kmers output',required=False,default='10')
	args = parser.parse_args()

	COG_id = args.report_file[:7]
	print(COG_id)

	true_pos_kmers = {}
	false_pos_kmers = {}

	with open(args.report_file) as opened_file:
		for line in opened_file:
			if line.startswith(COG_id):
				kmers = [x.strip() for x in line.split(',')[2:]]
				for kmer in kmers:
					true_pos_kmers.setdefault(kmer,0)
					true_pos_kmers[kmer] +=1
			elif not line.startswith(COG_id):
				kmers = [x.strip() for x in line.split(',')[2:]]
				for kmer in kmers:
					false_pos_kmers.setdefault(kmer,0)
					false_pos_kmers[kmer] +=1

	i = 0
	top_true_kmers = {}
	ttk_list = []
	top_false_kmers = {}
	tfk_list = []
	while i < int(args.top_results):
		
		top_true_kmer = ''
		for kmer in true_pos_kmers:
			if top_true_kmer == '':
				top_true_kmer = kmer
			elif true_pos_kmers[kmer] >= true_pos_kmers[top_true_kmer]:
				top_true_kmer = kmer
				
		top_true_kmers[top_true_kmer] = true_pos_kmers[top_true_kmer]
		ttk_list.append(top_true_kmer)
		del true_pos_kmers[top_true_kmer]

		top_false_kmer = ''
		for kmer in false_pos_kmers:
			if top_false_kmer == '':
				top_false_kmer = kmer
			elif false_pos_kmers[kmer] >= false_pos_kmers[top_false_kmer]:
				top_false_kmer = kmer

		top_false_kmers[top_false_kmer] = false_pos_kmers[top_false_kmer]
		tfk_list.append(top_false_kmer)
		del false_pos_kmers[top_false_kmer]

		i += 1

	print("\n\nThe top represented kmers in captured " + COG_id + " sequences are as follows...\n")
	for kmer in ttk_list:
		print(kmer + " : " + str(top_true_kmers[kmer]))

	print("\n\nThe top represented kmers in captured non-" + COG_id + " sequences are as follows...\n\nKmer" + (" "*(len(ttk_list[0])-4)) + "    non-" + COG_id + "    " + COG_id)
	for kmer in tfk_list:
		if kmer in true_pos_kmers:
			print(kmer + "    " + str(top_false_kmers[kmer]) + (" "*(11 - len(str(top_false_kmers[kmer])))) + "    " + str(true_pos_kmers[kmer]))
		else:
			print(kmer + "    " + str(top_false_kmers[kmer]) + (" "*(11 - len(str(top_false_kmers[kmer])))) + "    0")			



main()