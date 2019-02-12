for p in range(0,10):
	file_name = "covjobp" + str(p) + ".sh"
	job = open(file_name, 'w')
	job.write('#!/bin/sh\npython gen_kmer_list_w_cov.py -f ../../data/recall_partitions/trial_' + str(p) + '_train.fasta -o trial_' + str(p) + '_kcovlist.csv')
	job.close()