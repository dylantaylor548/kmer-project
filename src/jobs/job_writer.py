for k in range(35,42,2):
	for p in range(0,10):
		file_name = "jobk" + str(k) + "p" + str(p) + ".sh"
		job = open(file_name, 'w')
		job.write('#!/bin/sh\npython ../kmer_project_master.py -f ../../data/recall_partitions/trial_' + str(p) + '_train.fasta -o ../../data/recall_partitions/recall_stats_k' + str(k) + '/trial_' + str(p) + '_replist.csv -k ' + str(k))
		job.close()