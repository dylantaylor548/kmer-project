import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction

def graphrecall(recalldict,title,out_dir):
	fig1, ax1 = plt.subplots(nrows=1,ncols=1)

	klens = list(recalldict.keys())
	test = klens[0]
	klensint = [int(x) for x in klens]
	
	for i in range(0,len(recalldict[test])):
		recall = []
		for klen in klens:
			recall.append(recalldict[klen][i])
		ax1.plot(klensint,recall,"o",color='lightgray')

	avgrecalls = []
	for klen in klens:
		tot_recall = 0
		count = 0
		for recall in recalldict[klen]:
			tot_recall += recall
			count += 1
		avgrecall = tot_recall/count
		avgrecalls.append(avgrecall)

	ax1.plot(klensint,avgrecalls,"o",color='magenta')

	ax1.set_title(title)
	ax1.set_xlabel("kmer length")
	ax1.set_ylabel("Recall (%)")
	plt.tight_layout()
	plt.savefig(out_dir + "/k_based_recall.png")
	plt.show()


################################################################

partition_dir = 'C:/Users/Dylan/Desktop/Pop_Lab/kmer_project/data/recall_partitions'
out_dir = 'C:/Users/Dylan/Desktop'
graph_name = '8000/2000 Training/Test Sequences'


recalldict = {}
covdict = {}

for folder in os.listdir(partition_dir):
	if folder.startswith('recall_stats_k'):
		klen = folder[14:]
		recalldict[klen] = []
		folderpath = partition_dir + '/' + folder
		for file in os.listdir(folderpath):
			if file.endswith('stats.csv'):
				filepath = folderpath + '/' + file
				with open(filepath) as opened_file:
					for line in opened_file:
						if line.startswith('T'):
							continue
						else:
							linelist = line.split(',')
							tot_seq = int(linelist[0])
							cov_seq = int(linelist[1])
							recall = 100*Fraction(cov_seq,tot_seq)
							recalldict[klen].append(recall)

graphrecall(recalldict,graph_name,out_dir)

    
################################################################
