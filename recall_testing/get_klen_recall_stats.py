import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction

def graphrecall(recalldict,title,out_dir):
	fig1, ax1 = plt.subplots(nrows=1,ncols=1)

	klens = recalldict.keys()
	klensint = [int(x) for x in klens]
	
	for i in range(0,len(recalldict[klens[0]])):
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
		avrrecalls.append(avgrecall)

	ax1.plot(klensint,avgrecalls,"o",color='magenta')

	ax1.set_title(title)
	ax1.set_xlabel("kmer length")
	ax1.set_ylabel("Recall (%)")
	plt.tight_layout()
	plt.savefig(out_dir + "/k_based_recall.png")
	plt.show()


################################################################


def main():
    parser = argparse.ArgumentParser(description="Finds representative kmers for a set of sequences")
    parser.add_argument("-pd","--partition_dir", help="A directory containing the recall statistics for different k-lengths",required=True)
    parser.add_argument("-o","--out_dir",help="The directory that the figures will be outputted to",required=True)
    parser.add_argument("-n","--graph_name",help="The desired title of your figures",required=False,default='8000/2000 Training/Test Sequences')
    args = parser.parse_args()

    recalldict = {}
    covdict = {}

    for folder in os.listdir(args.partition_dir):
    	if folder.startswith('recall_stats_k'):
    		klen = folder[14:]
    		recalldict[klen] = []
    		folderpath = args.partition_dir + '/' + folder
    		for file in os.listdir(folderpath):
    			if file.endswith('stats.csv'):
    				filepath = folderpath + '/' + file
    				with open(filepath) as opened_file:
    					for line in opened_file:
    						if line.startswith('T'):
    							continue
    						else:
    							linelist = line.split(',')
    							tot_seq = linelist[0]
    							cov_seq = linelist[1]
    							recall = 100*Fraction(cov_seq,tot_seq)
    							recalldict[klen].append(recall)

    graphrecall(recalldict,args.graph_name,args.out_dir)

    
################################################################



if __name__ == '__main__':
    main()