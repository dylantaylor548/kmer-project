import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
import os
from fractions import Fraction


def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
        yield header, seq


def fasta2dict(file_path):
    seqdict = {}

    fiter = fasta_iter(file_path)
    for ff in fiter:
        header = ff[0]
        sequence = ff[1]
        seqdict[header] = sequence
    return seqdict


kmerlist_dir = 'C:/Users/Dylan/Desktop/Pop_Lab/kmer_project/recall_testing/track_kmer_cov'
testingfastadir = 'C:/Users/Dylan/Desktop/Pop_Lab/kmer_project/data/recall_partitions'
out_dir = 'C:/Users/Dylan/Desktop'

Error = True

fig1, ax1 = plt.subplots(nrows=1,ncols=1)
maxy = 0

covtracker = {}

for filename in os.listdir(kmerlist_dir):
    if filename.endswith('.csv'):
        kmercov = []
        numcovd = []
        prefix = filename[0:8]
        kmercovdict = {}
        for line in open(kmerlist_dir + '/' + filename):
            elements = line.split(',')
            kmercovdict[elements[0]] = int(elements[1])
        testfastafile = testingfastadir + '/' + prefix + 'test.fasta'
        fastadict = fasta2dict(testfastafile)
        for kmer in kmercovdict:
            # The below line is used to selectively examine just those kmers with at or below coverage of 50 when they were selected
            """if kmercovdict[kmer] <= 50:"""
            if True:
                kmercov.append(kmercovdict[kmer])
                if not str(kmercovdict[kmer]) in covtracker:
                    covtracker[str(kmercovdict[kmer])] = []
                count = 0
                for seq in fastadict:
                    if kmer in fastadict[seq]:
                        count += 1
                numcovd.append(count)
                covtracker[str(kmercovdict[kmer])].append(count)
        if not Error:
	        if max(numcovd) > maxy:
	            maxy = max(numcovd)
	        ax1.plot(kmercov,numcovd,"o",color='lightblue',markersize=2)

if Error:
	kcovlist = []
	meanlist = []
	errlist = []
	for kcov in covtracker:
		kcovlist.append(int(kcov))
		meanlist.append(np.mean(covtracker[kcov]))
		errlist.append(np.std(covtracker[kcov]))
		if np.mean(covtracker[kcov]) + np.std(covtracker[kcov]) > maxy:
			maxy = np.mean(covtracker[kcov]) + np.std(covtracker[kcov])
	"""ax1.errorbar(kcovlist,meanlist,errlist,marker="o",linestyle='None',capsize=3,color='blue',markersize=2)"""
	ax1.plot(kcovlist,meanlist,marker="o",linestyle='None',color='blue',markersize=2)

ax1.set_title('Kmer quality based on coverage')
ax1.set_xlabel("kmer coverage")
ax1.set_ylabel("Num of sequences covered by kmer")
ax1.axhline(y=0,xmin=0,xmax=1,color='red',linewidth=0.5)
plt.yticks(np.arange(0, maxy+2, 20))
ax1.grid(b=True, which='both', color='lightgray', linestyle='-')
plt.minorticks_on()
plt.tight_layout()
plt.savefig(out_dir + "/cov_based_usefulness.png")
plt.show()