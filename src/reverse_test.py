from itertools import groupby
import argparse
import time

def kmerize(dna,kmer_size):
    kmers = set()
    kmers_filtered = []
    dna = dna.upper()
    if (kmer_size <= len(dna) and kmer_size >= 1):
        for start in range(0,len(dna)-kmer_size+1,1):
            kmer = dna[start:start+kmer_size]
            kmers.add(kmer)
        kmers_filtered = list(kmers)
        return kmers_filtered

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def kmerize_directory(file_path,kmer_size):
    seqdict = {}
    print("Generating " + str(kmer_size)+ "-mers\n...")

    fiter = fasta_iter(file_path)
    for ff in fiter:
        header = ff[0]
        sequence = ff[1]
        kmers_in_dna = kmerize(sequence,kmer_size)
        seqdict[header] = dict((x, 1) for x in kmers_in_dna)
    print("\nAll " + str(kmer_size)+ "-mers generated!\n")
    return seqdict

def reverse_dict(seqdict):
    rev_dict = {}
    for seq in seqdict:
        for kmer in seqdict[seq]:
            if kmer in rev_dict:
                rev_dict[kmer][seq] = 1
            else:
                rev_dict[kmer] = {}
                rev_dict[kmer][seq] = 1
    return rev_dict

def gen_leastcovdict(kmerdict,seqcovdict):
    leastcovdict = {}
    for kmer in kmerdict:
        leastcov = 0
        leastcount = 0
        for seq in kmerdict[kmer]:
            if leastcov == 0:
                leastcov = seqcovdict[seq]
                leastcount = 1
            elif leastcov < seqcovdict[seq]:
                leastcov = seqcovdict[seq]
                leastcount = 1
            elif leastcov == seqcovdict[seq]:
                leastcount +=1
        leastcovdict[kmer] = [leastcov,leastcount]
    return leastcovdict

def checkcomplete(kmercovdict,coverage):
    check = True
    for kmer in kmercovdict:
        if kmercovdict[kmer] > coverage:
            check = False
            break
    return check

def findkmin(kmercovdict,kmerdict):
    kmin = None
    for kmer in kmercovdict:
        if kmin == None:
            kmin = kmer
        elif kmercovdict[kmer][0] < kmercovdict[kmin][0]:
            kmin = kmer
        elif kmercovdict[kmer][0] == kmercovdict[kmin][0]:
            if kmercovdict[kmer][1] > kmercovdict[kmin][1]:
                kmin = kmer
            elif kmercovdict[kmer][1] == kmercovdict[kmin][1]:
                if len(kmerdict[kmer]) < len(kmerdict[kmin]):
                    kmin = kmer
                elif len(kmerdict[kmer]) == len(kmerdict[kmin]):
                    if kmer < kmin:
                        kmin = kmer
    return kmin

def reduce_kmerlist(seqdict,coverage):
    kmerdict = reverse_dict(seqdict)
    print("There are " + str(len(kmerdict)) + " kmers in this file. Whittling it down...")
    
    seqcovdict = {}
    for seq in seqdict:
        seqcovdict[seq] = len(seqdict[seq])

    kmercovdict = gen_leastcovdict(kmerdict,seqcovdict)

    remcount = 0
    while not checkcomplete(kmercovdict,coverage):
    	if remcount != 0:
	        if (remcount%10) == 0:
	    	    print("We've removed " + str(remcount) + " kmers")

        kmin = findkmin(kmercovdict,kmerdict)
        for seq in kmerdict[kmin]:
            seqcovdict[seq] -= 1
        del kmerdict[kmin]
        remcount += 1
        kmercovdict = gen_leastcovdict(kmerdict,seqcovdict)

    return kmerdict

################################################################ vvv This is the part that does the code vvv

def main():
    parser = argparse.ArgumentParser(description="Finds representative kmers for a set of sequences")
    parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required=True)
    parser.add_argument("-o","--out_file", help="Output file with representative kmers",required=True)
    parser.add_argument("-c","--cutoff_value", help="cutoff value - it is the number of times each sequence needs to be covered",required=False,default='1')
    parser.add_argument("-k","--kmer_len", help="kmer length",required=False,default='21')
    args = parser.parse_args()

    start_time = time.time()

    seqdict = kmerize_directory(args.fasta_file,int(args.kmer_len))

    reduceddict = reduce_kmerlist(seqdict,int(args.cutoff_value))

    print("There are " + str(len(reduceddict)) + " " + str(args.kmer_len) + "-mers in the representative list.")
    
    end_time = time.time()
    runtime = end_time - start_time
    seconds = int(runtime % 60)
    minutes = int(runtime/60)
    print("It took " + str(minutes) + " minutes and " + str(seconds) + " seconds to generate a list from " + str(len(seqdict)) + " sequences.")

################################################################ ^^^ This is the part that does the code ^^^



if __name__ == '__main__':
    main()
    



