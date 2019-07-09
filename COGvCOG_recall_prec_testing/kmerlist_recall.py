import argparse
from itertools import groupby
from fractions import Fraction

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


def fasta2dict(file_path):
    seqdict = {}

    fiter = fasta_iter(file_path)
    for ff in fiter:
        header = ff[0]
        sequence = ff[1]
        seqdict[header] = sequence
    return seqdict


def replist2list(replist_path):
    replist = []
    with open(replist_path) as opened_file:
        for line in opened_file:
            if not line:
                continue
            kmer = line.strip()
            replist.append(kmer)
    return replist

#####################################################

def main():
    parser = argparse.ArgumentParser(description="Calculates the recall of a given list of kmers on a set of sequences")
    parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required=True)
    parser.add_argument("-r","--replist_file", help="A csv file containing a representative list of kmers",required=True)
    parser.add_argument("-s","--stats_out", help="Output file containing data about the proportion of sequences in the fasta file covered by the kmers in the representative list",required=False,default='stats')
    parser.add_argument("-m","--sum_out", help="Output file detailing the coverage of each of the sequences in the fasta file",required=False,default='summary')
    args = parser.parse_args()

    fastadict = fasta2dict(args.fasta_file)
    replist = replist2list(args.replist_file)

    stats = open(args.stats_out + '.csv','w')
    stats.write("Total Sequences" + "," + "Covered Sequences" + "," + "Recall" + "\n")

    summary = open(args.sum_out + '.csv','w')
    summary.write("Identifier" + "," + "Covered?" + "," + "Coverage" + "\n")

    numseq = len(fastadict)
    seq_covered = 0
    tot_coverage = 0
    for sequence in fastadict:
        coverage = 0
        for kmer in replist:
            if kmer in fastadict[sequence]:
                coverage += 1
                tot_coverage += 1
        if coverage != 0:
            seq_covered += 1
            summary.write(str(sequence) + "," + "True" + "," + str(coverage) + "\n")
        elif coverage == 0:
            summary.write(str(sequence) + "," + "False" + "," + str(coverage) + "\n")
    
    print("There were " + str(seq_covered) + " sequences covered, a total of " + str(tot_coverage) + " times")
    avg_coverage = float(Fraction(tot_coverage,seq_covered))
    recall = (100*seq_covered)/numseq
    stats.write(str(numseq) + "," + str(seq_covered) + "," + str(recall)[:5] + "\n")

    stats.close()
    summary.close()
    print("\nThe recall was " + str(recall)[:5] + "%")
    print("Of the sequences that were covered, they were covered " + str(avg_coverage) + " times, on average.")

#####################################################

if __name__ == '__main__':
    main()
    
