import argparse
from itertools import groupby

def fasta2dict(fasta_path):
    seqdict = {}
    with open(fasta_path) as opened_file:
        for line in opened_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                sequence_name = line[1:]
                continue
            sequence = line
            seqdict[sequence_name] = sequence
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
    parser.add_argument("-s","--stats_out", help="Output file containing data about the proportion of sequences in the fasta file covered by the kmers in the representative list",required=True)
    parser.add_argument("-m","--sum_out", help="Output file detailing the coverage of each of the sequences in the fasta file",required=True)
    args = parser.parse_args()

    fastadict = fasta2dict(args.fasta_file)
    replist = replist2list(args.replist_file)

    stats = open(args.stats_out,'w')
    stats.write("Total Sequences" + "," + "Covered Sequences" + "," + "Recall" + "\n")

    summary = open(args.sum_out,'w')
    summary.write("Identifier" + "," + "Covered?" + "," + "Coverage" + "\n")

    numseq = len(fastadict)
    seq_covered = 0
    tot_coverage = 0
    for sequence in fastadict.keys():
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

    avg_coverage = tot_coverage/seq_covered
    recall = (100*seq_covered)/numseq
    stats.write(str(numseq) + "," + str(seq_covered) + "," + str(recall)[:5] + "\n")

    stats.close()
    summary.close()
    print("\nThe recall was " + str(recall)[:5] + "%")
    print("Of the sequences that were covered, they were covered " + str(avg_coverage)[:3] + " times, on average.")

#####################################################

if __name__ == '__main__':
    main()
    
