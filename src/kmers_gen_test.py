from itertools import groupby
import argparse

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
    kmerdict = {}
    print("Generating " + str(kmer_size)+ "-mers\n...")

    fiter = fasta_iter(file_path)
    for ff in fiter:
        header = ff[0]
        sequence = ff[1]
        kmers_in_dna = kmerize(sequence,kmer_size)
        kmerdict[header] = dict((x, 1) for x in kmers_in_dna)
        # kmerdict[header] = kmers_in_dna
    print("\nAll " + str(kmer_size)+ "-mers generated!\n")
    return kmerdict

#############################################################

def main():
    parser = argparse.ArgumentParser(description="Finds representative kmers for a set of sequences")
    parser.add_argument("-f","--fasta_file", help="A fasta file of sequences",required=True)
    parser.add_argument("-o","--out_file", help="Output file with representative kmers",required=True)
    parser.add_argument("-k","--max_kmer_len", help="The maximum kmer length you'd like to try",required=False,default='31')
    args = parser.parse_args()

    fw = open(args.out_file, 'w')

    for k in range(1,int(args.max_kmer_len)+1,2):
        kmerized_dir = {}
        kmerized_dir = kmerize_directory(args.fasta_file,k)

        kmers = set()
        for seq in kmerized_dir:
            for kmer in kmerized_dir[seq]:
                kmers.add(kmer)

        fw.write(str(k) + ',' + str(len(kmers)) + '\n')

    fw.close()


if __name__ == '__main__':
    main()