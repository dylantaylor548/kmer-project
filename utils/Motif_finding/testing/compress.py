from itertools import groupby


#Reads a fasta file
def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
        yield header, seq


# For a directory containing files each of one DNA sequence, creates a dictionary
# with keywords: file_name and values: list of all kmers in that file
# (This has been tested and works)
def dict_directory(file_path):
    seqdict = {}

    fiter = fasta_iter(file_path)
    for ff in fiter:
        header = ff[0]
        sequence = ff[1]
        if sequence not in list(seqdict.values()):
            seqdict[header] = sequence

    return seqdict

def write_seqdict(seqdict,file_path):
    f = open(file_path,'w')
    for sequence in seqdict:
        f.write('>' + sequence + '\n' + seqdict[sequence] + '\n')
    f.close()

file2compress = 'C:/Users/Dylan/Desktop/Pop_Lab/kmer_project/data/COG0088/COG0088.fna'
outpath = 'C:/Users/Dylan/Desktop/compressed_COG0088.fasta'

seqdict = dict_directory(file2compress)
write_seqdict(seqdict,outpath)