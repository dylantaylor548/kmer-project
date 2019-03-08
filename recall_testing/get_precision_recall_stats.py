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


def store_seqs_w_id(file_path,identifier):
    interest_dict = {}
    noise_dict = {}

    fiter = fasta_iter(file_path)
    for ff in fiter:
        header = ff[0]
        sequence = ff[1]
        if identifier in header:
            interest_dict[header] = sequence
        else:
            noise_dict[header] = sequence
    return interest_dict, noise_dict


def store_kmerlist(oligo_csv_path):
    oligo_list = []
    with open(oligo_csv_path) as opened_file:
        for line in opened_file:
            line = line.strip()
            oligo_list.append(line)
    return oligo_list


def get_complement(DNA_seq):
    DNA = DNA_seq.split()
    for i in range(0,len(DNA)):
        base = DNA[i]
        if base == 'A':
            DNA[i] = 't'
        if base == 'C':
            DNA[i] = 'g'
        if base == 'G':
            DNA[i] = 'c'
        if base == 'T':
            DNA[i] = 'a'
    DNA_comp = ''.join(DNA)
    DNA_comp = DNA_comp.upper()
    return DNA_comp


##################################################################################

def main():
    parser = argparse.ArgumentParser(description="Calculates recall and precision statistics for an oligo list against a metagenomic sample")
    parser.add_argument("-f","--sample_fasta",help='A .fasta file containing the sequences in your metagenomic "sample"',required=True)
    parser.add_argument("-id","--identifier",help='The identifier for the correct sequence within the fasta file (i.e. "COG0088")',required=True)
    parser.add_argument("-l","--oligo_list",help='A .csv file containing the representative list from your training data set\n(i.e. only the sequences you are interested in capturing)',required=True)
    parser.add_argument("-o","--output_file",help='An output .csv file to store data about the recall and precision statistics',required=False)
    args = parser.parse_args()

    oligo_list = store_kmerlist(args.oligo_list)

    interest_dict, noise_dict = store_seqs_w_id(args.sample_fasta,args.identifier)

    true_pos = 0
    true_neg = 0
    false_pos = 0
    false_neg = 0

    for seq_name in interest_dict:
        covered = False
        for oligo in oligo_list:
            complement = get_complement(oligo)
            if (oligo in interest_dict[seq_name]) or (complement in interest_dict[seq_name]):
                covered = True
                true_pos += 1
                break
        if not covered:
            false_neg += 1

    for seq_name in noise_dict:
        covered = False
        for oligo in oligo_list:
            complement = get_complement(oligo)
            if (oligo in noise_dict[seq_name]) or (complement in noise_dict[seq_name]):
                covered = True
                false_pos += 1
                break
        if not covered:
            true_neg += 1

    recall = None
    if true_pos != 0:
    	recall = float(Fraction(true_pos,(true_pos + false_neg))*100)
    else:
    	recall = 0
    precision = float(Fraction(true_pos,(true_pos + false_pos))*100)

    print("Recall: " + str(recall) + "%")
    print("Precision: " + str(precision) + "%")
    print("Of the " + str(len(noise_dict)) + " sequences that were not " + args.identifier + ", we identified " + str(true_neg) + " as true negatives")
##################################################################################

main()