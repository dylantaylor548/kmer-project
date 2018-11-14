import random
import os
import sys

def sample_number(number,file_path,destination_folder):
    opened_file = open(file_path)
    file_dict = {}
    sample_dict = {}
    
    for line in opened_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            sequence_name = line
            continue
        sequence = line
        file_dict[sequence_name] = sequence

    sequences = [key for key in file_dict.keys()]
    random.shuffle(sequences)

    sequences = sequences[1:number+1]

    for sequence in sequences:
        sample_dict[sequence] = file_dict[sequence]

    destination_path = destination_folder + '/' + str(number) + '_' + 'sample.fasta'
    f = open(destination_path, 'w')
    
    for sequence in sample_dict:
        f.write(sequence+'\n'+sample_dict[sequence]+'\n')

    f.close()



chosen_file = input('Choose a file to sample: ')
samp_num = int(input('Choose a number of sequences to sample from the file: '))
destination = input('Choose a directory to place the sample file in: ')

sample_number(samp_num,chosen_file,destination)
