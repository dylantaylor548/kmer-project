import os

jobs = 40

file_names = []
for i in range(0,jobs):
	file_name = 'job00'[:-len(str(i+1))] + str(i+1) + '.sh'
	file_names.append(file_name)
	job = open(file_name,'w')
	job.write('#!/bin/sh\npython fastq_generator.py -f TIGR02012_83.fasta -c test_chip.csv -s test_spikes.fasta -o simulated_reads/' + file_name[:-3] + '_reads.fastq -l 200')
	job.close()

for file_name in os.listdir(os.getcwd()):
	if file_name in file_names:
		os.system('sbatch --mem=12gb --qos=throughput --time=18:00:00 --nodes=1 --ntasks=4 ' + file_name)
		os.system('rm ' + file_name)