import os


def avg(values):
	total = 0
	for item in values:
		total += item
	average = float(total)/float(len(values))
	return average


def stddev(values):
	average = avg(values)
	ssd = 0
	for item in values:
		ssd += (item - average)**2
	dev = (ssd/(len(values) - 1))**(0.5)
	return dev



cwd = os.getcwd()
for file_name in os.listdir(cwd):
	if file_name.endswith('.fastq'):
		prefix = file_name[:-6]
		cmd = 'python ../../abundancy_determination/get_relative_abundance.py -q ' + file_name + ' -f ../TIGR02012_83.fasta -s ../test_spikes.fasta -c ../test_chip.csv -o ' + prefix + '_calc_abund.csv'
		os.system(cmd)

seq_values = {}

for file_name in os.listdir(os.getcwd()):
	if file_name.endswith('_calc_abund.csv'):
		fo = open(file_name)
		for line in fo:
			line =line.strip().split(',')
			seq_id = line[0]
			seq_values.setdefault(seq_id,[])
			seq_values[seq_id].append(float(line[1]))
		os.system('rm ' + file_name)

fw = open('tmp_abundacy_report.csv','w')
fw.write('seq_id,average,stddev\n')
for seq_id in seq_values:
	fw.write(seq_id + ',' + str(avg(seq_values[seq_id])) + ',' + str(stddev(seq_values[seq_id])) + '\n')
fw.close()

