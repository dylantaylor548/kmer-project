import os


def avg(stuff):
	return sum(stuff)/len(stuff)


def stddevp(stuff):
	total = 0 
	stuffavg = avg(stuff)

	for elem in stuff:
		total += (elem - stuffavg)**2

	return total/len(stuff)


def count(stuff, value):
	total = 0

	for elem in stuff:
		if elem == value:
			total += 1

	return total


#########################################################################

outfile = open(os.getcwd()+'/COG0012_listgen_analysis.csv','w')

run_times = {}
num_kmers = {} 
avgs = {}
stddevs = {}
num_useless = {}

for file_name in os.listdir(os.getcwd()):
	if file_name.endswith('.list_analysis'):
		k = file_name[-16:-14]
		if k.startswith('k'):
			k = k[1]
		print('Analyzing k = ' + str(k))
		kmer_worths = []
		with open(file_name) as opened_file:
			for line in opened_file:
				line = line.strip()
				if line.startswith('runtime'):
					run_times[k] = float(line.split(',')[1])
				else:
					if len(line.split(',')) == 2:
						kmer_worths.append(float(line.split(',')[1]))
		num_kmers[k] = len(kmer_worths)
		avgs[k] = avg(kmer_worths)
		stddevs[k] = stddevp(kmer_worths)
		num_useless[k] = count(kmer_worths,1)

outfile.write('k =,run time (secs),list length,avg kmer worth,stddev kmer worth,num worth=1 kmers\n')

for k_id in run_times:
	outfile.write(k_id + ',' + str(run_times[k_id]) + ',' + str(num_kmers[k_id]) + ',' + str(avgs[k_id]) + ',' + str(stddevs[k_id]) + ',' + str(num_useless[k_id]) + '\n')

outfile.close()
