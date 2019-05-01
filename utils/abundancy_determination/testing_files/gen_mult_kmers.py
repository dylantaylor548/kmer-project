import random

fw = open('multiple_pres_kmer.fasta','w')
seq_len = 1002
kmer = 'ATCGCGCTATAGAGCTCGAGA'

alphabet = ['A','C','G','T']


fw.write('>NC_00001\nATG')
for i in range(0,400):
	fw.write(random.choice(alphabet))
fw.write(kmer)
for i in range(424,seq_len):
	fw.write(random.choice(alphabet))
fw.write('\n>NC_00002\nATG')
for i in range(0,300):
	fw.write(random.choice(alphabet))
fw.write(kmer)
for i in range(324,800):
	fw.write(random.choice(alphabet))
fw.write(kmer)
for i in range(821,seq_len):
	fw.write(random.choice(alphabet))

fw.close()


