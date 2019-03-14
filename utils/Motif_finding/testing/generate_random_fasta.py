import random

def generate_fasta(out_path,N_seq,seq_length):
	bases = ['A','C','G','T']
	f = open(out_path,'w')
	for x in range(1,N_seq+1):
		seq = ''
		for bp in range(0,seq_length):
			seq += random.choice(bases)
		f.write('>seq_' + str(x) + '\n' + seq + '\n')
	f.close()

out_path = 'C:/Users/Dylan/Desktop/random_seqs.fasta'
N_seq = 20
seq_length = 100

generate_fasta(out_path,N_seq,seq_length)
