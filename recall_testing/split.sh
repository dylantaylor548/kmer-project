#!/bin/sh

module load samtools 
samtools faidx $1
num_seqs=`wc -l ${1}.fai | cut -f 1 -d ' '`
ratio_train=$(bc <<< "$num_seqs * $2")
ratio_test=$(bc <<< "$num_seqs - $ratio_train")
cut -f 1 ${1}.fai | sort -R >  ${3}.names
head -n ${ratio_train%.*} ${3}.names > ${3}.train_names
tail -n ${ratio_test%.*} ${3}.names > ${3}.test_names
xargs samtools faidx ${1} < ${3}.train_names > ${3}_train.fasta
xargs samtools faidx ${1} < ${3}.test_names > ${3}_test.fasta 
rm ${3}.names ${3}.train_names ${3}.test_names
