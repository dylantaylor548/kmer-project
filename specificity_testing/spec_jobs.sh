#!/bin/sh

#SBATCH --array=1-40
#SBATCH -J spec_jobs # Job name
#SBATCH -o spec_jobs.o%j # Name of output file
#SBATCH -e spec_jobs.e%j # Name of error file
#SBATCH --mail-user=nidhiishahdb@gmail.com # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=100:00:00                                         # how long you think your job will take to complete; format=hh:mm:ss
#SBATCH --qos=workstation                                             # set QOS, this will determine what resources can be requested
#SBATCH --nodes=1                                               # number of nodes to allocate for your job
#SBATCH --ntasks=4                                            # request 4 cpu cores be reserved for your node total
#SBATCH --mem=48gb

module load samtools
j=`cat $1 | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1`
# j=`cat $1 | head -n 1 | tail -n 1`
# echo $j
train_file=`echo "$j" | cut -f 1 -d $'\t'`
test_file=`echo "$j" | cut -f 2 -d $'\t'`
train_file_name=`basename $train_file | cut -f 1 -d '.'`
echo $train_file_name
echo "Running partitioning of train and test file step"
python partition.py -f ${train_file} -r 0.8 -t ${test_file} -sr 1.0 -o ${train_file_name}
echo "Generating representative kmers list"
python /fs/cbcb-scratch/dtaylo12/kmer_project/src/kmer_project_master.py -f ${train_file_name}_train_seqs.fasta -o ${train_file_name}.list -c 1 -k 21
echo "Getting precision recall stats"
python /fs/cbcb-scratch/dtaylo12/kmer_project/recall_testing/get_precision_recall_stats.py -f ${train_file_name}_sample_test.fasta -id ${train_file_name} -l ${train_file_name}.list -o ${train_file_name}.report > ${train_file_name}.out

