#!/bin/sh

#SBATCH --array=27,29,31
#SBATCH -J k_dep													# Job name
#SBATCH -o k_dep.o%j												# Name of output file
#SBATCH -e k_dep.e%j												# Name of error file
#SBATCH --mail-user=dylantaylor548@gmail							# Email for job info
#SBATCH --mail-type=all												# Get email for begin, end, and fail
#SBATCH --time=100:00:00											# how long you think your job will take to complete; format=hh:mm:ss
#SBATCH --qos=workstation											# set QOS, this will determine what resources can be requested
#SBATCH --nodes=1													# number of nodes to allocate for your job
#SBATCH --ntasks=4													# request 4 cpu cores be reserved for your node total
#SBATCH --mem=48gb
j=`echo $SLURM_ARRAY_TASK_ID`

python single_COG_fullsearch.py -l ../general_test_results/kmer_length_comp_COG0012/COG0012_comp_true_k${j}.list -d genomes_annotations/ -o ../general_test_results/kmer_length_comp_COG0012/COG0012_comp_true_k${j}_wholegenome.report