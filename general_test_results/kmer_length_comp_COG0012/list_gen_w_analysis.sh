#!/bin/sh

#SBATCH --array=7,9,11,13,15,17,19,21,23,25,27,29,31
#SBATCH -J list_gen 												# Job name
#SBATCH -o list_gen.o%j												# Name of output file
#SBATCH -e list_gen.e%j												# Name of error file
#SBATCH --mail-user=dylantaylor548@gmail.com						# Email for job info
#SBATCH --mail-type=all												# Get email for begin, end, and fail
#SBATCH --time=100:00:00											# how long you think your job will take to complete; format=hh:mm:ss
#SBATCH --qos=workstation											# set QOS, this will determine what resources can be requested
#SBATCH --nodes=1													# number of nodes to allocate for your job
#SBATCH --ntasks=1													# request 4 cpu cores be reserved for your node total
#SBATCH --mem=48gb
j=`echo $SLURM_ARRAY_TASK_ID`

python list_gen_w_analysis.py -f COG0012_comp_true.fasta -o COG0012_comp_true_k${j}.list -c 1 -k ${j}