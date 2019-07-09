#!/bin/sh

#SBATCH --array=1-40
#SBATCH -J read_gen # Job name
#SBATCH -o read_gen.o%j # Name of output file
#SBATCH -e read_gen.e%j # Name of error file
#SBATCH --mail-user=dylantaylor548@gmail							# Email for job info
#SBATCH --mail-type=all												# Get email for begin, end, and fail
#SBATCH --time=100:00:00											# how long you think your job will take to complete; format=hh:mm:ss
#SBATCH --qos=workstation											# set QOS, this will determine what resources can be requested
#SBATCH --nodes=1													# number of nodes to allocate for your job
#SBATCH --ntasks=4													# request 4 cpu cores be reserved for your node total
#SBATCH --mem=48gb
j=`echo $SLURM_ARRAY_TASK_ID`

python fastq_generator.py -f TIGR02012_83.fasta -c test_chip.csv -s test_spikes.fasta -o simulated_reads/read_gen${j}_reads.fastq -l 200