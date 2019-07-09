#!/bin/sh

#SBATCH --array=102,103,124,172,184,185,186,197,200,201,202,215,256,495,525,533,541,552
#SBATCH -J match_find 												# Job name
#SBATCH -o match_find.o%j											# Name of output file
#SBATCH -e match_find.e%j											# Name of error file
#SBATCH --mail-user=dylantaylor548@gmail							# Email for job info
#SBATCH --mail-type=all												# Get email for begin, end, and fail
#SBATCH --time=100:00:00											# how long you think your job will take to complete; format=hh:mm:ss
#SBATCH --qos=workstation											# set QOS, this will determine what resources can be requested
#SBATCH --nodes=1													# number of nodes to allocate for your job
#SBATCH --ntasks=4													# request 4 cpu cores be reserved for your node total
#SBATCH --mem=48gb
j=`echo $SLURM_ARRAY_TASK_ID`

python single_COG_fullsearch.py -l Whole_COGs/COG_lists/COG0${j}_deoutliered.list -d genomes_annotations -o Whole_COGs/COG_reports/COG0${j}_deoutliered.report