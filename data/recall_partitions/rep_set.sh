#!/bin/bash

#SBATCH --array=1-10
#SBATCH --partition default
#SBATCH --time=0-10:00:00
#SBATCH --mem-per-cpu=10GB

#load modules you want
module load samtools

file=`head -n 1 ${SLURM_ARRAY_TASK_ID} list_files.txt | tail -n 1`
#file=`head -n 1 ${1} | tail -n 1`
rep_basename=`basename ${file}`
filename="${rep_basename%.*}"
python /fs/cbcb-scratch/dtaylo12/kmer-project/utils/kmer_project_master.py -f ${file} -o ${filename}.rep_set
