#!/bin/bash
#SBATCH -J task3 # name for job array
#SBATCH -o all.out #Standard output
#SBATCH -e all.err #Standard error
#SBATCH -p general #Partition
#SBATCH -t 01:00:00 #Running time of 1 hour
#SBATCH --mem-per-cpu 3000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine

# 1st argument = job id
Rscript YoungGentili_ps_task3.R $SLURM_ARRAY_TASK_ID