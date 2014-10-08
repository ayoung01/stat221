#!/bin/bash
#SBATCH -J task3 # name for job array
#SBATCH -o all.out #Standard output
#SBATCH -e all.err #Standard error
#SBATCH -p serial_requeue #Partition
#SBATCH --mem-per-cpu 4096 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine

# 1st argument = M (number of Ys to generate per node)
# 2nd argument = job id
# sbatch --array=1-12 YoungGentili_ps2_task3.slurm
mkdir out
Rscript YoungGentili_ps2_task5.R 50 $SLURM_ARRAY_TASK_ID