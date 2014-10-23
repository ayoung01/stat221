#!/bin/sh
#SBATCH -J ps3_2cd # name for job array
#SBATCH -o all.out #Standard output
#SBATCH -e all.err #Standard error
#SBATCH -p serial_requeue #Partition
#SBATCH --mem-per-cpu 4096 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
Rscript ps3_YoungGentili_2cd.R 2500 10
