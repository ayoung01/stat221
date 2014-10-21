#!/bin/bash
#SBATCH -J ps3_3 # name for job array
#SBATCH -o all.out #Standard output
#SBATCH -e all.err #Standard error
#SBATCH -p serial_requeue #Partition
#SBATCH --mem-per-cpu 4096 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=albert.t.young@gmail.com  #Email to which notifications will be sent

mkdir out
mkdir cmds
ns=(100, 1000, 5000, 1e4, 1e5, 1e6)
ps=(100, 1000, 5000, 20000, 50000)
rhos=(0, 0.1, 0.2, 0.5, 0.9, 0.95)

for n in "${ns[@]}"
do
	for p in "${ps[@]}"
	do
		for rho in "${rhos[@]}"
		do
			echo "Rscript ps3_YoungGentili_3.R ${n}$ ${p} ${rho} 3" > cmds/${n}_${p}_${rho}.sh
			cat cmds/${n}_${p}_${rho}.sh | sbatch
		done
	done
done