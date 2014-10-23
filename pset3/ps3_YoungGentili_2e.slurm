#!/bin/bash

mkdir 2e
mkdir 2e_cmds
ns=(100 500 1000 2500 1e4 1e5 1e6)
alphas=(60 70 80 90 100 110 120 130 140 150)

for n in "${ns[@]}"
do
	for alpha in "${alphas[@]}"
	do
		echo "#!/bin/sh" > 2e_cmds/${n}_${alpha}.sh
		echo '#SBATCH -J ps3_2e # name for job array' >> 2e_cmds/${n}_${alpha}.sh
		echo '#SBATCH -o all.out #Standard output' >> 2e_cmds/${n}_${alpha}.sh
		echo '#SBATCH -e all.err #Standard error' >> 2e_cmds/${n}_${alpha}.sh
		echo '#SBATCH -p serial_requeue #Partition' >> 2e_cmds/${n}_${alpha}.sh
		echo '#SBATCH --mem-per-cpu 4096 #Memory request' >> 2e_cmds/${n}_${alpha}.sh
		echo '#SBATCH -n 1 #Number of cores' >> 2e_cmds/${n}_${alpha}.sh
		echo '#SBATCH -N 1 #All cores on one machine' >> 2e_cmds/${n}_${alpha}.sh
		echo "Rscript ps3_YoungGentili_2e.R ${n} ${alpha}" >> 2e_cmds/${n}_${alpha}.sh
		cat 2e_cmds/${n}_${alpha}.sh | sbatch
	done
done