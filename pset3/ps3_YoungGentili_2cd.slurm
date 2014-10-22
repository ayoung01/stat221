#!/bin/bash

mkdir 2cd
mkdir 2cd_cmds
ns=(100 500 1000 2500 1e4 1e5 1e6)
alphas=(1e-4 1e-3 1e-2 0.1 0.5 1 10 20 50 100)

for n in "${ns[@]}"
do
	for alpha in "${alphas[@]}"
	do
		echo "#!/bin/sh" > 2cd_cmds/${n}_${alpha}.sh
		echo '#SBATCH -J ps3_2cd # name for job array' >> 2cd_cmds/${n}_${alpha}.sh
		echo '#SBATCH -o all.out #Standard output' >> 2cd_cmds/${n}_${alpha}.sh
		echo '#SBATCH -e all.err #Standard error' >> 2cd_cmds/${n}_${alpha}.sh
		echo '#SBATCH -p serial_requeue #Partition' >> 2cd_cmds/${n}_${alpha}.sh
		echo '#SBATCH --mem-per-cpu 4096 #Memory request' >> 2cd_cmds/${n}_${alpha}.sh
		echo '#SBATCH -n 1 #Number of cores' >> 2cd_cmds/${n}_${alpha}.sh
		echo '#SBATCH -N 1 #All cores on one machine' >> 2cd_cmds/${n}_${alpha}.sh
		echo "Rscript ps3_YoungGentili_2cd.R ${n} ${alpha}" >> 2cd_cmds/${n}_${alpha}.sh
		cat 2cd_cmds/${n}_${alpha}.sh | sbatch
	done
done