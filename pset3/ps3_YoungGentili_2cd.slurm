#!/bin/bash

mkdir 2cd
mkdir 2cd_cmds

alphas=(0.01 0.1 0.5 1 2 10 20 50 100 200)
rep.nos=(1 2 3 4 5)

for alpha in "${alpha[@]}"
do
	for rep.no in "${rep.nos[@]}"
	do
		echo "#!/bin/sh" > 2cd_cmds/${alpha}_${rep.no}.sh
		echo '#SBATCH -J ps3_2cd # name for job array' >> 2cd_cmds/${alpha}_${rep.no}.sh
		echo '#SBATCH -o all.out #Standard output' >> 2cd_cmds/${alpha}_${rep.no}.sh
		echo '#SBATCH -e all.err #Standard error' >> 2cd_cmds/${alpha}_${rep.no}.sh
		echo '#SBATCH -p serial_requeue #Partition' >> 2cd_cmds/${alpha}_${rep.no}.sh
		echo '#SBATCH --mem-per-cpu 4096 #Memory request' >> 2cd_cmds/${alpha}_${rep.no}.sh
		echo '#SBATCH -n 1 #Number of cores' >> 2cd_cmds/${alpha}_${rep.no}.sh
		echo '#SBATCH -N 1 #All cores on one machine' >> 2cd_cmds/${alpha}_${rep.no}.sh
		echo "Rscript ps3_YoungGentili_2cd.R ${alpha} ${rep.no}" >> 2cd_cmds/${alpha}_${rep.no}.sh
		cat 2cd_cmds/${alpha}_${rep.no}.sh | sbatch
	done
done