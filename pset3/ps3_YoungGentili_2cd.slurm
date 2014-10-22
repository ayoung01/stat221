#!/bin/bash

mkdir 2cd
mkdir 2cd_cmds
ns=(100 1000 5000 1e4 1e5 1e6)
alphas=(1e-4 1e-3 1e-2 0.1 0.5 1 10 20 50 100)

for n in "${ns[@]}"
do
	for alpha in "${alphas[@]}"
	do
		echo "#!/bin/sh" > cmds/${n}_${p}_${rho}.sh
		echo '#SBATCH -J ps3_3 # name for job array' >> cmds/${n}_${p}_${rho}.sh
		echo '#SBATCH -o all.out #Standard output' >> cmds/${n}_${p}_${rho}.sh
		echo '#SBATCH -e all.err #Standard error' >> cmds/${n}_${p}_${rho}.sh
		echo '#SBATCH -p serial_requeue #Partition' >> cmds/${n}_${p}_${rho}.sh
		echo '#SBATCH --mem-per-cpu 4096 #Memory request' >> cmds/${n}_${p}_${rho}.sh
		echo '#SBATCH -n 1 #Number of cores' >> cmds/${n}_${p}_${rho}.sh
		echo '#SBATCH -N 1 #All cores on one machine' >> cmds/${n}_${p}_${rho}.sh
		echo "Rscript ps3_YoungGentili_3.R ${n} ${p} ${rho} 3" >> cmds/${n}_${p}_${rho}.sh
		cat cmds/${n}_${p}_${rho}.sh | sbatch
	done
done