#!/bin/bash

mkdir out
mkdir cmds
ns=(100 1000 5000 1e4 1e5 1e6)
ps=(100 1000 5000 20000 50000)
rhos=(0 0.1 0.2 0.5 0.9 0.95)

for n in "${ns[@]}"
do
	for p in "${ps[@]}"
	do
		for rho in "${rhos[@]}"
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
done