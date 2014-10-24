#!/bin/bash

rm -r 2e
rm -r 2e_cmds

mkdir 2e
mkdir 2e/sgd
mkdir 2e/implicit
mkdir 2e_cmds

alphas=(60 70 80 90 100 110 120 130 140 150)
ns=(1000 2000 3000 4000 5000 6000 7000 8000 9000 10000)

for alpha in "${alphas[@]}"
do
  for n in "${ns[@]}"
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