#!/bin/bash
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 6
#SBATCH -J HMseu
#SBATCH -p scavenge
#SBATCH -t 3-
#SBATCH --mem-per-cpu=30G
#SBATCH -o HMseu_out.txt
#SBATCH -e HMseu_err.txt


srun -N 1 -n 1 Rscript inte.HMdb.v2.R &
srun -N 1 -n 1 Rscript inte.HMdb.v3.R &
wait



