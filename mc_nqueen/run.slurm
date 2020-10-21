#!/bin/bash
#SBATCH --time=0:10:00
# nodes=min-max
#SBATCH --nodes=1
#SBATCH --mem=32G 
#SBATCH --job-name=1-cpt-1-npn 
#SBATCH --cpus-per-task=1 
#SBATCH --partition=snowy 
#SBATCH --ntasks-per-node=1
#SBATCH --output=script.out
# You need to load a compiler before openmpi.

module load gcc/8.3.0
module load openmpi/3.1.4 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

g++ -O3 mc_nqueen.cpp -o nqueen
./nqueen < ../nqueen_8_90.data > nqueen_8_90.out
./nqueen < ../nqueen_50_90.data > nqueen_50_90.out
./nqueen < ../nqueen_100_90.data > nqueen_100_90.out