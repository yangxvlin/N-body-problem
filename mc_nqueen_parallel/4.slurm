#!/bin/bash
#SBATCH --time=0:10:00
# nodes=min-max
#SBATCH --nodes=4
#SBATCH --mem=32G
#SBATCH --partition=snowy
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=2-cpt-1-npn
#SBATCH --output=script.out
# You need to load a compiler before openmpi.

module load gcc/8.3.0
module load openmpi/3.1.4 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpicxx -O3 -fopenmp mc_nqueen_parallel.cpp -o nqueen
mpirun nqueen < test1.data > 4.out
