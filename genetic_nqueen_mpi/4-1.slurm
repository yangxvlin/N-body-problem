#!/bin/bash
#SBATCH --time=0:10:00
# nodes=min-max
#SBATCH --nodes=4
#SBATCH --mem=32G
#SBATCH --partition=snowy
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=1-cpt-1-npn
#SBATCH --output=script.out
# You need to load a compiler before openmpi.

module load gcc/8.3.0
module load openmpi/3.1.4 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpicxx genetic_nqueen_mpi.cpp -o nqueen -O3
mpirun nqueen > 4-1.out
