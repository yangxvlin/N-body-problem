#!/bin/bash
#SBATCH --time=0:15:00
# nodes=min-max
#SBATCH --nodes=12
#SBATCH --mem=32G
#SBATCH --partition=snowy
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=13
#SBATCH --job-name=1-cpt-1-npn
#SBATCH --output=script.out
# You need to load a compiler before openmpi.

module load gcc/8.3.0
module load openmpi/3.1.4 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpicxx -std=c++14 -fopenmp -O3 -o nlogn_hybrid nlogn_hybrid.cpp
mpirun nlogn_hybrid < ../body_500.data > 12-13-500.out
mpirun nlogn_hybrid < ../body_1000.data > 12-13-1000.out
mpirun nlogn_hybrid < ../body_2000.data > 12-13-2000.out
