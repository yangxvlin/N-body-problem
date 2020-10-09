#!/bin/bash
#SBATCH --time=0:10:00
# nodes=min-max
#SBATCH --nodes=1
#SBATCH --mem=32G 
# You need to load a compiler before openmpi.

module load gcc/8.3.0
module load openmpi/3.1.4 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

g++ -std=c++14 -fopenmp -O3 -o n2_openmp n2_openmp.cpp
./n2_openmp 10 6.67e-11 1.0 < ../body_10000.data