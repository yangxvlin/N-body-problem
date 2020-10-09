#!/bin/bash
#SBATCH --time=0:10:00
# nodes=min-max
#SBATCH --nodes=1
#SBATCH --mem=32G 
#SBATCH --output=1node-1-cpt-1-npn-physical.out 
#SBATCH --job-name=1-cpt-1-npn 
#SBATCH --cpus-per-task=1 
#SBATCH --partition=physical 
#SBATCH --ntasks-per-node=1
# You need to load a compiler before openmpi.

module load gcc/8.3.0
module load openmpi/3.1.4 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

g++ -std=c++14 -O3 -o n2_sequential n2_sequential.cpp
./n2_sequential 10 6.67e-11 1.0 < ../body_1000.data