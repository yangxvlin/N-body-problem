#!/bin/bash
#SBATCH --time=0:10:00
# nodes=min-max
#SBATCH --nodes=1
#SBATCH --mem=32G 
#SBATCH --job-name=1-cpt-1-npn 
#SBATCH --cpus-per-task=1 
#SBATCH --partition=snowy 
#SBATCH --ntasks-per-node=1
# You need to load a compiler before openmpi.

module load gcc/8.3.0
module load openmpi/3.1.4 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

g++ -std=c++14 -O3 -o nlogn_sequential nlogn_sequential.cpp
./nlogn_sequential < ../body_10.data > 1-{$SLURM_CPUS_PER_TASK}-10.out
./nlogn_sequential < ../body_1000.data > 1-{$SLURM_CPUS_PER_TASK}-1000.out
./nlogn_sequential < ../body_10000.data > 1-{$SLURM_CPUS_PER_TASK}-10000.out