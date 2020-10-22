#!/bin/bash
#SBATCH --time=0:15:00
# nodes=min-max
#SBATCH --nodes=12
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

mpicxx -std=c++14 -O3 -o nlogn_openmpi nlogn_openmpi.cpp
mpirun nlogn_openmpi < ../body_10.data > 12-1-10.out
mpirun nlogn_openmpi < ../body_100.data > 12-1-100.out
mpirun nlogn_openmpi < ../body_500.data > 12-1-500.out
mpirun nlogn_openmpi < ../body_1000.data > 12-1-1000.out
mpirun nlogn_openmpi < ../body_2000.data > 12-1-2000.out
