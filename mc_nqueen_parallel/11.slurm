#!/bin/bash
#SBATCH --time=0:50:00
# nodes=min-max
#SBATCH --nodes=11
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
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_8_90_1.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_50_90_1.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_100_90_1.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_8_90_2.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_50_90_2.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_100_90_2.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_8_90_3.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_50_90_3.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_100_90_3.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_8_90_4.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_50_90_4.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_100_90_4.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_8_90_5.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_50_90_5.out
mpirun nqueen < ../nqueen_11_90.data > nqueen_11_100_90_5.out
