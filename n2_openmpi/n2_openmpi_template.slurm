#!/bin/bash
#SBATCH --time=0:10:00
# nodes=min-max
#SBATCH --nodes=1
#SBATCH --mem=32G 
# You need to load a compiler before openmpi.

module load gcc/8.3.0
module load openmpi/3.1.4 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpicxx -std=c++14 -O3 -o n2_openmpi n2_openmpi.cpp
mpirun -np $SLURM_NODES n2_openmpi < ../body_10.data