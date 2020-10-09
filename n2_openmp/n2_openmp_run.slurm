#!/bin/bash
#SBATCH --output=run.out
#SBATCH --time=0:01:00

# for i in {2..16}
# do
#    sbatch --output=1node-1-cpt-${i}-npn-physical.out --job-name=1-cpt-${i}-npn --cpus-per-task=1 --partition=physical --ntasks-per-node=${i} "n2_openmp_template.slurm"
# done

sbatch --output=1node-1-cpt-16-npn-physical.out --job-name=1-cpt-16-npn --cpus-per-task=1 --partition=physical --ntasks-per-node=16 "n2_openmp_template.slurm"
