#!/bin/bash
#SBATCH --output=run.out
#SBATCH --time=0:01:00

# for i in {2..16}
# do
#    sbatch --output=1node-1-cpt-${i}-npn-snowy.out --job-name=1-cpt-${i}-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=${i} "n2_openmpi_template.slurm"
# done

sbatch --output=4node-1-cpt-1-npn-snowy.out --error=4node-1-cpt-1-npn-snowy.out --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=4 "n2_openmpi_template.slurm"
