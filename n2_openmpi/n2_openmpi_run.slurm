#!/bin/bash
#SBATCH --output=run.out
#SBATCH --time=0:01:00

# for i in {2..16}
# do
#    sbatch --output=${i}node-1-cpt-1-npn-snowy.out --error=${i}node-1-cpt-1-npn-snowy.out --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=${i} "n2_openmpi_template.slurm"
# done

# remember to switch to 10 body test file
# sbatch --error=2node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=2 "n2_openmpi_template.slurm"
# sbatch --error=3node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=3 "n2_openmpi_template.slurm"
# sbatch --error=4node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=4 "n2_openmpi_template.slurm"
# sbatch --error=5node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=5 "n2_openmpi_template.slurm"
# sbatch --error=6node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=6 "n2_openmpi_template.slurm"
# sbatch --error=7node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=7 "n2_openmpi_template.slurm"
# sbatch --error=8node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=8 "n2_openmpi_template.slurm"
# sbatch --error=9node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=9 "n2_openmpi_template.slurm"
# sbatch --error=10node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=10 "n2_openmpi_template.slurm"
# sbatch --error=11node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=11 "n2_openmpi_template.slurm"
# sbatch --error=12node-1-cpt-1-npn-snowy.err --job-name=1-cpt-1-npn --cpus-per-task=1 --partition=snowy --ntasks-per-node=1 --nodes=12 "n2_openmpi_template.slurm"