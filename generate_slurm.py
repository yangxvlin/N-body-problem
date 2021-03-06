# Author: Xulin Yang, 904904
# generate slurm scripts to be run on Spartan

def generate_slurm(dierectory: str, nodes: list, file_format: str, bodies: list, minutes=15):
    for n in nodes:
        with open("./" + dierectory + "/" + file_format.format(n) + ".slurm", "w") as f:
            print("#!/bin/bash", file=f)
            print("#SBATCH --time=0:{}:00".format(minutes), file=f)
            print("# nodes=min-max", file=f)
            print("#SBATCH --nodes={}".format(n), file=f)
            print("#SBATCH --mem=32G", file=f)
            print("#SBATCH --partition=snowy", file=f)
            print("#SBATCH --ntasks-per-node=1", file=f)
            print("#SBATCH --cpus-per-task=1", file=f)
            print("#SBATCH --job-name=1-cpt-1-npn", file=f)
            print("#SBATCH --output=script.out", file=f)
            print("# You need to load a compiler before openmpi.", file=f)
            print("", file=f)
            print("module load gcc/8.3.0", file=f)
            print("module load openmpi/3.1.4 ", file=f)
            print("", file=f)
            print("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK", file=f)
            print("", file=f)
            print("mpicxx -std=c++14 -O3 -o {} {}.cpp".format(dierectory, dierectory), file=f)
            for i in bodies:
                print("mpirun {} < ../body_{}.data > {}-1-{}.out".format(dierectory, i, n, i), file=f)

def generate_slurm4(dierectory: str, nodes: list, cores: list, file_format: str, bodies: list, minutes=15):
    for n in nodes:
        for core in cores:
            with open("./" + dierectory + "/" + file_format.format(n, core) + ".slurm", "w") as f:
                print("#!/bin/bash", file=f)
                print("#SBATCH --time=0:{}:00".format(minutes), file=f)
                print("# nodes=min-max", file=f)
                print("#SBATCH --nodes={}".format(n), file=f)
                print("#SBATCH --mem=32G", file=f)
                print("#SBATCH --partition=snowy", file=f)
                print("#SBATCH --ntasks-per-node=1", file=f)
                print("#SBATCH --cpus-per-task={}".format(core), file=f)
                print("#SBATCH --job-name=1-cpt-1-npn", file=f)
                print("#SBATCH --output=script.out", file=f)
                print("# You need to load a compiler before openmpi.", file=f)
                print("", file=f)
                print("module load gcc/8.3.0", file=f)
                print("module load openmpi/3.1.4 ", file=f)
                print("", file=f)
                print("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK", file=f)
                print("", file=f)
                print("mpicxx -std=c++14 -fopenmp -O3 -o {} {}.cpp".format(dierectory, dierectory), file=f)
                for i in bodies:
                    print("mpirun {} < ../body_{}.data > {}-{}-{}.out".format(dierectory, i, n, core, i), file=f)

def generate_slurm5(dierectory: str, node: int, cores: list, file_format: str, bodies: list, minutes=15):
    print(node)
    for core in cores:
        with open("./" + dierectory + "/" + file_format.format(node, core) + ".slurm", "w") as f:
            print("#!/bin/bash", file=f)
            print("#SBATCH --time=0:{}:00".format(minutes), file=f)
            print("# nodes=min-max", file=f)
            print("#SBATCH --nodes={}".format(node), file=f)
            print("#SBATCH --mem=32G", file=f)
            print("#SBATCH --partition=snowy", file=f)
            print("#SBATCH --ntasks-per-node=1", file=f)
            print("#SBATCH --cpus-per-task={}".format(core), file=f)
            print("#SBATCH --job-name=1-cpt-1-npn", file=f)
            print("#SBATCH --output=script.out", file=f)
            print("# You need to load a compiler before openmpi.", file=f)
            print("", file=f)
            print("module load gcc/8.3.0", file=f)
            print("module load openmpi/3.1.4 ", file=f)
            print("", file=f)
            print("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK", file=f)
            print("", file=f)
            print("mpicxx -std=c++14 -fopenmp -O3 -o {} {}.cpp".format(dierectory, dierectory), file=f)
            for i in bodies:
                print("mpirun {} < ../body_{}.data > {}-{}-{}.out".format(dierectory, i, node, core, i), file=f)

if __name__ == "__main__":
    generate_slurm("n2_sequential", [1], "{}-1", [10, 100, 500, 1000, 2000], minutes=15)
    generate_slurm("n2_openmpi", [i for i in range(2, 13)], "{}-1", [10, 100, 500, 1000, 2000], minutes=15)
    generate_slurm("n2_openmpi_profile", [i for i in range(2, 13)], "{}-1", [2000], minutes=15)
    generate_slurm("nlogn_sequential", [1], "{}-1", [10, 100, 500, 1000, 2000, 5000], minutes=15)
    generate_slurm("nlogn_openmpi", [i for i in range(2, 13)], "{}-1", [10, 100, 500, 1000, 2000, 5000], minutes=15)
    generate_slurm("nlogn_openmpi_profile", [i for i in range(2, 13)], "{}-1", [2000], minutes=15)

    generate_slurm4("nlogn_hybrid", [12], [i for i in range(2, 17)], "{}-{}", [500, 1000, 2000], minutes=15)
    generate_slurm4("n2_hybrid", [12], [i for i in range(2, 17)], "{}-{}", [500, 1000, 2000], minutes=15)
    
    generate_slurm5("n2_hybrid_profile", 12, [i for i in range(2, 17)], "{}-{}", [2000], minutes=15)
    generate_slurm5("nlogn_hybrid_profile", 12, [i for i in range(2, 17)], "{}-{}", [2000], minutes=15)
