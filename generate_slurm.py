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

def generate_slurm2(dierectory: str, nodes: list, file_format: str, datas: list, multiple_run: int):
    for n in nodes:
        with open("./" + dierectory + "/" + file_format.format(n) + ".slurm", "w") as f:
            print("#!/bin/bash", file=f)
            print("#SBATCH --time=0:45:00", file=f)
            print("# nodes=min-max", file=f)
            print("#SBATCH --nodes={}".format(n), file=f)
            print("#SBATCH --mem=32G", file=f)
            print("#SBATCH --partition=snowy", file=f)
            print("#SBATCH --ntasks-per-node=1", file=f)
            print("#SBATCH --cpus-per-task=2", file=f)
            print("#SBATCH --job-name=2-cpt-1-npn", file=f)
            print("#SBATCH --output=script.out", file=f)
            print("# You need to load a compiler before openmpi.", file=f)
            print("", file=f)
            print("module load gcc/8.3.0", file=f)
            print("module load openmpi/3.1.4 ", file=f)
            print("", file=f)
            print("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK", file=f)
            print("", file=f)
            print("mpicxx -O3 -fopenmp mc_nqueen_parallel.cpp -o nqueen", file=f)
            for i in range(1, multiple_run+1):
                for nqueen, k in datas:
                    print("mpirun nqueen < ../nqueen_{}_{}.data > nqueen_{}_{}_{}_{}.out".format(nqueen, k, n, nqueen, k, i), file=f)

def generate_slurm3(dierectory: str, nodes: list, file_format: str, datas: list):
    for n in nodes:
        with open("./" + dierectory + "/" + file_format.format(n) + ".slurm", "w") as f:
            print("#!/bin/bash", file=f)
            print("#SBATCH --time=0:15:00", file=f)
            print("# nodes=min-max", file=f)
            print("#SBATCH --nodes={}".format(n), file=f)
            print("#SBATCH --mem=32G", file=f)
            print("#SBATCH --partition=snowy", file=f)
            print("#SBATCH --ntasks-per-node=1", file=f)
            print("#SBATCH --cpus-per-task=2", file=f)
            print("#SBATCH --job-name=2-cpt-1-npn", file=f)
            print("#SBATCH --output=script.out", file=f)
            print("# You need to load a compiler before openmpi.", file=f)
            print("", file=f)
            print("module load gcc/8.3.0", file=f)
            print("module load openmpi/3.1.4 ", file=f)
            print("", file=f)
            print("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK", file=f)
            print("", file=f)
            print("mpicxx -O3 -fopenmp mc_nqueen_parallel.cpp -o nqueen", file=f)
            
            for nqueen, k in datas:
                print("mpirun nqueen < ../nqueen_{}_{}.data > nqueen_{}_{}_{}.out".format(nqueen, k, n, nqueen, k), file=f)

if __name__ == "__main__":
    generate_slurm("n2_sequential", [1], "{}-1", [10, 100, 1000, 10000], minutes=30)
    generate_slurm("nlogn_sequential", [1], "{}-1", [10, 100, 1000, 10000])
    generate_slurm("n2_openmpi", [i for i in range(2, 13)], "{}-1", [10, 100, 1000, 10000])
    # datas = [(8, 90), (50, 90), (100, 90)]
    # generate_slurm2("mc_nqueen_parallel", [i for i in range(2, 13)], "{}", datas, 5)
    # generate_slurm3("mc_nqueen_parallel", [i for i in range(2, 13)], "{}", datas)
