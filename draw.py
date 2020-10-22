# Author: Xulin Yang, 904904
# draw analysis diagrams

from matplotlib import pyplot as plt
from math import log10

def draw_speedup(sequential_path: str, 
                 parallel_directory: str, 
                 parallel_output_template: str, 
                 parallel_outputs: list,
                 title="1",
                 xlabel="",
                 ylabel=""
                ):
    with open(sequential_path) as f:
        seq_time = int(f.readline().split(' = ')[1])

    parallel_time = []
    for i in parallel_outputs:
        with open(parallel_directory + str(i) + parallel_output_template) as f:
            parallel_time.append(int(f.readline().split(' = ')[1]))
    
    print(seq_time)
    print(parallel_time)
    speedup = [seq_time / i for i in parallel_time]
    print(speedup)

    plt.figure()
    plt.plot(parallel_outputs, speedup, 'bo-')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(title)

def draw_nqueen_speedp(sequential_directory: str, parallel_directory: str, nqueens: list, ks: list, nodes: list,
                xlabel, ylabel, title):
    plt.figure()

    sequential_time = []
    parallel_time = [[] for _ in range(0, len(nqueens))]

    for i, nqueen in enumerate(nqueens):
        for k in ks:
            with open(sequential_directory + "nqueen_" + str(nqueen) + "_" +  str(k) + ".out") as f:
                while True:
                    tmp = f.readline().split(' = ')
                    if len(tmp) > 1:
                        sequential_time.append(int(tmp[1]))
                        break
            for n in nodes:
                with open(parallel_directory + "nqueen_" + str(n) + "_" + str(nqueen) + "_" +  str(k) + ".out") as f:
                    while True:
                        tmp = f.readline().split(' = ')
                        if len(tmp) > 1:
                            parallel_time[i].append(int(tmp[1]))
                            break
    print(sequential_time)
    print(parallel_time)

    speedup = []
    for i, times in enumerate(parallel_time):
        speedup.append([])
        for t in times:
            speedup[-1].append(sequential_time[i] / t)
    print(speedup)

    colors = ["r", "g", "b", "h", "i"]

    plt.figure()
    for i in range(0, len(nqueens)):
        plt.plot(nodes, speedup[i], colors[i] + 'o-')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend([str(i) + " queens" for i in nqueens])
    plt.savefig(title)

def draw_nqueen_time(sequential_directory: str, parallel_directory: str, nqueens: list, ks: list, nodes: list,
                xlabel, ylabel, title):
    plt.figure()

    sequential_time = []
    parallel_time = [[] for _ in range(0, len(nqueens))]

    for i, nqueen in enumerate(nqueens):
        for k in ks:
            with open(sequential_directory + "nqueen_" + str(nqueen) + "_" +  str(k) + ".out") as f:
                while True:
                    tmp = f.readline().split(' = ')
                    if len(tmp) > 1:
                        sequential_time.append(int(tmp[1]))
                        break
            for n in nodes:
                with open(parallel_directory + "nqueen_" + str(n) + "_" + str(nqueen) + "_" +  str(k) + ".out") as f:
                    while True:
                        tmp = f.readline().split(' = ')
                        if len(tmp) > 1:
                            parallel_time[i].append(int(tmp[1]))
                            break
    print(sequential_time)
    print(parallel_time)

    colors = ["r", "g", "b", "h", "i"]

    plt.figure()
    for i in range(0, len(nqueens)):
        plt.plot(nodes, parallel_time[i], colors[i] + 'o-')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.yscale("log")
    plt.title(title)
    plt.legend([str(i) + " queens" for i in nqueens])
    plt.savefig(title)

if __name__ == "__main__":
    # draw_speedup("./n2_sequential/1node-1-cpt-1-npn-snowy.out",
    #              "./n2_openmpi/",
    #              "node-1-cpt-1-npn-snowy.out",
    #              [i for i in range(2, 13)],
    #              title="Speed of MPI O(n^2)",
    #              xlabel="n nodes",
    #              ylabel="speedup"
    # )
    
    # draw_nqueen_speedp("./mc_nqueen/",
    #             "./mc_nqueen_parallel/",
    #             [8, 50, 100],
    #             [90],
    #             [i for i in range(2, 13)],
    #             title="Speedup of Minimum conflicts n-queens",
    #             xlabel="n nodes",
    #             ylabel="speedup"
    # )

    draw_nqueen_time("./mc_nqueen/",
                "./mc_nqueen_parallel/",
                [8, 50, 100],
                [90],
                [i for i in range(2, 13)],
                title="Wall time of Minimum conflicts n-queens",
                xlabel="n nodes",
                ylabel="Wall time (ms)"
    )
