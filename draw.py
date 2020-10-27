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

def get_speedup(sequential_directory: str, 
                 parallel_directory: str, 
                 nodes: list,
                 bodies: list
                ):
    sequential_time = []
    parallel_time = [[] for _ in range(len(bodies))]
    for i, b in enumerate(bodies):

        with open(sequential_directory + "1-1-{}".format(b) + ".out") as f:
            seq_time = int(f.readline().split(' = ')[1])
        sequential_time.append(seq_time)

        for n in nodes:
            with open(parallel_directory + "{}-1-{}".format(n, b) + ".out") as f:
                para_time = int(f.readline().split(' = ')[1])
            parallel_time[i].append(para_time)
    
    speedup = [[] for _ in range(len(bodies))]
    for i, pts in enumerate(parallel_time):
        for pt in pts:
            speedup[i].append(sequential_time[i] / pt)
    print(sequential_time)
    print(parallel_time)
    print(speedup)
    return sequential_time, parallel_time, speedup

def get_speedup2(sequential_directory: str, 
                 parallel_directory: str, 
                 sequential_node: int,
                 parallel_node: int,
                 bodies: list,
                 cores: list
                ):
    sequential_time = []
    parallel_time = [[] for _ in range(len(bodies))]
    for i, b in enumerate(bodies):

        with open(sequential_directory + "{}-1-{}".format(sequential_node, b) + ".out") as f:
            seq_time = int(f.readline().split(' = ')[1])
        sequential_time.append(seq_time)

        for c in cores:
            with open(parallel_directory + "{}-{}-{}".format(parallel_node, c, b) + ".out") as f:
                para_time = int(f.readline().split(' = ')[1])
            parallel_time[i].append(para_time)
    
    speedup = [[] for _ in range(len(bodies))]
    for i, pts in enumerate(parallel_time):
        for pt in pts:
            speedup[i].append(sequential_time[i] / pt)
    print(sequential_time)
    print(parallel_time)
    print(speedup)
    return sequential_time, parallel_time, speedup

def get_profile(parallel_directory: str, 
                nodes: list,
                body: int,
                hasTree=False):
    if hasTree:
        n_input = 6
    else:
        n_input = 4
    
    profile = [[] for _ in range(n_input)]

    for n in nodes:
        with open(parallel_directory + "{}-1-{}".format(n, body) + ".out") as f:

            for j in range(0, n_input):
                tmp = int(f.readline().split(' = ')[1])
                profile[j].append(tmp)
    
    print(profile)
    return profile

def draw(xs, yss, title, xlabel, ylabel, legend):
    plt.figure()
    for ys in yss:
        plt.plot(xs, ys, 'o-')
    plt.legend(legend)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.yscale('log')
    plt.title(title)
    plt.savefig(title)

if __name__ == "__main__":
    n2_profile = get_profile("n2_openmpi_profile/",
                             [i for i in range(2, 13)],
                             2000,
                             False
                            )
    draw([i for i in range(2, 13)],
         n2_profile,
         "parallel O(n^2) Runtime profile",
         "nodes",
         "Runtime (ms)",
         ["communication", "force calculation", "body update", "other"]
    )
    nlogn_profile = get_profile("nlogn_openmpi_profile/",
                             [i for i in range(2, 13)],
                             2000,
                             True
                            )
    draw([i for i in range(2, 13)],
         nlogn_profile,
         "parallel O(n logn) Runtime profile",
         "nodes",
         "Runtime (ms)",
         ["communication", "tree construct", "force calculation", "body update", "tree delete", "other"]
    )

    # n2_seq, n2_para, n2_speedup = get_speedup("n2_sequential/",
    #             "n2_openmpi/",
    #             [i for i in range(2, 13)],
    #             [10, 100, 500, 1000, 2000]
    #             )
    # nlogn_seq, nlogn_para, nlogn_speedup = get_speedup("nlogn_sequential/",
    #             "nlogn_openmpi/",
    #             [i for i in range(2, 13)],
    #             [10, 100, 500, 1000, 2000]
    #             )
    # draw([i for i in range(2, 13)],
    #      n2_speedup,
    #      "O(n^2) speedup",
    #      "nodes",
    #      "speedup",
    #      list(map(lambda x: "{}-body".format(x), [10, 100, 500, 1000, 2000]))
    # )
    # draw([i for i in range(2, 13)],
    #      nlogn_speedup,
    #      "O(n logn) speedup",
    #      "nodes",
    #      "speedup",
    #      list(map(lambda x: "{}-body".format(x), [10, 100, 500, 1000, 2000]))
    # )
    # draw([i for i in range(2, 13)],
    #      n2_para,
    #      "O(n^2) Runtime",
    #      "nodes",
    #      "Runtime (ms)",
    #      list(map(lambda x: "{}-body".format(x), [10, 100, 500, 1000, 2000]))
    # )
    # draw([i for i in range(2, 13)],
    #      nlogn_para,
    #      "O(n logn) Runtime",
    #      "nodes",
    #      "Runtime (ms)",
    #      list(map(lambda x: "{}-body".format(x), [10, 100, 500, 1000, 2000]))
    # )
    # n2_seq, n2_para, n2_speedup = get_speedup2("n2_openmpi/",
    #         "n2_hybrid/",
    #         5,
    #         5,
    #         [500, 1000, 2000],
    #         [i for i in range(2, 17)],
    #         )
    # draw([i for i in range(2, 17)],
    #     n2_speedup,
    #     "O(n^2) speedup openmpi+openmp 5 nodes",
    #     "threads",
    #     "speedup",
    #     list(map(lambda x: "{}-body".format(x), [500, 1000, 2000]))
    # )
    # n2_seq, n2_para, n2_speedup = get_speedup2("nlogn_openmpi/",
    #         "nlogn_hybrid/",
    #         5,
    #         5,
    #         [500, 1000, 2000],
    #         [i for i in range(2, 17)],
    #         )
    # draw([i for i in range(2, 17)],
    #     n2_speedup,
    #     "O(n logn) speedup openmpi+openmp 5 nodes",
    #     "threads",
    #     "speedup",
    #     list(map(lambda x: "{}-body".format(x), [500, 1000, 2000]))
    # )
    # n2_seq, n2_para, n2_speedup = get_speedup2("n2_sequential/",
    #         "n2_hybrid/",
    #         1,
    #         5,
    #         [500, 1000, 2000],
    #         [i for i in range(2, 17)],
    #         )
    # draw([i for i in range(2, 17)],
    #     n2_speedup,
    #     "O(n^2) speedup openmpi+openmp 5 nodes vs sequential",
    #     "threads",
    #     "speedup",
    #     list(map(lambda x: "{}-body".format(x), [500, 1000, 2000]))
    # )
    # n2_seq, n2_para, n2_speedup = get_speedup2("nlogn_sequential/",
    #         "nlogn_hybrid/",
    #         1,
    #         5,
    #         [500, 1000, 2000],
    #         [i for i in range(2, 17)],
    #         )
    # draw([i for i in range(2, 17)],
    #     n2_speedup,
    #     "O(n logn) speedup openmpi+openmp 5 nodes vs sequential",
    #     "threads",
    #     "speedup",
    #     list(map(lambda x: "{}-body".format(x), [500, 1000, 2000]))
    # )

    # draw_speedup("./n2_sequential/1node-1-cpt-1-npn-snowy.out",
    #              "./n2_openmpi/",
    #              "node-1-cpt-1-npn-snowy.out",
    #              [i for i in range(2, 13)],
    #              title="Speed of MPI O(n^2)",
    #              xlabel="n nodes",
    #              ylabel="speedup"
    # )
