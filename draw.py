# Author: Xulin Yang, 904904
# draw analysis diagrams

from matplotlib import pyplot as plt

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
    plt.plot(parallel_outputs, speedup)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(title)

if __name__ == "__main__":
    draw_speedup("./n2_sequential/1node-1-cpt-1-npn-snowy.out",
                 "./n2_openmpi/",
                 "node-1-cpt-1-npn-snowy.out",
                 [i for i in range(2, 12)],
                 title="Speed of MPI O(n^2)",
                 xlabel="n nodes",
                 ylabel="speedup"
    )
