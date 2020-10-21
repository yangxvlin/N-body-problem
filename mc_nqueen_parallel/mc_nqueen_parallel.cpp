// modify from https://github.com/pancr9/nQueen/blob/master/nQueenSingleDimension.cpp
// Xulin Yang, 904904
// COMP90025 Assignment 03, 2020 S2

#include <stdio.h>
#include <sstream>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unordered_set>
#include <mpi.h>
#include <omp.h>

using namespace std;

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;

// Return current time, for performance measurement
// Adapt from previous assignment
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

class nQueen {
public:	
	// Checks if a solution found by Minimum Conflicts Algorithm
	bool minConflictsHasSolution;
	// n queen
	int n;
	// temp n clashes in a state
	int tmpNConflicts;
	// store location of Queen in 1d array representation
	int* state;
	// stores value of best possible next state row and column values
	int bestCol, bestRow;
	// stores value of currentConflicts & best Possible conflict
	int currentConflict, bestConflict;
	// stores values of row and column for current state
	int currentRow, currentCol;
	// function to calculate #clashes for a state
	int calcClashes(int[]);
	// function to start searching a solution of nQueen Problem
	void start();
	// function to display a state in 2d on screen
	void displayState();
	// function to generate a random state
	void randomState();
	// function to evaluate minimum conflicts
	void minConflicts();
	// function to display a state in 1d
	void displayStateSimple();
	// function to convert a state to comma separated string of int(s)
	string stateToString();
};
void nQueen::start() {
	minConflictsHasSolution = false;
	minConflicts();	//Minimum Conflicts Algorithm
}
string intsToString(int * nums, int n) {
    string ans = "";
	stringstream ss;  
  	
	for (int i = 0; i < n; ++i) {
		ss << nums[i] << ",";
	}

	ss >> ans;
	return ans;
}

string nQueen::stateToString() {
	return intsToString(state, n);
}
void nQueen::randomState() {
	for (int i = 0; i < n; i++) {
		state[i] = (rand() % n);
	}
}
void nQueen::displayState() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << " ";
			if (j != state[i]) {
				cout << " _ ";

			} else {
				cout << " Q ";
			}
		}
		cout << endl;
	}
}
void nQueen::displayStateSimple() {
	for (int i = 0; i < n; ++i) {
		cout << state[i] << " ";
	}
	cout << endl;
}
int nQueen::calcClashes(int a[]) {
	int clashes = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			// calculate clashes in row, column, diagonal direction
			if ((a[i] == a[j]) || (abs(i - j) == abs(a[i] - a[j]))) {
				clashes++;
			}
		}
	}
	return clashes;
}
void nQueen::minConflicts() {
	while (!minConflictsHasSolution) {
		int maxelements = 10 * n;
		
		randomState();
		
		currentConflict = calcClashes(state);
		
		// solution found
		if (currentConflict == 0) {
			minConflictsHasSolution = true;
		}				

		// until no conflicts exists
		while (currentConflict > 0) {
			bestConflict = currentConflict;
			
			currentRow = rand() % n;
			
			currentCol = state[currentRow];
			for (int j = 0; j < n; j++) {
				state[currentRow] = j;
				tmpNConflicts = calcClashes(state);
				// updating best successor
				if (tmpNConflicts <= bestConflict) {							
					bestConflict = tmpNConflicts;
					bestRow = j;
					bestCol = currentRow;
				}
			}
		
			// To confirm max iterations are performed.
			if (maxelements < 1) {
				break;
			
			} else {
				state[bestCol] = bestRow;
				currentConflict = bestConflict;
				maxelements--;
			}
			// Solution found
			if (bestConflict == 0) {
				minConflictsHasSolution = true;
				break;
			}
		}
	}
}

int main(int argc, char **argv) {
	int prov;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov);
    int rank;
	MPI_Comm_rank(comm, &rank);
    int size;
    MPI_Comm_size(comm, &size);

	unordered_set <string> answerSet;

	// number of solution expected
	int numExpected;
    int N;
	
    if (rank == root) {
	    uint64_t start = GetTimeStamp();
        
        cin >> N;
        cin >> numExpected;
        // bcast n to all nodes
        MPI_Bcast(&N, 1, MPI_INT, root, comm);

        #pragma omp parallel num_threads(2) 
        {
            int thread_rank = omp_get_thread_num();

            if (thread_rank == 0) {
				MPI_Status status;
                int buffer[N];
                int response;

				bool firstTime = true;

                // find k solutions
                int preSize, postSize;
                int nWorkersRemaining = size;
				while (nWorkersRemaining > 0) {
                    // recv from slave
					preSize = answerSet.size();
					MPI_Recv(buffer, N, MPI_INT, MPI_ANY_SOURCE, 0, comm, &status);
					
					// cout << ": " << r.stateToString() << endl;
					string soln = intsToString(buffer, N);
					answerSet.insert(soln);
					postSize = answerSet.size();
					if ((postSize != preSize) && (postSize <= numExpected)) {
						cout << postSize << ": " << soln << endl;
					}
                    
                    // send continue or stop
                    if (postSize >= numExpected) {
                        // stop
                        response = 0;
						--nWorkersRemaining;
						// cout << "master finished" << endl;
                    } else {
                        // continue
                        response = 1;
                    }

					if (postSize >= numExpected && firstTime) {
						cout << "time = " << GetTimeStamp() - start << endl;
						firstTime= false;
					}
					MPI_Send(&response, 1, MPI_INT, status.MPI_SOURCE, 1, comm);
                }
                // cout << "master finished" << endl;
            } else if (thread_rank == 1) {
                MPI_Bcast(&N, 1, MPI_INT, root, comm);
                nQueen r;
                r.n = N;
                r.state = new int[r.n];

                int working = 1;   
                while (working == 1) {
                    r.start();
                    MPI_Send(r.state, N, MPI_INT, root, 0, comm);
					MPI_Recv(&working, 1, MPI_INT, root, 1, comm, MPI_STATUS_IGNORE);
                }
                delete[] r.state;
                // cout << "master's worker finished" << endl;
            }
        }
    } else {
        MPI_Bcast(&N, 1, MPI_INT, root, comm);
        nQueen r;
        r.n = N;
        r.state = new int[r.n];

        int working = 1;   
        while (working == 1) {
            r.start();
            MPI_Send(r.state, N, MPI_INT, root, 0, comm);
            MPI_Recv(&working, 1, MPI_INT, root, 1, comm, MPI_STATUS_IGNORE);
        }
	    delete[] r.state;
        // cout << "worker finished" << endl;
    }
	
    MPI_Finalize();
	return 0;
}

// mpicxx -O3 -fopenmp mc_nqueen_parallel.cpp -o nqueen
// mpirun nqueen < test1.data
// mpirun -np 2 nqueen < test1.data