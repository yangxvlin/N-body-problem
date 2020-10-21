// modify from: https://github.com/BaseMax/N-QueenGenetic
#include <iostream>
#include <algorithm> // sort, random_shuffle
#include <sstream> // ostringstream
// #include <string>
#include <map>
#include <queue> // vector, .push_back()
#include <math.h>
#include <vector>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>

using namespace std;

typedef struct {
	string way;
	int cost;
} individual;

int n_threads=4;
const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;

vector<individual*> population;

// Return current time, for performance measurement
// Adapt from previous assignment
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int fitValue(string way, vector<pair<int, int>> position_ij, int all_positions, int chessBoardSize, int rank, int size) {
	int fitness=(chessBoardSize*(chessBoardSize-1))/2;

    int sum2 = 0;

    char str_buffer[chessBoardSize];
    memcpy(str_buffer, way.c_str(), chessBoardSize);
    MPI_Bcast(str_buffer, chessBoardSize, MPI_CHAR, root, comm);

    for (int i = 1; i < size; ++i) {
        MPI_Recv(&sum2, 1, MPI_INT, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
        fitness -= sum2;
    }
	return fitness;

	// for(int i=0; i<chessBoardSize; i++)
	// 	for(int j=i+1; j<chessBoardSize; j++)
	// 		if((way[i] == way[j]) ||  (i-way[i] == j-way[j]) || (i+way[i] == j+way[j]))
	// 			fitness--;
	// cout << fitness << endl; 
	// return fitness;
}

individual* reproduce(individual *a, individual *b, vector<pair<int, int>> position_ij, int all_positions, int chessBoardSize, int rank, int size) {
	individual *item = new individual;
	int c = rand()%chessBoardSize;
	item->way =
				(a->way).substr(0, c) +
				(b->way).substr(c, chessBoardSize-c+1)
	;
	item->cost = fitValue(item->way, position_ij, all_positions, chessBoardSize, rank, size);
	return item;
}

bool compare(individual *a, individual*b) {
	return(a->cost > b->cost);
}

individual* GA(vector<pair<int, int>> position_ij, int all_positions, int chessBoardSize, int rank, int size) {
	int random1, random2;
	individual *individual1, *individual2, *item;
	bool found=0;
	while(!found) {
		vector<individual*> populationNew;
		for(unsigned int i=0; i<population.size(); i++) {
			sort(population.begin(), population.end(), compare);

			random1 = rand()%population.size() %2;
			random2 = rand()%population.size() %2;

			individual1 = population[random1];
			individual2 = population[random2];

			item = reproduce(individual1, individual2, position_ij, all_positions, chessBoardSize, rank, size);

			if(rand()%2==0)
				item->way[rand()%(chessBoardSize)+1] = (rand()%(chessBoardSize)+1)+48;

			if(fitValue(item->way, position_ij, all_positions, chessBoardSize, rank, size)==((chessBoardSize*(chessBoardSize-1))/2)) {
				found=1;
				return item;
			}
			populationNew.push_back(item);
		}
		//
		population = populationNew;
	}
	return item;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(comm, &rank);
    int size;
    MPI_Comm_size(comm, &size);

    uint64_t start;
    start = GetTimeStamp();
    
    int chessBoardSize;
    chessBoardSize=8;
    
    vector<pair<int, int>> position_ij;
    int all_positions = 0;
    for (int i=0; i<chessBoardSize; ++i) {
        for (int j=i+1; j<chessBoardSize; ++j) {
            position_ij.push_back(make_pair(i, j));
            ++all_positions;
        }
    }
    int tasks_per_rank = ceil((1.0* all_positions) / (size-1)); // master slave

    if (rank == root) {
        srand(time(0));
        
        int maxSolutions;
        int initialPopulationCount = 10;

        maxSolutions=10;
        map<string, int> solutions;
        int numFound=0;


        while(numFound!=maxSolutions) {
            string tempWay="";
            for(int i=1; i<=chessBoardSize; i++) {
                ostringstream ostr;
                ostr<<i;
                tempWay += ostr.str();
            }

            // fitValue(tempWay, position_ij, all_positions, chessBoardSize, rank, size);
            // break;

            individual *temp;
            for(int i=0; i<initialPopulationCount; i++) {
                random_shuffle(tempWay.begin(), tempWay.end());
                temp   = new individual;
                temp->way = tempWay;
                temp->cost = fitValue(tempWay, position_ij, all_positions, chessBoardSize, rank, size);
                population.push_back(temp);
            }

            individual *solution = GA(position_ij, all_positions, chessBoardSize, rank, size);
            if(!solutions[solution->way]) {
                solutions[solution->way]=1;
                cout<<"Possible Solution #"<<(++numFound)<<":\t"<<solution->way<<endl;
            }

        }
        cout << "master finished" << endl;
        char str_buffer[chessBoardSize];
        str_buffer[0] = '#';
        MPI_Bcast(str_buffer, chessBoardSize, MPI_CHAR, root, comm);

        cout << "time = " << GetTimeStamp() - start << endl;
    } else {
        char way[chessBoardSize+1];

        while (true) {
            MPI_Bcast(way, chessBoardSize, MPI_CHAR, root, comm);
            way[chessBoardSize] = '\0';
            // cout << "rank[" << rank << "] " << "str received: " << way[0] << way[1] << way[2] << way[3] << way[4] << endl;

            if (way[0] == '#') {
                // cout << "stop" << endl;
                break;
            }

            int sum2 = 0;

            for (int z = 0; z < tasks_per_rank; ++z) {
            int index = (rank-1) * tasks_per_rank + z;
			
                if (index < all_positions) {
                    pair<int, int> ij = position_ij.at(index);
                    
                    int i = ij.first;
                    int j = ij.second;
                    // cout << "rank[" << rank << "] " << index << "(" << i << " " << j << ")" << endl;
                    
                    if((way[i] == way[j]) || (i-way[i] == j-way[j]) || (i+way[i] == j+way[j])) {
                        sum2 += 1;
                    }
                }
            }

            // cout << "rank[" << rank << "] sum2: " << sum2 << endl;
            MPI_Send(&sum2, 1, MPI_INT, root, 0, comm);
        }

    }

    MPI_Finalize();
	return 0;
}

// mpicxx genetic_nqueen_mpi.cpp -o nqueen -O3
// mpirun -np 2 nqueen
// time mpirun -np 4 nqueen > 8queen.out