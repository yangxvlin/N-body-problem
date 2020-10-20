// modify from: https://github.com/BaseMax/N-QueenGenetic
#include <iostream>
#include <algorithm> // sort, random_shuffle
#include <sstream> // ostringstream
// #include <string>
#include <map>
#include <queue> // vector, .push_back()
#include <omp.h>
#include <math.h>
#include <vector>

using namespace std;

typedef struct {
	string way;
	int cost;
} individual;

vector<individual*> population;
int chessBoardSize;
int maxSolutions;
int initialPopulationCount = 10;

int n_threads=4;

int fitValue(string way, vector<pair<int, int>> position_ij, int all_positions) {
	int fitness=(chessBoardSize*(chessBoardSize-1))/2;
    
	omp_set_num_threads(n_threads);
    int sum2 = 0;
    int tasks_per_thread = ceil((1.0* all_positions) / n_threads);
    #pragma omp parallel num_threads(n_threads) reduction(+:sum2)
    {
        int thread_rank = omp_get_thread_num();
        for (int z = 0; z < tasks_per_thread; ++z) {
            int index = thread_rank * tasks_per_thread + z;
			#pragma omp critical
			
			
			if (index < all_positions) {
				pair<int, int> ij = position_ij.at(index);
				
				int i = ij.first;
				int j = ij.second;
				
				if((way[i] == way[j]) || (i-way[i] == j-way[j]) || (i+way[i] == j+way[j])) {
					sum2 += 1;
				}
			}
        }

    }
	return fitness - sum2;

	// for(int i=0; i<chessBoardSize; i++)
	// 	for(int j=i+1; j<chessBoardSize; j++)
	// 		if((way[i] == way[j]) ||  (i-way[i] == j-way[j]) || (i+way[i] == j+way[j]))
	// 			fitness--;
	// cout << fitness << endl; 
	// return fitness;
}

individual* reproduce(individual *a, individual *b, vector<pair<int, int>> position_ij, int all_positions) {
	individual *item = new individual;
	int c = rand()%chessBoardSize;
	item->way =
				(a->way).substr(0, c) +
				(b->way).substr(c, chessBoardSize-c+1)
	;
	item->cost = fitValue(item->way, position_ij, all_positions);
	return item;
}

bool compare(individual *a, individual*b) {
	return(a->cost > b->cost);
}

individual* GA(vector<pair<int, int>> position_ij, int all_positions) {
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

			item = reproduce(individual1, individual2, position_ij, all_positions);

			if(rand()%2==0)
				item->way[rand()%(chessBoardSize)+1] = (rand()%(chessBoardSize)+1)+48;

			if(fitValue(item->way, position_ij, all_positions)==((chessBoardSize*(chessBoardSize-1))/2)) {
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

int main() {
	srand(time(0));
	chessBoardSize=8;
	maxSolutions=92;
	map<string, int> solutions;
	int numFound=0;

    vector<pair<int, int>> position_ij;
    int all_positions = 0;
    for (int i=0; i<chessBoardSize; ++i) {
		for (int j=i+1; j<chessBoardSize; ++j) {
            position_ij.push_back(make_pair(i, j));
            ++all_positions;
        }
    }

	while(numFound!=maxSolutions) {
		string tempWay="";
		for(int i=1; i<=chessBoardSize; i++) {
			ostringstream ostr;
			ostr<<i;
			tempWay += ostr.str();
		}

		// fitValue(tempWay, position_ij, all_positions);
		// break;

		individual *temp;
		for(int i=0; i<initialPopulationCount; i++) {
			random_shuffle(tempWay.begin(), tempWay.end());
			temp   = new individual;
			temp->way = tempWay;
			temp->cost = fitValue(tempWay, position_ij, all_positions);
			population.push_back(temp);
		}

		individual *solution = GA(position_ij, all_positions);
		if(!solutions[solution->way]) {
			solutions[solution->way]=1;
			cout<<"Possible Solution #"<<(++numFound)<<":\t"<<solution->way<<endl;
		}
	}
	return 0;
}

// g++ genetic_nqueen_parallel.cpp -o nqueen -O3 -fopenmp
// ./nqueen