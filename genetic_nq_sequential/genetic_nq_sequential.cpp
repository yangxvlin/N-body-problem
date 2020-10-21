// modify from: https://kushalvyas.github.io/gen_8Q.html
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
#include <unordered_set>

using namespace std;

typedef struct {
	vector<int> queen_pos;
	float fitness;
} Individual;

float cal_fitness(vector<int> queen_pos, int N) {
    int clashes = 0;
    // calculate row and column clashes
    unordered_set<int> s( queen_pos.begin(), queen_pos.end() );
    int row_col_clashes = abs(N - (int) s.size());

    int dx, dy;
    // calculate diagonal clashes
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                dx = abs(i-j);
                dy = queen_pos.at(i) - queen_pos.at(j);
                if (dx == dy) {
                    clashes++;
                }
            }
        }
    }

    // fitness = 1 / (1 + clashes) from https://aljrico.github.io/blog/genetic-algorithms/ for n queens
    return 1.0 / (1.0 + clashes);
}

Individual* reproduce(Individual *a, Individual *b, int N) {
	Individual *item = new Individual;
	int c = rand()%N;

    (item->queen_pos).insert(item->queen_pos.end(), a->queen_pos.begin(), a->queen_pos.begin()+c);
    (item->queen_pos).insert(item->queen_pos.end(), b->queen_pos.begin()+c, b->queen_pos.end());

	item->fitness = cal_fitness(item->queen_pos, N);
	return item;
}

bool compare(Individual *a, Individual*b) {
	return(a->fitness > b->fitness);
}

int main(int argc, char** argv) {
    srand(time(NULL));

    int N;
    cin >> N;
    int populationSize;
    cin >> populationSize;
    int crossoverRate;
    cin >> crossoverRate;
    int mutationRate;
    cin >> mutationRate;
    int nSolnExpected;
    cin >> nSolnExpected;

    // int tournamentSize = populationSize / 5;

    int initialPopulationSize=10;

    int nSolnFound = 0;
    vector<Individual*> population;
    while (nSolnFound < nSolnExpected) {
        // randomize start
        Individual *temp;
        for (int i = 0; i < initialPopulationSize; i++) {
            temp   = new Individual;
            for (int j = 0; j < N; j++) { 
                temp->queen_pos.push_back(j);
            }
            random_shuffle(temp->queen_pos.begin(), temp->queen_pos.end());
            temp->fitness=cal_fitness(temp->queen_pos, N);
            population.push_back(temp);
        }

        int random1, random2;
        Individual *individual1, *individual2, *item;
        bool found=0;
        while(!found) {
            vector<Individual*> populationNew;
            for(unsigned int i=0; i<population.size(); i++) {
                // sort by fitness
                sort(population.begin(), population.end(), compare);

                random1 = rand()%population.size() %2;
                random2 = rand()%population.size() %2;

                individual1 = population.at(random1);
                individual2 = population.at(random2);

                item = reproduce(individual1, individual2, N);

                if(rand()%2==0)
                    item->queen_pos[rand()%(N)+1] = (rand()%(N)+1)+48;

                if(cal_fitness(item->queen_pos, N)==1.0) {
                    found=1;
                    return item;
                }
                populationNew.push_back(item);
            }
            population = populationNew;
        }
    }

    return 0;
}

// g++ genetic_nq_sequential.cpp -o nqueen -O3
// time ./nqueen test1.data