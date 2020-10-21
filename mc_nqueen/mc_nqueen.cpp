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

using namespace std;

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

int main() {
	uint64_t start = GetTimeStamp();
	unordered_set <string> answerSet;

	// number of solution expected
	int numExpected;
	
	// generate random time seeds
	srand(time(NULL));
	nQueen r;
	cin >> r.n;
	r.state = new int[r.n];
	cin >> numExpected;
	
	// find k solutions
	int pre_ans_size = answerSet.size(), post_ans_size;
	while (answerSet.size() < numExpected) {
		r.start();
		pre_ans_size = answerSet.size();
		answerSet.insert(r.stateToString());
		post_ans_size = answerSet.size();
		if (post_ans_size != pre_ans_size) {
			cout << post_ans_size << ": " << r.stateToString() << endl;
		}
	}
	delete[] r.state;

	cout << "time = " << GetTimeStamp() - start << endl;
	return 0;
}

// g++ -O3 mc_nqueen.cpp -o nqueen
// ./nqueen < test1.data