/*
 * genrate random 2D body
 * Author: Xulin Yang, 904904 
 */

#include <stdio.h>
#include <iostream>
#include <cstdlib>

using namespace std;

double rand_double( double low, double high ) {
    return ( ( double )rand() * ( high - low ) ) / ( double )RAND_MAX + low;
}

int main(int argc, char **argv){
    if (argc < 2) {
        cout << "Must has a positive integer N" << endl;
        return(-1);
    }
    int N = atoi(argv[1]); // number of body to be generated
    if (N <= 0) {
        cout << "N must be a positive integer!" << endl;
    }

    cout << N << endl;

    double mass, px, py, pz, vx, vy, vz;

    for (int i = 0; i < N; i++) {
        mass = rand_double(0.0, 1.0);
        px   = rand_double(0.0, 1.0);
        py   = rand_double(0.0, 1.0);
        pz   = rand_double(0.0, 1.0);
        vx   = rand_double(0.0, 1.0);
        vy   = rand_double(0.0, 1.0);
        vz   = rand_double(0.0, 1.0); 
        printf("%lf %lf %lf %lf %lf %lf %lf\n", mass, px, py, pz, vx, vy, vz);
    }
    
    return 0;
}

// g++ -std=c++14 -o random_body random_body.cpp
// ./random_body 10 > body_10.data