/*
 * genrate random 3D body
 * Author: Xulin Yang, 904904 
 */

#include <stdio.h>
#include <iostream>
#include <cstdlib>

using namespace std;

constexpr double X_BOUND = 5.9e12;      // Width of solar system
constexpr double Y_BOUND = 5.9e12;      // Height of solar system
constexpr double Z_BOUND = 5.9e12;      // Depth of solar system
constexpr double MAX_SPEED = 50000;     // orbit speed of mercury https://nssdc.gsfc.nasa.gov/planetary/factsheet/
constexpr double MAX_MASS = 1.0e30;     // Mass of the sun

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
    int T = 1440; // 24h * 60 minutes
    cout << T << endl;
    double G = 6.67e-11 ;
    cout << G << endl;
    double TIME_DELTA = 60; // 1 mintes
    cout << TIME_DELTA << endl;

    double mass, px, py, pz, vx, vy, vz;

    for (int i = 0; i < N; i++) {
        mass = rand_double(0.0, 1.0) * MAX_MASS;
        px   = rand_double(0.0, 1.0) * X_BOUND;
        py   = rand_double(0.0, 1.0) * Y_BOUND;
        pz   = rand_double(0.0, 1.0) * Z_BOUND;
        vx   = rand_double(0.0, 1.0) * MAX_SPEED;
        vy   = rand_double(0.0, 1.0) * MAX_SPEED;
        vz   = rand_double(0.0, 1.0) * MAX_SPEED; 
        printf("%lf %lf %lf %lf %lf %lf %lf\n", mass, px, py, pz, vx, vy, vz);
    }
    
    return 0;
}

// g++ -std=c++14 -o random_body random_body.cpp
// ./random_body 10 > body_10.data
// ./random_body 100 > body_100.data
// ./random_body 1000 > body_1000.data
// ./random_body 2000 > body_2000.data
// ./random_body 5000 > body_5000.data