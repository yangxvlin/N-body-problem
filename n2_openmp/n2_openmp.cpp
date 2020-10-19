/*
 * O(n^2) openmp paralleled
 * Author: Xulin Yang, 904904 
 * 
 * A body is in 3D x-y-z coordinates with mass and velocity
 */

#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>

#include <omp.h>

using namespace std;

// Return current time, for performance measurement
// Adapt from previous assignment
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

constexpr double EPSILON = 0.000001;
constexpr double X_BOUND = 1.0e6;      // Width of space
constexpr double Y_BOUND = 1.0e6;      // Height of space
constexpr double Z_BOUND = 1.0e6;      // Depth of space
constexpr int N_THREADS = 16;

// Body related calculation
struct Body {
    double mass, px, py, pz, vx, vy, vz;
};
// Overloaded operator for '<<' for struct output
ostream& operator << (ostream& os, const Body& body) {
    os << "body[" << body.px << "][" << body.py << "][" << body.pz << "] velocity: (" << body.vx << ", " << body.vy << ", " << body.vz << ") mass: " << body.mass;
    return os;
}

struct Force {
    double fx, fy, fz;
};

// compute force for body i based on body j where j != i
// i: body i
inline void compute_force(int i, int N, double G, Body *n_bodies, Force *n_bodies_forces){
    // reset force
    n_bodies_forces[i].fx = 0;
    n_bodies_forces[i].fy = 0;
    n_bodies_forces[i].fz = 0;

    double px_diff, py_diff, pz_diff, factor, euclidean_distance;
    for (int j = 0; j < N; j++) {
        if (i != j) {
            // distance in x direction
            px_diff = n_bodies[j].px - n_bodies[i].px;
            // distance in y direction
            py_diff = n_bodies[j].py - n_bodies[i].py;
            // distance in z direction
            pz_diff = n_bodies[j].pz - n_bodies[i].pz;

            // ||p_j - p_i||
            euclidean_distance = sqrt(pow(px_diff, 2) + pow(py_diff, 2) + pow(pz_diff, 2)) + EPSILON;  // add epsilon to avoid zero division

            // G * m_i * m_j / (||p_j - p_i||)^3
            factor = G * n_bodies[j].mass * n_bodies[i].mass / pow(euclidean_distance, 3);
            // f_ij = factor * (p_j - p_i)
            n_bodies_forces[i].fx += px_diff * factor; // force in x direction
            n_bodies_forces[i].fy += py_diff * factor; // force in y direction
            n_bodies_forces[i].fz += pz_diff * factor; // force in z direction
        }
    }
}

// update position and velocity of a body
inline void update_body(Body * body_next, int N, double G, double TIME_DELTA, Body body_cur, Force body_force) {
    // factor = dt / m
    double factor = TIME_DELTA / body_cur.mass;

    body_next->px += TIME_DELTA * body_cur.vx;
    body_next->py += TIME_DELTA * body_cur.vy;
    body_next->pz += TIME_DELTA * body_cur.vz;

    // v = dt * a = dt * (F / m) = (dt / m) * F = factor * F
    body_next->vx += factor * body_force.fx;
    body_next->vy += factor * body_force.fy;
    body_next->vz += factor * body_force.fz;

    // wrap the position if out of bound
    if (body_next->px >= X_BOUND) {
        body_next->px = fmod(body_next->px, X_BOUND);
    } else if (body_next->px <= 0) {
        while (body_next->px < 0) {
            body_next->px += X_BOUND;
        }
    }
    if (body_next->py >= Y_BOUND) {
        body_next->py = fmod(body_next->py, Y_BOUND);
    } else if (body_next->py <= 0) {
        while (body_next->py < 0) {
            body_next->py += Y_BOUND;
        }
    }
    if (body_next->pz >= Z_BOUND) {
        body_next->pz = fmod(body_next->pz, Z_BOUND);
    } else if (body_next->pz <= 0) {
        while (body_next->pz < 0) {
            body_next->pz += Z_BOUND;
        }
    }
}

inline void calculate(int N, int T, double G, double TIME_DELTA, Body *n_bodies) {
    Body n_bodies_next[N];
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        n_bodies_next[i] = n_bodies[i];
    }

    Force n_bodies_forces[N];
    for (int z = 0; z < T; ++z) {
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < N; ++i) {
                compute_force(i, N, G, n_bodies, n_bodies_forces);
            }
            #pragma omp for
            for (int i = 0; i < N; ++i) {
                update_body(&(n_bodies_next[i]), N, G, TIME_DELTA, n_bodies[i], n_bodies_forces[i]);
            }
            #pragma omp for
            for (int i = 0; i < N; i++) {
                n_bodies[i] = n_bodies_next[i];
            }
        }
    }
}

int main(int argc, char **argv) {
    uint64_t start, end;
    if (argc < 4) {
        cout << "at least three numbers for n iterations, gravity constant and time delta" << endl;
        return -1;
    }
    // n iterations
    int T = atoi(argv[1]);
    if (T <= 0) {
        cout << "Input an positive iteration number" << endl;
        return -1;
    }
    // grativity constant
    double G = atoi(argv[2]);
    if (G <= 0) {
        cout << "Input an positive grativity constant" << endl;
        return -1;
    }

    double TIME_DELTA = atoi(argv[3]);
    if (TIME_DELTA <= 0) {
        cout << "Input an positive time delta" << endl;
        return -1;
    }

    int n_threads = omp_get_max_threads();
    cout << "n_threads = " << n_threads << endl;
    omp_set_num_threads(n_threads);

    // n bodies
    int N;
    cin >> N;
    
    Body n_bodies[N];
    for (int i = 0; i < N; ++i) {
        cin >> n_bodies[i].mass;
        cin >> n_bodies[i].px;
        cin >> n_bodies[i].py;
        cin >> n_bodies[i].pz;
        cin >> n_bodies[i].vx;
        cin >> n_bodies[i].vy;
        cin >> n_bodies[i].vz;
    }

    start = GetTimeStamp();
    calculate(N, T, G, TIME_DELTA, n_bodies);
    cout << "time = " << GetTimeStamp() - start << endl;

    cout << N << endl;
    for (int i = 0; i < N; ++i) {
        cout << i << ": " << n_bodies[i] << endl;
    }

    return 0;
}

// g++ -std=c++14 -fopenmp -O3 -o n2_openmp n2_openmp.cpp
// ./n2_openmp 10 6.67e-11 1.0 < ../body_10.data
