/*
 * O(n^2) sequential
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
#include <mpi.h>

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

constexpr double EPSILON = 0.000001;

// Body related calculation
struct Body {
    double mass, px, py, pz, vx, vy, vz;
};
// Overloaded operator for '<<' for struct output
ostream& operator << (ostream& os, const Body& body) {
    os << "body[" << body.px << "][" << body.py << "][" << body.pz << "] velocity: (" << body.vx << ", " << body.vy << ", " << body.vz << ") mass: " << body.mass;
    return os;
}
inline MPI_Datatype create_MPI_Body() {
    MPI_Datatype MPI_Body;
    MPI_Type_contiguous(7, MPI_DOUBLE, &MPI_Body);
    MPI_Type_commit(&MPI_Body);
    return MPI_Body;
}

struct Force {
    double fx, fy, fz;
};
inline MPI_Datatype create_MPI_Force() {
    MPI_Datatype MPI_Force;
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_Force);
    MPI_Type_commit(&MPI_Force);
    return MPI_Force;
}

// compute force for body i based on body j where j != i
// i: body i
inline void compute_force(int i, int N, double G, Body *n_bodies, Force * force){
    // reset force
    force->fx = 0;
    force->fy = 0;
    force->fz = 0;

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
            force->fx += px_diff * factor; // force in x direction
            force->fy += py_diff * factor; // force in y direction
            force->fz += pz_diff * factor; // force in z direction
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
}

inline void calculate(int N, int T, double G, double TIME_DELTA, Body *n_bodies) {
    MPI_Bcast(&N, 1, MPI_INT, root, comm);
    MPI_Bcast(&T, 1, MPI_INT, root, comm);
    MPI_Bcast(&G, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&TIME_DELTA, 1, MPI_DOUBLE, root, comm);

    MPI_Datatype MPI_Body = create_MPI_Body();
    MPI_Datatype MPI_Force = create_MPI_Force();
    int rank;
    MPI_Comm_rank(comm, &rank);
    int size;
    MPI_Comm_size(comm, &size);
    
    int n_per_rank = (int) ceil((1.0 * N) / size);
    // workload for each rank calculation
    int n_start = rank * n_per_rank;
    int n_end   = min(N, (rank+1) * n_per_rank);
    int workload = n_end - n_start;
    Force tmp_forces[workload];
    Body  tmp_n_bodies[workload];
    for (int i = n_start; i < n_end; ++i) {
        tmp_n_bodies[i - n_start].mass = n_bodies[i].mass;
    }
    
    Force n_bodies_forces[N];

    if (rank == root) {
        cout << "test 0 " << rank << endl;
    }

    for (int z = 0; z < T; ++z) {
        MPI_Bcast(n_bodies, N, MPI_Body, root, comm);
        for (int i = n_start; i < n_end; ++i) {
            compute_force(i, N, G, n_bodies, &(tmp_forces[i - n_start]));
        }
        if (rank == root) {
            cout << z << " force computed" << endl;
        }
        MPI_Gather(tmp_forces,      workload, MPI_Force,
                   n_bodies_forces, N,        MPI_Force,
                   root, comm);
        if (rank == root) {
            cout << z << " force gathered" << endl;
        }
        for (int i = n_start; i < n_end; ++i) {
            update_body(&(tmp_n_bodies[i - n_start]), N, G, TIME_DELTA, n_bodies[i], n_bodies_forces[i]);
        }
        if (rank == root) {
            cout << z << " body updated" << endl;
        }
        MPI_Gather(tmp_n_bodies, workload, MPI_Body,
                   n_bodies,     N,        MPI_Body,
                   root, comm);
        if (rank == root) {
            cout << z << " iter finished" << endl;
            cout << "0: " << n_bodies[0] << endl;
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(comm, &rank);

    // n bodies
    int N;
    // n iterations
    int T;
    // frativity const
    double G;
    // time delta
    double TIME_DELTA;

    uint64_t start, end;

    if (rank == root) {
        cin >> N;
        cin >> T;
        cin >> G;
        cin >> TIME_DELTA;
    }

    Body n_bodies[N];
    if (rank == root) {
        for (int i = 0; i < N; ++i) {
            cin >> n_bodies[i].mass;
            cin >> n_bodies[i].px;
            cin >> n_bodies[i].py;
            cin >> n_bodies[i].pz;
            cin >> n_bodies[i].vx;
            cin >> n_bodies[i].vy;
            cin >> n_bodies[i].vz;
        }
    }

    cout << "start0 " << rank << endl;
    start = GetTimeStamp();
    calculate(N, T, G, TIME_DELTA, n_bodies);
    if (rank == root) {
        cout << "time = " << GetTimeStamp() - start << endl;

        cout << N << endl;
        for (int i = 0; i < N; ++i) {
            cout << i << ": " << n_bodies[i] << endl;
        }
    }

    MPI_Finalize();
    return 0;
}

// mpicxx -std=c++14 -O3 -o n2_openmpi n2_openmpi.cpp
// mpirun -np 4 n2_openmpi < ../body_10.data
