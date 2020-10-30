/*
 * proposed hybrid O(n^2) algorithm's runtime profile
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
            euclidean_distance = sqrt(pow(px_diff, 2) + pow(py_diff, 2) + pow(pz_diff, 2));

            // G * m_i * m_j / (||p_j - p_i||)^3
            factor = G * n_bodies[j].mass * n_bodies[i].mass / (pow(euclidean_distance, 3) + EPSILON); // + epsilon to avoid zero division
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
    uint64_t comm_time = 0, compute_force_time = 0, update_time = 0, tmp1;

    tmp1 = GetTimeStamp();
    MPI_Bcast(&T, 1, MPI_INT, root, comm);
    MPI_Bcast(&G, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&TIME_DELTA, 1, MPI_DOUBLE, root, comm);

    MPI_Datatype MPI_Body = create_MPI_Body();
    MPI_Datatype MPI_Force = create_MPI_Force();
    comm_time += GetTimeStamp() - tmp1;

    int rank;
    MPI_Comm_rank(comm, &rank);
    int size;
    MPI_Comm_size(comm, &size);
    
    int n_per_rank = (int) ceil((1.0 * N) / size);
    // workload for each rank calculation
    int n_start = rank * n_per_rank;
    int n_end   = min(N, (rank+1) * n_per_rank);
    int workload = n_end - n_start;
    Force tmp_forces[n_per_rank];
    Body  tmp_n_bodies[n_per_rank];
    
    Force n_bodies_forces[n_per_rank * size];
    Body  n_nodies_padded[n_per_rank * size];
    if (rank == root) {
        for (int i = 0; i < N; ++i) {
            n_nodies_padded[i] = n_bodies[i];
        }
    }

    int n_threads = omp_get_max_threads();
    omp_set_num_threads(n_threads);

    MPI_Bcast(n_nodies_padded, N, MPI_Body, root, comm);

    for (int z = 0; z < T; ++z) {
        #pragma omp parallel for
        for (int i = n_start; i < n_end; ++i) {
            tmp_n_bodies[i - n_start] = n_nodies_padded[i];
        }

        tmp1 = GetTimeStamp();
        #pragma omp parallel for
        for (int i = n_start; i < n_end; ++i) {
            compute_force(i, N, G, n_nodies_padded, &(tmp_forces[i - n_start]));
        }
        compute_force_time += GetTimeStamp() - tmp1;

        tmp1 = GetTimeStamp();
        MPI_Allgather(tmp_forces,      n_per_rank, MPI_Force,
                      n_bodies_forces, n_per_rank, MPI_Force,
                      comm);
        comm_time += GetTimeStamp() - tmp1;
        tmp1 = GetTimeStamp();
        #pragma omp parallel for
        for (int i = n_start; i < n_end; ++i) {
            update_body(&(tmp_n_bodies[i - n_start]), N, G, TIME_DELTA, n_nodies_padded[i], n_bodies_forces[i]);
        }
        update_time += GetTimeStamp() - tmp1;
        tmp1 = GetTimeStamp();
        MPI_Allgather(tmp_n_bodies,    n_per_rank, MPI_Body,
                      n_nodies_padded, n_per_rank, MPI_Body,
                      comm);
        comm_time += GetTimeStamp() - tmp1;
    }

    if (rank == root) {
        cout << "comm_time = " << comm_time << endl;
        cout << "cal_force_time = " << compute_force_time << endl;
        cout << "body_update_time = " << update_time << endl;
    }

    if (rank == root) {
        for (int i = 0; i < N; ++i) {
            n_bodies[i] = n_nodies_padded[i];
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
    MPI_Bcast(&N, 1, MPI_INT, root, comm);
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

// mpicxx -std=c++14 -fopenmp -O3 -o n2_openmpi n2_openmpi.cpp
// mpirun -np 4 n2_openmpi < ../body_10.data
