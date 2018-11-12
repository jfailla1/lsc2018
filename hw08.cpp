/*
 * hw11.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: fuji
 */
/* 
Mean Velocity = (-24291.7,-24291.7)
Time cost = 858.876(sec)
2 cores on my computer ^^^
*/
/*
1 Core
CPU0:0~100000
Mean Velocity = (-24291.7,-24291.7)
Time cost = 1359.87(sec)

1.58 Speed up
*/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <mpi.h>

// Dimension
#define DIM 2
// Blob Size
#define EPSILON 0.005

// timing method
double tsecond() {
	struct timeval tm;
	double t;
	static int base_sec = 0, base_usec = 0;

	gettimeofday(&tm, NULL);
	if (base_sec == 0 && base_usec == 0) {
		base_sec = tm.tv_sec;
		base_usec = tm.tv_usec;
		t = 0.0;
	} else {
		t = (double) (tm.tv_sec - base_sec) + ((double) (tm.tv_usec - base_usec)) / 1.0e6;
	}
	return t;
}

// function term 1
double term1(double r, double ep) {
#ifdef VALGRIND
	return 1.0;
#else
	double sq;

	sq = sqrt(r * r + ep * ep);

	return log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
#endif
}

// function term 2
double term2(double r, double ep) {
#ifdef VALGRIND
	return 1.0;
#else
	double sq;

	sq = sqrt(r * r + ep * ep);

	return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
#endif
}

// Main Routine
int main(int argc, char **argv) {
	
	// Initialize
  	MPI::Init(argc, argv);
  
  	// get myid and # of processors 
  	int numproc = MPI::COMM_WORLD.Get_size();
  	int myid = MPI::COMM_WORLD.Get_rank();

	const int numOfParticles = 100000;
	// Allocate space for position array
	double *loc = new double[numOfParticles * DIM];

	// Allocate space for force vector array
	double *foc = new double[numOfParticles * DIM];

	// Allocate space for velocity vector array
	double *vel = new double[numOfParticles * DIM];
	
	if (myid == 0){
	// Make Distribute particles and set forces
		for (int i = 0; i < numOfParticles; i++) {
			loc[i * DIM] = rand() / RAND_MAX;
			loc[i * DIM + 1] = rand() / RAND_MAX;
			foc[i * DIM] = rand() / RAND_MAX - 0.5;
			foc[i * DIM + 1] = rand() / RAND_MAX - 0.5;
		}
	}
	
	MPI::COMM_WORLD.Bcast(loc, numOfParticles*DIM, MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(foc, numOfParticles*DIM, MPI::DOUBLE, 0);

	/* divide loop */
        int mystart = (numOfParticles / numproc) * myid;
        int myend;
        if (numOfParticles % numproc > myid) {
                mystart += myid;
                myend = mystart + (numOfParticles / numproc) + 1;
        } else {
                mystart += numOfParticles % numproc;
                myend = mystart + (numOfParticles / numproc);
        }
        std::cout << "CPU" << myid << ":" << mystart << "~" << myend << std::endl;
	double mysize = myend - mystart;
	/// Compute Velocities
	double st = tsecond();
	for (int p = mystart; p < myend; p++) {
		/* zeros */
		vel[p * DIM] = 0.0;
		vel[p * DIM + 1] = 0.0;

		/* loop for particles  */
		for (int i = mystart; i < myend; i++) {
			double dx = loc[p * DIM] - loc[i * DIM];
			double dy = loc[p * DIM + 1] - loc[i * DIM + 1];
			double r = sqrt(dx * dx + dy * dy);

			double tr1 = term1(r, EPSILON) / (4.0 * M_PI);
			double tr2 = term2(r, EPSILON) / (4.0 * M_PI);

			tr2 *= foc[i * DIM] * dx + foc[i * DIM + 1] * dy;

			vel[p * DIM] += -foc[i * DIM] * tr1 + tr2 * dx;
			vel[p * DIM + 1] += -foc[i * DIM + 1] * tr1 + tr2 * dy;
		}
	}
	double et = tsecond();

	// Compute Average Velocity
	double vx = 0.0;
	double vy = 0.0;
	for (int i = mystart; i < myend; i++) {
		vx += vel[i * DIM];
		vy += vel[i * DIM + 1];
	}
	vx /= mysize;
	vy /= mysize;
	
	MPI::COMM_WORLD.Barrier();
	double vx_t, vy_t, et_t, st_t;

	// Send the vars to rank 0
	MPI::COMM_WORLD.Reduce(&vx, &vx_t, 1, MPI::DOUBLE, MPI::SUM, 0);
	MPI::COMM_WORLD.Reduce(&vy, &vy_t, 1, MPI::DOUBLE, MPI::SUM, 0);
	MPI::COMM_WORLD.Reduce(&et, &et_t, 1, MPI::DOUBLE, MPI::SUM, 0);
	MPI::COMM_WORLD.Reduce(&st, &st_t, 1, MPI::DOUBLE, MPI::SUM, 0);

	if (myid == 0){
	// Show Results
	std::cout << "Mean Velocity = (" << vx_t << "," << vy_t << ")\n";
	std::cout << "Time cost = " << et_t - st_t << "(sec)\n";
	}

	// cleanup
	delete [] loc;
	delete [] vel;
	delete [] foc;

	 
  	MPI::Finalize();
	return 0;
}

