/*
 * hw11.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: fuji
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

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
	const int numOfParticles = 100000;
	// Allocate space for position array
	double *loc = new double[numOfParticles * DIM];

	// Allocate space for force vector array
	double *foc = new double[numOfParticles * DIM];

	// Allocate space for velocity vector array
	double *vel = new double[numOfParticles * DIM];

	// Make Distribute particles and set forces
	for (int i = 0; i < numOfParticles; i++) {
		loc[i * DIM] = rand() / RAND_MAX;
		loc[i * DIM + 1] = rand() / RAND_MAX;
		foc[i * DIM] = rand() / RAND_MAX - 0.5;
		foc[i * DIM + 1] = rand() / RAND_MAX - 0.5;
	}

	/// Compute Velocities
	double st = tsecond();
	for (int p = 0; p < numOfParticles; p++) {
		/* zeros */
		vel[p * DIM] = 0.0; //x component
		vel[p * DIM + 1] = 0.0; //y component

		/* loop for particles  */
		int bsize = 1000;
		/* Loop tiling and parallelization */
		for (int ii = 0; ii < numOfParticles; ii+=bsize) {
			#pragma omp parallel for
			for (int j = 0; j < numOfParticles; j++) {
				for (int i = ii; i < std::min(numOfParticles,ii+ bsize); i++){
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
			}
	}
	double et = tsecond();

	// Compute Average Velocity
	double vx = 0.0;
	double vy = 0.0;
	for (int i = 0; i < numOfParticles; i++) {
		vx += vel[i * DIM];
		vy += vel[i * DIM + 1];
	}
	vx /= numOfParticles;
	vy /= numOfParticles;

	// Show Results
	std::cout << "Mean Velocity = (" << vx << "," << vy << ")\n";
	std::cout << "Time cost = " << et - st << "(sec)\n";

	// cleanup
//	free(loc);
	free(vel);
	free(foc);
	return 0;
}

