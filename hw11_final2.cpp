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
#pragma omp declare target
double term1(double r, double ep) {
#ifdef VALGRIND
	return 1.0;
#else
	double sq;

	sq = sqrt(r * r + ep * ep);

	return log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
#endif
}
#pragma omp end declare target

// function term 2
#pragma omp declare target
double term2(double r, double ep) {
#ifdef VALGRIND
	return 1.0;
#else
	double sq;

	sq = sqrt(r * r + ep * ep);

	return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
#endif
}
#pragma omp end declare target

// Main Routine
int main(int argc, char **argv) {
	
	// Gets Number of Devices
	int num_dev = omp_get_num_devices();
  	std::cout << "number of devices " << num_dev << std::endl;

	const int numOfParticles = 100000;
	const int num = numOfParticles;
	// Allocate space for position array
	double *loc = new double[numOfParticles * DIM];

	// Allocate space for force vector array
	double *foc = new double[numOfParticles * DIM];

	// Allocate space for velocity vector array
	double *vel = new double[numOfParticles * DIM];

	// Make Distribute particles and set forces
	for (int i = 0; i < numOfParticles; i++) {
		srand(time(NULL));
		loc[i * DIM] = static_cast<double>(rand()) / RAND_MAX;
		loc[i * DIM + 1] = static_cast<double>(rand()) / RAND_MAX;
		foc[i * DIM] = static_cast<double>(rand()) / RAND_MAX - 0.5;
		foc[i * DIM + 1] = static_cast<double>(rand()) / RAND_MAX - 0.5;
		vel[i * DIM] = 0.0;
		vel[i * DIM + 1] = 0.0;
	}
	double st = tsecond();
	
	int numproc = num_dev + 1;
	#pragma omp parallel num_threads(numproc)
	#pragma omp single
	{
	for (int dev = 0; dev < numproc; dev++) {
	#pragma omp task firstprivate(dev)		
	{
	/* divide loop */
	int mystart = (num / numproc) * dev;
	int myend;
	if (num % numproc > dev) {
	  mystart += dev;
	  myend = mystart + (num / numproc) + 1;
	} else {
	  mystart += num % numproc;
	  myend = mystart + (num / numproc);
	}
	int mysize = myend - mystart;
	double *vel_dev = new double[mysize * DIM];

        #pragma omp target if(dev != num_dev) device(dev) map(to:foc[0:numOfParticles * DIM]) map(to:loc[0:numOfParticles * DIM]) map(from:vel_dev[0:mysize * DIM])
        {// offload begins Transfer loc[mystart:myend] foc[mystart:myend] from host to device.
          std::cout << "Thread" << dev << ":" << mystart << "~" << myend << std::endl;
	
	//Computing Velocities
	#pragma omp parallel for	
	for (int p = mystart; p < myend; p++) {
		/* zeros */
		vel_dev[(p - mystart) * DIM] = 0.0;
		vel_dev[(p - mystart) * DIM + 1] = 0.0;

		/* loop for particles  */
		for (int i = mystart; i < myend; i++) {
			double dx = loc[p * DIM] - loc[i * DIM];
			double dy = loc[p * DIM + 1] - loc[i * DIM + 1];
			double r = sqrt(dx * dx + dy * dy);
			
			//double sq;
			//double ep = EPSILON;
       			//sq = sqrt(r * r + EPSILON * EPSILON);

			double tr1 = term1(r,EPSILON) / (4.0 * M_PI);
			double tr2 = term2(r,EPSILON) / (4.0 * M_PI);
			tr2 *= foc[i * DIM] * dx + foc[i * DIM + 1] * dy;

			vel_dev[(p - mystart) * DIM] += -foc[i * DIM] * tr1 + tr2 * dx;
			vel_dev[(p - mystart) * DIM + 1] += -foc[i * DIM + 1] * tr1 + tr2 * dy;
		}
	}
	}
	for (int i = mystart ; i < myend ; i++){
	  vel[i * DIM] = vel_dev[(i - mystart)* DIM];
	  vel[i * DIM + 1] = vel_dev[(i - mystart) * DIM + 1];
	}
	delete [] vel_dev;
	}}}

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
	delete [] loc;
	delete [] vel;
	delete [] foc;
	return 0;
}

