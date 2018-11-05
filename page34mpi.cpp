#include <iostream>
#include <math.h>
#include <omp.h>
#include <mpi.h>


int main(int argc, char **argv){
	//Initialize
	MPI::Init(argc, argv);

        // get myid and # of processors 
        int numproc = MPI::COMM_WORLD.Get_size();
        int myid = MPI::COMM_WORLD.Get_rank();
	
	//Let num = 10000
	const int num = 10000;
	double *aArray = new double[num];
	double *bArray = new double[num];
	double *cArray = new double[num];
	
	//Rank 0 assigns random values to those arrays.
	if (myid ==0) {
	for (int i = 0; i < num; i++){
		aArray[i] = static_cast<double>(rand()) / RAND_MAX;
		bArray[i] = static_cast<double>(rand()) / RAND_MAX;
	}}
	
	//Rank 0 sends array elements to all other ranks.
	/* Broadcast */
        MPI::COMM_WORLD.Bcast(aArray, num, MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(bArray, num, MPI::DOUBLE, 0);
	
	//Compute LaTeX /* divide loop */
        int mystart = (n / numproc) * myid;
        int myend;
        if (n % numproc > myid) {
                mystart += myid;
                myend = mystart + (n / numproc) + 1;
        } else {
                mystart += n % numproc;
                myend = mystart + (n / numproc);
        }
        std::cout << "CPU" << myid << ":" << mystart << "~" << myend << std::endl;

	 /* divide loop */
        int mystart = (num / numproc) * myid;
        int myend;
        if (num % numproc > myid) {
                mystart += myid;
                myend = mystart + (num / numproc) + 1;
        } else {
                mystart += num % numproc;
                myend = mystart + (num / numproc);
        }
        std::cout << "CPU" << myid << ":" << mystart << "~" << myend << std::endl;
	
	int mysize = myend - mystart;
	double *cArray = new double[mysize];

	for (int i = mystart; i < myend; i++){
		cArray[i-mystart] = 0.0;
		for (int j = 0; j < num; j++){
			cArray[i-mystart]+=aArray[i] * bArray[j];
		}
	}
	
	//Compute LaTeX: \left|C\right| | C | . Rank 0 gets the results. 
	double cNorm = 0.0;
	for (int i = 0; i < mysize; i++){
		cNorm += cArray[i] * cArray[i];
	}

	double tcNorm;
	MPI::COMM_WORLD.Reduce(&cNorm, &tcNorm, 1, MPI::DOUBLE, MPI::SUM, 0);
	if (myid ==0){
		tcNorm = std::sqrt(tcNorm);
		std::cout << "|C|=" << tcNorm << std::endl;
	}
	delete [] cArray;
	delete [] aArray;
	delete [] bArray;

	MPI::Finalize();	
	return 0;
}
