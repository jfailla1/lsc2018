#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv){
	//Initialize
	MPI::Init(argc,argv);
	
	int numproc = MPI::COMM_WORLD.Get_size();
	int myid = MPI::COMM_WORLD.Get_rank(); 

	if (argc < 2) {
		std::cout << argv[0] << " number\n";
		return 0;
	}
	srand (time(NULL));
	int num;	
	num =atoi(argv[1]);
	int count = 0;
	for (int i = 0; i < num; i++){
		double x = static_cast<double>(rand()) / RAND_MAX;
		double y = static_cast<double>(rand()) / RAND_MAX;
		if (x * x + y * y <= 1.0) {
			count++;
		}
	}

	int total_count;
	MPI::COMM_WORLD.Reduce(&count, &total_count, 1, MPI::INT, MPI::Sum, 0);
	if (myid == 0){
    	double res = 4 * static_cast<double>(count) / numproc * num;
	std::cout << "Result=" << res << std::endl;
	}
	
	MPI::Finalize();
	return 0;
}
