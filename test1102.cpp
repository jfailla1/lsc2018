#include <iostream>
#include <iomanip>
#include <omp.h>
#include <mpi.h>

double myFunction(double x){
	return 4/ (1 + x * x);
}

int main(int argc, char **argv){
	//Initialize
	MPI::Init(argc, argv);

	if (argc < 2){
		std::cout << argv[0] << "[N]\n";
		MPI::Finalize();
		return 0;
	}
	int n = atoi(argv[1]);
	
	// get myid and # of processors
	int numproc = MPI::COMM_WORLD.Get_size();
	int myid = MPI::COMM_WORLD.Get_rank();

	int mystart = (n/numproc)*myid;
	int myend;
	
	if (n % numproc > myid) {
		mystart += myid;
		myend = mystart + (n /numproc) + 1;
	} else {
		mystart += n % numproc;
		myend = mystart + (n / numproc);
	}
	
	double p = 0;
	//#pragma omp parallel
	//{
	//#pragma omp parallel for reduction(+:p)
	for (int i = 0; i < n ; i++){
		p+= myFunction((i + 0.5) / n) / n;
	}
	double psum;
	//}
	MPI::COMM_WORLD.Reduce(&p, &psum, 1 , MPI::DOUBLE, MPI::SUM, 0); 
	if (myid == 0) {
	std::cout << "P = " << std::setprecision(15) << p << std::endl;
	}
	MPI::Finalize();
	return 0;
}
