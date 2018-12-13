//
//  hellp_omp.cpp: display a message on the screen
//
#include <iostream>
#include <omp.h>

int main () {
  int id;
  std::cout << "C++ Start" << std::endl;
#pragma omp parallel private(id)
  {
    id = omp_get_thread_num();
#pragma omp critical
    std::cout <<  "hello from " << id << std::endl;
#pragma omp barrier
    if ( id == 0 ) {
      int nthreads = omp_get_num_threads();
      std::cout <<  nthreads <<	" threads said hello!" << std::endl;
    }
  }
  std::cout << "End" << std::endl;
}
