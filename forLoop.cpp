//
//  forLoops.cpp
//
#include <iostream>
#include <omp.
// array size
static const int arraySize = 1000;

int main(int argc, const char * argv[]) {
    double arrayA[arraySize];
    double arrayB[arraySize];
    double arrayC[arraySize];
    
    // Loop 1
#pragma omp parallel for 
    for (int i = 0 ; i < arraySize ; i++){
        arrayA[i] = i;
        arrayB[i] = arraySize - i;
    }
    
    // Loop 2
#pragma omp parallel for
    for (int i = 1 ; i < arraySize - 1; i++){
        arrayC[i] = arrayA[i-1] + arrayA[i+1];
    }
  
    // Loop 3 Nope Because its possible to access the wrong value of A out of order
    for (int i = 1 ; i < arraySize - 1; i++){
        arrayA[i] = arrayA[i-1] + arrayA[i+1];
    }
    
    // Loop 4
#pragma omp parallel for
    for (int i = 1 ; i < arraySize - 1; i+=2){
        arrayB[i] = arrayB[i-1] + arrayB[i+1];
    }
    
    // Loop 5
#pragma omp parallel for
    for (int i = 0 ; i < arraySize ; i++){
        arrayC[i] = 0.0;
        for (int j = 0 ; j < arraySize ; j++){
            arrayC[i] += arrayA[i] * arrayB[j];
        }
    }
    
    // Loop 6 NOPE will not preserve the sequential behaviour if parallelized because of the break
    for (int i = 0 ; i < arraySize ; i++){
        arrayA[i] -= arrayC[i];
        if (arrayA[i] > 0.0) break;
    }
    
    // Loop 7 
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0 ; i < arraySize ; i++){
        sum += arrayA[i];
    }
    
    std::cout << "End...\n";
    return 0;
}
