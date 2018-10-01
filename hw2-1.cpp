//
//  Triangle Matrix
//
//  Created by Hideki Fujioka on 9/11/18.
//
#include <iostream>
#include <cmath>
#include <algorithm>

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << argv[0] << " [size]\n";
        return 0;
    }
    
    int size = std::atoi(argv[1]);
    std::cout << "Size=" << size << std::endl;
    
    // allocation
    double *xvec = new double[size];
    double *bvec = new double[size];
    // Make 2D array
    double *amatBlock = new double[size * size];
    double **amat = new double*[size];
    for (int i = 0 ; i < size ; i++){
        amat[i] = amatBlock + size * i;
    }
    // set zeros
    std::fill(amatBlock,amatBlock+size * size,0.0);
    
    // setup matrix A lower elements
    for (int i = 0; i < size; i++) {
        double diag = 1.0;
        for (int j = 0; j < i; j++) {
            double elm = (rand() % 100) * 0.01;
            amat[i][j] = elm;
            diag += elm;
        }
        // diagonal element = -sum(off-diagonal elements)
        amat[i][i] = -diag;
    }
    
    // setup vector x elements
    for (int i = 0; i < size; i++) {
        xvec[i] = i;
    }

    // compute rhs
    std::fill(bvec,bvec+size,0.0);
    for (int i = 0; i < size ; i++) {
        for (int j = 0; j < size ; j++) {
            bvec[i] += amat[i][j] * xvec[j];
        }
    }
    
    // compute sol = A^-1 b
    double *sol = new double[size];
/*

WRITE YOUR CODE HERE

*/
    
   for (int i = 0; i < size; i++) {
	  double sumAs = 0.0;
	   for (int j = 0; j < i; j++) {
		sumAs += amat[i][j] *  sol[j];
	  }
	  sol[i] = (bvec[i] - sumAs) / amat[i][i];
    }	  
   
    // check solution
    double r = 0.0;
    for (int i = 0; i < size; i++){
        double dr = sol[i] - xvec[i];
        r += dr * dr;
    }
    r = std::sqrt(r);
    std::cout << "|x - x0|=" << r << std::endl;
    
    return 0;
}
