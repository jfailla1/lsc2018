// Let num = 10000
//Allocate two arrays of size num, LaTeX: A_i A i  and LaTeX: B_i B i .
//Rank 0 assigns random values to those arrays.
//Rank 0 sends array elements to all other ranks.
//Compute LaTeX: C_i=\sum_j^{num}{A_iB_j} C i = ∑ j n u m A i B j  in parallel.
//Compute LaTeX: \left|C\right| | C | . Rank 0 gets the results. 

#include <iostream>
#include <math.h>

int main(int argc, char **argv){
	
	//Let num = 10000
	const int num = 10000;
	double *aArray = new double[num];
	double *bArray = new double[num];
	double *cArray = new double[num];

	for (int i = 0; i < num; i++){
		aArray[i] = static_cast<double>(rand()) / RAND_MAX;
		bArray[i] = static_cast<double>(rand()) / RAND_MAX;
		
	}
	
	for (int i = 0; i < num; i++){
		cArray[i] = 0.0;
		for (int j = 0; j < num; j++){
			cArray[i]+=aArray[i] * bArray[j];
		}
	}
	
	double cNorm = 0.0;
	for (int i = 0; i < num; i++){
		cNorm += cArray[i] * cArray[i];
	}
	cNorm = std::sqrt(cNorm);

	std::cout << "|C|=" << cNorm << std::endl;

	return 0;
}
