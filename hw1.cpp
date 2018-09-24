#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

int main()
{
	int n;
	int* res = new int[n];
	time_t time0;
	time_t time1;
	time_t time2;
	time_t time3;

	cout << "Please Enter N: " << endl;
	cin >> n;
	if (n <=1 || n>=20000){
		cout << "Invalid Number" << endl;}
	else{
	srand (time(NULL));
		
	// Dynamic 2D array
	int** A = new int*[n];
	for (int i = 0; i < n; ++i)
		A[i] = new int[n];

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			A[i][j] = rand() % 10 + 1;


	int** trans = new int*[n];
	for (int i = 0; i < n; ++i)
		trans[i] = new int[n];

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			trans[i][j] = 0;


	//for (int i = 0; i < n; ++i)
	//	for (int j = 0; j < n; ++j)
	//		cout << A[i][j] << "\n";

	// Dynamic 1D array
	int* x = new int[n];
	for (int i = 0; i < n; ++i)
		x[i] = rand() % 10 + 1;
	
	//for (int i = 0; i < n; ++i)
	//	cout << x[i] << endl;

	for (int i = 0; i < n; ++i)
		res[i] = 0;

	// Multiplication before Transpose
	time(&time0);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			res[i] += A[i][j]*x[j];

	/*cout << "Original Result: " << endl;
	for (int i = 0; i < n; ++i)
		cout << res[i] << " ";
	cout << "\n";*/

	time(&time1);

	cout << "Time Elapsed: " << time1 - time0 << " second(s) for |Ax|" << endl;
	
	//double test = time1;
	//cout << test << endl;
	//double test2 = time0;
	//cout << test << endl;
	// Multiplaction with Transpose
	
	time(&time2);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			trans[j][i]= A[i][j];

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			res[i] += trans[i][j]*x[j];

	/*cout << "Transposed Result: " << endl;
	for (int i = 0; i < n; ++i)
		cout << res[i] << " ";
	cout << "\n";
	*/

	time(&time3);
	cout << "Time Elapsed: " << time3 - time2 << " second(s) for |Atx|" << endl;

	// Deletes Array from Memory
	for (int i = 0; i < n; ++i)
		delete [] A[i];
			delete A;
	for (int i = 0; i < n; ++i)
		delete [] trans[i];
			delete trans;
	delete x;
	delete res;
	}

	return 0;
}
		

