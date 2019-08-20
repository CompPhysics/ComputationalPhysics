// This program uses its own function for allocating and freeing memory for matrices
// It can be seen as an intermediate step towards the construction of a more general 
// matrix vector class

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "time.h"


using namespace std; // note use of namespace
int main (int argc, char* argv[])
{
  // read in dimension of square matrix
  int n = atoi(argv[1]);
  double s = 1.0/sqrt( (double) n);
  // Start timing
  clock_t start, finish;
  // Allocate space for the three matrices as one-dimensional vectors
  vector <double> A(n*n);
  vector <double> B(n*n);
  vector <double> C(n*n);
  // Set up values for matrix A
  for (auto i = 0; i < n; i++){
    for (auto j = 0; j < n; j++) {
      double angle = 2.0*M_PI*i*j/ (( double ) n);
      A[i*n+j] = s * ( sin ( angle ) + cos ( angle ) );
      B[i*n+j] = A[i*n+j];  
    }
  }
  // Then perform the matrix-matrix multiplication
  start = clock();
  for(auto i=0;i<n;i++){
    for(auto j=0;j<n;j++){
      C[i*n+j] = 0.0;
      for(auto k=0;k<n;k++){
	C[i*n+j]+=A[k*n+j]*B[k+n*i];
      }
    }
  }
  finish = clock();
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(10) << setw(20) << "Time used  for matrix-matrix multiplication using vectors=" << timeused  << endl;
  return 0;
}









