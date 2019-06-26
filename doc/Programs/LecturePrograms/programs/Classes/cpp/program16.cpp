// This program uses its own function for allocating and freeing memory for matrices
// It can be seen as an intermediate step towards the construction of a more general 
// matrix vector class

#include <cstdlib>
#include <armadillo>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "time.h"

using namespace std; // note use of namespace
using namespace arma; // note use of namespace
int main (int argc, char* argv[])
{
  // read in dimension of square matrix
  int n = atoi(argv[1]);
  double s = 1.0/sqrt( (double) n);
  // Start timing
  clock_t start, finish;
  // Allocate space for the three matrices
  mat A(n,n), B(n,n), C(n,n);
  // Set up values for matrix A and B and zero matrix C
  for (auto i = 0; i < n; i++){
    for (auto j = 0; j < n; j++) {
      double angle = 2.0*M_PI*i*j/ (( double ) n);
      A(i,j) = s * ( sin ( angle ) + cos ( angle ) );
      B(j,i) =  A(i,j);
    }
  }
  // Then perform the matrix-matrix multiplication
  start = clock();
  C = A*B;
  finish = clock();
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(10) << setw(20) << "Time used  for matrix-matrix multiplication=" << timeused  << endl;

  return 0;
}















