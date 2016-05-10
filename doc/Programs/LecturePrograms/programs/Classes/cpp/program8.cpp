#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "time.h"

using namespace std; // note use of namespace
int main (int argc, char* argv[])
{
  // read in dimension of square matrix
  int n = atoi(argv[1]);
  double **A, **B;
  // Allocate space for the two matrices
  A = new double*[n]; B = new double*[n];
  for (int i = 0; i < n; i++){
    A[i] = new double[n];
    B[i] = new double[n];
  }
  // Set up values for matrix A
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++) {
      A[i][j] =  cos(i*1.0)*sin(j*3.0);
    }
  }
  clock_t start, finish;
  start = clock();
  // Then compute the transpose
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++) {
      B[i][j]= A[j][i];
    }
  }

  finish = clock();
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(10) << setw(20) << "Time used  for setting up transpose of matrix=" << timeused  << endl;

  // Free up space
  for (int i = 0; i < n; i++){
    delete[] A[i];
    delete[] B[i];
  }
  delete[] A;
  delete[] B;
  return 0;
}
