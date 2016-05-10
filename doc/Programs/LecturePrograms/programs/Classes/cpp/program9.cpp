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
  double **A, **B, **C;
  // Allocate space for the two matrices
  A = new double*[n]; B = new double*[n]; C = new double*[n];
  for (int i = 0; i < n; i++){
    A[i] = new double[n];
    B[i] = new double[n];
    C[i] = new double[n];
  }
  // Set up values for matrix A and B and zero matrix C
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++) {
      A[i][j] =  cos(i*1.0)*sin(j*3.0);
      B[i][j] =  cos(i*5.0)*sin(j*4.0);
      C[i][j] =  0.0;    
    }
  }
  clock_t start, finish;
  start = clock();
  // Then perform the matrix-matrix multiplication
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++) {
       double sum = 0.0;
       for (int k = 0; j < n; j++) {
            sum += B[i][k]*A[k][j];
       }
       C[i][j] = sum;
    }
  }
  finish = clock();
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(10) << setw(20) << "Time used  for matrix-matrix multiplication=" << timeused  << endl;

  // Free up space
  for (int i = 0; i < n; i++){
    delete[] A[i];
    delete[] B[i];
    delete[] C[i];
  }
  delete[] A;
  delete[] B;
  delete[] C;
  return 0;
}
