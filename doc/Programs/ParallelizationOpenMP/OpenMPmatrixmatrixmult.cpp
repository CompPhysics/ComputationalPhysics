//  Matrix-matrix multiplication and Frobenius norm of a matrix with OpenMP
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include  <omp.h>
# include <ctime>

using namespace std; // note use of namespace
int main (int argc, char* argv[])
{
  // read in dimension of square matrix
  int n = atoi(argv[1]);
  double **A, **B, **C;
  int i, j, k;
  int thread_num;
  double wtime, Fsum, s, angle;
  cout << "  Compute matrix product C = A * B and Frobenius norm." << endl;
  omp_set_num_threads(4);
  thread_num = omp_get_max_threads ();
  cout << "  The number of processors available = " << omp_get_num_procs () << endl ;
  cout << "  The number of threads available    = " << thread_num <<  endl;
  cout << "  The matrix order n                 = " << n << endl;

  s = 1.0/sqrt( (double) n);
  wtime = omp_get_wtime ( );
  // Allocate space for the two matrices
  A = new double*[n]; B = new double*[n]; C = new double*[n];
  for (i = 0; i < n; i++){
    A[i] = new double[n];
    B[i] = new double[n];
    C[i] = new double[n];
  }
  Fsum = 0.0;
  // Define parallel region
  // private (angle, i, j, k)
# pragma omp parallel default(shared)  reduction(+:Fsum)
  {
  // Set up values for matrix A and B and zero matrix C
  # pragma omp for
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++) {
      angle = 2.0*M_PI*i*j/ (( double ) n);
      A[i][j] = s * ( sin ( angle ) + cos ( angle ) );
      B[j][i] =  A[i][j];
    }
  }
  // Then perform the matrix-matrix multiplication
  # pragma omp for
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++) {
       C[i][j] =  0.0;    
       for (k = 0; k < n; k++) {
            C[i][j] += A[i][k]*B[k][j];
       }
    }
  }
  // Compute now the Frobenius norm
  # pragma omp for
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++) {
      Fsum += C[i][j]*C[i][j];
    }
  }
  } // end parallel region and letting only one thread perform I/O
  Fsum = sqrt(Fsum);
  wtime = omp_get_wtime ( ) - wtime;
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(10) << setw(20) << "Time used  for matrix-matrix multiplication=" << wtime  << endl;
  cout << "  Frobenius norm  = " << Fsum << endl;
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

