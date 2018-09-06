// This program uses its own function for allocating and freeing memory for matrices
// It can be seen as an intermediate step towards the construction of a more general 
// matrix vector class

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "time.h"
#include "lib.h"

//  Declaring functions to allocate and free space for a matrix

double ** CreateMatrix(int m, int n);
void DestroyMatrix(double ** mat, int m, int n); 

int ** ICreateMatrix(int m, int n);
void IDestroyMatrix(int ** mat, int m, int n);

// function to compute the inverse of a matrix
void inverse(double **, int);

using namespace std; // note use of namespace
int main (int argc, char* argv[])
{
  // read in dimension of square matrix
  int n = atoi(argv[1]);
  double **A;
  // Start timing
  clock_t start, finish;
  start = clock();
  // Allocate space for the a matrix
  A = CreateMatrix(n, n);
  // Set up values for matrix A and B and zero matrix C
  for (int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      A[i][j] = 1.0+(double) i+j;
      //  symmetric matrix here 
      A[j][i] = A[i][j];
    }
  }
  cout << " Matrix a[][] before LU decomposition: " << endl;
  for(int i = 0; i < n; i++) {
    cout << endl;
    for(int j = 0; j < n; j++) {
       printf(" a[%2d][%2d] = %12.4E",i, j, A[i][j]);
    }
  }
  cout << endl;
  inverse(A,n);
  finish = clock();
  cout << " Matrix a[][] after inversion: " << endl;
  for(int i = 0; i < n; i++) {
    cout << endl;
    for(int j = 0; j < n; j++) {
      printf(" a[%2d][%2d] = %12.4E",i, j, A[i][j]);
    }
  }
  cout << endl;
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(10) << setw(20) << "Time used  for matrix inversion=" << timeused  << endl;
  // Free up space
  DestroyMatrix(A, n, n); 
  return 0;
}


//  Allocating space for a double type matrix
double ** CreateMatrix(int m, int n){
  double ** mat;
  mat = new double*[m];
  for(int i=0;i<m;i++){
    mat[i] = new double[n];
    for(int j=0;j<m;j++)
      mat[i][j] = 0.0;
  }
  return mat;
}

//  Allocating space for an integer type  matrix
int ** ICreateMatrix(int m, int n){
  int ** mat;
  mat = new int*[m];
  for(int i=0;i<m;i++){
    mat[i] = new int[n];
    for(int j=0;j<m;j++)
      mat[i][j] = 0;
  }
  return mat;
}

// Freeing space for a double type matrix 
void DestroyMatrix(double ** mat, int m, int n){
  for(int i=0;i<m;i++)
    delete[] mat[i];
  delete[] mat;
}

// Freeing space for an integer type matrix 
void IDestroyMatrix(int ** mat, int m, int n){
  for(int i=0;i<m;i++)
    delete[] mat[i];
  delete[] mat;
}

// Function to compute the inverse of a matrix
void inverse(double **a, int n)
{
  int          i,j, *indx;
  double       d, *col, **y;
  // allocate space in memory
  indx = new int[n];
  col  = new double[n];
  y    = CreateMatrix(n, n);
  ludcmp(a, n, indx, &d);   // LU decompose  a[][]
  printf("\n\nLU form of matrix of a[][]:\n");
  for(i = 0; i < n; i++) {
    printf("\n");
    for(j = 0; j < n; j++) {
      printf(" a[%2d][%2d] = %12.4E",i, j, a[i][j]);
    }
  }
  cout << endl;
  // find inverse of a[][] by columns
  for(j = 0; j < n; j++) {
    // initialize right-side of linear equations
    for(i = 0; i < n; i++) col[i] = 0.0;
    col[j] = 1.0;
    lubksb(a, n, indx, col);
    // save result in y[][]
    for(i = 0; i < n; i++) y[i][j] = col[i];
  }   //j-loop over columns
  // return the inverse matrix in a[][]
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) a[i][j] = y[i][j];
  }
  DestroyMatrix(y, n, n); 
  delete [] col;
  delete []indx;
}  // End: function inverse()
