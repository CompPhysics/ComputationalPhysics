//  Simple matrix inversion example some test for Ben

#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include    "lib.h"

using namespace std;

/* function declarations */

void inverse(double **, int);

#define   TRUE        1
#define   FALSE       0

/*
** This program sets up a simple 3x3 symmetric matrix
** and finds its determinant and inverse
*/

int main()
{
  int          i, j, k, result, n = 3;
  double       **matr, sum,  
    a[3][3]   = { {1.0, 3.0, 4.0},
		  {3.0, 4.0, 6.0},
		  {4.0, 6.0, 8.0}};

  // memory for  inverse matrix 
  matr = (double **) matrix(n, n, sizeof(double));   

  printf("\n\nInitial matrix a[][]:");      // print matrix a[][] 
  for(i = 0; i < n; i++) {
    printf("\n");
    for(j = 0; j < n; j++) {
      printf("  a[%2d][%2d] = %12.4E",i, j, a[i][j]);
    }
  }
  for(i = 0; i < n; i++) {                // transfer a[][] to matr[][] 
    for(j = 0; j < n; j++)  matr[i][j] = a[i][j];

  }
  printf("\n\n     RESULTS of LU inversion process:\n");

  inverse(matr, n);     // calculate and return inverse matrix  

  printf("\nInverse  matrix of a[][]:\n");
  for(i = 0; i < n; i++) {
    printf("\n");
    for(j = 0; j < n; j++) {
      printf("  matr[%2d][%2d] = %12.4E",i, j, matr[i][j]);
    }
  }
  result = TRUE;          // test matrix inversion 
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      sum = 0.0;
      for(k = 0; k < n; k++) sum += a[i][k] * matr[k][j];
      if((i != j) && (fabs(sum) > ZERO)) {
	printf("\nError in matrix inversion: ");
	printf(" product[%d][%d] = %10.4E\n", i, j, sum);
	result = FALSE;
      } 
      else if((i == j) && (fabs(sum - 1.0) > ZERO))  {
	printf("\nError in matrix inversion: ");
	printf(" product[%d][%d] = %10.4E\n", i, j, sum);
	result = FALSE;
      }
    }
  }   
  if(result == TRUE) {
    printf("\n\nMatrix inversion is ok!!!\n");
  }
  return 0;
} // End: function main() 

/* The function
**                inverse()
** perform a mtx inversion of the input matrix a[][] with
** dimension n. The method is described in Numerical Recipes
** sect. 2.3, page 48.
*/

void inverse(double **a, int n)
{        
  int          i,j, *indx;
  double       d, *col, **y;

  // allocate space in memory
  indx = new int[n];
  col  = new double[n];
  y    = (double **) matrix(n, n, sizeof(double)); 
   
  ludcmp(a, n, indx, &d);   // LU decompose  a[][] 

  printf("\n\nLU form of matrix of a[][]:\n");
  for(i = 0; i < n; i++) {
    printf("\n");
    for(j = 0; j < n; j++) {
      printf(" a[%2d][%2d] = %12.4E",i, j, a[i][j]);
    }
  }
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
  free_matrix((void **) y);     // release local memory 
  delete [] col;
  delete []indx;

}  // End: function inverse()
