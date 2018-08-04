/* matrix-vector product implemented in C++:

   let i,j=1,n:
   A(i,j) = 2 + i/j; x(j) = j/2;
   b(i) = sum_j (2 + i/j)*j/2 = sum_j j + i/2*sum_j 
        = n*(n+1)/2 + n*i/2
*/

#include <iostream>
#include <cstdio>     // C++'s view of stdio.h
#include "MatrixVector.h"

double sqr (double a)
{
  return a*a;
}

int main (int argc, char* argv[])
{
  MatVec<double> A, x, b;
  int n;
  if (argc >= 2) {
    n = atoi(argv[1]);
  } else {
    n = 5;
  }
  A.redim(n,n);  x.redim(n);  b.redim(n);

  int i,j;
  for (j=1; j<=n; j++) {
    x(j) = j/2.0;  /* or completely safe: double(j)/2.0 */
    for (i=1; i<=n; i++) {
      A(i,j) = 2.0 + double(i)/double(j);

      if (n < 10) { printf("A(%d,%d)=%g\t", i,j,A(i,j)); }
    }
    if (n < 10) { printf("  x(%d)=%g\n", j,x(j)); }
  }
  /* note: MatVec stores the array entries as in Fortran, that is
     why the j-loop is the outer loop here (different from the
     C version)
  */

  /* matrix-vector product: */
  double sum;
  for (i=1; i<=n; i++) {
    sum = 0.0;
    for (j=1; j<=n; j++) {
      sum += A(i,j)*x(j);
    }
    b(i) = sum;
  }
  /* check that b is correct: */
  sum = 0.0;
  for (i=1; i<=n; i++)
    sum += sqr(b(i) - (n*(n+1)/2.0 + n*i/2.0));
  std::cout << "square error in b = " << sum << std::endl;

  return 0;  /* success */
}

    
  
