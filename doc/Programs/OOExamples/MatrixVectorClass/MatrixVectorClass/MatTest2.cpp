/* matrix-vector product implemented in C++:

   let i,j=1,n:
   A(i,j) = 2 + i/j; x(j) = j/2;
   b(i) = sum_j (2 + i/j)*j/2 = sum_j j + i/2*sum_j 
        = n*(n+1)/2 + n*i/2
*/

#include <iostream>
#include <cstdio>     // C++'s view of stdio.h
#include "MatrixVector.h"

/* the signature of a function must be known before we can call it: */
double sqr (double a);
void init (MatVec<double>& A, MatVec<double>& x);
void mvproduct (const MatVec<double>& A, 
		const MatVec<double>& x, 
		MatVec<double>& b);

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

  init(A, x);
  mvproduct(A, x, b);

  /* check that b is correct: */
  double sum = 0.0; int i;
  for (i=1; i<=n; i++)
    sum += sqr(b(i) - (n*(n+1)/2.0 + n*i/2.0));
  std::cout << "square error in b = " << sum << std::endl;

  return 0;  /* success */
}


/* here we implement the various functions: */

inline double sqr (double a)
{
  return a*a;
}
    
void init (MatVec<double>& A, MatVec<double>& x)
{
  const int n = x.size();
  int i,j;
  for (j=1; j<=n; j++) {
    x(j) = j/2.0;  /* or completely safe: double(j)/2.0 */
    for (i=1; i<=n; i++) {
      A(i,j) = 2.0 + double(i)/double(j);

      if (n < 10) { printf("A(%d,%d)=%g\t", i,j,A(i,j)); }
    }
    if (n < 10) { printf("  x(%d)=%g\n", j,x(j)); }
  }
}

void mvproduct (const MatVec<double>& A, const MatVec<double>& x, 
		MatVec<double>& b)
{
  const int n = x.size();
  double sum; int i,j;
  for (i=1; i<=n; i++) {
    sum = 0.0;
    for (j=1; j<=n; j++) {
      sum += A(i,j)*x(j);
    }
    b(i) = sum;
  }
}
