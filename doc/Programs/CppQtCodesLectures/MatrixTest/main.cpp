#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

int main()
  {
   mat A = randu<mat>(100,100);
   vec b = randu<vec>(100);

  A.print("A =");
  b.print("b=");
  // solve Ax = b
  vec x = solve(A,b);
  // print x
  x.print("x=");
  // find LU decomp of A, P is the permutation matrix
  mat L, U, P;
  lu(L,U,P,A);
  // print l
  L.print("L=");
  // print U
  U.print("U=");
  //Check that A = LU
  (A-P*L*U).print("test");
    return 0;
  }

