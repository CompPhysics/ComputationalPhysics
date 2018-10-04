#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

int main()
  {

  mat A = randu<mat>(5,5);
  mat B = A.t()*A;  // generate a symmetric matrix

  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, B);
  eigval.print("Eigenvalues=");

  eigvec.print("Print eigenvalues");
  

  return 0;
  }

