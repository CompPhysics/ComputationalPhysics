//  Diagonalizing tridiagonal Toeplitz matrix  with Lapack functions
//  Compile as c++ -O3 -o Tridiag.x TridiagToeplitz.cpp -larmadillo -llapack -lblas
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "Tridiag.h" 
using namespace  std;
using namespace  arma;

vec GetEigenvalues(int Dim)
 {
   int       i, j;
   double    RMin, RMax, Step, DiagConst, NondiagConst; 
   RMin = 0.0; RMax = 1.0;
   mat Hamiltonian = zeros<mat>(Dim,Dim);
   // Integration step length
   Step    = RMax/ Dim;
   DiagConst = 2.0 / (Step*Step);
   NondiagConst =  -1.0 / (Step*Step);
   // Setting up tridiagonal matrix and diagonalization using Armadillo
   Hamiltonian(0,0) = DiagConst;
   Hamiltonian(0,1) = NondiagConst;
   for(i = 1; i < Dim-1; i++) {
     Hamiltonian(i,i-1)    = NondiagConst;
     Hamiltonian(i,i)    = DiagConst;
     Hamiltonian(i,i+1)    = NondiagConst;
   }
   Hamiltonian(Dim-1,Dim-2) = NondiagConst;
   Hamiltonian(Dim-1,Dim-1) = DiagConst;
   // diagonalize and obtain eigenvalues
   vec Eigval(Dim);
   eig_sym(Eigval, Hamiltonian);
   return Eigval;
}  //  end of eigenvalue function


