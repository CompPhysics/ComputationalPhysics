//  Diagonalizing tridiagonal Toeplitz matrix  with Lapack functions
//  Compile as c++ -O3 -o Tridiag.x TridiagToeplitz.cpp -larmadillo -llapack -lblas
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace  std;
using namespace  arma;

// Begin of main program   

int main(int argc, char* argv[])
{
  int       i, j, Dim;
  double    RMin, RMax, Step, DiagConst, NondiagConst; 
  RMin = 0.0; RMax = 1.0; Dim =20;  
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
  double pi = acos(-1.0);
  cout << "RESULTS:" << endl;
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout <<"Number of Eigenvalues = " << setw(15) << Dim << endl;  
  cout << "Exact versus numerical eigenvalues:" << endl;
  for(int i = 0; i < Dim; i++) {
    double Exact = DiagConst+2*NondiagConst*cos((i+1)*pi/(Dim+1));
    cout << setw(15) << setprecision(8) << fabs(Eigval[i]-Exact) << endl;
  }
  return 0;
}  //  end of main function

