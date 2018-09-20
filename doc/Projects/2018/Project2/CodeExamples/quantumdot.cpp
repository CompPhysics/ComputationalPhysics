//  Diagonalizing tridiagonal Toeplitz matrix  with Lapack functions
//  Compile as c++ -O3 -o Tridiag.x TridiagToeplitz.cpp -larmadillo -llapack -lblas
//  Here we add the harmonic oscillator to the tridiagonal Toeplitz matrix and 
//  use Richardson's extrapolation scheme to obtain results to high precision
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace  std;
using namespace  arma;
double potential(double);
void PolynomialInterpolation(vec, vec, int, double, double *, double *);
// Begin of main program   

int main(int argc, char* argv[])
{
  int       i, j, StartDim, Dim;
  double    RMin, RMax, Step, DiagConst, NondiagConst; 
  double error, FinalE;
  RMin = 0.0; RMax = 10.0; StartDim =400;
  vec StepSizes = zeros<vec>(5);  
  mat EigvaluesInterpol = zeros<mat>(5,20);  
  for (int interpolation = 0; interpolation < 5; interpolation++){
    if (interpolation == 1) {
      Dim = StartDim;
    }
    else{
      Dim = (interpolation-1)*2*StartDim;
    } 
    mat Hamiltonian = zeros<mat>(Dim,Dim);
    // Integration step length
    Step    = RMax/Dim;
    DiagConst = 2.0 / (Step*Step);
    NondiagConst =  -1.0 / (Step*Step);
    // local memory for r and the potential w[r] 
    double r;
    vec Potential = zeros<vec>(Dim);
    for(i = 0; i < Dim; i++) {
      r = RMin + (i+1) * Step;
      Potential[i] = potential(r);
    }
    StepSizes(interpolation) = Step;
    // Setting up tridiagonal matrix and diagonalization using Armadillo
    Hamiltonian(0,0) = DiagConst+Potential(0);
    Hamiltonian(0,1) = NondiagConst;
    for(i = 1; i < Dim-1; i++) {
      Hamiltonian(i,i-1)    = NondiagConst;
      Hamiltonian(i,i)    = DiagConst+Potential(i);
      Hamiltonian(i,i+1)    = NondiagConst;
    }
    Hamiltonian(Dim-1,Dim-2) = NondiagConst;
    Hamiltonian(Dim-1,Dim-1) = DiagConst+Potential(Dim-1);
    // diagonalize and obtain eigenvalues
    vec Eigval(Dim);
    //    eig_sym(Eigval, Hamiltonian);
    for (int l = 0; l < 20; l++){
      EigvaluesInterpol(interpolation,l) = Eigval(l);
      vec EigSend(interpolation);
           for (int k = 0; k < interpolation; k++) EigSend(k) = EigvaluesInterpol(k,l);
      if (interpolation > 0 ){
       	  PolynomialInterpolation(StepSizes, EigSend, interpolation, 0.0, &FinalE, &error);
	  cout << setiosflags(ios::showpoint | ios::uppercase);
	cout <<"Eigenvalue = " << setw(15) << l << endl;  
	cout << setw(15) << setprecision(8) << FinalE << endl;
      }
    }
  }
  return 0;
 }  //  end of main function


  /*
    The function potential()
    calculates and return the value of the 
    potential for a given argument x.
    The potential here is for the hydrogen atom
  */        

  double potential(double x)
  {
    return x*x;

  } // End: function potential()  


  void PolynomialInterpolation(vec xa, vec ya, int n, double x, double *y, double *dy)
  {
    int i,m,ns=1;
    double den,dif,dift,ho,hp,w;

    dif=fabs(x-xa(0));
    vec c = zeros<vec>(n);
    vec d = zeros<vec>(n);
    for (i=0;i < n;i++) {
      if ( (dift=fabs(x-xa(i))) < dif) {
	ns=i;
	dif=dift;
      }
      c(i)=ya(i);
      d(i)=ya(i);
    }
    *y=ya(ns--);
    for (m=0; m< n;m++) {
      for (i = 0;i < n-m;i++) {
	ho=xa(i)-x;
	hp=xa(i+m)-x;
	w=c(i+1)-d(i);
	//			if ( (den=ho-hp) == 0.0) nrerror("Error ");
	den=w/den;
	d(i)=hp*den;
	c(i)=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c(ns+1) : d(ns--)));
    }
  }

