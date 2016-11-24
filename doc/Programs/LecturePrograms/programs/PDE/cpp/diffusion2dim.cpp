/* Simple program for solving the two-dimensional diffusion 
   equation or Poisson equation using Jacobi's iterative method
   Note that this program does not contain a loop over the time 
   dependence.
*/

#include <iostream>
#include <iomanip>
#include <armadillo>
using namespace std;
using namespace arma;

int JacobiSolver(int, double, double, mat &, mat &, double);

int main(int argc, char * argv[]){
  int Npoints = 40;
  double ExactSolution;
  double dx = 1.0/(Npoints-1);
  double dt = 0.25*dx*dx;
  double tolerance = 1.0e-14;
  mat A = zeros<mat>(Npoints,Npoints);
  mat q = zeros<mat>(Npoints,Npoints);

  // setting up an additional source term
  for(int i = 0; i < Npoints; i++)
    for(int j = 0; j < Npoints; j++)
      q(i,j) = -2.0*M_PI*M_PI*sin(M_PI*dx*i)*sin(M_PI*dx*j);
    
  int itcount = JacobiSolver(Npoints,dx,dt,A,q,tolerance);
 
  // Testing against exact solution
  double sum = 0.0;
  for(int i = 0; i < Npoints; i++){
    for(int j=0;j < Npoints; j++){
      ExactSolution = -sin(M_PI*dx*i)*sin(M_PI*dx*j);
      sum += fabs((A(i,j) - ExactSolution));
    }
  }
  cout << setprecision(5) << setiosflags(ios::scientific);
  cout << "Jacobi: L2 Error is " << sum/Npoints << " in " << itcount << " iterations" << endl;
}


// Function for setting up the iterative Jacobi solver
int JacobiSolver(int N, double dx, double dt, mat &A, mat &q, double abstol)
{
  int MaxIterations = 100000;
  mat Aold = zeros<mat>(N,N);
  
  double D = dt/(dx*dx);
  
  for(int i=1;  i < N-1; i++)
    for(int j=1; j < N-1; j++)
      Aold(i,j) = 1.0;
  
  // Boundary Conditions -- all zeros
  for(int i=0; i < N; i++){
    A(0,i) = 0.0;
    A(N-1,i) = 0.0;
    A(i,0) = 0.0;
    A(i,N-1) = 0.0;
  }
  // Start the iterative solver
  for(int k = 0; k < MaxIterations; k++){
    for(int i = 1; i < N-1; i++){
      for(int j=1; j < N-1; j++){
	A(i,j) = dt*q(i,j) + Aold(i,j) +
	  D*(Aold(i+1,j) + Aold(i,j+1) - 4.0*Aold(i,j) + 
	     Aold(i-1,j) + Aold(i,j-1));
      }
    }
    double sum = 0.0;
    for(int i = 0; i < N;i++){
      for(int j = 0; j < N;j++){
	sum += (Aold(i,j)-A(i,j))*(Aold(i,j)-A(i,j));
	Aold(i,j) = A(i,j);
      }
    }
    if(sqrt (sum) <abstol){
      return k;
    }
  }
  cerr << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
  return MaxIterations;
}













