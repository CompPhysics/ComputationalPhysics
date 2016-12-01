/* Simple program for solving the two-dimensional diffusion 
   equation or Poisson equation using Jacobi's iterative method
   Note that this program does not contain a loop over the time 
   dependence. It uses OpenMP to parallelize the spatial part.
*/

#include <iostream>
#include <iomanip>
#include <armadillo>
#include <omp.h>
using namespace std;
using namespace arma;

int JacobiSolver(int, double, double, mat &, mat &, double);

int main(int argc, char * argv[]){
  int Npoints = 100;
  double ExactSolution;
  double dx = 1.0/(Npoints-1);
  double dt = 0.25*dx*dx;
  double tolerance = 1.0e-8;
  mat A = zeros<mat>(Npoints,Npoints);
  mat q = zeros<mat>(Npoints,Npoints);

  int thread_num;
  omp_set_num_threads(4);
  thread_num = omp_get_max_threads ();
  cout << "  The number of processors available = " << omp_get_num_procs () << endl ;
  cout << "  The number of threads available    = " << thread_num <<  endl;
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
  cout << "Jacobi error is " << sum/Npoints << " in " << itcount << " iterations" << endl;
}


// Function for setting up the iterative Jacobi solver
int JacobiSolver(int N, double dx, double dt, mat &A, mat &q, double abstol)
{
  int MaxIterations = 100000;

  
  double D = dt/(dx*dx);
  // initial guess
  mat Aold = randu<mat>(N,N);  
  // Boundary conditions, all zeros
  for(int i=0; i < N; i++){
    A(0,i) = 0.0;
    A(N-1,i) = 0.0;
    A(i,0) = 0.0;
    A(i,N-1) = 0.0;
  }
  double sum = 1.0;
  int k = 0;
  // Start the iterative solver
  while (k < MaxIterations && sum > abstol){
  int i, j;
  sum = 0.0;
  // Define parallel region
# pragma omp parallel default(shared) private (i, j) reduction(+:sum)
  {
     # pragma omp for
    for(i = 1; i < N-1; i++){
      for(j = 1; j < N-1; j++){
	A(i,j) = dt*q(i,j) + Aold(i,j) +
	  D*(Aold(i+1,j) + Aold(i,j+1) - 4.0*Aold(i,j) + 
	     Aold(i-1,j) + Aold(i,j-1));
      }
    }
    for(i = 0; i < N;i++){
      for(j = 0; j < N;j++){
	sum += fabs(Aold(i,j)-A(i,j));
	Aold(i,j) = A(i,j);
      }
    }
    sum /= (N*N); 
  } //end parallel region
  k++;
  } //end while loop
  return k;
}













