#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "vectormatrixclass.h"

using namespace  std;

  Vector ConjugateGradient(Matrix A, Vector b, Vector x0){
  int dim = x0.Dimension();
  const double tolerance = 1.0e-14;
  Vector x(dim),r(dim),v(dim),z(dim);
  double c,t,d;

  x = x0;
  r = b - A*x;
  v = r;
  c = dot(r,r);
  for(int i=0;i<dim;i++){
    if(sqrt(dot(v,v))<tolerance){
      cerr << "An error has occurred in ConjugateGradient: execution of function terminated" << endl;
      break;
    }
    z = A*v;
    t = c/dot(v,z);
    x = x + t*v;
    r = r - t*z;
    d = dot(r,r);
    if(sqrt(d) < tolerance)
      break;
    v = r + (d/c)*v;
    c = d;
  }
  return x;
} 


Vector SteepestDescent(Matrix A, Vector b, Vector x0){
  int IterMax, i;
  int dim = x0.Dimension();
  const double tolerance = 1.0e-14;
  Vector x(dim),f(dim),z(dim);
  double c,alpha,d;
  IterMax = 30;
  x = x0;
  f = A*x-b;
  i = 0;
  while (i <= IterMax || sqrt(dot(f,f)) < tolerance ){
    if(sqrt(dot(f,f))<tolerance){
       cerr << "An error has occurred: execution of function terminated" << endl;
       break;
    }
    z = A*f;
    c = dot(f,f);
    alpha = c/dot(f,z);
    x = x - alpha*f;
    f =  A*x-b;
    if(sqrt(dot(f,f)) < tolerance) break;
    i++;
  }
  return x;
} 



//   Main function begins here
int main(int  argc, char * argv[]){
  int dim = 2;
  Vector x(dim),xsd(dim), b(dim),x0(dim);
  Matrix A(dim,dim);
  
  // Set our initial guess
  x0(0) = x0(1) = 0;
  // Set the matrix  
  A(0,0) =  3;    A(1,0) =  2;   A(0,1) =  2;   A(1,1) =  6; 
  
  cout << "The Matrix A that we are using: " << endl;
  A.Print();
  cout << endl;

  Vector y(dim);
  y(0) = 2.;
  y(1) = -2.;

  cout << "The exact solution is: " << endl;
  y.Print();
  cout << endl;
  b = A*y;

  cout << "The right hand side, b, of the expression Ax=b: " << endl;
  b.Print();
  cout << endl;

  x = ConjugateGradient(A,b,x0);

  xsd = SteepestDescent(A,b,x0);
  
  cout << "The approximate solution using Conjugate Gradient is: " << endl;
  x.Print();
  cout << endl;

  cout << "The approximate solution using Steepest Descent is: " << endl;
  xsd.Print();
  cout << endl;


 
}




