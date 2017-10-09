#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "vector.h"

using namespace  std;


//   Main function begins here
int main(int  argc, char * argv[]){
  int dim = 2;
  Vector x(dim),xsd(dim), b(dim),x0(dim);
  
  // Set our initial guess
  x0(0) = x0(1) = 0;


  Vector y(dim);
  y(0) = 2.;
  y(1) = -2.;

  cout << "The exact solution is: " << endl;
  y.Print();
  cout << endl;
  b = y;

  cout << "The right hand side, b, of the expression Ax=b: " << endl;
  b.Print();
  cout << endl;

 
}




