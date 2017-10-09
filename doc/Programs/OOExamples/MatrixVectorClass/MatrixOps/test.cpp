#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "vectortemplates.h"

using namespace  std;


//   Main function begins here
int main(int  argc, char * argv[]){
  int dim = 2;
  Vector<double> x(dim), b(dim),x0(dim);
  
  // Set our initial guess
  x0(0) = x0(1) = 0.0;


  Vector<double> y(dim);
  y(0) = 2.;
  y(1) = -2.;

  cout << "The vector y: " << endl;
  y.Print();
  cout << endl;

 
}




