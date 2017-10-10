#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "vectortemplates.h"

using namespace  std;


//   Main function begins here
int main(int  argc, char * argv[]){
  int dim = 2;
  Vector<double> x(dim), b(dim);
  // Set values for x and y
  x(0) = x(1) = 10.0;
  Vector<double> y(dim);
  y(0) = 2.;
  y(1) = -2.;
  b = x+y;
  cout << "The vector b: " << endl;
  b.Print();
  cout << endl;
  cout << "The norm of  b: " << endl;
  cout<<  b.VectorNorm2() << endl;
  Vector<int> z(dim);
  z(0) = 2;
  z(1) = -2;
  cout << "The vector z: " << endl;
  z.Print();
  cout << endl;
  return 0;
}




