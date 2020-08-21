#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
// Note: not using namespace for std

// Begin of main program

int main(int argc, char* argv[])
{
  double derivative;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 2 ){
    std::cout << "Bad Usage: " << argv[0] <<
      " read number of integration points" << std::endl;
    exit(1);
  }
  // extracting number of mesh points
  int i = atoi(argv[1]);
  double x = atof(argv[2]);  // reading x-value
  double h = 1.0/((double) i); // setting up step size
  double Derivative = (exp(x+h)-2.*exp(x)+exp(x-h))/(h*h);
  double RelativeError = log10(fabs(Derivative-exp(x))/exp(x));
  std::cout << std::setw(15) << std::setprecision(8) << "relative error=" << RelativeError << std::endl;
  return 0;
}

