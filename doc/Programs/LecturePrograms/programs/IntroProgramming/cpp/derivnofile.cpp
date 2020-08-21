#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
// Note: not using namespace for std

// Begin of main program

int main(int argc, char* argv[])
{
  if( argc <= 2 ){
    std::cout << "Bad Usage: " << argv[0] <<
      " read step size and value of x" << std::endl;
    exit(1);
  }
  // extracting number of mesh points
  double h = atof(argv[1]);
  double x = atof(argv[2]);  // reading x-value
  double Derivative = (exp(x+h)-2.*exp(x)+exp(x-h))/(h*h);
  double RelativeError = log10(fabs(Derivative-exp(x))/exp(x));
  std::cout << std::setw(15) << std::setprecision(8) << "Log10 of step size=" << log10(h) << " Log10 of relative error=" << RelativeError << std::endl;
  return 0;
}

