#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
// Note: not using namespace for std
// output file as global variable

std::ofstream ofile;

// Begin of main program

int main(int argc, char* argv[])
{
  char *outfilename;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 3 ){
    std::cout << "Bad Usage: " << argv[0] <<
      " read also output file and number of elements on same line" << std::endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  //  opening a file for the program
  ofile.open(outfilename);
  // extracting number of mesh points
  int i = atoi(argv[2]);
  double x = atof(argv[3]);  // reading x-value
  double h = 1.0/((double) i); // setting up step size
  double Derivative = (exp(x+h)-2.*exp(x)+exp(x-h))/(h*h);
  double RelativeError = log10(fabs(Derivative-exp(x))/exp(x));
  ofile << std::setw(15) << std::setprecision(8) << "relative error=" << RelativeError << std::endl;
  ofile.close();  // close output file
  return 0;
}

