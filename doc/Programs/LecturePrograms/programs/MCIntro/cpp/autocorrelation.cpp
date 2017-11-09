//  This function computes the autocorrelation function for 
//  the standard c++ random number generator

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <string>
using namespace std;
// output file as global variable
ofstream ofile;  

//     Main function begins here     
int main(int argc, char* argv[])
{
  string filename;
  filename=argv[1];
  int n = atoi(argv[2]);  // number of Monte Carlo samples
  string fileout = filename;
  string argument = to_string(n);
  fileout.append(argument);
  ofile.open(fileout);

  double MCint = 0.;      double MCintsqr2=0.;
  double invers_period = 1./RAND_MAX; // initialise the random number generator
  double *X;  
  X = new double[n];
  for (int i = 0;  i < n; i++){
    double x = double(rand())*invers_period; 
    X[i] = x;
    MCint += x;
    MCintsqr2 += x*x;
  }
  double Mean = MCint/((double) n );
  MCintsqr2 = MCintsqr2/((double) n );
  double STDev = sqrt(MCintsqr2-Mean*Mean);
  double Variance = MCintsqr2-Mean*Mean;
  //   Write mean value and standard deviation 
  cout << " Standard deviation= " << STDev << " Integral = " << Mean << endl;

  // Now we compute the autocorrelation function
  double *autocor;  autocor = new double[n];
  for (int j = 0; j < n; j++){
    double sum = 0.0;
    for (int k = 0; k < (n-j); k++){
      sum  += (X[k]-Mean)*(X[k+j]-Mean); 
    }
    autocor[j] = sum/Variance/((double) n );
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << j;
    ofile << setw(15) << setprecision(8) << autocor[j] << endl;
  }
  ofile.close();  // close output file
  return 0;
}  // end of main program 
