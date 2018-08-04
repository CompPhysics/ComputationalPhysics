//  This function computes the autocorrelation function for 
//  the Mersenne random number generator with a uniform distribution
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include <cmath>
using namespace  std;
using namespace arma;
// output file
ofstream ofile;

//     Main function begins here     
int main(int argc, char* argv[])
{
  int MonteCarloCycles;
  string filename;
  if (argc > 1) {
    filename=argv[1];
    MonteCarloCycles = atoi(argv[2]);
    string fileout = filename;
    string argument = to_string(MonteCarloCycles);
    fileout.append(argument);
    ofile.open(fileout);
  }

  // Compute the variance and the mean value of the uniform distribution
  // Compute also the specific values x for each cycle in order to be able to
  // compute the covariance and the correlation function  

  vec X  = zeros<vec>(MonteCarloCycles);
  double MCint = 0.;      double MCintsqr2=0.;
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  for (int i = 0;  i < MonteCarloCycles; i++){
    double x =   RandomNumberGenerator(gen); 
    X(i) = x;
    MCint += x;
    MCintsqr2 += x*x;
  }
  double Mean = MCint/((double) MonteCarloCycles );
  MCintsqr2 = MCintsqr2/((double) MonteCarloCycles );
  double STDev = sqrt(MCintsqr2-Mean*Mean);
  double Variance = MCintsqr2-Mean*Mean;
  //   Write mean value and variance
  cout << " Sample variance= " << Variance  << " Mean value = " << Mean << endl;
  // Now we compute the autocorrelation function
  vec autocorrelation = zeros<vec>(MonteCarloCycles);
  for (int j = 0; j < MonteCarloCycles; j++){
    double sum = 0.0;
    for (int k = 0; k < (MonteCarloCycles-j); k++){
      sum  += (X(k)-Mean)*(X(k+j)-Mean); 
    }
    autocorrelation(j) = sum/Variance/((double) MonteCarloCycles );
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << j;
    ofile << setw(15) << setprecision(8) << autocorrelation(j) << endl;
  }
  // Now compute the exact covariance using the autocorrelation function
  double Covariance = 0.0;
  for (int j = 0; j < MonteCarloCycles; j++){
    Covariance  += autocorrelation(j);
  }
  Covariance *=  2.0/((double) MonteCarloCycles);
  // Compute now the total variance, including the covariance, and obtain the standard deviation
  double TotalVariance = (Variance/((double) MonteCarloCycles ))+Covariance;
  cout << "Covariance =" << Covariance << "Totalvariance= " << TotalVariance << "Sample Variance/n= " << (Variance/((double) MonteCarloCycles )) << endl;
  cout << " STD from sample variance= " << sqrt(Variance/((double) MonteCarloCycles )) << " STD with covariance = " << sqrt(TotalVariance) << endl;

  // Finally we start the blocking procedure and compute only the sample variance and the mean (to check that it ispreserved)
  int n = MonteCarloCycles;
  for (int i = 1; i <= 14; i++) {
     n /= 2; 
     cout << n << endl;
     vec Xp = zeros<vec>(n);
     Xp(0) = 0.5*(X(0)+X(1));
     for (int j = 1; j < n; j++){
       Xp(j) = 0.5*(X(2*j-1)+X(2*j));     
     }
     double MCint = 0.0; double MCintsqr2 =0.0;
     for (int i = 0;  i < n; i++){
       MCint += Xp(i);
       MCintsqr2 += Xp(i)*Xp(i);
     }
     double Mean = MCint/((double) n );
     MCintsqr2 = MCintsqr2/((double) n );
     double STDev = sqrt(MCintsqr2-Mean*Mean);
     cout << " " << Mean << " " << STDev/sqrt((double) n) << endl;
     vec X =  zeros<vec>(n); X = Xp;
  }

  ofile.close();  // close output file
  return 0;
}  // end of main program 

