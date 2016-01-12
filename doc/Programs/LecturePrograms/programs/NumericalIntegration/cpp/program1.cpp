#include <iostream>
#include "lib.h"
using namespace std;
//     Here we define various functions called by the main program

double int_function(double x);
//   Main function begins here
int main()
{
     int n;
     double a, b;
     cout << "Read in the number of integration points" << endl;
     cin >> n;
     cout << "Read in integration limits" << endl;
     cin >> a >> b;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
     double *x = new double [n];
     double *w = new double [n];
//   set up the mesh points and weights
     gauleg(a, b,x,w, n);
//   evaluate the integral with the Gauss-Legendre method
//   Note that we initialize the sum
     double int_gauss = 0.;
     for ( int i = 0;  i < n; i++){
        int_gauss+=w[i]*int_function(x[i]);
     }
//    final output
      cout << "Trapez-rule = " << trapezoidal_rule(a, b,n, &int_function)
           << endl;
      cout << "Simpson's rule = " << simpson(a, b,n, &int_function) 
           << endl;
      cout << "Gaussian quad = " << int_gauss << endl;
      delete [] x;
      delete [] w;
      return 0;
}  // end of main program
//  this function defines the function to integrate
double int_function(double x)
{
  double value = 4./(1.+x*x);
  return value;
} // end of function to evaluate
