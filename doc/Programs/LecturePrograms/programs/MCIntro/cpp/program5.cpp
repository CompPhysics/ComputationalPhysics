// importance sampling with gaussian deviates
#include <cmath>
#include <iostream>
#include "lib.h"
using namespace std;

double gaussian_MC(double *);
double gaussian_deviate(long *);
//     Main function begins here     
int main()
{
     int n;
     double x[6], y, fx; 
     printf("Read in the number of Monte-Carlo samples\n");
     scanf("%d", &n);

     double int_mc = 0.;  double variance = 0.;
     double sum_sigma= 0. ; long idum=-1 ;  
     double length=5.; // we fix the max size of the box to L=5
     double volume=pow(acos(-1.),3.);
     double sqrt2 = 1./sqrt(2.);
//   evaluate the integral with importance sampling    
     for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
       for (int j = 0; j < 6; j++) {
	 x[j] = gaussian_deviate(&idum)*sqrt2;
       }
       fx=gaussian_MC(x); 
       int_mc += fx;
       sum_sigma += fx*fx;
     }
     int_mc = int_mc/((double) n );
     sum_sigma = sum_sigma/((double) n );
     variance=sum_sigma-int_mc*int_mc;
//   final output 
     printf("%d standard deviation= %12.5E Inum= %12.5E", n,volume*sqrt(variance/n), volume*int_mc); 
     return 0;
}  // end of main program 

// this function defines the integrand to integrate 
 
double  gaussian_MC(double *x) 
{
   double a = 0.5; 
// evaluate the different terms of the exponential
   double xy=pow((x[0]-x[3]),2)+pow((x[1]-x[4]),2)+pow((x[2]-x[5]),2);
   return exp(-a*xy);
} // end function for the integrand

// random numbers with gaussian distribution
double gaussian_deviate(long * idum)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if ( idum < 0) iset =0;
  if (iset == 0) {
    do {
      v1 = 2.*ran0(idum) -1.0;
      v2 = 2.*ran0(idum) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset =0;
    return gset;
  }
} // end function for gaussian deviates
  





