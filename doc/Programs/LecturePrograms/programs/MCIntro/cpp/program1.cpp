#include <cmath>
#include <iostream>
#include "lib.h"
using namespace std;

//     Here we define various functions called by the main program  
//     this function defines the function to integrate  
double func(double x);
//     Main function begins here     
int main()
{
     int i, n;
     long idum;
     double crude_mc, x, sum_sigma, fx, variance, exact; 
     printf("Read in the number of Monte-Carlo samples\n");
     scanf("%d", &n);
     crude_mc = sum_sigma=0. ; idum=-1 ;  exact=acos(-1.);
//   evaluate the integral with a crude Monte-Carlo method    
     for ( i = 1;  i <= n; i++){
           x=ran0(&idum);
           fx=func(x);
           crude_mc += fx;
           sum_sigma += fx*fx;
      }
      crude_mc = crude_mc/((double) n );
      sum_sigma = sum_sigma/((double) n );
      variance=sum_sigma-crude_mc*crude_mc;
//    final output 
      printf("%d variance= %12.5E Inum= %12.5E pi= %12.5E\n", n, variance, crude_mc, exact); 
      return 0;
}  // end of main program 

// this function defines the function to integrate 

double func(double x)
{
  double value;
  value = 4/(1.+x*x);
  return value;
} // end of function to evaluate 
