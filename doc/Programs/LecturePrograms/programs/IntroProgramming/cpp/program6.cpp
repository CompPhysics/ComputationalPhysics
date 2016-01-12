// program to calculate the mean and standard deviation of
// a user created data set stored in array x[]
using namespace std;
#include <cmath>
#include <iostream>
int main()
{
  int      i;
  float    sum, sumsq2, xbar, sigma1, sigma2;
  // array declaration with fixed dimension
  float   x[127];
  //  initialise the data set   
  for ( i=0; i < 127 ; i++){
    x[i] = i + 100000.;
  }
  //  The variable sum is just the sum over all elements  
  //  The variable sumsq2 is the sum over x^2    
  sum=0.; 
  sumsq2=0.;
  //  Now we use the text book algorithm 
  for ( i=0; i < 127; i++){
    sum += x[i];
    sumsq2 += pow((double) x[i],2.);
  }
  //  calculate the average and sigma              
  xbar=sum/127.;
  sigma1=sqrt((sumsq2-sum*xbar)/126.);
  /*
  **  Here comes the cruder algorithm where we evaluate 
  **  separately first the average and thereafter the   
  **  sum which defines the standard deviation. The average 
  **  has already been evaluated through xbar          
  */
  sumsq2=0.;
  for ( i=0; i < 127; i++){
    sumsq2 += pow( (double) (x[i]-xbar),2.);
  }
  sigma2=sqrt(sumsq2/126.);
  cout << "xbar = " << xbar << "sigma1 = " << sigma1 << "sigma2 = " <<  sigma2;
  cout << endl;
  return 0;
}// End: function main() 
