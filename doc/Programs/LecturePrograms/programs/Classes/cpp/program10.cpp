#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "time.h"

using namespace std; // note use of namespace
int main (int argc, char* argv[])
{
  // read in dimension of square matrix
  int n = atoi(argv[1]);
  double s = 1.0/sqrt( (double) n);
  double *a, *b, *c;
  // Start timing
  clock_t start, finish;
  start = clock();
// Allocate space for the vectors to be used
    a = new double [n]; b = new double [n]; c = new double [n];
  // Define parallel region
  // Set up values for vectors  a and b
  for (int i = 0; i < n; i++){
    double angle = 2.0*M_PI*i/ (( double ) n);
    a[i] = s*(sin(angle) + cos(angle));
    b[i] =  s*sin(2.0*angle);
    c[i] = 0.0;
  }
  // Then perform the vector addition
  for (int i = 0; i < n; i++){
    c[i] += a[i]+b[i];
  }
  // Compute now the norm-2
  double Norm2 = 0.0;
  for (int i = 0; i < n; i++){
    Norm2  += c[i]*c[i];
  }
  finish = clock();
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(10) << setw(20) << "Time used  for norm computation=" << timeused  << endl;
  cout << "  Norm-2  = " << Norm2 << endl;
  // Free up space
  delete[] a;
  delete[] b;
  delete[] c;
  return 0;
}



