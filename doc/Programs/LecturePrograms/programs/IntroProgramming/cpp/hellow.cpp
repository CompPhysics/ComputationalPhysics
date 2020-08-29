#include <iostream>
#include <cmath>
//#include <cstdio>
double func(double x); 
using namespace std;
int main (int argc, char* argv[])
{
  // double r, s;        // declare variables
  double r = atof(argv[1]);  // convert the text argv[1] to double
  double t = atof(argv[2]);  // convert the text argv[1] to double 
  double s = sin(r);
  double v = cos(t);
  double z = func(v);
  cout << "Hello, World! sin(" << r << ") =" << s << endl;
  cout << "Hello, World! cos(" << t << ") =" << v << endl;
  cout << "Hello, World! func(" << v << ") =" << z << endl;    
  // return 0;           // success execution of the program
}

//  This returns y*y
double func(double y)
{
  return y*y;
}  
