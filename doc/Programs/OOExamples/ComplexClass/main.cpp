#include <iostream>
#include <cmath>
#include "Complex.h"
using namespace std;
int main()
{
  Complex<double> a(0.1,1.3);    // we declare a complex variable a
  Complex<double> b(3.0), c(5.0,-2.3);  // we declare  complex variables b and c
  Complex<double> d = a;         //  we declare  a new complex variable d
  Complex<double> e = d;         //  we declare  a new complex variable e
  Complex<double> f =0.0;
  d = a +c;
  e = a*c - d/b;  //   we subtract, multiply and divide two complex numbers
  f = e-d;
  cout << "Re(d)=" << d.Re() << ", Im(d)=" << d.Im() << endl;  // write out of the real and imaginary parts
  cout << "Re(d)=" << e.Re() << ", Im(d)=" << e.Im() << endl;  // write out of the real and imaginary parts
  cout << "Abs(d)=" << d.abs() << endl;  // write out absolute value
  cout << "Abs(e)=" << e.abs() << endl;  // write out absolute value
  return 0;
}
