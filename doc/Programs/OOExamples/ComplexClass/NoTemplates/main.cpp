#include <iostream>
#include <cmath>
#include "mycomplex.h"
using namespace std;
int main()
{
  Complex a(0.1,1.3);    // we declare a complex variable a
  Complex b(3.0), c(5.0,-2.3);  // we declare  complex variables b and c
  Complex d = a;         //  we declare  a new complex variable d
  Complex e = d;         //  we declare  a new complex variable e
  d = a +c;
  e = a*c - d/b;  //   we subtract, multiply and divide two complex numbers
  cout << "Re(d)=" << d.Re() << ", Im(d)=" << d.Im() << endl;  // write out of the real and imaginary parts
  cout << "Re(d)=" << e.Re() << ", Im(d)=" << e.Im() << endl;  // write out of the real and imaginary parts
  cout << "Abs(d)=" << d.abs() << endl;  // write out absolute value
  cout << "Abs(e)=" << e.abs() << endl;  // write out absolute value
  return 0;
}
