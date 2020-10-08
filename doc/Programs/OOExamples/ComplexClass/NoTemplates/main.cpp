#include <iostream>
#include <cmath>
#include "mycomplex.h"
using namespace std;
int main()
{
  // we declare a complex variable a
  Complex a(0.1,1.3);    
  // we declare  complex variables b and c
  Complex b(3.0), c(5.0,-2.3);  
  // using the copy constructor to define a new complex variable z=c
  //  Complex z(c);   // Could use C++11 as z{c}
  Complex z{c};   // Could use C++11 as z{c}
  // C++11 way of declaring compile with c++ -std=c++11
  Complex g{3,4};
  cout << g.Re() << " " << g.Im() << endl;
  //  we declare  a new complex variable d using the assignment operator
  Complex d = z;         
  //  we declare  a new complex variable e using the assignment operator
  Complex e = d;         
  d = a +c;
  e = a*c - d/b;  
  // write out of the real and imaginary parts
  cout << "Re(d)=" << d.Re() << ", Im(d)=" << d.Im() << endl;  
  cout << "Re(d)=" << e.Re() << ", Im(d)=" << e.Im() << endl;
  // write out absolute value
  cout << "Abs(d)=" << d.abs() << endl;  
  cout << "Abs(e)=" << e.abs() << endl;
  return 0;
}
