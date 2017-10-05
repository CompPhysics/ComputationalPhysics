#include <iostream>
#include <cmath>
#include "mycomplex.h"
using namespace std;
int main()
{
  Complex<double> a(0.1,1.3);   
  Complex<double> b(3.0,5.0), c(5.0,-2.3);  
  cout << "Re(a)=" << a.Re() << ", Im(a)=" << a.Im() << endl;  
  cout << "Re(b)=" << b.Re() << ", Im(b)=" << b.Im() << endl;  
  cout << "Abs(a)=" << a.abs() << endl;  
  cout << "Abs(c)=" << c.abs() << endl;  
  cout << "Re(a*b)=" << (a*b).Re() << ", Im(a*b)=" << (a*b).Im() << endl; 
  return 0;
}
