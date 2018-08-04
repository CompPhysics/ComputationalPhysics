#ifndef Complex_H
#define Complex_H
//   various include statements and definitions
#include <iostream>          // Standard ANSI-C++ include files
#include <new>

//  My own Complex class
class Complex
{
private:
  double re, im; // real and imaginary part
public:
  // Constructors, default, extended and copy constructor
  Complex ();                              
  Complex (double re = 0.0, double im = 0.0); 
  // copy constructor
  Complex(const Complex& c) : re(c.re), im(c.im) {}
  // destructor 
  ~Complex () {}         
  double   Re () const;        // T real_part = a.Re();
  double   Im () const;        // T imag_part = a.Im();
  double   abs () const;       // T m = a.abs(); // modulus
  // assignment operator c = a; 
  Complex& operator= (const Complex& c);  
  friend Complex operator+ (const Complex& a, const Complex& b);
  friend Complex  operator- (const Complex& a, const Complex& b);
  friend Complex operator* (const Complex& a, const Complex& b);
  friend Complex operator/ (const Complex& a, const Complex& b);
};

#endif
