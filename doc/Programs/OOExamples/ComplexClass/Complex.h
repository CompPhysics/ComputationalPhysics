#ifndef Complex_H
#define Complex_H
//   various include statements and definitions
#include <iostream>          // Standard ANSI-C++ include files
#include <new>

//  My own Complex class
template <class T>
class Complex
{
private:
  T re, im; // real and imaginary part
public:
  Complex ();                              // Complex c;
  Complex (T re = 0.0, T im = 0.0); // Definition of a complex variable;
  Complex& operator= (const Complex& c); // c = a;   //  equate two complex variables
  ~Complex () {}                        // destructor
  //  T   Re () const;        // T real_part = a.Re();
  //T   Im () const;        // T imag_part = a.Im();
  //T   abs () const;       // T m = a.abs(); // modulus
  //friend Complex operator+ (const Complex& a, const Complex& b);
  //friend Complex  operator- (const Complex& a, const Complex& b);
  //friend Complex operator* (const Complex& a, const Complex& b);
  //friend Complex operator/ (const Complex& a, const Complex& b);
};

#endif
