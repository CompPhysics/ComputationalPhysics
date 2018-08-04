#ifndef Complex_H
#define Complex_H
//   various include statements and definitions
#include <iostream>          // Standard ANSI-C++ include files
#include <new>

//  My own Complex class, now with templates but using only inline functions, not need to have functions which are longer
template<class T>
class Complex
{
private:
  T re, im; // real and imaginary part
public:
  inline Complex ();                              // Complex c;
  inline Complex (T re_a, T im_a) {re = re_a; im = im_a; }; // Definition of a complex variable;
  ~Complex () {}                        // destructor
  inline T   Re () const { return re; } //  getting the real part
  inline T   Im () const { return im; }        // T imag_part = a.Im();
  inline T   abs () const { return sqrt(re*re + im*im); }       // T m = a.abs(); // modulus
  inline friend Complex operator+ (const Complex& a, const Complex& b) { return Complex (a.re + b.re, a.im + b.im); }
  inline friend Complex  operator- (const Complex& a, const Complex& b) { return Complex (a.re - b.re, a.im - b.im); }
  inline friend Complex operator* (const Complex& a, const Complex& b) {return Complex(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);}
  inline friend Complex operator/ (const Complex& a, const Complex& b) {return Complex( (a.re*b.re +a.im*b.im)/(b.re*b.re+b.im*b.im), (-a.re*b.im + a.im*b.re)/(b.re*b.re+b.im*b.im));}
};



#endif

