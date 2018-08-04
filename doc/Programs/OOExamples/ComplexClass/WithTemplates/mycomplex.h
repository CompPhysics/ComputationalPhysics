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
  Complex ();                              // Complex c;
  Complex (T re_a, T im_a) {re = re_a; im = im_a; }; // Definition of a complex variable;
  // copy constructor
  Complex(const Complex& c) : re(c.re), im(c.im) {}
  // destructor
  ~Complex () {}                        // destructor
  Complex& operator= (const Complex& c)
    {
      re = c.re;
      im = c.im;
      return *this;
    }
  T   Re () const { return re; } //  getting the real part
  T   Im () const { return im; }        // T imag_part = a.Im();
  T   abs () const { return sqrt(re*re + im*im); }       // T m = a.abs(); // modulus
  friend Complex operator+ (const Complex& a, const Complex& b) { return Complex (a.re + b.re, a.im + b.im); }
  friend Complex  operator- (const Complex& a, const Complex& b) { return Complex (a.re - b.re, a.im - b.im); }
  friend Complex operator* (const Complex& a, const Complex& b) {return Complex(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);}
  friend Complex operator/ (const Complex& a, const Complex& b) {return Complex( (a.re*b.re +a.im*b.im)/(b.re*b.re+b.im*b.im), (-a.re*b.im + a.im*b.re)/(b.re*b.re+b.im*b.im));}
};



#endif

