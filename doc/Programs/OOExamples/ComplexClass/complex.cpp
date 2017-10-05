#include <cmath>
#include "Complex.h"

template<class T>
Complex<T>:: Complex () { re = im = 0.0; }

template<class T>
Complex<T>:: Complex (T re_a, T im_a) {re = re_a; im = im_a; }

template<class T> Complex<T>:: Re () const { return re; } //  getting the real part
template<class T> Complex<T>:: Im () const { return im; }  //   and the imaginary part
template<class T> Complex<T>:: abs () const { return sqrt(re*re + im*im); }
/*
Complex& Complex:: operator= (const Complex& c)
{
  re = c.re;
  im = c.im;
  return *this;
}
Complex operator+ (const Complex& a, const Complex& b) { return Complex (a.re + b.re, a.im + b.im); }
Complex operator- (const Complex& a, const Complex& b) { return Complex (a.re - b.re, a.im - b.im); }

Complex operator* (const Complex& a, const Complex& b) {return Complex(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);}
Complex operator/ (const Complex& a, const Complex& b) {return Complex( (a.re*b.re +a.im*b.im)/(b.re*b.re+b.im*b.im), (-a.re*b.im + a.im*b.re)/(b.re*b.re+b.im*b.im));}
*/
