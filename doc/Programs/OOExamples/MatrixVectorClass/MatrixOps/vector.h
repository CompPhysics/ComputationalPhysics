#ifndef _vectorclass
#define _vectorclass


#include <cmath>
#include <iostream>
using namespace  std;



class Vector{
 private:
  int   dimension;
  double *data;
  
 public:
  Vector();
  Vector(int dim);
  Vector(const Vector& v);
  ~Vector();
  
  void Initialize(int dim);
  int    Dimension() const;
  double Length();     /* Euclidean Norm of the Vector */
  void   Normalize();

  double Norm_l1();
  double Norm_l2();
  double Norm_linf();
  double MaxMod();
  double ElementofMaxMod();
  int MaxModindex();
  
  //************************
  // User Defined Operators
  //************************
  int operator==(const Vector& v) const;
  int operator!=(const Vector& v) const;
  Vector & operator=(const Vector& v);

  double  operator()(const int i) const;
  double& operator()(const int i);

  void Print() const;
  void Initialize(double a);
  void Initialize(double *v);
};


// Unitary operator -
Vector operator-(const Vector& v);

// Binary operator +,-
Vector operator+(const Vector& v1, const Vector& v2);
Vector operator-(const Vector& v1, const Vector& v2);

// Vector Scaling (multiplication by a scalar : defined commutatively)
Vector operator*(const double s, const Vector& v);
Vector operator*(const Vector& v, const double s);

// Vector Scaling (division by a scalar)
Vector operator/(const Vector& v, const double s);


/********************************/
/*   Function Declarations      */
/********************************/

int min_dimension(const Vector& u, const Vector& v);
double dot(const Vector& u, const Vector& v); 
double dot(int N, double *a, double *b);
double dot(int N, const Vector &u, const Vector &v); 
void Swap(double &a, double &b);
double Sign(double x);



#endif


