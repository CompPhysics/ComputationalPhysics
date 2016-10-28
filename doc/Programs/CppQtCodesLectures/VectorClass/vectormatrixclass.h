#ifndef _vectormatrixclass
#define _vectormatrixclass


#include <cmath>
#include <iostream>
using namespace  std;



class Point;
class Vector;
class Matrix;


/********************************/
/*        Point Class        */
/********************************/

class Point{
 private:
  int   dimension;
  double *data;
  
 public:
  Point(int dim);
  Point(const Point& v);
  ~Point();
  
  int    Dimension() const;

  //************************
  // User Defined Operators
  //************************
  int operator==(const Point& v) const;
  int operator!=(const Point& v) const;
  Point & operator=(const Point& v);

  double  operator()(const int i) const;
  double& operator()(const int i);

  void Print() const;
};



/********************************/
/*        Vector Class        */
/********************************/

class Vector{
 private:
  int   dimension;
  double *data;
  
 public:
  Vector();
  Vector(int dim);
  Vector(const Vector& v);
  Vector(int col, const Matrix &A);
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



/********************************/
/*        Matrix Class        */
/********************************/

class Matrix {
private:
  int rows, columns;
  double **data;
  
public:

  Matrix(int dim);
  Matrix(int rows1, int columns1);
  Matrix(const Matrix& m);
  Matrix(int num_vectors, const Vector * q);
  Matrix(int rows1, int columns1, double **rowptrs);
  ~Matrix();

  int Rows() const;
  int Columns() const;
  double ** GetPointer();
  void GetColumn(int col, Vector &x);
  void GetColumn(int col, Vector &x, int rowoffset);
  void PutColumn(int col, const Vector &x);
  double Norm_l1();
  double Norm_linf();

  //************************
  // User Defined Operators
  //************************
  Matrix& operator=(const Matrix& m);
  double operator()(const int i, const int j) const;
  double& operator()(const int i, const int j);

  double MaxModInRow(int row);
  double MaxModInRow(int row, int starting_column);
  int MaxModInRowindex(int row);
  int MaxModInRowindex(int row, int starting_column);
 
  double MaxModInColumn(int column);
  double MaxModInColumn(int column, int starting_row);
  int MaxModInColumnindex(int column);
  int MaxModInColumnindex(int column, int starting_row);

  void RowSwap(int row1, int row2);

  void Print() const;

};


/********************************/
/*   Operator Declarations      */
/********************************/

// Unitary operator -
Vector operator-(const Vector& v);

// Binary operator +,-
Vector operator+(const Vector& v1, const Vector& v2);
Vector operator-(const Vector& v1, const Vector& v2);

// Vector Scaling (multiplication by a scaler : defined commutatively)
Vector operator*(const double s, const Vector& v);
Vector operator*(const Vector& v, const double s);

// Vector Scaling (division by a scaler)
Vector operator/(const Vector& v, const double s);

Vector operator*(const Matrix& A, const Vector& x); 


/********************************/
/*   Function Declarations      */
/********************************/

int min_dimension(const Vector& u, const Vector& v);
double dot(const Vector& u, const Vector& v); 
double dot(int N, double *a, double *b);
double dot(int N, const Vector &u, const Vector &v); 
void Swap(double &a, double &b);
double Sign(double x);

/* Misc. useful functions to have */
double log2(double x);
double GammaF(double x);
int Factorial(int n);
double ** CreateMatrix(int m, int n);
void DestroyMatrix(double ** mat, int m, int n); 

int ** ICreateMatrix(int m, int n);
void IDestroyMatrix(int ** mat, int m, int n);

#endif


