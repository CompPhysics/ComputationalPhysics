#ifndef _vectorclass
#define _vectorclass
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
using namespace std;

template<class T>
class Vector{
 private:
  int   dimension;
  T *data;
  
 public:
  Vector();
  Vector(int dim);
  Vector(const Vector<T>& v);
  ~Vector();
  
  int    Dimension() const;
  //  T Length();     /* Euclidean Norm of the Vector */
  void   Normalize();

  T Norm_l1();
  T Norm_l2();
  T Norm_linf();
  T MaxMod();
  T ElementofMaxMod();
  int MaxModindex();
  T operator()(const int i) const;
  T& operator()(const int i);
  void Print() const;
  T Length();     /* Euclidean Norm of the Vector */
  Vector<T>& operator=(const Vector<T>& v);

// Unitary operator -
Vector<T> operator-(const Vector<T>& v);

// Binary operator +,-
Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2);
Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2);

// Vector Scaling (multiplication by a scalar : defined commutatively)
Vector<T> operator*(const T s, const Vector<T>& v);
Vector<T> operator*(const Vector<T>& v, const T s);
// Vector Scaling (division by a scalar)
Vector<T> operator/(const Vector<T>& v, const T s);
};


#endif


template <class T>
void Vector<T>::Print() const{
  cout << endl;
  cout << "[ ";
  if(dimension>0)
    cout << data[0];
  for(int i=1;i<dimension;i++)
    cout << "; " << data[i];
  cout << " ]" << endl;
}

template <class T>
Vector<T>::Vector(){
  dimension = 0;
  data = NULL;
}

template <class T>
Vector<T>::Vector(int dim){
  dimension = dim;
  data = new T[dimension];

  for(int i=0;i<dimension;i++)
    data[i] = 0.0;
}

template <class T>
Vector<T>::Vector(const Vector &v){
  dimension = v.Dimension();
  data = new T[dimension];

  for(int i=0;i<dimension;i++)
    data[i] = v.data[i];
}


template <class T>
Vector<T>::~Vector(){
  dimension = 0;
  delete[] data;
  data = NULL;
}



template <class T>
int Vector<T>::Dimension() const{
  return(dimension);
}



template <class T>
T Vector<T>::Norm_l1(){
  T sum = 0.0;
  for(int i=0; i < dimension;i++)
    sum += fabs(data[i]);
  return(sum);
}

template <class T>
T Vector<T>::Norm_l2(){
  T sum = 0.0;
  for(int i=0;i < dimension;i++)
    sum += data[i]*data[i];
  return(sqrt(sum));
}

template <class T>
void Vector<T>::Normalize(){
  T tmp = 1.0/Norm_l2();
  for(int i=0;i<dimension;i++)
    data[i] = data[i]*tmp;
}

template <class T>
T Vector<T>::Norm_linf(){
  T maxval = 0.0,tmp;
  
  for(int i=0;i<dimension;i++){
    tmp = fabs(data[i]);
    maxval = (maxval > tmp)?maxval:tmp;
  }
  return(maxval);
}

template <class T>
T Vector<T>::MaxMod(){
  T maxm = -1.0e+10;

  for(int i=0; i<dimension; i++)
    maxm = (maxm > fabs(data[i]))?maxm:fabs(data[i]);
  
  return maxm;
}

template <class T>
T Vector<T>::ElementofMaxMod(){
  return(data[MaxModindex()]);
}

template <class T>
int Vector<T>::MaxModindex(){
  T maxm = -1.0e+10;
  int maxmindex = 0;

  for(int i=0; i<dimension; i++){
    if(maxm<fabs(data[i])){
      maxm = fabs(data[i]);
      maxmindex = i;
    }
  }
  
  return maxmindex;
}


template <class T>
int min_dimension(const Vector<T>& v1, const Vector<T>& v2){
  int min_dim = (v1.Dimension()<v2.Dimension())?v1.Dimension():v2.Dimension();
  return(min_dim);
}

template <class T>
T dot(const Vector<T>& u, const Vector<T>& v){
  T sum = 0.0;
  int min_dim = min_dimension(u,v);

  for(int i=0;i<min_dim;i++)
    sum += u(i)*v(i);
  
  return sum; 
}

template <class T>
T dot(int N, const Vector<T>& u, const Vector<T>& v){
  T sum = 0.0;

  for(int i=0;i<N;i++)
    sum += u(i)*v(i);
  
  return sum;
}

template <class T>
T dot(int N, double *a, double *b){
  T sum = 0.0;
  
  for(int i=0;i<N;i++)
    sum += a[i]*b[i];

  return sum;
}


template <class T>
T Vector<T>::operator()(const int i) const{
  if(i>=0 && i<dimension)
    return data[i];

  cerr << "Vector::Invalid index " << i << " for Vector of dimension " << dimension << endl;
  return(0);
}


template <class T>
T& Vector<T>::operator()(const int i){
  if(i>=0 && i<dimension)
    return data[i];

  cerr << "Vector::Invalid index " << i << " for Vector of dimension " << dimension << endl;
  return(data[0]);
}





