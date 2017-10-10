//  Vector class with templates

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
  void   Normalize();
  T VectorNorm1();
  T VectorNorm2();
  T VectorNormInf();
  T MaxMod();
  T ElementofMaxMod();
  int MaxModindex();
  T operator()(const int i) const;
  T& operator()(const int i);
  void Print() const;
  // overloading the equality  operator, here as a template function defned below
  const Vector<T>  operator=(const Vector<T>& v);
  // overloading the binary  operator, here as a friend function
  friend  Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2){
    int min_dim = (v1.Dimension()<v2.Dimension())?v1.Dimension():v2.Dimension();
    Vector x(min_dim);
    for(int i=0;i < min_dim;i++)
      x(i) = v1(i) + v2(i);
    return x;
  };
  //  You can now define many more functions
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
  data = NULL;  // change to c++11 nullptr later
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
  data = NULL;   //  change to c++11 nullptr
}

template <class T>
int Vector<T>::Dimension() const{
  return(dimension);
}




template <class T>
T Vector<T>::VectorNorm1(){
  T sum = 0.0;
  for(int i=0; i < dimension;i++)
    sum += fabs(data[i]);
  return(sum);
}

template <class T>
T Vector<T>::VectorNorm2(){
  T sum = 0.0;
  for(int i=0;i < dimension;i++)
    sum += data[i]*data[i];
  return(sqrt(sum));
}

template <class T>
void Vector<T>::Normalize(){
  T tmp = 1.0/VectorNorm2();
  for(int i=0;i<dimension;i++)
    data[i] = data[i]*tmp;
}

template <class T>
T Vector<T>::VectorNormInf(){
  T maxval = 0.0, tmp;
  
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


template <class T>
const Vector<T> Vector<T>::operator=(const Vector<T>& v){
  dimension = v.Dimension();
  for(int i=0;i<dimension;i++)
    data[i] = v.data[i];
  return *this;
};

