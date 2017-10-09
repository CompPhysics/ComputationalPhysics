#include "vector.h"


Vector::Vector(){
  dimension = 0;
  data = NULL;
}


Vector::Vector(int dim){
  dimension = dim;
  data = new double[dimension];

  for(int i=0;i<dimension;i++)
    data[i] = 0.0;
}


Vector::Vector(const Vector &v){
  dimension = v.Dimension();
  data = new double[dimension];

  for(int i=0;i<dimension;i++)
    data[i] = v.data[i];
}



Vector::~Vector(){
  dimension = 0;
  delete[] data;
  data = NULL;
}


void Vector::Initialize(int dim){
  if(dimension!=0)
    delete[] data;

  dimension = dim;
  data = new double[dimension];
  
  for(int i=0;i<dimension;i++)
    data[i] = 0.0;
}


int Vector::Dimension() const{
  return(dimension);
}


double Vector::operator()(const int i) const{
  if(i>=0 && i<dimension)
    return data[i];

  cerr << "Vector::Invalid index " << i << " for Vector of dimension " << dimension << endl;
  return(0);
}



double& Vector::operator()(const int i){
  if(i>=0 && i<dimension)
    return data[i];

  cerr << "Vector::Invalid index " << i << " for Vector of dimension " << dimension << endl;
  return(data[0]);
}


Vector& Vector::operator=(const Vector &v) {
  dimension = v.Dimension();
  for(int i=0;i<dimension;i++)
    data[i] = v.data[i];
  return *this;
};

void Vector::Print() const{
  cout << endl;
  cout << "[ ";
  if(dimension>0)
    cout << data[0];
  for(int i=1;i<dimension;i++)
    cout << "; " << data[i];
  cout << " ]" << endl;
}


double Vector::Norm_l1(){
  double sum = 0.0;
  for(int i=0;i<dimension;i++)
    sum += fabs(data[i]);
  return(sum);
}


double Vector::Norm_l2(){
  double sum = 0.0;
  for(int i=0;i<dimension;i++)
    sum += data[i]*data[i];
  return(sqrt(sum));
}

void Vector::Normalize(){
  double tmp = 1.0/Norm_l2();
  for(int i=0;i<dimension;i++)
    data[i] = data[i]*tmp;
}


double Vector::Norm_linf(){
  double maxval = 0.0,tmp;
  
  for(int i=0;i<dimension;i++){
    tmp = fabs(data[i]);
    maxval = (maxval > tmp)?maxval:tmp;
  }
  return(maxval);
}

double Vector::MaxMod(){
  double maxm = -1.0e+10;

  for(int i=0; i<dimension; i++)
    maxm = (maxm > fabs(data[i]))?maxm:fabs(data[i]);
  
  return maxm;
}

double Vector::ElementofMaxMod(){
  return(data[MaxModindex()]);
}


int Vector::MaxModindex(){
  double maxm = -1.0e+10;
  int maxmindex = 0;

  for(int i=0; i<dimension; i++){
    if(maxm<fabs(data[i])){
      maxm = fabs(data[i]);
      maxmindex = i;
    }
  }
  
  return maxmindex;
}

void Vector::Initialize(double a){
  for(int i=0; i<dimension; i++)
    data[i] = a;
}

void Vector::Initialize(double *v){
  for(int i=0; i<dimension; i++)
    data[i] = v[i];
}


/****************************************************************/
/*                 Operator Definitions                         */
/****************************************************************/


Vector operator-(const Vector& v){
  Vector x(v.Dimension());
  for(int i=0;i<v.Dimension();i++)
    x(i) = -v(i);
  return x;
}


Vector operator+(const Vector& v1, const Vector& v2){
  int min_dim = min_dimension(v1,v2);
  Vector x(min_dim);
  for(int i=0;i<min_dim;i++)
    x(i) = v1(i) + v2(i);
  return x;
}


Vector operator-(const Vector& v1, const Vector& v2){
  int min_dim = min_dimension(v1,v2);
  Vector x(min_dim);
  for(int i=0;i<min_dim;i++)
    x(i) = v1(i) - v2(i);
  return x;
}


Vector operator/(const Vector& v, const double s) {
  Vector x(v.Dimension());
  for(int i=0;i<v.Dimension();i++)
    x(i) = v(i)/s;
  return x;
}



Vector operator*(const double s, const Vector &v) {
  Vector x(v.Dimension());
  for(int i=0;i<v.Dimension();i++)
    x(i) = s*v(i);
  return x;
}


Vector operator*(const Vector& v, const double s) {
  Vector x(v.Dimension());
  for(int i=0;i<v.Dimension();i++)
    x(i) = s*v(i);
  return x;
}


/****************************************************************/
/*                 Function Definitions                         */
/****************************************************************/

int min_dimension(const Vector& v1, const Vector& v2){
  int min_dim = (v1.Dimension()<v2.Dimension())?v1.Dimension():v2.Dimension();
  return(min_dim);
}


double dot(const Vector& u, const Vector& v){
  double sum = 0.0;
  int min_dim = min_dimension(u,v);

  for(int i=0;i<min_dim;i++)
    sum += u(i)*v(i);
  
  return sum; 
}


double dot(int N, const Vector& u, const Vector& v){
  double sum = 0.0;

  for(int i=0;i<N;i++)
    sum += u(i)*v(i);
  
  return sum;
}


double dot(int N, double *a, double *b){
  double sum = 0.0;
  
  for(int i=0;i<N;i++)
    sum += a[i]*b[i];

  return sum;
}



#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

template<class T>
class Array{
		//! Copy constructor
    Array(const Array<T>& array);    
    

    Array<T>& operator=(const Array<T>& array);

    Array<T> operator+(const Array<T>& array);
      
		Array<T> operator-(const Array<T>& array)const; /// w=u-v;
		
		
		Array<T>& operator+=(const Array<T>& w);


    Array<T>& operator-=(const Array<T>& w);
    
    
    Array<T>& operator*=(double scalar);
		
    Array<T>& operator/=(double scalar);	
		
    //! Index operators
    const T& operator()(int i)const;	
    const T& operator()(int i, int j)const;	
		const T& operator()(int i, int j, int k)const;	
		const T& operator()(int i, int j, int k, int l)const;		
		const T& operator()(int i, int j, int k, int l, int m)const;		
		const T& operator()(int i, int j, int k, int l, int m, int n)const;		
		
    T& operator()(int i);
		T& operator()(int i, int j);
		T& operator()(int i, int j, int k);	
		T& operator()(int i, int j, int k, int l);		
		T& operator()(int i, int j, int k, int l, int m);		
		T& operator()(int i, int j, int k, int l, int m, int n);		
		
		
   		
		/**************************************************************/
		/*								FRIEND FUNCTIONS 														*/
		/**************************************************************/
 		//! Unary operator +
		template <class T2>
    friend Array<T> operator+ (const Array<T>&);                 // u = + v
 
		//! Unary operator -
		template <class T2>
    friend Array<T> operator-(const Array<T>&);                 // u = - v
    
		template <class T2>
    friend Array<T> operator/(const Array<T>&, double);         // u = v/a 
		
		
// Destructor
template <class T>
inline Array<T>::~Array(){delete[] data;}

// Index operators
template <class T>
inline Array<T> operator+(const Array<T>& v){     // u = + v
	return v;
}



template <class T>
inline Array<T> operator-(const Array<T>& v){      // u = - v
	return Array<T>(v.size[0],v.size[1]) -v;
}



template <class T>
inline Array<T> operator*(const Array<T>& v, double scalar){   // u = v*a
  return Array<T>(v) *= scalar;
}


template <class T>
inline Array<T> operator*(double scalar, const Array<T>& v){   // u = a*v
  return v*scalar;  // Note the call to postmultiplication operator defined above
}


template <class T>
inline Array<T> operator/(const Array<T>& v, double scalar){ 
  if(!scalar) std::cout << "Division by zero!" << std::endl;
  return (1.0/scalar)*v;
}

