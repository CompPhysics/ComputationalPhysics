#ifndef MatrixVector_CPP
#define MatrixVector_CPP

// file        : MyArray.cpp
// description : implentation of class MyArray

#include <memory>

#include "MatrixVector.h"


template< typename T >
MatVec<T>:: MatVec() 
    : A(0)
{}


template< typename T >
MatVec<T>:: MatVec(int n1)
{
    A = allocate(n1);
}


template< typename T >
MatVec<T>:: MatVec(int n1, int n2)
{
    A = allocate(n1, n2);
}


template< typename T >
MatVec<T>:: MatVec(const MatVec<T>& array)
{
    if (array.dimensions == 2) {
	A = allocate(array.length[0], array.length[1]);

	for (int i = 0; i != (length[0] * length[1]); ++i)
	    A[i] = array.A[i];
    }
    else if (array.dimensions == 1) {
	A = allocate(array.length[0]);
	
	for (int i = 0; i != length[0]; ++i)
	    A[i] = array.A[i];
    }
    else {
	std::cerr << "MatVec::MatVec(const MatVec<T>&) -- "
		  << "illegal dimension " << array.dimension << std::endl;
	exit(1);
    }
}


template< typename T >
MatVec<T>:: ~MatVec()
{
    deallocate();
}


template< typename T >
void MatVec<T>:: redim(int n1)
{
    if (length[0] == n1)
	return;

    deallocate();
    A = allocate(n1);
}


template< typename T >
void MatVec<T>:: redim(int n1, int n2)
{
    if (length[0] == n1 && length[1] == n2)
	return;
    
    deallocate();
    A = allocate(n1, n2);
}


template< typename T >
int MatVec<T>:: size() const
{ return length[0]; }


template< typename T >
int MatVec<T>:: size(int dimension) const
{
    return (dimension == 1) ? length[0] : length[1];
}


template< typename T >
T* MatVec<T>:: allocate(int n1)
{
    T* ptr = 0;
    length[0] = n1;
    dimensions = 1;
    
    if (n1 < 0) {
	std::cerr << "MatVec::allocate -- illegal length " << n1
		  << std::endl;

	exit(1);
    }
    else if (n1 == 0) {
	ptr = 0;
	return ptr;
    }
    
    try {
	ptr = new T[n1];
    } 
    catch (std::bad_alloc&) {
	std::cerr << "MatVec::allocate -- unable to allocate array of length "
		  << n1 << std::endl;
	exit(1);
    }    
    return ptr;
}

template< typename T >
T* MatVec<T>:: allocate(int n1, int n2)
{
    T* ptr = 0;
    length[0] = n1;
    length[1] = n2;
    length0 = n1;       // special variable for efficient indexing
    length0p1 = n1+1;   // special variable for efficient indexing
    dimensions = 2;

    if (n1 < 0 || n2 < 0) {
	std::cerr << "MatVec::allocate -- illegal length " << n1
		  << ", " << n2 << std::endl;
	
	return 0;
    }
    
    if (n1 == 0) {
	ptr = 0;
	return ptr;
    }

    try {
	ptr = new T [n1 * n2];
    }
    catch (std::bad_alloc&) {
	std::cerr << "MatVec::allocate -- unable to allocate array of length "
		  << length[0] << ", " << length[1] << std::endl;
	exit(1);
    }
    
    return ptr;
}


template< typename T >
void MatVec<T>:: deallocate()
{
    delete [] A;
	
    length[0] = 0;
    length[1] = 0;
}


template< typename T >
T* MatVec<T>:: getPtr()
{ return A; }


template< typename T >
bool MatVec<T>:: indexOk(int i) const
{
  if (A == 0) {
      std::cerr << "vector index check; current array is of length 0,"
		<< " call redim!" << std::endl;
      
      return false;
  }
  else if (i < 1 || i > length[0]) {
      std::cerr << "vector index check; index " << i 
		<< "out of bounds in array (1:" << length[0] << ")" << std::endl;
      
      return false;
  }
  else
      return true;  // valid index!
}


template< typename T >
bool MatVec<T>:: indexOk(int i, int j) const
{
    if (A == 0) {
	std::cerr << "vector index check; current array is of length 0,"
		  << " call redim!" << std::endl;
	
	return false;
    }
    else if (i < 1 || i > length[0]
	     || j < 1 || j > length[1]) {
	std::cerr << "vector index check; index " << i << "," << j
		  << "out of bounds in array (" << 1 << ":"
		  << 1 + length[0] << ")(" << 1 << ":"
		  << 1 + length[1] << std::endl;
	
	return false;
    }
    else
        return true;  // valid index!
}


template< typename T >
void MatVec<T>:: print(std::ostream& os)
{
    for (int i = 1; i != size(1); ++i) {
	for (int j = 1; j != size(2); ++j)
	    os << (*this)(i, j) << " ";
	os << "\n";
    }
}


template<typename T> 
std::ostream& operator<<(std::ostream& os, MatVec<T>& arr)
{
    arr.print(os);
    return os;
}


#endif
