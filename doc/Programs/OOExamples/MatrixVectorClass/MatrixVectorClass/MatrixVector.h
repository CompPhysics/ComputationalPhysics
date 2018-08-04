#if !defined(MatrixVector_H)
#define MatrixVector_H

// file        : MyArrah.h
// description : the class MyArray is a templated 1D or 2D array, specialized 
//               for containing numerical data.

#include <iostream>

template< typename T >
class MatVec
{
    T* A;                                   // the data
    int dimensions;                         // number of dimensions (max is 2)
    int length[2];                          // lengths of each dimension

    // needed for optimizations:
    int length0;
    int length0p1;

    // these three functions are used by constructors, destructors and redim:
    T* allocate(int length1);
    T* allocate(int length1, int length2);
    void deallocate();

    // transform fake indicies to real:
    int transform_index(int i, int j) const;

public:
    MatVec();
    MatVec(int n1);                   // constructor for 1D array
    MatVec(int n1, int n2);           // constructor for 2D array
    MatVec(const MatVec<T>& array);  // copy constructor
    ~MatVec();                        // desctructor

    void redim(int n1);                // redim 1D array
    void redim(int n1, int n2);        // redim 2D array

    // return the size of the arrays dimensions
    int size() const;                  // length of 1D array
    int size(int dimension) const;

    bool indexOk(int i) const;         // check if index is ok to use
    bool indexOk(int i, int j) const;  // check if indicies are ok to use

    // operator() for 1D array
    const T& operator()(int i) const;
    T& operator()(int i);

    // operator() for 2D array
    const T& operator()(int i, int j) const;
    T& operator()(int i, int j);

    MatVec& operator=(const MatVec& v);  // assigment operator
                                           // (not implemented yet)
    // returns pointers to the data
    const T* getPtr() const;
    T* getPtr();

    void print(std::ostream& os);
};


template< typename T >
inline int MatVec<T>::  transform_index(int i, int j) const
{ return i + j*length0 - length0p1; }


template< typename T >
inline T& MatVec<T>:: operator()(int i)
{
#ifdef SAFETY_CHECKS
    indexOk(i);
#endif

    return  A[i - 1];
}


template< typename T >
inline const T& MatVec<T>:: operator()(int i) const
{
#ifdef SAFETY_CHECKS
    indexOk(i);
#endif
    
    return  A[i - 1];
}


template< typename T >
inline T& MatVec<T>:: operator()(int i, int j)
{
#ifdef SAFETY_CHECKS
    indexOk(i, j);
#endif
    
#ifdef NO_NESTED_INLINES
    return A[i + j*length0 - length0p1];
#else
    return A[transform_index(i, j)];
#endif
}


template< typename T >
inline const T& MatVec<T>:: operator()(int i, int j) const
{
#ifdef SAFETY_CHECKS
    indexOk(i, j);
#endif
    
#ifdef NO_NESTED_INLINES
    return A[i + j*length0 - length0p1];
#else
    return  A[transform_index(i, j)];
#endif
}

// since this is a template class, we must make the body of the
// functions available in this header file:
#include "MatrixVector.cpp"

#endif
