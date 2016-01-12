%module pylib_cpp

%{
#define SWIG_FILE_WITH_INIT
#include "lib.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%include "typemaps.i"


/* Wrap ludcmp */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* a, int a1, int a2)};
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* indx, int indx1)};
%apply (double* INOUT) {(double* d)};

%rename (ludcmp) ludcmp_wrapper;

%exception ludcmp_wrapper {
	$action
	if(PyErr_Occurred()) SWIG_fail;
}

%inline %{
void ludcmp_wrapper(double* a, int a1, int a2, int* indx, int indx1, double* d) {

	if (a1 != a2 || a1 != indx1) {
		PyErr_Format(PyExc_ValueError, "Sizes dont match: a is a (%d x %d) array, indx length %d", a1,a2,indx1);
		return;
	}

	//Convert single-pointer-array a into double-pointer-array A
	double** A = new double* [a1];
	for (int i = 0; i < a1; i++) {
		A[i] = a + a2*i;
	}

	ludcmp(A, a1, indx, d);
	return;
}
%}

%clear (double* a, int a1, int a2);
%clear (int* indx, int indx1);
%clear (double* d);

/* Wrap lubksb */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* a, int a1, int a2)};
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* indx, int indx1)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* d, int d1)};

%rename (luBackSubst) lubksb_wrapper;

%exception lubksb_wrapper {
	$action
	if(PyErr_Occurred()) SWIG_fail;
}
%inline %{
void lubksb_wrapper(double* a, int a1, int a2, int* indx, int indx1, double* d, int d1) {

	if (a1 != a2 || a1 != indx1 || a1 != d1) {
		PyErr_Format(PyExc_ValueError, "Sizes dont match: a is a (%d x %d) array, indx length %d, d length %d", a1,a2,indx1, d1);
		return;
	}

	//Convert single-pointer-array a into double-pointer-array A
	double** A = new double* [a1];
	for (int i = 0; i < a1; i++) {
		A[i] = a + a2*i;
	}

	lubksb(A, a1, indx, d);
	return;
}
%}

%clear (double* a, int a1, int a2);
%clear (int* indx, int indx1);
%clear (double* d, int d1);

/* Wrap gauleg */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* x, int xn)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* w, int wn)};

%rename (gausLegendre) gauleg_wrapper;

%exception gauleg_wrapper {
	$action
	if(PyErr_Occurred()) SWIG_fail;
}
%inline %{
void gauleg_wrapper(double x1, double x2, double* x, int xn, double* w, int wn, int n) {
	if (n != xn || n != wn) {
		PyErr_Format(PyExc_ValueError, "Sizes dont match: xn = %d, wn = %d, but n = %n", xn,wn,n);
		return;
	}
	gauleg(x1,x2,x,w,n);
	return;
}
%}

%clear (double* x, int xn);
%clear (double* w, int xw);

/* Wrap tred2 */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* a, int a1, int a2)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* d, int d1)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* e, int e1)};

%rename (tred2) tred2_wrapper;

%exception tred2_wrapper {
	$action
	if(PyErr_Occurred()) SWIG_fail;
}
%inline %{
void tred2_wrapper(double* a, int a1, int a2, int n, double* d, int d1, double* e, int e1) {
	if (n != a1 || n != a2 || n != d1 || n != e1) {
		PyErr_Format(PyExc_ValueError, "Sizes dont match: n = %d, a1 = %d, a2 = %d, d1 = %d, e1 = %d", n,a1,a2,d1,e1);
		return;
	}

	//Convert single-pointer-array a into double-pointer-array A
	double** A = new double* [a1];
	for (int i = 0; i < a1; i++) {
		A[i] = a + a2*i;
	}

	tred2(A,n,d,e);
	return;
}
%}

%clear (double* a, int a1, int a2);
%clear (double* d, int d1);
%clear (double* e, int e1);

/* Wrap tqli */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* d, int d1)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* e, int e1)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* z, int z1, int z2)};

%rename (tqli) tqli_wrapper;

%exception tqli_wrapper {
	$action
	if(PyErr_Occurred()) SWIG_fail;
}
%inline %{
void tqli_wrapper(double* d, int d1, double* e, int e1, int n, double* z, int z1, int z2) {
	if (n != d1 || n != e1 || n != z1 || n != z2) {
		PyErr_Format(PyExc_ValueError, "Sizes dont match: n = %d, d1 = %d, e1 = %d, z1 = %d, z2 = %d", n,d1,e1,z1,z2);
		return;
	}

	//Convert single-pointer-array z into double-pointer-array Z
	double** Z = new double* [z1];
	for (int i = 0; i < z1; i++) {
		Z[i] = z + z2*i;
	}

	tqli(d,e,n,Z);
	return;
}
%}

%clear (double* d, int d1);
%clear (double* e, int e1);
%clear (double* z, int z1, int z2);

/* Wrap jacobi */
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* a, int a1, int a2)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* d, int d1)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* v, int v1, int v2)};
%apply (double* OUTPUT) {(int &nrot)};

%rename (jacobi) jacobi_wrapper;

%exception jacobi_wrapper {
	$action
	if(PyErr_Occurred()) SWIG_fail;
}

%inline %{
void jacobi_wrapper(double* a, int a1, int a2, double* d, int d1, double* v, int v1, int v2, int n, int &nrot) {
	if (a1 != a2 || a1 != d1 || v1 != v2 || d1 != v1) {
		PyErr_Format(PyExc_ValueError, "Sizes dont match: a is a (%d x %d) array, d length %d, v is v (%d x %d) array, ", a1,a2,d,v1,v2);
		return;
	}

	//Convert single-pointer-array a into double-pointer-array A
	double** A = new double* [a1];
	for (int i = 0; i < a1; i++) {
		A[i] = a + a2*i;
	}
	double** V = new double* [v1];
	for(int i = 0; i < v1; i ++){
		V[i] = v + v2*i;
	}
	jacobi(A, d, V, n, nrot);
	return;
}
%}

%clear (double* a, int a1, int a2);
%clear (int* d, int d);
%clear (double* v, int v1, int v2);
%clear (int& d);