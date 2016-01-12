    /*
     * The definition module                              
     *                      lib.h                    
     * for the library function common for all C programs.
     */

     // Standard ANSI-C++ include files 


#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;


#define   NULL_PTR   (void *) 0
#define   ZERO       1.0E-10
//#define   INFINITY   1.0E15
#define   UL         unsigned long

         /* a macro used in function pythag() */

static float sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)


     /* Macro definitions for integer arguments only */

#define   SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

    // ******   data declaration  ******* 

typedef   struct  {   // structure definition for execution time   
  unsigned long long int
                          tick,
                           sec,
                           min,
                          hour;
} TID; 

    // Function declarations

TID time_step(int num);
void  **matrix(int, int, int);
void free_matrix(void **);
void rk4(double *, double *, int, double, double, double  *,
	           void (*derivs)(double, double *, double *));
void ludcmp(double **, int, int *, double*);
void lubksb(double **, int, int *, double *);
void tqli(double *, double *, int, double **);
void tred2(double **, int, double *, double *);
double pythag(double, double);
void gauleg(double, double, double *, double *, int);
void jacobi(double** a, double* d, double** v, int n, int& nrot);
double rectangle_rule(double, double, int, double (*func)(double));
double trapezoidal_rule(double, double, int, double (*func)(double));
void spline(double *, double *, int, double, double, double *);
void splint(double *, double *, double *, int, double, double *);
void polint(double *, double *, int, double, double *, double *);
double rtbis(double(*func)(double), double, double, double);
double rtsec(double( *func)(double), double, double, double);
double rtnewt(void ( *funcd)(double, double *, double *), double, double, double);
double zbrent(double( *func)(double), double, double, double);
double ran0(long *);
double ran1(long *);
double ran2(long *);
double ran3(long *);



