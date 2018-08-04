/*
  Solves the one-particle Schrodinger equation
  for a potential specified in function
  potential(). This example is for the harmonic oscillator in 3d
  It uses Francis algorithm for finding the eigenvalues of a tridiagonal matrix
  Note that the matrix is already tridiagonal
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#define   UL         unsigned long
/* a macro used in function pythag() */
static float sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)
/* Macro definitions for integer arguments only */
#define   SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

using namespace  std;

double potential(double);
void output(double, double, int, double *);
int comp(const double *, const double *);
double ** AllocateMatrix(int, int);
void DeallocateMatrix(double **, int, int); 
double pythag(double, double);
void tqli(double *, double *, int, double **);

// Begin of main program   

int main(int argc, char* argv[])
{
  int       i, j, Dim, lOrbital;
  //  double    RMin, RMax, Step, DiagConst, NondiagConst, OrbitalFactor; 
  // With spherical coordinates RMin = 0 always
  double RMin = 0.0;

  double RMax = 8.0;  lOrbital = 0;  Dim =1000;  
  // Integration step length
  double Step    = RMax/ Dim;
  double DiagConst = 2.0 / (Step*Step);
  double NondiagConst =  -1.0 / (Step*Step);
  double OrbitalFactor = lOrbital * (lOrbital + 1.0);
  
  // local memory for r and the potential w[r] 
  double *r = new double[Dim];
  double *Vpotential = new double[Dim];
  for(i = 0; i <  Dim; i++) {
    r[i] = RMin + (i+1) * Step;
    Vpotential[i] = potential(r[i]) + OrbitalFactor / (r[i] * r[i]);
  }
  // local memory for the diagonalization process 
  double *d = new double[Dim];    // diagonal elements 
  double *e = new double[Dim];    // tri-diagonal off-diagonal elements 
  double **Z = AllocateMatrix(Dim, Dim);
  //  z = (double **) matrix(max_step, max_step, sizeof(double));
  for(i = 0; i < Dim; i++) {
    d[i]    = DiagConst + Vpotential[i];
    e[i]    = NondiagConst;
    Z[i][i] = 1.0;
    for(j = i + 1; j < Dim; j++)  {
      Z[i][j] = 0.0;
    }
  }
  // diagonalize and obtain eigenvalues using Francis' algo
  tqli(d, e, Dim, Z);      
  // Sort eigenvalues as an ascending series 
  qsort(d,(UL) Dim,sizeof(double),
         (int(*)(const void *,const void *))comp);
  // send results to ouput file
  output(RMin , RMax, Dim, d);
  // Free memory
  delete [] r; delete [] Vpotential; delete [] e; delete [] d; 
  DeallocateMatrix(Z, Dim, Dim);
  return 0;
}  //  end of main function


/*
  The function potential()
  calculates and return the value of the 
  potential for a given argument x.
  The potential here is for the harmonic oscillator
*/        

double potential(double x)
{
  return x*x;

} // End: function potential()  


void output(double RMin , double RMax, int Dim, double *d)
{
  int i;
  cout << "RESULTS:" << endl;
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout <<"Rmin = " << setw(15) << setprecision(8) << RMin << endl;  
  cout <<"Rmax = " << setw(15) << setprecision(8) << RMax << endl;  
  cout <<"Number of steps = " << setw(15) << Dim << endl;  
  cout << "Five lowest eigenvalues:" << endl;
  for(i = 0; i < 5; i++) {
    cout << setw(15) << setprecision(8) << d[i] << endl;
  }
}  // end of function output         


int comp(const double *val_1, const double *val_2)
{
  if((*val_1) <= (*val_2))       return -1;
  else  if((*val_1) > (*val_2))  return +1;
  else                     return  0; 
} // End: function comp() 


// Allocate memory for a matrix and initialize the elements to zero

double ** AllocateMatrix(int m, int n){
  double ** Matrix;
  Matrix = new double*[m];
  for(int i=0;i<m;i++){
    Matrix[i] = new double[n];
    for(int j=0;j<m;j++)
      Matrix[i][j] = 0.0;
  }
  return Matrix;
}

// Free memory

void DeallocateMatrix(double ** Matrix, int m, int n){
  for(int i=0;i<m;i++)
    delete[] Matrix[i];
  delete[] Matrix;
}


    /*
    ** The function
    **                 tqli()
    ** determine eigenvalues and eigenvectors of a real symmetric
    ** tri-diagonal matrix, or a real, symmetric matrix previously
    ** reduced by function tred2[] to tri-diagonal form. On input,
    ** d[] contains the diagonal element and e[] the sub-diagonal
    ** of the tri-diagonal matrix. On output d[] contains the
    ** eigenvalues and  e[] is destroyed. If eigenvectors are
    ** desired z[][] on input contains the identity matrix. If
    ** eigenvectors of a matrix reduced by tred2() are required,
    ** then z[][] on input is the matrix output from tred2().
    ** On output, the k'th column returns the normalized eigenvector
    ** corresponding to d[k]. 
    ** The function is modified from the version in Numerical recipe.
    */

void tqli(double *d, double *e, int n, double **z)
{
   register int   m,l,iter,i,k;
   double         s,r,p,g,f,dd,c,b;

   for(i = 1; i < n; i++) e[i-1] = e[i];
     e[n] = 0.0;
   for(l = 0; l < n; l++) {
      iter = 0;
      do {
         for(m = l; m < n-1; m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
            if((double)(fabs(e[m])+dd) == dd) break;
         }
         if(m != l) {
            if(iter++ == 30) {
               printf("\n\nToo many iterations in tqli.\n");
               exit(1);
            }
            g = (d[l+1] - d[l])/(2.0 * e[l]);
            r = pythag(g,1.0);
            g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
            s = c = 1.0;
            p = 0.0;
            for(i = m-1; i >= l; i--) {
               f      = s * e[i];
               b      = c*e[i];
               e[i+1] = (r=pythag(f,g));
               if(r == 0.0) {
                  d[i+1] -= p;
                  e[m]    = 0.0;
                  break;
               }
               s      = f/r;
               c      = g/r;
               g      = d[i+1] - p;
               r      = (d[i] - g) * s + 2.0 * c * b;
               d[i+1] = g + (p = s * r);
               g      = c * r - b;
               for(k = 0; k < n; k++) {
                  f         = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i]   = c * z[k][i] - s * f;
               } /* end k-loop */
            } /* end i-loop */
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
         } /* end if-loop for m != 1 */
      } while(m != l);
   } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */
   

double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
