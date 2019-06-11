/*
  Solves the one-particle Schrodinger equation
  for a potential specified in function
  potential(). This example is for the hydrogenatom
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace  std;
// output file as global variable
ofstream ofile;  

// function declarations 

void initialise(double&, double&, int&, int&) ;
double potential(double);
int comp(const double *, const double *);
void output(double, double, int, double *);

int main(int argc, char* argv[])
{
  int       i, j, max_step, orb_l;
  double    r_min, r_max, step, const_1, const_2, orb_factor, 
            *e, *d, *w, *r, **z;
  char *outfilename;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] << 
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  ofile.open(outfilename); 
  //   Read in data 
  initialise(r_min, r_max, orb_l, max_step);
  // initialise constants
  step    = (r_max - r_min) / max_step; 
  const_2 = -1.0 / (step * step);
  const_1 =  2.0 * const_2;
  orb_factor = orb_l * (orb_l + 1);
  
  // local memory for r and the potential w[r] 
  r = new double[max_step + 1];
  w = new double[max_step + 1];
  for(i = 0; i <= max_step; i++) {
    r[i] = r_min + i * step;
    w[i] = potential(r[i]) + orb_factor / (r[i] * r[i]);
  }
  // local memory for the diagonalization process 
  d = new double[max_step];    // diagonal elements 
  e = new double[max_step];    // tri-diagonal off-diagonal elements 
  z = (double **) matrix(max_step, max_step, sizeof(double));
  for(i = 0; i < max_step; i++) {
    d[i]    = const_1 + w[i + 1];
    e[i]    = const_2;
    z[i][i] = 1.0;
    for(j = i + 1; j < max_step; j++)  {
      z[i][j] = 0.0;
    }
  }
  // diagonalize and obtain eigenvalues
  tqli(d, e, max_step - 1, z);      
  // Sort eigenvalues as an ascending series 
  qsort(d,(UL) max_step - 1,sizeof(double),
         (int(*)(const void *,const void *))comp);
  // send results to ouput file
  output(r_min , r_max, max_step, d);
  delete [] r; delete [] w; delete [] e; delete [] d; 
  free_matrix((void **) z); // free memory
  ofile.close();  // close output file
  return 0;
} // End: function main() 

/*
  The function potential()
  calculates and return the value of the 
  potential for a given argument x.
  The potential here is for the hydrogen atom
*/        

double potential(double x)
{
   return x*x;

} // End: function potential()  

/*
  The function   int comp()                  
  is a utility function for the library function qsort()
  to sort double numbers after increasing values.
*/       

int comp(const double *val_1, const double *val_2)
{
  if((*val_1) <= (*val_2))       return -1;
  else  if((*val_1) > (*val_2))  return +1;
  else                     return  0; 
} // End: function comp() 


void initialise(double& r_min, double& r_max, int& orb_l, int& max_step) 
{
  cout << "Min vakues of R = ";
  cin >> r_min;
  cout << "Max value of R = ";
  cin >> r_max;
  cout << "Orbital momentum = ";
  cin >> orb_l;
  cout << "Number of steps = ";
  cin >> max_step;
}  // end of function initialise   




void output(double r_min , double r_max, int max_step, double *d)
{
  int i;
  ofile << "RESULTS:" << endl;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile <<"R_min = " << setw(15) << setprecision(8) << r_min << endl;  
  ofile <<"R_max = " << setw(15) << setprecision(8) << r_max << endl;  
  ofile <<"Number of steps = " << setw(15) << max_step << endl;  
  ofile << "Five lowest eigenvalues:" << endl;
  for(i = 0; i < 5; i++) {
    ofile << setw(15) << setprecision(8) << d[i] << endl;
  }
}  // end of function output         


