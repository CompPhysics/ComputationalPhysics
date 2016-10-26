/* 
   Program to solve the two-dimensional Ising model 
   with zero external field using MPI
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include "vectormatrixclass.h"
using namespace  std;

// output file
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void initialize(int, Matrix, double&, double&);
// The metropolis algorithm 
void Metropolis(int, long&, Matrix, double&, double&, Vector);
// prints to file the results of the calculations  
void output(int, int, double, Vector);

// Main program begins here

int main(int argc, char* argv[])
{
  char *outfilename;
  double w[17], InitialTemp, FinalTemp, Energy, MagneticMoment, TempStep;

  if (argc <= 1) {
    cout << "Bad Usage: " << argv[0] << 
      " read output file" << endl;
    exit(1);
  }
  if (argc > 1) {
    outfilename=argv[1];
    ofile.open(outfilename); 
    int NSpins = argv[2];
    int MCcycles = argv[3];    
  }
  double InitialTemp = 2.1;
  double FinalTemp = 2.3;
  double TempStep = 0.01;
  //  Allocate memory for spin matrix
  Matrix SpinMatrix(NSpins,Npsins);
        // Start Monte Carlo sampling by looping over T first
  for ( double temperature = InitialTemp; temperature <= FinalTemp; temperature+=TempStep){
    //    initialise energy and magnetization 
    double Energy = MagneticMoment = 0.;
    // initialise array for expectation values
    initialize(NSpins, SpinMatrix, Energy, MagneticMoment);
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    for( int i = 0; i < 5; i++) average[i] = 0.;
    for( int i = 0; i < 5; i++) total_average[i] = 0.;
    // start Monte Carlo computation
    for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
      Metropolis(n_spins, idum, spin_matrix, E, M, w);
      // update expectation values  for local node
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);
    }
        output(n_spins, mcs, temperature, total_average);
    }
  }
  free_matrix((void **) spin_matrix); // free memory
  ofile.close();  // close output file
  // End MPI
  MPI_Finalize (); 
  return 0;
}

// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, int **spin_matrix, 
		double& E, double& M)
{
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      spin_matrix[y][x] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[y][x];
    }
  }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
	(spin_matrix[periodic(y,n_spins,-1)][x] +
	 spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
{
  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran2(&idum)*(double)n_spins);
      int iy = (int) (ran2(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
	(spin_matrix[iy][periodic(ix,n_spins,-1)]+
	 spin_matrix[periodic(iy,n_spins,-1)][ix] +
	 spin_matrix[iy][periodic(ix,n_spins,1)] +
	 spin_matrix[periodic(iy,n_spins,1)][ix]);
      if ( ran2(&idum) <= w[deltaE+8] ) {
	spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
      }
    }
  }
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temperature, double *total_average)
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles 
  double Etotal_average = total_average[0]*norm;
  double E2total_average = total_average[1]*norm;
  double Mtotal_average = total_average[2]*norm;
  double M2total_average = total_average[3]*norm;
  double Mabstotal_average = total_average[4]*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2total_average- Etotal_average*Etotal_average)/n_spins/n_spins;
  double Mvariance = (M2total_average - Mtotal_average*Mtotal_average)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Etotal_average/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins << endl;
} // end output function

/*
** The function 
**         ran2()
** is a long periode (> 2 x 10^18) random number generator of 
** L'Ecuyer and Bays-Durham shuffle and added safeguards.
** Call with idum a negative integer to initialize; thereafter,
** do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int            j;
  long           k;
  static long    idum2 = 123456789;
  static long    iy=0;
  static long    iv[NTAB];
  double         temp;

  if(*idum <= 0) {
    if(-(*idum) < 1) *idum = 1;
    else             *idum = -(*idum);
    idum2 = (*idum);
    for(j = NTAB + 7; j >= 0; j--) {
      k     = (*idum)/IQ1;
      *idum = IA1*(*idum - k*IQ1) - k*IR1;
      if(*idum < 0) *idum +=  IM1;
      if(j < NTAB)  iv[j]  = *idum;
    }
    iy=iv[0];
  }
  k     = (*idum)/IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if(*idum < 0) *idum += IM1;
  k     = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j     = iy/NDIV;
  iy    = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1) iy += IMM1;
  if((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran2()


  /*
   * The function                             
   *      void  **matrix()                    
   * reserves dynamic memory for a two-dimensional matrix 
   * using the C++ command new . No initialization of the elements. 
   * Input data:                      
   *  int row      - number of  rows          
   *  int col      - number of columns        
   *  int num_bytes- number of bytes for each 
   *                 element                  
   * Returns a void  **pointer to the reserved memory location.                                
   */

void **matrix(int row, int col, int num_bytes)
  {
  int      i, num;
  char     **pointer, *ptr;

  pointer = new(nothrow) char* [row];
  if(!pointer) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for "<< row << "row addresses !" << endl;
    return NULL;
  }
  i = (row * col * num_bytes)/sizeof(char);
  pointer[0] = new(nothrow) char [i];
  if(!pointer[0]) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for address to " << i << " characters !" << endl;
    return NULL;
  }
  ptr = pointer[0];
  num = col * num_bytes;
  for(i = 0; i < row; i++, ptr += num )   {
    pointer[i] = ptr; 
  }

  return  (void **)pointer;

  } // end: function void **matrix()

    /*
     * The function                         
     *      void free_matrix()              
     * releases the memory reserved by the function matrix() 
     *for the two-dimensional matrix[][] 
     * Input data:                          
     *  void far **matr - pointer to the matrix
     */

void free_matrix(void **matr)
{

  delete [] (char *) matr[0];
  delete [] matr;

}  // End:  function free_matrix() 



