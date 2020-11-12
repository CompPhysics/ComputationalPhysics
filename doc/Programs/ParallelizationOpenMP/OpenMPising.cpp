/* 
   Program to solve the two-dimensional Ising model 
   with zero external field using MPI
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
*/
#include "omp.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#define NUM_THREADS 4
using namespace  std;

// output file
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void initialize(int, int **, double&, double&);
// The metropolis algorithm 
void Metropolis(int, int **, double&, double&, double *);
// prints to file the results of the calculations  
void output(int, int, double, double, double, double, double, double);
//  Matrix memory allocation
//  allocate space for a matrix
void  **matrix(int, int, int);
//  free space for  a matrix
void free_matrix(void **);

minstd_rand0 generator;

inline double ran(){
        //return ((double) generator())/2147483647;
        return ((double) rand()) / RAND_MAX;
}

// Main program begins here

int main(int argc, char* argv[])
{
  char *outfilename;

  int **spin_matrix, n_spins, mcs;
  double w[17], initial_temp, final_temp, E, M, temp_step;

  if (argc <= 1) {
    cout << "Bad Usage: " << argv[0] << 
      " read output file" << endl;
    exit(1);
  }
  if (argc > 1) {
    outfilename=argv[1];
    ofile.open(outfilename); 
  }
  n_spins = 2; mcs = 100000;  initial_temp = 1.0; final_temp = 2.4; temp_step =0.1;
  cout << "  C++/OpenMP version" << endl;
  cout << "  Ising model with OpenMP" << endl;
  omp_set_num_threads(NUM_THREADS);
    //  int thread_num = omp_get_max_threads ( );
  cout << "  The number of processors available = " << omp_get_num_procs ( ) << endl;
  // defining time and various variables needed for the integration
  double wtime = omp_get_wtime ( );
  //  Allocate memory for spin matrix
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  // every node has its own seed for the random numbers, this is important else
  // if one starts with the same seed, one ends with the same random numbers

  srand(time(NULL));
  generator.seed(time(NULL));
  // Start Monte Carlo sampling by looping over T first
  for ( double temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){
    //    initialise energy and magnetization 

    // initialise array for expectation values
    initialize(n_spins, spin_matrix, E, M);
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    int cycles; 
    // start Monte Carlo computation and parallel region
    double totalE, totalM, Mabs, totalE2, totalM2;    
    totalE = 0.0; totalM = Mabs = totalE2 = totalM2 = 0.0;
# pragma omp parallel for default(shared) private (cycles, E, M, spin_matrix) reduction(+:totalE,totalM,totalE2,totalM2,Mabs)
    E = M = 0.;
    for (cycles = 1; cycles <= mcs; cycles++){
      Metropolis(n_spins, spin_matrix, E, M, w);
      // update expectation values
      totalE += E; totalE2 += E*E;
      totalM += M;    totalM2 += M*M; Mabs += fabs(M);
    }
    output(n_spins, mcs, temperature, totalE, totalM, Mabs, totalE2, totalM2);
  }
  free_matrix((void **) spin_matrix); // free memory
  wtime = omp_get_wtime ( ) - wtime;
  cout << "  Elapsed time in seconds = " << wtime << endl;
  ofile.close();  // close output file
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

void Metropolis(int n_spins, int **spin_matrix, double& E, double&M, double *w)
{
  long idum;
  idum = -1;  // random starting point
  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran()*(double)n_spins);
      int iy = (int) (ran()*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
	(spin_matrix[iy][periodic(ix,n_spins,-1)]+
	 spin_matrix[periodic(iy,n_spins,-1)][ix] +
	 spin_matrix[iy][periodic(ix,n_spins,1)] +
	 spin_matrix[periodic(iy,n_spins,1)][ix]);
      if ( ran() <= w[deltaE+8] ) {
	spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
      }
    }
  }
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temperature, double totalE, double totalM, double Mabs, double totalE2, double totalM2)
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles 
  double Etotal_average = totalE*norm;
  double E2total_average = totalE2*norm;
  double Mtotal_average = totalM*norm;
  double M2total_average = totalM2*norm;
  double Mabstotal_average = Mabs*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2total_average- Etotal_average*Etotal_average)/n_spins/n_spins;
  double Mvariance = (M2total_average - Mabstotal_average*Mabstotal_average)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Etotal_average/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins << endl;
} // end output function


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





