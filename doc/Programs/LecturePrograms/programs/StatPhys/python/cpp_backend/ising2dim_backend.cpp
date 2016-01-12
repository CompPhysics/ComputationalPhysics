#include <cmath>
#include "lib.h"

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}

//Practically identical to the python version
void monteCarlo(double temp, int size, int trials,
		double& E_av, double& E_variance,
		double& M_av, double& M_variance,
		double& Mabs_av) {

  //Setup and initialize spin matrix to ground state
  int** spin_matrix = new int*[size];
  for (int y = 0; y < size; y++) {
    spin_matrix[y] = new int[size];
    for (int x = 0; x < size; x++) {
      spin_matrix[y][x] = 1; //Note C-style double-array (y,x)
    }
  }

  //Setup (fixed) array for possible energy changes
  double w[17];
  for (int de = -8; de <= 8; de++   ) w[de+8] = 0.0;
  for (int de = -8; de <= 8; de += 4) w[de+8] = exp(-de/temp);
  
  //Calculate initial magnetization
  double M = 0.0; //State of magnetization at a given time
  for (int y = 0; y < size; y++) {
    for (int x = 0; x < size; x++) {
      M += spin_matrix[y][x];
    }
  }

  //Calculate intial energy
  double E = 0.0;
  for (int y = 0; y < size; y++) {
    for (int x = 0; x < size; x++) {
      E -= spin_matrix[y][x]*
	(spin_matrix[periodic(y,size,-1)][x] +
	 spin_matrix[y][periodic(x,size,-1)]);
    }
  }
  //Initialize random number generator
  long idum = -1;

  //Initialize other variables
  double E2_av, M2_av;
  E2_av = M2_av = 0.0;

  //Start Metropolis Montecarlo computation
  for (int i = 0; i < trials; i ++) {
    //Metropolis
    //Loop over all spins, pick a random spin each time
    for (int s = 0; s < size*size; s++) {
      int x = (int) (ran1(&idum)*size);
      int y = (int) (ran1(&idum)*size);
      int deltaE =  2*spin_matrix[y][x]*
	(spin_matrix[y][periodic(x,size,-1)]+
	 spin_matrix[periodic(y,size,-1)][x] +
	 spin_matrix[y][periodic(x,size,1)] +
	 spin_matrix[periodic(y,size,1)][x]);
      if ( ran1(&idum) <= w[deltaE+8] ) {
	spin_matrix[y][x] *= -1;  // flip one spin and accept new spin config
        M += (double) 2*spin_matrix[y][x];
        E += (double) deltaE;
      }
    }
    //Update expectation values
    E_av    += E;
    E2_av   += E*E;
    M_av    += M;
    M2_av   += M*M;
    Mabs_av += fabs(M);
  }

  //Normalize average values
  E_av       /= (double) trials;
  E2_av      /= (double) trials;
  M_av       /= (double) trials;
  M2_av      /= (double) trials;
  Mabs_av    /= (double) trials;
  //Calculate variance and normalize to per-point and temp
  E_variance  = (E2_av-E_av*E_av)/(double)(size*size*temp*temp);
  M_variance  = (M2_av-M_av*M_av)/(double)(size*size*temp);
  //Normalize returned averages to per-point
  E_av       /= (double)size*size;
  M_av       /= (double)size*size;
  Mabs_av    /= (double)size*size;
}
