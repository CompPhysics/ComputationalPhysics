/*
  This program computes the entropy for random walks
  based on a computed probability histogram.
  Many walkers (input) make several trials steps with
  a given number of walks per trial. Periodic boundary conditions
  have been implemented since the walkers are forced to move on a
  one-dimensional lattice of size -L  <= x <= L. This means that if the walker
  goes beyond L(-L) it is shifted to -L(L).
*/
#include <iostream>
#include <iomanip>
#include "lib.h"
using namespace  std;

// Function to read in data from screen, note call by reference
void initialise(int&, int&, int&, double&) ;
// The Mc sampling for random walks
void  mc_sampling(int, int, int, double, int *, int *);
// prints to file the results of the calculations 
void  output(int, int,  int, int *);

int main()
{
  int walkers, time_steps, length; 
  double move_probability;
  // Read in data 
  initialise(walkers, time_steps, length, move_probability);
  // allocate memory
  int *x = new int [walkers+1];
  int *probability = new int [2*(length+1)];
  // set all arrays equal zero
  for (int i = 1; i <= walkers; i++){   
    x[i] = 0;
  }
  for (int j = 0; j <= 2*length; j++){  
    probability[j] = 0;
  } // end initialization of vectors

  // Do the mc sampling  
  mc_sampling(walkers, time_steps, length, move_probability, x, probability);
  // Print out results 
  output(length, time_steps, walkers, probability);
  delete [] x; 
  delete [] probability;
  return 0; 
} // end main function


//  This is the sampling part
void mc_sampling(int walkers, int time_steps, int length, 
                 double move_probability, int *x, int *probability)
{
  long idum;
  idum=-1;  // initialise random number generator
  // loop over all time steps
  for (int step=1; step <= time_steps; step++){
    // move all walkers with periodic boundary conditions
    for (int walks = 1; walks <= walkers; walks++){   
      if (ran0(&idum) <= move_probability) {
	if ( x[walks] +1 > length) {
	  x[walks] = -length;
	}
	else{
	  x[walks] += 1; 
	}
      } 
      else {
	if ( x[walks] -1 < -length) {
	  x[walks] = length;
        }
	else{
	  x[walks] -= 1; 
	}
      }
    }  // end of loop over walks
  } // end of loop over trials
  // at the final time step we compute the probability
  // by counting the number of walkers at every position
  for ( int i = -length; i <= length; i++){
    int count = 0; 
    for( int j = 1; j <= walkers; j++){
      if ( x[j] == i ) {
        count += 1;
      }
    }
    probability[i+length] = count;
  }
}
// Writes the results to screen
void output(int length, int time_steps, int walkers, int *probability)
{
  double entropy, histogram;
  // find norm of probability
  double norm = 1.0/walkers;
  // compute the entropy
  entropy = 0.; histogram = 0.;
  for( int  i = -length; i <=  length; i++){
    histogram = (double) probability[i+length]*norm;
    if ( histogram > 0.0) {
    entropy -= histogram*log(histogram);
    }
  }
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setw(6) << time_steps;
    cout << setw(15) << setprecision(8) << entropy << endl;
}  // end of function output 

// Reads in data from screen
void initialise(int& walkers, int& time_steps, int& length, double& move_probability) 
{
  cout << "Number of walkers ="; 
  cin >> walkers;
  cout << "Number of time steps=";
  cin >> time_steps;
  cout << "Length of one dimensional lattice=";
  cin >> length;
  cout << "Move probability=";
  cin >> move_probability;
}  // end of function initialise   






