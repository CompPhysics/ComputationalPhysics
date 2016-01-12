
// Radioactive decay of nuclei 
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace  std;

ofstream ofile;

// Function to read in data from screen  
void initialise(int&, int&, int&, double& ) ;
// The Mc sampling for nuclear decay 
void  mc_sampling(int, int, int, double, int*);
// prints to screen the results of the calculations 
void  output(int, int, int *);
int main(int argc, char* argv[])
{
  char *outfilename;
  int initial_n_particles, max_time, number_cycles; 
  double decay_probability;
  int *ncumulative;
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
  // Read in data 
  initialise(initial_n_particles, max_time, number_cycles, 
	     decay_probability) ;
  ncumulative = new int [max_time+1];
  // Do the mc sampling  
  mc_sampling(initial_n_particles, max_time, number_cycles, 
	      decay_probability, ncumulative);
  // Print out results 
  output(max_time, number_cycles, ncumulative);
  delete [] ncumulative;
  return 0; 
}

void initialise(int& initial_n_particles, int& max_time, int& number_cycles, double& decay_probability) 
{
  cout << "Initial number of particles = " << endl ;
  cin >> initial_n_particles;
  cout << "maximum time = " << endl;
  cin >> max_time;
  cout << "# MC steps= " << endl;
  cin >> number_cycles;
  cout << "# Decay probability= " << endl;
  cin >> decay_probability;
}  // end of function initialise   




void output(int max_time, int number_cycles, int* ncumulative)
{
  int i;
  for( i=0; i <= max_time; i++){
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << i;
    ofile << setw(15) << setprecision(8);
    ofile  << ncumulative[i]/((double) number_cycles) << endl;
  }
}  // end of function output 


void mc_sampling(int initial_n_particles, int max_time, int number_cycles, double decay_probability, int *ncumulative)
{
  int cycles, time, np, n_unstable, particle_limit;
  long idum;

  idum=-1;  // initialise random number generator
  // loop over monte carlo cycles 
  // One monte carlo loop is one sample  
  for (cycles = 1; cycles <= number_cycles; cycles++){   
    n_unstable = initial_n_particles;
    //  accumulate the number of particles per time step per trial 
    ncumulative[0] += initial_n_particles;
    // loop over each time step 
    for (time=1; time <= max_time; time++){
      // for each time step, we check each particle
      particle_limit = n_unstable;
      for ( np = 1; np <=  particle_limit; np++) {
        if( ran0(&idum) <= decay_probability) {
          n_unstable=n_unstable-1;
        }
      }  // end of loop over particles 
      ncumulative[time] += n_unstable;
    }  // end of loop over time steps 
  }    // end of loop over MC trials 
}   // end mc_sampling function  


