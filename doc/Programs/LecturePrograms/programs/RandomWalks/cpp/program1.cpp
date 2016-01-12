/*
  1-dim random walk program. 
  A walker makes several trials steps with
  a given number of walks per trial
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace  std;

// Function to read in data from screen, note call by reference
void initialise(int&, int&, double&) ;
// The Mc sampling for random walks
void  mc_sampling(int, int, double, int *, int *, int *);
// prints to screen the results of the calculations 
void  output(int, int, int *, int *, int *);

int main()
{
  int max_trials, number_walks; 
  double move_probability;
  // Read in data 
  initialise(max_trials, number_walks, move_probability) ;
  int *walk_cumulative = new int [number_walks+1];
  int *walk2_cumulative = new int [number_walks+1];
  int *probability = new int [2*(number_walks+1)];
  for (int walks = 1; walks <= number_walks; walks++){   
    walk_cumulative[walks] = walk2_cumulative[walks] = 0;
  }
  for (int walks = 0; walks <= 2*number_walks; walks++){   
    probability[walks] = 0;
  } // end initialization of vectors
  // Do the mc sampling  
  mc_sampling(max_trials, number_walks, move_probability, 
              walk_cumulative, walk2_cumulative, probability);
  // Print out results 
  output(max_trials, number_walks, walk_cumulative, 
         walk2_cumulative, probability);
  delete [] walk_cumulative; // free memory
  delete [] walk2_cumulative; delete [] probability;
  return 0; 
} // end main function

// note call by reference, differs from C++
void initialise(int& max_trials, int& number_walks, double& move_probability) 
{
  cout << "Number of Monte Carlo trials ="; 
  cin >> max_trials;
  cout << "Number of attempted walks=";
  cin >> number_walks;
  cout << "Move probability=";
  cin >> move_probability;
}  // end of function initialise   


void output(int max_trials, int number_walks, 
            int *walk_cumulative, int *walk2_cumulative, int * probability)
{
  ofstream ofile("testwalkers.dat");
  ofstream probfile("probability.dat");
  for( int  i = 1; i <=  number_walks; i++){
    double xaverage = walk_cumulative[i]/((double) max_trials);
    double x2average = walk2_cumulative[i]/((double) max_trials);
    double variance = x2average - xaverage*xaverage;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(6) << i;
    ofile << setw(15) << setprecision(8) << xaverage;
    ofile << setw(15) << setprecision(8) << variance << endl;
  }
  ofile.close();
  // find norm of probability
  double norm = 0.;
  for( int  i = -number_walks; i <=  number_walks; i++){
    norm += (double) probability[i+number_walks];
  }
  // write probability
  for( int  i = -number_walks; i <=  number_walks; i++){
    double histogram = probability[i+number_walks]/norm;
    probfile << setiosflags(ios::showpoint | ios::uppercase);
    probfile << setw(6) << i;
    probfile << setw(15) << setprecision(8) << histogram << endl;
  }
  probfile.close();

}  // end of function output 


void mc_sampling(int max_trials, int number_walks, 
                 double move_probability, int *walk_cumulative, 
                 int *walk2_cumulative, int *probability)
{
  long idum;
  idum=-1;  // initialise random number generator
  for (int trial=1; trial <= max_trials; trial++){
    int position = 0;
    for (int walks = 1; walks <= number_walks; walks++){   
      if (ran0(&idum) <= move_probability) {
	position += 1;
      } 
      else {
	position -= 1;
      }
      walk_cumulative[walks] += position;
      walk2_cumulative[walks] += position*position;
      probability[position+number_walks] += 1;
    }  // end of loop over walks
  } // end of loop over trials
}   // end mc_sampling function  
