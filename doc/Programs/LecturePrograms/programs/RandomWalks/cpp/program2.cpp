
/*
  Two dimensional random walk program. 
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
void  mc_sampling(int, int, double, int *, int *, int *, int *);
// prints to screen the results of the calculations 
void  output(int, int, int *, int *, int *, int *);

int main()
{
  int max_trials, number_walks; 
  double move_probability;
  // Read in data 
  initialise(max_trials, number_walks, move_probability) ;
  int *x_cum = new int [number_walks+1];
  int *x2_cum = new int [number_walks+1];
  int *y_cum = new int [number_walks+1];
  int *y2_cum = new int [number_walks+1];
  for (int walks = 1; walks <= number_walks; walks++){   
    x_cum[walks] = x2_cum[walks] = y_cum[walks] = y2_cum[walks] = 0;
  } // end initialization of vectors
  // Do the mc sampling  
  mc_sampling(max_trials, number_walks, move_probability, 
              x_cum, x2_cum, y_cum, y2_cum);
  // Print out results 
  output(max_trials, number_walks, x_cum, x2_cum, y_cum, y2_cum); 
  delete [] x_cum; // free memory
  delete [] x2_cum; delete [] y_cum; delete [] y2_cum;
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
            int *x_cum, int *x2_cum, int *y_cum, int *y2_cum)
{
  ofstream ofile("twodim_walkers.dat");
  for( int  i = 1; i <=  number_walks; i++){
    double xaverage = x_cum[i]/((double) max_trials);
    double x2average = x2_cum[i]/((double) max_trials);
    double xvariance = x2average - xaverage*xaverage;
    double yaverage = y_cum[i]/((double) max_trials);
    double y2average = y2_cum[i]/((double) max_trials);
    double yvariance = y2average - yaverage*yaverage;
    double total_average = xaverage+yaverage;
    double total_variance = xvariance+yvariance;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(6) << i;
    ofile << setw(15) << setprecision(8) << xaverage;
    ofile << setw(15) << setprecision(8) << xvariance; 
    ofile << setw(15) << setprecision(8) << yaverage;
    ofile << setw(15) << setprecision(8) << yvariance; 
    ofile << setw(15) << setprecision(8) << total_average;
    ofile << setw(15) << setprecision(8) << total_variance << endl;
  }
  ofile.close();

}  // end of function output 


void mc_sampling(int max_trials, int number_walks, 
                 double move_probability, int *x_cum, 
                 int *x2_cum, int *y_cum, int *y2_cum)
{
  long idum;
  idum=-1;  // initialise random number generator
  for (int trial=1; trial <= max_trials; trial++){
    int x = 0; int y = 0;
    for (int walks = 1; walks <= number_walks; walks++){   
      double rantest = ran0(&idum);
      if (rantest <= move_probability) {
	x += 1;
      } 
      else if (rantest <= 2*move_probability){
	x -= 1;
      }
      else if (rantest <= 3*move_probability){
        y += 1;
      }
      else {
        y -= 1;
      }
      x_cum[walks] += x;
      x2_cum[walks] += x*x;
      y_cum[walks] += y;
      y2_cum[walks] += y*y;
    }  // end of loop over walks
  } // end of loop over trials
}   // end mc_sampling function  




