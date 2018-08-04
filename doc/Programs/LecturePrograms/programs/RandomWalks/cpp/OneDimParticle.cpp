

        /*
         The program simulates a process with a single particle in
	equilibrium with the surroindings at temperatur T.
        The model describes an ideal gass in one dimension and
	should reproduce the Boltzman distribution function
	using the Metropolis algorithm
	*/


#include <iostream>
#include <fstream>
#include <iomanip>
#include <new>
#include <cmath>
#include <stdio.h>
#include "lib.h"

using namespace std;

#define   ONE_LINE    100

           //   function declarations

void input_data(double *, int *, double *, double *, int *, char *);
void print_output(double, int, double, int, double, double, 
                                   double, int, int *, char *);

int main()
{
  char   
         file_name[ONE_LINE];
  int
         i, mc_cycles, num_intervals, accepted_steps, 
         *distribution_func;

  double
         kT, initial_v,  v_range, curr_v, v_step, test_v,
         beta, deltaE, energy, average_v, average_energy, average_energy2;
  long
         idum;

  for( ;  ; ) {

    input_data(&kT, &mc_cycles, &initial_v, &v_range, &num_intervals, file_name);

    distribution_func  = new(nothrow) int [num_intervals];    // number of velocity intervals
    if(!distribution_func) {
      printf("\n\nError in function main() :");
      printf("\n Not enough memory for the velcity intervals N = \n", num_intervals);
      return 1;
    }

    for(i = 0; i < num_intervals; i++) {               // initialization
      distribution_func[i] = 0;
    }

    idum             = -1;
    beta             = 1.0/kT;
    curr_v           = initial_v;
    average_v        = curr_v;
    energy           = 0.5 * curr_v * curr_v;
    average_energy   = energy;
    average_energy2  = energy * energy;
    accepted_steps  = 0;
    v_step          = 2.0 * v_range/((double) num_intervals);

    distribution_func[num_intervals/2 + (int) (curr_v/v_step)]++;      // storing initial point

    for(i = 0; i < mc_cycles; i++) {
      test_v   = curr_v + (2.0 * ran1(&idum) - 1.0) * v_range;
      deltaE = 0.5 * (test_v * test_v - curr_v * curr_v); 
      if(ran1(&idum) <= exp(-beta * deltaE)) {
	accepted_steps++;
        curr_v  = test_v;
        energy += deltaE;
      }
      distribution_func[num_intervals/2 + (int) (curr_v/v_step)]++; 
      average_v       += curr_v;
      average_energy  += energy;
      average_energy2 += energy * energy;
    
    } // end Monte Carl cycles
      
      // normalize to the total number of Monte Carlo cycles
      
    average_v       /= (mc_cycles + 1);
    average_energy  /= (mc_cycles + 1);
    average_energy2 /= (mc_cycles + 1);

    print_output(kT, mc_cycles, v_range, num_intervals, 
                average_v,average_energy, 
                 sqrt(average_energy2 - average_energy * average_energy),
                 num_intervals, distribution_func, file_name);

  } // process terminates

  return 0;

} // End: function main()


       /* 
       ** The function 
       **      input_data()
       ** reads all necessary starting data
       ** from standard input
       */

void input_data(double *kT, int *mc_cycles, double *initial_v, 
                double *v_range, int *num_intervals, char *file_name)
{
  double    value;

  *kT           = 0.0;                              // initialization
  *mc_cycles    = 0;
  *initial_v    = 0.0;
  value         = 0;
  *num_intervals = 0;


  cout << "\n\nSimulation of ideal gass system in one dimension:";
  cout <<"\nNumber of Monte Carlo cycles = ";
  cin >> *mc_cycles;

  if(*mc_cycles <= 0) exit(1);        // program terminates

  cout <<"\nTemperatur of the system: kT = ";
  cin >> *kT;
  cout <<"\nInitial velocity = ";
  cin >> *initial_v;
  cout <<"\nMaximum velocity range(value*sqrt(T)): value  = ";
  cin >> value;
  *v_range = value * sqrt(*kT);
  cout << "\nNumber of intervals (stops at N < 0) N = ";
  cin >> *num_intervals;

  cout << "\nName of data file for gnuplot = ";
  cin >> file_name;
  

} // End: function input_data()

     /* 
     ** The function 
     **        print_output()
     ** prints to standard output the cumulative
     ** distribution of a radioactive nuclear
     ** decay process
     */

void print_output(double kT, int mc_cycles, double v_range, int N, 
            double average_v, double average_energy, 
             double standard_deviation, int num_intervals,
	     int *distribution_func, char *file_name)
{
  int
             i;
  ofstream 
             ofile;


  cout << "\n\nSimulation of ideal gass system in one dimension";
  cout << "\n using the Metropolis algorithm:\n";
  cout << "\n       Result:";
  cout << "\nTemperatur of the system: kT = " <<  kT;
  cout << "\nNumber of Monte Carlo cycles = " <<  mc_cycles;
  cout << "\nMaximum velocity range(value*sqrt(T)) = " << v_range;
  cout << "\nNumber of intervals N = " << N;

  cout << "\n\n Average velocity = " << average_v;
  cout << "\n Average energy = " << average_energy; 
  cout << "\n Standard deviation =  " << standard_deviation;


  ofile.open(file_name);

  for(i = 0; i < N; i++) {
    ofile << setw(15) << setprecision(6) << (i - N/2);
    ofile << setw(15) << setprecision(6) <<distribution_func[i] << endl;
  }


} // End: function print_output()
