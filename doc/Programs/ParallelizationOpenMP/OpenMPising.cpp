/* 
   Program to solve the two-dimensional Ising model 
   with zero external field and no parallelization
   Parallel version using MPI
   The coupling constant J is set to J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis aolgorithm  is used as well as periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
   Run as
   ./executable Outputfile numberof spins number of MC cycles initial temp final temp tempstep
   ./test.x Lattice 100 10000000 2.1 2.4 0.01
   Compile and link as 
   c++ -O3 -std=c++11 -Rpass=loop-vectorize -o Ising.x <codename> -lomp
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
using namespace  std;
using namespace arma;
// output file
ofstream ofile;

// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void InitializeLattice(int, mat &, double&, double&);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &);
// prints to file the results of the calculations  
void WriteResultstoFile(int, int, double, vec);

// Main program begins here

int main(int argc, char* argv[])
{
  string filename;
  int NSpins, MonteCarloCycles;
  double InitialTemp, FinalTemp, TempStep;
  int NProcesses, RankProcess;
  if (argc <= 5) {
    cout << "Bad Usage: " << argv[0] << 
      " read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
    exit(1);
  }
    filename=argv[1];
    NSpins = atoi(argv[2]);
    MonteCarloCycles = atoi(argv[3]);    
    InitialTemp = atof(argv[4]);
    FinalTemp = atof(argv[5]);
    TempStep = atof(argv[6]);
    string fileout = filename;
    string argument = to_string(NSpins);
    fileout.append(argument);
    ofile.open(fileout);
  // Start Monte Carlo sampling by looping over the selected Temperatures
  for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
    vec LocalExpectationValues = zeros<mat>(5);
    // Start Monte Carlo computation and get local expectation values
    MetropolisSampling(NSpins, MonteCarloCycles, Temperature, LocalExpectationValues);
    // Find total average
    vec TotalExpectationValues = zeros<mat>(5);
    WriteResultstoFile(NSpins, MonteCarloCycles*NProcesses, Temperature, TotalExpectationValues);
  }  
  ofile.close();  // close output file
  return 0;
}


// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, int MonteCarloCycles, double Temperature, vec &ExpectationValues)
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // Initialize the lattice spin values
  mat SpinMatrix = zeros<mat>(NSpins,NSpins);
  //    initialize energy and magnetization 
  double Energy = 0.;     double MagneticMoment = 0.;
  // initialize array for expectation values
  InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);
  // setup array for possible energy changes
  vec EnergyDifference = zeros<mat>(17); 
  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  // Start Monte Carlo experiments
  int AllSpins = NSpins*NSpins;
  for (int cycles = 1; cycles <= MonteCarloCycles; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(int Spins =0; Spins < AllSpins; Spins++) {
      int ix = (int) (RandomNumberGenerator(gen)*NSpins);
      int iy = (int) (RandomNumberGenerator(gen)*NSpins);
      int deltaE =  2*SpinMatrix(ix,iy)*
	(SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
	 SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
	 SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
	 SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
      if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
	SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
	MagneticMoment += 2.0*SpinMatrix(ix,iy);
	Energy += (double) deltaE;
      }
    }
    // update expectation values  for local node after a sweep through the lattice
    ExpectationValues(0) += Energy;    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;    
    ExpectationValues(3) += MagneticMoment*MagneticMoment; 
    ExpectationValues(4) += fabs(MagneticMoment);
  }
} // end of Metropolis sampling over spins

// function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment)
{
  // setup spin matrix and initial magnetization using cold start, all spins pointing up or down
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      SpinMatrix(x,y) = 1.0; // spin orientation for the ground state
      MagneticMoment +=  (double) SpinMatrix(x,y);
    }
  }
  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
	(SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
	 SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialize



void WriteResultstoFile(int NSpins, int MonteCarloCycles, double temperature, vec ExpectationValues)
{
  double norm = 1.0/((double) (MonteCarloCycles));  // divided by  number of cycles 
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double AllSpins = 1.0/((double) NSpins*NSpins);
  double HeatCapacity = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)*AllSpins/temperature/temperature;
  double MagneticSusceptibility = (M2_ExpectationValues - M_ExpectationValues*M_ExpectationValues)*AllSpins/temperature;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << E_ExpectationValues*AllSpins;
  ofile << setw(15) << setprecision(8) << HeatCapacity;
  ofile << setw(15) << setprecision(8) << M_ExpectationValues*AllSpins;
  ofile << setw(15) << setprecision(8) << MagneticSusceptibility;
  ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues*AllSpins << endl;
} // end output function



    




