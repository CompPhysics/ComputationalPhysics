#include <math.h>
#include "lib.cpp"

class data { //An object of this class is returned to Python, and contains all relevant variables

public:
  double sum;
  double squaresum;
  int N;
  int accepted;
  long idum;
};

long seed() { //Returns a number used to initialize random number generator

  return time(NULL);
}

//Check for singularity
bool hasSingularity(double R[2][3]) {

  double r1 = sqrt(R[0][0]*R[0][0]+R[0][1]*R[0][1]+R[0][2]*R[0][2]);
  double r2 = sqrt(R[1][0]*R[1][0]+R[1][1]*R[1][1]+R[1][2]*R[1][2]);
  double r12 = sqrt((R[0][0]-R[1][0])*(R[0][0]-R[1][0]) +
		    (R[0][1]-R[1][1])*(R[0][1]-R[1][1]) +
		    (R[0][2]-R[1][2])*(R[0][2]-R[1][2]));

  if (r1 < 1e-10 || r2 < 1e-10 || r12 < 1e-10) {
    return true;
  }
  else {
    return false;
  }
}

//Local energy
double E_local(double R[2][3], double alpha) {

  //Absolute values
  double r1 = sqrt(R[0][0]*R[0][0]+R[0][1]*R[0][1]+R[0][2]*R[0][2]);
  double r2 = sqrt(R[1][0]*R[1][0]+R[1][1]*R[1][1]+R[1][2]*R[1][2]);
  double r12 = sqrt((R[0][0]-R[1][0])*(R[0][0]-R[1][0]) +
		    (R[0][1]-R[1][1])*(R[0][1]-R[1][1]) +
		    (R[0][2]-R[1][2])*(R[0][2]-R[1][2]));

  return (alpha - 2)*(1/r1 + 1/r2) + 1/r12 - alpha*alpha;
}

//Trial wavefunction
double Psi_trial(double R[2][3], double alpha) {

  //Absolute values
  double r1 = sqrt(R[0][0]*R[0][0]+R[0][1]*R[0][1]+R[0][2]*R[0][2]);
  double r2 = sqrt(R[1][0]*R[1][0]+R[1][1]*R[1][1]+R[1][2]*R[1][2]);

  return exp(-alpha*(r1+r2));
}

//The function which does the MC calculation of <H> (energy exp. value)
data runMC(int MCcycles, double delta, long idum, double alpha) {

  double sum = 0;
  double squaresum = 0;
  int N = 0; //Number of cycles that "count" (cycles without singularities)
  int accepted = 0; //Number of accepted moves

  //Initialize position
  double R[2][3]; //2 particles, 3 dimensions
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      R[i][j] = .5 * (ran0(&idum)*2 - 1);
    }
  }

  double EL; //Local energy

  //For Metropolis algo:
  double R_trial[2][3]; //Trial posistion
  double P_trial; //Probability density at trial position
  double P; //Probability density at current position
  double r; //Random number for Metropolis algo
  double w; //Ratio of probabilities

  //Loop over MC-cycles
  for (int k = 0; k < MCcycles; k++) {

    //Find R_trial (trial position for Metropolis algo)
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
	R_trial[i][j] = R[i][j] + delta * (ran0(&idum)*2 - 1);
      }
    }

    P_trial = Psi_trial(R_trial, alpha); //Find probability densities ...
    P_trial = P_trial * P_trial; //... at R_trial ...
    P = Psi_trial(R, alpha);
    P = P * P; //... and at R

    if (P_trial > P) {

      //Accept if larger probability density
      for (int l = 0; l < 2; l++) {
	for (int m = 0; m < 3; m++) {
	  R[l][m] = R_trial[l][m]; 
	}
      }
      accepted++;
    } 
    else { //If not, we use a random number
      
      r = ran0(&idum);
      w = P_trial / P;
      
      if (r < w) {
	for (int l = 0; l < 2; l++) {
	  for (int m = 0; m < 3; m++) {
	    R[l][m] = R_trial[l][m]; 
	  }
	}
	accepted++;
      }
    }
    
    if (!hasSingularity(R) && k > MCcycles/20) { //Calculate contribution unless singularity/thermalization

      EL = E_local(R, alpha);
      
      sum += EL;
      squaresum += EL * EL;

      N++;
    }
  }

  //Return sum, squaresum, N, accepted and idum to Python as a data object
  data x;
  x.sum = sum;
  x.squaresum = squaresum;
  x.N = N;
  x.accepted = accepted;
  x.idum = idum;
  return x;
}

