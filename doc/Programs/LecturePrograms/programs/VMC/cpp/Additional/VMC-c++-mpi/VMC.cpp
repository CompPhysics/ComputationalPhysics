//VMC-implementation for the Helium atom in pure C++ with MPI
//Written by Magnar K. Bugge
#include <mpi.h>
#include <iostream>
#include <math.h>
#include "lib.cpp"
#include <fstream>


using namespace std;

//Sjekker om vi får singularitetsproblemer med R
bool harSingularitet(double R[2][3]) {

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

//Lokal energi
double E_local(double R[2][3], double alpha) {

  //Absoluttverdier av r:
  double r1 = sqrt(R[0][0]*R[0][0]+R[0][1]*R[0][1]+R[0][2]*R[0][2]);
  double r2 = sqrt(R[1][0]*R[1][0]+R[1][1]*R[1][1]+R[1][2]*R[1][2]);
  double r12 = sqrt((R[0][0]-R[1][0])*(R[0][0]-R[1][0]) +
		    (R[0][1]-R[1][1])*(R[0][1]-R[1][1]) +
		    (R[0][2]-R[1][2])*(R[0][2]-R[1][2]));

  return (alpha - 2)*(1/r1 + 1/r2) + 1/r12 - alpha*alpha;
}

//Proeve-boelgefunksjon
double Psi_trial(double R[2][3], double alpha) {

  //Absoluttverdier av r:
  double r1 = sqrt(R[0][0]*R[0][0]+R[0][1]*R[0][1]+R[0][2]*R[0][2]);
  double r2 = sqrt(R[1][0]*R[1][0]+R[1][1]*R[1][1]+R[1][2]*R[1][2]);

  return exp(-alpha*(r1+r2));
}

//Funksjonen for MC-beregning av <H>
void kjoerMC(double& sum, double& kvadratsum, int antMC, 
	     double delta, long& idum, double alpha, int& N, int& accepted) {

  //Initialiser posisjon
  double R[2][3]; //2 partikler, 3 dimensjoner
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      R[i][j] = 0.2;
    }
  }

  double EL; //Lokal energi

  //For Metropolis algo:
  double R_trial[2][3]; //Trial-posisjon
  double P_trial; //Sannsynlighetstetthet i trial-posisjon
  double P; //Sannsynlighetstetthet i nåværende posisjon
  double r; //Tilfeldig tall for Metropolis algo
  double w; //Forhold mellom sannsynligheter

  //Løkke over MC-sykluser
  for (int k = 0; k < antMC; k++) {

    //Beregn R_trial (prøveposisjon for Metropolis algo)
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
	R_trial[i][j] = R[i][j] + delta * (ran0(&idum)*2 - 1);
      }
    }

    P_trial = Psi_trial(R_trial, alpha); //Beregner sannsynlighetstetthet
    P_trial = P_trial * P_trial; //i R_trial ...
    P = Psi_trial(R, alpha);
    P = P * P; //... og i R

    if (P_trial > P) {

      //Aksepter hvis større sannstetthet
      for (int l = 0; l < 2; l++) {
	for (int m = 0; m < 3; m++) {
	  R[l][m] = R_trial[l][m]; 
	}
      }
      accepted++;
    } 
    else { //Hvis ikke, må vi bruke tilfeldig tall
      
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
    
    if (!harSingularitet(R) && k > antMC/20) { //Beregn bidrag hvis vi ikke har singularitet + termalisering

      EL = E_local(R, alpha);
      
      sum += EL;
      kvadratsum += EL * EL;

      N++;
    }
  }  
}

//Funksjon som skal være nær 0 for optimal delta (50% aksepterte moves)
double difference(double sum, double kvadratsum, int antMC2, double delta, long idum, double alpha, int N, int accepted) {

   N = accepted = 0;
   kjoerMC(sum, kvadratsum, antMC2, delta, idum, alpha, N, accepted);
   return (double)accepted/antMC2 - .5;
}

int main() {
  
  //MPI-initialisering
  int numprocs;
  int minrank;
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &minrank);

  int antMC = 10000000; //Antall MC-sykluser
  int antMC2 = 10000; //Antall MC-sykluser for bestemmelse av delta
  double E; //Forventningsverdien for energien
  double E2; //Forventningsverdien for kvadratet av energien
  double sigma; //Standardavvik
  double alpha; //Variasjonsparameter
  double delta; //Lengde av Metropolis-hopp
  int N; //Antall sykluser uten singulartiteter
  int accepted = 0; //Antall aksepterte moves
  double sum = 0; //Sum i MC-beregning
  double kvadratsum = 0; //Kvadratsum i MC-beregning
  long idum = time(NULL);
  double delta_max = 2.0;//Størrelser for beregning av beste
  double delta_min = .01; //delta-verdi
  double minimum, maximum;
  double tolerance = .01;
  double error;

  if (idum == 0) {cout << "FEIL";}

  ofstream utfil("data");
  
  //Allokering av lokale størrelser
  double lokal_sum;
  double lokal_kvadratsum;
  int lokal_N;
  int lokal_accepted;
  int lokal_antMC = antMC / numprocs;

  //Løkke over alpha-verdier
  for (alpha = 1.4; alpha <= 2.51; alpha += .05) {

    //Her finnes en god verdi av delta ved at nullpunktet til difference
    //finnes ved halveringsmetoden
    minimum = delta_min;
    maximum = delta_max;
    while ((maximum - minimum) > tolerance) {

      if (difference(sum, kvadratsum, antMC2, minimum, idum, alpha, N, accepted) * 
	  difference(sum, kvadratsum, antMC2, (minimum+maximum)/2, idum, alpha, N, accepted) < 0)
	maximum = (minimum + maximum) / 2;
      else
	minimum = (minimum + maximum) / 2;
	 
    }
    delta = (minimum + maximum) / 2;
    
    //Initialiser lokale verdier
    lokal_sum = 0;
    lokal_kvadratsum = 0;
    lokal_N = 0;
    lokal_accepted = 0;

    kjoerMC(lokal_sum, lokal_kvadratsum, lokal_antMC, delta, idum, alpha, lokal_N, lokal_accepted); //Utfør MC-beregning
    
    //Summer opp fra alle maskiner
    MPI_Reduce(&lokal_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&lokal_kvadratsum, &kvadratsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&lokal_accepted, &accepted, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&lokal_N, &N, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    //Beregn størrelsene vi er ute etter
    E = sum / N;
    E2 = kvadratsum / N;
    sigma = sqrt(E2 - E*E);
    error = sigma / sqrt(N);

    //Skriv ut verdier inkl. andelen av aksepterte punkter
    if (minrank == 0) {
      utfil << alpha << " " << E << " " << sigma << " " << error << " " << (double)accepted/antMC << endl;
    }
  }
 
  utfil.close();
  
  //Avslutt MPI
  MPI_Finalize();

  return 0;
}
