#include <iostream>
#include <iomanip>
#include "ising2dim_backend.h"
using namespace std;

int main(int argc, char* argv) {
  double E_av, E_variance, M_av, M_variance, Mabs_av;
  E_av = E_variance = M_av = M_variance = Mabs_av = 0.0;
  double temp = 2.5;

  monteCarlo(temp,5,50000, E_av, E_variance, M_av, M_variance, Mabs_av);
  
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << "temp:       " << setw(15) << setprecision(8) << temp       << endl;
  cout << "E_av:       " << setw(15) << setprecision(8) << E_av       << endl;
  cout << "E_variance: " << setw(15) << setprecision(8) << E_variance << endl;
  cout << "M_av:       " << setw(15) << setprecision(8) << M_av       << endl;
  cout << "M_variance: " << setw(15) << setprecision(8) << M_variance << endl;
  cout << "Mabs_av:    " << setw(15) << setprecision(8) << Mabs_av    << endl;
}
