#include <iostream>
//#include "solarsystem.h"
#include "planet.h"
#include <cmath>
#include <armadillo>


using namespace arma;
using namespace std;



int main()
{
  solarsystem mysystem;

planet  Sun(1,0,0,0,0,0,0);
planet  Mercury(1.2e-7, 0.39, 0, 0,0,9.96,0);
planet  Venus(2.4e-6, 0.72, 0, 0,0,7.36,0);
planet  Earth(1.5e-6,1,0,0, 0, 6.26, 0);
planet  Mars(3.3e-7, 1.52, 0, 0,0,5.06,0);
planet  Jupiter(9.5e-4, 5.20, 0,0,0,2.75,0);
planet  Saturn(2.75e-4, 9.54, 0, 0,0,2.04,0);
planet  Uranus(4.4e-5, 19.19, 0, 0,0,1.43,0);
planet  Neptune(5.1e-5, 30.06, 0, 0,0,1.14,0);
planet  Pluto(5.6e-9, 39.53, 0, 0,0,0.99,0);


   mysystem.add(Sun);
   mysystem.add(Mercury);
   mysystem.add(Venus);
   mysystem.add(Earth);
   mysystem.add(Mars);
   mysystem.add(Jupiter);
   mysystem.add(Saturn);
   mysystem.add(Uranus);
   mysystem.add(Neptune);
   mysystem.add(Pluto);

    int elements = mysystem.number_planets;
    cout << "number of element" << elements<< endl;

    /*
   if(RK4) {
    mysystem.solverRK4(mysystem.all_planets, 0.001, 100 );
    cout << "RK4";
   }else{
     mysystem.solverVERLET(mysystem.all_planets, 0.001, 100 );
    cout << "Vertel";
   }
*/
}

