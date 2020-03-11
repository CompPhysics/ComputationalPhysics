#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string>   
#include <sstream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <random>


#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-8
#define NUM_THREADS 8


using namespace std;

minstd_rand0 generator;

inline double ran(){
	//return ((double) generator())/2147483647;
	return ((double) rand()) / RAND_MAX;
}



//THE INTEGRAND FUNCTION IN POLAR COORDINATES FOR THE IMPORTANCE SAMPLING MONTE CARLO
double func_polar_mc(double r1, double t1, double p1, double r2, double t2, double p2){
	double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
	double f = r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
		return f;
	else 
		return 0;
}





//Monte Carlo with polar coordinates and importance sampling
void Polar_MonteCarlo_Importance(int n, double  &integral, double  &std){
	double * x = new double [n];
	double r1, r2, t1, t2, p1, p2, f,rr1,rr2;
	double mc = 0.0;
	double sigma = 0.0;
	double jacob = 4*pow(M_PI,4)/16;
	int i;
	
	#pragma omp parallel for reduction(+:mc)  private (i, r1, r2, t1, t2, p1, p2, rr1, rr2, f)
	for (i = 0; i < n; i++){
		rr1=ran();
		rr2=ran();
		r1=-0.25*log(1-rr1);
		r2=-0.25*log(1-rr2);
		t1=ran()*M_PI;
		t2=ran()*M_PI;
		p1=ran()*2*M_PI;
		p2=ran()*2*M_PI;
		f=func_polar_mc(r1, t1, p1, r2, t2, p2);
		mc += f;
		x[i] = f;
	}
	mc = mc/((double) n);
	#pragma omp parallel for reduction(+:sigma)  private (i)
	for (i = 0; i < n; i++){
		sigma += (x[i] - mc)*(x[i] - mc); 
	}
	sigma = sigma*jacob/((double) n );
	std = sqrt(sigma)/sqrt(n);
	integral = mc*jacob;
	delete [] x;
}


int main(){
    omp_set_num_threads(NUM_THREADS);  
	int n_mc;
	double a,b;
	
	printf("EXACT RESULT:\t%.8f\t\n", 5*M_PI*M_PI/256);
	
	n_mc = 1000000;
	srand(time(NULL));
	generator.seed(time(NULL));
	
	double polar_imp_mc, polar_imp_std;
	Polar_MonteCarlo_Importance(n_mc, polar_imp_mc, polar_imp_std);
	printf("Polar MC IS: \t%.8f\t%.8f\n", polar_imp_mc, polar_imp_std);
}



