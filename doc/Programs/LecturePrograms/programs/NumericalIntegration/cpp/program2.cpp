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

//THE INTEGRAND FUNCTION IN CARTESIAN COORDINATES
double func_cart(double x1, double y1, double z1, double x2, double y2, double z2){
	if  ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) != 0)
		return exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))) 
		          / sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	else 
		return 0;
}

//THE INTEGRAND FUNCTION IN POLAR COORDINATES
double func_polar(double r1, double t1, double p1, double r2, double t2, double p2){
	double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
	double f = exp(-4*(r1+r2))*r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
		return f;
	else 
		return 0;
}

//THE INTEGRAND FUNCTION IN POLAR COORDINATES REDUCED FOR GAUSSIAN LAGUERRE
double func_polar_lag(double r1, double t1, double p1, double r2, double t2, double p2){
	double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
	double f = exp(-3*(r1+r2))*r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
		return f;
	else 
		return 0;
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

       /*
       ** The function 
       **              gauleg()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

void gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359; 
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */
 
	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /* 
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()




//FUNCTIONS FOR COMPUTING THE LAGUERRE POLYNOMIALS WEIGHTS
double gammln( double xx){
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
void gaulag(double *x, double *w, int n, double alf){
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}

//Plain Gauss-Legendre
void Gauss_Legendre(int n, double a, double b, double  &integral){
	double *x = new double [n];
	double *w = new double [n];
	gauleg(a,b,x,w,n);
	
	double int_gauss = 0.0;
	int i,j,k,l,f,t;
	#pragma omp parallel for reduction(+:int_gauss)  private (i,j,k,l,f,t)
	for (i = 0;  i < n; i++){
	for (j = 0;  j < n; j++){
	for (k = 0;  k < n; k++){
	for (l = 0;  l < n; l++){
	for (f = 0;  f < n; f++){
	for (t = 0;  t < n; t++){
		int_gauss+=w[i]*w[j]*w[k]*w[l]*w[f]*w[t]*func_cart(x[i],x[j],x[k],x[l],x[f],x[t]);
	}}}}}}
	integral = int_gauss;
	delete [] x;
	delete [] w;
}

//Gauss-Legendre Gauss-Laguerre in polar coordinates
void Gauss_Laguerre(int n_lag, int n_leg, double  &integral){
	double *xlag = new double [n_lag + 1];
	double *wlag = new double [n_lag + 1];
	double *xt = new double [n_leg];
	double *wt = new double [n_leg];
	double *xp = new double [n_leg];
	double *wp = new double [n_leg];

	gaulag(xlag,wlag,n_lag,0.0);
	gauleg(0,M_PI,xt,wt,n_leg);
	gauleg(0,2*M_PI,xp,wp,n_leg);
	double int_gauss = 0.0;
	
	int i,j,k,l,f,t;
	#pragma omp parallel for reduction(+:int_gauss)  private (i,j,k,l,f,t)
	for (i = 1;  i <= n_lag;  i++){    //r1
	for (j = 0;  j <  n_leg;  j++){    //t1
	for (k = 0;  k <  n_leg;  k++){    //p1
	for (l = 1;  l <= n_lag;  l++){    //r2
	for (f = 0;  f <  n_leg;  f++){    //t2
	for (t = 0;  t <  n_leg;  t++){    //p2
		int_gauss += wlag[i]*wlag[l]*wt[j]*wp[k]*wt[f]*wp[t]*func_polar_lag(xlag[i],xt[j],xp[k],xlag[l],xt[f],xp[t]);
	}}}}}}
	integral = int_gauss;
	
	delete [] xt;
	delete [] wt;
	delete [] xp;
	delete [] wp;
	delete [] xlag;
	delete [] wlag;
}

//Monte Carlo with finite domain
void Brute_MonteCarlo(int n, double a, double b, double  &integral, double  &std){
	double * x = new double [n];
	double x1, x2, y1, y2, z1, z2, f;
	double mc = 0.0;
	double sigma = 0.0;
	int i;
	double jacob = pow((b-a),6);
	
	#pragma omp parallel for reduction(+:mc)  private (i, x1, x2, y1, y2, z1, z2, f)
	for (i = 0; i < n; i++){
		x1=ran()*(b-a)+a;
		x2=ran()*(b-a)+a;
		y1=ran()*(b-a)+a;
		y2=ran()*(b-a)+a;
		z1=ran()*(b-a)+a;
		z2=ran()*(b-a)+a;
		f=func_cart(x1, x2, y1, y2, z1, z2);
		mc += f;
		x[i] = f;
	}
	mc = mc/((double) n );
	#pragma omp parallel for reduction(+:sigma)  private (i)
	for (i = 0; i < n; i++){
		sigma += (x[i] - mc)*(x[i] - mc); 
	}
	sigma = sigma*jacob/((double) n );
	//	std = sqrt(sigma)/sqrt(n);
        std = sigma;
	integral = mc*jacob;
	delete [] x;
}

//Monte Carlo with polar coordinates
void Brute_Polar_MonteCarlo(int n, double a, double  &integral, double  &std){
	double * x = new double [n];
	double r1, r2, t1, t2, p1, p2, f;
	double mc = 0.0;
	double sigma = 0.0;
	double jacob = a*a*4*pow(M_PI,4);
	int i;
	
	#pragma omp parallel for reduction(+:mc)  private (i, r1, r2, t1, t2, p1, p2, f)
	for (i = 0; i < n; i++){
		r1=ran()*a;
		r2=ran()*a;
		t1=ran()*M_PI;
		t2=ran()*M_PI;
		p1=ran()*2*M_PI;
		p2=ran()*2*M_PI;
		f=func_polar(r1, t1, p1, r2, t2, p2);
		mc += f;
		x[i] = f;
	}
	mc = mc/((double) n );
	#pragma omp parallel for reduction(+:sigma)  private (i)
	for (i = 0; i < n; i++){
		sigma += (x[i] - mc)*(x[i] - mc); 
	}
	sigma = sigma*jacob/((double) n );
        std = sigma;
	//	std = sqrt(sigma)/sqrt(n);
	integral = mc*jacob;
	delete [] x;
}

//Monte Carlo with polar coordinates and change of variables
void Polar_MonteCarlo(int n, double  &integral, double  &std){
	double * x = new double [n];
	double r1, r2, t1, t2, p1, p2, f;
	double mc = 0.0;
	double sigma = 0.0;
	double jacob = 4*pow(M_PI,4)/16;
	double rr1,rr2;
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
		f=func_polar(r1, t1, p1, r2, t2, p2)/((1-rr1)*(1-rr2));
		mc += f;
		x[i] = f;
	}
	mc = mc/((double) n );
	#pragma omp parallel for reduction(+:sigma)  private (i)
	for (i = 0; i < n; i++){
		sigma += (x[i] - mc)*(x[i] - mc); 
	}
	sigma = sigma*jacob/((double) n );	
	//	std = sqrt(sigma)/sqrt(n);
	std = sigma;
	integral = mc*jacob;
	delete [] x;
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
	std = sigma;
	//	std = sqrt(sigma)/sqrt(n);
	integral = mc*jacob;
	delete [] x;
}


int main(){
    omp_set_num_threads(NUM_THREADS);  
	int n,n_mc;
	double a,b;
	
	printf("EXACT RESULT:\t%.8f\t\n", 5*M_PI*M_PI/256);
	
	n = 25;
	a = -3;
	b =  3;
	
	double int_leg;
	Gauss_Legendre(n,a,b,int_leg);
	printf("Gau-Leg:     \t%.8f\t\n", int_leg);
	
	int n_lag = 25;
	int n_leg = 25;
	double int_lag;
	Gauss_Laguerre(n_lag, n_leg, int_lag);
	printf("Gau-Lag:     \t%.8f\t\n", int_lag);
	
	n_mc = 1000000;
	a = -3;
	b =  3;
	srand(time(NULL));
	generator.seed(time(NULL));
	
	double brute_mc, brute_std;
	Brute_MonteCarlo(n_mc, a, b, brute_mc, brute_std);
	printf("BF MC:       \t%.8f\t%.8f\n", brute_mc, brute_std);
	
	double brute_polar_mc, brute_polar_std;
	Brute_Polar_MonteCarlo(n_mc, 4, brute_polar_mc, brute_polar_std);
	printf("BF Pol MC:   \t%.8f\t%.8f\n", brute_polar_mc, brute_polar_std);
	
	double polar_mc, polar_std;
	Polar_MonteCarlo(n_mc, polar_mc, polar_std);
	printf("Polar MC:    \t%.8f\t%.8f\n", polar_mc, polar_std);
	
	double polar_imp_mc, polar_imp_std;
	Polar_MonteCarlo_Importance(n_mc, polar_imp_mc, polar_imp_std);
	printf("Polar MC IS: \t%.8f\t%.8f\n", polar_imp_mc, polar_imp_std);
}



