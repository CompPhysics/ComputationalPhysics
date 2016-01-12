#include <stdio.h>
#include <iostream.h>
#include <math.h>
#include <fstream.h>

/*

Different methods for solving ODEs are presented
We are solving the following eqation:

m*l*(phi)'' + reib*(phi)' + m*g*sin(phi) = A*cos(omega*t)

If you want to solve similar equations with other values you have to
rewrite the methods 'derivatives' and 'initialise' and change the variables in the private
part of the class Pendel

At first we rewrite the equation using the following definitions:

omega_0 = sqrt(g*l)
t_roof = omega_0*t
omega_roof = omega/omega_0
Q = (m*g)/(omega_0*reib)
A_roof = A/(m*g)

and we get a dimensionless equation

(phi)'' + 1/Q*(phi)' + sin(phi) = A_roof*cos(omega_roof*t_roof)

This equation can be written as two equations of first order:

(phi)' = v
(v)' = -v/Q - sin(phi) +A_roof*cos(omega_roof*t_roof)

All numerical methods are applied to the last two equations.
The algorithms are taken from the book "An introduction to computer simulation methods"
*/

class pendelum
 {
 private:
   double Q, A_roof, omega_0, omega_roof,g; //
   double y[2];          //for the initial-values of phi and v
   int n;                // how many steps
   double delta_t,delta_t_roof;

 public:
   void derivatives(double,double*,double*);
   void initialise();
   void euler();
   void euler_cromer();
   void midpoint();
   void euler_richardson();
   void half_step();
   void rk2(); //runge-kutta-second-order
   void rk4_step(double,double*,double*,double); // we need it in function rk4() and asc()
   void rk4(); //runge-kutta-fourth-order
   void asc(); //runge-kutta-fourth-order with adaptive stepsize control
 };

void pendelum::derivatives(double t, double* in, double* out)
{ /* Here we are calculating the derivatives at (dimensionless) time t
     'in' are the values of phi and v, which are used for the calculation
     The results are given to 'out' */
  
  out[0]=in[1];             //out[0] = (phi)'  = v
  if(Q)
    out[1]=-in[1]/((double)Q)-sin(in[0])+A_roof*cos(omega_roof*t);  //out[1] = (phi)''
  else
    out[1]=-sin(in[0])+A_roof*cos(omega_roof*t);  //out[1] = (phi)''
}


void pendelum::initialise()
{
  double m,l,omega,A,viscosity,phi_0,v_0,t_end;
  cout<<"Solving the differential eqation of the pendulum!\n";
  cout<<"We have a pendulum with mass m, length l. Then we have a periodic force with amplitude A and omega\n";
  cout<<"Furthermore there is a viscous drag coefficient.\n";
  cout<<"The initial conditions at t=0 are phi_0 and v_0\n";
  cout<<"Mass m: ";
  cin>>m;
  cout<<"length l: ";
  cin>>l;
  cout<<"omega of the force: ";
  cin>>omega;
  cout<<"amplitude of the force: ";
  cin>>A;
  cout<<"The value of the viscous drag constant (viscosity): ";
  cin>>viscosity;
  cout<<"phi_0: ";
  cin>>y[0];
  cout<<"v_0: ";
  cin>>y[1];
  cout<<"Number of time steps or integration steps:";
  cin>>n;
  cout<<"Final time steps as multiplum of pi:";
  cin>>t_end;
  t_end *= acos(-1.);
  g=9.81;
  // We need the following values:
  omega_0=sqrt(g/((double)l));      // omega of the pendulum
  if (viscosity)  Q= m*g/((double)omega_0*viscosity);
  else Q=0; //calculating Q
  A_roof=A/((double)m*g);
  omega_roof=omega/((double)omega_0);
  delta_t_roof=omega_0*t_end/((double)n);    //delta_t without dimension
  delta_t=t_end/((double)n);
}

void pendelum::euler()
{ //using simple euler-method
  int i;
  double yout[2],y_h[2];
  double t_h;

  y_h[0]=y[0];
  y_h[1]=y[1];
  t_h=0;
  ofstream fout("euler.out");
  fout.setf(ios::scientific);
  fout.precision(20);
  for(i=1;i<=n;i++){
    derivatives(t_h,y_h,yout);
    yout[1]=y_h[1]+yout[1]*delta_t_roof;
    yout[0]=y_h[0]+yout[0]*delta_t_roof;
    // Calculation with dimensionless values    
    fout<<i*delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
    t_h+=delta_t_roof;
    y_h[1]=yout[1];
    y_h[0]=yout[0];
  }
  fout.close;
}

void pendelum::euler_cromer()
{
  int i;
  double t_h;
  double yout[2],y_h[2];

  t_h=0;
  y_h[0]=y[0];  //phi
  y_h[1]=y[1];  //v
  ofstream fout("ec.out");
  fout.setf(ios::scientific);
  fout.precision(20);
  for(i=1; i<=n; i++){
    derivatives(t_h,y_h,yout);
    yout[1]=y_h[1]+yout[1]*delta_t_roof;
    yout[0]=y_h[0]+yout[1]*delta_t_roof;
    // The new calculated value of v is used for calculating phi 
    fout<<i*delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
    t_h+=delta_t_roof;
    y_h[0]=yout[0];
    y_h[1]=yout[1];
  }
  fout.close;
}

void pendelum::midpoint()
{
  int i;
  double t_h;
  double yout[2],y_h[2];
  
  t_h=0;
  y_h[0]=y[0];  //phi
  y_h[1]=y[1];  //v
  ofstream fout("midpoint.out");
  fout.setf(ios::scientific);
  fout.precision(20);
  for(i=1; i<=n; i++){
    derivatives(t_h,y_h,yout);
    yout[1]=y_h[1]+yout[1]*delta_t_roof;
    yout[0]=y_h[0]+0.5*(yout[1]+y_h[1])*delta_t_roof;
    fout<<i*delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
    t_h+=delta_t_roof;
    y_h[0]=yout[0];
    y_h[1]=yout[1];
  }
  fout.close;
}


void pendelum::euler_richardson()
{
  int i;
  double t_h,t_m;
  double yout[2],y_h[2],y_m[2];

  t_h=0;
  y_h[0]=y[0];  //phi
  y_h[1]=y[1];  //v
  ofstream fout("er.out");
  fout.setf(ios::scientific);
  fout.precision(20);
  for(i=1; i<=n; i++){
    derivatives(t_h,y_h,yout);
    y_m[1]=y_h[1]+0.5*yout[1]*delta_t_roof;
    y_m[0]=y_h[0]+0.5*y_h[1]*delta_t_roof;
    t_m=t_h+0.5*delta_t_roof;
    derivatives(t_m,y_m,yout);
    yout[1]=y_h[1]+yout[1]*delta_t_roof;
    yout[0]=y_h[0]+y_m[1]*delta_t_roof;
    fout<<i*delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
    t_h+=delta_t_roof;
    y_h[0]=yout[0];
    y_h[1]=yout[1];
  }
  fout.close;
}

void pendelum::half_step()
{
  /*We are using the half_step_algorith.
    The algorithm is not self-starting, so we calculate
    v_1/2 by using the Euler algorithm. */

  int i;
  double t_h;
  double yout[2],y_h[2];

  t_h=0;
  y_h[0]=y[0];  //phi
  y_h[1]=y[1];  //v
  ofstream fout("half_step.out");
  fout.setf(ios::scientific);
  fout.precision(20);
  /*At first we have to calculate v_1/2
    For this we use Euler's method:
    v_`1/2 = v_0 + 1/2*a_0*delta_t_roof
    For calculating a_0 we have to start derivatives
  */
  derivatives(t_h,y_h,yout);
  yout[1]=y_h[1]+0.5*yout[1]*delta_t_roof;
  yout[0]=y_h[0]+yout[1]*delta_t_roof;
  fout<<delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
  y_h[0]=yout[0];
  y_h[1]=yout[1];
  for(i=2; i<=n; i++){
    derivatives(t_h,y_h,yout);
    yout[1]=y_h[1]+yout[1]*delta_t_roof;
    yout[0]=y_h[0]+yout[1]*delta_t_roof;
    fout<<i*delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
    t_h+=delta_t_roof;
    y_h[0]=yout[0];
    y_h[1]=yout[1];
  }
  fout.close;
}

void pendelum::rk2()
{
  /*We are using the second-order-Runge-Kutta-algorithm
    We have to calculate the parameters k1 and k2 for v and phi,
    so we use to arrays k1[2] and k2[2] for this
    k1[0], k2[0] are the parameters for phi,
    k1[1], k2[1] are the parameters for v
  */

  int i;
  double t_h;
  double yout[2],y_h[2],k1[2],k2[2],y_k[2];
  
  t_h=0;
  y_h[0]=y[0];  //phi
  y_h[1]=y[1];  //v
  ofstream fout("rk2.out");
  fout.setf(ios::scientific);
  fout.precision(20);
  for(i=1; i<=n; i++){
    /*Calculation of k1 */
    derivatives(t_h,y_h,yout);
    k1[1]=yout[1]*delta_t_roof;
    k1[0]=yout[0]*delta_t_roof;
    y_k[0]=y_h[0]+k1[0]*0.5;
    y_k[1]=y_h[1]+k2[1]*0.5;
    /*Calculation of k2 */
    derivatives(t_h+delta_t_roof*0.5,y_k,yout);
    k2[1]=yout[1]*delta_t_roof;
    k2[0]=yout[0]*delta_t_roof;
    yout[1]=y_h[1]+k2[1];
    yout[0]=y_h[0]+k2[0];
    fout<<i*delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
    t_h+=delta_t_roof;
    y_h[0]=yout[0];
    y_h[1]=yout[1];
  }
  fout.close;
}

void pendelum::rk4_step(double t,double *yin,double *yout,double delta_t)
{
  /*
    The function calculates one step of fourth-order-runge-kutta-method
    We will need it for the normal fourth-order-Runge-Kutta-method and
    for RK-method with adaptive stepsize control

    The function calculates the value of y(t + delta_t) using fourth-order-RK-method
    Input: time t and the stepsize delta_t, yin (values of phi and v at time t)
    Output: yout (values of phi and v at time t+delta_t)
    
  */
  double k1[2],k2[2],k3[2],k4[2],y_k[2];
  // Calculation of k1 
  derivatives(t,yin,yout);
  k1[1]=yout[1]*delta_t;
  k1[0]=yout[0]*delta_t;
  y_k[0]=yin[0]+k1[0]*0.5;
  y_k[1]=yin[1]+k1[1]*0.5;
  /*Calculation of k2 */
  derivatives(t+delta_t*0.5,y_k,yout);
  k2[1]=yout[1]*delta_t;
  k2[0]=yout[0]*delta_t;
  y_k[0]=yin[0]+k2[0]*0.5;
  y_k[1]=yin[1]+k2[1]*0.5;
  /* Calculation of k3 */
  derivatives(t+delta_t*0.5,y_k,yout);
  k3[1]=yout[1]*delta_t;
  k3[0]=yout[0]*delta_t;
  y_k[0]=yin[0]+k3[0];
  y_k[1]=yin[1]+k3[1];
  /*Calculation of k4 */
  derivatives(t+delta_t,y_k,yout);
  k4[1]=yout[1]*delta_t;
  k4[0]=yout[0]*delta_t;
  /*Calculation of new values of phi and v */
  yout[0]=yin[0]+1.0/6.0*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
  yout[1]=yin[1]+1.0/6.0*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
}

void pendelum::rk4()
{
  /*We are using the fourth-order-Runge-Kutta-algorithm
    We have to calculate the parameters k1, k2, k3, k4 for v and phi,
    so we use to arrays k1[2] and k2[2] for this
    k1[0], k2[0] are the parameters for phi,
    k1[1], k2[1] are the parameters for v
  */

  int i;
  double t_h;
  double yout[2],y_h[2]; //k1[2],k2[2],k3[2],k4[2],y_k[2];
  
  t_h=0;
  y_h[0]=y[0];  //phi
  y_h[1]=y[1];  //v
  ofstream fout("rk4.out");
  fout.setf(ios::scientific);
  fout.precision(20);
  for(i=1; i<=n; i++){
    rk4_step(t_h,y_h,yout,delta_t_roof);
    fout<<i*delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
    t_h+=delta_t_roof;
    y_h[0]=yout[0];
    y_h[1]=yout[1];
  }
  fout.close;
}



void pendelum::asc()
{
  /*
    We are using the Runge-Kutta-algorithm with adaptive stepsize control
    according to "Numerical Recipes in C", S. 574 ff.
    
    At first we calculate y(x+h) using rk4-method  => y1
    Then we calculate y(x+h) using two times rk4-method at x+h/2 and x+h  => y2
    
    The difference between these values is called "delta" If it is smaller than a given value,
    we calculate y(x+h) by  y2 + (delta)/15 (page 575, Numerical R.)
    
    If delta is not smaller than ... we calculate a new stepsize using
    h_new=(Safety)*h_old*(.../delta)^(0.25) where "Safety" is constant (page 577 N.R.)
    and start again with calculating y(x+h)...
   */
  int i;
  double t_h,h_alt,h_neu,hh,errmax;
  double yout[2],y_h[2],y_m[2],y1[2],y2[2], delta[2], yscal[2];
  
  const double eps=1.0e-6;
  const double safety=0.9;
  const double errcon=6.0e-4;
  const double tiny=1.0e-30;

  t_h=0;
  y_h[0]=y[0];  //phi
  y_h[1]=y[1];  //v
  h_neu=delta_t_roof;
  ofstream fout("asc.out");
  fout.setf(ios::scientific);
  fout.precision(20);
  for(i=0;i<=n;i++){
    /* The error is scaled against yscal
       We use a yscal of the form yscal = fabs(y[i]) + fabs(h*derivatives[i])
       (N.R. page 567)
    */
    derivatives(t_h,y_h,yout);
    yscal[0]=fabs(y[0])+fabs(h_neu*yout[0])+tiny;
    yscal[1]=fabs(y[1])+fabs(h_neu*yout[1])+tiny;
    /* the do-while-loop is used until the */
    do{
      /* Calculating y2 by two half steps */
      h_alt=h_neu;
      hh=h_alt*0.5;
      rk4_step(t_h, y_h, y_m, hh);
      rk4_step(t_h+hh,y_m,y2,hh);
      /* Calculating y1 by one normal step */
      rk4_step(t_h,y_h,y1,h_alt);
      /* Now we have two values for phi and v at the time t_h + h  in y2 and y1
	 We can now calculate the delta for phi and v
      */
      delta[0]=fabs(y1[0]-y2[0]);
      delta[1]=fabs(y1[1]-y2[1]);
      errmax=(delta[0]/yscal[0] > delta[1]/yscal[1] ? delta[0]/yscal[0] : delta[1]/yscal[1]);
      
      /*We scale delta against the constant yscal
	Then we take the biggest one and call it errmax */
      errmax=(double)errmax/eps;
      /*We divide errmax by eps and have only   */
      h_neu=safety*h_alt*exp(-0.25*log(errmax));
    }while(errmax>1.0);
    /*Now we are outside the do-while-loop and have a delta which is small enough
      So we can calculate the new values of phi and v
    */
    yout[0]=y_h[0]+delta[0]/15.0;
    yout[1]=y_h[1]+delta[1]/15.0;
    fout<<(double)(t_h+h_alt)/omega_0<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";
    // Calculating of the new stepsize
    h_neu=(errmax > errcon ? safety*h_alt*exp(-0.20*log(errmax)) : 4.0*h_alt);
    y_h[0]=yout[0];
    y_h[1]=yout[1];
    t_h+=h_neu;
  }
}


int main()
{
  pendelum testcase;
  testcase.initialise();
  testcase.euler();
  testcase.euler_cromer();
  testcase.midpoint();
  testcase.euler_richardson();
  testcase.half_step();
  testcase.rk2();
  testcase.rk4();
  return 0;
}  // end of main function
