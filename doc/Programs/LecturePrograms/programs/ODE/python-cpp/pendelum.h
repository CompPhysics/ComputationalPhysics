#ifndef PENDELUM_H
#define PENDELUM_H

class pendelum
 {
 private:
   double Q, A_roof, omega_0, omega_roof,g; //
   double y[2];          //for the initial-values of phi and v
   int n;                // how many steps
   double delta_t,delta_t_roof;

 public:
   pendelum(); //Initialize from keyboard
   pendelum(double m, double l, double omega, double A,
	    double viscosity, double phi_0, double v_0,
	    int N, double t_end);
   void derivatives(double,double*,double*);
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

#endif
