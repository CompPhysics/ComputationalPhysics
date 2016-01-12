#ifndef ising2dim_backend
#define ising2dim_backend
void monteCarlo(double temp, int size, int trials,
	       double& E_av, double& E_variance,
	       double& M_av, double& M_variance,
	       double& Mabs_av);
#endif
