//    Trapezoidal rule and numerical integration using OpenMP

# include <cstdlib>
# include <iostream>
# include <ctime>
#include <omp.h>
using namespace std;

//     Here we define various functions called by the main program
double intfunction(double );
void timestamp ( );

//   Main function begins here
int main (int nargs, char* args[])
{
  cout << "  C++/OpenMP version" << endl;
  cout << "  Compute the integral of pi using the Trapezoidal rule" << endl;
  int thread_num = omp_get_max_threads ( );
  cout << "  The number of processors available = " << omp_get_num_procs ( ) << endl;
  cout << "  The number of threads available    = " << thread_num <<  endl;
  // defining time and various variables needed for the integration
  double wtime = omp_get_wtime ( );
  double a = 0.0 ; double b = 1.0;   int n = 10000;
  double h = (b-a)/n;    // h is the same for all processes 
  double fa = intfunction(a)/2.0;
  double fb = intfunction(b)/2.0;
  int j; 
  // Now we implement in parallel the trapezoidal rule
  double Trapezsum=0.;
  // Begin of parallel region 
# pragma omp parallel for default(shared) private (j) reduction(+:Trapezsum)
  for (j= 1; j <= n-1; j++){
    double x=j*h+a;
    Trapezsum+=intfunction(x);
  }
  Trapezsum=(Trapezsum+fb+fa)*h;
  wtime = omp_get_wtime ( ) - wtime;
  cout << "  Elapsed time in seconds = " << wtime << endl;
  cout << "Trapezoidal rule = " << Trapezsum << endl;
  timestamp ( );
  return 0;
}  // end of main program

//  this function defines the function to integrate
double intfunction(double x)
{
  double value = 4./(1.+x*x);
  return value;
} // end of function to evaluate


void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
