/*
**     Program to compute the second derivative of exp(x).
**     In this version we use C++ options for reading and
**     writing files and data. The rest of the code is as in
**     programs/chapter3/program1.cpp 
**     Three calling functions are included
**     in this version. In one function we read in the data from screen,
**     the next function computes the second derivative
**     while the last function prints out data to screen.
*/
using namespace std;
# include  <iostream> 
# include <fstream>
# include <iomanip>
# include <cmath>
void initialise (double *, double *, int *);
void second_derivative( int, double, double, double *, double *);
void output( double *, double *, double, int);

ofstream ofile;

int main(int argc, char* argv[])
{
  // declarations of variables 
  char *outfilename;
  int number_of_steps;
  double x, initial_step;
  double *h_step, *computed_derivative;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  ofile.open(outfilename);
  //  read in input data from screen 
  initialise (&initial_step, &x, &number_of_steps);
  //  allocate space in memory for the one-dimensional arrays  
  //  h_step and computed_derivative                           
  h_step =  new double[number_of_steps];
  computed_derivative = new double[number_of_steps];
  //  compute the second derivative of exp(x) 
  second_derivative( number_of_steps, x, initial_step, h_step, 
		     computed_derivative);        
  //  Then we print the results to file  
  output(h_step, computed_derivative, x, number_of_steps );
  // free memory
  delete [] h_step;
  delete [] computed_derivative; 
  // close output file
  ofile.close();
  return 0;
}   // end main program 

//     Read in from screen the initial step, the number of steps 
//     and the value of x 

void initialise (double *initial_step,  double *x, int *number_of_steps)
{
  printf("Read in from screen initial step, x and number of steps\n");
  scanf("%lf %lf %d",initial_step, x, number_of_steps);
  return;
}  // end of function initialise 

//  This function computes the second derivative 

void second_derivative( int number_of_steps, double x, 
                        double initial_step, double *h_step, 
                        double *computed_derivative)
{
  int counter;
  double h;
  //     calculate the step size  
  //     initialise the derivative, y and x (in minutes) 
  //     and iteration counter 
  h = initial_step;
  //  start computing for different step sizes 
  for (counter=0; counter < number_of_steps; counter++ ) {
      //  setup arrays with derivatives and step sizes
      h_step[counter] = h;
      computed_derivative[counter] = 
	(exp(x+h)-2.*exp(x)+exp(x-h))/(h*h);
      h = h*0.5;
    } // end of do loop 
  return;
}   // end of function second derivative 


//    function to write out the final results  
void output(double *h_step, double *computed_derivative, double x, int number_of_steps )
{
  int i;
  double epsilon;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for( i=0; i < number_of_steps; i++){
      epsilon = log10(fabs(computed_derivative[i]-exp(x))/exp(x));      
      ofile << setw(15) << setprecision(8) << log10(h_step[i]);
      ofile << setw(15) << setprecision(8) << epsilon << endl;
}
}  // end of function output

