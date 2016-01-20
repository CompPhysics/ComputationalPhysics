#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std; // note use of namespace                      \


// output file as global variable

ofstream ofile;

// Begin of main program

int main(int argc, char* argv[])
{
  char *outfilename;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 2 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file and number of elements on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }

//  opening a file for the program
 ofile.open(outfilename);
 int i = atoi(argv[2]);
 //  int *a;
 //a = new int[i];
 double *a = new double[i];
 cout << " bytes for i=" << sizeof(i) << endl;
 for (int j = 0; j < i; j++) {
   a[j] = j*exp(2.0);
   // ofile instead of cout
   ofile << setw(15) << setprecision(8) << "a=" << a[j] << endl;
 }
 delete [] a; // free memory
 ofile.close();  // close output file
 return 0;
}
