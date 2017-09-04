#include <iostream>
#include <cmath>
//#include <cstdio>
using namespace std;
int main (int argc, char* argv[])
{
  double r, s;        // declare variables
  r = atof(argv[1]);  // convert the text argv[1] to double 
  s = sin(r);
  cout << "Hello, World! sin(" << r << ") =" << s << endl;
  return 0;           // success execution of the program
}
