#include <iostream>
#include <cmath>
using namespace std; // note use of namespace                      \

int main (int argc, char* argv[])
{
  int i = atoi(argv[1]);
  //  Dynamic memory allocation: need tp declare -a- as a pointer
  //  You can use double *a = new double[i];  or
  double *a;
  a = new double[i];
  //  double a[10];
  // the first of element of a, a[0], and its address is the
  // value of the pointer.
  /* This is a longer comment
     if we want a static memory allocation
     this is the way to do it
  */
  cout << " bytes for i=" << sizeof(i) << endl;
  for (int j = 0; j < i; j++) {
    a[j] = j*exp(2.0);
    cout << "a=" << a[j] << endl;
  }
  // freeing memory
  delete [] a;
  // to check for memory leaks, use the software called -valgrind-
  return 0;           /* success execution of the program */
}
