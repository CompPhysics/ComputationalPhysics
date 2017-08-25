#include <iostream>
using namespace std;
//  Declare functions before main
int main(int argc, char *argv[])
{
  int a =112, b = -1;
  double c = 3.14;
  int *d = &a;
  double  *e = &c;

  cout << "Address of the integer variable a :" << &a <<endl;
  cout << "Value of the integer pointer variable d:" << d << endl;
  cout << "Value which pointer d is pointing at :" << *d << endl;

  cout << "Address of the double variable c :" << &c <<endl;
  cout << "Value of the integer pointer variable e:" << e << endl;
  cout << "Value which pointer e is pointing at :" << *e << endl;


  return 0;
} // End: function main()
