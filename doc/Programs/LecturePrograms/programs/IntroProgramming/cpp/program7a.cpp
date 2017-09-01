
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;
int main()  
{  
  int var;                    
  int *pointer;  
   
  pointer = &var;  
  var  = 421;  
  cout << "Address of the integer variable var :" << &var << endl ;
  cout << "Value of var" << var << endl;
  cout << "Value of the integer pointer variable:" << pointer << endl;
  cout << "Value which pointer is pointing at :" << *pointer << endl;
  cout<< "Address of the pointer variable :" << &pointer << endl;
  return 0;
}
