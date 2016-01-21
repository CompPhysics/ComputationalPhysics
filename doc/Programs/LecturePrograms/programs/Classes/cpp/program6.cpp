#include <iostream>
using namespace std;
//  Declare functions before main
void func(int, int*);
int main(int argc, char *argv[])
{
  int a;
  int *b;
  a = 10;
  b = new int[10];
  for(int i = 0; i < 10; i++) {
    b[i] = i;
    cout <<  b[i] << endl;
  }
  // the variable a is transferred by call by value. This means
  //  that the function func cannot change a in the calling function
  func( a,b);

  delete [] b ;
  return 0;
} // End: function main()

void func( int x, int *y)
{
  // a becomes locally x  and it can be changed locally
  x+=7;
  //  func gets the address of the first element of y (b)
   // it changes y[0] to 10 and when returning control to main
  // it changes also b[0]. Call by reference
  *y += 10;  //  *y = *y+10;
  //  explicit element
  y[6] += 10;
  //   in this function y[0]  and y[6] have been changed and when returning
  // control to main  this means that b[0] and b[6] are changed.
  return;
} // End: function func()
