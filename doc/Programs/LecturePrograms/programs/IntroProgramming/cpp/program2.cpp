#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

int main (int argc, char* argv[])
{
  int i; 
  int terms[32]; // storage of a0, a1, etc, up to 32 bits
  int number = atoi(argv[1]); 
  // initialise the term a0, a1 etc
  for (i=0; i < 32 ; i++){ terms[i] = 0;}
  for (i=0; i < 32 ; i++){ 
    terms[i] = number%2;
    number /= 2;
  }
  // write out results
  cout << "Number of bytes used= " << sizeof(number) << endl;
  for (i=0; i < 32 ; i++){ 
    cout << " Term nr: " << i << "Value= " << terms[i];
    cout << endl;
  }
  return 0;  
}
