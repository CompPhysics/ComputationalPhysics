#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;
int main()  
{   
  int matr[2];  
  int *pointer;  
  pointer = &matr[0];  
  matr[0] = 321;  
  matr[1] = 322;  
  printf("\nAddress of the matrix element matr[1]: %p",&matr[0]);
  printf("\nValue of the matrix element  matr[1]; %d",matr[0]);
  printf("\nAddress of the matrix element matr[2]: %p",&matr[1]);
  printf("\nValue of the matrix element  matr[2]: %d\n", matr[1]);
  printf("\nValue of the pointer : %p",pointer);
  printf("\nValue which pointer points at  : %d",*pointer);
  printf("\nValue which  (pointer+1) points at: %d\n",*(pointer+1));
  printf("\nAddress of the pointer variable: %p\n",&pointer);
  return 0;
}
