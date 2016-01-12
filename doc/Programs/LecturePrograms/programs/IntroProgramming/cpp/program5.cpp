// program to compute exp(-x) without factorials
using namespace std;
#include <iostream>
#include <cmath> 
#define  TRUNCATION     1.0E-10

int main()
{
  int       loop, n;
  double    x, term, sum;

  for(loop = 0; loop <= 100; loop += 10){
    x    = (double) loop;          // initialization 
    sum  = 1.0;
    term = 1;
    n    = 1;
    while(fabs(term) > TRUNCATION){
      term *= -x/((double) n);
      sum  += term;
      n++;
    } // end while loop 
    cout << " x =" << x << " exp = " << exp(-x) << "series = " << sum;
    cout  << " number of terms = " << n << endl;
  } // end of for loop 
}  //    End: function main() 
