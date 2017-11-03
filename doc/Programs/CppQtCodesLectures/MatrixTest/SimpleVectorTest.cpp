#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

int main()
{
  vec b = randu<vec>(5);
  b.print();
  vec c(5);
  c.fill(3);
  c = 3;  // c(0) =3
  (b.subvec(1,3)).print();
  return 0;
}
