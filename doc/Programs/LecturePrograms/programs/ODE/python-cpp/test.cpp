#include "pendelum.h"

int main()
{
  pendelum testcase;
  testcase.euler();
  testcase.euler_cromer();
  testcase.midpoint();
  testcase.euler_richardson();
  testcase.half_step();
  testcase.rk2();
  testcase.rk4();
  testcase.asc();
  return 0;
}
