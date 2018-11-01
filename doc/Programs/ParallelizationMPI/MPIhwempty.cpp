//    First C++ example of MPI Hello world
using namespace std;
#include <mpi.h>
#include <iostream>
  
int main (int nargs, char* args[])
{
  int numprocs, my_rank;
  //   MPI initializations
  // FILL in MPI statements


  cout << "Hello world, I have  rank " << my_rank << " out of " << numprocs << endl;
  //  End MPI

  return 0;
}
