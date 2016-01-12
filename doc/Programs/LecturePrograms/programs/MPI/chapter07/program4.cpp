//    Third C++ example of MPI Hello world
using namespace std;
#include <mpi.h>
#include <iostream>
  
int main (int nargs, char* args[])
{
     int numprocs, my_rank, flag;
//   MPI initializations
     MPI_Status status;
     MPI_Init (&nargs, &args);
     MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
     MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
     //   Send and Receive example
     if (my_rank > 0)
       MPI_Recv (&flag, 1, MPI_INT, my_rank-1, 100, MPI_COMM_WORLD, &status);
       cout << "Hello world, I have  rank " << my_rank << " out of " << numprocs << endl;
     if (my_rank < numprocs-1)
	 MPI_Send (&my_rank, 1, MPI_INT, my_rank+1, 100, MPI_COMM_WORLD);
//  End MPI
      MPI_Finalize ();
    return 0;
 }
