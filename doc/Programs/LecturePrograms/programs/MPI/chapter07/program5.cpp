//    Reactangle rule and numerical integration usign MPI, example 1, send and Receive
using namespace std;
#include <mpi.h>
#include <iostream>
  
int main (int nargs, char* args[])
{
  int numprocs, my_rank, i, n = 1000; 
     double local_sum, rectangle_sum, x, h;
     //   MPI initializations
     MPI_Init (&nargs, &args);
     MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
     MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
     //   Read from screen a possible new vaue of n
     if (my_rank == 0 && nargs > 1) {
        n = atoi(args[1]); 
     }
     h = 1.0/n;
     //  Broadcast n and h to all processes
     MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast (&h, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     //  Every process sets up its contribution to the integral
     local_sum = 0.;
     for (i = my_rank; i < n; i += numprocs) {  
       x = (i+0.5)*h;
       local_sum += 4.0/(1.0+x*x);
     }
     local_sum *= h;
     if (my_rank == 0) {
       MPI_Status status;  
       rectangle_sum = local_sum;  
       for (i=1; i < numprocs; i++) {
	 MPI_Recv(&local_sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,500,MPI_COMM_WORLD,&status);  
	 rectangle_sum += local_sum;
       }
       cout << "Result: " << rectangle_sum  << endl;  
     }  else
       MPI_Send(&local_sum,1,MPI_DOUBLE,0,500,MPI_COMM_WORLD);
     // End MPI
     MPI_Finalize ();  
     return 0;
 }
