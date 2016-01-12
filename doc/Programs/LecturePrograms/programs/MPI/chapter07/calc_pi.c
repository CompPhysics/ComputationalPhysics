#include "mpi.h"
#include <stdio.h>

int main (int nargs, char** args)
{
  int size, my_rank, i, n = 1000;
  double l_sum, g_sum, x, h;

  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  if (my_rank==0 && nargs>1)
    n = atoi(args[1]);

  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  h = 1.0/n;

  l_sum = 0.;
  for (i=my_rank; i<n; i+=size) {
    x = (i+0.5)*h;
    l_sum += 4.0/(1.0+x*x);
  }
  l_sum *= h;

  if (my_rank==0) {
    MPI_Status status;
    g_sum = l_sum;
    for (i=1; i<size; i++) {
      MPI_Recv(&l_sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,500,MPI_COMM_WORLD,&status);
      g_sum += l_sum;
    }
    printf("result=%g\n",g_sum);
  }
  else
    MPI_Send(&l_sum,1,MPI_DOUBLE,0,500,MPI_COMM_WORLD);

  MPI_Finalize ();
  return 0;
}
