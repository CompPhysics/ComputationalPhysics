#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <armadillo>
using namespace std;
using namespace arma;

double ** CreateMatrix(int m, int n){
  double ** mat;
  mat = new double*[m];
  for(int i=0;i<m;i++){
    mat[i] = new double[n];
    for(int j=0;j<m;j++)
      mat[i][j] = 0.0;
  }
  return mat;
}

void DestroyMatrix(double ** mat, int m, int n){
  for(int i=0;i<m;i++)
    delete[] mat[i];
  delete[] mat;
}

int Jacobi_P(int, int, int, double **, double *, double *, double);

int main(int argc, char * argv[]){
  int i,j, N = 20;
  double **A,*x,*q;
  int totalnodes,mynode;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  if(mynode==0){
    A = CreateMatrix(N,N);
    x = new double[N];
    q = new double[N];

    for(i=0;i<N;i++){
      q[i] = i+1;
      A[i][i] = -2.0;
      if(i<N-1){
        A[i][i+1] = 1.0;
        A[i+1][i] = 1.0;
      }
    }
  }
  Jacobi_P(mynode,totalnodes,N,A,x,q,1.0e-14);
  if(mynode==0){
    for(i=0;i<N;i++)
      cout << x[i] << endl;
    DestroyMatrix(A,N,N);
    delete[] x;
    delete[] q;
  }
  MPI_Finalize();
}


int Jacobi_P(int mynode, int numnodes, int N, double **A, double *x, double *b, double abstol){
  int i,j,k,i_global;
  int maxit = 100000;
  int rows_local,local_offset,last_rows_local,*count,*displacements;
  double sum1,sum2,*xold;
  double error_sum_local, error_sum_global;
  MPI_Status status;

  rows_local = (int) floor((double)N/numnodes);
  local_offset = mynode*rows_local;
  if(mynode == (numnodes-1)) 
    rows_local = N - rows_local*(numnodes-1);

  /*Distribute the Matrix and R.H.S. among the processors */
  if(mynode == 0){
    for(i=1;i<numnodes-1;i++){
      for(j=0;j<rows_local;j++)
        MPI_Send(A[i*rows_local+j],N,MPI_DOUBLE,i,j,MPI_COMM_WORLD);
      MPI_Send(b+i*rows_local,rows_local,MPI_DOUBLE,i,rows_local,
               MPI_COMM_WORLD);
    }
    last_rows_local = N-rows_local*(numnodes-1);
    for(j=0;j<last_rows_local;j++)
      MPI_Send(A[(numnodes-1)*rows_local+j],N,MPI_DOUBLE,numnodes-1,j,
               MPI_COMM_WORLD);
    MPI_Send(b+(numnodes-1)*rows_local,last_rows_local,MPI_DOUBLE,numnodes-1,
             last_rows_local,MPI_COMM_WORLD);
  }
  else{
    A = CreateMatrix(rows_local,N);
    x = new double[rows_local];    
    b = new double[rows_local];
    for(i=0;i<rows_local;i++)
      MPI_Recv(A[i],N,MPI_DOUBLE,0,i,MPI_COMM_WORLD,&status);
    MPI_Recv(b,rows_local,MPI_DOUBLE,0,rows_local,MPI_COMM_WORLD,&status);
  }


  xold = new double[N];
  count = new int[numnodes];
  displacements = new int[numnodes];


  //set initial guess to all 1.0
  for(i=0; i<N; i++){
    xold[i] = 1.0;
  }

  for(i=0;i<numnodes;i++){
    count[i] = (int) floor((double)N/numnodes);
    displacements[i] = i*count[i];
  }
  count[numnodes-1] = N - ((int)floor((double)N/numnodes))*(numnodes-1);
  
  for(k=0; k<maxit; k++){
    error_sum_local = 0.0;
    for(i = 0; i<rows_local; i++){
      i_global = local_offset+i;
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i_global; j++)
        sum1 = sum1 + A[i][j]*xold[j];
      for(j=i_global+1; j < N; j++)
        sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (-sum1 - sum2 + b[i])/A[i][i_global];
      error_sum_local += (x[i]-xold[i_global])*(x[i]-xold[i_global]);
    }
    
    MPI_Allreduce(&error_sum_local,&error_sum_global,1,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);
    MPI_Allgatherv(x,rows_local,MPI_DOUBLE,xold,count,displacements,
                   MPI_DOUBLE,MPI_COMM_WORLD);
    
    if(sqrt(error_sum_global)<abstol){
      if(mynode == 0){
        for(i=0;i<N;i++)
          x[i] = xold[i];
      }
      else{
        DestroyMatrix(A,rows_local,N);
        delete[] x;
        delete[] b;
      }
      delete[] xold;
      delete[] count;
      delete[] displacements;
      return k;
    }
  }

  cerr << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
  if(mynode == 0){
    for(i=0;i<N;i++)
      x[i] = xold[i];
  }
  else{
    DestroyMatrix(A,rows_local,N);
    delete[] x;
    delete[] b;
  }
  delete[] xold;
  delete[] count;
  delete[] displacements;
  
  return maxit;
}



