#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

int main()
{
  int dim = 5;
  mat W = zeros<mat>(dim,dim);
  vec wold = zeros<mat>(dim);
  vec wnew = zeros<mat>(dim);
  vec eigenvector = zeros<mat>(dim);
  // Initializing the first vector
  wold(0) = 1.0;
  // Setting up the stochastic matrix W
  W(0,0) = 0.; W(0,1) = 0.; W(0,2) = 0.25; W(0,3) = 0.0;   W(0,4) = 0.;
  W(1,0) = 0.; W(1,1) = 0.; W(1,2) = 0.25; W(1,3) = 0.; W(1,4) = 0.0;
  W(2,0) = 0.5; W(2,1) = 1.0; W(2,2) = 0.; W(2,3) = 0.5; W(2,4) = 0.;
  W(3,0) = 0.0; W(3,1) = 0.; W(3,2) = 0.25; W(3,3) = 0.;  W(3,4) = 0.;
  W(4,0) = 0.5; W(4,1) = 0.; W(4,2) = 0.25; W(4,3) = 0.5;  W(4,4) = 1.0;
  double eps = 1.0E-10;
  W.print("W =");
  double difference  = norm(wold-wnew, 2);
  int count = 0;
  do{
    // Multiplying the old vector with the transition probability
    count += 1;
    wnew = W*wold;
    difference  = norm(wold-wnew, 2);
    wold = wnew;
    cout << "Iteration number = " << count << endl;
    wnew.print("New vector =");
  } while(difference > eps);

  // Getting the eigenvectors and eigenvalues of the stochastic matrix
  cx_vec eigval;
  eig_gen(eigval, W);
  eigval.print("Eigenvalues=");
  return 0;
}
