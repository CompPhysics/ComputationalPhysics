%module MC
%{
#include "MC.cpp"
%}

class data {

public:
  double sum;
  double squaresum;
  int N;
  int accepted;
  long idum;
};

long seed();

data runMC(int MCcycles, double delta, long idum, double alpha);

