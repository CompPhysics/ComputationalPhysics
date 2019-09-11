#include "catch.hpp"
#include "Tridiag.h"


TEST_CASE("Testing eigenvalues of Toeplitz matrix"){
    int Dim = 10;
    //    Set up the exact eigenvalues
    vec Exact(Dim);
    double pi = acos(-1.0);
    // Integration step length
    double Step    = 1.0/Dim;
    double DiagConst = 2.0 / (Step*Step);
    double NondiagConst =  -1.0 / (Step*Step);
    for(int i = 0; i < Dim; i++) {
      Exact(i) = DiagConst+2*NondiagConst*cos((i+1)*pi/(Dim+1));
    }
    //get numerical eigenvalues
    vec EigvalueNum(Dim);
    EigvalueNum = GetEigenvalues(Dim);
    
    REQUIRE(EigvalueNum(0)==Approx(Exact(0)).epsilon(0.000000000001));
    REQUIRE(EigvalueNum[1]==Approx(Exact(1)).epsilon(0.000000000001));
    REQUIRE(EigvalueNum[2]==Approx(Exact(2)).epsilon(0.0000000000000001));
}


