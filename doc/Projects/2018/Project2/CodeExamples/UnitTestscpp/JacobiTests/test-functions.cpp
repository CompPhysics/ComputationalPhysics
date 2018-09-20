#include "catch.hpp"
#include "jacobi.h"


TEST_CASE("Testing max a(i,j)"){
    int n=3;
    double pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
    initialize(n,h,a,r,v,0,0);
    int p=0;
    int q=0;
    double apq=0;
    //find maximum matrix element
    find_max(a,p,q,apq,n);
    
    REQUIRE(p==2);
    REQUIRE(q==1);
    REQUIRE(apq==Approx(-0.09));
}

TEST_CASE("Testing eigenvalues"){
    int n=4,interact=0;
    double conv=0.001,wr=0.01, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    n=n-1;
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
    initialize(n,h,a,r,v,0,0);
    //do jacobi algorithm until convergence
    jacobi(n,interact,conv,wr,a,r,v);
    //get eigenvalue vector
    vector<double>eigen=get_eigenvals(a,n);
    
    REQUIRE(eigen[0]==Approx(6.56863));
    REQUIRE(eigen[1]==Approx(25.32055));
    REQUIRE(eigen[2]==Approx(56.57082));
}
TEST_CASE("Testing eigenvector orthogonality"){
    int n=4,interact=0;
    double conv=0.001,wr=0.01, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    n=n-1;
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
    initialize(n,h,a,r,v,0,0);
    //do jacobi algorithm until convergence
    jacobi(n,interact,conv,wr,a,r,v);
    mat eigenvec=get_eigenvecs(a,v,n);

    //test eigen vectors
    REQUIRE(eigenvec(0,0)==Approx(0.9996).epsilon(0.001));
    REQUIRE(eigenvec(0,1)==Approx(0.00853));
    REQUIRE(eigenvec(0,2)==Approx(0.0000).epsilon(0.001));
    REQUIRE(eigenvec(1,0)==Approx(-0.00853));
    REQUIRE(eigenvec(1,1)==Approx(0.99995).epsilon(0.001));
    REQUIRE(eigenvec(1,2)==Approx(0.00512).epsilon(0.001));
    REQUIRE(eigenvec(2,0)==Approx(0.0000).epsilon(0.001));
    REQUIRE(eigenvec(2,1)==Approx(-0.00512).epsilon(0.001));
    REQUIRE(eigenvec(2,2)==Approx(0.9999).epsilon(0.001));
    
    //test eigen vector orthogonality
    //dot1=v0*v1=0
    double dot1=eigenvec(0,0)*eigenvec(1,0)+eigenvec(0,1)*eigenvec(1,1)
        +eigenvec(0,2)*eigenvec(1,2);
    //dot2=v0*v0=1
    double dot2=eigenvec(0,0)*eigenvec(0,0)+eigenvec(0,1)*eigenvec(0,1)
        +eigenvec(0,2)*eigenvec(0,2);
    REQUIRE(dot1==Approx(0.000).epsilon(0.01));
    REQUIRE(dot2==Approx(1.000));
}
