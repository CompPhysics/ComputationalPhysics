!     Crude Monte-Carlo integration

PROGRAM crude_mc_integration
  USE constants
  USE f90library
  IMPLICIT NONE
  INTEGER :: i, idum, samples
  REAL(DP) :: x, sumf, sum2f, fx, sigma

  INTERFACE
     DOUBLE PRECISION FUNCTION  func(x)
       IMPLICIT NONE
       DOUBLE PRECISION , INTENT(IN) :: x

     END FUNCTION func
  END INTERFACE
  WRITE(*,*) ' Enter number of MC sampling points'
  READ(*,*) samples
  !     initialise the total sum, standard deviation and seed idum 
  sumf=0. ; sum2f=0. ; idum=-1
  DO i=1, samples
     x=ran1(idum)
     fx=func(x)
     sumf=sumf+fx
     sum2f=sum2f+fx*fx
  ENDDO
  sumf=sumf/FLOAT(samples)
  sum2f=sum2f/FLOAT(samples)
  sigma=SQRT((sum2f-sumf*sumf)/FLOAT(samples))
  WRITE(*,*) samples, sumf, sigma

END PROGRAM crude_mc_integration

!    The explicit function to be evaluated

REAL(DP)  FUNCTION  func(x)
  USE constants
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  func=4/(1.+x*x)  

END FUNCTION func


