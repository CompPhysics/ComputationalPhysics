!  Simple test program for numerical integration which
!  calls three methods in the program library (you must link
!  the library as well).


PROGRAM int_test
  USE constants
  USE F90library
  IMPLICIT NONE
  INTEGER :: i, n
  REAL(DP) :: a, b, int_trapez, int_gauss, int_simpson
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: x, w

  INTERFACE
     DOUBLE PRECISION FUNCTION  func(x)
       USE constants
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x

     END FUNCTION func
  END INTERFACE

  WRITE(*,*) ' Read in number of mesh points'
  READ(*,*) n
  WRITE(*,*) ' Read in integration limits [a,b]'
  READ(*,*) a, b

  !    reserve space in memory for vectors containing the mesh points
  !    weights and function values for the use of the Gauss-Legendre
  !    method

  ALLOCATE ( x(n), w(n) )

  !    Set up the mesh points and weights for Gauss Legendre

  CALL gauleg(a,b,x,w,n)

  !    Integrate using the trapezoidal rule and simpson's method
  !    Note the transfer of a function name to the methods
  CALL trapezoidal_rule(a,b,int_trapez,n,func)
  CALL simpson(a,b,int_simpson,n,func)
  !    Gaussian quadrature

  int_gauss=0.
  DO i=1,n
     int_gauss=int_gauss+w(i)*func(x(i))
  ENDDO

  !    final output
  WRITE (*, *) n, int_trapez, int_simpson, int_gauss
  DEALLOCATE ( x, w)

END PROGRAM int_test

!    The explicit function to be evaluated

DOUBLE PRECISION FUNCTION  func(x)
  USE constants
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  func=exp(-x)/x

END FUNCTION func
