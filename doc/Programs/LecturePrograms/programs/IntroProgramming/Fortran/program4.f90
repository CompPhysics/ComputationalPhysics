! In this module you can define for example global constants

MODULE constants
  ! definition of variables for double precisions and complex variables 
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  ! Global Truncation parameter
  REAL(DP), PARAMETER, PUBLIC ::  truncation=1.0E-10
END MODULE constants

! Here you can include specific functions which can be used by
! many subroutines or functions

MODULE functions

CONTAINS
  REAL(DP) FUNCTION factorial(n)
    USE CONSTANTS 
    INTEGER, INTENT(IN) :: n
    INTEGER  :: loop

    factorial = 1.0_dp
    IF ( n > 1 ) THEN
       DO loop = 2, n
          factorial=factorial*loop
       ENDDO
    ENDIF

  END FUNCTION factorial

END MODULE functions


PROGRAM exp_prog
  USE constants
  USE functions
  IMPLICIT NONE  
  REAL (DP) :: x, term, final_sum
  INTEGER :: n, loop_over_x

  !  loop over x-values
  DO loop_over_x=0, 100, 10
     x=loop_over_x
     !  initialize the EXP sum
     final_sum= 0.0_dp; term = 1.0_dp; n = 0 
     DO WHILE ( ABS(term) > truncation)
        term = ((-1.0_dp)**n)*(x**n)/ factorial(n)
        final_sum=final_sum+term
        n=n+1
     ENDDO
     !  write the argument x, the exact value, the computed value and n
     WRITE(*,*) x ,EXP(-x), final_sum, n
  ENDDO

END PROGRAM exp_prog


