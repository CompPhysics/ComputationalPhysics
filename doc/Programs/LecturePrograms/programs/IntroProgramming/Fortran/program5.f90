! In this module you can define for example global constants

MODULE constants
  ! definition of variables for double precisions and complex variables 
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  ! Global Truncation parameter
  REAL(DP), PARAMETER, PUBLIC ::  truncation=1.0E-10
END MODULE constants

PROGRAM improved_exp
  USE constants
  IMPLICIT NONE  
  REAL (dp) :: x, term, final_sum
  INTEGER  :: n, loop_over_x

  !  loop over x-values, no floats as loop variables
  DO loop_over_x=0, 100, 10
     x=loop_over_x
     !  initialize the EXP sum
     final_sum=1.0 ; term=1.0 ; n = 1
     DO WHILE ( ABS(term) > truncation)
        term = -term*x/FLOAT(n)
        final_sum=final_sum+term
        n=n+1
     ENDDO
     !  write the argument x, the exact value, the computed value and n
     WRITE(*,*) x ,EXP(-x), final_sum, n
  ENDDO

END PROGRAM improved_exp
