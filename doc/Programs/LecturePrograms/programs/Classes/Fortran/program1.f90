!     Program to compute the second derivative of exp(x). 
!     Three calling functions are included
!     in this version. In one function we read in the data from screen,
!     the next function computes the second derivative
!     while the last function prints out data to screen.

!     The variable h is the step size. We also fix the total number
!     of divisions by 2 of h. The total number of steps is read from
!     screen 


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
USE constants
IMPLICIT NONE
CONTAINS
  SUBROUTINE derivative(number_of_steps, x, initial_step, h_step, &
       computed_derivative)
    USE constants
    INTEGER, INTENT(IN) :: number_of_steps
    INTEGER  :: loop
    REAL(DP), DIMENSION(number_of_steps), INTENT(INOUT) :: &
         computed_derivative, h_step
    REAL(DP), INTENT(IN) :: initial_step, x 
    REAL(DP) :: h
    !     calculate the step size  
    !     initialise the derivative, y and x (in minutes) 
    !     and iteration counter 
    h = initial_step
    ! start computing for different step sizes 
    DO loop=1,  number_of_steps
       !  setup arrays with derivatives and step sizes
       h_step(loop) = h
       computed_derivative(loop) = (EXP(x+h)-2.*EXP(x)+EXP(x-h))/(h*h)
       h = h*0.5
    ENDDO
  END SUBROUTINE derivative

END MODULE functions

PROGRAM second_derivative
  USE constants
  USE functions
  IMPLICIT NONE
  ! declarations of variables 
  INTEGER :: number_of_steps, loop
  REAL(DP) :: x, initial_step
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: h_step, computed_derivative
  !  read in input data from screen 
  WRITE(*,*) 'Read in initial step, x value and number of steps'
  READ(*,*) initial_step, x, number_of_steps
  ! open file to write results on
  OPEN(UNIT=7,FILE='out.dat')
  !  allocate space in memory for the one-dimensional arrays  
  !  h_step and computed_derivative                           
  ALLOCATE(h_step(number_of_steps),computed_derivative(number_of_steps))
  ! compute the second derivative of exp(x)
  h_step = 0.0_dp; computed_derivative = 0.0_dp 
  CALL  derivative(number_of_steps,x,initial_step,h_step,computed_derivative)

  !  Then we print the results to file  
  DO loop=1,  number_of_steps
     WRITE(7,'(E16.10,2X,E16.10)') LOG10(h_step(loop)),&
     LOG10 ( ABS ( (computed_derivative(loop)-EXP(x))/EXP(x)))
     ! free memory
  ENDDO
  DEALLOCATE( h_step, computed_derivative) 
END PROGRAM second_derivative


