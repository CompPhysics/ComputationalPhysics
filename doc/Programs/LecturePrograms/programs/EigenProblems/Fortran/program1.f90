!     This version solves the Schroedinger equation for the 
!     one-dimensional 
!     harmonic oscillator using matrix diagonalization

!     declaration of constants used in the calculations
!     define variables, where we have chosen to use atomic units for 
!     m=c=hbar=k=1.  

MODULE parameters     
  USE constants
  !  the step
  REAL(DP), PUBLIC :: step
  !  r min
  REAL(DP), PUBLIC, PARAMETER :: rmin=-10.0_dp
  !  r max
  REAL(DP), PUBLIC, PARAMETER :: rmax=10.0_dp
  !  number of steps from rmin to rmin to rmax
  INTEGER, PUBLIC, PARAMETER ::  max_step=1000
  !  size of matrix to be diagonalized, do not need points at rmin and rmax
  INTEGER, PUBLIC, PARAMETER ::  n=998
END MODULE parameters

!  main program starts here

PROGRAM schroedinger_equation
  USE constants
  USE parameters
  USE f90library
  IMPLICIT NONE
  REAL(DP) :: potential, const1, const2
  REAL(DP), DIMENSION(:), ALLOCATABLE :: w, r
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: e, d
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: z
  INTEGER :: i

  !  reserve place in memory for the arrays w, r, e, d and z

  ALLOCATE ( e(n) , d(n))
  ALLOCATE ( w(0:max_step), r(0:max_step))

  !  define the step size

  step=(rmax-rmin)/FLOAT(max_step)

  !  define constants for the matrix to be diagonalized

  const1=2./(step*step)
  const2=-1./(step*step)

  !     set up r and the function w for energy =0
  !     w corresponds then to the potential
  !  values at 

  DO i=0, max_step
     r(i) = rmin+i*step
     w(i) = potential(r(i))
  ENDDO

  !     setup the diagonal d and the non-diagonal part e  of
  !     the  tri-diagonal matrix to be diagonalized

  d(1:n)=const1+w(1:n)  ;  e(1:n)=const2

  !  allocate space for eigenvector info

  ALLOCATE ( z(n,n) )

  !  obtain the eigenvalues

  CALL tqli(d,e,n,z)

  !  sort eigenvalues as an ascending series 

  CALL eigenvalue_sort(d,n)

  !  write out 

  WRITE(6,*) ' Rmin and Rmax ', rmin, rmax
  WRITE(6,*) ' Number of steps -1 ', n
  WRITE(6,*) ' The 5 lowest eigenvalues '
  WRITE(6,'(5F12.6)') d(1:5)

  !  find the wave function for the lowest lying state

  CALL plot(d(1), w, r) 

  !  deallocate space in memory

  DEALLOCATE (z)
  DEALLOCATE ( w, r, d, e)

END PROGRAM schroedinger_equation

!  Harmonic oscillator potential

REAL(DP) FUNCTION potential(x)
  USE constants
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x

  potential=x*x

END FUNCTION  potential


!
!           Sort the eigenvalues as
!           ascending series 
!

SUBROUTINE eigenvalue_sort(e,n)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  REAL(DP), DIMENSION(n), INTENT(INOUT) :: e
  REAL(DP) :: temp

  DO i = 1, n
     DO j =i, n
        IF( ABS(e(i)) >   ABS(e(j)) ) THEN
           ! change the energy  
           temp=e(i)
           e(i)=e(j)
           e(j)=temp
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE eigenvalue_sort



!  This function plots the wave function at the final energy

SUBROUTINE plot(energy, w ,r )
  USE constants
  USE parameters
  IMPLICIT NONE
  INTEGER :: i, match
  REAL(DP), DIMENSION(n) :: w, r
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wf
  REAL(DP) :: energy, wwf, fac, norm, pi

  ALLOCATE ( wf(0:max_step))
  pi=ACOS(-1.) ; pi=SQRT(pi)

  !  add the chosen energy to the potential 

  w=(w-energy)*step*step+2

  !  Boundary conditions
  !  value of the wave function at max step

  wf(max_step)=0.

  !  here we use the analytical form for the wave function for max step -1
  !  (taylor expansion of EXP(-r^2/2))

  wf(max_step-1)=step**2/2.

  !  start from right and move towards rmin

  OUTER_TO_INNER : DO
     DO i=max_step-2, 1, -1
        wf(i)=w(i+1)*wf(i+1)-wf(i+2)
        !  find point were the outer part starts turning in order to determine matching point
        IF ( wf(i) <=  wf(i+1) ) EXIT OUTER_TO_INNER  
     ENDDO
  ENDDO OUTER_TO_INNER
  match = i+1
  wwf=wf(match)

  !  start from rmin  and move towards the matching point

  wf(0)=0.
  wf(1)=step**2/2.
  DO i=2, match
     wf(i)=w(i-1)*wf(i-1)-wf(i-2)
  ENDDO

  !  scale the wave function below the matching point

  fac=wwf/wf(match)
  wf(0:match)=wf(0:match)*fac

  !  normalize the wave function

  norm=step*SUM(wf*wf)
  IF ( norm > 0. ) norm = 1./SQRT(norm)
  wf=wf*norm
  WRITE(6,*) norm 
  !  write out of numerically calculated and analytical expression for the wave functions

  WRITE(6,'(E12.6,2X,E12.6,2X,E12.6)')( r(i), wf(i), EXP(-(r(i)**2)/2.0_dp)/SQRT(pi),i=1, max_step)
  DEALLOCATE ( wf) 

END SUBROUTINE plot












