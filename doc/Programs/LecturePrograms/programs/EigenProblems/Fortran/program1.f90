!     This version solves the Schroedinger equation for the 
!     one-dimensional 
!     harmonic oscillator using matrix diagonalization

!     declaration of constants used in the calculations
!     define variables, where we have chosen to use atomic units for 
!     m=c=hbar=k=1.  

!    This module contains all constants and declarations 
!    of variables read in by the function read_data. These
!    variables are used by many functions.

MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
END MODULE constants

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
!  USE f90library
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







  REAL (DP) FUNCTION pythag(a,b)
    USE constants
    REAL(DP)  :: a,b
    REAL(DP)  :: absa,absb
    absa=ABS(a)
    absb=ABS(b)
    IF(absa > absb) THEN
       pythag=absa*sqrt(1.+(absb/absa)**2)
    ELSE
       IF(absb == 0.) THEN
          pythag=0.
       ELSE
          pythag=absb*sqrt(1.+(absa/absb)**2)
       ENDIF
    ENDIF

  END FUNCTION pythag


  !     determine eigenvalues and eigenvectors of a real symmetric
  !     tri-diagonal matrix, or a real, symmetric matrix previously
  !     reduced by function tred2 to tri-diagonal form. On input,
  !     d[] contains the diagonal element and e[] the sub-diagonal
  !     of the tri-diagonal matrix. On output d[] contains the
  !     eigenvalues and  e[] is destroyed. If eigenvectors are
  !     desired z[][] on input contains the identity matrix. If
  !     eigenvectors of a matrix reduced by tred2() are required,
  !     then z[][] on input is the matrix output from tred2().
  !     On output, the k'th column returns the normalized eigenvector
  !     corresponding to d[k]. 
  !     The function is modified from the version in Numerical recipe.

  SUBROUTINE tqli(d,e,n,z)
    USE constants
    IMPLICIT NONE
    INTEGER :: n 
    REAL(DP)  :: d(n),e(n),z(n,n)
    INTEGER :: i,iter,k,l,m
    REAL(DP)  :: b,c,dd,f,g,p,r,s,pythag,one
    DO i=2,n
       e(i-1)=e(i)
    ENDDO
    one=1.
    DO l=1,n
       iter=0
       ITERATE : DO
          DO m=l,n-1
             dd=ABS(d(m))+ABS(d(m+1))
             IF (ABS(e(m))+dd == dd) EXIT
          ENDDO
          IF(m == l) EXIT ITERATE
          IF(iter == 30) STOP 'too many iterations in tqli'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.*e(l))
          r=pythag(g,one)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.
          c=1.
          p=0.
          DO i=m-1,l,-1
             f=s*e(i)
             b=c*e(i)
             r=pythag(f,g)
             e(i+1)=r
             IF(r == 0.) THEN
                d(i+1)=d(i+1)-p
                e(m)=0.
                CYCLE ITERATE
             ENDIF
             s=f/r
             c=g/r
             g=d(i+1)-p
             r=(d(i)-g)*s+2.*c*b
             p=s*r
             d(i+1)=g+p
             g=c*r-b
             !     Omit lines from here ...
             DO k=1,n
                f=z(k,i+1)
                z(k,i+1)=s*z(k,i)+c*f
                z(k,i)=c*z(k,i)-s*f
             ENDDO
             !     ... to here when finding only eigenvalues.
          ENDDO
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.
       ENDDO ITERATE
    ENDDO

  END SUBROUTINE tqli

