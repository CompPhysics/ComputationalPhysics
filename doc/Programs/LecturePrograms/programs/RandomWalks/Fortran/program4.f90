!  This program computes the entropy for random walks
!  based on a computed probability histogram.
!  Many walkers (input) make several trials steps with
!  a given number of walks per trial. Periodic boundary conditions
!  have been implemented since the walkers are forced to move on a
!  one-dimensional lattice of size -L  <= x <= L. This means that if the walker
!  goes beyond L(-L) it is shifted to -L(L).

PROGRAM random_entropy 
  USE constants
  USE f90library
  IMPLICIT NONE
  INTEGER::  walkers, tsteps, length
  REAL(dp):: move_probability
  INTEGER, ALLOCATABLE, DIMENSION(:) :: x, probability

  ! get time steps, number of walkers and the probability to move
  CALL  initialise(walkers, tsteps, length, move_probability) 
  ! allocate memory for positions of all walkers
  ALLOCATE (x(walkers))
  ALLOCATE(probability(-length:length))
  x = 0; probability = 0
  ! do the Monte Carlo move
  CALL  mc_sampling(tsteps, walkers, length, move_probability, x, probability)
  ! Print out results 
  CALL  output(tsteps, walkers, length, probability) 
  ! free memory
  DEALLOCATE(x); DEALLOCATE(probability)

END PROGRAM random_entropy
! Read in data from screen

SUBROUTINE initialise(walkers, tsteps, length, move_probability) 
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(inout) :: walkers, tsteps, length
  REAL(DP) :: move_probability
  WRITE(*,*)'Number of walkers =' 
  READ(*,*)  walkers
  WRITE(*,*) 'Number of time steps='
  READ(*,*)  tsteps
  WRITE(*,*) 'Length of systems, from -L to +L='
  READ(*,*)  length
  WRITE(*,*) 'Move probability='
  READ(*,*)  move_probability
END SUBROUTINE initialise

! Output of expectation values, <x>, <y> and their variances
! <x^2> - <x>^2  and <y^2> - <y>^2

SUBROUTINE output(tsteps, walkers, length, probability) 
  USE constants
  IMPLICIT NONE
  INTEGER :: i
  INTEGER, INTENT(in) :: tsteps, walkers, length
  INTEGER, INTENT(in) :: probability(-length:length)
  REAL(dp) ::norm, histogram, entropy

  !find norm of probability
  OPEN(UNIT=6, FILE='prob.dat')    
  norm = 1.0_dp/walkers
  entropy = 0.0_dp; histogram = 0.0_dp
  DO i =-length, length, 1 
     histogram = probability(i)*norm;
     IF ( histogram > 0.0_dp) THEN
        entropy = entropy-histogram*LOG(histogram)
     ENDIF
  ENDDO
  !write entropy
  WRITE(*,'(I8,2X,F12.6)') tsteps, entropy

END SUBROUTINE output

! Here we perform the moves of the walkers as function of time steps with 
! periodic boundary conditions

SUBROUTINE mc_sampling(tsteps, walkers, length, move_probability, x, probability)
  USE constants
  USE f90library 
  IMPLICIT NONE
  INTEGER, INTENT(in) :: tsteps, walkers, length
  INTEGER, INTENT(inout) :: x(walkers)
  INTEGER, INTENT(inout) :: probability(-length:length)
  INTEGER :: idum, i, j, count
  REAL(dp) :: rantest
  REAL(dp), INTENT(in) :: move_probability
  idum=-1 
  ! loop over Monte Carlo trials
  DO i = 1, tsteps
     DO j=1, walkers
        IF (ran0(idum) <= move_probability) THEN
           IF ( x(j) +1 > length) THEN
              x(j) = -length;
           ELSE
              x(j) = x(j) +1;
           ENDIF
        ELSE
           IF ( x(j) -1 < -length) THEN
              x(j) = length;
           ELSE
              x(j) = x(j)+1;
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  ! After the final time step we compute the probabil
  ! by counting the number of walkers at every position
  DO i = -length, length, 1
     count = 0
     DO j = 1, walkers
        IF ( x(j) == i ) count = count + 1;
     ENDDO
     probability(i) = count;
  ENDDO

END SUBROUTINE mc_sampling


