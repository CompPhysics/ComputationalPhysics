!  Program to perform a random walk in 2 dimensions
!  with equal probability for moves to the letf, right, up and down
!  Only one walker present

PROGRAM random_walk_2dim
  USE constants
  USE f90library
  IMPLICIT NONE
  INTEGER::  max_trials, number_walks
  REAL(dp):: move_probability
  INTEGER, ALLOCATABLE, DIMENSION(:) :: x_cum, x2_cum, y_cum, y2_cum

  ! get max trials, number of walks and the probability to move
  CALL  initialise(max_trials, number_walks, move_probability) 
  ! allocate memory
  ALLOCATE (x_cum(number_walks), x2_cum(number_walks), &
       y_cum(number_walks), y2_cum(number_walks))
  x_cum = 0; y_cum = 0; x2_cum = 0; y2_cum = 0  
  ! do the Monte Carlo move
  CALL  mc_sampling(max_trials, number_walks, move_probability, &
       x_cum, x2_cum, y_cum, y2_cum)
  ! Print out results 
  CALL  output(max_trials, number_walks, x_cum, x2_cum, y_cum, y2_cum) 
  ! free memory
  DEALLOCATE(x_cum, x2_cum, y_cum, y2_cum) 
END PROGRAM random_walk_2dim

! Read in data from screen

SUBROUTINE initialise(max_trials, number_walks, move_probability) 
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(inout) ::max_trials, number_walks
  REAL(DP) :: move_probability
  WRITE(*,*)'Number of Monte Carlo trials =' 
  READ(*,*)  max_trials
  WRITE(*,*) 'Number of attempted walks='
  READ(*,*) number_walks
  WRITE(*,*) 'Move probability='
  READ(*,*)  move_probability
END SUBROUTINE initialise

! Output of expectation values, <x>, <y> and their variances
! <x^2> - <x>^2  and <y^2> - <y>^2

SUBROUTINE output(max_trials, number_walks, x_cum, x2_cum, y_cum, y2_cum)
  USE constants
  IMPLICIT NONE
  INTEGER :: i
  INTEGER, INTENT(in) :: max_trials, number_walks
  INTEGER, INTENT(in) :: x_cum(number_walks), x2_cum(number_walks), &
       y_cum(number_walks), y2_cum(number_walks)
  REAL(dp) :: xaverage, yaverage , x2average, y2average, total_average, &
       total_variance, yvariance, xvariance

  OPEN(UNIT=6, FILE='out.dat')
  DO i =1, number_walks 
     xaverage = x_cum(i)/FLOAT(max_trials)
     x2average = x2_cum(i)/FLOAT(max_trials)
     xvariance = x2average - xaverage*xaverage
     yaverage = y_cum(i)/FLOAT(max_trials)
     y2average = y2_cum(i)/FLOAT(max_trials)
     yvariance = y2average - yaverage*yaverage
     total_average = xaverage+yaverage
     total_variance = xvariance+yvariance
     WRITE(6,'(I3,2X,6F12.6)') i, xaverage, yaverage ,xvariance, yvariance, &
          total_average, total_variance
  ENDDO
END SUBROUTINE output

! Here we perform the Monte Carlo samplings and moves in 2 dim

SUBROUTINE mc_sampling(max_trials,number_walks, move_probability, x_cum, &
     x2_cum, y_cum, y2_cum)
  USE constants
  USE f90library
  IMPLICIT NONE
  INTEGER, INTENT(in) :: max_trials,number_walks
  INTEGER, INTENT(inout) :: x_cum(number_walks), x2_cum(number_walks), &
       y_cum(number_walks), y2_cum(number_walks)
  INTEGER :: idum, x, y, trial, walks
  REAL(dp) :: rantest
  REAL(dp), INTENT(in) :: move_probability
  idum=-1 
  ! loop over Monte Carlo trials
  DO trial = 1, max_trials
     x = 0;  y = 0
     ! loop over attempted walks in 2 dimensions
     DO walks = 1, number_walks
        !  note code taylored to equal probability for up, down, left and right
        rantest = ran0(idum)
        IF (rantest <= move_probability) THEN
           x = x+1
        ELSEIF(rantest <= 2*move_probability)THEN
           x = x-1
        ELSEIF(rantest <= 3*move_probability)THEN
           y = y+1
        ELSE
           y = y-1
        ENDIF
        x_cum(walks) = x_cum(walks)+x
        x2_cum(walks) = x2_cum(walks)+x*x
        y_cum(walks) = y_cum(walks)+y
        y2_cum(walks) = y2_cum(walks)+y*y
     ENDDO
  ENDDO
END SUBROUTINE mc_sampling


