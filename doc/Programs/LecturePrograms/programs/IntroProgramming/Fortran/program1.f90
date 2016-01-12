PROGRAM shw
  IMPLICIT NONE
  REAL (KIND =8) :: r         ! Input number
  REAL (KIND=8)  :: s         ! Result

  !  Get a number from user
  WRITE(*,*) 'Input a number: '
  READ(*,*) r
  !  Calculate the sine of the number
  s = SIN(r)
  !  Write result to screen
  WRITE(*,*) 'Hello World! SINE of ', r, ' =', s

END PROGRAM shw
