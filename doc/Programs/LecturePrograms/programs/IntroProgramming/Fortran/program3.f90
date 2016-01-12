PROGRAM integer_exp
  IMPLICIT NONE
  INTEGER (KIND=4) :: int1, int2, int3

  ! This is the begin of a comment line in Fortran 90
  ! Now we read from screen the variable int2

  WRITE(*,*) 'Read in the number to be exponentiated'   
  READ(*,*) int2 
  int1=2**int2
  WRITE(*,*) '2^N*2^N', int1*int1
  int3=int1-1
  WRITE(*,*) '2^N*(2^N-1)', int1*int3
  WRITE(*,*) '2^N-1', int3

END PROGRAM integer_exp
