PROGRAM binary_integer
IMPLICIT NONE
  INTEGER  i, number, terms(32) ! storage of a0, a1, etc, up to 32 bits
  
  WRITE(*,*) 'Give a number to transform to binary notation' 
  READ(*,*) number
! Initialise the terms a0, a1 etc
  terms = 0
! Fortran takes only integer loop variables
  DO i=0, 31
     terms(i) = MOD(number,2)
     number = number/2
  ENDDO
! write out results
  WRITE(*,*) 'Binary representation '
  DO i=0, 31
    WRITE(*,*)' Term nr and value', i, terms(i)
  ENDDO

END PROGRAM binary_integer
