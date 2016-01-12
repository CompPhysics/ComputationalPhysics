

PROGRAM standard_deviation
  IMPLICIT NONE
  REAL (KIND = 4) :: sum, sumsq2, xbar
  REAL (KIND = 4) :: sigma1, sigma2
  REAL (KIND = 4), DIMENSION (127) :: x
  INTEGER :: i

  x=0;
  DO i=1, 127
     x(i) = i + 100000.
  ENDDO
  sum=0.; sumsq2=0.
  !      standard deviation calculated with text book algorithm
  DO i=1, 127
     sum = sum +x(i)

     sumsq2 = sumsq2+x(i)**2
  ENDDO
  !      average
  xbar=sum/127.
  sigma1=SQRT((sumsq2-sum*xbar)/126.)
  !      second method to evaluate the standard deviation
  sumsq2=0.
  DO i=1, 127
     sumsq2=sumsq2+(x(i)-xbar)**2
  ENDDO
  sigma2=SQRT(sumsq2/126.)
  WRITE(*,*) xbar, sigma1, sigma2

END PROGRAM standard_deviation
