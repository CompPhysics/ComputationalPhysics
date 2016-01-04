!  Sets up integration points for
!  Gaussian quadrature using Laguerre polynomials
SUBROUTINE gauss_laguerre(x,w,n,alf)
  IMPLICIT NONE
  INTEGER :: MAXIT
  INTEGER, INTENT(IN) :: n
  REAL(KIND=8) :: EPS
  REAL(KIND=8), INTENT(IN) :: alf
  REAL(KIND=8), INTENT(INOUT) :: w(n),x(n)
  INTEGER :: i,its,j
  REAL(KIND=8) :: ai,gammln, nn
  REAL(KIND=8) :: p1,p2,p3,pp,z,z1

  maxit = 10; eps = 3.E-14
  nn = n
  DO i=1,n
     IF(i == 1)THEN
        z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
     ELSE IF(i == 2)THEN
        z=z+(15.+6.25*alf)/(1.+9.0*alf+2.5*n)
     ELSE
        ai=i-2
        z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))* &
          (z-x(i-2))/(1.+.3*alf)
     ENDIF
     DO its=1,MAXIT
        p1=1.d0
        p2=0.d0
        DO  j=1,n
           p3=p2
           p2=p1
           p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
        ENDDO
        pp=(n*p1-(n+alf)*p2)/z
        z1=z
        z=z1-p1/pp
        IF(ABS(z-z1) <= eps) GOTO 1
     ENDDO
     PAUSE 'too many iterations in gaulag'
1    x(i)=z
     w(i)=-EXP(gammln(alf+nn)-gammln(nn))/(pp*n*p2)
     write(*,*) x(i), w(i)
  ENDDO

END SUBROUTINE gauss_laguerre


REAL(KIND =8) FUNCTION gammln(xx)
  REAL(KIND=8), INTENT(IN) :: xx
  INTEGER :: j
  REAL(KIND=8) :: ser,tmp,x,y
  REAL(KIND=8), DIMENSION(6) :: cof(6)
  cof(1)=76.18009172947146; cof(2)=-86.50532032941677 
  cof(3)=24.01409824083091; cof(4)=-1.231739572450155
  cof(5)= .1208650973866179E-2; cof(6) = -.5395239384953E-5
  x=xx
  y=x
  tmp=x+5.5
  tmp=(x+0.5)*LOG(tmp)-tmp
  ser=1.000000000190015
  DO j=1,6
     y=y+1.0
     ser=ser+cof(j)/y
  ENDDO
  gammln = tmp+LOG(2.5066282746310005*ser/x)

END FUNCTION gammln

