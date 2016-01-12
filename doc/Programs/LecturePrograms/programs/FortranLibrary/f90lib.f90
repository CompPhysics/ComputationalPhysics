!    This module contains all constants and declarations 
!    of variables read in by the function read_data. These
!    variables are used by many functions.

MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
END MODULE constants


!    F90 program library, adapted from Numerical Recipes
!    All functions have been translated to F90 from F77

MODULE F90library
  USE constants
CONTAINS

  !
  !            Routines to do mtx inversion, from Numerical
  !            Recipes, Teukolsky et al. Routines included
  !            below are MATINV, LUDCMP and LUBKSB. See chap 2
  !            of Numerical Recipes for further details
  !
  SUBROUTINE matinv(a,n,d)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j
    REAL(DP), DIMENSION(n,n), INTENT(INOUT)  :: a
    REAL(DP), ALLOCATABLE :: y(:,:)
    REAL(DP) :: d
    INTEGER, ALLOCATABLE :: indx(:)

    ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
    y=0.
    !     setup identity matrix
    DO i=1,n
       y(i,i)=1.
    ENDDO
    !     LU decompose the matrix just once
    CALL  lu_decompose(a,n,indx,d)

    !     Find inverse by columns
    DO j=1,n
       CALL lu_linear_equation(a,n,indx,y(:,j))
    ENDDO
    !     The original matrix a was destroyed, now we equate it with the inverse y 
    a=y

    DEALLOCATE ( y ); DEALLOCATE ( indx )

  END SUBROUTINE matinv

  !     Given an NxN matrix A(N,N), this routine replaces it by the LU 
  !     decomposed one, where the matrix elements are stored in the same 
  !     matrix A. The array indx is  an output vector which records the row
  !     permutation effected by the partial pivoting. d is the determinant
  !
  SUBROUTINE lu_decompose(a,n,indx,d)
    IMPLICIT NONE
    INTEGER :: n, i, j, k, imax
    REAL(DP) :: sum , tiny, aamax, dum, d
    REAL(DP), DIMENSION(n,n) :: a
    INTEGER, DIMENSION(n) :: indx
    REAL(DP), ALLOCATABLE :: vv(:)

    tiny=1.0e-20
    ALLOCATE ( vv(n) )
    D=1.
    DO i=1,n
       aamax=0.
       DO j=1,n
          IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
       ENDDO
       !     Zero is the largest element
       IF (aamax == 0.) STOP 'Singular matrix.'
       !     No nonzero largest element
       vv(i)=1./aamax
    ENDDO
    !     loop over columns
    DO j=1,n
       !     solves equation 2.3.12 except for i=j of Numerical Recipes
       IF (j > 1) THEN
          DO i=1,j-1
             sum=a(i,j)
             IF (i > 1)THEN
                DO k=1,i-1
                   sum=sum-a(i,k)*a(k,j)
                ENDDO
                a(i,j)=sum
             ENDIF
          ENDDO
       ENDIF
       !    start searching for largest pivot element
       aamax=0.
       DO i=j,n
          sum=a(i,j)
          IF (j > 1)THEN
             DO k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             ENDDO
             a(i,j)=sum
          ENDIF
          dum=vv(i)*ABS(sum)
          IF (dum >= aamax) THEN
             imax=i
             aamax=dum
          ENDIF
       ENDDO
       !    interchange of rows
       IF (j /= imax)THEN
          DO k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          ENDDO
          !    change of parity for determinant
          d=-d
          vv(imax)=vv(j)
       ENDIF
       indx(j)=imax
       IF(j /= n) THEN
          IF(a(j,j) == 0.) a(j,j)=tiny
          dum=1./a(j,j)
          DO i=j+1,n
             a(i,j)=a(i,j)*dum
          ENDDO
       ENDIF
       !    set up determinant
       d=d*a(j,j)
    ENDDO
    IF(a(n,n) == 0.)  a(n,n)=tiny
    DEALLOCATE ( vv)

  END SUBROUTINE lu_decompose

  !     Solves set of linear equations Ax=b, A is input as an LU decompomsed
  !     matrix and indx keeps track of the permutations of the rows. b is input
  !     as the right-hand side vector b and returns the solution x. A, n and indx
  !     are not modified by this routine. This function takes into that b can contain
  !     many zeros and is therefore suitable for matrix inversion


  SUBROUTINE lu_linear_equation(a,n,indx,b)
    IMPLICIT NONE
    INTEGER :: n, ii, ll, i, j
    REAL(DP) :: sum 
    REAL(DP), DIMENSION(n,n) :: a
    REAL(DP), DIMENSION(n) :: b
    INTEGER, DIMENSION(n) :: indx

    ii=0
    !     First we solve equation 2.3.6 of numerical recipes 
    DO i=1,n
       ll=indx(i)
       sum=b(ll)
       b(ll)=b(i)
       IF (ii /= 0)THEN
          DO j=ii,i-1
             sum=sum-a(i,j)*b(j)
          ENDDO
       ELSEIF (sum /= 0.) THEN
          ii=i
       ENDIF
       b(i)=sum
    ENDDO
    !     then we solve equation 2.3.7
    DO i=n,1,-1
       sum=b(i)
       IF (i < n) THEN
          DO j=i+1,n
             sum=sum-a(i,j)*b(j)
          ENDDO
       ENDIF
       !     store a component of the solution x in the same place as b
       b(i)=sum/a(i,i)
    ENDDO

  END SUBROUTINE lu_linear_equation

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


  REAL(DP) FUNCTION pythag(a,b)
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

  !    perform a Housholder reduction of a real symmetric matrix
  !    a[][]. On output a[][] is replaced by the orthogonal matrix 
  !    effecting the transformation. d[] returns the diagonal elements
  !    of the tri-diagonal matrix, and e[] the off-diagonal elements, 
  !    with e[0] = 0.
  !    The function is modified from the version in Numerical recipes.

  SUBROUTINE tred2(a,n,d,e)
    IMPLICIT NONE
    INTEGER :: n
    REAL(DP) :: a(n,n),d(n),e(n)
    INTEGER :: i,j,k,l
    REAL(DP) ::  f,g,h,hh,scale

    DO i=n,2,-1
       l=i-1
       h=0.
       scale=0.
       IF (l > 1) THEN
          scale=SUM(abs(a(i,1:l)))
          IF (scale == 0.) THEN
             e(i)=a(i,l)
          ELSE
             a(i,1:l)=a(i,1:l)/scale
             h=sum(a(i,1:l)**2)
             f=a(i,l)
             g=-sign(sqrt(h),f)
             e(i)=scale*g
             h=h-f*g
             a(i,l)=f-g
             f=0.
             DO j=1,l
                !     Omit following line if finding only eigenvalues
                a(j,i)=a(i,j)/h
                g=0.
                DO k=1,j
                   g=g+a(j,k)*a(i,k)
                ENDDO
                DO  k=j+1,l
                   g=g+a(k,j)*a(i,k)
                ENDDO
                e(j)=g/h
                f=f+e(j)*a(i,j)
             ENDDO
             hh=f/(h+h)
             DO j=1,l
                f=a(i,j)
                g=e(j)-hh*f
                e(j)=g
                DO k=1,j
                   a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          e(i)=a(i,l)
       ENDIF
       d(i)=h
    ENDDO
    !     Omit following line if finding only eigenvalues.
    d(1)=0.
    e(1)=0.
    DO i=1,n
       !     Delete lines from here ...
       l=i-1
       IF (d(i) /= 0.) THEN
          DO j=1,l
             g=0.
             DO k=1,l
                g=g+a(i,k)*a(k,j)
             ENDDO
             DO k=1,l
                a(k,j)=a(k,j)-g*a(k,i)
             ENDDO
          ENDDO
       endif
       !     ... to here when finding only eigenvalues.
       d(i)=a(i,i)
       !     Also delete lines from here ...
       a(i,i)=1.
       DO j=1,l
          a(i,j)=0.
          a(j,i)=0.
       ENDDO
       !     ... to here when finding only eigenvalues.
    ENDDO

  END SUBROUTINE tred2

  ! takes as input xa[1,..,n] and ya[1,..,n] together with 
  ! a given value of x and returns a value y and an error 
  ! estimate dy. If P(x) is a polynomial of degree N - 1 such 
  ! that P(xa_i) = ya_i, i = 0,..,n-1, then the returned 
  ! value is y = P(x). 

  SUBROUTINE polint(xa,ya,n,x,y,dy)
    IMPLICIT NONE
    INTEGER :: i, ns, m, n
    REAL(DP), DIMENSION(n) :: xa,ya
    REAL(DP), DIMENSION(:), ALLOCATABLE :: c , d
    REAL(DP) :: den ,dif, dift, ho, hp, w, x, y, dy

    ALLOCATE ( c(n), d(n) )
    ns=1 
    dif=ABS(x-xa(1)) 
    DO i=1,n 
       dift=ABS(x-xa(i)) 
       IF (dift < dif) THEN 
          ns=i 
          dif=dift 
       ENDIF
       c(i)=ya(i) 
       d(i)=ya(i) 
    ENDDO
    y=ya(ns) 
    ns=ns-1 
    DO m=1,n-1 
       DO i=1,n-m 
          ho=xa(i)-x 
          hp=xa(i+m)-x 
          w=c(i+1)-d(i) 
          den=ho-hp 
          IF(den == 0.) WRITE(6,*)  'Error in POLINT, den =0'
          den=w/den 
          d(i)=hp*den 
          c(i)=ho*den 
       ENDDO
       IF (2*ns < n-m)THEN 
          dy=c(ns+1) 
       ELSE 
          dy=d(ns) 
          ns=ns-1 
       ENDIF
       y=y+dy 
    ENDDO
    DEALLOCATE ( c, d)

  END SUBROUTINE polint

  ! takes as input x[1,..,n] and y[1,..,n] containing a tabulation
  ! y_i = f(x_i) with x_0 < x_1 < .. < x_(n - 1) 
  ! together with yp_1 and yp2 for first derivatives  f(x) at x_0 
  ! and x_(n-1), respectively. Then the
  ! function returns y2[1,..,n] which contains the second 
  ! derivatives of f(x_i)at each point x_i. If yp1 and/or yp2 
  ! is larger than the constant INFINITY the function will 
  ! put corresponding second derivatives to zero.

  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
    IMPLICIT NONE
    INTEGER :: i, k, n
    REAL(DP), DIMENSION(n) :: x, y, y2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: u
    REAL(DP) :: p, qn, sig, un, ypn, yp1 

    ALLOCATE ( u (n) )
    IF (yp1 > .99E30) THEN
       y2(1)=0.
       u(1)=0.
    ELSE
       y2(1)=-0.5
       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    ENDIF
    DO i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.
       y2(i)=(sig-1.)/p
       u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
            /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    ENDDO
    IF (ypn > .99E30) THEN
       qn=0.
       un=0.
    ELSE
       qn=0.5
       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    ENDIF
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    DO k=n-1,1,-1
       y2(k)=y2(k)*y2(k+1)+u(k)
    ENDDO
    DEALLOCATE ( u )

  END SUBROUTINE spline


  ! takes xa[1,..,n] and y[1,..,n] which tabulates a function 
  ! (with the xa[i]'s in order) and given ya[0,..,n - 1], 
  ! which is the output from function spline() and with 
  ! given value of x returns a cubic--spline interpolation value y.

  SUBROUTINE splint(xa,ya,y2a,n,x,y)
    IMPLICIT NONE
    INTEGER :: k, n, klo, khi
    REAL(DP), DIMENSION(n) :: xa, ya, y2a
    REAL(DP) ::  x, y, h, b, a

    klo=1
    khi=n
    DO WHILE  (khi-klo > 1)
       k=(khi+klo)/2
       IF(xa(k)  > x)THEN
          khi=k
       ELSE
          klo=k
       ENDIF
    ENDDO
    h=xa(khi)-xa(klo)
    IF (h == 0.) WRITE (6,*) 'Bad XA input in SPLINT.F.'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+ &
         ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

  END SUBROUTINE splint


  !      This routine calculates gauss-legendre mesh points and weights
  !      input: 
  !      x1   : lower limit of the integration interval
  !      x2   : upper limit ---------- "" -------------
  !      n    : the desired number of mesh points
  !      output :
  !      x     : gauss-legendre mesh points on the interval (x1,x2)   
  !      w     : the corresponding weights

  SUBROUTINE gauleg(x1,x2,x,w,n)
    IMPLICIT NONE
    INTEGER :: i, j, m, n
    REAL(DP) :: eps, x1, x2, x, w 
    DIMENSION :: x(n), w(n) 
    PARAMETER (eps=3.D-14)
    REAL(DP) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    DO i=1,m
       z1=0.
       z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.
          p2=0.
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1.)
          z1=z
          z=z-p1/pp
       ENDDO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.*xl/((1.-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO

  END SUBROUTINE gauleg

  !     Function to integrate a function func over the 
  !     interval [a,b] with input a, b, and the number of steps
  !     n.  It returns the sum as the variable trapez_sum
  !     The trapezoidal rule is used


  SUBROUTINE trapezoidal_rule(a,b,trapez_sum,n,func)
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: a,b
    REAL(DP), INTENT(INOUT) :: trapez_sum
    REAL(DP) fa, fb, x, step
    INTEGER :: j
    INTERFACE
       DOUBLE PRECISION FUNCTION  func(x)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: x 

       END FUNCTION func
    END INTERFACE


    step=(b-a)/FLOAT(n)
    fa=func(a)/2. ; fb=func(b)/2. ; trapez_sum=0.
    DO j=1,n-1
       x=j*step+a
       trapez_sum=trapez_sum+func(x)
    ENDDO
    trapez_sum=(trapez_sum+fb+fa)*step

  END  SUBROUTINE trapezoidal_rule


  !     Function to integrate a function func over the 
  !     interval [a,b] with input a, b, and the number of steps
  !     n.  It returns the sum as the variable simpson_sum
  !     Simpson's method is used


  SUBROUTINE simpson(a,b,simpson_sum,n,func)
    REAL(DP), INTENT(IN) :: a,b
    REAL(DP), INTENT(INOUT) :: simpson_sum
    REAL(DP) fa, fb, x, step, fac
    INTEGER, INTENT(IN) :: n 
    INTEGER :: j

    INTERFACE
       DOUBLE PRECISION FUNCTION  func(x)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: x 

       END FUNCTION func
    END INTERFACE

    step=(b-a)/FLOAT(n)
    fa=func(a) ; fb=func(b) ; simpson_sum=fa ; fac=2.
    DO j=1,n-1
       IF ( fac == 2.) THEN
          fac = 4.
       ELSE
          fac = 2.
       ENDIF
       x=j*step+a
       simpson_sum=simpson_sum+func(x)*fac
    ENDDO
    simpson_sum=(simpson_sum+fb)*step/3.

  END SUBROUTINE simpson

  !
  !  4th-Runge-Kutta solution of coupled equations     
  !
  SUBROUTINE rk4(x,y,yout,dydx,step,n)
    IMPLICIT NONE
    INTEGER :: n
    REAL(DP), DIMENSION(n) :: yt, dyt, dym
    REAL(DP), DIMENSION(n), INTENT(IN) :: y, dydx
    REAL(DP), DIMENSION(n), INTENT(OUT) :: yout
    REAL(DP) :: hh, h6, xh, step
    REAL(DP), INTENT(IN) :: x 

    hh=step*0.5; h6=step/6. ; xh=x+hh
    !     first rk-step
    yt=y+hh*dydx
    !      CALL derivs(xh,yt,dyt)
    !     second rk-step
    yt=y+hh*dyt
    !      CALL derivs(xh,yt,dym)      
    !     third rk-step
    yt=y+step*dym;  dym=dyt+dym
    !      CALL derivs(x+step,yt,dyt)
    !     fourth rk-step
    yout=y+h6*(dydx+dyt+2.*dym)

  END SUBROUTINE rk4


  !  This function sets up the recursive relation
  !  for the associated Legendre polynomials

  REAL(DP) FUNCTION legendre_polynomials(l, m, x)
    IMPLICIT NONE
    REAL(DP) ::  fact,pll,pmm,pmmp1,somx2
    REAL(DP), INTENT(IN)  :: x
    INTEGER ::  i,ll
    INTEGER, INTENT(IN) :: l, m

    !  check whether m, l and x are ok

    IF((M < 0).OR.(M > L).OR.(ABS(X) > 1.)) THEN
       WRITE(6,*) 'bad arguments', m, l, x; RETURN
    ENDIF

    !  calculate now pmm as starting point for iterations

    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0-x)*(1.0+x))
       fact=1.0;
       DO i=1, m
          pmm = -fact*somx2*pmm
          fact = fact+2.0
       ENDDO
    ENDIF

    !  if l == m we do not need to use recursion relation

    IF (l == m) THEN
       legendre_polynomials=pmm

       !  recursive relation for associated Legendre polynomials

    ELSE
       pmmp1=x*(2*m+1)*pmm

       !  analytical formula for the case l == m+1

       IF (l == (m+1)) THEN
          legendre_polynomials=pmmp1
       ELSE 
          DO ll=m+2, l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          legendre_polynomials= pll
       ENDIF
    ENDIF

  END FUNCTION legendre_polynomials

  !     
  !     This function uses Newton-Raphson's 
  !     method to find the root of a function
  !     known to lie between x1 and x2. The root is returned as rtnewt
  !     and is refined till its accuracy is xacc. You need to provide
  !     this value when calling the function.
  !     You must also provide an external function called funcd which
  !     calculates both the function f and its derivative df.
  !
  REAL(DP) FUNCTION rtnewt(funcd,x1,x2,xacc)
    IMPLICIT NONE
    INTEGER :: jmax, j
    REAL(DP) :: x1,x2,xacc
    PARAMETER (jmax=20)
    REAL(DP) :: df,dx,f
    INTERFACE
       SUBROUTINE funcd(x,y,z)
         IMPLICIT NONE
         DOUBLE PRECISION :: x, y, z 

       END SUBROUTINE  funcd
    END INTERFACE


    rtnewt=.5*(x1+x2)
    DO j=1,jmax
       CALL funcd(rtnewt,f,df)
       dx=f/df
       rtnewt=rtnewt-dx
       IF ((x1-rtnewt)*(rtnewt-x2) < 0.) THEN
          WRITE(6,*) 'Error in rtnewt, jumped out of brackets'
          STOP
       ENDIF
       IF (ABS (dx) < xacc) RETURN
    ENDDO

  END FUNCTION rtnewt

  !     
  !     This function uses the bisection  method to find the 
  !     root of a function
  !     known to lie between x1 and x2. The root is returned as rtbis
  !     and is refined till its accuracy is xacc. 
  !     You must provide an external function called func.
  !
  REAL(DP) FUNCTION rtbis(func,x1,x2,xacc)
    IMPLICIT NONE
    INTEGER :: jmax, j
    REAL(DP) :: x1,x2,xacc
    PARAMETER (jmax=40)
    REAL(DP) :: dx,f,fmid,xmid

    INTERFACE
       DOUBLE PRECISION FUNCTION  func(x)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: x

       END FUNCTION func
    END INTERFACE

    fmid=func(x2)
    f=func(x1)
    IF(f*fmid >= 0.) STOP 'root must be bracketed in rtbis'
    IF (f < 0.) THEN
       rtbis=x1
       dx=x2-x1
    ELSE
       rtbis=x2
       dx=x1-x2
    ENDIF
    DO j=1,jmax
       dx=dx*.5
       xmid=rtbis+dx
       fmid=func(xmid)
       IF (fmid <= 0.) rtbis=xmid
       IF (ABS(dx) < xacc .or. fmid == 0.) RETURN
    ENDDO

  END FUNCTION rtbis
  !     
  !     This function uses the secant method 
  !     to find the root of a function
  !     known to lie between x1 and x2. The root is returned as rtsec
  !     and is refined till its accuracy is xacc
  !     You must provide an external function called func.
  !
  REAL(DP) FUNCTION rtsec(func,x1,x2,xacc)
    IMPLICIT NONE
    INTEGER :: maxit, j
    real(DP) :: x1,x2,xacc
    PARAMETER (maxit=30)
    REAL(DP) :: dx,f,fl,swap,xl
    INTERFACE
       DOUBLE PRECISION FUNCTION  func(x)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: x

       END FUNCTION func
    END INTERFACE


    fl=func(x1)
    f=func(x2)
    IF(ABS(fl) < ABS(f)) THEN
       rtsec=x1
       xl=x2
       swap=fl
       fl=f
       f=swap
    ELSE
       xl=x1
       rtsec=x2
    ENDIF
    DO j=1,maxit
       dx=(xl-rtsec)*f/(f-fl)
       xl=rtsec
       fl=f
       rtsec=rtsec+dx
       f=func(rtsec)
       IF((ABS(dx) < xacc).or.(f == 0.)) RETURN
    ENDDO

  END FUNCTION rtsec


  !     
  !     This function uses Brent's method to find the root of a function
  !     known to lie between x1 and x2. The root is returned as zbrent
  !     and is refined till its accuracy is tol
  !     You must provide an external function called func.
  !
  REAL(DP) FUNCTION zbrent(func,x1,x2,tol)
    IMPLICIT NONE
    REAL(DP) :: tol, x1, x2
    REAL(DP), PARAMETER :: EPS= epsilon(x1)
    INTEGER, PARAMETER :: ITMAX=100
    INTEGER :: iter
    REAL(DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    INTERFACE
       DOUBLE PRECISION FUNCTION  func(x)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: x

       END FUNCTION func
    END INTERFACE

    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    !      IF((fa > 0.0.AND.fb > 0.0).OR.(fa < 0.0.AND.fb < 0.0)) THEN
    !          WRITE(*,*) 'root must be bracketed for zbrent' ; STOP 
    !      ENDIF
    c=b
    fc=fb
    DO iter=1,ITMAX
       IF ((fb > 0.0.AND.fc > 0.).OR.(fb < 0.0.AND.fc < 0.0))THEN
          c=a
          fc=fa
          d=b-a
          e=d
       ENDIF
       IF(ABS(fc) < ABS(fb)) THEN
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       ENDIF
       tol1=2.*eps*ABS(b)+0.5*tol
       xm=.5*(c-b)
       IF((ABS(xm) <= tol1) .OR. (fb == 0.))THEN
          zbrent=b
          RETURN
       ENDIF
       IF((ABS(e) >= tol1) .AND.( ABS(fa) > ABS(fb)) ) THEN
          s=fb/fa
          IF(a == c) THEN
             p=2.*xm*s
             q=1.-s
          ELSE
             q=fa/fc
             r=fb/fc
             p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
             q=(q-1.)*(r-1.)*(s-1.)
          ENDIF
          IF(p > 0.) q=-q
          p=ABS(p)
          IF(2.*p < MIN(3.*xm*q-ABS(tol1*q),ABS(e*q))) THEN
             e=d
             d=p/q
          ELSE
             d=xm
             e=d
          ENDIF
       ELSE
          d=xm
          e=d
       ENDIF
       a=b
       fa=fb
       b=b+MERGE(d,SIGN(tol1,xm), ABS(d) > tol1 )
       fb=func(b)
    ENDDO
    WRITE (*,*) 'zbrent exceeding maximum iterations'; STOP
    zbrent=b

  END FUNCTION zbrent

  !
  !     The function
  !           ran0()
  !     is an "Minimal" random number generator of Park and Miller
  !     (see Numerical recipe page 279). Set or reset the input value
  !     idum to any integer value (except the unlikely value MASK)
  !     to initialize the sequence; idum must not be altered between
  !     calls for sucessive deviates in a sequence.
  !     The function returns a uniform deviate between 0.0 and 1.0.
  !


  REAL(DP) FUNCTION ran0(idum)
    IMPLICIT NONE
    INTEGER :: idum,ia,im,iq,ir,mask,k
    REAL(DP) :: am
    PARAMETER (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,mask=123459876)

    idum=ieor(idum,MASK)
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    IF  (idum < 0) idum=idum+IM
    ran0=am*idum
    idum=ieor(idum,MASK)

  END FUNCTION ran0


  !
  !      The function
  !           ran1()
  !      is an "Minimal" random number generator of Park and Miller
  !      (see Numerical recipe page 280) with Bays-Durham shuffle and
  !      added safeguards. Call with idum a negative integer to initialize;
  !      thereafter, do not alter idum between sucessive deviates in a
  !      sequence. RNMX should approximate the largest floating point value
  !      that is less than 1.
  !      The function returns a uniform deviate between 0.0 and 1.0
  !      (exclusive of end-point values).
  !



  REAL(DP) FUNCTION ran1(idum)
    IMPLICIT NONE
    INTEGER :: idum,ia,im,iq,ir,ntab,ndiv
    REAL(DP) :: am,eps,rnmx
    PARAMETER (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836, &
         ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
    INTEGER :: j,k,iv(ntab),iy
    SAVE iv,iy
    DATA iv /NTAB*0/, iy /0/

    IF ((idum == 0).OR.(iy==0)) THEN
       idum=MAX(-idum,1)
       DO j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          IF (idum < 0) idum=idum+im
          IF (j <= ntab) iv(j)=idum
       ENDDO
       iy=iv(1)
    ENDIF
    k=idum/IQ
    idum=ia*(idum-k*iq)-ir*k
    if (idum < 0) idum=idum+IM
    j=1+iy/ndiv
    iy=iv(j)
    iv(j)=idum
    ran1=MIN(am*iy,rnmx)

  END FUNCTION ran1


  REAL(DP) FUNCTION ran(idum)
    IMPLICIT NONE
    INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B), INTENT(INOUT) :: idum 
    INTEGER :: ia,im,iq,ir
    REAL(DP), SAVE :: am
    PARAMETER (ia=16807,im=2147483647,iq=127773,ir=2836)
    INTEGER(K4B), SAVE :: ix=-1, iy=-1, k 

    IF ((idum <= 0).OR.(iy < 0)) THEN
       am=NEAREST(1.,-1.)/im
       iy=IOR(IEOR(888889999,ABS(idum)),1)
       ix=IEOR(777755555,ABS(idum))
       idum=ABS(idum)+1
    ENDIF
    ix=IEOR(ix,ISHFT(ix,13))
    ix=IEOR(ix,ISHFT(ix,-17))
    ix=IEOR(ix,ISHFT(ix,5))         
    k=iy/iq
    iy=ia*(iy-k*iq)-ir*k
    IF ( iy < 0 ) iy=iy+im
    ran=am*IOR(IAND(im,IEOR(ix,iy)),1)

  END FUNCTION ran

  !
  !     The function 
  !         ran2
  !     is a long periode (> 2 x 10^18) random number generator of 
  !     L'Ecuyer and Bays-Durham shuffle and added safeguards.
  !     Call with idum a negative integer to initialize; thereafter,
  !     do not alter idum between sucessive deviates in a
  !     sequence. RNMX should approximate the largest floating point value
  !     that is less than 1.
  !     The function returns a uniform deviate between 0.0 and 1.0
  !     (exclusive of end-point values).
  !

  REAL(DP)  FUNCTION ran2(idum)
    IMPLICIT NONE
    INTEGER :: idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
    REAL(DP) :: am,eps,rnmx
    PARAMETER (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1, &
         ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791, &
         ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
    INTEGER :: idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    IF  (idum <= 0) then
       idum=MAX(-idum,1)
       idum2=idum
       DO j=NTAB+8,1,-1
          k=idum/IQ1
          idum=ia1*(idum-k*iq1)-k*ir1
          IF  (idum < 0) idum=idum+IM1
          IF (J <= NTAB) IV(J)=IDUM
       ENDDO
       IY=IV(1)
    ENDIF
    k=idum/iq1
    idum=ia1*(idum-k*iq1)-k*ir1
    if (idum < 0) idum=idum+im1
    k=idum2/iq2
    idum2=ia2*(idum2-k*iq2)-k*ir2
    IF (idum2 < 0) idum2=idum2+im2
    j=1+iy/ndiv
    iy=iv(j)-idum2
    iv(j)=idum
    IF (iy < 1)iy=iy+imm1
    ran2=MIN(am*iy,rnmx)

  END FUNCTION ran2


  !
  !     The function
  !        ran3
  !     returns a uniform random number deviate between 0.0 and 1.0. Set
  !     the idum to any negative value to initialize or reinitialize the
  !     sequence. Any large MBIG, and any small (but still large) MSEED
  !     can be substituted for the present values. 
  !

  REAL(DP) FUNCTION ran3(idum)
    IMPLICIT NONE
    INTEGER :: idum
    INTEGER :: mbig,mseed,mz
    REAL(DP) ::  fac
    PARAMETER (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
    INTEGER :: i,iff,ii,inext,inextp,k
    INTEGER :: mj,mk,ma(55)
    SAVE iff,inext,inextp,ma
    DATA iff /0/

    IF ( (idum < 0) .or. (iff == 0) ) THEN
       iff=1
       mj=mseed-IABS(idum)
       mj=MOD(mj,mbig)
       ma(55)=mj
       mk=1
       DO i=1,54
          ii=MOD(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          IF(mk < mz)mk=mk+mbig
          mj=ma(ii)
       ENDDO
       DO k=1,4
          DO i=1,55
             ma(i)=ma(i)-ma(1+MOD(i+30,55))
             IF (ma(i) < mz)ma(i)=ma(i)+mbig
          ENDDO
       ENDDO
       inext=0
       inextp=31
       idum=1
    ENDIF
    inext=inext+1
    IF (inext == 56) inext=1
    inextp=inextp+1
    IF (inextp == 56) inextp=1
    mj=ma(inext)-ma(inextp)
    IF (mj < mz) mj=mj+mbig
    ma(inext)=mj
    ran3=mj*fac

  END FUNCTION ran3

END MODULE F90library






