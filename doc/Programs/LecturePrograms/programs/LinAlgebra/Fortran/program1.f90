
!
!      This program sets up a simple 3x3 symmetric matrix  and finds its
!      determinant and inverse
!

PROGRAM matrix
  USE constants
  USE F90library
  IMPLICIT NONE
  !      The definition of the matrix, using dynamic allocation
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: a, ainv, unity
  !      the determinant
  REAL(DP) :: d
  !      The size of the matrix
  INTEGER :: n
  !      Here we set the dim n=3
  n=3
  !      Allocate now place in heap for a
  ALLOCATE ( a(n,n), ainv(n,n), unity(n,n) )
  a(1,1)=1. ; a(1,2)=3.; a(1,3)=4. 
  a(2,1)=3. ; a(2,2)=4.; a(2,3)=6. 
  a(3,1)=4. ; a(3,2)=6.; a(3,3)=8. 
  WRITE(6,*) ' The matrix before inversion'
  WRITE(6,'(3F12.6)') a
  ainv=a
  CALL matinv (ainv, n, d)
  WRITE(6,*) ' The determinant'
  WRITE(6,'(F12.6)') d
  WRITE(6,*) ' The matrix after inversion'
  WRITE(6,'(3F12.6)') ainv
  !      get the unity matrix
  unity=MATMUL(ainv,a)
  WRITE(6,*) ' The unity matrix'
  WRITE(6,'(3F12.6)') unity 
  !      deallocate all arrays
  DEALLOCATE (a, ainv, unity)

END PROGRAM matrix
