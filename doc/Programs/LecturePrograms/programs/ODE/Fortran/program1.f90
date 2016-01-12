!      This program solves Newton's equation for a block
!      sliding on a horizontal frictionless surface. The block
!      is tied  to a wall with a spring, and Newton's equation
!      takes the form
!           m d^2x/dt^2 =-kx
!      with k the spring tension and m the mass of the block.
!      The angular frequency is omega^2 = k/m and we set it equal
!      1 in this example program. 
!      Newton's equation is rewritten as two coupled differential
!      equations, one for the position x  and one for the velocity v
!           dx/dt = v    and
!           dv/dt = -x   when we set k/m=1
!      We use therefore a two-dimensional array to represent x and v
!      as functions of t
!      y[0] == x
!      y[1] == v
!      dy[0]/dt = v
!      dy[1]/dt = -x
!      The derivatives are calculated by the user defined function 
!      derivatives.
!      The user has to specify the initial velocity (usually v_0=0)
!      the number of steps and the initial position. In the programme
!      below we fix the time interval [a,b] to [0,2*pi].
!
!
!
!
!  this is the number of differential equations  as a global parameter

MODULE parameters
  INTEGER, PARAMETER, PUBLIC :: number_differential_eqs =2 
END MODULE parameters
!
!     Main function begins here 
!
PROGRAM diff_solver
  USE constants
  USE parameters
  IMPLICIT NONE
  REAL(DP), DIMENSION(number_differential_eqs) :: y, dydt, yout
  REAL(DP) ::  t, h, tmax, E0, initial_x, initial_v
  INTEGER :: i, number_of_steps

  !  read in the initial position, velocity and number of steps 
  CALL initialise (initial_x, initial_v, number_of_steps)
  !  setting initial values, step size and max time tmax  
  h = 2.0_dp*acos(-1.0_dp)/FLOAT(number_of_steps)   ! the step size     
  tmax = h*number_of_steps                  ! the final time    
  y(1) = initial_x                          ! initial position  
  y(2) = initial_v                          ! initial velocity  
  t=0.0_dp                                      ! initial time      
  E0 = 0.5_dp*(y(1)**2+y(2)**2)                ! the initial total energy
  ! now we start solving the differential equations using the RK4 method 
  OPEN(6,FILE='outf.dat')
  yout = 0.0_dp; dydt = 0.0_dp
  DO WHILE (t <= tmax)
     ! initial derivatives 
     CALL derivatives(t, y, dydt)                
     ! here we call the runge-kutta method and get the new y-value in yout                     
     CALL runge_kutta_4(y, t, h,yout,dydt)
     y = yout  
     t = t + h
     ! writing time, x, v, the exact solution, and the energy difference 
     WRITE(6,'(5(E12.6,1X))') t, y(1), y(2), cos(t),0.5*(y(1)**2+y(2)**2)-E0
  ENDDO
  CLOSE (6)

END PROGRAM diff_solver
!
!   this function sets up the derivatives for this special case  
!
SUBROUTINE derivatives(t, y, dydt)
  USE constants
  USE parameters
  IMPLICIT NONE
  REAL(DP), DIMENSION(number_differential_eqs) :: y, dydt
  REAL(DP) :: t

  dydt(1)=y(2);    ! derivative of x 
  dydt(2)=-y(1); ! derivative of v 

END SUBROUTINE derivatives
!
!     The function initialise
!     Reads in from screen the air temp, the number of steps
!     final time and the initial temperature
!
SUBROUTINE initialise(initial_x, initial_v, number_of_steps)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: number_of_steps
  REAL(DP), INTENT(OUT) :: initial_x, initial_v

  WRITE(*,*) ' Read in from screen intial x, initial v, and number of steps'
  READ(*,*) initial_x, initial_v, number_of_steps

END SUBROUTINE initialise
!
!     Runge-kutta procedure 
!
SUBROUTINE runge_kutta_4(y,x,diff_eq_step,yout,dydx)
  USE constants
  USE parameters
  IMPLICIT NONE
  REAL(DP), DIMENSION(number_differential_eqs) :: yt, dyt, dym
  REAL(DP), DIMENSION(number_differential_eqs), INTENT(IN) :: y, dydx
  REAL(DP), DIMENSION(number_differential_eqs), INTENT(OUT) :: yout
  REAL(DP) :: hh, h6, xh
  REAL(DP), INTENT(IN) :: x, diff_eq_step 

  hh=diff_eq_step*0.5; h6=diff_eq_step/6. ; xh=x+hh
  !     first rk-step
  yt=y+hh*dydx
  CALL derivatives(xh,yt,dyt)
  !     second rk-step
  yt=y+hh*dyt
  CALL derivatives(xh,yt,dym)
  !     third rk-step
  yt=y+diff_eq_step*dym;  dym=dyt+dym
  CALL derivatives(x+diff_eq_step,yt,dyt)
  !     fourth rk-step
  yout=y+h6*(dydx+dyt+2.*dym)

END SUBROUTINE runge_kutta_4

