
!   Program to solve the two-dimensional Ising model 
!   with zero external field.
!   The coupling constant J = 1
!   Boltzmann's constant = 1, temperature has thus dimension energy
!   Metropolis sampling is used. Periodic boundary conditions.

PROGRAM ising2dim
  USE constants
  USE f90library
  INTEGER :: X,Y, idum, n_spins, mcs, de, t, cycles, tstep, i, j 
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spin_matrix
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: m_matrix, phi
  REAL(dp) :: w(-8:8), average(5)
  REAL(DP) :: initial_temp, final_temp, E, M, temp_step, temperature, r,s,v,p

  CALL read_input(n_spins, mcs, initial_temp, final_temp, temp_step)
  ALLOCATE(spin_matrix(n_spins,n_spins))
  ALLOCATE(m_matrix(mcs),phi(mcs))
  ! We assume that we start with the lowest temperature and that all
  ! spins are ordered. These spins are then updated with the latest value
  ! for every temperature in the loop over temperatures. The next temp step
  ! receives as input the previous spin matrix
  DO y =1, n_spins
     DO x= 1, n_spins
        spin_matrix(x,y) = 1
!        IF ( ran1(idum) <= 0.5 )         spin_matrix(x,y) = -1 
     ENDDO
  ENDDO
  ! Find temperature step
!  tstep = INT((final_temp-initial_temp)/temp_step)
  tstep = 1; temperature = initial_temp
  ! random starting point
  idum = -1 
  m_matrix = 0.0_dp; phi = 0.0_dp
  OPEN(UNIT=6,FILE='mc_calc.dat');   OPEN(UNIT=7,FILE='phi.dat')
  DO t = 1, tstep
     !   initialise energy, magnetization and array for energy changes 
     E = 0.0_dp; M = 0.0_dp; w = 0.0_dp
     DO de =-8, 8, 4 
        w(de) = EXP(-de/temperature)
     ENDDO
     ! initialise array for expectation values
     average = 0.0_dp
     CALL initialize(n_spins, temperature, spin_matrix, E, M)
     ! start Monte Carlo computation
     DO cycles = 1, mcs
        CALL Metropolis(n_spins, idum, spin_matrix, E, M, w)
        ! update expectation values
        average(1) = average(1) + E;     average(2) = average(2) + E*E
        average(3) = average(3)+M;     average(4) = average(4)+M*M 
        average(5) = average(5)+ABS(M)
!        m_matrix(cycles) = abs(m)/(n_spins**2)
        m_matrix(cycles) = average(5)/(n_spins**2)/cycles
        IF ( cycles <= mcs*0.1) write(7,*) cycles, average(5)/(n_spins**2)/cycles, &
                                                   average(3)/(n_spins**2)/cycles
     ENDDO
     DO i = 1, mcs
        r = 1.0_dp/(mcs-i)
        s = 0.0_dp; v = 0.0_dp; p= 0.0_dp
        DO j = 1, mcs-i
           p = p+ m_matrix(j)*m_matrix(j+i)
           s = s+ m_matrix(j)
           v = v+ m_matrix(j+i)
        ENDDO
        phi(i)  = r*p-r*r*s*v
!        write(7,*) i, phi(i)/phi(1)
     ENDDO

     ! print results
!     CALL output(n_spins, mcs, temperature, average)
     ! update temp
     temperature = temperature + temp_step
  ENDDO
  DEALLOCATE(spin_matrix)
  DEALLOCATE(m_matrix)

END PROGRAM ising2dim
!
! read in input data
!
SUBROUTINE  read_input(n_spins, mcs, initial_temp, final_temp, temp_step)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(inout) :: n_spins, mcs
  REAL(dp), INTENT(inout) :: initial_temp, final_temp, temp_step

  WRITE(*,*) 'Number of Monte Carlo trials =' 
  READ(*,*) mcs
  WRITE(*,*) 'Lattice size or number of spins (x and y equal) ='
  READ(*,*) n_spins
  WRITE(*,*) 'Initial temperature with dimension energy='
  READ(*,*) initial_temp
  WRITE(*,*) 'Final temperature with dimension energy='
  READ(*,*) final_temp
  WRITE(*,*) 'Temperature step with dimension energy='
  READ(*,*) temp_step

END SUBROUTINE read_input
!
! function to initialise energy, spin matrix and magnetization
!
SUBROUTINE  initialize(n_spins, temperature, spin_matrix, E, M)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n_spins
  INTEGER, INTENT(inout) :: spin_matrix(n_spins, n_spins)
  REAL(dp), INTENT(in) :: temperature
  REAL(dp), INTENT(inout) :: E, M
  INTEGER :: x, y, right, left, up, down

  ! setup initial energy and magnetization
  DO y =1, n_spins
     DO x= 1, n_spins
        right = x+1 ; IF(x == n_spins  ) right = 1
        left = x-1 ; IF(x == 1  ) left = n_spins 
        up = y+1 ; IF(y == n_spins  ) up = 1 
        down = y-1 ; IF(y == 1  ) down = n_spins 
        e= e-spin_matrix(x,y)*(spin_matrix(right,y)+&
             spin_matrix(left,y)+spin_matrix(x,up)+ &
             spin_matrix(x,down) )
        m = m + spin_matrix(x,y) 
     ENDDO
  ENDDO
  e = e*0.5_dp

END SUBROUTINE  initialize
!
!  Here we perform the sampling
!
SUBROUTINE metropolis(n_spins, idum, spin_matrix, e, m, w)
  USE constants
  USE f90library
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n_spins, idum
  INTEGER :: x, y, ix, iy,  deltae, right, left, up, down
  INTEGER, INTENT(inout) :: spin_matrix(n_spins, n_spins)
  REAL(dp) :: e, m, w(-8:8)

  ! loop over all spins
  DO y =1, n_spins
     DO x= 1, n_spins
        ix = INT(ran1(idum)*n_spins)+1
        iy = INT(ran1(idum)*n_spins)+1
        right = ix+1 ; IF(ix == n_spins  ) right = 1 
        left = ix-1 ; IF(ix == 1  ) left = n_spins 
        up = iy+1 ; IF(iy == n_spins  ) up = 1 
        down = iy-1 ; IF(iy == 1  ) down = n_spins 
        deltae = 2*spin_matrix(ix,iy)*(spin_matrix(right,iy)+&
             spin_matrix(left,iy)+spin_matrix(ix,up)+ &
             spin_matrix(ix,down) )
        IF ( ran1(idum) <= w(deltae) ) THEN
           spin_matrix(ix,iy) = -spin_matrix(ix,iy)  !flip one spin and accept new spin config
           m = m+2*spin_matrix(ix,iy)
           e = e+deltaE
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE metropolis
!
!  Write out <E>, <M>, Heat capacity and susceptibility.
!
SUBROUTINE output(n_spins, mcs, temperature, average)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n_spins, mcs
  REAL(dp), INTENT(in) :: temperature, average(5)
  REAL(dp) :: norm, Eaverage, E2average, Maverage, M2average, Mabsaverage, &
       Evariance, Mvariance

  norm = 1.0_dp/mcs  ! divided by total number of cycles 
  Eaverage = average(1)*norm
  E2average = average(2)*norm
  Maverage = average(3)*norm
  M2average = average(4)*norm
  Mabsaverage = average(5)*norm
  ! all expectation values are per spin, divide by 1/n_spins/n_spins
  Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins
  Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins
  WRITE(6,'(6F12.6)')temperature, Eaverage/n_spins/n_spins, Evariance/temperature/temperature, &
       Maverage/n_spins/n_spins, Mvariance/temperature, Mabsaverage/n_spins/n_spins

END SUBROUTINE output

