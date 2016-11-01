!   Program to solve the two-dimensional Ising model 
!   with zero external field using MPI
!   The coupling constant J = 1
!   Boltzmann's constant = 1, temperature has thus dimension energy
!   Metropolis sampling is used. Periodic boundary conditions.
!
!   The program uses a single pseudorandom number or an array of pseudorandom
!   numbers from the uniform distribution over the range 0 \leq x < 1.
!   The runtime-library implements the xorshift1024* random number
!   generator (RNG). This generator has a period of 2^1024 - 1, and when
!   using multiple threads up to 2^512 threads can each generate 2^512
!   random numbers before any aliasing occurs.
!   Note that in a multi-threaded program (e.g. using OpenMP directives),
!   each thread will have its own random number state. For details of the
!   seeding procedure, see the documentation for the RANDOM_SEED
!   intrinsic.


!    This module contains all constants and declarations 
!    of variables read in by the function read_data. These
!    variables are used by many functions.

MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
END MODULE constants

!
! definition of MPI constants
!
MODULE mpi_constants
  INCLUDE 'mpif.h'
  INTEGER, PUBLIC :: no_intervalls, my_rank, myloop_begin, myloop_end, numprocs, ierror

END MODULE mpi_constants

!
!   Main program starts here, standard Metropolis algorithm
!

PROGRAM ising2dim
  USE constants
  USE mpi_constants
  INTEGER :: n_spins, de, t, tstep, i
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: spin_matrix
  REAL(dp) :: w(-8:8), local_average(5), total_average(5)
  REAL(DP) :: initial_temp, final_temp, E, M, temp_step, temperature
  REAL(DP) :: time_start, time_end

  ! initialize mpi
  !
  CALL MPI_INIT(ierror)
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierror )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, my_rank, ierror )
  ! Only master writes and reads to file
  IF ( my_rank == 0 ) THEN
     OPEN(UNIT=6,FILE='mc_calc.dat')
     OPEN(UNIT=5,FILE='input.dat')
     CALL read_input(n_spins, mcs, initial_temp, final_temp, temp_step)
  ENDIF
  ! Broadcats to all processes common values
  CALL MPI_Bcast(mcs, 1, MPI_INT, 0, MPI_COMM_WORLD,ierror)
  CALL MPI_Bcast(n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD,ierror)
  CALL MPI_Bcast(initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
  CALL MPI_Bcast(final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
  CALL MPI_Bcast(temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
  ! Find temperature step
  tstep = INT((final_temp-initial_temp)/temp_step); temperature = initial_temp
  ! Start look over time steps
  ALLOCATE(spin_matrix(n_spins,n_spins))
    ! All spins point up, cold start, we don't use temp here
  spin_matrix = 1
  time_start = MPI_Wtime()
  DO t = 1, tstep
     !   initialise energy, magnetization and array for energy changes 
     w = 0.0_dp
     DO de =-8, 8, 4 
        w(de) = EXP(-de/temperature)
     ENDDO
     ! initialise array for expectation values
     total_average = 0.0_dp; local_average  =0.0_dp
     CALL Metropolis(mcs, n_spins, my_rank, w, local_average, spin_matrix)
     ! collect all results
     CALL  MPI_REDUCE(local_average, total_average, 5, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
     ! print results
     IF (my_rank == 0) THEN
        CALL WritetoFile(n_spins, mcs*numprocs, temperature, total_average)
     ENDIF
     ! update temp
     temperature = temperature + temp_step
  ENDDO
  DEALLOCATE(spin_matrix)
  time_end = MPI_Wtime()
  IF (my_rank == 0) THEN
     WRITE(*,*) 'Total time =', time_end-time_start,' in seconds used  on number of processes = ', numprocs
  ENDIF   
  CALL MPI_FINALIZE(ierror)

END PROGRAM ising2dim

! read in input data
!
SUBROUTINE  read_input(n_spins, mcs, initial_temp, final_temp, temp_step)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(inout) :: n_spins, mcs
  REAL(dp), INTENT(inout) :: initial_temp, final_temp, temp_step

  ! 'Number of Monte Carlo trials, lattice, init temp, final temp and temp step
  READ(5,*) mcs, n_spins, initial_temp, final_temp, temp_step

END SUBROUTINE read_input
!
! function to initialise energy and magnetization
!
SUBROUTINE  initialize(n_spins, spin_matrix, E, M)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n_spins
  INTEGER, INTENT(IN) :: spin_matrix(n_spins, n_spins)
  REAL(dp), INTENT(INOUT) :: E, M
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
SUBROUTINE  Metropolis(mcs, n_spins, my_rank, w, average, spin_matrix)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n_spins, my_rank, mcs
  INTEGER :: spins, ix, iy,  deltae, right, left, up, down, allspins, cycles
  REAL(dp) :: e, m
  REAL(dp), INTENT(IN) :: w(-8:8)
  REAL(dp), INTENT(INOUT) :: average(5)
  INTEGER, INTENT(INOUT) :: spin_matrix(n_spins,n_spins)
  REAL(DP) :: RandomNumber(3)
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  INTEGER :: idum

  ! Initialize RNG
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)-my_rank
  CALL RANDOM_SEED(PUT = seed)
  ! initialize energy, magnetic moment
  e = 0.0_dp; m = 0.0_dp
  CALL initialize(n_spins, spin_matrix, e, m)
  ! start Monte Carlo computation
  allspins = n_spins*n_spins
  DO cycles = 1, mcs
     ! loop over all spins, that is a sweep over the lattice
     DO spins =1, allspins
        CALL RANDOM_NUMBER(RandomNumber)
        ix = INT(RandomNumber(1)*n_spins)+1
        iy = INT(RandomNumber(2)*n_spins)+1
        right = ix+1 ; IF(ix == n_spins  ) right = 1 
        left = ix-1 ; IF(ix == 1  ) left = n_spins 
        up = iy+1 ; IF(iy == n_spins  ) up = 1 
        down = iy-1 ; IF(iy == 1  ) down = n_spins 
        deltae = 2*spin_matrix(ix,iy)*(spin_matrix(right,iy)+&
             spin_matrix(left,iy)+spin_matrix(ix,up)+ &
             spin_matrix(ix,down) )
        IF ( RandomNumber(3) <= w(deltae) ) THEN
           spin_matrix(ix,iy) = -spin_matrix(ix,iy)  !flip one spin and accept new spin config
           m = m+2*spin_matrix(ix,iy)
           e = e+deltaE
        ENDIF
     ENDDO
     ! update expectation values
     average(1) = average(1) + e;     average(2) = average(2) + e*e
     average(3) = average(3)+m;     average(4) = average(4)+m*m 
     average(5) = average(5)+ABS(m)
  ENDDO
  DEALLOCATE(seed)

END SUBROUTINE metropolis
!
!  Write out <E>, <M>, Heat capacity and susceptibility.
!
SUBROUTINE WritetoFile(n_spins, mcs, temperature, average)
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
  Mvariance = (M2average - Maverage*Maverage)/n_spins/n_spins
  WRITE(6,'(6F12.6)')temperature, Eaverage/n_spins/n_spins, Evariance/temperature/temperature, &
       Maverage/n_spins/n_spins, Mvariance/temperature, Mabsaverage/n_spins/n_spins

END SUBROUTINE WritetoFile







