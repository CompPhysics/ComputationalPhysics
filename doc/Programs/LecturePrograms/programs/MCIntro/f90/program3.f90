
!   Note the structure of this module, it contains various
!   subroutines for initialisation of the problem, the metropolis
!   sampling itself and a function which deallocates the allocated
!   arrays. Note that all variables are PRIVATE, and belong only to
!   the module    

MODULE variables
  USE constants
  USE f90library
  INTEGER, PRIVATE :: initial_n_particles, max_time, number_cycles 
  REAL(DP), PRIVATE :: decay_probability
  INTEGER, ALLOCATABLE, DIMENSION(:), PRIVATE :: ncumulative 
CONTAINS

  SUBROUTINE initialise
    IMPLICIT NONE

    WRITE(*,*) ' Read in initial number of particles'
    READ(*,*) initial_n_particles
    WRITE(*,*) 'Read in number of Monte carlo cycles'
    READ(*,*) number_cycles
    WRITE(*,*) 'Read in decay probability'
    READ(*,*) decay_probability
    WRITE(*,*) 'Read in number of time intervals per cycle'
    READ(*,*) max_time

    ALLOCATE ( ncumulative( 0:max_time) )
    ncumulative=0

  END SUBROUTINE initialise

  !      This function performs the sampling for nuclear decays

  SUBROUTINE mc_sampling
    IMPLICIT NONE
    INTEGER :: idum, cycles, time, np, n_unstable, particle_limit

    OPEN ( UNIT=6, FILE='out.dat')
    idum=-1
    ! loop over Monte Carlo cycles
    DO cycles = 1, number_cycles   
       n_unstable = initial_n_particles
       ! accumulate the number of particles per time step per trial
       ncumulative(0) = ncumulative(0) + initial_n_particles
       ! loop over each time step 
       DO time = 1, max_time
          ! for each time step, we check the remaining particles
          particle_limit = n_unstable
          DO np= 1,  particle_limit
             IF( ran0(idum) <= decay_probability) THEN
                n_unstable=n_unstable-1
             ENDIF
          ENDDO
          ncumulative(time)=ncumulative(time)+n_unstable
       ENDDO
    ENDDO
    !      write time and mean number of unstable nuclei
    DO time = 0, max_time
       WRITE(6,*) time, ncumulative(time)/FLOAT(number_cycles)
    ENDDO
    CLOSE(6)
    DEALLOCATE ( ncumulative)

  END SUBROUTINE mc_sampling

END MODULE variables


PROGRAM nuclear_decay
  USE variables
  IMPLICIT NONE

  !   Read in data
  CALL initialise
  !   Do the mc sampling 
  CALL mc_sampling

END PROGRAM  nuclear_decay





