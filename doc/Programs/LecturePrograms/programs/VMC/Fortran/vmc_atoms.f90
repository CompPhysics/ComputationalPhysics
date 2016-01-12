! In this module we define the variables needed for the numerical
! calculation of the second derivative of the kinetic energy
! h is the step length while h2 is the squared inverse step length

MODULE derivative_step
USE constants
    REAL(DP), PUBLIC  :: h = 0.001_dp
    REAL(DP), PUBLIC  :: h2 = 1000000_dp
END MODULE derivative_step


! Begin of main program   

PROGRAM vmc_helium
USE constants
IMPLICIT NONE
INTEGER :: number_cycles, max_variations, thermalization, & 
           dimension, number_particles
REAL(DP) ::  step_length
REAL(DP), ALLOCATABLE, DIMENSION(:) :: cumulative_e, cumulative_e2

    !   Read in data 
    CALL initialise(dimension, number_particles, max_variations, &
                    number_cycles, thermalization, step_length) 
    ! allocate array for expectation values
    ALLOCATE(cumulative_e(max_variations), cumulative_e2(max_variations))
    cumulative_e=0.0_dp ;  cumulative_e2=0.0_dp
    !  Do the mc sampling  
    CALL mc_sampling(dimension, number_particles, max_variations, &
                     thermalization, number_cycles, step_length, &
                     cumulative_e, cumulative_e2)
    ! Print out results  
    CALL output(max_variations, number_cycles, cumulative_e, cumulative_e2)
END PROGRAM vmc_helium


! Monte Carlo sampling with the Metropolis algorithm  

SUBROUTINE mc_sampling(dimension, number_particles, max_variations,  &
                       thermalization, number_cycles, step_length,   &
                       cumulative_e, cumulative_e2)
USE constants
USE f90library
IMPLICIT NONE
REAL(DP), DIMENSION(max_variations), INTENT(INOUT) ::  cumulative_e, cumulative_e2
INTEGER, INTENT(IN) :: dimension, number_particles, max_variations, &
                       number_cycles, thermalization
REAL(DP), INTENT(IN) ::   step_length
INTEGER ::  cycles, variate, accept, dim, i, j
INTEGER idum
REAL(DP) :: wfnew, wfold, alpha, energy, energy2, &
                      delta_e, wave_function, local_energy
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: r_old, r_new
    alpha = 0.0_dp
    idum=-1
    ! allocate matrices which contain the position of the particles  
    ALLOCATE(r_new(dimension, number_particles))
    ALLOCATE(r_old(dimension, number_particles))
    ! loop over variational parameters  
    DO variate=1,  max_variations
       ! initialisations of variational parameters and energies 
       alpha = alpha + 0.1_dp  
       energy = 0.0_dp; energy2 = 0.0_dp; accept =0; delta_e=0.0_dp;
       !  initial trial position, note calling with alpha 
       !  and in three dimensions 
       DO i = 1,  number_particles
          DO j=1, dimension
             r_old(j,i) = step_length*(ran1(idum)-0.5_dp)/alpha
          ENDDO
       ENDDO
       wfold = wave_function(r_old, alpha, dimension, number_particles)
       ! loop over monte carlo cycles 
       DO cycles = 1, number_cycles+thermalization
          ! new position 
          DO i = 1,  number_particles
             DO j=1, dimension
                r_new(j,i)=r_old(j,i) + step_length*(ran1(idum)-0.5_dp)/alpha
             ENDDO
          ENDDO
          wfnew = wave_function(r_new, alpha, dimension, number_particles)
          ! Metropolis test 
          IF (ran1(idum) <= wfnew*wfnew/wfold/wfold ) THEN
             r_old=r_new
             wfold = wfnew
             accept = accept+1
          ENDIF
          ! compute local energy  
          IF ( cycles > thermalization ) THEN
             delta_e = local_energy(r_old, alpha, wfold, &
                                    dimension, number_particles)
             ! update energies  
             energy = energy+delta_e
             energy2 = energy2+delta_e*delta_e
          ENDIF
       ENDDO   ! end of loop over MC trials   
       WRITE(*,*)'variational parameter and accepted steps', alpha, accept
       ! update the energy average and its squared 
       cumulative_e(variate) = energy/number_cycles
       cumulative_e2(variate) = energy2/number_cycles
    ENDDO    ! end of loop over variational  steps 
    DEALLOCATE ( r_old, r_new)
END SUBROUTINE  mc_sampling  ! end mc_sampling function  


! Function to compute the squared wave function  

REAL(DP) FUNCTION wave_function(r, alpha,dimension, number_particles)
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: dimension, number_particles
REAL(DP), DIMENSION(dimension,number_particles), INTENT(IN) :: r
REAL(DP), INTENT(IN) :: alpha
INTEGER ::  i, j
REAL(DP) :: argument, r_single_particle
    argument = 0.0_dp
    DO i = 1, number_particles 
       r_single_particle = 0.0_dp
       DO j = 1, dimension 
          r_single_particle = r_single_particle+ r(j,i)*r(j,i)
       ENDDO
       argument = argument + SQRT(r_single_particle)
    ENDDO
    argument = argument*alpha
    wave_function = EXP(-argument) 
END FUNCTION wave_function

! Function to calculate the local energy 

REAL(DP) FUNCTION local_energy(r, alpha, wfold, dimension, &
                                       number_particles)
USE derivative_step
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: dimension, number_particles
REAL(DP), DIMENSION(dimension,number_particles), INTENT(IN) :: r
REAL(DP), INTENT(IN) :: alpha, wfold
INTEGER :: i, j , k
REAL(DP) :: e_local, wfminus, wfplus, e_kinetic, e_potential, &
                    r_12, r_single_particle, wave_function
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: r_plus, r_minus
    ! allocate matrices which contain the position of the particles  
    ! the function matrix is defined in the progam library 
    ALLOCATE(r_plus(dimension, number_particles))
    ALLOCATE(r_minus(dimension, number_particles))
    r_plus = r; r_minus = r
    ! compute the kinetic energy  
    e_kinetic = 0.0_dp
    DO i = 1, number_particles
       DO j = 1, dimension 
          r_plus(j,i) = r(j,i)+h
          r_minus(j,i) = r(j,i)-h
          wfminus = wave_function(r_minus, alpha, dimension, number_particles) 
          wfplus  = wave_function(r_plus, alpha, dimension, number_particles) 
          e_kinetic = e_kinetic-(wfminus+wfplus-2*wfold)
          r_plus(j,i) = r(j,i)
          r_minus(j,i) = r(j,i)
       ENDDO
    ENDDO
    ! include electron mass and hbar squared and divide by wave function 
    e_kinetic = 0.5_dp*h2*e_kinetic/wfold
    ! compute the potential energy 
    e_potential = 0.0_dp
    ! contribution from electron-proton potential  
    DO i = 1, number_particles
       r_single_particle = 0.0_dp
       DO j = 1, dimension
          r_single_particle = r_single_particle+r(j,i)*r(j,i)
       ENDDO
       e_potential = e_potential-number_particles/SQRT(r_single_particle)
    ENDDO
    ! contribution from electron-electron potential  
    DO i = 1, number_particles-1
       DO J= i+1, number_particles
	  r_12 = 0
          DO k = 1, dimension
             r_12 = r_12+(r(k,i)-r(k,j))*(r(k,i)-r(k,j))
          ENDDO
	  e_potential =  e_potential+1/SQRT(r_12)          
       ENDDO
    ENDDO
    local_energy = e_potential+e_kinetic
    DEALLOCATE (r_plus, r_minus)
END FUNCTION local_energy



SUBROUTINE initialise(dimension, number_particles, max_variations, &
                      number_cycles, thermalization, step_length) 
USE constants
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: dimension, number_particles, max_variations, &
                          number_cycles, thermalization
REAL(DP), INTENT(INOUT) ::   step_length
    WRITE(*,*)'number of particles = '
    READ(*,*) number_particles
    WRITE(*,*)'dimensionality = '
    READ(*,*) dimension
    WRITE(*,*)'maximum variational parameters = '
    READ(*,*) max_variations
    WRITE(*,*)'# Thermalization  steps= '
    READ(*,*) thermalization
    WRITE(*,*)'# MC steps= '
    READ(*,*)number_cycles
    WRITE(*,*)'# step length= '
    READ(*,*)step_length
END SUBROUTINE initialise



SUBROUTINE output(max_variations, number_cycles, cumulative_e, cumulative_e2)
USE constants
IMPLICIT NONE
REAL(DP), DIMENSION(max_variations), INTENT(IN) ::  cumulative_e, &
                                                            cumulative_e2
INTEGER, INTENT(IN) :: max_variations, number_cycles
INTEGER :: i
REAL(DP) :: alpha, variance, error
    ! open output file 
    OPEN(UNIT=6,FILE='output_file')
    alpha = 0
    DO i=1, max_variations
       alpha = alpha+0.1_dp  
       variance = cumulative_e2(i)-cumulative_e(i)*cumulative_e(i)
       error=SQRT(variance/number_cycles)
       WRITE(6,*)alpha,cumulative_e(i),variance, error
    ENDDO
    CLOSE(6)
END SUBROUTINE  output

