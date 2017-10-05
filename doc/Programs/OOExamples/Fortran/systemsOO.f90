!   This module contains some generic physical constants, such as hbarc, pi, masses of particles etc
MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  REAL(DP), PARAMETER, PUBLIC ::  bohr_r=0.05291772108_dp  ! Bohr radius in nm
  REAL(DP), PARAMETER, PUBLIC :: hbarc = 197.326968_dp    !  2006 value
  REAL(DP), PUBLIC, PARAMETER :: pi = 3.141592741012573_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_2 = 1.570796370506287_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_4 = 0.7853981852531433_dp
END MODULE constants

!   This module and its submodules allows the user to inherent quantum
!   numbers from lower and more basic objects
MODULE system
  USE constants
  ! This is the basis type used, and contains all quantum numbers necessary
  ! for fermions in one dimension
  TYPE SpQuantumNumbers
     !  n is the principal quantum number taken as number of nodes-1
     !  s is the spin and ms is the spin projection, and parity is obvious 
     INTEGER :: TotalOrbits
     INTEGER, DIMENSION(:), POINTER :: n, s, ms, parity
     REAL(DP), DIMENSION(:), POINTER :: masses, energy
  END TYPE SpQuantumNumbers
  ! We add then quantum numbers appropriate for two-dimensional systems,
  ! suitable  for electrons in quantum dots for example
  ! Use as TYPE(TwoDim) :: qdelectrons
  ! n => qdelectrons%n     
  TYPE, EXTENDS(SpQuantumNumbers) :: TwoDim
     INTEGER, DIMENSION(:), POINTER :: ml
  END TYPE TwoDim
  ! Then we extends to three dimensions, suitable for atoms and electrons in
  ! 3d traps
  ! Use as TYPE(ThreeDim) :: electrons
  ! n => electrons%n
  TYPE, EXTENDS(TwoDim) :: ThreeDim
     INTEGER, DIMENSION(:), POINTER :: l, j, mj
  END TYPE ThreeDim
  ! Then we extends to nucleons and isobars, note that the masses are in 
  ! SpQuantumNumbers. We add isospin and its projections
  ! Use as TYPE(nucleons) :: protons
  ! n => protons%n
  TYPE, EXTENDS(ThreeDim) :: nucleons
     INTEGER, DIMENSION(:), POINTER :: t, tz
  END TYPE nucleons
  ! Finally we allow for studies of hypernuclei, adding strangeness
  ! Use as TYPE(hyperons) :: sigma
  ! n => sigma%n; s => sigma%strange
  TYPE, EXTENDS(nucleons) :: hyperons
     INTEGER, DIMENSION(:), POINTER :: strange
  END TYPE hyperons
END MODULE system


