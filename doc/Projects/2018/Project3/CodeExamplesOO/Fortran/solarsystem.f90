program solarsystem
  use solver_class

  implicit none

  integer,parameter::numbodies=10 
  real(8)::Tmax=250,h
  integer::Num_Steps=10000,m
  real(8),dimension(10)::masses
  real(8),dimension(10,3)::position
  real(8),dimension(10,3)::velocity
  real(8),dimension(numbodies,numbodies,4)::relposition,updaterel
  real(8),dimension(numbodies,3)::relforce
  real(8),dimension(numbodies,3)::updatedforce
  real(8),dimension(numbodies+1)::kinetic,potential,angular

  type(solver)::solar
  
   
  solar%mass(1)=1.d0
  solar%mass(2)=3.d-6
  solar%mass(3)=9.5d-4
  solar%mass(4)=3.227d-7 
  solar%mass(5)=2.4478d-6
  solar%mass(6)=2.85812d-4
  solar%mass(7)=1.66012d-7 
  solar%mass(8)=4.36576d-5
  solar%mass(9)=5.1503d-5
  solar%mass(10)=6.583d-9 
  solar%position(1,1)=3.1647d-3
  solar%position(1,2)=4.4307d-3
  solar%position(1,3)=-1.51447d-4
  solar%velocity(1,1)=365.25*(-3.37957d-6)
  solar%velocity(1,2)=365.25*(6.606862d-6)
  solar%velocity(1,3)=365.25*(7.32297d-8)
  solar%position(2,1)=-9.8825d-1
  solar%position(2,2)=8.49978d-2
  solar%position(2,3)=-1.5199728d-4
  solar%velocity(2,1)=365.25*(-1.68024d-3)
  solar%velocity(2,2)=365.25*(-1.719988d-2)
  solar%velocity(2,3)=365.25*(4.34984d-7)
  solar%position(3,1)=-5.23294d0
  solar%position(3,2)=-1.52515d0
  solar%position(3,3)=1.233648d-1
  solar%velocity(3,1)=365.25*(2.022596d-3)
  solar%velocity(3,2)=365.25*(-6.88771645d-3)
  solar%velocity(3,3)=365.25*(-1.6694179d-5)
  solar%position(4,1)=7.78069d-1
  solar%position(4,2)=1.2797d0
  solar%position(4,3)=7.555377d-3
  solar%velocity(4,1)=365.25*(-1.143145d-2)
  solar%velocity(4,2)=365.25*(8.4664712d-3)
  solar%velocity(4,3)=365.25*(4.5782387d-4)
  solar%position(5,1)=-7.02894d-1
  solar%position(5,2)=1.359581d-1
  solar%position(5,3)=4.239547d-2
  solar%velocity(5,1)=365.25*(-3.8130624d-3)
  solar%velocity(5,2)=365.25*(-1.9968d-2)
  solar%velocity(5,3)=365.25*(-5.4012702e-5)
  solar%position(6,1)=-1.48071d0
  solar%position(6,2)=-9.935855d0
  solar%position(6,3)=2.31688d-1
  solar%velocity(6,1)=365.25*(5.212138d-3)
  solar%velocity(6,2)=365.25*(-8.3942195d-4)
  solar%velocity(6,3)=365.25*(-1.931769d-4)
  solar%position(7,1)=2.805339d-1
  solar%position(7,2)=1.7274317d-1
  solar%position(7,3)=-1.1844519d-2
  solar%velocity(7,1)=365.25*(-2.01015d-2)
  solar%velocity(7,2)=365.25*(2.5290758d-2)
  solar%velocity(7,3)=365.25*(3.9099347d-3)
  solar%position(8,1)=1.822435d1
  solar%position(8,2)=8.083455d0
  solar%position(8,3)=-2.0607748d-1
  solar%velocity(8,1)=365.25*(-1.623364d-3)
  solar%velocity(8,2)=365.25*(3.411947d-3)
  solar%velocity(8,3)=365.25*(3.381457d-5)
  solar%position(9,1)=2.8412218d1
  solar%position(9,2)=-9.4680088d0
  solar%position(9,3)=-4.598129d-1
  solar%velocity(9,1)=365.25*(9.7114038d-4)
  solar%velocity(9,2)=365.25*(2.99682d-3)
  solar%velocity(9,3)=365.25*(-8.375523d-5)
  solar%position(10,1)=9.890335d0
  solar%position(10,2)=-3.177864195d1
  solar%position(10,3)=5.3964748d-1
  solar%velocity(10,1)=365.25*(3.06860367d-3)
  solar%velocity(10,2)=365.25*(2.905788d-4)
  solar%velocity(10,3)=365.25*(-9.08614726d-4)

  open(3,file="Earth.dat")
  open(4,file="Jupiter.dat")
  open(5,file="energy.dat")
  open(7,file="momentum.dat")
  open(8,file="Mars.dat")
  open(9,file="Venus.dat")
  open(10,file="Saturn.dat")
  open(11,file="Mercury.dat")
  open(12,file="Uranus.dat")
  open(13,file="Neptune.dat")
  open(14,file="Pluto.dat")


  h=Tmax/Num_Steps
  do m=1,Num_Steps

     call relative_position(solar,numbodies,relposition)
     call forces(solar,Numbodies,relposition,relforce)
     call calc_position(solar,numbodies,relforce,h)
     call relative_position(solar,numbodies,updaterel)
     call forces(solar,numbodies,updaterel,updatedforce)
     call calc_velocities(solar,numbodies,relforce,updatedforce,h)
     call kinetic_energy(solar,numbodies,kinetic)
     call potential_energy(solar,numbodies,updaterel,potential)
     call angular_momentum(solar,numbodies,updaterel,angular)
     
!     write(5,*), kinetic(Numbodies+1), potential(numbodies+1)
!     write(7,*), angular(Numbodies+1)
  
     write(3,*),solar%position(2,1),solar%position(2,2),solar%position(2,3)
     write(4,*),solar%position(3,1),solar%position(3,2),solar%position(3,3)
     write(8,*),solar%position(4,1),solar%position(4,2),solar%position(4,3)
     write(9,*),solar%position(5,1),solar%position(5,2),solar%position(5,3)
     write(10,*),solar%position(6,1),solar%position(6,2),solar%position(6,3)
     write(11,*),solar%position(7,1),solar%position(7,2),solar%position(7,3)
     write(12,*),solar%position(8,1),solar%position(8,2),solar%position(8,3)
     write(13,*),solar%position(9,1),solar%position(9,2),solar%position(9,3)
     write(14,*),solar%position(10,1),solar%position(10,2),solar%position(10,3)

  end do

  
  close(3)
  close(4)
  close(5)
  close(7)
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
end program solarsystem
