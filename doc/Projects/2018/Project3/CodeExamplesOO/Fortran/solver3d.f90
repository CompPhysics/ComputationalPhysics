module solver_class
  
  type solver
     integer::Bodies
     real(8),dimension(10) :: mass
     real(8),dimension(10,3) :: position
     real(8),dimension(10,3) :: velocity
  contains
    procedure ::relative_position
    procedure ::forces
    procedure ::calc_position
    procedure ::calc_velocities
    procedure ::kinetic_energy
    procedure ::potential_energy
    procedure ::angular_momentum

 end type solver

contains
  subroutine relative_position(system,numbodies,relposition)
    implicit none

    class (solver),intent(in)::system
    integer :: i,j
    integer,intent(in)::NumBodies
    real(8),dimension(Numbodies,Numbodies,4):: relposition
    
    do i=1,NumBodies
       do j =1,NumBodies
          if (i.ne.j) then
             relposition(j,i,1)=system%position(j,1)-system%position(i,1)
             relposition(j,i,2)=system%position(j,2)-system%position(i,2)
             relposition(j,i,3)=system%position(j,3)-system%position(i,3)
             relposition(j,i,4)=sqrt(relposition(j,i,1)*relposition(j,i,1)+relposition(j,i,2)*relposition(j,i,2)+relposition(j,i,3)*relposition(j,i,3))
          end if
       end do
    end do
    
    return 
  end subroutine relative_position

  subroutine forces(system,Numbodies,relposition,relforce)
    implicit none

    class(solver),intent(in)::system
    integer :: j,i
    integer,intent(in):: Numbodies
    real(8),dimension(Numbodies,3)::relforce
    real(8),dimension(Numbodies,Numbodies,4):: relposition
    real(8)::rrr,Fourpi2
    Fourpi2 = 4.d0*3.14*3.14

    
    do i=1,numbodies
       do j=1,numbodies
          relforce(i,j)=0.d0
       end do
    end do
    

    do i=1,numbodies
       do j=1,numbodies
          if(j.ne.i) then
             rrr=(relposition(j,i,4)**3.d0)
             relforce(i,1) =relforce(i,1) - Fourpi2*system%mass(j)*relposition(j,i,1)/rrr
             relforce(i,2) =relforce(i,2) - Fourpi2*system%mass(j)*relposition(j,i,2)/rrr
             relforce(i,3) =relforce(i,3) - Fourpi2*system%mass(j)*relposition(j,i,3)/rrr
          end if
       end do
    end do
  end subroutine forces


  subroutine calc_position(system,Numbodies,relforce,h)
    implicit none
    class(solver), intent(inout)::system
    integer,intent(in)::Numbodies
    integer::i,j
    real(8),dimension(Numbodies,3)::relforce
    real(8)::h

    do i=1,numbodies
       system%position(i,1)=system%position(i,1)+h*system%velocity(i,1) - (h*h/2.d0)*relforce(i,1)
       system%position(i,2)=system%position(i,2)+h*system%velocity(i,2) - (h*h/2.d0)*relforce(i,2)
       system%position(i,3)=system%position(i,3)+h*system%velocity(i,3) - (h*h/2.d0)*relforce(i,3)
    end do

  end subroutine calc_position

  subroutine calc_velocities(system,numbodies,relforce,updatedforce,h)
    implicit none
    class(solver),intent(inout)::system
    integer,intent(in)::NumBodies
    integer::i,j
    real(8)::h
    real(8),dimension(Numbodies,3)::relforce
    real(8),dimension(Numbodies,3)::updatedforce

    do i=1, numbodies
       system%velocity(i,1)=system%velocity(i,1)-h*0.5*updatedforce(i,1)-h*0.5*relforce(i,1)
       system%velocity(i,2)=system%velocity(i,2)-h*0.5*updatedforce(i,2)-h*0.5*relforce(i,2)
       system%velocity(i,3)=system%velocity(i,3)-h*0.5*updatedforce(i,3)-h*0.5*relforce(i,3)
    end do

  end subroutine calc_velocities
    
  subroutine kinetic_energy(system,numbodies,kinetic)
    implicit none
    class(solver), intent(in)::system
    integer,intent(in)::Numbodies
    integer::i
    real(8),dimension(numbodies+1)::kinetic
    real(8)::totalke

    totalke=0
    
    do i =1,numbodies
       kinetic(i)=(1.d0/2.d0)*system%mass(i)*(system%velocity(i,1)*system%velocity(i,1)+system%velocity(i,2)*system%velocity(i,2)+system%velocity(i,3)*system%velocity(1,3))
       totalke=totalke+kinetic(i)
    end do

    kinetic(numbodies+1)=totalke
  end subroutine kinetic_energy

  subroutine potential_energy(system, numbodies,relposition,potential)
    implicit none
    class(solver), intent(in)::system
    integer,intent(in)::Numbodies
    integer::i,j
    real(8),dimension(numbodies+1)::potential
    real(8),dimension(numbodies,numbodies,4)::relposition
    real(8)::totalpe

    totalpe=0
    do i=1,numbodies+1
       potential(i)=0.d0
    end do

    do i=1,numbodies
      do j=1,NumBodies
         if (i .ne. j ) then
            potential(i)=potential(i)+(4*3.14*3.14*system%mass(i)*system%mass(j)/relposition(j,i,4))
         end if
      end do
      totalpe=totalpe+potential(i)
    end do

    potential(numbodies+1)=totalpe
  end subroutine potential_energy

  subroutine angular_momentum(system,numbodies,relposition,angular)
    implicit none
    class(solver),intent(in)::system
    integer,intent(in)::Numbodies
    integer::i,j
    real(8),dimension(numbodies+1)::angular
    real(8),dimension(numbodies,numbodies,3)::relposition
    real(8)::totalang

    totalang=0

    do i=1,numbodies
       angular(i)=system%mass(i)*relposition(1,i,4)*(sqrt(system%velocity(i,1)**2.d0 + system%velocity(i,2)**2.d0 + system%velocity(i,3)**2.d0))
       totalang=totalang+angular(i)
    end do

    angular(numbodies+1)=totalang
  end subroutine angular_momentum
    
  
end module solver_class
