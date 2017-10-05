MODULE class_Circle
  IMPLICIT NONE
  PRIVATE
  REAL :: pi = 3.1415926535897931d0 ! Class-wide private constant

  TYPE, PUBLIC :: Circle
     REAL :: radius
   CONTAINS
     PROCEDURE :: area => circle_area
     PROCEDURE :: PRINT => circle_print
  END TYPE Circle
CONTAINS
  FUNCTION circle_area(this) RESULT(area)
    CLASS(Circle), INTENT(in) :: this
    REAL :: area
    area = pi * this%radius**2
  END FUNCTION circle_area

  SUBROUTINE circle_print(this)
    CLASS(Circle), INTENT(in) :: this
    REAL :: area
    area = this%area()  ! Call the type-bound function
    PRINT *, 'Circle: r = ', this%radius, ' area = ', area
  END SUBROUTINE circle_print
END MODULE class_Circle


PROGRAM circle_test
  USE class_Circle
  IMPLICIT NONE

  TYPE(Circle) :: c     ! Declare a variable of type Circle.
  c = Circle(1.5)       ! Use the implicit constructor, radius = 1.5.
  CALL c%PRINT          ! Call the type-bound subroutine
END PROGRAM circle_test
