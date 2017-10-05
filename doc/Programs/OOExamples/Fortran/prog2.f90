MODULE shape_mod
PRIVATE    ! hide the type-bound procedure implementation procedures
PUBLIC :: shape, constructor   ! allow access to shape & constructor procedure
TYPE shape
   PRIVATE               ! hide the underlying details
        INTEGER :: color
        LOGICAL :: filled
        INTEGER :: x
        INTEGER :: y
   CONTAINS
   PRIVATE                 ! hide the type bound procedures by default
   PROCEDURE :: initShape  ! private type-bound procedure
   PROCEDURE, PUBLIC :: isFilled ! allow access to isFilled type-bound procedure
   PROCEDURE, PUBLIC :: PRINT ! allow access to print type-bound procedure
END TYPE shape
CONTAINS

LOGICAL FUNCTION isFilled(this)
CLASS(shape) :: this

isFilled = this%filled

END FUNCTION
 
FUNCTION constructor(color, filled, x, y)
  TYPE(shape) :: constructor
  INTEGER :: color
  LOGICAL :: filled
  INTEGER :: x
  INTEGER :: y
  CALL constructor%initShape(color, filled, x, y)
END FUNCTION constructor
SUBROUTINE initShape(this, color, filled, x, y)
  ! initialize shape objects
  CLASS(shape) :: this
  INTEGER :: color
  LOGICAL :: filled
  INTEGER :: x
  INTEGER :: y

  this%color = color
  this%filled = filled
  this%x = x
  this%y = y
END SUBROUTINE initShape

SUBROUTINE PRINT(this)
  CLASS(shape) :: this
  PRINT *, this%color, this%filled, this%x, this%y

END SUBROUTINE PRINT
END MODULE

PROGRAM shape_prg
USE shape_mod
TYPE(shape) :: sh
LOGICAL filled
sh = constructor(5, .TRUE., 100, 200)
CALL sh%PRINT()
END
