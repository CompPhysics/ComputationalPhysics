MODULE shape_mod
  TYPE shape
     INTEGER :: color
     LOGICAL :: filled
     INTEGER :: x
     INTEGER :: y
   CONTAINS
     PROCEDURE :: initialize
  END TYPE shape
  TYPE, EXTENDS(shape) :: rectangle
  INTEGER :: length
  INTEGER :: width
END TYPE rectangle
TYPE, EXTENDS(rectangle) :: square
END TYPE square
CONTAINS
SUBROUTINE initialize(sh, color, filled, x, y, length, width)
  ! initialize shape objects
CLASS(shape) :: sh
INTEGER :: color
LOGICAL :: filled
INTEGER :: x
INTEGER :: y
INTEGER, OPTIONAL :: length
INTEGER, OPTIONAL :: width

sh%color = color
sh%filled = filled
sh%x = x
sh%y = y
SELECT TYPE (sh)
TYPE is (shape)
   ! no further initialization required
 CLASS is (rectangle)
    ! rectangle or square specific initializations
 IF (PRESENT(length))  THEN
    sh%length = length
 ELSE
    sh%length = 0
 ENDIF
 IF (PRESENT(width)) THEN
    sh%width = width
 ELSE
    sh%width = 0
 ENDIF
 CLASS default
    ! give error for unexpected/unsupported type
 STOP 'initialize: unexpected type for sh object!'
END SELECT
END SUBROUTINE initialize
END MODULE



PROGRAM test

USE shape_mod
TYPE(shape) :: shp                       ! declare an instance of shape
CALL shp%initialize(1, .TRUE., 10, 20)   ! initialize shape

END PROGRAM test

