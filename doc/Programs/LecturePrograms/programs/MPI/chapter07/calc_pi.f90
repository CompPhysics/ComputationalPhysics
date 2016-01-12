!/////////////////////////////////////////////////////////////////////////
!//
!// MPI example program in FORTRAN 90
!//
!/////////////////////////////////////////////////////////////////////////

PROGRAM pi
   IMPLICIT NONE
   INCLUDE "mpif.h"

   INTEGER                                  :: psize    ! Number of procs
   INTEGER                                  :: my_rank  ! This process
   INTEGER                                  :: i        ! Index
   INTEGER                                  :: n = 1000 ! Used as constant
   INTEGER                                  :: rstat    ! Status variable
   INTEGER,DIMENSION(MPI_STATUS_SIZE)       :: status   ! MPI status

   DOUBLE PRECISION                         :: l_sum
   DOUBLE PRECISION                         :: g_sum
   DOUBLE PRECISION                         :: x
   DOUBLE PRECISION                         :: h


   CHARACTER (LEN=80)                       :: str      ! Text string

   !// Initialize MPI
   CALL MPI_Init(rstat)

   !// Get number of procs
   CALL MPI_Comm_size(MPI_COMM_WORLD, psize, rstat)

   !// Get process number for this process
   CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, rstat)

   !// Fortran has no default command line args into
   !// a process. Using fixed size of n

   CALL MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, rstat)

   h = 1.0/n
   l_sum = 0.

   DO i = my_rank, n
      x = (i+0.5)*h
      l_sum = l_sum + 4.0/(1.0+x*x)
   END DO

   l_sum = l_sum + h

   IF (my_rank == 0) THEN
      g_sum = l_sum
      DO i = 1, psize - 1
         CALL MPI_Recv(l_sum, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 500, &
                       MPI_COMM_WORLD, status, rstat)
         g_sum = g_sum + l_sum
      END DO
      PRINT *, 'Result = ', g_sum
   ELSE
      CALL MPI_Send(l_sum, 1, MPI_DOUBLE_PRECISION, 0, 500, MPI_COMM_WORLD)
   END IF

   CALL MPI_Finalize(rstat)

END PROGRAM pi
