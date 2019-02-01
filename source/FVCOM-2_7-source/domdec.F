!==============================================================================|
!  DECOMPOSE THE DOMAIN BY ELEMENTS USING METIS GRAPH PARTITIONING TOOL        |
!    RETURNS[EL_PID(1:MElemGL)]                                                    |
!    EL_PID(I) = Processor ID Of Owner of Element I                            |
!==============================================================================|

   SUBROUTINE DOMDEC(IUNIT,MElemGL,NPROCS,EL_PID,MSR)              
# if defined (MULTIPROCESSOR)

     use mpi
!==============================================================================|
   IMPLICIT NONE
!   include "mpif.h"
   INTEGER, INTENT(IN)  :: IUNIT,MElemGL,NPROCS
   INTEGER, INTENT(OUT) :: EL_PID(MElemGL)
   LOGICAL, INTENT(IN)  :: MSR
   INTEGER, ALLOCATABLE :: NVT(:)
   INTEGER :: I,NTEMP,IERR,ii
!==============================================================================|

!
!----------------READ IN NODE LIST FROM ***_grd.dat FILE-----------------------!
! 

   IF(MSR)THEN
   REWIND(IUNIT)
   ALLOCATE(NVT(3*MElemGL))
   DO I=1,MElemGL*3,3
     READ(IUNIT,*)NTEMP,NVT(I),NVT(I+1),NVT(I+2)
   END DO

!
!-------------DECOMPOSE ELEMENTS USING METIS GRAPH PARTITIONING ---------------!
!


   CALL PARTITION(NPROCS,MElemGL,MAXVAL(NVT),loc(NVT),loc(EL_PID))
   EL_PID = EL_PID + 1
   DEALLOCATE(NVT)
   END IF

!---------------------BROADCAST RESULT TO ALL PROCESSORS-----------------------!

!   do i=1,MElemGL
!   if(msr)write(500,'(2I20)')i,el_pid(i)
!   if(msr)read(500,*)ii,el_pid(i)
!   end do
!   call pstop


   CALL MPI_BCAST(EL_PID,MElemGL,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)


# endif
   END SUBROUTINE DOMDEC
!==============================================================================|
