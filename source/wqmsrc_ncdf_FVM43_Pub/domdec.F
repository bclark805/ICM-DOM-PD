!subroutine DOMDEC()

!==============================================================================|
!  DECOMPOSE THE DOMAIN BY ELEMENTS USING METIS GRAPH PARTITIONING TOOL        |
!    RETURNS[EL_PID(1:NGL)]                                                    |
!    EL_PID(I) = Processor ID Of Owner of Element I                            |
!==============================================================================|

   SUBROUTINE DOMDEC(NGL,NPROCS,EL_PID,MSR)              
# if defined (MULTIPROCESSOR)

   USE MOD_HYDROVARS, ONLY : NVG 
   
!==============================================================================|
   IMPLICIT NONE
   include "mpif.h"
   INTEGER, INTENT(IN)  :: NGL,NPROCS
   INTEGER, INTENT(OUT) :: EL_PID(NGL)
   LOGICAL, INTENT(IN)  :: MSR
   INTEGER, ALLOCATABLE :: NVT(:)
   INTEGER :: I,NTEMP,IERR,ii
!==============================================================================|

!
!----------------READ IN NODE LIST FROM ***_grd.dat FILE-----------------------!
!

   IF(MSR)THEN
   ALLOCATE(NVT(3*NGL))

! KURT GLAESEMANN - fix NVT creation
!  II = 0
!  DO I=1,NGL
!    II = II+1
!    NVT(II)   = NVG(I,1)    !TEMP(I,1)
!    NVT(II+1) = NVG(I,1)    !TEMP(I,3)
!    NVT(II+2) = NVG(I,1)    !TEMP(I,2)
!  ENDDO
   II = -2
   DO I=1,NGL
     II = II+3
     NVT(II)   = NVG(I,1)    !TEMP(I,1)
     NVT(II+1) = NVG(I,2)    !TEMP(I,3)
     NVT(II+2) = NVG(I,3)    !TEMP(I,2)
   ENDDO

!
!-------------DECOMPOSE ELEMENTS USING METIS GRAPH PARTITIONING ---------------!
!
   CALL PARTITION(NPROCS,NGL,maxval(NVT),loc(NVT),loc(EL_PID))

   EL_PID = EL_PID + 1
   DEALLOCATE(NVT)

   ENDIF

   CALL MPI_BCAST(EL_PID,NGL,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)


# endif
   END SUBROUTINE DOMDEC
!==============================================================================|
