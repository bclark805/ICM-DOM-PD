
!==============================================================================|
!     READ MESH AND DETERMINE NUMBER OF ELEMENTS AND NODES                     |
!==============================================================================|

   SUBROUTINE GETDIM

!==============================================================================|
   USE ALL_VARS
   USE MOD_OBCS
   IMPLICIT NONE
   INTEGER, ALLOCATABLE :: NVTMP(:,:)
   INTEGER              :: I,LM1,J,IOS,LMAX
!==============================================================================|

   IF(MSR)WRITE(IPT,*)'!                  GLOBAL INFORMATION                            !'
   IF(MSR)WRITE(IPT,*)'!                                                                !'
   IF(MSR)WRITE(IPT,*)'!  # OF PROCESSORS       :',NPROCS
!
!----------------Determine Number of Elements----------------------------------!
!
   LM1 = 1
   DO WHILE(.TRUE.)
     READ(INGRD,*,IOSTAT=IOS)J
     IF(IOS < 0)THEN
       IF(MSR)WRITE(IPT,*)'ERROR DETERMINING GRID DIMENSIONS'
       CALL PSTOP
     END IF
     IF(J == 1 .AND. LM1 /= 1)THEN
       MElemGL = LM1
       EXIT
     END IF
     LM1 = J
   END DO
   REWIND(INGRD)

!
!----------------Determine Number of Nodes-------------------------------------!
!
   ALLOCATE(NVTMP(MElemGL,3))
   DO I=1,MElemGL
     READ(INGRD,*)J,NVTMP(I,1),NVTMP(I,2),NVTMP(I,3)
   END DO
   NNodeGL = MAX(MAXVAL(NVTMP(:,1)) , MAXVAL(NVTMP(:,2)) , MAXVAL(NVTMP(:,3)))
   DEALLOCATE(NVTMP)
   REWIND(INGRD)
!
!----------------Determine Number of Nodes on Outer Boundary-------------------!
!
   READ(INOBC,*)
   LMAX = 0
   DO WHILE(.TRUE.)
     READ(INOBC,*,IOSTAT=IOS)J
     IF(IOS < 0)EXIT
     LMAX = J
   END DO
   REWIND(INOBC)
   IOBCN_GL = LMAX

!T.W. added for open boundary node#, per RGL
   IOBCN_RGL = IOBCN_GL
!Finish T.W.
!
!----------------Determine Number of Groundwater Flux Points ------------------!
!
   READ(INBFW,*)
   READ(INBFW,*) IBFW_GL
   REWIND(INBFW)
!
!----------------Determine Number of Freshwater Flux Points--------------------!
!
   READ(INRIV,*)
   READ(INRIV,*) NUMQBC_GL
   REWIND(INRIV)
   
!
!------Determine Number of Frictional Geostrophic Correction Inflow Nodes------!
!
   NOBCGEO_GL = 0
   IF(JMPOBC) READ(INJMP,*)NOBCGEO_GL
   REWIND(INJMP)
!
!----------------Report--------------------------------------------------------!
!
   IF(MSR)WRITE(IPT,*)'!  # OF NODES            :',NNodeGL
   IF(MSR)WRITE(IPT,*)'!  # OF ELEMENTS         :',MElemGL
   IF(MSR)WRITE(IPT,*)'!  # OF OPEN BNDRY NODES :',IOBCN_GL
   IF(MSR)WRITE(IPT,*)'!  # OF BOTTOM FLUX PTS  :',IBFW_GL
   IF(MSR)WRITE(IPT,*)'!  # OF FRESH WATER PTS  :',NUMQBC_GL
   IF(MSR)WRITE(IPT,*)'!  # OF JMP OBC NODES    :',NOBCGEO_GL

!
!-------------Set Local Number of Elements/Nodes for Serial Case---------------!
!
   MElem = MElemGL
   NNode = NNodeGL
   MTElem = MElem
   NTNode = NNode

!
!-------------Read In Global Node Numbering------------------------------------!
!
   ALLOCATE(NVG(0:MElemGL,4))
   DO I=1,MElemGL
     READ(INGRD,*) J,NVG(I,1),NVG(I,3),NVG(I,2)
     NVG(I,4) = NVG(I,1)
   END DO
   IF(MSR)WRITE(IPT,*)  '!  GLOBAL MESH READING   :    FINISHED'


   RETURN
   END SUBROUTINE GETDIM
!==============================================================================|
