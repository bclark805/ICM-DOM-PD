!==============================================================================!
!
!==============================================================================!

MODULE MOD_ATMOTIDE
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif   
   USE CONTROL
   IMPLICIT NONE
   SAVE


   REAL(SP), ALLOCATABLE :: ELF_ATMO(:),ELRK_ATMO(:),EL_ATMO(:),EGF_ATMO(:)
   REAL(SP), ALLOCATABLE :: PHI(:),LAMDA(:)
   
!--The amplitudes of atmospheric tides----------------------------------------------
   REAL(SP),PARAMETER :: APT_ATMO = 0.0113_SP            !S2
!--The frequncy of atmospheric tides------------------------------------------------
!   REAL(SP),PARAMETER :: FREQ_ATMO = 1.454410E-4_SP    !S2  
   REAL(SP) :: FREQ_ATMO                                !S2  
!-------------------------------------------------------------------------------
   REAL(SP) :: ALFA_ATMO = 112.0_SP                !S2   

   INTEGER  :: APT_FACT_ATMO
   CONTAINS


!==========================================================================|
!==========================================================================|
   SUBROUTINE ALLOCATE_ATMO
   USE ALL_VARS
   IMPLICIT NONE
   
   INTEGER :: I
   
   ALLOCATE(ELF_ATMO(0:NTNode)); ELF_ATMO = ZERO
   ALLOCATE(ELRK_ATMO(0:NTNode)); ELRK_ATMO = ZERO
   ALLOCATE(EL_ATMO(0:NTNode)); EL_ATMO = ZERO
   ALLOCATE(EGF_ATMO(0:NTNode)); EGF_ATMO = ZERO
   ALLOCATE(PHI(0:NTNode))   ; PHI     = ZERO
   ALLOCATE(LAMDA(0:NTNode))   ; LAMDA   = ZERO
   
#  if defined (SPHERICAL)
   IF(SERIAL)THEN
     PHI   = YG
     LAMDA = XG
   END IF  
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,NNode
       PHI(I)   = YG(NGID(I))
       LAMDA(I) = XG(NGID(I))
     END DO
     DO I=1,NHN
       PHI(I+NNode)   = YG(HN_LST(I))
       LAMDA(I+NNode) = XG(HN_LST(I))
     END DO
   END IF
#  endif
#  else
   IF(MSR) PRINT*,"THE EQUILIBRIUM TIDE HAS NOT BEEN ADDED IN THE ",     &
                  "NON SPHERICAL COORDINATE"
   CALL PSTOP		    
#  endif     

   FREQ_ATMO = PI2/PERIOD(1)                                !S2  
   ALFA_ATMO  = ALFA_ATMO*PI2/360.0_SP

   RETURN
   END SUBROUTINE ALLOCATE_ATMO
!==========================================================================|


!==========================================================================|
   SUBROUTINE ELEVATION_ATMO

!--------------------------------------------------------------------------|
!  Surface Elevation of ATMOSPHERIC TIDE                                   |
!--------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_OBCS
   IMPLICIT NONE

   INTEGER :: I,J
   REAL(SP):: TIME1
   REAL(SP):: FORCE,PHAI_IJ

   TIME1 = TIME * 86400.0_SP
!
!-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times-|
!
   IF(S_TYPE == 'julian')THEN
!  not finish yet
   END IF

!
!-Non-Julian: Set Elevation of Atmospheric Tide -----------------------------|
!

   IF(S_TYPE == 'non-julian')THEN
     DO I = 1, NTNode
       PHAI_IJ = LAMDA(I)*PI2/360.0_SP
       FORCE = APT_ATMO*APT_FACT_ATMO*COS(FREQ_ATMO*TIME1+2.0_SP*PHAI_IJ-ALFA_ATMO)
       ELF_ATMO(I) = FORCE * RAMP
     END DO

   END IF

   RETURN
   END SUBROUTINE ELEVATION_ATMO
!============================================================================|
!============================================================================|

END MODULE MOD_ATMOTIDE
