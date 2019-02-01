!==============================================================================!
!
!==============================================================================!

MODULE MOD_EQUITIDE
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif   
   USE CONTROL
   IMPLICIT NONE
   SAVE

   REAL(SP), ALLOCATABLE :: ELF_EQI(:),ELRK_EQI(:),EL_EQI(:),EGF_EQI(:)
   REAL(SP), ALLOCATABLE :: PHI(:),LAMDA(:)
   
!--The amplitudes of equilibrium tides----------------------------------------------
   REAL(SP),PARAMETER, DIMENSION(8) :: &

!                     S2           M2           N2              K2
     APT_EQI = (/0.112743_SP, 0.242334_SP, 0.046397_SP,  0.030684_SP,  &
!                     K1           P1           O1              Q1
                 0.141565_SP, 0.046848_SP, 0.100661_SP,  0.019273_SP/)

!--The frequncy of equilibrium tides------------------------------------------------
   REAL(SP),PARAMETER, DIMENSION(8) :: &
!                     S2              M2              N2              K2            
    FREQ_EQI = (/1.454441E-4_SP, 1.405189E-4_SP, 1.378797E-4_SP, 1.458423E-4_SP, &
!                     K1              P1              O1              Q1
                0.7292117E-4_SP,0.7252295E-4_SP,0.6759774E-4_SP,0.6495854E-4_SP/)  
!--The Love number k-----------------------------------------------------------------
!   REAL(SP),PARAMETER, DIMENSION(8) :: &
!                        S2        M2        N2        K1        P1        O1
!           K_LOVE = (/0.302_SP, 0.302_SP, 0.302_SP, 0.256_SP, 0.000_SP, 0.298_SP/)   
!--The Love number h -----------------------------------------------------------------
!   REAL(SP),PARAMETER, DIMENSION(8) :: &
!                        S2        M2        N2        K1        P1        O1
!           H_LOVE = (/0.602_SP, 0.602_SP, 0.602_SP, 0.520_SP, 0.000_SP, 0.603_SP/)   	   

!--The Love number BATA---------------------------------------------------------------
   REAL(SP),PARAMETER, DIMENSION(8) :: &
!                    S2        M2        N2      K2    
      BATA = (/0.693_SP, 0.693_SP, 0.693_SP, 0.693_SP, & 
!                    K1        P1        O1      Q1
               0.736_SP, 0.706_SP, 0.695_SP, 0.695_SP/)

   INTEGER   :: APT_FACT_EQUI(8)
   
   CONTAINS


!==========================================================================|
!==========================================================================|
   SUBROUTINE ALLOCATE_EQUI
   USE ALL_VARS
   IMPLICIT NONE
   
   INTEGER :: I
   
   ALLOCATE(ELF_EQI(0:NTNode)); ELF_EQI = ZERO
   ALLOCATE(ELRK_EQI(0:NTNode)); ELRK_EQI = ZERO
   ALLOCATE(EL_EQI(0:NTNode)); EL_EQI = ZERO
   ALLOCATE(EGF_EQI(0:NTNode)); EGF_EQI = ZERO
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

   RETURN
   END SUBROUTINE ALLOCATE_EQUI
!==========================================================================|


!==========================================================================|
   SUBROUTINE ELEVATION_EQUI

!--------------------------------------------------------------------------|
!  Surface Elevation of EQUILIBRIUM TIDE                                   |
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
!-Non-Julian: Set Elevation of Equilibrium Tide -----------------------------|
!

   IF(S_TYPE == 'non-julian')THEN
     DO I = 1, NTNode
       FORCE = 0.0_SP
       DO J = 1,8
         IF(J <= 4)THEN
           PHAI_IJ = LAMDA(I)*PI2/360.0_SP
           FORCE = BATA(J)*APT_EQI(J)*APT_FACT_EQUI(J)*COS(PHI(I)*PI2/360.0_SP)**2         &
                   *COS(FREQ_EQI(J)*TIME1+2.0_SP*PHAI_IJ) + FORCE
         ELSE
           PHAI_IJ = LAMDA(I)*PI2/360.0_SP
           FORCE = BATA(J)*APT_EQI(J)*APT_FACT_EQUI(J)*SIN(2.0_SP*PHI(I)*PI2/360.0_SP)     &
                  *COS(FREQ_EQI(J)*TIME1+PHAI_IJ) + FORCE
       END IF  
       END DO
       FORCE = FORCE
       ELF_EQI(I) = FORCE * RAMP
     END DO
   END IF

   RETURN
   END SUBROUTINE ELEVATION_EQUI
!============================================================================|
!============================================================================|

END MODULE MOD_EQUITIDE
