# if !defined (SEMI_IMPLICIT)
!==============================================================================!

   SUBROUTINE ADJUST2D3D(ADJUST_TYPE)
!==============================================================================|
!    ADJUST 3D VELOCITY USING DEFECT BETWEEN UPDATED AND CURRENT VERTICALLY    !
!    AVERAGED VELOCITIES						       !
! 									       !
!    FORMULA IS:							       !
!									       !
!      U_adjusted = U_orig + eps*(U_avg_new - U_avg_current)		       !
!      eps = 0 : no adjustment						       !
!      eps = 1 : full adjustment					       !
!==============================================================================|
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: ADJUST_TYPE
   INTEGER :: I,K
   REAL(SP), PARAMETER :: EPS = 1.0_SP
   REAL(SP) :: UAC,VAC,UTMP,VTMP
!==============================================================================!


   SELECT CASE(ADJUST_TYPE)

   CASE(1)
   DO I=1,MTElem
     UAC    = SUM(U(I,1:KBM1)*DZ1(I,1:KBM1))
     VAC    = SUM(V(I,1:KBM1)*DZ1(I,1:KBM1))
     U(I,1:KBM1) = U(I,1:KBM1) + EPS*(UA(I) - UAC) 
     V(I,1:KBM1) = V(I,1:KBM1) + EPS*(VA(I) - VAC) 
   END DO

   CASE(2)
   UARD = UARD/FLOAT(ISPLIT)
   VARD = VARD/FLOAT(ISPLIT)
!!#  if defined (WET_DRY)
!!   UARDS = UARDS/FLOAT(ISPLIT)
!!   VARDS = VARDS/FLOAT(ISPLIT)
!!#  endif

                                                                                                                         
   DO I=1,MTElem
#    if defined (WET_DRY)
     IF(ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#    endif
       UTMP = 0.0_SP ; VTMP = 0.0_SP
       DO K=1,KBM1
         UTMP = UTMP + U(I,K)*DZ1(I,K)
         VTMP = VTMP + V(I,K)*DZ1(I,K)
       END DO
       UTMP = UTMP*DT1(I)
       VTMP = VTMP*DT1(I)
       DO K=1,KBM1
         U(I,K) = U(I,K) - (UTMP-UARD(I))/DT1(I)
         V(I,K) = V(I,K) - (VTMP-VARD(I))/DT1(I)
       END DO
#    if defined (WET_DRY)
     END IF
#    endif
   END DO

!!#  if defined (WET_DRY)
!!   DO I=1,MTElem
!!     UTMP = 0.0_SP ; VTMP = 0.0_SP
!!     DO K=1,KBM1
!!       UTMP = UTMP + U(I,K)*DZ1(I,K)
!!       VTMP = VTMP + V(I,K)*DZ1(I,K)
!!     END DO
!!     UTMP = UTMP*DT1(I)
!!     VTMP = VTMP*DT1(I)
!!     DO K=1,KBM1
!!       US(I,K) = U(I,K) - (UTMP-UARDS(I))/DT1(I)
!!       VS(I,K) = V(I,K) - (VTMP-VARDS(I))/DT1(I)
!!     END DO
!!   END DO
!!#  endif

   END SELECT

   RETURN
   END SUBROUTINE ADJUST2D3D
!==============================================================================|
# endif
