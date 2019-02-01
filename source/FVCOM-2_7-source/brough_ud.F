!==============================================================================|
!  User Defined Bottom Drag Coefficient                                        !
!   Calculate Bottom Friction based on Regional Knowledge of Bottom Character  !
!   Current Configuration: GOM                                             ! 
!                                                                              ! 
!   CALCULATES:   WUBOT(MElem), WVBOT(MElem) : BOTTOM SHEAR STRESSES                   !
!==============================================================================|

   SUBROUTINE BROUGH_UD

!==============================================================================!
   USE ALL_VARS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif      
   IMPLICIT NONE
   INTEGER :: I
   REAL(SP), PARAMETER  :: VK2 = .160_SP  !!VON KARMAN LENGTH SCALE SQUARED
   REAL(SP)             :: CBCMIN,Z0,ZTEMP,BTPS
!==============================================================================!

!
!  SET CONSTANTS
!
   CBCMIN = BFRIC
   Z0     = Z0B

!
!  CALCULATE BOTTOM DRAG COEFFICIENT 
!

     DO I=1, MElem
       IF(D1(I) <= 40.0_SP) THEN
         Z0=3.E-3_SP
       ELSE IF(D1(I) > 40.0_SP.AND.D1(I) <= 70.0_SP) THEN
         Z0=3.E-3_SP*EXP(-(D1(I)-40.0_SP)/8.8204_SP)
       ELSE IF(D1(I) > 70.0_SP.AND.D1(I) <= 100.0_SP) THEN
         Z0=1.E-4_SP*EXP(-(D1(I)-70.0_SP)/13.0288_SP)
       ELSE
         Z0=1.E-5_SP
       END IF

       ZTEMP=(ZZ1(I,KBM1)-Z1(I,KB))*DT1(I)/Z0
       CBC(I) = VK2/(LOG(ZTEMP))**2
       IF(ZTEMP < 286.2468_SP) CBC(I)=0.005_SP
       END DO


!==============================================================================|
!  CALCULATE SHEAR STRESS ON BOTTOM  --> WUBOT/WVBOT                           |
!==============================================================================|
   DO I=1,MElem
     IF(D1(I) > 0.0_SP) THEN
       BTPS = CBC(I)*SQRT(U(I,KBM1)**2 + V(I,KBM1)**2)
       WUBOT(I) = -BTPS * U(I,KBM1)
       WVBOT(I) = -BTPS * V(I,KBM1)
     ELSE
       WUBOT(I) = 0.0_SP
       WVBOT(I) = 0.0_SP
     END IF
   END DO

   RETURN
   END SUBROUTINE BROUGH_UD
!==============================================================================|

