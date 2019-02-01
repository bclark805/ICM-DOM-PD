
!==============================================================================|
!     ADJUST TEMPERATURE NEAR RIVER MOUTHS USING ADJACENT NODES                |
!     ADJUST SALINITY AT RIVER MOUTHS
!==============================================================================|

   SUBROUTINE ADJUST_S

!==============================================================================|
   USE ALL_VARS
   use mod_par
   USE BCS
   IMPLICIT NONE
   REAL(SP) :: TAVE,TAVE1,TAVE2
   INTEGER :: I,K,JJ,I1,J,J1,J2,NUM_TAVE,NUM_TAVE1,NUM_TAVE2
!==============================================================================|

   
   IF(NUMQBC > 0)THEN   

     DO I=1,NUMQBC
       IF(INFLOW_TYPE == 'node')THEN
         J = INODEQ(I)
         DO K=1,KBM1
           S1(J,K) = MAX(S1(J,K),SDIS(I))
         END DO
       ELSE IF(INFLOW_TYPE == 'edge')THEN
         J1 = N_ICELLQ(I,1)
         J2 = N_ICELLQ(I,2)
         DO K=1,KBM1
           S1(J1,K) = MAX(S1(J1,K),SDIS(I))
           S1(J2,K) = MAX(S1(J2,K),SDIS(I))
         END DO
       END IF
     END DO

     CALL N2E3D(S1,S)

   END IF

   RETURN
   END SUBROUTINE ADJUST_S
!==============================================================================|

