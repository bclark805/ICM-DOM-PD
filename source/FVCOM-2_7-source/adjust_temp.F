
!==============================================================================|
!     ADJUST TEMPERATURE NEAR RIVER MOUTHS USING ADJACENT NODES                |
!==============================================================================|

   SUBROUTINE ADJUST_TEMP

!==============================================================================|
   USE ALL_VARS
   use mod_par
   USE BCS
   IMPLICIT NONE
   REAL(SP) :: TAVE
   INTEGER :: I,K,JJ,I1,J
!==============================================================================|

   
   IF(NUMQBC > 0)THEN   

   DO K=1,KBM1
     DO I=1,NUMQBC
       JJ=INODEQ(I)
       TAVE = 0.0_SP
       DO J=2,NTSN(JJ)-1
         I1=NBSN(JJ,J)
         TAVE = TAVE + T1(I1,K)
       END DO
       T1(JJ,K) = TAVE/FLOAT(NTSN(JJ)-2)
     END DO
   END DO
   
   END IF

   RETURN
   END SUBROUTINE ADJUST_TEMP
!==============================================================================|

