!==============================================================================|
!  AVERAGE THE TEMPERATURE                                                     |
!==============================================================================|

   SUBROUTINE AVER_T
#  if defined (WET_DRY)

!==============================================================================|
   USE ALL_VARS
   USE BCS
   IMPLICIT NONE
   REAL(SP):: AVE_T,TMAX,TMIN
   INTEGER :: I,J,K
!==============================================================================|

   IF(H_TYPE == 'body_h') GO TO 100
   DO I=1,NNode
     IF(NUMQBC > 0)THEN
       DO J=1,NUMQBC
         IF(INFLOW_TYPE == 'node' .AND. I == INODEQ(J))GO TO 200
         IF(INFLOW_TYPE == 'edge' .AND. &
           (I == N_ICELLQ(J,1) .OR. I == N_ICELLQ(J,2)))GO TO 200
       END DO
     END IF
     DO K=1,KBM1
       TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
       TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))

       IF(K == 1)THEN
         TMAX = MAX(TMAX,(T1(I,K)*DZ(K+1)+T1(I,K+1)*DZ(K))/(DZ(K)+DZ(K+1)))
         TMIN = MIN(TMIN,(T1(I,K)*DZ(K+1)+T1(I,K+1)*DZ(K))/(DZ(K)+DZ(K+1)))
       ELSE IF(K == KBM1)THEN
         TMAX = MAX(TMAX,(T1(I,K)*DZ(K-1)+T1(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1)))
         TMIN = MIN(TMIN,(T1(I,K)*DZ(K-1)+T1(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1)))
       ELSE
         TMAX = MAX(TMAX,(T1(I,K)*DZ(K-1)+T1(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1)), &
                 (T1(I,K)*DZ(K+1)+T1(I,K+1)*DZ(K))/(DZ(K)+DZ(K+1)))
         TMIN = MIN(TMIN,(T1(I,K)*DZ(K-1)+T1(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1)), &
                 (T1(I,K)*DZ(K+1)+T1(I,K+1)*DZ(K))/(DZ(K)+DZ(K+1)))
       END IF

       IF(TMIN-TF1(I,K) > 0.0_SP)TF1(I,K) = TMIN
       IF(TF1(I,K)-TMAX > 0.0_SP)TF1(I,K) = TMAX

     END DO
200 CONTINUE
   END DO

100 CONTINUE
   RETURN
#  endif
   END SUBROUTINE AVER_T
!==============================================================================|


