!==============================================================================|
!  FLUX CONTROL FOR TEMPERATURE                                                     |
!==============================================================================|

   SUBROUTINE FCT_T
!#  if defined (WET_DRY)

!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP):: TMAX,TMIN
   INTEGER :: I,J,K
!==============================================================================|

   IF(H_TYPE == 'body_h') GO TO 100
   DO I=1,NNode
     IF(IOBCN > 0)THEN
       DO J=1,IOBCN
         IF(I == I_OBC_N(J))GO TO 200
       END DO
     END IF  	 
     IF(NUMQBC > 0)THEN
       DO J=1,NUMQBC
         IF(INFLOW_TYPE == 'node')THEN
	   IF(I == INODEQ(J))GO TO 200
	 END IF  
         IF(INFLOW_TYPE == 'edge')THEN
	   IF(I == N_ICELLQ(J,1) .OR. I == N_ICELLQ(J,2))GO TO 200
	 END IF  
       END DO
     END IF
     DO K=1,KBM1
       TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
       TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))

       IF(K == 1)THEN
         TMAX = MAX(TMAX,(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
         TMIN = MIN(TMIN,(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
       ELSE IF(K == KBM1)THEN
         TMAX = MAX(TMAX,(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
         TMIN = MIN(TMIN,(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
       ELSE
         TMAX = MAX(TMAX,(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
         TMIN = MIN(TMIN,(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
       END IF

       IF(TMIN-TF1(I,K) > 0.0_SP)TF1(I,K) = TMIN
       IF(TF1(I,K)-TMAX > 0.0_SP)TF1(I,K) = TMAX

     END DO
200 CONTINUE
   END DO

100 CONTINUE
   RETURN
!#  endif
   END SUBROUTINE FCT_T
!==============================================================================|


