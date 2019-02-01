!==============================================================================|
!  FLUX CONTROL FOR SALINITY                                                        |
!==============================================================================|

   SUBROUTINE FCT_Q2
!#  if defined (WET_DRY)

!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP):: Q2MAX,Q2MIN
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
     DO K=2,KBM1
       Q2MAX = MAXVAL(Q2(NBSN(I,1:NTSN(I)),K))
       Q2MIN = MINVAL(Q2(NBSN(I,1:NTSN(I)),K))

       IF(K == 2)THEN
         Q2MAX = MAX(Q2MAX,(Q2(I,K)*DZZ(I,K)+Q2(I,K+1)*DZZ(I,K-1))/  &
                 (DZZ(I,K)+DZZ(I,K-1)))
         Q2MIN = MIN(Q2MIN,(Q2(I,K)*DZZ(I,K)+Q2(I,K+1)*DZZ(I,K-1))/  &
                 (DZZ(I,K)+DZZ(I,K-1)))
       ELSE IF(K == KBM1)THEN
         Q2MAX = MAX(Q2MAX,(Q2(I,K)*DZZ(I,K-2)+Q2(I,K-1)*DZZ(I,K-1))/  &
                 (DZZ(I,K-1)+DZZ(I,K-2)))
         Q2MIN = MIN(Q2MIN,(Q2(I,K)*DZZ(I,K-2)+Q2(I,K-1)*DZZ(I,K-1))/  &
                 (DZZ(I,K-1)+DZZ(I,K-2)))
       ELSE
         Q2MAX = MAX(Q2MAX,(Q2(I,K)*DZZ(I,K-2)+Q2(I,K-1)*DZZ(I,K-1))/  &
                 (DZZ(I,K-1)+DZZ(I,K-2)),                              &
                 (Q2(I,K)*DZZ(I,K)+Q2(I,K+1)*DZZ(I,K-1))/              &
                 (DZZ(I,K)+DZZ(I,K-1)))  
         Q2MIN = MIN(Q2MIN,(Q2(I,K)*DZZ(I,K-2)+Q2(I,K-1)*DZZ(I,K-1))/  &
                 (DZZ(I,K-1)+DZZ(I,K-2)),                              &
                 (Q2(I,K)*DZZ(I,K)+Q2(I,K+1)*DZZ(I,K-1))/              &
                 (DZZ(I,K)+DZZ(I,K-1)))
       END IF

       IF(Q2MIN-Q2F(I,K) > 0.0_SP)Q2F(I,K) = Q2MIN
       IF(Q2F(I,K)-Q2MAX > 0.0_SP)Q2F(I,K) = Q2MAX

     END DO
200 CONTINUE
   END DO

100 CONTINUE
   RETURN
!#  endif
   END SUBROUTINE FCT_Q2
!==============================================================================|


