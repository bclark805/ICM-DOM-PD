!==============================================================================|
!  FLUX CONTROL FOR SALINITY                                                        |
!==============================================================================|

   SUBROUTINE FCT_S
!#  if defined (WET_DRY)

!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP):: SMAX,SMIN
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
       SMAX = MAXVAL(S1(NBSN(I,1:NTSN(I)),K))
       SMIN = MINVAL(S1(NBSN(I,1:NTSN(I)),K))

       IF(K == 1)THEN
         SMAX = MAX(SMAX,(S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
         SMIN = MIN(SMIN,(S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
       ELSE IF(K == KBM1)THEN
         SMAX = MAX(SMAX,(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
         SMIN = MIN(SMIN,(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
       ELSE
         SMAX = MAX(SMAX,(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
         SMIN = MIN(SMIN,(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
       END IF

       IF(SMIN-SF1(I,K) > 0.0_SP)SF1(I,K) = SMIN
       IF(SF1(I,K)-SMAX > 0.0_SP)SF1(I,K) = SMAX

     END DO
200 CONTINUE
   END DO

100 CONTINUE
   RETURN
!#  endif
   END SUBROUTINE FCT_S
!==============================================================================|


