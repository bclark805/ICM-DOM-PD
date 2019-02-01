
!==============================================================================|
!     ADJUST TEMPERATURE NEAR RIVER MOUTHS USING ADJACENT NODES                |
!     ADJUST SALINITY AT RIVER MOUTHS
!==============================================================================|

   SUBROUTINE ADJUST_TS

!==============================================================================|
   USE ALL_VARS
   use mod_par
   USE BCS
   IMPLICIT NONE
   REAL(SP) :: TAVE,TAVE1,TAVE2
   INTEGER :: I,K,JJ,I1,J,J1,J2,NUM_TAVE,NUM_TAVE1,NUM_TAVE2
!==============================================================================|

!!T.W. comment out averaging adjacent nodes for temperature in this subroutine, minor fluctuations still can be seen, and not
!!sure if it is related to numerical error or code itself. For salinity, keep the original. Nov. 2011.
!! Limited Use of Tave in certain situations requested by T.K 11/8/11


   IF(NUMQBC > 0)THEN   

     IF(INFLOW_TYPE == 'node')THEN
       DO K=1,KBM1
         DO I=1,NUMQBC
           JJ=INODEQ(I)
           TAVE = 0.0_SP
           NUM_TAVE = 0
           DO J=2,NTSN(JJ)-1
             I1=NBSN(JJ,J)
	     IF(NUMQBC == 1)THEN
	       NUM_TAVE = NUM_TAVE + 1
	       TAVE = TAVE + T1(I1,K)
	     ELSE
	       IF(I == 1)THEN
	         IF(I1 /= INODEQ(I+1))THEN 
	         NUM_TAVE = NUM_TAVE + 1
	         TAVE = TAVE + T1(I1,K)
		 END IF
	       ELSE IF(I == NUMQBC)THEN
	         IF(I1 /= INODEQ(I-1))THEN	 
	         NUM_TAVE = NUM_TAVE + 1
	         TAVE = TAVE + T1(I1,K)
                 END IF
	       ELSE IF(I1 /= INODEQ(I-1) .AND. I1 /= INODEQ(I+1))THEN
	         NUM_TAVE = NUM_TAVE + 1
                 TAVE = TAVE + T1(I1,K)
	       END IF	 
	     END IF
           END DO
!           T1(JJ,K) = TAVE/FLOAT(NUM_TAVE)    !!T.W. ! Originally commented out by TW,  TK removed the comment 11/6/6 and tried - did not work
		! Limited Use of Tave in certain situations requested by T.K 11/8/11
		IF(TAVE/FLOAT(NUM_TAVE).LT.TDIS(I)) THEN		!!!T.W. helped add this 11/8/11 test
		T1(JJ,K) = MAX(T1(JJ,K),TAVE/FLOAT(NUM_TAVE))		!!!T.W. helped add this 11/8/11 test
		END IF
         END DO
       END DO
     ELSE IF(INFLOW_TYPE == 'edge')THEN
       DO K=1,KBM1
         DO I=1,NUMQBC
           J1=N_ICELLQ(I,1)
           J2=N_ICELLQ(I,2)
	   TAVE1 = 0.0_SP
	   TAVE2 = 0.0_SP
	   NUM_TAVE1 = 0
	   NUM_TAVE2 = 0

           DO J=2,NTSN(J1)-1
             I1=NBSN(J1,J)
	     IF(NUMQBC == 1)THEN
	       IF(I1 /= J2)THEN
	         NUM_TAVE1 = NUM_TAVE1 + 1
                 TAVE1 = TAVE1 + T1(I1,K)
	       END IF
	     ELSE IF(I == 1)THEN
	       IF(I1 /= J2 .AND. I1 /= N_ICELLQ(I+1,1) .AND. &
	         I1 /= N_ICELLQ(I+1,2))THEN
	         NUM_TAVE1 = NUM_TAVE1 + 1
                 TAVE1 = TAVE1 + T1(I1,K)
	       END IF
	     ELSE IF(I == NUMQBC)THEN
	       IF(I1 /= J2 .AND. I1 /= N_ICELLQ(I-1,1) .AND. &
	         I1 /= N_ICELLQ(I-1,2))THEN
	         NUM_TAVE1 = NUM_TAVE1 + 1
                 TAVE1 = TAVE1 + T1(I1,K)
	       END IF
             ELSE IF(I1 /= J2 .AND. &
	        I1 /= N_ICELLQ(I-1,1) .AND. I1 /= N_ICELLQ(I-1,2) .AND.  &
		I1 /= N_ICELLQ(I+1,1) .AND. I1 /= N_ICELLQ(I+1,2))THEN
	       NUM_TAVE1 = NUM_TAVE1 + 1
               TAVE1 = TAVE1 + T1(I1,K)
	     END IF
           END DO
 !          T1(J1,K) = TAVE1/FLOAT(NUM_TAVE1)    !!T.W. ! Originally commented out by TW,  TK removed the comment 11/6/6 and tried but did not work

           DO J=2,NTSN(J2)-1
             I1=NBSN(J2,J)
	     IF(NUMQBC == 1)THEN
	       IF(I1 /= J1)THEN
	         NUM_TAVE2 = NUM_TAVE2 + 1
                 TAVE2 = TAVE2 + T1(I1,K)
	       END IF
	     ELSE IF(I == 1)THEN
	       IF(I1 /= J1 .AND. I1 /= N_ICELLQ(I+1,1) .AND. &
	         I1 /= N_ICELLQ(I+1,2))THEN
	         NUM_TAVE2 = NUM_TAVE2 + 1
                 TAVE2 = TAVE2 + T1(I1,K)
	       END IF
	     ELSE IF(I == NUMQBC)THEN  
	       IF(I1 /= J1 .AND. I1 /= N_ICELLQ(I-1,1) .AND. &
	         I1 /= N_ICELLQ(I-1,2))THEN
	         NUM_TAVE2 = NUM_TAVE2 + 1
                 TAVE2 = TAVE2 + T1(I1,K)
	       END IF
	     ELSE IF(I1 /= J1 .AND. &
	        I1 /= N_ICELLQ(I-1,1) .AND. I1 /= N_ICELLQ(I-1,2) .AND.  &
		I1 /= N_ICELLQ(I+1,1) .AND. I1 /= N_ICELLQ(I+1,2))THEN
	       NUM_TAVE2 = NUM_TAVE2 + 1
               TAVE2 = TAVE2 + T1(I1,K)
	     END IF
           END DO
  !         T1(J2,K) = TAVE2/FLOAT(NUM_TAVE2)     !!T.W. ! Originally commented out by TW,  TK removed the comment 11/6/06 and tried but did not work
	   
         END DO
       END DO	 
     END IF
   
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

     CALL N2E3D(T1,T)
     CALL N2E3D(S1,S)

   END IF

   RETURN
   END SUBROUTINE ADJUST_TS
!==============================================================================|

