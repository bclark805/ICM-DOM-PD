
!==============================================================================|
!     ADJUST TEMPERATURE NEAR RIVER MOUTHS USING ADJACENT NODES                |
!     ADJUST SALINITY AT RIVER MOUTHS
!==============================================================================|

   SUBROUTINE ADJUST_T

!==============================================================================|
   USE ALL_VARS
   use mod_par
   USE BCS
   IMPLICIT NONE
   REAL(SP) :: TAVE,TAVE1,TAVE2
   INTEGER :: I,K,JJ,I1,J,J1,J2,NUM_TAVE,NUM_TAVE1,NUM_TAVE2
!==============================================================================|

   
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
           T1(JJ,K) = TAVE/FLOAT(NUM_TAVE)
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
           T1(J1,K) = TAVE1/FLOAT(NUM_TAVE1)

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
           T1(J2,K) = TAVE2/FLOAT(NUM_TAVE2)
	   
         END DO
       END DO	 
     END IF
  
     CALL N2E3D(T1,T)

   END IF

   RETURN
   END SUBROUTINE ADJUST_T
!==============================================================================|

