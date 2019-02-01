!==============================================================================|
!   Write Data Averages over Interval INT_AVGE To Assess if Model has          |
!   Achieved Quasi-Periodic Behavior                                           |
!==============================================================================|

   SUBROUTINE OUT_AVGE

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR 
#  endif
#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
#  endif

   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,WWTMP,KMTMP,KHTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: T1TMP,S1TMP,R1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP

#  if defined (BALANCE_2D)
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ADVUA2TMP,  ADVVA2TMP,  ADFX2TMP, ADFY2TMP  
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: DRX2D2TMP,  DRY2D2TMP,  CORX2TMP, CORY2TMP  
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: PSTX2TMP,   PSTY2TMP,   ADX2D2TMP,ADY2D2TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: WUSURBF2TMP,WVSURBF2TMP,DUDT2TMP, DVDT2TMP
#  endif

   CHARACTER(LEN=4)   :: FILENUMBER
   CHARACTER(LEN=120) :: DIR 
   CHARACTER(LEN=120) :: FNAME
   INTEGER  :: J1,J2,END_AVGE
   REAL(SP) :: TMP,FAC

!==============================================================================|
   
   END_AVGE = BEG_AVGE + NUM_AVGE*INT_AVGE - 1
 
!------------------------------------------------------------------------------!
!  ALLOCATE DATA FOR STORING AVERAGES                                          !
!------------------------------------------------------------------------------!

   IF(IINT == BEG_AVGE)THEN
     ALLOCATE(U_AVE(0:MTElem,KB))  ; U_AVE  = 0.0_SP
     ALLOCATE(V_AVE(0:MTElem,KB))  ; V_AVE  = 0.0_SP
     ALLOCATE(W_AVE(0:MTElem,KB))  ; W_AVE  = 0.0_SP
     ALLOCATE(KM_AVE(0:NTNode,KB)) ; KM_AVE = 0.0_SP
     ALLOCATE(KH_AVE(0:NTNode,KB)) ; KH_AVE = 0.0_SP
     ALLOCATE(S_AVE(0:NTNode,KB))  ; S_AVE  = 0.0_SP
     ALLOCATE(T_AVE(0:NTNode,KB))  ; T_AVE  = 0.0_SP
     ALLOCATE(R_AVE(0:NTNode,KB))  ; R_AVE  = 0.0_SP
     ALLOCATE(EL_AVE(0:NTNode))    ; EL_AVE = 0.0_SP

#  if defined (BALANCE_2D)
     ALLOCATE(ADVUA2_AVE(0:MTElem))   ; ADVUA2_AVE   = 0.0_SP
     ALLOCATE(ADVVA2_AVE(0:MTElem))   ; ADVVA2_AVE   = 0.0_SP
     ALLOCATE(ADFX2_AVE(0:MTElem))    ; ADFX2_AVE    = 0.0_SP  
     ALLOCATE(ADFY2_AVE(0:MTElem))    ; ADFY2_AVE    = 0.0_SP  
     ALLOCATE(DRX2D2_AVE(0:MTElem))   ; DRX2D2_AVE   = 0.0_SP
     ALLOCATE(DRY2D2_AVE(0:MTElem))   ; DRY2D2_AVE   = 0.0_SP
     ALLOCATE(CORX2_AVE(0:MTElem))    ; CORX2_AVE    = 0.0_SP  
     ALLOCATE(CORY2_AVE(0:MTElem))    ; CORY2_AVE    = 0.0_SP  
     ALLOCATE(PSTX2_AVE(0:MTElem))    ; PSTX2_AVE    = 0.0_SP
     ALLOCATE(PSTY2_AVE(0:MTElem))    ; PSTY2_AVE    = 0.0_SP
     ALLOCATE(ADX2D2_AVE(0:MTElem))   ; ADX2D2_AVE   = 0.0_SP
     ALLOCATE(ADY2D2_AVE(0:MTElem))   ; ADY2D2_AVE   = 0.0_SP
     ALLOCATE(WUSURBF2_AVE(0:MTElem)) ; WUSURBF2_AVE = 0.0_SP 
     ALLOCATE(WVSURBF2_AVE(0:MTElem)) ; WVSURBF2_AVE = 0.0_SP 
     ALLOCATE(DUDT2_AVE(0:MTElem))    ; DUDT2_AVE    = 0.0_SP    
     ALLOCATE(DVDT2_AVE(0:MTElem))    ; DVDT2_AVE    = 0.0_SP    
#  endif

   END IF

!------------------------------------------------------------------------------!
!  UPDATE AVERAGES                                                             !
!------------------------------------------------------------------------------!

   TMP = INT_AVGE
   FAC = 1.0_SP/TMP
   IF(IINT >= BEG_AVGE .AND. IINT <= END_AVGE)THEN 
     U_AVE  = U_AVE  + U*FAC
     V_AVE  = V_AVE  + V*FAC
     W_AVE  = W_AVE  + WW*FAC
     KM_AVE = KM_AVE + KM*FAC
     KH_AVE = KH_AVE + KH*FAC
     S_AVE  = S_AVE  + S1*FAC
     T_AVE  = T_AVE  + T1*FAC
     R_AVE  = R_AVE  + RHO1*FAC
     EL_AVE = EL_AVE + EL*FAC

#  if defined (BALANCE_2D)
     ADVUA2_AVE   = ADVUA2_AVE   + ADVUA2*FAC
     ADVVA2_AVE   = ADVVA2_AVE   + ADVVA2*FAC
     ADFX2_AVE    = ADFX2_AVE    + ADFX2*FAC
     ADFY2_AVE    = ADFY2_AVE    + ADFY2*FAC
     DRX2D2_AVE   = DRX2D2_AVE   + DRX2D2*FAC
     DRY2D2_AVE   = DRY2D2_AVE   + DRY2D2*FAC
     CORX2_AVE    = CORX2_AVE    + CORX2*FAC
     CORY2_AVE    = CORY2_AVE    + CORY2*FAC
     PSTX2_AVE    = PSTX2_AVE    + PSTX2*FAC
     PSTY2_AVE    = PSTY2_AVE    + PSTY2*FAC
     ADX2D2_AVE   = ADX2D2_AVE   + ADX2D2*FAC
     ADY2D2_AVE   = ADY2D2_AVE   + ADY2D2*FAC
     WUSURBF2_AVE = WUSURBF2_AVE + WUSURBF2*FAC
     WVSURBF2_AVE = WVSURBF2_AVE + WVSURBF2*FAC
     DUDT2_AVE    = DUDT2_AVE    + DUDT2*FAC
     DVDT2_AVE    = DVDT2_AVE    + DVDT2*FAC
#  endif

   END IF

!------------------------------------------------------------------------------!
!  OPEN FILE (Name Based on Iteration Number)                                  !
!------------------------------------------------------------------------------!

   J1 = MOD((IINT+1-BEG_AVGE),INT_AVGE)
   J2 = (IINT+1-BEG_AVGE)/INT_AVGE

   IF(IINT >= BEG_AVGE .AND.  J1 == 0 .AND. IINT <= END_AVGE)THEN

   
   IF(MSR)THEN
     WRITE(FILENUMBER,'(I4.4)') J2 
     FNAME = TRIM(CASENAME)//'_avge'//FILENUMBER//'.dat'
     DIR = TRIM(OUTDIR)//"/out"
     OPEN(1,FILE=TRIM(DIR)//"/"//TRIM(FNAME),STATUS='unknown',FORM='unformatted') 
     REWIND(1)
     WRITE(IPT,*)'DUMPING AVGES FILE: ',TRIM(FNAME)
   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Single Processor Case)                                !
!------------------------------------------------------------------------------!

   IF(SERIAL)THEN

     !! ELEMENT BASED VALUES
     WRITE(1) MElem,NNode,KB,IINT,IINT-INT_AVGE
     DO I=1,MElem
       WRITE(1) (U_AVE(I,K),V_AVE(I,K),W_AVE(I,K),K=1,KBM1)
     END DO
     
     !! NODE BASED VALUES
     DO I=1,NNode
       WRITE(1) (KM_AVE(I,K),KH_AVE(I,K),K=1,KBM1)
     END DO

     DO I=1,NNode
       WRITE(1) EL_AVE(I),(T_AVE(I,K),S_AVE(I,K),R_AVE(I,K),K=1,KBM1)
     END DO

#  if defined (BALANCE_2D)
     
     DO I=1,MElem
       WRITE(1) ADVUA2_AVE(I),   ADVVA2_AVE(I),  ADFX2_AVE(I),  ADFY2_AVE(I),&
                DRX2D2_AVE(I),   DRY2D2_AVE(I),  CORX2_AVE(I),  CORY2_AVE(I),&
		PSTX2_AVE(I),    PSTY2_AVE(I),   ADX2D2_AVE(I), ADY2D2_AVE(I),&
		WUSURBF2_AVE(I), WVSURBF2_AVE(I),DUDT2_AVE(I),  DVDT2_AVE(I)
     ENDDO
#   endif

   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Multi Processor Case)                                 !
!------------------------------------------------------------------------------!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN

     !!GATHER AND WRITE ELEMENT-BASED QUANTITIES (U,V,WW,KH,KM)
     ALLOCATE(UTMP(MElemGL,KB))
     ALLOCATE(VTMP(MElemGL,KB))
     ALLOCATE(WWTMP(MElemGL,KB))
     ALLOCATE(KMTMP(NNodeGL,KB))
     ALLOCATE(KHTMP(NNodeGL,KB))
     CALL GATHER(LBOUND(U,1), UBOUND(U,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP,U_AVE,  UTMP)
     CALL GATHER(LBOUND(V,1), UBOUND(V,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP,V_AVE,  VTMP)
     CALL GATHER(LBOUND(WW,1),UBOUND(WW,1),MElem,MElemGL,KB,MYID,NPROCS,EMAP,W_AVE,  WWTMP)
     CALL GATHER(LBOUND(KM,1),UBOUND(KM,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,KM_AVE, KMTMP)
     CALL GATHER(LBOUND(KH,1),UBOUND(KH,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,KH_AVE, KHTMP)
     IF(MSR)THEN
       WRITE(1) IINT,MElemGL,NNodeGL,THOUR
       DO I=1,MElemGL
         WRITE(1) (UTMP(I,K),VTMP(I,K),WWTMP(I,K),K=1,KBM1)
       END DO
       DO I=1,NNodeGL
         WRITE(1) (KMTMP(I,K),KHTMP(I,K),K=1,KBM1)
       END DO
     END IF
     DEALLOCATE(UTMP,VTMP,WWTMP,KMTMP,KHTMP)

     !!GATHER AND WRITE NODE-BASED QUANTITIES (EL,T1,S1,RHO1)
     ALLOCATE(ELTMP(NNodeGL))
     ALLOCATE(T1TMP(NNodeGL,KB))
     ALLOCATE(S1TMP(NNodeGL,KB))
     ALLOCATE(R1TMP(NNodeGL,KB))
     CALL GATHER(LBOUND(EL,1),  UBOUND(EL,1),  NNode,NNodeGL, 1,MYID,NPROCS,NMAP,EL_AVE,ELTMP)
     CALL GATHER(LBOUND(T1,1),  UBOUND(T1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,T_AVE, T1TMP)
     CALL GATHER(LBOUND(S1,1),  UBOUND(S1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,S_AVE, S1TMP)
     CALL GATHER(LBOUND(RHO1,1),UBOUND(RHO1,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,R_AVE, R1TMP)
     IF(MSR)THEN
       DO I=1,NNodeGL
       WRITE(1) ELTMP(I),(T1TMP(I,K),S1TMP(I,K),R1TMP(I,K),K=1,KBM1)
       END DO
     END IF
     DEALLOCATE(ELTMP,T1TMP,S1TMP,R1TMP)  

#  if defined (BALANCE_2D) 
     !!GATHER AND WRITE ELEMENT-BASED QUANTITIES (ADVUA2,ADVVA2,ADFX2,...)
     ALLOCATE(ADVUA2TMP(MElemGL))
     ALLOCATE(ADVVA2TMP(MElemGL))
     ALLOCATE(ADFX2TMP(MElemGL))
     ALLOCATE(ADFY2TMP(MElemGL))
     ALLOCATE(DRX2D2TMP(MElemGL))
     ALLOCATE(DRY2D2TMP(MElemGL))
     ALLOCATE(CORX2TMP(MElemGL))
     ALLOCATE(CORY2TMP(MElemGL))
     ALLOCATE(PSTX2TMP(MElemGL))
     ALLOCATE(PSTY2TMP(MElemGL))
     ALLOCATE(ADX2D2TMP(MElemGL))
     ALLOCATE(ADY2D2TMP(MElemGL))
     ALLOCATE(WUSURBF2TMP(MElemGL))
     ALLOCATE(WVSURBF2TMP(MElemGL)
     ALLOCATE(DUDT2TMP(MElemGL))
     ALLOCATE(DVDT2TMP(MElemGL))
     CALL GATHER(LBOUND(ADVUA2_AVE,1),  UBOUND(ADVUA2_AVE,1),   MElem,MElemGL,1,MYID,NPROCS,EMAP,ADVUA2_AVE,   ADVUA2TMP)
     CALL GATHER(LBOUND(ADVVA2_AVE,1),  UBOUND(ADVVA2_AVE,1),   MElem,MElemGL,1,MYID,NPROCS,EMAP,ADVVA2_AVE,   ADVVA2TMP)
     CALL GATHER(LBOUND(ADFX2_AVE,1),   UBOUND(ADFX2_AVE,1),    MElem,MElemGL,1,MYID,NPROCS,EMAP,ADFX2_AVE,    ADFX2TMP)
     CALL GATHER(LBOUND(ADFY2_AVE,1),   UBOUND(ADFY2_AVE,1),    MElem,MElemGL,1,MYID,NPROCS,EMAP,ADFY2_AVE,    ADFY2TMP)
     CALL GATHER(LBOUND(DRX2D2_AVE,1),  UBOUND(DRX2D2_AVE,1),   MElem,MElemGL,1,MYID,NPROCS,EMAP,DRX2D2_AVE,   DRX2D2TMP)
     CALL GATHER(LBOUND(DRY2D2_AVE,1),  UBOUND(DRY2D2_AVE,1),   MElem,MElemGL,1,MYID,NPROCS,EMAP,DRY2D2_AVE,   DRY2D2TMP)
     CALL GATHER(LBOUND(CORX2_AVE,1),   UBOUND(CORX2_AVE,1),    MElem,MElemGL,1,MYID,NPROCS,EMAP,CORX2_AVE,    CORX2TMP)
     CALL GATHER(LBOUND(CORY2_AVE,1),   UBOUND(CORY2_AVE,1),    MElem,MElemGL,1,MYID,NPROCS,EMAP,CORY2_AVE,    CORY2TMP)
     CALL GATHER(LBOUND(PSTX2_AVE,1),   UBOUND(PSTX2_AVE,1),    MElem,MElemGL,1,MYID,NPROCS,EMAP,PSTX2_AVE,    PSTX2TMP)
     CALL GATHER(LBOUND(PSTY2_AVE,1),   UBOUND(PSTY2_AVE,1),    MElem,MElemGL,1,MYID,NPROCS,EMAP,PSTY2_AVE,    PSTY2TMP)
     CALL GATHER(LBOUND(ADX2D2_AVE,1),  UBOUND(ADX2D2_AVE,1),   MElem,MElemGL,1,MYID,NPROCS,EMAP,ADX2D2_AVE,   ADX2D2TMP)
     CALL GATHER(LBOUND(ADY2D2_AVE,1),  UBOUND(ADY2D2_AVE,1),   MElem,MElemGL,1,MYID,NPROCS,EMAP,ADY2D2_AVE,   ADY2D2TMP)
     CALL GATHER(LBOUND(WUSURBF2_AVE,1),UBOUND(WUSURBF2_AVE,1), MElem,MElemGL,1,MYID,NPROCS,EMAP,WUSURBF2_AVE, WUSURBF2TMP)
     CALL GATHER(LBOUND(WVSURBF2_AVE,1),UBOUND(WVSURBF2_AVE,1), MElem,MElemGL,1,MYID,NPROCS,EMAP,WVSURBF2_AVE, WVSURBF2TMP)
     CALL GATHER(LBOUND(DUDT2_AVE,1),   UBOUND(DUDT2_AVE,1),    MElem,MElemGL,1,MYID,NPROCS,EMAP,DUDT2_AVE,    DUDT2TMP)
     CALL GATHER(LBOUND(DVDT2_AVE,1),   UBOUND(DVDT2_AVE,1),    MElem,MElemGL,1,MYID,NPROCS,EMAP,DVDT2_AVE,    DVDT2TMP)
     IF(MSR)THEN
       DO I=1,MElemGL
          WRITE(1) ADVUA2TMP(I),   ADVVA2TMP(I),  ADFX2TMP(I),  ADFY2TMP(I),&
                   DRX2D2TMP(I),   DRY2D2TMP(I),  CORX2TMP(I),  CORY2TMP(I),&
                   PSTX2TMP(I),    PSTY2TMP(I),   ADX2D2TMP(I), ADY2D2TMP(I),&
	           WUSURBF2TMP(I), WVSURBF2TMP(I),DUDT2TMP(I),  DVDT2TMP(I)
       ENDDO		   
     ENDIF     
     DEALLOCATE(ADVUA2TMP,  ADVVA2TMP,  ADFX2TMP, ADFY2TMP)  
     DEALLOCATE(DRX2D2TMP,  DRY2D2TMP,  CORX2TMP, CORY2TMP)  
     DEALLOCATE(PSTX2TMP,   PSTY2TMP,   ADX2D2TMP,ADY2D2TMP)  
     DEALLOCATE(WUSURBF2TMP,WVSURBF2TMP,DUDT2TMP, DVDT2TMP)  
#   endif

   END IF
#  endif

   IF(MSR) CLOSE(1)
!------------------------------------------------------------------------------!
!  REINITIALIZE AVERAGING ARRAYS                                               !
!------------------------------------------------------------------------------!

   U_AVE  = 0.0_SP
   V_AVE  = 0.0_SP
   W_AVE  = 0.0_SP
   KM_AVE = 0.0_SP
   KH_AVE = 0.0_SP
   S_AVE  = 0.0_SP
   T_AVE  = 0.0_SP
   R_AVE  = 0.0_SP
   EL_AVE = 0.0_SP
   
#  if defined (BALANCE_2D) 
     ADVUA2_AVE   = 0.0_SP
     ADVVA2_AVE   = 0.0_SP
     ADFX2_AVE    = 0.0_SP  
     ADFY2_AVE    = 0.0_SP  
     DRX2D2_AVE   = 0.0_SP
     DRY2D2_AVE   = 0.0_SP
     CORX2_AVE    = 0.0_SP  
     CORY2_AVE    = 0.0_SP  
     PSTX2_AVE    = 0.0_SP
     PSTY2_AVE    = 0.0_SP
     ADX2D2_AVE   = 0.0_SP
     ADY2D2_AVE   = 0.0_SP
     WUSURBF2_AVE = 0.0_SP 
     WVSURBF2_AVE = 0.0_SP 
     DUDT2_AVE    = 0.0_SP    
     DVDT2_AVE    = 0.0_SP    
#   endif

   END IF

   RETURN
   END SUBROUTINE OUT_AVGE   
!==============================================================================|
