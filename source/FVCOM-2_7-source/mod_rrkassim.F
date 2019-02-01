MODULE MOD_RRKA 
# if defined (RRK_ASSIM) || defined (RRK_PRE)
   USE CONTROL
   IMPLICIT NONE
   SAVE
   
   INTEGER NLOC
   INTEGER IDUM
   INTEGER NCYC
   INTEGER N1CYC
   INTEGER ICYC
   INTEGER I_INITIAL
   INTEGER INORRK
   INTEGER OUTERR
   INTEGER,ALLOCATABLE :: STLOC(:)

   REAL(DP),ALLOCATABLE :: STFCT(:)
   REAL(DP),ALLOCATABLE :: STAN1(:)
   REAL(DP),ALLOCATABLE :: CRSTATE(:)
   REAL(DP),ALLOCATABLE :: KAL(:,:)
   REAL(DP),ALLOCATABLE :: MODDATA(:)
   REAL(DP),ALLOCATABLE :: OBSDATA(:)
   REAL(DP),ALLOCATABLE :: ERRVEC(:)
   REAL(DP),ALLOCATABLE :: R(:,:)
   REAL(DP),ALLOCATABLE :: OBSERR(:)
   REAL(DP),ALLOCATABLE :: STTRUE(:,:)
   REAL(DP),ALLOCATABLE :: STTRUE1(:)
   REAL(DP),ALLOCATABLE :: STTR1(:)
   REAL(DP),ALLOCATABLE :: STTR0(:)
   REAL(DP),ALLOCATABLE :: SEOF(:,:)
   REAL(DP),ALLOCATABLE :: STTEMP1(:)
   REAL(DP),ALLOCATABLE :: STSD(:)
   REAL(DP),ALLOCATABLE :: RRKEL2(:)
   REAL(DP),ALLOCATABLE :: RRKU2(:,:)
   REAL(DP),ALLOCATABLE :: RRKV2(:,:)
   REAL(DP),ALLOCATABLE :: RRKT2(:,:)
   REAL(DP),ALLOCATABLE :: RRKS2(:,:)
   REAL(DP)             :: ELSD, USD, TSD, SSD
   REAL(DP)             :: RRKSUM
   REAL(DP)             :: ERR2_INN_FCT
   REAL(DP)             :: ERR2_INN_AN1
   REAL(DP)             :: ERR2_TOT_FCT
   REAL(DP)             :: ERR2_TOT_AN1   
   REAL                 :: RNOBS
   
   CHARACTER(LEN=6)     :: FCYC
   CHARACTER(LEN=20)    :: TEXT
   CHARACTER(LEN=80)    :: ERRFILE
   CHARACTER(LEN=80)    :: AN1FILE
   CHARACTER(LEN=80)    :: FILENAME   

   INTEGER  TIMEN
   REAL(DP),ALLOCATABLE  :: EL_SRS(:,:),SRS_TMP(:)
   REAL(DP),ALLOCATABLE  :: TIME_SER(:)
   REAL(DP) AMP(6),PHAS(6),PHAI_IJ,FORCE

   CONTAINS

   SUBROUTINE RRK_RRKF
!-------------------------------------------------------------------------------|
! START OF UPDATING THE FORCAST BY THE RRKF                                     | 
!-------------------------------------------------------------------------------|

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE RRKVAL
   USE MOD_RRK
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER IDUMMY
   INTEGER I,J,K
   REAL(DP) :: ERR1
   REAL(DP) :: ERR2
   REAL(DP) :: ERR3
   INTEGER SS_DIM
   INTEGER STDIM
   INTEGER IERR
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: UATMP,VATMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TTMP,STMP
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1

   SS_DIM = RRK_NVAR*RRK_NEOF
   
!rrkf   Set the output file name for the analysis
   WRITE(FCYC,'(I6.6)') ICYC*DELTA_ASS
   AN1FILE=TRIM(OUTDIR)//'/rrktemp/'//'analysis'//FCYC//'.cdf'

   IF (SERIAL) THEN
     IF(EL_ASSIM) THEN
       DO I=1, NNodeGL
         RRKEL2(I) = EL(I)
       ENDDO
     ENDIF
     IF(UV_ASSIM) THEN  
       DO K=1, KB
         DO I=1, MElemGL
           RRKU2(I,K)  = U(I,K)
           RRKV2(I,K)  = V(I,K)
         ENDDO
       ENDDO
!       DO I=1, MElemGL
!         RRKU2(I,1)  = UA(I)
!         RRKV2(I,1)  = VA(I)
!       ENDDO 
     ENDIF
     IF(T_ASSIM) THEN
       DO K=1, KB
         DO I=1, NNodeGL
           RRKT2(I,K)  = T1(I,K)
         ENDDO
       ENDDO
     ENDIF
     IF(S_ASSIM) THEN
       DO K=1, KB
         DO I=1, NNodeGL
           RRKS2(I,K)  = S1(I,K)
         ENDDO
       ENDDO
     ENDIF
   ELSE
#    if defined(MULTIPROCESSOR)  
       IF(PAR)THEN   
         ALLOCATE(ELTMP(NNodeGL))        ; ELTMP = 0.0_SP
         ALLOCATE(UATMP(MElemGL))        ; UATMP = 0.0_SP
         ALLOCATE(VATMP(MElemGL))        ; VATMP = 0.0_SP
         ALLOCATE(UTMP(MElemGL,KB))      ; UTMP  = 0.0_SP
         ALLOCATE(VTMP(MElemGL,KB))      ; VTMP  = 0.0_SP
	 ALLOCATE(TTMP(NNodeGL,KB))      ; TTMP  = 0.0_SP
	 ALLOCATE(STMP(NNodeGL,KB))      ; STMP  = 0.0_SP

         CALL GATHER(LBOUND(EL,1), UBOUND(EL,1), NNode,NNodeGL,1 ,MYID,NPROCS,NMAP,EL, ELTMP)
         CALL GATHER( LBOUND(U,1),  UBOUND(U,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP, U,  UTMP) 
         CALL GATHER( LBOUND(V,1),  UBOUND(V,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP, V,  VTMP) 
         CALL GATHER(LBOUND(UA,1), UBOUND(UA,1), MElem,MElemGL,1 ,MYID,NPROCS,EMAP,UA, UATMP) 
         CALL GATHER(LBOUND(VA,1), UBOUND(VA,1), MElem,MElemGL,1 ,MYID,NPROCS,EMAP,VA, VATMP)
	 CALL GATHER(LBOUND(T1,1), UBOUND(T1,1), NNode,NNodeGL,KB,MYID,NPROCS,EMAP,T1,  TTMP)
	 CALL GATHER(LBOUND(S1,1), UBOUND(S1,1), NNode,NNodeGL,KB,MYID,NPROCS,EMAP,S1,  STMP)
	 
	 IF(EL_ASSIM) THEN
	   DO I=1, NNodeGL
             RRKEL2(I) = ELTMP(I)
           ENDDO
           CALL MPI_BCAST(RRKEL2,NNodeGL,MPI_F,0,MPI_COMM_WORLD,IERR)
	 ENDIF
         IF(UV_ASSIM) THEN
	   DO K=1, KB
             DO I=1, MElemGL
               RRKU2(I,K)  = UTMP(I,K)
               RRKV2(I,K)  = VTMP(I,K)
             ENDDO
           ENDDO
!           DO I=1, MElemGL
!	     RRKU2(I,1) = UATMP(I)
!	     RRKV2(I,1) = VATMP(I)
!	   ENDDO
           CALL MPI_BCAST(RRKU2,MElemGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
           CALL MPI_BCAST(RRKV2,MElemGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
	 ENDIF
	 IF(T_ASSIM) THEN
	   DO K=1, KB
             DO I=1, NNodeGL
               RRKT2(I,K)  = TTMP(I,K)
             ENDDO
           ENDDO
           CALL MPI_BCAST(RRKT2,NNodeGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
	 ENDIF
	 IF(S_ASSIM) THEN
	   DO K=1, KB
             DO I=1, NNodeGL
               RRKS2(I,K)  = STMP(I,K)
             ENDDO
           ENDDO
	   CALL MPI_BCAST(RRKS2,NNodeGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
	 ENDIF	 
       
         DEALLOCATE(ELTMP,UTMP,VTMP,TTMP,STMP,UATMP,VATMP)
       END IF
#    endif
   ENDIF

!rrkf   Save the forecast state as 'stfct' 
    IDUMMY = 0
    IF(EL_ASSIM) THEN
      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        STFCT(IDUMMY)= RRKEL2(I)
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MElemGL 
           IDUMMY = IDUMMY + 1         
	   STFCT(IDUMMY) = RRKU2(J,I)
        ENDDO  
      ENDDO
      DO I=1, KBM1
        DO J=1, MElemGL
           IDUMMY = IDUMMY + 1
           STFCT(IDUMMY) = RRKV2(J,I)
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL 
           IDUMMY = IDUMMY + 1         
	   STFCT(IDUMMY) = RRKT2(J,I)
        ENDDO  
      ENDDO
    ENDIF
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL 
           IDUMMY = IDUMMY + 1         
	   STFCT(IDUMMY) = RRKS2(J,I)
        ENDDO  
      ENDDO
    ENDIF

!rrkf   Get the observations and model counterparts
! y=H(x^t) -> obsdata
!    CALL GETOBSDATA
!     CALL GETOBSLOC1
     CALL READMEDM2
     STTRUE1 = STTR0

     DO I=1, NLOC
      OBSDATA(I) = STTRUE1(STLOC(I))
     ENDDO

!rrkf   Perturb the observation with random errors of prescribed observatin error 'obserr'
!rrkf   Generate the random number of Gaussian deviation 
!rrkf   with mean of 0 and one standard deviation of 1
   
    DO I=1, NLOC
      CALL GASDEV(IDUM,RNOBS)
      OBSDATA(I) = OBSDATA(I)+ OBSERR(STLOC(I))*RNOBS
    ENDDO

!rrkf   H(x^f)   -> moddata
!rrkf  Get the forecast values at the observation locations
    DO J=1,NLOC
      MODDATA(J)=STFCT(STLOC(J))
    ENDDO

!rrkf  Calculate innovation vector y' = y - H(x^f) ->moddata
   DO I = 1, NLOC
      MODDATA(I) = OBSDATA(I) - MODDATA(I)
   ENDDO
   
   TEXT='FCT(innv)'
!rrkf  Calculate the rms error
   ERR1 = 0.0_DP
   ERR2 = 0.0_DP
   ERR3 = 0.0_DP
   DO K = 1, NLOC
     IF(ABS(MODDATA(K)) > ERR1) THEN
        ERR1 = ABS(MODDATA(K))
     ENDIF
     ERR2 = ERR2 + MODDATA(K) * MODDATA(K)
     ERR3 = ERR3 + ABS(MODDATA(K))
   ENDDO
     ERR2 = DSQRT(ERR2/NLOC)
     ERR3 = ERR3 / NLOC
     ERR2_INN_FCT = ERR2
     
!******************************ended by pengfei************************* 

!rrkf  Also size of total error (full state vector)
!rrkf  Read the true state from 'ref2.cdf'
!   CALL READREF2(ICYC+(RRK_START - REF_TIME2)/DELTA_ASS-1)

!rrkf  Subtract B from A
   DO I=1,STDIM
     ERRVEC(I)=STFCT(I)-STTRUE1(I)
   ENDDO
   text='FCT(full)'
!rrkf  Calculate the rms error
   ERR1 = 0.
   ERR2 = 0.
   ERR3 = 0.
   DO K = 1, STDIM
     IF(ABS(ERRVEC(K)) > ERR1) THEN
        ERR1 = ABS(ERRVEC(K))
     ENDIF
     ERR2 = ERR2 + ERRVEC(K) * ERRVEC(K)
     ERR3 = ERR3 + ABS(ERRVEC(K))
   ENDDO
     ERR2 = DSQRT(ERR2/STDIM)
     ERR3 = ERR3 / STDIM
     ERR2_TOT_FCT = ERR2
!     IF(MSR) WRITE(IPT,*) TEXT,' max error=', ERR1, ' rms error=', &
!                     ERR2,' mean error=', ERR3
!*************************************************************************

!rrkf  Compute correction by applying gain matrix to innovation vector: Kr * (y - H(x^f))
   DO I=1,SS_DIM
     CRSTATE(I)=0.0_DP
     DO K=1,NLOC
        CRSTATE(I)= CRSTATE(I)+ KAL(I,K)*MODDATA(K)
     ENDDO
   ENDDO

!rrkf  Project this increment to the full state space(Er * Kr * (y-Hx^f)) and add it to x^f to get x^a

   DO I=1,STDIM
     RRKSUM=0.0_DP
     DO J=1,SS_DIM
       RRKSUM = RRKSUM + STSD(I)*SEOF(I,J)*CRSTATE(J)
     ENDDO

     STAN1(I) = STFCT(I) + RRKSUM
   ENDDO


!  WRITE(IPT,*) 'DA 5 stanl(1)',ICYC,STAN1(1)

!rrkf  Transfer x^a to the model variable, El, U, and V
    IDUMMY=0
    IF(EL_ASSIM) THEN
      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        RRKEL2(I)  = STAN1(IDUMMY) 
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN 
      DO I=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1
          RRKU2(J,I) = STAN1(IDUMMY) 
        ENDDO
      ENDDO
      DO I=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1
          RRKV2(J,I) = STAN1(IDUMMY)
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN 
      DO I=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1
          RRKT2(J,I) = STAN1(IDUMMY) 
        ENDDO
      ENDDO
    ENDIF
    IF(S_ASSIM) THEN 
      DO I=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1
          RRKS2(J,I) = STAN1(IDUMMY) 
        ENDDO
      ENDDO
    ENDIF
!  WRITE(IPT,*) 'DA 6 Elanl(1)',ICYC,El(1)

!   IF(MSR) CALL RRKF_INVERSE

   IF (SERIAL) THEN
     IF(EL_ASSIM) THEN
       DO I=1, NNodeGL
         EL(I) = RRKEL2(I)
       ENDDO
       D  = H + EL
       ET = EL
       DT = D
       DO I=1, MElemGL
         EL1(I)=(EL(NVG(I,1)) + EL(NVG(I,2)) + EL(NVG(I,3)) )/3.0_DP
       ENDDO
       D1  = H1 + EL1
       ET1 = EL1
       DT1 = D1
     ENDIF
     IF(UV_ASSIM) THEN  
       DO K=1, KB
         DO I=1, MElemGL
           U(I,K)  =  RRKU2(I,K)
           V(I,K)  =  RRKV2(I,K) 
         ENDDO
       ENDDO
!       DO I=1, MElemGL
!         UA(I)  =  RRKU2(I,1)
!         VA(I)  =  RRKV2(I,1) 
!       ENDDO
     ENDIF
     IF(T_ASSIM) THEN  
       DO K=1, KB
         DO I=1, NNodeGL
           T1(I,K)  =  RRKT2(I,K)
         ENDDO
       ENDDO
     ENDIF
     IF(S_ASSIM) THEN  
       DO K=1, KB
         DO I=1, NNodeGL
           S1(I,K)  =  RRKS2(I,K)
         ENDDO
       ENDDO
     ENDIF
   ELSE
#  if defined (MULTIPROCESSOR)
     IF(PAR) THEN
       IF(EL_ASSIM) THEN
         DO I=1, NNode
           EL(I)=RRKEL2(NGID(I))
         ENDDO
         D  = H + EL
         ET = EL
         DT = D
         DO I=1, MElem
           EL1(I)=(EL(NV(I,1)) + EL(NV(I,2)) + EL(NV(I,3)) )/3.0_DP
         ENDDO
         D1  = H1 + EL1
         ET1 = EL1
         DT1 = D1
       ENDIF
       IF(UV_ASSIM) THEN
         DO I=1, MElem
           DO J=1, KB
             U(I,J)=RRKU2(EGID(I),J)
             V(I,J)=RRKV2(EGID(I),J)
           ENDDO
         ENDDO
         DO I=1, MElem
           UA(I)=RRKU2(EGID(I),1)
           VA(I)=RRKV2(EGID(I),1)
         ENDDO
       ENDIF
       IF(T_ASSIM) THEN
         DO I=1, NNode
           DO J=1, KB
	     T1(I,J)=RRKT2(NGID(I),J)
	   ENDDO
	 ENDDO  
       ENDIF
       IF(S_ASSIM) THEN
         DO I=1, NNode
           DO J=1, KB
	     S1(I,J)=RRKS2(NGID(I),J)
	   ENDDO
	 ENDDO  
       ENDIF
     ENDIF
#  endif
   ENDIF 

#  if defined (WET_DRY)   

     CALL WET_JUDGE_EL

#  endif
  
!rrkf  Compute diagnostic on analysis 
!rrkf  Compute the analysis error in the observation space 
   DO J=1,NLOC
     MODDATA(J)=STAN1(STLOC(J))
   ENDDO
   DO K = 1, NLOC
      ERRVEC(K) = OBSDATA(K) - MODDATA(K)
   ENDDO
   TEXT='ANL(innv)'
   ERR1 = 0.0_DP
   ERR2 = 0.0_DP
   ERR3 = 0.0_DP
   DO K = 1, NLOC
     IF(ABS(ERRVEC(K)) > ERR1) THEN
        ERR1 = ABS(ERRVEC(K))
     ENDIF
     ERR2 = ERR2 + ERRVEC(K) * ERRVEC(K)
     ERR3 = ERR3 + ABS(ERRVEC(K))
   ENDDO
     ERR2 = DSQRT(ERR2/NLOC)
     ERR3 = ERR3 / NLOC
     ERR2_INN_AN1 = ERR2
!     IF(MSR) WRITE(IPT,*) TEXT,' max error=', ERR1, ' rms error=', &
!                     ERR2,' mean error=', ERR3

!******************************ended by pengfei*************************   

!rrkf  Compute the analysis error in the full state space 
   DO I=1,STDIM
      ERRVEC(I)=STAN1(I)-STTRUE1(I)
   ENDDO
   TEXT='ANL(full)'
   ERR1 = 0.
   ERR2 = 0.
   ERR3 = 0.
   DO K = 1, STDIM
     IF(ABS(ERRVEC(K)) > ERR1) THEN
        ERR1 = ABS(ERRVEC(K))
     ENDIF
     ERR2 = ERR2 + ERRVEC(K) * ERRVEC(K)
     ERR3 = ERR3 + ABS(ERRVEC(K))
   ENDDO
     ERR2 = DSQRT(ERR2/STDIM)
     ERR3 = ERR3 / STDIM
     ERR2_TOT_AN1 = ERR2
!     IF(MSR) WRITE(IPT,*) TEXT,' max error=', ERR1, ' rms error=', &
!                     ERR2,' mean error=', ERR3
!*******************************************************************

!rrkf  Write out error diagnostics for forecast and analysis
   IF(MSR) WRITE(OUTERR,'(I5,4E15.7)') ICYC,ERR2_INN_FCT, ERR2_INN_AN1, &
                                            ERR2_TOT_FCT, ERR2_TOT_AN1

   RETURN
   END SUBROUTINE RRK_RRKF



   SUBROUTINE SETUP_RRKA 

!------------------------------------------------------------------------------|
!  SETUP RRK_ASSIMILATION RUN                                                  |
!------------------------------------------------------------------------------|
   
   USE LIMS
   USE CONTROL
   USE MOD_RRK
   USE MOD_OBCS
   IMPLICIT NONE

   INTEGER I,J
   INTEGER SS_DIM
   INTEGER STDIM
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1

   SS_DIM = RRK_NVAR*RRK_NEOF

!rrkf  Set the first seed integer for the random number generation in gasdev, 
!rrkf  which will be used for the measurement error. 
   IDUM = -31

!READ THE OBSERVATION NUMBER AND LOCATION
   CALL READOBS2
   
!rrkf  Read the one standard deviation and set observation error'
   FILENAME = TRIM(OUTDIR)//'/rrktemp/avgstd.dat'
   OPEN(INORRK,FILE=FILENAME,STATUS='OLD')
   READ(INORRK,*)
   READ(INORRK,*) ELSD, USD, TSD, SSD
   CLOSE(INORRK)   
   
!rrkf  Set the magnitue of the measurement error as 1% of the spatially averaged one
!rrkf  standard deviation.
   CALL MAKEOBSERR   
   
!rrkf  Set a state vector consisting of the spatially averaged one standard deviation for later use.
   CALL MAKESTSD   
  
!rrkf  Read the retained eofs from 'eof.cdf'
   DO I=1, SS_DIM
      CALL READEOF2(I)
      DO J=1,STDIM
         SEOF(J,I)= STTEMP1(J)
      ENDDO
   ENDDO   
 
!rrkf  Read the Kalman gain matrix (Kr), which is creasted by 'fvcomrrK'. 
   FILENAME = TRIM(OUTDIR)//'/rrktemp/rrK.dat'
   OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
   DO J=1, NLOC
     READ(INORRK) (KAL(I,J), I=1, SS_DIM)    
   ENDDO 
   CLOSE(INORRK)
   
!rrkf  Set the file name of the assimilation results, forecast and analysis error
   ERRFILE = './ErrOut.dat'
   OPEN(OUTERR,FILE=ERRFILE,STATUS='UNKNOWN')   

!   I_INITIAL = IINT
   I_INITIAL = RRK_START
   NCYC      = (IEND-I_INITIAL)/DELTA_ASS
   N1CYC     = 1

   IF(MSR) THEN
      IF(IBCN_GL(1) > 0) CALL PERT_BC
   ENDIF
   
   RETURN
   END SUBROUTINE SETUP_RRKA



  SUBROUTINE ALLOC_RRKA

   USE LIMS
   USE CONTROL
   USE MOD_RRK
   IMPLICIT NONE
   
   INTEGER STDIM
   INTEGER SS_DIM

   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   
   SS_DIM = RRK_NVAR*RRK_NEOF

! ALLOCATE ARRYS

   ALLOCATE(STLOC(RRK_NOBSMAX))            ;STLOC     = 0 
   ALLOCATE(STFCT(STDIM))                  ;STFCT     = ZERO
   ALLOCATE(STAN1(STDIM))                  ;STAN1     = ZERO
   ALLOCATE(CRSTATE(SS_DIM))               ;CRSTATE   = ZERO
   ALLOCATE(KAL(SS_DIM,RRK_NOBSMAX))       ;KAL       = ZERO
   ALLOCATE(MODDATA(RRK_NOBSMAX))          ;MODDATA   = ZERO
   ALLOCATE(OBSDATA(RRK_NOBSMAX))          ;OBSDATA   = ZERO
   ALLOCATE(ERRVEC(STDIM))                 ;ERRVEC    = ZERO
   ALLOCATE(R(RRK_NOBSMAX,RRK_NOBSMAX))    ;R         = ZERO
   ALLOCATE(OBSERR(STDIM))                 ;OBSERR    = ZERO
   ALLOCATE(STTRUE(STDIM,(RRK_END-RRK_START)/DELTA_ASS+1))         ; STTRUE  = ZERO
   ALLOCATE(STTRUE1(STDIM))                ;STTRUE1   = ZERO
   ALLOCATE(STTR1(STDIM))                  ;STTR1     = ZERO
   ALLOCATE(STTR0(STDIM))                  ;STTR0     = ZERO
   ALLOCATE(SEOF(STDIM,SS_DIM))            ;SEOF      = ZERO
   ALLOCATE(STTEMP1(STDIM))                ;STTEMP1   = ZERO
   ALLOCATE(STSD(STDIM))                   ;STSD      = ZERO   
   ALLOCATE(RRKU2(MElemGL,KB))                 ;RRKU2     = ZERO
   ALLOCATE(RRKV2(MElemGL,KB))                 ;RRKV2     = ZERO
   ALLOCATE(RRKEL2(NNodeGL))                   ;RRKEL2    = ZERO
   ALLOCATE(RRKT2(NNodeGL,KB))                 ;RRKT2     = ZERO
   ALLOCATE(RRKS2(NNodeGL,KB))                 ;RRKS2     = ZERO
   
! END ALLOCATION
   
  RETURN
  END SUBROUTINE ALLOC_RRKA


  SUBROUTINE RRK_FNLOUT

   USE LIMS
   USE CONTROL
   USE MOD_RRK
   IMPLICIT NONE
      
   CHARACTER(LEN=80) :: FILNAME
   INTEGER I
   INTEGER STDIM
  
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   
   FILNAME = TRIM(OUTDIR)//'/rrktemp/forecast.cdf'
!   CALL PLOTSTATE_CDF(FILNAME,STFCT) 

   FILNAME = TRIM(OUTDIR)//'/rrktemp/analysis.cdf'
!   CALL PLOTSTATE_CDF(FILNAME,STAN1)
   
   DO I=1, STDIM
     STTEMP1(I) = STAN1(I) - STFCT(I)
   ENDDO       
   FILNAME = TRIM(OUTDIR)//'/rrktemp/incstate.cdf'
!   CALL PLOTSTATE_CDF(FILNAME,STTEMP1)
   
  RETURN
  END SUBROUTINE RRK_FNLOUT


!=======================================================================|
!rrkf  Set the magnitue of the measurement error as 1% of the spatially |
!rrkf  averaged one standard deviation.                                 | 
!=======================================================================|

   SUBROUTINE MAKEOBSERR
      
    USE LIMS
    USE CONTROL
    USE MOD_RRK
    IMPLICIT NONE

    INTEGER IDUMMY
    INTEGER I,J
    
    IDUMMY = 0
    
    IF(EL_ASSIM) THEN
      DO I =1, NNodeGL
        IDUMMY = IDUMMY + 1
        OBSERR(IDUMMY) = ELSD*DBLE(RRK_RSCALE)
      ENDDO
    ENDIF
    
    IF(UV_ASSIM) THEN
      DO J=1, 2*KBM1
        DO I=1, MElemGL
          IDUMMY = IDUMMY + 1
          OBSERR(IDUMMY) = USD*DBLE(RRK_RSCALE)
        ENDDO
      ENDDO
    ENDIF
    
    IF(T_ASSIM) THEN
      DO J=1, KBM1
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          OBSERR(IDUMMY)= TSD*DBLE(RRK_RSCALE)
        ENDDO
      ENDDO
    ENDIF
    
    IF(S_ASSIM) THEN
      DO J=1, KBM1
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          OBSERR(IDUMMY)= SSD*DBLE(RRK_RSCALE)
        ENDDO
      ENDDO
    ENDIF  
        
   RETURN
   END SUBROUTINE MAKEOBSERR

!=======================================================================!
!rrkf  Set a state vector consisting of the spatially averaged one      !
!rrkf  standard deviation for later use.                                !
!=======================================================================!
   SUBROUTINE MAKESTSD

    USE LIMS
    USE CONTROL
    USE MOD_RRK
    IMPLICIT NONE

    INTEGER IDUMMY
    INTEGER I,J    
    
    IDUMMY = 0
    
    IF(EL_ASSIM) THEN    
      DO I =1, NNodeGL
        IDUMMY = IDUMMY + 1
        STSD(IDUMMY) = ELSD   
      ENDDO
    ENDIF

    IF(UV_ASSIM) THEN
      DO J=1, 2*KBM1
        DO I=1, MElemGL
          IDUMMY = IDUMMY + 1
          IF(RRK_OPTION == 0) THEN
	     STSD(IDUMMY) = USD
	  ELSE
	     STSD(IDUMMY) = USD*DSQRT(DBLE(KBM1)) 
	  ENDIF     
        ENDDO
      ENDDO
    ENDIF
    
    IF(T_ASSIM) THEN
      DO J=1, KBM1
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          STSD(IDUMMY)= TSD
        ENDDO
      ENDDO
    ENDIF
    
    IF(S_ASSIM) THEN
      DO J=1, KBM1
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          STSD(IDUMMY)= SSD
        ENDDO
      ENDDO
    ENDIF    

   RETURN
   END SUBROUTINE MAKESTSD

!======================================================================  
!rrkf   Get the observations  
!======================================================================
   SUBROUTINE GETOBSDATA

    USE LIMS
    USE CONTROL
    USE MOD_RRK
    IMPLICIT NONE
    
    INTEGER J
    
    DO J=1,NLOC
      OBSDATA(J)= STTRUE(STLOC(J),(IINT-1-RRK_START)/DELTA_ASS+1)
    ENDDO

   RETURN
   END SUBROUTINE GETOBSDATA

!==========================================================================
!rrkf   Generate the random number of Gaussian deviation 
!rrkf   with mean of 0 and one standard deviation of 1
!==========================================================================
   SUBROUTINE GASDEV(IDUM,REV)
    
    USE LIMS
    USE CONTROL
    IMPLICIT NONE

    INTEGER IDUM
    REAL REV
    REAL REV2
!U    USES ran2
    INTEGER :: iset=0 
    REAL FAC
    REAL GSET
    REAL RSQ
    REAL v1
    REAL v2
    SAVE iset,gset

    IF(ISET==0) THEN
1     CALL RAN2(IDUM,REV2)
      V1=2.*REV2-1.
      CALL RAN2(IDUM,REV2)
      V2=2.*REV2-1.
      RSQ=V1**2+V2**2
      IF(RSQ >=1. .OR. RSQ == 0.) GOTO 1
      FAC=SQRT(-2.*LOG(RSQ)/RSQ)
      GSET=V1*FAC
      REV=V2*FAC
      ISET=1
    ELSE
      REV=GSET
      ISET=0
    ENDIF

   RETURN
   END SUBROUTINE GASDEV

!=====================================================================
  SUBROUTINE RAN2(IDUM,REV2)

    USE LIMS
    USE CONTROL
    IMPLICIT NONE
    
    INTEGER IDUM,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL AM
    REAL EPS
    REAL RNMX
    REAL REV2
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
               IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
               NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    
    IF(IDUM<=0) THEN
      IDUM=MAX(-IDUM,1)
      IDUM2=IDUM
      DO 11 J= NTAB+8,1,-1
        K=IDUM/IQ1
        IDUM=IA1*(IDUM-K*IQ1)-K*IR1
        IF (IDUM<0) IDUM=IDUM+IM1
        IF (J<=NTAB) IV(J)=IDUM
11    CONTINUE
      IY=IV(1)
    ENDIF
    K=IDUM/IQ1
    IDUM=IA1*(IDUM-K*IQ1)-K*IR1
    IF(IDUM<0) IDUM=IDUM+IM1
    K=IDUM2/IQ2
    IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
    IF(IDUM2<0) IDUM2=IDUM2+IM2
    J=1+IY/NDIV
    IY=IV(j)-IDUM2
    IV(J)=IDUM
    IF(IY<1) IY=IY+IMM1
    REV2 = MIN(AM*IY,RNMX)
    
  RETURN
  END SUBROUTINE RAN2

!==========================================================================|
!rrkf  Store the state vector for record and postprocess                   |   
!==========================================================================|
  SUBROUTINE PLOTSTATE_CDF(FILENAME,TEMP)

    USE LIMS
    USE CONTROL
    USE MOD_RRK
    IMPLICIT NONE

#include "/hosts/salmon01/data00/medm/src/netcdf-3.6.0-p1/src/fortran/netcdf.inc"   
!#include "/usr/local/include/nedcdf.inc" 
     INTEGER I,J,II
     INTEGER RCODE
     INTEGER START(2)
     INTEGER COUNT(2)
     INTEGER DIMS(2)
     INTEGER IDUMMY
     INTEGER IDVQ1, IDVQ2, IDVQ3, IDVQ4, IDVQ5
     INTEGER STDIM
     CHARACTER(LEN=80) :: FILENAME
     REAL(DP),ALLOCATABLE :: TEMP(:)
     REAL(DP),ALLOCATABLE :: ELTMP(:)
     REAL(DP),ALLOCATABLE :: UTMP(:,:)
     REAL(DP),ALLOCATABLE :: VTMP(:,:)
     REAL(DP),ALLOCATABLE :: TTMP(:,:)
     REAL(DP),ALLOCATABLE :: STMP(:,:)
    
     STDIM = 0
     IF(EL_ASSIM) STDIM = STDIM + NNodeGL
     IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
     IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
     IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1

     ALLOCATE(UTMP(MElemGL,KBM1))                ; UTMP     = ZERO
     ALLOCATE(VTMP(MElemGL,KBM1))                ; VTMP     = ZERO
     ALLOCATE(ELTMP(NNodeGL))                    ; ELTMP    = ZERO
     ALLOCATE(TTMP(NNodeGL,KBM1))                ; TTMP     = ZERO
     ALLOCATE(STMP(NNodeGL,KBM1))                ; STMP     = ZERO    

     IDUMMY = 0
     IF(EL_ASSIM) THEN
       DO I=1,NNodeGL
         IDUMMY = IDUMMY + 1
         ELTMP(I)= TEMP(IDUMMY)
       ENDDO
     ENDIF
     IF(UV_ASSIM) THEN
       DO I=1,KBM1
         DO J=1,MElemGL
           IDUMMY = IDUMMY + 1
           UTMP(J,I) = TEMP(IDUMMY)
         ENDDO
       ENDDO
       DO I=1,KBM1
         DO J=1,MElemGL
           IDUMMY = IDUMMY + 1
           VTMP(J,I) = TEMP(IDUMMY)
         ENDDO
       ENDDO
     ENDIF   
     IF(T_ASSIM) THEN  
       DO I=1,KBM1
         DO J=1,NNodeGL
           IDUMMY = IDUMMY + 1
           TTMP(J,I) = TEMP(IDUMMY)
         ENDDO
       ENDDO
     ENDIF
     IF(S_ASSIM) THEN  
       DO I=1,KBM1
         DO J=1,NNodeGL
           IDUMMY = IDUMMY + 1
           STMP(J,I) = TEMP(IDUMMY)
         ENDDO
       ENDDO
     ENDIF     
! Output state in cdf file
      RCODE = nf_create(FILENAME,NF_CLOBBER,II)

      IF(EL_ASSIM) THEN
        RCODE = nf_def_dim(II,'n_el',NNodeGL,DIMS(1))
        RCODE = nf_def_var(II,'el',NF_DOUBLE,1,DIMS,IDVQ1)
      ENDIF
      IF(UV_ASSIM) THEN 
        RCODE = nf_def_dim(II,'c_u',MElemGL,DIMS(1))
        RCODE = nf_def_dim(II,'z_u',KBM1,DIMS(2))
        RCODE = nf_def_var(II,'u',NF_DOUBLE,2,DIMS,IDVQ2)
        RCODE = nf_def_var(II,'v',NF_DOUBLE,2,DIMS,IDVQ3)
      ENDIF 
      IF(T_ASSIM) THEN
        RCODE = nf_def_dim(II,'c_temp',NNodeGL,DIMS(1))
	RCODE = nf_def_dim(II,'z_temp',KBM1,DIMS(2))
        RCODE = nf_def_var(II,'temp',NF_DOUBLE,2,DIMS,IDVQ4)
      ENDIF
      IF(S_ASSIM) THEN
        RCODE = nf_def_dim(II,'c_sal',NNodeGL,DIMS(1))
	RCODE = nf_def_dim(II,'z_sal',KBM1,DIMS(2))
        RCODE = nf_def_var(II,'sal',NF_DOUBLE,2,DIMS,IDVQ5)
      ENDIF
      
      RCODE = nf_enddef(II)

      IF(EL_ASSIM) THEN 
        START(1)=1
        START(2)=1
        COUNT(1)= NNodeGL
        COUNT(2)= 1
        RCODE = nf_put_var_double(II,IDVQ1,ELTMP)
      ENDIF

      IF(UV_ASSIM) THEN
        START(1)=1
        COUNT(1)= MElemGL
        START(2)=1
        COUNT(2)= KBM1
        RCODE = nf_put_var_double(II,IDVQ2,UTMP)

        RCODE = nf_put_var_double(II,IDVQ3,VTMP)
      ENDIF
      
      IF(T_ASSIM) THEN
        START(1)=1
        COUNT(1)= NNodeGL
        START(2)=1
        COUNT(2)= KBM1
        RCODE = nf_put_var_double(II,IDVQ4,TTMP)
      ENDIF
      
      IF(S_ASSIM) THEN
        START(1)=1
        COUNT(1)= NNodeGL
        START(2)=1
        COUNT(2)= KBM1
        RCODE = nf_put_var_double(II,IDVQ5,STMP)
      ENDIF

      RCODE = nf_close(II)
  
    DEALLOCATE(UTMP,VTMP,ELTMP,TTMP,STMP)
  
  RETURN
  END SUBROUTINE PLOTSTATE_CDF

!=====================================================================================/
!  READ THE EOFs FROM eof.cdf                                                         /
!=====================================================================================/
  SUBROUTINE READEOF2(IEOF)
 
   USE LIMS
   USE MOD_RRK
   IMPLICIT NONE

#include "/hosts/salmon01/data00/medm/src/netcdf-3.6.0-p1/src/fortran/netcdf.inc"   
!#include "/usr/local/include/nedcdf.inc" 
    INTEGER I,J,K
    INTEGER RCODE
    INTEGER IEOF
    INTEGER START(3)
    INTEGER COUNT(3)
    INTEGER STATUS
    INTEGER NCID
    INTEGER VARID
    INTEGER IDUMMY
    CHARACTER(LEN=80) FNAME
    REAL(DP),ALLOCATABLE :: UTMP(:,:)
    REAL(DP),ALLOCATABLE :: VTMP(:,:)
    REAL(DP),ALLOCATABLE :: ELTMP(:)
    REAL(DP),ALLOCATABLE :: TTMP(:,:)
    REAL(DP),ALLOCATABLE :: STMP(:,:)

    ALLOCATE(UTMP(MElemGL,KBM1))         ; UTMP  = ZERO
    ALLOCATE(VTMP(MElemGL,KBM1))         ; VTMP  = ZERO
    ALLOCATE(ELTMP(NNodeGL))             ; ELTMP = ZERO 
    ALLOCATE(TTMP(NNodeGL,KBM1))         ; TTMP  = ZERO
    ALLOCATE(STMP(NNodeGL,KBM1))         ; STMP  = ZERO
    
    FNAME=TRIM(OUTDIR)//'/rrktemp/'//'eof.cdf'
    STATUS = nf_open(FNAME,NF_NOWRITE,NCID)
    
    IF(EL_ASSIM) THEN
      STATUS = nf_inq_varid(NCID,'eof_el',VARID)
    
      START(1) = 1
      COUNT(1) = NNodeGL
      START(2) = IEOF
      COUNT(2) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,ELTMP)
    
      IDUMMY = 0
      DO J=1, NNodeGL
        IDUMMY = IDUMMY +1
        STTEMP1(IDUMMY) = ELTMP(J)
      ENDDO
    ENDIF     
 
    IF(UV_ASSIM) THEN
      STATUS = nf_inq_varid(NCID,'eof_u',VARID)
     
      START(1) = 1
      COUNT(1) = MElemGL
      START(2) = 1
      COUNT(2) = KBM1
      START(3) = IEOF
      COUNT(3) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,UTMP)
    
      DO K=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1 
          STTEMP1(IDUMMY) = UTMP(J,K)
        ENDDO
      ENDDO

      STATUS = nf_inq_varid(NCID,'eof_v',VARID)
     
      START(1) = 1
      COUNT(1) = MElemGL
      START(2) = 1
      COUNT(2) = KBM1
      START(3) = IEOF
      COUNT(3) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,VTMP)
    
      DO K=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1 
          STTEMP1(IDUMMY) = VTMP(J,K)
        ENDDO
      ENDDO
    ENDIF

    IF(T_ASSIM) THEN
      STATUS = nf_inq_varid(NCID,'eof_temp',VARID)
     
      START(1) = 1
      COUNT(1) = NNodeGL
      START(2) = 1
      COUNT(2) = KBM1
      START(3) = IEOF
      COUNT(3) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,TTMP)
    
      DO K=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1 
          STTEMP1(IDUMMY) = TTMP(J,K)
        ENDDO
      ENDDO   
    ENDIF     

    IF(S_ASSIM) THEN
      STATUS = nf_inq_varid(NCID,'eof_sal',VARID)
     
      START(1) = 1
      COUNT(1) = NNodeGL
      START(2) = 1
      COUNT(2) = KBM1
      START(3) = IEOF
      COUNT(3) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,STMP)
    
      DO K=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1 
          STTEMP1(IDUMMY) = STMP(J,K)
        ENDDO
      ENDDO   
    ENDIF
    
    RCODE = nf_close(NCID)

    DEALLOCATE(UTMP,VTMP,ELTMP,TTMP,STMP)
    
    RETURN
  END SUBROUTINE READEOF2

!======================================================================!
!                                                                      !
!======================================================================!
  SUBROUTINE READOBS2
    
   USE LIMS
   USE CONTROL
   USE MOD_RRK
   IMPLICIT NONE
   
     INTEGER ::  NUM  = 0 
     INTEGER ::  SWITCH = 0
     INTEGER ::  J,K
     INTEGER ::  IDUMMY = 0
     INTEGER ::  TMP
     CHARACTER(LEN=80) FILENAME
     CHARACTER(LEN=24) HEADINFO
     INTEGER LAY(RRK_NOBSMAX)

     FILENAME = TRIM(INPDIR)//"/"//trim(casename)//"_assim_rrkf.dat"
    
     OPEN(INORRK,FILE=TRIM(FILENAME),FORM='FORMATTED')

     NLOC = 0

 100 READ(INORRK,'(A24)',END=200) HEADINFO
     IF(SWITCH/=1) THEN
       IF(HEADINFO=='!=== READ IN OBSERVATION') THEN
         SWITCH = 1
         GOTO 100
       ELSE
         GOTO 100
       ENDIF
     ENDIF 
     
     IF(TRIM(HEADINFO)=='!EL') THEN
       IF(EL_OBS) THEN
         READ(INORRK,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INORRK,*)  (STLOC(K), K=1,NLOC)	 
	 
       ENDIF

       IF(EL_ASSIM) THEN
	 IDUMMY = IDUMMY + NNodeGL
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!UV') THEN
       IF(UV_OBS) THEN
         READ(INORRK,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC+NUM>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc+num, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INORRK,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)
	 READ(INORRK,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MElemGL*(LAY(K)-1)
	 ENDDO   
	 
	 NLOC = NLOC + NUM
	 DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K-NUM)+MElemGL*KBM1+MElemGL*(LAY(K-NUM)-1)
	 ENDDO
	 
       ENDIF
       IF(UV_ASSIM) THEN  
          IDUMMY = IDUMMY + MElemGL*KBM1
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!T') THEN
       IF(T_OBS) THEN
         READ(INORRK,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INORRK,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)       
         READ(INORRK,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NNodeGL*(LAY(K)-1)
	 ENDDO   
         
       ENDIF   
       IF(T_ASSIM) THEN
         IDUMMY = IDUMMY + NNodeGL*KBM1
       ENDIF
     ENDIF

     IF(TRIM(HEADINFO)=='!S') THEN
       IF(S_OBS) THEN
         READ(INORRK,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INORRK,*)  (STLOC(K),K=NLOC-NUM+1,NLOC)        
         READ(INORRK,*)  (LAY(K),  K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NNodeGL*(LAY(K)-1)
	 ENDDO   
        
       ENDIF
       IF(S_ASSIM) THEN
         IDUMMY = IDUMMY + NNodeGL*KBM1 
       ENDIF 
     ENDIF
     
     GOTO 100
 200 CONTINUE

     DO J=1, NLOC-1
       DO K=2, NLOC
	 IF(STLOC(K)<STLOC(J)) THEN
	   TMP = STLOC(J)
	   STLOC(J) = STLOC(K)
	   STLOC(K) = TMP 
         ENDIF
       ENDDO
     ENDDO

     CLOSE(INORRK)
    
  RETURN 
  END SUBROUTINE READOBS2  

  SUBROUTINE RRKF_INVERSE

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   USE MOD_RRK
   IMPLICIT NONE
   
   INTEGER INUM,I,J,K
   INTEGER NCON
   INTEGER,ALLOCATABLE :: NODE_SBC(:)
   INTEGER :: NCON_FLG(6)
   
!UPDATE THE ELEVATION TIMESER

   OPEN(INORRK,FILE=TRIM(OUTDIR)//'/rrktemp/el_srs.dat')
   DO J=1,TIMEN 
     READ(INORRK,*) (EL_SRS(I,J), I=1,IBCN_GL(1))
   ENDDO
   CLOSE(INORRK)

   DO I=1,IBCN_GL(1)
     EL_SRS(I,ICYC) = RRKEL2(I_OBC_GL(I))             
   ENDDO

   OPEN(INORRK,FILE=TRIM(OUTDIR)//'/rrktemp/el_srs.dat') 
   DO J=1,TIMEN 
     WRITE(INORRK,'(1000(F13.5))') (EL_SRS(I,J),I=1,IBCN_GL(1))
   ENDDO	
   CLOSE(INORRK)

   OPEN(INORRK,FILE=TRIM(OUTDIR)//'/rrktemp/bc.dat',STATUS='UNKNOWN')
   DEALLOCATE(EMEAN,APT,PHAI)
   ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
   ALLOCATE(APT(IBCN_GL(1),6), PHAI(IBCN_GL(1),6))
   APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP
   DO I=1,IBCN_GL(1)
      READ(INORRK,*)  NODE_SBC(I),EMEAN(I)
      READ(INORRK,*) 
      READ(INORRK,*) 
   ENDDO
   CLOSE(INORRK)

!ADJUST AMPLITUDE
     OPEN(INORRK,FILE=TRIM(OUTDIR)//'/rrktemp/el_srs.dat',STATUS='OLD')
     DO J=1,TIMEN
       READ(INORRK,*) (EL_SRS(I,J),I=1,IBCN_GL(1))
     ENDDO	
     CLOSE(INORRK)
     
     NCON     = 0
     NCON_FLG = 0
     DO I=1, 6	
        IF(BC_AMP_ERR(I)>0.0001 .AND. BC_PHA_ERR(I)>0.0001 ) THEN
	   NCON = NCON + 1
	   NCON_FLG(I) = 1
        ENDIF
     ENDDO

     DO I=1,IBCN_GL(1)

       ALLOCATE(SRS_TMP(TIMEN))     ; SRS_TMP  = 0.0_DP 
       SRS_TMP = EL_SRS(I,:)*100.0_DP
       CALL LSQFIT(SRS_TMP,TIME_SER,TIMEN,NCON,NCON_FLG,AMP,PHAS)
       DEALLOCATE(SRS_TMP)

       K = 0
       DO J=1, 6
         IF (NCON_FLG(J)==1) THEN
            K = K + 1
	    APT(I,J)  = AMP(K)
	    PHAI(I,J) = PHAS(K)
         ENDIF
       ENDDO

     ENDDO
     PHAI = MOD(PHAI,360.0_SP)

!     OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/rrktemp/amppha1_'//FNUM//'.dat',STATUS='UNKNOWN',ACCESS='APPEND')
!     WRITE(IOBCKF,*) ICYC, APT(1,2), PHAI(1,2)
!     CLOSE(IOBCKF)

     OPEN(INORRK,FILE=TRIM(OUTDIR)//'/rrktemp/bc.dat',STATUS='UNKNOWN')
     DO I=1,IBCN_GL(1)
       WRITE(INORRK,'(I10,1000(F13.5))')  NODE_SBC(I),EMEAN(I)
       WRITE(INORRK,'(1000(F13.5))') (APT(I,J), J=1,6)
       WRITE(INORRK,'(1000(F13.5))') (PHAI(I,J), J=1,6)
     ENDDO
     CLOSE(INORRK)

!RE-BUILD ELEVATION TIME SERIES
     OPEN(INORRK,FILE=TRIM(OUTDIR)//'/rrktemp/el_srs.dat',STATUS='REPLACE')
     DO I =1, IBCN_GL(1)
       DO J = ICYC+1, TIMEN
         FORCE = 0.0_SP
         DO K = 1, 6
           PHAI_IJ= PHAI(I,K)*PI2/360.0_SP
           FORCE  = FORCE + APT(I,K)/100.0_DP * COS(PI2/PERIOD(K)*(DTI*FLOAT(RRK_START+(J-1)*DELTA_ASS)) -PHAI_IJ)
!          FORCE  = FORCE + APT(I,K)/100.0_DP * COS(PI2/PERIOD(K)*TIME_SER(J) -PHAI_IJ)    
         ENDDO
         EL_SRS(I,J)= FORCE + EMEAN(I)       
       ENDDO
     ENDDO
     DO J=1, TIMEN  
       WRITE(INORRK,'(1000(F13.5))') (EL_SRS(I,J),I=1,IBCN_GL(1))
     ENDDO	
     CLOSE(INORRK)

     DEALLOCATE(NODE_SBC)
     RETURN 
   END SUBROUTINE RRKF_INVERSE

   
   SUBROUTINE PERT_BC

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   USE MOD_RRK
   IMPLICIT NONE
   
   INTEGER I,J,K,KK
   CHARACTER(LEN=4) FNUM 
   CHARACTER(LEN=80)   :: ISTR
   CHARACTER(LEN=80)   :: COMT
   INTEGER,ALLOCATABLE :: NODE_SBC(:)
   REAL(DP),ALLOCATABLE :: APTTMP(:,:),PHAITMP(:,:)
   INTEGER   ISBCN1
   
   
   IF(IBCN_GL(1) > 0)THEN

   IF(S_TYPE == 'non-julian') THEN

     ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)
     OPEN(INORRK,FILE=TRIM(ISTR)//'_el_obc.dat')

     READ(INORRK ,1000) COMT
     READ(INORRK,*) ISBCN1
!
!-------ENSURE SAME NUMBER OF SPECIFIED OPEN BOUNDARY POINTS AS FILE-casename_obc.dat----|
!
     IF(ISBCN1 /= IBCN_GL(1))THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'NUMBER OF OPEN BOUNDARY POINTS IN OPEN BOUNDARY SURFACE'
       WRITE(IPT,*)'ELEVATION FILE IS LARGER THAN NUMBER OF OPEN BOUNDARY '
       WRITE(IPT,*)'POINTS OF PRESCRIBED ELEVATION TYPE IN CASENAME_obc.dat'
       WRITE(IPT,*) 'SEE SUBROUTINE BCS_FORCE'
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

!
!----READ IN BOUNDARY POINTS, AMPLITUDES, AND PHASES OF TIDE-------------------|
!
     DEALLOCATE(EMEAN,APT,PHAI)
     ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
     ALLOCATE(APT(IBCN_GL(1),6), PHAI(IBCN_GL(1),6))
     ALLOCATE(APTTMP(IBCN_GL(1),6), PHAITMP(IBCN_GL(1),6))
     APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP
     APTTMP = 0.0_SP; PHAITMP = 0.0_SP
     DO I=1,IBCN_GL(1)
       READ(INORRK,*)  NODE_SBC(I), EMEAN(I)
       READ(INORRK,*) (APT(I,J), J=1,6)
       READ(INORRK,*) (PHAI(I,J), J=1,6)
     END DO

     PHAI = MOD(PHAI,360.0_SP)
     APTTMP  = APT
     PHAITMP = PHAI 

     CLOSE(INORRK)

   ELSE
     WRITE(IPT,*) 'INVERSE METHOD CAN ONLY BE USED FOR NON-JULIAN TIDAL SIMULATION RIGHT NOW!'
     CALL PSTOP   
   ENDIF  
   
   ELSE
     WRITE(IPT,*) 'NO TIDAL B.C.s ARE SPECIFIED TO DO INVERSE, PLEASE CHECK AGAIN!'
     CALL PSTOP
   ENDIF 

   TIMEN = (RRK_END - RRK_START)/DELTA_ASS+1
   ALLOCATE(EL_SRS(IBCN_GL(1),TIMEN))    ; EL_SRS   = 0.0_DP
   ALLOCATE(TIME_SER(TIMEN))             ; TIME_SER = 0.0_DP
        
   DO I=1,TIMEN 
     TIME_SER(I)=(DTI*DBLE(RRK_START+(I-1)*DELTA_ASS)) 
   ENDDO

   DO I=1, IBCN_GL(1)
     DO J=1, 6
       APT(I,J)  =  APTTMP(I,J) + BC_AMP_ERR(J) 
       PHAI(I,J) = PHAITMP(I,J) + BC_PHA_ERR(J)
     ENDDO	   
   ENDDO

   PHAI = MOD(PHAI,360.0_SP)

   OPEN(INORRK,FILE=TRIM(OUTDIR)//'/rrktemp/'//'el_srs.dat',STATUS='REPLACE')
   DO I = 1,IBCN_GL(1)
      DO J = 1, TIMEN 
         FORCE = 0.0_SP
         DO KK = 1, 6
            PHAI_IJ= PHAI(I,KK)*PI2/360.0_SP 
            FORCE  = FORCE + APT(I,KK)/100._DP * COS(PI2/PERIOD(KK)*(DTI*FLOAT(RRK_START+(J-1)*DELTA_ASS)) - PHAI_IJ)
         ENDDO
         EL_SRS(I,J) = FORCE + EMEAN(I)
      END DO
   END DO
   DO J=1, TIMEN
      WRITE(INORRK,'(1000(F13.5))') (EL_SRS(I,J),I=1,IBCN_GL(1))
   ENDDO	
   CLOSE(INORRK)

   OPEN(INORRK,FILE=TRIM(OUTDIR)//'/rrktemp/bc.dat')
   DO I=1, IBCN_GL(1)
      WRITE(INORRK,'(I10,1000(F13.5))') NODE_SBC(I),EMEAN(I)
      WRITE(INORRK,'(1000(F13.5))')     ( APT(I,J),J=1,6)
      WRITE(INORRK,'(1000(F13.5))')     (PHAI(I,J),J=1,6)
   ENDDO
   CLOSE(INORRK) 

   DEALLOCATE(NODE_SBC,APTTMP,PHAITMP)
   1000 FORMAT(A80)
   RETURN   
   END SUBROUTINE PERT_BC

   SUBROUTINE READ_BC
   
   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER I,J,K
   CHARACTER(LEN=4) FNUM 
   CHARACTER(LEN=80)   :: ISTR
   CHARACTER(LEN=80)   :: COMT
   INTEGER,ALLOCATABLE :: NODE_SBC(:)
   INTEGER   ISBCN1,NCNT,JN
   REAL(SP), ALLOCATABLE :: RTEMP(:),RTEMP1(:,:),RTEMP2(:,:)
   INTEGER,  ALLOCATABLE :: TEMP2(:)
      
   ISTR = TRIM(OUTDIR)//'/rrktemp/bc.dat'
   OPEN(INORRK,FILE=TRIM(ISTR))
   
!----READ IN BOUNDARY POINTS, AMPLITUDES, AND PHASES OF TIDE-------------------|
!
     DEALLOCATE(EMEAN,APT,PHAI)
     ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
     ALLOCATE(APT(IBCN_GL(1),6), PHAI(IBCN_GL(1),6))
     APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP
     DO I=1,IBCN_GL(1)
       READ(INORRK,*)  NODE_SBC(I),EMEAN(I)
       READ(INORRK,*) (APT(I,J), J=1,6)
       READ(INORRK,*) (PHAI(I,J), J=1,6)
     END DO
   CLOSE(INORRK)  

#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
     ALLOCATE( TEMP2(IBCN_GL(1)) ,RTEMP(IBCN_GL(1)))
     ALLOCATE( RTEMP1(IBCN_GL(1),6) , RTEMP2(IBCN_GL(1),6))
     NCNT = 0
     DO I=1,IBCN_GL(1)
       IF(NLID(NODE_SBC(I)) /= 0)THEN
         NCNT = NCNT + 1
         TEMP2(NCNT)     = NLID(NODE_SBC(I))
         RTEMP(NCNT)     = EMEAN(I)
         RTEMP1(NCNT,1:6) = APT(I,1:6)
         RTEMP2(NCNT,1:6) = PHAI(I,1:6)
       END IF
     END DO

     IF(NCNT /= IBCN(1))THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'LOCAL OPEN BOUNDARY NODE COUNTS DIFFER BETWEEN TIDE'
       WRITE(IPT,*)'FORCING AND OPEN BOUNDARY NODE FILES'
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     DEALLOCATE(NODE_SBC,EMEAN,APT,PHAI)
     IF(IBCN(1) > 0)THEN
       ALLOCATE(NODE_SBC(IBCN(1)),EMEAN(IBCN(1)))
       ALLOCATE(APT(IBCN(1),6),PHAI(IBCN(1),6))
       NODE_SBC = TEMP2(1:NCNT)
       EMEAN    = RTEMP(1:NCNT)
       APT      = RTEMP1(1:NCNT,1:6)
       PHAI     = RTEMP2(1:NCNT,1:6)
     ELSE
       ALLOCATE(NODE_SBC(1),EMEAN(1))
       ALLOCATE(APT(1,6),PHAI(1,6))
       NODE_SBC = 0.0_SP ; EMEAN = 0.0_SP ; APT = 0.0_SP ; PHAI = 0.0_SP
     END IF

     DEALLOCATE(TEMP2,RTEMP,RTEMP1,RTEMP2)

     END IF !!PAR
#    endif

!
!----MAKE SURE LOCAL NODE NUMBERS OF SPECIFIED NODES MATCHES LOCAL NODE--------|
!----NUMBER OF SPECIFIED NODES IN obc.dat FILE---------------------------------|
!
     DO I=1,IBCN(1)
       JN = OBC_LST(1,I)
       IF(NODE_SBC(I) /= I_OBC_N(JN))THEN
         WRITE(IPT,*)'==================ERROR=================================='
         WRITE(IPT,*)'LOCAL OPEN BOUNDARY NODE LIST DIFFERS BETWEEN TIDE'
         WRITE(IPT,*)'FORCING AND OPEN BOUNDARY NODE (TYPE 1 OR 2) FILES'
         WRITE(IPT,*)'========================================================='
         WRITE(IPT,*)NODE_SBC(I),I_OBC_N(JN)
         CALL PSTOP
       END IF
     END DO

     APT = APT/100.0_SP
     PHAI = MOD(PHAI,360.0_SP)
   
     DEALLOCATE(NODE_SBC)
     RETURN
   END SUBROUTINE READ_BC

      SUBROUTINE LSQFIT(EL,T,NUM,NCOMP,IDX,AMP2,PHA2)
!-----------------------------------------------------------------------
!     THIS SUBROUTINE IS USED FOR ANALYSIS AMPLITUDE AND PHASE BY LSQFIT
!
!     INPUT: 
!         EL(NUM)---ELEVATION (m)
!         T(NUM) ---TIME (s)
!         NUM    ---NUMBER OF TIMESERIES
!     OUTPUT:
!         AMP(6) ---AMPPLITUDE (m)
!         PHA(6) ---PHASE (deg)
!     from QXU, modified by LZG 10/25/05, 
!-----------------------------------------------------------------------
      IMPLICIT NONE
!      INTEGER, PARAMETER :: NCOMP=6, NCOMP2 = NCOMP*2+1
      REAL, PARAMETER, DIMENSION(6) :: &
                   !  S2       M2       N2       K1       P1       O1 
      PERIOD2 = (/43200.0, 44712.0, 45570.0, 86164.0, 86637.0, 92950.0/) !(sec)
              != 12.0000  12.4200  12.6583  23.9344  24.0658  25.8194 hours
      REAL, PARAMETER :: PI = 3.1415926
      
      INTEGER NCOMP, NCOMP2, IDX(6)
      INTEGER NUM,N,I,J,K,I1,I2,J1,J2
      REAL*8, DIMENSION(NUM) :: EL
      REAL*8, DIMENSION(NUM) :: T
      REAL*8  STEL,AEL
      REAL*8 A(NCOMP*2+1,NCOMP*2+1),B(NCOMP*2+1),F(NCOMP)
      REAL*8 AMP1(NCOMP),PHA1(NCOMP)
      REAL*8 AMP2(NCOMP),PHA2(NCOMP)
      REAL*8 PERIOD(NCOMP)
      
!      F = 2.0*PI/(PERIOD/3600.0)        !(1/HOUR)
!      F = 2.0*PI/PERIOD                 !(1/s)

      NCOMP2 = NCOMP*2+1
      
      N=0
      DO I=1, 6
         IF(IDX(I)==1) THEN       
           N = N + 1
	   PERIOD(N) = PERIOD2(I)
	 ENDIF
      ENDDO      
      
      F = 2.0*PI/PERIOD                  !(1/s)
     
      AEL = 0.0
      DO N=1,NUM
         AEL = AEL +EL(N)
      ENDDO
      AEL = AEL/FLOAT(NUM) 	 
      DO N=1,NUM
         EL(N)=EL(N)-AEL
      ENDDO
      STEL=0.0
      DO N=1,NUM
         STEL=STEL+EL(N)*EL(N)	 
      ENDDO
      STEL = SQRT(STEL/FLOAT(NUM))
      DO N=1,NUM
         EL(N)=EL(N)/STEL
      ENDDO

    
      DO J = 1, NCOMP2
         DO K = 1, NCOMP2
            A(J,K) = 0.0
         ENDDO
      ENDDO
      DO J = 1, NCOMP2
         B(J) = 0.0
      ENDDO
           
      DO N = 1,NUM
         A(1,1)    = A(1,1)    + 1
	 DO I=1,NCOMP
	    I1 = I*2
	    I2 = I1+1
	    A(1,I1)= A(1,I1)   + COS(F(I)*T(N))
	    A(1,I2)= A(1,I2)   + SIN(F(i)*T(N))
	 ENDDO 
	 DO I=1,NCOMP  
	    I1 = I*2
	    I2 = I1+1
	    DO J=I,NCOMP
	       J1 = J*2
	       J2 = J1+1
	       A(I1,J1) = A(I1,J1) + COS(F(I)*T(N))* COS(F(J)*T(N))
	       A(I1,J2) = A(I1,J2) + COS(F(I)*T(N))* SIN(F(J)*T(N))
	       A(I2,J1) = A(I2,J1) + SIN(F(I)*T(N))* COS(F(J)*T(N))
	       A(I2,J2) = A(I2,J2) + SIN(F(I)*T(N))* SIN(F(J)*T(N))
	    ENDDO   
	 ENDDO      
         
         B(1) = B(1) + EL(N)
	 DO I=1,NCOMP
	    I1 = I*2
	    I2 = I1+1
	    B(I1) = B(I1) + EL(N)*COS(F(I)*T(N))
	    B(I2) = B(I2) + EL(N)*SIN(F(I)*T(N))
	 ENDDO   
      ENDDO
      DO I=2,NCOMP2
         DO J=1,I
            A(I,J)=A(J,I)
         ENDDO
      ENDDO
         
      CALL GAUSSJ_2(A,NCOMP2,NCOMP2,B,1,1)
      
      DO I=1,NCOMP
         I1 = I*2
	 I2 = I1+1
         AMP1(I) = SQRT(B(I1)*B(I1)+B(I2)*B(I2))*STEL
	 PHA1(I) = ATAN2(B(I2),B(I1))*180/PI
         IF(PHA1(I).LT.0) PHA1(I)=PHA1(I)+360.0
      ENDDO  
      AMP2 = AMP1
      PHA2 = PHA1 	   
          
      RETURN
      END SUBROUTINE LSQFIT
     

!--------------------------------------------------      
      SUBROUTINE GAUSSJ_2(A,N,NP,B,M,MP)
      IMPLICIT NONE
      INTEGER M,MP,N,NP,NMAX
      REAL*8 A(NP,NP),B(NP,MP)
      PARAMETER (NMAX=50)
      INTEGER I,ICOL,IROW,J,K,L,LL,INDXC(NMAX),INDXR(NMAX),&
             IPIV(NMAX)
      REAL*8 BIG,DUM,PIVINV
      DO J=1,N
         IPIV(J)=0
      ENDDO
      DO 22 I=1,N
         BIG=0.
         DO 13 J=1,N
            IF(IPIV(J).NE.1)THEN
              DO 12 K=1,N
                 IF (IPIV(K).EQ.0) THEN
                   IF (ABS(A(J,K)).GE.BIG)THEN
                     BIG=ABS(A(J,K))
                     IROW=J
                     ICOL=K
                   ENDIF

                 ELSE IF (IPIV(K).GT.1) THEN
                   PAUSE 'SINGULAR MATRIX IN GAUSSJ'
                 ENDIF
12            CONTINUE
            ENDIF
13       CONTINUE
         IPIV(ICOL)=IPIV(ICOL)+1
         IF (IROW.NE.ICOL) THEN
           DO 14 L=1,N
              DUM=A(IROW,L)
              A(IROW,L)=A(ICOL,L)
              A(ICOL,L)=DUM
14         CONTINUE
           DO 15 L=1,M
              DUM=B(IROW,L)
              B(IROW,L)=B(ICOL,L)
              B(ICOL,L)=DUM
15         CONTINUE
         ENDIF
         INDXR(I)=IROW
         INDXC(I)=ICOL
         IF (A(ICOL,ICOL).EQ.0.) PAUSE 'SINGULAR MATRIX IN GAUSSJ'

         PIVINV=1./A(ICOL,ICOL)
         A(ICOL,ICOL)=1.
         DO 16 L=1,N
            A(ICOL,L)=A(ICOL,L)*PIVINV
16       CONTINUE
         DO 17 L=1,M
            B(ICOL,L)=B(ICOL,L)*PIVINV
17       CONTINUE
         DO 21 LL=1,N
            IF(LL.NE.ICOL)THEN
              DUM=A(LL,ICOL)
              A(LL,ICOL)=0.
              DO 18 L=1,N
                 A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18            CONTINUE
              DO 19 L=1,M
                 B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19            CONTINUE
            ENDIF
21       CONTINUE
22    CONTINUE

      DO 24 L=N,1,-1
         IF(INDXR(L).NE.INDXC(L))THEN

           DO 23 K=1,N
              DUM=A(K,INDXR(L))
              A(K,INDXR(L))=A(K,INDXC(L))
              A(K,INDXC(L))=DUM
23         CONTINUE
         ENDIF
24    CONTINUE

      RETURN
      END SUBROUTINE GAUSSJ_2
      
!=====================================================================
!rrkf  Read the true state from 'ref2.cdf'
!=====================================================================
   SUBROUTINE READREF2(INDEX)
    
    USE LIMS
    USE CONTROL
    USE MOD_RRK
    IMPLICIT NONE

!#include "/usr/medm/include/netcdf.inc"
#include "/hosts/salmon01/data00/medm/src/netcdf-3.6.0-p1/src/fortran/netcdf.inc"   
     INTEGER INDEX
     INTEGER IDUMMY
     INTEGER J,K
     INTEGER NCINP
     INTEGER VINPID
     INTEGER VARID
     INTEGER STATUS
     INTEGER RCODE
     INTEGER START(3)
     INTEGER COUNT(3)
     CHARACTER(LEN=80) FNAME
     REAL(DP),ALLOCATABLE :: RKU(:,:)
     REAL(DP),ALLOCATABLE :: RKV(:,:)
     REAL(DP),ALLOCATABLE :: RKEL(:)
     REAL(DP),ALLOCATABLE :: RRKT(:,:)
     REAL(DP),ALLOCATABLE :: RRKS(:,:)
     
     ALLOCATE(RKU(MElemGL,1))              ; RKU  = 0.
     ALLOCATE(RKV(MElemGL,1))              ; RKV  = 0.
     ALLOCATE(RKEL(NNodeGL))               ; RKEL = 0.
     ALLOCATE(RRKT(NNodeGL,1))             ; RRKT = 0.
     ALLOCATE(RRKS(NNodeGL,1))             ; RRKS = 0.

     FNAME = TRIM(OUTDIR)//'/rrktemp/ref2.cdf'
     
     STATUS = nf_open(FNAME,NF_NOWRITE,NCINP)

     IF (STATUS /=0) THEN
       WRITE(IPT,*) 'Could not open cdf file:',FNAME
       STOP
     ENDIF
     
     IDUMMY = 0

! STORE THE ELEVATION IN STTRUE1

     IF(EL_ASSIM) THEN 
       START(1)=1
       START(2)=INDEX
       COUNT(1)=NNodeGL
       COUNT(2)=1
       STATUS = nf_inq_varid(NCINP,'el',VARID)
       RCODE = nf_get_vara_double(NCINP,VARID,START,COUNT,RKEL)
     
       DO J=1,NNodeGL
         IDUMMY = IDUMMY + 1
         STTRUE1(IDUMMY)=RKEL(J)
       ENDDO
       
     ENDIF

!     WRITE(IPT,*)'el', 1,RKEL(1)

! STORE THE SPEED U(MElemGL,KBM1) IN STTRUE1

     IF(UV_ASSIM) THEN
       START(1)=1
       START(2)=1
       START(3)=INDEX
       COUNT(1)=MElemGL
       COUNT(2)=KBM1
       COUNT(3)=1
       STATUS = nf_inq_varid(NCINP,'u',VARID)
       RCODE = nf_get_vara_double(NCINP,VARID,START,COUNT,RKU)

       DO K=1,KBM1
         DO J=1,MElemGL
          IDUMMY = IDUMMY + 1
          STTRUE1(IDUMMY) = RKU(J,K)
         ENDDO
       ENDDO

! STORE THE SPEED V(MElemGL,KBM1) IN STTRUE1

     STATUS = nf_inq_varid(NCINP,'v',VARID)
     RCODE = nf_get_vara_double(NCINP,VARID,START,COUNT,RKV)

       DO K=1,KBM1
         DO J=1,MElemGL
           IDUMMY = IDUMMY + 1
           STTRUE1(IDUMMY) = RKV(J,K)
         ENDDO
       ENDDO
     ENDIF
     
     IF(T_ASSIM) THEN
        START(1) = 1
	START(2) = 1
        START(3) = INDEX
        COUNT(1) = NNodeGL
        COUNT(2) = KBM1
        COUNT(3) = 1
        STATUS = nf_inq_varid(NCINP,'temp',VINPID)
        STATUS = nf_get_vara_double(NCINP,VINPID,START,COUNT,RRKT)
        DO K = 1, KBM1
          DO J = 1, NNodeGL
	     IDUMMY = IDUMMY + 1
	     STTRUE1(IDUMMY) = RRKT(J,K)
	  ENDDO
        ENDDO	
      ENDIF

      IF(S_ASSIM) THEN
        START(1) = 1
	START(2) = 1
        START(3) = INDEX
        COUNT(1) = NNodeGL
        COUNT(2) = KBM1
        COUNT(3) = 1
        STATUS = nf_inq_varid(NCINP,'sal',VINPID)
        STATUS = nf_get_vara_double(NCINP,VINPID,START,COUNT,RRKS)
	
	DO K = 1, KBM1
          DO J = 1, NNodeGL
	     IDUMMY = IDUMMY + 1
	     STTRUE1(IDUMMY) = RRKS(J,K)
	  ENDDO
        ENDDO	
      ENDIF

     RCODE = nf_close(NCINP)

     DEALLOCATE(RKU,RKV,RKEL,RRKT,RRKS)
   
   RETURN
   END SUBROUTINE READREF2
   

   
   SUBROUTINE READMEDM2
   
   USE ALL_VARS
   USE MOD_RRK
   IMPLICIT NONE

   INTEGER NUMS
   integer i,k,idump,IINTT,NGL9,MGL9
   REAL(SP) THOUR9
   CHARACTER DIR*120,FILE_NO*4,FILENAME*120
   REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::UGL,VGL,WGL,KMGL
   REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::S1GL,T1GL,RHO1GL
   REAL(SP), ALLOCATABLE, DIMENSION(:)     ::ELGL,UAGL,VAGL

   ALLOCATE(UGL(MElemGL,KBM1))     
   ALLOCATE(VGL(MElemGL,KBM1))   
   ALLOCATE(WGL(MElemGL,KBM1))     
   ALLOCATE(KMGL(MElemGL,KBM1))     
   ALLOCATE(ELGL(NNodeGL))
   ALLOCATE(UAGL(0:MElemGL))
   ALLOCATE(VAGL(0:MElemGL))
   ALLOCATE(T1GL(NNodeGL,KBM1))      
   ALLOCATE(S1GL(NNodeGL,KBM1))      
   ALLOCATE(RHO1GL(NNodeGL,KBM1))  
   if(MOD(IEND,DELTA_ASS) .ne. 0) then
      print*,'-------Error in read ture obs data'
      print*,'iint,int(iint/DELTA_ASS)*DELTA_ASS=',iint,int(iint/DELTA_ASS)*DELTA_ASS
      stop
   endif
   
!=== read true state of current time step to sttr0
   
   WRITE(FILE_NO,'(I4.4)') ISTART/DELTA_ASS
   DIR = '../medm_bck/'
   FILENAME='chn_sim'//FILE_NO//'.dat'
   OPEN(1,file=TRIM(DIR)//TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED')

   READ(1) IINTT,NGL9,MGL9,THOUR9  

     DO I=1,MElemGL
        READ(1) (UGL(I,K),VGL(I,K),WGL(I,K),KMGL(I,K), K=1,KBM1)
     ENDDO
     DO I=1,NNodeGL
        READ(1) ELGL(I),(T1GL(I,K),S1GL(I,K),RHO1GL(I,K),K=1,KBM1)
     ENDDO     

     IDUMP=0
       IF(EL_ASSIM) THEN
       DO I=1,NNodeGL
          IDUMP=IDUMP+1
          STTR0(IDUMP)=ELGL(I)
       ENDDO
     ENDIF
     
     IF(UV_ASSIM) THEN
       DO K=1, KBM1
	 DO I=1, MElemGL
           IDUMP = IDUMP + 1
           STTR0(IDUMP) = UGL(I,K) 
	 ENDDO
       ENDDO
       DO K=1, KBM1
	 DO I=1, MElemGL
           IDUMP = IDUMP + 1
           STTR0(IDUMP) = VGL(I,K)
	 ENDDO
       ENDDO
     ENDIF
     
     IF(T_ASSIM) THEN         ! NEED READ T ABOVE
     DO K=1, KBM1
       DO I=1, NNodeGL
         IDUMP = IDUMP + 1
         STTR0(IDUMP) = T1GL(I,K) 
       ENDDO
     ENDDO
     ENDIF
     
     IF(S_ASSIM) THEN
       DO K=1, KBM1
	 DO I=1, NNodeGL
           IDUMP = IDUMP + 1
           STTR0(IDUMP) = S1GL(I,K) 
	 ENDDO
       ENDDO
     ENDIF

!=== read true state of next time step to sttr1
     
   WRITE(FILE_NO,'(I4.4)') ISTART/DELTA_ASS+1
   DIR = '../medm_bck/'
   FILENAME='chn_sim'//FILE_NO//'.dat'
   OPEN(1,file=TRIM(DIR)//TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED')

   READ(1) IINTT,NGL9,MGL9,THOUR9  

     DO I=1,MElemGL
        READ(1) (UGL(I,K),VGL(I,K),WGL(I,K),KMGL(I,K), K=1,KBM1)
     ENDDO
     DO I=1,NNodeGL
        READ(1) ELGL(I),(T1GL(I,K),S1GL(I,K),RHO1GL(I,K),K=1,kbm1)
     ENDDO    

     IDUMP=0
       IF(EL_ASSIM) THEN
       DO I=1,NNodeGL
          IDUMP=IDUMP+1
          STTR1(IDUMP)=ELGL(I)
       ENDDO
     ENDIF
     
     IF(UV_ASSIM) THEN
       DO K=1, KBM1
	 DO I=1, MElemGL
           IDUMP = IDUMP + 1
           STTR1(IDUMP) = UGL(I,K) 
	 ENDDO
       ENDDO
       DO K=1, KBM1
	 DO I=1, MElemGL
           IDUMP = IDUMP + 1
           STTR1(IDUMP) = VGL(I,K)
	 ENDDO
       ENDDO
     ENDIF
     
     IF(T_ASSIM) THEN         ! NEED READ T ABOVE
     DO K=1, KBM1
       DO I=1, NNodeGL
         IDUMP = IDUMP + 1
         STTR1(IDUMP) = T1GL(I,K) 
       ENDDO
     ENDDO
     ENDIF
     
     IF(S_ASSIM) THEN
       DO K=1, KBM1
	 DO I=1, NNodeGL
           IDUMP = IDUMP + 1
           STTR1(IDUMP) = S1GL(I,K) 
	 ENDDO
       ENDDO
     ENDIF
     
     DEALLOCATE(UGL,VGL,WGL,KMGL)
     DEALLOCATE(ELGL,UAGL,VAGL)
     DEALLOCATE(S1GL,T1GL,RHO1GL)
     CLOSE(1)
     RETURN
     
     END SUBROUTINE READMEDM2
     
# endif
END MODULE MOD_RRKA
