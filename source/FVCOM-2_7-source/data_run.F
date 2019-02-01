!==============================================================================|
!   Input Parameters Which Control the Model Run                               |
!==============================================================================|

   SUBROUTINE DATA_RUN            

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_INP
!C. Lee, 12.21.2006, modified for restart T.W.?? for DTI_Ratio
   USE CONTROL
!C. Lee end T.W.

#if defined (TS_OBC)
   USE MOD_TSOBC, only : TSOBC_ON,ALPHA_SERIES_OBC ! B Clark Aug 2015 for obc nudging clarification
#endif
   
# if defined (VISIT)
   USE MOD_VISIT, only : VISIT_OPT
# endif
   IMPLICIT NONE
   REAL(SP) REALVEC(150)
   INTEGER  INTVEC(150),ISCAN,KTEMP
   REAL(SP) ZKUTMP,ZKLTMP
   CHARACTER(LEN=120) :: FNAME
   INTEGER I


!==============================================================================|
!   READ IN VARIABLES AND SET VALUES                                           |
!==============================================================================|

   FNAME = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------|
!     "INFO FILE"   !!
!------------------------------------------------------------------------------|

! David Changed default IPT to 6 to print error reading infofile to screen
   IPT =6
   ISCAN = SCAN_FILE(FNAME,"INFOFILE",CVAL = INFOFILE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING INFOFILE: ',ISCAN
     CALL PSTOP
   END IF

!
!-----------------OPEN RUNTIME INFO FILE---------------------------------------!
!
   IF(TRIM(INFOFILE) /= "screen")THEN
      IPT = 71
      CALL FOPEN(IPT, TRIM(INFOFILE) ,"ofr")
!   ELSE
!     IPT = 6
   END IF

!
!-----------------WRITE BANNER ------------------------------------------------!
!
   IF(NPROCS > 1 .AND. MSR)THEN
   CALL WRITE_BANNER(IPT)
   END IF


!------------------------------------------------------------------------------|
!   EXTERNAL TIME STEP (DTE) 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"DTE",FSCAL = DTE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DTE: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "ISPLIT"   -   RATIO OF EXTERNAL/INTERNAL MODE STEPS
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ISPLIT",ISCAL = ISPLIT)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ISPLIT: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "IRAMP"   -NUMBER OF INTEGRATIONS OVER WHICH TO RAMP UP MODEL
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IRAMP",ISCAL = IRAMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IRAMP: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "NSTEPS"   -NUMBER OF TIME STEPS TO RUN MODEL                  
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"NSTEPS",ISCAL = NSTEPS)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING NSTEPS: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "IRHO_MEAN"   -NUMBER OF INTEGRATIONS OVER WHICH TO call rho_mean
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IRHO_MEAN",ISCAL = IRHO_MEAN)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IRHO_MEAN: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "IRECORD"   -CONTROLS RECORD PRINTING                                  
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IRECORD",ISCAL = IRECORD)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IRECORD: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "IDMPSMS"   -CONTROLS SMS FILE DUMPING
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IDMPSMS",ISCAL = IDMPSMS)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IDMPSMS: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "RESTART -CONTROLS RESTART TYPE 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"RESTART",CVAL = RESTART)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RESTART: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "BFRIC"   -CONTROLS BOTTOM FRICTION
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BFRIC",FSCAL = BFRIC)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BFRIC: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "BROUGH_TYPE"   !!CONTROLS BOTTOM ROUGHNESS CALCULATION 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BROUGH_TYPE",CVAL = BROUGH_TYPE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BROUGH_TYPE: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "Z0B"   -CONTROLS BOTTOM FRICTION
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"Z0B",FSCAL = Z0B)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING Z0B: ',ISCAN
     CALL PSTOP 
   END IF

!------------------------------------------------------------------------------|
!     "HORZMIX"   -HORIZONTAL DIFFUSION TYPE
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"HORZMIX",CVAL = HORZMIX)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING HORZMIX: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "HORCON"   !!HORIZONTAL DIFFUSION COEFFICIENT
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"HORCON",FSCAL = HORCON)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING HORCON: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "HPRNU"   !!HORIZONTAL DIFFUSION COEFFICIENT
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"HPRNU",FSCAL = HPRNU)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING HPRNU: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "VERTMIX"   !!VERTICAL DIFFUSION TYPE
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"VERTMIX",CVAL = VERTMIX)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING VERTMIX: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "UMOL"   !!VERTICAL DIFFUSION COEFFICIENT
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"UMOL",FSCAL = UMOL)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING UMOL: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "VPRNU"   !!VERTICAL DIFFUSION COEFFICIENT
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"VPRNU",FSCAL = VPRNU)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING VPRNU: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "C_BAROPG"   !!CONTROLS BARO PRESSURE GRADIENT CALC
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"C_BAROPG",CVAL = C_BAROPG)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING C_BAROPG: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "CTRL_DEN"   !!CONTROLS BARO PRESSURE GRADIENT CALC
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"CTRL_DEN",CVAL = CTRL_DEN)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING CTRL_DEN: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "H_TYPE"   !!CONTROLS BARO PRESSURE GRADIENT CALC
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"H_TYPE",CVAL = H_TYPE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING H_TYPE: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "S_TYPE"   !!CONTROLS TIDAL FORCING 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"S_TYPE",CVAL = S_TYPE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING S_TYPE: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "M_TYPE"   !!CONTROLS METEROLOGICAL FORCING        
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"M_TYPE",CVAL = M_TYPE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING M_TYPE: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "WINDTYPE"   -CONTROLS BARO PRESSURE GRADIENT CALC
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"WINDTYPE",CVAL = WINDTYPE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING WINDTYPE: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "DJUST"   !!DEPTH ADJUSTMENT FACTOR        
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DJUST",FSCAL = DJUST)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DJUST: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "ZETA1"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ZETA1",FSCAL = ZETA1) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ZETA1: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "ZETA2"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ZETA2",FSCAL = ZETA2)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ZETA2: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "RHEAT"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"RHEAT",FSCAL = RHEAT)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RHEAT: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "THOUR_HS"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"THOUR_HS",FSCAL = THOUR_HS)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING THOUR_HS: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "MIN_DEPTH"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"MIN_DEPTH",FSCAL = MIN_DEPTH)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING MIN_DEPTH: ',ISCAN
     CALL PSTOP
   END IF
   
!------------------------------------------------------------------------------|
!   "AVGE_ON" TURNS ON OUTPUT OF FLOW FIELD AVERAGES 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"AVGE_ON",LVAL = AVGE_ON)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING AVGE_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "INT_AVGE"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"INT_AVGE",ISCAL = INT_AVGE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING INT_AVGE: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "NUM_AVGE"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"NUM_AVGE",ISCAL = NUM_AVGE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING NUM_AVGE: ',ISCAN
     CALL PSTOP
   END IF
   IF(INT_AVGE <= 0 .OR. NUM_AVGE <=0) AVGE_ON = .FALSE.
!------------------------------------------------------------------------------|
!     "BEG_AVGE"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BEG_AVGE",ISCAL = BEG_AVGE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BEG_AVGE: ',ISCAN
     CALL PSTOP
   END IF


!------------------------------------------------------------------------------|
!     "IRESTART"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IRESTART",ISCAL = IRESTART)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IRESTART: ',ISCAN
     CALL PSTOP 
   END IF


!C. Lee, 12.21.2006, added for restart, to control restarting time in main input file. T.W.
!------------------------------------------------------------------------------|
!     "DTI_RATIO"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DTI_RATIO",FSCAL = DTI_RATIO)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DTI_RATIO: ',ISCAN
     CALL PSTOP
   END IF
!C. Lee end, T.W.



!------------------------------------------------------------------------------|
!     "IREPORT"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IREPORT",ISCAL = IREPORT)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IREPORT: ',ISCAN
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!     "KSL"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"KSL",ISCAL = KSL)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING KSL: ',ISCAN
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!     "DPTHSL"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DPTHSL",FVEC = REALVEC,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DPTHSL: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= KSL)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED DEPTHS IN DEPTHSL IS NOT EQUAL TO KSL' 
       WRITE(*,*)'KSL: ',KSL
       WRITE(*,*)'DPTHSL: ',REALVEC       !DPTHSL
     END IF
   END IF
  
   ALLOCATE(DPTHSL(KSL))
   DPTHSL(1:KSL)= REALVEC(1:KSL)
!------------------------------------------------------------------------------|
!     "INDEX_VERCOR"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"INDEX_VERCOR",ISCAL = INDEX_VERCOR)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING INDEX_VERCOR: ',ISCAN
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!     "P_SIGMA"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"P_SIGMA",FSCAL = P_SIGMA)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING P_SIGMA: ',ISCAN
     CALL PSTOP
   END IF

! Added by B Clark Nov 2018 for setting minimum depth for hybrid coordinates
!  
!		HMIN2
!-------------------------------------------------   
    ISCAN = SCAN_FILE(FNAME,"HMIN2",FSCAL = HMIN2)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING HMIN2: ',ISCAN
     CALL PSTOP
   END IF    
!------------------------------------------------------------------------------|
!     "DU2"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DU2",FSCAL = DU2)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DU2: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "DL2"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DL2",FSCAL = DL2)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DL2: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "DUU"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DUU",FSCAL = DUU)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DUU: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "DLL"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DLL",FSCAL = DLL)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DLL: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "HMIN1"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"HMIN1",FSCAL = HMIN1)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING HMIN1: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "KU"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"KU",ISCAL = KU)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING KU: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "KL"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"KL",ISCAL = KL)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING KL: ',ISCAN
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "ZKU"   !! 
!------------------------------------------------------------------------------|
   IF(KU > 1)THEN
     ISCAN = SCAN_FILE(FNAME,"ZKU",FVEC = REALVEC,NSZE = KTEMP)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING ZKU: ',ISCAN
       CALL PSTOP
     END IF
     IF(MSR)THEN
       IF(KTEMP /= KU)THEN
         WRITE(*,*)'NUMBER OF SPECIFIED DEPTHS IN ZKU IS NOT EQUAL TO KU' 
         WRITE(*,*)'KU: ',KU
         WRITE(*,*)'ZKU: ',REALVEC      
       END IF
     END IF
  
     ALLOCATE(ZKU(KU)); ZKU=0.0_SP
     ZKU(1:KU)= REALVEC(1:KU)
   ELSE
     ISCAN = SCAN_FILE(FNAME,"ZKU",FSCAL = ZKUTMP)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING ZKU: ',ISCAN
       CALL PSTOP
     END IF
  
     ALLOCATE(ZKU(KU)); ZKU=0.0_SP
     IF(KU > 0)ZKU(1)= ZKUTMP
       
   END IF  
!------------------------------------------------------------------------------|
!     "ZKL"   !! 
!------------------------------------------------------------------------------|
   IF(KL > 1)THEN
     ISCAN = SCAN_FILE(FNAME,"ZKL",FVEC = REALVEC,NSZE = KTEMP)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING ZKL: ',ISCAN
       CALL PSTOP
     END IF
     IF(MSR)THEN
       IF(KTEMP /= KL)THEN
         WRITE(*,*)'NUMBER OF SPECIFIED DEPTHS IN ZKL IS NOT EQUAL TO KL' 
         WRITE(*,*)'KL: ',KL
         WRITE(*,*)'ZKL: ',REALVEC      
       END IF
     END IF
  
     ALLOCATE(ZKL(KL)); ZKL=0.0_SP
     ZKL(1:KL)= REALVEC(1:KL)
   ELSE IF(KL > 0)THEN
     ISCAN = SCAN_FILE(FNAME,"ZKL",FSCAL = ZKLTMP)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING ZKL: ',ISCAN
       CALL PSTOP
     END IF
  
     ALLOCATE(ZKL(KL)); ZKL=0.0_SP
     ZKL(1)= ZKLTMP
   ELSE IF(KL == 0)THEN
     ALLOCATE(ZKL(KL+1)); ZKL=0.0_SP
   ELSE 
     WRITE(IPT,*)'KL SHOULD BE >= 0: KL= ',KL
     CALL PSTOP
   END IF  
!------------------------------------------------------------------------------|
!     "KB"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"KB",ISCAL = KB)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING KB: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "DELTT"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DELTT",FSCAL = DELTT)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DELTT: ',ISCAN
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "CASETITLE"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"CASETITLE",CVAL = CASETITLE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING CASETITLE: ',ISCAN
     CALL PSTOP 
   END IF

!------------------------------------------------------------------------------|
!     "INPDIR"   !!DIRECTORY FOR INPUT FILES            
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"INPDIR",CVAL = INPDIR)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING INPDIR: ',ISCAN
     CALL PSTOP 
   END IF
   I = LEN_TRIM(INPDIR)
   IF(INPDIR(I:I) == "/") INPDIR(I:I) = " "

!------------------------------------------------------------------------------|
!     "OUTDIR"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"OUTDIR",CVAL = OUTDIR)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OUTDIR: ',ISCAN
     CALL PSTOP 
   END IF
   I = LEN_TRIM(OUTDIR)
   IF(OUTDIR(I:I) == "/") OUTDIR(I:I) = " "

!------------------------------------------------------------------------------|
!   CONVECTIVE OVERTURNING FLAG
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"VERT_STAB",LVAL = VERT_STAB)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING VERT_STAB: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!  TEMPERATURE/SALINITY AVERAGING FLAG
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"TS_FCT",LVAL = TS_FCT)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING TS_FCT: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "TS_FCT" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
                                                                                                                            
!------------------------------------------------------------------------------|
!  BAROTROPIC FLAG
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"BAROTROPIC",LVAL = BAROTROPIC)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BAROTROPIC: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "BAROTROPIC" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!  TEMP_ON        
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"TEMP_ON",LVAL = TEMP_ON)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING TEMP_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "TEMP_ON" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
                                                                                                                        
!------------------------------------------------------------------------------|
!  SALINITY_ON 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"SALINITY_ON",LVAL = SALINITY_ON)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING SALINITY_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "SALINITY_ON" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
                                                                                                                        
!------------------------------------------------------------------------------|
!  'SURFACEWAVE_MIX'  If T, Surface wave induced by wind   
!                     If F, No wind induced,considering Richardson # dep. 
!                           dissipation correction
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"SURFACEWAVE_MIX",LVAL = SURFACEWAVE_MIX)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING SURFACEWAVE_MIX: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "SURFACEWAVE_MIX" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF

!  GWC FOR TEMP/SALT NUDGING ON OBC
!------------------------------------------------------------------------------|
!  'TS_NUDGING_OBC'   If T, Nudge Temperature and Salinity on OBC
!                     If F, No Nuding on OBC     
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"TS_NUDGING_OBC",LVAL = TS_NUDGING_OBC)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING TS_NUDGING_OBC: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "TS_NUDGING_OBC" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!     "ALPHA_OBC"   !!
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ALPHA_OBC",FSCAL = ALPHA_OBC)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ALPHA_OBC: ',ISCAN
     CALL PSTOP
   END IF
   IF(.not.TS_NUDGING_OBC) ALPHA_OBC = 0.0_SP

!------------------------------------------------------------------------------|
!     TIDE_INITIAL and TIDE_INTERVAL
!------------------------------------------------------------------------------|
#  if defined (TIDE_OUTPUT)
   ISCAN = SCAN_FILE(FNAME,"TIDE_INITIAL",ISCAL = TIDE_INITIAL)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING TIDE_INITIAL: ',ISCAN
     CALL PSTOP 
   END IF
    
   ISCAN = SCAN_FILE(FNAME,"TIDE_INTERVAL",ISCAL = TIDE_INTERVAL)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING TIDE_INTERVAL: ',ISCAN
     CALL PSTOP 
   END IF                                                                                                                    
# endif
! SEDIMENT BEGIN

!------------------------------------------------------------------------------|
!  'SEDIMENT_ON'      If T, Use Sediment Model
!------------------------------------------------------------------------------|
   SEDIMENT_ON = .false.
   RESTART_SED = .false.
#  if defined(SEDIMENT)
   ISCAN = SCAN_FILE(TRIM(FNAME),"SEDIMENT_ON",LVAL = SEDIMENT_ON)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING SEDIMENT_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "SEDIMENT_ON" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
!------------------------------------------------------------------------------|
!  'RESTART_SED'      If T and Restart Case, Restart Sediment also
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RESTART_SED",LVAL = RESTART_SED)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RESTART_SED: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "RESTART_SED" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
   IF(TRIM(RESTART) == 'cold_start') RESTART_SED = .false.
#  endif 
!------------------------------------------------------------------------------|
! SEDIMENT END

!------------------------------------------------------------------------------|
!   "ADCOR_ON" TURNS ON OUTPUT OF SEMI-IMPLICIT CORIOLIS TERM
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"ADCOR_ON",LVAL = ADCOR_ON)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING AVGE_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   END IF
   
   
! "WIND_BLOCK logical TRUE if want to specify nodes that have 0.0 wind velocity on 
 ! uniformly forced Grid   B Clark Nov 2015
   
      ISCAN = SCAN_FILE(TRIM(FNAME),"WIND_BLOCK",LVAL = WIND_BLOCK )
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING WIND_BLOCK: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "WIND_BLOCK" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF


# if defined(VISIT)
!------------------------------------------------------------------------------|
!     "VISIT_OPT"   Controls complexity of Meta Data Returned to Visit 
!------------------------------------------------------------------------------|
   ! NOTE: Visit is set up not to crash the model when ever possible
   !       Missing run file option returns basic instead of pstop.
   ISCAN = SCAN_FILE(FNAME,"VISIT_OPT",CVAL = VISIT_OPT )
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING VISIT_OPT: ',ISCAN

     WRITE(IPT,*) 'VISIT_OPT run file parameter missing:'
     WRITE(IPT,*) '++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(IPT,*) 'SETTING VISIT_OPT TO BASIC'
     WRITE(IPT,*) '++++++++++++++++++++++++++++++++++++++++++++'

     VISIT_OPT = 'basic'
   END IF
# endif


!==============================================================================|
!            SET PHYSICAL PARAMETERS                                           !
!==============================================================================|

   KBM1 = KB-1 ; KBM2 = KB-2 
   DTI=DTE*FLOAT(ISPLIT)
   IEND=NSTEPS
   DAYS=NSTEPS*DTI/24.0_SP/3600.0_SP
   IINT = 0
   IF(IREPORT == 0) IREPORT = IEND+2

!==============================================================================|
!            ERROR CHECKING                                                    !
!==============================================================================|


   IF(RESTART /= 'cold_start' .AND. RESTART /= 'hot_cold_s' .AND. &
         RESTART /= 'hot_start')THEN
     IF(MSR)WRITE(IPT,*) 'RESTART NOT CORRECT --->',RESTART   
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "cold_start","hot_cold_s", or "hot_start"'
     STOP
   END IF
   IF(HORZMIX /= 'constant' .AND. HORZMIX /= 'closure') THEN
     IF(MSR)WRITE(IPT,*) 'HORZMIX NOT CORRECT --->',HORZMIX  
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "constant" or "closure"'
     STOP
   END IF
   IF(C_BAROPG /= 'sigma' .AND. C_BAROPG /= 's_levels') THEN
     IF(MSR)WRITE(IPT,*) 'C_BAROPG NOT CORRECT --->',C_BAROPG
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "sigma" or "s_levels"'
     STOP
   END IF
   IF(CTRL_DEN /= 'sigma-t' .AND. CTRL_DEN /= 'pdensity' .AND. &
      CTRL_DEN /= 'sigma-t_stp') THEN
     IF(MSR)WRITE(IPT,*) 'CTRL_DEN NOT CORRECT,--->',CTRL_DEN
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "sigma-t" , "pdensity", or "sigma-t_stp"'
     STOP
   END IF
   IF(WINDTYPE /= 'stress' .AND. WINDTYPE /= 'speed') THEN
     IF(MSR)WRITE(IPT,*) 'WINDTYPE NOT CORRECT,--->',WINDTYPE
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "stress" or "speed"'
     STOP 
   END IF
   IF(M_TYPE /= 'uniform' .AND. M_TYPE /= 'non-uniform') THEN
     IF(MSR)WRITE(IPT,*) 'M_TYPE NOT CORRECT,--->',M_TYPE
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "uniform" or "non-uniform"'
     STOP 
   END IF
   IF(S_TYPE /= 'julian' .AND. S_TYPE /= 'non-julian') THEN
     IF(MSR)WRITE(IPT,*) 'S_TYPE NOT CORRECT,--->',S_TYPE
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "julian" or "non-julian"'
     STOP 
   END IF
   IF(BROUGH_TYPE /= 'orig' .AND. BROUGH_TYPE /= 'gotm' .AND. &
      BROUGH_TYPE /= 'user_defined') THEN
     IF(MSR)WRITE(IPT,*) 'BROUGH_TYPE NOT CORRECT,--->',BROUGH_TYPE
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "orig" or "gotm" or "user_defined"'
     STOP 
   END IF
   IF(H_TYPE /= 'body_h' .AND. H_TYPE /= 'flux_h') THEN
     IF(MSR)WRITE(IPT,*) 'H_TYPE NOT CORRECT,--->',H_TYPE
     IF(MSR)WRITE(IPT,*) 'SHOULD BE "body_h" or "flux_h"'
     STOP 
   END IF

   IF(KB > 200)THEN
     WRITE(IPT,*)'KB EXCEEDS 200'
     WRITE(IPT,*)'THIS WILL CAUSE ERROR IN SEVERAL READ STATEMENTS SINCE'
     WRITE(IPT,*)'ASSOCIATED FORMAT STATEMENT ASSUMES MAX KB OF 200'
     WRITE(IPT,*)'GREP CODE FOR READ AND 200 TO SEE'
     call PSTOP
   END IF

# if defined (VISIT)
   IF (TRIM(VISIT_OPT) .NE. 'basic' .AND. TRIM(VISIT_OPT) .NE. 'advanced') then
     IF(MSR)WRITE(IPT,*) 'VISIT_OPT returned undefined option: ',TRIM(VISIT_OPT)
     IF(MSR)WRITE(IPT,*) 'Know options for VISIT_OPT are: advanced  OR  basic'
     IF(MSR)WRITE(IPT,*) '++++++++++++++++++++++++++++++++++++++++++++'
     IF(MSR)WRITE(IPT,*) 'SETTING VISIT_OPT TO BASIC'
     IF(MSR)WRITE(IPT,*) '++++++++++++++++++++++++++++++++++++++++++++'
     VISIT_OPT = 'basic'
  END IF
# endif 

!==============================================================================|
!            REPORTING                                                         !
!==============================================================================|
   IF(MSR)WRITE(IPT,*)'!  # STD SALINITY LEVELS :',KSL
   IF(MSR)WRITE(IPT,*)'!  # OF SIGMA LEVELS     :',KB

!==============================================================================|
!            SCREEN REPORT OF SET VARIABlES                                    !
!==============================================================================|
   IF(MSR)THEN
   WRITE(IPT,*)'!  # DTE                 :',DTE
   WRITE(IPT,*)'!  # ISPLIT              :',ISPLIT
   WRITE(IPT,*)'!  # IRAMP               :',IRAMP
   WRITE(IPT,*)'!  # NSTEPS              :',NSTEPS
   WRITE(IPT,*)'!  # IRHO_MEAN           :',IRHO_MEAN
   WRITE(IPT,*)'!  # RESTART             :',TRIM(RESTART)
   WRITE(IPT,*)'!  # BFRIC               :',BFRIC
   WRITE(IPT,*)'!  # MIN_DEPTH           :',MIN_DEPTH
   WRITE(IPT,*)'!  # Z0B                 :',Z0B   
   WRITE(IPT,*)'!  # HORZMIX             :',TRIM(HORZMIX)
   WRITE(IPT,*)'!  # HORCON              :',HORCON
   WRITE(IPT,*)'!  # HPRNU               :',HPRNU
   WRITE(IPT,*)'!  # VERTMIX             :',TRIM(VERTMIX)
   WRITE(IPT,*)'!  # UMOL                :',UMOL
   WRITE(IPT,*)'!  # VPRNU               :',VPRNU
   WRITE(IPT,*)'!  # C_BAROPG            :',TRIM(C_BAROPG)
   WRITE(IPT,*)'!  # CTRL_DEN            :',TRIM(CTRL_DEN)
   WRITE(IPT,*)'!  # H_TYPE              :',TRIM(H_TYPE   )
   WRITE(IPT,*)'!  # WINDTYPE            :',TRIM(WINDTYPE   )
   WRITE(IPT,*)'!  # DJUST               :',DJUST   
   WRITE(IPT,*)'!  # ZETA1               :',ZETA1   
   WRITE(IPT,*)'!  # ZETA2               :',ZETA2   
   WRITE(IPT,*)'!  # RHEAT               :',RHEAT   
   WRITE(IPT,*)'!  # THOUR_HS            :',THOUR_HS           
   WRITE(IPT,*)'!  # IRESTART            :',IRESTART   
!C. Lee,12.21.2006, modified for restart  T.W.
   WRITE(IPT,*)'!  # DTI_RATIO           :',DTI_RATIO
!C. Lee end T.W.
   WRITE(IPT,*)'!  # IREPORT             :',IREPORT    
   WRITE(IPT,*)'!  # IRECORD             :',IRECORD    
   WRITE(IPT,*)'!  # IDMPSMS             :',IDMPSMS    
   WRITE(IPT,*)'!  # KSL                 :',KSL      
   WRITE(IPT,*)'!  # DPTHSL              :',DPTHSL   
   WRITE(IPT,*)'!  # P_SIGMA             :',P_SIGMA   
   WRITE(IPT,*)'!  # KB                  :',KB   
   WRITE(IPT,*)'!  # DELTT               :',DELTT   
   WRITE(IPT,*)'!  # CASETITLE           :',TRIM(CASETITLE)
   WRITE(IPT,*)'!  # M_TYPE              :',TRIM(M_TYPE   )
   WRITE(IPT,*)'!  # S_TYPE              :',TRIM(S_TYPE   )
   WRITE(IPT,*)'!  # BROUGH_TYPE         :',TRIM(BROUGH_TYPE)
   WRITE(IPT,*)'!  # OUTDIR              :',TRIM(OUTDIR   )
   WRITE(IPT,*)'!  # INPDIR              :',TRIM(INPDIR   )
   WRITE(IPT,*)'!  # INFOFILE            :',TRIM(INFOFILE )

   IF(AVGE_ON)THEN
     WRITE(IPT,*)'!  # FLOW AVGES OUTPUT   :  ACTIVE'
     WRITE(IPT,*)'!  # START ITERATION     :',BEG_AVGE
     WRITE(IPT,*)'!  # AVGING INTERVAL     :',INT_AVGE
     WRITE(IPT,*)'!  # NUMBER OF INTERVALS :',NUM_AVGE
   ELSE
     WRITE(IPT,*)'!  # FLOW AVGES OUTPUT   :  INACTIVE'
   END IF
   IF(VERT_STAB)THEN
     WRITE(IPT,*)'!  # CONVECTIVE OVERTURN :  ACTIVE'
   ELSE
     WRITE(IPT,*)'!  # CONVECTIVE OVERTURN :  INACTIVE'
   END IF
   IF(TS_FCT)THEN
     WRITE(IPT,*)'!  # TEMP/SAL AVERAGING  :  ACTIVE'
   ELSE
     WRITE(IPT,*)'!  # TEMP/SAL AVERAGING  :  INACTIVE'
   END IF
   IF(BAROTROPIC)THEN
     WRITE(IPT,*)'!  # BAROTROPIC RUN      :  ACTIVE'
   END IF
   IF(TEMP_ON .AND. .NOT. BAROTROPIC)THEN
     WRITE(IPT,*)'!  # TEMPERATURE EQUATION:  ACTIVE'
   END IF
   IF(SALINITY_ON .AND. .NOT. BAROTROPIC)THEN
     WRITE(IPT,*)'!  # SALINITY EQUATION   :  ACTIVE'
   END IF
   IF(SURFACEWAVE_MIX)THEN
     WRITE(IPT,*)'!  # SURFACE WAVE MIXING :  ACTIVE'
   END IF
#if defined(TS_OBC) 
  IF(TSOBC_ON) THEN     ! Added by B Clark for TSOBC nudging clarification on run time.  
   write(ipt,*)'! # OBC TS SERIES NUDGING  : ACTIVE '
   write(ipt,*)'! # NUDGING COEFF          :',ALPHA_SERIES_OBC
  END IF
# else 
   IF(TS_NUDGING_OBC)THEN 
     WRITE(IPT,*)'!  # OBC TS NUDGING      :  ACTIVE'
     WRITE(IPT,*)'!  # NUDGING  COEFF      :',ALPHA_OBC
   END IF
#endif  
!  SEDIMENT BEGIN
   IF(SEDIMENT_ON)THEN
     WRITE(IPT,*)'!  # SEDIMENT MODEL      :  ACTIVE'
   END IF
   IF(RESTART_SED)THEN
     WRITE(IPT,*)'!  # SEDIMENT MODEL      :  HOT_STARTED'
  ELSE
     WRITE(IPT,*)'!  # SEDIMENT MODEL      :  COLD_STARTED'
  END IF
  !  END SEDIMENT
  
# if defined (VISIT)
  IF(VISIT_OPT .EQ. "advanced") then 
     WRITE(IPT,*)'!  # VISIT META DATA     :  ADVANCED'
  ELSEIF(VISIT_OPT .EQ. "basic") then 
     WRITE(IPT,*)'!  # VISIT META DATA     :  BASIC'
  END IF
# endif
  
  
  IF(WIND_BLOCK)THEN
  write(IPT,*)'!   # WIND BLOCK IS ON'
  END IF

END IF


!==============================================================================|
!            FORMATS                                                           |
!==============================================================================|
 101  FORMAT(A10," = ",F10.6)
 102  FORMAT(A10," = ",I10)
 103  FORMAT(A10," = ",A25)
1000  FORMAT (80a1)
4000  FORMAT (3i10,1x,a10)
5000  FORMAT (a10,2e10.3)
6000  FORMAT (3(2x,a8),4x,a6)

   RETURN
   END SUBROUTINE DATA_RUN    
!------------------------------------------------------------------------------|
