MODULE MOD_ENKF 
# if defined (ENKF_ASSIM)
   USE CONTROL
   IMPLICIT NONE
   SAVE

   CHARACTER(LEN=80)  ENKF_INIT   !!INITIAL PERTUBATION FIELD OPTION(usr/default)
   CHARACTER(LEN=80)  ENKF_TYPE   !!ENSEMBLE UPDATE OPTION(default/square root)

   INTEGER      ENKF_NENS         !!GLOBAL NUMBER OF ENSEMBLES
   INTEGER      DELTA_ASS         !!ASSIMILATON TIME INTERVAL IN SECONDS
   INTEGER      ENKF_INT          !!ASSIMILATION TIME INTERVAL/FILE OUTPUT INTERVAL (>=1) 
   INTEGER      ENKF_NOBSMAX      !!MAXIMUM NUMBER OF THE OBSERVATION STATIONS  
   INTEGER      ENKF_START        !!ASSIMILATION START TIME      
   INTEGER      ENKF_END          !!ASSIMILATION END TIME
   REAL(SP) ::  ENKF_CINF         !!MAX DISTANCE OF CORRELATIN  
   REAL(SP) ::  OBSERR_EL         !!EL OBSERVATION ERROR SPECIFIED   
   REAL(SP) ::  OBSERR_UV         !!UV OBSERVATION ERROR SPECIFIED  
   REAL(SP) ::  OBSERR_T          !!TEMPERATURE OBSERVATION ERROR SPECIFIED
   REAL(SP) ::  OBSERR_S          !!SALINITY OBSERVATION ERROR SPECIFIED

   LOGICAL  ::  EL_ASSIM          !!OPTION FOR CHOSING ELEVATION AS ASSIMILATION VARIABLES
   LOGICAL  ::  UV_ASSIM          !!OPTION FOR CHOSING CURRENT AS ASSIMILATION VARIABLES  
   LOGICAL  ::  T_ASSIM           !!OPTION FOR CHOSING TEMPERATURE AS ASSIMILATION VARIABLES
   LOGICAL  ::  S_ASSIM           !!OPTION FOR CHOSING SALINITY AS ASSIMILATION VARIABLES

   LOGICAL  ::  EL_OBS            !!OPTION FOR ELEVATION OBSERVATION DATA
   LOGICAL  ::  UV_OBS            !!OPTION FOR CURRENT OBSERVATION DATA  
   LOGICAL  ::  T_OBS             !!OPTION FOR TEMPERATURE OBSERVATION DATA
   LOGICAL  ::  S_OBS             !!OPTION FOR SALINITY OBSERVATION DATA

   INTEGER  ::  INOOB             !! 72 RESERVE FOR I/O INPUT OF OBSERVATION FILE
   INTEGER  ::  INOKF             !! 73 FILE I/O PIPE NUMBER
   INTEGER  ::  IOBCKF            !! 76 I/O PIPE ONLY FOR B.C. TREATMENT IN ENKF 

!-VARIABLES for dimension                 
   INTEGER  ::  NLOC              !!Number of observation locations
   INTEGER  ::  IENS              !!ensemble number index
   INTEGER  ::  NCYC              !!Number of cycles for assimilation
   INTEGER  ::  N1CYC
   INTEGER  ::  ICYC = 0          !!cycle number index
   CHARACTER(LEN=4)   FCYC 

   INTEGER I_INITIAL
   INTEGER,ALLOCATABLE  :: STLOC(:)     !counting number in the state vector of observation location.  
   REAL(DP),ALLOCATABLE :: STSD(:)      !for normalization use 
   REAL(DP),ALLOCATABLE :: WKTMP(:)     !TEMP ARRAY FOR WK() 
   REAL(SP)  BC_AMP_ERR(6)
   REAL(SP)  BC_PHA_ERR(6)

   REAL(DP)  :: INFLREF = 0.0_DP   ! 10.0d9

!-VARIABLES OF enkf     
   REAL(DP),PARAMETER :: ZEROD = 0.0_DP
  
   REAL(DP),ALLOCATABLE  :: WK(:,:)     !!observation error covariance matrix
                                        ! Nobs*Nobs
   REAL(DP),ALLOCATABLE  :: STFCT(:)    !!state vector of one ensemble forecast
                                        !!stfct(stdim)
   REAL(DP),ALLOCATABLE  :: AKF(:,:)    !!state vector of ensemble forecasts,
                                        !!and the difference whith it's mean
                                        !!Akf(stdim,Nens)
                                        !!Akf = Akf -Mkf   
   REAL(DP),ALLOCATABLE  :: MKF(:)      !!mean of state Vector forecast over 
                                        !!Nens ensembles, Mkf(stdim)
   REAL(DP),ALLOCATABLE  :: AFSTD(:)    !!standard difference of state vector Akf
                                        ! Afstd(stdim)
   REAL(DP),ALLOCATABLE  :: SFH(:,:)    !!inflation factor,SfH=H(Xf)/sqrt(R*(Ne-1)) 
                                        ! SfH(Nens,Nobs)
   REAL(DP),ALLOCATABLE  :: SFSF(:,:)   !!SfSf=SfH*SfH'=SfU*SfD*SfU',  SfSf(Nens,Nens)
   REAL(DP),ALLOCATABLE  :: SFD(:)      !!The singular values of of SfSf, SfD(Nens)
   REAL(DP),ALLOCATABLE  :: SFU(:,:)    !!SfU(Nens,Nens)				
   REAL(DP),ALLOCATABLE  :: SK(:,:)     !!H*(Af-Mf),   Sk(1:Nobs,1:Nens)
   REAL(DP),ALLOCATABLE  :: BK(:,:)     !!Bk=Sk*Akf^{T}/(Nens-1)
                                        !   =Sk*(Af-Mf)^{T}/(Nens-1)
                                        !   =H*(Af-Mf)*(Af-Mf)^{T}/(Nens-1)
                                        !   =H*Pf
                                        ! BK=BK*HBHT 
                                        !    the HBHT if function of distance, is 0 at 
                                        !    long distancd
                                        ! BK(1:Nobs,1:stdim)
   REAL(DP),ALLOCATABLE  :: RK(:,:)     !!Rk=H*Pf*H^{T}=H*Bk^{T},  :Rk(1:Nobs,1:Nobs)
   REAL(DP),ALLOCATABLE  :: TK(:,:)     !!Tk = Rk+ Wk , and
                                        ! Tk = Tk^-1,              :Tk(1:Nobs,1:Nobs)
   REAL(DP),ALLOCATABLE  :: DUMMY(:)    !!input for subroutine gaussj,dummy(1:Nobs)
   REAL(DP),ALLOCATABLE  :: KGAIN(:,:)  !!Kalman Gain Matrix,
                                        ! Kgain = Bk'*Tk
                                        !       = Pf*H'*(H*Pf*H'+Wk)^{-1}
                                        !    Kgain(1:stdim, 1:Nobs)
   REAL(DP),ALLOCATABLE  :: STTR(:)     !!sttr: true state vector getting from ref2.cdf
                                        !   sttr(stdim)
   REAL(DP),ALLOCATABLE  :: STTR1(:)
   REAL(DP),ALLOCATABLE  :: STTEMP(:)   !! = Mkf             :sttemp(stdim)
   REAL(DP),ALLOCATABLE  :: OBSDATA(:)  !!true Obs.          :obsdata(Nobs)
   REAL(DP),ALLOCATABLE  :: ERRVEC(:)   !!the difference between the mean of
                                        ! ensemble forecast and the true
                                        ! = H*(Xf-Xt)
                                        ! = Mkf(stloc(i))-sttr(stloc(i))
                                        !  :errvec(1:Nobs)
                                        ! and
                                        ! = Xf-Xt  :errvec(1:stdim)
   REAL(DP),ALLOCATABLE  :: RPERT(:,:)  !!random matrix, which mean equas zero
                                        !Rpert(Nens,Nobs)
   REAL(DP),ALLOCATABLE  :: RPAVE(:)    !! mean of perturbations,   :RpAve(Nobs)
   REAL(DP),ALLOCATABLE  :: MODDATA(:)  !!1.innovation vector y' = H(x_fct)-H(x_obs)
                                        ! 2.model data in observation location 
                                        ! moddata(1:Nloc)
   REAL(DP),ALLOCATABLE  :: STMEAN(:)   !!stmean = Mkf+Kgain*(obsdata(k) - H*Mkf)
                                        ! anlysis vector for mean of ensemble
                                        ! stmean(1:stdim)

   REAL(DP),ALLOCATABLE  :: DIFMD(:)    !!inovation(y(i)-H(xf(i))) for ensemble memeber
                                        ! difmd(Nobs)
   REAL(DP),ALLOCATABLE  :: SA(:,:)
   REAL(DP),ALLOCATABLE  :: SF(:,:)
   REAL(DP),ALLOCATABLE  :: SFH1(:)

   INTEGER  TIMEN
   REAL(DP),ALLOCATABLE  :: EL_SRS(:,:),SRS_TMP(:)
   REAL(DP),ALLOCATABLE  :: TIME_SER(:)
   REAL(SP),ALLOCATABLE  :: EL_INV(:)
   REAL(DP) AMP(6),PHAS(6),PHAI_IJ,FORCE

!-----------------------------------------------------------------------------|  

   CONTAINS !-----------------------------------------------------------------|
            !SET_ENKF_PARAM                                                   | 
            !                                                                 | 
            !                                                                 |
          

   SUBROUTINE SET_ENKF_PARAM
   
   USE CONTROL
   USE MOD_INP
   IMPLICIT NONE
   INTEGER  :: I, ISCAN, KTEMP
   CHARACTER(LEN=120) :: FNAME
   REAL(SP) REALVEC(150)      
   
!  initialize iens
   IENS = 0  

   FNAME = TRIM(INPDIR)//"/"//trim(casename)//"_assim_enkf.dat"
   INQUIRE(FILE=FNAME,EXIST=FEXIST)
   IF(.NOT.FEXIST)THEN
     WRITE(IPT,*)'FILE ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'STOPPING...'
     CALL PSTOP
   END IF
   
!----------------------------------------------------------------------------|
!     "ENKF_INIT"   !! 
!----------------------------------------------------------------------------|     
   ISCAN = SCAN_FILE(TRIM(FNAME),"ENKF_INIT",CVAL = ENKF_INIT)
    IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ENKF_INIT: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ENKF_INIT NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF

!----------------------------------------------------------------------------|
!     "ENKF_NENS"   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ENKF_NENS",ISCAL = ENKF_NENS)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ENKF_NENS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ENKF_NENS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF   
   
!----------------------------------------------------------------------------|
!     "DELTA_ASS "  !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"DELTA_ASS",ISCAL = DELTA_ASS)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DELTA_ASS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE DELTA_ASS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF   
   
!----------------------------------------------------------------------------|
!     "ENKF_INT "   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ENKF_INT",ISCAL = ENKF_INT)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ENKF_INT: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ENKF_INT NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF    
   
!----------------------------------------------------------------------------|
!     "ENKF_NOBSMAX "   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ENKF_NOBSMAX",ISCAL = ENKF_NOBSMAX)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ENKF_NOBSMAX: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ENKF_NOBSMAX NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF       
   
!----------------------------------------------------------------------------|
!     "ENKF_START "   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ENKF_START",ISCAL = ENKF_START)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ENKF_START: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ENKF_START NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF

!----------------------------------------------------------------------------|
!     "ENKF_END "   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ENKF_END",ISCAL = ENKF_END)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ENKF_END: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ENKF_END NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF

!------------------------------------------------------------------------------|
!     "BC_AMP_ERR"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BC_AMP_ERR",FVEC = REALVEC,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BC_AMP_ERR: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= 6)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED TIDAL COMSTITUENTS ERROR ARE NOT EQUAL TO 6' 
     END IF
   END IF
  
   BC_AMP_ERR(1:6)= REALVEC(1:6)
!------------------------------------------------------------------------------|
!     "BC_PHA_ERR"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BC_PHA_ERR",FVEC = REALVEC,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BC_PHA_ERR: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= 6)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED TIDAL COMSTITUENTS ERROR ARE NOT EQUAL TO 6' 
     END IF
   END IF
  
   BC_PHA_ERR(1:6)= REALVEC(1:6)

!----------------------------------------------------------------------------|
!   "ENKF_CINF" MAX LONG DISTANCE OF CORRELATIN 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"ENKF_CINF",FSCAL = ENKF_CINF)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ENKF_CINF: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ENKF_CINF NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   "OBSERR_EL" 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OBSERR_EL",FSCAL = OBSERR_EL)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OBSERR_EL: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE OBSERR_EL NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   "OBSERR_UV" 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OBSERR_UV",FSCAL = OBSERR_UV)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OBSERR_UV: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE OBSERR_UV NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   "OBSERR_T" 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OBSERR_T",FSCAL = OBSERR_T)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OBSERR_T: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE OBSERR_T NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   "OBSERR_S" 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OBSERR_S",FSCAL = OBSERR_S)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OBSERR_S: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE OBSERR_S NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF
   
!----------------------------------------------------------------------------|
!     "ENKF_TYPE"
!----------------------------------------------------------------------------|     
   ISCAN = SCAN_FILE(TRIM(FNAME),"ENKF_TYPE",CVAL = ENKF_TYPE)
    IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ENKF_TYPE: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ENKF_TYPE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
   
!----------------------------------------------------------------------------|
!   EL_ASSIM  -  OPTION FOR CHOSING ELEVATION AS ASSIMILATION VARIABLES       
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"EL_ASSIM",LVAL = EL_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING EL_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE EL_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF   
   
!----------------------------------------------------------------------------|
!   EL_OBS  -  ELEVATION OBSERVATION DATA OPTION                              
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"EL_OBS",LVAL = EL_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING EL_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE EL_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF   

!----------------------------------------------------------------------------|
!   UV_ASSIM  -  OPTION FOR CHOSING CURRENT AS ASSIMILATION VARIABLES         
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"UV_ASSIM",LVAL = UV_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING UV_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE UV_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!----------------------------------------------------------------------------|
!   UV_OBS  -  CURRENT OBSERVATION DATA OPTION                                
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"UV_OBS",LVAL = UV_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING UV_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE UV_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!----------------------------------------------------------------------------|
!   T_ASSIM  -   OPTION FOR CHOSING TEMPERATURE AS ASSIMILATION VARIABLES     
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"T_ASSIM",LVAL = T_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING T_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE T_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!----------------------------------------------------------------------------|
!   T_OBS  -  TEMPERATURE OBERVATION DATA OPTION                              
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"T_OBS",LVAL = T_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING T_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE T_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!----------------------------------------------------------------------------|
!   S_ASSIM  -  OPTION FOR CHOSING SALINITY AS ASSIMILATION VARIABLES         
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"S_ASSIM",LVAL = S_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING S_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE S_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF

!----------------------------------------------------------------------------|
!   S_OBS  -  SALINITY OBSERVATION DATA OPTION                                
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"S_OBS",LVAL = S_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING S_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE S_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF

!==============================================================================|
!            SCREEN REPORT OF SET ENKF VARIABlES                        !
!==============================================================================|
   IF(MSR) THEN  
     WRITE(IPT,*) '!                                                    !'     
     WRITE(IPT,*) '!------SPECIFY ENKF DATA ASSIMINATION PARAMETERS-----!'     
     WRITE(IPT,*) '!                                                    !'     
     WRITE(IPT,*) '!  # SPECIFICY INITIAL PERTUBATION FIELD             :',ENKF_INIT
     WRITE(IPT,*) '!  # GLOBAL NUMBER OF ENSEMBLES                      :',ENKF_NENS
     WRITE(IPT,*) '!  # ASSIMILATON TIME INTERVAL IN SECONDS            :',DELTA_ASS
     WRITE(IPT,*) '!  # ASSIMILATION TIME INTERVAL/FILE OUTPUT INTERVAL :',ENKF_INT
     WRITE(IPT,*) '!  # MAXIMUM NUMBER OF THE OBSERVATION STATIONS      :',ENKF_NOBSMAX
     WRITE(IPT,*) '!  # ASSIMILATION START TIME                         :',ENKF_START
     WRITE(IPT,*) '!  # ASSIMILATION END TIME                           :',ENKF_END
     WRITE(IPT,*) '!  # TIDAL AMPLITUDE ERROR RANGE SPECIFIED           :',(BC_AMP_ERR(I),I=1,6)
     WRITE(IPT,*) '!  # TIDAL PHASE ERROR RANGE SPECIFIED               :',(BC_PHA_ERR(I),I=1,6) 
     WRITE(IPT,*) '!  # MAX DISTANCE OF CORRELATIN                      :',ENKF_CINF
     WRITE(IPT,*) '!  # ELEVATION AS ASSIMILATION VARIABLES             :',EL_ASSIM
     WRITE(IPT,*) '!  # ELEVATION AS OBSERVATION DATA OPTION            :',EL_OBS
     WRITE(IPT,*) '!  # CURRENTS AS ASSIMILATION VARIABLES              :',UV_ASSIM
     WRITE(IPT,*) '!  # CURRENTS OBSERVATION DATA OPTION                :',UV_OBS
     WRITE(IPT,*) '!  # TEMPERATURE AS ASSIMILATION VARIABLES           :',T_ASSIM
     WRITE(IPT,*) '!  # TEMPERATURE OBSERVATION DATA OPTION             :',T_OBS
     WRITE(IPT,*) '!  # SALINITY AS ASSIMILATION VARIABLES              :',S_ASSIM
     WRITE(IPT,*) '!  # SALINITY OBSERVATION DATA OPTION                :',S_OBS
     IF(EL_OBS) WRITE(IPT,*) '!  # EL OBSERVATION ERROR SPECIFIED                  :',OBSERR_EL
     IF(UV_OBS) WRITE(IPT,*) '!  # UV OBSERVATION ERROR SPECIFIED                  :',OBSERR_UV
     IF(T_OBS)  WRITE(IPT,*) '!  # T OBSERVATION ERROR SPECIFIED                   :',OBSERR_T
     IF(S_OBS)  WRITE(IPT,*) '!  # S OBSERVATION ERROR SPECIFIED                   :',OBSERR_S
     WRITE(IPT,*) '!                                                   !'           
  ENDIF
   
   RETURN
   END SUBROUTINE SET_ENKF_PARAM

   SUBROUTINE ALLOC_VARS_ENKF
   
   USE LIMS
   USE CONTROL
   IMPLICIT NONE  
    
   INTEGER         STDIM               !!dimension of state vector,number of elements for all assimilation
   INTEGER         NDB9 
   REAL(DP)  ::    MEMCNT9 
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1

   MEMCNT9 = 0.0_DP
   NDB9    = 2
   
   ALLOCATE(MKF(STDIM))     ;MKF       = ZEROD   !!mean of state Vector forecast over Nens ensembles
   ALLOCATE(AFSTD(STDIM))   ;AFSTD     = ZEROD   !!standard difference of state vector Akf
   ALLOCATE(STTEMP(STDIM))  ;STTEMP    = ZEROD   !! = Mkf             
   ALLOCATE(ERRVEC(STDIM))  ;ERRVEC    = ZEROD   !!the difference between the mean of ensemble forecast and the true
   ALLOCATE(STMEAN(STDIM))  ;STMEAN    = ZEROD   !!stmean = Mkf+Kgain*(obsdata(k) - H*Mkf) 
                                               ! anlysis vector for mean of ensemble
               MEMCNT9 = MEMCNT9+STDIM*6*NDB9+NLOC

   ALLOCATE(SFD(ENKF_NENS)) ;SFD       = ZEROD   !!The singular values of of SfSf
   ALLOCATE(DUMMY(NLOC))    ;DUMMY     = ZEROD   !!input for subroutine gaussj
   ALLOCATE(OBSDATA(NLOC))  ;OBSDATA   = ZEROD   !!true Obs.         
   ALLOCATE(RPAVE(NLOC))    ;RPAVE     = ZEROD   !! mean of perturbations
   ALLOCATE(MODDATA(NLOC))  ;MODDATA   = ZEROD   !!1.innovation vector y' = H(x_fct)-H(x_obs)
                                               ! 2.model data in observation location 
   ALLOCATE(DIFMD(NLOC))    ;DIFMD     = ZEROD   !!inovation(y(i)-H(xf(i))) for ensemble memeber

               MEMCNT9 = MEMCNT9+(ENKF_NENS+NLOC*4+NLOC)*NDB9

   ALLOCATE(WK(NLOC,NLOC))  ;WK        = ZEROD   !!observation error covariance matrix
   ALLOCATE(AKF(STDIM,ENKF_NENS)) ;AKF = ZEROD   !!1.state vector of ensemble forecasts,
                                                 ! 2.and the difference whith it's mean
                                                 !   Akf = Akf -Mkf   
   ALLOCATE(SFH(ENKF_NENS,NLOC))  ;SFH = ZEROD   !!inflation factor,SfH=H(Xf)/sqrt(R*(Ne-1)) 
   ALLOCATE(SA(STDIM,ENKF_NENS+ENKF_NOBSMAX))   ;SA        = ZEROD
   ALLOCATE(SF(STDIM,ENKF_NENS+ENKF_NOBSMAX))   ;SF        = ZEROD
   ALLOCATE(SFH1(ENKF_NENS))          ;SFH1      = ZEROD
   
   
               MEMCNT9 = MEMCNT9+(NLOC*NLOC+STDIM*ENKF_NENS+ENKF_NENS*NLOC)*NDB9

   ALLOCATE(SFSF(ENKF_NENS,ENKF_NENS)) ;SFSF  = ZEROD   !!SfSf=SfH*SfH'=SfU*SfD*SfU'
   ALLOCATE(SFU(ENKF_NENS,ENKF_NENS))  ;SFU   = ZEROD   !!
   ALLOCATE(SK(NLOC,ENKF_NENS))        ;SK    = ZEROD   !!H*(Af-Mf)
   ALLOCATE(BK(NLOC,STDIM))            ;BK    = ZEROD   !!Bk=Sk*Akf^{T}/(Nens-1)
                                                        !   =Sk*(Af-Mf)^{T}/(Nens-1)
                                                        !   =H*(Af-Mf)*(Af-Mf)^{T}/(Nens-1)
                                                        !   =H*Pf
                                                        ! BK=BK*HBHT 
                                                        !    the HBHT if function of distance, is 0 at long distancd
 
               MEMCNT9 = MEMCNT9+(2*ENKF_NENS*ENKF_NENS+NLOC*ENKF_NENS+NLOC*STDIM)*NDB9
 
   ALLOCATE(RK(NLOC,NLOC))             ;RK    = ZEROD   !!Rk=H*Pf*H^{T}=H*Bk^{T}
   ALLOCATE(TK(NLOC,NLOC))             ;Tk    = ZEROD   !!Tk = Rk+ Wk , and
                                                        ! Tk = Tk^{-1}=(H*Pf*H'+Wk)^{-1}
   ALLOCATE(KGAIN(STDIM,NLOC))         ;KGAIN = ZEROD    !!Kalman Gain Matrix,
                                                        ! Kgain = Bk'*Tk
                                                        !       = Pf*H'*(H*Pf*H'+Wk)^{-1}
   
   ALLOCATE(RPERT(ENKF_NENS,NLOC))     ;RPERT = ZEROD   !!random matrix, which mean equas zero

       MEMCNT9 = MEMCNT9+(2*NLOC*NLOC+STDIM*NLOC+ENKF_NENS*NLOC)*NDB9
  
   
   RETURN
   END SUBROUTINE ALLOC_VARS_ENKF


   SUBROUTINE DEALLOC_VARS_ENKF
   USE LIMS

   DEALLOCATE(MKF,AFSTD,STTEMP,ERRVEC,STMEAN)
   DEALLOCATE(SFD,DUMMY,OBSDATA,RPAVE,MODDATA,DIFMD)
   DEALLOCATE(WK,AKF,SFH,SFSF,SFU,SK,BK,SA,SFH1,SF)
   DEALLOCATE(RK,TK,KGAIN,RPERT)
   RETURN
   END SUBROUTINE DEALLOC_VARS_ENKF
   
   SUBROUTINE SET_ASSIM_ENKF

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   use mod_obcs
   IMPLICIT NONE

   INTEGER I,J,K,IERR
   CHARACTER(LEN=100) MKANLDIR, MKFCTDIR, MKERRDIR, MKOUTDIR
   CHARACTER(LEN=120)  ::  FLNAME
   CHARACTER(LEN=4)    ::  IFIL
   CHARACTER(LEN=5)    ::  IFIL1
   
   IF(MSR) THEN 
#  if !defined (DOS)
      MKANLDIR = "mkdir -p "//TRIM(OUTDIR)//"/anl"
      MKFCTDIR = "mkdir -p "//TRIM(OUTDIR)//"/fct"
      MKERRDIR = "mkdir -p "//TRIM(OUTDIR)//"/out_err"  
      MKOUTDIR = "mkdir -p "//TRIM(OUTDIR)//"/flow"
#     if !defined (CRAY)
         CALL SYSTEM( TRIM(MKANLDIR) )  
         CALL SYSTEM( TRIM(MKFCTDIR) )
         CALL SYSTEM( TRIM(MKERRDIR) )
         CALL SYSTEM( TRIM(MKOUTDIR) )
#     endif
#     if defined (CRAY)
         CALL CRAY_SYSTEM_CALL(TRIM(MKANLDIR))
         CALL CRAY_SYSTEM_CALL(TRIM(MKFCTDIR))
         CALL CRAY_SYSTEM_CALL(TRIM(MKERRDIR))    
         CALL CRAY_SYSTEM_CALL(TRIM(MKOUTDIR))
#     endif             
#  endif         
   ENDIF

   ICYC = 0
   I_INITIAL = IINT            !! BE CAREFULE ABOUT IINT, IT SHOULD BE RELATED TO THE INITIAL FIELD YOU BUILD

   NCYC  = IEND/DELTA_ASS      !! GET THE ASSIMILATON LOOP NUMBERS           
   N1CYC = (I_INITIAL+DELTA_ASS)/DELTA_ASS

   IF(MSR) THEN
     OPEN(INOOB,FILE=TRIM(INPDIR)//"/"//trim(casename)//"_assim_enkf.dat",FORM='FORMATTED')

     OPEN(74,FILE=TRIM(OUTDIR)//'/out_err/EnSp.dat')
     WRITE(74,*) ' Icyc       aa        infl1        infl2      inflold    inflation' 
   
     OPEN(75,FILE='./err.dat')

   ENDIF
   IF(ENKF_INIT == 'default') CALL SET_INI_ENKF   ! only used for test case!
 
   IF(MSR) THEN
      IF(IBCN_GL(1) > 0) THEN
        CALL PERT_BC 
      ENDIF
   ENDIF
   
#  if defined (MULTIPROCESSOR)
      IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#  endif 

   RETURN
   END SUBROUTINE SET_ASSIM_ENKF

   SUBROUTINE ENKF_ASS

   USE LIMS
   USE CONTROL
   USE ALL_VARS
#  if defined(WET_DRY)
   USE MOD_WD
#  endif
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER I, J, K, JJ, IERR, ENKF_WD
   INTEGER STDIM
   CHARACTER(LEN=120) FNAM, GNAM    ! gnam : directory for get binary data
                                    ! fnam : directory for store binary data 
   CHARACTER(LEN=4)  FENS   
   

   INTEGER IDOBS, IDOBS2            ! idobs,idobs2 : input number for random number generator
   REAL(DP)  AA, BB, DELT           ! delt =distst :the distance between two location
   REAL*8    DISTST
   REAL(DP)  HBHT                   ! HBHT: for localization 
   REAL(DP)  RNOBS, GASDEV          ! RNobs= gasdev : random number
   REAL(DP)  AAA,SUM0,RSCALE
   REAL(DP)  ENKF_CINF2
   CHARACTER(LEN=80) TEXT 
   REAL(DP)  SUM9, AVGRMSERR,FCTRMSERR,ANLRMSERR,FCTOBSERR,ANLOBSERR
                                    ! avgrmserr : averaged rms error for all ensemble member
                                    ! fctrmserr : rms error of the mean of all ensemble memeber forcast
                                    ! anlrmserr : rms error of the mean of all ensemble memeber analysis
                                    ! fctobserr : rms error between forcast and observation
                                    ! anlobserr : rms error between analysis and boservation
   
   REAL(DP) INFLATION,INFL1,INFL2,INFL1_2,INFLOLD
   REAL(DP) VT
   REAL(DP) ELSTD, USTD, TSTD, SSTD
   INTEGER  LWORK4, LDVT, RCODE
   INTEGER  IDUMMY 
   REAL(DP),ALLOCATABLE    ::   WORK4(:)
   INTEGER ,ALLOCATABLE    ::   ISWDN(:)
   INTEGER ,ALLOCATABLE    ::   ISWDC(:)
   INTEGER ,ALLOCATABLE    ::   ISWD(:,:)


   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1

   ALLOCATE(STFCT(STDIM))   ;STFCT     = ZEROD   !!state vector of one ensemble forecast

!   IF (MSR) THEN  !  0

   LDVT   = 1
   LWORK4 = 5*ENKF_NENS
   ALLOCATE(WORK4(LWORK4))  ; WORK4 = ZEROD
   ALLOCATE(WKTMP(STDIM))   ; WKTMP = ZEROD
   ALLOCATE(STTR(STDIM))    ; STTR  = ZEROD   !!STTR : true state vector 
   ALLOCATE(STTR1(STDIM))   ; STTR1 = ZEROD 
   ALLOCATE(STLOC(ENKF_NOBSMAX))   ; STLOC = 0   !!STLOC: 
   ALLOCATE(EL_INV(1:NNodeGL))         ; EL_INV = 0.0_SP      
   ALLOCATE(ISWDN(NNodeGL))  ; ISWDN = 0
   ALLOCATE(ISWDC(MElemGL))  ; ISWDC = 0
   ALLOCATE(ISWD(STDIM,ENKF_NENS))  ; ISWD = 0


   IF(ICYC == ENKF_START/DELTA_ASS) THEN
      ALLOCATE(STSD(STDIM))        ; STSD  = ZEROD
   ENDIF
   
!   ENDIF ! MSR 0

#  if defined(WET_DRY)
   
   DO K= 1, ENKF_NENS

       WRITE(FENS,'(I4.4)') K
       GNAM= TRIM(OUTDIR)//'/fct/restart'//FENS//'_wd.dat'
       IF(WET_DRY_ON) CALL WD_READ(GNAM)   


!   IF(MSR)THEN  
    
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1, NNodeGL
       IDUMMY = IDUMMY + 1
       ISWD(IDUMMY,K) = ISWET_NODE_CURRENTSTEP(I) 
     ENDDO
   ENDIF
   
   IF(UV_ASSIM) THEN
     DO J=1, KBM1
       DO I=1, MElemGL
         IDUMMY = IDUMMY + 1
         ISWD(IDUMMY,K) = ISWET_CELL_CURRENTSTEP(I) 
       ENDDO
     ENDDO
     DO J=1, KBM1
       DO I=1, MElemGL
         IDUMMY = IDUMMY + 1
         ISWD(IDUMMY,K) = ISWET_CELL_CURRENTSTEP(I)
       ENDDO
     ENDDO
   ENDIF
   
   IF(T_ASSIM) THEN
     DO J=1, KBM1
       DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         ISWD(IDUMMY,K) = ISWET_NODE_CURRENTSTEP(I) 
       ENDDO
     ENDDO
   ENDIF
   
   IF(S_ASSIM) THEN
     DO J=1, KBM1
       DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         ISWD(IDUMMY,K) = ISWET_NODE_CURRENTSTEP(I) 
       ENDDO
     ENDDO
   ENDIF

!   END IF        ! msr

   ENDDO
   
#  else
   ISWD = 1    ! assume all wet
#  endif
   
   
   CALL GETOBSLOC1

   CALL ALLOC_VARS_ENKF

!CC----------------------------------------------------------CC
!CC Get the Observation Covariance Matrix: Wk Nobs * Nobs    CC
!CC----------------------------------------------------------CC

   WK = 0.0_DP
   DO I=1, NLOC
      WK(I,I) = WKTMP(STLOC(I))
      WK(I,I) = WK(I,I)**2.0_DP      
   ENDDO

   AKF   = ZEROD
   STFCT = ZEROD

   DO K=1,ENKF_NENS

      WRITE(FENS,'(I4.4)') K
             
      GNAM= TRIM(OUTDIR)//'/fct/restart'//FENS//'.dat'

      OPEN(INOKF,FILE=TRIM(GNAM),FORM='UNFORMATTED',STATUS='OLD')
      CALL GR2ST(INOKF)        ! RETURN TO STFCT
      DO I=1,STDIM
        AKF(I,K)=STFCT(I)
      ENDDO
   ENDDO
   
!CC--------Calculate the Ensemble Mean and Anomalies-------CC            
!CC the Notation is consistent with that by                CC
!CC Evensen and Van Leeuwen, MWR, 1996. P85--              CC 

   MKF   = ZEROD
   AFSTD = ZEROD
   DO I = 1, STDIM
     AA = 0.0_DP
    
     DO K= 1, ENKF_NENS
       AA = AA + AKF(I,K)      
     ENDDO

       MKF(I)= AA / DBLE(ENKF_NENS)

     BB = 0.0_DP
     DO K= 1, ENKF_NENS

        AKF(I,K) = AKF(I,K) - MKF(I)

     ENDDO
       
   ENDDO    
   
!
!  Calculate sum of eigenvalues for the inflation factor
!  according to Wang and Bishop (2003) in JAS
!
   SFH  = ZEROD
   SFSF = ZEROD
   DO J = 1, NLOC
     DO I=1,ENKF_NENS
        SFH(I,J)=AKF(STLOC(J),I)/DSQRT(WK(J,J)*DBLE((ENKF_NENS-1)))
     ENDDO
   ENDDO

   DO I=1,ENKF_NENS
     DO J=1,ENKF_NENS
        SUM9 = 0.0_DP
        DO K=1,NLOC
           SUM9 = SUM9 + SFH(I,K)*SFH(J,k)	   
        ENDDO
        SFSF(I,J) = SUM9
     ENDDO
   ENDDO

! SVD of SfSf: SfSf=(psi)*R^-1*(psi)'=SfH*SfH'=SfU*SfD*SfU'
!input:SFSF      , output SFU,SFD
   
   CALL DGESVD('A','N',ENKF_NENS,ENKF_NENS,SFSF,ENKF_NENS,SFD,SFU,ENKF_NENS,VT,LDVT,WORK4,LWORK4,RCODE)
   
!sum_of_lambda
   OPEN(INOKF,FILE=TRIM(OUTDIR)//'/out_err/Lambda'//fcyc//'.dat')
   INFL1 = 0.0_DP
   DO J=1,ENKF_NENS
     INFL1 = INFL1 + SFD(J)
     WRITE(INOKF,'(I5,E15.7)') ICYC,SFD(J)
   ENDDO
   CLOSE(INOKF)

!---------- Calculate Sk(1:Nobs,1:Nens)-------------------------C
!-----------   Sk=H_{k}*(Af_{k}-Mkf_{k})=Nobs*Nens     ---------C
   SK = ZEROD        
   DO K=1,ENKF_NENS
     DO I=1,NLOC
       SK(I,K) = AKF(STLOC(I),K)
     ENDDO
   ENDDO

!CC--------------- Calculate Bk(1:Nobs,1:stdim) -------------CC
!CC  Bk=Sk*Akf^{T}/(Nens-1)=Sk*(Af_{k}-Mf_{k})^{T}/(Nens-1)--CC
!CC    =H_{k}*(Af_{k}-Mkf_{k})*(Af_{k}-Mf_{k})^{T}/(Nens-1)
!CC    =H_{k}*Pf 
!CC---------- Bk=Nobs * stdim -------------------------------CC
   Bk = ZEROD  
   DO I=1,NLOC
     DO K=1,STDIM
        AA = 0.0_DP
	SUM0=0.0_DP 
        DO J=1,ENKF_NENS 
           IF(ISWD(K,J)==1) THEN
	     AA = AA + SK(I,J)*AKF(K,J)
	     SUM0=SUM0+1
	   ENDIF

        ENDDO
        IF(SUM0>0) BK(I,K)=AA/DBLE(SUM0)  
     ENDDO

   ENDDO
   
!CC---------------------------------------------------CC
!CC---- Correct the Covariance Function; get rid of   CC 
!CC--  false correlation at long distance          ---CC
!CC--------  Bk=Bk*HBHT;  ----------------------------CC
!CC-----See Eq. 4.10 in Gaspari and Cohn (1999) in Quart. J. Royal Meteor. Soc.
!CC

   IF(ABS(ENKF_CINF)>0.001) THEN
    ENKF_CINF2 = DBLE(ENKF_CINF)/2.0_DP

    DO I=1,NLOC

      DO K=1,STDIM
 
        DELT = DISTST(STLOC(I),K)        

        IF(DELT > ENKF_CINF) THEN
           HBHT = 0.
        ELSE IF(DELT > ENKF_CINF2) THEN
           HBHT = 1.0_DP/12.0_DP*(DELT/ENKF_CINF2)**5 &
                    -1.0_DP/2.0_DP*(DELT/ENKF_CINF2)**4 &
                    +5.0_DP/8.0_DP*(DELT/ENKF_CINF2)**3 &
                    +5.0_DP/3.0_DP*(DELT/ENKF_CINF2)**2 &
                    -5.0_DP*(DELT/ENKF_CINF2) &
                    +4.0_DP -2.0_DP/3.0_DP*(ENKF_CINF2/DELT)
        ELSE
           HBHT = -1.0_DP/4.0_DP*(DELT/ENKF_CINF2)**5 &
                    +1.0_DP/2.0_DP*(DELT/ENKF_CINF2)**4 &
                    +5.0_DP/8.0_DP*(DELT/ENKF_CINF2)**3 &
                    -5.0_DP/3.0_DP*(DELT/ENKF_CINF2)**2 &
                    +1
        ENDIF
    
        BK(I,K)=BK(I,K)*HBHT    

     ENDDO
   ENDDO
   ENDIF

!CC---------  Calculate Rk(1:Nobs,1:Nobs)  -----------CC
!CC---------     Note Rk=H*Pf*H^{T}=H*Bk^{T} ------------------CC

   Rk = ZEROD
   Tk = ZEROD
   DO I=1,NLOC
     DO J=1,NLOC
        RK(I,J) = BK(J,STLOC(I))               
     ENDDO
   ENDDO
   
!CC----------------------------------------------------CC
!CC                  Tk= Rk+ Wk                        CC
!CC----------------------------------------------------CC

   DO I=1,NLOC
     DO J=1,NLOC
        TK(I,J) = RK(I,J) + WK(I,J)
     ENDDO
   ENDDO

   DO I=1,NLOC
      DUMMY(I)= 0.0_DP
   ENDDO

   CALL GAUSSJ(TK,NLOC,NLOC,DUMMY,1,1)   

!C-----------------------------------------------------C
!CC- Calculate the Kalman Gain Matrix : Akf(Nxy, Nobs) C
!C-----------------------------------------------------C
   KGAIN = ZEROD
   DO I=1, STDIM
      DO J=1, NLOC
         AA = 0.0_DP
         DO K=1, NLOC
            AA = AA + BK(K,I)*TK(K,J)
         ENDDO
         KGAIN(I,J) = AA
      ENDDO
   ENDDO 
   
!C-----------------------------------------------------C
!CC----------- Get the Observations         ----------CC
!C-----------------------------------------------------C

    CALL READMEDM                       

   DO I=1,NLOC
      OBSDATA(I) = STTR1(STLOC(I))
   ENDDO
   
!C---------Calculate the difference between the mean of -C
!C-------- ensemble forecast and the true ---------------C
! 
   DO I=1,NLOC  
      ERRVEC(I) = MKF(STLOC(I)) - OBSDATA(I) 
   ENDDO
   TEXT='FctObserr'
   CALL PRINT_ERR(ERRVEC,NLOC,TEXT,FCTOBSERR)

   DO I=1,STDIM
      ERRVEC(I) = MKF(I) - STTR1(I)
   ENDDO
   TEXT='Fctrmserr'
   CALL PRINT_ERR(ERRVEC,STDIM,TEXT,FCTRMSERR)


!C-------------------------------------------------------C
!C---------        Assimilation  ------------------------C
!C-------------------------------------------------------C

!C--- Initialize the Random Number Generator for Obs. ---C
!C--- Initial value for IDobs must be negative        ---C

   IF (ICYC == ENKF_START/DELTA_ASS) THEN
      IDOBS  = -31
      IDOBS2 = -711
   ELSE
      OPEN(INOKF,FILE = TRIM(OUTDIR)//'/out_err/IDobstmp.dat')
      READ(INOKF,*) IDOBS
!        print *, 'IDobs= ', IDobs
      CLOSE(INOKF)

      OPEN(INOKF,FILE = TRIM(OUTDIR)//'/out_err/IDobs2tmp.dat')
      READ(INOKF,*) IDOBS2
!         print *, 'IDobs2= ', IDobs2
      CLOSE(INOKF)

   ENDIF

       
!CC----------- Generate  Observations  from True State-------------CC
!CC-----------   By adding noises to true state       -------------CC

   DO I=1, NLOC
      RNOBS = GASDEV(IDOBS)
      OBSDATA(I) = OBSDATA(I) + DSQRT(WK(I,I))*RNOBS
   ENDDO      
       
!C---------------------------------------------------------------------C
!C---  Remove the mean of perturbations added to the observations -----C
!C---------------------------------------------------------------------C

   DO K=1, ENKF_NENS
      DO I=1, NLOC
         RNOBS = GASDEV(IDOBS2)
         RPERT(K,I)= DSQRT(WK(I,I))*RNOBS
      ENDDO
   ENDDO

   DO I=1, NLOC
      AAA = 0.0_DP
      DO K=1, ENKF_NENS
         AAA = AAA + RPERT(K,I)
      ENDDO
      RPAVE(I) = AAA/DBLE(ENKF_NENS)
      DO K=1, ENKF_NENS
         RPERT(K,I) = RPERT(K,I) - RPAVE(I)
      ENDDO
   ENDDO

!C------------------------------------------------------C
!C------------------------------------------------------C
!C------------------------------------------------------C

!c
!c Calculate innovation vector y' = H(x_fct)-H(x_obs) ->moddata
!c
   INFL2 = 0.0_DP
   DO K = 1, NLOC
      MODDATA(K) = OBSDATA(K) - MKF(STLOC(K))
      INFL2 = INFL2 + MODDATA(K)**2/WK(K,K)
   ENDDO
   
   DO I=1, STDIM
      BB = 0.0_DP
      DO J=1,NLOC
         BB = BB + KGAIN(I,J)*MODDATA(J)
      ENDDO
      STMEAN(I) = MKF(I) + BB 
   ENDDO

   
!====================================   
!  CALCUATE XA-XABAR USING ENSRKF   |
!====================================
  
   ENKF_CINF2 = DBLE(ENKF_CINF)/2.0_DP
   
   DO I = 1, STDIM
     DO J = 1, ENKF_NENS
      SF(I,J) = AKF(I,J)/DSQRT(DBLE(ENKF_NENS-1))
     ENDDO
   ENDDO
 
   DO J = 1, NLOC

      DO I=1,ENKF_NENS  
         SFH1(I) = SF(STLOC(J),I)
      ENDDO
      SUM0 = 0.0D0
      DO I=1,ENKF_NENS
         SUM0 = SUM0 + SFH1(I)**2
      ENDDO
      SUM0 = SUM0 + WK(J,J)
      SUM0 = SUM0 + DSQRT(WK(J,J)*SUM0) 
      DO I = 1,ENKF_NENS
         DO K = 1, ENKF_NENS
            SFSF(I,K) = SFH1(I)*SFH1(K)
         ENDDO
      ENDDO
      DO I=1,ENKF_NENS
         DO K = 1, ENKF_NENS
            SFU(I,K)= 1.0_SP/SUM0*SFSF (I,K)
         ENDDO
      ENDDO

       DO I = 1, STDIM
         RSCALE = 1
       
       IF(ABS(ENKF_CINF)>0.001) THEN
         DELT = DISTST(STLOC(J),I)
                
         IF(DELT > ENKF_CINF) THEN
           RSCALE = 0.
         ELSE IF(DELT > ENKF_CINF2) THEN
           RSCALE = 1.0_DP/12.0_DP*(DELT/ENKF_CINF2)**5 &
                    -1.0_DP/2.0_DP*(DELT/ENKF_CINF2)**4 & 
                    +5.0_DP/8.0_DP*(DELT/ENKF_CINF2)**3 &
                    +5.0_DP/3.0_DP*(DELT/ENKF_CINF2)**2 &
                    -5.0_DP*(DELT/ENKF_CINF2) &
                    +4.0_DP -2.0_DP/3.0_DP* (ENKF_CINF2/DELT)
         ELSE
           RSCALE = -1.0_DP/4.0_DP*(DELT/ENKF_CINF2)**5 &
                    +1.0_DP/2.0_DP*(DELT/ENKF_CINF2)**4 &
                    +5.0_DP/8.0_DP*(DELT/ENKF_CINF2)**3 & 
                    -5.0_DP/3.0_DP*(DELT/ENKF_CINF2)**2 &
                    +1
         ENDIF
       ENDIF
	 
        DO JJ = 1, ENKF_NENS
           SUM0 = 0.0D0
           DO K = 1, ENKF_NENS
              SUM0 = SUM0 + SF(I,K)*SFU(K,JJ)
           ENDDO
           SA(I,JJ) = SF(I,JJ) - SUM0*RSCALE      

        ENDDO
       ENDDO
       
       DO I = 1,STDIM
         DO JJ = 1,ENKF_NENS
	   SF(I,JJ) = SA(I,JJ)
	 END DO
       END DO	   
      
   ENDDO    

! ========= END OF XA CALCUALATION ==============   

   
!  Calculate an inflation factor=[innovation'*innovation-Nloc]/sum_of_lambda

   INFLATION = (INFL2 - NLOC)/INFL1
  
   IF(INFLATION < INFLREF) THEN   
      INFLOLD = INFLATION
      INFLATION = 1.0_DP
   ELSE
      INFLOLD = INFLATION
   ENDIF

!   INFL1_2 = DSQRT(INFLATION)    ! TURN ON  INFLATION
   INFL1_2 = 1.0_DP               ! TURN OFF INFLATION

   DO K=1, ENKF_NENS

!!CC----------- Perturb the observations  -------------CC
      DO I = 1,NLOC
        OBSDATA(I) = STTR1(STLOC(I))
      END DO	
 
      DO I=1, NLOC
         OBSDATA(I) = OBSDATA(I) + RPERT(K,I)
      ENDDO

!CC----------------------------------------------------CC       
       
      DO I=1,STDIM
         AKF(I,K) = AKF(I,K) + MKF(I)       
      ENDDO
      
      DO I=1,NLOC
         MODDATA(I) = AKF(STLOC(I),K)
         DIFMD(I)   = OBSDATA(I) - MODDATA(I)
      ENDDO      
      
      DO I=1,STDIM
         BB = 0.0_DP
         DO J= 1, NLOC
            BB = BB + KGAIN(I,J)*DIFMD(J)
         ENDDO
	 
	 IF(ENKF_TYPE == 'default') THEN
           AKF(I,K) = (AKF(I,K)+BB-STMEAN(I))*INFL1_2 + STMEAN(I)             ! EnKF
         ENDIF

         IF(ENKF_TYPE == 'square root')	THEN 
           AKF(I,K) = SA(I,K)*DSQRT(DBLE(ENKF_NENS-1.))*INFL1_2 + STMEAN(I)   ! EnSR   
         ENDIF
	 
         STFCT(I) = AKF(I,K)	 

      ENDDO 
      	   
      WRITE(FENS,'(I4.4)') K
            
      FNAM = TRIM(OUTDIR)//'/fct/restart'//fens//'.dat'            
      GNAM = TRIM(OUTDIR)//'/anl/restart'//fens//'.dat'            
!---------------------------------------------------------------------------|
!    read fvcom forecast for every ensembles from fct/                      |
!    write EnKF anlysis to anl/                                             |
!---------------------------------------------------------------------------|

       CALL FCT2ANL(FNAM,GNAM)

#      if defined(WET_DRY)
       FNAM = TRIM(OUTDIR)//'/anl/restart'//fens//'_wd.dat' 
       CALL WET_JUDGE_EL     
       OPEN(INOKF,FILE=FNAM,FORM='FORMATTED')
       CALL WD_DUMP_EL(INOKF,IINT-1)
#      endif 

     ENDDO

     IINT=IINT-1
     
#  if defined(WET_DRY)

     DO K=1,ENKF_NENS
       WRITE(FENS,'(I4.4)') K
       GNAM= TRIM(OUTDIR)//'/anl/restart'//FENS//'_wd.dat'
       IF(WET_DRY_ON) CALL WD_READ(GNAM) 

!   IF(MSR)THEN  

   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1, NNodeGL
       IDUMMY = IDUMMY + 1
       ISWD(IDUMMY,K) = ISWET_NODE_CURRENTSTEP(I) 
     ENDDO
   ENDIF
   
   IF(UV_ASSIM) THEN
     DO J=1, KBM1
       DO I=1, MElemGL
         IDUMMY = IDUMMY + 1
         ISWD(IDUMMY,K) = ISWET_CELL_CURRENTSTEP(I) 
       ENDDO
     ENDDO
     DO J=1, KBM1
       DO I=1, MElemGL
         IDUMMY = IDUMMY + 1
         ISWD(IDUMMY,K) = ISWET_CELL_CURRENTSTEP(I)
       ENDDO
     ENDDO
   ENDIF
   
   IF(T_ASSIM) THEN
     DO J=1, KBM1
       DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         ISWD(IDUMMY,K) = ISWET_NODE_CURRENTSTEP(I) 
       ENDDO
     ENDDO
   ENDIF
   
   IF(S_ASSIM) THEN
     DO J=1, KBM1
       DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         ISWD(IDUMMY,K) = ISWET_NODE_CURRENTSTEP(I) 
       ENDDO
     ENDDO
   ENDIF

!   END IF   ! msr

   ENDDO   
   
#  else
   ISWD = 1    ! assume all wet
#  endif
   
      
   DO K=1, ENKF_NENS

      IDUMMY = 0
      DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         EL_INV(I) = AKF(IDUMMY,K)
      ENDDO

      IF(ENKF_INIT == 'inv_bc') THEN
        CALL ENKF_INVERSE(K)
      ENDIF
   ENDDO
   
!CC--------Calculate the average rms error of each ensemble member -------CC            

   TEXT = 'eachrmserr'
   SUM9 = 0.0_DP
   DO K=1, ENKF_NENS
      DO I=1,STDIM
         ERRVEC(I) = AKF(I,K) - STTR1(I)
      ENDDO
      CALL PRINT_ERR(ERRVEC,STDIM,TEXT,AVGRMSERR)
      SUM9 = SUM9 + AVGRMSERR
   ENDDO
   AVGRMSERR = SUM9/ENKF_NENS             

!CC--------Calculate the Ensemble Mean -------CC            

   DO I=1,STDIM
      AA = 0.0_DP
      
      DO K=1,ENKF_NENS
         AA = AA + AKF(I,K)
      ENDDO           

         MKF(I)= AA / DBLE(ENKF_NENS)

      BB = 0.0_DP
      DO K=1,ENKF_NENS
         AKF(I,K) = AKF(I,K) - MKF(I)
      ENDDO
  
   ENDDO   
   
!C---------Calculate the difference between the mean of -C
!C-------- ensemble forecast and the true ---------------C
         

!  CALCULATE ANALYSIS RMS ERROR
              
   DO I=1,NLOC
      ERRVEC(I) = MKF(STLOC(I)) - STTR1(STLOC(I))
   ENDDO
   TEXT = 'AnlObserr'
   CALL PRINT_ERR(ERRVEC,NLOC,TEXT,ANLOBSERR)

   ERRVEC = MKF - STTR1
   TEXT = 'Anlrmserr'
   CALL PRINT_ERR(ERRVEC,STDIM,TEXT,ANLRMSERR)
   
 
   WRITE(74,'(i5,5d13.5)') ICYC,AA,INFL1,INFL2,INFLOLD,INFLATION
   WRITE(75,'(i5,6d12.4)') ICYC,FCTOBSERR,FCTRMSERR,ANLOBSERR,ANLRMSERR
!C----------------------------------------------------------C

   OPEN(INOKF, FILE = TRIM(OUTDIR)//'/out_err/IDobstmp.dat',STATUS ='UNKNOWN')
   WRITE(INOKF, *) IDOBS
   CLOSE(INOKF)
!  PRINT *, 'IDobs= ', IDobs

   OPEN(INOKF, FILE= TRIM(OUTDIR)//'/out_err/IDobs2tmp.dat',STATUS ='UNKNOWN')
   WRITE(INOKF, *) IDOBS2
   CLOSE(INOKF)
!  PRINT *, 'IDobs2= ', IDobs2           
   
   CALL DEALLOC_VARS_ENKF
   DEALLOCATE(WORK4, STTR, STTR1, STLOC, WKTMP, EL_INV)


   DEALLOCATE(STFCT) 
   RETURN
   END SUBROUTINE ENKF_ASS


   SUBROUTINE GR2ST(INF)
   
   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
   IMPLICIT NONE
   
   INTEGER I,K,N1,IDUMMY 
   INTEGER INF
   INTEGER ITMP
   
   REAL(SP), ALLOCATABLE :: UTMP(:,:)
   REAL(SP), ALLOCATABLE :: VTMP(:,:)
   REAL(SP), ALLOCATABLE :: S1TMP(:,:)
   REAL(SP), ALLOCATABLE :: T1TMP(:,:)
   REAL(SP), ALLOCATABLE :: ELTMP(:)
   REAL(SP), ALLOCATABLE :: UATMP(:)
   REAL(SP), ALLOCATABLE :: VATMP(:)   
   REAL(SP), ALLOCATABLE :: TMP1(:,:)
   REAL(SP), ALLOCATABLE :: TMP2(:,:)
   REAL(SP), ALLOCATABLE :: TMP3(:,:)
   REAL(SP), ALLOCATABLE :: TMP4(:)
   REAL(SP), ALLOCATABLE :: TMP5(:)
   REAL(SP), ALLOCATABLE :: TMP6(:)
   REAL(SP), ALLOCATABLE :: TMP8(:,:)
#  if defined (WATER_QUALITY)
   REAL(SP), ALLOCATABLE :: TMP7(:,:,:)
#  endif

   ALLOCATE(UTMP(0:MElemGL,1:KB))      ; UTMP  = 0.0_SP 
   ALLOCATE(VTMP(0:MElemGL,1:KB))      ; VTMP  = 0.0_SP
   ALLOCATE(S1TMP(1:NNodeGL,1:KB))     ; S1TMP = 0.0_SP
   ALLOCATE(T1TMP(1:NNodeGL,1:KB))     ; T1TMP = 0.0_SP
   ALLOCATE(ELTMP(1:NNodeGL))          ; ELTMP = 0.0_SP
   ALLOCATE(UATMP(0:MElemGL))          ; UATMP = 0.0_SP 
   ALLOCATE(VATMP(0:MElemGL))          ; VATMP = 0.0_SP
   ALLOCATE(TMP1(0:MElemGL,1:KB))      ; TMP1  = 0.0_SP
   ALLOCATE(TMP2(1:NNodeGL,1:KB))      ; TMP2  = 0.0_SP
   ALLOCATE(TMP8(0:NNodeGL,1:KB))      ; TMP8  = 0.0_SP
   ALLOCATE(TMP3(1:MElemGL,1:KB))      ; TMP3  = 0.0_SP
   ALLOCATE(TMP4(0:MElemGL))           ; TMP4  = 0.0_SP
   ALLOCATE(TMP5(1:MElemGL))           ; TMP5  = 0.0_SP
   ALLOCATE(TMP6(1:NNodeGL))           ; TMP6  = 0.0_SP
#  if defined (WATER_QUALITY)
   ALLOCATE(TMP7(1:NNodeGL,1:KB,1:NB)) ; TMP7  = 0.0_SP 
#  endif

   REWIND(INF)
   READ(INF) ITMP
   READ(INF) ((UTMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INF) ((VTMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,MElemGL)
#  if defined (GOTM)
   READ(INF) ((TMP8(I,K),K=1,KB),I=0,NNodeGL)
   READ(INF) ((TMP8(I,K),K=1,KB),I=0,NNodeGL)
#  else
   READ(INF) ((TMP8(I,K),K=1,KB),I=0,NNodeGL)
   READ(INF) ((TMP8(I,K),K=1,KB),I=0,NNodeGL)
   READ(INF) ((TMP8(I,K),K=1,KB),I=0,NNodeGL)
#  endif
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,MElemGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,MElemGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,MElemGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,MElemGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,MElemGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,MElemGL)

   READ(INF) ((S1TMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INF) ((T1TMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,NNodeGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,NNodeGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,NNodeGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,NNodeGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,NNodeGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,NNodeGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,NNodeGL)

   READ(INF) (UATMP(I), I=0,MElemGL)
   READ(INF) (VATMP(I), I=0,MElemGL)
   READ(INF) (TMP5(I), I=1,MElemGL)
   READ(INF) (TMP5(I), I=1,MElemGL)
   READ(INF) (TMP5(I), I=1,MElemGL)
   READ(INF) (TMP5(I), I=1,MElemGL)
   READ(INF) (TMP5(I), I=1,MElemGL)
   READ(INF) (TMP5(I), I=1,MElemGL)

   READ(INF) (ELTMP(I), I=1,NNodeGL)
   READ(INF) (TMP6(I), I=1,NNodeGL)
   READ(INF) (TMP6(I), I=1,NNodeGL)
   READ(INF) (TMP6(I), I=1,NNodeGL)
   READ(INF) (TMP6(I), I=1,NNodeGL)

#  if defined (EQUI_TIDE)
     READ(INF) (TMP6(I), I=1,NNodeGL)
#  endif
#  if defined (ATMO_TIDE)
     READ(INF) (TMP6(I), I=1,NNodeGL)
#  endif

#  if defined (WATER_QUALITY)
   DO N1=1,NB
     READ(INF) ((TMP7(I,K,N1),K=1,KB),I=1,NNodeGL)
   END DO
#  endif

   CLOSE(INF) 
   
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1, NNodeGL
       IDUMMY = IDUMMY + 1
       STFCT(IDUMMY) = ELTMP(I) 
     ENDDO
   ENDIF
   
   IF(UV_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, MElemGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = UTMP(I,k) 
       ENDDO
     ENDDO
     DO K=1, KBM1
       DO I=1, MElemGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = VTMP(I,k)
       ENDDO
     ENDDO
   ENDIF
   
   IF(T_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = T1TMP(I,K) 
       ENDDO
     ENDDO
   ENDIF
   
   IF(S_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = S1TMP(I,K) 
       ENDDO
     ENDDO
   ENDIF
   
   DEALLOCATE(UTMP,VTMP,S1TMP,T1TMP,ELTMP,UATMP,VATMP,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP8)
#  if defined (WATER_QUALITY)
   DEALLOCATE(TMP7)
#  endif
   RETURN
   END SUBROUTINE GR2ST 


   SUBROUTINE GETOBSLOC1
   
   USE LIMS
   USE CONTROL
   IMPLICIT NONE

   INTEGER ::  NUM 
   INTEGER  SWITCH
   SAVE     SWITCH
   INTEGER ::  J,K,PNT
   INTEGER ::  IDUMMY
   INTEGER ::  TMP
   CHARACTER(LEN=80) FILENAME
   CHARACTER(LEN=24) HEADINFO
   INTEGER STLTMP(ENKF_NOBSMAX)
   INTEGER LAY(ENKF_NOBSMAX)

   NUM     = 0 
   IDUMMY  = 0
   NLOC    = 0
   PNT     = 0
   STLOC   = 0
   STLTMP  = 0 
   LAY     = 0
   STTR    = 0.0_DP
   WKTMP   = 0.0_DP

   SWITCH  = 0
   
   FILENAME = TRIM(INPDIR)//"/"//trim(casename)//"_assim_enkf.dat"
    
   OPEN(INOOB,FILE=TRIM(FILENAME),FORM='FORMATTED')
   
 100 READ(INOOB,*,END=200) HEADINFO
     IF(SWITCH/=1) THEN
       IF(HEADINFO=='!===READ') THEN
         SWITCH = 1
         GOTO 100
       ELSE
         GOTO 100
       ENDIF
     ENDIF 

     IF(TRIM(HEADINFO)=='!EL') THEN
       IF(EL_OBS) THEN
         READ(INOOB,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>ENKF_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INOOB,*)  (STLOC(K), K=1,NLOC)	
	 
	   DO K=1, NLOC   
              WKTMP(STLOC(K)) = OBSERR_EL      ! FIND THE NAME OF VARIABLE
           ENDDO 
	 
       ENDIF

       IF(EL_ASSIM) THEN
	 IDUMMY = IDUMMY + NNodeGL
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!UV') THEN
       IF(UV_OBS) THEN
         READ(INOOB,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC+NUM>ENKF_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc+num, 'Nobsmax=', ENKF_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INOOB,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)
	 READ(INOOB,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MElemGL*(LAY(K)-1)
	   WKTMP(STLOC(K)) = OBSERR_UV
	 ENDDO   
	 
	 NLOC = NLOC + NUM
	 DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K-NUM)+IDUMMY+MElemGL*KBM1
	   WKTMP(STLOC(K)) = OBSERR_UV
	 ENDDO
	 
       ENDIF
       IF(UV_ASSIM) THEN  
          IDUMMY = IDUMMY + 2*MElemGL*KBM1
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!T') THEN
       IF(T_OBS) THEN
         READ(INOob,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>ENKF_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INOOB,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)       
         READ(INOOB,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NNodeGL*(LAY(K)-1)
	   WKTMP(STLOC(K)) = OBSERR_T
	 ENDDO   
        
       ENDIF   
       IF(T_ASSIM) THEN
         IDUMMY = IDUMMY + NNodeGL*1
       ENDIF
     ENDIF

     IF(TRIM(HEADINFO)=='!S') THEN
       IF(S_OBS) THEN
         READ(INOOB,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>ENKF_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INOOB,*)  (STLOC(K),K=NLOC-NUM+1,NLOC)        
         READ(INOOB,*)  (LAY(K),  K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NNodeGL*(LAY(K)-1)
	   WKTMP(STLOC(K)) = OBSERR_S
	 ENDDO   
         
       ENDIF
       IF(S_ASSIM) THEN
         IDUMMY = IDUMMY + NNodeGL*KBM1 
       ENDIF 
     ENDIF

     GOTO 100
 200 CONTINUE
     REWIND(INOOB)
     
    IF(NLOC==0) THEN 
       WRITE(IPT,*) "!WARNING, NOT OBSERVATION DATA FOUND!"
       CALL PSTOP
    ELSE   
       IF(NLOC>1) THEN
         DO J=1, NLOC-1
           DO K=J+1, NLOC
             IF(STLOC(K)<STLOC(J)) THEN
               TMP = STLOC(J)
               STLOC(J) = STLOC(K)
               STLOC(K) = TMP 
             ENDIF
           ENDDO
         ENDDO
       ENDIF
    ENDIF

    RETURN
   END SUBROUTINE GETOBSLOC1


   SUBROUTINE gaussj(a,n,np,b,m,mp)

   USE MOD_PREC
   IMPLICIT NONE
     
   INTEGER m,mp,n,np,NMAX
   REAL(DP) a(np,np),b(np,mp)
   PARAMETER (NMAX=1000)
   INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
   REAL(DP) big,dum,pivinv

   do 11 j=1,n
      ipiv(j)=0
11 continue
   do 22 i=1,n
      big=0.
      do 13 j=1,n
        if(ipiv(j).ne.1)then
          do 12 k=1,n
            if (ipiv(k).eq.0) then
              if (abs(a(j,k)).ge.big)then
                big=abs(a(j,k))
                irow=j
                icol=k
              endif

            else if (ipiv(k).gt.1) then
              pause 'singular matrix in gaussj'
            endif
12       continue
       endif
13   continue
     
     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
       do 14 l=1,n
         dum=a(irow,l)
         a(irow,l)=a(icol,l)
         a(icol,l)=dum
14     continue
       do 15 l=1,m
         dum=b(irow,l)
         b(irow,l)=b(icol,l)
         b(icol,l)=dum
15     continue
       endif
       indxr(i)=irow
       indxc(i)=icol
       if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'

       pivinv=1./a(icol,icol)
       a(icol,icol)=1.
       do 16 l=1,n
         a(icol,l)=a(icol,l)*pivinv
16     continue
       do 17 l=1,m
         b(icol,l)=b(icol,l)*pivinv
17     continue
       do 21 ll=1,n
         if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0.
           do 18 l=1,n
             a(ll,l)=a(ll,l)-a(icol,l)*dum
18         continue
           do 19 l=1,m
             b(ll,l)=b(ll,l)-b(icol,l)*dum
19         continue
         endif
21     continue
22  continue
    do 24 l=n,1,-1
       if(indxr(l).ne.indxc(l))then

         do 23 k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
23       continue
       endif
24  continue
    
    RETURN
    END SUBROUTINE GAUSSJ

   FUNCTION GASDEV(IDUM)
   IMPLICIT NONE
   
   INTEGER IDUM
   REAL GASDEV
!CU    USES ran2
!      REAL ran2
   INTEGER iset
   REAL fac,gset,rsq,v1,v2,ran2
   SAVE iset,gset
   DATA iset/0/

   if (iset.eq.0) then
1     v1=2.*ran2(idum)-1.
      v2=2.*ran2(idum)-1.
      rsq=v1**2+v2**2
      if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
   RETURN
   END FUNCTION GASDEV

   FUNCTION RAN2(IDUM)
   IMPLICIT NONE
   
   INTEGER IDUM,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
   REAL ran2,AM,EPS,RNMX
   PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
              IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
              NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   INTEGER idum2,j,k,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/, iv/NTAB*0/, iy/0/
   
   if (idum.le.0) then
      idum=max(-idum,1)
      idum2=idum
      do 11 j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         if (idum.lt.0) idum=idum+IM1
         if (j.le.NTAB) iv(j)=idum
11    continue
      iy=iv(1)
   endif
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   if (idum.lt.0) idum=idum+IM1
   k=idum2/IQ2
   idum2=IA2*(idum2-k*IQ2)-k*IR2
   if (idum2.lt.0) idum2=idum2+IM2
   j=1+iy/NDIV
   iy=iv(j)-idum2
   iv(j)=idum
   if(iy.lt.1)iy=iy+IMM1
   ran2=min(AM*iy,RNMX)
   return
   END FUNCTION RAN2

   SUBROUTINE ST2GR(FLNAME)

   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif 
#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE

   INTEGER I,J,K,N1
   INTEGER IDUMMY
   INTEGER IERR, NPC, me
   REAL(SP), ALLOCATABLE :: RTP(:)
   CHARACTER(LEN=120)    :: FLNAME
#  if defined (MULTIPROCESSOR)   
     real(sp), allocatable ::  enkfel(:)   
     real(sp), allocatable ::  enkfu(:,:), enkfv(:,:)
     real(sp), allocatable ::  enkft(:,:), enkfs(:,:)
 
     allocate(enkfel(NNodeGL))        ; enkfel = 0.0_sp
     allocate(enkfu(MElemGL,kbm1))    ; enkfu  = 0.0_sp
     allocate(enkfv(MElemGL,kbm1))    ; enkfv  = 0.0_sp
     allocate(enkft(NNodeGL,kbm1))    ; enkft  = 0.0_sp
     allocate(enkfs(NNodeGL,kbm1))    ; enkfs  = 0.0_sp
#  endif 

   ME = MYID ; NPC = NPROCS
   ALLOCATE(RTP(MElemGL)) ; RTP = 0.0_SP
  
   IF(MSR) THEN     
     OPEN(1,FILE=TRIM(FLNAME), FORM='UNFORMATTED') 
     REWIND(1)
     WRITE(1) I_INITIAL
   ENDIF
   
   IF(SERIAL)THEN

   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         EL(I) = STFCT(IDUMMY)
      ENDDO
      DO I=1, MElemGL
         EL1(I)=(EL(NVG(I,1)) + EL(NVG(I,2)) + EL(NVG(I,3)) )/3.0_DP
      ENDDO
   ENDIF

   IF(UV_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MElemGL
          IDUMMY = IDUMMY + 1
          U(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  

      DO K=1, KBM1
        DO I=1, MElemGL
          IDUMMY = IDUMMY + 1
          V(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO
 
   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          T1(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          S1(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   WRITE(1) ((U(I,K),    K=1,KB),I=0,MElem)
   WRITE(1) ((V(I,K),    K=1,KB),I=0,MElem)
   WRITE(1) ((W(I,K),    K=1,KB),I=0,MElem)
#  if defined (GOTM)
   WRITE(1) ((TKE(I,K),   K=1,KB),I=0,NNode)
   WRITE(1) ((TEPS(I,K),  K=1,KB),I=0,NNode)
#  else
   WRITE(1) ((Q2(I,K),   K=1,KB),I=0,NNode)
   WRITE(1) ((Q2L(I,K),  K=1,KB),I=0,NNode)
   WRITE(1) ((L(I,K)  ,  K=1,KB),I=0,NNode)
#  endif
   WRITE(1) ((S(I,K),    K=1,KB),I=0,MElem)
   WRITE(1) ((T(I,K),    K=1,KB),I=0,MElem)
   WRITE(1) ((RHO(I,K),  K=1,KB),I=0,MElem)
   WRITE(1) ((TMEAN(I,K),K=1,KB),I=0,MElem)
   WRITE(1) ((SMEAN(I,K),K=1,KB),I=0,MElem)
   WRITE(1) ((RMEAN(I,K),K=1,KB),I=0,MElem)

   WRITE(1) ((S1(I,K),    K=1,KB),I=1,NNode)
   WRITE(1) ((T1(I,K),    K=1,KB),I=1,NNode)
   WRITE(1) ((RHO1(I,K),  K=1,KB),I=1,NNode)
   WRITE(1) ((TMEAN1(I,K),K=1,KB),I=1,NNode)
   WRITE(1) ((SMEAN1(I,K),K=1,KB),I=1,NNode)
   WRITE(1) ((RMEAN1(I,K),K=1,KB),I=1,NNode)

   WRITE(1) ((KM(I,K),K=1,KB),I=1,NNode)
   WRITE(1) ((KH(I,K),K=1,KB),I=1,NNode)
   WRITE(1) ((KQ(I,K),K=1,KB),I=1,NNode)

   WRITE(1) (UA(I), I=0,MElem)
   WRITE(1) (VA(I), I=0,MElem)

   WRITE(1) (EL1(I), I=1,MElem)
   WRITE(1) (ET1(I), I=1,MElem)
   WRITE(1) (H1(I),  I=1,MElem)
   WRITE(1) (D1(I),  I=1,MElem)
   WRITE(1) (DT1(I), I=1,MElem)
   WRITE(1) (RTP(I), I=1,MElem)

   WRITE(1) (EL(I), I=1,NNode)
   WRITE(1) (ET(I), I=1,NNode)
   WRITE(1) (H(I),  I=1,NNode)
   WRITE(1) (D(I),  I=1,NNode)
   WRITE(1) (DT(I), I=1,NNode)

#  if defined (EQUI_TIDE)
   WRITE(1) (EL_EQI(I), I=1,NNode)
#  endif
#  if defined (ATMO_TIDE)
   WRITE(1) (EL_ATMO(I), I=1,NNode)
#  endif

#  if defined (WATER_QUALITY)
   DO N1 = 1, NB
     WRITE(1) ((WQM(I,K,N1),K=1,KB),I=1,NNode)
   END DO
#  endif
#  if defined (DYE_RELEASE)
   IF(IINT.GT.IINT_SPE_DYE_B) THEN
   WRITE(1) ((DYE(I,K),K=1,KB),I=1,NNode)
   WRITE(1) ((DYEMEAN(I,K),K=1,KB),I=1,NNode)
   ENDIF
#  endif

   ELSE   
#    if defined (MULTIPROCESSOR)   
     IDUMMY=0
     IF(EL_ASSIM) THEN
       DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         ENKFEL(I)  = STFCT(IDUMMY)
       ENDDO
       DO I=1, NNode
         EL(I)=ENKFEL(NGID(I))
       ENDDO
       DO I=1, MElem
         EL1(I)=(EL(NV(I,1)) + EL(NV(I,2)) + EL(NV(I,3)) )/3.0_DP
       ENDDO
     ENDIF
      IF(UV_ASSIM) THEN
        DO I=1, KBM1
          DO J=1, MElemGL
            IDUMMY = IDUMMY + 1
            ENKFU(J,I) = STFCT(IDUMMY)
          ENDDO
        ENDDO
        DO I=1, KBM1
          DO J=1, MElemGL
            IDUMMY = IDUMMY + 1
            ENKFV(J,I) = STFCT(IDUMMY)
          ENDDO
        ENDDO

        DO I=1, MElem
          DO J=1, KBM1
            U(I,J)=ENKFU(EGID(I),J)
            V(I,J)=ENKFV(EGID(I),J)
          ENDDO
        ENDDO
!        DO I=1, MElem
!          UA(I)=ENKFU(EGID(I),1)
!          VA(I)=ENKFV(EGID(I),1)
!        ENDDO
      ENDIF
      IF(T_ASSIM) THEN
        DO I=1, KBM1
          DO J=1, NNodeGL
            IDUMMY = IDUMMY + 1
            ENKFT(J,I) = STFCT(IDUMMY)
          ENDDO
        ENDDO
        DO I=1, NNode
          DO J=1, KBM1
            T1(I,J)=ENKFT(NGID(I),J)
          ENDDO
        ENDDO
      ENDIF
      IF(S_ASSIM) THEN
        DO I=1, KBM1
          DO J=1, NNodeGL
            IDUMMY = IDUMMY + 1
            ENKFS(J,I) = STFCT(IDUMMY)
          ENDDO
        ENDDO
        DO I=1, NNode
          DO J=1, KBM1
            S1(I,J)=ENKFS(NGID(I),J)
          ENDDO
        ENDDO
      ENDIF

      CALL PWRITE(1,ME,NPC,U,    LBOUND(U,1),    UBOUND(U,1),    MElem,MElemGL,KB,EMAP,0,"U"    )
      CALL PWRITE(1,ME,NPC,V,    LBOUND(V,1),    UBOUND(V,1),    MElem,MElemGL,KB,EMAP,0,"V"    )
      CALL PWRITE(1,ME,NPC,W,    LBOUND(W,1),    UBOUND(W,1),    MElem,MElemGL,KB,EMAP,0,"W"    )
#     if defined (GOTM)
      CALL PWRITE(1,ME,NPC,TKE,  LBOUND(TKE,1),  UBOUND(TKE,1),  NNode,NNodeGL,KB,NMAP,0,"TKE"  )
      CALL PWRITE(1,ME,NPC,TEPS, LBOUND(TEPS,1), UBOUND(TEPS,1), NNode,NNodeGL,KB,NMAP,0,"TEPS" )
#     else
      CALL PWRITE(1,ME,NPC,Q2,   LBOUND(Q2,1),   UBOUND(Q2,1),   NNode,NNodeGL,KB,NMAP,0,"Q2"   )
      CALL PWRITE(1,ME,NPC,Q2L,  LBOUND(Q2L,1),  UBOUND(Q2L,1),  NNode,NNodeGL,KB,NMAP,0,"Q2L"  )
      CALL PWRITE(1,ME,NPC,L,    LBOUND(L,1),    UBOUND(L,1),    NNode,NNodeGL,KB,NMAP,0,"L"  )
#     endif
      CALL PWRITE(1,ME,NPC,S,    LBOUND(S,1),    UBOUND(S,1),    MElem,MElemGL,KB,EMAP,0,"S"    )
      CALL PWRITE(1,ME,NPC,T,    LBOUND(T,1),    UBOUND(T,1),    MElem,MElemGL,KB,EMAP,0,"T"    )
      CALL PWRITE(1,ME,NPC,RHO,  LBOUND(RHO,1),  UBOUND(RHO,1),  MElem,MElemGL,KB,EMAP,0,"RHO"  )
      CALL PWRITE(1,ME,NPC,TMEAN,LBOUND(TMEAN,1),UBOUND(TMEAN,1),MElem,MElemGL,KB,EMAP,0,"TMEAN")
      CALL PWRITE(1,ME,NPC,SMEAN,LBOUND(SMEAN,1),UBOUND(SMEAN,1),MElem,MElemGL,KB,EMAP,0,"SMEAN")
      CALL PWRITE(1,ME,NPC,RMEAN,LBOUND(RMEAN,1),UBOUND(RMEAN,1),MElem,MElemGL,KB,EMAP,0,"RMEAN")

      CALL PWRITE(1,ME,NPC,S1,    LBOUND(S1,1),    UBOUND(S1,1),    NNode,NNodeGL,KB,NMAP,1,"S1"    )
      CALL PWRITE(1,ME,NPC,T1,    LBOUND(T1,1),    UBOUND(T1,1),    NNode,NNodeGL,KB,NMAP,1,"T1"    )
      CALL PWRITE(1,ME,NPC,RHO1,  LBOUND(RHO1,1),  UBOUND(RHO1,1),  NNode,NNodeGL,KB,NMAP,1,"RHO1"  )
      CALL PWRITE(1,ME,NPC,TMEAN1,LBOUND(TMEAN1,1),UBOUND(TMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"TMEAN1")
      CALL PWRITE(1,ME,NPC,SMEAN1,LBOUND(SMEAN1,1),UBOUND(SMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"SMEAN1")
      CALL PWRITE(1,ME,NPC,RMEAN1,LBOUND(RMEAN1,1),UBOUND(RMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"RMEAN1")

      CALL PWRITE(1,ME,NPC,KM,LBOUND(KM,1),UBOUND(KM,1),NNode,NNodeGL,KB,NMAP,1,"KM")
      CALL PWRITE(1,ME,NPC,KH,LBOUND(KH,1),UBOUND(KH,1),NNode,NNodeGL,KB,NMAP,1,"KH")
      CALL PWRITE(1,ME,NPC,KQ,LBOUND(KQ,1),UBOUND(KQ,1),NNode,NNodeGL,KB,NMAP,1,"KQ")

      CALL PWRITE(1,ME,NPC,UA,LBOUND(UA,1),UBOUND(UA,1),MElem,MElemGL,1,EMAP,0,"UA")
      CALL PWRITE(1,ME,NPC,VA,LBOUND(VA,1),UBOUND(VA,1),MElem,MElemGL,1,EMAP,0,"VA")

      CALL PWRITE(1,ME,NPC,EL1,LBOUND(EL1,1),UBOUND(EL1,1),MElem,MElemGL,1,EMAP,1,"EL1")
      CALL PWRITE(1,ME,NPC,ET1,LBOUND(ET1,1),UBOUND(ET1,1),MElem,MElemGL,1,EMAP,1,"ET1")
      CALL PWRITE(1,ME,NPC,H1, LBOUND(H1,1), UBOUND(H1,1), MElem,MElemGL,1,EMAP,1,"H1" )
      CALL PWRITE(1,ME,NPC,D1, LBOUND(D1,1), UBOUND(D1,1), MElem,MElemGL,1,EMAP,1,"D1" )
      CALL PWRITE(1,ME,NPC,DT1,LBOUND(DT1,1),UBOUND(DT1,1),MElem,MElemGL,1,EMAP,1,"DT1")
      CALL PWRITE(1,ME,NPC,RTP,LBOUND(RTP,1),UBOUND(RTP,1),MElem,MElemGL,1,EMAP,1,"RTP")

      CALL PWRITE(1,ME,NPC,EL,LBOUND(EL,1),UBOUND(EL,1),NNode,NNodeGL,1,NMAP,1,"EL")
      CALL PWRITE(1,ME,NPC,ET,LBOUND(ET,1),UBOUND(ET,1),NNode,NNodeGL,1,NMAP,1,"ET")
      CALL PWRITE(1,ME,NPC,H, LBOUND(H,1), UBOUND(H,1), NNode,NNodeGL,1,NMAP,1,"H" )
      CALL PWRITE(1,ME,NPC,D, LBOUND(D,1), UBOUND(D,1), NNode,NNodeGL,1,NMAP,1,"D" )
      CALL PWRITE(1,ME,NPC,DT,LBOUND(DT,1),UBOUND(DT,1),NNode,NNodeGL,1,NMAP,1,"DT")

#     if defined (EQUI_TIDE)
      CALL PWRITE(1,ME,NPC,EL_EQI,LBOUND(EL_EQI,1),UBOUND(EL_EQI,1),NNode,NNodeGL,1,NMAP,1,"EL_EQI")
#     endif
#     if defined (ATMO_TIDE)
      CALL PWRITE(1,ME,NPC,EL_ATMO,LBOUND(EL_ATMO,1),UBOUND(EL_ATMO,1),NNode,NNodeGL,1,NMAP,1,"EL_ATMO")
#     endif

#     if defined (WATER_QUALITY)
      DO N1=1,NB
        CALL PWRITE(1,ME,NPC,WQM(1:NNode,1:KB,N1),LBOUND(WQM(1:NNode,1:KB,N1),1),     &
                      UBOUND(WQM(1:NNode,1:KB,N1),1),NNode,NNodeGL,KB,NMAP,1,"WQM")
      END DO
#     endif
#     if defined (DYE_RELEASE)
     IF(IINT.GT.IINT_SPE_DYE_B) THEN
      CALL PWRITE(1,ME,NPC,DYE,    LBOUND(DYE,1),    UBOUND(DYE,1),    NNode,NNodeGL,KB,NMAP,1,"DYE"    )
      CALL PWRITE(1,ME,NPC,DYEMEAN,LBOUND(DYEMEAN,1),UBOUND(DYEMEAN,1),NNode,NNodeGL,KB,NMAP,1,"DYEMEAN"    )
     ENDIF
#     endif

#     endif

   ENDIF

   IF(MSR) CLOSE(1)

   DEALLOCATE(RTP)
#  if defined (MULTIPROCESSOR)  
     deallocate(enkfel,enkfu,enkfv,enkft,enkfs)
#  endif   
   RETURN
   END SUBROUTINE ST2GR

   SUBROUTINE FCT2ANL(FLNAME,GLNAME)

   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
#  if defined (MULTIPROCESSOR)
   use mod_par
#  endif
   IMPLICIT NONE

   INTEGER I,K,J,N1,ITMP
   INTEGER IDUMMY
   CHARACTER(LEN=120)    :: GLNAME, FLNAME

#  if defined (MULTIPROCESSOR)   
     real(sp), allocatable ::  enkfel(:)   
     real(sp), allocatable ::  enkfu(:,:), enkfv(:,:)
     real(sp), allocatable ::  enkft(:,:), enkfs(:,:)
#  endif 
        

   REAL(SP), ALLOCATABLE :: UTMP(:,:)
   REAL(SP), ALLOCATABLE :: VTMP(:,:)
   REAL(SP), ALLOCATABLE :: WTMP(:,:)
#  if defined (GOTM)
   REAL(SP), ALLOCATABLE :: TKETMP(:,:)
   REAL(SP), ALLOCATABLE :: TEPSTMP(:,:)
#  else   
   REAL(SP), ALLOCATABLE :: Q2TMP(:,:)
   REAL(SP), ALLOCATABLE :: Q2LTMP(:,:)
   REAL(SP), ALLOCATABLE :: LTMP(:,:)
#  endif   
   REAL(SP), ALLOCATABLE :: STMP(:,:)
   REAL(SP), ALLOCATABLE :: TTMP(:,:)
   REAL(SP), ALLOCATABLE :: RHOTMP(:,:)
   REAL(SP), ALLOCATABLE :: TMEANTMP(:,:)
   REAL(SP), ALLOCATABLE :: SMEANTMP(:,:)
   REAL(SP), ALLOCATABLE :: RMEANTMP(:,:)   
   
   REAL(SP), ALLOCATABLE :: S1TMP(:,:)
   REAL(SP), ALLOCATABLE :: T1TMP(:,:)
   REAL(SP), ALLOCATABLE :: RHO1TMP(:,:)
   REAL(SP), ALLOCATABLE :: TMEAN1TMP(:,:)
   REAL(SP), ALLOCATABLE :: SMEAN1TMP(:,:)
   REAL(SP), ALLOCATABLE :: RMEAN1TMP(:,:)    
   REAL(SP), ALLOCATABLE :: KMTMP(:,:)
   REAL(SP), ALLOCATABLE :: KHTMP(:,:)
   REAL(SP), ALLOCATABLE :: KQTMP(:,:)   
   
   REAL(SP), ALLOCATABLE :: UATMP(:)
   REAL(SP), ALLOCATABLE :: VATMP(:)
   REAL(SP), ALLOCATABLE :: EL1TMP(:)
   REAL(SP), ALLOCATABLE :: ET1TMP(:)
   REAL(SP), ALLOCATABLE :: H1TMP(:)
   REAL(SP), ALLOCATABLE :: D1TMP(:)    
   REAL(SP), ALLOCATABLE :: DT1TMP(:)
   REAL(SP), ALLOCATABLE :: RTPTMP(:)   
   
   REAL(SP), ALLOCATABLE :: ELTMP(:)
   REAL(SP), ALLOCATABLE :: ETTMP(:)
   REAL(SP), ALLOCATABLE :: HTMP(:)    
   REAL(SP), ALLOCATABLE :: DTMP(:)
   REAL(SP), ALLOCATABLE :: DTTMP(:)   

#  if defined (EQUI_TIDE)
   REAL(SP), ALLOCATABLE :: EL_EQITMP(:)
#  endif
#  if defined (ATMO_TIDE)
   REAL(SP), ALLOCATABLE :: EL_ATMOTMP(:)
#  endif

#  if defined (WATER_QUALITY)
   REAL(SP), ALLOCATABLE :: WQMTMP(:,:,:)
#  endif

  
#  if defined (MULTIPROCESSOR)
   allocate(enkfel(NNodeGL))        ; enkfel = 0.0_sp
   allocate(enkfu(MElemGL,kbm1))    ; enkfu  = 0.0_sp
   allocate(enkfv(MElemGL,kbm1))    ; enkfv  = 0.0_sp
   allocate(enkft(NNodeGL,kbm1))    ; enkft  = 0.0_sp
   allocate(enkfs(NNodeGL,kbm1))    ; enkfs  = 0.0_sp
#  endif
   
   IF(MSR) THEN

   ALLOCATE(UTMP(0:MElemGL,1:KB))               ; UTMP = 0.0_SP
   ALLOCATE(VTMP(0:MElemGL,1:KB))               ; VTMP = 0.0_SP
   ALLOCATE(WTMP(0:MElemGL,1:KB))               ; WTMP = 0.0_SP  
#  if defined (GOTM)
   ALLOCATE(TKETMP(0:NNodeGL,1:KB))             ; TKETMP  = 0.0_SP
   ALLOCATE(TEPSTMP(0:NNodeGL,1:KB))            ; TEPSTMP = 0.0_SP
#  else   
   ALLOCATE(Q2TMP(0:NNodeGL,1:KB))              ; Q2TMP  = 0.0_SP
   ALLOCATE(Q2LTMP(0:NNodeGL,1:KB))             ; Q2LTMP = 0.0_SP
   ALLOCATE(LTMP(0:NNodeGL,1:KB))               ; LTMP   = 0.0_SP
#  endif   
   ALLOCATE(STMP(0:MElemGL,1:KB))               ; STMP      = 0.0_SP
   ALLOCATE(TTMP(0:MElemGL,1:KB))               ; TTMP      = 0.0_SP
   ALLOCATE(RHOTMP(0:MElemGL,1:KB))             ; RHOTMP    = 0.0_SP
   ALLOCATE(TMEANTMP(0:MElemGL,1:KB))           ; TMEANTMP  = 0.0_SP
   ALLOCATE(SMEANTMP(0:MElemGL,1:KB))           ; SMEANTMP  = 0.0_SP
   ALLOCATE(RMEANTMP(0:MElemGL,1:KB))           ; RMEANTMP  = 0.0_SP

   ALLOCATE(S1TMP(1:NNodeGL,1:KB))              ; S1TMP      = 0.0_SP
   ALLOCATE(T1TMP(1:NNodeGL,1:KB))              ; T1TMP      = 0.0_SP
   ALLOCATE(RHO1TMP(1:NNodeGL,1:KB))            ; RHO1TMP    = 0.0_SP
   ALLOCATE(TMEAN1TMP(1:NNodeGL,1:KB))          ; TMEAN1TMP  = 0.0_SP
   ALLOCATE(SMEAN1TMP(1:NNodeGL,1:KB))          ; SMEAN1TMP  = 0.0_SP
   ALLOCATE(RMEAN1TMP(1:NNodeGL,1:KB))          ; RMEAN1TMP  = 0.0_SP
   ALLOCATE(KMTMP(1:NNodeGL,1:KB))              ; KMTMP      = 0.0_SP
   ALLOCATE(KHTMP(1:NNodeGL,1:KB))              ; KHTMP      = 0.0_SP
   ALLOCATE(KQTMP(1:NNodeGL,1:KB))              ; KQTMP      = 0.0_SP  

   ALLOCATE(UATMP(0:MElemGL))                   ; UATMP      = 0.0_SP
   ALLOCATE(VATMP(0:MElemGL))                   ; VATMP      = 0.0_SP
   ALLOCATE(EL1TMP(1:MElemGL))                  ; EL1TMP     = 0.0_SP
   ALLOCATE(ET1TMP(1:MElemGL))                  ; ET1TMP     = 0.0_SP
   ALLOCATE(H1TMP(1:MElemGL))                   ; H1TMP      = 0.0_SP
   ALLOCATE(D1TMP(1:MElemGL))                   ; D1TMP      = 0.0_SP
   ALLOCATE(DT1TMP(1:MElemGL))                  ; DT1TMP     = 0.0_SP
   ALLOCATE(RTPTMP(1:MElemGL))                  ; RTPTMP     = 0.0_SP

   ALLOCATE(ELTMP(1:NNodeGL))                   ; ELTMP      = 0.0_SP
   ALLOCATE(ETTMP(1:NNodeGL))                   ; ETTMP      = 0.0_SP
   ALLOCATE(HTMP(1:NNodeGL))                    ; HTMP       = 0.0_SP
   ALLOCATE(DTMP(1:NNodeGL))                    ; DTMP       = 0.0_SP
   ALLOCATE(DTTMP(1:NNodeGL))                   ; DTTMP      = 0.0_SP

#  if defined (EQUI_TIDE)
   ALLOCATE(EL_EQITMP(1:NNodeGL))               ; EL_EQITMP  = 0.0_SP
#  endif
#  if defined (ATMO_TIDE)
   ALLOCATE(EL_ATMOTMP(1:NNodeGL))              ; EL_ATMOTMP = 0.0_SP
#  endif

#  if defined (WATER_QUALITY)
   ALLOCATE(WQMTMP(1:NNodeGL,1:KB,1:NB))        ; WQMTMP     = 0.0_SP
#  endif

   ENDIF
   
#  if defined (MULTIPROCESSOR)   
   IDUMMY=0
   IF(EL_ASSIM) THEN
     DO I=1, NNodeGL
       IDUMMY = IDUMMY + 1
       ENKFEL(I)  = STFCT(IDUMMY)
     ENDDO
     DO I=1, NNode
       EL(I)=ENKFEL(NGID(I))
     ENDDO
     DO I=1, MElem
       EL1(I)=(EL(NV(I,1)) + EL(NV(I,2)) + EL(NV(I,3)) )/3.0_DP
     ENDDO
   ENDIF
   IF(UV_ASSIM) THEN
    DO I=1, KBM1
       DO J=1, MElemGL
         IDUMMY = IDUMMY + 1
         ENKFU(J,I) = STFCT(IDUMMY)
       ENDDO
     ENDDO
     DO I=1, KBM1
       DO J=1, MElemGL
         IDUMMY = IDUMMY + 1
         ENKFV(J,I) = STFCT(IDUMMY)
       ENDDO
     ENDDO

     DO I=1, MElem
       DO J=1, KBM1
         U(I,J)=ENKFU(EGID(I),J)
         V(I,J)=ENKFV(EGID(I),J)
       ENDDO
     ENDDO
!     DO I=1, MElem
!       UA(I)=ENKFU(EGID(I),1)
!       VA(I)=ENKFV(EGID(I),1)
!     ENDDO
   ENDIF
   IF(T_ASSIM) THEN
     DO I=1, KBM1
       DO J=1, NNodeGL
         IDUMMY = IDUMMY + 1
         ENKFT(J,I) = STFCT(IDUMMY)
       ENDDO
     ENDDO
     DO I=1, NNode
       DO J=1, KBM1
         T1(I,J)=ENKFT(NGID(I),J)
       ENDDO
     ENDDO
   ENDIF
   IF(S_ASSIM) THEN
     DO I=1, KBM1
       DO J=1, NNodeGL
         IDUMMY = IDUMMY + 1
         ENKFS(J,I) = STFCT(IDUMMY)
       ENDDO
     ENDDO
     DO I=1, NNode
       DO J=1, KBM1
         S1(I,J)=ENKFS(NGID(I),J)
       ENDDO
     ENDDO
   ENDIF
#  endif  
   
   IF(MSR) THEN

   OPEN(INOKF,FILE=TRIM(FLNAME), FORM='UNFORMATTED') 

   REWIND(INOKF)
   READ(INOKF) ITMP
   READ(INOKF) ((UTMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INOKF) ((VTMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INOKF) ((WTMP(I,K),K=1,KB),I=0,MElemGL)
#  if defined (GOTM)
   READ(INOKF) ((TKETMP(I,K),K=1,KB),I=0,NNodeGL)
   READ(INOKF) ((TEPSTMP(I,K),K=1,KB),I=0,NNodeGL)
#  else
   READ(INOKF) ((Q2TMP(I,K),K=1,KB),I=0,NNodeGL)
   READ(INOKF) ((Q2LTMP(I,K),K=1,KB),I=0,NNodeGL)
   READ(INOKF) ((LTMP(I,K),K=1,KB),I=0,NNodeGL)
#  endif
   READ(INOKF) ((STMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INOKF) ((TTMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INOKF) ((RHOTMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INOKF) ((TMEANTMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INOKF) ((SMEANTMP(I,K),K=1,KB),I=0,MElemGL)
   READ(INOKF) ((RMEANTMP(I,K),K=1,KB),I=0,MElemGL)

   READ(INOKF) ((S1TMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INOKF) ((T1TMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INOKF) ((RHO1TMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INOKF) ((TMEAN1TMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INOKF) ((SMEAN1TMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INOKF) ((RMEAN1TMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INOKF) ((KMTMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INOKF) ((KHTMP(I,K),K=1,KB),I=1,NNodeGL)
   READ(INOKF) ((KQTMP(I,K),K=1,KB),I=1,NNodeGL)

   READ(INOKF) (UATMP(I), I=0,MElemGL)
   READ(INOKF) (VATMP(I), I=0,MElemGL)
   READ(INOKF) (EL1TMP(I), I=1,MElemGL)
   READ(INOKF) (ET1TMP(I), I=1,MElemGL)
   READ(INOKF) (H1TMP(I), I=1,MElemGL)
   READ(INOKF) (D1TMP(I), I=1,MElemGL)
   READ(INOKF) (DT1TMP(I), I=1,MElemGL)
   READ(INOKF) (RTPTMP(I), I=1,MElemGL)

   READ(INOKF) (ELTMP(I), I=1,NNodeGL)
   READ(INOKF) (ETTMP(I), I=1,NNodeGL)
   READ(INOKF) (HTMP(I), I=1,NNodeGL)
   READ(INOKF) (DTMP(I), I=1,NNodeGL)
   READ(INOKF) (DTTMP(I), I=1,NNodeGL)

#  if defined (EQUI_TIDE)
   READ(INOKF) (EL_EQITMP(I), I=1,NNodeGL)
#  endif
#  if defined (ATMO_TIDE)
   READ(INOKF) (EL_ATMOTMP(I), I=1,NNodeGL)
#  endif

#  if defined (WATER_QUALITY)
   DO N1=1,NB
     READ(INOKF) ((WQMTMP(I,K,N1),K=1,KB),I=1,NNodeGL)
   END DO
#  endif

   CLOSE(INOKF) 

   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         ELTMP(I) = STFCT(IDUMMY)
      ENDDO
      
!      DTMP = HTMP + ELTMP
!      ETTMP = ELTMP
!      DTTMP = DTMP
      
      DO I=1, MElemGL
         EL1TMP(I)=(ELTMP(NVG(I,1)) + ELTMP(NVG(I,2)) + ELTMP(NVG(I,3)) )/3.0_DP
      ENDDO
   
!      D1TMP = H1TMP + EL1TMP
!      ET1TMP = EL1TMP
!      DT1TMP = D1TMP  
   
   ENDIF

   IF(UV_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MElemGL
          IDUMMY = IDUMMY + 1
          UTMP(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  

      DO K=1, KBM1
        DO I=1, MElemGL
          IDUMMY = IDUMMY + 1
          VTMP(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
     
   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          T1TMP(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          S1TMP(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   OPEN(INOKF,FILE=TRIM(GLNAME), FORM='UNFORMATTED') 
 
   REWIND(INOKF)

   WRITE(INOKF) IEND
   WRITE(INOKF) ((UTMP(I,K),    K=1,KB),I=0,MElemGL)
   WRITE(INOKF) ((VTMP(I,K),    K=1,KB),I=0,MElemGL)
   WRITE(INOKF) ((WTMP(I,K),    K=1,KB),I=0,MElemGL)
#  if defined (GOTM)
   WRITE(INOKF) ((TKETMP(I,K),   K=1,KB),I=0,NNodeGL)
   WRITE(INOKF) ((TEPSTMP(I,K),  K=1,KB),I=0,NNodeGL)
#  else
   WRITE(INOKF) ((Q2TMP(I,K),   K=1,KB),I=0,NNodeGL)
   WRITE(INOKF) ((Q2LTMP(I,K),  K=1,KB),I=0,NNodeGL)
   WRITE(INOKF) ((LTMP(I,K)  ,  K=1,KB),I=0,NNodeGL)
#  endif
   WRITE(INOKF) ((STMP(I,K),    K=1,KB),I=0,MElemGL)
   WRITE(INOKF) ((TTMP(I,K),    K=1,KB),I=0,MElemGL)
   WRITE(INOKF) ((RHOTMP(I,K),  K=1,KB),I=0,MElemGL)
   WRITE(INOKF) ((TMEANTMP(I,K),K=1,KB),I=0,MElemGL)
   WRITE(INOKF) ((SMEANTMP(I,K),K=1,KB),I=0,MElemGL)
   WRITE(INOKF) ((RMEANTMP(I,K),K=1,KB),I=0,MElemGL)

   WRITE(INOKF) ((S1TMP(I,K),    K=1,KB),I=1,NNodeGL)
   WRITE(INOKF) ((T1TMP(I,K),    K=1,KB),I=1,NNodeGL)
   WRITE(INOKF) ((RHO1TMP(I,K),  K=1,KB),I=1,NNodeGL)
   WRITE(INOKF) ((TMEAN1TMP(I,K),K=1,KB),I=1,NNodeGL)
   WRITE(INOKF) ((SMEAN1TMP(I,K),K=1,KB),I=1,NNodeGL)
   WRITE(INOKF) ((RMEAN1TMP(I,K),K=1,KB),I=1,NNodeGL)

   WRITE(INOKF) ((KMTMP(I,K),K=1,KB),I=1,NNodeGL)
   WRITE(INOKF) ((KHTMP(I,K),K=1,KB),I=1,NNodeGL)
   WRITE(INOKF) ((KQTMP(I,K),K=1,KB),I=1,NNodeGL)

   WRITE(INOKF) (UATMP(I), I=0,MElemGL)
   WRITE(INOKF) (VATMP(I), I=0,MElemGL)

   WRITE(INOKF) (EL1TMP(I), I=1,MElemGL)
   WRITE(INOKF) (ET1TMP(I), I=1,MElemGL)
   WRITE(INOKF) (H1TMP(I),  I=1,MElemGL)
   WRITE(INOKF) (D1TMP(I),  I=1,MElemGL)
   WRITE(INOKF) (DT1TMP(I), I=1,MElemGL)
   WRITE(INOKF) (RTPTMP(I), I=1,MElemGL)

   WRITE(INOKF) (ELTMP(I), I=1,NNodeGL)
   WRITE(INOKF) (ETTMP(I), I=1,NNodeGL)
   WRITE(INOKF) (HTMP(I),  I=1,NNodeGL)
   WRITE(INOKF) (DTMP(I),  I=1,NNodeGL)
   WRITE(INOKF) (DTTMP(I), I=1,NNodeGL)

#  if defined (EQUI_TIDE)
   WRITE(INOKF) (EL_EQITMP(I), I=1,NNodeGL)
#  endif
#  if defined (ATMO_TIDE)
   WRITE(INOKF) (EL_ATMOTMP(I), I=1,NNodeGL)
#  endif

#  if defined (WATER_QUALITY)
   DO N1 = 1, NB
     WRITE(INOKF) ((WQMTMP(I,K,N1),K=1,KB),I=1,NNodeGL)
   END DO
#  endif

   CLOSE(INOKF)
   DEALLOCATE(UTMP,VTMP,WTMP,STMP,TTMP,RHOTMP,TMEANTMP,SMEANTMP,RMEANTMP,S1TMP,T1TMP,RHO1TMP)
   DEALLOCATE(TMEAN1TMP,SMEAN1TMP,RMEAN1TMP,KMTMP,KHTMP,KQTMP,UATMP,VATMP,EL1TMP,ET1TMP,H1TMP)
   DEALLOCATE(D1TMP,DT1TMP,RTPTMP,ELTMP,ETTMP,HTMP,DTMP,DTTMP)
#  if defined (EQUI_TIDE)
   DEALLOCATE(EL_EQITMP)
#  endif
#  if defined (ATMO_TIDE)
   DEALLOCATE(EL_ATMOTMP)
#  endif
#  if defined(GOTM)
   DEALLOCATE(TKETMP,TEPSTMP)
#  else
   DEALLOCATE(Q2TMP,Q2LTMP,LTMP)
#  endif    
#  if defined (WATER_QUALITY)
   DEALLOCATE(WQMTMP)
#  endif

   ENDIF

   RETURN
   END SUBROUTINE FCT2ANL

   SUBROUTINE PRINT_ERR(ERR,NDIM,TEXT,ERR2)
   USE MOD_PREC
   IMPLICIT NONE
      
   INTEGER NDIM,K
   REAL(DP) ERR1,ERR2,ERR3,ERR(NDIM)
   CHARACTER*15 TEXT

   ERR1 = 0.0_DP
   ERR2 = 0.0_DP
   ERR3 = 0.0_DP
   DO K = 1, NDIM
     IF (ABS(ERR(K)) > ERR1) THEN
        ERR1 = ABS(ERR(K))
     ENDIF
     ERR2 = ERR2 + ERR(K) * ERR(K)
     ERR3 = ERR3 + ABS(ERR(K))
   ENDDO
   ERR2 = DSQRT(ERR2/NDIM)
   ERR3 = ERR3 / NDIM
!   WRITE (6,'(a15,3(a10,e12.4))') text,' max err=', err1,' rms err=',err2,' mean err=', err3

   RETURN
   END SUBROUTINE PRINT_ERR
   
   SUBROUTINE SET_INI_ENKF
     
   USE LIMS
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE

   INTEGER I,J,K
   INTEGER STDIM
   INTEGER IERR
   REAL(DP) SUM9
   REAL(DP),ALLOCATABLE   ::  STREF(:)
   REAL(DP),ALLOCATABLE   ::  RPETMP(:,:)
   REAL(DP),ALLOCATABLE   ::  RPATMP(:)   
   CHARACTER(LEN=120)     ::  FLNAME, FLNAME2
   CHARACTER(LEN=4)       ::  IFIL
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1

   ALLOCATE(STFCT(STDIM))            ; STFCT   = ZEROD   !!state vector of one ensemble forecast

   IF(MSR) THEN
     ALLOCATE(STREF(STDIM))            ; STREF   = ZEROD  
     ALLOCATE(RPETMP(ENKF_NENS,STDIM)) ; RPETMP  = ZEROD
     ALLOCATE(RPATMP(STDIM))           ; RPATMP  = ZEROD

     DO I=1, ENKF_NENS
       WRITE(IFIL,'(I4.4)') ENKF_START/DELTA_ASS-1-(ENKF_NENS-I)*16
       print *, ifil
       FLNAME = '../restart/re_'//IFIL//'.dat' ! must output restart file at certain time step intervals
       OPEN(INOKF,FILE=FLNAME,FORM='UNFORMATTED')          ! and put them into an "out" folder.
       CALL GR2ST(INOKF)
       
       DO J=1, STDIM
          RPETMP(I,J) = STFCT(J) 
       ENDDO
     ENDDO
     
     WRITE(IFIL,'(I4.4)') (ENKF_START-DELTA_ASS)/DELTA_ASS/ENKF_INT
     FLNAME = '../restart/re_'//IFIL//'.dat'
     OPEN(INOKF,FILE=FLNAME,FORM='UNFORMATTED')
     CALL GR2ST(INOKF)
     STREF = STFCT

     RPATMP = 0.0_DP
     DO K=1, ENKF_NENS     
       DO I=1, STDIM                           
         RPATMP(I) = RPATMP(I) + RPETMP(K,I)/DBLE(ENKF_NENS)
       ENDDO      
     ENDDO  

     ENDIF

     DO K=1, ENKF_NENS
       IF(MSR) THEN                      
         DO I=1, STDIM
           STFCT(I) = RPETMP(K,I)  
         ENDDO
        ENDIF
       WRITE(IFIL,'(I4.4)') K
       FLNAME  = TRIM(OUTDIR)//'/anl/restart'//IFIL//'.dat'
#      if defined(WET_DRY)
       FLNAME2 = TRIM(OUTDIR)//'/anl/restart'//IFIL//'_wd.dat' 
#      endif

#      if defined (MULTIPROCESSOR)
          IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(STFCT,STDIM,MPI_F,0,MPI_COMM_WORLD,IERR)
#      endif

   
       CALL ST2GR(FLNAME)

#    if defined(WET_DRY)

     CALL WET_JUDGE_EL 
     OPEN(INOKF,FILE=FLNAME2,FORM='FORMATTED');
     CALL WD_DUMP_EL(INOKF,I_INITIAL)     

#    endif 
     ENDDO   
   
   DEALLOCATE(STFCT)
   IF(MSR) DEALLOCATE(STREF,RPETMP,RPATMP)

   RETURN 
   END SUBROUTINE SET_INI_ENKF

   SUBROUTINE INI_VALS
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
   IMPLICIT NONE   
   
   U   = 0.0_SP
   V   = 0.0_SP
   W   = 0.0_SP
   UA  = 0.0_SP
   VA  = 0.0_SP
   EL1 = 0.0_SP
   ET1 = 0.0_SP
   EL  = 0.0_SP
   ET  = 0.0_SP
#  if defined (EQUI_TIDE)
   EL_EQI  = 0.0_SP
#  endif
#  if defined (ATMO_TIDE)
   EL_ATMO = 0.0_SP
#  endif
   
   RETURN
   END SUBROUTINE INI_VALS
   
   SUBROUTINE READMEDM
   USE ALL_VARS
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
   WRITE(FILE_NO,'(I4.4)') IEND/DELTA_ASS
   DIR = '../medm_bck/'
   FILENAME=TRIM(CASENAME)//'_sim'//FILE_NO//'.dat'
   OPEN(INOKF,file=TRIM(DIR)//TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED')

   READ(INOKF) IINTT,NGL9,MGL9,THOUR9
!   print*,'medm FILENAME,IINTT,THOUR9=',TRIM(DIR)//TRIM(FILENAME),IINTT,THOUR9

     DO I=1,MElemGL
        READ(INOKF) (UGL(I,K),VGL(I,K),WGL(I,K),KMGL(I,K) , K = 1, KBM1)
     ENDDO

     DO I=1,NNodeGL
        READ(INOKF) ELGL(I),(T1GL(I,K),S1GL(I,K),RHO1GL(I,K),K=1,kbm1)
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
           STTR1(IDUMP) = UGL(I,k) 
	 ENDDO
       ENDDO
       DO K=1, KBM1
	 DO I=1, MElemGL
           IDUMP = IDUMP + 1
           STTR1(IDUMP) = VGL(I,k)
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
     CLOSE(INOKF)
     RETURN

     END SUBROUTINE READMEDM


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
      
!==============================================================================
! PERTUBATE THE B.C. FOR INVERSE METHOD IN NON-JULIAN TIDAL SIMULATION
!==============================================================================
   SUBROUTINE PERT_BC

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER I,J,K,KK,ISEED
   REAL(DP)  RNBC(ENKF_NENS),RNBC_AVG,TMP
   CHARACTER(LEN=4) FNUM 
   CHARACTER(LEN=80)   :: ISTR
   CHARACTER(LEN=80)   :: COMT
   INTEGER,ALLOCATABLE :: NODE_SBC(:)
   REAL(DP),ALLOCATABLE :: APTTMP(:,:),PHAITMP(:,:)
   INTEGER   ISBCN1
   
   
   IF(IBCN_GL(1) > 0)THEN

   IF(S_TYPE == 'non-julian') THEN

     ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)
     OPEN(INOEL,FILE=TRIM(ISTR)//'_el_obc.dat')

     READ(INOEL ,1000) COMT
     READ(INOEL,*) ISBCN1
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
       READ(INOEL,*)  NODE_SBC(I), EMEAN(I)
       READ (INOEL,*) (APT(I,J), J=1,6)
       READ (INOEL,*) (PHAI(I,J), J=1,6)
     END DO

     PHAI = MOD(PHAI,360.0_SP)
     APTTMP  = APT
     PHAITMP = PHAI 

     CLOSE(INOEL)

   ELSE
     WRITE(IPT,*) 'INVERSE METHOD CAN ONLY BE USED FOR NON-JULIAN TIDAL SIMULATION RIGHT NOW!'
     CALL PSTOP   
   ENDIF  
   
   ELSE
     WRITE(IPT,*) 'NO TIDAL B.C.s ARE SPECIFIED TO DO INVERSE, PLEASE CHECK AGAIN!'
     CALL PSTOP
   ENDIF 

   TIMEN = (ENKF_END - ENKF_START)/DELTA_ASS+1
   ALLOCATE(EL_SRS(IBCN_GL(1),TIMEN))    ; EL_SRS   = 0.0_DP
   ALLOCATE(TIME_SER(TIMEN))             ; TIME_SER = 0.0_DP
        
   DO I=1,TIMEN 
     TIME_SER(I)=(DTI*DBLE(ENKF_START+(I-1)*DELTA_ASS)) 
   ENDDO

!  CALL RANDOM_SEED
   RNBC_AVG = 0.0_DP
   TMP = 0.0_DP
   ISEED = -7111
   DO K=1, ENKF_NENS
      RNBC(K) = GASDEV(ISEED)
      RNBC_AVG = RNBC_AVG + RNBC(K)
   ENDDO
   RNBC_AVG = RNBC_AVG/DBLE(ENKF_NENS)
   DO K=1, ENKF_NENS 
      RNBC(K) = RNBC(K) - RNBC_AVG
      TMP = TMP + RNBC(K)
   ENDDO
   WRITE(100,*) TMP

   DO K=1, ENKF_NENS 
      WRITE(FNUM,'(I4.4)') K 

      DO I=1, IBCN_GL(1)
        DO J=1, 6
           APT(I,J)  =  APTTMP(I,J) + BC_AMP_ERR(J)*RNBC(K) 
           PHAI(I,J) = PHAITMP(I,J) + BC_PHA_ERR(J)*RNBC(K)
        ENDDO   
      ENDDO

      PHAI = MOD(PHAI,360.0_SP)

      OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/'//'el_srs'//FNUM//'.dat',STATUS='REPLACE')
      DO I = 1,IBCN_GL(1)
        DO J = 1, TIMEN 
          FORCE = 0.0_SP
          DO KK = 1, 6
            PHAI_IJ= PHAI(I,KK)*PI2/360.0_SP 
            FORCE  = FORCE + APT(I,KK)/100._DP * COS(PI2/PERIOD(KK)*(DTI*FLOAT(ENKF_START+(J-1)*DELTA_ASS)) - PHAI_IJ)
!           FORCE  = FORCE + APT(I,KK)/100._DP * COS(PI2/PERIOD(KK)*TIME_SER(J) - PHAI_IJ)
          ENDDO
          EL_SRS(I,J) = FORCE + EMEAN(I)
        END DO
      END DO
      DO J=1, TIMEN
        WRITE(IOBCKF,'(1000(F13.5))') (EL_SRS(I,J),I=1,IBCN_GL(1))
      ENDDO
      CLOSE(IOBCKF)

      OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/bc_'//FNUM//'.dat')
      DO I=1, IBCN_GL(1)
        WRITE(IOBCKF,'(I10,1000(F13.5))') NODE_SBC(I),EMEAN(I)
        WRITE(IOBCKF,'(1000(F13.5))') ( APT(I,J),J=1,6)
        WRITE(IOBCKF,'(1000(F13.5))') (PHAI(I,J),J=1,6)
      ENDDO
      CLOSE(IOBCKF) 
   ENDDO

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
      
   WRITE(FNUM,'(I4.4)') IENS 
   ISTR = TRIM(OUTDIR)//'/out_err/bc_'//FNUM//'.dat'
   OPEN(IOBCKF,FILE=TRIM(ISTR))
   
!----READ IN BOUNDARY POINTS, AMPLITUDES, AND PHASES OF TIDE-------------------|
!
     DEALLOCATE(EMEAN,APT,PHAI)
     ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
     ALLOCATE(APT(IBCN_GL(1),6), PHAI(IBCN_GL(1),6))
     APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP
     DO I=1,IBCN_GL(1)
       READ(IOBCKF,*)  NODE_SBC(I),EMEAN(I)
       READ (IOBCKF,*) (APT(I,J), J=1,6)
       READ (IOBCKF,*) (PHAI(I,J), J=1,6)
     END DO
   CLOSE(IOBCKF)  

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
   

   SUBROUTINE ENKF_INVERSE(INUM)

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER INUM,I,J,K
   INTEGER NCON
   CHARACTER(LEN=4) FNUM
   INTEGER,ALLOCATABLE :: NODE_SBC(:)
   INTEGER :: NCON_FLG(6)
   
   WRITE(FNUM,'(I4.4)') INUM 

!UPDATE THE ELEVATION TIMESER

   OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/el_srs'//FNUM//'.dat')
   DO J=1,TIMEN 
     READ(IOBCKF,*) (EL_SRS(I,J), I=1,IBCN_GL(1))
   ENDDO
   CLOSE(IOBCKF)

   OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/el_srs'//FNUM//'.dat')
   DO I=1,IBCN_GL(1)
     EL_SRS(I,(ICYC-ENKF_START/DELTA_ASS+1)) = EL_INV(I_OBC_GL(I))             
   ENDDO
   CLOSE(IOBCKF)

   OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/el_srs'//FNUM//'.dat') 
   DO J=1,TIMEN 
     WRITE(IOBCKF,'(1000(F13.5))') (EL_SRS(I,J),I=1,IBCN_GL(1))
   ENDDO
   CLOSE(IOBCKF)

   OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/bc_'//FNUM//'.dat',STATUS='UNKNOWN')
   DEALLOCATE(EMEAN,APT,PHAI)
   ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
   ALLOCATE(APT(IBCN_GL(1),6), PHAI(IBCN_GL(1),6))
   APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP
   DO I=1,IBCN_GL(1)
      READ(IOBCKF,*)  NODE_SBC(I),EMEAN(I)
      READ (IOBCKF,*) 
      READ (IOBCKF,*) 
   ENDDO
   CLOSE(IOBCKF)

!ADJUST AMPLITUDE
     OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/el_srs'//FNUM//'.dat',STATUS='OLD')
     DO J=1,TIMEN
       READ(IOBCKF,*) (EL_SRS(I,J),I=1,IBCN_GL(1))
     ENDDO
     CLOSE(IOBCKF)
     
     NCON     = 0
     NCON_FLG = 0
     DO I=1, 6
        IF(ABS(BC_AMP_ERR(I))>0.0001 .OR. ABS(BC_PHA_ERR(I))>0.0001 ) THEN
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

!     OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/amppha1_'//FNUM//'.dat',STATUS='UNKNOWN',ACCESS='APPEND')
!     WRITE(IOBCKF,*) ICYC, APT(1,2), PHAI(1,2)
!     CLOSE(IOBCKF)

     OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/bc_'//FNUM//'.dat',STATUS='UNKNOWN')
     DO I=1,IBCN_GL(1)
       WRITE(IOBCKF,'(I10,1000(F13.5))')  NODE_SBC(I),EMEAN(I)
       WRITE(IOBCKF,'(1000(F13.5))') (APT(I,J), J=1,6)
       WRITE(IOBCKF,'(1000(F13.5))') (PHAI(I,J), J=1,6)
     ENDDO
     CLOSE(IOBCKF)

!RE-BUILD ELEVATION TIME SERIES
     OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/el_srs'//FNUM//'.dat',STATUS='REPLACE')
     DO I =1, IBCN_GL(1)
       DO J = ICYC-ENKf_START/DELTA_ASS+2, TIMEN
         FORCE = 0.0_SP
         DO K = 1, 6
           PHAI_IJ= PHAI(I,K)*PI2/360.0_SP
           FORCE  = FORCE + APT(I,K)/100.0_DP * COS(PI2/PERIOD(K)*(DTI*FLOAT(ENKF_START+(J-1)*DELTA_ASS)) -PHAI_IJ)
!          FORCE  = FORCE + APT(I,K)/100.0_DP * COS(PI2/PERIOD(K)*TIME_SER(J) -PHAI_IJ)    
         ENDDO
         EL_SRS(I,J)= FORCE + EMEAN(I)       
       ENDDO
     ENDDO
     DO J=1, TIMEN  
       WRITE(IOBCKF,'(1000(F13.5))') (EL_SRS(I,J),I=1,IBCN_GL(1))
     ENDDO
     CLOSE(IOBCKF)

     DEALLOCATE(NODE_SBC)
     RETURN 
   END SUBROUTINE ENKF_INVERSE
   
!=====================================================================================/
!  DETERMINE IF NODES/ELEMENTS ARE WET OR DRY                                         /
!=====================================================================================/
  SUBROUTINE WET_JUDGE_EL

   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined(WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   REAL(SP) :: DTMP
   INTEGER  :: ITA_TEMP
   INTEGER  :: I,IL,IA,IB,K1,K2,K3,K4,K5,K6

#  if defined(WET_DRY)

!
!--Determine If Node Points Are Wet/Dry Based on Depth Threshold---------------!
!
   ISWET_NODE_CURRENTSTEP = 1
   DO I = 1, NNode
     DTMP = H(I) + EL(I)
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWET_NODE_CURRENTSTEP(I) = 0
   END DO

!
!--Determine if Cells are Wet/Dry Based on Depth Threshold---------------------!
!
   ISWET_CELL_CURRENTSTEP = 1
   DO I = 1, MElem
     DTMP =  MAX(EL(NV(I,1)),EL(NV(I,2)),EL(NV(I,3)))  + &
             MIN(  H(NV(I,1)),  H(NV(I,2)),  H(NV(I,3)))
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWET_CELL_CURRENTSTEP(I) = 0
   END DO

!
!--A Secondary Condition for Nodal Dryness-(All Elements Around Node Are Dry)--!
!
   DO I = 1, NNode
     IF(SUM(ISWET_CELL_CURRENTSTEP(NBVE(I,1:NTVE(I)))) == 0)  ISWET_NODE_CURRENTSTEP(I) = 0
   END DO

!
!--Adjust Water Surface So It Does Not Go Below Minimum Depth------------------!
!
   EL = MAX(EL,-H + MIN_DEPTH)

!
!--Recompute Element Based Depths----------------------------------------------!
!
   DO I = 1, MElem
     EL1(I) = ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))
   END DO

!
!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     FWET_N_N = FLOAT(ISWET_NODE_CURRENTSTEP)
     FWET_C_C = FLOAT(ISWET_CELL_CURRENTSTEP)
     CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,FWET_C_C)
     CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,FWET_N_N)
     ISWET_NODE_CURRENTSTEP = INT(FWET_N_N+.5)
     ISWET_CELL_CURRENTSTEP = INT(FWET_C_C+.5)
   END IF
#  endif

#  endif

   RETURN

  END SUBROUTINE WET_JUDGE_EL
  
!==============================================================================|
!   DUMP WET/DRY FLAG DATA FOR RESTART                                         |
!==============================================================================|

   SUBROUTINE WD_DUMP_EL(INF,I_START)

!------------------------------------------------------------------------------|
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined(WET_DRY)
   USE MOD_WD
#  endif

   IMPLICIT NONE
   INTEGER, ALLOCATABLE,DIMENSION(:) :: NTEMP1,NTEMP2
   INTEGER I, INF
   INTEGER I_START
!==============================================================================|

#  if defined(WET_DRY)

   IF(MSR)THEN
     REWIND(INF)
     WRITE(INF,*) I_START
     WRITE(INF,*) MElemGL,NNodeGL
   END IF

   IF(SERIAL)THEN
     WRITE(INF,*) (ISWET_CELL_CURRENTSTEP(I), I=1,MElem)
     WRITE(INF,*) (ISWET_NODE_CURRENTSTEP(I), I=1,NNode)
   ELSE
   ALLOCATE(NTEMP1(MElemGL),NTEMP2(NNodeGL))
#  if defined (MULTIPROCESSOR)
   CALL IGATHER(LBOUND(ISWET_CELL_CURRENTSTEP,1),UBOUND(ISWET_CELL_CURRENTSTEP,1),MElem,MElemGL,1,MYID,NPROCS,EMAP,ISWET_CELL_CURRENTSTEP,NTEMP1)
   CALL IGATHER(LBOUND(ISWET_NODE_CURRENTSTEP,1),UBOUND(ISWET_NODE_CURRENTSTEP,1),NNode,NNodeGL,1,MYID,NPROCS,NMAP,ISWET_NODE_CURRENTSTEP,NTEMP2)
   IF(MSR)THEN
     WRITE(INF,*) (NTEMP1(I), I=1,MElemGL)
     WRITE(INF,*) (NTEMP2(I), I=1,NNodeGL)
   END IF
   DEALLOCATE(NTEMP1,NTEMP2)
#  endif
   END IF

   CLOSE(INF)

#  endif

   RETURN
   END SUBROUTINE WD_DUMP_EL
   
   FUNCTION DISTST(STLOC1,STLOC2)

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   IMPLICIT NONE

   INTEGER IDUMMY, IDUMMY1
   INTEGER RNG(2)
   INTEGER STLOC1,STLOC2
   INTEGER STB,STE,LOC1,LOC2,K
   REAL(DP) BIG_DIST
   REAL*8   DISTST
   
   BIG_DIST = 9000000.0_DP    !1000 Km
   DISTST = BIG_DIST

   IDUMMY  = 0
   IDUMMY1 = 0
   IF(EL_ASSIM) THEN
      RNG(1)  = IDUMMY + 1
      RNG(2)  = IDUMMY + NNodeGL 
      IDUMMY  = NNodeGL
     
      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            DISTST = SQRT((XG(STLOC1)-XG(STLOC2))**2 + (YG(STLOC1)-YG(STLOC2))**2)
            RETURN
         ELSE
            RETURN
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
   ENDIF

   IF(UV_ASSIM) THEN
      RNG(1) = IDUMMY + 1
      RNG(2) = IDUMMY + KBM1*MElemGL
      IDUMMY = IDUMMY + KBM1*MElemGL

      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            IF( INT((STLOC1-IDUMMY1+1)/MElemGL) == INT((STLOC2-IDUMMY1+1)/MElemGL) ) THEN           
               LOC1 = MOD(STLOC1-IDUMMY1+1,MElemGL)
               LOC2 = MOD(STLOC2-IDUMMY1+1,MElemGL)     
               DISTST = SQRT((XCG(LOC1)-XCG(LOC2))**2 + (YCG(LOC1)-YCG(LOC2))**2)
               RETURN
            ELSE
               RETURN
            ENDIF
         ELSE
            RETURN 
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
      
      RNG(1) = IDUMMY + 1
      RNG(2) = IDUMMY + KBM1*MElemGL
      IDUMMY = IDUMMY + KBM1*MElemGL

      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            IF( INT((STLOC1-IDUMMY1+1)/MElemGL) == INT((STLOC2-IDUMMY1+1)/MElemGL) ) THEN           
               LOC1 = MOD(STLOC1-IDUMMY1+1,MElemGL)
               LOC2 = MOD(STLOC2-IDUMMY1+1,MElemGL)     
               DISTST = SQRT((XCG(LOC1)-XCG(LOC2))**2 + (YCG(LOC1)-YCG(LOC2))**2)
               RETURN
            ELSE
               RETURN
            ENDIF
         ELSE
            RETURN 
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
   ENDIF   

   IF(T_ASSIM) THEN
      RNG(1) = IDUMMY + 1
      RNG(2) = IDUMMY + KBM1*NNodeGL
      IDUMMY = IDUMMY + KBM1*NNodeGL
      
      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            IF( INT((STLOC1-IDUMMY1)/NNodeGL) == INT((STLOC2-IDUMMY1)/NNodeGL) ) THEN           
               LOC1 = MOD(STLOC1-IDUMMY1,NNodeGL)
               LOC2 = MOD(STLOC2-IDUMMY1,NNodeGL)     
               DISTST = SQRT((XG(LOC1)-XG(LOC2))**2 + (YG(LOC1)-YG(LOC2))**2)
               RETURN
            ELSE
	       DISTST = 9000000.0_DP
               RETURN
            ENDIF
         ELSE
	    DISTST = 9000000.0_DP
            RETURN 
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
   ENDIF
      
   IF(S_ASSIM) THEN
      RNG(1) = IDUMMY + 1
      RNG(2) = IDUMMY + KBM1*NNodeGL
      IDUMMY = IDUMMY + KBM1*NNodeGL

      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            IF( INT((STLOC1-IDUMMY1)/NNodeGL) == INT((STLOC2-IDUMMY1)/NNodeGL) ) THEN           
               LOC1 = MOD(STLOC1-IDUMMY1,NNodeGL)
               LOC2 = MOD(STLOC2-IDUMMY1,NNodeGL) 	        
               DISTST = SQRT((XG(LOC1)-XG(LOC2))**2 + (YG(LOC1)-YG(LOC2))**2)
	       RETURN
            ELSE
	       DISTST = 9000000.0_DP
               RETURN
            ENDIF
         ELSE
	    DISTST = 9000000.0_DP
            RETURN 
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
   ENDIF  
   
   WRITE(IPT,*) '!ERROR: COULD NOT CALCULATE THE DISTST!'
   CALL PSTOP 
   
   RETURN
   END FUNCTION DISTST
#  endif   
END MODULE MOD_ENKF
