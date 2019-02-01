MODULE MOD_ASSIM
#if defined (NG_OI_ASSIM)

   USE MOD_PREC
#  if defined (WET_DRY)
   USE MOD_WD
#  endif      
   IMPLICIT NONE
   SAVE
!
!--Current Assimilation Object Type
!
   TYPE ASSIM_OBJ_CUR
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS 
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
     INTEGER  :: N_LAYERS                                 !!NUMBER OF OBSERVATIONS IN THE VERTICAL
     INTEGER  :: N_CELL                                   !!CELL NUMBER OF OBSERVATION
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: DEPTH                                    !!DEPTH AT OBSERVATION STATION (MOORING)
     REAL(SP) :: SITA                                     !!LOCAL ISOBATH ANGLE AT OBSERVATION
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ODEPTH      !!OBSERVATION DEPTHS
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TIMES       !!DATA TIMES
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UO,VO       !!OBSERVATION DATA FOR X,Y VELOCITY COMPONENTS
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC 
     INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: S_WEIGHT    !!SIGMA WEIGHTING                               
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
   END TYPE ASSIM_OBJ_CUR

!
!--Temperature/Salinity Assimilation Object Type
!
   TYPE ASSIM_OBJ_TS  
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
     INTEGER  :: N_LAYERS                                 !!NUMBER OF OBSERVATIONS IN THE VERTICAL
     INTEGER  :: N_CELL                                   !!CELL NUMBER OF OBSERVATION
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: DEPTH                                    !!DEPTH AT OBSERVATION STATION (MOORING)
     REAL(SP) :: SITA                                     !!LOCAL ISOBATH ANGLE AT OBSERVATION
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ODEPTH      !!OBSERVATION DEPTHS
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TIMES       !!DATA TIMES
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP        !!OBSERVATION DATA FOR TEMPERATURE 
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: SAL         !!OBSERVATION DATA FOR SALINITY 
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC
     INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: S_WEIGHT    !!SIGMA WEIGHTING
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
   END TYPE ASSIM_OBJ_TS    
                                                                                                                       


!
!--SST Assimilation Object Type
!
   TYPE ASSIM_OBJ_SST
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS 
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR SST OBSERVATION WEIGHTING
     INTEGER  :: N_CELL                                   !!CELL NUMBER OF OBSERVATION
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: T_INT_AVGD                               !!INTERPOLATE SIM DATA AVERAGED OVER OBSERVATION PERIOD
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: T_INT_HOUR  !!INTERPOLATED SIM DATA OBTAINED AT EACH HOUR DURING SIM
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TIMES       !!DATA TIMES (DAYS)
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: SST         !!OBSERVATION DATA FOR SST 
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC 
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
   END TYPE ASSIM_OBJ_SST

     
!
!--Data Assimilation Parameters for SST Assimilation
!
   LOGICAL  :: SST_ASSIM                    !!TRUE IF SST ASSIMILATION ACTIVE
   CHARACTER(LEN=80) :: SST_METHOD           !!'NG' IF SST NUDGING ASSIMILATION ACTIVE
                                            !!'OI' IF SST OI ASSIMILATION ACTIVE
   REAL(SP) :: RAD_SST                      !!SEARCH RADIUS FOR INTERPOLATION POINTS
   REAL(SP) :: GAMA_SST                     
   REAL(SP) :: GALPHA_SST
   REAL(SP) :: ASTIME_WINDOW_SST            !!TIME WINDOW FOR OBSERVATION ASSIMILATION
   INTEGER  :: IAV_DAY
   INTEGER  :: N_INFLU_SST                  !!NUMBER OF INFLUENTIAL OBSERVATIONS FOR OI
   REAL(SP),ALLOCATABLE,DIMENSION(:,:)  :: PARAM_SST

!
!--Data Assimilation Parameters for Current Assimilation
!
   LOGICAL  :: CURRENT_ASSIM                 !!TRUE IF CURRENT ASSIMILATION ACTIVE
   CHARACTER(LEN=80) :: CURRENT_METHOD        !!'NG' IF CURRENT NUDGING ASSIMILATION ACTIVE
                                             !!'OI' IF CURRENT OI ASSIMILATION ACTIVE
   REAL(SP) :: RAD_CUR                       !!SEARCH RADIUS FOR INTERPOLATION POINTS 
   REAL(SP) :: GAMA_CUR
   REAL(SP) :: GALPHA_CUR
   REAL(SP) :: ASTIME_WINDOW_CUR
   INTEGER  :: MAX_LAYER_CUR                 !!MAXIMUM NUMBER OF VERTICAL DATA FROM ANY OBS POINT
   INTEGER  :: N_INFLU_CUR                   !!NUMBER OF INFLUENTIAL OBSERVATIONS FOR OI
   REAL(SP),ALLOCATABLE,DIMENSION(:,:)  :: PARAM_CUR

!
!--Data Assimilation Parameters for Temp/Salinity Data Assimilation
!
   LOGICAL  :: TS_ASSIM                     !!TRUE IF TEMP/SAL ASSIMILATION ACTIVE
   CHARACTER(LEN=80) :: TS_METHOD            !!'NG' IF TEMP/SAL NUDGING ASSIMILATION ACTIVE
                                            !!'OI' IF TEMP/SAL OI ASSIMILATION ACTIVE
   REAL(SP) :: RAD_TS                       !!SEARCH RADIUS FOR INTERPOLATION POINTS
   REAL(SP) :: GAMA_TS  
   REAL(SP) :: GALPHA_TS  
   REAL(SP) :: ASTIME_WINDOW_TS  
   INTEGER  :: MAX_LAYER_TS                 !!MAXIMUM NUMBER OF VERTICAL DATA FROM ANY OBS POINT
   INTEGER  :: N_INFLU_TS                   !!NUMBER OF INFLUENTIAL OBSERVATIONS FOT OI
   REAL(SP),ALLOCATABLE,DIMENSION(:,:)  :: PARAM_TS
   REAL(SP),ALLOCATABLE,DIMENSION(:,:)  :: AW0G,AWXG,AWYG
                                                                                                                       

!
!--Current Data Assimilation Variables
!
   INTEGER                             N_ASSIM_CUR !!NUMBER OF CURRENT OBSERVATIONS 
   TYPE(ASSIM_OBJ_CUR), ALLOCATABLE :: CUR_OBS(:)  !!CURRENT ASSIMILATION DATA OBJECTS
   INTEGER, ALLOCATABLE             :: DA_CUR(:)   !!FLAG IF ELEMENT IS USED FOR CURRENT DA INTERP 

!
!--Salinity/Temperature Data Assimilation Variables
!
   INTEGER                                N_ASSIM_TS   !!NUMBER OF TEMPERATURE OBSERVATIONS
   TYPE(ASSIM_OBJ_TS), ALLOCATABLE     :: TS_OBS(:)    !!TEMP ASSIMILATION DATA OBJECTS
   INTEGER, ALLOCATABLE                :: DA_TS(:)     !!FLAG IF NODE IS USED FOR CURRENT DA INTERP
                                                                                                                       
!
!--SST Data Assimilation Variables
!
   LOGICAL                             PURE_SIM    !!TRUE IF ONLY SIMULATING
   INTEGER                             SST_CYC   !!ITERATION NUMBER FOR SST SWEEPS
   INTEGER                             N_ASSIM_SST !!NUMBER OF SST  OBSERVATIONS 
   INTEGER                             N_TIMES_SST !!NUMBER OF SST OBSERVATIONS TIMES
   INTEGER                             N_DA_HOURS  !!NUMBER OF HOURS IN SIM/DA SWEEP (IAV_DAY*24)
   INTEGER                             DA_HOUR     !!DATA ASSIMALATION HOUR (1-->IAV_DAY*24)               
   INTEGER                             IINT_SST    !!SWEEP NUMBER OF SST ASSIMILATION 
   INTEGER                             INTERNAL_SST!!INNER LOOP OVER SIM/DA CYCLES    
   INTEGER                             ISTART_DAY  !!STARTING DAY FOR SST ASSIM SWEEP
   INTEGER                             IEND_DAY    !!STARTING DAY FOR SST ASSIM SWEEP
   INTEGER                             ASSIM_FLAG  !!TRUE IF ON ASSIMILATION SWEEP   
   INTEGER                             ISWEEP      !!NUMBER OF SWEEPS IN SST ASSIMILATION PROCEDURE
   TYPE(ASSIM_OBJ_SST), ALLOCATABLE :: SST_OBS(:)  !!SST ASSIMILATION DATA OBJECTS
   INTEGER, ALLOCATABLE             :: DA_SST(:)   !!FLAG IF NODE IS USED FOR SST DA INTERP 


   CONTAINS !------------------------------------------------------------------!
            ! SET_ASSIM_PARAM     :   READ ASSIMILATION PARAMETERS FROM INPUT  !
            ! SET_ASSIM_INTERVALS :   SET UP DATA ASSIMILATION SWEEPS          !
            ! SET_CUR_ASSIM_DATA  :   READ AND SET CURRENT ASSIMILATION DATA   ! 
            ! SET_TS_ASSIM_DATA   :   READ AND SET TEMP/SAL ASSIMILATION DATA  ! 
            ! SET_SST_ASSIM_DATA  :   READ AND SET SST ASSIMILATION DATA       ! 
            ! CURRENT_NUDGING     :   NUDGE CURRENT USING ASSIMILATION         !
            ! TEMP NUDGING        :   NUDGE TEMP USING ASSIMILATION            !
            ! SALT_NUDGING        :   NUDGE SALT USING ASSIMILATION            !
            ! SST_NUDGING         :   NUDGE SST USING ASSIMILATION             !
            ! SST_INT             :   INTERPOLATE HOURLY/AVGE SST DATA TO OBS  !
            ! HOT_START_SST       :   READ ICs FOR SST_ASSIMILATION STAGE      !
            ! ARC_SST             :   DUMP ICs FOR SST_ASSIMILATION STAGE      ! 
            ! -----------------------------------------------------------------!

!==============================================================================|
!==============================================================================|

   SUBROUTINE SET_ASSIM_PARAM 

!------------------------------------------------------------------------------|
!  READ IN PARAMETERS CONTROLLING NUDGING OR OI ASSIMILATION                                 |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   USE MOD_INP
   IMPLICIT NONE
   REAL(SP) REALVEC(150)
   INTEGER  INTVEC(150),ISCAN
   CHARACTER(LEN=120) :: FNAME
   INTEGER I


   SST_ASSIM         = .FALSE.
   SST_METHOD        = " "
   RAD_SST           = 0.0_SP
   GAMA_SST          = 0.0_SP
   GALPHA_SST        = 0.0_SP
   ASTIME_WINDOW_SST = 0.0_SP
   IAV_DAY           = 0
   N_INFLU_SST       = 0 

   CURRENT_ASSIM     = .FALSE.
   CURRENT_METHOD    = " "
   RAD_CUR           = 0.0_SP
   GAMA_CUR          = 0.0_SP
   GALPHA_CUR        = 0.0_SP
   ASTIME_WINDOW_CUR = 0.0_SP
   MAX_LAYER_CUR     = 0
   N_INFLU_CUR       = 0

   TS_ASSIM          = .FALSE.
   TS_METHOD         = " "
   RAD_TS            = 0.0_SP
   GAMA_TS           = 0.0_SP
   GALPHA_TS         = 0.0_SP
   ASTIME_WINDOW_TS  = 0.0_SP
   MAX_LAYER_TS      = 0
   N_INFLU_TS        = 0
   
!------------------------------------------------------------------------------|
!   READ IN VARIABLES AND SET VALUES                                           |
!------------------------------------------------------------------------------|

   FNAME = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------|
!   CURRENT ASSIMILATION FLAG
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CURRENT_ASSIM",LVAL = CURRENT_ASSIM)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING CURRENT_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF

   IF(CURRENT_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"CURRENT_METHOD",CVAL = CURRENT_METHOD)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING CURRENT_METHOD: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP 
     END IF
   END IF

!------------------------------------------------------------------------------|
!   TEMP/SALINITY ASSIMILATION FLAG
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"TS_ASSIM",LVAL = TS_ASSIM)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING TS_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   END IF

   IF(TS_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"TS_METHOD",CVAL = TS_METHOD)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING TS_METHOD: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF


!------------------------------------------------------------------------------|
!   SST ASSIMILATION FLAG     
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"SST_ASSIM",LVAL = SST_ASSIM)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING SST_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF

   IF(SST_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"SST_METHOD",CVAL = SST_METHOD)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING SST_METHOD: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP 
     END IF
   END IF

!------------------------------------------------------------------------------|
!   RAD_CUR: CURRENT ASSIMILATION INFLUENCE RADIUS
!------------------------------------------------------------------------------|
   IF(CURRENT_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"RAD_CUR",FSCAL = RAD_CUR)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING RAD_CUR: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP 
     END IF
   END IF

!------------------------------------------------------------------------------|
!   RAD_TS: TEMP/SALINITY ASSIMILATION INFLUENCE RADIUS
!------------------------------------------------------------------------------|
   IF(TS_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"RAD_TS",FSCAL = RAD_TS)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING RAD_TS: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   RAD_SST: SST ASSIMILATION INFLUENCE RADIUS
!------------------------------------------------------------------------------|
   IF(SST_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"RAD_SST",FSCAL = RAD_SST)
        IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING RAD_SST: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP 
     END IF
   END IF

!------------------------------------------------------------------------------|
!   GAMA_CUR: CURRENT ASSIMILATION COEFFICIENT 
!------------------------------------------------------------------------------|
   IF(CURRENT_ASSIM .AND. CURRENT_METHOD == "NG")THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"GAMA_CUR",FSCAL = GAMA_CUR)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING GAMA_CUR: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP 
     END IF
   END IF
   
!------------------------------------------------------------------------------|
!   GAMA_TS: TEMP/SALINITY ASSIMILATION COEFFICIENT
!------------------------------------------------------------------------------|
   IF(TS_ASSIM .AND. TS_METHOD == "NG")THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"GAMA_TS",FSCAL = GAMA_TS)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING GAMA_TS: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   GAMA_SST: SST ASSIMILATION COEFFICIENT 
!------------------------------------------------------------------------------|
   IF(SST_ASSIM .AND. SST_METHOD == "NG")THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"GAMA_SST",FSCAL = GAMA_SST)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING GAMA_SST: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   GALPHA_CUR: CURRENT ASSIMILATION COEFFICIENT 
!------------------------------------------------------------------------------|
   IF(CURRENT_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"GALPHA_CUR",FSCAL = GALPHA_CUR)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING GALPHA_CUR: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   GALPHA_TS: TEMP/SALINITY ASSIMILATION COEFFICIENT
!------------------------------------------------------------------------------|
   IF(TS_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"GALPHA_TS",FSCAL = GALPHA_TS)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING GALPHA_TS: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   GALPHA_SST: SST ASSIMILATION COEFFICIENT 
!------------------------------------------------------------------------------|
   IF(SST_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"GALPHA_SST",FSCAL = GALPHA_SST)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING GALPHA_SST: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF
!------------------------------------------------------------------------------|
!   N_INFLU_CUR: NUMBER OF INFLUENTIAL OBSERVATIONS
!------------------------------------------------------------------------------|
   IF(CURRENT_ASSIM .AND. CURRENT_METHOD == "OI")THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"N_INFLU_CUR",ISCAL = N_INFLU_CUR)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING N_INFLU_CUR: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP 
     END IF
   END IF

!------------------------------------------------------------------------------|
!   N_INFLU_TS: NUMBER OF INFLUENTIAL OBSERVATIONS
!------------------------------------------------------------------------------|
   IF(TS_ASSIM .AND. TS_METHOD == "OI")THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"N_INFLU_TS",ISCAL = N_INFLU_TS)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING N_INFLU_TS: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   N_INFLU_SST: NUMBER OF INFLUENTIAL OBSERVATIONS
!------------------------------------------------------------------------------|
   IF(SST_ASSIM .AND. SST_METHOD == "OI")THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"N_INFLU_SST",ISCAL = N_INFLU_SST)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING N_INFLU_SST: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   ASTIME_WINDOW_CUR: TIME WINDOW FOR CURRENT OBSERVATION INFLUENCE
!------------------------------------------------------------------------------|
   IF(CURRENT_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"ASTIME_WINDOW_CUR",FSCAL = ASTIME_WINDOW_CUR)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING ASTIME_WINDOW_CUR: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   ASTIME_WINDOW_TS: TIME WINDOW FOR TEMP/SALINITY OBSERVATION INFLUENCE
!------------------------------------------------------------------------------|
   IF(TS_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"ASTIME_WINDOW_TS",FSCAL = ASTIME_WINDOW_TS)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING ASTIME_WINDOW_TS: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   ASTIME_WINDOW_SST: TIME WINDOW FOR SST OBSERVATION INFLUENCE
!------------------------------------------------------------------------------|
   IF(SST_ASSIM)THEN
     ISCAN = SCAN_FILE(TRIM(FNAME),"ASTIME_WINDOW_SST",FSCAL = ASTIME_WINDOW_SST)
     IF(ISCAN /= 0)THEN
       WRITE(IPT,*)'ERROR READING ASTIME_WINDOW_SST: ',ISCAN
       IF(ISCAN == -2)THEN
         WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       END IF
       CALL PSTOP
     END IF
   END IF

!------------------------------------------------------------------------------|
!   IAV_DAY: NUMBER OF DAYS OVER WHICH TO AVERAGE SST DATA      
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"IAV_DAY",ISCAL = IAV_DAY)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IAV_DAY: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   END IF
   N_DA_HOURS = IAV_DAY*24

!------------------------------------------------------------------------------|
!   MODIFY VARIABLES TO CORRESPOND TO CORRECT UNITS              
!------------------------------------------------------------------------------|
   ASTIME_WINDOW_CUR = ASTIME_WINDOW_CUR*3600.   !!CONVERT HOURS --> SECONDS
   ASTIME_WINDOW_TS  = ASTIME_WINDOW_TS*3600.    !!CONVERT HOURS --> SECONDS
   ASTIME_WINDOW_SST = ASTIME_WINDOW_SST*3600.   !!CONVERT HOURS --> SECONDS

!------------------------------------------------------------------------------|
!            SCREEN REPORT OF SET VARIABlES                                    !
!------------------------------------------------------------------------------|
   IF(MSR)THEN
     WRITE(IPT,*)''
     WRITE(IPT,*)'!        DATA ASSIMILATION PARAMETERS       '
     IF(CURRENT_ASSIM)THEN
       WRITE(IPT,*)'!  # CURRENT_ASSIM       :  ACTIVE'
       WRITE(IPT,*)'!  # CURRENT_METHOD      :',TRIM(CURRENT_METHOD)
       WRITE(IPT,*)'!  # RAD_CUR             :',RAD_CUR
       WRITE(IPT,*)'!  # GALPHA_CUR          :',GALPHA_CUR
       IF(CURRENT_METHOD == 'NG')WRITE(IPT,*)'!  # GAMA_CUR            :',GAMA_CUR   
       IF(CURRENT_METHOD == 'OI')WRITE(IPT,*)'!  # N_INFLU_CUR         :',N_INFLU_CUR
       WRITE(IPT,*)'!  # ASTIME_WINDOW_CUR   :',ASTIME_WINDOW_CUR
     ELSE
       WRITE(IPT,*)'!  # CURRENT_ASSIM       :  NOT ACTIVE'
     END IF
     IF(SST_ASSIM)THEN
       WRITE(IPT,*)'!  # SST_ASSIM           :  ACTIVE'
       WRITE(IPT,*)'!  # SST_METHOD          :',TRIM(SST_METHOD)
       WRITE(IPT,*)'!  # RAD_SST             :',RAD_SST
       WRITE(IPT,*)'!  # GALPHA_SST          :',GALPHA_SST
       IF(SST_METHOD == 'NG')WRITE(IPT,*)'!  # GAMA_SST            :',GAMA_SST   
       IF(SST_METHOD == 'OI')WRITE(IPT,*)'!  # N_INFLU_SST         :',N_INFLU_SST
       WRITE(IPT,*)'!  # ASTIME_WINDOW_SST   :',ASTIME_WINDOW_SST
       WRITE(IPT,*)'!  # IAV_DAY             :',IAV_DAY          
       WRITE(IPT,*)'!  # N_DA_HOURS          :',N_DA_HOURS       
     ELSE
       WRITE(IPT,*)'!  # SST_ASSIM           :  NOT ACTIVE'
     END IF
     IF(TS_ASSIM)THEN
       WRITE(IPT,*)'!  # TS_ASSIM            :  ACTIVE'
       WRITE(IPT,*)'!  # TS_METHOD           :',TRIM(TS_METHOD)
       WRITE(IPT,*)'!  # RAD_TS              :',RAD_TS 
       WRITE(IPT,*)'!  # GALPHA_TS           :',GALPHA_TS 
       IF(TS_METHOD == 'NG')WRITE(IPT,*)'!  # GAMA_TS             :',GAMA_TS 
       IF(TS_METHOD == 'OI')WRITE(IPT,*)'!  # N_INFLU_TS          :',N_INFLU_TS
       WRITE(IPT,*)'!  # ASTIME_WINDOW_TS    :',ASTIME_WINDOW_TS 
     ELSE
       WRITE(IPT,*)'!  # TS_ASSIM            :  NOT ACTIVE'
     END IF
   END IF


   RETURN
   END SUBROUTINE SET_ASSIM_PARAM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
   SUBROUTINE SET_ASSIM_INTERVALS
!==============================================================================!
                                                                                                                            
!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION SWEEP INTERVALS                                         |
!------------------------------------------------------------------------------!
   USE ALL_VARS
   IMPLICIT NONE

   PURE_SIM   = .TRUE.
   ISTART_DAY = 1
   IEND_DAY   = 1
   ISWEEP     = 1
   IF(SST_ASSIM)THEN
     ISTART_DAY=INT(FLOAT(ISTART-1)*DTI/86400.+0.5)
     IEND_DAY  =INT(FLOAT(IEND)*DTI/86400.-1.+0.5)
     PURE_SIM  = .FALSE.
     ISWEEP    = 2
     IF(MSR)THEN
       WRITE(IPT,*)'!'
       WRITE(IPT,*)'!  STARTING DAY :  ',ISTART_DAY
       WRITE(IPT,*)'!  ENDING DAY   :  ',IEND_DAY
       WRITE(IPT,*)'!  SST DA INT   :  ',IAV_DAY
       WRITE(IPT,*)'!  #SST DA INT  :  ',(IEND_DAY-ISTART_DAY-1)/IAV_DAY
     END IF
   END IF
   RETURN
   END SUBROUTINE SET_ASSIM_INTERVALS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
   SUBROUTINE SET_CUR_ASSIM_DATA 
!==============================================================================!

!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR CURRENT OBSERVATIONS                           |
!------------------------------------------------------------------------------!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)   
   USE MOD_PAR  
# endif
   IMPLICIT NONE
   INTEGER I,J,K,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   CHARACTER(LEN= 2 ) :: NAC   
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
   REAL(SP), PARAMETER :: ALF = 0.05_SP
   LOGICAL :: FEXIST
   INTEGER :: MAXEL,NBD_CNT
   REAL(SP) :: LMIN 
   INTEGER :: JMIN,JJ
       
   REAL(SP), DIMENSION(1:MElemGL,1) :: RDLIST
   REAL(SP), DIMENSION(3) :: XTRI,YTRI
   REAL(SP) :: RDLAST
   INTEGER :: LOCIJ(2),MIN_LOC,IERR,Nsite_tmp
   INTEGER :: ND1,ND2,ND3
   REAL(SP) :: DELTA,COFA,COFB,COFC
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP
   REAL(SP), DIMENSION(KB) :: ZZ_OB
#  if defined (MULTIPROCESSOR)   
   REAL(SP), ALLOCATABLE :: ZZ_G(:,:)  
# endif
       
!------------------------------------------------------------------------------!
!  Read Number of Current Observations and Coordinates of Each                 !
!------------------------------------------------------------------------------!
       
   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_current.xy"
!
!--Make Sure Current Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number of Current Measurement Stations---------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_CUR
   ALLOCATE(CUR_OBS(N_ASSIM_CUR))

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!

   DO I=1,N_ASSIM_CUR
     READ(1,*)ITMP,CUR_OBS(I)%X,CUR_OBS(I)%Y,CUR_OBS(I)%DEPTH,NLAY,CUR_OBS(I)%SITA
     CUR_OBS(I)%N_LAYERS = NLAY
     ALLOCATE(CUR_OBS(I)%ODEPTH(NLAY))
     DO J=1,NLAY
       READ(1,*)CUR_OBS(I)%ODEPTH(J)
       IF(CUR_OBS(I)%ODEPTH(J) > CUR_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'OBSERVATION DEPTH',J,'OF CURRENT MOORING',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH'
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO
   MAX_LAYER_CUR = MAXVAL(CUR_OBS(1:N_ASSIM_CUR)%N_LAYERS)

!
!--Shift Coordinates-----------------------------------------------------------!
!
   CUR_OBS(:)%X = CUR_OBS(:)%X - VXMIN 
   CUR_OBS(:)%Y = CUR_OBS(:)%Y - VYMIN 
   
   IF(CURRENT_METHOD == 'OI')THEN
!
!--find the cell number (TS_OBS(:)%N_CELL) of Obs station---------------------
!
   RRTMP = 100000.0     !100km
   DO J= 1,N_ASSIM_CUR 
      X0 = CUR_OBS(J)%X
      Y0 = CUR_OBS(J)%Y
      DO I=1,MElemGL
         Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
	 if(Rtmp.LT.RRTMP) then
	    S11 = (XG(NVG(I,2))-X0)*(YG(NVG(I,3))-Y0)-&
	          (XG(NVG(I,3))-X0)*(YG(NVG(I,2))-Y0)
            S22 = (XG(NVG(I,3))-X0)*(YG(NVG(I,1))-Y0)-&
	          (XG(NVG(I,1))-X0)*(YG(NVG(I,3))-Y0)
            S33 = (XG(NVG(I,1))-X0)*(YG(NVG(I,2))-Y0)-&
	          (XG(NVG(I,2))-X0)*(YG(NVG(I,1))-Y0)
	    IF(S11.LE.0.AND.S22.LE.0.AND.S33.LE.0) THEN
	      CUR_OBS(J)%N_CELL = I
	      GOTO 300
	    ELSE
	      CUR_OBS(J)%N_CELL = 0
	    ENDIF
	 ELSE
	    CUR_OBS(J)%N_CELL = -1
	 ENDIF
      ENDDO	
300 CONTINUE
      IF(CUR_OBS(J)%N_CELL.LE.0) THEN
         IF(MSR) WRITE(IPT,*)'ERROR--CURRENT OBS SITE:',J,' OUT OF DOMAN',&
	         CUR_OBS(J)%N_CELL 
         CALL PSTOP
      ENDIF 	 

   ENDDO             
!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   ALLOCATE(AW0G(MElemGL,3))
   ALLOCATE(AWXG(MElemGL,3))
   ALLOCATE(AWYG(MElemGL,3))
   
   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF

#  if defined (MULTIPROCESSOR)     
   IF(PAR)THEN
     CALL GATHER(LBOUND(AW0,1),UBOUND(AW0,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AW0,AW0G)
     CALL GATHER(LBOUND(AWX,1),UBOUND(AWX,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AWX,AWXG)
     CALL GATHER(LBOUND(AWY,1),UBOUND(AWY,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AWY,AWYG)
   END IF
#  endif     
   END IF
!
!--Close Current Observation Global File---------------------------------------!
!
   CLOSE(1)
       
   IF(CURRENT_METHOD == 'OI')THEN
!------------------------------------------------------------------------------!
!  Read Correlation Length of Current Observations                             !
!------------------------------------------------------------------------------!
       
!JQI   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_radius_cur.dat"

!JQI   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
!JQI   IF(MSR .AND. .NOT.FEXIST)THEN
!JQI     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
!JQI     WRITE(IPT,*)'HALTING.....'
!JQI     CALL PSTOP
!JQI   END IF
     
!JQI   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')

   ALLOCATE(PARAM_CUR(2,N_ASSIM_CUR))

!JQI   DO I=1,N_ASSIM_CUR
!JQI     READ(1,*)PARAM_CUR(1,I),PARAM_CUR(2,I)
!JQI   END DO

!JQI   CLOSE(1)
   PARAM_CUR = 30000.0_SP
   END IF
       
!------------------------------------------------------------------------------!
!  Open Current Observation Files for Each Observation Point and Read Data     !
!------------------------------------------------------------------------------!
!----Make Sure Current Observation File Exists--------------------!
   ONAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//'_cur.dat'
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF

!----Open Current Observation File for Read------------------------------!
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)


   DO I=1,N_ASSIM_CUR
     READ(1,*,IOSTAT=IOS) nsite_tmp,CUR_OBS(I)%N_TIMES
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_cur.dat at site number:',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Current (UA,VA) and Time (TIME)-------------------!
     ALLOCATE(CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES))
     ALLOCATE(CUR_OBS(I)%UO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))
     ALLOCATE(CUR_OBS(I)%VO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))

!----Read in Current Data for Observation I------------------------------------!
     NLAY = CUR_OBS(I)%N_LAYERS
     DO J=1,CUR_OBS(I)%N_TIMES
!       READ(1,*)CUR_OBS(I)%TIMES(J),(CUR_OBS(I)%UO(J,K),CUR_OBS(I)%VO(J,K),K=1,NLAY)
       READ(1,*,IOSTAT=IOS)CUR_OBS(I)%TIMES(J),   &
                          (CUR_OBS(I)%UO(J,K),CUR_OBS(I)%VO(J,K),K=1,NLAY)
			  
       IF(IOS < 0)THEN         ! ios=0 if all goes ok.
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_cur.dat at site:',I,  &
                     'Time No:',J
         CALL PSTOP
       END IF  
     END DO

!----Convert Time to Seconds---------------------------------------------------!
!----Shift Jan 1 Based Time Data to Dec 1 Based Time Data-----CASESPECIFIC-----!
     IF(trim(CASENAME) == 'gom')THEN
       CUR_OBS(I)%TIMES = ((CUR_OBS(I)%TIMES-1.0_SP)*24.0_SP+744.0_SP)*3600.0_SP
     ELSE   
       CUR_OBS(I)%TIMES = CUR_OBS(I)%TIMES*3600.0_SP*24.0_SP
     END IF
 
!----Convert Current Data from cm/s to m/s-------------------------------------!
     CUR_OBS(I)%UO = CUR_OBS(I)%UO * .01_SP
     CUR_OBS(I)%VO = CUR_OBS(I)%VO * .01_SP
  END DO
  CLOSE(1)

!------------------------------------------------------------------------------!
!  Count Number of Points with Bad Data (UO = 0. + V0 = 0.)         
!------------------------------------------------------------------------------!
  NBD_CNT = 0
  DO I=1,N_ASSIM_CUR
    DO J=1,CUR_OBS(I)%N_TIMES
      DO K=1,CUR_OBS(I)%N_LAYERS
        IF(ABS(CUR_OBS(I)%UO(J,K))+ABS(CUR_OBS(I)%VO(J,K)) < .0001)THEN
          NBD_CNT = NBD_CNT + 1
        END IF
      END DO
    END DO
  END DO
 
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Mooring Location 
!------------------------------------------------------------------------------!
!   DO I=1,N_ASSIM_CUR
!   LMIN = 100000000.
!     X0 = CUR_OBS(I)%X
!     Y0 = CUR_OBS(I)%Y
!     DO J=1,NNodeGL
!       DX = ABS(XG(J)-X0)
!       DY = ABS(YG(J)-Y0)
!       IF(SQRT(DX**2 + DY**2) < LMIN)THEN
!         LMIN = SQRT(DX**2 + DY**2)
!         JMIN = J
!       END IF
!      END DO
!      CUR_OBS(I)%SITA = SITA_GD(JMIN) + 3.14159_SP/2.0_SP
!    END DO


   ALLOCATE(ITEMP(MElemGL),FTEMP(MElemGL),DA_CUR(MElemGL))     ; DA_CUR = 0
   DO I=1,N_ASSIM_CUR
     X0 = CUR_OBS(I)%X
     Y0 = CUR_OBS(I)%Y
     ISOBATH_ANGLE = CUR_OBS(I)%SITA/180.0_SP*3.1415926_SP
     ECNT = 0
     DO J=1,MElemGL
       DX = ABS(XCG(J)-X0)
       DY = ABS(YCG(J)-Y0)
       RD = SQRT(DX**2 + DY**2)
       IF(RD <= RAD_CUR)THEN
         DA_CUR(J)   = 1
         ECNT        = ECNT + 1      
         ITEMP(ECNT) =  J
         FTEMP(ECNT) = (RAD_CUR**2 - RD**2) / (RAD_CUR**2 + RD**2)
         ANG_OBS_SIM = ATAN2(DY,DX)
         D_ANGLE     = ANG_OBS_SIM - ISOBATH_ANGLE 
         D_ANGLE     = D_ANGLE - INT(D_ANGLE/3.1415926_SP)*3.1415926_SP
         D_ANGLE     = ABS(D_ANGLE)
         DIR_WEIGHT  = (ABS(D_ANGLE-0.5*3.1415926_SP)+ALF*3.1415926_SP)/ &
                       ((0.5_SP+ALF)*3.1415926_SP)
         FTEMP(ECNT) = FTEMP(ECNT)*DIR_WEIGHT
       END IF
     END DO
     IF(ECNT == 0)THEN
       WRITE(IPT,*)'ERROR SETTING UP CURRENT DATA ASSIMILATION'
       WRITE(IPT,*)'NO ELEMENTS LIE WITHIN RADIUS',RAD_CUR
       WRITE(IPT,*)'OF OBSERVATION POINT',I
       CALL PSTOP   
     ELSE
       CUR_OBS(I)%N_INTPTS = ECNT
       ALLOCATE(CUR_OBS(I)%INTPTS(ECNT))
       ALLOCATE(CUR_OBS(I)%X_WEIGHT(ECNT))
       CUR_OBS(I)%INTPTS(1:ECNT)  = ITEMP(1:ECNT)
       CUR_OBS(I)%X_WEIGHT(1:ECNT) = FTEMP(1:ECNT)
     END IF
   END DO
   DEALLOCATE(FTEMP,ITEMP)


     
!------------------------------------------------------------------------------!
!  Compute Sigma Layer Weights for Vertical Interpolation                                                                 
!------------------------------------------------------------------------------!
#  if defined (MULTIPROCESSOR)   
   ALLOCATE(ZZ_G(0:NNodeGL,KB))  
   IF(PAR)CALL GATHER(LBOUND(ZZ,1),UBOUND(ZZ,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,ZZ,ZZ_G)
   IF(PAR)CALL MPI_BCAST(ZZ_G,NNodeGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
#  endif


   DO I=1,N_ASSIM_CUR
     NLAY = CUR_OBS(I)%N_LAYERS
     ALLOCATE(CUR_OBS(I)%S_INT(NLAY,2))
     ALLOCATE(CUR_OBS(I)%S_WEIGHT(NLAY,2))

     X0 = CUR_OBS(I)%X
     Y0 = CUR_OBS(I)%Y
     RDLIST(1:MElemGL,1) = SQRT((XCG(1:MElemGL)-X0)**2 + (YCG(1:MElemGL)-Y0)**2)
     RDLAST = -1.0_SP
in:  DO WHILE(.TRUE.)
       LOCIJ = MINLOC(RDLIST,RDLIST>RDLAST)
       MIN_LOC = LOCIJ(1)
       IF(MIN_LOC == 0)THEN
         EXIT in
       END IF
       XTRI = XG(NVG(MIN_LOC,1:3))
       YTRI = YG(NVG(MIN_LOC,1:3))
       RDLAST = RDLIST(MIN_LOC,1)
       IF(ISINTRIANGLE1(XTRI,YTRI,X0,Y0))THEN
         JJ = MIN_LOC
	 EXIT IN
       END IF
       RDLAST = RDLIST(MIN_LOC,1)
     END DO IN  	 	      

     ND1=NVG(JJ,1)
     ND2=NVG(JJ,2)
     ND3=NVG(JJ,3)
     DELTA=(XG(ND2)-XG(ND1))*(YG(ND3)-YG(ND1))-     &
           (XG(ND3)-XG(ND1))*(YG(ND2)-YG(ND1))

     IF(SERIAL)THEN
     DO K=1,KBM1
       COFA=(YG(ND3)-YG(ND1))*(ZZ(ND2,K)-ZZ(ND1,K))-   &
            (YG(ND2)-YG(ND1))*(ZZ(ND3,K)-ZZ(ND1,K))
       COFB=(XG(ND2)-XG(ND1))*(ZZ(ND3,K)-ZZ(ND1,K))-   &
            (XG(ND3)-XG(ND1))*(ZZ(ND2,K)-ZZ(ND1,K))
       COFA=COFA/DELTA
       COFB=COFB/DELTA
       COFC=ZZ(ND1,K)-COFA*XG(ND1)-COFB*YG(ND1)
       ZZ_OB(K)=COFA*X0+COFB*Y0+COFC
     END DO  
     END IF

#  if defined (MULTIPROCESSOR)   
     IF(PAR)THEN
     DO K=1,KBM1
       COFA=(YG(ND3)-YG(ND1))*(ZZ_G(ND2,K)-ZZ_G(ND1,K))-   &
            (YG(ND2)-YG(ND1))*(ZZ_G(ND3,K)-ZZ_G(ND1,K))
       COFB=(XG(ND2)-XG(ND1))*(ZZ_G(ND3,K)-ZZ_G(ND1,K))-   &
            (XG(ND3)-XG(ND1))*(ZZ_G(ND2,K)-ZZ_G(ND1,K))
       COFA=COFA/DELTA
       COFB=COFB/DELTA
       COFC=ZZ_G(ND1,K)-COFA*XG(ND1)-COFB*YG(ND1)
       ZZ_OB(K)=COFA*X0+COFB*Y0+COFC
     END DO  
     END IF
#  endif
     
     DO J=1,NLAY
       SIGMA_C = -CUR_OBS(I)%ODEPTH(J)/CUR_OBS(I)%DEPTH
       DO K=2,KBM1
         IF(ZZ_OB(K) <= SIGMA_C .AND. ZZ_OB(K-1) > SIGMA_C)THEN 
           CUR_OBS(I)%S_INT(J,1) = K-1
           CUR_OBS(I)%S_INT(J,2) = K
           CUR_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ_OB(K))/(ZZ_OB(K-1)-ZZ_OB(K))
           CUR_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - CUR_OBS(I)%S_WEIGHT(J,1) 
         END IF  
       END DO
       IF(ZZ_OB(1) < SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
         CUR_OBS(I)%S_INT(J,1) = 1
         CUR_OBS(I)%S_INT(J,2) = 1
         CUR_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         CUR_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       IF(ZZ_OB(KBM1) > SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
         CUR_OBS(I)%S_INT(J,1) = KBM1
         CUR_OBS(I)%S_INT(J,2) = KBM1
         CUR_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         CUR_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF

     END DO
   END DO
#  if defined (MULTIPROCESSOR)   
   DEALLOCATE(ZZ_G)  
#  endif


!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN

   WRITE(IPT,*)
   WRITE(IPT,*)'!            CURRENT OBSERVATION DATA           '
   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_EL   SITA"
   DO I=1,N_ASSIM_CUR
     MAXEL = MAXLOC(CUR_OBS(I)%X_WEIGHT,DIM=1)
     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
     I,CUR_OBS(I)%X/1000.,CUR_OBS(I)%Y/1000., &
       CUR_OBS(I)%N_INTPTS,CUR_OBS(I)%N_TIMES,CUR_OBS(I)%INTPTS(MAXEL),&
       CUR_OBS(I)%SITA
   END DO
   WRITE(IPT,*)
   WRITE(IPT,*)'NUMBER OF BAD CURRENT DATA POINTS: ',NBD_CNT
   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
   DO I=1,N_ASSIM_CUR
   WRITE(IPT,*)I,CUR_OBS(I)%TIMES(1)/(24.*3600.),&
       CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES)/(24.*3600.)
   END DO

   RETURN
   END SUBROUTINE SET_CUR_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE CURRENT_ASSIMILATION
!==============================================================================|
!  USE CURRENT OBSERVATION DATA TO ADJUST VELOCITY COMPONENTS                  |
!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR 
#  endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UINT,VINT,UCORR,VCORR,UG,VG,TWGHT
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UCORR1,VCORR1
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
   INTEGER I,J,K,J1,K1,K2,NLAY,ITIME,NTIME,IERR
   INTRINSIC MINLOC

   INTEGER IP,JN1,JN2,JN3
   REAL(SP) XSC,YSC,COFT0,COFTX,COFTY
!==============================================================================|

       
!------------------------------------------------------------------------------!
!  Gather U and V Fields to Master Processor                                   ! 
!------------------------------------------------------------------------------!
   ALLOCATE(UG(MElemGL,KB))
   ALLOCATE(VG(MElemGL,KB))
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     CALL GATHER(LBOUND(UF,1),  UBOUND(UF,1),  MElem,MElemGL,KB,MYID,NPROCS,EMAP,UF,UG)
     CALL GATHER(LBOUND(VF,1),  UBOUND(VF,1),  MElem,MElemGL,KB,MYID,NPROCS,EMAP,VF,VG)
   END IF
#  endif
   IF(SERIAL)THEN
     UG(1:MElemGL,1:KBM1) = UF(1:MElemGL,1:KBM1)
     VG(1:MElemGL,1:KBM1) = VF(1:MElemGL,1:KBM1)
   END IF
!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!

   IF(MSR)THEN
   CUR_OBS%T_WEIGHT = 0. 
   T_THRESH         = ASTIME_WINDOW_CUR   
   SIMTIME          = TIME*86400

   DO I=1,N_ASSIM_CUR       
     NTIME = CUR_OBS(I)%N_TIMES
     ALLOCATE(FTEMP(NTIME)) 
     FTEMP(1:NTIME) = ABS(SIMTIME - CUR_OBS(I)%TIMES(1:NTIME))
     DT_MIN = MINVAL(FTEMP(1:NTIME))
     CUR_OBS(I)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)

     IF(DT_MIN < T_THRESH)THEN     
       IF(DT_MIN < .5_SP*T_THRESH) THEN
         CUR_OBS(I)%T_WEIGHT = 1.0_SP
       ELSE
         CUR_OBS(I)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
       END IF
     END IF

     DEALLOCATE(FTEMP)
   END DO
   
       
!------------------------------------------------------------------------------!
!  Interpolate Simulation Data to Local Observation Point                      ! 
!------------------------------------------------------------------------------!
       
   ALLOCATE(UINT(N_ASSIM_CUR,MAX_LAYER_CUR)) ; UINT = 0. 
   ALLOCATE(VINT(N_ASSIM_CUR,MAX_LAYER_CUR)) ; VINT = 0.

   IF(CURRENT_METHOD == 'NG')THEN
     DO I=1,N_ASSIM_CUR  
       DO J=1,CUR_OBS(I)%N_INTPTS
         J1        = CUR_OBS(I)%INTPTS(J)
         WGHT      = CUR_OBS(I)%X_WEIGHT(J)
         NLAY      = CUR_OBS(I)%N_LAYERS
         DO K=1,NLAY
           U1 = UG(J1,CUR_OBS(I)%S_INT(K,1))
           U2 = UG(J1,CUR_OBS(I)%S_INT(K,2))
           V1 = VG(J1,CUR_OBS(I)%S_INT(K,1))
           V2 = VG(J1,CUR_OBS(I)%S_INT(K,2))
           W1 = CUR_OBS(I)%S_WEIGHT(K,1)
           W2 = CUR_OBS(I)%S_WEIGHT(K,2)
           UINT(I,K) = UINT(I,K) + (U1*W1 + U2*W2)*WGHT 
           VINT(I,K) = VINT(I,K) + (V1*W1 + V2*W2)*WGHT 
         END DO
       END DO
       TOT_WGHT = SUM(CUR_OBS(I)%X_WEIGHT(1:CUR_OBS(I)%N_INTPTS))
       UINT(I,1:NLAY) = UINT(I,1:NLAY)/TOT_WGHT
       VINT(I,1:NLAY) = VINT(I,1:NLAY)/TOT_WGHT
     END DO
   ELSE IF(CURRENT_METHOD == 'OI')THEN
     DO I=1,N_ASSIM_CUR  
       IP = CUR_OBS(I)%N_CELL
       NLAY      = CUR_OBS(I)%N_LAYERS
       XSC = CUR_OBS(I)%X-XCG(IP)
       YSC = CUR_OBS(I)%Y-YCG(IP)
       JN1 = NVG(IP,1)
       JN2 = NVG(IP,2)
       JN3 = NVG(IP,3)
      
       DO K=1,NLAY
         K1 = CUR_OBS(I)%S_INT(K,1)
         K2 = CUR_OBS(I)%S_INT(K,2)
         COFT0 = AW0G(IP,1)*UG(JN1,K1)+AW0G(IP,2)*UG(JN2,K1)+AW0G(IP,3)*UG(JN3,K1)
         COFTX = AWXG(IP,1)*UG(JN1,K1)+AWXG(IP,2)*UG(JN2,K1)+AWXG(IP,3)*UG(JN3,K1)
         COFTY = AWYG(IP,1)*UG(JN1,K1)+AWYG(IP,2)*UG(JN2,K1)+AWYG(IP,3)*UG(JN3,K1)
         U1 = COFT0 + COFTX*XSC + COFTY*YSC
         COFT0 = AW0G(IP,1)*UG(JN1,K2)+AW0G(IP,2)*UG(JN2,K2)+AW0G(IP,3)*UG(JN3,K2)
         COFTX = AWXG(IP,1)*UG(JN1,K2)+AWXG(IP,2)*UG(JN2,K2)+AWXG(IP,3)*UG(JN3,K2)
         COFTY = AWYG(IP,1)*UG(JN1,K2)+AWYG(IP,2)*UG(JN2,K2)+AWYG(IP,3)*UG(JN3,K2)
         U2 = COFT0 + COFTX*XSC + COFTY*YSC
         COFT0 = AW0G(IP,1)*VG(JN1,K1)+AW0G(IP,2)*VG(JN2,K1)+AW0G(IP,3)*VG(JN3,K1)
         COFTX = AWXG(IP,1)*VG(JN1,K1)+AWXG(IP,2)*VG(JN2,K1)+AWXG(IP,3)*VG(JN3,K1)
         COFTY = AWYG(IP,1)*VG(JN1,K1)+AWYG(IP,2)*VG(JN2,K1)+AWYG(IP,3)*VG(JN3,K1)
         V1 = COFT0 + COFTX*XSC + COFTY*YSC
         COFT0 = AW0G(IP,1)*VG(JN1,K2)+AW0G(IP,2)*VG(JN2,K2)+AW0G(IP,3)*VG(JN3,K2)
         COFTX = AWXG(IP,1)*VG(JN1,K2)+AWXG(IP,2)*VG(JN2,K2)+AWXG(IP,3)*VG(JN3,K2)
         COFTY = AWYG(IP,1)*VG(JN1,K2)+AWYG(IP,2)*VG(JN2,K2)+AWYG(IP,3)*VG(JN3,K2)
         V2 = COFT0 + COFTX*XSC + COFTY*YSC
         W1 = CUR_OBS(I)%S_WEIGHT(K,1)
         W2 = CUR_OBS(I)%S_WEIGHT(K,2)
         UINT(I,K) = U1*W1 + U2*W2 
         VINT(I,K) = V1*W1 + V2*W2 
       END DO
     END DO
   ELSE
     WRITE(IPT,*)'CURRENT_METHOD SHOULD BE NG OR OI. BUT HERE CURRENT_METHOD=',  &
                  TRIM(CURRENT_METHOD)
     CALL PSTOP
   END IF
     
!------------------------------------------------------------------------------!
!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!------------------------------------------------------------------------------!

   ALLOCATE(TWGHT(MElemGL,KBM1))   ; TWGHT = 0.
   ALLOCATE(UCORR(MElemGL,KBM1))   ; UCORR   = 0.
   ALLOCATE(VCORR(MElemGL,KBM1))   ; VCORR   = 0.

   IF(CURRENT_METHOD == 'NG')THEN
     DO I=1,N_ASSIM_CUR
       DO J=1,CUR_OBS(I)%N_INTPTS
         J1     = CUR_OBS(I)%INTPTS(J)
         ITIME  = CUR_OBS(I)%N_T_WEIGHT
         NLAY   = CUR_OBS(I)%N_LAYERS
         DO K=1,NLAY
           K1           = CUR_OBS(I)%S_INT(K,1)
           K2           = CUR_OBS(I)%S_INT(K,2)
           W1           = CUR_OBS(I)%S_WEIGHT(K,1)
           W2           = CUR_OBS(I)%S_WEIGHT(K,2)
           WEIGHT1      = CUR_OBS(I)%T_WEIGHT*CUR_OBS(I)%X_WEIGHT(J)*W1
           WEIGHT2      = CUR_OBS(I)%T_WEIGHT*CUR_OBS(I)%X_WEIGHT(J)*W2
           TWGHT(J1,K1) = TWGHT(J1,K1) + WEIGHT1   
           TWGHT(J1,K2) = TWGHT(J1,K2) + WEIGHT2   
           DEFECT       = CUR_OBS(I)%UO(ITIME,K) - UINT(I,K)
           CORRECTION   = GAMA_CUR*DEFECT
!         UCORR(J1,K1) = UCORR(J1,K1) + CORRECTION*WEIGHT1
!         UCORR(J1,K2) = UCORR(J1,K2) + CORRECTION*WEIGHT2
           UCORR(J1,K1) = UCORR(J1,K1) + CORRECTION*WEIGHT1**2
           UCORR(J1,K2) = UCORR(J1,K2) + CORRECTION*WEIGHT2**2
           DEFECT       = CUR_OBS(I)%VO(ITIME,K) - VINT(I,K)
           CORRECTION   = GAMA_CUR*DEFECT
!         VCORR(J1,K1) = VCORR(J1,K1) + CORRECTION*WEIGHT1
!         VCORR(J1,K2) = VCORR(J1,K2) + CORRECTION*WEIGHT2
           VCORR(J1,K1) = VCORR(J1,K1) + CORRECTION*WEIGHT1**2
           VCORR(J1,K2) = VCORR(J1,K2) + CORRECTION*WEIGHT2**2
!        GEOFF NEW
           IF(ABS(CUR_OBS(I)%UO(ITIME,K)) + ABS(CUR_OBS(I)%VO(ITIME,K)) < .0001)THEN
             TWGHT(J1,K1) = 0.
             TWGHT(J1,K2) = 0.
           END IF
         END DO
       END DO
     END DO

!------------------------------------------------------------------------------!
!  Nudge Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     DO I=1,MElemGL
       DO K=1,KBM1
         IF(DA_CUR(I) == 1 .AND. TWGHT(I,K) > 1.0E-08)THEN
           UG(I,K) = UG(I,K) + DTI*GALPHA_CUR*UCORR(I,K)/TWGHT(I,K)
           VG(I,K) = VG(I,K) + DTI*GALPHA_CUR*VCORR(I,K)/TWGHT(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT,UCORR,VCORR,UINT,VINT)
   ELSE IF(CURRENT_METHOD == 'OI')THEN
     DO I=1,N_ASSIM_CUR
       ITIME  = CUR_OBS(I)%N_T_WEIGHT
       NLAY   = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         K1           = CUR_OBS(I)%S_INT(K,1)
         K2           = CUR_OBS(I)%S_INT(K,2)
         W1           = CUR_OBS(I)%S_WEIGHT(K,1)
         W2           = CUR_OBS(I)%S_WEIGHT(K,2)
         WEIGHT1      = CUR_OBS(I)%T_WEIGHT*W1
         WEIGHT2      = CUR_OBS(I)%T_WEIGHT*W2
         TWGHT(I,K1) = TWGHT(I,K1) + WEIGHT1   
         TWGHT(I,K2) = TWGHT(I,K2) + WEIGHT2   
         DEFECT       = CUR_OBS(I)%UO(ITIME,K) - UINT(I,K)
         CORRECTION   = DEFECT
         UCORR(I,K1) = UCORR(I,K1) + CORRECTION*WEIGHT1
         UCORR(I,K2) = UCORR(I,K2) + CORRECTION*WEIGHT2
         DEFECT       = CUR_OBS(I)%VO(ITIME,K) - VINT(I,K)
         CORRECTION   = DEFECT
         VCORR(I,K1) = VCORR(I,K1) + CORRECTION*WEIGHT1
         VCORR(I,K2) = VCORR(I,K2) + CORRECTION*WEIGHT2
         IF(ABS(CUR_OBS(I)%UO(ITIME,K)) + ABS(CUR_OBS(I)%VO(ITIME,K)) < .0001)THEN
           TWGHT(I,K1) = 0.
           TWGHT(I,K2) = 0.
         END IF
       END DO
     END DO

     DO I=1,N_ASSIM_CUR
       DO K=1,KBM1
         IF(TWGHT(I,K) > 1.0E-8)THEN
           UCORR(I,K)=UCORR(I,K)/TWGHT(I,K)
	   VCORR(I,K)=VCORR(I,K)/TWGHT(I,K)
         END IF
       END DO
     END DO
       	 
!------------------------------------------------------------------------------!
!  Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     ALLOCATE(UCORR1(MElemGL,KBM1));UCORR1=0.0_SP
     ALLOCATE(VCORR1(MElemGL,KBM1));VCORR1=0.0_SP
   
     DO K=1,KBM1
       CALL CUR_OPTIMINTERP(UCORR(:,K),VCORR(:,K),UCORR1(:,K),VCORR1(:,K))
     END DO

     DO I=1,MElemGL
       DO K=1,KBM1
         IF(DA_CUR(I) == 1)THEN
!          UG(I,K) = UG(I,K)+UCORR1(I,K)
!          VG(I,K) = VG(I,K)+VCORR1(I,K)
           UG(I,K) = UG(I,K)+ GALPHA_CUR*UCORR1(I,K)
           VG(I,K) = VG(I,K)+ GALPHA_CUR*VCORR1(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT,UCORR,VCORR,UINT,VINT)
     DEALLOCATE(UCORR1,VCORR1)
   ELSE
     WRITE(IPT,*)'CURRENT_METHOD SHOULD BE NG OR OI. BUT HERE CURRENT_METHOD=',  &
                  TRIM(CURRENT_METHOD)
     CALL PSTOP
   END IF

   END IF  !!MASTER
       
!------------------------------------------------------------------------------!
!  Disperse New Data Fields to Slave Processors                                ! 
!------------------------------------------------------------------------------!
   IF(SERIAL)THEN
     UF(1:MElem,1:KBM1) = UG(1:MElem,1:KBM1)
     VF(1:MElem,1:KBM1) = VG(1:MElem,1:KBM1)
   END IF
#  if defined (MULTIPROCESSOR) 
   CALL MPI_BCAST(UG,MElemGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(VG,MElemGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)THEN
     DO I=1,MElem
       UF(I,1:KBM1) = UG(EGID(I),1:KBM1)
       VF(I,1:KBM1) = VG(EGID(I),1:KBM1)
     END DO
   END IF
#  endif

   DEALLOCATE(UG,VG)
          

   RETURN
   END SUBROUTINE CURRENT_ASSIMILATION
!==============================================================================|
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
!==============================================================================|
   SUBROUTINE CUR_OPTIMINTERP(F1,F2,FI1,FI2)
   USE MOD_OPTIMAL_INTERPOLATION
   USE MOD_PREC
   USE ALL_VARS
   IMPLICIT NONE

!------------------------------------------------------------------------------|
!  xi(1,:) and xi(2,:) represent the x and y coordindate of the grid of the    |
!  interpolated field                                                          |
!  fi and vari are the interpolated field and its error variance resp.         |
!------------------------------------------------------------------------------|
   REAL(SP) :: XI(2,MElemGL),FI1(MElemGL),FI2(MElemGL),VARI(MElemGL)

!------------------------------------------------------------------------------|
!  x(1,:) and x(2,:) represent the x and y coordindate of the observations     |
!  f and var are observations and their error variance resp.                   |
!------------------------------------------------------------------------------|
   REAL(SP) :: X(2,N_ASSIM_CUR),VAR(N_ASSIM_CUR),F1(N_ASSIM_CUR),F2(N_ASSIM_CUR)

!------------------------------------------------------------------------------|
!  param: inverse of the correlation length                                    |
!------------------------------------------------------------------------------|
   REAL(SP) :: PARAM(2,N_ASSIM_CUR)

   INTEGER  :: I,J,MM

!------------------------------------------------------------------------------|
!  create a regular 2D grid                                                    |
!------------------------------------------------------------------------------|
   DO I=1,MElemGL
     XI(1,I) = XCG(I)
     XI(2,I) = YCG(I)
   END DO	
   
!------------------------------------------------------------------------------|   
!  param is the inverse of the correlation length                              |
!------------------------------------------------------------------------------|
   PARAM = 1.0_SP/PARAM_CUR 
   MM = N_INFLU_CUR
 
!------------------------------------------------------------------------------| 
!  the error variance of the observations                                      |
!------------------------------------------------------------------------------|
   VAR = 0.0_SP   

   DO I=1,N_ASSIM_CUR
     X(1,I) = CUR_OBS(I)%X
     X(2,I) = CUR_OBS(I)%Y
   END DO  

!------------------------------------------------------------------------------|
!  fi is the interpolated function and vari its error variance                 |
!------------------------------------------------------------------------------|
   CALL OPTIMINTERP(X,F1,VAR,PARAM,MM,XI,FI1,VARI)
   CALL OPTIMINTERP(X,F2,VAR,PARAM,MM,XI,FI2,VARI)

   RETURN
   END SUBROUTINE CUR_OPTIMINTERP
!==============================================================================|
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|


   SUBROUTINE SET_SST_ASSIM_DATA 

!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR SST OBSERVATIONS                               |
!------------------------------------------------------------------------------!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif      
   IMPLICIT NONE
   INTEGER I,J,K,ECNT,ITMP,NCNT,IOS
   CHARACTER(LEN=120) :: FNAME
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,TEMPF
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP
   LOGICAL :: FEXIST
       
!------------------------------------------------------------------------------!
       
   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_sst.dat"
!
!--Make Sure SST Assimilation Data File Exists---------------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'SST OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number SST Measurements and Data Set Size------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD') ; REWIND(1)
   READ(1,*) N_TIMES_SST 
   READ(1,*) TEMPF,N_ASSIM_SST 
   ALLOCATE(SST_OBS(N_ASSIM_SST))
   DO I=1,N_ASSIM_SST
     SST_OBS(I)%N_TIMES = N_TIMES_SST
     ALLOCATE(SST_OBS(I)%SST(N_TIMES_SST))
     ALLOCATE(SST_OBS(I)%TIMES(N_TIMES_SST))
     ALLOCATE(SST_OBS(I)%T_INT_HOUR(N_DA_HOURS)) ; SST_OBS(I)%T_INT_HOUR = 0.
   END DO
   REWIND(1)  ; READ(1,*)
   
!
!--Read SST Data---------------------------------------------------------------!
!
   DO I=1,N_TIMES_SST
     READ(1,*) TEMPF,N_ASSIM_SST 
     DO J=1,N_ASSIM_SST
       READ(1,*)SST_OBS(J)%X,SST_OBS(J)%Y,SST_OBS(J)%SST(I)
       SST_OBS(J)%TIMES(I) = TEMPF
     END DO
   END DO
        write(*,*)'x and y from input = ',sst_obs%x,sst_obs%y  ! B Clark Debug
!
!--Shift Coordinates-of SST Observation Locations------------------------------!
!       
   SST_OBS(:)%X = SST_OBS(:)%X - VXMIN 
   SST_OBS(:)%Y = SST_OBS(:)%Y - VYMIN 
        write(*,*) 'vxmin, vymin = ',vxmin,vymin ! B Clark debug  
        
   IF(SST_METHOD == 'OI')THEN
!
!--find the cell number (SST_OBS(:)%N_CELL) of Obs station---------------------
!
   RRTMP = 100000.0     !100km
   DO J= 1,N_ASSIM_SST 
      X0 = SST_OBS(J)%X
      Y0 = SST_OBS(J)%Y
      DO I=1,MElemGL
         Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
        if(Rtmp.LT.RRTMP) then
	    S11 = (XG(NVG(I,2))-X0)*(YG(NVG(I,3))-Y0)-&
	          (XG(NVG(I,3))-X0)*(YG(NVG(I,2))-Y0)
            S22 = (XG(NVG(I,3))-X0)*(YG(NVG(I,1))-Y0)-&
	          (XG(NVG(I,1))-X0)*(YG(NVG(I,3))-Y0)
            S33 = (XG(NVG(I,1))-X0)*(YG(NVG(I,2))-Y0)-&
	          (XG(NVG(I,2))-X0)*(YG(NVG(I,1))-Y0)
	    IF(S11.LE.0.AND.S22.LE.0.AND.S33.LE.0) THEN
	      SST_OBS(J)%N_CELL = I
	      GOTO 300
	    ELSE
	      SST_OBS(J)%N_CELL = 0
	    ENDIF
	 ELSE
	    SST_OBS(J)%N_CELL = -1
	 ENDIF
      ENDDO	
300 CONTINUE
      IF(SST_OBS(J)%N_CELL.LE.0) THEN
         IF(MSR) WRITE(IPT,*)'ERROR--SST OBS SITE:',J,' OUT OF DOMAN',&
	         SST_OBS(J)%N_CELL 
         CALL PSTOP
      ENDIF 	 

   ENDDO             
!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   ALLOCATE(AW0G(MElemGL,3))
   ALLOCATE(AWXG(MElemGL,3))
   ALLOCATE(AWYG(MElemGL,3))
   
   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN     
     CALL GATHER(LBOUND(AW0,1),UBOUND(AW0,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AW0,AW0G)
     CALL GATHER(LBOUND(AWX,1),UBOUND(AWX,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AWX,AWXG)
     CALL GATHER(LBOUND(AWY,1),UBOUND(AWY,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AWY,AWYG)
   END IF
#  endif     
   END IF
!
!--Close SST Observation Data File---------------------------------------------!
!
   CLOSE(1)
       
   IF(SST_METHOD == 'OI')THEN
!------------------------------------------------------------------------------!
!  Read Correlation Length of SST Observations                             !
!------------------------------------------------------------------------------!
       
!JQI   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_radius_sst.dat"

!JQI   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
!JQI   IF(MSR .AND. .NOT.FEXIST)THEN
!JQI     WRITE(IPT,*)'SST OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
!JQI     WRITE(IPT,*)'HALTING.....'
!JQI     CALL PSTOP
!JQI   END IF
     
!JQI   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')

   ALLOCATE(PARAM_SST(2,N_ASSIM_SST))

!JQI   DO I=1,N_ASSIM_SST
!JQI     READ(1,*)PARAM_SST(1,I),PARAM_SST(2,I)
!JQI   END DO

!JQI   CLOSE(1)

   PARAM_SST = 10000.0_SP
   END IF
       
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Observation Location 
!------------------------------------------------------------------------------!

   ALLOCATE(ITEMP(NNodeGL),FTEMP(NNodeGL),DA_SST(NNodeGL))     ; DA_SST = 0
   DO I=1,N_ASSIM_SST
     X0 = SST_OBS(I)%X
     Y0 = SST_OBS(I)%Y
     ECNT = 0
     DO J=1,NNodeGL
       DX = ABS(XG(J)-X0)
       DY = ABS(YG(J)-Y0)
       RD = SQRT(DX**2 + DY**2)
       IF(RD <= RAD_SST)THEN
         DA_SST(J) = 1
         ECNT = ECNT + 1      
         ITEMP(ECNT) =  J
         FTEMP(ECNT) = (RAD_SST**2 - RD**2) / (RAD_SST**2 + RD**2)
       END IF
     END DO
     IF(ECNT == 0)THEN
       SST_OBS(I)%N_INTPTS = ECNT
       ALLOCATE(SST_OBS(I)%INTPTS(1))
       ALLOCATE(SST_OBS(I)%X_WEIGHT(1))
     ELSE
       SST_OBS(I)%N_INTPTS = ECNT
       ALLOCATE(SST_OBS(I)%INTPTS(ECNT))
       ALLOCATE(SST_OBS(I)%X_WEIGHT(ECNT))
       SST_OBS(I)%INTPTS(1:ECNT)  = ITEMP(1:ECNT)
       SST_OBS(I)%X_WEIGHT(1:ECNT) = FTEMP(1:ECNT)
     END IF
   END DO
   DEALLOCATE(FTEMP,ITEMP)
     

!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN
   WRITE(IPT,*)
   WRITE(IPT,*)'!                SST OBSERVATION DATA           '
!   WRITE(IPT,*)"  OBS#      X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES"
!   DO I=1,N_ASSIM_SST
!     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X)') &
!       I,SST_OBS(I)%X/1000.,SST_OBS(I)%Y/1000.,SST_OBS(I)%N_INTPTS,SST_OBS(I)%N_TIMES
!   END DO
   WRITE(IPT,*)'NUMBER OF NODES WITHOUT NUDGING: ',NNodeGL-SUM(DA_SST(1:NNodeGL))
   WRITE(IPT,*)'MAXIMUM NUMBER OF INTERP POINTS: ',MAXVAL(SST_OBS(1:N_ASSIM_SST)%N_INTPTS)
   WRITE(IPT,*)
   RETURN
   END SUBROUTINE SET_SST_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE SST_ASSIMILATION
!==============================================================================|
!  USE SST OBSERVATION DATA TO ADJUST SURFACE TEMPERATURE WITH NUDGE OR OI     |
!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: T1G
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TINT,TCORR,TCORR1
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TWEIGHT,FTEMP
   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
   REAL(SP) :: TRUTH,ttemp,TIME_24,TMP,ASTIME,tmax,TSHIFT
   REAL(SP) :: U1,U2,V1,V2   !,W1,W2,WEIGHT1,WEIGHT2
   INTEGER I,J,J1,ITIME,NTIME,IERR,ITMP,jtmax
   INTRINSIC MINLOC

   INTEGER IP,JN1,JN2,JN3
   REAL(SP) XSC,YSC,COFT0,COFTX,COFTY
!==============================================================================|

   IF(MSR)THEN
     IF(SST_CYC > N_TIMES_SST)THEN
     WRITE(IPT,*)'NUMBER OF DA SWEEPS EXCEEDS NUMBER OF SST DATA INTERVALS'
     WRITE(IPT,*)'NUMBER OF SWEEPS: ',SST_CYC
     WRITE(IPT,*)'NUMBER OF DATA INTERVALS: ',N_TIMES_SST
     WRITE(IPT,*)'HALTING'
     CALL PSTOP
     END IF
   END IF
!------------------------------------------------------------------------------!
!  Calculate Hour Number of Current Assimilation Cycle                         ! 
!------------------------------------------------------------------------------!

   TSHIFT = THOUR-FLOAT(ISTART_DAY*24)-FLOAT((SST_CYC-1)*IAV_DAY*24)
   DA_HOUR = INT(TSHIFT+.5)
   SST_OBS(:)%N_T_WEIGHT = DA_HOUR
   IF(DA_HOUR == 0)DA_HOUR = 1
   ASTIME = ABS(TSHIFT-DA_HOUR)*3600.


!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!


   IF(ASTIME < 0.5*ASTIME_WINDOW_SST) THEN
     SST_OBS(:)%T_WEIGHT = 1.0
   ELSE IF(ASTIME < ASTIME_WINDOW_SST) THEN
     SST_OBS(:)%T_WEIGHT = (ASTIME_WINDOW_SST-ASTIME)/ASTIME_WINDOW_SST*2.0
   ELSE
     SST_OBS(:)%T_WEIGHT = 0.0
   END IF

       
!------------------------------------------------------------------------------!
!  Gather Surface Temperature Field to Master Processor                        ! 
!------------------------------------------------------------------------------!
   ALLOCATE(T1G(NNodeGL,1:KB))
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     CALL GATHER(LBOUND(T1,1),UBOUND(T1,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,T1,T1G)
   END IF
#  endif
   IF(SERIAL)THEN
     T1G(1:NNodeGL,1) = T1(1:NNodeGL,1)
   END IF
       
!------------------------------------------------------------------------------!
!  Interpolate Simulation Data to Local Observation Point                      ! 
!------------------------------------------------------------------------------!
       
   IF(MSR)THEN
   ALLOCATE(TINT(N_ASSIM_SST)) ; TINT = 0. 

   IF(SST_METHOD == 'NG')THEN
     DO I=1,N_ASSIM_SST  
       IF(SST_OBS(I)%N_INTPTS > 0)THEN
         DO J=1,SST_OBS(I)%N_INTPTS
           J1        = SST_OBS(I)%INTPTS(J)
           WGHT      = SST_OBS(I)%X_WEIGHT(J)
           TINT(I)   = TINT(I) + SST_OBS(I)%X_WEIGHT(J)*T1G(J1,1) 
         END DO
         TOT_WGHT = SUM(SST_OBS(I)%X_WEIGHT(1:SST_OBS(I)%N_INTPTS))
         TINT(I)  = TINT(I)/TOT_WGHT
       END IF 
     END DO
   ELSE IF(SST_METHOD == 'OI')THEN
     DO I=1,N_ASSIM_SST  
       IP   = SST_OBS(I)%N_CELL
       XSC  = SST_OBS(I)%X-XCG(IP)
       YSC  = SST_OBS(I)%Y-YCG(IP)
       JN1  = NVG(IP,1)
       JN2  = NVG(IP,2)
       JN3  = NVG(IP,3)
      
       COFT0 = AW0G(IP,1)*T1G(JN1,1)+AW0G(IP,2)*T1G(JN2,1)+AW0G(IP,3)*T1G(JN3,1)
       COFTX = AWXG(IP,1)*T1G(JN1,1)+AWXG(IP,2)*T1G(JN2,1)+AWXG(IP,3)*T1G(JN3,1)
       COFTY = AWYG(IP,1)*T1G(JN1,1)+AWYG(IP,2)*T1G(JN2,1)+AWYG(IP,3)*T1G(JN3,1)
       U1 = COFT0 + COFTX*XSC + COFTY*YSC
       TINT(I) = U1
     END DO
   ELSE
     WRITE(IPT,*)'SST_METHOD SHOULD BE NG OR OI. BUT HERE SST_METHOD=',TRIM(SST_METHOD)
     CALL PSTOP
   END IF
     
!------------------------------------------------------------------------------!
!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!------------------------------------------------------------------------------!

   ALLOCATE(TWEIGHT(NNodeGL))    ; TWEIGHT = 0.
   ALLOCATE(TCORR(NNodeGL))      ; TCORR   = 0.

   IF(SST_METHOD == 'NG')THEN
     DO I=1,N_ASSIM_SST
       DO J=1,SST_OBS(I)%N_INTPTS
         J1           = SST_OBS(I)%INTPTS(J)
         WEIGHT       = SST_OBS(I)%T_WEIGHT*SST_OBS(I)%X_WEIGHT(J)
         TWEIGHT(J1)  = TWEIGHT(J1) + WEIGHT
         TRUTH        = &
           SST_OBS(I)%T_INT_HOUR(DA_HOUR) + SST_OBS(I)%SST(SST_CYC) - SST_OBS(I)%T_INT_AVGD
         DEFECT       = TRUTH - TINT(I)
         CORRECTION   = GAMA_SST*DEFECT*WEIGHT**2
         TCORR(J1)    = TCORR(J1) + CORRECTION
       END DO
     END DO

!------------------------------------------------------------------------------!
!  Nudge Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     DO I=1,NNodeGL
       IF(DA_SST(I) == 1 .AND. TWEIGHT(I) > 1.0E-08)THEN
         T1G(I,1) = T1G(I,1) + DTI*GALPHA_SST*TCORR(I)/TWEIGHT(I)
       END IF
     END DO

     DEALLOCATE(TWEIGHT,TCORR,TINT)
   ELSE IF(SST_METHOD == 'OI')THEN
     DO I=1,N_ASSIM_SST
       WEIGHT     = SST_OBS(I)%T_WEIGHT
       TWEIGHT(I) = TWEIGHT(I) + WEIGHT
       TRUTH      = &
         SST_OBS(I)%T_INT_HOUR(DA_HOUR) + SST_OBS(I)%SST(SST_CYC) - SST_OBS(I)%T_INT_AVGD
       DEFECT       = TRUTH - TINT(I)
       CORRECTION   = DEFECT*WEIGHT
       TCORR(I)    = TCORR(I) + CORRECTION
     END DO

     DO I=1,N_ASSIM_SST
       IF(TWEIGHT(I) > 1.0E-8)THEN
         TCORR(I)=TCORR(I)/TWEIGHT(I)
       END IF
     END DO

!------------------------------------------------------------------------------!
!  'OI' Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     ALLOCATE(TCORR1(NNodeGL)); TCORR1=0.0_SP
     CALL SST_OPTIMINTERP(TCORR(:),TCORR1(:))

     DO I=1,NNodeGL
       IF(DA_SST(I) == 1)THEN
!        T1G(I,1) = T1G(I,1) + TCORR1(I)
         T1G(I,1) = T1G(I,1) + GALPHA_SST*TCORR1(I)
       END IF
     END DO

     DEALLOCATE(TWEIGHT,TCORR,TINT)
   ELSE
     WRITE(IPT,*)'SST_METHOD SHOULD BE NG OR OI. BUT HERE SST_METHOD=',TRIM(SST_METHOD)
     CALL PSTOP
   END IF
   
   END IF !!MSR
       
!------------------------------------------------------------------------------!
!  Disperse New Data Fields to Slave Processors                                ! 
!------------------------------------------------------------------------------!
   IF(SERIAL)THEN
     T1(1:NNode,1) = T1G(1:NNode,1)
   END IF
#  if defined (MULTIPROCESSOR) 
   CALL MPI_BCAST(T1G,NNodeGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)THEN
     DO I=1,NNode
       T1(I,1) = T1G(NGID(I),1)
     END DO
   END IF
#  endif

   DEALLOCATE(T1G)

   RETURN
   END SUBROUTINE SST_ASSIMILATION
!==============================================================================|
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
!==============================================================================|
   SUBROUTINE SST_OPTIMINTERP(F,FI)
   USE MOD_OPTIMAL_INTERPOLATION
   USE MOD_PREC
   USE ALL_VARS
   IMPLICIT NONE

!------------------------------------------------------------------------------!
! xi(1,:) and xi(2,:) represent the x and y coordindate of the grid of the     !
! interpolated field                                                           !
! fi and vari are the interpolated field and its error variance resp.          !
!------------------------------------------------------------------------------!
   REAL(SP) :: XI(2,NNodeGL),FI(NNodeGL),VARI(NNodeGL)

!------------------------------------------------------------------------------!
! x(1,:) and x(2,:) represent the x and y coordindate of the observations      !
! f and var are observations and their error variance resp.                    !
!------------------------------------------------------------------------------!
   REAL(SP) :: X(2,N_ASSIM_SST),VAR(N_ASSIM_SST),F(N_ASSIM_SST)

!------------------------------------------------------------------------------!
! param: inverse of the correlation length                                     !
!------------------------------------------------------------------------------!
   REAL(SP) :: PARAM(2,N_ASSIM_SST)

   INTEGER  :: I,J,MM 

!------------------------------------------------------------------------------!
! create a regular 2D grid                                                     !
!------------------------------------------------------------------------------!
   DO I=1,NNodeGL
     XI(1,I) = XG(I)
     XI(2,I) = YG(I)
   END DO	
   
!------------------------------------------------------------------------------!
! param is the inverse of the correlation length                               !
!------------------------------------------------------------------------------!
   PARAM = 1.0_SP/PARAM_SST            
   
   MM = N_INFLU_SST

!------------------------------------------------------------------------------|   
! the error variance of the observations                                       |
!------------------------------------------------------------------------------|
   VAR = 0.0_SP   

!------------------------------------------------------------------------------|
! location of observations                                                     |
!------------------------------------------------------------------------------|
   DO I=1,N_ASSIM_SST
     X(1,I) = SST_OBS(I)%X
     X(2,I) = SST_OBS(I)%Y
   END DO  

!------------------------------------------------------------------------------|
! fi is the interpolated function and vari its error variance                  |
!------------------------------------------------------------------------------|
   CALL OPTIMINTERP(X,F,VAR,PARAM,MM,XI,FI,VARI)

   RETURN
   END SUBROUTINE SST_OPTIMINTERP
!==============================================================================|
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE SST_INT
!==============================================================================|
!  INTERPOLATE HOURLY AND INTERVAL-AVERAGAED SST DATA TO OBSERVATIN POINTS     | 
!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: T_TMP_L,T_TMP_G 
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TWEIGHT,FTEMP
   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
   INTEGER I,J,J1,ITIME,NTIME,IERR
   INTRINSIC MINLOC
!==============================================================================|

!------------------------------------------------------------------------------!
!  Calculate Hour Number of Current Assimilation Cycle                         ! 
!------------------------------------------------------------------------------!

   IF(ABS(THOUR-INT(THOUR+0.5)) < 1.E-8) THEN
     DA_HOUR = MOD(INT(THOUR-FLOAT(ISTART_DAY*24)+0.5),24*IAV_DAY)
     IF(DA_HOUR == 0) DA_HOUR = 24*IAV_DAY
   ELSE
     RETURN
   END IF

!------------------------------------------------------------------------------!
!  Extract Sea Surface Temp and Gather to Global Array                         ! 
!------------------------------------------------------------------------------!

   ALLOCATE(T_TMP_G(NNodeGL))
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(T_TMP_L(0:NNode))
     T_TMP_L(1:NNode) = T1(1:NNode,1)
     CALL GATHER(0,NNode,NNode,NNodeGL,1,MYID,NPROCS,NMAP,T_TMP_L,T_TMP_G)
   END IF
#  endif
   IF(SERIAL)THEN
     T_TMP_G(1:NNodeGL) = T1(1:NNodeGL,1)
   END IF



!------------------------------------------------------------------------------!
!  Interpolate Hourly Data to Local Observation Point                          ! 
!------------------------------------------------------------------------------!
   IF(MSR)THEN
   DO I=1,N_ASSIM_SST
     SST_OBS(I)%T_INT_HOUR(DA_HOUR) = 0.
     IF(SST_OBS(I)%N_INTPTS /= 0)THEN
     DO J=1,SST_OBS(I)%N_INTPTS
       J1        = SST_OBS(I)%INTPTS(J)
       WGHT      = SST_OBS(I)%X_WEIGHT(J)
       SST_OBS(I)%T_INT_HOUR(DA_HOUR)   = &
               SST_OBS(I)%T_INT_HOUR(DA_HOUR) + WGHT*T_TMP_G(J1)
     END DO
     TOT_WGHT = SUM(SST_OBS(I)%X_WEIGHT(1:SST_OBS(I)%N_INTPTS))
     SST_OBS(I)%T_INT_HOUR(DA_HOUR)  = SST_OBS(I)%T_INT_HOUR(DA_HOUR)/TOT_WGHT
     END IF
   END DO
  
!   WRITE(IPT,*)'  INTERPOLATING HOURLY SST DATA',DA_HOUR,N_DA_HOURS
   END IF

   DEALLOCATE(T_TMP_G)
   IF(PAR) DEALLOCATE(T_TMP_L)

!------------------------------------------------------------------------------!
!  Average Hourly Data Over Integration Period                                 ! 
!------------------------------------------------------------------------------!

   IF(DA_HOUR /=  N_DA_HOURS) THEN
      write(*,*)'exiting subroutine SST_INT' !B Clark Debug
   ENDIF     

   IF(DA_HOUR /=  N_DA_HOURS) RETURN
   
   IF(MSR)THEN
   DO I=1,N_ASSIM_SST
     SST_OBS(I)%T_INT_AVGD = SUM(SST_OBS(I)%T_INT_HOUR(1:N_DA_HOURS))/FLOAT(N_DA_HOURS)
   END DO
   WRITE(IPT,*)'  COMPUTING AVERAGE INTERPOLATED SST DATA',DA_HOUR,N_DA_HOURS
   END IF
   
   RETURN
   END SUBROUTINE SST_INT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
!   READ IN DATA FROM SIMULATION STAGE AND FOR ASSIMILATION STAGE STARTUP      |
!==============================================================================|

   SUBROUTINE HOT_START_SST      

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   INTEGER :: I,K
   LOGICAL :: FEXIST


   IF(MSR)THEN
     INQUIRE(FILE='restart_sst',EXIST=FEXIST)
     IF(.NOT. FEXIST)THEN
       WRITE(IPT,*)'FILE restart_sst DOES NOT EXIST'
       WRITE(IPT,*)'CANNOT START ASSIMILATION STAGE'
       WRITE(IPT,*)'HALTING....'
       CALL PSTOP
     END IF
   END IF

   OPEN(1,FILE='restart_sst',FORM='UNFORMATTED')
   REWIND(1)

   IF(SERIAL)THEN
     REWIND(1)
     READ(1) IINT
     READ(1) ((U(I,K),K=1,KB),I=0,MElem)
     READ(1) ((V(I,K),K=1,KB),I=0,MElem)
     READ(1) ((W(I,K),K=1,KB),I=0,MElem)
     
     READ(1) ((Q2(I,K),K=1,KB),I=0,NNode)
     READ(1) ((Q2L(I,K),K=1,KB),I=0,NNode)
     
     READ(1) ((S(I,K),K=1,KB),I=0,MElem)
     READ(1) ((T(I,K),K=1,KB),I=0,MElem)
     READ(1) ((RHO(I,K),K=1,KB),I=0,MElem)
     READ(1) ((TMEAN(I,K),K=1,KB),I=0,MElem)
     READ(1) ((SMEAN(I,K),K=1,KB),I=0,MElem)
     READ(1) ((RMEAN(I,K),K=1,KB),I=0,MElem)

     READ(1) ((S1(I,K),K=1,KB),I=1,NNode)
     READ(1) ((T1(I,K),K=1,KB),I=1,NNode)
     READ(1) ((RHO1(I,K),K=1,KB),I=1,NNode)
     READ(1) ((TMEAN1(I,K),K=1,KB),I=1,NNode)
     READ(1) ((SMEAN1(I,K),K=1,KB),I=1,NNode)
     READ(1) ((RMEAN1(I,K),K=1,KB),I=1,NNode)
     
     READ(1) ((KM(I,K),K=1,KB),I=1,NNode)
     READ(1) ((KH(I,K),K=1,KB),I=1,NNode)
     READ(1) ((KQ(I,K),K=1,KB),I=1,NNode)
     READ(1)  ((L(I,K),K=1,KB),I=1,NNode)

     READ(1) (UA(I), I=0,MElem)
     READ(1) (VA(I), I=0,MElem)
     READ(1) (EL1(I), I=1,MElem)
     READ(1) (ET1(I), I=1,MElem)
     READ(1) (H1(I), I=1,MElem)
     READ(1) (D1(I), I=1,MElem)
     READ(1) (DT1(I), I=1,MElem)

     READ(1) (EL(I), I=1,NNode)
     READ(1) (ET(I), I=1,NNode)
     READ(1) (H(I), I=1,NNode)
     READ(1) (D(I), I=1,NNode)
     READ(1) (DT(I), I=1,NNode)
     CLOSE(1)
   ELSE
#  if defined (MULTIPROCESSOR)
     REWIND(1)
     READ(1) IINT
     CALL PREAD(1,U     ,LBOUND(U,1),    UBOUND(U,1),    MElem,MElemGL,KB,EGID(1),0,"U"     )
     CALL PREAD(1,V     ,LBOUND(V,1),    UBOUND(V,1),    MElem,MElemGL,KB,EGID(1),0,"V"     )
     CALL PREAD(1,W     ,LBOUND(W,1),    UBOUND(W,1),    MElem,MElemGL,KB,EGID(1),0,"W"     )
     
     CALL PREAD(1,Q2    ,LBOUND(Q2,1),   UBOUND(Q2,1),   NNode,NNodeGL,KB,NGID(1),0,"Q2"    )
     CALL PREAD(1,Q2L   ,LBOUND(Q2L,1),  UBOUND(Q2L,1),  NNode,NNodeGL,KB,NGID(1),0,"Q2L"   )
     
     CALL PREAD(1,S     ,LBOUND(S,1),    UBOUND(S,1),    MElem,MElemGL,KB,EGID(1),0,"S"     )
     CALL PREAD(1,T     ,LBOUND(S,1),    UBOUND(T,1),    MElem,MElemGL,KB,EGID(1),0,"T"     )
     CALL PREAD(1,RHO   ,LBOUND(RHO,1),  UBOUND(RHO,1),  MElem,MElemGL,KB,EGID(1),0,"RHO"   )
     CALL PREAD(1,TMEAN ,LBOUND(TMEAN,1),UBOUND(TMEAN,1),MElem,MElemGL,KB,EGID(1),0,"TMEAN" )
     CALL PREAD(1,SMEAN ,LBOUND(SMEAN,1),UBOUND(SMEAN,1),MElem,MElemGL,KB,EGID(1),0,"SMEAN" )
     CALL PREAD(1,RMEAN ,LBOUND(RMEAN,1),UBOUND(RMEAN,1),MElem,MElemGL,KB,EGID(1),0,"RMEAN" )

     CALL PREAD(1,S1    ,LBOUND(S1,1),    UBOUND(S1,1),    NNode,NNodeGL,KB,NGID,1,"S1"     )
     CALL PREAD(1,T1    ,LBOUND(T1,1),    UBOUND(T1,1),    NNode,NNodeGL,KB,NGID,1,"T1"     )
     CALL PREAD(1,RHO1  ,LBOUND(RHO1,1),  UBOUND(RHO1,1),  NNode,NNodeGL,KB,NGID,1,"RHO1"   )
     CALL PREAD(1,TMEAN1,LBOUND(TMEAN1,1),UBOUND(TMEAN1,1),NNode,NNodeGL,KB,NGID,1,"TMEAN1" )
     CALL PREAD(1,SMEAN1,LBOUND(SMEAN1,1),UBOUND(SMEAN1,1),NNode,NNodeGL,KB,NGID,1,"SMEAN1" )
     CALL PREAD(1,RMEAN1,LBOUND(RMEAN1,1),UBOUND(RMEAN1,1),NNode,NNodeGL,KB,NGID,1,"RMEAN1" )
  
     CALL PREAD(1,KM  ,LBOUND(KM,1),UBOUND(KM,1),NNode,NNodeGL,KB,NGID(1),1,"KM" )
     CALL PREAD(1,KH  ,LBOUND(KH,1),UBOUND(KH,1),NNode,NNodeGL,KB,NGID(1),1,"KH" )
     CALL PREAD(1,KQ  ,LBOUND(KQ,1),UBOUND(KQ,1),NNode,NNodeGL,KB,NGID(1),1,"KQ" )
     CALL PREAD(1,L   ,LBOUND( L,1),UBOUND( L,1),NNode,NNodeGL,KB,NGID(1),1,"L" )

     CALL PREAD(1,UA  ,LBOUND(UA,1), UBOUND(UA,1), MElem,MElemGL,1 ,EGID(1),0,"UA"  )
     CALL PREAD(1,VA  ,LBOUND(VA,1), UBOUND(VA,1), MElem,MElemGL,1 ,EGID(1),0,"VA"  )
     CALL PREAD(1,EL1 ,LBOUND(EL1,1),UBOUND(EL1,1),MElem,MElemGL,1 ,EGID(1),1,"EL1" )
     CALL PREAD(1,ET1 ,LBOUND(ET1,1),UBOUND(ET1,1),MElem,MElemGL,1 ,EGID(1),1,"ET1" )
     CALL PREAD(1,H1  ,LBOUND(H1,1), UBOUND(H1,1), MElem,MElemGL,1 ,EGID(1),1,"H1"  )
     CALL PREAD(1,D1  ,LBOUND(D1,1), UBOUND(D1,1), MElem,MElemGL,1 ,EGID(1),1,"D1"  )
     CALL PREAD(1,DT1 ,LBOUND(DT1,1),UBOUND(DT1,1),MElem,MElemGL,1 ,EGID(1),1,"DT1" )

     CALL PREAD(1,EL  ,LBOUND(EL,1),UBOUND(EL,1),NNode,NNodeGL,1 ,NGID,1,"EL"   )
     CALL PREAD(1,ET  ,LBOUND(ET,1),UBOUND(ET,1),NNode,NNodeGL,1 ,NGID,1,"ET"   )
     CALL PREAD(1,H   ,LBOUND(H,1), UBOUND(H,1), NNode,NNodeGL,1 ,NGID,1,"H"    )
     CALL PREAD(1,D   ,LBOUND(D,1), UBOUND(D,1), NNode,NNodeGL,1 ,NGID,1,"D"    )
     CALL PREAD(1,DT  ,LBOUND(DT,1),UBOUND(DT,1),NNode,NNodeGL,1 ,NGID,1,"DT"   )
     CLOSE(1)
#    endif
   END IF

   RETURN
   END SUBROUTINE HOT_START_SST 
!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
!   DUMP DATA FILE FOR ASSIMILATION RESTART                                    |
!==============================================================================|

   SUBROUTINE ARC_SST           

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   INTEGER I,K,ME,NPC
!==============================================================================|
   
   ME = MYID ; NPC = NPROCS 

   IF(MSR)THEN
     OPEN(1,FILE='restart_sst',FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(1)
     WRITE(1) IINT
   END IF

   IF(SERIAL)THEN
     WRITE(1) ((U(I,K),    K=1,KB),I=0,MElem)
     WRITE(1) ((V(I,K),    K=1,KB),I=0,MElem)
     WRITE(1) ((W(I,K),    K=1,KB),I=0,MElem)
     
     WRITE(1) ((Q2(I,K),   K=1,KB),I=0,NNode)
     WRITE(1) ((Q2L(I,K),  K=1,KB),I=0,NNode)
     
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
     WRITE(1) (( L(I,K),K=1,KB),I=1,NNode)

     WRITE(1) (UA(I), I=0,MElem)
     WRITE(1) (VA(I), I=0,MElem)

     WRITE(1) (EL1(I), I=1,MElem)
     WRITE(1) (ET1(I), I=1,MElem)
     WRITE(1) (H1(I),  I=1,MElem)
     WRITE(1) (D1(I),  I=1,MElem)
     WRITE(1) (DT1(I), I=1,MElem)

     WRITE(1) (EL(I), I=1,NNode)
     WRITE(1) (ET(I), I=1,NNode)
     WRITE(1) (H(I),  I=1,NNode)
     WRITE(1) (D(I),  I=1,NNode)
     WRITE(1) (DT(I), I=1,NNode)
   ELSE
#     if defined (MULTIPROCESSOR)
      CALL PWRITE(1,ME,NPC,U,    LBOUND(U,1),    UBOUND(U,1),    MElem,MElemGL,KB,EMAP,0,"U"    )
      CALL PWRITE(1,ME,NPC,V,    LBOUND(V,1),    UBOUND(V,1),    MElem,MElemGL,KB,EMAP,0,"V"    )
      CALL PWRITE(1,ME,NPC,W,    LBOUND(W,1),    UBOUND(W,1),    MElem,MElemGL,KB,EMAP,0,"W"    )
      
      CALL PWRITE(1,ME,NPC,Q2,   LBOUND(Q2,1),   UBOUND(Q2,1),   NNode,NNodeGL,KB,NMAP,0,"Q2"   )
      CALL PWRITE(1,ME,NPC,Q2L,  LBOUND(Q2L,1),  UBOUND(Q2L,1),  NNode,NNodeGL,KB,NMAP,0,"Q2L"  )
      
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
      CALL PWRITE(1,ME,NPC, L,LBOUND( L,1),UBOUND( L,1),NNode,NNodeGL,KB,NMAP,1,"L")

      CALL PWRITE(1,ME,NPC,UA,LBOUND(UA,1),UBOUND(UA,1),MElem,MElemGL,1,EMAP,0,"UA")
      CALL PWRITE(1,ME,NPC,VA,LBOUND(VA,1),UBOUND(VA,1),MElem,MElemGL,1,EMAP,0,"VA")

      CALL PWRITE(1,ME,NPC,EL1,LBOUND(EL1,1),UBOUND(EL1,1),MElem,MElemGL,1,EMAP,1,"EL1")
      CALL PWRITE(1,ME,NPC,ET1,LBOUND(ET1,1),UBOUND(ET1,1),MElem,MElemGL,1,EMAP,1,"ET1")
      CALL PWRITE(1,ME,NPC,H1, LBOUND(H1,1), UBOUND(H1,1), MElem,MElemGL,1,EMAP,1,"H1" )
      CALL PWRITE(1,ME,NPC,D1, LBOUND(D1,1), UBOUND(D1,1), MElem,MElemGL,1,EMAP,1,"D1" )
      CALL PWRITE(1,ME,NPC,DT1,LBOUND(DT1,1),UBOUND(DT1,1),MElem,MElemGL,1,EMAP,1,"DT1")

      CALL PWRITE(1,ME,NPC,EL,LBOUND(EL,1),UBOUND(EL,1),NNode,NNodeGL,1,NMAP,1,"EL")
      CALL PWRITE(1,ME,NPC,ET,LBOUND(ET,1),UBOUND(ET,1),NNode,NNodeGL,1,NMAP,1,"ET")
      CALL PWRITE(1,ME,NPC,H, LBOUND(H,1), UBOUND(H,1), NNode,NNodeGL,1,NMAP,1,"H" )
      CALL PWRITE(1,ME,NPC,D, LBOUND(D,1), UBOUND(D,1), NNode,NNodeGL,1,NMAP,1,"D" )
      CALL PWRITE(1,ME,NPC,DT,LBOUND(DT,1),UBOUND(DT,1),NNode,NNodeGL,1,NMAP,1,"DT")
#     endif
   END IF
   IF(MSR) CLOSE(1)

   RETURN
   END SUBROUTINE ARC_SST
!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|


   SUBROUTINE SET_TS_ASSIM_DATA 

!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR TEMP/SAL OBSERVATIONS                          |
!------------------------------------------------------------------------------!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)   
   USE MOD_PAR  
# endif
   IMPLICIT NONE
   INTEGER I,J,K,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   CHARACTER(LEN= 2 ) :: NAC   
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
   REAL(SP), PARAMETER :: ALF = 0.05_SP
   LOGICAL :: FEXIST
   INTEGER :: MAXEL,NBD_CNT
   INTEGER :: JMIN 
   REAL(SP):: LMIN

   REAL(SP), DIMENSION(1:MElemGL,1) :: RDLIST
   REAL(SP), DIMENSION(3) :: XTRI,YTRI
   REAL(SP) :: RDLAST
   INTEGER :: LOCIJ(2),MIN_LOC,JJ,IERR,Nsite_tmp
   INTEGER :: ND1,ND2,ND3
   REAL(SP) :: DELTA,COFA,COFB,COFC
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP
   REAL(SP), DIMENSION(KB) :: ZZ_OB
#  if defined (MULTIPROCESSOR)   
   REAL(SP), ALLOCATABLE :: ZZ_G(:,:)  
# endif
       
!------------------------------------------------------------------------------!
!  Read Number of Scalar Observations and Coordinates of Each                  !
!------------------------------------------------------------------------------!
       
   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_ts.xy"
!
!--Make Sure Temperature and Salinity Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMP/SALINITY OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number of T/S Measurement Stations---------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_TS                  !nomber of TS Obs station
   ALLOCATE(TS_OBS(N_ASSIM_TS))          !Type for TS_OBS

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!

   DO I=1,N_ASSIM_TS  
     READ(1,*)ITMP,TS_OBS(I)%X,TS_OBS(I)%Y,TS_OBS(I)%DEPTH,NLAY,TS_OBS(I)%SITA
     TS_OBS(I)%N_LAYERS = NLAY
     ALLOCATE(TS_OBS(I)%ODEPTH(NLAY))
     DO J=1,NLAY
       READ(1,*)TS_OBS(I)%ODEPTH(J)
       IF(TS_OBS(I)%ODEPTH(J) > TS_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'OBSERVATION DEPTH',J,'OF TEMP/SALINITY OBS',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH'
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO
   MAX_LAYER_TS = MAXVAL(TS_OBS(1:N_ASSIM_TS)%N_LAYERS)

!
!--Shift Coordinates-----------------------------------------------------------!
!
   TS_OBS(:)%X = TS_OBS(:)%X - VXMIN 
   TS_OBS(:)%Y = TS_OBS(:)%Y - VYMIN 
   
   IF(TS_METHOD == 'OI')THEN
!
!--find the cell number (TS_OBS(:)%N_CELL) of Obs station---------------------
!
   RRTMP = 100000.0     !100km
   DO J= 1,N_ASSIM_TS 
      X0 = TS_OBS(J)%X
      Y0 = TS_OBS(J)%Y
      DO I=1,MElemGL
         Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
	 if(Rtmp.LT.RRTMP) then
	    S11 = (XG(NVG(I,2))-X0)*(YG(NVG(I,3))-Y0)-&
	          (XG(NVG(I,3))-X0)*(YG(NVG(I,2))-Y0)
            S22 = (XG(NVG(I,3))-X0)*(YG(NVG(I,1))-Y0)-&
	          (XG(NVG(I,1))-X0)*(YG(NVG(I,3))-Y0)
            S33 = (XG(NVG(I,1))-X0)*(YG(NVG(I,2))-Y0)-&
	          (XG(NVG(I,2))-X0)*(YG(NVG(I,1))-Y0)
	    IF(S11.LE.0.AND.S22.LE.0.AND.S33.LE.0) THEN
	      TS_OBS(J)%N_CELL = I
	      GOTO 300
	    ELSE
	      TS_OBS(J)%N_CELL = 0
	    ENDIF
	 ELSE
	    TS_OBS(J)%N_CELL = -1
	 ENDIF
      ENDDO	
300 CONTINUE
      IF(TS_OBS(J)%N_CELL.LE.0) THEN
         IF(MSR) WRITE(IPT,*)'ERROR--T/S OBS SITE:',J,' OUT OF DOMAN',&
	         TS_OBS(J)%N_CELL 
         CALL PSTOP
      ENDIF 	 

   ENDDO             
!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   ALLOCATE(AW0G(MElemGL,3))
   ALLOCATE(AWXG(MElemGL,3))
   ALLOCATE(AWYG(MElemGL,3))
   
   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN     
     CALL GATHER(LBOUND(AW0,1),UBOUND(AW0,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AW0,AW0G)
     CALL GATHER(LBOUND(AWX,1),UBOUND(AWX,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AWX,AWXG)
     CALL GATHER(LBOUND(AWY,1),UBOUND(AWY,1),MElem,MElemGL,3,MYID,NPROCS,EMAP,AWY,AWYG)
   END IF
#  endif     
   END IF
!
!--Close Current Observation Global File---------------------------------------!
!
   CLOSE(1)
       
   IF(TS_METHOD == 'OI')THEN
!------------------------------------------------------------------------------!
!  Read Correlation Length of Scalar Observations                             !
!------------------------------------------------------------------------------!
       
!JQI   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_radius_ts.dat"

!JQI   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
!JQI   IF(MSR .AND. .NOT.FEXIST)THEN
!JQI     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
!JQI     WRITE(IPT,*)'HALTING.....'
!JQI     CALL PSTOP
!JQI   END IF
     
!JQI   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')

   ALLOCATE(PARAM_TS(2,N_ASSIM_TS))

!JQI   DO I=1,N_ASSIM_TS
!JQI     READ(1,*)PARAM_TS(1,I),PARAM_TS(2,I)
!JQI   END DO

!JQI   CLOSE(1)
     PARAM_TS = 30000.0_SP
   END IF
   
!------------------------------------------------------------------------------!
!  Open Temp/Sal Observation Files for Each Observation Point and Read Data    !
!------------------------------------------------------------------------------!

!----Make Sure Temperature/Salinity Observation File Exists--------------------!
   ONAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//'_ts.dat'
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMP/SALINITY OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!----Open Temp/Salinity Observation File for Read------------------------------!
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_TS
      READ(1,*,IOSTAT=IOS) nsite_tmp,TS_OBS(I)%N_TIMES
      IF(IOS<0) then
        WRITE(IPT,*) 'ERROR in read ',trim(casename),'_ts.dat at site number:',I
	CALL PSTOP
      ENDIF	

!----Allocate Arrays to Hold Temp/Salinity (TEMP/SAL) and Time (TIME)----------!
     ALLOCATE(TS_OBS(I)%TIMES(TS_OBS(I)%N_TIMES))
     ALLOCATE(TS_OBS(I)%TEMP( TS_OBS(I)%N_TIMES , TS_OBS(I)%N_LAYERS ))
     ALLOCATE(TS_OBS(I)%SAL(  TS_OBS(I)%N_TIMES , TS_OBS(I)%N_LAYERS ))

!----Read in Current Data for Observation I------------------------------------!
     NLAY = TS_OBS(I)%N_LAYERS
     DO J=1,TS_OBS(I)%N_TIMES
!       READ(1,*)TS_OBS(I)%TIMES(J),(TS_OBS(I)%TEMP(J,K),TS_OBS(I)%SAL(J,K),K=1,NLAY)
       READ(1,*,IOSTAT=IOS)TS_OBS(I)%TIMES(J),(TS_OBS(I)%TEMP(J,K),TS_OBS(I)%SAL(J,K),K=1,NLAY)
       IF(IOS < 0) then         ! ios=0 if all goes ok.
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_ts.dat at site:',I,  &
                     'Time No:',J
         CALL PSTOP
       ENDIF  
     END DO
     

!----Convert Time to Seconds---------------------------------------------------!
!----Shift Jan 1 Based Time Data to Dec 1 Based Time Data-----CASESPECIFIC-----!
     IF(trim(CASENAME) == 'gom')THEN
       TS_OBS(I)%TIMES = ((TS_OBS(I)%TIMES-1.0_SP)*24.0_SP+744.0_SP)*3600.0_SP
     ELSE   
       TS_OBS(I)%TIMES = TS_OBS(I)%TIMES*3600.0_SP*24.0_SP
     END IF
 
!----Convert Temperature and Salinity to PSU/Celsius-(If Necessary)------------!
     TS_OBS(I)%TEMP = TS_OBS(I)%TEMP 
     TS_OBS(I)%SAL  = TS_OBS(I)%SAL  
  END DO
  CLOSE(1)

!------------------------------------------------------------------------------!
!  Count Number of Points with Bad Data (TEMP = 0. .OR. SAL = 0.)         
!------------------------------------------------------------------------------!
  NBD_CNT = 0
  DO I=1,N_ASSIM_TS
    DO J=1,TS_OBS(I)%N_TIMES
      DO K=1,TS_OBS(I)%N_LAYERS
        IF(TS_OBS(I)%TEMP(J,K) < -90. .OR. TS_OBS(I)%SAL(J,K) < -90.) THEN
          NBD_CNT = NBD_CNT + 1
        END IF
      END DO
    END DO
  END DO
 
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Mooring Location 
!------------------------------------------------------------------------------!
!   LMIN = 100000000.
!   DO I=1,N_ASSIM_TS
!     X0 = TS_OBS(I)%X
!     Y0 = TS_OBS(I)%Y
!     DO J=1,NNodeGL
!       DX = ABS(XG(J)-X0)
!       DY = ABS(YG(J)-Y0)
!       IF(SQRT(DX**2 + DY**2) < LMIN)THEN
!         LMIN = SQRT(DX**2 + DY**2)
!         JMIN = J
!       END IF
!      END DO
!      TS_OBS(I)%SITA = SITA_GD(JMIN) + 3.14159_SP/2.0_SP
!    end do


   ALLOCATE(ITEMP(NNodeGL),FTEMP(NNodeGL),DA_TS(NNodeGL))     ; DA_TS = 0

   DO I=1,N_ASSIM_TS
     X0 = TS_OBS(I)%X
     Y0 = TS_OBS(I)%Y
     ISOBATH_ANGLE = TS_OBS(I)%SITA/180.0_SP*3.1415926_SP
     ECNT = 0
     DO J=1,NNodeGL
       DX = ABS(XG(J)-X0)
       DY = ABS(YG(J)-Y0)
       RD = SQRT(DX**2 + DY**2)
       IF(RD <= RAD_TS)THEN
         DA_TS(J)   = 1
         ECNT        = ECNT + 1      
         ITEMP(ECNT) =  J
         FTEMP(ECNT) = (RAD_TS**2 - RD**2) / (RAD_TS**2 + RD**2)
         ANG_OBS_SIM = ATAN2(DY,DX)
         D_ANGLE     = ANG_OBS_SIM - ISOBATH_ANGLE
         D_ANGLE     = D_ANGLE - INT(D_ANGLE/3.1415926_SP)*3.1415926_SP
         D_ANGLE     = ABS(D_ANGLE)
         DIR_WEIGHT  = (ABS(D_ANGLE-0.5*3.1415926_SP)+ALF*3.1415926_SP)/ &
                       ((0.5_SP+ALF)*3.1415926_SP)
         FTEMP(ECNT) = FTEMP(ECNT)*DIR_WEIGHT
       END IF
     END DO
     IF(ECNT == 0)THEN
       WRITE(IPT,*)'ERROR SETTING UP TEMP/SAL DATA ASSIMILATION'
       WRITE(IPT,*)'NO ELEMENTS LIE WITHIN RADIUS',RAD_TS
       WRITE(IPT,*)'OF OBSERVATION POINT',I
       CALL PSTOP   
     ELSE
       TS_OBS(I)%N_INTPTS = ECNT
       ALLOCATE(TS_OBS(I)%INTPTS(ECNT))
       ALLOCATE(TS_OBS(I)%X_WEIGHT(ECNT))
       TS_OBS(I)%INTPTS(1:ECNT)  = ITEMP(1:ECNT)
       TS_OBS(I)%X_WEIGHT(1:ECNT) = FTEMP(1:ECNT)
     END IF
   END DO
   DEALLOCATE(FTEMP,ITEMP)

!------------------------------------------------------------------------------!
!  Compute Sigma Layer Weights for Vertical Interpolation            
!------------------------------------------------------------------------------!
#  if defined (MULTIPROCESSOR)   
   ALLOCATE(ZZ_G(0:NNodeGL,KB))  
   IF(PAR)CALL GATHER(LBOUND(ZZ,1),UBOUND(ZZ,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,ZZ,ZZ_G)
   IF(PAR)CALL MPI_BCAST(ZZ_G,NNodeGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
#  endif

   DO I=1,N_ASSIM_TS
     NLAY = TS_OBS(I)%N_LAYERS
     ALLOCATE(TS_OBS(I)%S_INT(NLAY,2))
     ALLOCATE(TS_OBS(I)%S_WEIGHT(NLAY,2))
     
     X0 = TS_OBS(I)%X
     Y0 = TS_OBS(I)%Y
     RDLIST(1:MElemGL,1) = SQRT((XCG(1:MElemGL)-X0)**2 + (YCG(1:MElemGL)-Y0)**2)
     RDLAST = -1.0_SP
in:  DO WHILE(.TRUE.)
       LOCIJ = MINLOC(RDLIST,RDLIST>RDLAST)
       MIN_LOC = LOCIJ(1)
       IF(MIN_LOC == 0)THEN
         EXIT in
       END IF
       XTRI = XG(NVG(MIN_LOC,1:3))
       YTRI = YG(NVG(MIN_LOC,1:3))
       RDLAST = RDLIST(MIN_LOC,1)
       IF(ISINTRIANGLE1(XTRI,YTRI,X0,Y0))THEN
         JJ = MIN_LOC
	 EXIT IN
       END IF
       RDLAST = RDLIST(MIN_LOC,1)
     END DO IN  	 	      

     ND1=NVG(JJ,1)
     ND2=NVG(JJ,2)
     ND3=NVG(JJ,3)
     DELTA=(XG(ND2)-XG(ND1))*(YG(ND3)-YG(ND1))-     &
           (XG(ND3)-XG(ND1))*(YG(ND2)-YG(ND1))

     IF(SERIAL)THEN
     DO K=1,KBM1
       COFA=(YG(ND3)-YG(ND1))*(ZZ(ND2,K)-ZZ(ND1,K))-   &
            (YG(ND2)-YG(ND1))*(ZZ(ND3,K)-ZZ(ND1,K))
       COFB=(XG(ND2)-XG(ND1))*(ZZ(ND3,K)-ZZ(ND1,K))-   &
            (XG(ND3)-XG(ND1))*(ZZ(ND2,K)-ZZ(ND1,K))
       COFA=COFA/DELTA
       COFB=COFB/DELTA
       COFC=ZZ(ND1,K)-COFA*XG(ND1)-COFB*YG(ND1)
       ZZ_OB(K)=COFA*X0+COFB*Y0+COFC
     END DO  
     END IF

#  if defined (MULTIPROCESSOR)   
     IF(PAR)THEN
     DO K=1,KBM1
       COFA=(YG(ND3)-YG(ND1))*(ZZ_G(ND2,K)-ZZ_G(ND1,K))-   &
            (YG(ND2)-YG(ND1))*(ZZ_G(ND3,K)-ZZ_G(ND1,K))
       COFB=(XG(ND2)-XG(ND1))*(ZZ_G(ND3,K)-ZZ_G(ND1,K))-   &
            (XG(ND3)-XG(ND1))*(ZZ_G(ND2,K)-ZZ_G(ND1,K))
       COFA=COFA/DELTA
       COFB=COFB/DELTA
       COFC=ZZ_G(ND1,K)-COFA*XG(ND1)-COFB*YG(ND1)
       ZZ_OB(K)=COFA*X0+COFB*Y0+COFC
     END DO  
     END IF
#  endif
     
     DO J=1,NLAY
       SIGMA_C = -TS_OBS(I)%ODEPTH(J)/TS_OBS(I)%DEPTH
       DO K=2,KBM1
         IF(ZZ_OB(K) <= SIGMA_C .AND. ZZ_OB(K-1) > SIGMA_C)THEN 
           TS_OBS(I)%S_INT(J,1) = K-1
           TS_OBS(I)%S_INT(J,2) = K
           TS_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ_OB(K))/(ZZ_OB(K-1)-ZZ_OB(K))
           TS_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - TS_OBS(I)%S_WEIGHT(J,1) 
         END IF  
       END DO
       IF(ZZ_OB(1) <= SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
         TS_OBS(I)%S_INT(J,1) = 1
         TS_OBS(I)%S_INT(J,2) = 1
         TS_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         TS_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       IF(ZZ_OB(KBM1) > SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
         TS_OBS(I)%S_INT(J,1) = KBM1
         TS_OBS(I)%S_INT(J,2) = KBM1
         TS_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         TS_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF

     END DO
   END DO
#  if defined (MULTIPROCESSOR)   
   DEALLOCATE(ZZ_G)  
#  endif

!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN

   WRITE(IPT,*)
   WRITE(IPT,*)'!            TEMP/SALINITY OBSERVATION DATA           '
   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_NODE  SITA"
   DO I=1,N_ASSIM_TS
     MAXEL = MAXLOC(TS_OBS(I)%X_WEIGHT,DIM=1)
     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
     I,TS_OBS(I)%X/1000.,TS_OBS(I)%Y/1000., &
       TS_OBS(I)%N_INTPTS,TS_OBS(I)%N_TIMES,TS_OBS(I)%INTPTS(MAXEL),&
       TS_OBS(I)%SITA
   END DO
   WRITE(IPT,*)
   WRITE(IPT,*)'NUMBER OF BAD TS DATA POINTS: ',NBD_CNT
   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
   DO I=1,N_ASSIM_TS
   WRITE(IPT,*)I,TS_OBS(I)%TIMES(1)/(24.*3600.),&
       TS_OBS(I)%TIMES(TS_OBS(I)%N_TIMES)/(24.*3600.)
   END DO

   RETURN
   END SUBROUTINE SET_TS_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE TEMP_ASSIMILATION
!==============================================================================|
!  USE TEMP OBSERVATION DATA TO ADJUST TEMP FIELD                              |
!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR 
#  endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TINT,TCORR,TG,TCORR1
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TWGHT_T
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
   INTEGER I,J,K,J1,K1,K2,NLAY,ITIME,NTIME,IERR
   INTRINSIC MINLOC

   INTEGER IP,JN1,JN2,JN3
   REAL(SP) XSC,YSC,COFT0,COFTX,COFTY
!==============================================================================|

       
!------------------------------------------------------------------------------!
!  Gather T and S Fields to Master Processor                                   ! 
!------------------------------------------------------------------------------!
   ALLOCATE(TG(NNodeGL,KB))
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     CALL GATHER(LBOUND(TF1,1),UBOUND(TF1,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,TF1,TG)
   END IF
#  endif
   IF(SERIAL)THEN
     TG(1:NNodeGL,1:KBM1) = TF1(1:NNodeGL,1:KBM1)
   END IF
!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!

   IF(MSR)THEN
   TS_OBS%T_WEIGHT = 0. 
   T_THRESH         = ASTIME_WINDOW_TS    
   SIMTIME          = TIME*86400

   DO I=1,N_ASSIM_TS
     NTIME = TS_OBS(I)%N_TIMES
     ALLOCATE(FTEMP(NTIME)) 
     FTEMP(1:NTIME) = ABS(SIMTIME - TS_OBS(I)%TIMES(1:NTIME))
     DT_MIN = MINVAL(FTEMP(1:NTIME))
     TS_OBS(I)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)

     IF(DT_MIN < T_THRESH)THEN     
       IF(DT_MIN < .5_SP*T_THRESH) THEN
         TS_OBS(I)%T_WEIGHT = 1.0_SP
       ELSE
         TS_OBS(I)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
       END IF
     END IF

     DEALLOCATE(FTEMP)
   END DO
   
       
!------------------------------------------------------------------------------!
!  Interpolate Simulation Data to Local Observation Point                      ! 
!------------------------------------------------------------------------------!
       
   ALLOCATE(TINT(N_ASSIM_TS,MAX_LAYER_TS)) ; TINT = 0. 

   IF(TS_METHOD == 'NG')THEN
     DO I=1,N_ASSIM_TS   
       DO J=1,TS_OBS(I)%N_INTPTS
         J1        = TS_OBS(I)%INTPTS(J)
         WGHT      = TS_OBS(I)%X_WEIGHT(J)
         NLAY      = TS_OBS(I)%N_LAYERS
         DO K=1,NLAY
           U1 = TG(J1,TS_OBS(I)%S_INT(K,1))
           U2 = TG(J1,TS_OBS(I)%S_INT(K,2))
           W1 = TS_OBS(I)%S_WEIGHT(K,1)
           W2 = TS_OBS(I)%S_WEIGHT(K,2)
           TINT(I,K) = TINT(I,K) + (U1*W1 + U2*W2)*WGHT 
         END DO
       END DO
       TOT_WGHT = SUM(TS_OBS(I)%X_WEIGHT(1:TS_OBS(I)%N_INTPTS))
       TINT(I,1:NLAY) = TINT(I,1:NLAY)/TOT_WGHT
     END DO

   ELSE IF(TS_METHOD == 'OI')THEN
     DO I=1,N_ASSIM_TS
       IP = TS_OBS(I)%N_CELL
       NLAY      = TS_OBS(I)%N_LAYERS
       XSC = TS_OBS(I)%X-XCG(IP)
       YSC = TS_OBS(I)%Y-YCG(IP)
       JN1 = NVG(IP,1)
       JN2 = NVG(IP,2)
       JN3 = NVG(IP,3)
      
       DO K=1,NLAY
         K1 = TS_OBS(I)%S_INT(K,1)
         K2 = TS_OBS(I)%S_INT(K,2)
         COFT0 = AW0G(IP,1)*TG(JN1,K1)+AW0G(IP,2)*TG(JN2,K1)+AW0G(IP,3)*TG(JN3,K1)
         COFTX = AWXG(IP,1)*TG(JN1,K1)+AWXG(IP,2)*TG(JN2,K1)+AWXG(IP,3)*TG(JN3,K1)
         COFTY = AWYG(IP,1)*TG(JN1,K1)+AWYG(IP,2)*TG(JN2,K1)+AWYG(IP,3)*TG(JN3,K1)
         U1 = COFT0 + COFTX*XSC + COFTY*YSC
         COFT0 = AW0G(IP,1)*TG(JN1,K2)+AW0G(IP,2)*TG(JN2,K2)+AW0G(IP,3)*TG(JN3,K2)
         COFTX = AWXG(IP,1)*TG(JN1,K2)+AWXG(IP,2)*TG(JN2,K2)+AWXG(IP,3)*TG(JN3,K2)
         COFTY = AWYG(IP,1)*TG(JN1,K2)+AWYG(IP,2)*TG(JN2,K2)+AWYG(IP,3)*TG(JN3,K2)
         U2 = COFT0 + COFTX*XSC + COFTY*YSC
         W1 = TS_OBS(I)%S_WEIGHT(K,1)
         W2 = TS_OBS(I)%S_WEIGHT(K,2)
	 TINT(I,K) = (U1*W1 + U2*W2)
       ENDDO
     ENDDO   	 
   ELSE
     WRITE(IPT,*)'TS_METHOD SHOULD BE NG OR OI. BUT HERE TS_METHOD=',TRIM(TS_METHOD)
     CALL PSTOP
   END IF  
!------------------------------------------------------------------------------!
!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!------------------------------------------------------------------------------!

   ALLOCATE(TWGHT_T(NNodeGL,KBM1))   ; TWGHT_T   = 0.
   ALLOCATE(TCORR(NNodeGL,KBM1))   ; TCORR   = 0.

   IF(TS_METHOD == 'NG')THEN
     DO I=1,N_ASSIM_TS 
       DO J=1,TS_OBS(I)%N_INTPTS
         J1     = TS_OBS(I)%INTPTS(J)
         ITIME  = TS_OBS(I)%N_T_WEIGHT
         NLAY   = TS_OBS(I)%N_LAYERS
         DO K=1,NLAY
           K1             = TS_OBS(I)%S_INT(K,1)
           K2             = TS_OBS(I)%S_INT(K,2)
           W1             = TS_OBS(I)%S_WEIGHT(K,1)
           W2             = TS_OBS(I)%S_WEIGHT(K,2)

           DEFECT         = TS_OBS(I)%TEMP(ITIME,K) - TINT(I,K)
           IF(ABS(DEFECT) < 20.0)THEN     !quality control
             WEIGHT1        = TS_OBS(I)%T_WEIGHT*TS_OBS(I)%X_WEIGHT(J)*W1
             WEIGHT2        = TS_OBS(I)%T_WEIGHT*TS_OBS(I)%X_WEIGHT(J)*W2
             TWGHT_T(J1,K1) = TWGHT_T(J1,K1) + WEIGHT1   
             TWGHT_T(J1,K2) = TWGHT_T(J1,K2) + WEIGHT2   

             CORRECTION     = GAMA_TS*DEFECT
!           TCORR(J1,K1)   = TCORR(J1,K1) + CORRECTION*WEIGHT1
!           TCORR(J1,K2)   = TCORR(J1,K2) + CORRECTION*WEIGHT2
             TCORR(J1,K1)   = TCORR(J1,K1) + CORRECTION*WEIGHT1**2
             TCORR(J1,K2)   = TCORR(J1,K2) + CORRECTION*WEIGHT2**2
           ENDIF
         END DO
       END DO
     END DO

!------------------------------------------------------------------------------!
!  Nudge Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     DO I=1,NNodeGL
       DO K=1,KBM1
         IF(DA_TS(I) == 1 .AND. TWGHT_T(I,K) > 1.0E-08)THEN
           TG(I,K) = TG(I,K) + DTI*GALPHA_TS*TCORR(I,K)/TWGHT_T(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT_T,TCORR,TINT)
   ELSE IF(TS_METHOD == 'OI')THEN
     DO I=1,N_ASSIM_TS 
       ITIME  = TS_OBS(I)%N_T_WEIGHT
       NLAY   = TS_OBS(I)%N_LAYERS
       DO K=1,NLAY
         K1             = TS_OBS(I)%S_INT(K,1)
         K2             = TS_OBS(I)%S_INT(K,2)
         W1             = TS_OBS(I)%S_WEIGHT(K,1)
         W2             = TS_OBS(I)%S_WEIGHT(K,2)

         DEFECT         = TS_OBS(I)%TEMP(ITIME,K) - TINT(I,K)
         IF(ABS(DEFECT) < 20.0)THEN     !quality control
           WEIGHT1        = TS_OBS(I)%T_WEIGHT*W1
           WEIGHT2        = TS_OBS(I)%T_WEIGHT*W2
           TWGHT_T(I,K1) = TWGHT_T(I,K1) + WEIGHT1   
           TWGHT_T(I,K2) = TWGHT_T(I,K2) + WEIGHT2   
           CORRECTION     = DEFECT
           TCORR(I,K1)   = TCORR(I,K1) + CORRECTION*WEIGHT1
           TCORR(I,K2)   = TCORR(I,K2) + CORRECTION*WEIGHT2
         ENDIF
       END DO
     END DO

     DO I=1,N_ASSIM_TS
       DO K=1,KBM1
         IF(TWGHT_T(I,K) > 1.0E-8)THEN
           TCORR(I,K)=TCORR(I,K)/TWGHT_T(I,K)
         END IF
       END DO
     END DO
       	 
!------------------------------------------------------------------------------!
!  'OI' Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     ALLOCATE(TCORR1(NNodeGL,KBM1));TCORR1=0.0_SP

     DO K=1,KBM1
       CALL TS_OPTIMINTERP(TCORR(:,K),TCORR1(:,K))
     END DO

     DO I=1,NNodeGL
       DO K=1,KBM1
         IF(DA_TS(I) == 1)THEN
!QXU         TG(I,K) = TG(I,K) + TCORR1(I,K)
           TG(I,K) = TG(I,K) + GALPHA_TS*TCORR1(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT_T,TCORR,TINT,TCORR1)
   ELSE
     WRITE(IPT,*)'TS_METHOD SHOULD BE NG OR OI. BUT HERE TS_METHOD=',TRIM(TS_METHOD)
     CALL PSTOP
   END IF

   END IF  !!MASTER


!------------------------------------------------------------------------------!
!  Disperse New Data Fields to Slave Processors                                ! 
!------------------------------------------------------------------------------!
   IF(SERIAL)THEN
     TF1(1:NNode,1:KBM1) = TG(1:NNode,1:KBM1)
   END IF
#  if defined (MULTIPROCESSOR) 
   CALL MPI_BCAST(TG,NNodeGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)THEN
     DO I=1,NNode
       TF1(I,1:KBM1) = TG(NGID(I),1:KBM1)
     END DO
   END IF
#  endif

   DEALLOCATE(TG)
          

   RETURN
   END SUBROUTINE TEMP_ASSIMILATION
!==============================================================================|
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE SALT_ASSIMILATION
!==============================================================================|
!  USE SALT OBSERVATION DATA TO ADJUST SALINITY FIELD                          |
!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR 
#  endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: SINT,SCORR,SG,SCORR1
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TWGHT_S
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
   INTEGER I,J,K,J1,K1,K2,NLAY,ITIME,NTIME,IERR
   INTRINSIC MINLOC

   INTEGER IP,JN1,JN2,JN3
   REAL(SP) XSC,YSC,COFT0,COFTX,COFTY
!==============================================================================|

       
!------------------------------------------------------------------------------!
!  Gather S Fields to Master Processor                                         ! 
!------------------------------------------------------------------------------!
   ALLOCATE(SG(NNodeGL,KB))
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     CALL GATHER(LBOUND(SF1,1),UBOUND(SF1,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,SF1,SG)
   END IF
#  endif
   IF(SERIAL)THEN
     SG(1:NNodeGL,1:KBM1) = SF1(1:NNodeGL,1:KBM1)
   END IF
!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!

   IF(MSR)THEN
   TS_OBS%T_WEIGHT = 0. 
   T_THRESH         = ASTIME_WINDOW_TS    
   SIMTIME          = TIME*86400

   DO I=1,N_ASSIM_TS
     NTIME = TS_OBS(I)%N_TIMES
     ALLOCATE(FTEMP(NTIME)) 
     FTEMP(1:NTIME) = ABS(SIMTIME - TS_OBS(I)%TIMES(1:NTIME))
     DT_MIN = MINVAL(FTEMP(1:NTIME))
     TS_OBS(I)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)

     IF(DT_MIN < T_THRESH)THEN     
       IF(DT_MIN < .5_SP*T_THRESH) THEN
         TS_OBS(I)%T_WEIGHT = 1.0_SP
       ELSE
         TS_OBS(I)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
       END IF
     END IF

     DEALLOCATE(FTEMP)
   END DO
   
       
!------------------------------------------------------------------------------!
!  Interpolate Simulation Data to Local Observation Point                      ! 
!------------------------------------------------------------------------------!
       
   ALLOCATE(SINT(N_ASSIM_TS,MAX_LAYER_TS)) ; SINT = 0.

   IF(TS_METHOD == 'NG')THEN
     DO I=1,N_ASSIM_TS   
       DO J=1,TS_OBS(I)%N_INTPTS
         J1        = TS_OBS(I)%INTPTS(J)
         WGHT      = TS_OBS(I)%X_WEIGHT(J)
         NLAY      = TS_OBS(I)%N_LAYERS
         DO K=1,NLAY
           V1 = SG(J1,TS_OBS(I)%S_INT(K,1))
           V2 = SG(J1,TS_OBS(I)%S_INT(K,2))
           W1 = TS_OBS(I)%S_WEIGHT(K,1)
           W2 = TS_OBS(I)%S_WEIGHT(K,2)
           SINT(I,K) = SINT(I,K) + (V1*W1 + V2*W2)*WGHT 
         END DO
       END DO
       TOT_WGHT = SUM(TS_OBS(I)%X_WEIGHT(1:TS_OBS(I)%N_INTPTS))
       SINT(I,1:NLAY) = SINT(I,1:NLAY)/TOT_WGHT
     END DO
   ELSE IF(TS_METHOD == 'OI')THEN
     DO I=1,N_ASSIM_TS
       IP = TS_OBS(I)%N_CELL
       NLAY      = TS_OBS(I)%N_LAYERS
       XSC = TS_OBS(I)%X-XCG(IP)
       YSC = TS_OBS(I)%Y-YCG(IP)
       JN1 = NVG(IP,1)
       JN2 = NVG(IP,2)
       JN3 = NVG(IP,3)
      
       DO K=1,NLAY
         K1 = TS_OBS(I)%S_INT(K,1)
         K2 = TS_OBS(I)%S_INT(K,2)
         COFT0 = AW0G(IP,1)*SG(JN1,K1)+AW0G(IP,2)*SG(JN2,K1)+AW0G(IP,3)*SG(JN3,K1)
         COFTX = AWXG(IP,1)*SG(JN1,K1)+AWXG(IP,2)*SG(JN2,K1)+AWXG(IP,3)*SG(JN3,K1)
         COFTY = AWYG(IP,1)*SG(JN1,K1)+AWYG(IP,2)*SG(JN2,K1)+AWYG(IP,3)*SG(JN3,K1)
         U1 = COFT0 + COFTX*XSC + COFTY*YSC
         COFT0 = AW0G(IP,1)*SG(JN1,K2)+AW0G(IP,2)*SG(JN2,K2)+AW0G(IP,3)*SG(JN3,K2)
         COFTX = AWXG(IP,1)*SG(JN1,K2)+AWXG(IP,2)*SG(JN2,K2)+AWXG(IP,3)*SG(JN3,K2)
         COFTY = AWYG(IP,1)*SG(JN1,K2)+AWYG(IP,2)*SG(JN2,K2)+AWYG(IP,3)*SG(JN3,K2)
         U2 = COFT0 + COFTX*XSC + COFTY*YSC
         W1 = TS_OBS(I)%S_WEIGHT(K,1)
         W2 = TS_OBS(I)%S_WEIGHT(K,2)
	 SINT(I,K) = (U1*W1 + U2*W2)
       ENDDO
     ENDDO   	 
   ELSE
     WRITE(IPT,*)'TS_METHOD SHOULD BE NG OR OI. BUT HERE TS_METHOD=',TRIM(TS_METHOD)
     CALL PSTOP
   END IF  
!------------------------------------------------------------------------------!
!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!------------------------------------------------------------------------------!

   ALLOCATE(TWGHT_S(NNodeGL,KBM1))   ; TWGHT_S   = 0.
   ALLOCATE(SCORR(NNodeGL,KBM1))   ; SCORR   = 0.

   IF(TS_METHOD == 'NG')THEN
     DO I=1,N_ASSIM_TS 
       DO J=1,TS_OBS(I)%N_INTPTS
         J1     = TS_OBS(I)%INTPTS(J)
         ITIME  = TS_OBS(I)%N_T_WEIGHT
         NLAY   = TS_OBS(I)%N_LAYERS
         DO K=1,NLAY
           K1             = TS_OBS(I)%S_INT(K,1)
           K2             = TS_OBS(I)%S_INT(K,2)
           W1             = TS_OBS(I)%S_WEIGHT(K,1)
           W2             = TS_OBS(I)%S_WEIGHT(K,2)

           DEFECT         = TS_OBS(I)%SAL(ITIME,K) - SINT(I,K)
	   IF(ABS(DEFECT).LT.20.0) THEN             !quality control

             WEIGHT1        = TS_OBS(I)%T_WEIGHT*TS_OBS(I)%X_WEIGHT(J)*W1
             WEIGHT2        = TS_OBS(I)%T_WEIGHT*TS_OBS(I)%X_WEIGHT(J)*W2
             TWGHT_S(J1,K1) = TWGHT_S(J1,K1) + WEIGHT1   
             TWGHT_S(J1,K2) = TWGHT_S(J1,K2) + WEIGHT2   

             CORRECTION     = GAMA_TS*DEFECT

!          SCORR(J1,K1)   = SCORR(J1,K1) + CORRECTION*WEIGHT1
!          SCORR(J1,K2)   = SCORR(J1,K2) + CORRECTION*WEIGHT2
             SCORR(J1,K1)   = SCORR(J1,K1) + CORRECTION*WEIGHT1**2
             SCORR(J1,K2)   = SCORR(J1,K2) + CORRECTION*WEIGHT2**2

           ENDIF
         END DO
       END DO
     END DO

!------------------------------------------------------------------------------!
!  Nudge Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     DO I=1,NNodeGL
       DO K=1,KBM1
         IF(DA_TS(I) == 1 .AND. TWGHT_S(I,K) > 1.0E-08)THEN
           SG(I,K) = SG(I,K) + DTI*GALPHA_TS*SCORR(I,K)/TWGHT_S(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT_S,SCORR,SINT)
   ELSE IF(TS_METHOD == 'OI')THEN
     DO I=1,N_ASSIM_TS 
       ITIME  = TS_OBS(I)%N_T_WEIGHT
       NLAY   = TS_OBS(I)%N_LAYERS
       DO K=1,NLAY
         K1             = TS_OBS(I)%S_INT(K,1)
         K2             = TS_OBS(I)%S_INT(K,2)
         W1             = TS_OBS(I)%S_WEIGHT(K,1)
         W2             = TS_OBS(I)%S_WEIGHT(K,2)

         DEFECT         = TS_OBS(I)%SAL(ITIME,K) - SINT(I,K)
         IF(ABS(DEFECT).LT.20.0) THEN             !quality control
           WEIGHT1        = TS_OBS(I)%T_WEIGHT*W1
           WEIGHT2        = TS_OBS(I)%T_WEIGHT*W2
           TWGHT_S(I,K1) = TWGHT_S(I,K1) + WEIGHT1   
           TWGHT_S(I,K2) = TWGHT_S(I,K2) + WEIGHT2   
           CORRECTION     = DEFECT
           SCORR(I,K1)   = SCORR(I,K1) + CORRECTION*WEIGHT1
           SCORR(I,K2)   = SCORR(I,K2) + CORRECTION*WEIGHT2
         ENDIF
       END DO
     END DO

     DO I=1,N_ASSIM_TS
       DO K=1,KBM1
         IF(TWGHT_S(I,K) > 1.0E-8)THEN
           SCORR(I,K)=SCORR(I,K)/TWGHT_S(I,K)
         END IF
       END DO
     END DO

!------------------------------------------------------------------------------!
!  'OI' Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     ALLOCATE(SCORR1(NNodeGL,KBM1));SCORR1=0.0_SP

     DO K=1,KBM1
       CALL TS_OPTIMINTERP(SCORR(:,K),SCORR1(:,K))
     END DO
   
     DO I=1,NNodeGL
       DO K=1,KBM1
         IF(DA_TS(I) == 1)THEN
!QXU       SG(I,K)= SG(I,K) + SCORR1(I,K)
           SG(I,K) = SG(I,K) + GALPHA_TS*SCORR1(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT_S,SCORR,SINT,SCORR1)
   ELSE
     WRITE(IPT,*)'TS_METHOD SHOULD BE NG OR OI. BUT HERE TS_METHOD=',TRIM(TS_METHOD)
     CALL PSTOP
   END IF

   END IF  !!MASTER
   
!------------------------------------------------------------------------------!
!  Disperse New Data Fields to Slave Processors                                ! 
!------------------------------------------------------------------------------!
   IF(SERIAL)THEN
     SF1(1:NNode,1:KBM1) = SG(1:NNode,1:KBM1)
   END IF
#  if defined (MULTIPROCESSOR) 
   CALL MPI_BCAST(SG,NNodeGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)THEN
     DO I=1,NNode
#  if defined (WET_DRY)
       IF(ISWET_NODE_LASTSTEP(I)*ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif       
       SF1(I,1:KBM1) = SG(NGID(I),1:KBM1)
#  if defined (WET_DRY)
       END IF
#  endif              
     END DO
   END IF
#  endif

   DEALLOCATE(SG)
          
   RETURN
   END SUBROUTINE SALT_ASSIMILATION
!==============================================================================|
!==============================================================================|

!==============================================================================|
!==============================================================================|
   SUBROUTINE TS_OPTIMINTERP(F,FI)
   USE MOD_OPTIMAL_INTERPOLATION
   USE MOD_PREC
   USE ALL_VARS
   IMPLICIT NONE

!------------------------------------------------------------------------------|
! xi(1,:) and xi(2,:) represent the x and y coordindate of the grid of the     |
! interpolated field                                                           |
! fi and vari are the interpolated field and its error variance resp.          |
!------------------------------------------------------------------------------|
   REAL(SP) :: XI(2,NNodeGL),FI(NNodeGL),VARI(NNodeGL)

!------------------------------------------------------------------------------|
! x(1,:) and x(2,:) represent the x and y coordindate of the observations      |
! f and var are observations and their error variance resp.                    |
!------------------------------------------------------------------------------|
   REAL(SP) :: X(2,N_ASSIM_TS),VAR(N_ASSIM_TS),F(N_ASSIM_TS)

!------------------------------------------------------------------------------|
! param: inverse of the correlation length                                     |
!------------------------------------------------------------------------------|
   REAL(SP) :: PARAM(2,N_ASSIM_TS)

   INTEGER  :: I,J,MM   

!------------------------------------------------------------------------------|
! create a regular 2D grid                                                     |
!------------------------------------------------------------------------------|
   DO I=1,NNodeGL
     XI(1,I) = XG(I)
     XI(2,I) = YG(I)
   END DO	

!------------------------------------------------------------------------------|   
! param is the inverse of the correlation length                               |
!------------------------------------------------------------------------------|
   PARAM = 1.0_SP/PARAM_TS             
 
   MM = N_INFLU_TS

!------------------------------------------------------------------------------|   
! the error variance of the observations                                       |
!------------------------------------------------------------------------------|
   VAR = 0.0_SP   

!------------------------------------------------------------------------------|
! location of observations                                                     |
!------------------------------------------------------------------------------|
   DO I=1,N_ASSIM_TS
     X(1,I) = TS_OBS(I)%X
     X(2,I) = TS_OBS(I)%Y
   END DO  

!------------------------------------------------------------------------------|
! fi is the interpolated function and vari its error variance                  |
!------------------------------------------------------------------------------|
   CALL OPTIMINTERP(X,F,VAR,PARAM,MM,XI,FI,VARI)

   RETURN
   END SUBROUTINE TS_OPTIMINTERP
!==============================================================================|
!==============================================================================|

!==============================================================================|
   LOGICAL FUNCTION ISINTRIANGLE1(XT,YT,X0,Y0) 
!==============================================================================|
!  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    |
!  using algorithm used for scene rendering in computer graphics               |
!  algorithm works well unless particle happens to lie in a line parallel      |
!  to the edge of a triangle.                                                  |
!  This can cause problems if you use a regular grid, say for idealized        |
!  modelling and you happen to see particles right on edges or parallel to     |
!  edges.                                                                      |
!==============================================================================|

   USE MOD_PREC
   IMPLICIT NONE
   REAL(SP), INTENT(IN) :: X0,Y0
   REAL(SP), INTENT(IN) :: XT(3),YT(3)
   REAL(SP) :: F1,F2,F3
   REAL(SP) :: X1(2)
   REAL(SP) :: X2(2)
   REAL(SP) :: X3(2)
   REAL(SP) :: P(2)

!------------------------------------------------------------------------------|

   ISINTRIANGLE1 = .FALSE.  

   IF(Y0 < MINVAL(YT) .OR. Y0 > MAXVAL(YT))THEN
     ISINTRIANGLE1 = .FALSE.
     RETURN
   END IF
   IF(X0 < MINVAL(XT) .OR. X0 > MAXVAL(XT))THEN
     ISINTRIANGLE1 = .FALSE.
     RETURN
   END IF

   F1 = (Y0-YT(1))*(XT(2)-XT(1)) - (X0-XT(1))*(YT(2)-YT(1))
   F2 = (Y0-YT(3))*(XT(1)-XT(3)) - (X0-XT(3))*(YT(1)-YT(3))
   F3 = (Y0-YT(2))*(XT(3)-XT(2)) - (X0-XT(2))*(YT(3)-YT(2))
   IF(F1*F3 >= 0.0_SP .AND. F3*F2 >= 0.0_SP) ISINTRIANGLE1 = .TRUE.

   RETURN
   END FUNCTION ISINTRIANGLE1
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
#  endif
   END MODULE MOD_ASSIM
