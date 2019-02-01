!==============================================================================|
!   GLOBAL LIMITS AND ARRAY SIZING PARAMETERS                                  !
!==============================================================================|

MODULE LIMS
   USE MOD_PREC
   IMPLICIT NONE
   SAVE
 
   INTEGER MElemGL                !!GLOBAL NUMBER OF ELEMENTS
   INTEGER NNodeGL                !!GLOBAL NUMBER OF NODES
   INTEGER IBFW_GL            !!GLOBAL NUMBER OF GROUNDWATER NODES
   INTEGER NUMQBC_GL          !!GLOBAL NUMBER OF FRESHWATER INFLOW NODES
   INTEGER NOBCGEO_GL         !!GLOBAL NUMBER OF OPEN BOUNDARY 
   INTEGER NDRFT_GL           !!GLOBAL NUMBER OF LAGRANGIAN TRACKING PARTICLES
   INTEGER IOBCN_RGL          !!T.W., added to save open boundary node#, per RGL
 

   INTEGER MElem                  !!LOCAL NUMBER OF ELEMENTS 
   INTEGER NNode                  !!LOCAL NUMBER OF NODES
   INTEGER IBFW               !!LOCAL NUMBER OF GROUNDWATER NODES
   INTEGER NUMQBC             !!LOCAL NUMBER OF FRESHWATER INFLOW NODES
   INTEGER NOBCGEO            !!LOCAL NUMBER OF OPEN BOUNDARY 
   INTEGER NDRFT              !!LOCAL NUMBER OF LAGRANGIAN TRACKING PARTICLES
   INTEGER NISBCE_1           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 1
   INTEGER NISBCE_2           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 2
   INTEGER NISBCE_3           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 3

   INTEGER KB                 !!NUMBER OF SIGMA LEVELS
   INTEGER KBM1               !!NUMBER OF SIGMA LEVELS-1
   INTEGER KBM2               !!NUMBER OF SIGMA LEVELS-2
   INTEGER MYID               !!UNIQUE PROCESSOR ID (1 => NPROCS)
   INTEGER KSL                !!NUMBER OF STANDARD SEA LEVELS 
   INTEGER NPROCS             !!NUMBER OF PROCESSORS
   INTEGER NE                 !!NUMBER OF UNIQUE EDGES (LOCAL DOMAIN ONLY)
   INTEGER NCV                !!NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)
   
   INTEGER NCV_I              !!NUMBER OF INTERNAL CONTROL VOLUMES (LOCAL ONLY)
   INTEGER MTElem                 !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
   INTEGER NTNode                 !!TOTAL OF LOCAL INTERNAL + HALO NODES
   INTEGER MX_NBR_ELEM        !!MAX NUMBER OF ELEMENTS SURROUNDING A NODE

END MODULE LIMS



!==============================================================================|
!   CONTROL VARIABLES                                                          |
!==============================================================================|

MODULE CONTROL
   USE MOD_PREC
   IMPLICIT NONE
   SAVE

   LOGICAL SERIAL                  !!TRUE IF SINGLE PROCESSOR
   LOGICAL MSR                     !!TRUE IF MASTER PROCESSOR (MYID==1)
   LOGICAL PAR                     !!TRUE IF MULTIPROCESSOR RUN
   LOGICAL JMPOBC                  !!TRUE IF USING JS MODIFIED OPEN BOUNDARY CONDITION 
   LOGICAL EVP_FLAG                !!TRUE IF INPUT EVAPOTATION AND PRECIPITATION 
   LOGICAL :: WIND_BLOCK           ! wind block logical, read from run file input
   
   
   CHARACTER(LEN=80) FVCOM_VERSION !!STRING DESCRIBING VERSION
   CHARACTER(LEN=80) FVCOM_WEBSITE !!STRING DESCRIBING WEBSITE FOR FVCOM INFO 
   CHARACTER(LEN=80) CASENAME      !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
   CHARACTER(LEN=120) CASETITLE    !!CASE TITLE                                 
   CHARACTER(LEN=80) OUTDIR        !!PARENT OUTPUT DIRECTORY
   CHARACTER(LEN=80) INPDIR        !!MAIN   INPUT DIRECTORY
   CHARACTER(LEN=80) INFOFILE      !!INFO  FILE            
   CHARACTER(LEN=80) M_TYPE        !!METEOROLOGICAL FORCING TYPE (uniform/non-uniform)
   CHARACTER(LEN=80) S_TYPE        !!(julian/non-julian)  
   CHARACTER(LEN=80) BROUGH_TYPE   !!(orig/gotm/user_defined)  
   CHARACTER(LEN=80) WINDTYPE      !!WIND FORCING (stress/speed)  
   CHARACTER(LEN=80) INFLOW_TYPE   !!SPECIFIED RIVER INFLOW TYPE (edge/node) 
   CHARACTER(LEN=80) POINT_ST_TYPE !!(calculated/specified)
   CHARACTER(LEN=80) CTRL_DEN      !!CONTROLS DENSITY CALCULATION (sigma-t/pdensity/sigma-t_stp)
   CHARACTER(LEN=80) H_TYPE        !!SURFACE HEATING (body_h/flux_h)
   CHARACTER(LEN=80) C_BAROPG      !!CONTROL BAROCLINIC PRES GRADIENT CALC (sigma/s_levels) 
   CHARACTER(LEN=80) HORZMIX       !!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
   CHARACTER(LEN=80) VERTMIX       !!CONTROLS VERTICAL DIFFUSION COEF CALC [UNUSED]
   CHARACTER(LEN=80) RESTART       !!CONTROLS RESTART TYPE (cold_start/hot_cold_s/hot_start)


!--Constants-------------------------------------------------------------------!
   REAL(SP), PARAMETER, DIMENSION(4) :: ALPHA_RK = (/0.2500_SP,0.333333_SP,0.5000_SP,1.0_SP/)
   REAL(SP), PARAMETER :: GRAV      = 9.81_SP
   REAL(SP), PARAMETER :: PI        = 3.141592653_SP
   REAL(SP), PARAMETER :: PI2       = 6.283185307_SP
   REAL(SP), PARAMETER :: ZERO      = 0.0_SP 
   REAL(SP), PARAMETER :: ONE_THIRD = 1.0_SP/3.0_SP 
   REAL(SP), PARAMETER :: REARTH    = 6371.0E03_SP   !!Earth Radius in Meters
   REAL(SP), PARAMETER :: DEG2RAD   = PI2/360.0_SP   !!Radians/Degree
   REAL(SP), PARAMETER :: TPI       = DEG2RAD*REARTH !TPI=pi*rearth/180.=3.14159265/180.0*6371.*1000.
   REAL(SP), PARAMETER :: ROFVROS   = 0.9775171065_SP!!RATIO OF THE DENSITY OF FRESH AND SEA WATER 1000./1023.   

   REAL(SP), PARAMETER :: SLP0      = 101325.0_SP    !! mean sea surface pressure (Pa)

!--Gravity (vary with latitute) -----------------------------------------------!
   REAL(SP),ALLOCATABLE :: GRAV_N(:),GRAV_E(:)
   
!--Tidal Periods------------------:    
   REAL(SP), PARAMETER, DIMENSION(8) :: &
                    !     s2         m2           n2          k2          k1          p1         o1         q1
           PERIOD = (/43200.0_SP, 44712.0_SP, 45570.0_SP, 43082.0_SP, 86164.0_SP, 86637.0_SP, 92950.0_SP, 96726.0_SP/)

!--File Unit Specifiers (see iofiles.f for details)----------------------------!
   INTEGER  INRUN,INDEP,INGRD,INOBC,INCWH,INRIV,INITS,INOEL, &
            INJUL,INHFX,INWND,INEVP,INELF,INUVF,INLAG,INSPO, &
            INCOR,INBFW,IOPRT,IOPLT,IOTSR, &
            IOSMSD,IOSMSV,IOSMST,INRES,INJMP,IREST,IPT,INRIVW,IGOTM, &
			INWNDBLK  ! Wind block file added by B Clark
			

!--Parameter Controlling Vertical Coordinate Distribution----------------------!
   INTEGER  :: INDEX_VERCOR 
!--Sigma Level Parameters------------------------------------------------------!
   REAL(SP) :: P_SIGMA      !!PARAMETER CONTROLLING SIGMA LEVEL DISTRIBUTION
   REAL(SP):: HMIN2
!--General Vertical Level Parameters for case 2--------------------------------!
   REAL(SP) :: DU2          !!PARAMETER CONTROLLING LEVEL DISTRIBUTION OF SURFACE
   REAL(SP) :: DL2          !!PARAMETER CONTROLLING LEVEL DISTRIBUTION OF BOTTOM

!--General Vertical Level Parameters for case 3--------------------------------!
   REAL(SP) :: DUU          !!THE UPPER BOUNDARY OF PARALLEL COORDINATE
   REAL(SP) :: DLL          !!THE LOWER BOUNDARY OF PARALLEL COORDINATE
   REAL(SP) :: HMIN1        !!THE MIN DEPTH AT WHICH THE LAYERS ARE CONSTANT
   
   INTEGER  :: KU           !!THE NUMBERS OF LAYERS ABOVE UPPER BOUNDARY
   INTEGER  :: KL           !!THE NUMBVER OF LAYERS BELOW LOWER BOUNDARY
   
   REAL(SP), ALLOCATABLE :: ZKU(:)       !!THE DEPTHS OF PARALLEL LAYERS ABOVE UPPER BOUNDARY
   REAL(SP), ALLOCATABLE :: ZKL(:)       !!THE DEPTHS OF PARALLEL LAYERS BELOW LOWER BOUNDARY

!--Bottom Friction Parameters--------------------------------------------------!
   REAL(SP) :: BFRIC        !!MINIMUM BOTTOM DRAG COEFFICIENT
   REAL(SP) :: Z0B          !!BOTTOM FRICTION DEPTH LENGTH SCALE

!--Flow Stability Parameters---------------------------------------------------!
   LOGICAL  :: VERT_STAB    !!TRUE IF VERTICAL ADJUSTMENT PERFORMED

!--Temp/Salinity Flux Enhancement Parameters-----------------------------------!
   LOGICAL  :: TS_FCT       !!TRUE IF TEMP/SALINITY ADJUSTMENT IS ON

!--Run Mode Parameters---------------------------------------------------------!
   LOGICAL  :: BAROTROPIC   !!TRUE CODE IS TO RUN IN BAROTROPIC MODE
   LOGICAL  :: SALINITY_ON  !!TRUE CODE IS TO UPDATE SALINITY EQUATION
   LOGICAL  :: TEMP_ON      !!TRUE CODE IS TO UPDATE TEMPERATURE EQUATION

!--Bathymetric Parameters------------------------------------------------------!
   REAL(SP) :: DJUST        !!WATER LEVEL CONSTANT ADJUSTMENT
   REAL(SP) :: MIN_DEPTH    !!MINIMUM ALLOWABLE DEPTH         
   REAL(SP) :: HMAX         !!GLOBAL MAXIMUM DEPTH
   REAL(SP) :: HMIN         !!GLOBAL MINIMUM DEPTH

!--Diffusion Coefficients/Viscosities------------------------------------------!
   REAL(SP) :: UMOL         !!VERTICAL DIFFUSION COEFFICIENT
   REAL(SP) :: HORCON       !!HORIZONTAL DIFFUSION COEFFICIENT
   REAL(SP) :: VPRNU        !!VERTICAL KINEMATIC VISCOSITY
   REAL(SP) :: HPRNU        !!HORIZONTAL KINEMATIC VISCOSITY

!--Water Heating Parameters----------------------------------------------------!
   REAL(SP) :: RHEAT        !!VERTICAL HEATING SWITCH 
   REAL(SP) :: ZETA1        !!VERTICAL HEATING LENGTH SCALE
   REAL(SP) :: ZETA2        !!SECONDARY VERTICAL HEATING LENGTH SCALE

!--Parameters Controlling Time/Time Stepping-----------------------------------!
   REAL(SP) :: DTE        !!EXTERNAL TIME STEP (Seconds)
   REAL(SP) :: DTI        !!INTERNAL TIME STEP (Seconds)
   REAL(SP) :: TIME       !!SIMULATION TIME IN DAYS
   REAL(SP) :: TIMERK     !!SIMULATION TIME AT CURRENT RK STAGE [UNUSED] 
   REAL(SP) :: THOUR      !!SIMULATION TIME AT CURRENT INTERNAL STEP (IINT) IN HOURS
   REAL(SP) :: THOUR1     !!SIMULATION TIME AT END OF CURRENT EXTERNAL STEP (IEXT) IN HOURS
   REAL(SP) :: RAMP       !!RAMP FACTOR USED TO EASE STARTUP = f(IINT)
   REAL(SP) :: DELTT      !!TIME BETWEEN JULIAN TIDE DATASETS 
   REAL(SP) :: THOUR_HS   !!TIME IN HOURS AT WHICH FLUX/BODY HEATING INITIATED
   REAL(SP) :: DAYS       !!TOTAL NUMBER OF DAYS IN SIMULATION
   INTEGER IINT           !!INTERNAL TIME STEP ITERATION NUMBER (ISTART => IEND)
   INTEGER IEXT           !!EXTERNAL TIME STEP ITERATION NUMBER (1 => ISPLIT)
   INTEGER ISTART         !!STARTING INTERNAL TIME STEP ITERATION NUMBER
   INTEGER IEND           !!ENDING INTERNAL TIME STEP ITERATION NUMBER
   INTEGER ISPLIT         !!NUMBER OF ITERATIONS OF EXTERNAL MODE/INTERNAL STEP
   INTEGER IRAMP          !!NUMBER OF INTERNAL ITERATIONS OVER WHICH TO RAMP UP
   INTEGER NSTEPS         !!ITERATION NUMBER AT END OF SIMULATION = IEND
!C. Lee, 12.21.2006,modified for restart T.W.
   REAL(SP) ::DTI_RATIO   !!INTERAL TIME STEP RATIO (DTIold/DTIrestart)
                          !!IF TIME STEP IS CHANGED WHEN RESTART, YOU NEED TO MAKE SURE MODEL
						  !STARTS AT SAME PHYSICAL TIME. HOWEVER, IT CAN BE USED TO RESTART AT ANY
						  !TIME BY SETTING THE PARAMETER VALUES.
!C.  Lee end T.W.
!--Parameters Controlling call rho_mean-----------------------------------!
   INTEGER IRHO_MEAN      !!NUMBER OF INTERNAL ITERATIONS OVER WHICH TO call rho_mean

!--Parameters Controlling Output
   INTEGER IRECORD        !!INTERVAL (IINT) FOR DUMPING BINARY OUTPUT
   INTEGER IDMPSMS        !!INTERVAL (IINT) FOR DUMPING SMS OUTPUT
   INTEGER IREPORT        !!INTERVAL (IINT) FOR REPORTING OF FLOWFIELD STATISTICS
   INTEGER IRESTART       !!INTERVAL (IINT) FOR RESTART FILE DUMPS

!--Parameters Related to Archiving of Flow Field Averages----------------------!
   LOGICAL AVGE_ON
   INTEGER BEG_AVGE
   INTEGER INT_AVGE
   INTEGER NUM_AVGE

!--Parameters Pertaining to Archiving------------------------------------------!
   REAL(SP) :: T1EL,T2EL
   
!--Parameter Related to Surface Wave Breaking and Ri # dissipation correction--!
   LOGICAL  :: SURFACEWAVE_MIX    !!TRUE IF SURFACE WIND INDUCED

!--Parameter Related to Open Boundary Condition Nudging------------------------!
   LOGICAL  :: TS_NUDGING_OBC     !!TRUE IF OBC Nudging Activated 
   REAL(SP) :: ALPHA_OBC          !!OBC NUDGING COEFFICIENT         

!--Parameter Related to Tidal Open Boundary Output-----------------------------!
#  if defined (TIDE_OUTPUT)
   INTEGER  :: TIDE_INITIAL       !!STARTING TIME OF OUTPUT (in time steps)
   INTEGER  :: TIDE_INTERVAL      !!OUTPUT INTERVAL (in time steps)  
#  endif

!--Parameter Related to Inclusion of Sediment Model----------------------------!
   LOGICAL  :: SEDIMENT_ON        !!TRUE IF SEDIMENT MODEL ACTIVE
   LOGICAL  :: RESTART_SED        !!TRUE IF SEDIMENT SHOULD BE RESTARTED

!--Option for semi-implicit coriolis force-------------------------------------!
   LOGICAL  :: ADCOR_ON

END MODULE CONTROL

!==============================================================================|

MODULE ALL_VARS
   USE MOD_PREC
   USE LIMS
   USE CONTROL 
   IMPLICIT NONE
   SAVE

!--------------------------Temporary Array------------------------------------------!

  INTEGER, ALLOCATABLE :: NVG(:,:)

!--------------------------Global Grid Variables------------------------------------!

  REAL(SP), ALLOCATABLE :: XG(:)               !!GLOBAL X-COORD AT NODE 
  REAL(SP), ALLOCATABLE :: YG(:)               !!GLOBAL X-COORD AT NODE 
  REAL(SP), ALLOCATABLE :: HG(:)               !!GLOBAL DEPTH AT NODE 
  REAL(SP), ALLOCATABLE :: XCG(:)              !!GLOBAL X-COORD AT FACE CENTER 
  REAL(SP), ALLOCATABLE :: YCG(:)              !!GLOBAL X-COORD AT FACE CENTER 

!--------------------------Grid Metrics---------------------------------------------!

   REAL(SP)              :: VXMIN,VYMIN,VXMAX,VYMAX
   REAL(SP), ALLOCATABLE :: XC(:)               !!X-COORD AT FACE CENTER 
   REAL(SP), ALLOCATABLE :: YC(:)               !!Y-COORD AT FACE CENTER
   REAL(SP), ALLOCATABLE :: VX(:)               !!X-COORD AT GRID POINT
   REAL(SP), ALLOCATABLE :: VY(:)               !!Y-COORD AT GRID POINT
   REAL(SP), ALLOCATABLE :: ART(:)              !!AREA OF ELEMENT
   REAL(SP), ALLOCATABLE :: ART1(:)             !!AREA OF NODE-BASE CONTROl VOLUME
   REAL(SP), ALLOCATABLE :: ART2(:)             !!AREA OF ELEMENTS AROUND NODE
   
!----------------Node, Boundary Condition, and Control Volume-----------------------!

   INTEGER, ALLOCATABLE :: NV(:,:)             !!NODE NUMBERING FOR ELEMENTS
   INTEGER, ALLOCATABLE :: NBE(:,:)            !!INDICES OF ELMNT NEIGHBORS
   INTEGER, ALLOCATABLE :: NTVE(:)      
   INTEGER, ALLOCATABLE :: NTSN(:)      
   INTEGER, ALLOCATABLE :: ISONB(:)            !!NODE MARKER = 0,1,2
   INTEGER, ALLOCATABLE :: ISBC(:)     
   INTEGER, ALLOCATABLE :: ISBCE(:)     
   INTEGER, ALLOCATABLE :: IEC(:,:)
   INTEGER, ALLOCATABLE :: IENODE(:,:)
   INTEGER, ALLOCATABLE :: NBSN(:,:)
   INTEGER, ALLOCATABLE :: NIEC(:,:)
   INTEGER, ALLOCATABLE :: NTRG(:)
   INTEGER, ALLOCATABLE :: NBVE(:,:)
   INTEGER, ALLOCATABLE :: NBVT(:,:)
   INTEGER, ALLOCATABLE :: LISBCE_1(:)          !!LIST OF ELEMENTS WITH ISBCE=1
   INTEGER, ALLOCATABLE :: LISBCE_2(:)          !!LIST OF ELEMENTS WITH ISBCE=2
   INTEGER, ALLOCATABLE :: LISBCE_3(:)          !!LIST OF ELEMENTS WITH ISBCE=3
   REAL(SP),ALLOCATABLE :: DLTXC(:)
   REAL(SP),ALLOCATABLE :: DLTYC(:)
   REAL(SP),ALLOCATABLE :: DLTXYC(:)
   REAL(SP),ALLOCATABLE :: DLTXE(:)
   REAL(SP),ALLOCATABLE :: DLTYE(:)
   REAL(SP),ALLOCATABLE :: DLTXYE(:)
   REAL(SP),ALLOCATABLE :: SITAC(:) 
   REAL(SP),ALLOCATABLE :: SITAE(:) 
   REAL(SP),ALLOCATABLE :: XIJC(:) 
   REAL(SP),ALLOCATABLE :: YIJC(:)
   REAL(SP),ALLOCATABLE :: XIJE(:,:) 
   REAL(SP),ALLOCATABLE :: YIJE(:,:) 
   REAL(SP),ALLOCATABLE :: EPOR(:)            !!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
   INTEGER, ALLOCATABLE :: IBCGEO(:)        !!LOCAL GEOSTROPHIC FRICTION CORRECTION NODES
   INTEGER, ALLOCATABLE :: N_ICELLQ(:,:)    !!FLUX ANGLE 

!----------------2-d arrays for the general vertical coordinate -------------------------------!

   REAL(SP), ALLOCATABLE :: Z(:,:)                    !!SIGMA COORDINATE VALUE 
   REAL(SP), ALLOCATABLE :: ZZ(:,:)                   !!INTRA LEVEL SIGMA VALUE
   REAL(SP), ALLOCATABLE :: DZ(:,:)                   !!DELTA-SIGMA VALUE
   REAL(SP), ALLOCATABLE :: DZZ(:,:)                  !!DELTA OF INTRA LEVEL SIGMA 
   REAL(SP), ALLOCATABLE :: Z1(:,:)                   !!SIGMA COORDINATE VALUE 
   REAL(SP), ALLOCATABLE :: ZZ1(:,:)                  !!INTRA LEVEL SIGMA VALUE
   REAL(SP), ALLOCATABLE :: DZ1(:,:)                  !!DELTA-SIGMA VALUE
   REAL(SP), ALLOCATABLE :: DZZ1(:,:)                 !!DELTA OF INTRA LEVEL SIGMA 
   REAL(SP), ALLOCATABLE :: DPTHSL(:)                 !!Z-DEPTHS FOR SALINITY/TEMP ICs
   REAL(SP), ALLOCATABLE :: S_HYBRID(:,:)		      !!Hyrbid sigma coordinates for output B clark 2015

!---------------2-d flow variable arrays at elements-------------------------------!

   REAL(SP), ALLOCATABLE :: UA(:)            !!VERTICALLY AVERAGED X-VELOC
   REAL(SP), ALLOCATABLE :: VA(:)            !!VERTICALLY AVERAGED Y-VELOC
   REAL(SP), ALLOCATABLE :: UAF(:)           !!UA FROM PREVIOUS RK STAGE 
   REAL(SP), ALLOCATABLE :: VAF(:)           !!VA FROM PREVIOUS RK STAGE 
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: UARK(:)          !!UA FROM PREVIOUS TIMESTEP 
   REAL(SP), ALLOCATABLE :: VARK(:)          !!VA FROM PREVIOUS TIMESTEP 
   REAL(SP), ALLOCATABLE :: UARD(:)          !!UA AVERAGED OVER EXTERNAL INT
   REAL(SP), ALLOCATABLE :: VARD(:)          !!VA AVERAGED OVER EXTERNAL INT
# endif
   REAL(SP), ALLOCATABLE :: COR(:)           !!CORIOLIS PARAMETER
   REAL(SP), ALLOCATABLE :: H1(:)            !!BATHYMETRIC DEPTH   
   REAL(SP), ALLOCATABLE :: D1(:)            !!CURRENT DEPTH
   REAL(SP), ALLOCATABLE :: DT1(:)           !!DEPTH AT PREVIOUS TIME STEP
   REAL(SP), ALLOCATABLE :: EL1(:)           !!CURRENT SURFACE ELEVATION
   REAL(SP), ALLOCATABLE :: ET1(:)           !!SURFACE ELEVATION AT PREVIOUS TIME STEP
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: ELRK1(:)         !!SURFACE ELEVATION AT BEGINNING OF RK INT 
# endif
   REAL(SP), ALLOCATABLE :: ELF1(:)          !!SURFACE ELEVATION STORAGE FOR RK INT
   REAL(SP), ALLOCATABLE :: DTFA(:)          !!ADJUSTED DEPTH FOR MASS CONSERVATION


   REAL(SP), ALLOCATABLE :: CC_SPONGE(:)     !!SPONGE DAMPING COEFFICIENT FOR MOMENTUM
   
!---------------2-d flow variable arrays at nodes----------------------------------!

   REAL(SP), ALLOCATABLE :: H(:)             !!BATHYMETRIC DEPTH   
   REAL(SP), ALLOCATABLE :: D(:)             !!CURRENT DEPTH   
   REAL(SP), ALLOCATABLE :: DT(:)            !!DEPTH AT PREVIOUS TIME STEP
   REAL(SP), ALLOCATABLE :: EL(:)            !!CURRENT SURFACE ELEVATION
   REAL(SP), ALLOCATABLE :: ET(:)            !!SURFACE ELEVATION AT PREVIOUS TIME STEP
   REAL(SP), ALLOCATABLE :: EGF(:)           !!AVERAGE SURFACE ELEVATION OVER EXTERNAL INT
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: ELRK(:)          !!SURFACE ELEVATION AT BEGINNING OF RK INT
# endif
   REAL(SP), ALLOCATABLE :: ELF(:)           !!SURFACE ELEVATION STORAGE FOR RK INT

!---------------surface/bottom boundary conditions---------------------------------!

   REAL(SP), ALLOCATABLE :: CBC(:)           !!BOTTOM FRICTION     
   REAL(SP), ALLOCATABLE :: SWRAD(:)         !!SURFACE INCIDENT RADIATION
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: WUSURF2(:)       !!SURFACE FRICTION FOR EXT
   REAL(SP), ALLOCATABLE :: WVSURF2(:)       !!SURFACE FRICTION FOR EXT
# endif
   REAL(SP), ALLOCATABLE :: WUBOT(:)         !!BOTTOM FRICTION
   REAL(SP), ALLOCATABLE :: WVBOT(:)         !!BOTTOM FRICTION
   REAL(SP), ALLOCATABLE :: WUSURF(:)        !!SURFACE FRICTION FOR INT
   REAL(SP), ALLOCATABLE :: WVSURF(:)        !!SURFACE FRICTION FOR INT
   REAL(SP), ALLOCATABLE :: WTSURF(:)        
   REAL(SP), ALLOCATABLE :: BFWDIS(:)        !!GROUNDWATER FLUX AT CURRENT TIME ! pengfei marked no need for ground water
   REAL(SP), ALLOCATABLE :: BFWTDIS(:)        !!GROUNDWATER TMEPERATURE  AT CURRENT TIME
   REAL(SP), ALLOCATABLE :: BFWSDIS(:)        !!GROUNDWATER SALINITY AT CURRENT TIME
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: BFWDIS2(:)        !!GROUNDWATER FLUX AT CURRENT TIME FOR EXT
   REAL(SP), ALLOCATABLE :: BFWTDIS2(:)       !!GROUNDWATER TMEPERATURE  AT CURRENT TIME FOR EXT
   REAL(SP), ALLOCATABLE :: BFWSDIS2(:)       !!GROUNDWATER SALINITY AT CURRENT TIME FOR EXT
# endif
   REAL(SP), ALLOCATABLE :: BFWDIS3(:)        !!GROUNDWATER FLUX AT CURRENT TIME FOR INT
   REAL(SP), ALLOCATABLE :: BFWTDIS3(:)       !!GROUNDWATER TMEPERATURE  AT CURRENT TIME FOR INT
   REAL(SP), ALLOCATABLE :: BFWSDIS3(:)       !!GROUNDWATER SALINITY AT CURRENT TIME FOR INT

   REAL(SP), ALLOCATABLE :: QDIS(:)          !!RIVER FLUX AT CURRENT TIME
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: QDIS2(:)         !!RIVER FLUX (EXT MODE, NOT USED) 
# endif
   REAL(SP), ALLOCATABLE :: TDIS(:)          !!RIVER WATER TEMP AT CURRENT TIME
   REAL(SP), ALLOCATABLE :: SDIS(:)          !!RIVER WATER SLNT AT CURRENT TIME
   REAL(SP), ALLOCATABLE :: QAREA(:)         !!AREA OF RIVER DISCHARGE 
   REAL(SP), ALLOCATABLE :: RDISQ(:,:)       !!AREA OF FLUX                    
   REAL(SP), ALLOCATABLE :: ANGLEQ(:)        !!RIVER DISCHARGE ANGLE           
   REAL(SP), ALLOCATABLE :: VLCTYQ(:)        !!RIVER DISCHARGE VELOCITY
   REAL(SP), ALLOCATABLE :: UUWIND(:)        !!SURFACE X-WIND 
   REAL(SP), ALLOCATABLE :: VVWIND(:)        !!SURFACE Y-WIND 
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: QPREC2(:)        !!SURFACE PRECIPITATION FOR EXT
   REAL(SP), ALLOCATABLE :: QEVAP2(:)        !!SURFACE EVAPORATION FOR EXT
# endif
   REAL(SP), ALLOCATABLE :: QPREC3(:)        !!SURFACE PRECIPITATION FOR INT
   REAL(SP), ALLOCATABLE :: QEVAP3(:)        !!SURFACE EVAPORATION FOR INT
!----------------boundary conditions: meteo conditions-----------------!
!# if defined (ICE) 
   REAL(SP), ALLOCATABLE  :: T_AIR(:)         !!BULK AIR TEMPERATURE( deg C ) AT HEIGHT ?(m)
   REAL(SP), ALLOCATABLE  :: RH_AIR(:)        !!RELATIVE HUMIDITY
   REAL(SP), ALLOCATABLE  :: QA_AIR(:)        !!SPECIFIC HUMIDITY  (KG/KG)
   REAL(SP), ALLOCATABLE  :: PA_AIR(:)        !!SURFACE PRESSURE
   REAL(SP), ALLOCATABLE  :: DLW_AIR(:)       !!DOWNWARD LONGWAVE RADIATION
   REAL(SP), ALLOCATABLE  :: DSW_AIR(:)       !!DOWNWARD SHORTWAVE RADIATION
                                              !!Total Short Wave Solar Radiation
   REAL(SP), ALLOCATABLE  :: CLOUD(:)         !!Total Cloud Cover
   REAL(SP), ALLOCATABLE  :: WSSURF(:)        !!  sea surface salinity    flux
   REAL(SP), ALLOCATABLE  :: WHSURF(:)        !!  sea surface heat        flux
   REAL(SP), ALLOCATABLE  :: FWtSURF(:)       !!  sea surface fresh water flux
!# endif


!!By T.W. for local array of kelp bed, April 2013

#  if defined (KELP)

   INTEGER,  ALLOCATABLE  :: n_of_kelp(:,:)         !!local array 
   REAL(SP), ALLOCATABLE  :: drag_coef(:,:)         !!local drag
   REAL(SP), ALLOCATABLE  :: kelp_area(:,:)         !!kelp area

   REAL(SP), ALLOCATABLE  :: EMS_X(:,:)             !! momentum sink east
   REAL(SP), ALLOCATABLE  :: EMS_Y(:,:)             !! momentum sink north

#  endif
!!finished by T.W.
   
! Added by B Clark for wind velocity block, Nov 2015
   
   INTEGER, allocatable :: WBTEMP(:)!,BLOCKED_WIND(:) ! Wind Block read in variable and 
   INTEGER :: Nblocked  ! Number of blocked nodes, read from input 

!-----------------------2-d flow fluxes--------------------------------------------!
   REAL(SP), ALLOCATABLE :: PSTX(:)           !!EXT MODE BAROTROPIC TERMS
   REAL(SP), ALLOCATABLE :: PSTY(:)           !!EXT MODE BAROTROPIC TERMS
   REAL(SP), ALLOCATABLE :: ADVUA(:) 
   REAL(SP), ALLOCATABLE :: ADVVA(:) 
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: ADX2D(:) 
   REAL(SP), ALLOCATABLE :: ADY2D(:) 
   REAL(SP), ALLOCATABLE :: DRX2D(:) 
   REAL(SP), ALLOCATABLE :: DRY2D(:) 
# endif
   REAL(SP), ALLOCATABLE :: TPS(:)            !!WORKING ARRAY             
# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE :: ADVX(:,:)      
   REAL(SP), ALLOCATABLE :: ADVY(:,:)     
# endif
!---------------- internal mode   arrays-(element based)----------------------------!
   
   REAL(SP), ALLOCATABLE :: U(:,:)         !X-VELOCITY
   REAL(SP), ALLOCATABLE :: V(:,:)         !Y-VELOCITY
   REAL(SP), ALLOCATABLE :: UBETA(:,:)     !X-VELOCITY temp time step
   REAL(SP), ALLOCATABLE :: VBETA(:,:)     !Y-VELOCITY temp time step
   REAL(SP), ALLOCATABLE :: UBETA2D(:)
   REAL(SP), ALLOCATABLE :: VBETA2D(:)

   REAL(SP), ALLOCATABLE :: W(:,:)         !VERTICAL VELOCITY IN SIGMA SYSTEM
   REAL(SP), ALLOCATABLE :: WW(:,:)        !Z-VELOCITY
   REAL(SP), ALLOCATABLE :: UF(:,:)        !X-VELOCITY FROM PREVIOUS TIMESTEP
   REAL(SP), ALLOCATABLE :: VF(:,:)        !Y-VELOCITY FROM PREVIOUS TIMESTEP
   REAL(SP), ALLOCATABLE :: WT(:,:)        !Z-VELOCITY FROM PREVIOUS TIMESTEP
   REAL(SP), ALLOCATABLE :: RHO(:,:)       !DENSITY AT ELEMENTS
   REAL(SP), ALLOCATABLE :: RMEAN(:,:)     !INITIAL DENSITY AT ELEMENTS
   REAL(SP), ALLOCATABLE :: T(:,:)         !TEMPERATURE AT ELEMENTS
   REAL(SP), ALLOCATABLE :: TMEAN(:,:)     !INITIAL TEMPERATURE AT ELEMENTS
   REAL(SP), ALLOCATABLE :: S(:,:)         !SALINITY AT ELEMENTS
   REAL(SP), ALLOCATABLE :: SMEAN(:,:)     !INITIAL SALINITY AT ELEMENTS
   REAL(SP), ALLOCATABLE :: Q2(:,:)        !2 X TURBULENT KINETIC ENERGY AT NODES
   REAL(SP), ALLOCATABLE :: L(:,:)         !TURBULENT LENGTH MACROSCALE 
   REAL(SP), ALLOCATABLE :: Q2L(:,:)       !2 X TURBULENT KE X LENGTH AT NODES
#  if defined (GOTM)
  REAL(SP), ALLOCATABLE :: TKE(:,:)       !TURBULENT KINETIC ENERGY AT NODES
  REAL(SP), ALLOCATABLE :: TKEF(:,:)      !TURBULENT KINETIC ENERGY AT NODES
  REAL(SP), ALLOCATABLE :: TEPS(:,:)      !TURBULENT DISSIPATION AT NODES
  REAL(SP), ALLOCATABLE :: TEPSF(:,:)     !TURBULENT DISSIPATION AT NODES
#  endif
  REAL(SP), ALLOCATABLE :: KM(:,:)        !TURBULENT EDDY VISCOSITY FOR MOMENTUM 
  REAL(SP), ALLOCATABLE :: KH(:,:)        !TURBULENT DIFFUSIVITY FOR SALINITY/TEMP 
  REAL(SP), ALLOCATABLE :: KQ(:,:)        !TURBULENT DIFFUSIVITY FOR Q2/Q2L 
  REAL(SP), ALLOCATABLE :: AAM(:,:)       !STORAGE FOR OUTPUT OF HORIZONTAL VISCOSITY
  REAL(SP), ALLOCATABLE :: Q2F(:,:)       !WORKING ARRAY FOR UPDATING Q2
  REAL(SP), ALLOCATABLE :: Q2LF(:,:)      !WORKING ARRAY FOR UPDATING Q2F 
  REAL(SP), ALLOCATABLE :: KM1(:,:)       !TURBULENT EDDY VISCOSITY FOR MOMENTUM 

  !-----------------------3d variable arrays-(node based)-----------------------------!

  REAL(SP), ALLOCATABLE :: T1(:,:)         !!TEMPERATURE AT NODES               
  REAL(SP), ALLOCATABLE :: S1(:,:)         !!SALINITY AT NODES               
  REAL(SP), ALLOCATABLE :: RHO1(:,:)       !!DENSITY AT NODES               
  REAL(SP), ALLOCATABLE :: TF1(:,:)        !!TEMPERATURE FROM PREVIOUS TIME
  REAL(SP), ALLOCATABLE :: SF1(:,:)        !!SALINITY FROM PREVIOUS TIME 
  REAL(SP), ALLOCATABLE :: TMEAN1(:,:)     !!MEAN INITIAL TEMP
  REAL(SP), ALLOCATABLE :: SMEAN1(:,:)     !!MEAN INITIAL SALINITY 
  REAL(SP), ALLOCATABLE :: RMEAN1(:,:)     !!MEAN INITIAL DENSITY 
  REAL(SP), ALLOCATABLE :: WTS(:,:)        !!VERTICAL VELOCITY IN SIGMA SYSTEM
  REAL(SP), ALLOCATABLE :: WTTS(:,:)       !!WTS FROM PREVIOUS TIMESTEP        

  !---------------------------internal mode fluxes------------------------------------!

  REAL(SP), ALLOCATABLE :: DRHOX(:,:)      !!BAROCLINIC PG IN X DIRECTION
  REAL(SP), ALLOCATABLE :: DRHOY(:,:)      !!BAROCLINIC PG IN Y DIRECTION

  !------------shape coefficient arrays and control volume metrics--------------------!

  REAL(SP), ALLOCATABLE :: A1U(:,:)      
  REAL(SP), ALLOCATABLE :: A2U(:,:)     
  REAL(SP), ALLOCATABLE :: AWX(:,:)   
  REAL(SP), ALLOCATABLE :: AWY(:,:)  
  REAL(SP), ALLOCATABLE :: AW0(:,:) 
  REAL(SP), ALLOCATABLE :: ALPHA(:)


  !-----salinity and temperature bottom diffusion condition/bottom depth gradients----!

  REAL(SP), ALLOCATABLE :: PHPN(:)   
  REAL(SP), ALLOCATABLE :: PFPXB(:)   
  REAL(SP), ALLOCATABLE :: PFPYB(:)   
  REAL(SP), ALLOCATABLE :: SITA_GD(:)  
  REAL(SP), ALLOCATABLE :: AH_BOTTOM(:)  

  !-----arrays used for averaging of flow quantities for output-----------------------!

  REAL(SP), ALLOCATABLE :: U_AVE(:,:)   !U AVERAGED OVER INT_AVGE ITERATIONS
  REAL(SP), ALLOCATABLE :: V_AVE(:,:)   !V AVERAGED OVER INT_AVGE ITERATIONS
  REAL(SP), ALLOCATABLE :: W_AVE(:,:)   !WW AVERAGED OVER INT_AVGE ITERATIONS
  REAL(SP), ALLOCATABLE :: KM_AVE(:,:)  !KM AVERAGED OVER INT_AVGE ITERATIONS
  REAL(SP), ALLOCATABLE :: KH_AVE(:,:)  !KH AVERAGED OVER INT_AVGE ITERATIONS
  REAL(SP), ALLOCATABLE :: T_AVE(:,:)   !T1 AVERAGED OVER INT_AVGE ITERATIONS
  REAL(SP), ALLOCATABLE :: S_AVE(:,:)   !S1 AVERAGED OVER INT_AVGE ITERATIONS
  REAL(SP), ALLOCATABLE :: R_AVE(:,:)   !RHO1 AVERAGED OVER INT_AVGE ITERATIONS
  REAL(SP), ALLOCATABLE :: EL_AVE(:)    !EL AVERAGED OVER INT_AVGE ITERATIONS

  REAL(SP), ALLOCATABLE :: VISCOFH(:,:)


  REAL(DP), ALLOCATABLE :: XCA(:)       !Longitude of center point of arc
  REAL(DP), ALLOCATABLE :: YCA(:)       !Latitude of center point of arc
  REAL(DP), ALLOCATABLE :: XCB(:)       !Longitude of center point of arc
  REAL(DP), ALLOCATABLE :: YCB(:)       !Latitude of center point of arc

  REAL(DP), ALLOCATABLE :: XCC(:,:)       !Longitude of center point of arc
  REAL(DP), ALLOCATABLE :: YCC(:,:)       !Latitude of center point of arc
  REAL(DP), ALLOCATABLE :: XCD(:,:)       !Longitude of center point of arc
  REAL(DP), ALLOCATABLE :: YCD(:,:)       !Latitude of center point of arc

  REAL(DP), ALLOCATABLE :: XCE(:)       !Longitude of center point of arc
  REAL(DP), ALLOCATABLE :: YCE(:)       !Latitude of center point of arc
  REAL(DP), ALLOCATABLE :: XCF(:)       !Longitude of center point of arc
  REAL(DP), ALLOCATABLE :: YCF(:)       !Latitude of center point of arc

  REAL(DP), ALLOCATABLE :: XCG2(:)       !Longitude of center point of arc
  REAL(DP), ALLOCATABLE :: YCG2(:)       !Latitude of center point of arc
  
  REAL(SP), ALLOCATABLE :: VAL_COS_VY(:) !THE COSINE VALUE OF LATITUDE
  
END MODULE ALL_VARS

MODULE BCS
   USE MOD_TYPES
   USE MOD_PREC
   IMPLICIT NONE
   SAVE


!----------------boundary conditions: ground water----------------------------------!

   INTEGER, ALLOCATABLE :: NODE_BFW(:)      !!LOCAL GROUNDWATER NODES
   INTEGER, ALLOCATABLE :: BFW_GL2LOC(:)    !!GLOBAL TO LOCAL MAPPING OF GWATER NODES
   REAL(SP),  ALLOCATABLE :: BFWQDIS(:,:)   !!GROUNDWATER FRESH WATER FLUX DATA
   REAL(SP),  ALLOCATABLE :: BFWQTDIS(:,:)   !!GROUNDWATER WATER  temperature  DATA
   REAL(SP),  ALLOCATABLE :: BFWQSDIS(:,:)   !!GROUNDWATER WATER salinity DATA

   TYPE(BC)      :: BFW_TM                  !!TIME MAP FOR GROUNDWATER DATA

!----------------boundary conditions: non-Julian tidal forcing----------------------!
 
   REAL(SP),  ALLOCATABLE :: APT(:,:)      !!TIDE AMPLITUDE
   REAL(SP),  ALLOCATABLE :: PHAI(:,:)     !!TIDE PHASE 
   REAL(SP),  ALLOCATABLE :: EMEAN(:)      !!MEAN SURFACE ELEVATION

!----------------boundary conditions: Julian tidal forcing--------------------------!

   REAL(SP), ALLOCATABLE    :: ELSBC(:,:)    !!INPUT SURFACE ELEVATION
   TYPE(BC)                 :: ELO_TM        !!TIME MAP FOR SURFACE ELEVATION DATA

!----------------boundary conditions: fresh water discharge-------------------------!

   INTEGER, ALLOCATABLE  :: INODEQ(:)        !!LOCAL FRESH WATER INFLOW NODES
   INTEGER, ALLOCATABLE  :: ICELLQ(:)        !!LOCAL FRESH WATER INFLOW ELEMENTS 
   INTEGER, ALLOCATABLE  :: RIV_GL2LOC(:)    !!GLOBAL TO LOCAL MAP OF FW NODES
   TYPE(BC)              :: QBC_TM           !!TIME MAP FOR RIVER DISCHARGE DATA 
   REAL(SP), ALLOCATABLE  :: VQDIST(:,:)     !!DISCHARGE VERTICAL DISTRIBUTION
   REAL(SP), ALLOCATABLE  :: DQDIS(:,:)      !!WATER FLUX DISCHARGE DATA
   REAL(SP), ALLOCATABLE  :: DTDIS(:,:)      !!WATER TEMPERATURE DISCHARGE DATA
   REAL(SP), ALLOCATABLE  :: DSDIS(:,:)      !!WATER SALINITY DISCHARGE DATA

!----------------surface forcing: heat flux and radiation---------------------------!

   TYPE(BC)               :: HFX_TM           !!TIME MAP FOR TIME DEPENDENT HF/RAD DATA
   REAL(SP), ALLOCATABLE  :: DHFLUX(:,:)      !!HEAT FLUX DATA
   REAL(SP), ALLOCATABLE  :: DHSHORT(:,:)     !!SHORT WAVE RADIATION DATA
   
!----------------surface forcing: wind field----------------------------------------!

   TYPE(BC)               :: WND_TM           !!TIME MAPPING FOR WIND FIELD DATA
   REAL(SP), ALLOCATABLE  :: DTX(:,:)         !!SURFACE WIND (X)
   REAL(SP), ALLOCATABLE  :: DTY(:,:)         !!SURFACE WIND (Y)

!----------------surface forcing: evaporation and precipitation---------------------!

   TYPE(BC)               :: EVP_TM           !!TIME MAPPING FOR EVAP AND PREC DATA
   REAL(SP), ALLOCATABLE  :: DQPREC(:,:)      !!PRECIPITATION
   REAL(SP), ALLOCATABLE  :: DQEVAP(:,:)      !!EVAPOTATION

!----------------boundary conditions: uniform meteo conditions----------------------!

   TYPE(BC)               :: UMF_TM           !!TIME MAPPING FOR UNIFORM METEOS
   REAL(SP), ALLOCATABLE  :: UQPREC(:)        !!PRECIPITATION
   REAL(SP), ALLOCATABLE  :: UQEVAP(:)        !!EVAPORATION   
   REAL(SP), ALLOCATABLE  :: UWIND(:)         !!SURFACE WIND (X)
   REAL(SP), ALLOCATABLE  :: VWIND(:)         !!SURFACE WIND (Y)
   REAL(SP), ALLOCATABLE  :: UHFLUX(:)        !!UNIFORM HEAT FLUX 
   REAL(SP), ALLOCATABLE  :: UHSHORT(:)       !!UNIFORM SHORT WAVE RADIATION

!==============================================================================|

END MODULE  BCS
