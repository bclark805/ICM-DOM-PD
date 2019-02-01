!==============================================================================!
!  VERSION 2.7.1 
!==============================================================================!

   PROGRAM US_FVCOM

!==============================================================================!
!                                                                              !
!                            USG-FVCOM                                         !
!    The Unstructured Grid Finite Volume Coastal Ocean Model                   !
!                                                                              !
!    The USG-FVCOM (publically called FVCOM) was developed by Drs. Changsheng  !
!  Chen and Hedong Liu at the Marine Ecosystem Dynamics Modeling Laboratory    !
!  at the School of Marine Science and Technology (SMAST), University of       !
!  Massachusetts-Dartmouth (UMASSD) and Dr. Robert C. Beardsley at the         !
!  Department of Physical Oceanography, Woods Hole Oceanographic Institution   !
!  (WHOI). This code was rewritten in Fortran 90/2K, modularized, rendered     ! 
!  somewhat understandable, and parallelized by Geoff Cowles at SMAST/UMASSD.  !
!                                                                              !
!    The Development was initially supported by the Georgia Sea Grant College  !
!  Program for the study of the complex dynamic system in Georgia estuaries.   !
!  The code improvement has been supported by Dr. Chen's research grants       !
!  received from NSF and NOAA Coastal Ocean research grants and SMAST-NASA     !
!  fishery research grants. The University of Massachusetts-Dartmouth          !
!                                                                              !
!    FVCOM is a three dimensional,time dependent,primitive equations,          !
!  coastal ocean circulation model. The model computes the momentum,           !
!  continuity, temperature, salinity, and density equations and is closed      !
!  physically and mathematically using the Mellor and Yamada level-2.5         !
!  turbulent closure submodel. The irregular bottom slope is represented       !
!  using a sigma-coordinate transformation, and the horizontal grids           !
!  comprise unstructured triangular cells. The finite-volume method (FVM)      !
!  used in this model combines the advantages of the finite-element            !
!  method (FEM) for geometric flexibility and the finite-difference            !
!  method (FDM) for simple discrete computation. Current, temperature,         !
!  and salinity in the model are computed in the integral form of the          !
!  equations, which provides a better representation of the conservation       !
!  laws for mass, momentum, and heat in the coastal region with complex        !
!  geometry.                                                                   !
!                                                                              !
!    All users should read this agreement carefully.  A user, who receives any !  
!  version of the source code of FVCOM, must accept all the terms and          !
!  conditions of this agreement and also agree that this agreement is like any !
!  written negotiated agreement signed by you. You may be required to have     !
!  another written agreement directly with Dr. Changsheng Chen at SMAST/UMASS-D!
!  or Dr. Brian Rothschild, Director of the SMAST/UMASS-D that supplements     !
!  all or portions of this agreement. Dr. Changsheng Chen, leader of the       !
!  FVCOM development team, owns all intellectual property rights to the        !
!  software. The University of Massachusetts-Dartmouth and the Georgia Sea     !
!  Grant Program share the copyright of the software. All copyrights are       !
!  reserved. Unauthorized reproduction and re-distribution of this program     !
!  are expressly prohibited. This program is only permitted for use in         !
!  non-commercial academic research and education.  Commercial use must be     !
!  approved by Dr. Chen and is subject to a license fee. Registration is       !
!  required for all new users.  Users should realize that this model software  !
!  is a research product without any warranty. Users must take full            !
!  responsibility for any mistakes that might be caused by any inappropriate   !
!  modification of the source code.  Modification is not encouraged for users  !
!  who do not have a deep understanding of the finite-volume numerical methods !
!  used in FVCOM. Contributions made to correcting and modifying the programs  !
!  will be credited, but will not affect copyrights. No duplicate              !
!  configurations of FVCOM are allowed in the same geographical region,        !
!  particularly in the regions where FVCOM has been already been applied.      !
!  Users who want to use FVCOM in a region that the SMAST/UMASS Marine         !
!  Ecosystem Dynamics Modeling (MEDM) group (led by Dr. Chen) is working on    !
!  must request permission from Dr. Chen. No competition is allowed in the     !
!  same region using FVCOM, especially with Dr. Chen's group. FVCOM has been   !
!  validated for many standard model test cases.  Users are welcome to do any  !
!  further model validation experiments. These experiments shall be carried    !
!  out in collaboration with the SMAST/UMASSD model development team. To avoid !
!  or reduce deriving any incorrect conclusions due to an inappropriate use of !
!  FVCOM, users are required to contact the scientific leaders of the FVCOM    !
!  development team (Dr. Chen at SMAST/UMASS-D and Dr. Beardsley at WHOI)      !
!  before any formal publications are prepared for model validation.           !
!                                                                              !
!    For public use, all users should name this model as "FVCOM". In any       !
!  publications with the use of FVCOM, acknowledgement must be included. The   !
!  rationale behind this FVCOM distribution policy is straightforward.  New    !
!  researchers and educators who want to use FVCOM and agree to the above      !
!  requirements get free access to the latest version of FVCOM and the         !
!  collective wisdom and experience of the FVCOM development team and existing !
!  users. Problems arising in new FVCOM applications, both related to          !
!  conceptual as well as numerical and coding issues, can be shared with the   !
!  development team and other users who can work together on physics and code  !
!  improvements that over time will lead to a better FVCOM.                    !
!                                                                              !
!    FVCOM has been developed to date with state and federal funding with the  !
!  idea that FVCOM will become a community model that new users start to use   !
!  the model and its scientific usefulness and numerical accuracy and          !
!  efficiency continue to improve.  The FVCOM distribution policy is designed  !
!  to encourage this transition while maintaining a central core group         !
!  responsible for overall FVCOM development direction, implementing official!
!  code improvements, and maintaining well tested and documented updated code  !
!  versions.                                                                   !       
!                                                                              !
!                                                                              !
!  External forces used to drive this model:                                   !
!                                                                              !
!  1) Tidal amplitudes and phases at open boundaries (initial designs          !
!         include 6 tidal consituents, more can be added as needed);           !
!  2) Wind Stress [3 ways: a) uniform wind speed and direction, b)             !
!         spatially distributed wind velocity field, and c)the MM5 model-out   !
!         wind fields]                                                         !
!  3) Surface heat flux [3 ways: a) uniform heat flux, b) spatially            !
!         distributed heat flux, and c) the MM5-output heat flux fields        !
!         All the surface net heat flux and short-wave radiation are needed    !
!         in the input file                                                    ! 
!  4) River discharges: specify the location and discharge volume,             !
!         temperature, and salinity                                            !
!  5) Groundwater input: currently diffused bottom flux only                   !
!                                                                              !
!  Initial conditions:                                                         !
!                                                                              !
!  The model can be prognostically run for both barotropic and baroclinic      !
!  cases.                                                                      !
!                                                                              !
!  Tidal forcing can be added into the system with zero velocity               !
!  field at initial or specified the 3-D tidal initial velocity field          !
!  using the model-predicted harmonic tidal currents.                          !
!                                                                              !
!  Initial fields of temperature and salinity needed to be specified           !
!  by using either climatological field, real-time observed field or           !
!  idealized functions. The model has included Gregorian time for the          !
!  time simulation for tidal currents.                                         !
!                                                                              !
!  For the purpose of interdisciplinary studies, biological, chemical, and     !
!  sediment suspension models are available for FVCOM.  These submodels are    !
!  directly driven by the FVCOM physical model. A description of these         !
!  submodels follows.                                                          !
!                                                                              !
!  Generalized biological modules-a software platform that allows users to     !
!             build his own biological model in FVCOM                          !
!                                                                              !
!  NPZ model--a 3 component nutrient-phytoplankton-zooplankton model           !
!                                                                              !
!  NPZD model--an 8 component nutrient-phytolankton-zooplankton-detritus       !
!              model;                                                          !
!                                                                              !
!  NPZDB-model-a 9 phosphorus-controlled component nutrient-                   !
!               phytoplankton-zooplankton-detritus-bacteria model;             !
!                                                                              !
!  Water quality model with inclusion of the benthic process                   !
!                                                                              !
!  Sediment model--A new module that was developed by Dr. Geoff Cowles         !
!                                                                              !
!  Lagrangian particle tracking:                                               !
!                                                                              !
!  A bilinear interpolation scheme is used to determine the particle           !
!  velocity for the Lagrangian particle tracking. A random walk process        !
!  also could be included with a specified function related to horizontal      !
!  and vertical diffusion coefficients                                         !
!                                                                              !
!  Key reference:                                                              !
!                                                                              !
!   Chen, C., H. Liu, and R. C. Beardsley, 2003. An unstructured grid,         !
!       finite-volume, three-dimensional, primitive equations ocean            !
!       model: application to coastal ocean and estuaries, Journal             !
!       of Atmospheric and Oceanic Technology,  20, 159-186.                   !
!                                                                              !
!                                                                              !
!                                                                              !
!  Please direct criticisms and suggestions to                                 !
!                                                                              !
!               Changsheng Chen                                                !
!               School for Marine Science and Technology                       !
!               University of Massachusetts-Dartmouth                          !
!               New Bedford, MA 02742                                          !
!               Phone: 508-910-6388, Fax: 508-910-6371                         !
!               E-mail: c1chen@umassd.edu                                      !
!               Web: http://fvcom.smast.umassd.edu                             !
!                                                                              !
! What are new for version 2.7.1?                                           !
!    1) Multiple choice to set up the coordinate in the vertical is developed  !
!       by Dr. Qi                                                              !
!       a) sigma levels                                                        !
!       b) general vertical levels                                             ! 
!       c) constant layer transformation                                       !
!                                                                              !
! What are new for version 2.6?                                                !
!    1) A new Lagrangian particle tracking module is added with multiprocessor !
!       formulation, restart capability and general tracking et al.            !
!       (by Dr. Cowles)                                                        !
!                                                                              !
! What are new for version 2.5?                                                !
!    1) A new spherical coordinate version is added with an accurate treatment !
!       of the north pole (Arctic Ocedan) by Dr. Chen et al. at UMASSD         !
!    2) Spherical and x-y coordinate versions was merged into a single code    !
!       with users' choice for either coordiante for their application         !
!    3) Multiple choices of open radiation boundary conditions are added       !
!    4) General turbulence modules were linked by Dr. Cowles                   !
!    5) The selection for a 2-D barotrophic application is added               !
!    6) bugs in paralleziation and wet/dry point treatments are corrected      !
! For more detailed information, please read the upgrade user manaul           !
!                                                                              !
! What will be included in version 2.4 (under validation tests)                !
!    1) Generalized Biological Modules (developed by Dr. Tian and Chen)        !
!    2) 3-D sediment model(developed by Dr. Cowles)                            !
!    3) Reduced Kalman Filter, Ensemble Kalman Filter, and Ensemble Tranistion !
!       Kalman Filter (implented into FVCOM with UMASSD and MIT groups led     !
!       by Chen and Roziili                                                    !
!                                                                              !   
! Enjoy!                                                                       !
!==============================================================================!

!==============================================================================!
!  INCLUDE MODULES                                                             !
!==============================================================================!


   USE ALL_VARS
# if defined (MULTIPROCESSOR)
   USE MOD_PAR  
# endif

# if defined (VISIT)
  USE MOD_VISIT, only : VISIT_OPT, VISIT_CMD_DUMP, visit_time_ext,&
       & visit_time_int, visit_cycle
# endif

   USE MOD_CLOCK
#  if defined (NG_OI_ASSIM)
   USE MOD_ASSIM
#  endif
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  if defined (NORTHPOLE)
   USE MOD_NORTHPOLE
#  endif
#  endif

#  if defined (WET_DRY)
   USE MOD_WD
#  endif
#  if defined (NETCDF_IO)
   USE MOD_NCDIO
   USE MOD_NCDAVE
#  endif
   USE BCS
   USE PROBES     
   USE MOD_LAG
#  if defined (GOTM)
   USE MOD_GOTM
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
!  SEDIMENT BEGIN
#  if defined (SEDIMENT)
   USE MOD_SED
#  endif
!  SEDIMENT END
#  if defined (BioGen)
   USE MOD_BIO_3D
#  endif
#  if defined (ONE_D_MODEL)
   USE MOD_ONEDTIDE
#  endif   

   USE MOD_OBCS

! New Open Boundary Condition ----1
#  if defined (MEAN_FLOW)
   USE MOD_MEANFLOW
   USE MOD_OBCS2
   USE MOD_OBCS3
#  endif

#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif 
#  if defined (TS_OBC)
   USE MOD_TSOBC
#  endif   
#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
#  endif

#  if defined (ENKF_ASSIM)
   USE MOD_ENKF
#  endif
#  if defined (ETKF_ASSIM)
   USE MOD_ETKF
#  endif
#  if defined (RRK_PRE)
   USE MOD_RRK
#  endif
#  if defined (RRK_ASSIM)
   USE MOD_RRK
   USE MOD_RRKA 
#  endif 
#  if defined (HEAT_FLUX) || (ICE) 
   USE MOD_HEATFLUX
#  endif

#  if defined (ICE)
   USE MOD_ICE
#  endif

#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif
#  if defined (SEMI_IMPLICIT)
   USE MOD_PETSC
#  endif

! TW, to use Module velocity block
#  if defined (V_BLOCK)
   USE VELOCITY_BLOCK
#  endif

! added by T.W. for Kelp bed drag
# if defined (KELP)
  USE MOD_KELP
# endif

! finished by T.W. April 2013

 !added by Wen Long for outputing flux
# if defined (OUTPUT_FLUX)  
   USE MOD_OUTPUT_FLUX
# endif   
  
!------------------------------------------------------------------------------|
   IMPLICIT NONE

   LOGICAL  :: FEXIST
   REAL(SP) :: TMP1,TMP,UTMP,VTMP,TTIME
   INTEGER  :: I,K,J,IERR,N1,J1,J2,JN
   REAL(SP), ALLOCATABLE :: FTEMP(:),FTEMP2(:)
   CHARACTER(LEN=13) :: TSTRING
!   INTEGER, ALLOCATABLE :: kount(:)   !Added by TW for counting the total number of blocked vertical layers
!   LOGICAL  :: block_count            !added by TW 

# if defined (RRK_PRE)
   CHARACTER(LEN=120):: RRKFILES
   CHARACTER(LEN=8)  :: RRKINP1
   CHARACTER(LEN=4)  :: RRKINP2
# endif
   character(len=4)  :: ch4
   character(len=8)  :: ch8
   CHARACTER(LEN=120) :: FNAME

# if defined (ICE) 
    real(SP) :: DUVI,ice_ocnDX,ice_ocnDY
# endif

  integer temp
  integer :: irecn(4), nmj1(201)
  real(sp), allocatable :: el_tmp(:), u_tmp(:), v_tmp(:)

!------------------------------------------------------------------------------|

!------------------------------------------------------------------------------|
!   block_count = .TRUE.   !added by TW

!==============================================================================!
!  FVCOM VERSION                                                               !
!==============================================================================!

   FVCOM_VERSION     = 'FVCOM_2.7'
   FVCOM_WEBSITE     = 'http://fvcom.smast.umassd.edu'
#  if defined (NETCDF_IO)
   INSTITUTION       = 'School for Marine Science and Technology'
   NETCDF_TIMESTRING = 'seconds after 00:00:00'
#  endif

!==============================================================================!
!   SETUP PARALLEL ENVIRONMENT                                                 !
!==============================================================================!

   SERIAL = .TRUE. 
   PAR    = .FALSE. 
   MSR    = .TRUE.
   MYID   = 1
   NPROCS = 1
#  if defined (MULTIPROCESSOR)
   CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR)
#  endif

!==============================================================================!
!   IMPORT CASENAME FROM COMMAND LINE                                          !
!==============================================================================!

   CALL GET_CASENAME(CASENAME)

!==============================================================================!
!   SETUP MODEL RUN                                                            !
!==============================================================================!


!  READ PARAMETERS CONTROLLING T/S OBC Series NUDGING
!  
! B Clark moved first in sequence so can report in log file if time series nudging is active and the coefficient
#  if defined (TS_OBC) 
   CALL SET_TSOBC_PARAM  
#  endif   
!
!  READ PARAMETERS CONTROLLING MODEL RUN
!
!
   CALL DATA_RUN
#  if defined (HEAT_FLUX) || (ICE)
   CALL DATA_RUN_HFX
#  endif

! T.W., April 2013, for kelp
#  if defined (KELP)
   CALL DATA_RUN_KELP
#  endif

! W. Long, April 2014, for output transect flux  
#  if defined (OUTPUT_FLUX)
   CALL DATA_RUN_FLUX
#  endif
   
# if defined (VISIT)
   If(MSR) WRITE(IPT,*)"INITILIZING VISIT"
   Call INIT_VISIT
# endif

#  if defined (ONE_D_MODEL)
   CALL GET_UMAX
#  endif      
!
!  READ PARAMETERS CONTROLLING OBC FOR TEMPERATURE AND SALINITY
!
   CALL TSOBC_TYPE
!
!  READ PARAMETERS CONTROLLING DATA ASSIMILATION
!
#  if defined (NG_OI_ASSIM)
   CALL SET_ASSIM_PARAM
#  endif
!
!  READ PARAMETERS CONTROLLING WET_DRY TREATMENT
!
#  if defined (WET_DRY)
   CALL SET_WD_PARAM
#  endif

#  if defined (BALANCE_2D)
   CALL SET_BALANCE_PARAM
#  endif

!  READ PARAMETERS CONTROLLING ENSEMBLE KALMAN FILTERS ASSIMILATION
!
#  if defined (ENKF_ASSIM)
   CALL SET_ENKF_PARAM
#  endif
#  if defined (ETKF_ASSIM)
   CALL SET_ETKF_PARAM
#  endif
!
!  READ PARAMETERS CONTROLLING KALMAN FILTERS ASSIMILATION
!
#  if defined (RRK_PRE)
   CALL SET_RRK_PARAM
#  endif
#  if defined (RRK_ASSIM)
   CALL SET_RRK_PARAM
#  endif

!
!  READ PARAMETERS CONTROLLING DYE RELEASE
!
#  if defined (DYE_RELEASE)
   CALL SET_DYE_PARAM
#  endif   

!
!  OPEN INPUT/OUTPUT FILES
!
   CALL IOFILES
#  if defined (HEAT_FLUX) || (ICE)
   IF(C_HFX)CALL IOFILES_HFX  
!   write(*,*)'opened the hfx input files'
#  endif

!
!  DETERMINE NUMBER OF ELEMENTS AND NODES IN THE MODEL
!
   CALL GETDIM
!
!  READ PARAMETERS CONTROLLING WATER QUALITY MODEL 
!
#  if defined (WATER_QUALITY)
   CALL GET_WQMPAR
#  endif
!
!  READ PARAMETERS CONTROLLING NETCDF OUTPUT
!

#  if defined (NETCDF_IO)
   CALL SET_NCD_IO
   IF(AVGE_ON) CALL SET_NCD_AVE
#  endif

!
!  DECOMPOSE DOMAIN BY ELEMENTS USING METIS
!


#  if defined (MULTIPROCESSOR) 
   ALLOCATE(EL_PID(MElemGL))  ;  EL_PID = 1
   IF(PAR)CALL DOMDEC(INGRD,MElemGL,NPROCS,EL_PID,MSR)

#  endif 

!
!  GENERATE GLOBAL<==>LOCAL ELEMENT/NODE MAPPING
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL GENMAP
#  endif
 
!
!  MAP OPEN BOUNDARY CONDITION NODES TO LOCAL DOMAIN
!
   CALL BCMAP

!
!  INPUT AND SETUP BOUNDARY FORCING (HEAT/RIVERS/WIND/etc)
!
   CALL BCS_FORCE

#  if defined (HEAT_FLUX) || (ICE)
   IF(C_HFX)CALL BCS_FORCE_HFX
#  endif

!
!  INPUT WATER QUALITY MODEL VARIABLES
!
#  if defined (WATER_QUALITY)
   IF(WQM_ON)CALL BCS_FORCE_WQM
#  endif

#  if defined (TS_OBC)
!
!  INPUT T/S OBC Series NUDGING VARIABLES
!
   IF(TSOBC_ON) CALL READ_TSOBC
#  endif 

!
!  ALLOCATE FLOWFIELD VARIABLES
!
   CALL ALLOC_VARS

!
!  ALLOCATE ELM1 AND ELM2 FOR ORLANSKI RADIATION OPEN BOUNDARY CONDITION
!
   CALL ALLOC_OBC_DATA

!
!  ALLOCATE AND WET/DRY CONTROL ARRAYS
!
#  if defined (WET_DRY)
   CALL ALLOC_WD_DATA
#  endif
!
!  ALLOCATE WATER QUALITY MODEL VARIABLES
!
#  if defined (WATER_QUALITY)
   IF(WQM_ON)CALL ALLOC_WQM_VARS
#  endif

!
!  ALLOCATE SPHERICAL COORDINATE SYSTEM VARS
!
#  if defined (SPHERICAL)
   CALL ALLOC_SPHERE_VARS
#  endif

!
!  ALLOCATE MOMENTUM BALANCE CHEKING VARS
!
#  if defined (BALANCE_2D)
   CALL ALLOC_BALANCE_VARS
#  endif
!
!  ALLOCATE DYE VARS
!
#  if defined (DYE_RELEASE)
   CALL ALLOC_VARS_DYE
#  endif

#  if defined (BioGen)
   KBV=KB
   CALL GET_PARAMETER
   POINT_ST_TYPE='NONE' !No bio river yet
#  endif
!
!  SHIFT GRID/CORIOLIS/BATHYMETRY TO LOCAL DOMAIN
!
   CALL PDOMDEC
!
!  ALLOCATE EQUILIBRIUM AND ATMOSPHERIC TIDE VARS
!
# if defined (EQUI_TIDE)
  CALL ALLOCATE_EQUI
# endif  
# if defined (ATMO_TIDE)
  CALL ALLOCATE_ATMO
# endif  

!
!  SET UP GRID METRICS (FLUX EDGES/CONTROL VOLUMES/ETC)
!

!!add by TW for velocity BLOCK
# if defined (V_BLOCK)
  CALL ALLOC_VAR_BLOCK   !provide default values
  CALL SET_BLOCK_PARAM   !further update if BLOCK = T
# endif

   CALL TRIANGLE_GRID_EDGE      !Set up fluxes and control Volumes
   CALL SET_SIGMA(INDEX_VERCOR) !Build Vertical Coordinate
   CALL CELL_AREA               !Calculate Element and Control Volume Areas
# if defined (GCN)
   CALL SHAPE_COEF_GCN          !Calc Shape Coefficients for Flux Construction
# else
   CALL SHAPE_COEF_GCY          !Calc Shape Coefficients for Flux Construction
# endif

!T.W., for Kelp, April 2013
#  if defined (KELP)
   IF(C_KELP) CALL IOFILES_KELP
#  endif
   
!W. Long, for output of transect flux, April, 2014
   
# if defined (OUTPUT_FLUX)
   IF(WRITE_TCE_FLUX) CALL IOFILES_FLUX
# endif

#  if defined (SPHERICAL) && (NORTHPOLE)
   CALL SHAPE_COEF_XY
!  FIND THE NODE NUMBER OF NORTH POLE IF EXIST
   CALL FIND_NORTHPOLE
   CALL FIND_CELLSIDE
#  endif   

!  SET ISBCE AND ISONB CORRECTLY IN HALO CELLS/NODES
# if defined (MULTIPROCESSOR)
  ALLOCATE(FTEMP(0:MTElem)) ; FTEMP = ISBCE
  IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,FTEMP)
  ISBCE = FTEMP
  DEALLOCATE(FTEMP)
  ALLOCATE(FTEMP(0:NTNode)) ; FTEMP = ISONB
  IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,FTEMP)
  ISONB = FTEMP
  DEALLOCATE(FTEMP)
# endif


!below is commented out and moved to front, for the purpose of sigma coordiate modification proposed by Tarang
!!add by TW for velocity BLOCK
!# if defined (V_BLOCK)
!  CALL ALLOC_VAR_BLOCK   !provide default values
!  CALL SET_BLOCK_PARAM   !further update if BLOCK = T
!# endif

!
!  SETUP OPEN BOUNDARY METRICS                                   
!
   CALL SETUP_OBC 

   CALL SET_BNDRY               !Boundary Condition Metrics

#  if defined (SEMI_IMPLICIT)
   CALL SET_IMPLICIT_PARAM
   CALL ALLOC_VARS_SEMI
#  endif
#  if defined (SEMI_IMPLICIT)
   CALL PETSc_SET
#  endif

!
!  INITIALIZE FLOWFIELD  --> [T,S,U,V,EL,D,Q2,Q2L]
!
   CALL STARTUP

 !  CALL PDOMDEC

# if defined (SPHERICAL)              
     CALL CAL_CENTER      
# endif
!================================================================================

!  ALLOCATE VARIABLES AND SETUP KALMAN FILTR ASSIMILATION
#  if defined(RRK_ASSIM)
   CALL ALLOC_RRKA
   CALL SETUP_RRKA
#  endif

!
!  INITIALIZE GOTM
!
#  if defined (GOTM)
   CALL INIT_GOTM
#  endif
!
!  INITIALIZE SEDIMENT MODEL
!
#  if defined (SEDIMENT)
   IF(SEDIMENT_ON)CALL SETUP_SED(RESTART_SED)
#  endif


!
!  CALCULATE DEPTH HORIZONTAL DERIVATIVES
!
   CALL DEPTH_GRADIENT
! 
!  GENERATE SECTION INFORMATION FILES
!
   CALL SECTINF
!
!  EXCHANGE SHAPE FACTOR INFORMATION
!
# if defined (MULTIPROCESSOR)
  IF(PAR)CALL EXCHANGE(EC,MTElem,4,MYID,NPROCS,A1U,A2U) 
  IF(PAR)CALL EXCHANGE(EC,MTElem,3,MYID,NPROCS,AWX,AWY,AW0) 
  IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,ALPHA) 
  IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,ART)
# endif

# if defined (BioGen)
  IF(RESTART == 'hot_start')  CALL BIO_HOT_START
  IF(RESTART == 'cold_start') CALL BIO_INITIAL
!        CALL BIO_INITIAL
# endif


!
!  READ PARAMETERS CONTROLLING LAGRANGIAN TRACKING 
!
#  if defined (NETCDF_IO)
   CALL SET_LAG
#  endif
!
!  INITIALIZE VISUALIZATION SERVER
!

#  if defined(PV3) 
   CALL VIZ_INIT(1)
#  endif

!
!  INITIALIZE DATA ASSIMILATION VARIABLES
!
#  if defined (NG_OI_ASSIM)
   IF(CURRENT_ASSIM) CALL SET_CUR_ASSIM_DATA
   IF(SST_ASSIM)     CALL SET_SST_ASSIM_DATA
   IF(TS_ASSIM)      CALL SET_TS_ASSIM_DATA
#  endif

!
!  INITIALIZE TIME SERIES OBJECTS 
!
   CALL SET_PROBES      

!
!  INITIALIZE ICE DYNAMIC/THERMODYNAMIC
!
# if defined (ICE)
   IF(RESTART == 'cold_start') THEN
     IINT =1
     IF(C_HFX)CALL BCOND_HFX
     IINT =0
   END IF
   call ice_init_0
#endif
!
!  CALCULATE INTERVALS FOR SST DATA ASSIMILATION 
!
#  if defined (NG_OI_ASSIM)
   CALL SET_ASSIM_INTERVALS
#  endif

!===============================================================================
!  MAIN LOOP OVER ENSEMBLE KALMAN FILTER DATA ASSMILATION 
!===============================================================================
#  if defined (ENKF_ASSIM)
  
   CALL SET_ASSIM_ENKF

   DO ICYC = N1CYC, NCYC
     IINT = (ICYC-1)*DELTA_ASS
     IEND = IINT + DELTA_ASS
     DAYS = IEND*DTI/24.0_SP/3600.0_SP

     WRITE(FCYC,'(I4.4)') ICYC
     IF(ICYC>N1CYC) RESTART = 'hot_start'

     DO IENS = 1, ENKF_NENS 

       IINT = (ICYC-1)*DELTA_ASS
       IF(IBCN_GL(1) > 0) CALL READ_BC
       CALL STARTUP
       IF(ENKF_INIT == 'inv_bc'.AND. ICYC==N1CYC) CALL INI_VALS 
#  endif
!================================================================
! END OF ENSEMBLE KALMAN FILTER PART
!================================================================

!===============================================================================
!  MAIN LOOP OVER ENSEMBLE TRANSFORM KALMAN FILTER DATA ASSMILATION
!===============================================================================
#  if defined (ETKF_ASSIM)

   CALL SET_ASSIM_ETKF

   DO ICYC = N1CYC, NCYC
     IINT = (ICYC-1)*DELTA_ASS
     IEND = IINT + DELTA_ASS
     DAYS = IEND*DTI/24.0_SP/3600.0_SP

     WRITE(FCYC,'(I4.4)') ICYC
     IF(ICYC>N1CYC) RESTART = 'hot_start'

    IF(ETKF_CTR) THEN
       ETKF_RUNS = ETKF_NENS + 1
     ELSE
       ETKF_RUNS = ETKF_NENS
     ENDIF

     DO IENS = 1, ETKF_RUNS

       IINT = (ICYC-1)*DELTA_ASS
       IF(IBCN_GL(1) > 0) CALL READ_BC
       CALL STARTUP
       IF(ETKF_INIT == 'inv_bc'.AND. ICYC==N1CYC) CALL INI_VAL
#  endif
!================================================================
! END OF ENSEMBLE TRANSFORM KALMAN FILTER PART
!================================================================

!================================================================
! REDUCED KALMAN FILTER PREPARATION LOOP
!================================================================
#  if defined(RRK_PRE)
       IF (RRK_RUNCONTR > 1) THEN
2000     CONTINUE
         IEND = REF_TIME1+ DELTA_ASS
         IRESTART = REF_TIME1+ DELTA_ASS
         IRECORD  = 0
         RESTART = 'hot_start'
  
!
!  OPEN INPUT/OUTPUT FILES
!
         CALL IOFILES

!
!  INITIALIZE FLOWFIELD  --> [T,S,U,V,EL,D,Q2,Q2L]
!
         CALL STARTUP

       ENDIF
#  endif
!================================================================
! END OF REDUCED KALMAN FILTER PREPARATION LOOP
!================================================================

!
!  REPORT STATISTICS ON INITIAL VALUES
!
   CALL REPORT('INITIAL VALUE INFORMATION')



!
!  CALCULATE INTERNAL TIME STEP AND SET INTEGRATION LIMITS
!
   ISTART=IINT+1

   CALL START_CLOCK


#  if defined(RRK_PRE)   
    IF(RRK_RUNCONTR==1) GOTO 4000 
#  endif   

!
!  REPORT INTEGRATION INITIAL TIMES 
!
   IF(MSR)CALL REPORT_SIMTIME

!==================================================================
!  MAIN LOOP OVER KALMAN FILTER DATA ASSMILATION 
!==================================================================
#  if defined (RRK_ASSIM)

   DO ICYC = N1CYC, NCYC

!rrkf   Set the start and end time step within one data assimilation cycle
   ISTART = I_INITIAL + (ICYC-N1CYC)*DELTA_ASS + 1
   IEND   = ISTART + DELTA_ASS - 1

   IF(IBCN_GL(1) > 0) CALL READ_BC

#  endif
!=================================================================
! END OF REDUCED KALMAN FILTER DATA ASSIMILATION MAIN LOOP
!=================================================================

!
!  CALCULATE INTERVALS FOR SST DATA ASSIMILATION 
!
#  if defined (NG_OI_ASSIM)
   CALL SET_ASSIM_INTERVALS

   ASSIM_FLAG = 0 
   SST_CYC    = 0
!==============================================================================!
!  MAIN LOOP OVER DATA ASSIMILATION SWEEPS                                     !
!==============================================================================!
   DO IINT_SST=ISTART_DAY,IEND_DAY,IAV_DAY

   SST_CYC = SST_CYC + 1
   DO INTERNAL_SST= 1,ISWEEP
     IF(SST_ASSIM) CALL START_CLOCK
     IF(SST_ASSIM)THEN
       ISTART=INT(FLOAT(IINT_SST)*(86400./DTI)+0.5)+1
       IEND=ISTART-1+INT(IAV_DAY*86400./DTI+0.5)
       IF(MSR)THEN
         WRITE(IPT,*  )'!'
         IF(INTERNAL_SST == 1) WRITE(IPT,*)'!   SIMULATION STAGE!                 '
         IF(INTERNAL_SST == 2) WRITE(IPT,*)'!   ASSIMILATION STAGE!               '
         WRITE(IPT,7001)'! DA INTERVAL NUMBER     :  ',SST_CYC 
         WRITE(IPT,*)'!'
         CALL GETTIME(TSTRING,INT(DTI)*(ISTART-1))
         WRITE(IPT,7000)'! DA INTERVAL TIME BEGIN :  ',TSTRING
         CALL GETTIME(TSTRING,INT(DTI)*IEND)
         WRITE(IPT,7000)'! DA INTERVAL TIME ENDS  :  ',TSTRING
         WRITE(IPT,*)'!'
       END IF
     END IF

     ASSIM_FLAG = INTERNAL_SST-1 

     IF(SST_ASSIM)THEN
     IF(ASSIM_FLAG ==1)THEN
       IF(MSR)WRITE(IPT,*)'RESTARTING FROM SST RESTART FILE'
       CALL HOT_START_SST
       IF(PAR) CALL EXCHANGE_ALL
     ELSE
       IF(MSR)WRITE(IPT,*)'WRITING SST RESTART FILE'
       CALL ARC_SST
     END IF
     END IF

#    endif
       
! New Open Boundary Condition ----2
#  if defined (MEAN_FLOW)
     CALL FIND_OBSIDE
     CALL ALLOC_OBC3_DATA
     CALL SETUP_OBC3
     CALL ALLOC_OBC2_DATA
#    endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  BEGIN MAIN LOOP OVER PHYSICAL TIME STEPS                                    |
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   DO IINT=ISTART,IEND
     TIME =DTI*FLOAT(IINT)/86400.0_SP
     THOUR=DTI*FLOAT(IINT)/3600.0_SP

!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
     RAMP=1.0_SP
     IF(IRAMP /= 0) RAMP=TANH(FLOAT(IINT)/FLOAT(IRAMP))
!     IF(IRAMP /= 0) RAMP=TANH(FLOAT(IINT-ISTART+1)/FLOAT(IRAMP))

!----SET UP WATER QUALITY MODEL COEFFICIENTS-----------------------------------!
#    if defined (WATER_QUALITY)
     IF(WQM_ON)THEN
        TIME_R=MOD(IINT*DTI/3600.0_SP-14.743_SP, 24.0_SP)+6.0_SP
       CALL WQMCONST
     END IF
#    endif

#    if !defined (TWO_D_MODEL)

# if !defined(SEMI_IMPLICIT)
!----ADJUST CONSISTENCY BETWEEN 3-D VELOCITY AND VERT-AVERAGED VELOCITIES------!
# if !defined (ONE_D_MODEL)
     CALL ADJUST2D3D(1)
# endif
!----SPECIFY THE SOLID BOUNDARY CONDITIONS OF U&V INTERNAL MODES---------------!
# if defined (GCN)
     CALL BCOND_GCN(5,0)
# else
     CALL BCOND_GCY(5,0)
# endif
#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,U,V)
#    endif
# endif
! endif defined semi-implicit

!----SPECIFY THE SURFACE FORCING OF INTERNAL MODES-----------------------------!
# if defined (GCN)
	! 	 write(ipt,*)'OPENING BCOND_GCN for INTERNAL MODE'
     CALL BCOND_GCN(8,0)

# else
     CALL BCOND_GCY(8,0)
# endif

! New Open Boundary Condition ----3
#  if defined (MEAN_FLOW)
     CALL BCOND_MEANFLOW
     CALL BCOND_TIDE_3D
!     CALL BCOND_NG_3D       !change BKI
!     CALL BCOND_NG_2D       !change BKI
     CALL BCOND_BKI_3D(1)  
!     CALL BCOND_BKI_2D(4)

     CALL FLUX_OBC3D_2
     CALL FLUX_OBC2D
# endif

# endif    
!end !defined (TWO_D_MODEL)

!----SPECIFY BOTTOM FRESH WATER INPUT BOUNDARY CONDITION-----------------------!
     CALL BCOND_BFW(1)

!----SPECIFY THE BOTTOM ROUGHNESS AND CALCULATE THE BOTTOM STRESSES------------!
     CALL BOTTOM_ROUGHNESS

# if defined (HEAT_FLUX) || (ICE)
     IF(C_HFX) CALL BCOND_HFX
!     write(*,*)' called the hfx calculation' ! B Clark debug
# endif

!added by T.W. to calculate momentum sink by kelp, April, 2013
# if defined (KELP)
  if(C_KELP) call MS_KELP
# endif


# if defined (KELP)
# if defined (MULTIPROCESSOR)
     IF(PAR .AND. C_KELP) CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,EMS_X)
     IF(PAR .AND. C_KELP) CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,EMS_Y)
# endif
# endif
!finish by T.W.

# if defined (ICE)
!----ITERATE THE ICE DYNAMIC/THERMODYNAMIC MODEL ------------!

  IF(MSR)   WRITE(IPT,*)IINT, 'Calculating ICE'
!!  the thermodynamics part are modified from CICE
!!  the dynamics and advection part are use FVCOM-ice

      DO i=1,NNode
         DSW_AIR(I)=- SWRAD(I)*4.2174E3_SP*1.023E3_SP
      END DO
#    if defined(MULTIPROCESSOR)
      IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,DSW_AIR)
#    endif

     CALL ICE_1D_2D
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!
#endif

!
!==============================================================================!
!  CALCULATE DISPERSION (GX/GY) AND BAROCLINIC PRESSURE GRADIENT TERMS         !
!==============================================================================!
# if !defined (TWO_D_MODEL)

# if !defined (SEMI_IMPLICIT)
# if defined (GCN)
     CALL ADVECTION_EDGE_GCN(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !
# else
     CALL ADVECTION_EDGE_GCY(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !
# endif
# endif
! if defined semi-implicit

     IF(.NOT. BAROTROPIC)THEN                    !Barotropic Flow ?            !
!Calculate the rho mean
     IF(IRHO_MEAN > 0) THEN
       IF(MOD(IINT,IRHO_MEAN) == 0) CALL RHO_MEAN    
     END IF
     IF(C_BAROPG == 'sigma')    CALL BAROPG      !Sigma Level Pressure Gradient!
     IF(C_BAROPG == 's_levels') CALL PHY_BAROPG  !Z Level Pressure Gradient    !
     END IF                                      !                             !
                                                 !                             !
# if !defined (SEMI_IMPLICIT)

     ADX2D = 0.0_SP ; ADY2D = 0.0_SP             !Initialize GX/GY Terms       !
     DRX2D = 0.0_SP ; DRY2D = 0.0_SP             !Initialize BCPG for Ext Mode !

     DO K=1,KBM1
       DO I=1,MElem
          ADX2D(I)=ADX2D(I)+ADVX(I,K)   !*DZ1(I,K)
          ADY2D(I)=ADY2D(I)+ADVY(I,K)   !*DZ1(I,K)
          DRX2D(I)=DRX2D(I)+DRHOX(I,K)  !*DZ1(I,K)
          DRY2D(I)=DRY2D(I)+DRHOY(I,K)  !*DZ1(I,K)
        END DO
     END DO

#    if defined (GCN)
     CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
#    else
     CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
#    endif
     ADX2D = ADX2D - ADVUA                       !Subtract to Form GX
     ADY2D = ADY2D - ADVVA                       !Subtract to Form GY


!----INITIALIZE ARRAYS USED TO CALCULATE AVERAGE UA/E  OVER EXTERNAL STEPS-----!
     UARD = 0.0_SP
     VARD = 0.0_SP
     EGF  = 0.0_SP

#  if defined (EQUI_TIDE)
     EGF_EQI = 0.0_SP
#  endif
#  if defined (ATMO_TIDE)
     EGF_ATMO = 0.0_SP
#  endif       
# if defined (ICE) || (HEAT_FLUX)
     EGF_AIR = 0.0_SP
# endif

!!#    if defined (WET_DRY)
!!     UARDS = 0.0_SP
!!     VARDS = 0.0_SP
!!#    endif


     IF(IOBCN > 0) THEN
         UARD_OBCN(1:IOBCN)=0.0_SP
     END IF
# endif
! if defined semi-implicit

# endif    
!end defined (TWO_D_MODEL)

#  if defined (BALANCE_2D)
    ADVUA2    = 0.0_SP
    ADVVA2    = 0.0_SP
    ADFX2     = 0.0_SP
    ADFY2     = 0.0_SP
    DRX2D2    = 0.0_SP
    DRY2D2    = 0.0_SP
    CORX2     = 0.0_SP
    CORY2     = 0.0_SP
    PSTX2     = 0.0_SP
    PSTY2     = 0.0_SP
    ADX2D2    = 0.0_SP
    ADY2D2    = 0.0_SP
    WUSURBF2  = 0.0_SP
    WVSURBF2  = 0.0_SP 
    DUDT2     = 0.0_SP 
    DVDT2     = 0.0_SP
     
    DIVX2D2   = ZERO 
    DIVY2D2   = ZERO 
    DEDT2     = ZERO  
#  endif

! New Open Boundary Condition ----4
#  if defined (MEAN_FLOW)
     CALL ZERO_OBC3
# endif   

#  if !defined(ONE_D_MODEL) && !defined (SEMI_IMPLICIT)
!==============================================================================!
!  LOOP OVER EXTERNAL TIME STEPS                                               !
!==============================================================================!
     DO IEXT=1,ISPLIT

       !! David for VISIT
       !!=======================================================
# if defined (VISIT)
        visit_cycle=IINT-1 
        visit_time_int=DTI*FLOAT(IINT-1)/3600.0_SP
        visit_time_ext=(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT-1))/3600.0_SP
        Call VisitCheck
# endif
       !!=======================================================

       TIME  =(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT))/86400.0_SP
       THOUR1=(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT))/3600.0_SP

!
!----- USE RAMP VARIABLE TO EASE MODEL STARTUP---------------------------------!
!
!       TMP1 = FLOAT(IINT-ISTART)+FLOAT(IEXT)/FLOAT(ISPLIT)
       TMP1 = FLOAT(IINT-1)+FLOAT(IEXT)/FLOAT(ISPLIT)
       RAMP = 1.0_SP
       IF(IRAMP /= 0) RAMP = TANH(TMP1/FLOAT(IRAMP))
!
!------SURFACE BOUNDARY CONDITIONS FOR EXTERNAL MODEL--------------------------!
!
#      if defined (GCN)
	!   write(IPT,*)'Calling BCOND_GCN, external mode'
       CALL BCOND_GCN(9,0)
#      else
       CALL BCOND_GCY(9,0)
#      endif

# if defined (ICE) 
       DO I=1,MElem
         IF(ISICEC(I) == 1)THEN
       !! (magnitude of relative ocean current)*rhow*drag*aice
           DUVI = DRAGW*SQRT((U(I,1)-UICE2(I))**2+(V(I,1)-VICE2(I))**2)  ! m/s
       !! ice/ocean stress
           ice_ocnDX = DUVI*((U(I,1)*COSW-V(I,1)*SINW)-(uice2(i)*cosw-vice2(i)*sinw))
           ice_ocnDY = DUVI*((V(I,1)*COSW+U(I,1)*SINW)-(vice2(i)*cosw+uice2(i)*sinw))
           wusurf2(I)=wusurf2(I)*(1.0_SP-AIU(I))-ice_ocnDX*aiu(I)*1.0E-3
           wvsurf2(I)=wvsurf2(I)*(1.0_SP-AIU(I))-ice_ocnDY*aiu(I)*1.0E-3
         ENDIF
       END DO
#      if defined (MULTIPROCESSOR)
       IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,WUSURF2,WVSURF2)
#      endif

#endif

       CALL BCOND_BFW(2)
!
!------SAVE VALUES FROM CURRENT TIME STEP--------------------------------------!
!
       ELRK1 = EL1
       ELRK  = EL
       UARK  = UA
       VARK  = VA

#      if defined (EQUI_TIDE)
       ELRK_EQI = EL_EQI
#      endif
#      if defined (ATMO_TIDE)
       ELRK_ATMO = EL_ATMO
#      endif 
#      if defined (HEAT_FLUX)
       ELRK_AIR = EL_AIR
#      endif

! New Open Boundary Condition ----5
#  if defined (MEAN_FLOW)
       ELRKT  = ELT
       ELRKP  = ELP
       UARKNT = UANT   !change BKI
       VARKNT = VANT   !change BKI
       UARKN  = UAN    !change BKI
       VARKN  = VAN    !change BKI
#  endif 

!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!
       DO K=1,4
         TIMERK = TIME + (ALPHA_RK(K)-1.)*DTE/86400.0_SP

! New Open Boundary Condition ----6
#  if defined (MEAN_FLOW)
       CALL BCOND_TIDE_2D
!       CALL BCOND_NG_2D        ! change BKI
       CALL FLUX_OBN2D(K)
       CALL FLUX_OBC2D
#  endif 

         !FREE SURFACE AMPLITUDE UPDATE  --> ELF
         CALL EXTEL_EDGE(K)
# if defined (MULTIPROCESSOR)
         IF(PAR) CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,ELF)
# endif

#        if defined (EQUI_TIDE)
         CALL ELEVATION_EQUI
         ELF_EQI = ELRK_EQI +ALPHA_RK(K)*(ELF_EQI-ELRK_EQI) 
#        endif
#        if defined (ATMO_TIDE)
         CALL ELEVATION_ATMO
         ELF_ATMO = ELRK_ATMO +ALPHA_RK(K)*(ELF_ATMO-ELRK_ATMO) 
#        endif

#       if defined (HEAT_FLUX)
        CALL BCOND_PA_AIR 
        ELF_AIR = ELRK_AIR +ALPHA_RK(K)*(ELF_AIR-ELRK_AIR) 
#       endif


! New Open Boundary Condition ----7
#  if defined (MEAN_FLOW)
       IF (ntidenode > 0) THEN
          DO I = 1, ntidenode     ! need to calculate NEXT_OBC column of ELPF
             J = I_TIDENODE_N(I)
             ELPF(I) = ELF(J) - ELTF(I)
          END DO
       END IF

       CALL EXTELPF_EDGE(K)

       IF (IOBCN > 0) THEN
          DO I = 1, IOBCN
             J = I_OBC_N(I)
             J1= I_OBC_NODE(J)
             ELF(J) = ELTF(J1) + ELPF(J1)
          END DO
       END IF
#  else
#        if defined (GCN)
         CALL BCOND_GCN(1,K)
#        else
         CALL BCOND_GCY(1,K)
#        endif

         DO I=1,IBCN(1)
           JN = OBC_LST(1,I)
           J=I_OBC_N(JN)
           ELF(J)=ELRK(J)+ALPHA_RK(K)*(ELF(J)-ELRK(J))
         END DO
#  endif

         CALL N2E2D(ELF,ELF1)


#        if defined (WET_DRY)
         IF(WET_DRY_ON)CALL WET_JUDGE
#        endif

         CALL FLUX_OBN(K)

         !CALCULATE ADVECTIVE, DIFFUSIVE, AND BAROCLINIC MODES --> UAF ,VAF
# if defined (GCN)
         CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
# else
         CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
# endif
         CALL EXTUV_EDGE(K)
# if defined (GCN)
         CALL BCOND_GCN(2,K)
# else
         CALL BCOND_GCY(2,K)
# endif

#        if defined (MULTIPROCESSOR)
         IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,1,MYID,NPROCS,ELF)
#        endif

         !UPDATE WATER SURFACE ELEVATION
         CALL ASSIGN_ELM1_TO_ELM2

         EL  = ELF
         EL1 = ELF1

#        if defined (EQUI_TIDE)
         EL_EQI = ELF_EQI
#        endif
#        if defined (ATMO_TIDE)
         EL_ATMO = ELF_ATMO
#        endif       
#        if defined (HEAT_FLUX)
         EL_AIR = ELF_AIR
#        endif       

         !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
         CALL N2E2D(EL,EL1)

         !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
         D   = H + EL
         D1  = H1 + EL1
         UA  = UAF
         VA  = VAF
         DTFA = D

! New Open Boundary Condition ----8
#  if defined (MEAN_FLOW)
         ELT = ELTF
         UAT = UATF
         VAT = VATF
! The part below is equivalent to do both NODE MATCH and EXCHANGE for ELPF
# if defined (MULTIPROCESSOR)
         IF (ntidenode > 0) THEN
            DO I = 1, ntidenode
               J = I_TIDENODE_N(I)
               ELPF(I) = ELF(J) - ELTF(I)
            END DO
         END IF
# endif
         ELP = ELPF
         IF(K == 4 .and. nmfcell > 0)THEN
           DO I = 1, nmfcell
              J = I_MFCELL_N(I)
              OBC2D_X_TIDE(I) = OBC2D_X_TIDE(I) + UANT(I) * D1(J)
              OBC2D_Y_TIDE(I) = OBC2D_Y_TIDE(I) + VANT(I) * D1(J)
           ENDDO
         END IF
         CALL BCOND_BKI_2D(K)            ! change BKI
#  endif


         !!ENSURE ALL CELLS ARE WET IN NO FLOOD/DRY CASE  
#        if !defined (WET_DRY)
         CALL DEPTH_CHECK
#        endif

         !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
#        if defined (MULTIPROCESSOR)
         IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,UA,VA,D1)
!!#        if defined (WET_DRY)
!!         IF(PAR .AND. K==3)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,UAS,VAS)
!!#        endif
#        endif


# if !defined (TWO_D_MODEL)
         !SAVE VALUES FOR 3D MOMENTUM CORRECTION AND UPDATE
         IF(K == 3)THEN
           UARD = UARD + UA*D1
           VARD = VARD + VA*D1
           EGF  = EGF  + EL/ISPLIT

#          if defined (EQUI_TIDE)
           EGF_EQI = EGF_EQI + EL_EQI/ISPLIT
#          endif
#          if defined (ATMO_TIDE)
           EGF_ATMO = EGF_ATMO + EL_ATMO/ISPLIT
#          endif       
#          if defined (ICE) || (HEAT_FLUX)
           EGF_AIR = EGF_AIR + EL_AIR/ISPLIT
#          endif

!!#          if defined (WET_DRY)
!!           UARDS = UARDS + UAS*D1
!!           VARDS = VARDS + VAS*D1
!!#          endif
         END IF

         !CALCULATE VALUES USED FOR SALINITY/TEMP BOUNDARY CONDITIONS
         IF(K == 4.AND.IOBCN > 0) THEN
           DO I=1,IOBCN
             J=I_OBC_N(I)
             TMP=-(ELF(J)-ELRK(J))*ART1(J)/DTE-XFLUX_OBCN(I)
             UARD_OBCN(I)=UARD_OBCN(I)+TMP/FLOAT(ISPLIT)
            END DO
         END IF
# endif    
!end !defined (TWO_D_MODEL)

          !UPDATE WET/DRY FACTORS
#        if defined (WET_DRY)
         IF(WET_DRY_ON)CALL WD_UPDATE(1)
#        endif

       END DO     !! END RUNGE-KUTTA LOOP

     END DO     !! EXTERNAL MODE LOOP
!==============================================================================!
!  END LOOP OVER EXTERN STEPS                                                  !
!==============================================================================!
#  endif
!end !defined (ONE_D_MODEL) or !defined (SEMI_IMPLICIT)

#  if defined (ONE_D_MODEL)
   EL = EL - DTI*(QEVAP3-QPREC3)*ROFVROS
   CALL N2E2D(EL,EL1)
   D    = H + EL
   D1   = H1 + EL1
   DTFA = D
#  endif   

# if defined (SEMI_IMPLICIT)
   !! David for VISIT
   !!=======================================================
#  if defined (VISIT)
   visit_cycle=IINT-1
   visit_time_int=DTI*FLOAT(IINT-1)/3600.0_SP
!  visit_time_ext=(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT-1))/3600.0_SP
   Call VisitCheck
#  endif
   !!=======================================================

#  if defined (EQUI_TIDE)
   CALL ELEVATION_EQUI
#  endif
#  if defined (ATMO_TIDE)
   CALL ELEVATION_ATMO
#  endif
#  if defined (HEAT_FLUX)
   CALL BCOND_PA_AIR 
#  endif

#  if defined (GCN)
   CALL BCOND_GCN(1,0)
#  else
   CALL BCOND_GCY(1,0)
#  endif

# endif
! if defined semi-implicit

#  if !defined (TWO_D_MODEL)
!==============================================================================!
!    ADJUST INTERNAL VELOCITY FIELD TO CORRESPOND TO EXTERNAL                  !
!==============================================================================!
# if !defined (SEMI_IMPLICIT)
# if !defined (ONE_D_MODEL)
     CALL ADJUST2D3D(2)
# endif
# endif
! if defined semi-implicit

! New Open Boundary Condition ----9
#  if defined (MEAN_FLOW)
!     CALL BCOND_NG_3D       ! change BWI
     CALL BCOND_BKI_3D(2)      
     CALL FLUX_OBC3D
#  endif

# if defined (ICE)
       !  ggao 0521 for ice ocean couple
       DO I=1,MElem
         IF(ISICEC(I)==1) THEN
       !! (magnitude of relative ocean current)*rhow*drag*aice
           DUVI = DRAGW*SQRT((U(I,1)-UICE2(I))**2+(V(I,1)-VICE2(I))**2)  ! m/s
       !! ice/ocean stress
           ice_ocnDX = DUVI*((U(I,1)*COSW-V(I,1)*SINW)-(uice2(i)*cosw-vice2(i)*sinw))
           ice_ocnDY = DUVI*((V(I,1)*COSW+U(I,1)*SINW)-(vice2(i)*cosw+uice2(i)*sinw))
           wusurf(I) =wusurf(I) *(1.0_SP-AIU(I)) +ice_ocnDX*aiu(I)*1.0E-3
           wvsurf(I) =wvsurf(I) *(1.0_SP-AIU(I)) +ice_ocnDY*aiu(I)*1.0E-3
         ENDIF
       END DO
#      if defined (MULTIPROCESSOR)
       IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,WUSURF,WVSURF)
#      endif
# endif

# endif
!if !defined (TWO_D_MODEL)

!==============================================================================!
!     CALCULATE INTERNAL VELOCITY FLUXES                                       |
!==============================================================================!
# if !defined (SEMI_IMPLICIT)

# if !defined (TWO_D_MODEL)
                          !                                                    !
     CALL VERTVL_EDGE     ! Calculate/Update Sigma Vertical Velocity (Omega)   !
#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL WD_UPDATE(2)
#    endif

#    if defined (GCN)
     CALL ADV_UV_EDGE_GCN ! Horizontal Advect/Diff + Vertical Advection        !
#    else
     CALL ADV_UV_EDGE_GCY ! Horizontal Advect/Diff + Vertical Advection        !
#    endif
     CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !

     IF(ADCOR_ON) THEN
       CALL ADCOR
       CALL VDIF_UV       ! Implicit Integration of Vertical Diffusion of U/V  !
     ENDIF

#    if defined (WET_DRY)
     DO I=1,MElem
       IF(H1(I) <= DJUST ) THEN
         DO K=1,KBM1
           UF(I,K)=UA(I)
           VF(I,K)=VA(I)
         END DO
       END IF
     END DO
#    endif

#    if defined (GCN)
     CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input           !
#    else
     CALL BCOND_GCY(3,0)    ! Boundary Condition on U/V At River Input           !
#    endif

# endif
! if !defined (TWO_D_MODEL)

# else

#    if defined (TWO_D_MODEL) 

#    if defined (GCN)
     CALL BCOND_GCN(8,0)
#    else
     CALL BCOND_GCY(8,0)
#    endif

#    if defined (GCN)
     CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
#    else
     CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
#    endif     
     CALL UV2D_SBD
     CALL SEMI_IMPLICIT_EL
     IF(ADCOR_ON) THEN
       CALL ADCOR
       CALL SEMI_IMPLICIT_EL
     ENDIF

     UA  = UAF
     VA  = VAF
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,UA,VA)
#    endif

#    else

#    if defined (GCN)
     CALL ADV_UV_EDGE_GCN ! Horizontal Advect/Diff + Vertical Advection        !
#    else
     CALL ADV_UV_EDGE_GCY ! Horizontal Advect/Diff + Vertical Advection        !
#    endif
     CALL UV3D_SBD
     CALL SEMI_IMPLICIT_EL
     CALL VDIF_UV
     IF(ADCOR_ON) THEN
       CALL ADCOR
       CALL SEMI_IMPLICIT_EL
       CALL VDIF_UV
     ENDIF
#    if defined (SPHERICAL) && (NORTHPOLE)
     DO J=1, NP
       I = NP_LST(J)
       IF(CELL_NORTHAREA(I) == 1)THEN
         DO K=1, KBM1
           UTMP=UF(I,K)
           VTMP=VF(I,K)
           UF(I,K) = VTMP*COS(XC(I)*DEG2RAD)-UTMP*SIN(XC(I)*DEG2RAD)
           VF(I,K) = -( UTMP*COS(XC(I)*DEG2RAD)+VTMP*SIN(XC(I)*DEG2RAD) )
         ENDDO
       ENDIF
     ENDDO
#    endif
     CALL UPDATE_EL

#    if defined (WET_DRY)
     DO I=1, MTElem
       IF(H1(I) <= DJUST ) THEN
         DO K=1, KBM1
           UF(I,K)=UA(I)
           VF(I,K)=VA(I)
         ENDDO
       ENDIF
     ENDDO
!!     US = U
!!     VS = V
#    endif
     
     CALL VERTVL_EDGE
#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL WD_UPDATE(2)
#    endif

     DO I=1,IOBCN            ! do linear T/S case, must redo this part, carefule!
       J=I_OBC_N(I)
       UARD_OBCN(I)=-(EL(J)-ET(J))*ART1(J)/DTI-XFLUX_OBCN(I)
     END DO

#    endif
!if defined (TWO_D_MODEL) 

#    if defined (EQUI_TIDE)
     EL_EQI = ELF_EQI
#    endif
#    if defined (ATMO_TIDE)
     EL_ATMO = ELF_ATMO
#    endif
#    if defined (HEAT_FLUX)
     EL_AIR = ELF_AIR
#    endif

#    if !defined (WET_DRY)
     CALL DEPTH_CHECK
#    endif

# endif
!if !defined (SEMI_IMPLICIT)

# if !defined (TWO_D_MODEL)
           
     CALL WREAL           ! Calculate True Vertical Velocity (W)               !
                          !                                                    !
     CALL VISCOF_H        ! Calculate horizontal diffusion coefficient for     !
                          ! the scalar                                         !
!==============================================================================!
!    TURBULENCE MODEL SECTION                                                  |
!==============================================================================!
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,WUSURF,WVSURF)
     IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,WUBOT,WVBOT)
#    endif

     IF(VERTMIX == 'closure')THEN
     !=================General Ocean Turbulence Model==========================!
#    if defined (GOTM)

        ! There is a bug in the advection of turbulent kinetic energy
        ! when using the GOTM MODULE

!     CALL ADV_Q(TKE,TKEF)     !!Advection of Tubulent Kinetic Energy
!     CALL ADV_Q(TEPS,TEPSF)   !!Advection of Turbulent Dissipation Rate 
!#    if defined(MULTIPROCESSOR)
!     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,TKEF)
!     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,TEPSF)
!#    endif
!     TKE  = TKEF
!     TEPS = TEPSF
     CALL ADVANCE_GOTM            !!Solve TKE/TEPS eqns for KH/KM/KQ

#    else     
     !===================Original FVCOM MY-2.5/Galperin 1988 Model=============!
#    if !defined (ONE_D_MODEL)
     CALL ADV_Q(Q2,Q2F)       !!Advection of Q2 
     CALL ADV_Q(Q2L,Q2LF) 
#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,Q2F)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,Q2LF)
#    endif
     IF(TS_FCT) CALL FCT_Q2                              !Conservation Correction   !
     IF(TS_FCT) CALL FCT_Q2L                             !Conservation Correction   !
#    endif
!end !defined (ONE_D_MODEL)
     CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange   !
#    endif
     Q2  = Q2F
     Q2L = Q2LF

#    endif
     ELSE
       KM = UMOL
       KH = UMOL*VPRNU
     END IF  

#   if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,KM,KQ,KH)
#   endif
    CALL N2E3D(KM,KM1)
!==============================================================================!
!    SEDIMENT MODEL SECTION                                                    |
!==============================================================================!
#    if defined (SEDIMENT)
     ALLOCATE(FTEMP(0:NTNode),FTEMP2(0:MTElem))
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,WUBOT,WVBOT)
#    endif
     FTEMP2 = SQRT(WUBOT**2 + WVBOT**2)
     CALL E2N2D(FTEMP2,FTEMP)
     IF(SEDIMENT_ON)CALL ADVANCE_SED(DTI,THOUR*3600,FTEMP)
     DEALLOCATE(FTEMP,FTEMP2)
#    endif
#    if defined (BioGen)
     CALL BIO_3D1D                          
#    endif
!==============================================================================!
!    UPDATE TEMPERATURE IN NON-BAROTROPIC CASE                                 !
!==============================================================================!
     IF(TEMP_ON)THEN 
     
     CALL ADV_T                                     !Advection                 !
#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,TF1)
#    endif

!#                                                   if !defined (DOUBLE_PRECISION)
     IF(TS_FCT) CALL FCT_T                             !Conservation Correction   !
!#                                                   endif

     IF(CASENAME == 'gom')THEN
       CALL VDIF_TS_GOM(1,TF1)
     ELSE  
       CALL VDIF_TS(1,TF1)                            !Vertical Diffusion        !
     END IF

#                                                   if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,TF1) !Interprocessor Exchange   !
#                                                   endif
     CALL BCOND_TS(1)                               !Boundary Conditions       !
#    if defined (NG_OI_ASSIM)
!    CALL TEMP_NUDGING                              !Nudge Temperature         !
     IF(TS_ASSIM)THEN
      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL TEMP_ASSIMILATION
      IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL TEMP_ASSIMILATION
     END IF
#    endif

     T1 = TF1                                       !Update to new time level  !
     CALL N2E3D(T1,T)                               !Shift to Elements         !

     END IF                                         !                          !
!==============================================================================!
!    UPDATE SALINITY IN NON-BAROTROPIC CASE                                    !
!==============================================================================!
     IF(SALINITY_ON)THEN                            !                          !   

     CALL ADV_S                                     !Advection                 !
#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,SF1)
#    endif

!#                                                   if !defined (DOUBLE_PRECISION)
     IF(TS_FCT) CALL FCT_S                             !Conservation Correction   !
!#                                                   endif

     CALL VDIF_TS(2,SF1)                            !Vertical Diffusion        !
#                                                   if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,SF1) !Interprocessor Exchange   !
#                                                   endif
     CALL BCOND_TS(2)                               !Boundary Conditions       !
#    if defined (NG_OI_ASSIM)
!    CALL SALT_NUDGING                              !Nudge Salinity            !
     IF(TS_ASSIM)THEN
      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL SALT_ASSIMILATION
      IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL SALT_ASSIMILATION
     END IF

#    endif

     S1 = SF1                                       !Update to new time level  !
     CALL N2E3D(S1,S)                               !Shift to Elements         !

     END IF                                         !                          !
!==============================================================================!
#  if defined (DYE_RELEASE)

!==============================================================================!
!    UPDATE DYE IN NON-BAROTROPIC CASE                                         !
!==============================================================================!
!     IF(DYE_ON.AND.IINT.GE.IINT_SPE_DYE_B) THEN     !                          !                          !
     IF(DYE_ON) THEN     !                          !                          !
     CALL ADV_DYE                                   !Advection                 !

#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,DYEF)
#    endif

#                                                   if !defined (DOUBLE_PRECISION)
!     IF(TS_FCT) CALL AVER_DYE                      !Conservation Correction   !
#                                                   endif
     CALL VDIF_DYE(DYEF)                            !Vertical Diffusion        !

!check
!!     DYE = DYEF                                     !Update to new time level  !
!!     IF(IINT.GE.IINT_SPE_DYE_B) CALL ARCHIVE
!check    
#                                                   if defined (MULTIPROCESSOR)
    IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,DYEF) !Interprocessor Exchange   !
#                                                   endif
     CALL BCOND_DYE                                 !Boundary Conditions       !
     DYE = DYEF                                     !Update to new time level  !

!!     IF(MSR) WRITE(IPT,*) 'CALL Dye_on--iint=',iint,IINT_SPE_DYE_B
     END IF                                         !                          !
!==================================================================================!
#  endif

#    if !defined (MPDATA)
     IF(POINT_ST_TYPE == 'calculated')THEN
!    ADJUST TEMPERATURE AND SALINITY AT RIVER MOUTHS
       CALL ADJUST_TS
     END IF  
#    endif   

#    if defined (WATER_QUALITY)
!==============================================================================!
!    CALCULATE WATER QUALITY VARIABLES CONCENTRATIONS                          |
!==============================================================================!
                                                    !                          !
     IF(WQM_ON)THEN                                 !Water Quality Active?     !
     CALL ADV_WQM                                   !Advection                 !
     CALL VDIF_WQM(WQM_F)                           !Vertical Diffusion        !
     CALL EXCHANGE_WQM                              !Interprocessor Exchange   !
     CALL BCOND_WQM                                 !Boundary Conditions       !
     WQM(1:NNode,1:KBM2,1:NB) = WQM_F(1:NNode,1:KBM2,1:NB)  !Update                    !
     END IF                                         !                          !
!==============================================================================!
#    endif

#    if defined (NG_OI_ASSIM)
!==============================================================================!
!     ASSIMILATE SEA SURFACE TEMPERATURE DATA                                  |
!==============================================================================!
     IF(SST_ASSIM .AND. ASSIM_FLAG == 1) CALL SST_ASSIMILATION !NUDGING OR OI SST !
     IF(SST_ASSIM .AND. ASSIM_FLAG == 0) CALL SST_INT     !STORE HOURLY SST    !
     IF(SST_ASSIM .AND. ASSIM_FLAG == 1) CALL N2E3D(T1,T) !RECALCULATE T       !

!
!----Recalculate Element Based Temperatures Based on Nudged or OIed Temp Field T1------!
!
     IF(SST_ASSIM .AND. ASSIM_FLAG == 1)CALL N2E3D(T1,T)

!==============================================================================!
#    endif

!==============================================================================!
!     UPDATE THE DENSITY IN NON-BAROTROPIC CASE                                |
!==============================================================================!
     IF(.NOT.BAROTROPIC)THEN
       IF(CTRL_DEN == 'pdensity'   ) CALL DENS
       IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
       IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3
     END IF
!==============================================================================!
!     MIMIC CONVECTIVE OVERTURNING TO STABILIZE VERTICAL DENSITY PROFILE       |
!==============================================================================!

     IF(VERT_STAB)THEN
       CALL CONV_OVER
       IF(.NOT.BAROTROPIC)THEN
         IF(CTRL_DEN == 'pdensity'   ) CALL DENS
         IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
         IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3
       END IF
     END IF  
#  endif     
! end if !defined (TWO_D_MODEL)

#    if defined (NG_OI_ASSIM)
!==============================================================================!
!     DATA ASSIMILATION FOR CURRENT FIELD                                      |
!==============================================================================!
     IF(CURRENT_ASSIM)THEN
       IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL CURRENT_ASSIMILATION
       IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL CURRENT_ASSIMILATION
# if defined (GCN)
       CALL BCOND_GCN(5,0)
# else
       CALL BCOND_GCY(5,0)
# endif
     END IF
!==============================================================================!
#    endif

!==============================================================================!
!    LAGRANGIAN PARTICLE TRACKING                                              |
!==============================================================================!
#  if defined (NETCDF_IO)
     CALL LAG_UPDATE
#  endif     
#    if !defined (TWO_D_MODEL)
!==============================================================================!
!     UPDATE VELOCITY FIELD (NEEDED TO WAIT FOR SALINITY/TEMP/TURB/TRACER)     |
!==============================================================================!
     U = UF
     V = VF

!===== by TW to update velocity field by blocking designated cells ============!
#if defined V_BLOCK	
!===== this method is a shortcut and thus not rigorous ========================!
!     IF(BLOCK) THEN
!       if(block_count) then
!          block_count = .FALSE.
!          allocate(kount(0:MTElem))
!          kount = 0 
!          do i=1,MTElem
!             do j=1,kbm1
!                if(CELL_BLOCK(I,J)==1) kount(i) = kount(i) + 1      !count the blocked vertical layer#
!             end do
!          end do
!        end if

!       do i=1,MElem   !update current total water depth for external mode vertical integration, may be wrong
!          D1(i)=D1(i)*(KBM1-kount(i))/KBM1
!       end do
!     END IF

     IF(BLOCK) THEN
       DO J=1,KBM1
          DO I=1,MElem
             IF(CELL_BLOCK(I,J)==1) THEN
               U(I,J) = 0.
               V(I,J) = 0.

!               W(I,J) = 0.
!               W(I,J+1) = 0.
             END IF
          END DO
       END DO
     END IF
#endif
	 
!===== end velocity blockage ==================================================!
!==============================================================================!
!    PERFORM DATA EXCHANGE FOR ELEMENT BASED INFORMATION AT PROC BNDRIES       |
!==============================================================================!

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,U,V)
#      if defined (GOTM)
       CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,TKE,TEPS)
#      else 
       CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,Q2,Q2L)
#      endif
       CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,RHO,T,S)
       CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,S1,T1,RHO1)

#   if defined (DYE_RELEASE)
       CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,DYE)
#   endif       

     END IF
#   endif
	 
#   if defined (OUTPUT_FLUX)
	  IF(WRITE_TCE_FLUX .AND. Flux_INT /=0 )THEN
	     IF(MOD(IINT,Flux_INT)==0) CALL WRITE_EDGELINE_FLUX			
	  ENDIF
#   endif	 
	 
	 
!==============================================================================!
!     PERFORM DATA EXCHANGE FOR WATER QUALITY VARIABLES                        |
!==============================================================================!
#    if defined (WATER_QUALITY)
     CALL EXCHANGE_WQM
#    endif
#    endif      
! end if !defined (TWO_D_MODEL)

!
!----SHIFT SEA SURFACE ELEVATION AND DEPTH TO CURRENT TIME LEVEL---------------!
!
     ET  = EL  
     DT  = D 
     ET1 = EL1
     DT1 = D1
     
#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL WD_UPDATE(3)
#    endif

! New Open Boundary Condition ----10
#  if defined (MEAN_FLOW)
     ELTDT = ELT
#  endif

!==============================================================================!
!    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
!==============================================================================!

      IF(MSR)CALL REPORT_TIME(IINT,ISTART,IEND,TIME*86400,IPT) 
     IF(MOD(IINT,IREPORT)==0) THEN
       CALL REPORT("FLOW FIELD STATS")
     !  IF(MSR) CALL REPORT_TIME(IINT,ISTART,IEND,TIME*86400,IPT)     
     ENDIF
      CALL DUMP_PROBE_DATA 

#    if !defined (AIX) && !defined(ABSOFT)
     CALL FLUSH(6)
#    endif
!
!-------------UPDATE THE VISUALIZATION SERVER----------------------------------!
!
#    if defined (PV3)
     CALL PV_UPDATE(THOUR)
#    endif

#    if !defined (NG_OI_ASSIM)
     CALL ARCHIVE
#    endif

! New Open Boundary Condition ----11
#    if defined (TIDE_OUTPUT) && !defined (MEAN_FLOW) 
       CALL PRINT_VALS  !Wen Long un-commented this line 
#    endif

# if defined (VISIT)

     IF (VISIT_CMD_DUMP) THEN
        VISIT_CMD_DUMP=.FALSE. ! Turn it off for next time
        
        IF(MSR)WRITE(IPT,*)  '!  VISIT CMD => DUMPING  :    RESTART FILE'
        CALL ARCRST
        
# if defined(WET_DRY)
        
        write(CH8,'(I8.8)') IINT
        FNAME = 're_'//trim(CH8)//'_wd'
        
        IF(WET_DRY_ON) CALL WD_DUMP(FNAME) 
        
# endif 
        
# if defined(SEDIMENT)
        IF(SEDIMENT_ON) CALL ARCHIVE_SED
# endif       
        
     END IF
     
# endif
  
  
!
!---------------WRITE OUTPUT FILES---------------------------------------------!
!
#    if defined (NG_OI_ASSIM)
     IF((SST_ASSIM .AND. ASSIM_FLAG==1) .OR. PURE_SIM) CALL  ARCHIVE
#    endif
  
     CALL SHUTDOWN_CHECK

   END DO !!MAIN LOOP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    END MAIN LOOP OVER PHYSICAL TIME STEPS                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#  if defined (NG_OI_ASSIM)
   IF(SST_ASSIM .AND. ASSIM_FLAG==1) CALL ARCRST
#  endif

#  if defined (NG_OI_ASSIM)
   END DO
   END DO
#  endif
!==============================================================================!
!  END MAIN LOOP OVER DATA ASSIMILATION INTERVALS AND SWEEP #                  !
!==============================================================================!

!================================================================================!
!  END OF PROPAGATING OF THE PREVIOUS ANALYSIS TO GET THE FORCAST AT ASSIMILATION
!  START OF UPDATING THE FORCAST BY THE RRKF
!================================================================================!
#  if defined (RRK_ASSIM) 
 
     IF(IEND>=RRK_START .AND. IEND<=RRK_END) CALL RRK_RRKF

   ENDDO

!   IF(MSR) CALL RRK_FNLOUT  
#  endif

!===============================================================================
!  MAIN LOOP OVER ENSEMBLE KALMAN FILTERS DATA ASSMILATION 
!===============================================================================
#  if defined (ENKF_ASSIM)

        CALL ARCRST
# if defined(WET_DRY)
        write(CH4,'(I4.4)') iens
        FNAME = TRIM(OUTDIR)//'/fct/restart'//trim(CH4)//'_wd.dat'
        IF(WET_DRY_ON) CALL WD_DUMP(FNAME) 
# endif
     ENDDO


!     IF(MSR) THEN
       IF(IEND>=ENKF_START .AND. IEND<=ENKF_END) THEN
          CALL ENKF_ASS      
       ENDIF 
!     ENDIF

#    if defined (MULTIPROCESSOR)
        IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#    endif

   ENDDO

#  endif

!==============================================================================!
!  CLOSE UP COMPUTATION                                                        !
!==============================================================================!
   IINT = IEND
   CALL CLOSEFILES
   
!W. Long close files and deallocate variable for edgelien fux
# if defined (OUTPUT_FLUX)
  IF(WRITE_TCE_FLUX)THEN
	CALL DEALLOC_EDGELINE_FLUX
  ENDIF
# endif
   
#  if defined(ENKF_ASSIM) 

#  elif defined(ETKF_ASSIM)

#  elif defined(RRK_PRE)

#  elif defined(RRK_ASSIM)

#  else
   CALL ARCRST
   IF(MSR)THEN
     WRITE(IPT,*)'DUMPING RESTART'
   END IF
#  if defined (WET_DRY)
   IF(WET_DRY_ON) THEN
      write(CH8,'(I8.8)') IINT
      FNAME = 're_'//trim(CH8)//'_wd'
      CALL WD_DUMP(FNAME)
   ENDIF
#  endif
#  endif

   CALL REPORT('FINAL VALUES INFORMATION')

   IF(MSR)THEN
     WRITE(IPT,*) ; WRITE(IPT,*)'Computation completed, congratulations!'
     CALL GET_CLOCK 
     CALL GETTIME(TSTRING,INT(TCURRENT-TINIT))
   END IF

4000 continue

#  if defined (RRK_PRE)
   IF(MSR) THEN

    IF(RRK_RUNCONTR == 2) THEN  

      WRITE(RRKINP1,'(I8.8)') REF_TIME1+ DELTA_ASS
#     if !defined (DOS)
        RRKFILES = "mv -f "//"re_"//trim(RRKINP1)//" "//TRIM(OUTDIR)//"/rrktemp/basfct"  
#       if !defined (CRAY)
           CALL SYSTEM( TRIM(RRKFILES) )  
#       endif
#       if defined (CRAY)
           CALL CRAY_SYSTEM_CALL(TRIM(RRKFILES)) 
#       endif             
#     endif         

#     if defined(WET_DRY)
        IF(WET_DRY_ON) THEN
#         if !defined (DOS)
          RRKFILES = "mv -f "//"re_"//trim(RRKINP1)//"_wd "//TRIM(OUTDIR)//"/rrktemp/basfct_wd"  
#         if !defined (CRAY)
             CALL SYSTEM( TRIM(RRKFILES) )  
#         endif
#         if defined (CRAY)
             CALL CRAY_SYSTEM_CALL(TRIM(RRKFILES)) 
#         endif             
#         endif
        ENDIF
#     endif 
      RRK_RUNCONTR = RRK_RUNCONTR + 1
      GOTO 3000

    ELSE IF(RRK_RUNCONTR > 2 .AND. RRK_RUNCONTR <= RRK_NEOF+2) THEN

      WRITE(RRKINP1,'(I8.8)') REF_TIME1+ DELTA_ASS
      WRITE(RRKINP2,'(I4.4)') RRK_RUNCONTR -2
#     if !defined (DOS)
        RRKFILES = "mv -f "//"re_"//trim(RRKINP1)//" "//TRIM(OUTDIR)//"/rrktemp/eoffct"//TRIM(RRKINP2)  
#       if !defined (CRAY)
           CALL SYSTEM( TRIM(RRKFILES) )  
#       endif
#       if defined (CRAY)
           CALL CRAY_SYSTEM_CALL(TRIM(RRKFILES)) 
#       endif             
#     endif         

#     if defined(WET_DRY)
        IF(WET_DRY_ON) THEN
#         if !defined (DOS)
          RRKFILES = "mv -f "//"re_"//trim(RRKINP1)//"_wd "//TRIM(OUTDIR)//"/rrktemp/eoffct"//TRIM(RRKINP2)//"_wd"  
#         if !defined (CRAY)
             CALL SYSTEM( TRIM(RRKFILES) )  
#         endif
#         if defined (CRAY)
             CALL CRAY_SYSTEM_CALL(TRIM(RRKFILES)) 
#         endif             
#         endif
        ENDIF
#     endif 

      RRK_RUNCONTR = RRK_RUNCONTR + 1
      GOTO 3000
    
    ENDIF

   ENDIF

     RRK_RUNCONTR = RRK_RUNCONTR + 1
     
     IF(MSR) WRITE(IPT,*) 
     IF(MSR) WRITE(IPT,*) 'Starting Kalman filter input data processing......'

!     IF(RRK_RUNCONTR >=2) goto 6000

!  CALL ref.F PROGRAM----
     IF(MSR) CALL RRK_REF(1)
     IF(MSR) CALL RRK_REF(2)

!  CALL eof.F PROGRAM----
     IF(MSR) CALL RRK_EOF

3000 CONTINUE

#    if defined (MULTIPROCESSOR)
        IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#    endif

!  CALL fvcomrrK.F PROGRAM-------

     IF(RRK_RUNCONTR ==2) CALL RRK_RRK(1) 

     IF(RRK_RUNCONTR > 1 .AND. RRK_RUNCONTR <= 2+ RRK_NEOF) GOTO 2000

     CALL RRK_RRK(2)

     IF(MSR) CALL RRK_RRK(4)

     IF(MSR) WRITE(IPT,*)
     IF(MSR) WRITE(IPT,*) 'Finish Kalman filter data processing, congratulation! '

#  endif   

#  if defined (SEMI_IMPLICIT)
   CALL PETSc_CLEANUP
#  endif

   CALL PSTOP

!
!----------------------FORMAT STATEMENTS---------------------------------------|
!
7000 FORMAT(1X,A28,A13)  
7001 FORMAT(1X,A28,I8)  

   END PROGRAM US_FVCOM
!==============================================================================!
