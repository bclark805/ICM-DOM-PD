!==============================================================================|
!    Allocate and Initialize Most Arrays                                       !
!==============================================================================|

   SUBROUTINE ALLOC_VARS

!==============================================================================!
   USE ALL_VARS
   IMPLICIT NONE

   INTEGER NCT,NDB,IERR
   REAL(SP) :: MEMTOT,MEMCNT
!==============================================================================!

   NDB = 1       !!GWC BASE THIS ON KIND
   NCT = MTElem*3
   MEMCNT = 0.

   IF(MSR)WRITE(IPT,*  )'!'
   IF(MSR)WRITE(IPT,*)'!                ALLOCATING MEMORY    '
   IF(MSR)WRITE(IPT,*  )'!'

!==============================================================================!
!  ALLOCATE:                                                                   !
!==============================================================================!

!--------------------------Grid Metrics---------------------------------------------!

   ALLOCATE(XC(0:MTElem))            ;XC   = ZERO   !!X-COORD AT FACE CENTER 
   ALLOCATE(YC(0:MTElem))            ;YC   = ZERO   !!Y-COORD AT FACE CENTER
   ALLOCATE(VX(0:NTNode))            ;VX   = ZERO   !!X-COORD AT GRID POINT
   ALLOCATE(VY(0:NTNode))            ;VY   = ZERO   !!Y-COORD AT GRID POINT
   ALLOCATE(ART(0:MTElem))           ;ART  = ZERO   !!AREA OF ELEMENT
   ALLOCATE(ART1(0:NTNode))          ;ART1 = ZERO   !!AREA OF NODE-BASE CONTROl VOLUME
   ALLOCATE(ART2(NTNode))            ;ART2 = ZERO   !!AREA OF ELEMENTS AROUND NODE

   MEMCNT = NTNode*7*NDB + MEMCNT

!----------------Node, Boundary Condition, and Control Volume-----------------------!

   ALLOCATE(NV(0:MTElem,4))          ;NV       = 0  !!NODE NUMBERING FOR ELEMENTS
   ALLOCATE(NBE(0:MTElem,3))         ;NBE      = 0  !!INDICES OF ELMNT NEIGHBORS
   ALLOCATE(NTVE(0:NTNode))          ;NTVE     = 0 
   ALLOCATE(NTSN(NTNode))            ;NTSN     = 0 
   ALLOCATE(ISONB(0:NTNode))         ;ISONB    = 0  !!NODE MARKER = 0,1,2
   ALLOCATE(ISBCE(0:MTElem))         ;ISBCE    = 0 
   ALLOCATE(NIEC(NCT,2))         ;NIEC     = 0
   ALLOCATE(NTRG(NCT))           ;NTRG     = 0
   ALLOCATE(DLTXE(NCT))          ;DLTXE    = ZERO
   ALLOCATE(DLTYE(NCT))          ;DLTYE    = ZERO
   ALLOCATE(DLTXYE(NCT))         ;DLTXYE   = ZERO
   ALLOCATE(SITAE(NCT))          ;SITAE    = ZERO
   ALLOCATE(XIJE(NCT,2))         ;XIJE     = ZERO
   ALLOCATE(YIJE(NCT,2))         ;YIJE     = ZERO 
   ALLOCATE(N_ICELLQ(NUMQBC,2))  ;N_ICELLQ = 0

                 MEMCNT = MTElem*8 + NTNode*3  + NCT*3 + NCT*7*NDB  + MEMCNT

!----------------2-d arrays for the general vertical coordinate -------------------------------!

   ALLOCATE(Z(0:NTNode,KB))               ; Z    = ZERO    !!SIGMA COORDINATE VALUE 
   ALLOCATE(ZZ(0:NTNode,KB))              ; ZZ   = ZERO    !!INTRA LEVEL SIGMA VALUE
   ALLOCATE(DZ(0:NTNode,KB))              ; DZ   = ZERO    !!DELTA-SIGMA VALUE
   ALLOCATE(DZZ(0:NTNode,KB))             ; DZZ  = ZERO    !!DELTA OF INTRA LEVEL SIGMA 
   ALLOCATE(Z1(0:MTElem,KB))              ; Z1   = ZERO    !!SIGMA COORDINATE VALUE 
   ALLOCATE(ZZ1(0:MTElem,KB))             ; ZZ1  = ZERO    !!INTRA LEVEL SIGMA VALUE
   ALLOCATE(DZ1(0:MTElem,KB))             ; DZ1  = ZERO    !!DELTA-SIGMA VALUE
   ALLOCATE(DZZ1(0:MTElem,KB))            ; DZZ1 = ZERO    !!DELTA OF INTRA LEVEL SIGMA 
   ALLOCATE(S_HYBRID(0:NTNode,KB))	  ; S_HYBRID = 0.0  !! Hybrid sigma coordinate output, B Clark 2015  	

                   MEMCNT = NTNode*KB*4*NDB + MTElem*KB*4*NDB +MEMCNT

!---------------2-d flow variable arrays at elements-------------------------------!

   ALLOCATE(UA(0:MTElem))            ;UA        = ZERO  !!VERTICALLY AVERAGED X-VELOC
   ALLOCATE(VA(0:MTElem))            ;VA        = ZERO  !!VERTICALLY AVERAGED Y-VELOC
   ALLOCATE(UAF(0:MTElem))           ;UAF       = ZERO  !!UA FROM PREVIOUS RK STAGE 
   ALLOCATE(VAF(0:MTElem))           ;VAF       = ZERO  !!VA FROM PREVIOUS RK STAGE 
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(UARK(0:MTElem))          ;UARK      = ZERO  !!UA FROM PREVIOUS TIMESTEP 
   ALLOCATE(VARK(0:MTElem))          ;VARK      = ZERO  !!VA FROM PREVIOUS TIMESTEP 
   ALLOCATE(UARD(0:MTElem))          ;UARD      = ZERO  !!UA AVERAGED OVER EXTERNAL INT
   ALLOCATE(VARD(0:MTElem))          ;VARD      = ZERO  !!VA AVERAGED OVER EXTERNAL INT
#  endif
   ALLOCATE(COR(0:MTElem))           ;COR       = ZERO  !!CORIOLIS PARAMETER
   ALLOCATE(H1(0:MTElem))            ;H1        = ZERO  !!BATHYMETRIC DEPTH   
   ALLOCATE(D1(0:MTElem))            ;D1        = ZERO  !!DEPTH
   ALLOCATE(DT1(0:MTElem))           ;DT1       = ZERO  !!DEPTH  
   ALLOCATE(EL1(0:MTElem))           ;EL1       = ZERO  !!SURFACE ELEVATION
   ALLOCATE(ELF1(0:MTElem))          ;ELF1      = ZERO  !!SURFACE ELEVATION
   ALLOCATE(DTFA(0:NTNode))          ;DTFA      = ZERO  !!ADJUSTED DEPTH FOR MASS CONSERVATION
   ALLOCATE(ET1(0:MTElem))           ;ET1       = ZERO  !!SURFACE ELEVATION
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK1(0:MTElem))         ;ELRK1     = ZERO  !!SURFACE ELEVATION
#  endif
   ALLOCATE(CC_SPONGE(0:MTElem))     ;CC_SPONGE = ZERO  !!SPONGE DAMPING COEFFICIENT FOR MOMENTUM
                 MEMCNT = MTElem*17*NDB + MEMCNT

!---------------2-d flow variable arrays at nodes----------------------------------!

   ALLOCATE(H(0:NTNode))             ;H    = ZERO       !!BATHYMETRIC DEPTH   
   ALLOCATE(D(0:NTNode))             ;D    = ZERO       !!DEPTH   
   ALLOCATE(DT(0:NTNode))            ;DT   = ZERO       !!DEPTH   
   ALLOCATE(EL(0:NTNode))            ;EL   = ZERO       !!SURFACE ELEVATION
   ALLOCATE(ELF(0:NTNode))           ;ELF  = ZERO       !!SURFACE ELEVATION
   ALLOCATE(ET(0:NTNode))            ;ET   = ZERO       !!SURFACE ELEVATION
   ALLOCATE(EGF(0:NTNode))           ;EGF  = ZERO       !!SURFACE ELEVATION
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK(0:NTNode))          ;ELRK = ZERO       !!SURFACE ELEVATION
#  endif
               MEMCNT = NTNode*8*NDB + MEMCNT

!---------------surface/bottom/edge boundary conditions-----------------------------!

   ALLOCATE(CBC(0:MTElem))           ;CBC     = ZERO     !!BOTTOM FRICTION     
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(WUSURF2(0:MTElem))       ;WUSURF2 = ZERO     !!SURFACE FRICTION FOR EXT
   ALLOCATE(WVSURF2(0:MTElem))       ;WVSURF2 = ZERO     !!SURFACE FRICTION FOR EXT
# endif
   ALLOCATE(WUBOT(0:MTElem))         ;WUBOT   = ZERO     !!BOTTOM FRICTION
   ALLOCATE(WVBOT(0:MTElem))         ;WVBOT   = ZERO     !!BOTTOM FRICTION
   ALLOCATE(WUSURF(0:MTElem))        ;WUSURF  = ZERO     !!SURFACE FRICTION FOR INT
   ALLOCATE(WVSURF(0:MTElem))        ;WVSURF  = ZERO     !!SURFACE FRICTION FOR INT
   ALLOCATE(UUWIND(0:MTElem))        ;UUWIND  = ZERO     !!SURFACE X-WIND 
   ALLOCATE(VVWIND(0:MTElem))        ;VVWIND  = ZERO     !!SURFACE Y-WIND 
   ALLOCATE(SWRAD(0:NTNode))         ;SWRAD   = ZERO     !!SURFACE INCIDENT RADIATION
   ALLOCATE(WTSURF(0:NTNode))        ;WTSURF  = ZERO 
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(QPREC2(0:NTNode))        ;QPREC2  = ZERO     !!SURFACE PRECIPITATION FOR EXT
   ALLOCATE(QEVAP2(0:NTNode))        ;QEVAP2  = ZERO     !!SURFACE EVAPORATION FOR EXT
#  endif
   ALLOCATE(QPREC3(0:NTNode))        ;QPREC3  = ZERO     !!SURFACE PRECIPITATION FOR INT
   ALLOCATE(QEVAP3(0:NTNode))        ;QEVAP3  = ZERO     !!SURFACE EVAPORATION FOR INT

!#  if defined (ICE)
   ALLOCATE(T_AIR(0:NTNode))         ;T_AIR   = ZERO     !!BULK AIR TEMPERATURE
   ALLOCATE(RH_AIR(0:NTNode))        ;RH_AIR  = ZERO     !!RELATIVE HUMIDITY
   ALLOCATE(QA_AIR(0:NTNode))        ;QA_AIR  = ZERO     !!SPECIFIC HUMIDITY
   ALLOCATE(PA_AIR(0:NTNode))        ;PA_AIR  = ZERO     !!SURFACE PRESSURE
   ALLOCATE(DLW_AIR(0:NTNode))       ;DLW_AIR = ZERO     !!DOWNWARD LONGWAVE RADIATION
   ALLOCATE(DSW_AIR(0:NTNode))       ;DSW_AIR = ZERO     !!DOWNWARD SHORTWAVE RADIATION
   ALLOCATE(CLOUD(0:NTNode))         ;CLOUD   = ZERO     !!Total Cloud Cover
   ALLOCATE(WSSURF(0:NTNode))        ;WSSURF  = ZERO     !!sea surface salinity    flux
   ALLOCATE(WHSURF(0:NTNode))        ;WHSURF  = ZERO     !!sea surface heat        flux
   ALLOCATE(FWtSURF(0:NTNode))       ;FWtSURF = ZERO     !!sea surface fresh water flux
!#  endif

!!By T.W. for initializing kelp arrays, April 2013
#  if defined (KELP)

   ALLOCATE(n_of_kelp(0:MTElem,kbm1)); n_of_kelp = 0     !!local array for # of kelps
   ALLOCATE(drag_coef(0:MTElem,kbm1)); drag_coef = ZERO  !!energy extraction coef.
   ALLOCATE(kelp_area(0:MTElem,kbm1)); kelp_area = ZERO  !!local array for turbine facing angle

   ALLOCATE(EMS_X(0:MTElem,kb)); EMS_X = ZERO
   ALLOCATE(EMS_Y(0:MTElem,kb)); EMS_Y = ZERO

#  endif
!! finished by T.W. for kelp module

!------------------------ GROUNDWATER T,S FLUX AT CURRENT TIM-------------------
   ALLOCATE(BFWDIS(IBFW))        ;BFWDIS  = ZERO     !!GROUNDWATER FLUX AT CURRENT TIME
   ALLOCATE(BFWTDIS(IBFW))       ;BFWTDIS = ZERO     !!GROUNDWATER TEMPERATURE FLUX AT CURRENT TIME
   ALLOCATE(BFWSDIS(IBFW))       ;BFWSDIS = ZERO     !!GROUNDWATER SALINITY FLUX AT CURRENT TIME
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(BFWDIS2(IBFW))       ;BFWDIS2  = ZERO     !!GROUNDWATER FLUX FOR EXT
!   ALLOCATE(BFWTDIS2(IBFW))      ;BFWTDIS2 = ZERO     !!GROUNDWATER TEMPERATURE FLUX AT CURRENT TIME FOR EXT
!   ALLOCATE(BFWSDIS2(IBFW))      ;BFWSDIS2 = ZERO     !!GROUNDWATER SALINITY FLUX AT CURRENT TIME FOR EXT
#  endif
   ALLOCATE(BFWDIS3(IBFW))       ;BFWDIS3  = ZERO     !!GROUNDWATER FLUX FOR INT
   ALLOCATE(BFWTDIS3(IBFW))      ;BFWTDIS3 = ZERO     !!GROUNDWATER TEMPERATURE FLUX AT CURRENT TIME FOR INT
   ALLOCATE(BFWSDIS3(IBFW))      ;BFWSDIS3 = ZERO     !!GROUNDWATER SALINITY FLUX AT CURRENT TIME FOR INT

!------------------------------------------------------------------------------
   ALLOCATE(QDIS(NUMQBC))        ;QDIS    = ZERO     !!FRESH WATER FLUX AT CURRENT TIME
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(QDIS2(NUMQBC))       ;QDIS2   = ZERO     !!FRESH WATER (EXT MODE NOT USED )  
#  endif
   ALLOCATE(TDIS(NUMQBC))        ;TDIS    = ZERO     !!FRESH WATER TEMP AT CURRENT TIME
   ALLOCATE(SDIS(NUMQBC))        ;SDIS    = ZERO     !!FRESH WATER SLNT AT CURRENT TIME
   ALLOCATE(QAREA(NUMQBC))       ;QAREA   = ZERO     !!AREA OF FLUX                    
   ALLOCATE(RDISQ(NUMQBC,2))     ;RDISQ   = ZERO     !!AREA OF FLUX                    
   ALLOCATE(ANGLEQ(NUMQBC))      ;ANGLEQ  = ZERO     !!FLUX VELOCITY                   
   ALLOCATE(VLCTYQ(NUMQBC))      ;VLCTYQ  = ZERO     !!FLUX ANGLE 
          MEMCNT = MTElem*9*NDB + NTNode*4*NDB + NUMQBC*9*NDB + MEMCNT

!-----------------------2-d flow fluxes---------------------------------------------!
   ALLOCATE(PSTX(0:MTElem))          ;PSTX  = ZERO       !!EXT MODE BAROTROPIC TERMS
   ALLOCATE(PSTY(0:MTElem))          ;PSTY  = ZERO       !!EXT MODE BAROTROPIC TERMS
   ALLOCATE(ADVUA(0:MTElem))         ;ADVUA = ZERO 
   ALLOCATE(ADVVA(0:MTElem))         ;ADVVA = ZERO
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ADX2D(0:MTElem))         ;ADX2D = ZERO
   ALLOCATE(ADY2D(0:MTElem))         ;ADY2D = ZERO
   ALLOCATE(DRX2D(0:MTElem))         ;DRX2D = ZERO
   ALLOCATE(DRY2D(0:MTElem))         ;DRY2D = ZERO
# endif
   ALLOCATE(TPS(0:MTElem))           ;TPS   = ZERO      !!WORKING ARRAY             
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ADVX(0:MTElem,KB))       ;ADVX  = ZERO 
   ALLOCATE(ADVY(0:MTElem,KB))       ;ADVY  = ZERO 
# endif
               MEMCNT = MTElem*9*NDB + MTElem*KB*2*NDB + MEMCNT


!---------------- internal mode   arrays-(element based)----------------------------!

   ALLOCATE(U(0:MTElem,KB))          ;U       = ZERO   !!X-VELOCITY
   ALLOCATE(V(0:MTElem,KB))          ;V       = ZERO   !!Y-VELOCITY
   ALLOCATE(UBETA(0:MTElem,KB))      ;UBETA   = ZERO   !X-VELOCITY temp time step without Coriolis
   ALLOCATE(VBETA(0:MTElem,KB))      ;VBETA   = ZERO   !Y-VELOCITY temp time step without Coriolis
   ALLOCATE(UBETA2D(0:MTElem))       ;UBETA2D = ZERO
   ALLOCATE(VBETA2D(0:MTElem))       ;VBETA2D = ZERO

   ALLOCATE(W(0:MTElem,KB))          ;W     = ZERO   !!VERTICAL VELOCITY IN SIGMA SYSTEM
   ALLOCATE(WW(0:MTElem,KB))         ;WW    = ZERO   !!Z-VELOCITY
   ALLOCATE(UF(0:MTElem,KB))         ;UF    = ZERO   !!X-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(VF(0:MTElem,KB))         ;VF    = ZERO   !!Y-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(WT(0:MTElem,KB))         ;WT    = ZERO   !!Z-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(RHO(0:MTElem,KB))        ;RHO   = ZERO   !!DENSITY AT ELEMENTS
   ALLOCATE(RMEAN(0:MTElem,KB))      ;RMEAN = ZERO   !!MEAN INITIAL DENSITY AT ELEMENTS
   ALLOCATE(T(0:MTElem,KB))          ;T     = ZERO   !!TEMPERATURE AT ELEMENTS
   ALLOCATE(TMEAN(0:MTElem,KB))      ;TMEAN = ZERO   !!MEAN INITIAL TEMPERATURE AT ELEMENTS
   ALLOCATE(S(0:MTElem,KB))          ;S     = ZERO   !!SALINITY AT ELEMENTS
   ALLOCATE(SMEAN(0:MTElem,KB))      ;SMEAN = ZERO   !!MEAN INITIAL SALINITY AT ELEMENTS
               MEMCNT = MTElem*KB*13*NDB + MEMCNT

!-----------------------3d variable arrays-(node based)-----------------------------!

   ALLOCATE(T1(0:NTNode,KB))         ;T1     = ZERO  !!TEMPERATURE AT NODES               
   ALLOCATE(S1(0:NTNode,KB))         ;S1     = ZERO  !!SALINITY AT NODES               
   ALLOCATE(RHO1(0:NTNode,KB))       ;RHO1   = ZERO  !!DENSITY AT NODES               
   ALLOCATE(TF1(0:NTNode,KB))        ;TF1    = ZERO  !!TEMPERATURE FROM PREVIOUS TIME
   ALLOCATE(SF1(0:NTNode,KB))        ;SF1    = ZERO  !!SALINITY FROM PREVIOUS TIME 
   ALLOCATE(TMEAN1(0:NTNode,KB))     ;TMEAN1 = ZERO  !!MEAN INITIAL TEMP
   ALLOCATE(SMEAN1(0:NTNode,KB))     ;SMEAN1 = ZERO  !!MEAN INITIAL SALINITY 
   ALLOCATE(RMEAN1(0:NTNode,KB))     ;RMEAN1 = ZERO  !!MEAN INITIAL DENSITY 
   ALLOCATE(WTS(0:NTNode,KB))        ;WTS    = ZERO  !!VERTICAL VELOCITY IN SIGMA SYSTEM
   ALLOCATE(WTTS(0:NTNode,KB))       ;WTTS   = ZERO  !!WTS FROM PREVIOUS TIMESTEP        
   ALLOCATE(Q2(0:NTNode,KB))         ;Q2     = ZERO   !!TURBULENT KINETIC ENERGY AT NODES
   ALLOCATE(Q2L(0:NTNode,KB))        ;Q2L    = ZERO   !!TURBULENT KE*LENGTH AT NODES
#  if defined (GOTM)
   ALLOCATE(TKE(0:NTNode,KB))        ;TKE   = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TKEF(0:NTNode,KB))       ;TKEF  = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TEPS(0:NTNode,KB))       ;TEPS  = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TEPSF(0:NTNode,KB))      ;TEPSF = ZERO   !!TURBULENT KE*LENGTH AT ELEMENTS
#  endif
   ALLOCATE(L(0:NTNode,KB))          ;L     = ZERO   !!TURBULENT LENGTH SCALE AT ELEMENTS
   ALLOCATE(KM(0:NTNode,KB))         ;KM    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(KH(0:NTNode,KB))         ;KH    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(KQ(0:NTNode,KB))         ;KQ    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(AAM(0:NTNode,KB))        ;AAM   = ZERO   !!??
   ALLOCATE(KM1(0:MTElem,KB))        ;KM1   = ZERO   !!TURBULENT QUANTITY AT ELEMENTS
               MEMCNT = NTNode*KB*18*NDB + MTElem*KB*NDB + MEMCNT
  
!---------------------------internal mode fluxes------------------------------------!

   ALLOCATE(DRHOX(0:MTElem,KB))      ;DRHOX = ZERO 
   ALLOCATE(DRHOY(0:MTElem,KB))      ;DRHOY = ZERO 
   ALLOCATE(Q2F(0:NTNode,KB))        ;Q2F   = ZERO 
   ALLOCATE(Q2LF(0:NTNode,KB))       ;Q2LF  = ZERO
               MEMCNT = MTElem*KB*2*NDB + NTNode*KB*2*NDB + MEMCNT

!------------shape coefficient arrays and control volume metrics--------------------!

   ALLOCATE(A1U(0:MTElem,4))         ;A1U   = ZERO
   ALLOCATE(A2U(0:MTElem,4))         ;A2U   = ZERO 
   ALLOCATE(AWX(0:MTElem,3))         ;AWX   = ZERO 
   ALLOCATE(AWY(0:MTElem,3))         ;AWY   = ZERO 
   ALLOCATE(AW0(0:MTElem,3))         ;AW0   = ZERO 
   ALLOCATE(ALPHA(0:MTElem))         ;ALPHA = ZERO
               MEMCNT = MTElem*4*2 + MTElem*3*3 + MTElem + MEMCNT

!-----salinity and temperature bottom diffusion condition/bottom depth gradients----!

   ALLOCATE(PHPN(NTNode))            ;PHPN      = ZERO 
   ALLOCATE(PFPXB(NTNode))           ;PFPXB     = ZERO
   ALLOCATE(PFPYB(NTNode))           ;PFPYB     = ZERO
   ALLOCATE(SITA_GD(NTNode))         ;SITA_GD   = ZERO 
   ALLOCATE(AH_BOTTOM(NTNode))       ;AH_BOTTOM = ZERO 
               MEMCNT = NTNode*5*NDB + MEMCNT

   ALLOCATE(VISCOFH(0:NTNode,KB))    ;VISCOFH = ZERO
               MEMCNT = NTNode*KB*NDB + MEMCNT

!---------------Coordinates of Center Pionts around the Nodes-----------------------!
   ALLOCATE(XCA(NNode))             ;XCA        = ZERO
   ALLOCATE(YCA(NNode))             ;YCA        = ZERO
   ALLOCATE(XCB(NNode))             ;XCB        = ZERO
   ALLOCATE(YCB(NNode))             ;YCB        = ZERO
   ALLOCATE(XCC(NNode,10))          ;XCC        = ZERO !ASSUMING THE MAXIUM OF NEIGHBORING ELEMENTS IS NOT MORE THAN 10
   ALLOCATE(YCC(NNode,10))          ;YCC        = ZERO !ASSUMING THE MAXIUM OF NEIGHBORING ELEMENTS IS NOT MORE THAN 10
   ALLOCATE(XCD(NNode,10))          ;XCD        = ZERO !ASSUMING THE MAXIUM OF NEIGHBORING ELEMENTS IS NOT MORE THAN 10
   ALLOCATE(YCD(NNode,10))          ;YCD        = ZERO !ASSUMING THE MAXIUM OF NEIGHBORING ELEMENTS IS NOT MORE THAN 10
   ALLOCATE(XCE(NNode))             ;XCE        = ZERO
   ALLOCATE(YCE(NNode))             ;YCE        = ZERO
   ALLOCATE(XCF(NNode))             ;XCF        = ZERO
   ALLOCATE(YCF(NNode))             ;YCF        = ZERO
   ALLOCATE(VAL_COS_VY(NNode))      ;VAL_COS_VY = ZERO

!-----special initialization which probably do nothing------------------------------!

   DT1(0) = 100.0_SP

!---------------report approximate memory usage-------------------------------------!

   MEMTOT = MEMCNT*4
#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(MEMCNT,MEMTOT,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
#  endif
   IF(MSR)WRITE(IPT,*)'!  # MBYTES TOTAL MEM    :',MEMTOT/1E+6
   IF(MSR .AND. .NOT.SERIAL )WRITE(IPT,*)'!  # AVERAGE MBYTES/PROC :',MEMTOT/(1E+6*NPROCS)
   

   RETURN
   END SUBROUTINE ALLOC_VARS
!==============================================================================|
