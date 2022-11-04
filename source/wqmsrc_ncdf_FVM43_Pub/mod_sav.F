MODULE MOD_SAV

!************************************************************************
!**                                                                    **
!**       Submerged Aquatic Vegetation Subroutines for FVCOM-ICM       **
!**               Wen Long, Seattle, WA, Nov 2014					   **
!**                                                                    **
!************************************************************************

	USE MOD_PREC, ONLY: SP
	

     CHARACTER(LEN=8) SPNAM,    &  !species name of the SAV
                      EPINAM       !name of the epiphyte

     INTEGER :: NSAVM            !global maximum number of SAV species in all of tracer control elements (TCE's) - polygon
                                 !surrounding a given node
     
     REAL(SP) :: KHSTOX,         &   !Sulfide toxicity coefficient in calculating HS^- limitation on SAV growth (m^3/mole-S2)
             KESAV,          &   !Effect on light attenuation by SAV leaf and stem, KESAV (gC/m^3)  
             KEEPI,          &   !Shading effect on light on SAV shoots (leaf) due to Epiphytes ((m^2-leaf-surface/gDW-epiphyte)
             KHEP,           &   !half saturation desnity of epiphyte (gC-Epiphyte/gC-leaf) for epiphyte growth on SAV leaf
             NLIMEPI             !Never used ???                      
             
     REAL(SP) :: KHNLEAF,         &  !Half saturation constant for N uptake by leaf (shoots) (gN/m^3)
             KHNROOT,         &  !Half saturation constant for N uptake by root (gN/m^3)
             KHPLEAF,         &  !Half saturation constant for P uptake by leaf (gP/m^3)
             KHPROOT,         &  !Half saturation constant for P uptake by root (gP/m^3)
             KHNEPI,          &  !Half saturation constant for N uptake by epiphytes (gN/m^3)
             KHPEPI,          &  !Half saturation constant for P uptake by epiphytes (gN/m^3)
             KHNPSAV,         &  !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for SAV
             KHNPEPI             !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for epiphytes

     REAL(SP) :: PMSAV,       &   !maximum production rate of SAV, gC/gDW/day  where DW is dry weight
             PRSPSAVB,    &   !fraction of leaf production rate that is counted as photochemical respiration at salt << saltmax
                              !at salt >> saltmax, PRSPSAV approaches 1, and there will be no production of SAV leaf hence
                              !no photosynthesis for SAV growth
             BMSAV,       &   !Basal metabolism rate of SAV (1/day)
             BMTBRREF,    &   !Reference basal metabolism rate of tube [1/day] BMTUBE=BMTBRREF*f(T), i.e. function of temperature
             FDOSR,       &   !Fraction of DO production discounted due to route of Leaf production to root
             SALMAX,      &   !Salinity at which PRSPSAV is approaching 1 (maximum) from PRSPSAVB, [psu]
             ALPHSAV,     &   !Initial slope of P-I curve for SAV [ (gC/gDW)/(E/m^2) ]
             ALAC,        &   !m^2-leaf-area/gC-leaf (leaf area to leaf carbon ratio) 
                              !Note in 2002 version, it is the inverse ACLA: gC-leaf/m^2-leaf-area
             SLSAV,       &   !???? Never used !Sloughing rate of SAV ?  1/day
             ACDWSAV,     &   !Carbon to Dry Weight ratio of SAV (gC/gDW)
             ANDWSAV,     &   !Nitrogen to Dry Weight ratio of SAV (gN/gDW)
             APDWSAV,     &   !Phosphorus to Dry Weight ratio of Sav (gP/gDW)
             ANCSAV,      &   !Nitrogen to Carbon ratio of SAV (gN/gC)             
             APCSAV,      &   !P to C ratio of SAV (gP/gC)
             HCAN,        &   !Canopy height (m), HCAN = ACAN + BCAN *( leaf+stem )
             ACAN,        &   !basic canopy height (m)
             BCAN             !coefficient for calculating leaf (gC/m^2) and stem (gC/m^2) 
                              !contribution to canopy height,  [m^3/gC]

   !Distribution of C, N,P on SAV mortality      
     REAL(SP) :: FNISAV,   &      !Fraction of NH4 Nitrogen generated due to SAV (leaf, stem etc) respiration (mortality)
             FNLDSAV,  &      !Fraction of Nitrogen routed to LDON on SAV mortality
             FNRDSAV,  &      !Fraction of Nitrogen routed to RDON on SAV mortality
             FNLPSAV,  &      !Fraction of Nitrogen routed to LPON on SAV mortality
             FNRPSAV,  &      !Fraction of Nitrogen routed to RPON on SAV mortality
             FPISAV,   &      !Fraction of Phosphorus routed to PO4 in water column due to SAV (leaf,stem) respiration
             FPLDSAV,  &      !Fraction of P routed to LDOP in water column due to SAV mortality
             FPRDSAV,  &      !Fraction of P routed to RDOP in water column due to SAV mortality
             FPLPSAV,  &      !Fraction of P routed to LPOP in water column due to SAV mortality
             FPRPSAV,  &      !Fraction of P routed to RPOP in water column due to SAV mortality
             FCLDSAV,  &      !Fraction of Carbon routed to LDOC in water column due to SAV mortality
             FCRDSAV,  &      !Fraction of Carbon routed to RDOC in water column due to SAV mortality
             FCLPSAV,  &      !Fraction of Carbon routed to LPOC in water column due to SAV mortality
             FCRPSAV,  &      !Fraction of Carbon routed to RPOC in water column due to SAV mortality
             FDOSAV           !Fraction of Carbon routed to CO2 (oxygen consumption) in water column due to SAV mortality
                              !Note FDOSAV + FCLDSAV+ FCRDSAV + FCLPSAV + FCRPSAV = 1
                              !Wen Long: we really need to rename FDOSAV to FCO2SAV and make it equivalent to FDOSAV!

     REAL(SP) :: WSSSAV,   &      !Increase of fixed solids (ISS) settling by SAV leaf and stem, (m/day)/(gC/m^2)==>m^3/gC/day
             WSLSAV,   &      !Increase of LPOM settling by SAV leaf and stem, m^3/gC/day
             WSRSAV,   &      !Increase of RPOM settling by SAV leaf and stem, m^3/gC/day
             WS1SAV,   &      !Increase of Algae group 1 settling by SAV leaf and stem, m^3/gC/day
             WS2SAV,   &      !                        2
             WS3SAV,   &      !                        3
             WSUSAV           !Increase of particulate biogenic silica settling rate by SAV leaf and stem, m^3/gC/day
                     
  INTEGER  :: NDDI,     & !number of degree days for entries of look up table of FPLEAF etc for Senescence of leaf etc
              NSAVCELL, & !total number of grid cells (ICM boxes) that have SAV --> need to be changed to nodes of FVCOM
              NSAVSPCM    !maximum number of species in a cell
              
     REAL(SP) :: PMEPI,        & !maximum photosynthesis production rate (gC/gCHLA/day) of epiphytes
             BMEPI,        & !basal metabolism loss rate of epiphyte (1/day)
             PREPI,        & !predation rate on epiphytes (gC-leaf/gC-epiphyte/day)
             PRSPEPI         !Photochemical respiration rate of epiphyte --> BMEP

     REAL(SP) :: ALPHEPI,      & !slope of PI curve for epiphyte photothenesis (gC/gDW)/(E/m^2)
             CCHLEPI,      & !Carbon to Chlorophyll ratio of epiphyte (gC/gCHLA)
             ANCEPI,       & !N to C ratio of epiphytes (gN/gC)
             APCEPI,       & !P to C ratio of epiphytes (gP/gC)
             ADWCEPI         !ratio of total epiphyte dry weight to viable carbon (gDW/gC)
             
     REAL(SP) :: FCLDEPI,     &   !Fraction of C routed to LDOC due to basal metabolism of epiphyte
             FCRDEPI,     &   !Fraction of C routed to RDOC due to basal metabolism of epiphyte
                              !WLong: missing FCLPEPI and FCRPEPI to routt to LPOC and RPOC ??? (seems not coded)
                              !FCLPEPI and FCRPEPI should be added (Wen Long)
                              !also missing FDOEPI compared to FDOSAV for SAV
                              !code assumes that no carbon will be routed to LPOP and RPOP due to metabolism
                              !and also FDOEPI=1-FCLDEPI-FCRDEPI, and FDOEPI is only implicitly used in calculating 
                              !DOEPIW
             
             FNIEPI,      &   !Fraction of N routed to NH4  due to epiphyte basal metabolism
             FNLDEPI,     &   !Fraction of N routed to LDON due to epiphyte basal metabolism
             FNRDEPI,     &   !Fraction of N routed to RDON due to epiphyte basal metabolism
             FNLPEPI,     &   !Fraction of N routed to LPON due to epiphyte basal metabolism
             FNRPEPI,     &   !Fraction of N routed to RPON due to epiphyte basal metabolism
             FPIEPI,      &   !Fraction of P routed to  PO4 due to basal metabolism of epiphyte
             FPLDEPI,     &   !Fraction of P routed to LDOP due to basal metabolism of epiphyte
             FPRDEPI,     &   !Fraction of P routed to RPOP due to basal metabolism of epiphyte
             FPLPEPI,     &   !Fraction of P routed to LPOP due to basal metabolism of epiphyte
             FPRPEPI          !Fraction of P routed to RPOP due to basal metabolism of epiphyte             

     !Distribution of C, N, P upon epiphyte mortality (predation)
     REAL(SP) :: FCLDPEP,     &   !Fraction of C routed to LDOC due to predation and sloughing on epiphytes
             FCRDPEP,     &   !Fraction of C routed to RDOC due to predation and sloughing on epiphytes
             FCLPPEP,     &   !Fraction of C routed to LPOC due to predation and sloughing on epiphytes
             FCRPPEP,     &   !Fraction of C routed to RPOC due to predation and sloughing on epiphytes
             FNIPEP,      &   !Fraction of Nitrogen routed to NH4  due to predation on epiphytes and Slough
             FNLDPEP,     &   !Fraction of Nitrogen routed to LDON due to predation on epiphytes and Slough
             FNRDPEP,     &   !Fraction of Nitrogen routed to RDON due to predation on epiphytes and Slough
             FNLPPEP,     &   !Fraction of Nitrogen routed to LPON due to predation on epiphytes and Slough
             FNRPPEP,     &   !Fraction of Nitrogen routed to RPON due to predation on epiphytes and Slough
             FPIPEP,      &   !Fraction of P routed to PO4  due to predation on epiphytes and slough
             FPLDPEP,     &   !Fraction of P routed to LDOP due to predation on epiphytes and slough
             FPRDPEP,     &   !Fraction of P routed to RDOP due to predation on epiphytes and slough
             FPLPPEP,     &   !Fraction of P routed to LPOP due to predation on epiphytes and slough
             FPRPPEP          !Fraction of P routed to RPOP due to predation on epiphytes and slough            

     REAL(SP) :: FPSR      !Never used!!!! This is supposed to be the fraction of leaf production routed to root
                       !(replaced by FPROOT)

     !Wen Long: these fractions should change due to aging (senescence)
     REAL(SP),DIMENSION(0:1000) :: FPLEAF, & !fraction of production assigned to leaf growth
                               FPROOT, & !fraction of production assigned to root growth
                               FPSTEM, & !fraction of production assigned to stem growth
                               FPTUBER,& !fraction of production assigned to tuber growth
                                         !WLong: FPLEAF +FPROOT+FPSTEM+FPTUBER = 1
                                         !all growth is based on photosynthesis by leaf
                                         !(so stem and tuber have no photosynthesis?? I thought they look green )
                               TRTBRLF   !transfer rate from tuber to leaf (gC-Leaf/day/gC-Tuber)
                               
     !Look up table of temperature control in a year on growth and mortality rates of SAV and epiphytes
     REAL(SP),DIMENSION(-50:400) ::  FTPSAV, & !temperature control on SAV growth
                                 FTRSAV, & !temperature control on SAV resipration 
                                 FTPEP,  & !temperature control on epiphyte growth rate
                                 FTREP,  & !temperature control on epiphyte resipration (basal metabolism)
                                 FTPREP    !temperature control on predation on epiphyte
                               
     !fluxes to sediments by SAV mortality
     REAL(SP),DIMENSION(3) :: FRPOCSAV, FRPONSAV, FRPOPSAV    !3G partitioning of POC, PON, POP sources to sediments from SAV (non-dimensional)
  
     !State variables of SAV model 
     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: LEAF,   & !Leaf biomass (gC/m^2) 
                                        ROOT,   & !Root biomass (gC/m^2)
                                        STEM,   & !Stem biomass (gC/m^2)
                                        TUBER,  & !Tuber biomass (gC/m^2)
                                        EP        !Epiphytes biomass (gC-epiphytes/gm^2-leaf)
                                                  !Note in 2002 version of ICM, EP's unit was (gC-epiphytes/gC-leaf)
                                 
     REAL(SP),ALLOCATABLE,DIMENSION(:)    :: SAVAREA    !area of SAV coverage for each surface cell (m^2) --> needs to be changed for FVCOM grid                                 
     INTEGER,ALLOCATABLE,DIMENSION(:) :: SAVCELL    !grid cell number for all cells that have SAV in it
                                                    !size 1xNSAVCELL --> need to change to element number or node number!!
     REAL(SP),ALLOCATABLE :: NSAVSPC(:)                 !number of SAV species in each grid cell
                                                    !size 1xNSAVCELL

     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: SAVDPH,  &  !Depth of SAV (mean depth of the plant) for all species in each grid cell
                                                    !size, NSAVCELL x NSAVSPC(B) where B is cell number, B=1,...,NSAVCELL
                                        SAVFRAC     !coverage for each SAV cell and each vertical deph increment
                                                    !i.e. coverage is function of grid cell and SAV species

!Growth and mortality rates of SAV
     REAL(SP),ALLOCATABLE,DIMENSION(:,:) ::  PLEAF    !Leaf growth rate (1/day), P (SAV shoots growth) in C-M 2001 paper's SH equation
     REAL(SP),ALLOCATABLE,DIMENSION(:) ::  BMLEAF, &  !respiration rate of leaf of SAV (1/day) with temperature control 
                                       BMSTEM, &  !resipration rate of stem of SAV (1/day) with temperature control
                                       BMROOT, &  !respiration rate of root of SAV (1/day) with temperature control
                                       BMTUBER    !respiration rate of tuber of SAV (1/day) with temperature control                                      
     
    !Growth and mortality rates of SAV                                      
     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: PEP,   &  !epiphytes growth rate (1/day)  = PMEPI*light_limitation*deinsity_limitation/CCHLEPI
                                        BMEP,  &  !epiphytes basal metabolism rate (1/day) (also including photochemical respiration)
                                        PREP      !predation rate on epiphytes (1/day) = PREPI*EP
                                                  !where PREPI has unit gC-shoot/gC-epiphyte/day

     REAL(SP), ALLOCATABLE :: SLSH(:)                 !epiphyte sloughing rate (1/day)
     
    !Nutrient limitations on SAV and epiphytes growth
    REAL(SP),ALLOCATABLE,DIMENSION(:) ::   NLSAV,    &  !Nitrogen limiation on SAV growth (non-dimensional)
                                       PLSAV,    &  !Phosphorus limitation on SAV growth (non-dimensional)
                                       FNSEDSAV, &  !fraction of nitrogen uptake by SAV due to nitrogen in sediment pore water
                                                    !Wen Long: I think we should further split into NH4 and NO3 
                                                    !          uptake by SAV from sediments
                                       FPSEDSAV, &  !fraction of phosphorus uptake by SAV due to P in sediment pore water
                                       NLEPI,    &  !Nitrogen limitation on epiphytes growth (non-dimensional)
                                       PLEPI,    &  !Phosphorus limitation on epiphytes growth (non-dimensional)
                                       FHS          !fresh water SAV species growth limitation due to sulfide HS^- 
                                                    !(in sediment anerobic layer) toxicity (dimensionless)

     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: EPATN,  & !Light attenuation by epiphyte
                                                  !e^(-Kep*Acla*Adwcep*EP) in Cerco-Moore 2001 paper for Ish eqn
                                                  !e^(-Kep*Adwcep*EP) in new version as EP changed to gC-Epi/m^2-leaf unit
                                        SAVATN, & !Light attenuation by SAV (shoots)
                                                  ![1-e^(-Ksh*SH+KESS*HCAN)]/(Ksh*SH+KESS*HCAN) in C-M 2001 Iwc eqn
                                        WATATN, & !water column attenuation of light(dimensionless) at canopy top
                                                  != exp(-KESS*ZTOCNPY) where ZTOCNPY is distance from surface
                                                  !                           to canopy top
                                        FISH,   & !light limitation function for LEAF (shoots) f(I)=Ish/sqrt(Ish^2+Ik^2)
                                        FIEP,   & !light limitation function for epiphytes, f(I)=Iep/sqrt(Iep^2+Ik^2)
                                        NPPSAV, & !net leaf production per day (growth-basal metabolism) (gC/m^2/day)
                                        NPPEPI    !net epiphyte production per day (growth - basal metabolism) (gC-epi/gC-leaf/day)

     !fluxes to water column
     REAL(SP),ALLOCATABLE,DIMENSION(:) :: DOSAVW, &  !DO source (flux) due to SAV (source-sink) (gO2/m^2/day ) 
                                                 !DOSAVW --> DTDO (gO2/m^3/sec)--> DTC (gO2/m^3/sec)
                                                 !-->DTM(gO2/m^3/sec)-->DTM*DLT-->C1(:,:,27)(gO2/m^3) !27'th constituents is DO
                                    LDOCSAVW, &  !LDOC flux due to SAV (gC/m^2/day)
!WL                                              !***ALSO used by TVDS.F ****
                                                 !LDOCSAVW --> DTDOC (gC/m^3/sec) -> DTC --> DTM --> DTM*DLT --> C1(:,:,9)
                                    RDOCSAVW, &  !RDOC flux due to SAV (gC/m^2/day)                                 
                                    LPOCSAVW, &  !LPOC flux due to SAV (gC/m^2/day)
!WL                                              !***ALSO used by TVDS.F ****
                                    RPOCSAVW, &  !RPOC flux due to SAV (gC/m^2/day)
!WL                                              !***ALSO used by TVDS.F ****
                                    NH4SAVW,  &  !NH4 flux to water column (respiration source - uptake sink by SAV growth) (gN/m^2/day)
                                                 !NH4SAVW (gN/m^2/day) --> DTNH4  (gN/m^3/sec) --> DTC --> DTM --> C1(:,:,13)
                                    NO3SAVW,  &  !NO3 flux to water column (due to uptake by SAV, positive increasing NO3 in water) (gN/m^2/day)
                                    LDONSAVW, &  !LDON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf, basal metabolism of stem) (gN/m^2/day)
                                    RDONSAVW, &  !RDON flux                                     
                                    LPONSAVW, &  !LPON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem (gN/m^2/day)
                                    RPONSAVW, &  !RPON flux 
                                    PO4SAVW,  &  !PO4 flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem and uptake by leaf) (gP/m^2/day)
                                    LDOPSAVW, &  !LDOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day)
                                    RDOPSAVW, &  !RDOP flux (gP/m^2/day)
                                    LPOPSAVW, &  !LPOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day) 
                                    RPOPSAVW, &  !RPOP flux (gP/m^2/day)
                                      DOEPIW, &  !DO   flux (positive into water column) due to epiphytes (gO2/m^2/day)
                                    LDOCEPIW, &  !LDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
                                    RDOCEPIW, &  !RDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
                                    LPOCEPIW, &  !LPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
                                    RPOCEPIW, &  !RPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
                                     NH4EPIW, &  !NH4  flux (positive into water column) due to epiphytes (gN/m^2/day) (basal meabolism + predation + sloughing - uptake by photosynthesis)
                                     NO3EPIW, &  !NO3  flux
                                    LDONEPIW, &  !LDON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
                                    RDONEPIW, &  !RDON flux (positive into water column) (gN/m^2/day)
                                    LPONEPIW, &  !LPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
                                    RPONEPIW, &  !RPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
                                    PO4EPIW,  &  !PO4  flux (positive into water column) due to epiphytes mortality (basal metabolism, salt toxicity, sloughing) and uptake (sink) (growth) (gP/m^2/day)                                    
                                    LDOPEPIW, &  !LDOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
                                    RDOPEPIW, &  !RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
                                    LPOPEPIW, &  !LPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
                                    RPOPEPIW     !RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
  !fluxes to sediments                                                                         
  REAL(SP), ALLOCATABLE, DIMENSION(:):: SEDPOCSAV, &  !POC flux to sediments (positive into sediments) due to SAV (gC/m^2/day)
                                                  !SEDPOCSAV (gC/m^2/day) --> JPOC (mgC/m^2/day) -> POC1, POC2, POC3 (source of particulate organic carbon in sediments)
                                    SEDPONSAV, &  !PON flux to sediments (positive into sediments) due to SAV (gN/m^2/day)
                                                  !SEDPONSAV (gN/m^2/day) --> JPON (mgN/m^2/day) -> PON1, PON2, PON3 (source of particulate organic nitrogen in sediments)
                                    SEDPOPSAV, &  !POP flux to sedimetns (positive into sediments) due to SAV (gP/m^2/day)
                                                  !SEDPOPSAV (gP/m^2/day) --> JPON (mgP/m^2/day) -> POP1, POP2, POP3 (source of particulate organic phosphorus in sediments)

                                    SEDDOSAV,  &  !DO flux to water column (positive into water) due to SAV growth in sediments (source) and mortality of SAV in sediments (tuber and root basal metabolism)(sink)
                                                  !(gO2/m^2/day) (WLong: does this have to be to water ?? I thought the sinks should be in sediments directly for root metabolism)
									
												  !WLong : we might need to split into DO production in water column
												  !due to leaf growht 
												  ! and Sediment Oxygen Demand due to root respiration
									
                                    SEDNH4SAV, &  !NH4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
                                                  !(gN/m^2/day), SEDNH4SAV --> NH4T2TM1 (mg/N/m^3) decerase in sediments 
                                    SEDNO3SAV, &  !NO3 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
                                                  !(gN/m^2/day), SEDNO3SAV --> NO3T2TM1 (mg/N/m^3) decerase in sediments 
                                    SEDPO4SAV     !PO4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
                                    
   !misc
   REAL(SP), ALLOCATABLE, DIMENSION(:):: DGRDAYS     !Degree day calculated from temperature of water

   !
   !time average values 
   !
     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: ALEAF,   & !Average leaf biomass (gC/m^2)
                                        AROOT,   & !Average root biomass (gC/m^2)
                                        ASTEM,   & !Average stem biomass (gC/m^2)       
                                       ATUBER,   & !Average tuber biomass (gC/m^2)       
                                          AEP      !Average epiphyte biomass (gC-epi/gC-leaf)

     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: APLEAF         !Averaged leaf growth rate (1/day)
     REAL(SP),ALLOCATABLE,DIMENSION(:)   :: ABMLEAF,   &   !Averaged basal metabolism rate of leaf (1/day)
                                        ABMTUBER       !Average basal metabolism rate of tuber (1/day)

     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: APEP,    &     !Average growth rate  for epiphytes (1/day)
                                        ABMEP,   &     !Average basal metabolism rate for epiphytes (1/day)
                                        APREP          !Average predation rate on epiphytes (1/day) = APREP*EP

     REAL(SP),ALLOCATABLE,DIMENSION(:) ::  ASLSH           !Average epiphute sloughing rate (1/day)
                                        
     REAL(SP),ALLOCATABLE,DIMENSION(:) ::  ANLSAV,    &     !Averaged N limitation on SAV
                                       APLSAV,    &     !Averaged P limitation on SAV
                                       ANLEPI,    &     !Averaged N limitation on epiphytes
                                       APLEPI,    &     !Averaged P limitation on epiphytes
                                       AFNSED,    &     !Averaged fraction of uptake of N from sediment
                                       AFPSED,    &     !Averaged fraction of uptake of P from sediment
                                         AFHS           !Averaged SAV growth limiation due to sulfide HS^-            
                    
     !attenuation and light limitation
     REAL(SP),ALLOCATABLE,DIMENSION(:,:) ::  AEPATN,   &    !Averaged attenuation by epiphyte
                                        AWATATN,   &    !Averaged water column attenuation of light at canopy top
                                          AFISH,   &    !Averaged light limitation for LEAF (shoots)
                                          AFIEP         !Averaged light limitation for Epiphytes      
     !net production                                          
     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: ANPPSAV,   &    !Averaged net leaf production per day (gC/m^2/day)
                                        ANPPEPI         !Averaged net epiphyte production per day (gC-epi/gC-leaf/day)

     !time averaged flux to water column
     REAL(SP),ALLOCATABLE,DIMENSION(:) ::  ADOSAVW,  &  !average DO  flux due to SAV gO2/m^2/day
                                      ADOCSAVW,  &  !average DOC flux due to SAV
                                      APOCSAVW,  &  !average POC flux due to SAV
                                      ANH4SAVW,  &  !average NH4 flux due to SAV
                                      ANO3SAVW,  &  !average NO3 flux due to SAV
                                      ADONSAVW,  &  !average DON flux due to SAV              
                                      APONSAVW,  &  !average PON flux due to SAV
                                      APO4SAVW,  &  !average PO4 flux due to SAV
                                      ADOPSAVW,  &  !average DOP flux due to SAV
                                      APOPSAVW,  &  !average POL flux due to SAV
                                       ADOEPIW,  &  !DO  flux (positive into water column ) due to epiphytes (gO2/m2/day)
                                      ANH4EPIW,  &  !NH4 flux (positive into water column) due to epiphytes (gN/m2/day)
                                      ANO3EPIW,  &  !NO3 flux (positive into water column) due to epiphytes (gN/m2/day)
                                      APO4EPIW,  &  !PO4 flux (positive into water column) due to epiphytes (gP/m2/day)
                     !following not realy used???
                                      ADOCEPIW,  &  ! DOC flux to water column due to epiphytes (gC/m^2/day)
                                      APOCEPIW,  &  ! POC flux to water column due to epiphytes (gC/m^2/day)
                                      ADONEPIW,  &  ! DON flux to water column due to epiphytes (gN/m^2/day)
                                      APONEPIW,  &  ! PON flux to water column due to epiphytes (gN/m^2/day)                  
                                      ADOPEPIW,  &  ! DOP flux to water column due to epiphytes (gP/m^2/day)
                                      APOPEPIW      ! POP flux to water column due to epiphytes (gP/m^2/day)

    ! time averaged flux to sediments
     REAL(SP),ALLOCATABLE,DIMENSION(:) :: ASEDDOSAV,  &   !DO  flux due to SAV,  positive leaving sediments
                                      ASEDPOCSAV,  &  !POC flux due to SAV, positive into sediments
                                      ASEDPONSAV,  &  !PON flux due to SAV, positive into sediments
                                      ASEDPOPSAV,  &  !POP flux due to SAV, positive into sediments
                                      ASEDNH4SAV,  &  !NH4 flux due to SAV, positive leaving sediemnts
                                      ASEDPO4SAV      !PO4 flux due to SAV, positive leaving sediments

CONTAINS
		!Subroutines:
		!	subroutine SAV_ALLOC()
		!	subroutine SAV_DEALLOC()
		!	subroutine SAV_LOADS_ALLOC()
		!	subroutine SAV_LOADS_DEALLOC()
		!	subroutine SAV_READ()
		!	subroutine SAV_COMP()

!************************************************************************
!**             S U B R O U T I N E   SAV_ALLOC                        **
!************************************************************************

   SUBROUTINE SAV_ALLOC
   !Allocate all variables related to SAV
   !much of this section is moved from wqm_modules.F
   USE MOD_LIMS, ONLY : MTLOC
   

   ALLOCATE(LEAF(MTLOC,NSAVM));         LEAF   = 0.0          !1
   ALLOCATE(ROOT(MTLOC,NSAVM));         ROOT   = 0.0          !2
   ALLOCATE(STEM(MTLOC,NSAVM));         STEM   = 0.0          !3
   ALLOCATE(TUBER(MTLOC,NSAVM));        TUBER  = 0.0          !4
   ALLOCATE(EP(MTLOC,NSAVM));           EP     = 0.0          !5

   ALLOCATE(SAVAREA(MTLOC));            SAVAREA = 0.0         !6
   ALLOCATE(SAVCELL(MTLOC));            SAVCELL = 0           !7

   ALLOCATE(NSAVSPC(MTLOC));            NSAVSPC = 0.0         !8
!  ALLOCATE(SAVWDH(MTLOC,NSAVM));       SAVWDH  = 0.0         !9
   ALLOCATE(SAVDPH(MTLOC,NSAVM));       SAVDPH  = 0.0         !10
   ALLOCATE(SAVFRAC(MTLOC,NSAVM));      SAVFRAC = 0.0         !11

!Production of SAV and Basal Metabolism of SAV   
   ALLOCATE(PLEAF(MTLOC,NSAVM));        PLEAF  = 0.0          !12
   ALLOCATE(BMLEAF(MTLOC));             BMLEAF    = 0.0       !13
   ALLOCATE(BMSTEM(MTLOC));             BMSTEM    = 0.0       !14
   ALLOCATE(BMROOT(MTLOC));             BMROOT    = 0.0       !15
   ALLOCATE(BMTUBER(MTLOC));            BMTUBER   = 0.0       !16

!Production, Basal Metabolism and Predation of Epiphytes   
   ALLOCATE(PEP(MTLOC,NSAVM));          PEP  = 0.0            !17
   ALLOCATE(BMEP(MTLOC,NSAVM));         BMEP = 0.0            !18
   ALLOCATE(PREP(MTLOC,NSAVM));         PREP = 0.0            !19

!Slough rate of Epiphytes   
   ALLOCATE(SLSH(MTLOC));               SLSH    = 0.0         !20

!Nitrogen and Phosphorus limitation to SAV growth, Fraction Nitrogen and Phosphorus uptake fron sediments by SAV
!Nitrogen and Phosphorus limitation on Epiphyte growth
!Sulfide HS^- in sediments toxicity to limited SAV growth
   ALLOCATE(NLSAV(MTLOC));              NLSAV     = 0.0       !21
   ALLOCATE(PLSAV(MTLOC));              PLSAV     = 0.0       !22
   ALLOCATE(FNSEDSAV(MTLOC));           FNSEDSAV  = 0.0       !23
   ALLOCATE(FPSEDSAV(MTLOC));           FPSEDSAV  = 0.0       !24
   ALLOCATE(NLEPI(MTLOC));              NLEPI     = 0.0       !25
   ALLOCATE(PLEPI(MTLOC));              PLEPI     = 0.0       !26
   ALLOCATE(FHS(MTLOC));                FHS       = 0.0       !27

!Light attenuation by SAV, Epiphytes, Water
!Light limitation to SAV (shoots) growth, Epiphyte growth 
!Net leaf, epiphytes production per day
   ALLOCATE(SAVATN(MTLOC,NSAVM));       SAVATN = 0.0         !28
   ALLOCATE(EPATN(MTLOC,NSAVM));        EPATN  = 0.0         !29
   ALLOCATE(WATATN(MTLOC,NSAVM));       WATATN = 0.0         !30
   ALLOCATE(FISH(MTLOC,NSAVM));         FISH   = 0.0         !31
   ALLOCATE(FIEP(MTLOC,NSAVM));         FIEP   = 0.0         !32
   ALLOCATE(NPPSAV(MTLOC,NSAVM));       NPPSAV = 0.0         !33
   ALLOCATE(NPPEPI(MTLOC,NSAVM));       NPPEPI = 0.0         !34
   
!DO, LDOC, RDOC, LPOC, RPOC flux to water column due to SAV   
   ALLOCATE(DOSAVW(MTLOC));               DOSAVW  = 0.0      !35
   ALLOCATE(LDOCSAVW(MTLOC));            LDOCSAVW = 0.0      !36
   ALLOCATE(RDOCSAVW(MTLOC));            RDOCSAVW = 0.0      !37
   ALLOCATE(LPOCSAVW(MTLOC));            LPOCSAVW = 0.0      !38
   ALLOCATE(RPOCSAVW(MTLOC));            RPOCSAVW = 0.0      !39
!NH4, NO4, LDON, RDON, LPON, RPON flux to water column due to SAV
   ALLOCATE(NH4SAVW(MTLOC));            NH4SAVW   = 0.0      !40
   ALLOCATE(NO3SAVW(MTLOC));            NO3SAVW   = 0.0      !41
   ALLOCATE(LDONSAVW(MTLOC));           LDONSAVW  = 0.0      !42
   ALLOCATE(RDONSAVW(MTLOC));           RDONSAVW  = 0.0      !43
   ALLOCATE(LPONSAVW(MTLOC));           LPONSAVW  = 0.0      !44
   ALLOCATE(RPONSAVW(MTLOC));           RPONSAVW  = 0.0      !45
!PO4, LDOP, RDOP, LPOP, RPOP flux to water column due to SAV
   ALLOCATE(PO4SAVW(MTLOC));            PO4SAVW   = 0.0      !46
   ALLOCATE(LDOPSAVW(MTLOC));           LDOPSAVW  = 0.0      !47
   ALLOCATE(RDOPSAVW(MTLOC));           RDOPSAVW  = 0.0      !48
   ALLOCATE(LPOPSAVW(MTLOC));           LPOPSAVW  = 0.0      !49
   ALLOCATE(RPOPSAVW(MTLOC));           RPOPSAVW  = 0.0      !50

!DO, LDOC, RDOC, RPOC flux to water colum due to epiphytes  
   ALLOCATE(DOEPIW(MTLOC));              DOEPIW   = 0.0      !51
   ALLOCATE(LDOCEPIW(MTLOC));            LDOCEPIW = 0.0      !52
   ALLOCATE(RDOCEPIW(MTLOC));            RDOCEPIW = 0.0      !53
   ALLOCATE(LPOCEPIW(MTLOC));            LPOCEPIW = 0.0      !54
   ALLOCATE(RPOCEPIW(MTLOC));            RPOCEPIW = 0.0      !55
!NH4, NO3, LDON, RDON, LPON, RPON flux to water column due to epiphytes   
   ALLOCATE(NH4EPIW(MTLOC));             NH4EPIW    = 0.0    !56
   ALLOCATE(NO3EPIW(MTLOC));             NO3EPIW    = 0.0    !57
   ALLOCATE(LDONEPIW(MTLOC));            LDONEPIW   = 0.0    !58
   ALLOCATE(RDONEPIW(MTLOC));            RDONEPIW   = 0.0    !59
   ALLOCATE(LPONEPIW(MTLOC));            LPONEPIW   = 0.0    !60
   ALLOCATE(RPONEPIW(MTLOC));            RPONEPIW   = 0.0    !61
!PO4, LDOP, RDOP, LPOP, RPOP flux to water column due to epiphytes   
   ALLOCATE(PO4EPIW(MTLOC));             PO4EPIW    = 0.0    !62
   ALLOCATE(LDOPEPIW(MTLOC));            LDOPEPIW   = 0.0    !63
   ALLOCATE(RDOPEPIW(MTLOC));            RDOPEPIW   = 0.0    !64
   ALLOCATE(LPOPEPIW(MTLOC));            LPOPEPIW   = 0.0    !65
   ALLOCATE(RPOPEPIW(MTLOC));            RPOPEPIW   = 0.0    !66

!loading flux : POC, PON, POP (positive into sediments)      to sediments
!return  flux : DO, NH4, PO4  (positive leaving sediments) from sediments 
!applied to sediments due to SAV
   
   ALLOCATE(SEDPOCSAV(MTLOC));          SEDPOCSAV = 0.0   !67
   ALLOCATE(SEDPONSAV(MTLOC));          SEDPONSAV = 0.0   !68
   ALLOCATE(SEDPOPSAV(MTLOC));          SEDPOPSAV = 0.0   !69
   ALLOCATE(SEDDOSAV(MTLOC));           SEDDOSAV  = 0.0   !70
   ALLOCATE(SEDNH4SAV(MTLOC));          SEDNH4SAV = 0.0   !71
   ALLOCATE(SEDNO3SAV(MTLOC));          SEDNO3SAV = 0.0   !71-1
   ALLOCATE(SEDPO4SAV(MTLOC));          SEDPO4SAV = 0.0   !72

   !misc
   ALLOCATE(DGRDAYS(MTLOC));            DGRDAYS = 0.0     !73

   
   !SAV model time averaged quantities
   
   !state varaiables
   ALLOCATE(ALEAF(MTLOC,NSAVM));        ALEAF   = 0.0   
   ALLOCATE(AROOT(MTLOC,NSAVM));        AROOT   = 0.0     
   ALLOCATE(ASTEM(MTLOC,NSAVM));        ASTEM   = 0.0   
   ALLOCATE(ATUBER(MTLOC,NSAVM));       ATUBER  = 0.0
   ALLOCATE(AEP(MTLOC,NSAVM));          AEP     = 0.0

   !growth and mortality of SAV
   ALLOCATE(APLEAF(MTLOC,NSAVM));       APLEAF  = 0.0
   ALLOCATE(ABMLEAF(MTLOC));           ABMLEAF  = 0.0
   ALLOCATE(ABMTUBER(MTLOC));         ABMTUBER  = 0.0   
   
   !growth and mortality of Epiphytes
   ALLOCATE(APEP(MTLOC,NSAVM));         APEP  = 0.0   
   ALLOCATE(ABMEP(MTLOC,NSAVM));        ABMEP = 0.0
   ALLOCATE(APREP(MTLOC,NSAVM));        APREP = 0.0
   ALLOCATE(ASLSH(MTLOC));              ASLSH = 0.0
   
   !nutrient, light, HS^- limitation on SAV and Epiphytes growth
   ALLOCATE(ANLSAV(MTLOC));             ANLSAV = 0.0       
   ALLOCATE(APLSAV(MTLOC));             APLSAV = 0.0
   ALLOCATE(ANLEPI(MTLOC));             ANLEPI = 0.0
   ALLOCATE(APLEPI(MTLOC));             APLEPI = 0.0
   ALLOCATE(AFNSED(MTLOC));             AFNSED = 0.0 !fraction uptake of N from sediment
   ALLOCATE(AFPSED(MTLOC));             AFPSED = 0.0 !fraction uptake of P from sediment
   ALLOCATE(AFHS(MTLOC));               AFHS   = 0.0         

   !Attenuation, Light limitation and net production 
   ALLOCATE(AEPATN(MTLOC,NSAVM));       AEPATN  = 0.0     
   ALLOCATE(AWATATN(MTLOC,NSAVM));      AWATATN = 0.0 
   ALLOCATE(AFISH(MTLOC,NSAVM));        AFISH   = 0.0
   ALLOCATE(AFIEP(MTLOC,NSAVM));        AFIEP   = 0.0   
   ALLOCATE(ANPPSAV(MTLOC,NSAVM));      ANPPSAV = 0.0   !net production of SAV
   ALLOCATE(ANPPEPI(MTLOC,NSAVM));      ANPPEPI  = 0.0  !net production of Epiphytes
   
 !time averaged flux to water column due to SAV and epiphytes
   ALLOCATE(ADOCSAVW(MTLOC));      ADOCSAVW = 0.0   !average of LDOCSAVW and RDOCSAVW accumulated over time
   ALLOCATE(APOCSAVW(MTLOC));      APOCSAVW = 0.0   
   ALLOCATE(ADOCEPIW(MTLOC));      ADOCEPIW = 0.0
   ALLOCATE(APOCEPIW(MTLOC));      APOCEPIW = 0.0

   ALLOCATE(ANH4SAVW(MTLOC));      ANH4SAVW = 0.0
   ALLOCATE(ANO3SAVW(MTLOC));      ANO3SAVW = 0.0
   ALLOCATE(ADONSAVW(MTLOC));      ADONSAVW = 0.0
   ALLOCATE(APONSAVW(MTLOC));      APONSAVW = 0.0
   ALLOCATE(ANH4EPIW(MTLOC));      ANH4EPIW = 0.0
   ALLOCATE(ANO3EPIW(MTLOC));      ANO3EPIW = 0.0
   ALLOCATE(ADONEPIW(MTLOC));      ADONEPIW = 0.0
   ALLOCATE(APONEPIW(MTLOC));      APONEPIW = 0.0

   ALLOCATE(APO4SAVW(MTLOC));      APO4SAVW = 0.0
   ALLOCATE(ADOPSAVW(MTLOC));     ADOPSAVW  = 0.0
   ALLOCATE(APOPSAVW(MTLOC));      APOPSAVW = 0.0
   ALLOCATE(APO4EPIW(MTLOC));    APO4EPIW   = 0.0
   ALLOCATE(ADOPEPIW(MTLOC));      ADOPEPIW = 0.0
   ALLOCATE(APOPEPIW(MTLOC));     APOPEPIW  = 0.0

   ALLOCATE(ADOSAVW(MTLOC));     ADOSAVW    = 0.0
   ALLOCATE(ADOEPIW(MTLOC));     ADOEPIW    = 0.0

   !time averaged flux to sediment due to SAV and epiphytes
   ALLOCATE(ASEDDOSAV(MTLOC));    ASEDDOSAV = 0.0
   ALLOCATE(ASEDPOCSAV(MTLOC));  ASEDPOCSAV = 0.0   
   ALLOCATE(ASEDPONSAV(MTLOC));  ASEDPONSAV = 0.0
   ALLOCATE(ASEDNH4SAV(MTLOC));  ASEDNH4SAV = 0.0
   ALLOCATE(ASEDPOPSAV(MTLOC));  ASEDPOPSAV = 0.0
   ALLOCATE(ASEDPO4SAV(MTLOC));  ASEDPO4SAV = 0.0
      
   END SUBROUTINE SAV_ALLOC

   
!************************************************************************
!**             S U B R O U T I N E   SAV_LOADS_ALLOC                  **
!************************************************************************

   SUBROUTINE SAV_LOADS_ALLOC
   !Allocate all variables related to SAV
   !much of this section is moved from wqm_modules.F
   USE MOD_LIMS, ONLY : MTLOC
   
        ALLOCATE(  DOSAVW(MTLOC));     DOSAVW  = 0.0       !35
        ALLOCATE(LDOCSAVW(MTLOC));     LDOCSAVW = 0.0      !36
        ALLOCATE(RDOCSAVW(MTLOC));     RDOCSAVW = 0.0      !37
        ALLOCATE(LPOCSAVW(MTLOC));     LPOCSAVW = 0.0      !38
        ALLOCATE(RPOCSAVW(MTLOC));     RPOCSAVW = 0.0      !39
   END SUBROUTINE SAV_LOADS_ALLOC
   
!************************************************************************
!**             S U B R O U T I N E   SAV_LOADS_DEALLOC                **
!************************************************************************

   SUBROUTINE SAV_LOADS_DEALLOC
   !Allocate all variables related to SAV
   !much of this section is moved from wqm_modules.F
     IF(  ALLOCATED(DOSAVW)) DEALLOCATE (DOSAVW)      !35
     IF(ALLOCATED(LDOCSAVW)) DEALLOCATE (LDOCSAVW)    !36
     IF(ALLOCATED(RDOCSAVW)) DEALLOCATE (RDOCSAVW)    !37
     IF(ALLOCATED(LPOCSAVW)) DEALLOCATE (LPOCSAVW)    !38
     IF(ALLOCATED(RPOCSAVW)) DEALLOCATE (RPOCSAVW)    !39
   END SUBROUTINE SAV_LOADS_DEALLOC
   

!************************************************************************
!**             S U B R O U T I N E   SAV_DEALLOC                      **
!************************************************************************
   
   SUBROUTINE SAV_DEALLOC
   
   !Deallocation of all variables related to SAV
   !much of this section is moved from wqm_main.F
   !
     IF(ALLOCATED(LEAF)) DEALLOCATE (LEAF)            !1
     IF(ALLOCATED(ROOT)) DEALLOCATE (ROOT)            !2
     IF(ALLOCATED(STEM)) DEALLOCATE (STEM)            !3
     IF(ALLOCATED(TUBER)) DEALLOCATE (TUBER)          !4
     IF(ALLOCATED(EP)) DEALLOCATE (EP)                !5

     IF(ALLOCATED(SAVAREA)) DEALLOCATE (SAVAREA)      !6
     IF(ALLOCATED(SAVCELL)) DEALLOCATE (SAVCELL)      !7

     IF(ALLOCATED(NSAVSPC)) DEALLOCATE (NSAVSPC)      !8  
    !IF(ALLOCATED(SAVWDH)) DEALLOCATE (SAVWDH)        !9 !this is now completely deprecated and replaced by SAVFRAC   
     IF(ALLOCATED(SAVDPH)) DEALLOCATE (SAVDPH)        !10
     IF(ALLOCATED(SAVFRAC)) DEALLOCATE (SAVFRAC)      !11

     IF(ALLOCATED(PLEAF)) DEALLOCATE (PLEAF)          !12
     IF(ALLOCATED(BMLEAF)) DEALLOCATE (BMLEAF)        !13
     IF(ALLOCATED(BMSTEM)) DEALLOCATE (BMSTEM)        !14
     IF(ALLOCATED(BMROOT)) DEALLOCATE (BMROOT)        !15
     IF(ALLOCATED(BMTUBER)) DEALLOCATE (BMTUBER)      !16

     IF(ALLOCATED(PEP)) DEALLOCATE (PEP)              !17
     IF(ALLOCATED(BMEP)) DEALLOCATE (BMEP)            !18
     IF(ALLOCATED(PREP)) DEALLOCATE (PREP)            !19
     IF(ALLOCATED(SLSH)) DEALLOCATE (SLSH)            !20
     
     IF(ALLOCATED(NLSAV)) DEALLOCATE (NLSAV)          !21
     IF(ALLOCATED(PLSAV)) DEALLOCATE (PLSAV)          !22
     IF(ALLOCATED(FNSEDSAV)) DEALLOCATE (FNSEDSAV)    !23
     IF(ALLOCATED(FPSEDSAV)) DEALLOCATE (FPSEDSAV)    !24
     IF(ALLOCATED(NLEPI)) DEALLOCATE (NLEPI)          !25
     IF(ALLOCATED(PLEPI)) DEALLOCATE (PLEPI)          !26
     IF(ALLOCATED(FHS)) DEALLOCATE (FHS)              !27
   
     IF(ALLOCATED(EPATN)) DEALLOCATE (EPATN)          !28
     IF(ALLOCATED(SAVATN)) DEALLOCATE (SAVATN)        !29
     IF(ALLOCATED(WATATN)) DEALLOCATE (WATATN)        !30
     IF(ALLOCATED(FISH)) DEALLOCATE (FISH)            !31
     IF(ALLOCATED(FIEP)) DEALLOCATE (FIEP)            !32
     IF(ALLOCATED(NPPSAV)) DEALLOCATE (NPPSAV)        !33
     IF(ALLOCATED(NPPEPI)) DEALLOCATE (NPPEPI)        !34

     IF(ALLOCATED(DOSAVW)) DEALLOCATE (DOSAVW)        !35
     IF(ALLOCATED(LDOCSAVW)) DEALLOCATE (LDOCSAVW)    !36
     IF(ALLOCATED(RDOCSAVW)) DEALLOCATE (RDOCSAVW)    !37
     IF(ALLOCATED(LPOCSAVW)) DEALLOCATE (LPOCSAVW)    !38
     IF(ALLOCATED(RPOCSAVW)) DEALLOCATE (RPOCSAVW)    !39
     IF(ALLOCATED(NH4SAVW)) DEALLOCATE (NH4SAVW)      !40
     IF(ALLOCATED(NO3SAVW)) DEALLOCATE (NO3SAVW)      !41
     IF(ALLOCATED(LDONSAVW)) DEALLOCATE (LDONSAVW)    !42
     IF(ALLOCATED(RDONSAVW)) DEALLOCATE (RDONSAVW)    !43
     IF(ALLOCATED(LPONSAVW)) DEALLOCATE (LPONSAVW)    !44
     IF(ALLOCATED(RPONSAVW)) DEALLOCATE (RPONSAVW)    !45
     IF(ALLOCATED(PO4SAVW)) DEALLOCATE (PO4SAVW)      !46
     IF(ALLOCATED(LDOPSAVW)) DEALLOCATE (LDOPSAVW)    !47
     IF(ALLOCATED(RDOPSAVW)) DEALLOCATE (RDOPSAVW)    !48
     IF(ALLOCATED(LPOPSAVW)) DEALLOCATE (LPOPSAVW)    !49
     IF(ALLOCATED(RPOPSAVW)) DEALLOCATE (RPOPSAVW)    !50

     IF(ALLOCATED(DOEPIW)) DEALLOCATE (DOEPIW)        !51
     IF(ALLOCATED(LDOCEPIW)) DEALLOCATE (LDOCEPIW)    !52
     IF(ALLOCATED(RDOCEPIW)) DEALLOCATE (RDOCEPIW)    !53
     IF(ALLOCATED(LPOCEPIW)) DEALLOCATE (LPOCEPIW)    !54
     IF(ALLOCATED(RPOCEPIW)) DEALLOCATE (RPOCEPIW)    !55
     IF(ALLOCATED(NH4EPIW)) DEALLOCATE (NH4EPIW)      !56
     IF(ALLOCATED(NO3EPIW)) DEALLOCATE (NO3EPIW)      !57
     IF(ALLOCATED(LDONEPIW)) DEALLOCATE (LDONEPIW)    !58
     IF(ALLOCATED(RDONEPIW)) DEALLOCATE (RDONEPIW)    !59
     IF(ALLOCATED(LPONEPIW)) DEALLOCATE (LPONEPIW)    !60
     IF(ALLOCATED(RPONEPIW)) DEALLOCATE (RPONEPIW)    !61
     IF(ALLOCATED(PO4EPIW)) DEALLOCATE (PO4EPIW)      !62
     IF(ALLOCATED(LDOPEPIW)) DEALLOCATE (LDOPEPIW)    !63
     IF(ALLOCATED(RDOPEPIW)) DEALLOCATE (RDOPEPIW)    !64
     IF(ALLOCATED(LPOPEPIW)) DEALLOCATE (LPOPEPIW)    !65
     IF(ALLOCATED(RPOPEPIW)) DEALLOCATE (RPOPEPIW)    !66

     
     IF(ALLOCATED(SEDPOCSAV)) DEALLOCATE (SEDPOCSAV)  !67
     IF(ALLOCATED(SEDDOSAV)) DEALLOCATE (SEDDOSAV)    !68
     IF(ALLOCATED(SEDNH4SAV)) DEALLOCATE (SEDNH4SAV)  !69
     IF(ALLOCATED(SEDNO3SAV)) DEALLOCATE (SEDNO3SAV)  !69-1
     IF(ALLOCATED(SEDPO4SAV)) DEALLOCATE (SEDPO4SAV)  !70
     IF(ALLOCATED(SEDPONSAV)) DEALLOCATE (SEDPONSAV)  !71
     IF(ALLOCATED(SEDPOPSAV)) DEALLOCATE (SEDPOPSAV)  !72

     IF(ALLOCATED(DGRDAYS)) DEALLOCATE (DGRDAYS)      !73

! SAV time averaged quantities
     
     IF(ALLOCATED(AEPATN))  DEALLOCATE (AEPATN)  
     IF(ALLOCATED(AWATATN)) DEALLOCATE (AWATATN)
     IF(ALLOCATED(AFISH))   DEALLOCATE (AFISH)
     IF(ALLOCATED(AFIEP))   DEALLOCATE (AFIEP)   
     IF(ALLOCATED(ANPPSAV)) DEALLOCATE (ANPPSAV)
     IF(ALLOCATED(ANPPEPI))  DEALLOCATE (ANPPEPI)
    
     IF(ALLOCATED(ADOSAVW)) DEALLOCATE (ADOSAVW)     
     IF(ALLOCATED(ADOCSAVW)) DEALLOCATE (ADOCSAVW)
     IF(ALLOCATED(APOCSAVW)) DEALLOCATE (APOCSAVW)
     IF(ALLOCATED(ANH4SAVW)) DEALLOCATE (ANH4SAVW)
     IF(ALLOCATED(ANO3SAVW)) DEALLOCATE (ANO3SAVW)
     IF(ALLOCATED(ADONSAVW)) DEALLOCATE (ADONSAVW)
     IF(ALLOCATED(APONSAVW)) DEALLOCATE (APONSAVW)
     IF(ALLOCATED(APO4SAVW)) DEALLOCATE (APO4SAVW)
     IF(ALLOCATED(ADOPSAVW)) DEALLOCATE (ADOPSAVW)
     IF(ALLOCATED(APOPSAVW)) DEALLOCATE (APOPSAVW)   

     IF(ALLOCATED(ADOEPIW)) DEALLOCATE (ADOEPIW)
     IF(ALLOCATED(ADOCEPIW)) DEALLOCATE (ADOCEPIW)   
     IF(ALLOCATED(APOCEPIW)) DEALLOCATE (APOCEPIW)     
     IF(ALLOCATED(ANH4EPIW)) DEALLOCATE (ANH4EPIW)
     IF(ALLOCATED(ANO3EPIW)) DEALLOCATE (ANO3EPIW)
     IF(ALLOCATED(ADONEPIW)) DEALLOCATE (ADONEPIW)
     IF(ALLOCATED(APONEPIW)) DEALLOCATE (APONEPIW)
     IF(ALLOCATED(APO4EPIW)) DEALLOCATE (APO4EPIW)
     IF(ALLOCATED(ADOPEPIW)) DEALLOCATE (ADOPEPIW)
     IF(ALLOCATED(APOPEPIW)) DEALLOCATE (APOPEPIW)

     IF(ALLOCATED(ASEDDOSAV))  DEALLOCATE (ASEDDOSAV)
     IF(ALLOCATED(ASEDPOCSAV)) DEALLOCATE (ASEDPOCSAV)
     IF(ALLOCATED(ASEDPONSAV)) DEALLOCATE (ASEDPONSAV)
     IF(ALLOCATED(ASEDPOPSAV)) DEALLOCATE (ASEDPOPSAV)
     IF(ALLOCATED(ASEDNH4SAV)) DEALLOCATE (ASEDNH4SAV)   
     IF(ALLOCATED(ASEDPO4SAV)) DEALLOCATE (ASEDPO4SAV)

     
   END SUBROUTINE SAV_DEALLOC


!************************************************************************
!**             S U B R O U T I N E   S A V _ R E A D                  **
!************************************************************************

   SUBROUTINE SAV_READ
   
	   USE MOD_SIZES, ONLY : MGL
	   USE MOD_WQM, ONLY :B, SVOFN
	   
       USE MOD_FILEINFO, ONLY : 		&
				!,DIA			& !
				!,CBC 			& !
				!,S1			& !
				!,S2			& !
				!,S3			& !                
				!,BFI			& !
	            !,BAI           & !
				!,MET			& !
				!,BFO			& !
				!,KEI			& !
				!,ATM			& !
				!,STL			& !
				!,AGR			& !
				 SVI			& !
				,SVO	!		& !
				!,KFL			& !
				!,ZOO			& !
				!,ZFO			& !
				!,ALO      		&!
				!,RSO			&!
				!,SNP			&!
				!,PLT			&!
				!,APL 			&!
				!,TFL 			&!
				!,OPL			&!
				!,SFI			&!
				!,SFO			&!
                !,MAP 			&!
				!,ICI 			&!
				!,ICO			&!
				!,MRL			&!
				!,MBL			&!
				!,RSI			&!
                !,UNIT_LINKAGE	&!
				!,UNIT_STN		&!
				!,UNIT_HIS		&!           
				!,CNAME			&!
				!,CONFN
				
		USE MOD_LIMS, ONLY : MLOC

   !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
   USE MOD_CONTROL, ONLY : 		&
			SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			,MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			,PAR        !	&           !!TRUE IF MULTIPROCESSOR RUN
			!,CASENAME  	&   		!!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
			!,CASETITLE  	&  			!!CASE TITLE                                 
			!,HMAX       	&  			!!GLOBAL MAXIMUM DEPTH
			!,HMIN       	&  			!!GLOBAL MINIMUM DEPTH
			!,UMOL       	&  			!!VERTICAL DIFFUSION COEFFICIENT
			!,HORCON     	&  			!!HORIZONTAL DIFFUSION COEFFICIENT
			!,DTI        	&  			!!internal time step
			!,HORZMIX    	&   		!!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
			!,FILENUMBER	&			!!
			!,PREFF			&			!!
			!,INPDIR		&			!!
			!,GEOAREA		&			!!
			!,RIV_FILENUMBER	&			!!
            !,INFLOW_TYPE   	&			!!SPECIFIED RIVER INFLOW TYPE (edge/node) 
            !,POINT_ST_TYPE 	&			!!(calculated/specified)
            !,PNT_SOURCE    	&			!!point_source
            !,DAY				&
			!,in_jday		
		     
#  if defined (MULTIPROCESSOR)   
   USE MOD_PAR, ONLY : NLID
#  endif
   IMPLICIT NONE
   CHARACTER(LEN=72) :: TITLE(6)
   CHARACTER(LEN=8)  :: TIMVAR, PRINTIT
   REAL(SP) :: KTPS1,        &  !temperature control coefficenet for growth rate of SAV [(1/degc)^2] at T < T_optimal
           KTPS2,        &  !temperature control coefficenet for growth rate of SAV [(1/degc)^2] at T > T_optimal
           KTBMSAV,      &  !temperature control coefficient for basal metabolism (basal respiration for root, leaf, stem, tube [1/degC]
           KTPE1,        &  !temperature control coefficient on photosynthesis rate of SAV species 1 [(1/degC)^2]
           KTPE2,        &  !temperature control coefficient on photosynthesis rate of SAV species 2 [(1/degC)^2]
           KTBME,        &  !temperature control coefficient for epiphyte basal metabolism [1/degC]
           KTPRE            !temperature control coefficient for epiphyte predation loss rate [1/degC]

! vjp 8/23/05   added these declarations
   INTEGER :: I, J  &
              , N         ! Wen Long: N is never used here

   INTEGER:: NSAV           !Number of SAV species 

   REAL(SP)    :: TRPMSAV,   &  !Optimal (reference) temperature for SAV growth (degC)
              TRBMSAV,   &  !Reference temperature for SAV basal metabolism (degC)
              TRPMEPI,   &  !Optimal (reference) temperature for epiphyte growth (degC)
              TRBMEPI,   &  !Reference temperature for epiphyte metabolism (degC)
              TRPREPI,   &  !Reference temperature for epiphyte predation (degC)
              TLOOK         !Tempoary variable for temperature

   !Wen Long added following temporary variables              
   INTEGER, ALLOCATABLE :: SAVCELL_GL(:)  &    !global TCE ID
                          ,NSAVSPC_GL(:)       !number of species in each global TCE
                             
   REAL(SP), ALLOCATABLE :: SAVFRAC_GL(:,:)        !global SAV coverage in global TCE
   REAL(SP), ALLOCATABLE ::  SAVDPH_GL(:,:)        !global SAV depth in global TCE

   INTEGER :: NSAVCELL_GL                      !global total number of SAV TCE's
   INTEGER NSAVM_TMP
   
! TITLE CARDS

   READ(SVI,1000)
!1000 FORMAT(///)          !go to the end of 4th line
   
   READ(SVI,1010) (TITLE(J),J=1,6)  !read one line
!1010 FORMAT(A72)
! READ SAV PROPERTIES
   
   READ(SVI, '(//8X,I8)') NSAVM_TMP     !go down two lines and read the third line
   
!1040 FORMAT(//8X,2A8)    !go down two lines and read the third line
   READ(SVI,1040) SPNAM   !read name of SAV species  
     !growth rate, respiration ratel, basal metabolism, DO consumption fraction, sulfide toxicity coeffient on SAV
     !cut-off salinity for salinity toxicity on freshwater species of SAV. (ramp up respiration (mortality) with high salt)
   READ(SVI,1030) PMSAV,PRSPSAVB,BMSAV,BMTBRREF,FDOSR,KHSTOX,SALMAX  
     !C, N, P to SAV dryweight ratios (gC/gDW), (gN/gDW), (gP/gDW)
   READ(SVI,1030) ACDWSAV,ANDWSAV,APDWSAV
!1030 FORMAT(//(8X,9F8.0))  !go down two lines and read the third line
   
! COMPUTE ADDITIONAL STOICHIOMETRY OF SAV

   ANCSAV = ANDWSAV/ACDWSAV   !N to C ratio  gN/gC in SAV (WL: should it be different for leaf, root and epiphytes ?)
   APCSAV = APDWSAV/ACDWSAV   !P to C ratio  gP/gC in SAV (WL: should it be different for leaf, root and epiphytes ?)

   !Half saturation constant for N uptake by leaf (shoots) (gN/m^3)
   !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) 
   !Half saturation constant for N uptake by root (gN/m^3)
   !Half saturation constant for P uptake by leaf (KHPLEAF) and root (KHPROOT) (gP/m^3)
   
   READ(SVI,1030) KHNLEAF,KHNPSAV,KHNROOT,KHPLEAF,KHPROOT   !go down two lines and read the third line

   !PI curve slope alpha for SAV, ALPHSAV, [(gC/gDW)/(Einstein/m^2)]
   !Light attenuation due to SAV leaf and stem, KESAV (gC/m^3)  
   !ALAC = carbon weight per unit leaf area (m^2-leaf-area/gC-leaf)
   !ACAN = basic canopy height (m)
   !BCAN = effect of SAV leaf and stem biomass on canopy height (m^3/gC)
   
   READ(SVI,1030) ALPHSAV,KESAV,ALAC,ACAN,BCAN  !go down two lines and read the third line
   
   !Optimal temperature for SAV growth (degC) and reference temperature for SAV basal metabolism
   READ(SVI,1030) TRPMSAV,TRBMSAV      !go down two lines and read the third line
   
   !effect of temperature below optimal temperature and above optimal temperature for SAV growth 
   !effect of temperature on SAV basal metabolism
   READ(SVI,1030) KTPS1,KTPS2,KTBMSAV  !go down two lines and read the third line

! READ PIECEWISE FUNCTIONS DESCRIBING PARTITIONING OF PRODUCTION,
! TRANSFER FROM TUBERS TO LEAVES
! VALUES ARE ENTERED AT INTERVALS OF 10 DEGREE DAYS FROM ZERO TO
! SENESCENCE

   !number of data entries for FPLEAF, FPSTEM, FPROOT, FPTUBER, TRTBRLF
   !if these change over time of the year (temperature). It is specified every 1/10 of DGRDAYS
   !where DGRDAYS is caluclated as TDIFF*DLTDY, and TDIFF = max(T-3.0,0)  !Wen Long: I don't like this
   !                                                    !            conversion of temperature to days
   !                                                    !      spring and autumn have same temperature
   !                                                    !but may not have same entries for FPLEAF etc here
   
   READ(SVI,1042) NDDI, PRINTIT
!1042 FORMAT(//8X,I8,A8)                   !go down two lines and read the third
   READ(SVI,1032) (FPLEAF(J), FPSTEM(J),                     &
                   FPROOT(J), FPTUBER(J), TRTBRLF(J),J=0,NDDI-1)
!1032 FORMAT(//(8X,5F8.0))  !go down two lines and read the third on ward until finish 
! READ DISTRIBUTION OF SAV BIOMASS AFTER MORTALITY

   !Nitrogen distribution to NH4, LDON, RDON, LPON, RPON 
   READ(SVI,1030) FNISAV,FNLDSAV,FNRDSAV,FNLPSAV,FNRPSAV ! go down two lines and read the third line
   
   !3G partitioning of PON flux into sediments
   READ(SVI,1030) FRPONSAV(1),FRPONSAV(2),FRPONSAV(3)    ! go down two lines and read the third line
   
   !Phosphorus distribution to PO4, LDOP, RDOP, LPOL, RPOP
   READ(SVI,1030) FPISAV,FPLDSAV,FPRDSAV,FPLPSAV,FPRPSAV ! go down two lines and read the third line
   
   !3G partitioning of POP flux into sediments
   READ(SVI,1030) FRPOPSAV(1),FRPOPSAV(2),FRPOPSAV(3)    ! go down two lines and read the third line
   
   !carbon distribution into DO consumption (CO2), LDOC, RDOC, LPOC, RPOC
   READ(SVI,1030) FDOSAV,FCLDSAV,FCRDSAV,FCLPSAV,FCRPSAV ! go down two lines and read the third line

   !3G partitioning of POC flux into sediments
   READ(SVI,1030) FRPOCSAV(1),FRPOCSAV(2),FRPOCSAV(3)    ! go down two lines and read the third line

! READ EFFECT OF SAV ON NET SETTLING TO SEDIMENTS

   !fixed solids, LPOM, RPOM, algae1, algae 2, algae 3, and particulate biogenic silica
   READ(SVI,1030) WSSSAV,WSLSAV,WSRSAV,WS1SAV,WS2SAV,WS3SAV,WSUSAV  ! go down two lines and read the third line

! READ EPIPHYTE NAME AND PROPERTIES
   !Name of epipytes
!1040 FORMAT(//8X,2A8)   
   READ(SVI,1040) EPINAM       !go down two lines and read the third
   
   !growth rate, photoresprication rate, basal metabolism and predation rate
   READ(SVI,1030) PMEPI,PRSPEPI,BMEPI,PREPI  ! go down two lines and read the third
   
   !C to CHLA ratio, DW to C ratio, N to C ratio, P to C ratio of Epiphytes
   READ(SVI,1030) CCHLEPI,ADWCEPI,ANCEPI,APCEPI  !go down two lines and read the third
   
   !Nitrogen uptake half constant, NH4 preference constant, Phosphors uptake half constant
   READ(SVI,1030) KHNEPI,KHNPEPI,KHPEPI           !go down two lines and read the third
   
   !PI curve slope for Epiphyte, Light attenuation of epiphyte
   READ(SVI,1030) ALPHEPI,KEEPI                   !go down two lines and read the third
   
   !primary production optimal temperature, basal metabolism reference termperature, predation for epiphytes
   READ(SVI,1030) TRPMEPI,TRBMEPI,TRPREPI         !go down two lines and read the third
   
   !temperature effect of growth (production) below and above optimal temperature TRPMEPI, 
   !temperature effect on basal metabolism and predation for epiphytes
   READ(SVI,1030) KTPE1,KTPE2,KTBME,KTPRE         !go down two lines and read the third
   
   !half saturation density of epiphyte for epiphyte growth
   READ(SVI,1030) KHEP                            !go down two lines and read the third

! READ DISTRIBUTION OF EPIPHYTE BIOMASS AFTER MORTALITY

   !nitrogen distribution to inoragnic N (NH4), LDON, RDON, LPON, RPON upon basal metabolism loss
   READ(SVI,1030) FNIEPI,FNLDEPI,FNRDEPI,FNLPEPI,FNRPEPI    !go down two lines and read the third
   
   !nitrogen distribution to inorganic N (NH4), LDON, RDON, LPON, RPON upon predation loss
   READ(SVI,1030) FNIPEP,FNLDPEP,FNRDPEP,FNLPPEP,FNRPPEP    !go down two lines and read the third
   
   !phosphorus distribution to inorganic P (PO4), LDOP, RDOP, LPOP, RPOP upon basal metabolism loss
   READ(SVI,1030) FPIEPI,FPLDEPI,FPRDEPI,FPLPEPI,FPRPEPI    !go down two lines and read the third
   
   !phosphorus distribution to inorganic P (PO4), LDOP, RDOP, LPOP, RPOP upon predation loss
   READ(SVI,1030) FPIPEP,FPLDPEP,FPRDPEP,FPLPPEP,FPRPPEP    !go down two lines and read the third
   
   !Carbon distribution to LDOC, RDOC up on basal metabolism 
   READ(SVI,1030) FCLDEPI,FCRDEPI                            !go down two lines and read the third
                              !WLong: missing FCLPEPI and FCRPEPI to routt to LPOC and RPOC ??? (seems not coded)
                              !FCLPEPI and FCRPEPI should be added (Wen Long),as these are coded for nitrogen
                              !and there is no reason not to do the same for carbon. Also below LPOP and RPOP 
                              !have their share of carbon from predation loss. No reason that metabolism should be 
                              !different
   
                              !also missing FDOEPI compared to FDOSAV for SAV
                              !code assumes that no carbon will be routed to LPOP and RPOP due to metabolism
                              !and also FDOEPI=1-FCLDEPI-FCRDEPI, and FDOEPI is only implicitly used in calculating 
                              !DOEPIW

   !Carbon distribution to LDOP, RDOP, LPOP, RPOP up on predation loss for epiphytes
   READ(SVI,1030) FCLDPEP,FCRDPEP,FCLPPEP,FCRPPEP    !go down two lines and read the third
      
! READ NUMBER OF CELLS WITH SAV, MAXIMUM NUMBER OF DEPTH INCREMENTS

   READ(SVI,1020) NSAVCELL_GL, NSAVSPCM          !go down two lines and read the third
                                                 !Wen Long: Note that NSAVSPCM is never used
!1020 FORMAT(//(8X,9I8))   
      
! READ NUMBER OF DEPTH INCREMENTS, CHARACTERISTIC WIDTHS, AND DEPTHS

   DO I=1, MLOC             
     SAVCELL(I) = 0
     NSAVSPC(I) = 0
   ENDDO
   DO J=1, NSAVM
     DO I=1, MLOC           
!       SAVWDH(I,J) = 0.0
       SAVFRAC(I,J) = 0.0   !Wen Long: replace SAVWDH with SAVFRAC
        SAVDPH(I,J) = 0.0
     ENDDO
   ENDDO

   READ(SVI,1006)
!1006 FORMAT(/)   

!
!WLong: Parallized the code here by distributing global quantity to each subdomain that an individual computer node 
!       is assigned.  NSAVCELL is total number of SAV cells in a subdomain of parallel computer node
!
   
   ALLOCATE(SAVCELL_GL(NSAVCELL_GL)); SAVCELL_GL=0;        !global TCE ID
   ALLOCATE(NSAVSPC_GL(MGL));         NSAVSPC_GL=0;        !number of species in each global TCE
   ALLOCATE(SAVFRAC_GL(MGL,NSAVM));   SAVFRAC_GL=0.0;      !global SAV coverage in global TCE
   ALLOCATE( SAVDPH_GL(MGL,NSAVM));    SAVDPH_GL=0.0;      !global SAV depth in global TCE

   DO I=1,NSAVCELL_GL             !Loop through all SAV cells
     READ(SVI,*) NSAV, B
     SAVCELL_GL(I) = B            !The overall box number of I'th SAV cell
     NSAVSPC_GL(B) = NSAV         !Number of SAV species in box B, what a bad name!!! should cal NSAV directly
                                  !
     READ(SVI,*) (SAVFRAC_GL(B,J),J=1,NSAVSPC_GL(B)) !Fraction of SAV coverage in the projected surface area of the
                                                     !bottom grid cell 
     READ(SVI,*) (SAVDPH_GL(B,J) ,J=1,NSAVSPC_GL(B)) !Depth of SAV (mean depth of the plant) for all 
                                                     !species in each grid cell (Wen Long: should replace with
                                                     !still water depth of the SAV cell as we no longer
                                                     !use fixed z grid that was in the old ICM)
   ENDDO

   IF(SERIAL) THEN
      NSAVCELL = NSAVCELL_GL
      NSAVSPC  =  NSAVSPC_GL
      SAVFRAC  =  SAVFRAC_GL
      SAVDPH   =   SAVDPH_GL
   ENDIF

# if defined (MULTIPROCESSOR)
   IF(PAR)THEN
         NSAVCELL=0
         DO I=1,NSAVCELL_GL
            !find global node
             B=SAVCELL_GL(I)          !global TCE (node) ID       
             J=NLID(B)                !local TCE (node) ID
             IF ( J .NE. 0 ) THEN     !find if B is in local domain
                !accumulate number of local NSAVCELL
                 NSAVCELL          =  NSAVCELL+1 
                 !local TCE ID  
                 SAVCELL(NSAVCELL) =  J
                !set the local SAVFRAC and SAVDPH from the global value
                 NSAVSPC(J)  =  NSAVSPC_GL(B) !number of species
                 SAVFRAC(J,1:NSAVSPC(J)) =  SAVFRAC_GL(B,1:NSAVSPC_GL(B))  !fraction of coverage in TCE J
                 SAVDPH(J,1:NSAVSPC(J)) =   SAVDPH_GL(B,1:NSAVSPC_GL(B))  !depth of SAV in TCE J
             ENDIF 
         ENDDO
         
         !Wen Long:
         !Should we deal with halo nodes using NGLD_X ?
         !I guess not
         
   ENDIF   
# endif
   
!***** Input FORMAT statements

1000 FORMAT(///)
1006 FORMAT(/)
1010 FORMAT(A72)
1020 FORMAT(//(8X,9I8))
1030 FORMAT(//(8X,9F8.0))
1032 FORMAT(//(8X,5F8.0))
1040 FORMAT(//8X,2A8)
1042 FORMAT(//8X,I8,A8)
1050 FORMAT(//(8X,A8,F8.0))

! OUTPUT WHAT WAS INPUT

   IF(MSR) THEN
  
     OPEN (SVO,FILE=SVOFN)
     WRITE(SVO,2010) (TITLE(J),J=1,6)

! WRITE SAV SPECIES NAME AND PROPERTIES

     WRITE(SVO,2040) SPNAM
     WRITE(SVO,2050) PMSAV,TRPMSAV,KTPS1,KTPS2
     WRITE(SVO,2052) PMSAV/ACDWSAV,TRPMSAV
     WRITE(SVO,2056) ACDWSAV,ANDWSAV,APDWSAV
     WRITE(SVO,4010) KHNLEAF,KHNPSAV,KHNROOT,KHPLEAF,KHPROOT
     WRITE(SVO,2060) BMSAV,TRBMSAV,KTBMSAV
     WRITE(SVO,2058) PRSPSAVB,SALMAX
     WRITE(SVO,3070) FDOSR
     WRITE(SVO,2090) KHSTOX
     WRITE(SVO,3000) ALPHSAV
     WRITE(SVO,3010) KESAV,ALAC,ACAN,BCAN
     IF (PRINTIT == '     ALL') THEN
       WRITE(SVO,3080)
       WRITE(SVO,3090) (J,FPLEAF(J),FPSTEM(J),FPROOT(J),              &
                        FPTUBER(J),TRTBRLF(J),J=1,NDDI)
     ENDIF

! WRITE DISTRIBUTION OF BIOMASS UPON MORTALITY

     WRITE(SVO,4020)
     WRITE(SVO,4030) FNISAV,FNLDSAV,FNRDSAV,FNLPSAV,FNRPSAV
     WRITE(SVO,4040) FPISAV,FPLDSAV,FPRDSAV,FPLPSAV,FPRPSAV
     WRITE(SVO,4050) FDOSAV,FCLDSAV,FCRDSAV,FCLPSAV,FCRPSAV
     WRITE(SVO,4080)
     WRITE(SVO,4090) FRPONSAV(1),FRPONSAV(2),FRPONSAV(3)
     WRITE(SVO,5000) FRPOPSAV(1),FRPOPSAV(2),FRPOPSAV(3)
     WRITE(SVO,5010) FRPOCSAV(1),FRPOCSAV(2),FRPOCSAV(3)

! WRITE EFFECT OF SAV ON NET SETTLING TO SEDIMENTS
     WRITE(SVO,5060) WSSSAV,WSLSAV,WSRSAV,WS1SAV,WS2SAV,WS3SAV,WSUSAV      

! WRITE EPIPHYTE PROPERTIES

     WRITE(SVO,2040) EPINAM
     WRITE(SVO,2051) PMEPI,TRPMEPI,KTPE1,KTPE2
     WRITE(SVO,2052) PMEPI/CCHLEPI,TRPMEPI
     WRITE(SVO,2062) CCHLEPI,ADWCEPI,ANCEPI,APCEPI
     WRITE(SVO,5020) KHNEPI,KHNPEPI,KHPEPI
     WRITE(SVO,2065) BMEPI,TRBMEPI,KTBME
     WRITE(SVO,3020) PREPI,TRPREPI,KTPRE
     WRITE(SVO,3002) ALPHEPI
     WRITE(SVO,3030) KEEPI
     WRITE(SVO,4000) KHEP

! WRITE DISTRIBUTION OF BIOMASS UPON MORTALITY

     WRITE(SVO,5030)
     WRITE(SVO,4030) FNIEPI,FNLDEPI,FNRDEPI,FNLPEPI,FNRPEPI
     WRITE(SVO,4040) FPIEPI,FPLDEPI,FPRDEPI,FPLPEPI,FPRPEPI
     WRITE(SVO,5050) FCLDEPI,FCRDEPI
     WRITE(SVO,5040)
     WRITE(SVO,4030) FNIPEP,FNLDPEP,FNRDPEP,FNLPPEP,FNRPPEP
     WRITE(SVO,4040) FPIPEP,FPLDPEP,FPRDPEP,FPLPPEP,FPRPPEP
     WRITE(SVO,4070) FCLDPEP,FCRDPEP,FCLPPEP,FCRPPEP
      
! WRITE NUMBER OF SAV CELLS, SAV WIDTH AND DEPTH INCREMENTS

     WRITE(SVO,5070) NSAVCELL_GL, NSAVSPCM !Wen Long: note that  NSAVSPCM is never used
     IF (PRINTIT == '     ALL') THEN
       DO I=1,NSAVCELL_GL
         B = SAVCELL_GL(I)
         WRITE(SVO,5080) I, B, NSAVSPC_GL(B)
         WRITE(SVO,5082) (SAVFRAC_GL(B,J),J=1,NSAVSPC_GL(B))  !Wen Long       
         WRITE(SVO,5090) ( SAVDPH_GL(B,J),J=1,NSAVSPC_GL(B))
       ENDDO
     ENDIF

  ENDIF  !END of MSR block

  DEALLOCATE(SAVCELL_GL,  &
             NSAVSPC_GL,  &
             SAVFRAC_GL,  &
              SAVDPH_GL)
  
! COMPUTE FRACTIONAL AREA OF EACH SAV DEPTH INCREMENT

!Wen Long    !The whole purpose of SAVWDH is to calculate
!Wen Long    !fraction of SAV occupation in box B, hence I deprecated
!Wen Long    !deprecate SAVWDH and replaced with SAVFRAC directly read from input

!Wen Long    
!Wen Long     DO I=1,NSAVCELL
!Wen Long       B=SAVCELL(I)                 !Box number of the I'th SAV cell
!Wen Long       DO N=1,NSAVSPC(B)            !Fraction of area SAV occupied
!Wen Long                                    !for each vertical increment step 
!Wen Long
!Wen Long deprecated SAVWDH.  Array SAVFRAC is now read in from input file directly
!Wen Long
!Wen Long !JQI         SAVFRAC(B,N)=SAVWDH(B,N)/BL(B,2)  !Wen Long: BL(B,2) is bottom box's y dimension length (horizontal)
!Wen Long                                                !was used in original structured grid version of ICM
!Wen Long                                                !to indicate extent of SAV growth away from shoreline
!Wen Long                                                !within alongshore grids and when the grids are too large
!Wen Long                                                !to resolve area of SAV beds
!Wen Long       
!Wen Long       
!Wen Long       ENDDO
!Wen Long     ENDDO
!Wen Long

! CREATE LOOKUP TABLE OF TEMPERATURE EFFECTS

     DO J = -50,400
       TLOOK = FLOAT(J)/10.

       !temperature control on SAV growth rate
       IF (TLOOK < TRPMSAV) THEN
         FTPSAV(J) = EXP(-KTPS1*(TLOOK-TRPMSAV)**2)
       ELSE
         FTPSAV(J) = EXP(-KTPS2*(TLOOK-TRPMSAV)**2)
       ENDIF

       !temperature control on SAV respiration
       FTRSAV(J) = EXP(KTBMSAV*(TLOOK-TRBMSAV))

     ENDDO

     DO J = -50,400
       TLOOK = FLOAT(J)/10.  !TLOOK ranges from -5 to 40 degree C

       !temperature control on epiphyte grwth
       IF (TLOOK < TRPMEPI) THEN
         FTPEP(J) = EXP(-KTPE1*(TLOOK-TRPMEPI)**2)
       ELSE
         FTPEP(J) = EXP(-KTPE2*(TLOOK-TRPMEPI)**2)
       ENDIF
       
       !temperature control on epiphytes respiration
       FTREP(J)  = EXP(KTBME*(TLOOK-TRBMEPI))  
       
       !temperature control on predation on epiphytes
       FTPREP(J) = EXP(KTPRE*(TLOOK-TRPREPI))  

     ENDDO
      
! INITIALIZE DEGREE DAYS TO ZERO

     DO I=1,MLOC
       DGRDAYS(I)=0
     ENDDO
      
! ZERO OUT ALL COUNTERS FOR EFECTS ON WATER COLUMN AND SEDIMENTS

     DO I=1,MLOC

       DOSAVW(I)   = 0.
       LDOCSAVW(I) = 0.
       RDOCSAVW(I) = 0.
       LPOCSAVW(I) = 0.
       RPOCSAVW(I) = 0.

       DOEPIW(I)   = 0.
       LDOCEPIW(I) = 0.
       RDOCEPIW(I) = 0.
       LPOCEPIW(I) = 0.
       RPOCEPIW(I) = 0.

       NH4SAVW(I)  = 0.
       NO3SAVW(I)  = 0.
       LDONSAVW(I) = 0.
       RDONSAVW(I) = 0.
       LPONSAVW(I) = 0.
       RPONSAVW(I) = 0.

       NH4EPIW(I)  = 0.
       NO3EPIW(I)  = 0.
       RDONEPIW(I) = 0.
       LDONEPIW(I) = 0.
       LPONEPIW(I) = 0.
       RPONEPIW(I) = 0.

       PO4SAVW(I)  = 0.
       LDOPSAVW(I) = 0.
       RDOPSAVW(I) = 0.
       LPOPSAVW(I) = 0.
       RPOPSAVW(I) = 0.
        
       PO4EPIW(I)  = 0.
       LDOPEPIW(I) = 0.
       RDOPEPIW(I) = 0.
       LPOPEPIW(I) = 0.
       RPOPEPIW(I) = 0.
        
     ENDDO

     DO I=1,MLOC

       SEDDOSAV(I) = 0.      !Dissolved oxygen flux leaving sediments due to SAV and epiphytes
       SEDNH4SAV(I) = 0.     !NH4 flux leaving sediments due to SAV and epiphytes
       SEDNO3SAV(I) = 0.     !NO3 flux (Wen Long: added NO3 here) leaving sediments due to SAV and epiphytes      
       SEDPO4SAV(I) = 0.     !PO4 flux leaving sediments due to SAV

       SEDPOCSAV(I)  = 0.    !POC flux into sediments due to SAV and epiphytes
       SEDPONSAV(I)   = 0.   !PON flux into sediments due to SAV and epiphytes
       SEDPOPSAV(I)   = 0.   !POP flux into sediments due to SAV and epiphytes

     ENDDO
    
!***** Output FORMAT statements

2010 FORMAT(1X,A72)
2040 FORMAT(/1X,A8)
2050 FORMAT(' PRODUCTION = ',F8.3,' GM C/GM DW/DAY AT ',F8.2,' C.'/         &
            ' KT1 = ',F8.3,' KT2 = ',F8.3,' PER DEGREE**2')
2051 FORMAT(' PRODUCTION = ',F8.3,' GM C/GM CHL/DAY AT ',F8.2,' C.'/        &
            ' KT1 = ',F8.3,' KT2 = ',F8.3,' PER DEGREE**2')
2052 FORMAT(' CARBON SPECIFIC GROWTH RATE = ',F8.3,' PER DAY AT ',F8.2,' C.')
2056 FORMAT(' CARBON TO DRY WEIGHT RATIO = ',F8.3/                          &
            ' NITROGEN TO DRY WEIGHT RATIO = ',F8.3/                        &
            ' PHOSPHORUS TO DRY WEIGHT RATIO = ',F8.3)
2058 FORMAT(' PHOTORESPIRATION = ',F8.3,' * PRODUCTION @ ZERO SALT'/        &
            ' PHOTORESPIRATION = 1.0 WHEN SALINITY > ',F8.3)     
2060 FORMAT(/' PLANT RESPIRATION = ',F8.3,' PER DAY AT ',F8.2,              &
            ' C.  KT = ',F8.3,' PER DEGREE')
2062 FORMAT(' CARBON TO CHLOROPHYLL RATIO = ',F8.3/                         &
            ' DRY WEIGHT TO CARBON RATIO = ',F8.3/                          &
            ' NITROGEN TO CARBON RATIO = ',F8.3/                            &
            ' PHOSPHORUS TO CARBON RATIO = ',F8.3)
2065 FORMAT(/' RESPIRATION = ',F8.3,' PER DAY AT ',F8.2,                    &
            ' C.  KT = ',F8.3,' PER DEGREE')
2090 FORMAT(' SULFIDE TOXICITY EFFECT = ',F8.3,' PER MOL HS')
3000 FORMAT(/' ALPHA = ',F8.3,' (GM C/GM DW) / (E/M**2)')
3002 FORMAT(/' ALPHA = ',F8.3,' (GM C/GM CHL) / (E/M**2)')
3010 FORMAT(' LIGHT ATTENUATION = ',F8.3,' M**2/GM SHOOT C',/               &
            ' M**2 LEAF AREA PER GM LEAF CARBON  = ',F8.3,/                 &
            ' CANOPY HEIGHT = ',F8.3,' + ',F8.3,' * SHOOTS (M)')
3020 FORMAT(' PREDATION = ',F8.3,' PER DAY AT ',F8.2,                       &
            ' C.  KT = ',F8.3,' PER DEGREE')
3030 FORMAT(' LIGHT ATTENUATION = ',F8.3,' M**2 LEAF AREA PER GM EPIPHYTE C'/)
3040 FORMAT(' SPECIES SPECIFIED INCORRECTLY IN CELL ',I5,' AS ',A8)
3050 FORMAT(/3(' CELL SPECIES CODE BLSAV'))
3060 FORMAT(3(I5,1X,A8,I4,F8.2))
3070 FORMAT(' FRACTION DO TRANSFERRED FROM SHOOTS TO ROOTS =',F8.3)
3080 FORMAT(/'  DGR DAY  FPLEAF   FPSTEM   FPROOT  FPTUBER TRTBRLF'/        &
             '           0<  <1   0<  <1   0<  <1  0<  <1  1/DAY')
3090 FORMAT(I7,3X,5F10.3)
4000 FORMAT(/' EPIPHYTE DENSITY AT WHICH GROWTH IS HALVED = ',F8.3,' GM C/GM C')
4010 FORMAT(/' NITROGEN HALF-SATURATION CONC (SHOOTS) = ',F8.3,             &
            ' GM/M**3'/                                                     &
            ' NH4 HALF-SATURATION CONC (PREFERENCE) = ',F8.3,' GM/M**3'/    &
            ' NITROGEN HALF-SATURATION CONC (ROOTS) = ',F8.3,' GM/M**3'/    &
            ' PHOSPHORUS HALF-SATURATION CONC (SHOOTS) = ',F8.3,            &
            ' GM/M**3'/                                                     &
            ' PHOSPHORUS HALF-SATURATION CONC (ROOTS) = ',F8.3,             &
            ' GM/M**3')
4020 FORMAT(/' DISTRIBUTION OF SAV UPON MORTALITY'/                         &
            ' SHOOT RESPIRATION  DIS INORG  LAB DISS  REF DISS  ',          &
            'LAB PART  REF PART')
4030 FORMAT(' NITROGEN          ',5F10.3)
4040 FORMAT(' PHOSPHORUS        ',5F10.3)
4050 FORMAT(' CARBON            ',5F10.3)
4070 FORMAT(' CARBON            ',10X,4F10.3)
4080 FORMAT(/' ROOT MORTALITY   LAB PART  REF PART     INERT')
4090 FORMAT(' NITROGEN       ',3F10.3)
5000 FORMAT(' PHOSPHORUS     ',3F10.3)
5010 FORMAT(' CARBON         ',3F10.3)
5020 FORMAT(/' NITROGEN HALF-SATURATION CONC = ',F8.3,' GM/M**3'/           &
            ' NH4 HALF-SATURATION CONC (PREFERENCE) = ',F8.3,' GM/M**3'/    &
            ' PHOSPHORUS HALF-SATURATION CONC = ',F8.3,' GM/M**3')
5030 FORMAT(/' DISTRIBUTION OF EPIPHYTES UPON MORTALITY'/                   &
            ' RESPIRATION        DIS INORG  LAB DISS  REF DISS  ',          &
            'LAB PART  REF PART')
5040 FORMAT(' PREDATION')
5050 FORMAT(' CARBON            ',10X,2F10.3)
5060 FORMAT(/'EFFECT OF SAV ON NET SETTLING TO SEDIMENTS ',                 &
            '(M**3/GM C/DAY)'/'  SOLIDS  LABILE REFRACT      B1      B2',   &
            '      B3     PBS'/7F8.2)
5070 FORMAT(/' NUMBER OF SAV CELLS = ',I8,/' MAXIMUM NUMBER OF DEPTHS = ',I8)
!5080 FORMAT(/' SAV CELL ',I6,', WQM CELL ',I6, ', # DEPTH INCR ',I6)
5080 FORMAT(/' SAV CELL ',I6,', WQM CELL ',I6, ', # Species ',I6)
5082 FORMAT(' Coverage Fraction ',10F8.1)  !Wen Long changed from Width (m) to Coverage (dimensionless)
5090 FORMAT(' DEPTH (M) ',10F8.2)
 
     RETURN
     END SUBROUTINE SAV_READ

!************************************************************************
!**                    S U B R O U T I N E   S A V                     **
!************************************************************************

   SUBROUTINE SAV_COMP
   USE MOD_LIMS, ONLY: MLOC, KBM1
   USE MOD_WQM, ONLY:  			&
                AOCR,           &!
                B,              &!
                DLT,            &!
                I0,             &!
                JDAY,           &!
                KADPO4,         &!
                NH4,            &!
                NO3,            &!
                PO4,            &!
                SALT,           &!
                SREFLECT,       &!
                SSI,            &!
                T!,             &!

   USE MOD_OWQ, ONLY : KDLAMB_int
   USE MOD_SED_SAV_EXCHANGE_VARS, ONLY: NH4T2TM1S_SHARE,PO4T2TM1S_SHARE,HST2TM1S_SHARE,M2_SHARE,PIE2HS_SHARE
   IMPLICIT NONE
   INTEGER :: IDGD(MLOC)        !Degree day (basically infering day of year from temperature)
   
   REAL(SP) :: DENLIM(MLOC,NSAVM),    &
           PRSPSAV(MLOC)        !Photochemical respiration and salinity toxicity to salinity 

   REAL(SP) :: IK,       &    !half constant in light limitation  equation (E/m^2/day)
           IATCNPY,  &    !light level at canopy top (E/m^2/day)
           IAVEP,    &    !light level at SAV epiphytes (E/m^2/day)        
           IAVSH,    &    !light level at SAV leaf (E/m^2/day) (light reaches epiphytes before SAV leaf)
           NSED,     &    !NH4 concentration in sediment pore water (anaerobic layer) (gN/m^3)
           NWAT,     &    !NH4 and NO3 concentration in water column (gN/m^3) exsposed to SAV
           PSED,     &    !PO4 concentration in sediment pore water (gP/m^3)
           PWAT,     &    !PO4 concentration in water column (gP/m^3) exposed to SAV
           LFOLD          !SAV leaf biomass (gC/m^2) of previous time step
           
   REAL(SP) :: NPRSAV(MLOC),   & !NH4 preference on nitrogen uptake by SAV (dimensionless)
           NPREPI(MLOC)      !NH4 preference on nitrogen uptake by Epiphytes (dimensionless)
!vjp 9/20/05 added declarations
   INTEGER :: ITEMP, I, J, N
   REAL(SP)    :: ATN,     &
              DF,      & !fraction of dissolved PO4 concentration relative to total PO4 cocentration in water column
              HSMOL,   & !HS (sulfide) concentration in sediment pore water (mole-S2/m^3)
!             STEST,   & !salinity used in calculating salt toxicity (Wen Long changed this to STOXI below) (psu)
              STOXI,   & !salinity used in calculating salt toxicity (psu)
              GLIMSAV, & !overall (global) limitation to SAV growth (N, P, light and sulfide HS toxicity) 
              GLIMEPI, & !overall (global) limitation to SAV growth (N, P, light)
              DELLF      !delta leaf, i.e. leaf biomass change (gC/m^2/day) 

   REAL(SP)    :: DLTDY,  &  !time step in days (day) for updating SAV and epiphyte biomass 
              TDIFF,  &  !temerature - 3 (degC), used to calculate DGRDAYS and then IDGD for discretizing
                         !leaf,stem,root,tuber share of the SAV growth etc as a function of temperature to
                         !account for change of functions due to aging of SAV (senescence)
              ZMEAN,  &  !SAVDPH (vertical coordinate of SAV distribution discretized vertically) (distance from surface) (m)
              ZTOCNPY    !Ztc in Cerco and Moore 2001 paper, depth (m) at the top of the canopy (distance from surface to 
                         !canopy top

   INTEGER:: KWC         !water column layer index    
   
! CONVERT MODEL TIME STEP TO DAYS, DETERMINE JULIAN DAY

   DLTDY = DLT/86400.
   J     = 1.0 + AMOD(JDAY,365.25)
      
! COMPUTE DEGREE DAYS FOR EACH CELL, DISCRETE ARRAY ELEMENT FOR
! TRANSFERS ETC

   DO I=1,NSAVCELL
     B = SAVCELL(I)
     TDIFF = MAX(T(B,KBM1)-3.,0.)              !Wen Long changed to KBM1
                                               !when T< 3, TDIFF is bounded from below by 0
                                               !so that when T changes over time, DGRDAYS keeps increasing
                                               !
     DGRDAYS(B) = DGRDAYS(B)+TDIFF*DLTDY
     IDGD(B) = DGRDAYS(B)/10.                  !Index in FPLEAF, FPSTEM, FPROOT, FPTUBER, TRTBRLF
                                               !which are functions of DGRDAYS/10
                                               !Here DGRDAYS mean integration of tempreature in time
                                               ! i.e. the accumulative history of temperature that the 
                                               ! plant has experienced and endured/enjoyed, which units [degC day]
                                               !Wen Long: I do not like the hard-wired 1/10 here.
                                               !
                                               !The look up table FPLEAF, FPSTEM, FPROOT, FPTUBER, TRTBRLF
                                               !has to be large enough to contain the whole life history
                                               !in [degC*day] units of the plant
                                               !when DGRDAYS --> infinity as time increases, the look up table
                                               !should asymptotically approach a set of limiting values corresponding
                                               !to the dying phase of the plant.
                                               !
                                               !Wen Long: what about the regeneration (newly born) of plants?
                                               ! When should we introduce new generations ? and how to keep 
                                               ! different SAV's with in different aging phase simutaneously existing
                                               ! in the model?
     
                                               !Wen Long: We ought to keep a profile distribution of plants that are in different
                                               !life stages at the same time, some of them really yong, some of them  
                                               !being adult and reproductive, some of them being very old and not as 
                                               !energetic and active.  We also ought to keep track of 
                                               !several species. And then we have to keep track of their height of above
                                               !ground part, and their volume of occupancty of their roots (below ground). 
                                               !And we ought to sort the height of all classes(species, life stage) 
                                               !to determine the amount of light that each class is receiving. 
                                               !The amount of nutrients being uptaken by each class is 
                                               !determined by its height or biomass of leaves (photosynthesis)
                                               !and its root biomass and root volume occupancy
                                               !when root density exceeds certain value, we shall assume roots grow
                                               !horizontally (rhizome growth) into neighboring territory. The above ground height is composed
                                               !by leaf, stem and tuber. Reference: Keiko Aioi, Seasonal Change in the
                                               !Standing Crop of Eelgrass (Zostera Marina L. ) in Odawa Bay, Central Japan,
                                               !Aquatic Botany, Vol. 8, pp 343-354, 1980. Orth and Moore, Seasonal and 
                                               !year to year variations in the growth of zosteria marina L. in the lower
                                               !Chesapeake Bay, Aquatic Botany, 24, pp. 335-341, 1986. 
                                               !
                                               !Laugier T. et. al., Seasonal dynamics in mixed eelgrass beds,
                                               !Zostera marina L. and Z. noltii Hornem.,in a Mediterranean coastal lagoon
                                               !(Thau lagoon, France) Aquatic Botany, 63, pp 51-69, 1999 gives log distribution
                                               !of mortality as function of shoot and root density
                                               !
                                               !Shoot recruitment should also be a function of age and a parameter
                                               !to the model.  see Laugier paper Table 5 above. 
                                               !Shoot turn over (1/shoot-life [1/year]) should also be modeled
                                               !separately from stem and root
                                               !
                                               !Tuber is an interesting component, it stores nutrients for next
                                               !year or when nutrients are depleted. So its function should 
                                               !aslo be modeled differently.  See reference : Harwell and Orth, Influence
                                               ! of a tube-dwelling polychaete on the dispersal of fragmented 
                                               ! reproductive shoots of eelgrass, Aquatic Botany, 70, pp. 1-7, 2001
                                               !
                                               !When rhizomes become over crowded, the plant stops flowering too.
                                               !
                                               !Alternatively, we model mortality difference through size of patches
                                               !
     
   ENDDO

! ZERO OUT AFFECTS OF SAV CELLS ON WATER COLUMN AND SEDIMENTS

   DO I=1,NSAVCELL
     B = SAVCELL(I)      
     DOSAVW(B)   = 0.
     LDOCSAVW(B) = 0.
     RDOCSAVW(B) = 0.
     LPOCSAVW(B) = 0.
     RPOCSAVW(B) = 0.

     DOEPIW(B)   = 0.
     LDOCEPIW(B) = 0.
     RDOCEPIW(B) = 0.
     LPOCEPIW(B) = 0.
     RPOCEPIW(B) = 0.

     NH4SAVW(B)  = 0.
     NO3SAVW(B)  = 0.
     LDONSAVW(B) = 0.
     RDONSAVW(B) = 0.
     LPONSAVW(B) = 0.
     RPONSAVW(B) = 0.

     NH4EPIW(B)  = 0.
     NO3EPIW(B)  = 0.
     RDONEPIW(B) = 0.
     LDONEPIW(B) = 0.
     LPONEPIW(B) = 0.
     RPONEPIW(B) = 0.

     PO4SAVW(B)  = 0.
     LDOPSAVW(B) = 0.
     RDOPSAVW(B) = 0.
     LPOPSAVW(B) = 0.
     RPOPSAVW(B) = 0.
        
     PO4EPIW(B)  = 0.
     LDOPEPIW(B) = 0.
     RDOPEPIW(B) = 0.
     LPOPEPIW(B) = 0.
     RPOPEPIW(B) = 0.
        
   ENDDO

   DO I=1,NSAVCELL
     B = SAVCELL(I)      
     SEDDOSAV(B)    = 0.
     SEDPOCSAV(B)   = 0. 
     SEDNH4SAV(B)   = 0.
     SEDPO4SAV(B)   = 0.
     SEDPONSAV(B)   = 0.
     SEDPOPSAV(B)   = 0.
   ENDDO

! LIGHT EFFECTS

   DO I=1,NSAVCELL
     B = SAVCELL(I)

     ITEMP = 10. * T(B,KBM1) + 0.05   !Here T should be T(B,KBM1), i.e. bottom temperature   
     
     DO N=1,NSAVSPC(B)                !What is N here? species index
        
! COMPUTE IRRADIANCE AT CANOPY HEIGHT 
!
!    Titus and Adams 1979, "Coexistence and comparative light relations of 
!          the submersed macrophytes Myriophyllum spicatum L. and Vallisneria americana Michx",
!          Oecologia (Berl.) 40, 273-286)
!    Titus, J.E., Goldstein, R.A., Adams, M.S., Mankin, J.B., O'Neill, R.V., Weiler, P.R. Jr., Shugart,
!           H.H., Booth, R.S. : A production model for Myriophyllum spicatum L. Ecology 56, 1129-1138

       ZMEAN   = SAVDPH(B,N)                        !depth of SAV veritcal increment (mean depth level of canopy (m from surface)?)
       HCAN = ACAN + BCAN*(LEAF(B,N)+STEM(B,N))     !canopy height increases (m) with leaf and stem growth
       HCAN = MIN(ZMEAN,HCAN)                       !make sure ZMEAN > HCAN              
                                                    !so that ZTOCNPY will be positive below
       ZTOCNPY = ZMEAN-HCAN                         !depth level of canopy top (m from surface), always positive
!       WATATN(B,N) = EXP(-KESS(B,1)*ZTOCNPY)       !attenuation due to water column 
                                                    !Wen Long: KESS(B) should now be KESS(B,k), where k is
                                                    !layer index, and B is bottom layer
                                                    !also we should do it from surface of water
                                                    !repeatdly until reaching ZTOCNPY depth to get IATCNPY below
       
       WATATN(B,N) = EXP(-KDLAMB_int(B,KBM1)*ZTOCNPY)     !attenuation due to water column 
       
       IATCNPY = I0*(1.-SREFLECT/100.)*WATATN(B,N)  !Ic in Cerco and Moore 2001, light at canopy top

! COMPUTE ATTENUATION BY EPIPHYTES AND SELF SHADING

! EPATN(B,N) = EXP(-KEEPI*ACLA*ADWCEPI*EP(B,N))  !Old version with ACLA= gC-leaf/m^2-leaf, EP=gC-epi/gC-leaf      
                                                     !Cerco & Moore,
                                                     !Estuaries Vol. 24, No. 4, p. 522-534 August 2001
                                                     !equation for Ish
                                                     !KEEPI unit:           m^2 leaf/gDW-epi 
                                                     !ADWCEPI unit:         gDW-epi/gC-epi
                                                     !EP :                   gC-epi/m^2-leaf
                                                     !---------------------------------------------------
                                                     !KEEPI*ADWCEP*EP :      1    (nondimensional)  !check ! 
 
       EPATN(B,N) = EXP(-KEEPI*ADWCEPI*EP(B,N))      !New version now EP is in unit gC-epi/m^2-leaf

       ATN = KDLAMB_int(B,KBM1)*HCAN+KESAV*(LEAF(B,N)+STEM(B,N))
       SAVATN(B,N) = (1.-EXP(-ATN))/(ATN+1.0E-6)  

       !mean light within canopy
       IAVEP = IATCNPY*SAVATN(B,N)                !Iwc in Cerco and Moore 2001 (mean light within canopy)

       !mean light on shoots 
       IAVSH = IAVEP*EPATN(B,N)                   !Ish in Cerco and Moore 2001

! COMPUTE LIGHT LIMITATIONS TO GROWTH

       IK = PMSAV*FTPSAV(ITEMP)/ALPHSAV
       FISH(B,N) = IAVSH/(SQRT(IK*IK+IAVSH*IAVSH)+1.0E-30)
       IK = PMEPI*FTPEP(ITEMP)/ALPHEPI
       FIEP(B,N) = IAVEP/(SQRT(IK*IK+IAVEP*IAVEP)+1.0E-30)
          
     ENDDO

   ENDDO

! COMPUTE NUTRIENT LIMITATIONS TO GROWTH, SAV FIRST

   KWC=KBM1   !water column layer (bottom layer)   

   DO I=1,NSAVCELL
      B = SAVCELL(I)
    
      NSED = MAX(0.,NH4T2TM1S_SHARE(B)/1000.)    !NH4 in sediment anaerobic (2nd) layer  mgN/m^3 --> gN/m^3  ==  mgN/L
                                               !WL: why not also take NO3 in sediemnt top layer ?
      NWAT = MAX(0.,NH4(B,KWC)+NO3(B,KWC))         !NH4+NO3 in bottom layer of water column (gN/m^3)  (WAT means water)
      PSED = MAX(0.,PO4T2TM1S_SHARE(B)/1000.)    !PO4 in pore water of sediment anaerobic layer (gP/m^3)
                                               !WL: why not use pore water, the upper calculation is total concentration
                                               !
      DF = 1./(1.+KADPO4*SSI(B,KWC))       !fraction of dissolved PO4 in water column
                                           ! where KADPO4 is PO4 partitioning coefficient m^3/gSolids, as SSI has unit gSolids/m^3
      PWAT = MAX(0.,DF*PO4(B,KWC))         !Dissolved PO4 (gP/m^3), where PO4 is total concentration (gP/m^3)

      !Nitrogen limitation on SAV growth
      NLSAV(B) = (KHNROOT*NWAT+KHNLEAF*NSED)/                 &
                 (KHNLEAF*KHNROOT+(KHNROOT*NWAT+KHNLEAF*NSED))

      !fraction of nitrogen uptake by SAV due to inorganic nitrogen in sediment pore water
      FNSEDSAV(B) = (KHNLEAF*NSED)/                           &
                   (KHNROOT*NWAT+KHNLEAF*NSED+1.0E-6)

      !preference of NH4 over NO3 on N uptake by SAV, controled by threshould KHNPSAV (gN/m^3)
    
      NPRSAV(B) = NH4(B,KWC)/(KHNPSAV+NO3(B,KWC)+1.E-6)*              &   
                 (NO3(B,KWC)/(KHNPSAV+NH4(B,KWC)+1.E-6)               &
                 +KHNPSAV/(1.E-6+NH4(B,KWC)+NO3(B,KWC)))

      !Phosphorus limitation on SAV growth 
      PLSAV(B) = (KHPROOT*PWAT+KHPLEAF*PSED)/                 &
                 (KHPLEAF*KHPROOT+KHPROOT*PWAT+KHPLEAF*PSED)  

      !fraction of PO4 uptake by SAV due to PO4 in sediment pore water
      FPSEDSAV(B) = (KHPLEAF*PSED)/(KHPROOT*PWAT+KHPLEAF*PSED+1.0E-6)

! EPIPHYTES

      NLEPI(B) = NWAT/(KHNEPI+NWAT+1.0E-6)
      NPREPI(B) = NH4(B,KWC)/(KHNPEPI+NO3(B,KWC)+1.E-6)*              &
                 (NO3(B,KWC)/(KHNPEPI+NH4(B,KWC)+1.E-6)               &
                 +KHNPEPI/(1.E-6+NH4(B,KWC)+NO3(B,KWC)))    !Same as eqn (6-17) of 2002 Chesapeake Bay Model report
                                                            !"The 2002 Chesapeake Bay Eutrophication Model", 
                                                            !Cerco C.F. and Noel M. R., July 2004, EPA 903-R-04-004
     PLEPI(B) = PWAT/(KHPEPI+PWAT+1.0E-6)
      
! COMPUTE SULFIDE TOXICITY EFECT ON SAV, BASICALLY LIMITING GROWTH ON SAV
  
      HSMOL = HST2TM1S_SHARE(B)/(1.+M2_SHARE*PIE2HS_SHARE)/64000.    !HS concetration in sediment pore water
                                                  !where PIE2HS is Sulfide partitiong coefficient  (L/Kg)
                                                  !M2 is solids concentration (kg/L)
                                                  !HST2TM1S is total concentation of HS in sediments (mg/m^3)
                                                  !==> HSMOL is in units of mole-S2/m^3

      FHS(B) = EXP(-KHSTOX*HSMOL)                 !KHSTOX is of unit (m^3/mole-S2)

! COMPUTE EPIPHYTE DENSITY LIMITATION

      DO N=1,NSAVSPC(B)                   
         DENLIM(B,N) = KHEP/(KHEP+EP(B,N)+1.0E-6) !densitity limitation of species N
      ENDDO
    ENDDO

    KWC=KBM1
    
    DO I=1,NSAVCELL
       B = SAVCELL(I)

       ITEMP = 10. * T(B,KBM1) + 0.05  !WLong, changed to KBM1
        
! COMPUTE RESPIRATION

       !Photosynshetic respiration with toxicity to SAV enhanced by salinity      
       STOXI   = MAX(SALT(B,KWC),0.0)
       PRSPSAV(B) = PRSPSAVB + (1.-PRSPSAVB)   &
                              * 0.5 * (1.+TANH(STOXI-SALMAX)) !Transition from PRSPSAVB (value < 1) to 1 
                                                              !using a tanh function with transition
                                                              !centered at salinity = SALMAX
       BMLEAF(B)= BMSAV*FTRSAV(ITEMP)
       BMSTEM(B)= BMSAV*FTRSAV(ITEMP)
       BMROOT(B)= BMSAV*FTRSAV(ITEMP)
       BMTUBER(B)= BMTBRREF*FTRSAV(ITEMP)
       SLSH(B)  = BMSAV*FTRSAV(ITEMP)*(FCLPSAV+FCRPSAV)

! DETERMINE GROWTH FOR INDIVIDUAL DEPTH INCREMENTS

      DO N=1,NSAVSPC(B)  !Loop over all species

        GLIMSAV = MIN(NLSAV(B),PLSAV(B),FISH(B,N),FHS(B)) !global limitation due to N, P, light and HS toxicity
        GLIMEPI = MIN(NLEPI(B),PLEPI(B),FIEP(B,N))

      !leaf growth rate (1/day)
        PLEAF(B,N) = PMSAV              &  !gC/gDW/day
                    *FTPSAV(ITEMP)      &  !temperature control
                    *GLIMSAV            &  !light limitation
                    /ACDWSAV               !carbon to dry weight ratio gC/gDW
                                           !==> PLEAF has unit 1/day
     !net production of SAV (gC/m^2/day)       
        NPPSAV(B,N) = (PLEAF(B,N)-BMLEAF(B))*LEAF(B,N) !!net production of SAV (gC/m^2/day)     

! EPIPHYTE GROWTH AND RESPIRATION
! PHOTORESPIRATION IS NOT CONSISTENT WITH SAV
! FOR NOW LEAVE IT ALONE

      !epiphytes growth rate (1/day)
         PEP(B,N) = PMEPI             &  !gC/gCHLA/day
                *FTPEP(ITEMP)         &  !temperature control
                *GLIMEPI              &  !light   limitation
                *DENLIM(B,N)          &  !density limitation
                /CCHLEPI                 !gC/gCHLA
                                         !==> PEP has unit 1/day 
!WL   !respiration rate 1/day
       BMEP(B,N) = BMEPI*FTREP(ITEMP)+PRSPEPI*PEP(B,N)  !basal metabolism and photochemical respiration (1/day)
!WL   !predation 1/day
       PREP(B,N) = PREPI*FTPREP(ITEMP)*EP(B,N)    !predation rate m^2-leaf/gC-Epi/day * (gC-Epi/m^2-leaf) ==> (1/day)
!WL   !net production gC-Epi/gC-leaf/day (growth - respiration)
       NPPEPI(B,N) = (PEP(B,N)-BMEP(B,N))*EP(B,N) !net epiphytes growth [gC-epiphyte/m^2-leaf)/day] 

     ENDDO

   ENDDO
!
! INTERACTIONS WITH WATER COLUMN
!

   DO I=1,NSAVCELL
     B = SAVCELL(I)
     DO N=1,NSAVSPC(B)  !Loop over all species
      
      !DISSOLVED OXYGEN AND CARBON
      !DO flux due to SAV (gO2/m^2/day), positive generating DO in the bottom layer (SAV cell)

       DOSAVW(B) = DOSAVW(B)+SAVFRAC(B,N)*                         &
                 (AOCR*(1.-FDOSR)*PLEAF(B,N)*LEAF(B,N)             & !DO production discounted (reduceed) by FDOSR 
                 -AOCR*((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))          & !for leaf's production routed to root 
                 *LEAF(B,N)+BMSTEM(B)*STEM(B,N))*FDOSAV)             !Discount (reduce) by (1-FDOSAV) because
                                                                     !oxygen consumption of  STEM and LEAF are
                                                                     !routed to ROOT and Tuber ??
     
       LDOCSAVW(B) = LDOCSAVW(B)+SAVFRAC(B,N)*                     &
                  ((LEAF(B,N)*(BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))    &  
                  +STEM(B,N)*BMSTEM(B))*FCLDSAV)  !WLong: PLEAF has unit 1/day, PRSPSAV has no units
     
       RDOCSAVW(B) = RDOCSAVW(B)+SAVFRAC(B,N)*                     &
                   ((LEAF(B,N)*(BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))   &
                   +STEM(B,N)*BMSTEM(B))*FCRDSAV)
       LPOCSAVW(B) = LPOCSAVW(B)+SAVFRAC(B,N)*                     &
                   ((LEAF(B,N)*(BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))   &
                   +STEM(B,N)*BMSTEM(B))*FCLPSAV)
       RPOCSAVW(B) = RPOCSAVW(B)+SAVFRAC(B,N)*                     &
                   ((LEAF(B,N)*(BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))   &
                   +STEM(B,N)*BMSTEM(B))*FCRPSAV)

       !flux of oxygen to water column due to epiphytes growth (photosynthesis) - respiration consumption of O2
       !gO2/m^2/day (positive into water column)
       DOEPIW(B) = DOEPIW(B)               & 
                   + SAVFRAC(B,N)          &     !             1
                     *EP(B,N)              &     !       gC-Epi/m^2-leaf
                     *LEAF(B,N)            &     !      gC-leaf/m^2
                     *ALAC                 &     !     m^2-leaf/gC-leaf     !(2002 version: ACLA: gC-leaf/m^2-leaf)
                     *AOCR                 &     !          gO2/gC
                     *(PEP(B,N)-BMEP(B,N)  &     !            1/day
                     *(1.-FCLDEPI-FCRDEPI))      !             1
                                                 !---------------------------
                                                 !          gO2/m^2/day

       !flux to LDOC by epi mortality (gC/m^2/day), positive to water column
       LDOCEPIW(B) = LDOCEPIW(B)                                             &
                    +SAVFRAC(B,N)                                            &
                    *(  EP(B,N)*LEAF(B,N)*ALAC*BMEP(B,N)          *FCLDEPI   & !basal metabolism routed to LDOC
                      + EP(B,N)*LEAF(B,N)*ALAC*(PREP(B,N)+SLSH(B))*FCLDPEP   & !predation and sloughing routed to LDOC
                     )
                    
       !flux to RDOC by epi mortality (gC/m^2/day), positive to water column
       RDOCEPIW(B) = RDOCEPIW(B)                                             &
                    +SAVFRAC(B,N)                                            &
                    *(  EP(B,N)*LEAF(B,N)*ALAC*BMEP(B,N)          *FCRDEPI   & !basal metabolism routed to RDOC
                      + EP(B,N)*LEAF(B,N)*ALAC*(PREP(B,N)+SLSH(B))*FCRDPEP   & !predation and sloughing routed to RDOC
                     )
       !flux to LPOC by epi mortality (gC/m^2/day), positive to water column
       LPOCEPIW(B) = LPOCEPIW(B)                                             &
                    +SAVFRAC(B,N)                                            &
                    *(EP(B,N)*LEAF(B,N)*ALAC*(PREP(B,N)+SLSH(B))*FCLPPEP     &  !predation and sloughing
                     )  

       !flux to RPOC by epi mortality (gC/m^2/day), positive to water column 
       RPOCEPIW(B) = RPOCEPIW(B)                                             &
                    +SAVFRAC(B,N)                                            &
                    *( EP(B,N)*LEAF(B,N)*ALAC*(PREP(B,N)+SLSH(B))            &  !predation and sloughing
                      *(1.-FCLDPEP-FCRDPEP-FCLPPEP)                          &  !routed to RPOC
                     ) !Note that 1 - FCLDPEP - FCRDPEP - FCLPPEP <==> FCRPPEP

! NITROGEN AND PHOSPHORUS

!WLong
       NH4SAVW(B) = NH4SAVW(B)+SAVFRAC(B,N)*                     &
                   ( ANCSAV*FNISAV*(                             &
                                    (  BMLEAF(B)                 &   !basal metabolism of leaf
                                      +PLEAF(B,N)*PRSPSAV(B)     &   !mortality due to salt toxicity and photochemical respiration
                                     )*LEAF(B,N)                 &   
                                    +BMSTEM(B)*STEM(B,N)         &   !basal metabolism of stem
                                   )                             &
                     - (1.-FNSEDSAV(B))                          &
                       *NPRSAV(B)*PLEAF(B,N)*ANCSAV*LEAF(B,N)    &
                    )
     
       NO3SAVW(B) = NO3SAVW(B)-SAVFRAC(B,N)*                     &
                   (1.-FNSEDSAV(B))*(1.-NPRSAV(B))*PLEAF(B,N)*   &
                    ANCSAV*LEAF(B,N)
                    
       !
       !Wen Long: above we have NPRSAV used in calculating NH4SAVW and NO3SAVW
       !          I think it should also play a role in SEDNH4SAV and SEDNO3SAV
       !          It seems that we are missing the split to SEDNO3SAV
       !
       
       !LDON flux (gN/m^2/day)
       LDONSAVW(B) = LDONSAVW(B)+SAVFRAC(B,N)*FNLDSAV*ANCSAV     &
                    *((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))          &
                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N))
       !RDON flux (gN/m^2/day)
       RDONSAVW(B) = RDONSAVW(B)+SAVFRAC(B,N)*FNRDSAV*ANCSAV     &
                    *((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))          &
                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N))
                    
       !LPON flux (gN/m^2/day)
       LPONSAVW(B) = LPONSAVW(B)+SAVFRAC(B,N)*FNLPSAV*ANCSAV     &
                    *((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))          &
                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N))

       !RPON flux (gN/m^2/day)
       RPONSAVW(B) = RPONSAVW(B)+SAVFRAC(B,N)*FNRPSAV*ANCSAV     &
                    *((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))          &
                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N))

       !NH4 flux (gN/m^2/day)
       NH4EPIW(B)  = NH4EPIW(B)                            &
                      +SAVFRAC(B,N)*                       &
                       (  FNIEPI*BMEP(B,N)                 &  !Basal metabolism and photochemical respiration
                         +FNIPEP*(PREP(B,N)+SLSH(B))       &
                         -NPREPI(B)*PEP(B,N)               &
                       )                                   &
                       *ANCEPI*EP(B,N)*LEAF(B,N)*ALAC         

       !NO3 flux (gN/m^2/day)
       NO3EPIW(B) = NO3EPIW(B)-SAVFRAC(B,N)*(1.-NPREPI(B))*PEP(B,N) &
                   *ANCEPI*EP(B,N)*LEAF(B,N)*ALAC

       !LDON flux due to Epiphyte (gN/m^2/day)
       LDONEPIW(B) = LDONEPIW(B)+SAVFRAC(B,N)*                       &
                     (FNLDEPI*BMEP(B,N)+FNLDPEP*(PREP(B,N)+SLSH(B))) &
                    *ANCEPI*EP(B,N)*LEAF(B,N)*ALAC                  

       !RDON flux due to Epiphyte (gN/m^2/day)
       RDONEPIW(B) = RDONEPIW(B)+SAVFRAC(B,N)*                       &
                    (FNRDEPI*BMEP(B,N)+FNRDPEP*(PREP(B,N)+SLSH(B)))  &
                    *ANCEPI*EP(B,N)*LEAF(B,N)*ALAC

       !LPON flux (gN/m^2/day)
       LPONEPIW(B) = LPONEPIW(B)+SAVFRAC(B,N)*                       &
                    (FNLPEPI*BMEP(B,N)+FNLPPEP*(PREP(B,N)+SLSH(B)))  &
                   *ANCEPI*EP(B,N)*LEAF(B,N)*ALAC
                   
       !RPON flux (gN/m^2/day)
       RPONEPIW(B) = RPONEPIW(B)+SAVFRAC(B,N)*                       &
                    (FNRPEPI*BMEP(B,N)+FNRPPEP*(PREP(B,N)+SLSH(B)))  &
                    *ANCEPI*EP(B,N)*LEAF(B,N)*ALAC

      !P fluxes to water column (gP/m^2/day)

       !PO4 flux (gP/m^2/day)
       PO4SAVW(B) = PO4SAVW(B)+                                         &
                   SAVFRAC(B,N)*(APCSAV*FPISAV*                        &
                                   ((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))  &
                                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N)     &
                                   )                                   &
                                 -(1.-FPSEDSAV(B))                     &
                                  *PLEAF(B,N)*APCSAV*LEAF(B,N)         &
                                 )
    
       !LDOP flux (gP/m^2/day)
       LDOPSAVW(B) = LDOPSAVW(B)+SAVFRAC(B,N)*FPLDSAV*APCSAV        &
                    *((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))             &
                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N))

       !RDOP flux (gP/m^2/day)
       RDOPSAVW(B) = RDOPSAVW(B)+SAVFRAC(B,N)*FPRDSAV*APCSAV        &
                    *((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))             &
                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N))

       !LPOP flux (gP/m^2/day)
       LPOPSAVW(B) = LPOPSAVW(B)+SAVFRAC(B,N)*FPLPSAV*APCSAV        &
                    *((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))             &
                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N))
   
       !RPOP flux (gP/m^2/day)
       RPOPSAVW(B) = RPOPSAVW(B)+SAVFRAC(B,N)*FPRPSAV*APCSAV        &
                    *((BMLEAF(B)+PLEAF(B,N)*PRSPSAV(B))             &
                    *LEAF(B,N)+BMSTEM(B)*STEM(B,N))

       !PO4 flux due to Epiphyte
       PO4EPIW(B)  = PO4EPIW(B)+SAVFRAC(B,N)*                       &
                  (FPIEPI*BMEP(B,N)+FPIPEP*(PREP(B,N)+SLSH(B))      &
                  -PEP(B,N))*APCEPI*EP(B,N)*LEAF(B,N)*ALAC 
                  

       LDOPEPIW(B) = LDOPEPIW(B)+SAVFRAC(B,N)*                        &
                    (FPLDEPI*BMEP(B,N)+FPLDPEP*(PREP(B,N)+SLSH(B)))   &
                    *APCEPI*EP(B,N)*LEAF(B,N)*ALAC
       RDOPEPIW(B) = RDOPEPIW(B)+SAVFRAC(B,N)*                        &
                   (FPRDEPI*BMEP(B,N)+FPRDPEP*(PREP(B,N)+SLSH(B)))    &
                    *APCEPI*EP(B,N)*LEAF(B,N)*ALAC
       LPOPEPIW(B) = LPOPEPIW(B)+SAVFRAC(B,N)*                        &
                    (FPLPEPI*BMEP(B,N)+FPLPPEP*(PREP(B,N)+SLSH(B)))   &
                   *APCEPI*EP(B,N)*LEAF(B,N)*ALAC
       RPOPEPIW(B) = RPOPEPIW(B)+SAVFRAC(B,N)*                        &
                   (FPRPEPI*BMEP(B,N)+FPRPPEP*(PREP(B,N)+SLSH(B)))    &
                    *APCEPI*EP(B,N)*LEAF(B,N)*ALAC

     ENDDO
   ENDDO

! INTERACTIONS WITH SEDIMENTS

   DO I=1,NSAVCELL
     B = SAVCELL(I)
     DO N=1,NSAVSPC(B)  !loop over all species
      
! DISSOLVED OXYGEN AND CARBON
  
       SEDDOSAV(B) = SEDDOSAV(B)+SAVFRAC(B,N)*                    &    !Wen Long: why positive here? Should be negative
                    (LEAF(B,N)*AOCR*FDOSR*PLEAF(B,N)              &    !Leaf mortality related DO consumption
                                                                       !(oxidation of C to CO2)
                                                                       !
                    -( ROOT(B,N)*BMROOT(B)                        &    !Root respiration and tuber respiration
                      +TUBER(B,N)*BMTUBER(B)                      &    !cost of oxygen (mortality loss of 
                     )*AOCR*FDOSAV                                &    !Root and Tuber carbon converted to oxygen
                    )                                                  !consumption and oxydation of C to CO2 )
                                                                       !FDOSAV is the portion of carbon loss
                                                                       !assocatied with oxydation of Organic Carbon 
                                                                       !to CO2 (DIC) 

		!WLong : we might need to split the above into DO production in water column
		!due to leaf growth and Sediment Oxygen Demand due to root respiration
					
        !SEDDOSAV should be added to SOD directly and let it affect sediment
        !structure (such as first layer thickness)
        						
                                     
       SEDPOCSAV(B) = SEDPOCSAV(B)+SAVFRAC(B,N)*                  &    !the rest part (1-FDOSAV) that remains to be 
                   (ROOT(B,N)*BMROOT(B)+TUBER(B,N)*BMTUBER(B))*   &    !POC and settles to sediment up on 
                   (1.-FDOSAV)                                         !Root and Tuber mortality

! NITROGEN AND PHOSPHORUS

!WL       SEDNH4SAV(B) = SEDNH4SAV(B)+SAVFRAC(B,N)*                    &
!WL                      FNSEDSAV(B)*PLEAF(B,N)*ANCSAV*LEAF(B,N)

      !                   
      !Wen Long: what about NO3 flux with sediments due to SAV uptake?
      ! Here I split it using NPRSAV(B) (preference on NH4 over NO3)
                   
       SEDNH4SAV(B) = SEDNH4SAV(B)+SAVFRAC(B,N)*                    &
                      FNSEDSAV(B)*NPRSAV(B)*PLEAF(B,N)*ANCSAV*LEAF(B,N)

       SEDNO3SAV(B) = SEDNO3SAV(B)+SAVFRAC(B,N)*                    &
                      FNSEDSAV(B)*(1.-NPRSAV(B))*PLEAF(B,N)*ANCSAV*LEAF(B,N)
                      
       SEDPO4SAV(B) = SEDPO4SAV(B)+SAVFRAC(B,N)*                    &
                      FPSEDSAV(B)*PLEAF(B,N)*APCSAV*LEAF(B,N)

       SEDPONSAV(B) = SEDPONSAV(B)+SAVFRAC(B,N)*                    &
                     (ROOT(B,N)*BMROOT(B)+TUBER(B,N)*BMTUBER(B))*ANCSAV

       SEDPOPSAV(B) = SEDPOPSAV(B)+SAVFRAC(B,N)*                    &
                     (ROOT(B,N)*BMROOT(B)+TUBER(B,N)*BMTUBER(B))*APCSAV

     ENDDO
   ENDDO

! COMPUTE NEW SHOOT (leaf, tube, stem), ROOT, AND EPIPHYTE BIOMASS

   DO I=1,NSAVCELL
     B = SAVCELL(I)

! DETERMINE DEGREE DAY INTERVAL.  SET MAXIMUM TO OCCURANCE OF SENESCENCE
! TO PREVENT ARRAY OVERFLOW      
     J = MIN(IDGD(B),NDDI-1)   !entry to the loop up table in FPROOT, FPLEAF, etc for senescence of SAV (aging)
          
     DO N=1,NSAVSPC(B) !loop over all species

       LFOLD = LEAF(B,N)  !previous leaf biomass gC/m^2

       !leaf biomass change rate(gC/m^2/day)     
       DELLF = (  PLEAF(B,N)                      & !leaf growth 1/day
              *(1.-PRSPSAV(B))                    & !discount mortality due to salt toxicity 
                                                    !when salt> SALTMAX, PRSPSAV --> 1 and there will be no growth
              *FPLEAF(J)                          & !
              -BMLEAF(B)                          & !discount basal metabolism
             )*LFOLD                              & 
            + TRTBRLF(J)*TUBER(B,N)                 !transfer from tuber to leaf

       !update root biomass (gC/m^2)
       ROOT(B,N) = ROOT(B,N) +                    &
                     (   FPROOT(J)                &
                        *PLEAF(B,N)               & !leaf growth 1/day
                        *(1.-PRSPSAV(B))          & !discount mortality due to salt toxicity
                        *LFOLD                    &
                       - BMROOT(B)*ROOT(B,N)      & !discount basal metabolism of root
                     )*DLTDY

       !update stem biomass (gC/m^2)
       STEM(B,N) = STEM(B,N) +                    & 
               (   FPSTEM(J)                      &
                  *PLEAF(B,N)                     &
                     *(1.-PRSPSAV(B))             &
                  *LFOLD                          &
                 - BMSTEM(B)*STEM(B,N)            &
                )*DLTDY
     
       TUBER(B,N) = TUBER(B,N) +                  &
                  (  FPTUBER(J)                   &
                    *PLEAF(B,N)                   &
                    *(1.-PRSPSAV(B))              &
                    *LFOLD                        &
                      -(BMTUBER(B)+TRTBRLF(J))    &
                    *TUBER(B,N)                   &
                  )*DLTDY

       !update leaf biomass (gC/m^2)
       LEAF(B,N) = LEAF(B,N) + DLTDY*DELLF   !leaf mass increase (gC/m^2) in DLTDY (days) time step
                                             !Note LEAF is updated before EP!

       !update epiphyte biomass (gC-Epiphyte)/(gm^2-Leaf)
       EP(B,N) = EP(B,N)*LFOLD*           &
              (1.0                        &
                  + DLTDY                 & 
                   *( PEP(B,N)            & !epiphyte growth rate  (1/day)
                     -BMEP(B,N)           & !epiphyte basal metabolism (1/day)
                     -PREP(B,N)           & !epiphyte loss due to predation (1/day)
                     -SLSH(B)             & !epiphyte loss due to sloughing (1/day)
                    )                     & 
               )  /(LEAF(B,N)+1.0E-20)      !normalize by new LEAF to get EP  
                                            !Note LFOLD is LEAF mass in previous time step
     ENDDO
   ENDDO

! TEMPORARILY SET MINIMUM BIOMASS

   DO I=1,NSAVCELL
     B = SAVCELL(I)
     DO N=1,NSAVSPC(B)  !loop over all species
       ROOT(B,N) = MAX(ROOT(B,N),0.01)
       LEAF(B,N) = MAX(LEAF(B,N),0.01)
       STEM(B,N) = MAX(STEM(B,N),0.01)
       TUBER(B,N) = MAX(TUBER(B,N),0.01)
       EP(B,N) = MAX(EP(B,N),0.01)
     ENDDO
   ENDDO
      
   RETURN
   END SUBROUTINE SAV_COMP

END MODULE MOD_SAV

