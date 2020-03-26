MODULE MOD_WQM   
            
            USE MOD_SIZES, ONLY :       &    !
                NCP,            &  !
                NQFP,           &  !
                !NHQP,           &  !
                NS1P,           &  !
                NS2P,           &  !
                NS3P,           &  !
                NBCP,           &  !
                !NDP,            &  !
                NFLP,           &  !
                !NOIP,           &  !
                NSSFP,             &  !
                MGL,            &  !
                NGL ,            &  !
                !OBCGL,          &  !
                NOBTY  

           USE MOD_PREC, ONLY : SP, CDF_PREC
           
           USE MOD_LIMS, ONLY : MTLOC, NTLOC, KBM1,KB


           
#if defined (MULTIPROCESSOR)
           USE MOD_BUFFERS , ONLY: SNDBUF,RCVBUF
#endif
           IMPLICIT NONE
           SAVE
        !********************************************************************************
        !**                                                                            **
        !**                        COMMON file for CE-QUAL-IC                          **
        !**                                                                            **
        !**                                Version 1.1                                 **
        !**                               July 28, 1992                                **
        !**                                                                            **
        !**                        Water Quality Modeling Group                        **
        !**                        U.S. Army Corps of Engineers                        **
        !**                        Waterways Experiment Station                        **
        !**                        Vicksburg, Mississippi 39180                        **
        !**                                                                            **
        !********************************************************************************

        !********************************************************************************
        !**                                                                            **
        !**                           PARAMETER definitions                            **
        !**                                                                            **
        !**    NBP   - Number of boxes                                                 **
        !**    NQFP  - Number of horizontal and vertical flow faces                    **
        !**    NHQP  - Number of horizontal flow faces                                 **
        !**    NSBP  - Number of boxes in the surface layer                            **
        !**    NLP   - Maximum number of layers                                        **
        !**    NS1P  - Number of source one inputs                                     **
        !**    NS2P  - Number of source two inputs                                     **
        !**    NS3P  - Number of source three inputs                                   **
        !**    NBCP  - Number of boundary concentration inputs                         **
        !**    NMP   - Number of modifications to initial box concentrations           **
        !**    NDP   - Maximum number of days for any output                           **
        !**    NFLP  - Number of files for each type of time-varying data              **
        !**    NCP   - Number of constituent state variables                           ** 
        !**    NSAVM - Maximum Number of SAV Species                                   ** !Wen Long: 
        !**    NSSFP - Number of suspension-feeding species                            **
        !**                                                                            **
        !********************************************************************************


        !********************************************************************************
        !**                        Water Quality Model Setup                           **
        !********************************************************************************


        !
        ! index variables (typically used for do loops)
        !

             INTEGER :: B           !index for overall  boxes (i.e. nodes in fvcom) 
             INTEGER :: BB          !index for surface or bottom boxes (now equivalent to B)
         
             INTEGER :: SAVB        !array storing boxes indices where there is SAV in it in SAV model
                                    !this should be globalized and then broadcast to individual processors
                                    !******!!Never used!!*****

        ! water column constituents 
                                   
                 REAL(SP), ALLOCATABLE, DIMENSION (:,:)  :: &
                                    T,    & !1 water column temperature (degC)
                                    SALT, & !2 salinity in water column (psu)
                                    SSI,  & !3 inorganic solids (g/m^3)
                                    B1,   & !4 algal group 1 biomass concentration (gC/m^3)
                                    B2,   & !5 algal group 2 biomass concentration (gC/m^3)
                                    B3,   & !6 algal group 3 biomass concentration (gC/m^3)
!                                    SZ,   & !7 micro(small) zooplankton concentration (gC/m^3)  !WLong moved to mod_zoop
!                                    LZ,   & !8 meso(large) zooplankton concentration (gC/m^3)     !WLong moved to mod_zoop
                           !         LDOC, & !9 labile dissolved organic carbon concentration (gC/m^3)  ! Moved to mod_wc_doc B Clark Sep 2015
                             !       RDOC, & !10 refractory dissolved organic carbon concentration (gC/m^3)  ! Moved to mod_wc_doc B Clark Sep 2015                           
                                    LPOC, & !11 labile particulate organic carbon concentration (gC/m^3)
                                    RPOC, & !12 refractory particulate organic carbon concentration (gC/m^3)
                                    NH4,  & !13 ammonia concentration in water column (gN/m^3)
                                    NO3,  & !14 nitrate concentration in water column (gN/m^3)
									
									DTNH4_out, DTNO3_out, &!Integrated change in Nh4 and NO3 g N
									
                                    UREA, & !15 urea concentration (gN/m^3)    !Wen Long: this term is never used !!!
                               !     LDON, & !16 labile dissolved organic nitrogen concentration (gN/m^3)  ! ! Moved to mod_wc_doc B Clark Sep 2015
                                !    RDON, & !17 labile dissolved organic nitrogen concentration (gN/m^3)  ! Moved to mod_wc_doc B Clark Sep 2015
                                    LPON, & !18 labile particulate organic nitrogen concentration (gN/m^3)
                                    RPON, & !19 refractory particulate organic nitrogen concentration (gN/m^3)
                                    PO4,  & !20 phosphate concentation in water column (gP/m^3) 
                                 !   LDOP, & !21 labile dissolved organic phosphorus concentration (gP/m^3) ! Moved to mod_wc_doc B Clark Sep 2015
                                 !   RDOP, & !22 refractory dissolved organic phosphorus concentration (gP/m^3)  ! Moved to mod_wc_doc B Clark Sep 2015
                                    LPOP, & !23 labile particulate organic phosphorus concentration (gP/m^3)
                                    RPOP, & !24 refractory particulate organic phosphorus concentration (gP/m^3) 
                                    PIP,  & !25 particulate inorganic phosphorus (adsorbed to SSI) (gP/m^3)
                                    COD,  & !26 chemical oxygen demand in water column (gO2-equiv/m^3)
                                    DOXG, & !27 dissolved oxygen in water column (gO2/m^3)
                                    SIUPB,& !28 particulate biogenic silica (unavaiable ) concentration in water column (gSi/m^3) (why not call it SIUPB ??) !Wen Long changed SU to PSIW to SIUPB
                                    SIAT, & !29 total silica (available) concentration in water column (gSi/m^3) (Wen Long changed SA to DSI and to SIAT)
                                    PIB1, & !Phosphorus in algae 1 (gP/m^3)
                                    PIB2, & !Phosphorus in algae 2 (gP/m^3)
                                    PIB3    !Phosphorus in algae 3 (gP/m^3)
									
									
			REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: DTLPOC_out,DTRPOC_out, &
												   DTLPON_out,DTRPON_out   !  B Clark added for Budget Calculation, Feb 2019

            !Water column settling rates
            !3D Spatially varying (horizontally and vertically) settling velocities (read from STLFN)
            REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: WSS, & !Fixed solids (m/day)
                                               WSL, & !Labile POM (m/day)
                                               WSR, & !Refractory POM (m/day)
                                               WS1, & !Algal group 1 (m/day)
                                               WS2, & !Algal group 2 (m/day)
                                               WS3, & !Algal group 3 (m/day)
                                               WSU, & !Particulate biogenic sillica (m/day)                                                                   
                                               WSSHI  !fixed solids settling rate (m/day) when solids concentration is high
                                                      ! Set to WSS when SSI (fixed solids concentration ) <100 g/m^3
                                                      ! Set to 5 (m/day) when fixed solids concentration > 100 g/m^3                           

             ! KHR1,   &!!Half constant of respiration of algal group 1 to release DOC for oxygen (gO2/m^3) ! B Clark moved to wc_dom Sep 2015
             !        KHR2,   &    !Half constant of respiration of algal group 2 to release DOC for oxygen (gO2/m^3)
             !        KHR3,   &    !Half constant of respiration of algal group 3 to release DOC for oxygen (gO2/m^3)
         REAL(SP) :: KHN1,   &    !Algal group 1 half saturation concentration for nitrogen uptake (gN/m^3)
                     KHN2,   &    !Algal group 2 half saturation concentration for nitrogen uptake (gN/m^3)
                     KHN3,   &    !Algal group 3 half saturation concentration for nitrogen uptake (gN/m^3)
                     KHNAVG, &    !average of KHN1, KHN2 and KHN3                         
                     KHP1,   &    !Algal group 1 half saturation concentration for phosphorus uptake (gP/m^3)
                     KHP2,   &    !Algal group 2 half saturation concentration for phosphorus uptake (gP/m^3)
                     KHP3,   &    !Algal group 3 half saturation concentration for phosphorus uptake (gP/m^3)
                     KHPAVG, &    !Average of KHP1, KHP2 and KHP3
                     KHONT,  &    !Half saturation constant of dissolved oxygen required for nitrification (gO2/m^3)
                     KHOCOD, &    !Half saturation constant of dissolved oxygen required for excertion of chemical oxygen demand (gO2m^3)
                     KHODOC, &    !Half saturation constant of dissolved oxygen required for oxic respiration (gO2/m^3)
                     KHNNT,  &    !Half saturation constant of NH4 required for nitrification (gN/m^3)
                     KHNDN,  &    !Half saturation concentration of nitrate required for denitrification (0.1gN/m^3)
                     KHNFIX, &    !Half constant of DIN (NH4+NO3) at which nitrogen limitation = 2KHNFIX/(2KHNFIX+KHN1) (Not in the origianl ICM model!!!)
                     FNFIX,  &    !Fraction of N fixers (not in the original ICM model!)
                     KHS1,   &    !Half saturation concentration for silica uptake (gSi/m^3) by algal group 1
                     KHS2,   &    !Half saturation concentration for silica uptake (gSi/m^3) by algal group 2
                     KHS3,   &    !Half saturation concentration for silica uptake (gSi/m^3) by algal group 3
                     KHST1,  &    !Half saturation constant of salinity at which mortality is half the maximum value for algal group 1 (ppt)
                     KHST2,  &    !Half saturation constant of salinity at which mortality is half the maximum value for algal group 2 (ppt)
                     KHST3        !Half saturation constant of salinity at which mortality is half the maximum value for algal group 3 (ppt)

        !
        !parameters for specifying benthic flux for the water column model
        !original water column model was not coupled with sediment flux model
        !and benthic fluxes are specified using these parameters
        !                 
             REAL(SP) :: KSO,        &  !effect of temperature on sediment oxygen consumption (SOD) (1/degC)
                     KSNH4,      &  !effect of temperature on ammonium flux (1/degC)
                     KSNO3,      &  !effect of temperature on nitrate flux (1/degC)
                     KSPO4,      &  !effect of temperature on phosphate flux (1/degC)
                     KSSA,       &  !effect of temperature on dissolved silica flux (1/degC)
                     KSDOC,      &  !effect of temperature on dissolved organic carbon flux (1/degC)
                     TRSO,       &  !reference temperature for specification of SOD (degC)                          
                     TRSNH4,     &  !reference tempreature for ammonium flux (degC)
                     TRSNO3,     &  !reference temperature for nitrate flux (degC)
                     TRSPO4,     &  !reference temperature for phosphate flux (degC)
                     TRSSA,      &  !reference temperature for dissolved silica flux (degC)
                     TRSDOC,     &  !reference temperature for dissolved organic carbon flux          
                     MTCNO3,     &  !Nitrate mass transfer coefficient (m/day)
                     SEDNO3,     &  !Sediment interstitial nitrate concentration gN/m^3
                     KHSO           !half saturation constant of SOD for oxygen (gO2/m^3)

             REAL(SP),ALLOCATABLE :: BENDOCB(:),  &  !benthic DOC flux (gC/m^2/day) at reference temperature
                                 BENNH4B(:),  &  !benthic NH4 flux (gN/m^2/day) at reference temperature
                                 BENNO3B(:),  &  !benthic NO3 flux (gN/m^2/day) at reference temperature
                                 BENPO4B(:),  &  !benthic PO4 flux (gN/m^2/day) at reference temperature
                                 BENCODB(:),  &  !benthic COD flux (gN/m^2/day) at reference temperature
                                 BENDOB(:),   &  !benthic SOD flux (gO2/m^2/day) at reference temperature
                                 BENSAB(:)       !benthic dissolved silica flux (gSi/m^2/day) at reference temperature

             !Temperature control on nitriciation
             REAL(SP) :: TMNT,  &  !optimal temperature for nitriciation
                     KTNT1, &  !effect of temperature below TMNT on nitrificatin (1/degC^2)
                     KTNT2     !effect of temperature above TMNT on nitrifcation (1/degC^2)
                
             REAL(SP) :: KTCOD, &  !effect of temperature on excertion of COD (1/day)
                         KTMNL, &  !effect of temperture on minerailation rate (1/degC)
                         KTHDR, &  !effect of temperature on hydrolysis rates (1/degC)
                         KTSUA     !effect of temperature on particulate biogenic silica dissolution rates (1/degC)
             
             REAL(SP):: KHCOAG !salinity at which coagulation of RDOC and RDON is half the maxium value (ppt)                                                                                             
             
             INTEGER :: HYDPTR, METPTR, CBCPTR, S1PTR,  S2PTR,  S3PTR,  BFIPTR,     &
                        ATMPTR, SAVPTR, &
                        BAIPTR, &!WLong: pointer to benthic algae input
                        WCLPTR, &!WLong: pointer to sediment flux model water column input
                        NTVWCLF   !number of sediment flux model water column input files

             REAL(SP) :: KHTIS ! Never used

             REAL(SP) :: NXCBC,  JDAY,   JDAYMBL

             REAL(SP) :: KT  !Coef of heat exchange (watts/m^2/C)

             REAL(SP) :: NXMET,  NXS1,   NXS2,   NXS3,   NXBFI,  NXATM,  NXSAV

             REAL(SP) :: I0,   & !irradiance (E/m^2/sec)
                     IT,   & !irradiance (E/m^2/sec)
                     ITNX, & !irradiance (E/m^2/sec)
                     ISMIN   !NOT used in the code

             REAL(SP) :: KADPO4,  & !Partition coefficient for soprtion of phosphate on inorganic solids  (m^3/gSolids)
                     KADSA      !Partition coefficient for sorption of silica on inorganic solids (m^3/gSolids)
                                !SA-- available silica
                                !SUA -- unavailable particulate biogenic silica
                     
                     
        !Wen Long moved all these to the mod_zoop.F
        !!Zooplankton parameters
        !                         
        !     REAL(SP) :: CTSZ,    &  !carbon threshold for grazing (gC/m^3) for microzooplankton
        !             CTLZ,    &  !carbon threshold for grazing (gC/m^3) for mesozooplankton
        !             KHCSZ,   &  !prey density at which grazing is halved (gC/m^3) for microzooplankton
        !             KHCLZ,   &  !prey density at which grazing is halved (gC/m^3) for mesozooplankton
        !             MZEROSZ, &  !mortality at zero dissolved oxygen (1/day) for microzooplankton
        !             MZEROLZ     !mortality at zero dissolved oxygen (1/day) for mesozooplankton

        !     REAL(SP) :: UB1SZ, &   !Utilization of algal group 1 by microzooplankton, range [0,1]
        !             UB2SZ, &   !Utilization of algal group 2 by microzooplankton, range [0,1]
        !             UB3SZ, &   !Utilization of algal group 3 by microzooplankton, range [0,1]
        !             UB1LZ, &   !Utilization of algal group 1 by mesozooplankton, range [0,1]
        !             UB2LZ, &   !Utilization of algal group 2 by mesozooplankton, range [0,1]
        !             UB3LZ, &   !Utilization of algal group 3 by mesozooplankton, range [0,1]                         
        !             UDSZ,  &   !Utilization of dissolved organic carbon by microzooplankton, range [0,1]
        !             ULDSZ, &   !Utilization of labile dissolved organic carbon by microzooplankton, range [0,1]
        !             URDSZ, &   !Utilization of refractory dissolved organic carbon by microzooplankton, range [0,1]
        !             ULPSZ, &   !microzooplankton utilization of labile particulate organic carbon, range [0,1]
        !             URPSZ, &   !microzooplankton utilization of refractory particulate organic carbon, range [0,1]
        !             ULLZ,  &   !mesozooplankton utilization of labile particulate organic carbon, range [0,1]
        !             URLZ,  &   !mesozooplankton utilization of refractory particulate organic carbon, range [0,1]
        !             USZLZ, &   !mesozooplankton utilization of microzooplankton, range [0,1], (~1.0)
        !             TRSZ,  &   !reference temperature for microzooplankton metabolism (degC)
        !             TRLZ,  &   !reference temperature for mesozooplankton metabolism (degC)
        !             DOCRITSZ, &!dissolved oxygen below which mortality occurs for microzooplankton (gO2/m^3)
        !             DOCRITLZ, &!dissolved oxygen below which mortality occurs for mesozooplankton (gO2/m^3)
        !             ANCSZ,  &  !microzooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
        !             ANCLZ,  &  !mesozooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
        !             APCSZ,  &  !microzooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
        !             APCLZ,  &  !mesozooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
        !             AOCRSZ, &  !ratio of oxygen consumed to microzooplankton carbon metabolized (gO2/gC) (~2.67)
        !             AOCRLZ, &  !ratio of oxygen consumed to mesozooplankton carbon metabolized (gO2/gC) (~2.67)
        !             FRSASZ, &  !fraction of microzooplankton silica recycled to dissolved silica pool, range [0,1]
        !             FRSALZ     !fraction of mesozooplankton silica recycled to dissolved silica pool, range [0,1]
        !                         
        !     REAL(SP) :: FLDOCSZ,&    !fraction of microzooplankton carbon released to LDOC, range [0,1]
        !             FRDOCSZ,&    !fraction of microzooplankton carbon released to RDOC, range [0,1]
        !             FLPOCSZ,&    !fraction of microzooplankton carbon released to LPOC, range [0,1]
        !             FRPOCSZ,&    !fraction of microzooplankton carbon released to RPOC, range [0,1]
        !             FLDONSZ,&    !fraction of microzooplankton nitrogen released to LDON, range [0,1]
        !             FRDONSZ,&    !fraction of microzooplankton nitrogen released to RDON, range [0,1]
        !             FLPONSZ,&    !fraction of microzooplankton nitrogen released to LPON, range [0,1]
        !             FRPONSZ,&    !fraction of microzooplankton nitrogen released to RPON, range [0,1]
        !             FLDOPSZ,&    !fraction of microzooplankton phosphorus released to LDOP range [0,1]    
        !             FRDOPSZ,&    !fraction of microzooplankton phosphorus released to RDOP range [0,1]    
        !             FLPOPSZ,&    !fraction of microzooplankton phosphorus released to LPOP range [0,1]    
        !             FRPOPSZ,&    !fraction of microzooplankton phosphorus released to RPOP range [0,1]    
        !             FNH4SZ, &    !fraction of microzooplankton nitrogen recycled to DIN as NH4, range [0,1]
        !             FPO4SZ, &    !fraction of microzooplankton phosphorus recycled to DIP as PO4, range [0,1]
        !             FUREASZ,&    !Never used !!!                         
        !             FLDOCLZ,&    !fraction of mesozooplankton carbon released to LDOC, range [0,1]
        !             FRDOCLZ,&    !fraction of mesozooplankton carbon released to RDOC, range [0,1]
        !             FLPOCLZ,&    !fraction of mesozooplankton carbon released to LPOC, range [0,1]
        !             FRPOCLZ,&    !fraction of mesozooplankton carbon released to RPOC, range [0,1]
        !             FLDONLZ,&    !fraction of mesozooplankton nitrogen released to LDON, range [0,1]
        !             FRDONLZ,&    !fraction of mesozooplankton nitrogen released to RDON, range [0,1]
        !             FLPONLZ,&    !fraction of mesozooplankton nitrogen released to LPON, range [0,1]
        !             FRPONLZ,&    !fraction of mesozooplankton nitrogen released to RPON, range [0,1]
        !             FLDOPLZ,&    !fraction of mesozooplankton phosphorus released to LDOP range [0,1]    
        !             FRDOPLZ,&    !fraction of mesozooplankton phosphorus released to RDOP range [0,1]    
        !             FLPOPLZ,&    !fraction of mesozooplankton phosphorus released to LPOP range [0,1]    
        !             FRPOPLZ,&    !fraction of mesozooplankton phosphorus released to RPOP range [0,1]    
        !             FNH4LZ, &    !fraction of mesozooplankton nitrogen recycled to DIN as NH4, range [0,1]
        !             FPO4LZ, &    !fraction of mesozooplankton phosphorus recycled to DIP as PO4, range [0,1]
        !             FUREALZ      !Never used!!!

        !     REAL(SP),ALLOCATABLE::  B1ASZ(:,:),   & 
        !                         B2ASZ(:,:),   &
        !                         B3ASZ(:,:),   &
        !                         LPOCASZ(:,:), &
        !                         RPOCASZ(:,:), &
        !                         PRASZ(:,:),   &
        !                         B1ALZ(:,:),   &
        !                         B2ALZ(:,:),   &
        !                         B3ALZ(:,:),   &
        !                         SZALZ(:,:),   &
        !                         LPOCALZ(:,:), &
        !                         RPOCALZ(:,:), &
        !                         PRALZ(:,:),   &
        !                         CLSZ(:,:),    & 
        !                         CLLZ(:,:),    &
        !                         RSZ(:,:),     &
        !                         RLZ(:,:),     &
        !                         RMAXSZ(:,:),  &
        !                         RMAXLZ(:,:),  &
        !                         BMSZ(:,:),    &
        !                         BMLZ(:,:),    &
        !                         BMRSZ(:,:),   &
        !                         BMRLZ(:,:),   &
        !                         MSZ(:,:),     &
        !                         MLZ(:,:),     &
        !                         PRSZLZ(:,:),  &
        !                         GSZ(:,:),     &
        !                         GLZ(:,:),     &
        !                         ESZ(:,:),     &
        !                         ELZ(:,:),     &
        !                         RFSZ(:,:),    &
        !                         RFLZ(:,:),    &
        !                         PRSZ(:,:),    &
        !                         PRLZ(:,:),    &
        !                         LDOCASZ(:,:), &
        !                         BPRSZ(:,:),   &
        !                         BPRLZ(:,:),   &
        !                         RDOCASZ(:,:), &
        !                         LDOCSZ(:,:),  &
        !                         LPOCSZ(:,:),  &
        !                         RPOCSZ(:,:),  &
        !                         LDOCLZ(:,:),  &
        !                         LPOCLZ(:,:),  &
        !                         RPOCLZ(:,:),  &
        !                         NH4SZ(:,:),   &
        !                         LDONSZ(:,:),  &
        !                         LPONSZ(:,:),  &
        !                         RPONSZ(:,:),  &
        !                         NH4LZ(:,:),   &
        !                         LDONLZ(:,:),  &
        !                         LPONLZ(:,:),  &
        !                         RPONLZ(:,:),  &
        !                         PO4SZ(:,:),   &
        !                         LDOPSZ(:,:),  &
        !                         LPOPSZ(:,:),  &
        !                         RPOPSZ(:,:),  &
        !                         PO4LZ(:,:),   &
        !                         LDOPLZ(:,:),  &
        !                         LPOPLZ(:,:),  &
        !                         RPOPLZ(:,:),  &
        !                         RDOCSZ(:,:),  &
        !                         RDONSZ(:,:),  &
        !                         RDOPSZ(:,:),  &
        !                         RDOCLZ(:,:),  &
        !                         RDONLZ(:,:),  &
        !                         RDOPLZ(:,:),  &
        !                         PIB1SZ(:,:),  &
        !                         PIB2SZ(:,:),  &
        !                         PIB3SZ(:,:),  &
        !                         PIB1LZ(:,:),  &
        !                         PIB2LZ(:,:),  &
        !                         PIB3LZ(:,:),  &
        !                         B1SZ(:,:),    &
        !                         B2SZ(:,:),    &
        !                         B3SZ(:,:),    &
        !                         B1LZ(:,:),    &
        !                         B2LZ(:,:),    &
        !                         B3LZ(:,:),    &
        !                         DOSZ(:,:),    &
        !                         DOLZ(:,:),    &
        !                         SASZ(:,:),    &
        !                         SUSZ(:,:),    &
        !                         SALZ(:,:),    &
        !                         SULZ(:,:)  
        !     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: ACLSZ,   ACLLZ,    ARSZ,    ARLZ,     &
        !               ABMSZ,   ABMLZ,    AMSZ,    AMLZ,     &
        !               APRSZLZ, AGSZ,     AGLZ,    ADOCSZ,   &
        !               APOCSZ,  ADOCLZ,   APOCLZ,                 &
        !               ANH4SZ,  ADONSZ,   APONSZ,                 &
        !               ANH4LZ,  ADONLZ,   APONLZ,                 &
        !               APO4SZ,  ADOPSZ,   APOPSZ,                 &
        !               APO4LZ,  ADOPLZ,   APOPLZ,                 &
        !               APRSZ,   APRLZ,    APISZ,  APILZ

        !     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: AB1SZ,   AB2SZ,    AB3SZ,  AB1LZ,     &
        !               AB2LZ,   AB3LZ,    ADOSZ,  ADOLZ,     &
        !               ASASZ,   ASUSZ,    ASALZ,  ASULZ

        !     REAL(SP),DIMENSION(-50:400) :: FTLZ,  FTSZ,  FTBMSZ,          &
        !               FTBMLZ,FTPRSZ,FTPRLZ,          

        
        !!Tempreature control on algae, should move to mod_algal.F   !!LB: commented because we now use get_ft1,get_ftbm1,etc
        !       REAL(SP),DIMENSION(-50:400) ::  FT1,            &
        !                                   FT2,            &
        !                                   FT3,         &
        !                                FTBM1,            &
        !                                FTBM2,            &
        !                                FTBM3,          &
        !                                FTPR
                     
                     
        !  variables needed by autostep and main   VJP 10/11/04
             INTEGER,SAVE :: NXHYD = 0
             REAL(SP),SAVE :: COURMX, COURQS, COURVS
             REAL(SP),SAVE :: DIFFMX, DIFFAS, DIFFDS, FNDLT

             REAL(SP),ALLOCATABLE :: V1SINGLE(:,:)

             REAL(SP),SAVE :: DLT8, MAXDLTDP, MXDLTDP, NXHYDDP, ELTMS1

             REAL(SP) :: ELTMS 

             LOGICAL :: BOUNDARY_CONC,             &    !
                        SOURCE_ONE,               &    !
                        SOURCE_TWO,             &    !
                        SOURCE_THR,             &    !
                        SOLIDS_CALC,              &    !
                        PIP_CALC,               &    !
                        BENTHIC_FLUXES,         &    !
                        SEDIMENT_CALC,            &    !
                        LIGHT_EXTINCTION,       &    !
                        ATMOS_LOADS,            &    !
                        TRANSPORT_FLUXES,          &    !
                        AVERAGE_PLOTS,          &    !
                        QUALITY_DIAG,             &    !
                        KINETIC_FLUXES,         &    !
                        SEDIMENT_DIAG,          &    !
                        DIAGNOSTICS,              &    !
                        CONSERVE_MASS,          &    !
                        SETTLING,               &    !
                        STEP_BOUNDARY,            &    !
                        MASS_BALANCE      
                        
             LOGICAL :: SAV_LOADS,              &
                        SAV_CALC,                 &
                        SAV_PLOTS                
                        
             LOGICAL :: FLOW,                     &
                        XY_DIFFUSION,             &
                        Z_DIFFUSION,             &
                        BINARY_HYDRO,             &
                        ASCII_HYDRO,              &    
                        DEPTH_AVG_HYDRO,         &
                        SIGMA_HYDRO             
             CHARACTER(LEN=8) ::  SPVARKE, PRINTKE
             CHARACTER(LEN=72) :: KEIFN,         &
                                  SVOFN,         &
                                  BFOFN,         &
                                  ZOOFN,         &
                                  ZFOFN,         &
                                  ALOFN,         &
                                  BAOFN,        &!WLong Benthic Algae Output
                                  DFOFN             !WLong, Deposition feeder output
                                  
             CHARACTER(LEN=8) :: MINPROD

        !***** Dimension declarations

             REAL(SP),ALLOCATABLE :: DTM(:,:,:)

        !RGL
             REAL(SP), ALLOCATABLE :: DEPTHR(:,:)    !Never used
             REAL(SP), ALLOCATABLE,DIMENSION(:,:,:) ::     C1,     &
                                                        C2,     &
                                                        C2F,     &
														C3,    & ! Added by B Clark for advection model 
                                                        CSTAR,     &
                                                        DTC,     &!g/m^3/sec 
                                                        AC1
														
		    REAL(SP), ALLOCATABLE, DIMENSION(:,:,:)  :: XFLUX,XFLUX_ADV
			
        ! KURT GLAESEMANN add C2_GL
             REAL(SP), ALLOCATABLE,DIMENSION(:,:,:) :: C2_GL
             REAL(SP), ALLOCATABLE,DIMENSION(:,:)   :: CCHL1_GL,T_GL,S_GL ! RGL added 2nd
             REAL(SP), ALLOCATABLE,DIMENSION(:)     :: total_netPP_GL
             REAL(SP), ALLOCATABLE,DIMENSION(:)     :: D_GL,H_GL,EL_GL    !Wen Long, D_GL : total depth (global)
                                                                     !          H_GL : still depth (global)
                                                                     !         EL_GL : surface elevation (global)
        !Wen Long moved following BFLUX terms from tvds.F to here
             REAL(SP), ALLOCATABLE,DIMENSION(:,:)     :: BFLUX_GL
             REAL(SP), ALLOCATABLE,DIMENSION(:,:)     :: BFLUXB_GL
             REAL(SP), ALLOCATABLE,DIMENSION(:,:)     :: BFLUXNX, BFLUXNX_GL, BFLUX


        !Wen Long moved the following variables from ncdio.F to here

             REAL(CDF_PREC), ALLOCATABLE, DIMENSION(:,:)  :: UL_GL,VL_GL                 !u and v
			 REAL(CDF_PREC), ALLOCATABLE, DIMENSION(:,:)  :: WW_GL                      !ww
             REAL(CDF_PREC), ALLOCATABLE, DIMENSION(:,:)  :: WTSL_GL,KHL_GL              !wts, kh
             REAL(CDF_PREC), ALLOCATABLE, DIMENSION(:,:)  :: SL_GL,TL_GL                 !salinity,temp
             REAL(CDF_PREC), ALLOCATABLE, DIMENSION(:)    :: ELL_GL,DTFAL_GL             !zeta,dtfa
             REAL(CDF_PREC), ALLOCATABLE, DIMENSION(:)    :: UARD_OBCN_GL                !uard_obcn
             REAL(CDF_PREC), ALLOCATABLE, DIMENSION(:,:)  :: XFLUX_OBC_GL                !xflux_obc
			 
			 REAL(SP),ALLOCATABLE,DIMENSION(:,:)::WW_OUT ! B clark Dec 2018 vertical velocity output for netcdf


             REAL(SP),DIMENSION(NCP) :: CMASS    !Total mass of all constiutents (kg)

             REAL(SP) :: S1L(NS1P,NCP),    S2L(NS2P,NCP),    S3L(NS3P,NCP)
             INTEGER :: S1LB(NS1P,NCP),   S2LB(NS2P,NCP),   S3LB(NS3P,NCP)

             REAL(SP) :: CB(NBCP,NCP),CBNX(NBCP,NCP), CBOLD(NBCP,NCP)
             
             REAL(SP),ALLOCATABLE,DIMENSION(:) ::    CTEMP      !sediment temperature (degC)
           
            !!The following needs to be moved to SED_INIT 
            ! REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: CPOP,   &  !Sediment POP (mgP/gSediment)
            !                                         CPON,   &  !Sediment PON (mgN/gSediment)
            !                                         CPOC       !Sediment POC (mgN/gSediment)

            !                                    
           
             
            !REAL(SP),ALLOCATABLE,DIMENSION(:) ::  CPOS,    &  !Sediment particulate organic
            !                                              !silica (mgSi/m^3)
            !                                  CPO4,    &  !Sediment inorganic phosphorus (mgP/m^3))
            !                                  CNO3,    &  !sediment NO3 concentration
            !                                              !(2nd layer) (mgC/m^3)
            
            !                                  !CDTEMP, &  !never used 
            !                                  !CNH4,    &  !sediment NH4 concentration    !Noved to mod_sed.F
            !                                  CCH4,    &  !sediment CH4 concentration
            !                                              !mgO2/^3
            !                                  CSO4,    &  !sediment SO4 concentration 
            !                                              !mgO2/m^3
            !                                  CHS,     &  !Sediment HS concentration 
            !                                              !(mgO2/m^3)
            !                                  CSI         !Silicate concentration
            !                                              !(mgSi/m^3)

             
             !REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: JPOP, & !POP flux (mgP/m^2/day)
             !                                             !positive into sediments
             !                                   JPON, & !PON flux (mgN/m^2/day)
             !                                             !positive into sediments
             !                                   JPOC    !POC flux (mgC/m^2/day)
             !                                             !positive into sediments

             !REAL(SP),ALLOCATABLE,DIMENSION(:) :: JPOS         !Particulate organic silica
             !                                                !flux, (mgSi/m^2/day) positive
             !                                                !into sediments 

             !REAL(SP),ALLOCATABLE,DIMENSION(:) :: HSED         !sediment layer thickness (H2 (m))
             REAL(SP),ALLOCATABLE,DIMENSION(:) :: BSVOL          !WLong: never used

            !water column nitrogen balance check  (positive into water column)
             REAL(SP),ALLOCATABLE,DIMENSION(:) :: ATMFLXNB,    & !atmospheric flux of N
                                                             !(kgN/day), positive into
                                                             !water
                                              ATMFLXPB,    & !atmospheric flux of P
                                                             !(kgP/day), positive into
                                                             !water
                                              ATMFLXCB       !atmospheric flux of carbon
                                                             !(kgC/day), positive into
                                                             !water
            REAL(SP),ALLOCATABLE,DIMENSION(:)::  BENFLXPNB,   & !benthic PON flux (kgN/day)
                                                             !positive into water
                                              BENFLXDNB,   & !benthic  dissolved
                                                             !inorganic nitrogen  flux
                                                             !(NH4 and NO3) (kgN/day)
                                              BENFLXPPB,   & !benthic POP flux (kgP/day)
                                              BENFLXDPB,   & !benthic dissolved PO4 flux
                                                             !(kgP/day)
                                              BENFLXPCB      !benthic POC flux (kgC/day)

            !sediment nitrogen check (positive into sediments)
            REAL(SP),ALLOCATABLE,DIMENSION(:) ::    DLSEDKNB,    &     !sediment denitrification
                                                                    !kinetic flux (kgN/day)
                                                    DLSEDKCB,    &     !sediment carbib diagenesis
                                                                    !kinetic flux (kgC/day)
                                                    BURIALFLXNB, &     !sediment nitrogen flux due
                                                                    !to burial (kgN/day)
                                                    BURIALFLXPB, &     !sediment P flux due to
                                                                    !burial (kgP/day) positive
                                                                    !into sediments
                                                    BURIALFLXCB        !sediment C flux due to
                                                                    !buria (kgC/day), positive

             REAL(SP) :: FLUXT(0:NQFP,NCP), &      !fluxes through all faces in the model (this is old ICM) 
                        AFLUX(NQFP,13)            !fluxes through all faces of 13 chosen species

             !might be used for calculating flxues through all TCE faces in the water
             !column, WLong: need to separate horizontal face fluxes and
             !surface/bottom fluxes
             !NQFP should be replaced by 2D arrays  instead
             
            REAL(SP),DIMENSION(NQFP) ::  FLXTTEM,        & !1
                                        FLXTSAL,        & !2
                                        FLXTSSI,        & !3
                                        FLXT1,          & !4
                                        FLXT2,          & !5
                                        FLXT3,          & !6
                                        FLXTLDOC,       & !7
                                        FLXTLPOC,       & !8
                                        FLXTRPOC,       & !9
                                        FLXTNH4,        & !10
                                        FLXTNO3,        & !11
                                        FLXTLDON,       & !12
                                        FLXTLPON,       & !13
                                        FLXTRPON,       & !14
                                        FLXTPO4,        & !15
                                        FLXTLDOP,       & !16
                                        FLXTLPOP,       & !17
                                        FLXTRPOP,       & !18
                                        FLXTCOD,        & !19
                                        FLXTDO,         & !20
                                        FLXTSU,         & !21
                                        FLXTSA,         & !22
                                        FLXTSZ,         & !23
                                        FLXTLZ,         & !24
                                        FLXTRDOC,       & !25
                                        FLXTRDON,       & !26
                                        FLXTRDOP,       & !27
                                        FLXTUREA,       & !28
                                        FLXTPIP,        & !29
                                        FLXTPIB1,       & !30
                                        FLXTPIB2,       & !31
                                        FLXTPIB3          !32

             !settling flux of all particulate constituents in water column through
             !bottom face of each layer for each TCE 
             REAL(SP),ALLOCATABLE :: FLUXS(:,:,:)  !size MTLOC x KBM1 x 13
                                               !13 settling fluxes above :
                                               !1  SSI    
                                               !2  B1
                                               !3  B2
                                               !4  B3
                                               !5  POC 
                                               !6  PON calculated by NITROG() in wqm_main.F
                                               !       which is dummy FLXSPON in
                                               !       NITROG() method with size
                                               !       MTLOC x KBM1, need to be careful
                                               !       here passing parts of FLUXS to
                                               !       FLXSPON(:,:)
                                               !7  PO4
                                               !8  POP
                                               !9  PIP
                                               !10 Si
                                               !11 PIB1
                                               !12 PIB2
                                               !13 PIB3
             !never used
             REAL(SP),ALLOCATABLE,DIMENSION(:,:) ::        FLXSPIB1,    &  !11 of FLUXS above
                                                        FLXSPIB2,    &  !12 of FLUXS above
                                                        FLXSPIB3        !13 of FLUXS above

             REAL(SP),ALLOCATABLE,DIMENSION(:) ::     ABENDOC,        &
                                                    ABENNH4,        &
                                                    ABENNO3,        &
                                                    ABENPO4,        &
                                                    ABENCOD,        &
                                                    ABENDO,         &
                                                    ABENSA,         &
                                                    ABENCH4G,        &
                                                    ABENCH4A,        &
                                                    ASSFWS,            &
                                                    APCFWS,         &
                                                    APNFWS,            &
                                                    APPFWS,            &
                                                    APSFWS,         &
                                                    ACPIP,            &
                                                    ACPOS,          &
                                                    ADFEED

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: ACPOC,            &
                                                    ACPON,            &
                                                    ACPOP,            &
                                                    ASFEED            
             
             REAL(SP),ALLOCATABLE,DIMENSION(:) ::    AJNSF,             &
                                                    AJPSF,            &
                                                    ASODSF,         &
                                                    ASASF,            &
                                                    ASUSF,             &
                                                    ASFGCIN,        &
                                                    ASFCFEC,        &
                                                    ASFCPSF,        &
                                                    AFLXCSF,        &
                                                    AFLXNSF,        &
                                                    AFLXPSF,        &
                                                    ARPOCSF,        &
                                                    ARPONSF,        &
                                                    ARPOPSF,        &
                                                    ASSISF,         &
                                                    ASSISASF,        &
                                                    ASSISUSF,        &
                                                    ASSIPSF         

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: ACCHL1,          &
                                                    ACCHL2,            &
                                                    ACCHL3

             INTEGER :: S1LN(NCP),        &
                        S2LN(NCP),        &
                        S3LN(NCP),        &
                        AC(NCP),        &
                        NCB(NCP)        !J    ,IBC(NCP,NBCP)
             
             REAL(SP),DIMENSION(0:NQFP) :: Q, A, DIFF

!WLong removed the following variables
!             REAL(SP),ALLOCATABLE ::                     BL(:,:,:)
!             REAL(SP),ALLOCATABLE,DIMENSION(:,:) ::     V1S,    &
!                                                        HMV,    &
!                                                        ZD            
!             REAL(SP),ALLOCATABLE,DIMENSION(:) :: HMBV    

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) ::     V1,    &
                                                        V2

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: PN1,     &
                                                    PN2,    &
                                                    PN3,    &
                                                    BM1,      &
                                                    BM2,    &
                                                    BM3,      &
                                                    FTCOD,  &
                                                    NPP,     &
                                                    P1,        &
                                                    P2,        &
                                                    P3,     &
                                                    NT,     &
													NT_OUT, &!
													NT_GL,   & ! Global loss of Nh4 due to nitrification, B Clark July 2018
                                                    PR1,     &
                                                    PR2,    &
                                                    PR3,    &
                                                    RATOX,  &
                                                    GPP,    &
                                                    P1NNF
                                                
             REAL(SP),ALLOCATABLE,DIMENSION(:) ::     ASRAT,            &
                                                total_netPP

             REAL(SP),ALLOCATABLE,DIMENSION(:) ::     CFIX,            &
                                                    SNFIX
             
             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: AFI1,        &
                                                    ANL1,        &
                                                    APL1,        &
                                                    AFI2,          &
                                                    ANL2,        &
                                                    APL2,         &
                                                    ASL2,         &
                                                    AFI3,        &
                                                    ANL3,        &
                                                    APL3,        &
                                                    ANPP,        &
                                                    ARESP,         &
                                                    AKE,        &
                                                    ASL1,        &
                                                    ASL3,        &
                                                    AGPP        

             REAL(SP),ALLOCATABLE,DIMENSION(:) ::     AASRAT,         &
                                                    ACFIX

             REAL(SP),ALLOCATABLE,DIMENSION(:) ::     AFIB,   &
                                                    ANLB,   &
                                                    APLB,   &
                                                    ANPPB,   &
                                                    ABBM,   &
                                                    ABLITE,   &
                                                    ABMB,   &
                                                    APB,   &
                                                    APRB,   &
                                                    ABADOC,   &
                                                    ABAPOC,   &
                                                    ABANH4,   &
                                                    ABANO3,   &
                                                    ABAPON,   &
                                                    ABAPO4,   &
                                                    ABAPOP,   &
                                                    ABADO
        
         
             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: PM1,   & !Maximum photosynthetic rate for algal group 1 (gC/gCHLA/day)
                                                    PM2,   & !Maximum photosynthetic rate for algal group 2 (gC/gCHLA/day)
                                                    PM3,   & !Maximum photosynthetic rate for algal group 3 (gC/gCHLA/day)
                                                    BMR1,  & !Algal group 1 basal metabolism rate at reference temperature (1/day)
                                                    BMR2,  & !Algal group 2 basal metabolism rate at reference temperature (1/day)
                                                    BMR3,  & !Algal group 3 basal metabolism rate at reference temperature (1/day)
                                                    BPR1,  & !First order predation rate on algal group 1 (1/day)
                                                    BPR2,  & !First order predation rate on algal group 2 (1/day)
                                                    BPR3,  & !First order predation rate on algal group 3 (1/day)
                                                    CCHL1, & !Carbon to chlorophyll ratio (gC/gCHLA) for algal group 1
                                                    CCHL2, & !Carbon to chlorophyll ratio (gC/gCHLA) for algal group 2
                                                    CCHL3, & !Carbon to chlorophyll ratio (gC/gCHLA) for algal group 3  
                                                    Q1,    & !Phosphorus to Carbon ratio (gP/gC) for algal group 1 
                                                    Q2,    & !Phosphorus to Carbon ratio (gP/gC) for algal group 2 
                                                    Q3       !Phosphorus to Carbon ratio (gP/gC) for algal group 3 


        !     REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: DTT, DTSSI, DTB1, DTB2, DTB3,      &
        !      REAL(SP),ALLOCATABLE,DIMENSION(:,:) ::   &
        !               DTLDOC, DTLPOC, DTRPOC, DTNH4, DTNO3, DTLDON, DTLPON,       &
        !               DTRPON, DTPO4, DTLDOP, DTLPOP, DTRPOP, DTCOD, DTDO, DTSIUPB,   &
        !               DTSIAT, DTSZ, DTLZ, DTRDOC, DTRDON, DTRDOP, DTUREA, DTPIP,    &
        !                DTUREA
        !               DTPIB1, DTPIB2, DTPIB3

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: FI1,         &
                                                    FI2,         &
                                                    FI3,         &
                                                    NL1,         &
                                                    NL2,         &
                                                    NL3,          &
                                                    PL1,         &
                                                    PL2,         &
                                                    PL3,         &
                                                    SL2,         &
                                                    RESP,         &
    !                                                KESS,         &! WLong mooved this to mod_owq.F
                                                    SL1,         &
                                                    SL3,         &
                                                    !IAVG,         &!    WLong moved this to mod_owq.F
                                                    !IATBOT,    &!     Wlong moved this to mod_owq.F
                                                    KE,         &!    Not used, should remove
                                                    KEISS,         &!    Not used, should remove
                                                    KEVSS,         &!  Not used, should remove
                                                    KEDOC         !  Not used, should remove

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: FTMNL,       &
                                                    FTHDR

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: RESPC,       &
                                                    DLALGC

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: IWCMNB,      &
                                                    IWCMPB,      &
                                                    IWCMCB,      &
                                                    IWCMSB,        &
                                                    WCMNB,       &
                                                    WCMPB,       &
                                                    WCMCB,       &
                                                    WCMSB,      &
                                                    DLWCMNB,     &
                                                    DLWCMPB,     &
                                                    DLWCKMNB,      &
                                                    DLWCKMCB,    &
                                                    S1FLXNB,     &
                                                    S1FLXPB,     &
                                                    S1FLXCB,    &
                                                    S2FLXNB,     &
                                                    S2FLXPB,     &
                                                    S2FLXCB,    &
                                                    S3FLXNB,     &
                                                    S3FLXPB,     &
                                                    S3FLXCB

             REAL(SP),DIMENSION(NCP) :: DLWCMCB

             REAL(SP),ALLOCATABLE,DIMENSION(:) ::     BENDOC,   & !dissolved organic carbon flux (LDOC) (gC/m^2/day)
                                                    BENNH4,   &    !gN/m^2/day
                                                    BENNO3,   &    !gN/m^2/day
                                                    BENDON,   &    !gN/m^2/day
                                                    BENDOP,   & !gP/m^2/day
                                                    BENPO4,   & !gP/m^2/day
                                                    BENCOD,   & !gO2/m^2/day
                                                    BENDO,    & !gO2/m^2/day (positive adding oxygen to water column)
                                                    BENSA,    & !benthic flux of avaiable dissolved silicate (gSi/m^2/day) (positive into water column)
                                                    BENDEN,   & !bebtguc denitrification flux (gN/m^2/day)
                                                    BENCH4G,  &
                                                    BENCH4A
             
             INTEGER,ALLOCATABLE,DIMENSION(:) ::     SBN    &
                                                    ,BBN    &
                                                    ,HMSBV  
            !Wen Long                               ,SFA      !  Wen Long deprecated SFA, replaced by ART1

             REAL(SP),ALLOCATABLE,DIMENSION(:) ::     WSSNET,  & !settling rate of suspended solids (m/day)
                                                    WSLNET,  & !settling rate of LPOM (m/d)
                                                    WSRNET,  & !settling rate of RPOM (m/d)
                                                    WS1NET,  & !settling rate of alg l (m/d)
                                                    WS2NET,  & !settling rate of alg 2 (m/d)
                                                    WS3NET,  & !settling rate of alg 3 (m/d)
                                                    WSUNET     !settling rate of particulate biogenic (unavaiable) silica (m/d)
             !WLong moved these to wqm_sed.F          
             !REAL(SP), ALLOCATABLE,DIMENSION (:)::    VSED,  &
             !                                       VPMIX, &
            !                                       VDMIX


!            !WLong moved these to mod_sed.F
!             REAL(SP), ALLOCATABLE,DIMENSION (:) ::DIAGP,        &!  
!                                                    DIAGN,        &!
!                                                    DIAGC,        &!CFC
!                                                    DIAGS         !,&
!                                                !MTVEL            !WL moved this to mod_sed

            !Wen Long moved these to wqm_sed.F
            ! REAL(SP),ALLOCATABLE,DIMENSION(:) ::     WSSBNET,  &
            !                                          WSLBNET,  &
            !                                          WSRBNET,  &
            !                                          WS1BNET,  &
            !                                          WS2BNET,  &
            !                                          WS3BNET,  &
            !                                          WSUBNET

             REAL(SP),ALLOCATABLE,DIMENSION(:) ::    PPFWS,  &    !POP flux at bottom of water column (gP/m^2/day)
                                                    PNFWS,  &   !PON flux at bottom of water column (gN/m^2/day)
                                                    PCFWS,  &   !POC flux at bottom of water column (gC/m^2/day)
                                                    PSFWS,  &    !Particulate biogenic silicate flux (gSi/m^2/day)
                                                    SSFWS       !suspended solids flux (gSolids/m^2/dat)

             CHARACTER(LEN=72),DIMENSION(NFLP) ::     METFN,     &
                                                    S1FN,    &
                                                    S2FN,    &
                                                    S3FN,    &
                                                    HYDFN,  &
                                                    CBCFN,    &
                                                    ATMFN,    &
                                                    SVIFN,    &
                                                    BFIFN        !Benthic Flux data file or Sediment Diagenesis module input file
                       
             CHARACTER(LEN=72)BAIFN  !benthic algae input file name
             
             CHARACTER(LEN=72)DFIFN  !deposition feeder input file name
             


             !
             ! moved to sediment module
             !
             !REAL,DIMENSION(3) :: FRPALG1,    & !Algal group 1 Phosphorus fraction for
             !                                  !G1,G2,G3
             !                     FRPALG2,    & !Algal group 2 Phosphorus fraction for 
             !                                  !G1,G2,G3
             !                     FRPALG3,    & !Algal group 3 Phosphorus fraction for 
             !                                  !G1,G2,G3
             !                     FRNALG1,    & !Algal group 1 Nitrogen fractions for
             !                                  !G1,G2,G3
             !                     FRNALG2,    & !Algal group 2 Nitrogen fractions for 
             !                                  !G1,G2,G3
             !                     FRNALG3,    & !Algal group 3 Nitrogen fractions for G1
             !                                  !G2 G3
             !                     FRCALG1,    & !Algal group 1 Carbon fractions for G1,
             !                                  !G2, G3
             !                     FRCALG2,    & !Algal group 2 Carbon fractions for G1,G2
             !                                  !G3
             !                     FRCALG3,    & !Algal group 3 Carbon fractions for G1,
             !                                  !G2, G3
             !                     FRCPHB,    & !Fraction of POC generated for G1,G2,G3
             !                                  !due to predation on benthic algae
             !                     FRNPHB,    & !Fraction of PON genereated for G1,G2,G3
             !                                  !due to predation on benthic algae
             !                     FRPPHB       !Fraction of POP generated for G1,G2,G3
             !                                  !due to predation on benthic algae

             !!moved to wqm_sed.F
             !REAL,ALLOCATABLE,DIMENSION(:,:) :: FRPOP, & !Fraction of non-algal particulate organic phosphorus (dimensionless)
             !                                   FRPON, & !Fraction of non-algal particulate organic nitrogen (dimensionless)
             !                                   FRPOC    !Fraction of non-algal particulate organic carbon (dimensionless)


             REAL(SP),ALLOCATABLE,DIMENSION(:,:) ::& !KLDC,    &  !labile dissolved organic carbin mineralization rate (1/day) ! B CLark moved to WC_DOM
                                               ! KRDC,    &  !refractory dissolved organic carbin mineralization rate (1/day)
                                               ! KLDN,    &  !labile dissolved organic nitrogen mineralization rate (1/day)
                                               ! KRDN,    &  !refractory dissolved organic nitrogen mineralization rate (1/day)
                                               ! KLDP,    &  !labile dissolved organic phosphorus mineralization rate (1/day)
                                                !KRDP,    &  !refractory dissolved organic phosphorus mineralization rate (1/day)
                                                KLPC,    &  !labile particulate organic carbon mineralization rate (1/day)
                                                KRPC,    &  !refractory particulate organic carbon mineralization rate (1/day)
                                                KLPN,    &  !labile particulate organic nitrogen remin rate (1/day)
                                                KRPN,    &  !refractory particulate organic nitrogen remin rate (1/day)
                                                KLPP,    &  !labile particulate organic phosphorus remin rate (1/day)
                                                KRPP,    &  !refractory particulate organic phosphorus remin rate (1/day)
                                                KDCALG,  &  !constant relating DOC respiration rate to total algae biomass C (1/day)/(gC/m^3)
                                                KDNALG,  &  !constant relating DON respiration rate to total algae biomass N (1/day)/(gN/m^3)
                                                KDPALG,  &  !constant relating DOP respiration rate to total algae biomass P (1/day)/(gP/m^3)
                                                KLCALG,  &  !constant relating LPOC hydrolysis rate to total algae biomass C (1/day)/(gC/m^3)
                                                KLNALG,  &  !constant relating LPON hydrolysis rate to total algae biomass N (1/day)/(gN/m^3)
                                                KLPALG,  &  !constant relating LPOP hydrolysis rate to total algae biomass P (1/day)/(gP/m^3)
                                                KSUA,    &  !particulate silica dissolution rate (1/day)
                                                KCOD,    &  !oxidation rate of chemicao oxygen demand (1/day)
                                                KRCOAG,  &  !dissolved organic carbon coagulation rate (1/day)
                                                            ! Coagulation rate increases as salinity increases
                                                NTM         !maximum nitrication nrate at optimal temperature (gN/m^3/day from mannual) !Laura: should be 1/d
                                                            ! nitrification is limited by NH4 and DOXG availability

            REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: ALGDOC,  &      !Dissolvoed organic carbon produced due to algae (all groups)
                                                                !mortality (respiration, metabolism and predation)
                                                                !labile and refractory DOC's together (gC/m^3/day)
												ALGDOC_GL,   &   ! Global algal production of DOC, B Clark July 2018
												ALGDOC_OUT, &!
                                                ALGPOC,  &      !Particulate organic carbon produced due to algae (all groups)
                                                                !mortaility (respiration, metabolism and predation)
                                                                !(gC/m^3/day)
												ALGPOC_OUT, &!
												ALGPOC_GL,   &   ! Global algal production of DOC, B Clark July 2018
												
																								
												ALGPredDOC,  &  ! Added by B Clark to output the algal predation produced DOC
												ALGPredDOC_GL, & !
												ALGPPDOC_GL,  & ! B clark added to output the Primary Production DOC exudation 
												ALGPPdOC, & !
												ALGPredDON,  &  ! Added by B Clark to output the algal predation produced DON
												ALGPredDON_GL, & !
												ALGPPDON_GL,  &  ! same as above
												ALGPPdON, & !
												
                                                ALGDON,  &      !Dissolvoed organic nitrogen produced due to algae g N m^-3 d^-1
												ALGDON_OUT, &!
												ALGDON_GL,   &  ! Globala algal production of DON, added by B Clark July 2018 for output
												
                                                ALGPON,  &      !Particulate organic nitrogen produced due to algae (all groups)
                                                                   !mortaility (respiration, metabolism and predation) 
												ALGPON_OUT, &!
												
												ALGPON_GL,   &  ! Globala algal production of PON, added by B Clark July 2018 for output
                                                                   !(gN/m^3/day) 
												
                                                ALGDOP,  &      !Dissolvoed organic phosphorus produced due to algae (all groups)
                                                                   !mortality (respiration, metabolism and predation)
                                                                   !labile and refractory DON's together
                                                                   !(gP/m^3/day)
                                                ALGPOP,  &      !Particulate organic phosphorus produced due to algae (all groups)
                                                                   !mortality (respiration, metabolism and predation)
                                                                   !(gP/m^3/day)
                                                ALGNH4,  &      !Soure+sink term for water column NH4 due to algae uptake (gN/m^3/day)
                                                                !   positive as source of NH4		
												ALGNH4_OUT, & !
												ALGNH4_GL,  & ! Global delta NH4 for algae
												
                                                ALGNO3,  &      !Source+sink term for water column NO3 due to algae uptake (gN/m^3/day)
                                                                !   positive as source of NO3
												ALGNO3_OUT,& !
												ALGNO3_GL ,   & ! Global delta NO3 for algae
												
                                                ALGPO4,  &      !Source+sink term for water column PO4 due to algae uptake (gP/m^3/day)
                                                                !   positive as source of PO4
  !                                             MNLLDOC, &      !Mineralization rate of labile DOC (positive losing LDOC)     ! B Clark moved all mineralization rates to wc_dom Sep 2015
  !                                                              !  gC/m^3/day
 !                                               MNLRDOC, &      !Mineralization rate of refractory DOC (positive losing RDOC)
 !                                                              !  gC/m^3/day
 !                                               MNLLDON, &      !Mineralization rate of labile DON (positive losing LDON)
 !                                                               !  gN/m^3/day                                        
  !                                              MNLRDON, &      !Mineralization rate of refractory DON (positive losing RDON)
 !                                                               !  gN/m^3/day                                        
  !                                              MNLLDOP, &      !Mineralization rate of labile DOP (positive losing LDOP)
 !                                                               !  gP/m^3/day
 !                                               MNLRDOP, &      !Mineralization rate of refractory DOC (positive losing RDOP)
!                                                                !  gP/m^3/day
												
                                                HDRLPOC, &      !Hydrolysis loss rate of labile POC (positive losing LPOC)
                                                                !   gaining LDOC (gC/m^3/day)
												HDRLPOC_OUT, &!
												HDRLPOC_GL,  & ! Global hydrolysis of LPOC
												
                                                HDRRPOC, &      !Hydrolysis loss rate of refractory POC (positive losing RPOC)
                                                                !   gaining RDOC (gC/m^3/day)
												HDRRPOC_OUT, &!
												HDRRPOC_GL,  & ! Global hydrolysis of LPOC
												
                                                HDRLPON, &      !Hydrolysis loss rate of labile PON (positive losing LPON)
                                                                !   gaining LDON (gN/m^3/day)
												HDRLPON_OUT, &!
												HDRLPON_GL,  & ! Global hydrolysis of LPON
												
                                                HDRRPON, &      !Hydrolysis loss rate of refractory PON (positive losing RPON)
                                                                !   gaining RDOC (gN/m^3/day)	
												HDRRPON_OUT, &!
												HDRRPON_GL,  &  ! Global Hydrolysis of RPON
												
                                                HDRLPOP, &      !Hydrolysis loss rate of labile POP (positive losing LPOP)
                                                                !   gaining LDOP (gP/m^3/day)
                                                HDRRPOP, &      !Hydrolysis loss rate of refractory POP (positive losing RPOP)
                                                                !   gaining LDOC (gP/m^3/day)
												
                                                DENNO3,  &      !Denitrification loss rate of NO3 (gN/m^3/day)
                                                                !    Positive losing NO3
												DENNO3_OUT, &!
												
												DENNO3_GL,  &! Global denitrification loss of NO3
												
                                                DENIT,   &      !Denitrication loss rate ot DOC (gC/m^3/day)
                                                                !    positive losing DOC (nitrate limited reaction)
												DENIT_OUT, &!
												DENIT_GL,  & ! Global denitrification loss of DOC
												
                                                DENITN,   &      !Denitrication loss rate ot DON (gN/m^3/day)
                                                                !    positive losing DON (nitrate limited reaction)	
												DENITN_OUT, &!
												
												DENITN_GL,  & !Global loss rate of DON from denitrification
												
                                                DENITP,   &      !Denitrication loss rate ot DOP (gC/m^3/day)
                                                                !    positive losing DOP (nitrate limited reaction)												
                                                NFIX            !Nitrogen fixation rate (gN/m^3/day)


             
             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: PSD,           &
                                                SAP,        &
                                                ALGUP,        &
                                                ALGRES

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: DORALG,        &
                                                DOPR,        &
                                                DCOD,        &
                                                DDOC,        &
                                                NITRIF

             REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: A_T,         &
                                                AP1,         &
                                                ABM1,         &
                                                APR1,         &
                                                AP2,         &
                                                ABM2,         &
                                                APR2,         &
                                                AP3,         &
                                                ABM3,         &
                                                APR3,       &
                                                AALGDOC,     &
                                                AALGPOC,     &
                                                ADENIT,     &
                                                AMNLDOC,    &
                                                AHDRPOC,    &
                                                AALGNH4,    &
                                                AALGNO3,    &
                                                AALGDON,    &
                                                AALGPON,    &
                                                ANT,         &
                                                ADENNO3,    &
                                                AMNLDON,    &
                                                AHDRPON,    &
                                                ANFIX,      &  !WLong inconsistent with ACFIX and not used
                                                AALGPO4,    &
                                                AALGDOP,    &
                                                AALGPOP,    &
                                                AMNLDOP,    &
                                                AHDRPOP,    &
                                                APSD,        &
                                                ASAP,        &
                                                AALGUP,        &
                                                AALGRES,    &
                                                ADO,        &
                                                ADORALG,    &
                                                ADOPR,        &
                                                ADCOD,        &
                                                ADDOC,        &
                                                ANITRIF    

            !WLong moved to mod_zoop.F
            !REAL,DIMENSION(366) :: TVPRSZ,            &
            !                        TVPRLZ
                                    
             REAL(SP),DIMENSION(366) :: TVPR,            &        !A multiplier to predation rate
                                    B2GR

             REAL(SP),ALLOCATABLE ::     QLIT(:,:),       &
                                    CONLIT(:,:,:) 

            !RGL ADDED LENGTH2 AND CTR
                                    
             INTEGER,SAVE ::      LENGTH,        & 
                                IJDAY,        &
                                LENGTH2,     &
                                IWRIT 
                                
             INTEGER,SAVE ::      NAPL,         &
                                NDLT,         &
                                NSNP,         &
                                NPLT,         &
                                NTFL
                                
             INTEGER,SAVE ::      NKFL,         &
                                NOPL,         &
                                NOINT,         &
                                NMBL
                                
             INTEGER,SAVE ::      KFLDP,         &
                                MBLDP,         &
                                INFLOW
                                
             INTEGER,SAVE ::      NDIA,         &
                                NRSO,         &
                                NHYDF,         &
                                NTVDF,         &
                                KEIPTR!, NUMBOX

             INTEGER,SAVE ::      JCB,         &
                                IT1,         &
                                IT2,         &
                                IT3,        & 
                                JT1,        & 
                                IVTMP,         &
                                JVTMP        
                                
             INTEGER,SAVE ::      JCS1MAX,     &
                                JCS2MAX,     &
                                JCS3MAX,     &
                                JDAYVB

             INTEGER ::     BAOPTR 
             
             REAL(SP),SAVE ::      TH,         &
                            TMSTRT,            &
                            XYDF!,            &
!                            CTEMPI            !WLong moved CTEMPI to mod_sed.F
                            
!             REAL(SP),SAVE ::      BENSTI        !WLong moved to mod_sed.F
             
             
             REAL(SP), SAVE::     BBMI,            &
                                DFEEDI,            &
                                DUMMY
                            
             REAL(SP),SAVE ::      ELTMSPLT,        &
                                ELTMSTFL,        &
                                ELTMSKFL,        &
                                EPI!,             &                            
!                                HST2I,             &    !WLong moved to mod_sed.F
!                                PO4T2I
!                                CPOSI,            &
!                                CH4T2I,            &
!                                CH41TI            
!                                SIT2I,             &
!                                SO4T2I
                            
             REAL(SP),SAVE ::      RTI,             &
                                SFATOT,         &
                                SHI!,             &
!                                SUMI             !WLong moved this to mod_owq.F
                            
             REAL(SP),SAVE ::      RATIO,             &
                                RELTMS,            &
                                TODS,             &
                                DSSR,             &
                                SAVENT,         &
                                SEDENT,         &
                                SFENT            
                            
             REAL(SP),SAVE ::      SFATMP,         &
                                ALFAS,            &
                                CONC1,             &
                                CONC2,             &
                                CONC3,             &
                                CADV,             &
                                OMTH             
                            
             REAL(SP),SAVE ::      BNDC,            &
                                BNDN,             &
                                BNDP
                            
                            
             REAL(SP),SAVE ::      DOS,             &
                                DLTAV,             &
                                ELTMJD,         &
                                HOURVB,         &
                                HMTV2
                            
             REAL(SP),SAVE ::      WCMN,             &
                                WCMP,             &
                                WCMC,             &
                                WCMS,             &
                                WQMDLV,         &
                                WQMTV2            
                            
             REAL(SP),SAVE ::      DLWCMC,            &
                                DLWCMN,         &
                                DLWCMP,         &
                                DLWCSSMC,         &
                                DLWCSSMN,         &
                                DLWCSSMP
                            
             REAL(SP),SAVE ::    DTMP1
                                         
             REAL(SP),SAVE ::      AVGINT,         &
                                DTMP2,             &
                                DFTMP  
                            
             REAL(SP),SAVE ::      DLSEDMC,         &!kgC
                                DLSEDMN,         &!kgN
                                DLSEDMP,         &!kgP
                                DLSEDSSMC,         &
                                DLSEDSSMN,         &
                                DLSEDSSMP 
                            
            REAL(SP),SAVE ::      ERRWCC,         &
                                ERRWCN,         &
                                ERRWCP,            &
                                ERRSEDC,        &
                                ERRSEDN,         &
                                ERRSEDP
                            
            REAL(SP), SAVE::    ELTMMBL


        !***** Equivalence declarations

        !     EQUIVALENCE (C1(0,1),    CSTAR(0,1))
        !     EQUIVALENCE (DTT(1),     DTC(1,1)),                                &
        !                 (DTSSI(1),   DTC(1,3)),    (DTB1(1),     DTC(1,4)),    &
        !                 (DTB2(1),    DTC(1,5)),    (DTB3(1),     DTC(1,6)),    &
        !                 (DTSZ(1),    DTC(1,7)),    (DTLZ(1),     DTC(1,8)),    &
        !                 (DTLDOC(1),  DTC(1,9)),    (DTRDOC(1),   DTC(1,10)),   &
        !                 (DTLPOC(1),  DTC(1,11)),   (DTRPOC(1),   DTC(1,12)),   &
        !                 (DTNH4(1),   DTC(1,13)),   (DTNO3(1),    DTC(1,14)),   &
        !                 (DTUREA(1),  DTC(1,15)),   (DTLDON(1),   DTC(1,16)),   &
        !                 (DTRDON(1),  DTC(1,17)),   (DTLPON(1),   DTC(1,18)),   &
        !                 (DTRPON(1),  DTC(1,19)),   (DTPO4(1),    DTC(1,20)),   &
        !                 (DTLDOP(1),  DTC(1,21)),   (DTRDOP(1),   DTC(1,22)),   &
        !                 (DTLPOP(1),  DTC(1,23)),   (DTRPOP(1),   DTC(1,24)),   &
        !                 (DTPIP(1),   DTC(1,25)),   (DTCOD(1),    DTC(1,26)),   &
        !                 (DTDO(1),    DTC(1,27)),   (DTSIUPB(1),   DTC(1,28)),   &
        !                 (DTSIAT(1),    DTC(1,29)),   (DTPIB1(1),  DTC(1,30)),   &
        !                 (DTPIB2(1),  DTC(1,31)),   (DTPIB3(1),   DTC(1,32))


        !     EQUIVALENCE (DTC(0,1), DTM(0,1))

        !Water column fluxes
        !     EQUIVALENCE (FLXTTEM(1), FLUXT(1,1)),  (FLXTSAL(1),  FLUXT(1,2)),  &
        !                 (FLXTSSI(1), FLUXT(1,3)),  (FLXT1(1),    FLUXT(1,4)),  &
        !                 (FLXT2(1),   FLUXT(1,5)),  (FLXT3(1),    FLUXT(1,6)),  &
        !                 (FLXTSZ(1),  FLUXT(1,7)),  (FLXTLZ(1),   FLUXT(1,8)),  &
        !                 (FLXTLDOC(1),FLUXT(1,9)),  (FLXTRDOC(1), FLUXT(1,10)), &
        !                 (FLXTLPOC(1),FLUXT(1,11)), (FLXTRPOC(1), FLUXT(1,12)), &
        !                 (FLXTNH4(1), FLUXT(1,13)), (FLXTNO3(1),  FLUXT(1,14)), &
        !                 (FLXTUREA(1),FLUXT(1,15)), (FLXTLDON(1), FLUXT(1,16)), &
        !                 (FLXTRDON(1),FLUXT(1,17)), (FLXTLPON(1), FLUXT(1,18)), &
        !                 (FLXTRPON(1),FLUXT(1,19)), (FLXTPO4(1),  FLUXT(1,20)), &
        !                 (FLXTLDOP(1),FLUXT(1,21)), (FLXTRDOP(1), FLUXT(1,22)), &
        !                 (FLXTLPOP(1),FLUXT(1,23)), (FLXTRPOP(1), FLUXT(1,24)), &
        !                 (FLXTPIP(1), FLUXT(1,25)), (FLXTCOD(1),  FLUXT(1,26)), &
        !                 (FLXTDO(1),  FLUXT(1,27)), (FLXTSU(1),   FLUXT(1,28)), &
        !                 (FLXTSA(1),  FLUXT(1,29)), (FLXTPIB1(1), FLUXT(1,30)), &
        !                 (FLXTPIB2(1),FLUXT(1,31)), (FLXTPIB3(1), FLUXT(1,32)) 

        !fluxes that affect sediment 
        !     EQUIVALENCE (FLXSSSI(1), FLUXS(1,1)),  (FLXS1(1),    FLUXS(1,2)),  &
        !                 (FLXS2(1),   FLUXS(1,3)),  (FLXS3(1),    FLUXS(1,4)),  &
        !                 (FLXSPOC(1), FLUXS(1,5)),  (FLXSPON(1),  FLUXS(1,6)),  &
        !                 (FLXSPO4(1), FLUXS(1,7)),  (FLXSPOP(1),  FLUXS(1,8)),  &
        !                 (FLXSPIP(1), FLUXS(1,9)),  (FLXSSI(1),   FLUXS(1,10)), &
        !                 (FLXSPIB1(1),FLUXS(1,11)),                             &
        !                 (FLXSPIB2(1),FLUXS(1,12)), (FLXSPIB3(1), FLUXS(1,13))
             
        !     EQUIVALENCE (AFLXPOC(1), AFLUX(1,1)),  (AFLXDOC(1),  AFLUX(1,2)),  &
        !                 (AFLXTC(1),  AFLUX(1,3)),  (AFLXDON(1),  AFLUX(1,4)),  &
        !                 (AFLXDIN(1), AFLUX(1,5)),  (AFLXPON(1),  AFLUX(1,6)),  &
        !                 (AFLXTN(1),  AFLUX(1,7)),  (AFLXDOP(1),  AFLUX(1,8)),  &
        !                 (AFLXDIP(1), AFLUX(1,9)),  (AFLXPOP(1),  AFLUX(1,10)), &
        !                 (AFLXPIP(1), AFLUX(1,11)), (AFLXTP(1),   AFLUX(1,12)), &
        !                 (AFLXSSI(1), AFLUX(1,13))
        !     EQUIVALENCE (BENDOC(1),  BFLUX(1,1)),  (BENNH4(1),   BFLUX(1,2)),  &
        !                 (BENNO3(1),  BFLUX(1,3)),  (BENDON(1),   BFLUX(1,4)),  &
        !                 (BENPO4(1),  BFLUX(1,5)),  (BENDOP(1),   BFLUX(1,6)),  &
        !                 (BENCOD(1),  BFLUX(1,7)),  (BENDO(1),    BFLUX(1,8)),  &
        !                 (BENSA(1),   BFLUX(1,9)),  (BENCH4A(1),  BFLUX(1,10)), &
        !                 (BENCH4G(1), BFLUX(1,11))
        !     EQUIVALENCE (BENDOCB(1), BFLUXB(1,1)), (BENNH4B(1),  BFLUXB(1,2)), &
        !                 (BENNO3B(1), BFLUXB(1,3)), (BENPO4B(1),  BFLUXB(1,4)), &
        !                 (BENCODB(1), BFLUXB(1,5)), (BENDOB(1),   BFLUXB(1,6)), &
        !                 (BENSAB(1),  BFLUXB(1,7))

        !***** Common declarations

             REAL(SP) ::     TMEND,        &
                        OLDNXCBC,    &
                        TE,            &
                        WMS        
                        
             REAL(SP) ::     AOCR,        &
                        AONT,        &
                        FDOP,        &
                        AREAR,        &
                        BREAR,        &
                        CREAR
                        
             INTEGER  ::     NAC,    &
                            NB,        &
                            NQF,    &
                            NL,        &
                            NSB,    &
                            NBB,    &!
                            NHQF,    &!
                            NSQF
                            
             REAL(SP) ::     TRCOD
             
             REAL(SP) ::     ASC1,        &
                        ASC2,        &
                        ASC3,        &
                        FSAP,        &
                        TRSUA
                        
             REAL(SP) :: FR,SCTOX
             REAL(SP) :: PRSP1,  &  !Fraction of production consumed by predation for Algae 1
                     PRSP2,  &  !Fraction of production consumed by predation for Algae 2
                     PRSP3,  &  !Fraction of production consumed by predation for Algae 3
                     STF1,   &  !Salinity toxicity coefficient for algae 1
                     STF2,   &  !Salinity toxicity coefficient for algae 2
                     STF3,   &  !Salinity toxicity coefficient for algae 3
                     VMAX1,  &  !maximum P:C ratio for algae 1? (gP/gC)
                     VMAX2,  &  !maximum P:C ratio for algae 2? (gP/gC)
                     VMAX3      !maximum P:C ratio for algae 3? (gP/gC)
                     
             REAL(SP) ::     ANC1,    &
                        ANC2,    &
                        ANC3,    &
                        FNI1,    &
                        FNI2,    &
                        FNI3,    &                        
                        FNLP1,    &
                        FNLP2,    &
                        FNLP3,    &
                        FNRP1,    &
                        FNRP2,    &
                        FNRP3,  &
                        FNIP,    &
                        FNLPP,    &
                        FNRPP,    &
                        FNUP    
                        ! FNLDP,    &
                      !  FNRDP
						
						
             REAL(SP) ::     ALPHMIN1,    &
                        ALPHMIN2,    &
                        ALPHMIN3,    &
                        ALPHRAT1,    &
                        ALPHRAT2,    &
                        ALPHRAT3,      &    
                        CHLCMN1,    &
                        CHLCMN2,    &
                        CHLCMN3,    &
                        ACHLC1,        &
                        ACHLC2,        &
                        ACHLC3,        &
                        BCHLC1,        &
                        BCHLC2,        &
                        BCHLC3,        &
                        CCHLC1,        &
                        CCHLC2,        &
                        CCHLC3
                        
             REAL(SP) :: ANDC
             REAL(SP):: Q01, &  !Previously called APC1, Algal group 1 Phosphorus-to-carbon ratio (gP/gC)
                     Q02, &     !Previously called APC2, Algal group 2 Phosphorus-to-carbon ratio (gP/gC)
                     Q03      !Previously called APC3, Algal group 3 Phosphorus-to-carbon ratio (gP/gC)
                     
            REAL(SP) :: FPI1,        &   ! 
                    FPI2,        &
                    FPI3,        &
                    FPLP1,        &
                    FPLP2,        &
                    FPLP3,        &
                    FPRP1,        &
                    FPRP2,        &
                    FPRP3,        &
                    FPIP,          &
				    FPLPP,        &
                    FPRPP    
                    !FPLDP,        &
                    !FPRDP,        &
    
             !       
            ! REAL(SP) :: FCD1 !FCLD1,   &     ! USED IN MOD_WC_DOM FOR ALGAL CONTRIBUTION TO DOM POOLS, renamed to FCD1 because DOM is fractionated to the different pools after it is produced by algae
            !             FCD2!,        FCLD2! &   
			!			! FCD3
            !    !    FCLD3,        &
            !     !   FCRD1,        &
            !    !    FCRD2,        &
            !     !   FCRD3,  
			!		 
			!REAL(SP):: FPD1,        &    ! Same as above for Phosphorus
            !           FPD2!,        &!
			!		  ! FPD3
            !           !FPLD3,        &
            !           !FPRD1,         &
            !           !FPRD2,        &
            !           !FPRD3,     
			!		
			!REAL(SP):: FND1,    &  ! same as above for Nitrogen
            !            FND2  !,    &!
			!			!FND3
            !     !      FNLD3,  &
            !    !        FNRD1,    &
            !     !       FNRD2,    &
            !      !      FNRD3    
          REAL(SP) ::	 FCDP ,&      		! Fraction of dissolved carbon production from algal predation, same for all algal groups
						 FNDP, &
						 FPDP
             REAL(SP) ::       FCLP1,        &
                    FCLP2,        &
                    FCLP3,      &
                    FCRP1,        &
                    FCRP2,        &
                    FCRP3,        &
                   
                    FCLPP,        &
                    FCRPP,        &
                    TRMNL,        &
                    TRHDR,      &
                    AANOX        

             REAL(SP) ::     FD,            &
                        TTSS,        &
                        SREFLECT,    &
                        PRECIP
                        
            REAL(SP) ::        ATMNH4,        &
                        ATMNO3,        &
                        ATMLDON,    &
                        ATMRDON,      &
                        ATMPO4,        &
                        ATMLDOP,    &
                        ATMRDOP 

             REAL(SP) ::     REDS1C,        &
                        REDS1N,        &
                        REDS1P,        &
                        REDS2C,        &
                        REDS2N,        &
                        REDS2P,        &
                        REDCBC,     &
                        REDCBN,        &
                        REDCBP,        &
                        REDS3C,        &
                        REDS3N,        &
                        REDS3P       

             REAL(SP) ::     PS1,        &
                        PS2,        &
                        PS3
                        
             REAL(SP) ::     Z1TP,        &
                        Z2TP
        
             REAL(SP) ::     DLT,        &
                        AHMDLT,        &
                        FILGTH,        &
                        ZDFMUL,        &
                        ZDFBCK,        &
                        WQ_DT
                        
            !WLong moved to mod_zoop.F            
             !INTEGER  ::     KTBSZ,        &
             !                KTBLZ
            INTEGER::        NHMDLT,        &
                            NWQMR,        &
                            NHMR,        &
                            NIT                    !count of integration time steps


        !Wen Long added the following for debugging algae production and mortality

             REAL(SP),ALLOCATABLE :: P1_GL(:,:),  &  !gross production of ALG1 (1/day)
                                 P2_GL(:,:),  &  !gross production of ALG2 (1/day)
                                 BM1_GL(:,:), &  !basal metabolism of ALG1 (1/day)
                                 BM2_GL(:,:), &  !                       2
                                 NL1_GL(:,:), &  !Nitrogen limitation of ALG1 (non-dimensional)
                                 NL2_GL(:,:), &  !                          2
                                 PL1_GL(:,:), &  !Phosphorus                1
                                 PL2_GL(:,:), &  !                          2
                                 FI1_GL(:,:), &  !Light                     1
                                 FI2_GL(:,:), &  !                          2
                                 B1SZ_GL(:,:),&  !micro zoop consumption (mgC/L/day) for ALG1
                                 B1LZ_GL(:,:),&  !macro                                  ALG1
                                 B2SZ_GL(:,:),&  !micro                                  ALG2
                                 B2LZ_GL(:,:),&  !macro                                  ALG2
                                 PR1_GL(:,:), &  !predation loss of ALG1  (mgC/L/day)     
                                 PR2_GL(:,:), &  !                  ALG2  (mgC/L/day)
								 
                                 IAVG_GL(:,:),&  !light intensity (E/m^2/day)
								 KDLAMB_int_GL(:,:) ,& !integrated Kd for PAR (m^-1)
								 SSI_GL(:,:),  & ! Added By B Clark June 2018 for outputting ssi
                                 IK1_GL(:,:),  &  !saturation light intensity (E/m^2/day) based on PI curve initial slope
                                 IK2_GL(:,:),  &  !saturation light intensity (E/m^2/day) based on PI curve initial slope
                                 IK1(:,:),    &  ! local IK for alg 1
                                 IK2(:,:)        ! local IK for alg 2
								 

        !Wen Long added the following for debugging DOXG source and sinks in water column

             REAL(SP), ALLOCATABLE, DIMENSION (:,:) :: & 
                                  DDOC_GL,       &   !oxygen loss due to remin of DOC (LDOC and RDOC)          (mgO2/L/day)  (sink)
                                  DCOD_GL,       &   !oxygen loss due to oxidation of COD in water column      (mgO2/L/day)  (sink)
                                  NITRIF_GL,     &   !oxygen loss due to nitrification                         (mgO2/L/day)  (sink)
                                  DOSZ_GL,       &   !oxygen loss due to microzooplankton respiration          (mgO2/L/day)  (sink)
                                  DOLZ_GL,       &   !oxygen loss due to macrozooplankton respiration          (mgO2/L/day)  (sink)
                                  DOPR_GL,       &   !sum of DOP1 and DOP2 (total loss due to predation)       (mgO2/L/day)  (sink)
                                  DORALG_GL,     &   !sum of DOR1 and DOR2 (total production)                  (mgO2/L/day)  (source)
                                  RESP_GL,       &   !total oxygen loss due to respiration (ALGDO, DOPR, DDOC, DCOD, NITRIF, DOSZ, DOLZ)  (mgO2/L/day)
                                  FTCOD_GL,      &   !                  where ALGDO is oxygen loss due to algal resipration and basal metabolism  (sink)
                                    COD_GL,      &   !COD decay rate after applying temperature control (m/day)
                                  REAERDO_GL,     &      !COD (mgO2-equiv/L)
                                                    !areation of DOXG (mg-O2/L/day)  (global)
                                  REAERDO            !local variable of DOXG areation


! B Clark for sediment Resuspension, June 2018
	REAL(SP), ALLOCATABLE :: ShearStress(:), RHO1(:,:), RHO1_GL(:,:),ShearStress_GL(:) ! Added by B Clark to calculate shear stress for sediment resuspension, June 2018, Pascals and kg m^3

    
        !WLong moved these to mod_sed.F
        !     REAL(SP) :: SEDMN,      & !total sediment organic nitrogen  (kgN)
        !             SEDMP,      & !total sediment organic phosphorus (kgP)
        !             SEDMC         !total sediment organic carbon  (kgC)

             REAL(SP) :: S1FLXN,     & !Source one load of N flux (gN/sec or gN/day ?)
                     S1FLXP,     & !Source one load of P flux (gP/sec or gP/day ?)
                     S1FLXC,     & !Source one load of C flux (gC/sec or gC/day ?)
                     S3FLXN,     & !Source three load of N flux (gN/sec or gN/day ?)
                     S3FLXP,     & !Source three load of P flux (gP/sec or gN/day ?)
                     S3FLXC,     & !source three load of C flux (gC/sec or gN/day ?)
                     S2FLXN,     & !source two load N (gN/day)
                     S2FLXP,     & !source two load P (gP/day)
                     S2FLXC       !source two load C (gC/day)

             REAL(SP) :: ATMFLXN,    & !Atmospheric nitrogen flux (kgN) in DLT (sec) time
                     ATMFLXP       !Atmospheric phosphorus flux  (kgP) in DLT (sec) time
             
             REAL(SP) :: BENFLXPN,   &  !benthic particulate organic nitrogen flux  (kgN)
                                    !positive into water column, leaving sediments
                     BENFLXPP,   &  !benthic particulate organic phosphorus flux (kgP)
                     BENFLXPC,   &  !benthic particulate organic Carbon flx (kgC)
                     BENFLXDN,   &  !benthic dissolved nitrogen (NH4+NO3) flux (kgN)
                     BENFLXDP,   &  !benthic dissolved phosphorus (PO4) flux (kgP) 
                     DLWCKMN,    &  !nitogen change (kgN) due to denitrification kinetics in water column
                     DLWCKMC,    &  !carbon change (kgC) due to kinetics in water column
                     BNDFLXN,    &  !latetal boundary nitrogen flux due to flow (advection) and diffusion
                                    !in water column  (kgN)
                     BNDFLXP,    &  !lateral boundary phosphorus flux due to flow and diffusion
                                    !in water colum (kgP)
                     BNDFLXC,    &  !lateral boundary carbon flux due to flow (advction)
                                    !and diffusion  in water column (kgC)
                     DLSEDKN,    &  !nitrogen change (kgN) due to denitrification in
                                    !sediments (layer 1 and layer 2)
                     DLSEDKC,    &  !carbon change (kgC) due to sediment carbon
                                    !diagenesis
                     BURIALFLXN, &  !nitrogen change (kgN) due to sediment burial flux 
                                    !positive into sediments
                     BURIALFLXP, &  !phosphorus change (kgP) due to sediment burial flux
                                    !positive into sediments
                     BURIALFLXC     !carbon change (kgC) due to sediment burial flux
                                    !positive into sediments

             REAL(SP)     :: IWCMS, &  !total initial silica in water column
                         IWCMN, &  !total initial nitrogen in water column 
                         IWCMP, &  !total initial phosphorus in water column 
                         IWCMC     !total initial carbon in water column 
                       
        CONTAINS
           
           !Subroutine WQM_ALLOC()
           !Subroutine WQM_DEALLOC()
           
        SUBROUTINE WQM_ALLOC

           IMPLICIT NONE

           ALLOCATE(T(0:MTLOC,KBM1));     T = 0.0
           ALLOCATE(SSI(0:MTLOC,KBM1));   SSI = 0.0
           ALLOCATE(B1(0:MTLOC,KBM1));    B1  = 0.0
           ALLOCATE(B2(0:MTLOC,KBM1));    B2  = 0.0
           ALLOCATE(B3(0:MTLOC,KBM1));    B3  = 0.0
           !ALLOCATE(LDOC(0:MTLOC,KBM1));  LDOC = 0.0
           ALLOCATE(LPOC(0:MTLOC,KBM1));  LPOC = 0.0
           ALLOCATE(RPOC(0:MTLOC,KBM1));  RPOC = 0.0
	
		   !Below added by B Clark Feb 2019 For Budget Calculation
		   ALLOCATE(DTRPOC_out(MTLOC,KBM1)); DTRPOC_out=0.
		   ALLOCATE(DTLPOC_out(MTLOC,KBM1)); DTLPOC_out=0.
		   ALLOCATE(DTRPON_out(MTLOC,KBM1)); DTRPON_out=0.
		   ALLOCATE(DTLPON_out(MTLOC,KBM1)); DTLPON_out=0.
		   
           ALLOCATE(NH4(0:MTLOC,KBM1));   NH4  = 0.0
           ALLOCATE(NO3(0:MTLOC,KBM1));   NO3  = 0.0

           ALLOCATE(DTNH4_out(1:MTLOC,KBM1));   DTNH4_out  = 0.0
           ALLOCATE(DTNO3_out(1:MTLOC,KBM1));   DTNO3_out  = 0.0		   
		   
		   
           !ALLOCATE(LDON(0:MTLOC,KBM1));  LDON = 0.0
           ALLOCATE(LPON(0:MTLOC,KBM1));  LPON = 0.0
           ALLOCATE(RPON(0:MTLOC,KBM1));  RPON = 0.0
           ALLOCATE(PO4(0:MTLOC,KBM1));   PO4  = 0.0
           !ALLOCATE(LDOP(0:MTLOC,KBM1));  LDOP = 0.0
           ALLOCATE(LPOP(0:MTLOC,KBM1));  LPOP = 0.0
           ALLOCATE(RPOP(0:MTLOC,KBM1));  RPOP = 0.0
           ALLOCATE(COD(0:MTLOC,KBM1));   COD  = 0.0
           ALLOCATE(DOXG(0:MTLOC,KBM1));    DOXG   = 0.0
           ALLOCATE(SIUPB(0:MTLOC,KBM1));    SIUPB   = 0.0
           ALLOCATE(SIAT(0:MTLOC,KBM1));    SIAT   = 0.0
           ALLOCATE(SALT(0:MTLOC,KBM1));  SALT = 0.0

           !ALLOCATE(RDOC(0:MTLOC,KBM1));  RDOC = 0.0
           !ALLOCATE(RDON(0:MTLOC,KBM1));  RDON = 0.0
           !ALLOCATE(RDOP(0:MTLOC,KBM1));  RDOP = 0.0
           ALLOCATE(UREA(0:MTLOC,KBM1));  UREA = 0.0
           ALLOCATE(PIP(0:MTLOC,KBM1));   PIP  = 0.0
           ALLOCATE(PIB1(0:MTLOC,KBM1));  PIB1 = 0.0
           ALLOCATE(PIB2(0:MTLOC,KBM1));  PIB2 = 0.0 
           ALLOCATE(PIB3(0:MTLOC,KBM1));  PIB3 = 0.0

           ALLOCATE(PN1(MTLOC,KBM1));            PN1     = 0.0
           ALLOCATE(PN2(MTLOC,KBM1));            PN2     = 0.0
           ALLOCATE(PN3(MTLOC,KBM1));            PN3     = 0.0
           ALLOCATE(BM1(MTLOC,KBM1));            BM1     = 0.0
           ALLOCATE(BM2(MTLOC,KBM1));            BM2     = 0.0
           ALLOCATE(BM3(MTLOC,KBM1));            BM3     = 0.0
    !       ALLOCATE(MNLLDOC(MTLOC,KBM1));        MNLLDOC = 0.0
           ALLOCATE(FTCOD(MTLOC,KBM1));          FTCOD   = 0.0     
           ALLOCATE(NPP(MTLOC,KBM1));            NPP     = 0.0
           ALLOCATE(P1(MTLOC,KBM1));             P1      = 0.0
           ALLOCATE(P2(MTLOC,KBM1));             P2      = 0.0
           ALLOCATE(P3(MTLOC,KBM1));             P3      = 0.0
		   
           ALLOCATE(NT(MTLOC,KBM1));             NT      = 0.0
           ALLOCATE(NT_OUT(MTLOC,KBM1));             NT_OUT      = 0.0		   
		   ALLOCATE(NT_GL(MGL,KBM1));             NT_GL      = 0.0
		   
           ALLOCATE(PR1(MTLOC,KBM1));            PR1     = 0.0
           ALLOCATE(PR2(MTLOC,KBM1));            PR2     = 0.0
           ALLOCATE(PR3(MTLOC,KBM1));            PR3     = 0.0
		   
           ALLOCATE(DENIT(MTLOC,KBM1));          DENIT   = 0.0 
           ALLOCATE(DENIT_OUT(MTLOC,KBM1));          DENIT_OUT   = 0.0 		   
           ALLOCATE(DENIT_GL(MGL,KBM1));          DENIT_GL   = 0.0 
		   
           ALLOCATE(DENITN(MTLOC,KBM1));          DENITN   = 0.0 
           ALLOCATE(DENITN_OUT(MTLOC,KBM1));          DENITN_OUT   = 0.0 		   
           ALLOCATE(DENITN_GL(MGL,KBM1));          DENITN_GL   = 0.0 	
		   
		   ALLOCATE(DENITP(MTLOC,KBM1));          DENITP   = 0.0 
           ALLOCATE(RATOX(MTLOC,KBM1));          RATOX   = 0.0
           ALLOCATE(GPP(MTLOC,KBM1));            GPP     = 0.0
      !     ALLOCATE(MNLRDOC(MTLOC,KBM1));        MNLRDOC = 0.0
           ALLOCATE(P1NNF(MTLOC,KBM1));          P1NNF   = 0.0
           
           ALLOCATE(ASRAT(MTLOC));             ASRAT   = 0.0
           ALLOCATE(total_netPP(MTLOC));       total_netPP = 0.0

           ALLOCATE(KE(MTLOC,KBM1));           KE     = 0.0  
           ALLOCATE(FI1(MTLOC,KBM1));          FI1    = 0.0
           ALLOCATE(FI2(MTLOC,KBM1));          FI2    = 0.0
           ALLOCATE(FI3(MTLOC,KBM1));          FI3    = 0.0
           ALLOCATE(NL1(MTLOC,KBM1));          NL1    = 0.0
           ALLOCATE(NL2(MTLOC,KBM1));          NL2    = 0.0
           ALLOCATE(NL3(MTLOC,KBM1));          NL3    = 0.0
           ALLOCATE(PL1(MTLOC,KBM1));          PL1    = 0.0
           ALLOCATE(PL2(MTLOC,KBM1));          PL2    = 0.0 
           ALLOCATE(PL3(MTLOC,KBM1));          PL3    = 0.0

           ALLOCATE(RESP(MTLOC,KBM1));         RESP   = 0.0

           ALLOCATE(SL1(MTLOC,KBM1));          SL1    = 0.0
           ALLOCATE(SL2(MTLOC,KBM1));          SL2    = 0.0           
           ALLOCATE(SL3(MTLOC,KBM1));          SL3    = 0.0
           
           ALLOCATE(KEISS(MTLOC,KBM1));        KEISS  = 0.0
           ALLOCATE(KEVSS(MTLOC,KBM1));        KEVSS  = 0.0
           ALLOCATE(KEDOC(MTLOC,KBM1));        KEDOC  = 0.0

           ALLOCATE(PM1(MTLOC,KBM1));          PM1   = 0.0   
           ALLOCATE(PM2(MTLOC,KBM1));          PM2   = 0.0
           ALLOCATE(PM3(MTLOC,KBM1));          PM3   = 0.0
           
           ALLOCATE(BMR1(MTLOC,KBM1));         BMR1  = 0.0
           ALLOCATE(BMR2(MTLOC,KBM1));         BMR2  = 0.0
           ALLOCATE(BMR3(MTLOC,KBM1));         BMR3  = 0.0
           
           ALLOCATE(BPR1(MTLOC,KBM1));         BPR1  = 0.0
           ALLOCATE(BPR2(MTLOC,KBM1));         BPR2  = 0.0
           ALLOCATE(BPR3(MTLOC,KBM1));         BPR3  = 0.0
           
           ALLOCATE(DEPTHR(MTLOC,KBM1));       DEPTHR= 0.0!Never used
           ALLOCATE(CCHL1(MTLOC,KBM1));        CCHL1 = 0.0           
           ALLOCATE(CCHL2(MTLOC,KBM1));        CCHL2 = 0.0
           ALLOCATE(CCHL3(MTLOC,KBM1));        CCHL3 = 0.0
           ALLOCATE(Q1(MTLOC,KBM1));           Q1    = 0.0
           ALLOCATE(Q2(MTLOC,KBM1));           Q2    = 0.0
           ALLOCATE(Q3(MTLOC,KBM1));           Q3    = 0.0

!           ALLOCATE(DTUREA(MTLOC,KBM1));       DTUREA = 0.0 !Never used


           ALLOCATE(CFIX(MTLOC));              CFIX  = 0.0
           ALLOCATE(SNFIX(MTLOC));             SNFIX = 0.0

           ALLOCATE(WSS(0:MTLOC,KBM1));        WSS   = 0.0
           ALLOCATE(WSL(0:MTLOC,KBM1));        WSL   = 0.0
           ALLOCATE(WSR(0:MTLOC,KBM1));        WSR   = 0.0
           ALLOCATE(WS1(0:MTLOC,KBM1));        WS1   = 0.0
           ALLOCATE(WS2(0:MTLOC,KBM1));        WS2   = 0.0
           ALLOCATE(WS3(0:MTLOC,KBM1));        WS3   = 0.0
           ALLOCATE(WSU(0:MTLOC,KBM1));        WSU   = 0.0
           ALLOCATE(WSSHI(0:MTLOC,KBM1));      WSSHI = 0.0

           ALLOCATE(FLXSPIB1(MTLOC,KBM1));     FLXSPIB1 = 0.0
           ALLOCATE(FLXSPIB2(MTLOC,KBM1));     FLXSPIB2 = 0.0
           ALLOCATE(FLXSPIB3(MTLOC,KBM1));     FLXSPIB3 = 0.0

           ALLOCATE(V1(0:MTLOC,KBM1));         V1 = 0.0    
           ALLOCATE(V2(0:MTLOC,KBM1));         V2 = 0.0
           
           ALLOCATE(WSSNET(MTLOC));            WSSNET = 0.0
           ALLOCATE(WSLNET(MTLOC));            WSLNET = 0.0
           ALLOCATE(WSRNET(MTLOC));            WSRNET = 0.0
           ALLOCATE(WS1NET(MTLOC));            WS1NET = 0.0
           ALLOCATE(WS2NET(MTLOC));            WS2NET = 0.0
           ALLOCATE(WS3NET(MTLOC));            WS3NET = 0.0
           
           
           !ALLOCATE(DIAGP(MTLOC));             DIAGP  = 0.0
           !ALLOCATE(DIAGN(MTLOC));             DIAGN  = 0.0
           !ALLOCATE(DIAGC(MTLOC));             DIAGC  = 0.0
           !ALLOCATE(DIAGS(MTLOC));             DIAGS  = 0.0
           !ALLOCATE(MTVEL(MTLOC));             MTVEL  = 0.0    !WL moved to mod_sed.F
           ALLOCATE(WSUNET(MTLOC));            WSUNET = 0.0
           
!           ALLOCATE(BL(0:MTLOC,KBM1,3));       BL   = 0.0
!           ALLOCATE(V1S(0:MTLOC,KBM1));        V1S  = 0.0
!           ALLOCATE(HMV(0:MTLOC,KBM1));        HMV  = 0.0
!           ALLOCATE(HMBV(MTLOC));              HMBV = 0.0
!           ALLOCATE(ZD(0:MTLOC,KBM1));         ZD   = 0.0

           ALLOCATE(DTM(0:MTLOC,KBM1,NCP));    DTM = 0.0

           
        
           ALLOCATE(ALGDOC(MTLOC,KBM1));       ALGDOC  = 0.0 
           ALLOCATE(ALGDOC_OUT(MTLOC,KBM1));       ALGDOC_OUT  = 0.0 		   
           ALLOCATE(ALGDOC_GL(MGL,KBM1));       ALGDOC_GL  = 0.0		   
           ALLOCATE(ALGPOC(MTLOC,KBM1));       ALGPOC  = 0.0
           ALLOCATE(ALGPOC_OUT(MTLOC,KBM1));       ALGPOC_OUT  = 0.0		   
           ALLOCATE(ALGPOC_GL(MGL,KBM1));       ALGPOC_GL  = 0.0
		   
		   ALLOCATE(ALGPPDOC(MTLOC,KBM1));       ALGDOC  = 0.0 
           ALLOCATE(ALGPPDOC_GL(MGL,KBM1));       ALGDOC_GL  = 0.0	
		   ALLOCATE(ALGPredDOC(MTLOC,KBM1));       ALGDOC  = 0.0 
           ALLOCATE(ALGPredDOC_GL(MGL,KBM1));       ALGDOC_GL  = 0.0	
           
		   ALLOCATE(HDRLPOC(MTLOC,KBM1));      HDRLPOC = 0.0
		   ALLOCATE(HDRLPOC_OUT(MTLOC,KBM1));      HDRLPOC_OUT = 0.0
           ALLOCATE(HDRLPOC_GL(MGL,KBM1));      HDRLPOC_GL = 0.0	
		   
           ALLOCATE(HDRRPOC(MTLOC,KBM1));      HDRRPOC = 0.0
           ALLOCATE(HDRRPOC_OUT(MTLOC,KBM1));      HDRRPOC_OUT = 0.0		   
           ALLOCATE(HDRRPOC_GL(MGL,KBM1));      HDRRPOC_GL = 0.0
		   
           !ALLOCATE(KLDC(MTLOC,KBM1));         KLDC = 0.0   
           !ALLOCATE(KRDC(MTLOC,KBM1));         KRDC = 0.0
           ALLOCATE(KDCALG(MTLOC,KBM1));       KDCALG = 0.0
           ALLOCATE(KLPC(MTLOC,KBM1));         KLPC = 0.0
           ALLOCATE(KLCALG(MTLOC,KBM1));       KLCALG = 0.0
           ALLOCATE(KRPC(MTLOC,KBM1));         KRPC = 0.0
           ALLOCATE(KSUA(MTLOC,KBM1));         KSUA = 0.0
           ALLOCATE(KCOD(MTLOC,KBM1));         KCOD = 0.0
           !ALLOCATE(KLDN(MTLOC,KBM1));         KLDN = 0.0
           !ALLOCATE(KRDN(MTLOC,KBM1));         KRDN = 0.0
           ALLOCATE(KDNALG(MTLOC,KBM1));       KDNALG = 0.0
           !ALLOCATE(KLDP(MTLOC,KBM1));         KLDP = 0.0
           !ALLOCATE(KRDP(MTLOC,KBM1));         KRDP = 0.0
           ALLOCATE(KDPALG(MTLOC,KBM1));       KDPALG = 0.0
           ALLOCATE(KLPN(MTLOC,KBM1));         KLPN = 0.0
           ALLOCATE(KLNALG(MTLOC,KBM1));       KLNALG = 0.0
           ALLOCATE(KLPP(MTLOC,KBM1));         KLPP = 0.0
           ALLOCATE(KLPALG(MTLOC,KBM1));       KLPALG = 0.0
           ALLOCATE(KRPN(MTLOC,KBM1));         KRPN = 0.0
           ALLOCATE(KRPP(MTLOC,KBM1));         KRPP = 0.0
           ALLOCATE(NTM(MTLOC,KBM1));          NTM = 0.0
           ALLOCATE(KRCOAG(MTLOC,KBM1));       KRCOAG = 0.0

           ALLOCATE(FTMNL(MTLOC,KBM1));        FTMNL = 0.0  
           ALLOCATE(FTHDR(MTLOC,KBM1));        FTHDR = 0.0
           ALLOCATE(RESPC(MTLOC,KBM1));        RESPC = 0.0   
           ALLOCATE(DLALGC(MTLOC,KBM1));       DLALGC = 0.0
           
           ALLOCATE(ALGNH4(MTLOC,KBM1));       ALGNH4 = 0.0
           ALLOCATE(ALGNH4_OUT(MTLOC,KBM1));       ALGNH4_OUT = 0.0		   
           ALLOCATE(ALGNO3(MTLOC,KBM1));       ALGNO3 = 0.0
          ALLOCATE(ALGNO3_OUT(MTLOC,KBM1));       ALGNO3_OUT = 0.0		   
           ALLOCATE(ALGNH4_GL(MGL,KBM1));       ALGNH4_GL = 0.0
           ALLOCATE(ALGNO3_GL(MGL,KBM1));       ALGNO3_GL = 0.0
		   
		   
           ALLOCATE(ALGDON(MTLOC,KBM1));       ALGDON = 0.0
           ALLOCATE(ALGDON_OUT(MTLOC,KBM1));       ALGDON_OUT = 0.0		   
		   ALLOCATE(ALGDON_GL(MGL,KBM1));       ALGDON_GL = 0.0
		   
		   ALLOCATE(ALGPPDON(MTLOC,KBM1));       ALGDOC  = 0.0 
           ALLOCATE(ALGPPDON_GL(MGL,KBM1));       ALGDOC_GL  = 0.0	
		   ALLOCATE(ALGPredDON(MTLOC,KBM1));       ALGDOC  = 0.0 
           ALLOCATE(ALGPredDON_GL(MGL,KBM1));       ALGDOC_GL  = 0.0	
		   
           ALLOCATE(ALGPON(MTLOC,KBM1));       ALGPON = 0.0
           ALLOCATE(ALGPON_OUT(MTLOC,KBM1));       ALGPON_OUT = 0.0		   
           ALLOCATE(ALGPON_GL(MGL,KBM1));       ALGPON_GL = 0.0		   
		   
           ALLOCATE(DENNO3(MTLOC,KBM1));       DENNO3 = 0.0
           ALLOCATE(DENNO3_OUT(MTLOC,KBM1));       DENNO3_OUT = 0.0		   
           ALLOCATE(DENNO3_GL(MGL,KBM1));       DENNO3_GL = 0.0		   
		   
           !ALLOCATE(MNLLDON(MTLOC,KBM1));      MNLLDON = 0.0
		   
           ALLOCATE(HDRLPON(MTLOC,KBM1));      HDRLPON = 0.0
           ALLOCATE(HDRLPON_OUT(MTLOC,KBM1));      HDRLPON = 0.0		   
           ALLOCATE(HDRLPON_GL(MGL,KBM1));      HDRLPON_GL = 0.0	
		   
           ALLOCATE(HDRRPON(MTLOC,KBM1));      HDRRPON = 0.0
           ALLOCATE(HDRRPON_OUT(MTLOC,KBM1));      HDRRPON_OUT = 0.0		   
           ALLOCATE(HDRRPON_GL(MGL,KBM1));      HDRRPON_GL = 0.0	
		   
           !ALLOCATE(MNLRDON(MTLOC,KBM1));      MNLRDON = 0.0
           ALLOCATE(NFIX(MTLOC,KBM1));         NFIX = 0.0

           ALLOCATE(ALGPO4(MTLOC,KBM1));       ALGPO4 = 0.0
           ALLOCATE(ALGDOP(MTLOC,KBM1));       ALGDOP = 0.0
           ALLOCATE(ALGPOP(MTLOC,KBM1));       ALGPOP = 0.0
           !ALLOCATE(MNLLDOP(MTLOC,KBM1));      MNLLDOP = 0.0
           ALLOCATE(HDRLPOP(MTLOC,KBM1));      HDRLPOP = 0.0
           ALLOCATE(HDRRPOP(MTLOC,KBM1));      HDRRPOP = 0.0
           !ALLOCATE(MNLRDOP(MTLOC,KBM1));      MNLRDOP = 0.0
           
           ALLOCATE(DORALG(MTLOC,KBM1));       DORALG = 0.0
           ALLOCATE(DOPR(MTLOC,KBM1));         DOPR = 0.0
           ALLOCATE(DCOD(MTLOC,KBM1));         DCOD = 0.0
           ALLOCATE(DDOC(MTLOC,KBM1));         DDOC = 0.0
           ALLOCATE(NITRIF(MTLOC,KBM1));       NITRIF = 0.0
           
           ALLOCATE(PSD(MTLOC,KBM1));          PSD = 0.0
           ALLOCATE(SAP(MTLOC,KBM1));          SAP = 0.0
           ALLOCATE(ALGUP(MTLOC,KBM1));        ALGUP = 0.0
           ALLOCATE(ALGRES(MTLOC,KBM1));       ALGRES = 0.0

           ALLOCATE(BENDOC(MTLOC));            BENDOC = 0.0
           ALLOCATE(BENNH4(MTLOC));            BENNH4 = 0.0
           ALLOCATE(BENNO3(MTLOC));            BENNO3 = 0.0
           ALLOCATE(BENDON(MTLOC));            BENDON = 0.0
           ALLOCATE(BENDOP(MTLOC));            BENDOP = 0.0
           ALLOCATE(BENPO4(MTLOC));            BENPO4 = 0.0
           ALLOCATE(BENCOD(MTLOC));            BENCOD = 0.0
           ALLOCATE(BENDO(MTLOC));             BENDO = 0.0
           ALLOCATE(BENSA(MTLOC));             BENSA = 0.0
           ALLOCATE(BENDEN(MTLOC));            BENDEN = 0.0
           ALLOCATE(BENCH4G(MTLOC));           BENCH4G = 0.0
           ALLOCATE(BENCH4A(MTLOC));           BENCH4A = 0.0

           
           ALLOCATE(C1(0:MTLOC,KBM1,NCP));         C1    = 0.0
           ALLOCATE(C2(0:MTLOC,KBM1,NCP));         C2    = 0.0
		!    ALLOCATE(C3(0:MTLOC,KBM1,NAC));         C3    = 0.0 ! B Clark added so that the work array in advection model is smaller than c2, should speed up a bit Sep 29 2015
           ALLOCATE(C2F(0:MTLOC,KBM1,NCP));        C2F    = 0.0
           ALLOCATE(CSTAR(0:MTLOC,KBM1,NCP));      CSTAR = 0.0
           ALLOCATE(DTC(0:MTLOC,KBM1,NCP));        DTC   = 0.0 
           ALLOCATE(AC1(0:MTLOC,KBM1,NCP));        AC1   = 0.0
		   
           ALLOCATE(XFLUX(0:MTLOC,KBM1,NCP))
           ALLOCATE(XFLUX_ADV(0:MTLOC,KBM1,NCP))
           
           ALLOCATE(QLIT(MTLOC,KBM1));             QLIT = 0.0
           ALLOCATE(CONLIT(MTLOC,KBM1,NCP));       CONLIT = 0.0

           ALLOCATE(AFI1(MTLOC,KBM1));       AFI1 = 0.0 
           ALLOCATE(ANL1(MTLOC,KBM1));       ANL1 = 0.0
           ALLOCATE(APL1(MTLOC,KBM1));       APL1 = 0.0
           ALLOCATE(AFI2(MTLOC,KBM1));       AFI2 = 0.0
           ALLOCATE(ANL2(MTLOC,KBM1));       ANL2 = 0.0
           ALLOCATE(APL2(MTLOC,KBM1));       APL2 = 0.0
           ALLOCATE(ASL2(MTLOC,KBM1));       ASL2 = 0.0
           ALLOCATE(AFI3(MTLOC,KBM1));       AFI3 = 0.0
           ALLOCATE(ANL3(MTLOC,KBM1));       ANL3 = 0.0
           ALLOCATE(APL3(MTLOC,KBM1));       APL3 = 0.0
           ALLOCATE(ANPP(MTLOC,KBM1));       ANPP = 0.0
           ALLOCATE(ARESP(MTLOC,KBM1));      ARESP = 0.0
           ALLOCATE(AKE(MTLOC,KBM1));        AKE = 0.0
           ALLOCATE(ASL1(MTLOC,KBM1));       ASL1 = 0.0
           ALLOCATE(ASL3(MTLOC,KBM1));       ASL3 = 0.0
           ALLOCATE(AGPP(MTLOC,KBM1));       AGPP = 0.0
           
           ALLOCATE(ACCHL1(MTLOC,KBM1));     ACCHL1 = 0.0  
           ALLOCATE(ACCHL2(MTLOC,KBM1));     ACCHL2 = 0.0
           ALLOCATE(ACCHL3(MTLOC,KBM1));     ACCHL3 = 0.0

           
           ALLOCATE(ABENDOC(MTLOC));          ABENDOC = 0.0  
           ALLOCATE(ABENNH4(MTLOC));          ABENNH4 = 0.0
           ALLOCATE(ABENNO3(MTLOC));          ABENNO3 = 0.0
           ALLOCATE(ABENPO4(MTLOC));          ABENPO4 = 0.0
           ALLOCATE(ABENCOD(MTLOC));          ABENCOD = 0.0
           ALLOCATE(ABENDO(MTLOC));           ABENDO = 0.0
           ALLOCATE(ABENSA(MTLOC));           ABENSA = 0.0
           ALLOCATE(ASSFWS(MTLOC));           ASSFWS = 0.0
           ALLOCATE(APCFWS(MTLOC));           APCFWS = 0.0
           ALLOCATE(APNFWS(MTLOC));           APNFWS = 0.0
           ALLOCATE(APPFWS(MTLOC));           APPFWS = 0.0
           ALLOCATE(APSFWS(MTLOC));           APSFWS = 0.0
           ALLOCATE(ACPIP(MTLOC));            ACPIP = 0.0
           ALLOCATE(ACPOS(MTLOC));            ACPOS = 0.0
           ALLOCATE(ADFEED(MTLOC));           ADFEED = 0.0
           ALLOCATE(ABENCH4G(MTLOC));         ABENCH4G = 0.0
           ALLOCATE(ABENCH4A(MTLOC));         ABENCH4A = 0.0

           ALLOCATE(AFIB(MTLOC));             AFIB = 0.0
           ALLOCATE(ANLB(MTLOC));             ANLB = 0.0
           ALLOCATE(APLB(MTLOC));             APLB = 0.0
           ALLOCATE(ANPPB(MTLOC));            ANPPB = 0.0
           ALLOCATE(ABBM(MTLOC));             ABBM = 0.0
           ALLOCATE(ABLITE(MTLOC));           ABLITE = 0.0
           ALLOCATE(ABMB(MTLOC));             ABMB = 0.0
           ALLOCATE(APB(MTLOC));              APB = 0.0
           ALLOCATE(APRB(MTLOC));             APRB = 0.0
           ALLOCATE(ABADOC(MTLOC));           ABADOC = 0.0
           ALLOCATE(ABAPOC(MTLOC));           ABAPOC = 0.0
           ALLOCATE(ABANH4(MTLOC));           ABANH4 = 0.0
           ALLOCATE(ABANO3(MTLOC));           ABANO3 = 0.0
           ALLOCATE(ABAPON(MTLOC));           ABAPON = 0.0
           ALLOCATE(ABAPO4(MTLOC));           ABAPO4 = 0.0
           ALLOCATE(ABAPOP(MTLOC));           ABAPOP = 0.0
           ALLOCATE(ABADO(MTLOC));            ABADO = 0.0
           
           ALLOCATE(ACPOC(MTLOC,3));          ACPOC = 0.0
           ALLOCATE(ACPON(MTLOC,3));          ACPON = 0.0
           ALLOCATE(ACPOP(MTLOC,3));          ACPOP = 0.0
           ALLOCATE(ASFEED(MTLOC,NSSFP));     ASFEED = 0.0
              
           ALLOCATE(AJNSF(MTLOC));            AJNSF = 0.0    
           ALLOCATE(AJPSF(MTLOC));            AJPSF = 0.0
           ALLOCATE(ASODSF(MTLOC));           ASODSF = 0.0
           ALLOCATE(ASASF(MTLOC));            ASASF = 0.0
           ALLOCATE(ASUSF(MTLOC));            ASUSF = 0.0
           ALLOCATE(ASFGCIN(MTLOC));          ASFGCIN = 0.0
           ALLOCATE(ASFCFEC(MTLOC));          ASFCFEC = 0.0
           ALLOCATE(ASFCPSF(MTLOC));          ASFCPSF = 0.0
           ALLOCATE(AFLXCSF(MTLOC));          AFLXCSF = 0.0
           ALLOCATE(AFLXNSF(MTLOC));          AFLXNSF = 0.0
           ALLOCATE(AFLXPSF(MTLOC));          AFLXPSF = 0.0
           ALLOCATE(ARPOCSF(MTLOC));          ARPOCSF = 0.0
           ALLOCATE(ARPONSF(MTLOC));          ARPONSF = 0.0
           ALLOCATE(ARPOPSF(MTLOC));          ARPOPSF = 0.0
           ALLOCATE(ASSISF(MTLOC));           ASSISF = 0.0
           ALLOCATE(ASSISASF(MTLOC));         ASSISASF = 0.0
           ALLOCATE(ASSISUSF(MTLOC));         ASSISUSF = 0.0
           ALLOCATE(ASSIPSF(MTLOC));          ASSIPSF = 0.0
           
           ALLOCATE(A_T(MTLOC,KBM1));              A_T = 0.0
           ALLOCATE(AP1(MTLOC,KBM1));              AP1 = 0.0
           ALLOCATE(ABM1(MTLOC,KBM1));             ABM1 = 0.0
           ALLOCATE(APR1(MTLOC,KBM1));             APR1 = 0.0
           ALLOCATE(AP2(MTLOC,KBM1));              AP2 = 0.0
           ALLOCATE(ABM2(MTLOC,KBM1));             ABM2 = 0.0
           ALLOCATE(APR2(MTLOC,KBM1));             APR2 = 0.0
           ALLOCATE(AP3(MTLOC,KBM1));              AP3 = 0.0
           ALLOCATE(ABM3(MTLOC,KBM1));             ABM3  = 0.0
           ALLOCATE(APR3(MTLOC,KBM1));             APR3 = 0.0
           ALLOCATE(AALGDOC(MTLOC,KBM1));          AALGDOC = 0.0
           ALLOCATE(AALGPOC(MTLOC,KBM1));          AALGPOC = 0.0
           ALLOCATE(ADENIT(MTLOC,KBM1));           ADENIT = 0.0
           ALLOCATE(AMNLDOC(MTLOC,KBM1));          AMNLDOC = 0.0
           ALLOCATE(AHDRPOC(MTLOC,KBM1));          AHDRPOC = 0.0
           ALLOCATE(AALGNH4(MTLOC,KBM1));          AALGNH4 = 0.0
           ALLOCATE(AALGNO3(MTLOC,KBM1));          AALGNO3 = 0.0
           ALLOCATE(AALGDON(MTLOC,KBM1));          AALGDON = 0.0
           ALLOCATE(AALGPON(MTLOC,KBM1));          AALGPON = 0.0
           ALLOCATE(ANT(MTLOC,KBM1));              ANT = 0.0
           ALLOCATE(ADENNO3(MTLOC,KBM1));          ADENNO3  = 0.0
           ALLOCATE(AMNLDON(MTLOC,KBM1));          AMNLDON = 0.0
           ALLOCATE(AHDRPON(MTLOC,KBM1));          AHDRPON = 0.0
           

           ALLOCATE(AALGPO4(MTLOC,KBM1));          AALGPO4 = 0.0
           ALLOCATE(AALGDOP(MTLOC,KBM1));          AALGDOP = 0.0
           ALLOCATE(AALGPOP(MTLOC,KBM1));          AALGPOP = 0.0
           ALLOCATE(AMNLDOP(MTLOC,KBM1));          AMNLDOP = 0.0
           ALLOCATE(AHDRPOP(MTLOC,KBM1));          AHDRPOP = 0.0
           ALLOCATE(APSD(MTLOC,KBM1));             APSD = 0.0
           ALLOCATE(ASAP(MTLOC,KBM1));             ASAP = 0.0
           ALLOCATE(AALGUP(MTLOC,KBM1));           AALGUP = 0.0
           ALLOCATE(AALGRES(MTLOC,KBM1));          AALGRES = 0.0
           ALLOCATE(ADO(MTLOC,KBM1));              ADO = 0.0
           ALLOCATE(ADORALG(MTLOC,KBM1));          ADORALG = 0.0
           ALLOCATE(ADOPR(MTLOC,KBM1));            ADOPR = 0.0
           ALLOCATE(ADCOD(MTLOC,KBM1));            ADCOD = 0.0
           ALLOCATE(ADDOC(MTLOC,KBM1));            ADDOC = 0.0
           ALLOCATE(ANITRIF(MTLOC,KBM1));          ANITRIF = 0.0

           ALLOCATE(V1SINGLE(0:MTLOC,KBM1)); V1SINGLE = 0.0

           ALLOCATE(SBN(MTLOC));          SBN = 0  
           ALLOCATE(BBN(MTLOC));          BBN = 0
           ALLOCATE(HMSBV(MTLOC));        HMSBV = 0 
            !Wen Long   ALLOCATE(SFA(MTLOC));          SFA = 0   !Wen Long deprecated SFA, replaced by ART1    

           ALLOCATE(FLUXS(0:MTLOC,KBM1,13));         FLUXS = 0.0

           ALLOCATE(IWCMNB(MTLOC,KBM1));        IWCMNB = 0.0
           ALLOCATE(IWCMPB(MTLOC,KBM1));        IWCMPB = 0.0
           ALLOCATE(IWCMCB(MTLOC,KBM1));        IWCMCB = 0.0
           ALLOCATE(IWCMSB(MTLOC,KBM1));        IWCMSB = 0.0
           ALLOCATE(WCMNB(MTLOC,KBM1));         WCMNB = 0.0
           ALLOCATE(WCMPB(MTLOC,KBM1));         WCMPB = 0.0
           ALLOCATE(WCMCB(MTLOC,KBM1));         WCMCB = 0.0
           ALLOCATE(WCMSB(MTLOC,KBM1));         WCMSB = 0.0
           ALLOCATE(DLWCMNB(MTLOC,KBM1));       DLWCMNB = 0.0
           ALLOCATE(DLWCMPB(MTLOC,KBM1));       DLWCMPB = 0.0
           ALLOCATE(DLWCKMNB(MTLOC,KBM1));      DLWCKMNB = 0.0
           ALLOCATE(DLWCKMCB(MTLOC,KBM1));      DLWCKMCB = 0.0
           ALLOCATE(S1FLXNB(MTLOC,KBM1));       S1FLXNB = 0.0
           ALLOCATE(S1FLXPB(MTLOC,KBM1));       S1FLXPB = 0.0
           ALLOCATE(S1FLXCB(MTLOC,KBM1));       S1FLXCB = 0.0
           ALLOCATE(S2FLXNB(MTLOC,KBM1));       S2FLXNB = 0.0
           ALLOCATE(S2FLXPB(MTLOC,KBM1));       S2FLXPB = 0.0
           ALLOCATE(S2FLXCB(MTLOC,KBM1));       S2FLXCB = 0.0
           ALLOCATE(S3FLXNB(MTLOC,KBM1));       S3FLXNB = 0.0
           ALLOCATE(S3FLXPB(MTLOC,KBM1));       S3FLXPB = 0.0
           ALLOCATE(S3FLXCB(MTLOC,KBM1));       S3FLXCB = 0.0
           

!           ALLOCATE(CPOC(MTLOC,3));             CPOC   = 0.0
!           ALLOCATE(CPOP(MTLOC,3));             CPOP   = 0.0
!           ALLOCATE(CPON(MTLOC,3));             CPON   = 0.0
!           ALLOCATE(CPOS(MTLOC));               CPOS   = 0.0      
!           
!           ALLOCATE(JPOC(MTLOC,3));           JPOC = 0.0
!           ALLOCATE(JPOP(MTLOC,3));           JPOP = 0.0 
!           ALLOCATE(JPON(MTLOC,3));           JPON = 0.0
!           ALLOCATE(JPOS(MTLOC));             JPOS = 0.0     
!           
!           ALLOCATE(CPO4(MTLOC));               CPO4   = 0.0      
!           ALLOCATE(CNO3(MTLOC));               CNO3   = 0.0
!           ALLOCATE(CNH4(MTLOC));               CNH4   = 0.0     

!!          ALLOCATE(HSED(MTLOC));               HSED   = 0.0 !moved to mod_sed.F
!           ALLOCATE(CCH4(MTLOC));               CCH4   = 0.0      
!           ALLOCATE(CSO4(MTLOC));               CSO4   = 0.0     
!           ALLOCATE(CHS(MTLOC));                CHS    = 0.0
!           ALLOCATE(CSI(MTLOC));                CSI    = 0.0

           ALLOCATE(CTEMP(MTLOC));              CTEMP  = 0.0
           
           ALLOCATE(BSVOL(MTLOC));              BSVOL  = 0.0     

           ALLOCATE(PCFWS(MTLOC));              PCFWS = 0.0            
           ALLOCATE(PNFWS(MTLOC));              PNFWS = 0.0            
           ALLOCATE(PPFWS(MTLOC));              PPFWS = 0.0
           ALLOCATE(PSFWS(MTLOC));              PSFWS = 0.0 
           ALLOCATE(SSFWS(MTLOC));              SSFWS = 0.0
                
           ALLOCATE(AASRAT(MTLOC));             AASRAT = 0.0
           
           ALLOCATE(ACFIX(MTLOC));              ACFIX  = 0.0
           ALLOCATE(ANFIX(MTLOC,KBM1));            ANFIX = 0.0

           !----
             
           ALLOCATE(BENDOCB(MTLOC));            BENDOCB = 0.0
           ALLOCATE(BENNH4B(MTLOC));            BENNH4B = 0.0
           ALLOCATE(BENNO3B(MTLOC));            BENNO3B = 0.0
           ALLOCATE(BENPO4B(MTLOC));            BENPO4B = 0.0
           ALLOCATE(BENCODB(MTLOC));            BENCODB = 0.0
           ALLOCATE(BENDOB(MTLOC));             BENDOB  = 0.0
           ALLOCATE(BENSAB(MTLOC));             BENSAB  = 0.0
     
           ALLOCATE(ATMFLXNB(MTLOC));           ATMFLXNB    = 0.0  
           ALLOCATE(ATMFLXPB(MTLOC));           ATMFLXPB    = 0.0 
           ALLOCATE(ATMFLXCB(MTLOC));           ATMFLXCB    = 0.0
           ALLOCATE(BENFLXPNB(MTLOC));          BENFLXPNB   = 0.0  
           ALLOCATE(BENFLXDNB(MTLOC));          BENFLXDNB   = 0.0
           ALLOCATE(BENFLXPPB(MTLOC));          BENFLXPPB   = 0.0
           ALLOCATE(BENFLXDPB(MTLOC));          BENFLXDPB   = 0.0
           ALLOCATE(BENFLXPCB(MTLOC));          BENFLXPCB   = 0.0
           ALLOCATE(DLSEDKNB(MTLOC));           DLSEDKNB    = 0.0
           ALLOCATE(DLSEDKCB(MTLOC));           DLSEDKCB    = 0.0 
           ALLOCATE(BURIALFLXNB(MTLOC));        BURIALFLXNB = 0.0
           ALLOCATE(BURIALFLXPB(MTLOC));        BURIALFLXPB = 0.0
           ALLOCATE(BURIALFLXCB(MTLOC));        BURIALFLXCB = 0.0


        ! KURT GLAESEMANN add C2_GL
           ALLOCATE(C2_GL(MGL,KBM1,NCP));   C2_GL = 0.0
        ! RGl added below for output of CCHL1_GL
           ALLOCATE(CCHL1_GL(MGL,KBM1));   CCHL1_GL = 0.0
           ALLOCATE(T_GL(MGL,KBM1));        T_GL = 0.0
           ALLOCATE(S_GL(MGL,KBM1));        S_GL = 0.0
           ALLOCATE(total_netPP_GL(MGL));   total_netPP_GL = 0.0

        !Wen Long added D_GL, EL_GL, H_GL (total depth, surface elevation and bathymetric depth)
           ALLOCATE(D_GL(MGL));        D_GL=0.0;
           ALLOCATE(H_GL(MGL));        H_GL=0.0;
           ALLOCATE(EL_GL(MGL));       EL_GL=0.0;

        !--Wen Long debugging benthic fluxes ---
           ALLOCATE(BFLUX_GL(MGL,9))  ;    BFLUX_GL   =0.0;
           ALLOCATE(BFLUXB_GL(MGL,9)) ;    BFLUXB_GL  =0.0; 
           ALLOCATE(BFLUXNX(MTLOC,9)) ;    BFLUXNX    =0.0;
           ALLOCATE(BFLUXNX_GL(MGL,9));    BFLUXNX_GL =0.0;
           ALLOCATE(BFLUX(MTLOC,9))   ;    BFLUX      =0.0;

        !--Wen Long allocating global hydro variables here ---
           
           ALLOCATE(UL_GL(NGL,KBM1))          ;  UL_GL=0.0;  !LB: had added "0:" because parallel mode gave error "Subscript #1 of the array UL_GL has value 0 which is less than the lower bound of 1"
           ALLOCATE(VL_GL(NGL,KBM1))          ;  VL_GL=0.0;  !but then, PSM did not run.  In revision 117, we need it from 1 for PSM to run!!!
           ALLOCATE(WW_GL(NGL,KBM1))          ;  WW_GL=0.0;
		   ALLOCATE(WW_OUT(NTLOC,KBM1))			  ; WW_OUT=0
		   ALLOCATE(WTSL_GL(MGL,KB))          ;  WTSL_GL=0.0;
           ALLOCATE(KHL_GL(MGL,KB))           ;  KHL_GL=0.0;
           ALLOCATE(SL_GL(MGL,KBM1))          ;  SL_GL=0.0;
           ALLOCATE(TL_GL(MGL,KBM1))          ;  TL_GL=0.0;
           ALLOCATE(ELL_GL(1:MGL))            ;  ELL_GL=0.0;
           ALLOCATE(DTFAL_GL(MGL))            ;  DTFAL_GL=0.0;
           ALLOCATE(UARD_OBCN_GL(NOBTY+1))    ;  UARD_OBCN_GL=0.0;
           ALLOCATE(XFLUX_OBC_GL(NOBTY,KBM1)) ;  XFLUX_OBC_GL=0.0;
         

        !--Wen Long  added the following algae debugging global variables here ---

           ALLOCATE(P1_GL(MGL,KBM1))     ;    P1_GL   = 0.0
           ALLOCATE(P2_GL(MGL,KBM1))     ;    P2_GL   = 0.0
           ALLOCATE(BM1_GL(MGL,KBM1))    ;    BM1_GL  = 0.0
           ALLOCATE(BM2_GL(MGL,KBM1))    ;    BM2_GL  = 0.0
           ALLOCATE(NL1_GL(MGL,KBM1))    ;    NL1_GL  = 0.0
           ALLOCATE(NL2_GL(MGL,KBM1))    ;    NL2_GL  = 0.0
           ALLOCATE(PL1_GL(MGL,KBM1))    ;    PL1_GL  = 0.0
           ALLOCATE(PL2_GL(MGL,KBM1))    ;    PL2_GL  = 0.0
           ALLOCATE(B1SZ_GL(MGL,KBM1))   ;    B1SZ_GL = 0.0
           ALLOCATE(B2SZ_GL(MGL,KBM1))   ;    B2SZ_GL = 0.0
           ALLOCATE(B1LZ_GL(MGL,KBM1))   ;    B1LZ_GL = 0.0
           ALLOCATE(B2LZ_GL(MGL,KBM1))   ;    B2LZ_GL = 0.0
           ALLOCATE(PR1_GL(MGL,KBM1))    ;    PR1_GL  = 0.0
           ALLOCATE(PR2_GL(MGL,KBM1))    ;    PR2_GL  = 0.0
           
           IF(LIGHT_EXTINCTION)THEN
                ALLOCATE(IAVG_GL(MGL,KBM1))   ;    IAVG_GL = 0.0
                ALLOCATE(FI1_GL(MGL,KBM1))    ;    FI1_GL  = 0.0
				ALLOCATE(KDLAMB_int_GL(MGL,KBM1))   ;    KDLAMB_int_GL = 0.0
                ALLOCATE(FI2_GL(MGL,KBM1))    ;    FI2_GL  = 0.0
           ENDIF
		   
		   ALLOCATE(SSI_GL(MGL,KBM1))   ;    SSI_GL = 0.0
		   
           ALLOCATE(IK1_GL(MGL,KBM1))    ;     IK1_GL = 0.0
           ALLOCATE(IK2_GL(MGL,KBM1))    ;     IK2_GL = 0.0
           ALLOCATE(IK1(MTLOC,KBM1))     ;       IK1  = 0.0
           ALLOCATE(IK2(MTLOC,KBM1))     ;       IK2  = 0.0

        !---Wen Long added the following for debugging DOXG

            ALLOCATE(DDOC_GL(MGL,KBM1))   ;  DDOC_GL=0.0
            ALLOCATE(DCOD_GL(MGL,KBM1))   ;  DCOD_GL=0.0
            ALLOCATE(NITRIF_GL(MGL,KBM1)) ;  NITRIF_GL=0.0
            ALLOCATE(DOSZ_GL(MGL,KBM1))   ;  DOSZ_GL=0.0
            ALLOCATE(DOLZ_GL(MGL,KBM1))   ;  DOLZ_GL=0.0
            ALLOCATE(DOPR_GL(MGL,KBM1))   ;  DOPR_GL=0.0
            ALLOCATE(DORALG_GL(MGL,KBM1)) ;  DORALG_GL=0.0
            ALLOCATE(RESP_GL(MGL,KBM1))   ;  RESP_GL=0.0
            ALLOCATE(FTCOD_GL(MGL,KBM1))  ;  FTCOD_GL=0.0
            ALLOCATE(COD_GL(MGL,KBM1))    ;  COD_GL=0.0 
            ALLOCATE(REAERDO_GL(MGL,KBM1));  REAERDO_GL=0.0
            
            ALLOCATE(REAERDO(MTLOC,KBM1)) ;  REAERDO   =0.0   ! local variabl

        !---Wen Long declared a large buffer (so that we do not have to dynamically allocate and deallocate buffers 
        !   all the time--
#if defined (MULTIPROCESSOR)
            IF(.NOT.ALLOCATED(RCVBUF))ALLOCATE(RCVBUF(MAX(MGL,NGL)*3*KB));RCVBUF=0.0;
            IF(.NOT.ALLOCATED(SNDBUF))ALLOCATE(SNDBUF(MAX(MGL,NGL)*3*KB));SNDBUF=0.0;
														  									  
														  
#endif
														  
	! B Clark for SSI Resuspension
			    ALLOCATE(ShearStress(MTLOC))
				ALLOCATE(ShearStress_GL(MGL))
				ALLOCATE(RHO1(MTLOC,KBM1))
				ALLOCATE(RHO1_GL(MGL,KBM1))

														  

           RETURN

        END SUBROUTINE WQM_ALLOC

        SUBROUTINE WQM_DEALLOC

            IMPLICIT NONE
            
           IF(ALLOCATED(T))DEALLOCATE(T)
           IF(ALLOCATED(SALT))DEALLOCATE(SALT)           
           IF(ALLOCATED(SSI))DEALLOCATE(SSI)
           IF(ALLOCATED(B1))DEALLOCATE(B1)
           IF(ALLOCATED(B2))DEALLOCATE(B2)
           IF(ALLOCATED(B3))DEALLOCATE(B3)
           !IF(ALLOCATED(LDOC))DEALLOCATE(LDOC)
           IF(ALLOCATED(LPOC))DEALLOCATE(LPOC)
           IF(ALLOCATED(RPOC))DEALLOCATE(RPOC)
		   
           IF(ALLOCATED(DTLPOC_out))DEALLOCATE(DTLPOC_out)
           IF(ALLOCATED(DTRPOC_out))DEALLOCATE(DTRPOC_out)
           IF(ALLOCATED(DTLPON_out))DEALLOCATE(DTLPON_out)
           IF(ALLOCATED(DTRPON_out))DEALLOCATE(DTRPON_out)
		   
           IF(ALLOCATED(NH4))DEALLOCATE(NH4)
           IF(ALLOCATED(NO3))DEALLOCATE(NO3)
		   
           IF(ALLOCATED(DTNH4_out))DEALLOCATE(DTNH4_out)
           IF(ALLOCATED(DTNO3_out))DEALLOCATE(DTNO3_out)		   
		   
           !IF(ALLOCATED(LDON))DEALLOCATE(LDON)
           IF(ALLOCATED(LPON))DEALLOCATE(LPON)
           IF(ALLOCATED(RPON))DEALLOCATE(RPON)
           IF(ALLOCATED(PO4))DEALLOCATE(PO4)
           !IF(ALLOCATED(LDOP))DEALLOCATE(LDOP)
           IF(ALLOCATED(LPOP))DEALLOCATE(LPOP)
           IF(ALLOCATED(RPOP))DEALLOCATE(RPOP)
           IF(ALLOCATED(COD))DEALLOCATE(COD)
           IF(ALLOCATED(DOXG))DEALLOCATE(DOXG)
           IF(ALLOCATED(SIUPB))DEALLOCATE(SIUPB)
           IF(ALLOCATED(SIAT))DEALLOCATE(SIAT)
           

           !IF(ALLOCATED(RDOC))DEALLOCATE(RDOC)
           !IF(ALLOCATED(RDON))DEALLOCATE(RDON)
           !IF(ALLOCATED(RDOP))DEALLOCATE(RDOP)
           
           IF(ALLOCATED(UREA))DEALLOCATE(UREA)
           
           IF(ALLOCATED(PIP))DEALLOCATE(PIP)
           IF(ALLOCATED(PIB1))DEALLOCATE(PIB1)
           IF(ALLOCATED(PIB2))DEALLOCATE(PIB2)
           IF(ALLOCATED(PIB3))DEALLOCATE(PIB3)

           IF(ALLOCATED(PN1))DEALLOCATE(PN1)
           IF(ALLOCATED(PN2))DEALLOCATE(PN2)
           IF(ALLOCATED(PN3))DEALLOCATE(PN3)
           
           IF(ALLOCATED(BM1))DEALLOCATE(BM1)
           IF(ALLOCATED(BM2))DEALLOCATE(BM2)
           IF(ALLOCATED(BM3))DEALLOCATE(BM3)
           
    !       IF(ALLOCATED(MNLLDOC))DEALLOCATE(MNLLDOC)
    !       IF(ALLOCATED(MNLRDOC))DEALLOCATE(MNLRDOC)
           
           IF(ALLOCATED(FTCOD))DEALLOCATE(FTCOD)
           
           IF(ALLOCATED(NPP))DEALLOCATE(NPP)
           IF(ALLOCATED(P1))DEALLOCATE(P1)
           IF(ALLOCATED(P2))DEALLOCATE(P2)
           IF(ALLOCATED(P3))DEALLOCATE(P3)
           
           IF(ALLOCATED(NT))DEALLOCATE(NT)
           IF(ALLOCATED(NT_OUT))DEALLOCATE(NT_OUT)		   
           IF(ALLOCATED(NT_GL))DEALLOCATE(NT_GL)	
		   
           IF(ALLOCATED(PR1))DEALLOCATE(PR1)
           IF(ALLOCATED(PR2))DEALLOCATE(PR2)
           IF(ALLOCATED(PR3))DEALLOCATE(PR3)
           
           IF(ALLOCATED(DENIT))DEALLOCATE(DENIT)
           IF(ALLOCATED(DENIT_OUT))DEALLOCATE(DENIT_OUT)		   
           IF(ALLOCATED(DENIT_GL))DEALLOCATE(DENIT_GL)
		   
           IF(ALLOCATED(DENITN))DEALLOCATE(DENITN)
           IF(ALLOCATED(DENITN_OUT))DEALLOCATE(DENITN_OUT)		   
           IF(ALLOCATED(DENITN_GL))DEALLOCATE(DENITN_GL)
		   
           IF(ALLOCATED(DENITP))DEALLOCATE(DENITP)

           IF(ALLOCATED(RATOX))DEALLOCATE(RATOX)
           IF(ALLOCATED(GPP))DEALLOCATE(GPP)

           IF(ALLOCATED(P1NNF))DEALLOCATE(P1NNF)
           
           IF(ALLOCATED(ASRAT))DEALLOCATE(ASRAT)
           IF(ALLOCATED(total_netPP))DEALLOCATE(total_netPP)

           IF(ALLOCATED(KE))DEALLOCATE(KE)
           IF(ALLOCATED(KEISS))DEALLOCATE(KEISS)
           IF(ALLOCATED(KEVSS))DEALLOCATE(KEVSS)
           IF(ALLOCATED(KEDOC))DEALLOCATE(KEDOC)
           
           IF(ALLOCATED(FI1))DEALLOCATE(FI1)
           IF(ALLOCATED(FI2))DEALLOCATE(FI2)
           IF(ALLOCATED(FI3))DEALLOCATE(FI3)
           IF(ALLOCATED(NL1))DEALLOCATE(NL1)
           IF(ALLOCATED(NL2))DEALLOCATE(NL2)
           IF(ALLOCATED(NL3))DEALLOCATE(NL3)
           IF(ALLOCATED(PL1))DEALLOCATE(PL1)
           IF(ALLOCATED(PL2))DEALLOCATE(PL2)
           IF(ALLOCATED(PL3))DEALLOCATE(PL3)

           IF(ALLOCATED(RESP))DEALLOCATE(RESP)

           IF(ALLOCATED(SL1))DEALLOCATE(SL1)
           IF(ALLOCATED(SL2))DEALLOCATE(SL2)           
           IF(ALLOCATED(SL3))DEALLOCATE(SL3)
           

           IF(ALLOCATED(PM1))DEALLOCATE(PM1)
           IF(ALLOCATED(PM2))DEALLOCATE(PM2)
           IF(ALLOCATED(PM3))DEALLOCATE(PM3)
           
           IF(ALLOCATED(BMR1))DEALLOCATE(BMR1)
           IF(ALLOCATED(BMR2))DEALLOCATE(BMR2)
           IF(ALLOCATED(BMR3))DEALLOCATE(BMR3)
           
           IF(ALLOCATED(BPR1))DEALLOCATE(BPR1)
           IF(ALLOCATED(BPR2))DEALLOCATE(BPR2)
           IF(ALLOCATED(BPR3))DEALLOCATE(BPR3)
           
           IF(ALLOCATED(DEPTHR))DEALLOCATE(DEPTHR)!Never used
           
           IF(ALLOCATED(CCHL1))DEALLOCATE(CCHL1)               
           IF(ALLOCATED(CCHL2))DEALLOCATE(CCHL2)
           IF(ALLOCATED(CCHL3))DEALLOCATE(CCHL3)
           
           IF(ALLOCATED(Q1))DEALLOCATE(Q1)
           IF(ALLOCATED(Q2))DEALLOCATE(Q2)
           IF(ALLOCATED(Q3))DEALLOCATE(Q3)

!          IF(ALLOCATED(DTUREA))DEALLOCATE()

           IF(ALLOCATED(CFIX))DEALLOCATE(CFIX)
           IF(ALLOCATED(SNFIX))DEALLOCATE(SNFIX)

           IF(ALLOCATED(WSS))DEALLOCATE(WSS)
           IF(ALLOCATED(WSL))DEALLOCATE(WSL)
           IF(ALLOCATED(WSR))DEALLOCATE(WSR)
           
           IF(ALLOCATED(WS1))DEALLOCATE(WS1)
           IF(ALLOCATED(WS2))DEALLOCATE(WS2)
           IF(ALLOCATED(WS3))DEALLOCATE(WS3)
           IF(ALLOCATED(WSU))DEALLOCATE(WSU)
           IF(ALLOCATED(WSSHI))DEALLOCATE(WSSHI)

           IF(ALLOCATED(FLXSPIB1))DEALLOCATE(FLXSPIB1)
           IF(ALLOCATED(FLXSPIB2))DEALLOCATE(FLXSPIB2)
           IF(ALLOCATED(FLXSPIB3))DEALLOCATE(FLXSPIB3)

           IF(ALLOCATED(V1))DEALLOCATE(V1)
           IF(ALLOCATED(V2))DEALLOCATE(V2)
           
           IF(ALLOCATED(WSSNET))DEALLOCATE(WSSNET)
           IF(ALLOCATED(WSLNET))DEALLOCATE(WSLNET)
           IF(ALLOCATED(WSRNET))DEALLOCATE(WSRNET)
           IF(ALLOCATED(WS1NET))DEALLOCATE(WS1NET)
           IF(ALLOCATED(WS2NET))DEALLOCATE(WS2NET)
           IF(ALLOCATED(WS3NET))DEALLOCATE(WS3NET)
           IF(ALLOCATED(WSUNET))DEALLOCATE(WSUNET)           
           
           !IF(ALLOCATED(DIAGC))DEALLOCATE(DIAGC)           
           !IF(ALLOCATED(DIAGN))DEALLOCATE(DIAGN)           
           !IF(ALLOCATED(DIAGP))DEALLOCATE(DIAGP)
           !IF(ALLOCATED(DIAGS))DEALLOCATE(DIAGS)

           !IF(ALLOCATED(MTVEL))DEALLOCATE(MTVEL)    !WLong moved to mod_sed.F

           
!           IF(ALLOCATED(BL))DEALLOCATE(BL)
!           IF(ALLOCATED(V1S))DEALLOCATE(V1S)
!           IF(ALLOCATED(HMV))DEALLOCATE(HMV)
!           IF(ALLOCATED(HMBV))DEALLOCATE(HMBV)
!           IF(ALLOCATED(ZD))DEALLOCATE(ZD)

           IF(ALLOCATED(DTM))DEALLOCATE(DTM)

        
           IF(ALLOCATED(ALGDOC))DEALLOCATE(ALGDOC)
           IF(ALLOCATED(ALGDOC_OUT))DEALLOCATE(ALGDOC_OUT)		   
           IF(ALLOCATED(ALGDOC_GL))DEALLOCATE(ALGDOC_GL)		   
           IF(ALLOCATED(ALGPOC))DEALLOCATE(ALGPOC)
           IF(ALLOCATED(ALGPOC_OUT))DEALLOCATE(ALGPOC_OUT)		   
           IF(ALLOCATED(ALGPOC_GL))DEALLOCATE(ALGPOC_GL)
		   
		   IF(ALLOCATED(ALGPredDOC))DEALLOCATE(ALGPredDOC)
           IF(ALLOCATED(ALGPredDOC_GL))DEALLOCATE(ALGPredDOC_GL)		   
           IF(ALLOCATED(ALGPPDOC))DEALLOCATE(ALGPPDOC)
           IF(ALLOCATED(ALGPPDOC_GL))DEALLOCATE(ALGPPDOC_GL)	
           
           IF(ALLOCATED(HDRLPOC))DEALLOCATE(HDRLPOC)
           IF(ALLOCATED(HDRRPOC))DEALLOCATE(HDRRPOC)
           IF(ALLOCATED(HDRLPOC_OUT))DEALLOCATE(HDRLPOC_OUT)
           IF(ALLOCATED(HDRRPOC_OUT))DEALLOCATE(HDRRPOC_OUT)		   
           IF(ALLOCATED(HDRLPOC_GL))DEALLOCATE(HDRLPOC_GL)
           IF(ALLOCATED(HDRRPOC_GL))DEALLOCATE(HDRRPOC_GL)
           
           !IF(ALLOCATED(KLDC))DEALLOCATE(KLDC)
           !IF(ALLOCATED(KRDC))DEALLOCATE(KRDC)

           IF(ALLOCATED(KDCALG))DEALLOCATE(KDCALG)
           IF(ALLOCATED(KLCALG))DEALLOCATE(KLCALG)
           
           IF(ALLOCATED(KLPC))DEALLOCATE(KLPC)
           IF(ALLOCATED(KRPC))DEALLOCATE(KRPC)           
  

           IF(ALLOCATED(KSUA))DEALLOCATE(KSUA)
           IF(ALLOCATED(KCOD))DEALLOCATE(KCOD)
           
           !IF(ALLOCATED(KLDN))DEALLOCATE(KLDN)
           !IF(ALLOCATED(KRDN))DEALLOCATE(KRDN)
           
           !IF(ALLOCATED(KLDP))DEALLOCATE(KLDP)
           !IF(ALLOCATED(KRDP))DEALLOCATE(KRDP)

           IF(ALLOCATED(KDNALG))DEALLOCATE(KDNALG)
           IF(ALLOCATED(KLNALG))DEALLOCATE(KLNALG)           
           
           IF(ALLOCATED(KLPN))DEALLOCATE(KLPN)
           IF(ALLOCATED(KRPN))DEALLOCATE(KRPN)
           

           IF(ALLOCATED(KLPP))DEALLOCATE(KLPP)
           IF(ALLOCATED(KRPP))DEALLOCATE(KRPP)
           
           IF(ALLOCATED(KLPALG))DEALLOCATE(KLPALG)
           IF(ALLOCATED(KDPALG))DEALLOCATE(KDPALG)
           
           
           IF(ALLOCATED(NTM))DEALLOCATE(NTM)
           
           IF(ALLOCATED(KRCOAG))DEALLOCATE(KRCOAG) !coagulation

           IF(ALLOCATED(FTMNL))DEALLOCATE(FTMNL)
           IF(ALLOCATED(FTHDR))DEALLOCATE(FTHDR)
           
           IF(ALLOCATED(RESPC))DEALLOCATE(RESPC)

           IF(ALLOCATED(DLALGC))DEALLOCATE(DLALGC)
           
           IF(ALLOCATED(ALGNH4))DEALLOCATE(ALGNH4)
           IF(ALLOCATED(ALGNO3))DEALLOCATE(ALGNO3)
		   
           IF(ALLOCATED(ALGDON))DEALLOCATE(ALGDON)
		   IF(ALLOCATED(ALGDON_OUT))DEALLOCATE(ALGDON_OUT)
		   IF(ALLOCATED(ALGDON_GL))DEALLOCATE(ALGDON_GL)
		   
		   IF(ALLOCATED(ALGPredDON))DEALLOCATE(ALGPredDON)
           IF(ALLOCATED(ALGPredDON_GL))DEALLOCATE(ALGPredDON_GL)		   
           IF(ALLOCATED(ALGPPDON))DEALLOCATE(ALGPPDON)
           IF(ALLOCATED(ALGPPDON_GL))DEALLOCATE(ALGPPDON_GL)	
		   
           IF(ALLOCATED(ALGPON))DEALLOCATE(ALGPON)
		   IF(ALLOCATED(ALGPON_OUT))DEALLOCATE(ALGPON_OUT)
           IF(ALLOCATED(ALGPON_GL))DEALLOCATE(ALGPON_GL)		   
		   
           IF(ALLOCATED(DENNO3))DEALLOCATE(DENNO3)
		   IF(ALLOCATED(DENNO3_OUT))DEALLOCATE(DENNO3_OUT)
           IF(ALLOCATED(DENNO3_GL))DEALLOCATE(DENNO3_GL)
		   
           !IF(ALLOCATED(MNLLDON))DEALLOCATE(MNLLDON)
           !IF(ALLOCATED(MNLRDON))DEALLOCATE(MNLRDON)           
           IF(ALLOCATED(HDRLPON))DEALLOCATE(HDRLPON)
           IF(ALLOCATED(HDRRPON))DEALLOCATE(HDRRPON)
           IF(ALLOCATED(HDRLPON_OUT))DEALLOCATE(HDRLPON_OUT)
           IF(ALLOCATED(HDRRPON_OUT))DEALLOCATE(HDRRPON_OUT)
           IF(ALLOCATED(HDRLPON_GL))DEALLOCATE(HDRLPON_GL)
           IF(ALLOCATED(HDRRPON_GL))DEALLOCATE(HDRRPON_GL)
           
           IF(ALLOCATED(NFIX))DEALLOCATE(NFIX)

           IF(ALLOCATED(ALGPO4))DEALLOCATE(ALGPO4)
           IF(ALLOCATED(ALGDOP))DEALLOCATE(ALGDOP)
           IF(ALLOCATED(ALGPOP))DEALLOCATE(ALGPOP)
           !IF(ALLOCATED(MNLLDOP))DEALLOCATE(MNLLDOP)
           !IF(ALLOCATED(MNLRDOP))DEALLOCATE(MNLRDOP)           
           IF(ALLOCATED(HDRLPOP))DEALLOCATE(HDRLPOP)
           IF(ALLOCATED(HDRRPOP))DEALLOCATE(HDRRPOP)
           
           IF(ALLOCATED(DORALG))DEALLOCATE(DORALG)
           IF(ALLOCATED(DOPR))DEALLOCATE(DOPR)
           IF(ALLOCATED(DCOD))DEALLOCATE(DCOD)
           IF(ALLOCATED(DDOC))DEALLOCATE(DDOC)
           
           IF(ALLOCATED(NITRIF))DEALLOCATE(NITRIF)
           
           IF(ALLOCATED(PSD))DEALLOCATE(PSD)
           IF(ALLOCATED(SAP))DEALLOCATE(SAP)

           IF(ALLOCATED(ALGUP))DEALLOCATE(ALGUP)
           IF(ALLOCATED(ALGRES))DEALLOCATE(ALGRES)

           IF(ALLOCATED(BENDOC))DEALLOCATE(BENDOC)
           IF(ALLOCATED(BENNH4))DEALLOCATE(BENNH4)
           IF(ALLOCATED(BENNO3))DEALLOCATE(BENNO3)
           IF(ALLOCATED(BENDON))DEALLOCATE(BENDON)
           IF(ALLOCATED(BENDOP))DEALLOCATE(BENDOP)
           IF(ALLOCATED(BENPO4))DEALLOCATE(BENPO4)
           IF(ALLOCATED(BENCOD))DEALLOCATE(BENCOD)
           IF(ALLOCATED(BENDO))DEALLOCATE(BENDO)
           IF(ALLOCATED(BENSA))DEALLOCATE(BENSA)
           IF(ALLOCATED(BENDEN))DEALLOCATE(BENDEN)
           IF(ALLOCATED(BENCH4G))DEALLOCATE(BENCH4G)
           IF(ALLOCATED(BENCH4A))DEALLOCATE(BENCH4A)

           
           IF(ALLOCATED(C1))DEALLOCATE(C1)
           IF(ALLOCATED(C2))DEALLOCATE(C2)
           IF(ALLOCATED(C2F))DEALLOCATE(C2F)
           IF(ALLOCATED(CSTAR))DEALLOCATE(CSTAR)
           IF(ALLOCATED(DTC))DEALLOCATE(DTC)
           IF(ALLOCATED(AC1))DEALLOCATE(AC1)
           
           IF(ALLOCATED(XFLUX)) DEALLOCATE (XFLUX)
           IF(ALLOCATED(XFLUX_ADV)) DEALLOCATE (XFLUX_ADV)
		   
           IF(ALLOCATED(QLIT))DEALLOCATE(QLIT)
           IF(ALLOCATED(CONLIT))DEALLOCATE(CONLIT)

           IF(ALLOCATED(AFI1))DEALLOCATE(AFI1)
           IF(ALLOCATED(AFI2))DEALLOCATE(AFI2)
           IF(ALLOCATED(AFI3))DEALLOCATE(AFI3)
           
           IF(ALLOCATED(ANL1))DEALLOCATE(ANL1)           
           IF(ALLOCATED(ANL2))DEALLOCATE(ANL2)
           IF(ALLOCATED(ANL3))DEALLOCATE(ANL3)
           
           IF(ALLOCATED(APL1))DEALLOCATE(APL1)
           IF(ALLOCATED(APL2))DEALLOCATE(APL2)
           IF(ALLOCATED(APL3))DEALLOCATE(APL3)

           IF(ALLOCATED(ASL1))DEALLOCATE(ASL1)
           IF(ALLOCATED(ASL2))DEALLOCATE(ASL2)
           IF(ALLOCATED(ASL3))DEALLOCATE(ASL3)

           IF(ALLOCATED(ACCHL1))DEALLOCATE(ACCHL1)
           IF(ALLOCATED(ACCHL2))DEALLOCATE(ACCHL2)
           IF(ALLOCATED(ACCHL3))DEALLOCATE(ACCHL3)
           
           IF(ALLOCATED(ANPP))DEALLOCATE(ANPP)
           
           IF(ALLOCATED(ARESP))DEALLOCATE(ARESP)
           IF(ALLOCATED(AKE))DEALLOCATE(AKE)
           IF(ALLOCATED(AGPP))DEALLOCATE(AGPP)
           

           
           IF(ALLOCATED(ABENDOC))DEALLOCATE(ABENDOC)
           IF(ALLOCATED(ABENNH4))DEALLOCATE(ABENNH4)
           IF(ALLOCATED(ABENNO3))DEALLOCATE(ABENNO3)
           IF(ALLOCATED(ABENPO4))DEALLOCATE(ABENPO4)
           IF(ALLOCATED(ABENCOD))DEALLOCATE(ABENCOD)
           IF(ALLOCATED(ABENDO))DEALLOCATE(ABENDO)
           IF(ALLOCATED(ABENSA))DEALLOCATE(ABENSA)
           
           IF(ALLOCATED(ASSFWS))DEALLOCATE(ASSFWS)
           IF(ALLOCATED(APCFWS))DEALLOCATE(APCFWS)
           IF(ALLOCATED(APNFWS))DEALLOCATE(APNFWS)
           IF(ALLOCATED(APPFWS))DEALLOCATE(APPFWS)
           IF(ALLOCATED(APSFWS))DEALLOCATE(APSFWS)
           IF(ALLOCATED(ACPIP))DEALLOCATE(ACPIP)
           IF(ALLOCATED(ACPOS))DEALLOCATE(ACPOS)
           
           IF(ALLOCATED(ADFEED))DEALLOCATE(ADFEED)
           
           IF(ALLOCATED(ABENCH4G))DEALLOCATE(ABENCH4G)
           IF(ALLOCATED(ABENCH4A))DEALLOCATE(ABENCH4A)

            !benthic algae
           IF(ALLOCATED(AFIB))DEALLOCATE(AFIB)
           IF(ALLOCATED(ANLB))DEALLOCATE(ANLB)
           IF(ALLOCATED(APLB))DEALLOCATE(APLB)
           
           IF(ALLOCATED(ANPPB))DEALLOCATE(ANPPB)
           
           IF(ALLOCATED(ABBM))DEALLOCATE(ABBM)!benthic algae
           IF(ALLOCATED(ABLITE))DEALLOCATE(ABLITE)
           IF(ALLOCATED(ABMB))DEALLOCATE(ABMB)
           IF(ALLOCATED(APB))DEALLOCATE(APB)
           IF(ALLOCATED(APRB))DEALLOCATE(APRB)
            
           !absolute values 
           IF(ALLOCATED(ABADOC))DEALLOCATE(ABADOC)
           IF(ALLOCATED(ABAPOC))DEALLOCATE(ABAPOC)
           IF(ALLOCATED(ABANH4))DEALLOCATE(ABANH4)
           IF(ALLOCATED(ABANO3))DEALLOCATE(ABANO3)
           IF(ALLOCATED(ABAPON))DEALLOCATE(ABAPON)
           IF(ALLOCATED(ABAPO4))DEALLOCATE(ABAPO4)
           IF(ALLOCATED(ABAPOP))DEALLOCATE(ABAPOP)
           IF(ALLOCATED(ABADO))DEALLOCATE(ABADO)
           
           IF(ALLOCATED(ACPOC))DEALLOCATE(ACPOC)
           IF(ALLOCATED(ACPON))DEALLOCATE(ACPON)
           IF(ALLOCATED(ACPOP))DEALLOCATE(ACPOP)
                              
            !Suspension feeder
           IF(ALLOCATED(ASFEED))DEALLOCATE(ASFEED)              
           IF(ALLOCATED(AJNSF))DEALLOCATE(AJNSF)
           IF(ALLOCATED(AJPSF))DEALLOCATE(AJPSF)
           IF(ALLOCATED(ASODSF))DEALLOCATE(ASODSF)
           IF(ALLOCATED(ASASF))DEALLOCATE(ASASF)
           IF(ALLOCATED(ASUSF))DEALLOCATE(ASUSF)
           
           IF(ALLOCATED(ASFGCIN))DEALLOCATE(ASFGCIN)
           IF(ALLOCATED(ASFCFEC))DEALLOCATE(ASFCFEC)
           IF(ALLOCATED(ASFCPSF))DEALLOCATE(ASFCPSF)
           
           IF(ALLOCATED(AFLXCSF))DEALLOCATE(AFLXCSF)
           IF(ALLOCATED(AFLXNSF))DEALLOCATE(AFLXNSF)
           IF(ALLOCATED(AFLXPSF))DEALLOCATE(AFLXPSF)
           
           IF(ALLOCATED(ARPOCSF))DEALLOCATE(ARPOCSF)
           IF(ALLOCATED(ARPONSF))DEALLOCATE(ARPONSF)
           IF(ALLOCATED(ARPOPSF))DEALLOCATE(ARPOPSF)
           
           IF(ALLOCATED(ASSISF))DEALLOCATE(ASSISF)
           IF(ALLOCATED(ASSISASF))DEALLOCATE(ASSISASF)
           IF(ALLOCATED(ASSISUSF))DEALLOCATE(ASSISUSF)
           IF(ALLOCATED(ASSIPSF))DEALLOCATE(ASSIPSF)
           
                                
           IF(ALLOCATED(A_T))DEALLOCATE(A_T)
           
           IF(ALLOCATED(AP1))DEALLOCATE(AP1)
           IF(ALLOCATED(ABM1))DEALLOCATE(ABM1)
           IF(ALLOCATED(APR1))DEALLOCATE(APR1)
           
           IF(ALLOCATED(AP2))DEALLOCATE(AP2)
           IF(ALLOCATED(ABM2))DEALLOCATE(ABM2)
           IF(ALLOCATED(APR2))DEALLOCATE(APR2)
           
           IF(ALLOCATED(AP3))DEALLOCATE(AP3)
           IF(ALLOCATED(ABM3))DEALLOCATE(ABM3)
           IF(ALLOCATED(APR3))DEALLOCATE(APR3)
           
           IF(ALLOCATED(AALGDOC))DEALLOCATE(AALGDOC)
           IF(ALLOCATED(AALGPOC))DEALLOCATE(AALGPOC)
           IF(ALLOCATED(ADENIT))DEALLOCATE(ADENIT)
           IF(ALLOCATED(AMNLDOC))DEALLOCATE(AMNLDOC)
           IF(ALLOCATED(AHDRPOC))DEALLOCATE(AHDRPOC)
           IF(ALLOCATED(AALGNH4))DEALLOCATE(AALGNH4)
           IF(ALLOCATED(AALGNO3))DEALLOCATE(AALGNO3)
           IF(ALLOCATED(AALGDON))DEALLOCATE(AALGDON)
           IF(ALLOCATED(AALGPON))DEALLOCATE(AALGPON)
           IF(ALLOCATED(ANT))DEALLOCATE(ANT)
           IF(ALLOCATED(ADENNO3))DEALLOCATE(ADENNO3)
           IF(ALLOCATED(AMNLDON))DEALLOCATE(AMNLDON)
           IF(ALLOCATED(AHDRPON))DEALLOCATE(AHDRPON)
           
           IF(ALLOCATED(ACFIX))DEALLOCATE(ACFIX)
           IF(ALLOCATED(ANFIX))DEALLOCATE(ANFIX)
           
           IF(ALLOCATED(AALGPO4))DEALLOCATE(AALGPO4)
           IF(ALLOCATED(AALGDOP))DEALLOCATE(AALGDOP)
           IF(ALLOCATED(AALGPOP))DEALLOCATE(AALGPOP)
           IF(ALLOCATED(AMNLDOP))DEALLOCATE(AMNLDOP)
           IF(ALLOCATED(AHDRPOP))DEALLOCATE(AHDRPOP)
           
           IF(ALLOCATED(APSD))DEALLOCATE(APSD)
           IF(ALLOCATED(ASAP))DEALLOCATE(ASAP)
           
           IF(ALLOCATED(AALGUP))DEALLOCATE(AALGUP)
           IF(ALLOCATED(AALGRES))DEALLOCATE(AALGRES)
           
           IF(ALLOCATED(ADO))DEALLOCATE(ADO)
           IF(ALLOCATED(ADORALG))DEALLOCATE(ADORALG)
           IF(ALLOCATED(ADOPR))DEALLOCATE(ADOPR)
           
           IF(ALLOCATED(ADCOD))DEALLOCATE(ADCOD)
           IF(ALLOCATED(ADDOC))DEALLOCATE(ADDOC)
           
           IF(ALLOCATED(ANITRIF))DEALLOCATE(ANITRIF)

           IF(ALLOCATED(V1SINGLE))DEALLOCATE(V1SINGLE)

                                 
           IF(ALLOCATED(SBN))DEALLOCATE(SBN)
           IF(ALLOCATED(BBN))DEALLOCATE(BBN)
           IF(ALLOCATED(HMSBV))DEALLOCATE(HMSBV)
           
            !Wen Long   IF(ALLOCATED(SFA))DEALLOCATE(SFA)

           IF(ALLOCATED(IWCMCB))DEALLOCATE(IWCMCB)       
           IF(ALLOCATED(IWCMNB))DEALLOCATE(IWCMNB)
           IF(ALLOCATED(IWCMPB))DEALLOCATE(IWCMPB)
           IF(ALLOCATED(IWCMSB))DEALLOCATE(IWCMSB)
           
           IF(ALLOCATED(WCMNB))DEALLOCATE(WCMNB)
           IF(ALLOCATED(WCMPB))DEALLOCATE(WCMPB)
           IF(ALLOCATED(WCMCB))DEALLOCATE(WCMCB)
           IF(ALLOCATED(WCMSB))DEALLOCATE(WCMSB)
           
           IF(ALLOCATED(DLWCMNB))DEALLOCATE(DLWCMNB)
           IF(ALLOCATED(DLWCMPB))DEALLOCATE(DLWCMPB)
           IF(ALLOCATED(DLWCKMNB))DEALLOCATE(DLWCKMNB)
           IF(ALLOCATED(DLWCKMCB))DEALLOCATE(DLWCKMCB)

           IF(ALLOCATED(FLUXS))DEALLOCATE(FLUXS)
             
           IF(ALLOCATED(S1FLXCB))DEALLOCATE(S1FLXCB)
           IF(ALLOCATED(S1FLXNB))DEALLOCATE(S1FLXNB)
           IF(ALLOCATED(S1FLXPB))DEALLOCATE(S1FLXPB)
           
           IF(ALLOCATED(S2FLXCB))DEALLOCATE(S2FLXCB)
           IF(ALLOCATED(S2FLXNB))DEALLOCATE(S2FLXNB)
           IF(ALLOCATED(S2FLXPB))DEALLOCATE(S2FLXPB)
           
           IF(ALLOCATED(S3FLXCB))DEALLOCATE(S3FLXCB)
           IF(ALLOCATED(S3FLXNB))DEALLOCATE(S3FLXNB)
           IF(ALLOCATED(S3FLXPB))DEALLOCATE(S3FLXPB)

           
            
           !IF(ALLOCATED(JPOC))DEALLOCATE(JPOC)           
           !IF(ALLOCATED(JPON))DEALLOCATE(JPON)
           !IF(ALLOCATED(JPOP))DEALLOCATE(JPOP)
           !IF(ALLOCATED(JPOS))DEALLOCATE(JPOS)
           
           !IF(ALLOCATED(CPOC))DEALLOCATE(CPOC)
           !IF(ALLOCATED(CPON))DEALLOCATE(CPON)             
           !IF(ALLOCATED(CPOP))DEALLOCATE(CPOP)

           !IF(ALLOCATED(CPOS))DEALLOCATE(CPOS)
           !IF(ALLOCATED(CPO4))DEALLOCATE(CPO4)
           !IF(ALLOCATED(CNO3))DEALLOCATE(CNO3)
           !IF(ALLOCATED(CNH4))DEALLOCATE(CNH4)
           

           !IF(ALLOCATED(CCH4))DEALLOCATE(CCH4)
           !IF(ALLOCATED(CSO4))DEALLOCATE(CSO4)
           !IF(ALLOCATED(CHS))DEALLOCATE(CHS)
           !IF(ALLOCATED(CSI))DEALLOCATE(CSI)           

           IF(ALLOCATED(BSVOL))DEALLOCATE(BSVOL)
!           IF(ALLOCATED(HSED))DEALLOCATE(HSED)

            IF(ALLOCATED(CTEMP))DEALLOCATE(CTEMP)
           
           IF(ALLOCATED(PCFWS))DEALLOCATE(PCFWS)      
           IF(ALLOCATED(PNFWS))DEALLOCATE(PNFWS)
           IF(ALLOCATED(PPFWS))DEALLOCATE(PPFWS)
           IF(ALLOCATED(PSFWS))DEALLOCATE(PSFWS)
           IF(ALLOCATED(SSFWS))DEALLOCATE(SSFWS)
       
           IF(ALLOCATED(AASRAT))DEALLOCATE(AASRAT)


           !----
             
           IF(ALLOCATED(BENDOCB))DEALLOCATE(BENDOCB)
           IF(ALLOCATED(BENNH4B))DEALLOCATE(BENNH4B)
           IF(ALLOCATED(BENNO3B))DEALLOCATE(BENNO3B)
           IF(ALLOCATED(BENPO4B))DEALLOCATE(BENPO4B)
           IF(ALLOCATED(BENCODB))DEALLOCATE(BENCODB)
           IF(ALLOCATED(BENDOB))DEALLOCATE(BENDOB)
           IF(ALLOCATED(BENSAB))DEALLOCATE(BENSAB)

           IF(ALLOCATED(ATMFLXCB))DEALLOCATE(ATMFLXCB)           
           IF(ALLOCATED(ATMFLXNB))DEALLOCATE(ATMFLXNB)
           IF(ALLOCATED(ATMFLXPB))DEALLOCATE(ATMFLXPB)

           IF(ALLOCATED(BENFLXPCB))DEALLOCATE(BENFLXPCB)

           IF(ALLOCATED(BENFLXDNB))DEALLOCATE(BENFLXDNB)
           IF(ALLOCATED(BENFLXPNB))DEALLOCATE(BENFLXPNB)
           
           
           IF(ALLOCATED(BENFLXDPB))DEALLOCATE(BENFLXDPB)           
           IF(ALLOCATED(BENFLXPPB))DEALLOCATE(BENFLXPPB)

           IF(ALLOCATED(DLSEDKNB))DEALLOCATE(DLSEDKNB)
           IF(ALLOCATED(DLSEDKCB))DEALLOCATE(DLSEDKCB)
           
           IF(ALLOCATED(BURIALFLXCB))DEALLOCATE(BURIALFLXCB)           
           IF(ALLOCATED(BURIALFLXNB))DEALLOCATE(BURIALFLXNB)
           IF(ALLOCATED(BURIALFLXPB))DEALLOCATE(BURIALFLXPB)


        ! KURT GLAESEMANN add C2_GL
           IF(ALLOCATED(C2_GL))DEALLOCATE(C2_GL)
        ! RGl added below for output of CCHL1_GL
           IF(ALLOCATED(CCHL1_GL))DEALLOCATE(CCHL1_GL)
           IF(ALLOCATED(T_GL))DEALLOCATE(T_GL)
           IF(ALLOCATED(S_GL))DEALLOCATE(S_GL)
           IF(ALLOCATED(total_netPP_GL))DEALLOCATE(total_netPP_GL)

        !Wen Long added D_GL, EL_GL, H_GL (total depth, surface elevation and bathymetric depth)
           IF(ALLOCATED(D_GL))DEALLOCATE(D_GL)
           IF(ALLOCATED(H_GL))DEALLOCATE(H_GL)
           IF(ALLOCATED(EL_GL))DEALLOCATE(EL_GL)

        !--Wen Long debugging benthic fluxes ---
           IF(ALLOCATED(BFLUX))DEALLOCATE(BFLUX)           
           IF(ALLOCATED(BFLUX_GL))DEALLOCATE(BFLUX_GL)

           IF(ALLOCATED(BFLUXNX))DEALLOCATE(BFLUXNX)
           IF(ALLOCATED(BFLUXNX_GL))DEALLOCATE(BFLUXNX_GL)

           IF(ALLOCATED(BFLUXB_GL))DEALLOCATE(BFLUXB_GL)
           
           
        !--Wen Long deallocating global hydro variables here ---
           
           IF(ALLOCATED(UL_GL))DEALLOCATE(UL_GL)
		   IF(ALLOCATED(WW_GL))DEALLOCATE(WW_GL)
		   IF(ALLOCATED(WW_OUT))DEALLOCATE(WW_OUT)
           IF(ALLOCATED(VL_GL))DEALLOCATE(VL_GL)
           IF(ALLOCATED(WTSL_GL))DEALLOCATE(WTSL_GL)
           IF(ALLOCATED(KHL_GL))DEALLOCATE(KHL_GL)
           IF(ALLOCATED(SL_GL))DEALLOCATE(SL_GL)
           IF(ALLOCATED(TL_GL))DEALLOCATE(TL_GL)
           IF(ALLOCATED(ELL_GL))DEALLOCATE(ELL_GL)
           IF(ALLOCATED(DTFAL_GL))DEALLOCATE(DTFAL_GL)
           IF(ALLOCATED(UARD_OBCN_GL))DEALLOCATE(UARD_OBCN_GL)
           IF(ALLOCATED(XFLUX_OBC_GL))DEALLOCATE(XFLUX_OBC_GL)
         

        !--Wen Long  added the following algae debugging global variables here ---

           IF(ALLOCATED(P1_GL))DEALLOCATE(P1_GL)
           IF(ALLOCATED(P2_GL))DEALLOCATE(P2_GL)
           IF(ALLOCATED(BM1_GL))DEALLOCATE(BM1_GL)
           IF(ALLOCATED(BM2_GL))DEALLOCATE(BM2_GL)
           IF(ALLOCATED(NL1_GL))DEALLOCATE(NL1_GL)
           IF(ALLOCATED(NL2_GL))DEALLOCATE(NL2_GL)
           IF(ALLOCATED(PL1_GL))DEALLOCATE(PL1_GL)
           IF(ALLOCATED(PL2_GL))DEALLOCATE(PL2_GL)
           IF(ALLOCATED(B1SZ_GL))DEALLOCATE(B1SZ_GL)
           IF(ALLOCATED(B2SZ_GL))DEALLOCATE(B2SZ_GL)
           IF(ALLOCATED(B1LZ_GL))DEALLOCATE(B1LZ_GL)
           IF(ALLOCATED(B2LZ_GL))DEALLOCATE(B2LZ_GL)
           IF(ALLOCATED(PR1_GL))DEALLOCATE(PR1_GL)
           IF(ALLOCATED(PR2_GL))DEALLOCATE(PR2_GL)
           
           IF(LIGHT_EXTINCTION)THEN
                IF(ALLOCATED(IAVG_GL))DEALLOCATE(IAVG_GL)
				IF(ALLOCATED(KDLAMB_int_GL)) DEALLOCATE(KDLAMB_int_GL)
                IF(ALLOCATED(FI1_GL))DEALLOCATE(FI1_GL)
                IF(ALLOCATED(FI2_GL))DEALLOCATE(FI2_GL)
           ENDIF
           
		   IF(ALLOCATED(SSI_GL)) DEALLOCATE(SSI_GL)
		   
           IF(ALLOCATED(IK1_GL))DEALLOCATE(IK1_GL)
           IF(ALLOCATED(IK2_GL))DEALLOCATE(IK2_GL)
           IF(ALLOCATED(IK1))DEALLOCATE(IK1)
           IF(ALLOCATED(IK2))DEALLOCATE(IK2)
		   
		   
		   ! B CLark for SSI respuspension
		   
		   IF(ALLOCATED(ShearStress)) Deallocate(ShearStress)
			IF(ALLOCATED(ShearStress_GL)) Deallocate(ShearStress_GL)
			IF(ALLOCATED(RHO1)) Deallocate(RHO1)
			IF(ALLOCATED(RHO1_GL)) Deallocate(RHO1_GL)
		   

        !---Wen Long added the following for debugging DOXG

            IF(ALLOCATED(DDOC_GL))DEALLOCATE(DDOC_GL)
            IF(ALLOCATED(DCOD_GL))DEALLOCATE(DCOD_GL)
            
            IF(ALLOCATED(NITRIF_GL))DEALLOCATE(NITRIF_GL)
            IF(ALLOCATED(DOSZ_GL))DEALLOCATE(DOSZ_GL)
            IF(ALLOCATED(DOLZ_GL))DEALLOCATE(DOLZ_GL)
            IF(ALLOCATED(DOPR_GL))DEALLOCATE(DOPR_GL)
            IF(ALLOCATED(DORALG_GL))DEALLOCATE(DORALG_GL)
            IF(ALLOCATED(RESP_GL))DEALLOCATE(RESP_GL)
            IF(ALLOCATED(FTCOD_GL))DEALLOCATE(FTCOD_GL)
            IF(ALLOCATED(COD_GL))DEALLOCATE(COD_GL)
            IF(ALLOCATED(REAERDO_GL))DEALLOCATE(REAERDO_GL)
            IF(ALLOCATED(REAERDO))DEALLOCATE(REAERDO)

        !---Wen Long declared a large buffer (so that we do not have to dynamically allocate and deallocate buffers 
        !   all the time--
#if defined (MULTIPROCESSOR)
            IF(ALLOCATED(RCVBUF))DEALLOCATE(RCVBUF)
            IF(ALLOCATED(SNDBUF))DEALLOCATE(SNDBUF)
#endif
          Write(*,*)'Allocated WQM vars' ! B Clark Debug

            RETURN

        END SUBROUTINE WQM_DEALLOC
        
END MODULE MOD_WQM


