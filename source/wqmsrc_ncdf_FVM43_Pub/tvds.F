!Subroutine TVDS()  !Time varying data series

!************************************************************************
!**                  S U B R O U T I N E   T V D S                     **
!************************************************************************

    SUBROUTINE TVDS (NXTVD)
	
	USE MOD_PREC, ONLY: SP
	USE MOD_LIMS, ONLY: MLOC
	
	USE MOD_SIZES, ONLY:	&
				MGL,		&!
				NCP!,       &!
						
    USE MOD_WQM, ONLY: 			&
				BB,             &!
                BENDOC,         &!
                BENTHIC_FLUXES, &!

                DIAGNOSTICS,    &!
                FD,             &!
                FILGTH,         &!
                IT,             &!
                ITNX,           &!
                JDAY,           &!
                KT,             &!
                METFN,          &!
                METPTR,         &!
                NXBFI,          &!
                NXMET,          &!
                TE,             &!
                TMEND,          &!
                TTSS,           &!
                WMS,            &!
				BENCOD,         &!
                BENCODB,        &!
                BENDO,          &!
                BENDOB,         &!
                BENDOCB,        &!
                BENDON,         &!
                BENDOP,         &!
                BENNH4,         &!
                BENNH4B,        &!
                BENNO3,         &!
                BENNO3B,        &!
                BENPO4,         &!
                BENPO4B,        &!
                BENSA,          &!
                BENSAB,         &!
                BFIFN,          &!
                BFIPTR,         &!

                KSDOC,          &!
				ATMFN,          &!
                ATMLDON,        &!
                ATMLDOP,        &!
                ATMNH4,         &!
                ATMNO3,         &!
                ATMOS_LOADS,    &!
                ATMPO4,         &!
                ATMPTR,         &!
                ATMRDON,        &!
                ATMRDOP,        &!
                KHSO,           &!
                KSNH4,          &!
                KSNO3,          &!
                KSO,            &!
                KSPO4,          &!
                KSSA,           &!
                MTCNO3,         &!
                NXATM,          &!
                PRECIP,         &!
                SAV_LOADS,      &!
                SEDNO3,         &!
                TRSDOC,         &!
                TRSNH4,         &!
                TRSNO3,         &!
                TRSO,           &!
                TRSPO4,         &!
                TRSSA,          &!
				 B,              &!
                NBB,            &!
                NS1P,           &!
                NS2P,           &!
                NS3P,           &!
                NXSAV,          &!
                SAVPTR,         &!
                SVIFN!,          &!

	USE MOD_HYDROVARS, ONLY: ART1
    USE MOD_CONTROL, ONLY : 		&
			SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			!,MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
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
			
        USE MOD_FILEINFO, ONLY : 	&
				DIA				&!
				!,CBC 			&!
				!,S1			&!
				!,S2			&!
				!,S3			& !                 
				,BFI			&!
				,BAI            &!
				,MET			&!
				!,BFO			& !
				!,KEI			&!
				,ATM			&!
				!,STL			& !
				!,AGR			& !
				,SVI	!		& !
				!,SVO			& !
				!,KFL			& !
				!,ZOO			& !
				!,ZFO			& !
				!,ALO      		&!
				!,CON			&!
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

    USE MOD_SAV, ONLY	:	&
			SPNAM,   		&	!species name of the SAV
            EPINAM,  		&   !name of the epiphyte
			NSAVM,   		&  	!global maximum number of SAV species in all of tracer control elements (TCE's) - polygon
                                !surrounding a given node
			
			KHSTOX,         &   !Sulfide toxicity coefficient in calculating HS^- limitation on SAV growth (m^3/mole-S2)
            KESAV,      &   !Effect on light attenuation by SAV leaf and stem, KESAV (gC/m^3)  
            KEEPI,      &   !Shading effect on light on SAV shoots (leaf) due to Epiphytes ((m^2-leaf-surface/gDW-epiphyte)
            KHEP,       &  !half saturation desnity of epiphyte (gC-Epiphyte/gC-leaf) for epiphyte growth on SAV leaf
            NLIMEPI,	&  !Never used ???                      
             
			KHNLEAF,    &  !Half saturation constant for N uptake by leaf (shoots) (gN/m^3)
            KHNROOT,    &  !Half saturation constant for N uptake by root (gN/m^3)
            KHPLEAF,    &  !Half saturation constant for P uptake by leaf (gP/m^3)
            KHPROOT,    &  !Half saturation constant for P uptake by root (gP/m^3)
            KHNEPI,     &  !Half saturation constant for N uptake by epiphytes (gN/m^3)
            KHPEPI,     &  !Half saturation constant for P uptake by epiphytes (gN/m^3)
            KHNPSAV,    &  !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for SAV
            KHNPEPI,    &  !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for epiphytes

			PMSAV,      &   !maximum production rate of SAV, gC/gDW/day  where DW is dry weight
            PRSPSAVB,   &   !fraction of leaf production rate that is counted as photochemical respiration at salt << saltmax
                            !at salt >> saltmax, PRSPSAV approaches 1, and there will be no production of SAV leaf hence
                            !no photosynthesis for SAV growth
            BMSAV,      &   !Basal metabolism rate of SAV (1/day)
            BMTBRREF,   &   !Reference basal metabolism rate of tube [1/day] BMTUBE=BMTBRREF*f(T), i.e. function of temperature
            FDOSR,      &   !Fraction of DO production discounted due to route of Leaf production to root
            SALMAX,     &   !Salinity at which PRSPSAV is approaching 1 (maximum) from PRSPSAVB, [psu]
            ALPHSAV,    &   !Initial slope of P-I curve for SAV [ (gC/gDW)/(E/m^2) ]
            ALAC,       &   !m^2-leaf-area/gC-leaf (leaf area to leaf carbon ratio) 
                            !Note in 2002 version, it is the inverse ACLA: gC-leaf/m^2-leaf-area
            SLSAV,      &   !???? Never used !Sloughing rate of SAV ?  1/day
            ACDWSAV,    &   !Carbon to Dry Weight ratio of SAV (gC/gDW)
            ANDWSAV,    &   !Nitrogen to Dry Weight ratio of SAV (gN/gDW)
            APDWSAV,    &   !Phosphorus to Dry Weight ratio of Sav (gP/gDW)
            ANCSAV,     &   !Nitrogen to Carbon ratio of SAV (gN/gC)             
            APCSAV,     &   !P to C ratio of SAV (gP/gC)
            HCAN,       &   !Canopy height (m), HCAN = ACAN + BCAN *( leaf+stem )
            ACAN,       &   !basic canopy height (m)
            BCAN,       &	 !coefficient for calculating leaf (gC/m^2) and stem (gC/m^2) 
                             !contribution to canopy height,  [m^3/gC]

			FNISAV,   	&   !Fraction of NH4 Nitrogen generated due to SAV (leaf, stem etc) respiration (mortality)
            FNLDSAV,  	&   !Fraction of Nitrogen routed to LDON on SAV mortality
            FNRDSAV,  	&   !Fraction of Nitrogen routed to RDON on SAV mortality
            FNLPSAV,  	&   !Fraction of Nitrogen routed to LPON on SAV mortality
            FNRPSAV,  	&   !Fraction of Nitrogen routed to RPON on SAV mortality
            FPISAV,   	&   !Fraction of Phosphorus routed to PO4 in water column due to SAV (leaf,stem) respiration
            FPLDSAV,  	&   !Fraction of P routed to LDOP in water column due to SAV mortality
            FPRDSAV,  	&   !Fraction of P routed to RDOP in water column due to SAV mortality
            FPLPSAV,  	&   !Fraction of P routed to LPOP in water column due to SAV mortality
            FPRPSAV,  	&   !Fraction of P routed to RPOP in water column due to SAV mortality
            FCLDSAV,  	&   !Fraction of Carbon routed to LDOC in water column due to SAV mortality
            FCRDSAV,  	&   !Fraction of Carbon routed to RDOC in water column due to SAV mortality
            FCLPSAV,  	&   !Fraction of Carbon routed to LPOC in water column due to SAV mortality
            FCRPSAV,  	&   !Fraction of Carbon routed to RPOC in water column due to SAV mortality
            FDOSAV,     &   !Fraction of Carbon routed to CO2 (oxygen consumption) in water column due to SAV mortality
                            !Note FDOSAV + FCLDSAV+ FCRDSAV + FCLPSAV + FCRPSAV = 1
                            !Wen Long: we really need to rename FDOSAV to FCO2SAV and make it equivalent to FDOSAV!

			WSSSAV,   	&   !Increase of fixed solids (ISS) settling by SAV leaf and stem, (m/day)/(gC/m^2)==>m^3/gC/day
            WSLSAV,   	&   !Increase of LPOM settling by SAV leaf and stem, m^3/gC/day
            WSRSAV,   	&   !Increase of RPOM settling by SAV leaf and stem, m^3/gC/day
            WS1SAV,   	&   !Increase of Algae group 1 settling by SAV leaf and stem, m^3/gC/day
            WS2SAV,   	&   !                        2
            WS3SAV,   	&   !                        3
            WSUSAV,   	&   !Increase of particulate biogenic silica settling rate by SAV leaf and stem, m^3/gC/day
                     
			NDDI,     	& 	!number of degree days for entries of look up table of FPLEAF etc for Senescence of leaf etc
            NSAVCELL, 	& 	!total number of grid cells (ICM boxes) that have SAV --> need to be changed to nodes of FVCOM
            NSAVSPCM, 	&   !maximum number of species in a cell
              
			PMEPI,      & 	!maximum photosynthesis production rate (gC/gCHLA/day) of epiphytes
            BMEPI,      & 	!basal metabolism loss rate of epiphyte (1/day)
            PREPI,      & 	!predation rate on epiphytes (gC-leaf/gC-epiphyte/day)
            PRSPEPI,    & 	!Photochemical respiration rate of epiphyte --> BMEP

			ALPHEPI,    & 	!slope of PI curve for epiphyte photothenesis (gC/gDW)/(E/m^2)
            CCHLEPI,    & 	!Carbon to Chlorophyll ratio of epiphyte (gC/gCHLA)
            ANCEPI,     & 	!N to C ratio of epiphytes (gN/gC)
            APCEPI,     & 	!P to C ratio of epiphytes (gP/gC)
            ADWCEPI,    &  	!ratio of total epiphyte dry weight to viable carbon (gDW/gC)
             
			FCLDEPI,    &   !Fraction of C routed to LDOC due to basal metabolism of epiphyte
            FCRDEPI,    &   !Fraction of C routed to RDOC due to basal metabolism of epiphyte
                            !WLong: missing FCLPEPI and FCRPEPI to routt to LPOC and RPOC ??? (seems not coded)
                            !FCLPEPI and FCRPEPI should be added (Wen Long)
                            !also missing FDOEPI compared to FDOSAV for SAV
                            !code assumes that no carbon will be routed to LPOP and RPOP due to metabolism
                            !and also FDOEPI=1-FCLDEPI-FCRDEPI, and FDOEPI is only implicitly used in calculating 
                            !DOEPIW
             
            FNIEPI,     &   !Fraction of N routed to NH4  due to epiphyte basal metabolism
            FNLDEPI,    &   !Fraction of N routed to LDON due to epiphyte basal metabolism
            FNRDEPI,    &   !Fraction of N routed to RDON due to epiphyte basal metabolism
            FNLPEPI,    &   !Fraction of N routed to LPON due to epiphyte basal metabolism
            FNRPEPI,    &   !Fraction of N routed to RPON due to epiphyte basal metabolism
            FPIEPI,     &   !Fraction of P routed to  PO4 due to basal metabolism of epiphyte
            FPLDEPI,    &   !Fraction of P routed to LDOP due to basal metabolism of epiphyte
            FPRDEPI,    &   !Fraction of P routed to RPOP due to basal metabolism of epiphyte
            FPLPEPI,    &   !Fraction of P routed to LPOP due to basal metabolism of epiphyte
            FPRPEPI,	&   !Fraction of P routed to RPOP due to basal metabolism of epiphyte             

			FCLDPEP,    &   !Fraction of C routed to LDOC due to predation and sloughing on epiphytes
            FCRDPEP,    &   !Fraction of C routed to RDOC due to predation and sloughing on epiphytes
            FCLPPEP,    &   !Fraction of C routed to LPOC due to predation and sloughing on epiphytes
            FCRPPEP,    &   !Fraction of C routed to RPOC due to predation and sloughing on epiphytes
            FNIPEP,     &   !Fraction of Nitrogen routed to NH4  due to predation on epiphytes and Slough
            FNLDPEP,    &   !Fraction of Nitrogen routed to LDON due to predation on epiphytes and Slough
            FNRDPEP,    &   !Fraction of Nitrogen routed to RDON due to predation on epiphytes and Slough
            FNLPPEP,    &   !Fraction of Nitrogen routed to LPON due to predation on epiphytes and Slough
            FNRPPEP,    &   !Fraction of Nitrogen routed to RPON due to predation on epiphytes and Slough
            FPIPEP,     &   !Fraction of P routed to PO4  due to predation on epiphytes and slough
            FPLDPEP,    &   !Fraction of P routed to LDOP due to predation on epiphytes and slough
            FPRDPEP,    &   !Fraction of P routed to RDOP due to predation on epiphytes and slough
            FPLPPEP,    &   !Fraction of P routed to LPOP due to predation on epiphytes and slough
            FPRPPEP,    &   !Fraction of P routed to RPOP due to predation on epiphytes and slough            

			FPSR,      	&	!Never used!!!! This is supposed to be the fraction of leaf production routed to root
							!(replaced by FPROOT)

     
			FPLEAF, 	& 	!fraction of production assigned to leaf growth
            FPROOT, 	& 	!fraction of production assigned to root growth
            FPSTEM, 	& 	!fraction of production assigned to stem growth
            FPTUBER,	& 	!fraction of production assigned to tuber growth
                            !WLong: FPLEAF +FPROOT+FPSTEM+FPTUBER = 1
                            !all growth is based on photosynthesis by leaf
                            !(so stem and tuber have no photosynthesis?? I thought they look green )
            TRTBRLF,	&   !transfer rate from tuber to leaf (gC-Leaf/day/gC-Tuber)
                               
			!Look up table of temperature control in a year on growth and mortality rates of SAV and epiphytes
			FTPSAV, 	& 	!temperature control on SAV growth
            FTRSAV, 	& 	!temperature control on SAV resipration 
			FTPEP,  	& 	!temperature control on epiphyte growth rate
            FTREP,  	& 	!temperature control on epiphyte resipration (basal metabolism)
            FTPREP,		&   !temperature control on predation on epiphyte
                               
			!fluxes to sediments by SAV mortality
			FRPOCSAV, 	&	!
			FRPONSAV, 	&	!
			FRPOPSAV,   &	!3G partitioning of POC, PON, POP sources to sediments from SAV (non-dimensional)
  
			!State variables of SAV model 
			LEAF,   	& 	!Leaf biomass (gC/m^2) 
            ROOT,   	& 	!Root biomass (gC/m^2)
            STEM,   	& 	!Stem biomass (gC/m^2)
            TUBER,  	& 	!Tuber biomass (gC/m^2)
            EP,			&   !Epiphytes biomass (gC-epiphytes/gm^2-leaf)
							!Note in 2002 version of ICM, EP's unit was (gC-epiphytes/gC-leaf)
                                 
			SAVAREA,	&   !area of SAV coverage for each surface cell (m^2) --> needs to be changed for FVCOM grid                                 
			SAVCELL,	&   !grid cell number for all cells that have SAV in it
							!size 1xNSAVCELL --> need to change to element number or node number!!
			
			NSAVSPC,	&	!(:) !number of SAV species in each grid cell
							!size 1xNSAVCELL

			SAVDPH,  	&  	!Depth of SAV (mean depth of the plant) for all species in each grid cell
							!size, NSAVCELL x NSAVSPC(B) where B is cell number, B=1,...,NSAVCELL
            SAVFRAC,	&  	!coverage for each SAV cell and each vertical deph increment
							!i.e. coverage is function of grid cell and SAV species

							!Growth and mortality rates of SAV
            PLEAF,		&	!Leaf growth rate (1/day), P (SAV shoots growth) in C-M 2001 paper's SH equation
            BMLEAF, 	&  	!respiration rate of leaf of SAV (1/day) with temperature control 
            BMSTEM, 	&  	!resipration rate of stem of SAV (1/day) with temperature control
            BMROOT, 	&  	!respiration rate of root of SAV (1/day) with temperature control
            BMTUBER,    &	!respiration rate of tuber of SAV (1/day) with temperature control                                      
     
			!Growth and mortality rates of SAV                                      
            PEP,   		&  	!epiphytes growth rate (1/day)  = PMEPI*light_limitation*deinsity_limitation/CCHLEPI
            BMEP,  		&  	!epiphytes basal metabolism rate (1/day) (also including photochemical respiration)
            PREP,  		&   !predation rate on epiphytes (1/day) = PREPI*EP
                            !where PREPI has unit gC-shoot/gC-epiphyte/day

            SLSH,		&	!(:)      !epiphyte sloughing rate (1/day)
     
			!Nutrient limitations on SAV and epiphytes growth
            NLSAV,    	&  	!Nitrogen limiation on SAV growth (non-dimensional)
            PLSAV,    	&  	!Phosphorus limitation on SAV growth (non-dimensional)
            FNSEDSAV, 	&  	!fraction of nitrogen uptake by SAV due to nitrogen in sediment pore water
							!Wen Long: I think we should further split into NH4 and NO3 
							!          uptake by SAV from sediments
            FPSEDSAV, 	&  	!fraction of phosphorus uptake by SAV due to P in sediment pore water
            NLEPI,    	&  	!Nitrogen limitation on epiphytes growth (non-dimensional)
            PLEPI,    	&  	!Phosphorus limitation on epiphytes growth (non-dimensional)
            FHS,		&   !fresh water SAV species growth limitation due to sulfide HS^- 
							!(in sediment anerobic layer) toxicity (dimensionless)

			EPATN,  	& 	!Light attenuation by epiphyte
							!e^(-Kep*Acla*Adwcep*EP) in Cerco-Moore 2001 paper for Ish eqn
							!e^(-Kep*Adwcep*EP) in new version as EP changed to gC-Epi/m^2-leaf unit
			SAVATN, 	& 	!Light attenuation by SAV (shoots)
							![1-e^(-Ksh*SH+KESS*HCAN)]/(Ksh*SH+KESS*HCAN) in C-M 2001 Iwc eqn
			WATATN, 	& 	!water column attenuation of light(dimensionless) at canopy top
							!= exp(-KESS*ZTOCNPY) where ZTOCNPY is distance from surface
							!                           to canopy top
			FISH,   	& 	!light limitation function for LEAF (shoots) f(I)=Ish/sqrt(Ish^2+Ik^2)
            FIEP,   	& 	!light limitation function for epiphytes, f(I)=Iep/sqrt(Iep^2+Ik^2)
            NPPSAV, 	& 	!net leaf production per day (growth-basal metabolism) (gC/m^2/day)
            NPPEPI,		&   !net epiphyte production per day (growth - basal metabolism) (gC-epi/gC-leaf/day)

			!fluxes to water column
			DOSAVW, 	& 	!DO source (flux) due to SAV (source-sink) (gO2/m^2/day ) 
							!DOSAVW --> DTDO (gO2/m^3/sec)--> DTC (gO2/m^3/sec)
							!-->DTM(gO2/m^3/sec)-->DTM*DLT-->C1(:,:,27)(gO2/m^3) !27'th constituents is DO
			LDOCSAVW, 	& 	!LDOC flux due to SAV (gC/m^2/day)
			
!WL 		!***ALSO used by TVDS.F ****
							!LDOCSAVW --> DTDOC (gC/m^3/sec) -> DTC --> DTM --> DTM*DLT --> C1(:,:,9)
			RDOCSAVW, 	&  	!RDOC flux due to SAV (gC/m^2/day)                                 
			LPOCSAVW, 	&  	!LPOC flux due to SAV (gC/m^2/day) 
!WL                                              !***ALSO used by TVDS.F ****
		    RPOCSAVW, 	&  	!RPOC flux due to SAV (gC/m^2/day)
!WL                                              !***ALSO used by TVDS.F ****
		    NH4SAVW,  	&  	!NH4 flux to water column (respiration source - uptake sink by SAV growth) (gN/m^2/day)
							!NH4SAVW (gN/m^2/day) --> DTNH4  (gN/m^3/sec) --> DTC --> DTM --> C1(:,:,13)
		    NO3SAVW,  	&  	!NO3 flux to water column (due to uptake by SAV, positive increasing NO3 in water) (gN/m^2/day)
		    LDONSAVW, 	&  	!LDON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf, basal metabolism of stem) (gN/m^2/day)
		    RDONSAVW, 	&  	!RDON flux                                     
		    LPONSAVW, 	&  	!LPON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem (gN/m^2/day)
		    RPONSAVW, 	&  	!RPON flux 
		    PO4SAVW,  	&  	!PO4 flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem and uptake by leaf) (gP/m^2/day)
		    LDOPSAVW, 	&  	!LDOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day)
		    RDOPSAVW, 	&  	!RDOP flux (gP/m^2/day)
		    LPOPSAVW, 	&  	!LPOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day) 
		    RPOPSAVW, 	&  	!RPOP flux (gP/m^2/day)
		      DOEPIW, 	&  	!DO   flux (positive into water column) due to epiphytes (gO2/m^2/day)
		    LDOCEPIW, 	&  	!LDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
		    RDOCEPIW, 	&  	!RDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
		    LPOCEPIW, 	&  	!LPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
		    RPOCEPIW, 	&  	!RPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
		     NH4EPIW, 	&  	!NH4  flux (positive into water column) due to epiphytes (gN/m^2/day) (basal meabolism + predation + sloughing - uptake by photosynthesis)
		     NO3EPIW, 	&  	!NO3  flux
		    LDONEPIW, 	&  	!LDON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
		    RDONEPIW, 	&  	!RDON flux (positive into water column) (gN/m^2/day)
		    LPONEPIW, 	&  	!LPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
		    RPONEPIW, 	&  	!RPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
		    PO4EPIW,  	&  	!PO4  flux (positive into water column) due to epiphytes mortality (basal metabolism, salt toxicity, sloughing) and uptake (sink) (growth) (gP/m^2/day)                                    
		    LDOPEPIW, 	&  	!LDOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
		    RDOPEPIW, 	&  	!RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
		    LPOPEPIW, 	&  	!LPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
		    RPOPEPIW, 	&	!RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)

			!fluxes to sediments                                                                         
		    SEDPOCSAV, 	&  	!POC flux to sediments (positive into sediments) due to SAV (gC/m^2/day)
							!SEDPOCSAV (gC/m^2/day) --> JPOC (mgC/m^2/day) -> POC1, POC2, POC3 (source of particulate organic carbon in sediments)
		    SEDPONSAV, 	&  	!PON flux to sediments (positive into sediments) due to SAV (gN/m^2/day)
							!SEDPONSAV (gN/m^2/day) --> JPON (mgN/m^2/day) -> PON1, PON2, PON3 (source of particulate organic nitrogen in sediments)
		    SEDPOPSAV, 	&  	!POP flux to sedimetns (positive into sediments) due to SAV (gP/m^2/day)
							!SEDPOPSAV (gP/m^2/day) --> JPOP (mgP/m^2/day) -> POP1, POP2, POP3 (source of particulate organic phosphorus in sediments)
		    SEDDOSAV,  	&  	!DO flux to water column (positive into water) due to SAV growth in sediments (source) and mortality of SAV in sediments (tuber and root basal metabolism)(sink)
							!(gO2/m^2/day) (WLong: does this have to be to water ?? I thought the sinks should be in sediments directly for root metabolism)
		    SEDNH4SAV, 	&  	!NH4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
							!(gN/m^2/day), SEDNH4SAV --> NH4T2TM1 (mg/N/m^3) decerase in sediments 
		    SEDNO3SAV, 	&  	!NO3 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
							!(gN/m^2/day), SEDNO3SAV --> NO3T2TM1 (mg/N/m^3) decerase in sediments 
		    SEDPO4SAV, 	&   !PO4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
		    
			!misc
		    DGRDAYS,   	&  	!Degree day calculated from temperature of water

			!
			!time average values 
			!
			
			ALEAF,   	& 	!Average leaf biomass (gC/m^2)
			AROOT,   	& 	!Average root biomass (gC/m^2)
			ASTEM,   	& 	!Average stem biomass (gC/m^2)       
			ATUBER,   	& 	!Average tuber biomass (gC/m^2)       
			AEP,		&   !Average epiphyte biomass (gC-epi/gC-leaf)

			APLEAF,   	&   !Averaged leaf growth rate (1/day)
			ABMLEAF,   	&   !Averaged basal metabolism rate of leaf (1/day)
			ABMTUBER,   &   !Average basal metabolism rate of tuber (1/day)

			APEP,    	&   !Average growth rate  for epiphytes (1/day)
			ABMEP,   	&   !Average basal metabolism rate for epiphytes (1/day)
			APREP,		&	!Average predation rate on epiphytes (1/day) = APREP*EP

			ASLSH,		&	!Average epiphute sloughing rate (1/day)
				
	       ANLSAV,    	&   !Averaged N limitation on SAV
	       APLSAV,    	&   !Averaged P limitation on SAV
	       ANLEPI,    	&   !Averaged N limitation on epiphytes
	       APLEPI,    	&   !Averaged P limitation on epiphytes
	       AFNSED,    	&   !Averaged fraction of uptake of N from sediment
	       AFPSED,    	&   !Averaged fraction of uptake of P from sediment
			 AFHS,   	& 	!Averaged SAV growth limiation due to sulfide HS^-            
    
		!attenuation and light limitation
		 AEPATN,   		&   !Averaged attenuation by epiphyte
		AWATATN,   		&   !Averaged water column attenuation of light at canopy top
		  AFISH,   		&   !Averaged light limitation for LEAF (shoots)
		  AFIEP,        &	!Averaged light limitation for Epiphytes      

		!net production                                          
		ANPPSAV,   		&   !Averaged net leaf production per day (gC/m^2/day)
		ANPPEPI,		&	!Averaged net epiphyte production per day (gC-epi/gC-leaf/day)

		!time averaged flux to water column
	       ADOSAVW,  	&  	!average DO  flux due to SAV gO2/m^2/day
	      ADOCSAVW,  	&  	!average DOC flux due to SAV
	      APOCSAVW,  	&  	!average POC flux due to SAV
	      ANH4SAVW,  	&  	!average NH4 flux due to SAV
	      ANO3SAVW,  	&  	!average NO3 flux due to SAV
	      ADONSAVW,  	&  	!average DON flux due to SAV              
	      APONSAVW,  	&  	!average PON flux due to SAV
	      APO4SAVW,  	&  	!average PO4 flux due to SAV
	      ADOPSAVW,  	&  	!average DOP flux due to SAV
	      APOPSAVW,  	&  	!average POL flux due to SAV
	       ADOEPIW,  	&  	!DO  flux (positive into water column ) due to epiphytes (gO2/m2/day)
	      ANH4EPIW,  	&  	!NH4 flux (positive into water column) due to epiphytes (gN/m2/day)
	      ANO3EPIW,  	&  	!NO3 flux (positive into water column) due to epiphytes (gN/m2/day)
	      APO4EPIW,  	&  	!PO4 flux (positive into water column) due to epiphytes (gP/m2/day)

	      ADOCEPIW,  	&  	! DOC flux to water column due to epiphytes (gC/m^2/day)
	      APOCEPIW,  	&  	! POC flux to water column due to epiphytes (gC/m^2/day)
	      ADONEPIW,  	&  	! DON flux to water column due to epiphytes (gN/m^2/day)
	      APONEPIW,  	&  	! PON flux to water column due to epiphytes (gN/m^2/day)                  
	      ADOPEPIW,  	&  	! DOP flux to water column due to epiphytes (gP/m^2/day)
	      APOPEPIW,		&  	! POP flux to water column due to epiphytes (gP/m^2/day)

		  ASEDDOSAV,   	&   !DO  flux due to SAV,  positive leaving sediments
	      ASEDPOCSAV,  	&  	!POC flux due to SAV, positive into sediments
	      ASEDPONSAV,  	&  	!PON flux due to SAV, positive into sediments
	      ASEDPOPSAV,  	&  	!POP flux due to SAV, positive into sediments
	      ASEDNH4SAV,  	&  	!NH4 flux due to SAV, positive leaving sediemnts
	      ASEDPO4SAV!,  &	!PO4 flux due to SAV, positive leaving sediments

#   if defined (MULTIPROCESSOR)
		USE MOD_PAR, ONLY:   &	  

			 NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 

			NHN,  		&	!			!!NUMBER OF HALO NODES
			HN_LST!, 	&	!			!!LIST OF HALO NODES 

#   endif        
    IMPLICIT NONE
    SAVE
    INTEGER :: MAX, I,	&	!WLong: This MAX should be removed, as it is a built-in max() function for float types
				JCON		!Index for constituents
				

    REAL(SP) :: FACTOR,FDNX,TENX,WMSNX,PRECIPNX, &
            ANH4NX,ANO3NX,ALDONNX,ARDONNX,APO4NX,ALDOPNX,ARDOPNX, &
            SLDOCNX,  &  !LDOC flux to water column due to SAV  (gC/m^2/day)
            SRDOCNX,  &  !RDOC flux to water column due to SAV  (gC/m^2/day)
            SLPOCNX,  &  !LPOC flux to water column due to SAV  (gC/m^2/day)
            SRPOCNX,  &  !RPOC flux to water column due to SAV  (gC/m^2/day)
            SDONX        !DO flux to water column due to SAV    (gC/m^2/day)
    REAL(SP) :: S2LNX, KTNX,  NXTVD,S1LNX, S3LNX
!
!Wen Long moved BFLUXNX, BFLUXNX_GL, BFLUX into wqm_module.f module name MOD_WQM!!!

     DIMENSION S1LNX(NS1P,NCP), S2LNX(NS2P,NCP), S3LNX(NS3P,NCP)


    NXTVD = TMEND
	! FDNX is fixed because the light is scaled by the time of day and reflectance in mod_owq
  FDNX = 1 ! Kurt Glaesemann 13 April 2015
    TENX = 0 ! Kurt Glaesemann 13 April 2015
10  CONTINUE
    DO WHILE(JDAY >= NXMET)
      KT = KTNX
      FD = FDNX

      IT = ITNX


      TTSS = 86400.*(1-FD)/2.
          
      TE = MAX(TENX,0.0)
      WMS = WMSNX
      READ (MET,*,END=11)NXMET,	&!	JDAY (day)
							KTNX,	&!	Coef of heat exchange (Watt/m^2/degC)
							TENX,	&!	Equilibrium temperature (degC)
							ITNX,	&!  I0 (PAR) (E/m2/day)
							FDNX,	&!	Fraction of day
							WMSNX	 !  Wind speed (m/s)
							
			
      NXMET = (METPTR-1)*FILGTH+NXMET
    ENDDO 
    GO TO 12

!******* Open next data file

11  CONTINUE
    METPTR = METPTR+1
    IF (DIAGNOSTICS) WRITE (DIA,*) 'Opening meteorologic file ',    &
                                    METPTR,' at day ',JDAY
    CLOSE (MET)
    OPEN  (MET,FILE=METFN(METPTR),STATUS='OLD')
    READ  (MET,1000)


    READ  (MET,*)NXMET,KTNX,TENX,ITNX,FDNX,WMSNX

    NXMET = (METPTR-1)*FILGTH+NXMET
    IF (JDAY >= NXMET) GOTO 10
12  CONTINUE
    NXTVD = MIN(NXTVD,NXMET)



!******* Atmospheric Loads

    IF (ATMOS_LOADS) THEN
133   CONTINUE
      DO WHILE(JDAY >= NXATM)
        PRECIP = PRECIPNX/8640000.
        ATMNH4 = ANH4NX
        ATMNO3 = ANO3NX
        ATMLDON = ALDONNX
        ATMRDON = ARDONNX
        ATMPO4 = APO4NX
        ATMLDOP = ALDOPNX
        ATMRDOP = ARDOPNX
        READ (ATM,1010,END=134) NXATM,PRECIPNX,ANH4NX,ANO3NX,      &
              ALDONNX,ARDONNX,APO4NX,ALDOPNX,ARDOPNX
        NXATM = (ATMPTR-1)*FILGTH+NXATM
      ENDDO
      GO TO 135

!********* Open next data file

134   CONTINUE
      ATMPTR = ATMPTR+1
      IF (DIAGNOSTICS) WRITE (DIA,*) 'Opening atmospheric load',   &
                                     'ing file ',ATMPTR,' at ',    &
                                     'day ',JDAY
      CLOSE (ATM)
      OPEN  (ATM,FILE=ATMFN(ATMPTR),STATUS='OLD')
      READ  (ATM,1000)
      READ  (ATM,1010) NXATM,PRECIPNX,ANH4NX,ANO3NX,ALDONNX,       &
                       ARDONNX,APO4NX,ALDOPNX,ARDOPNX
      NXATM = (ATMPTR-1)*FILGTH+NXATM
      IF (JDAY >= NXATM) GOTO 133
135   CONTINUE
      NXTVD = MIN(NXTVD,NXATM)
    ENDIF
	

!******* Submerged Aquatic Vegetation
!Wen Long, this part should be moved to wqm_sav.F as part of MOD_SAV
!   
    IF (SAV_LOADS) THEN
136   continue
      DO WHILE(JDAY >= NXSAV)
        DO B=1,MLOC
!WL            FACTOR = SAVAREA(B)/SFA(B)            !Wen Long: SFA is not calculated at all!
               FACTOR = SAVAREA(B)/ART1(B)           !replaced SFA by ART1 by Wen Long
          LDOCSAVW(B) = SLDOCNX  * FACTOR       
          RDOCSAVW(B) = SRDOCNX  * FACTOR       
          LPOCSAVW(B) = SLPOCNX  * FACTOR 
          RPOCSAVW(B) = SRPOCNX  * FACTOR
            DOSAVW(B) = SDONX    * FACTOR
        ENDDO  
        READ (SVI,1010,END=137) NXSAV, SLDOCNX, SRDOCNX, SLPOCNX, SRPOCNX, SDONX
        NXSAV = (SAVPTR-1)*FILGTH+NXSAV
      ENDDO
      GO TO 138

!********* Open next data file
137   CONTINUE
      SAVPTR = SAVPTR+1
      IF (DIAGNOSTICS) WRITE (DIA,*) 'Opening aquatic vegetation',   &
                                     ' file ',SAVPTR,' at day ',     &
                                      JDAY
      CLOSE (SVI)
      OPEN  (SVI,FILE=SVIFN(SAVPTR),STATUS='OLD')
      READ  (SVI,1080)
      READ  (SVI,1026) (SAVAREA(B),B=1,NBB)
      READ  (SVI,1080)
      READ  (SVI,1010) NXSAV, SLDOCNX, SRDOCNX, SLPOCNX, SRPOCNX, SDONX
      NXSAV = (SAVPTR-1)*FILGTH+NXSAV
      IF(JDAY >= NXSAV) GOTO 136
138   CONTINUE
      NXTVD = MIN(NXTVD,NXSAV)
    ENDIF


!******* Input FORMAT's

1000 FORMAT(///)
1010 FORMAT(10F8.0,:/(:8X,9F8.0))
1020 FORMAT(8X,9F8.0,:/(:16X,8F8.0))
1025 FORMAT(16X,F8.0)
1026 FORMAT(//(8X,9F8.0))
1030 FORMAT(//(8X,9I8))
1040 FORMAT(/)
1041 FORMAT(//)            !advance two lines
1050 FORMAT(://8X,9F8.0)   !go down 2 lines and read
                           !but do not advance to next line when all data 
                           !items are processed
1051 FORMAT(://8X,9F8.4)   !go down 2 lines and write, but do not
                           !continue going down when writing is finished
1060 FORMAT(8X,9F8.0)
1070 FORMAT(//(:8X,6F8.0))
1080 FORMAT(/)
1082 FORMAT(//8X,2A8//)
1100 FORMAT(/)
1101 FORMAT(1F8.0)
1102 FORMAT(1F8.0)
1103 FORMAT(10F8.0)
1104 FORMAT(10F8.4)

   END SUBROUTINE TVDS
 
