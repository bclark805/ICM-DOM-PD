MODULE MOD_SED_DOM

USE MOD_SED_DOM_EXCHANGE, ONLY: DLTS_SHARE,  &
								 W2_SHARE,    &
								 SED_TEMP_SHARE,   &
								 ZSED_SHARE,     &
								 JWCDOC1,      &
								 JWCDOC2,      &
								 JWCDOC3,      &
								 JWNCDOC1,      &
								 JWNCDOC2,      &
								 JWNCDOC3,      &
								 JWCDON1,      &
								 JWCDON2,      &
								 JWCDON3,      &
								 JWNCDON1,      &
								 JWNCDON2,      &
								 JWNCDON3,      &
								 JWCDOP1,      &
								 JWCDOP2,      &
								 JWCDOP3,      &
								 JWNCDOP1,      &
								 JWNCDOP2,      &
								 JWNCDOP3,  &
								 PERCOLATION_FLAG, & S_SHARE, KL12_SHARE, HSED1_SHARE
								 
USE WC_DOM, ONLY: WC_CDOC1,WC_CDOC2,WC_CDOC3, &
					  WC_NCDOC1,WC_NCDOC2,WC_NCDOC3, &
					  WC_CDON1,WC_CDON2,WC_CDON3, &
					  WC_NCDON1,WC_NCDON2,WC_NCDON3, &
					  WC_CDOP1,WC_CDOP2,WC_CDOP3, &
					  WC_NCDOP1,WC_NCDOP2,WC_NCDOP3
					  
USE MOD_HYDROVARS, ONLY : Nstationnum_gl ,   &  ! station output number	
						  Nstation,&             ! number of stations
						  ELNC2, ELNC1  ,  &  ! water surface elevation at time t and time t -1 
						  ART1,D,DZ ! area of node based control volume, to calculate total mass fluxed at each node
				   
USE MOD_SIZES, ONLY : MGL,  &  !global node and element indices
					  NGL    

USE MOD_LIMS, ONLY :  NTLOC, &              !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
                      MTLOC, &             !!TOTAL OF LOCAL INTERNAL + HALO NODES
		              MLOC,  &			   
					  KBM1,	 &
					  MYID,  &
					  NPROCS

#if defined (MULTIPROCESSOR)	
					  
USE MOD_PAR, ONLY : NGID ,   &         ! globale node ID of local node
					HN_LST,  &            ! list of halo nodes
					NHN,     &           ! Number of halo nodes
					GATHER,   &
					NMAP
# endif
					
USE MOD_CONTROL, ONLY: MSR, SERIAL, PAR

USE MOD_PREC, ONLY : SP
					  
USE MOD_WQM, ONLY : SAV_CALC, C2, ShearStress, Rho1!,     &
				!    LDOC, RDOC,       &
				!	LDON, RDON,       &
				!	LDOP, RDOP
					
			   				 				   
IMPLICIT NONE
!
!This module will read in the kelp definition file from 
! mod_kelp to get the density of plants.  The density of plants will than be converted to
!a biomass that will be used to input carbon into our marsh sees.

!the carbon from the plants will be converted to colored DOM and non-colored DOM
!with sinks for bacterial remineralization, sedimentation and diagenesis.

!there will be a diffusive flux that transfers the DOM into the water column as the 
!marsh is wet

!Started B Clark April 4th, 2015

!SUBROUTINE SED_DOM_DEFINE_VARS


INTEGER :: I,B,K,  &
		   NMARSHTOTAL  !total number of nodes in the marsh with marsh plants, from marsh file
!********************** LOGICALS ***************************

LOGICAL  :: SED_DOM = .TRUE. 
LOGICAL  :: MARSH_SEDIMENTS = .TRUE.
LOGICAL  :: HYDRO_FLAG = .TRUE.
LOGICAL :: WETLAND_FLAG = .TRUE.

!LOGICAL  :: SAV_CALC 

!********************* END LOGICALS **********************

!REAL                   ::  HSED,                  &    !depth of marsh sediments, m
!			               ZOLW			     !depth of the overlying water column, m

! **** temp outside vars for compiling seperate from mod_sed

!REAL :: DT_GL

!INTEGER		       :: MTLOC=8903
REAL(SP)		       :: DT
!INTEGER	  :: DTIN !the time step the model is on 
!INTEGER  		   ::   DTINT ! time step in this module
!!REAL		       :: DLTS
REAL(SP),ALLOCATABLE,DIMENSION(:)  :: ZSED ! sediment depth for this module

REAL(SP),ALLOCATABLE,DIMENSION(:)  ::SED_TEMP ! sediment temperature for this module



!*********************** PARAMETERS **************************
REAL(SP)      :: W2_Z   ! settling velocity divided by depth, used for integration

REAL(SP)      :: DOM_diff != 0.0104   ! THE DOM DIFFUSION RATE, CAN BE SPECIFIED OR TAKEN FROM SEDIMENT MODULE, 	m^2/day
							    !! DOM diffusion coefficient 1.2x10^-5 cm/s from Hesslein 1980
REAL(SP)      :: DDOM_TCON != 1.00  ! DOM diffusion rate temperature control
REAL(SP)      :: DDOM

REAL(SP), ALLOCATABLE :: PSTY(:), &  !sediment porosity unitless <1 if sediment is more porous than basal DDOM > 1 if sediment is less porous than DDOM
							DNSTY(:),&   !  and density, gC cm^-1
							SED_ADV(:)  !advection of water via percolation (m/day), f
REAL(SP)      :: SEDPERC_COEFF =  1.143,  &  ! SEdiment percolation coefficnet, multiplied by the sediment density to yield the sqrt of sediment percolation velocity, added to the mixing velocity
                  SEDPERC_INT = -0.478!,&     ! intercept of linear relationship between density and velocity
				  
				  
REAL(SP) :: TEMP_ABSC = 0.0, TEMP_ABSN =0.0		! add in a temporary irreversible absorption sink for testing purposes


				   ! above is from Bradley and Morris 1990 MEPS
				  
				  !SED_ADV = SEDPERC_COEFF*DNSTY(I)+SEDPERC_INT)/100*24 (convert from cm/hr to m/day) ! sediment vertical advection of porewater
				  
				!******FALLING TIDE*************
				  !CDOC1 is LOST out of bottom via w2 - SED_ADV /ZSED*CDOC1
				  !if wet gain CDOC1 from water vi +SED_ADV/Z(KBM1)*WC_DOC1
				  
				!******RISING TIDE*************
				  !CDOC1 is GAINED out of bottom via w2 + SED_ADV /ZSED*CDOC1
				  !if wet LOSE CDOC1 to water vi -SED_ADV/Z(KBM1)*CDOC1                                                        
				  
				  
                  
REAL(SP)	  :: DD,DD_Z  ! THE MIXING VELOCITY (m/day) and depth integrated diffusion rate (1/day)

REAL(SP)	  :: ALPHA1, ALPHA2    ! the reaction coefficients for the sediment DOM solver for layer 1 and layer 2
									!, 1/day


REAL(SP)                   ::       NCDOCFRAC ,             &    !fraction of MDOM to NCDOM  These should add up to 1 for conservation
							         CDOCFRAC                     !fraction of MDOM to CDOM

REAL(SP)	   ::   W2M 			        ! loss velocity of DOC (m/day) in this module

REAL(SP)                  ::   M_NtoC   ,	        & !the marsh N:C ratio (gN/gC) 
							M_PtoC   		  !the marsh P:C rati  (gP/gC)


REAL(SP)		   ::           ANCM     ,              & ! conversion factor to local N:C 
							APCM	 		  ! conversion factor to local P:C

REAL(SP)      	           ::   MFDOC1    , 		& ! the marsh G1 DOC fractionation coefficient
							    MFDOC2	  ,		& ! the marsh G2 DOC fractionation coefficient
							    MFDOC3			  ! the marsh G3 DOC fractionation coefficient
								
REAL(SP)      	           ::   WFDOM1    , 		& ! THE FRACTION OF LDOM FROM WQM THAT GOES TO G2 IN SEDIMENTS
							    WFDOM2   ,          & ! THE FRACTION OF RDOM FROM WQM THAT GOES TO G2 IN SEDIMENTS
								WFDOM3                !THE FRACTION OF G2 CARBON THAT IS LABILE FROM SEDS ------> WATER COLUMN	
								
!REAL(SP) 	               :: WC_FRAC = 0.25         ! THE FRACTION OF WATER COLUMN DOM THAT IS COLORED, this is temporary until the colored
													! non-colored carbon kinetics are coded into the WQM

REAL(SP)                  ::   KDOC1    ,              & ! the marsh doc1 decay rate  (1/day)
							KDOC2    ,		& ! the marsh doc2 decay rate  (1/day)
							KDOC3	  ,		& ! the marsh doc3 decay rate  (1/day)	
							KDON1    ,		& ! the marsh DON1 decay rate  (1/day) 
							KDON2    , 		& ! the marsh DON2 decay rate  (1/day)
							KDON3    ,		& ! the marsh DON3 decay rate  (1/day)	
							KDOP1    ,		& ! the marsh DOP1 decay rate  (1/day)
							KDOP2    ,              & ! the marsh DOP2 decay rate  (1/day)
							KDOP3			 ! the marsh DOP3 decay rate  (1/day)

			
REAL(SP)       	   ::   thetaDOC1,		& ! temperature theta for diagenesis reaction DOC1
					thetaDOC2,		& ! temperature theta for diagenesis reaction DOC2
					thetaDOC3,		& ! temperature theta for diagenesis reaction DOC3
					thetaDON1,		& ! temperature theta for diagenesis reaction DON1
					thetaDON2,		& ! temperature theta for diagenesis reaction DON2
					thetaDON3,		& ! temperature theta for diagenesis reaction DON3
					thetaDOP1,		& ! temperature theta for diagenesis reaction DOP1
					thetaDOP2, 		& ! temperature theta for diagenesis reaction DOP2  
					thetaDOP3		  ! temperature theta for diagenesis reaction DOP3

					
					
real(sp)		::  POCCF1,POCCF2,POCCF3   !fraction of colored hydrolyzed POC	


!********************** END PARAMETERS **************************

!************************* IO VARS ******************************

character(len=120)::FSEDDOM_GRID,FSEDDOM_CON !input file names
character(len=120)::SEDDOM_GRIDIN,SEDDOM_CONIN 
character(len=120)::SEDDOM_OUT,FSEDDOM_OUT !output file
character(len=120)::WDOM_OUT,FWDOM_OUT !output file FOR WATER DOM
character(len=120)::JWDOM_OUT,FJWDOM_OUT ! output for fluxes

!temporary variables for reading the file in
INTEGER, ALLOCATABLE 			:: NMNTEMP(:)
REAL(SP), ALLOCATABLE :: MASSTEMP(:),  &
					   FRACTEMP(:),  &
					   EX_TEMP(:), &!
					   PSTY_TEMP(:), DNSTY_TEMP(:),  &! 
					   DOM_INITIAL_TEMP(:),DOM_INITIAL(:)  ! B Clark,
					   !temporary solution to see if a spatially variant settling rate helps with the NH4 fluxes out of sediments

INTEGER, DIMENSION(:),ALLOCATABLE  ::   NMARSHNODE           ! node number from input file



REAL(SP), DIMENSION(:), ALLOCATABLE  ::	MARSH_FRAC,          &   !fraction of the TCE that is marsh plants  
				MARSH_MASS,          &   !plant biomass in each marsh cell, (specified)(g C/m^2)
				P_EX_RATE,  &                !dom exudation rate (specified)(1/day)
				Plant_NH4_UPTAKE   ! Uptake flux of NH4 due to wetland plant growth				
				
REAL(SP) :: EXRATE_TEMPCON ! temperature control of marsh plant exudation
REAL(SP) :: EXRATE ! exudation rate after temperature affects of marsh plants DOM

REAL(SP) :: VMAX_MARSH_NH4 = 1.344 ! Maximum marsh NH4 uptake rate mg N g root^-1 day^-1 before temperature control; Morris 1990
REAL(SP) :: VMAX_MARSH_NH4_out = 0.0! marsh Nh4 uptake rate after temperature adjustement

REAL(SP) :: KMNH4_MARSH = 70.0  ! half saturation rate of Michaelis-Menton marsh NH4 uptake kinetics (mg N m^-3); Morris 1990  at 20 ppt salinity

REAL(SP) :: TAU_crit, Resusp_mass, RoughHeight ! Critical shear stress, Pascals, added by B Clark June 2018



!************************ END IO VARS ***************************

!************************* CALC VARS ****************************

REAL(SP) 	                         ::	JMDOM      	             !marsh plant dom exudation (calculated) (g C/m^2/day)’	
				
REAL(SP),   DIMENSION(:),  ALLOCATABLE  ::	JNCMDOC,	        &   !Non-colored marsh DOM flux (gC/m^2/day)
											JCMDOC   	           !colored marsh DOM flux	(gC/m^2/day)

!Anaerobic Layer Carbon
REAL(SP)     		   ::   CDOC1TM1 ,               & ! colored [DOC1] at time t-dt  ( gC/m^3)   
                                CDOC2TM1 , 		& ! colored [DOC2] at time t-dt  (gC/m^3)
                                CDOC3TM1 , 		& ! colored [DOC3] at time t-dt
                                NCDOC1TM1 ,		& ! non colored [DOC1] at time t-dt
                                NCDOC2TM1 ,		& ! non colored [DOC2] at time t-dt
                                NCDOC3TM1	          ! non colored [DOC3] at time t-dt		

!Aerobic Layer Carbon
REAL(SP)     		   ::   CDOC11TM1 ,               & ! colored [DOC1]  at time t-dt  ( gC/m^3)   
                                CDOC21TM1 , 		& ! colored [DOC2] at time t-dt  (gC/m^3)
                                CDOC31TM1 , 		& ! colored [DOC3] at time t-dt
                                NCDOC11TM1 ,		& ! non colored [DOC1] at time t-dt
                                NCDOC21TM1 ,		& ! non colored [DOC2] at time t-dt
                                NCDOC31TM1	          ! non colored [DOC3] at time t-dt		
                                        
!Anaerobic layer Nitrogen                               
REAL(SP)                   ::   CDON1TM1 ,		& ! colored [DON1] at time t-dt   (g N/m^3)
                                CDON2TM1 ,		& ! colored [DON2] at time t-dt
                                CDON3TM1 ,		& ! colored [DON3] at time t-dt
                                NCDON1TM1 ,		& ! non colored [DON1] at time t-dt
                                NCDON2TM1 ,		& ! non colored [DON2] at time t-dt
                                NCDON3TM1 		  ! non colored [DON3] at time t-dt
								
!Aerobic layer Nitrogen                               
REAL(SP)                   ::   CDON11TM1 ,		& ! colored [DON1] at time t-dt   (g N/m^3)
                                CDON21TM1 ,		& ! colored [DON2] at time t-dt
                                CDON31TM1 ,		& ! colored [DON3] at time t-dt
                                NCDON11TM1 ,		& ! non colored [DON1] at time t-dt
                                NCDON21TM1 ,		& ! non colored [DON2] at time t-dt
                                NCDON31TM1 		  ! non colored [DON3] at time t-dt
                                        
!Aerobic layer P	                               
REAL(SP)                   ::   CDOP1TM1 , 		& ! colored [DOP1] at time t-dt (g P/m^3)
                                CDOP2TM1 ,		& ! colored [DOP2] at time t-dt 
                                CDOP3TM1 ,		& ! colored [DOP3] at time t-dt
                                NCDOP1TM1 ,		& ! non colored [DOP1] at time t-dt
                                NCDOP2TM1 ,		& ! non colored [DOP2] at time t-dt
                                NCDOP3TM1 		  ! non colored [DOP3] at time t-dt
!Anaerobic layer P								
REAL(SP)                   ::   CDOP11TM1 , 		& ! colored [DOP1] at time t-dt (g P/m^3)
                                CDOP21TM1 ,		& ! colored [DOP2] at time t-dt 
                                CDOP31TM1 ,		& ! colored [DOP3] at time t-dt
                                NCDOP11TM1 ,		& ! non colored [DOP1] at time t-dt
                                NCDOP21TM1 ,		& ! non colored [DOP2] at time t-dt
                                NCDOP31TM1 		  ! non colored [DOP3] at time t-dt

REAL(SP)                   ::   XKDOC1,			& !temperature adjusted decay rate (1/day)
				XKDOC2,			&
			        XKDOC3,			&
				XKDON1,			&
				XKDON2,			&
				XKDON3,			&
				XKDOP1,			&
				XKDOP2,			&
				XKDOP3

!Aerobic layer inputs from overlying water column, plants etc
REAL(SP)		::   JinCDOC11,	        & ! New CDOC1 input  (g C/m^3)
				JinCDOC21,		& ! New	CDOC2 input  (g C/m^3)
				JinCDOC31,		& ! New CDOC3 input  (g C/m^3)
				JinNCDOC11,		& ! NEW NCDOC1 input (g C/m^3)
				JinNCDOC21,		& ! New NCDOC2 input (g C/m^3)
				JinNCDOC31		  ! New NCDOC3 input (g C/m^3)

REAL(SP)		   ::	JinCDON11,		& ! New CDON1 input  (g N/m^3)
				JinCDON21,		& ! New CDON2 input  (g N/m^3) 
				JinCDON31,		& ! New CDON3 input  (g N/m^3)
				JinNCDON11,		& ! NEW NCDON1 input (g N/m^3)
				JinNCDON21,		& ! New NCDON2 input (g N/m^3)
				JinNCDON31		  ! New NCDON3 input (g N/m^3)

REAL(SP)		   ::	JinCDOP11,		& ! New CDOP1 input  (g N/m^3)
				JinCDOP21,		& ! New CDOP2 input  (g N/m^3) 
				JinCDOP31,		& ! New CDOP3 input  (g N/m^3)
				JinNCDOP11,		& ! NEW NCDOP1 input (g N/m^3)
				JinNCDOP21,		& ! New NCDOP2 input (g N/m^3)
				JinNCDOP31		  ! New NCDOP3 input (g N/m^3)				
				
				
				
!! Anaerobic Layer INputs from plants, etc				
REAL(SP)		   ::   JinCDOC1,	        & ! New CDOC1 input  (g C/m^3)
				JinCDOC2,		& ! New	CDOC2 input  (g C/m^3)
				JinCDOC3,		& ! New CDOC3 input  (g C/m^3)
				JinNCDOC1,		& ! NEW NCDOC1 input (g C/m^3)
				JinNCDOC2,		& ! New NCDOC2 input (g C/m^3)
				JinNCDOC3		  ! New NCDOC3 input (g C/m^3)

REAL(SP)		   ::	JinCDON1,		& ! New CDON1 input  (g N/m^3)
				JinCDON2,		& ! New CDON2 input  (g N/m^3) 
				JinCDON3,		& ! New CDON3 input  (g N/m^3)
				JinNCDON1,		& ! NEW NCDON1 input (g N/m^3)
				JinNCDON2,		& ! New NCDON2 input (g N/m^3)
				JinNCDON3		  ! New NCDON3 input (g N/m^3)

REAL(SP)		   ::	JinCDOP1,		& ! New CDOP1 input  (g N/m^3)
				JinCDOP2,		& ! New CDOP2 input  (g N/m^3) 
				JinCDOP3,		& ! New CDOP3 input  (g N/m^3)
				JinNCDOP1,		& ! NEW NCDOP1 input (g N/m^3)
				JinNCDOP2,		& ! New NCDOP2 input (g N/m^3)
				JinNCDOP3		  ! New NCDOP3 input (g N/m^3)
				
				
				
				

REAL(SP)				:: C1TEMP,    &  ! SEDIMENT LAYER TEMPORARY CONCENTRATION FROM DIFFUSION
						   C2TEMP	 ! WATER COLUMN CONCENTRATION FROM DIFFUSION, UNUSED,
				
!Anaerbioc Layer Carbon
REAL(SP),DIMENSION(:), ALLOCATABLE  ::   CDOC1     ,                & !the Colored G1 DOC  gC/m^3
										 CDOC2     ,		& !the colored G2 DOC  gC/m^3
										 CDOC3	 	          !the colored G3 DOC  gC/m^3

REAL(SP),DIMENSION(:), ALLOCATABLE  ::   NCDOC1    ,		& !the non-colored G1 DOC gC/m^3
										 NCDOC2    ,		& !the non-colored G2 DOC gC/m^3
										 NCDOC3	                  !the non-colored G3 DOC gC/m^3
! Aerobic Layer carbon										 
REAL(SP),DIMENSION(:), ALLOCATABLE  ::   CDOC11     ,                & !the Colored G1 DOC  gC/m^3
										 CDOC21     ,		& !the colored G2 DOC  gC/m^3
										 CDOC31	 	          !the colored G3 DOC  gC/m^3

REAL(SP),DIMENSION(:), ALLOCATABLE  ::   NCDOC11    ,		& !the non-colored G1 DOC gC/m^3
										 NCDOC21   ,		& !the non-colored G2 DOC gC/m^3
										 NCDOC31	                  !the non-colored G3 DOC gC/m^3										 
! Anaerobic layer nitrogen
REAL(SP), DIMENSION(:),  ALLOCATABLE  ::   CDON1     ,              & !the Colored G1 DON  gN/m^3
										   CDON2     ,	         	& !the colored G2 DON  gN/m^3
										   CDON3 	 	                  !the colored G3 DON  gN/m^3

REAL(SP),DIMENSION(:), ALLOCATABLE  ::   NCDON1    ,		& !the non-colored G1 DON gN/m^3
										 NCDON2    ,	        	& !the non-colored G2 DON gN/m^3
										 NCDON3	                          !the non-colored G3 DON gN/m^3

!Aerobic layer nitrogen
REAL(SP), DIMENSION(:),  ALLOCATABLE  ::   CDON11     ,              & !the Colored G1 DON  gN/m^3
										   CDON21     ,	         	& !the colored G2 DON  gN/m^3
										   CDON31 	 	                  !the colored G3 DON  gN/m^3

REAL(SP),DIMENSION(:), ALLOCATABLE  ::   NCDON11    ,		& !the non-colored G1 DON gN/m^3
										 NCDON21   ,	        	& !the non-colored G2 DON gN/m^3
										 NCDON31	                          !the non-colored G3 DON gN/m^3
!aerobic layer P										 
REAL(SP),DIMENSION(:), ALLOCATABLE  ::   CDOP11     ,                & !the Colored G1 DOP  gP/m^3
										 CDOP21          ,		& !the colored G2 DOP  gP/m^3
										 CDOP31	 	                  !the colored G3 DOP  gP/m^3

REAL(SP),DIMENSION(:), ALLOCATABLE  ::   NCDOP11    ,		& !the non-colored G1 DOP gP/m^3
										 NCDOP21       ,		& !the non-colored G2 DOP gP/m^3
										 NCDOP31	  	                  !the non-colored G3 DOP gP/m^3
										 
!anaerobic layer p										 
REAL(SP),DIMENSION(:), ALLOCATABLE  ::   CDOP1     ,                & !the Colored G1 DOP  gP/m^3
										 CDOP2          ,		& !the colored G2 DOP  gP/m^3
										 CDOP3	 	                  !the colored G3 DOP  gP/m^3

REAL(SP),DIMENSION(:), ALLOCATABLE  ::   NCDOP1    ,		& !the non-colored G1 DOP gP/m^3
										 NCDOP2       ,		& !the non-colored G2 DOP gP/m^3
										 NCDOP3	  	                  !the non-colored G3 DOP gP/m^3


REAL(SP), ALLOCATABLE  ::          VAR_OUT_GL(:,:), &              ! the output SEDDOM variable array, size is (MGL,18)
                                   !LDOC_OUT_GL(:,:),&         ! water dom output size is (MGL,KBM1)
                                   !RDOC_OUT_GL(:,:),&
                                   !LDON_OUT_GL(:,:),&         ! water dom output
                                   !RDON_OUT_GL(:,:),&
                                   !LDOP_OUT_GL(:,:),&         ! water dom output
                                   !RDOP_OUT_GL(:,:), &!
									JDOM_OUT(:,:)  ,&           ! Flux between seds and water size is (MGL,6)
								   JDOM_OUT_GL(:,:)  ,&           ! Flux between seds and water size is (MGL,6)
								   WCDOM_GL(:,:,:)!,  &
								   
								   
REAL(SP), ALLOCATABLE :: TIDE_STAGE(:)      ! Difference in height at time t and t-1 to see if the tide is ebb or flow, calculated	


!************************* END CALC VARS ************************
										 
							 
!========================== WATER COLUMN ORGANICS fluxes =============================	
										 
! temp variables used to calculate the diffusive flux to the water column from the sediments
									

!									
!REAL 							::   WCDOC1     ,       & !the Colored G1 DOC  gC/m^3
!									 WCDOC2     ,		& !the colored G2 DOC  gC/m^3
!									 WCDOC3	 	          !the colored G3 DOC  gC/m^3

!REAL						    ::   WNCDOC1    ,		& !the non-colored G1 DOC gC/m^3
!									 WNCDOC2    ,		& !the non-colored G2 DOC gC/m^3
!									 WNCDOC3	                  !the non-colored G3 DOC gC/m^3

!REAL						    ::   WCDON1     ,              & !the Colored G1 DON  gN/m^3
!									   WCDON2     ,	         	& !the colored G2 DON  gN/m^3
!									   WCDON3 	 	                  !the colored G3 DON  gN/m^3

!REAL	                        ::   WNCDON1    ,		& !the non-colored G1 DON gN/m^3
!									 WNCDON2    ,	        	& !the non-colored G2 DON gN/m^3
!									 WNCDON3	                          !the non-colored G3 DON gN/m^3

!REAL						    ::   WCDOP1     ,                & !the Colored G1 DOP  gP/m^3
!								     WCDOP2          ,		& !the colored G2 DOP  gP/m^3
!									 WCDOP3	 	                  !the colored G3 DOP  gP/m^3

!REAL						    ::   WNCDOP1    ,		& !the non-colored G1 DOP gP/m^3
!									 WNCDOP2       ,		& !the non-colored G2 DOP gP/m^3
!									 WNCDOP3	  	                  !the non-colored G3 DOP gP/m^3
!									 
									 									

!*************************** FLUX VARS ***************************

REAL(SP)		           ::   JDIAGDOP1    ,		& ! the flux of G1 DOP lost to diagenesis (g P/m^2/day)
								JDIAGDOP2    ,		& ! the flux of G2 DOP lost to diagenesis (g P/m^2/day)
								JDIAGDOP3		  ! the flux of G3 DOP lost to diagenesis (g P/m^2/day)

	
REAL(SP)                   ::   JDIAGDOC1     ,         & ! the flux of G1 DOC lost to diagenisis (gC/m^2/day)
								JDIAGDOC2     ,         & ! the flux of G2 DOC lost to diagenisis (gC/m^2/day)
								JDIAGDOC3       	  ! the flux of G3 DOC lost to diagensis  (gC/m^2/day)	


REAL(SP)		           ::   JDIAGDON1   ,		& ! the flux of G1 DON lost to diagenesis (g N/m^2/day)
								JDIAGDON2   ,		& ! the flux of G2 DON lost to diagenesis (g N/m^2/day)
								JDIAGDON3		  ! the flux of G3 DON lost to diagenesis (g N/m^2/day)


REAL(SP),DIMENSION(:), ALLOCATABLE  ::   JDOCX_OUT1   , 		& ! the flux of DOC that is passed to diagenesis in mod_sed (g C/m^2/day)
										  JDOCX_OUT2   , 		& ! the flux of DOC that is passed to diagenesis in mod_sed (g C/m^2/day) 
										JDONX_OUT1   ,		& ! the flux of DON that is passed to diagenesis in mod_sed (g N/m^2/day) in layer 1
										 JDONX_OUT2   ,		& ! the flux of DON that is passed to diagenesis in mod_sed (g N/m^2/day) in layer 2
										 JDOPX_OUT1	 , &     	  ! the flux of DOP that is passed to diagenesis in mod_sed (g P/m^2/day)
										JDOPX_OUT2
!************************ END FLUX VARS **************************
!==========*=================*===============*=============*=====*


CONTAINS

SUBROUTINE SED_DOM_ALLOCATE_VARS
IMPLICIT NONE
!temporary
ALLOCATE(ZSED(MTLOC))
ALLOCATE(SED_TEMP(MTLOC))
!our I/O vars

!our run vars

  ALLOCATE(Plant_NH4_UPTAKE(MTLOC))

 ALLOCATE(JNCMDOC(MTLOC))
 ALLOCATE(JCMDOC(MTLOC))

 ALLOCATE(CDOC1(MTLOC))
 ALLOCATE(CDOC2(MTLOC))
 ALLOCATE(CDOC3(MTLOC))

 ALLOCATE(NCDOC1(MTLOC))
 ALLOCATE(NCDOC2(MTLOC))
 ALLOCATE(NCDOC3(MTLOC))
 
 ALLOCATE(CDOC11(MTLOC))
 ALLOCATE(CDOC21(MTLOC))
 ALLOCATE(CDOC31(MTLOC))

 ALLOCATE(NCDOC11(MTLOC))
 ALLOCATE(NCDOC21(MTLOC))
 ALLOCATE(NCDOC31(MTLOC))
 
 ALLOCATE(CDON1(MTLOC))
 ALLOCATE(CDON2(MTLOC))
 ALLOCATE(CDON3(MTLOC))

 ALLOCATE(NCDON1(MTLOC))
 ALLOCATE(NCDON2(MTLOC))
 ALLOCATE(NCDON3(MTLOC))
 
 ALLOCATE(CDON11(MTLOC))
 ALLOCATE(CDON21(MTLOC))
 ALLOCATE(CDON31(MTLOC))

 ALLOCATE(NCDON11(MTLOC))
 ALLOCATE(NCDON21(MTLOC))
 ALLOCATE(NCDON31(MTLOC))
 
  ALLOCATE(CDOP11(MTLOC))
 ALLOCATE(CDOP21(MTLOC))
 ALLOCATE(CDOP31(MTLOC))

 ALLOCATE(NCDOP11(MTLOC))
 ALLOCATE(NCDOP21(MTLOC))
 ALLOCATE(NCDOP31(MTLOC))
 
 
 ALLOCATE(CDOP1(MTLOC))
 ALLOCATE(CDOP2(MTLOC))
 ALLOCATE(CDOP3(MTLOC))

 ALLOCATE(NCDOP1(MTLOC))
 ALLOCATE(NCDOP2(MTLOC))
 ALLOCATE(NCDOP3(MTLOC))

 ALLOCATE(JDOCX_OUT1(MTLOC))
  ALLOCATE(JDOCX_OUT2(MTLOC))
 ALLOCATE(JDONX_OUT1(MTLOC))
 ALLOCATE(JDONX_OUT2(MTLOC))
 ALLOCATE(JDOPX_OUT1(MTLOC))
  ALLOCATE(JDOPX_OUT2(MTLOC))
 
 ALLOCATE(VAR_OUT_GL(MGL,1:18)) ! allocate to the global array size and the 18 variables
 ALLOCATE(WCDOM_GL(MGL,KBM1,18))
 
 ALLOCATE(JDOM_OUT(MTLOC,1:18))     ! For Flux output
 ALLOCATE(JDOM_OUT_GL(MGL,1:18))     ! For Flux output
 
 ALLOCATE(TIDE_STAGE(MTLOC))

END SUBROUTINE SED_DOM_ALLOCATE_VARS

!====*==========*================*================*================*

SUBROUTINE SED_DOM_DEALLOC

!deallocate all the vars!!
!***** temporary for development

IF(ALLOCATED(ZSED)) DEALLOCATE(ZSED)
IF(ALLOCATED(SED_TEMP)) DEALLOCATE(SED_TEMP)

IF(ALLOCATED(Plant_NH4_uptake)) Deallocate(Plant_NH4_Uptake)

! I/O vars

IF(ALLOCATED(NMARSHNODE)) DEALLOCATE(NMARSHNODE)
IF(ALLOCATED(MARSH_FRAC)) DEALLOCATE(MARSH_FRAC)
IF(ALLOCATED(MARSH_MASS)) DEALLOCATE(MARSH_MASS)
IF(ALLOCATED(P_EX_RATE))  DEALLOCATE(P_EX_RATE)
IF(ALLOCATED(DNSTY))  DEALLOCATE(DNSTY)
IF(ALLOCATED(PSTY))  DEALLOCATE(PSTY)
IF(ALLOCATED(DOM_INITIAL)) DEALLOCATE(DOM_INITIAL)

! run vars

IF(ALLOCATED(JNCMDOC)) DEALLOCATE(JNCMDOC)
IF(ALLOCATED(JCMDOC))  DEALLOCATE(JCMDOC)

IF(ALLOCATED(CDOC1))   DEALLOCATE(CDOC1)
IF(ALLOCATED(CDOC2))   DEALLOCATE(CDOC2)
IF(ALLOCATED(CDOC3))   DEALLOCATE(CDOC3)

IF(ALLOCATED(NCDOC1))  DEALLOCATE(NCDOC1)
IF(ALLOCATED(NCDOC2))  DEALLOCATE(NCDOC2)
IF(ALLOCATED(NCDOC3))  DEALLOCATE(NCDOC3)

IF(ALLOCATED(CDOC11))   DEALLOCATE(CDOC11)
IF(ALLOCATED(CDOC21))   DEALLOCATE(CDOC21)
IF(ALLOCATED(CDOC31))   DEALLOCATE(CDOC31)

IF(ALLOCATED(NCDOC11))  DEALLOCATE(NCDOC11)
IF(ALLOCATED(NCDOC21))  DEALLOCATE(NCDOC21)
IF(ALLOCATED(NCDOC31))  DEALLOCATE(NCDOC31)

IF(ALLOCATED(CDON1))   DEALLOCATE(CDON1)
IF(ALLOCATED(CDON2))   DEALLOCATE(CDON2)
IF(ALLOCATED(CDON3))   DEALLOCATE(CDON3)

IF(ALLOCATED(NCDON1))  DEALLOCATE(NCDON1)
IF(ALLOCATED(NCDON2))  DEALLOCATE(NCDON2)
IF(ALLOCATED(NCDON3))  DEALLOCATE(NCDON3)

IF(ALLOCATED(CDON11))   DEALLOCATE(CDON11)
IF(ALLOCATED(CDON21))   DEALLOCATE(CDON21)
IF(ALLOCATED(CDON31))   DEALLOCATE(CDON31)

IF(ALLOCATED(NCDON11))  DEALLOCATE(NCDON11)
IF(ALLOCATED(NCDON21))  DEALLOCATE(NCDON21)
IF(ALLOCATED(NCDON31))  DEALLOCATE(NCDON31)

IF(ALLOCATED(CDOP1))   DEALLOCATE(CDOP1)
IF(ALLOCATED(CDOP2))   DEALLOCATE(CDOP2)
IF(ALLOCATED(CDOP3))   DEALLOCATE(CDOP3)

IF(ALLOCATED(NCDOP1))  DEALLOCATE(NCDOP1)
IF(ALLOCATED(NCDOP2))  DEALLOCATE(NCDOP2)
IF(ALLOCATED(NCDOP3))  DEALLOCATE(NCDOP3)

IF(ALLOCATED(CDOP11))   DEALLOCATE(CDOP11)
IF(ALLOCATED(CDOP21))   DEALLOCATE(CDOP21)
IF(ALLOCATED(CDOP31))   DEALLOCATE(CDOP31)

IF(ALLOCATED(NCDOP11))  DEALLOCATE(NCDOP11)
IF(ALLOCATED(NCDOP21))  DEALLOCATE(NCDOP21)
IF(ALLOCATED(NCDOP31))  DEALLOCATE(NCDOP31)

IF(ALLOCATED(VAR_OUT_GL)) DEALLOCATE(VAR_OUT_GL)
IF(ALLOCATED(WCDOM_GL)) DEALLOCATE(WCDOM_GL)

IF(ALLOCATED(JDOM_OUT)) DEALLOCATE(JDOM_OUT)
IF(ALLOCATED(JDOM_OUT_GL)) DEALLOCATE(JDOM_OUT_GL)


IF(ALLOCATED(TIDE_STAGE)) DEALLOCATE(TIDE_STAGE)
IF(ALLOCATED(SED_ADV)) DEALLOCATE(SED_ADV)



! ========================= Flux to Sediments ============================

IF(ALLOCATED(JDOCX_OUT1))  DEALLOCATE(JDOCX_OUT1)
IF(ALLOCATED(JDOCX_OUT2))  DEALLOCATE(JDOCX_OUT2)
IF(ALLOCATED(JDONX_OUT1))  DEALLOCATE(JDONX_OUT1)
IF(ALLOCATED(JDONX_OUT2))  DEALLOCATE(JDONX_OUT2)
IF(ALLOCATED(JDOPX_OUT1))  DEALLOCATE(JDOPX_OUT1)
IF(ALLOCATED(JDOPX_OUT1))  DEALLOCATE(JDOPX_OUT2)

END SUBROUTINE SED_DOM_DEALLOC

!=========*=================*=================*=================*

!DOM from roots (mg/L) in sediments (totalmasss) = DOC + DOP +DON + others 

!DOC is mgC//L from roots (roots also is in carbon units) 
!DON is DOC * a N/C ratio 
!DOP is DOC * a P/C ratio

!DOC  ==> NCDOC + CDOC =(NCDOC1 + NCDOC2 +NCDOC3)

!all plants(roots, leaf, shoot) are in carbon units
!DOM in carbon units first(DOC), then multiply by alpha_NC to get DON, alpha_PC to get DOP

!======================*============================*============================*


SUBROUTINE SED_DOM_INITIALIZE
IMPLICIT NONE

!MARSH_SEDIMENTS = .TRUE.
SAV_CALC = .FALSE.

!DT = 1
!ZSED = 0.1
!CTEMP = 20.0

!ZSED = 0.1      
!ZOLW = 2	  
           
           
 !NCDOCFRAC  = 0.4     ! NCDOCFRAC+CDOCFRAC = 1  
 !CDOCFRAC   = 0.6
 !                    
 !M_NtoC     = 0.083333 ! 1:12 N:C by weight
 !M_PtoC     = 0.0125! 1:80 P:C by weight
           
           
 ANCM     = 0.0
 APCM	  = 0.0
           
           
           
 !MFDOC1 = 0.65 ! marsh DOC G1 frac from Di Toro et al pg 253 
 !MFDOC2	= 0.20  ! marsh DOC G2 frac 
 !MFDOC3	= 0.15  ! marsh DOC G3 frac
 
 !KDOC1    = 0.035  ! 1/day from Di Toro pg 253, need to find for DOC not POC!!!!!
 !KDOC2    = 0.0018 ! " "
 !KDOC3	  = 0.00    ! "  "
 !KDON1    = 0.035  ! use the same as DOC 
 !KDON2    = 0.0018
 !KDON3    = 0.00
 !KDOP1    = 0.035
 !KDOP2    = 0.0018
 !KDOP3	  = 0.00
           
           
 !thetaDOC1 = 1.1 ! from Di Toro pg 253, same issue as above
 !thetaDOC2 = 1.150
 !thetaDOC3 = 0.0
 !thetaDON1 = 1.1
 !thetaDON2 = 1.150
 !thetaDON3 = 0.0
 !thetaDOP1 = 1.1
 !thetaDOP2 = 1.150
 !thetaDOP3 = 0.0
 !          
           
 !NMARSHNODE(:) = 0.0
           
 JMDOM     = 0.0    
           
! MARSH_FRAC = 1.0  ! fraction of each element that is covered in grass
! MARSH_MASS = 200 ! biomass in the marss g C/m^2
! 
!P_EX_RATE(:) = 0.01 ! fraction of the plant biomass that is exuded as DOM into the seds (1/day)
           
 JNCMDOC(:) = 0.0
 JCMDOC(:) =0.0
 
 !Anaerobic Layer carbon
 CDOC1TM1 = 20
 CDOC2TM1 = 15
 CDOC3TM1 = 3
 NCDOC1TM1 = 20 
 NCDOC2TM1 = 15
 NCDOC3TM1 =3

 !Aerobic Layer carbon
 CDOC11TM1 = 1
 CDOC21TM1 = 1
 CDOC31TM1 = 1
 NCDOC11TM1 = 1 
 NCDOC21TM1 = 1
 NCDOC31TM1 =1
         
!Anaerobic Layer Nitrogen
 CDON1TM1 = 0.2334 ! (CDOC1TM1*M_NtoC)
 CDON2TM1 = 0.1750 
 CDON3TM1 = 0.0350
 NCDON1TM1 = 0.2334
 NCDON2TM1 = 0.1750
 NCDON3TM1 = 0.0350
 
!Aerobic Layer carbon
 CDON11TM1 = 0.1
 CDON11TM1 = 0.1
 CDON31TM1 = 0.1
 NCDON11TM1 = 0.1 
 NCDON21TM1 = 0.1
 NCDON31TM1 = 0.1         

 CDOP11TM1 = 0.02334
 CDOP21TM1 = 0.01750
 CDOP31TM1 = 0.0350
 NCDOP11TM1 = 0.02334
 NCDOP21TM1 = 0.01750
 NCDOP31TM1 = 0.0350
           
 CDOP1TM1 = 0.02334
 CDOP2TM1 = 0.01750
 CDOP3TM1 = 0.0350
 NCDOP1TM1 = 0.02334
 NCDOP2TM1 = 0.01750
 NCDOP3TM1 = 0.0350
           
 !XKDOC1 = 0.1  
 !XKDOC2 = 0.1	  
 !XKDOC3 = 0.1	  
 !XKDON1 = 0.1	  
 !XKDON2 = 0.1	  
 !XKDON3 = 0.1	  
 !XKDOP1 = 0.1	  
 !XKDOP2 = 0.1   
 !XKDOP3 = 0.1
           
           
 JinCDOC1 = 0.0  
 JinCDOC2 = 0.0  
 JinCDOC3 = 0.0  
 JinNCDOC1 = 0.0  
 JinNCDOC2 = 0.0  
 JinNCDOC3 = 0.0   
           
 JinCDON1 = 0.0  
 JinCDON2  = 0.0 
 JinCDON3  = 0.0 
 JinNCDON1 = 0.0 
 JinNCDON2 = 0.0 
 JinNCDON3  = 0.0 
           
 JinCDOP1 = 0.0   
 JinCDOP2 = 0.0  
 JinCDOP3  = 0.0 
 JinNCDOP1 = 0.0 
 JinNCDOP2 = 0.0 
 JinNCDOP3  = 0.0 
           
 CDOC1(:)   = 0.0 
 CDOC2(:)   = 0.0 
 CDOC3(:)   = 0.0  
           
 NCDOC1(:)  = 0.0 
 NCDOC2(:)  = 0.0 
 NCDOC3(:)  = 0.0 
 
 CDOC11(:)   = 0.0 
 CDOC21(:)   = 0.0 
 CDOC31(:)   = 0.0  
           
 NCDOC11(:)  = 0.0 
 NCDOC21(:)  = 0.0 
 NCDOC31(:)  = 0.0 
           
 CDON1(:)   = 0.0 
 CDON2(:)   = 0.0 
 CDON3(:)   = 0.0 
           
 NCDON1(:)  = 0.0 
 NCDON2(:)  = 0.0 
 NCDON3(:)  = 0.0 
 
 CDON11(:)   = 0.0 
 CDON21(:)   = 0.0 
 CDON31(:)   = 0.0 
           
 NCDON11(:)  = 0.0 
 NCDON21(:)  = 0.0 
 NCDON31(:)  = 0.0 
           
 CDOP1(:)   = 0.0 
 CDOP2(:)   = 0.0
 CDOP3(:)   = 0.0
           
 NCDOP1(:)  = 0.0
 NCDOP2(:)  = 0.0 
 NCDOP3(:)  = 0.0 
 
 CDOP11(:)   = 0.0 
 CDOP21(:)   = 0.0
 CDOP31(:)   = 0.0
           
 NCDOP11(:)  = 0.0
 NCDOP21(:)  = 0.0 
 NCDOP31(:)  = 0.0 
           
           
 JDIAGDOP1  = 0.1 
 JDIAGDOP2  = 0.1 
 JDIAGDOP3  = 0.0
           
           
 JDIAGDOC1  = 0.1 
 JDIAGDOC2  = 0.1 
 JDIAGDOC3  = 0.0 
           
           
 JDIAGDON1  = 0.1 
 JDIAGDON2  = 0.1 
 JDIAGDON3  = 0.0 
           
           
 JDOCX_OUT1(:) = 0.0 
  JDOCX_OUT2(:) = 0.0 
 JDONX_OUT1(:) = 0.0 
 JDONX_OUT2(:) = 0.0 
 JDOPX_OUT1(:) = 0.0 
  JDOPX_OUT2(:) = 0.0 

END SUBROUTINE SED_DOM_INITIALIZE

!=========================*============================*==========================*

SUBROUTINE SED_DOM_INPUT
IMPLICIT NONE	

  seddom_gridin='inputs/seddom_grid.dat'
  FSEDDOM_GRID = SEDDOM_GRIDIN  !marsh input file

    OPEN(unit=1111,file = FSEDDOM_GRID)
    write(*,*)'opened the marshgrass grid file'

    READ(1111,*)NMARSHTOTAL,EXRATE_TEMPCON,TAU_crit, Resusp_mass, RoughHeight !get the total number of marsh cells and the temperature ramp for plant DOC exudation

    write(*,*)'Total Sed DOM node numbers:  ',NMARSHTOTAL
	write(*,*)'Critical Shear Stress = ',TAU_crit
	write(*,*)'Resuspension Flux of TSS = ', Resusp_mass
	write(*,*)'Roughness Height = ',RoughHeight
    !write(*,*)'read in the first line'

    WRITE(*,*) 'MTLOC = ', MTLOC
    WRITE(*,*) 'MGL = ', MGL

    IF (MGL /= NMARSHTOTAL) THEN
       WRITE(*,*)'THE NUMBER OF NODES FROM THE DOM DEFINITION FILE MUST COVER THE ENTIRE DOMAIN'
	   write(*,*)'Stopping Program'
	   STOP
    END IF
       ALLOCATE(NMNTEMP(MGL))
       ALLOCATE(MASSTEMP(MGL))
       ALLOCATE(FRACTEMP(MGL))
       ALLOCATE(EX_TEMP(MGL))
       ALLOCATE(PSTY_TEMP(MGL))
	   ALLOCATE(DNSTY_TEMP(MGL))
	   ALLOCATE(DOM_INITIAL_TEMP(MGL))
	   
	   ALLOCATE(NMARSHNODE(MTLOC))
       ALLOCATE(MARSH_FRAC(MTLOC))
       ALLOCATE(MARSH_MASS(MTLOC))
       ALLOCATE(P_EX_RATE(MTLOC))
	   ALLOCATE(PSTY(MTLOC))
	   ALLOCATE(DNSTY(MTLOC))
	    ALLOCATE(DOM_INITIAL(MTLOC))
	   
	   ALLOCATE(SED_ADV(MTLOC))
	   
read(1111,*)  ! skip two lines
read(1111,*)

    DO I=1,NMARSHTOTAL
       READ(1111,*) NMNTEMP(I),MASSTEMP(I),FRACTEMP(I),EX_TEMP(I),PSTY_TEMP(I),DNSTY_TEMP(I),DOM_INITIAL_TEMP(I)
    END DO



   IF(SERIAL) THEN
      NMARSHNODE= NMNTEMP
      MARSH_MASS = MASSTEMP
      MARSH_FRAC = FRACTEMP
      P_EX_RATE = EX_TEMP
	  PSTY = PSTY_TEMP
	  DNSTY=DNSTY_TEMP
	  DOM_INITIAL=DOM_INITIAL_TEMP
	  write(*,*)'Read the sed_dom grid file'
   ENDIF


#if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I = 1,MLOC !loop over local nodes
        NMARSHNODE(I)=NMNTEMP(NGID(I))
        MARSH_MASS(I)=MASSTEMP(NGID(I))
        MARSH_FRAC(I)=FRACTEMP(NGID(I))
        P_EX_RATE(I) = EX_TEMP(NGID(I))
		PSTY(I) = PSTY_TEMP(NGID(I))
		DNSTY(I) = DNSTY_TEMP(NGID(I))
		DOM_INITIAL(I) = DOM_INITIAL_TEMP(NGID(I))  ! DOM initial Condition

        !IF(MARSH_FRAC(I) > 0)THEN
        !      write(*,*)'Exudation Tempcon = ',EXRATE_TEMPCON
        !      write(*,*) 'NODE number', NGID(I)
		!	  WRITE(*,*)'NMARSHNODE1 = ', NMARSHNODE(I)
		!	  WRITE(*,*)'MARSH_MASS1 = ', MARSH_MASS(I)
		!	  WRITE(*,*)'MARSH_FRAC1 = ', MARSH_FRAC(I)
		!	  write(*,*)'Exudation rate =',P_Ex_RATE(I)
	    !endif

      ENDDO

      DO I = 1, NHN
         NMARSHNODE(I+MLOC) = NMNTEMP(HN_LST(I))
         MARSH_MASS(I+MLOC) = MASSTEMP(HN_LST(I))
         MARSH_FRAC(I+MLOC) = FRACTEMP(HN_LST(I))
         P_EX_RATE(I+MLOC)  = EX_TEMP(HN_LST(I))
		 PSTY(I+MLOC)=PSTY_TEMP(HN_LST(I))  ! The sediment porosity  should be < 1 if sediment is more porous than basal rate, > 1 if less porous than basal rate
		 DNSTY(I+MLOC)=DNSTY_TEMP(HN_LST(I))
		 DOM_INITIAL(I+MLOC)=DOM_INITIAL_TEMP(HN_LST(I))
      ENDDO
	 ENDIF
	 IF(MSR) THEN
	    WRITE (*,*) ' Read the sed_dom grid file'
     ENDIF !Serial loop
	 
!	MARSH_MASS=0.0 
!	MARSH_FRAC=0.0 
!       DOM_INITIAL=20.
# endif
	 

	 
    ! get the sediment water percolation velocity, based on the density of the sediments
DO I = 1,MLOC	

	 SED_ADV(I) = DNSTY(I)*SEDPERC_COEFF + SEDPERC_INT   ! SEdiment percolation coefficnet, multiplied by the sediment density to yield the sqrt of sediment percolation velocity
	 SED_ADV(I) = SED_ADV(I)**2/100.*24.   ! convert from sqrt(cm/hr) to m/day
	 SED_ADV(I) = SED_ADV(I)*PSTY(I)  ! correct for the porosity of the sediments, i.e. if highly dense and porous (gravel) high percolation, if highly dens and non-porous (marine sediments) percolation is low
										! 0<PSTY<1  0 is no pores, 1 is empty space 
 !        write(*,*)'sed adv = ,',NGID(i), SED_ADV(I)
     
END DO                 
!MARSH_MASS = 00.
!MARSH_FRAC = 1.0


       DEALLOCATE(NMNTEMP)
       DEALLOCATE(MASSTEMP)
       DEALLOCATE(FRACTEMP)
       DEALLOCATE(EX_TEMP)
	   DEALLOCATE(DNSTY_TEMP)

    !   P_EX_RATE(:) = 0.1  ! temperorarily have here due to bad read in file
      !EXRATE_TEMPCON = 1.15  ! new parameter for temperature control of plant exudation
	  
    !100 format(I8,F8.3,F)

	345 format(I4,3(F8.4))
 
	CLOSE (1111)

!Read in the marsh control file with the parameters
		seddom_conin='inputs/seddom_control.in'
        FSEDDOM_CON =SEDDOM_CONIN !Sediment DOM control file


	OPEN(unit=22,file=FSEDDOM_CON)

   write(*,*)'opened the control file';									      
	READ(22,222) thetaDOC1,thetaDOC2,thetaDOC3,KDOC1,KDOC2,KDOC3, &
 	             thetaDON1,thetaDON2,thetaDON3,KDON1,KDON2,KDON3, &
	             thetaDOP1,thetaDOP2,thetaDOP3,KDOP1,KDOP2,KDOP3, &
		         MFDOC1,MFDOC2,MFDOC3,	                     &!
				
			!	 WFDOM1,WFDOM2,WFDOM3,  &!
				 M_NtoC,M_PtoC,         &
		         NCDOCFRAC,CDOCFRAC,    &
                 DOM_diff,DDOM_TCON,  &
				 POCCF1,POCCF2,POCCF3
	
	REWIND(unit=22)
    READ(22,221) SED_DOM
	READ(22,221) HYDRO_FLAG
	CLOSE(22)
	
IF(MSR) THEN
 	WRITE(*,*)    'SED_DOM = ', SED_DOM
	write(*,*)    ' HYDROLysis FLAG = ',HYDRO_FLAG
	   WRITE(*,145) 'thetaDOC1 = ',thetaDOC1,    &
					'thetaDOC2 = ',thetaDOC2,    &
					'thetaDOC3 = ',thetaDOC3,    &
					'KDOC1 = ',KDOC1,             &
					'KDOC2 = ',KDOC2,			   &
					'KDOC3 = ',KDOC3, 			   &
 	                'thetaDON1 = ',thetaDON1,     & 
				    'thetaDON2 = ',thetaDON2,     &
					'thetaDON3 = ',thetaDON3,     &
					'KDON1 = ',KDON1,			   &
					'KDON2 = ',KDON2,			   &
					'KDON3 = ',KDON3, 			   &
					'thetaDOP1 = ',thetaDOP1,	   &
					'thetaDOP2 = ',thetaDOP2,	   &
					'thetaDOP3 = ',thetaDOP3,     &
				    'KDOP1 = ',KDOP1,			   &
					'KDOP2 = ',KDOP2,			   &
					'KDOP3 = ',KDOP3, 			   &
					'MFDOC1 = ',MFDOC1,		   &
					'MFDOC2 = ',MFDOC2,           &
					'MFDOC3 = ',MFDOC3,           &
					!'WFDOM1 = ',WFDOM1,			 &!
					!'WFDOM2 = ',WFDOM2,			&!
					!'WFDOM3 = ',WFDOM3,		   &!
					'M_NtoC = ',M_NtoC,		   &
					'M_PtoC = ',M_PtoC,     	   &
				    'NCDOCFRAC = ',NCDOCFRAC,     &
					'CDOCFRAC = ',CDOCFRAC,      &
					'DOM_diff = ',DOM_diff,      &
                    'DDOM_TCON = ',DDOM_TCON,   &
					'POCCF1 =  ',POCCF1,   &
					'POCCF2 =  ',POCCF2,   &
					'POCCF3 = ',POCCF3
ENDIF

	 ! Fractionate out the initial bulk DOM from the input file
!DOM_INITIAL=20.0	 
	 
DO I = 1,MLOC
       ! B Clark hardwired the initial dom to always be 20 
      DOM_INITIAL(I)=DOM_INITIAL(I)!*MARSH_FRAC(I) + 20.   ! Make the subtidal DOM
                                                        !always start at 20 mg C
                                                        !L^-1 

	CDOC1(I) = DOM_INITIAL(I) * CDOCFRAC* 0.1!hardwire that 0.1 percent of the DOM is labile, 0.2 is semilabile and 0.7 is refractory, can make it input later
	CDOC2(I) = DOM_INITIAL(I) * CDOCFRAC* 0.2!
	CDOC3(I) = DOM_INITIAL(I) * CDOCFRAC* 0.7!
    NCDOC1(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.1!
	NCDOC2(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.2!
	NCDOC3(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.7!

!NITROGEN
	CDON1(I) = DOM_INITIAL(I) * CDOCFRAC* 0.1 * M_NTOC
	CDON2(I) = DOM_INITIAL(I) * CDOCFRAC* 0.2 * M_NTOC
    CDON3(I) = DOM_INITIAL(I) * CDOCFRAC* 0.7 * M_NTOC

	NCDON1(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.1 * M_NTOC
	NCDON2(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.2 * M_NTOC
	NCDON3(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.7 * M_NTOC

!PHOSPHORUS
	CDOP1(I) = DOM_INITIAL(I) * CDOCFRAC* 0.1 * M_PTOC	
	CDOP2(I) = DOM_INITIAL(I) * CDOCFRAC* 0.2 * M_PTOC
	CDOP3(I) = DOM_INITIAL(I) * CDOCFRAC* 0.7 * M_PTOC
	
	NCDOP1(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.1 * M_PTOC
	NCDOP2(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.2 * M_PTOC
	NCDOP3(I) = DOM_INITIAL(I) * NCDOCFRAC* 0.7 * M_PTOC	
         
! set first layer same as second layer
  !      CDOC11(I)=CDOC1(I)
  !      CDOC21(I)=CDOC2(I)
  !      CDOC31(I)=CDOC3(I)
  !      NCDOC11(I)=NCDOC1(I)
  !      NCDOC21(I)=NCDOC2(I)
  !      NCDOC31(I)=NCDOC3(I)
  !      CDON11(I)=CDON1(I)
  !      CDON21(I)=CDON2(I)
  !      CDON31(I)=CDON3(I)
  !      NCDON11(I)=NCDON1(I)
  !      NCDON21(I)=NCDON2(I)
  !      NCDON31(I)=NCDON3(I)

   !     CDOP11(I)=CDOP1(I)
   !     CDOP21(I)=CDOP2(I)
   !     CDOP31(I)=CDOP3(I)
   !     NCDOP11(I)=NCDOP1(I)
   !     NCDOP21(I)=NCDOP2(I)
   !     NCDOP31(I)=NCDOP3(I)

ENDDO

!DOM_INITIAL=0.0

write(*,*)'DOM inital is ',DOM_INITIAL(352)


    221 format(6/,L8) 
	
	222 format(17/,7(F8.4,3x,F8.4,3x,F8.4,3/), &
 		          F8.4,3x,F8.4,3/,	       &  
	              F8.4,3x,F8.4,3/,   &  
	              F12.6,3x,F8.4,3/,  &
				  F8.4,3x,F8.4,3x,F8.4)     !,3/,         &
!				  F8.4)		
					 
	145 format(30(A12,F12.6,/)) 
 !open write the headers for the output file 
	
	    FSEDDOM_OUT='./outputs/sediment_DOM_stations.dat'
	    OPEN(UNIT=33,file = FSEDDOM_OUT)
		WRITE(33,300)'NODE','CC','NC',		&
						 'CN','NN',		&
						 'CP','NP'

	
        FWDOM_OUT = './outputs/watercolumn_DOM_stations.dat'  
        OPEN(UNIT = 34, file = FWDOM_OUT)
		WRITE(34,300)'NODE','CC','NC',		&
						 'CN','NN',		&
						 'CP','NP'
						 
		FJWDOM_OUT = './outputs/fluxes_DOM_stations.dat'
        OPEN(UNIT = 355, file = FJWDOM_OUT)
		WRITE(355,300)'NODE','CC','NC',		&
						 'CN','NN',		&
						 'CP','NP'
						 
	    300 format(A,2x,5(A,36X),A)
		301 format(A,6(7X,A))
		
END SUBROUTINE SED_DOM_INPUT

!!=====================*============================*=============================*


SUBROUTINE SED_DOM_OUTPUT

# if defined (MULTIPROCESSOR)
USE MOD_PAR, ONLY : NGID ,   &         ! globale node ID of local node
					HN_LST,  &            ! list of halo nodes
					NHN,     &           ! Number of halo nodes
					GATHER,   &
					NMAP
# endif
USE MOD_LIMS, ONLY :  NTLOC, &              !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
                      MTLOC, &             !!TOTAL OF LOCAL INTERNAL + HALO NODES
		              MLOC,  &			   
					  KBM1,	 &
					  MYID,  &
					  NPROCS

USE MOD_SED_DOM_EXCHANGE, ONLY : JNH4_GL_SHARE,JNO3_GL_SHARE, NH41_GL_SHARE, NH42_GL_SHARE, NO31_GL_SHARE, NO32_GL_SHARE

IMPLICIT NONE

  IF(SERIAL) THEN	

	DO I = 1,nstation
	
	   WRITE(33,404)Nstationnum_gl(I),CDOC1(Nstationnum_gl(I)),CDOC2(Nstationnum_gl(I)),CDOC3(Nstationnum_gl(I)),  &
					  NCDOC1(Nstationnum_gl(I)),NCDOC2(Nstationnum_gl(I)),NCDOC3(Nstationnum_gl(I)),    &
				      CDON1(Nstationnum_gl(I)),CDON2(Nstationnum_gl(I)),CDON3(Nstationnum_gl(I)),      &
					  NCDON1(Nstationnum_gl(I)),NCDON2(Nstationnum_gl(I)),NCDON3(Nstationnum_gl(I)),    &
					  CDOP1(Nstationnum_gl(I)),CDOP2(Nstationnum_gl(I)),CDOP3(Nstationnum_gl(I)),      &
					  NCDOP1(Nstationnum_gl(I)),NCDOP2(Nstationnum_gl(I)),NCDOP3(Nstationnum_gl(I)),  &
					  JNH4_GL_SHARE(Nstationnum_gl(I)),JNO3_GL_SHARE(Nstationnum_gl(I))
					  
	END DO
	
  ELSE
  
  !write(*,*)' writing SED DOM output'
# if defined (MULTIPROCESSOR)
  
! OUTPUT CARBON
     !write(*,*)'gathering the variable into a big array'
  
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDOC1(1),          VAR_OUT_GL(1,1))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDOC2(1)  ,         VAR_OUT_GL(1,2))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDOC3(1)  ,        VAR_OUT_GL(1,3))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDOC1(1)  ,        VAR_OUT_GL(1,4))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDOC2(1)  ,        VAR_OUT_GL(1,5))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDOC3(1)  ,        VAR_OUT_GL(1,6))
! OUTPUT NITROGEN
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDON1(1)  ,        VAR_OUT_GL(1,7))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDON2(1)  ,        VAR_OUT_GL(1,8))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDON3(1)  ,        VAR_OUT_GL(1,9))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDON1(1)  ,        VAR_OUT_GL(1,10))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDON2(1)  ,        VAR_OUT_GL(1,11))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDON3(1)  ,        VAR_OUT_GL(1,12))
!OUTPUT PHOSPHORUS
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDOP1(1)  ,        VAR_OUT_GL(1,13))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDOP2(1)  ,        VAR_OUT_GL(1,14))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CDOP3(1)  ,        VAR_OUT_GL(1,15))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDOP1(1)  ,        VAR_OUT_GL(1,16))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDOP2(1)  ,        VAR_OUT_GL(1,17))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NCDOP3(1)  ,        VAR_OUT_GL(1,18))
		
! NOW WRITE OUT THE WATER COLUMN DOM AT OUR STATION, makes it a little bit more accessible for the NASA project
	DO K=1,KBM1 
	!    !CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,LDOC(1,K),          LDOC_OUT_GL(1,K))
	!    !CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,RDOC(1,K),          RDOC_OUT_GL(1,K))
	!    !CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,LDON(1,K),          LDON_OUT_GL(1,K))
	!    !CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,RDON(1,K),          RDON_OUT_GL(1,K))
	!    !CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,LDOP(1,K),          LDOP_OUT_GL(1,K))
	!    !CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,RDOP(1,K),          RDOP_OUT_GL(1,K))	
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDOC1(1,K),          WCDOM_GL(1,K,1))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDOC2(1,K),          WCDOM_GL(1,K,2))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDOC3(1,K),          WCDOM_GL(1,K,3))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDOC1(1,K),          WCDOM_GL(1,K,4))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDOC2(1,K),          WCDOM_GL(1,K,5))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDOC3(1,K),          WCDOM_GL(1,K,6))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDON1(1,K),          WCDOM_GL(1,K,7))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDON2(1,K),          WCDOM_GL(1,K,8))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDON3(1,K),          WCDOM_GL(1,K,9))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDON1(1,K),          WCDOM_GL(1,K,10))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDON2(1,K),          WCDOM_GL(1,K,11))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDON3(1,K),          WCDOM_GL(1,K,12))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDOP1(1,K),          WCDOM_GL(1,K,13))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDOP2(1,K),          WCDOM_GL(1,K,14))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_CDOP3(1,K),          WCDOM_GL(1,K,15))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDOP1(1,K),          WCDOM_GL(1,K,16))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDOP2(1,K),          WCDOM_GL(1,K,17))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,WC_NCDOP3(1,K),          WCDOM_GL(1,K,18))
    END DO
	
	!! Gather the fluxes for output, can be used for comparison to literature values for fluxes in CBAY
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,1),          JDOM_OUT_GL(1,1))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,2),          JDOM_OUT_GL(1,2))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,3),          JDOM_OUT_GL(1,3))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,4),          JDOM_OUT_GL(1,4))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,5),          JDOM_OUT_GL(1,5))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,6),          JDOM_OUT_GL(1,6))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,7),          JDOM_OUT_GL(1,7))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,8),          JDOM_OUT_GL(1,8))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,9),          JDOM_OUT_GL(1,9))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,10),          JDOM_OUT_GL(1,10))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,11),          JDOM_OUT_GL(1,11))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,12),          JDOM_OUT_GL(1,12))
		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,13),          JDOM_OUT_GL(1,13))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,14),          JDOM_OUT_GL(1,14))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,15),          JDOM_OUT_GL(1,15))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,16),          JDOM_OUT_GL(1,16))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,17),          JDOM_OUT_GL(1,17))
	    CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JDOM_OUT(1,18),          JDOM_OUT_GL(1,18))
		
		
		! B Clark add in sum up of JDOM_OUT_GL  to get the time integrated flux in the output
		
		!JDOM_OUT_GL = JDOM_OUT_GL + JDOM_OUT_GL

	IF(MSR) THEN

		DO I = 1,nstation
		   WRITE(33,400)Nstationnum_gl(I), VAR_OUT_GL(Nstationnum_gl(I),1), VAR_OUT_GL(Nstationnum_gl(I),2), VAR_OUT_GL(Nstationnum_gl(I),3),  &
					   VAR_OUT_GL(Nstationnum_gl(I),4), VAR_OUT_GL(Nstationnum_gl(I),5), VAR_OUT_GL(Nstationnum_gl(I),6),    &
				       VAR_OUT_GL(Nstationnum_gl(I),7), VAR_OUT_GL(Nstationnum_gl(I),8), VAR_OUT_GL(Nstationnum_gl(I),9),      &
					   VAR_OUT_GL(Nstationnum_gl(I),10), VAR_OUT_GL(Nstationnum_gl(I),11), VAR_OUT_GL(Nstationnum_gl(I),12),    &
					   VAR_OUT_GL(Nstationnum_gl(I),13), VAR_OUT_GL(Nstationnum_gl(I),14), VAR_OUT_GL(Nstationnum_gl(I),15),      &
					   VAR_OUT_GL(Nstationnum_gl(I),16), VAR_OUT_GL(Nstationnum_gl(I),17), VAR_OUT_GL(Nstationnum_gl(I),18),   &
                       NH42_GL_SHARE(Nstationnum_gl(I)), NO32_GL_SHARE(Nstationnum_gl(I))    !,            &
   !  NO31_GL,            &


			!Write the output for water column DOM.
		  DO K = 1,KBM1
          !   !write(34,401)Nstationnum_gl(I),K,LDOC_OUT_GL(Nstationnum_gl(I),K),RDOC_OUT_GL(Nstationnum_gl(I),K),& 
          !   !                                 LDON_OUT_GL(Nstationnum_gl(I),K),RDON_OUT_GL(Nstationnum_gl(I),K), &
          !   !                                 LDOP_OUT_GL(Nstationnum_gl(I),K),RDOP_OUT_GL(Nstationnum_gl(I),K)
			 write(34,401)Nstationnum_gl(I),K,WCDOM_GL(Nstationnum_gl(I),K,1),WCDOM_GL(Nstationnum_gl(I),K,2),& !
                                              WCDOM_GL(Nstationnum_gl(I),K,3),WCDOM_GL(Nstationnum_gl(I),K,4), &!
                                              WCDOM_GL(Nstationnum_gl(I),K,5),WCDOM_GL(Nstationnum_gl(I),K,6), &!
											  WCDOM_GL(Nstationnum_gl(I),K,7),WCDOM_GL(Nstationnum_gl(I),K,8),& !
                                              WCDOM_GL(Nstationnum_gl(I),K,9),WCDOM_GL(Nstationnum_gl(I),K,10), &!
                                              WCDOM_GL(Nstationnum_gl(I),K,11),WCDOM_GL(Nstationnum_gl(I),K,12), &!
											  WCDOM_GL(Nstationnum_gl(I),K,13),WCDOM_GL(Nstationnum_gl(I),K,14),& !
                                              WCDOM_GL(Nstationnum_gl(I),K,15),WCDOM_GL(Nstationnum_gl(I),K,16), &!
                                              WCDOM_GL(Nstationnum_gl(I),K,17),WCDOM_GL(Nstationnum_gl(I),K,18)!, &!
          ENDDO
	
	       !write the output for fluxes between sediments and water column
			WRITE(355,400)Nstationnum_gl(I), JDOM_OUT_GL(Nstationnum_gl(I),1),JDOM_OUT_GL(Nstationnum_gl(I),2),JDOM_OUT_GL(Nstationnum_gl(I),3), &
											JDOM_OUT_GL(Nstationnum_gl(I),4),JDOM_OUT_GL(Nstationnum_gl(I),5),JDOM_OUT_GL(Nstationnum_gl(I),6), &	
											JDOM_OUT_GL(Nstationnum_gl(I),7),JDOM_OUT_GL(Nstationnum_gl(I),8),JDOM_OUT_GL(Nstationnum_gl(I),9), &
											JDOM_OUT_GL(Nstationnum_gl(I),10),JDOM_OUT_GL(Nstationnum_gl(I),11),JDOM_OUT_GL(Nstationnum_gl(I),12), &	
											JDOM_OUT_GL(Nstationnum_gl(I),13),JDOM_OUT_GL(Nstationnum_gl(I),14),JDOM_OUT_GL(Nstationnum_gl(I),15), &
											JDOM_OUT_GL(Nstationnum_gl(I),16),JDOM_OUT_GL(Nstationnum_gl(I),17),JDOM_OUT_GL(Nstationnum_gl(I),18), &!	
											JNH4_GL_SHARE(Nstationnum_gl(I)),JNO3_GL_SHARE(Nstationnum_gl(I))
       ENDDO
	ENDIF	
 
#endif
		
ENDIF ! SERIAL CHECK LOOP
	
	400 format (I,20(F16.8)) 
    401 format (2(I),3x,18(F16.8))
	402 format (I,3x,6(F12.4))
	404 format (I,20(F16.8))
	444 format (I4,F12.4)
	
	!OPEN(UNIT=44,file = 'fluxes.out')
	!    !WRITE(44,300)'CC','NC',		&
	!	!				 'CN','NN',		&
	!	!				 'CP','NP'
	!  	
	!	WRITE(44,500)JDIAGDOC1,JDIAGDOC2,JDIAGDOC3,JDIAGDON1,JDIAGDON2,JDIAGDON3,JDIAGDOP1,JDIAGDOP2,JDIAGDOP3
	!	
		500 format(F8.4,3x,F8.4,3x,F8.4,3x,F8.4,3x,F8.4,3x,F8.4,3x,      &
			   F8.4,3x,F8.4,3x,F8.4)
		
    !OPEN(unit=55,file='sed_DOM_tempcontrol.out')
	
	!WRITE(55,500)XKDOC1,XKDOC2,XKDOC3,XKDON2,XKDON2,XKDON3,XKDOP1,XKDOP2,XKDOP3
	!	
END SUBROUTINE SED_DOM_OUTPUT	

!!=====================*===========================*==============================*
	
SUBROUTINE SED_DOM_CALC

USE MOD_SED_DOM_EXCHANGE, ONLY: JDOCX_SHARE1,    &  !Diagenesis DOC flux  (mg C/m^2)
								 JDOCX_SHARE2,    &  !Diagenesis DOC flux  (mg C/m^2)
								JDONX_SHARE1,    &  !Diagenesis DON flux layer 1
		        				JDONX_SHARE2,    &  !Diagenesis DON flux layer 2
								JDOPX_SHARE1,    &   !Diagenesis DOP flux
								JDOPX_SHARE2,    &   !Diagenesis DOP flux
								SED_TEMP_SHARE,    &
								ZSED_SHARE,     &
								DLTS_SHARE,     &
								W2_SHARE ,   &
								FX_POC1TM1,  &   !  Added by B Clark Nov 2015 for hydrolysis
	                        	FX_POC2TM1,   &
	                         	FX_POC3TM1,  &
	                        	FX_PON1TM1,  &
	                          	FX_PON2TM1,   &
		                        FX_PON3TM1,    &	
		                        FX_POP1TM1,  &
		                        FX_POP2TM1,   &
	                        	FX_POP3TM1
								

IMPLICIT NONE


Real(SP), ALLOCATABLE, DIMENSION(:,:)::A! matrix for passing to LAPACK routines
INTEGER, ALLOCATABLE, DIMENSION(:,:)::IPIV! matrix for passing to LAPACK routines
!REAL(SP) :: A11, A12, A21, A22  ! matrix coefficients
REAL(SP) :: Diffusion_Distance ! distance between middle of bottom layer in
! in water column and sed layer 1
REAL(SP),ALLOCATABLE,DIMENSION(:,:):: B! input and output from matrix solver for rhs

INTEGER ::  INFO ! stuff for LU factorization output
allocate(A(2,2))
allocate(IPIV(2,2))
allocate(B(1,2))




!write(*,*)'W2_share = ',w2_share
!write(*,*)'Zsed_share = ',zsed_share
!write(*,*)'sed temp share = ',sed_temp_share
 IF(MARSH_SEDIMENTS) THEN

!!!ADD IN ALL OF OUR MARSH DOM PARAMETERIZATIONS HERE

	APCM=M_PtoC  !alpha P:C ratio 
	ANCM=M_NtoC  !alpha N:C ratio

 ELSE

!!!! ADD IN OUR ESTUARY PARAMETERIZATIONS HERE


 END IF

IF (SAV_CALC) THEN

   !MARSH_MASS=ROOT+LEAF+xxxxx  ! we will just assign the marsh mass to the root biomass, this makes it easy to start, we can modify to make 
                    !more complicated with leaf litter etc in the future

END IF

!***** get our global variables used for the calculations

  SED_TEMP = SED_TEMP_SHARE   !  =CTEMP(I)  !temperature of sediments
  ZSED = ZSED_SHARE     !depth of sediments
  DT = DLTS_SHARE          !time step in days
  
!  write(*,*)'DT = ',DT
!PERCOLATION_FLAG = .FALSE.  
  
  !!IF(PERCOLATION_FLAG) THEN
  !   TIDE_STAGE = ELNC1 - ELNC2
 ! ENDIF

  DO I=1,MLOC   ! loop over the local domain + NODES  !NMARSHTOTAL !number of marsh nodes

!get the concentration from the previous time step

!CARBON	, ANaerobic Layer
	CDOC1TM1  = CDOC1(I)  !get the DOC concentration from the previous time step (CDOC1(t-dt))
	CDOC2TM1  = CDOC2(I)
	CDOC3TM1  = CDOC3(I)
	NCDOC1TM1 = NCDOC1(I)
	NCDOC2TM1 = NCDOC2(I)
	NCDOC3TM1 = NCDOC3(I)
!CARBON	in aerobic layer
	CDOC11TM1  = CDOC11(I)  !get the DOC concentration from the previous time step (CDOC1(t-dt))
	CDOC21TM1  = CDOC21(I)
	CDOC31TM1  = CDOC31(I)
	NCDOC11TM1 = NCDOC11(I)
	NCDOC21TM1 = NCDOC21(I)
	NCDOC31TM1 = NCDOC31(I)
	

!NITROGEN in Anaerobic Layer
	CDON1TM1 = CDON1(I)
	CDON2TM1 = CDON2(I)
    CDON3TM1 = CDON3(I)
!	write(*,*)'CDON3TM1 = ',CDON3TM1
	NCDON1TM1 = NCDON1(I)
	NCDON2TM1 = NCDON2(I)
	NCDON3TM1 = NCDON3(I)
	
!Nitrogen	in aerobic layer
	CDON11TM1  = CDON11(I)  !get the DOC concentration from the previous time step (CDOC1(t-dt))
	CDON21TM1  = CDON21(I)
	CDON31TM1  = CDON31(I)
	NCDON11TM1 = NCDON11(I)
	NCDON21TM1 = NCDON21(I)
	NCDON31TM1 = NCDON31(I)	

!PHOSPHORUS
	
	!PHOSPHORUS in aerobic layer
	CDOP11TM1 =  CDOP11(I)	
	CDOP21TM1 =  CDOP21(I)
	CDOP31TM1 =  CDOP31(I)
	
	NCDOP11TM1 = NCDOP11(I)
	NCDOP21TM1 = NCDOP21(I)
	NCDOP31TM1 = NCDOP31(I)
	
	!P in anaerobic layer
	CDOP1TM1 =  CDOP1(I)	
	CDOP2TM1 =  CDOP2(I)
	CDOP3TM1 =  CDOP3(I)
	
	NCDOP1TM1 = NCDOP1(I)
	NCDOP2TM1 = NCDOP2(I)
	NCDOP3TM1 = NCDOP3(I)	

!NEW FLUX INTO THE MARSH FROM PLANTS
	
    EXRATE = GET_ZHTA(P_EX_RATE(I),SED_TEMP(I),EXRATE_TEMPCON)	
!EXRATE = 0.0  ! B Clark temp adtiion to get more DOm into marsh sediments	
	JMDOM = MARSH_MASS(I)*MARSH_FRAC(I)*EXRATE!*0.1  ! we will only loop over the marsh nodes, but will have MDOM be defined globally
!write(*,*)'JMDOM = ',JMDOM
!FRACTIONATION OF OUR DOM INTO COLORED AND NON COLORED PARTS
	JCMDOC(I)  = CDOCFRAC*JMDOM       !flux of cdom and ncdom from the marshes
	JNCMDOC(I) = NCDOCFRAC*JMDOM
!JMDOM=0.0
   ! if(ngid(i) ==  7656) THEN
!	write(*,*)'Marsh JMDOM =',JMDOM
!	endif
! now we put our colored and non-colored DOM into the labile (G1), semilabile (G2), refractory (G3)

!THE new input of DOC from our marsh grasses /sources  (for estuary seds this will be zero...)

!CARBON 
 	JinCDOC1 = JCMDOC(I)*MFDOC1/ZSED(I)
	JinCDOC2 = JCMDOC(I)*MFDOC2/ZSED(I)
	JinCDOC3 = JCMDOC(I)*MFDOC3/ZSED(I)

	JinNCDOC1 = JNCMDOC(I)*MFDOC1/ZSED(I)
	JinNCDOC2 = JNCMDOC(I)*MFDOC2/ZSED(I)
	JinNCDOC3 = JNCMDOC(I)*MFDOC3/ZSED(I)


!NITROGEN	
	JinCDON1 = JinCDOC1*ANCM
	JinCDON2 = JinCDOC2*ANCM
	JinCDON3 = JinCDOC3*ANCM
	
	JinNCDON1 = JinNCDOC1*ANCM
	JinNCDON2 = JinNCDOC2*ANCM
	JinNCDON3 = JinNCDOC3*ANCM


!PHOSPHORUS
	JinCDOP1 = JinCDOC1*APCM
	JinCDOP2 = JinCDOC2*APCM
	JinCDOP3 = JinCDOC3*APCM
	
	JinNCDOP1 = JinNCDOC1*APCM
	JinNCDOP2 = JinNCDOC2*APCM
	JinNCDOP3 = JinNCDOC3*APCM
	
!====== New Addition, B Clark April 2017, Two layer model diffusion from OLW into layer 1
	
	JinCDOC11 = JinCDOC1 + S_SHARE(I)*WC_CDOC1(I,KBM1)  ! plant input (g C m^-3 d^-1) + Diffusive flux (g m^-3 d^-1) using S = SOD/O20 from mod_sed
	JinCDOC21 = JinCDOC1 + S_SHARE(I)*WC_CDOC2(I,KBM1)
	JinCDOC31 = JinCDOC1 + S_SHARE(I)*WC_CDOC3(I,KBM1)	
	JinNCDOC11 = JinNCDOC1 + S_SHARE(I)*WC_NCDOC1(I,KBM1)  ! plant input (g C m^-3 d^-1) + Diffusive flux (g m^-3 d^-1) using S = SOD/O20 from mod_sed
	JinNCDOC21 = JinNCDOC1 + S_SHARE(I)*WC_NCDOC2(I,KBM1)
	JinNCDOC31 = JinNCDOC1 + S_SHARE(I)*WC_NCDOC3(I,KBM1)		

! Nitrogen
	JinCDON11 = JinCDON1 + S_SHARE(I)*WC_CDON1(I,KBM1)  ! plant input (g N m^-3 d^-1) + Diffusive flux (g m^-3 d^-1) using S = SOD/O20 from mod_sed
	JinCDON21 = JinCDON1 + S_SHARE(I)*WC_CDON2(I,KBM1)
	JinCDON31 = JinCDON1 + S_SHARE(I)*WC_CDON3(I,KBM1)	
	JinNCDON11 = JinNCDON1 + S_SHARE(I)*WC_NCDON1(I,KBM1)  ! plant input (g N m^-3 d^-1) + Diffusive flux (g m^-3 d^-1) using S = SOD/O20 from mod_sed
	JinNCDON21 = JinNCDON1 + S_SHARE(I)*WC_NCDON2(I,KBM1)
	JinNCDON31 = JinNCDON1 + S_SHARE(I)*WC_NCDON3(I,KBM1)		
	
! Phosphorus
	JinCDOP11 = JinCDOP1 + S_SHARE(I)*WC_CDOP1(I,KBM1)  ! plant input (g P m^-3 d^-1) + Diffusive flux (g m^-3 d^-1) using S = SOD/O20 from mod_sed
	JinCDOP21 = JinCDOP1 + S_SHARE(I)*WC_CDOP2(I,KBM1)
	JinCDOP31 = JinCDOP1 + S_SHARE(I)*WC_CDOP3(I,KBM1)	
	JinNCDOP11 = JinNCDOP1 + S_SHARE(I)*WC_NCDOP1(I,KBM1)  ! plant input (g P m^-3 d^-1) + Diffusive flux (g m^-3 d^-1) using S = SOD/O20 from mod_sed
	JinNCDOP21 = JinNCDOP1 + S_SHARE(I)*WC_NCDOP2(I,KBM1)
	JinNCDOP31 = JinNCDOP1 + S_SHARE(I)*WC_NCDOP3(I,KBM1)			

!Get the temperature adjusted decay rate from GET_ZHTA function
!CARBON	

	XKDOC1=GET_ZHTA(KDOC1,SED_TEMP(I),thetaDOC1)
	XKDOC2=GET_ZHTA(KDOC2,SED_TEMP(I),thetaDOC2)
	XKDOC3=GET_ZHTA(KDOC3,SED_TEMP(I),thetaDOC3)

!NITROGEN

	XKDON1=GET_ZHTA(KDON1,SED_TEMP(I),thetaDON1)
	XKDON2=GET_ZHTA(KDON2,SED_TEMP(I),thetaDON2)
	XKDON3=GET_ZHTA(KDON3,SED_TEMP(I),thetaDON3)

!PHOSPHORUS

	XKDOP1=GET_ZHTA(KDOP1,SED_TEMP(I),thetaDOP1)
	XKDOP2=GET_ZHTA(KDOP2,SED_TEMP(I),thetaDOP2)
	XKDOP3=GET_ZHTA(KDOP3,SED_TEMP(I),thetaDOP3)
	
	!!!*********************** ADD IN HYDROLYSIS FLUX FROM POC****************************
	!! Only occurs in the anaerobic layer
	
	IF (HYDRO_FLAG) THEN
	
	    CDOC1TM1=CDOC1TM1+POCCF1 * FX_POC1TM1(I)*DT
		CDOC2TM1=CDOC2TM1+POCCF2 * FX_POC2TM1(I)*DT
		CDOC3TM1=CDOC3TM1+POCCF3 * FX_POC3TM1(I)*DT
		
		NCDOC1TM1=NCDOC1TM1+(1-POCCF1) * FX_POC1TM1(I)*DT
		NCDOC2TM1=NCDOC2TM1+(1-POCCF2) * FX_POC2TM1(I)*DT
		NCDOC3TM1=NCDOC3TM1+(1-POCCF3) * FX_POC3TM1(I)*DT

			
		CDON1TM1=CDON1TM1+POCCF1 * FX_PON1TM1(I)*DT
		CDON2TM1=CDON2TM1+POCCF2 * FX_PON2TM1(I)*DT
		CDON3TM1=CDON3TM1+POCCF3 * FX_PON3TM1(I)*DT
		NCDON1TM1=NCDON1TM1+(1-POCCF1) * FX_PON1TM1(I)*DT
		NCDON2TM1=NCDON2TM1+(1-POCCF2) * FX_PON2TM1(I)*DT
		NCDON3TM1=NCDON3TM1+(1-POCCF3) * FX_PON3TM1(I)*DT
		
		CDOP1TM1=CDOP1TM1+POCCF1 * FX_POP1TM1(I)*DT
		CDOP2TM1=CDOP2TM1+POCCF2 * FX_POP2TM1(I)*DT
		CDOP3TM1=CDOP3TM1+POCCF3 * FX_POP3TM1(I)*DT
		NCDOP1TM1=NCDOP1TM1+(1-POCCF1) * FX_POP1TM1(I)*DT
		NCDOP2TM1=NCDOP2TM1+(1-POCCF2) * FX_POP2TM1(I)*DT
		NCDOP3TM1=NCDOP3TM1+(1-POCCF3) * FX_POP3TM1(I)*DT

	
	ENDIF
	
 !*********************Sediment Water percolation************************** 
! Added by B Clark Dec 2015
!	write(*,*)'Tide stage = ',TIDE_STAGE(I) 
  ! IF(PERCOLATION_FLAG .AND. (TIDE_STAGE(I) > 0.0001 .OR. TIDE_STAGE(I) < -0.0001) )THEN  ! put in a tolerance level for numerical elevation changes that can occur in the intertidal

 !         SED_ADV(I) = SED_ADV(I)*TIDE_STAGE(I)/abs(TIDE_STAGE(I))  ! changes the sign of the advection, where if the tide is falling (tide_stage < 0) advection is into sediments, if tide is rising, advection is out of the sediments)
!	write(*,*)'SED_ADV =',SED_ADV(I)	  
  ! ELSE	 
   
   !     SED_ADV(I) = 0.0
		 
  ! END IF
!write(*,*)'SED ADV = ',SED_ADV(I)
!*******************end Sediment Water percolation************************** 

!!************************ UPDATE CONCENTRATION *************************
	
!write(*,*)'WC_CDOC1 before = ',WC_CDOC1(I,KBM1)
	 W2M = W2_SHARE(I)           !DOC loss velocity is m/day
	 
     W2_Z=W2M/ZSED(I)

Diffusion_distance=(HSED1_SHARE(I)+D(I)*DZ(KBM1))*0.5

!Diffusion_distance=(HSED1_SHARE(I))+D(I)*DZ(KBM1))*0.5

	! write(*,*) 'ZSED = ',ZSED(I)
	 
      !DDOM = GET_ZHTA(DOM_diff,SED_TEMP(I),DDOM_TCON) ! temperature affects of diffusion
	 !  DDOM = DOM_diff
	  !DDOM = 0.0005
	!  write(*,*)'DDOM = ',DDOM
	!  DD = DDOM/ZSED(I) !+ SED_ADV(I)*0.1  ! add in the percolation to the diffusion coefficient to get total mixing  !diffusion is usually m^2/d, mixing rate = diffusion/depth = m/d
	 
	!  write(*,*) 'Diffusion = ',DD
	!  write(*,*)'DDOM =',DDOM
	  
	!  write(*,*)'DD =',DD
!#if defined(SEDIMENT_DOM) 
   SED_ADV=0.0 
       DD = 0.1*(KL12_SHARE(I)+SED_ADV(I))  /ZSED(I)   ! Mas Transfer velocity from mod_sed m/d,divided by depth of sediment
!	     ! decrease by order of magnitude to account of molecular weight of DOM
	  DD_Z=DD!/ZSED(I)   ! DEPTH INTEGRATED DIFFUSION RATE FOR SOLVING FOR SEDIMENT CONCS  1/DAY
!	  
!#endif	  
	  !write(*,*)'KL12_SHARE = ',KL12_SHARE
	 ! DD_Z = S_SHARE/HSED1!S_SHARE/ZSED(I)*0.1
	  
	  !S_SHARE = DD_Z

!COlORED CARBON	  
	!			   
	  !W2_Z = 0.0
	  !XKDOC1 = 0.3.
	  !S_SHARE = 0.0
	  !DD = 0.005
	  !DD_Z = 0.005
	  !JINCDOC11 = 100.
	  !JINCDOC1 = 100.

!CALL SGETRF(1,1,ALPHA1,1)

	   ALPHA1 = -(XKDOC1  + S_SHARE(I) )!+temp_absC)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	  
	   ALPHA2 = -(XKDOC1 + W2_Z + DD_Z +temp_absC)
	   
!write(*,*)'DD_Z = ',DD_Z
!   !A11 =  1 - ALPHA1*DT      ! time integrate everything
!   !A12 = - DD*DT
!   !A21=  - DD*DT
!   !A22 = 1  - ALPHA2*DT
!   	   !write(*,*)'A12 = ',A12
!	   write(*,*)'DD = ',DD
!	   write(*,*)'DT = ',DT
!   A(1,1) =  1 - ALPHA1*DT
!   A(1,2) = - DD*DT
!   A(2,1) = - DD*DT
!   A(2,2) = 1  - ALPHA2*DT
!   write(*,*) 'A before =',A

!   CALL SGETRF(2,2,A,2,IPIV,INFO)
!    write(*,*)'INFO = ',INFO
!   write(*,*)'A after = ',A
!  
!   B(1,1) = CDOC11TM1 + JinCDOC11*DT
!   B(1,2) = CDOC1TM1 + JinCDOC1*DT

!   write(*,*)'B Before = ',B 
!   CALL SGETRS('N',2,2,A,2,IPIV,B,2,INFO)
!  
!   write(*,*)'B after = ',B
!  write(*,*)'INFO = ',INFO
!CDOC11(I) = B(1,1)				   
!CDOC1(I) = B(1,2)

 CALL DDOM_SOLVER(CDOC11(I),CDOC1(I),CDOC11TM1, CDOC1TM1,ALPHA1, ALPHA2, DD, DT, JinCDOC11, JinCDOC1)
 
	 JWCDOC1(I) = -S_SHARE(I)*(WC_CDOC1(I,KBM1) - CDOC11(I))*Diffusion_distance

	   ALPHA1 = -( XKDOC2 + S_SHARE(I) )!+temp_absC)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOC2 + W2_Z + DD_Z + temp_absC)
		   
 CALL DDOM_SOLVER(CDOC21(I),CDOC2(I),CDOC21TM1, CDOC2TM1,ALPHA1, ALPHA2, DD, DT, JinCDOC21, JinCDOC2)

	 JWCDOC2(I) = -S_SHARE(I)*(WC_CDOC2(I,KBM1) - CDOC21(I))*Diffusion_distance
 
	   ALPHA1 = -( XKDOC3  + S_SHARE(I) )!+temp_absC)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOC3 + W2_Z + DD_Z + temp_absC)
		   
 CALL DDOM_SOLVER(CDOC31(I),CDOC3(I),CDOC31TM1, CDOC3TM1,ALPHA1, ALPHA2, DD,DT, JinCDOC31, JinCDOC3)

	 JWCDOC3(I) = -S_SHARE(I)*(WC_CDOC3(I,KBM1) - CDOC31(I))*Diffusion_distance

!Non-COlORED CARBON	  
	   ALPHA1 = -( XKDOC1  + S_SHARE(I))!+temp_absC )  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOC1 + W2_Z + DD_Z + temp_absC)
	
 CALL DDOM_SOLVER(NCDOC11(I),NCDOC1(I),NCDOC11TM1, NCDOC1TM1,ALPHA1, ALPHA2,DD, DT, JinNCDOC11, JinNCDOC1)

	 JWNCDOC1(I) = -S_SHARE(I)*(WC_NCDOC1(I,KBM1) - NCDOC11(I))*Diffusion_distance

	   ALPHA1 = -(XKDOC2  + S_SHARE(I))! +temp_absC)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOC2 + W2_Z + DD_Z + temp_absC)
		   
 CALL DDOM_SOLVER(NCDOC21(I),NCDOC2(I),NCDOC21TM1, NCDOC2TM1,ALPHA1, ALPHA2, DD,DT, JinNCDOC21, JinNCDOC2)
	 
	 JWNCDOC2(I) = -S_SHARE(I)*(WC_NCDOC2(I,KBM1) - NCDOC21(I))*Diffusion_distance

	   ALPHA1 = -(XKDOC3  + S_SHARE(I))! +temp_absC)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOC3 + W2_Z + DD_Z + temp_absC)
		   
 CALL DDOM_SOLVER(NCDOC31(I),NCDOC3(I),NCDOC31TM1, NCDOC3TM1,ALPHA1, ALPHA2, DD,DT, JinNCDOC31, JinNCDOC3)

	 JWNCDOC3(I) = -S_SHARE(I)*(WC_NCDOC3(I,KBM1) - NCDOC31(I))*Diffusion_distance

!COlORED Nitrogen	  
	   ALPHA1 = -(XKDON1 + S_SHARE(I) )!+temp_absN)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDON1 + W2_Z + DD_Z + temp_absN)

 CALL DDOM_SOLVER(CDON11(I),CDON1(I),CDON11TM1,CDON1TM1,ALPHA1, ALPHA2, DD, DT, JinCDON11, JinCDON1)

	 JWCDON1(I) = -S_SHARE(I)*(WC_CDON1(I,KBM1) - CDON11(I))*Diffusion_distance
	  
	   ALPHA1 = -(XKDON2  + S_SHARE(I) )!+temp_absN)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDON2 + W2_Z + DD_Z + temp_absN)
		   
 CALL DDOM_SOLVER(CDON21(I),CDON2(I),CDON21TM1, CDON2TM1,ALPHA1, ALPHA2, DD, DT, JinCDON21, JinCDON2)
	 
	 JWCDON2(I) = -S_SHARE(I)*(WC_CDON2(I,KBM1) - CDON21(I))*Diffusion_distance

	   ALPHA1 = -(XKDON3  + S_SHARE(I) )!+temp_absN)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDON3 + W2_Z + DD_Z + temp_absN)
		   
 CALL DDOM_SOLVER(CDON31(I),CDON3(I),CDON31TM1, CDON3TM1,ALPHA1, ALPHA2, DD,DT, JinCDON31, JinCDON3)
	 
	 JWCDON3(I) = -S_SHARE(I)*(WC_CDON3(I,KBM1) - CDON31(I))*Diffusion_distance

!Non-COlORED Nitrogen	  
	   ALPHA1 = -(XKDON1  + S_SHARE(I) )!+temp_absN)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDON1 + W2_Z + DD_Z + temp_absN)
	
 CALL DDOM_SOLVER(NCDON11(I),NCDON1(I),NCDON11TM1, NCDON1TM1,ALPHA1, ALPHA2,DD, DT, JinNCDON11, JinNCDON1)

	 JWNCDON1(I) = -S_SHARE(I)*(WC_NCDON1(I,KBM1) - NCDON11(I))*Diffusion_distance
	 
	   ALPHA1 = -(XKDON2 + S_SHARE(I) )!+temp_absN)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDON2 + W2_Z + DD_Z + temp_absN)
		   
 CALL DDOM_SOLVER(NCDON21(I),NCDON2(I),NCDON21TM1, NCDON2TM1,ALPHA1, ALPHA2, DD,DT, JinNCDON21, JinNCDON2)
	 
	 JWNCDON2(I) = -S_SHARE(I)*(WC_NCDON2(I,KBM1) - NCDON21(I))*Diffusion_distance
	 
	   ALPHA1 = -(XKDON3  + S_SHARE(I))! +temp_absN)  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDON3 + W2_Z + DD_Z + temp_absN)
		   
 CALL DDOM_SOLVER(NCDON31(I),NCDON3(I),NCDON31TM1, NCDON3TM1,ALPHA1, ALPHA2, DD,DT, JinNCDON31, JinNCDON3)
	 
	 JWNCDON3(I) = -S_SHARE(I)*(WC_NCDON3(I,KBM1) - NCDON31(I))*Diffusion_distance
	  	 
!!COlORED Phosphorous
!	 
	   ALPHA1 = -(XKDOP1  + S_SHARE(I) )  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOP1 + W2_Z + DD_Z )

 CALL DDOM_SOLVER(CDOP11(I),CDOP1(I),CDOP11TM1,CDOP1TM1,ALPHA1, ALPHA2, DD, DT, JinCDOP11, JinCDOP1)

	 JWCDOP1(I) = -S_SHARE(I)*(WC_CDOP1(I,KBM1) - CDOP11(I))*Diffusion_distance
	  
	   ALPHA1 = -(XKDOP2  + S_SHARE(I) )  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOP2 + W2_Z + DD_Z )	
	   
 CALL DDOM_SOLVER(CDOP21(I),CDOP2(I),CDOP21TM1, CDOP2TM1,ALPHA1, ALPHA2, DD, DT, JinCDOP21, JinCDOP2)
	 
	 JWCDOP2(I) = -S_SHARE(I)*(WC_CDOP2(I,KBM1) - CDOP21(I))*Diffusion_distance

	   ALPHA1 = -(XKDOP3 + S_SHARE(I) )  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOP3 + W2_Z + DD_Z )	
	   
 CALL DDOM_SOLVER(CDOP31(I),CDOP3(I),CDOP31TM1, CDOP3TM1,ALPHA1, ALPHA2, DD,DT, JinCDOP31, JinCDOP3)
	 
	 JWCDOP3(I) = -S_SHARE(I)*(WC_CDOP3(I,KBM1) - CDOP31(I))*Diffusion_distance

!Non-COlORED Phosphorous	  
	   ALPHA1 = -(XKDOP1  + S_SHARE(I) )  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOP1 + W2_Z + DD_Z )
	
 CALL DDOM_SOLVER(NCDOP11(I),NCDOP1(I),NCDOP11TM1, NCDOP1TM1,ALPHA1, ALPHA2,DD, DT, JinNCDOP11, JinNCDOP1)

	 JWNCDOP1(I) = -S_SHARE(I)*(WC_NCDOP1(I,KBM1) - NCDOP11(I))*Diffusion_distance
	 
	   ALPHA1 = -(XKDOP2  + S_SHARE(I) )  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOP2 + W2_Z + DD_Z )
		   
 CALL DDOM_SOLVER(NCDOP21(I),NCDOP2(I),NCDOP21TM1, NCDOP2TM1,ALPHA1, ALPHA2, DD,DT, JinNCDOP21, JinNCDOP2)
	 
	 JWNCDOP2(I) = -S_SHARE(I)*(WC_NCDOP2(I,KBM1) - NCDOP21(I))*Diffusion_distance
	 
	   ALPHA1 = -(XKDOP3  + S_SHARE(I) )  ! these are all the coefficients for C1(t+dt) units = 1/day, make negative because they are loss terms	
	   ALPHA2 = -(XKDOP3 + W2_Z + DD_Z )
	
 CALL DDOM_SOLVER(NCDOP31(I),NCDOP3(I),NCDOP31TM1, NCDOP3TM1,ALPHA1, ALPHA2, DD,DT, JinNCDOP31, JinNCDOP3)
	 
	 JWNCDOP3(I) = -S_SHARE(I)*(WC_NCDOP3(I,KBM1) - NCDOP31(I))*Diffusion_distance

!!********************** END UPDATE CONCENTRATION ***********************
	 
	 
!!********************** Acummulate fluxes over time **************************

!!#if(ACCUMULATE_FLUXES)	 
!	 	 
!	JDOM_OUT(I,1)= JDOM_OUT(I,1)+JWCDOC1(I) * DT * ART1(I)   ! time and area integrate the flux to get g m^-2 d^-1 ---> g 
!	JDOM_OUT(I,2)= JDOM_OUT(I,2)+JWCDOC2(I) * DT * ART1(I) 	 
!	JDOM_OUT(I,3)= JDOM_OUT(I,3)+JWCDOC3(I) * DT * ART1(I) 	 
!	JDOM_OUT(I,4)= JDOM_OUT(I,4)+JWNCDOC1(I) * DT * ART1(I) 	 
!	JDOM_OUT(I,5)= JDOM_OUT(I,5)+JWNCDOC2(I) * DT * ART1(I) 
!	JDOM_OUT(I,6)= JDOM_OUT(I,6)+JWNCDOC3(I) * DT * ART1(I) 
!	JDOM_OUT(I,7)= JDOM_OUT(I,7)+JWCDON1(I) * DT * ART1(I) 
!	JDOM_OUT(I,8)= JDOM_OUT(I,8)+JWCDON2(I) * DT * ART1(I) 
!	JDOM_OUT(I,9)= JDOM_OUT(I,9)+JWCDON3(I) * DT * ART1(I) 
!	JDOM_OUT(I,10)= JDOM_OUT(I,10)+JWNCDON1(I) * DT * ART1(I) 
!	JDOM_OUT(I,11)= JDOM_OUT(I,11)+JWNCDON2(I) * DT * ART1(I) 
!	JDOM_OUT(I,12)= JDOM_OUT(I,12)+JWNCDON3(I) * DT * ART1(I) 
!	JDOM_OUT(I,13)= JDOM_OUT(I,13)+JWCDOP1(I) * DT * ART1(I) 
!	JDOM_OUT(I,14)= JDOM_OUT(I,14)+JWCDOP2(I) * DT * ART1(I) 
!	JDOM_OUT(I,15)= JDOM_OUT(I,15)+JWCDOP3(I) * DT * ART1(I) 
!	JDOM_OUT(I,16)= JDOM_OUT(I,16)+JWNCDOP1(I) * DT * ART1(I) 
!	JDOM_OUT(I,17)= JDOM_OUT(I,17)+JWNCDOP2(I) * DT * ART1(I) 
!	JDOM_OUT(I,18)= JDOM_OUT(I,18)+JWNCDOP3(I) * DT * ART1(I) 
!!	write(*,*) 'JDOM_Out = ',JDOM_OUT(I,1)
!	
!#else	
	
	JDOM_OUT(I,1)= JWCDOC1(I) !* DLTS * ART1   !get the instantaneous flux accross the interface for output
	JDOM_OUT(I,2)= JWCDOC2(I)! * DLTS * ART1 	 
	JDOM_OUT(I,3)= JWCDOC3(I) !* DLTS * ART1 	 
	JDOM_OUT(I,4)= JWNCDOC1(I) !* DLTS * ART1 	 
	JDOM_OUT(I,5)=  JWNCDOC2(I)! * DLTS * ART1 
	JDOM_OUT(I,6)=  JWNCDOC3(I) !* DLTS * ART1 
	JDOM_OUT(I,7)= JWCDON1(I)! * DLTS * ART1 
	JDOM_OUT(I,8)= JWCDON2(I)! * DLTS * ART1 
	JDOM_OUT(I,9)= JWCDON3(I)! * DLTS * ART1 
	JDOM_OUT(I,10)= JWNCDON1(I)! * DLTS * ART1 
	JDOM_OUT(I,11)= JWNCDON2(I) !* DLTS * ART1 
	JDOM_OUT(I,12)= JWNCDON3(I)! * DLTS * ART1 
	JDOM_OUT(I,13)= JWCDOP1(I) !* DLTS * ART1 
	JDOM_OUT(I,14)= JWCDOP2(I) !* DLTS * ART1 
	JDOM_OUT(I,15)= JWCDOP3(I) !* DLTS * ART1 
	JDOM_OUT(I,16)= JWNCDOP1(I) !* DLTS * ART1 
	JDOM_OUT(I,17)= JWNCDOP2(I)! * DLTS * ART1 
	JDOM_OUT(I,18)= JWNCDOP3(I) !* DLTS * ART1 

!#endif

!**********************COLLECT WATER COLUMN FLUXES **********************
	 
! Gather the fluxes that are going to the water column via diffusion 
! into labile and refractory portions
! THERE WILL EVENTUALLY BE COLORED AND NON-COLORED CARBON IN THE WATER COLUMN
! BUT FOR NOW ALL OF THE DOM WILL BE 'NON-COLORED'

! this is obsolete, now we can pass the flux directly to the water column	 
!!UNITS = g /m^2/day
!	 
!J_LDOC_SHARE(I)=(JWCDOC1+JWNCDOC1) +   &    ! ALL OF THE G1 CARBON IS LABILE IN THE WATER COLUMN
!				(JWCDOC2+JWNCDOC2)*WFDOM3!,   & ! SOME FRACTION OF G2 IS LABILE
!				
!J_RDOC_SHARE(I)=(JWCDOC2+JWNCDOC2)*(1-WFDOM3) +  &   ! SOME FRACTION OF G2 CARBON IS REFRACTORY
!				(JWCDOC3+JWNCDOC3)                   !IN THE WATER COLUMN	
!				
!J_LDON_SHARE(I)=(JWCDON1+JWNCDON1) +   &    !ORGANIC NITROGEN, SAME AS CARBON
!				(JWCDON2+JWNCDON2)*WFDOM3   
!				
!J_RDON_SHARE(I)=(JWCDON2+JWNCDON2)*(1-WFDOM3) +  &   
!				(JWCDON3+JWNCDON3)                   

!J_LDOP_SHARE(I)=(JWCDOP1+JWNCDOP1) +   &    !ORGANIC PHOSPHORUS, SAME AS CARBON
!				(JWCDOP2+JWNCDOP2)*WFDOM3   
!				
!J_RDOP_SHARE(I)=(JWCDOP2+JWNCDOP2)*(1-WFDOM3) +  &   
!				(JWCDOP3+JWNCDOP3)       


!******************* END COLLECT WATER COLUMN FLUXES **********************				

!!************************** FLUX CALCULATION ***************************

!The diagenetic flux, we pass out MARSH_JDOCX_OUT to the sediment moduel to get the 	
!*************** layer 1 Aerobic Layer ammonium flux (no carbon flux in aerobic layer in model)
	! we will multiply by layer thickness in mod_sed where we can have the updated sediment layer 1 depth

!NITROGEN remineralization

	JDIAGDON1 = XKDON1 * (CDON11TM1 + NCDON11TM1)  !pass to G1 nitrogen diagenesis in layer 1 g N m^-3 d^-1
	JDIAGDON2 = XKDON2 * (CDON21TM1 + NCDON21TM1) !pass to G2 nitrogen diagenesis
	JDIAGDON3 = XKDON3 * (CDON31TM1 + NCDON31TM1)   !pass to G3 nitrogen diagenesis

	JDONX_OUT1(I)=JDIAGDON1+JDIAGDON2+JDIAGDON3 ! this goes to mod_sed (g N/m^3/day)
! Carbon remineralization
	JDIAGDOC1 = XKDOC1 * (CDOC11TM1 + NCDOC11TM1)  !pass to G1 carbon diagenesis
	JDIAGDOC2 = XKDOC2 * (CDOC21TM1 + NCDOC21TM1)  !pass to G2 carbon diagenesis
	JDIAGDOC3 = XKDOC3 * (CDOC31TM1 + NCDOC31TM1)   !pass to G3 carbon diagenesis
		
	JDOCX_OUT1(I)=JDIAGDOC1+JDIAGDOC2+JDIAGDOC3 ! this goes to mod_sed (g C/m^3/day)

!Phosphorus remineralization
	JDIAGDOP1 = XKDOP1 * (CDOP11TM1 + NCDOP11TM1)  !pass to G1 phosphorus diagenesis
	JDIAGDOP2 = XKDOP2 * (CDOP21TM1 + NCDOP21TM1)   !pass to G2 phosphorus diagenesis
	JDIAGDOP3 = XKDOP3 * (CDOP31TM1 + NCDOP31TM1)   !pass to G3 phosphorus diagenesis

	JDOPX_OUT1(I)=JDIAGDOP1+JDIAGDOP2+JDIAGDOP3 ! this goes to mod_sed  (g P/m^3/day)
		
!************ Layer 2 Anaerobic Layer	
!CARBON
	JDIAGDOC1 = XKDOC1 * (CDOC1TM1 + NCDOC1TM1) * ZSED(I)  !pass to G1 carbon diagenesis
	JDIAGDOC2 = XKDOC2 * (CDOC2TM1 + NCDOC2TM1) * ZSED(I)  !pass to G2 carbon diagenesis
	JDIAGDOC3 = XKDOC3 * (CDOC3TM1 + NCDOC3TM1) * ZSED(I)  !pass to G3 carbon diagenesis

	JDOCX_OUT2(I)=JDIAGDOC1+JDIAGDOC2+JDIAGDOC3 ! this goes to mod_sed (g C/m^2/day)

!NITROGEN

	JDIAGDON1 = XKDON1 * (CDON1TM1 + NCDON1TM1) * ZSED(I)  !pass to G1 nitrogen diagenesis
	JDIAGDON2 = XKDON2 * (CDON2TM1 + NCDON2TM1) * ZSED(I)  !pass to G2 nitrogen diagenesis
	JDIAGDON3 = XKDON3 * (CDON3TM1 + NCDON3TM1) * ZSED(I)  !pass to G3 nitrogen diagenesis

	JDONX_OUT2(I)=JDIAGDON1+JDIAGDON2+JDIAGDON3 ! this goes to mod_sed (g N/m^2/day)	

!PHOSPHORUS

	JDIAGDOP1 = XKDOP1 * (CDOP1TM1 + NCDOP1TM1) * ZSED(I)  !pass to G1 phosphorus diagenesis
	JDIAGDOP2 = XKDOP2 * (CDOP2TM1 + NCDOP2TM1) * ZSED(I)  !pass to G2 phosphorus diagenesis
	JDIAGDOP3 = XKDOP3 * (CDOP3TM1 + NCDOP3TM1) * ZSED(I)  !pass to G3 phosphorus diagenesis

	JDOPX_OUT2(I)=JDIAGDOP1+JDIAGDOP2+JDIAGDOP3 ! this goes to mod_sed  (g P/m^2/day)
	
!!***********************END FLUX CALCULATION ****************************
	
ENDDO



   ! UPDATE THE SHARED ARRAYS FOR PASSING TO THE BIG SEDIMENT MODULE
! Take note of unit difference between layer 1 and layer 2
! Layer 1 is depth integrated in mod_sed with the depth of the aerobic layer
	JDOCX_SHARE1 = JDOCX_OUT1    !Diagenesis DOC flux  (g C m^-3 d^-1)
	JDOCX_SHARE2 = JDOCX_OUT2    !Diagenesis DOC flux  (g C m^-2 d^-1)
    JDONX_SHARE1 = JDONX_OUT1	   !Diagenesis DON flux layer 1 (g N m^-3 d-1
    JDONX_SHARE2 = JDONX_OUT2	   !Diagenesis DON flux (g N m-2 d^-1)
	JDOPX_SHARE1 = JDOPX_OUT1    !Diagenesis DOP flux (g P m^-3 d^-1)
	JDOPX_SHARE2 = JDOPX_OUT2    !Diagenesis DOP flux (g P m^-2 d^-1)

END SUBROUTINE SED_DOM_CALC

!!!!!=========*===========***** FUNCTIONS *****===========*============*

!*************************DDOM SOLVER ******************************

SUBROUTINE DDOM_SOLVER(C1,C2,C1TM1,C2TM1,CC1,CC2,DD,DT_IN,J_IN1, J_IN2)
 use mod_prec, only: SP
IMPLICIT none

  !Implicitly Solve the sediment and water concentration  for DOM
  !initially start with a constant diffusion rate (non- temperature dependent) and a
  ! formulation that doesn't factor in porosity of the sediments

   REAL(SP), INTENT(OUT) :: C1, C2 ! concentrations(t+dt) of sediment(2) and sediment(1)
   REAL(SP), INTENT(IN)  :: C1TM1, C2TM1  !concentrations of sediment(2) and sediment(1) at the previous time step
   REAL(SP), INTENT(IN)	 :: J_IN1, J_IN2 ! flux into the seds from outside sources in layer 1 and layer 2 (plants, OLW...), g /m^3/day,
   REAL(SP), INTENT(IN)  :: DT_IN ! The time step, days
   
   REAL(SP) :: TDOC_IN1, &  ! flux into layer 1 from plant input and diffusion from overlying water column, 
								!g /m^3/d, + DOC from previous time step
				
				TDOC_IN2 ! flux from plants plus the DOC from the previous time step
				
   REAL(SP) :: CC1,CC2 ! reaction coefficients for sediment layer 1 and 2 concentration
   REAL(SP) :: DD  ! diffusion rate between layer 1 and lyaer 2 (1/d) = D/(Hsed^2)       

   REAL(SP)	::   GAMMA    ! Gamma = Determinant of systems of equations,used to find A^-1 and solve
   REAL(SP) :: A11, A12, A21, A22  ! matrix coefficients

   A11 =  1 - CC1*DT_IN       ! time integrate everything
   A12 = - DD*DT_IN
   A21=  - DD*DT_in
   A22 = 1  - CC2*DT_IN
   
   TDOC_IN1 = J_IN1*DT_IN + C1TM1! the total DOC into the solver from previous step and outside fluxes
   TDOC_IN2 = J_IN2*DT_IN + C2TM1! the total DOC into the solver from previous step and outside fluxes

   GAMMA = A11*A22 - A12*A21  ! the determinant
   
   C1 = (TDOC_IN1 * A22  -  TDOC_IN2 * A12 ) / GAMMA    ! = Flux into layer 1 * A11  - 

   C2 = (- A21 * TDOC_IN1 + A11 * TDOC_IN2) / GAMMA
   
   !  write(*,*)'C2tm1 in solver = ',C2TM1
   ! write(*,*)'C1tm1 in solver = ',C1TM1
!	
!	write(*,*) 'A21 = ',A21
!	write(*,*) 'TDOC_In = ',TDOC_in1
!	write(*,*) 'A11 = ',A11
!	write(*,*)  'C2TM1 = ',C2TM1
!	write(*,*)  'GAMMA = ',GAMMA
   !write(*,*) ' C2 in solver  = ',C2
   !write(*,*) ' C1 in solver  = ',C1

END SUBROUTINE DDOM_SOLVER

!*************************TEMPERATURE SOLVER ******************************

REAL FUNCTION GET_ZHTA(K_IN,TEMP,THETA_IN)
use mod_prec, only: SP
implicit none
	! this function solves the temperature reactoin dependence for the diagenetic flux
	! K_IN = REACTION RATE
	! TEMP = TEMPERATURE
	! THETA_IN = INITIAL THETA


	REAL(SP) :: K_IN 
	REAL(SP) :: THETA_IN  
	REAL(SP) :: TEMP 
    REAL(SP) :: TEMP20,TEMP202	

     TEMP20=TEMP-20.d0
     TEMP202=TEMP20/2.d0
	 !write(*,*) 'K_in = ',k_in
	 !write(*,*) 'theta in = ',theta_in
	 !write(*,*)'Temp20 = ',temp20   ! B Clark Debug
	 !write(*,*)'temp202 =',temp202
		
     GET_ZHTA=K_IN*THETA_IN**TEMP20
	!write(*,*)'temp function = ',get_zhta  ! B Clark Debug
   RETURN
  END FUNCTION 

END



