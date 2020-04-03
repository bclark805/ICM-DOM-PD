!************************************************************************
!**                                                                    **
!** Benthic Algage Module                                              **
!**                                                                    **
!**    Wen Long,9/19/2014                                              **
!**                                                                    **
!************************************************************************
!**                                                                    **
!** Inputs:                                                            **
!**                                                                    **
!**            Required inputs for benthic algae model                 **
!**                                                                    **
!** Outpus:                                                            **
!**                                                                    **
!**            Required inputs for benthic algae model                 **
!**                                                                    **
!************************************************************************

    MODULE MOD_BA
		USE MOD_PREC, ONLY : SP 
		USE MOD_LIMS, ONLY: MTLOC, MLOC
		USE MOD_WQM, ONLY : &
                AOCR,           &!
                BAOFN,          &!
                BB,             &!
                BENDO,          &!
                BENDOC,         &!
                BENNH4,         &!
                BENNO3,         &!
                BENPO4,         &!
                DLT,            &!
                DOXG,           &!
                KADPO4,         &!
                NH4,            &!
                NO3,            &!
                PO4,            &!
                SSI,            &!
                T!,             &!
		
        USE MOD_FILEINFO, ONLY : 		& !
				!,DIA					&!
				!,CBC 					&!
				!,S1					&!
				!,S2					&!
				!,S3					& !                 
				BAI				&!
				!,MET			&!
				,BAO	!		& !
				!,KEI			&!
				!,ATM			&!
				!,STL			& !
				!,AGR			& !
				!,SVI			& !
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
				!,UNIT_HIS		& !           
				!,CNAME			&!
				!,CONFN
				
USE MOD_HYDROVARS, ONLY:   &
   		!GRAV	&		!
		!,PI	&		!
		!,PI2	&		!
		!,ZERO	&		!
		!,ONE_THIRD	&	!
		!,NVG	&		!
		!,XG	&		!GLOBAL X-COORD AT NODE 
		!,YG	&		!GLOBAL X-COORD AT NODE 
		!,HG	&		!GLOBAL DEPTH AT NODE 
		!,XCG	&		!GLOBAL X-COORD AT FACE CENTER 
		!,YCG	&		!GLOBAL X-COORD AT FACE CENTER 
		!,VXMIN	&		!
		!,VYMIN	&		!
		!,VXMAX	&		!
		!,VYMAX	&		!
		!,XC	&		!X-COORD AT FACE CENTER 
		!,YC	&		!Y-COORD AT FACE CENTER
		!,VX	&		!X-COORD AT GRID POINT
		!,VY	&		!Y-COORD AT GRID POINT
		!,ART	&		!AREA OF ELEMENT
		!ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		!,ART2	&		!AREA OF ELEMENTS AROUND NODE
		!,NV	&		!NODE NUMBERING FOR ELEMENTS
		!,NBE	&		!INDICES OF ELMNT NEIGHBORS
		!,NTVE	&		!
		!,NTSN	&		!
		!,ISONB	&		!NODE MARKER = 0,1,2 
		!,ISBC	&		!
		!,ISBCE	&		!
		!,IEC	&		!
		!,IENODE &		!
		!,NBSN	&		!
		!,NIEC	&		!
		!,NTRG	&		!
		!,NBVE	&		!
		!,NBVT	&		!
		!,LISBCE_1	&	!LIST OF ELEMENTS WITH ISBCE=1
		!,LISBCE_2	&	!LIST OF ELEMENTS WITH ISBCE=2
		!,LISBCE_3	&	!LIST OF ELEMENTS WITH ISBCE=3
		!,DLTXC	&		!
		!,DLTYC	&		!
		!,DLTXYC	&	!
		!,DLTXE	&		!
		!,DLTYE	&		!
		!,DLTXYE	&	!
		!,SITAC	&		!
		!,SITAE	&		!
		!,XIJC	&		!
		!,YIJC	&		!
		!,XIJE	&		!
		!,YIJE	&		!
		!,EPOR	&		!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
		!,IBCGEO	&	!LOCAL GEOSTROPHIC FRICTION CORRECTION NODES

		!,Z	&			!SIGMA COORDINATE VALUE 
		!,ZZ	&		!INTRA LEVEL SIGMA VALUE
		DZ	&		!DELTA-SIGMA VALUE
		!,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
		!,H1	&		!BATHYMETRIC DEPTH 
		!,H	&			!BATHYMETRIC DEPTH 
		,D	!&			!CURRENT DEPTH 
		!,DT	&		!DEPTH AT PREVIOUS TIME STEP
		!,DT1	&		!DEPTH AT PREVIOUS TIME STEP
		!,EL	&		!CURRENT SURFACE ELEVATION
		!,ET	&		!SURFACE ELEVATION AT PREVIOUS TIME STEP
		!,DTFA	&		!ADJUSTED DEPTH FOR MASS CONSERVATION
		!,UU	&		!X-VELOCITY
		!,VV	&		!Y-VELOCITY
		!,UUT	&		!X-VELOCITY FROM PREVIOUS TIMESTEP
		!,VVT	&		!Y-VELOCITY FROM PREVIOUS TIMESTEP
		!,WWT	&		!Z-VELOCITY FROM PREVIOUS TIMESTEP
		!,WTST	&		!Vertical velocity in sigma from PREVIOUS TIMESTEP
		!,UARD_OBCNT	&!tykim
		!,XFLUX_OBCT	&!tykim
		!,DTFAT	&		!tykim
		!,TT_T	&		!tykim
		!,SALTT	&		!tykim
		!,WTS	&		!VERTICAL VELOCITY IN SIGMA SYSTEM
		!,UARD_OBCN	&	! tykim 
		!,XFLUX_OBC	&	! tykim 
		!,WTTS	&		!VERTICAL VELOCITY IN SIGMA SYSTEM 
		!,KH	&		!TURBULENT DIFFUSIVITY
		!,A1U	&		!
		!,A2U	&		!
		!,AWX	&		!
		!,AWY	&		!
		!,AW0	&		!
		!,VISCOFH	&	!
		!,UNC1	&		!
		!,VNC1	&		!
		!,WNC1	&		!
		!,WTSNC1	&		!
		!,UARD_OBCNNC1	&	!
		!,XFLUX_OBCNC1	&	!
		!,DTFANC1	&		!
		!,KHNC1	&		!
		!,TNC1	&		!
		!,SNC1	&		!
		!,ELNC1	&		!
		!,UNC2	&		!
		!,VNC2	&		!
		!,WNC2	&		!
		!,WTSNC2	&	!
		!,UARD_OBCNNC2	&!
		!,XFLUX_OBCNC2	&!
		!,DTFANC2	&	!
		!,KHNC2	&		!
		!,TNC2	&		!
		!,SNC2	&		!
		!,ELNC2	&		!
		!,num_hyd_ints	&!number of records in each hydrodynamics netcdf file
		!,TIME_MAP	&	!
		!,THOUR1	&	!SIMULATION TIME AT END OF CURRENT EXTERNAL STEP (IEXT) IN HOURS
		!,THOUR	&		!
		!,NCFILE_DIR	&!
		!,NCFILE_PREFIX	&!
		!,NCFILE_SUFFIX	&!
		!,NCFILE_NUMBER	&!
		!,FORMAT_STR	&!
		!,hydro_dir, 	&	! directory name where hydrodynamics results (netcdf) files are stored
		!,hydro_prefix, &	! prefix of file name, e.g. 'psm_'
		!,hydro_suffix	&	! suffix of filename, e.g. '.nc'
		!,hydro_filenumwidth, &	! number of digits in filename following hydro_prefix, e.g. 4 for psm_0002.nc
		!,hydro_filenumstart, &	! starting number of the file name in the digital part of the file name, e.g. 185 for psm_0185.nc
		!,hydro_Nrec	&		! number of records in each of hydrodynamics file
		!,hydro_dlt	&			! time step in hydrodynamics file (in seconds), e.g. 100 for 100sec
		!,t_his_start	&		!
		!,t_his_end	&			!
		!,t_his_dlt	&			!starting time, ending time, and interval of history outputs (days)
		!,Nstation	&			!
		!,NstationNum_GL	&	!maximum number of station is NstationMax!
		!,t_stn_start	&		!
		!,t_stn_end	&			!
		!,t_stn_dlt	&			!starting time, ending time, and interval of station outputs (days)
		!,STNFN	&				!file name for station output
		!,HISFN	&				!file name for history output
		!,HISFN_PREFIX	&		!prefix of history output file
		!,HISFN_EXT	&			!extention name of history output file
		!,HISFN_FINAL	&		! 
		!,HISFN_SPLIT_BYLEVEL	&!True or False for splitting history output in files level by level (default is .FALSE.)
		!,hydro_netcdf	&		!
		!,wqm_history	&		!
		!,wqm_stations	&		!
		!,IFNC	&				!file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise 
		!,NTRECNC	&			!time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file. 
		!,NTHYDRO				!overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read

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


   USE MOD_OWQ, ONLY:   				&!
						!INITKE,		&	!
						!INTKE,			&	!
						!NXKEI,			&	!
						!LAT, 			&	!
						!LONG,			&	!
						!GSTAR440,		&	!
						!PHISTAR676,	&	!
						!BNSTARNTU,		&	!
						!BNSTARCHL,     &	!
						!PSTARINT,		&	!
						!PSTAR440,		&	!
						!PSTARCOLOR,	&	!
						!SCDOM,			&	!
						!EZEROINT,    	&	!
						!TCHL2CHL,		&	!
						!DOFFSET,		&	!
						!DEPTHKE,		&	!
						!NWAVEL,		&	!
						!COLOR, 		&	!
						!TURB,			&	!
						!OPZ,		 	&	!
						!G1,			&	!
						!G2,			&	!
						!ALAMB, 		&	!
						!GLAMB,			&	!
						!PLAMB,			&	!
						!PHILAMB,       &	!
                        !BLAMB,			&	!
						!WAVEL,			&	!
						!EZERO,			&	!
						!GMUZERO,		&	!
						!IAVG,			&	!
						IATBOT_out !,		&	!
						!OWQ_ALLOC,		&	!
						!OWQ_DEALLOC,	&	!
						!OWQ_READ!,		&	!
						!LGHT_ATTN
						
	IMPLICIT NONE
    SAVE
   
	REAL(SP), ALLOCATABLE ::	DTPOC1_BA(:), &  !change of POC in sediments mgC/m^3 dring a sediment time step dlt_ba
							DTPOC2_BA(:), &	 !
							DTPOC3_BA(:), &  !
							DTPON1_BA(:), &  !change of PON in sediments mgN/m^3 dring a sediment time step dlt_ba
							DTPON2_BA(:), &  !
							DTPON3_BA(:), &  !
							DTPOP1_BA(:), &  !change of POP in sediments mgP/m^3 dring a sediment time step dlt_ba
							DTPOP2_BA(:), &
							DTPOP3_BA(:) 
							
	REAL(SP), ALLOCATABLE ::   	DTAG3CFL_BA(:),  & !chaneg of accumulated G3 Carbon flux 
												!during a sediment time step dlt_ba
							DTAG3NFL_BA(:),	 & !change of accumulated G3 Nitrogen flux 
												!during a sediment time step dlt_ba
							DTAG3PFL_BA(:)     !change of accumulated G3 Phosphorus flux 
												!during a sediment time step dlt_ba
	
    REAL(SP), ALLOCATABLE :: HBA(:)	!thickness of sediment layer that grows benthic algae (H1+H2 (m))

!
!variables that are used for as flags of model simulation 
!


      LOGICAL :: BALGAE_CALC            !flag for activating benthic algae
                                        !calculation

    
!
! variables that are local to benthic algae module
!
    REAL(SP) :: DLT_BA                     !benthic algae biomass simulation time step (day)
                                       !==DLT/86400, where DLT is in seconds
    	
	REAL(SP),DIMENSION(3) :: FRCPHB,    & !Fraction of POC generated for G1,G2,G3
	                                  !due to predation on benthic algae
	                     FRNPHB,    & !Fraction of PON genereated for G1,G2,G3
                                      !due to predation on benthic algae
                         FRPPHB       !Fraction of POP generated for G1,G2,G3
	                                  !due to predation on benthic algae

  
	
	
	!WLong moved from subroutine SED_READ
	CHARACTER(LEN=3)  :: BALC
						   
	REAL(SP) :: HBAALL    !depth of sediments with benthic algae (m)
						      
	INTEGER :: INTBAC !steady state flag of bentic algae simulation 
					  !(1 for steady state simulation, 0 for time varying)
	LOGICAL :: STEADY_STATE_BA !flag for steady state simulation of benthic algae
	
	LOGICAL :: BA_INITIALIZED !flag for if benthic algae model is initialized
	
	!benthic algae growth kinetics	
	!WLong moved from wqm_modules.F
	REAL(SP) ::  PMB, 		&!
			 ANCB, 		&!
             APCB, 		&!
             KTGB1, 	&!
             KTGB2, 	&!
             TMB,   	&!
             KESED, 	&!
             KEBALG, 	&!
             KHNB, 		&!
             KHPB, 		&!
             KHRB, 		&!
             BMRB, 		&!
             BPRB, 		&!
             KTBB,		&!
             TRB,  		&!
             BALGMIN, 	&!1
			 FNIB,      &!
             FPIB  		 

	
	REAL(SP) :: ALPHB,    &  !alpha of P-I curve in light calculation for benthic algae growth
	        CCHLB        !C/CHLA ratio of benthic algae gC/gCHLA
			
	LOGICAL ::   BA_OUTPUT    !Flag for benthic algae output
	
    !WLong has checked the units of these settling velocities
!	REAL(SP),ALLOCATABLE,DIMENSION(:) ::! WSSBNET,  & !netsettling rate of inorganic suspended sediments (m/d) 
									! WSLBNET,  & !net settling rate of LPOM (m/d)
									! WSRBNET,  & !net settling rate of RPOM (m/d)
									! WS1BNET,  & !net settling rate of alg 1 (m/d) (normal settling discounted by resuspension rate)
									! WS2BNET,  & !net settling rate of alg 2 (m/d)
									! WS3BNET,  & !net settling rate of alg 3 (m/d)
									! WSUBNET     !net settling rate of particulate biogenic silica (m/d)
							  
	!WLong moved from subroutine SED_READ
!	CHARACTER(LEN=8)  ::	SPVARB, &
!							PRINTB
							

	!WLong moved from subroutine SED_READ
!	CHARACTER(LEN=8)  :: SPVARLR,  &
!						 PRINTLR
	
!
!other local variables that are used in the benthic algae module
!

   

  	INTEGER ::    KWC!, !& 
				! IERR!, &
				! ITEMP
  
	CHARACTER(LEN=20) :: SSNAME_BA(17)  
	
	REAL(SP) :: NH4AVL,   &
			NO3AVL !,   &
!WLong			KETOTL    !never used
				
	REAL(SP)     :: DF, 		&          !dissolved fraction of a constituent 
				PO4AVL!, 	&          !available PO4 for benthic algae growth (=DF*PO4) in bottom layer of water column (gP/L ?)
				


    REAL(SP)     :: & !LOGICT, & 			   !
				IK
				
    REAL(SP)     :: & !PF,		&!
					!PPO4,	&!
 				FTB,	&!
 				PRNB,	&!
 				FRDOB	 !

	!Benthic algae growth limitation and kinetic fluxes
	!WLong moved this to wqm_sed.F.
	!benthic algae nutrient limitation and kinetic fluxes
	
	REAL(SP),ALLOCATABLE,DIMENSION(:) :: FIB,      &  !benthic algae growth light limitation (dimensionless)
	                                 NLB,      &  !benthic algae growth nitrogen limitation ***need to check NH4 limitation *** dimensionless)
	                                 PLB,      &  !benthic algae growth phosphorus limitation (dimensionless)
	                                 NPPB,     &  !net primary production of benthic algae (gC/m^2/day ?)
	                                 BBM,      &  !benthic algae biomass (gC/m^2 ?)
	                                 BLITE,    &  !benthic algae light (light at bottom of water column) (W/m^2)
	                                 BMB,      &  !benthic algae basal metabolism  after temperature control (1/day)
	                                 PB,       &  !primary production rate of benthic algae
	                                 PRB,      &  !predation rate of benthic algae (1/day)
	                                 BANH4,    &  !adjustment to benthic NH4 flux gN/m^2/day by benthic algae
	                                 BANO3,    &  !adjustment to benthic NO3 flux gN/m^2/day by benthic algae
	                                 BAPO4,    &  !adjustment to benthic PO4 flux gP/m^2/day by benthic algae
	                                 BADOC,    &  !adjustment to benthic DOC flux by benthic algae (gC/m^2/day)
	                                 BADO,     &  !adjustment to benthic DO  flux by benthic algae (gO2/m^2/day)
	                                 BAPOC,    &  !adjustment to benthic POC flux by benthic algae (gC/m^2/day)
	                                 BAPON,    &  !adjustment to benthic PON flux by benthic algae (gN/m^2/day)
	                                 BAPOP        !adjustment to benthic POP flux by benthic algae (gP/m^2/day)

	!WLong moved here from wqm_modules.F	
	REAL(SP) :: TINTIM_BA !time (days) for steady state simulation of sediments

   
	REAL(SP), ALLOCATABLE ::   BBM_GL(:)
  
   
   CONTAINS

   !subroutines:
		!subroutine BA_INIT()
		!subroutine BA_ALLOC()
		!subroutine BA_DEALLOC()
		!subroutine BA_READ()
		!subroutine BA_CALC()
		!subroutine BA_INT()
   
   
!********************************************************************************
!**                    S U B R O U T I N E   BA_INIT                           **
!********************************************************************************
  SUBROUTINE BA_INIT

    INTEGER :: I
	!WLong moved from wqm_modules.F
   
    SSNAME_BA(15) = 'Benthic Algae       '

	!allocate variables related to sediment diagenesis
    CALL BA_ALLOC
	
	!Set initial benthic algae sediment depth
	 DO I=1,MLOC
		HBA(I) = HBAALL*0.01        !WLong convert from cm to meter
     ENDDO
	 
	 BA_INITIALIZED = .TRUE. 
	 
  END SUBROUTINE BA_INIT

!********************************************************************************
!**                    S U B R O U T I N E   SED_ALLOC                         **
!********************************************************************************
  SUBROUTINE BA_ALLOC

	!WLong moved this here from wqm_modules.F
     ALLOCATE(FIB(MTLOC));                FIB   = 0.0
     ALLOCATE(NLB(MTLOC));                NLB   = 0.0
     ALLOCATE(PLB(MTLOC));                PLB   = 0.0
     ALLOCATE(NPPB(MTLOC));               NPPB  = 0.0
     ALLOCATE(BBM(MTLOC));                BBM   = 0.0
     ALLOCATE(BLITE(MTLOC));              BLITE = 0.0
     ALLOCATE(BMB(MTLOC));                BMB   = 0.0
     ALLOCATE(PB(MTLOC));                 PB    = 0.0
     ALLOCATE(PRB(MTLOC));                PRB   = 0.0
     ALLOCATE(BANH4(MTLOC));              BANH4 = 0.0
     ALLOCATE(BANO3(MTLOC));              BANO3 = 0.0
     ALLOCATE(BAPO4(MTLOC));              BAPO4 = 0.0
     ALLOCATE(BADOC(MTLOC));              BADOC = 0.0
     ALLOCATE(BADO(MTLOC));               BADO  = 0.0
     ALLOCATE(BAPOC(MTLOC));              BAPOC = 0.0
     ALLOCATE(BAPON(MTLOC));              BAPON = 0.0
     ALLOCATE(BAPOP(MTLOC));              BAPOP = 0.0
	 
	 ALLOCATE(DTPOC1_BA(MTLOC)); 	DTPOC1_BA=0.0;
	 ALLOCATE(DTPOC2_BA(MTLOC));	DTPOC2_BA=0.0;
	 ALLOCATE(DTPOC3_BA(MTLOC)); 	DTPOC3_BA=0.0;
	 ALLOCATE(DTPON1_BA(MTLOC));	DTPON1_BA=0.0;
	 ALLOCATE(DTPON2_BA(MTLOC));	DTPON2_BA=0.0;
	 ALLOCATE(DTPON3_BA(MTLOC));	DTPON3_BA=0.0;
	 ALLOCATE(DTPOP1_BA(MTLOC));	DTPOP1_BA=0.0;
	 ALLOCATE(DTPOP2_BA(MTLOC));	DTPOP2_BA=0.0;
	 ALLOCATE(DTPOP3_BA(MTLOC));	DTPOP3_BA=0.0;
	 
	 ALLOCATE(DTAG3CFL_BA(MTLOC)); DTAG3CFL_BA=0.0;
	 ALLOCATE(DTAG3NFL_BA(MTLOC)); DTAG3NFL_BA=0.0;
	 ALLOCATE(DTAG3PFL_BA(MTLOC)); DTAG3PFL_BA=0.0;
	 
	 ALLOCATE(HBA(MTLOC));	HBA=0.0;
			  
	 
END SUBROUTINE BA_ALLOC

!********************************************************************************
!**                    S U B R O U T I N E   BA_DEALLOC                        **
!********************************************************************************

SUBROUTINE BA_DEALLOC

     !WLong moved here from wqm_main.F
	 IF(ALLOCATED(FIB)) DEALLOCATE (FIB)
	 IF(ALLOCATED(NLB)) DEALLOCATE (NLB)
	 IF(ALLOCATED(PLB)) DEALLOCATE (PLB)
	 IF(ALLOCATED(NPPB)) DEALLOCATE (NPPB)
	 IF(ALLOCATED(BBM)) DEALLOCATE (BBM)
	 IF(ALLOCATED(BLITE)) DEALLOCATE (BLITE)
	 IF(ALLOCATED(BMB)) DEALLOCATE (BMB)
	 IF(ALLOCATED(PB)) DEALLOCATE (PB)
	 IF(ALLOCATED(PRB)) DEALLOCATE (PRB)
	 IF(ALLOCATED(BANH4)) DEALLOCATE (BANH4)
	 IF(ALLOCATED(BANO3)) DEALLOCATE (BANO3)
	 IF(ALLOCATED(BAPO4)) DEALLOCATE (BAPO4)
	 IF(ALLOCATED(BADOC)) DEALLOCATE (BADOC)
	 IF(ALLOCATED(BADO)) DEALLOCATE (BADO)
	 IF(ALLOCATED(BAPOC)) DEALLOCATE (BAPOC)
	 IF(ALLOCATED(BAPON)) DEALLOCATE (BAPON)
	 IF(ALLOCATED(BAPOP)) DEALLOCATE (BAPOP)

	 IF(ALLOCATED(DTPOC1_BA))DEALLOCATE(DTPOC1_BA)
	 IF(ALLOCATED(DTPOC2_BA))DEALLOCATE(DTPOC2_BA)
	 IF(ALLOCATED(DTPOC3_BA))DEALLOCATE(DTPOC3_BA)
	 IF(ALLOCATED(DTPON1_BA))DEALLOCATE(DTPON1_BA)
	 IF(ALLOCATED(DTPON2_BA))DEALLOCATE(DTPON2_BA)
	 IF(ALLOCATED(DTPON3_BA))DEALLOCATE(DTPON3_BA)
	 IF(ALLOCATED(DTPOP1_BA))DEALLOCATE(DTPOP1_BA)
	 IF(ALLOCATED(DTPOP2_BA))DEALLOCATE(DTPOP2_BA)
	 IF(ALLOCATED(DTPOP3_BA))DEALLOCATE(DTPOP3_BA)
	 
	 IF(ALLOCATED(DTAG3CFL_BA))DEALLOCATE(DTAG3CFL_BA)
	 IF(ALLOCATED(DTAG3NFL_BA))DEALLOCATE(DTAG3NFL_BA)
	 IF(ALLOCATED(DTAG3PFL_BA))DEALLOCATE(DTAG3PFL_BA)
	 
	 IF(ALLOCATED(HBA))DEALLOCATE(HBA)
	 
     IF(ALLOCATED(BBM_GL))        DEALLOCATE(BBM_GL)
	
	 
END SUBROUTINE BA_DEALLOC

!********************************************************************************
!**                    S U B R O U T I N E   B A _ R E A D                     **
!********************************************************************************
   SUBROUTINE BA_READ
#  if defined (MULTIPROCESSOR)
		USE MOD_PAR, ONLY:   &	  
			!EL_PID,	&	!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			!ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			!NLID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			!ELID_X,	&	!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			!NLID_X,	&	!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			!EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
			!NHE,       &   !			!NUMBER OF HALO ELEMENTS
			!HE_LST,	&	!(:)       	!!GLOBAL IDENTITIES OF HALO ELEMENTS
			!HE_OWN,	&	!(:)       	!!OWNER OF HALO ELEMENTS
			!NBN ,      &   !          	!!NUMBER OF BOUNDARY NODES
			!MX_MLT,    &   !          	!!MAX MULTIPLICITY OF BOUNDARY NODES
			!BN_LST, 	&	!	   		!!GLOBAL IDENTITY OF BOUNDARY NODES
			!BN_LOC, 	&	!			!!LOCAL IDENTITY OF BOUNDARY NODES
			!BN_MLT, 	&	!			!!MULTIPLICITY OF BOUNDARY NODES
			!BN_NEY, 	&	!			!!NODE OWNER LIST
			!NDE_ID, 	&	!			!! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
			NHN,  		&	!			!!NUMBER OF HALO NODES
			HN_LST!, 	&	!			!!LIST OF HALO NODES 
			!HN_OWN, 	&	!			!!PRIMARY OWNER OF HALO NODES
			!EC,		&	!
			!NC,		&	!
			!BNC, 		&	!
			!EMAP,		&	!
			!NMAP, 		&	!
			!PNE,		&	!(:)        !!NUMBER OF ELEMENTS IN EACH PROC
			!PNN,		&	!(:)        !!NUMBER OF NODES IN EACH PROC
			!PNHE,		&	!(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
			!PNBN,		&	!(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
			!PMBM,		&	!(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
			!PNHN,		&	!(:)       !!NUMBER OF HALO NODES IN EACH PROC
			!FILELOCK,   &
			!MYFILELOCK, &
			!FILELOCKED , &
			!NODE_MATCH, &
			!NODE_MATCH_ALL, &
			!EXCHANGE, &
			!EXCHANGE_ALL, &
			!GATHER, &
			!IGATHER, &
			!INIT_MPI_ENV,&!			
			!SORT,	&
			!PREAD,	&
			!PWRITE,&
			!PPRINT, &
			!LOCK, &
			!UNLOCK, &
			!GETLOC
#  endif  


   IMPLICIT NONE
   SAVE    

!***** Variable declarations

   CHARACTER(LEN=24) :: FRNAME_BA(3)
   INTEGER :: I, JG! , JT
  
   DATA FRNAME_BA                               &
        /'Benthic algal carbon    ',            &     !1   Benthic algae C
         'Benthic algal nitrogen  ',            &     !2   Benthic algae N
         'Benthic algal phosphorus'/                  !3   Benthic algae P

		 
   !open file for diagnostics on input format check
		 
   IF (BA_OUTPUT) THEN
     IF(MSR)THEN
		OPEN (BAO,FILE=BAOFN)
	 ENDIF
   ENDIF	 
!********************************************************************************
!**                                  Inputs                                    **
!********************************************************************************

 
! benthic algae 
   READ(BAI,1000,ERR=10100)  HBAALL, INTBAC   !depth of sediments with benthic algae,
											   !steady state flag
   READ(BAI,1015,ERR=10100)  PMB, ANCB, APCB, KTGB1, KTGB2, TMB
   READ(BAI,1015,ERR=10100)  ALPHB, CCHLB, KESED, KEBALG, KHNB, KHPB, KHRB
   READ(BAI,1015,ERR=10100)  BMRB, BPRB, KTBB, TRB, BALGMIN
   READ(BAI,1015,ERR=10100)  FNIB, FPIB
   READ(BAI,1010,ERR=10100)  FRPPHB  !Fractions of P in benthic algae predation to POM for G1, G2, G3
   READ(BAI,1010,ERR=10100)  FRNPHB  !Fractions of N in benthic algae predation to POM for G1, G2, G3
   READ(BAI,1010,ERR=10100)  FRCPHB  !Fractions of C in benthic algae predation to POM for G1, G2, G3
   
   CLOSE(BAI)
   
   STEADY_STATE_BA = INTBAC == 1
   
   BA_INITIALIZED = .FALSE.
   
!********************************************************************************
!**                                 Outputs                                    **
!********************************************************************************

   IF (BA_OUTPUT) THEN
     IF(MSR)THEN
			WRITE(BAO,2020) HBAALL
			IF (STEADY_STATE_BA) THEN
				WRITE(BAO,2022)
			ELSE
				WRITE(BAO,2025)
			ENDIF
			WRITE(BAO,2030)
			WRITE(BAO,2110)  FRNAME_BA(3),FRPPHB   
			WRITE(BAO,2110)  FRNAME_BA(2),FRNPHB
			WRITE(BAO,2110)  FRNAME_BA(1),FRCPHB
	
			WRITE(BAO,2252)  BALC
			WRITE(BAO,2242)  PMB, ANCB, APCB, KTGB1, KTGB2, TMB
			WRITE(BAO,2244)  ALPHB, CCHLB, KESED, KEBALG, KHNB, KHPB, KHRB
			WRITE(BAO,2246)  BMRB, BPRB, KTBB, TRB, BALGMIN
			WRITE(BAO,2248)  FNIB, FPIB
	 
		CLOSE(BAO)
     ENDIF
   ENDIF


!***** Input FORMAT'S

1000 FORMAT(:////8X,F8.0,I8)
1010 FORMAT(8F10.0)

1015 FORMAT(//8X,8F8.1)
 
!***** Output FORMAT'S
2020 FORMAT(/' ACTIVE LAYER DEPTH ',F8.3,' CM')
2022 FORMAT(/' STEADY-STATE VALUES OF BENTHIC ALGAE  COMPUTED'/)
2025 FORMAT(/' NON STEADY-STATE VALUES OF BENTHIC ALGAE COMPUTED'/)

2030 FORMAT(////33X,' BENTHIC ALGAE    I N I T I A L   C O N D I T ',   &
            'I O N S'/)
2110 FORMAT(6X,A24,11X,3F7.2)

2242 FORMAT(' PMB = ',F10.3/' ANCB = ',F10.3/' APCB = ',F10.3/          &
            ' KTGB1 = ',F10.3/' KTGB2 = ',F10.3/' TMB = ',F10.3)
2244 FORMAT(' ALPHB = ',F10.3/' CCHLB = ',F10.3/                        &
            ' KESED = ',F10.3/' KEBALG = ',F10.3,                       &
            ' KHNB = ',F10.3/' KHPB = ',F10.3/' KHRB = ',F10.3)
2246 FORMAT(' BMRB = ',F10.3/' BPRB = ',F10.3/' KTBB = ', F10.3/        &
            ' TRB = ',F10.3/ ' BALGMIN = ',F10.3)
2248 FORMAT(' FNIB = ',F10.3/' FPIB = ',F10.3)
2252 FORMAT(/' BENTHIC ALGAE ARE ',A3)
   RETURN


10100 CONTINUE

   IF(MSR)THEN
		IF (BA_OUTPUT)THEN
			WRITE(BAO,3010)			
			CLOSE(BAO)
		ENDIF
   ENDIF
3010 FORMAT(/' Read error in benthic algae input deck')
   STOP 'BA_READ'

   RETURN
   END SUBROUTINE BA_READ


!********************************************************************************
!**                    S U B R O U T I N E   B A _ C A L C                     **
!**                           Benthic Algae Calculations                       **
!********************************************************************************

   SUBROUTINE BA_CALC
   
   IMPLICIT NONE                     
   SAVE    
   INTEGER :: I
   !, J
   !, JSF, N, IERR
!   REAL(SP) :: SODMIN, SODMAX, SAVEFCT
!   REAL(SP) :: SOD,ZBOUT

   IF(.NOT.BA_INITIALIZED)THEN
		CALL BA_INIT()
   ENDIF
!***** Pass MOD_WQM time-step (in days) to sediment subr

   DLT_BA = DLT/86400.
   IF (STEADY_STATE_BA) TINTIM_BA = TINTIM_BA+DLT_BA  


!***** Update sediment concentrations

   DO I=1,MLOC


!WLong     TEMPD = CTEMP(I)
!WLong     STP20 = TEMPD-20.

     IF (BALGAE_CALC) THEN

!******* Benthic algae algorithms start here        

!******* Calculate mean light in algal mat

       BLITE(I) = IATBOT_out(I)*EXP(-KESED)/(KEBALG+1.0E-8)/BBM(I)   &
                   *(1. - EXP(-(KEBALG+1.0E-8)*BBM(I)))
        
!******* Temperature effects

       IF (T(I,KWC) < TMB) THEN
         FTB = EXP(-KTGB1*(T(I,KWC)-TMB)**2)
       ELSE
         FTB = EXP(-KTGB2*(TMB-T(I,KWC))**2)
       ENDIF
        
!******* Light effects

       IK = PMB*FTB/ALPHB
!WLong       FIB(I) = BLITE(I)/SQRT(IK*IK+BLITE(I)*BLITE(I)+1.0E-8)
        
!******* Nutrient limitations

! COMPUTE AVAILABLE AMMONIUM AND NITRATE

       NH4AVL = BENNH4(I)*DLT_BA + NH4(I,KWC)*D(I)*DZ(KWC)
       NH4AVL = MAX(0.,NH4AVL)
       NO3AVL = BENNO3(I)*DLT_BA + NO3(I,KWC)*D(I)*DZ(KWC)
       NO3AVL = MAX(0.,NO3AVL)

! COMPUTE NITROGEN LIMITATION 

       NLB(I) = (NH4AVL+NO3AVL)/(KHNB+NH4AVL+NO3AVL)

! COMPUTE NITROGEN PREFERENCE

       PRNB = NH4AVL*NO3AVL/((KHNB+NH4AVL)*(KHNB+NO3AVL))         &
             +NH4AVL*KHNB/((1.E-30+NH4AVL+NO3AVL)*(KHNB+NO3AVL))

!******* Phosphorus available for algal growth

       DF     = 1./(1.+KADPO4*SSI(I,KWC))
       PO4AVL = DF*PO4(I,KWC)*D(I)*DZ(KWC)
       PO4AVL = PO4AVL + BENPO4(I)*DLT_BA
       PO4AVL = MAX(0.,PO4AVL)
       PLB(I) = PO4AVL/(KHPB+PO4AVL)

!******* Base metabolism

! IF BIOMASS IS LESS THAN ALLOWED MINIMUM, SET METABOLISM TO ZERO

       IF (BBM(I) > BALGMIN) THEN
         BMB(I) = BMRB*EXP(KTBB*(T(I,KWC)-TRB))
       ELSE
         BMB(I) = 0.
       ENDIF

!******* Production

       PB(I) = PMB*FTB*AMIN1(FIB(I),NLB(I),PLB(I))/CCHLB

!******* Net primary production
      
       NPPB(I) = (PB(I)-BMB(I))*BBM(I)                     
 
!******* Predation

! IF BIOMASS IS LESS THAN ALLOWED MINIMUM, SET PREDATION TO ZERO

       IF (BBM(I) > BALGMIN) THEN
         PRB(I) = BBM(I)*BPRB*EXP(KTBB*(T(I,KWC)-TRB)) !WLong and Laura: should not have BBM here 
													   !for PRB should have units 1/d
       ELSE
         PRB(I) = 0.
       ENDIF

! ADJUST PREDATION SO BIOMASS DOESN'T GO NEGATIVE

       PRB(I) = MIN(PRB(I),PB(I)-BMB(I)+0.99/DLT_BA)              

!******* Compute effects of algal activity on benthic flux

       BANH4(I) = ANCB *(BMB(I)*FNIB - PRNB*PB(I)+ PRB(I)*FNIB) * BBM(I) 
       BANO3(I) = -(1. - PRNB) * PB(I) * ANCB * BBM(I)
       BAPO4(I) = APCB *(BMB(I)*FPIB - PB(I)+ PRB(I)*FPIB) * BBM(I) 
       FRDOB     = 1.-KHRB/(DOXG(I,KWC)+KHRB)
       BADO(I)  = ((1.3-0.3*PRNB)*PB(I)-FRDOB*BMB(I)) * AOCR*BBM(I)
       BADOC(I) = (1.-FRDOB)*BMB(I)*BBM(I)

!stuff to pass to water column due to change to NH4, NO3, PO4, DOC, DOXG in water column
       BENNH4(I) = BENNH4(I) + BANH4(I)	!gN/m^2/day
       BENNO3(I) = BENNO3(I) + BANO3(I) !gN/m^2/day
       BENPO4(I) = BENPO4(I) + BAPO4(I)	!gP/m^2/day
       BENDOC(I) = BENDOC(I) + BADOC(I)	!gC/m^2/day !LDOC
       BENDO(I)  = BENDO(I)  + BADO(I)  !gO2/m^2/day

! COMPUTE EFFECTS OF ALGAL ACTIVITY ON ORGANIC PARTICULATES (MG/M**3)
	!stuff to pass to sediment module
       BAPOC(I) = PRB(I)*BBM(I)
       BAPON(I) = ANCB*(1.-FNIB)*(BMB(I)+PRB(I))*BBM(I)
       BAPOP(I) = APCB*(1.-FPIB)*(BMB(I)+PRB(I))*BBM(I)
          
       !POC1 = POC1 + 1000. * BAPOC(I)*FRCPHB(1)*DLT_BA/H2   !1000*(gC/m^2/day) * (day) /(m) ==> mgC/m^3
       !POC2 = POC2 + 1000. * BAPOC(I)*FRCPHB(2)*DLT_BA/H2
       !POC3 = POC3 + 1000. * BAPOC(I)*FRCPHB(3)*DLT_BA/H2
       !PON1 = PON1 + 1000. * BAPON(I)*FRNPHB(1)*DLT_BA/H2
       !PON2 = PON2 + 1000. * BAPON(I)*FRNPHB(2)*DLT_BA/H2
       !PON3 = PON3 + 1000. * BAPON(I)*FRNPHB(3)*DLT_BA/H2
       !POP1 = POP1 + 1000. * BAPOP(I)*FRPPHB(1)*DLT_BA/H2
       !POP2 = POP2 + 1000. * BAPOP(I)*FRPPHB(2)*DLT_BA/H2
       !POP3 = POP3 + 1000. * BAPOP(I)*FRPPHB(3)*DLT_BA/H2

	   DTPOC1_BA(I) =  1000. * BAPOC(I)*FRCPHB(1)*DLT_BA/HBA(I)   !1000*(gC/m^2/day) * (day) /(m) ==> mgC/m^3
       DTPOC2_BA(I) =  1000. * BAPOC(I)*FRCPHB(2)*DLT_BA/HBA(I)
       DTPOC3_BA(I) =  1000. * BAPOC(I)*FRCPHB(3)*DLT_BA/HBA(I)
       DTPON1_BA(I) =  1000. * BAPON(I)*FRNPHB(1)*DLT_BA/HBA(I)
       DTPON2_BA(I) =  1000. * BAPON(I)*FRNPHB(2)*DLT_BA/HBA(I)
       DTPON3_BA(I) =  1000. * BAPON(I)*FRNPHB(3)*DLT_BA/HBA(I)
       DTPOP1_BA(I) =  1000. * BAPOP(I)*FRPPHB(1)*DLT_BA/HBA(I)
       DTPOP2_BA(I) =  1000. * BAPOP(I)*FRPPHB(2)*DLT_BA/HBA(I)
       DTPOP3_BA(I) =  1000. * BAPOP(I)*FRPPHB(3)*DLT_BA/HBA(I)

!******* Accumulate fluxes for steady-state computation

       IF (STEADY_STATE_BA) THEN
        
         !AG3CFL(I) = AG3CFL(I)+1000.*PRB(I)*FRCPHB(3)*BBM(I)*DLT_BA
         !AG3NFL(I) = AG3NFL(I)+1000.*PRB(I)*FRNPHB(3)*ANCB*BBM(I)*DLT_BA
         !AG3PFL(I) = AG3PFL(I)+1000.*PRB(I)*FRPPHB(3)*APCB*BBM(I)*DLT_BA
	   
         DTAG3CFL_BA(I) = 1000.*PRB(I)*FRCPHB(3)*BBM(I)*DLT_BA
         DTAG3NFL_BA(I) = 1000.*PRB(I)*FRNPHB(3)*ANCB*BBM(I)*DLT_BA
         DTAG3PFL_BA(I) = 1000.*PRB(I)*FRPPHB(3)*APCB*BBM(I)*DLT_BA
       ENDIF
 
!******* Change in benthic algal biomass

       BBM(I) = BBM(I) * (1. + DLT_BA*(PB(I)-BMB(I)-PRB(I)))

     ENDIF
   
   ENDDO


   RETURN
   END SUBROUTINE BA_CALC
   
   SUBROUTINE BA_INT()

   INTEGER :: I
   
   RETURN
   END SUBROUTINE BA_INT


   END MODULE MOD_BA
