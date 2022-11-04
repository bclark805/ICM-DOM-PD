MODULE MOD_KIN
	USE MOD_PREC, ONLY: SP
#if defined (MULTIPROCESSOR)      
 USE MOD_PAR, ONLY : NGID	
#endif
#if defined(ZOOPLANKTON)
	USE MOD_ZOOP, ONLY: 		&
					 CTSZ,    	&!carbon threshold for grazing (gC/m^3) for microzooplankton
                     CTLZ,    	&!carbon threshold for grazing (gC/m^3) for mesozooplankton
                     KHCSZ,   	&!prey density at which grazing is halved (gC/m^3) for microzooplankton
                     KHCLZ,   	&!prey density at which grazing is halved (gC/m^3) for mesozooplankton
                     MZEROSZ, 	&!mortality at zero dissolved oxygen (1/day) for microzooplankton
                     MZEROLZ,   &!mortality at zero dissolved oxygen (1/day) for mesozooplankton
					 UB1SZ, 	&!Utilization of algal group 1 by microzooplankton, range [0,1]
                     UB2SZ, 	&!Utilization of algal group 2 by microzooplankton, range [0,1]
                     UB3SZ, 	&!Utilization of algal group 3 by microzooplankton, range [0,1]
                     UB1LZ, 	&!Utilization of algal group 1 by mesozooplankton, range [0,1]
                     UB2LZ, 	&!Utilization of algal group 2 by mesozooplankton, range [0,1]
                     UB3LZ, 	&!Utilization of algal group 3 by mesozooplankton, range [0,1]                         
                     !UDSZ,  	&!Never used
                     ULDSZ, 	&!Utilization of labile dissolved organic carbon by microzooplankton, range [0,1]
                     URDSZ, 	&!Utilization of refractory dissolved organic carbon by microzooplankton, range [0,1]
                     ULPSZ, 	&!microzooplankton utilization of labile particulate organic carbon, range [0,1]
                     URPSZ, 	&!microzooplankton utilization of refractory particulate organic carbon, range [0,1]
                     ULLZ,  	&!mesozooplankton utilization of labile particulate organic carbon, range [0,1]
                     URLZ,  	&!mesozooplankton utilization of refractory particulate organic carbon, range [0,1]
                     USZLZ, 	&!mesozooplankton utilization of microzooplankton, range [0,1], (~1.0)
                     !TRSZ,  	&!reference temperature for microzooplankton metabolism (degC) 
                     !TRLZ,  	&!reference temperature for mesozooplankton metabolism (degC)
                     DOCRITSZ, 	&!dissolved oxygen below which mortality occurs for microzooplankton (gO2/m^3)
                     DOCRITLZ, 	&!dissolved oxygen below which mortality occurs for mesozooplankton (gO2/m^3)
                     ANCSZ,  	&!microzooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     ANCLZ,  	&!mesozooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     APCSZ,  	&!microzooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
                     APCLZ,  	&!mesozooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
                     AOCRSZ, 	&!ratio of oxygen consumed to microzooplankton carbon metabolized (gO2/gC) (~2.67)
                     AOCRLZ, 	&!ratio of oxygen consumed to mesozooplankton carbon metabolized (gO2/gC) (~2.67)
                     FRSASZ, 	&!fraction of microzooplankton silica recycled to dissolved silica pool, range [0,1]
                     FRSALZ,    &!fraction of mesozooplankton silica recycled to dissolved silica pool, range [0,1]
					 !FLDOCSZ,	&!fraction of microzooplankton carbon released to LDOC, range [0,1]
                     !FRDOCSZ,	&!fraction of microzooplankton carbon released to RDOC, range [0,1]
                     FLPOCSZ,	&!fraction of microzooplankton carbon released to LPOC, range [0,1]
                     FRPOCSZ,	&!fraction of microzooplankton carbon released to RPOC, range [0,1]
                     !FLDONSZ,	&!fraction of microzooplankton nitrogen released to LDON, range [0,1]
                     !FRDONSZ,	&!fraction of microzooplankton nitrogen released to RDON, range [0,1]
                     FLPONSZ,	&!fraction of microzooplankton nitrogen released to LPON, range [0,1]
                     FRPONSZ,	&!fraction of microzooplankton nitrogen released to RPON, range [0,1]
                     !FLDOPSZ,	&!fraction of microzooplankton phosphorus released to LDOP range [0,1]    
                     !FRDOPSZ,	&!fraction of microzooplankton phosphorus released to RDOP range [0,1]    
                     FLPOPSZ,	&!fraction of microzooplankton phosphorus released to LPOP range [0,1]    
                     FRPOPSZ,	&!fraction of microzooplankton phosphorus released to RPOP range [0,1]    
                     FNH4SZ, 	&!fraction of microzooplankton nitrogen recycled to DIN as NH4, range [0,1]
                     FPO4SZ, 	&!fraction of microzooplankton phosphorus recycled to DIP as PO4, range [0,1]
                     !FUREASZ,	&!Never used !!!                         
                   !  FLDOCLZ,	&!fraction of mesozooplankton carbon released to LDOC, range [0,1]
                    ! FRDOCLZ,	&!fraction of mesozooplankton carbon released to RDOC, range [0,1]
                     FLPOCLZ,	&!fraction of mesozooplankton carbon released to LPOC, range [0,1]
                     FRPOCLZ,	&!fraction of mesozooplankton carbon released to RPOC, range [0,1]
                   !  FLDONLZ,	&!fraction of mesozooplankton nitrogen released to LDON, range [0,1]
                   !  FRDONLZ,	&!fraction of mesozooplankton nitrogen released to RDON, range [0,1]
                     FLPONLZ,	&!fraction of mesozooplankton nitrogen released to LPON, range [0,1]
                     FRPONLZ,	&!fraction of mesozooplankton nitrogen released to RPON, range [0,1]
                   !  FLDOPLZ,	&!fraction of mesozooplankton phosphorus released to LDOP range [0,1]    
                    ! FRDOPLZ,	&!fraction of mesozooplankton phosphorus released to RDOP range [0,1]    
                     FLPOPLZ,	&!fraction of mesozooplankton phosphorus released to LPOP range [0,1]    
                     FRPOPLZ,	&!fraction of mesozooplankton phosphorus released to RPOP range [0,1]    
                     FNH4LZ,	&!fraction of mesozooplankton nitrogen recycled to DIN as NH4, range [0,1]
                     FPO4LZ,	&!fraction of mesozooplankton phosphorus recycled to DIP as PO4, range [0,1]
                     !FUREALZ,	&!Never used!!!
					 B1ASZ,   	&! 
                     B2ASZ,   	&!
                     B3ASZ,   	&!
					 LPOCASZ,	&!
                     RPOCASZ,	&!
                     PRASZ,  	&!
					 B1ALZ,   	&!
					 B2ALZ,   	&!
					 B3ALZ,   	&!
					 SZALZ,   	&!
					 LPOCALZ, 	&!
					 RPOCALZ, 	&!
					 PRALZ,   	&!
					 CLSZ,    	&! 
					 CLLZ,    	&!
					 RSZ,     	&!
					 RLZ,     	&!
					 RMAXSZ,  	&!
					 RMAXLZ,  	&!
					 BMSZ,    	&!
					 BMLZ,    	&!
					 BMRSZ,   	&!
					 BMRLZ,   	&!
					 MSZ,     	&!
					 MLZ,     	&!
					 PRSZLZ,  	&!
					 GSZ,     	&!
					 GLZ,     	&!
					 ESZ,     	&!
					 ELZ,     	&!
					 RFSZ,    	&!
					 RFLZ,    	&!
					 PRSZ,    	&!
					 PRLZ,    	&!
					 DOC1ASZ, 	&!
					 BPRSZ,   	&!
					 BPRLZ,   	&!
					 DOC2ASZ, 	&!
					 !LDOCSZ,  	&!
					 LPOCSZ,  	&!
					 RPOCSZ,  	&!
					 !LDOCLZ,  	&!
					 LPOCLZ,  	&!
					 RPOCLZ,  	&!
					 NH4SZ,   	&!
					 !LDONSZ,  	&!
					 LPONSZ,  	&!
					 RPONSZ,  	&!
					 NH4LZ,   	&!
					 !LDONLZ,  	&!
					 LPONLZ,  	&!
					 RPONLZ,  	&!
					 PO4SZ,   	&!
					 !LDOPSZ,  	&!
					 LPOPSZ,  	&!
					 RPOPSZ,  	&!
					 PO4LZ,   	&!
					 !LDOPLZ,  	&!
					 LPOPLZ,  	&!
					 RPOPLZ,  	&!
					 !RDOCSZ,  	&!
					 !RDONSZ,  	&!
					 !RDOPSZ,  	&!
					 !RDOCLZ,  	&!
					 !!RDONLZ,  	&!
					 !RDOPLZ,  	&!
					 PIB1SZ,  	&!
					 PIB2SZ,  	&!
					 PIB3SZ,  	&!
					 PIB1LZ,  	&!
					 PIB2LZ,  	&!
					 PIB3LZ,  	&!
					 B1SZ,    	&!
					 B2SZ,    	&!
					 B3SZ,    	&!
					 B1LZ,    	&!
					 B2LZ,    	&!
					 B3LZ,    	&!
					 DOSZ,    	&!
					 DOLZ,    	&!
					 SASZ,    	&!
					 SUSZ,    	&!
					 SALZ,    	&!
					 SULZ,    	&!
					 
					 !ACLSZ,	&!Accumulation of stuff for diagnostics
					 !ACLLZ,	&!
					 !ARSZ,		&!
					 !ARLZ,   	&!
					 !ABMSZ,	&!
					 !ABMLZ,	&!
					 !AMSZ,		&!
					 !AMLZ,   	&!
					 !APRSZLZ,	&!
					 !AGSZ,		&!
					 !AGLZ, 	&!
					 !ADOCSZ,	&!		
					 !APOCSZ,	&!
				 	 !ADOCLZ,	&!
					 !APOCLZ, 	&!
					 !ANH4SZ,	&!
					 !ADONSZ,	&!
					 !APONSZ,	&!
					 !ANH4LZ,	&!
					 !ADONLZ,	&!
					 !APONLZ,	&!
					 !APO4SZ,	&!
					 !ADOPSZ,	&!
					 !APOPSZ,	&!
					 !APO4LZ,	&!
					 !ADOPLZ,	&!
					 !APOPLZ,	&!
					 !APRSZ,	&!
					 !APRLZ,	&!
					 !APISZ,	&!
					 !APILZ,	&!
					 
					 !AB1SZ, 	&!
					 !AB2SZ,	&!
					 !AB3SZ,	&!
					 !AB1LZ,	&!
					 !AB2LZ,	&!
					 !AB3LZ,	&!
					 !ADOSZ,	&!
					 !ADOLZ,	&!
					 !ASASZ,	&!
					 !ASUSZ,	&!
					 !ASALZ,	&!
					 !ASULZ,	&!
					 
					 FTLZ,		&!
					 FTSZ,		&!
					 FTBMSZ,	&!
					 FTBMLZ,	&!
					 FTPRSZ,	&!
					 FTPRLZ,	&!
					 SZ,		&!
					 LZ,		&!
					 TVPRSZ,	&!
					 TVPRLZ
#endif					 
!USE MOD_SED_DOM_EXCHANGE, ONLY:  J_LDOC_SHARE,  &  ! THIS ARE THE FLUXES OF DOM FROM THE SEDIMENTS TO THE BOTTOM LAYER OF WATER COLUMN
!								   J_RDOC_SHARE,  &  ! g /m^2/day
!								   J_LDON_SHARE,  &
!								   J_RDON_SHARE,  &
!								   J_LDOP_SHARE,  &
!								   J_RDOP_SHARE
								   
USE MOD_SED_DOM, ONLY:  SED_DOM, TAU_crit, Resusp_mass
USE MOD_WQM, ONLY: ShearStress

		
										
!************************************************************************
!** additional paramters were made spatially variable these include:   **
!**     KDC,    KLC,    KRC,    KDN,    KLN,    KRN,    KDP,           **
!**     KLP,    KRP,   KSUA,   KCOD, KDCALG, KLCALG, KRCALG,           **
!**  KDNALG, KLNALG, KRNALG, KDPALG, KLPALG, KRPALG,    NTM,           **
!**      DL,      R        ***** MNOEL 2-20-93 *****                   **
!**  **** took out kfl_%%% subroutines that performed writing          **
!**       operations and included them in main program so they could   **
!**       be averaged over time periods instead of instantaneous       **
!**       output.  This process included vectorizing some variables.   **
!**                        ***** MNOEL 7-28-93 *****                   **
!**                                                                    **
!**                Kinetics Subroutines for CE-QUAL-IC                 **
!**                                                                    **
!**                            Version 1.0                             **
!**                         February 24, 1989                          **
!**                                                                    **
!**                    Water Quality Modeling Group                    **
!**                    U.S. Army Corps of Engineers                    **
!**                    Waterways Experiment Station                    **
!**                    Vicksburg, Mississippi 39180                    **
!**                                                                    **
!************************************************************************

	CONTAINS
		!Subroutines:
		!	Subroutine TEMPER() 
		!	Subroutine SOLIDS()
		!	Subroutine ZOOPL()
		!	Subroutine CARBON()
		!	Subroutine NITROG()
		!	Subroutine PHOSPH()
		!	Subroutine CODMND()
		!	Subroutine OXYGEN()
		!	Subroutine SILICA()
		!	Subroutine BEN_FLUX()
		!

!************************************************************************
!**                 S U B R O U T I N E   T E M P E R                  **
!************************************************************************

   SUBROUTINE TEMPER(DTT) 
   
   USE MOD_WQM, ONLY :  		&		 
                KT,             &!  heat transfer coefficient from METFN
                T,              &!	Temperature (degC)
                TE!,             &!	Equilibrium temperature (degC)

   USE MOD_LIMS, ONLY : &
				MLOC, 	&
				KBM1, 	&
				MTLOC
				
   USE MOD_HYDROVARS, ONLY: &
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
		!,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
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


   IMPLICIT NONE
   REAL(SP) :: RHO,CP
   REAL(SP) :: DTT(1:MTLOC,KBM1)    
   INTEGER :: I
   DATA     RHO /1.0E6/, CP /4.1796/

   DTT = 0.0

   DO I=1,MLOC
     DTT(I,1) = KT/(RHO*CP*D(I)*DZ(1))*(TE-T(I,1))
   ENDDO  

   RETURN
   END SUBROUTINE TEMPER

!************************************************************************
!**                 S U B R O U T I N E   S O L I D S                  **
!************************************************************************

   SUBROUTINE SOLIDS(DTSSI)
   
   USE MOD_WQM, ONLY :			&!
				SEDIMENT_CALC,  &!
                SSI,            &!
                V2,             &!
                WSS,            &!
                WSS1,            &!
                WSS2,            &!
                WSSHI,          &!
                WSSNET,WSSLIM,        &!
                Ksed, S_vmax,    &! B Clark April 20201 Ksed is half sat. conc for sediment sinking, S_vmax is maximum settling velocity
                SedHyp,SedOriginal !Added by B Clark April 2021 for sediment sinking velocity functionality

   USE MOD_LIMS, ONLY: MLOC, KBM1, MTLOC
   
   USE MOD_HYDROVARS, ONLY:  DZ	,D	
	
   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: DTSSI  

!******* Set special settling rate for high solids concs
   DTSSI = 0.0


!******* Settle particulate fraction

!Suspended sediment settling option for the velocity
   !first loop is originial with a baseline rate, and a high rate at concs > WSSLIM
   !second loop is a michaelis menten function that saturates velocity to the maximum sediment value based 
   !off of a half saturation coefficient (Kmax) which here is set to the mean observed sediment concentration
IF (SedOriginal) THEN
   DO I=1,MLOC
	  IF (SSI(I,1) > WSSLIM) THEN
		WSS=WSSHI 
	  ENDIF
     DTSSI(I,1)   = -WSS*SSI(I,1)/(D(I)*DZ(1)*86400.)
	 
   ENDDO  

   DO K=2,KBM1
	 DO I=1,MLOC
		IF (SSI(I,K) > WSSLIM) THEN
			WSS=WSSHI 
		ENDIF
       DTSSI(I,K)   = (WSS*SSI(I,K-1)-WSS*SSI(I,K))/(D(I)*DZ(K)*86400.)
     ENDDO  
   ENDDO


   ELSEIF(SedHyp) THEN
 !new exponential sinking velocity formulation based on 520 g m^-3 maximum concentration and the velocities from 
  !the above citation

  
  DO I=1,MLOC
    WSS1=S_vmax*SSI(I,1)/(Ksed+SSI(I,1))
     DTSSI(I,1)   = -WSS1*SSI(I,1)/(D(I)*DZ(1)*86400.)
   ENDDO  

   DO K=2,KBM1
    DO I=1,MLOC

       WSS1=S_vmax*SSI(I,K-1)/(Ksed+SSI(I,K-1))
       WSS2=S_vmax*SSI(I,K)/(Ksed+SSI(I,K))
       DTSSI(I,K)   = (WSS1*SSI(I,K-1)-WSS2*SSI(I,K))/(D(I)*DZ(K)*86400.)
    ENDDO  
   ENDDO



ENDIF
!******  Resuspension

! IF (SEDIMENT_CALC) THEN

  DO I = 1, MLOC	

    TAU_crit=0.005;
    Resusp_Mass=1.e-5;

    ! write(*,*)'Shear STress = ',ShearStress(i)
    ! write(*,*)'Resusp_mass = ',Resusp_mass
    ! write(*,*)'Tau Crit = ',TAU_crit

		IF(ShearStress(I) > TAU_crit) THEN
!		
			DTSSI(I,KBM1)=DTSSI(I,KBM1)+Resusp_mass*(ShearStress(I)-TAU_crit)/(D(I)*DZ(KBM1))! g m^-2 s^-1 Pa^-1 * Pa / m --> g m^-3 s^-1
!			
		ENDIF
  ENDDO
! ENDIF

   RETURN
   END SUBROUTINE SOLIDS


!************************************************************************
!**                  S U B R O U T I N E   Z O O P L                   **
!************************************************************************
#if defined(ZOOPLANKTON)
   SUBROUTINE ZOOPL(DTSZ,DTLZ)
   
   
   USE MOD_LIMS, ONLY : MLOC, KBM1, MTLOC

   USE MOD_WQM, ONLY:			&
                ANC1,           &!
                ANC2,           &!
                ANC3,           &!
                ASC1,           &!
                ASC2,           &!
                ASC3,           &!
                B1,             &!
                B2,             &!
                B3,             &!
                DOXG,           &!
                JDAY,           &!
            !    LDOC,           &!
                LPOC,           &!
                LPON,           &!
                LPOP,           &!
                Q1,             &!
                Q2,             &!
                Q3,             &!
              !  RDOC,           &!
                RPOC,           &!
                RPON,           &!
                RPOP,           &!
                T!,              &!
   USE WC_DOM, ONLY:         &
				CRRATE_SZ,CRRATE_LZ,  &
				NRRATE_SZ, NRRATE_LZ,  &
				PRRATE_SZ, PRRATE_LZ,   &
				FRD1,FRD2,  & !   
				 WC_NCDOC1, WC_CDOC1, &  !FRP
				WC_NCDOC2, WC_CDOC2
				
				
     
   IMPLICIT NONE  
   REAL(SP) :: DOREF,FRLP,FRRP,ANCPSZ,   &  !NRRATE, CRRATE,PRRATE ! B clark pull out C,N and P zooplankton rates to make allocatable for passing into mod_wc_doc
           ANCPLZ,APCPSZ,APCPLZ,FRB1SZ,FRB2SZ,FRB3SZ,FRB1LZ,  &  !FRD1,FRD2, NOW DEFINED IN WC_DOM
		   FRB2LZ,FRB3LZ
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: DTSZ, DTLZ
   INTEGER :: I,J,K,ITEMP

! DETERMINE JULIAN DAY
   DTSZ = 0.0
   DTLZ = 0.0
   
   J = 1.0 + AMOD(JDAY,365.25)

   DO K=1,KBM1
     DO I=1,MLOC

! COMPUTE PREY AVAILABLE TO ZOOPLANKTON

       B1ASZ(I,K) = MAX(0.,B1(I,K)-CTSZ)          
       B2ASZ(I,K) = MAX(0.,B2(I,K)-CTSZ)          
       B3ASZ(I,K) = MAX(0.,B3(I,K)-CTSZ)
	   
       DOC1ASZ(I,K) = MAX(0.,WC_NCDOC1(I,K)+WC_CDOC1(I,K)-CTSZ)  ! B Clark  This is the amount of DOC1 available to micro zooplankton        
       DOC2ASZ(I,K) = MAX(0.,WC_NCDOC2(I,K)+WC_CDOC2(I,K)-CTSZ)   ! B CLark this is amount of     DOC 2 available to micro zooplankton, micro zooplankton can't eat refractory DOC  
	   
       LPOCASZ(I,K) = MAX(0.,LPOC(I,K)-CTSZ)          
       RPOCASZ(I,K) = MAX(0.,RPOC(I,K)-CTSZ)
	   
       PRASZ(I,K) = UB1SZ*B1ASZ(I,K)+UB2SZ*B2ASZ(I,K)+UB3SZ*B3ASZ(I,K)      &
                + ULDSZ*DOC1ASZ(I,K)+URDSZ*DOC2ASZ(I,K)+ULPSZ*LPOCASZ(I,K)  &
                + URPSZ*RPOCASZ(I,K)        

       B1ALZ(I,K) = MAX(0.,B1(I,K)-CTLZ)          
       B2ALZ(I,K) = MAX(0.,B2(I,K)-CTLZ)          
       B3ALZ(I,K) = MAX(0.,B3(I,K)-CTLZ)
       SZALZ(I,K) = MAX(0.,SZ(I,K)-CTLZ)
       LPOCALZ(I,K) = MAX(0.,LPOC(I,K)-CTLZ)          
       RPOCALZ(I,K) = MAX(0.,RPOC(I,K)-CTLZ)
	   
       PRALZ(I,K) = UB1LZ*B1ALZ(I,K)+UB2LZ*B2ALZ(I,K)+UB3LZ*B3ALZ(I,K)      &
                + USZLZ*SZALZ(I,K)+ULLZ*LPOCALZ(I,K)+URLZ*RPOCALZ(I,K)

     ENDDO
   ENDDO  

! COMPUTE ZOOPLANKTON RATION AND METABOLISM

   DO K=1,KBM1
     DO I = 1,MLOC

! NUTRIENT LIMITATION

       CLSZ(I,K) = PRASZ(I,K)/(KHCSZ+PRASZ(I,K)+1.0E-30)
       CLLZ(I,K) = PRALZ(I,K)/(KHCLZ+PRALZ(I,K)+1.0E-30)

! RATION (GM PREY CARBON/GM ZOOPL CARBON/DAY)

       ITEMP = 10.*(T(I,K)+0.05)
       RSZ(I,K) = FTSZ(ITEMP)*CLSZ(I,K)*RMAXSZ(I,K)
       RLZ(I,K) = FTLZ(ITEMP)*CLLZ(I,K)*RMAXLZ(I,K)

! BASAL METABOLISM (1/DAY)

       BMSZ(I,K) = BMRSZ(I,K)*FTBMSZ(ITEMP)
       BMLZ(I,K) = BMRLZ(I,K)*FTBMLZ(ITEMP)

! LOSS TO PREDATION (GM C/M**3/DAY)

       PRSZLZ(I,K) = USZLZ*SZALZ(I,K)/(PRALZ(I,K)+1.0E-30)*RLZ(I,K)*LZ(I,K)
       PRSZ(I,K) = BPRSZ(I,K)*SZ(I,K)*SZ(I,K)*FTPRSZ(ITEMP)*TVPRSZ(J)
       PRLZ(I,K) = BPRLZ(I,K)*LZ(I,K)*LZ(I,K)*FTPRLZ(ITEMP)*TVPRLZ(J)
        
     ENDDO
   ENDDO  

! COMPUTE MORTALITY FACTORS

   DO K=1,KBM1
     DO I=1,MLOC

! DISSOLVED OXYGEN EFFECTS (1/DAY)

       DOREF = MIN(DOCRITSZ,DOXG(I,K))
       MSZ(I,K) = MZEROSZ*(1.-DOREF/DOCRITSZ)
       DOREF = MIN(DOCRITLZ,DOXG(I,K))
       MLZ(I,K) = MZEROLZ*(1.-DOREF/DOCRITLZ)

     ENDDO
   ENDDO  

! RATE OF CHANGE OF ZOOPLANKTON BIOMASS (GM C/M**3/DAY)

   DO K=1,KBM1
     DO I=1,MLOC

       GSZ(I,K) = ESZ(I,K)*(1-RFSZ(I,K))*RSZ(I,K)
       GLZ(I,K) = ELZ(I,K)*(1-RFLZ(I,K))*RLZ(I,K)
       DTSZ(I,K) = ((GSZ(I,K)-BMSZ(I,K)-MSZ(I,K))*SZ(I,K)-PRSZLZ(I,K)    &
                   -PRSZ(I,K))/86400.
       DTLZ(I,K) = ((GLZ(I,K)-BMLZ(I,K)-MLZ(I,K))*LZ(I,K)-PRLZ(I,K))/86400.

     ENDDO
   ENDDO  

! EFFECT OF ZOOPLANKTON ON ALGAE (GM C/M**3/DAY)

   DO K=1,KBM1
     DO I=1,MLOC

       B1SZ(I,K) = UB1SZ*B1ASZ(I,K)/(PRASZ(I,K)+1.0E-30)*RSZ(I,K)*SZ(I,K)
       B2SZ(I,K) = UB2SZ*B2ASZ(I,K)/(PRASZ(I,K)+1.0E-30)*RSZ(I,K)*SZ(I,K)
       B3SZ(I,K) = UB3SZ*B3ASZ(I,K)/(PRASZ(I,K)+1.0E-30)*RSZ(I,K)*SZ(I,K)

       B1LZ(I,K) = UB1LZ*B1ALZ(I,K)/(PRALZ(I,K)+1.0E-30)*RLZ(I,K)*LZ(I,K)
       B2LZ(I,K) = UB2LZ*B2ALZ(I,K)/(PRALZ(I,K)+1.0E-30)*RLZ(I,K)*LZ(I,K)
       B3LZ(I,K) = UB3LZ*B3ALZ(I,K)/(PRALZ(I,K)+1.0E-30)*RLZ(I,K)*LZ(I,K)

     ENDDO
   ENDDO  

! EFECT OF ZOOPLANKTON ON CARBON (GM C/M**3/DAY)

   DO K=1,KBM1   ! DOC_flag  all below is used to get the zooplankton doc contribution
     DO I=1,MLOC

       CRRATE_SZ(I,K) = ((1.-ESZ(I,K))*RSZ(I,K)+MSZ(I,K))*SZ(I,K)+PRSZ(I,K)  ! Need CRRATE for DOC calcs in mod_wc_doc, therefore expand to be an array, and make crrate_sz and crrate_lz
	   
       FRD1(I,K) = ULDSZ*DOC1ASZ(I,K)/(PRASZ(I,K)+1.0E-30)     ! B_Clark made CRRATE_SZ and LZ an array Sep 2015
       FRD2(I,K) = URDSZ*DOC2ASZ(I,K)/(PRASZ(I,K)+1.0E-30)     ! only small zooplankton consume DOC
	   
       FRLP = ULPSZ*LPOCASZ(I,K)/(PRASZ(I,K)+1.0E-30)
       FRRP = URPSZ*RPOCASZ(I,K)/(PRASZ(I,K)+1.0E-30)
       !LDOCSZ(I,K)  = -FRD1*RSZ(I,K)*SZ(I,K)+CRRATE*FLDOCSZ  ! DOC_flag  ! B Clark DOM migration, Sep 2015
       !RDOCSZ(I,K)  = -FRD2*RSZ(I,K)*SZ(I,K)+CRRATE*FRDOCSZ  ! DOC_flag  
       LPOCSZ(I,K) = -FRLP*RSZ(I,K)*SZ(I,K)+CRRATE_SZ(I,K)*FLPOCSZ
       RPOCSZ(I,K) = -FRRP*RSZ(I,K)*SZ(I,K)+CRRATE_SZ(I,K)*FRPOCSZ

       CRRATE_LZ(I,K) = ((1.-ELZ(I,K))*RLZ(I,K)+MLZ(I,K))*LZ(I,K)+PRLZ(I,K)
	   
       FRLP = ULLZ*LPOCALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       FRRP = URLZ*RPOCALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       !LDOCLZ(I,K)  = CRRATE*FLDOCLZ   ! DOC_flag  ! B Clark DOM migration, Sep 2015
       !RDOCLZ(I,K)  = CRRATE*FRDOCLZ   ! DOC_flag
       LPOCLZ(I,K) = -FRLP*RLZ(I,K)*LZ(I,K)+CRRATE_LZ(I,K)*FLPOCLZ  
       RPOCLZ(I,K) = -FRRP*RLZ(I,K)*LZ(I,K)+CRRATE_LZ(I,K)*FRPOCLZ

     ENDDO
   ENDDO  

! EFFECT OF ZOOPLANKTON ON NITROGEN (GM N/M**3/DAY)

   DO K=1,KBM1
     DO I=1,MLOC

       ANCPSZ = (ANC1*UB1SZ*B1ASZ(I,K)+ANC2*UB2SZ*B2ASZ(I,K)              &
                +ANC3*UB3SZ*B3ASZ(I,K)                                    &
                +ULPSZ*LPON(I,K)*LPOCASZ(I,K)/(LPOC(I,K)+1.0E-30)         & 
                +URPSZ*RPON(I,K)*RPOCASZ(I,K)/(RPOC(I,K)+1.0E-30))        &
                /(PRASZ(I,K)+1.0E-30)
       NRRATE_SZ(I,K) = ((ANCPSZ-ANCSZ*ESZ(I,K)*(1.-RFSZ(I,K)))*RSZ(I,K)          &
						+(BMSZ(I,K)+MSZ(I,K))*ANCSZ)*SZ(I,K)                      &
                        +PRSZ(I,K)*ANCSZ
       FRLP = ULPSZ*LPOCASZ(I,K)/(PRASZ(I,K)+1.0E-30)  ! ADDDED IN ARRAY SETUP FOR PASSING TO MOD_WC_DOM
       FRRP = URPSZ*RPOCASZ(I,K)/(PRASZ(I,K)+1.0E-30)
	   
       NH4SZ(I,K) =  NRRATE_SZ(I,K)*FNH4SZ
       !LDONSZ(I,K) = NRRATE*FLDONSZ   ! DON_flag   ! B Clark DOM migration Sep 2015, made NRRATE_SZ,LZ array 
       !RDONSZ(I,K) = NRRATE*FRDONSZ   ! DON_flag
       LPONSZ(I,K) = -FRLP*LPON(I,K)/(LPOC(I,K)+1.0E-30)                   &   !if negative, the Zooplankton are net-eating POC
                   *RSZ(I,K)*SZ(I,K)+NRRATE_SZ(I,K)*FLPONSZ         
       RPONSZ(I,K) = -FRRP*RPON(I,K)/(RPOC(I,K)+1.0E-30)                   &
                   *RSZ(I,K)*SZ(I,K)+NRRATE_SZ(I,K)*FRPONSZ         

       ANCPLZ = (ANC1*UB1LZ*B1ALZ(I,K)+ANC2*UB2LZ*B2ALZ(I,K)              &
                +ANC3*UB3LZ*B3ALZ(I,K)+ANCSZ*USZLZ*SZALZ(I,K)             &
                +ULLZ*LPON(I,K)*LPOCALZ(I,K)/(LPOC(I,K)+1.0E-30)          &
                +URLZ*RPON(I,K)*RPOCALZ(I,K)/(RPOC(I,K)+1.0E-30))         &
                /(PRALZ(I,K)+1.0E-30)
       NRRATE_LZ(I,K) = ((ANCPLZ-ANCLZ*ELZ(I,K)*(1.-RFLZ(I,K)))*RLZ(I,K)          &
                +(BMLZ(I,K)+MLZ(I,K))*ANCLZ)*LZ(I,K)                      &
                +PRLZ(I,K)*ANCLZ
       FRLP = ULLZ*LPOCALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       FRRP = URLZ*RPOCALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       NH4LZ(I,K) = NRRATE_LZ(I,K)*FNH4LZ
       !LDONLZ(I,K) = NRRATE*FLDONLZ     ! DON_flag     ! for some reason, there is no consumption term from Zooplankton for DON and DOP, maybe the assumption is that the zooplankton only create DOC.  
       !RDONLZ(I,K) = NRRATE*FRDONLZ     ! DON_flag  ! What about UREA  ???
       LPONLZ(I,K) = -FRLP*LPON(I,K)/(LPOC(I,K)+1.0E-30)*RLZ(I,K)*LZ(I,K)  &
                   +NRRATE_LZ(I,K)*FLPONLZ         
       RPONLZ(I,K) = -FRRP*RPON(I,K)/(RPOC(I,K)+1.0E-30)*RLZ(I,K)*LZ(I,K)  &
                   +NRRATE_LZ(I,K)*FRPONLZ

     ENDDO
   ENDDO  

! EFFECT OF ZOOPLANKTON ON PHOSPHORUS (GM P/M**3/DAY)

   DO K=1,KBM1
     DO I=1,MLOC

       APCPSZ = (Q1(I,K)*UB1SZ*B1ASZ(I,K)+Q2(I,K)*UB2SZ*B2ASZ(I,K)        &
                +Q3(I,K)*UB3SZ*B3ASZ(I,K)                                 &
                +ULPSZ*LPOP(I,K)*LPOCASZ(I,K)/(LPOC(I,K)+1.0E-30)         &
                +URPSZ*RPOP(I,K)*RPOCASZ(I,K)/(RPOC(I,K)+1.0E-30))        &
                /(PRASZ(I,K)+1.0E-30)
       PRRATE_SZ(I,K) = ((APCPSZ-APCSZ*ESZ(I,K)*(1.-RFSZ(I,K)))*RSZ(I,K)          &
						+(BMSZ(I,K)+MSZ(I,K))*APCSZ)*SZ(I,K)                      &
                       +PRSZ(I,K)*APCSZ
					   
       FRLP = ULPSZ*LPOCASZ(I,K)/(PRASZ(I,K)+1.0E-30)
       FRRP = URPSZ*RPOCASZ(I,K)/(PRASZ(I,K)+1.0E-30)
	   
       PO4SZ(I,K)  = PRRATE_LZ(I,K)*FPO4SZ
       !LDOPSZ(I,K) = PRRATE*FLDOPSZ       ! DOP_flag    ! B CLark DOM migration Sep 2015, made PRRATE_SZ, LZ array
       !RDOPSZ(I,K) = PRRATE*FRDOPSZ       ! DOP_flag  ! no loss term of phosphorus from zooplankton... is this correct? BC
       LPOPSZ(I,K) = -FRLP*LPOP(I,K)/(LPOC(I,K)+1.0E-30)                   &
                   *RSZ(I,K)*SZ(I,K)+PRRATE_SZ(I,K)*FLPOPSZ         
       RPOPSZ(I,K) = -FRRP*RPOP(I,K)/(RPOC(I,K)+1.0E-30)                   &
                   *RSZ(I,K)*SZ(I,K)+PRRATE_SZ(I,K)*FRPOPSZ 
       PIB1SZ(I,K) = Q1(I,K)*B1SZ(I,K)        
       PIB2SZ(I,K) = Q2(I,K)*B2SZ(I,K)        
       PIB3SZ(I,K) = Q3(I,K)*B3SZ(I,K)        

       APCPLZ = (Q1(I,K)*UB1LZ*B1ALZ(I,K)+Q2(I,K)*UB2LZ*B2ALZ(I,K)        &
                +Q3(I,K)*UB3LZ*B3ALZ(I,K)+USZLZ*APCSZ*SZALZ(I,K)          &
                +ULLZ*LPOP(I,K)*LPOCALZ(I,K)/(LPOC(I,K)+1.0E-30)          &
                +URLZ*RPOP(I,K)*RPOCALZ(I,K)/(RPOC(I,K)+1.0E-30))         &
                /(PRALZ(I,K)+1.0E-30)
       PRRATE_LZ(I,K) = ((APCPLZ-APCLZ*ELZ(I,K)*(1.-RFLZ(I,K)))*RLZ(I,K)          &
                +(BMLZ(I,K)+MLZ(I,K))*APCLZ)*LZ(I,K)                      &
                +PRLZ(I,K)*APCLZ
       FRLP = ULLZ*LPOCALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       FRRP = URLZ*RPOCALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       PO4LZ(I,K)  = PRRATE_LZ(I,K)*FPO4LZ
       !LDOPLZ(I,K) = PRRATE*FLDOPLZ    ! DOP_flag
       !RDOPLZ(I,K) = PRRATE*FRDOPLZ    ! DOP_flag
       LPOPLZ(I,K) = -FRLP*LPOP(I,K)/(LPOC(I,K)+1.0E-30)                   &
                   *RLZ(I,K)*LZ(I,K)+PRRATE_LZ(I,K)*FLPOPLZ         
       RPOPLZ(I,K) = -FRRP*RPOP(I,K)/(RPOC(I,K)+1.0E-30)                   &
                   *RLZ(I,K)*LZ(I,K)+PRRATE_LZ(I,K)*FRPOPLZ
       PIB1LZ(I,K) = Q1(I,K)*B1LZ(I,K)        
       PIB2LZ(I,K) = Q2(I,K)*B2LZ(I,K)        
       PIB3LZ(I,K) = Q3(I,K)*B3LZ(I,K)        

     ENDDO
   ENDDO  

! EFFECT OF ZOOPLANKTON ON DISSOLVED OXYGEN (GM DO/M**3/DAY)

   DO K=1,KBM1
     DO I=1,MLOC

       DOSZ(I,K) = (ESZ(I,K)*RFSZ(I,K)*RSZ(I,K)+BMSZ(I,K))*AOCRSZ*SZ(I,K)
       DOLZ(I,K) = (ELZ(I,K)*RFLZ(I,K)*RLZ(I,K)+BMLZ(I,K))*AOCRLZ*LZ(I,K)

     ENDDO
   ENDDO  

! EFFECT OF ZOOPLANKTON ON SILICA (GM SI/M**3/DAY)

   DO K=1,KBM1
     DO I=1,MLOC

       FRB1SZ = UB1SZ*B1ASZ(I,K)/(PRASZ(I,K)+1.0E-30)
       FRB2SZ = UB2SZ*B2ASZ(I,K)/(PRASZ(I,K)+1.0E-30)
       FRB3SZ = UB3SZ*B3ASZ(I,K)/(PRASZ(I,K)+1.0E-30)
       SASZ(I,K) = (FRB1SZ*ASC1+FRB2SZ*ASC2+FRB3SZ*ASC3)*RSZ(I,K)*SZ(I,K)*FRSASZ
       SUSZ(I,K) = (FRB1SZ*ASC1+FRB2SZ*ASC2+FRB3SZ*ASC3)*RSZ(I,K)*SZ(I,K)*(1.-FRSASZ)

       FRB1LZ = UB1LZ*B1ALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       FRB2LZ = UB2LZ*B2ALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       FRB3LZ = UB3LZ*B3ALZ(I,K)/(PRALZ(I,K)+1.0E-30)
       SALZ(I,K) = (FRB1LZ*ASC1+FRB2LZ*ASC2+FRB3LZ*ASC3)*RLZ(I,K)*LZ(I,K)*FRSALZ
       SULZ(I,K) = (FRB1LZ*ASC1+FRB2LZ*ASC2+FRB3LZ*ASC3)*RLZ(I,K)*LZ(I,K)*(1.-FRSALZ)
         
     ENDDO
   ENDDO  

   RETURN
   END SUBROUTINE ZOOPL
#endif         
        
!************************************************************************
!**                 S U B R O U T I N E   C A R B O N                  **
!************************************************************************

   SUBROUTINE PARTICULATE_CARBON(DTLPOC,DTRPOC)!DTLDOC,DTRDOC ! renamed to particulate carbon, all DOM is now in mod_wc_dom  
                                                                      ! B Clark Sep 2015 DOM migration														
   USE MOD_LIMS, ONLY : MLOC, KBM1, MTLOC
   
   USE MOD_WQM, ONLY :			&!
				!ALGDOC,         &!
                ALGPOC,         &!
                B1,             &!
                B2,             &!
                BM1,            &!
                BM2,            &!
                DOXG,           &!
                !FCLD1,          &! ! all moved to algdoc
                !FCLD2,          &!
                !FCLDP,          &!
                FCLP1,          &!
                FCLP2,          &!
                FCLPP,          &!
               !FCRD1,          &!
               ! FCRD2,          &! FCRDP,          &!
                
                FCRP1,          &!
                FCRP2,          &!
                FCRPP,          &!
                FDOP,           &!
                KDCALG,         &!
              !  KHR1,           &!
              !  KHR2,           &!
!                KLDC,           &!
                P1,             &!
                P2,             &!
                PR1,            &!
                PR2,            &!
                PRSP1,          &!
                PRSP2,          &!

				AANOX,          &!
                DENIT,          &!

                FTHDR,          &!
                FTMNL,          &!
                HDRLPOC,        &!
                HDRRPOC,        &!
                KHCOAG,         &!
                KHNDN,          &!
                KHODOC,         &!
                KLCALG,         &!
                KLPC,           &!
                ! KRCOAG,         &!
              !  KRDC,           &!
                KRPC,           &!
                KTHDR,          &!
                KTMNL,          &!
              !  LDOC,           &!
                LPOC,           &!

      !          MNLLDOC,        &!
       !         MNLRDOC,        &!
                NO3,            &!
                RATOX,          &!
              !  RDOC,           &!
                RPOC,           &!
                SALT,           &!
                T,              &!
                TRHDR,          &!
                TRMNL,          &!

				B,              &!
                BENDOC,         &!
                ! DLALGC,         &!
                
                SAV_CALC,       &!
                SAV_LOADS,      &!
                SEDIMENT_CALC,  &!
                V2,             &!
                WSL,            &!
                WSLNET,         &!
                WSR,            &!
                WSRNET!  &!,        &!
			!    SALTC

   USE MOD_HYDROVARS, ONLY:  DZ, D	
	
   USE MOD_SAV,ONLY: NSAVCELL, SAVCELL                                  &
                   ,LDOCSAVW, RDOCSAVW, LPOCSAVW, RPOCSAVW              &
                   ,LDOCEPIW, RDOCEPIW, LPOCEPIW, RPOCEPIW 
				   
   USE WC_DOM, ONLY : DOCBM1,DOCBM2,CP1C,CP2C,CP3C, &   !,DOCBM3
			                KDOC1, KDOC2, KDOC3,           &
							COAGC,                         &
							WC_NCDOC1, WC_CDOC1, ALGCAR,   &
							KHR1, KHR2,FCD1,FCD2!, KHR3

   IMPLICIT NONE
   REAL(SP) :: LPOC1, LPOC2, LPOC3, LPSETL, KLDOC, KLPOC    !, !LDOC1, LDOC2, LDOC3
   INTEGER :: I,K
   REAL(SP) ::  RPOC1, RPOC2, RPOC3 
   REAL(SP) ::  RPSETL,FTN
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: DTLPOC, DTRPOC!, ALGCAR!, COAG  !  DTLDOC, DTRDOC B Clark moved ALGCAR here to make array  

 !  DTLDOC = 0.0
 !  DTRDOC = 0.0
   DTLPOC = 0.0
   DTRPOC = 0.0

 !RGL modified to take out all algae 3 calculations
   
   DO K=1,KBM1
     DO I=1,MLOC

!********* Algal contribution

       DOCBM1(I,K) = KHR1/(KHR1+DOXG(I,K))    ! DOC_flag    ! hyperbolic saturation for DOC excretion from algae
       DOCBM2(I,K) = KHR2/(KHR2+DOXG(I,K))    ! DOC_flag

       CP1C(I,K) = (P1(I,K)*PRSP1+BM1(I,K))*B1(I,K)    ! DOC_flag   ! P1 = temperature dependent reaction velocity PRSP1
       CP2C(I,K) = (P2(I,K)*PRSP2+BM2(I,K))*B2(I,K)    ! DOC_flag
!
       LPOC1      = FCLP1*CP1C(I,k)+FCLPP*PR1(I,K)
       LPOC2      = FCLP2*CP2C(I,K)+FCLPP*PR2(I,K)

       RPOC1      = FCRP1*CP1C(I,K)+FCRPP*PR1(I,K)
       RPOC2      = FCRP2*CP2C(I,K)+FCRPP*PR2(I,K)

      ALGPOC(I,K)  = LPOC1+LPOC2+RPOC1+RPOC2
!********* Mineralization and hydrolysis

       ALGCAR(I,K)  = B1(I,K)+B2(I,K)!+B3(I,K)       !DOC_flag
	   
       !KLDOC       = KLDC(I,K)+KDCALG(I,K)*ALGCAR   !DOC_flag ! B Clark moved to wc_dom
       KLPOC       = KLPC+KLCALG*ALGCAR(I,K)   !DOC_flag
	   
       !SALTC       = MAX(0., SALT(I,K))          
       
	!   COAG(I,K)        = 0.5*(1.+TANH(SALTC-KHCOAG))*KRCOAG(I,K)*RDOC(I,K)   !DOC_flag
	   
       FTMNL(I,K)  = EXP(KTMNL*(T(I,K)-TRMNL))
       FTHDR(I,K)  = EXP(KTHDR*(T(I,K)-TRHDR))
       RATOX(I,K)  = (DOXG(I,K)+AANOX*KHODOC)/(KHODOC+DOXG(I,K))
	   
       DENIT(I,K)  = KDOC1(I,K)*FTMNL(I,K)*AANOX*KHODOC/(KHODOC+DOXG(I,K))        &    !DOC_flag
                    *NO3(I,K)/(KHNDN+NO3(I,K))*(WC_NCDOC1(I,K)+WC_CDOC1(I,K))

       HDRLPOC(I,K)= RATOX(I,K)*KLPOC*FTHDR(I,K)*LPOC(I,K)
       HDRRPOC(I,K)= RATOX(I,K)*KRPC*FTHDR(I,K)*RPOC(I,K)

 										 ! DOC_flag
	   DTLPOC(I,K) = (LPOC1+LPOC2-HDRLPOC(I,K))/86400.
	   DTRPOC(I,K) = (RPOC1+RPOC2-HDRRPOC(I,K)+COAGC(I,K))/86400. 

#if defined(ZOOPLANKTON)
	   DTLPOC(I,K) = DTLPOC(I,K)+(LPOCSZ(I,K)+LPOCLZ(I,K))/86400.  
	   DTRPOC(I,K)= DTRPOC(I,K)+(RPOCSZ(I,K)+RPOCLZ(I,K))/86400.
#endif
     ENDDO
   ENDDO  



   DO I=1,MLOC

!********* Settling

     LPSETL = -WSL*LPOC(I,1)/(D(I)*DZ(1))
     RPSETL = -WSR*RPOC(I,1)/(D(I)*DZ(1))
     DTLPOC(I,1) = DTLPOC(I,1)+LPSETL/86400.
     DTRPOC(I,1) = DTRPOC(I,1)+RPSETL/86400.

   ENDDO  

   DO K=2,KBM1
     DO I=1,MLOC

!********* Settling

       LPSETL = (WSL*LPOC(I,K-1)-WSL*LPOC(I,K))/(D(I)*DZ(K))
       RPSETL = (WSR*RPOC(I,K-1)-WSR*RPOC(I,K))/(D(I)*DZ(K))
       DTLPOC(I,K) = DTLPOC(I,K)+LPSETL/86400.
       DTRPOC(I,K) = DTRPOC(I,K)+RPSETL/86400.

!********* Settling flux


     ENDDO  
   ENDDO  
!Uncommented B CLark July 2020 for use in the Yukon River

   TAU_crit=0.005;
   Resusp_Mass=1e-5;
   DO I=1,MLOC    

   IF(ShearStress(I) > TAU_crit)  THEN ! Above commented bY B Clark June 2018 for below formulation that resuspends by critical shear stress
    DTLPOC(I,KBM1)=DTLPOC(I,KBM1)+Resusp_mass*(ShearStress(I)-TAU_crit)/(D(I)*DZ(KBM1))
    DTRPOC(I,KBM1)=DTRPOC(I,KBM1)+Resusp_mass*(ShearStress(I)-TAU_crit)/(D(I)*DZ(KBM1))
	
  ENDIF
   ENDDO
!******* Resuspension

   IF (SEDIMENT_CALC) THEN

     DO I=1,MLOC
       DTLPOC(I,KBM1) = DTLPOC(I,KBM1)+(WSL-WSLNET)   &
                        *LPOC(I,KBM1)/(D(I)*DZ(KBM1))/86400.
       DTRPOC(I,KBM1) = DTRPOC(I,KBM1)+(WSR-WSRNET)      &
                        *RPOC(I,KBM1)/(D(I)*DZ(KBM1))/86400.
     ENDDO           
   ENDIF


   END SUBROUTINE PARTICULATE_CARBON
   
!************************************************************************
!**                 S U B R O U T I N E   N I T R O G                  **
!************************************************************************

   SUBROUTINE NITROG(DTNH4,DTNO3,DTLPON,DTRPON)  ! ,DTLDON,DTRDON B CLark old DON
   
   USE MOD_LIMS, ONLY: MLOC, KBM1, MTLOC
   
   USE MOD_WQM, ONLY: 			&!
				ANC1,           &!
                ANC2,           &!
                B1,             &!
                B2,             &!
                BM1,            &!
                BM2,            &!
                DOXG,           &!

                FNI1,           &!
                FNIP,           &!
                KHN1,           &!
                KHN2,           &!
                KHNNT,          &!
                KHONT,          &!
                KTNT1,          &!
                KTNT2,          &!
                NFIX,           &!
                NH4,            &!
                NO3,            &!
                NT,             &!
                NTM,            &!
                P1,             &!
                P1NNF,          &!
                P2,             &!
                PN1,            &!
				NL1,            &!
				NL2,            &!
                PN2,            &!
                PRSP1,          &!
                PRSP2,          &!
                T,              &!
                TMNT,           &!

                ATMOS_LOADS,    &!
				
				ALGDON,         &!
                ALGNH4,         &!
                ALGNO3,         &!
                ALGPON,         &!
                B3,             &!
                
                FNI2,           &!
           !     FNLD1,          &! ! kinetics moved to mod_wc_dom B Clark Sep 2015
              !  FNLD2,          &!
              !  FNLDP,          &!
                FNLP1,          &!
                FNLP2,          &!
                FNLPP,          &!
            !    FNRD1,          &!
              !  FNRD2,          &!
              !  FNRDP,          &!
                FNRP1,          &!
                FNRP2,          &!
                FNRPP,          &!
                KDNALG,         &!
                KHCOAG,         &!
                KHNAVG,         &!
              !  KLDN,           &!
                KLNALG,         &!
                KLPN,           &!
                ! KRCOAG,         &!

                PR1,            &!
                PR2,            &!
             !   RDON,           &!
                SALT,           &!
				
				ANDC,           &!
                ATMNH4,         &!
                B,              &!
                BENDON,         &!
                BENNH4,         &!
                BENNO3,         &!
				
				AANOX,   &
				KHODOC,  &
				KHNDN,   &
				DOXG,   &
                DENIT,          &!
				                
			    DENITN,          &!! Added by B Clark to account of N and P under denitrification

                DENNO3,         &!

                FTHDR,          &!
                FTMNL,          &!
                HDRLPON,        &!
                HDRRPON,        &!
              !  KRDN,           &!
                KRPN,           &!
              !  LDON,           &!
                LPON,           &!
              !  MNLLDON,        &!
               ! MNLRDON,        &!
                PRECIP,         &!
                RATOX,          &!
                RPON,           &!
                SAV_CALC,       &!
                SEDIMENT_CALC,  &!
                V2,             &!
                WSL,            &!
                WSLNET,         &!
                WSR,            &!
                WSRNET,         &!
				
				ATMLDON,        &!
                ATMNO3,         &!
                ATMRDON!,        &!


   
   USE MOD_HYDROVARS, ONLY: DZ,D 
	
   USE MOD_SAV, ONLY: NSAVCELL, SAVCELL                                        &
                    ,NH4SAVW, NO3SAVW, LDONSAVW, RDONSAVW, LPONSAVW, RPONSAVW  &
                    ,NH4EPIW, NO3EPIW, LDONEPIW, RDONEPIW, LPONEPIW, RPONEPIW
   
   USE WC_DOM, ONLY : COAGN, &
					      NP1, NP2, NP3,  &
						  KDON1, KDON2, KDON3,   &
						  KDON1_IN, KDON2_IN, KDON3_IN, &
						  ALGCAR, MNLDON,  &   !
						  DPD30N,DPD20N,DPD10N, &  ! inorganic production of NH4 from photochemistry (I,K) dimensions
						  PHOTODEG, WC_CDON1, WC_NCDON1
   IMPLICIT NONE
   REAL(SP) ::   KLPON, LPON1, LPON2, LPON3,                         &!  KLDON,
              NH4A1, NH4A2,  NH4A3,  NO3A1,  NO3A2,  NO3A3,  LPSETL,    &
              NL1MOD ! NP1,  NP2,   NP3,   LDON1, LDON2, LDON3, !!! Now defined in mod_wc_dom B clark sep 2015
   REAL(SP) :: FTN,  RPON1, RPON2, RPON3    		! RDON1,RDON2, RDON3,
   REAL(SP):: DIN! B Clark added, Nh4+NO3
   REAL(SP) :: RPSETL! , SALTC!,  ALGCAR
   INTEGER :: I,K
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: 	DTNH4, 		&!gN/m^3/sec
										DTNO3, 		&!gN/m^3/sec
										!DTLDON, 	&!gN/m^3/sec     ! DON_flag
										!DTRDON, 	&!gN/m^3/sec      ! DON_flag
										DTLPON, 	&!gN/m^3/sec
										DTRPON       !gN/m^3/sec


! NITRIFICATION

   DO K=1,KBM1
     DO I=1,MLOC
        
!********* Temperature effect

       IF (T(I,K) < TMNT) THEN
         FTN = EXP(-KTNT1*(T(I,K)-TMNT)**2)
       ELSE
         FTN = EXP(-KTNT2*(TMNT-T(I,K))**2)
       ENDIF

!	          NT(I,K) = DOXG(I,K)/(KHONT+DOXG(I,K))*NH4(I,K)/(KHNNT+NH4(I,K))*FTN*NTM(I,K)
		! B Clark updated nitrification rate to include ammonium as a limiting reacting, now NTM is a first order rate rather than maximal concentration loss (weird logic before)
	    ! Now the function is a simple Michaelis mMeton curve with an ammonium and oxygen control, as well as Gaussian Temperature
       NT(I,K) = DOXG(I,K)/(KHONT+DOXG(I,K))*NH4(I,K)/(KHNNT+NH4(I,K))*FTN*NTM*NH4(I,K)
	   
     ENDDO
   ENDDO  

   DO K=1,KBM1
     DO I=1,MLOC
	 
	 DIN=NH4(I,K)+NO3(I,K)

!********* Algal nitrogen preference
       PN1(I,K) = NH4(I,K)*(NO3(I,K)/((KHN1+NH4(I,K))*(KHN1+NO3(I,K)))+KHN1  &
                /((1.E-30+NH4(I,K)+NO3(I,K))*(KHN1+NO3(I,K))))

	   
       PN2(I,K) = NH4(I,K)*(NO3(I,K)/((KHN2+NH4(I,K))*(KHN2+NO3(I,K)))+KHN2  &
                /((1.E-30+NH4(I,K)+NO3(I,K))*(KHN2+NO3(I,K))))
  
	   !NL1 added by B Clark May 2018 for consistency between mod_algal limitation and mod_kin uptake
       NL1(I,K) = ( NH4(I,K) + DIN - PN1(I,K)*DIN)   &
                     /(KHN1+NH4(I,K) + DIN - PN1(I,K)*DIN+1E-30)
					 
       NL2(I,K) = (NH4(I,K) + DIN - PN2(I,K)*DIN)   &
                     /(KHN2+NH4(I,K) + DIN - PN2(I,K)*DIN+1E-30)	   
	   
!********* Nitrogen Fixation

       NFIX(I,K) = (P1(I,K)-P1NNF(I,K))*ANC1*B1(I,K)

!********* Algal sources/sinks
	   
		! NP1 = loss of algal 1 N to different sinks
       NP1(I,K) = ANC1*(P1(I,K)*PRSP1+BM1(I,K))*B1(I,K) ! NP1 AND NP2 GET PASSED INTO MOD_WC_DOM FOR DON FORMULATIONS ! Algal p
       NP2(I,K) = ANC2*(P2(I,K)*PRSP2+BM2(I,K))*B2(I,K)          

! B Clark removed algal production of NH4 from net primary production, NP1 (shouldn't occur)	   
	   NH4A1 = -PN1(I,K)*P1(I,K)*ANC1*B1(I,K)+FNIP*PR1(I,K)*ANC1  !FNI1*NP1(I,K)
       NH4A2 = -PN2(I,K)*P2(I,K)*ANC2*B2(I,K)+FNIP*PR2(I,K)*ANC2 !FNI2*NP2(I,K)

       NO3A1      = (PN1(I,K)-1.)*P1(I,K)*ANC1*B1(I,K)	   
       NO3A2      = (PN2(I,K)-1.)*P2(I,K)*ANC2*B2(I,K)

       LPON1     = FNLP1*NP1(I,K)+FNLPP*PR1(I,K)*ANC1
       LPON2     = FNLP2*NP2(I,K)+FNLPP*PR2(I,K)*ANC2

       RPON1     = FNRP1*NP1(I,K)+FNRPP*PR1(I,K)*ANC1
       RPON2     = FNRP2*NP2(I,K)+FNRPP*PR2(I,K)*ANC2

       ALGNH4(I,K) = NH4A1+NH4A2!+NH4A3
       ALGNO3(I,K) = NO3A1+NO3A2!+NO3A3
     ALGPON(I,K) = LPON1+LPON2+RPON1+RPON2
!********* Mineralization and hydrolysis

	   ! CALCULATE REMINERALIZATION RATES. LEAVE THIS IN HERE BECAUSE OF DEPENDENCY ON NUTRIENTS
       KDON1(I,K)      = KDON1_IN + KDNALG*ALGCAR(I,K)*KHNAVG/(KHNAVG+NH4(I,K)+NO3(I,K)) ! DON_flag
	     KDON2(I,K)      = KDON2_IN + KDNALG*ALGCAR(I,K)*KHNAVG/(KHNAVG+NH4(I,K)+NO3(I,K)) ! DON_flag
	     KDON3(I,K)      = KDON3_IN + KDNALG*ALGCAR(I,K)*KHNAVG/(KHNAVG+NH4(I,K)+NO3(I,K)) ! DON_flag
	   
       KLPON      = KLPN+KLNALG*ALGCAR(I,K)*KHNAVG/(KHNAVG+NH4(I,K)+NO3(I,K)) ! DON_flag
	   
     !  SALTC       = MAX(0., SALT(I,K))          
   !    COAG(I,K)        = 0.5*(1.+TANH(SALTC-KHCOAG))*KRCOAG(I,K)*RDON(I,K)     ! DON_flag

       HDRLPON(I,K) = RATOX(I,K)*KLPON*FTHDR(I,K)*LPON(I,K)    ! DON_flag
       HDRRPON(I,K) = RATOX(I,K)*KRPN*FTHDR(I,K)*RPON(I,K)    ! DON_flag
       DENNO3(I,K)  = -ANDC*DENIT(I,K)    ! DON_flag

	     DENITN(I,K)  = KDON1(I,K)*FTMNL(I,K)*AANOX*KHODOC/(KHODOC+DOXG(I,K))        &    !Added by B CLark to account of DON loss and NH4 production from Denitrification
                    *NO3(I,K)/(KHNDN+NO3(I,K))*(WC_NCDON1(I,K)+WC_CDON1(I,K))

	     DTNH4(I,K) = (NH4A1+NH4A2+MNLDON(I,K)-NT(I,K) + DENITN(I,K) )/86400.
	   
       DTNO3(I,K) = (NT(I,K)-ANDC*DENIT(I,K)+NO3A1+NO3A2)/86400.
											 ! DON_flag
       DTLPON(I,K) = (LPON1+LPON2-HDRLPON(I,K))/86400.      
     
       DTRPON(I,K) = (RPON1+RPON2-HDRRPON(I,K) +COAGN(I,K))/86400.  
                     

#if defined(ZOOPLANKTON)
    DTNH4(I,K)=DTNH4(I,K)+ (NH4SZ(I,K)+NH4LZ(I,K))/86400.
	  DTLPON(I,K) =DTLPON(I,K) + (LPONSZ(I,K)+LPONLZ(I,K))/86400.      
  	DTRPON(I,K) =DTRPON(I,K) + (RPONSZ(I,K)+RPONLZ(I,K))/86400.
#endif
     ENDDO
   ENDDO  

   
 ! Photochemical addition of NH4 (and NO3??)  ! B Clark added April 2018
   	   
	   
#if defined(PHOTODEGRADATION)
	IF(PHOTODEG) THEN  
	
	DO I = 1,MLOC
	    DO K=1,KBM1
	
	      DTNH4(I,K)=DTNH4(I,K) + (DPD30N(I,K)+DPD20N(I,K)+DPD10N(I,K)) /(D(I)*DZ(K))   ! Addition of photochemical flux of NH4
												 !(mg N m^-2 s^-1) / layer thickness (m)  ----> mg N m^-3 s^-1
	    END DO
	 END DO  
    ENDIF  	   
#endif
   
   DO I=1,MLOC

!********* Settling

     LPSETL = -WSL*LPON(I,1)/(D(I)*DZ(1))
     RPSETL = -WSR*RPON(I,1)/(D(I)*DZ(1))
     DTLPON(I,1) = DTLPON(I,1)+LPSETL/86400.
     DTRPON(I,1) = DTRPON(I,1)+RPSETL/86400.

   ENDDO
   
   DO K=2,KBM1
     DO I=1,MLOC



!********* Settling

       LPSETL = (WSL*LPON(I,K-1)-WSL*LPON(I,K))/(D(I)*DZ(K))
       RPSETL = (WSR*RPON(I,K-1)-WSR*RPON(I,K))/(D(I)*DZ(K))
       DTLPON(I,K) = DTLPON(I,K)+LPSETL/86400.
       DTRPON(I,K) = DTRPON(I,K)+RPSETL/86400.

     ENDDO
   ENDDO  

!******* Benthic fluxes  

   DO I=1,MLOC
     DTNH4(I,KBM1)  = DTNH4(I,KBM1)+BENNH4(I)/(D(I)*DZ(KBM1))/86400.		!gN/m^3/sec
     DTNO3(I,KBM1)  = DTNO3(I,KBM1)+BENNO3(I)/(D(I)*DZ(KBM1))/86400.		!This is the JNO3 going to watercolumn bottom layer of nitrate 
	 
																	!JNO3 =S(NO31-NO30) has unit gN/m^2/day  which is flux of Nitrate to water column = cocentration * m/d 
																	!If you can get JLDOC in gC/m^2/day you can add this ti DTLDOC equation, same for DTRDOC
	 
   ENDDO  

!******* Resuspension
   !uncommented B Clark, July 2020 for use in Yukon River
       
   TAU_crit=0.005;
   Resusp_Mass=0.0000;
	 	IF(ShearStress(I) > TAU_crit) THEN ! Above commented bY B Clark June 2018 for below formulation that resuspends by critical shear stress
			DTLPON(I,KBM1)=DTLPON(I,KBM1)+Resusp_mass*(ShearStress(I)-TAU_crit)/(D(I)*DZ(KBM1)) ! scaled by factor of 10 from POC and iss
			DTRPON(I,KBM1)=DTRPON(I,KBM1)+Resusp_mass*(ShearStress(I)-TAU_crit)/(D(I)*DZ(KBM1))			
		ENDIF
   IF (SEDIMENT_CALC) THEN

     DO I=1,MLOC
       DTLPON(I,KBM1) = DTLPON(I,KBM1)+(WSL-WSLNET)       &
                        *LPON(I,KBM1)/(D(I)*DZ(KBM1))/86400.
       DTRPON(I,KBM1) = DTRPON(I,KBM1)+(WSR-WSRNET)       &
                        *RPON(I,KBM1)/(D(I)*DZ(KBM1))/86400.
	 
     ENDDO
   ENDIF

!******* SAV and epiphytes

   IF (SAV_CALC) THEN
     DO I=1,NSAVCELL
       B=SAVCELL(I)
       DTNH4(B,KBM1)   = DTNH4(B,KBM1) +( NH4SAVW(B)+ NH4EPIW(B))/(D(B)*DZ(KBM1))/86400.
       DTNO3(B,KBM1)   = DTNO3(B,KBM1) +( NO3SAVW(B)+ NO3EPIW(B))/(D(B)*DZ(KBM1))/86400.
       !DTLDON(B,KBM1)  = DTLDON(B,KBM1)+(LDONSAVW(B)+LDONEPIW(B))/(D(B)*DZ(KBM1))/86400.
       !DTRDON(B,KBM1)  = DTRDON(B,KBM1)+(RDONSAVW(B)+RDONEPIW(B))/(D(B)*DZ(KBM1))/86400.
       DTLPON(B,KBM1)  = DTLPON(B,KBM1)+(LPONSAVW(B)+LPONEPIW(B))/(D(B)*DZ(KBM1))/86400.
       DTRPON(B,KBM1)  = DTRPON(B,KBM1)+(RPONSAVW(B)+RPONEPIW(B))/(D(B)*DZ(KBM1))/86400.
     ENDDO
   ENDIF

!******* Atmospheric Loads
 
   IF(ATMOS_LOADS) THEN ! Kurt Glaesemann 13 April 2015       ! DON_flag
    DO I=1,MLOC  											  ! DON_flag
        DTNH4(I,1) = DTNH4(I,1) +PRECIP*ATMNH4/(D(I)*DZ(1))    
        DTNO3(I,1) = DTNO3(I,1) +PRECIP*ATMNO3/(D(I)*DZ(1))
        !DTLDON(I,1)= DTLDON(I,1)+PRECIP*ATMLDON/(D(I)*DZ(1))     ! DON_flag
        !DTRDON(I,1)= DTRDON(I,1)+PRECIP*ATMRDON/(D(I)*DZ(1))      ! DON_flag
    ENDDO
   ENDIF

 !  write(*,*)'DTNH4 = ',DTNH4
   !IF(SED_DOM)  THEN  ! B Clark SEDIMENT DOM update  July 2015  ! MOVED TO MOD_WC_DOM SEP 2015
   ! DO I = 1,MLOC					! DON_flag
   !	 DTLDON(I,KBM1) = DTLDON(I,KBM1)+J_LDON_SHARE(I)/(D(I)*DZ(KBM1))/86400.         ! DON_flag
!	 DTRDON(I,KBM1) = DTRDON(I,KBM1)+J_RDON_SHARE(I)/(D(I)*DZ(KBM1))/86400.          ! DON_flag
!	ENDDO  ! DON_flag
   !ENDIF   ! DON_flag
!	   
   RETURN
   END SUBROUTINE NITROG

!************************************************************************
!**                 S U B R O U T I N E   P H O S P H                  **
!************************************************************************
#if defined(PHOSPHORUS)
   SUBROUTINE PHOSPH(DTPO4,DTLPOP,DTRPOP,DTPIP,           &   ! DOP_flag ,DTLDOP,DTRDOP B Clark old DOP
                     DTPIB1,DTPIB2,DTPIB3)

   USE MOD_LIMS, ONLY: MLOC, KBM1, MTLOC
   
   USE MOD_WQM, ONLY : 			&!
				B1,             &!
                B2,             &!
                BM1,            &!
                BM2,            &!
                FPI1,           &!
                FPI2,           &!
                FPIP,           &!
            !    FPLD1,          &!
             !   FPLD2,          &!
             !   FPLDP,          &!
                FPLP1,          &!
                FPLPP,          &!
                !FPRD1,          &!
               ! FPRD2,          &!
             !   FPRDP,          &!
                KHP1,           &!
                KHP2,           &!
                P1,             &!
                P2,             &!
                PO4,            &!
                PR1,            &!
                PR2,            &!
                PRSP1,          &!
                PRSP2,          &!
                Q01,            &!
                Q02,            &!
                Q1,             &!
                Q2,             &!
                VMAX1,          &!
                VMAX2,          &!

				ALGDOP,         &!
                ALGPO4,         &!
                ALGPOP,         &!
                B3,             &!
				
				
				AANOX,   &
				KHODOC,  &
				KHNDN,   &
				DOXG,   &
               DENITP,          &! ! B Clark for denitirfication phosphorus, April 2018
				NO3,   &
				
                ATMOS_LOADS,    &!

                FPLP2,          &!
                FPRP1,          &!
                FPRP2,          &!
                FPRPP,          &!
                FTHDR,          &!
                FTMNL,          &!
                HDRLPOP,        &!
                HDRRPOP,        &!
                KADPO4,         &!
                KDPALG,         &!
                KHPAVG,         &!
              !  KLDP,           &!
                KLPALG,         &!
                KLPP,           &!
              !  KRDP,           &!
                KRPP,           &!
              !  LDOP,           &!
                LPOP,           &!
              !  MNLLDOP,        &!
              !  MNLRDOP,        &!
                RATOX,          &!
              !  RDOP,           &!
                RPOP,           &!
                SSI,            &!
                WSL,            &!
                WSS,            &!
                WSS1,           &!
                WSS2,           &!
				        ATMLDOP,        &!
                ATMPO4,         &!
                ATMRDOP,        &!
                B,              &!
                BENDOP,         &!
                BENPO4,         &!
                PIB1,           &!
                PIB2,           &!
                PIB3,           &!
                PIP,            &!
                PIP_CALC,       &!
                PRECIP,         &!
                SAV_CALC,       &!
                SEDIMENT_CALC,  &!
                SOLIDS_CALC,    &!
                ! V2,             &!
                WS1,            &!
                WS1NET,         &!
                WS2,            &!
                WS2NET,         &!
                WS3,            &!
                WSLNET,         &!
                WSR,            &!
                WSRNET,         &!
                WSSHI,          &!
                WSSNET!,         &!

   
    USE MOD_HYDROVARS, ONLY:  &
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
		!,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
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
		 ,D	 !&			!CURRENT DEPTH 
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

   USE MOD_SAV,ONLY: NSAVCELL, SAVCELL                                  &
                   ,PO4SAVW, LDOPSAVW, RDOPSAVW, LPOPSAVW, RPOPSAVW     &
                   ,PO4EPIW, LDOPEPIW, RDOPEPIW, LPOPEPIW, RPOPEPIW

   USE  WC_DOM, ONLY : ALGCAR, PP1, PP2,    &
                       KDOP1,KDOP2,KDOP3,   &
                       KDOP1_IN,KDOP2_IN,KDOP3_IN,   &
					   MNLDOP, WC_CDOP1, WC_NCDOP1
					   
					   
  
   IMPLICIT NONE
   
   !Wen Long added PO41, PO42
   
   REAL(SP) PO41, PO42
   
   
   REAL(SP) ::  KLPOP, LPOP1, LPOP2, LPOP3, LPSETL !LDOP1, LDOP2, LDOP3  ! KLDOP
   REAL(SP) :: DF, PO4AVL, PUP1, PUP2, PUP3,  PO43,     &  ! PP1, PP2, PP3, RDOP1, RDOP2, RDOP3,
            RPOP1, RPOP2, RPOP3
   REAL(SP) ::PFU, PPU, PFD, PPD, PO4SET, RPSETL, PIB1SETL,    &   ! ALGCAR, 
           PIB2SETL, PIB3SETL, PF, PP 
   INTEGER :: I, J, K
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: DTPO4,  DTLPOP, DTRPOP,  &  ! DTLDOP, DTRDOP,
        DTPIP, DTPIB1, DTPIB2, DTPIB3

   DTPO4 = 0.0
   !DTLDOP = 0.0
   !DTRDOP = 0.0
   DTLPOP = 0.0
   DTRPOP = 0.0
   DTPIP = 0.0
   DTPIB1 = 0.0
   DTPIB2 = 0.0
   DTPIB3 = 0.0

   DO K=1,KBM1
   
     DO I=1,MLOC

!********* Algae sources/sinks

!********* Uptake via Monod formulation

       DF     = 1. !1./(1.+KADPO4*SSI(I,K))
       PO4AVL = MAX(0.,DF*PO4(I,K))
       PUP1   = VMAX1*B1(I,K)*PO4AVL/(KHP1+PO4AVL+1.0E-6)
       PUP2   = VMAX2*B2(I,K)*PO4AVL/(KHP2+PO4AVL+1.0E-6)
!       PUP3   = VMAX3*B3(I,K)*PO4AVL/(KHP3+PO4AVL+1.0E-6)
          
       PP1(I,K)  = Q01*(P1(I,K)*PRSP1+BM1(I,K))*B1(I,K)
       PP2(I,K)  = Q02*(P2(I,K)*PRSP2+BM2(I,K))*B2(I,K)
!       PP3  = Q03*(P3(I,K)*PRSP3+BM3(I,K))*B3(I,K)
       PO41 = FPI1*PP1(I,K)-PUP1+FPIP*PR1(I,K)*Q1(I,K)
       PO42 = FPI2*PP2(I,K)-PUP2+FPIP*PR2(I,K)*Q2(I,K)

       LPOP1     = FPLP1*PP1(I,K)+FPLPP*PR1(I,K)*Q1(I,K)
       LPOP2     = FPLP2*PP2(I,K)+FPLPP*PR2(I,K)*Q2(I,K)

       RPOP1     = FPRP1*PP1(I,K)+FPRPP*PR1(I,K)*Q1(I,K)
       RPOP2     = FPRP2*PP2(I,K)+FPRPP*PR2(I,K)*Q2(I,K)
  
       ALGPO4(I,K) = PO41+PO42!+PO43
	   ALGPOP(I,K) = LPOP1+LPOP2+RPOP1+RPOP2

!********* Mineralization and hydrolysis

 
       DF      = 1. !1./(1.+KADPO4*SSI(I,K))
       PO4AVL  = MAX(1.E-6,DF*PO4(I,K))
       KDOP1(I,K)   = KDOP1_IN +KDPALG*ALGCAR(I,K)*KHPAVG/(KHPAVG+PO4AVL)
	   KDOP2(I,K)   = KDOP2_IN +KDPALG*ALGCAR(I,K)*KHPAVG/(KHPAVG+PO4AVL)
	   KDOP3(I,K)   = KDOP3_IN +KDPALG*ALGCAR(I,K)*KHPAVG/(KHPAVG+PO4AVL)
	   
       KLPOP   = KLPP+KLPALG*ALGCAR(I,K)*KHPAVG/(KHPAVG+PO4AVL)           

       HDRLPOP(I,K) = RATOX(I,K)*KLPOP*FTHDR(I,K)*LPOP(I,K)
       HDRRPOP(I,K) = RATOX(I,K)*KRPP*FTHDR(I,K)*RPOP(I,K)
	   
	   DENITP(I,K)  = KDOP1(I,K)*FTMNL(I,K)*AANOX*KHODOC/(KHODOC+DOXG(I,K))        &    !Added by B CLark to account of DOP loss and PO4 production from Denitrification
                    *NO3(I,K)/(KHNDN+NO3(I,K))*(WC_NCDOP1(I,K)+WC_CDOP1(I,K))
					   
	   
       DTPO4(I,K) = (PO41+PO42+MNLDOP(I,K)      &   !+MNLRDOP(I,K)
                    +PO4SZ(I,K)+PO4LZ(I,K) + DENITP(I,K))/86400.
         
       DTLPOP(I,K) = (LPOP1+LPOP2-HDRLPOP(I,K)+LPOPSZ(I,K)+LPOPLZ(I,K))/86400.
       DTRPOP(I,K) = (RPOP1+RPOP2-HDRRPOP(I,K)+RPOPSZ(I,K)+RPOPLZ(I,K))/86400.
       DTPIB1(I,K) = (PUP1-PP1(I,K)-PR1(I,K)*Q1(I,K)-PIB1SZ(I,K)-PIB1LZ(I,K))/86400.
       DTPIB2(I,K) = (PUP2-PP2(I,K)-PR2(I,K)*Q2(I,K)-PIB2SZ(I,K)-PIB2LZ(I,K))/86400.

     ENDDO
   ENDDO  

   DO I=1,MLOC

!********* Settling

     PFD    = KADPO4*SSI(I,1)/(1.+KADPO4*SSI(I,1))
     PPD    = PFD*PO4(I,1)
     PO4SET = -WSS*PPD/(D(I)*DZ(1))
     LPSETL = -WSL*LPOP(I,1)/(D(I)*DZ(1))
     RPSETL = -WSR*RPOP(I,1)/(D(I)*DZ(1))
     PIB1SETL = -WS1*PIB1(I,1)/(D(I)*DZ(1))
     PIB2SETL = -WS2*PIB2(I,1)/(D(I)*DZ(1))
!     PIB3SETL = -WS3(I,1)*PIB3(I,1)/(D(I)*DZ(1))
     DTPO4(I,1)  = DTPO4(I,1)+PO4SET/86400.
     DTLPOP(I,1) = DTLPOP(I,1)+LPSETL/86400.
     DTRPOP(I,1) = DTRPOP(I,1)+RPSETL/86400.
     DTPIB1(I,1) = DTPIB1(I,1)+PIB1SETL/86400.
     DTPIB2(I,1) = DTPIB2(I,1)+PIB2SETL/86400.
!     DTPIB3(I,1) = DTPIB3(I,1)+PIB3SETL/86400.

   ENDDO

   DO K=2,KBM1
     DO I=1,MLOC

!********* Settling below top depth

       PFU    = KADPO4*SSI(I,K-1)/(1.+KADPO4*SSI(I,K-1))  !Particle fraction of P in upperward grid (K-1)
       PPU    = PFU*PO4(I,K-1)                            !Particulate Inorganic P (PO4 in particles) (gP/m^3) at K-1
       PFD    = KADPO4*SSI(I,K)/(1.+KADPO4*SSI(I,K))      !Particle fraction of P in downward grid(K)
       PPD    = PFD*PO4(I,K)                              !Particulate Inorganic P (PO4 in particles) (gP/m^3) at K
       
       PO4SET = (WSS*PPU-WSS*PPD)/(D(I)*DZ(K))
       
       LPSETL = (WSL*LPOP(I,K-1)-WSL*LPOP(I,K))/(D(I)*DZ(K))
       RPSETL = (WSR*RPOP(I,K-1)-WSR*RPOP(I,K))/(D(I)*DZ(K))
       PIB1SETL = (WS1*PIB1(I,K-1)-WS1*PIB1(I,K))/(D(I)*DZ(K))
       PIB2SETL = (WS2*PIB2(I,K-1)-WS2*PIB2(I,K))/(D(I)*DZ(K))
       PIB3SETL = (WS3*PIB3(I,K-1)-WS3*PIB3(I,K))/(D(I)*DZ(K))
       DTPO4(I,K)  = DTPO4(I,K)+PO4SET/86400.
       DTLPOP(I,K) = DTLPOP(I,K)+LPSETL/86400.
       DTRPOP(I,K) = DTRPOP(I,K)+RPSETL/86400.
       DTPIB1(I,K) = DTPIB1(I,K)+PIB1SETL/86400.
       DTPIB2(I,K) = DTPIB2(I,K)+PIB2SETL/86400.
       DTPIB3(I,K) = DTPIB3(I,K)+PIB3SETL/86400.

     
     ENDDO
   ENDDO  

!******* Benthic fluxes

   DO I=1,MLOC
     DTPO4(I,KBM1)  = DTPO4(I,KBM1)+BENPO4(I)/(D(I)*DZ(KBM1))/86400.
     !DTLDOP(I,KBM1) = DTLDOP(I,KBM1)+BENDOP(I)/(D(I)*DZ(KBM1))/86400.
   ENDDO

!******* Resuspension
       
   TAU_crit=0.005;
   Resusp_Mass=0.00000;
						
		IF(ShearStress(I) > TAU_crit) THEN ! Above commented bY B Clark June 2018 for below formulation that resuspends by critical shear stress
			DTLPOP(I,KBM1)=DTLPOP(I,KBM1)+Resusp_mass*(ShearStress(I)-TAU_crit)/(D(I)*DZ(KBM1)) ! scaled from Carbon by a factor of 100
			DTRPOP(I,KBM1)=DTRPOP(I,KBM1)+Resusp_mass*(ShearStress(I)-TAU_crit)/(D(I)*DZ(KBM1))			
    ENDIF
    
   IF (SEDIMENT_CALC) THEN

     DO I=1,MLOC
       PF = KADPO4*SSI(I,KBM1)/(1.+KADPO4*SSI(I,KBM1))
       PP = PF*PO4(I,KBM1)
       DTPO4(I,KBM1)  = DTPO4(I,KBM1)+(WSS-WSSNET)*PP     &
                        /(D(I)*DZ(KBM1))/86400.
       DTLPOP(I,KBM1) = DTLPOP(I,KBM1)+(WSL-WSLNET)       &
                        *LPOP(I,KBM1)/(D(I)*DZ(KBM1))/86400.
       DTRPOP(I,KBM1) = DTRPOP(I,KBM1)+(WSR-WSRNET)       &
                        *RPOP(I,KBM1)/(D(I)*DZ(KBM1))/86400.

       DTPIB1(I,KBM1) = DTPIB1(I,KBM1)+(WS1-WS1NET)*      &
                        PIB1(I,KBM1)/(D(I)*DZ(KBM1))/86400.
       DTPIB2(I,KBM1) = DTPIB2(I,KBM1)+(WS2-WS2NET)*      &
                        PIB2(I,KBM1)/(D(I)*DZ(KBM1))/86400.

     ENDDO
   ENDIF

!******* SAV and epiphytes

     IF (SAV_CALC) THEN
       DO I=1,NSAVCELL
         B=SAVCELL(I)
          ! 
         !Wen Long: There is a potential problem here, if 
         !          SAV's height is larger than the thickness of the bottom layer
         !          we should actually distribute the flux into being added multiple layers
         !          in the water column
         !
         DTPO4(B,KBM1)   = DTPO4(B,KBM1)+(PO4SAVW(B)+PO4EPIW(B))                  &
                      /(D(B)*DZ(KBM1))/86400.

         DTLPOP(B,KBM1)  = DTLPOP(B,KBM1)+(LPOPSAVW(B)+LPOPEPIW(B))               &
                      /(D(B)*DZ(KBM1))/86400.
         DTRPOP(B,KBM1)  = DTRPOP(B,KBM1)+(RPOPSAVW(B)+RPOPEPIW(B))               &
                      /(D(B)*DZ(KBM1))/86400.
       ENDDO
     ENDIF

!******* Atmospheric loads

     IF(ATMOS_LOADS) THEN ! Kurt Glaesemann 13 April 2015
        DO I=1,MLOC
            DTPO4(I,1) = DTPO4(I,1)+PRECIP*ATMPO4/(D(I)*DZ(1))
        ENDDO  
    ENDIF

!******* Particulate Inorganic Phosphorus

   IF (SOLIDS_CALC .AND. PIP_CALC) THEN

!******  Settling
     DO I=1,MLOC
       DTPIP(I,1)   = WSS*PIP(I,1)/(D(I)*DZ(1)*86400.)
     ENDDO
 
     DO K=2,KBM1
       DO I=1,MLOC
         DTPIP(I,K)   = (WSS*PIP(I,K-1)-WSS*PIP(I,K))       &
                      /(D(I)*DZ(K)*86400.)
       ENDDO
     ENDDO  

!******  Resuspension

     IF (SEDIMENT_CALC) THEN
       DO I=1,MLOC
         DTPIP(I,KBM1) = DTPIP(I,KBM1)+(WSS-WSSNET)    &
                         *PIP(I,KBM1)/(D(I)*DZ(KBM1))/86400.
						 
						 
       ENDDO
     ENDIF
          
	 
	 ! SEDIMENT DOP diffusion
	 ! B  Clark sediment DOM update July 2015
   !IF(SED_DOM)  THEN   ! DOP_flag
   ! DO I = 1,MLOC     ! DOP_flag
   !	 DTLDOP(I,KBM1) = DTLDOP(I,KBM1)+J_LDOP_SHARE(I)/(D(I)*DZ(KBM1))/86400. ! DOP_flag
!	 DTRDOP(I,KBM1) = DTRDOP(I,KBM1)+J_RDOP_SHARE(I)/(D(I)*DZ(KBM1))/86400.  ! DOP_flag
!	ENDDO   ! DOP_flag
   !ENDIF     ! DOP_flag
!	 
	 
   ENDIF
        
   RETURN
   END SUBROUTINE PHOSPH
      
#endif

!************************************************************************
!**                S U B R O U T I N E   C O D M N D                   **
!************************************************************************

   SUBROUTINE CODMND(DTCOD_out)
   
   USE MOD_LIMS, ONLY: MLOC, MTLOC, KBM1
   
   USE MOD_WQM, ONLY: 			&!
				COD,            &!
                DOXG,           &!

                FTCOD,          &!
                KCOD,           &!
                KHOCOD,         &!
                KTCOD,          &!
                T,              &!
                TRCOD, BENCOD!,          &!

   USE MOD_HYDROVARS, ONLY:  &
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
		!,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
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
		 ,D	 !&			!CURRENT DEPTH 
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

     
   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP) :: DTCOD_out(1:MTLOC,KBM1)
   
   DTCOD_out  = 0.0

!******* Change in chemical oxygen demand

   DO K=1,KBM1
     DO I=1,MLOC
!! Analytical solution
!       FTCOD(I,K) = KCOD(I,K)!*EXP(-KTCOD)
!RGL need to check below
!write(*,*)'FTCOD = ',FTCOD(I,K)       
        FTCOD(I,K) = KCOD*EXP(KTCOD*(T(I,K)-TRCOD))
        DTCOD_out(I,K) = -(DOXG(I,K)/(KHOCOD+DOXG(I,K))*FTCOD(I,K)*COD(I,K))/86400.
!       DTCOD(I,K) = -FTCOD(I,K)*COD(I,K)/86400.
     ENDDO
   ENDDO  

!!******* Sediment demand
!! why is sediment not included?
!   
!   !oxidation of COD is negative (loss) while Benthic flux ia a gain into the bottom layer, loss term is also
!   ! added into the loss of O2 in oxygen subroutine below
!   ! B Clark uncommented BENCOD TERM, March 2018
!   DO I=1,MLOC
!     DTCOD_out(I,KBM1) = DTCOD_out(I,KBM1)+BENCOD(I)/(D(I)*DZ(KBM1))/86400.  !Wen Long: Wrong not to include
!                                                                             !benthic COD contribution!!!
!                                                                            !Original code has BENCOD term here
!    IF(NGID(I)==352)THEN
!          write(*,*)'BENCOD in mod_kin = ',BENCOD(I)
!    ENDIF


!   ENDDO

   RETURN
   END SUBROUTINE CODMND

!************************************************************************
!**                 S U B R O U T I N E   O X Y G E N                  **
!************************************************************************

   SUBROUTINE OXYGEN(DTDO)

   USE MOD_LIMS, ONLY : MLOC, MTLOC, KBM1
   
   USE MOD_WQM, ONLY: 			&!
				AOCR,           &!
                AONT,           &!
                B1,             &!
                B2,             &!
                BM1,            &!
                BM2,            &!
                DCOD,           &!
                DOPR,           &!
                DOXG,           &!

              !  FCLD1,          &!
              !  FCLD2,          &!
                FCLP1,          &!
                FCLP2,          &!
              !  FCRD1,          &!
              !  FCRD2,          &!
                FCRP1,          &!
                FCRP2,          &!
                FDOP,           &!
              !  KHR1,           &!
              !  KHR2,           &!
                NITRIF,         &!
                NT,             &!
                P1,             &!
                P2,             &!
                PN1,            &!
                PN2,            &!
                PR1,            &!
                PR2,            &!
                PRSP1,          &!
                PRSP2,          &!

				AREAR,          &!
                B,              &!
                BENDO,          &!
                BREAR,          &!
                COD,            &!
                CREAR,          &!
                DDOC,           &!
                ! DORALG,         &!
                DOS,            &!
                FTCOD,          &!
                KHOCOD,         &!
     !           MNLLDOC,                &!
         !       MNLRDOC,                &!
                REAERDO,                &!
                RESP,           &!
                SALT,           &!
                SAV_CALC,               &!
                SAV_LOADS,              &!
                T,              &!
                WMS!,            &!

   
   USE MOD_HYDROVARS, ONLY:  &
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
		!,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
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
		 
   USE WC_DOM, ONLY : MNLDOC,  &
                      KHR1, KHR2  ,  &
					  FCD1, FCD2
   USE MOD_SAV, ONLY: NSAVCELL, SAVCELL         &
                    ,DOSAVW,DOEPIW              &
                    ,SEDDOSAV
                    
   USE MOD_WQMINIT, ONLY : COD_CALC
                    
   IMPLICIT NONE
   REAL(SP) :: KRDO
   INTEGER :: I,K
   REAL(SP) :: FRDO1, FRDO2, FRDO3, CP1, CP2, CP3, ALGDO, DOR1, DOR2, DOR3,  &
           DOP1, DOP2, DOP3
   REAL(SP) :: FACTOR, CHLORN, TDOS, RNU 
   REAL(SP) :: DTDO(1:MTLOC,KBM1)
   
   DTDO   = 0.0

   DO K=1,KBM1
     DO I=1,MLOC

!********* Nitrification

       NITRIF = AONT*NT(I,K)

!********* Uptake/respiration

	   ! B CLARK UPDATED TO ONLY USE THE FRACTION THAT IS LOST TO DOC, PARTICULATE AND DISSOLVED
	   ! April 2018
	   FRDO1  = (1.-FCD1-FCLP1-FCRP1)!*DOXG(I,K)/(KHR1+DOXG(I,K))  
       FRDO2  = (1.-FCD2-FCLP2-FCRP2)!*DOXG(I,K)/(KHR2+DOXG(I,K))   
	   
       CP1    = P1(I,K)*PRSP1+BM1(I,K)
       CP2    = P2(I,K)*PRSP2+BM2(I,K)

       ALGDO = (FRDO1*CP1*B1(I,K)+FRDO2*CP2*B2(I,K))*AOCR
       DOR1   = ((1.3-0.3*PN1(I,K))*P1(I,K)-FRDO1*CP1)*AOCR*B1(I,K)
       DOR2   = ((1.3-0.3*PN2(I,K))*P2(I,K)-FRDO2*CP2)*AOCR*B2(I,K)          

	   ! B Clark updated so just use a straing DIC production rate (FDOP) from predation (PR1) and an O:C ratio 
	   ! For predation April 2018
	     DOP1   = FDOP*PR1(I,K)*AOCR!*DOXG(I,K)/(KHR1+DOXG(I,K))    
       DOP2   = FDOP*PR2(I,K)*AOCR!*DOXG(I,K)/(KHR2+DOXG(I,K))
	   
       DOPR = DOP1+DOP2!+DOP3
       
       DDOC  = AOCR*(MNLDOC(I,K)) ! CHANGE IN DISSOLVED OXYGEN DUE TO DOC RESPIRATION, MNLDOC is from mod_wc_dom

!********* Change in dissolved oxygen

!       DTDO(I,K)  = (DOR1+DOR2+DOR3-DOP1-DOP2-DOP3-DDOC(I,K)-DCOD(I,K)       &
!                  -NITRIF(I,K)-DOSZ(I,K)-DOLZ(I,K))/86400.
        DTDO(I,K) = (DOR1+DOR2-DOP1-DOP2-DDOC -NITRIF)/86400.    
        
        ! RESP(I,K)  =(ALGDO + DOPR+DDOC+NITRIF)/86400.
		
IF(COD_CALC)THEN
	   
            DCOD  = DOXG(I,K)/(KHOCOD+DOXG(I,K))*FTCOD(I,K)*COD(I,K)
			
            DTDO(I,K) = DTDO(I,K) - DCOD/86400.
		  
            ! RESP(I,K) = RESP(I,K) +DCOD
			
ENDIF
       
#if defined(ZOOPLANKTON)
               DTDO(I,K)=DTDO(I,K)-(DOSZ(I,K)-DOLZ(I,K))/86400.
#endif
    
        ENDDO
   ENDDO  

	
!******* Reaeration

   FACTOR = AREAR*(BREAR*WMS)**CREAR
   DO I=1,MLOC
     CHLORN    = SALT(I,1)/1.80655
     TDOS      = T(I,1)
     RNU       = 0.54 + 0.7*TDOS/30 - 0.07*SALT(I,1)/35.
!     KRDO      = 0.157*RNU*(1.5*WMS)**1.5
     KRDO      = FACTOR*RNU
     DOS       = 14.5532+TDOS*(0.0054258*TDOS-0.38217)-CHLORN     &
                    *(0.1665+TDOS*(9.796E-5*TDOS-5.866E-3))
     DTDO(I,1) = DTDO(I,1)+KRDO/(D(I)*DZ(1))*(DOS-DOXG(I,1))/86400.

     !Wen Long debugging areation
     REAERDO(i,1) = KRDO/(D(I)*DZ(1))*(DOS-DOXG(I,1))  !gO2/m^3/day rate of areation in surface layer
                                                     !Stumm and Morgan (third edition) pp242
   ENDDO

!******* Sediment oxygen demand

   DO I=1,MLOC
     DTDO(I,KBM1) = DTDO(I,KBM1)+BENDO(I)/(D(I)*DZ(KBM1))/86400.  !gO2/m^3/sec
   ENDDO  
								!BENDO needs to be gO2/m2/day
   
   
   RETURN
   END SUBROUTINE OXYGEN

!************************************************************************
!**                  S U B R O U T I N E   S I L I C A                 **
!************************************************************************

   SUBROUTINE SILICA(DTSIUPB,DTSIAT)
   
   USE MOD_LIMS, ONLY: MLOC, MTLOC, KBM1
   
   USE MOD_WQM, ONLY: 			&!
				ALGRES,         &!
                ALGUP,          &!
                ASC1,           &!
                ASC2,           &!
                ASC3,           &!
                B1,             &!
                B2,             &!
                B3,             &!
                BM1,            &!
                BM2,            &!
                BM3,            &!
                
                FSAP,           &!
                KADSA,          &!
                KSUA,           &!
                KTSUA,          &!
                P1,             &!
                P2,             &!
                P3,             &!
                PR1,            &!
                PR2,            &!
                PR3,            &!
                PRSP1,          &!
                PRSP2,          &!
                PRSP3,          &!
                PSD,            &!
                SIUPB,           &!
                SAP,            &!
                T,              &!
                TRSUA,          &!

				BENSA,          &!
                SIAT,            &!
                SEDIMENT_CALC,  &!
                SSI,            &!
                V2,             &!
                WSS,            &!
                WSSNET,         &!
                WSU,            &!
                WSUNET!,         &!

   
   USE MOD_HYDROVARS, ONLY:  DZ	,D	 
	
     
   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP) :: SUP, SP1, SP2, SP3, PFU, PF, SASETL,SUSETL 
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: 	DTSIUPB,		&  !rate of change of particulate biogenic silica (unavaiable, gSi/m^3/sec)
									DTSIAT	!rate of change of avaiable total silica (avaiable, gSi/m^3/sec)


   DTSIUPB = 0.0
   DTSIAT = 0.0


   DO K=1,KBM1
     DO I=1,MLOC

!********* Particulate silica dissolution

       PSD = KSUA*SIUPB(I,K)*EXP(KTSUA*(T(I,K)-TRSUA))
 
!********* Predation

       SAP = FSAP*(PR1(I,K)*ASC1+PR2(I,K)*ASC2+PR3(I,K)*ASC3)
       SUP = (1.-FSAP)*(PR1(I,K)*ASC1+PR2(I,K)*ASC2+PR3(I,K)*ASC3)

!********* Algal uptake and respiration

       ALGUP = P1(I,K)*ASC1*B1(I,K)+P2(I,K)*ASC2*B2(I,K)+P3(I,K)*ASC3*B3(I,K)
	   
       SP1 = ASC1*(P1(I,K)*PRSP1+BM1(I,K))*B1(I,K)
       SP2 = ASC2*(P2(I,K)*PRSP2+BM2(I,K))*B2(I,K)
    !    SP3 = ASC3*(P3(I,K)*PRSP3+BM3(I,K))*B3(I,K)
       ALGRES = SP1+SP2!+SP3

!********* Change in silica

       DTSIUPB(I,K) = (SUP-PSD+ALGRES)/86400.
       DTSIAT(I,K)  = (SAP+PSD-ALGUP)/86400.

! #if defined(ZOOPLANKTON)	
!        DTSIUPUB(I,K)=DTSIUPUB(I,K)+(SUSZ(I,K)+SULZ(I,K))/86400.
!        DTSIAT(I,K)=DTSIAT(I,K) + (SASZ(I,K)+SALZ(I,K))/86400.
! #endif
     ENDDO
   ENDDO  

   DO I=1,MLOC

!********* Settling

     PF     = KADSA*SSI(I,1)/(1.+KADSA*SSI(I,1))
     SASETL = -WSS*PF*SIAT(I,1)/(D(I)*DZ(1))
     SUSETL = -WSU*SIUPB(I,1)/(D(I)*DZ(1))

!********* Change in silica

     DTSIUPB(I,1) = DTSIUPB(I,1)+SUSETL/86400.
     DTSIAT(I,1)  = DTSIAT(I,1) +SASETL/86400.



   ENDDO

   DO K=2,KBM1
     DO I=1,MLOC

!********* Settling (source term of layer k due to settling)

       PFU    = KADSA*SSI(I,K-1)/(1.+KADSA*SSI(I,K-1))
       PF     = KADSA*SSI(I,K)/(1.+KADSA*SSI(I,K))
       SASETL = (WSS*PFU*SIAT(I,K-1)-WSS*PF*SIAT(I,K))/(D(I)*DZ(K))
       SUSETL = (WSU*SIUPB(I,K-1)-WSU*SIUPB(I,K))/(D(I)*DZ(K))

!********* Change in silica

       DTSIUPB(I,K) = DTSIUPB(I,K)+SUSETL/86400.
       DTSIAT(I,K) = DTSIAT(I,K)+SASETL/86400.


     ENDDO
   ENDDO  

!******* Benthic fluxes

   DO I=1,MLOC
     DTSIAT(I,KBM1) = DTSIAT(I,KBM1)+BENSA(I)/(D(I)*DZ(KBM1))/86400.
   ENDDO

!******* Resuspension

   IF (SEDIMENT_CALC) THEN

     DO I=1,MLOC
       PF           = KADSA*SSI(I,KBM1)/(1.+KADSA*SSI(I,KBM1))
       DTSIAT(I,KBM1) = DTSIAT(I,KBM1)+(WSS-WSSNET)*PF        &
                      *SIAT(I,KBM1)/(D(I)*DZ(KBM1))/86400.
       DTSIUPB(I,KBM1) = DTSIUPB(I,KBM1)+(WSU-WSUNET)           &
                      *SIUPB(I,KBM1)/(D(I)*DZ(KBM1))/86400.
     ENDDO         
   ENDIF

   RETURN
   END SUBROUTINE SILICA

!************************************************************************
!**               S U B R O U T I N E   B E N C O M P                  **
!************************************************************************

   SUBROUTINE BEN_FLUX  !simple bottom boundary flux calcualtin (parallel to sediment flux model)
						!modulate the input benthic flux to give some temperature control
   
   USE MOD_LIMS, ONLY : MLOC, MTLOC, KBM1
   
   USE MOD_WQM, ONLY : 			&
				BENCOD,         &!
                BENCODB,        &!
                BENDO,          &!
                BENDOB,         &!
                BENNH4,         &!
                BENNH4B,        &!
                BENNO3,         &!
                BENNO3B,        &!
                BENPO4,         &!
                BENPO4B,        &!
                BENSA,          &!
                DOXG,           &!

                KHSO,           &!
                KSDOC,          &!
                KSNH4,          &!
                KSNO3,          &!
                KSO,            &!
                KSPO4,          &!
                KSSA,           &!
                MTCNO3,         &!
                NO3,            &!
                SEDNO3,         &!
                T,              &!
                TRSDOC,         &!
                TRSNH4,         &!
                TRSNO3,         &!
                TRSO,           &!
                TRSPO4,         &!
                TRSSA,          &!

				B1,             &!
                B2,             &!
                BENDOC,         &!
                BENDOCB,        &!
                BENSAB,         &!

                LPOC,           &!
                LPON,           &!
                LPOP,           &!
                PCFWS,          &!
                PNFWS,          &!
                PPFWS,          &!
                PSFWS,          &!
				SSFWS,			&!
                SIUPB,           &!
                Q1,             &!
                Q2,             &!
                RPOC,           &!
                RPON,           &!
                RPOP,           &!
                WS1,            &!
                WS2,            &!
                WSL,            &!
                WSR,            &!
                WSU,            &!
				
				ANC1,           &!
                ANC2,           &!
                ASC1,           &!
                ASC2,           &!
                SIAT,            &!
                KADPO4,         &!
                KADSA,          &!
                PO4,            &!
                SSI,            &!
                WSS!,            &!

   
   USE MOD_HYDROVARS, ONLY: DZ ,D	
	
     
   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP) :: FTSOD, FTNH4, FTNO3, FTPO4, FTSA, FTDOC 
   REAL(SP) :: PF, PPO4

   DO I=1,MLOC

!******* Calculate temperature effects

     FTSOD    = EXP(KSO*(T(I,KBM1)-TRSO))
     FTNH4    = EXP(KSNH4*(T(I,KBM1)-TRSNH4))
     FTNO3    = EXP(KSNO3*(T(I,KBM1)-TRSNO3))
     FTPO4    = EXP(KSPO4*(T(I,KBM1)-TRSPO4))
     FTSA     = EXP(KSSA*(T(I,KBM1)-TRSSA))
     FTDOC    = EXP(KSDOC*(T(I,KBM1)-TRSDOC))

!********* Dissolved oxygen effect on sediment oxygen demand

     BENDO(I) = BENDOB(I)*FTSOD*DOXG(I,KBM1)/(KHSO+DOXG(I,KBM1))
     BENCOD(I) = BENCODB(I)-BENDOB(I)*FTSOD*KHSO/(KHSO+DOXG(I,KBM1)) !WLong, here we should have BENCODB on RHS instead of BENDOB ??

!********* Balance of nutrients (g/m2/day)

     BENNH4(I) = BENNH4B(I)*FTNH4
     BENNO3(I) = BENNO3B(I)+FTNO3*MTCNO3*(SEDNO3-NO3(I,KBM1))
     BENPO4(I) = BENPO4B(I)*FTPO4
     BENSA(I)  = BENSAB(I)*FTSA
     BENDOC(I) = BENDOCB(I)*FTDOC

   ENDDO

!******* Compute particle flux to sediments for mass-balance purposes

   DO I=1,MLOC

!******* First zero flux accumlators

     PCFWS(I) = 0.
     PNFWS(I) = 0.
     PPFWS(I) = 0.
     PSFWS(I) = 0.
	 SSFWS(I) = 0.

!******* Now accumulate fluxes of labile and refractory particles as well as suspended solids

	 !!!!!WLong: Seems these were already calculated in mod_sed.F!!!!!
     PPFWS(I) = PPFWS(I)-WSL*LPOP(I,KBM1)-WSR*RPOP(I,KBM1)
     PNFWS(I) = PNFWS(I)-WSL*LPON(I,KBM1)-WSR*RPON(I,KBM1)
     PCFWS(I) = PCFWS(I)-WSL*LPOC(I,KBM1)-WSR*RPOC(I,KBM1)
     PSFWS(I) = PSFWS(I)-WSU*SIUPB(I,KBM1)
     SSFWS(I) = SSFWS(I)-WSS*SSI(I,KBM1)
	 
!******* Now accumulate fluxes of algal biomass

				
	!RGL fixed bug  in PCFWS - twp WS2's instead of WS2 and WS3
     PCFWS(I) = PCFWS(I)-WS1*B1(I,KBM1)                          &
                        -WS2*B2(I,KBM1)!-WS3(I,KBM1)*B3(I,KBM1)				

! below per tykim replaced ANC3 with ANC2 in second equation
				
     PNFWS(I) = PNFWS(I)-WS1*ANC1*B1(I,KBM1)                     &
                        -WS2*ANC2*B2(I,KBM1)
!				        -WS3(I,KBM1)*ANC3*B3(I,KBM1)
				
     PPFWS(I) = PPFWS(I)-Q1(I,KBM1)*WS1*B1(I,KBM1)               &
                        -Q2(I,KBM1)*WS2*B2(I,KBM1)
!                       -Q3(I,KBM1)*WS3*B3(I,KBM1)

				
     PSFWS(I) = PSFWS(I)-WS1*ASC1*B1(I,KBM1)                     &
                        -WS2*ASC2*B2(I,KBM1)
!                       -WS3*ASC3*B3(I,KBM1)

!******* Accumulate adsorbed phosphate and silica flux (things that move with suspended solids****

     PF           = KADPO4*SSI(I,KBM1)/(1.+KADPO4*SSI(I,KBM1))
     PPO4         = PF*PO4(I,KBM1)
     PPFWS(I)     = PPFWS(I)-WSS*PPO4
     PF           = KADSA*SSI(I,KBM1)/(1.+KADSA*SSI(I,KBM1))
     PSFWS(I)     = PSFWS(I)-WSS*PF*SIAT(I,KBM1)

   ENDDO

   RETURN
   END SUBROUTINE BEN_FLUX

END MODULE MOD_KIN
