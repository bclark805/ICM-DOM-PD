!************************************************************************
!** This is a new version of algal input and kinetics.  Created for    **
!** The Tributary Refinements portion of the Chesapeake Bay study.     **
!**                                                                    **
!** Refinements include creation of an algal input routine, a new      **
!** P vs. I curve, a silica fraction for all algal groups, option to   **
!** spatially vary inputs, lookup table for temperature functions,     **
!** switch to uE from langleys, option to specify minimum or product   **
!** function for growth limitations.                                   **
!**                                                                    **
!**                  Algal  Subroutines for CE-QUAL-ICM                **
!**                                                                    **
!**                  Chesapeake Bay Tributary Refinements              **
!**                         February 7, 1996                           **
!**                                                                    **
!**                    Water Quality Modeling Group                    **
!**                    U.S. Army Corps of Engineers                    **
!**                    Waterways Experiment Station                    **
!**                    Vicksburg, Mississippi 39180                    **
!**                                                                    **
!************************************************************************

   MODULE MOD_ALGAL
   
	   USE MOD_PREC, ONLY: SP
	   USE MOD_SIZES, ONLY : MGL
	   USE MOD_LIMS, ONLY: KBM1, MTLOC, MLOC
	   
	   USE MOD_WQM, ONLY :		&!				
                ACHLC1,         &!
                ACHLC2,         &!
                ACHLC3,         &!
                ALOFN,          &!
                ALPHMIN1,       &!
                ALPHMIN2,       &!
                ALPHMIN3,       &!
                ALPHRAT1,       &!
                ALPHRAT2,       &!
                ALPHRAT3,       &!
                ANC1,           &!
                ANC2,           &!
                ANC3,           &!
                ASC1,           &!
                ASC2,           &!
                ASC3,           &!
                ASRAT,          &!
                B1,             &!
                B2,             &!
                B2GR,           &!
                B3,             &!
                BCHLC1,         &!
                BCHLC2,         &!
                BCHLC3,         &!
                BM1,            &!
                BM2,            &!
                BMR1,           &!
                BMR2,           &!
                BMR3,           &!
                BPR1,           &!
                BPR2,           &!
                BPR3,           &!
                CCHL1,          &!
                CCHL2,          &!
                CCHL3,          &!
                CCHLC1,         &!
                CCHLC2,         &!
                CCHLC3,         &!
                CFIX,           &!
                CHLCMN1,        &!
                CHLCMN2,        &!
                CHLCMN3,        &!
                !FCLD1,          &!now FCD1
                !FCLD2,          &!now FCD2
                !FCLD3,          &!now FCD3
               ! FCLDP,          &!!now FCDP
                FCLP1,          &!
                FCLP2,          &!
                FCLP3,          &!
                FCLPP,          &!
                !FCRD1,          &!now FCD1
                !FCRD2,          &!!now FCD2
                !FCRD3,          &!!now FCD3
                !FCRDP,          &!now FCDP
                FCRP1,          &!
                FCRP2,          &!
                FCRP3,          &!
                FCRPP,          &!
                FDOP,           &!
                FI1,            &!
                FI2,            &!
                FNFIX,          &!
                FNI1,           &!
                FNI2,           &!
                FNI3,           &!
                FNIP,           &!
                !FNLD1,          &!
                !FNLD2,          &!
                !FNLD3,          &!
             !   FNLDP,          &! ! now FNDP
                FNLP1,          &!
                FNLP2,          &!
                FNLP3,          &!
                FNLPP,          &!
                !FNRD1,          &! Now FND1
                !FNRD2,          &! Now FND2
                !FNRD3,          &!Now FND3
               ! FNRDP,          &!now FNDP
                FNRP1,          &!
                FNRP2,          &!
                FNRP3,          &!
                FNRPP,          &!
                FNUP,           &!  !Not used at all
                FPI1,           &!
                FPI2,           &!
                FPI3,           &!
                FPIP,           &!
                !FPLD1,          &! now FPD1
                !FPLD2,          &!now FPD2
                !FPLD3,          &!now FPD3
                !FPLDP,          &! now FPDP
                FPLP1,          &!
                FPLP2,          &!
                FPLP3,          &!
                FPLPP,          &!
                !FPRD1,          &!!now FPD1
                !FPRD2,          &!!now FPD2
                !FPRD3,          &!!now FPD3
               ! FPRDP,          &!! now FPDP
                FPRP1,          &!
                FPRP2,          &!
                FPRP3,          &!
                FPRPP,          &!
                FSAP,           &!
                !FT1,            &!
                !FT2,            &!
                !FT3,            &!
                !FTBM1,          &!
                !FTBM2,          &!
                !FTBM3,          &!
                !FTPR,           &!
                GPP,            &!
                IK1,            &!
                IK2,            &!
                IT,             &!
                JDAY,           &!
                KHN1,           &!
                KHN2,           &!
                KHN3,           &!
                KHNAVG,         &!
                KHNFIX,         &!
                KHP1,           &!
                KHP2,           &!
                KHP3,           &!
                KHPAVG,         &!
                !KHR1,           &!
                !KHR2,           &!
                !KHR3,           &!
                KHS1,           &!
                KHS2,           &!
                KHS3,           &!
                KHST1,          &!
                KHST2,          &!
                KHST3,          &!
                MINPROD,        &!
                NH4,            &!
                NL1,            &!
                NL2,            &!
                NO3,            &!
                NPP,            &!
                P1,             &!
                P1NNF,          &!
                P2,             &!
                PL1,            &!
                PL2,            &!
                PM1,            &!
                PM2,            &!
                PM3,            &!
                PN1,            &!
                PN2,            &!
                PO4,            &!
                PR1,            &!
                PR2,            &!
                PRSP1,          &!
                PRSP2,          &!
                PRSP3,          &!
                Q01,            &!
                Q02,            &!
                Q03,            &!
                SEDIMENT_CALC,  &!
                SNFIX,          &!
                STF1,           &!
                STF2,           &!
                STF3,           &!
                T,              &!
                TOTAL_NETPP,    &!
                TVPR,           &!
                V2,             &!
                VMAX1,          &!
                VMAX2,          &!
                VMAX3,          &!
                WS1,            &!
                WS1NET,         &!
                WS2,            &!
                WS2NET,         &!
                WS3,            &!
                WS3NET
				
       USE MOD_FILEINFO, ONLY : 	&!
				!,DIA			&!
				!,CBC 			&!
				!,S1			&!
				!,S2			&!
				!,S3			& !                 
				!,BFI			&!
	            !,BAI           &!
				!,MET			&!
				!,BFO			& !
				!,KEI			&!
				!,ATM			&!
				!,STL			& !
				AGR			& !
				!,SVI			& !
				!,SVO			& !
				!,KFL			& !
				!,ZOO			& !
				!,ZFO			& !
				,ALO      	!	&!
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

	!Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
		USE MOD_CONTROL, ONLY : 	&
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


   USE MOD_OWQ, ONLY:   				&
						!INITKE,			&	!
						!INTKE,			&	!
						!NXKEI,			&	!
						!LAT, 			&	!
						!LONG,			&	!
						!GSTAR440,		&	!
						!PHISTAR676,		&	!
						!BNSTARNTU,		&	!
						!BNSTARCHL,      &	!
						!PSTARINT,		&	!
						!PSTAR440,		&	!
						!PSTARCOLOR,		&	!
						!SCDOM,			&	!
						!EZEROINT,    	&	!
						!TCHL2CHL,		&	!
						!DOFFSET,		&	!
						!DEPTHKE,		&	!
						!NWAVEL,			&	!
						!COLOR, 			&	!
						!TURB,			&	!
						!OPZ,		 	&	!
						!G1,				&	!
						!G2,				&	!
						!ALAMB, 			&	!
						!GLAMB,			&	!
						!PLAMB,			&	!
						!PHILAMB,        &	!
                        !BLAMB,			&	!
						!WAVEL,			&	!
						!EZERO,			&	!
						!GMUZERO,		&	!
						IAVG,			&	!
						KESS!,			&	!
						!IATBOT,		&	!
						!OWQ_ALLOC,		&	!
						!OWQ_DEALLOC,	&	!
						!OWQ_READ!,		&	!
						!LGHT_ATTN

	USE MOD_ZOOP, ONLY:			&!
					 B1SZ,    	&!
					 B2SZ,    	&!
					 !B3SZ,    	&!
					 B1LZ,    	&!
					 B2LZ!,    	&!
					 !B3LZ,    	&!
					 !DOSZ,    	&!
					 !DOLZ,    	&!
					 !SASZ,    	&!
					 !SUSZ,    	&!
					 !SALZ,    	&!
					 !SULZ,    	&!
					 !ACLSZ,	&!
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
					 !FTLZ,		&!
					 !FTSZ,		&!
					 !FTBMSZ,	&!
					 !FTBMLZ,	&!
					 !FTPRSZ,	&!
					 !FTPRLZ,	&!
					 !SZ,		&!
					 !LZ				
						 
   IMPLICIT NONE
   
   !LB: moved the definitions below to this location to have acces to TRPR, TMP1, TR1, TMP2, TR2, TMP3, TR3 withing get_f*** functions
   REAL(SP) :: KTG11, KTG12, KTB1, KTG21, KTG22, KTB2, KTG31
   REAL(SP) :: KTG32, KTB3,  KTPR
   REAL(SP) :: TRPR, TMP1, TR1, TMP2, TR2, TMP3, TR3 
   
		!
		!   KURT GLAESEMANN 8 SEPT 2009 MOVE ALL VARIABLES OUT OF GLOBAL (HERE)
		!   AND DOWN TO WHERE THEY BELONG
		!
   
   CONTAINS
   
	!Subroutines:
	! subroutine	ALG_READ()
	! subroutine	ALGAE()
    !
   
!************************************************************************
!**             S U B R O U T I N E   A L G _ R E A D                  **
!************************************************************************

   SUBROUTINE ALG_READ
#  if defined (MULTIPROCESSOR)
	USE MOD_PAR, ONLY:   &	  
			!EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			!ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			!NLID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			!ELID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			!NLID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,	&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			!EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
			!NHE,        &   !			!NUMBER OF HALO ELEMENTS
			!HE_LST,		&	!(:)       	!!GLOBAL IDENTITIES OF HALO ELEMENTS
			!HE_OWN,		&	!(:)       	!!OWNER OF HALO ELEMENTS
			!NBN ,       &   !          	!!NUMBER OF BOUNDARY NODES
			!MX_MLT,     &   !          	!!MAX MULTIPLICITY OF BOUNDARY NODES
			!BN_LST, 	&	!	   		!!GLOBAL IDENTITY OF BOUNDARY NODES
			!BN_LOC, 	&	!			!!LOCAL IDENTITY OF BOUNDARY NODES
			!BN_MLT, 	&	!			!!MULTIPLICITY OF BOUNDARY NODES
			!BN_NEY, 	&	!			!!NODE OWNER LIST
			!NDE_ID, 	&	!			!! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
			NHN,  		&	!			!!NUMBER OF HALO NODES
			HN_LST!, 	&	!			!!LIST OF HALO NODES 
			!HN_OWN, 	&	!			!!PRIMARY OWNER OF HALO NODES
			!EC,			&	!
			!NC,			&	!
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
   use wc_dom, only: sanity_check,sant, FND1, FND2, FPD1, FPD2, FCD1, FCD2   
	USE MOD_WQM, only : FNDP, FCDP, FPDP
	
   IMPLICIT NONE
   real(sp) :: sanity
   INTEGER :: I, J, F, K
   REAL(SP), ALLOCATABLE :: RTMP21(:,:),RTMP22(:,:),RTMP23(:,:)
   CHARACTER(LEN=72) :: ALG_TITLE(6)
   CHARACTER(LEN=8)  :: SPVAR1, SPVAR2, SPVAR3
   CHARACTER(LEN=8)  :: PRINT1, PRINT2, PRINT3
   CHARACTER(LEN=8)  :: TPVAR,  TPRINT, TB2GR,  PRINTB2
   !REAL(SP) :: KTG11, KTG12, KTB1, KTG21, KTG22, KTB2, KTG31               !LB moved these 3 rows to the MOD_ALGAL definitions
   !REAL(SP) :: KTG32, KTB3,  KTPR
   !REAL(SP) :: TRPR, TMP1, TR1, TMP2, TR2, TMP3, TR3  !, TLOOK  !LB: TLOOK not needed anymore

!vjp  OPEN(UNIT=150,FILE='CCHL.OUT',STATUS='UNKNOWN')

! TITLE CARDS

     READ(AGR,1010) (ALG_TITLE(J),J=1,6)

! READ PARAMETERS WHICH ARE IDENTICAL FOR ALL GROUPS

     READ (AGR,1040)  MINPROD
     READ (AGR,1030)  TRPR,   KTPR
     READ (AGR,1030)  FNIP,   FNUP, FNDP,   FNLPP,   FNRPP  !need to add up to 1, FNUP not used at all, FNLDP,  FNRDP,
	 
	 
	 sanity = FNIP +  FNUP + FNDP + FNLPP + FNRPP
	 sant='Algal Predation Nitrogen'
	 CALL sanity_check(sant,sanity,1.0,0.011)
	 
	 
     READ (AGR,1030)  FPIP,   FPDP,  FPLPP,  FPRPP !FPLDP,  FPRDP
	  sanity = FPIP +  FPDP + FPLPP + FPRPP
	  sant='Algal phoshorous'
	 CALL sanity_check(sant,sanity,1.0,0.011)
	 
     READ (AGR,1030)  FDOP,   FCDP, FCLPP,  FCRPP  !FCLDP,  FCRDP, moved to wc_dom into one variable, FCDP 
	 sanity =   FDOP +  FCDP + FCLPP + FCRPP
	 sant='Algal Carbon'
	 CALL sanity_check(sant,sanity,1.0,0.011)
	 
     READ (AGR,1030)  FSAP
! READ SPATIALLY-INVARIANT PARAMETERS FOR ALGAL GROUP 1

     READ(AGR,1030) ANC1, Q01, ASC1, STF1                    !
     READ(AGR,1030) CHLCMN1, ACHLC1, BCHLC1, CCHLC1          !CHL to C ratios etc
     READ(AGR,1030) KHN1,KHP1,KHS1,KHST1,KHNFIX,FNFIX   !KH half constants   !,KHR1 B Clark moved to wc_dom
     READ(AGR,1030) ALPHMIN1,ALPHRAT1,PRSP1,VMAX1
     READ(AGR,1030) TMP1,TR1                !TR is temperature reference, TMP1 is optimal temperature for grwoth
     READ(AGR,1030) KTG11,KTG12,KTB1
     READ(AGR,1030) FNI1,FNLP1,FNRP1 !WLong: should check if they add up to 1 after reading  
												 !for FNLD1, FNRD1, FNLP1 , FNRP1  !!!!,FNLD1,FNRD1 B Clark moved to wc_dom
	 sanity = FNI1 +  FNLP1 + FNRP1 + FND1
	 sant='Algal1 Resp Nitrogen'
     READ(AGR,1030) FPI1,FPLP1,FPRP1  !,FPLD1,FPRD1
   !  READ(AGR,1030) FCLD1,FCRD1,FCLP1,FCRP1  ! B Clark change this to read in new paramaters for Fractions of DOM production
	 
	 sanity = FPI1 +  FPLP1 + FPRP1 + FPD1
	 sant='Algal1 Resp Phosphorus'
	! CALL sanity_check(sant,sanity,1.0,0.011)
	 
	 
     READ(AGR,1030)FCLP1,FCRP1   ! algal carbon from respiration
! READ SPATIALLY-INVARIANT PARAMETERS FOR ALGAL GROUP 2
	 sanity = FCLP1 +  FCRP1 + FCD1 ! algae don't lose carbon from respiration currently, this can have a range of values, as the balance all goes to "DIC"
	 sant='Algal1 Resp Carbon'
	! CALL sanity_check(sant,sanity,1.0,1.)
	 

     READ(AGR,1030) ANC2, Q02, ASC2, STF2
     READ(AGR,1030) CHLCMN2, ACHLC2, BCHLC2, CCHLC2      
     READ(AGR,1030) KHN2,KHP2,KHS2,KHST2   !KHR2 B clark moved to wc_dom
     READ(AGR,1030) ALPHMIN2,ALPHRAT2,PRSP2,VMAX2
     READ(AGR,1030) TMP2,TR2                !TR2 : temperaure reference   !TMP2: optimal temperature for growth
     READ(AGR,1030) KTG21,KTG22,KTB2
	 
     READ(AGR,1030) FNI2, FNLP2,FNRP2  !,FNLD2,FNRD2,
	 												 !for FNLD1, FNRD1, FNLP1 , FNRP1  !!!!,FNLD1,FNRD1 B Clark moved to wc_dom
	 sanity = FNI2 +  FNLP2 + FNRP2 + FND2! same as above for algae 2
	 sant='Algal2 Resp Nitrogen'
!	 CALL sanity_check(sant,sanity,1.0,0.011)
	 
     READ(AGR,1030) FPI2,FPLP2,FPRP2   !,FPLD2,FPRD2
	 sanity = FPI2 +  FPLP2 + FPRP2 + FPD2
	 sant='Algal1 Resp Phosphorus'
!	 CALL sanity_check(sant,sanity,1.0,0.011)
	 
     READ(AGR,1030) FCLP2,FCRP2  !FCLD2,FCRD2,
	 sanity = FCLP2 +  FCRP2 + FCD2 
	 sant='Algal1 Resp Carbon'
!	 CALL sanity_check(sant,sanity,1.0,1.)

! READ SPATIALLY-INVARIANT PARAMETERS FOR ALGAL GROUP 3

     READ(AGR,1030) ANC3, Q03, ASC3, STF3
     READ(AGR,1030) CHLCMN3, ACHLC3, BCHLC3, CCHLC3      
     READ(AGR,1030) KHN3,KHP3,KHS3,KHST3   ! :, B Clark moved to wc_dom
     READ(AGR,1030) ALPHMIN3,ALPHRAT3,PRSP3,VMAX3
     READ(AGR,1030) TMP3,TR3
     READ(AGR,1030) KTG31,KTG32,KTB3
     READ(AGR,1030) FNI3,FNLP3,FNRP3  !FNLD3,FNRD3,
     READ(AGR,1030) FPI3,FPLP3,FPRP3  !,FPLD3,FPRD3
     READ(AGR,1030) FCLP3,FCRP3  !FCLD3,FCRD3

!! CREATE LOOKUP TABLE FOR TEMPERATURE EFFECTS  !!LB: COMMENTED BECAUSE NOW WE USE FUNCTIONS GET_FT1,GET_FTBM1,ETC
!     DO I = -50,400
!       TLOOK = REAL(I)/10.

!       IF (TLOOK < TMP1) THEN                  !shape is different for < TMP1 and > TMP1
!         FT1(I) = EXP(-KTG11*(TLOOK-TMP1)**2)  !TMP1 is optimal temperature 
!       ELSE                                    !KTG11 is the shape parameter (bell shape) for T < TMP1
!         FT1(I) = EXP(-KTG12*(TMP1-TLOOK)**2)  !KTG12 is for T > TMP1
!       ENDIF

!       FTBM1(I)=EXP(KTB1*(TLOOK-TR1))          !Metabolism's temperature dependence


!       IF (TLOOK < TMP2) THEN                  !Look up table for ALG 2 growth rate
!         FT2(I) = EXP(-KTG21*(TLOOK-TMP2)**2)  !KTG21 is for shape with T < TMP2 
!       ELSE
!         FT2(I) = EXP(-KTG22*(TMP2-TLOOK)**2)  !KTG22 is for shape with T > TMP2
!       ENDIF

!       FTBM2(I)=EXP(KTB2*(TLOOK-TR2))          !
!       IF (TLOOK < TMP3) THEN                  !ALG 3
!         FT3(I) = EXP(-KTG31*(TLOOK-TMP3)**2)
!       ELSE
!         FT3(I) = EXP(-KTG32*(TMP3-TLOOK)**2)
!       ENDIF
!       FTBM3(I)=EXP(KTB3*(TLOOK-TR3))
!       FTPR(I)=EXP(KTPR*(TLOOK-TRPR))
!     ENDDO

     KHNAVG = (KHN1+KHN2+KHN3)/3.
     KHPAVG = (KHP1+KHP2+KHP3)/3.

! ARE REMAINING GROUP 1 PARAMETERS SPATIALLY VARYING?

     READ(AGR,1040) SPVAR1, PRINT1
     IF (SPVAR1 == 'CONSTANT') THEN
       READ(AGR,1030) PM1(1,1), BMR1(1,1), BPR1(1,1)
       DO K=1,KBM1
         DO I=1,MTLOC
           PM1(I,K) = PM1(1,1)
           BMR1(I,K) = BMR1(1,1)
           BPR1(I,K) = BPR1(1,1)
     ENDDO  
       ENDDO
     ELSE
       ALLOCATE(RTMP21(MGL,KBM1));     RTMP21 = 0.0
       ALLOCATE(RTMP22(MGL,KBM1));     RTMP22 = 0.0
       ALLOCATE(RTMP23(MGL,KBM1));     RTMP23 = 0.0
       DO K=1,KBM1
         DO I=1,MGL
           READ(AGR,1032) RTMP21(I,K), RTMP22(I,K), RTMP23(I,K)
     ENDDO  
       ENDDO
       IF(SERIAL)THEN
         PM1  = RTMP21
         BMR1 = RTMP22
         BPR1 = RTMP23
       ENDIF
#      if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         DO K=1,KBM1
           DO I=1,MLOC
             PM1(I,K)  = RTMP21(NGID(I),K)
             BMR1(I,K) = RTMP22(NGID(I),K)
             BPR1(I,K) = RTMP23(NGID(I),K)
           ENDDO
           DO I=1,NHN
             PM1(I+MLOC,K)  = RTMP21(HN_LST(I),K)
             BMR1(I+MLOC,K) = RTMP22(HN_LST(I),K)
             BPR1(I+MLOC,K) = RTMP23(HN_LST(I),K)
           ENDDO
         ENDDO
       ENDIF
#      endif
       DEALLOCATE(RTMP21,RTMP22,RTMP23)
     ENDIF      

! ARE REMAINING GROUP 2 PARAMETERS SPATIALLY VARYING?

     READ(AGR,1040) SPVAR2, PRINT2
     IF (SPVAR2 == 'CONSTANT') THEN
       READ(AGR,1030) PM2(1,1), BMR2(1,1), BPR2(1,1)
       DO K=1,KBM1
         DO I=1,MTLOC
           PM2(I,K) = PM2(1,1)
           BMR2(I,K) = BMR2(1,1)
           BPR2(I,K) = BPR2(1,1)
     ENDDO  
       ENDDO
     ELSE
       ALLOCATE(RTMP21(MGL,KBM1));    RTMP21 = 0.0
       ALLOCATE(RTMP22(MGL,KBM1));    RTMP22 = 0.0
       ALLOCATE(RTMP23(MGL,KBM1));    RTMP23 = 0.0
       DO K=1,KBM1
         DO I=1,MGL
           READ(AGR,1032) RTMP21(I,K), RTMP22(I,K), RTMP23(I,K)
     ENDDO  
       ENDDO
       IF(SERIAL)THEN
         PM2  = RTMP21
         BMR2 = RTMP22
         BPR2 = RTMP23
       ENDIF
#      if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         DO K=1,KBM1
           DO I=1,MLOC
             PM2(I,K)  = RTMP21(NGID(I),K)
             BMR2(I,K) = RTMP22(NGID(I),K)
             BPR2(I,K) = RTMP23(NGID(I),K)
           ENDDO
           DO I=1,NHN
             PM2(I+MLOC,K)  = RTMP21(HN_LST(I),K)
             BMR2(I+MLOC,K) = RTMP22(HN_LST(I),K)
             BPR2(I+MLOC,K) = RTMP23(HN_LST(I),K)
           ENDDO
         ENDDO
       ENDIF
#      endif
       DEALLOCATE(RTMP21,RTMP22,RTMP23)
     ENDIF      

! ARE REMAINING GROUP 3 PARAMETERS SPATIALLY VARYING?

     READ(AGR,1040) SPVAR3, PRINT3
     IF (SPVAR3 == 'CONSTANT') THEN
       READ(AGR,1030) PM3(1,1), BMR3(1,1), BPR3(1,1)
       DO K=1,KBM1
         DO I=1,MTLOC
           PM3(I,K) = PM3(1,1)
           BMR3(I,K) = BMR3(1,1)
           BPR3(I,K) = BPR3(1,1)
     ENDDO  
       ENDDO
     ELSE
       ALLOCATE(RTMP21(MGL,KBM1));    RTMP21 = 0.0
       ALLOCATE(RTMP22(MGL,KBM1));    RTMP22 = 0.0
       ALLOCATE(RTMP23(MGL,KBM1));    RTMP23 = 0.0
       DO K=1,KBM1
         DO I=1,MGL
           READ(AGR,1032) RTMP21(I,K), RTMP22(I,K), RTMP23(I,K)
     ENDDO  
       ENDDO
       IF(SERIAL)THEN
         PM3  = RTMP21
         BMR3 = RTMP22
         BPR3 = RTMP23
       ENDIF
#      if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         DO K=1,KBM1
           DO I=1,MLOC
             PM3(I,K)  = RTMP21(NGID(I),K)
             BMR3(I,K) = RTMP22(NGID(I),K)
             BPR3(I,K) = RTMP23(NGID(I),K)
           ENDDO
           DO I=1,NHN
             PM3(I+MLOC,K)  = RTMP21(HN_LST(I),K)
             BMR3(I+MLOC,K) = RTMP22(HN_LST(I),K)
             BPR3(I+MLOC,K) = RTMP23(HN_LST(I),K)
           ENDDO
         ENDDO
       ENDIF
#      endif
       DEALLOCATE(RTMP21,RTMP22,RTMP23)
     ENDIF

! TIME DEPENDENCE OF PREDATION BY HIGHER TROPHIC LEVELS

     READ(AGR,1040) TPVAR, TPRINT   
     IF (TPVAR == 'CONSTANT') THEN
       READ(AGR,1060) TVPR(1)
       DO J=2,366
         TVPR(J) = TVPR(1)
       ENDDO
     ELSE
       READ(AGR,1060) (TVPR(J), J=1,366)
     ENDIF 

! TIME DEPENDENCE OF GROWTH BY SPRING ALGAL GROUP

     READ(AGR,1040) TB2GR, PRINTB2
     IF (TB2GR == 'CONSTANT') THEN
       READ(AGR,1060) B2GR(1)
       DO J=2,366
         B2GR(J) = B2GR(1)
       ENDDO
     ELSE
       READ(AGR,1060) (B2GR(J), J=1,366)   !daily series of B2GR 
     ENDIF 

!***** Input FORMAT statements

1010 FORMAT(A72)
1030 FORMAT(//(8X,9F8.0))
1032 FORMAT(8X,9F8.0)
1040 FORMAT(//8X,8A8)
1060 FORMAT(//(16X,F8.0))

! OUTPUT WHAT WAS INPUT
     OPEN(ALO,FILE=ALOFN)
     IF(MSR)WRITE(ALO,2002) (ALG_TITLE(J),J=1,6)

! WRITE ALGAL PROPERTIES WHICH ARE NOT SPATIALLY VARYING

     IF(MSR)WRITE(ALO,3000)
     IF(MSR)WRITE(ALO,3010) CHLCMN1, CHLCMN2, CHLCMN3
     IF(MSR)WRITE(ALO,3012) ACHLC1, ACHLC2, ACHLC3
     IF(MSR)WRITE(ALO,3014) BCHLC1, BCHLC2, BCHLC3
     IF(MSR)WRITE(ALO,3016) CCHLC1, CCHLC2, CCHLC3
     IF(MSR)WRITE(ALO,3020) ANC1, ANC2, ANC3
     IF(MSR)WRITE(ALO,3022) Q01,  Q02,  Q03
     IF(MSR)WRITE(ALO,3030) ASC1, ASC2, ASC3
     IF(MSR)WRITE(ALO,3072) STF1, STF2, STF3
     IF(MSR)WRITE(ALO,3040) KHN1, KHN2, KHN3
     IF(MSR)WRITE(ALO,3050) KHP1, KHP2, KHP3
     IF(MSR)WRITE(ALO,3042) KHNFIX
     IF(MSR)WRITE(ALO,3043) FNFIX
     IF(MSR)WRITE(ALO,3060) KHS1, KHS2, KHS3
   !  IF(MSR)WRITE(ALO,3070) KHR1, KHR2, KHR3
     IF(MSR)WRITE(ALO,3074) KHST1,KHST2,KHST3
     IF(MSR)WRITE(ALO,3076) VMAX1, VMAX2, VMAX3
     IF(MSR)WRITE(ALO,3080) ALPHMIN1, ALPHMIN2, ALPHMIN3
     IF(MSR)WRITE(ALO,3090) ALPHRAT1, ALPHRAT2, ALPHRAT3
     IF(MSR)WRITE(ALO,4000) PRSP1, PRSP2, PRSP3
     IF(MSR)WRITE(ALO,4010) TMP1, TMP2, TMP3
     IF(MSR)WRITE(ALO,5020) TR1, TR2, TR3
     IF(MSR)WRITE(ALO,5030) KTG11, KTG21, KTG31
     IF(MSR)WRITE(ALO,5040) KTG12, KTG22, KTG32
     IF(MSR)WRITE(ALO,5050) KTB1, KTB2, KTB3

! WRITE DISTRIBUTION OF BIOMASS UPON MORTALITY

     IF(MSR)WRITE(ALO,4020)
     IF(MSR)WRITE(ALO,4030) FNI1,FNLP1,FNRP1 !,FNLD1,FNRD1   ! FNI1,FNLP1 etc + FND1 must add up to 1!!!!!
     IF(MSR)WRITE(ALO,4040) FPI1,FPLP1,FPRP1 !,FPLD1,FPRD1
     IF(MSR)WRITE(ALO,4050) FCLP1,FCRP1  ! FCLD1,FCRD1,
     IF(MSR)WRITE(ALO,4060)
     IF(MSR)WRITE(ALO,4030) FNI2,FNLP2,FNRP2 !,FNLD2,FNRD2
     IF(MSR)WRITE(ALO,4040) FPI2,FPLP2,FPRP2 !,FPLD2,FPRD2
     IF(MSR)WRITE(ALO,4050) FCLP2,FCRP2 !FCLD2,FCRD2
     IF(MSR)WRITE(ALO,4070)
     IF(MSR)WRITE(ALO,4030) FNI3,FNLP3,FNRP3 !,FNLD3,FNRD3
     IF(MSR)WRITE(ALO,4040) FPI3,FPLP3,FPRP3 !,FPLD3,FPRD3
     IF(MSR)WRITE(ALO,4050) FCLP3,FCRP3 !FCLD3,FCRD3
     IF(MSR)WRITE(ALO,4090)
     IF(MSR)WRITE(ALO,4030) FNIP,FNDP,FNLPP,FNRPP
     IF(MSR)WRITE(ALO,4040) FPIP,FPDP,FPLPP,FPRPP
     IF(MSR)WRITE(ALO,4080) FDOP,FCDP,FCLPP,FCRPP
     IF(MSR)WRITE(ALO,4082) FSAP,(1.-FSAP)
 
     IF (MINPROD == ' MINIMUM' .OR. MINPROD == ' PRODUCT') THEN
       IF(MSR)WRITE(ALO,2010) MINPROD
     ELSE
       IF(MSR)WRITE(ALO,2020) MINPROD
       STOP
     ENDIF

! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 1

     IF(MSR)WRITE(ALO,2000)
     IF (SPVAR1 == 'CONSTANT') THEN
       IF(MSR)WRITE(ALO,2085)
     ELSE
       IF(MSR)WRITE(ALO,2086)
     ENDIF
     IF(MSR)WRITE(ALO,5060)
     IF (PRINT1 /= '     ALL') THEN
       IF(MSR)WRITE(ALO,5070) PM1(1,1), BMR1(1,1), BPR1(1,1)
     ELSE
       DO K=1,KBM1
         DO I=1,MLOC
       IF(MSR)WRITE(ALO,5080) I, K, PM1(I,K), BMR1(I,K), BPR1(I,K)
     ENDDO
       ENDDO      
     ENDIF
      
! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 2

     IF(MSR)WRITE(ALO,2005)
     IF (SPVAR2 == 'CONSTANT') THEN
       IF(MSR)WRITE(ALO,2085)
     ELSE
       IF(MSR)WRITE(ALO,2086)
     ENDIF
     IF(MSR)WRITE(ALO,5060)
     IF (PRINT2 /= '     ALL') THEN
       IF(MSR)WRITE(ALO,5070) PM2(1,1), BMR2(1,1), BPR2(1,1)
     ELSE
       DO K=1,KBM1
         DO I=1,MLOC
           IF(MSR)WRITE(ALO,5080) I, K, PM2(I,K), BMR2(I,K), BPR2(I,K)
     ENDDO
       ENDDO      
     ENDIF

! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 3

     IF(MSR)WRITE(ALO,2006)
     IF (SPVAR3 == 'CONSTANT') THEN
       IF(MSR)WRITE(ALO,2085)
     ELSE
       IF(MSR)WRITE(ALO,2086)
     ENDIF
     IF(MSR)WRITE(ALO,5060)
     IF (PRINT3 /= '     ALL') THEN
       IF(MSR)WRITE(ALO,5070) PM3(1,1), BMR3(1,1), BPR3(1,1)
     ELSE
       DO K=1,KBM1
         DO I=1,MLOC
           IF(MSR)WRITE(ALO,5080) I, K, PM3(I,K), BMR3(I,K), BPR3(I,K)
     ENDDO
       ENDDO      
     ENDIF

! WRITE TERMS FOR PREDATION BY HIGHER TROPHIC LEVELS

     IF (TPVAR == 'CONSTANT') THEN
       IF(MSR)WRITE(ALO,2087)
     ELSE
       IF(MSR)WRITE(ALO,2088)
     ENDIF 
     IF(MSR)WRITE(ALO,2040) TRPR, KTPR
     IF(MSR)WRITE(ALO,2092)
     IF (TPRINT /= '     ALL') THEN
       IF(MSR)WRITE(ALO,2090) (J, TVPR(J), J=1,1)
     ELSE
       IF(MSR)WRITE(ALO,2090) (J, TVPR(J), J=1,366)
     ENDIF
      
! WRITE TIME DEPENDENCE OF GROWTH OF SPRING ALGAL GROUP

     IF (TB2GR == 'CONSTANT') THEN
       IF(MSR)WRITE(ALO,2094)
     ELSE
       IF(MSR)WRITE(ALO,2095)
     ENDIF 
     IF(MSR)WRITE(ALO,2093)
     IF (PRINTB2 /= '     ALL') THEN
       IF(MSR)WRITE(ALO,2090) (J, B2GR(J), J=1,1)
     ELSE
       IF(MSR)WRITE(ALO,2090) (J, B2GR(J), J=1,366)
     ENDIF
      
!***** Output FORMAT statements

2000 FORMAT(/' ALGAL GROUP 1')
2002 FORMAT(1X,A72)
2005 FORMAT(/' ALGAL GROUP 2')
2006 FORMAT(/' ALGAL GROUP 3')
2010 FORMAT(/A8,' FORMULATION SELECTED FOR LIGHT AND NUTRIENT',' LIMITATIONS')
2020 FORMAT(/' MINPROD INCORRECTLY SPECIFIED AS ',A8)
2040 FORMAT(' PREDATION SPECIFIED AT ',F8.2,' C.'/              &
            ' TEMPERATURE EFFECT = ',F8.3,' PER DEGREE')
2085 FORMAT(/' REMAINING PARAMETERS ARE SPATIALLY-INVARIANT')
2086 FORMAT(/' REMAINING PARAMETERS ARE SPATIALLY-VARYING')
2087 FORMAT(/' HIGHER-LEVEL PREDATION IS TEMPORALLY-INVARIANT')
2088 FORMAT(/' HIGHER-LEVEL PREDATION VARIES TEMPORALLY')
2090 FORMAT(I8,F8.3)
2092 FORMAT(/'    DAY     TVPR')
2093 FORMAT(/'    DAY     B2GR')
2094 FORMAT(/' BASE GROUP 2 GROWTH IS TEMPORALLY-INVARIANT')
2095 FORMAT(/' BASE GROUP 2 VARIES TEMPORALLY')
3000 FORMAT(/' ALGAL PROPERTIES',30X,'GROUP 1   GROUP 2   GROUP 3')
3010 FORMAT(/' MINIMUM GM CHL/GM C    ',20X,3F10.4)
3012 FORMAT(' ACHCL (GM CHL/GM C)    ',20X,3F10.4)
3014 FORMAT(' BCHCL (1/DEGREE C)     ',20X,3F10.4)
3016 FORMAT(' CCHCL (M**2 DAY/E)     ',20X,3F10.4)
3020 FORMAT(' GM N/GM C              ',20X,3F10.4)
3022 FORMAT(' MINIMUM GM P/GM C      ',20X,3F10.4)
3030 FORMAT(' GM SI/GM C             ',20X,3F10.4)
3040 FORMAT(' KHN (GM N/M**3)        ',20X,3F10.4)
3042 FORMAT(' KHNFIX (GM N/M**3)     ',20X,F10.4)
3043 FORMAT(' FRACTION N FIXERS      ',20X,F10.4)
3050 FORMAT(' KHP (GM P/M**3)        ',20X,3F10.4)
3060 FORMAT(' KHS (GM SI/M**3)       ',20X,3F10.4)
3070 FORMAT(' KHR (GM DO/M**3)       ',20X,3F10.4)
3072 FORMAT(' SALT TOXICITY (1/DAY)  ',20X,3F10.4)
3074 FORMAT(' KHST (ppt)             ',20X,3F10.4)
3076 FORMAT(' VMAX (GM P/GM C/DAY)   ',20X,3F10.4)
3080 FORMAT(' ALPHA (GM C/GM CHL/DAY)/(uE/M**2/SEC)      ',3F10.3)
3090 FORMAT(' ALPHA TO P RATIO (1/(uE/M**2/SEC))',9X,3F10.3)
4000 FORMAT(' PHOTORESPIRATION FRACTION     ',13X,3F10.3)
4010 FORMAT(' OPTIMAL TEMPERATURE FOR PRODUCTION (C)',5X,3F10.3) 
4020 FORMAT(/' DISTRIBUTION OF ALGAE UPON MORTALITY'//                       &
     ' GROUP 1 RESPIRATION  DIS INORG  LAB DISS  REF DISS  LAB PART  REF PART')
4030 FORMAT(' NITROGEN            ',5F10.3)
4040 FORMAT(' PHOSPHORUS          ',5F10.3)
4050 FORMAT(' CARBON              ',10X,4F10.3)
4060 FORMAT(/' GROUP 2 RESPIRATION')
4070 FORMAT(/' GROUP 3 RESPIRATION')
4080 FORMAT(' CARBON              ',5F10.3)
4082 FORMAT(' SILICA              ',F10.3,10X,F10.3)
4090 FORMAT(/' PREDATION')
5000 FORMAT(' PHOSPHORUS     ',4F10.3)
5010 FORMAT(' CARBON         ',4F10.3)
5020 FORMAT(' REFERENCE TEMPERATURE FOR RESPIRATION (C)  ',3F10.3)
5030 FORMAT(' EFFECT OF SUBOPTIMAL TEMP ON PRODUCTION    ',3F10.4) 
5040 FORMAT(' EFFECT OF SUPEROPTIMAL TEMP ON PRODUCTION  ',3F10.4) 
5050 FORMAT(' EFFECT OF TEMPERATURE ON RESPIRATION       ',3F10.4)
5060 FORMAT('       CELL     PMAX     METAB     PRDTN',/                   &
            '             C/CHL/DAY   1/DAY     1/DAY'/)
5070 FORMAT(10X,3F10.3)
5080 FORMAT(2I10,3F10.3)

     RETURN
     END SUBROUTINE ALG_READ
      


!************************************************************************
!**                  S U B R O U T I N E   A L G A E                   **
!************************************************************************

   SUBROUTINE ALGAE(DTB1,DTB2,DTB3,FLXS1,FLXS2,FLXS3)
   IMPLICIT NONE
   INTEGER  ::  I, J, F, K, ITEMP
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: DTB1, DTB2,DTB3
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: FLXS1,FLXS2,FLXS3
   REAL(SP),ALLOCATABLE :: XL(:), STOX1(:,:), STOX2(:,:), STOX3(:,:)
   REAL(SP),ALLOCATABLE :: ITBOT(:,:), NL1NNF(:,:)
   REAL(SP) :: FI01, FI02, FI03, GPP1, GPP2, GPP3, TREC
   REAL(SP) :: NETP1, NETP2, NETP3, IK, ITTOP, ITAVG, NUTLIM
   REAL(SP) :: SALTOX, DIN, DSIL, ALPHA, OPTDEPTH

   
   ALLOCATE(XL(MLOC));          XL = 0.0
   ALLOCATE(STOX1(MLOC,KBM1));  STOX1 = 0.0
   ALLOCATE(STOX2(MLOC,KBM1));  STOX2 = 0.0
   ALLOCATE(STOX3(MLOC,KBM1));  STOX3 = 0.0
   ALLOCATE(NL1NNF(MLOC,KBM1)); NL1NNF = 0.0
   ALLOCATE(ITBOT(MLOC,KBM1));  ITBOT = 0.0

! DETERMINE JULIAN DAY

   J = 1.0 + AMOD(JDAY,365.25)
!     if (jday .gt. 730.) fnfix = 0.9

!!! SALINITY TOXICITY
!!!RGL commented this out for runs in the sound - take this out for Skagit
!!   DO K=1,KBM1
!!     DO I=1,MLOC
!!       SALTOX = MAX(0., SALT(I,K))
!!       STOX1(I,K)=STF1*0.5*(1.+TANH(SALTOX-KHST1))
!!       STOX2(I,K)=STF2*(1.-0.5*(1.+TANH(SALTOX-KHST2)))
!!       STOX3(I,K)=STF3*0.5*(1.+TANH(SALTOX-KHST3))
!!     ENDDO
!!   ENDDO  
!! SALINITY TOXICITY- use code below instead
!
!   DO K=1,KBM1
!     DO I=1,MLOC
!       SALTOX = MAX(0., SALT(I,K))
!! 14 SEPT 2009 KURT GLAESEMANN, ADD IF, SINCE TANH CAN BE COSTLY, BUT ZERO IS EASY
!       if(STF1 .ne. 0) STOX1(I,K)=STF1*0.5*(1.+TANH(SALTOX-KHST1))
!       if(STF2 .ne. 0) STOX2(I,K)=STF2*(1.-0.5*(1.+TANH(SALTOX-KHST2)))
!       if(STF3 .ne. 0) STOX3(I,K)=STF3*0.5*(1.+TANH(SALTOX-KHST3))
!     ENDDO
!   ENDDO  
! TEMPERATURE EFFECTS ON PRODUCTION AND RESPIRATION

   DO K=1,KBM1
     DO I=1,MLOC
       !ITEMP = 10.* (max(T(I,K),0.0)) +0.05         !LB: added max( T, 0.) to ensure positive temperature
!       P1(I,K) = PM1(I,K)*GET_FT1(T(I,K))            !!LB: replace FT1(ITEMP)-lookup-table with "get_FT1(T)"
!       P2(I,K) = PM2(I,K)*GET_FT2(T(I,K))*B2GR(J)
!!       P3(I,K) = PM3(I,K)*GET_FT3(T(I,K))

!       BM1(I,K) = BMR1(I,K)*GET_FTBM1(T(I,K))+STOX1(I,K)   !ALG1 metabolism 
!       BM2(I,K) = BMR2(I,K)*GET_FTBM2(T(I,K))+STOX2(I,K)   !ALG2 metabolism
!!       BM3(I,K) = BMR3(I,K)*GET_FTBM3(T(I,K))+STOX3(I,K)
!       
!	   PR1(I,K) = BPR1(I,K)*B1(I,K)*B1(I,K)*GET_FTPR(T(I,K))*TVPR(J)  !gC/m^3/day ??  !Laura: should not repeat B1	   
!       PR2(I,K) = BPR2(I,K)*B2(I,K)*B2(I,K)*GET_FTPR(T(I,K))*TVPR(J)	 !Laura: should not repeat B2
!!       PR3(I,K) = BPR3(I,K)*B3(I,K)*B3(I,K)*GET_FTPR(T(I,K))*TVPR(J) !Laura: should not repeat B3
       
       
       P1(I,K) = PM1(I,K)*GET_FT(T(I,K),TMP1,KTG11,KTG21)            !!LB: replace FT1(ITEMP)-lookup-table with "get_FT1(T)" -NOW get_FT(T,TMP,KTG1,KTG2)
       P2(I,K) = PM2(I,K)*GET_FT(T(I,K),TMP2,KTG21,KTG22)*B2GR(J)   ! B2Gr = 1
!       P3(I,K) = PM3(I,K)*GET_FT(T(I,K),TMP3,KTG31,KTG32)

       BM1(I,K) = BMR1(I,K)*GET_FTBM(T(I,K),TR1,KTB1)+STOX1(I,K)   !ALG1 metabolism 
       BM2(I,K) = BMR2(I,K)*GET_FTBM(T(I,K),TR2,KTB2)+STOX2(I,K)   !ALG2 metabolism
!       BM3(I,K) = BMR3(I,K)*GET_FTBM(T(I,K),TR3,KTB3)+STOX3(I,K)
       
	   PR1(I,K) = BPR1(I,K)*B1(I,K)*B1(I,K)*GET_FTPR(T(I,K),TRPR,KTPR)*TVPR(J)  !gC/m^3/day ??  !Laura: should not repeat B1	   
        PR2(I,K) = BPR2(I,K)*B2(I,K)*B2(I,K)*GET_FTPR(T(I,K),TRPR,KTPR)*TVPR(J)	 !Laura: should not repeat B2
!       PR3(I,K) = BPR3(I,K)*B3(I,K)*B3(I,K)*GET_FTPR(T(I,K),TRPR,KTPR)*TVPR(J) !Laura: should not repeat B3
     ENDDO 
   ENDDO  
   
! NUTRIENT LIMITATIONS - NITROGEN

   DO K=1,KBM1
     DO I=1,MLOC
       DIN = NH4(I,K)+NO3(I,K)
!RGL       NL1(I,K) = (1.-FNFIX)*DIN/(KHN1+DIN) + FNFIX*(DIN+KHNFIX)/(KHN1+DIN+KHNFIX)

!WLONG       NL1(I,K) = DIN/(KHN1+DIN)
!WLONG       NL1NNF(I,K) = DIN/(KHN1+DIN) 
!WLONG       NL2(I,K) = DIN/(KHN2+DIN)
!JQI       NL3(I,K) = DIN/(KHN3+DIN)

!Wen Long : fix bug in nitrogen limiation (one should consider ammonia preference PN)
!           when PN -->1,       (i.e. KHN becomes small), NL should be NH4/(KHN+NH4) 
!           when PN --> NH4/DIN,(i.e. KHN becomes large), NL should be DIN/(KHN+DIN)

           !ammonia preference  (moved from wqm_kin.F to here)  
           PN1(I,K) = NH4(I,K)*(NO3(I,K)/((KHN1+NH4(I,K))*(KHN1+NO3(I,K)))+KHN1  &
                /((1.E-30+NH4(I,K)+NO3(I,K))*(KHN1+NO3(I,K))))
           PN2(I,K) = NH4(I,K)*(NO3(I,K)/((KHN2+NH4(I,K))*(KHN2+NO3(I,K)))+KHN2  &
                /((1.E-30+NH4(I,K)+NO3(I,K))*(KHN2+NO3(I,K))))
!WLONG       PN3(I,K) = NH4(I,K)*(NO3(I,K)/((KHN3+NH4(I,K))*(KHN3+NO3(I,K)))+KHN3  &
!WLONG                /((1.E-30+NH4(I,K)+NO3(I,K))*(KHN3+NO3(I,K))))

           NL1(I,K) = ( NH4(I,K) + DIN - PN1(I,K)*DIN)   &
                     /(KHN1+NH4(I,K) + DIN - PN1(I,K)*DIN+1E-30)
           NL1NNF(I,K) = (NH4(I,K) + DIN - PN1(I,K)*DIN)   &
                        /(KHN1+NH4(I,K) + DIN - PN1(I,K)*DIN+1E-30)
           NL2(I,K) = (NH4(I,K) + DIN - PN2(I,K)*DIN)   &
                     /(KHN2+NH4(I,K) + DIN - PN2(I,K)*DIN+1E-30)

!WLONG     NL3(I,K) = (     NH4(I,K) + DIN - PN(I,K)*DIN)   &
!WLONG               /(KHN3+NH4(I,K) + DIN - PN(I,K)*DIN+1E-30)

!     ENDDO 
!   ENDDO  


!! PHOSPHORUS
!!! DROOP FORMULA FOR P LIMIT
!!! RGL need to check
!   DO K=1,KBM1
!     DO I=1,MLOC
!below is Droop - RGL converted to Michaelis-Menten b/c Droop issues
!       PL1(I,K) = MAX((Q1(I,K)-Q01)/(Q1(I,K)+1.0E-30),0.0)
!       PL2(I,K) = MAX((Q2(I,K)-Q02)/(Q2(I,K)+1.0E-30),0.0)
!!      PL3(I,K) = MAX((Q3(I,K)-Q03)/(Q3(I,K)+1.0E-30),0.0)
!Michaelis-Menten here...
        PL1(I,K)= PO4(I,K)/(KHP1+PO4(I,K))
        PL2(I,K)= PO4(I,K)/(KHP2+PO4(I,K))
!        PL3(I,K)= PO4(I,K)/(KHP3+PO4(I,K))
     ENDDO 
   ENDDO  

!! SILICA 

!   DO K=1,KBM1
!     DO I=1,MLOC
!       DSIL = SIAT(I,K)/(1.+KADSA*SSI(I,K))             !Wen Long: SIAT is total concentration of silca
!       SL1(I,K) = (DSIL+1.0E-30)/(KHS1+DSIL+1.0E-30)
!       SL2(I,K) = (DSIL+1.0E-30)/(KHS2+DSIL+1.0E-30)
!       SL3(I,K) = (DSIL+1.0E-30)/(KHS3+DSIL+1.0E-30)
!     ENDDO
!   ENDDO  

! PHOTOSYNTHESIS VS IRRADIANCE

!vjp 10/07/04  modified sqrt formulae in next loop
   DO K=1,KBM1
     DO I=1,MLOC
       ALPHA = ALPHMIN1+ALPHRAT1*P1(I,K)
! RGL why are they using lowered Pmax instead of true Pmax?
! try this with true pmax and redo with Ek instead of alpha/Pmax
!       IK = P1(I,K)*AMIN1(NL1(I,K),PL1(I,K),SL1(I,K))/(ALPHA+1.0E-10)
!below without SiLim
!       IK = P1(I,K)*AMIN1(NL1(I,K),PL1(I,K))/(ALPHA+1.0E-10)

!       IK = P1(I,K)*NL1(I,K)/(ALPHA+1.0E-10)
!RGL changed IK to reflect Pmax/alpha - not modified P
       IK= P1(I,K)/(ALPHA+1.0E-10)        !Wen Long, I believe this is correct
                                          !I agree with RGL that IK calculation should not have 
                                          !nitrogen limitation in it
       IK1(I,K)=IK
!IF ( IAVG(I,K) > 1.) THEN
   ! write(*,*)'IAVG in algal mod= ',IAVG(I,K)
       FI1(I,K) = IAVG(I,K)/(SQRT(IK*IK+IAVG(I,K)*IAVG(I,K)+1.0E-10))
!ENDIF
!      FI1(B) = (IAVG(B)+1.0E-10)/(SQRT(IK*IK+IAVG(B)*IAVG(B)+1.0E-10))

       ALPHA = ALPHMIN2+ALPHRAT2*P2(I,K)
!       IK = P2(I,K)*AMIN1(NL2(I,K),PL2(I,K),SL2(I,K))/(ALPHA+1.0E-10)
!RGL below without P, Si limitation
!       IK = P2(I,K)*NL2(I,K)/(ALPHA+1.0E-10)
       IK= P2(I,K)/(ALPHA+1.0E-10)      !Wen Long, I believe this is correct 
	   
       IK2(I,K)=IK

       FI2(I,K) = IAVG(I,K)/(SQRT(IK*IK+IAVG(I,K)*IAVG(I,K)+1.0E-10))
!      FI2(B) = (IAVG(B)+1.0E-10)/(SQRT(IK*IK+IAVG(B)*IAVG(B)+1.0E-10))

!       ALPHA = ALPHMIN3+ALPHRAT3*P3(I,K)
!!       IK = P3(I,K)*AMIN1(NL3(I,K),PL3(I,K),SL3(I,K))/(ALPHA+1.0E-10)
!!RGL below without P, Si limitation
!       IK = P3(I,K)*NL3(I,K)/(ALPHA+1.0E-10)
!       FI3(I,K) = IAVG(I,K)/(SQRT(IK*IK+IAVG(I,K)*IAVG(I,K)+1.0E-10))
     ENDDO
    ENDDO  

!     *** Compute carbon to chlorophyll ratio in each cell 
! 14 SEPT 2009 KURT GLAESEMANN - THESE ARE SOMETIME ALL ZERO - MUCH SIMPLER CODE
   IF (ACHLC1 .EQ. 0 .AND. ACHLC2 .EQ. 0 .AND. ACHLC3 .EQ. 0) THEN
   DO K=1,KBM1
     DO I=1,MLOC
       CCHL1(I,K) = 1./CHLCMN1
       CCHL2(I,K) = 1./CHLCMN2
       CCHL3(I,K) = 1./CHLCMN3
     ENDDO
   ENDDO

   ELSE ! 14 SEPT 2009
   
   DO I=1,MLOC
   

#if defined(PHOTODEGRADATION) 
     ITAVG  = IAVG(I,1)
	! write(*,*)'ITAVG in mod_algal = ',ITAVG
#else
	 ITTOP = IT
     OPTDEPTH = KESS(I,1)*D(I)*DZ(1)
     ITBOT(I,1) = ITTOP*EXP(-OPTDEPTH)
     ITAVG = (ITTOP-ITBOT(I,1))/OPTDEPTH
#endif	 
	 
     NUTLIM = AMIN1(NL1(I,1),PL1(I,1))
! below without Plim
!     NUTLIM = NL1(I,1)
     CCHL1(I,1) = 1./(CHLCMN1+ACHLC1*EXP(BCHLC1*T(I,1))*NUTLIM/EXP(CCHLC1*ITAVG))
     NUTLIM = AMIN1(NL2(I,1),PL2(I,1))
!     NUTLIM = NL2(I,1)
     CCHL2(I,1) = 1./(CHLCMN2+ACHLC2*EXP(BCHLC2*T(I,1))*NUTLIM/EXP(CCHLC2*ITAVG))
!     NUTLIM = AMIN1(NL3(I,1),PL3(I,1))
!     NUTLIM = NL3(I,1)
!     CCHL3(I,1) = 1./(CHLCMN3+ACHLC3*EXP(BCHLC3*T(I,1))*NUTLIM/EXP(CCHLC3*ITAVG))
   ENDDO 

   DO K=2,KBM1
     DO I=1,MLOC
#if defined(PHOTODEGRADATION) 
       ITAVG  = IAVG(I,K)
#else
	   ITTOP = IT
       OPTDEPTH = KESS(I,K)*D(I)*DZ(K)
       ITBOT(I,K) = ITTOP*EXP(-OPTDEPTH)
       ITAVG = (ITTOP-ITBOT(I,1))/OPTDEPTH
#endif	 
       NUTLIM = AMIN1(NL1(I,K),PL1(I,K))
!RGL below without P, Si limitation
!       NUTLIM = NL1(I,K)
       CCHL1(I,K) = 1./(CHLCMN1+ACHLC1*EXP(BCHLC1*T(I,K))*NUTLIM/EXP(CCHLC1*ITAVG))
       NUTLIM = AMIN1(NL2(I,K),PL2(I,K))
!       NUTLIM = NL2(I,K)
       CCHL2(I,K) = 1./(CHLCMN2+ACHLC2*EXP(BCHLC2*T(I,K))*NUTLIM/EXP(CCHLC2*ITAVG))
!       NUTLIM = AMIN1(NL3(I,K),PL3(I,K))
!       CCHL3(I,K) = 1./(CHLCMN3+ACHLC3*EXP(BCHLC3*T(I,K))*NUTLIM/EXP(CCHLC3*ITAVG))
     ENDDO
   ENDDO 
   ENDIF ! ENDIF 14 SEPT 2009

!END RGL output section
! COMPUTE ASSIMILATION RATIO AT WATER SURFACE (GM C/GM CHL/DAY), 
! CARBON FIXATION (GM C/M**3/DAY)  COPIED FROM CHES BAY CODE MARCH 14, 2007

   
!   DO I=1,MLOC
!!
!!     ASRAT(I) = (P1(I,1)*AMIN1(NL1(I,1),PL1(I,1),SL1(I,1))*B1(I,1)                  &
!!               + P2(I,1)*AMIN1(NL2(I,1),PL2(I,1),SL2(I,1))*B2(I,1)        &
!!               + P3(I,1)*AMIN1(NL3(I,1),PL3(I,1),SL3(I,1))*B3(I,1))       &
!!               / (B1(I,1)+B2(I,1)+B3(I,1)+1.0e-6)
!!RGL below without P or Si limitation
!!     ASRAT(I) = (P1(I,1)*NL1(I,1)*B1(I,1)                  &
!!               + P2(I,1)*NL2(I,1)*B2(I,1)        &
!!               + P3(I,1)*NL3(I,1)*B3(I,1))       &
!!               / (B1(I,1)+B2(I,1)+B3(I,1)+1.0e-6)
!      ASRAT(I) = (P1(I,1)*AMIN1(NL1(I,1),PL1(I,1))*B1(I,1)       &
!                + P2(I,1)*AMIN1(NL2(I,1),PL2(I,1))*B2(I,1)) /     &
!                 (B1(I,1)+B2(I,1)+1.0e-6)

!!     NETP1 = (P1(I,1)*AMIN1(NL1(I,1),PL1(I,1),SL1(I,1))*(1.-PRSP1)        &
!!             /CCHL1(I,1)-BM1(I,1))*B1(I,1)
!!     NETP2 = (P2(I,1)*AMIN1(NL2(I,1),PL2(I,1),SL2(I,1))*(1.-PRSP2)        &
!!             /CCHL2(I,1)-BM2(I,1))*B2(I,1)
!!     NETP3 = (P3(I,1)*AMIN1(NL3(I,1),PL3(I,1),SL3(I,1))*(1.-PRSP3)        &
!!             /CCHL3(I,1)-BM3(I,1))*B3(I,1)
!! RGL below without P or Si limitation
!     NETP1 = (P1(I,1)*AMIN1(NL1(I,1),PL1(I,1))*(1.-PRSP1)        &
!             /CCHL1(I,1)-BM1(I,1))*B1(I,1)
!     NETP2 = (P2(I,1)*AMIN1(NL2(I,1),PL2(I,1))*(1.-PRSP2)        &
!             /CCHL2(I,1)-BM2(I,1))*B2(I,1)
!!     NETP3 = (P3(I,1)*NL3(I,1)*(1.-PRSP3)        &
!!             /CCHL3(I,1)-BM3(I,1))*B3(I,1)

!     CFIX(I) = NETP1+NETP2!+NETP3
!   ENDDO

! EFFECTS OF LIGHT AND NUTRIENTS ON PRODUCTION
! CONVERT FROM GM C/GM CHL/DAY TO SPECIFIC GROWTH RATE

   IF (MINPROD == ' MINIMUM') THEN

     DO K=1,KBM1
       DO I=1,MLOC
!         P1NNF(I,K) = P1(I,K)*AMIN1(FI1(I,K),NL1NNF(I,K),PL1(I,K),SL1(I,K))  &
!                /CCHL1(I,K)            
!         P1(I,K) = P1(I,K)*AMIN1(FI1(I,K),NL1(I,K),PL1(I,K),SL1(I,K))/CCHL1(I,K)
!         P2(I,K) = P2(I,K)*AMIN1(FI2(I,K),NL2(I,K),PL2(I,K),SL2(I,K))/CCHL2(I,K)
!         P3(I,K) = P3(I,K)*AMIN1(FI3(I,K),NL3(I,K),PL3(I,K),SL3(I,K))/CCHL3(I,K)
!RGL below without P or Si limitation
         P1NNF(I,K) = P1(I,K)*AMIN1(FI1(I,K),NL1NNF(I,K),PL1(I,K))/CCHL1(I,K)
         P1(I,K)    = P1(I,K)*AMIN1(FI1(I,K),NL1(I,K),PL1(I,K))/CCHL1(I,K)
         P2(I,K)    = P2(I,K)*AMIN1(FI2(I,K),NL2(I,K),PL2(I,K))/CCHL2(I,K)
!         P3(I,K) = P3(I,K)*AMIN1(FI3(I,K),NL3(I,K))/CCHL3(I,K)
       ENDDO    
     ENDDO

   ELSE

     DO K=1,KBM1
       DO I=1,MLOC
!         P1NNF(I,K) = P1(I,K)*FI1(I,K)*AMIN1(NL1NNF(I,K),PL1(I,K),SL1(I,K))  &
!                /CCHL1(I,K)
!         P1(I,K) = P1(I,K)*FI1(I,K)*AMIN1(NL1(I,K),PL1(I,K),SL1(I,K))/CCHL1(I,K)
!         P2(I,K) = P2(I,K)*FI2(I,K)*AMIN1(NL2(I,K),PL2(I,K),SL2(I,K))/CCHL2(I,K)
!         P3(I,K) = P3(I,K)*FI3(I,K)*AMIN1(NL3(I,K),PL3(I,K),SL3(I,K))/CCHL3(I,K)
!RGL below without P or Si limitation
         P1NNF(I,K) = P1(I,K)*FI1(I,K)*AMIN1(NL1NNF(I,K),PL1(I,K))/CCHL1(I,K)
         P1(I,K)    = P1(I,K)*FI1(I,K)*AMIN1(   NL1(I,K),PL1(I,K))/CCHL1(I,K)
         P2(I,K)    = P2(I,K)*FI2(I,K)*AMIN1(   NL2(I,K),PL2(I,K))/CCHL2(I,K)
!         P3(I,K) = P3(I,K)*FI3(I,K)*NL3(I,K)/CCHL3(I,K)
       ENDDO
     ENDDO  

   ENDIF
! RATE OF CHANGE DUE TO PRODUCTION, RESPIRATION, PREDATION

   DO K=1,KBM1
     DO I=1,MLOC
       
       NETP1 = (P1(I,K)*(1.-PRSP1)-BM1(I,K))*B1(I,K)  ! predation - metabolism
       NETP2 = (P2(I,K)*(1.-PRSP2)-BM2(I,K))*B2(I,K)
	   
!       NETP3 = (P3(I,K)*(1.-PRSP3)-BM3(I,K))*B3(I,K)
       GPP1 = P1(I,K)*B1(I,K)
       GPP2 = P2(I,K)*B2(I,K)
!       GPP3 = P3(I,K)*B3(I,K)

! PRIMARY PRODUCTION IN GM C/M**2/DAY

!       GPP(I,K) = (GPP1+GPP2+GPP3)*D(I)*DZ(K)
!       NPP(I,K) = (NETP1+NETP2+NETP3)*D(I)*DZ(K)
!
!not integrating it over total depth? looks like just calc ea depth  
!We should intergrate verticall here for GPP and NPP --Wen Long, now it is for each layer
!
        GPP(I,K) = (GPP1+GPP2)*D(I)*DZ(K)
        NPP(I,K) = (NETP1+NETP2)*D(I)*DZ(K)
! is this a rate or a total? looks like a rate. when added back to B1?       
       DTB1(I,K) = (NETP1-PR1(I,K)-B1SZ(I,K)-B1LZ(I,K))/86400.
	   
	   !dB1/Dt 			= 		(NEPT1/86400) - PR1/86400
	   
	   !gC/m^3/sec 		=		gC/m^3/d/86400	- 
	   
	   !==> PR1 must have unit gC/m^3/day 
	   
	   
       DTB2(I,K) = (NETP2-PR2(I,K)-B2SZ(I,K)-B2LZ(I,K))/86400.
!       DTB3(I,K) = (NETP3-PR3(I,K)-B3SZ(I,K)-B3LZ(I,K))/86400.

     ENDDO
   ENDDO 
      total_netPP(:)=0.0 
      DO K=1,KBM1
       DO I=1,MLOC
         total_netPP(I)=total_netPP(I)+NPP(I,K)
       ENDDO
      ENDDO
! NITROGEN FIXATION IN SURFACE BOXES (G N/G CHL/DAY)

   DO I=1,MLOC
     SNFIX(I) = (P1(I,1)-P1NNF(I,1))*ANC1*CCHL1(I,1)
   ENDDO

! RATE OF CHANGE DUE TO SETTLING
   DO I=1,MLOC
! 14 SEPT 2009 - KURT GLAESEMANN CHANGED X/Y/86400 TO X/(Y*86400)
     DTB1(I,1) = DTB1(I,1)-WS1(I,1)*B1(I,1)/(D(I)*DZ(1)*86400.)
     DTB2(I,1) = DTB2(I,1)-WS2(I,1)*B2(I,1)/(D(I)*DZ(1)*86400.)
     DTB3(I,1) = DTB3(I,1)-WS3(I,1)*B3(I,1)/(D(I)*DZ(1)*86400.)
   ENDDO          
!   DO I=1,MLOC
!     DTB1(I,1) = DTB1(I,1)-WS1(I,1)*B1(I,1)/(D(I)*DZ(1))/86400.
!     DTB2(I,1) = DTB2(I,1)-WS2(I,1)*B2(I,1)/(D(I)*DZ(1))/86400.
!!     DTB3(I,1) = DTB3(I,1)-WS3(I,1)*B3(I,1)/(D(I)*DZ(1))/86400.
!   ENDDO

   DO K=2,KBM1
     DO I=1,MLOC
! 14 SEPT 2009 - KURT GLAESEMANN CHANGED X/Y/86400 TO X/(Y*86400)
       DTB1(I,K) = DTB1(I,K)+(WS1(I,K-1)*B1(I,K-1)-WS1(I,K)*B1(I,K))  &
                  /(D(I)*DZ(K)*86400.)
       DTB2(I,K) = DTB2(I,K)+(WS2(I,K-1)*B2(I,K-1)-WS2(I,K)*B2(I,K))  &
                  /(D(I)*DZ(K)*86400.)                                    
       DTB3(I,K) = DTB3(I,K)+(WS3(I,K-1)*B3(I,K-1)-WS3(I,K)*B3(I,K))  &
                  /(D(I)*DZ(K)*86400.)
     ENDDO
   ENDDO  
!   DO K=2,KBM1
!     DO I=1,MLOC
!       DTB1(I,K) = DTB1(I,K)+(WS1(I,K-1)*B1(I,K-1)-WS1(I,K)*B1(I,K))  &
!                  /(D(I)*DZ(K))/86400.
!       DTB2(I,K) = DTB2(I,K)+(WS2(I,K-1)*B2(I,K-1)-WS2(I,K)*B2(I,K))  &
!                  /(D(I)*DZ(K))/86400.
!!       DTB3(I,K) = DTB3(I,K)+(WS3(I,K-1)*B3(I,K-1)-WS3(I,K)*B3(I,K))  &
!!                  /(D(I)*DZ(K))/86400.
!     ENDDO
!   ENDDO  

! SETTLING FLUX FOR MASS BALANCE

   DO K=1,KBM1
     DO I=1,MLOC
       FLXS1(I,K) = WS1(I,K)*B1(I,K)*V2(I,K)/(D(I)*DZ(K)*86400.)
       FLXS2(I,K) = WS2(I,K)*B2(I,K)*V2(I,K)/(D(I)*DZ(K)*86400.)
!       FLXS3(I,K) = WS3(I,K)*B3(I,K)*V2(I,K)/(D(I)*DZ(K)*86400.)
     ENDDO
   ENDDO  
        
!  RESUSPENSION

   IF (SEDIMENT_CALC) THEN 
! 14 SEPT 2009 - KURT GLAESEMANN MAKE THIS SIMPLER LOOP
   
   !WLong and Laura Bianucci comments:
   !     resuspension rate should be (WS1-WS1NET ) in m/d 
   !  and resuspension is counter-acting on the sinking of WS1 
   ! so the net sinking is WS1NET = WS1 - (WS1-WS1NET) , i.e. sinking - resuspension 
   ! 
   
     	DO I=1,MLOC
		  XL(I)=1.0/(D(I)*DZ(KBM1)*86400)
		  DTB1(I,KBM1) = DTB1(I,KBM1)+(WS1(I,KBM1)-WS1NET(I))*B1(I,KBM1)*XL(I)
		  DTB2(I,KBM1) = DTB2(I,KBM1)+(WS2(I,KBM1)-WS2NET(I))*B2(I,KBM1)*XL(I)
		  DTB3(I,KBM1) = DTB3(I,KBM1)+(WS3(I,KBM1)-WS3NET(I))*B3(I,KBM1)*XL(I)
		ENDDO
    

!     TREC=1.0/86400.
!     DO I=1,MLOC 
!       XL(I)=TREC/(D(I)*DZ(KBM1))     
!     ENDDO
!
!     DO I=1,MLOC 
!       DTB1(I,KBM1) = DTB1(I,KBM1)+(WS1(I,KBM1)-WS1NET(I))*B1(I,KBM1)*XL(I)
!     ENDDO
!     DO I=1,MLOC 
!       DTB2(I,KBM1) = DTB2(I,KBM1)+(WS2(I,KBM1)-WS2NET(I))*B2(I,KBM1)*XL(I)
!     ENDDO
!!     DO I=1,MLOC 
!!       DTB3(I,KBM1) = DTB3(I,KBM1)+(WS3(I,KBM1)-WS3NET(I))*B3(I,KBM1)*XL(I)
!!     ENDDO
   ENDIF 

   DEALLOCATE(XL,STOX1,STOX2,STOX3,NL1NNF,ITBOT)


   RETURN
   END SUBROUTINE ALGAE
   
   
   
   FUNCTION GET_FT(TEMPVAL,TMP,KTG1,KTG2)     !m/d  
   REAL(SP) ::     TEMPVAL,TMP,KTG1,KTG2             !temperature (degC)
   REAL(SP) ::    GET_FT       
       IF (TEMPVAL < TMP) THEN                  !shape is different for < TMP1 and > TMP1
         GET_FT = EXP(-KTG1*(TEMPVAL-TMP)**2)  !TMP1 is optimal temperature 
       ELSE                                    !KTG11 is the shape parameter (bell shape) for T < TMP1
         GET_FT = EXP(-KTG2*(TMP-TEMPVAL)**2)  !KTG12 is for T > TMP1
       ENDIF
   RETURN
   END FUNCTION GET_FT
   
   FUNCTION GET_FTBM(TEMPVAL,TR,KTB)     !m/d  
   REAL(SP) ::     TEMPVAL,TR,KTB             !temperature (degC)
   REAL(SP) ::    GET_FTBM      
       GET_FTBM=EXP(KTB*(TEMPVAL-TR))
   RETURN
   END FUNCTION GET_FTBM
   
   FUNCTION GET_FTPR(TEMPVAL,TRPR,KTPR)     !m/d  
   REAL(SP) ::     TEMPVAL,TRPR,KTPR             !temperature (degC)
   REAL(SP) ::    GET_FTPR    
       GET_FTPR=EXP(KTPR*(TEMPVAL-TRPR))
   RETURN
   END FUNCTION GET_FTPR
   
   !FUNCTION GET_FT1(TEMPVAL)     !m/d  
   !REAL(SP) ::     TEMPVAL             !temperature (degC)
   !REAL(SP) ::    GET_FT1       
   !    IF (TEMPVAL < TMP1) THEN                  !shape is different for < TMP1 and > TMP1
   !      GET_FT1 = EXP(-KTG11*(TEMPVAL-TMP1)**2)  !TMP1 is optimal temperature 
   !    ELSE                                    !KTG11 is the shape parameter (bell shape) for T < TMP1
   !      GET_FT1 = EXP(-KTG12*(TMP1-TEMPVAL)**2)  !KTG12 is for T > TMP1
   !    ENDIF
   !RETURN
   !END FUNCTION GET_FT1
   
   !FUNCTION GET_FTBM1(TEMPVAL)     !m/d  
   !REAL(SP) ::     TEMPVAL             !temperature (degC)
   !REAL(SP) ::    GET_FTBM1      
   !    GET_FTBM1=EXP(KTB1*(TEMPVAL-TR1))
   !RETURN
   !END FUNCTION GET_FTBM1
   
   
   !FUNCTION GET_FT2(TEMPVAL)     !m/d  
   !REAL(SP) ::     TEMPVAL             !temperature (degC)
   !REAL(SP) ::    GET_FT2       
   !    IF (TEMPVAL < TMP2) THEN                  !shape is different for < TMP1 and > TMP1
   !      GET_FT2 = EXP(-KTG21*(TEMPVAL-TMP2)**2)  !TMP1 is optimal temperature 
   !    ELSE                                    !KTG11 is the shape parameter (bell shape) for T < TMP1
   !      GET_FT2 = EXP(-KTG22*(TMP2-TEMPVAL)**2)  !KTG12 is for T > TMP1
   !    ENDIF
   !RETURN
   !END FUNCTION GET_FT2
   
   !FUNCTION GET_FTBM2(TEMPVAL)     !m/d  
   !REAL(SP) ::     TEMPVAL             !temperature (degC)
   !REAL(SP) ::    GET_FTBM2      
   !    GET_FTBM2=EXP(KTB2*(TEMPVAL-TR2))
   !RETURN
   !END FUNCTION GET_FTBM2
   
   
   !FUNCTION GET_FT3(TEMPVAL)     !m/d  
   !REAL(SP) ::     TEMPVAL             !temperature (degC)
   !REAL(SP) ::    GET_FT3       
   !    IF (TEMPVAL < TMP3) THEN                  !shape is different for < TMP1 and > TMP1
   !      GET_FT3 = EXP(-KTG31*(TEMPVAL-TMP3)**2)  !TMP1 is optimal temperature 
   !    ELSE                                    !KTG11 is the shape parameter (bell shape) for T < TMP1
   !      GET_FT3 = EXP(-KTG32*(TMP3-TEMPVAL)**2)  !KTG12 is for T > TMP1
   !    ENDIF
   !RETURN
   !END FUNCTION GET_FT3
   
   !FUNCTION GET_FTBM3(TEMPVAL)     !m/d  
   !REAL(SP) ::     TEMPVAL             !temperature (degC)
   !REAL(SP) ::    GET_FTBM3      
   !    GET_FTBM3=EXP(KTB3*(TEMPVAL-TR3))
   !RETURN
   !END FUNCTION GET_FTBM3
   
   !FUNCTION GET_FTPR(TEMPVAL)     !m/d  
   !REAL(SP) ::     TEMPVAL             !temperature (degC)
   !REAL(SP) ::    GET_FTPR    
   !    GET_FTPR=EXP(KTPR*(TEMPVAL-TRPR))
   !RETURN
   !END FUNCTION GET_FTPR
   
   
   
   END MODULE MOD_ALGAL
