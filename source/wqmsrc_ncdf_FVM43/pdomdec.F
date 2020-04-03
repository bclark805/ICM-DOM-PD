!Subroutine PDOMDEC()

!==============================================================================|
!  SET UP LOCAL PHYSICAL DOMAIN (CONNECTIVITY/MESH)                            |
!==============================================================================|

   SUBROUTINE PDOMDEC

!==============================================================================!
	USE MOD_LIMS, ONLY: NLOC, NTLOC, MLOC, MTLOC &
				, KBM1, NPROCS
	USE MOD_PREC, ONLY: SP			
	USE MOD_TGE, ONLY: NV
	USE MOD_HYDROVARS, ONLY: &
   		!GRAV	&		!
		!,PI	&		!
		!,PI2	&		!
		ZERO	&		!
		!,ONE_THIRD	&	!
		  ,NVG	&		!
		 ,XG	&		!GLOBAL X-COORD AT NODE 
		 ,YG	&		!GLOBAL X-COORD AT NODE 
		,HG	&		!GLOBAL DEPTH AT NODE 
		,XCG	&		!GLOBAL X-COORD AT FACE CENTER 
		,YCG	&		!GLOBAL X-COORD AT FACE CENTER 
		,VXMIN	&		!
		,VYMIN	&		!
		,VXMAX	&		!
		,VYMAX	&		!
		,XC	&		!X-COORD AT FACE CENTER 
		,YC	&		!Y-COORD AT FACE CENTER
		,VX	&		!X-COORD AT GRID POINT
		,VY	&		!Y-COORD AT GRID POINT
		!,ART	&		!AREA OF ELEMENT
		!,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		!,ART2	&		!AREA OF ELEMENTS AROUND NODE
!		 ,NV	&		!NODE NUMBERING FOR ELEMENTS
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
		!,DZ	&		!DELTA-SIGMA VALUE
		!,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
		,H1	&		!BATHYMETRIC DEPTH 
		,H	!&			!BATHYMETRIC DEPTH 
		!,D	&			!CURRENT DEPTH 
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
			!,MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			,PAR        	&           !!TRUE IF MULTIPROCESSOR RUN
			!,CASENAME  	&   		!!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
			!,CASETITLE  	&  			!!CASE TITLE                                 
			,HMAX       	&  			!!GLOBAL MAXIMUM DEPTH
			,HMIN   !    	&  			!!GLOBAL MINIMUM DEPTH
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
		USE MOD_PAR, ONLY:   &	  
			!EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			!ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			NLID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			!ELID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			NLID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,	&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
			NHE,        &   !			!NUMBER OF HALO ELEMENTS
			HE_LST,		&	!(:)       	!!GLOBAL IDENTITIES OF HALO ELEMENTS
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
	
		 USE MOD_SIZES, ONLY :        &	!
			!NCP,            &  !
			!NQFP,           &  !
            !NHQP,           &  !
			!NS1P,           &  !
            !NS2P,           &  !
            !NS3P,           &  !
            !NBCP,           &  !
            !NDP,            &  !
            !NFLP,           &  !
            !NOIP,           &  !
            !NSSFP,			 &  !
            MGL,            &  !
			NGL!,            &  !
			!OBCGL,          &  !
			!NOBTY  
  
   IMPLICIT NONE
   INTEGER I,EGL,J,IERR,I1,I2  

!==============================================================================|
!  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
!  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
!  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
!==============================================================================|
   !WLong moved this to HYDRO_ALLOC()
   !ALLOCATE(NV(0:NTLOC,4));    NV = 0  !!NODE NUMBERING FOR ELEMENTS
   IF(SERIAL) NV = NVG

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,NLOC
       EGL = EGID(I)
       NV(I,1:4) = NLID(NVG(EGID(I),1:4))
     ENDDO
     DO I=1,NHE
       EGL = HE_LST(I)
       NV(I+NLOC,1:4) = NLID_X(NVG(EGL,1:4))
     ENDDO
   ENDIF
#  endif

!==============================================================================|
!   SET UP LOCAL MESH (HORIZONTAL COORDINATES)                                 |
!==============================================================================|

!--------------CALCULATE GLOBAL MINIMUMS AND MAXIMUMS--------------------------!

   VXMIN = MINVAL(XG(1:MGL)) ; VXMAX = MAXVAL(XG(1:MGL))
   VYMIN = MINVAL(YG(1:MGL)) ; VYMAX = MAXVAL(YG(1:MGL))

!--------------SHIFT GRID TO UPPER RIGHT CARTESIAN-----------------------------!

   XG = XG - VXMIN
   YG = YG - VYMIN
   XG(0) = 0.0_SP ; YG(0) = 0.0_SP

!--------------CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES----------------!
	!WLong moved this to HYDRO_GEOM_ALLOC
    !ALLOCATE(XCG(0:NGL),YCG(0:NGL)) ; XCG = 0.0_SP ; YCG = 0.0_SP
   
   DO I=1,NGL   
     XCG(I)  = (XG(NVG(I,1)) + XG(NVG(I,2)) + XG(NVG(I,3)))/3.0_SP
     YCG(I)  = (YG(NVG(I,1)) + YG(NVG(I,2)) + YG(NVG(I,3)))/3.0_SP
   ENDDO

   XCG(0) = 0.0_SP ; YCG(0) = 0.0_SP


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
   !WLong moved this to HYDO_ALLOC()
   !ALLOCATE(VX(0:MTLOC));        VX   = ZERO   !!X-COORD AT GRID POINT   
   !ALLOCATE(VY(0:MTLOC));        VY   = ZERO   !!X-COORD AT GRID POINT   
   IF(SERIAL)THEN
     VX = XG
     VY = YG
   ENDIF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,MLOC
       VX(I) = XG(NGID(I))
       VY(I) = YG(NGID(I))
     ENDDO

     DO I=1,NHN
       VX(I+MLOC) = XG(HN_LST(I))
       VY(I+MLOC) = YG(HN_LST(I))
     ENDDO
   ENDIF
# endif

!==============================================================================|
!   SET UP LOCAL MESH (BATHYMETRIC DEPTH)                                      |
!==============================================================================|

!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
   !WLong moved this to mod_hydrovars.F HYDRO_ALLOC
   !ALLOCATE(H(0:MTLOC));       H = ZERO       !!BATHYMETRIC DEPTH
   IF(SERIAL) H = HG

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,MLOC
       H(I)   = HG(NGID(I))
     ENDDO
     DO I=1,NHN
       H(I+MLOC) = HG(HN_LST(I))
     ENDDO
   ENDIF
#  endif

!--------------CALCULATE EXTREMUMS---------------------------------------------!

   HMAX = MAXVAL(ABS(HG(1:MGL)))
   HMIN = MINVAL(HG(1:MGL))

!==============================================================================|
!   COMPUTE FACE CENTER VALUES FOR GRID, DEPTH, AND CORIOLIS PARAMETER         |
!==============================================================================|

   !WLong moved these to HYDRO_ALLOC
   !ALLOCATE(XC(0:NTLOC))       ;XC   = ZERO   !!X-COORD AT FACE CENTER
   !ALLOCATE(YC(0:NTLOC))       ;YC   = ZERO   !!X-COORD AT FACE CENTER
   !ALLOCATE(H1(0:NTLOC))       ;H1   = ZERO   !!BATHYMETRIC DEPTH
   
   
   DO I=1,NTLOC
     XC(I)  = (VX(NV(I,1)) + VX(NV(I,2)) + VX(NV(I,3)))/3.0_SP
     YC(I)  = (VY(NV(I,1)) + VY(NV(I,2)) + VY(NV(I,3)))/3.0_SP
     H1(I)  = SUM( H(NV(I,1:3)))/3.0_SP
   ENDDO

   RETURN
   END SUBROUTINE PDOMDEC
!==============================================================================|
