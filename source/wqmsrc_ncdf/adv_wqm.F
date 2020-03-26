!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!=============================================================================!

!Subroutine ADV_WQM()


   SUBROUTINE ADV_WQM  
   USE MOD_LIMS, ONLY: KB, MYID, NPROCS, NCV,NUMPNT &
				  ,MLOC, MTLOC, NTLOC, KBM1, NCV_I
   USE MOD_PREC, ONLY : SP
#if defined(XFLUX_FLAG)
   USE MOD_OUTPUT_FLUX, ONLY: NEL_LOC, I_EDGELINE_ELE_EDGE_IA, &
								  I_EDGELINE_ELE_EDGE_IB, &
								  EDGELINE_FLUX, I_EDGELINE_FLUX_SIGN
#endif								  
   USE MOD_TGE, ONLY : 		&!
				 ISONB		&!
				,NTSN		&!
				,NBSN		&!
				,NIEC		&!
				,NTRG		&!
				,XIJE		&!
				,YIJE		&!
				,DLTXE		&!
				,DLTYE   &!
				,NTVE, NBVE, NV! Added by B Clark June 2018

   USE MOD_HYDROVARS, ONLY: &
   		!GRAV	&		!
		!,PI	&		!
		!,PI2	&		!
		ZERO	&		!
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
		,VX	&		!X-COORD AT GRID POINT
		,VY	&		!Y-COORD AT GRID POINT
		!,ART	&		!AREA OF ELEMENT
		,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		,ART2	&		!AREA OF ELEMENTS AROUND NODE
		!,NV	&		!NODE NUMBERING FOR ELEMENTS
		!,NBE	&		!INDICES OF ELMNT NEIGHBORS
		! ,NTVE	&		!! Uncommented by B Clark June 2018 for calculated average nodal velocity for sheer stress on bottom
		!,NTSN	&		!
		!,ISONB	&		!NODE MARKER = 0,1,2 
		!,ISBC	&		!
		!,ISBCE	&		!
		!,IEC	&		!
		!,IENODE &		!
!		,NBSN	&		!
!		,NIEC	&		!
!		,NTRG	&		!
	!	,NBVE	&		! ! Uncommented by B Clark June 2018 for calculated average nodal velocity for sheer stress on bottom
		!,NBVT	&		!
		!,LISBCE_1	&	!LIST OF ELEMENTS WITH ISBCE=1
		!,LISBCE_2	&	!LIST OF ELEMENTS WITH ISBCE=2
		!,LISBCE_3	&	!LIST OF ELEMENTS WITH ISBCE=3
		!,DLTXC	&		!
		!,DLTYC	&		!
		!,DLTXYC	&	!
		!,DLTXE	&		!
!		,DLTYE	&		!
		!,DLTXYE	&	!
		!,SITAC	&		!
		!,SITAE	&		!
		!,XIJC	&		!
		!,YIJC	&		!
!		,XIJE	&		!
!		,YIJE	&		!
		!,EPOR	&		!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
		!,IBCGEO	&	!LOCAL GEOSTROPHIC FRICTION CORRECTION NODES

		!,Z	&			!SIGMA COORDINATE VALUE 
		!,ZZ	&		!INTRA LEVEL SIGMA VALUE
		 ,DZ	&		!DELTA-SIGMA VALUE
		!,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
		!,H1	&		!BATHYMETRIC DEPTH 
		!,H	&			!BATHYMETRIC DEPTH 
		,D	&			!CURRENT DEPTH 
		,DT	&		!DEPTH AT PREVIOUS TIME STEP
		,DT1	&		!DEPTH AT PREVIOUS TIME STEP
		!,EL	&		!CURRENT SURFACE ELEVATION
		!,ET	&		!SURFACE ELEVATION AT PREVIOUS TIME STEP
		,DTFA	&		!ADJUSTED DEPTH FOR MASS CONSERVATION
		,UU	&		!X-VELOCITY
		,VV	&		!Y-VELOCITY
		!,WWL & ! updward velocity
		!,UUT	&		!X-VELOCITY FROM PREVIOUS TIMESTEP
		!,VVT	&		!Y-VELOCITY FROM PREVIOUS TIMESTEP
		!,WWT	&		!Z-VELOCITY FROM PREVIOUS TIMESTEP
		!,WTST	&		!Vertical velocity in sigma from PREVIOUS TIMESTEP
		!,UARD_OBCNT	&!tykim
		!,XFLUX_OBCT	&!tykim
		!,DTFAT	&		!tykim
		!,TT_T	&		!tykim
		!,SALTT	&		!tykim
		 ,WTS	&		!VERTICAL VELOCITY IN SIGMA SYSTEM
		!,UARD_OBCN	&	! tykim 
		!,XFLUX_OBC	&	! tykim 
		!,WTTS	&		!VERTICAL VELOCITY IN SIGMA SYSTEM 
		!,KH	&		!TURBULENT DIFFUSIVITY
		!,A1U	&		!
		!,A2U	&		!
		!,AWX	&		!
		!,AWY	&		!
		!,AW0	&		!
		,VISCOFH	&	!
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
		,THOUR	!&		!
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
			!SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
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
   USE MOD_WQM, ONLY : NAC,AC,C2,C1,C2F,DLT,WQ_DT,DTM,NCP,XYDF,ELTMS,XFLUX,XFLUX_ADV, ShearStress,RHO1, WW_OUT ! Added by B Clark June 2018 
   USE MOD_SED_DOM, ONLY:  RoughHeight
   USE MOD_WQMINIT, ONLY : XYDFU
   USE MOD_BCS, ONLY : 		&
!         INODEQ		    &  !!LOCAL FRESH WATER INFLOW NODES
		INOPNT      	&  !!LOCAL NON-POINT SOURCE NODES
		,PNT_TM         &  !!TIME MAP FOR NON POINT SOURCE DATA
		,WVQDIST   		&  !!DISCHARGE VERTICAL DISTRIBUTION for point source
		,PWQDIS    		&  !!Current WATER QUALITY at non-point source node
		,PDWQDIS 		&  !!WATER QUALITY DATA at non-point source node
		,PQDIS       	&  !!Current DISCHARGE at non-point source node
		,PDQDIS      	   !!DISCHARGE at non-point source node

# if defined (MULTIPROCESSOR)
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
			NBN ,       &   !          	!!NUMBER OF BOUNDARY NODES
			!MX_MLT,     &   !          	!!MAX MULTIPLICITY OF BOUNDARY NODES
			!BN_LST, 	&	!	   		!!GLOBAL IDENTITY OF BOUNDARY NODES
			BN_LOC, 	&	!			!!LOCAL IDENTITY OF BOUNDARY NODES
			BN_MLT, 	&	!			!!MULTIPLICITY OF BOUNDARY NODES
			!BN_NEY, 	&	!			!!NODE OWNER LIST
			!NDE_ID, 	&	!			!! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
			!NHN,  		&	!			!!NUMBER OF HALO NODES
			!HN_LST!, 	&	!			!!LIST OF HALO NODES 
			!HN_OWN, 	&	!			!!PRIMARY OWNER OF HALO NODES
			!EC,			&	!
			!NC,			&	!
			BNC, 		&	!
			!EMAP,		&	!
			!NMAP, 		&	!
			!PNE,		&	!(:)        !!NUMBER OF ELEMENTS IN EACH PROC
			!PNN,		&	!(:)        !!NUMBER OF NODES IN EACH PROC
			!PNHE,		&	!(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
			!PNBN,		&	!(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
			!PMBM,		&	!(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
			!PNHN,		&	!(:)       !!NUMBER OF HALO NODES IN EACH PROC
			!FILELOCK,   & !
			!MYFILELOCK, & !
			!FILELOCKED , & !
			!NODE_MATCH, & !
			 NODE_MATCH_ALL  !, & !
			!EXCHANGE, & !
			!EXCHANGE_ALL, & !
			!GATHER, & !
			!IGATHER, & !
			!INIT_MPI_ENV,&!
			!INIT_MPI_ENV,&!
			!SORT,	& !
			!PREAD,	& !
			!PWRITE,& !
			!PPRINT, & !
			!LOCK, & !
			!UNLOCK, & !
			!GETLOC !
		 
# endif

   IMPLICIT NONE
!   REAL(SP), ALLOCATABLE, DIMENSION(:,:,:)  :: XFLUX,XFLUX_ADV
   REAL(SP), DIMENSION(MLOC)           :: PUPX,PUPY,PVPX,PVPY
   REAL(SP), DIMENSION(MLOC)           :: PFPX,PFPY,PFPXD,PFPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NTLOC))      :: DTIJ
   REAL(SP), DIMENSION(3*(NTLOC),KBM1) :: UVN
   REAL(SP) :: FFD,FF1,X11,Y11,X22,Y22,X33,Y33,XI,YI
   REAL(SP) :: DXA,DYA,DXB_A(NCV_I),DYB_A(NCV_I),DXB_B(NCV_I),DYB_B(NCV_I),FIJ1,FIJ2,UN,DYB,DXB
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP
   REAL(SP) :: FACT,FM1
   REAL(SP) :: TT,TTIME,STPOINT,STPOINT1,STPOINT2
   INTEGER  :: I,I1,I2,IA,IB,IC,J,J1,J2,JTMP,K,JJ,		&
				JCON,									&	!Index for constituents
				II,IX
   REAL(SP) :: C21MIN, C21MAX, C22MIN, C22MAX
! KURT GLAESEMANN, fix C2MEAN size
!  REAL(SP) :: C2MEAN(MLOC,KBM1,NCP)
  REAL(SP) :: C2MEAN(MTLOC,KBM1,NCP)
! KURT GLAESEMANN 10 SEPT 2009 - precalculate things
  LOGICAL XYDFUISON
  REAL(SP) :: TMP1(MLOC), TMP2(MLOC)

  INTEGER :: NJ1, NJ2, NJ3, EXTRA
!TYKIM for nonpoint sources
  INTEGER :: L1, L2, IERR
  REAL(SP)    :: UFACT, DFACT
    
  REAL(SP) :: Nodal_V, Nodal_U, USTAR
  REAL(SP) :: Kappa = 0.41 ! Von Karman's constants
!------------------------------------------------------------------------------
!  KURT GLAESEMANN 11 SEPT 2009 C2MEAN IS ZERO
!   C2MEAN = 0.0

!  KURT GLAESEMANN 17 SEPT 2009, make XFLUX allocated
!   ALLOCATE(XFLUX(0:MTLOC,KB,NCP))
!   ALLOCATE(XFLUX_ADV(0:MTLOC,KB,NCP))
!  KURT GLAESEMANN 10 SEPT 2009 - convert this to logic below
!   FACT = 0.0
!   FM1  = 1.0
!   IF(XYDFU == ' ON')THEN
!     FACT = 1.0
!     FM1  = 0.0
!   ENDIF
    XYDFUISON = .false.
    IF(XYDFU == ' ON') XYDFUISON = .true.

! KURT GLAESEMANN 10 SEPT 2009 - precalculate stuff
    DO I = 1,MLOC
      TMP1(I) = DLT/DT(I)/ART1(I)
      TMP2(I) = DT(I)/DTFA(I)
   ENDDO

! KURT GLAESEMANN 17 SEPT 2009 precalculate more stuff
   DO I=1,NCV_I
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     XI = 0.5_SP*(XIJE(I,1)+XIJE(I,2))  ! middle point of element edge x coords
     YI = 0.5_SP*(YIJE(I,1)+YIJE(I,2))   ! middle point of element edge y coords
     DXB_A(I)=XI-VX(IA)   ! distance from node to midpoint of element edge to the right (m)
     DYB_A(I)=YI-VY(IA)  
     DXB_B(I)=XI-VX(IB)
     DYB_B(I)=YI-VY(IB)
   ENDDO

!
!--Initialize Fluxes-----------------------------------------------------------
!
   XFLUX = 0.0
   XFLUX_ADV = 0.0
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------
!
   DO I=1,NCV ! loop over the internal control volumes
     I1=NTRG(I) ! The element associated with this control volume
     DTIJ(I)=DT1(I1)  ! Previous time step depth
     DO K=1,KBM1
       UVN(I,K)=VV(I1,K)*DLTXE(I) - UU(I1,K)*DLTYE(I)   ! normal velocity
     ENDDO
   ENDDO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------
!
   DO JCON=1,NAC
     II = AC(JCON)
!write(*,*) 'JCON = ',JCON
!write(*,*) 'II= ',II
     DO K=1,KBM1
       PFPX  = 0.0
       PFPY  = 0.0
! KURT GLAESEMANN 10 SEPT 2009 - C2MEAN IS ZERO
!      PFPXD = 0.0
!      PFPYD = 0.0

       DO I=1,MLOC
         DO J=1,NTSN(I)-1  ! total number of surrounding nodes
           I1=NBSN(I,J)  ! ID number of surrouding nodes
           I2=NBSN(I,J+1)
! KURT GLAESEMANN 10 SEPT 2009 - C2MEAN IS ZERO
!          FFD=0.5_SP*(C2(I1,K,II)+C2(I2,K,II)-C2MEAN(I1,K,II)-C2MEAN(I2,K,II))
           FF1=0.5_SP*(C2(I1,K,II)+C2(I2,K,II))

           PFPX(I)=PFPX(I)+FF1*(VY(I1)-VY(I2)) ! Get total average mixing between the nodes surrounding each node g / m^3 * m (Conc * dx) =  g / m^2
           PFPY(I)=PFPY(I)+FF1*(VX(I2)-VX(I1))
! KURT GLAESEMANN 10 SEPT 2009 - C2MEAN IS ZERO
!          PFPXD(I)=PFPXD(I)+FFD*(VY(I1)-VY(I2))
!          PFPYD(I)=PFPYD(I)+FFD*(VX(I2)-VX(I1))
		 ENDDO

         PFPX(I)=PFPX(I)/ART2(I)  ! di  !   ! total area of surrounding elements
         PFPY(I)=PFPY(I)/ART2(I)   ! g  / m^4
! KURT GLAESEMANN 10 SEPT 2009 - C2MEAN IS ZERO
!        PFPXD(I)=PFPXD(I)/ART2(I)
!        PFPYD(I)=PFPYD(I)/ART2(I)
       ENDDO

       DO I=1,MLOC
	   
         VISCOFF(I)=VISCOFH(I,K)
		 		 	 
       ENDDO
	   
	!    write(*,*)'NCV_I = ',NCV_I

	   DO I=1,NCV_I ! loop over total number of control volumes (local)
		
		 IA=NIEC(I,1) ! counting number of left and right elements
		 IB=NIEC(I,2)
		 
         UN=UVN(I,K)  ! normal velocity
		 
!  KURT GLAESEMANN 10 SEPT 2009 - convert FACT and FM1 to logic
         IF(XYDFUISON)THEN
            VISCOF=XYDF*((VISCOFF(IA)+VISCOFF(IB))*0.5)
         ELSE
            VISCOF=XYDF   
         ENDIF
!         VISCOF=XYDF*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5 + FM1)

! KURT GLAESEMANN 10 SEPT 2009 - C2MEAN IS ZERO
!        TXX=0.5*(PFPXD(IA)+PFPXD(IB))*VISCOF
!        TYY=0.5*(PFPYD(IA)+PFPYD(IB))*VISCOF
         TXX=0.5*(PFPX (IA)+PFPX (IB))*VISCOF  ! Take the average of both left and right nodes multiplied by the diffusion coeff to get a flux  ( g /m ^4 * m^2/s ---> g /m^2/s)
         TYY=0.5*(PFPY (IA)+PFPY (IB))*VISCOF

         FXX=-DTIJ(I)*TXX*DLTYE(I)  ! Previous time step depth * flux * y distance on this edge = g /s
         FYY= DTIJ(I)*TYY*DLTXE(I)

         if(UN .ge. 0) then ! If normal velocity is positive than use the left node  ---->
           ic = ib
           DYB=DYB_B(I)
           DXB=DXB_B(I)
         else
           ic = ia ! if normal velocity is negative than use the right node   <----
           DYB=DYB_A(I)  ! the distance bet
           DXB=DXB_A(I)
         endif
         C22MIN=C2(IC,K,II)
         C22MAX=C2(IC,K,II)
         do j=1,NTSN(IC)-1
            temp=C2(NBSN(IC,j),K,II)
            C22MIN=min(C22MIN,temp)
            c22max=max(c22max,temp)
         enddo
         FIJ2=C2(IC,K,II)+DXB*PFPX(IC)+DYB*PFPY(IC)
         IF(FIJ2 < C22MIN) FIJ2=C22MIN
         IF(FIJ2 > C22MAX) FIJ2=C22MAX
!caps above for minmax criteria
		 
		 EXFLUX=-UN*DTIJ(I)*FIJ2+FXX+FYY
         
! END KURT
		 !Below added by B Clark, NASA GSCF March 2020
		 ! to get xflux across
		 !control volume faces we care about
		 !Edgeline_IA and Edgline_IB are the two nodes we care about
		 !positive flux is out of the transect
		 !negative flux is into the transect
#if defined(XFLUX_FLAG)		 
      IF(NEL_LOC>0)THEN
		 DO IX=1,NEL_LOC

			IF(IA==I_EDGELINE_ELE_EDGE_IA(IX) .AND. & 
			  (IB) == I_EDGELINE_ELE_EDGE_IB(IX)) THEN

			   EDGELINE_FLUX(IX,K,II)=EXFLUX*i_edgeline_flux_sign(IX)

		    ENDIF
		 ENDDO
	 ENDIF
			 
#endif
         XFLUX(IA,K,II)=XFLUX(IA,K,II)+EXFLUX
         XFLUX(IB,K,II)=XFLUX(IB,K,II)-EXFLUX

         XFLUX_ADV(IA,K,II)=XFLUX_ADV(IA,K,II)+(EXFLUX-FXX-FYY)
		 XFLUX_ADV(IB,K,II)=XFLUX_ADV(IB,K,II)-(EXFLUX-FXX-FYY)
		 
	   ENDDO

     ENDDO
   ENDDO

!
!-Accumulate Fluxes at Boundary Nodes
!
   NODAL_V=0.0
   NODAL_U=0.0
   DO I = 1,MLOC
   
   		 ! B Clark add in calculation for bottom sheer stress to apply sediment resuspension
		 ! June 2018
		 
		 !Average Velocity
		 
!		 write(*,*)'Calculating nodal vel'
		
		 Nodal_V=SUM(VV(NBVE(I,1:NTVE(I)),KBM1))/NTVE(I) ! average V velocity
		 Nodal_U=SUM(UU(NBVE(I,1:NTVE(I)),KBM1))/NTVE(I) ! Average U velocity

!		 write(*,*)'Calculating Ustar'		 
		 
		 USTAR=(Nodal_V**2+Nodal_U**2)**0.5*kappa/LOG((D(I)*DZ(KBM1))/RoughHeight) ! ustar over the bed, common formula 
!		 		 write(*,*)'Calculating shear stress'
		 
		 
		 ShearStress(I)=RHO1(I,KBM1)*USTAR**2 ! sheer stress from density and ustar
		 
		 
   ENDDO
		 
		
# if defined (MULTIPROCESSOR)
   IF(PAR)THEN
! KURT GLAESEMANN 17 SEPT 2009 DO MULTIPLE VARIABLES AT ONCE
!    DO JCON=1,NAC
!      II=AC(JCON)
!      CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MTLOC,KB,MYID,NPROCS,XFLUX(:,:,II))
!    ENDDO 

#if(0)
     EXTRA = MOD(NAC,3)
     IF (EXTRA .eq. 1) THEN
       NJ1= AC(1)
       CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MTLOC,KB,MYID,NPROCS,XFLUX(:,:,NJ1))
     ELSEIF (EXTRA .eq. 2) THEN
       NJ1= AC(1)
       NJ2= AC(2)
       CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MTLOC,KB,MYID,NPROCS,XFLUX(:,:,NJ1),XFLUX(:,:,NJ2))
     ENDIF

     DO JCON=EXTRA+1,NAC-2,3
       NJ1= AC(JCON+0)
       NJ2= AC(JCON+1)
       NJ3= AC(JCON+2)
       CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MTLOC,KB,MYID,NPROCS,XFLUX(:,:,NJ1),XFLUX(:,:,NJ2),XFLUX(:,:,NJ3))
     ENDDO
#else
!  KURT GLAESEMANN 22 SEPT 2009 - DO ALL AT ONCE 
     CALL NODE_MATCH_ALL(0,NBN,BN_MLT,BN_LOC,BNC,MTLOC,KB,NCP,MYID,NPROCS,XFLUX)
#endif
   ENDIF  
# endif
!     DO JCON=1,NAC
!       II=AC(JCON)
!       DO K=1,KBM1
!          IF(NUMPNT > 0) THEN
!             DO J=1,NUMPNT
!                JJ=INOPNT(J)
!                XFLUX(JJ,K,II)=XFLUX_ADV(JJ,K,II)
!             ENDDO
!          ENDIF
!       ENDDO  
!     ENDDO

!TYKIM
!NON-POINT SOURCE INTERPOLATION
      IF(NUMPNT > 0) THEN
         ALLOCATE(PWQDIS(NUMPNT,NCP))  ;PWQDIS = ZERO
         ALLOCATE(PQDIS(NUMPNT))       ;PQDIS  = ZERO
         CALL BRACKET(PNT_TM,THOUR,L1,L2,DFACT,UFACT,IERR)
         PWQDIS(:,:)=UFACT*PDWQDIS(:,:,L1)+DFACT*PDWQDIS(:,:,L2)
         PQDIS(:)   =UFACT*PDQDIS(:,L1)   +DFACT*PDQDIS(:,L2)
         PQDIS(:)=PQDIS(:)*TANH(ELTMS/FLOAT(86400))
      ENDIF

   DO JCON=1,NAC
     II = AC(JCON)

!
!--Calculate the Vertical Terms------------------------------------------------
!
     DO K=1,KBM1
       DO I=1,MLOC
         IF(K == 1) THEN
           TEMP=-WTS(I,K+1)*(C2(I,K,II)*DZ(K+1)+C2(I,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1))
         ELSEIF(K == KBM1) THEN
           TEMP=WTS(I,K)*(C2(I,K,II)*DZ(K-1)+C2(I,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1))
         ELSE
           TEMP=WTS(I,K)*(C2(I,K,II)*DZ(K-1)+C2(I,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1))-  &
                WTS(I,K+1)*(C2(I,K,II)*DZ(K+1)+C2(I,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1))
         ENDIF

!--Total Fluxes ---------------------------------------------------------------
!
         IF(ISONB(I) == 2) THEN
           XFLUX(I,K,II)=TEMP*ART1(I)/DZ(K)
         ELSE
           XFLUX(I,K,II)=XFLUX(I,K,II)+TEMP*ART1(I)/DZ(K)
         ENDIF
! KURT GLAESEMANN 10 SEPT 2009 MERGE TWO LOOPS- rgl uncommented
! to allow nonpoint
       ENDDO
     ENDDO   
	 
	 ! ! Get average vertical velocity at element center
	 
	 DO K = 1,KBM1
	   DO I=1,NTLOC 
		 WW_OUT(I,K)=(WTS(NV(I,1),K)+WTS(NV(I,2),K)+WTS(NV(I,3),K))/3
	   END DO
     END DO
		 
	 

!TYKIM
!===============================================================================
!POINT SOURCE FLUX
!===============================================================================
    IF(NUMPNT > 0) THEN
      DO J=1,NUMPNT
        JJ=INOPNT(J)
        DO K=1,KBM1
         STPOINT = PWQDIS(J,II)
!         STPOINT = C2(J,K,II)
         XFLUX(JJ,K,II)=XFLUX(JJ,K,II)-PQDIS(J)*WVQDIST(J,K)*STPOINT/DZ(K)
        ENDDO
      ENDDO
    ENDIF
!
!
!--Update Water Quality Variables--------------------------------
!
     DO I = 1,MLOC
         DO K = 1, KBM1
! KURT GLAESEMANN 10 SEPT 2009 PRECALCULATE A BUNCH OF STUFF
!           C1(I,K,II)=(C2(I,K,II)-XFLUX(I,K,II)/ART1(I)*(DLT/DT(I)))*   &
!                         (DT(I)/DTFA(I))+DTM(I,K,II)*DLT
!Here is where wq variables updated! Hooray
! RGL added below to to water quality sub-cycling
!      IF(MOD(ELTMS,WQ_DT).eq.0) THEN
          
			C1(I,K,II)=(C2(I,K,II)-XFLUX(I,K,II)*TMP1(I))*TMP2(I)!+DTM(I,K,II)*DLT ! See below, DTM is now added after flux limiter
		    !wen long:
			!DTM should have unit of concentration/sec	
			!DLT should have unit of sec, which is the time step
			!cocnentration for each constituent is different,in general g/m^3
			!
		  
!      ELSE
!          	C1(I,K,II)=(C2(I,K,II)-XFLUX(I,K,II)*TMP1(I))*TMP2(I)
!      ENDIF
          
!ORI        C1(I,K,II)=(C2(I,K,II)-XFLUX(I,K,II)/ART1(I)*(DLT/DT(I)))*   &
!                         (DT(I)/DTFA(I))+DTM(I,K,II)*DLT
				! IF(II==27 .or. II==35)then
				! 	if(K==3)then
				! 		if(NGID(I)==7849) THEN
				! 			write(*,*)'II = ',II
				! 			write(*,*)'C1 before assignment = ',C1(I,K,II)
				! 			write(*,*)'XFLUX before assignment = ',XFLUX(I,K,II)

				! 		endif 
				! 	endif
				! endif
		    ! C2F(I,K,II)=C1(I,K,II)
			  C2F(I,K,II)=MAX(C1(I,K,II),1.e-15)
			! IF(II==27 .or. II==35)then
			! 	if(K==3)then
			! 		if(NGID(I)==7849) THEN
			! 			write(*,*)'II = ',II
			! 	        write(*,*)'C2F after assignment = ',C2F(I,K,II)
			! 		endif 
			! 	endif
			! endif
            !   C2F(I,K,II) = MAX( C1(I,K,II),  10E-15 )
			! if(C2F(I,K,II)<=0.)THEN
			! 	C2F(I,K,II) = 1E-10;
			! endif
        !    IF(II.EQ.4.OR.II.EQ.5.OR.II.EQ.6)C2F(I,K,II)=MAX( C1(I,K,II), 0.003)  !make sure there are seeds for ALG1,2 and 3

         ENDDO
     ENDDO

   ENDDO
   
   
     CALL  FCT_NUT  ! Moved here by B Clark to fix problem with reaction/ source sink terms being smoothed out by FCT_NUT
                  ! Now, the change in concentration is updated after the flux control is put in place
   
    DO JCON=1,NAC
     II = AC(JCON)
       DO K=1,KBM1
		 DO I=1,MLOC
			! IF(II==27 .or. II==35)then
			! 	if(K==3)then
			! 		if(NGID(I)==7849) THEN
			! 			write(*,*)'II = ',II
			! 			write(*,*)'C2F before DTM = ',C2F(I,K,II)
			! 			write(*,*)'DTM = ',DTM(I,K,II)
			! 		endif 
			! 	endif
			! endif
			 C2F(I,K,II)=C2F(I,K,II) + DTM(I,K,II)*DLT
			   C2F(I,K,II) = DMAX1 ( C2F(I,K,II), 1.e-15 )
			! if(C2F(I,K,II)<=0.)THEN
			! 	C2F(I,K,II) = 1E-10;
			! endif
			! IF(II==27 .or. II==35)then
			! 	if(K==3)then
			! 		if(NGID(I)==7849) THEN
			! 			write(*,*)'II = ',II
            !             write(*,*)'C2F AFTER DTM = ',C2F(I,K,II)
			! 		endif 
			! 	endif
			! endif
            !       C2F(I,K,II) = MAX( C2F(I,K,II), 10E-15 )
			! IF(II.EQ.4.OR.II.EQ.5.OR.II.EQ.6)C2F(I,K,II)=MAX( C2F(I,K,II), 0.003)  !make sure there are seeds for ALG1,2 and 3
               !write(*,*)'C2F after= ',C2F(I,K,II)

         ENDDO
       ENDDO
     ENDDO
!	 IF(MSR)THEN
!		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
!		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," FCTNUT"
!	 ENDIF
	 
!	 IF(MSR)THEN
!		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
!		NPROFILE=NPROFILE+1
!	 ENDIF
   IF(NUMPNT > 0) DEALLOCATE(PQDIS, PWQDIS) 
!   DEALLOCATE (XFLUX)
!   DEALLOCATE (XFLUX_ADV)
   RETURN
   END SUBROUTINE ADV_WQM
!==============================================================================!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 
