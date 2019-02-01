!==============================================================================|
!    This subroutine is used to calculate the area of individual               !
!    triangle based on the three vertex coordinates and also calculate         !
!    the sigma-surface area of individual control volume consisted of          !
!    triangles with a common node point                                        !
!                                          !
! calculates: art(ntri)   = area of element (triangle)                         !
! calculates: art1(nnode) = area of interior cv (for node value integration)   !
! calculates: art2(nnode) = sum area of all cells around node                  !
!==============================================================================|

!subroutine CELL_AREA()

	SUBROUTINE CELL_AREA 

!==============================================================================!
	USE MOD_LIMS, ONLY: NTLOC, MTLOC, NLOC, MLOC, MYID, NPROCS
	USE MOD_PREC, ONLY: SP
	
    USE MOD_TGE, ONLY : 	&!
				ISONB		&!
				,NTVE		&!
				,NBVE		&!		
				,NBVT		&!
				,NV			

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
		,XC	&		!X-COORD AT FACE CENTER 
		,YC	&		!Y-COORD AT FACE CENTER
		,VX	&		!X-COORD AT GRID POINT
		,VY	&		!Y-COORD AT GRID POINT
		,ART	&		!AREA OF ELEMENT
		,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		,ART2 !	&		!AREA OF ELEMENTS AROUND NODE
!		,NV	&		!NODE NUMBERING FOR ELEMENTS
		!,NBE	&		!INDICES OF ELMNT NEIGHBORS
!		,NTVE	&		!
		!,NTSN	&		!
!		,ISONB	&		!NODE MARKER = 0,1,2 
		!,ISBC	&		!
		!,ISBCE	&		!
		!,IEC	&		!
		!,IENODE &		!
		!,NBSN	&		!
		!,NIEC	&		!
		!,NTRG	&		!
!		,NBVE	&		!
!		,NBVT!	&		!
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
		!,H1	&		!BATHYMETRIC DEPTH 
		!,H	&			!BATHYMETRIC DEPTH 
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
			!SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			!,MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			PAR        !	&           !!TRUE IF MULTIPROCESSOR RUN
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
			 NBN,       &   !          	!!NUMBER OF BOUNDARY NODES
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
			!MYFILELOCK, &!
			!FILELOCKED, &!
			 NODE_MATCH!, &!
			!NODE_MATCH_ALL, &!
			!EXCHANGE, &!
			!EXCHANGE_ALL, &!
			!GATHER, &!
			!IGATHER, &!
			!INIT_MPI_ENV,&!			 
			!SORT,	&!
			!PREAD,	&!
			!PWRITE,&!
			!PPRINT, &!
			!LOCK, &!
			!UNLOCK, &!
			!GETLOC
			 
# endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE :: XX(:),YY(:) 
   REAL(SP) :: ARTMAX,ARTTOT,ARTMIN
   INTEGER  :: I,J,II,J1,J2,MAX_NBRE
!==============================================================================!

!---------------INITIALIZE ARRAYS----------------------------------------------!

   !WLong moved these to HYDRO_ALLOC()
   !ALLOCATE(ART(0:NTLOC))           ;ART  = ZERO   !!AREA OF ELEMENT
   !ALLOCATE(ART1(0:MTLOC))          ;ART1 = ZERO   !!AREA OF NODE-BASE CONTROl VOLUME
   !ALLOCATE(ART2(0:MTLOC))            ;ART2 = ZERO   !!AREA OF ELEMENTSAROUND NODE
   
   MAX_NBRE = MAXVAL(NTVE)+1
   ALLOCATE(XX(2*MAX_NBRE+1),YY(2*MAX_NBRE+1))
   XX = 0.0_SP ; YY = 0.0_SP

!---------------COMPUTE AREA OF TRIANGLES USING CROSS PRODUCT------------------!
  
    DO I=1,NTLOC   !surface area of all elements
    ART(I) = (VX(NV(I,2)) - VX(NV(I,1))) * (VY(NV(I,3)) - VY(NV(I,1))) - &
             (VX(NV(I,3)) - VX(NV(I,1))) * (VY(NV(I,2)) - VY(NV(I,1)))
    ENDDO
   ART    = ABS(.5_SP*ART)

!---------------COMPUTE MESH STATISTICS----------------------------------------!

   ARTMIN = MINVAL(ART(1:NLOC))
   ARTMAX = MAXVAL(ART(1:NLOC))
   ARTTOT =    SUM(ART(1:NLOC))

!-------COMPUTE CONTROL VOLUME ART1: CV FOR FLUXES OF NODAL BASED VALUES-------!

   DO I=1,MLOC
     IF(ISONB(I) == 0) THEN
       DO J=1,NTVE(I)
         II=NBVE(I,J)
         J1=NBVT(I,J)
         J2=J1+1-INT((J1+1)/4)*3
         XX(2*J-1)=(VX(NV(II,J1))+VX(NV(II,J2)))*0.5_SP-VX(I)
         YY(2*J-1)=(VY(NV(II,J1))+VY(NV(II,J2)))*0.5_SP-VY(I)
         XX(2*J)=XC(II)-VX(I)
         YY(2*J)=YC(II)-VY(I)
       ENDDO
       XX(2*NTVE(I)+1)=XX(1)
       YY(2*NTVE(I)+1)=YY(1)

       DO J=1,2*NTVE(I)
          ART1(I)=ART1(I)+0.5_SP*(XX(J+1)*YY(J)-XX(J)*YY(J+1))
       ENDDO
       ART1(I)=ABS(ART1(I))
     ELSE
       DO J=1,NTVE(I)
         II=NBVE(I,J)
         J1=NBVT(I,J)
         J2=J1+1-INT((J1+1)/4)*3
         XX(2*J-1)=(VX(NV(II,J1))+VX(NV(II,J2)))*0.5_SP-VX(I)
         YY(2*J-1)=(VY(NV(II,J1))+VY(NV(II,J2)))*0.5_SP-VY(I)
         XX(2*J)=XC(II)-VX(I)
         YY(2*J)=YC(II)-VY(I)
       ENDDO
       J=NTVE(I)+1
       II=NBVE(I,J-1)
       J1=NBVT(I,NTVE(I))
       J2=J1+2-INT((J1+2)/4)*3

       XX(2*J-1)=(VX(NV(II,J1))+VX(NV(II,J2)))*0.5_SP-VX(I)
       YY(2*J-1)=(VY(NV(II,J1))+VY(NV(II,J2)))*0.5_SP-VY(I)

       XX(2*J)=VX(I)-VX(I)
       YY(2*J)=VY(I)-VY(I)

       XX(2*J+1)=XX(1)
       YY(2*J+1)=YY(1)

       DO J=1,2*NTVE(I)+2
        ART1(I)=ART1(I)+0.5_SP*(XX(J+1)*YY(J)-XX(J)*YY(J+1))
       ENDDO
       ART1(I)=ABS(ART1(I))
     ENDIF
   ENDDO

!---COMPUTE AREA OF CONTROL VOLUME ART2(I) = SUM(ALL TRIS SURROUNDING NODE I)--!

   DO I=1,MLOC
     ART2(I) = SUM(ART(NBVE(I,1:NTVE(I))))  !summation of surrounding triangular elements area 
                                            !
   ENDDO

   ART(0) = ART(1) 
   ART1(0) = ART1(1) 
!   IF(NTLOC > NLOC)ART(NLOC+1:NTLOC) = ART(NLOC)
   IF(MTLOC > MLOC)ART2(MLOC+1:MTLOC) = ART2(MLOC)
   IF(MTLOC > MLOC)ART1(MLOC+1:MTLOC) = ART1(MLOC)
   DEALLOCATE(XX,YY)

# if defined (MULTIPROCESSOR)
   if(par)call node_match(1,NBN,BN_MLT,BN_LOC,BNC,MTLOC,1,MYID,NPROCS,art1)
   if(par)call node_match(1,NBN,BN_MLT,BN_LOC,BNC,MTLOC,1,MYID,NPROCS,art2)
# endif

	 

!		DO I = 1,MLOC
!			write(*,1245) I,ART1(I)
!		enddo
!		
   ! !  endif
   
   !1245 format(I4,F16.8)
   
   RETURN
   

   
   write(*,*)'Computed CELL area'
   
   
   
   
   END SUBROUTINE CELL_AREA
!==============================================================================|
