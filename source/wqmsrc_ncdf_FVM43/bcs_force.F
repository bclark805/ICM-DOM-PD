!subroutine BCS_FORCE()

!==============================================================================|
!   Set Up the Following Boundary Conditions:                                  |
!     Freshwater River Water Quality Input                                     |
!==============================================================================|

   SUBROUTINE BCS_FORCE           

!------------------------------------------------------------------------------|
    USE MOD_SIZES, ONLY :        &	!
			 NCP !,            &  !
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
            !MGL,            &  !
			!NGL,            &  !
			!OBCGL,          &  !
			!NOBTY  
  
   
	USE MOD_LIMS, ONLY: NUMQBC_GL,NUMQBC, NUMPNT_GL, NUMPNT, KBM1 &
				, MLOC
	USE MOD_PREC, ONLY: SP				
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
    USE MOD_CONTROL, ONLY : 		&
			SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			,MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			,PAR        	&           !!TRUE IF MULTIPROCESSOR RUN
			,CASENAME  	&   		!!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
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
			,RIV_FILENUMBER	&			!!
            ,INFLOW_TYPE   	&			!!SPECIFIED RIVER INFLOW TYPE (edge/node) 
            ,POINT_ST_TYPE 	&			!!(calculated/specified)
            ,PNT_SOURCE !   	&			!!point_source
            !,DAY				&
			!,in_jday		
			
   USE MOD_BCS, ONLY :      &
         INODEQ		    &  !!LOCAL FRESH WATER INFLOW NODES
		,ICELLQ     	&  !!LOCAL FRESH WATER INFLOW ELEMENTS 
		,RIV_GL2LOC  	&  !!GLOBAL TO LOCAL MAP OF FW NODES
		,INOPNT      	&  !!LOCAL NON-POINT SOURCE NODES
		,QBC_TM         &  !!TIME MAP FOR RIVER DISCHARGE DATA 
		,PNT_TM         &  !!TIME MAP FOR NON POINT SOURCE DATA

		,WVQDIST   		&  !!DISCHARGE VERTICAL DISTRIBUTION for point source
		,DWDIS   		&  !!WATER QUALITY DISCHARGE DATA
		,PDWQDIS 		&  !!WATER QUALITY DATA at non-point source node
		,PDQDIS     	&  !!DISCHARGE at non-point source node
		,N_ICELLQ		   !node number of the edge for element edge source

	USE MOD_UTILS, ONLY: &
		PERROR,   		&!
		!GET_TIMESAMP,	&!
		FOPEN!,			&!
		!CRAY_SYSTEM_CALL
		
   
   USE MOD_WQM, ONLY: JDAY, IJDAY
# if defined (MULTIPROCESSOR)
	USE MOD_PAR, ONLY:   &	  
			!EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			NLID!,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			!ELID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			!NLID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,	&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			!EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			!NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
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
			!NHN,  		&	!			!!NUMBER OF HALO NODES
			!HN_LST!, 	&	!			!!LIST OF HALO NODES 
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
			!FILELOCK,   &		!
			!MYFILELOCK, &		!
			!FILELOCKED , &		!
			!NODE_MATCH, &		!
			!NODE_MATCH_ALL, &	!
			!EXCHANGE, &		!
			!EXCHANGE_ALL, &	!
			!GATHER, &			!
			!IGATHER, &			!
			!INIT_MPI_ENV,&!			
			!SORT,	&			!
			!PREAD,	&			!
			!PWRITE,&			!
			!PPRINT, &			!
			!LOCK, &			!
			!UNLOCK, &			
			!GETLOC 

		USE MOD_PREC, ONLY: MPI_F
# endif

   IMPLICIT NONE
# if defined (MULTIPROCESSOR)
   include "mpif.h"
# endif
   CHARACTER(LEN=80) :: COMT
   REAL(SP) :: TTIME
   REAL(SP) :: FTEMP1,FTEMP2,FTEMP3,RBUF1,RBUF2,RBUF3
   REAL(SP), ALLOCATABLE :: RTEMP(:,:),RTEMP1(:,:,:),RTEMP2(:,:)
   INTEGER,  ALLOCATABLE :: RIV_2_LOC_PNT(:)
   INTEGER,  ALLOCATABLE :: TEMP(:),TEMP2(:),TEMP3(:),TEMP4(:),NODE_SBC(:)
!   INTEGER,  ALLOCATABLE :: TEMPD(:,:),TEMP2D(:,:),TEMP3D(:,:) !WLong not used at all
   INTEGER   I,J,K,NQTIME,NBFWTIME,ISBCN1,INMAX,IOS,NCNT,IGL,IERR,JN
   INTEGER   PNQTIME
   CHARACTER(LEN=13) :: TSTRING
   CHARACTER(LEN=100) :: ISTR,ISTPT
   INTEGER :: INRIV,NJ,INPT
   
   INRIV = 61
   INPT = 63

!------------------------------------------------------------------------------|


!----------------------------REPORT--------------------------------------------!
   IF(MSR)WRITE(*,*  )'!'
   IF(MSR)WRITE(*,*)'!           SETTING UP PRESCRIBED BOUNDARY CONDITIONS   '
   IF(MSR)WRITE(*,*  )'!'

!==============================================================================|
!   Input River/Dam/Intake/Outfall Boundary Values                             |
!==============================================================================|

!
!-------Check Selected Combination for Validity--------------------------------!
!
   IJDAY=AINT(JDAY)
   WRITE(RIV_FILENUMBER,'(I4.4)') IJDAY
   ISTR = "./inputs/"//trim(casename)//trim(RIV_FILENUMBER)
    write(*,*) ISTR
   CALL FOPEN(INRIV, TRIM(ISTR)//'_riv.dat',"cfr")
!
   REWIND(INRIV)
   READ(INRIV,'(A4,2X,A10)') INFLOW_TYPE,POINT_ST_TYPE

   IF(INFLOW_TYPE /= 'edge' .AND. INFLOW_TYPE /= 'node') THEN
     CALL PERROR(6,"INFLOW TYPE NOT CORRECT","SHOULD BE edge or node")
   ENDIF

   IF(POINT_ST_TYPE /= 'calculated' .AND. POINT_ST_TYPE /= 'specified') THEN
     CALL PERROR(6,"POINT_ST TYPE NOT CORRECT","SHOULD BE calculated or specified")
   ENDIF
   write(*,*)'opened the river fle' ! B Clark Debug
!
!--Read in Number of Discharge Nodes/Edges-------------------------------------!
!
   IF(MSR)THEN
     READ(INRIV,*) NUMQBC_GL
   ENDIF

#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_BCAST(NUMQBC_GL,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#  endif

   NUMQBC = 0
   IF(NUMQBC_GL > 0)THEN

     ALLOCATE(TEMP(NUMQBC_GL),TEMP2(NUMQBC_GL),TEMP3(NUMQBC_GL))
     IF(MSR)THEN
       DO I=1,NUMQBC_GL
         READ(INRIV,*) TEMP(I)

       ENDDO
     ENDIF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(TEMP,NUMQBC_GL,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

!
!--Determine Global--> Local Mapping of Freshwater Discharge Nodes
!
     IF(SERIAL)THEN
       NUMQBC = NUMQBC_GL
       IF(INFLOW_TYPE == 'node') THEN
         ALLOCATE(INODEQ(NUMQBC))
         INODEQ = TEMP
       ELSEIF(INFLOW_TYPE == 'edge') THEN
         ALLOCATE(ICELLQ(NUMQBC))
         ICELLQ = TEMP(1:NUMQBC)
		 ALLOCATE(N_ICELLQ(NUMQBC,2))  !WLong added this to record the node numbes of the edge
									   !which should be calculated using set_bndry.F
		 N_ICELLQ=0
       ENDIF
     ENDIF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       NCNT = 0
       IF(INFLOW_TYPE == 'node') THEN
         DO I=1,NUMQBC_GL
           IF(NLID(TEMP(I)) /= 0)THEN
             NCNT = NCNT + 1
             TEMP2(NCNT) = NLID(TEMP(I))
             TEMP3(NCNT) = I
           ENDIF
         ENDDO
         NUMQBC = NCNT
         ALLOCATE(INODEQ(NUMQBC),RIV_GL2LOC(NUMQBC))
         INODEQ = TEMP2(1:NCNT)
         RIV_GL2LOC = TEMP3(1:NCNT)
       ELSE IF(INFLOW_TYPE == 'edge') THEN
         DO I=1,NUMQBC_GL
           IF(ELID(TEMP(I)) /= 0)THEN
             NCNT = NCNT + 1
             TEMP2(NCNT) = ELID(TEMP(I))
             TEMP3(NCNT) = I
           ENDIF
         ENDDO
         NUMQBC = NCNT
         ALLOCATE(ICELLQ(NUMQBC),RIV_GL2LOC(NUMQBC))
		 ALLOCATE(N_ICELLQ(NUMQBC,2))  !WLong added this to record the node numbes of the edge
									   !which should be calculated using set_bndry.F
		 N_ICELLQ=0
		 
         ICELLQ = TEMP2(1:NCNT)
         RIV_GL2LOC = TEMP3(1:NCNT)
       ENDIF
     ENDIF
#   endif

     DEALLOCATE(TEMP,TEMP2,TEMP3)

!
!----Read in Time Dependent DataSets (DQDIS,DSDIS,DTDIS)------------------------!
!
     IF(MSR)READ(INRIV,*) NQTIME

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(NQTIME,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

     QBC_TM%NTIMES = NQTIME
     QBC_TM%LABEL  = "Freshwater Discharge"
     ALLOCATE(QBC_TM%TIMES(NQTIME))
     ALLOCATE(RTEMP1(NUMQBC_GL,NCP,NQTIME))
!    KURT GLAESEMANN added MSR to print out
!     IF(MSR)WRITE(*,*)'NQTIME', NQTIME
     IF(MSR)THEN
       DO I = 1, NQTIME
         READ(INRIV,*) TTIME
         QBC_TM%TIMES(I) = TTIME
!RGL this is where wq data read in for rivers
!need to create software to generate input file
         DO NJ=1,NCP
           READ(INRIV,*) (RTEMP1(J,NJ,I),J = 1,NUMQBC_GL)
         ENDDO

       ENDDO
     ENDIF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(QBC_TM%TIMES,NQTIME,MPI_F,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_BCAST(RTEMP1(:,:,:),NUMQBC_GL*NQTIME*NCP,MPI_F,0,MPI_COMM_WORLD,IERR)
#    endif

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!

     IF(NUMQBC > 0)THEN
       ALLOCATE(DWDIS(NUMQBC,NCP,NQTIME))

       IF(SERIAL)THEN
         DWDIS(1:NUMQBC_GL,:,:) = RTEMP1(1:NUMQBC_GL,:,:)
       ENDIF

#     if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         DO I=1,NQTIME
           DO NJ = 1,NCP
             DWDIS(1:NUMQBC,NJ,I) = RTEMP1(RIV_GL2LOC(1:NUMQBC),NJ,I)
           ENDDO
         ENDDO
       ENDIF
#     endif
	   

     ENDIF

     DEALLOCATE(RTEMP1)

     CLOSE(INRIV)

   ENDIF !! NUMQBC_GL > 0
!TYKIM
!==============================================================================
! NON-POINT source Values                             |
!==============================================================================
!
!-------Check Selected Combination for Validity--------------------------------
!
   ISTPT = "./inputs/"//trim(casename)
   CALL FOPEN(INPT, TRIM(ISTPT)//'_pnt_wq.dat',"cfr")!
   REWIND(INPT)

   READ(INPT,'(A5)') PNT_SOURCE

   IF(PNT_SOURCE /= 'point') THEN
     CALL PERROR(6,"POINT_SOURCE INPUT FILE IS  NOT CORRECT")
   ENDIF
!
!--Read in Number of Discharge Nodes/Edges-------------------------------------

   READ(INPT,*) NUMPNT_GL

   NUMPNT = 0
   IF(NUMPNT_GL > 0)THEN

     ALLOCATE(TEMP(NUMPNT_GL),TEMP2(NUMPNT_GL),TEMP3(NUMPNT_GL))
     DO I=1,NUMPNT_GL
       READ(INPT,*) TEMP(I)
     ENDDO

!--Determine Global--> Local Mapping of Freshwater Discharge Nodes
!
     IF(SERIAL)THEN
       NUMPNT = NUMPNT_GL
       ALLOCATE(INOPNT(NUMPNT))
       INOPNT(:) = TEMP(:)
     ENDIF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       NCNT = 0
       DO I=1,NUMPNT_GL
         IF(NLID(TEMP(I)) /= 0)THEN
           NCNT = NCNT + 1
           TEMP2(NCNT) = NLID(TEMP(I))
           TEMP3(NCNT) = I
         ENDIF
       ENDDO
       NUMPNT = NCNT
       ALLOCATE(RIV_2_LOC_PNT(NUMPNT),INOPNT(NUMPNT))
       INOPNT(1:NCNT) = TEMP2(1:NCNT) ! The point locally
       RIV_2_LOC_PNT(1:NCNT) = TEMP3(1:NCNT) ! The point globally
     ENDIF
#   endif

     DEALLOCATE(TEMP,TEMP2,TEMP3)

!
!----Read in Time Dependent Data and Flow Distribution Rate!
!
     ALLOCATE(RTEMP(NUMPNT_GL,KBM1))
     DO I = 1, NUMPNT_GL
       READ(INPT,*) J,(RTEMP(I,K),K = 1,KBM1)
     ENDDO

     ALLOCATE(WVQDIST(NUMPNT,KBM1))
     IF(SERIAL) WVQDIST(1:NUMPNT,:) = RTEMP(1:NUMPNT_GL,:)
#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       WVQDIST(1:NUMPNT,:) = RTEMP(RIV_2_LOC_PNT(1:NUMPNT),:)
     ENDIF
#   endif

!    READ DATA
     READ(INPT,*) PNQTIME

     PNT_TM%NTIMES = PNQTIME
     PNT_TM%LABEL  = "Point source"
     ALLOCATE(PNT_TM%TIMES(PNQTIME))
     ALLOCATE(RTEMP1(NUMPNT_GL,NCP,PNQTIME))
     ALLOCATE(RTEMP2(NUMPNT_GL,PNQTIME))

     DO I = 1, PNQTIME
       READ(INPT,*) TTIME
       PNT_TM%TIMES(I) = TTIME
       READ(INPT,*) (RTEMP2(J,I),J = 1,NUMPNT_GL)
       DO NJ=1,NCP
         READ(INPT,*) (RTEMP1(J,NJ,I),J = 1,NUMPNT_GL)
       ENDDO
     ENDDO
!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
       ALLOCATE(PDQDIS(NUMPNT,PNQTIME))
       ALLOCATE(PDWQDIS(NUMPNT,NCP,PNQTIME))
!
       IF(SERIAL)THEN
         PDQDIS(1:NUMPNT,:) = RTEMP2(1:NUMPNT_GL,:)
         PDWQDIS(1:NUMPNT,:,:) = RTEMP1(1:NUMPNT_GL,:,:)
       ENDIF

#     if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         PDQDIS(1:NUMPNT,:) = RTEMP2(RIV_2_LOC_PNT(1:NUMPNT),:)
         PDWQDIS(1:NUMPNT,:,:) = RTEMP1(RIV_2_LOC_PNT(1:NUMPNT),:,:)
         DEALLOCATE(RIV_2_LOC_PNT)
       ENDIF
#     endif


     DEALLOCATE(RTEMP,RTEMP1,RTEMP2)

     CLOSE(INPT)

   ENDIF !! NUMQBC_GL > 0

   RETURN
   END SUBROUTINE BCS_FORCE
!==============================================================================|
