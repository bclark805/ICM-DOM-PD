!Subroutines :
!	Subroutine NCD_READ_GRID()
!	Subroutine NCD_READ_SHAPE()
!	Subroutine NCD_READ_OPEN()
!	Subroutine NCD_READ()

SUBROUTINE NCD_READ_GRID(INFILE)
  !---------------------------------------------------------------------
  ! READ DIMENSIONS IN A NETCDF FILES
  !---------------------------------------------------------------------
  USE MOD_NCD, ONLY : NC_FID, GETDIM, GETSVAR, GETDVAR
  USE NETCDF
  USE MOD_LIMS, ONLY : MLOC, NLOC, NTLOC, KBM2, KB, KBM1, MTLOC
  USE MOD_SIZES, ONLY: NGL, MGL
  USE MOD_PREC, ONLY : SP
  	!Wen Long took CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY : 		&
			!SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			MSR  !      	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			!,PAR        !	&           !!TRUE IF MULTIPROCESSOR RUN
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

  	
  IMPLICIT NONE
  !----------------------------------------------------------------------------!
  CHARACTER(LEN=1024), INTENT(IN) :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: N_ELEMS,N_NODES,N_SIG_M1,N_SIG
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP
  !----------------------------------------------------------------------------!

  !--Open NetCDF DATA FILE
  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
  if(MSR)WRITE(*,*)'opening netcdf file, INFILE, NC_FID, NGL '//TRIM(INFILE),NC_FID,NGL

  IF(IERR /=NF90_NOERR)THEN
     WRITE(*,*)'ERROR READING ',TRIM(INFILE)
     WRITE(*,*)TRIM(NF90_STRERROR(IERR))
     STOP
  ENDIF

  !--Get Model Dimensions
  N_ELEMS   = GETDIM(NC_FID,LEN_TRIM('nele'),'nele')
  N_NODES   = GETDIM(NC_FID,LEN_TRIM('node'),'node')
  N_SIG_M1  = GETDIM(NC_FID,LEN_TRIM('siglay'),'siglay')

  N_SIG     = GETDIM(NC_FID,LEN_TRIM('siglev'),'siglev')
  IF(N_ELEMS /= NGL)THEN
    IF(MSR)THEN 
      WRITE(*,*) 'NGL is inconsistent with the # of ELEMENTS in NetCDF file'
      WRITE(*,*) 'NGL should be equal to',N_ELEMS,NGL
     STOP
    ENDIF
  ENDIF    
  IF(N_NODES /= MGL)THEN
    IF(MSR)THEN 
      WRITE(*,*) 'MGL IS inconsistent with the # of NODES in NetCDF file'
      WRITE(*,*) 'MGL should be equal to',N_NODES
      STOP
    ENDIF
  ENDIF    
  IF(N_SIG /= KB)THEN
    IF(MSR)THEN 
      WRITE(*,*) 'KB IS inconsistent with the # of SIGMA LAYERS in NetCDF file'
      WRITE(*,*) 'KB should be equal to',N_SIG
     STOP
    ENDIF
  ENDIF    
  
  MLOC=N_NODES
  NLOC=N_ELEMS
  MTLOC = MLOC
  NTLOC = NLOC
  
  KB=N_SIG    
  KBM1=N_SIG_M1
  KBM2=KB-2

  !--close file
  IERR = NF90_CLOSE(NC_FID)

  RETURN
END SUBROUTINE NCD_READ_GRID

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE NCD_READ_SHAPE(INFILE)
  !---------------------------------------------------------------------
  ! READ BATHYMETRY, SIGMA LEVELS AND GRID COEFFICIENTS IN A NETCDF FILES
  !---------------------------------------------------------------------
  USE MOD_NCD, ONLY : NC_FID, GETSVAR
  USE NETCDF
  USE MOD_LIMS, ONLY: KBM2,MTLOC,MLOC
  USE MOD_HYDROVARS, ONLY: &
   		!GRAV	&		!
		!,PI	&		!
		!,PI2	&		!
		!,ZERO	&		!
		!,ONE_THIRD	&	!
		NVG	&		!
		,XG	&		!GLOBAL X-COORD AT NODE 
		,YG	&		!GLOBAL X-COORD AT NODE 
		,HG	&		!GLOBAL DEPTH AT NODE 

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
, SIGLEV  &
, SIGLAY &
		,Z	&			!SIGMA COORDINATE VALUE 
		,ZZ	&		!INTRA LEVEL SIGMA VALUE
		,DZ	&		!DELTA-SIGMA VALUE
		,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
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
		,num_hyd_ints  !	&!number of records in each hydrodynamics netcdf file
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

	!Wen Long took CONTROL out of MOD_HYDROVARS and put the used variables here
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
			
  USE MOD_LIMS, ONLY :  KB, KBM1
  USE MOD_SIZES, ONLY: NGL, MGL
  USE MOD_PREC, ONLY : SP
 
  IMPLICIT NONE
  !----------------------------------------------------------------------------!
  CHARACTER(LEN=1024), INTENT(IN) :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: I,K
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP

  !--OPEN NETCDF DATA FILE
  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
  if(MSR)WRITE(*,*)'opening netcdf file, INFILE, NC_FID'//TRIM(INFILE),NC_FID
  IF(IERR /= NF90_NOERR)THEN
    WRITE(*,*)'ERROR READING ',TRIM(INFILE)
    WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    STOP
  ENDIF

  !--Get Node Coordinates
  !WLong moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F 
  !ALLOCATE(XG(0:MGL),YG(0:MGL)) ; XG = 0.0_SP ; YG = 0.0_SP
  
  ALLOCATE(TEMP(MGL,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('x'),'x',MGL,1,TEMP)
  XG(1:MGL) = TEMP(1:MGL,1)
  DEALLOCATE(TEMP)

  ALLOCATE(TEMP(MGL,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('y'),'y',MGL,1,TEMP)
  YG(1:MGL) = TEMP(1:MGL,1)
  DEALLOCATE(TEMP)

  !--Get Node Numbering
  !WLong moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F 
  !ALLOCATE(NVG(0:NGL,4)); NVG = 0

  ALLOCATE(TEMP(NGL,3))
  CALL GETSVAR(NC_FID,LEN_TRIM('nv'),'nv',NGL,3,TEMP)
  NVG(1:NGL,1:3) = TEMP(1:NGL,1:3)
  DEALLOCATE(TEMP)
  NVG(:,4) = NVG(:,1)

  !--Get Bathymetry
  !!WLong moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F 
  !ALLOCATE(HG(0:MGL))  ; HG = 0.0_SP

  ALLOCATE(TEMP(MGL,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('h'),'h',MGL,1,TEMP)
  HG(1:MGL) = TEMP(1:MGL,1)
  DEALLOCATE(TEMP)
!   DO I=1,MGL
!	IF(HG(I)< 1.0)THEN
!		HG(I)=1.
 !       END IF
!   END DO
 !WLong moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F 
 ! ALLOCATE(Z(KB))  ; Z   = 0.0_SP    !!SIGMA COORDINATE VALUE 
 ! ALLOCATE(ZZ(KB)) ; ZZ  = 0.0_SP    !!INTRA LEVEL SIGMA VALUE
 ! ALLOCATE(DZ(KB)) ; DZ  = 0.0_SP    !!DELTA-SIGMA VALUE
 ! ALLOCATE(DZZ(KB)); DZZ = 0.0_SP    !!DELTA OF INTRA LEVEL SIGMA 

  !--Get Sigma levels
  ALLOCATE(TEMP(KB,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('siglev'),'siglev',1,KB,TEMP)
  Z(1:KB) = TEMP(1:KB,1)

  DEALLOCATE(TEMP)
    !  do K=1,KB
!tykim00
!         Z(K)=0.0-(K-1)*(1.0/(KB*1.0-1.0))

!
!WLong: this vertical sigma coordinate formulation is hardired
!       we should have this information in the netcdf file itself
!       or we should get this **1.5 as a parameter from input control
!
        !Z(K) = -((K-1)/FLOAT(KB-1))**1.5  =
      !enddo
                 
  !--Compute derivative and intra-sigma levels
  DO K=1,KBM1
    ZZ(K)=0.5_SP*(Z(K)+Z(K+1))
    DZ(K)=Z(K)-Z(K+1)
  ENDDO
  ZZ(KB)=2.0_SP*ZZ(KBM1)-ZZ(KBM2)

  DO K=1,KBM2
    DZZ(K)=ZZ(K)-ZZ(K+1)
  ENDDO
  DZZ(KB-1)=0.0
  DZ(KB)=0.0
  	   DO I =1,MTLOC
	    DO K = 1,KB
	      SIGLEV(I,K)=Z(K)
		END DO
		DO K=1,KBM1
		  siglay(I,K)=ZZ(K)
		ENDDO
     END DO
	!    write(*,*) 'Z = ',Z
  ! write(*,*)'Siglev = ', SIGLEV(100,:)
  ! write(*,*)'Siglay = ',SIGLAY(1090,:)

  !--Close file
  IERR = NF90_CLOSE(NC_FID)

  RETURN
END SUBROUTINE NCD_READ_SHAPE

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
! 15 SEPT 2009, KURT GLAESEMANN, REMOVE ALL THE CODE FROM THIS ROUTINE AND CALL DUPLICATE CODE IN NCD_READ
SUBROUTINE NCD_READ_OPEN(INFILE,UL,VL,WTSL,UARD_OBCNL,XFLUX_OBCL,DTFAL,KHL,ELL,TL,SL,HO)
  !---------------------------------------------------------------------
  ! READ DATA FROM DAILY NETCDF FILES
  !---------------------------------------------------------------------

  USE MOD_NCD, ONLY : NF90_OPEN,NF90_NOWRITE,NF90_STRERROR,NF90_NOERR,NC_FID
  USE MOD_SIZES, ONLY : NOBTY
  USE MOD_LIMS, ONLY: MTLOC, NTLOC, KBM1, KB
  USE MOD_PREC, ONLY: SP
  
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
  IMPLICIT NONE
  !----------------------------------------------------------------------------!

  REAL(SP), DIMENSION(0:NTLOC,KB),INTENT(OUT)   :: UL,VL     !,WWL
  REAL(SP), DIMENSION(0:MTLOC,KB),INTENT(OUT)   :: KHL,WTSL
  REAL(SP), DIMENSION(0:MTLOC,KBM1),INTENT(OUT) :: TL,SL
  REAL(SP), DIMENSION(0:MTLOC),INTENT(OUT)      :: ELL,DTFAL
  REAL(SP), DIMENSION(0:NOBTY+1),INTENT(OUT)    :: UARD_OBCNL
  REAL(SP), DIMENSION(0:NOBTY,KBM1),INTENT(OUT) :: XFLUX_OBCL
  INTEGER, INTENT(IN)                       :: HO
  CHARACTER(LEN=1024), INTENT(IN)           :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
!  CHARACTER(LEN=300) TASK

  !--Open NetCDF Datafile 
  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
! KURT GLAESEMANN added MSR to print out
  

  if(MSR)WRITE(*,*)'opening netcdf file, INFILE, NC_FID'//TRIM(INFILE),NC_FID

  IF(IERR /=NF90_NOERR)THEN
     WRITE(*,*)'ERROR READING ',TRIM(INFILE)
     WRITE(*,*)TRIM(NF90_STRERROR(IERR))
     STOP
  ENDIF

  CALL NCD_READ(INFILE,UL,VL,WTSL,UARD_OBCNL,XFLUX_OBCL,DTFAL,KHL,ELL,TL,SL,HO)

  RETURN
END SUBROUTINE NCD_READ_OPEN


!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE NCD_READ(INFILE,UL,VL,WTSL,UARD_OBCNL,XFLUX_OBCL,DTFAL,KHL,ELL,TL,SL,HO)
  !---------------------------------------------------------------------
  ! READ DATA FROM DAILY NETCDF FILES
  !---------------------------------------------------------------------
  USE MOD_PREC, ONLY: SP
  USE MOD_NCD, ONLY : NC_FID, GETDVAR
  USE NETCDF
  
  USE MOD_WQM, ONLY: 		&!
			DTFAL_GL,       &!
            ELL_GL,         &!
            KHL_GL,         &!
            SL_GL,          &!
            TL_GL,          &!
            UARD_OBCN_GL,   &!
            UL_GL,          &!
            VL_GL,          &!
            WTSL_GL,        &!
            XFLUX_OBC_GL

  USE MOD_LIMS, ONLY: MTLOC, NTLOC, MLOC, NLOC, KB,KBM1
  USE MOD_SIZES, ONLY: MGL,NGL, NOBTY
  
  
  
  USE MOD_HYDROVARS, ONLY : NUM_HYD_INTS, SIGLEV, SIGLAY
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
		
  USE MOD_BCMAP, ONLY: 	&
		IOBCN,			&
		I_OBC_GL		

   
# if defined (MULTIPROCESSOR)
		USE MOD_PAR, ONLY:   &	  
			!EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			!ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			 NLID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			!ELID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			!NLID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,	&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			 EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			 NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
			 NHE,        &   !			!NUMBER OF HALO ELEMENTS
			 HE_LST,		&	!(:)       	!!GLOBAL IDENTITIES OF HALO ELEMENTS
			!HE_OWN,		&	!(:)       	!!OWNER OF HALO ELEMENTS
			 NBN ,       &   !          	!!NUMBER OF BOUNDARY NODES
			!MX_MLT,     &   !          	!!MAX MULTIPLICITY OF BOUNDARY NODES
			!BN_LST, 	&	!	   		!!GLOBAL IDENTITY OF BOUNDARY NODES
			 BN_LOC, 	&	!			!!LOCAL IDENTITY OF BOUNDARY NODES
			 BN_MLT, 	&	!			!!MULTIPLICITY OF BOUNDARY NODES
			!BN_NEY, 	&	!			!!NODE OWNER LIST
			!NDE_ID, 	&	!			!! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
			 NHN,  		&	!			!!NUMBER OF HALO NODES
			 HN_LST, 	&	!			!!LIST OF HALO NODES 
			!HN_OWN, 	&	!			!!PRIMARY OWNER OF HALO NODES
			!EC,			&	!
			!NC,			&	!
			 BNC!, 		&	!
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

# endif  
  IMPLICIT NONE

!----------------------------------------------------------------------------!
!
!  REAL(CDF_PREC), DIMENSION(NGL,KBM1)   :: UL_GL,VL_GL                 !u and v
!  REAL(CDF_PREC), DIMENSION(MGL,KB)     :: WTSL_GL,KHL_GL              !wts, kh
!  REAL(CDF_PREC), DIMENSION(MGL,KBM1)   :: SL_GL,TL_GL                 !salinity,temp
!  REAL(CDF_PREC), DIMENSION(MGL)        :: ELL_GL,DTFAL_GL             !zeta,dtfa
!  REAL(CDF_PREC), DIMENSION(NOBTY+1)    :: UARD_OBCN_GL                !uard_obcn
!  REAL(CDF_PREC), DIMENSION(NOBTY,KBM1) :: XFLUX_OBC_GL                !xflux_obc
!

  REAL(SP), DIMENSION(0:NTLOC,KB),INTENT(OUT)   :: UL,VL     !,WWL
  REAL(SP), DIMENSION(0:MTLOC,KB),INTENT(OUT)   :: KHL,WTSL
  REAL(SP), DIMENSION(0:MTLOC,KBM1),INTENT(OUT) :: TL,SL
  REAL(SP), DIMENSION(0:MTLOC),INTENT(OUT)      :: ELL,DTFAL
  REAL(SP), DIMENSION(0:NOBTY+1),INTENT(OUT)    :: UARD_OBCNL
  REAL(SP), DIMENSION(0:NOBTY,KBM1),INTENT(OUT) :: XFLUX_OBCL
  INTEGER, INTENT(IN)                       :: HO
  CHARACTER(LEN=1024), INTENT(IN)           :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: HT
  INTEGER            :: I,K,J

!--Adjustement to read in Netcdf file
!  HT=HO+1
   HT=HO 

  !---------------------------------------------------------------------
  ! Read Data from file INFILE at time level ht
  !---------------------------------------------------------------------

! KURT GLAESEMANN SEPT 24 2009 - REORDER DATA READS TO MATCH ORDER IN FILE - FASTER
  !--U velocity
  CALL GETDVAR(NC_FID,LEN_TRIM('u'),'u',NGL,KBM1,UL_GL,HT)

  !--V velocity
  CALL GETDVAR(NC_FID,LEN_TRIM('v'),'v',NGL,KBM1,VL_GL,HT)

 !--W Velocity 
 ! CALL GETDVAR(NC_FID,LEN_TRIM('ww'),'ww',NGL,KBM1,WW_GL,HT)
!--Wen Long debugging----
!  IF(MSR)THEN
!     WRITE(*,'(A32,1x,F20.10,1x,F20.10,1x,F20.10)')'Wen Long (ELTMS,U,V) Surface == ',ELTMS,UL_GL(565,1),VL_GL(565,1)
!  ENDIF
!------------------------

  !--WTS velocity
  CALL GETDVAR(NC_FID,LEN_TRIM('omega'),'omega',MGL,KB,WTSL_GL,HT)

  !--UARD_OBCN
  CALL GETDVAR(NC_FID,LEN_TRIM('uard_obcn'),'uard_obcn',NOBTY,1,UARD_OBCN_GL,HT)

  !--XFLUX_OBCN
  CALL GETDVAR(NC_FID,LEN_TRIM('xflux_obc'),'xflux_obc',NOBTY,KBM1,XFLUX_OBC_GL,HT)

  !--free surface elevation
  CALL GETDVAR(NC_FID,LEN_TRIM('dtfa'),'dtfa',MGL,1,DTFAL_GL,HT)

! KURT GLAESEMAN 23 SEPT 2009 - fix dimensions of KH (KB, not KBM1)- here and above and below
  !--KH
  CALL GETDVAR(NC_FID,LEN_TRIM('kh'),'kh',MGL,KB,KHL_GL,HT)

  !--free surface elevation
  CALL GETDVAR(NC_FID,LEN_TRIM('zeta'),'zeta',MGL,1,ELL_GL,HT)

  !--salinity
  CALL GETDVAR(NC_FID,LEN_TRIM('salinity'),'salinity',MGL,KBM1,SL_GL,HT)

  !--temperature
  CALL GETDVAR(NC_FID,LEN_TRIM('temp'),'temp',MGL,KBM1,TL_GL,HT)
  
  ! ! B Clark add sanity check for bad temperatures to correct to reasonable value;
  !   DO I = 1,MGL
	!   DO K = 1,KBM1
	!     IF (TL_GL(I,K) < -5.0) THEN
	!         TL_GL(I,K) = -5.0
	! 		write(*,*)'corrected temp at',MGL
	! 		write(*,*)'back to a value of',TL_GL(I,K)
  !       ENDIF
  !     END DO
  !   END DO
! 15 SEPT KURT GLAESEMANN - MOVED ALL COPIES TOGETHER
  IF(SERIAL) THEN
     ELL(1:MGL) = ELL_GL(1:MGL)
     SL(1:MGL,1:KBM1) = SL_GL(1:MGL,1:KBM1)
     TL(1:MGL,1:KBM1) = TL_GL(1:MGL,1:KBM1)
     UL(1:NGL,1:KBM1) = UL_GL(1:NGL,1:KBM1)
     VL(1:NGL,1:KBM1) = VL_GL(1:NGL,1:KBM1)
     WTSL(1:MGL,1:KB) = WTSL_GL(1:MGL,1:KB)
     UARD_OBCNL(1:NOBTY) = UARD_OBCN_GL(1:NOBTY)
     XFLUX_OBCL(1:NOBTY,1:KBM1) = XFLUX_OBC_GL(1:NOBTY,1:KBM1)
     DTFAL(1:MGL) = DTFAL_GL(1:MGL)
     KHL(1:MGL,1:KB) = KHL_GL(1:MGL,1:KB)
  ENDIF

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    DO I=1,MLOC
	  !WRITE(*,*)'ELL_GL Debug here here 05 '
      ELL(I) = ELL_GL(NGID(I))
	  !WRITE(*,*)'ELL_GL Debug here here 5 '
    ENDDO

    DO I=1,NHN
	  !WRITE(*,*)'ELL_GL Debug here here 06 '
      ELL(I+MLOC) = ELL_GL(HN_LST(I))
	  !WRITE(*,*)'ELL_GL Debug here here 6 '
    ENDDO
!
    DO K=1,KBM1
      DO I=1,MLOC
        SL(I,K) = SL_GL(NGID(I),K)
      ENDDO

      DO I=1,NHN
        SL(I+MLOC,K) = SL_GL(HN_LST(I),K)
      ENDDO
    ENDDO

!
    DO K=1,KBM1
      DO I=1,MLOC
        TL(I,K) = TL_GL(NGID(I),K)
      ENDDO

      DO I=1,NHN
        TL(I+MLOC,K) = TL_GL(HN_LST(I),K)
      ENDDO
    ENDDO
!
    DO K=1,KBM1
      DO I=1,NLOC
        UL(I,K) = UL_GL(EGID(I),K)
      ENDDO

      DO I=1,NHE
        UL(I+NLOC,K) = UL_GL(HE_LST(I),K)
      ENDDO
    ENDDO
!
    DO K=1,KBM1
      DO I=1,NLOC
        VL(I,K) = VL_GL(EGID(I),K)
      ENDDO

      DO I=1,NHE
        VL(I+NLOC,K) = VL_GL(HE_LST(I),K)
      ENDDO
    ENDDO
!
    DO K=1,KB
      DO I=1,MLOC
        WTSL(I,K) = WTSL_GL(NGID(I),K)
      ENDDO

      DO I=1,NHN
        WTSL(I+MLOC,K) = WTSL_GL(HN_LST(I),K)
      ENDDO
    ENDDO

       UARD_OBCNL = -99999
       K = 0
       DO I = 1, NOBTY
         J = NLID(I_OBC_GL(I))
         IF ( J .NE. 0) THEN
            K = K + 1
            UARD_OBCNL(K) = UARD_OBCN_GL(I)
         ENDIF
       ENDDO
       IF (K .NE. IOBCN) THEN
         write(*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
         call pstop
       endif
       XFLUX_OBCL = -99999
       K = 0
       DO I = 1, NOBTY
         J = NLID(I_OBC_GL(I))
         IF ( J .NE. 0) THEN
            K = K + 1
            XFLUX_OBCL(K,:) = XFLUX_OBC_GL(I,:)
         ENDIF
       ENDDO
       IF (K .NE. IOBCN) THEN
         write(*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
         call pstop
       endif
!
    DO I=1,MLOC
      DTFAL(I) = DTFAL_GL(NGID(I))
    ENDDO

    DO I=1,NHN
      DTFAL(I+MLOC) = DTFAL_GL(HN_LST(I))
    ENDDO
!
    DO K=1,KB
      DO I=1,MLOC
        KHL(I,K) = KHL_GL(NGID(I),K)
      ENDDO

      DO I=1,NHN
        KHL(I+MLOC,K) = KHL_GL(HN_LST(I),K)
      ENDDO
    ENDDO

  ENDIF
# endif

  IF(HT.eq.num_hyd_ints) THEN
   IERR = NF90_CLOSE(NC_FID)
  ENDIF 

  RETURN
END SUBROUTINE NCD_READ
