!subroutine GENMAP()

!==============================================================================!
!  SET UP GLOBAL <--> LOCAL NODE AND ELEMENT NUMBERING MAPS                    !
!  SET UP INTERPROCESSOR COMMUNICATION                                         !
!==============================================================================!

   SUBROUTINE GENMAP 
# if defined (MULTIPROCESSOR)
!==============================================================================!
   USE MOD_PREC, ONLY: SP
   USE MOD_LIMS, ONLY: NLOC, MLOC, MYID, NPROCS,KB, NTLOC, MTLOC
   USE MOD_HYDROVARS, ONLY:  &
   		!GRAV	&		!
		!,PI	&		!
		!,PI2	&		!
		!,ZERO	&		!
		!,ONE_THIRD	&	!
		NVG	!&		!
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

   
   !Wen Long commented out the following CONTROL module for it is already included in MOD_HYDROVARS
    USE MOD_CONTROL, ONLY : 		&
			SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			,MSR  !      	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			!,PAR        	&           !!TRUE IF MULTIPROCESSOR RUN
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
			 
	USE MOD_PAR, ONLY:   &	  
			EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			NLID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			ELID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			NLID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,	&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
			NHE,        &   !			!NUMBER OF HALO ELEMENTS
			HE_LST,		&	!(:)       	!!GLOBAL IDENTITIES OF HALO ELEMENTS
			HE_OWN,		&	!(:)       	!!OWNER OF HALO ELEMENTS
			NBN ,       &   !          	!!NUMBER OF BOUNDARY NODES
			MX_MLT,     &   !          	!!MAX MULTIPLICITY OF BOUNDARY NODES
			BN_LST, 	&	!	   		!!GLOBAL IDENTITY OF BOUNDARY NODES
			BN_LOC, 	&	!			!!LOCAL IDENTITY OF BOUNDARY NODES
			BN_MLT, 	&	!			!!MULTIPLICITY OF BOUNDARY NODES
			BN_NEY, 	&	!			!!NODE OWNER LIST
			NDE_ID, 	&	!			!! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
			NHN,  		&	!			!!NUMBER OF HALO NODES
			HN_LST, 	&	!			!!LIST OF HALO NODES 
			HN_OWN, 	&	!			!!PRIMARY OWNER OF HALO NODES
			EC,			&	!
			NC,			&	!
			BNC, 		&	!
			EMAP,		&	!
			NMAP, 		&	!
			PNE,		&	!(:)        !!NUMBER OF ELEMENTS IN EACH PROC
			PNN,		&	!(:)        !!NUMBER OF NODES IN EACH PROC
			PNHE,		&	!(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
			PNBN,		&	!(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
			PMBM,		&	!(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
			PNHN,		&	!(:)       !!NUMBER OF HALO NODES IN EACH PROC
			!FILELOCK,   &  !
			!MYFILELOCK, &  !
			!FILELOCKED , &  !
			 NODE_MATCH, &  !
			! NODE_MATCH_ALL, &  !
			 EXCHANGE, & !
			!EXCHANGE_ALL, & !
			GATHER, & !
			!IGATHER, & !
			!INIT_MPI_ENV,&!			 
			 SORT!,	& !
			!PREAD,	& !
			!PWRITE,& !
			!PPRINT, & !
			!LOCK, & !
			!UNLOCK, & !
			!GETLOC 

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
# if defined (MULTIPROCESSOR)
   include "mpif.h"
# endif
   INTEGER I,J,K,N1 ,N2,N3,IERR,CHECK,NGLOB,NN,I1,I2,I3,NSZE,ISND,IRCV,OWNER
   INTEGER, ALLOCATABLE :: NP(:,:),NTEMP(:),NTEMP2(:),NTEMP3(:,:),NTEMP4(:,:)
   INTEGER, ALLOCATABLE :: NODEMARK(:)
   INTEGER :: NNODES_MINE
   INTEGER IGLID,IC1,IC2,TERCV
   INTEGER STAT(MPI_STATUS_SIZE)
   LOGICAL :: NEW,INELEM,ISMYN
   REAL(SP), ALLOCATABLE :: TESTA(:,:),TESTB(:,:),TESTC(:,:)
   REAL(SP)  DIFF,VAL,TIME_PAR,NEAVE,NNAVE
   REAL(SP) :: FCHECK
   INTRINSIC SYSTEM_CLOCK


   IF(MSR)WRITE(*,*  )'!'
   IF(MSR)WRITE(*,*)'!            SETTING UP MAPPING/MESSAGE PASSING '
   IF(MSR)WRITE(*,*  )'!'

!==============================================================================|
!   COUNT NUMBER OF LOCAL DOMAIN ELEMENTS                        NLOC             |
!   CREATE LOCAL--> GLOBAL ELEMENT MAP                           EGID(1:NLOC)     |
!   EGID(I) = GLOBAL ELEMENT NUMBER OF LOCAL ELEMENT I                         |
!==============================================================================|
   ALLOCATE(NTEMP(NGL)) ; NTEMP = 0
   NLOC = 0
   DO I=1,NGL
     IF(EL_PID(I) == MYID) THEN
       NLOC = NLOC + 1
       NTEMP(NLOC) = I
     ENDIF
   ENDDO
   ALLOCATE(EGID(0:NLOC)) ; EGID = 0
   EGID(1:NLOC) = NTEMP(1:NLOC)
   DEALLOCATE(NTEMP)
   IF(MSR)WRITE(*,*)  '!  ELEMENT MAP           :    COMPLETE'
!==============================================================================|
!   COUNT NUMBER OF LOCAL NODES                                  MLOC             |
!   CREATE LOCAL--> GLOBAL NODE MAP                              NGID(1:MLOC)     |
!   NGID(I) = GLOBAL NODE NUMBER OF LOCAL NODE I                               |
!==============================================================================|
   ALLOCATE(NODEMARK(MGL)) ; NODEMARK = 0
   ALLOCATE(NTEMP(MGL))
   DO I=1,NGL
     IF(EL_PID(I) == MYID)THEN
       NODEMARK(NVG(I,1:3)) = 1
     ENDIF
   ENDDO

   MLOC = 0
   DO I=1,MGL
     IF(NODEMARK(I) == 1)THEN
       MLOC = MLOC + 1
       NTEMP(MLOC) = I
     ENDIF
   ENDDO
!   DEALLOCATE(NODEMARK)


   ALLOCATE(NGID(MLOC))
   NGID(1:MLOC) = NTEMP(1:MLOC)
   DEALLOCATE(NTEMP)
   
   IF(MSR)WRITE(*,*)  '!  NODE MAP              :    COMPLETE'

!==============================================================================|
!   DETERMINE NUMBER OF HALO ELEMENTS                            NHE           |
!   DETERMINE GLOBAL MAPPING OF HALO ELEMENTS                    HE_LST(NHE)   | 
!   DETERMINE LOCAL MAPPING OF HALO ELEMENTS                     HE_LOC(NHE)   | 
!==============================================================================|
   
   ALLOCATE(NTEMP(NGL))
   ALLOCATE(NTEMP2(NGL))
   NHE = 0
   DO I=1,NGL
     NNODES_MINE = SUM(NODEMARK(NVG(I,1:3)))
!     IF(NNODES_MINE == 1 .OR. NNODES_MINE ==2)THEN
     IF(NNODES_MINE > 0 .AND. EL_PID(I) /= MYID)THEN
       NHE         = NHE + 1
       NTEMP(NHE)  = I
       NTEMP2(NHE) = EL_PID(I)
     ENDIF
   ENDDO
   DEALLOCATE(NODEMARK)

   ALLOCATE(HE_LST(NHE),HE_OWN(NHE))
   HE_LST(1:NHE) =  NTEMP(1:NHE)
   HE_OWN(1:NHE) = NTEMP2(1:NHE)
   DEALLOCATE(NTEMP,NTEMP2)
   IF(MSR)WRITE(*,*)  '!  HALO ELEMENTS         :    COMPLETE'

!==============================================================================|
!   DETERMINE LOCAL MAPPING OF GLOBAL ELMNTS (INTERNAL ONLY)     ELID          |
!   ELID(I) = LOCAL ELEMENT ID OF GLOBAL ELEMENT I                             |
!             RETURNS 0 IF I IS NOT AN INTERIOR ELEMENT                        |
!                                                                              |
!   DETERMINE LOCAL MAPPING OF GLOBAL ELMNTS (INTERNAL + HALO)   ELID_X        |
!   ELID_X(I) = LOCAL ELEMENT ID OF GLOBAL ELEMENT I                           |
!             RETURNS 0 IF I IS NOT AN INTERIOR OR HALO ELEMENT                |
!==============================================================================|

   ALLOCATE(ELID(0:NGL)) ; ELID = 0
   ALLOCATE(ELID_X(0:NGL)) ; ELID = 0

   DO I=1,NLOC
     ELID(EGID(I)) = I
   ENDDO
   ELID_X = ELID

   DO I=1,NHE
     ELID_X(HE_LST(I)) = I+NLOC
   ENDDO
   
   
!==============================================================================|
!   DETERMINE NUMBER OF INTERNAL BOUNDARY NODES             NBN                |
!   DETERMINE GLOBAL MAPPING OF INTERNAL BOUNDARY NODES     BN_LST(NBN)        | 
!   DETERMINE LOCAL MAPPING OF INTERNAL BOUNDARY NODES      BN_LOC(NBN)        |
!   DETERMINE MULTIPLICITY (HOW MANY OWNERS) OF IBNs        BN_MLT(NBN)        |
!   MARK OTHER OWNERS OF IBNs                               BN_NEY(NBN,NPROCS) | 
!==============================================================================|

   ALLOCATE(NTEMP(MLOC)) ; NTEMP = 0
   ALLOCATE(NTEMP3(MLOC,NPROCS)) ; NTEMP3 = 0
   ALLOCATE(NDE_ID(MLOC)) ; NDE_ID = 0

   DO I=1,MLOC  
     DO J=1,NGL
       NGLOB = NGID(I)
       FCHECK= DBLE(NVG(J,1)-NGLOB)*DBLE(NVG(J,2)-NGLOB)*DBLE(NVG(J,3)-NGLOB)
       IF(ABS(FCHECK) < 1.0 .AND. MYID /= EL_PID(J))THEN
          NTEMP(I) = 1
          NTEMP3(I,EL_PID(J)) = 1
          NDE_ID(I) = 1
       ENDIF
     ENDDO
   ENDDO

   NBN = SUM(NTEMP)
   ALLOCATE(BN_LST(NBN))
   ALLOCATE(BN_MLT(NBN))
   ALLOCATE(BN_LOC(NBN))
   ALLOCATE(BN_NEY(NBN,NPROCS))
   BN_NEY = 0 ; BN_MLT = 0 ; BN_LST = 0
   NN = 1
   DO I=1,MLOC
     IF(NTEMP(I)==1)THEN
       BN_LST(NN) = NGID(I) 
       BN_LOC(NN) = I
       BN_NEY(NN,MYID) = 1
       DO J=1,NPROCS
         IF(NTEMP3(I,J)==1) BN_NEY(NN,J) = 1
       ENDDO
       BN_MLT(NN) = SUM(NTEMP3(I,:))
       NN = NN + 1
     ENDIF
   ENDDO
   BN_MLT = BN_MLT + 1
   MX_MLT = 1 
   IF(NBN > 0)MX_MLT = MAXVAL(BN_MLT)
   DEALLOCATE(NTEMP,NTEMP3)
   IF(MSR)WRITE(*,*)  '!  BOUNDARY NODES        :    COMPLETE'
   
!==============================================================================|
!   DETERMINE NUMBER OF HALO NODES                               NHN           |
!   DETERMINE GLOBAL NUMBERING OF HALO NODES                     HN_LST(NHN)   |
!   DETERMINE PRIMARY OWNER OF HALO NODE                         HN_OWN(NHN)   !
!==============================================================================|

   ALLOCATE(NTEMP(MGL)) ; NTEMP = 0
   NHN = 0
   DO I=1,MGL
     DO J=1,NHE
       I1 = NVG(HE_LST(J),1)  
       I2 = NVG(HE_LST(J),2)  
       I3 = NVG(HE_LST(J),3)  
       IF(I1 == I .OR. I2 == I .OR. I3 == I)THEN
         INELEM = .TRUE.
         ISMYN  = .FALSE.
         DO K=1,MLOC
           NGLOB=NGID(K)
           IF(I == NGLOB) ISMYN = .TRUE. 
         ENDDO
         IF(INELEM .AND. .NOT.ISMYN)THEN
           NHN = NHN + 1
           NTEMP(NHN) = I
           EXIT
         ENDIF
       ENDIF
     ENDDO
   ENDDO
   ALLOCATE(HN_LST(NHN)) ; HN_LST = NTEMP(1:NHN)
   DEALLOCATE(NTEMP)

   ALLOCATE(HN_OWN(NHN))
   HN_OWN = 0
   DO I=1,NHN
     IGLID = HN_LST(I)
     DO J=1,NHE
       I1 = NVG(HE_LST(J),1)  
       I2 = NVG(HE_LST(J),2)  
       I3 = NVG(HE_LST(J),3)  
       IF(I1 == IGLID .OR. I2 == IGLID .OR. I3 == IGLID)THEN
         HN_OWN(I) = HE_OWN(J)
         EXIT
       ENDIF
     ENDDO
     IF(HN_OWN(I) == 0)THEN
       WRITE(*,*)'NO OWNER FOUND FOR HALO NODE',I,'OF DOMAIN',MYID
     ENDIF
   ENDDO

   IF(MSR)WRITE(*,*)  '!  HALO NODES            :    COMPLETE'

!==============================================================================|
!   DETERMINE LOCAL MAPPING OF GLOBAL NODES  (INTERNAL ONLY)     NLID          |
!   DETERMINE LOCAL MAPPING OF GLOBAL NODES  (INTERNAL + HALO)   NLID_X        |
!==============================================================================|

   ALLOCATE(NLID(0:MGL)) ; NLID = 0
   ALLOCATE(NLID_X(0:MGL)) ; NLID_X = 0

   DO I=1,MLOC
     NLID(NGID(I)) = I
   ENDDO
   NLID_X = NLID

   DO I=1,NHN
     NLID_X(HN_LST(I)) = I+MLOC
   ENDDO
   
         
!==============================================================================|
!   SEND INFORMATION TO PROCESSOR ZERO FOR STATISTICS OUTPUT                   |
!   PNE  :: NUMBER OF ELEMENTS IN EACH DOMAIN                                  |
!   PNN  :: NUMBER OF NODES IN EACH DOMAIN                                     |
!   PNHE :: NUMBER OF HALO ELEMENTS IN EACH DOMAIN                             |
!   PNBN :: NUMBER OF BOUNDARY NODES IN EACH DOMAIN                            |
!   PMBM :: MAXIMUM MULTIPLICITY OF BOUNDARY NODES IN EACH DOMAIN              |
!   PNHN :: NUMBER OF HALO NODES IN EACH DOMAIN                                |
!==============================================================================|

   ALLOCATE(PNE(NPROCS))
   ALLOCATE(PNN(NPROCS))
   ALLOCATE(PNHE(NPROCS))
   ALLOCATE(PNBN(NPROCS))
   ALLOCATE(PMBM(NPROCS))
   ALLOCATE(PNHN(NPROCS))

   CALL MPI_GATHER(NLOC,1,MPI_INTEGER,PNE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(MLOC,1,MPI_INTEGER,PNN,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   IF(MSR)THEN
   WRITE(*,*)'         PROC     #ELEMENTS   IMBALANCE      #NODES  IMBALANCE'
   NEAVE = FLOAT(NGL)/FLOAT(NPROCS)
   NNAVE = FLOAT(SUM(PNN(1:NPROCS)))/FLOAT(NPROCS) !average number of local nodes per proc
   DO I=1,NPROCS
     WRITE(*,'(I12,I14,F12.6,I14,F12.6)')I,PNE(I),FLOAT(PNE(I))/NEAVE,PNN(I),FLOAT(PNN(I))/NNAVE
   ENDDO
   ENDIF
   CALL MPI_GATHER(NHE,1,MPI_INTEGER,PNHE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(NBN,1,MPI_INTEGER,PNBN,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(MX_MLT,1,MPI_INTEGER,PMBM,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(NHN,1,MPI_INTEGER,PNHN,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   IF(MSR)WRITE(*,*)'         PROC     #HALO ELMNTS   #BNDRY NODES  MAX MLTPLCTY  #HALO NODES'
   DO I=1,NPROCS
     IF(MSR)WRITE(*,'(I12,4I14)')I,PNHE(I),PNBN(I),PMBM(I),PNHN(I)
   ENDDO


!==============================================================================|
!   SET UP ELEMENT MAPPING FOR GLOBAL REASSEMBLY OF ARRAYS (MASTER ONLY)       | 
!   ELEMENT MAP :: EMAP(NPROCS)                                                |
!     EMAP(1-->NPROCS)%NSIZE  :: NUMBER OF ELEMENTS IN EACH DOM                |
!     EMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM  |
!==============================================================================|

   ALLOCATE(EMAP(NPROCS))  ; EMAP(:)%NSIZE = 0
   IF(MSR)THEN

!--Determine Number of Elements for Each Processor
   DO I=1,NGL    
     EMAP(EL_PID(I))%NSIZE = EMAP(EL_PID(I))%NSIZE + 1
   ENDDO

!--Allocate Mapping Array for Each Processor
   DO I=1,NPROCS
     ALLOCATE(EMAP(I)%LOC_2_GL(0:EMAP(I)%NSIZE))
     EMAP(I)%LOC_2_GL(0) = 0
   ENDDO
   
!--Construct Mapping Array for Each Processor 
   ALLOCATE(NTEMP(NPROCS))
   NTEMP = 0
   DO I=1,NGL
     I1 = EL_PID(I)
     NTEMP(I1) = NTEMP(I1) + 1
     EMAP(I1)%LOC_2_GL(NTEMP(I1)) = I
   ENDDO
   DEALLOCATE(NTEMP)
   ENDIF

!==============================================================================|
!   TEST ELEMENT MAPPING USING GATHER                                          | 
!==============================================================================|
   
   IF(MSR)THEN
     ALLOCATE(TESTA(NLOC,KB))
     ALLOCATE(TESTB(NGL,KB))
   ELSE
     ALLOCATE(TESTA(NLOC,KB))
     ALLOCATE(TESTB(1,1))
   ENDIF

   TESTA = FLOAT(MYID)
   CALL GATHER(LBOUND(TESTA,1),UBOUND(TESTA,1),NLOC,NGL,KB,MYID,NPROCS,EMAP,TESTA,TESTB) 
   IF(MSR)THEN
     DO J=1,KB
       DO I=1,NGL
         IF(TESTB(I,J) /= EL_PID(I))THEN
           WRITE(*,*)'GLOBAL ARRAY NOT CORRECT',I,J,TESTB(I,J),EL_PID(I)
         ENDIF
       ENDDO
     ENDDO
   ENDIF
   DEALLOCATE(TESTA,TESTB)

    
   IF(MSR)WRITE(*,*)  '!  ELEMENT GATHER TEST   :    PASSED'
       
   
!==============================================================================|
!   SET UP NODE MAPPING FOR GLOBAL REASSEMBLY OF ARRAYS (MASTER ONLY)          | 
!   NODAL   MAP :: NMAP(NPROCS)                                                |
!     NMAP(1-->NPROCS)%NSIZE  :: NUMBER OF NODES IN EACH DOM                   |
!     NMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM  |
!==============================================================================|

   ALLOCATE(NMAP(NPROCS))  ; NMAP(:)%NSIZE = 0
   IF(MSR)THEN
   NMAP(:)%NSIZE = PNN(:)

   DO I=1,NPROCS
     ALLOCATE(NMAP(I)%LOC_2_GL(0:NMAP(I)%NSIZE))
     NMAP(I)%LOC_2_GL(0) = 0
   ENDDO
   ENDIF
  
   IF(.NOT. MSR)THEN
     ISND = MYID+1000
     CALL MPI_SEND(NGID,MLOC,MPI_INTEGER,0,ISND,MPI_COMM_WORLD,IERR)
   ENDIF

   IF(MSR)THEN
     NMAP(1)%LOC_2_GL(1:MLOC) = NGID(1:MLOC)
     DO I=2,NPROCS
       NSZE = NMAP(I)%NSIZE
       ALLOCATE(NTEMP(NSZE))
       IRCV = I+1000
       CALL MPI_RECV(NTEMP,NSZE,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
       NMAP(I)%LOC_2_GL(1:NSZE) = NTEMP
       DEALLOCATE(NTEMP)
     ENDDO
   ENDIF
!WL   debug3602 = -1 ! not found yet
!WL   do j=1,MLOC
!WL      if(NGID(j) .eq. 3602) debug3602=j
!WL   enddo

!==============================================================================|
!   TEST NODE MAPPING USING GATHER                                             | 
!==============================================================================|

   ALLOCATE(TESTA(MLOC,KB))
   ALLOCATE(TESTB(MGL,KB))


   TESTA = 100.0
   CALL GATHER(LBOUND(TESTA,1),UBOUND(TESTA,1),MLOC,MGL,KB,MYID,NPROCS,NMAP,TESTA,TESTB)
   IF(MSR)THEN
     DO J=1,KB
       DO I=1,MGL
         IF(TESTB(I,J) /= 100.0)THEN
           WRITE(*,*)'GLOBAL ARRAY NOT CORRECT',I,J,TESTB(I,J)
         ENDIF
       ENDDO
     ENDDO
   ENDIF
   DEALLOCATE(TESTA,TESTB)

   IF(MSR)WRITE(*,*)  '!  NODE GATHER TEST      :    PASSED'
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)



!==============================================================================|
!   SET UP COMMUNICATION: HALO ELEMENTS                                        |
!==============================================================================|

   ALLOCATE(EC(NPROCS))


   EC(:)%NRCV = 0
   ALLOCATE(NTEMP3(NPROCS,NGL))

!--Set up Recieve Information First
   DO I=1,NHE
     OWNER = HE_OWN(I)
     EC(OWNER)%NRCV = EC(OWNER)%NRCV + 1
     NTEMP3(OWNER,EC(OWNER)%NRCV) = I + NLOC
   ENDDO
   
   NN = 0
   DO I=1,NPROCS
     EC(I)%RCPT = NN
     NSZE = EC(I)%NRCV 
     NN = NN + NSZE
!     ALLOCATE(EC(I)%RCVP(NSZE))
     IF(NSZE > 0)ALLOCATE(EC(I)%RCVP(NSZE))
     IF(NSZE == 0)ALLOCATE(EC(I)%RCVP(1))
     IF(NSZE > 0)EC(I)%RCVP = NTEMP3(I,1:NSZE)
   ENDDO
   DEALLOCATE(NTEMP3)
   TERCV = SUM(EC(:)%NRCV)
   IF(TERCV /= NHE)THEN
     WRITE(*,*)'TOTAL NUMBER OF ELEMENTS SET UP TO RECEIVE DATA'
     WRITE(*,*)'NOT EQUAL TO TOTAL NUMBER OF HALO ELEMENTS'
     WRITE(*,*)TERCV,NHE
     CALL MPI_FINALIZE(IERR)
     CALL PSTOP
   ENDIF

!--Set up Send
   EC(:)%NSND = 0

   DO I=1,NPROCS
   IF(MYID /=I)THEN
     ISND = MYID + 1000*I
     CALL MPI_SEND(EC(I)%NRCV,1,MPI_INTEGER,I-1,ISND,MPI_COMM_WORLD,IERR)
   ENDIF
   ENDDO
   DO I=1,NPROCS
   IF(MYID /=I)THEN
     IRCV = 1000*MYID + I
     CALL MPI_RECV(EC(I)%NSND,1,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
   ENDIF
!   ALLOCATE(EC(I)%SNDP(EC(I)%NSND))
   IF(EC(I)%NSND > 0)ALLOCATE(EC(I)%SNDP(EC(I)%NSND))
   IF(EC(I)%NSND == 0)ALLOCATE(EC(I)%SNDP(1))
   ENDDO
   DO I=1,NPROCS
   IF(MYID /=I .AND. EC(I)%NRCV > 0)THEN
     ALLOCATE(NTEMP(EC(I)%NRCV))
     DO J=1,EC(I)%NRCV
       NTEMP(J) = HE_LST(EC(I)%RCVP(J)-NLOC)
     ENDDO
     ISND = MYID + 2000*I
     CALL MPI_SEND(NTEMP,EC(I)%NRCV,MPI_INTEGER,I-1,ISND,MPI_COMM_WORLD,IERR)
     DEALLOCATE(NTEMP)
   ENDIF
   ENDDO
   DO I=1,NPROCS
   IF(MYID /=I .AND. EC(I)%NSND > 0)THEN
     ALLOCATE(NTEMP(EC(I)%NSND))
     IRCV = 2000*MYID + I
     CALL MPI_RECV(NTEMP,EC(I)%NSND,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
     N1 = 0
     DO J=1,EC(I)%NSND
     DO K=1,NLOC
       IF(EGID(K)==NTEMP(J))THEN
         N1 = N1 + 1
         EC(I)%SNDP(N1) = K
       ENDIF
     ENDDO
     ENDDO
     DEALLOCATE(NTEMP)
   ENDIF
   ENDDO
   

!==============================================================================|
!   CHECK ELEMENT COMMUNICATION                                                |
!==============================================================================|
!--Set Up Testing Array 
   ALLOCATE(TESTA(0:NLOC+NHE,KB))  ; TESTA = 0.0
   ALLOCATE(TESTB(0:NLOC+NHE,KB))  ; TESTB = 0.0
   ALLOCATE(TESTC(0:NLOC+NHE,KB))  ; TESTC = 0.0

!--Initialize with Function of Global Element Index in Interior, Zero in Halo
   DO I=1,KB
   TESTA(1:NLOC,I) =  (FLOAT(EGID(1:NLOC))  +10000*I)
   TESTB(1:NLOC,I) = -(FLOAT(EGID(1:NLOC))  +10000*I)
   TESTC(1:NLOC,I) =  (FLOAT(EGID(1:NLOC))  +10000*I) + .5
   ENDDO
   

!--Perform Communication

   CALL EXCHANGE(EC,NLOC+NHE,KB,MYID,NPROCS,TESTA,TESTB,TESTC) 

!--Check Results For Interior on Array A
   DO I=1,NLOC
   DO J=1,KB
   VAL  = FLOAT(EGID(I)+10000*J)  
   DIFF = ABS(TESTA(I,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'INTERNAL ELEMENT VALUE CHANGED on A',I,J,TESTA(I,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO

!--Check Results For Halo on Array A 
   DO I=1,NHE
   DO J=1,KB
   VAL = FLOAT(HE_LST(I)+10000*J)
   DIFF = ABS(TESTA(I+NLOC,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'HALO ELEMENT CHANGED on A',I,J,TESTA(I+NLOC,J),VAL,I+NLOC
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO
   

!--Check Results For Interior on Array B
   DO I=1,NLOC
   DO J=1,KB
   VAL  = -FLOAT(EGID(I)+10000*J)  
   DIFF = ABS(TESTB(I,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'INTERNAL ELEMENT VALUE CHANGED on B',I,J,TESTB(I,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO

!--Check Results For Halo on Array B 
   DO I=1,NHE
   DO J=1,KB
   VAL = -FLOAT(HE_LST(I)+10000*J)
   DIFF = ABS(TESTB(I+NLOC,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'HALO ELEMENT CHANGED on B',I,J,TESTB(I+NLOC,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO

!--Check Results For Interior on Array C
   DO I=1,NLOC
   DO J=1,KB
   VAL  = FLOAT(EGID(I)+10000*J)+.5
   DIFF = ABS(TESTC(I,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'INTERNAL ELEMENT VALUE CHANGED on C',I,J,TESTC(I,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO

!--Check Results For Halo on Array C 
   DO I=1,NHE
   DO J=1,KB
   VAL = FLOAT(HE_LST(I)+10000*J)+.5
   DIFF = ABS(TESTC(I+NLOC,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'HALO ELEMENT CHANGED on C',I,J,TESTC(I+NLOC,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO
  
   DEALLOCATE(TESTA,TESTB,TESTC)
   
   IF(MSR)WRITE(*,*)  '!  ELEMENT COMM TEST     :    PASSED'
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

!==============================================================================|
!   SET UP COMMUNICATION: BOUNDARY NODES                                       |
!==============================================================================|

   ALLOCATE(NTEMP3(NBN,NPROCS),NTEMP4(NBN,NPROCS))
   ALLOCATE(BNC(NPROCS)) 
   BNC%NRCV = 0      ; BNC%NSND = 0

!--Count Receives from Each Processor
   DO I=1,NBN
     DO J=1,NPROCS
       IF(BN_NEY(I,J)==1 .AND. MYID /= J) THEN
         BNC(J)%NRCV = BNC(J)%NRCV + 1
         NTEMP3(BNC(J)%NRCV,J) = BN_LOC(I) 
         NTEMP4(BNC(J)%NRCV,J) = BN_MLT(I) 
       ENDIF
     ENDDO
   ENDDO

!--Sent up Indices for Receive Buffer
   NN = 0
   DO I=1,NPROCS
     BNC(I)%RCPT = NN
     NSZE = BNC(I)%NRCV 
     NN = NN + NSZE
   ENDDO

!--Sends = Receives for Boundary Nodes
   BNC%NSND = BNC%NRCV

!--Array of Receives/Sends/Multiplicities for Each Processor
   DO I=1,NPROCS
     NSZE = BNC(I)%NRCV
     ALLOCATE(BNC(I)%RCVP(NSZE))
     ALLOCATE(BNC(I)%SNDP(NSZE))
     ALLOCATE(BNC(I)%MLTP(NSZE))
     BNC(I)%RCVP = NTEMP3(1:NSZE,I)
     BNC(I)%SNDP = NTEMP3(1:NSZE,I)
     BNC(I)%MLTP = NTEMP4(1:NSZE,I)
   ENDDO
   DEALLOCATE(NTEMP3)


!--Sort Send/Recive Arrays According to Global Identity to Match Order    
   
   DO I=1,NPROCS
     IF(BNC(I)%NRCV > 0)THEN
       NSZE     = BNC(I)%NRCV
       ALLOCATE(NTEMP(NSZE),NTEMP2(NSZE))
       NTEMP(:) = NGID(BNC(I)%RCVP(:))
       CALL SORT(NTEMP,NTEMP2,NSZE)
       NTEMP = BNC(I)%RCVP
       DO J=1,NSZE
         BNC(I)%RCVP(J) = NTEMP(NTEMP2(J))
       ENDDO
       NTEMP = BNC(I)%MLTP
       DO J=1,NSZE
         BNC(I)%MLTP(J) = NTEMP(NTEMP2(J))
       ENDDO
       BNC(I)%SNDP = BNC(I)%RCVP
       DEALLOCATE(NTEMP,NTEMP2)
     ENDIF
   ENDDO
 
 
!==============================================================================|
!   TEST COMMUNICATION: HALO NODES                                             |
!==============================================================================|

!--Set Up Testing Array 
   ALLOCATE(TESTA(0:MLOC+NHN,KB))  
   ALLOCATE(TESTB(0:MLOC+NHN,KB))  
   ALLOCATE(TESTC(0:MLOC+NHN,KB))  
   TESTA = FLOAT(MYID)  
   TESTB = FLOAT(MYID)  * 100.0
   TESTC = FLOAT(MYID)  * 10000.0
   DO J=1,KB
   DO I=1,NBN
     TESTA(BN_LOC(I),J) = J
     TESTB(BN_LOC(I),J) = 100*J
     TESTC(BN_LOC(I),J) = 10000*J
   ENDDO
   ENDDO

!--Perform Communication
   
   CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MLOC+NHN,KB,MYID,NPROCS,TESTA,TESTB,TESTC) 
!

!--Check Results For Interior and Boundary on Array A
   DO J=1,KB
   DO I=1,MLOC
     IF(NDE_ID(I)==0)THEN
     IF( ABS(TESTA(I,J) -  FLOAT(MYID))> .01 )THEN
       WRITE(*,*)'INTERNAL NODE CHANGED',I,J,TESTA(I,J),FLOAT(MYID)
       CALL PSTOP
     ENDIF
     ELSE
     IF(ABS(TESTA(I,J) - J) > .01)THEN
       WRITE(*,*)'BOUNDARY PROBLEM',MYID,I,J,TESTA(I,J)
       WRITE(*,*)'GLOBAL NODE: ',NGID(I)
       CALL PSTOP
     ENDIF
     ENDIF
   ENDDO
   ENDDO
   
!--Check Results For Interior and Boundary on Array B
   DO J=1,KB
   DO I=1,MLOC
     IF(NDE_ID(I)==0)THEN
     IF( ABS(TESTB(I,J) -  FLOAT(MYID)*100.0)> .01 )THEN
       WRITE(*,*)'INTERNAL NODE CHANGED',I,J,TESTB(I,J),FLOAT(MYID)
       CALL PSTOP
     ENDIF
     ELSE
     IF(ABS(TESTB(I,J) - 100*J) > .01)THEN
       WRITE(*,*)'BOUNDARY PROBLEM',MYID,I,J,TESTB(I,J)
       CALL PSTOP
     ENDIF
     ENDIF
   ENDDO
   ENDDO

!--Check Results For Interior and Boundary on Array C
   DO J=1,KB
   DO I=1,MLOC
     IF(NDE_ID(I)==0)THEN
     IF( ABS(TESTC(I,J) -  FLOAT(MYID)*10000.0)> .01 )THEN
       WRITE(*,*)'INTERNAL NODE CHANGED',I,J,TESTC(I,J),FLOAT(MYID)
       CALL PSTOP
     ENDIF
     ELSE
     IF(ABS(TESTC(I,J) - 10000.0*J) > .01)THEN
       WRITE(*,*)'BOUNDARY PROBLEM',MYID,I,J,TESTC(I,J)
       CALL PSTOP
     ENDIF
     ENDIF
   ENDDO
   ENDDO

   DEALLOCATE(TESTA,TESTB,TESTC)


   IF(MSR)WRITE(*,*)'!  BNDRY NODE COMM TEST  :    PASSED'

!==============================================================================|
!   SET UP COMMUNICATION: HALO NODES                                           |
!==============================================================================|
! NHN,HN_LST --HN_OWN (owner of halo node i, like he_own)

   ALLOCATE(NC(NPROCS))
   NC(:)%NRCV = 0

   ALLOCATE(NTEMP3(NPROCS,MGL))

!--Set up Recieve Information First
   DO I=1,NHN
     OWNER = HN_OWN(I)
     NC(OWNER)%NRCV = NC(OWNER)%NRCV + 1
     NTEMP3(OWNER,NC(OWNER)%NRCV) = I + MLOC
   ENDDO
   
   NN = 0
   DO I=1,NPROCS
     NC(I)%RCPT = NN
     NSZE = NC(I)%NRCV 
     NN = NN + NSZE
     ALLOCATE(NC(I)%RCVP(NSZE))
     NC(I)%RCVP = NTEMP3(I,1:NSZE)
   ENDDO
   DEALLOCATE(NTEMP3)
   TERCV = SUM(NC(:)%NRCV)
   IF(TERCV /= NHN)THEN
     WRITE(*,*)'TOTAL NUMBER OF NODES SET UP TO RECEIVE DATA'
     WRITE(*,*)'NOT EQUAL TO TOTAL NUMBER OF HALO NODES'
     WRITE(*,*)TERCV,NHN
     CALL MPI_FINALIZE(IERR)
     CALL PSTOP
   ENDIF

!--Set up Send
   NC(:)%NSND = 0

   DO I=1,NPROCS
   IF(MYID /=I)THEN
     ISND = MYID + 1000*I
     CALL MPI_SEND(NC(I)%NRCV,1,MPI_INTEGER,I-1,ISND,MPI_COMM_WORLD,IERR)
   ENDIF
   ENDDO
   DO I=1,NPROCS
   IF(MYID /=I)THEN
     IRCV = 1000*MYID + I
     CALL MPI_RECV(NC(I)%NSND,1,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
   ENDIF
   ALLOCATE(NC(I)%SNDP(NC(I)%NSND))
   ENDDO

   DO I=1,NPROCS
   IF(MYID /=I .AND. NC(I)%NRCV > 0)THEN
     ALLOCATE(NTEMP(NC(I)%NRCV))
     DO J=1,NC(I)%NRCV
       NTEMP(J) = HN_LST(NC(I)%RCVP(J)-MLOC)
     ENDDO
     ISND = MYID + 1000*I
     CALL MPI_SEND(NTEMP,NC(I)%NRCV,MPI_INTEGER,I-1,ISND,MPI_COMM_WORLD,IERR)
     DEALLOCATE(NTEMP)
   ENDIF
   ENDDO

   DO I=1,NPROCS
   IF(MYID /=I .AND. NC(I)%NSND > 0)THEN
     ALLOCATE(NTEMP(NC(I)%NSND))
     IRCV = 1000*MYID + I
     CALL MPI_RECV(NTEMP,NC(I)%NSND,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
     N1 = 0
     DO J=1,NC(I)%NSND
     DO K=1,MLOC
       IF(NGID(K)==NTEMP(J))THEN
         N1 = N1 + 1
         NC(I)%SNDP(N1) = K
       ENDIF
     ENDDO
     ENDDO
     DEALLOCATE(NTEMP)
   ENDIF
   ENDDO
   
!==============================================================================|
!   CHECK HALO NODE COMMUNICATION                                              |
!==============================================================================|

!--Set Up Testing Array 
   ALLOCATE(TESTA(0:MLOC+NHN,KB))  ; TESTA = 0.0
   ALLOCATE(TESTB(0:MLOC+NHN,KB))  ; TESTB = 0.0
   ALLOCATE(TESTC(0:MLOC+NHN,KB))  ; TESTC = 0.0

!--Initialize with Global Nodal ID in Interior, Zero in Halo
   DO I=1,KB
   TESTA(1:MLOC,I) =   FLOAT(NGID)  +10000*I
   TESTB(1:MLOC,I) = -(FLOAT(NGID)  +10000*I)
   TESTC(1:MLOC,I) = (FLOAT(NGID)  +10000*I) + .50
   ENDDO

!--Perform Communication
   CALL EXCHANGE(NC,MLOC+NHN,KB,MYID,NPROCS,TESTA,TESTB,TESTC) 
!--Check Results For Interior on Array A
   DO I=1,MLOC
   DO J=1,KB
   VAL  = FLOAT(NGID(I)+10000*J)  
   DIFF = ABS(TESTA(I,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'INTERNAL NODE VALUE CHANGED on A',I,J,TESTA(I,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO

!--Check Results For Halo on Array A 
   DO I=1,NHN
   DO J=1,KB
   VAL = FLOAT(HN_LST(I)+10000*J)
   DIFF = ABS(TESTA(I+MLOC,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'HALO NODE CHANGED on A',I,J,TESTA(I+MLOC,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO
   

!--Check Results For Interior on Array B
   DO I=1,MLOC
   DO J=1,KB
   VAL  = -FLOAT(NGID(I)+10000*J)  
   DIFF = ABS(TESTB(I,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'INTERNAL NODE VALUE CHANGED on B',I,J,TESTB(I,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO

!--Check Results For Halo on Array B 
   DO I=1,NHN
   DO J=1,KB
   VAL = -FLOAT(HN_LST(I)+10000*J)
   DIFF = ABS(TESTB(I+MLOC,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'HALO NODE CHANGED on B',I,J,TESTB(I+MLOC,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO

!--Check Results For Interior on Array C
   DO I=1,MLOC
   DO J=1,KB
   VAL  = FLOAT(NGID(I)+10000*J)+.5
   DIFF = ABS(TESTC(I,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'INTERNAL NODE VALUE CHANGED on C',I,J,TESTC(I,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO

!--Check Results For Halo on Array C 
   DO I=1,NHN
   DO J=1,KB
   VAL = FLOAT(HN_LST(I)+10000*J)+.5
   DIFF = ABS(TESTC(I+MLOC,J)-VAL)
   IF(DIFF > .01)THEN
   PRINT *,MYID,'HALO ELEMENT CHANGED on C',I,J,TESTC(I+MLOC,J),VAL
   CALL PSTOP
   ENDIF
   ENDDO
   ENDDO
  
   DEALLOCATE(TESTA,TESTB,TESTC)
   
   IF(MSR)WRITE(*,*)  '!  NODE COMM TEST        :    PASSED'

!==============================================================================|
!   SET TOTAL FOR EACH PROCESSOR                                               |
!   NTLOC:   TOTAL NUMBER OF ELEMENTS = INTERIOR (NLOC) + HALO (NHE)                 |
!   MTLOC:   TOTAL NUMBER OF NODES    = INTERIOR (MLOC) + HALO (NHN)                 |
!==============================================================================|
   NTLOC = NLOC + NHE
   MTLOC = MLOC + NHN
  
   RETURN
# endif
   END SUBROUTINE GENMAP
!==============================================================================!
