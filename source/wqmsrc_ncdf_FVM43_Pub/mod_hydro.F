MODULE MOD_HYDRO

        CONTAINS
        !subroutine HYDRO ()

        !WLong: This should be combined with mod_hydrovars.F

        !************************************************************************
        !**                  S U B R O U T I N E   H Y D R O                   **
        !************************************************************************

            SUBROUTINE HYDRO ()
            USE MOD_PREC, ONLY: SP
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
                          UNC2	&		!
                         ,VNC2	&		!
                        !,WNC2	&		!
                         ,WTSNC2	&	!
                         ,UARD_OBCNNC2	&!
                         ,XFLUX_OBCNC2	&!
                         ,DTFANC2	&	!
                         ,KHNC2	&		!
                         ,TNC2	&		!
                         ,SNC2	&		!
                         ,ELNC2	&		!
                         ,num_hyd_ints	&!number of records in each hydrodynamics netcdf file
                        !,TIME_MAP	&	!
                        !,THOUR1	&	!SIMULATION TIME AT END OF CURRENT EXTERNAL STEP (IEXT) IN HOURS
                        !,THOUR	&		!
                        !,NCFILE_DIR	&!
                        ,NCFILE_PREFIX	&!
                        ,NCFILE_SUFFIX	&!
                        ,NCFILE_NUMBER	&!
                        ,FORMAT_STR	&!
                        ,hydro_dir 	&	! directory name where hydrodynamics results (netcdf) files are stored
                        !,hydro_prefix  &	! prefix of file name, e.g. 'psm_'
                        !,hydro_suffix	&	! suffix of filename, e.g. '.nc'
                        ,hydro_filenumwidth &	! number of digits in filename following hydro_prefix, e.g. 4 for psm_0002.nc
                        ,hydro_filenumstart &	! starting number of the file name in the digital part of the file name, e.g. 185 for psm_0185.nc
                        ,hydro_Nrec	&		! number of records in each of hydrodynamics file
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
                        ,IFNC	&				!file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise 
                        ,NTRECNC	&			!time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file. 
                        ,NTHYDRO				!overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read

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

                USE MOD_LIMS, ONLY: MTLOC, KBM1
            USE MOD_SIZES, ONLY :        &	!
                                NCP!,            &  !
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
                                !NGL!,            &  !
                                !OBCGL,          &  !
                                !NOBTY  
          
              USE MOD_NCDIO
            IMPLICIT NONE
#if defined (MULTIPROCESSOR)
                include "mpif.h"
#endif
            SAVE
            INTEGER   F, SB, L, JCON, I, J, NC_ID2
            LOGICAL   END_OF_FILE
            REAL(SP)   NXDAY, TDUM
            REAL(SP)   MASS(0:MTLOC,KBM1,NCP)  !WLong: this is unncessary variable
            CHARACTER(LEN=1024) :: NCFILE

        !************************************************************************
        !**                              Inputs                                **
        !************************************************************************

        ! Wen Long moved these to the main program wam_main.F 
        !    UNC1  = UNC2
        !    VNC1  = VNC2
        !    WNC1  = WNC2
        !    WTSNC1 = WTSNC2
        !    UARD_OBCNNC1=UARD_OBCNNC2
        !    XFLUX_OBCNC1=XFLUX_OBCNC2
        !    DTFANC1=DTFANC2
        !    KHNC1 = KHNC2
        !    ELNC1 = ELNC2
        !    TNC1  = TNC2
        !    SNC1  = SNC2

        !**************RGL section begins below**********************************
        ! RGl need to change below harddwire if change number of tsteps/file netcdf
          num_hyd_ints = hydro_Nrec

        ! find indices for the next record to read 
        !         index of file:                        IFNC 
        !         time record index within the file:    NTRECNC
        !         overwall hydrodynamics record index:  NTHYDRO

          IF(NTRECNC.eq.num_hyd_ints) THEN
             NTRECNC = 0
             IFNC = IFNC +1       
          ENDIF
          NTRECNC=NTRECNC+1     !NTRECNC is time record # within one hydro file
          NTHYDRO=NTHYDRO+1     !NTHYDRO is overall record # of all hydro files that have been used
         
              IF(NTRECNC.eq.1)THEN  !open new file (IFNC) when trying to read the first record
                                    
                 NCFILE_NUMBER=''
                 WRITE(NCFILE_NUMBER(1:hydro_filenumwidth),TRIM(FORMAT_STR))(IFNC+hydro_filenumstart-1)
                 NCFILE=TRIM(hydro_dir)//TRIM(NCFILE_PREFIX)//TRIM(NCFILE_NUMBER)//TRIM(NCFILE_SUFFIX)
                 IF(MSR) WRITE(*,*)'NCFILE=',TRIM(NCFILE)

                                                           !        variable             varilable                         unit       dimensions (e.g. :
                                                           !        name in              meaning                                      siglay=10, siglev=11
                                                           !        NCFILE                                                            obc=11,obc2=11,node=9013
                                                           !                                                                          nele=13941 )
            !                                    ----------------------------------------------------------------------------------------------------------
                 CALL NCD_READ_OPEN(NCFILE,             &  !        --NetCDF file name 
                                      UNC2,             &  !        --'u'                 "Eastward Water Velocity"        m/s        (time,siglay,nele)
                                      VNC2,             &  !        --'v'                 "Northward Water Velocity"       m/s        (time,siglay,nele)
                                    WTSNC2,             &  !        --'wts'               "Upward Water Velocity at node"  m/s        (time,siglev,node)
                              UARD_OBCNNC2,             &  !        --'uard_obcn'         "UARD at OBC"                    m/s        (time,obc)
                              XFLUX_OBCNC2,             &  !        --'xflux_obc'         "Xflux at OBC"                   m/s***     (time,siglay,obc2)
                                   DTFANC2,             &  !        --'dtfa'              "FSH + Water Height"              m         (time,node)
                                     KHNC2,             &  !        --'kh'                "Turbulent Eddy Diffusivity"     m^2/s      (time,siglev,node)
                                     ELNC2,             &  !        --'zeta'              "Water Surface Elevation"         m         (time,node)
                                      TNC2,             &  !        --'temp'              "temperature"                    deg-C      (time,siglay,node)
                                      SNC2,             &  !        --'salinity'          "salinity"                        psu       (time,siglay,node)
                                   NTRECNC)                !        -- time record number corresponding to index of 'time' variable in NCFILE
            !                                    -----------------------------------------------------------------------------------------------------------
            !Wen Long                           !*** the unit of xflux_obc should be [m][m/s][concentration]  i.e. [m^2/s][concentration] rather than [m/s]!!
            !                                   !    where concentration could be for dye, deg-C for temperature, psu for salinity
            !                                   !    FVCOM is not clear on what unit xflux_obc has, but the calculation gives (velocity * depth * concentration)

              ELSE

                 CALL NCD_READ(NCFILE,UNC2,VNC2,WTSNC2,UARD_OBCNNC2,XFLUX_OBCNC2,DTFANC2,&
                      KHNC2,ELNC2,TNC2,SNC2,NTRECNC)

              ENDIF

            RETURN
            END SUBROUTINE HYDRO


#if defined (MULTIPROCESSOR)
            SUBROUTINE BROADCAST_HYDRO (MPI_SRC,UL,VL,WTSL,UARD_OBCNL,XFLUX_OBCL,DTFAL,KHL,ELL,TL,SL)

            !
            !Subroutine to broadcast hydrodynamic fields to child processes from master process
            !
              

                  USE MOD_LIMS, ONLY :  &
                  NLOC,         &       !!NUMBER OF ELEMENTS 
                  MLOC,         &       !!NUMBER OF NODES
                  !NISBCE_1,     &      !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 1
                  !NISBCE_2,     &      !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 2
                  !NISBCE_3,     &      !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 3
                   KB,           &      !!NUMBER OF SIGMA LEVELS
                  KBM1,         &      !!NUMBER OF SIGMA LEVELS-1
                  !KBM2 ,        &      !!NUMBER OF SIGMA LEVELS-2
                  !MYID ,        &      !!UNIQUE PROCESSOR ID (1 => NPROCS)
                  !NPROCS,       &      !!NUMBER OF PROCESSORS
                  !NE   ,        &      !!NUMBER OF UNIQUE EDGES
                  !NCV  ,        &      !!NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)
                  !IINT  ,       &      !!TYKIM added for nudging
                  !NCV_I ,       &      !!NUMBER OF INTERNAL CONTROL VOLUMES (LOCAL ONLY)
                  NTLOC ,       &      !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
                  MTLOC! ,       &      !!TOTAL OF LOCAL INTERNAL + HALO NODES
                  !NCT   ,       &
                  !MX_NBR_ELEM , &      !!MAX NUMBER OF ELEMENTS SURROUNDING A NODE
                  !NUMQBC_GL ,   &
                          !NUMPNT_GL ,   &
                  !NUMQBC,       &
                          !NUMPNT,       &
                  !NstationMax         !Maximum number of stations 
                             
              USE MOD_WQM, ONLY :  &
                                  UL_GL,   &
                                  VL_GL,   &
								  WW_GL,&
								  WW_OUT, &
                                                  WTSL_GL, &
                                                  UARD_OBCN_GL,  &
                                                  XFLUX_OBC_GL,  &
                                                  DTFAL_GL,   &
                                                  KHL_GL, &
                                                  ELL_GL, &
                                                  SL_GL, &
                                                  TL_GL
                                                
                                                
                  
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
                                NGL,            &  !
                                !OBCGL,          &  !
                                NOBTY  
          
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

                  USE MOD_BCMAP, ONLY: IOBCN, I_OBC_GL
                 
                  
                  USE MOD_PREC, ONLY: MPI_F, SP, MPI_CDF
              IMPLICIT NONE
                  include "mpif.h"
                  
            !passed variables
              INTEGER,INTENT(IN) :: MPI_SRC  ! ID of the source process for broadcasting
              REAL(SP), DIMENSION(0:NTLOC,KB),INTENT(OUT)   :: UL,VL!, WWL
              REAL(SP), DIMENSION(0:MTLOC,KB),INTENT(OUT)   :: KHL,WTSL
              REAL(SP), DIMENSION(0:MTLOC,KBM1),INTENT(OUT) :: TL,SL       !Wen Long changed KB to KBM1 for TL,SL
              REAL(SP), DIMENSION(0:MTLOC),INTENT(OUT)      :: ELL,DTFAL
              REAL(SP), DIMENSION(0:NOBTY+1),INTENT(OUT)    :: UARD_OBCNL
              REAL(SP), DIMENSION(0:NOBTY,KBM1),INTENT(OUT) :: XFLUX_OBCL

            !local variables

               INTEGER :: MPI_COUNT_TMP=0                !count to send
               INTEGER :: I,J,K
               INTEGER :: IERR
            !Broadcast global variables
             IF(PAR)THEN  
        !   KURT GLAESEMANN 14 APRIL 2015 - these are MPI_CDF, not MPI_F precision. Always real*4
                      MPI_COUNT_TMP=NGL*KBM1  ; CALL MPI_BCAST(       UL_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=NGL*KBM1  ; CALL MPI_BCAST(       VL_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=MGL*KB    ; CALL MPI_BCAST(     WTSL_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=NOBTY+1   ; CALL MPI_BCAST(  UARD_OBCN_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=NOBTY*KBM1; CALL MPI_BCAST(XFLUX_OBC_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=MGL       ; CALL MPI_BCAST(      DTFAL_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=MGL*KB    ; CALL MPI_BCAST(      KHL_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=MGL       ; CALL MPI_BCAST(        ELL_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=MGL*KBM1  ; CALL MPI_BCAST(       SL_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=MGL*KBM1  ; CALL MPI_BCAST(       TL_GL,MPI_COUNT_TMP,MPI_CDF,MPI_SRC,MPI_COMM_WORLD,IERR)
             ENDIF

            !---then push the global arrays into local arays including halo nodes and halo elements
             IF(PAR)THEN
                DO I=1,MLOC
                         !WRITE(*,*)'ELL_GL Debug here here 01 '
                  ELL(I) = ELL_GL(NGID(I))
                         ! WRITE(*,*)'ELL_GL Debug here here 1 '
                          
                ENDDO

                DO I=1,NHN

                  !WRITE(*,*)'ELL_GL Debug here here 02'
                          !WRITE(*,*)'I, NHN, MLOC, HN_LST(I) ', I, NHN, MLOC, HN_LST(I)
                  ELL(I+MLOC) = ELL_GL(HN_LST(I))
                          !WRITE(*,*)'ELL_GL Debug here here 2 '
                ENDDO

                DO K=1,KBM1
                  DO I=1,MLOC
                    SL(I,K) = SL_GL(NGID(I),K)
                  ENDDO

                  DO I=1,NHN
                    SL(I+MLOC,K) = SL_GL(HN_LST(I),K)
                  ENDDO
                ENDDO

                DO K=1,KBM1
                  DO I=1,MLOC
                    TL(I,K) = TL_GL(NGID(I),K)
                  ENDDO

                  DO I=1,NHN
                    TL(I+MLOC,K) = TL_GL(HN_LST(I),K)
                  ENDDO
                ENDDO

                DO K=1,KBM1
                  DO I=1,NLOC
                    UL(I,K) = UL_GL(EGID(I),K)
				!	WWL(I,K) = WW_GL(EGID(I),K)
                  ENDDO

                  DO I=1,NHE
                    UL(I+NLOC,K) = UL_GL(HE_LST(I),K)
					!WWL(I+NLOC,K) = WW_GL(EGID(I),K)
                  ENDDO
                ENDDO
				!WW_OUT=WWL

                DO K=1,KBM1
                  DO I=1,NLOC
                    VL(I,K) = VL_GL(EGID(I),K)
                  ENDDO

                  DO I=1,NHE
                    VL(I+NLOC,K) = VL_GL(HE_LST(I),K)
                  ENDDO
                ENDDO

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
                     CALL pstop
                ENDIF

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
                     CALL PSTOP
                ENDIF

                DO I=1,MLOC
                  DTFAL(I) = DTFAL_GL(NGID(I))
                ENDDO

                DO I=1,NHN
                  DTFAL(I+MLOC) = DTFAL_GL(HN_LST(I))
                ENDDO

                DO K=1,KB
                  DO I=1,MLOC
                    KHL(I,K) = KHL_GL(NGID(I),K)
                  ENDDO

                  DO I=1,NHN
                    KHL(I+MLOC,K) = KHL_GL(HN_LST(I),K)
                  ENDDO
                ENDDO

              ENDIF

            RETURN
            END SUBROUTINE BROADCAST_HYDRO
#endif 

#if defined (MULTIPROCESSOR)
            SUBROUTINE BROADCAST_HYDRO_REDUCED(MPI_SRC,UARD_OBCNL,XFLUX_OBCL,DTFAL,ELL)
            !
            !Subroutine to broadcast hydrodynamic fields to child processes from master process
            !

              USE MOD_LIMS, ONLY:  NTLOC, MLOC, MTLOC, KBM1
                  
                  USE MOD_PREC, ONLY : MPI_F, SP
                  
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
                    MGL	,            &  !
                                !NGL,            &  !
                                !OBCGL,          &  !
                                NOBTY  
                  
              USE MOD_WQM, ONLY:		&
                                        UARD_OBCN_GL, 	&
                                        XFLUX_OBC_GL,	&
                                        DTFAL_GL,		&
                                        ELL_GL!,		
                                        
                                        
              
                        !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
                  USE MOD_CONTROL, ONLY : 		&
                                !SERIAL  				&           !!TRUE IF SINGLE PROCESSOR
                                !,MSR        			&           !!TRUE IF MASTER PROCESSOR (MYID==1)
                                PAR        !			&           !!TRUE IF MULTIPROCESSOR RUN
                                !,CASENAME  			&   		!!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
                                !,CASETITLE  			&  			!!CASE TITLE                                 
                                !,HMAX       			&  			!!GLOBAL MAXIMUM DEPTH
                                !,HMIN       			&  			!!GLOBAL MINIMUM DEPTH
                                !,UMOL       			&  			!!VERTICAL DIFFUSION COEFFICIENT
                                !,HORCON     			&  			!!HORIZONTAL DIFFUSION COEFFICIENT
                                !,DTI        			&  			!!internal time step
                                !,HORZMIX    			&   		!!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
                                !,FILENUMBER			&			!!
                                !,PREFF					&			!!
                                !,INPDIR				&			!!
                                !,GEOAREA				&			!!
                                !,RIV_FILENUMBER		&			!!
                    !,INFLOW_TYPE   		&			!!SPECIFIED RIVER INFLOW TYPE (edge/node) 
                    !,POINT_ST_TYPE 		&			!!(calculated/specified)
                    !,PNT_SOURCE    		&			!!point_source
                    !,DAY					&
                                !,in_jday		
                        
                  USE MOD_BCMAP, ONLY : 	&
                                IOBCN,				&!
                                I_OBC_GL
                        
           
                        USE MOD_PAR, ONLY:   &	  
                                !EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
                                !ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
                                NLID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
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

                  
                  IMPLICIT NONE
                  include "mpif.h"	  
                  
            !passed variables
              INTEGER,INTENT(IN) :: MPI_SRC  ! ID of the source process for broadcasting
              REAL(SP), DIMENSION(0:MTLOC),INTENT(OUT)      :: ELL,DTFAL
              REAL(SP), DIMENSION(0:NOBTY+1),INTENT(OUT)    :: UARD_OBCNL
              REAL(SP), DIMENSION(0:NOBTY,KBM1),INTENT(OUT) :: XFLUX_OBCL

            !local variables
              INTEGER ::  MPI_COUNT_TMP=0                !count to send
              INTEGER :: I,J,K
                  INTEGER :: IERR
            !Broadcast global variables
             IF(PAR)THEN  
                      MPI_COUNT_TMP=NOBTY+1   ; CALL MPI_BCAST(  UARD_OBCN_GL(:),MPI_COUNT_TMP,MPI_F,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=NOBTY*KBM1; CALL MPI_BCAST(XFLUX_OBC_GL(:,:),MPI_COUNT_TMP,MPI_F,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=MGL       ; CALL MPI_BCAST(      DTFAL_GL(:),MPI_COUNT_TMP,MPI_F,MPI_SRC,MPI_COMM_WORLD,IERR)
                      MPI_COUNT_TMP=MGL       ; CALL MPI_BCAST(        ELL_GL(:),MPI_COUNT_TMP,MPI_F,MPI_SRC,MPI_COMM_WORLD,IERR)
             ENDIF

            !---then push the global arrays into local arays including halo nodes and halo elements
             IF(PAR)THEN
                DO I=1,MLOC
                          !WRITE(*,*)'ELL_GL Debug here here 03 '
                  ELL(I) = ELL_GL(NGID(I))
                          !WRITE(*,*)'ELL_GL Debug here here 3 '
                ENDDO

                DO I=1,NHN
                          !WRITE(*,*)'ELL_GL Debug here here 04 '
                  ELL(I+MLOC) = ELL_GL(HN_LST(I))
                          !WRITE(*,*)'ELL_GL Debug here here 4 '
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
                     CALL pstop
                ENDIF

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
                     CALL PSTOP
                ENDIF

                DO I=1,MLOC
                  DTFAL(I) = DTFAL_GL(NGID(I))
                ENDDO

                DO I=1,NHN
                  DTFAL(I+MLOC) = DTFAL_GL(HN_LST(I))
                ENDDO

              ENDIF

            RETURN
                
            END SUBROUTINE BROADCAST_HYDRO_REDUCED
#endif
                
#if defined (MULTIPROCESSOR)
              SUBROUTINE SCATTER(N1,NL,NT,NTG,KT,MYID,NPROCS,GM,CM,NID,NID_X,AG,AL,BG,BL,CG,CL)
              !======================================================================!
              ! PASS global array AG from master proc to child procs local array AL  !
              ! Also can do the same for BG -->BL, CG->CL if present                 !
              !-----------------------------------------------------------------------

               USE MOD_TYPES , ONLY:   & 
                                GMAP, 		&	! grid mapping type
                                !NSIZE, 		&	! size
                                !LOC_2_GL, 	&	! conversion from local to global grid for a mapping
                                        !
                                COMM!, 		&	! communication world type
                                !NSND, 		&	! sending buffer size
                                !NRCV, 		&	! receiving buffer size
                                !RCPT, 		&	! pointer in receiving buffer
                                !SNDP,       &    ! sending buffer
                                !RCVP, 		&	! receiving buffer
                                !MLTP, 		&	! muliplicity of each processor
                                 !
                                !BC!, 		&	!boundary type
                                !NTIMES,		&	!
                                !TIMES,		&	!
                                !LABEL

               USE MOD_PAR, ONLY :  EXCHANGE
               USE MOD_CONTROL, ONLY :  PAR!,SERIAL,MSR 
               USE MOD_PREC, ONLY : MPI_F,CDF_PREC,SP,DP,MPI_DP
                  
               USE MOD_BUFFERS, ONLY : SNDBUF, RCVBUF
                   
               IMPLICIT NONE
                   include "mpif.h"
                   
                   
               INTEGER, INTENT(IN)    :: N1,  &  !start index in local array
                                         NL,  &  !size of local array AL without halos
                                         NT,  &  !size of local array AL with halos 
                                         NTG, &  !size of global array
                                         KT,  &  !size of 2nd dimension of both global array AG and local array AL
                                         MYID,&  !ID of the calling proc
                                         NPROCS  !total number of procs

               TYPE(GMAP), INTENT(IN)   :: GM(NPROCS)                 !mapping (element or node)
               REAL(CDF_PREC),   INTENT(IN)   :: AG(1:NTG,KT)         !global arry AG always starts at 1
               REAL(SP),   INTENT(OUT)  :: AL(N1:NT,KT)               !local array AL, size includs halos, NT-NL = number of halo elemetns/nodes
               REAL(CDF_PREC),   INTENT(IN),  OPTIONAL   :: BG(1:NTG,KT)
               REAL(SP),   INTENT(OUT), OPTIONAL   :: BL(N1:NT,KT)
               REAL(CDF_PREC),   INTENT(IN),  OPTIONAL   :: CG(1:NTG,KT)
               REAL(SP),   INTENT(OUT), OPTIONAL   :: CL(N1:NT,KT)
               TYPE(COMM), INTENT(IN)              :: CM(NPROCS)      !communicator (for elements or nodes) 
               INTEGER,    INTENT(IN)              :: NID(0:NTG)      !local index of a global element/node, 
                                                                      !      0 if out of current proc
               INTEGER,    INTENT(IN)              :: NID_X(0:NTG)    !local index of a global element/node including halos,
                                                                      !      0 if out of current proca
              !local variables 
               INTEGER ::  IERR
               INTEGER ::  MPI_COUNT_TMP=0  !count to send
               INTEGER ::  I,J,K
               INTEGER ::  IP
               INTEGER ::  TSND
               INTEGER ::  TRCV
               INTEGER ::  ISNDTAG
               INTEGER ::  IRCVTAG
               INTEGER ::  ISND
               INTEGER ::  IRCV
        !       REAL(SP), ALLOCATABLE, DIMENSION(:) :: RCVBUF,SNDBUF
               INTEGER ::  NMSG
               INTEGER ::  IREQR(NPROCS),IREQS(NPROCS)
               INTEGER ::  STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS)
               INTEGER ::  INDX
               INTEGER ::  NVAR_PASS = 1         !number of variables to pass

               LOGICAL ::  BYES = .FALSE.,  &    !Yes if BG and BL are passed in
                           CYES = .FALSE.        !Yes if CG and CL are passed in 

#if(0)  /* skip debugging after testing */
               !debugging local variables
                REAL(SP) :: SLEEPTIME
                CHARACTER(LEN=1024)::COMMAND_STR,SLEEPTIME_STR
#endif

                IF(PAR)THEN

                 !master proc knows every child procs' global ID (element or nodes) through mapping GM
                 !GM contains the following attributes:
                 !        NSIZE             --- size of the footprint on a proc
                 !        LOC_2_GL          --- global index of a local arraya
                 !Note that GM is only avaiable on master proc

                 !on the other hand, all procs know the correspoonding global ID through NID and NID_X
                 !NID(I)   gives local value of global element or nodes  (0:NTG) (RETURNS 0 IF I IS NOT AN INTERIOR element or node)
                 !NID_X(I) gives local value of global element or nodes including halos  (0:NTG) (RETURNS 0 IF I IS NOT AN INTERIOR OR HALO element/node)
                 

                  !sanity check on optional argments
                  IF(PRESENT(BG).XOR.PRESENT(BL))THEN
                     WRITE(*,*)'ERROR calling SCATTER(), BG and BL must be passed at the same time'
                     CALL PSTOP
                  ENDIF
                  IF(PRESENT(CG).XOR.PRESENT(CL))THEN
                     WRITE(*,*)'ERROR calling SCATTER(), CG and CL must be passed at the same time'
                     CALL PSTOP
                  ENDIF

                  IF(PRESENT(BG).AND.PRESENT(BL))BYES=.TRUE. 
                  IF(PRESENT(CG).AND.PRESENT(CL))CYES=.TRUE. 

                  IF(CYES.EQ..TRUE..AND.BYES.EQ..FALSE.)THEN
                     WRITE(*,*)'ERROR calling SCATTER(), BG and BL must be avaiable when CG and CL are used'
                     CALL PSTOP 
                  ENDIF

                  IF(BYES) NVAR_PASS=NVAR_PASS+1
                  IF(CYES) NVAR_PASS=NVAR_PASS+1

                  IF(.NOT.ALLOCATED(RCVBUF))ALLOCATE(RCVBUF(1:NTG*KT*NVAR_PASS))  !actually can find maximum length across all procs so do not 
                                                                              !really need size 1:NTG, however child does not have GM
                                                                              !so this buffer size is hard for child to allocate
                                                                              !hence we just use NTG as an overkill but simple
                  IF(.NOT.ALLOCATED(SNDBUF))ALLOCATE(SNDBUF(1:NTG*KT*NVAR_PASS))  

                  !Child receives  data
                   IF(MYID/=1)THEN
                         TRCV=0
                         MPI_COUNT_TMP=NL*NVAR_PASS*KT    ! number of data points to receive (locally)
                         RCVBUF(:)=0.0
                         TRCV=TRCV+1
                         IRCVTAG= 1000+MYID
                        !non-blocking receiving method with receiving status in IREQR to querry later
                         CALL MPI_IRECV(RCVBUF(1),MPI_COUNT_TMP,MPI_F,0,IRCVTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
                    ENDIF


                  !Master sends out data to child processes
                   IF(MYID==1)THEN
                       TSND=0  !count number of sends by master proc
                       DO IP=2,NPROCS  !send data to child procs (2:NPROCS)
                           !Master finds chunk of data for each child
                           !find size of buffer to send
                            SNDBUF(:)=0.0
                           !fill out buffer
                            ISND=0
                            DO I=1,GM(IP)%NSIZE    !loop through all values for proc IP
                               DO K=1,KT
                                  ISND=ISND+1
                                  SNDBUF(ISND)=AG(GM(IP)%LOC_2_GL(I),K)   !master fetches global data for this child
                                                                        !and puts it in the sending buffer
                               ENDDO
                            ENDDO

                            IF(BYES)THEN  !pass array BG
                               DO I=1,GM(IP)%NSIZE    !loop through all values for proc IP
                               DO K=1,KT
                                  ISND=ISND+1
                                  SNDBUF(ISND)=BG(GM(IP)%LOC_2_GL(I),K)   !master fetches global data for this child
                                                                        !and puts it in the sending buffer
                               ENDDO
                               ENDDO
                            ENDIF

                            IF(CYES)THEN  !pass array CG
                               DO I=1,GM(IP)%NSIZE    !loop through all values for proc IP
                               DO K=1,1
                                  ISND=ISND+1
                                  SNDBUF(ISND)=CG(GM(IP)%LOC_2_GL(I),K)   !master fetches global data for this child
                                                                        !and puts it in the sending buffer
                               ENDDO
                               ENDDO
                            ENDIF

                            MPI_COUNT_TMP=ISND         !number of data points to send to proc IP. 
                            !Master sends chunk of data to each child
                            ISNDTAG=IP+1000      !Note possible bug here with 1000, we need to make sure ISNDTAG is unique for each process IP
                            TSND=TSND+1

                            !Send with non-blocking method, with request in IREQS to querry later
                            CALL MPI_ISEND(SNDBUF(1),MPI_COUNT_TMP, MPI_F, IP-1, ISNDTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
                             
                            !or use the blocking send
                            !CALL MPI_SEND(SNDBUF(1),MPI_COUNT_TMP, MPI_F, IP-1, ISNDTAG,MPI_COMM_WORLD,IERR)
                       ENDDO
                   ENDIF 

        !           !Child receives data 
        !           IF(MYID/=1)THEN
        !                 TRCV=0
        !                 MPI_COUNT_TMP=NL*NVAR_PASS*KT    ! number of data points to receive (locally)
        !                 RCVBUF(:)=0.0
        !                 TRCV=TRCV+1
        !                 IRCVTAG= 1000+MYID
        !                !non-blocking receiving method with receiving status in IREQR to querry later
        !                 CALL MPI_IRECV(RCVBUF(1),MPI_COUNT_TMP,MPI_F,0,IRCVTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
        !
        !                !or blocking receve
        !                ! CALL MPI_RECV(RCVBUF(1),MPI_COUNT_TMP,MPI_F,0,IRCVTAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,IERR)
        !           ENDIF

                  !Loop over procs until a message is received and unpack the data
                   IF(MYID/=1)THEN
                      DO NMSG=1,TRCV
                         CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR,RCVBUF,SNDBUF) !use this when use the non-blocking recieve method
                         !Unpack the message
                         IRCV=0                    !index in the buffer
                         DO I=1,NL                 !NL must be equal to GM(MYID)%NSIZE although GM is only available in master
                            DO K=1,KT
                               IRCV=IRCV+1
                               AL(I,K)=RCVBUF(IRCV)      !Child retrieves the data from the receiving buffer
                            ENDDO
                         ENDDO

                         IF(BYES)THEN   !retrive array BL
                            DO I=1,NL                 !NL must be equal to GM(MYID)%NSIZE although GM is onlya vailable in master
                            DO K=1,KT
                               IRCV=IRCV+1
                               BL(I,K)=RCVBUF(IRCV)      !Child retrieves the data from the receiving buffer
                            ENDDO
                            ENDDO
                         ENDIF

                         IF(CYES)THEN   !retrieve array CL
                            DO I=1,NL                 !NL must be equal to GM(MYID)%NSIZE although GM is onlya vailable in master
                            DO K=1,KT
                               IRCV=IRCV+1
                               CL(I,K)=RCVBUF(IRCV)      !Child retrieves the data from the receiving buffer
                            ENDDO
                            ENDDO
                         ENDIF

                      ENDDO
                   ENDIF

                   !WAIT FOR COMPLETION OF NON-BLOCKING SENDS (turn on this when using non-blocking sends)
                   IF(MYID==1)THEN
                      CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR,RCVBUF,SNDBUF) ! Add RBUF,SBUF to fool compiler into NOT optimizing into a bug
                   ENDIF

                   IF(ALLOCATED(RCVBUF))DEALLOCATE(RCVBUF)
                   IF(ALLOCATED(SNDBUF))DEALLOCATE(SNDBUF)

                   !Master gets its own copy of the global data into local array
                   IF(MYID==1)THEN
                      DO K=1,KT
                         DO I=1,GM(MYID)%NSIZE
                            AL(I,K) = AG(GM(MYID)%LOC_2_GL(I),K)
                            IF(BYES)BL(I,K) = BG(GM(MYID)%LOC_2_GL(I),K)
                            IF(CYES)CL(I,K) = CG(GM(MYID)%LOC_2_GL(I),K)
                         ENDDO
                      ENDDO
                   ENDIF

                   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

                   !Note here CM is type COMM for elements and nodes, i.e. element communicator and node communicator
                   !They have the followinng properties:

                       !SND: TRUE if you are to send to processor
                       !RCV: TRUE IF you are to receive from processor
                       !NSND: Number of data to send to processor
                       !NRCV: Number of data to receive from processor
                       !SNDP: Array pointing to locations to send to processor
                       !RCVP: Array pointing to locations received from process
                       !RCPT: Pointer to location in receive buffer

                   !The EXCHANGE() subroutine will prepare these properties and then send or receive
                   !and then set the halo elements/nodes with values in the received buffer
                 
                   !
                   !Call the EXCHANGE() subroutine after scattering to make sure
                   !each process gets its halo node/element information from neighbors
                   !

                   CALL EXCHANGE(CM, NT, KT, MYID, NPROCS, AL)

                   IF(BYES)CALL EXCHANGE(CM, NT, KT, MYID, NPROCS, BL)
                   IF(CYES)CALL EXCHANGE(CM, NT, KT, MYID, NPROCS, CL)


#if(0) /*skip debugging and testing */ 
                   !
                   !Test : master gather's the same data that was scattered to child processes
                   !       and check if the gathered data is the same as what was sent out

                   !Master prints the global data
                   DO IP=1,NPROCS
                      IF(IP.eq.1.AND.MYID.eq.1)THEN
                            DO K=1,KT
                            DO I=1,NTG
                                WRITE(*,        '(A9,1x,I10,1x,A2,1x,I10,1x,A8,F20.5)')'M-A,GIND=',I,'K=',K,'AG(I,K)=',AG(I,K)
                                IF(BYES)WRITE(*,'(A9,1x,I10,1x,A2,1x,I10,1x,A8,F20.5)')'M-B,GIND=',I,'K=',K,'BG(I,K)=',BG(I,K)
                                IF(CYES)WRITE(*,'(A9,1x,I10,1x,A2,1x,I10,1x,A8,F20.5)')'M-C,GIND=',I,'K=',K,'CG(I,K)=',CG(I,K)
                            ENDDO
                            ENDDO
                      ENDIF
                      SLEEPTIME=MOD(MYID,2)                       !(MYID+1)*1  ! sleep time in microseconds
                      SLEEPTIME_STR=''
                      WRITE(SLEEPTIME_STR,'(I0)')INT(SLEEPTIME)    ! put number of microseconds in integer in a string
                      COMMAND_STR='usleep'//' '//SLEEPTIME_STR     ! command line string for system command
                      CALL SYSTEM(TRIM(COMMAND_STR))
                      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
                   ENDDO

                   !every one prints local data 
                   DO IP=1,NPROCS                 !printing is serialized by nesting MPI_BARRIER in this IP loop
                       IF(IP.eq.MYID)THEN
                         DO K=1,KT
                         DO I=1,NTG               !loop globally
                            IF(NID(I).NE.0)THEN
                               WRITE(*,        '(A9,1x,I10,1x,A2,1x,I10,1x,A8,1x,F20.5)')'C-A,GIND=',I,'K=',K,'AL(I,K)=',AL(NID(I),K),MYID
                               IF(BYES)WRITE(*,'(A9,1x,I10,1x,A2,1x,I10,1x,A8,1x,F20.5)')'C-B,GIND=',I,'K=',K,'BL(I,K)=',BL(NID(I),K),MYID
                               IF(CYES)WRITE(*,'(A9,1x,I10,1x,A2,1x,I10,1x,A8,1x,F20.5)')'C-C,GIND=',I,'K=',K,'CL(I,K)=',CL(NID(I),K),MYID
                            ENDIF
                            IF(NID_X(I).NE.0.AND.NID(I).EQ.0)THEN  !print the HALO data only, this is to make sure HALO data are correct
                               WRITE(*,        '(A10,1x,I10,1x,A2,1x,I10,1x,A8,1x,F20.5)')'CH-A,GIND=',I,'K=',K,'AL(I,K)=',AL(NID_X(I),K),MYID
                               IF(BYES)WRITE(*,'(A10,1x,I10,1x,A2,1x,I10,1x,A8,1x,F20.5)')'CH-B,GIND=',I,'K=',K,'BL(I,K)=',BL(NID_X(I),K),MYID
                               IF(CYES)WRITE(*,'(A10,1x,I10,1x,A2,1x,I10,1x,A8,1x,F20.5)')'CH-C,GIND=',I,'K=',K,'CL(I,K)=',CL(NID_X(I),K),MYID
                            ENDIF
                         ENDDO
                         ENDDO
                      ENDIF

                      SLEEPTIME=MOD(MYID,2)    !Sleep shortly so that different procs spend different amount time in IP loop                   
                      SLEEPTIME_STR=''         !to make sure this printing is serialized
                      WRITE(SLEEPTIME_STR,'(I0)')INT(SLEEPTIME)    ! put number of microseconds in integer in a string
                      COMMAND_STR='usleep'//' '//SLEEPTIME_STR     ! command line string for system command
                      CALL SYSTEM(TRIM(COMMAND_STR))
                      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
                   ENDDO
                   CALL PSTOP
#endif
              ENDIF
              RETURN
            END SUBROUTINE SCATTER
#endif

END MODULE MOD_HYDRO
