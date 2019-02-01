
        !==============================================================================|
        !  VARIABLES                                                                   |
        !==============================================================================|

        !WL MODULE ALL_VARS  !Wen Long: this is very confusing. Modify to MOD_HYDROVARS

        MODULE MOD_HYDROVARS  !Wen Long: this is very confusing. Modify to MOD_HYDROVARS
		
          USE MOD_PREC, ONLY: SP
          USE MOD_LIMS, ONLY: NSTATIONMAX, KB,KBM1,NTLOC, MTLOC
		  USE MOD_TYPES , ONLY:  & 
			!GMAP, 		&	! grid mapping type
			!NSIZE, 		&	! size
			!LOC_2_GL, 	&	! conversion from local to global grid for a mapping
		                !
			!COMM, 		&	! communication world type
			!NSND, 		&	! sending buffer size
			!NRCV, 		&	! receiving buffer size
			!RCPT, 		&	! pointer in receiving buffer
			!SNDP,       &    ! sending buffer
			!RCVP, 		&	! receiving buffer
			!MLTP, 		&	! muliplicity of each processor
                         !
			 BC!, 		&	!boundary type
			!NTIMES,		&	!
			!TIMES,		&	!
			!LABEL

		  USE MOD_SIZES, ONLY : MGL, NGL, NOBTY
			 
          IMPLICIT NONE
          SAVE

         !--Constants-------------------------------------------------------------------!
          REAL(SP), PARAMETER :: GRAV      = 9.81
          REAL(SP), PARAMETER :: PI        = 3.141592653
          REAL(SP), PARAMETER :: PI2       = 6.283185307
          REAL(SP), PARAMETER :: ZERO      = 0.0 
          REAL(SP), PARAMETER :: ONE_THIRD = 1.0/3.0 

        !--------------------------Temporary Array------------------------------------------!

          INTEGER, ALLOCATABLE :: NVG(:,:)

        !--------------------------Global Grid Variables------------------------------------!

          REAL(SP), ALLOCATABLE :: XG(:)               !!GLOBAL X-COORD AT NODE 
          REAL(SP), ALLOCATABLE :: YG(:)               !!GLOBAL X-COORD AT NODE 
          REAL(SP), ALLOCATABLE :: HG(:)               !!GLOBAL DEPTH AT NODE 
          REAL(SP), ALLOCATABLE :: XCG(:)              !!GLOBAL X-COORD AT FACE CENTER 
          REAL(SP), ALLOCATABLE :: YCG(:)              !!GLOBAL X-COORD AT FACE CENTER 

          !--------------------------Grid Metrics---------------------------------------------!

          REAL(SP)              :: VXMIN,VYMIN,VXMAX,VYMAX
          REAL(SP), ALLOCATABLE :: XC(:)               !!X-COORD AT FACE CENTER 
          REAL(SP), ALLOCATABLE :: YC(:)               !!Y-COORD AT FACE CENTER
          REAL(SP), ALLOCATABLE :: VX(:)               !!X-COORD AT GRID POINT
          REAL(SP), ALLOCATABLE :: VY(:)               !!Y-COORD AT GRID POINT
          REAL(SP), ALLOCATABLE :: ART(:)              !!AREA OF ELEMENT
          REAL(SP), ALLOCATABLE :: ART1(:),ART1_GL(:)             !!AREA OF NODE-BASE CONTROl VOLUME, and global
          REAL(SP), ALLOCATABLE :: ART2(:)             !!AREA OF ELEMENTS AROUND NODE

          !----------------Node, Boundary Condition, and Control Volume-----------------------!
		   
		   !The following are moved to tge.F
		   
		   !INTEGER, ALLOCATABLE :: NV(:,:)             !!NODE NUMBERING FOR ELEMENTS
		  
           !INTEGER, ALLOCATABLE :: ISBCE(:)     
		   !INTEGER, ALLOCATABLE :: ISONB(:)  !!NODE MARKER = 0,1,2   
		   !INTEGER, ALLOCATABLE :: NBE(:,:)  !!INDICES OF ELMNT NEIGHBORS
           !INTEGER, ALLOCATABLE :: NTVE(:)         		
           !INTEGER, ALLOCATABLE :: NTSN(:)      

          ! INTEGER, ALLOCATABLE :: IEC(:,:)     	!1
          ! INTEGER, ALLOCATABLE :: IENODE(:,:)  	!2
		  ! REAL(SP),ALLOCATABLE :: XIJC(:) 			!3
          ! REAL(SP),ALLOCATABLE :: YIJC(:)				!4
		  ! REAL(SP),ALLOCATABLE :: DLTXC(:)			!5
          ! REAL(SP),ALLOCATABLE :: DLTYC(:)			!6
          ! REAL(SP),ALLOCATABLE :: DLTXYC(:)			!7

 		  !REAL(SP),ALLOCATABLE :: SITAC(:) 			!8
		  !INTEGER, ALLOCATABLE :: ISBC(:) 			!9    
          !INTEGER, ALLOCATABLE :: NBVE(:,:)		!10
          !INTEGER, ALLOCATABLE :: NBVT(:,:)		!11
		  !INTEGER, ALLOCATABLE :: NBSN(:,:)		!12
          !INTEGER, ALLOCATABLE :: NIEC(:,:)		!13
          !INTEGER, ALLOCATABLE :: NTRG(:)			!14
		  
		  !REAL(SP),ALLOCATABLE :: XIJE(:,:) 			!15
          !REAL(SP),ALLOCATABLE :: YIJE(:,:) 			!16
          !REAL(SP),ALLOCATABLE :: DLTXE(:)				!17
          !REAL(SP),ALLOCATABLE :: DLTYE(:)				!18
          !REAL(SP),ALLOCATABLE :: DLTXYE(:)			!19 
          !REAL(SP),ALLOCATABLE :: SITAE(:) 			!20

          !INTEGER, ALLOCATABLE :: LISBCE_1(:)  	!21   !!LIST OF ELEMENTS WITH ISBCE=1
          !INTEGER, ALLOCATABLE :: LISBCE_2(:)  	!22   !!LIST OF ELEMENTS WITH ISBCE=2
          !INTEGER, ALLOCATABLE :: LISBCE_3(:)  	!23   !!LIST OF ELEMENTS WITH ISBCE=3
           
          !REAL(SP),ALLOCATABLE :: EPOR(:)      		!24   !!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)


!
!WLong commented out these, they are never used
!          INTEGER, ALLOCATABLE :: IBCGEO(:)        !!LOCAL GEOSTROPHIC FRICTION CORRECTION NODES  !Never used
		 
!WLong moved this to mod_bcs.F
!          INTEGER, ALLOCATABLE :: N_ICELLQ(:,:)    !!node number of the edge for element edge source, 
							
!
!		  !------------shape coefficient arrays and control volume metrics--------------------!
!
!          REAL(SP), ALLOCATABLE :: A1U(:,:)      
!          REAL(SP), ALLOCATABLE :: A2U(:,:)     
!          REAL(SP), ALLOCATABLE :: AWX(:,:)   
!          REAL(SP), ALLOCATABLE :: AWY(:,:)  
!          REAL(SP), ALLOCATABLE :: AW0(:,:) 
!
         !----------------1-d arrays for the sigma coordinate -------------------------------!

          REAL(SP), ALLOCATABLE :: Z(:)                    !!SIGMA COORDINATE VALUE 
          REAL(SP), ALLOCATABLE :: ZZ(:)                   !!INTRA LEVEL SIGMA VALUE
          REAL(SP), ALLOCATABLE :: DZ(:)                   !!DELTA-SIGMA VALUE
          REAL(SP), ALLOCATABLE :: DZZ(:)                  !!DELTA OF INTRA LEVEL SIGMA 
		  
		  !2d arrays from FVCOM input for spatial data, used for netcdf output
		  !added by B Clark Dec. 2018
		  REAL(SP), ALLOCATABLE :: SIGLEV(:,:)                    !!SIGMA COORDINATE VALUE 
          REAL(SP), ALLOCATABLE :: SIGLAY(:,:)                   !!Sigma Layer value 

        !---------------2-d flow variable arrays at elements-------------------------------!

          REAL(SP), ALLOCATABLE :: H1(:)            !!BATHYMETRIC DEPTH   

        !---------------2-d flow variable arrays at nodes----------------------------------!

          REAL(SP), ALLOCATABLE :: H(:)            !!BATHYMETRIC DEPTH   
          REAL(SP), ALLOCATABLE :: D(:)            !!CURRENT DEPTH   
          REAL(SP), ALLOCATABLE :: DT(:)            !!DEPTH AT PREVIOUS TIME STEP
          REAL(SP), ALLOCATABLE :: DT1(:)           !!DEPTH AT PREVIOUS TIME STEP
          REAL(SP), ALLOCATABLE :: EL(:)           !!CURRENT SURFACE ELEVATION
          REAL(SP), ALLOCATABLE :: ET(:)           !!SURFACE ELEVATION AT PREVIOUS TIME STEP
          REAL(SP), ALLOCATABLE :: DTFA(:)          !!ADJUSTED DEPTH FOR MASS CONSERVATION

        !---------------- internal mode   arrays-(element based)----------------------------!

          REAL(SP), ALLOCATABLE :: UU(:,:)         !X-VELOCITY
          REAL(SP), ALLOCATABLE :: VV(:,:)         !Y-VELOCITY
          REAL(SP), ALLOCATABLE :: UUT(:,:)        !X-VELOCITY FROM PREVIOUS TIMESTEP
          REAL(SP), ALLOCATABLE :: VVT(:,:)        !Y-VELOCITY FROM PREVIOUS TIMESTEP
!         REAL(SP), ALLOCATABLE :: WWT(:,:)        !Z-VELOCITY FROM PREVIOUS TIMESTEP !WLong never used
          REAL(SP), ALLOCATABLE :: WTST(:,:)       !Vertical velocity in sigma from PREVIOUS TIMESTEP
          REAL(SP), ALLOCATABLE :: UARD_OBCNT(:)   !tykim
          REAL(SP), ALLOCATABLE :: XFLUX_OBCT(:,:) !tykim
          REAL(SP), ALLOCATABLE :: DTFAT(:)        !tykim
!         REAL(SP), ALLOCATABLE :: TT_T(:,:)       !tykim !WLong, never used
!         REAL(SP), ALLOCATABLE :: SALTT(:,:)      !tykim !WLong, never used

          !-----------------------3d variable arrays-(node based)-----------------------------!

          REAL(SP), ALLOCATABLE :: WTS(:,:)        !!VERTICAL VELOCITY IN SIGMA SYSTEM
          REAL(SP), ALLOCATABLE :: UARD_OBCN(:)    !! tykim  
          REAL(SP), ALLOCATABLE :: XFLUX_OBC(:,:)  !! tykim     
!         REAL(SP), ALLOCATABLE :: WTTS(:,:)       !!VERTICAL VELOCITY IN SIGMA SYSTEM     !WLong, never used
          REAL(SP), ALLOCATABLE :: KH(:,:)         !!TURBULENT DIFFUSIVITY
       
		  !--------------------hydrodynamics-------------------------------------------------

          REAL(SP), ALLOCATABLE :: VISCOFH(:,:)

          REAL(SP), ALLOCATABLE :: UNC1(:,:)
          REAL(SP), ALLOCATABLE :: VNC1(:,:)
          REAL(SP), ALLOCATABLE :: WNC1(:,:)
          REAL(SP), ALLOCATABLE :: WTSNC1(:,:)
          REAL(SP), ALLOCATABLE :: UARD_OBCNNC1(:)
          REAL(SP), ALLOCATABLE :: XFLUX_OBCNC1(:,:)
          REAL(SP), ALLOCATABLE :: DTFANC1(:)
          REAL(SP), ALLOCATABLE :: KHNC1(:,:)
          REAL(SP), ALLOCATABLE :: TNC1(:,:)
          REAL(SP), ALLOCATABLE :: SNC1(:,:)
          REAL(SP), ALLOCATABLE :: ELNC1(:)
		  
          REAL(SP), ALLOCATABLE :: UNC2(:,:)
          REAL(SP), ALLOCATABLE :: VNC2(:,:)
          REAL(SP), ALLOCATABLE :: WNC2(:,:)
          REAL(SP), ALLOCATABLE :: WTSNC2(:,:)
          REAL(SP), ALLOCATABLE :: UARD_OBCNNC2(:)
          REAL(SP), ALLOCATABLE :: XFLUX_OBCNC2(:,:)
          REAL(SP), ALLOCATABLE :: DTFANC2(:)
          REAL(SP), ALLOCATABLE :: KHNC2(:,:)
          REAL(SP), ALLOCATABLE :: TNC2(:,:)
          REAL(SP), ALLOCATABLE :: SNC2(:,:)
          REAL(SP), ALLOCATABLE :: ELNC2(:)
          
          INTEGER(4) ::  num_hyd_ints   !number of records in each hydrodynamics netcdf file

        !WL  INTEGER :: debug3602 = 3602
          
          TYPE(BC) :: TIME_MAP
          REAL(SP) :: THOUR1     !!SIMULATION TIME AT END OF CURRENT EXTERNAL STEP (IEXT) IN HOURS
          REAL(SP) :: THOUR

          !Wen Long added following for hydrodynamics input
           CHARACTER(LEN=1024) :: NCFILE_DIR,NCFILE_PREFIX,NCFILE_SUFFIX, NCFILE_NUMBER
           CHARACTER(LEN=1024) :: FORMAT_STR
           CHARACTER(LEN=1024) :: hydro_dir,              & ! directory name where hydrodynamics results (netcdf) files are stored
                                  hydro_prefix,           & ! prefix of file name, e.g. 'psm_'
                                  hydro_suffix              ! suffix of filename, e.g. '.nc'
           INTEGER(4)          :: hydro_filenumwidth,     & ! number of digits in filename following hydro_prefix, e.g. 4 for psm_0002.nc
                                  hydro_filenumstart,     & ! starting number of the file name in the digital part of the file name, e.g. 185 for psm_0185.nc
                                  hydro_Nrec                ! number of records in each of hydrodynamics file
           REAL(SP) :: hydro_dlt                                ! time step in hydrodynamics file (in seconds), e.g. 100 for 100sec


           REAL(SP) :: t_his_start,t_his_end,t_his_dlt  !starting time, ending time, and interval of history outputs (days)

           INTEGER(4) :: Nstation, NstationNum_GL(NstationMax) !maximum number of station is NstationMax!
           REAL(SP)   :: t_stn_start,t_stn_end,t_stn_dlt       !starting time, ending time, and interval of station outputs (days)

           CHARACTER(LEN=1024) :: STNFN                        !file name for station output
           CHARACTER(LEN=1024) :: HISFN                        !file name for history output
           CHARACTER(LEN=1024) :: HISFN_PREFIX                 !prefix of history output file
           CHARACTER(LEN=1024) :: HISFN_EXT                    !extention name of history output file
           CHARACTER(LEN=1024) :: HISFN_FINAL                   
           LOGICAL :: HISFN_SPLIT_BYLEVEL =.FALSE.             !Ture or False for splitting history output in files level by level (default is .FALSE.)

           NAMELIST /hydro_netcdf/hydro_dir,hydro_prefix,hydro_suffix, &
                           hydro_filenumwidth,hydro_filenumstart, &
                           hydro_Nrec,hydro_dlt
           NAMELIST /wqm_history/HISFN,t_his_start,t_his_end,t_his_dlt,HISFN_SPLIT_BYLEVEL
           NAMELIST /wqm_stations/STNFN,Nstation,NstationNum_GL,t_stn_start,t_stn_end,t_stn_dlt


           INTEGER(4) :: IFNC     !file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise 
                                  !for restart
           INTEGER(4) :: NTRECNC=0  !time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file. 
           INTEGER(4) :: NTHYDRO  !overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read

		   
	   CONTAINS
	   
			!Subroutines:
			!	SUBROUTINE HYDRO_GEOM_ALLOC()
			!	SUBROUTINE HYDRO_GEOM_DEALLOC()
			! 	SUBROUTINE HYDRO_ALLOC()
			! 	SUBROUTINE HYDRO_DEALLOC()
	   
			!Subroutine to allocate geometry related variables from hydrodyanmics (netcdf)
			SUBROUTINE HYDRO_GEOM_ALLOC  
					ALLOCATE(XG(0:MGL)) 		; XG = 0.0_SP 
					ALLOCATE(YG(0:MGL))			; YG = 0.0_SP
					ALLOCATE(NVG(0:NGL,4))		; NVG = 0	  
					ALLOCATE(HG(0:MGL))  		; HG = 0.0_SP
					ALLOCATE(Z(KB))  			; Z   = 0.0_SP    !!SIGMA COORDINATE VALUE 
					ALLOCATE(ZZ(KB)) 			; ZZ  = 0.0_SP    !!INTRA LEVEL SIGMA VALUE
					ALLOCATE(DZ(KB)) 			; DZ  = 0.0_SP    !!DELTA-SIGMA VALUE
					ALLOCATE(DZZ(KB))			; DZZ = 0.0_SP    !!DELTA OF INTRA LEVEL SIGMA 
					
				    ALLOCATE(XCG(0:NGL)) ; XCG = 0.0_SP 
					ALLOCATE(YCG(0:NGL)) ; YCG = 0.0_SP
					
					ALLOCATE(SIGLEV(0:MTLOC,KB)); SIGLEV=0.0
					ALLOCATE(SIGLAY(0:MTLOC,KB)); SIGLAY=0.
					                                   
			END SUBROUTINE HYDRO_GEOM_ALLOC
			
			!Subroutine to allocate geometry related variables from hydrodyanmics (netcdf)
			SUBROUTINE HYDRO_GEOM_DEALLOC
			
					IF(ALLOCATED(XG))DEALLOCATE(XG)
					IF(ALLOCATED(YG))DEALLOCATE(YG)
					IF(ALLOCATED(NVG))DEALLOCATE(NVG)
					IF(ALLOCATED(HG))DEALLOCATE(HG)
					IF(ALLOCATED(Z))DEALLOCATE(Z)
					IF(ALLOCATED(ZZ))DEALLOCATE(ZZ)
					IF(ALLOCATED(DZ))DEALLOCATE(DZ)
					IF(ALLOCATED(DZZ))DEALLOCATE(DZZ)
					IF(ALLOCATED(XCG))DEALLOCATE(XCG)
					IF(ALLOCATED(YCG))DEALLOCATE(YCG)

			END SUBROUTINE HYDRO_GEOM_DEALLOC
			
			SUBROUTINE HYDRO_ALLOC
			
			    !local geometry
				ALLOCATE(VX(0:MTLOC));      VX  = 0.0_SP  !!X-COORD AT GRID POINT   
				ALLOCATE(VY(0:MTLOC));      VY  = 0.0_SP  !!X-COORD AT GRID POINT   
				ALLOCATE(XC(0:NTLOC));		XC  = 0.0_SP  !!X-COORD AT FACE CENTER
				ALLOCATE(YC(0:NTLOC));		YC  = 0.0_SP  !!X-COORD AT FACE CENTER
				ALLOCATE(H(0:MTLOC));  		H 	= 0.0_SP  !!BATHYMETRIC DEPTH
				ALLOCATE(H1(0:NTLOC));		H1  = 0.0_SP  !!BATHYMETRIC DEPTH

				!WLong moved from cell_area.F to here
				ALLOCATE(ART(0:NTLOC))      ;ART  = 0.0_SP !!AREA OF ELEMENT
				ALLOCATE(ART1(0:MTLOC))     ;ART1 = 0.0_SP !!AREA OF NODE-BASE CONTROl VOLUME
				ALLOCATE(ART1_GL(MGL))   ; ART1_GL = 0.0
				
				ALLOCATE(ART2(0:MTLOC))     ;ART2 = 0.0_SP !!AREA OF ELEMENTSAROUND NODE
    				
				ALLOCATE(UNC1(0:NTLOC,KB));  UNC1  = 0.0
				ALLOCATE(VNC1(0:NTLOC,KB));  VNC1  = 0.0
				ALLOCATE(WNC1(0:MTLOC,KB));  WNC1  = 0.0
				ALLOCATE(WTSNC1(0:MTLOC,KB));  WTSNC1  = 0.0
				ALLOCATE(UARD_OBCNNC1(0:NOBTY+1))   ;  UARD_OBCNNC1  = 0.0
				ALLOCATE(XFLUX_OBCNC1(0:NOBTY,KBM1));  XFLUX_OBCNC1  = 0.0
				ALLOCATE(DTFANC1(0:MTLOC))   ;  DTFANC1  = 0.0
				ALLOCATE(KHNC1(0:MTLOC,KB)); KHNC1 = 0.0
				ALLOCATE(TNC1(0:MTLOC,KBM1));  TNC1 = 0.0
				ALLOCATE(SNC1(0:MTLOC,KBM1));  SNC1 = 0.0
				ALLOCATE(ELNC1(0:MTLOC)); ELNC1 = 0.0
 
				ALLOCATE(UNC2(0:NTLOC,KB));  UNC2 = 0.0
				ALLOCATE(VNC2(0:NTLOC,KB));  VNC2 = 0.0
				ALLOCATE(WNC2(0:MTLOC,KB));  WNC2  = 0.0
				ALLOCATE(WTSNC2(0:MTLOC,KB));  WTSNC2  = 0.0
				ALLOCATE(UARD_OBCNNC2(0:NOBTY+1))   ;  UARD_OBCNNC2  = 0.0
				ALLOCATE(XFLUX_OBCNC2(0:NOBTY,KBM1));  XFLUX_OBCNC2  = 0.0
				ALLOCATE(DTFANC2(0:MTLOC))   ;  DTFANC2  = 0.0
				ALLOCATE(KHNC2(0:MTLOC,KB)); KHNC2 = 0.0
				ALLOCATE(TNC2(0:MTLOC,KBM1));  TNC2 = 0.0
				ALLOCATE(SNC2(0:MTLOC,KBM1));  SNC2 = 0.0
				ALLOCATE(ELNC2(0:MTLOC)); ELNC2 = 0.0

				ALLOCATE(UU(0:NTLOC,KB));  UU  = 0.0
				ALLOCATE(VV(0:NTLOC,KB));  VV  = 0.0
				ALLOCATE(UUT(0:NTLOC,KB));  UUT  = 0.0
				ALLOCATE(VVT(0:NTLOC,KB));  VVT  = 0.0				
				ALLOCATE(WTST(0:MTLOC,KB));  WTST  = 0.0
				ALLOCATE(UARD_OBCNT(0:NOBTY+1))   ;  UARD_OBCNT  = 0.0
				ALLOCATE(XFLUX_OBCT(0:NOBTY,KBM1));  XFLUX_OBCT  = 0.0
				ALLOCATE(DTFAT(0:MTLOC)); DTFAT = 0.0
		  
				ALLOCATE(WTS(0:MTLOC,KB));  WTS  = 0.0
				ALLOCATE(UARD_OBCN(0:NOBTY+1))   ;  UARD_OBCN  = 0.0
				ALLOCATE(XFLUX_OBC(0:NOBTY,KBM1));  XFLUX_OBC  = 0.0
				ALLOCATE(KH(0:MTLOC,KB)); KH = 0.0

				ALLOCATE(EL(0:MTLOC)); EL = 0.0
				ALLOCATE(ET(0:MTLOC)); ET = 0.0
				ALLOCATE(D(0:MTLOC)); D = 0.0
				ALLOCATE(DT(0:MTLOC)); DT = 0.0
				ALLOCATE(DT1(0:NTLOC)); DT1 = 0.0
				ALLOCATE(DTFA(0:MTLOC)); DTFA = 0.0
	
				ALLOCATE(VISCOFH(0:NTLOC,KB));  VISCOFH  = 0.0
				
			END SUBROUTINE HYDRO_ALLOC
			
			SUBROUTINE HYDRO_DEALLOC

				!local geometry
				IF(ALLOCATED(VX))DEALLOCATE(VX)			!1
				IF(ALLOCATED(VY))DEALLOCATE(VY)			!2
				IF(ALLOCATED(XC))DEALLOCATE(XC)			!3
				IF(ALLOCATED(YC))DEALLOCATE(YC)			!4
				IF(ALLOCATED(H))DEALLOCATE(H)			!5
				IF(ALLOCATED(H1))DEALLOCATE(H1)			!6

				IF(ALLOCATED(ART))DEALLOCATE(ART)		!7
				IF(ALLOCATED(ART1))DEALLOCATE(ART1)		!8
				IF(ALLOCATED(ART1_GL))DEALLOCATE(ART1_GL)		!8
				IF(ALLOCATED(ART2))DEALLOCATE(ART2)		!9
								
				IF(ALLOCATED(UNC1))DEALLOCATE(UNC1)		!10
				IF(ALLOCATED(VNC1))DEALLOCATE(VNC1)		!11
				IF(ALLOCATED(WNC1))DEALLOCATE(WNC1)		!12
				IF(ALLOCATED(WTSNC1))DEALLOCATE(WTSNC1)	!13
				IF(ALLOCATED(UARD_OBCNNC1))DEALLOCATE(UARD_OBCNNC1)	!14
				IF(ALLOCATED(XFLUX_OBCNC1))DEALLOCATE(XFLUX_OBCNC1)	!15
				IF(ALLOCATED(DTFANC1))DEALLOCATE(DTFANC1)			!16
				IF(ALLOCATED(KHNC1))DEALLOCATE(KHNC1)				!17
				IF(ALLOCATED(TNC1))DEALLOCATE(TNC1)					!18
				IF(ALLOCATED(SNC1))DEALLOCATE(SNC1)					!19
				IF(ALLOCATED(ELNC1))DEALLOCATE(ELNC1)				!20
 
				IF(ALLOCATED(UNC2))DEALLOCATE(UNC2)					!21
				IF(ALLOCATED(VNC2))DEALLOCATE(VNC2)					!22
				IF(ALLOCATED(WNC2))DEALLOCATE(WNC2)					!23
				IF(ALLOCATED(WTSNC2))DEALLOCATE(WTSNC2)				!24
				IF(ALLOCATED(UARD_OBCNNC2))DEALLOCATE(UARD_OBCNNC2)	!25
				IF(ALLOCATED(XFLUX_OBCNC2))DEALLOCATE(XFLUX_OBCNC2)	!26
				IF(ALLOCATED(DTFANC2))DEALLOCATE(DTFANC2)			!27
				IF(ALLOCATED(KHNC2))DEALLOCATE(KHNC2)				!28
				IF(ALLOCATED(TNC2))DEALLOCATE(TNC2)					!29
				IF(ALLOCATED(SNC2))DEALLOCATE(SNC2)					!30
				IF(ALLOCATED(ELNC2))DEALLOCATE(ELNC2)				!31

				IF(ALLOCATED(UU))DEALLOCATE(UU)						!32
				IF(ALLOCATED(VV))DEALLOCATE(VV)						!33
				
				IF(ALLOCATED(UUT))DEALLOCATE(UUT)					!34	
				IF(ALLOCATED(VVT))DEALLOCATE(VVT)					!35
				IF(ALLOCATED(WTST))DEALLOCATE(WTST)					!36
				IF(ALLOCATED(UARD_OBCNT))DEALLOCATE(UARD_OBCNT)		!37
				IF(ALLOCATED(XFLUX_OBCT))DEALLOCATE(XFLUX_OBCT)		!38
				IF(ALLOCATED(DTFAT))DEALLOCATE(DTFAT)				!39
			
				IF(ALLOCATED(WTS))DEALLOCATE(WTS)					!40
				IF(ALLOCATED(UARD_OBCN))DEALLOCATE(UARD_OBCN)		!41
				IF(ALLOCATED(XFLUX_OBC))DEALLOCATE(XFLUX_OBC)		!42
				IF(ALLOCATED(KH))DEALLOCATE(KH)						!43
		
				IF(ALLOCATED(EL))DEALLOCATE(EL)						!44
				IF(ALLOCATED(ET))DEALLOCATE(ET)						!45
				
				IF(ALLOCATED(D))DEALLOCATE(D)						!46
				IF(ALLOCATED(DT))DEALLOCATE(DT)						!47
				IF(ALLOCATED(DT1))DEALLOCATE(DT1)					!48
				IF(ALLOCATED(DTFA))DEALLOCATE(DTFA)					!49

				IF(ALLOCATED(VISCOFH))DEALLOCATE(VISCOFH)			!50
				

				
			END SUBROUTINE HYDRO_DEALLOC 

			
END MODULE MOD_HYDROVARS
