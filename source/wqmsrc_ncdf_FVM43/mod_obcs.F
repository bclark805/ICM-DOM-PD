MODULE MOD_OBCS
   
	USE MOD_PREC, ONLY: SP!
   
	USE MOD_BCMAP, ONLY :  	&
		IOBCN,				&!
		IOBCN_GL,			&!
		I_OBC_N				!
		
   IMPLICIT NONE
   SAVE

!--Open Boundary Types, Lists, Pointers
   !WLong moved the following to bcmap.F
   !INTEGER               :: IOBCN_GL         !!GLOBAL NUMBER OF OPEN BOUNDARY NODES   !move to bcmap.F
   !INTEGER               :: IOBCN            !!LOCAL NUMBER OF OPEN BOUNDARY NODES    !move to bcmap.F
   !INTEGER,  ALLOCATABLE :: I_OBC_GL(:)      !!GLOBAL ID OF OPEN BOUNDARY NODES		  !move to bcmap.F
   !INTEGER,  ALLOCATABLE :: I_OBC_N(:)       !!OPEN BOUNDARY NODE LIST				  !move to bcmap.F
   
   INTEGER,  ALLOCATABLE :: NEXT_OBC(:)      !!INTERIOR NEIGHBOR OF OPEN BOUNDARY NODE
   INTEGER,  ALLOCATABLE :: NEXT_OBC2(:)     !!INTERIOR NEIGHBOR OF NEXT_OBC
   INTEGER,  ALLOCATABLE :: TYPE_OBC(:)      !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION) !never allocated
   INTEGER,  ALLOCATABLE :: TYPE_OBC_GL(:)   !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION) !never allocated
   INTEGER               :: IBCN(5)          !!NUMBER OF EACH TYPE OF OBN IN LOCAL  DOM
   INTEGER               :: IBCN_GL(5)       !!NUMBER OF EACH TYPE OF OBN IN GLOBAL DOM
   INTEGER,  ALLOCATABLE :: OBC_LST(:,:)     !!MAPPING OF OPEN BOUNDARY ARRAYS TO EACH TYPE
   INTEGER,  ALLOCATABLE :: NADJN_OBC(:)     !!NUMBER OF ADJACENT OPEN BOUNDARY NODES TO OBN
   INTEGER,  ALLOCATABLE :: ADJN_OBC(:,:)    !!ADJACENT OBN's of OBN
   INTEGER,  ALLOCATABLE :: NADJC_OBC(:)     !!NUMBER OF ADJACENT OPEN BOUNDARY CELLS TO OBN
   INTEGER,  ALLOCATABLE :: ADJC_OBC(:,:)    !!ADJACENT OPEN BOUNDARY CELLS

!--Open Boundary Metrics
   INTEGER,  ALLOCATABLE :: NFLUXF_OBC(:)    !!NUMBER OF FLUX SEGMENTS TO OBN   !never used
   REAL(SP), ALLOCATABLE :: FLUXF_OBC(:,:)   	 !!FLUX FRACTION ON EACH SIDE OF OBN	!never used
   REAL(SP), ALLOCATABLE :: NXOBC(:)         	 !!INWARD POINTING X-NORMAL OF OBN
   REAL(SP), ALLOCATABLE :: NYOBC(:)         	 !!INWARD POINTING Y-NORMAL OF OBN
   REAL(SP), ALLOCATABLE :: DLTN_OBC(:)      	 !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY
! tykim
   REAL(SP), ALLOCATABLE :: NUT_OBC(:,:)     	 !!NUTRIENT AT OPEN BOUNDARY							!never used
   REAL(SP), ALLOCATABLE :: NUT_OBC_GL(:,:)  	 !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY	!never used

		   
   CONTAINS

	!subroutines:

	!	subroutine SEPARATE_OBC()
	!	subroutine SETUP_OBC()
	!	subroutine OBCS_DEALLOC()
	!functions:
	!	
   
!==========================================================================|

!==========================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
   SUBROUTINE SEPARATE_OBC  !Function never used in this code

!------------------------------------------------------------------------------|
! Accumulate separately the amounts of nodes for 11 types of open boundaries   |
!------------------------------------------------------------------------------|

   IMPLICIT NONE

   INTEGER :: I,I1,I2,I3,I4,I5,II,J

   IBCN = 0
   IBCN_GL = 0

   DO I = 1, IOBCN_GL
     IF(TYPE_OBC_GL(I) == 1 .OR. TYPE_OBC_GL(I) == 2)  IBCN_GL(1) = IBCN_GL(1) + 1
     IF(TYPE_OBC_GL(I) == 3 .OR. TYPE_OBC_GL(I) == 4)  IBCN_GL(2) = IBCN_GL(2) + 1
     IF(TYPE_OBC_GL(I) == 5 .OR. TYPE_OBC_GL(I) == 6)  IBCN_GL(3) = IBCN_GL(3) + 1
     IF(TYPE_OBC_GL(I) == 7 .OR. TYPE_OBC_GL(I) == 8)  IBCN_GL(4) = IBCN_GL(4) + 1
     IF(TYPE_OBC_GL(I) == 9 .OR. TYPE_OBC_GL(I) == 10) IBCN_GL(5) = IBCN_GL(5) + 1
   ENDDO

   DO I = 1, IOBCN
     IF(TYPE_OBC(I) == 1 .OR. TYPE_OBC(I) == 2)  IBCN(1) = IBCN(1) + 1
     IF(TYPE_OBC(I) == 3 .OR. TYPE_OBC(I) == 4)  IBCN(2) = IBCN(2) + 1
     IF(TYPE_OBC(I) == 5 .OR. TYPE_OBC(I) == 6)  IBCN(3) = IBCN(3) + 1
     IF(TYPE_OBC(I) == 7 .OR. TYPE_OBC(I) == 8)  IBCN(4) = IBCN(4) + 1
     IF(TYPE_OBC(I) == 9 .OR. TYPE_OBC(I) == 10) IBCN(5) = IBCN(5) + 1
   ENDDO

   I1 = 0
   I2 = 0
   I3 = 0
   I4 = 0
   I5 = 0


   ALLOCATE(OBC_LST(5,MAXVAL(IBCN))) ; OBC_LST = 0

   DO I=1,IOBCN
     IF(TYPE_OBC(I) == 1 .OR. TYPE_OBC(I) == 2)THEN
       I1 = I1 + 1
       OBC_LST(1,I1) = I
     ELSEIF(TYPE_OBC(I) == 3 .OR. TYPE_OBC(I) == 4)THEN
       I2 = I2 + 1
       OBC_LST(2,I2) = I
     ELSEIF(TYPE_OBC(I) == 5 .OR. TYPE_OBC(I) == 6)THEN
       I3 = I3 + 1
       OBC_LST(3,I3) = I
     ELSEIF(TYPE_OBC(I) == 7 .OR. TYPE_OBC(I) == 8)THEN
       I4 = I4 + 1
       OBC_LST(4,I4) = I
     ELSEIF(TYPE_OBC(I) == 9 .OR. TYPE_OBC(I) == 10)THEN
       I5 = I5 + 1
       OBC_LST(5,I5) = I
     ENDIF
   ENDDO

   RETURN
   END SUBROUTINE SEPARATE_OBC
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
   SUBROUTINE SETUP_OBC

	USE MOD_LIMS, ONLY: MYID
	
    USE MOD_TGE, ONLY : &!
		ISONB,			&!
		NTVE,			&!
		NTSN,			&!
		NBSN,			&!
		NBVE,			&!
		NV

!------------------------------------------------------------------------------!
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
		 XC	&		!X-COORD AT FACE CENTER 
		,YC	&		!Y-COORD AT FACE CENTER
		,VX	&		!X-COORD AT GRID POINT
		,VY			!Y-COORD AT GRID POINT
		!,ART	&		!AREA OF ELEMENT
		!,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		!,ART2	&		!AREA OF ELEMENTS AROUND NODE
!		,NV	&		!NODE NUMBERING FOR ELEMENTS
		!,NBE	&		!INDICES OF ELMNT NEIGHBORS
!		,NTVE	&		!
		!,NTSN	&		!
!		,ISONB	&		!NODE MARKER = 0,1,2 
		!,ISBC	&		!
		!,ISBCE	&		!
		!,IEC	&		!
		!,IENODE &		!
!		 ,NBSN	&		!
		!,NIEC	&		!
		!,NTRG	&		!
!		,NBVE	!&		!
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

   
   IMPLICIT NONE

   REAL(SP) :: DXC,DYC,DXN,DYN,CROSS,E1,E2,DOTMAX,DOT,DX,DY,DXN_TMP,DYN_TMP
   INTEGER  :: I,J,JJ,INODE,JNODE,I1,I2,IC,N1,N2,N3
   LOGICAL  :: DEBUG

   REAL(SP), ALLOCATABLE :: NXOBC_TMP(:),NYOBC_TMP(:)

!------------------------------------------------------------------------------!

!--Determine Adjacent Open Boundary Points-------------------------------------!

   ALLOCATE(NADJN_OBC(IOBCN))  ; NADJN_OBC = 0
   ALLOCATE(ADJN_OBC(IOBCN,2)) ; ADJN_OBC = 0

   ALLOCATE(NEXT_OBC(IOBCN))           ;NEXT_OBC   = 0
   ALLOCATE(NEXT_OBC2(IOBCN))          ;NEXT_OBC2  = 0


   DO I=1,IOBCN
     INODE = I_OBC_N(I)
     DO J=1,NTSN(INODE)-1
       JNODE = NBSN(INODE,J)
       IF(ISONB(JNODE) == 2 .AND. INODE /= JNODE)THEN
         NADJN_OBC(I) = NADJN_OBC(I) + 1
         ADJN_OBC(I,NADJN_OBC(I)) = JNODE
       ENDIF
     ENDDO
   ENDDO


   DO I=1,IOBCN
     IF(NADJN_OBC(I) == 0)THEN
       WRITE(*,*)'NO ADJACENT NODE FOUND FOR BOUNDARY NODE',I
       WRITE(*,*)'IN PROCESSOR',MYID
       CALL PSTOP
     ENDIF
   ENDDO

!--Determine Adjacent Cells-(Nonlinear Only)-----------------------------------!
!--Simultaneously Determine INWARD Pointing Normal NXOBC,NYOBC                 !

   ALLOCATE(NADJC_OBC(IOBCN))  ; NADJC_OBC = 0
   ALLOCATE(ADJC_OBC(IOBCN,2)) ; ADJC_OBC = 0
   ALLOCATE(NXOBC(IOBCN)) ; NXOBC = 0
   ALLOCATE(NYOBC(IOBCN)) ; NYOBC = 0
   ALLOCATE(NXOBC_TMP(IOBCN)) ; NXOBC_TMP = 0
   ALLOCATE(NYOBC_TMP(IOBCN)) ; NYOBC_TMP = 0

   DO I=1,IOBCN
     I1 = I_OBC_N(I)

     !!Mark First Cell on Boundary Edge Adjacent to Node I
     I2 = ADJN_OBC(I,1)
     DO J = 1, NTVE(I1)
       IC = NBVE(I1,J)
       N1 = NV(IC,1) ; N2 = NV(IC,2) ; N3 = NV(IC,3)
       IF( N1-I2 == 0 .OR. N2-I2 == 0 .OR. N3-I2 == 0)THEN
         DXN = VX(I2)-VX(I1) ; DYN = VY(I2)-VY(I1)
         DXC = XC(IC)-VX(I1) ; DYC = YC(IC)-VY(I1)

         CROSS = SIGN(1.0_SP,DXC*DYN - DYC*DXN)
         NXOBC_TMP(I) =  CROSS*DYN/SQRT(DXN**2 +DYN**2)
         NYOBC_TMP(I) = -CROSS*DXN/SQRT(DXN**2 +DYN**2)
         NXOBC(I) = NXOBC_TMP(I)
         NYOBC(I) = NYOBC_TMP(I)
           NADJC_OBC(I) = NADJC_OBC(I) + 1
           ADJC_OBC(I,NADJC_OBC(I)) = IC
       ENDIF
     ENDDO

     IF(NADJN_OBC(I) > 1)THEN
     I2 = ADJN_OBC(I,2)
     DO J = 1, NTVE(I1)
       IC = NBVE(I1,J)
       N1 = NV(IC,1) ; N2 = NV(IC,2) ; N3 = NV(IC,3)
       IF( N1-I2 == 0 .OR. N2-I2 == 0 .OR. N3-I2 == 0)THEN
         DXN = VX(I2)-VX(I1) ; DYN = VY(I2)-VY(I1)
         DXC = XC(IC)-VX(I1) ; DYC = YC(IC)-VY(I1)

         CROSS = SIGN(1.0_SP,DXC*DYN - DYC*DXN)
         NXOBC_TMP(I) = NXOBC_TMP(I) + CROSS*DYN/SQRT(DXN**2 +DYN**2)
         NYOBC_TMP(I) = NYOBC_TMP(I) - CROSS*DXN/SQRT(DXN**2 +DYN**2)
         NXOBC(I) = NXOBC_TMP(I)/SQRT(NXOBC_TMP(I)**2 + NYOBC_TMP(I)**2)
         NYOBC(I) = NYOBC_TMP(I)/SQRT(NXOBC_TMP(I)**2 + NYOBC_TMP(I)**2)

           NADJC_OBC(I) = NADJC_OBC(I) + 1
           ADJC_OBC(I,NADJC_OBC(I)) = IC
       ENDIF
     ENDDO
     ENDIF
   ENDDO

   DEALLOCATE(NXOBC_TMP,NYOBC_TMP)

!--Determine 1st Layer Neighbor for Open Boundary Points-----------------------!
!--Node Chosen is Node That is Connected to OBC Node and is Oriented           !
!--Most Normal to the Boundary.  It is not Necessarily the Closest Node        !
!--Determine also DLTN_OBC, the normal component of the distance between       !
!--Next_obc and the open boundary node                                         !

   DO I=1,IOBCN
     DOTMAX =  -1.0
     INODE = I_OBC_N(I)
     DO J=1,NTSN(INODE)-1
       JNODE = NBSN(INODE,J)
       IF(ISONB(JNODE) /= 2 .AND. INODE /= JNODE)THEN
         DXN_TMP = VX(JNODE)-VX(INODE)
         DYN_TMP = VY(JNODE)-VY(INODE)

         DXN = DXN_TMP/SQRT(DXN_TMP**2 + DYN_TMP**2)
         DYN = DYN_TMP/SQRT(DXN_TMP**2 + DYN_TMP**2)
         DOT = DXN*NXOBC(I) + DYN*NYOBC(I)
         IF(DOT > DOTMAX)THEN
           DOTMAX = DOT
           NEXT_OBC(I) = JNODE
         ENDIF
       ENDIF
     ENDDO
   ENDDO

!--Determine DLTN_OBC----------------------------------------------------------!
   ALLOCATE(DLTN_OBC(IOBCN))
   DO I=1,IOBCN
     I1 = I_OBC_N(I)
     I2 = NEXT_OBC(I)

     DX = VX(I2)-VX(I1)
     DY = VY(I2)-VY(I1)
     DLTN_OBC(I) = ABS(DX*NXOBC(I) + DY*NYOBC(I))
   ENDDO

   RETURN
   END SUBROUTINE SETUP_OBC
   
   
   SUBROUTINE OBCS_DEALLOC
   
    !WLong moved these to mod_obcs.F
     IF(ALLOCATED(TYPE_OBC_GL)) DEALLOCATE(TYPE_OBC_GL)  !never allocated
     IF(ALLOCATED(NUT_OBC_GL)) DEALLOCATE(NUT_OBC_GL)	 !never allocated
     IF(ALLOCATED(NUT_OBC)) DEALLOCATE(NUT_OBC)			 !never allocated
    
     IF(ALLOCATED(OBC_LST)) DEALLOCATE(OBC_LST)			 !mod_obcs
	 
     IF(ALLOCATED(NADJN_OBC)) DEALLOCATE(NADJN_OBC)		 !mod_obcs
     IF(ALLOCATED(ADJN_OBC)) DEALLOCATE(ADJN_OBC)		 !mod_obcs
     IF(ALLOCATED(NEXT_OBC)) DEALLOCATE(NEXT_OBC)		 !mod_obcs
     IF(ALLOCATED(NEXT_OBC2)) DEALLOCATE(NEXT_OBC2)		 !mod_obcs
     IF(ALLOCATED(NADJC_OBC)) DEALLOCATE(NADJC_OBC)		 !mod_obcs
     IF(ALLOCATED(ADJC_OBC)) DEALLOCATE(ADJC_OBC)		 !mod_obcs
     IF(ALLOCATED(NXOBC)) DEALLOCATE(NXOBC)				 !mod_obcs
     IF(ALLOCATED(NYOBC)) DEALLOCATE(NYOBC)				 !mod_obcs
     IF(ALLOCATED(DLTN_OBC)) DEALLOCATE(DLTN_OBC)		 !mod_obcs
	 
   
   END SUBROUTINE OBCS_DEALLOC
   
!==============================================================================|
END MODULE MOD_OBCS
