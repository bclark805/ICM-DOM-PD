!Subroutine VERTVL()

!==============================================================================|
!   CALCULATE THE SIGMA COORDINATE VERTICAL VELOCITY FOR THE 3D MODE (omega)   |
!                                                  |
!   DETERMINED FROM EQUATION:                              |
!                                              !
!   d/dt(D) + d/dx(uD) + d/dy(uD) = d/sigma(omega)                             !
!==============================================================================|

   SUBROUTINE VERTVL         

!------------------------------------------------------------------------------|
    USE MOD_PREC, ONLY: SP
	USE MOD_LIMS, ONLY: NCV, MLOC,MTLOC, KBM1, KB
	
	USE MOD_TGE, ONLY : 	&!
				ISONB		&!		
				,NIEC		&!
				,NTRG		&!
				,DLTXE		&!
				,DLTYE
					
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
		 ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		!,ART2	&		!AREA OF ELEMENTS AROUND NODE
		!,NV	&		!NODE NUMBERING FOR ELEMENTS
		!,NBE	&		!INDICES OF ELMNT NEIGHBORS
		!,NTVE	&		!
		!,NTSN	&		!
!		 ,ISONB	&		!NODE MARKER = 0,1,2 
		!,ISBC	&		!
		!,ISBCE	&		!
		!,IEC	&		!
		!,IENODE &		!
		!,NBSN	&		!
!		,NIEC	&		!
!		,NTRG	&		!
		!,NBVE	&		!
		!,NBVT	&		!
		!,LISBCE_1	&	!LIST OF ELEMENTS WITH ISBCE=1
		!,LISBCE_2	&	!LIST OF ELEMENTS WITH ISBCE=2
		!,LISBCE_3	&	!LIST OF ELEMENTS WITH ISBCE=3
		!,DLTXC	&		!
		!,DLTYC	&		!
		!,DLTXYC	&	!
!		,DLTXE	&		!
!		,DLTYE	&		!
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
		 ,DZ	&		!DELTA-SIGMA VALUE
		!,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
		!,H1	&		!BATHYMETRIC DEPTH 
		,H	&			!BATHYMETRIC DEPTH 
		!,D	&			!CURRENT DEPTH 
		!,DT	&		!DEPTH AT PREVIOUS TIME STEP
		,DT1	&		!DEPTH AT PREVIOUS TIME STEP
		,EL	&		!CURRENT SURFACE ELEVATION
		,ET	&		!SURFACE ELEVATION AT PREVIOUS TIME STEP
		,DTFA	&		!ADJUSTED DEPTH FOR MASS CONSERVATION
		,UU	&		!X-VELOCITY
		,VV	&		!Y-VELOCITY
		!,UUT	&		!X-VELOCITY FROM PREVIOUS TIMESTEP
		!,VVT	&		!Y-VELOCITY FROM PREVIOUS TIMESTEP
		!,WWT	&		!Z-VELOCITY FROM PREVIOUS TIMESTEP
		!,WTST	&		!Vertical velocity in sigma from PREVIOUS TIMESTEP
		!,UARD_OBCNT	&!tykim
		!,XFLUX_OBCT	&!tykim
		!,DTFAT	&		!tykim
		!,TT_T	&		!tykim
		!,SALTT	&		!tykim
		 ,WTS!	&		!VERTICAL VELOCITY IN SIGMA SYSTEM
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

   USE MOD_WQM, ONLY : DLT

   IMPLICIT NONE 
   REAL(SP) :: XFLUX(MTLOC,KBM1)
   REAL(SP) :: DIJ,UIJ,VIJ,UN,EXFLUX,TMP1
   INTEGER  :: I,K,IA,IB,I1 ,J,JJ,J1,J2
!------------------------------------------------------------------------------|

!----------------------INITIALIZE FLUX-----------------------------------------!

   XFLUX = 0.0_SP

!----------------------ACCUMULATE FLUX-----------------------------------------!

   DO I=1,NCV
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     DIJ=DT1(I1)
     DO K=1,KBM1
       UIJ=UU(I1,K)
       VIJ=VV(I1,K)
       EXFLUX=DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
       XFLUX(IA,K)=XFLUX(IA,K)-EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)+EXFLUX
     ENDDO
   ENDDO

!-----------------------NULLIFY BOUNDARY FLUX----------------------------------!

      DO I=1,MLOC
!        KURT GLAESEMANN - 16 SEPT 2009 - Pull IF out of loop
!        DO K=1,KBM1
!          IF(ISONB(I) == 2) XFLUX(I,K)=0.0_SP  
!        ENDDO
         IF(ISONB(I) == 2) XFLUX(I,1:KBM1)=0.0_SP
      ENDDO

!---IF NO FRESH WATER INFLOW, OMEGA IS ZERO AT FREE SURFACE AND BOTTOM---------!

! KURT GLAESEMANN - 16 SEPT 2009 - WTS(*,KB) is set below
!   DO I=1,MLOC
!     WTS(I,1)=0.0_SP
!     WTS(I,KB)=0.0_SP
!   ENDDO
    WTS(1:MLOC,1)=0.0_SP

!--------------------------CALCULATE OMEGA-------------------------------------!

   DO I=1,MLOC
     DO K=1,KBM1
       WTS(I,K+1)=WTS(I,K)+DZ(K)*(XFLUX(I,K)/ART1(I)+(EL(I)-ET(I))/DLT)
     ENDDO
   ENDDO

!--------------------------ADJUST OMEGA----------------------------------------!
! IMPROVES MASS CONSERVATION

   DO I=1,MLOC
     IF(ABS(WTS(I,KB)) > 1.0E-8_SP)THEN
       IF(ISONB(I) /= 2)THEN

         TMP1=EL(I)*FLOAT(KBM1)-WTS(I,KB)*DLT/DZ(1)  !change of layer thickness due to vertical velocity
         
         TMP1=TMP1/FLOAT(KBM1)
         DTFA(I)=TMP1+H(I)  
         DO K=2,KB
           WTS(I,K)=WTS(I,K)-FLOAT(K-1)/FLOAT(KBM1)*WTS(I,KB)
         ENDDO
       ENDIF    
     ENDIF
   ENDDO

   RETURN
   END SUBROUTINE VERTVL
!==============================================================================|
