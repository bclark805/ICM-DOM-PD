!Subroutine VISCOF_H() 

!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|

   SUBROUTINE VISCOF_H               

!------------------------------------------------------------------------------|
	USE MOD_LIMS, ONLY : KBM1, MTLOC, KB, MLOC
	
	USE MOD_TGE, ONLY : 		&!
			ISONB				&!
			,NTVE				&!		
			,NBVT				&!
			,NBVE				&!
			,NV

	USE MOD_HYDROVARS, ONLY:   &
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
		 ,VY	&		!Y-COORD AT GRID POINT
		!,ART	&		!AREA OF ELEMENT
		 ,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		!,ART2	&		!AREA OF ELEMENTS AROUND NODE
!		 ,NV	&		!NODE NUMBERING FOR ELEMENTS
		!,NBE	&		!INDICES OF ELMNT NEIGHBORS
		!,NTVE	&		!
		!,NTSN	&		!
!		 ,ISONB	&		!NODE MARKER = 0,1,2 
		!,ISBC	&		!
		!,ISBCE	&		!
		!,IEC	&		!
		!,IENODE &		!
		!,NBSN	&		!
		!,NIEC	&		!
		!,NTRG	&		!
!		 ,NBVE	&		!
!		 ,NBVT	&		!
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
		 ,VISCOFH	!&	!
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

   USE MOD_PREC, ONLY: SP		

   IMPLICIT NONE
! KURT GLAESEMANN 8 SEPT 2009, REMOVE MLOC DIMENSION - NOT NEEDED FOR PUPX,PUPY,PVPX,PVPY & VISCOFF
   REAL(SP) :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP) :: VISCOFF
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
   INTEGER  :: I,I1,IA,IB,J,J1,J2,K,JTMP

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!
   DO K=1,KBM1
     DO I=1,MLOC
       PUPX=0.0_SP
       PUPY=0.0_SP
       PVPX=0.0_SP
       PVPY=0.0_SP

! KURT GLAESEMANN 8 SEPT 2009 FOLD IN SPECIAL CASE FOR J=1 AND J=NTVE(I)
! HAVE TO LEAVE ISONB(I) /= 0 BELOW
       DO J=1,NTVE(I)
         I1=NBVE(I,J)
         JTMP=NBVT(I,J)
         J1=JTMP+1-(JTMP+1)/4*3
         J2=JTMP+2-(JTMP+2)/4*3
         X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
         Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
         X22=XC(I1)
         Y22=YC(I1)
         X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
         Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

         PUPX=PUPX+UU(I1,K)*(Y11-Y33)
         PUPY=PUPY+UU(I1,K)*(X33-X11)
         PVPX=PVPX+VV(I1,K)*(Y11-Y33)
         PVPY=PVPY+VV(I1,K)*(X33-X11)
         IF(ISONB(I) /= 0) THEN
           IF(J .EQ. 1) THEN
             PUPX=PUPX+UU(I1,K)*(VY(I)-Y11)
             PUPY=PUPY+UU(I1,K)*(X11-VX(I))
             PVPX=PVPX+VV(I1,K)*(VY(I)-Y11)
             PVPY=PVPY+VV(I1,K)*(X11-VX(I))
           ENDIF
           IF(J .EQ. NTVE(I)) THEN
             PUPX=PUPX+UU(I1,K)*(Y11-VY(I))
             PUPY=PUPY+UU(I1,K)*(VX(I)-X11)
             PVPX=PVPX+VV(I1,K)*(Y11-VY(I))
             PVPY=PVPY+VV(I1,K)*(VX(I)-X11)
           ENDIF
         ENDIF
       ENDDO

       PUPX=PUPX/ART1(I)
       PUPY=PUPY/ART1(I)
       PVPX=PVPX/ART1(I)
       PVPY=PVPY/ART1(I)
       TMP1=PUPX**2+PVPY**2
       TMP2=0.5_SP*(PUPY+PVPX)**2
       VISCOFF=SQRT(TMP1+TMP2)*ART1(I)
       
       VISCOFH(I,K) = VISCOFF

     ENDDO
   ENDDO  
    
   RETURN
   END SUBROUTINE VISCOF_H
!==============================================================================|
