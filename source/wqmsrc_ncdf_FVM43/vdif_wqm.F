!Subroutine
! Subroutine VDIF_WQM()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   SUBROUTINE VDIF_WQM(F)
   
	USE MOD_PREC,ONLY: SP
	USE MOD_LIMS, ONLY : MTLOC, KBM1, MLOC, KB, KBM2,KBM1
	
    USE MOD_TGE, ONLY : 		&!
					ISONB		!	

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
!		 ISONB	&		!NODE MARKER = 0,1,2 
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
		,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
		!,H1	&		!BATHYMETRIC DEPTH 
		!,H	&			!BATHYMETRIC DEPTH 
		 ,D	&			!CURRENT DEPTH 
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
		 ,KH	!&		!TURBULENT DIFFUSIVITY
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

   USE MOD_SIZES, ONLY : NCP
   USE MOD_WQM,ONLY : NAC,AC,DLT,ZDFBCK
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MTLOC,KBM1,NCP) :: F
   REAL(SP), DIMENSION(0:MTLOC,KBM1) :: FF
   REAL(SP), DIMENSION(MLOC,KB)    :: VHF, VHPF
!  KURT GLAESEMANN 8 SEP 2009 - Fix dimensions on AF, CF, RAD KB!=NCP
   REAL(SP), DIMENSION(MLOC,KB)    :: AF, CF, RAD
   REAL(SP), DIMENSION(MLOC)       :: BENFLUX,WFSURF
   REAL(SP), DIMENSION(MLOC)          :: SOURCE1,SOURCE2,SOURCE3
   REAL(SP), DIMENSION(MLOC)          :: TBOT
   REAL(SP) :: FKH,UMOLPR
   REAL(SP) :: TEMPWUVBOT,TMP
   INTEGER  :: I,K,J,KI,JCON,II

   UMOLPR = ZDFBCK*1.E0
!  KURT GLAESEMANN - this was not initalized
  BENFLUX = 0

!----------------------------------------------------------------
!                                                                
!  the following section solves the equation               
!  dti*(kh*f')' -f=-fb
!                                                                
!----------------------------------------------------------------

   DO K = 2, KBM1
     DO I = 1, MLOC
       IF(D(I) > 0.0)THEN
         FKH=KH(I,K)
         AF(I,K-1)=-DLT*(FKH+UMOLPR)/(DZ(K-1)*DZZ(K-1)*D(I)*D(I))
         CF(I,K)=-DLT*(FKH+UMOLPR)/(DZ(K)*DZZ(K-1)*D(I)*D(I))
       ENDIF
     ENDDO
   ENDDO

   WFSURF = 0.0

!------------------------------------------------
!  Surface BCs; WFSURF
!----------------------------------------------- 

!  KURT GLAESEMANN 8 SEP 2009 -  MERGE ALL JCON LOOPS INTO ONE BIG ONE - FASTER
!  A BUNCH OF VARIABLES LOST A NCP DIMENSION
   DO JCON=1,NAC
     II = AC(JCON)
     DO I = 1, MLOC
       IF (D(I) > 0.0) THEN
         VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
         VHPF(I,1) = -DLT * WFSURF(I) / (-DZ(1)*D(I)) - F(I,1,II)
         VHPF(I,1) = VHPF(I,1) / (AF(I,1)-1.)
       ENDIF
     ENDDO
       
     DO K = 2, KBM2
       DO I = 1, MLOC
         IF(D(I) > 0.0) THEN
           VHPF(I,K)=1./ (AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
           VHF(I,K) = AF(I,K) * VHPF(I,K)
           VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-F(I,K,II))*VHPF(I,K)
         ENDIF
       ENDDO
     ENDDO

     DO K = 1, KBM1 
       DO I = 1, MLOC
         If (D(I) > 0.0)THEN
           FF(I,K) = F(I,K,II)
         ENDIF
       ENDDO
     ENDDO

     DO I = 1, MLOC
       IF (D(I) > 0.0 .AND. ISONB(I) /= 2) THEN
         FF(I,KBM1) = (CF(I,KBM1)*VHPF(I,KBM2)-FF(I,KBM1)   &
                         -DLT*BENFLUX(I)/(D(I)*DZ(KBM1)))/  &
                         (CF(I,KBM1)*(1.-VHF(I,KBM2))-1.)
       ENDIF
     ENDDO

     DO K = 2, KBM1
       KI = KB - K
       DO I = 1, MLOC
         IF (D(I) > 0.0 .AND. ISONB(I) /= 2) THEN
           FF(I,KI) = (VHF(I,KI)*FF(I,KI+1)+VHPF(I,KI))
         ENDIF
       ENDDO
     ENDDO

     DO I = 1, MLOC
       IF(D(I) > 0.0)THEN
         DO K = 1, KBM1
           F(I,K,II) = FF(I,K)
         ENDDO
       ENDIF
     ENDDO
   ENDDO

   RETURN
   END SUBROUTINE VDIF_WQM
!==============================================================================!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 
