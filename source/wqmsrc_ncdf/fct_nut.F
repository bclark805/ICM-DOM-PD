!subroutine FCT_NUT()

!==============================================================================|
!  FLUX CONTROL FOR NUTRIENTS                                                        |
!==============================================================================|

   SUBROUTINE FCT_NUT
!#  if defined (WET_DRY)

!==============================================================================|
    USE MOD_LIMS, ONLY: MLOC, KBM1, NUMQBC, NUMPNT
		
	USE MOD_TGE, ONLY : 	&	!
		 NTSN				&	!
		,NBSN
							
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
!		 ,NBSN	&		!
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
		 DZ 	&		!DELTA-SIGMA VALUE
		!,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
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
		,THOUR	!&		!
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
            INFLOW_TYPE   	!&			!!SPECIFIED RIVER INFLOW TYPE (edge/node) 
            !,POINT_ST_TYPE 	&			!!(calculated/specified)
            !,PNT_SOURCE    	&			!!point_source
            !,DAY				&
			!,in_jday		
   USE MOD_BCS, ONLY :  &
         INODEQ		    &  !!LOCAL FRESH WATER INFLOW NODES
		,INOPNT      	   !!LOCAL NON-POINT SOURCE NODES

   USE MOD_WQM, ONLY : C2,C2F,NAC,AC
   USE MOD_PREC, ONLY: SP	
   USE MOD_PAR, ONLY : NGID

   USE MOD_BCMAP, ONLY: IOBCN, I_OBC_N
   
   IMPLICIT NONE
   REAL(SP):: NUTMAX,NUTMIN
   INTEGER :: I,J,K,JCON,II

!==============================================================================|
  DO JCON=1,NAC
   II=AC(JCON)
   DO I=1,MLOC
      IF(IOBCN > 0)THEN
        DO J=1,IOBCN
          IF(I == I_OBC_N(J)) GO TO 200
       ENDDO
      ENDIF

      IF(NUMQBC > 0 )then
        DO J=1,NUMQBC
           IF(INFLOW_TYPE == 'node') THEN
             IF (I == INODEQ(J)) GO TO 200
           ENDIF
        ENDDO
      ENDIF

      IF(NUMPNT > 0) THEN
       DO J=1,NUMPNT
         IF (I == INOPNT(J)) GO TO 200
       ENDDO
      ENDIF

	  DO K=1,KBM1

       NUTMAX = MAXVAL(C2(NBSN(I,1:NTSN(I)),K,II))
       NUTMIN = MINVAL(C2(NBSN(I,1:NTSN(I)),K,II))

       IF(K == 1)THEN
         NUTMAX = MAX(NUTMAX,(C2(I,K,II)*DZ(K+1)+C2(I,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1)))
         NUTMIN = MIN(NUTMIN,(C2(I,K,II)*DZ(K+1)+C2(I,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1)))
       ELSE IF(K == KBM1)THEN
         NUTMAX = MAX(NUTMAX,(C2(I,K,II)*DZ(K-1)+C2(I,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1)))
         NUTMIN = MIN(NUTMIN,(C2(I,K,II)*DZ(K-1)+C2(I,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1)))
       ELSE
         NUTMAX = MAX(NUTMAX,(C2(I,K,II)*DZ(K-1)+C2(I,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1)), &
                 (C2(I,K,II)*DZ(K+1)+C2(I,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1)))
         NUTMIN = MIN(NUTMIN,(C2(I,K,II)*DZ(K-1)+C2(I,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1)), &
                 (C2(I,K,II)*DZ(K+1)+C2(I,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1)))
       ENDIF

!tykim commented out two lines below
	    IF(NUTMIN-C2F(I,K,II) >  0.0_SP ) THEN

		 C2F(I,K,II) = NUTMIN

		
       ENDIF
	    IF(C2F(I,K,II)-NUTMAX >  0.0_SP ) THEN

!
		C2F(I,K,II) = NUTMAX

		
	   ENDIF

      ENDDO

 200 CONTINUE
	
   ENDDO
   ENDDO
   RETURN
!#  endif
   END SUBROUTINE FCT_NUT
!==============================================================================|


