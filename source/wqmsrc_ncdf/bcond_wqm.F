!subroutine BCOND_WQM()
!subroutine BCOND_NUT_PERTURBATION()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   SUBROUTINE BCOND_WQM     
!==============================================================================|
!   Set Boundary Conditions on Water Quality                                   |
!==============================================================================|

!------------------------------------------------------------------------------|
  USE MOD_SIZES, ONLY :        &	!
			 NCP,            &  !
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
			 NOBTY  
  
   
	USE MOD_PREC , ONLY : SP
	USE MOD_LIMS, ONLY : KBM1, NUMQBC, IINT
	
	USE MOD_TGE, ONLY : NTSN, NBSN
	
	USE MOD_BCMAP, ONLY: IOBCN, I_OBC_N, WQOBC, NUT_TM
								
	USE MOD_HYDROVARS, ONLY: &
   		 GRAV	&		!
		!,PI	&		!
		!,PI2	&		!
		 ,ZERO	&		!
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
		 ,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
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
		,DZ	&		!DELTA-SIGMA VALUE
		!,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
		!,H1	&		!BATHYMETRIC DEPTH 
		 ,H	&			!BATHYMETRIC DEPTH 
		 ,D	&			!CURRENT DEPTH 
		 ,DT	&		!DEPTH AT PREVIOUS TIME STEP
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
		 ,UARD_OBCN	&	! tykim 
		 ,XFLUX_OBC	&	! tykim 
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
			! DTI        	&  			!!internal time step
			!,HORZMIX    	&   		!!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
			!,FILENUMBER	&			!!
			!,PREFF			&			!!
			!,INPDIR		&			!!
			!,GEOAREA		&			!!
			!,RIV_FILENUMBER	&			!!
            INFLOW_TYPE   	&			!!SPECIFIED RIVER INFLOW TYPE (edge/node) 
            ,POINT_ST_TYPE !	&			!!(calculated/specified)
            !,PNT_SOURCE    	&			!!point_source
            !,DAY				&
			!,in_jday		
   USE MOD_BCS, ONLY : &
         INODEQ		    &  !!LOCAL FRESH WATER INFLOW NODES
		,ICELLQ     	&  !!LOCAL FRESH WATER INFLOW ELEMENTS 
		,QBC_TM         &  !!TIME MAP FOR RIVER DISCHARGE DATA 

		,DWDIS   		&  !!WATER QUALITY DISCHARGE DATA
		,WDIS      		&  !!FRESH WATER QUALITY AT CURRENT TIME

		,N_ICELLQ			!node number of the edge for element edge source

		
	USE MOD_WQM, ONLY : C2,C2F,NAC,AC,JDAY,DLT
    USE MOD_PREC, ONLY : SP
    USE MOD_OBCS, ONLY:  &
		 NEXT_OBC ,       & !!INTERIOR NEIGHBOR OF OPEN BOUNDARY NODE
		!NEXT_OBC2,       & !!INTERIOR NEIGHBOR OF NEXT_OBC
		!TYPE_OBC ,       & !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
		!TYPE_OBC_GL,     & !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
		!IBCN       ,     & !!NUMBER OF EACH TYPE OF OBN IN LOCAL  DOM
		!IBCN_GL    ,     & !!NUMBER OF EACH TYPE OF OBN IN GLOBAL DOM
		!OBC_LST    ,     & !!MAPPING OF OPEN BOUNDARY ARRAYS TO EACH TYPE
		!NADJN_OBC  ,     & !!NUMBER OF ADJACENT OPEN BOUNDARY NODES TO OBN
		!ADJN_OBC   ,     & !!ADJACENT OBN's of OBN
		!NADJC_OBC  ,     & !!NUMBER OF ADJACENT OPEN BOUNDARY CELLS TO OBN
		!ADJC_OBC   ,     & !!ADJACENT OPEN BOUNDARY CELLS
		!NFLUXF_OBC ,     & !!NUMBER OF FLUX SEGMENTS TO OBN
		!FLUXF_OBC  ,     & !!FLUX FRACTION ON EACH SIDE OF OBN
		!NXOBC      ,     & !!INWARD POINTING X-NORMAL OF OBN
		!NYOBC      ,     & !!INWARD POINTING Y-NORMAL OF OBN
		 DLTN_OBC  ! ,     & !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY
		!NUT_OBC    ,     & !!NUTRIENT AT OPEN BOUNDARY
		!NUT_OBC_GL ,     & !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY
		
		!SEPARATE_OBC, &
		!SETUP_OBC 
   
   IMPLICIT NONE
   INTEGER  :: I,J,K,J1,J11,J22,NCON2,		&
				JCON,						&	!Index for activated constituents
				II
   REAL(SP)     :: UFACT,FACT
   REAL(SP) :: T2D,T2D_NEXT,T2D_OBC,XFLUX2D,TMP,TMAX,TMIN,RAMP_TS
   INTEGER  :: L1,L2,IERR
   REAL(SP) :: CC,CP,PERT,PERT_NEXT
   REAL(SP), ALLOCATABLE :: TTMP(:,:)
   REAL(SP) :: TY2D_N,TYC2F
   REAL(SP), ALLOCATABLE :: WQOBC_TMP(:,:,:)
!------------------------------------------------------------------------------|


!
!--SET CONDITIONS FOR FRESH WATER INFLOW---------------------------------------|
!
   IF(NUMQBC > 0)THEN
     ALLOCATE(WDIS(NUMQBC,NCP))     ;WDIS   = ZERO
     CALL BRACKET(QBC_TM,JDAY,L1,L2,FACT,UFACT,IERR)
      WDIS(:,:) = DWDIS(:,:,L1)      !UFACT*DWDIS(:,L1) + FACT*DWDIS(:,L2)
   ENDIF

   IF(POINT_ST_TYPE == 'specified') THEN
     IF(NUMQBC > 0) THEN
       IF(INFLOW_TYPE == 'node') THEN
         DO JCON =1, NAC
	          II = AC(JCON)
	          DO K=1,KBM1
	          DO I=1,NUMQBC
               J11=INODEQ(I)
               C2F(J11,K,II) = WDIS(I,II)
            ENDDO
            ENDDO
         ENDDO
       ELSEIF(INFLOW_TYPE == 'edge') THEN
         DO JCON = 1,NAC
            II = AC(JCON)
            DO K=1,KBM1
            DO I=1,NUMQBC
               J11=N_ICELLQ(I,1)  		!first node of the element edge
               J22=N_ICELLQ(I,2)  		!second node of the element edge
               C2F(J11,K,II)=WDIS(I,II) !Wen Long has doubt here, shouldn' it be divided by 2 so that 
										!	it is distirbuted evenly by the two nodes connecting the edge?
               C2F(J22,K,II)=WDIS(I,II) ! !Wen Long has doubt here, shouldn' it be divided by 2 so that 
										!	it is distirbuted evenly by the two nodes connecting the edge?
            ENDDO
            ENDDO
         ENDDO
       ENDIF
     ENDIF
   ENDIF

   IF(NUMQBC > 0)DEALLOCATE(WDIS)
      
   IF(IOBCN > 0) THEN
!
!  SET WATER QUALITY CONDITIONS ON OUTER BOUNDARY
!
     ALLOCATE(WQOBC_TMP(1:IOBCN,1:KBM1,1:NCP))
     CALL BRACKET(NUT_TM,THOUR,L1,L2,FACT,UFACT,IERR)
     IF(IERR==-1)THEN
       WQOBC_TMP = 0.0_SP
     ELSE

       WQOBC_TMP(1:IOBCN,1:KBM1,1:NCP) =  &
               UFACT*WQOBC(1:IOBCN,1:KBM1,L1,1:NCP)&
            +  FACT*WQOBC(1:IOBCN,1:KBM1,L2,1:NCP)
     ENDIF

     RAMP_TS=TANH(FLOAT(IINT)/FLOAT(8640+1))
!
!
     DO JCON=1,NAC
       II = AC(JCON)
       ALLOCATE(TTMP(NOBTY,KBM1)); TTMP=0.0_SP
       DO I=1,IOBCN
         J=I_OBC_N(I)
         J1=NEXT_OBC(I)
         T2D=0.0_SP
         T2D_NEXT=0.0_SP
         XFLUX2D=0.0_SP
         DO K=1,KBM1
            T2D=T2D+C2(J,K,II)*DZ(K)
            T2D_NEXT=T2D_NEXT+C2F(J1,K,II)*DZ(K)
            XFLUX2D=XFLUX2D+XFLUX_OBC(I,K)*DZ(K)
!                IF(II.eq.2) write(2236,*)T2D_NEXT,C2F(J1,K,II),DZ(K)
         ENDDO
      
         IF(UARD_OBCN(I) > 0.0) THEN
            TMP=XFLUX2D+T2D*UARD_OBCN(I)
!tykim
            T2D_OBC=(T2D*DT(J)-TMP*DLT/ART1(J))/D(J) 
!             DTI=DTE*FLOAT(ISPLIT) -- In FVCOM 

!            CALL BCOND_NUT_PERTURBATION(T2D_NEXT,T2D,TTMP,I,J,J1)
! Perturbation case 1 TTMP(I,K) = TF1(J1,K) - T2D_NEXT
!            DO K=1,KBM1
!               if(II.eq.27) write(2238,*)'before',C2F(J,K,II),T2D_OBC,T2D_NEXT
!                C2F(J,K,II)=T2D_OBC+(C2F(J1,K,II)-T2D_NEXT)
!                if(II.eq.27) write(2238,*)'after',C2F(J,K,II)
!            ENDDO
! Perturbation case 3
            CC = SQRT(GRAV*H(J))*DLT/DLTN_OBC(I)
            CP = CC + 1.0_SP
            DO K=1,KBM1
               PERT_NEXT = C2F(J1,K,II) - T2D_NEXT
               TYC2F=C2F(J1,K,II)
               TY2D_N=T2D_NEXT
               PERT      = C2(J,K,II)  - T2D

               C2F(J,K,II)= T2D_OBC + (CC*PERT_NEXT + PERT*(1.0_SP - DLT/10800.0_SP))/CP

! Perturbation case 2
!               C2F(J,K,II)= T2D_OBC + (CC*PERT_NEXT + PERT)/CP
            ENDDO
! The end of Perturbation
      
            DO K=1,KBM1
               TMAX = MAXVAL(C2(NBSN(J,1:NTSN(J)),K,II))
               TMIN = MINVAL(C2(NBSN(J,1:NTSN(J)),K,II))

           IF(K == 1)THEN
            TMAX = MAX(TMAX,(C2(J,K,II)*DZ(K+1)+C2(J,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1)))
            TMIN = MIN(TMIN,(C2(J,K,II)*DZ(K+1)+C2(J,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1)))
           ELSEIF(K == KBM1)THEN
            TMAX = MAX(TMAX,(C2(J,K,II)*DZ(K-1)+C2(J,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1)))
            TMIN = MIN(TMIN,(C2(J,K,II)*DZ(K-1)+C2(J,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1)))
           ELSE
            TMAX = MAX(TMAX,(C2(J,K,II)*DZ(K-1)+C2(J,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1)), &
                   (C2(J,K,II)*DZ(K+1)+C2(J,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1)))
            TMIN = MIN(TMIN,(C2(J,K,II)*DZ(K-1)+C2(J,K-1,II)*DZ(K))/(DZ(K)+DZ(K-1)), &
                   (C2(J,K,II)*DZ(K+1)+C2(J,K+1,II)*DZ(K))/(DZ(K)+DZ(K+1)))
           ENDIF

           IF(TMIN-C2F(J,K,II) > 0.0_SP) C2F(J,K,II) = TMIN
           IF(C2F(J,K,II)-TMAX > 0.0_SP) C2F(J,K,II) = TMAX

!           IF(TMIN-C2F(J,K,II) > C2F(J,K,II)*0.0001) C2F(J,K,II) = TMIN
!           IF(C2F(J,K,II)-TMAX > C2F(J,K,II)*0.0001) C2F(J,K,II) = TMAX

           ENDDO
 
        ELSE

           DO K=1,KBM1

!            C2F(J,K,II) = C2(J,K,II)-0.5*RAMP_TS*(C2(J,K,II)-NUT_OBC_GL(II,I))

!           ENDDO

              IF(IERR.NE.-1)THEN
!         if(II.eq.4.and.I.eq.1)write(666,*)'tykim',I,K,II,WQOBC_TMP(I,K,II),C2(J,K,II)
                 C2F(J,K,II)=C2(J,K,II)-0.014*RAMP_TS*(C2(J,K,II)-WQOBC_TMP(I,K,II))

              ELSE
                 C2F(J,K,II) = C2(J,K,II)
              ENDIF

           ENDDO

        ENDIF
       ENDDO
       DEALLOCATE(TTMP)
      ENDDO !!OUTER LOOP OVER WQ VARIABLES
       DEALLOCATE(WQOBC_TMP)
   ENDIF

!
!--SET BOUNDARY CONDITIONS-----------------------------------------------------|
!
   C2F(0,:,:) = 0.0

   RETURN
   END SUBROUTINE BCOND_WQM

!==============================================================================|
 
!==============================================================================|
   SUBROUTINE BCOND_NUT_PERTURBATION(T2D_NEXT,T2D,TTMP,I,J,J1)
!==============================================================================|
! Calculate the OBC for temperature perturbation                               |
!==============================================================================|
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
            !MGL,            &  !
			!NGL,            &  !
			!OBCGL,          &  !
			 NOBTY  
  
	USE MOD_PREC, ONLY : SP
	USE MOD_LIMS, ONLY : KBM1, IINT
	IMPLICIT NONE

!   INTEGER :: I1,I2,J,JN
   INTEGER :: I,J,J1,J2,K
   REAL(SP):: CC,CP,MU,CL
   REAL(SP):: PERT_NEXT,PERT,T2D_NEXT,T2D
   REAL(SP):: T2D_NEXT1,TM12D_NEXT2,TM12D_NEXT1,TM22D_NEXT1
   REAL(SP):: TTMP(NOBTY,KBM1)

!   SELECT CASE(TYPE_TSOBC)
   
!   CASE(1)
!     DO K=1,KBM1
!       TTMP(I,K) = TF1(J1,K) - T2D_NEXT
!     ENDDO  
!   CASE(2)
!     CC = SQRT(GRAV*H(J))*DTI/DLTN_OBC(I)
!     CP = CC + 1.0_SP
!     DO K=1,KBM1
!       PERT_NEXT = TF1(J1,K) - T2D_NEXT
!       PERT      = T1(J,K) - T2D
!       TTMP(I,K) = (CC*PERT_NEXT + PERT)/CP
!     ENDDO
!   CASE(3)
!     CC = SQRT(GRAV*H(J))*DTI/DLTN_OBC(I)
!     CP = CC + 1.0_SP
!     DO K=1,KBM1
!       PERT_NEXT = TF1(J1,K) - T2D_NEXT
!       PERT      = T1(J,K) - T2D
!       TTMP(I,K) = (CC*PERT_NEXT + PERT*(1.0_SP - DTI/10800.0_SP))/CP
!     ENDDO
!   CASE(4)
!     J2 = NEXT_OBC2(I)
!     T2D_NEXT1  =0.0_SP
!     TM12D_NEXT2=0.0_SP
!     TM12D_NEXT1=0.0_SP
!     TM22D_NEXT1=0.0_SP
!     DO K=1,KBM1
!       T2D_NEXT1  =T2D_NEXT1  +T1(J1,K)*DZ(K)
!       TM12D_NEXT2=TM12D_NEXT2+T1M1(J2,K)*DZ(K)
!       TM12D_NEXT1=TM12D_NEXT1+T1M1(J,K)*DZ(K)
!       TM22D_NEXT1=TM22D_NEXT1+T1M2(J1,K)*DZ(K)
!     ENDDO
!
!     DO K=1,KBM1
!       CL = ((T1M2(J1,K)-TM22D_NEXT1)-(T1(J1,K)-T2D_NEXT1))/   &
!            ((T1(J1,K)-T2D_NEXT1)+(T1M2(J1,K)-TM22D_NEXT1)     &
!	           -2.0*(T1M1(J2,K)-TM12D_NEXT2))
!       IF(CL >= 1.0)THEN
!         MU = 1.0
!       ELSEIF(CL > 0.0 .AND. CL < 1.0)THEN
!         MU = CL
!       ELSE
!         MU = 0.0
!       ENDIF

!       TTMP(I,K)=((T1M1(J,K)-TM12D_NEXT1)*(1.0-MU)     &
!                 +2.0*MU*(T1(J1,K)-T2D_NEXT1))/(1.0+MU)
!     ENDDO

!   END SELECT
     
   RETURN
   END SUBROUTINE BCOND_NUT_PERTURBATION
!========================================================================
