MODULE MOD_TGE

	USE MOD_LIMS, ONLY: MX_NBR_ELEM,	&!
						MTLOC, 			&!
						NTLOC, 			&!
						MYID, 			&!
						MLOC, 			&!
						NLOC, 			&!
						NE,				&!
						NCT, 			&!
						NCV_I, 			&!
						NCV, 			&!
						NISBCE_1, 		&!
						NISBCE_2, 		&!
						NISBCE_3 

    USE MOD_HYDROVARS, ONLY: & 
   		!GRAV	&		!
		!,PI	&		!
		!,PI2	&		!
		ZERO	&		!
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
		,XC	&		!X-COORD AT FACE CENTER 
		,YC	&		!Y-COORD AT FACE CENTER
		,VX	&		!X-COORD AT GRID POINT
		,VY	!&		!Y-COORD AT GRID POINT
		!,ART	&		!AREA OF ELEMENT
		!,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		!,ART2	&		!AREA OF ELEMENTS AROUND NODE
!		,NV	&		!NODE NUMBERING FOR ELEMENTS
!		,NBE	&		!INDICES OF ELMNT NEIGHBORS
!		,NTVE	&		!
!		,NTSN	&		!
!		,ISONB	&		!NODE MARKER = 0,1,2 
!		,ISBCE 	&		!
!		,ISBC	&		!	
!		,IEC	&		!
!		,IENODE &		!
!		,NBSN	&		!
!		,NIEC	&		!
!		,NTRG	&		!
!		,NBVE	&		!
!		,NBVT	&		!
!		,LISBCE_1	&	!LIST OF ELEMENTS WITH ISBCE=1
!		,LISBCE_2	&	!LIST OF ELEMENTS WITH ISBCE=2
!		,LISBCE_3	&	!LIST OF ELEMENTS WITH ISBCE=3
!		,DLTXC	&		!
!		,DLTYC	&		!
!		,DLTXYC	&		!
!		,DLTXE	&		!
!		,DLTYE	&		!
!		,DLTXYE	&		!
!		,SITAC	&		!
!		,SITAE	&		!
!		,XIJC	&		!
!		,YIJC 	&		!
!		,XIJE	&		!
!		,YIJE 	&		!
!		,EPOR	!&		!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
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

	!Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
	USE MOD_CONTROL, ONLY : 		&
			!SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			MSR       ! 	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
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
#  if defined (MULTIPROCESSOR)
	USE MOD_PAR, only : NGID
# endif

	USE MOD_BCMAP, ONLY: 	&
			IOBCN, 			&!
			I_OBC_N
			
	USE MOD_PREC, ONLY: SP	

   
	INTEGER, ALLOCATABLE :: ISBCE(:)     
	INTEGER, ALLOCATABLE :: ISONB(:)  !!NODE MARKER = 0,1,2   
	INTEGER, ALLOCATABLE :: NBE(:,:)  !!INDICES OF ELMNT NEIGHBORS
    INTEGER, ALLOCATABLE :: NTVE(:)         		
    INTEGER, ALLOCATABLE :: NTSN(:)      
	
	INTEGER, ALLOCATABLE :: NV(:,:)         !!NODE NUMBERING FOR ELEMENTS
	INTEGER, ALLOCATABLE :: IEC(:,:)     	!1
    INTEGER, ALLOCATABLE :: IENODE(:,:)  	!2
	
	REAL(SP),ALLOCATABLE :: 	XIJC(:) 		!3
    REAL(SP),ALLOCATABLE :: 	YIJC(:)			!4
	REAL(SP),ALLOCATABLE :: 	DLTXC(:)		!5
	
    REAL(SP),ALLOCATABLE :: 	DLTYC(:)		!6
    REAL(SP),ALLOCATABLE :: 	DLTXYC(:)		!7
	
 	REAL(SP),ALLOCATABLE :: 	SITAC(:) 		!8
	INTEGER, ALLOCATABLE :: ISBC(:) 		!9    
	
    INTEGER, ALLOCATABLE :: NBVE(:,:)		!10
    INTEGER, ALLOCATABLE :: NBVT(:,:)		!11
	INTEGER, ALLOCATABLE :: NBSN(:,:)		!12
	
    INTEGER, ALLOCATABLE :: NIEC(:,:)		!13
    INTEGER, ALLOCATABLE :: NTRG(:)			!14
		  
	REAL(SP),ALLOCATABLE :: 	XIJE(:,:) 		!15
    REAL(SP),ALLOCATABLE :: 	YIJE(:,:) 		!16
    REAL(SP),ALLOCATABLE :: 	DLTXE(:)		!17
    REAL(SP),ALLOCATABLE :: 	DLTYE(:)		!18
    REAL(SP),ALLOCATABLE :: 	DLTXYE(:)		!19 
    REAL(SP),ALLOCATABLE :: 	SITAE(:) 		!20
	
	INTEGER, ALLOCATABLE :: LISBCE_1(:)  	!21   !!LIST OF ELEMENTS WITH ISBCE=1
	INTEGER, ALLOCATABLE :: LISBCE_2(:)  	!22   !!LIST OF ELEMENTS WITH ISBCE=2
	INTEGER, ALLOCATABLE :: LISBCE_3(:)  	!23   !!LIST OF ELEMENTS WITH ISBCE=3
 
	REAL(SP),ALLOCATABLE :: EPOR(:)      		!24   !!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
	
	
	CONTAINS
	
        !Subroutines
		!	Subroutine TRIANGLE_GRID_EDGE()
		!	Subroutine TGE_DEALLOC()

	SUBROUTINE TRIANGLE_GRID_EDGE

!==============================================================================!
!  This program is used to define the non-overlapped, unstructured             !
!  triangular meshes used for flux computations. The mesh could be             !
!  created using the commerical software called "sms8.0" or other              !
!  mesh generation programs. The mesh file generated by sms8.0 can             !
!  be directly used for this subroutine, while the mesh file                   !
!  generated using other programs must be coverted the data format             !
!  to meet the required format used here.                                      !
!==============================================================================!
!     variable list:							       !
!  vx(m)    :: vx(i) = x-coordinate of node i (input from mesh)	               !
!  vy(m)    :: vy(i) = y-coordinate of node i (input from mesh)	               !
!  nv(n,3)  :: nv(i:1-3) = 3 node numbers of element i                         !
!  xc(n)    :: xc(i) = x-coordinate of element i (calculated from vx)          !
!  yc(n)    :: yc(i) = y-coordinate of element i (calculated from vy)          !
!  cor(n)   :: cor(i) = f plane coriolis at element i                          !
!                                                                              !
!  nbe(n,3) :: nbe(i,1->3) = element index of 1->3 neighbors of element i      !
!  isbce(n) :: flag if element is on the boundary, see below for values        !
!  isonb(m) :: flag is node is on the boundary, see below for values           !
!                                                                              !
!  ntve(m)  :: the number of neighboring elements of node m                    !
!  nbve(m,ntve(m)) :: nbve(i,1->ntve(i)) = ntve elements containing node i     !
!  nbvt(m,ntve(m)) :: nbvt(i,j) = the node number of node i in element         ! 
!                     nbve(i,j) (has a value of 1,2,or 3)                      ! 
!                                                                              !
!   ne       ::   number of unique element edges                               !
!  iec(ne,2) :: nbe(i,1->2) cell number of cell(s) connected to edge i         !
!   isbc(ne) :: flag marking edge property                                     !
!              isbc(i) = 0:  element edge i is in the interior                 !
!              isbc(i) = 1:  element edge i is on the boundary                 !
!ienode(ne,2):: ienode(i,1->2) node numbers at each end of element edge i      ! 
!  xijc(ne)  :: xijc(i) = x-coordinate of mid point of element edge i          ! 
!  yijc(ne)  :: yijc(i) = y-coordinate of mid point of element edge i          ! 
!  dltxyc(ne):: dltxyc(i) = length of element edge i                           !
!  dltxc(ne) :: dltxc(i) = deltax (x-projection) of element edge i             !
!  dltyc(ne) :: dltyc(i) = deltay (y-projection) of element edge i             !
!  sitac(ne) :: sitac(i) =  arctg(dltyc,dltxc) (angle of inclination of edge)  !
!                                                                              !
!==============================================================================!
!     classification of the triangles nodes, and edges                         !
!                                                                              !
!     isonb(i)=0:  node in the interior computational domain                   !
!     isonb(i)=1:  node on the solid boundary                                  !
!     isonb(i)=2:  node on the open boundary                                   !
!                                                                              !
!     isbce(i)=0:  element in the interior computational domain                !
!     isbce(i)=1:  element on the solid boundary                               !
!     isbce(i)=2:  element on the open boundary                                !
!     isbce(i)=3:  element with 2 solid boundary edges                         !
!                                                                              !
!      isbc(i)=0:  element edge in interior                                    !
!      isbc(i)=1:  element edge on boundary                                    !
!==============================================================================!


!==============================================================================|
!   FIND NEIGHBORING ELEMENTS, MARK BOUNDARY NODES AND ELEMENTS                |
!									       |
!   NBE(N,3) :  NBE(I,J) = IDENTITY OF NBOR ELMNT TO TRI I ON EDGE J           |
!   IBCE(N)  :  DESCRIBED IN SUBROUTINE HEADING			               |	
!   ISONB(M):  DESCRIBED IN SUBROUTINE HEADING			               |	
!==============================================================================|
   
   
   IMPLICIT NONE
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TEMP,TEMP2,NB_TMP,ISET
   INTEGER I,J,DIF1,DIF2,DIF3,II,JJ,NTMP,NCNT,INEY,NFLAG
   INTEGER ITMP1,ITMP2,ITMP3,JN,JJB,IBCETMP,NCTMP,NCETMP,NPT
   INTEGER, ALLOCATABLE :: CELLS(:,:),CELLCNT(:),NBET(:,:)
   REAL(SP)  DTMP
   INTEGER N1,N2,N3,J1,J2,J3

!----------------------------REPORT--------------------------------------------!
   IF(MSR)WRITE(*,*  )'!'
   IF(MSR)WRITE(*,*)'!           SETTING UP TRIS/ELEMENTS/EDGES/CVS          '
   IF(MSR)WRITE(*,*  )'!'
!----------------------------INITIALIZE----------------------------------------!
   
   ALLOCATE(ISBCE(0:NTLOC))         ;ISBCE = 0
   ALLOCATE(ISONB(0:MTLOC))         ;ISONB = 0  !!NODE MARKER = 0,1,2   
   ALLOCATE(NBE(0:NTLOC,3))         ;NBE   = 0  !!INDICES OF ELMNTNEIGHBORS

!
!----DETERMINE NBE(i=1:n,j=1:3): INDEX OF 1 to 3 NEIGHBORING ELEMENTS----------!
!
   ALLOCATE(NBET(NTLOC,3)) ; NBET = 0
   ALLOCATE(CELLS(MTLOC,50)) ; CELLS = 0
   ALLOCATE(CELLCNT(MTLOC))  ; CELLCNT = 0
   DO I=1,NTLOC
     N1 = NV(I,1) ; CELLCNT(N1) = CELLCNT(N1)+1
     N2 = NV(I,2) ; CELLCNT(N2) = CELLCNT(N2)+1
     N3 = NV(I,3) ; CELLCNT(N3) = CELLCNT(N3)+1
     CELLS(NV(I,1),CELLCNT(N1)) = I
     CELLS(NV(I,2),CELLCNT(N2)) = I
     CELLS(NV(I,3),CELLCNT(N3)) = I
   ENDDO
   if(maxval(cellcnt) > 50)write(*,*)'bad',maxval(cellcnt)
   DO I=1,NTLOC
     N1 = NV(I,1)
     N2 = NV(I,2)
     N3 = NV(I,3)
     DO J1 = 1,CELLCNT(N1) 
     DO J2 = 1,CELLCNT(N2) 
       IF((CELLS(N1,J1) == CELLS(N2,J2)).AND. CELLS(N1,J1) /= I)NBE(I,3) = CELLS(N1,J1)
     ENDDO
     ENDDO
     DO J2 = 1,CELLCNT(N2) 
     DO J3 = 1,CELLCNT(N3) 
       IF((CELLS(N2,J2) == CELLS(N3,J3)).AND. CELLS(N2,J2) /= I)NBE(I,1) = CELLS(N2,J2)
     ENDDO
     ENDDO
     DO J1 = 1,CELLCNT(N1) 
     DO J3 = 1,CELLCNT(N3) 
       IF((CELLS(N1,J1) == CELLS(N3,J3)).AND. CELLS(N1,J1) /= I)NBE(I,2) = CELLS(N3,J3)
     ENDDO
     ENDDO
   ENDDO
   DEALLOCATE(CELLS,CELLCNT,NBET)  !LB added NBET to the deallocation statement
   IF(MSR)WRITE(*,*)  '!  NEIGHBOR FINDING      :    COMPLETE'
!
!--ENSURE ALL ELEMENTS HAVE AT LEAST ONE NEIGHBOR------------------------------!
!
   NFLAG = 0
   DO I=1,NTLOC
     IF(SUM(NBE(I,1:3))==0)THEN 
       NFLAG = 1
       WRITE(*,*)'ELEMENT ',I,' AT ',XC(I),YC(I),' HAS NO NEIGHBORS'
       CALL PSTOP
     ENDIF
   ENDDO
   IF(NFLAG == 1) CALL PSTOP
     
!
!----IF ELEMENT ON BOUNDARY SET ISBCE(I)=1 AND ISONB(J)=1 FOR BOUNDARY NODES J-!
!

   DO I=1,NTLOC 
     IF(MIN(NBE(I,1),NBE(I,2),NBE(I,3))==0)THEN    !!ELEMENT ON BOUNDARY
       ISBCE(I) = 1
       IF(NBE(I,1) == 0)THEN 
         ISONB(NV(I,2)) = 1 ; ISONB(NV(I,3)) = 1
       ENDIF
       IF(NBE(I,2) ==0) THEN
         ISONB(NV(I,1)) = 1 ; ISONB(NV(I,3)) = 1
       ENDIF
       IF(NBE(I,3) ==0) THEN
         ISONB(NV(I,1)) = 1 ; ISONB(NV(I,2)) = 1
       ENDIF
     ENDIF
   ENDDO
   IF(MSR)WRITE(*,*)  '!  ISONB SETTING         :    COMPLETE'
        
!==============================================================================|
!             DEFINE NTVE, NBVE, NBVT                                          !
!                                                                              !
! ntve(1:m):           total number of the surrounding triangles               !
!                      connected to the given node                             !                                        
! nbve(1:m, 1:ntve+1): the identification number of surrounding                !
!                      triangles with a common node (counted clockwise)        !
! nbvt(1:m,ntve(1:m)): the idenfication number of a given node over            !
!                      each individual surrounding triangle(counted            !
!                      clockwise)                                              !
! ntsn(1:m):           total number of surrounding nodes                       !
! nbsn(1:m, ntsn):     the identification number of surrounding nodes          !
!                      (counted clockwise)                                     !
! nbse(1:m,2*ntsn):    the identification number of control volume s           !
!                      edges between two neighbor nodes                        !
!==============================================================================|

!
!----DETERMINE MAX NUMBER OF SURROUNDING ELEMENTS------------------------------!
!
   MX_NBR_ELEM = 0
   DO I=1,MLOC
     NCNT = 0
     DO J=1,NTLOC
!       IF( (NV(J,1)-I)*(NV(J,2)-I)*(NV(J,3)-I)==0) NCNT = NCNT + 1
!       IF( (NV(J,1)-I) ==0 .OR. (NV(J,2)-I) ==0 .OR. (NV(J,3)-I)==0) NCNT = NCNT + 1
       IF( FLOAT(NV(J,1)-I)*FLOAT(NV(J,2)-I)*FLOAT(NV(J,3)-I) == 0.0_SP) &
         NCNT = NCNT + 1
     ENDDO
     MX_NBR_ELEM = MAX(MX_NBR_ELEM,NCNT)
   ENDDO
!   WRITE(*,*) 'MAXIMUM NUMBER OF NEIGHBOR ELEMENTS',MX_NBR_ELEM

!
!----ALLOCATE ARRAYS BASED ON MX_NBR_ELEM--------------------------------------!
! 
   
   ALLOCATE(NBVE(MLOC,MX_NBR_ELEM+1))
   ALLOCATE(NBVT(MLOC,MX_NBR_ELEM+1))
   ALLOCATE(NBSN(MLOC,MX_NBR_ELEM+3)) !!MHB
!
!--DETERMINE NUMBER OF SURROUNDING ELEMENTS FOR NODE I = NTVE(I)---------------!
!--DETERMINE NBVE - INDICES OF NEIGHBORING ELEMENTS OF NODE I------------------!
!--DETERMINE NBVT - INDEX (1,2, or 3) OF NODE I IN NEIGHBORING ELEMENT---------!
!

   
   ALLOCATE(NTVE(0:MTLOC))       ;NTVE = 0       
   DO I=1,MLOC
     NCNT=0
     DO J=1,NTLOC
!       IF( (NV(J,1)-I) == 0 .OR.  (NV(J,2)-I) == 0 .OR. (NV(J,3)-I) == 0)THEN 
        IF (FLOAT(NV(J,1)-I)*FLOAT(NV(J,2)-I)*FLOAT(NV(J,3)-I) == 0.0_SP)THEN
         NCNT = NCNT+1
         NBVE(I,NCNT)=J
         IF((NV(J,1)-I) == 0) NBVT(I,NCNT)=1
         IF((NV(J,2)-I) == 0) NBVT(I,NCNT)=2
         IF((NV(J,3)-I) == 0) NBVT(I,NCNT)=3
       ENDIF
     ENDDO
     NTVE(I)=NCNT
   ENDDO

!
!--Reorder Order Elements Surrounding a Node to Go in a Cyclical Procession----!
!--Determine NTSN  = Number of Nodes Surrounding a Node (+1)-------------------!
!--Determine NBSN  = Node Numbers of Nodes Surrounding a Node------------------!
!
  
   ALLOCATE(NTSN(MTLOC))         ;NTSN     = 0   
   
   ALLOCATE(NB_TMP(MLOC,MX_NBR_ELEM+1))
   
   DO I=1,MLOC
     IF(ISONB(I) == 0) THEN
       NB_TMP(1,1)=NBVE(I,1)
       NB_TMP(1,2)=NBVT(I,1)
       DO J=2,NTVE(I)+1
         II=NB_TMP(J-1,1)
         JJ=NB_TMP(J-1,2)
         NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
         JJ=NB_TMP(J,1)
         IF((NV(JJ,1)-I) == 0) NB_TMP(J,2)=1
         IF((NV(JJ,2)-I) == 0) NB_TMP(J,2)=2
         IF((NV(JJ,3)-I) == 0) NB_TMP(J,2)=3
       ENDDO

       DO J=2,NTVE(I)+1
         NBVE(I,J)=NB_TMP(J,1)
       ENDDO

       DO J=2,NTVE(I)+1
         NBVT(I,J)=NB_TMP(J,2)
       ENDDO

       NTMP=NTVE(I)+1
       IF(NBVE(I,1) /= NBVE(I,NTMP)) THEN
          WRITE(*,*) MYID,I,'NBVE(I) NOT CORRECT!!'
          CALL PSTOP
       ENDIF
       IF(NBVT(I,1) /= NBVT(I,NTMP)) THEN
          WRITE(*,*) I,'NBVT(I) NOT CORRECT!!'
          CALL PSTOP
       ENDIF

       NTSN(I)=NTVE(I)

       DO J=1,NTSN(I)
         II=NBVE(I,J)
         JJ=NBVT(I,J)
         NBSN(I,J)=NV(II,JJ+1-INT((JJ+1)/4)*3)
       ENDDO

       NTSN(I)=NTSN(I)+1
       NBSN(I,NTSN(I))=NBSN(I,1)

     ELSE 
           JJB=0

       DO J=1,NTVE(I)
         JJ=NBVT(I,J)
         IF(NBE(NBVE(I,J),JJ+2-INT((JJ+2)/4)*3) == 0) THEN
           JJB=JJB+1
           NB_TMP(JJB,1)=NBVE(I,J)
           NB_TMP(JJB,2)=NBVT(I,J)
         ENDIF
       ENDDO

       IF(JJB /= 1) THEN
         WRITE(*,*), 'ERROR IN ISONB !,I,J', I,J
         PAUSE
       ENDIF

       DO J=2,NTVE(I)
         II=NB_TMP(J-1,1)
         JJ=NB_TMP(J-1,2)
         NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
         JJ=NB_TMP(J,1)
         IF((NV(JJ,1)-I) == 0) NB_TMP(J,2)=1
         IF((NV(JJ,2)-I) == 0) NB_TMP(J,2)=2
         IF((NV(JJ,3)-I) == 0) NB_TMP(J,2)=3
       ENDDO

       DO J=1,NTVE(I)
         NBVE(I,J)=NB_TMP(J,1)
         NBVT(I,J)=NB_TMP(J,2)
       ENDDO

       NBVE(I,NTVE(I)+1)=0
       NTSN(I)=NTVE(I)+1
       NBSN(I,1)=I

       DO J=1,NTSN(I)-1
         II=NBVE(I,J)
         JJ=NBVT(I,J)
         NBSN(I,J+1)=NV(II,JJ+1-INT((JJ+1)/4)*3)
       ENDDO

       J=NTSN(I)
       II=NBVE(I,J-1)
       JJ=NBVT(I,J-1)
       NBSN(I,J+1)=NV(II,JJ+2-INT((JJ+2)/4)*3)
       NTSN(I)=NTSN(I)+2
       NBSN(I,NTSN(I))=I
     ENDIF
   ENDDO
   DEALLOCATE(NB_TMP)
   IF(MX_NBR_ELEM+3 -MAXVAL(NTSN) < 0)THEN
      WRITE(*,*)'CHECK NTSN/NBSN',MAXVAL(NTSN),MX_NBR_ELEM+3
   !   CALL PSTOP  ! B Clark commented to run ideal test case, dont need NTSN
   ENDIF
   IF(MSR)WRITE(*,*)  '!  NBVE/NBVT             :    COMPLETE'
     

!==============================================================================!
!  Define the parameters of each triangular edge                               !
!                                                                              !
!  ne           :    number of unique element edges                            !
!  iec(1:ne,1:2):    counting number identifying two connected cells           !
!  isbc(1:ne):       0: triangle s edge in the interior                        !
!                    1: triangle s edge on the boundary                        !
!  ienode(1:ne,1:2): the identification number of two end points of a          !
!                    edge                                                      !
!  xijc(1:ne):       the x-coordinate location of the middle points            !
!                    of a edge                                                 !
!  yijc(1:ne):       the y-coordinate location of the middle points            !
!                    of a edge                                                 !
!  dltxyc(1:ne):     length of the edge                                        !
!  dltxc(1:ne):      vx(ienode(i,2))-vx(idnode(i,1))                           !
!  dltyc(1:ne):      vy(ienode(i,2))-vy(idnode(i,1))                           !
!  sitac(1:ne):      arctg(dltyc,dltxc)                                        !
!==============================================================================!

   ALLOCATE(ISET(NTLOC,3),TEMP((NTLOC)*3,2),TEMP2((NTLOC)*3,2))
   ISET = 0
   NE = 0
   TEMP = 0 
   TEMP2 = 0
   DO I=1,NTLOC
     DO J=1,3
       IF(ISET(I,J) == 0)THEN
         NE   = NE + 1
         INEY = NBE(I,J)
         ISET(I,J) = 1
         DO JN=1,3
           IF(I == NBE(INEY,JN)) ISET(INEY,JN) = 1
         ENDDO
         TEMP(NE,1) = I ; TEMP(NE,2) = INEY
         TEMP2(NE,1) = NV(I,J+1-INT((J+1)/4)*3)
         TEMP2(NE,2) = NV(I,J+2-INT((J+2)/4)*3)
       ENDIF
     ENDDO
   ENDDO
   DEALLOCATE(ISET)
!
!--ALLOCATE ARRAYS REQUIRING NUMBER OF EDGES-----------------------------------!
!
   ALLOCATE(IEC(NE,2))
   ALLOCATE(IENODE(NE,2))
   ALLOCATE(XIJC(NE))
   ALLOCATE(YIJC(NE))
   ALLOCATE(DLTXC(NE))
   ALLOCATE(DLTYC(NE))
   ALLOCATE(DLTXYC(NE))
   ALLOCATE(SITAC(NE))
   ALLOCATE(ISBC(NE))


   IEC(:,1) = TEMP(1:NE,1)
   IEC(:,2) = TEMP(1:NE,2)
   IENODE(:,1) = TEMP2(1:NE,1)
   IENODE(:,2) = TEMP2(1:NE,2)


   DEALLOCATE(TEMP,TEMP2)
         
         
!
!------MARK ELEMENT EDGES THAT ARE ON THE BOUNDARY-----------------------------!
!
   ISBC = 0
   DO I=1,NE
     IF((IEC(I,1) == 0) .OR. (IEC(I,2) == 0)) ISBC(I) = 1 
   ENDDO

!
!------CALCULATE ELEMENT EDGE METRICS------------------------------------------!
!
   DO I=1,NE
     DLTXC(I) =  VX(IENODE(I,2))-VX(IENODE(I,1))
     DLTYC(I) =  VY(IENODE(I,2))-VY(IENODE(I,1))
     XIJC(I)  = (VX(IENODE(I,1))+VX(IENODE(I,2)))/2.0
     YIJC(I)  = (VY(IENODE(I,1))+VY(IENODE(I,2)))/2.0
     DLTXYC(I)= SQRT(DLTXC(I)**2+DLTYC(I)**2)
     SITAC(I) = ATAN2(DLTYC(I),DLTXC(I))
   ENDDO

   IF(MSR)WRITE(*,*)  '!  EDGE SETUP            :    COMPLETE'


!==============================================================================!
!  read triangular mesh parameters on open boundary :                          !
!  iobce:   number of open boundary cells.                                     !
!  isbcn:   number of open boundary nodes.                                     !
!  i_obc_e: counter number of open boundary cells   !never used                !
!  i_obc_n: counter number of open boundary nodes                              !
!==============================================================================!

!
!----TRAVERSE  BOUNDARY NODE NUMBERS AND SET ISONB(NODE)=2---------------------!
!
   DO I=1,IOBCN
     ISONB(I_OBC_N(I))=2
   ENDDO

!
!----DETERMINE IF ELEMENT IS ON OPEN BOUNDARY (CONTAINS EDGE ON OPEN BOUNDARY)-!
!
   IBCETMP=0
   DO I=1,NLOC
     ITMP1=ISONB(NV(I,1))
     ITMP2=ISONB(NV(I,2))
     ITMP3=ISONB(NV(I,3))

     IF(SUM(ISONB(NV(I,1:3))) == 4) THEN
       ISBCE(I)=2
       IBCETMP =IBCETMP+1
     ELSE IF(SUM(ISONB(NV(I,1:3))) > 4) THEN
#  if defined (MULTIPROCESSOR)
       WRITE(*,*)'SORRY, THE BOUNDARY CELL',NGID(I),'IS NOT GOOD FOR MODEL.'
#  else
       WRITE(*,*)'SORRY, THE BOUNDARY CELL',I,'IS NOT GOOD FOR MODEL.'
#  endif       
       WRITE(*,*)'IT HAS EITHER TWO SIDES OF OPEN BOUNDARY OR ONE OPEN BOUNDARY'
       WRITE(*,*),'AND ONE SOLID BOUNDARY. PLEASE CHECK AND MODIFIED IT.'
       WRITE(*,*),'THIS MESSAGE IS IN SUBROUTINE TRIANGLE_GRID_EDGE (MOD_TGE.F)'
       WRITE(*,*),'STOP RUNNING...'
       CALL PSTOP
     ENDIF
   ENDDO

   DO I=1,NTLOC
     IF((NBE(I,1)+NBE(I,2)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,1)+NBE(I,2) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,2)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,1)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
   ENDDO

!==============================================================================!
!  xije(1:nc,1:2):  the x coordinate locations of starting and ending          !
!                   points of the control volumes edges                        !
!  yije(1:nc,1:2):  the y coordinate locations of starting and ending          !
!                   points of the control volumes edges                        !
!  niec(1:nc,1:2):  the counting number of left and right nodes                !
!                   conected to this control volumes edge from                 !
!                   starting point to ending point                             !
!  dltxe(1:nc):     the x distance of individual edges                         !
!  dltye(1:nc)      the y distance of individual edges                         !
!  dltxye(1:nc):    the length of individual edges                             !
!  ntrg(1:nc)  :    element associated with this control volume edge           !
!==============================================================================!

   NCT = NTLOC*3
   ALLOCATE(XIJE(NCT,2))         ;XIJE     = ZERO
   ALLOCATE(YIJE(NCT,2))         ;YIJE     = ZERO
   ALLOCATE(NIEC(NCT,2))         ;NIEC     = 0
   ALLOCATE(NTRG(NCT))           ;NTRG     = 0
   ALLOCATE(DLTXE(NCT))          ;DLTXE    = ZERO
   ALLOCATE(DLTYE(NCT))          ;DLTYE    = ZERO
   ALLOCATE(DLTXYE(NCT))         ;DLTXYE   = ZERO
   ALLOCATE(SITAE(NCT))          ;SITAE    = ZERO
   NCTMP  = 0
   NCETMP = 0 

   DO I=1,NE
     IF(ISBC(I) == 0) THEN
       IF(IEC(I,1) <= NLOC)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*NLOC)
       ENDIF
       XIJE(NPT,1) = XC(IEC(I,1))
       YIJE(NPT,1) = YC(IEC(I,1))
       XIJE(NPT,2) = XIJC(I)
       YIJE(NPT,2) = YIJC(I)
       NIEC(NPT,1) = IENODE(I,1)
       NIEC(NPT,2) = IENODE(I,2)
       NTRG(NPT)   = IEC(I,1)
       DLTXE(NPT)  = XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)  = YIJE(NPT,2)-YIJE(NPT,1)
       DTMP        = DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT) = SQRT(DTMP)
       SITAE(NPT)  = ATAN2(DLTYE(NPT),DLTXE(NPT))

       IF(IEC(I,2) <= NLOC)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*NLOC)
       ENDIF
       XIJE(NPT,1)=XC(IEC(I,2))
       YIJE(NPT,1)=YC(IEC(I,2))
       XIJE(NPT,2)=XIJC(I)
       YIJE(NPT,2)=YIJC(I)
       NIEC(NPT,1)=IENODE(I,2)
       NIEC(NPT,2)=IENODE(I,1)
       NTRG(NPT)=IEC(I,2)
       DLTXE(NPT)=XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)=YIJE(NPT,2)-YIJE(NPT,1)
       DTMP=DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT)=SQRT(DTMP)
       SITAE(NPT)=ATAN2(DLTYE(NPT),DLTXE(NPT))
     ELSE IF(ISBC(I) == 1) THEN
       IF(IEC(I,1) <= NLOC)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*NLOC)
       ENDIF
       IF(IEC(I,1) == 0) THEN
         WRITE(*,*)I,'IEC(I,1)===0'
         CALL PSTOP
       ENDIF
       XIJE(NPT,1)=XC(IEC(I,1))
       YIJE(NPT,1)=YC(IEC(I,1))
       XIJE(NPT,2)=XIJC(I)
       YIJE(NPT,2)=YIJC(I)
       NIEC(NPT,1)=IENODE(I,1)
       NIEC(NPT,2)=IENODE(I,2)
       NTRG(NPT)=IEC(I,1)
       DLTXE(NPT)=XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)=YIJE(NPT,2)-YIJE(NPT,1)
       DTMP=DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT)=SQRT(DTMP)
       SITAE(NPT)=ATAN2(DLTYE(NPT),DLTXE(NPT))
     ELSE
       WRITE(*,*) 'ISBC(I) NOT CORRECT, I==',I
       CALL PSTOP
     ENDIF
   ENDDO

   NCV_I  = NCTMP
   NCV    = NCETMP+NCTMP
      
   IF(NCV /= 3*(NTLOC)) THEN
      WRITE(*,*)'NCV IS NOT CORRECT, PLEASE CHECK THE SETUP'
      CALL PSTOP
   ENDIF
   IF(NCV_I /= 3*NLOC) THEN
      WRITE(*,*)'NCV_I IS NOT CORRECT, PLEASE CHECK THE SETUP'
      CALL PSTOP
   ENDIF

   DO I=1,NCV_I
     IF(NIEC(I,1) > MLOC .OR. NIEC(I,2) > MLOC)THEN
       WRITE(*,*)'PROBLEMS:',niec(i,1),niec(i,2),mloc
       CALL PSTOP
     ENDIF
   ENDDO
   
!==============================================================================!
!  nisbce_1/nisbce_2/nisbce_3:  number of elements with isbce of 1,2,3         !
!  lisbce_1/lisbce_2/lisbce_3:  list of elements with isbce of 1,2,3           !
!  epor                      :  element porosity (=0 if isbce = 2)             !
!==============================================================================!

!
!--COUNT NUMBER OF ELEMENTS OF EACH TYPE (ISBCE=1,2,3)-------------------------!
!
   NISBCE_1 = 0 ; NISBCE_2 = 0 ; NISBCE_3 = 0
   DO I=1,NLOC 
     IF(ISBCE(I) == 1) NISBCE_1 = NISBCE_1 + 1
     IF(ISBCE(I) == 2) NISBCE_2 = NISBCE_2 + 1
     IF(ISBCE(I) == 3) NISBCE_3 = NISBCE_3 + 1
   ENDDO

!
!--ALLOCATE ELEMENT TYPE ARRAYS LISBCE_1,LISBCE_2,LISBCE_3---------------------!
!
   IF(NISBCE_1 > 0)THEN
     ALLOCATE( LISBCE_1(NISBCE_1) )
   ELSE
!     WRITE(IPT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=1'
   ENDIF
     
   IF(NISBCE_2 > 0)THEN
     ALLOCATE( LISBCE_2(NISBCE_2) )
   ELSE
!     WRITE(IPT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=2'
   ENDIF
     
   IF(NISBCE_3 > 0)THEN
     ALLOCATE( LISBCE_3(NISBCE_3) )
   ELSE
!     WRITE(IPT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=3'
   ENDIF

!
!--LOAD ELEMENT TYPE ARRAYS LISBCE_1,LISBCE_2,LISBCE_3--------------------------!
!
   NISBCE_1 = 0 ; NISBCE_2 = 0 ; NISBCE_3 = 0
   DO I=1,NLOC 
     IF(ISBCE(I) == 1) THEN
       NISBCE_1 = NISBCE_1 + 1
       LISBCE_1(NISBCE_1) = I
     ENDIF
     IF(ISBCE(I) == 2) THEN
       NISBCE_2 = NISBCE_2 + 1
       LISBCE_2(NISBCE_2) = I
     ENDIF
     IF(ISBCE(I) == 3) THEN
       NISBCE_3 = NISBCE_3 + 1
       LISBCE_3(NISBCE_3) = I
     ENDIF
   ENDDO

!
!--SET ELEMENT POROSITY---------------------------------------------------------!
!
   ALLOCATE(EPOR(0:NTLOC)) ; EPOR = 1.0_SP
   DO I=1,NLOC
     IF(ISBCE(I) == 2) EPOR(I) = 0.0_SP
   ENDDO
     
   IF(MSR)WRITE(*,*)  '!  NISBCE/LISBCE/EPOR    :    COMPLETE'

   IF(MSR)WRITE(*,*)  '!  TRIS/EDGES/CVOLS      :    COMPLETE'

   RETURN
   END SUBROUTINE TRIANGLE_GRID_EDGE
   
   SUBROUTINE TGE_ALLOC
		!WLong moved here from pdomdec.F
		ALLOCATE(NV(0:NTLOC,4));    NV = 0  !!NODE NUMBERING FOR ELEMENTS
				
   END SUBROUTINE TGE_ALLOC
   
   SUBROUTINE TGE_DEALLOC
   
	!Deallocate variables defined in MOD_TGE
   
   	 IF (ALLOCATED(ISBCE)) 	DEALLOCATE(ISBCE)
     IF (ALLOCATED(ISONB)) 	DEALLOCATE(ISONB)
     IF (ALLOCATED(NBE)) 	DEALLOCATE(NBE) 
	 IF (ALLOCATED(NTVE)) 	DEALLOCATE(NTVE)
	 IF (ALLOCATED(NTSN)) 	DEALLOCATE(NTSN)
	 
	 IF	(ALLOCATED(NV)) 	DEALLOCATE(NV)
	 IF (ALLOCATED(IEC)) 	DEALLOCATE(IEC)
	 IF (ALLOCATED(IENODE)) DEALLOCATE(IENODE)
	 
	 IF (ALLOCATED(XIJC)) 	DEALLOCATE(XIJC)
	 IF (ALLOCATED(YIJC)) 	DEALLOCATE(YIJC)
	 IF (ALLOCATED(DLTXC)) 	DEALLOCATE(DLTXC)
	 IF (ALLOCATED(DLTYC)) 	DEALLOCATE(DLTYC)
	 IF (ALLOCATED(DLTXYC)) DEALLOCATE(DLTXYC)
	 IF (ALLOCATED(SITAC)) 	DEALLOCATE(SITAC)
	 IF (ALLOCATED(ISBC)) 	DEALLOCATE(ISBC)

	 IF (ALLOCATED(NBVE)) 	DEALLOCATE(NBVE)
	 IF (ALLOCATED(NBVT)) 	DEALLOCATE(NBVT)
	 IF (ALLOCATED(NBSN)) 	DEALLOCATE(NBSN)

	 IF (ALLOCATED(NIEC)) 	DEALLOCATE(NIEC)
     IF (ALLOCATED(NTRG)) 	DEALLOCATE(NTRG)
	 
	 IF (ALLOCATED(XIJE)) 	DEALLOCATE(XIJE)
     IF (ALLOCATED(YIJE)) 	DEALLOCATE(YIJE)
     IF (ALLOCATED(DLTXE)) 	DEALLOCATE(DLTXE)
     IF (ALLOCATED(DLTYE)) 	DEALLOCATE(DLTYE)
     IF (ALLOCATED(DLTXYE)) DEALLOCATE(DLTXYE)
     IF (ALLOCATED(SITAE)) 	DEALLOCATE(SITAE)
	 IF (ALLOCATED(EPOR)) 	DEALLOCATE(EPOR)	 
     
     !LB: these ones were missing
     IF (ALLOCATED(LISBCE_1))   DEALLOCATE(LISBCE_1)
     IF (ALLOCATED(LISBCE_2))   DEALLOCATE(LISBCE_2)
     IF (ALLOCATED(LISBCE_3))   DEALLOCATE(LISBCE_3)

     

   RETURN

   END SUBROUTINE TGE_DEALLOC

!==============================================================================!

END MODULE MOD_TGE
