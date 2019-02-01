MODULE MOD_BCMAP

	USE MOD_PREC, ONLY: SP
	USE MOD_TYPES, ONLY :  BC
	
   !WLong moved the following here from mod_obcs.F
   INTEGER               :: IOBCN_GL         !!GLOBAL NUMBER OF OPEN BOUNDARY NODES   
   INTEGER               :: IOBCN            !!LOCAL NUMBER OF OPEN BOUNDARY NODES    
   INTEGER,  ALLOCATABLE :: I_OBC_GL(:)      !!GLOBAL ID OF OPEN BOUNDARY NODES		  
   INTEGER,  ALLOCATABLE :: I_OBC_N(:)       !!OPEN BOUNDARY NODE LIST				  
   
	!WLong moved here from mod_bcs.F
   REAL(SP), ALLOCATABLE     :: WQOBC(:,:,:,:)   !!WATER QUALITY DATA AT BOUNDARY	
   TYPE(BC)              :: NUT_TM           !!TIME MAPPING FOR NUTRIENT OBC  !should be 

   CONTAINS

		!Subroutine BCMAP()
		!Subroutine BCMAP_DEALLOC

!--------------------------------------------
   SUBROUTINE  BCMAP
!--------------------------------------------
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
             MGL!,            &  !
			!NGL,            &  !
			!OBCGL,          &  !
			!NOBTY  
   
   USE MOD_LIMS, ONLY: KBM1, NPROCS
	!Wen Long took CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY : 		&
			SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			,MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			,PAR        	&           !!TRUE IF MULTIPROCESSOR RUN
			,CASENAME ! 	&   		!!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
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
			
   USE MOD_PREC, ONLY: SP
     
   USE MOD_UTILS, ONLY: &
		!PERROR,   		&!
		!GET_TIMESAMP,	&!
		FOPEN!,			&!
		!CRAY_SYSTEM_CALL
		
# if defined (MULTIPROCESSOR)
		USE MOD_PAR, ONLY:   &	  
			!EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			!ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			NLID!,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			!ELID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			!NLID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,	&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			!EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			!NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
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
			!NHN,  		&	!			!!NUMBER OF HALO NODES
			!HN_LST!, 	&	!			!!LIST OF HALO NODES 
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

# endif
	IMPLICIT NONE	
	

#if defined (MULTIPROCESSOR)
   include "mpif.h"
#endif
   INTEGER              :: I,I1,I2,NCNT,IERR,J,JCON,JJ,JT,K
   INTEGER, ALLOCATABLE :: TEMP1(:),TEMP2(:),TEMP3(:),TEMP4(:),&
                           TEMP5(:),TEMP6(:),TEMP7(:),ITEMP(:)
   REAL(SP), ALLOCATABLE :: TEMP9(:,:,:,:)  ! KURT GLAESEMANN 3-MAR-2010
   REAL(SP) :: FTEMP1
   REAL(SP), ALLOCATABLE  :: WQOBC_GL(:,:,:,:)   !!NUT Series at OB sigma layer

   INTEGER :: INOBC,NTMP 			   
   CHARACTER(LEN=100) ::ISTR

   INOBC = 13
!------------------------------------------------------------------------------!

!==============================================================================|
!   OPEN BOUNDARY CONDITION NODES                                              |
!==============================================================================|

!----------------------------REPORT--------------------------------------------!

   IF(MSR)WRITE(*,*  )'!'
   IF(MSR)WRITE(*,*)'!           SETTING UP OPEN BOUNDARY NODES  '
   IF(MSR)WRITE(*,*  )'!'

   IOBCN = 0
!   IBCN  = 0
 
   ISTR = "./inputs/"//trim(casename)
    CALL FOPEN(INOBC, TRIM(ISTR)//'_obc_wq.dat',"cfr")
    REWIND(INOBC)

      READ(INOBC,*)IOBCN_GL
   
   IF(IOBCN_GL > 0)THEN

!------------Read in Open Boundary Nodes and Temperature/Salinity Conditions---!

!
!----Input Number of Data Times for T/S at Open Boundary Every Sigma Layers---!
!
     ALLOCATE(I_OBC_GL(IOBCN_GL))

     DO I=1,IOBCN_GL
       READ(INOBC,*) I1,I_OBC_GL(I)
     ENDDO
      NCNT = 0
      DO WHILE(.TRUE.)
        READ(INOBC,*,END=10)FTEMP1
        DO JJ=1,NCP
           DO J=1,IOBCN_GL
              READ(INOBC,*)
           ENDDO
        ENDDO
        NCNT = NCNT + 1
      ENDDO
 10   CONTINUE
      REWIND(INOBC)

      IF(NCNT == 0) write(*,*)'NO WQ DATA at BND'

      NUT_TM%NTIMES = NCNT

!    READ nutrient open boundary values

     ALLOCATE(NUT_TM%TIMES(NUT_TM%NTIMES))         ; NUT_TM%TIMES = 0.
     ALLOCATE(WQOBC_GL(IOBCN_GL,KBM1,NUT_TM%NTIMES,NCP)) ; WQOBC_GL = 0.

     READ(INOBC,*)
     DO I=1,IOBCN_GL
       READ(INOBC,*)
     ENDDO
      DO JT=1,NUT_TM%NTIMES
         READ(INOBC,*) NUT_TM%TIMES(JT)               !time(hours)
	!	 write(*,*)'read the time',JT
         DO JJ=1,NCP
            DO I=1,IOBCN_GL            
               READ(INOBC,*) ntmp,(WQOBC_GL(I,K,JT,JJ),K=1,KBM1) !temp in OB sigma layer
			          
		!	write(*,*)'Variable number at obc read'
		!	write(*,*) JJ, i
            ENDDO

          ENDDO
!          NUT_TM%TIMES = 24.0*(NUT_TM%TIMES) !shift to model hours
      ENDDO
	  write(*,*)'Read in all OBC data'
      CLOSE(INOBC)

!----------------------Make Sure It Is In Global Domain------------------------!

     DO I=1,IOBCN_GL
       IF((I_OBC_GL(I) > MGL))THEN
         WRITE(*,*)'==================ERROR=================================='
         WRITE(*,*)'OPEN BOUNDARY NODE NUMBER',I,'IS NOT IN THE'
         WRITE(*,*)'GLOBAL DOMAIN'
         WRITE(*,*)'CHECK INPUT FILE AND ENSURE OPEN BOUNDARY NODES <= ',MGL
         WRITE(*,*)'========================================================='
         CALL PSTOP
       ENDIF
     ENDDO

!----------Shift Open Boundary Node List,Type,Salt,and Temp to Local-----------!

     IF(SERIAL)THEN

        IOBCN    = IOBCN_GL
        ALLOCATE(I_OBC_N(IOBCN))
        I_OBC_N = I_OBC_GL

        ALLOCATE(WQOBC(IOBCN,KBM1,NUT_TM%NTIMES,NCP)) ; WQOBC=0.0
        WQOBC(:,:,:,:) = WQOBC_GL(:,:,:,:) ! KURT GLAESEMANN 19-MAR-2010

     ENDIF

#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       ALLOCATE(TEMP1(IOBCN_GL))
       ALLOCATE(TEMP9(IOBCN_GL,KBM1,NUT_TM%NTIMES,NCP)) ; TEMP9 = 0 ! KURT GLAESEMANN 3-MAR-2010
       NCNT = 0
       !!SET UP LOCAL OPEN BOUNDARY NODES
       DO I=1,IOBCN_GL
         I1 = NLID( I_OBC_GL(I) )
         IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP9(NCNT,:,:,:) = WQOBC_GL(I,:,:,:) ! KURT GLAESEMANN 3-MAR-2010
         ENDIF
       ENDDO
       IOBCN = NCNT

       IF(NCNT > 0)THEN
         ALLOCATE(I_OBC_N(NCNT))
         I_OBC_N  = TEMP1(1:NCNT)
         ALLOCATE(WQOBC(NCNT,KBM1,NUT_TM%NTIMES,NCP)) ; WQOBC=0.0 ! KURT GLAESEMANN 3-MAR-2010
         WQOBC(1:NCNT,:,:,:) = TEMP9(1:NCNT,:,:,:)  ! KURT GLAESEMANN 3-MAR-2010
       ENDIF

       DEALLOCATE(TEMP9) ! KURT GLAESEMANN 3-MAR-2010
       DEALLOCATE(TEMP1)
     ENDIF     
#    endif

      DEALLOCATE(WQOBC_GL)
      WRITE(*,*)'!  TEMP/SALT OBC READ      :    COMPLETE'
   ENDIF !!IOBCN_GL > 0

!==============================================================================|
!   REPORT AND CHECK RESULTS                                                   |
!==============================================================================|
   ALLOCATE(TEMP1(NPROCS))
   TEMP1(1)  = IOBCN

# if defined (MULTIPROCESSOR)
   CALL MPI_GATHER(IOBCN,  1,MPI_INTEGER,TEMP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
# endif

   
   IF(MSR)WRITE(*,100)'!  IOBCN                 :',IOBCN_GL,   (TEMP1(I),I=1,NPROCS)
   DEALLOCATE(TEMP1)

   RETURN
   100 FORMAT(1X,A26,I6," =>",2X,4(I5,1H,))
   END SUBROUTINE BCMAP
!==============================================================================|
   
   SUBROUTINE  BCMAP_DEALLOC
   
	!WLong moved here from wqm_main.F
		IF(ALLOCATED(I_OBC_N)) DEALLOCATE(I_OBC_N)
		IF(ALLOCATED(I_OBC_GL)) DEALLOCATE(I_OBC_GL)
		
	!WLong added deallocation WQOBC
		
		IF(ALLOCATED(WQOBC)) DEALLOCATE(WQOBC)
        
     !LB: these deallocations were missing 
        IF(ALLOCATED(NUT_TM%TIMES)) DEALLOCATE(NUT_TM%TIMES)

        
	 
   END SUBROUTINE BCMAP_DEALLOC
   
END MODULE MOD_BCMAP
