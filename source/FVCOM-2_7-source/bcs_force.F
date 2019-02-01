!==============================================================================|
!   Set Up the Following Boundary Conditions:                                  |
!     Bottom Freshwater (Groundwater) Info			               |
!     Tidal Forcing at Open Boundary					       |
!     Freshwater River Discharge					       |
!     Meteorological Forcing						       |
!==============================================================================|

   SUBROUTINE BCS_FORCE           

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE BCS
   USE MOD_CLOCK
   USE MOD_UTILS
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   USE MOD_OBCS
# if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif

#if defined(HEAT_FLUX)   
   USE MOD_HEATFLUX, ONLY : C_HFX
#endif    
   IMPLICIT NONE
   CHARACTER(LEN=80) :: COMT
   REAL(SP) :: QPREC,QEVAP,WDS,WDD,HFLUX,HSHORT,RBUF
   REAL(SP) :: TX,TY,CD,TTIME,BFWTIME
   REAL(SP) :: FTEMP1,FTEMP2,FTEMP3,RBUF1,RBUF2,RBUF3
   REAL(SP), ALLOCATABLE :: RTEMP(:),RTEMP1(:,:),RTEMP2(:,:),RTEMP3(:,:)
   REAL(SP), ALLOCATABLE :: RTEMP11(:),RTEMP22(:)
   INTEGER,  ALLOCATABLE :: TEMP(:),TEMP2(:),TEMP3(:),TEMP4(:),NODE_SBC(:)
   INTEGER,  ALLOCATABLE :: TEMPD(:,:),TEMP2D(:,:),TEMP3D(:,:)
   INTEGER   I,J,K,NQTIME,NBFWTIME,ISBCN1,INMAX,IOS,NCNT,IGL,IERR,JN
   CHARACTER(LEN=13) :: TSTRING

   !real(sp), allocatable :: WBTEMP(:),BLOCKED_WIND(:) ! Wind Block read in variable and 
   !INTEGER :: Nblocked  ! Number of blocked nodes, read from input 
   

!------------------------------------------------------------------------------|


!----------------------------REPORT--------------------------------------------!
   IF(MSR)WRITE(IPT,*  )'!'
   IF(MSR)WRITE(IPT,*)'!           SETTING UP PRESCRIBED BOUNDARY CONDITIONS   '
   IF(MSR)WRITE(IPT,*  )'!'

!==============================================================================|
!   Ground Water Information  BFWQDIS: m^3/s                                                 |
!==============================================================================|

!
!--------------determine global number of groundwater points and bcast---------!
!
   IF(MSR)THEN
     WRITE(IOPRT,*)'GROUNDWATER INFORMATION'
     READ(INBFW ,1000) COMT
     WRITE(IOPRT,1000) COMT
     READ(INBFW ,*) IBFW_GL
     WRITE(IOPRT,*) IBFW_GL
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_BCAST(IBFW_GL,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#  endif

   IBFW = 0
   IF(IBFW_GL > 0) THEN
     NCNT = 0

!
!--------------input node numbers for ground water inflow----------------------!
!
     ALLOCATE( NODE_BFW(IBFW_GL) )

     IF(MSR)THEN
     READ(INBFW ,*) (NODE_BFW(I),I=1,IBFW_GL)
     WRITE(IOPRT,*) (NODE_BFW(I),I=1,IBFW_GL)
!
!-----------------ensure all nodes exist in global domain----------------------!
!
     DO I=1,IBFW_GL
       IF(NODE_BFW(I) > NNodeGL)THEN
         WRITE(IPT,*)'==================ERROR=================================='
         WRITE(IPT,*)'GROUND WATER NODE NUMBER',I,'IS NOT IN THE GLOBAL DOMAIN'
         WRITE(IPT,*)'ENSURE GROUNDWATER NODES <= ',NNodeGL
         WRITE(IPT,*)'========================================================='
         CALL PSTOP
       END IF
     END DO
     END IF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(NODE_BFW,IBFW_GL,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif
!
!----------------------Shift To Local Domain If Parallel-----------------------!
!

     IF(SERIAL) IBFW = IBFW_GL

# if defined (MULTIPROCESSOR)
     IF(PAR)THEN
!       ALLOCATE(TEMP(IBFW_GL) , TEMP2(IBFW_GL) , TEMP3(IBFW_GL))
       ALLOCATE(TEMP2(IBFW_GL) , TEMP3(IBFW_GL))
       DO I=1,IBFW_GL
!         IF(NLID(TEMP(I)) /= 0)THEN
         IF(NLID(NODE_BFW(I)) /= 0)THEN
           NCNT = NCNT + 1
!           TEMP2(NCNT) = NLID(TEMP(I))
           TEMP2(NCNT) = NLID(NODE_BFW(I))
           TEMP3(NCNT) = I
         END IF
       END DO
       IBFW = NCNT

       DEALLOCATE(NODE_BFW)
       ALLOCATE(NODE_BFW(IBFW),BFW_GL2LOC(IBFW))
       NODE_BFW   = TEMP2(1:IBFW)
       BFW_GL2LOC = TEMP3(1:IBFW)
!       DEALLOCATE(TEMP,TEMP2,TEMP3)
       DEALLOCATE(TEMP2,TEMP3)
     END IF
#   endif


!
!----INPUT NUMBER OF DATA TIMES FOR GROUNDWATER DATA---------------------------!
!
     BFW_TM%LABEL = "Groundwater"
     IF(MSR)THEN
       READ(INBFW ,*) NBFWTIME
       WRITE(IOPRT,*) NBFWTIME
     END IF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(NBFWTIME,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

     BFW_TM%NTIMES = NBFWTIME
     ALLOCATE(BFW_TM%TIMES(NBFWTIME))

!
!----READ IN FRESH WATER FLUX AT EACH TIME=BFWTIME-----------------------------!
!
     ALLOCATE(RTEMP(IBFW_GL))
     ALLOCATE(RTEMP11(IBFW_GL))  
     ALLOCATE(RTEMP22(IBFW_GL)) 

     ALLOCATE(BFWQDIS(IBFW,NBFWTIME))  ; BFWQDIS = 0.0_SP
     ALLOCATE(BFWQTDIS(IBFW,NBFWTIME)) ; BFWQTDIS = 0.0_SP 
     ALLOCATE(BFWQSDIS(IBFW,NBFWTIME)) ; BFWQSDIS = 0.0_SP 

     DO I=1,NBFWTIME
       IF(MSR)THEN
         READ(INBFW,*) BFWTIME
         WRITE(IOPRT,5000) BFWTIME
         BFW_TM%TIMES(I) = BFWTIME
         READ(INBFW,*) (RTEMP(J),J = 1,IBFW_GL)
         READ(INBFW,*) (RTEMP11(J),J = 1,IBFW_GL)
         READ(INBFW,*) (RTEMP22(J),J = 1,IBFW_GL)
       END IF

#      if defined (MULTIPROCESSOR)
       IF(PAR)CALL MPI_BCAST(RTEMP,IBFW_GL,MPI_F,0,MPI_COMM_WORLD,IERR)
       IF(PAR)CALL MPI_BCAST(RTEMP11,IBFW_GL,MPI_F,0,MPI_COMM_WORLD,IERR)
       IF(PAR)CALL MPI_BCAST(RTEMP22,IBFW_GL,MPI_F,0,MPI_COMM_WORLD,IERR)
#      endif

       IF(SERIAL)BFWQDIS(1:IBFW_GL,I)  = RTEMP(1:IBFW_GL)
       IF(SERIAL)BFWQTDIS(1:IBFW_GL,I) = RTEMP11(1:IBFW_GL)
       IF(SERIAL)BFWQSDIS(1:IBFW_GL,I) = RTEMP22(1:IBFW_GL)

#      if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         DO J=1,IBFW
           BFWQDIS(J,I)  = RTEMP(BFW_GL2LOC(J))
           BFWQTDIS(J,I) = RTEMP11(BFW_GL2LOC(J))
           BFWQSDIS(J,I) = RTEMP22(BFW_GL2LOC(J))
         END DO
       END IF
#      endif

       IF(MSR)WRITE(IOPRT,5000) (RTEMP(J),J = 1,IBFW_GL) 
       IF(MSR)WRITE(IOPRT,5000) (RTEMP11(J),J = 1,IBFW_GL) 
       IF(MSR)WRITE(IOPRT,5000) (RTEMP22(J),J = 1,IBFW_GL) 
     END DO
     DEALLOCATE(RTEMP,RTEMP11,RTEMP22)


#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(BFW_TM%TIMES,NBFWTIME,MPI_F,0,MPI_COMM_WORLD,IERR)
#    endif

   END IF !!IBFW_GL > 0

!
!--REPORT RESULTS--------------------------------------------------------------!
!
   ALLOCATE(TEMP(NPROCS))
   TEMP(1)  = IBFW

# if defined (MULTIPROCESSOR)
   CALL MPI_GATHER(IBFW,1,MPI_INTEGER,TEMP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
# endif

   IF(IBFW_GL == 0)THEN
     IF(MSR)WRITE(IPT,*)'!  GROUNDWATER FLUX      :    NONE'
   ELSE
   IF(MSR)WRITE(IPT,*)'!'
   IF(MSR)WRITE(IPT,100)'!  GROUNDWATER POINTS    :',IBFW_GL, (TEMP(I),I=1,NPROCS)
   IF(NBFWTIME > 0)THEN
     IF(MSR)WRITE(IPT,101)'!  GWATER DATA BEGIN     :',BFW_TM%TIMES(1)
     IF(MSR)WRITE(IPT,101)'!  GWATER DATA END       :',BFW_TM%TIMES(NBFWTIME)
   END IF
   END IF
   DEALLOCATE(TEMP)

!==============================================================================|
!   Input Non-Julian Open Boundary Tidal Forcing                               |
!==============================================================================|

# if !defined (MEAN_FLOW)
   IF(IBCN_GL(1) > 0)THEN
# endif

   IF(S_TYPE == 'non-julian') THEN
     READ(INOEL ,1000) COMT
     IF(MSR)WRITE(IOPRT,*) 'Non-Julian Tide Information'
     IF(MSR)WRITE(IOPRT,1000) COMT

     READ(INOEL,*) ISBCN1
     IF(MSR)WRITE(IOPRT,*) ISBCN1

!
!-------ENSURE SAME NUMBER OF SPECIFIED OPEN BOUNDARY POINTS AS FILE-casename_obc.dat----|
!
     IF(ISBCN1 /= IBCN_GL(1))THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'NUMBER OF OPEN BOUNDARY POINTS IN OPEN BOUNDARY SURFACE'
       WRITE(IPT,*)'ELEVATION FILE IS LARGER THAN NUMBER OF OPEN BOUNDARY '
       WRITE(IPT,*)'POINTS OF PRESCRIBED ELEVATION TYPE IN CASENAME_obc.dat'
       WRITE(IPT,*) 'SEE SUBROUTINE BCS_FORCE'
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

!
!----READ IN BOUNDARY POINTS, AMPLITUDES, AND PHASES OF TIDE-------------------|
!
     ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
     ALLOCATE(APT(IBCN_GL(1),8), PHAI(IBCN_GL(1),8))
     APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP
     NCNT = 0
#   if defined (EQUI_TIDE)       
     APT_FACT_EQUI = 0
#   endif       
#   if defined (ATMO_TIDE)       
     APT_FACT_ATMO = 0
#   endif       
     DO I=1,IBCN_GL(1)
       READ(INOEL,*)  NODE_SBC(I),EMEAN(I)
       READ (INOEL,*) (APT(I,J), J=1,8)
       READ (INOEL,*) (PHAI(I,J), J=1,8)
#   if defined (EQUI_TIDE)       
       IF(APT(I,1) > 0.0_SP)APT_FACT_EQUI(1) = 1
       IF(APT(I,2) > 0.0_SP)APT_FACT_EQUI(2) = 1
       IF(APT(I,3) > 0.0_SP)APT_FACT_EQUI(3) = 1
       IF(APT(I,4) > 0.0_SP)APT_FACT_EQUI(4) = 1
       IF(APT(I,5) > 0.0_SP)APT_FACT_EQUI(5) = 1
       IF(APT(I,6) > 0.0_SP)APT_FACT_EQUI(6) = 1
       IF(APT(I,7) > 0.0_SP)APT_FACT_EQUI(7) = 1
       IF(APT(I,8) > 0.0_SP)APT_FACT_EQUI(8) = 1
#   endif       
#   if defined (ATMO_TIDE)       
       IF(APT(I,1) > 0.0_SP)APT_FACT_ATMO = 1
#   endif       
       IF(MSR)WRITE(IOPRT,*) NODE_SBC(I),EMEAN(I)
       IF(MSR)WRITE(IOPRT,*) (APT(I,J), J=1,8)
       IF(MSR)WRITE(IOPRT,*) (PHAI(I,J), J=1,8)
     END DO

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!

#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
     ALLOCATE( TEMP2(IBCN_GL(1)) ,RTEMP(IBCN_GL(1)))
     ALLOCATE( RTEMP1(IBCN_GL(1),8) , RTEMP2(IBCN_GL(1),8))
     NCNT = 0
     DO I=1,IBCN_GL(1)
       IF(NLID(NODE_SBC(I)) /= 0)THEN
         NCNT = NCNT + 1
         TEMP2(NCNT)     = NLID(NODE_SBC(I))
         RTEMP(NCNT)     = EMEAN(I)
         RTEMP1(NCNT,1:8) = APT(I,1:8)
         RTEMP2(NCNT,1:8) = PHAI(I,1:8)
       END IF
     END DO

     IF(NCNT /= IBCN(1))THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'LOCAL OPEN BOUNDARY NODE COUNTS DIFFER BETWEEN TIDE'
       WRITE(IPT,*)'FORCING AND OPEN BOUNDARY NODE FILES'
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     DEALLOCATE(NODE_SBC,EMEAN,APT,PHAI)
     IF(IBCN(1) > 0)THEN
       ALLOCATE(NODE_SBC(IBCN(1)),EMEAN(IBCN(1)))
       ALLOCATE(APT(IBCN(1),8),PHAI(IBCN(1),8))
       NODE_SBC = TEMP2(1:NCNT)
       EMEAN    = RTEMP(1:NCNT)
       APT      = RTEMP1(1:NCNT,1:8)
       PHAI     = RTEMP2(1:NCNT,1:8)
     ELSE
       ALLOCATE(NODE_SBC(1),EMEAN(1))
       ALLOCATE(APT(1,8),PHAI(1,8))
       NODE_SBC = 0.0_SP ; EMEAN = 0.0_SP ; APT = 0.0_SP ; PHAI = 0.0_SP
     END IF

     DEALLOCATE(TEMP2,RTEMP,RTEMP1,RTEMP2)

     END IF !!PAR
#    endif

!
!----MAKE SURE LOCAL NODE NUMBERS OF SPECIFIED NODES MATCHES LOCAL NODE--------|
!----NUMBER OF SPECIFIED NODES IN obc.dat FILE---------------------------------|
!
     DO I=1,IBCN(1)
       JN = OBC_LST(1,I)
       IF(NODE_SBC(I) /= I_OBC_N(JN))THEN
         WRITE(IPT,*)'==================ERROR=================================='
         WRITE(IPT,*)'LOCAL OPEN BOUNDARY NODE LIST DIFFERS BETWEEN TIDE'
         WRITE(IPT,*)'FORCING AND OPEN BOUNDARY NODE (TYPE 1 OR 2) FILES'
         WRITE(IPT,*)'========================================================='
         WRITE(IPT,*)NODE_SBC(I),I_OBC_N(JN)
         CALL PSTOP
       END IF
     END DO

!
!----MODIFY AMPLITUDE AND PHASE------------------------------------------------|
!
     APT = APT/100.0_SP
     PHAI = MOD(PHAI,360.0_SP)

     CLOSE(INOEL)
!
!--REPORT RESULTS--------------------------------------------------------------!
!
   RBUF = MAXVAL(APT)
#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(MAXVAL(APT),RBUF,1,MPI_F,MPI_MAX,0,MPI_COMM_WORLD,IERR)
#  endif
   IF(MSR)WRITE(IPT,*)'!'
   IF(MSR)WRITE(IPT,*  )'!  NON-JULIAN TIDE       :    SET'
   IF(MSR)WRITE(IPT,101)'!  MAX TIDE AMPLITUDE    : ',RBUF

!==============================================================================|
!   Input Julian Open Boundary Tidal Forcing                                   |
!==============================================================================|


   ELSE IF(S_TYPE == 'julian')THEN

!
!----Count Number of Data Series Available-------------------------------------!
!
     IF(MSR)THEN
       NCNT = 0
       DO WHILE(.TRUE.)
         READ(INJUL,*,IOSTAT=IOS)
         IF(IOS < 0)EXIT
         NCNT = NCNT + 1
       END DO

       IF(NCNT == 0)CALL PERROR(6,"JULIAN TIDE SELECTED BUT NO DATA IN FILE")
     END IF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(NCNT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

     ELO_TM%NTIMES = NCNT

     IF(IBCN(1) > 0)THEN
       ALLOCATE(ELSBC(IBCN(1),NCNT))
       ALLOCATE(ELO_TM%TIMES(NCNT))
     ELSE
       ALLOCATE(ELSBC(1,NCNT))
       ALLOCATE(ELO_TM%TIMES(NCNT))
     END IF
     ELSBC = 0.0_SP

     TTIME = 0.0_SP
     DO I=1,NCNT
       ELO_TM%TIMES(I) = TTIME
       TTIME = TTIME + DELTT
     END DO

!
!----Read in Data Series ------------------------------------------------------!
!

     REWIND(INJUL)
!     IF(IBCN_GL(1) > 200) THEN
!       WRITE(IPT,*)'CHANGE FORMAT STATEMENT BELOW TO ACCOMODATE'
!       WRITE(IPT,*)'IBCN_GL(1) NUMBER OF NODES AND RECOMPILE'
!       CALL PSTOP
!     END IF

     ALLOCATE(RTEMP1(IBCN_GL(1),ELO_TM%NTIMES))

     IF(MSR)THEN
       DO I=1,ELO_TM%NTIMES
         READ(INJUL,*) (RTEMP1(J,I),J=1,IBCN_GL(1))
       END DO
     END IF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(RTEMP1,IBCN_GL(1)*ELO_TM%NTIMES,MPI_F,0,MPI_COMM_WORLD,IERR)
#    endif

!
!---Map to Local Domain--------------------------------------------------------!
!

     IF(SERIAL) ELSBC = RTEMP1

     IF(IBCN(1) > 0)THEN
#    if defined (MULTIPROCESSOR)
     NCNT = 0
     IF(PAR)THEN
       DO  J=1,IBCN_GL(1)
!        IF(NLID(I_OBC_GL(J)) /= 0 .AND. (TYPE_OBC_GL(J) == 1 .OR. TYPE_OBC_GL(J) ==2))THEN
         JN=I_OBC_GL(OBC_LST_GL(1,J))
         IF(NLID(JN) /= 0 .AND. (TYPE_OBC_GL(OBC_LST_GL(1,J)) == 1 .OR.  &
	    TYPE_OBC_GL(OBC_LST_GL(1,J)) == 2))THEN
           NCNT = NCNT +1
           ELSBC(NCNT,:) = RTEMP1(J,:)
         END IF
       END DO
     END IF
#   endif
     END IF

     DEALLOCATE(RTEMP1)

!
!--REPORT RESULTS--------------------------------------------------------------!
!

   IF(MSR)WRITE(IPT,*)'!'
   IF(MSR)WRITE(IPT,*  )'!  JULIAN TIDE           :    SET'
   RBUF = 0.
   IF(IBCN(1) > 0)RBUF = MAXVAL(ELSBC)
#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(MAXVAL(ELSBC),RBUF,1,MPI_F,MPI_MAX,0,MPI_COMM_WORLD,IERR)
#  endif
   IF(MSR)WRITE(IPT,101)'!  MAX TIDE AMPLITUDE    : ',RBUF
   IF(ELO_TM%NTIMES > 0)THEN
     CALL GETTIME(TSTRING,INT(ELO_TM%TIMES(1)))
     IF(MSR)WRITE(IPT,102)'!  TIDAL DATA BEGIN      :  ',TSTRING
     CALL GETTIME(TSTRING,INT(ELO_TM%TIMES(ELO_TM%NTIMES)))
     IF(MSR)WRITE(IPT,102)'!  TIDAL DATA END        :  ',TSTRING
   END IF

   END IF    !!JULIAN
# if !defined (MEAN_FLOW)
   END IF    !!IBCN_GL(1) > 0
# endif

#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#  endif
!==============================================================================|
!   Input River/Dam/Intake/Outfall Boundary Values                             |
!==============================================================================|

!
!-------Check Selected Combination for Validity--------------------------------!
!
   REWIND(INRIV)
   READ(INRIV,'(A4,2X,A10)') INFLOW_TYPE,POINT_ST_TYPE
   IF(MSR)WRITE(IOPRT,*) 'River Inflow Information'
   IF(MSR)WRITE(IOPRT,*) 'INFLOW_TYPE==',INFLOW_TYPE
   IF(MSR)WRITE(IOPRT,*) 'POINT_ST_TYPE==',POINT_ST_TYPE

   IF(INFLOW_TYPE /= 'edge' .AND. INFLOW_TYPE /= 'node') THEN
     CALL PERROR(6,"INFLOW TYPE NOT CORRECT","SHOULD BE edge or node")
   END IF

   IF(POINT_ST_TYPE /= 'calculated' .AND. POINT_ST_TYPE /= 'specified') THEN
     CALL PERROR(6,"POINT_ST TYPE NOT CORRECT","SHOULD BE calculated or specified")
   END IF

!
!--Read in Number of Discharge Nodes/Edges-------------------------------------!
!
   IF(MSR)THEN
     READ(INRIV,*) NUMQBC_GL
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_BCAST(NUMQBC_GL,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#  endif

   NUMQBC = 0
   IF(NUMQBC_GL > 0)THEN
!
!--Shut off Temp/Salinity Averaging if River Flux is of type "specified"
!
!   IF(POINT_ST_TYPE == 'specified' .AND. TS_FCT)THEN
!     IF(MSR)THEN
!       WRITE(IPT,*)'=========WARNING================'
!       WRITE(IPT,*)'RIVER QUANTITIES ARE "specified"'
!       WRITE(IPT,*)'DEACTIVATING TS_FCT'
!       WRITE(IPT,*)'================================'
!     END IF
!     TS_FCT = .FALSE.
!   END IF
!
!--Read in Freshwater Discharge Nodes------------------------------------------!
!
     ALLOCATE(TEMP(NUMQBC_GL),TEMP2(NUMQBC_GL),TEMP3(NUMQBC_GL))
     IF(MSR)THEN
       DO I=1,NUMQBC_GL
         READ(INRIV,*) TEMP(I)
       END DO
     END IF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(TEMP,NUMQBC_GL,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

!
!--Determine Global--> Local Mapping of Freshwater Discharge Nodes
!
     IF(SERIAL)THEN
       NUMQBC = NUMQBC_GL
       IF(INFLOW_TYPE == 'node') THEN
         ALLOCATE(INODEQ(NUMQBC))
         INODEQ = TEMP
       ELSE IF(INFLOW_TYPE == 'edge') THEN
         ALLOCATE(ICELLQ(NUMQBC))
!         ICELLQ = TEMP(1:NCNT)
         ICELLQ = TEMP(1:NUMQBC)
       END IF
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       NCNT = 0
       IF(INFLOW_TYPE == 'node') THEN
         DO I=1,NUMQBC_GL
           IF(NLID(TEMP(I)) /= 0)THEN
             NCNT = NCNT + 1
             TEMP2(NCNT) = NLID(TEMP(I))
             TEMP3(NCNT) = I
           END IF
         END DO
         NUMQBC = NCNT
         ALLOCATE(INODEQ(NUMQBC),RIV_GL2LOC(NUMQBC))
         INODEQ = TEMP2(1:NCNT)
         RIV_GL2LOC = TEMP3(1:NCNT)
       ELSE IF(INFLOW_TYPE == 'edge') THEN
         DO I=1,NUMQBC_GL
           IF(ELID(TEMP(I)) /= 0)THEN
             NCNT = NCNT + 1
             TEMP2(NCNT) = ELID(TEMP(I))
             TEMP3(NCNT) = I
           END IF
         END DO
         NUMQBC = NCNT
         ALLOCATE(ICELLQ(NUMQBC),RIV_GL2LOC(NUMQBC))
         ICELLQ = TEMP2(1:NCNT)
         RIV_GL2LOC = TEMP3(1:NCNT)
       END IF
     END IF
#   endif

     DEALLOCATE(TEMP,TEMP2,TEMP3)


!
!----Read in Freshwater Flux Vertical Distribution-----------------------------!
!

     ALLOCATE(RTEMP1(NUMQBC_GL,KBM1))
     IF(MSR)THEN
       DO I = 1, NUMQBC_GL
         READ(INRIV ,*) J,(RTEMP1(I,K),K = 1,KBM1)
         WRITE(IOPRT,*) J,(RTEMP1(I,K),K = 1,KBM1)
       END DO
     END IF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(RTEMP1,NUMQBC_GL*KBM1,MPI_F,0,MPI_COMM_WORLD,IERR)
#    endif

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     IF(NUMQBC > 0)THEN
     ALLOCATE(VQDIST(NUMQBC,KBM1))

     IF(SERIAL) VQDIST = RTEMP1

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,NUMQBC
         DO K=1,KBM1
           VQDIST(I,K) = RTEMP1(RIV_GL2LOC(I),K)
         END DO
       END DO
     END IF
#   endif
     END IF

     DEALLOCATE(RTEMP1)

!
!----Read in Time Dependent DataSets (DQDIS,DSDIS,DTDIS)------------------------!
!
     IF(MSR)READ(INRIV,*) NQTIME

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(NQTIME,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

     QBC_TM%NTIMES = NQTIME
     QBC_TM%LABEL  = "Freshwater Discharge"
     ALLOCATE(QBC_TM%TIMES(NQTIME))
     ALLOCATE(RTEMP1(NUMQBC_GL,NQTIME))
     ALLOCATE(RTEMP2(NUMQBC_GL,NQTIME))
     ALLOCATE(RTEMP3(NUMQBC_GL,NQTIME))

     IF(MSR)THEN
       DO I = 1, NQTIME
         READ(INRIV,*) TTIME
         QBC_TM%TIMES(I) = TTIME
         READ(INRIV,*) (RTEMP1(J,I),J = 1,NUMQBC_GL)
         READ(INRIV,*) (RTEMP2(J,I),J = 1,NUMQBC_GL)
         READ(INRIV,*) (RTEMP3(J,I),J = 1,NUMQBC_GL)
         WRITE(IOPRT,5000) TTIME
         WRITE(IOPRT,5000) (RTEMP1(J,I),J = 1,NUMQBC_GL)
         WRITE(IOPRT,5000) (RTEMP2(J,I),J = 1,NUMQBC_GL)
         WRITE(IOPRT,5000) (RTEMP3(J,I),J = 1,NUMQBC_GL)
       END DO
     END IF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(QBC_TM%TIMES,NQTIME,MPI_F,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_BCAST(RTEMP1,NUMQBC_GL*NQTIME,MPI_F,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_BCAST(RTEMP2,NUMQBC_GL*NQTIME,MPI_F,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_BCAST(RTEMP3,NUMQBC_GL*NQTIME,MPI_F,0,MPI_COMM_WORLD,IERR)
#    endif

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     IF(NUMQBC > 0)THEN
       ALLOCATE(DQDIS(NUMQBC,NQTIME))
       ALLOCATE(DTDIS(NUMQBC,NQTIME))
       ALLOCATE(DSDIS(NUMQBC,NQTIME))

       IF(SERIAL)THEN
         DQDIS(1:NUMQBC_GL,:) = RTEMP1(1:NUMQBC_GL,:)
         DTDIS(1:NUMQBC_GL,:) = RTEMP2(1:NUMQBC_GL,:)
         DSDIS(1:NUMQBC_GL,:) = RTEMP3(1:NUMQBC_GL,:)
       END IF

#     if defined (MULTIPROCESSOR)
       IF(PAR)THEN
       DO I=1,NQTIME
         DQDIS(1:NUMQBC,I) = RTEMP1(RIV_GL2LOC(1:NUMQBC),I)
         DTDIS(1:NUMQBC,I) = RTEMP2(RIV_GL2LOC(1:NUMQBC),I)
         DSDIS(1:NUMQBC,I) = RTEMP3(RIV_GL2LOC(1:NUMQBC),I)
  

       END DO
       END IF
#     endif



!write(*,*)'Discharge that is scaled is '
 !       write(*,*) DQDIS(:,:)!, DQDIS(1,1:10) 
 
!DQDIS=DQDIS*2.
! DTDIS= 0.0
  !DSDIS=0.0 
! B Clark

     END IF

     DEALLOCATE(RTEMP1,RTEMP2,RTEMP3)

   CLOSE(INRIV)
!
!--REPORT RESULTS--------------------------------------------------------------!
!
   ALLOCATE(TEMP(NPROCS))
   TEMP(1)  = NUMQBC
   FTEMP1 = 0.0_SP; FTEMP2 = 0.0_SP; FTEMP3 = 0.0_SP;
   IF(NUMQBC > 0) FTEMP1 = MAXVAL(DQDIS)
   IF(NUMQBC > 0) FTEMP2 = MAXVAL(DTDIS)
   IF(NUMQBC > 0) FTEMP3 = MAXVAL(DSDIS)
   RBUF1 = FTEMP1 ; RBUF2 = FTEMP2 ; RBUF3 = FTEMP3

# if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_GATHER(NUMQBC,1,MPI_INTEGER,TEMP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_REDUCE(FTEMP1,RBUF1,1,MPI_F,MPI_MAX,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_REDUCE(FTEMP2,RBUF2,1,MPI_F,MPI_MAX,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_REDUCE(FTEMP3,RBUF3,1,MPI_F,MPI_MAX,0,MPI_COMM_WORLD,IERR)
# endif

   END IF !! NUMQBC_GL > 0

   IF(MSR)WRITE(IPT,*)'!'
   IF(NUMQBC_GL == 0)THEN
     IF(MSR)WRITE(IPT,*)'!  FRESHWATER FLUX       :    NONE'
   ELSE
     IF(MSR)WRITE(IPT,100)'!  FRESHWATER POINTS     :',NUMQBC_GL, (TEMP(I),I=1,NPROCS)
     IF(MSR)CALL GETTIME(TSTRING,3600*INT(QBC_TM%TIMES(1)))
     IF(MSR)WRITE(IPT,102)'!  FWATER DATA BEGIN     :  ',TSTRING
     IF(MSR)CALL GETTIME(TSTRING,3600*INT(QBC_TM%TIMES(QBC_TM%NTIMES)))
     IF(MSR)WRITE(IPT,102)'!  FWATER DATA END       :  ',TSTRING
     IF(MSR)WRITE(IPT,101)'!  MAX DQDIS             :',RBUF1
     IF(MSR)WRITE(IPT,101)'!  MAX DTDIS             :',RBUF2
     IF(MSR)WRITE(IPT,101)'!  MAX DSDIS             :',RBUF3
     DEALLOCATE(TEMP)
   END IF


!==============================================================================|
!   Input Meteorological Boundary Conditions                                   |
!==============================================================================|
!    precipitation: mm/s       "qprec"                                         |
!    evaporation:   mm/s       "qevap"                                         |
!    wind:          wds (speed) wdd (direction)       			       |
!    heat flux:     w/m^2                              			       |
!==============================================================================|

   IF(M_TYPE == 'uniform')THEN

!==============================================================================|
!   UNIFORM METEOLOGICAL CONDITIONS                                            |
!==============================================================================|

     READ(INCWH,1000) COMT
     IF(MSR)WRITE(IOPRT,*)'Meteorological Forcing Info'
     IF(MSR)WRITE(IOPRT,1000) COMT

!
!----Determine Number of Data Times--------------------------------------------!
!
     NCNT = 0
     DO WHILE(.TRUE.)
       READ(INCWH,*,END=15,IOSTAT=IOS)
       READ(INCWH,*,END=15,IOSTAT=IOS)
       IF(IOS < 0)EXIT
       NCNT = NCNT + 1
     END DO
 15  CONTINUE
     IF(NCNT == 0)CALL PERROR(6,"NO UNIFORM METEO DATA PROVIDED")

     REWIND(INCWH) ; READ(INCWH,*)

!
!----Read in Precipitation/Evap/Wind/Heat Flux/Radiation Data at Each Time-----!
!

     UMF_TM%NTIMES = NCNT
     ALLOCATE(UMF_TM%TIMES(NCNT))
     ALLOCATE(UQPREC(NCNT),UQEVAP(NCNT))
     ALLOCATE(UWIND(NCNT),VWIND(NCNT))
     ALLOCATE(UHFLUX(NCNT),UHSHORT(NCNT))

     DO I=1,NCNT
       READ(INCWH ,*) TTIME
       IF(MSR)WRITE(IOPRT,*) TTIME
       UMF_TM%TIMES(I) = TTIME

       READ(INCWH ,*) QPREC, QEVAP, WDS, WDD, HFLUX,HSHORT

       IF(MSR)WRITE(IOPRT,5000) QPREC, QEVAP, WDS, WDD, HFLUX,HSHORT

!       WRITE(*,*) QPREC, QEVAP, WDS, WDD, HFLUX,HSHORT
!       UQPREC(I) = QPREC / (86400.0_SP*365.0_SP)
!       UQEVAP(I) = QEVAP / (86400.0_SP*365.0_SP)
!       UQPREC(I) = QPREC / 1000.0_SP
!       UQEVAP(I) = QEVAP / 1000.0_SP
       UQPREC(I) =  QPREC
       UQEVAP(I) = QEVAP

       WDD = MOD(WDD,360.0_SP)
!       UWIND(I) = WDS * COS(6.28319_SP*WDD/360.0_SP)
!       VWIND(I) = WDS * SIN(6.28319_SP*WDD/360.0_SP)
       UWIND(I) = -WDS * SIN(WDD*DEG2RAD)
       VWIND(I) = -WDS * COS(WDD*DEG2RAD)
! B Clark
       UHFLUX(I)  = HFLUX
       UHSHORT(I) = HSHORT
     END DO

     IF(WINDTYPE /= 'speed' .AND. WINDTYPE /='stress') THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'NO UNIFORM METEO DATA PROVIDED'
       WRITE(IPT,*)'WINDTYPE IS NOT CORRECT, --->',WINDTYPE
       WRITE(IPT,*)'MUST BE "speed" or "stress"'
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

     CLOSE(INCWH)

!UWIND=0.0
!VWIND=0.0
!QPREC=0.0
!QEVAP=0.0


!
!--REPORT RESULTS--------------------------------------------------------------!
!

   IF(MSR)THEN
     WRITE(IPT,*)'!'
     WRITE(IPT,*    )'!  UNIFORM METEO         :    SET'
      IF(UMF_TM%NTIMES > 0)THEN
        CALL GETTIME(TSTRING,INT(3600.*UMF_TM%TIMES(1)))
        WRITE(IPT,102)'!  METEO DATA BEGIN      :  ',TSTRING
        CALL GETTIME(TSTRING,INT(3600.*UMF_TM%TIMES(UMF_TM%NTIMES)))
        WRITE(IPT,102)'!  METEO DATA END        :  ',TSTRING
      END IF
    END IF
	
	
! Wind Block on a uniformly forced model domain, for non uniform, use can specify nodes
	!that have 0 velocity during pre-process
! B Clark Nov 2015	
	IF(WIND_BLOCK) THEN
	   READ(INWNDBLK,*)Nblocked
     ALLOCATE(WBTEMP(Nblocked))
       DO J=1,Nblocked
         READ(INWNDBLK,*) WBTEMP(J)   ! THe elements to block wind for
!		 write(IPT,*)'WBTEMP = ',WBTEMP(J)
       END DO

!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|

     IF(SERIAL)THEN
        DO J = 1, Nblocked
		  UWIND(WBTEMP(j)) = 0.0   ! Now all elements with wind blocked have 0 windspeed
		  VWIND(WBTEMP(j)) = 0.0
       END DO
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,MElem
	       DO j = 1,Nblocked
               IF(EGID(I) == WBTEMP(J))THEN
			      UWIND(I) = 0.0
				  VWIND(I) = 0.0
			    END IF
			END DO
       END DO
     END IF
#   endif
   END IF
   
!   DEALLOCATE(WBTEMP)


!==============================================================================|
!   NON-UNIFORM METEOLOGICAL CONDITIONS                                        |
!==============================================================================|

   ELSE IF (M_TYPE == 'non-uniform')THEN 


!=====================HEAT FLUX/SHORT WAVE RADIATION===========================!
#if defined(HEAT_FLUX)
IF(.NOT. C_HFX)THEN
#endif
     REWIND(INHFX)
!
!----Input Number of Data Times for Heat Flux and Short Wave Radiation---------!
!
     IF(MSR)THEN
     NCNT = 0
     DO WHILE(.TRUE.)
       READ(INHFX,*,END=10)FTEMP1
       READ(INHFX,*)
       NCNT = NCNT + 1
     END DO
 10  CONTINUE
     REWIND(INHFX)

     IF(NCNT == 0)CALL PERROR(6,"NO DATA PROVIDED FOR HEAT FLUX AND SHORT WAVE RAD")
     END IF

     HFX_TM%NTIMES = NCNT 

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(HFX_TM%NTIMES,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

!
!----Read in Data Times and Global Heat Flux/Short Wave Radiation Data---------!
!

     ALLOCATE(HFX_TM%TIMES(HFX_TM%NTIMES))
     ALLOCATE(RTEMP1(NNodeGL,HFX_TM%NTIMES),RTEMP2(NNodeGL,HFX_TM%NTIMES))

     IF(MSR)THEN
       DO J=1,HFX_TM%NTIMES
         READ(INHFX,*) HFX_TM%TIMES(J) 
         READ(INHFX,*) (RTEMP1(I,J),RTEMP2(I,J),I=1,NNodeGL)
       END DO
     END IF
!
!----Broadcast Data------------------------------------------------------------!
!

#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       CALL MPI_BCAST(HFX_TM%TIMES,HFX_TM%NTIMES,MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP1,HFX_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP2,HFX_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
     END IF
#    endif

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     ALLOCATE(DHFLUX(NNode,HFX_TM%NTIMES),DHSHORT(NNode,HFX_TM%NTIMES))

     IF(SERIAL)THEN
       DHFLUX(1:NNodeGL,:)  = RTEMP1(1:NNodeGL,:)
       DHSHORT(1:NNodeGL,:) = RTEMP2(1:NNodeGL,:)
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,NNode
         DHFLUX(I,:)  = RTEMP1(NGID(I),:) 
         DHSHORT(I,:) = RTEMP2(NGID(I),:) 
       END DO
     END IF
#   endif

     DEALLOCATE(RTEMP1,RTEMP2)
     IF(MSR)WRITE(IPT,101)'!  HFLUX/SWRAD READ      :    COMPLETE'
#if defined(HEAT_FLUX)      
ENDIF
#endif
!=====================TIME DEPENDENT WIND FIELD================================!


     REWIND(INWND)
!
!----Input Number of Data Times for Wind Field---------------------------------!
!
     IF(MSR)THEN
     NCNT = 0
     DO WHILE(.TRUE.)
       READ(INWND,*,END=20)FTEMP1
       READ(INWND,*) 
       NCNT = NCNT + 1       
     END DO
 20  CONTINUE
     REWIND(INWND)

     IF(NCNT == 0)CALL PERROR(6,"NO DATA PROVIDED FOR SURFACE WIND FIELD")
     END IF

     WND_TM%NTIMES = NCNT 

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(WND_TM%NTIMES,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

!
!----Read in Data Times and Global Wind Data-----------------------------------!
!

     ALLOCATE(WND_TM%TIMES(WND_TM%NTIMES))
     ALLOCATE(RTEMP1(MElemGL,WND_TM%NTIMES),RTEMP2(MElemGL,WND_TM%NTIMES))

     IF(MSR)THEN
       DO J=1,WND_TM%NTIMES
         READ(INWND,*) WND_TM%TIMES(J) 
         READ(INWND,*) (RTEMP1(I,J),RTEMP2(I,J),I=1,MElemGL)
       END DO
     END IF
!
!----Send/Recv Data------------------------------------------------------------!
!

#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       CALL MPI_BCAST(WND_TM%TIMES,WND_TM%NTIMES,MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP1,WND_TM%NTIMES*MElemGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP2,WND_TM%NTIMES*MElemGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
     END IF
#    endif

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     ALLOCATE(DTX(MElem,WND_TM%NTIMES),DTY(MElem,WND_TM%NTIMES))

     IF(SERIAL)THEN
       DTX(1:MElemGL,:)  = RTEMP1(1:MElemGL,:)
       DTY(1:MElemGL,:) = RTEMP2(1:MElemGL,:)
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,MElem
         DTX(I,:) = RTEMP1(EGID(I),:) 
         DTY(I,:) = RTEMP2(EGID(I),:) 
       END DO
     END IF
#   endif

     DEALLOCATE(RTEMP1,RTEMP2)
     IF(MSR)WRITE(IPT,101)'!  WIND FIELD READ       :    COMPLETE'
      


!=====================TIME DEPENDENT EVAPORATION AND PRECIPITATION=============!

     IF(EVP_FLAG)THEN
       REWIND(INEVP)
!
!----Input Number of Data Times for Evaporation and Precipitation--------------!
!
       IF(MSR)THEN
         NCNT = 0
         DO WHILE(.TRUE.)
           READ(INEVP,*,END=30)FTEMP1
           READ(INEVP,*) 
           NCNT = NCNT + 1       
         END DO
 30      CONTINUE
         REWIND(INEVP)

         IF(NCNT == 0)CALL PERROR(6,"NO DATA PROVIDED FOR EVAPORATION AND PRECIPITATION")
       END IF

       EVP_TM%NTIMES = NCNT 

#    if defined (MULTIPROCESSOR)
       IF(PAR)CALL MPI_BCAST(EVP_TM%NTIMES,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

!
!----Read in Data Times and Global Evaporation and Precipitation---------------!
!

       ALLOCATE(EVP_TM%TIMES(EVP_TM%NTIMES))
       ALLOCATE(RTEMP1(NNodeGL,EVP_TM%NTIMES),RTEMP2(NNodeGL,EVP_TM%NTIMES))

       IF(MSR)THEN
         DO J=1,EVP_TM%NTIMES
           READ(INEVP,*) EVP_TM%TIMES(J) 
           READ(INEVP,*) (RTEMP1(I,J),RTEMP2(I,J),I=1,NNodeGL)
         END DO
       END IF
!
!----Send/Recv Data------------------------------------------------------------!
!

#    if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         CALL MPI_BCAST(EVP_TM%TIMES,EVP_TM%NTIMES,MPI_F,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RTEMP1,EVP_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RTEMP2,EVP_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       END IF
#    endif

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
       ALLOCATE(DQEVAP(NNode,EVP_TM%NTIMES),DQPREC(NNode,EVP_TM%NTIMES))

       IF(SERIAL)THEN
         DQEVAP(1:NNodeGL,:) = RTEMP1(1:NNodeGL,:)
         DQPREC(1:NNodeGL,:) = RTEMP2(1:NNodeGL,:)
       END IF

#   if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         DO I=1,NNode
           DQEVAP(I,:) = RTEMP1(NGID(I),:) 
           DQPREC(I,:) = RTEMP2(NGID(I),:) 
         END DO
       END IF
#   endif

    !DQEVAP=0.0!Added by B Clark so no evaporation 
       DEALLOCATE(RTEMP1,RTEMP2)
       IF(MSR)WRITE(IPT,101)'!  EVAPORATION AND PRECIPITATION READ : COMPLETE'
     END IF 

!
!--REPORT RESULTS--------------------------------------------------------------!
!

     IF(MSR)WRITE(IPT,*)'!'
     IF(MSR)WRITE(IPT,*    )'!  NON-UNIFORM METEO     :    SET'
     IF(HFX_TM%NTIMES > 0 )THEN
       CALL GETTIME(TSTRING,3600*INT(HFX_TM%TIMES(1)))
       IF(MSR)WRITE(IPT,102)'!  HEAT/RAD DATA BEGIN   :  ',TSTRING        
       CALL GETTIME(TSTRING,3600*INT(HFX_TM%TIMES(HFX_TM%NTIMES)))
       IF(MSR)WRITE(IPT,102)'!  HEAT/RAD DATA END     :  ',TSTRING
     END IF
     IF(WND_TM%NTIMES > 0)THEN
       CALL GETTIME(TSTRING,3600*INT(WND_TM%TIMES(1)))
       IF(MSR)WRITE(IPT,102)'!  WIND DATA BEGIN       :  ',TSTRING          
       CALL GETTIME(TSTRING,3600*INT(WND_TM%TIMES(WND_TM%NTIMES)))
       IF(MSR)WRITE(IPT,102)'!  WIND DATA END         :  ',TSTRING
     END IF
     IF(EVP_FLAG)THEN
       IF(EVP_TM%NTIMES > 0)THEN
         CALL GETTIME(TSTRING,3600*INT(EVP_TM%TIMES(1)))
         IF(MSR)WRITE(IPT,102)'!  EVAP/PREC DATA BEGIN       :  ',TSTRING          
         CALL GETTIME(TSTRING,3600*INT(EVP_TM%TIMES(EVP_TM%NTIMES)))
         IF(MSR)WRITE(IPT,102)'!  EVAP/PREC DATA END         :  ',TSTRING
       END IF
     END IF
     
     FTEMP1 = SUM(DHFLUX/FLOAT(NNode))/FLOAT(HFX_TM%NTIMES)
     FTEMP2 = SUM(DHSHORT)/FLOAT(NNode*HFX_TM%NTIMES)
     IF(SERIAL)THEN
       RBUF1 = FTEMP1
       RBUF2 = FTEMP2
     END IF
#  if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_REDUCE(FTEMP1,RBUF1,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_REDUCE(FTEMP2,RBUF2,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
#  endif
     IF(MSR)WRITE(IPT,101)'!  AVE HEAT FLUX         :',RBUF1/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,101)'!  AVE SHORT WAVE RAD    :',RBUF2/FLOAT(NPROCS)

!     FTEMP1 = SUM(DTX)/FLOAT(NNode*WND_TM%NTIMES)
!     FTEMP2 = SUM(DTY)/FLOAT(NNode*WND_TM%NTIMES)
     FTEMP1 = SUM(DTX)/FLOAT(MElem*WND_TM%NTIMES)
     FTEMP2 = SUM(DTY)/FLOAT(MElem*WND_TM%NTIMES)
     IF(SERIAL)THEN
       RBUF1 = FTEMP1
       RBUF2 = FTEMP2
     END IF
#  if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_REDUCE(FTEMP1,RBUF1,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_REDUCE(FTEMP2,RBUF2,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
#  endif
     IF(MSR)WRITE(IPT,101)'!  AVE WIND X-COMP       :',RBUF1/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,101)'!  AVE WIND Y-COMP       :',RBUF2/FLOAT(NPROCS)

     IF(EVP_FLAG)THEN
       FTEMP1 = SUM(DQEVAP)/FLOAT(NNode*EVP_TM%NTIMES)
       FTEMP2 = SUM(DQPREC)/FLOAT(NNode*EVP_TM%NTIMES)
       IF(SERIAL)THEN
         RBUF1 = FTEMP1
         RBUF2 = FTEMP2
       END IF
#  if defined (MULTIPROCESSOR)
       IF(PAR)CALL MPI_REDUCE(FTEMP1,RBUF1,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
       IF(PAR)CALL MPI_REDUCE(FTEMP2,RBUF2,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
#  endif
       IF(MSR)WRITE(IPT,101)'!  AVE EVAPORATION       :',RBUF1/FLOAT(NPROCS)
       IF(MSR)WRITE(IPT,101)'!  AVE PRECIPITATION     :',RBUF2/FLOAT(NPROCS)
     END IF
   ELSE
     WRITE(IPT,*)'==================ERROR=================================='
     WRITE(IPT,*)'M_TYPE NOT CORRECT, --->',M_TYPE
     WRITE(IPT,*)'MUST BE "uniform" or "non-uniform"'
     WRITE(IPT,*)'========================================================='
     CALL PSTOP
   END IF
      
  
!
!--Format Statements-----------------------------------------------------------!
!

   100  FORMAT(1X,A26,I6," =>",2X,4(I5,1H,))
   101  FORMAT(1X,A26,F10.4)  
   102  FORMAT(1X,A28,A13)  
   1000 FORMAT(A80)
   5000 FORMAT(8E14.5)

   RETURN
   END SUBROUTINE BCS_FORCE
!==============================================================================|
