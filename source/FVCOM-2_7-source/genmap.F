!==============================================================================!
!  SET UP GLOBAL <--> LOCAL NODE AND ELEMENT NUMBERING MAPS                    !
!  SET UP INTERPROCESSOR COMMUNICATION                                         !
!==============================================================================!

   SUBROUTINE GENMAP 
# if defined (MULTIPROCESSOR)

!==============================================================================!
   USE ALL_VARS
   USE CONTROL 
   USE MOD_PAR  
   IMPLICIT NONE
   INTEGER I,J,K,N1 ,N2,N3,IERR,CHECK,NGLOB,NN,I1,I2,I3,NSZE,ISND,IRCV,OWNER
   INTEGER, ALLOCATABLE :: NP(:,:),NTEMP(:),NTEMP2(:),NTEMP3(:,:),NTEMP4(:,:)
   INTEGER, ALLOCATABLE :: NODEMARK(:)
   INTEGER :: NNODES_MINE
   INTEGER IGLID,IC1,IC2,TERCV
   INTEGER STAT(MPI_STATUS_SIZE)
   LOGICAL :: NEW,INELEM,ISMYN
   REAL(SP), ALLOCATABLE :: TESTA(:,:),TESTB(:,:),TESTC(:,:)
   REAL(SP)  DIFF,VAL,TIME_PAR,NEAVE,NNAVE
   REAL(DP) :: FCHECK
   INTRINSIC SYSTEM_CLOCK


   IF(MSR)WRITE(IPT,*  )'!'
   IF(MSR)WRITE(IPT,*)'!            SETTING UP MAPPING/MESSAGE PASSING '
   IF(MSR)WRITE(IPT,*  )'!'

!==============================================================================|
!   COUNT NUMBER OF LOCAL DOMAIN ELEMENTS                        MElem             |
!   CREATE LOCAL--> GLOBAL ELEMENT MAP                           EGID(1:MElem)     |
!   EGID(I) = GLOBAL ELEMENT NUMBER OF LOCAL ELEMENT I                         |
!==============================================================================|
   ALLOCATE(NTEMP(MElemGL)) ; NTEMP = 0
   MElem = 0
   DO I=1,MElemGL
     IF(EL_PID(I) == MYID) THEN
       MElem = MElem + 1
       NTEMP(MElem) = I
     END IF
   END DO
   ALLOCATE(EGID(0:MElem)) ; EGID = 0
   EGID(1:MElem) = NTEMP(1:MElem)
   DEALLOCATE(NTEMP)
   IF(MSR)WRITE(IPT,*)  '!  ELEMENT MAP           :    COMPLETE'
!==============================================================================|
!   COUNT NUMBER OF LOCAL NODES                                  NNode         |
!   CREATE LOCAL--> GLOBAL NODE MAP                              NGID(1:NNode) |
!   NGID(I) = GLOBAL NODE NUMBER OF LOCAL NODE I                               |
!==============================================================================|
   ALLOCATE(NODEMARK(NNodeGL)) ; NODEMARK = 0
   ALLOCATE(NTEMP(NNodeGL))
   DO I=1,MElemGL
     IF(EL_PID(I) == MYID)THEN
       NODEMARK(NVG(I,1:3)) = 1
     END IF
   END DO

   NNode = 0
   DO I=1,NNodeGL
     IF(NODEMARK(I) == 1)THEN
       NNode = NNode + 1
       NTEMP(NNode) = I
     END IF
   END DO
!   DEALLOCATE(NODEMARK)


   ALLOCATE(NGID(NNode))
   NGID(1:NNode) = NTEMP(1:NNode)
   DEALLOCATE(NTEMP)
   
   IF(MSR)WRITE(IPT,*)  '!  NODE MAP              :    COMPLETE'

!==============================================================================|
!   DETERMINE NUMBER OF HALO ELEMENTS                            MHE           |
!   DETERMINE GLOBAL MAPPING OF HALO ELEMENTS                    HE_LST(MHE)   | 
!   DETERMINE LOCAL MAPPING OF HALO ELEMENTS                     HE_LOC(MHE)   | 
!==============================================================================|
   
   ALLOCATE(NTEMP(MElemGL))
   ALLOCATE(NTEMP2(MElemGL))
   MHE = 0
   DO I=1,MElemGL
     NNODES_MINE = SUM(NODEMARK(NVG(I,1:3)))
!     IF(NNODES_MINE == 1 .OR. NNODES_MINE ==2)THEN
     IF(NNODES_MINE > 0 .AND. EL_PID(I) /= MYID)THEN
       MHE         = MHE + 1
       NTEMP(MHE)  = I
       NTEMP2(MHE) = EL_PID(I)
     END IF
   END DO
   DEALLOCATE(NODEMARK)

   ALLOCATE(HE_LST(MHE),HE_OWN(MHE))
   HE_LST(1:MHE) =  NTEMP(1:MHE)
   HE_OWN(1:MHE) = NTEMP2(1:MHE)
   DEALLOCATE(NTEMP,NTEMP2)
   IF(MSR)WRITE(IPT,*)  '!  HALO ELEMENTS         :    COMPLETE'

!==============================================================================|
!   DETERMINE LOCAL MAPPING OF GLOBAL ELMNTS (INTERNAL ONLY)     ELID          |
!   ELID(I) = LOCAL ELEMENT ID OF GLOBAL ELEMENT I                             |
!             RETURNS 0 IF I IS NOT AN INTERIOR ELEMENT                        |
!                                                                              |
!   DETERMINE LOCAL MAPPING OF GLOBAL ELMNTS (INTERNAL + HALO)   ELID_X        |
!   ELID_X(I) = LOCAL ELEMENT ID OF GLOBAL ELEMENT I                           |
!             RETURNS 0 IF I IS NOT AN INTERIOR OR HALO ELEMENT                |
!==============================================================================|

   ALLOCATE(ELID(0:MElemGL)) ; ELID = 0
   ALLOCATE(ELID_X(0:MElemGL)) ; ELID = 0

   DO I=1,MElem
     ELID(EGID(I)) = I
   END DO
   ELID_X = ELID

   DO I=1,MHE
     ELID_X(HE_LST(I)) = I+MElem
   END DO
   
   
!==============================================================================|
!   DETERMINE NUMBER OF INTERNAL BOUNDARY NODES             NBN                |
!   DETERMINE GLOBAL MAPPING OF INTERNAL BOUNDARY NODES     BN_LST(NBN)        | 
!   DETERMINE LOCAL MAPPING OF INTERNAL BOUNDARY NODES      BN_LOC(NBN)        |
!   DETERMINE MULTIPLICITY (HOW MANY OWNERS) OF IBNs        BN_MLT(NBN)        |
!   MARK OTHER OWNERS OF IBNs                               BN_NEY(NBN,NPROCS) | 
!==============================================================================|

   ALLOCATE(NTEMP(NNode)) ; NTEMP = 0
   ALLOCATE(NTEMP3(NNode,NPROCS)) ; NTEMP3 = 0
   ALLOCATE(NDE_ID(NNode)) ; NDE_ID = 0

!!$   DO I=1,NNode  
!!$     DO J=1,MElemGL
!!$       NGLOB = NGID(I)
!!$       FCHECK= DBLE(NVG(J,1)-NGLOB)*DBLE(NVG(J,2)-NGLOB)*DBLE(NVG(J,3)-NGLOB)
!!$       IF(ABS(FCHECK) < 1.0 .AND. MYID /= EL_PID(J))THEN
!!$          NTEMP(I) = 1
!!$          NTEMP3(I,EL_PID(J)) = 1
!!$          NDE_ID(I) = 1
!!$       END IF
!!$     END DO
!!$   END DO

   DO I=1,NNode          !loop through all local nodes
      NGLOB = NGID(I)    !get global ID of node
      DO K=1,MHE         !loop number of halo elements 
         J = HE_LST(K)   !get the halo element k's global element number J
		 
         FCHECK= DBLE(NVG(J,1)-NGLOB)*DBLE(NVG(J,2)-NGLOB)*DBLE(NVG(J,3)-NGLOB)
		                 !if any of global node of halo global element J matches global node NGLOB, then this 
		                 !is a boundary node (this nodes belong to one of the halo elements)
		 
         IF(ABS(FCHECK) < 1.0)THEN
            NTEMP(I) = 1
            NTEMP3(I,EL_PID(J)) = 1
            NDE_ID(I) = 1
         END IF
      END DO
   END DO

   NBN = SUM(NTEMP)  !total number of boundary nodes (all nodes that involve in halo elements)
   ALLOCATE(BN_LST(NBN))
   ALLOCATE(BN_MLT(NBN))
   ALLOCATE(BN_LOC(NBN))
   ALLOCATE(BN_NEY(NBN,NPROCS))
   BN_NEY = 0 ; BN_MLT = 0 ; BN_LST = 0
   NN = 1
   DO I=1,NNode
     IF(NTEMP(I)==1)THEN      !if this is boundary node
       BN_LST(NN) = NGID(I)   !record the boundary node's global node ID in BN_LST
       BN_LOC(NN) = I         !record the boundary node's local node ID in BN_LOC
       BN_NEY(NN,MYID) = 1    !mark other owners of this boundary node
       DO J=1,NPROCS
         IF(NTEMP3(I,J)==1) BN_NEY(NN,J) = 1
       END DO
       BN_MLT(NN) = SUM(NTEMP3(I,:))
       NN = NN + 1
     END IF
   END DO
   BN_MLT = BN_MLT + 1         !number of boundary node owners for boundary node (include self)
   MX_MLT = 1 
   IF(NBN > 0)MX_MLT = MAXVAL(BN_MLT)  !maximum number of owners for all boundary nodes for this processor
                                       !not all boundary nodes ahve the same number of owners
   
   DEALLOCATE(NTEMP,NTEMP3)
   IF(MSR)WRITE(IPT,*)  '!  BOUNDARY NODES        :    COMPLETE'
   
!==============================================================================|
!   DETERMINE NUMBER OF HALO NODES                               NHN           |
!   DETERMINE GLOBAL NUMBERING OF HALO NODES                     HN_LST(NHN)   |
!   DETERMINE PRIMARY OWNER OF HALO NODE                         HN_OWN(NHN)   !
!==============================================================================|

   ALLOCATE(NTEMP(NNodeGL)) ; NTEMP = 0
   NHN = 0
   DO I=1,NNodeGL
     DO J=1,MHE                  !loop through all halo elements
       I1 = NVG(HE_LST(J),1)     !fetch global node numbers of halo elements
       I2 = NVG(HE_LST(J),2)     
       I3 = NVG(HE_LST(J),3)  
       IF(I1 == I .OR. I2 == I .OR. I3 == I)THEN  !any node of the halo element which is equal to this global node I
         INELEM = .TRUE.
         ISMYN  = .FALSE.
         DO K=1,NNode                      !loop through local nodes
           NGLOB=NGID(K)                   !get global id of local nodes
           IF(I == NGLOB) ISMYN = .TRUE.   !if I is also a local node
         END DO
         IF(INELEM .AND. .NOT.ISMYN)THEN   !exclude local nodes (this global node is part of a halo element
           NHN = NHN + 1                   !                     yet not part of local nodes)
           NTEMP(NHN) = I
           EXIT
         END IF
       END IF
     END DO
   END DO
   ALLOCATE(HN_LST(NHN)) ; HN_LST = NTEMP(1:NHN)
   DEALLOCATE(NTEMP)

   ALLOCATE(HN_OWN(NHN))
   HN_OWN = 0
   DO I=1,NHN
     IGLID = HN_LST(I)
     DO J=1,MHE
       I1 = NVG(HE_LST(J),1)  
       I2 = NVG(HE_LST(J),2)  
       I3 = NVG(HE_LST(J),3)  
       IF(I1 == IGLID .OR. I2 == IGLID .OR. I3 == IGLID)THEN
         HN_OWN(I) = HE_OWN(J)
         EXIT
       END IF
     END DO
     IF(HN_OWN(I) == 0)THEN
       WRITE(IPT,*)'NO OWNER FOUND FOR HALO NODE',I,'OF DOMAIN',MYID
     END IF
   END DO

   IF(MSR)WRITE(IPT,*)  '!  HALO NODES            :    COMPLETE'

!==============================================================================|
!   DETERMINE LOCAL MAPPING OF GLOBAL NODES  (INTERNAL ONLY)     NLID          |
!   DETERMINE LOCAL MAPPING OF GLOBAL NODES  (INTERNAL + HALO)   NLID_X        |
!==============================================================================|

   ALLOCATE(NLID(0:NNodeGL)) ; NLID = 0
   ALLOCATE(NLID_X(0:NNodeGL)) ; NLID_X = 0

   DO I=1,NNode
     NLID(NGID(I)) = I
   END DO
   NLID_X = NLID

   DO I=1,NHN
     NLID_X(HN_LST(I)) = I+NNode
   END DO
   
         
!==============================================================================|
!   SEND INFORMATION TO PROCESSOR ZERO FOR STATISTICS OUTPUT                   |
!   PNE  :: NUMBER OF ELEMENTS IN EACH DOMAIN                                  |
!   PNN  :: NUMBER OF NODES IN EACH DOMAIN                                     |
!   PNHE :: NUMBER OF HALO ELEMENTS IN EACH DOMAIN                             |
!   PNBN :: NUMBER OF BOUNDARY NODES IN EACH DOMAIN                            |
!   PMBM :: MAXIMUM MULTIPLICITY OF BOUNDARY NODES IN EACH DOMAIN              |
!   PNHN :: NUMBER OF HALO NODES IN EACH DOMAIN                                |
!==============================================================================|
   
   ALLOCATE(PNE(NPROCS))
   ALLOCATE(PNN(NPROCS))
   ALLOCATE(PNHE(NPROCS))
   ALLOCATE(PNBN(NPROCS))
   ALLOCATE(PMBM(NPROCS))
   ALLOCATE(PNHN(NPROCS))
   
   CALL MPI_GATHER(MElem,1,MPI_INTEGER,PNE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(NNode,1,MPI_INTEGER,PNN,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   IF(MSR)THEN
   WRITE(IPT,*)'         PROC     #ELEMENTS   IMBALANCE      #NODES  IMBALANCE'
   NEAVE = FLOAT(MElemGL)/FLOAT(NPROCS)
   NNAVE = FLOAT(SUM(PNN(1:NPROCS)))/FLOAT(NPROCS)
   DO I=1,NPROCS
     WRITE(IPT,'(I12,I14,F12.6,I14,F12.6)')I,PNE(I),FLOAT(PNE(I))/NEAVE,PNN(I),FLOAT(PNN(I))/NNAVE
   END DO
   END IF
   CALL MPI_GATHER(MHE,1,MPI_INTEGER,PNHE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(NBN,1,MPI_INTEGER,PNBN,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(MX_MLT,1,MPI_INTEGER,PMBM,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(NHN,1,MPI_INTEGER,PNHN,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   IF(MSR)WRITE(IPT,*)'         PROC     #HALO ELMNTS   #BNDRY NODES  MAX MLTPLCTY  #HALO NODES'
   DO I=1,NPROCS
     IF(MSR)WRITE(IPT,'(I12,4I14)')I,PNHE(I),PNBN(I),PMBM(I),PNHN(I)
   END DO


!==============================================================================|
!   SET UP ELEMENT MAPPING FOR GLOBAL REASSEMBLY OF ARRAYS (MASTER ONLY)       | 
!   ELEMENT MAP :: EMAP(NPROCS)                                                |
!     EMAP(1-->NPROCS)%NSIZE  :: NUMBER OF ELEMENTS IN EACH DOM                |
!     EMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM  |
!==============================================================================|

   ALLOCATE(EMAP(NPROCS))  ; EMAP(:)%NSIZE = 0
   IF(MSR)THEN

!--Determine Number of Elements for Each Processor
   DO I=1,MElemGL    
     EMAP(EL_PID(I))%NSIZE = EMAP(EL_PID(I))%NSIZE + 1
   END DO

!--Allocate Mapping Array for Each Processor
   DO I=1,NPROCS
     ALLOCATE(EMAP(I)%LOC_2_GL(0:EMAP(I)%NSIZE))
     EMAP(I)%LOC_2_GL(0) = 0
   END DO
   
!--Construct Mapping Array for Each Processor 
   ALLOCATE(NTEMP(NPROCS))
   NTEMP = 0
   DO I=1,MElemGL
     I1 = EL_PID(I)
     NTEMP(I1) = NTEMP(I1) + 1
     EMAP(I1)%LOC_2_GL(NTEMP(I1)) = I
   END DO
   DEALLOCATE(NTEMP)
   END IF

!==============================================================================|
!   TEST ELEMENT MAPPING USING GATHER                                          | 
!==============================================================================|
   
   IF(MSR)THEN
     ALLOCATE(TESTA(MElem,KB))
     ALLOCATE(TESTB(MElemGL,KB))
   ELSE
     ALLOCATE(TESTA(MElem,KB))
     ALLOCATE(TESTB(1,1))
   END IF

   TESTA = FLOAT(MYID)
   CALL GATHER(LBOUND(TESTA,1),UBOUND(TESTA,1),MElem,MElemGL,KB,MYID,NPROCS,EMAP,TESTA,TESTB) 
   IF(MSR)THEN
     DO J=1,KB
       DO I=1,MElemGL
         IF(TESTB(I,J) /= EL_PID(I))THEN
           WRITE(IPT,*)'GLOBAL ARRAY NOT CORRECT',I,J,TESTB(I,J),EL_PID(I)
         END IF
       END DO
     END DO
   END IF
   DEALLOCATE(TESTA,TESTB)

    
   IF(MSR)WRITE(IPT,*)  '!  ELEMENT GATHER TEST   :    PASSED'
       
   
!==============================================================================|
!   SET UP NODE MAPPING FOR GLOBAL REASSEMBLY OF ARRAYS (MASTER ONLY)          | 
!   NODAL   MAP :: NMAP(NPROCS)                                                |
!     NMAP(1-->NPROCS)%NSIZE  :: NUMBER OF NODES IN EACH DOM                   |
!     NMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM  |
!==============================================================================|

   ALLOCATE(NMAP(NPROCS))  ; NMAP(:)%NSIZE = 0
   IF(MSR)THEN
   NMAP(:)%NSIZE = PNN(:)

   DO I=1,NPROCS
     ALLOCATE(NMAP(I)%LOC_2_GL(0:NMAP(I)%NSIZE))
     NMAP(I)%LOC_2_GL(0) = 0
   END DO
   END IF
  
   IF(.NOT. MSR)THEN
     ISND = MYID+1000
     CALL MPI_SEND(NGID,NNode,MPI_INTEGER,0,ISND,MPI_COMM_WORLD,IERR)
   END IF

   IF(MSR)THEN
     NMAP(1)%LOC_2_GL(1:NNode) = NGID(1:NNode)
     DO I=2,NPROCS
       NSZE = NMAP(I)%NSIZE
       ALLOCATE(NTEMP(NSZE))
       IRCV = I+1000
       CALL MPI_RECV(NTEMP,NSZE,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
       NMAP(I)%LOC_2_GL(1:NSZE) = NTEMP
       DEALLOCATE(NTEMP)
     END DO
   END IF

!==============================================================================|
!   TEST NODE MAPPING USING GATHER                                             | 
!==============================================================================|

   ALLOCATE(TESTA(NNode,KB))
   ALLOCATE(TESTB(NNodeGL,KB))

   TESTA = 100.0_SP
   CALL GATHER(LBOUND(TESTA,1),UBOUND(TESTA,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,TESTA,TESTB)
   IF(MSR)THEN
     DO J=1,KB
       DO I=1,NNodeGL
         IF(TESTB(I,J) /= 100.0_SP)THEN
           WRITE(IPT,*)'GLOBAL ARRAY NOT CORRECT',I,J,TESTB(I,J)
         END IF
       END DO
     END DO
   END IF
   DEALLOCATE(TESTA,TESTB)

   IF(MSR)WRITE(IPT,*)  '!  NODE GATHER TEST      :    PASSED'
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)



!==============================================================================|
!   SET UP COMMUNICATION: HALO ELEMENTS                                        |
!==============================================================================|

   ALLOCATE(EC(NPROCS))


   EC(:)%NRCV = 0
   ALLOCATE(NTEMP3(NPROCS,MElemGL))

!--Set up Recieve Information First
   DO I=1,MHE
     OWNER = HE_OWN(I)
     EC(OWNER)%NRCV = EC(OWNER)%NRCV + 1
     NTEMP3(OWNER,EC(OWNER)%NRCV) = I + MElem
   END DO
   
   NN = 0
   DO I=1,NPROCS
     EC(I)%RCPT = NN
     NSZE = EC(I)%NRCV 
     NN = NN + NSZE
!     ALLOCATE(EC(I)%RCVP(NSZE))
     IF(NSZE > 0)ALLOCATE(EC(I)%RCVP(NSZE))
     IF(NSZE == 0)ALLOCATE(EC(I)%RCVP(1))
     IF(NSZE > 0)EC(I)%RCVP = NTEMP3(I,1:NSZE)
   END DO
   DEALLOCATE(NTEMP3)
   TERCV = SUM(EC(:)%NRCV)
   IF(TERCV /= MHE)THEN
     WRITE(IPT,*)'TOTAL NUMBER OF ELEMENTS SET UP TO RECEIVE DATA'
     WRITE(IPT,*)'NOT EQUAL TO TOTAL NUMBER OF HALO ELEMENTS'
     WRITE(IPT,*)TERCV,MHE
     CALL MPI_FINALIZE(IERR)
     CALL PSTOP
   END IF

!--Set up Send
   EC(:)%NSND = 0

   DO I=1,NPROCS
   IF(MYID /=I)THEN
     ISND = MYID + 1000*I
     CALL MPI_SEND(EC(I)%NRCV,1,MPI_INTEGER,I-1,ISND,MPI_COMM_WORLD,IERR)
   END IF
   END DO
   DO I=1,NPROCS
   IF(MYID /=I)THEN
     IRCV = 1000*MYID + I
     CALL MPI_RECV(EC(I)%NSND,1,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
   END IF
!   ALLOCATE(EC(I)%SNDP(EC(I)%NSND))
   IF(EC(I)%NSND > 0)ALLOCATE(EC(I)%SNDP(EC(I)%NSND))
   IF(EC(I)%NSND == 0)ALLOCATE(EC(I)%SNDP(1))
   END DO
   DO I=1,NPROCS
   IF(MYID /=I .AND. EC(I)%NRCV > 0)THEN
     ALLOCATE(NTEMP(EC(I)%NRCV))
     DO J=1,EC(I)%NRCV
       NTEMP(J) = HE_LST(EC(I)%RCVP(J)-MElem)
     END DO
     ISND = MYID + 2000*I
     CALL MPI_SEND(NTEMP,EC(I)%NRCV,MPI_INTEGER,I-1,ISND,MPI_COMM_WORLD,IERR)
     DEALLOCATE(NTEMP)
   END IF
   END DO
   DO I=1,NPROCS
   IF(MYID /=I .AND. EC(I)%NSND > 0)THEN
     ALLOCATE(NTEMP(EC(I)%NSND))
     IRCV = 2000*MYID + I
     CALL MPI_RECV(NTEMP,EC(I)%NSND,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
     N1 = 0
     DO J=1,EC(I)%NSND
     DO K=1,MElem
       IF(EGID(K)==NTEMP(J))THEN
         N1 = N1 + 1
         EC(I)%SNDP(N1) = K
       END IF
     END DO
     END DO
     DEALLOCATE(NTEMP)
   END IF
   END DO
   

!==============================================================================|
!   CHECK ELEMENT COMMUNICATION                                                |
!==============================================================================|
!--Set Up Testing Array 
   ALLOCATE(TESTA(0:MElem+MHE,KB))  ; TESTA = 0.0_SP
   ALLOCATE(TESTB(0:MElem+MHE,KB))  ; TESTB = 0.0_SP
   ALLOCATE(TESTC(0:MElem+MHE,KB))  ; TESTC = 0.0_SP

!--Initialize with Function of Global Element Index in Interior, Zero in Halo
   DO I=1,KB
   TESTA(1:MElem,I) =  (FLOAT(EGID(1:MElem))  +10000*I)
   TESTB(1:MElem,I) = -(FLOAT(EGID(1:MElem))  +10000*I)
   TESTC(1:MElem,I) =  (FLOAT(EGID(1:MElem))  +10000*I) + .5_SP 
   END DO
   

!--Perform Communication

   CALL EXCHANGE(EC,MElem+MHE,KB,MYID,NPROCS,TESTA,TESTB,TESTC) 

!--Check Results For Interior on Array A
   DO I=1,MElem
   DO J=1,KB
   VAL  = FLOAT(EGID(I)+10000*J)  
   DIFF = ABS(TESTA(I,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'INTERNAL ELEMENT VALUE CHANGED on A',I,J,TESTA(I,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO

!--Check Results For Halo on Array A 
   DO I=1,MHE
   DO J=1,KB
   VAL = FLOAT(HE_LST(I)+10000*J)
   DIFF = ABS(TESTA(I+MElem,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'HALO ELEMENT CHANGED on A',I,J,TESTA(I+MElem,J),VAL,I+MElem
   CALL PSTOP
   END IF
   END DO
   END DO
   

!--Check Results For Interior on Array B
   DO I=1,MElem
   DO J=1,KB
   VAL  = -FLOAT(EGID(I)+10000*J)  
   DIFF = ABS(TESTB(I,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'INTERNAL ELEMENT VALUE CHANGED on B',I,J,TESTB(I,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO

!--Check Results For Halo on Array B 
   DO I=1,MHE
   DO J=1,KB
   VAL = -FLOAT(HE_LST(I)+10000*J)
   DIFF = ABS(TESTB(I+MElem,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'HALO ELEMENT CHANGED on B',I,J,TESTB(I+MElem,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO

!--Check Results For Interior on Array C
   DO I=1,MElem
   DO J=1,KB
   VAL  = FLOAT(EGID(I)+10000*J)+.5_SP
   DIFF = ABS(TESTC(I,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'INTERNAL ELEMENT VALUE CHANGED on C',I,J,TESTC(I,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO

!--Check Results For Halo on Array C 
   DO I=1,MHE
   DO J=1,KB
   VAL = FLOAT(HE_LST(I)+10000*J)+.5_SP
   DIFF = ABS(TESTC(I+MElem,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'HALO ELEMENT CHANGED on C',I,J,TESTC(I+MElem,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO
  
   DEALLOCATE(TESTA,TESTB,TESTC)
   
   IF(MSR)WRITE(IPT,*)  '!  ELEMENT COMM TEST     :    PASSED'
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

!==============================================================================|
!   SET UP COMMUNICATION: BOUNDARY NODES                                       |
!==============================================================================|

   ALLOCATE(NTEMP3(NBN,NPROCS),NTEMP4(NBN,NPROCS))
   ALLOCATE(BNC(NPROCS)) 
   BNC%NRCV = 0      ; BNC%NSND = 0

!--Count Receives from Each Processor
   DO I=1,NBN
     DO J=1,NPROCS
       IF(BN_NEY(I,J)==1 .AND. MYID /= J) THEN
         BNC(J)%NRCV = BNC(J)%NRCV + 1
         NTEMP3(BNC(J)%NRCV,J) = BN_LOC(I) 
         NTEMP4(BNC(J)%NRCV,J) = BN_MLT(I) 
       END IF
     END DO
   END DO

!--Sent up Indices for Receive Buffer
   NN = 0
   DO I=1,NPROCS
     BNC(I)%RCPT = NN
     NSZE = BNC(I)%NRCV 
     NN = NN + NSZE
   END DO

!--Sends = Receives for Boundary Nodes
   BNC%NSND = BNC%NRCV

!--Array of Receives/Sends/Multiplicities for Each Processor
   DO I=1,NPROCS
     NSZE = BNC(I)%NRCV
     ALLOCATE(BNC(I)%RCVP(NSZE))
     ALLOCATE(BNC(I)%SNDP(NSZE))
     ALLOCATE(BNC(I)%MLTP(NSZE))
     BNC(I)%RCVP = NTEMP3(1:NSZE,I)
     BNC(I)%SNDP = NTEMP3(1:NSZE,I)
     BNC(I)%MLTP = NTEMP4(1:NSZE,I)
   END DO
   DEALLOCATE(NTEMP3)


!--Sort Send/Recive Arrays According to Global Identity to Match Order    
   
   DO I=1,NPROCS
     IF(BNC(I)%NRCV > 0)THEN
       NSZE     = BNC(I)%NRCV
       ALLOCATE(NTEMP(NSZE),NTEMP2(NSZE))
       NTEMP(:) = NGID(BNC(I)%RCVP(:))
       CALL SORT(NTEMP,NTEMP2,NSZE)
       NTEMP = BNC(I)%RCVP
       DO J=1,NSZE
         BNC(I)%RCVP(J) = NTEMP(NTEMP2(J))
       END DO
       NTEMP = BNC(I)%MLTP
       DO J=1,NSZE
         BNC(I)%MLTP(J) = NTEMP(NTEMP2(J))
       END DO
       BNC(I)%SNDP = BNC(I)%RCVP
       DEALLOCATE(NTEMP,NTEMP2)
     END IF
   END DO
 
 
!==============================================================================|
!   TEST COMMUNICATION: HALO NODES                                             |
!==============================================================================|

!--Set Up Testing Array 
   ALLOCATE(TESTA(0:NNode+NHN,KB))  
   ALLOCATE(TESTB(0:NNode+NHN,KB))  
   ALLOCATE(TESTC(0:NNode+NHN,KB))  
   TESTA = FLOAT(MYID)  
   TESTB = FLOAT(MYID)  * 100.0_SP
   TESTC = FLOAT(MYID)  * 10000.0_SP
   DO J=1,KB
   DO I=1,NBN
     TESTA(BN_LOC(I),J) = J
     TESTB(BN_LOC(I),J) = 100*J
     TESTC(BN_LOC(I),J) = 10000*J
   END DO
   END DO

!--Perform Communication
   
   CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NNode+NHN,KB,MYID,NPROCS,TESTA,TESTB,TESTC) 
!

!--Check Results For Interior and Boundary on Array A
   DO J=1,KB
   DO I=1,NNode
     IF(NDE_ID(I)==0)THEN
     IF( ABS(TESTA(I,J) -  FLOAT(MYID))> .01_SP )THEN
       WRITE(IPT,*)'INTERNAL NODE CHANGED',I,J,TESTA(I,J),FLOAT(MYID)
       CALL PSTOP
     END IF
     ELSE
     IF(ABS(TESTA(I,J) - J) > .01_SP)THEN
       WRITE(IPT,*)'BOUNDARY PROBLEM',MYID,I,J,TESTA(I,J)
       WRITE(*,*)'GLOBAL NODE: ',NGID(I)
       CALL PSTOP
     END IF
     END IF
   END DO
   END DO
   
!--Check Results For Interior and Boundary on Array B
   DO J=1,KB
   DO I=1,NNode
     IF(NDE_ID(I)==0)THEN
     IF( ABS(TESTB(I,J) -  FLOAT(MYID)*100.0_SP)> .01_SP )THEN
       WRITE(IPT,*)'INTERNAL NODE CHANGED',I,J,TESTB(I,J),FLOAT(MYID)
       CALL PSTOP
     END IF
     ELSE
     IF(ABS(TESTB(I,J) - 100*J) > .01_SP)THEN
       WRITE(IPT,*)'BOUNDARY PROBLEM',MYID,I,J,TESTB(I,J)
       CALL PSTOP
     END IF
     END IF
   END DO
   END DO

!--Check Results For Interior and Boundary on Array C
   DO J=1,KB
   DO I=1,NNode
     IF(NDE_ID(I)==0)THEN
     IF( ABS(TESTC(I,J) -  FLOAT(MYID)*10000.0_SP)> .01_SP )THEN
       WRITE(IPT,*)'INTERNAL NODE CHANGED',I,J,TESTC(I,J),FLOAT(MYID)
       CALL PSTOP
     END IF
     ELSE
     IF(ABS(TESTC(I,J) - 10000.0_SP*J) > .01_SP)THEN
       WRITE(IPT,*)'BOUNDARY PROBLEM',MYID,I,J,TESTC(I,J)
       CALL PSTOP
     END IF
     END IF
   END DO
   END DO

   DEALLOCATE(TESTA,TESTB,TESTC)


   IF(MSR)WRITE(IPT,*)'!  BNDRY NODE COMM TEST  :    PASSED'

!==============================================================================|
!   SET UP COMMUNICATION: HALO NODES                                           |
!==============================================================================|
! NHN,HN_LST --HN_OWN (owner of halo node i, like he_own)

   ALLOCATE(NC(NPROCS))
   NC(:)%NRCV = 0

   ALLOCATE(NTEMP3(NPROCS,NNodeGL))

!--Set up Recieve Information First
   DO I=1,NHN
     OWNER = HN_OWN(I)
     NC(OWNER)%NRCV = NC(OWNER)%NRCV + 1
     NTEMP3(OWNER,NC(OWNER)%NRCV) = I + NNode
   END DO
   
   NN = 0
   DO I=1,NPROCS
     NC(I)%RCPT = NN
     NSZE = NC(I)%NRCV 
     NN = NN + NSZE
     ALLOCATE(NC(I)%RCVP(NSZE))
     NC(I)%RCVP = NTEMP3(I,1:NSZE)
   END DO
   DEALLOCATE(NTEMP3)
   TERCV = SUM(NC(:)%NRCV)
   IF(TERCV /= NHN)THEN
     WRITE(IPT,*)'TOTAL NUMBER OF NODES SET UP TO RECEIVE DATA'
     WRITE(IPT,*)'NOT EQUAL TO TOTAL NUMBER OF HALO NODES'
     WRITE(IPT,*)TERCV,NHN
     CALL MPI_FINALIZE(IERR)
     CALL PSTOP
   END IF

!--Set up Send
   NC(:)%NSND = 0

   DO I=1,NPROCS
   IF(MYID /=I)THEN
     ISND = MYID + 1000*I
     CALL MPI_SEND(NC(I)%NRCV,1,MPI_INTEGER,I-1,ISND,MPI_COMM_WORLD,IERR)
   END IF
   END DO
   DO I=1,NPROCS
   IF(MYID /=I)THEN
     IRCV = 1000*MYID + I
     CALL MPI_RECV(NC(I)%NSND,1,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
   END IF
   ALLOCATE(NC(I)%SNDP(NC(I)%NSND))
   END DO

   DO I=1,NPROCS
   IF(MYID /=I .AND. NC(I)%NRCV > 0)THEN
     ALLOCATE(NTEMP(NC(I)%NRCV))
     DO J=1,NC(I)%NRCV
       NTEMP(J) = HN_LST(NC(I)%RCVP(J)-NNode)
     END DO
     ISND = MYID + 1000*I
     CALL MPI_SEND(NTEMP,NC(I)%NRCV,MPI_INTEGER,I-1,ISND,MPI_COMM_WORLD,IERR)
     DEALLOCATE(NTEMP)
   END IF
   END DO

   DO I=1,NPROCS
   IF(MYID /=I .AND. NC(I)%NSND > 0)THEN
     ALLOCATE(NTEMP(NC(I)%NSND))
     IRCV = 1000*MYID + I
     CALL MPI_RECV(NTEMP,NC(I)%NSND,MPI_INTEGER,I-1,IRCV,MPI_COMM_WORLD,STAT,IERR)
     N1 = 0
     DO J=1,NC(I)%NSND
     DO K=1,NNode
       IF(NGID(K)==NTEMP(J))THEN
         N1 = N1 + 1
         NC(I)%SNDP(N1) = K
       END IF
     END DO
     END DO
     DEALLOCATE(NTEMP)
   END IF
   END DO
   
!==============================================================================|
!   CHECK HALO NODE COMMUNICATION                                              |
!==============================================================================|

!--Set Up Testing Array 
   ALLOCATE(TESTA(0:NNode+NHN,KB))  ; TESTA = 0.0_SP
   ALLOCATE(TESTB(0:NNode+NHN,KB))  ; TESTB = 0.0_SP
   ALLOCATE(TESTC(0:NNode+NHN,KB))  ; TESTC = 0.0_SP

!--Initialize with Global Nodal ID in Interior, Zero in Halo
   DO I=1,KB
   TESTA(1:NNode,I) =   FLOAT(NGID)  +10000*I
   TESTB(1:NNode,I) = -(FLOAT(NGID)  +10000*I)
   TESTC(1:NNode,I) =  (FLOAT(NGID)  +10000*I) + .50_SP
   END DO

!--Perform Communication
   CALL EXCHANGE(NC,NNode+NHN,KB,MYID,NPROCS,TESTA,TESTB,TESTC) 
!--Check Results For Interior on Array A
   DO I=1,NNode
   DO J=1,KB
   VAL  = FLOAT(NGID(I)+10000*J)  
   DIFF = ABS(TESTA(I,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'INTERNAL NODE VALUE CHANGED on A',I,J,TESTA(I,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO

!--Check Results For Halo on Array A 
   DO I=1,NHN
   DO J=1,KB
   VAL = FLOAT(HN_LST(I)+10000*J)
   DIFF = ABS(TESTA(I+NNode,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'HALO NODE CHANGED on A',I,J,TESTA(I+NNode,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO
   

!--Check Results For Interior on Array B
   DO I=1,NNode
   DO J=1,KB
   VAL  = -FLOAT(NGID(I)+10000*J)  
   DIFF = ABS(TESTB(I,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'INTERNAL NODE VALUE CHANGED on B',I,J,TESTB(I,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO

!--Check Results For Halo on Array B 
   DO I=1,NHN
   DO J=1,KB
   VAL = -FLOAT(HN_LST(I)+10000*J)
   DIFF = ABS(TESTB(I+NNode,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'HALO NODE CHANGED on B',I,J,TESTB(I+NNode,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO

!--Check Results For Interior on Array C
   DO I=1,NNode
   DO J=1,KB
   VAL  = FLOAT(NGID(I)+10000*J)+.5_SP
   DIFF = ABS(TESTC(I,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'INTERNAL NODE VALUE CHANGED on C',I,J,TESTC(I,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO

!--Check Results For Halo on Array C 
   DO I=1,NHN
   DO J=1,KB
   VAL = FLOAT(HN_LST(I)+10000*J)+.5_SP
   DIFF = ABS(TESTC(I+NNode,J)-VAL)
   IF(DIFF > .01_SP)THEN
   PRINT *,MYID,'HALO ELEMENT CHANGED on C',I,J,TESTC(I+NNode,J),VAL
   CALL PSTOP
   END IF
   END DO
   END DO
  
   DEALLOCATE(TESTA,TESTB,TESTC)
   
   IF(MSR)WRITE(IPT,*)  '!  NODE COMM TEST        :    PASSED'

!==============================================================================|
!   SET TOTAL FOR EACH PROCESSOR                                               |
!   MTElem:   TOTAL NUMBER OF ELEMENTS = INTERIOR (MElem) + HALO (MHE)                 |
!   NTNode:   TOTAL NUMBER OF NODES    = INTERIOR (NNode) + HALO (NHN)                 |
!==============================================================================|
   MTElem = MElem + MHE
   NTNode = NNode + NHN
  
   RETURN
# endif
   END SUBROUTINE GENMAP
!==============================================================================!


