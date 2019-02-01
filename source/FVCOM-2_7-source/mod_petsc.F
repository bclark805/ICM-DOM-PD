MODULE MOD_PETSC
# if defined (SEMI_IMPLICIT) && (MULTIPROCESSOR)
  IMPLICIT NONE

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  REAL*8, POINTER          :: XVALS_EL(:)

!======== PETSc VARIABLE BLOCK ==============================================
!common variables
!  PetscInt    :: IERR
!  PetscScalar :: NEG_ONE = -1.0D0
!  PetscScalar :: ONE     =  1.0D0
!  PetscScalar :: ZERO    =  0.0D0
!  PetscViewer :: viewer

!petsc 2d data structures
  Mat         :: A_EL
  Vec         :: X_EL
  Vec         :: B_EL
  Vec         :: XL_EL
  Vec         :: BL_EL
  PC          :: Pc_EL
  KSP         :: Ksp_EL
  PetscInt    :: N_VERTS
  PetscInt    :: N_HALO
  PetscInt    :: PSIZE_EL_HALO
  PetscReal   :: NORM_EL
  IS          :: ISLOCAL_EL
  IS          :: ISGLOBAL_EL
  VecScatter  :: L2G_EL
  VecScatter  :: G2L_EL
  PetscInt    :: ITS_EL
  ISLocalToGlobalMapping :: ISL2G_EL

!integer mappings
  INTEGER,ALLOCATABLE  :: PLO_2_ALO_NODE(:)  !maps PLO to ALO on a node (surface vertex) basis
  INTEGER,ALLOCATABLE  :: ALO_2_PLO_NODE(:)  !maps ALO to PLO on a node (surface vertex) basis

  INTEGER,PARAMETER    :: MAX_NZ_ROW = 30    !max nonzero values in any row

  INTEGER :: PUNIT
  INTEGER :: PSCRN
  LOGICAL :: PMSR  = .TRUE.                  !petsc master processor

!-------------------------------------------------------------------------------------
! The relative residual reduction required for convergence. this can be overriden
! at runtime from the command line.
!-------------------------------------------------------------------------------------
  PetscReal :: RTOL = 1.0E-10

!-------------------------------------------------------------------------------------
! AR Aspect ratio. In the fake poisson discretization, the d^2/dz^2 terms are set to
! this number instead of order 1.
!-------------------------------------------------------------------------------------
  PetscInt  :: AR = 1

!-------------------------------------------------------------------------------------
! CHECK_EXACT If this is  true, we set up a RHS using a solution vector u = 1 (Au=b)
! then we solve Ax=b and make sure x and u are the same (within tol).
!-------------------------------------------------------------------------------------
  LOGICAL   :: CHECK_EXACT = .FALSE.

!-------------------------------------------------------------------------------------
! If this is true. we set the x vector from qnh and this is used by petsc as an
! initial condition for the Solve.  If this is not true we skip setting the initial
! condition (to save time) and Petsc zeroes out the x vector on its own.
!-------------------------------------------------------------------------------------
  LOGICAL   :: USE_LAST = .TRUE.

!-------------------------------------------------------------------------------------
! If this is true, extra debug statements are printed to file unit numbers myid+800
! where myid is the processor id.
!-------------------------------------------------------------------------------------
  LOGICAL   :: PDEBUG   = .FALSE.

!============ END PETSc BLOCK ========================================================

  !----------------------------------------------------------------------------
  CONTAINS
  ! PETSc_SET
  ! PETSc_AllocA
  ! free surface operation
  !   --> PETSc_SetICS_EL    :
  !   --> PETSc_Solve_EL     :
  !
  ! PETSc_CleanUp            : clear petsc vars in memory
  ! PLO_2_ALO                : map [i] PLO to [i,k] ALO
  ! ALO_2_PLO                : map [i,k] ALO to [i] PLO
  !----------------------------------------------------------------------------

  SUBROUTINE PETSc_SET
  !============================================================================
  !  ALO:  Application Local Node Ordering
  !  PLO:  Petsc Local Node Ordering
  !  AGO:  Application Global Node Ordering
  !  PGO:  Petsc Global Node Ordering
  !============================================================================
  USE MOD_PAR,  ONLY : EL_PID,NLID,NLID_X
  USE ALL_VARS, ONLY : NVG, MSR, MElemGL, NNodeGL, NTNode, NPROCS, MYID, NNode, IPT
  IMPLICIT NONE

  PetscInt:: IERR
  INTEGER :: I,L,J,K,ICNT,PID,N1,PGL_NODE
  INTEGER :: MINDOMSIZE, MINDOM
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: CNT
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: ACCUM
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PARTID
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PART_NODES
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: OWNER_COUNT
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PLO_2_AGO_NODE
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: AGO_2_PGO_NODE
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PLO_2_PGO_NODE
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: OWNERS
  CHARACTER(LEN=20)                  :: SUBNAME = 'PETSc_SETUP'
  LOGICAL                            :: NEWOWN

  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PTMP
  !-----------------------------------------------------------------------------
  !setup Petsc vectors and matrices (first call only)
  !we will assume that the non-zero structure of matrix DOES NOT CHANGE
  !during iteration process.  This will allow us to preallocate properly
  !-----------------------------------------------------------------------------

  PSCRN = IPT             ! set petsc unit number for dumping regular information
  PUNIT = 800+MYID        ! set petsc unit number for dumping debug information
  PMSR = MSR              ! set petsc master (primarily for printing) to FVCOM master

  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ', TRIM(SUBNAME)
  IF(PMSR)   WRITE(PSCRN,*) 'PETSc SETUP        : BEGINNING'

  CALL PETScInitialize(PETSC_NULL_CHARACTER,IERR) ; CHKERRQ(IERR)
  IF(PMSR)   WRITE(PSCRN,*) 'PETSc INITIALIZE   : COMPLETE'

  !---------------------------------------------------------
  ! set up AO<-->PO mappings, sizes, etc.
  !---------------------------------------------------------

  !-----> simple partitioning does not give great load balance partition the domain by
  ! nodes such that each local node in Petsc Partition is also local in Application partition
  ! But, we do not have overlapping domains What were originally boundary nodes in FVCOM decomp
  ! and had multiple partition owners are assigned to the adjacent partition of lower number
  ! petsc partition id for each node in global application order ==> [partid(1:NNodeGL)]
  ! we should really consider a load balancing procedure where we first assign node that
  ! arent boundary nodes and then assign boundary nodes to the adjacent partition which is smaller.
  ! Ultimately, we want N_verts to be as equal as possible in each domain.  For right now,
  ! use this simple partitioning strategy

!  allocate(partid(NNodeGL)) ; partid = nprocs+1
!  do i=1,MElemGL
!  do l=1,3
!    partid(nvg(i,l)) = min( partid(nvg(i,l)) , el_pid(i))
!  end do
!  end do

  !-----> more complex partitioning for better load balance partition the domain by nodes
  ! such that each local node in Petsc Partition is also local in Application partition
  ! we will do two sweeps. In a first sweep, any node that was not an interprocessor boundary
  ! node will be assigned to the owner that had the node in the FVCOM partition.
  ! then we will count nodes in each Petsc Partition then we will loop over boundary nodes
  ! (which have multiple owners in the FVCOM partition) and decide who gets it based
  ! on who has the lowest current node count this should result in a nearly optimal loadbalance

  !--sweep 1, determine owners, owner count, assume max # owners = 10
  ALLOCATE(OWNER_COUNT(NNodeGL)) ; OWNER_COUNT = 0
  ALLOCATE(OWNERS(NNodeGL,0:10)) ; OWNERS      = 0
  DO I =1, MElemGL
    DO L =1, 3
      !check if owner is already in list --> can use loc!
      NEWOWN = .TRUE.
      DO J=1,OWNER_COUNT(NVG(I,L))
        IF(OWNERS(NVG(I,L),J) == EL_PID(I) ) NEWOWN = .FALSE.
      ENDDO
      IF(NEWOWN)THEN
        OWNER_COUNT(NVG(I,L)) = OWNER_COUNT(NVG(I,L)) + 1
        IF(OWNER_COUNT(NVG(I,L)) > 10) THEN
          WRITE(PSCRN,*) 'ERROR IN PETSc PARTITIONING'
          WRITE(PSCRN,*) 'INCREASE OWNER COUNT FROM 10'
          CALL MPI_FINALIZE(IERR)
          STOP
        ENDIF
        OWNERS( NVG(I,L),OWNER_COUNT(NVG(I,L)) ) = EL_PID(I)
      ENDIF
    ENDDO
  ENDDO

  !--now sweep and petsc partition the non-interprocessor boundary nodes
  ALLOCATE(PART_NODES(0:NPROCS)) ; PART_NODES = 0
  ALLOCATE(PARTID(NNodeGL)) ; PARTID = -1
  DO I= 1, NNodeGL
    IF(OWNER_COUNT(I) == 1) THEN
      PART_NODES(OWNERS(I,1)) = PART_NODES(OWNERS(I,1)) + 1
      PARTID(I) = OWNERS(I,1)
    ENDIF
  ENDDO

  !--now sweep boundary nodes, assigning a node to owner with least count
  DO I= 1, NNodeGL
    IF(OWNER_COUNT(I) > 1) THEN
      MINDOMSIZE = HUGE(I)
      DO J= 1, OWNER_COUNT(I)
        !find which adjacent dom has min count
        IF (PART_NODES(OWNERS(I,J)) < MINDOMSIZE) THEN
          MINDOMSIZE = PART_NODES(OWNERS(I,J) )
          MINDOM     = OWNERS(I,J)
        ENDIF
      ENDDO
      !assign boundary node to that partition
      PARTID(I) = MINDOM
      PART_NODES(MINDOM) = PART_NODES(MINDOM) + 1
    ENDIF
  ENDDO

  !--make sure everybody got assigned
  DO I = 1, NNodeGL
    IF(PARTID(I) == -1) THEN
      WRITE(PSCRN,*) 'ERROR IN PETSc PARTITIONING'
      WRITE(PSCRN,*) 'GLOBAL NODE: ',I,' NEVER GOT ASSIGNED A PETSc PARTITION!'
      WRITE(PSCRN,*) 'FIX IT!'
      CALL MPI_FINALIZE(IERR)
      STOP
    ENDIF
  ENDDO

  !--deallocate the data we used for 'fancy' partition
  DEALLOCATE(OWNERS)
  DEALLOCATE(OWNER_COUNT)

  ! lets count (or recount) nodes in each Petsc partition ==> [part_nodes(0:nprocs)]
  ! we need to start from zero so we can access partition (myid-1) later the total will equal NNodeGL,
  ! unlike the sum of nodes over FVCOM partitions which is > NNodeGL
  IF(.NOT.ALLOCATED(PART_NODES) ) ALLOCATE(PART_NODES(0:NPROCS))
  PART_NODES = 0
  DO I= 1, NNodeGL
    PART_NODES( PARTID(I)) = PART_NODES( PARTID(I) ) + 1
  ENDDO

  !count number of nodes in my local Petsc partition ==> [N_verts] halt if partition is empty
  N_VERTS = PART_NODES(MYID)

  print*,'N_VERTS and NNode=',N_VERTS,NNode,NTNode,MYID
  IF(N_VERTS == 0) THEN
    WRITE(PUNIT,*) 'MY PETSc PARTITION DOSE NOT CONTAIN ANY UNKNOWNS'
    CALL MPI_FINALIZE(IERR)
    STOP
  ENDIF
  IF(PDEBUG) WRITE(PUNIT,*) 'NODES IN PETSc PARTION: ',N_VERTS,' APPLICATION',NNode

  !report all sizes(master only)
  IF(PMSR) WRITE(PSCRN,*)'PARTITIONING       : COMPLETE'
  IF(PMSR) THEN
    DO I= 1, NPROCS
      WRITE(PSCRN,*)'PETSc PARTITION: ',I,'  SIZE  ',PART_NODES(I)
    ENDDO
  ENDIF

  !make sure total sum of nodes in all petsc partitions is equal to NNodeGL
  IF(SUM(PART_NODES) /= NNodeGL) THEN
    WRITE(PSCRN,*)'TOTAL SUM OF NODES IN PETSc PARTIONS DOES NOT EQUAL GLOBAL NUMBER OF NODES'
    WRITE(PSCRN,*)'EXITING!'
    CALL MPI_FINALIZE(IERR)
    STOP
  ENDIF

  !determine which global nodes (in application order) are in my Petsc partition
  !gives us a mapping (for interior nodes) of PLO to AGO
  !allocate it to size NNodeGL because we do not know size of petsc partition halo
  !gives us ==> [plo_2_ago_node(1:N_verts)]
  ALLOCATE(PLO_2_AGO_NODE(NNodeGL)) ; PLO_2_AGO_NODE = 0
  ICNT = 0
  DO I= 1, NNodeGL
    IF(PARTID(I) == MYID) THEN
      ICNT = ICNT + 1
      PLO_2_AGO_NODE(ICNT) = I
    ENDIF
  ENDDO

  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'PLO_2_AGO_NODE'
    DO I= 1, N_VERTS
      WRITE(PUNIT,*) I, PLO_2_AGO_NODE(I)
    ENDDO
  ENDIF

  !get a mapping from petsc local nodes to fvcom local nodes ==> plo_2_alo_node(1:N_verts)
  !we already have the global application order number of each node, we will
  !use nlid to convert this to local application order of each node
  !this will be used to in our plo_2_alo function which is ultimately used to
  !load and unload FVCOM vectors (q and source term)
  !gives us ==> [plo_2_alo_node(N_verts)]
  ALLOCATE(PLO_2_ALO_NODE(N_VERTS)); PLO_2_ALO_NODE = 0
  IF(NPROCS > 1) THEN
    DO I=1,N_VERTS
      PLO_2_ALO_NODE(I) = NLID(PLO_2_AGO_NODE(I))
    ENDDO
  ELSE
    PLO_2_ALO_NODE = PLO_2_AGO_NODE
  ENDIF
  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'PLO_2_ALO_NODE:'
    DO I= 1, N_VERTS
      WRITE(PUNIT,*) I, PLO_2_ALO_NODE(I)
    ENDDO
  ENDIF

  !get a mapping from application local nodes to petsc local nodes
  !first we map without setting petsc halos
  !gives us ==> [alo_2_plo_node(1:NNode)]
  ALLOCATE(ALO_2_PLO_NODE(NTNode)) ; ALO_2_PLO_NODE = 0
  IF(NPROCS > 1)THEN
    ICNT = 0
    DO I= 1, NNodeGL
      IF(NLID_X(I) /= 0 .AND. PARTID(I) == MYID) THEN
        ICNT = ICNT + 1
        ALO_2_PLO_NODE(NLID_X(I)) = ICNT
      ENDIF
    ENDDO
  ELSE
    ALO_2_PLO_NODE = PLO_2_ALO_NODE
  ENDIF
  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'ALO_2_PLO_NODE:'
    DO I= 1, N_VERTS
      WRITE(PUNIT,*) I, ALO_2_PLO_NODE(I)
    ENDDO
  ENDIF

  IF(ICNT /= N_VERTS)THEN
    WRITE(PUNIT,*) 'PROBLEM IN ALO to PLO MAPPING!'
    CALL MPI_FINALIZE(IERR)
    STOP
  ENDIF

  !we continue mapping ALO to PLO but now we accumulate halo nodes in PLO
  !we also extend our plo_2_ago map to include a mapping from petsc local
  !order to application global order in the petsc partition halos
  !gives us ==> [n_halo]
  !gives us ==> [alo_2_plo_node(NNode+1:NTNode)]
  !gives us ==> [plo_2_alo_node(n_verts+1:n_verts+n_halo)]
  N_HALO = 0
  IF(NPROCS > 1)THEN
    DO I= 1, NNodeGL
      IF(NLID_X(I) /= 0 .AND. PARTID(I) /= MYID) THEN
        ICNT = ICNT + 1
        N_HALO = N_HALO + 1
        ALO_2_PLO_NODE(NLID_X(I)) = ICNT
        PLO_2_AGO_NODE(ICNT     ) = I
      ENDIF
    ENDDO
  ENDIF
  IF(PDEBUG) THEN
    WRITE(PUNIT,*) 'ALO_2_PLO_NODE WITH HALO:'
    DO I= 1, NTNode
      WRITE(PUNIT,*) I, ALO_2_PLO_NODE(I)
    ENDDO
  ENDIF

  !we now have a map from all all local application order nodes, including
  !halos to petsc local nodes, including  halos.  We also have a map from all local
  !petsc partition nodes, including the halo to application global order.
  !what we need is a map from petsc partition nodes including halo to petsc
  !global order, so that we can insert matrix values in local order and they will
  !get mapped accordingly.
  !we will construct a map of application global order to petsc global order and
  !use that to map our palo_2_ago into what we really need, a plo_2_pgo
  !gives us ==> [ago_2_pgo_node(1:NNodeGL)]
  ALLOCATE(AGO_2_PGO_NODE(NNodeGL)) ; AGO_2_PGO_NODE = 0
  ALLOCATE(CNT(NPROCS))         ; CNT   = 0
  ALLOCATE(ACCUM(NPROCS) )      ; ACCUM = 0
  DO I= 2, NPROCS
    ACCUM(I)  = ACCUM(I-1) + PART_NODES(I-1)
  ENDDO
  DO I= 1, NNodeGL
    PID = PARTID(I)
    CNT(PID) = CNT(PID) + 1
    AGO_2_PGO_NODE(I) = ACCUM(PID) + CNT(PID)
  ENDDO
  IF(PDEBUG)THEN
    WRITE(PUNIT,*)'AGO_2_PGO_NODE:'
    DO I= 1, NNodeGL
      WRITE(PUNIT,*) I, AGO_2_PGO_NODE(I)
    ENDDO
  ENDIF

  !now, use our ago_2_pgo_node to construct a plo_2_pgo_node
  !gives us ==> [plo_2_pgo_node(1:n_verts+n_halo)]
  ALLOCATE(PLO_2_PGO_NODE(N_VERTS+N_HALO)) ; PLO_2_PGO_NODE = 0
  DO I= 1, N_VERTS + N_HALO
    PLO_2_PGO_NODE(I) = AGO_2_PGO_NODE( PLO_2_AGO_NODE(I) )
  ENDDO
  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'COUNT PLO_2_PGO_NODE:'
    DO I=1, N_VERTS+N_HALO
      WRITE(PUNIT,*) I, PLO_2_PGO_NODE(I)
    ENDDO
  ENDIF

  !set a total problem size including petsc partition halo nodes
  ! ==> [Psize_whalo]
  PSIZE_EL_HALO = N_VERTS + N_HALO

  !---------------------------------------------------------
  ! set up the petsc matrix A
  !---------------------------------------------------------
  IF(PDEBUG) WRITE(PUNIT,*) 'SETTING UP THE MATRIX'

  !setup petsc local to global mapping (isl2g) to all for local indexing in assembly of matrix
  ALLOCATE(PTMP(N_VERTS+N_HALO)) ; PTMP = 0
  DO I= 1, N_VERTS + N_HALO
    PTMP(I) = PLO_2_PGO_NODE(I)-1
  ENDDO
  CALL ISLocalToGlobalMappingCreate(MPI_COMM_SELF,PSIZE_EL_HALO,PTMP,ISL2G_EL,IERR);CHKERRQ(IERR)

  !preallocate data
  CALL PETSc_AllocA

  !allow user to set options from command line
  CALL MatSetFromOptions(A_EL,IERR);CHKERRQ(IERR)

  !allowing local indexing for assembly of matrix
  CALL MatSetLocalToGlobalMapping(A_EL,ISL2G_EL,IERR);CHKERRQ(IERR)

  IF(PDEBUG) WRITE(PUNIT,*) 'MatSetLocalToGlobalMapping COMPLETE!'
  IF(PMSR)   WRITE(PSCRN,*) 'MATRIX PREALLOC    : COMPLETE!'

  !---------------------------------------------------------
  ! set up vector index sets used for scatter ops
  !---------------------------------------------------------

  !setup local vector index sets  (Petsc Ordering) -->  these are zero based
  CALL ISCreateStride(MPI_COMM_SELF,N_VERTS,0,1,ISLOCAL_EL,IERR);CHKERRQ(IERR)

  !setup global vector index set (Petsc Ordering) --> zero based
  CALL ISCreateStride(MPI_COMM_SELF,N_VERTS,ACCUM(MYID),1,ISGLOBAL_EL,IERR);CHKERRQ(IERR)
  IF(PMSR) WRITE(PSCRN,*)'VECTOR INDEX SETS  : COMPLETE!'

  !---------------------------------------------------------
  ! set up vector contexts: local vecs (xlocal,blocal)
  !---------------------------------------------------------
  CALL VecCreateSeq(MPI_COMM_SELF,N_VERTS,BL_EL,IERR);CHKERRQ(IERR)
  CALL VecCreateSeq(MPI_COMM_SELF,N_VERTS,XL_EL,IERR);CHKERRQ(IERR)

  !---------------------------------------------------------
  ! set up vector contexts: global parallel vecs (x,u,b)
  !---------------------------------------------------------

  CALL VecCreate(MPI_COMM_WORLD,X_EL,IERR);CHKERRQ(IERR)
  CALL VecSetSizes(X_EL,N_VERTS,PETSC_DECIDE,IERR);CHKERRQ(IERR)
  CALL VecSetFromOptions(X_EL,IERR);CHKERRQ(IERR)
  CALL VecDuplicate(X_EL,B_EL,IERR);CHKERRQ(IERR)
  IF(PMSR) WRITE(PSCRN,*) 'VECTOR SETUP       : COMPLETE!'

  !---------------------------------------------------------
  ! set up vector scatter operations
  !---------------------------------------------------------
  !for insertion (gather):  petsc local --> petsc global
  CALL VecScatterCreate(XL_EL,ISLOCAL_EL,X_EL,ISGLOBAL_EL,L2G_EL,IERR);CHKERRQ(IERR)

  !for a true scatter: petsc global --> petsc local
  CALL VecScatterCreate(X_EL,ISGLOBAL_EL,XL_EL,ISLOCAL_EL,G2L_EL,IERR);CHKERRQ(IERR)
  IF(PMSR) WRITE(PSCRN,*) 'VECTOR SCATTER OPS : COMPLETE!'

  !---------------------------------------------------------
  ! create the linear solver and setup the options
  !   set type
  !   set initial guess is zero or nonzero
  !   set convergence tolerance
  !   use set from options to allow override of these
  !   and defaults from the command line at runtime
  !---------------------------------------------------------

  CALL KSPCreate(MPI_COMM_WORLD,Ksp_EL,IERR);CHKERRQ(IERR)
  CALL KSPSetOperators(Ksp_EL,A_EL,A_EL,SAME_NONZERO_PATTERN,IERR);CHKERRQ(IERR)
  CALL KSPSetType(Ksp_EL,KSPGMRES,IERR);CHKERRQ(IERR)
  CALL KSPGetPC(Ksp_EL,Pc_EL,IERR);CHKERRQ(IERR)
  CALL PCSetType(Pc_EL,PCHYPRE,IERR);CHKERRQ(IERR)
  CALL PCHYPRESetType(Pc_EL,"boomeramg",IERR);CHKERRQ(IERR)
  IF(USE_LAST) CALL KSPSetInitialGuessNonzero(Ksp_EL,PETSC_TRUE,IERR);CHKERRQ(IERR)
  CALL KSPSetTolerances(Ksp_EL,RTOL,PETSC_DEFAULT_DOUBLE_PRECISION,  &
         PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,IERR);CHKERRQ(IERR)
  CALL KSPSetFromOptions(Ksp_EL,IERR);CHKERRQ(IERR)

  IF(PMSR) WRITE(PSCRN,*) 'SOLVER CONTEXT     : COMPLETE!'

  !---------------------------------------------------------
  ! deallocate unneeded arrays
  ! we need to keep the scatters (l2g,g2l)
  ! we need to keep the local to global matrix map (isl2g)
  ! we need to keep both local-local maps: alo_2_plo_node
  !                                        plo_2_alo_node
  ! we may or may need the index sets (islocal,isglobal)
  !---------------------------------------------------------
  DEALLOCATE(ACCUM)
  DEALLOCATE(PARTID)
  DEALLOCATE(PART_NODES)
  DEALLOCATE(PLO_2_AGO_NODE)
  DEALLOCATE(AGO_2_PGO_NODE)
  DEALLOCATE(CNT)
  DEALLOCATE(PLO_2_PGO_NODE)

  !report and reset logical first_entry
  IF(PMSR) WRITE(PSCRN,*) 'PETSc SETUP        : COMPLETE!'
  IF(PDEBUG) WRITE(PUNIT,*) 'FINISHING: ',TRIM(SUBNAME)

  RETURN
  END SUBROUTINE PETSc_SET

!=========================================================================================
  SUBROUTINE PETSc_AllocA
!============================================================================
!                         PETSc_AllocA
!     preallocate matrix A
!     assume our nonzero pattern will  not change
!     this will make the implementation much more efficient
!     malloc-ing on the fly can be very slow
!
!     currently assumes stencil is:
!       -node       + surrounding nodes
!       -node above + surrounding nodes
!       -node below + surrounding nodes
!
!     typical number of nonzero entries:  21
!
!     dependencies:  all internal to this module
!          except:  NNode,kbm1
!============================================================================
  USE ALL_VARS , ONLY : NNode, NTSN, NBSN, NNodeGL
  IMPLICIT NONE
  INTEGER :: NODE,PETSc_POS,NCOL2,NEYNUM,NEY,COLCOUNT,I
  INTEGER :: COL2(MAX_NZ_ROW)
  INTEGER,ALLOCATABLE  :: ONNZ2(:)
  INTEGER,ALLOCATABLE  :: DNNZ2(:)
  CHARACTER(LEN=20)    :: SUBNAME = 'PETSc_AllocA'
  INTEGER :: KBM1

  PetscInt :: IERR

  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ', TRIM(SUBNAME)
  !-----------------------------------------------------------------------------
  ! assemble matrix, use local vars (requires that MatSetLocalToGlobalMapping is setup)
  ! note: row and column are C order (starts from 0)
  ! see page 51 of petsc manual, bottom paragraph
  ! we will set entries row by row
  !-----------------------------------------------------------------------------
  ALLOCATE(ONNZ2(N_VERTS)) ; ONNZ2 = 0
  ALLOCATE(DNNZ2(N_VERTS)) ; DNNZ2 = 0

  DO NODE=1, NNode

    COL2  = 0
    NCOL2 = 0
    IF( ALO_2_PLO_NODE(NODE) <= N_VERTS) THEN
      NCOL2 = 1
      COL2(1) = ALO_2_PLO_NODE(NODE)

      DO NEYNUM =1, NTSN(NODE)-1
        NEY = NBSN(NODE,NEYNUM)
        IF(NEY == NODE) CYCLE
        NCOL2 = NCOL2 + 1
        !make sure we havent exceeded max_nz_row
        IF(NCOL2 > 10) THEN
          WRITE(PUNIT,*) 'WE HAVE EXCEEDED MAX_NZ_ROW==10'
          WRITE(PUNIT,*) 'INCREASE AND RECOMPILE'
          CALL MPI_FINALIZE(IERR)
          STOP
        ENDIF
        COL2(NCOL2) = ALO_2_PLO_NODE(NEY)
      ENDDO

      DO COLCOUNT =1, NCOL2
        IF(COL2(COLCOUNT) <= N_VERTS) THEN
          ONNZ2(ALO_2_PLO_NODE(NODE)) = ONNZ2(ALO_2_PLO_NODE(NODE)) + 1
        ELSE
          DNNZ2(ALO_2_PLO_NODE(NODE)) = DNNZ2(ALO_2_PLO_NODE(NODE)) + 1
        ENDIF
      ENDDO

    ENDIF

  ENDDO

  CALL MatCreateMPIAIJ(MPI_COMM_WORLD,N_VERTS,N_VERTS,NNodeGL,NNodeGL,0,DNNZ2,0,ONNZ2,A_EL,IERR);CHKERRQ(IERR)

  !deallocate non-zero counters
  DEALLOCATE(ONNZ2,DNNZ2)

  IF(PDEBUG) WRITE(PUNIT,*)'FINISHING: ',TRIM(SUBNAME)

  RETURN
  END SUBROUTINE PETSc_AllocA
!=============================================================================

  SUBROUTINE PETSc_CLEANUP
  !============================================================================
  !                         petsc_cleanup
  !      cleanup memory from petsc vars (called at end of FVCOM run)
  !
  ! external deps:  ksp,u,x,b,blocal,xlocal,A,ierr,pmsr from this module
  !============================================================================
  IMPLICIT NONE
  PetscInt :: IERR

  CALL KSPDestroy(Ksp_EL,IERR)   ;CHKERRQ(IERR)
  CALL VecDestroy(X_EL,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(B_EL,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(BL_EL,IERR)    ;CHKERRQ(IERR)
  CALL VecDestroy(XL_EL,IERR)    ;CHKERRQ(IERR)
  CALL MatDestroy(A_EL,IERR)     ;CHKERRQ(IERR)

  CALL petscfinalize(IERR)       ;CHKERRQ(IERR)

  IF(PMSR) WRITE(PSCRN,*) '    '
  IF(PMSR) WRITE(PSCRN,*) 'PETSc CLEANUP      : COMPLETE!'

  RETURN
  END SUBROUTINE PETSc_CLEANUP

  SUBROUTINE PETSc_SETICS_EL
  USE ALL_VARS
  IMPLICIT NONE
  INTEGER :: I, IK
  PetscReal :: QTERM
  PetscInt  :: IERR
  CHARACTER(LEN=20) :: SUBNAME = 'PETSc_SETICS'

  IF(.NOT.USE_LAST) THEN
    CALL VecSet(X_EL,ZERO,IERR);CHKERRQ(IERR)
    RETURN
  ENDIF

  DO I=1,N_VERTS
    IK = PLO_2_ALO_NODE(I)
    QTERM = EL(IK)
    CALL VecSetValues(XL_EL,1,I-1,QTERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

  CALL VecScatterBegin(XL_EL,X_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(XL_EL,X_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)

  RETURN
  END SUBROUTINE PETSc_SETICS_EL

  SUBROUTINE PETSc_SOLVER_EL
  IMPLICIT NONE
  PetscInt :: REASON
  PetscInt :: IERR

  CALL KSPSolve(Ksp_EL,B_EL,X_EL,IERR);CHKERRQ(IERR)

  CALL KSPGetConvergedReason(Ksp_EL,REASON,IERR)
  IF(REASON < 0) THEN
    IF(PMSR) THEN
      WRITE (PSCRN,*) 'PETSc SOLVER HAS DIVERGED IN EL CAL: STOPPING...'
    ENDIF
    CALL MPI_FINALIZE(IERR)
    STOP
  ELSE
    CALL KSPGetIterationNumber(Ksp_EL,ITS_EL,IERR);CHKERRQ(IERR)
    IF(PMSR) WRITE(PSCRN,*)'PETSc CONVERGED IN : ',ITS_EL,'AT RATE ',RTOL**(1.0d0/FLOAT(ITS_EL))
  ENDIF

  CALL VecNorm(X_EL,NORM_2,NORM_EL,IERR);CHKERRQ(IERR)
  IF(PMSR) WRITE(PSCRN,*) 'SOLUTION OF EL NORM      : ',NORM_EL

  RETURN
  END SUBROUTINE PETSc_SOLVER_EL

# endif
END MODULE MOD_PETSC

