
!!==================================================================================
!!==================================================================================
!!
!! US FVCOM VISIT SIMULATION INTERFACE 
!! PROGRAMER: David StUebe 
!!
!!==================================================================================
!!==================================================================================


!  MODULE INTERFACE DOES NOT WORK WITH VISIT LIBRARY NAMING SCHEME
!
!==================================================================================
!==================================================================================

# ifdef VISIT

MODULE MOD_VISIT
  USE CONTROL, only : MSR, IPT
  USE MOD_PREC
#  if defined (NETCDF_IO)
  USE MOD_LAG, only : LAG_ON, particle_list
# endif

  IMPLICIT NONE
  ! Variables Used in us_fvcom.F and run_data.F go here

   CHARACTER(LEN=80) :: VISIT_OPT   !!(advanced/basic)  
   LOGICAL :: VISIT_CMD_DUMP

   real*4 :: visit_time_ext
   real*4 :: visit_time_int
   integer :: visit_cycle

! Variables only used in visitsim.F and mod_visit.F

  LOGICAL :: VisitOneTStep

  LOGICAL :: VisitHalt

  LOGICAL :: VisitAnimate
  
  integer :: VisitRunFlag 
  
  integer :: VisitStep
  INTEGER :: VisitStepCount
    
  integer :: VisitParRank
  
!  integer :: loopcount
  
  INTEGER, parameter :: VISIT_COMMAND_PROCESS = 0
  INTEGER, parameter :: VISIT_COMMAND_SUCCESS = 1
  INTEGER, parameter :: VISIT_COMMAND_FAILURE = 2

  INTEGER, parameter :: VISIT_STEP_EXT = 0
  INTEGER, parameter :: VISIT_STEP_INT = 1
  INTEGER, parameter :: VISIT_STEP_10XINT = 2
  INTEGER, parameter :: VISIT_STEP_100XINT = 3

  
  INTEGER, parameter :: VISIT_TWODMESH = 1
  INTEGER, parameter :: VISIT_BATHYMESH = 2
  INTEGER, parameter :: VISIT_SSHMESH = 3
  INTEGER, parameter :: VISIT_LAYERMESH = 4
  INTEGER, parameter :: VISIT_LEVELMESH = 5
  
  
  ! Debug stuff!
  INTEGER :: VDB

# ifdef VISIT_DEBUG
  integer :: BROADCASTINTCOUNT
  
  integer :: BROADCASTSTRCOUNT
  
  integer :: SLAVECALLBACKCOUNT
# endif
  
  TYPE VisitMeshType
     LOGICAL :: MESHALLOCATED
     REAL(SP) :: Updated_Time ! the time at which the mesh was last updated  
     INTEGER :: DATAOWNER
     INTEGER :: NDIMS
     INTEGER ::Nodes
     INTEGER :: Zones
     INTEGER :: LCONN 
     INTEGER :: MESH
     INTEGER :: LYRS
     INTEGER :: LVLS
     INTEGER, POINTER, DIMENSION(:,:) :: NV
     INTEGER, POINTER, DIMENSION(:,:) :: GHOSTZONES
     REAL*4, POINTER, DIMENSION(:,:) :: VX,VY,VZ
  END TYPE VisitMeshType
  
  
  TYPE VisitLag
     REAL(SP) :: Updated_Time ! the time at which the mesh was last updated  
     INTEGER :: NDIMS
     INTEGER ::Nodes
     REAL*4, POINTER, DIMENSION(:) :: VX,VY,VZ,VU,VW,VV,VS, VD
  END TYPE VisitLag

  TYPE VisitSphereVel
     REAL*4, POINTER, DIMENSION(:,:) :: VelX,VelY,VelZ
     REAL(SP) :: Updated_Time ! the time at which the velocity was last updated       
  END TYPE VisitSphereVel



  TYPE(VisitMeshType), TARGET, DIMENSION(5) :: VISIT_MESH

  TYPE(VisitLag), TARGET :: VISIT_LAGDATA
  
# if defined (SPHERICAL)
  TYPE(VisitSphereVel), target :: VisitSphericalVel,&
       & VisitSphericalAVel, VisitSphericalWindVel
# if defined (ICE)
  TYPE(VisitSphereVel), target :: VisitSphericalIceVel
# endif

# endif

  SAVE
  
  
  !===================================================================================|
CONTAINS   !!INCLUDED SUBROUTINES FOLLOW
  !===================================================================================|
  !===================================================================================
  


  !============================================================================
  !============================================================================
  !============================================================================
  !
  ! Adding fvcom to visit data handeling functions:
  ! FVCOM data is in 2 dimensional arrays with rows and columns padding
  ! the idicies. These functions make passing the data to visit as
  ! transparent as possible to the user by making the mesh data complex.
  !
  !============================================================================
  !============================================================================
  !============================================================================
  
  SUBROUTINE PRINT_MESH(vmp)
    implicit none
    TYPE(VISITMESHTYPE), POINTER :: vmp
    
    if (.NOT. ASSOCIATED(VMP)) then
       WRITE(VDB,*)"Mesh data pointer passed to UPDATE MESH was not associa&
            &ted"
       return
    end if
    

    WRITE(VDB,*)"==========================================="
    WRITE(VDB,*)"VMP%MESH=",VMP%MESH  
    if(VMP%MESHALLOCATED) then
       WRITE(VDB,*)"VMP has benn allocated" 
    else
       WRITE(VDB,*)"VMP has NOT benn allocated" 
    end if
    WRITE(VDB,*)"VMP%NDIMS=",VMP%NDIMS
    WRITE(VDB,*)"VMP%LYRS=",VMP%LYRS
    WRITE(VDB,*)"VMP%LVLS=",VMP%LVLS
    WRITE(VDB,*)"VMP%Nodes=",VMP%Nodes
    WRITE(VDB,*)"VMP%Zones=",VMP%Zones
    WRITE(VDB,*)"VMP%LCONN=",VMP%LCONN
    WRITE(VDB,*)"VMP%DataOwner=",VMP%DataOwner
    WRITE(VDB,*)"VMP%Updated_Time=",VMP%Updated_Time
    WRITE(VDB,*)"==========================================="
    
    
  END SUBROUTINE PRINT_MESH
  
  SUBROUTINE UPDATE_MESH(vmp,error)
    USE ALL_VARS
    implicit none
    include "visitfortransiminterface.inc"
    integer, intent(out) :: error
    integer :: ind, VnCP, lind, uind
    TYPE(VISITMESHTYPE), POINTER :: vmp
    REAL*4, pointer, Dimension(:) :: PX,PY,PZ

    error =-1 ! BAD RESULT!
    
    if (.NOT. ASSOCIATED(VMP)) then
       if(MSR) WRITE(VDB,*) "Mesh data pointer passed to UPDATE MESH was not associa&
            &ted"
       return
    end if    


# ifdef VISIT_DEBUG
    WRITE(VDB,*)"Starting update_mesh"
    CALL PRINT_MESH(VMP)
# endif

    IF (.NOT.VMP%MESHALLOCATED) then
       
       SELECT CASE(VMP%MESH)
          
       CASE ( VISIT_TWODMESH )
          VnCP=4  ! The number of integers to define one cell
          VMP%NDIMS=2 ! Spatial dimension of the mesh
          VMP%LYRS=1 ! The number of layers
          VMP%LVLS=1 ! The number of levels
          VMP%Nodes=VMP%LVLS * NTNode ! Total number of nodes
          VMP%Zones=VMP%LYRS * MTElem ! Total number of zones or cells
          VMP%LCONN=VnCP*VMP%Zones ! The number of integer in the
          ! connectivity array passed to visit
          VMP%DATAOWNER = VISIT_OWNER_SIM
         ! SINCE THE TWOD MESH IS STATIC SET THE TIME TO ZERO AND LEAVE IT
          VMP%Updated_Time=-1
       CASE ( VISIT_BATHYMESH )
          VnCP=4  ! The number of integers to define one cell
          VMP%NDIMS=3
          VMP%LYRS=1
          VMP%LVLS=1
          VMP%Nodes= VMP%LVLS * NTNode
          VMP%Zones= VMP%LYRS * MTElem
          VMP%LCONN=VnCP*VMP%Zones
          VMP%DATAOWNER = VISIT_OWNER_SIM
           ! SINCE THE BATHYMETRY MESH IS STATIC SET THE TIME TO ZERO AND LEAVE IT
          VMP%Updated_Time=-1
       CASE ( VISIT_SSHMESH )
          VnCP=4  ! The number of integers to define one cell
          VMP%NDIMS=3
          VMP%LYRS=1
          VMP%LVLS=1
          VMP%Nodes= VMP%LVLS * NTNode
          VMP%Zones= VMP%LYRS * MTElem
          VMP%LCONN=VnCP*VMP%Zones
          VMP%Nodes=NTNode
          VMP%Zones=MTElem
          VMP%LCONN=VnCP*MTElem
          VMP%DATAOWNER = VISIT_OWNER_SIM
          ! SINCE THE SSH MESH IS NOT STATIC SET THE TIME TO -1
          VMP%Updated_Time=-1
       CASE ( VISIT_LAYERMESH )

          VnCP=7  ! The number of integers to define one cell
          VMP%NDIMS=3
          VMP%LYRS=KBM2
          VMP%LVLS=KBM1
          VMP%Nodes= VMP%LVLS * NTNode
          VMP%Zones= VMP%LYRS * MTElem
          VMP%LCONN=VnCP*VMP%Zones
          VMP%Nodes= KBM1 * NTNode
          VMP%Zones= KBM2 * MTElem
          VMP%LCONN=VnCP*MTElem*KBM2
          VMP%DATAOWNER = VISIT_OWNER_SIM
          ! SINCE THE LAYER MESH IS NOT STATIC SET THE TIME TO -1
          VMP%Updated_Time=-1
       CASE ( VISIT_LEVELMESH )
          VnCP=7  ! The number of integers to define one cell
          VMP%NDIMS=3 
          VMP%LYRS=KBM1
          VMP%LVLS=KB
          VMP%Nodes= VMP%LVLS * NTNode
          VMP%Zones= VMP%LYRS * MTElem
          VMP%LCONN=VnCP*VMP%Zones
          VMP%DATAOWNER = VISIT_OWNER_SIM
          ! SINCE THE LEVEL MESH IS NOT STATIC SET THE TIME TO -1
          VMP%Updated_Time=-1
          
       CASE DEFAULT
          if(MSR) WRITE(VDB,*)"CASE(DEFAULT) in UPDATE_MESH: bad mesh! Calling PSTOP:"
          call PSTOP
       END SELECT

       VMP%MESHALLOCATED=.TRUE.
       
       
       ! NV must be allocated by the total number of zones(MTElem * layers)
       ! by VnCP
       
       ! All other variables will be zones by layer or nodes by layers
       ALLOCATE(VMP%NV(VnCP,VMP%Zones),STAT=error)
       if ( error  .ne. 0 ) then
           if(MSR) WRITE(VDB,*)"Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
               &TE NV, error= ",error,"; MYID= ",MYID
          return
       end if

       if (VnCP == 4) then
          VMP%NV(1,:)=VISIT_CELL_TRI
          VMP%NV(2:4,:)=Transpose(NV(1:MTElem,1:3)) -1
          
       else if (VnCP == 7) then
          VMP%NV(1,:)=VISIT_CELL_WEDGE
          DO ind =0,(VMP%LYRS-1)

             lind=ind*MTElem+1
             uind=(ind+1)*MTElem
             !           1    MTElem
             !          MTElem+1 2*MTElem
             !
             ! (LYRS-1)*MTElem+1 LYRS*MTElem
             VMP%NV(2:4,lind:uind)=ind*NTNode + TRANSPOSE(NV(1:MTElem,1:3)) -1
             
             VMP%NV(5:7,lind:uind)=(ind+1)*NTNode + TRANSPOSE(NV(1:MTElem,1:3)) -1
          END DO
          
       else 
          if(MSR)WRITE(VDB,*) "Unknown Cell type: VnCP=", VnCP
          if(MSR)WRITE(VDB,*) "known types: 4 => Triangle"
          if(MSR)WRITE(VDB,*) "             7 => Wedge"
          if(MSR)WRITE(VDB,*)"Calling PSTOP"
          Call PSTOP
       end if
       
       ALLOCATE(VMP%GHOSTZONES(MTElem,VMP%LYRS),STAT=error)
       if ( error  .ne. 0 ) then
          if(MSR)WRITE(VDB,*) "Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
               &TE NV, error= ",error,"; MYID= ",MYID
          return
       end if
       
       DO ind = 1,VMP%LYRS
          VMP%GHOSTZONES(1:MElem,ind)=0
          if(MTElem .GT. MElem) VMP%GHOSTZONES((MElem+1):MTElem,ind)=1
       END DO
       
       
       ALLOCATE(VMP%VX(NTNode,VMP%LVLS),STAT=error)
       if (error  .ne. 0 ) then
          if(MSR)WRITE(VDB,*) "Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
             &TE VX, MYID=",MYID
          return
       end if
       
       
       
       ALLOCATE(VMP%VY(NTNode,VMP%LVLS),STAT=error)
       if (error  .ne. 0 ) then
          if(MSR)WRITE(VDB,*) "Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
               &TE VY, MYID=",MYID
          return
       end if
       
       
       
       ALLOCATE(VMP%VZ(NTNode,VMP%LVLS),STAT=error)
       if (error  .ne. 0 ) then
          if(MSR)WRITE(VDB,*) "Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
               &TE VZ, MYID=",MYID
          return
       end if
       
       ! Add static data
       DO ind=1,VMP%LVLS
          VMP%VX(1:NTNode,ind)=VX(1:NTNode)
          VMP%VY(1:NTNode,ind)=VY(1:NTNode)
       END DO
       
       if (VMP%MESH == VISIT_TWODMESH) VMP%VZ(1:NTNode,1)=0
       
       if (VMP%MESH == VISIT_BATHYMESH) VMP%VZ(1:NTNode,1)=-H(1:NTNode)

# if defined (SPHERICAL)
       if (VMP%MESH == VISIT_TWODMESH .OR. &
            &VMP%MESH == VISIT_BATHYMESH) then

          PX=>VMP%VX(:,1)
          PY=>VMP%VY(:,1)
          PZ=>VMP%VZ(:,1)
          Call Sphere2Cart(PX,PY,PZ,NTNode)
          

       end if
# endif
      
    END IF  ! 
    
    
    SELECT CASE(VMP%MESH)
       
    CASE ( VISIT_TWODMESH )
       ! Do NOTHING
    CASE ( VISIT_BATHYMESH )
     ! Do NOTHING
    CASE ( VISIT_SSHMESH )
       
       if  (VMP%Updated_Time .NE. VISIT_TIME_EXT) then
          
# if defined (SPHERICAL)
          VMP%VX(1:NTNode,1)=VX(1:NTNode)
          VMP%VY(1:NTNode,1)=VY(1:NTNode)
          VMP%VZ(1:NTNode,1)=EL(1:NTNode)
          
          PX=>VMP%VX(:,1)
          PY=>VMP%VY(:,1)
          PZ=>VMP%VZ(:,1)
          Call Sphere2Cart(PX,PY,PZ,NTNode)
# else
          VMP%VZ(1:NTNode,1)=EL(1:NTNode)
          VMP%Updated_Time=VISIT_TIME_EXT
# endif
!          if (MSR) WRITE(VDB,*) "TIME CHANGED: Updated SSH"
       end if
       
    CASE ( VISIT_LAYERMESH )
       
       if  (VMP%Updated_Time .NE. VISIT_TIME_INT) then
# if defined (SPHERICAL)
          Do ind = 1,VMP%LVLS
             VMP%VX(1:NTNode,ind)=VX(1:NTNode)
             VMP%VY(1:NTNode,ind)=VY(1:NTNode)
             VMP%VZ(1:NTNode,ind)=EL(1:NTNode)+ ZZ(1:NTNode,ind)*D(1:NTNode)
             
             PX=>VMP%VX(:,ind)
             PY=>VMP%VY(:,ind)
             PZ=>VMP%VZ(:,ind)
             Call Sphere2Cart(PX,PY,PZ,NTNode)
          End Do
# else
          Do ind = 1,VMP%LVLS
             VMP%VZ(1:NTNode,ind)=EL(1:NTNode)+ ZZ(1:NTNode,ind)*D(1:NTNode)
          End Do
# endif       
!          if (MSR) WRITE(VDB,*) "TIME CHANGED: Updated SigmaLayers"
          VMP%Updated_Time=VISIT_TIME_INT
       end if
       
    CASE ( VISIT_LEVELMESH )

       if  (VMP%Updated_Time .NE. VISIT_TIME_INT) then
# if defined (SPHERICAL)
          Do ind = 1,VMP%LVLS
             VMP%VX(1:NTNode,ind)=VX(1:NTNode)
             VMP%VY(1:NTNode,ind)=VY(1:NTNode)
             VMP%VZ(1:NTNode,ind)=EL(1:NTNode)+ Z(1:NTNode,ind)*D(1:NTNode)

             PX=>VMP%VX(:,ind)
             PY=>VMP%VY(:,ind)
             PZ=>VMP%VZ(:,ind)
             Call Sphere2Cart(PX,PY,PZ,NTNode)
          End Do

# else
          Do ind = 1,VMP%LVLS
             VMP%VZ(1:NTNode,ind)=EL(1:NTNode)+ Z(1:NTNode,ind)*D(1:NTNode)
          End Do
# endif      
!          if (MSR) WRITE(VDB,*) "TIME CHANGED: Updated SigmaLevels"
          VMP%Updated_Time=VISIT_TIME_INT
       end if
       
    CASE DEFAULT
       if(MSR) WRITE(VDB,*)"CASE(DEFAULT) in UPDATE_MESH setting values: bad mesh! Calling PSTOP:"
       call PSTOP
    END SELECT
    


!    if (MSR)   WRITE(VDB,*) "FINISHED MESH UPDATE"
!    if(MSR)    CALL  PRINT_MESH(VMP)

    
    error = 0
  END SUBROUTINE UPDATE_MESH

  
# if defined (SPHERICAL)
  SUBROUTINE SPHERE2CART(PX,PY,PZ,n)
    IMPLICIT NONE  
    integer, intent(in) :: n
    REAL*4, pointer, Dimension(:), intent(inout) :: PX,PY,PZ
    
    real*4, parameter :: rad=100000.0;
    real*4, parameter :: d2r=3.14159/180.0;
    
    real*4, dimension(n) :: tx, ty, tz
    
    tx = (rad+PZ)*cos(PX*d2r)*cos(PY*d2r);
    ty = (rad+PZ)*sin(PX*d2r)*cos(PY*d2r);
    tz = (rad+PZ)*sin(PY*d2r);
    
    PX=tx
    PY=ty
    PZ=tz
    
  END SUBROUTINE SPHERE2CART
  
  
  SUBROUTINE UpdateSphereVel(VSV)
    USE ALL_VARS, only : MTElem, KBM1, U, V, WW
    implicit none
    TYPE(VisitSphereVel), POINTER :: VSV
    REAL*4, pointer, Dimension(:) :: PU,PV,PW
    integer :: ind

    if (.NOT.associated(VSV%Velx)) allocate(VSV%Velx(MTElem,KBM1))
    if (.NOT.associated(VSV%Vely)) allocate(VSV%Vely(MTElem,KBM1))
    if (.NOT.associated(VSV%Velz)) allocate(VSV%Velz(MTElem,KBM1))
    
    VSV%Velx=U(1:MTElem,1:KBM1)
    VSV%Vely=V(1:MTElem,1:KBM1)
    VSV%Velz=WW(1:MTElem,1:KBM1)

    Do ind=1,KBM1
       PU=>VSV%VelX(:,ind)
       PV=>VSV%VelY(:,ind)
       PW=>VSV%VelZ(:,ind)
       CALL SphereVel2Cart(PU,PV,PW)
    End Do

    Nullify(PU)
    Nullify(PV)
    Nullify(PW)

    VSV%Updated_Time=Visit_Time_INT
    
  END SUBROUTINE UpdateSphereVel

  SUBROUTINE UpdateSphereAVel(VSV)
    USE ALL_VARS, only : MTElem, UA, VA
    implicit none
    TYPE(VisitSphereVel), POINTER :: VSV
    REAL*4, pointer, Dimension(:) :: PU,PV,PW
    
    if (.NOT.associated(VSV%Velx)) allocate(VSV%Velx(MTElem,1))
    if (.NOT.associated(VSV%Vely)) allocate(VSV%Vely(MTElem,1))
    if (.NOT.associated(VSV%Velz)) allocate(VSV%Velz(MTElem,1))
    
    VSV%Velx(1:MTElem,1)=UA(1:MTElem)
    VSV%Vely(1:MTElem,1)=VA(1:MTElem)
    VSV%Velz=0
    
    PU=>VSV%VelX(:,1)
    PV=>VSV%VelY(:,1)
    PW=>VSV%VelZ(:,1)
    CALL SphereVel2Cart(PU,PV,PW)
    
    Nullify(PU)
    Nullify(PV)
    Nullify(PW)
    
    VSV%Updated_Time=Visit_Time_EXT
    
  END SUBROUTINE UpdateSphereAVel

  SUBROUTINE UpdateSphereWindVel(VSV)
    USE ALL_VARS, only : MTElem, UUWIND, VVWIND
    implicit none
    TYPE(VisitSphereVel), POINTER :: VSV
    REAL*4, pointer, Dimension(:) :: PU,PV,PW
    
    if (.NOT.associated(VSV%Velx)) allocate(VSV%Velx(MTElem,1))
    if (.NOT.associated(VSV%Vely)) allocate(VSV%Vely(MTElem,1))
    if (.NOT.associated(VSV%Velz)) allocate(VSV%Velz(MTElem,1))
    
    VSV%Velx(1:MTElem,1)=UUWIND(1:MTElem)
    VSV%Vely(1:MTElem,1)=VVWIND(1:MTElem)
    VSV%Velz=0
    
    PU=>VSV%VelX(:,1)
    PV=>VSV%VelY(:,1)
    PW=>VSV%VelZ(:,1)
    CALL SphereVel2Cart(PU,PV,PW)
    
    Nullify(PU)
    Nullify(PV)
    Nullify(PW)
    
    VSV%Updated_Time=Visit_Time_INT
    
  END SUBROUTINE UpdateSphereWindVel
  
#if defined (ICE)
  SUBROUTINE UpdateSphereIceVel(VSV)
    USE mod_ice2d, only: UICE2, VICE2
    USE ALL_VARS, only : MTElem
    implicit none
    TYPE(VisitSphereVel), POINTER :: VSV
    REAL*4, pointer, Dimension(:) :: PU,PV,PW
    
    if (.NOT.associated(VSV%Velx)) allocate(VSV%Velx(MTElem,1))
    if (.NOT.associated(VSV%Vely)) allocate(VSV%Vely(MTElem,1))
    if (.NOT.associated(VSV%Velz)) allocate(VSV%Velz(MTElem,1))
    
    VSV%Velx(1:MTElem,1)=UICE2(1:MTElem)
    VSV%Vely(1:MTElem,1)=VICE2(1:MTElem)
    VSV%Velz=0
    
    PU=>VSV%VelX(:,1)
    PV=>VSV%VelY(:,1)
    PW=>VSV%VelZ(:,1)
    CALL SphereVel2Cart(PU,PV,PW)
    
    Nullify(PU)
    Nullify(PV)
    Nullify(PW)
    
    VSV%Updated_Time=Visit_Time_INT
    
  END SUBROUTINE UpdateSphereIceVel
# endif
  
  SUBROUTINE SPHEREVEL2CART(PU,PV,PW)
    USE ALL_VARS, only : XC,YC, MTElem
    IMPLICIT NONE  
    REAL*4, pointer, Dimension(:), intent(inout) :: PU,PV,PW
    real*4, parameter :: d2r=3.14159/180.0;
    real*4, dimension(MTElem) :: tu, tv, tw


    tu=-sin(XC(1:MTElem)*d2r)*PU &
         & -sin(YC(1:MTElem)*d2r)*cos(XC(1:MTElem)*d2r)*PV &
         & + cos(XC(1:MTElem)*d2r)*cos(YC(1:MTElem)*d2r)*PW

    tv=cos(XC(1:MTElem)*d2r)*PU &
         & -sin(YC(1:MTElem)*d2r)*sin(XC(1:MTElem)*d2r)*PV &
         & + sin(XC(1:MTElem)*d2r)*cos(YC(1:MTElem)*d2r)*PW

    tw= 0.0 &
         & + cos(YC(1:MTElem)*d2r)*PV &
         & + sin(YC(1:MTElem)*d2r)*PW


    PU = tu
    PV = tv
    PW = tw

  END SUBROUTINE SPHEREVEL2CART

  SUBROUTINE DeAllocate_SphereVel(VSV,error)
    implicit none
    integer, intent(out) :: error
    TYPE(VisitSphereVel), POINTER :: VSV
    error =-1 ! BAD RESULT!
    
    if (associated(VSV%VELX)) deallocate(VSV%VELX)
    if (associated(VSV%VELY)) deallocate(VSV%VELY)
    if (associated(VSV%VELZ)) deallocate(VSV%VELZ)
    

    error = 0
    
  END SUBROUTINE DeAllocate_SphereVel


# endif

  SUBROUTINE DeAllocate_MESH(vmp,error)
    implicit none
    integer, intent(out) :: error
    TYPE(VISITMESHTYPE), POINTER :: vmp
    error =-1 ! BAD RESULT!
    
    if (associated(VMP%VX)) deallocate(VMP%VX)
    if (associated(VMP%VY)) deallocate(VMP%VY)
    if (associated(VMP%VZ)) deallocate(VMP%VZ)
    
    if (associated(VMP%NV)) deallocate(VMP%NV)
    
    if (associated(VMP%GHOSTZONES)) deallocate(VMP%GHOSTZONES)
  
    VMP%MESHALLOCATED=.FALSE.
    error = 0
    
  END SUBROUTINE DeAllocate_MESH
  

#  if defined (NETCDF_IO)
! UPDATE FOR LAGRANGIAN DATA
  SUBROUTINE UPDATE_LAG(LMP,error)
    USE LINKED_LIST
    implicit none
    integer :: ind, error
    ! FOR LAG TRACKING MESH DATA
    type(link_node), pointer :: lp
    type(VISITLAG), pointer  :: LMP
    
    error = -1

    if  (LMP%Updated_Time .NE. visit_time_int .OR. .NOT. associated(LMP%VX)) then
       
       if (associated(LMP%VX)) deallocate(LMP%VX)
       if (associated(LMP%VY)) deallocate(LMP%VY)
       if (associated(LMP%VZ)) deallocate(LMP%VZ)

       if (associated(LMP%VU)) deallocate(LMP%VU)
       if (associated(LMP%VV)) deallocate(LMP%VV)
       if (associated(LMP%VW)) deallocate(LMP%VW)

       if (associated(LMP%VS)) deallocate(LMP%VS)
       if (associated(LMP%VD)) deallocate(LMP%VD)

       
       
       LMP%NODES = listsize(particle_list)
       allocate(LMP%VX(LMP%NODES))
       allocate(LMP%VY(LMP%NODES))
       allocate(LMP%VZ(LMP%NODES))

       allocate(LMP%VU(LMP%NODES))
       allocate(LMP%VV(LMP%NODES))
       allocate(LMP%VW(LMP%NODES))
       
       allocate(LMP%VS(LMP%NODES))
       allocate(LMP%VD(LMP%NODES))
       
       ind = 1
       lp  => particle_list%first%next
       do
          if(.not. associated(lp) ) exit  !end of list, exit
          
          LMP%VX(ind)=lp%v%x(1)
          LMP%VY(ind)=lp%v%x(2)
          LMP%VZ(ind)=lp%v%zloc

          LMP%VU(ind)=lp%v%U
          LMP%VV(ind)=lp%v%V
          LMP%VW(ind)=lp%v%W

          LMP%VS(ind)=lp%v%S

          LMP%VD(ind)=lp%v%pathlength
          
          ind=ind+1
          lp => lp%next                          !set object
       end do
       
       
       if (MSR) WRITE(VDB,*) "TIME CHANGED: Updated LAGRANGIAN MESH"
       LMP%Updated_Time=visit_time_int
       error = 0

    else
       error = 0
    end if
    
  END SUBROUTINE UPDATE_LAG

  SUBROUTINE DeAllocate_LAG(LMP,error)
    USE LINKED_LIST
    implicit none
    integer, intent(out) :: error
    ! FOR LAG TRACKING MESH DATA
    type(VISITLAG), pointer  :: LMP
    
    error = -1

    if (associated(LMP%VX)) deallocate(LMP%VX)
    if (associated(LMP%VY)) deallocate(LMP%VY)
    if (associated(LMP%VZ)) deallocate(LMP%VZ)
    
    if (associated(LMP%VU)) deallocate(LMP%VU)
    if (associated(LMP%VV)) deallocate(LMP%VV)
    if (associated(LMP%VW)) deallocate(LMP%VW)
    
    if (associated(LMP%VS)) deallocate(LMP%VS)
    if (associated(LMP%VD)) deallocate(LMP%VD)
    
    error = 0
        
  END SUBROUTINE DeAllocate_LAG

# endif

END MODULE MOD_VISIT


#else
! if visit is not defined compile a dummy subroutine!

subroutine visit_dummy_mod

implicit none

end subroutine visit_dummy_mod

#endif

!end MODULE MOD_VISIT
