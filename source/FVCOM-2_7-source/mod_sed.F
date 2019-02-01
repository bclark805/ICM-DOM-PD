!=======================================================================
! FVCOM Sediment Module 
!
! Copyright:    2005(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Authors:      G. Cowles 
!               School for Marine Science and Technology, Umass-Dartmouth
!
! Based on the Community Sediment Transport Model (CSTM) as implemented
!     in ROMS by J. Warner (USGS) 
!
! Comments:     Sediment Dynamics Module 
!
! Current FVCOM (main program) dependency
!   archive.F:   - calls archive_sed
!   init_sed.F:  - user defined sediment model initial conditions
!   mod_ncdio.F: - netcdf output includes concentration fields
!   us_fvcom.F:  - main calls sediment setup and sediment advance subs
!
! History
!   Feb 7, 2008: added initialization of bottom(:,:) to 0 (w/ T. Hamada)
!              : fixed loop bounds in hot start and archive for conc (w/ T. Hamada)
!              : added comments describing theoretical bases of dynamics
!   Feb 14,2008: added non-constant settling velocity for cohesive sediments (w/ T. Hamada) 
!              : updated vertical flux routine to handle non-constant vertical velocity (w/ T. Hamada)
!              : added a user-defined routine to calculate settling velocity based on concentration (w/ T. Hamada)
!              : added a user-defined routine to calculate erosion for a general case (w/ T. Hamada)
!
! ToDo with Hamada
!   1.) Add active layer thickness constraint on erosive flux
!   2.) Add potential for infinite sediment supply through inf_bed
!
!  Later
!   1.) Modify vertical flux routines to work with general vertical coordinate
!   2.) Add divergence term for bedload transport calc 
!   3.) Add ripple roughness calculation
!   4.) Add morphological change (bathymetry + vertical velocity condition) 
!   5.) Eliminate excess divisions and recalcs
!   
!=======================================================================
Module Mod_Sed  
# if defined (SEDIMENT)
Use Mod_Prec 
Use Mod_Types

# if defined (MULTIPROCESSOR)
  Use MOD_PAR
# endif

implicit none 

!--------------------------------------------------
!Sediment Type                         
!
! sname        => sediment name (silt,clay,etc)     
! stype        => sediment type: 'cohesive'/'non-cohesive'
! Sd50         => sediment mean diameter ()
! Wset         => mean sediment settling velocity
! tau_ce       => critical shear stress for erosion
! tau_cd       => critical shear stress for deposition 
! Srho         => sediment density
! Spor         => sediment porosity
! erate        => surface erosion mass flux         [kg m^-2 s^-1] 
! conc         => sed concentration in water column [kg m^-3]
! cnew         => sed concentration during update   [kg m^-3]
! mass         => sediment mass in bed layers       [kg m^-2]
! frac         => sediment fraction in bed layers   [-]
! bflx         => bedload sediment flux (kg/m^2)  (+out of bed) 
! eflx         => suspended sediment erosive flux (+out of bed)
! dflx         => suspended sed depositional flux (+into bed) 
! cdis         => concentration at river source 
! cflx         => store advective flux at open bndry 
! cobc         => user specd open bndry concentration 
! depm         => store deposited mass
! arraysize    => spatial size of sediment arrays
!--------------------------------------------------
type sed_type
  character(len=20) :: sname    
  character(len=20) :: stype 
  real(sp)          :: Sd50
  real(sp)          :: Wset  
  real(sp)          :: tau_ce
  real(sp)          :: tau_cd
  real(sp)          :: Srho 
  real(sp)          :: Spor 
  real(sp)          :: erate 
  real(sp)          :: cmax 
  real(sp)          :: cmin 
  real(sp)          :: crms 
  real(sp), pointer :: conc(:,:)
  real(sp), pointer :: cnew(:,:)
  real(sp), pointer :: mass(:,:)
  real(sp), pointer :: frac(:,:)
  real(sp), pointer :: bflx(:)        
  real(sp), pointer :: eflx(:)        
  real(sp), pointer :: dflx(:)        
  real(sp), pointer :: cdis(:)
  real(sp), pointer :: cflx(:,:)
  real(sp), pointer :: cobc(:)
  real(sp), pointer :: depm(:)     
  integer           :: arraysize
end type
 
public

!--------------------------------------------------
!Global Model Parameters
!
! sedfile   : sediment input parameter control file
! nsed      : number of sediment classes
! nbed      : number of layers in sediment bed
! min_Srho  : minimum Sediment density
! inf_bed   : true if bed has infinite sediment supply
! bedload   : true if bedload is to be considered
! susload   : true if suspended load is to be considered
! DTsed     : sediment model time step (seconds)
! T_model   : model time (seconds)
! taub      : array to hold bottom shear stress
! rho0      : mean density parameter used to convert tau/rho to tau
! tau_max   : max bottom stress (N-m^-2)
! tau_min   : min bottom stress (N-m^-2)
! thck_cr   : critical thickness for initiating new surface layer (m)
! n_report  : iteration interval for statistics printing
! sed_start : start interval for sed model
! sed_its   : sediment model iteration counter
! sed_nudge : flag for activiating sediment nudging on obc
! sed_alpha : sediment nudging relaxation factor 
! sed_ramp  : number of iterations over which to ramp sed nudging
! sed_source: flag for activiating sediment point sources
!
!--------------------------------------------------
character(len=120) :: sedfile
integer  :: nsed
integer  :: nbed
real(sp) :: min_Srho
logical  :: inf_bed
logical  :: sed_nudge
real(sp) :: sed_alpha
integer  :: sed_ramp  
logical  :: sed_source
logical  :: bedload
logical  :: susload
integer  :: n_report
integer  :: sed_start
integer  :: sed_its
real(sp) :: DTsed
real(sp) :: T_model
real(sp) :: tau_max 
real(sp) :: tau_min 
real(sp), allocatable :: taub(:)
real(sp), parameter   :: rho0 = 1025.
real(sp), parameter   :: thck_cr = .005 
logical, parameter    :: debug_sed = .false.

!--------------------------------------------------
!Bedload 
!
!Meyer-Peter Muller Bedload Formulation Parameters
!
! Shield_Cr_MPM => Effective Critical Shields number
! Gamma_MPM     => MPM Power Law Coefficient
! k_MPM         => MPM Multiplier
!
!--------------------------------------------------

real(sp) :: Shield_Cr_MPM    !  default: 0.047 
real(sp) :: Gamma_MPM        !  default: 1.5
real(sp) :: k_MPM            !  default: 8.0

!--------------------------------------------------
!Bed 
! data --> bed(horizontal,vertical,proptype)
! size --> bed(0:NTNode , nbed , n_bed_chars))
!
! with proptype=:
!   ithck  => layer thickness
!   iaged  => layer age 
!   iporo  => layer porosity  [volume of voids/total volume]
!
! global bed characteristics 
!
! n_bed_chars => number of bed characteristics 
!-------------------------------------------------

integer, parameter    :: n_bed_chars = 3
integer, parameter    :: ithck = 1
integer, parameter    :: iaged = 2
integer, parameter    :: iporo = 3
real(sp), allocatable :: bed(:,:,:)


!--------------------------------------------------
!Bottom (Exposed Sediment Layer)
! data --> bed(0:NTNode,n_bot_vars)
!
! with proptype=:
!   isd50  => mean grain diameter 
!   idens  => mean grain density     (not currently used)
!   iwset  => mean settling velocity (not currently used)
!   itauc  => mean critical erosion stress     
!   iactv  => active layer thickness            
!   nthck  => new total thickness of sediment layer
!   lthck  => last thickness of sediment layer 
!   dthck  => accumulated delta of layer thickness [m]
!   tmass  => total mass in sediment layer [kg/m^2]
!
! global bottom characteristics 
!   n_bot_vars => number of bottom characteristics
!
!-------------------------------------------------

integer, parameter    :: n_bot_vars  = 10
integer, parameter    :: isd50  = 1  
integer, parameter    :: idens  = 2
integer, parameter    :: iwset  = 3
integer, parameter    :: itauc  = 5
integer, parameter    :: iactv  = 6 
integer, parameter    :: nthck  = 7 
integer, parameter    :: lthck  = 8 
integer, parameter    :: dthck  = 9 
integer, parameter    :: tmass  = 10
real(sp), allocatable :: bottom(:,:)

!--------------------------------------------------
!Sediment Point Source Data 
!   sbc_tm:  time map for source data
!   seddis:  source data
!--------------------------------------------------
type(bc)              :: sbc_tm
real(sp), allocatable :: seddis(:,:,:)

!--------------------------------------------------
!Sediment Array                  
!--------------------------------------------------

type(sed_type), allocatable :: sed(:)

contains


!==========================================================================
! Allocate Data, Initialize Concentration Fields and Bed Parameters       
!   restart_yes: initialization type ('hot_start','cold_start')  
!==========================================================================
  Subroutine Setup_Sed(restart_yes) 
  use control, only : msr,casename,inpdir
  use lims,    only : NNode,numqbc_gl
  implicit none
  logical           :: restart_yes
  !-Local--------------------------
  logical           :: fexist
  integer           :: ised
  integer           :: i,k
  character(len=120):: fname

!----------------------------------------------------------------
! check arguments and ensure necessary files exist
!----------------------------------------------------------------

!Z Yang check point_1
	write(*,*)'mod_sed.F ZY_1'
	
  !initialize sediment model iteration counter
  sed_its = 0

  !ensure sediment parameter file exists
  sedfile = "./"//trim(inpdir)//"/"//trim(casename)//'_sediment.inp'
  inquire(file=trim(sedfile),exist=fexist)
  print*,'T.W., test1'

  if(.not.fexist)then
    write(*,*)'sediment parameter file: ',trim(sedfile),' does not exist'
    write(*,*)'stopping'
    call pstop
  end if
  print*,'T.W.,test2'

  !if restart case, ensure sediment restart file exists
  if(restart_yes)then
  fname = "./"//trim(inpdir)//"/"//trim(casename)//"_restart_sed.dat"
  inquire(file=trim(fname),exist=fexist)
  if(.not.fexist)then
    write(*,*)'sediment restart file: ',trim(fname),' does not exist'
    write(*,*)'halting.....'
    call pstop
  end if
  endif
  
  print*, 'T.W., test3' 
!----------------------------------------------------------------
! read in sediment parameter file                    
!----------------------------------------------------------------
  call read_sed_params

  print*, 'T.W., test4'
!----------------------------------------------------------------
! allocate data                                      
!----------------------------------------------------------------
  call alloc_sed_vars

  print*, 'T.W., test5'

!----------------------------------------------------------------
! initialize fields                                  
!----------------------------------------------------------------
  if(restart_yes)then
    call hot_start_sed
  else
    call init_sed 
  endif

!----------------------------------------------------------------
! setup open boundary condition nudging values       
!----------------------------------------------------------------
  if(sed_nudge) call setup_sed_obc 

!----------------------------------------------------------------
! setup point source river forcing data 
!----------------------------------------------------------------
  if(sed_source) call setup_sed_PTsource

  !--------------------------------------------------
  !Initialize Bed_Mass properties
  !--------------------------------------------------
   Do k=1,Nbed
     Do i=1,NNode
       Do ised=1,Nsed
         sed(ised)%mass(i,k) = bed(i,k,ithck)*       &
                              (1.0-bed(i,k,iporo))*  &
                               sed(ised)%Srho*       &
                               sed(ised)%frac(i,k)
       end do
     end do
  end do

  !--------------------------------------------------
  !Calculate Bottom Thicknesses
  !--------------------------------------------------

  do i=1,NNode
    do k=1,nbed
      bottom(i,nthck) = bottom(i,nthck) + bed(i,k,ithck)
    end do
  end do


!----------------------------------------------------------------
! report sediment parameters and statistics to screen
!----------------------------------------------------------------
  call report_sed_setup

!----------------------------------------------------------------
! convert parameters to MKS 
!----------------------------------------------------------------
  call convert_sed_params 

!----------------------------------------------------------------
! update bottom properties                           
!----------------------------------------------------------------
  call update_bottom_properties

End Subroutine Setup_Sed

!=======================================================================
!Convert Sed Params from nonstandard (mm, etc) to mks
!=======================================================================
  Subroutine Convert_Sed_Params 
  Implicit None
  integer :: i

  !convert settling velocity to m/s from input mm/s
  do i=1,nsed
    sed(i)%Wset = sed(i)%Wset*.001
  end do

  !convert mean grain diameter from mm to m 
  do i=1,nsed
    sed(i)%Sd50 = sed(i)%Sd50*.001
  end do

  End Subroutine Convert_Sed_Params

!=======================================================================
!Read Sediment Parameters From Sediment Input File
!=======================================================================
  Subroutine Read_Sed_Params
  Use Input_Util
  Implicit None
  integer linenum,i,k1
  real(sp)           :: ftemp
  character(len=120) :: stemp


  !read in number of sediment classes 
  Call Get_Val(nsed,sedfile,'NSED',line=linenum,echo=.false.)  

  !read in start interval for sed model 
  Call Get_Val(sed_start,sedfile,'SED_START',line=linenum,echo=.false.)  

  !read in interation interval for reporting 
  Call Get_Val(n_report,sedfile,'N_REPORT',line=linenum,echo=.false.)  

  !read in number of bed layers 
  Call Get_Val(nbed,sedfile,'NBED',line=linenum,echo=.false.)

  !read in logical for infinite sediment supply 
  Call Get_Val(inf_bed,sedfile,'INF_BED',line=linenum,echo=.false.)

  !read in logical for bedload calculation 
  Call Get_Val(bedload,sedfile,'BEDLOAD',line=linenum,echo=.false.)

  !read in logical for suspended load calculation 
  Call Get_Val(susload,sedfile,'SUSLOAD',line=linenum,echo=.false.)

  !read in minumum sediment density
  Call Get_Val(min_Srho,sedfile,'MIN_SRHO',line=linenum,echo=.false.)

  !read in nudging switch 
  Call Get_Val(sed_nudge,sedfile,'SED_NUDGE',line=linenum,echo=.false.)

  !read in nudging relaxation factor 
  if(sed_nudge) then
  Call Get_Val(sed_alpha,sedfile,'SED_ALPHA',line=linenum,echo=.false.)
  Call Get_Val(sed_ramp ,sedfile,'SED_RAMP ',line=linenum,echo=.false.)
  endif

  !read in point source switch 
  Call Get_Val(sed_source,sedfile,'SED_PTSOURCE',line=linenum,echo=.false.)

  !check values
  if(nsed < 1)then
    write(*,*)'sediment input file must have at least one sediment class'
    write(*,*)'currently nsed = ',nsed
    stop
  endif
  if(nbed < 1)then
    write(*,*)'sediment input file must have at least one bed layer'
    write(*,*)'currently nbed = ',nbed
    stop
  endif

  !allocate sediment data space
  Allocate(sed(nsed))

  !read in sediment parameters
  k1 = 1

  do i=1,nsed

    !read SED_NAME and mark position
    Call Get_Val(stemp,sedfile,'SED_NAME',line=linenum,echo=.false.,start=k1)
    sed(i)%sname = stemp
    k1        = linenum+1

    !read type
    Call Get_Val(stemp,sedfile,'SED_TYPE',line=linenum,echo=.false.,start=k1)
    sed(i)%stype = stemp

    !read mean diameter
    Call Get_Val(ftemp,sedfile,'SED_SD50',line=linenum,echo=.false.,start=k1)
    sed(i)%Sd50 = ftemp

    !read sediment density 
    Call Get_Val(ftemp,sedfile,'SED_SRHO',line=linenum,echo=.false.,start=k1)
    sed(i)%Srho = ftemp

    !read sediment settling rate 
    Call Get_Val(ftemp,sedfile,'SED_WSET',line=linenum,echo=.false.,start=k1)
    sed(i)%Wset = ftemp

    !read sediment surface erosion rate 
    Call Get_Val(ftemp,sedfile,'SED_ERAT',line=linenum,echo=.false.,start=k1)
    sed(i)%erate = ftemp

    !read sediment critical erosive shear stress 
    Call Get_Val(ftemp,sedfile,'SED_TAUE',line=linenum,echo=.false.,start=k1)
    sed(i)%tau_ce = ftemp

    !read sediment critical depositional shear stress 
    Call Get_Val(ftemp,sedfile,'SED_TAUD',line=linenum,echo=.false.,start=k1)
    sed(i)%tau_cd = ftemp

    !read sediment porosity
    Call Get_Val(ftemp,sedfile,'SED_PORS',line=linenum,echo=.false.,start=k1)
    sed(i)%Spor = ftemp

  end do

  ! read in bedload function parameters
  Call Get_Val(Shield_Cr_MPM,sedfile,'MPM_CS',line=linenum,echo=.false.)
  Call Get_Val(    Gamma_MPM,sedfile,'MPM_GM',line=linenum,echo=.false.)
  Call Get_Val(        k_MPM,sedfile,'MPM_K ',line=linenum,echo=.false.)
 
  End Subroutine Read_Sed_Params 
     
!=======================================================================
! Advance Sediment Model by time step DTSED 
!=======================================================================
  Subroutine Advance_Sed(DTin,Tin,Taub_in)
  Use Input_Util
  Use Scalar
  Use Control, only : ireport,iint,msr,par
  Use Lims,    only : NNode,NTNode,kbm1,numqbc,kb,nprocs,myid
  Use Mod_OBCS,only : iobcn
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif

  Implicit None
  real(sp), intent(in) :: DTin,Tin
  real(sp), intent(in) :: Taub_in(0:NTNode) 
  integer :: i,k,ised,l1,l2,ierr,d_cdis,d_cflx
  character(len=4) :: fnum
  real(sp) :: fact,ufact,temp

  !------------------------------------------------------
  !Check for Initialization
  !------------------------------------------------------
  if(iint < sed_start) return

  !------------------------------------------------------
  !Increment Sed Model Counter and Report Sed Model Init 
  !------------------------------------------------------
  sed_its = sed_its + 1
  if(msr .and.  sed_its == 1)then
    write(*,*)'========Sed Model Initiated======='
  endif

  !------------------------------------------------------
  !Set up Model Forcing  -> time step and bottom shear 
  !Convert Model tau/rho to tau [N/m^2]
  !------------------------------------------------------
  DTsed = DTin
  T_model = Tin
  taub = Taub_in*rho0
  tau_max = maxval(taub(1:NNode))
  tau_min = minval(taub(1:NNode))

  !------------------------------------------------------
  !Calculate Active Layer Thickness --> bottom(:,iactv)
  !------------------------------------------------------
  call calc_active_layer

  !------------------------------------------------------
  !Vertical Sediment Dynamics --> Sediment 'Model'
  !------------------------------------------------------
  if(bedload)then
    !calculate bedload

    !calculate bedload flux => sed(:)%bflx 

    !adjust surface bed properties
    call update_surface_bed('bedload')

    !adjust bottom properties
    call update_bottom_properties

    !update bed layer stratigraphy
    call update_stratigraphy
  endif

  if(susload)then
    !calculated depositional flux from settling
    call calc_deposition

    !calculate erosive flux
    call calc_erosion

    !adjust surface bed properties
    call update_surface_bed('susload')

    !add new surface layer if necessary
    call add_new_layer

    !update bed layer stratigraphy
    call update_stratigraphy

    !adjust bottom properties
    call update_bottom_properties
  endif
   
  !------------------------------------------------------
  !Horizontal Advection of Sediment Concentration 
  !------------------------------------------------------
  d_cdis = 1
  if(numqbc > 0) d_cdis = numqbc
  d_cflx = 1
  if(iobcn > 0) d_cflx = iobcn

  do i=1,nsed

    !set discharge sediment concentrations
    if(numqbc > 0 .and. sed_source)then
      call bracket(sbc_tm,T_model/3600.,l1,l2,fact,ufact,ierr)
      sed(i)%cdis(:) = ufact*seddis(:,i,l1) + fact*seddis(:,i,l2)
    endif

    call adv_scal(  sed(i)%conc, &
                    sed(i)%cnew, &
                    d_cdis,      &
                    sed(i)%cdis, &
                    d_cflx,      &
                    sed(i)%cflx, &
                    DTsed ,sed_source)
  end do

  !------------------------------------------------------
  !Vertical Diffusion of Sediment Concentration
  !------------------------------------------------------
  do i=1,nsed
    call vdif_scal(sed(i)%cnew,DTsed) 
  end do

  !------------------------------------------------------
  !Exchange Concentration on Processor Boundaries
  !------------------------------------------------------
# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,nsed
      call exchange(nc,NTNode,kb,myid,nprocs,sed(i)%cnew) 
    end do
  endif
# endif

  !------------------------------------------------------
  ! Open Boundary Conditions 
  !------------------------------------------------------
  if(sed_nudge)then
    if(sed_ramp == 0) sed_ramp = 1
    temp = min(float(sed_its/sed_ramp), 1.0)*sed_alpha
    do i=1,nsed
      call bcond_scal_OBC(sed(i)%conc, &
                          sed(i)%cnew, &
                          sed(i)%cflx, &
                          sed(i)%cobc, &
                          DTsed,       &
                          temp) 
    end do
  endif

  !------------------------------------------------------
  ! Point Source Boundary Conditions 
  !------------------------------------------------------
  if(sed_source)then !geoff???
    do i=1,nsed
      call bcond_scal_PTsource(sed(i)%conc, &
                               sed(i)%cnew, &
                               sed(i)%cdis) 
    end do
  endif

  !------------------------------------------------------
  ! Update Concentration Variables      
  !------------------------------------------------------
  do i=1,nsed
    sed(i)%conc = sed(i)%cnew
  end do

  !------------------------------------------------------
  ! Limit to Non-Negative (Until Venkat Limiter Debug) 
  !------------------------------------------------------
  do i=1,nsed
    sed(i)%conc = max(sed(i)%cnew,0.0)
  end do

  !------------------------------------------------------
  !Exchange Concentration on Processor Boundaries
  !------------------------------------------------------
# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,nsed
      call exchange(nc,NTNode,kb,myid,nprocs,sed(i)%cnew)
    end do
  endif
# endif

  !------------------------------------------------------
  ! Update Depth Delta
  !------------------------------------------------------
  call update_thickness_delta

  !------------------------------------------------------
  ! Report 
  !------------------------------------------------------
  if(mod(sed_its,n_report)==0)call sed_report    

!  call layer_report(1)

  End Subroutine Advance_Sed

!=======================================================================
! Allocate Sediment Variables                      
!=======================================================================
  Subroutine Alloc_Sed_Vars
  use lims,     only: MTElem,NTNode,kb,numqbc,kbm1
  use mod_obcs, only: iobcn
  implicit none
  integer i,tmp1,tmp2,tmp3

  !ensure arrays have nonzero dimension
  tmp1 = max(numqbc ,1)
  tmp2 = max(iobcn+1,1)
  tmp3 = max(iobcn  ,1)

  !allocate suspended sediment arrays
  do i=1,nsed
    sed(i)%arraysize = NTNode+1
    allocate(sed(i)%conc(0:NTNode,kb  ))      ; sed(i)%conc = 0.0
    allocate(sed(i)%cnew(0:NTNode,kb  ))      ; sed(i)%conc = 0.0
    allocate(sed(i)%mass(0:NTNode,nbed))      ; sed(i)%mass = 0.0
    allocate(sed(i)%frac(0:NTNode,nbed))      ; sed(i)%frac = 0.0
    allocate(sed(i)%bflx(0:NTNode     ))      ; sed(i)%bflx = 0.0
    allocate(sed(i)%eflx(0:NTNode     ))      ; sed(i)%eflx = 0.0
    allocate(sed(i)%dflx(0:NTNode     ))      ; sed(i)%dflx = 0.0
    allocate(sed(i)%cdis(tmp1     ))      ; sed(i)%cdis = 0.0
    allocate(sed(i)%cflx(tmp2,kbm1))      ; sed(i)%cflx = 0.0
    allocate(sed(i)%cobc(tmp3     ))      ; sed(i)%cobc = 0.0
    allocate(sed(i)%depm(0:NTNode     ))      ; sed(i)%depm = 0.0
  end do

  !allocate bottom shear stress array 
  allocate(taub(0:NTNode)) ; taub = 0.0

  !allocate bed data
  allocate(bed(0:NTNode , nbed , n_bed_chars)) ; bed = 0.0

  !allocate bottom data
  allocate(bottom(0:NTNode , n_bot_vars)) ; bottom = 0.0

  End Subroutine Alloc_Sed_Vars


!==========================================================================
! Report Sediment Setup To Screen 
!==========================================================================

  Subroutine Report_Sed_Setup
  use control, only: msr 
  implicit none
  integer :: i

  ! echo sediment model parameters to screen
  if(msr)then

  write(*,*)
  write(*,*)'------------------ Sediment Model Setup --------------------'
  write(*,*)'!  nbed                  :',nbed    
  write(*,*)'!  nsed                  :',nsed    
  write(*,*)'!  n_report              :',n_report
  if(sed_start > 0)then
    write(*,*)'!  sed_start             :',sed_start
  else
    write(*,*)'!  sed_start             : no delay'  
  endif
  write(*,*)'!  min_Srho              :',min_Srho
  if(sed_nudge)then
    write(*,*)'!  open boundary nudging :  active'
    write(*,*)'!  nudging relax factor  :',sed_alpha 
    write(*,*)'!  # nudging ramp its    :',sed_ramp  
  else
    write(*,*)'!  open boundary nudging :  not active'
  endif
  if(sed_source)then
    write(*,*)'!  point source forcing  :  active'
  else
    write(*,*)'!  point source forcing  :  not active'
  endif
  if(bedload)then
    write(*,*)'!  bedload dynamics      :  active'
  else
    write(*,*)'!  bedload dynamics      :  not active'
  endif
  if(susload)then
    write(*,*)'!  susload dynamics      :  active'
  else
    write(*,*)'!  susload dynamics      :  not active'
  endif
  if(inf_bed)then
    write(*,*)'!  sediment supply       :  infinite'
  else
    write(*,*)'!  sediment supply       :  finite'
  endif
  write(*,*)'!  Critical Shields      :',Shield_Cr_MPM
  write(*,*)'!  Bedload Exponent      :',Gamma_MPM     
  write(*,*)'!  Bedload Constant      :',k_MPM     
  write(*,*)'!'
  write(*,*)'!   class        type   Sd50    Wset   tau_ce   '//&
            'tau_cd  Erate   Spor    Srho'
  do i=1,nsed
    write(*,'(1X,A1,1X,A12,1X,A6,6(F7.4,1X),F8.2)')'!', &
              trim(sed(i)%sname),trim(sed(i)%stype(1:6)), &
              sed(i)%Sd50,sed(i)%Wset, &
              sed(i)%tau_ce,sed(i)%tau_cd,sed(i)%erate,sed(i)%Spor, &
              sed(i)%Srho  
  end do

  end if

  ! echo sediment model initial conditions to screen
  call sed_report  

  write(*,*)'------------------------------------------------------------'
  write(*,*)

  End Subroutine Report_Sed_Setup

!==========================================================================
! Calculate Sediment Statistics
!==========================================================================

  Subroutine Calc_Sed_Stats 
  use control, only: msr 
  use lims,    only: kbm1,NTNode
  implicit none
  integer :: i,dim1,dim2

  !set up limits
  dim1 = NTNode 
  dim2 = kbm1

  !calculate max/min/rms concentrations
  do i=1,nsed
    sed(i)%cmax = maxval(sed(i)%conc(1:dim1,1:dim2))
    sed(i)%cmin = minval(sed(i)%conc(1:dim1,1:dim2))
    sed(i)%crms = sum(dble(abs(sed(i)%conc(1:dim1,1:dim2))))
    sed(i)%crms = sed(i)%crms/(dim1*dim2)
  end do

  End Subroutine Calc_Sed_Stats

!==========================================================================
! Monitor Sediment Model: Write Statistics to Screen
!==========================================================================

  Subroutine Sed_Report  
  use control, only: msr 
  use lims   , only: NNode
  implicit none
  integer :: i

  ! calculate statistics
  call calc_sed_stats

  ! write statistics to screen
  if(.not.msr)return
  write(*,*  )'====================Sediment Model Stats======================'
  write(*,*  )'!  quantity              :     avg           max        '
  do i=1,nsed
    write(*,100)'!',trim(sed(i)%sname),'    :',sed(i)%crms,sed(i)%cmax
  end do
  write(*,*  )'!  max sed thick change  :     ',maxval(abs(bottom(1:NNode,dthck)))
  write(*,*  )'!  max bottom stress     :     ',tau_max             
  write(*,*  )'!  min bottom stress     :     ',tau_min             
  

  100 format(1x,a1,a20,a5,3f12.6)
  101 format(1x,a26,f12.6)

  End Subroutine Sed_Report  

!==========================================================================
! Calculate Bed Load using Dimensionless Meyer-Peter Muller Formulation   
!==========================================================================

  Real(sp) Function Bedload_T(R,rho0,D_50,tau_b) 
  use control, only: grav
  implicit none
  real(sp), intent(in) :: R
  real(sp), intent(in) :: rho0
  real(sp), intent(in) :: D_50 
  real(sp), intent(in) :: tau_b  
  !--------------------------------
  real(sp) :: tau_nd,bedld_nd

  !nondimensional bottom stress
  write(*,*)'make sure dimensions (tau/rho or tau) are correct'
  stop
  tau_nd = tau_b/(rho0*R*grav*D_50) 

  !compute nondimensional bedload
  bedld_nd = k_MPM*(MAX((tau_nd - Shield_Cr_MPM),0.0)**Gamma_MPM)

  !compute bedload (kg/(m s))
  bedload_t = sqrt(R*grav*D_50)*D_50*bedld_nd

  End Function Bedload_T 


!==========================================================================
! Calculate Erosion Rate (kg/m^2) 
! Calculate erosion rate using the method of Ariathurai and Arulanandan (1978)
!    -Erosion rates of cohesive soils. J. Hydraulics Division, ASCE, 104(2),279-282.-
!    surface_mass_flux[i] = erate[i]*(1-porosity)*bfrac[i]*(tau_w/tau_crit[i] - 1)
!      where
!      i               = index of sediment type
!      bfrac[i]        = fraction of bed composed by sediment i
!      erate[i]        = bed erodibility constant
!      porosity        = volume of voids/total volume in the top layer
!      tau_w           = shear stress on the bed
!      tau_crit[i]     = critical shear stress of 
!==========================================================================

  Subroutine Calc_Erosion
  use all_vars
  implicit none
  integer  :: ised ,i
!Z Yang 8/2/2010
! exclude erosion as a real variable as it is called as a function
! real(sp) :: bed_por,bed_frac,tau_w,tau_ce,erate,dep,erosion,active
  real(sp) :: bed_por,bed_frac,tau_w,tau_ce,erate,dep,active

  !Compute Erosive Flux (kg/m^2), limited by available material in active layer 
  do ised=1,nsed
    erate  = sed(ised)%erate
    tau_ce = sed(ised)%tau_ce 
    do i=1,NNode
      bed_por  = bed(i,1,iporo)
      bed_frac = sed(ised)%frac(i,1)
      tau_w    = taub(i)
      dep      = sed(ised)%dflx(i)
      active   = (1.0-bed_por)*Bed_Frac*(tau_w/tau_ce-1.0)*bottom(i,iactv) + dep
      sed(ised)%eflx(i) = DTsed*erosion(erate,bed_por,bed_frac,tau_w,tau_ce) !min(erosion,active) !geoff
    end do
  end do

  !Update Concentration in Bottom of Water Column
  do ised=1,nsed
    do i=1,NNode
      sed(ised)%conc(i,kbm1) = sed(ised)%conc(i,kbm1) + sed(ised)%eflx(i)/(d(i)*dz(i,kbm1))
    end do
  end do
  
  End Subroutine Calc_Erosion 


!==========================================================================
! Calculate Erosion (kgm^-2s^-1) using user-defined formula 
!==========================================================================
  Real(sp) Function Erosion(erate,bed_por,bed_frac,tau_w,tau_ce)
  implicit none
  real(sp), intent(in) ::  erate 
  real(sp), intent(in) ::  bed_por 
  real(sp), intent(in) ::  bed_frac 
  real(sp), intent(in) ::  tau_w 
  real(sp), intent(in) ::  tau_ce 

  !standard CSTM formulation 
  Erosion = MAX(Erate*(1.0-bed_por)*Bed_Frac*(tau_w/tau_ce-1.0),0.0)

  !Winterwerp
  !Erosion =  MAX(Erate*Bed_Frac*(tau_w/tau_ce-1.0),0.0)
  
  End Function Erosion
  
  
!==========================================================================
! Calculate Depositional Flux and Update Concentration from Settling
!==========================================================================
  Subroutine Calc_Deposition
  use all_vars
  implicit none
  integer  :: ised ,i,mcyc,ncyc,k
  real(sp) :: wset(kbm1) ,c(kbm1),dx(kbm1),flux,DTmax,eps,DTdep
  
  eps = epsilon(eps)

  !Loop over Sediment Class       
  do ised=1,nsed
    do i=1,NNode
      !initialize flux
      flux = 0.0

      !setup 1D concentration and grid arrays
      c(1:kbm1) = sed(ised)%conc(i,1:kbm1)
      dx(1:kbm1) = dz(i,1:kbm1)*d(i)

      !calculate the settling velocity if sediment is cohesive
      if(sed(ised)%stype == 'cohesive')then
        call calc_wset(kbm1,c,wset,sed(ised)%Wset)
      else
        wset = sed(ised)%Wset 
      endif

      !set up cycles (use max(CFL) == 1)
      DTmax = minval(dx(1:kbm1)/wset(1:kbm1))
      mcyc = int(DTsed/DTmax + 1. - eps)
      DTdep = DTsed/float(mcyc)

      !call flux-limited settling equation
      do ncyc = 1,mcyc
        call settle_flux(kbm1,c,dx,wset,flux,DTdep)
      end do

      !store solution in 3D array
      sed(ised)%conc(i,1:kbm1) = c(1:kbm1)

      !store depositional flux
      sed(ised)%dflx(i) = flux 
    end do
  end do

  End Subroutine Calc_Deposition

!==========================================================================
! Calculate settling velocity of cohesive sediment using the method of
!    Mehta, 1986 or
!    Hindered settlement strategy
!    or
!    User Defined
!==========================================================================
  Subroutine Calc_Wset(MElem,c,wset,Wset_Mean) 
  implicit none
  integer , intent(in   ) :: MElem
  real(sp), intent(in   ) ::  c(MElem)
  real(sp), intent(inout) :: wset(MElem)
  real(sp), intent(in   ) :: Wset_Mean
  integer :: i

  !setup for Ariake Bay
  real(sp) :: a = 0.5 !Mehta proportional constant
  real(sp) :: b = 2.0 !Mehta exponential constant
  real(sp) :: cmax = 1.0 !Concentration for hindered settling

  do i=1,MElem
    if(c(i) < cmax)then
       wset(i) = a*(c(i)**b)  
    else
       !hindered settling here
    endif
  end do

  End Subroutine Calc_Wset

!==========================================================================
! Calculate Settling Flux and Update Sediment Concentration in Water Column 
!   use the SLIP limiter
!   Jameson, A., Analysis and Design of Numerical Schemes for Gas Dynamics I
!      Artificial Diffusion, Upwind Biasing, Limiters, and their effects
!      on Accuracy and Multigrid Convergence, International Journal of 
!      Computational Fluid Dynamics, Vol 4, 171-218, 1995.
!
!      -second order accurate in smooth regions
!      -guaranteed local extremum diminishing
!      
!   conv:  convective flux
!   diss:  dissipative flux
!   
!==========================================================================
  Subroutine Settle_Flux(MElem,c,dx,wset,flux,deltat) 
  implicit none
  integer , intent(in   ) :: MElem
  real(sp), intent(inout) ::  c(MElem)
  real(sp), intent(in   ) :: dx(MElem)
  real(sp), intent(in   ) :: wset(MElem)
  real(sp), intent(inout) :: flux 
  real(sp), intent(in   ) :: deltat 
  real(sp) :: conv(MElem+1),diss(MElem+1)
  real(sp) :: cin(-1:MElem+2)
  real(sp) :: win(-1:MElem+2)
  real(sp) :: dis4,wvel
  integer  :: i

  !transfer to working array
  cin(1:MElem) = c(1:MElem)
  win(1:MElem) = wset(1:MElem) 

  !surface bcs (no flux)
  cin(0)  =  -cin(1) 
  cin(-1) =  -cin(1)
  win(0)  =  -wset(1)
  win(-1) =  -wset(1)
  
  !bottom bcs (extrapolate)
  cin(MElem+1) = cin(MElem) 
  cin(MElem+2) = cin(MElem)
  win(MElem+1) = win(MElem)
  win(MElem+2) = win(MElem)

  !flux computation
  do i=1,MElem+1
    wvel    = .5*(win(i)+win(i-1))  !settle velocity at interface
    dis4    = wvel/2.
    conv(i) = wvel*(cin(i)+cin(i-1))/2. 
    diss(i) = dis4*(cin(i)-cin(i-1)-lim(cin(i+1)-cin(i),cin(i-1)-cin(i-2))) 
  end do

  !zero out surface flux
  conv(1) = 0. ; diss(1) = 0.

  !update
  do i=1,MElem
    c(i) = cin(i) + (deltat/dx(i))*(-conv(i+1)+conv(i) + diss(i+1)-diss(i)) 
  end do

  !set bottom flux (> 0 = deposition)
  flux = flux + deltat*(conv(MElem+1)-diss(MElem+1))

  End Subroutine Settle_Flux
  
!==========================================================================
! Calculate LED Limiter L(u,v)  
!==========================================================================
  Function Lim(a,b)
  real(sp) lim,a,b
  real(sp) q,R
  real(sp) eps
  eps = epsilon(eps)
  
 ! exponent
 ! q = 0. !1st order
 ! q = 1. !minmod
   q = 2. !van leer

  R = abs(   (a-b)/(abs(a)+abs(b)+eps) )**q
  lim = .5*(1-R)*(a+b)

  End Function Lim

!==========================================================================
! Update Surface Bed Properties Due to Bedload/SusLoad Flux 
! Change:
!   sed(:)%mass(:,1)   to reflect increased/decrease mass (kg)
!   bed(:,1,itchk)     to reflect changed bed thickness properties
!==========================================================================

  Subroutine Update_Surface_Bed(fluxtype)
  use lims, only: NNode
  implicit none
  character(len=*) :: fluxtype
  real(sp) :: ddz,accum,flux,t_mass,eps
  integer  :: i,ised

  eps = epsilon(eps)

  !======= adjust for bedload flux ===================
  if(fluxtype == 'bedload')then

  !top layer mass and thickness
  do i=1,NNode
    do ised=1,nsed
      flux = sed(ised)%bflx(i)

      !change top layer bed mass of each sed type 
      accum  = sed(ised)%mass(i,1) - flux 
      sed(ised)%mass(i,1) = max(accum,0.0) 

      !change top layer bed thickness 
      ddz = flux/(sed(ised)%srho*(1.0-bed(i,1,iporo)))
      bed(i,1,ithck) = max(bed(i,1,ithck)-ddz,0.0)

    end do
  end do

  !top layer fractions 
  do i=1,NNode
    t_mass = 0
    do ised=1,nsed
      t_mass = t_mass + sed(ised)%mass(i,1)  
    end do
    do ised=1,nsed
      sed(ised)%frac(i,1) = sed(ised)%mass(i,1)/MAX(t_mass,eps)
    end do
  end do

  !======= adjust for susload flux ===================
  elseif(fluxtype =='susload')then
  
  do i=1,NNode
    do ised=1,nsed
      flux = sed(ised)%eflx(i)-sed(ised)%dflx(i)
  
      !if depositional + first deposit time step, save mass and loop  
      !http://woodshole.er.usgs.gov/project-pages/sediment-transport/
      sed(ised)%depm(i) = 0.0
      if(flux < 0)then                                  !!depositional
        if(T_model > bed(i,1,iaged)+1.1*DTsed .and. &   !!first deposit time step
           bed(i,1,ithck) > thck_cr)then                !!thickness surpasses critical
          sed(ised)%depm(i) = -flux
          cycle
        else                                            !!update age of surface layer
          bed(i,1,iaged) = T_model
        endif
      endif

      !change top layer bed mass of each sed type 
      accum  = sed(ised)%mass(i,1) - flux 
      sed(ised)%mass(i,1) = max(accum,0.0) 

      !change top layer bed thickness 
      ddz = flux/(sed(ised)%srho*(1.0-bed(i,1,iporo)))
      bed(i,1,ithck) = max(bed(i,1,ithck)-ddz,0.0)

    end do
  end do

  !top layer fractions 
  do i=1,NNode
    t_mass = 0
    do ised=1,nsed
      t_mass = t_mass + sed(ised)%mass(i,1)  
    end do
    do ised=1,nsed
      sed(ised)%frac(i,1) = sed(ised)%mass(i,1)/MAX(t_mass,eps)
    end do
  end do

  !======= error ======================
  else
    write(*,*)'argument to Update_Surface_Bed must be: "susload" or "bedload"'
    stop
  endif !fluxtype

  End Subroutine Update_Surface_Bed


!==========================================================================
! Update Bottom Properties due to changes in surface layer sed fractions 
!    Use geometric mean, good for log normal distribution
!==========================================================================
  Subroutine Update_Bottom_Properties
  use lims, only: NNode
  implicit none
  integer  :: i,ised
  real(sp) :: temp,sum1,sum2,sum3,sum4,eps

  ! initialize
  eps = epsilon(temp) 
  
  ! update bottom properties using new bed fractions
  do i=1,NNode
    sum1 = 1.0 ; sum2 = 1.0 ; sum3 = 1.0 ; sum4 = 1.0
    do ised=1,nsed
      sum1 = sum1*(sed(ised)%tau_ce)**sed(ised)%frac(i,1) 
      sum2 = sum2*(sed(ised)%Sd50  )**sed(ised)%frac(i,1)
      sum3 = sum3*(sed(ised)%wset  )**sed(ised)%frac(i,1)
      sum4 = sum4*(sed(ised)%Srho  )**sed(ised)%frac(i,1)
    end do
    bottom(i,itauc) = sum1
    bottom(i,isd50) = sum2
    bottom(i,iwset) = sum3
    bottom(i,idens) = max(sum4,min_Srho) 
  end do 
  
  End Subroutine Update_Bottom_Properties 

!==========================================================================
! Update Active Layer Thickness 
!   Use method of Harris and Wiberg (1997):
!   -Approaches to quantifying long-term continental shelf sediment transport
!   -with an example from the northern California STRESS mid-shelf site.
!   [Continental Shelf REsearch, 17, 1389-1418]
!
!   z_a = max[k_1*(tau_bot - tau_ce)*rho_0 , 0] + k_2 * D_50
!       where
!       k_1 = .007 (empirical constant)
!       k_2 = 6.0  (empirical constant)
!       tau_bot = shear stress on the bed
!       tau_ce  = critical shear stress for erosion of the bed (averaged over
!                 sediment classes)
!       D_50    = median grain diameter of the surface sediment
!       note: our shear stresses are nondimensionalized by rho, so
!             rho does not appear in actual calculation (below)
!==========================================================================

  Subroutine Calc_Active_Layer 
  use lims,     only: NNode
  implicit none
  integer :: i

  !calculate active layer thickness 
  do i=1,NNode
    bottom(i,iactv) = max(0.,.007*(taub(i)-bottom(i,itauc)))+6.0*bottom(i,isd50)
  end do

  End Subroutine Calc_Active_Layer
  
!==========================================================================
! Update Stratigraphy       
!   -Calculate empirical active layer thickness from bottom stress
!   -Expand top layer to this thickness by shifting necessary mass up
!    through the bed.
!==========================================================================

  Subroutine Update_Stratigraphy
  use lims,     only: NNode
  implicit none
  integer :: i,k,ksed,ks,ised
  real(sp) :: thck_avail,thck_to_add,eps,tmp1,tmp2,tmp3,top_layer_mass
  real(sp) :: delta(NNode) 
 
  !initialize
  eps = epsilon(tmp1) 
  delta = 0.0

  !calculate deficit between active layer thickness and surface layer thickness
  do i=1,NNode
    delta(i) = bottom(i,iactv)-bed(i,1,ithck)
  end do
  
  ! (single layer system) ensure top layer > active layer thickness 
  if(nbed == 1)then 
    do i=1,NNode
       if(delta(i) > 0) bottom(i,iactv) = bed(i,1,ithck)
    end do
    return
  endif
   
  ! (multi layer system) ensure top layer > active layer thickness 
  do i=1,NNode
    thck_to_add = delta(i)
    if(thck_to_add > 0.0)then !must redistribute layers 
      thck_avail=0.0

      !find fractional layer (below which there will be no mass after forming active layer) 
      Ksed=1 
      do k=2,nbed
        if (thck_avail < thck_to_add) then
          thck_avail=thck_avail+bed(i,k,ithck)
          Ksed=k
        end if
      end do

      !if not enough material in bed to satisfy active layer req, use all available 
      if (thck_avail < thck_to_add) then
        bottom(i,iactv)=bed(i,1,ithck)+thck_avail
        thck_to_add=thck_avail
      end if
  
      !update the bed mass of surface and fractional layers
      tmp2=MAX(thck_avail-thck_to_add,0.0)/MAX(bed(i,Ksed,ithck),eps)
      do ised=1,nsed
        tmp1=0.0 
        do k=1,Ksed
          tmp1=tmp1+sed(ised)%mass(i,k)
        end do
        sed(ised)%mass(i,1   ) = tmp1-sed(ised)%mass(i,Ksed)*tmp2
        sed(ised)%mass(i,Ksed) = sed(ised)%mass(i,Ksed)*tmp2
      end do

      !update thickness of fractional layer Ksed
      bed(i,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0)

      !update top layer bed fraction 
      top_layer_mass = 0.0
      do ised=1,nsed
        top_layer_mass = top_layer_mass + sed(ised)%mass(i,1)
      end do
      do ised=1,nsed
        sed(ised)%frac(i,1)=sed(ised)%mass(i,1)/MAX(top_layer_mass,eps) 
      end do

      !update bed thickness of top layer
      bed(i,1,ithck)=bottom(i,iactv)

      !pull layers Ksed to Bottom up to fill layer 2 down 
      do k=Ksed,Nbed
        ks=Ksed-2
        bed(i,k-ks,ithck)=bed(i,k,ithck)
        bed(i,k-ks,iporo)=bed(i,k,iporo)
        bed(i,k-ks,iaged)=bed(i,k,iaged)
        do ised=1,nsed
          sed(ised)%frac(i,k-ks) = sed(ised)%frac(i,k)
          sed(ised)%mass(i,k-ks) = sed(ised)%mass(i,k)
        end do
      end do
 
      !split what was in the bottom layer to fill empty bottom cells
      !note:  porosity of nbed (bed(i,nbed,iporo)) does not change.
      ks=Ksed-2
      tmp3=1.0/real(ks+1)
      do k=Nbed,Nbed-ks,-1
        bed(i,k,ithck)=bed(i,Nbed-ks,ithck)*tmp3
        bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
        do ised=1,nsed
          sed(ised)%frac(i,k)=sed(ised)%frac(i,Nbed-ks)
          sed(ised)%mass(i,k)=sed(ised)%mass(i,Nbed-ks)*tmp3
        end do
      end do

     end if !thck_to_add > 0
  end do !loop over nodes

    
  End Subroutine Update_Stratigraphy

!==========================================================================
! Add New Layer                 
!   If first deposit time step:  
!       a.) combine bottom layers 
!       b.) shift upper layers down
!       c.) create new top layer
!==========================================================================

  Subroutine Add_New_Layer 
  use lims,     only: NNode
  implicit none
  integer :: i,k,ksed,ised
  real(sp) :: eps,t_mass,cnt
 
  !initialize
  eps = epsilon(eps) 

  do i=1,NNode
  
    cnt = 0.0
    do ised=1,nsed
      cnt = cnt + sed(ised)%depm(i)
    end do

    if(cnt > 0)then !!need new surface layer

      if(nbed > 1)then!--->
      !combine bottom two layers
      bed(i,nbed,ithck) =     bed(i,nbed-1,ithck) + bed(i,nbed,ithck)
      bed(i,nbed,iporo) = .5*(bed(i,nbed-1,iporo) + bed(i,nbed,iporo))
      bed(i,nbed,iaged) = .5*(bed(i,nbed-1,iaged) + bed(i,nbed,iaged))

      t_mass =0.0
      do ised=1,nsed
        sed(ised)%mass(i,nbed) = sed(ised)%mass(i,nbed-1)+sed(ised)%mass(i,nbed) 
        t_mass = t_mass +sed(ised)%mass(i,nbed)
      end do
      do ised=1,nsed
        sed(ised)%frac(i,nbed) = sed(ised)%mass(i,nbed)/MAX(t_mass,eps)
      end do

      !push layers down
      do k=Nbed-1,2,-1
        bed(i,k,ithck) = bed(i,k-1,ithck)
        bed(i,k,iporo) = bed(i,k-1,iporo)
        bed(i,k,iaged) = bed(i,k-1,iaged)
        do ised =1,nsed
          sed(ised)%frac(i,k) = sed(ised)%frac(i,k-1)
          sed(ised)%mass(i,k) = sed(ised)%mass(i,k-1)
        end do
      end do

      !refresh top layer in multi layer bed
      bed(i,1,ithck)=0.0
      do ised=1,nsed
        sed(ised)%mass(i,1)=0.0
      end do
      end if ! <---

     

      !update surface layer properties
      t_mass = 0.0
      do ised=1,nsed
        sed(ised)%mass(i,1) = sed(ised)%mass(i,1)+ sed(ised)%depm(i)
        t_mass = t_mass + sed(ised)%mass(i,1)
        bed(i,1,ithck)=bed(i,1,ithck)+sed(ised)%depm(i)/    &
                         (sed(ised)%Srho*(1.0-bed(i,1,iporo)))
      end do
      bed(i,1,iaged) = T_model  
      do ised=1,nsed
        sed(ised)%frac(i,1) = sed(ised)%mass(i,1)/MAX(t_mass,eps) 
        sed(ised)%depm(i  ) = 0.0
      end do
 
    end if  !need new surface layer
  end do !node loop

  End Subroutine Add_New_Layer

!==========================================================================
! Update Total Sediment Thickness for Reporting
!==========================================================================

  Subroutine Update_Thickness_Delta
  use lims,     only: NNode
  implicit none
  integer :: i,k,ised
 
  !shift last thickness to lthck
  bottom(:,lthck) = bottom(:,nthck)
  bottom(:,nthck) = 0.0

  !calculate new thickness
  do i=1,NNode
    do k=1,nbed
      bottom(i,nthck) = bottom(i,nthck) + bed(i,k,ithck)
    end do
  end do

  !calculate new delta
  do i=1,NNode
    bottom(i,dthck) = bottom(i,dthck) + bottom(i,nthck) - bottom(i,lthck)
  end do

  !calculate new total mass
  do i=1,NNode
    bottom(i,tmass) = 0.0 
    do ised=1,nsed
      do k=1,nbed
        bottom(i,tmass) = bottom(i,tmass) + sed(ised)%mass(i,k)
      end do
    end do
  end do
  
  End Subroutine Update_Thickness_Delta

!==========================================================================
! Report Layer Stats for Debug   
!==========================================================================

  Subroutine Layer_Report(i)
  implicit none
  integer, intent(in) :: i
  integer :: j,ised

  write(*,*)'----------bottom props'
  write(*,*)'active thickness',bottom(i,iactv)
  write(*,*)'total  thickness',bottom(i,nthck)
  write(*,*)'mean grain size ',bottom(i,isd50)
  write(*,*)'mean density    ',bottom(i,idens)
  write(*,*)'mean settle V   ',bottom(i,iwset)
  write(*,*)'mean E stress   ',bottom(i,itauc)

  write(*,*)'--------- bed props'
  write(*,*)'total  thickness',bottom(i,nthck)
  write(*,*)'total  mass'     ,bottom(i,tmass)
  if(bottom(i,nthck) <= 0.)stop

  End Subroutine Layer_Report

!==========================================================================
! Dump sediment restart file             
!==========================================================================
  Subroutine Archive_Sed
  Use Mod_Prec 
  Use Lims, only: NNode,NNodeGL,kbm1,myid,nprocs
  Use Control, only: msr,serial,par
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif
  implicit none
  integer :: i,j,k,ised
  character(len=120) :: fname
  integer, parameter :: sedrestart = 81
  real(sp), allocatable :: gtemp(:,:),temp(:,:)
  real(sp), allocatable :: vgtemp(:),vtemp(:)

  !--------------------------------------------------
  !Open Sediment Restart File for Writing
  !--------------------------------------------------
  if(msr)then
    fname = "./restart_sed"
    open(unit=sedrestart,file=fname,form='unformatted')
  endif

  !--------------------------------------------------
  !Dump Sediment Model Info
  !--------------------------------------------------
  if(msr)then
    write(sedrestart) NNodeGL
    write(sedrestart) nsed
    write(sedrestart) nbed
  endif

  !--------------------------------------------------
  !Dump Bed Properties
  !--------------------------------------------------
  if(serial)then
    write(sedrestart) ((bed(i,j,iaged),i=1,NNodeGL),j=1,Nbed) !iaged
    write(sedrestart) ((bed(i,j,ithck),i=1,NNodeGL),j=1,Nbed) !ithck
    write(sedrestart) ((bed(i,j,iporo),i=1,NNodeGL),j=1,Nbed) !iporo
  endif
  
# if defined (MULTIPROCESSOR)
  if(par)then
  allocate(gtemp(NNodeGL,Nbed))
  allocate(temp(NNode,Nbed))
  
  !iaged
  temp(1:NNode,1:Nbed) = bed(1:NNode,1:Nbed,iaged)
  call gather(lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,myid,nprocs,nmap,temp,gtemp)
  if(msr) write(sedrestart) ((gtemp(i,j),i=1,NNodeGL),j=1,Nbed)

  !ithck
  temp(1:NNode,1:Nbed) = bed(1:NNode,1:Nbed,ithck)
  call gather(lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,myid,nprocs,nmap,temp,gtemp)
  if(msr) write(sedrestart) ((gtemp(i,j),i=1,NNodeGL),j=1,Nbed)

  !iporo 
  temp(1:NNode,1:Nbed) = bed(1:NNode,1:Nbed,iporo)
  call gather(lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,myid,nprocs,nmap,temp,gtemp)
  if(msr) write(sedrestart) ((gtemp(i,j),i=1,NNodeGL),j=1,Nbed)
  deallocate(gtemp,temp)
  endif

# endif

  !--------------------------------------------------
  !Dump Bed_Frac properties
  !--------------------------------------------------
  allocate(temp(NNode,Nbed))
  if(serial)then
    do ised=1,Nsed
      temp(1:NNode,1:Nbed) = sed(ised)%frac(1:NNode,1:Nbed)
      write(sedrestart) ((temp(i,j),i=1,NNode),j=1,Nbed) 
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    allocate(gtemp(NNodeGL,Nbed))
    do ised=1,Nsed
      temp(1:NNode,1:Nbed) = sed(ised)%frac(1:NNode,1:Nbed)
      call gather(lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,myid,nprocs,nmap,temp,gtemp)
      write(sedrestart) ((gtemp(i,j),i=1,NNodeGL),j=1,Nbed) 
    end do
    deallocate(gtemp)
  endif
# endif
  
  deallocate(temp)

  !--------------------------------------------------
  !Dump Sediment Concentrations
  !--------------------------------------------------
!  allocate(temp(NNode,Nbed))
  allocate(temp(NNode,kbm1))
  if(serial)then
    do ised=1,Nsed
!      temp(1:NNode,1:Nbed) = sed(ised)%cnew(1:NNode,1:kbm1)
      temp(1:NNode,1:kbm1) = sed(ised)%cnew(1:NNode,1:kbm1)
      write(sedrestart) ((temp(i,j),i=1,NNode),j=1,kbm1)
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
!    allocate(gtemp(NNodeGL,Nbed))
    allocate(gtemp(NNodeGL,kbm1))

    do ised=1,Nsed
!      temp(1:NNode,1:Nbed) = sed(ised)%cnew(1:NNode,1:kbm1)
      temp(1:NNode,1:kbm1) = sed(ised)%cnew(1:NNode,1:kbm1)
!      call gather(lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,myid,nprocs,nmap,temp,gtemp)
      call gather(lbound(temp,1),ubound(temp,1),NNode,NNodeGL,kbm1,myid,nprocs,nmap,temp,gtemp)
      write(sedrestart) ((gtemp(i,j),i=1,NNodeGL),j=1,kbm1)
    end do
    deallocate(gtemp)
  endif
# endif

  deallocate(temp)

  !--------------------------------------------------
  !Dump Bottom Thickness Change 
  !--------------------------------------------------
  if(serial)then
    write(sedrestart) (bottom(i,dthck),i=1,NNodeGL)
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    allocate(vgtemp(NNodeGL))
    allocate(vtemp(NNode))

    vtemp(1:NNode) = bottom(1:NNode,dthck)
    call gather(lbound(vtemp,1),ubound(vtemp,1),NNode,NNodeGL,Nbed,myid,nprocs,nmap,vtemp,vgtemp)
    if(msr) write(sedrestart) (vgtemp(i),i=1,NNodeGL)
    deallocate(vgtemp,vtemp)
  endif
# endif

  close(sedrestart)

  End Subroutine Archive_Sed
  

!==========================================================================
! Read sediment restart file             
!==========================================================================
  Subroutine Hot_Start_Sed
  Use Mod_Prec 
  Use Lims, only: NNode,kbm1,myid,nprocs,NNodeGL
  Use Control, only: msr,serial,par,casename,inpdir
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif
  implicit none
  integer :: i,j,k,ised,tmp
  character(len=120) :: fname
  integer, parameter :: sedrestart = 81
  real(sp), allocatable :: gtemp(:,:),temp(:,:)
  real(sp), allocatable :: vgtemp(:),vtemp(:)
  logical :: fexist

  !--------------------------------------------------
  !Open Sediment Restart File for Writing
  !--------------------------------------------------
  fname = "./"//trim(inpdir)//"/"//trim(casename)//"_restart_sed.dat"
  open(unit=sedrestart,file=trim(fname),form='unformatted')
  rewind(sedrestart)

  !----------------------------------------------------
  !Read Sediment Model Info, Check for Case Consistency
  !----------------------------------------------------

  read(sedrestart) tmp 
  if(msr)then
    if(tmp /= NNodeGL)then
      write(*,*)'number of nodes in sed restart file: ',tmp
      write(*,*)'number of nodes in current model:       ',NNodeGL  
      write(*,*)'inconsistent'
      write(*,*)'halting'
      call pstop
    endif
  endif

  read(sedrestart) tmp 
  if(msr)then
    if(tmp /= nsed)then
      write(*,*)'number of sed types in sed restart file: ',tmp
      write(*,*)'number of sed types current model:       ',nsed
      write(*,*)'inconsistent'
      write(*,*)'halting'
      call pstop
    endif
  endif

  read(sedrestart) tmp
  if(msr)then
    if(tmp /= nbed)then
      write(*,*)'number of bed layers in sed restart file: ',tmp
      write(*,*)'number of bed layers current model:       ',nbed  
      write(*,*)'inconsistent'
      write(*,*)'halting'
      call pstop
    endif
  endif

  !--------------------------------------------------
  !Read Bed Properties
  !--------------------------------------------------
  if(serial)then
    read(sedrestart) ((bed(i,j,iaged),i=1,NNodeGL),j=1,Nbed) !iaged
    read(sedrestart) ((bed(i,j,ithck),i=1,NNodeGL),j=1,Nbed) !ithck
    read(sedrestart) ((bed(i,j,iporo),i=1,NNodeGL),j=1,Nbed) !iporo
  endif
  
# if defined (MULTIPROCESSOR)
  if(par)then
    allocate(temp(NNode,Nbed))

    !iaged
    call pread(sedrestart,temp,lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,ngid(1),1,"bed")
    bed(1:NNode,1:Nbed,iaged) = temp(1:NNode,1:Nbed) 
  
    !ithck
    call pread(sedrestart,temp,lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,ngid(1),1,"bed")
    bed(1:NNode,1:Nbed,ithck) = temp(1:NNode,1:Nbed) 

    !iporo
    call pread(sedrestart,temp,lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,ngid(1),1,"bed")
    bed(1:NNode,1:Nbed,iporo) = temp(1:NNode,1:Nbed) 

    deallocate(temp)
  endif

# endif

  !--------------------------------------------------
  !Read Bed_Frac properties
  !--------------------------------------------------
  allocate(temp(NNode,Nbed))
  if(serial)then
    do ised=1,Nsed
      read(sedrestart) ((temp(i,j),i=1,NNode),j=1,Nbed) 
      sed(ised)%frac(1:NNode,1:Nbed) = temp(1:NNode,1:Nbed) 
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    do ised=1,Nsed
      call pread(sedrestart,temp,lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,ngid(1),1,"bed_frac")
      sed(ised)%frac(1:NNode,1:Nbed) = temp(1:NNode,1:Nbed) 
    end do
  endif
# endif
  
  deallocate(temp)

  !--------------------------------------------------
  !Read Sediment Concentrations
  !--------------------------------------------------
  allocate(temp(NNode,Nbed))
  if(serial)then
    do ised=1,Nsed
      read(sedrestart) ((temp(i,j),i=1,NNode),j=1,Nbed) 
      sed(ised)%conc(1:NNode,1:Nbed) = temp(1:NNode,1:kbm1) 
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    do ised=1,Nsed
      call pread(sedrestart,temp,lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,ngid(1),1,"sed_conc")
      sed(ised)%conc(1:NNode,1:Nbed) = temp(1:NNode,1:kbm1) 
    end do
  endif
# endif

  deallocate(temp)

  !--------------------------------------------------
  !Read Bottom Thickness  Delta
  !--------------------------------------------------
  if(serial)then
    read(sedrestart) (bottom(i,dthck),i=1,NNodeGL)
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    allocate(vtemp(NNode))
    vtemp(1:NNode) = bottom(1:NNode,dthck)
    call pread(sedrestart,vtemp,lbound(vtemp,1),ubound(vtemp,1),NNode,NNodeGL,1,ngid(1),1,"sed_conc")
    bottom(1:NNode,dthck) = vtemp(1:NNode) 
    deallocate(vtemp)
  endif
# endif

  End Subroutine Hot_Start_Sed
  
!==========================================================================
!  Setup open boundary nudging for sediment obc
!==========================================================================
  Subroutine Setup_Sed_OBC
  Use Control,  only : casename,inpdir,msr,serial,par
  Use Mod_OBCs, only : iobcn,i_obc_gl,iobcn_gl
# if defined (MULTIPROCESSOR)
  Use Mod_Par,  only : nlid
# endif

  integer :: i,j,tmp,ncnt,i1
  integer,  allocatable :: vtmp(:),atemp(:,:)
  real(sp), allocatable :: sed_obc_gl(:,:)
  integer, parameter    :: sedobc = 81
  character(len=120)    :: fname
  logical :: fexist
  

  !------------------------------------------------------
  ! make sure sediment nudging data file exists 
  !------------------------------------------------------
  fname = "./"//trim(inpdir)//"/"//trim(casename)//"_sed_nudge.dat"

  inquire(file=trim(fname),exist=fexist)
  if(msr .and. .not.fexist)then
    write(*,*)'sediment nudging data: ',fname,' does not exist'
    write(*,*)'halting.....'
    call pstop
  end if

  open(unit=sedobc,file=fname,form='formatted')
  !------------------------------------------------------
  ! allocate global data
  !------------------------------------------------------
  allocate(sed_obc_gl(iobcn_gl,nsed)) ; sed_obc_gl = 0.

  !header
  read(sedobc,*)
  read(sedobc,*)
  read(sedobc,*)

  !# obc nodes
  read(sedobc,*) tmp 

  if(tmp /= iobcn_gl .and. msr)then
    write(*,*)'number of open boundary nodes for sediment nudging: ',tmp
    write(*,*)'this is not consistent with number of open boundary nodes in'
    write(*,*)'the model: ',iobcn_gl
    write(*,*)'halting'
    call pstop
  endif

  allocate(vtmp(iobcn_gl))
  do i=1,iobcn_gl
    read(sedobc,*) i1,vtmp(i),(sed_obc_gl(i,j),j=1,nsed)
  end do
  close(sedobc)

  !------------------------------------------------------
  ! make sure global obc node list is consistent
  !------------------------------------------------------
  
  if(msr)then
  do i=1,iobcn_gl
    if(i_obc_gl(i) /= vtmp(i))then
      write(*,*)'obc node list for sediment nudging not consistent with model'
      write(*,*)'obc node list' 
      write(*,*)'halting'
      call pstop
    end if
  end do
  endif

  !---------------------------------------------------------
  ! shift sediment nudging data to local open boundary nodes
  !---------------------------------------------------------

  if(serial)then
    do i=1,iobcn_gl
      do j=1,nsed
        sed(j)%cobc(i) = sed_obc_gl(i,j)
      end do
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    allocate(atemp(iobcn_gl,nsed))
    ncnt = 0
    !!set up local open boundary nodes
    do i=1,iobcn_gl
      i1 = nlid( i_obc_gl(i) )
      if(i1 /= 0)then
        ncnt = ncnt + 1
        do j=1,nsed
          atemp(ncnt,j) = sed_obc_gl(i,j)
        end do
      end if
    end do

    !transfer sed nudging data global --> local
    if(ncnt > 0)then
      do j=1,nsed
        sed(j)%cobc(1:ncnt) = atemp(1:ncnt,j)
      end do
    end if

    deallocate(atemp) 
  end if
# endif

  if(msr)write(*,*)'obc nuding setup for sediment model: complete'

  !---------------------------------------------------------
  ! cleanup 
  !---------------------------------------------------------
  deallocate(sed_obc_gl)

  End Subroutine Setup_Sed_OBC


!==========================================================================
! Setup point source (river) loading of sediment 
! Only Called if Number of Rivers in Model Domain > 0
!==========================================================================
  Subroutine Setup_Sed_PTsource
  Use Control, only : casename,inpdir,msr,serial,par
  Use Lims,    only : numqbc_gl,numqbc
  Use BCS,     only : riv_gl2loc
  Use Mod_Clock
# if defined (MULTIPROCESSER)
  Use Mod_Par
# endif

  integer :: i,j,ns,ierr
  integer :: ntime,nriv
  integer :: sed_riv
  logical :: fexist
  integer, parameter :: sedriv = 81
  character(len=120) :: fname
  character(len=13)  :: tstring
  real(sp), allocatable :: temp(:,:,:)
  real(sp) :: ttime

  !-------------------------------------------------------
  ! Read Sediment River Data                             
  !-------------------------------------------------------

  fname = "./"//trim(inpdir)//"/"//trim(casename)//"_sed_ptsource.dat"

  inquire(file=trim(fname),exist=fexist)
  if(msr .and. .not.fexist)then
    write(*,*)'sediment point source data: ',fname,' does not exist'
    write(*,*)'halting.....'
    call pstop
  end if

  open(unit=sedriv,file=fname,form='formatted')
  if(msr)read(sedriv,*) ntime,nriv

  !------------------------------------------------------
  ! make sure number of rivers is consistent 
  !------------------------------------------------------
  if(msr)then
    if(nriv /= numqbc_gl)then
      write(*,*)'number of rivers in sediment point source data file: ',nriv   
      write(*,*)'number of rivers in model: ',numqbc_gl 
      write(*,*)'inconsistent'
      write(*,*)'halting'
      call pstop
    end if
  endif
    

# if defined (MULTIPROCESSOR)
  if(par)call mpi_bcast(ntime,1,mpi_integer,0,mpi_comm_world,ierr)
# endif

  !------------------------------------------------------
  ! read data times and data and broadcast
  !------------------------------------------------------
  sbc_tm%ntimes = ntime  
  sbc_tm%label  = 'Freshwater Discharge of Sediment'
  allocate(sbc_tm%times(ntime))
  allocate(temp(numqbc_gl,nsed,ntime)) 

  if(msr)then
    do i = 1, ntime
      read(sedriv,*) ttime
      sbc_tm%times(i) = ttime
      do ns=1,nsed
        read(sedriv,*) (temp(j,ns,i),j=1,numqbc_gl)  
      end do
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)call mpi_bcast(sbc_tm%times,ntime,mpi_f,0,mpi_comm_world,ierr)
  if(par)call mpi_bcast(temp,nsed*numqbc_gl*ntime,mpi_f,0,mpi_comm_world,ierr)
# endif

  !--------------------------------------------------------
  ! allocate local data array and transform global--> local 
  !--------------------------------------------------------

 if (numqbc > 0) then 

  allocate(seddis(numqbc,nsed,ntime))

  if(serial) seddis = temp

# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,ntime
      do j=1,nsed
         seddis(1:numqbc,j,i) = temp(riv_gl2loc(1:numqbc),j,i) 
      end do
    end do
  end if 
# endif

 end if

  !--------------------------------------------------------
  ! report results 
  !--------------------------------------------------------
  if(msr)then
    write(*,*)'!'
    call gettime(tstring,3600*int(sbc_tm%times(1)))
    write(*,102)'!  sed ptsource data beg :  ',tstring           
    call gettime(tstring,3600*int(sbc_tm%times(sbc_tm%ntimes)))
    write(*,102)'!  sed ptsource data end :  ',tstring           
    write(*,101)'!  max conc               :',maxval(temp)
  end if
  
  !---------------------------------------------------------
  ! cleanup 
  !---------------------------------------------------------
  deallocate(temp)
  close(sedriv)

  !---------------------------------------------------------
  ! format statements  
  !---------------------------------------------------------
   101  format(1x,a26,f10.4)  
   102  format(1x,a28,a13)  

 End Subroutine Setup_Sed_PTsource

!comment by T.W., copied from ZY's mod_sed.F in the sediment model

!==========================================================================
! Dump sediment results to medm binary file format             
!==========================================================================
  Subroutine Out_Binary_Sed(IINTT)
!  USE Mod_main
  Use Mod_Prec 
  Use Lims, only: NNode,NNodeGL,kbm1,myid,nprocs
  Use Control
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif
  implicit none
  integer :: i,j,k,ised
!  character(len=120) :: fname
!  integer, parameter :: sedrestart = 81
  real(sp), allocatable :: gtemp(:,:),temp(:,:)
  real(sp), allocatable :: vgtemp(:),vtemp(:)
  INTEGER, INTENT(IN) :: IINTT
  CHARACTER(LEN=4) :: FILENUMBER
  CHARACTER(LEN=120) :: DIR
  CHARACTER(LEN=120) :: FNAME

  !--------------------------------------------------
  !Dump Sediment Concentrations
  !--------------------------------------------------

!open medm file for sediment
   IF(MSR)THEN
     WRITE(FILENUMBER,'(I4.4)') IINTT
     FNAME = TRIM(CASENAME)//'_sed'//FILENUMBER//'.dat'
     DIR = TRIM(OUTDIR)//"/medm"
     OPEN(817,FILE=TRIM(DIR)//"/"//TRIM(FNAME),STATUS='unknown',FORM='unformatted')
     REWIND(817)
     WRITE(IPT,*)'DUMPING SEDIMENT MEDM FILE: ',TRIM(FNAME)
   END IF
!end open, start writing

 !--------------------------------------------------
  !Dump Sediment Model Info
  !--------------------------------------------------

!output heading info.
  if(msr)then
    WRITE(817) IINTT,THOUR, NNodeGL, nsed, nbed
  endif
!----------------------------------------------------------!
!--output water column sediment concentration--
!----------------------------------------------------------!
  allocate(temp(NNode,kbm1))
  if(serial)then
    do ised=1,Nsed
      temp(1:NNode,1:kbm1) = sed(ised)%cnew(1:NNode,1:kbm1)
      write(817) ((temp(i,j),i=1,NNode),j=1,kbm1)
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    allocate(gtemp(NNodeGL,kbm1))
    do ised=1,Nsed
      temp(1:NNode,1:kbm1) = sed(ised)%cnew(1:NNode,1:kbm1)
      call gather(lbound(temp,1),ubound(temp,1),NNode,NNodeGL,kbm1,myid,nprocs,nmap,temp,gtemp)
     if(msr)  then
       write(817) ((gtemp(i,j),i=1,NNodeGL),j=1,kbm1)
     end if
    end do
    deallocate(gtemp)
  endif
# endif
  deallocate(temp)
!-------------------------------------------------------!

  !--------------------------------------------------
  !--Dump Bed Properties
  !--------------------------------------------------
  if(serial)then
    write(817) ((bed(i,j,ithck),i=1,NNodeGL),j=1,Nbed) !ithck = 1
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
  allocate(gtemp(NNodeGL,Nbed))
  allocate(temp(NNode,Nbed))

  !ithck
  temp(1:NNode,1:Nbed) = bed(1:NNode,1:Nbed,ithck)
  call gather(lbound(temp,1),ubound(temp,1),NNode,NNodeGL,Nbed,myid,nprocs,nmap,temp,gtemp)
  if(msr) write(817) ((gtemp(i,j),i=1,NNodeGL),j=1,Nbed)

  deallocate(gtemp,temp)
  endif
# endif

IF(MSR) close(817)

  End Subroutine Out_Binary_Sed

# endif
End Module Mod_Sed 
