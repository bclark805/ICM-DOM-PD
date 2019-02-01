module mod_lag
#  if defined (NETCDF_IO)
!------------------------------------------------------------------------------
! Lagrangian particle tracking using linked lists
!
!  Current Capabilities (September 08, 2006)
!     1.) Multiprocessor Formulation:
!           Subdomains will only traject particles contained in their interior
!           Particles passing across interprocessor boundaries are exchanged 
!           between neighboring processors.
!     2.) Restart Capability:
!           The code can be restarted and the particle output file will
!           act as if the code had run with intermission.  Read
!           below for details.
!     3.) General tracking:  particles are trajected at each internal time step
!         (dti) using a 4 stage RK scheme and the (u,v,sigma) velocity field.
!
!     4.) Output File Contents
!         particle position, (x,y,z, + sigma) as well as an associated scalar
!         variable at the particle position.  Current options include salinity
!         temperature, densities, and eddy viscosities but can be modified 
!         easily to provide user defined fields such as a biological quantitiy
!         or tracer concentration. 
!
!     5.) initial particle position file.  This file is in netcdf format.
!         It must contain information about number of particles, location
!         (x,y,z) as well as global parameters including an information
!         string, time of last output, number of outputs made.  The last
!         two parameters are used when restarting the model.  To create
!         an initial position file, start from the create_initial_lagpos.f90
!         utility supplied in the utilities folder of your FVCOM distribution.
!         You can use this script to easily convert an ascii particle file
!         to the FVCOM netcdf format or you can generate a grid of particles
!         within a user defined box.
!
!  Future Capabilities
!     1.) Horizontal Random Walk.  
!         To mimic the affects of turbulent eddies, we have use these algorithms. 
!         They are quite varied, a bit hacky, and have the dangerous property of
!         perhaps violating local CFL conditions which will mean that particles
!         passing through domain boundaries will be forever lost.  A particle
!         must first be detected inside a domain's halo element to know that it should
!         be sent to the neighboring processor who contains that particular element
!         in their interior.  If it jumps right into the other processors interior.
!         It will be lost.  There are ways of dealing with this, of course, but I don't 
!         suspect it will happen anytime soon.
!
!     2.) Vertical random walk. 
!         This will be easier to implement, but the wide variety of models (read: 
!         a lot of parameters) will mean a separate inputfile may be necessary
!         for their specification.  Probably coming soon.
!    
!  =======Notes====
!  Note:  there is no specific array holding the particles because
!         they may traverse from one processor subdomain to another.
!         The particles are stored in the linked list 'particle_list'.
!         Particles may be added or deleted from the list (see linklist.F)
!         To access a particle property, for example, id, loop through the
!         linked list as follows:
!
! 1.)  lp  => particle_list%first%next
! 2.)  do
! 3.)    if(.not. associated(lp) ) exit      
! 4.)      write(*,*)'this is the current particles id: ',lp%v%id   
! 5.)     lp => lp%next                     
! 6.)  end do
!
!             A list is comprised of 0 or more nodes
! Line 1 sets the pointer to the first node in the list. (lp)
! Line 2 starts an infinite loop
! Line 3 exits when the end of the list is reached
! Line 4 shows example usage:  note lp is the node and lp%v is the particle
! Line 5 move to the next node, without this line you have an infinite loop 
! Line 6 is the end of the loop
!
! Reference: Object Oriented Programming via Fortran 90/95  , Ed Akin
!
! Note, you must include NETCDF_IO in the Makefile to build the Lagrangian library
! This will require you to build the netcdf librarie for your system.
! This is a result of a code modernization effort and is part of conspiracy 
!   against 'old school' programmers.  Well, maybe it is a conspiracy, but
!   tough luck.
!
! Runtime Parameters (to be placed in file:  casename_run.dat):
! Must be capital letters.
!
!  LAG_INTERVAL:  output interval in model hours
!  LAG_INPFILE:   input file name (will read INPDIR directory) 
!  LAG_RESFILE:   restart file (will dump to executable directory)
!  LAG_OUTFILE:   output file  (will write to OUTDIR/netcdf directory)
!  LAG_SCAL:      scalar variable to include with output
!  LAG_COLD_START logical (T/F) giving startup type
!
!  LAG_COLD_START = T
!       --> a new outputfile will be created
!       --> initial positions will be output to outputfile 
!       --> dumps are made every subsequent LAG_INTERVAL to outputfile
!
!  LAG_COLD_START = F (assumes input file generated by restart dump)
!       --> new outputfile will not be created
!       --> time of last dump is obtained from input file
!       --> stack # of last dump obtained from input file 
!       --> next dump made to stack# + 1 at time of last dump + LAG_INTERVAL 
!       --> initial positions not dumped to file (the assumed initial position
!           would be the first output in the file.
!
!  LAG_SCAL
!       current choices
!       s1   (salinity at nodes)
!       t1   (temperature at nodes)
!       rho1 (density at nodes)
!       km   (eddy viscosity)
!       kh   (vertical eddy diffusivity)
!       user_defined
!          --> users must set variable names in get_ud_scalnames
!          --> users must supply 3D field in add_scalar 
!          --> recompile after modification!
!       search the code for "users" to see where to modify
!
!  General information on restarting the code.
!    If you restart the code but don't want to interrupt your particle output 
!    file:
!
!    1.) run the code.  It will dump a particle restart file (designated by
!        LAG_RESFILE in the casename_run.dat file), every IRESTART iterations.
!        It will also dump one at the end of the run.
!    2.) move the LAG_RESFILE to your INPDIR directory.  Rename if you want.
!        Let's say you name it brainless_coding.nc
!    3.) modify the following parameters in your runtime file 
!           RESTART = hot_start  (may already be set to this)
!           LAG_COLD_START = F
!           LAG_INPFILE = brainless_coding.nc 
!    4.) restart the code.
!
!    The restart file you used as your LAG_INPFILE contains information about
!    when the last output was made and the stack# of the output so that the  
!    Lagrangian tracking module knows when to make the last dump and what
!    stack# to use in the netcdf file.  By setting LAG_COLD_START=F, the code
!    knows not to create a new outputfile, but rather use the othe one.
!
!    If you do something like restart the code but restart with a different
!    number of particles but use LAG_COLD_START = F, I cannot be responsible
!    for what might happen.  In simulations performed during the  extensive 
!    testing procedure used for every line of FVCOM code, we found that results
!    can range from crashed code to broken keyboards.
!
!------------------------------------------------------------------------------
  use mod_prec
  use linked_list
  use particle_class
  use netcdf

  implicit none

  private

  public :: set_lag
  public :: lag_update
  public :: dump_lag_restart
  
  !===========================================================
  ! ADDED FOR VISIT
  public :: LAG_ON
  public :: particle_list
  !===========================================================


  !global variables
  type(link_list) :: particle_list    !linked list holding particles
  integer :: nlag_gl = 0              !# global particles
  integer :: nlag = 0                 !# particles in local domain
  integer, parameter :: inlag  = 81   !initial position file unit #
  integer, parameter :: outlag = 82   !lagrangian output file #
  integer, allocatable :: ney_list(:) !list of domain neighbor processor ids
  logical lag_on                      !true if Lagrangian tracking active
  character(len=80 ) :: lag_inpfile   !particle location file name
  character(len=80 ) :: lag_resfile   !particle restart file name  
  character(len=120) :: lag_info_strng!information string from Lagrangian input file
  logical            :: cold_start    !true if cold starting

  !output control
  real(sp) :: lag_interval            !output interval in hours
  real(sp) :: t_last_dump = 0.0       !time (in hours) of last particle location dump
  character(len=80 ) :: lag_outfile   !particle output file name 

  !netcdf variables - input file
  integer :: lagin_nc_fid
  integer :: xin_vid
  integer :: yin_vid
  integer :: zin_vid
  integer :: groupin_vid
  integer :: markin_vid
  integer :: pathlengthin_vid
  integer :: tbegin_vid
  integer :: tendin_vid

  !netcdf variables - output file
  integer :: lag_nc_fid
  integer :: lag_nlag_did
  integer :: lag_time_did
  integer :: lag_time_vid
  integer :: lag_x_vid
  integer :: lag_y_vid
  integer :: lag_z_vid
  integer :: lag_sig_vid
  integer :: lag_s_vid
  integer :: dump_cnt = 0
  character(len=80) :: scalar_name 
  character(len=80) :: scalar_long_name 
  character(len=80) :: scalar_units 

  !Lagrangian particle associated scalar
  character(len=80) :: lag_scal_choice
  integer, parameter:: n_scal_choice = 6
  character(len=80) :: scal_choices(n_scal_choice) = (/"s1","t1","rho1","km","kh","user_defined"/)
  logical           :: scal_choice_valid



   contains !------------------------------------------------------------------!
            ! set_lag             :   read in initial particle locations       !
            !                     :   determine initial host elements          !
            ! lag_update          :   update particle location and write data  !
            ! traject             :   use 4 stage rk scheme to track particle  !
            ! interp_v            :   interpolate velocity at particle position!
            ! interp_elh          :   interpolate eta/h at particle position   !
            ! fhe                 :   fine element containint particle (driver)!
            ! fhe_robust          :   find element containing particle (robust)!
            ! fhe_quick           :   find element containing particle (quick) !
            ! exchange_particles  :   pass particles among processors          !
            ! isintriangle        :   true if particle is in given triangle    !
            ! output_lag          :   gather particles and call to output      !
            ! dump_particles_ncd  :   dump particle information to netcdf file !
            ! gather_particles    :   gather particles to master processor     !
            ! and many more all for a guaranteed low price!!                   !        
            ! -----------------------------------------------------------------!


!==============================================================================|
   subroutine set_lag 
!------------------------------------------------------------------------------|
!  read control parameters from input file
!  read initial positions from LAG_INPFILE
!  set remaining particle parameters (sigma location, pathlength, etc) 
!  report
!------------------------------------------------------------------------------|

  use all_vars
  use mod_inp
  use netcdf
  implicit none
  integer              :: i,ii,icnt,iscan,ierr
  logical              :: fexist,check
  character(len=120)   :: fname,temp
  type(particle)       :: p
  integer              :: dims(1)
  integer, allocatable :: tmp(:)
  type(link_node), pointer      :: lp,nlp


!----------------------------------------------------------
! initialize particle list (linked list structure)        
!----------------------------------------------------------
  particle_list = new_list()

!------------------------------------------------------------
!  check if Lagrangian particle tracking is activated 
!------------------------------------------------------------
                                                                                                             
   fname = "./"//trim(casename)//"_run.dat"
                                                                                                             
   iscan = scan_file(trim(fname),"LAG_ON",lval = lag_on)
   if(iscan /= 0)then
     write(ipt,*)'error reading lag_on: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!-------------------------------------------------------------
!  Read particle tracking parameters from runtime file       
!-------------------------------------------------------------
   if(lag_on)then

     !scan Lagrangian output interval (hours)
     iscan = scan_file(trim(fname),"LAG_INTERVAL",fscal = lag_interval)
     if(iscan /= 0)then
       write(ipt,*)'error reading lag_interval: ',iscan
       if(iscan == -2)then
         write(ipt,*)'variable not found in input file: ',trim(fname)
       end if
       call pstop
     end if

     !read in scalar string
     iscan = scan_file(trim(fname),"LAG_SCAL",cval = lag_scal_choice)
     if(iscan /= 0)then
       if(msr)write(ipt,*)'error reading lag_scal: ',iscan
       if(iscan == -2 .and. msr)then
         write(ipt,*)'variable not found in input file: ',trim(fname)
       end if
       call pstop
     end if
 
     !check validity of lag_scal_choice
     scal_choice_valid = .false.
     do i=1,n_scal_choice
       if(trim(lag_scal_choice) == trim(scal_choices(i))) scal_choice_valid = .true.
     end do
     if(msr .and. .not. scal_choice_valid)then
       write(ipt,*)'particle scalar choice "LAG_SCAL" not valid'
       write(ipt,*)'your choice: ',trim(lag_scal_choice)
       write(ipt,*)'valid options'
       do i=1,n_scal_choice
         write(ipt,*)i,' ',trim(scal_choices(i))
       end do
       call pstop
     end if

     !read startup type (cold or hot)
     iscan = scan_file(trim(fname),"LAG_COLD_START",lval = cold_start)
     if(iscan /= 0)then
       if(msr)write(ipt,*)'error reading LAG_COLD_START: ',iscan
       if(iscan == -2 .and. msr)then
         write(ipt,*)'variable not found in input file: ',trim(fname)
       end if
       call pstop
     end if



     !read in initial particle file name
     iscan = scan_file(trim(fname),"LAG_INPFILE",cval = lag_inpfile)   
     if(iscan /= 0)then
       if(msr)write(ipt,*)'error reading lag_inpfile: ',iscan
       if(iscan == -2 .and. msr)then
         write(ipt,*)'variable not found in input file: ',trim(fname)
       end if
       call pstop
     end if

     !read in particle restart file name 
     iscan = scan_file(trim(fname),"LAG_RESFILE",cval = lag_resfile)
     if(iscan /= 0)then
       if(msr)write(ipt,*)'error reading lag_resfile: ',iscan
       if(iscan == -2 .and. msr)then
         write(ipt,*)'variable not found in input file: ',trim(fname)
       end if
       call pstop
     end if

     !read in particle output file name 
     iscan = scan_file(trim(fname),"LAG_OUTFILE",cval = lag_outfile)
     if(iscan /= 0)then
       if(msr)write(ipt,*)'error reading lag_outfile: ',iscan
       if(iscan == -2 .and. msr)then
         write(ipt,*)'variable not found in input file: ',trim(fname)
       end if
       call pstop
     end if

     !make sure initial particle location file file exists
     fname = "./"//trim(inpdir)//"/"//trim(lag_inpfile)   
     inquire(file=trim(fname),exist=fexist)
     if(msr .and. .not.fexist)then
       write(ipt,*)'Lagrangian particle initial position file: '
       write(ipt,*) fname,' does not exist'
       write(ipt,*)'halting.....'
       call pstop
     end if

    !open file 
    ierr = nf90_open(trim(fname),nf90_nowrite,lagin_nc_fid)
    if(ierr /= nf90_noerr .and. msr)then
       write(ipt,*)'error opening', trim(fname)
       write(ipt,*)trim(nf90_strerror(ierr))
       call pstop
     endif

    !read dimensions
    temp = ""
    ierr = nf90_get_att(lagin_nc_fid,nf90_global,"info_string",temp)
    lag_info_strng = trim(temp)
    ierr = nf90_get_att(lagin_nc_fid,nf90_global,"number_particles",nlag_gl)
    ierr = nf90_get_att(lagin_nc_fid,nf90_global,"dump_counter",dump_cnt)    
    ierr = nf90_get_att(lagin_nc_fid,nf90_global,"t_last_dump",t_last_dump)    

    !sanity check
    if(nlag_gl < 1)then
      write(ipt,*)'# particles initial position file'
      write(ipt,*)'must be > 1'
      write(ipt,*)'is: ',nlag_gl
      stop
    endif

    !make sure output file already exists for hot start case
    if(.not. cold_start)then
      fname = "./"//trim(outdir)//"/"//"netcdf/"//trim(lag_outfile)
      inquire(exist=fexist,file=trim(fname))
      if(msr .and. .not.fexist)then
        write(ipt,*)'Lagrangian output file: ',trim(fname)        
        write(ipt,*)'does not exist'
        write(ipt,*)'for LAG_COLD_START = F case'
        write(ipt,*)'file must already exist'
        write(ipt,*)'halting.....'
        call pstop
      end if
    endif

  

    !report
    if(msr)then
    write(ipt,*  )'!'
    write(ipt,*)'!       Lagrangian tracking info'
    write(ipt,*  )'!'
    write(ipt,*)'!  # particle tracking   :  active'
    write(ipt,*)'!  # out interval (hrs)  :',lag_interval
    write(ipt,*)'!  # lag points          :',nlag_gl
    write(ipt,*)'!  input file            :',trim(lag_inpfile)
    write(ipt,*)'!  restart file          :',trim(lag_resfile)
    write(ipt,*)'!  output file           :',trim(lag_outfile)
    write(ipt,*)'!  information           :',trim(lag_info_strng)
    if(cold_start)then
      write(ipt,*)'!  startup type          :  cold start'
      write(ipt,*)'!  outputfile            :  new'
     else
      write(ipt,*)'!  startup type          :  hot start'
      write(ipt,*)'!  outputfile            :  old'
      write(ipt,*)'!  next data dump        :',dump_cnt +1       
      write(ipt,*)'!  time of last dump     :',t_last_dump  
    endif
    endif
 
    icnt = 0
    !read data for each particle and insert into link list 
    !reading:  x,y,z,mark,group,tbeg,tend,pathlength
    !setting:  [x,y,z,mark,group,tbeg,tend,pathlength]
    ierr = nf90_inq_varid(lagin_nc_fid,'x',xin_vid)
    ierr = nf90_inq_varid(lagin_nc_fid,'y',yin_vid)
    ierr = nf90_inq_varid(lagin_nc_fid,'z',zin_vid)
    ierr = nf90_inq_varid(lagin_nc_fid,'tbeg',tbegin_vid)
    ierr = nf90_inq_varid(lagin_nc_fid,'tend',tendin_vid)
    ierr = nf90_inq_varid(lagin_nc_fid,'group',groupin_vid)
    ierr = nf90_inq_varid(lagin_nc_fid,'mark',markin_vid)
    ierr = nf90_inq_varid(lagin_nc_fid,'pathlength',pathlengthin_vid)

    do i=1,nlag_gl
  
      !read particle (i) data
      dims(1) = i
      ierr = nf90_get_var(lagin_nc_fid,xin_vid,p%x(1),dims)
      ierr = nf90_get_var(lagin_nc_fid,yin_vid,p%x(2),dims)
      ierr = nf90_get_var(lagin_nc_fid,zin_vid,p%zloc,dims)
      ierr = nf90_get_var(lagin_nc_fid,groupin_vid,p%group,dims)
      ierr = nf90_get_var(lagin_nc_fid,markin_vid,p%mark,dims)
      ierr = nf90_get_var(lagin_nc_fid,pathlengthin_vid,p%pathlength,dims)
      ierr = nf90_get_var(lagin_nc_fid,tbegin_vid,p%tbeg,dims)
      ierr = nf90_get_var(lagin_nc_fid,tendin_vid,p%tend,dims)
  
      !adjust coordinate system
      p%x(1) = p%x(1) - vxmin
      p%x(2) = p%x(2) - vymin

      !set id
      p%id   = i

      !insert into linked list
      call node_insert(particle_list,p)
    end do
    ierr = nf90_close(lagin_nc_fid)

!------------------------------------------------------------
!   find element containint particle [elem, mark] 
!------------------------------------------------------------
    !find the element containing each particle: set [elem]
    if(msr)then
      write(ipt,*)'finding elements containing particles....'
    endif
    call fhe

    lp  => particle_list%first%next
    do
      if(.not. associated(lp) ) exit  !end of list, exit
        if(lp%v%elem > MElem .or. lp%v%elem == 0)then
          call node_delete(particle_list,lp%v,check)
          lp => particle_list%first%next
          else
          lp => lp%next
        endif
    end do

   
!------------------------------------------------------------
!  modify particle positions in z coordinate [zloc] 
!  ensure [-h < z > el]
!  if cold start, assume z = 0 means z = el if el /= 0
!------------------------------------------------------------

   !get surface elevation and bathymetry at particles (x,y)
   !set [el,h]
   call interp_elh

   !now calculate z from sigma,x,y
   lp  => particle_list%first%next
   do  
     if(.not. associated(lp) ) exit       !end of list, exit
      if(cold_start .and. lp%v%zloc == 0.0) lp%v%zloc = lp%v%el
      if(lp%v%zloc < -lp%v%h) lp%v%zloc = -lp%v%h
      lp => lp%next                              !set object
   end do

!------------------------------------------------------------
!  set sigma value of particle positions [x(3)]
!------------------------------------------------------------

   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit       !end of list, exit
      lp%v%x(3) = (lp%v%zloc - lp%v%el)/(lp%v%h + lp%v%el)
     lp => lp%next                              !set object
   end do

!------------------------------------------------------------
!   particle position at last time step [xn]       
!------------------------------------------------------------
    call shift_pos_list(particle_list)

!------------------------------------------------------------
!  set particle processor id  => [pid]                        
!------------------------------------------------------------

   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit   
      lp%v%pid = myid  
      lp => lp%next                    
   end do


!-------------------------------------------------------------
!  screen reporting of particle processor distribution 
!-------------------------------------------------------------


   !use mpi_gather to get domain distribution to master 
   allocate(tmp(nprocs))
   nlag = listsize(particle_list) !my # particles
   tmp(1) = nlag
   
# if defined (MULTIPROCESSOR)
   call mpi_gather(nlag,  1,mpi_integer,tmp,1,mpi_integer,0,mpi_comm_world,ierr)
#  endif

   !check for problems (sum of domain particles > global # particles)
   if(msr)then
   if(sum(tmp) > nlag_gl)then
     write(ipt,*)'some particles are considered to be in multiple domains'
     write(ipt,*)'nlag_gl',nlag_gl 
     write(ipt,*)'sum of particles from all domains',sum(tmp)
     write(ipt,*)'processor       # particles'
     do i=1,nprocs
       write(ipt,*)i,tmp(i)
     end do
     call pstop
    endif 
   endif


!-------------------------------------------------------------
!  set up scalar_name info for netcdf output
!-------------------------------------------------------------

    do i=1,n_scal_choice
      if(lag_scal_choice == scal_choices(i))then
      select case(i)

      case(1) !salinity
        scalar_name = 'salinity'
        scalar_long_name = 'salinity'
        scalar_units = 'ppt'
      case(2) !temperature
        scalar_name = 'temperature'
        scalar_long_name = 'temperature'
        scalar_units = 'C'
      case(3) !density
        scalar_name = 'density'
        scalar_long_name = 'density'
        scalar_units = 'kgm^-3'
      case(4) !vertical eddy viscosity - km
        scalar_name = 'km'
        scalar_long_name = 'Turbulent Eddy Viscosity'
        scalar_units = 'm^2 s^-1'
      case(5) !vertical eddy viscosity - kh
        scalar_name = 'kh'
        scalar_long_name = 'Turbulent Eddy Diffusivity'
        scalar_units = 'm^2 s^-1'
      case(6) !user_defined --> must set in get_ud_scalnames
        call get_ud_scalnames
      end select
      endif
   end do

   !output initial positions to file
   if(cold_start) call output_lag

   else  !lag tracking inactive
    if(msr) write(ipt,*)'!  # Lagrangian tracking :  inactive'
   end if

   end subroutine set_lag

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   subroutine lag_update
!==============================================================================|
!  update particle positions, calculate scalar fields and particle velocities  |
!==============================================================================|

   use all_vars
   use mod_inp
   implicit none
   integer i,ii,ierr
   type(link_node), pointer    :: lp 

!---------------------------------------------------
!  exit if Lagrangian tracking not active          
!---------------------------------------------------
   if(.not.lag_on) return


!---------------------------------------------------
!  update particle position                       
!---------------------------------------------------

   call traject(dti,thour,u,uf,v,vf,wts,wtts,h,el) 

!--------------------------------------------------
!  add scalar value at particle position
!--------------------------------------------------
   call add_scalar

!--------------------------------------------------
!  write particle position data to file            
!--------------------------------------------------
   call output_lag

!--------------------------------------------------
!  write particle position data to screen         
!   this will generate a lot of info, but might
!   be useful for debugging
!--------------------------------------------------
!   call print_list(particle_list)

   return
   end subroutine lag_update


!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   subroutine traject(deltat,time,u1,u2,v1,v2,w1,w2,hl,el) 
!==============================================================================|
!  integrate particle position from x0 to xn using velocity fields at time     |
!  t0 and time tn                                                              |
!  deltat:   time step between calls to traject (usually dti)                  |
!  time:     model time in hours
!  u1/v1/w1: velocity field (u,v,omega) at time t0                             |
!  u2/v2/w2: velocity field (u,v,omega) at time tn                             |
!  hl:       bathymetry                                                        | 
!  el:       free surface elevation at time tn                                 | 
!==============================================================================|
   use mod_prec
   use lims
   use control, only : iint
   implicit none
   real(sp),                     intent(in)    :: deltat
   real(sp),                     intent(in)    :: time 
   real(sp), dimension(0:MTElem,kb), intent(in)    :: u1,u2,v1,v2
   real(sp), dimension(0:NTNode,kb), intent(in)    :: w1,w2
   real(sp), dimension(0:NTNode),    intent(in)    :: hl,el 
!------------------------------------------------------------------------------|
!  n    : stage counter                                                        |
!  ns   : number of stages in explicit runga-kutta                             |
!  chix : stage function evaluation for x-velocity                             | 
!  pdx  : stage particle x-position                                            |
!  ul   : stage u velocity                                                     |
!  eps  : parameter defining depth of dry element                              |
!  dmax : maximum sigma depth                                                  |
!  a_rk : erk coefficients (a)                                                 |
!  b_rk : erk coefficients (b)                                                 |
!  c_rk : erk_coefficients (c)                                                 |
!------------------------------------------------------------------------------|
   integer                            :: ns
   integer,  parameter                :: mstage = 4
   real(sp), dimension(0:MTElem,kb)       :: ul,vl 
   real(sp), dimension(0:NTNode,kb)       :: wl 
   real(sp), parameter                :: eps  = 1.0e-5
   real(sp), parameter                :: dmax = -1.0_dp
   real(sp), dimension(3)             :: xtmp
   real(sp), parameter, dimension(mstage) :: a_rk = (/0.0_dp,0.5_dp,0.5_dp,1.0_dp/) 
   real(sp), parameter, dimension(mstage) :: b_rk = (/1.0_dp/6.0_dp,1.0_dp/3.0_dp, &
                                                      1.0_dp/3.0_dp,1.0_dp/6.0_dp/) 
   real(sp), parameter, dimension(mstage) :: c_rk = (/0.0_dp,0.5_dp,0.5_dp,1.0_dp/) 
   type(link_node), pointer :: p
!------------------------------------------------------------------------------|


!--set particle time step (0 if particle is not released or exceeded T_end)
   p  => particle_list%first%next
   do
     if(.not. associated(p) ) exit  !end of list, exit
      if(time >= p%v%tbeg .and. time <= p%v%tend .and. p%v%mark == 0)then
        p%v%deltat = deltat
      else
        p%v%deltat = 0.0_sp
      endif
      p => p%next                          !set object
   end do

!--initialize stage functional evaluations
   p  => particle_list%first%next
   do
     if(.not. associated(p) ) exit  !end of list, exit
     p%v%chi = 0.0_sp
     p => p%next                          !set object
   end do

!--loop over Runge-Kutta stages and calculate stage velocities 
   do ns=1,mstage

     !!particle position at stage n (x,y,sigma)
     p  => particle_list%first%next
     do
       if(.not. associated(p) ) exit  !end of list, exit
        p%v%x(1) = p%v%xn(1) + a_rk(ns)*p%v%deltat*p%v%chi(1,ns-1)  
        p%v%x(2) = p%v%xn(2) + a_rk(ns)*p%v%deltat*p%v%chi(2,ns-1)  
        p%v%x(3) = p%v%xn(3) + a_rk(ns)*p%v%deltat*p%v%chi(3,ns-1)  
        p => p%next                          !set object
     end do

     !!determine element location 
     call fhe

     !!shift particles to new domain if they have crossed boundary
#    if defined (MULTIPROCESSOR)
     call exchange_particles
#    endif

     !!determine el/h at stage position to convert to sigma
     call interp_elh

     !!adjust sigma position of particle at boundaries 
     p  => particle_list%first%next
     do
       if(.not. associated(p) ) exit  !end of list, exit
        p%v%x(3) = max(p%v%x(3),-(2.0+p%v%x(3))) !mirror bottom 
        p%v%x(3) = min(p%v%x(3),0.0_sp)          !don't penetrate free surf
        p => p%next                         !set object
     end do

     !!calculate velocity field for stage n using c_rk coefficients
     ul  = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2 
     vl  = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2 
     wl  = (1.0_sp-c_rk(ns))*w1 + c_rk(ns)*w2 


     call interp_v(ul,vl,wl)

     !!evaluate bathymetry and free surface height at stage ns particle position
     call interp_elh

     !!calculate u,v,w at stage ns
     p  => particle_list%first%next
     do
       if(.not. associated(p) ) exit  !end of list, exit
        p%v%chi(1,ns) = p%v%u
        p%v%chi(2,ns) = p%v%v
        p%v%chi(3,ns) = p%v%w/(p%v%h+p%v%el)         !delta_sigma/deltat = omega/d
        p => p%next                         !set object
     end do

     !!do not allow vertical motion in very shallow water
     p  => particle_list%first%next
     do
       if(.not. associated(p) ) exit  !end of list, exit
        if((p%v%h + p%v%el) < eps) p%v%chi(3,:) = 0.0_sp
        p => p%next                         !set object
     end do

  end do !stage loop

!-sum stage contributions to get updated particle positions-------------------!
!-store accumulated pathlength
   p  => particle_list%first%next
   do
     if(.not. associated(p) ) exit  !end of list, exit
      xtmp = p%v%xn
      do ns=1,mstage
        p%v%xn(1) = p%v%xn(1) + p%v%deltat*p%v%chi(1,ns)*b_rk(ns)
        p%v%xn(2) = p%v%xn(2) + p%v%deltat*p%v%chi(2,ns)*b_rk(ns)
        p%v%xn(3) = p%v%xn(3) + p%v%deltat*p%v%chi(3,ns)*b_rk(ns)
      end do
      p%v%pathlength = p%v%pathlength + sqrt( (p%v%xn(1)-xtmp(1))**2 + (p%v%xn(2)-xtmp(2))**2)
      p => p%next                         !set object
   end do

!-adjust sigma locs near top/bottom-------------------------------------------!
   p  => particle_list%first%next
   do
     if(.not. associated(p) ) exit  !end of list, exit
     p%v%xn(3) = max(p%v%xn(3),-(2.0+p%v%xn(3)))      !mirror bottom 
     p%v%xn(3) = min(p%v%xn(3),0.0_sp)               !don't penetrate free surf
     p => p%next                         !set object
   end do
 
!-update x location-----------------------------------------------------------!
   p  => particle_list%first%next
   do
     if(.not. associated(p) ) exit  !end of list, exit
     p%v%x(1) = p%v%xn(1)
     p%v%x(2) = p%v%xn(2)
     p%v%x(3) = p%v%xn(3)
     p => p%next                         !set object
   end do
    
!--evaluate the element
   call fhe

!--evaluate bathymetry and free surface height at updated particle position----!
   call interp_elh

!--adjust depth of updated particle positions and calculate depth--------------!
     p  => particle_list%first%next
     do
       if(.not. associated(p) ) exit       !end of list, exit
        if(p%v%mark == 0)then
          p%v%zloc = (p%v%x(3)*(p%v%h+p%v%el) + p%v%el) !particle depth (m)
        endif
        p => p%next                              !set object
     end do

   return
   end subroutine traject
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   subroutine interp_v(uin,vin,win) 
!==============================================================================|
!  obtain a linear interpolation of velocity field uin,vin,win at particle     |
!  locations in particle list 'particle_list'                                  |
!  returns:                                                                    |
!     p%u,p%v,p%w, velocity field at particle locations for each particle      |
!     in list                                                                  |
!									       |
!==============================================================================|
                                                                                                         
!------------------------------------------------------------------------------|
                                                                                                         
   use all_vars
   implicit none
   real(sp), intent(in) :: uin(0:MTElem,1:kb),vin(0:MTElem,1:kb),win(0:NTNode,1:kb)
!------------------------------------------------------------------------------!
   integer  :: np,i,e1,e2,e3,n1,n2,n3,k1,k2,k
   real(sp) :: w01,wx1,wy1,w02,wx2,wy2
   real(sp) :: x0c,y0c ,zp
   real(sp) :: dudx,dudy,dvdx,dvdy 
   real(sp) :: ue01,ue02,ve01,ve02,we01,we02
   real(sp) :: zf1,zf2
   logical all_found
   type(link_node), pointer :: p
   
!==============================================================================|


!===============================================================================!
!  determine velocity, bathymetry, and free surface height at (xp,yp,zp)        !
!===============================================================================!

   p  => particle_list%first%next

   do !loop over objects in list
     if(.not. associated(p) ) exit  !end of list, exit

     if(p%v%elem == 0)cycle                 !particle not in domain
     i   = p%v%elem 
     e1  = nbe(i,1)
     e2  = nbe(i,2)
     e3  = nbe(i,3)
     n1  = nv(i,1)
     n2  = nv(i,2)
     n3  = nv(i,3)
     x0c = p%v%x(1) - xc(i)
     y0c = p%v%x(2) - yc(i)
     zp  = p%v%x(3)

!----determine sigma layers above and below particle (for u/v velocities)--------!
     if(zp > zz1(i,1))then     !!particle near surface
       k1  = 1
       k2  = 1
       zf1 = 1.0_sp 
       zf2 = 0.0_sp
     else if(zp < zz1(i,kbm1)) then !!particle near bottom
       k1 = kbm1
       k2 = kbm1
       zf1 = 1.0_sp
       zf2 = 0.0_sp
     else
       k1  = int( (zz1(i,1)-zp)/dz1(i,1) ) + 1
       k2  = k1 + 1
       zf1 = (zp-zz1(i,k2))/dz1(i,1)
       zf2 = (zz1(i,k1)-zp)/dz1(i,1)
     end if
       
!----linear interpolation of particle velocity in sigma level above particle-----!
     k = k1 

     dudx = a1u(i,1)*uin(i,k)+a1u(i,2)*uin(e1,k)+a1u(i,3)*uin(e2,k)+a1u(i,4)*uin(e3,k)
     dudy = a2u(i,1)*uin(i,k)+a2u(i,2)*uin(e1,k)+a2u(i,3)*uin(e2,k)+a2u(i,4)*uin(e3,k)
     dvdx = a1u(i,1)*vin(i,k)+a1u(i,2)*vin(e1,k)+a1u(i,3)*vin(e2,k)+a1u(i,4)*vin(e3,k)
     dvdy = a2u(i,1)*vin(i,k)+a2u(i,2)*vin(e1,k)+a2u(i,3)*vin(e2,k)+a2u(i,4)*vin(e3,k)
     ue01 = uin(i,k) + dudx*x0c + dudy*y0c
     ve01 = vin(i,k) + dvdx*x0c + dvdy*y0c

!----linear interpolation of particle position in sigma level below particle-----!
     k = k2 

     dudx = a1u(i,1)*uin(i,k)+a1u(i,2)*uin(e1,k)+a1u(i,3)*uin(e2,k)+a1u(i,4)*uin(e3,k)
     dudy = a2u(i,1)*uin(i,k)+a2u(i,2)*uin(e1,k)+a2u(i,3)*uin(e2,k)+a2u(i,4)*uin(e3,k)
     dvdx = a1u(i,1)*vin(i,k)+a1u(i,2)*vin(e1,k)+a1u(i,3)*vin(e2,k)+a1u(i,4)*vin(e3,k)
     dvdy = a2u(i,1)*vin(i,k)+a2u(i,2)*vin(e1,k)+a2u(i,3)*vin(e2,k)+a2u(i,4)*vin(e3,k)
     ue02 = uin(i,k) + dudx*x0c + dudy*y0c
     ve02 = vin(i,k) + dvdx*x0c + dvdy*y0c

!----interpolate particle velocity between two sigma layers----------------------!

     p%v%u = ue01*zf1 + ue02*zf2
     p%v%v = ve01*zf1 + ve02*zf2 

!----determine sigma layers above and below particle (for node-based omega)------!
     k1 = max(int( abs(zp)/dz1(i,1) ) + 1,kbm1)
     k2 = k1 + 1
     zf1 = (z1(i,k1)-zp)/dz1(i,1)
     zf2 = (zp-z1(i,k2))/dz1(i,1)

!----linear interpolation of sigma velocity in sigma layer above particle--------!
     k = k1
     w01  = aw0(i,1)*win(n1,k)+aw0(i,2)*win(n2,k)+aw0(i,3)*win(n3,k)
     wx1  = awx(i,1)*win(n1,k)+awx(i,2)*win(n2,k)+awx(i,3)*win(n3,k)
     wy1  = awy(i,1)*win(n1,k)+awy(i,2)*win(n2,k)+awy(i,3)*win(n3,k)
     we01 = w01 + wx1*x0c + wy1*y0c

!----linear interpolation of sigma velocity in sigma layer below particle--------!
     k = k2
     w02  = aw0(i,1)*win(n1,k)+aw0(i,2)*win(n2,k)+aw0(i,3)*win(n3,k)
     wx2  = awx(i,1)*win(n1,k)+awx(i,2)*win(n2,k)+awx(i,3)*win(n3,k)
     wy2  = awy(i,1)*win(n1,k)+awy(i,2)*win(n2,k)+awy(i,3)*win(n3,k)
     we02 = w02 + wx2*x0c + wy2*y0c

!----vertical interpolation of sigma velocity------------------------------------!
     p%v%w = we01*zf1 + we02*zf2 
     p => p%next                          !set object

   end do !!loop over particles


   return
   end subroutine interp_v
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   subroutine fhe
!==============================================================================|
!  determine element containing point (xp,yp,zp)                                !
!    first try a quick algorithm searaching only last known element and 
!    surrouding elements.  Then, search the whole domain by first searching
!    nearest element (based on centroid) and moving out
!===============================================================================!

   use all_vars, only : myid
   implicit none 
   logical ::  all_found  = .false.

   call set_not_found(particle_list)
   call fhe_quick(all_found)
   if(.not. all_found) call fhe_robust

   end subroutine fhe

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   subroutine interp_elh
!==============================================================================|
!  given a set of pts (xp,yp,zp) of size npts, obtain a linear interpolation   |
!  of the provided free surface/bathymetry (hin,ein) at these points           |
!                                                                              |
!  sets:                                                                       |
!     h(bathymetry at xp,yp) and el (free surface elevation at xp,yp)          |
!                                                                              |
!==============================================================================|

!------------------------------------------------------------------------------|
                                                                                                         
   use all_vars
   implicit none
!------------------------------------------------------------------------------!
   integer  :: np,i,n1,n2,n3,k1,k2,k
   real(sp) :: h0,hx,hy,e0,ex,ey
   real(sp) :: x0c,y0c 
   type(link_node), pointer :: p
   
!==============================================================================|


!===============================================================================!
!  linearly interpolate free surface height and bathymetry                      !
!===============================================================================!

   p  => particle_list%first%next

   do !loop over objects in list
     if(.not. associated(p) ) exit  !end of list, exit

     if(p%v%elem == 0)cycle
     i  = p%v%elem      
     n1  = nv(i,1)
     n2  = nv(i,2)
     n3  = nv(i,3)
     x0c = p%v%x(1) - xc(i)
     y0c = p%v%x(2) - yc(i)

!----linear interpolation of bathymetry------------------------------------------!
     h0 = aw0(i,1)*h(n1)+aw0(i,2)*h(n2)+aw0(i,3)*h(n3)
     hx = awx(i,1)*h(n1)+awx(i,2)*h(n2)+awx(i,3)*h(n3)
     hy = awy(i,1)*h(n1)+awy(i,2)*h(n2)+awy(i,3)*h(n3)
     p%v%h = h0 + hx*x0c + hy*y0c
                                                                                                                          
!----linear interpolation of free surface height---------------------------------!
     e0 = aw0(i,1)*el(n1)+aw0(i,2)*el(n2)+aw0(i,3)*el(n3)
     ex = awx(i,1)*el(n1)+awx(i,2)*el(n2)+awx(i,3)*el(n3)
     ey = awy(i,1)*el(n1)+awy(i,2)*el(n2)+awy(i,3)*el(n3)
     p%v%el = e0 + ex*x0c + ey*y0c

     p => p%next                          !set object
   end do !!loop over objects in linked list 


   return
   end subroutine interp_elh
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!==============================================================================|
   subroutine interp_scalar(n_nodes,n_lay,field)
!==============================================================================|
! interpolate a scalar field p%s at the current particle location              |
!==============================================================================|

!------------------------------------------------------------------------------|
                                                                                                         
   use all_vars, only : nv,aw0,awx,awy,xc,yc,kbm1,zz1
   implicit none
   integer , intent(in) :: n_nodes
   integer , intent(in) :: n_lay
   real(sp), intent(in) :: field(n_nodes,n_lay)
   type(link_node), pointer :: p
   integer                  :: i,n1,n2,n3,k1,k2,k
   real(sp)                 :: x0c,y0c,s0,sx,sy,alpha,s_lower,s_upper
   real(sp)                 :: dsig,sigloc
   
!==============================================================================|


!===============================================================================!
!  linearly interpolate free surface height and bathymetry                      !
!===============================================================================!

   p  => particle_list%first%next

   do !loop over objects in list
     if(.not. associated(p) ) exit  !end of list, exit

     if(p%v%elem == 0)cycle

     !element location (i) and surrounding nodes (n1,n2,n3)
     i  = p%v%elem      
     n1  = nv(i,1)
     n2  = nv(i,2)
     n3  = nv(i,3)

     !offset from element center
     x0c = p%v%x(1) - xc(i)
     y0c = p%v%x(2) - yc(i)

     !determine location in vertical grid and interpolation coefficients
     sigloc = p%v%x(3)

     !initialize for surface particle
     k1 = 1
     k2 = 1
     alpha = .5

     !bottom
     if(sigloc < zz1(i,kbm1))then
       k1 = kbm1
       k2 = kbm1
       alpha = .5
     else !intermediate 
       do k=1,kbm1-1
         if(sigloc  < zz1(i,k) .and. sigloc >= zz1(i,k+1) )then
           k1 = k
           k2 = k+1 
           dsig = zz1(i,k1)-zz1(i,k2)
           alpha = (zz1(i,k)-sigloc)/dsig
         endif
       end do
     endif

     !interpolate to xp,yp in layer below 
     s0 = aw0(i,1)*field(n1,k1)+aw0(i,2)*field(n2,k1)+aw0(i,3)*field(n3,k1)
     sx = awx(i,1)*field(n1,k1)+awx(i,2)*field(n2,k1)+awx(i,3)*field(n3,k1)
     sy = awy(i,1)*field(n1,k1)+awy(i,2)*field(n2,k1)+awy(i,3)*field(n3,k1)
     s_upper = s0 + sx*x0c + sy*y0c

     !interpolate to xp,yp in layer above
     s0 = aw0(i,1)*field(n1,k2)+aw0(i,2)*field(n2,k2)+aw0(i,3)*field(n3,k2)
     sx = awx(i,1)*field(n1,k2)+awx(i,2)*field(n2,k2)+awx(i,3)*field(n3,k2)
     sy = awy(i,1)*field(n1,k2)+awy(i,2)*field(n2,k2)+awy(i,3)*field(n3,k2)
     s_lower = s0 + sx*x0c + sy*y0c

     !interpolate between layers
     p%v%s = (alpha)*s_lower + (1.0-alpha)*s_upper

     p => p%next                          !set object
   end do !!loop over objects in linked list 

   return
   end subroutine interp_scalar
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   subroutine fhe_robust
!==============================================================================|
!  find home element for points (x,y)                                          |
!  search nearest element to progressively further elements. updates lagrangian|  
!  component "el_host" and marks lagrangian component "ifound"     with 1 if 
!  found.  returns logical variable "all_found" if all lagrangian variables    |
!  have a known host element.  the host element may have been found prior to   |
!  entry in this routine                                                       |
!==============================================================================|


!------------------------------------------------------------------------------|

   use all_vars
   implicit none
   
!------------------------------------------------------------------------------|
   integer i,min_loc
   real(sp), dimension(1:MTElem,1) :: radlist
   real(sp), dimension(3) :: xtri,ytri
   real(sp) :: xlag,ylag,radlast
   integer  :: locij(2)
   type(link_node), pointer :: p

!==============================================================================|
   
   p  => particle_list%first%next

   do !loop over objects in list

     !set object
     if(.not. associated(p) ) exit
! David added check mark!
     if(.not.p%v%found .and. p%v%mark == 0 )then   
     xlag  = p%v%x(1) 
     ylag  = p%v%x(2)
     radlist(1:MTElem,1) = sqrt((xc(1:MTElem)-xlag)**2 + (yc(1:MTElem)-ylag)**2)
     radlast = -1.0_sp
in:  do while(.true.)
       locij   = minloc(radlist,radlist>radlast)
       min_loc = locij(1)
       if(min_loc == 0) then
         exit in
       end if
       xtri    = vx(nv(min_loc,1:3)) 
       ytri    = vy(nv(min_loc,1:3)) 
       radlast = radlist(min_loc,1)
       if(isintriangle(xtri,ytri,xlag,ylag))then
         p%v%found  = .true. 
         p%v%elem   = min_loc
         exit in 
       end if
       radlast = radlist(min_loc,1)
     end do in

     endif  !.not. found
     if(.not.p%v%found) p%v%mark = 1

     p => p%next
   end do


   return
   end subroutine fhe_robust
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   subroutine fhe_quick(all_found)
!==============================================================================|
!  determine which element a list of particles reside in by searching          |
!  neighboring elements.  updates "inelem" component of lagrangian particle    |  
!  type and updates logical array "elem_found" flagging whether the host       |
!  has been found							       |
!==============================================================================|


!------------------------------------------------------------------------------|

   use all_vars
   implicit none
   logical, intent(inout) :: all_found

   integer i,j,k,ilast,iney,ncheck
   real(sp), dimension(3) :: xlast,ylast,xney,yney
   real(sp) :: xlag,ylag
   type(link_node), pointer :: p

!==============================================================================|
   

   p  => particle_list%first%next

   do !loop over objects in list

     !set object
     if(.not. associated(p) ) exit
      
     ilast = p%v%elem
! David added check mark!
!     if(.not. p%v%found .and. ilast /= 0)then 
     if(.not. p%v%found .and. ilast /= 0 .and. p%v%mark ==0 )then 
     xlag  = p%v%x(1) 
     ylag  = p%v%x(2) 
     xlast = vx(nv(ilast,1:3))
     ylast = vy(nv(ilast,1:3))
     if(isintriangle(xlast,ylast,xlag,ylag))then       !!particle remains in element
       p%v%found = .true.
     else                                             !!check neighbors
outer: do j=1,3
         ncheck = nv(ilast,j)
         do k=1,ntve(nv(ilast,j))
           iney = nbve(ncheck,k) 
           xney = vx(nv(iney,1:3))
           yney = vy(nv(iney,1:3))
           if(isintriangle(xney,yney,xlag,ylag))then
             p%v%found  =  .true. 
             p%v%elem   = iney 
             exit outer
           end if
         end do
       end do outer
     end if
     endif !.not. found

     p => p%next
   end do


   return
   end subroutine fhe_quick 
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   logical function isintriangle(xt,yt,x0,y0) 
!==============================================================================|
!  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    |
!  using algorithm used for scene rendering in computer graphics               |
!  algorithm works well unless particle happens to lie in a line parallel      |
!  to the edge of a triangle.                                                  |
!  This can cause problems if you use a regular grid, say for idealized        |
!  modelling and you happen to see particles right on edges or parallel to     |
!  edges.                                                                      |
!==============================================================================|

   use mod_prec
   implicit none
   real(sp), intent(in) :: x0,y0
   real(sp), intent(in) :: xt(3),yt(3)
   real(sp) :: f1,f2,f3
   real(sp) :: x1(2)
   real(sp) :: x2(2)
   real(sp) :: x3(2)
   real(sp) :: p(2)

!------------------------------------------------------------------------------|

   isintriangle = .false.  
!   p(1) = x0
!   p(2) = x0
!   x1(1) = xt(1)
!   x1(2) = yt(1)
!   x2(1) = xt(2)
!   x2(2) = yt(2)
!   x3(1) = xt(3)
!   x3(2) = yt(3)
!
!   if(sameside(p,x1,x2,x3).and.sameside(p,x2,x1,x3).and. &
!      sameside(p,x3,x1,x2)) isintriangle = .true. 
   if(y0 < minval(yt) .or. y0 > maxval(yt)) then
     isintriangle = .false.
     return
   endif
   if(x0 < minval(xt) .or. x0 > maxval(xt)) then
     isintriangle = .false.
     return
   endif

   f1 = (y0-yt(1))*(xt(2)-xt(1)) - (x0-xt(1))*(yt(2)-yt(1))
   f2 = (y0-yt(3))*(xt(1)-xt(3)) - (x0-xt(3))*(yt(1)-yt(3))
   f3 = (y0-yt(2))*(xt(3)-xt(2)) - (x0-xt(2))*(yt(3)-yt(2))
   if(f1*f3 >= 0.0_sp .and. f3*f2 >= 0.0_sp) isintriangle = .true.

   return
   end function isintriangle
!==============================================================================|

  function sameside(p1,p2,a,b) result(value)
     real(sp), intent(in) :: p1(2)
     real(sp), intent(in) :: p2(2)
     real(sp), intent(in) :: a(2)
     real(sp), intent(in) :: b(2)
     logical value
     real(sp) :: cp1,cp2
  
     cp1 = (b(1)-a(1))*(p1(2)-a(2)) - (b(2)-a(2))*(p1(1)-a(1))
     cp2 = (b(1)-a(1))*(p2(2)-a(2)) - (b(2)-a(2))*(p2(1)-a(1))
  
     value = .false.
     if(cp1*cp2 >= 0) value = .true.

  end function sameside

!==============================================================================|
! exchange particles across processor domains                                  |
!==============================================================================|
# if defined (MULTIPROCESSOR)
  subroutine exchange_particles
!------------------------------------------------------------------------------|
  use lims,    only : myid ,nprocs, MElem
  use mod_par, only : el_pid, he_own, he_lst, elid
  use control, only : iint
  use all_vars, only : vx,vy,nv
  implicit none

  integer               ::ireqr(nprocs),ireqs(nprocs)
  real(sp), allocatable :: rbuf(:),sbuf(:)
  integer   stat(mpi_status_size),istatr(mpi_status_size,nprocs),ierr,j,n1,n2,ncnt
  integer   i,ifrom,ito,istag,irtag,trcv,tsnd,nvars,lbuf,nmsg,indx,lproc,nsze,rp
  integer   ibuf,recv_buf_size,send_buf_size
  integer   tri_local,tri_global

  integer, allocatable :: nrcv(:)
  integer, allocatable :: nsnd(:)
  type(particle)                :: p
  type(link_node), pointer      :: lp

  logical check
!------------------------------------------------------------------------------

  !--------------------------------------------
  !find particles in halo
  !set processor id to neighboring processor
  !we will send particle to that processor
  !--------------------------------------------
  lp  => particle_list%first%next
  do
    if(.not. associated(lp) ) exit 
      tri_local = lp%v%elem
      if(tri_local > MElem)then
       lp%v%pid  = he_own(tri_local-MElem) 
       lp%v%elem = he_lst(tri_local-MElem)   
      endif
      lp => lp%next
  end do

  !-------------------------------------------
  !determine # to send to each processor 
  !-------------------------------------------
  allocate(nsnd(nprocs)) ; nsnd = 0
  allocate(nrcv(nprocs)) ; nrcv = 0
  lp  => particle_list%first%next
  do
    if(.not. associated(lp) ) exit  
      if(lp%v%pid /= myid) nsnd(lp%v%pid) = nsnd(lp%v%pid) + 1
      lp => lp%next
  end do
 

  !-------------------------------------------
  !determine # to recv from each processor 
  !-------------------------------------------
  call recv_count(nsnd,nrcv)

  !-------------------------------------------
  !allocate buffer space      
  !-------------------------------------------
  send_buf_size = max(par_type_size*sum(nsnd) , 1)
  recv_buf_size = max(par_type_size*sum(nrcv) , 1)
  allocate(sbuf(send_buf_size)) ; sbuf = 0.0_sp  
  allocate(rbuf(recv_buf_size)) ; rbuf = 0.0_sp  

  !-------------------------------------------
  !Post non-blocking receives from neighbors
  !-------------------------------------------

  trcv = 0
  rp   = 1
  do i=1,nprocs
    if(nrcv(i) > 0)then
      ifrom = i-1
      irtag = i*1000
      trcv  = trcv + 1
      lbuf  = par_type_size*nrcv(i) 
      call mpi_irecv(rbuf(rp),lbuf,mpi_f,ifrom,irtag,mpi_comm_world,ireqr(trcv),ierr)
    end if
    rp = rp + nrcv(i)*par_type_size
  end do

  !-------------------------------------------
  !send data to neighbors
  !-------------------------------------------
  tsnd = 0
  ibuf = 1
  do i=1,nprocs
    lbuf = par_type_size*nsnd(i) 
    if(lbuf == 0 .or. myid == i)cycle

    !load send array
    ncnt = ibuf
    lp   => particle_list%first%next
    do
      if(.not. associated(lp) ) exit  !end of list, exit
        if(lp%v%pid == i)then
          call load_float(sbuf(ibuf),lp%v)
          ibuf = ibuf + par_type_size
        endif
        lp => lp%next                          !set object
    end do

    !send
    tsnd  = tsnd + 1
    ito   = i-1
    istag = myid*1000
    call mpi_isend(sbuf(ncnt),lbuf,mpi_f,ito,istag,mpi_comm_world,ireqs(tsnd),ierr)
  end do

  !--------------------------------------------
  !loop over procs until a message is received
  !--------------------------------------------
  do nmsg = 1,trcv 
    rp = 1
    call mpi_waitany(trcv,ireqr,indx,stat,ierr)
    lproc = stat(mpi_source) +1 
    do i=1,lproc-1
      rp = rp + nrcv(i)*par_type_size
    end do
    do i=1,nrcv(lproc)
      call load_particle(rbuf(rp),p)
      p%elem = elid(p%elem) 
      call node_insert(particle_list,p)
      rp = rp + par_type_size
    end do
  end do

  !--------------------------------------------
  !wait for completion of non-blocking sends  
  !--------------------------------------------

  call mpi_waitall(tsnd,ireqs,istatr,ierr)
  call mpi_waitall(trcv,ireqr,istatr,ierr)
  deallocate(sbuf)
  deallocate(rbuf)

  !--------------------------------------------
  !delete particles which have exited domain
  !--------------------------------------------
!  lp  => particle_list%first%next
!  do
!    if(.not. associated(lp) ) exit  !end of list, exit
!      if(lp%v%pid /= myid)then
!        call node_delete(particle_list,lp%v,check)  
!      endif
!      lp => lp%next                       
!  end do
   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit  !end of list, exit
       if(lp%v%pid /= myid)then
         call node_delete(particle_list,lp%v,check)
         lp => particle_list%first%next
         else
         lp => lp%next
       endif
   end do


  !--------------------------------------------
  !reset home element
  !--------------------------------------------
  call fhe


  end subroutine exchange_particles


!==============================================================================|
! count particle receives across processor boundaries                          |
!==============================================================================|
  subroutine recv_count(nsnd,nrcv)
!------------------------------------------------------------------------------|
# if defined (MULTIPROCESSOR)
  use mod_par, only : MHE, he_own
# endif
  use lims,    only : myid ,nprocs, MElem
  implicit none
  integer, intent(in ) :: nsnd(nprocs)
  integer, intent(out) :: nrcv(nprocs)

  integer               ::ireqr(nprocs),ireqs(nprocs)
  real(sp), allocatable :: rbuf(:),sbuf(:)
  integer   stat(mpi_status_size),istatr(mpi_status_size,nprocs),ierr,j,n1,n2,ncnt
  integer   i,ifrom,ito,istag,irtag,trcv,tsnd,nvars,lbuf,nmsg,indx,lproc,nsze,rp
  integer   ibuf,recv_buf_size,send_buf_size
  integer   tri_local,tri_global

  type(particle)                :: p
  type(link_node), pointer      :: lp
  logical :: first_entry = .true.

  logical check
!------------------------------------------------------------------------------

  !---------------------------------------------
  !if first entry, determine adjacent processors
  !--------------------------------------------
  if(first_entry)then
    allocate(ney_list(nprocs)) ; ney_list = 0
    do i=1,MHE
      ney_list(he_own(i)) = 1
    end do
    first_entry = .false.
  endif
    
  !-------------------------------------------
  !Post non-blocking receives from neighbors
  !-------------------------------------------

  trcv = 0
  do i=1,nprocs
    if(ney_list(i) > 0)then
      ifrom = i-1
      irtag = i*1000
      trcv  = trcv + 1
      call mpi_irecv(nrcv(i),1,mpi_integer,ifrom,irtag,mpi_comm_world,ireqr(trcv),ierr)
    end if
  end do

  !-------------------------------------------
  !send data to neighbors
  !-------------------------------------------
  tsnd = 0
  do i=1,nprocs
    if(ney_list(i) > 0)then
      tsnd  = tsnd + 1
      ito   = i-1
      istag = myid*1000
      call mpi_isend(nsnd(i),1,mpi_integer,ito,istag,mpi_comm_world,ireqs(tsnd),ierr)
    endif
  end do

  !--------------------------------------------
  !loop over procs until a message is received
  !--------------------------------------------
  do nmsg = 1,trcv 
    call mpi_waitany(trcv,ireqr,indx,stat,ierr)
  end do

  !--------------------------------------------
  !wait for completion of non-blocking sends  
  !--------------------------------------------

  call mpi_waitall(tsnd,ireqs,istatr,ierr)


  end subroutine recv_count         
# endif

!==============================================================================|
! output_lag  
!==============================================================================|
  subroutine output_lag 
   use all_vars
   use mod_ncdio, only : institution
   use mod_utils, only : get_timestamp
   implicit none
   logical, save :: first = .true.
   character(len=120) :: fname
   character(len=100) :: timestamp
   character(len=4  ) :: proc_string 
   integer            :: ierr
   integer            :: dynm2d(2)
   integer            :: dynmtime(1)


   !--------------------------------------------------------
   !open file
   !   if netcdf, write header information
   !--------------------------------------------------------

   if(first .and. cold_start)then
     first = .false.
     t_last_dump = thour - 1.1*LAG_INTERVAL
     dump_cnt = 0

     !open netcdf file
     if(msr)then
     fname = "./"//trim(outdir)//"/"//"netcdf/"//trim(lag_outfile)
     ierr = nf90_create(trim(fname),nf90_64bit_offset,lag_nc_fid)
     if(ierr /= nf90_noerr)then
       write(ipt,*)'error creating', trim(fname)
       write(ipt,*)trim(nf90_strerror(ierr))
     endif

     !global attributes
     call get_timestamp(timestamp)
     ierr = nf90_put_att(lag_nc_fid,nf90_global,"title"      ,trim(casetitle))
     ierr = nf90_put_att(lag_nc_fid,nf90_global,"institution",trim(institution))
     ierr = nf90_put_att(lag_nc_fid,nf90_global,"source"     ,'Particle'//trim(fvcom_version))
     ierr = nf90_put_att(lag_nc_fid,nf90_global,"history"    ,trim(timestamp))
     ierr = nf90_put_att(lag_nc_fid,nf90_global,"references" ,trim(fvcom_website))
     write(proc_string,'(I4.4)') nprocs
     ierr = nf90_put_att(lag_nc_fid,nf90_global,"num_processes" ,nprocs)

     !dimensions
     ierr = nf90_def_dim(lag_nc_fid,"nlag",nlag_gl,       lag_nlag_did)
     ierr = nf90_def_dim(lag_nc_fid,"time",nf90_unlimited,lag_time_did)
     dynm2d   = (/lag_nlag_did,lag_time_did/)
     dynmtime = (/lag_time_did/)

     !-----time------
     ierr = nf90_def_var(lag_nc_fid,"time",nf90_float,dynmtime,lag_time_vid)
     ierr = nf90_put_att(lag_nc_fid,lag_time_vid,"long_name","time")
     ierr = nf90_put_att(lag_nc_fid,lag_time_vid,"units","days")
     !-----x---------
     ierr = nf90_def_var(lag_nc_fid,"x",nf90_float,dynm2d,lag_x_vid)
     ierr = nf90_put_att(lag_nc_fid,lag_x_vid,"long_name","particle x position")
     ierr = nf90_put_att(lag_nc_fid,lag_x_vid,"units","m")
     !-----y---------
     ierr = nf90_def_var(lag_nc_fid,"y",nf90_float,dynm2d,lag_y_vid)
     ierr = nf90_put_att(lag_nc_fid,lag_y_vid,"long_name","particle y position")
     ierr = nf90_put_att(lag_nc_fid,lag_y_vid,"units","m")
     !-----z---------
     ierr = nf90_def_var(lag_nc_fid,"z",nf90_float,dynm2d,lag_z_vid)
     ierr = nf90_put_att(lag_nc_fid,lag_z_vid,"long_name","particle z position")
     ierr = nf90_put_att(lag_nc_fid,lag_z_vid,"units","m")
     !-----sig-------
     ierr = nf90_def_var(lag_nc_fid,"sigma",nf90_float,dynm2d,lag_sig_vid)
     ierr = nf90_put_att(lag_nc_fid,lag_sig_vid,"long_name","particle sigma position")
     ierr = nf90_put_att(lag_nc_fid,lag_sig_vid,"units","-")
     !-----scal------
     ierr = nf90_def_var(lag_nc_fid,trim(scalar_name),nf90_float,dynm2d,lag_s_vid)
     ierr = nf90_put_att(lag_nc_fid,lag_s_vid,"long_name",trim(scalar_long_name))
     ierr = nf90_put_att(lag_nc_fid,lag_s_vid,"units",trim(scalar_units))

     !enddef and close
     ierr = nf90_enddef(lag_nc_fid)
     ierr = nf90_close(lag_nc_fid)
     endif !msr
   endif


   !if we have proceeded LAG_INTERVAL hours since last output:
   if((thour - t_last_dump) > LAG_INTERVAL)then 

     !reset t_last_dump
     t_last_dump = thour

     !gather particles to global
#    if defined (MULTIPROCESSOR)
     call gather_particles
#    endif

     !output to file 
     call dump_particles_ncd     
   endif
 
  end subroutine output_lag

!==============================================================================|
! collect particles to global list in master only 
!     needed to dump particles to file
!==============================================================================|
# if defined (MULTIPROCESSOR)
  subroutine gather_particles
   use all_vars
   implicit none
   integer, allocatable :: tmp(:) 
   real(sp),allocatable :: send_space(:),rcve_space(:)
   integer   stat(mpi_status_size)
   integer :: ierr,isnd,ircv,ip,nsze,ibuf,i
   type(particle)                :: p
   type(link_node), pointer      :: lp

  
   !count number of particles in local domain => nlag
   nlag = listsize(particle_list)

   !send to master
   allocate(tmp(nprocs))
   tmp(1) = nlag

#  if defined (MULTIPROCESSOR)
   call mpi_gather(nlag,  1,mpi_integer,tmp,1,mpi_integer,0,mpi_comm_world,ierr)
#  endif

   !report processor particle distribution to screen
   if(msr)then
     write(ipt,*)'   particle/processor distribution'
     write(ipt,*)'processor       # particles'
     do i=1,nprocs
       write(ipt,*)i,tmp(i)
     end do
     write(ipt,*)'max: ',maxval(tmp),' total: ',sum(tmp),' min ',minval(tmp)
   endif

   !loop through other processors, receive particles, add to linked list
#  if defined (MULTIPROCESSOR)
   if(PAR)then

      !----------------------------------------------------
      !send data to master
      !----------------------------------------------------
      if(myid /= 1)then
        !allocate send space
        allocate(send_space(par_type_size*nlag)) ; send_space = 0
  
        !load send array
        ibuf = 1
        lp  => particle_list%first%next
        do
          if(.not. associated(lp) ) exit  !end of list, exit
            call load_float(send_space(ibuf),lp%v)  
            ibuf = ibuf + par_type_size
            lp => lp%next                          !set object
        end do
  
        !send to master processor
        isnd = myid+100
        call mpi_send(send_space,par_type_size*nlag,mpi_f,0,isnd,mpi_comm_world,ierr)
        deallocate(send_space)
        return
      endif

      !----------------------------------------------------
      !recv data and add to linked list
      !----------------------------------------------------

      do ip=2,nprocs
        nsze = (tmp(ip)*par_type_size)
        allocate(rcve_space(nsze))
        ircv = ip+100
        call mpi_recv(rcve_space,nsze,mpi_f,ip-1,ircv,mpi_comm_world,stat,ierr)
        ibuf = 1
        do i=1,tmp(ip)
          call load_particle(rcve_space(ibuf),p)
          call node_insert(particle_list,p)
          ibuf = ibuf + par_type_size
        end do
        deallocate(rcve_space)
      end do

   endif !PAR
#  endif
   deallocate(tmp)

  end subroutine gather_particles
# endif

!==============================================================================|
! dump particles to netcdf file 
!      dumps:  time, position (x,y,z) and associated scalar value
!==============================================================================|
  subroutine dump_particles_ncd
   use all_vars, only : vxmin,vymin,msr,iint,thour,time,outdir,casename,ipt
   implicit none
   type(link_node), pointer :: lp,nlp
   real(sp), allocatable, dimension(:) :: xtmp
   real(sp), allocatable, dimension(:) :: ytmp
   real(sp), allocatable, dimension(:) :: ztmp
   real(sp), allocatable, dimension(:) :: sigtmp
   real(sp), allocatable, dimension(:) :: stmp
   integer , allocatable, dimension(:) :: missing
   integer , allocatable, dimension(:) :: tmp 
   integer                             :: i,ii,ierr,nmissing
   integer                             :: varid
   type(particle)                      :: p
   integer, dimension(1)               :: dims
   integer, dimension(2)               :: dims2
   logical                             :: check
   character(len=120)                  :: fname


   if(.not.msr) return 

   !look for missing particles
   allocate(tmp(nlag_gl)) ; tmp = 0
   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit  !end of list, exit
     tmp(lp%v%id) = 1
     lp => lp%next                          !set object
   end do

   nmissing = nlag_gl - sum(tmp)  

   !find missing particle id's and add to blank particles to linked list
   if(nmissing > 0)then

     allocate(missing(nmissing))
     ii = 0
     do i=1,nlag_gl
       if(tmp(i) == 0)then
         ii = ii + 1
         missing(ii) = i
       endif
     end do

     do i=1,nmissing
       call zero_out(p)
       p%id = missing(i)
       p%x(1) = -100.
       p%x(2) = -100.
       call node_insert(particle_list,p)
     end do

   endif

   !update dump counter
   dump_cnt = dump_cnt + 1

   !check # particles
   if(listsize(particle_list) /= nlag_gl)then
     call print_list(particle_list)
     write(ipt,*)'number of particles being tracked: ',listsize(particle_list)
     write(ipt,*)'does not equal original number of particles: ',nlag_gl
     write(ipt,*)'stopping..........'
     call pstop
   endif

   !put particle data into temporary 1d arrays
   allocate(xtmp(nlag_gl)) ; xtmp = 0.0_sp
   allocate(ytmp(nlag_gl)) ; ytmp = 0.0_sp
   allocate(ztmp(nlag_gl)) ; ztmp = 0.0_sp
   allocate(sigtmp(nlag_gl)) ; sigtmp = 0.0_sp
   allocate(stmp(nlag_gl)) ; stmp = 0.0_sp


   ii = 0
   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit  !end of list, exit
     ii = ii + 1
     xtmp(ii) = lp%v%x(1) + vxmin
     ytmp(ii) = lp%v%x(2) + vymin
     sigtmp(ii) = lp%v%x(3) 
     ztmp(ii) = lp%v%zloc
     stmp(ii) = lp%v%s
     lp => lp%next                          !set object
   end do

   !dump to netcdf file
   fname = "./"//trim(outdir)//"/"//"netcdf/"//trim(lag_outfile)  
   ierr = nf90_open(fname,nf90_write,lag_nc_fid)
   if(ierr /= nf90_noerr)then
     write(ipt,*)'error opening', trim(fname)
     write(ipt,*)trim(nf90_strerror(ierr))
   endif
   dims(1)   = dump_cnt 
   
   ierr = nf90_inq_varid(lag_nc_fid,'time',varid)
   ierr = nf90_put_var(lag_nc_fid,varid,time,START=dims)
   if(ierr /= nf90_noerr)then
     write(ipt,*)'error dumping time to', trim(fname)
     write(ipt,*)trim(nf90_strerror(ierr))
     call pstop
   endif
   dims2(1) = 1
   dims2(2) = dump_cnt 
   ierr = nf90_inq_varid(lag_nc_fid,'x',varid)
   ierr = nf90_put_var(lag_nc_fid, varid, xtmp,START=dims2)
   ierr = nf90_inq_varid(lag_nc_fid,'y',varid)
   ierr = nf90_put_var(lag_nc_fid, varid, ytmp,START=dims2)
   ierr = nf90_inq_varid(lag_nc_fid,'z',varid)
   ierr = nf90_put_var(lag_nc_fid, varid, ztmp,START=dims2)
   ierr = nf90_inq_varid(lag_nc_fid,'sigma',varid)
   ierr = nf90_put_var(lag_nc_fid, varid, sigtmp,START=dims2)
   ierr = nf90_inq_varid(lag_nc_fid,trim(scalar_name),varid)
   ierr = nf90_put_var(lag_nc_fid, varid, stmp,START=dims2)

   !close up
   ierr = nf90_close(lag_nc_fid)


   !deallocate temporary arrays
   deallocate( xtmp )
   deallocate( ytmp )
   deallocate( ztmp )
   deallocate( sigtmp )
   deallocate( stmp )

   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit  !end of list, exit
       if(lp%v%pid /= 1)then
         call node_delete(particle_list,lp%v,check)
         lp => particle_list%first%next
         else
         lp => lp%next
       endif
   end do

  end subroutine dump_particles_ncd

  subroutine get_ud_scalnames
    implicit none

    !users, define your names, long_names, and units for the scalar field here
    scalar_name = 'generic_name'              !generally the variable name
    scalar_long_name = 'generic_long_name'    !generally the oceanographic name
    scalar_units = 'generic_units'            !units in mks: example:  m^2s^-1
  end subroutine

  subroutine add_scalar 
    use all_vars, only : s1,t1,rho1,km,kh,NNode,kbm1
!   users, you must include your module containig the scalar field here
!   use my_module, only : my_var  !add user module here
    implicit none
    real(sp) :: field(NNode,kbm1)
    integer :: i,k,ns

    do ns=1,n_scal_choice
      if(lag_scal_choice == scal_choices(ns))then
      select case(ns)

      case(1) !salinity
        do k=1,kbm1
        do i=1,NNode
          field(i,k) = s1(i,k)
        end do
        end do
      case(2) !temperature
        do k=1,kbm1
        do i=1,NNode
          field(i,k) = t1(i,k)
        end do
        end do
      case(3) !density
        do k=1,kbm1
        do i=1,NNode
          field(i,k) = rho1(i,k)
        end do
        end do
      case(4) !vertical eddy viscosity - km
        do k=1,kbm1
        do i=1,NNode
          field(i,k) = .5*(km(i,k)+km(i,k+1))
        end do
        end do
      case(5) !vertical eddy viscosity - kh
        do k=1,kbm1
        do i=1,NNode
          field(i,k) = .5*(kh(i,k)+kh(i,k+1))
        end do
        end do
      case(6) !user-defined --> users , you must include your data module and supply array 
        do k=1,kbm1
        do i=1,NNode
!          field(i,k) = my_var(i,k)
        end do
        end do
      end select

      endif
   end do

   call interp_scalar(NNode,kbm1,field)

  end subroutine add_scalar

  subroutine dump_lag_restart

   use all_vars, only : vxmin,vymin,msr,ipt
   implicit none
   type(link_node), pointer :: lp,nlp
   real(sp), allocatable, dimension(:) :: xtmp
   real(sp), allocatable, dimension(:) :: ytmp
   real(sp), allocatable, dimension(:) :: ztmp
   real(sp), allocatable, dimension(:) :: ptmp
   real(sp), allocatable, dimension(:) :: btmp
   real(sp), allocatable, dimension(:) :: etmp
   integer , allocatable, dimension(:) :: gtmp
   integer , allocatable, dimension(:) :: mtmp
   integer , allocatable, dimension(:) :: missing
   integer , allocatable, dimension(:) :: tmp 
   integer                             :: i,ii,ierr,nmissing
   type(particle)                      :: p
   integer, dimension(1)               :: dims
   logical                             :: check
   character(len=120)                  :: fname
   integer                             :: fid
   integer                             :: lagres_did 
   integer                             :: lagres_x_vid,lagres_y_vid,lagres_z_vid 
   integer                             :: lagres_p_vid,lagres_b_vid,lagres_e_vid 
   integer                             :: lagres_g_vid,lagres_m_vid


    !collect particles
#   if defined (MULTIPROCESSOR)
    call gather_particles
#   endif

   if(.not.msr) return 

   !look for missing particles
   allocate(tmp(nlag_gl)) ; tmp = 0
   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit  !end of list, exit
     tmp(lp%v%id) = 1
     lp => lp%next                          !set object
   end do

   nmissing = nlag_gl - sum(tmp)  

   !find missing particle id's and add to blank particles to linked list
   if(nmissing > 0)then

     allocate(missing(nmissing))
     ii = 0
     do i=1,nlag_gl
       if(tmp(i) == 0)then
         ii = ii + 1
         missing(ii) = i
       endif
     end do

     do i=1,nmissing
       call zero_out(p)
       p%id = missing(i)
       p%x(1) = -100.
       p%x(2) = -100.
       call node_insert(particle_list,p)
     end do

   endif

   !put particle data into temporary 1d arrays
   allocate(xtmp(nlag_gl))       ; xtmp = 0.0_sp
   allocate(ytmp(nlag_gl))       ; ytmp = 0.0_sp
   allocate(ztmp(nlag_gl))       ; ztmp = 0.0_sp
   allocate(ptmp(nlag_gl))       ; ptmp = 0.0_sp
   allocate(btmp(nlag_gl))       ; btmp = 0.0_sp
   allocate(etmp(nlag_gl))       ; etmp = 0.0_sp
   allocate(gtmp(nlag_gl))       ; gtmp = 0
   allocate(mtmp(nlag_gl))       ; mtmp = 0

   ii = 0
   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit  !end of list, exit
     ii = ii + 1
     xtmp(ii) = lp%v%x(1) + vxmin
     ytmp(ii) = lp%v%x(2) + vymin
     ztmp(ii) = lp%v%zloc
     ptmp(ii) = lp%v%pathlength
     btmp(ii) = lp%v%tbeg
     etmp(ii) = lp%v%tend
     gtmp(ii) = lp%v%group
     mtmp(ii) = lp%v%mark  
     lp => lp%next                          !set object
   end do

   !open netcdf file
   fname = trim(lag_resfile)
   write(ipt,*)'dumping Lagrangian data to restartfile: ',trim(lag_resfile)
   ierr = nf90_create(trim(fname),nf90_64bit_offset,fid)

   !global attributes 
   ierr = nf90_put_att(fid,nf90_global,"info_string",trim(lag_info_strng))
   ierr = nf90_put_att(fid,nf90_global,"dump_counter",dump_cnt)
   ierr = nf90_put_att(fid,nf90_global,"t_last_dump",t_last_dump)
   ierr = nf90_put_att(fid,nf90_global,"number_particles",nlag_gl)

   !dimensions
   ierr = nf90_def_dim(fid,"nparticles",nlag_gl, lagres_did)

   !-----x---------
   ierr = nf90_def_var(fid,"x",nf90_float,lagres_did,lagres_x_vid)
   ierr = nf90_put_att(fid,lagres_x_vid,"long_name","particle x position")
   ierr = nf90_put_att(fid,lagres_x_vid,"units","m")
   !-----y---------
   ierr = nf90_def_var(fid,"y",nf90_float,lagres_did,lagres_y_vid)
   ierr = nf90_put_att(fid,lagres_y_vid,"long_name","particle y position")
   ierr = nf90_put_att(fid,lagres_y_vid,"units","m")
   !-----z---------
   ierr = nf90_def_var(fid,"z",nf90_float,lagres_did,lagres_z_vid)
   ierr = nf90_put_att(fid,lagres_z_vid,"long_name","particle z position")
   ierr = nf90_put_att(fid,lagres_z_vid,"units","m")
   !-----tbeg------
   ierr = nf90_def_var(fid,"tbeg",nf90_float,lagres_did,lagres_b_vid)
   ierr = nf90_put_att(fid,lagres_b_vid,"long_name","particle release time")
   ierr = nf90_put_att(fid,lagres_b_vid,"units","hours")
   !-----tend------
   ierr = nf90_def_var(fid,"tend",nf90_float,lagres_did,lagres_e_vid)
   ierr = nf90_put_att(fid,lagres_e_vid,"long_name","particle freeze time")
   ierr = nf90_put_att(fid,lagres_e_vid,"units","hours")
   !-----pathlength
   ierr = nf90_def_var(fid,"pathlength",nf90_float,lagres_did,lagres_p_vid)
   ierr = nf90_put_att(fid,lagres_p_vid,"long_name","particle integrated path length")
   ierr = nf90_put_att(fid,lagres_p_vid,"units","m")
   !-----group      
   ierr = nf90_def_var(fid,"group",nf90_int,lagres_did,lagres_g_vid)
   ierr = nf90_put_att(fid,lagres_g_vid,"long_name","particle group")
   ierr = nf90_put_att(fid,lagres_g_vid,"units","-")
   !-----mark       
   ierr = nf90_def_var(fid,"mark",nf90_int,lagres_did,lagres_m_vid)
   ierr = nf90_put_att(fid,lagres_m_vid,"long_name","particle marker (0=in domain)")
   ierr = nf90_put_att(fid,lagres_m_vid,"units","-")

   !enddef and close
   ierr = nf90_enddef(fid)
   ierr = nf90_close(fid)

   !reopen and write data
   ierr = nf90_open(fname,nf90_write,fid)
   if(ierr /= nf90_noerr)then
     write(ipt,*)'error opening', trim(fname)
     write(ipt,*)trim(nf90_strerror(ierr))
   endif

   dims(1)   = 1 
   ierr = nf90_put_var(fid, lagres_x_vid, xtmp,START=dims)
   ierr = nf90_put_var(fid, lagres_y_vid, ytmp,START=dims)
   ierr = nf90_put_var(fid, lagres_z_vid, ztmp,START=dims)
   ierr = nf90_put_var(fid, lagres_p_vid, ptmp,START=dims)
   ierr = nf90_put_var(fid, lagres_b_vid, btmp,START=dims)
   ierr = nf90_put_var(fid, lagres_e_vid, etmp,START=dims)
   ierr = nf90_put_var(fid, lagres_g_vid, gtmp,START=dims)
   ierr = nf90_put_var(fid, lagres_m_vid, mtmp,START=dims)

   !close up
   ierr = nf90_close(fid)

   !deallocate temporary arrays
   deallocate( xtmp )
   deallocate( ytmp )
   deallocate( ztmp )
   deallocate( ptmp )
   deallocate( btmp )
   deallocate( etmp )
   deallocate( gtmp )
   deallocate( mtmp )

   !delete nodes that belong to other processors
   lp  => particle_list%first%next
   do
     if(.not. associated(lp) ) exit  !end of list, exit
       if(lp%v%pid /= 1)then
         call node_delete(particle_list,lp%v,check)
         lp => particle_list%first%next
         else
         lp => lp%next
       endif
   end do

  end subroutine dump_lag_restart   

# endif
end module mod_lag

