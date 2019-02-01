!==============================================================================!
!                LAGRANGIAN PARTICLE CLASS                                     !
!==============================================================================!
  
module particle_class
  use mod_prec ,only : sp
  integer, parameter :: par_type_size = 34  !id => deltat 

  type particle 
    integer  :: id            =0 !!global particle number
    integer  :: pid           =0 !!processor id number        
    integer  :: elem          =0 !!element containing particle
    integer  :: group         =0 !!element group id
    integer  :: mark          =0 !!=1 out of domain               
    real(sp) :: tbeg          =0.!!start time in seconds
    real(sp) :: tend          =0.!!end time in seconds
    real(sp) :: x(3)          =0.!!particle position
    real(sp) :: xn(3)         =0.!!last time step position 
    real(sp) :: el            =0.!!surface elevation at particle point
    real(sp) :: h             =0.!!bathymetry at particle point
    real(sp) :: u             =0.!!x-velocity at particle location
    real(sp) :: v             =0.!!y-velocity at particle location
    real(sp) :: w             =0.!!sigma-velocity at particle location
    real(sp) :: zloc          =0.!!z position of particle [m]            
    real(sp) :: s             =0.!!particle scalar
    real(sp) :: chi(3,4)      =0.!!Runge-Kutta stage contributions
    real(sp) :: deltat        =0.!!particle time step
    real(sp) :: pathlength    =0.!!particle integrated pathlength [m]
    logical  :: found         = .false.
  end type particle 


  interface operator (<)    !for sorting, insert
    module procedure less_than_particle ; end interface
  interface operator (==)   !for sorting or delete
    module procedure equal_to_particle  ; end interface

  contains

  function less_than_particle (particle_1,particle_2) result(Boolean)
    type(particle), intent(in) :: particle_1
    type(particle), intent(in) :: particle_2
    logical                    :: Boolean
    Boolean = particle_1%id < particle_2%id
  end function less_than_particle

  function equal_to_particle (particle_1,particle_2) result(Boolean)
    type(particle), intent(in) :: particle_1
    type(particle), intent(in) :: particle_2
    logical                    :: Boolean
    Boolean = particle_1%id == particle_2%id
  end function equal_to_particle

  subroutine screen_write(p,hprint)
    type(particle), intent(in) :: p
    logical,        intent(in) :: hprint
    if(hprint)then
      write(*,*)'id   group    x           y         z      elem  mark'
    endif
    write(*,'(2I5,3F10.2,3I5)')p%id,p%group,p%x(1),p%x(2),p%x(3),p%elem,p%mark,p%pid
  end subroutine screen_write

  subroutine shift_pos(p)
    type(particle), intent(inout) :: p
    p%xn = p%x
  end subroutine shift_pos 

  subroutine set_pathlength(p)
    type(particle), intent(inout) :: p
    p%pathlength = sqrt( (p%x(1)-p%xn(1))**2 + (p%x(2)-p%xn(2))**2)
  end subroutine set_pathlength

  !convert particle type data (p) to a vector of floats (a)
  subroutine load_float(a,p)
    real(sp),       intent(out) :: a(par_type_size)
    type(particle), intent(in ) :: p
    a(1 ) = float(p%id)    !!global particle number
    a(2 ) = float(p%pid )  !!processor id number
    a(3 ) = float(p%elem)  !!element containing particle
    a(4 ) = float(p%group) !!element group id
    a(5 ) = float(p%mark)  !!=1 out of domain
    a(6 ) = p%tbeg         !!start time in seconds
    a(7 ) = p%tend         !!end time in seconds
    a(8 ) = p%x(1)         !!x particle position
    a(9 ) = p%x(2)         !!y particle position
    a(10) = p%x(3)         !!z particle position
    a(11) = p%xn(1)        !!last x particle position
    a(12) = p%xn(2)        !!last y particle position
    a(13) = p%xn(3)        !!last z particle position
    a(14) = p%el           !!surface elevation at particle point
    a(15) = p%h            !!bathymetry at particle point
    a(16) = p%u            !!x-velocity at particle location
    a(17) = p%v            !!y-velocity at particle location
    a(18) = p%w            !!sigma-velocity at particle location
    a(19) = p%zloc         !!z location of particle [m]
    a(20) = p%s            !!particle associated scalar
    a(21) = p%chi(1,1)     !!stage contributions
    a(22) = p%chi(2,1)
    a(23) = p%chi(3,1)
    a(24) = p%chi(1,2)
    a(25) = p%chi(2,2)
    a(26) = p%chi(3,2)
    a(27) = p%chi(1,3)
    a(28) = p%chi(2,3)
    a(29) = p%chi(3,3)
    a(30) = p%chi(1,4)
    a(31) = p%chi(2,4)
    a(32) = p%chi(3,4)
    a(33) = p%deltat   
    a(34) = p%pathlength
  end subroutine load_float

  !convert a vector of floats (a) to particle type data (p)
  subroutine load_particle(a,p)
    real(sp),       intent(in ) :: a(par_type_size)
    type(particle), intent(out) :: p
    p%id       = int(a(1)) !!global particle number
    p%pid      = int(a(2)) !!processor id number
    p%elem     = int(a(3)) !!element containing particle
    p%group    = int(a(4)) !!element group id
    p%mark     = int(a(5)) !!=1 out of domain
    p%tbeg     = a(6 )     !!start time in seconds
    p%tend     = a(7 )     !!end time in seconds
    p%x(1)     = a(8 )     !!x particle position
    p%x(2)     = a(9 )     !!y particle position
    p%x(3)     = a(10)     !!z particle position
    p%xn(1)    = a(11)     !!last x particle position
    p%xn(2)    = a(12)     !!last y particle position
    p%xn(3)    = a(13)     !!last z particle position
    p%el       = a(14)     !!surface elevation at particle point
    p%h        = a(15)     !!bathymetry at particle point
    p%u        = a(16)     !!x-velocity at particle location
    p%v        = a(17)     !!y-velocity at particle location
    p%w        = a(18)     !!sigma-velocity at particle location
    p%zloc     = a(19)     !!z value of particle [m]  
    p%s        = a(20)     !!particle associated scalar           
    p%chi(1,1) = a(21)     !!stage contributions
    p%chi(2,1) = a(22) 
    p%chi(3,1) = a(23) 
    p%chi(1,2) = a(24) 
    p%chi(2,2) = a(25) 
    p%chi(3,2) = a(26) 
    p%chi(1,3) = a(27) 
    p%chi(2,3) = a(28) 
    p%chi(3,3) = a(29) 
    p%chi(1,4) = a(30) 
    p%chi(2,4) = a(31) 
    p%chi(3,4) = a(32) 
    p%deltat   = a(33) 
    p%pathlength = a(34) 
  end subroutine load_particle

  !dump particle info to screen
  subroutine pprint(p)
    type(particle), intent(in) :: p
    write(*,*)
    write(*,*)'id:        ',p%id
    write(*,*)'processor: ',p%pid
    write(*,*)'element:   ',p%elem
    write(*,*)'group:     ',p%group
    write(*,*)'mark:      ',p%mark
    write(*,*)'tbeg:      ',p%tbeg
    write(*,*)'tend:      ',p%tend
    write(*,*)'x:         ',p%x 
    write(*,*)'xn:        ',p%xn
    write(*,*)'el:        ',p%el
    write(*,*)'u:         ',p%u
    write(*,*)'v:         ',p%v
    write(*,*)'w:         ',p%w
    write(*,*)'zloc:      ',p%zloc
    write(*,*)'s:         ',p%s
    write(*,*)'chi1       ',p%chi(:,1)
    write(*,*)'chi2       ',p%chi(:,2)
    write(*,*)'chi3       ',p%chi(:,3)
    write(*,*)'chi4       ',p%chi(:,4)
    write(*,*)'deltat     ',p%deltat   
    write(*,*)'pathlength ',p%pathlength
  end subroutine pprint

  !initialize values
  subroutine zero_out(p)
    type(particle), intent(inout) :: p
    p%id     = 0
    p%pid    = 0 
    p%group  = 0
    p%mark   = 0
    p%tbeg   = 0.0
    p%tend   = 0.0
    p%x      = 0.0
    p%xn     = 0.0
    p%el     = 0.0
    p%u      = 0.0
    p%v      = 0.0
    p%w      = 0.0
    p%zloc   = 0.0
    p%s      = 0.0
    p%chi    = 0.0
    p%deltat = 0.0
    p%pathlength = 0.0
    p%found  = .false.
  end subroutine zero_out 


end module particle_class
