MODULE mod_ncdave
!==============================================================================!
!  Save Time Averaged Fields and Dump to NetCDF File                           !
!                                                                              !
!    see: http://www.cgd.ucar.edu/cms/eaton/cf-metadata/ for info              !
!                                                                              !
!    current time dependent variables set up                                   !
!         el:    surface elevation                                             !
!          u:    x-velocity                                                    !
!          v:    y-velocity                                                    !
!         ww:    z-velocity                                                    !
!         kh:    turbulent diffusivity                                         !
!         km:    turbulent viscosity                                           !
!         t1:    temperature                                                   !
!         s1:    salinity                                                      !
!         ua:    vertically-averaged x-velocity                                !
!         va:    vertically-averaged y-velocity                                !
!          d:    depth at nodes                                                !
!         ua:    vertically-averaged x-velocity                                !
!                In spherical coordinate,vertically-averaged lon-velocity      !
!         va:    vertically-averaged y-velocity                                !
!                In spherical coordinate,vertically-averaged lat-velocity      !
!          d:    depth at nodes                                                !
!        uca:    In spherical coordinate,x-velocity                            !
!                (Polar Stereographic projection)                              !
!        vca:    In spherical coordinate,y-velocity                            !
!                (Polar Stereographic projection)                              !
!       uaca:    In spherical coordinate,vertically-averaged x-velocity        !
!                (Polar Stereographic projection)                              !
!       vaca:    In spherical coordinate,vertically-averaged y-velocity        !
!                (Polar Stereographic projection)                              !
!                                                                              !
!    to add additional variables:                                              !
!      1.) add to list above                                                   !
!      2.) add *_vid to variables vid in section "new variable vid"            !
!      3.) go to definition section "new variable definition"                  !
!      4.) add output section "new variable output"                            !
!==============================================================================!

#if defined (NETCDF_IO)

   USE mod_prec
!JQI{
#  if defined (SPHERICAL) && (CARTESIAN_OUTPUT)
   use mod_sp2xy
#  endif
!JQI}
   USE netcdf
   USE mod_ncdio, only : institution,netcdf_timestring
   implicit none
   save

!--Control Variables----------------------------------------------!
   logical,public  :: cdf_out_ave
   integer,private :: n_cnt              !!counts files written
   integer,private :: n_start            !!start iteration for averaging
   integer,private :: n_its              !!number of iterations in averaged period
   integer,private :: n_per              !!number of periods to average
   character(len=120),private :: cdfname !!netcdf file name

   integer,private :: nout_vars          !!number of variables to average and dump
   character(len=80),private,allocatable, dimension(:) :: cdf_vdp

!--NetCDF IDs----------------------------------------------------!

   !--NetCDF File 
   integer,private :: nc_ofid

   !--Dimensions
   integer,private :: node_did,nele_did
   integer,private :: scl_did,siglay_did,siglev_did
   integer,private :: three_did,four_did
   integer,private :: time_did

   !--Grid Variables
   integer,private :: nprocs_vid,partition_vid
   integer,private :: idens_vid
   integer,private :: x_vid,y_vid,lat_vid,lon_vid
!JQI{   
#  if defined (SPHERICAL) && (CARTESIAN_OUTPUT)             
!  Polar Stereographic projection Cartesian coordinate       
   integer :: xca_vid,yca_vid
#  endif   
!JQI}

   integer,private :: nv_vid
   integer,private :: aw0_vid,awx_vid,awy_vid
   integer,private :: a1u_vid,a2u_vid
   integer,private :: siglay_vid,siglev_vid,siglay_shift_vid
  
   !--Flow Variables 
   integer,private :: time_vid
   integer,private :: iint_vid
   integer,private :: u_vid
   integer,private :: v_vid
   integer,private :: ww_vid
   integer,private :: s1_vid
   integer,private :: t1_vid
   integer,private :: el_vid
   integer,private :: h_vid
   integer,private :: km_vid
   integer,private :: kh_vid
   integer,private :: ua_vid
   integer,private :: va_vid
   integer,private :: d_vid
!  JQI{
#  if defined (SPHERICAL) && (CARTESIAN_OUTPUT)
   integer :: uca_vid
   integer :: vca_vid
   integer :: uaca_vid
   integer :: vaca_vid
#  endif
!JQI}
   ! new variable vid  
   ! add *_vid here, e.g. for rho, add rho_vid


   !Averaged Flow Variable Arrays
   real(sp), private, allocatable, dimension(:,:) :: s1_a
   real(sp), private, allocatable, dimension(:,:) :: t1_a
   real(sp), private, allocatable, dimension(:,:) :: u_a
   real(sp), private, allocatable, dimension(:,:) :: v_a
   real(sp), private, allocatable, dimension(:,:) :: ww_a
   real(sp), private, allocatable, dimension(:,:) :: km_a
   real(sp), private, allocatable, dimension(:,:) :: kh_a
   real(sp), private, allocatable, dimension(:)   :: el_a
   real(sp), private, allocatable, dimension(:)   :: ua_a
   real(sp), private, allocatable, dimension(:)   :: va_a
   real(sp), private, allocatable, dimension(:)   :: d_a

   contains !------------------------------------------------------------------!
            ! handle_ncerr_ave    :   deal with netcdf error                   !
            ! avge_fields         :   average fields and call to write file    !
            ! write_netcdf_ave    :   dumped averaged file in netcdf format    !
            ! putvar2             :   collect variable to global and dump      ! 
            ! -----------------------------------------------------------------!

!==============================================================================|
!==============================================================================|

!------------------------------------------------------------------------------|
!  CHECK NETCDF ERROR AND WRITE MESSAGE TO FILE UNIT IPT                       |
!------------------------------------------------------------------------------|
   SUBROUTINE handle_ncerr_ave(status,errmsge,ipt)
   integer, intent(in) :: status,ipt
   character(len=*)    :: errmsge
   if(status /=nf90_noerr)then
     write(ipt,*)trim(errmsge)
     write(ipt,*)trim(nf90_strerror(status))
     call pstop
   end if
   END SUBROUTINE handle_ncerr_ave

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
!  READ IN PARAMETERS CONTROLLING INPUT/OUTPUT FROM RUNTIME PARAM FILE         |
!==============================================================================|
   SUBROUTINE set_ncd_ave  
   use mod_prec
   use all_vars
   use mod_inp
   use netcdf
   implicit none
!--Local Vars----------------------------|
   real(sp)           :: realvec(150)
   integer            :: intvec(150)
   integer            :: iscan
   character(len=120) :: fname
   character(len=80), dimension(100) :: charvec
   integer            :: i
!----------------------------------------|

   fname = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------|
!   cdf_out_ave: netcdf activation flag        
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_OUT_AVE",LVAL = CDF_OUT_AVE)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_out_ave: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  n_start: start iteration for netcdf averaging         
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"BEG_AVGE",ISCAL = N_START)
   if(iscan /= 0)then
     write(ipt,*)'error reading n_start: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  n_its: number of iterations to average over 
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"INT_AVGE",ISCAL = N_ITS)
   if(iscan /= 0)then
     write(ipt,*)'error reading n_its: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  n_per: number of averaging periods (one file / period) 
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"NUM_AVGE",ISCAL = N_PER)
   if(iscan /= 0)then
     write(ipt,*)'error reading n_per: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!     cdf_vdp: list of variables to write to output file
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_VDP_AVE",CVEC = CHARVEC,NSZE = NOUT_VARS)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_vdp_ave: ',iscan
     call pstop
   end if
   if(nout_vars <= 0)then
     write(ipt,*)'incorrect number of netcdf cdf_vdp_ave variables specified'
     write(ipt,*)'in input file',nout_vars
     call pstop
   end if

   allocate(cdf_vdp(nout_vars))
   cdf_vdp(1:nout_vars)= charvec(1:nout_vars)

!------------------------------------------------------------------------------|
!            SCREEN REPORT OF SET VARIABLES                                    !
!------------------------------------------------------------------------------|
   if(msr)then
     write(ipt,*)''
     write(ipt,*)'!     netcdf averaging parameters           '
     if(cdf_out_ave)then
       write(ipt,*)'!  netcdf i/o            :  active'
       write(ipt,*)'!  start iteration       : ',n_start
       write(ipt,*)'!  # iterations/average  : ',n_its   
       write(ipt,*)'!  # averaging periods   : ',n_per   
       write(ipt,*)'!  # variables to write  : ',nout_vars
       do i=1,nout_vars
         write(ipt,999)i,trim(cdf_vdp(i))
       end do
     else
       write(ipt,*)'!  # netcdf averaging    :  not active'
     end if
   end if


   return
   999 format(' !  variable #',i4,'        :',a13)
   END SUBROUTINE set_ncd_ave 
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
!  Average Fields                                                              |
!   itnum:  model iteration number                                             |
!==============================================================================|
   SUBROUTINE avge_fields(itnum) 

   use all_vars
   implicit none
   integer, intent(in) :: itnum
   logical  :: firstit
   logical  :: lastit
   real(sp) :: fac_ave
   integer  :: i

   firstit = .false.
   lastit  = .false.  
   
   fac_ave = 1.0_sp/float(n_its)
 
   !--------------------------------------------------------
   !process logicals based on iteration number
   !--------------------------------------------------------

   !check if we are before or after averaging area
   if(itnum < n_start .or. itnum >=n_start+n_per*n_its)return  

   !first iteration in averaging process?
   if( mod(itnum-n_start,n_its) == 0)then  
     firstit = .true.
   endif

   !last iteration in process
   if( mod(itnum-n_start+1,n_its) == 0)then  
     lastit = .true.
     n_cnt = (itnum+n_its-n_start-1)/n_its
   endif
   if(firstit) write(*,*)'firstit',itnum
   if(lastit) write(*,*)'lastit',itnum,n_cnt
 

   do i=1,nout_vars

     select case(trim(cdf_vdp(i)))
     case("u")  !!=================u=======================================!
       if(.not.allocated(u_a))then
         allocate(u_a(0:MTElem,kb)) ; u_a = 0.0
       endif
       if(firstit) u_a = 0.0
       u_a = u_a + u 
       if(lastit) u_a = u_a*fac_ave
     case("v")  !!===============v=======================================!
       if(.not.allocated(v_a))then
         allocate(v_a(0:MTElem,kb)) ; v_a = 0.0
       endif
       if(firstit) v_a = 0.0
       v_a = v_a + v 
       if(lastit) v_a = v_a*fac_ave
     case("ww") !!===============ww======================================!
       if(.not.allocated(ww_a))then
         allocate(ww_a(0:MTElem,kb)) ; ww_a = 0.0
       endif
       if(firstit) ww_a = 0.0
       ww_a = ww_a + ww
       if(lastit) ww_a = ww_a*fac_ave
     case("km") !!===============km======================================!
       if(.not.allocated(km_a))then
         allocate(km_a(0:MTElem,kb)) ; km_a = 0.0
       endif
       if(firstit) km_a = 0.0
       km_a = km_a + km 
       if(lastit) km_a = km_a*fac_ave
     case("kh") !!===============kh======================================!
       if(.not.allocated(kh_a))then
         allocate(kh_a(0:MTElem,kb)) ; kh_a = 0.0
       endif
       if(firstit) kh_a = 0.0
       kh_a = kh_a + kh
       if(lastit) kh_a = kh_a*fac_ave
     case("t1") !!===============t1======================================!
       if(.not.allocated(t1_a))then
         allocate(t1_a(0:NTNode,kb)) ; t1_a = 0.0
       endif
       if(firstit) t1_a = 0.0
       t1_a = t1_a + t1
       if(lastit) t1_a = t1_a*fac_ave
     case("s1") !!===============s1======================================!
       if(.not.allocated(s1_a))then
         allocate(s1_a(0:NTNode,kb)) ; s1_a = 0.0
       endif
       if(firstit) s1_a = 0.0
       s1_a = s1_a + s1
       if(lastit) s1_a = s1_a*fac_ave
     case("el") !!===============el======================================!
       if(.not.allocated(el_a))then
         allocate(el_a(0:NTNode)) ; el_a = 0.0
       endif
       if(firstit) el_a = 0.0
       el_a = el_a + el 
       if(lastit) el_a = el_a*fac_ave
     case("d") !!===============d=======================================!
       if(.not.allocated(d_a))then
         allocate(d_a(0:NTNode)) ; d_a = 0.0
       endif
       if(firstit) d_a = 0.0
       d_a = d_a + d 
       if(lastit) d_a = d_a*fac_ave
     case("ua") !!===============ua======================================!
       if(.not.allocated(ua_a))then
         allocate(ua_a(0:MTElem)) ; ua_a = 0.0
       endif
       if(firstit) ua_a = 0.0
       ua_a = ua_a + ua
       if(lastit) ua_a = ua_a*fac_ave
     case("va") !!===============va======================================!
       if(.not.allocated(va_a))then
         allocate(va_a(0:MTElem)) ; va_a = 0.0
       endif
       if(firstit) va_a = 0.0
       va_a = va_a + va
       if(lastit) va_a = va_a*fac_ave
!ex  case("var") !!===============var====================================!
     case default 

     end select 

   end do

   !call to write netcdf file
   if(lastit) call write_netcdf_ave


   END SUBROUTINE avge_fields



!==============================================================================|
!  Write NetCDF Averaging File                                                 |
!==============================================================================|
   SUBROUTINE write_netcdf_ave

   use all_vars
#  if defined (MULTIPROCESSOR)
   use mod_par 
#  endif
   use netcdf
   use mod_types
   use mod_utils
   implicit none
   integer, dimension(2) :: stat3de_lev,stat3de_lay 
   integer, dimension(2) :: stat3dn_lev,stat3dn_lay 
   integer, dimension(2) :: specdim
   integer, dimension(1) :: stat2de,stat2dn
   integer, dimension(1) :: stat_lev,stat_lay,stat_scl
   character(len=100)    :: netcdf_convention
   character(len=100)    :: timestamp ,temp
   integer               :: i,j,ierr,i1,i2
   integer               :: maxnode,maxnodep,maxelem,maxelemp
   real(sp), allocatable :: tmp(:,:),tvec(:)
   character(len=4)      :: nchar
   integer :: dims(1)

!==============================================================================|

!==============================================================================|
!  Set up Constants and Initialize Counters                                    |
!==============================================================================|

!--NetCDF Convention String
   netcdf_convention = 'CF-1.0'

!--Time Stamp for History
   call get_timestamp(temp)
   timestamp = 'model started at: '//trim(temp)

!==============================================================================|
!  OPEN FILE AND DEFINE VARIABLES                                              |
!==============================================================================|
   if(msr)then

!--Define NetCDF Output Filename 
   write(nchar,'(I4)')n_cnt
   if(n_cnt < 10)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_ave_000'//trim(adjustl(nchar))//'.nc'
   elseif(n_cnt < 100)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_ave_00'//trim(adjustl(nchar))//'.nc'
   elseif(n_cnt < 1000)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_ave_0'//trim(adjustl(nchar))//'.nc'
   elseif(n_cnt < 10000)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_ave_'//trim(adjustl(nchar))//'.nc'
   else
     write(*,*)'error in netcdf module'
     write(*,*)'# history files > 10000'
     stop
   endif

!--Create File 
   ierr = nf90_create(path=cdfname,cmode=nf90_clobber,ncid=nc_ofid)
   if(ierr /= nf90_eexist)then
     call handle_ncerr_ave(ierr,"file creation error",ipt)
   else
     write(ipt,*)'file :',cdfname,' already exists'
     write(ipt,*)'exiting'
     stop
   end if

!--Description of File Contents
   ierr = nf90_put_att(nc_ofid,nf90_global,"title"      ,trim(casetitle))
   ierr = nf90_put_att(nc_ofid,nf90_global,"institution",trim(institution))
   ierr = nf90_put_att(nc_ofid,nf90_global,"source"     ,trim(fvcom_version))
   ierr = nf90_put_att(nc_ofid,nf90_global,"history"    ,trim(timestamp))
   ierr = nf90_put_att(nc_ofid,nf90_global,"references" ,trim(fvcom_website))
   ierr = nf90_put_att(nc_ofid,nf90_global,"Conventions",trim(netcdf_convention))

!--Define Fixed Model Dimensions 
   ierr = nf90_def_dim(nc_ofid,"scalar" ,1      ,scl_did    )        
   ierr = nf90_def_dim(nc_ofid,"node"   ,NNodeGL    ,node_did   )        
   ierr = nf90_def_dim(nc_ofid,"nele"   ,MElemGL    ,nele_did   )
   ierr = nf90_def_dim(nc_ofid,"siglay" ,kbm1   ,siglay_did )
   ierr = nf90_def_dim(nc_ofid,"siglev" ,kb     ,siglev_did )
   ierr = nf90_def_dim(nc_ofid,"three"  ,3      ,three_did  )
   ierr = nf90_def_dim(nc_ofid,"four"   ,4      ,four_did   )

!--Define Unlimited Model Dimension
   ierr = nf90_def_dim(nc_ofid,"time"   ,nf90_unlimited,time_did)

!--Set Up Data Dimensioning - Static Vars
   stat_scl     = (/scl_did/)             !!scalar variable               
   stat_lay     = (/siglay_did/)          !!vertical variables at layers
   stat_lev     = (/siglev_did/)          !!vertical variables at levels
   stat2de      = (/nele_did/)            !!2d element vars
   stat2dn      = (/node_did/)            !!2d nodal vars
   stat3de_lay  = (/nele_did,siglay_did/) !!3d element vars at layers
   stat3de_lev  = (/nele_did,siglev_did/) !!3d element vars at levels
   stat3dn_lay  = (/node_did,siglay_did/) !!3d node    vars at layers
   stat3dn_lev  = (/node_did,siglev_did/) !!3d node    vars at levels

!--Define Coordinate Variables and Attributes

   !!====NPROCS: Number of Processors=======================!
   ierr = nf90_def_var(nc_ofid,"nprocs",nf90_int,stat_scl,nprocs_vid)
   ierr = nf90_put_att(nc_ofid,nprocs_vid,"long_name","number of processors")

   !!====PARTITION: Partion Number of Element===============!
   ierr = nf90_def_var(nc_ofid,"partition",nf90_int,stat2de,partition_vid)
   ierr = nf90_put_att(nc_ofid,partition_vid,"long_name","partition")

   !!====Initial Density (Used for Constructing 3D Domain)==!
   ierr = nf90_def_var(nc_ofid,"Initial_Density",nf90_float,stat3dn_lay,idens_vid)
   ierr = nf90_put_att(nc_ofid,idens_vid,"long_name","Initial Density")

   !!====X Grid Coordinate at Nodes (VX) (Meters)===========!
   ierr = nf90_def_var(nc_ofid,"x",nf90_float,stat2dn,x_vid)
   ierr = nf90_put_att(nc_ofid,x_vid,"long_name","nodal x-coordinate")
   ierr = nf90_put_att(nc_ofid,x_vid,"units","meters")

   !!====Y Grid Coordinate at Nodes (VY) (Meters)===========!
   ierr = nf90_def_var(nc_ofid,"y",nf90_float,stat2dn,y_vid)
   ierr = nf90_put_att(nc_ofid,y_vid,"long_name","nodal y-coordinate")
   ierr = nf90_put_att(nc_ofid,y_vid,"units","meters")

   !!====Longitudinal Coordinate at Nodes (LON) (degrees)===!
   ierr = nf90_def_var(nc_ofid,"lon",nf90_float,stat2dn,lon_vid)
   ierr = nf90_put_att(nc_ofid,lon_vid,"long_name","Longitude")
   ierr = nf90_put_att(nc_ofid,lon_vid,"standard_name","longitude")
   ierr = nf90_put_att(nc_ofid,lon_vid,"units","degrees_east")

   !!====Latitudinal  Coordinate at Nodes (LAT) (degrees)===!
   ierr = nf90_def_var(nc_ofid,"lat",nf90_float,stat2dn,lat_vid)
   ierr = nf90_put_att(nc_ofid,lat_vid,"long_name","Latitude")
   ierr = nf90_put_att(nc_ofid,lat_vid,"standard_name","latitude")
   ierr = nf90_put_att(nc_ofid,lat_vid,"units","degrees_north")
!JQI{
#  if defined (SPHERICAL) && (CARTESIAN_OUTPUT)
   !!====X Grid Coordinate at Nodes (VXCA) (Meters) in Polar Stereographic projection==!
   ierr = nf90_def_var(nc_ofid,"xca",nf90_float,stat2dn,xca_vid)
   ierr = nf90_put_att(nc_ofid,xca_vid,"long_name","nodal x-coordinate in Polar Stereographic projection")
   ierr = nf90_put_att(nc_ofid,xca_vid,"units","meters")
   
   !!====Y Grid Coordinate at Nodes (VYCA) (Meters) in Polar Stereographic projection==!
   ierr = nf90_def_var(nc_ofid,"yca",nf90_float,stat2dn,yca_vid)
   ierr = nf90_put_att(nc_ofid,yca_vid,"long_name","nodal y-coordinate in Polar Stereographic projection")
   ierr = nf90_put_att(nc_ofid,yca_vid,"units","meters")
#  endif
!JQI}


   !!====Sigma Coordinate for Sigma Layers (ZZ)  (-)========!
   ierr = nf90_def_var(nc_ofid,"siglay",nf90_float,stat_lay,siglay_vid)
   ierr = nf90_put_att(nc_ofid,siglay_vid,"long_name","Sigma Layers")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"formula_terms","siglay:siglay eta:zeta depth:depth")

   !!====Shifted Sigma Layer Coordinate for Viz ============!
   ierr = nf90_def_var(nc_ofid,"siglay_shift",nf90_float,stat_lay,siglay_shift_vid)
   ierr = nf90_put_att(nc_ofid,siglay_shift_vid,"long_name","Shifted Sigma Layers")

   !!====Sigma Coordinate for Sigma Levels (Z)   (-)========!
   ierr = nf90_def_var(nc_ofid,"siglev",nf90_float,stat_lev,siglev_vid)
   ierr = nf90_put_att(nc_ofid,siglev_vid,"long_name","Sigma Levels")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"formula_terms","siglev:siglev eta:zeta depth:depth")



!--Define Mesh Relevant Variables and Attributes

   !!====Bathymetry at Nodes (H) (meters)===================!
   ierr = nf90_def_var(nc_ofid,"h",nf90_float,stat2dn,h_vid)
   ierr = nf90_put_att(nc_ofid,h_vid,"long_name","Bathymetry")   
   ierr = nf90_put_att(nc_ofid,h_vid,"units","meters")
   ierr = nf90_put_att(nc_ofid,h_vid,"positive","down")
   ierr = nf90_put_att(nc_ofid,h_vid,"standard_name","depth")
   ierr = nf90_put_att(nc_ofid,h_vid,"grid","fvcom_grid")

   !!====Nodes surrounding each Element (NV)================!
   specdim = (/nele_did,three_did/) 
   ierr = nf90_def_var(nc_ofid,"nv",nf90_float,specdim,nv_vid)
   ierr = nf90_put_att(nc_ofid,nv_vid,"long_name","nodes surrounding element")     

   !!====Momentum Stencil Interpolation Coefficients========!
   specdim = (/nele_did,four_did/) 
   ierr = nf90_def_var(nc_ofid,"a1u",nf90_float,specdim,a1u_vid)
   ierr = nf90_put_att(nc_ofid,a1u_vid,"long_name","a1u")
   ierr = nf90_def_var(nc_ofid,"a2u",nf90_float,specdim,a2u_vid)
   ierr = nf90_put_att(nc_ofid,a2u_vid,"long_name","a2u")

   !!====Element Based Interpolation Coefficients===========!
   specdim = (/nele_did,three_did/) 
   ierr = nf90_def_var(nc_ofid,"aw0",nf90_float,specdim,aw0_vid)
   ierr = nf90_put_att(nc_ofid,aw0_vid,"long_name","aw0")
   ierr = nf90_def_var(nc_ofid,"awx",nf90_float,specdim,awx_vid)
   ierr = nf90_put_att(nc_ofid,awx_vid,"long_name","awx")
   ierr = nf90_def_var(nc_ofid,"awy",nf90_float,specdim,awy_vid)
   ierr = nf90_put_att(nc_ofid,awy_vid,"long_name","awy")

!--Define Model Time Variables and Attributes    
   ierr = nf90_def_var(nc_ofid,"time",nf90_float,stat_scl,time_vid)
   ierr = nf90_put_att(nc_ofid,time_vid,"long_name","Time")
   ierr = nf90_put_att(nc_ofid,time_vid,"units",trim(netcdf_timestring))
   ierr = nf90_put_att(nc_ofid,time_vid,"calendar","none")
   ierr = nf90_def_var(nc_ofid,"iint",nf90_float,stat_scl,iint_vid)
   ierr = nf90_put_att(nc_ofid,iint_vid,"long_name","internal mode iteration number")

!--Define Time Dependent Flow Variables (selected by user from input file)
   do i=1,nout_vars

     select case(trim(cdf_vdp(i)))

     case("u")  !!===============u=======================================!
     ierr = nf90_def_var(nc_ofid,"u",nf90_float,stat3de_lay,u_vid)
     ierr = nf90_put_att(nc_ofid,u_vid,"long_name","Eastward Water Velocity")
     ierr = nf90_put_att(nc_ofid,u_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,u_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,u_vid,"type","data")
       
     case("v")  !!===============v=======================================!
     ierr = nf90_def_var(nc_ofid,"v",nf90_float,stat3de_lay,v_vid)
     ierr = nf90_put_att(nc_ofid,v_vid,"long_name","Northward Water Velocity")
     ierr = nf90_put_att(nc_ofid,v_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,v_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,v_vid,"type","data")

     case("ww") !!===============ww======================================!
     ierr = nf90_def_var(nc_ofid,"ww",nf90_float,stat3de_lay,ww_vid)
     ierr = nf90_put_att(nc_ofid,ww_vid,"long_name","Upward Water Velocity")
     ierr = nf90_put_att(nc_ofid,ww_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,ww_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,ww_vid,"type","data")

     case("km") !!===============km======================================!
     ierr = nf90_def_var(nc_ofid,"km",nf90_float,stat3de_lay,km_vid)
     ierr = nf90_put_att(nc_ofid,km_vid,"long_name","Turbulent Eddy Viscosity")
     ierr = nf90_put_att(nc_ofid,km_vid,"units","meters2 s-1")
     ierr = nf90_put_att(nc_ofid,km_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,km_vid,"type","data")

     case("kh") !!===============kh======================================!
     ierr = nf90_def_var(nc_ofid,"kh",nf90_float,stat3de_lay,kh_vid)
     ierr = nf90_put_att(nc_ofid,kh_vid,"long_name","Turbulent Eddy Diffusivity")
     ierr = nf90_put_att(nc_ofid,kh_vid,"units","meters2 s-1")
     ierr = nf90_put_att(nc_ofid,kh_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,kh_vid,"type","data")

     case("t1") !!===============t1======================================!
     ierr = nf90_def_var(nc_ofid,"temp",nf90_float,stat3dn_lay,t1_vid)
     ierr = nf90_put_att(nc_ofid,t1_vid,"long_name","temperature")
     ierr = nf90_put_att(nc_ofid,t1_vid,"standard_name","sea_water_temperature")
     ierr = nf90_put_att(nc_ofid,t1_vid,"units","degrees_C")
     ierr = nf90_put_att(nc_ofid,t1_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,t1_vid,"type","data")

     case("s1") !!===============s1======================================!
     ierr = nf90_def_var(nc_ofid,"salinity",nf90_float,stat3dn_lay,s1_vid)
     ierr = nf90_put_att(nc_ofid,s1_vid,"long_name","salinity")
     ierr = nf90_put_att(nc_ofid,s1_vid,"standard_name","sea_water_salinity")
     ierr = nf90_put_att(nc_ofid,s1_vid,"units","1e-3")
     ierr = nf90_put_att(nc_ofid,s1_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,s1_vid,"type","data")

     case("el") !!===============el======================================!
     ierr = nf90_def_var(nc_ofid,"zeta",nf90_float,stat2dn,el_vid)
     ierr = nf90_put_att(nc_ofid,el_vid,"long_name","Water Surface Elevation")
     ierr = nf90_put_att(nc_ofid,el_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,el_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,el_vid,"standard_name","sea_surface_elevation")
     ierr = nf90_put_att(nc_ofid,el_vid,"type","data")

     case("d") !!===============d=======================================!
     ierr = nf90_def_var(nc_ofid,"depth",nf90_float,stat2dn,d_vid)
     ierr = nf90_put_att(nc_ofid,d_vid,"long_name","Water Depth")
     ierr = nf90_put_att(nc_ofid,d_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,d_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,d_vid,"type","data")

     case("ua") !!===============ua======================================!
     ierr = nf90_def_var(nc_ofid,"ua",nf90_float,stat2de,ua_vid)
     ierr = nf90_put_att(nc_ofid,ua_vid,"long_name","Vertically Averaged x-velocity")
     ierr = nf90_put_att(nc_ofid,ua_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,ua_vid,"type","data")

     case("va") !!===============va======================================!
     ierr = nf90_def_var(nc_ofid,"va",nf90_float,stat2de,va_vid)
     ierr = nf90_put_att(nc_ofid,va_vid,"long_name","Vertically Averaged y-velocity")
     ierr = nf90_put_att(nc_ofid,va_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,va_vid,"type","data")

!JQI{
#    if defined (SPHERICAL) && (CARTESIAN_OUTPUT)
!    In Polar Stereographic projection
     case("uca")!!===============uca=====================================!
     ierr = nf90_def_var(nc_ofid,"uca",nf90_float,dynm3de_lay,uca_vid)
     ierr = nf90_put_att(nc_ofid,uca_vid,"long_name","X Velocity in Polar Stereographic projection")
     ierr = nf90_put_att(nc_ofid,uca_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,uca_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,uca_vid,"type","data")

     case("vca")!!===============vca=====================================!
     ierr = nf90_def_var(nc_ofid,"vca",nf90_float,dynm3de_lay,vca_vid)
     ierr = nf90_put_att(nc_ofid,vca_vid,"long_name","Y Velocity in Polar Stereographic projection")
     ierr = nf90_put_att(nc_ofid,vca_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,vca_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,vca_vid,"type","data")

     case("uaca")!!===============uaca===================================!
     ierr = nf90_def_var(nc_ofid,"uaca",nf90_float,dynm2de,uaca_vid)
     ierr = nf90_put_att(nc_ofid,uaca_vid,"long_name","Vertically Averaged x-velocity in PS projection")
     ierr = nf90_put_att(nc_ofid,uaca_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,uaca_vid,"type","data")

     case("vaca")!!===============vaca======================================!
     ierr = nf90_def_var(nc_ofid,"vaca",nf90_float,dynm2de,vaca_vid)
     ierr = nf90_put_att(nc_ofid,vaca_vid,"long_name","Vertically Averaged y-velocity in PS projection")
     ierr = nf90_put_att(nc_ofid,vaca_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,vaca_vid,"type","data")
#    endif
!JQI}



!ex  case("var") !!===============var====================================!
!ex  ierr = nf90_def_var(nc_ofid,"truevar",nf90_float,dimensions,var_vid)
!ex  ierr = nf90_put_att(nc_ofid,var_vid,"long_name","A Good Descriptive Name")
!ex  ierr = nf90_put_att(nc_ofid,var_vid,"units","UDUNITS compatible units")
!ex  ierr = nf90_put_att(nc_ofid,var_vid,"standard_name","CF-convention standard name")
!ex  ierr = nf90_put_att(nc_ofid,var_vid,"type","data")

     !    new variable definition
     !1.) add new definition above here by copying example above and modifying
     !2.) copy dimensions from variable which has same dimensions as var 
     !3.) change variable name if necessary to something more descriptive
     !   e.g. model name for temperature is t1, use temp instead 
     !4.) give the variable a reasonable "long_name"
     !5.) look up the variables standard_name from the cf-convention standard_name list
     !   http://www.cgd.ucar.edu/cms/eaton/cf-metadata/standard_name.html
     !   if it does not exist, do not provide a standard name attribute
     !6.) set variable units conforming to udunits standard   
     !   http://my.unidata.ucar.edu/content/software/udunits/index.html
      


     case default 
       write(ipt,*)'variable',cdf_vdp(i),' not set up for netcdf output'
       write(ipt,*)'modify module mod_ncdio.f' 
       call pstop
     end select

   end do

!--Exit Define Mode
   ierr = nf90_enddef(nc_ofid)
   ierr = nf90_close(nc_ofid)

   end if !(msr)

!==============================================================================|
!  WRITE VARIABLES TO FILE                                                     |
!==============================================================================|
   if(msr)then
     ierr = nf90_open(cdfname,nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr_ave(ierr,"file open error",ipt)
     end if
   end if
   
   !!====Longitude at Nodes (LON) ==========================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar2(i1,i2,NNode,NNodeGL,1,1,"n",vx+vxmin,nc_ofid,lon_vid,myid,nprocs,ipt)

   !!====Latitude  at Nodes (LAT) ==========================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar2(i1,i2,NNode,NNodeGL,1,1,"n",vx+vxmin,nc_ofid,lat_vid,myid,nprocs,ipt)

   !!====Number of Processors (NPROCS) =====================!
   if(msr)then 
   ierr = nf90_put_var(nc_ofid,nprocs_vid,nprocs)
   if(ierr /= nf90_noerr)then
     call handle_ncerr_ave(ierr,"error writing nprocs variable to netcdf",ipt)
   end if
#  if defined (MULTIPROCESSOR)
   ierr = nf90_put_var(nc_ofid,partition_vid,el_pid)
   if(ierr /= nf90_noerr)then
     call handle_ncerr_ave(ierr,"error writing el_pid variable to netcdf",ipt)
   end if
#  endif
   end if
!JQI{
#  if defined (SPHERICAL) && (CARTESIAN_OUTPUT)
!  Polar Stereographic projection 
   !!====X Grid Coordinate at Nodes in Polar Stereographic projection(VXCA)====================!
   i1 = lbound(vxca,1) ; i2 = ubound(vxca,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vxca+vxmin,nc_ofid,xca_vid,myid,nprocs,ipt)

   !!====Y Grid Coordinate at Nodes in Polar Stereographic projection(VYCA)====================!
   i1 = lbound(vyca,1) ; i2 = ubound(vyca,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vyca+vymin,nc_ofid,yca_vid,myid,nprocs,ipt)
#  endif
!JQI}


   !!====Initial Density Field==============================!
   i1 = lbound(rho1,1) ; i2 = ubound(rho1,1)
   call putvar2(i1,i2,NNode,NNodeGL,kb,kb-1,"n",rho1,nc_ofid,idens_vid,myid,nprocs,ipt)


   !!====X Grid Coordinate at Nodes (VX)====================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar2(i1,i2,NNode,NNodeGL,1,1,"n",vx+vxmin,nc_ofid,x_vid,myid,nprocs,ipt)

   !!====Y Grid Coordinate at Nodes (VY)====================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar2(i1,i2,NNode,NNodeGL,1,1,"n",vy+vymin,nc_ofid,y_vid,myid,nprocs,ipt)

   !!====Bathymetry at Nodes (H)============================!
   i1 = lbound(h,1) ; i2 = ubound(h,1)
   call putvar2(i1,i2,NNode,NNodeGL,1,1,"n",h,nc_ofid,h_vid,myid,nprocs,ipt)

   !!====Nodes surrounding each Element (NV)================!
   allocate(tmp(0:MTElem,3))
   if(serial)then
     tmp(0:MTElem,1:3) = real(nv(0:MTElem,1:3),sp) 
   end if
#  if defined (MULTIPROCESSOR)
   if(par)then
   do j=1,3
   do i=1,MElem
     tmp(i,j) = real(ngid(nv(i,j)),sp)
   end do
   end do
   end if
#  endif
   i1 = lbound(tmp,1) ; i2 = ubound(tmp,1)
   call putvar2(i1,i2,MElem,MElemGL,3,3,"e",tmp,nc_ofid,nv_vid,myid,nprocs,ipt)
   deallocate(tmp)

   !!====Momentum Stencil Interpolation Coefficients========!
   i1 = lbound(a1u,1) ; i2 = ubound(a1u,1)
   call putvar2(i1,i2,MElem,MElemGL,4,4,"e",a1u,nc_ofid,a1u_vid,myid,nprocs,ipt)
   i1 = lbound(a2u,1) ; i2 = ubound(a2u,1)
   call putvar2(i1,i2,MElem,MElemGL,4,4,"e",a2u,nc_ofid,a2u_vid,myid,nprocs,ipt)

   !!====Element Based Interpolation Coefficients===========!
   i1 = lbound(aw0,1) ; i2 = ubound(aw0,1)
   call putvar2(i1,i2,MElem,MElemGL,3,3,"e",aw0,nc_ofid,aw0_vid,myid,nprocs,ipt)
   i1 = lbound(awx,1) ; i2 = ubound(awx,1)
   call putvar2(i1,i2,MElem,MElemGL,3,3,"e",awx,nc_ofid,awx_vid,myid,nprocs,ipt)
   i1 = lbound(awy,1) ; i2 = ubound(awy,1)
   call putvar2(i1,i2,MElem,MElemGL,3,3,"e",awy,nc_ofid,awy_vid,myid,nprocs,ipt)

   !!====Sigma Layers (zz)==================================!
   if(msr)then 
   allocate(tvec(kbm1))
!   tvec(1:kbm1) = zz(1:kbm1)
   tvec(1:kbm1) = zz(1,1:kbm1)  !T.W. added
   ierr = nf90_put_var(nc_ofid,siglay_vid,tvec)
   if(ierr /= nf90_noerr)then
     call handle_ncerr_ave(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(tvec)

   allocate(tvec(kbm1))
!   tvec(1:kbm1) = z(2:kb)
   tvec(1:kbm1) = z(1,2:kb)  !T.W. added
   ierr = nf90_put_var(nc_ofid,siglay_shift_vid,tvec)
   if(ierr /= nf90_noerr)then
     call handle_ncerr_ave(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(tvec)

   allocate(tvec(kb))
!   tvec(1:kb) = z(1:kb)
   tvec(1:kb) = z(1,1:kb)  !T.W. added
   ierr = nf90_put_var(nc_ofid,siglev_vid,tvec)
   if(ierr /= nf90_noerr)then
     call handle_ncerr_ave(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(tvec)
   endif
   
!--Dump Time/IINT to File
   dims(1) = 1 
   ierr    = nf90_put_var(nc_ofid,iint_vid,float(iint),START=dims)
   ierr    = nf90_put_var(nc_ofid,time_vid,thour*3600.,START=dims)


!--Write Variables to File
   if(msr) write(ipt,*)'dumping averaged netcdf file: ',trim(cdfname)

   do i=1,nout_vars

     select case(trim(cdf_vdp(i)))

     case("u")  !!===============U=======================================!
       i1 = lbound(u,1) ; i2 = ubound(u,1)
       call putvar2(i1,i2,MElem,MElemGL,kb,kb-1,"e",u_a,nc_ofid,u_vid,myid,nprocs,ipt)

     case("v")  !!===============V=======================================!
       i1 = lbound(v,1) ; i2 = ubound(v,1)
       call putvar2(i1,i2,MElem,MElemGL,kb,kb-1,"e",v_a,nc_ofid,v_vid,myid,nprocs,ipt)

     case("ww") !!===============WW======================================!
       i1 = lbound(ww,1) ; i2 = ubound(ww,1)
       call putvar2(i1,i2,MElem,MElemGL,kb,kb-1,"e",ww_a,nc_ofid,ww_vid,myid,nprocs,ipt)

     case("km") !!===============KM======================================!
       i1 = lbound(km,1) ; i2 = ubound(km,1)
       call putvar2(i1,i2,NNode,NNodeGL,kb,kb-1,"n",km_a,nc_ofid,km_vid,myid,nprocs,ipt)

     case("kh") !!===============KH======================================!
       i1 = lbound(kh,1) ; i2 = ubound(kh,1)
       call putvar2(i1,i2,NNode,NNodeGL,kb,kb-1,"n",kh_a,nc_ofid,kh_vid,myid,nprocs,ipt)

     case("el") !!===============EL======================================!
       i1 = lbound(el,1) ; i2 = ubound(el,1)
       call putvar2(i1,i2,NNode,NNodeGL,1,1,"n",el_a,nc_ofid,el_vid,myid,nprocs,ipt)

     case("d") !!===============D=======================================!
       i1 = lbound(d,1) ; i2 = ubound(d,1)
       call putvar2(i1,i2,NNode,NNodeGL,1,1,"n",d_a,nc_ofid,d_vid,myid,nprocs,ipt)

     case("t1") !!===============T1======================================!
       i1 = lbound(t1,1) ; i2 = ubound(t1,1)
       call putvar2(i1,i2,NNode,NNodeGL,kb,kb-1,"n",t1_a,nc_ofid,t1_vid,myid,nprocs,ipt)

     case("s1") !!===============S1======================================!
       i1 = lbound(s1,1) ; i2 = ubound(s1,1)
       call putvar2(i1,i2,NNode,NNodeGL,kb,kb-1,"n",s1_a,nc_ofid,s1_vid,myid,nprocs,ipt)

     case("ua") !!===============UA======================================!
       i1 = lbound(ua,1) ; i2 = ubound(ua,1)
       call putvar2(i1,i2,MElem,MElemGL,1,1,"e",ua_a,nc_ofid,ua_vid,myid,nprocs,ipt)

     case("va") !!===============VA======================================!
       i1 = lbound(va,1) ; i2 = ubound(va,1)
       call putvar2(i1,i2,MElem,MElemGL,1,1,"e",va_a,nc_ofid,va_vid,myid,nprocs,ipt)

!JQI{
#    if defined (SPHERICAL) && (CARTESIAN_OUTPUT)
     case("uca")!!===============UCA=======================================!
       i1 = lbound(uca,1) ; i2 = ubound(uca,1)
       call putvar(i1,i2,MElem,MElemGL,kb,kb-1,"e",uca,nc_ofid,uca_vid,myid,nprocs,ipt)
     
     case("vca")!!===============VCA=======================================!
       i1 = lbound(vca,1) ; i2 = ubound(vca,1)
       call putvar(i1,i2,MElem,MElemGL,kb,kb-1,"e",vca,nc_ofid,vca_vid,myid,nprocs,ipt)

     
     case("uaca")!!===============UACA======================================!
       i1 = lbound(uaca,1) ; i2 = ubound(uaca,1)
       call putvar(i1,i2,MElem,MElemGL,1,1,"e",uaca,nc_ofid,uaca_vid,myid,nprocs,ipt)
     
     case("vaca")!!===============VACA======================================!
       i1 = lbound(vaca,1) ; i2 = ubound(vaca,1)
       call putvar(i1,i2,MElem,MElemGL,1,1,"e",vaca,nc_ofid,vaca_vid,myid,nprocs,ipt)

#    endif
!JQI}     


!ex     case("s1") !!===============S1======================================!
!ex       i1 = lbound(s1,1) ; i2 = ubound(s1,1)
!ex       call putvar2(i1,i2,NNode,NNodeGL,kb,kb-1,"n",s1,nc_ofid,s1_vid,myid,nprocs,ipt)

     !new variable output - add a new variable (e.g. 'var') to output
     !1.) copy example section above
     !2.) modify case for your variable 'case("var")' 
     !3.) modify bounds for your variable
     !4.) modify putvar for your variable by finding a putvar for a variable
     !    with same dimensions and type ("e" or "n")
     !5.) modify variable vid with your variables vid (e.g. "var_vid")

     case default
       if(msr)then
         write(ipt,*)'variable',cdf_vdp(i),' not set up for netcdf output'
         write(ipt,*)'modify module MOD_NCDIO.f' 
         call pstop
       end if
     end select

   end do

!==============================================================================|
!  CLOSE THE FILE                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)

   return
   end subroutine write_netcdf_ave


!==============================================================================|
!  Collect Data to Global Array and Write to Netcdf File                       |
!==============================================================================|
                                                                                                                  
   subroutine putvar2(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid,myid,nprocs,ipt)
 
!------------------------------------------------------------------------------|

#  if defined (MULTIPROCESSOR)
   use mod_par
#  endif
   use mod_types
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs,ipt
   character(len=*),intent(in)   :: map_type
   real(sp), dimension(i1:i2,kt) :: var

   real(sp), allocatable, dimension(:,:) :: temp,gtemp
   integer :: ierr,k1m1
   integer, allocatable :: dims(:)
   

   k1m1 = k1 
   if(k1m1 == 1)then
     allocate(dims(2))
     dims(1) = 1 
     dims(2) = 1 
   else
     allocate(dims(3))
     dims(1) = 1 
     dims(2) = 1 
     dims(3) = 1 
   end if
     

   if(map_type(1:1) /= "e" .and. map_type(1:1) /= "n")then
     write(ipt,*)'map_type input to putvar should be "e" OR "n"'
     call pstop
   end if

   if(nprocs==1)then
     allocate(temp(n1,k1m1))  ; temp(1:n1,1:k1m1) = var(1:n1,1:k1m1)
   end if

#  if defined (MULTIPROCESSOR)
   if(nprocs > 1)then
     allocate(gtemp(n1gl,kt))
     if(map_type(1:1) == "e")then
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,emap,var,gtemp)
     else 
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,nmap,var,gtemp)
     end if
     allocate(temp(n1gl,k1m1))  ; temp(1:n1gl,1:k1m1) = gtemp(1:n1gl,1:k1m1)
     deallocate(gtemp)
   end if
#  endif

   if(myid /= 1) return
   ierr = nf90_put_var(nc_fid,vid,temp,START=dims)
   if(ierr /= nf90_noerr)then
     call handle_ncerr_ave(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(dims)

   return
   end subroutine putvar2
!==============================================================================|

#  endif

   END MODULE mod_ncdave
