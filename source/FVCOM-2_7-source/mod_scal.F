!=======================================================================
! FVCOM Scalar Module  
!
!    contains methods:
!        Adv_Scal            => Advect a Scalar Quantity 
!        Vdif_Scal           => Vertical Diffusion of Scalar Quantity
!        Bcond_Scal_OBC      => Open Boundary Condition for Scalar
!        Bcond_Scal_PTsource => Point Sources of Scalar
!=======================================================================
Module Scalar

  logical, parameter :: debug = .true. 

  contains
!==============================================================================|
! Calculate Horizontal Advection and Diffusion For Scalar (f)                  |
!==============================================================================|
  Subroutine Adv_Scal(f,fn,d_fdis,fdis,d_fflux,fflux_obc,deltat,source)
!------------------------------------------------------------------------------|

  use all_vars
  use lims, only: NNode,NTNode,MElem,MTElem,kbm1,kb
  use bcs
  use mod_obcs
# if defined (MULTIPROCESSOR)
  use mod_par
# endif
# if defined (WET_DRY)
  use mod_wd
# endif

  implicit none
  real(sp), intent(in ), dimension(0:NTNode,kb)      :: f 
  real(sp), intent(out), dimension(0:NTNode,kb)      :: fn
  integer , intent(in )                          :: d_fdis
  real(sp), intent(in ), dimension(d_fdis)       :: fdis
  integer , intent(in )                          :: d_fflux
  real(sp), intent(out), dimension(d_fflux,kbm1) :: fflux_obc 
  real(sp), intent(in )                          :: deltat
  logical , intent(in )                          :: source

  !----------------local--------------------------------------
  real(sp), dimension(0:NTNode,kb)   :: xflux,xflux_adv
  real(sp), dimension(NNode)         :: pupx,pupy,pvpx,pvpy  
  real(sp), dimension(NNode)         :: pfpx,pfpy,pfpxd,pfpyd,viscoff
  real(sp), dimension(3*MTElem)      :: dtij 
  real(sp), dimension(3*MTElem,kbm1) :: uvn
  real(sp), dimension(kb)        :: vflux
  real(sp) :: utmp,vtmp,sitai,ffd,ff1,x11,y11,x22,y22,x33,y33
  real(sp) :: tmp1,tmp2,xi,yi
  real(sp) :: dxa,dya,dxb,dyb,fij1,fij2,un
  real(sp) :: txx,tyy,fxx,fyy,viscof,exflux,temp,fpoint
  real(sp) :: fact,fm1,fmean
  integer  :: i,i1,i2,ia,ib,j,j1,j2,k,jtmp,jj
# if defined (SPHERICAL)
  real(sp) :: ty,txpi,typi
# endif

!------------------------------------------------------------------------------!

!-------------------------------------------------------
!Calculate Mean Values
!-------------------------------------------------------

  fmean = sum(f(1:NNode,1:kbm1))/float(NNode*kbm1)

!-------------------------------------------------------
!Initialize Multipliers to Control Horizontal Diff
!-------------------------------------------------------

  fact = 0.0_sp
  fm1  = 1.0_sp
  if(horzmix == 'closure') then
    fact = 1.0_sp
    fm1  = 0.0_sp
  end if
     
!-------------------------------------------------------
!Initialize Fluxes
!-------------------------------------------------------
  xflux     = 0.0_sp
  xflux_adv = 0.0_sp

!-------------------------------------------------------
!Calculate Normal Velocity on Control Volume Edges
!-------------------------------------------------------
!!# if !defined (WET_DRY)
  do i=1,ncv
    i1=ntrg(i)
    dtij(i)=dt1(i1)
    do k=1,kbm1
      uvn(i,k) = v(i1,k)*dltxe(i) - u(i1,k)*dltye(i)
    end do
  end do
!!# else
!!  do i=1,ncv
!!    i1=ntrg(i)
!!    dtij(i)=dt1(i1)
!!    do k=1,kbm1
!!      uvn(i,k) = vs(i1,k)*dltxe(i) - us(i1,k)*dltye(i)
!!    end do
!!  end do
!!# endif

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   do k=1,kbm1
      pfpx  = 0.0_sp 
      pfpy  = 0.0_sp 
      pfpxd = 0.0_sp 
      pfpyd = 0.0_sp
     do i=1,NNode
       do j=1,ntsn(i)-1
         i1=nbsn(i,j)
         i2=nbsn(i,j+1)

#    if defined (WET_DRY)
         IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 1)THEN
          FFD=0.5_SP*(f(I,K)+f(I2,K))
          FF1=0.5_SP*(f(I,K)+f(I2,K))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 1 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(f(I1,K)+f(I,K))
          FF1=0.5_SP*(f(I1,K)+f(I,K))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(f(I,K)+f(I,K))
          FF1=0.5_SP*(f(I,K)+f(I,K))
	 ELSE
          FFD=0.5_SP*(f(I1,K)+f(I2,K))
          FF1=0.5_SP*(f(I1,K)+f(I2,K))
	 END IF 
#    else	 
         ffd=0.5_sp*(f(i1,k)+f(i2,k)) !-fmean1(i1,k)-fmean1(i2,k))
         ff1=0.5_sp*(f(i1,k)+f(i2,k))
#    endif	 
	 
#        if defined (SPHERICAL)
         ty=0.5_sp*(vy(i1)+vy(i2))
         txpi=(vx(i2)-vx(i1))*tpi*cos(deg2rad*ty)
         typi=(vy(i1)-vy(i2))*tpi
         pfpx(i)=pfpx(i)+ff1*typi
         pfpy(i)=pfpy(i)+ff1*txpi
         pfpxd(i)=pfpxd(i)+ffd*typi
         pfpyd(i)=pfpyd(i)+ffd*txpi
#        else
         pfpx(i) = pfpx(i) +ff1*(vy(i1)-vy(i2))
         pfpy(i) = pfpy(i) +ff1*(vx(i2)-vx(i1))
         pfpxd(i)= pfpxd(i)+ffd*(vy(i1)-vy(i2))
         pfpyd(i)= pfpyd(i)+ffd*(vx(i2)-vx(i1))
#        endif
       end do
       pfpx(i)  =pfpx(i )/art2(i)
       pfpy(i)  =pfpy(i )/art2(i)
       pfpxd(i) =pfpxd(i)/art2(i)
       pfpyd(i) =pfpyd(i)/art2(i)
     end do
          
     if(k == kbm1)then
       do i=1,NNode
         pfpxb(i) = pfpx(i)
         pfpyb(i) = pfpy(i)
       end do
     end if

     do i=1,NNode
       pupx(i)=0.0_sp
       pupy(i)=0.0_sp
       pvpx(i)=0.0_sp
       pvpy(i)=0.0_sp
       j=1
       i1=nbve(i,j)
       jtmp=nbvt(i,j)
       j1=jtmp+1-(jtmp+1)/4*3
       j2=jtmp+2-(jtmp+2)/4*3
       x11=0.5_sp*(vx(i)+vx(nv(i1,j1)))
       y11=0.5_sp*(vy(i)+vy(nv(i1,j1)))
       x22=xc(i1)
       y22=yc(i1)
       x33=0.5_sp*(vx(i)+vx(nv(i1,j2)))
       y33=0.5_sp*(vy(i)+vy(nv(i1,j2)))

#      if defined (SPHERICAL)
       ty  =0.5_sp*(y11+y33)
       txpi=(x33-x11)*tpi*cos(deg2rad*ty)
       typi=(y11-y33)*tpi
       pupx(i)=pupx(i)+u(i1,k)*typi 
       pupy(i)=pupy(i)+u(i1,k)*txpi
       pvpx(i)=pvpx(i)+v(i1,k)*typi
       pvpy(i)=pvpy(i)+v(i1,k)*txpi
#      else
       pupx(i)=pupx(i)+u(i1,k)*(y11-y33)
       pupy(i)=pupy(i)+u(i1,k)*(x33-x11)
       pvpx(i)=pvpx(i)+v(i1,k)*(y11-y33)
       pvpy(i)=pvpy(i)+v(i1,k)*(x33-x11)
#      endif

       if(isonb(i) /= 0) then
#        if defined (SPHERICAL)
         ty=0.5_sp*(vy(i)+y11)
         txpi=(x11-vx(i))*tpi*cos(deg2rad*ty)
         typi=(vy(i)-y11)*tpi
         pupx(i)=pupx(i)+u(i1,k)*typi
         pupy(i)=pupy(i)+u(i1,k)*txpi
         pvpx(i)=pvpx(i)+v(i1,k)*typi
         pvpy(i)=pvpy(i)+v(i1,k)*txpi
#        else
         pupx(i)=pupx(i)+u(i1,k)*(vy(i)-y11)
         pupy(i)=pupy(i)+u(i1,k)*(x11-vx(i))
         pvpx(i)=pvpx(i)+v(i1,k)*(vy(i)-y11)
         pvpy(i)=pvpy(i)+v(i1,k)*(x11-vx(i))
#        endif
       end if

       do j=2,ntve(i)-1
         i1=nbve(i,j)
         jtmp=nbvt(i,j)
         j1=jtmp+1-(jtmp+1)/4*3
         j2=jtmp+2-(jtmp+2)/4*3
         x11=0.5_sp*(vx(i)+vx(nv(i1,j1)))
         y11=0.5_sp*(vy(i)+vy(nv(i1,j1)))
         x22=xc(i1)
         y22=yc(i1)
         x33=0.5_sp*(vx(i)+vx(nv(i1,j2)))
         y33=0.5_sp*(vy(i)+vy(nv(i1,j2)))

#        if defined (SPHERICAL)
         ty=0.5_sp*(y11+y33)
         txpi=(x33-x11)*tpi*COS(deg2rad*TY)
         typi=(y11-y33)*tpi
         pupx(i)=pupx(i)+u(i1,k)*typi
         pupy(i)=pupy(i)+u(i1,k)*txpi
         pvpx(i)=pvpx(i)+v(i1,k)*typi
         pvpy(i)=pvpy(i)+v(i1,k)*txpi
#        else
         pupx(i)=pupx(i)+u(i1,k)*(y11-y33)
         pupy(i)=pupy(i)+u(i1,k)*(x33-x11)
         pvpx(i)=pvpx(i)+v(i1,k)*(y11-y33)
         pvpy(i)=pvpy(i)+v(i1,k)*(x33-x11)
#        endif
       end do
       j=ntve(i)
       i1=nbve(i,j)
       jtmp=nbvt(i,j)
       j1=jtmp+1-(jtmp+1)/4*3
       j2=jtmp+2-(jtmp+2)/4*3
       x11=0.5_sp*(vx(i)+vx(nv(i1,j1)))
       y11=0.5_sp*(vy(i)+vy(nv(i1,j1)))
       x22=xc(i1)
       y22=yc(i1)
       x33=0.5_sp*(vx(i)+vx(nv(i1,j2)))
       y33=0.5_sp*(vy(i)+vy(nv(i1,j2)))

#      if defined (SPHERICAL)
       ty=0.5*(Y11+Y33)
       txpi=(x33-x11)*tpi*cos(deg2rad*TY)
       typi=(y11-y33)*tpi
       pupx(i)=pupx(i)+u(i1,k)*typi
       pupy(i)=pupy(i)+u(i1,k)*txpi
       pvpx(i)=pvpx(i)+v(i1,k)*typi
       pvpy(i)=pvpy(i)+v(i1,k)*txpi
#      else
       pupx(i)=pupx(i)+u(i1,k)*(y11-y33)
       pupy(i)=pupy(i)+u(i1,k)*(x33-x11)
       pvpx(i)=pvpx(i)+v(i1,k)*(y11-y33)
       pvpy(i)=pvpy(i)+v(i1,k)*(x33-x11)
#      endif

       if(isonb(i) /= 0) then
#      if defined (SPHERICAL)
         ty=0.5*(Y11+VY(I))
         txpi=(VX(I)-X11)*tpi*COS(deg2rad*ty)
         typi=(Y11-VY(I))*tpi
         pupx(i)=pupx(i)+u(i1,k)*typi
         pupy(i)=pupy(i)+u(i1,k)*txpi
         pvpx(i)=pvpx(i)+v(i1,k)*typi
         pvpy(i)=pvpy(i)+v(i1,k)*txpi
#        else
         pupx(i)=pupx(i)+u(i1,k)*(y11-vy(i))
         pupy(i)=pupy(i)+u(i1,k)*(vx(i)-x11)
         pvpx(i)=pvpx(i)+v(i1,k)*(y11-vy(i))
         pvpy(i)=pvpy(i)+v(i1,k)*(vx(i)-x11)
#        endif
       end if
       pupx(i)=pupx(i)/art1(i)
       pupy(i)=pupy(i)/art1(i)
       pvpx(i)=pvpx(i)/art1(i)
       pvpy(i)=pvpy(i)/art1(i)
       tmp1=pupx(i)**2+pvpy(i)**2
       tmp2=0.5_sp*(pupy(i)+pvpx(i))**2
       viscoff(i)=sqrt(tmp1+tmp2)*art1(i)
     end do
!     if(k == kbm1) then
!       ah_bottom(1:NNode) = horcon*(fact*viscoff(1:NNode) + fm1)
!     end if


     do i=1,ncv_i
       ia=niec(i,1)
       ib=niec(i,2)
       xi=0.5_sp*(xije(i,1)+xije(i,2))
       yi=0.5_sp*(yije(i,1)+yije(i,2))
#      if defined (SPHERICAL)
       ty=0.5_sp*(yi+vy(ia))
       dxa=(xi-vx(ia))*tpi*cos(deg2rad*ty)
       dya=(yi-vy(ia))*tpi
       ty=0.5*(YI+VY(IB))
       DXB=(XI-VX(IB))*tpi*COS(deg2rad*ty)
       DYB=(YI-VY(IB))*tpi
#      else
       dxa=xi-vx(ia)
       dya=yi-vy(ia)
       dxb=xi-vx(ib)
       dyb=yi-vy(ib)
#      endif
       fij1=f(ia,k)+dxa*pfpx(ia)+dya*pfpy(ia)
       fij2=f(ib,k)+dxb*pfpx(ib)+dyb*pfpy(ib)
       un=uvn(i,k)

       viscof=horcon*(fact*(viscoff(ia)+viscoff(ib))*0.5_sp + fm1)

       txx=0.5_sp*(pfpxd(ia)+pfpxd(ib))*viscof
       tyy=0.5_sp*(pfpyd(ia)+pfpyd(ib))*viscof

       fxx=-dtij(i)*txx*dltye(i)
       fyy= dtij(i)*tyy*dltxe(i)

       exflux=-un*dtij(i)* &
          ((1.0_sp+sign(1.0_sp,un))*fij2+(1.0_sp-sign(1.0_sp,un))*fij1)*0.5_sp+fxx+fyy

       xflux(ia,k)=xflux(ia,k)+exflux
       xflux(ib,k)=xflux(ib,k)-exflux

       xflux_adv(ia,k)=xflux_adv(ia,k)+(exflux-fxx-fyy)
       xflux_adv(ib,k)=xflux_adv(ib,k)-(exflux-fxx-fyy)
     end do
  end do !!sigma loop

!---------------------------------------------------------------------------------
! Accumulate Fluxes at Boundary Nodes
!---------------------------------------------------------------------------------
 
# if defined (MULTIPROCESSOR)
  if(par)call node_match(0,nbn,bn_mlt,bn_loc,bnc,NTNode,kb,myid,nprocs,xflux,xflux_adv)
# endif

!---------------------------------------------------------------------------------
! Store Advective Fluxes at the Boundary
!---------------------------------------------------------------------------------
  do k=1,kbm1
     if(iobcn > 0) then
       do i=1,iobcn
         i1=i_obc_n(i)
         fflux_obc(i,k)=xflux_adv(i1,k)
       end do
     end if
  end do

!---------------------------------------------------------------------------------
! Calculate Vertical Advection Terms 
!---------------------------------------------------------------------------------

   do i=1,NNode 
#    if defined (WET_DRY)
     if(ISWET_NODE_CURRENTSTEP(i)*ISWET_NODE_LASTSTEP(i) == 1) then
#    endif
     call calc_vflux(kbm1,f(i,1:kbm1),wts(i,1:kb),vflux)
     do k=1,kbm1
       if(isonb(i) == 2) then
         xflux(i,k)= (vflux(k)-vflux(k+1))*art1(i)/dz(i,k)
       else
         xflux(i,k)=xflux(i,k)+ (vflux(k)-vflux(k+1))*art1(i)/dz(i,k)
       end if
     end do
#    if defined (WET_DRY)
     end if
#    endif
   end do

!-------------------------------------------------------
!Point Source                                      
!-------------------------------------------------------
  if(source)then  !!user specified

  if(point_st_type == 'calculated') then
    if(inflow_type == 'node') then
        do j=1,numqbc
          jj=inodeq(j)
          fpoint=fdis(j)
          do k=1,kbm1
            xflux(jj,k)=xflux(jj,k) - qdis(j)*vqdist(j,k)*fpoint/dz(jj,k)
          end do
        end do
    else if(inflow_type == 'edge') then
      write(*,*)'scalar advection not setup for "edge" point source'
      stop
    end if
  end if

  else

  if(point_st_type == 'calculated') then
    if(inflow_type == 'node') then
        do j=1,numqbc
          jj=inodeq(j)
          do k=1,kbm1
            fpoint = f(jj,k)
            xflux(jj,k)=xflux(jj,k) - qdis(j)*vqdist(j,k)*fpoint/dz(jj,k)
          end do
        end do
    else if(inflow_type == 'edge') then
      write(*,*)'scalar advection not setup for "edge" point source'
      stop
    end if
  end if

  endif
!------------------------------------------------------------------------
!Update Scalar Quantity
!------------------------------------------------------------------------

  do i=1,NNode
#   if defined (WET_DRY)
    if(ISWET_NODE_CURRENTSTEP(i)*ISWET_NODE_LASTSTEP(i) == 1 )then
#   endif
    do k=1,kbm1
      fn(i,k)=(f(i,k)-xflux(i,k)/art1(i)*(deltat/dt(i)))*(dt(i)/dtfa(i))
    end do
#   if defined (WET_DRY)
    else
    do k=1,kbm1
      fn(i,k)=f(i,k)
    end do
    end if
#   endif
  end do

  return
  End Subroutine Adv_Scal
!==============================================================================|

!==============================================================================|
! Vertical Diffusion of Scalar                                                 |
!==============================================================================|
  Subroutine Vdif_Scal(f,deltat)

  use mTridiagonal_scal
  use all_vars 
  Implicit None 
  Real(sp), intent(inout) :: f(0:NTNode,kb)
  Real(sp), intent(in   ) :: deltat
  !--local--------------------
  integer  :: i,k,ll
  real(sp) :: dsqrd,dfdz,visb
  real(sp) :: fsol(0:kb)

  call init_tridiagonal_scal(kb)

  Do i=1,NNode
     dsqrd = d(i)*d(i)

    !----------------------------------------------------------------
    !  Set up Diagonals of Matrix (lower=au_m,diag=bu_m,upper=cu_m)
    !----------------------------------------------------------------
    

    !Surface
    au_m(1) = 0.0
    cu_m(1)=      - deltat*(kh(i,2)+umol)/(dzz(i,1)*dz(i,1)*dsqrd)
    bu_m(1)=  1.0 - cu_m(1) 

    !Interior
    do k=2,kbm1-1
      au_m(k) =     - deltat*(kh(i,k  )+umol)/(dzz(i,k-1)*dz(i,k)*dsqrd)
      cu_m(k) =     - deltat*(kh(i,k+1)+umol)/(dzz(i,k  )*dz(i,k)*dsqrd)
      bu_m(k) = 1.0 - cu_m(k) - au_m(k) 
    end do

    !Bottom
     au_m(kbm1) =     - deltat*(kh(i,kbm1)+umol)/(dzz(i,kbm1-1)*dz(i,kbm1)*dsqrd)
     cu_m(kbm1) = 0.0
     bu_m(kbm1) = 1.0 - au_m(kbm1) 

    !----------------------------------------------------------------
    ! Set up RHS forcing vector and boundary conditions 
    !----------------------------------------------------------------
    do k=1,kbm1
      du_m(k) = f(i,k)
    end do

    !Free Surface: No flux

    !Bottom: No flux
      

    !----------------------------------------------------------------
    ! Solve 
    !----------------------------------------------------------------

     call tridiagonal_scal(kb,1,kbm1,fsol)
    
     !Transfer
     f(i,1:kbm1) = fsol(1:kbm1)

  End Do

  End Subroutine Vdif_Scal


!==============================================================================|
! Set Point Source Conditions for Scalar Function                              |
!==============================================================================|

  Subroutine Bcond_Scal_PTsource(f,fn,fdis)

!------------------------------------------------------------------------------|
  use all_vars
  use bcs
  use mod_obcs
  implicit none
  real(sp), intent(in ), dimension(0:NTNode,kb)      :: f 
  real(sp), intent(out), dimension(0:NTNode,kb)      :: fn
  real(sp), intent(in ), dimension(numqbc )      :: fdis
!--local-------------------------------------------
  integer  :: i,j,k,j1,j11,j22
!------------------------------------------------------------------------------|


!--------------------------------------------
! Set Source Terms
!--------------------------------------------
  if(point_st_type == 'specified') then
    if(numqbc > 0) then
      if(inflow_type == 'node') then
        do i=1,numqbc
          j11=inodeq(i)
          do k=1,kbm1
            fn(j11,k)=fdis(i)
          end do
        end do
      else if(inflow_type == 'edge') then
        do i=1,numqbc
          j11=n_icellq(i,1)
          j22=n_icellq(i,2)
          do k=1,kbm1
            fn(j11,k)=fdis(i)
            fn(j22,k)=fdis(i)
          end do
        end do
      end if
    end if
  end if

  return
  End Subroutine Bcond_Scal_PTSource 
!==============================================================================|
!==============================================================================|

!==============================================================================|
!   Set Boundary Conditions for Scalar Function on Open Boundary               |
!==============================================================================|

  Subroutine Bcond_Scal_OBC(f,fn,fflux_obc,f_obc,deltat,alpha_nudge)

!------------------------------------------------------------------------------|
  use all_vars
  use bcs
  use mod_obcs
  implicit none
  real(sp), intent(in   ), dimension(0:NTNode,kb)      :: f 
  real(sp), intent(inout), dimension(0:NTNode,kb)      :: fn
  real(sp), intent(in   ), dimension(iobcn+1,kbm1) :: fflux_obc 
  real(sp), intent(in   ), dimension(iobcn       ) :: f_obc 
  real(sp), intent(in   )                          :: deltat
  real(sp), intent(in   )                          :: alpha_nudge 
!--local-------------------------------------------
  real(sp) :: f2d,f2d_next,f2d_obc,xflux2d,tmp
  integer  :: i,j,k,j1,j11,j22
!------------------------------------------------------------------------------|
       
!--------------------------------------------
! Set Scalar Value on Open Boundary
!--------------------------------------------
  if(iobcn > 0) then
    do i=1,iobcn
      j=i_obc_n(i)
      j1=next_obc(i)
      f2d=0.0_sp
      f2d_next=0.0_sp
      xflux2d=0.0_sp
      do k=1,kbm1
        f2d=f2d+f(j,k)*dz(j,k)
        f2d_next=f2d_next+fn(j1,k)*dz(j1,k)
        xflux2d=xflux2d+fflux_obc(i,k)*dz(j,k)
      end do
  
      if(uard_obcn(i) > 0.0_sp) then
        tmp=xflux2d+f2d*uard_obcn(i)
        f2d_obc=(f2d*dt(j)-tmp*deltat/art1(j))/d(j)
        do k=1,kbm1
          fn(j,k)=f2d_obc+(fn(j1,k)-f2d_next)
        end do
      else
        do k=1,kbm1
          fn(j,k) = f(j,k)-alpha_nudge*(f(j,k)-f_obc(i))
        end do
      end if
    end do
  endif

  return
  End Subroutine Bcond_Scal_OBC 
!==============================================================================|
!==============================================================================|


!==========================================================================
! Calculate Fluxes for Vertical Advection Equation                            
! MElem: number of cells
! c: scalar variable (1:MElem)
! w: velocity field at cell interfaces (1:MElem+1)
! note:  we dont use face normals to construct inflow/outflow
!        thus we add dissipation term instead of subtracting because 
!        positive velocity is up while computational coordinates increase
!        down towards bottom.  
!==========================================================================
  Subroutine Calc_VFlux(MElem,c,w,flux) 
  use mod_prec
  implicit none
  integer , intent(in ) :: MElem
  real(sp), intent(in ) ::  c(MElem)
  real(sp), intent(in ) ::  w(MElem+1) 
  real(sp), intent(out) ::  flux(MElem+1)
  real(sp) :: conv(MElem+1),diss(MElem+1)
  real(sp) :: cin(-1:MElem+2)
  real(sp) :: dis4
  integer  :: i

  !transfer to working array
  cin(1:MElem) = c(1:MElem)

  !surface bcs (no flux)
  cin(0)  =  -cin(1) 
  cin(-1) =  -cin(2)
  
  !bottom bcs (no flux)
  cin(MElem+1) = -cin(MElem) 
  cin(MElem+2) = -cin(MElem-1)

  !flux computation
  do i=1,MElem+1
    dis4    = .5*abs(w(i))
    conv(i) = w(i)*(cin(i)+cin(i-1))/2. 
    diss(i) = dis4*(cin(i)-cin(i-1)-lim(cin(i+1)-cin(i),cin(i-1)-cin(i-2))) 
    flux(i) = conv(i)+diss(i)
  end do

  End Subroutine Calc_VFlux
  
!==========================================================================
! Calculate LED Limiter L(u,v)  
!==========================================================================
  Function Lim(a,b)
  use mod_prec
  real(sp) lim,a,b
  real(sp) q,R
  real(sp) eps
  eps = epsilon(eps)
  
  q = 0. !1st order
  q = 1. !minmod
  q = 2. !van leer

  R = abs(   (a-b)/(abs(a)+abs(b)+eps) )**q
  lim = .5*(1-R)*(a+b)

  End Function Lim


End Module Scalar
