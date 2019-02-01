!!This mod is used for release dye concentration.

#  if defined (DYE_RELEASE)
MODULE MOD_DYE

   USE MOD_PREC
   USE MOD_INP
   USE CONTROL
   USE MOD_PAR
   IMPLICIT NONE
   SAVE
!
!--VARIABLES for SPECIFY DYE RELEASE                 
     LOGICAL  :: DYE_ON                        !!RELEASE DYE ACTIVE
     INTEGER  :: IINT_SPE_DYE_B                !!INTERNAL TIME STEP OF BEGIN RELEASE DYE
     INTEGER  :: IINT_SPE_DYE_E                !!INTERNAL TIME STEP END OF RELEASE DYE
     INTEGER  :: KSPE_DYE                      !!NUMBER OF SIGMA LAYER FOR SPECIFY DYE RELEASE 
     INTEGER  :: MSPE_DYE                      !!NUMBER OF NODE FOR SPECIFY DYE  RELEASE
     INTEGER, ALLOCATABLE  :: K_SPECIFY(:)    !!NO of sigma layer for specify dye release
     INTEGER, ALLOCATABLE  :: M_SPECIFY(:)    !!NO of node for specify dye release
     REAL(SP) :: DYE_SOURCE_TERM               !!Specify source term value of dye releasing
!--VARIABLES OF DYE     
     REAL(SP), ALLOCATABLE :: DYE(:,:)       !!DYE CONCENTRATION AT NODE
     REAL(SP), ALLOCATABLE :: DYEF(:,:)      !!DYE CONCENTRATION FROM PREVIOUS TIME
     REAL(SP), ALLOCATABLE :: DYEMEAN(:,:)   !!MEAN INITIAL DYE
     REAL(SP), ALLOCATABLE :: DYE_S(:,:)    !! temporary dye in modified upwind
     REAL(SP), ALLOCATABLE :: DYE_SF(:,:)    !! temporary dye in modified upwind

     REAL(SP), ALLOCATABLE :: WWWS(:,:)    !! temporary dye in modified upwind
     REAL(SP), ALLOCATABLE :: WWWSF(:,:)    !! temporary dye in modified upwind
     REAL(SP), ALLOCATABLE :: DTWWWS(:)    !! temporary dye in modified upwind
     REAL(SP), ALLOCATABLE :: zzzflux(:,:)    !! temporary total flux in corrected part
     REAL(SP), ALLOCATABLE :: beta(:,:)    !! temporary beta coefficient in corrected part

     REAL(SP), ALLOCATABLE :: betain(:,:)    !! temporary beta coefficient in corrected part
     REAL(SP), ALLOCATABLE :: betaout(:,:)    !! temporary beta coefficient in corrected part


    integer itera,ntera
     real ssss
   CONTAINS !------------------------------------------------------------------!
            ! ALLOC_VARS_DYE  : Allocate and Initialize Arrays of dye
            ! SET_DYE_PARA  : specify dye soerce tterm and parameter
            ! ADV_DYE       : Horizontal Advection/Diffusion of dye Variables  !
            ! VDIF_DYE      : Vertical Diffusion of dye Variables              !
            ! INITIAL_DYE   : Initialize for dye Variables                     !
            ! BCOND_DYE     : Boundary Conditions (River Flux) of DYE Variables!

!==============================================================================|
!    Allocate and Initialize Arrays of Dye                                     !
!==============================================================================|

   SUBROUTINE ALLOC_VARS_DYE

!==============================================================================!
   USE ALL_VARS
   IMPLICIT NONE
!#  if defined (MULTIPROCESSOR)
!   include "mpif.h"
!#  endif
!==============================================================================!


!==============================================================================!
!  ALLOCATE:                                                                   !
!==============================================================================!
   ALLOCATE(DYE(0:NTNode,KB))           ;DYE    = ZERO
   ALLOCATE(DYEF(0:NTNode,KB))          ;DYEF    = ZERO
   ALLOCATE(DYEMEAN(0:NTNode,KB))       ;DYEMEAN    = ZERO
   ALLOCATE(DYE_S(0:NTNode,KB))          ;DYE_S    = ZERO
   ALLOCATE(DYE_SF(0:NTNode,KB))          ;DYE_SF    = ZERO
   ALLOCATE(WWWS(0:NTNode,KB))          ;WWWS    = ZERO
   ALLOCATE(WWWSF(0:NTNode,KB))          ;WWWSF    = ZERO
   ALLOCATE(DTWWWS(0:NTNode))       ;DTWWWS    = ZERO

   ALLOCATE(zzzflux(0:NTNode,KB))          ;zzzflux    = ZERO 
   ALLOCATE(beta(0:NTNode,KB))       ;beta    = ZERO  
   ALLOCATE(betain(0:NTNode,KB))       ;betain    = ZERO
   ALLOCATE(betaout(0:NTNode,KB))       ;betaout    = ZERO

   RETURN
   END SUBROUTINE ALLOC_VARS_DYE


!==============================================================================|
!   Specify source term                |
!==============================================================================|
   SUBROUTINE SET_DYE_PARAM
   USE MOD_PREC
   USE CONTROL
   IMPLICIT NONE
   INTEGER  INTVEC(150),ISCAN,KTEMP
   CHARACTER(LEN=120) :: FNAME
   FNAME = "./"//trim(casename)//"_run.dat"
!------------------------------------------------------------------------------|
!     "DYE_ON"   !! 
!------------------------------------------------------------------------------|     
   ISCAN = SCAN_FILE(TRIM(FNAME),"DYE_ON",LVAL = DYE_ON)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DYE_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "IINT_SPE_DYE_B"   !! 
!------------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"IINT_SPE_DYE_B",ISCAL = IINT_SPE_DYE_B)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IINT_SPE_DYE_B: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "IINT_SPE_DYE_E"   !! 
!------------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"IINT_SPE_DYE_E",ISCAL = IINT_SPE_DYE_E)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IINT_SPE_DYE_E: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "KSPE_DYE"   !! 
!------------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"KSPE_DYE",ISCAL = KSPE_DYE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING KSPE_DYE: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "MSPE_DYE"   !! 
!------------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"MSPE_DYE",ISCAL = MSPE_DYE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING MSPE_DYE: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "K_SPECIFY"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"K_SPECIFY",IVEC =INTVEC ,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING K_SPECIFY: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= KSPE_DYE)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED K_SPECIFY IS NOT EQUAL TO KSPE_DYE' 
       WRITE(*,*)'KSPE_DYE: ',KSPE_DYE
       WRITE(*,*)'K_SPECIFY: ',INTVEC
     END IF
   END IF
  
   ALLOCATE(K_SPECIFY(KSPE_DYE)) ; K_SPECIFY=0
   K_SPECIFY(1:KSPE_DYE)= INTVEC(1:KSPE_DYE)
!------------------------------------------------------------------------------|
!     "M_SPECIFY"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"M_SPECIFY",IVEC =INTVEC ,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING M_SPECIFY: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= MSPE_DYE)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED M_SPECIFY IS NOT EQUAL TO MSPE_DYE' 
       WRITE(*,*)'MSPE_DYE: ',MSPE_DYE
       WRITE(*,*)'M_SPECIFY: ',INTVEC
     END IF
   END IF
  
   ALLOCATE(M_SPECIFY(MSPE_DYE)) ; M_SPECIFY=0
   M_SPECIFY(1:MSPE_DYE)= INTVEC(1:MSPE_DYE)
   
!------------------------------------------------------------------------------|
!     "DYE_SOURCE_TERM"   -CONCENTRATION OF DYE SOURCE
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"DYE_SOURCE_TERM",FSCAL = DYE_SOURCE_TERM)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DYE_SOURCE_TERM: ',ISCAN
     CALL PSTOP 
   END IF

!==============================================================================|
!            SCREEN REPORT OF SET DYE RELEASE VARIABlES                        !
!==============================================================================|
   IF(MSR) THEN  
     WRITE(IPT,*) '!                                                   !'     
     WRITE(IPT,*) '!------SPECIFY DYE RELEASE VARIABlES----------------!'     
     WRITE(IPT,*) '!                                                   !'     
     WRITE(IPT,*) '!  # DYE_ON              :',DYE_ON
     WRITE(IPT,*) '!  # IINT_SPE_DYE_B      :',IINT_SPE_DYE_B
     WRITE(IPT,*) '!  # IINT_SPE_DYE_E      :',IINT_SPE_DYE_E
     WRITE(IPT,*) '!  # KSPE_DYE            :',KSPE_DYE
     WRITE(IPT,*) '!  # K_SPECIFY           :',K_SPECIFY
     WRITE(IPT,*) '!  # MSPE_DYE            :',MSPE_DYE
     WRITE(IPT,*) '!  # M_SPECIFY           :',M_SPECIFY
   END IF
   RETURN
   END SUBROUTINE SET_DYE_PARAM


!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for DYE                 |
!==============================================================================|

   SUBROUTINE ADV_DYE               

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif   
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:NTNode,KB)     :: XFLUX,XFLUX_ADV
   REAL(SP), DIMENSION(NNode)           :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(NNode)           :: PDYEPX,PDYEPY,PDYEPXD,PDYEPYD,VISCOFF
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: DTIJ 
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
   REAL(SP) :: FACT,FM1
   REAL(SP) :: s1min, s1max, s2min, s2max,SMIN,SMAX,xxxx
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,KK,IP,II
# if defined (SPHERICAL)
   REAL(SP) :: ty,txpi,typi
   REAL(DP) :: XTMP1,XTMP
   REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
# endif
#  if defined (WET_DRY)
   INTEGER :: N_NTVE
#  endif

!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

     
!
!--Initialize Fluxes and dyef----------------------------------------------------------!
!
   XFLUX     = 0.0_SP
   XFLUX_ADV = 0.0_SP
   DYEF      = 0.0_SP
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
!!#  if !defined (WET_DRY)
   DO I=1,NCV
     I1=NTRG(I)
!     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
       DTIJ(I,K)=DT1(I1)*DZ1(I1,K)
       UVN(I,K) = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I)
     END DO
   END DO
!!#  else
!!   DO I=1,NCV
!!     I1=NTRG(I)
!!!     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!       DTIJ(I,K)=DT1(I1)*DZ1(I1,K)
!!       UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!     END DO
!!   END DO
!!#  endif

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=1,KBM1
      PDYEPX  = 0.0_SP 
      PDYEPY  = 0.0_SP 
      PDYEPXD = 0.0_SP 
      PDYEPYD = 0.0_SP
     DO I=1,NNode
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)

#    if defined (WET_DRY)
         IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 1)THEN
          FFD=0.5_SP*(DYE(I,K)+DYE(I2,K)-DYEMEAN(I,K)-DYEMEAN(I2,K))
          FF1=0.5_SP*(DYE(I,K)+DYE(I2,K))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 1 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(DYE(I1,K)+DYE(I,K)-DYEMEAN(I1,K)-DYEMEAN(I,K))
          FF1=0.5_SP*(DYE(I1,K)+DYE(I,K))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(DYE(I,K)+DYE(I,K)-DYEMEAN(I,K)-DYEMEAN(I,K))
          FF1=0.5_SP*(DYE(I,K)+DYE(I,K))
	 ELSE
          FFD=0.5_SP*(DYE(I1,K)+DYE(I2,K)-DYEMEAN(I1,K)-DYEMEAN(I2,K))
          FF1=0.5_SP*(DYE(I1,K)+DYE(I2,K))
	 END IF 
#    else	 
         FFD=0.5_SP*(DYE(I1,K)+DYE(I2,K)-DYEMEAN(I1,K)-DYEMEAN(I2,K))
         FF1=0.5_SP*(DYE(I1,K)+DYE(I2,K))
#    endif	 
	 
#        if defined (SPHERICAL)
         XTMP  = VX(I2)*TPI-VX(I1)*TPI
	 XTMP1 = VX(I2)-VX(I1)
	 IF(XTMP1 >  180.0_SP)THEN
	   XTMP = -360.0_SP*TPI+XTMP
	 ELSE IF(XTMP1 < -180.0_SP)THEN
	   XTMP =  360.0_SP*TPI+XTMP
	 END IF  
         TXPI=XTMP*COS(DEG2RAD*VY(I))
         TYPI=(VY(I1)-VY(I2))*TPI

         PDYEPX(I)=PDYEPX(I)+FF1*TYPI
         PDYEPY(I)=PDYEPY(I)+FF1*TXPI
         PDYEPXD(I)=PDYEPXD(I)+FFD*TYPI
         PDYEPYD(I)=PDYEPYD(I)+FFD*TXPI
#        else
         PDYEPX(I)=PDYEPX(I)+FF1*(VY(I1)-VY(I2))
         PDYEPY(I)=PDYEPY(I)+FF1*(VX(I2)-VX(I1))
         PDYEPXD(I)=PDYEPXD(I)+FFD*(VY(I1)-VY(I2))
         PDYEPYD(I)=PDYEPYD(I)+FFD*(VX(I2)-VX(I1))
#        endif
       END DO
       PDYEPX(I)=PDYEPX(I)/ART2(I)
       PDYEPY(I)=PDYEPY(I)/ART2(I)
       PDYEPXD(I)=PDYEPXD(I)/ART2(I)
       PDYEPYD(I)=PDYEPYD(I)/ART2(I)
     END DO
          
     IF(K == KBM1)THEN
       DO I=1,NNode
         PFPXB(I) = PDYEPX(I)
         PFPYB(I) = PDYEPY(I)
       END DO
     END IF

     DO I=1,NNode
!       PUPX(I)=0.0_SP
!       PUPY(I)=0.0_SP
!       PVPX(I)=0.0_SP
!       PVPY(I)=0.0_SP
!       J=1
!       I1=NBVE(I,J)
!       JTMP=NBVT(I,J)
!       J1=JTMP+1-(JTMP+1)/4*3
!       J2=JTMP+2-(JTMP+2)/4*3
!       X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
!       Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!       X22=XC(I1)
!       Y22=YC(I1)
!       X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
!       Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

!#      if defined (SPHERICAL)
!       TY  =0.5_SP*(Y11+Y33)
!       TXPI=(X33-X11)*TPI*COS(DEG2RAD*TY)
!       TYPI=(Y11-Y33)*TPI
!       PUPX(I)=PUPX(I)+U(I1,K)*typi
!       PUPY(I)=PUPY(I)+U(I1,K)*txpi
!       PVPX(I)=PVPX(I)+V(I1,K)*typi
!       PVPY(I)=PVPY(I)+V(I1,K)*txpi
!#      else
!       PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
!       PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
!       PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
!       PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
!#      endif

!       IF(ISONB(I) /= 0) THEN
!#        if defined (SPHERICAL)
!         TY=0.5_SP*(VY(I)+Y11)
!         TXPI=(X11-VX(I))*TPI*COS(DEG2RAD*TY)
!         TYPI=(VY(I)-Y11)*TPI
!         PUPX(I)=PUPX(I)+U(I1,K)*typi
!         PUPY(I)=PUPY(I)+U(I1,K)*txpi
!         PVPX(I)=PVPX(I)+V(I1,K)*typi
!         PVPY(I)=PVPY(I)+V(I1,K)*txpi
!#        else
!         PUPX(I)=PUPX(I)+U(I1,K)*(VY(I)-Y11)
!         PUPY(I)=PUPY(I)+U(I1,K)*(X11-VX(I))
!         PVPX(I)=PVPX(I)+V(I1,K)*(VY(I)-Y11)
!         PVPY(I)=PVPY(I)+V(I1,K)*(X11-VX(I))
!#        endif
!       END IF

!       DO J=2,NTVE(I)-1
!         I1=NBVE(I,J)
!         JTMP=NBVT(I,J)
!         J1=JTMP+1-(JTMP+1)/4*3
!         J2=JTMP+2-(JTMP+2)/4*3
!         X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
!         Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!         X22=XC(I1)
!         Y22=YC(I1)
!         X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
!         Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

!#        if defined (SPHERICAL)
!         TY=0.5_SP*(Y11+Y33)
!         txpi=(x33-x11)*TPI*COS(DEG2RAD*TY)
!         typi=(y11-y33)*TPI
!         PUPX(I)=PUPX(I)+U(I1,K)*typi
!         PUPY(I)=PUPY(I)+U(I1,K)*txpi
!         PVPX(I)=PVPX(I)+V(I1,K)*typi
!         PVPY(I)=PVPY(I)+V(I1,K)*txpi
!#        else
!         PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
!         PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
!         PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
!         PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
!#        endif
!       END DO
!       J=NTVE(I)
!       I1=NBVE(I,J)
!       JTMP=NBVT(I,J)
!       J1=JTMP+1-(JTMP+1)/4*3
!       J2=JTMP+2-(JTMP+2)/4*3
!       X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
!       Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!       X22=XC(I1)
!       Y22=YC(I1)
!       X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
!       Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

!#      if defined (SPHERICAL)
!       TY=0.5*(y11+y33)
!       txpi=(x33-x11)*TPI*cos(DEG2RAD*TY)
!       typi=(y11-y33)*TPI
!       pupx(i)=pupx(i)+u(i1,k)*typi
!       pupy(i)=pupy(i)+u(i1,k)*txpi
!       pvpx(i)=pvpx(i)+v(i1,k)*typi
!       pvpy(i)=pvpy(i)+v(i1,k)*txpi
!#      else
!       PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
!       PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
!       PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
!       PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
!#      endif

!       IF(ISONB(I) /= 0) THEN
!#      if defined (SPHERICAL)
!         TY=0.5*(y11+vy(i))
!         txpi=(vx(i)-x11)*TPI*cos(DEG2RAD*TY)
!         typi=(y11-vy(i))*TPI
!         PUPX(I)=PUPX(I)+U(I1,K)*typi
!         PUPY(I)=PUPY(I)+U(I1,K)*txpi
!         PVPX(I)=PVPX(I)+V(I1,K)*typi
!         PVPY(I)=PVPY(I)+V(I1,K)*txpi
!#        else
!         PUPX(I)=PUPX(I)+U(I1,K)*(Y11-VY(I))
!         PUPY(I)=PUPY(I)+U(I1,K)*(VX(I)-X11)
!         PVPX(I)=PVPX(I)+V(I1,K)*(Y11-VY(I))
!         PVPY(I)=PVPY(I)+V(I1,K)*(VX(I)-X11)
!#        endif
!       END IF
!       PUPX(I)=PUPX(I)/ART1(I)
!       PUPY(I)=PUPY(I)/ART1(I)
!       PVPX(I)=PVPX(I)/ART1(I)
!       PVPY(I)=PVPY(I)/ART1(I)
!       TMP1=PUPX(I)**2+PVPY(I)**2
!       TMP2=0.5_SP*(PUPY(I)+PVPX(I))**2
!       VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)

       VISCOFF(I)=VISCOFH(I,K)
       
     END DO
     IF(K == KBM1) THEN
       AH_BOTTOM(1:NNode) = HORCON*(FACT*VISCOFF(1:NNode) + FM1)
    END IF


     DO I=1,NCV_I
       IA=NIEC(I,1)
       IB=NIEC(I,2)
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
#      if defined (SPHERICAL)
       X1_DP=XIJE(I,1)
       Y1_DP=YIJE(I,1)
       X2_DP=XIJE(I,2)
       Y2_DP=YIJE(I,2)
       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
       XI=XII		
       XTMP  = XI*TPI-VX(IA)*TPI
       XTMP1 = XI-VX(IA)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF	 

       DXA=XTMP*COS(DEG2RAD*VY(IA))    
       DYA=(YI-VY(IA))*TPI
       XTMP  = XI*TPI-VX(IB)*TPI
       XTMP1 = XI-VX(IB)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF	 

       DXB=XTMP*COS(DEG2RAD*VY(IB)) 
       DYB=(YI-VY(IB))*TPI
#      else
       DXA=XI-VX(IA)
       DYA=YI-VY(IA)
       DXB=XI-VX(IB)
       DYB=YI-VY(IB)
#      endif
       FIJ1=DYE(IA,K)+DXA*PDYEPX(IA)+DYA*PDYEPY(IA)
       FIJ2=DYE(IB,K)+DXB*PDYEPX(IB)+DYB*PDYEPY(IB)

       s1min=MINVAL(DYE(nbsn(IA,1:NTSN(IA)-1),k))
       s1min=MIN(s1min, DYE(IA,K))
       s1max=MAXVAL(DYE(nbsn(IA,1:NTSN(IA)-1),k))
       s1max=MAX(s1max, DYE(IA,K))
       s2min=MINVAL(DYE(nbsn(IB,1:NTSN(IB)-1),k))
       s2min=MIN(s2min, DYE(IB,K))
       s2max=MAXVAL(DYE(nbsn(IB,1:NTSN(IB)-1),k))
       s2max=MAX(s2max, DYE(IB,K))
       if (FIJ1 < s1min) FIJ1=s1min
       if (FIJ1 > s1max) FIJ1=s1max
       if (FIJ2 < s2min) FIJ2=s2min
       if (FIJ2 > s2max) FIJ2=s2max

       UN=UVN(I,K)

       VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

       TXX=0.5_SP*(PDYEPXD(IA)+PDYEPXD(IB))*VISCOF
       TYY=0.5_SP*(PDYEPYD(IA)+PDYEPYD(IB))*VISCOF

       FXX=-DTIJ(I,K)*TXX*DLTYE(I)
       FYY= DTIJ(I,K)*TYY*DLTXE(I)

       EXFLUX=-UN*DTIJ(I,K)* &
          ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

       XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
       XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)
     END DO
  END DO !!SIGMA LOOP


!
!-Accumulate Fluxes at Boundary Nodes
!
# if defined (MULTIPROCESSOR)
  IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,XFLUX,XFLUX_ADV)
# endif

  DO K=1,KBM1
     IF(IOBCN > 0) THEN
       DO I=1,IOBCN
         I1=I_OBC_N(I)
         XFLUX_OBC(I,K)=XFLUX_ADV(I1,K)
       END DO
     END IF
  ENDDO
 

!   --------------------------------------------------------------------
!   The central difference scheme in vertical advection
!   --------------------------------------------------------------------
!   DO K=1,KBM1
!     DO I=1,NNode
!#    if defined (WET_DRY)
!       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
!#    endif
!       IF(K == 1) THEN
!         TEMP=-WTS(I,K+1)*(DYE(I,K)*DZ(K+1)+DYE(I,K+1)*DZ(K))/(DZ(K)+DZ(K+1))
!       ELSE IF(K == KBM1) THEN
!         TEMP= WTS(I,K)*(DYE(I,K)*DZ(K-1)+DYE(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1))
!       ELSE
!         TEMP= WTS(I,K)*(DYE(I,K)*DZ(K-1)+DYE(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1))-&
!               WTS(I,K+1)*(DYE(I,K)*DZ(K+1)+DYE(I,K+1)*DZ(K))/(DZ(K)+DZ(K+1))
!       END IF

!       IF(ISONB(I) == 2) THEN
!         XFLUX(I,K)=TEMP*ART1(I)/DZ(K)
!       ELSE
!         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/DZ(K)
!       END IF
!#    if defined (WET_DRY)
!       END IF
!#    endif
!     END DO
!   END DO  !! SIGMA LOOP

!   -------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------

!#   if defined (MPDATA)

!   --------------------------------------------------------------------------------
!   ---------------------------------------------------------------------------------
!   Songhu {
!   Using smolarkiewicz, P. K; A fully multidimensional positive definite advection
!   transport algorithm with small implicit diffusion, Journal of Computational
!   Physics, 54, 325-362, 1984
        
! The horizontal term of advection is neglected here

   DO K=1,KBM1
     DO I=1,NNode
       IF(ISONB(I) == 2) THEN
         XFLUX(I,K)=0.
       ENDIF
     ENDDO
   ENDDO

! Initialize variables of MPDATA
   DYE_S=0._SP
   DYE_SF=0._SP
   WWWS=0._SP
   WWWSF=0._SP
   DTWWWS=0._SP
   ZZZFLUX=0._SP
   BETA=0._SP
   BETAIN=0._SP
   BETAOUT=0._SP

!!   first loop for vertical upwind
!!   flux including horizontal and vertical upwind
   DO K=1,KBM1
     DO I=1,NNode
#    if defined (WET_DRY)
       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
#    endif
         IF(K == 1) THEN
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*DYE(I,K)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*DYE(I,K+1) &
                  +(WTS(I,K)+ABS(WTS(I,K)))*DYE(I,K)    
         ELSE IF(K == KBM1) THEN
           TEMP =  +(WTS(I,K)-ABS(WTS(I,K)))*DYE(I,K-1)     &
                   +(WTS(I,K)+ABS(WTS(I,K)))*DYE(I,K)
         ELSE
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*DYE(I,K)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*DYE(I,K+1) &
                  +(WTS(I,K)-ABS(WTS(I,K)))*DYE(I,K-1)     &
                  +(WTS(I,K)+ABS(WTS(I,K)))*DYE(I,K)
         END IF
         TEMP = 0.5_SP*TEMP 

         IF(K == 1)THEN
           SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMAX = MAX(SMAX,DYE(I,K+1),DYE(I,K))
           SMIN = MIN(SMIN,DYE(I,K+1),DYE(I,K))
         ELSE IF(K == KBM1)THEN
           SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMAX = MAX(SMAX,DYE(I,K-1),DYE(I,K))
           SMIN = MIN(SMIN,DYE(I,K-1),DYE(I,K))
         ELSE
           SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMAX = MAX(SMAX,DYE(I,K+1),DYE(I,K-1),DYE(I,K))
           SMIN = MIN(SMIN,DYE(I,K+1),DYE(I,K-1),DYE(I,K))
         END IF

        ZZZFLUX(I,K)= TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K)/ART1(I)*(DTI/DT(I))/DZ(I,K) 
        XXXX= ZZZFLUX(I,K)*DT(I)/DTFA(I)+DYE(I,K)-DYE(I,K)*DT(I)/DTFA(I) 

        BETA(I,K)=0.5*(1.-SIGN(1.,XXXX)) * (SMAX-DYE(I,K))/(ABS(XXXX)+1.E-10) &
                 +0.5*(1.-SIGN(1.,-XXXX)) * (DYE(I,K)-SMIN)/(ABS(XXXX)+1.E-10)

         DYE_SF(I,K)=DYE(I,K)-MIN(1.,BETA(I,K))*XXXX

#    if defined (WET_DRY)
       END IF
#    endif
     END DO
   END DO  !! SIGMA LOOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   NTERA = 4
   DO ITERA=1,NTERA   !! Smolaricizw Loop 
     IF(ITERA == 1)THEN
       WWWSF = WTS
       DYE_S =DYE_SF
       DTWWWS = DT
     ELSE
       WWWSF = WWWS
       DYE_S =DYE_SF
       DTWWWS = DTFA
     END IF
     DO K=2,KBM1
       DO I=1,NNode
         TEMP=ABS(WWWSF(I,K))-DTI*(WWWSF(I,K))*(WWWSF(I,K))/DZ(I,K)/DTWWWS(I)
         WWWS(I,K)=TEMP*(DYE_S(I,K-1)-DYE_S(I,K))/(ABS(DYE_S(I,K-1))+     &
	           ABS(DYE_S(I,K))+1.E-14)

         IF(TEMP < 0.0_SP .OR. DYE_S(I,K) == 0.0_SP)THEN 
           WWWS(I,K)=0.0_SP 
         END IF
       END DO
     END DO
     DO I=1,NNode
       WWWS(I,1)=0.0_SP
     END DO


     DO I=1,NNode
       SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),1))
       SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),1))
       SMAX = MAX(SMAX,DYE(I,2),DYE(I,1))
       SMIN = MIN(SMIN,DYE(I,2),DYE(I,1))
 
       TEMP=0.5*((WWWS(I,2)+ABS(WWWS(I,2)))*DYE_S(I,2))*(DTI/DTFA(I))/DZ(I,1)
       BETAIN(I,1)=(SMAX-DYE_S(I,1))/(TEMP+1.E-10)

       TEMP=0.5*((WWWS(I,1)+ABS(WWWS(I,1)))*DYE_S(I,1)-    &
            (WWWS(I,2)-ABS(WWWS(I,2)))*DYE_S(I,1))*(DTI/DTFA(I))/DZ(I,1)
       BETAOUT(I,1)=(DYE_S(I,1)-SMIN)/(TEMP+1.E-10)

       WWWSF(I,1)=0.5*MIN(1.,BETAOUT(I,1))*(WWWS(I,1)+ABS(WWWS(I,1))) + &
                  0.5*MIN(1.,BETAIN(I,1))*(WWWS(I,1)-ABS(WWWS(I,1)))
     END DO

     DO K=2,KBM1-1
       DO I=1,NNode
         SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
         SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
         SMAX = MAX(SMAX,DYE(I,K+1),DYE(I,K-1),DYE(I,K))
         SMIN = MIN(SMIN,DYE(I,K+1),DYE(I,K-1),DYE(I,K))
 
         TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*DYE_S(I,K+1)-   &
	      (WWWS(I,K)-ABS(WWWS(I,K)))*DYE_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
         BETAIN(I,K)=(SMAX-DYE_S(I,K))/(TEMP+1.E-10)

         TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)-         &
	      (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*DYE_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
         BETAOUT(I,K)=(DYE_S(I,K)-SMIN)/(TEMP+1.E-10)

         WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
                    0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
       END DO
     END DO
     
     K=KBM1
     DO I=1,NNode
       SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
       SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
       SMAX = MAX(SMAX,DYE(I,K-1),DYE(I,K))
       SMIN = MIN(SMIN,DYE(I,K-1),DYE(I,K))
 
       TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*DYE_S(I,K+1)-   &
            (WWWS(I,K)-ABS(WWWS(I,K)))*DYE_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
       BETAIN(I,K)=(SMAX-DYE_S(I,K))/(TEMP+1.E-10)

       TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)-         &
	    (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*DYE_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
       BETAOUT(I,K)=(DYE_S(I,K)-SMIN)/(TEMP+1.E-10)

       WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
                  0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
     END DO
     WWWS=WWWSF 

     DO K=1,KBM1
       DO I=1,NNode
#      if defined (WET_DRY)
         IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
#      endif
           IF(K == 1) THEN
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*DYE_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*DYE_S(I,K+1) &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)
           ELSE IF(K == KBM1) THEN
             TEMP =  +(WWWS(I,K)-ABS(WWWS(I,K)))*DYE_S(I,K-1)     &
                     +(WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)
           ELSE
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*DYE_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*DYE_S(I,K+1) &
                    +(WWWS(I,K)-ABS(WWWS(I,K)))*DYE_S(I,K-1)     &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)
           END IF
           TEMP = 0.5_SP*TEMP
           DYE_SF(I,K)=(   DYE_S(I,K)-TEMP*(DTI/DTFA(I))/DZ(I,K)   ) 
#      if defined (WET_DRY)
         END IF
#      endif
       END DO
     END DO  !! SIGMA LOOP
   END DO  !! Smolarvizw Loop


!--------------------------------------------------------------------------
! End of smolarkiewicz upwind loop
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!#  endif

!# if ! defined(MPDATA)
!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------
!   DO K=1,KBM1
!     DO I=1,NNode
!#    if defined (WET_DRY)
!       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
!#    endif
!       IF(K == 1) THEN
!         TEMP=-WTS(I,K+1)*(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
!	      (DZ(I,K)+DZ(I,K+1))
!       ELSE IF(K == KBM1) THEN
!         TEMP= WTS(I,K)*(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/    &
!	      (DZ(I,K)+DZ(I,K-1))
!       ELSE
!         TEMP= WTS(I,K)*(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/    &
!	      (DZ(I,K)+DZ(I,K-1))-&
!               WTS(I,K+1)*(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
!	      (DZ(I,K)+DZ(I,K+1))
!       END IF

!       IF(ISONB(I) == 2) THEN
!         XFLUX(I,K)=TEMP*ART1(I)    !/DZ(K)
!       ELSE
!         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)    !/DZ(K)
!       END IF
!#    if defined (WET_DRY)
!       END IF
!#    endif
!     END DO
!   END DO  !! SIGMA LOOP


!
!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
!   IF(POINT_ST_TYPE == 'calculated') THEN
!     IF(INFLOW_TYPE == 'node') THEN
!       IF(NUMQBC > 0) THEN
!         DO J=1,NUMQBC
!           JJ=INODEQ(J)
!           STPOINT=SDIS(J)
!           DO K=1,KBM1
!!             XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT/DZ(K)
!             XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT
!           END DO
!         END DO
!       END IF
!     ELSE IF(INFLOW_TYPE == 'edge') THEN
!       IF(NUMQBC > 0) THEN
!         DO J=1,NUMQBC
!           J1=N_ICELLQ(J,1)
!           J2=N_ICELLQ(J,2)
!           STPOINT=SDIS(J) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
!           DO K=1,KBM1
!!             XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT/DZ(K)
!!             XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT/DZ(K)
!             XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
!             XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
!           END DO
!         END DO
!       END IF
!     END IF
!   END IF

!#  endif
!
!--Update Dye-------------------------------------------------------------!
!

   DO I=1,NNode
#  if defined (WET_DRY)
     IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1 )THEN
#  endif
     DO K=1,KBM1
!#  if defined (MPDATA)
     DYEF(I,K)=DYE_SF(I,K)
!#  else     
!     DYEF(I,K)=(DYE(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
!#  endif
     END DO
#  if defined (WET_DRY)
     ELSE
     DO K=1,KBM1
       DYEF(I,K)=DYE(I,K)
     END DO
     END IF
#  endif
   END DO
!   
!---- specify the source term------------------------------------------|
!
!   IF(IINT.GE.IINT_SPE_DYE_B.AND.IINT.LE.IINT_SPE_DYE_E) THEN
!     DO KK=1, KSPE_DYE
!        K=K_SPECIFY(KK)
!        DO J=1,MSPE_DYE
!           I = M_SPECIFY(J)
!           DYEF(I,K)= DYE_SOURCE_TERM
!        ENDDO
!     ENDDO
!   ENDIF

   RETURN
   END SUBROUTINE ADV_DYE
!==============================================================================|



!==============================================================================|
!     this subroutine is used to calculate the dye                             !
!     by including vertical diffusion implicitly.                              !
!==============================================================================|

   SUBROUTINE VDIF_DYE(F)                

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE BCS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   INTEGER :: I,K,J,KI,KK
   REAL(DP) :: TMP,TMP1,TMP2,TMP3,QTMP,GW,ZDEP,FKH,UMOLPR
   REAL(SP), DIMENSION(0:NTNode,KB)  :: F
   REAL(DP), DIMENSION(NNode,KB)     :: FF,AF,CF,VHF,VHPF,RAD
   REAL(DP), DIMENSION(NNode)        :: KHBOTTOM,WFSURF,SWRADF


   UMOLPR = UMOL*1.E0_SP

!------------------------------------------------------------------------------!
!                                                                              !
!        the following section solves the equation                             !
!         dti*(kh*f')'-f=-fb                                                   !
!                                                                              !
!------------------------------------------------------------------------------!


   DO K = 2, KBM1
     DO I = 1, NNode
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
!         FKH=0.0_SP
!         DO J=1,NTVE(I)
!           FKH=FKH+KH(NBVE(I,J),K)
!         END DO
!         FKH=FKH/FLOAT(NTVE(I))
         FKH = KH(I,K)

         IF(K == KBM1) THEN
           KHBOTTOM(I)=FKH
         END IF

!         AF(I,K-1)=-DTI*(FKH+UMOLPR)/(DZ(K-1)*DZZ(K-1)*D(I)*D(I))
!         CF(I,K)=-DTI*(FKH+UMOLPR)/(DZ(K)*DZZ(K-1)*D(I)*D(I))
         AF(I,K-1)=-DTI*(FKH+UMOLPR)/(DZ(I,K-1)*DZZ(I,K-1)*D(I)*D(I))
         CF(I,K)=-DTI*(FKH+UMOLPR)/(DZ(I,K)*DZZ(I,K-1)*D(I)*D(I))
       END IF
     END DO
   END DO


!------------------------------------------------------------------------------!
!     the net heat flux input.                                                 !
!     the method shown below can be used when we turn off the                  !
!     body force in subroutine advt. be sure this method could                 !
!     cause the surface overheated if the vertical resolution                  !
!     is not high enough.                                                      !
!------------------------------------------------------------------------------!

     DO I = 1, NNode
       SWRADF(I)= 0.0_SP
       WFSURF(I)=0.0_SP
       DO K=1,KB
         RAD(I,K)=0.0_SP
       END DO
     END DO


!------------------------------------------------------------------------------!
!   surface bcs; wfsurf                                                        !
!------------------------------------------------------------------------------!

   DO I = 1, NNode
#  if !defined (WET_DRY)
     IF (D(I) > 0.0_SP) THEN
#  else
     IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
       VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
       VHPF(I,1) = -DTI *(WFSURF(I)-SWRADF(I) &
                   +RAD(I,1)-RAD(I,2)) / (-DZ(I,1)*D(I)) - F(I,1)
       VHPF(I,1) = VHPF(I,1) / (AF(I,1)-1.)
     END IF
   END DO

   DO K = 2, KBM2
     DO I = 1, NNode
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
         VHPF(I,K)=1./ (AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
         VHF(I,K) = AF(I,K) * VHPF(I,K)
         VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-DBLE(F(I,K)) &
                     +DTI*(RAD(I,K)-RAD(I,K+1))/(D(I)*DZ(I,K)))*VHPF(I,K)
       END IF
     END DO
   END DO


   DO  K = 1, KBM1
     DO  I = 1, NNode
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
         FF(I,K) = F(I,K)
       END IF
     END DO
   END DO

   DO I = 1, NNode
#  if !defined (WET_DRY)
     IF (D(I) > 0.0_SP .AND.ISONB(I) /= 2) THEN
#  else
     IF(ISWET_NODE_CURRENTSTEP(I) == 1 .AND.ISONB(I) /= 2)THEN
#  endif
       TMP1=PFPXB(I)*COS(SITA_GD(I))+PFPYB(I)*SIN(SITA_GD(I))
       TMP2=AH_BOTTOM(I)*PHPN(I)
       TMP3=KHBOTTOM(I)+UMOLPR+AH_BOTTOM(I)*PHPN(I)*PHPN(I)
       TMP=TMP1*TMP2/TMP3*(KHBOTTOM(I)+UMOLPR)
! --- Huang change
!       IF (TMP1 > 0.0_SP) TMP=0.0_SP
       TMP=0.0_SP
! change end
       GW=0.0_SP
!       IF(IBFW > 0) THEN
!         DO J=1,IBFW
!           IF(I == NODE_BFW(J)) THEN
!!             QTMP=-(F(I,KBM1)*D1(I)*DZ(KBM1)*BFWDIS(J))/ &
!!                   (D1(I)*DZ(KBM1)*ART1(I)+BFWDIS(J))
!!             GW=DTI/D1(I)/DZ(KBM1)*QTMP
!             QTMP=-(F(I,KBM1)*D(I)*DZ(I,KBM1)*BFWDIS(J))/ &
!                   (D(I)*DZ(I,KBM1)*ART1(I)+BFWDIS(J))
!             GW=DTI/D(I)/DZ(I,KBM1)*QTMP
!             TMP=0.0_SP
!           END IF
!         END DO
!       END IF

       FF(I,KBM1) = ((CF(I,KBM1)*VHPF(I,KBM2)-FF(I,KBM1)-GW &
               +DTI*(RAD(I,KBM1)-RAD(I,KB)-TMP)/(D(I)*DZ(I,KBM1))) &
                /(CF(I,KBM1)*(1.-VHF(I,KBM2))-1.))
     END IF
   END DO

   DO  K = 2, KBM1
     KI = KB - K
     DO  I = 1, NNode
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP .AND.ISONB(I) /= 2) THEN
#  else
       IF(ISWET_NODE_CURRENTSTEP(I) == 1 .AND.ISONB(I) /= 2)THEN
#  endif
         FF(I,KI) = (VHF(I,KI)*FF(I,KI+1)+VHPF(I,KI))
       END IF
     END DO
   END DO

   DO I = 1, NNode
#  if defined (WET_DRY)
     IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1 )then
#  endif
       DO K = 1, KBM1
         F(I,K) = FF(I,K)
       END DO
#  if defined (WET_DRY)
     END IF
#  endif
   END DO

!   
!---- specify the source term------------------------------------------|
!
   IF(IINT.GE.IINT_SPE_DYE_B.AND.IINT.LE.IINT_SPE_DYE_E) THEN

    IF(SERIAL)THEN
     DO KK=1, KSPE_DYE
       K=K_SPECIFY(KK)
       DO J=1,MSPE_DYE
         I = M_SPECIFY(J)
         F(I,K)= DYE_SOURCE_TERM
       END DO
     END DO
    END IF 

# if defined (MULTIPROCESSOR)
    IF(PAR)THEN
     DO KK=1, KSPE_DYE
        K=K_SPECIFY(KK)
        DO J=1,MSPE_DYE
           I = NLID(M_SPECIFY(J))
           F(I,K)= DYE_SOURCE_TERM
        ENDDO
     ENDDO
    END IF 
#  endif
   ENDIF

   RETURN
   END SUBROUTINE VDIF_DYE
!==============================================================================|
!==============================================================================|
!  AVERAGE THE dye                                                     |
!==============================================================================|

   SUBROUTINE AVER_DYE

!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP):: AVE_DYE,SMAX,SMIN
   INTEGER :: I,J,K
!==============================================================================|

   IF(H_TYPE == 'body_h') GO TO 100
   DO I=1,NNode
     IF(IOBCN > 0)THEN
       DO J=1,IOBCN
         IF(I == I_OBC_N(J))GO TO 200
       END DO
     END IF  	 
     IF(NUMQBC > 0)THEN
       DO J=1,NUMQBC
         IF(INFLOW_TYPE == 'node' .AND. I == INODEQ(J))GO TO 200
         IF(INFLOW_TYPE == 'edge' .AND. &
           (I == N_ICELLQ(J,1) .OR. I == N_ICELLQ(J,2)))GO TO 200
       END DO
     END IF
     DO K=1,KBM1
       SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
       SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))

       IF(K == 1)THEN
         SMAX = MAX(SMAX,(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K+1)))
         SMIN = MIN(SMIN,(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K+1)))
       ELSE IF(K == KBM1)THEN
         SMAX = MAX(SMAX,(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K-1)))
         SMIN = MIN(SMIN,(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K-1)))
       ELSE
         SMAX = MAX(SMAX,(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K-1)), &
                 (DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/   &
		 (DZ(I,K)+DZ(I,K+1)))
         SMIN = MIN(SMIN,(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K-1)), &
                 (DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/   &
		 (DZ(I,K)+DZ(I,K+1)))
       END IF

       IF(SMIN-DYEF(I,K) > 0.0_SP)DYEF(I,K) = SMIN
       IF(DYEF(I,K)-SMAX > 0.0_SP)DYEF(I,K) = SMAX

     END DO
200 CONTINUE
   END DO

100 CONTINUE
   RETURN
   END SUBROUTINE AVER_DYE
!==============================================================================|


!==============================================================================|
!    Initialize dye fields (dye)                                               !
!    Calculate Mean Fields (DYEMEAN)                                           !
!==============================================================================|

   SUBROUTINE INITIAL_DYE

!==============================================================================!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
!# if defined (MULTIPROCESSOR)
!   include "mpif.h"
!# endif
   INTEGER :: I,K,ierr
   real(sp), allocatable :: temp_dye(:,:),temp_dyemean(:,:)

!==============================================================================!




   DO I=1,NNode
     DO K=1,KB
        DYE(I,K)=0.0_SP
        DYEMEAN(I,K)=0.0_SP
     END DO
   END DO

!! read global dye data into temporary arrays
!   allocate(temp_dye(0:NNodeGL,kb))
!   allocate(temp_dyemean(0:NNodeGL,kb))
!   open(78,file='gom_bioini.dat')
!   do i=1,NNodeGL
!      read(78,*)(temp_dye(i,k),k=1,kb)
!   enddo
!!   do i=1,NNodeGL
!!      read(78,*)(temp_dyemean(i,k),k=1,kb)
!!   enddo
!   temp_dyemean = temp_dye
!   close(78)

!! broadcast to all processors
!   call mpi_bcast(temp_dye,kb*(NNodeGL+1), mpi_f,0,mpi_comm_world,ierr)
!   call mpi_bcast(temp_dyemean,kb*(NNodeGL+1), mpi_f,0,mpi_comm_world,ierr)

!! transform to local arrays 
!   if(serial)then
!     dye = temp_dye
!     dyemean = temp_dyemean
!   end if

!# if defined (MULTIPROCESSOR)
!   if(par)then
!     do i=1,NNode
!        dye(i,:) = temp_dye(ngid(i),:)
!        dyemean(i,:) = temp_dyemean(ngid(i),:)
!     end do
!     if(par)call exchange(nc,NTNode,kb,myid,nprocs,dye,dyemean)
!   end if

!   deallocate(temp_dye)
!   deallocate(temp_dyemean)
!# endif
    
   RETURN
   END SUBROUTINE INITIAL_DYE
!==============================================================================|
!==============================================================================|
!   Set Boundary Conditions on DYE                                             |
!==============================================================================|

   SUBROUTINE BCOND_DYE     

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: S2D,S2D_NEXT,S2D_OBC,T2D,T2D_NEXT,T2D_OBC,XFLUX2D,TMP
   INTEGER  :: I,J,K,J1,J11,J22
   REAL(SP) :: SMAX,SMIN
!------------------------------------------------------------------------------|


       
   IF(IOBCN > 0) THEN


!
!  SET dye CONDITIONS ON OUTER BOUNDARY
!
     DO I=1,IOBCN
       J=I_OBC_N(I)
       J1=NEXT_OBC(I)
       S2D=0.0_SP
       S2D_NEXT=0.0_SP
       XFLUX2D=0.0_SP
       DO K=1,KBM1
         S2D=S2D+DYE(J,K)*DZ(J,K)
         S2D_NEXT=S2D_NEXT+DYEF(J1,K)*DZ(J1,K)
         XFLUX2D=XFLUX2D+XFLUX_OBC(I,K)                 !*DZ(K)
       END DO
 
       IF(UARD_OBCN(I) > 0.0_SP) THEN
         TMP=XFLUX2D+S2D*UARD_OBCN(I)
         S2D_OBC=(S2D*DT(J)-TMP*DTI/ART1(J))/D(J)
         DO K=1,KBM1
!           DYEF(J,K)=S2D_OBC+(DYEF(J1,K)-S2D_NEXT)  !!bug 2 
           DYEF(J,K)=DYEF(J1,K)
          END DO

         DO K=1,KBM1
           SMAX = MAXVAL(DYE(NBSN(J,1:NTSN(J)),K))
           SMIN = MINVAL(DYE(NBSN(J,1:NTSN(J)),K))

           IF(K == 1)THEN
            SMAX = MAX(SMAX,(DYE(J,K)*DZ(J,K+1)+DYE(J,K+1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K+1)))
            SMIN = MIN(SMIN,(DYE(J,K)*DZ(J,K+1)+DYE(J,K+1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K+1)))
           ELSE IF(K == KBM1)THEN
            SMAX = MAX(SMAX,(DYE(J,K)*DZ(J,K-1)+DYE(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)))
            SMIN = MIN(SMIN,(DYE(J,K)*DZ(J,K-1)+DYE(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)))
           ELSE
            SMAX = MAX(SMAX,(DYE(J,K)*DZ(J,K-1)+DYE(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)),                             &
                   (DYE(J,K)*DZ(J,K+1)+DYE(J,K+1)*DZ(J,K))/           &
                   (DZ(J,K)+DZ(J,K+1)))
            SMIN = MIN(SMIN,(DYE(J,K)*DZ(J,K-1)+DYE(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)),                             &
                   (DYE(J,K)*DZ(J,K+1)+DYE(J,K+1)*DZ(J,K))/           &
                   (DZ(J,K)+DZ(J,K+1)))
           END IF

           IF(SMIN-DYEF(J,K) > 0.0_SP) DYEF(J,K) = SMIN
           IF(DYEF(J,K)-SMAX > 0.0_SP) DYEF(J,K) = SMAX

         END DO

        ELSE
         DO K=1,KBM1
           DYEF(J,K)=DYE(J,K)
         END DO
       END IF
     END DO
   END IF

!
!--SET BOUNDARY CONDITIONS-----------------------------------------------------|
!
!   DO K=1,KBM1
!     DYE(0,K)=0.0_SP
!   END DO

   RETURN
   END SUBROUTINE BCOND_DYE
!==============================================================================|
     
END MODULE MOD_DYE
# else
   SUBROUTINE DUM_DYE
   END SUBROUTINE DUM_DYE

# endif

