!   THIS PROGRAM LINKS 1D BIOLOGICAL CALCULATION TO FVCOM 3D COMPUTAITION
MODULE MOD_BIO_3D 
#  if defined (BioGen)
   USE ALL_VARS
   USE MOD_NCDIO
   USE BCS
   USE MOD_OBCS
   USE LIMS
   USE CONTROL
# if defined (MULTIPROCESSOR)
   USE MOD_PAR   
# endif
#  if defined (WET_DRY) 
   USE MOD_WD
#  endif
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
#if defined (NETCDF_IO)
   use netcdf
#endif
   use mod_types
   use mod_utils
   USE MOD_1D
   USE MOD_PHYTOPLANKTON !,   ONLY: BIO_P,NNP,INP,IRRAD0,PARFRAC
   USE MOD_ZOOPLANKTON !,     ONLY: BIO_Z,NNZ,INZ
   USE MOD_BACTERIA !,        ONLY: BIO_B,NNB,INB
   USE MOD_DETRITUS !,        ONLY: BIO_D,NND,IND
   USE MOD_DOM !,             ONLY: BIO_DOM,NNM,INM
   USE MOD_NUTRIENT !,        ONLY: BIO_N,NNN,INN
   USE MOD_PARAMETER

   IMPLICIT NONE
   SAVE
   REAL(SP), ALLOCATABLE ::  BIO_ALL(:,:,:)        !3D BIO_VARIABLES
   REAL(SP), ALLOCATABLE ::  BIO_F(:,:,:)          !FORECASTED VARIABLES
   REAL(SP), ALLOCATABLE ::  BIO_MEAN(:,:,:)       !MEAN VARIABLES
   REAL(SP), ALLOCATABLE ::  XFLUX_OBCB(:,:,:)     !OPEN BOUNDARY FLUX
   REAL(SP), ALLOCATABLE ::  BIO_MEANN(:,:,:)      !MEAN IN CELLS
!************   FOR NETCDF OUTPUT
   CHARACTER(LEN=120) :: fldnam,CDFNAME_BIO   !BioGen
   INTEGER :: bio_ofid,ii,mm,stck_bio          !BioGen
   INTEGER :: ntve_vid,nbve_vid

   integer :: node_didb,nele_didb
   integer :: scl_didb,siglay_didb,siglev_didb
   integer :: three_didb,four_didb
   integer :: time_didb

   !--Grid Variables
   integer :: nprocs_vidb,partition_vidb
   integer :: idens_vidb
   integer :: x_vidb,y_vidb,lat_vidb,lon_vidb

   integer :: nv_vidb,nbe_vidb
   integer :: aw0_vidb,awx_vidb,awy_vidb
   integer :: a1u_vidb,a2u_vidb
   integer :: siglay_vidb,siglev_vidb,siglay_shift_vidb
  
   !--Flow Variables 
   integer :: time_vidb
   integer :: iint_vidb
   integer :: u_vidb
   integer :: v_vidb
   integer :: wd_vidb
   integer :: ww_vidb
   integer :: s1_vidb
   integer :: t1_vidb
   integer :: el_vidb
   integer :: h_vidb
   integer :: km_vidb
   integer :: kh_vidb
   integer :: ua_vidb
   integer :: va_vidb
   integer :: d_vidb

   INTEGER, ALLOCATABLE :: trcsid(:)
   REAL(SP),ALLOCATABLE :: BIO_VAR_MEAN(:,:,:)

   CONTAINS !------------------------------------------------------------------!
            ! BIO_3D1D             :ADVANCE TMODEL USING GOTM LIBRARIES        !
            ! BIO_ADV              :ADVECTION OF BIOLOGICAL STATE VARIABLES    !
            ! BIO_BCOND            :BOUNDARY CONDITION                         !
            ! BIO_NETCDF_HEADER    :NETCDF HEADER WRITER                       !
            ! BIO_OUT_NETCDF       :NETCDF OUTPUT                              !
            ! BIO_EXCHANGE         :INTERPROCESSOR EXCHANGE                    !
            ! BIO_INITIAL          :INITIALIZATION                             !
            ! BIO_HOT_START        :HOT_START FROM BIO_RESTART.NC              !
            !------------------------------------------------------------------!
SUBROUTINE BIO_3D1D
     IMPLICIT NONE
     SAVE
     INTEGER  :: I,J,K,L
     REAL(SP) :: SPCP,ROSEA,SPRO !,BIO_VAR_MEAN(NNode,KBM1,NTT)
     REAL(SP) :: DEPTH_Z(KB)
     SPCP  = 4.2174E3_SP                        !HEAT SPECIFIC CAPACITY
     ROSEA = 1.023E3_SP                         !RHO OF SEA WATER
     SPRO=SPCP*ROSEA
!     BIO_VAR_MEAN = 0.0_SP
!---------------------------                    !MAIN LOOP OVER ELEMENTS

#    if defined (ONE_D_MODEL)
     DO I=NNode,NNode
#    else
     DO I=1,NNode
#    endif
       DO K=1,KBM1                              !3D TO 1D FIELD
         DO J=1,NNN
            BIO_N(K,J)=BIO_ALL(I,K,J+INN-1)
         END DO
         DO J=1,NNP
            BIO_P(K,J)=BIO_ALL(I,K,J+INP-1)
         END DO
         DO J=1,NNZ
            BIO_Z(K,J)=BIO_ALL(I,K,J+INZ-1)
         END DO
         DO J=1,NNM
            BIO_DOM(K,J)=BIO_ALL(I,K,J+INM-1)
         END DO
         DO J=1,NNB
            BIO_B(K,J)=BIO_ALL(I,K,J+INB-1)
         END DO
         DO J=1,NND
            BIO_D(K,J)=BIO_ALL(I,K,J+IND-1)
         END DO
         DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS
         DELTA_Z(K)=DZZ(I,K)*D(I)                  !DISTANCE BETWEEN LAYERS
         DEPTH_Z(K)=Z(I,K)*D(I)                    !LAYER CENTER DEPTH
         IRRAD0=-SWRAD(I)*PARFRAC*SPRO/RAMP      !PAR FRACTION
         L_NH4N=30._SPP                          !NITRIFICATION USE	 
         T_BIO(K)=T1(I,K)
       END DO                                    !K=1,KB
       T_STEP=DTI
       CALL ZOOPLANKTON
       CALL PHYTOPLANKTON
       CALL BACTERIA
       CALL DETRITUS
       CALL DOM
       CALL NUTRIENT
        DO K=1,KBM1                                !1D TO 3D FIELD
          DO J=1,NNN
             BIO_ALL(I,K,J+INN-1)=BIO_N(K,J)
          END DO
          DO J=1,NNP
             BIO_ALL(I,K,J+INP-1)=BIO_P(K,J)
          END DO
          DO J=1,NNZ
             BIO_ALL(I,K,J+INZ-1)=BIO_Z(K,J)
          END DO
          DO J=1,NNM
             BIO_ALL(I,K,J+INM-1)=BIO_DOM(K,J)
          END DO
          DO J=1,NNB
             BIO_ALL(I,K,J+INB-1)=BIO_B(K,J)
          END DO
          DO J=1,NND
             BIO_ALL(I,K,J+IND-1)=BIO_D(K,J)
          END DO
        END DO 
      KM_BIO(:)=KH(I,:)
      BIO_VAR(1:KBV,1:NTT)=BIO_ALL(I,1:KBV,1:NTT)
      CALL BIO_MIXING 
      BIO_ALL(I,1:KBV,1:NTT)=BIO_VAR(1:KBV,1:NTT)
      END DO !I=1,NNode 
#    if !defined (ONE_D_MODEL)
#    if defined (MULTIPROCESSOR)
       CALL BIO_EXCHANGE
#    endif
      CALL BIO_ADV
#    if defined (MULTIPROCESSOR)
      CALL BIO_EXCHANGE
#    endif
      CALL BIO_BCOND
      BIO_ALL=BIO_F                                   !UPDATE
#    endif
!    end if defined 1D
      WHERE (BIO_ALL < 0.001) BIO_ALL=0.001

#if defined (NETCDF_IO)
      IF(MOD(IINT-1,CDF_INT)==0) CALL BIO_OUT_NETCDF
#endif
!endif if defined (NETCDF_IO)
END SUBROUTINE BIO_3D1D

!=============================================================================!
   SUBROUTINE BIO_ADV  
!=============================================================================!
!                                                                             !
!   This subroutine is used to calculate the horizontal advection and         !
!   and diffusion terms for the state variables of the adjustable biomodel    !
!=============================================================================!

   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
# if defined (MULTIPROCESSOR)
   USE MOD_PAR   
# endif
#  if defined (WET_DRY) 
   USE MOD_WD
#  endif
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
   IMPLICIT NONE

   REAL(SP), DIMENSION(0:NTNode,KB,ntt)  :: XFLUX,RF,XFLUX_ADV
   REAL(SP), DIMENSION(NNode)           :: PUPX,PUPY,PVPX,PVPY
   REAL(SP), DIMENSION(NNode)           :: PFPX,PFPY,PFPXD,PFPYD,VISCOFF
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: DTIJ
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: UVN
   REAL(SP) :: FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP
   REAL(SP) :: FACT,FM1
   REAL(SP) :: TT,TTIME,STPOINT
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,JTMP,K,JJ,N1
   REAL(SP) :: WQM1MIN, WQM1MAX, WQM2MIN, WQM2MAX

   REAL(SP), ALLOCATABLE :: DWDIS(:,:,:)  !!WATER QUALITY DISCHARGE DATA
   REAL(SP), ALLOCATABLE :: WDIS(:,:)     !!FRESH WATER QUALITY AT CURRENT TIME

# if defined (SPHERICAL)
   REAL(SP) :: ty,txpi,typi
   REAL(DP) :: XTMP,XTMP1
   REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
   REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
# endif
#  if defined (MPDATA)
   REAL(SP) :: WQMMIN,WQMMAX,XXXX
   REAL(SP), DIMENSION(0:NTNode,KB)     :: WQM_S    !! temporary salinity in modified upwind
   REAL(SP), DIMENSION(0:NTNode,KB)     :: WQM_SF   !! temporary salinity in modified upwind
   REAL(SP), DIMENSION(0:NTNode,KB)     :: WWWS     
   REAL(SP), DIMENSION(0:NTNode,KB)     :: WWWSF   
   REAL(SP), DIMENSION(0:NTNode)        :: DTWWWS  
   REAL(SP), DIMENSION(0:NTNode,KB)     :: ZZZFLUX !! temporary total flux in corrected part
   REAL(SP), DIMENSION(0:NTNode,KB)     :: BETA    !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:NTNode,KB)     :: BETAIN  !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:NTNode,KB)     :: BETAOUT !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:NTNode,KB)     :: BIO_FRESH    !! for source term

   INTEGER ITERA, NTERA
#  endif

   ALLOCATE(WDIS(NUMQBC,NTT))     ;WDIS      = ZERO
 !  ALLOCATE(DWDIS(NUMQBC,NB,NQTIME))

!------------------------------------------------------------------------------
   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure')THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Fluxes-----------------------------------------------------------
!
   XFLUX = 0.0_SP
   XFLUX_ADV = 0.0_SP
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------
!
!!#  if !defined (WET_DRY)
   DO I=1,NCV
     I1=NTRG(I)
!     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
       UVN(I,K)=V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I) 
     END DO
   END DO
!!#  else
!!   DO I=1,NCV
!!     I1=NTRG(I)
!!!     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
!!       UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!     END DO
!!   END DO
!!#  endif

   TTIME=THOUR

   RF = 0.0_SP

!--Calculate the Advection and Horizontal Diffusion Terms----------------------

   DO N1=1,NTT
     DO K=1,KBM1
       PFPX  = 0.0_SP
       PFPY  = 0.0_SP
       PFPXD = 0.0_SP
       PFPYD = 0.0_SP

       DO I=1,NNode
         DO J=1,NTSN(I)-1
           I1=NBSN(I,J)
           I2=NBSN(I,J+1)

#    if defined (WET_DRY)
         IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 1)THEN
          FFD=0.5_SP*(BIO_ALL(I,K,N1)+BIO_ALL(I2,K,N1)           &
	      -BIO_MEAN(I,K,N1)-BIO_MEAN(I2,K,N1))
          FF1=0.5_SP*(BIO_ALL(I,K,N1)+BIO_ALL(I2,K,N1))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 1 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I,K,N1)           &
	      -BIO_MEAN(I1,K,N1)-BIO_MEAN(I,K,N1))
          FF1=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I,K,N1))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=BIO_ALL(I,K,N1)-BIO_MEAN(I,K,N1)
          FF1=BIO_ALL(I,K,N1)
	 ELSE
          FFD=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I2,K,N1)          &
	      -BIO_MEAN(I1,K,N1)-BIO_MEAN(I2,K,N1))
          FF1=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I2,K,N1))
	 END IF 
#    else	 
           FFD=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I2,K,N1)          &
               -BIO_MEAN(I1,K,N1)-BIO_MEAN(I2,K,N1))
           FF1=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I2,K,N1))
#    endif	 
	 
#          if defined (SPHERICAL)
           XTMP  = VX(I2)*TPI-VX(I1)*TPI
           XTMP1 = VX(I2)-VX(I1)
	   IF(XTMP1 >  180.0_SP)THEN
	     XTMP = -360.0_SP*TPI+XTMP
	   ELSE IF(XTMP1 < -180.0_SP)THEN
	     XTMP =  360.0_SP*TPI+XTMP
	   END IF  
           TXPI=XTMP*COS(DEG2RAD*VY(I))
           TYPI=(VY(I1)-VY(I2))*TPI
           PFPX(I)=PFPX(I)+FF1*TYPI
           PFPY(I)=PFPY(I)+FF1*TXPI
           PFPXD(I)=PFPXD(I)+FFD*TYPI
           PFPYD(I)=PFPYD(I)+FFD*TXPI
#          else
           PFPX(I)=PFPX(I)+FF1*(VY(I1)-VY(I2))
           PFPY(I)=PFPY(I)+FF1*(VX(I2)-VX(I1))
           PFPXD(I)=PFPXD(I)+FFD*(VY(I1)-VY(I2))
           PFPYD(I)=PFPYD(I)+FFD*(VX(I2)-VX(I1))
#          endif
         END DO
         PFPX(I)=PFPX(I)/ART2(I)
         PFPY(I)=PFPY(I)/ART2(I)
         PFPXD(I)=PFPXD(I)/ART2(I)
         PFPYD(I)=PFPYD(I)/ART2(I)
       END DO

      IF(K == KBM1)THEN
        DO I=1,NNode
          PFPXB(I) = PFPX(I)
          PFPYB(I) = PFPY(I)
        END DO
      END IF


       DO I=1,NNode
!         PUPX(I) = 0.0_SP 
!         PUPY(I) = 0.0_SP
!         PVPX(I) = 0.0_SP
!         PVPY(I) = 0.0_SP
!         J=1
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
!         X1_DP=VX(I)
!         Y1_DP=VY(I)
!         X2_DP=VX(NV(I1,J1))
!         Y2_DP=VY(NV(I1,J1))
!         CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
!         X11=X11_TMP
!         Y11=Y11_TMP
!         X2_DP=VX(NV(I1,J2))
!         Y2_DP=VY(NV(I1,J2))
!         CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
!         X33=X33_TMP
!         Y33=Y33_TMP

!         XTMP  = X33*TPI-X11*TPI
!         XTMP1 = X33-X11
!         IF(XTMP1 >  180.0_SP)THEN
!           XTMP = -360.0_SP*TPI+XTMP
!         ELSE IF(XTMP1 < -180.0_SP)THEN
!           XTMP =  360.0_SP*TPI+XTMP
!         END IF	 
!         TXPI=XTMP*COS(DEG2RAD*VY(I))
!         TYPI=(Y11-Y33)*TPI
!         PUPX(I)=PUPX(I)+U(I1,K)*typi
!         PUPY(I)=PUPY(I)+U(I1,K)*txpi
!         PVPX(I)=PVPX(I)+V(I1,K)*typi
!         PVPY(I)=PVPY(I)+V(I1,K)*txpi
!#      else
!         PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
!         PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
!         PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
!         PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
!#      endif

!         IF(ISONB(I) /= 0) THEN
!#          if defined (SPHERICAL)
!           XTMP  = X11*TPI-VX(I)*TPI
!           XTMP1 = X11-VX(I)
!           IF(XTMP1 >  180.0_SP)THEN
!	     XTMP = -360.0_SP*TPI+XTMP
!           ELSE IF(XTMP1 < -180.0_SP)THEN
!	     XTMP =  360.0_SP*TPI+XTMP
!           END IF  
!           TXPI=XTMP*COS(DEG2RAD*VY(I))
!           TYPI=(VY(I)-Y11)*TPI
!           PUPX(I)=PUPX(I)+U(I1,K)*typi
!           PUPY(I)=PUPY(I)+U(I1,K)*txpi
!           PVPX(I)=PVPX(I)+V(I1,K)*typi
!           PVPY(I)=PVPY(I)+V(I1,K)*txpi
!#          else
!           PUPX(I)=PUPX(I)+U(I1,K)*(VY(I)-Y11)
!           PUPY(I)=PUPY(I)+U(I1,K)*(X11-VX(I))
!           PVPX(I)=PVPX(I)+V(I1,K)*(VY(I)-Y11)
!           PVPY(I)=PVPY(I)+V(I1,K)*(X11-VX(I))
!#          endif
!         END IF

!         DO J=2,NTVE(I)-1
!           I1=NBVE(I,J)
!           JTMP=NBVT(I,J)
!           J1=JTMP+1-(JTMP+1)/4*3
!           J2=JTMP+2-(JTMP+2)/4*3
!           X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
!           Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!           X22=XC(I1)
!           Y22=YC(I1)
!           X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
!           Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

!#          if defined (SPHERICAL)
!           X1_DP=VX(I)
!           Y1_DP=VY(I)
!           X2_DP=VX(NV(I1,J1))
!           Y2_DP=VY(NV(I1,J1))
!           CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
!  	   X11=X11_TMP
!	   Y11=Y11_TMP
!           X2_DP=VX(NV(I1,J2))
!           Y2_DP=VY(NV(I1,J2))
!           CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
!	   X33=X33_TMP
!	   Y33=Y33_TMP

!           XTMP  = X33*TPI-X11*TPI
!           XTMP1 = X33-X11
!           IF(XTMP1 >  180.0_SP)THEN
!	     XTMP = -360.0_SP*TPI+XTMP
!           ELSE IF(XTMP1 < -180.0_SP)THEN
!	     XTMP =  360.0_SP*TPI+XTMP
!  	   END IF  
!           TXPI=XTMP*COS(DEG2RAD*VY(I))
!           TYPI=(Y11-Y33)*TPI
!           PUPX(I)=PUPX(I)+U(I1,K)*typi
!           PUPY(I)=PUPY(I)+U(I1,K)*txpi
!           PVPX(I)=PVPX(I)+V(I1,K)*typi
!           PVPY(I)=PVPY(I)+V(I1,K)*txpi
!#          else
!           PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
!           PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
!           PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
!           PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
!#          endif
!         END DO
!         J=NTVE(I)
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
!         X1_DP=VX(I)
!         Y1_DP=VY(I)
!         X2_DP=VX(NV(I1,J1))
!         Y2_DP=VY(NV(I1,J1))
!         CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
!         X11=X11_TMP
!         Y11=Y11_TMP
!         X2_DP=VX(NV(I1,J2))
!         Y2_DP=VY(NV(I1,J2))
!         CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
!         X33=X33_TMP
!         Y33=Y33_TMP

!         XTMP  = X33*TPI-X11*TPI
!         XTMP1 = X33-X11
!         IF(XTMP1 >  180.0_SP)THEN
!           XTMP = -360.0_SP*TPI+XTMP
!         ELSE IF(XTMP1 < -180.0_SP)THEN
!           XTMP =  360.0_SP*TPI+XTMP
!         END IF	 
!         TXPI=XTMP*COS(DEG2RAD*VY(I))
!         typi=(y11-y33)*TPI
!         pupx(i)=pupx(i)+u(i1,k)*typi
!         pupy(i)=pupy(i)+u(i1,k)*txpi
!         pvpx(i)=pvpx(i)+v(i1,k)*typi
!         pvpy(i)=pvpy(i)+v(i1,k)*txpi
!#        else
!         PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
!         PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
!         PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
!         PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
!#        endif

!         IF(ISONB(I) /= 0) THEN
!#          if defined (SPHERICAL)
!           TY=0.5*(y11+vy(i))  !Not in adv_s.F
!           XTMP  = VX(I)*TPI-X11*TPI
!           XTMP1 = VX(I)-X11
!           IF(XTMP1 >  180.0_SP)THEN
!	     XTMP = -360.0_SP*TPI+XTMP
!           ELSE IF(XTMP1 < -180.0_SP)THEN
!	     XTMP =  360.0_SP*TPI+XTMP
!	   END IF  
!           TXPI=XTMP*COS(DEG2RAD*VY(I))
!           typi=(y11-vy(i))*TPI
!           PUPX(I)=PUPX(I)+U(I1,K)*typi
!           PUPY(I)=PUPY(I)+U(I1,K)*txpi
!           PVPX(I)=PVPX(I)+V(I1,K)*typi
!           PVPY(I)=PVPY(I)+V(I1,K)*txpi
!#          else
!           PUPX(I)=PUPX(I)+U(I1,K)*(Y11-VY(I))
!           PUPY(I)=PUPY(I)+U(I1,K)*(VX(I)-X11)
!           PVPX(I)=PVPX(I)+V(I1,K)*(Y11-VY(I))
!           PVPY(I)=PVPY(I)+V(I1,K)*(VX(I)-X11)
!#          endif
!         END IF
!         PUPX(I)=PUPX(I)/ART1(I)
!         PUPY(I)=PUPY(I)/ART1(I)
!         PVPX(I)=PVPX(I)/ART1(I)
!         PVPY(I)=PVPY(I)/ART1(I)
!         TMP1=PUPX(I)**2+PVPY(I)**2
!         TMP2=0.5_SP*(PUPY(I)+PVPX(I))**2
!         VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)

         VISCOFF(I)=VISCOFH(I,K)  !CALCULATED IN viscofh.F
       END DO
     IF(K == KBM1) THEN
       AH_BOTTOM(1:NNode) = HORCON*(FACT*VISCOFF(1:NNode) + FM1)
     END IF

       DO I=1,NCV_I
         IA=NIEC(I,1)
         IB=NIEC(I,2)
         XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
         YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
#        if defined (SPHERICAL)
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
!         dya=(yi-vy(ia))*TPI
         DYA=(YI-VY(IA))*TPI
!         TY=0.5*(yi+vy(ib))
         XTMP  = XI*TPI-VX(IB)*TPI
         XTMP1 = XI-VX(IB)
         IF(XTMP1 >  180.0_SP)THEN
           XTMP = -360.0_SP*TPI+XTMP
         ELSE IF(XTMP1 < -180.0_SP)THEN
           XTMP =  360.0_SP*TPI+XTMP
         END IF	 
         DXB=XTMP*COS(DEG2RAD*VY(IB))  
!         dyb=(yi-vy(ib))*TPI
         DYB=(YI-VY(IB))*TPI
#        else
         DXA=XI-VX(IA)
         DYA=YI-VY(IA)
         DXB=XI-VX(IB)
         DYB=YI-VY(IB)
#        endif
         FIJ1=BIO_ALL(IA,K,N1)+DXA*PFPX(IA)+DYA*PFPY(IA)
         FIJ2=BIO_ALL(IB,K,N1)+DXB*PFPX(IB)+DYB*PFPY(IB)

         WQM1MIN=MINVAL(BIO_ALL(NBSN(IA,1:NTSN(IA)-1),K,N1))
         WQM1MIN=MIN(WQM1MIN, BIO_ALL(IA,K,N1))
         WQM1MAX=MAXVAL(BIO_ALL(NBSN(IA,1:NTSN(IA)-1),K,N1))
         WQM1MAX=MAX(WQM1MAX, BIO_ALL(IA,K,N1))
         WQM2MIN=MINVAL(BIO_ALL(NBSN(IB,1:NTSN(IB)-1),K,N1))
         WQM2MIN=MIN(WQM2MIN, BIO_ALL(IB,K,N1))
         WQM2MAX=MAXVAL(BIO_ALL(NBSN(IB,1:NTSN(IB)-1),K,N1))
         WQM2MAX=MAX(WQM2MAX, BIO_ALL(IB,K,N1))
         IF(FIJ1 < WQM1MIN) FIJ1=WQM1MIN
         IF(FIJ1 > WQM1MAX) FIJ1=WQM1MAX
         IF(FIJ2 < WQM2MIN) FIJ2=WQM2MIN
         IF(FIJ2 > WQM2MAX) FIJ2=WQM2MAX

         UN=UVN(I,K)
  
         VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

         TXX=0.5_SP*(PFPXD(IA)+PFPXD(IB))*VISCOF
         TYY=0.5_SP*(PFPYD(IA)+PFPYD(IB))*VISCOF

         FXX=-DTIJ(I,K)*TXX*DLTYE(I)
         FYY= DTIJ(I,K)*TYY*DLTXE(I)

         EXFLUX=-UN*DTIJ(I,K)*                           &
                ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+          &
                 (1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP   &
                +FXX+FYY
 
         XFLUX(IA,K,N1)=XFLUX(IA,K,N1)+EXFLUX
         XFLUX(IB,K,N1)=XFLUX(IB,K,N1)-EXFLUX

       XFLUX_ADV(IA,K,N1)=XFLUX_ADV(IA,K,N1)+(EXFLUX-FXX-FYY)
       XFLUX_ADV(IB,K,N1)=XFLUX_ADV(IB,K,N1)-(EXFLUX-FXX-FYY)

       END DO !to NNode

#    if defined (SPHERICAL)
#    if defined (NORTHPOLE)
!     CALL ADV_T_XY(XFLUX(:,:,N1),XFLUX_ADV(:,:,N1),PTPX,PTPY,PTPXD,PTPYD,VISCOFF,K)
#    endif
#    endif  

     END DO !to KBM1
   END DO !to nnt

!
!-Accumulate Fluxes at Boundary Nodes
!
# if defined (MULTIPROCESSOR)
      DO N1=1,NTT
!      IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,       &
!                            XFLUX(:,:,I))    
   IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,       &
                            XFLUX(:,:,N1),XFLUX_ADV(:,:,N1))
      END DO
# endif

!#  if !defined (MPDATA)
   DO N1=1,NTT
     DO K=1,KBM1
        IF(IOBCN > 0) THEN
          DO I=1,IOBCN
            I1=I_OBC_N(I)
            XFLUX_OBCB(I,K,N1)=XFLUX_ADV(I1,K,N1)
          END DO
        END IF
      END DO
    END DO
!# endif


!#  if !defined (MPDATA)
   DO N1=1,ntt
#  if !defined (MPDATA)
!
!--Calculate the Vertical Terms------------------------------------------------
!
     DO K=1,KBM1
       DO I=1,NNode
#      if defined (WET_DRY)
       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
#      endif
         IF(K == 1) THEN  !Is there any violation ?
           TEMP=-WTS(I,K+1)*(BIO_ALL(I,K,N1)*DZ(I,K+1)+BIO_ALL(I,K+1,N1)*DZ(I,K))/   &
	         (DZ(I,K)+DZ(I,K+1))
         ELSE IF(K == KBM1) THEN
           TEMP=WTS(I,K)*(BIO_ALL(I,K,N1)*DZ(I,K-1)+BIO_ALL(I,K-1,N1)*DZ(I,K))/      &
	         (DZ(I,K)+DZ(I,K-1))
         ELSE
           TEMP=WTS(I,K)*(BIO_ALL(I,K,N1)*DZ(I,K-1)+BIO_ALL(I,K-1,N1)*DZ(I,K))/      &
	         (DZ(I,K)+DZ(I,K-1))-  &
                WTS(I,K+1)*(BIO_ALL(I,K,N1)*DZ(I,K+1)+BIO_ALL(I,K+1,N1)*DZ(I,K))/    &
		 (DZ(I,K)+DZ(I,K+1))
         END IF

!
!--Total Fluxes ---------------------------------------------------------------
!
         IF(ISONB(I) == 2) THEN
           XFLUX(I,K,N1)=TEMP*ART1(I)
         ELSE
           XFLUX(I,K,N1)=XFLUX(I,K,N1)+TEMP*ART1(I)
         END IF
#    if defined (WET_DRY)
       END IF
#    endif
       END DO  !i=1,NNode
     END DO    !k=1,kbm1

!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
     IF(POINT_ST_TYPE == 'calculated') THEN
       IF(INFLOW_TYPE == 'node') THEN
         IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
             JJ=INODEQ(J)
             STPOINT=WDIS(J,N1)
             DO K=1,KBM1
               XFLUX(JJ,K,N1)=XFLUX(JJ,K,N1) - QDIS(J)*VQDIST(J,K)*STPOINT
             END DO
           END DO
         END IF
       ELSE IF(INFLOW_TYPE == 'edge') THEN
         IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
             J1=N_ICELLQ(J,1)
             J2=N_ICELLQ(J,2)
             STPOINT=WDIS(J,N1) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
             DO K=1,KBM1
               XFLUX(J1,K,N1)=XFLUX(J1,K,N1)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
               XFLUX(J2,K,N1)=XFLUX(J2,K,N1)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
             END DO
           END DO
         END IF
       END IF
     END IF

#  else
!--------------------------------------------------------------------------------
!   S. HU
!   Using smolarkiewicz, P. K; A fully multidimensional positive definite advection
!   TEMPport algorithm with small implicit diffusion, Journal of Computational
!   Physics, 54, 325-362, 1984
!-----------------------------------------------------------------        

	BIO_FRESH=BIO_ALL(:,:,N1)


   IF(POINT_ST_TYPE == 'calculated') THEN
     IF(INFLOW_TYPE == 'node') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           JJ=INODEQ(J)
           STPOINT=SDIS(J)
           DO K=1,KBM1
!	    S1_FRESH(JJ,K)=SDIS(J)   !NEED TO CHANGE FROM S TO BIO
            XFLUX(JJ,K,N1)=XFLUX(JJ,K,N1) - QDIS(J)*VQDIST(J,K)*STPOINT
           END DO
         END DO
       END IF
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           J1=N_ICELLQ(J,1)
           J2=N_ICELLQ(J,2)
           STPOINT=SDIS(J) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
           DO K=1,KBM1
!             S1_FRESH(J1,K)=SDIS(J)  !NEED BIO CONCENTRATION
!             S1_FRESH(J1,K)=SDIS(J)
             XFLUX(J1,K,N1)=XFLUX(J1,K,N1)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
             XFLUX(J2,K,N1)=XFLUX(J2,K,N1)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
           END DO
         END DO
       END IF
     END IF
   END IF
!

! The horizontal term of advection is neglected here
   DO K=1,KBM1
     DO I=1,NNode
       IF(ISONB(I) == 2) THEN
         XFLUX(I,K,N1)=0.
       ENDIF
     END DO
   END DO

! Initialize variables of MPDATA
   WQM_S=0._SP
   WQM_SF=0._SP
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
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*BIO_ALL(I,K,N1)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*BIO_ALL(I,K+1,N1) &
                  +(WTS(I,K)+ABS(WTS(I,K)))*BIO_ALL(I,K,N1)    
         ELSE IF(K == KBM1) THEN
           TEMP = +(WTS(I,K)-ABS(WTS(I,K)))*BIO_ALL(I,K-1,N1)     &
                  +(WTS(I,K)+ABS(WTS(I,K)))*BIO_ALL(I,K,N1)
         ELSE
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*BIO_ALL(I,K,N1)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*BIO_ALL(I,K+1,N1) &
                  +(WTS(I,K)-ABS(WTS(I,K)))*BIO_ALL(I,K-1,N1)     &
                  +(WTS(I,K)+ABS(WTS(I,K)))*BIO_ALL(I,K,N1)
         END IF
         TEMP = 0.5_SP*TEMP 

         IF(K /= 1)THEN
           WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
           WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
!           WQMMAX = MAX(WQMMAX,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_FRESH(I,K,N1))
!           WQMMIN = MIN(WQMMIN,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_FRESH(I,K,N1))
           WQMMAX = MAX(WQMMAX,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_FRESH(I,K))
           WQMMIN = MIN(WQMMIN,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_FRESH(I,K))
         ELSE
           WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
           WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
!           WQMMAX = MAX(WQMMAX,BIO_ALL(I,K+1,N1),BIO_FRESH(I,K,N1))
!           WQMMIN = MIN(WQMMIN,BIO_ALL(I,K+1,N1),BIO_FRESH(I,K,N1))
           WQMMAX = MAX(WQMMAX,BIO_ALL(I,K+1,N1),BIO_FRESH(I,K))
           WQMMIN = MIN(WQMMIN,BIO_ALL(I,K+1,N1),BIO_FRESH(I,K))
         END IF

         ZZZFLUX(I,K) = TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K,N1)/ART1(I)*(DTI/DT(I))/DZ(I,K) 
         XXXX = ZZZFLUX(I,K)*DT(I)/DTFA(I)+BIO_ALL(I,K,N1)-BIO_ALL(I,K,N1)*DT(I)/DTFA(I) 

         BETA(I,K)=0.5*(1.-SIGN(1.,XXXX)) * (WQMMAX-BIO_ALL(I,K,N1))/(ABS(XXXX)+1.E-10) &
                  +0.5*(1.-SIGN(1.,-XXXX)) * (BIO_ALL(I,K,N1)-WQMMIN)/(ABS(XXXX)+1.E-10)

         WQM_SF(I,K)=BIO_ALL(I,K,N1)-MIN(1.,BETA(I,K))*XXXX

#    if defined (WET_DRY)
       END IF
#    endif
     END DO
   END DO  !! SIGMA LOOP

!----------------------------------------------------------------------------------------
   NTERA = 4
   DO ITERA=1,NTERA   !! Smolaricizw Loop 
     IF(ITERA == 1)THEN
       WWWSF  = WTS
       WQM_S   = WQM_SF
       DTWWWS = DT
     ELSE
       WWWSF  = WWWS
       WQM_S   = WQM_SF
       DTWWWS = DTFA
     END IF
     DO K=2,KBM1
       DO I=1,NNode
         TEMP=ABS(WWWSF(I,K))-DTI*(WWWSF(I,K))*(WWWSF(I,K))/DZ(I,K)/DTWWWS(I)
         WWWS(I,K)=TEMP*(WQM_S(I,K-1)-WQM_S(I,K))/(ABS(WQM_S(I,K-1))+ABS(WQM_S(I,K))+1.E-14)
 
         IF(TEMP < 0.0_SP .OR. WQM_S(I,K) == 0.0_SP)THEN 
           WWWS(I,K)=0. 
         END IF
       END DO 
     END DO

     DO I=1,NNode
       WWWS(I,1)=0.
     END DO

     DO I=1,NNode
       WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),1,N1))
       WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),1,N1))
       WQMMAX = MAX(WQMMAX,BIO_ALL(I,2,N1),BIO_ALL(I,1,N1),BIO_FRESH(I,1))
       WQMMIN = MIN(WQMMIN,BIO_ALL(I,2,N1),BIO_ALL(I,1,N1),BIO_FRESH(I,1))
 
       TEMP=0.5*((WWWS(I,2)+ABS(WWWS(I,2)))*WQM_S(I,2))*(DTI/DTFA(I))/DZ(I,1)
       BETAIN(I,1)=(WQMMAX-WQM_S(I,1))/(TEMP+1.E-10)

       TEMP=0.5*((WWWS(I,1)+ABS(WWWS(I,1)))*WQM_S(I,1)-        &
	           (WWWS(I,2)-ABS(WWWS(I,2)))*WQM_S(I,1))*(DTI/DTFA(I))/DZ(I,1)
       BETAOUT(I,1)=(WQM_S(I,1)-WQMMIN)/(TEMP+1.E-10)

       WWWSF(I,1)=0.5*MIN(1.,BETAOUT(I,1))*(WWWS(I,1)+ABS(WWWS(I,1))) + &
                    0.5*MIN(1.,BETAIN(I,1))*(WWWS(I,1)-ABS(WWWS(I,1)))
     END DO

     DO K=2,KBM1
       DO I=1,NNode
         WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
         WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
         WQMMAX = MAX(WQMMAX,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_FRESH(I,K))
         WQMMIN = MIN(WQMMIN,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_FRESH(I,K))
 
         TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1)-  &
	           (WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
         BETAIN(I,K)=(WQMMAX-WQM_S(I,K))/(TEMP+1.E-10)

         TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)-        &
	           (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
         BETAOUT(I,K)=(WQM_S(I,K)-WQMMIN)/(TEMP+1.E-10)

         WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
                    0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
       END DO
     END DO

     WWWS=WWWSF 

     DO K=1,KBM1
       DO I=1,NNode
#      if defined (WET_DRY)
         IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
#      endif
           IF(K == 1) THEN
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1) &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
           ELSE IF(K == KBM1) THEN
             TEMP = +(WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1)     &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
           ELSE
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1) &
                    +(WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1)     &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
           END IF
           TEMP = 0.5_SP*TEMP
           WQM_SF(I,K)=(WQM_S(I,K)-TEMP*(DTI/DTFA(I))/DZ(I,K)) 
#      if defined (WET_DRY)
         END IF
#      endif
       END DO
     END DO  !! SIGMA LOOP
   END DO  !! Smolarvizw Loop
!--------------------------------------------------------------------------
! End of smolarkiewicz upwind loop
!--------------------------------------------------------------------------
#  endif

!--Update Variables--------------------------------
!

     DO I = 1,NNode
#    if defined (WET_DRY)
       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1 )THEN
#      endif
         DO K = 1, KBM1
       XFLUX(I,K,N1) = XFLUX(I,K,N1) - RF(I,K,N1)*ART1(I)     !/DZ(K)
#        if !defined (MPDATA) 
           BIO_F(I,K,N1)=(BIO_ALL(I,K,N1)-XFLUX(I,K,N1)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*   &
                         (DT(I)/D(I)) 
#        else
           BIO_F(I,K,N1)=WQM_SF(I,K)
#    endif  		 			 
         END DO 
#      if defined (WET_DRY)
       ELSE
         DO K=1,KBM1
           BIO_F(I,K,N1)=BIO_ALL(I,K,N1)
         END DO
       END IF
#      endif
     END DO

   END DO !do N1=1,ntt
   RETURN
   END SUBROUTINE BIO_ADV



   SUBROUTINE BIO_BCOND
!==============================================================================|
!   Set Boundary Conditions for BioGen                                         |
!==============================================================================|

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: T2D,T2D_NEXT,T2D_OBC,XFLUX2D,TMP,RAMP_TS
   INTEGER  :: I,J,K,J1,J11,J22,NCON2,N1
   REAL(SP), ALLOCATABLE :: WDIS(:,:)     !!FRESH WATER QUALITY AT CURRENT TIME
   REAL(SP) ::WQMMAX,WQMMIN
   ALLOCATE(WDIS(NUMQBC,ntt))     ;WDIS      = ZERO
!------------------------------------------------------------------------------|


!
!--SET CONDITIONS FOR FRESH WATER INFLOW---------------------------------------|
!
   POINT_ST_TYPE = 'NONE'  !TEMPORAL, NO RIVER BIOLOGY YET
   IF(POINT_ST_TYPE == 'specified') THEN
     IF(NUMQBC > 0) THEN
       IF(INFLOW_TYPE == 'node') THEN
         DO I=1,NUMQBC
           J11=INODEQ(I)
           DO K=1,KBM1
             DO N1=1,NTT
               BIO_F(J11,K,N1) = WDIS(I,N1)
             END DO
           END DO
         END DO
       ELSE IF(INFLOW_TYPE == 'edge') THEN
         DO I=1,NUMQBC
           J11=N_ICELLQ(I,1)
           J22=N_ICELLQ(I,2)
           DO K=1,KBM1
             DO N1=1,NTT
               BIO_F(J11,K,N1)=WDIS(I,N1)
               BIO_F(J22,K,N1)=WDIS(I,N1)
             END DO
           END DO
         END DO
       END IF
     END IF
   END IF

       
   IF(IOBCN > 0) THEN
!
!  SET CONDITIONS ON OUTER BOUNDARY
!
   RAMP_TS = TANH(FLOAT(IINT)/FLOAT(IRAMP+1))
     DO N1=1,NTT
       DO I=1,IOBCN
         J=I_OBC_N(I)
         J1=NEXT_OBC(I)
         T2D=0.0_SP
         T2D_NEXT=0.0_SP
         XFLUX2D=0.0_SP
         DO K=1,KBM1
           T2D=T2D+BIO_ALL(J,K,N1)*DZ(J,K)
           T2D_NEXT=T2D_NEXT+BIO_F(J1,K,N1)*DZ(J1,K)
           XFLUX2D=XFLUX2D+XFLUX_OBCB(I,K,N1)           !*DZ(K)
         END DO
         IF(UARD_OBCN(I) > 0.0_SP) THEN
           TMP=XFLUX2D+T2D*UARD_OBCN(I)
           T2D_OBC=(T2D*DT(J)-TMP*DTI/ART1(J))/D(J)
           DO K=1,KBM1
            BIO_ALL(J,K,N1)=T2D_OBC+(BIO_ALL(J1,K,N1)-T2D_NEXT)
!#  if !defined (MPDATA)	
!            taf(J,K,N1)=T2D_OBC+(taf(J1,K,N1)-T2D_NEXT)
!!               taf(J,K,N1) = taf(J1,K,N1)
!#  else
!               taf(J,K,N1) = taf(J1,K,N1) 
!#  endif
           END DO
         DO K=1,KBM1
           WQMMAX = MAXVAL(BIO_ALL(NBSN(J,1:NTSN(J)),K,N1))
           WQMMIN = MINVAL(BIO_ALL(NBSN(J,1:NTSN(J)),K,N1))
         
           IF(K == 1)THEN
            WQMMAX = MAX(WQMMAX,(BIO_ALL(J,K,N1)*DZ(J,K+1)+BIO_ALL(J,K+1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K+1)))
            WQMMIN = MIN(WQMMIN,(BIO_ALL(J,K,N1)*DZ(J,K+1)+BIO_ALL(J,K+1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K+1)))
           ELSE IF(K == KBM1)THEN
            WQMMAX = MAX(WQMMAX,(BIO_ALL(J,K,N1)*DZ(J,K-1)+BIO_ALL(J,K-1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K-1)))
            WQMMIN = MIN(WQMMIN,(BIO_ALL(J,K,N1)*DZ(J,K-1)+BIO_ALL(J,K-1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K-1)))
           ELSE
            WQMMAX = MAX(WQMMAX,(BIO_ALL(J,K,N1)*DZ(J,K-1)+BIO_ALL(J,K-1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K-1)), &
                    (BIO_ALL(J,K,N1)*DZ(J,K+1)+BIO_ALL(J,K+1,N1)*DZ(J,K))/  &
		     (DZ(J,K)+DZ(J,K+1)))
            WQMMIN = MIN(WQMMIN,(BIO_ALL(J,K,N1)*DZ(J,K-1)+BIO_ALL(J,K-1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K-1)), &
                    (BIO_ALL(J,K,N1)*DZ(J,K+1)+BIO_ALL(J,K+1,N1)*DZ(J,K))/  &
		     (DZ(J,K)+DZ(J,K+1)))
           END IF
 
           IF(WQMMIN-BIO_F(J,K,N1) > 0.0_SP)BIO_F(J,K,N1) = WQMMIN
           IF(BIO_F(J,K,N1)-WQMMAX > 0.0_SP)BIO_F(J,K,N1) = WQMMAX

         END DO

          ELSE
           DO K=1,KBM1
               BIO_F(J,K,N1)=BIO_ALL(J,K,N1)
!           taf(J,K,N1) = tb(J,K,N1) - ALPHA_OBC*RAMP_TS*(tb(J,K,N1)-TB_OBC(I,N1))
           END DO
         END IF
       END DO
     END DO !!OUTER LOOP OVER BIO-VARIABLES

   END IF

!
!--SET BOUNDARY CONDITIONS-----------------------------------------------------|
!
       BIO_ALL(0,:,:)=ZERO
   RETURN
   END SUBROUTINE BIO_BCOND



SUBROUTINE BIO_NETCDF_HEADER
#if defined (NETCDF_IO)
!==============================================================================!
! THIS PROGRAM OUTPUT BIOLOGICAL VARIABLE, ADOPTED FROM mod_ncdio.F            !
!  NetCDF Output for FVCOM using CF Metadata Convention                        !
!                                                                              !
!    see: http://www.cgd.ucar.edu/cms/eaton/cf-metadata/ for info              !
!                                                                              !
!    current time dependent variables set up                                   !
!         el:    surface elevation                                             !
!          u:    x-velocity. In spherical coordinate,lon-velocity              !
!          v:    y-velocity. In spherical coordinate,lat-velocity              !
!         ww:    z-velocity                                                    !
!         kh:    turbulent diffusivity                                         !
!         km:    turbulent viscosity                                           !
!         t1:    temperature                                                   !
!         s1:    salinity                                                      !
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
!       wd:      wet/dry flag (0 or 1)
!                                                                              !
!    to add additional variables:                                              !
!      1.) add to list above                                                   !
!      2.) add *_vid to variables vid in section "new variable vid"            !
!      3.) go to definition section "new variable definition"                  !
!      4.) add output section "new variable output"                            !
!==============================================================================!
   USE mod_ncdio
!   use netcdf
   use mod_types
   use mod_utils
   use all_vars
   IMPLICIT NONE
   integer, dimension(3) :: dynm3de_lev,dynm3de_lay
   integer, dimension(3) :: dynm3dn_lev,dynm3dn_lay
   integer, dimension(2) :: stat3de_lev,stat3de_lay 
   integer, dimension(2) :: stat3dn_lev,stat3dn_lay 
   integer, dimension(2) :: specdim
   integer, dimension(2) :: dynm2de,dynm2dn
   integer, dimension(1) :: stat2de,stat2dn
   integer, dimension(1) :: stat_lev,stat_lay,dynmtime ,stat_scl
   character(len=100)    :: netcdf_convention
   character(len=100)    :: timestamp ,time_bio
   integer               :: i,j,ierr,i1,i2
   integer               :: maxnode,maxnodep,maxelem,maxelemp,itmp
   real(sp), allocatable :: tmp(:,:),tvec(:)
   character(len=4)      :: nchar
   integer               :: ii,mm
   ALLOCATE (trcsid(ntt))
!--Initialize Stack Count

   stck_cnt = 1
   stck_bio = 0
!--NetCDF Convention String
   netcdf_convention = 'CF-1.0'

!--Time Stamp for History
 !  call get_timestamp(time_bio)
 !  timestamp = 'model started at: '//trim(time_bio)


!==============================================================================|
!  OPEN FILE AND DEFINE VARIABLES                                              |
!==============================================================================|
   IF(MSR)THEN

!--Define NetCDF Output Filename 
     cdfname_bio = trim(outdir)//"/netcdf/"//'bio_netcdf.nc'
   ierr = nf90_create(path=cdfname_bio,cmode=nf90_clobber,ncid=bio_ofid)
!--Description of File Contents
   ierr = nf90_put_att(bio_ofid,nf90_global,"title"      ,trim(casetitle))
   ierr = nf90_put_att(bio_ofid,nf90_global,"institution",trim(institution))
   ierr = nf90_put_att(bio_ofid,nf90_global,"source"     ,trim(fvcom_version))
   ierr = nf90_put_att(bio_ofid,nf90_global,"history"    ,trim(timestamp))
   ierr = nf90_put_att(bio_ofid,nf90_global,"references" ,trim(fvcom_website))
   ierr = nf90_put_att(bio_ofid,nf90_global,"Conventions",trim(netcdf_convention))

!--Define Fixed Model Dimensions 
   ierr = nf90_def_dim(bio_ofid,"scalar" ,1      ,scl_didb    )        
   ierr = nf90_def_dim(bio_ofid,"node"   ,NNodeGL    ,node_didb   )        
   ierr = nf90_def_dim(bio_ofid,"nele"   ,MElemGL    ,nele_didb   )
   ierr = nf90_def_dim(bio_ofid,"siglay" ,kbm1   ,siglay_didb )
   ierr = nf90_def_dim(bio_ofid,"siglev" ,kb     ,siglev_didb )
   ierr = nf90_def_dim(bio_ofid,"three"  ,3      ,three_didb  )
   ierr = nf90_def_dim(bio_ofid,"four"   ,4      ,four_didb   )

!--Define Unlimited Model Dimension
   ierr = nf90_def_dim(bio_ofid,"time"   ,nf90_unlimited,time_didb)

!--Set Up Data Dimensioning - Static Vars
   stat_scl     = (/scl_didb/)             !!scalar variable               
   stat_lay     = (/siglay_didb/)          !!vertical variables at layers
   stat_lev     = (/siglev_didb/)          !!vertical variables at levels
   stat2de      = (/nele_didb/)            !!2d element vars
   stat2dn      = (/node_didb/)            !!2d nodal vars
   stat3de_lay  = (/nele_didb,siglay_didb/) !!3d element vars at layers
   stat3de_lev  = (/nele_didb,siglev_didb/) !!3d element vars at levels
   stat3dn_lay  = (/node_didb,siglay_didb/) !!3d node    vars at layers
   stat3dn_lev  = (/node_didb,siglev_didb/) !!3d node    vars at levels

!--Set Up Data Dimensioning - Dynamic Vars 
   dynm2de      = (/nele_didb,time_didb/)            !!2d element vars
   dynm2dn      = (/node_didb,time_didb/)            !!2d nodal vars
   dynm3de_lay  = (/nele_didb,siglay_didb,time_didb/) !!3d elem vars at layers
   dynm3de_lev  = (/nele_didb,siglev_didb,time_didb/) !!3d elem vars at levels
   dynm3dn_lay  = (/node_didb,siglay_didb,time_didb/) !!3d node vars at layers
   dynm3dn_lev  = (/node_didb,siglev_didb,time_didb/) !!3d node vars at levels
   dynmtime     = (/time_didb/)   

!--Define Coordinate Variables and Attributes

   !!====NPROCS: Number of Processors=======================!
   ierr = nf90_def_var(bio_ofid,"nprocs",nf90_int,stat_scl,nprocs_vidb)
   ierr = nf90_put_att(bio_ofid,nprocs_vidb,"long_name","number of processors")

   !!====PARTITION: Partion Number of Element===============!
   ierr = nf90_def_var(bio_ofid,"partition",nf90_int,stat2de,partition_vidb)
   ierr = nf90_put_att(bio_ofid,partition_vidb,"long_name","partition")

   !!====Initial Density (Used for Constructing 3D Domain)==!
!   ierr = nf90_def_var(bio_ofid,"Initial_Density",nf90_float,stat3dn_lay,idens_vidb)
!   ierr = nf90_put_att(bio_ofid,idens_vidb,"long_name","Initial Density")

   !!====X Grid Coordinate at Nodes (VX) (Meters)===========!
   ierr = nf90_def_var(bio_ofid,"x",nf90_float,stat2dn,x_vidb)
   ierr = nf90_put_att(bio_ofid,x_vidb,"long_name","nodal x-coordinate")
   ierr = nf90_put_att(bio_ofid,x_vidb,"units","meters")

   !!====Y Grid Coordinate at Nodes (VY) (Meters)===========!
   ierr = nf90_def_var(bio_ofid,"y",nf90_float,stat2dn,y_vidb)
   ierr = nf90_put_att(bio_ofid,y_vidb,"long_name","nodal y-coordinate")
   ierr = nf90_put_att(bio_ofid,y_vidb,"units","meters")

   !!====Longitudinal Coordinate at Nodes (LON) (degrees)===!
   ierr = nf90_def_var(bio_ofid,"lon",nf90_float,stat2dn,lon_vidb)
   ierr = nf90_put_att(bio_ofid,lon_vidb,"long_name","Longitude")
   ierr = nf90_put_att(bio_ofid,lon_vidb,"standard_name","longitude")
   ierr = nf90_put_att(bio_ofid,lon_vidb,"units","degrees_east")

   !!====Latitudinal  Coordinate at Nodes (LAT) (degrees)===!
   ierr = nf90_def_var(bio_ofid,"lat",nf90_float,stat2dn,lat_vidb)
   ierr = nf90_put_att(bio_ofid,lat_vidb,"long_name","Latitude")
   ierr = nf90_put_att(bio_ofid,lat_vidb,"standard_name","latitude")
   ierr = nf90_put_att(bio_ofid,lat_vidb,"units","degrees_north")

   !!====Sigma Coordinate for Sigma Layers (ZZ)  (-)========!
   ierr = nf90_def_var(bio_ofid,"siglay",nf90_float,stat3dn_lay,siglay_vidb)
   ierr = nf90_put_att(bio_ofid,siglay_vidb,"long_name","Sigma Layers")
   ierr = nf90_put_att(bio_ofid,siglay_vidb,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(bio_ofid,siglay_vidb,"positive","up")
   ierr = nf90_put_att(bio_ofid,siglay_vidb,"valid_min","-1")
   ierr = nf90_put_att(bio_ofid,siglay_vidb,"valid_max","0")
   ierr = nf90_put_att(bio_ofid,siglay_vidb,"formula_terms","siglay:siglay eta:zeta depth:depth")

   !!====Shifted Sigma Layer Coordinate for Viz ============!
   ierr = nf90_def_var(bio_ofid,"siglay_shift",nf90_float,stat3dn_lay,siglay_shift_vidb)
   ierr = nf90_put_att(bio_ofid,siglay_shift_vidb,"long_name","Shifted Sigma Layers")

   !!====Sigma Coordinate for Sigma Levels (Z)   (-)========!
   ierr = nf90_def_var(bio_ofid,"siglev",nf90_float,stat3dn_lev,siglev_vidb)
   ierr = nf90_put_att(bio_ofid,siglev_vidb,"long_name","Sigma Levels")
   ierr = nf90_put_att(bio_ofid,siglev_vidb,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(bio_ofid,siglev_vidb,"positive","up")
   ierr = nf90_put_att(bio_ofid,siglev_vidb,"valid_min","-1")
   ierr = nf90_put_att(bio_ofid,siglev_vidb,"valid_max","0")
   ierr = nf90_put_att(bio_ofid,siglev_vidb,"formula_terms","siglev:siglev eta:zeta depth:depth")



!--Define Mesh Relevant Variables and Attributes

   !!====Bathymetry at Nodes (H) (meters)===================!
   ierr = nf90_def_var(bio_ofid,"h",nf90_float,stat2dn,h_vidb)
   ierr = nf90_put_att(bio_ofid,h_vidb,"long_name","Bathymetry")   
   ierr = nf90_put_att(bio_ofid,h_vidb,"units","meters")
   ierr = nf90_put_att(bio_ofid,h_vidb,"positive","down")
   ierr = nf90_put_att(bio_ofid,h_vidb,"standard_name","depth")
   ierr = nf90_put_att(bio_ofid,h_vidb,"grid","fvcom_grid")

   !!====Nodes surrounding each Element (NV)================!
   specdim = (/nele_didb,three_didb/) 
   ierr = nf90_def_var(bio_ofid,"nv",nf90_float,specdim,nv_vidb)
   ierr = nf90_put_att(bio_ofid,nv_vidb,"long_name","nodes surrounding element")     

   !!====Momentum Stencil Interpolation Coefficients========!
   specdim = (/nele_didb,four_didb/) 
   ierr = nf90_def_var(bio_ofid,"a1u",nf90_float,specdim,a1u_vidb)
   ierr = nf90_put_att(bio_ofid,a1u_vidb,"long_name","a1u")
   ierr = nf90_def_var(bio_ofid,"a2u",nf90_float,specdim,a2u_vidb)
   ierr = nf90_put_att(bio_ofid,a2u_vidb,"long_name","a2u")

   !!====Element Based Interpolation Coefficients===========!
   specdim = (/nele_didb,three_didb/) 
   ierr = nf90_def_var(bio_ofid,"aw0",nf90_float,specdim,aw0_vidb)
   ierr = nf90_put_att(bio_ofid,aw0_vidb,"long_name","aw0")
   ierr = nf90_def_var(bio_ofid,"awx",nf90_float,specdim,awx_vidb)
   ierr = nf90_put_att(bio_ofid,awx_vidb,"long_name","awx")
   ierr = nf90_def_var(bio_ofid,"awy",nf90_float,specdim,awy_vidb)
   ierr = nf90_put_att(bio_ofid,awy_vidb,"long_name","awy")

!--Define Model Time Variables and Attributes    
   ierr = nf90_def_var(bio_ofid,"time",nf90_float,dynmtime,time_vidb)
   ierr = nf90_put_att(bio_ofid,time_vidb,"long_name","Time")
   ierr = nf90_put_att(bio_ofid,time_vidb,"units",trim(netcdf_timestring))
   ierr = nf90_put_att(bio_ofid,time_vidb,"calendar","none")
   ierr = nf90_def_var(bio_ofid,"iint",nf90_float,dynmtime,iint_vidb)
   ierr = nf90_put_att(bio_ofid,iint_vidb,"long_name","internal mode iteration number")

!--Define Time Dependent Flow Variables (selected by user from input file)

      do mm=1,ntt
          fldnam=TRIM(BIO_NAME(mm,1))
     ierr = nf90_def_var(bio_ofid,fldnam,nf90_float,dynm3dn_lay,ii)	  
!          IERR = NF90_DEF_VAR(bio_ofid,fldnam,NF90_FLOAT,DIMS3DN,ii)
        trcsid(mm)=ii
     ierr = nf90_put_att(bio_ofid,ii,"unit",TRIM(BIO_NAME(mm,2)))
     ierr = nf90_put_att(bio_ofid,ii,"long_name",TRIM(BIO_NAME(mm,3)))
!     ierr = nf90_put_att(bio_ofid,ii,"grid","fvcom_grid")
!     ierr = nf90_put_att(bio_ofid,ii,"type","data")  
      end do

!--Exit Define Mode
   ierr = nf90_enddef(bio_ofid)
   ierr = nf90_close(bio_ofid)

   END IF !(MSR)

!==============================================================================|
!  WRITE VARIABLES TO FILE                                                     |
!==============================================================================|

   IF (MSR) THEN
     ierr = nf90_open(cdfname_bio,nf90_write,bio_ofid)
!     if(ierr /= nf90_noerr)then
!       call handle_ncerr(ierr,"file open error",ipt)
!     end if
   END IF

   !!====Longitude at Nodes (LON) ==========================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vx+vxmin,bio_ofid,lon_vidb,myid,nprocs,ipt)
!   ierr = nf90_put_var(bio_ofid,lon_vidb,vx(i1:i2))
!   ierr = nf90_put_var(bio_ofid,lon_vidb,vx,START=dims)
!   ierr = nf90_put_var(nc_fid,vid,temp,START=dims)
   !!====Latitude  at Nodes (LAT) ==========================!

   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vy+vymin,bio_ofid,lat_vidb,myid,nprocs,ipt)
!   ierr = nf90_put_var(bio_ofid,lon_vidb,vy(i1:i2))
   !!====Number of Processors (NPROCS) =====================!
   if(msr)then 
   ierr = nf90_put_var(bio_ofid,nprocs_vidb,nprocs)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing nprocs variable to netcdf",ipt)
   end if
#  if defined (MULTIPROCESSOR)
!   ierr = nf90_put_var(bio_ofid,partition_vidb,el_pidb)
!   if(ierr /= nf90_noerr)then
!     call handle_ncerr(ierr,"error writing el_pid variable to netcdf",ipt)
!   end if
#  endif
   end if

   !!====Initial Density Field==============================!
!   i1 = lbound(rho1,1) ; i2 = ubound(rho1,1)
!   call putvar(i1,i2,NNode,NNodeGL,kb,kb-1,"n",rho1,bio_ofid,idens_vidb,myid,nprocs,ipt)


   !!====X Grid Coordinate at Nodes (VX)====================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vx+vxmin,bio_ofid,x_vidb,myid,nprocs,ipt)

   !!====Y Grid Coordinate at Nodes (VY)====================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vy+vymin,bio_ofid,y_vidb,myid,nprocs,ipt)

   !!====Bathymetry at Nodes (H)============================!
   i1 = lbound(h,1) ; i2 = ubound(h,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",h,bio_ofid,h_vidb,myid,nprocs,ipt)

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
   call putvar(i1,i2,MElem,MElemGL,3,3,"e",tmp,bio_ofid,nv_vidb,myid,nprocs,ipt)
   deallocate(tmp)

   !!====Momentum Stencil Interpolation Coefficients========!
   i1 = lbound(a1u,1) ; i2 = ubound(a1u,1)
   call putvar(i1,i2,MElem,MElemGL,4,4,"e",a1u,bio_ofid,a1u_vidb,myid,nprocs,ipt)
   i1 = lbound(a2u,1) ; i2 = ubound(a2u,1)
   call putvar(i1,i2,MElem,MElemGL,4,4,"e",a2u,bio_ofid,a2u_vidb,myid,nprocs,ipt)

   !!====Element Based Interpolation Coefficients===========!
   i1 = lbound(aw0,1) ; i2 = ubound(aw0,1)
   call putvar(i1,i2,MElem,MElemGL,3,3,"e",aw0,bio_ofid,aw0_vidb,myid,nprocs,ipt)
   i1 = lbound(awx,1) ; i2 = ubound(awx,1)
   call putvar(i1,i2,MElem,MElemGL,3,3,"e",awx,bio_ofid,awx_vidb,myid,nprocs,ipt)
   i1 = lbound(awy,1) ; i2 = ubound(awy,1)
   call putvar(i1,i2,MElem,MElemGL,3,3,"e",awy,bio_ofid,awy_vidb,myid,nprocs,ipt)

   IF(MSR)THEN
   !!====Sigma Layers (zz)==================================!
   i1 = lbound(zz,1) ; i2 = ubound(zz,1)
   call putvar(i1,i2,NNode,NNodeGL,kb-1,kb-1,"n",zz,bio_ofid,siglay_vid,myid,nprocs,ipt)

!   if(msr)then 
!  allocate(tvec(kbm1))
!   tvec(1:kbm1) = zz(1:kbm1)
!  ierr = nf90_put_var(bio_ofid,siglay_vidb,tvec)
!   if(ierr /= nf90_noerr)then
!     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
!   end if
!   deallocate(tvec)

   allocate(tmp(0:NTNode,kbm1))
   tmp(:,1:kbm1) = z(:,2:kb)
   i1 = lbound(tmp,1) ; i2 = ubound(tmp,1)
   call putvar(i1,i2,NNode,NNodeGL,kb-1,kb-1,"n",tmp,bio_ofid,siglay_shift_vid,myid,nprocs,ipt)
   deallocate(tmp)

!   allocate(tvec(kbm1))
!   tvec(1:kbm1) = z(2:kb)
!   ierr = nf90_put_var(bio_ofid,siglay_shift_vidb,tvec)
!   if(ierr /= nf90_noerr)then
!     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
!   end if
!   deallocate(tvec)

   i1 = lbound(z,1) ; i2 = ubound(z,1)
   call putvar(i1,i2,NNode,NNodeGL,kb,kb,"n",z,bio_ofid,siglev_vid,myid,nprocs,ipt)

!   allocate(tvec(kb))
!   tvec(1:kb) = z(1:kb)
!   ierr = nf90_put_var(bio_ofid,siglev_vidb,tvec)
!   if(ierr /= nf90_noerr)then
!     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
!   end if
!   deallocate(tvec)
   ierr = nf90_close(bio_ofid)  
   END IF !msr?

   RETURN
#endif
END SUBROUTINE BIO_NETCDF_HEADER
!================================================

   subroutine BIO_OUT_NETCDF 
#if defined (NETCDF_IO)
!==============================================================================|
!   Write Time Dependent NetCDF Data to File                                   |
!==============================================================================|

   use all_vars
   use netcdf
#  if defined (WET_DRY)
   use mod_wd
#  endif
   implicit none
   integer :: i,ierr,i1,i2,k,icheck
   integer :: dims(1)
   real*4, allocatable :: ftemp(:)

!==============================================================================|
   
!--Update Counter
!   out_cnt = out_cnt + 1
!   stck_cnt = stck_cnt + 1 
    stck_bio = stck_bio + 1
!--Open File
   if(msr)then
     ierr = nf90_open(cdfname_bio,nf90_write,bio_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error for bio_netcdf",ipt)
     end if
   end if
!--Dump Time/IINT to File
   dims(1) = stck_cnt
!   dims(1) = stck_bio
   if(msr)then
   ierr    = nf90_put_var(bio_ofid,iint_vidb,float(iint),START=dims)
    if(ierr /= nf90_noerr) print*,'error writing time steps to netcdf'
   ierr    = nf90_put_var(bio_ofid,time_vidb,thour*3600.,START=dims)
    if(ierr /= nf90_noerr) print*,'error writing time to netcdf'
   end if
!--Write Variables to File
!   if(msr) write(ipt,*)'dumping to netcdf file: ',trim(cdfname),stck_cnt

!      IF (MSR) THEN
      do ii=1,ntt
      I1 = LBOUND(BIO_ALL,1) ; I2 = UBOUND(BIO_ALL,1)
      mm=trcsid(ii)
 CALL putvar_bio(i1,i2,NNode,NNodeGL,kb,kb-1,"n",BIO_ALL(I1:I2,:,ii),bio_ofid,mm,MYID,NPROCS,IPT)
      end do
!      END IF !MSR
!==============================================================================|
!  CLOSE THE FILE                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(bio_ofid)
#endif
!end if defined netcdf-io
   return
   end subroutine BIO_OUT_NETCDF


!==============================================================================|
!  Collect Data to Global Array and Write to Netcdf File                       |
!==============================================================================|
                                                                                                                  
   subroutine putvar_bio(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid,myid,nprocs,ipt)

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

#if defined (NETCDF_IO)
   k1m1 = k1 
   if(k1m1 == 1)then
     allocate(dims(2))
     dims(1) = 1 
!     dims(2) = stck_cnt
    dims(2) = stck_bio
   else
     allocate(dims(3))
     dims(1) = 1 
     dims(2) = 1 
!     dims(3) = stck_cnt
     dims(3) = stck_bio
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
!   ierr = nf90_put_var(bio_ofid,vid,temp,START=dims)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(dims)
   return
#endif
!end if defined netcdf-io
   end subroutine putvar_bio


SUBROUTINE BIO_EXCHANGE
!==============================================================================!
!     PERFORM DATA EXCHANGE FOR the Generalized biological model               |
!==============================================================================!
#if defined (MULTIPROCESSOR)
!     USE ALL_VARS
     USE MOD_PAR
     USE LIMS
     USE CONTROL
     IMPLICIT NONE
     INTEGER :: I3
     DO I3=1,NTT
      IF(PAR) CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,BIO_ALL(:,:,I3))
      IF(PAR) CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,BIO_MEAN(:,:,I3))
      IF(PAR) CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,BIO_F(:,:,I3))
     END DO 
   RETURN
#endif
  END SUBROUTINE BIO_EXCHANGE

   SUBROUTINE BIO_INITIAL
!=============================================================================!
! THS PROGRAM INITIALIZES THE 3D BIOLOGICAL FIELD FOR THE GENERALIZED         !
! BIOLOGICAL MODEL: BIO_ALL(I,K,Nl),N1=1,NTT),AND MEAN VALUES BIO_MEAN        !
! EACH BIOLOGICAL STATE VARIABLE HAS AN INDEPENDENT INITIAL CONDITION FILE    !
! PLACED IN INPDIR. THEY SHOULD BE NAME AS "NUTRIENT_INI_1", "NUTRIENT_INI_2",!
! "PHYTOPLANKTON_INI_1", "ZOOPLANKTON_INI_1", "BACTERIA_INI_1", 'DETRITUS_    !
! INI_1", "DOM_INI_1" AND SO FORTH. THREE TYPES OF INITIAL CONDITIONS WERE    !
! IMPLEMENTED: (1) 'CONSTANT': A SINGLE VALUE; (2) 'LINEAR':WITH AT LEAST TWO !
! PAIRS OF VALUES WITH DEPTH. VARIABLE VALUES WILL BE LINEARLY INTERPOLATED   !
! BETWEEN THE VALUES GIVEN), (3) "DATA": OBSERVATION DATA SHOULD BE INTER-    !
! POLATED ONTO THE GRID POINTS AT STANDARD LEVELS. VARIABLE VALUES WILL BE    !
! INTERPOLATED AT EACH GRID POINT FROM THE DATA. THE TYPE OF INITIAL CONDI-   !
! TION SHOULD BE PUT ON THE FIRST LINE OF EACH INITIAL FILE                   !
!=============================================================================!
      IMPLICIT NONE
      INTEGER :: I,J,K,LL,N_DATA
      CHARACTER(LEN=80) :: ISTR
      CHARACTER(LEN=1)  :: BIO_NUMBER
      CHARACTER(LEN=10) :: INI_TYPE
      REAL(SP), DIMENSION(KBM1)    :: ZM           !GRID DEPTH
      REAL(SP), DIMENSION(500)     :: DEPTH_STD    !STANDARD DEPTH OF DATA
      REAL(SP), DIMENSION(500)     :: DATA_BIO     !STANDARD DATA FOR LINEAR INTERPOLATION
      REAL(SP), DIMENSION(KB)      :: DATA_INT     !INTERPOLDATED VALUES
      REAL(SP), DIMENSION(NNodeGL,KSL) :: TEMPB        !TEMPERAL FOR DATA INPUT
      REAL(SP), DIMENSION(NNode,KSL)   :: DATA_3D      !3D OBSERVATION DATA
      ALLOCATE(BIO_ALL(0:NTNode,KB,NTT))    ; BIO_ALL     =  0.001_SP
      ALLOCATE(BIO_F(0:NTNode,KB,NTT))      ; BIO_F       =  0.001_SP
      ALLOCATE(BIO_MEAN(0:NTNode,KB,NTT))   ; BIO_MEAN    =  0.001_SP
      ALLOCATE(XFLUX_OBCB(0:NTNode,KB,NTT)) ; XFLUX_OBCB  =  0.0_SP
      ALLOCATE(BIO_MEANN(0:MTElem,KB,NTT))  ; BIO_MEANN   =  0.001_SP
      ALLOCATE(BIO_VAR_MEAN(0:NTNode,KB,NTT)) ; BIO_VAR_MEAN   =  0.0_SP


!*******************************************************************!
! *********     PRINT OUT MODEL SETUP AND PARAMETER VALUES      ****!
!*******************************************************************!
       IF (MSR) THEN
       PRINT*
       PRINT*,'*****************************************************'
       PRINT*,'**  STRUCTURE AND FUNCTION OF THE BIOLOGICAL MODEL **'
       PRINT*,'*****************************************************'

          PRINT*
          PRINT*,'MODEL STRUCTURE        : ',  MODEL_STRUCTURE
         DO I=1,NNN
          PRINT*,'                         ',  NUTRIENT_NAME(I)
         END DO
         DO I=1,NNP
          PRINT*,'                         ',  PHYTO_NAME(I)
         END DO
         DO I=1,NNZ
          PRINT*,'                         ',  ZOO_NAME(I)
         END DO
         DO I=1,NND
          PRINT*,'                         ',  DETRITUS_NAME(I)
         END DO
         DO I=1,NNM
          PRINT*,'                         ', DOM_NAME(I)
         END DO
         DO I=1,NNB
          PRINT*,'                         ', BACTERIA_NAME(I)
         END DO
          WRITE(*,'(A26,A20)')' LIGHT FUNCTION         : ', L_FUNCTION
          WRITE(*,'(A26,A20)')' GRAZING FUNCTION       : ', G_FUNCTION
          PRINT*
      PRINT*,'*********    PHYTOPLANTON PARAMETERS    **************'
        PRINT*
        IF(L_FUNCTION.NE.'EXP_LIGHT'.AND.L_FUNCTION.NE.'SL62_LIGHT') THEN
          PRINT*,'ALPHA                  : '  , (ALPHAP(I),I=1,NNP)
        END IF
          PRINT*,'L_N COMBINE            : ' , (ALPHA_U(I),I=1,NNP)
          PRINT*,'T FORCING              : ' , (A_TP(I),I=1,NNP)
          PRINT*,'CHL ATTANUATION        : ' , ATANU_C
          PRINT*,'D ATTANUATION          : ', ATANU_D
          PRINT*,'WATER ATTANUATION      : ', ATANU_W
        IF(L_FUNCTION.EQ.'PGH80_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT'.OR. &
           L_FUNCTION.EQ.'BWDC9_LIGHT') THEN       
          PRINT*,'BETAP                  : ', (BETAP(I),I=1,NNP)
        END IF
          PRINT*,'CHL:C                  : ',(CHL2C(I),I=1,NNP)
          PRINT*,'ACTIVE DOM EXUD.       : ',(D_DOM(I),I=1,NNP)
          PRINT*,'PASSIVE DOM EXUD.      : ', (DPDOM(I),I=1,NNP)
        IF(L_FUNCTION.EQ.'SL62_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT'.OR. &
           L_FUNCTION.EQ.'PE78_LIGHT') THEN       
          PRINT*,'OPTIMAL LIGHT          : ',(I_OPT(I),I=1,NNP)
        END IF
        IF(L_FUNCTION.EQ.'MM_LIGHT'.OR.L_FUNCTION.EQ.'LB_LIGHT') THEN 
          PRINT*,'LIGHT HALF SATURATTION : ', (K_LIGHT(I),I=1,NNP)
        END IF
          PRINT*,'MORTALITY              : ', (MPD(I),I=1,NNP)
          PRINT*,'MORTALITY POWER        : ', (M_P(I),I=1,NNP)
        IF(L_FUNCTION.EQ.'LB_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT') THEN 
          PRINT*,'POWER OF LIGHT         : ', (N_P(I),I=1,NNP)
        END IF
          PRINT*,'THRESHOLD              : ', (P_0(I),I=1,NNP)
          PRINT*,'T ON RESPIRATION       : ', RP_T
          PRINT*,'RESPIRATION            : ', (R_P(I),I=1,NNP)
          PRINT*,'OPTIMAL T              : ', (T_OPTP(I),I=1,NNP)
          PRINT*,'MAXIMUM GROWTH         : ', (UMAX(I),I=1,NNP)
          PRINT*,'SINKING VELOCITY       : ', (W_P(I),I=1,NNP)
          PRINT*

          PRINT*,'********    ZOOPLANKTON PARAMETERS   ****************'
          PRINT*
          PRINT*,'ACTIVE RESPIRATION     : ', (ACTIVE_R(I),I=1,NNZ)
          PRINT*,'T FORCING EXPONENTIAL  : ', (A_TZ(I),I=1,NNZ)
        IF (NNB.GE.1) THEN
          PRINT*,'EFFICIENCY ON BACTERIA : ',((EFFIB(I,J),I=1,NNB),J=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'EFFICIENCY ON DETRITUS : ',((EFFID(I,J),I=1,NND),J=1,NNZ)
        END IF
          PRINT*,'EFFICIENCY ON PHYTO    : ',((EFFIP(I,J),I=1,NNP),J=1,NNZ)
          PRINT*,'EFFICIENCY ON ZOO      : ',((EFFIZ(I,J),I=1,NNZ),J=1,NNZ)
          PRINT*,'MAX GRAIZING RATE      : ',(G_MAX(I),I=1,NNZ)
        IF(G_FUNCTION.EQ.'RECTI_G'.OR.G_FUNCTION.EQ.'MM1_G' .AND. &
           G_FUNCTION.EQ.'MM2_G') THEN 
          PRINT*,'HALF SATURATION        : ', (K_ZG(I),I=1,NNZ)
        END IF
          PRINT*,'GRAZING POWER          : ',(M_G(I),I=1,NNZ)
          PRINT*,'MORALITY               : ',(MZD(I),I=1,NNZ)
          PRINT*,'MORTALITY POWER        : ',(M_Z(I),I=1,NNZ)
        IF(G_FUNCTION.EQ.'MM2_G') THEN 
          PRINT*,'GRAZING THRESHOLD      : ',(P_C(I),I=1,NNZ)
        END IF
        IF (NNZ.GE.1) THEN
          PRINT*,'RECRUITMENT            : ',(R_RECRUIT(I),I=1,NNZ)
        END IF
          PRINT*,'RESPIRATION            : ',(R_Z(I),I=1,NNZ)
        IF (NNB.GE.1) THEN
          PRINT*,'PREFERENCE ON BACTERIA : ',((SIGMA_B(I,J),I=1,NNB),J=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'PREFERENCE ON DETRITUS : ',((SIGMA_D(I,J),I=1,NND),J=1,NNZ)
        END IF
          PRINT*,'PREFERENCE ON PHYTO    : ',((SIGMA_P(I,J),I=1,NNP),J=1,NNZ)
          PRINT*,'PREFERENCE ON ZOO      : ',((SIGMA_Z(I,J),I=1,NNZ),J=1,NNZ)
          PRINT*,'OPTIMAL T              : ',(T_OPTZ(I),I=1,NNZ)
          PRINT*,'ZOO THRESHOLD          : ',(Z_0(I),I=1,NNZ)
          PRINT*
          PRINT*,'*********    NUTRIENT PARAMETERS    ***********'
          PRINT*
          PRINT*,'HALF-SATURATION',((KSN(I,J),I=1,NNN),J=1,NNP)
        IF (NNB.GE.1) THEN
          PRINT*,'ELEMENT RATIO IN BAC.  : ',((N2CB(I,J),I=1,NNN),J=1,NNB)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'ELEMENT RATIO IN D     : ',((N2CD(I,J),I=1,NNN),J=1,NND)
        END IF
          PRINT*,'ELEMENT RATIO IN PHYTO : ',((N2CP(I,J),I=1,NNN),J=1,NNP)
          PRINT*,'ELEMENT RATIO IN ZOO   : ',((N2CZ(I,J),I=1,NNN),J=1,NNZ)
        IF (NNM.GE.1)  THEN
          PRINT*,'ELEMENT RATIO IN DOM   : ', ((N2CDOM(I,J),I=1,NNN),J=1,NNM)
        END IF
          PRINT*,'THRESHOLD              : ', (N_0(I),I=1,NNN)
        IF (NO3_ON)    THEN
          PRINT*,'NITRIFICATION RATE     : ', R_AN
        END IF
          PRINT*
        IF (NNB.GE.1) THEN
          PRINT*,'*********  BACTERIA PARAMETERS  ************'
          PRINT*
          PRINT*,'T_FORCING              : ',(A_TB(I),I=1,NNB)
          PRINT*,'THRESHOLD              : ', (B_0(I),I=1,NNB)
          PRINT*,'RATIO OF NH4 VS DON    : ',(DELTA_B(I),I=1,NNB)
          PRINT*,'EFFICIENCY OF DETRITUS : ',((EFFIBD(I,J),I=1,NND),J=1,NNB)
          PRINT*,'EFFICIENCY OF DOM      : ',((EFFIDOM(I,J),I=1,NNM),J=1,NNB)
          PRINT*,'EFFICIENCY OF NUTRIENT : ', ((EFFIN(I,J),I=1,NNN),J=1,NNB)
          PRINT*,'RESPIRATION            : ',(R_B(I),I=1,NNB)
          PRINT*,'PREFERENCE ON DETRITUS : ',((SIGMA_BD(I,J),I=1,NND),J=1,NNB)
          PRINT*,'PREFERENCE ON DOM      : ',((SIGMA_DOM(I,J),I=1,NNM),J=1,NNB)
          PRINT*,'PREFERENCE ON NUTRIENT : ',((SIGMA_N(I,J),I=1,NNN),J=1,NNB)
          PRINT*,'OPTIMAL TEMPERATURE    : ',(T_OPTB(I),I=1,NNB)
          PRINT*,'MAXIMUM GROWTH RATE    : ',(UBMAX(I),I=1,NNB)
          PRINT*
        END IF
        IF (NND.GE.1) THEN
          PRINT*, '********  DETRITUS PARAMETERS  *************'
          PRINT*
        IF (NNB.GE.1) THEN
          PRINT*,'GRAZING ON B TO D      : ', (((ALPHA_BD(I,J,K),I=1,NND),J=1,NNB),K=1,NNZ)
        END IF
          PRINT*,'AGGREGATION            : ', (ALPHA_DAG(I),I=1,NND)
          PRINT*,'GRAZING LOSS ON D TO D : ', (((ALPHA_DD(I,J,K),I=1,NND),J=1,NND),K=1,NNZ)
          PRINT*,'DISAGGREGATION         : ', (ALPHA_DAG(I),I=1,NND)
          PRINT*,'GRAZING LOSS ON P TO D : ', (((ALPHA_PD(I,J,K),I=1,NND),J=1,NNP),K=1,NNZ)
          PRINT*,'GRAZING LOSS ON Z TO D : ', (((ALPHA_ZD(I,J,K),I=1,NND),J=1,NNZ),K=1,NNZ)
        IF (NNM.GE.1) THEN
          PRINT*,'DISSOLUTION            : ', (D_DOM(I),I=1,NND)
        END IF
          PRINT*,'REMINERALIZATION       : ', (D_RN(I),I=1,NND)
          PRINT*,'THRESHOLD              : ', (D_0(I),I=1,NND)
          PRINT*,'P MORTALITY TO DETRITUS: ', ((EPSILON_PD(I,J),I=1,NND),J=1,NNP)
          PRINT*,'Z MORTALITY TO DETRITUS: ', ((EPSILON_ZD(I,J),I=1,NND),J=1,NNZ)
          PRINT*,'SINING VELOCITY        : ', (W_D(I),I=1,NND)
          PRINT*
        END IF
        IF (NNM.GE.1) THEN
          PRINT*,'*********     DOM PARAMETERS     *************'
          PRINT*
          PRINT*,'DOM AGEING COEFFICIENT : ', (ALPHA_DOM(I),I=1,NNM)
          PRINT*,'PHYTO EXUDATION        : ', ((ALPHA_PDOM(I,J),I=1,NNM),J=1,NNP)
          PRINT*,'DETRITUS DISSOLUTION   : ', ((ALPHA_DDOM(I,J),I=1,NNM),J=1,NND)
        IF (NNB.GE.1) THEN
          PRINT*,'GRAZING LOSS ON B > DOM: ', (((ALPHA_ZBDOM(I,J,K),I=1,NNM),J=1,NNB),K=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'GRAZING LOSS ON D > DOM: ', (((ALPHA_ZDDOM(I,J,K),I=1,NNM),J=1,NND),K=1,NNZ)
        END IF
          PRINT*,'GRAZING LOSS ON P > DOM: ', (((ALPHA_ZPDOM(I,J,K),I=1,NNM),J=1,NNP),K=1,NNZ)
          PRINT*,'GRAZING LOSS ON Z > DOM: ', (((ALPHA_ZZDOM(I,J,K),I=1,NNM),J=1,NNZ),K=1,NNZ)
          PRINT*,'THRESHOLD              : ', (DOM_0(I),I=1,NNM)
          PRINT*
        END IF
       PRINT*
       PRINT*,'*****************************************************'
       PRINT*,'*********    END OF BIOLOGICAL MODEL     ************'
       PRINT*,'*****************************************************'
       END IF !(MSR)


!*********      NUTRIENT INITIAL CONDITIONS   *************
      DO LL=1,NNN
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE='INPDIR/NUTRIENT_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,NNode
            DO K=1,KB
              BIO_ALL(I,K,LL+INN-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 11       READ(1,*,END=12)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 11
 12       N_DATA=N_DATA-1
          DO I=1,NNode
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*(D(I)+EL(I))
            END DO
            CALL SINTER(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INN-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR
        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
          DO I=1,NNodeGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,NNode
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,NNode
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INN-1) = DATA_INT(K)
            END DO
          END DO
        END IF !DATA
      END DO !L=1,NNN; NUTRIENT INITIALIZATION

!*********      PHYTOPLANKTON INITIAL CONDITIONS   *************
      DO LL=1,NNP
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE='INPDIR/PHYTOPLANKTON_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,NNode
            DO K=1,KB
              BIO_ALL(I,K,LL+INP-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 21       READ(1,*,END=22)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 21
 22       N_DATA=N_DATA-1
          DO I=1,NNode
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INP-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
          DO I=1,NNodeGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,NNode
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,NNode
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INP-1) = DATA_INT(K)
            END DO
          END DO
        END IF !DATA
      END DO !L=1,NNP; PHYTOPLANKTON INITIALIZATION

!*********      ZOOPLANKTON INITIAL CONDITIONS   *************
      DO LL=1,NNZ
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE='INPDIR/ZOOPLANKTON_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,NNode
            DO K=1,KB
              BIO_ALL(I,K,LL+INZ-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 31       READ(1,*,END=32)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 31
 32       N_DATA=N_DATA-1
          DO I=1,NNode
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INZ-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
          DO I=1,NNodeGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,NNode
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,NNode
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INZ-1) = DATA_INT(K)
            END DO
          END DO
        END IF !DATA
      END DO !L=1,NNZ; ZOOPLANKTON INITIALIZATION


!*********      BACTERIA INITIAL CONDITIONS   *************
      DO LL=1,NNB
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE='INPDIR/BACTERIA_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,NNode
            DO K=1,KB
              BIO_ALL(I,K,LL+INB-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 41       READ(1,*,END=42)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 41
 42       N_DATA=N_DATA-1
          DO I=1,NNode
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INB-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
          DO I=1,NNodeGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,NNode
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,NNode
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INB-1) = DATA_INT(K)
            END DO
          END DO
        END IF !DATA
      END DO !L=1,NNB; BACTERIA INITIALIZATION


!*********      DOM INITIAL CONDITIONS   *************
      DO LL=1,NNM
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE='INPDIR/DOM_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,NNode
            DO K=1,KB
              BIO_ALL(I,K,LL+INM-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 51       READ(1,*,END=52)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 51
 52       N_DATA=N_DATA-1
          DO I=1,NNode
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INM-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
          DO I=1,NNodeGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,NNode
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,NNode
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INM-1) = DATA_INT(K)
            END DO
          END DO
        END IF !DATA
      END DO !L=1,NNM; DOM INITIALIZATION

!*********      DETRITUS INITIAL CONDITIONS   *************
      DO LL=1,NND
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE='INPDIR/DETRITUS_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,NNode
            DO K=1,KB
              BIO_ALL(I,K,LL+IND-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 61       READ(1,*,END=62)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 61
 62       N_DATA=N_DATA-1
          DO I=1,NNode
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+IND-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
          DO I=1,NNodeGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,NNode
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,NNode
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+IND-1) = DATA_INT(K)
            END DO
          END DO
        END IF !DATA
      END DO !L=1,NND; DETRITUS INITIALIZATION
      WHERE (BIO_ALL < 0.001) BIO_ALL=0.001
      BIO_MEAN=BIO_ALL   !3D ASSIGNMENT
#    if defined (MULTIPROCESSOR)
      IF (PAR) CALL BIO_EXCHANGE
#    endif
      CALL BIO_NETCDF_HEADER
   RETURN
END SUBROUTINE BIO_INITIAL


SUBROUTINE BIO_HOT_START
! THIS SUBROUTINE READS IN RESTART BIOLOGICAL DATA FROM THE NETCDF FILE restart_bio.nc  !
! IT QUERIES THE DIMENSION AND TAKES THE LAST TIME STEP AS THE RESTART DATA.            !

#if defined (NETCDF_IO)
      USE NETCDF  
#endif      
      IMPLICIT NONE
      SAVE
      INTEGER  ::  I,J,K,IERR,NC_FID,N_START,VARID,DIMID,DIMS(3)
      REAL(SP) ::  TEMPB(NNodeGL,KBM1)
      CHARACTER(LEN=20) ::  TEMPNAME,time
      ALLOCATE(BIO_ALL(0:NTNode,KB,NTT))    ; BIO_ALL     =  0.001_SP
      ALLOCATE(BIO_F(0:NTNode,KB,NTT))      ; BIO_F       =  0.001_SP
      ALLOCATE(BIO_MEAN(00:NTNode,KB,NTT))  ; BIO_MEAN    =  0.001_SP
      ALLOCATE(XFLUX_OBCB(0:NTNode,KB,NTT)) ; XFLUX_OBCB  =  0.0_SP
      ALLOCATE(BIO_MEANN(0:MTElem,KB,NTT))  ; BIO_MEANN   =  0.001_SP
      ALLOCATE(BIO_VAR_MEAN(0:NTNode,KB,NTT)) ; BIO_VAR_MEAN   =  0.0_SP

!*******************   START EXECUTABLE      *******************!
#if defined (NETCDF_IO)
!      IF(MSR) THEN
        IERR = NF90_OPEN('restart_bio.nc',NF90_NOWRITE,NC_FID)
        IF(IERR /=NF90_NOERR)THEN
           WRITE(*,*)' ERROR IN OPENNING restart_bio.nc '
           STOP
        END IF
!      END IF
      IERR = NF90_INQ_DIMID(NC_FID,'time',DIMID)   !GET time DIMENSION ID
      IF(IERR /=NF90_NOERR)THEN
         WRITE(*,*)' ERROR GETTING time ID: '
         STOP
      END IF
      IERR = NF90_INQUIRE_DIMENSION(NC_FID,DIMID,TEMPNAME,N_START)
      IF(IERR /=NF90_NOERR)THEN
         WRITE(*,*)' ERROR GETTING TIME STEPS '
         STOP
      END IF
!**************   DETERMINE START POINT
      DIMS(1) = 1
      DIMS(2) = 1
      DIMS(3) = N_START
      DO I=1,NTT
         IERR = nf90_inq_varid(NC_FID,TRIM(BIO_NAME(I,1)),VARID)
         IF(IERR /=NF90_NOERR)THEN
           WRITE(*,*)'ERROR GETTING THE ID OF ',TRIM(BIO_NAME(I,1))
           STOP
         END IF

         IERR = nf90_get_var(NC_FID,VARID,TEMPB,START=DIMS)
         IF(IERR /=NF90_NOERR)THEN
            WRITE(*,*)'ERROR GETTING RESTART DATA OF ',TRIM(BIO_NAME(I,1))
            STOP
         END IF

      IF(SERIAL)THEN
        BIO_ALL(:,:,I)=TEMPB(:,:)
      END IF
#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO J=1,NNode
              BIO_ALL(J,1:KBM1,I)=TEMPB(NGID(J),1:KBM1)
            END DO
          END IF
#    endif
      END DO  !DO I=1,NNT
      BIO_MEAN=BIO_ALL
#    if defined (MULTIPROCESSOR)
      IF (PAR) CALL BIO_EXCHANGE
#    endif
      CALL BIO_NETCDF_HEADER
      IERR=NF90_CLOSE(NC_FID)
# endif      
END SUBROUTINE BIO_HOT_START

#    endif
!    END IF DEFINED BioGen AT THE BEGINNING
END MODULE MOD_BIO_3D
