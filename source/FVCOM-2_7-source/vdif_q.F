!==============================================================================|
!   this subroutine is used to calculate the q2 and q2l by solving             !
!   the vertical diffusion equation implicitly.                                !
!==============================================================================|

   SUBROUTINE VDIF_Q              

!------------------------------------------------------------------------------|
   USE ALL_VARS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   use mod_par
   IMPLICIT NONE
   INTEGER :: I,J,K,KI
   REAL(SP)  :: CONST1,COEF1,COEF2,COEF3,COEF4,COEF5,LMAX
   REAL(SP), DIMENSION(0:NTNode,KB) :: GH,SM,SH,PROD,DTEF,KN,BOYGR,A,C,VHP,VH,STF
   REAL(SP), DIMENSION(0:NTNode) :: L0
   REAL(SP)  :: UTAU2
   REAL(SP)  :: WUSURF_NODE,WVSURF_NODE,WUBOT_NODE,WVBOT_NODE
   REAL(SP)  :: UU_NODE_K,VV_NODE_K,UU_NODE_KM1,VV_NODE_KM1

!--Stability Function Coefficients-------------
   REAL(SP), PARAMETER :: A1    = 0.92_SP
   REAL(SP), PARAMETER :: B1    = 16.6_SP
   REAL(SP), PARAMETER :: A2    = 0.74_SP
   REAL(SP), PARAMETER :: B2    = 10.1_SP
   REAL(SP), PARAMETER :: C1    = 0.08_SP

!--Source/Sink Term Coefficients---------------
   REAL(SP), PARAMETER :: E1    = 1.80_SP
   REAL(SP), PARAMETER :: E2    = 1.33_SP
   REAL(SP), PARAMETER :: SEF   = 1.00_SP

!--Von Karmans Constant------------------------
   REAL(SP), PARAMETER :: KAPPA = 0.40_SP

!--Gravitational Constant----------------------
   REAL(SP), PARAMETER :: GEE   = 9.806_SP

!--Limiting Values of Q2/Q2L-------------------
   REAL(SP), PARAMETER :: SMALL = 1.E-8_SP

!--Coefficient of buoyancy length scale--------
   REAL(SP), PARAMETER :: CB    = 1.E-8_SP

!--Denominator small number threshold----------
   REAL(SP), PARAMETER :: DEN_S = 1.E-10_SP

!--Upper bound of GH for unstable flow---------
   REAL(SP), PARAMETER :: GH_MAX = .0233_SP

!--Lower bound for GH for stable flow----------
   REAL(SP), PARAMETER :: GH_MIN = -.281_SP

!----------------------------------------------
   REAL(SP), PARAMETER :: CBCNST = 100.0_SP
   REAL(SP), PARAMETER :: SURFL = 2.E+5_SP
   REAL(SP), PARAMETER :: SHIW = 0.0_SP
   REAL(SP), PARAMETER :: GHC = -6.0_SP

#  if defined (WET_DRY)
   REAL(SP)  :: KQ_TMP,KM_TMP,KH_TMP
#  endif

!------------------------------------------------------------------------------|

!------------------------------------------------------------------------------|
!  place a threshold on the turbulence quantities following advection          |
!------------------------------------------------------------------------------|

   DO K = 2, KBM1
     DO I = 1, NNode
       IF (Q2F(I,K) <= small .OR. Q2LF(I,K) <= small) THEN
!         Q2F(I,K)  = 0.0_SP
!         Q2LF(I,K) = 0.0_SP
         Q2F(I,K)  = small
         Q2LF(I,K) = small
       END IF
     END DO
   END DO

!------------------------------------------------------------------------------!
!  set up coefficients for implicit calculation of vertical diffusion          !
!------------------------------------------------------------------------------!

   DO I=1,NNode
#  if defined (WET_DRY)
       IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
     DO K=2,KBM1
       A(I,K) = -0.5_SP* DTI * (KQ(I,K+1)+KQ(I,K)+2.0_SP*UMOL)/ &
             (DZZ(I,K-1)*DZ(I,K)*D(I)*D(I))
       C(I,K) = -0.5_SP* DTI*(KQ(I,K-1)+KQ(I,K)+2.0_SP*UMOL)/ &
             (DZZ(I,K-1)*DZ(I,K-1)*D(I)*D(I))
     END DO
#  if defined (WET_DRY)
     END IF
#  endif
   END DO

!------------------------------------------------------------------------------!
!    the following section solves the equation                                 !
!    dti*(kq*q2')' - q2*(2.*dti*dtef+1.) = -q2b                                !
!------------------------------------------------------------------------------!

   CONST1 = 16.6_SP ** .6666667_SP * SEF
   DO I = 1, NNode
!       VHP(I,1) = SQRT(WUSURF(I)**2+WVSURF(I)**2) * CONST1
!      VH(I,1)  = 0.0_SP

     WUSURF_NODE = 0.0_SP
     WVSURF_NODE = 0.0_SP
     DO J=1,NTVE(I)
       WUSURF_NODE = WUSURF_NODE + WUSURF(NBVE(I,J))
       WVSURF_NODE = WVSURF_NODE + WVSURF(NBVE(I,J))
     END DO
     WUSURF_NODE = WUSURF_NODE/FLOAT(NTVE(I))
     WVSURF_NODE = WVSURF_NODE/FLOAT(NTVE(I))

     UTAU2 = SQRT(WUSURF_NODE**2+WVSURF_NODE**2)
     VHP(I,1) = (15.8_SP*CBCNST)**0.6666667_SP*UTAU2
     VH(I,1)  = 0.0_SP
!     L0(I)    = SURFL*UTAU2/GEE
     L0(I)    = SURFL*UTAU2/GRAV_N(I)

     WUBOT_NODE = 0.0_SP
     WVBOT_NODE = 0.0_SP
     DO J=1,NTVE(I)
       WUBOT_NODE = WUBOT_NODE + WUBOT(NBVE(I,J))
       WVBOT_NODE = WVBOT_NODE + WVBOT(NBVE(I,J))
     END DO
     WUBOT_NODE = WUBOT_NODE/FLOAT(NTVE(I))
     WVBOT_NODE = WVBOT_NODE/FLOAT(NTVE(I))

     Q2F(I,KB) = SQRT(WUBOT_NODE**2+WVBOT_NODE**2) * CONST1
   END DO

   Q2  = ABS(Q2)
   Q2L = ABS(Q2L)
 
!------------------------------------------------------------------------------!
!  calculate boygr = -Brunt Vaisala Frequency^2                                !
!  calculate internal wave shear energy contribution: prod = 200*(b.v.freq)**3 !
!------------------------------------------------------------------------------!

   DO K=2,KBM1
     DO I=1,NNode
#  if defined (WET_DRY)
      IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
!       BOYGR(I,K) = GEE * (RHO1(I,K-1)-RHO1(I,K))/(DZZ(I,K-1)*D(I))
       BOYGR(I,K) = GRAV_N(I) * (RHO1(I,K-1)-RHO1(I,K))/(DZZ(I,K-1)*D(I))
       PROD(I,K)  = KM(I,K) * 0.0_SP * (.5_SP*(-BOYGR(I,K)+ABS(BOYGR(I,K)))) ** 1.5_SP
#  if defined (WET_DRY)
      END IF
#  endif
     END DO
   END DO

!------------------------------------------------------------------------------!
!  calculate shear production source term = prod                               !
!  calculate buoyancy production source term = kh*boygr                        !
!------------------------------------------------------------------------------!
   DO  K = 2, KBM1
     DO  I = 1, NNode
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif

         UU_NODE_K   = 0.0_SP
         VV_NODE_K   = 0.0_SP
         UU_NODE_KM1 = 0.0_SP
         VV_NODE_KM1 = 0.0_SP
         DO J=1,NTVE(I)
           UU_NODE_K   = UU_NODE_K + U(NBVE(I,J),K)
           VV_NODE_K   = VV_NODE_K + V(NBVE(I,J),K)
           UU_NODE_KM1 = UU_NODE_KM1 + U(NBVE(I,J),K-1)
           VV_NODE_KM1 = VV_NODE_KM1 + V(NBVE(I,J),K-1)
         END DO
         UU_NODE_K   = UU_NODE_K/FLOAT(NTVE(I))
         VV_NODE_K   = VV_NODE_K/FLOAT(NTVE(I))
         UU_NODE_KM1 = UU_NODE_KM1/FLOAT(NTVE(I))
         VV_NODE_KM1 = VV_NODE_KM1/FLOAT(NTVE(I))

         PROD(I,K) = PROD(I,K) + KM(I,K) * SEF * (( &
                     UU_NODE_K-UU_NODE_KM1)**2+(VV_NODE_K-VV_NODE_KM1)**2)/ &
                     (DZZ(I,K-1)*D(I))**2
         PROD(I,K) = PROD(I,K) + KH(I,K) * BOYGR(I,K)
       END IF
     END DO
   END DO
!------------------------------------------------------------------------------!
!  solve for turbulent length scale l = q2l/q2                                 !
!------------------------------------------------------------------------------!

!   L = Q2L/Q2
    DO K=2,KBM1
      DO I=1,NNode
        L(I,K) = Q2L(I,K)/Q2(I,K)
      END DO
    END DO

!------------------------------------------------------------------------------!
!  in stably stratified regions, length scale is limited by buoyancy length    !
!  scale, l_b = c_b*(k^.5/N).  length scale limitation also puts an effective  !
!  lower bound of -.281 on gh in stably stratified regions                     !
!  length scale near surface is modified with wind mixed length scale          !
!------------------------------------------------------------------------------!

   DO K=2,KBM1
     DO I=1,NNode
       IF(Z(I,K) > -0.5_SP)L(I,K) = MAX(L(I,K),KAPPA*L0(I)) 
       IF(BOYGR(I,K) < 0.0_SP)THEN
         LMAX      = SQRT(ABS(GH_MIN) * Q2(I,K) / (-BOYGR(I,K) + SMALL) )
         L(I,K)    = MIN(L(I,K),LMAX)
         Q2L(I,K) = Q2(I,K)*L(I,K)
       END IF
     END DO
   END DO

!------------------------------------------------------------------------------!
!  solve for gh = (-l^2*N^2/q^2)                                               !
!------------------------------------------------------------------------------!

   GH = (L**2/Q2)*BOYGR

!------------------------------------------------------------------------------!
!  limit maximum of GH in unstable regions to reasonable physical value        !
!  limiting GH also constrains stability functions SH/SM to finite values      !
!------------------------------------------------------------------------------!

   GH = MIN(GH,GH_MAX)

   L(:,1)   = KAPPA*L0(:)         !0.0_SP
   L(:,KB)  = 0.0_SP
   GH(:,1)  = 0.0_SP
   GH(:,KB) = 0.0_SP


!------------------------------------------------------------------------------!
!  calculate eddy kinetic energy dissipation rate: dtef                        !
!------------------------------------------------------------------------------!

   STF = 1.0_SP

   IF(.not.SURFACEWAVE_MIX)THEN
   DO K = 1,KB
     DO I = 1,NNode
       IF(GH(I,K) < 0.0_SP) STF(I,K)=1.0_SP-0.9_SP*(GH(I,K)/GHC)**1.5_SP
       IF(GH(I,K) < GHC) STF(I,K) = 0.1_SP
     END DO
   END DO
   END IF

!   DTEF = Q2*SQRT(Q2)/(B1*Q2L + small)*STF
   DTEF = SQRT(Q2)/(B1*L + small)*STF

   DO K = 2, KBM1
     DO I = 1, NNode
       VHP(I,K) = 1. / (A(I,K)+C(I,K)*(1.-VH(I,K-1))-(2.*DTI*DTEF(I,K)+1.))
       VH(I,K) = A(I,K) * VHP(I,K)
       VHP(I,K) = (-2.*DTI*PROD(I,K)+C(I,K)*VHP(I,K-1)-Q2F(I,K))*VHP(I,K)
     END DO
   END DO

   DO K=1,KBM1
     KI = KB-K
     DO I = 1, NNode
#  if defined (WET_DRY)
     IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
       Q2F(I,KI) = VH(I,KI) * Q2F(I,KI+1) + VHP(I,KI)
#  if defined (WET_DRY)
     END IF
#  endif
     END DO
   END DO

!------------------------------------------------------------------------------!
!      the following section solves the equation                               !
!      dti*(kq*q2l')' - q2l*(dti*dtef+1.) = -q2lb                              !
!------------------------------------------------------------------------------!

   VH(:,1)  = 0.0_SP
   VHP(:,1) = 0.0_SP

   DO  K = 2, KBM1
     DO  I = 1, NNode
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
         DTEF(I,K) = DTEF(I,K) * (1.+E2*((1./ABS(Z(I,K)-Z(I,1))+1./ &
                     ABS(Z(I,K)-Z(I,KB)))*L(I,K)/(D(I)*KAPPA))**2)
         VHP(I,K)  = 1. / (A(I,K)+C(I,K)*(1.-VH(I,K-1))- &
                     (DTI*DTEF(I,K)+1.))
         VH(I,K)   = A(I,K) * VHP(I,K)
         VHP(I,K)  = (DTI*(-PROD(I,K)*L(I,K)*E1)+C(I,K)* &
                     VHP(I,K-1)-Q2LF(I,K)) * VHP(I,K)
       END IF
     END DO
   END DO


   DO K=1,KBM1
     KI = KB - K
     DO I = 1, NNode
#  if defined (WET_DRY)
     IF(ISWET_NODE_CURRENTSTEP(I) == 1)THEN
#  endif
       Q2LF(I,KI) = VH(I,KI) * Q2LF(I,KI+1) + VHP(I,KI)
#  if defined (WET_DRY)
     END IF
#  endif
     END DO
   END DO

!------------------------------------------------------------------------------|
!  place a threshold on the turbulence quantities following vertical diffusion |
!------------------------------------------------------------------------------|

   DO K = 2, KBM1
     DO I = 1, NNode
       IF (Q2F(I,K) <= small .OR. Q2LF(I,K) <= small) THEN
         Q2F(I,K)  = small
         Q2LF(I,K) = small
       END IF
     END DO
   END DO

!===========THRESHOLD ON LENGTH SCALE ALT STYLE (Active in Barotropic Case)====
!   DO K=2,KBM1
!     DO I=1,MElem
!       LMAX      = SQRT(ABS(GH_MIN) * Q2F(I,K) / MAX(0.0_SP,-BOYGR(I,K)+SMALL))
!       L(I,K)    = MIN(L(I,K),LMAX)
!       Q2LF(I,K) = Q2F(I,K)*L(I,K)
!     END DO
!   END DO
!------------------------------------------------------------------------------!

!   L(:,1)   = 0.0_SP
!   L(:,KB)  = 0.0_SP
!   GH(:,1)  = 0.0_SP
!   GH(:,KB) = 0.0_SP

!------------------------------------------------------------------------------!
!  calculate stability functions sh + sm using Galperin 1988 formulation       !
!------------------------------------------------------------------------------!
   COEF4 = 18.0_SP * A1 * A1 + 9.0_SP * A1 * A2
   COEF5 = 9.0_SP * A1 * A2

   DO K = 1,KB
     DO I = 1,NNode
       COEF1 = A2 * (1.0_SP-6.0_SP*A1/B1*STF(I,K))
       COEF2 = 3.0_SP * A2 * B2/STF(I,K) + 18.0_SP * A1 * A2
       COEF3 = A1 * (1.0_SP-3.0_SP*C1-6.0_SP*A1/B1*STF(I,K))

       SH(I,K) = COEF1/(1.0_SP-COEF2*GH(I,K))
       SM(I,K) = (COEF3+SH(I,K)*COEF4*GH(I,K))/(1.-COEF5*GH(I,K))
     END DO
   END DO

!------------------------------------------------------------------------------!
!  calculate turbulent eddy viscosities/diffusivities kq,km,kh                 !
!------------------------------------------------------------------------------!

   KN = L*SQRT(Q2)
   KQ = 0.5_SP*(KN*KAPPA*SM+KQ)
   KM = 0.5_SP*(KN*SM+KM)
!   KH = 0.5_SP*(KN*SH+KH)
   KH = 0.5_SP*(KN*SH*VPRNU+KH)

#  if defined (WET_DRY)
    DO I=1,NNode
     IF(H(I) <= DJUST)THEN
      KQ_TMP = 0.0_SP
      KM_TMP = 0.0_SP
      KH_TMP = 0.0_SP
      DO K=2,KBM1
        KQ_TMP = KQ_TMP + KQ(I,K)
        KM_TMP = KM_TMP + KM(I,K)
        KH_TMP = KH_TMP + KH(I,K)
      END DO
        KQ_TMP = KQ_TMP/(KBM1-1)
        KM_TMP = KM_TMP/(KBM1-1)
        KH_TMP = KH_TMP/(KBM1-1)
      DO K=2,KBM1
        KQ(I,K) = KQ_TMP
        KM(I,K) = KM_TMP
        KH(I,K) = KH_TMP
      END DO
     END IF
    END DO
#  endif

   RETURN
   END SUBROUTINE VDIF_Q
!==============================================================================|
