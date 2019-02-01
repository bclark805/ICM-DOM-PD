!==============================================================================|
!   Calculate the Turbulent Kinetic Energy and Mixing Length Based on          |
!   The Mellor-Yamada Level 2.5 Turbulent Closure Model                        |
!==============================================================================|

   SUBROUTINE ADV_Q_GCN(Q,QF)

!------------------------------------------------------------------------------|
   USE ALL_VARS
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MTElem,KB) :: Q,QF,XFLUX 
   REAL(SP) :: TMPA,TMPB,SITA,DIJ
   REAL(SP) :: UA1,UA2,UA3,UA4,VA1,VA2,VA3,VA4
   REAL(SP) :: QA1,QA2,QA3,QA4,QB1,QB2,QB3,QB4
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UB1,UB2,UB3,UB4,VB1,VB2,VB3,VB4
   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,UIJ,VIJ
   REAL(SP) :: FXX,FYY,XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP
   REAL(SP) :: TMP,TPA,TPB
   REAL(SP) :: COFAQ1,COFAQ2,COFAQ3,COFAQ4,QIJ1,QIJ2
   REAL(SP) :: FACT,FM1,EXFLUX,ISWETTMP
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K

!------------------------------------------------------------------------------|

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Fluxes-----------------------------------------------------------!               
!
   QF    = 0.0_SP
   XFLUX = 0.0_SP

!
!-Loop Over Edges and Accumulate Advection and Diffusion Terms at each Element-!       
!

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
#  if defined (WET_DRY)
     IF(ISWET_CELL_LAST_INT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#  endif
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     DIJ=0.5_SP*(DT(J1)+DT(J2))
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)
#    if defined (SPHERICAL)
     XIJA=DLTXNE(I,1)
     YIJA=DLTYNE(I,1)
     XIJB=DLTXNE(I,2)
     YIJB=DLTYNE(I,2)
#    else
     XIJA=XIJC(I)-XC(IA)
     YIJA=YIJC(I)-YC(IA)
     XIJB=XIJC(I)-XC(IB)
     YIJB=YIJC(I)-YC(IB)
#    endif

     DO K=2,KBM1
       !!LEFT VELOCITY
       UA1=0.5_SP*(U(IA,K)+U(IA,K-1))
       UA2=0.5_SP*(U(K1,K)+U(K1,K-1))
       UA3=0.5_SP*(U(K2,K)+U(K2,K-1))
       UA4=0.5_SP*(U(K3,K)+U(K3,K-1))
       VA1=0.5_SP*(V(IA,K)+V(IA,K-1))
       VA2=0.5_SP*(V(K1,K)+V(K1,K-1))
       VA3=0.5_SP*(V(K2,K)+V(K2,K-1))
       VA4=0.5_SP*(V(K3,K)+V(K3,K-1))
       COFA1=A1U(IA,1)*UA1+A1U(IA,2)*UA2+A1U(IA,3)*UA3+A1U(IA,4)*UA4
       COFA2=A2U(IA,1)*UA1+A2U(IA,2)*UA2+A2U(IA,3)*UA3+A2U(IA,4)*UA4
       COFA5=A1U(IA,1)*VA1+A1U(IA,2)*VA2+A1U(IA,3)*VA3+A1U(IA,4)*VA4
       COFA6=A2U(IA,1)*VA1+A2U(IA,2)*VA2+A2U(IA,3)*VA3+A2U(IA,4)*VA4
       UIJ1=UA1+COFA1*XIJA+COFA2*YIJA
       VIJ1=VA1+COFA5*XIJA+COFA6*YIJA

       !!RIGHT VELOCITY
       UB1=0.5_SP*(U(IB,K)+U(IB,K-1))
       UB2=0.5_SP*(U(K4,K)+U(K4,K-1))
       UB3=0.5_SP*(U(K5,K)+U(K5,K-1))
       UB4=0.5_SP*(U(K6,K)+U(K6,K-1))
       VB1=0.5_SP*(V(IB,K)+V(IB,K-1))
       VB2=0.5_SP*(V(K4,K)+V(K4,K-1))
       VB3=0.5_SP*(V(K5,K)+V(K5,K-1))
       VB4=0.5_SP*(V(K6,K)+V(K6,K-1))
       COFA3=A1U(IB,1)*UB1+A1U(IB,2)*UB2+A1U(IB,3)*UB3+A1U(IB,4)*UB4
       COFA4=A2U(IB,1)*UB1+A2U(IB,2)*UB2+A2U(IB,3)*UB3+A2U(IB,4)*UB4
       COFA7=A1U(IB,1)*VB1+A1U(IB,2)*VB2+A1U(IB,3)*VB3+A1U(IB,4)*VB4
       COFA8=A2U(IB,1)*VB1+A2U(IB,2)*VB2+A2U(IB,3)*VB3+A2U(IB,4)*VB4
       UIJ2=UB1+COFA3*XIJB+COFA4*YIJB
       VIJ2=VB1+COFA7*XIJB+COFA8*YIJB

       !!NORMAL VELOCITY ACROSS THE EDGE
       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))    

       QA1=Q(IA,K)
       QA2=Q(K1,K)
       QA3=Q(K2,K)
       QA4=Q(K3,K)
       COFAQ1=A1U(IA,1)*QA1+A1U(IA,2)*QA2+A1U(IA,3)*QA3+A1U(IA,4)*QA4
       COFAQ2=A2U(IA,1)*QA1+A2U(IA,2)*QA2+A2U(IA,3)*QA3+A2U(IA,4)*QA4
       QIJ1=QA1+COFAQ1*XIJA+COFAQ2*YIJA

       QB1=Q(IB,K)
       QB2=Q(K4,K)
       QB3=Q(K5,K)
       QB4=Q(K6,K)
       COFAQ3=A1U(IB,1)*QB1+A1U(IB,2)*QB2+A1U(IB,3)*QB3+A1U(IB,4)*QB4
       COFAQ4=A2U(IB,1)*QB1+A2U(IB,2)*QB2+A2U(IB,3)*QB3+A2U(IB,4)*QB4
       QIJ2=QB1+COFAQ3*XIJB+COFAQ4*YIJB

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON

       TXXIJ=0.5_SP*(COFAQ1+COFAQ3)*VISCOF
       TYYIJ=0.5_SP*(COFAQ2+COFAQ4)*VISCOF
       FXX= DIJ*TXXIJ*DLTYC(I)
       FYY=-DIJ*TYYIJ*DLTXC(I)

       !!UPWIND THE ADVECTION
       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*QIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*QIJ1)*0.5_SP

       !!AUGMENT FLUX AT THE BOUNDARY 
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)

       !!ACCUMULATE FLUX
       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
     END DO
#  if defined (WET_DRY)
    END IF
#  endif
   END DO

#  if defined (WET_DRY)
   DO I=1,MElem
     ISWETTMP = ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I)
     DO K=1,KBM1
       XFLUX(I,K) = XFLUX(I,K)*ISWETTMP
     END DO
   END DO    	 
#  endif   


   DO I=1,MElem
#  if defined (WET_DRY)
    IF(ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif
     DO K=2,KBM1
       TEMP=W(I,K-1)*Q(I,K-1)-W(I,K+1)*Q(I,K+1)
       XFLUX(I,K)=XFLUX(I,K)+TEMP*ART(I)/(DZ(K-1)+DZ(K))
     END DO
#  if defined (WET_DRY)
    END IF
#  endif
   END DO


   DO  I=1,MElem
#  if !defined (WET_DRY)
    IF(D1(I) > 0.0_SP)THEN
#  else
    IF(ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif
      DO  K=2,KBM1
         QF(I,K)=Q(I,K)*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/D1(I)
     END DO
    END IF 
   END DO

   RETURN
   END SUBROUTINE ADV_Q_GCN
!==============================================================================|




