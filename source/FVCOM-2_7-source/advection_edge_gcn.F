#  if !defined (SEMI_IMPLICIT)
!==============================================================================|
   SUBROUTINE ADVECTION_EDGE_GCN(XFLUX,YFLUX)
!==============================================================================|
!   Calculate the Advection and Diffusion Terms of 3D Velocity Field           |
!   These Terms will be vertically integrated to form the Mean Terms in        |
!   the Gx and Gy Terms of the External Mode Equation                          |
!==============================================================================|

   USE ALL_VARS
   USE BCS
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  if defined (NORTHPOLE)
   USE MOD_NORTHPOLE
#  endif   
#  endif
#  if defined (WET_DRY)
   USE MOD_WD
#  endif

#  if defined (MEAN_FLOW)
   USE MOD_MEANFLOW
   USE MOD_OBCS2
   USE MOD_OBCS3
#  endif

   IMPLICIT NONE
   REAL(SP), INTENT(OUT), DIMENSION(0:MTElem,KB) :: XFLUX,YFLUX
   REAL(SP) :: DIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN_TMP
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2
   REAL(SP) :: ISWETTMP

#  if defined (LIMITED_NO)
   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
#  else
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UALFA,VALFA
   REAL(SP) :: UALFA_TMP,VALFA_TMP
   INTEGER  :: ERROR
   REAL(SP) :: EPS
#  endif
!------------------------------------------------------------------------------|

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Variables--------------------------------------------------------|
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP

!
!--Loop Over Edges and Accumulate Fluxes-For Each Element----------------------|
!
#  if !defined (LIMITED_NO)
   ALLOCATE(UIJ1(NE),VIJ1(NE),UIJ2(NE),VIJ2(NE))
   ALLOCATE(UALFA(0:MTElem),VALFA(0:MTElem))
   ALLOCATE(FXX(NE),FYY(NE))

   DO K=1,KBM1
     UIJ1=0.0_SP;  VIJ1=0.0_SP;  UIJ2=0.0_SP;  VIJ2=0.0_SP
     UALFA=1.0_SP; VALFA=1.0_SP
     FXX=0.0_SP;   FYY=0.0_SP

     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)
#      if defined (WET_DRY)
       IF(ISWET_CELL_LAST_INT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#      endif
       J1=IENODE(I,1)
       J2=IENODE(I,2)
!       DIJ= 0.5_SP*(DT(J1)+DT(J2))

       K1=NBE(IA,1)
       K2=NBE(IA,2)
       K3=NBE(IA,3)
       K4=NBE(IB,1)
       K5=NBE(IB,2)
       K6=NBE(IB,3)
#      if defined (SPHERICAL)
       XIJA=DLTXNE(I,1)
       YIJA=DLTYNE(I,1)
       XIJB=DLTXNE(I,2)
       YIJB=DLTYNE(I,2)
#      else
       XIJA=XIJC(I)-XC(IA)
       YIJA=YIJC(I)-YC(IA)
       XIJB=XIJC(I)-XC(IB)
       YIJB=YIJC(I)-YC(IB)
#      endif

       DIJ= 0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
       !!FORM THE LEFT FLUX
       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)
       UIJ1(I)=COFA1*XIJA+COFA2*YIJA
       VIJ1(I)=COFA5*XIJA+COFA6*YIJA
       UALFA_TMP=ABS(U(IA,K)-U(IB,K))/ABS(UIJ1(I)+EPSILON(EPS))
       VALFA_TMP=ABS(V(IA,K)-V(IB,K))/ABS(VIJ1(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1.0_SP)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1.0_SP)VALFA_TMP = 1.0_SP
       UALFA(IA)=MIN(UALFA(IA),UALFA_TMP)
       VALFA(IA)=MIN(VALFA(IA),VALFA_TMP)
     
       !!FORM THE RIGHT FLUX
       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)
       UIJ2(I)=COFA3*XIJB+COFA4*YIJB
       VIJ2(I)=COFA7*XIJB+COFA8*YIJB
       UALFA_TMP=ABS(U(IA,K)-U(IB,K))/ABS(UIJ2(I)+EPSILON(EPS))
       VALFA_TMP=ABS(V(IA,K)-V(IB,K))/ABS(VIJ2(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1.0_SP)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1.0_SP)VALFA_TMP = 1.0_SP
       UALFA(IB)=MIN(UALFA(IB),UALFA_TMP)
       VALFA(IB)=MIN(VALFA(IB),VALFA_TMP)
     
       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
        
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX(I)=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY(I)=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))
#      if defined (WET_DRY)
       END IF
#      endif
     END DO

     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)
#      if defined (WET_DRY)
       IF(ISWET_CELL_LAST_INT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#      endif
       J1=IENODE(I,1)
       J2=IENODE(I,2)
       DIJ= 0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
       UIJ1(I)=U(IA,K)+UALFA(IA)*UIJ1(I)
       VIJ1(I)=V(IA,K)+VALFA(IA)*VIJ1(I)
       UIJ2(I)=U(IB,K)+UALFA(IB)*UIJ2(I)
       VIJ2(I)=V(IB,K)+VALFA(IB)*VIJ2(I)

#      if defined (LIMITED_1)
       IF(UIJ1(I)> MAX(U(IA,K),U(IB,K)) .OR. UIJ1(I) < MIN(U(IA,K),U(IB,K)) .OR.   &
          UIJ2(I)> MAX(U(IA,K),U(IB,K)) .OR. UIJ2(I) < MIN(U(IA,K),U(IB,K)))THEN
         UIJ1(I)=U(IA,K)
         UIJ2(I)=U(IB,K)
       END IF

       IF(VIJ1(I) > MAX(V(IA,K),V(IB,K)) .OR. VIJ1(I) < MIN(V(IA,K),V(IB,K)) .OR.  &
          VIJ2(I) > MAX(V(IA,K),V(IB,K)) .OR. VIJ2(I) < MIN(V(IA,K),V(IB,K)))THEN
         VIJ1(I)=V(IA,K)
         VIJ2(I)=V(IB,K)
       END IF
#      endif       
       
       !!COMPUTE THE NORMAL VELOCITY ACROSS THE EDGE
       UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
       VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
       UN_TMP=VIJ*DLTXC(I) - UIJ*DLTYC(I)

       !!UPWIND THE ADVECTIVE FLUX
       XADV=DIJ*UN_TMP*((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1(I))*0.5_SP
       YADV=DIJ*UN_TMP*((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1(I))*0.5_SP
        

       !!COMPUTE BOUNDARY FLUX AUGMENTERS   
#      if !defined (MEAN_FLOW)
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)
       !!ACCUMULATE THE FLUX
!       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
!       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
!       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
!       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB

        XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I)))*EPOR(IA)
        YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I)))*EPOR(IA)
        XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I)))*EPOR(IB)
        YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I)))*EPOR(IB)
#       else
        TPA = FLOAT(1-ISBC(I))
        TPB = FLOAT(1-ISBC(I))
        XFLUX(IA,K)=XFLUX(IA,K)+(XADV*TPA+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I))))*IUCP(IA)
        YFLUX(IA,K)=YFLUX(IA,K)+(YADV*TPA+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I))))*IUCP(IA)
        XFLUX(IB,K)=XFLUX(IB,K)-(XADV*TPB+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I))))*IUCP(IB)
        YFLUX(IB,K)=YFLUX(IB,K)-(YADV*TPB+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I))))*IUCP(IB)
#       endif

#       if defined (WET_DRY)
        END IF
#       endif
     END DO
   END DO

   DEALLOCATE(UIJ1,VIJ1,UIJ2,VIJ2)
   DEALLOCATE(UALFA,VALFA)
   DEALLOCATE(FXX,FYY)

#  else

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
#  if defined (WET_DRY)
    IF(ISWET_CELL_LAST_INT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#  endif
     J1=IENODE(I,1)
     J2=IENODE(I,2)
!     DIJ= 0.5_SP*(DT(J1)+DT(J2))


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

     DO K=1,KBM1

       DIJ= 0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))

       !!FORM THE LEFT FLUX
       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)
       UIJ1=U(IA,K)+COFA1*XIJA+COFA2*YIJA
       VIJ1=V(IA,K)+COFA5*XIJA+COFA6*YIJA

       !!FORM THE RIGHT FLUX
       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)
       UIJ2=U(IB,K)+COFA3*XIJB+COFA4*YIJB
       VIJ2=V(IB,K)+COFA7*XIJB+COFA8*YIJB

       !!COMPUTE THE NORMAL VELOCITY ACROSS THE EDGE
       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       UN_TMP=VIJ*DLTXC(I) - UIJ*DLTYC(I)

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
        
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

       !!UPWIND THE ADVECTIVE FLUX
       XADV=DIJ*UN_TMP*((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1)*0.5_SP
       YADV=DIJ*UN_TMP*((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1)*0.5_SP
        

       !!COMPUTE BOUNDARY FLUX AUGMENTERS   
#  if !defined (MEAN_FLOW)
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)
       !!ACCUMULATE THE FLUX
!       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
!       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
!       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
!       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB

        XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
        YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
        XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
        YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
#  else
        TPA = FLOAT(1-ISBC(I))
        TPB = FLOAT(1-ISBC(I))
        XFLUX(IA,K)=XFLUX(IA,K)+(XADV*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I))))*IUCP(IA)
        YFLUX(IA,K)=YFLUX(IA,K)+(YADV*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I))))*IUCP(IA)
        XFLUX(IB,K)=XFLUX(IB,K)-(XADV*TPB+(FXX+3.0_SP*FXX*FLOAT(ISBC(I))))*IUCP(IB)
        YFLUX(IB,K)=YFLUX(IB,K)-(YADV*TPB+(FYY+3.0_SP*FYY*FLOAT(ISBC(I))))*IUCP(IB)
#  endif

     END DO
#  if defined (WET_DRY)
    END IF
#  endif
   END DO
#  endif

#  if defined (SPHERICAL) && (NORTHPOLE)
   CALL ADVECTION_EDGE_XY(XFLUX,YFLUX)
#  endif 

#  if defined (WET_DRY)
   DO I=1,MElem
     ISWETTMP = ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I)
     DO K=1,KBM1
       XFLUX(I,K) = XFLUX(I,K)*ISWETTMP
       YFLUX(I,K) = YFLUX(I,K)*ISWETTMP
     END DO
   END DO
#  endif       	 


!
!--Boundary Conditions on Flux-------------------------------------------------|
!
#  if !defined (MEAN_FLOW)
      DO I=1,MElem
        IF(ISBCE(I) == 2)THEN
           DO K=1,KBM1
             XFLUX(I,K)=0.0_SP
             YFLUX(I,K)=0.0_SP
           END DO
        END IF
      END DO
#  else
   IF (nmfcell_i > 0) THEN
     DO II=1,nmfcell_i
        I1=I_MFCELL_N(II)
        DO K=1,KBM1
           XFLUX(I1,K) = XFLUX(I1,K) + FLUXOBC3D_X_2(II,K)*IUCP(I1)
           YFLUX(I1,K) = YFLUX(I1,K) + FLUXOBC3D_Y_2(II,K)*IUCP(I1)
        END DO
     END DO
   END IF
#  endif

!
!--Adjust Flux for Fresh Water Inflow------------------------------------------|
!

   IF(NUMQBC > 0) THEN
     IF(INFLOW_TYPE == 'node') THEN
       DO II=1,NUMQBC
         J=INODEQ(II)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
!           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ(J,K)
!           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ(J,K)*COS(ANGLEQ(II))
!           XFLUX(I2,K)=XFLUX(I2,K)-TEMP/DZ(J,K)*COS(ANGLEQ(II))
!           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ(J,K)*SIN(ANGLEQ(II))
!           YFLUX(I2,K)=YFLUX(I2,K)-TEMP/DZ(J,K)*SIN(ANGLEQ(II))
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           XFLUX(I2,K)=XFLUX(I2,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
           YFLUX(I2,K)=YFLUX(I2,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO II=1,NUMQBC
         I1=ICELLQ(II)
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
!           TEMP=QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           TEMP=QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ1(I1,K)
!           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ1(I1,K)*COS(ANGLEQ(II))
!           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ1(I1,K)*SIN(ANGLEQ(II))
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     END IF
   END IF

!
!--Adjust Flux for Open Boundary Inflow/Outflow------------------------------------------|
!
#  if defined (MEAN_FLOW)
   IF(nmfcell_i > 0) THEN
     DO II=1,nmfcell_i
       I1=I_MFCELL_N(II)
       DO K=1,KBM1
         VLCTYMF(II)=MFQDIS(II)/MFAREA(II)
!         TEMP=MFQDIS(II)*MFDIST(II,K)*VLCTYMF(II)
         TEMP=MFQDIS(II)*MFDIST(II,K)*MFDIST(II,K)*VLCTYMF(II)/DZ1(I1,K)
!         XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ1(I1,K)*COS(ANGLEMF(II))
!         YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ1(I1,K)*SIN(ANGLEMF(II))
         XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEMF(II))
         YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEMF(II))
       END DO
     END DO
   END IF
#  endif

   RETURN
   END SUBROUTINE ADVECTION_EDGE_GCN
#  endif
!==============================================================================|
