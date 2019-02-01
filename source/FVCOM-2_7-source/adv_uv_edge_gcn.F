!==============================================================================!

   SUBROUTINE ADV_UV_EDGE_GCN

!==============================================================================!
! this subroutine calculate advective, coriolis, pressure gradient, etc in     !
! x and y momentum equations except vertical diffusion terms for internal mode ! 
!==============================================================================!

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
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif

#  if defined (MEAN_FLOW)
   USE MOD_MEANFLOW
   USE MOD_OBCS2
   USE MOD_OBCS3
#  endif

#  if defined (HEAT_FLUX)
   USE MOD_HEATFLUX
#  endif   

#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif

!TW, added to update velocity block
#  if defined (V_BLOCK)
   USE VELOCITY_BLOCK
#  endif

#  if defined (KELP)
   USE MOD_KELP
#  endif
!! finish addition, T.W., April 2013


   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:MTElem,KB),YFLUX(0:MTElem,KB)
   REAL(SP) :: PSTX_TM(0:MTElem,KB),PSTY_TM(0:MTElem,KB)
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: DIJ,ELIJ,TMPA,TMPB,TMP,XFLUXV,YFLUXV
   REAL(SP) :: FACT,FM1,EXFLUX,ISWETTMP
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2
#  if defined (SPHERICAL)
   REAL(DP) :: XTMP,XTMP1
#  if defined (NORTHPOLE)
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP
#  endif   
#  endif

#  if defined (LIMITED_NO)
   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
#  else
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UALFA,VALFA
   REAL(SP) :: UALFA_TMP,VALFA_TMP
   INTEGER :: ERROR
   REAL(SP) :: EPS
#  endif
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!-----Initialize Flux Variables------------------------------------------------!
!
   XFLUX  = 0.0_SP
   YFLUX  = 0.0_SP
   PSTX_TM = 0.0_SP
   PSTY_TM = 0.0_SP
# if defined (SEMI_IMPLICIT)
   XFLUX3= 0.0_SP
   YFLUX3= 0.0_SP
# endif

!
!-----Loop Over Edges and Accumulate Flux--------------------------------------!
!
#  if !defined (LIMITED_NO)
   ALLOCATE(UIJ1(NE),VIJ1(NE),UIJ2(NE),VIJ2(NE))
   ALLOCATE(UALFA(0:MTElem),VALFA(0:MTElem))
   ALLOCATE(FXX(NE),FYY(NE))
   
   DO K=1,KBM1
     UIJ1=0.0_SP;VIJ1=0.0_SP;UIJ2=0.0_SP;VIJ2=0.0_SP
     UALFA=1.0_SP;VALFA=1.0_SP
     FXX=0.0_SP;FYY=0.0_SP
     
     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)
       J1=IENODE(I,1)
       J2=IENODE(I,2)
!       DIJ=0.5_SP*(DT(J1)+DT(J2))

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

       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
#      if !defined (SEMI_IMPLICIT)
       IF(ISWET_CELL_LAST_INT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#      else
       IF(ISWET_CELL_LAST_INT_STEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB) == 1)THEN
#      endif
#      endif
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
       
!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON/HPRNU

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
       J1=IENODE(I,1)
       J2=IENODE(I,2)
!       DIJ=0.5_SP*(DT(J1)+DT(J2))

# if !defined (SEMI_IMPLICIT)
       ELIJ=0.5_SP*(EGF(J1)+EGF(J2))

#      if defined (HEAT_FLUX)
	   if(C_HFX) then
       ELIJ=ELIJ-0.5_SP*(EGF_AIR(J1)+EGF_AIR(J2))
	   endif
#      endif

#      if defined (EQUI_TIDE)
       ELIJ=ELIJ-0.5_SP*(EGF_EQI(J1)+EGF_EQI(J2))
#      endif
#      if defined (ATMO_TIDE)
       ELIJ=ELIJ-0.5_SP*(EGF_ATMO(J1)+EGF_ATMO(J2))
#      endif       

# else

       ELIJ=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
#      if defined (HEAT_FLUX)
	   if(C_HFX) then
       ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) ) 
	   endif
#      endif
#      if defined (EQUI_TIDE)
       ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#      endif
#      if defined (ATMO_TIDE)
       ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#      endif       

# endif

       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
#      if !defined (SEMI_IMPLICIT)
       IF(ISWET_CELL_LAST_INT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#      else
       IF(ISWET_CELL_LAST_INT_STEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB) == 1)THEN
#      endif
#      endif       
       UIJ1(I)=U(IA,K)+UALFA(IA)*UIJ1(I)
       VIJ1(I)=V(IA,K)+VALFA(IA)*VIJ1(I)
       UIJ2(I)=U(IB,K)+UALFA(IB)*UIJ2(I)
       VIJ2(I)=V(IB,K)+VALFA(IB)*VIJ2(I)

#      if defined (LIMITED_1)
       IF(UIJ1(I) > MAX(U(IA,K),U(IB,K)) .OR. UIJ1(I) < MIN(U(IA,K),U(IB,K)) .OR.   &
          UIJ2(I) > MAX(U(IA,K),U(IB,K)) .OR. UIJ2(I) < MIN(U(IA,K),U(IB,K)))THEN
         UIJ1(I)=U(IA,K)
         UIJ2(I)=U(IB,K)
       END IF

       IF(VIJ1(I) > MAX(V(IA,K),V(IB,K)) .OR. VIJ1(I) < MIN(V(IA,K),V(IB,K)) .OR.   &
          VIJ2(I) > MAX(V(IA,K),V(IB,K)) .OR. VIJ2(I) < MIN(V(IA,K),V(IB,K)))THEN
         VIJ1(I)=V(IA,K)
         VIJ2(I)=V(IB,K)
       END IF
#      endif

!      NORMAL VELOCITY              
       UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
       VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))

       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,EXFLUX))*UIJ1(I))*0.5_SP
       YADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,EXFLUX))*VIJ1(I))*0.5_SP

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
#  if !defined (MEAN_FLOW)
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)

       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
!       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
!       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
!       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
!       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB
       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I)))*EPOR(IA)
       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I)))*EPOR(IA)
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I)))*EPOR(IB)
       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I)))*EPOR(IB)
#  else
       TPA = FLOAT(1-ISBC(I))
       TPB = FLOAT(1-ISBC(I))
       XFLUX(IA,K)=XFLUX(IA,K)+(XADV*TPA+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I))))*IUCP(IA)
       YFLUX(IA,K)=YFLUX(IA,K)+(YADV*TPA+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I))))*IUCP(IA)
       XFLUX(IB,K)=XFLUX(IB,K)-(XADV*TPB+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I))))*IUCP(IB)
       YFLUX(IB,K)=YFLUX(IB,K)-(YADV*TPB+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I))))*IUCP(IB)
#  endif

#  if defined (WET_DRY)
    END IF
#  endif

! for spherical coordinator and domain across 360^o longitude
#    if defined (SPHERICAL)
        XTMP  = VX(J2)*TPI-VX(J1)*TPI
        XTMP1 = VX(J2)-VX(J1)
        IF(XTMP1 >  180.0_SP)THEN
          XTMP = -360.0_SP*TPI+XTMP
        ELSE IF(XTMP1 < -180.0_SP)THEN
          XTMP =  360.0_SP*TPI+XTMP
        END IF
  
!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
#    else
!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*DLTXC(I)
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*DLTXC(I)
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)
#    endif

     END DO
   END DO

   DEALLOCATE(UIJ1,VIJ1,UIJ2,VIJ2)
   DEALLOCATE(UALFA,VALFA)
   DEALLOCATE(FXX,FYY)

#  else

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
!     DIJ=0.5_SP*(DT(J1)+DT(J2))

#    if !defined (SEMI_IMPLICIT)

     ELIJ=0.5_SP*(EGF(J1)+EGF(J2))

#    if defined (HEAT_FLUX)
	 if(C_HFX) then
     ELIJ=ELIJ-0.5_SP*(EGF_AIR(J1)+EGF_AIR(J2))
	endif
#    endif

#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-0.5_SP*(EGF_EQI(J1)+EGF_EQI(J2))
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-0.5_SP*(EGF_ATMO(J1)+EGF_ATMO(J2))
#    endif       

#    else

     ELIJ=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
#    if defined (HEAT_FLUX)
	 if(C_HFX) then
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) ) 
	 endif
#    endif
#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#    endif       

#    endif

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
       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#     if defined (WET_DRY)
#     if !defined (SEMI_IMPLICIT)
      IF(ISWET_CELL_LAST_INT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#     else
      IF(ISWET_CELL_LAST_INT_STEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB) == 1)THEN
#     endif
#     endif
       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)

       UIJ1=U(IA,K)+COFA1*XIJA+COFA2*YIJA
       VIJ1=V(IA,K)+COFA5*XIJA+COFA6*YIJA

       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)

       UIJ2=U(IB,K)+COFA3*XIJB+COFA4*YIJB
       VIJ2=V(IB,K)+COFA7*XIJB+COFA8*YIJB

       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))

!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON/HPRNU

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*UIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*UIJ1)*0.5_SP
       YADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*VIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*VIJ1)*0.5_SP

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
#  if !defined (MEAN_FLOW)
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)

       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
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

#  if defined (WET_DRY)
    END IF
#  endif

! for spherical coordinator and domain across 360^o longitude
#    if defined (SPHERICAL)
        XTMP  = VX(J2)*TPI-VX(J1)*TPI
        XTMP1 = VX(J2)-VX(J1)
        IF(XTMP1 >  180.0_SP)THEN
          XTMP = -360.0_SP*TPI+XTMP
        ELSE IF(XTMP1 < -180.0_SP)THEN
          XTMP =  360.0_SP*TPI+XTMP
        END IF
  
!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
#    else
!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*DLTXC(I)
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*DLTXC(I)
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)
#    endif     

     END DO
   END DO
#  endif

      DO I=1,MElem
#     if defined (WET_DRY)
#       if !defined (SEMI_IMPLICIT)
        ISWETTMP = ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I)
#       else
        ISWETTMP = ISWET_CELL_LAST_INT_STEP(I)
#       endif
        DO K=1,KBM1
	 XFLUX(I,K)  = XFLUX(I,K)*ISWETTMP
	 YFLUX(I,K)  = YFLUX(I,K)*ISWETTMP
         PSTX_TM(I,K)= PSTX_TM(I,K)*ISWETTMP
         PSTY_TM(I,K)= PSTY_TM(I,K)*ISWETTMP
        END DO
#     endif
       DO K=1,KBM1
        XFLUX(I,K)=XFLUX(I,K)+PSTX_TM(I,K)
        YFLUX(I,K)=YFLUX(I,K)+PSTY_TM(I,K)
       END DO
      END DO

!
!-------ADD VERTICAL CONVECTIVE FLUX, CORIOLIS TERM AND BAROCLINIC PG TERM----!
!
   DO I=1,MElem
#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#    else
     IF(ISWET_CELL_LAST_INT_STEP(I) == 1)THEN
#    endif
#    endif
     DO K=1,KBM1
       IF(K == 1) THEN
         XFLUXV=-W(I,K+1)*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV=-W(I,K+1)*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       ELSE IF(K == KBM1) THEN
         XFLUXV= W(I,K)*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
         YFLUXV= W(I,K)*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
       ELSE
         XFLUXV= W(I,K)*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 W(I,K+1)*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV= W(I,K)*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 W(I,K+1)*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       END IF
#      if defined (SPHERICAL)
!       XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
!                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*ART(I)&
!                 -U(I,K)*V(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*ART(I)&
!                 +0.5_SP*U(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*ART(I)
!       YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
!                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*ART(I)&
!                 +U(I,K)*U(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*ART(I)&
!                 +0.5_SP*V(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*ART(I)
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)&
                 -U(I,K)*V(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*DZ1(I,K)*ART(I)&
                 +0.5_SP*U(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*DZ1(I,K)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)&
                 +U(I,K)*U(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*DZ1(I,K)*ART(I)&
                 +0.5_SP*V(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*DZ1(I,K)*ART(I)
#      else
!       XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
!                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*ART(I)
!       YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
!                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*ART(I)
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)

!update momentum sink by kelp bed, T.W., April 2013
#  if defined (KELP)
   IF(C_KELP) THEN
     XFLUX(I,K) = XFLUX(I,K) - EMS_X(I,K)
     YFLUX(I,K) = YFLUX(I,K) - EMS_Y(I,K)
   END IF
#  endif


#      endif

     END DO
#  if defined (WET_DRY)
    END IF
#  endif
   END DO

#  if defined (SPHERICAL)  && (NORTHPOLE)
#  if !defined (SEMI_IMPLICIT)
     CALL ADV_UV_EDGE_XY(XFLUX,YFLUX,0.0_SP)
#  else
     CALL ADV_UV_EDGE_XY(XFLUX,YFLUX,IFCETA)
#  endif
#  endif

#  if !defined (MEAN_FLOW)
      DO I=1,MElem
         IF(ISBCE(I) == 2) THEN
            DO K=1,KBM1
#              if !defined (SEMI_IMPLICIT)
               XFLUX(I,K)=0.0_SP
               YFLUX(I,K)=0.0_SP
#              else
               XFLUX(I,K)=PSTX_TM(I,K)
               YFLUX(I,K)=PSTY_TM(I,K)
#              endif
            END DO
         END IF
      END DO
#  else
   IF(nmfcell_i > 0) THEN
     DO II=1,nmfcell_i
       I1=I_MFCELL_N(II)
       DO K=1,KBM1
         XFLUX(I1,K) = XFLUX(I1,K) + FLUXOBC3D_X(II,K)*IUCP(I1)
         YFLUX(I1,K) = YFLUX(I1,K) + FLUXOBC3D_Y(II,K)*IUCP(I1)
       END DO
     END DO
   END IF
#  endif

   !ADJUST FLUX AT RIVER INFLOWS
   IF(NUMQBC >= 1) THEN
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
     ELSE
       PRINT*,'INFLOW_TYPE NOT CORRECT'
       CALL PSTOP
     END IF
   END IF

   !ADJUST FLUX AT OPEN BOUNDARY MEAN FLOW
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

   DO I=1,MElem
#  if defined (WET_DRY)
#   if !defined (SEMI_IMPLICIT)
    IF(ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#   else
    IF(ISWET_CELL_LAST_INT_STEP(I) == 1)THEN
#   endif
#  endif

#  if defined (SPHERICAL) && (NORTHPOLE)
   IF(CELL_NORTHAREA(I) == 1)THEN
     DO K=1,KBM1
       U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
       V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)

#      if !defined (SEMI_IMPLICIT)
!       UF_TMP=U_TMP*DT1(I,K)/D1(I,K)-DTI*XFLUX(I,K)/ART(I)/D1(I,K)
!       VF_TMP=V_TMP*DT1(I,K)/D1(I,K)-DTI*YFLUX(I,K)/ART(I)/D1(I,K)
       UF_TMP=U_TMP*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
       VF_TMP=V_TMP*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))

       IF (ADCOR_ON) THEN
         UF(I,K)  = UF_TMP
         VF(I,K)  = VF_TMP
         UBETA(I,K)=XFLUX(I,K) +COR(I)*V_TMP*DT1(I)*DZ1(I,K)*ART(I)
         VBETA(I,K)=YFLUX(I,K) -COR(I)*U_TMP*DT1(I)*DZ1(I,K)*ART(I)
       ELSE
         UF(I,K)  = VF_TMP*COS(XC(I)*DEG2RAD)-UF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = UF_TMP*COS(XC(I)*DEG2RAD)+VF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = -VF(I,K)
       ENDIF
#      else
       XFLUX3(I,K)=XFLUX(I,K)
       YFLUX3(I,K)=YFLUX(I,K)
       IF (ADCOR_ON) THEN
         UBETA(I,K)=XFLUX(I,K) +COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)
         VBETA(I,K)=YFLUX(I,K) -COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)
         UBETA_NP(I,K)=XFLUX3_NP(I,K) +COR(I)*V_TMP*DT1(I)*DZ1(I,K)*ART(I)
         VBETA_NP(I,K)=YFLUX3_NP(I,K) -COR(I)*U_TMP*DT1(I)*DZ1(I,K)*ART(I)
       ENDIF
#      endif

     END DO
   ELSE
#  endif   

     DO K=1,KBM1
#      if !defined (SEMI_IMPLICIT)
!       UF(I,K)=U(I,K)*DT1(I,K)/D1(I,K)-DTI*XFLUX(I,K)/ART(I)/D1(I,K)
!       VF(I,K)=V(I,K)*DT1(I,K)/D1(I,K)-DTI*YFLUX(I,K)/ART(I)/D1(I,K)
       UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
       VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))

       IF(ADCOR_ON) THEN
         UBETA(I,K)=XFLUX(I,K) +COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)
         VBETA(I,K)=YFLUX(I,K) -COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)
       ENDIF
#      else
       XFLUX3(I,K)=XFLUX(I,K)
       YFLUX3(I,K)=YFLUX(I,K)
       IF(ADCOR_ON) THEN
         UBETA(I,K)=XFLUX(I,K) +COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
         VBETA(I,K)=YFLUX(I,K) -COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
       ENDIF
#      endif
     END DO

#  if defined (SPHERICAL) && (NORTHPOLE)
   END IF
#  endif        

#  if defined (WET_DRY)
    ELSE
     DO K=1,KBM1
       UF(I,K)=0.0_SP
       VF(I,K)=0.0_SP
     END DO
    END IF
#  endif
   END DO

   RETURN
   END SUBROUTINE ADV_UV_EDGE_GCN
!==============================================================================!
