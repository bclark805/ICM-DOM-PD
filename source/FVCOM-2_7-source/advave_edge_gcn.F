!==============================================================================|
!   CALCULATE CONVECTION AND DIFFUSION FLUXES FOR EXTERNAL MODE                !
!==============================================================================|
   SUBROUTINE ADVAVE_EDGE_GCN(XFLUX,YFLUX)
!==============================================================================|

   USE ALL_VARS
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  if defined (NORTHPOLE)
   USE MOD_NORTHPOLE
#  endif   
#  endif
   USE BCS
   USE MOD_OBCS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
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

!add by TW to update velocity block
#  if defined (V_BLOCK)
   USE VELOCITY_BLOCK
#  endif

   IMPLICIT NONE
   INTEGER  :: I,J,K,IA,IB,J1,J2,K1,K2,K3,I1,I2
   REAL(SP) :: DIJ,ELIJ,XIJ,YIJ,UIJ,VIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN_TMP
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP
   REAL(SP) :: XFLUX(0:MTElem),YFLUX(0:MTElem)
   REAL(SP) :: FACT,FM1,ISWETTMP

   REAL(SP) :: TPA,TPB

#  if defined (SPHERICAL)
   REAL(DP) :: XTMP,XTMP1
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

   REAL(SP) :: BTPS
   REAL(SP) :: U_TMP,V_TMP,UAC_TMP,VAC_TMP,WUSURF_TMP,WVSURF_TMP
   REAL(SP) :: WUBOT_TMP,WVBOT_TMP,UAF_TMP,VAF_TMP 
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!-------------------------INITIALIZE FLUXES------------------------------------!
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP
   PSTX  = 0.0_SP
   PSTY  = 0.0_SP

#  if defined (BALANCE_2D)
   ADFXA =0.0_SP
   ADFYA =0.0_SP
#  endif

!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!
#  if !defined (LIMITED_NO)
   ALLOCATE(UIJ1(NE),VIJ1(NE),UIJ2(NE),VIJ2(NE))
   UIJ1=0.0_SP;VIJ1=0.0_SP;UIJ2=0.0_SP;VIJ2=0.0_SP
   
   ALLOCATE(UALFA(0:MTElem),VALFA(0:MTElem))
   UALFA=1.0_SP;VALFA=1.0_SP
   
   ALLOCATE(FXX(NE),FYY(NE))
   FXX=0.0_SP;FYY=0.0_SP

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if !defined (SEMI_IMPLICIT)
     DIJ=0.5_SP*(D(J1)+D(J2))
#    else
     DIJ= 0.5_SP*(DT(J1)+DT(J2))     
#    endif

#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWET_CELL_LAST_EXT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_EXT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#    else
     IF(ISWET_CELL_LAST_INT_STEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB) == 1)THEN
#    endif
#    endif
!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
         
     COFA1=A1U(IA,1)*UA(IA)+A1U(IA,2)*UA(K1)+A1U(IA,3)*UA(K2)+A1U(IA,4)*UA(K3)
     COFA2=A2U(IA,1)*UA(IA)+A2U(IA,2)*UA(K1)+A2U(IA,3)*UA(K2)+A2U(IA,4)*UA(K3)
     COFA5=A1U(IA,1)*VA(IA)+A1U(IA,2)*VA(K1)+A1U(IA,3)*VA(K2)+A1U(IA,4)*VA(K3)
     COFA6=A2U(IA,1)*VA(IA)+A2U(IA,2)*VA(K1)+A2U(IA,3)*VA(K2)+A2U(IA,4)*VA(K3)
     
#    if defined (SPHERICAL)
     UIJ1(I)=COFA1*DLTXNE(I,1)+COFA2*DLTYNE(I,1)
     VIJ1(I)=COFA5*DLTXNE(I,1)+COFA6*DLTYNE(I,1)
#    else
     XIJ=XIJC(I)-XC(IA)
     YIJ=YIJC(I)-YC(IA)
     UIJ1(I)=COFA1*XIJ+COFA2*YIJ
     VIJ1(I)=COFA5*XIJ+COFA6*YIJ
#    endif
     UALFA_TMP=ABS(UA(IA)-UA(IB))/ABS(UIJ1(I)+EPSILON(EPS))
     VALFA_TMP=ABS(VA(IA)-VA(IB))/ABS(VIJ1(I)+EPSILON(EPS))
     IF(UALFA_TMP > 1.0_SP)UALFA_TMP = 1.0_SP
     IF(VALFA_TMP > 1.0_SP)VALFA_TMP = 1.0_SP
     UALFA(IA)=MIN(UALFA(IA),UALFA_TMP)
     VALFA(IA)=MIN(VALFA(IA),VALFA_TMP)

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)
          
     COFA3=A1U(IB,1)*UA(IB)+A1U(IB,2)*UA(K1)+A1U(IB,3)*UA(K2)+A1U(IB,4)*UA(K3)
     COFA4=A2U(IB,1)*UA(IB)+A2U(IB,2)*UA(K1)+A2U(IB,3)*UA(K2)+A2U(IB,4)*UA(K3)
     COFA7=A1U(IB,1)*VA(IB)+A1U(IB,2)*VA(K1)+A1U(IB,3)*VA(K2)+A1U(IB,4)*VA(K3)
     COFA8=A2U(IB,1)*VA(IB)+A2U(IB,2)*VA(K1)+A2U(IB,3)*VA(K2)+A2U(IB,4)*VA(K3)
     
#    if defined (SPHERICAL)
     UIJ2(I)=COFA3*DLTXNE(I,2)+COFA4*DLTYNE(I,2)
     VIJ2(I)=COFA7*DLTXNE(I,2)+COFA8*DLTYNE(I,2)
#    else
     XIJ=XIJC(I)-XC(IB)
     YIJ=YIJC(I)-YC(IB)
     UIJ2(I)=COFA3*XIJ+COFA4*YIJ
     VIJ2(I)=COFA7*XIJ+COFA8*YIJ
#    endif
     UALFA_TMP=ABS(UA(IA)-UA(IB))/ABS(UIJ2(I)+EPSILON(EPS))
     VALFA_TMP=ABS(VA(IA)-VA(IB))/ABS(VIJ2(I)+EPSILON(EPS))
     IF(UALFA_TMP > 1.0_SP)UALFA_TMP = 1.0_SP
     IF(VALFA_TMP > 1.0_SP)VALFA_TMP = 1.0_SP
     UALFA(IB)=MIN(UALFA(IB),UALFA_TMP)
     VALFA(IB)=MIN(VALFA(IB),VALFA_TMP)

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU

!    SHEAR STRESSES
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX(I)=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
     FYY(I)=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))
#    if defined (WET_DRY)
     ENDIF
#    endif
   END DO

   DO I=1, NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if !defined (SEMI_IMPLICIT)
     ELIJ=0.5_SP*(EL(J1)+EL(J2))
     DIJ=0.5_SP*(D(J1)+D(J2))

#    if defined (HEAT_FLUX)
	 if(C_HFX) then
     ELIJ = ELIJ-0.5_SP*(EL_AIR(J1)+EL_AIR(J2))*RAMP
	 endif
#    endif

#    if defined (EQUI_TIDE)
     ELIJ= ELIJ-0.5_SP*(EL_EQI(J1)+EL_EQI(J2))
#    endif
#    if defined (ATMO_TIDE)
     ELIJ= ELIJ-0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))
#    endif       

#    else

     DIJ= 0.5_SP*(DT(J1)+DT(J2))     
     ELIJ=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
#    if defined (HEAT_FLUX)
	 if(C_HFX) then
     ELIJ=ELIJ-((1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) )*RAMP
	 endif
#    endif
#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-((1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-((1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#    endif

#    endif

     UIJ1(I)=UA(IA)+UALFA(IA)*UIJ1(I)
     VIJ1(I)=VA(IA)+VALFA(IA)*VIJ1(I)
     UIJ2(I)=UA(IB)+UALFA(IB)*UIJ2(I)
     VIJ2(I)=VA(IB)+VALFA(IB)*VIJ2(I)

#    if defined (LIMITED_1)
     IF(UIJ1(I) > MAX(UA(IA),UA(IB)) .OR. UIJ1(I) < MIN(UA(IA),UA(IB))  .OR.  &
        UIJ2(I) > MAX(UA(IA),UA(IB)) .OR. UIJ2(I) < MIN(UA(IA),UA(IB)))THEN
       UIJ1(I)=UA(IA)
       UIJ2(I)=UA(IB)
     END IF
      
     IF(VIJ1(I) > MAX(VA(IA),VA(IB)) .OR. VIJ1(I) < MIN(VA(IA),VA(IB))  .OR.  &
        VIJ2(I) > MAX(VA(IA),VA(IB)) .OR. VIJ2(I) < MIN(VA(IA),VA(IB)))THEN
       VIJ1(I)=VA(IA)
       VIJ2(I)=VA(IB)
     END IF
#    endif

!    NORMAL VELOCITY
     UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
     VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
     UN_TMP=-UIJ*DLTYC(I) + VIJ*DLTXC(I)
          
#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWET_CELL_LAST_EXT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_EXT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#    else
     IF(ISWET_CELL_LAST_INT_STEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB) == 1)THEN
#    endif
#    endif
!    ADD CONVECTIVE AND VISCOUS FLUXES
     XADV=DIJ*UN_TMP*&
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1(I))*0.5_SP
     YADV=DIJ*UN_TMP* &
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1(I))*0.5_SP

!    ACCUMULATE FLUX
#  if !defined (MEAN_FLOW)
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX(I)*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY(I)*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX(I)*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY(I)*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
#  else
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX(I))*(1.0_SP-ISBC(I))*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY(I))*(1.0_SP-ISBC(I))*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX(I))*(1.0_SP-ISBC(I))*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY(I))*(1.0_SP-ISBC(I))*IUCP(IB)
#    endif

#    if defined (WET_DRY)
     END IF
#    endif

#  if defined (BALANCE_2D)
     ADFXA(IA)=ADFXA(IA)+FXX(I)
     ADFYA(IA)=ADFYA(IA)+FYY(I)
     ADFXA(IB)=ADFXA(IB)-FXX(I)
     ADFYA(IB)=ADFYA(IB)-FYY(I)
#  endif


!    ACCUMULATE BAROTROPIC FLUX
!for spherical coordinator and domain across 360^o longitude         
#    if defined (SPHERICAL)
     XTMP  = VX(J2)*TPI-VX(J1)*TPI
     XTMP1 = VX(J2)-VX(J1)
     IF(XTMP1 >  180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP =  360.0_SP*TPI+XTMP
     END IF  

#    if !defined (SEMI_IMPLICIT)
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))

     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#    else
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))

     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#    endif
#    else
#    if !defined (SEMI_IMPLICIT)
#     if defined (V_BLOC)	 
		PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I) * (1.0 - bwt(ia))    !depth adjustment added by TW for block
		PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC(I) * (1.0 - bwt(ia))
		PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I) * (1.0 - bwt(ib))
		PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC(I) * (1.0 - bwt(ib))
#     else
		PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
		PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC(I)
		PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
		PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC(I)
#     endif
#    else
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*DLTXC(I)
     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*DLTXC(I)
#    endif
#    endif     

   END DO

#  else
   
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if !defined (SEMI_IMPLICIT)
     ELIJ = 0.5_SP*(EL(J1)+EL(J2))
     DIJ  = 0.5_SP*(D(J1)+D(J2))

#    if defined (HEAT_FLUX)
	 if(C_HFX) then
     ELIJ = ELIJ-0.5_SP*(EL_AIR(J1)+EL_AIR(J2))*RAMP
	 endif
#    endif

#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-0.5_SP*(EL_EQI(J1)+EL_EQI(J2))
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))
#    endif       

#    else

     DIJ= 0.5_SP*(DT(J1)+DT(J2))     
     ELIJ=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
#    if defined (HEAT_FLUX)
	 if(C_HFX) then
     ELIJ=ELIJ-((1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) )*RAMP
	 endif
#    endif
#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-((1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-((1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#    endif

#    endif

#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWET_CELL_LAST_EXT_STEP(IA)*ISWET_CELL_CURRENTSTEP(IA) == 1 .OR. ISWET_CELL_LAST_EXT_STEP(IB)*ISWET_CELL_CURRENTSTEP(IB) == 1)THEN
#    else
     IF(ISWET_CELL_LAST_INT_STEP(IA) == 1 .OR. ISWET_CELL_LAST_INT_STEP(IB) == 1)THEN
#    endif
#    endif
!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
         
     COFA1=A1U(IA,1)*UA(IA)+A1U(IA,2)*UA(K1)+A1U(IA,3)*UA(K2)+A1U(IA,4)*UA(K3)
     COFA2=A2U(IA,1)*UA(IA)+A2U(IA,2)*UA(K1)+A2U(IA,3)*UA(K2)+A2U(IA,4)*UA(K3)
     COFA5=A1U(IA,1)*VA(IA)+A1U(IA,2)*VA(K1)+A1U(IA,3)*VA(K2)+A1U(IA,4)*VA(K3)
     COFA6=A2U(IA,1)*VA(IA)+A2U(IA,2)*VA(K1)+A2U(IA,3)*VA(K2)+A2U(IA,4)*VA(K3)
     
#    if defined (SPHERICAL)
     UIJ1=UA(IA)+COFA1*DLTXNE(I,1)+COFA2*DLTYNE(I,1)
     VIJ1=VA(IA)+COFA5*DLTXNE(I,1)+COFA6*DLTYNE(I,1)
#    else
     XIJ=XIJC(I)-XC(IA)
     YIJ=YIJC(I)-YC(IA)
     UIJ1=UA(IA)+COFA1*XIJ+COFA2*YIJ
     VIJ1=VA(IA)+COFA5*XIJ+COFA6*YIJ
#    endif

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)
          
     COFA3=A1U(IB,1)*UA(IB)+A1U(IB,2)*UA(K1)+A1U(IB,3)*UA(K2)+A1U(IB,4)*UA(K3)
     COFA4=A2U(IB,1)*UA(IB)+A2U(IB,2)*UA(K1)+A2U(IB,3)*UA(K2)+A2U(IB,4)*UA(K3)
     COFA7=A1U(IB,1)*VA(IB)+A1U(IB,2)*VA(K1)+A1U(IB,3)*VA(K2)+A1U(IB,4)*VA(K3)
     COFA8=A2U(IB,1)*VA(IB)+A2U(IB,2)*VA(K1)+A2U(IB,3)*VA(K2)+A2U(IB,4)*VA(K3)
     
#    if defined (SPHERICAL)
     UIJ2=UA(IB)+COFA3*DLTXNE(I,2)+COFA4*DLTYNE(I,2)
     VIJ2=VA(IB)+COFA7*DLTXNE(I,2)+COFA8*DLTYNE(I,2)
#    else
     XIJ=XIJC(I)-XC(IB)
     YIJ=YIJC(I)-YC(IB)
     UIJ2=UA(IB)+COFA3*XIJ+COFA4*YIJ
     VIJ2=VA(IB)+COFA7*XIJ+COFA8*YIJ
#    endif

!    NORMAL VELOCITY
     UIJ=0.5_SP*(UIJ1+UIJ2)
     VIJ=0.5_SP*(VIJ1+VIJ2)
     UN_TMP=-UIJ*DLTYC(I) + VIJ*DLTXC(I)

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU

!    SHEAR STRESSES
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
     FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

!    ADD CONVECTIVE AND VISCOUS FLUXES
     XADV=DIJ*UN_TMP*&
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1)*0.5_SP
     YADV=DIJ*UN_TMP* &
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1)*0.5_SP

!    ACCUMULATE FLUX
#  if !defined (MEAN_FLOW)
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
#  else
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX)*(1.0_SP-ISBC(I))*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY)*(1.0_SP-ISBC(I))*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX)*(1.0_SP-ISBC(I))*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY)*(1.0_SP-ISBC(I))*IUCP(IB)
#    endif

#    if defined (WET_DRY)
     END IF
#    endif

#  if defined (BALANCE_2D)
     ADFXA(IA)=ADFXA(IA)+FXX
     ADFYA(IA)=ADFYA(IA)+FYY
     ADFXA(IB)=ADFXA(IB)-FXX
     ADFYA(IB)=ADFYA(IB)-FYY
#  endif


!    ACCUMULATE BAROTROPIC FLUX
!for spherical coordinator and domain across 360^o longitude         
#    if defined (SPHERICAL)
     XTMP  = VX(J2)*TPI-VX(J1)*TPI
     XTMP1 = VX(J2)-VX(J1)
     IF(XTMP1 >  180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP =  360.0_SP*TPI+XTMP
     END IF  

#    if !defined (SEMI_IMPLICIT)
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))

     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#    else
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))

     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#    endif
#    else
#    if !defined (SEMI_IMPLICIT)
#	 	if defined (V_BLOCK)
			PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I) * (1.0 - bwt(ia))    !depth adjustment added by TW for block
			PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC(I) * (1.0 - bwt(ia))
			PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I) * (1.0 - bwt(ib))
			PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC(I) * (1.0 - bwt(ib))
#		else
			PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
			PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC(I)
			PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
			PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC(I)
#		endif
#    else
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*DLTXC(I)
     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*DLTXC(I)
#    endif
#    endif     

   END DO
#  endif

#  if !defined (SEMI_IMPLICIT)
#  if defined (SPHERICAL)  && (NORTHPOLE)
   CALL ADVAVE_EDGE_XY(XFLUX,YFLUX)
#  endif  
#  endif 

#  if defined (WET_DRY)
   DO I = 1,MElem
#    if !defined (SEMI_IMPLICIT)
     ISWETTMP = ISWET_CELL_LAST_EXT_STEP(I)*ISWET_CELL_CURRENTSTEP(I)
#    else
     ISWETTMP = ISWET_CELL_LAST_INT_STEP(I)
#    endif
     XFLUX(I) = XFLUX(I)*ISWETTMP
     YFLUX(I) = YFLUX(I)*ISWETTMP
   END DO
#  endif   


!
!-------------------------SET BOUNDARY VALUES----------------------------------!
!

!  MODIFY BOUNDARY FLUX
#  if !defined (MEAN_FLOW)
      DO I=1,MElem
        IF(ISBCE(I) == 2) THEN
#         if !defined (SEMI_IMPLICIT)
          XFLUX(I)=(XFLUX(I)+Fluxobn(I)*UA(I))*IUCP(I)
          YFLUX(I)=(YFLUX(I)+Fluxobn(I)*VA(I))*IUCP(I)
#         else
          XFLUX(I)=0.0_SP
          YFLUX(I)=0.0_SP
#         endif
        ENDIF
      END DO
#  else
   IF (nmfcell_i > 0) THEN
     DO K=1,nmfcell_i
       I1=I_MFCELL_N(K)
       XFLUX(I1) = XFLUX(I1) + FLUXOBC2D_X(K)*IUCP(I1)
       YFLUX(I1) = YFLUX(I1) + FLUXOBC2D_Y(K)*IUCP(I1)
     END DO
   END IF
#  endif 


!  ADJUST FLUX FOR RIVER INFLOW
   IF(NUMQBC > 0) THEN
     IF(INFLOW_TYPE == 'node')THEN
       DO K=1,NUMQBC
         J=INODEQ(K)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         VLCTYQ(K)=QDIS(K)/QAREA(K)
         XFLUX(I1)=XFLUX(I1)-0.5_SP*QDIS(K)*VLCTYQ(K)*COS(ANGLEQ(K))
         YFLUX(I1)=YFLUX(I1)-0.5_SP*QDIS(K)*VLCTYQ(K)*SIN(ANGLEQ(K))
         XFLUX(I2)=XFLUX(I2)-0.5_SP*QDIS(K)*VLCTYQ(K)*COS(ANGLEQ(K))
         YFLUX(I2)=YFLUX(I2)-0.5_SP*QDIS(K)*VLCTYQ(K)*SIN(ANGLEQ(K))
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO K=1,NUMQBC
         I1=ICELLQ(K)
         VLCTYQ(K)=QDIS(K)/QAREA(K)
         TEMP=QDIS(K)*VLCTYQ(K)
         XFLUX(I1)=XFLUX(I1)-TEMP*COS(ANGLEQ(K))
         YFLUX(I1)=YFLUX(I1)-TEMP*SIN(ANGLEQ(K))
       END DO
     END IF
   END IF

!  ADJUST FLUX FOR OPEN BOUNDARY MEAN FLOW
#  if defined (MEAN_FLOW)
   IF(nmfcell_i > 0) THEN
     DO K=1,nmfcell_i
       I1=I_MFCELL_N(K)
       VLCTYMF(K)=MFQDIS(K)/MFAREA(K)
       TEMP=MFQDIS(K)*VLCTYMF(K)
       XFLUX(I1)=XFLUX(I1)-TEMP*COS(ANGLEMF(K))
       YFLUX(I1)=YFLUX(I1)-TEMP*SIN(ANGLEMF(K))
     END DO
   END IF
#  endif

#  if defined (SEMI_IMPLICIT) && (TWO_D_MODEL) 
   DO I=1, MTElem
     XFLUX(I) = XFLUX(I) + PSTX(I) - COR(I)*VA(I)*DT1(I)*ART(I)*EPOR(I)
     YFLUX(I) = YFLUX(I) + PSTY(I) + COR(I)*UA(I)*DT1(I)*ART(I)*EPOR(I)
!     XFLUX(I) = PSTX(I)
!     YFLUX(I) = PSTY(I)     
#    if defined (SPHERICAL)   
     XFLUX(I) = XFLUX(I) - UA(I)*VA(I)/REARTH*TAN(DEG2RAD*YC(I))*DT1(I)*ART(I)*EPOR(I)
     YFLUX(I) = YFLUX(I) + UA(I)*UA(I)/REARTH*TAN(DEG2RAD*YC(I))*DT1(I)*ART(I)*EPOR(I)
#    endif
     IF(ADCOR_ON) THEN
       UBETA2D(I) = XFLUX(I) + COR(I)*VA(I)*DT1(I)*ART(I)*EPOR(I)
       VBETA2D(I) = YFLUX(I) - COR(I)*UA(I)*DT1(I)*ART(I)*EPOR(I)
     ENDIF
   ENDDO
#  endif

#  if !defined (LIMITED_NO) 
   DEALLOCATE(UIJ1,VIJ1,UIJ2,VIJ2)
   DEALLOCATE(UALFA,VALFA)
   DEALLOCATE(FXX,FYY)
#  endif

   RETURN
   END SUBROUTINE ADVAVE_EDGE_GCN
!==============================================================================|
