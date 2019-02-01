!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|

   SUBROUTINE ADV_T               

!------------------------------------------------------------------------------|

   USE ALL_VARS
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   USE BCS
   USE MOD_OBCS
#  if defined (WET_DRY)
   USE MOD_WD 
#  endif
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  if defined (NORTHPOLE)
   USE MOD_NORTHPOLE
#  endif   
#  endif
#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:NTNode,KB)      :: XFLUX,XFLUX_ADV,RF
   REAL(SP), DIMENSION(0:NTNode)         :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(0:NTNode)         :: PTPX,PTPY,PTPXD,PTPYD,VISCOFF
   REAL(SP), DIMENSION(3*(MTElem),KBM1)  :: DTIJ 
   REAL(SP), DIMENSION(3*(MTElem),KBM1)  :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN,TTIME,ZDEP
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT,STPOINT1,STPOINT2
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: T1MIN, T1MAX, T2MIN, T2MAX
#  if defined (SPHERICAL)
   REAL(DP) TY,TXPI,TYPI
   REAL(DP) :: XTMP1,XTMP
   REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
   REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
#  if defined (NORTHPOLE)
   REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
   REAL(DP) :: TXPI_TMP,TYPI_TMP
#  endif
#  endif
#  if defined (SEMI_IMPLICIT)
   REAL(SP) :: UN1
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: UVN1
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: DTIJ1 
#  endif
#  if defined (MPDATA)
   REAL(SP) :: TMIN,TMAX,XXXX
   REAL(SP), DIMENSION(0:NTNode,KB)     :: T1_S    !! temporary salinity in modified upwind
   REAL(SP), DIMENSION(0:NTNode,KB)     :: T1_SF   !! temporary salinity in modified upwind
   REAL(SP), DIMENSION(0:NTNode,KB)     :: WWWS     
   REAL(SP), DIMENSION(0:NTNode,KB)     :: WWWSF   
   REAL(SP), DIMENSION(0:NTNode)        :: DTWWWS  
   REAL(SP), DIMENSION(0:NTNode,KB)     :: ZZZFLUX !! temporary total flux in corrected part
   REAL(SP), DIMENSION(0:NTNode,KB)     :: BETA    !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:NTNode,KB)     :: BETAIN  !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:NTNode,KB)     :: BETAOUT !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:NTNode,KB)     :: T1_FRESH    !! for source term which also bring mass volume

   INTEGER ITERA, NTERA
#  endif
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Fluxes-----------------------------------------------------------!
!
   XFLUX     = 0.0_SP
   XFLUX_ADV = 0.0_SP

!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
!!#  if !defined (WET_DRY)
   DO I=1,NCV
     I1=NTRG(I)
!     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
       UVN(I,K)  = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I)
#      if defined (SEMI_IMPLICIT)
       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
#      endif
     END DO
   END DO
!!#  else
!!   DO I=1,NCV
!!     I1=NTRG(I)
!!!     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
!!       UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!#      if defined (SEMI_IMPLICIT)
!!       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
!!       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
!!#      endif
!!     END DO
!!   END DO
!!#  endif

!
!--Add the Shortwave Radiation Body Force--------------------------------------!
!
     TTIME=THOUR
   IF(H_TYPE == 'body_h')THEN
     TTIME=THOUR
     DO  K=1,KBM1
       DO  I=1,NNode
         IF(TTIME < THOUR_HS)THEN
           RF(I,K)=0.0_SP
         ELSE
           IF(DT(I) > 0.0_SP) THEN
             ZDEP=0.5_SP*(Z(I,K)+Z(I,K+1))*DT(I)
           END IF
           RF(I,K)=-SWRAD(I)*((RHEAT/ZETA1)*EXP(ZDEP/ZETA1) &
                  +((1-RHEAT)/ZETA2)*EXP(ZDEP/ZETA2))*DT(I)
         END IF
       END DO
     END DO
   ELSE !! H_TYPE = 'flux_h'
     RF = 0.0_SP
   ENDIF

   DO  I=1,NNode
     IF(TTIME < THOUR_HS)THEN
       RF(I,1)=0.0_SP
     ELSE
       RF(I,1)=RF(I,1)-ROFVROS*(QEVAP3(I)-QPREC3(I))*T1(I,1)    
     END IF
   END DO
!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=1,KBM1
      PTPX  = 0.0_SP
      PTPY  = 0.0_SP
      PTPXD = 0.0_SP
      PTPYD = 0.0_SP
     DO I=1,NNode
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)
	 
#    if defined (WET_DRY)
         IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 1)THEN
          FFD=0.5_SP*(T1(I,K)+T1(I2,K)-TMEAN1(I,K)-TMEAN1(I2,K))
          FF1=0.5_SP*(T1(I,K)+T1(I2,K))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 1 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(T1(I1,K)+T1(I,K)-TMEAN1(I1,K)-TMEAN1(I,K))
          FF1=0.5_SP*(T1(I1,K)+T1(I,K))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(T1(I,K)+T1(I,K)-TMEAN1(I,K)-TMEAN1(I,K))
          FF1=0.5_SP*(T1(I,K)+T1(I,K))
	 ELSE
          FFD=0.5_SP*(T1(I1,K)+T1(I2,K)-TMEAN1(I1,K)-TMEAN1(I2,K))
          FF1=0.5_SP*(T1(I1,K)+T1(I2,K))
	 END IF 
#    else	 
         FFD=0.5_SP*(T1(I1,K)+T1(I2,K)-TMEAN1(I1,K)-TMEAN1(I2,K))
         FF1=0.5_SP*(T1(I1,K)+T1(I2,K))
#    endif	 
	 
#        if defined (SPHERICAL)
         XTMP  = VX(I2)*TPI-VX(I1)*TPI
         XTMP1 = VX(I2)-VX(I1)
         IF(XTMP1 >  180.0_SP)THEN
           XTMP = -360.0_SP*TPI+XTMP
         ELSE IF(XTMP1 < -180.0_SP)THEN
           XTMP =  360.0_SP*TPI+XTMP
         END IF  

!         TXPI=XTMP*COS(DEG2RAD*VY(I))
         TXPI=XTMP*VAL_COS_VY(I)
         TYPI=(VY(I1)-VY(I2))*tpi
#    if defined (NORTHPOLE)
         IF(NODE_NORTHAREA(I) == 1)THEN
           VX1_TMP = REARTH * COS(VY(I1)*DEG2RAD) * COS(VX(I1)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I1)*DEG2RAD))
           VY1_TMP = REARTH * COS(VY(I1)*DEG2RAD) * SIN(VX(I1)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I1)*DEG2RAD))

           VX2_TMP = REARTH * COS(VY(I2)*DEG2RAD) * COS(VX(I2)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I2)*DEG2RAD))
           VY2_TMP = REARTH * COS(VY(I2)*DEG2RAD) * SIN(VX(I2)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I2)*DEG2RAD))

           TXPI = (VX2_TMP-VX1_TMP)/(2._SP /(1._SP+SIN(VY(I)*DEG2RAD)))
           TYPI = (VY1_TMP-VY2_TMP)/(2._SP /(1._SP+SIN(VY(I)*DEG2RAD)))
  	   IF(I /= NODE_NORTHPOLE)THEN
	     TXPI_TMP = TYPI*COS(VX(I)*DEG2RAD)-TXPI*SIN(VX(I)*DEG2RAD)
	     TYPI_TMP = TXPI*COS(VX(I)*DEG2RAD)+TYPI*SIN(VX(I)*DEG2RAD)
	     TYPI_TMP = -TYPI_TMP
	    
	     TXPI = TXPI_TMP
	     TYPI = TYPI_TMP
	   END IF  
	 END IF 
#    endif	 
         PTPX(I)=PTPX(I)+FF1*TYPI
         PTPY(I)=PTPY(I)+FF1*TXPI
         PTPXD(I)=PTPXD(I)+FFD*TYPI
         PTPYD(I)=PTPYD(I)+FFD*TXPI
#        else
         PTPX(I)=PTPX(I)+FF1*(VY(I1)-VY(I2))
         PTPY(I)=PTPY(I)+FF1*(VX(I2)-VX(I1))
         PTPXD(I)=PTPXD(I)+FFD*(VY(I1)-VY(I2))
         PTPYD(I)=PTPYD(I)+FFD*(VX(I2)-VX(I1))
#        endif
       END DO
       PTPX(I)=PTPX(I)/ART2(I)
       PTPY(I)=PTPY(I)/ART2(I)
       PTPXD(I)=PTPXD(I)/ART2(I)
       PTPYD(I)=PTPYD(I)/ART2(I)
     END DO
          
     IF(K == KBM1)THEN
       DO I=1,NNode
         PFPXB(I) = PTPX(I)
         PFPYB(I) = PTPY(I)
       END DO
     END IF

     DO I=1,NNode
       VISCOFF(I) = VISCOFH(I,K)
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

!       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
!       XII=XCG2(I)
!       YII=YCG2(I)

       XI=XCG2(I)
       XTMP  = XI*TPI-VX(IA)*TPI
       XTMP1 = XI-VX(IA)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF 

!       DXA=XTMP*COS(DEG2RAD*VY(IA))  
       DXA=XTMP*VAL_COS_VY(IA)      
       DYA=(YI-VY(IA))*TPI

       XTMP  = XI*TPI-VX(IB)*TPI
       XTMP1 = XI-VX(IB)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF 

!       DXB=XTMP*COS(DEG2RAD*VY(IB))
       DXB=XTMP*VAL_COS_VY(IB)     
       DYB=(YI-VY(IB))*TPI
#      else
       DXA=XI-VX(IA)
       DYA=YI-VY(IA)
       DXB=XI-VX(IB)
       DYB=YI-VY(IB)
#      endif
       FIJ1=T1(IA,K)+DXA*PTPX(IA)+DYA*PTPY(IA)
       FIJ2=T1(IB,K)+DXB*PTPX(IB)+DYB*PTPY(IB)

       T1MIN=MINVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
       T1MIN=MIN(T1MIN, T1(IA,K))
       T1MAX=MAXVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
       T1MAX=MAX(T1MAX, T1(IA,K))
       T2MIN=MINVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
       T2MIN=MIN(T2MIN, T1(IB,K))
       T2MAX=MAXVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
       T2MAX=MAX(T2MAX, T1(IB,K))
       IF(FIJ1 < T1MIN) FIJ1=T1MIN
       IF(FIJ1 > T1MAX) FIJ1=T1MAX
       IF(FIJ2 < T2MIN) FIJ2=T2MIN
       IF(FIJ2 > T2MAX) FIJ2=T2MAX
    
       UN=UVN(I,K)
#      if defined (SEMI_IMPLICIT)
       UN1=UVN1(I,K)
#      endif

       VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

       TXX=0.5_SP*(PTPXD(IA)+PTPXD(IB))*VISCOF
       TYY=0.5_SP*(PTPYD(IA)+PTPYD(IB))*VISCOF

       FXX=-DTIJ(I,K)*TXX*DLTYE(I)
       FYY= DTIJ(I,K)*TYY*DLTXE(I)

#      if !defined (SEMI_IMPLICIT)
       EXFLUX=-UN*DTIJ(I,K)* &
           ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY
#      else
       EXFLUX=-UN*DTIJ(I,K)* &
           ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP
       EXFLUX=(1.0_SP-IFCETA)*EXFLUX+IFCETA*(-UN1*DTIJ1(I,K)*((1.0_SP+SIGN(1.0_SP,UN1))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN1))*FIJ1)*0.5_SP)+FXX+FYY
#      endif

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

       XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
       XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)

     END DO

#    if defined (SPHERICAL)  && (NORTHPOLE)
#    if !defined (SEMI_IMPLICIT)
       CALL ADV_T_XY(XFLUX,XFLUX_ADV,PTPX,PTPY,PTPXD,PTPYD,VISCOFF,K,0.0_SP)
#    else
       CALL ADV_T_XY(XFLUX,XFLUX_ADV,PTPX,PTPY,PTPXD,PTPYD,VISCOFF,K,IFCETA)
#    endif
#    endif

   END DO !! K LOOP

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
   END DO



!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
# if defined (MPDATA)
!--------------------------------------------------------------------------------
!   S. HU
!   Using smolarkiewicz, P. K; A fully multidimensional positive definite advection
!   transport algorithm with small implicit diffusion, Journal of Computational
!   Physics, 54, 325-362, 1984
!-----------------------------------------------------------------

!-----combine all the horizontal flux first-----------------------------------

!-------fresh water part-----------------------------------------------------
	T1_FRESH=T1

   IF(POINT_ST_TYPE == 'calculated') THEN
     IF(INFLOW_TYPE == 'node') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           JJ=INODEQ(J)
           STPOINT = TDIS(J)   !add by T.W. to correct T riv_BC bug. Sept. 2011
!           STPOINT=TDIS(J)  
           DO K=1,KBM1
!             STPOINT    = T1(JJ,K)   !comment out by T.W. to correct the riv_bc bug
             T1_FRESH(JJ,K)=TDIS(J)            !add by T.W. to correct riv_BC bug, Sept. 2011
!        T1_FRESH(JJ,K)=TDIS(J)            ! for gauge guess in the following
             XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT     !/DZ(JJ,K)
           END DO
         END DO
       END IF
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           J1=N_ICELLQ(J,1)
           J2=N_ICELLQ(J,2)
           DO K=1,KBM1
             STPOINT1 = T1(J1,K)
             STPOINT2 = T1(J2,K)
             XFLUX(J1,K)=XFLUX(J1,K)-  &
                         QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT1     !/DZ(J1,K)
             XFLUX(J2,K)=XFLUX(J2,K)-  &      
                         QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT2     !/DZ(J2,K)
           END DO
         END DO
       END IF
     END IF
   END IF
!------------------------------------------------------------------------------

! The horizontal term of advection is neglected here
   DO K=1,KBM1
     DO I=1,NNode
       IF(ISONB(I) == 2) THEN
         XFLUX(I,K)=0.
       ENDIF
     END DO
   END DO
   
! Initialize variables of MPDATA
   T1_S=0._SP
   T1_SF=0._SP
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
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*T1(I,K)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*T1(I,K+1) &
                  +(WTS(I,K)+ABS(WTS(I,K)))*T1(I,K)    
         ELSE IF(K == KBM1) THEN
           TEMP = +(WTS(I,K)-ABS(WTS(I,K)))*T1(I,K-1)     &
                  +(WTS(I,K)+ABS(WTS(I,K)))*T1(I,K)
         ELSE
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*T1(I,K)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*T1(I,K+1) &
                  +(WTS(I,K)-ABS(WTS(I,K)))*T1(I,K-1)     &
                  +(WTS(I,K)+ABS(WTS(I,K)))*T1(I,K)
         END IF
         TEMP = 0.5_SP*TEMP 

         IF(K == 1)THEN
           TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
           TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
           TMAX = MAX(TMAX,T1(I,K+1),T1(I,K),T1_FRESH(I,K))
           TMIN = MIN(TMIN,T1(I,K+1),T1(I,K),T1_FRESH(I,K))
         ELSEIF(K == KBM1) THEN
           TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
           TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
           TMAX = MAX(TMAX,T1(I,K-1),T1(I,K),T1_FRESH(I,K))
           TMIN = MIN(TMIN,T1(I,K-1),T1(I,K),T1_FRESH(I,K))
	ELSE
           TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
           TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
           TMAX = MAX(TMAX,T1(I,K+1),T1(I,K-1),T1(I,K),T1_FRESH(I,K))
           TMIN = MIN(TMIN,T1(I,K+1),T1(I,K-1),T1(I,K),T1_FRESH(I,K))
         END IF



         ZZZFLUX(I,K) = TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K)/ART1(I)*(DTI/DT(I))/DZ(I,K) 
         XXXX = ZZZFLUX(I,K)*DT(I)/DTFA(I)+T1(I,K)-T1(I,K)*DT(I)/DTFA(I) 

         BETA(I,K)=0.5*(1.-SIGN(1.,XXXX)) * (TMAX-T1(I,K))/(ABS(XXXX)+1.E-10) &
                  +0.5*(1.-SIGN(1.,-XXXX)) * (T1(I,K)-TMIN)/(ABS(XXXX)+1.E-10)

         T1_SF(I,K)=T1(I,K)-MIN(1.,BETA(I,K))*XXXX

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
       T1_S   = T1_SF
       DTWWWS = DT
     ELSE
       WWWSF  = WWWS
       T1_S   = T1_SF
       DTWWWS = DTFA
     END IF
     DO K=2,KBM1
       DO I=1,NNode
         TEMP=ABS(WWWSF(I,K))-DTI*(WWWSF(I,K))*(WWWSF(I,K))/DZ(I,K)/DTWWWS(I)
         WWWS(I,K)=TEMP*(T1_S(I,K-1)-T1_S(I,K))/(ABS(T1_S(I,K-1))+ABS(T1_S(I,K))+1.E-14)
 
         IF(TEMP < 0.0_SP .OR. T1_S(I,K) == 0.0_SP)THEN 
           WWWS(I,K)=0. 
         END IF
       END DO 
     END DO
     DO I=1,NNode
       WWWS(I,1)=0.
     END DO

     DO I=1,NNode
       TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),1))
       TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),1))
       TMAX = MAX(TMAX,T1(I,2),T1(I,1),T1_FRESH(I,1))
       TMIN = MIN(TMIN,T1(I,2),T1(I,1),T1_FRESH(I,1))
 
       TEMP=0.5*((WWWS(I,2)+ABS(WWWS(I,2)))*T1_S(I,2))*(DTI/DTFA(I))/DZ(I,1)
       BETAIN(I,1)=(TMAX-T1_S(I,1))/(TEMP+1.E-10)

       TEMP=0.5*((WWWS(I,1)+ABS(WWWS(I,1)))*T1_S(I,1)-        &
	           (WWWS(I,2)-ABS(WWWS(I,2)))*T1_S(I,1))*(DTI/DTFA(I))/DZ(I,1)
       BETAOUT(I,1)=(T1_S(I,1)-TMIN)/(TEMP+1.E-10)

       WWWSF(I,1)=0.5*MIN(1.,BETAOUT(I,1))*(WWWS(I,1)+ABS(WWWS(I,1))) + &
                    0.5*MIN(1.,BETAIN(I,1))*(WWWS(I,1)-ABS(WWWS(I,1)))
     END DO

     DO K=2,KBM1-1
       DO I=1,NNode
         TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
         TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
         TMAX = MAX(TMAX,T1(I,K+1),T1(I,K-1),T1(I,K),T1_FRESH(I,K))
         TMIN = MIN(TMIN,T1(I,K+1),T1(I,K-1),T1(I,K),T1_FRESH(I,K))
 
         TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*T1_S(I,K+1)-  &
	           (WWWS(I,K)-ABS(WWWS(I,K)))*T1_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
         BETAIN(I,K)=(TMAX-T1_S(I,K))/(TEMP+1.E-10)

         TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)-        &
	           (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*T1_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
         BETAOUT(I,K)=(T1_S(I,K)-TMIN)/(TEMP+1.E-10)

         WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
                    0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
       END DO
     END DO


      K=KBM1
       DO I=1,NNode
         TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
         TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
         TMAX = MAX(TMAX,T1(I,K-1),T1(I,K),T1_FRESH(I,K))
         TMIN = MIN(TMIN,T1(I,K-1),T1(I,K),T1_FRESH(I,K))
 
         TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*T1_S(I,K+1)-  &
	           (WWWS(I,K)-ABS(WWWS(I,K)))*T1_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
         BETAIN(I,K)=(TMAX-T1_S(I,K))/(TEMP+1.E-10)

         TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)-        &
	           (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*T1_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
         BETAOUT(I,K)=(T1_S(I,K)-TMIN)/(TEMP+1.E-10)

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
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*T1_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*T1_S(I,K+1) &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)
           ELSE IF(K == KBM1) THEN
             TEMP = +(WWWS(I,K)-ABS(WWWS(I,K)))*T1_S(I,K-1)     &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)
           ELSE
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*T1_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*T1_S(I,K+1) &
                    +(WWWS(I,K)-ABS(WWWS(I,K)))*T1_S(I,K-1)     &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)
           END IF
           TEMP = 0.5_SP*TEMP
           T1_SF(I,K)=(T1_S(I,K)-TEMP*(DTI/DTFA(I))/DZ(I,K)) 
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



# if !defined (MPDATA)
# if defined (ONE_D_MODEL)
    XFLUX = 0.0_SP
# endif    

!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------
   DO K=1,KBM1
     DO I=1,NNode
#    if defined (WET_DRY)
       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
#    endif
       IF(K == 1) THEN
         TEMP=-WTS(I,K+1)*(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
       ELSE IF(K == KBM1) THEN
         TEMP= WTS(I,K)*(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
       ELSE
         TEMP= WTS(I,K)*(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
               WTS(I,K+1)*(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
       END IF

       IF(ISONB(I) == 2) THEN
!         XFLUX(I,K)=TEMP*ART1(I)/DZ(I,K)
         XFLUX(I,K)=TEMP*ART1(I)
       ELSE
!         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/DZ(I,K)
         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)
       END IF
#    if defined (WET_DRY)
       END IF
#    endif
     END DO
   END DO  !!K LOOP
!
!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
   IF(POINT_ST_TYPE == 'calculated') THEN
     IF(INFLOW_TYPE == 'node') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           JJ=INODEQ(J)
           STPOINT=TDIS(J)   !add by T.W., to correct riv_bc bug, july 2011 
!           STPOINT=TDIS(J)  
           DO K=1,KBM1
!             STPOINT    = T1(JJ,K)  ! comment out by T.W., july 2011, to correct riv_bc bug, july 2011
!             XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT/DZ(JJ,K)
             XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT
           END DO
         END DO
       END IF
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           J1=N_ICELLQ(J,1)
           J2=N_ICELLQ(J,2)
           DO K=1,KBM1
             STPOINT1 = T1(J1,K)
             STPOINT2 = T1(J2,K)
!             XFLUX(J1,K)=XFLUX(J1,K)-  &
!                         QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT1/DZ(J1,K)
!             XFLUX(J2,K)=XFLUX(J2,K)-  &
!                         QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT2/DZ(J2,K)
             XFLUX(J1,K)=XFLUX(J1,K)-  &
                         QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT1
             XFLUX(J2,K)=XFLUX(J2,K)-  &
                         QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT2
           END DO
         END DO
       END IF
     END IF
   END IF
!---------------------------------------------------------------------

# endif

!---------------- heat flux by ground water------------------
   IF(IBFW > 0)THEN
     DO I=1,NNode
       DO J=1,IBFW
         IF(I == NODE_BFW(J))THEN
           XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS3(J)*BFWTDIS3(J)      !/DZ(I,KBM1)
         END IF
       END DO
     END DO
   END IF

!
!--Update Temperature----------------------------------------------------------!
!
   DO I=1,NNode
#  if defined (WET_DRY)
     IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1 )THEN
#  endif
     DO K=1,KBM1
       XFLUX(I,K) = XFLUX(I,K) - RF(I,K)*ART1(I)    !/DZ(I,K)
#    if !defined (MPDATA)     
!       TF1(I,K)=(T1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/DT(I)))*(DT(I)/DTFA(I))
       TF1(I,K)=(T1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#    else
       TF1(I,K)=T1_SF(I,K) 
!	TF1(I,K)=T1_SF(I,K) - RF(I,K)*DTI/DTFA(I)/DZ(I,K)  ! considering about heat flux source term, not tested yet
#    endif              
     END DO
#  if defined (WET_DRY)
     ELSE
     DO K=1,KBM1
       TF1(I,K)=T1(I,K)
     END DO
     END IF
#  endif
   END DO
    
   RETURN
   END SUBROUTINE ADV_T
!==============================================================================|
