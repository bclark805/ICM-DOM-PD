# if !defined (SEMI_IMPLICIT)
!==============================================================================|
!     ACCUMLATE FLUXES FOR EXTERNAL MODE                                       |
!==============================================================================|

   SUBROUTINE EXTUV_EDGE(K)       

!==============================================================================|
   USE ALL_VARS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif

#  if defined (SPHERICAL)
#  if defined (NORTHPOLE)
   USE MOD_NORTHPOLE
#  endif   
#  endif   

#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
#  endif

!TW, update to include velocity block
#  if defined (V_BLOCK)
   USE VELOCITY_BLOCK
#  endif

#  if defined (KELP)
   USE MOD_KELP
#  endif
!! finished for Kelp module, T.W, April 2013


   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:MTElem) :: RESX,RESY,TMP
   REAL(SP) :: UAFT,VAFT
   INTEGER  :: I

!==============================================================================|

!
!--ACCUMULATE RESIDUALS FOR EXTERNAL MODE EQUATIONS----------------------------|
!
   UAFT = UAF(0)
   VAFT = VAF(0)
   H1(0)= H1(1)

!!#  if defined (WET_DRY)
!!   IF(K == 3)THEN
!!     RESX = ADX2D + ADVUA + DRX2D + PSTX - COR*VA*D1*ART  &
!!            -(WUSURF2 + WUBOT)*ART
!!     RESY = ADY2D + ADVVA + DRY2D + PSTY + COR*UA*D1*ART  &
!!            -(WVSURF2 + WVBOT)*ART
!!#  if defined (SPHERICAL)
!!     RESX = RESX -UA*VA/REARTH*TAN(DEG2RAD*YC)*D1*ART
!!     RESY = RESY +UA*UA/REARTH*TAN(DEG2RAD*YC)*D1*ART
!!#  endif

!!!
!!!--UPDATE----------------------------------------------------------------------|
!!!

!!     UAF = (UARK*(H1+ELRK1)-ALPHA_RK(K)*DTE*RESX/ART)/(H1+ELF1)
!!     VAF = (VARK*(H1+ELRK1)-ALPHA_RK(K)*DTE*RESY/ART)/(H1+ELF1)
!!     UAS = UAF
!!     VAS = VAF
!!   END IF
!!#  endif

   DO I=1,MTElem
#  if defined (WET_DRY)
     IF(ISWET_CELL_LAST_EXT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif
       RESX(I) = ADX2D(I)+ADVUA(I)+DRX2D(I)+PSTX(I)-COR(I)*VA(I)*D1(I)*ART(I) &
#  if defined V_BLOC	   
					* (1.0 - bwt(i)) &    !depth adjust added by TW for block
#  endif						
                 -(WUSURF2(I)+WUBOT(I))*ART(I)
       RESY(I) = ADY2D(I)+ADVVA(I)+DRY2D(I)+PSTY(I)+COR(I)*UA(I)*D1(I)*ART(I) &
#  if defined V_BLOCK	   
					* (1.0 - bwt(i)) &
#  endif						
                 -(WVSURF2(I)+WVBOT(I))*ART(I)

!update for MHK device momentum removal, by T.W.

#  if defined (KELP)
      IF(C_KELP) THEN
        RESX(I) = RESX(I) - SUM(EMS_X(I,:))
        RESY(I) = RESY(I) - SUM(EMS_Y(I,:))
      END IF
#  endif
! finished addition, T.W.


#  if defined (SPHERICAL)
       RESX(I) = RESX(I)     &
                 -UA(I)*VA(I)/REARTH*TAN(DEG2RAD*YC(I))*D1(I)*ART(I)
       RESY(I) = RESY(I)     &
                 +UA(I)*UA(I)/REARTH*TAN(DEG2RAD*YC(I))*D1(I)*ART(I)
#  endif
!
!--UPDATE----------------------------------------------------------------------|
!

       UAF(I)  =  (UARK(I)*(H1(I)+ELRK1(I))-ALPHA_RK(K)*DTE*RESX(I)/ART(I))/(H1(I)+ELF1(I))
       VAF(I)  =  (VARK(I)*(H1(I)+ELRK1(I))-ALPHA_RK(K)*DTE*RESY(I)/ART(I))/(H1(I)+ELF1(I))
#  if defined (WET_DRY)
     ELSE
       UAF(I) = 0.0_SP
       VAF(I) = 0.0_SP
     END IF
#  endif
   END DO

#  if defined (SPHERICAL)
#  if defined (NORTHPOLE)
   CALL EXTUV_EDGE_XY(K)
#  endif
#  endif
   
   VAF(0) = VAFT
   UAF(0) = UAFT

!
!--ADJUST EXTERNAL VELOCITY IN SPONGE REGION-----------------------------------|
!
   UAF = UAF-CC_SPONGE*UAF
   VAF = VAF-CC_SPONGE*VAF

!
!--STORE VARIABLES FOR MOMENTUM BALANCE CHECK----------------------------------|
!
#  if defined (BALANCE_2D)
   IF(K == 4) THEN
     TMP=ART*(H1+ELF1)
     ADVUA2 = ADVUA2 + (ADVUA-ADFXA)/TMP/FLOAT(ISPLIT)    !X- HORIZONTAL ADVECTION (m/s^2)
     ADVVA2 = ADVVA2 + (ADVVA-ADFYA)/TMP/FLOAT(ISPLIT)    !X- HORIZONTAL ADVECTION (m/s^2)
     ADFX2  = ADFX2 + ADFXA/TMP/FLOAT(ISPLIT)
     ADFY2  = ADFY2 + ADFYA/TMP/FLOAT(ISPLIT)
     DRX2D2 = DRX2D2 + DRX2D/TMP/FLOAT(ISPLIT)            !X- BAROCLINIC PRESURE GRADIENT FORCE
     DRY2D2 = DRY2D2 + DRY2D/TMP/FLOAT(ISPLIT)            !Y- BAROCLINIC PRESURE GRADIENT FORCE
     CORX2  = CORX2  - COR*VA*D1/(H1+ELF1)/FLOAT(ISPLIT)  !X- CORIOLIS FORCE
     CORY2  = CORY2  + COR*UA*D1/(H1+ELF1)/FLOAT(ISPLIT)  !Y- CORIOLIS FORCE
     PSTX2  = PSTX2  + PSTX/TMP/FLOAT(ISPLIT)             !X- BAROTROPIC PRESURE GRADIENT FORCE
     PSTY2  = PSTY2  + PSTY/TMP/FLOAT(ISPLIT)             !Y- BAROTROPIC PRESURE GRADIENT FORCE
     ADX2D2 = ADX2D2 + ADX2D/TMP/FLOAT(ISPLIT)            !GX (m/s^2)
     ADY2D2 = ADY2D2 + ADY2D/TMP/FLOAT(ISPLIT)            !GY (m/s^2)
     WUSURBF2=WUSURBF2-(WUSURF2+WUBOT)/(H1+ELF1)/FLOAT(ISPLIT) !X-SURFACE & BOTTOM FRICTION
     WVSURBF2=WVSURBF2-(WVSURF2+WVBOT)/(H1+ELF1)/FLOAT(ISPLIT) !Y-SURFACE & BOTTOM FRICTION
     DUDT2  = DUDT2 + (UAF-UARK*(H1+ELRK1)/(H1+ELF1))/DTE/FLOAT(ISPLIT)
     DVDT2  = DVDT2 + (VAF-VARK*(H1+ELRK1)/(H1+ELF1))/DTE/FLOAT(ISPLIT)
     
   END IF     
#  endif

   RETURN
   END SUBROUTINE EXTUV_EDGE
!==============================================================================|
# endif
