   SUBROUTINE ADCOR

   USE ALL_VARS
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
   REAL(SP) :: UFC(0:MTElem,KB),VFC(0:MTElem,KB)
   REAL(SP),PARAMETER :: BETA0=0.5
   REAL(SP) ::CURCOR,PRECOR
   INTEGER :: I,K  
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP  
#  if defined (SEMI_IMPLICIT)
#  if defined (SPHERICAL) && (NORTHPOLE)
   REAL(SP) :: UFC_NP(0:MTElem,KB),VFC_NP(0:MTElem,KB)
   REAL(SP) :: UU_TMP1, VV_TMP1, UU_TMP2, VV_TMP2
#  endif
#  endif

#  if !defined (TWO_D_MODEL)
   UFC=0.0
   VFC=0.0
#  if defined (SEMI_IMPLICIT)
#  if defined (SPHERICAL) && (NORTHPOLE)
   UFC_NP=0.0
   VFC_NP=0.0
#  endif
#  endif
#  endif

   DO I = 1, MElem
#  if defined (SEMI_IMPLICIT) && (TWO_D_MODEL)
     CURCOR=BETA0*COR(I)*VAF(I)
     PRECOR=(1.0_SP-BETA0)*COR(I)*VA(I)
     ADVUA(I)=UBETA2D(I)-(CURCOR+PRECOR)*DT1(I)*ART(I)*EPOR(I)
#  else

#  if defined (SPHERICAL) && (NORTHPOLE)
     IF(CELL_NORTHAREA(I) == 1)THEN
       DO K = 1, KBM1
#        if !defined (SEMI_IMPLICIT)
         V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)

         CURCOR=BETA0*COR(I)*VF(I,K)
         PRECOR=(1.-BETA0)*COR(I)*V_TMP
         UFC(I,K)=UBETA(I,K)-(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)
#        else
         UU_TMP2 = UF(I,K)
         VV_TMP2 = VF(I,K)
         UU_TMP1 = VV_TMP2*COS(XC(I)*DEG2RAD)-UU_TMP2*SIN(XC(I)*DEG2RAD) 
         VV_TMP1 = -( UU_TMP2*COS(XC(I)*DEG2RAD)+VV_TMP2*SIN(XC(I)*DEG2RAD) )
         CURCOR=BETA0*COR(I)*VV_TMP1
         PRECOR=(1.-BETA0)*COR(I)*V(I,K)
         UFC(I,K)=UBETA(I,K)-(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)

         V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)
         CURCOR=BETA0*COR(I)*VF(I,K)
         PRECOR=(1.-BETA0)*COR(I)*V_TMP
         UFC_NP(I,K)=UBETA_NP(I,K)-(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)
#        endif
       END DO
     ELSE
#  endif   
       DO K = 1, KBM1
         CURCOR=BETA0*COR(I)*VF(I,K)
         PRECOR=(1.-BETA0)*COR(I)*V(I,K)
         UFC(I,K)=UBETA(I,K)-(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
       END DO
#  if defined (SPHERICAL) && (NORTHPOLE)
     END IF
#  endif        
   
#  endif   
   END DO

   DO I = 1, MElem
#  if defined (SEMI_IMPLICIT) && (TWO_D_MODEL)
     CURCOR=BETA0*COR(I)*UAF(I)
     PRECOR=(1.0_SP-BETA0)*COR(I)*UA(I)
     ADVVA(I)=VBETA2D(I)+(CURCOR+PRECOR)*DT1(I)*ART(I)*EPOR(I)
#  else

#  if defined (SPHERICAL) && (NORTHPOLE)
     IF(CELL_NORTHAREA(I) == 1)THEN
       DO K = 1, KBM1
#        if !defined (SEMI_IMPLICIT)
         U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)

         CURCOR=BETA0*COR(I)*UF(I,K)
         PRECOR=(1.-BETA0)*COR(I)*U_TMP
         VFC(I,K)=VBETA(I,K)+(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)
#        else
         UU_TMP2 = UF(I,K)
         VV_TMP2 = VF(I,K)
         UU_TMP1 = VV_TMP2*COS(XC(I)*DEG2RAD)-UU_TMP2*SIN(XC(I)*DEG2RAD)
         VV_TMP1 = -( UU_TMP2*COS(XC(I)*DEG2RAD)+VV_TMP2*SIN(XC(I)*DEG2RAD) )
         CURCOR=BETA0*COR(I)*UU_TMP1
         PRECOR=(1.-BETA0)*COR(I)*U(I,K)
         VFC(I,K)=VBETA(I,K)+(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)

         U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
         CURCOR=BETA0*COR(I)*UF(I,K)
         PRECOR=(1.-BETA0)*COR(I)*U_TMP
         VFC_NP(I,K)=VBETA_NP(I,K)+(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)
#        endif
       END DO
     ELSE
#  endif   
       DO K = 1, KBM1
         CURCOR=BETA0*COR(I)*UF(I,K)
         PRECOR=(1.-BETA0)*COR(I)*U(I,K)
         VFC(I,K)=VBETA(I,K)+(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
       END DO
#  if defined (SPHERICAL) && (NORTHPOLE)
     END IF  
#  endif        

#  endif
   END DO

#  if !defined (TWO_D_MODEL)
   DO I=1,MElem
#  if defined (SPHERICAL) && (NORTHPOLE)
     IF(CELL_NORTHAREA(I) == 1)THEN
       DO K=1,KBM1
#        if !defined (SEMI_IMPLICIT)
         U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
         V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)

         UF_TMP=U_TMP*DT1(I)/D1(I)-DTI*UFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         VF_TMP=V_TMP*DT1(I)/D1(I)-DTI*VFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))

         UF(I,K)  = VF_TMP*COS(XC(I)*DEG2RAD)-UF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = UF_TMP*COS(XC(I)*DEG2RAD)+VF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = -VF(I,K)    
#        else
         XFLUX3(I,K) = UFC(I,K)
         YFLUX3(I,K) = VFC(I,K)
         XFLUX3_NP(I,K) = UFC_NP(I,K)
         YFLUX3_NP(I,K) = VFC_NP(I,K)
#        endif
       END DO
     ELSE
#  endif   
       DO K=1,KBM1
#        if !defined (SEMI_IMPLICIT)
         UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*UFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*VFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
#        else
         XFLUX3(I,K) = UFC(I,K)
         YFLUX3(I,K) = VFC(I,K) 
#        endif
       END DO
#  if defined (SPHERICAL) && (NORTHPOLE)
     END IF  
#  endif        
   END DO
#  endif

   RETURN
   END SUBROUTINE ADCOR
