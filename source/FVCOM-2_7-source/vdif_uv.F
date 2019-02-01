!==============================================================================|
! this subroutine uses an implicit method to calculate the vertical            !
! diffusion term in the x and y momentum equations                             !
!==============================================================================|

   SUBROUTINE VDIF_UV             

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   INTEGER I,K,ITMP1,ITMP2,ITMP3,KI
   REAL(SP), DIMENSION(0:MElem,KB) :: C,A,VHU,VHPU,VHV,VHPV
   REAL(SP) :: TMP

   C = KM1

   DO K = 2, KBM1
     DO I = 1, MElem
#  if !defined (WET_DRY)
       IF (D1(I) > 0.0_SP) THEN
#  else
       IF(ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif

         A(I,K-1) = -DTI *(C(I,K)+UMOL) / (DZ1(I,K-1)*DZZ1(I,K-1)* D1(I)*D1(I))
         C(I,K) = -DTI *(C(I,K)+UMOL) / (DZ1(I,K)*DZZ1(I,K-1)* D1(I)*D1(I))
       END IF
     END DO
   END DO

   DO I = 1, MElem
#  if !defined (WET_DRY)
     IF (D1(I) > 0.0_SP) THEN
#  else
     IF(ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif
       ITMP1=ISONB(NV(I,1))
       ITMP2=ISONB(NV(I,2))
       ITMP3=ISONB(NV(I,3))
       IF(ITMP1 == 2 .OR. ITMP2 == 2 .OR. ITMP3 == 2) THEN
         TMP=0.0_SP
       ELSE
         TMP=1.0_SP
       END IF
       VHU(I,1) = A(I,1) / (A(I,1)-1.)
       VHV(I,1) = A(I,1) / (A(I,1)-1.)
       VHPU(I,1) = (-DTI*WUSURF(I)*TMP/(-DZ1(I,1)*D1(I))-UF(I,1)) / (A(I,1)-1.)
       VHPV(I,1) = (-DTI*WVSURF(I)*TMP/(-DZ1(I,1)*D1(I))-VF(I,1)) / (A(I,1)-1.)
     END IF
   END DO

   DO  K = 2, KBM2
     DO I = 1, MElem
#  if !defined (WET_DRY)
       IF (D1(I) > 0.0_SP) THEN
#  else
       IF(ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif
         VHPU(I,K) = 1.0_SP / (A(I,K)+C(I,K)*(1.-VHU(I,K-1))-1.)
         VHPV(I,K) = 1.0_SP / (A(I,K)+C(I,K)*(1.-VHV(I,K-1))-1.)
         VHU(I,K)  = A(I,K) * VHPU(I,K)
         VHV(I,K)  = A(I,K) * VHPV(I,K)
         VHPU(I,K) = (C(I,K)*VHPU(I,K-1)-UF(I,K))*VHPU(I,K)
         VHPV(I,K) = (C(I,K)*VHPV(I,K-1)-VF(I,K))*VHPV(I,K)
       END IF
     END DO
   END DO

   DO  I = 1, MElem
#  if !defined (WET_DRY)
     IF (D1(I) > 0.0_SP) THEN
#  else
     IF(ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif
       TPS(I) = CBC(I) * SQRT(U(I,KBM1)**2+V(I,KBM1)**2)
       UF(I,KBM1) = (C(I,KBM1)*VHPU(I,KBM2)-UF(I,KBM1))/ &
       (TPS(I)*DTI/(-DZ1(I,KBM1)*D1(I))-1.-(VHU(I,KBM2)-1.)*C(I,KBM1))
       VF(I,KBM1) = (C(I,KBM1)*VHPV(I,KBM2)-VF(I,KBM1))/ &
       (TPS(I)*DTI/(-DZ1(I,KBM1)*D1(I))-1.-(VHV(I,KBM2)-1.)*C(I,KBM1))
     END IF
   END DO

   DO  K = 2, KBM1
     KI = KB - K
     DO  I = 1, MElem
#  if !defined (WET_DRY)
       IF (D1(I) > 0.0_SP) THEN
#  else
       IF(ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif
         UF(I,KI) = (VHU(I,KI)*UF(I,KI+1)+VHPU(I,KI))
         VF(I,KI) = (VHV(I,KI)*VF(I,KI+1)+VHPV(I,KI))
#  if defined (WET_DRY)
       ELSE
         UF(I,KI) = 0.0_SP
         VF(I,KI) = 0.0_SP
#  endif
       END IF
     END DO
   END DO

!
!--Damp Velocity In Sponge Region----------------------------------------------!
!
   DO K = 1, KBM1
     DO I = 1, MElem
       UF(I,K) = UF(I,K)-CC_SPONGE(I)*UF(I,K)
       VF(I,K) = VF(I,K)-CC_SPONGE(I)*VF(I,K)
     END DO
   END DO

   DO  I = 1, MElem
#  if !defined (WET_DRY)
     IF (D1(I) > 0.0_SP) THEN
#  else
     IF(ISWET_CELL_CURRENTSTEP(I) == 1)THEN
#  endif
       WUBOT(I) = -TPS(I) * U(I,KBM1)
       WVBOT(I) = -TPS(I) * V(I,KBM1)
#  if defined (WET_DRY)
     ELSE
       WUBOT(I) = 0.0_SP
       WVBOT(I) = 0.0_SP
#  endif
     END IF
   END DO

   RETURN
   END SUBROUTINE VDIF_UV
!==============================================================================|
