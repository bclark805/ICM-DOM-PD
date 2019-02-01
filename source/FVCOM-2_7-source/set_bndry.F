!==============================================================================|
!  SET METRICS FOR THE BOUNDARY CONDITIONS       			       |
!==============================================================================|

   SUBROUTINE SET_BNDRY       

!------------------------------------------------------------------------------!
   USE ALL_VARS
   USE BCS
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
   USE MOD_OBCS

   IMPLICIT NONE
   REAL(DP)  DX12,DY12,DX32,DY32,ATMP1,ATMP2,DXYTMP,HTMP,AREATMP
   REAL(DP)  XNORM,YNORM,XP,YP,XN,YN,XI,YI,FAC,XNEXT,YNEXT,MODNR
   INTEGER I,J,I1,I2,I3,J1,J2,II,ITMP,JTMP,INODE,JNODE,KNODE,NNORM
# if defined (SPHERICAL)
   REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,XDIS1,YDIS1,XDIS2,YDIS2,SIDE
# endif
!------------------------------------------------------------------------------!


   IF(NUMQBC > 0)THEN
     IF(INFLOW_TYPE == 'node')THEN
       DO I=1,NUMQBC
         J=INODEQ(I)
         I1=NBSN(J,2)
         I2=NBSN(J,1)
         I3=NBSN(J,NTSN(J)-1)
         DY12=VY(I1)-VY(I2)
         DY32=VY(I3)-VY(I2)
#        if defined (SPHERICAL)
         X1_DP = VX(I2)
         Y1_DP = VY(I2)
         X2_DP = VX(I1)
         Y2_DP = VY(I1)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         DX12 = SIDE

         X2_DP = VX(I3)
         Y2_DP = VY(I3)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         DX32 = SIDE

         DY12 = TPI*DY12
         DY32 = TPI*DY32
#        else
         DX12=VX(I1)-VX(I2)
         DX32=VX(I3)-VX(I2)
#        endif

         ATMP1=ATAN2(DY12,DX12)
         ATMP2=ATAN2(DY32,DX32)
!         IF(ATMP1 < 0.0_SP) ATMP1=ATMP1+2.0_SP*3.1415927_SP
!         IF(ATMP2 < 0.0_SP) ATMP2=ATMP2+2.0_SP*3.1415927_SP
         IF(ATMP1 < ATMP2) ATMP1=ATMP1+2.0_SP*3.1415927_SP
         DXYTMP=SQRT(DX12**2+DY12**2)+SQRT(DX32**2+DY32**2)
         QAREA(I)=0.5_SP*DXYTMP*H(INODEQ(I))
         ANGLEQ(I)=(ATMP1-ATMP2)/2.+ATMP2
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO I=1,NUMQBC
         II=ICELLQ(I)
         IF(ISBCE(II) /= 1) THEN
            PRINT*, 'NO.',I,'FRESHW INFLOW CELL'
            PRINT*, 'IS NOT A SLLID BOUNDARY ONE'
            CALL PSTOP
         END IF
         ITMP=0
         DO J=1,3
            IF(NBE(II,J) == 0) THEN
              JTMP=J
              ITMP=ITMP+1
            END IF
         END DO
         IF(ITMP /= 1) THEN
            PRINT*, 'NO SOLID BOUNDARY OR MORE THAN ONE SOLID BOUNDARY'
            PRINT*, 'IN NO.',I,'FRESHW INFLOW CELL'
            CALL PSTOP
         END IF
         J1=JTMP+1-INT((JTMP+1)/4)*3
         J2=JTMP+2-INT((JTMP+2)/4)*3
         I1=NV(II,J1)
         I2=NV(II,J2)
         N_ICELLQ(I,1)=I1
         N_ICELLQ(I,2)=I2
         HTMP=0.5_SP*(H(I1)+H(I2))
         DY12=VY(I1)-VY(I2)
#        if defined (SPHERICAL)
         X1_DP = VX(I2)
         Y1_DP = VY(I2)
         X2_DP = VX(I1)
         Y2_DP = VY(I1)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         DX12 = SIDE

         DY12 = TPI*DY12
#        else
         DX12=VX(I1)-VX(I2)
#        endif
         ATMP1=ATAN2(DY12,DX12)
         QAREA(I)=SQRT(DX12**2+DY12**2)*HTMP
         ANGLEQ(I)=ATMP1+3.1415927/2.0
         AREATMP=ART1(I1)+ART1(I2)
         RDISQ(I,1)=ART1(I1)/AREATMP
         RDISQ(I,2)=ART1(I2)/AREATMP
       END DO
     ELSE
       PRINT*, 'INFLOW_TYPE NOT CORRECT!!'
       CALL PSTOP
     END IF
   END IF

   RETURN
   END SUBROUTINE SET_BNDRY
!==============================================================================|
