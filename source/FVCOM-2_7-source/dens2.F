!==============================================================================|
!     COMPUTE DENSITY USING SALINITY AND POTENTIAL TEMP                        |
!									       | 
!  CALCULATES: RHO1(NNode) DENSITY AT NODES				               |
!  CALCULATES: RHO (MElem) DENSITY AT ELEMENTS				       |
!==============================================================================|

   SUBROUTINE DENS2               

!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   REAL(SP), DIMENSION(NNode,KB) :: RHOF,TF,SF
   INTEGER :: I,K
!==============================================================================|

!
!  CALCULATE DENSITY FROM EQUATION OF STATE
!
   DO I=1,NNode
     DO K=1,KBM1
!       TF(I,K) = T1(I,K)
        SF(I,K) = S1(I,K)

       TF(I,K) = max(T1(I,K),1.0) !T.K 1/31/2013 added this for stability in intertidal regions

       RHOF(I,K) = SF(I,K) * SF(I,K) * SF(I,K) * &
           6.76786136E-6_SP - SF(I,K) * SF(I,K) * 4.8249614E-4_SP + &
           SF(I,K) * 8.14876577E-1_SP - 0.22584586E0_SP

       RHOF(I,K) = RHOF(I,K) * (TF(I,K)*TF(I,K)* &
           TF(I,K)*1.667E-8_SP-TF(I,K)*TF(I,K)*8.164E-7_SP+ &
           TF(I,K)*1.803E-5_SP)

       RHOF(I,K) = RHOF(I,K) + 1. - TF(I,K) * TF(I,K) * &
           TF(I,K) * 1.0843E-6_SP + TF(I,K) * TF(I,K) * &
           9.8185E-5_SP - TF(I,K) * 4.786E-3_SP

       RHOF(I,K) = RHOF(I,K) * (SF(I,K)*SF(I,K)* &
           SF(I,K)*6.76786136E-6_SP-SF(I,K)*SF(I,K)* &
           4.8249614E-4_SP+SF(I,K)*8.14876577E-1_SP+3.895414E-2_SP)

       RHOF(I,K) = RHOF(I,K) - (TF(I,K)-3.98_SP) ** 2 * ( &
           TF(I,K)+283.0_SP) / (503.57_SP*(TF(I,K)+67.26_SP))
     END DO
   END DO

!
!  CALCULATE RHO1
!
   DO I=1,NNode
     IF (D(I) > 0.0_SP)THEN
       DO K=1,KBM1
         RHO1(I,K) =  RHOF(I,K)*1.e-3_SP
       END DO
     END IF
   END DO

!
!  AVERAGE FROM NODES TO FACE CENTERS
!
   CALL N2E3D(RHO1,RHO)

   RETURN
   END SUBROUTINE DENS2
!==============================================================================|
