!==============================================================================|
!     COMPUTE DENSITY USING SALINITY AND POTENTIAL TEMP                        |
!									       | 
!  CALCULATES: RHO1(NNode) DENSITY AT NODES				               |
!  CALCULATES: RHO (MElem) DENSITY AT ELEMENTS				       |
!==============================================================================|

   SUBROUTINE DENS2               

!==============================================================================|
   !USE ALL_VARS
   USE MOD_LIMS
   !USE ADV_WQM
   USE MOD_WQM, ONLY: T, SALT
   USE MOD_HYDROVARS, ONLY: D
   USE MOD_PREC, ONLY : SP
   USE MOD_SED_DOM, ONLY : RHO1
   
   IMPLICIT NONE
   REAL(SP), DIMENSION(MTLOC,KBM1) :: RHOF,TF,SF
   
   INTEGER :: I,K
!==============================================================================|

!
!  CALCULATE DENSITY FROM EQUATION OF STATE
!
   DO I=1,MLOC
     DO K=1,KBM1
!       TF(I,K) = T1(I,K)
        SF(I,K) = SALT(I,K)

       TF(I,K) = max(T(I,K),1.0) !T.K 1/31/2013 added this for stability in intertidal regions

       RHOF(I,K) = SF(I,K) * SF(I,K) * SF(I,K) * &
           6.76786136E-6 - SF(I,K) * SF(I,K) * 4.8249614E-4 + &
           SF(I,K) * 8.14876577E-1 - 0.22584586E0

       RHOF(I,K) = RHOF(I,K) * (TF(I,K)*TF(I,K)* &
           TF(I,K)*1.667E-8-TF(I,K)*TF(I,K)*8.164E-7+ &
           TF(I,K)*1.803E-5)

       RHOF(I,K) = RHOF(I,K) + 1. - TF(I,K) * TF(I,K) * &
           TF(I,K) * 1.0843E-6 + TF(I,K) * TF(I,K) * &
           9.8185E-5 - TF(I,K) * 4.786E-3

       RHOF(I,K) = RHOF(I,K) * (SF(I,K)*SF(I,K)* &
           SF(I,K)*6.76786136E-6-SF(I,K)*SF(I,K)* &
           4.8249614E-4+SF(I,K)*8.14876577E-1+3.895414E-2)

       RHOF(I,K) = RHOF(I,K) - (TF(I,K)-3.98_SP) ** 2 * ( &
           TF(I,K)+283.0) / (503.57*(TF(I,K)+67.26))
     END DO
   END DO

!
!  CALCULATE RHO1
!
   DO I=1,MLOC
     IF (D(I) > 0.0)THEN
       DO K=1,KBM1
         RHO1(I,K) =  RHOF(I,K)+1000.!*1.e-3
       END DO
     END IF
   END DO

!
!  AVERAGE FROM NODES TO FACE CENTERS
!
   !CALL N2E3D(RHO1,RHO)

   RETURN
   END SUBROUTINE DENS2
!==============================================================================|
