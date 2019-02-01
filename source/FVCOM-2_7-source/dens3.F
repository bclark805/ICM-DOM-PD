!==============================================================================|
!     COMPUTE IN SITU DENSITY - 1000  USING SALINITY, POTENTIAL TEMP,          |
!     AND PRESSURE FROM A POLYNOMIAL EXPRESSION (JACKETT & MCDOUGALL,          |
!     1995). IT ASSUMES  NO  PRESSURE  VARIATION  ALONG GEOPOTENTIAL           |
!     SURFACES, THAT IS, DEPTH (METERS; NEGATIVE) AND PRESSURE (DBAR           |
!     ASSUMED NEGATIVE HERE) ARE INTERCHANGEABLE.                              |
!                                                                              |
!     check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                             |
!        RHOF  = 1050.3639165364     (kg/m3)                                   |
!        DEN1  = 1028.2845117925     (kg/m3)                                   |
!                                                                              |
!  Reference:                                                                  |
!                                                                              |
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of             |
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.            |
!    and Oceanic Techn., vol. 12, pp. 381-389.                                 |
!									       | 
!    CALCULATES: RHO1(NNode) DENSITY AT NODES		                       |
!    CALCULATES: RHO (MElem) DENSITY AT ELEMENTS				       |
!==============================================================================|

   SUBROUTINE DENS3        

!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   REAL(SP), DIMENSION(NNode,KB) :: RHOF, DEN1
   REAL(SP) :: TF,SF,sqrtSF,PZ,PBAR,TEMP(10),BULK,BULK0,BULK1,BULK2
   INTEGER :: I,K
!==============================================================================|
!  Polynomial  expansion  coefficients for the computation of in situ          |
!  density  via  the  nonlinear  equation of state  for seawater as a          |
!  function of potential temperature, salinity, and pressure (Jackett          |
!  and McDougall, 1995).                                                       |
   REAL(SP), PARAMETER :: A00 = +1.965933e+04_SP
   REAL(SP), PARAMETER :: A01 = +1.444304e+02_SP
   REAL(SP), PARAMETER :: A02 = -1.706103e+00_SP
   REAL(SP), PARAMETER :: A03 = +9.648704e-03_SP
   REAL(SP), PARAMETER :: A04 = -4.190253e-05_SP
   REAL(SP), PARAMETER :: B00 = +5.284855e+01_SP
   REAL(SP), PARAMETER :: B01 = -3.101089e-01_SP
   REAL(SP), PARAMETER :: B02 = +6.283263e-03_SP
   REAL(SP), PARAMETER :: B03 = -5.084188e-05_SP
   REAL(SP), PARAMETER :: D00 = +3.886640e-01_SP
   REAL(SP), PARAMETER :: D01 = +9.085835e-03_SP
   REAL(SP), PARAMETER :: D02 = -4.619924e-04_SP
   REAL(SP), PARAMETER :: E00 = +3.186519e+00_SP
   REAL(SP), PARAMETER :: E01 = +2.212276e-02_SP
   REAL(SP), PARAMETER :: E02 = -2.984642e-04_SP
   REAL(SP), PARAMETER :: E03 = +1.956415e-06_SP
   REAL(SP), PARAMETER :: F00 = +6.704388e-03_SP
   REAL(SP), PARAMETER :: F01 = -1.847318e-04_SP
   REAL(SP), PARAMETER :: F02 = +2.059331e-07_SP
   REAL(SP), PARAMETER :: G00 = +1.480266e-04_SP
   REAL(SP), PARAMETER :: G01 = +2.102898e-04_SP
   REAL(SP), PARAMETER :: G02 = -1.202016e-05_SP
   REAL(SP), PARAMETER :: G03 = +1.394680e-07_SP
   REAL(SP), PARAMETER :: H00 = -2.040237e-06_SP
   REAL(SP), PARAMETER :: H01 = +6.128773e-08_SP
   REAL(SP), PARAMETER :: H02 = +6.207323e-10_SP

   REAL(SP), PARAMETER :: Q00 = +9.99842594e+02_SP
   REAL(SP), PARAMETER :: Q01 = +6.793952e-02_SP
   REAL(SP), PARAMETER :: Q02 = -9.095290e-03_SP
   REAL(SP), PARAMETER :: Q03 = +1.001685e-04_SP
   REAL(SP), PARAMETER :: Q04 = -1.120083e-06_SP
   REAL(SP), PARAMETER :: Q05 = +6.536332e-09_SP
   REAL(SP), PARAMETER :: U00 = +8.24493e-01_SP
   REAL(SP), PARAMETER :: U01 = -4.08990e-03_SP
   REAL(SP), PARAMETER :: U02 = +7.64380e-05_SP
   REAL(SP), PARAMETER :: U03 = -8.24670e-07_SP
   REAL(SP), PARAMETER :: U04 = +5.38750e-09_SP
   REAL(SP), PARAMETER :: V00 = -5.72466e-03_SP
   REAL(SP), PARAMETER :: V01 = +1.02270e-04_SP
   REAL(SP), PARAMETER :: V02 = -1.65460e-06_SP
   REAL(SP), PARAMETER :: W00 = +4.8314e-04_SP

!
!  CALCULATE DENSITY FROM EQUATION OF STATE
!
   DO I=1,NNode
     DO K=1,KBM1
       TF = T1(I,K)
       SF = S1(I,K)
       sqrtSF = sqrt(SF)
       PZ = -ZZ(I,K)*D(I)
       PBAR = GRAV_N(I)*1.025_SP*PZ*0.01_SP

!  Compute density (kg/m3) at standard one atmosphere pressure
       TEMP(1)=Q00+TF*(Q01+TF*(Q02+TF*(Q03+TF*(Q04+TF*Q05))))
       TEMP(2)=U00+TF*(U01+TF*(U02+TF*(U03+TF*U04)))
       TEMP(3)=V00+TF*(V01+TF*V02)
       DEN1(I,K)=TEMP(1)+SF*(TEMP(2)+sqrtSF*TEMP(3)+SF*W00)

!  Compute secant bulk modulus (BULK = BULK0 + BULK1*PBAR + BULK2*PBAR*PBAR)
       TEMP(4)=A00+TF*(A01+TF*(A02+TF*(A03+TF*A04)))
       TEMP(5)=B00+TF*(B01+TF*(B02+TF*B03))
       TEMP(6)=D00+TF*(D01+TF*D02)
       TEMP(7)=E00+TF*(E01+TF*(E02+TF*E03))
       TEMP(8)=F00+TF*(F01+TF*F02)
       TEMP(9)=G01+TF*(G02+TF*G03)
       TEMP(10)=H00+TF*(H01+TF*H02)

       BULK0=TEMP(4)+SF*(TEMP(5)+sqrtSF*TEMP(6))
       BULK1=TEMP(7)+SF*(TEMP(8)+sqrtSF*G00)
       BULK2=TEMP(9)+SF*TEMP(10)
       BULK = BULK0 + PBAR * (BULK1 + PBAR * BULK2)

!  Compute "in situ" density anomaly (kg/m3)
       RHOF(I,K)=(DEN1(I,K)*BULK)/(BULK-PBAR)
       RHOF(I,K)= RHOF(I,K)-1000.0_SP
     END DO
   END DO

!
!  CALCULATE RHO1
!
   DO I=1,NNode
     DO K=1,KBM1
       RHO1(I,K) =  RHOF(I,K)*1.e-3_SP
     END DO
   END DO

!
!  AVERAGE FROM NODES TO FACE CENTERS
!
   CALL N2E3D(RHO1,RHO)

   RETURN
   END SUBROUTINE DENS3
!==============================================================================|
