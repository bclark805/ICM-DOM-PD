!==============================================================================!
!==============================================================================!
   FUNCTION SVAN(S4,T4,P04,SIGMA)
!==============================================================================!
! specific volume anomaly (steric anomaly) based on 1980 equation              |
! of state for seawater and 1978 practerical salinity scale.                   |
! references:                                                                  |
! millero, et al (1980) deep-sea res.,27a,255-264                              |
! millero and poisson 1981,deep-sea res.,28a pp 625-629.                       |
! both above references are also found in unesco report 38 (1981)              |
!                                                                              |
! units:                                                                       |
!       pressure        p04       decibars                                     |
!       temperature     t4        deg celsius (ipts-68)                        |
!       salinity        s4        (ipss-78)                                    |
!       spec. vol. ana. svan     m**3/kg *1.0e-8                               |
!       density ana.    sigma    kg/m**3                                       |
!                                                                              |
! check value: svan=981.3021 e-8 m**3/kg. for s = 40 (ipss-78),                |
! t = 40 deg c, p0= 10000 decibars.                                            |
! check value: sigma = 59.82037  kg/m**3. for s = 40 (ipss-78) ,               |
! t = 40 deg c, p0= 10000 decibars.                                            |
!==============================================================================!


   USE MOD_PREC
   IMPLICIT NONE
   REAL(SP) :: SVAN
   REAL(SP), INTENT(IN)  :: S4,T4,P04
   REAL(SP), INTENT(OUT) :: SIGMA
   REAL(SP) P4,SIG,SR,RR1,RR2,RR3,V350P,DK
   REAL(SP) A4,B4,C4,D4,E4,AA1,BB1,AW,BW,KK,K0,KW,K35,SVA
   REAL(SP) GAM,PK,DVAN,DR35P

   REAL(SP), PARAMETER :: R3500 = 1028.1063_SP
   REAL(SP), PARAMETER :: RR4   = 4.8314E-4_SP
   REAL(SP), PARAMETER :: DR350 = 28.106331_SP


!   rr4 is refered to as  c  in millero and poisson 1981
! convert pressure to bars and take square root salinity.

   P4=P04/10.0_SP
   SR = SQRT(ABS(S4))

! pure water density at atmospheric pressure
!   bigg p.h.,(1967) br. j. applied physics 8 pp 521-537.
!

   RR1=((((6.536332E-9_SP*T4-1.120083E-6_SP)*T4+1.001685E-4_SP)*T4 &
          -9.095290E-3_SP)*T4+6.793952E-2_SP)*T4-28.263737_SP


! seawater density atm press.
!  coefficients involving salinity
!  rr2 = a   in notation of millero and poisson 1981

   RR2=(((5.3875E-9_SP*T4-8.2467E-7_SP)*T4+7.6438E-5_SP)*T4-4.0899E-3_SP)*T4 &
          +8.24493E-1_SP

!  rr3 = b4  in notation of millero and poisson 1981

   RR3=(-1.6546E-6_SP*T4+1.0227E-4_SP)*T4-5.72466E-3_SP

!  international one-atmosphere equation of state of seawater

   SIG=(RR4*S4+RR3*SR+RR2)*S4+RR1

! specific volume at atmospheric pressure

   V350P = 1.0_SP/R3500
   SVA = -SIG*V350P/(R3500+SIG)
   SIGMA=SIG+DR350

!  scale specific vol. anamoly to normally reported units

   SVAN=SVA*1.0E+8_SP
   IF(P4 == 0.0_SP) RETURN

!-------------------------------------------------------------|
!    new high pressure equation of sate for seawater          |
!                                                             |
!        millero, el al., 1980 dsr 27a, pp 255-264            |
!        constant notation follows article                    |
!-------------------------------------------------------------|
! compute compression terms

   E4  = (9.1697E-10*T4+2.0816E-8_SP)*T4-9.9348E-7_SP
   BW  = (5.2787E-8_SP*T4-6.12293E-6_SP)*T4+3.47718E-5_SP
   B4  = BW + E4*S4
 
   D4  = 1.91075E-4
   C4  = (-1.6078E-6_SP*T4-1.0981E-5_SP)*T4+2.2838E-3_SP
   AW  = ((-5.77905E-7_SP*T4+1.16092E-4_SP)*T4+1.43713E-3_SP)*T4 &
          -0.1194975_SP
   A4  = (D4*SR + C4)*S4 + AW
 
   BB1 = (-5.3009E-4_SP*T4+1.6483E-2_SP)*T4+7.944E-2_SP
   AA1 = ((-6.1670E-5_SP*T4+1.09987E-2_SP)*T4-0.603459_SP)*T4+54.6746
   KW  = (((-5.155288E-5_SP*T4+1.360477E-2_SP)*T4-2.327105_SP)*T4 &
            +148.4206_SP)*T4-1930.06_SP
   K0  = (BB1*SR + AA1)*S4 + KW

! evaluate pressure polynomial
!-----------------------------------------------------|
!   k equals the secant bulk modulus of seawater      |
!   dk=k(s,t,p)-k(35,0,p)                             |
!   k35=k(35,0,p)                                     |
!-----------------------------------------------------|

   DK = (B4*P4 + A4)*P4 + K0
   K35  = (5.03217E-5_SP*P4+3.359406_SP)*P4+21582.27_SP
   GAM=P4/K35
   PK = 1.0_SP - GAM
   SVA = SVA*PK + (V350P+SVA)*P4*DK/(K35*(K35+DK))

!  scale specific vol. anamoly to normally reported units

   SVAN=SVA*1.0E+8_SP
   V350P = V350P*PK

!----------------------------------------------------------|
! compute density anamoly with respect to 1000.0 kg/m**3   |
!  1) dr350: density anamoly at 35 (ipss-78),              |
!                               0 deg. c and 0 decibars    |
!  2) dr35p: density anamoly at 35 (ipss-78),              |
!                               0 deg. c, pres. variation  |
!  3) dvan : density anamoly variations involving specific |
!            volume anamoly                                |
!                                                          |
! check values: sigma = 59.82037 kg/m**3                   |
! for s = 40 (ipss-78), t = 40 deg c, p0= 10000 decibars.  |
!----------------------------------------------------------|

   DR35P=GAM/V350P
   DVAN=SVA/(V350P*(V350P+SVA))
   SIGMA=DR350+DR35P-DVAN

   RETURN
   END FUNCTION SVAN
!==============================================================================!
!==============================================================================!
