!==============================================================================|
   FUNCTION THETA(S4,T04,P04,PR)
!==============================================================================|
! to compute local potential temperature at pr using                           |
! bryden 1973 polynomial for adiabatic lapse rate and                          |
! runge-kutta 4th order integration algorithm.                                 |
! ref: bryden,h.,1973,deep-sea res.,20,401-408;                                |
! fofonoff,n.,1977,deep-sea res.,24,489-491                                    |
!                                                                              |
! units:                                                                       |
!       pressure        p04       decibars                                     |
!       temperature     t04       deg celsius (ipts-68)                        |
!       salinity         s4        (ipss-78)                                   |
!       reference prs    pr       decibars                                     |
!       potential tmp.  theta     deg celsius                                  |
! checkvalue:                                                                  |
!             theta= 36.89073 c,s=40 (ipss-78),                                |
!             t0=40 deg c,p0=10000 decibars,pr=0 decibars                      |
!                                                                              |
!                                                                              |
! set up intermediate temperature and pressure variables.                      |
!==============================================================================|

   USE MOD_PREC
   IMPLICIT NONE
   REAL(SP) :: THETA
   REAL(SP), INTENT(IN) :: S4,T04,P04,PR
   REAL(SP) :: P4,T4,H4,Q4,XK
   REAL(SP) :: ATG

!==============================================================================|

   P4 = P04
   T4 = T04
   H4 = PR - P4
   XK = H4*ATG(S4,T4,P4)
   T4 = T4 + .5_SP*XK
   Q4 = XK
   P4 = P4 + 0.5_SP*H4
   XK = H4*ATG(S4,T4,P4)
   T4 = T4 + 0.29289322_SP*(XK-Q4)
   Q4 = 0.58578644_SP*XK + 0.121320344_SP*Q4
   XK = H4*ATG(S4,T4,P4)
   T4 = T4 + 1.707106781_SP*(XK-Q4)
   Q4 = 3.414213562_SP*XK - 4.121320344_SP*Q4
   P4 = P4 + 0.5_SP*H4
   XK = H4*ATG(S4,T4,P4)
   THETA = T4 + (XK-2.0_SP*Q4)/6.0_SP

   RETURN
   END FUNCTION THETA
!==============================================================================|
