!==============================================================================|
! adiabatic temperature gradient deg c per decibar    			       |
! ref: bryden, h., 1973,deep-sea res.,20,401-408                               |
!                                                                              |
! units:                                                                       |
!       pressure        P4        decibars                                     |
!       temperature     T4        deg celsius(ipts-68)                         |
!       salinity        s4        (ipss-78)                                    |
!       adiabatic      atg        deg. c/decibar                               |
! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),                       |
! t=40 deg c,p0=10000 decibars                                                 |
!==============================================================================|

   REAL FUNCTION ATG(S4,T4,P4)
   USE MOD_PREC

!------------------------------------------------------------------------------|

   IMPLICIT NONE
   REAL(SP), INTENT(IN) :: S4,T4,P4
   REAL(SP)  :: DS

!==============================================================================|
   
   DS  = S4 - 35.0_SP
   ATG = (((-2.1687e-16_SP*T4+1.8676e-14_SP)*T4-4.6206e-13_SP)*P4 &
        +((2.7759e-12_SP*T4-1.1351e-10_SP)*DS+((-5.4481e-14_SP*T4 &
        +8.733e-12_SP)*T4-6.7795e-10_SP)*T4+1.8741e-8_SP))*P4 &
        +(-4.2393e-8_SP*T4+1.8932e-6_SP)*DS &
        +((6.6228e-10_SP*T4-6.836e-8_SP)*T4+8.5258e-6_SP)*T4+3.5803e-5_SP

   RETURN
   END FUNCTION ATG
!==============================================================================|
