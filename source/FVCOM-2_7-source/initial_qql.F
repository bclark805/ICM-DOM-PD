
!==============================================================================|
!   Initialize Turbulent Kinetic Energy and Length Scale                       |
!==============================================================================|

   SUBROUTINE INITIAL_QQL         

!------------------------------------------------------------------------------|

   USE ALL_VARS
   IMPLICIT NONE
   INTEGER :: I,K
!==============================================================================|
   
   DO K=1,KBM1
     DO I=1,NTNode
       IF(D(I) > 0.0_SP) AAM(I,K) = HORCON
     END DO
   END DO

!
!------------------------BOUNDARY VALUES---------------------------------------!
!

   DO I = 1, NTNode
     KM(I,1)   = 0.0_SP
     KM(I,KB)  = 0.0_SP
     KH(I,1)   = 0.0_SP
     KH(I,KB)  = 0.0_SP
     KQ(I,1)   = 0.0_SP
     KQ(I,KB)  = 0.0_SP
     L(I,1)    = 0.0_SP
     L(I,KB)   = 0.0_SP
     Q2(I,1)   = 0.0_SP
     Q2(I,KB)  = 0.0_SP
     Q2L(I,1)  = 0.0_SP
     Q2L(I,KB) = 0.0_SP
#    if defined (GOTM)
     TKE(I,1)  = 0.0_SP
    TEPS(I,1)  = 0.0_SP
     TKE(I,KB) = 0.0_SP
    TEPS(I,KB) = 0.0_SP
#    endif
   END DO


!
!------------------------INTERNAL VALUES---------------------------------------!
!
   DO  K = 2, KBM1
     DO I = 1, NTNode
       IF (D(I) > 0.0_SP) THEN
         Q2(I,K)  = 1.E-8
         Q2L(I,K) = 1.E-8
         L(I,K)   = 1.
#        if defined (GOTM)
         TKE(I,K) = Q2(I,K)/2.
         TEPS(I,K)= (.1704*TKE(I,K)**1.5)/L(I,K)
#        endif
         KM(I,K)  = 2.*UMOL
         KQ(I,K)  = 2.*UMOL
         KH(I,K)  = 2.*UMOL / VPRNU
       END IF
     END DO
   END DO

   RETURN
   END SUBROUTINE INITIAL_QQL
!==============================================================================|




