!==============================================================================|
!   Calculate Potential Density Based on Potential Temp and Salinity           |
!     Pressure effects are incorported (Can Model Fresh Water < 4 Deg C)       |
!     Ref:  algorithms for computation of fundamental properties of            |
!         seawater , Fofonoff and Millard.				       |
!                                                                              |
!  calculates: rho1(nnode) density at nodes				       |
!  calculates: rho (ncell) density at elements				       |
!==============================================================================|

   SUBROUTINE DENS                

!------------------------------------------------------------------------------|

   USE ALL_VARS 
   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP)  :: PR,RZU,PT,SVA,SIGMA
   REAL(SP) :: THETA,SVAN
   

!==============================================================================|
   
   PR = 0.0_SP

   DO K=1,KBM1
     DO I=1,NNode    
       IF (D(I)>0.0_SP)THEN
!       RZU = GRAV*1.025_SP*(ZZ(I,K)*D(I))*0.1_SP
       RZU = -GRAV_N(I)*1.025_SP*(ZZ(I,K)*D(I))*0.1_SP
       PT  = THETA(S1(I,K),T1(I,K),RZU,PR)
       SVA = SVAN(S1(I,K),PT,RZU,SIGMA)
       RHO1(I,K) = SIGMA*1.e-3_SP
       END IF
     END DO
   END DO

!----------------TRANSFORM TO FACE CENTER--------------------------------------

   CALL N2E3D(RHO1,RHO)

   RETURN
   END SUBROUTINE DENS
!==============================================================================|
