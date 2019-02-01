!==============================================================================|
!  ADJUST THE VERTICAL WATER COLUMN WHEN UNSTABLE                              |
!==============================================================================|

   SUBROUTINE CONV_OVER           

!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   REAL(SP):: AVE_T,AVE_S,AVE_R 
   INTEGER :: I,K,KK,J1,J2,J3
!==============================================================================|

!--APPROXIMATE CONVECTIVE OVERTURNING------------------------------------------!

   DO I=1,NNode
     DO K=KBM1,2,-1
       DO KK=K-1,1,-1
         IF(RHO1(I,K) < RHO1(I,KK)) THEN
           AVE_T = SUM(  T1(I,KK:K))/FLOAT(K-KK+1) 
           AVE_S = SUM(  S1(I,KK:K))/FLOAT(K-KK+1) 
           AVE_R = SUM(RHO1(I,KK:K))/FLOAT(K-KK+1) 
           T1(I,KK:K)   = AVE_T
           S1(I,KK:K)   = AVE_S
           RHO1(I,KK:K) = AVE_R
         END IF
       END DO
     END DO
   END DO
       
!-----RECALCULATE ELEMENT-BASED VALUES OF SALINITY/TEMP/DENSITY----------------!

   DO I=1,MElem
     J1=NV(I,1) ; J2 = NV(I,2) ; J3 = NV(I,3)
     DO K=1,KBM1
       T(I,K)  = ONE_THIRD*(  T1(J1,K)+  T1(J2,K)+  T1(J3,K))
       S(I,K)  = ONE_THIRD*(  S1(J1,K)+  S1(J2,K)+  S1(J3,K))
     END DO
   END DO   
          
       
   RETURN
   END SUBROUTINE CONV_OVER       
!==============================================================================|

