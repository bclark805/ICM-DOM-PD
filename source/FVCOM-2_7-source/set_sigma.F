!==============================================================================|
! This program is used to set up the coordinate in the vertical.               !
!								               !
! case(1): sigma levels                                                        !
! sigma levels are determined by a formula of                                  !
!                      sigma(k)=-[(k-1)/(kb-1)]^k11                            !
!    p_sigma=1: uniform sigma layers                                           !
!    p_sigma=2: layers satisfying a parabolic function with high               !
!               vertical resolution near the surface and bottom.               !
!    p_sigma can be used any real number                                       !
!									       !
! case(2): general vertical level                                              !
! vertical levels are determined by the formula                                !
!                tanh[(dl+du)((kbm1-k)/kbm1)-dl]+tanh(dl)                      !
!        z(k)= ------------------------------------------  - 1                 !
!                      tanh(dl) + tanh(du)                                     !
!                                                                              !
! case(3): constant layer transformation                                       !
! four values need to be specified:                                            !
!  DUU the upper boundaries, up to which the co-ordinates are parallel must be !
!      defined.                                                                !
!  DLL the lower boundaries, up to which the co-ordinates are parallel must be !
!      defined.                                                                !
!  HMIN1 the minimum water depth at which the layers are constant. If H < HMIN1!
!      then sigma co-ordinates are used.                                       !
!                                                                              !
! Reference of case(2), case(3) and case(4):                                   !
! Pietrzak, J.D., Jan B. Jakobson, Hans Burchard, Hans Jacob Vested, Ole       !
! Petersen , 2002. A three-dimensional hydrostatic model for coastal and ocean !
! modelling using a generalised topography following co-ordinate system. Ocean !
! Modelling 4, 173-205                                                         !
!                                                                              !
!  calculates: z(NNode,kb) vertical levels					       !
!  calculates: dz(NNode,kb-1) delta between vertical levels		               !
!  calculates: zz(NNode,kb-1) intra-vertical levels				       !
!  calculates: dzz(NNode,kb-2) delta between intra-vertical levels		       !
!==============================================================================|

   SUBROUTINE SET_SIGMA(IDX)           

!==============================================================================|
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif

!TW, added to update velocity block
#  if defined (V_BLOCK)
   USE VELOCITY_BLOCK
#  endif

   IMPLICIT NONE
   INTEGER :: K,KK
   INTEGER :: IDX,I
   REAL(SP):: ZTMP(KB)
   REAL(SP):: X1,X2,X3   !,DL,DU
!   REAL(SP):: HMIN1,DUU,DLL
   REAL(SP):: DR,RCL,RCU
!   INTEGER :: KU,KL
!   REAL(SP):: ZKU(KB),ZKL(KB)
    REAL(SP):: P_SIGMA_SHAL,P_SIGMA_DEP  ! for Hybrid sigma system B Clark 2015 HMIN2,
!==============================================================================|

   SELECT CASE(IDX)

   CASE(1)

   IF(P_SIGMA > 1 .AND. MOD(KB,2) == 0)THEN
     IF(MSR) PRINT*, 'kb shoude be an odd number,stop ....'
     CALL PSTOP
   END IF

!--------  SET SIGMA LEVELS  --------------------------------------------------!  
!orginal formula to set sigma   
!   IF(P_SIGMA == 1)THEN
!     DO K=1,KB
!       ZTMP(K) = -((K-1)/FLOAT(KB-1))**P_SIGMA 
!     END DO
!   ELSE
!     DO K=1,(KB+1)/2
!       ZTMP(K) = -((K-1)/FLOAT((KB+1)/2-1))**P_SIGMA/2 
!     END DO
!     DO K=(KB+1)/2+1,KB
!       ZTMP(K) = ((KB-K)/FLOAT((KB+1)/2-1))**P_SIGMA/2-1.0
!     END DO
!   END IF

!=======================B Clark, updated hybrid sigma system ===========================  

!HMIN2=1.E6

P_SIGMA_SHAL=1.5
P_SIGMA_DEP=1.5

      DO K=1,KB
          ZTMP(K) = -((K-1)/FLOAT(KB-1))**P_SIGMA 
        END DO
  DO I = 1,NNode
     IF (H(I) <= HMIN2) THEN
 	 	
        DO K=1,KB
          ZTMP(K) = -((K-1)/FLOAT(KB-1))**P_SIGMA_SHAL 
        END DO
   
     ELSE
 	 
        DO K=1,(KB+1)/2
 	       ZTMP(K) = -((K-1)/FLOAT((KB+1)/2-1))**P_SIGMA_DEP/2 
        END DO
 
        DO K=(KB+1)/2+1,KB
          ZTMP(K) = ((KB-K)/FLOAT((KB+1)/2-1))**P_SIGMA_DEP/2-1.0
        END DO
 	
      END IF

     DO K=1,KB
       Z(I,K)=ZTMP(K)
!update sigma if block activated, TW
#if defined V_BLOCK
       if(BLOCK) then
         Z(I,K) = ZTMP(K) * (1.0 + dep_block(I)) + dep_block(I)
       end if
#endif
!end update sigma, TW
    
! Added by B Clark for hybrid sigma output variable
 
    S_HYBRID(I,K)=Z(I,K)

!! end Update Hybrid Sigma B Clark 2015

     END DO
  END DO    

   DO I=1,MElem
     DO K=1,KB
       Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0_SP
     END DO
   END DO    

   CASE(2)
   
   Z=0.0;Z1=0.0
   DO K=1,KBM1
     X1=DL2+DU2
     X1=X1*(KBM1-K)/KBM1
     X1=X1-DL2
     X1=TANH(X1)
     X2=TANH(DL2)
     X3=X2+TANH(DU2)
     DO I=1,NNode
       Z(I,K+1)=(X1+X2)/X3-1.0_SP
     END DO
     DO I=1,MElem
       Z1(I,K+1)=(X1+X2)/X3-1.0_SP
     END DO
   END DO

   CASE(3)

   DO I=1,NNode
     IF(H(I) < HMIN1)THEN
       Z(I,1)=0.0
       DL2=0.001;DU2=0.001
       DO K=1,KBM1
         X1=DL2+DU2
         X1=X1*(KBM1-K)/KBM1
         X1=X1-DL2
         X1=TANH(X1)
         X2=TANH(DL2)
         X3=X2+TANH(DU2)

         Z(I,K+1)=(X1+X2)/X3-1.0_SP
       END DO
     ELSE
       DR=(H(I)-DUU-DLL)/H(I)/(KB-KU-KL-1)

       Z(I,1)=0.0_SP

       DO K=2,KU+1
         Z(I,K)=Z(I,K-1)-ZKU(K-1)/H(I)
       END DO

       DO K=KU+2,KB-KL
         Z(I,K)=Z(I,K-1)-DR
       END DO

       KK=0
       DO K=KB-KL+1,KB
         KK=KK+1
         Z(I,K)=Z(I,K-1)-ZKL(KK)/H(I)
       END DO
     END IF
   END DO
  S_HYBRID=Z
 
   DO I=1,MElem
     DO K=1,KB
       Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0_SP
     END DO
   END DO    

   CASE(4)

   DO I=1,NNode
     IF(H(I) < HMIN1)THEN
       RCU=-DUU/HMIN1
       RCL=DLL/HMIN1-1
       DR=(RCL-RCU)/(KB-KU-KL-1)

       DO K=1,KU
         ZKU(K)=RCU/KU
       END DO
       DO K=1,KL
         ZKL(K)=(-1.0_SP-RCL)/KL
       END DO

       Z(I,1)=0.0_SP

       DO K=2,KU+1
         Z(I,K)=Z(I,K-1)+ZKU(K-1)
       END DO

       DO K=KU+2,KB-KL
         Z(I,K)=Z(I,K-1)+DR
       END DO

       KK=0
       DO K=KB-KL+1,KB
         KK=KK+1
         Z(I,K)=Z(I,K-1)+ZKL(KK)
       END DO

     ELSE
       DR=(H(I)-DUU-DLL)/H(I)/(KB-KU-KL-1)

       Z(I,1)=0.0_SP

       DO K=2,KU+1
         Z(I,K)=Z(I,K-1)-ZKU(K-1)/H(I)
       END DO

       DO K=KU+2,KB-KL
         Z(I,K)=Z(I,K-1)-DR
       END DO

       KK=0
       DO K=KB-KL+1,KB
         KK=KK+1
         Z(I,K)=Z(I,K-1)-ZKL(KK)/H(I)
       END DO
     END IF
   END DO
   
   DO I=1,MElem
     DO K=1,KB
       Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0_SP
     END DO
   END DO    

   END SELECT    

!---------COMPUTE SIGMA DERIVATIVES AND INTRA SIGMA LEVELS---------------------!
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,Z1)
     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,Z)
#    endif
   
   DO K=1,KB-1
     DO I=1,NTNode
       DZ(I,K)  = Z(I,K)-Z(I,K+1)
       ZZ(I,K)  = .5_SP*(Z(I,K)+Z(I,K+1))
     END DO
     DO I=1,MTElem
       DZ1(I,K)  = Z1(I,K)-Z1(I,K+1)
       ZZ1(I,K)  = .5_SP*(Z1(I,K)+Z1(I,K+1))
     END DO
   END DO

   DO I=1,NTNode
     ZZ(I,KB) = 2.0_SP*ZZ(I,KB-1)-ZZ(I,KB-2)
   END DO
   DO I=1,MTElem
     ZZ1(I,KB) = 2.0_SP*ZZ1(I,KB-1)-ZZ1(I,KB-2)
   END DO

   DO K=1,KBM2
     DO I=1,NTNode
       DZZ(I,K) = ZZ(I,K)-ZZ(I,K+1)
     END DO
     DO I=1,MTElem
       DZZ1(I,K) = ZZ1(I,K)-ZZ1(I,K+1)
     END DO
   END DO
   
   DZZ(:,KBM1) = 0.0_SP
   DZ(:,KB)    = 0.0_SP
   DZZ1(:,KBM1) = 0.0_SP
   DZ1(:,KB)    = 0.0_SP


  
!----------OUTPUT VALUES-TO INFOFILE-------------------------------------------!

   IF(MSR)THEN
     WRITE(IPT,*  )'!'
     WRITE(IPT,*  )'!'
     WRITE(IPT,*)'!                SIGMA LAYER INFO     '
     WRITE(IPT,70)
     SELECT CASE(IDX)
     CASE(1)
     DO K=1,KB
       WRITE(IPT,80) K,Z(1,K),ZZ(1,K),DZ(1,K),DZZ(1,K)
     END DO
     CASE(2)
     DO K=1,KB
       WRITE(IPT,80) K,Z(1,K),ZZ(1,K),DZ(1,K),DZZ(1,K)
     END DO
     CASE(3)
     DO K=1,KB
!       WRITE(IPT,80) K,Z(30762,K),ZZ(30762,K),DZ(30762,K),DZZ(30762,K)
     END DO
     DO K=1,KB
!       WRITE(IPT,80) K,Z(564,K),ZZ(564,K),DZ(564,K),DZZ(564,K)
     END DO
     CASE(4)
     DO K=1,KB
!       WRITE(IPT,80) K,Z(1,K),ZZ(1,K),DZ(1,K),DZZ(1,K)
     END DO
     DO K=1,KB
!       WRITE(IPT,80) K,Z(158,K),ZZ(158,K),DZ(158,K),DZZ(158,K)
     END DO
     END SELECT    
     WRITE(IPT,*  )'!'
   END IF

!----------FORMAT STATEMENTS---------------------------------------------------!

70 FORMAT(2x,'k',13x,'z',11x,'zz',11x,'dz',11x,'dzz')
80 FORMAT(' ',i5,4f13.8)

   RETURN
   END SUBROUTINE SET_SIGMA
!==============================================================================|
