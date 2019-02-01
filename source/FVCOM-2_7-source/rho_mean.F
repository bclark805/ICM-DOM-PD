!==============================================================================|

   SUBROUTINE RHO_MEAN
!==============================================================================|
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)  
   USE MOD_PAR
#  endif
   IMPLICIT NONE

   INTEGER, PARAMETER :: KBB=601
   INTEGER, PARAMETER :: KBBM1=KBB-1

   REAL(SP) DELTZ
   REAL(SP), DIMENSION(KBBM1)   :: PHY_Z   !Depth(m) in every standary Z levels 
   REAL(SP), DIMENSION(KBBM1)   :: RHOZTMP !density in standary Z levels 
   REAL(SP), DIMENSION(NNode,KBBM1) :: RHOZ    !density in standary Z levals,in Nodes 
   REAL(SP), DIMENSION(KBBM1)   :: RHOA    !density mean in standary Z levals
   
   REAL(SP), DIMENSION(NNodeGL)   :: HTMP
   REAL(SP), DIMENSION(NNodeGL,KBBM1)   :: RHOZZTMP
   
   REAL(SP), DIMENSION(KBM1)    :: ZM      !Depth (m) in every sigma levels for giving node
   REAL(SP), DIMENSION(KBM1)    :: RHOS    !density for giving node
   INTEGER ICOUNT,I,K,IERR

!--CALCULATE Z-LEVELS TO MAX DEPTH---------------------------------------------|
   
   DELTZ=HMAX/FLOAT(KBBM1)

   DO K=1,KBBM1
     PHY_Z(K)=(0.5_SP-FLOAT(K))*DELTZ
   END DO

!--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT Z LEVELS-------------------|
   
   DO I=1,NNode
     DO K=1,KBM1
       ZM(K)=ZZ(I,K)*DT(I)+EL(I)
       RHOS(K)= RHO1(I,K)                   !=DBLE(RHO1(I,K))
     END DO
     
     CALL SINTER_P(ZM,RHOS,PHY_Z,RHOZTMP,KBM1,KBBM1)

     DO K=1,KBBM1
       RHOZ(I,K)=RHOZTMP(K)
     END DO
   END DO

!--DO THE AVERAGE OVER Z_levels 

   IF(SERIAL)THEN
     RHOA = 0.0_SP
     DO K=1,KBBM1
       ICOUNT=0
       DO I=1,NNode
         IF(-H(I) <= PHY_Z(K))THEN
           ICOUNT = ICOUNT + 1
           RHOA(K)=RHOA(K)+RHOZ(I,K)
         END IF
       ENDDO

       IF(ICOUNT >= 1)THEN
         RHOA(K) = RHOA(K)/FLOAT(ICOUNT)
       ELSE 
         RHOA(K) = RHOA(K-1) 
       END IF
     END DO
   END IF  

#  if defined(MULTIPROCESSOR)
   IF(PAR)THEN
     CALL GATHER(LBOUND(RHOZ,1),UBOUND(RHOZ,1),        &
                  NNode,NNodeGL,KBBM1,MYID,NPROCS,NMAP,RHOZ,RHOZZTMP)
     CALL GATHER(LBOUND(H,1),UBOUND(H,1),        &
                  NNode,NNodeGL,1,MYID,NPROCS,NMAP,H,HTMP)
     
     IF(MSR)THEN
       RHOA = 0.0_SP
       DO K=1,KBBM1
         ICOUNT = 0.0_SP
         DO I=1,NNodeGL
           IF(-HTMP(I) <= PHY_Z(K))THEN
             ICOUNT = ICOUNT + 1.0_SP
             RHOA(K)=RHOA(K)+RHOZZTMP(I,K)
           END IF
         ENDDO

         IF(ICOUNT >= 1)THEN
           RHOA(K) = RHOA(K)/FLOAT(ICOUNT)
         ELSE 
           RHOA(K) = RHOA(K-1) 
         END IF
       END DO
     END IF  
     
     CALL MPI_BCAST(RHOA,KBBM1,MPI_F,0,MPI_COMM_WORLD,IERR)
   END IF  
#  endif

!--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT SIGMA LEVELS-------------------|

   DO I=1,NNode
     DO K=1,KBM1
       ZM(K)=ZZ(I,K)*DT(I)+EL(I)
     END DO

     CALL SINTER(PHY_Z,RHOA,ZM,RHOS,KBBM1,KBM1)
     RMEAN1(I,1:KBM1)=RHOS
   END DO
   RMEAN1(:,KB) = 0.0_SP
   
   CALL N2E3D(RMEAN1,RMEAN)
   
   RETURN
   END
   
