
!==============================================================================|
!    Initialize Temperature, Salinity, Density Fields (T1,T,S1,S,RHO1,RHO)     !
!    Calculate Mean Fields (RMEAN,RMEAN1,TMEAN,TMEAN1,SMEAN,SMEAN1)            !
!==============================================================================|

   SUBROUTINE INITIAL_TS

!==============================================================================!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
!# if defined (MULTIPROCESSOR)
!   include "mpif.h"
!# endif
   REAL(SP), DIMENSION(KBM1)  :: ZM
   REAL(SP), DIMENSION(KB )   :: TI,SI  !initial temperature and salinity
   REAL(SP), DIMENSION(KSL)   :: TA,SA
   REAL(SP), DIMENSION(NNode,KSL) :: TS1,SS1 !specified temperature and salinity at KSL levels on nodes
!   REAL(SP), DIMENSION(MElem,KSL) :: TS,SS   !temperature and salinity at KSL levels on elements (not used)
   CHARACTER(LEN=80) :: COMT
   REAL(SP) :: T0,S0,T11,T22,S11,S22,DMAX,DMIN,SCOUNT,SCTOT,FAC,TATOT,SATOT
   REAL(SP), ALLOCATABLE :: TEMP1(:,:),TEMP2(:,:)
   CHARACTER(LEN=10) TSTYPE
   INTEGER :: I,II,IW,J1,J2,J3,K,IB,IERR

!==============================================================================!



!==============================================================================|
!    READ IN TEMP/SALINITY PARAMETERS FROM: 'casename_its.dat'                 !
!==============================================================================|

!
!  READ HEADER INFORMATION
!

   READ(INITS,*)!HEADER
   READ(INITS,'(A10)')   TSTYPE

!
!  MAKE SURE TEMP IS CONSTANT IN BAROTROPIC CASE
!
   IF(BAROTROPIC .AND. TSTYPE /= 'constant')THEN
     WRITE(IPT,*)'ERROR USING BAROTROPIC MODE'
     WRITE(IPT,*)'TEMPERATURE AND SALINITY MUST BE CONSTANT'
     CALL PSTOP
   END IF

!
!  CONSTANT SALINITY AND TEMPERATURE
!

   IF(TSTYPE == 'constant')THEN
     READ(INITS,*) T0,S0
     TS1 = T0
     SS1 = S0

!
!  LINEAR PROFILE OF SALINITY AND TEMPERATURE
!
   ELSE IF(TSTYPE == 'linear')THEN  !Linear in depth direction
     READ(INITS,*) IW
!    t11:  bottom temperature;
!    t22:  surface temperature
     DMIN=-1
     DO II=1,IW
       READ(INITS,*) T11, T22, S11, S22, DMAX
       DO K= 1, KSL
         IF(-DPTHSL(K) > DMIN.AND.-DPTHSL(K) <= DMAX) THEN
           T0=T22+(T22-T11)*DPTHSL(K)/DMAX
           S0=S22+(S22-S11)*DPTHSL(K)/DMAX
           DO I = 1,NNode
             TS1(I,K) = T0
             SS1(I,K) = S0
           END DO
         ENDIF
       END DO
     DMIN=DMAX
   ENDDO
!
!   OBSERVED SALINITY AND TEMPERATURE
!
   ELSE IF(TSTYPE == 'observed')THEN  !observed in KSL levels for all nodes
     ALLOCATE(TEMP1(NNodeGL,KSL),TEMP2(NNodeGL,KSL))
     DO I=1,NNodeGL 
       READ(INITS,*) (TEMP1(I,K), K=1,KSL)
       READ(INITS,*) (TEMP2(I,K), K=1,KSL)
     END DO

     IF(SERIAL)THEN
       TS1 = TEMP1
       SS1 = TEMP2
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,NNode
         TS1(I,1:KSL) = TEMP1(NGID(I),1:KSL)
         SS1(I,1:KSL) = TEMP2(NGID(I),1:KSL)
       END DO
     END IF
#    endif

     DEALLOCATE(TEMP1,TEMP2)
!
!   OTHER: REPORT FAILURE
!
   ELSE
     WRITE(IPT,*) 'SECOND LINE OF TEMP/SALINITY INIT FILE HAS AN ERROR'
     WRITE(IPT,*) 'OPTION SHOULD BE (constant/linear/observed)'
     WRITE(IPT,*) 'SEE initial_ts.f90'
     CALL PSTOP
   END IF

!==============================================================================|
!    CALCULATE MEAN VALUES OF TEMPERATURE/SALINITY/DENSITY                     !
!==============================================================================|

   TA = 0.0_SP
   SA = 0.0_SP
   IF(TSTYPE == 'constant')THEN  !constant 
     TA = T0
     SA = S0
	 
   ELSE  !Linear or observed, which is domain wide and covers until DPTHSL(KSL) from the surface
 
     IF(SERIAL)THEN
       DO K = 1, KSL  !loop all depth
         SCOUNT = 0.0_SP
         DO I = 1, NNode  !loop all nodes
           IF(-H(I) <= DPTHSL(K)) THEN
             SCOUNT = SCOUNT + 1.0_SP
             TA(K) = TA(K) + TS1(I,K)
             SA(K) = SA(K) + SS1(I,K)
           END IF
         END DO
         IF(SCOUNT >= 1.0_SP)THEN
           TA(K)=TA(K)/SCOUNT
           SA(K)=SA(K)/SCOUNT
         END IF
       END DO
     END IF


# if defined(MULTIPROCESSOR)
     IF(PAR)THEN
     DO K = 1, KSL
       SCOUNT = 0.0_SP
       DO I = 1, NNode
         IF(-H(I) <= DPTHSL(K)) THEN
           IF(NDE_ID(I) == 0)THEN !!INTERNAL NODE
             SCOUNT = SCOUNT + 1.0_SP
             TA(K) = TA(K) + TS1(I,K)
             SA(K) = SA(K) + SS1(I,K)
           ELSE  !!BOUNDARY NODE, ACCUMULATE FRACTION ONLY
             DO IB = 1,NBN
               IF(BN_LOC(IB) == I)THEN
                 FAC = 1./FLOAT(BN_MLT(IB))
                 SCOUNT = SCOUNT + FAC 
                 TA(K) = TA(K) + FAC*TS1(I,K)
                 SA(K) = SA(K) + FAC*SS1(I,K)
               END IF
             END DO
           END IF
         END IF
       END DO
       CALL MPI_REDUCE(TA(K),TATOT,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
       CALL MPI_REDUCE(SA(K),SATOT,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
       CALL MPI_REDUCE(SCOUNT,SCTOT,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
       SCOUNT = SCTOT
       TA(K)  = TATOT
       SA(K)  = SATOT
       CALL MPI_BCAST(SA(K),1, MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(TA(K),1, MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(SCOUNT,1, MPI_F,0,MPI_COMM_WORLD,IERR)
       IF(SCOUNT >= 1.0_SP)THEN
         TA(K)=TA(K)/SCOUNT
         SA(K)=SA(K)/SCOUNT
       END IF
     END DO !!LOOP OVER KSL
     END IF  !!PARALLEL
# endif
   END IF

!
!Vertical interpolation of TA, SA from data depth DPTHSL to model depth ZM
!with interpolated value in TI and SI for each node
!

   DO I=1,NNode
     DO K=1,KBM1
       ZM(K)=ZZ(I,K)*D(I)+EL(I)
     END DO

     CALL SINTER(DPTHSL,TA,ZM,TI,KSL,KBM1)
     CALL SINTER(DPTHSL,SA,ZM,SI,KSL,KBM1)

     DO K =1,KBM1
       T1(I,K) = TI(K)
       S1(I,K) = SI(K)
     END DO

   END DO

   IF(CTRL_DEN == 'pdensity'   ) CALL DENS
   IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
   IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3

!
!horizontal layer mean of T, S, RHO for nodes. Note TMEAN1, SMEAN1, RMEAN1 have no horizontal gradients in original (x,y,z) coordinates
!but they do have horizontal gradients in (x,y,sigma) coordinates
!

   DO K=1,KBM1
     DO I=1,NNode
       TMEAN1(I,K)=T1(I,K)
       SMEAN1(I,K)=S1(I,K)
       RMEAN1(I,K)=RHO1(I,K)
     END DO
   END DO

!
!horizontal layer mean of T, S, RHO for elements
!   
   DO I=1,MElem
     DO K=1,KBM1
       J1=NV(I,1)
       J2=NV(I,2)
       J3=NV(I,3)
       TMEAN(I,K)=(TMEAN1(J1,K)+TMEAN1(J2,K)+TMEAN1(J3,K))/3.0_SP
       SMEAN(I,K)=(SMEAN1(J1,K)+SMEAN1(J2,K)+SMEAN1(J3,K))/3.0_SP
       RMEAN(I,K)=(RMEAN1(J1,K)+RMEAN1(J2,K)+RMEAN1(J3,K))/3.0_SP
     END DO
     TMEAN(I,KB)=0.0_SP
     SMEAN(I,KB)=0.0_SP
     RMEAN(I,KB)=0.0_SP
   END DO

!==============================================================================|
!    CALCULATE TRUE VALUES OF TEMPERATURE/SALINITY/DENSITY                     !
!==============================================================================|

   DO I=1,NNode
     DO K=1,KSL
       TA(K)=TS1(I,K)
       SA(K)=SS1(I,K)
     END DO

     DO K= 1,KBM1
       ZM(K) =ZZ(I,K)*D(I)+EL(I)
     END DO

     CALL SINTER(DPTHSL,TA,ZM,TI,KSL,KBM1)
     CALL SINTER(DPTHSL,SA,ZM,SI,KSL,KBM1)

     DO K =1,KBM1
       T1(I,K) = TI(K)
       S1(I,K) = SI(K)
     END DO
   END DO


   IF(CTRL_DEN == 'sigma-t') THEN
     CALL DENS2
   ELSE IF(CTRL_DEN == 'pdensity') THEN
     CALL DENS
   ELSE IF(CTRL_DEN == 'sigma-t_stp') THEN
     CALL DENS3
   END IF


   DO I=1,MElem
     DO K=1,KBM1
       J1=NV(I,1)
       J2=NV(I,2)
       J3=NV(I,3)
       T(I,K)=(T1(J1,K)+T1(J2,K)+T1(J3,K))/3.0_SP
       S(I,K)=(S1(J1,K)+S1(J2,K)+S1(J3,K))/3.0_SP
       RHO(I,K)=(RHO1(J1,K)+RHO1(J2,K)+RHO1(J3,K))/3.0_SP
     END DO
     T(I,KB)=0.0_SP
     S(I,KB)=0.0_SP
     RHO(I,KB)=0.0_SP
   END DO

   RETURN
   END SUBROUTINE INITIAL_TS
!==============================================================================|
