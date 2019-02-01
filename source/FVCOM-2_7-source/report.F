
!==============================================================================|
!     REPORT INITIAL INFORMATION                                               |
!==============================================================================|

   SUBROUTINE REPORT(INFO_STRING)

!==============================================================================|
   USE ALL_VARS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   CHARACTER(LEN=*) :: INFO_STRING     !!INFORMATION STRING
   INTEGER :: E3TOT,ESTOT,IERR
   REAL(DP), DIMENSION(17) :: SBUF,RBUF1,RBUF2,RBUF3
   
   REAL(SP), ALLOCATABLE :: Q21(:,:),Q2L1(:,:),L1(:,:)
   REAL(SP), ALLOCATABLE :: KH1(:,:),KQ1(:,:)
#  if defined (GOTM)
   REAL(SP), ALLOCATABLE :: TKE1(:,:),TEPS1(:,:)
#  endif   
   INTEGER :: I,J,K
!==============================================================================|

   ALLOCATE(Q21(1:MElem,KBM1));   Q21   = 0.0_SP
   ALLOCATE(Q2L1(1:MElem,KBM1));  Q2L1  = 0.0_SP
   ALLOCATE(L1(1:MElem,KBM1));    L1    = 0.0_SP
   ALLOCATE(KH1(1:MElem,KBM1));   KH1   = 0.0_SP
   ALLOCATE(KQ1(1:MElem,KBM1));   KQ1   = 0.0_SP
#  if defined (GOTM)
   ALLOCATE(TKE1(1:MElem,KBM1));  TKE1  = 0.0_SP
   ALLOCATE(TEPS1(1:MElem,KBM1)); TEPS1 = 0.0_SP
#  endif   
   
   DO K=1,KBM1
     DO I=1,MElem
       DO J=1,3
         Q21(I,K)  = Q21(I,K)+Q2(NV(I,J),K)
         Q2L1(I,K) = Q2L1(I,K)+Q2L(NV(I,J),K)
         L1(I,K)   = L1(I,K)+L(NV(I,J),K)
         KH1(I,K)  = KH1(I,K)+KH(NV(I,J),K)
         KQ1(I,K)  = KQ1(I,K)+KQ(NV(I,J),K)
#  if defined (GOTM)
         TKE1(I,K) = TKE1(I,K)+TKE(NV(I,J),K)
	 TEPS1(I,K)= TEPS1(I,K)+TEPS(NV(I,J),K)
#  endif   
       END DO
       Q21(I,K)  = Q21(I,K)/3.0_SP
       Q2L1(I,K) = Q2L1(I,K)/3.0_SP
       L1(I,K)   = L1(I,K)/3.0_SP
       KH1(I,K)  = KH1(I,K)/3.0_SP
       KQ1(I,K)  = KQ1(I,K)/3.0_SP
#  if defined (GOTM)
         TKE1(I,K) = TKE1(I,K)/3.0_SP
	 TEPS1(I,K)= TEPS1(I,K)/3.0_SP
#  endif   
     END DO
   END DO    	    
      
   SBUF(1)  = SUM(DBLE(UA(1:MElem)))
   SBUF(2)  = SUM(DBLE(VA(1:MElem)))
   SBUF(3)  = SUM(DBLE(EL1(1:MElem)))
   SBUF(4)  = SUM(DBLE(H1(1:MElem)))
   SBUF(5)  = SUM(DBLE(U(1:MElem,1:KBM1)))
   SBUF(6)  = SUM(DBLE(V(1:MElem,1:KBM1)))
   SBUF(7)  = SUM(DBLE(S(1:MElem,1:KBM1)))
   SBUF(8)  = SUM(DBLE(T(1:MElem,1:KBM1)))
#  if defined (GOTM)
   SBUF(9)  = SUM(DBLE(TKE1(1:MElem,2:KBM1)))
   SBUF(10) = SUM(DBLE(TEPS1(1:MElem,2:KBM1)))
#  else
   SBUF(9)  = SUM(DBLE(Q21(1:MElem,2:KBM1)))
   SBUF(10) = SUM(DBLE(Q2L1(1:MElem,2:KBM1)))
#  endif
   SBUF(11) = SUM(DBLE(L1(1:MElem,1:KBM1)))
   SBUF(12) = SUM(DBLE(KM1(1:MElem,1:KBM1)))
   SBUF(13) = SUM(DBLE(KQ1(1:MElem,1:KBM1)))
   SBUF(14) = SUM(DBLE(KH1(1:MElem,1:KBM1)))
   SBUF(15) = SUM(DBLE(RHO(1:MElem,1:KBM1)))
   SBUF(16) = SUM(DBLE(D1(1:MElem)))
#  if !defined (WET_DRY)
   SBUF(17) = FLOAT(MElem)
#  else
   SBUF(17) = SUM(ISWET_CELL_CURRENTSTEP(1:MElem))
#  endif


   RBUF1 = SBUF
# if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(SBUF,RBUF1,17,MPI_DP,MPI_SUM,0,MPI_COMM_WORLD,IERR)
# endif

   SBUF(1)  = MAXVAL(UA(1:MElem))
   SBUF(2)  = MAXVAL(VA(1:MElem))
   SBUF(3)  = MAXVAL(EL(1:NNode))
   SBUF(4)  = MAXVAL(H(1:NNode))
   SBUF(5)  = MAXVAL(U(1:MElem,1:KBM1))
   SBUF(6)  = MAXVAL(V(1:MElem,1:KBM1))
   SBUF(7)  = MAXVAL(S1(1:NNode,1:KBM1))
   SBUF(8)  = MAXVAL(T1(1:NNode,1:KBM1))
#  if defined (GOTM)
   SBUF(9)  = MAXVAL(TKE(1:NNode,2:KBM1))
   SBUF(10) = MAXVAL(TEPS(1:NNode,2:KBM1))
#  else
   SBUF(9)  = MAXVAL(Q2(1:NNode,1:KBM1))
   SBUF(10) = MAXVAL(Q2L(1:NNode,1:KBM1))
#  endif
   SBUF(11) = MAXVAL(L(1:NNode,1:KBM1))
   SBUF(12) = MAXVAL(KM(1:NNode,1:KBM1))
   SBUF(13) = MAXVAL(KQ(1:NNode,1:KBM1))
   SBUF(14) = MAXVAL(KH(1:NNode,1:KBM1))
   SBUF(15) = MAXVAL(RHO1(1:NNode,1:KBM1))
   SBUF(16) = MAXVAL(D(1:NNode))

   RBUF2 = SBUF
# if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(SBUF,RBUF2,16,MPI_DP,MPI_MAX,0,MPI_COMM_WORLD,IERR)
# endif

   SBUF(1)  = MINVAL(UA(1:MElem))
   SBUF(2)  = MINVAL(VA(1:MElem))
   SBUF(3)  = MINVAL(EL(1:NNode))
   SBUF(4)  = MINVAL(H(1:NNode))
   SBUF(5)  = MINVAL(U(1:MElem,1:KBM1))
   SBUF(6)  = MINVAL(V(1:MElem,1:KBM1))
   SBUF(7)  = MINVAL(S1(1:NNode,1:KBM1))
   SBUF(8)  = MINVAL(T1(1:NNode,1:KBM1))
#  if defined (GOTM)
   SBUF(9)  = MINVAL(TKE(1:NNode,2:KBM1))
   SBUF(10)  = MINVAL(TEPS(1:NNode,2:KBM1))
#  else
   SBUF(9)  = MINVAL(Q2(1:NNode,1:KBM1))
   SBUF(10)  = MINVAL(Q2L(1:NNode,1:KBM1))
#  endif
   SBUF(11) = MINVAL(L(1:NNode,1:KBM1))
   SBUF(12) = MINVAL(KM(1:NNode,1:KBM1))
   SBUF(13) = MINVAL(KQ(1:NNode,1:KBM1))
   SBUF(14) = MINVAL(KH(1:NNode,1:KBM1))
   SBUF(15) = MINVAL(RHO1(1:NNode,1:KBM1))
   SBUF(16) = MINVAL(D(1:NNode))

   RBUF3 = SBUF
# if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(SBUF,RBUF3,16,MPI_DP,MPI_MIN,0,MPI_COMM_WORLD,IERR)
# endif

   IF(MSR)THEN
   IF(LEN_TRIM(INFO_STRING) /= 0)THEN
     WRITE(IPT,*  )'===================',TRIM(INFO_STRING),'======================'
   END IF
   RBUF1(15) = (RBUF1(15)+MElemGL*KBM1)*1000.
   RBUF2(15) = (RBUF2(15)+1.)*1000.
   RBUF3(15) = (RBUF3(15)+1.)*1000.
   E3TOT = DBLE(MElemGL*KBM1)
   ESTOT = DBLE(MElemGL)
   WRITE(IPT,*  )'!  QUANTITY              :     AVG           MAX         MIN'
   WRITE(IPT,100)'!  EXTERNAL UVEL         :',RBUF1(1)/ESTOT,RBUF2(1),RBUF3(1)
   WRITE(IPT,100)'!  EXTERNAL VVEL         :',RBUF1(2)/ESTOT,RBUF2(2),RBUF3(2)
   WRITE(IPT,100)'!  FREE SURFACE          :',RBUF1(3)/ESTOT,RBUF2(3),RBUF3(3)
   WRITE(IPT,100)'!  BATH DEPTH            :',RBUF1(4)/ESTOT,RBUF2(4),RBUF3(4)
   WRITE(IPT,100)'!  INTERNAL UVEL         :',RBUF1(5)/E3TOT,RBUF2(5),RBUF3(5)
   WRITE(IPT,100)'!  INTERNAL VVEL         :',RBUF1(6)/E3TOT,RBUF2(6),RBUF3(6)
   WRITE(IPT,100)'!  SALINITY              :',RBUF1(7)/E3TOT,RBUF2(7),RBUF3(7)
   WRITE(IPT,100)'!  TEMPERATURE           :',RBUF1(8)/E3TOT,RBUF2(8),RBUF3(8)
#  if defined (GOTM)
   WRITE(IPT,100)'!  TURBULENT KE          :',RBUF1(9)/E3TOT,RBUF2(9),RBUF3(9)
   WRITE(IPT,100)'!  TURBULENT DISSIPATION :',RBUF1(10)/E3TOT,RBUF2(10),RBUF3(10)
#  else
   WRITE(IPT,100)'!  TURBULENT KE          :',RBUF1(9)/E3TOT,RBUF2(9),RBUF3(9)
   WRITE(IPT,100)'!  TURB KE*L             :',RBUF1(10)/E3TOT,RBUF2(10),RBUF3(10)
#  endif
   WRITE(IPT,100)'!  TURB LENGTH SCALE     :',RBUF1(11)/E3TOT,RBUF2(11),RBUF3(11)
   WRITE(IPT,100)'!  KM                    :',RBUF1(12)/E3TOT,RBUF2(12),RBUF3(12)
   WRITE(IPT,100)'!  KQ                    :',RBUF1(13)/E3TOT,RBUF2(13),RBUF3(13)
   WRITE(IPT,100)'!  KH                    :',RBUF1(14)/E3TOT,RBUF2(14),RBUF3(14)
   WRITE(IPT,100)'!  DENSITY               :',RBUF1(15)/E3TOT,RBUF2(15),RBUF3(15)
   WRITE(IPT,100)'!  DEPTH                 :',RBUF1(16)/ESTOT,RBUF2(16),RBUF3(16)
#  if defined (WET_DRY)
   WRITE(IPT,*  )'!  WET/DRY INFO          :   #WET       #DRY             %WET'
   IF(RBUF1(17) == FLOAT(MElemGL))THEN
   WRITE(IPT,*)'!  NO DRY POINTS          '
   ELSE
   WRITE(IPT,101)'!  WET/DRY DATA          :',INT(RBUF1(17)),MElemGL-INT(RBUF1(17)),100.*RBUF1(17)/FLOAT(MElemGL)
   END IF
#  endif

   END IF

   DEALLOCATE(Q21,Q2L1,L1)
   DEALLOCATE(KH1,KQ1)
#  if defined (GOTM)
   DEALLOCATE(TKE1,TEPS1)
#  endif
   
   RETURN
   100 FORMAT(1X,A26,3F12.6)
   101 FORMAT(1X,A26,2I12,F12.6)
   END SUBROUTINE REPORT 
!==============================================================================|

