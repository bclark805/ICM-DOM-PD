!==============================================================================|
!   READ IN SURFACE ELEVATION AND FLOW VELOCITIES			       |
!   SET DEPTH                                                                  |
!==============================================================================|

   SUBROUTINE INITIAL_UVEL        

!==============================================================================|
   USE ALL_VARS
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP)  :: XTEMP,YTEMP
   REAL(SP), PARAMETER :: FAC = 1.0_SP/3.0_SP
   REAL(SP), ALLOCATABLE :: TEMP(:),TEMP1(:,:),TEMP2(:,:) 
!==============================================================================|

!==============================================================================|
!          INITIALIZE FREE SURFACE ELEVATION AND RELATED QUANTITIES            !
!==============================================================================|


!----------READ IN THE ELEVATION-----------------------------------------------!
   
   ALLOCATE(TEMP(0:NNodeGL)) ; TEMP = 0.0_SP
   DO I=1,NNodeGL
!     READ(INELF,'(E14.5)') TEMP(I) !!EL(I)
     READ(INELF,*) TEMP(I) !!EL(I)
   END DO
   CLOSE(INELF)

!----------TRANSFORM TO LOCAL ARRAY IF PARALLEL--------------------------------!

   IF(SERIAL) EL = TEMP

# if defined (MULTIPROCESSOR)
   IF(PAR) THEN
   DO I=1,NNode
     EL(I) = TEMP(NGID(I))
   END DO
   DO I=1,NHN
     EL(I+NNode) = TEMP(HN_LST(I))
   END DO
   END IF
# endif
   DEALLOCATE(TEMP)
 
!-------- TRANSFER ELEVATION FIELD TO DEPTH AND OLD TIME LEVELS----------------!

   ET = EL 
   D  = EL + H
   DT = D

!-------- TRANSFER ELEVATIONS FROM NODES TO ELEMENTS---------------------------!

   DO I=1,MElem                                            !! GWC SHOULD HAVE TRANSFER SUB
     EL1(I) = FAC*SUM(EL(NV(I,1:3)))
      D1(I) = FAC*SUM( D(NV(I,1:3)))
   END DO

   ET1 = EL1
   DT1 = D1  

!==============================================================================|
!          INITIALIZE INTERNAL AND EXTERNAL VELOCITIY FIELDS                   !
!==============================================================================|

!---------READ IN INTERNAL VELOCITY FIELD--------------------------------------!

   ALLOCATE(TEMP1(0:MElemGL,KB),TEMP2(0:MElemGL,KB))
   DO I=1,MElemGL
     READ(INUVF,'(200e14.5)') (TEMP1(I,K),TEMP2(I,K),K=1,KBM1)   !!U,V
   END DO
   CLOSE(INUVF)

!----------TRANSFORM TO LOCAL ARRAY IF PARALLEL--------------------------------!

   IF(SERIAL)THEN
     U = TEMP1
     V = TEMP2
   END IF

# if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,MElem
       U(I,1:KBM1) = TEMP1(EGID(I),1:KBM1) 
       V(I,1:KBM1) = TEMP2(EGID(I),1:KBM1) 
     END DO
     DO I=1,MHE
       U(I+MElem,1:KBM1) = TEMP1(HE_LST(I),1:KBM1) 
       V(I+MElem,1:KBM1) = TEMP2(HE_LST(I),1:KBM1) 
     END DO
   END IF
#  endif
   DEALLOCATE(TEMP1,TEMP2)


!-------- COMPUTE VERTICALLY AVERAGED VELOCITIES-------------------------------!

   UA = 0.0_SP ; VA = 0.0_SP

   DO I=1,MTElem
!     UA(I) = SUM(U(I,1:KBM1))/FLOAT(KBM1)
!     VA(I) = SUM(V(I,1:KBM1))/FLOAT(KBM1)
     UA(I) = SUM(U(I,1:KBM1)*DZ1(I,1:KBM1))
     VA(I) = SUM(V(I,1:KBM1)*DZ1(I,1:KBM1))
   END DO

   RETURN
   END SUBROUTINE INITIAL_UVEL
!==============================================================================|
