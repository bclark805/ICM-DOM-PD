!==============================================================================|
!   Implement Fresh Water Boundary Condition on Bottom (Groundwater)           |
!==============================================================================|

   SUBROUTINE BCOND_BFW(IDX)           

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IDX
   REAL(SP) :: FACT,UFACT
   INTEGER  :: L1,L2,IERR 
!------------------------------------------------------------------------------|

   SELECT CASE(IDX)
   
   CASE(1)
   
   IF(IBFW <= 0)RETURN

   CALL BRACKET(BFW_TM,THOUR,L1,L2,FACT,UFACT,IERR)

   IF(IERR == -1)THEN
     BFWDIS3 = 0.0_SP
     BFWTDIS3 = 0.0_SP
     BFWSDIS3 = 0.0_SP
   ELSE
     BFWDIS3(:)=UFACT*BFWQDIS(:,L1)+FACT*BFWQDIS(:,L2)
     BFWDIS3   =BFWDIS3*RAMP

     BFWTDIS3(:)=UFACT*BFWQTDIS(:,L1)+FACT*BFWQTDIS(:,L2)
     BFWTDIS3   =BFWTDIS3*RAMP

     BFWSDIS3(:)=UFACT*BFWQSDIS(:,L1)+FACT*BFWQSDIS(:,L2)
     BFWSDIS3   =BFWSDIS3*RAMP
   END IF

   CASE(2)
# if !defined (SEMI_IMPLICIT)   
   IF(IBFW <= 0)RETURN

   CALL BRACKET(BFW_TM,THOUR1,L1,L2,FACT,UFACT,IERR)

   IF(IERR == -1)THEN
     BFWDIS2 = 0.0_SP
     BFWTDIS2 = 0.0_SP
     BFWSDIS2 = 0.0_SP
   ELSE
     BFWDIS2(:)=UFACT*BFWQDIS(:,L1)+FACT*BFWQDIS(:,L2)
     BFWDIS2   =BFWDIS2*RAMP

!     BFWTDIS2(:)=UFACT*BFWQTDIS(:,L1)+FACT*BFWQTDIS(:,L2)
!     BFWTDIS2   =BFWTDIS2*RAMP

!     BFWSDIS2(:)=UFACT*BFWQSDIS(:,L1)+FACT*BFWQSDIS(:,L2)
!     BFWSDIS2   =BFWSDIS2*RAMP
   END IF
# endif
   END SELECT
   
   RETURN
   END SUBROUTINE BCOND_BFW
!==============================================================================|
