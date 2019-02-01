MODULE MOD_SEMI_IMPLICIT
# if defined (SEMI_IMPLICIT) && (MULTIPROCESSOR)

  USE MOD_PREC, ONLY : SP, DP
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: SET_IMPLICIT_PARAM
  PUBLIC :: ALLOC_VARS_SEMI
  PUBLIC :: SEMI_IMPLICIT_EL
  PUBLIC :: UPDATE_EL
  PUBLIC :: UV2D_SBD
  PUBLIC :: UV3D_SBD

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  REAL(SP), PUBLIC  ::  IFCETA, BEDF
  REAL(SP), PUBLIC, ALLOCATABLE  :: XFLUX3(:,:), YFLUX3(:,:)
  REAL(SP), ALLOCATABLE :: UA_N(:), VA_N(:) 
  REAL(SP), ALLOCATABLE :: U_N(:,:), V_N(:,:)

  PetscViewer :: viewer

  CONTAINS

  SUBROUTINE SET_IMPLICIT_PARAM

  USE CONTROL
  USE MOD_INP
  IMPLICIT NONE
  INTEGER  :: ISCAN
  CHARACTER(LEN=120) :: FNAME

  FNAME = "./"//trim(casename)//"_run.dat"

!----------------------------------------------------------------------------|
!   implicit factor ceta 
!----------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"IFCETA",FSCAL = IFCETA)
  IF(ISCAN /= 0)THEN
    WRITE(IPT,*)'ERROR READING IFCETA: ',ISCAN
    CALL PSTOP
  END IF

!----------------------------------------------------------------------------|
!   implicit factor ceta
!----------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"BEDF",FSCAL = BEDF)
  IF(ISCAN /= 0)THEN
    WRITE(IPT,*)'ERROR READING BEDF: ',ISCAN
    CALL PSTOP
  END IF

!==============================================================================|
!            SCREEN REPORT                                                     !
!==============================================================================|
  IF(MSR) THEN
    WRITE(IPT,*) '!                                                    !'
    WRITE(IPT,*) '!-------------- SPECIFY IMPLICIT FACTOR -------------!'
    WRITE(IPT,*) '!                                                    !'
    WRITE(IPT,*) '!  # IMPLICIT FACTOR SELECTED                        :',IFCETA
    WRITE(IPT,*) '!  # SOLID BOUNDARY ENERGY DAMPING FACTOR SELECTED   :',BEDF
    WRITE(IPT,*) '!                                                    !'
  ENDIF

  RETURN
  END SUBROUTINE SET_IMPLICIT_PARAM

  SUBROUTINE ALLOC_VARS_SEMI

  USE MOD_PREC, ONLY : SP, DP
  USE ALL_VARS, ONLY : MTElem,KBM1,KB
# if defined (SPHERICAL) && (NORTHPOLE)  
  USE MOD_NORTHPOLE, ONLY : XFLUX3_NP, YFLUX3_NP, UBETA_NP, VBETA_NP
# endif
  IMPLICIT NONE

  ALLOCATE(XFLUX3(0:MTElem,KBM1))         ; XFLUX3    = 0.0_SP
  ALLOCATE(YFLUX3(0:MTElem,KBM1))         ; YFLUX3    = 0.0_SP
# if defined (SPHERICAL) && (NORTHPOLE)
  ALLOCATE(XFLUX3_NP(0:MTElem,KBM1))      ; XFLUX3_NP = 0.0_SP
  ALLOCATE(YFLUX3_NP(0:MTElem,KBM1))      ; YFLUX3_NP = 0.0_SP
  ALLOCATE(UBETA_NP(0:MTElem,KB))         ; UBETA_NP  = 0.0_SP   
  ALLOCATE(VBETA_NP(0:MTElem,KB))         ; VBETA_NP  = 0.0_SP 
# endif
  ALLOCATE(UA_N(0:MTElem))                ; UA_N      = 0.0_SP
  ALLOCATE(VA_N(0:MTElem))                ; VA_N      = 0.0_SP
  ALLOCATE(U_N(0:MTElem,KBM1))            ; U_N       = 0.0_SP
  ALLOCATE(V_N(0:MTElem,KBM1))            ; V_N       = 0.0_SP

  RETURN
  END SUBROUTINE ALLOC_VARS_SEMI

  SUBROUTINE SEMI_IMPLICIT_EL

  USE LIMS
  USE CONTROL
  USE MOD_OBCS
  USE ALL_VARS
  USE BCS
  USE MOD_PETSC, ONLY: A_EL,X_EL,B_EL,BL_EL,XL_EL,N_VERTS,L2G_EL,G2L_EL,XVALS_EL,ALO_2_PLO_NODE,PLO_2_ALO_NODE,PETSc_SETICS_EL,PETSc_SOLVER_EL
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
# if defined (SPHERICAL)   
  USE MOD_SPHERICAL
# if defined(NORTHPOLE)
  USE MOD_NORTHPOLE,   ONLY :  NPE, NPEDGE_LST, NPCV,  NCEDGE_LST, NODE_NORTHPOLE, NP, NP_LST, CELL_NORTHAREA, XFLUX3_NP, YFLUX3_NP
# endif
# endif   
# if defined (WET_DRY)
  USE MOD_WD, ONLY : WET_DRY_ON, ISWET_CELL_CURRENTSTEP, ISWET_CELL_LAST_INT_STEP, WET_JUDGE
# endif

  IMPLICIT NONE 
  INTEGER II, I, J, JJ, J1, J2, IK, K, I1, IA, IB, TMP2, TMP3
  INTEGER  NNZ2, LN, LOC
  INTEGER  PETSc_POS, NODE, PROW1, PROW2, PCOL1, PCOL2, PCOL3
  REAL(SP) DIJ, UIJ, VIJ, EXFLUX, COEF
  REAL(SP) Y3Y2,Y1Y3,Y2Y1,X3X2,X1X3,X2X1
  REAL(SP) :: XFLUX(0:NTNode), XFLUX_TMP(0:MTElem), YFLUX_TMP(0:MTElem), XFLUX_TMP_NP(0:MTElem), YFLUX_TMP_NP(0:MTElem)

  REAL(SP) :: UA_TMP, VA_TMP, WUSURF_TMP, WVSURF_TMP, WUBOT_TMP, WVBOT_TMP
  REAL(SP) :: U_TMP, V_TMP, UF_TMP, VF_TMP
  REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP,FLUX_TMP, ELIJ
  REAL(SP) :: DLTXE_TMP,DLTYE_TMP,XTMP, XTMP1, DLTXC_TMP, DLTYC_TMP
  REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP,VX3_TMP,VY3_TMP
  REAL(SP) :: XXX1, XXX2, XXX3, YYY1, YYY2, YYY3

  PetscReal :: STERM, VCOL1,VCOL2,VCOL3
  PetscInt  :: IERR
  PetscScalar :: ZERO1   =  0.0D0  

# if defined (TWO_D_MODEL)
  XFLUX_TMP = 0.0_SP
  YFLUX_TMP = 0.0_SP
  DO I=1, MElem
#   if defined (WET_DRY)
    IF(ISWET_CELL_LAST_INT_STEP(I) == 1) THEN
#   endif
    XFLUX_TMP(I) = UA(I)*DT1(I) - DTI*IFCETA*ADVUA(I)/ART(I) + DTI*IFCETA*( -WUSURF(I) + WUBOT(I) ) - IFCETA*BEDF*UA_N(I)*DT1(I)
    YFLUX_TMP(I) = VA(I)*DT1(I) - DTI*IFCETA*ADVVA(I)/ART(I) + DTI*IFCETA*( -WVSURF(I) + WVBOT(I) ) - IFCETA*BEDF*VA_N(I)*DT1(I)
!!!    XFLUX_TMP(I) = UA(I)*H1(I) - DTI*IFCETA*ADVUA(I)/ART(I)
!!!    YFLUX_TMP(I) = VA(I)*H1(I) - DTI*IFCETA*ADVVA(I)/ART(I)
#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,XFLUX_TMP,YFLUX_TMP)
# endif

# else
  XFLUX_TMP = 0.0_SP
  YFLUX_TMP = 0.0_SP
  DO I=1, MElem
#   if defined (WET_DRY)
    IF(ISWET_CELL_LAST_INT_STEP(I) == 1) THEN
#   endif
    DO K=1, KBM1
      XFLUX_TMP(I) = XFLUX_TMP(I) - DTI*IFCETA*XFLUX3(I,K)/ART(I) - IFCETA*BEDF*U_N(I,K)*DT1(I)*DZ1(I,K)
      YFLUX_TMP(I) = YFLUX_TMP(I) - DTI*IFCETA*YFLUX3(I,K)/ART(I) - IFCETA*BEDF*V_N(I,K)*DT1(I)*DZ1(I,K)
    ENDDO
    XFLUX_TMP(I) = XFLUX_TMP(I) + UA(I)*DT1(I) + DTI*IFCETA*( -WUSURF(I) + WUBOT(I) )
    YFLUX_TMP(I) = YFLUX_TMP(I) + VA(I)*DT1(I) + DTI*IFCETA*( -WVSURF(I) + WVBOT(I) )
!!!    XFLUX_TMP(I) = XFLUX_TMP(I) + UA(I)*H1(I)
!!!    YFLUX_TMP(I) = YFLUX_TMP(I) + VA(I)*H1(I)
#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,XFLUX_TMP,YFLUX_TMP)
# endif

# if defined (SPHERICAL) && (NORTHPOLE)
  XFLUX_TMP_NP = 0.0_SP
  YFLUX_TMP_NP = 0.0_SP
  DO II=1, NP
    I = NP_LST(II)
    IF(CELL_NORTHAREA(I) ==1 ) THEN
      DO K=1, KBM1
        XFLUX_TMP_NP(I) = XFLUX_TMP_NP(I) - DTI*IFCETA*XFLUX3_NP(I,K)/ART(I)
        YFLUX_TMP_NP(I) = YFLUX_TMP_NP(I) - DTI*IFCETA*YFLUX3_NP(I,K)/ART(I)
      ENDDO
      UA_TMP = -VA(I)*COS(XC(I)*PI/180.0_SP)-UA(I)*SIN(XC(I)*PI/180.0_SP)
      VA_TMP = -VA(I)*SIN(XC(I)*PI/180.0_SP)+UA(I)*COS(XC(I)*PI/180.0_SP)
      WUSURF_TMP = -WVSURF(I)*COS(XC(I)*PI/180.0_SP)-WUSURF(I)*SIN(XC(I)*PI/180.0_SP)
      WVSURF_TMP = -WVSURF(I)*SIN(XC(I)*PI/180.0_SP)+WUSURF(I)*COS(XC(I)*PI/180.0_SP)
      WUBOT_TMP  = -WVBOT(I)*COS(XC(I)*PI/180.0_SP)-WUBOT(I)*SIN(XC(I)*PI/180.0_SP)
      WVBOT_TMP  = -WVBOT(I)*SIN(XC(I)*PI/180.0_SP)+WUBOT(I)*COS(XC(I)*PI/180.0_SP)
      XFLUX_TMP_NP(I) = XFLUX_TMP_NP(I) + UA_TMP*DT1(I) + DTI*IFCETA*( -WUSURF_TMP + WUBOT_TMP )
      YFLUX_TMP_NP(I) = YFLUX_TMP_NP(I) + VA_TMP*DT1(I) + DTI*IFCETA*( -WVSURF_TMP + WVBOT_TMP )
!!!      XFLUX_TMP_NP(I) = XFLUX_TMP_NP(I) + UA_TMP*H1(I)
!!!      YFLUX_TMP_NP(I) = YFLUX_TMP_NP(I) + VA_TMP*H1(I)
    ENDIF
  ENDDO
# endif
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,XFLUX_TMP_NP,YFLUX_TMP_NP)
# endif

# endif

  CALL VecSet(BL_EL,ZERO1,IERR);CHKERRQ(IERR)
  CALL VecSet(B_EL,ZERO1,IERR);CHKERRQ(IERR)

  XFLUX   = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)

    UIJ=XFLUX_TMP(I1)
    VIJ=YFLUX_TMP(I1)
    EXFLUX= -UIJ*DLTYE(I)+VIJ*DLTXE(I)
    XFLUX(IA)=XFLUX(IA)-EXFLUX
    XFLUX(IB)=XFLUX(IB)+EXFLUX
  ENDDO

# if defined (SPHERICAL) && (NORTHPOLE)
  IF(NODE_NORTHPOLE/=0) XFLUX(NODE_NORTHPOLE) = 0.0_SP   

  DO II=1,NPCV
    I = NCEDGE_LST(II)
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)

    UIJ = XFLUX_TMP_NP(I1)
    VIJ = YFLUX_TMP_NP(I1)

    IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
      VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)!&
!                * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))
      VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)!&
!                * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))

      VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)!&
!                * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))
      VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)!&
!                * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))

      DLTXE_TMP = VX2_TMP-VX1_TMP
      DLTYE_TMP = VY2_TMP-VY1_TMP

      EXFLUX_TMP = -UIJ*DLTYE_TMP+VIJ*DLTXE_TMP
    END IF

    IF(IA == NODE_NORTHPOLE) THEN
      XFLUX(IA) = XFLUX(IA)-EXFLUX_TMP
    ELSE IF(IB == NODE_NORTHPOLE)THEN
      XFLUX(IB) = XFLUX(IB)+EXFLUX_TMP
    END IF

  END DO 
# endif

! this part may need be changed to semi-implicit form
  XFLUX = XFLUX+(QEVAP3-QPREC3)*ROFVROS*ART1

  IF(IBFW > 0)THEN
    DO I=1,NNode
      DO J=1,IBFW
        IF(I == NODE_BFW(J))THEN
          XFLUX(I)=XFLUX(I)+BFWDIS3(J)*ROFVROS*ART1(I)
        END IF
      END DO
    END DO
  END IF

  IF(NUMQBC >= 1) THEN
    IF(INFLOW_TYPE == 'node') THEN
      DO J=1,NUMQBC
        JJ=INODEQ(J)
        XFLUX(JJ)=XFLUX(JJ)-QDIS(J)
      END DO
    ELSE IF(INFLOW_TYPE == 'edge') THEN
      DO J=1,NUMQBC
        J1=N_ICELLQ(J,1)
        J2=N_ICELLQ(J,2)
        XFLUX(J1)=XFLUX(J1)-QDIS(J)*RDISQ(J,1)
        XFLUX(J2)=XFLUX(J2)-QDIS(J)*RDISQ(J,2)
      END DO
    END IF
  END IF

  DO I=1, NNode
    PETSc_POS = ALO_2_PLO_NODE(I)
    IF (PETSc_POS > N_VERTS) CYCLE
    IF( ISONB(I)/=2 ) THEN
      STERM = ET(I) - DTI/ART1(I)*XFLUX(I)
    ELSE
      STERM = ELF(I)
    ENDIF
    CALL VecSetValues(BL_EL,1,PETSc_POS-1,STERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

  CALL VecScatterBegin(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)

  CALL MatZeroEntries(A_EL,IERR);CHKERRQ(IERR)

  COEF  = DTI**2*IFCETA**2
  DO I=1, NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)

    IF((ISONB(IA)+ISONB(IB))<4) THEN

#     if defined (WET_DRY)
      IF(ISWET_CELL_LAST_INT_STEP(I1) == 1) THEN
#     endif

      DO J=1, 3
        IF(NV(I1,J)==IA) THEN
          TMP2 = NV(I1,J+1-INT((J+1)/4)*3)
          TMP3 = NV(I1,J+2-INT((J+2)/4)*3)
        ENDIF
      ENDDO

      PROW1 = ALO_2_PLO_NODE(IA)
      PROW2 = ALO_2_PLO_NODE(IB)

      PCOL1 = ALO_2_PLO_NODE(IA)
      PCOL2 = ALO_2_PLO_NODE(IB)
      PCOL3 = ALO_2_PLO_NODE(TMP3)

#     if defined (SPHERICAL)

#     if defined (NORTHPOLE)
      IF (CELL_NORTHAREA(I1)==1) THEN

        VX1_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * COS(VX(IA)*PI/180.0_SP)&
                * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))
        VY1_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * SIN(VX(IA)*PI/180.0_SP)&
                * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))
    
        VX2_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * COS(VX(IB)*PI/180.0_SP)&
                * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))
        VY2_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * SIN(VX(IB)*PI/180.0_SP)&
                * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))
    
        VX3_TMP = REARTH * COS(VY(TMP3)*PI/180.0_SP) * COS(VX(TMP3)*PI/180.0_SP)&
                * 2._SP /(1._SP+sin(VY(TMP3)*PI/180.0_SP))
        VY3_TMP = REARTH * COS(VY(TMP3)*PI/180.0_SP) * SIN(VX(TMP3)*PI/180.0_SP)&
                * 2._SP /(1._SP+sin(VY(TMP3)*PI/180.0_SP))

        X3X2  = VX3_TMP-VX2_TMP
        X1X3  = VX1_TMP-VX3_TMP
        X2X1  = VX2_TMP-VX1_TMP
        Y3Y2  = VY3_TMP-VY2_TMP
        Y1Y3  = VY1_TMP-VY3_TMP
        Y2Y1  = VY2_TMP-VY1_TMP
      
        XXX1  =( X3X2*COS(XC(I1)*PI/180.0_SP)+Y3Y2*SIN(XC(I1)*PI/180.0_SP) )/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
        XXX2  =( X1X3*COS(XC(I1)*PI/180.0_SP)+Y1Y3*SIN(XC(I1)*PI/180.0_SP) )/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
        XXX3  =( X2X1*COS(XC(I1)*PI/180.0_SP)+Y2Y1*SIN(XC(I1)*PI/180.0_SP) )/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
        YYY1  =(-X3X2*SIN(XC(I1)*PI/180.0_SP)+Y3Y2*COS(XC(I1)*PI/180.0_SP) )/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
        YYY2  =(-X1X3*SIN(XC(I1)*PI/180.0_SP)+Y1Y3*COS(XC(I1)*PI/180.0_SP) )/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
        YYY3  =(-X2X1*SIN(XC(I1)*PI/180.0_SP)+Y2Y1*COS(XC(I1)*PI/180.0_SP) )/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))

        VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(XXX1*DLTYE(I)-YYY1*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IA)
        VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(XXX2*DLTYE(I)-YYY2*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IA)
        VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(XXX3*DLTYE(I)-YYY3*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IA)
!!!        VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(XXX1*DLTYE(I)-YYY1*DLTXE(I))*H1(I1)/ART(I1)/ART1(IA)
!!!        VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(XXX2*DLTYE(I)-YYY2*DLTXE(I))*H1(I1)/ART(I1)/ART1(IA)
!!!        VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(XXX3*DLTYE(I)-YYY3*DLTXE(I))*H1(I1)/ART(I1)/ART1(IA)

      ELSE
#     endif
      XTMP  = VX(TMP3)-VX(TMP2)
      IF(XTMP >  180.0_SP)THEN
          X3X2 = -360.0_SP*TPI+XTMP*TPI
        ELSE IF(XTMP < -180.0_SP)THEN
          X3X2 =  360.0_SP*TPI+XTMP*TPI
        ELSE
          X3X2 = XTMP*TPI
      END IF
      XTMP  = VX(IA)  -VX(TMP3)
      IF(XTMP >  180.0_SP)THEN
          X1X3 = -360.0_SP*TPI+XTMP*TPI
        ELSE IF(XTMP < -180.0_SP)THEN
          X1X3 =  360.0_SP*TPI+XTMP*TPI
        ELSE
          X1X3 = XTMP*TPI
      END IF
      XTMP  = VX(TMP2)-VX(IA)
      IF(XTMP >  180.0_SP)THEN
          X2X1 = -360.0_SP*TPI+XTMP*TPI
        ELSE IF(XTMP < -180.0_SP)THEN
          X2X1 =  360.0_SP*TPI+XTMP*TPI
        ELSE
          X2X1 = XTMP*TPI
      END IF
      Y3Y2  = (VY(TMP3)-VY(TMP2))*TPI
      Y1Y3  = (VY(IA)  -VY(TMP3))*TPI
      Y2Y1  = (VY(TMP2)-VY(IA))*TPI

        VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y3Y2*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
        VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y1Y3*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
        VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y2Y1*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
!!!        VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y3Y2*DLTYE(I))*H1(I1)/ART(I1)/ART1(IA)
!!!        VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y1Y3*DLTYE(I))*H1(I1)/ART(I1)/ART1(IA)
!!!        VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y2Y1*DLTYE(I))*H1(I1)/ART(I1)/ART1(IA)
#     if defined (NORTHPOLE)
      ENDIF
#     endif      

#     else
      X3X2  = VX(TMP3)-VX(TMP2)
      X1X3  = VX(IA)  -VX(TMP3)
      X2X1  = VX(TMP2)-VX(IA)
      Y3Y2  = VY(TMP3)-VY(TMP2)
      Y1Y3  = VY(IA)  -VY(TMP3)
      Y2Y1  = VY(TMP2)-VY(IA)
      VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE(I)-Y3Y2*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
      VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE(I)-Y1Y3*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
      VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE(I)-Y2Y1*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
!!!      VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE(I)-Y3Y2*DLTYE(I))*H1(I1)/ART(I1)/ART1(IA)
!!!      VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE(I)-Y1Y3*DLTYE(I))*H1(I1)/ART(I1)/ART1(IA)
!!!      VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE(I)-Y2Y1*DLTYE(I))*H1(I1)/ART(I1)/ART1(IA)
#     endif
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL3,ADD_VALUES,IERR);CHKERRQ(IERR)

#     if defined (SPHERICAL)

#     if defined (NORTHPOLE)
      IF (CELL_NORTHAREA(I1)==1) THEN

        VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(XXX1*DLTYE(I)-YYY1*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IB)
        VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(XXX2*DLTYE(I)-YYY2*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IB)
        VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(XXX3*DLTYE(I)-YYY3*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IB)
!!!        VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(XXX1*DLTYE(I)-YYY1*DLTXE(I))*H1(I1)/ART(I1)/ART1(IB)
!!!        VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(XXX2*DLTYE(I)-YYY2*DLTXE(I))*H1(I1)/ART(I1)/ART1(IB)
!!!        VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(XXX3*DLTYE(I)-YYY3*DLTXE(I))*H1(I1)/ART(I1)/ART1(IB)

      ELSE
#     endif
        VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y3Y2*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
        VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y1Y3*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
        VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y2Y1*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
!!!        VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y3Y2*DLTYE(I))*H1(I1)/ART(I1)/ART1(IB)
!!!        VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y1Y3*DLTYE(I))*H1(I1)/ART(I1)/ART1(IB)
!!!        VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y2Y1*DLTYE(I))*H1(I1)/ART(I1)/ART1(IB)
#     if defined (NORTHPOLE)
      ENDIF
#     endif

#     else
      VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE(I)-Y3Y2*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
      VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE(I)-Y1Y3*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
      VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE(I)-Y2Y1*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
!!!      VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE(I)-Y3Y2*DLTYE(I))*H1(I1)/ART(I1)/ART1(IB)
!!!      VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE(I)-Y1Y3*DLTYE(I))*H1(I1)/ART(I1)/ART1(IB)
!!!      VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE(I)-Y2Y1*DLTYE(I))*H1(I1)/ART(I1)/ART1(IB)
#     endif
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL3-1,VCOL3,ADD_VALUES,IERR);CHKERRQ(IERR)

#     if defined (WET_DRY)
      ENDIF
#     endif
    ENDIF
  ENDDO

  DO I=1, NNode
    PROW1  = ALO_2_PLO_NODE(I)
    IF(PROW1 > N_VERTS) CYCLE
    PCOL1  = ALO_2_PLO_NODE(I)
    VCOL1 = 1.0D0
    CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
# endif

  DO I=1, IOBCN
    PROW1 = ALO_2_PLO_NODE(I_OBC_N(I))
    IF(PROW1>N_VERTS) CYCLE
    VCOL1 = 0.0D0
    NODE  = I_OBC_N(I)
    DO J=1, NTSN(NODE)-1
      PCOL1 = ALO_2_PLO_NODE(NBSN(NODE,J))
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
    ENDDO
    VCOL1 = 1.0D0
    PROW1 = ALO_2_PLO_NODE(NODE)
    CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

# if defined(SPHERICAL) && (NORTHPOLE)
  IF(NODE_NORTHPOLE/=0) THEN
    PROW1 = ALO_2_PLO_NODE(NODE_NORTHPOLE)
    IF(PROW1<=N_VERTS) THEN
      VCOL1 = 0.0D0
      NODE  = NODE_NORTHPOLE
      DO J=1, NTSN(NODE)-1
        PCOL1 = ALO_2_PLO_NODE(NBSN(NODE,J))
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
      ENDDO
      PROW1 = ALO_2_PLO_NODE(NODE)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
    ENDIF
  ENDIF
# endif

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)

# if defined(SPHERICAL) && (NORTHPOLE)
  COEF  = DTI**2*IFCETA**2
  DO II=1,NPCV
    I = NCEDGE_LST(II)
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)

    IF(IA==NODE_NORTHPOLE .OR. IB==NODE_NORTHPOLE) THEN
      DO J=1, 3
        IF(NV(I1,J)==IA) THEN
          TMP2 = NV(I1,J+1-INT((J+1)/4)*3)
          TMP3 = NV(I1,J+2-INT((J+2)/4)*3)
        ENDIF
      ENDDO

      PROW1 = ALO_2_PLO_NODE(IA)
      PROW2 = ALO_2_PLO_NODE(IB)

      PCOL1 = ALO_2_PLO_NODE(IA)
      PCOL2 = ALO_2_PLO_NODE(IB)
      PCOL3 = ALO_2_PLO_NODE(TMP3)

      VX1_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * COS(VX(IA)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP)) 
      VY1_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * SIN(VX(IA)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))

      VX2_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * COS(VX(IB)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP)) 
      VY2_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * SIN(VX(IB)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))

      VX3_TMP = REARTH * COS(VY(TMP3)*PI/180.0_SP) * COS(VX(TMP3)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(VY(TMP3)*PI/180.0_SP))
      VY3_TMP = REARTH * COS(VY(TMP3)*PI/180.0_SP) * SIN(VX(TMP3)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(VY(TMP3)*PI/180.0_SP))

      X3X2  = VX3_TMP-VX2_TMP
      X1X3  = VX1_TMP-VX3_TMP
      X2X1  = VX2_TMP-VX1_TMP
      Y3Y2  = VY3_TMP-VY2_TMP
      Y1Y3  = VY1_TMP-VY3_TMP
      Y2Y1  = VY2_TMP-VY1_TMP

      VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)!&
!              * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))
      VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)!&
!              * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))

      VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)!&
!              * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))
      VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)!&
!              * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))

      DLTXE_TMP = VX2_TMP-VX1_TMP
      DLTYE_TMP = VY2_TMP-VY1_TMP

      IF(IA == NODE_NORTHPOLE) THEN
       VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE_TMP-Y3Y2*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
       VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE_TMP-Y1Y3*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
       VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE_TMP-Y2Y1*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
!!!       VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE_TMP-Y3Y2*DLTYE_TMP)*H1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
!!!       VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE_TMP-Y1Y3*DLTYE_TMP)*H1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
!!!       VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE_TMP-Y2Y1*DLTYE_TMP)*H1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL3,ADD_VALUES,IERR);CHKERRQ(IERR)
      ELSE IF (IB == NODE_NORTHPOLE) THEN
       VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE_TMP-Y3Y2*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
       VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE_TMP-Y1Y3*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
       VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE_TMP-Y2Y1*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
!!!       VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE_TMP-Y3Y2*DLTYE_TMP)*H1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
!!!       VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE_TMP-Y1Y3*DLTYE_TMP)*H1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
!!!       VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE_TMP-Y2Y1*DLTYE_TMP)*H1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*PI/180.0_SP)))
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL3-1,VCOL3,ADD_VALUES,IERR);CHKERRQ(IERR)
      ENDIF
    ENDIF

  ENDDO

  IF(NODE_NORTHPOLE/=0) THEN
    PROW1 = ALO_2_PLO_NODE(NODE_NORTHPOLE)
    IF(PROW1<=N_VERTS) THEN
      VCOL1 = 1.0D0
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
    ENDIF
  ENDIF
 
  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
# endif

  CALL PETSc_SETICS_EL
  CALL PETSc_SOLVER_EL

  CALL VecScatterBegin(X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,G2L_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,G2L_EL,IERR);CHKERRQ(IERR)
  CALL VecGetArrayF90(XL_EL,XVALS_EL,IERR);CHKERRQ(IERR)

  EL = 0.0_SP
  DO I=1,N_VERTS
    IK = PLO_2_ALO_NODE(I)
    EL(IK) = XVALS_EL(I)
  ENDDO

# if defined (MULTIPROCESSOR)
    IF(PAR) CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,NTNode,1,MYID,NPROCS,EL)
    IF(PAR) CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,EL)
# endif

# if defined (WET_DRY)
  IF(WET_DRY_ON) THEN
    ELF = EL
    CALL WET_JUDGE
    IF(PAR) CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,ELF)
    EL = ELF
  ENDIF
# endif

  DO I = 1, MTElem
    EL1(I) = (EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))/3.0_SP
  END DO
  D   = H + EL
  D1  = H1 + EL1

  DO I=1, NE
    IA=IEC(I,1)
    IB=IEC(I,2)
    J1=IENODE(I,1)
    J2=IENODE(I,2)
    ELIJ=IFCETA*0.5_SP*(EL(J1)+EL(J2))

#   if defined (TWO_D_MODEL)

#   if defined (SPHERICAL)
    XTMP  = VX(J2)*TPI-VX(J1)*TPI
    XTMP1 = VX(J2)-VX(J1)
    IF(XTMP1 >  180.0_SP)THEN
      XTMP = -360.0_SP*TPI+XTMP
    ELSE IF(XTMP1 < -180.0_SP)THEN
      XTMP =  360.0_SP*TPI+XTMP
    END IF

    ADVUA(IA)=ADVUA(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
    ADVVA(IA)=ADVVA(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
    ADVUA(IB)=ADVUA(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
    ADVVA(IB)=ADVVA(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
!!!    ADVUA(IA)=ADVUA(IA)-GRAV_E(IA)*H1(IA)*ELIJ*DLTYC(I)
!!!    ADVVA(IA)=ADVVA(IA)+GRAV_E(IA)*H1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
!!!    ADVUA(IB)=ADVUA(IB)+GRAV_E(IB)*H1(IB)*ELIJ*DLTYC(I)
!!!    ADVVA(IB)=ADVVA(IB)-GRAV_E(IB)*H1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#   else
    ADVUA(IA)=ADVUA(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
    ADVVA(IA)=ADVVA(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*DLTXC(I)
    ADVUA(IB)=ADVUA(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
    ADVVA(IB)=ADVVA(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*DLTXC(I)
!!!    ADVUA(IA)=ADVUA(IA)-GRAV_E(IA)*H1(IA)*ELIJ*DLTYC(I)
!!!    ADVVA(IA)=ADVVA(IA)+GRAV_E(IA)*H1(IA)*ELIJ*DLTXC(I)
!!!    ADVUA(IB)=ADVUA(IB)+GRAV_E(IB)*H1(IB)*ELIJ*DLTYC(I)
!!!    ADVVA(IB)=ADVVA(IB)-GRAV_E(IB)*H1(IB)*ELIJ*DLTXC(I)
#   endif

#   else
    
    DO K=1, KBM1
#     if defined (SPHERICAL)
      XTMP  = VX(J2)*TPI-VX(J1)*TPI
      XTMP1 = VX(J2)-VX(J1)
      IF(XTMP1 >  180.0_SP)THEN
        XTMP = -360.0_SP*TPI+XTMP
      ELSE IF(XTMP1 < -180.0_SP)THEN
        XTMP =  360.0_SP*TPI+XTMP
      END IF

      XFLUX3(IA,K)=XFLUX3(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
      YFLUX3(IA,K)=YFLUX3(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
      XFLUX3(IB,K)=XFLUX3(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
      YFLUX3(IB,K)=YFLUX3(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
!!!      XFLUX3(IA,K)=XFLUX3(IA,K)-GRAV_E(IA)*H1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
!!!      YFLUX3(IA,K)=YFLUX3(IA,K)+GRAV_E(IA)*H1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
!!!      XFLUX3(IB,K)=XFLUX3(IB,K)+GRAV_E(IB)*H1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
!!!      YFLUX3(IB,K)=YFLUX3(IB,K)-GRAV_E(IB)*H1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#     else
      XFLUX3(IA,K)=XFLUX3(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
      YFLUX3(IA,K)=YFLUX3(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
      XFLUX3(IB,K)=XFLUX3(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
      YFLUX3(IB,K)=YFLUX3(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)
!!!      XFLUX3(IA,K)=XFLUX3(IA,K)-GRAV_E(IA)*H1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
!!!      YFLUX3(IA,K)=YFLUX3(IA,K)+GRAV_E(IA)*H1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
!!!      XFLUX3(IB,K)=XFLUX3(IB,K)+GRAV_E(IB)*H1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
!!!      YFLUX3(IB,K)=YFLUX3(IB,K)-GRAV_E(IB)*H1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)
#     endif
    ENDDO

#   endif
  ENDDO

# if defined (SPHERICAL) && (NORTHPOLE)
  DO II=1,NPE
    I=NPEDGE_LST(II)
    IA=IEC(I,1)
    IB=IEC(I,2)
    J1=IENODE(I,1)
    J2=IENODE(I,2)
    ELIJ=IFCETA*0.5_SP*(EL(J1)+EL(J2))

    DO K=1, KBM1
      VX1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * COS(VX(IENODE(I,1))*PI/180.0_SP) &
                * 2._SP /(1._SP+sin(VY(IENODE(I,1))*PI/180.0_SP))
      VY1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * SIN(VX(IENODE(I,1))*PI/180.0_SP) &
                * 2._SP /(1._SP+sin(VY(IENODE(I,1))*PI/180.0_SP))

      VX2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * COS(VX(IENODE(I,2))*PI/180.0_SP) &
                * 2._SP /(1._SP+sin(VY(IENODE(I,2))*PI/180.0_SP))
      VY2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * SIN(VX(IENODE(I,2))*PI/180.0_SP) &
                * 2._SP /(1._SP+sin(VY(IENODE(I,2))*PI/180.0_SP))

      DLTXC_TMP = VX2_TMP-VX1_TMP
      DLTYC_TMP = VY2_TMP-VY1_TMP

      IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
        XFLUX3_NP(IA,K)=XFLUX3_NP(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC_TMP/  &
                       (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
        YFLUX3_NP(IA,K)=YFLUX3_NP(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC_TMP/  &
                       (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
        XFLUX3_NP(IB,K)=XFLUX3_NP(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC_TMP/  &
                       (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
        YFLUX3_NP(IB,K)=YFLUX3_NP(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC_TMP/  &
                       (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
!!!        XFLUX3_NP(IA,K)=XFLUX3_NP(IA,K)-GRAV_E(IA)*H1(IA)*DZ1(IA,K)*ELIJ*DLTYC_TMP/  &
!!!                      (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
!!!        YFLUX3_NP(IA,K)=YFLUX3_NP(IA,K)+GRAV_E(IA)*H1(IA)*DZ1(IA,K)*ELIJ*DLTXC_TMP/  &
!!!                      (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
!!!        XFLUX3_NP(IB,K)=XFLUX3_NP(IB,K)+GRAV_E(IB)*H1(IB)*DZ1(IB,K)*ELIJ*DLTYC_TMP/  &
!!!                      (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
!!!        YFLUX3_NP(IB,K)=YFLUX3_NP(IB,K)-GRAV_E(IB)*H1(IB)*DZ1(IB,K)*ELIJ*DLTXC_TMP/  &
!!!                      (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
      ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
        XFLUX3_NP(IA,K)=XFLUX3_NP(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC_TMP/   &
                      (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
        YFLUX3_NP(IA,K)=YFLUX3_NP(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC_TMP/   &
                      (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
!!!        XFLUX3_NP(IA,K)=XFLUX3_NP(IA,K)-GRAV_E(IA)*H1(IA)*DZ1(IA,K)*ELIJ*DLTYC_TMP/   &
!!!                      (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
!!!        YFLUX3_NP(IA,K)=YFLUX3_NP(IA,K)+GRAV_E(IA)*H1(IA)*DZ1(IA,K)*ELIJ*DLTXC_TMP/   &
!!!                      (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
      ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN
        XFLUX3_NP(IB,K)=XFLUX3_NP(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC_TMP/  &
                      (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
        YFLUX3_NP(IB,K)=YFLUX3_NP(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC_TMP/  &
                      (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
!!!        XFLUX3_NP(IB,K)=XFLUX3_NP(IB,K)+GRAV_E(IB)*H1(IB)*DZ1(IB,K)*ELIJ*DLTYC_TMP/  &
!!!                      (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
!!!        YFLUX3_NP(IB,K)=YFLUX3_NP(IB,K)-GRAV_E(IB)*H1(IB)*DZ1(IB,K)*ELIJ*DLTXC_TMP/  &
!!!                      (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
      END IF
    ENDDO

  ENDDO
# endif

  DO I=1, MElem

# if defined (WET_DRY)
  IF(ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1) THEN
# endif
    
# if defined (TWO_D_MODEL)

#   if defined (SPHERICAL)
    UAF(I) = UA(I)*DT1(I)/D1(I) - DTI*ADVUA(I)/ART(I)/D1(I) + DTI*( -WUSURF(I) + WUBOT(I) )/D1(I) - BEDF*UA_N(I)*DT1(I)/D1(I) 
    VAF(I) = VA(I)*DT1(I)/D1(I) - DTI*ADVVA(I)/ART(I)/D1(I) + DTI*( -WVSURF(I) + WVBOT(I) )/D1(I) - BEDF*VA_N(I)*DT1(I)/D1(I)
!!!    UAF(I) = UA(I) - DTI*ADVUA(I)/ART(I)/H1(I) 
!!!    VAF(I) = VA(I) - DTI*ADVVA(I)/ART(I)/H1(I) 
#   else
    UAF(I) = UA(I)*DT1(I)/D1(I) - DTI*ADVUA(I)/ART(I)/D1(I) + DTI*( -WUSURF(I) + WUBOT(I) )/D1(I) - BEDF*UA_N(I)*DT1(I)/D1(I)
    VAF(I) = VA(I)*DT1(I)/D1(I) - DTI*ADVVA(I)/ART(I)/D1(I) + DTI*( -WVSURF(I) + WVBOT(I) )/D1(I) - BEDF*VA_N(I)*DT1(I)/D1(I)
!!!    UAF(I) = UA(I) - DTI*ADVUA(I)/ART(I)/H1(I)
!!!    VAF(I) = VA(I) - DTI*ADVVA(I)/ART(I)/H1(I)
#   endif

# else

    DO K=1, KBM1
#     if defined (SPHERICAL)
#     if defined (NORTHPOLE)
      IF(CELL_NORTHAREA(I)==1) THEN      
        U_TMP = -V(I,K)*COS(XC(I)*PI/180.0_SP)-U(I,K)*SIN(XC(I)*PI/180.0_SP)
        V_TMP = -V(I,K)*SIN(XC(I)*PI/180.0_SP)+U(I,K)*COS(XC(I)*PI/180.0_SP)  
        UF_TMP = U_TMP*DT1(I)/D1(I) - DTI*XFLUX3_NP(I,K)/ART(I)/D1(I)/DZ1(I,K)
        VF_TMP = V_TMP*DT1(I)/D1(I) - DTI*YFLUX3_NP(I,K)/ART(I)/D1(I)/DZ1(I,K)
!!!        UF_TMP = U_TMP - DTI*XFLUX3_NP(I,K)/ART(I)/H1(I)/DZ1(I,K)
!!!        VF_TMP = V_TMP - DTI*YFLUX3_NP(I,K)/ART(I)/H1(I)/DZ1(I,K)

        UF(I,K) = UF_TMP
        VF(I,K) = VF_TMP  
      ELSE
#     endif
       UF(I,K) = U(I,K)*DT1(I)/D1(I) - DTI*XFLUX3(I,K)/ART(I)/D1(I)/DZ1(I,K) - BEDF*U_N(I,K)*DT1(I)/D1(I)
       VF(I,K) = V(I,K)*DT1(I)/D1(I) - DTI*YFLUX3(I,K)/ART(I)/D1(I)/DZ1(I,K) - BEDF*V_N(I,K)*DT1(I)/D1(I)
!!!       UF(I,K) = U(I,K) - DTI*XFLUX3(I,K)/ART(I)/H1(I)/DZ1(I,K) 
!!!       VF(I,K) = V(I,K) - DTI*YFLUX3(I,K)/ART(I)/H1(I)/DZ1(I,K) 
#     if defined (NORTHPOLE)
      ENDIF
#     endif

#     else
      UF(I,K) = U(I,K)*DT1(I)/D1(I) - DTI*XFLUX3(I,K)/ART(I)/D1(I)/DZ1(I,K) - BEDF*U_N(I,K)*DT1(I)/D1(I)
      VF(I,K) = V(I,K)*DT1(I)/D1(I) - DTI*YFLUX3(I,K)/ART(I)/D1(I)/DZ1(I,K) - BEDF*V_N(I,K)*DT1(I)/D1(I)
!!!      UF(I,K) = U(I,K) - DTI*XFLUX3(I,K)/ART(I)/H1(I)/DZ1(I,K)
!!!      VF(I,K) = V(I,K) - DTI*YFLUX3(I,K)/ART(I)/H1(I)/DZ1(I,K)
#     endif
    ENDDO
    
# endif

# if defined (WET_DRY)
  ELSE
#   if defined (TWO_D_MODEL)
    UAF(I) = 0.0_SP
    VAF(I) = 0.0_SP
#   else
    DO K=1, KBM1
      UF(I,K) = 0.0_SP
      VF(I,K) = 0.0_SP
    ENDDO  
#   endif
  ENDIF
# endif

  ENDDO

# if defined (TWO_D_MODEL)
  UAF = UAF-CC_SPONGE*UAF
  VAF = VAF-CC_SPONGE*VAF
# endif

  RETURN
  END SUBROUTINE SEMI_IMPLICIT_EL

  SUBROUTINE UPDATE_EL

  USE LIMS
  USE CONTROL 
  USE MOD_OBCS
  USE ALL_VARS
  USE BCS
  USE MOD_PETSC, ONLY: A_EL,X_EL,B_EL,BL_EL,XL_EL,N_VERTS,L2G_EL,G2L_EL,XVALS_EL,ALO_2_PLO_NODE,PLO_2_ALO_NODE,PETSc_SETICS_EL,PETSc_SOLVER_EL
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
# if defined (SPHERICAL)   
  USE MOD_SPHERICAL
# if defined(NORTHPOLE)
  USE MOD_NORTHPOLE,   ONLY :  NPCV,  NCEDGE_LST, NODE_NORTHPOLE
# endif
# endif   

  IMPLICIT NONE
  INTEGER  II, I, J, JJ, J1, J2, IK, K, I1, IA, IB, TMP2, TMP3
  INTEGER  NNZ2, LN, LOC
  INTEGER  PETSc_POS, NODE, PROW1, PROW2, PCOL1, PCOL2, PCOL3
  REAL(SP) DIJ, DIJ1, UIJ, VIJ, EXFLUX, COEF
  REAL(SP) Y3Y2,Y1Y3,Y2Y1,X3X2,X1X3,X2X1
  REAL(SP) :: XFLUX(0:NTNode), XFLUX_TMP(0:MTElem), YFLUX_TMP(0:MTElem)
  REAL(SP) :: UF_AVG(0:MTElem), VF_AVG(0:MTElem)

  REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP,FLUX_TMP
  REAL(SP) :: DLTXE_TMP,DLTYE_TMP
  REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP

  PetscReal :: STERM, VCOL1,VCOL2,VCOL3
  PetscInt  :: IERR
  PetscScalar :: ZERO1   =  0.0D0

  UF_AVG = 0.0_SP
  VF_AVG = 0.0_SP
# if defined(MULTIPROCESSOR)
   IF(PAR) CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,UF,VF)
# endif
  DO I=1,MTElem
    UF_AVG(I) = SUM(UF(I,1:KBM1)*DZ1(I,1:KBM1))
    VF_AVG(I) = SUM(VF(I,1:KBM1)*DZ1(I,1:KBM1))
  END DO

  CALL VecSet(BL_EL,ZERO,IERR);CHKERRQ(IERR)
  CALL VecSet(B_EL,ZERO,IERR);CHKERRQ(IERR)

  XFLUX   = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)
    DIJ=(DT(NV(I1,1))+DT(NV(I1,2))+DT(NV(I1,3)))/3.0_SP
    DIJ1=(D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP

    UIJ=UF_AVG(I1)
    VIJ=VF_AVG(I1)
    EXFLUX=(1.0_SP-IFCETA)*DIJ*(-UA(I1)*DLTYE(I)+VA(I1)*DLTXE(I))+IFCETA*DIJ1*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
    XFLUX(IA)=XFLUX(IA)-EXFLUX
    XFLUX(IB)=XFLUX(IB)+EXFLUX
  ENDDO
# if defined (SPHERICAL) && (NORTHPOLE)
  IF(NODE_NORTHPOLE/=0) XFLUX(NODE_NORTHPOLE) = 0.0_SP

  DO II=1,NPCV
    I = NCEDGE_LST(II)
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)

    DIJ = (DT(NV(I1,1))+DT(NV(I1,2))+DT(NV(I1,3)))/3.0_SP
    DIJ1= (D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP
    UIJ = -VA(I1)*COS(XC(I1)*PI/180.0_SP)-UA(I1)*SIN(XC(I1)*PI/180.0_SP)
    VIJ = -VA(I1)*SIN(XC(I1)*PI/180.0_SP)+UA(I1)*COS(XC(I1)*PI/180.0_SP)
    UIJ_TMP = -VF_AVG(I1)*COS(XC(I1)*PI/180.0_SP)-UF_AVG(I1)*SIN(XC(I1)*PI/180.0_SP)
    VIJ_TMP = -VF_AVG(I1)*SIN(XC(I1)*PI/180.0_SP)+UF_AVG(I1)*COS(XC(I1)*PI/180.0_SP)

    IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
      VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)!&
!                * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))
      VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)!&
!                * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))

      VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)!&
!                * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))
      VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)!&
!                * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))

      DLTXE_TMP = VX2_TMP-VX1_TMP
      DLTYE_TMP = VY2_TMP-VY1_TMP

      EXFLUX_TMP =(1.0_SP-IFCETA)*DIJ*(-UIJ*DLTYE_TMP+VIJ*DLTXE_TMP)+IFCETA*DIJ1*(-UIJ_TMP*DLTYE_TMP+VIJ_TMP*DLTXE_TMP)
    END IF

    IF(IA == NODE_NORTHPOLE) THEN
      XFLUX(IA) = XFLUX(IA)-EXFLUX_TMP
    ELSE IF(IB == NODE_NORTHPOLE)THEN
      XFLUX(IB) = XFLUX(IB)+EXFLUX_TMP
    END IF

  END DO
# endif

! this part may need be changed to semi-implicit form 
  XFLUX = XFLUX+(QEVAP3-QPREC3)*ROFVROS*ART1

  IF(IBFW > 0)THEN
    DO I=1,NNode
      DO J=1,IBFW
        IF(I == NODE_BFW(J))THEN
          XFLUX(I)=XFLUX(I)+BFWDIS3(J)*ROFVROS*ART1(I)
        END IF
      END DO
    END DO
  END IF

  IF(NUMQBC >= 1) THEN
    IF(INFLOW_TYPE == 'node') THEN
      DO J=1,NUMQBC
        JJ=INODEQ(J)
        XFLUX(JJ)=XFLUX(JJ)-QDIS(J)
      END DO
    ELSE IF(INFLOW_TYPE == 'edge') THEN
      DO J=1,NUMQBC
        J1=N_ICELLQ(J,1)
        J2=N_ICELLQ(J,2)
        XFLUX(J1)=XFLUX(J1)-QDIS(J)*RDISQ(J,1)
        XFLUX(J2)=XFLUX(J2)-QDIS(J)*RDISQ(J,2)
      END DO
    END IF
  END IF

  IF(IOBCN > 0) THEN
    DO I=1,IOBCN
      XFLUX_OBCN(I)=XFLUX(I_OBC_N(I))
      XFLUX(I_OBC_N(I)) = 0.0_SP
    END DO
  END IF

  if(.false.) then
  DO I=1, NNode
    PETSc_POS = ALO_2_PLO_NODE(I)
    IF (PETSc_POS > N_VERTS) CYCLE
    IF( ISONB(I)/=2 ) THEN
      STERM = ET(I) - DTI/ART1(I)*XFLUX(I)
    ELSE
      STERM = ELF(I)
    ENDIF
    CALL VecSetValues(BL_EL,1,PETSc_POS-1,STERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

  CALL VecScatterBegin(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)

  CALL MatZeroEntries(A_EL,IERR);CHKERRQ(IERR)

  DO I=1, NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)

    IF((ISONB(IA)+ISONB(IB))<4) THEN
      DO J=1, 3
        IF(NV(I1,J)==IA) THEN
          TMP2 = NV(I1,J+1-INT((J+1)/4)*3)
          TMP3 = NV(I1,J+2-INT((J+2)/4)*3)
        ENDIF
      ENDDO

      PROW1 = ALO_2_PLO_NODE(IA)
      PROW2 = ALO_2_PLO_NODE(IB)

      PCOL1 = ALO_2_PLO_NODE(IA)
      PCOL2 = ALO_2_PLO_NODE(IB)
      PCOL3 = ALO_2_PLO_NODE(TMP3)

      VCOL1 = DTI*IFCETA*(-VF_AVG(I1)*DLTXE(I)+UF_AVG(I1)*DLTYE(I))/ART1(IA)/3.0_SP
      VCOL2 = -DTI*IFCETA*(-VF_AVG(I1)*DLTXE(I)+UF_AVG(I1)*DLTYE(I))/ART1(IB)/3.0_SP
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL1-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)

      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)

      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL3-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
    ENDIF
  ENDDO

  DO I=1, NNode
    PROW1  = ALO_2_PLO_NODE(I)
    IF(PROW1 > N_VERTS) CYCLE
    PCOL1  = ALO_2_PLO_NODE(I)
    VCOL1 = 1.0D0
    CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
# if defined (MULTIPROCESSOR)
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
# endif

  DO I=1, IOBCN
    PROW1 = ALO_2_PLO_NODE(I_OBC_N(I))
    IF(PROW1>N_VERTS) CYCLE
    VCOL1 = 0.0D0
    NODE  = I_OBC_N(I)
    DO J=1, NTSN(NODE)-1
      PCOL1 = ALO_2_PLO_NODE(NBSN(NODE,J))
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
    ENDDO
    VCOL1 = 1.0D0
    PROW1 = ALO_2_PLO_NODE(NODE)
    CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

# if defined(SPHERICAL) && (NORTHPOLE)
  IF (NODE_NORTHPOLE/=0) THEN
    PROW1 = ALO_2_PLO_NODE(NODE_NORTHPOLE)
    IF(PROW1<=N_VERTS) THEN
      VCOL1 = 0.0D0
      NODE  = NODE_NORTHPOLE
      DO J=1, NTSN(NODE)-1
        PCOL1 = ALO_2_PLO_NODE(NBSN(NODE,J))
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
      ENDDO
      PROW1 = ALO_2_PLO_NODE(NODE)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
    ENDIF
  ENDIF
# endif

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)

# if defined(SPHERICAL) && (NORTHPOLE) 
  DO II=1,NPCV
    I = NCEDGE_LST(II)
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)

    IF(IA==NODE_NORTHPOLE .OR. IB==NODE_NORTHPOLE) THEN
      DO J=1, 3
        IF(NV(I1,J)==IA) THEN
          TMP2 = NV(I1,J+1-INT((J+1)/4)*3)
          TMP3 = NV(I1,J+2-INT((J+2)/4)*3)
        ENDIF
      ENDDO

      PROW1 = ALO_2_PLO_NODE(IA)
      PROW2 = ALO_2_PLO_NODE(IB)

      PCOL1 = ALO_2_PLO_NODE(IA)
      PCOL2 = ALO_2_PLO_NODE(IB)
      PCOL3 = ALO_2_PLO_NODE(TMP3)

      UIJ_TMP = -VF_AVG(I1)*COS(XC(I1)*PI/180.0_SP)-UF_AVG(I1)*SIN(XC(I1)*PI/180.0_SP)
      VIJ_TMP = -VF_AVG(I1)*SIN(XC(I1)*PI/180.0_SP)+UF_AVG(I1)*COS(XC(I1)*PI/180.0_SP)

      VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))
      VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))

      VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))
      VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)&
              * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))

      DLTXE_TMP = VX2_TMP-VX1_TMP
      DLTYE_TMP = VY2_TMP-VY1_TMP

      IF(IA == NODE_NORTHPOLE) THEN
        VCOL1 = DTI*IFCETA*(-VIJ_TMP*DLTXE_TMP+UIJ_TMP*DLTYE_TMP)/ART1(IA)/3.0_SP
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      ELSE IF (IB == NODE_NORTHPOLE) THEN
        VCOL2 = -DTI*IFCETA*(-VIJ_TMP*DLTXE_TMP+UIJ_TMP*DLTYE_TMP)/ART1(IB)/3.0_SP
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL1-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL3-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
      ENDIF
    ENDIF

  ENDDO

  IF(NODE_NORTHPOLE/=0) THEN
    PROW1 = ALO_2_PLO_NODE(NODE_NORTHPOLE)
    IF(PROW1<=N_VERTS) THEN
      VCOL1 = 1.0D0
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
    ENDIF
  ENDIF

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
# endif

  CALL PETSc_SETICS_EL
  CALL PETSc_SOLVER_EL

  CALL VecScatterBegin(X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,G2L_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,G2L_EL,IERR);CHKERRQ(IERR)
  CALL VecGetArrayF90(XL_EL,XVALS_EL,IERR);CHKERRQ(IERR)

  EL = 0.0_SP
  DO I=1,N_VERTS
    IK = PLO_2_ALO_NODE(I)
    EL(IK) = XVALS_EL(I)
  ENDDO

# if defined (MULTIPROCESSOR)
    IF(PAR) CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,NTNode,1,MYID,NPROCS,EL)
    IF(PAR) CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,EL)
# endif
  endif

  UA = UF_AVG
  VA = VF_AVG
!  DO I = 1, MTElem
!    EL1(I) = (EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))/3.0_SP
!  END DO
!  D   = H + EL
!  D1  = H1 + EL1
  DTFA = D

  RETURN
  END SUBROUTINE UPDATE_EL

  SUBROUTINE UV2D_SBD
  USE ALL_VARS
  USE BCS
  USE MOD_OBCS
  IMPLICIT NONE

  INTEGER  :: I,J,I1,J1,J2
  REAL(SP) :: ALPHA1,UNTMP,VNTMP,UI,VI

  UA_N = 0.0_SP
  VA_N = 0.0_SP

  DO I=1, MElem

    IF(ISBCE(I) == 1) THEN
      ALPHA1=ALPHA(I)
      IF(NUMQBC > 0) THEN
        IF(INFLOW_TYPE == 'node') THEN
          DO J=1,NUMQBC
            I1=INODEQ(J)
            J1=NBVE(I1,1)
            J2=NBVE(I1,NTVE(I1))
            IF((I == J1).OR.(I == J2)) THEN
              UNTMP=UA(I)*COS(ANGLEQ(J))+VA(I)*SIN(ANGLEQ(J))
              IF(UNTMP<0.0_SP) THEN
                UA_N(I) = UNTMP*COS(ANGLEQ(J))
                VA_N(I) = UNTMP*SIN(ANGLEQ(J))
              ENDIF
              GOTO 21
            END IF
          END DO
        ELSE IF(INFLOW_TYPE == 'edge') THEN
          DO J=1,NUMQBC
            J1=ICELLQ(J)
            IF(I == J1) THEN
              UNTMP=UA(I)*COS(ANGLEQ(J))+VA(I)*SIN(ANGLEQ(J))
              IF(UNTMP<0.0_SP) THEN
                UA_N(I) = UNTMP*COS(ANGLEQ(J))
                VA_N(I) = UNTMP*SIN(ANGLEQ(J))
              ENDIF
              GOTO 21
            END IF
          END DO
        END IF
      END IF

      UNTMP = UA(I)*COS(ALPHA1)+VA(I)*SIN(ALPHA1)
      IF(UNTMP>0.0_SP) THEN
        UA_N(I)= UNTMP*COS(ALPHA1)
        VA_N(I)= UNTMP*SIN(ALPHA1)
      ENDIF

21  CONTINUE
    END IF
  END DO

  RETURN
  END SUBROUTINE UV2D_SBD

  SUBROUTINE UV3D_SBD
  USE ALL_VARS
  USE BCS
  USE MOD_OBCS
  IMPLICIT NONE

  INTEGER  :: I,J,K,I1,J1,J2
  REAL(SP) :: ALPHA1,UNTMP,VNTMP,UI,VI

  U_N = 0.0_SP
  V_N = 0.0_SP

  DO I=1, MElem
    DO K=1, KBM1
      IF(ISBCE(I).EQ.1) THEN
        ALPHA1=ALPHA(I)
        IF(NUMQBC.GE.1) THEN
          IF(INFLOW_TYPE .EQ. 'node') THEN
            DO J=1,NUMQBC
              I1=INODEQ(J)
              J1=NBVE(I1,1)
              J2=NBVE(I1,NTVE(I1))
              IF((I.EQ.J1).OR.(I.EQ.J2)) THEN
                UNTMP=U(I,K)*COS(ANGLEQ(J))+V(I,K)*SIN(ANGLEQ(J))
                IF(UNTMP<0.0_SP) THEN
                  U_N(I,K) = UNTMP*COS(ANGLEQ(J))
                  V_N(I,K) = UNTMP*SIN(ANGLEQ(J))
                ENDIF
                GOTO 31
              ENDIF
            ENDDO
          ELSE IF(INFLOW_TYPE .EQ. 'edge') THEN
            DO J=1,NUMQBC
              J1=ICELLQ(J)
              IF(I.EQ.J1) THEN
                UNTMP=U(I,K)*COS(ANGLEQ(J))+V(I,K)*SIN(ANGLEQ(J))
                IF(UNTMP<0.0_SP) THEN
                  U_N(I,K) = UNTMP*COS(ANGLEQ(J))
                  V_N(I,K) = UNTMP*SIN(ANGLEQ(J))
                ENDIF
                GOTO 31
              ENDIF
            ENDDO
          ELSE
             PRINT*, 'inflow_type not correct'
             CALL PSTOP
          ENDIF
        ENDIF

        UNTMP=U(I,K)*COS(ALPHA1)+V(I,K)*SIN(ALPHA1)
        IF(UNTMP>0.0_SP) THEN
          U_N(I,K) = UNTMP*COS(ALPHA1)
          V_N(I,K) = UNTMP*SIN(ALPHA1)
        ENDIF

31    CONTINUE
      ENDIF
    ENDDO
  ENDDO

  RETURN
  END SUBROUTINE UV3D_SBD

# endif
END MODULE MOD_SEMI_IMPLICIT
