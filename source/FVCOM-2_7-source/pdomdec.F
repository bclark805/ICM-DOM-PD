!==============================================================================|
!  SET UP LOCAL PHYSICAL DOMAIN (CONNECTIVITY/MESH)                            |
!==============================================================================|

   SUBROUTINE PDOMDEC

!==============================================================================!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR  
#  endif
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
   IMPLICIT NONE
   INTEGER I,EGL,J,IERR,I1,I2,N_SPONGE
   REAL(SP), ALLOCATABLE :: CORRG(:),CORR(:)
   REAL(SP), ALLOCATABLE :: R_SPG(:),C_SPG(:) 
   INTEGER, ALLOCATABLE  :: N_SPG(:)
   REAL(SP)  TEMP,DTMP,C_SPONGE
#  if defined (SPHERICAL)
   INTEGER K,ITMP
   REAL(DP) VX1,VY1,VX2,VY2,VX3,VY3,EVX12,EVX13,EVX23,&
            EVY12,EVY13,EVY23,EVXY,VX12,VY12,VX23,VY23,VX31,VY31,&
	    X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP
#  endif

!==============================================================================|
!  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
!  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
!  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
!==============================================================================|

   IF(SERIAL) NV = NVG

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,MElem
       EGL = EGID(I)
       NV(I,1:4) = NLID(NVG(EGID(I),1:4))
     END DO
     DO I=1,MHE
       EGL = HE_LST(I)
       NV(I+MElem,1:4) = NLID_X(NVG(EGL,1:4))
     END DO
   END IF
#  endif

!==============================================================================|
!   SET UP LOCAL MESH (HORIZONTAL COORDINATES)                                 |
!==============================================================================|


!--------------READ IN X AND Y GLOBAL COORDINATES AT NODES---------------------!

   ALLOCATE(XG(0:NNodeGL),YG(0:NNodeGL)) ; XG = 0.0_SP ; YG = 0.0_SP
   DO I=1,NNodeGL
     READ(INGRD,*)J,XG(I),YG(I)
#  if defined (SPHERICAL)
     IF(XG(I) < 0.0) XG(I) = XG(I) + 360.0
!     IF(YG(I) < 0.0) YG(I) = YG(I) + 360.0
#  endif
   END DO
   CLOSE(INGRD)

!--------------CALCULATE GLOBAL MINIMUMS AND MAXIMUMS--------------------------!

#  if defined (SPHERICAL)
   VXMIN = 0.0_SP ; VXMAX = MAXVAL(XG(1:NNodeGL))
   VYMIN = 0.0_SP ; VYMAX = MAXVAL(YG(1:NNodeGL))
#  else
   VXMIN = MINVAL(XG(1:NNodeGL)) ; VXMAX = MAXVAL(XG(1:NNodeGL))
   VYMIN = MINVAL(YG(1:NNodeGL)) ; VYMAX = MAXVAL(YG(1:NNodeGL))
#  endif

!--------------SHIFT GRID TO UPPER RIGHT CARTESIAN-----------------------------!

   XG = XG - VXMIN
   YG = YG - VYMIN
   XG(0) = 0.0_SP ; YG(0) = 0.0_SP

!--------------CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES----------------!

   ALLOCATE(XCG(0:MElemGL),YCG(0:MElemGL)) ; XCG = 0.0_SP ; YCG = 0.0_SP
   DO I=1,MElemGL   
     XCG(I)  = (XG(NVG(I,1)) + XG(NVG(I,2)) + XG(NVG(I,3)))/3.0_SP
     YCG(I)  = (YG(NVG(I,1)) + YG(NVG(I,2)) + YG(NVG(I,3)))/3.0_SP
   END DO
# if defined (SPHERICAL)
   DO I=1,MElemGL
     VX1=XG(NVG(I,1))
     VY1=YG(NVG(I,1))
     VX2=XG(NVG(I,2))
     VY2=YG(NVG(I,2))
     VX3=XG(NVG(I,3))
     VY3=YG(NVG(I,3))

     DO 56 K=1,1000000

       EVX12=VX2-VX1
       EVX13=VX3-VX1
       EVX23=VX3-VX2

       IF(EVX12 >  180.0_SP)THEN
         EVX12 = -360.0_SP+EVX12
       ELSE IF(EVX12 < -180.0_SP)THEN
         EVX12 =  360.0_SP+EVX12
       END IF
       IF(EVX13 >  180.0_SP)THEN
	 EVX13 = -360.0_SP+EVX13
       ELSE IF(EVX13 < -180.0_SP)THEN
	 EVX13 =  360.0_SP+EVX13
       END IF
       IF(EVX23 >  180.0_SP)THEN
         EVX23 = -360.0_SP+EVX23
       ELSE IF(EVX23 < -180.0_SP)THEN
         EVX23 =  360.0_SP+EVX23
       END IF

       EVX12=ABS(EVX12)
       EVX13=ABS(EVX13)
       EVX23=ABS(EVX23)

       EVY12=ABS(VY2-VY1)
       EVY13=ABS(VY3-VY1)
       EVY23=ABS(VY3-VY2)

       EVXY=1.E-10_SP

       IF((EVX12 < EVXY) .AND.(EVX13 < EVXY) .AND. (EVX23 < EVXY) &
          .AND.(EVY12 < EVXY) .AND. (EVY13 < EVXY)                &
          .AND.(EVY23 < EVXY))THEN
         XCG(I)=VX1
         YCG(I)=VY1
         GOTO 57
       ELSE
         CALL ARCC(VX1,VY1,VX2,VY2,VX12,VY12)
         CALL ARCC(VX2,VY2,VX3,VY3,VX23,VY23)
         CALL ARCC(VX3,VY3,VX1,VY1,VX31,VY31)

         VX1=VX12
         VY1=VY12
         VX2=VX23
         VY2=VY23
         VX3=VX31
         VY3=VY31
       END IF
56     CONTINUE
57     CONTINUE
     END DO
# endif

     XCG(0) = 0.0_SP ; YCG(0) = 0.0_SP


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!

     IF(SERIAL)THEN
       VX = XG
       VY = YG
     END IF

#  if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,NNode
         VX(I) = XG(NGID(I))
         VY(I) = YG(NGID(I))
       END DO

       DO I=1,NHN
         VX(I+NNode) = XG(HN_LST(I))
         VY(I+NNode) = YG(HN_LST(I))
       END DO
     END IF
# endif

!==============================================================================|
!   SET UP LOCAL MESH (BATHYMETRIC DEPTH)                                      |
!==============================================================================|

!--------------READ IN BATHYMETRY----------------------------------------------!

     ALLOCATE(HG(0:NNodeGL))  ; HG = 0.0_SP
     DO I=1,NNodeGL
       READ(INDEP,*) TEMP,TEMP,HG(I)
     END DO
     CLOSE(INDEP)


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!

     IF(SERIAL) H = HG

#  if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,NNode
         H(I)   = HG(NGID(I))
! B Clark, get rid of places high in elevation
    !    IF(H(I) < -3.0) THEN
     !          H(I) = -3.0
      !  END IF
        END DO
       DO I=1,NHN
         H(I+NNode) = HG(HN_LST(I))
       END DO
     END IF

#  endif

!--------------CALCULATE EXTREMUMS---------------------------------------------!

     HMAX = MAXVAL(ABS(HG(1:NNodeGL)))
     HMIN = MINVAL(HG(1:NNodeGL))

!==============================================================================|
!   SET UP LOCAL CORIOLIS FORCE                                                |
!==============================================================================|

!--------------READ IN CORIOLIS PARAMETER--------------------------------------!

     ALLOCATE(CORRG(0:NNodeGL))  ; CORRG = 0.0_SP
!  MHB:  ADJUST FOR DIFFERENT CORIOLIS FILE FORMAT
#  if defined (SPHERICAL)
     CORRG=YG
#  else
     IF(CASENAME == "mhb")THEN
       DO I=1,NNodeGL
         READ(INCOR,*) TEMP,CORRG(I)
       END DO
     ELSE
       DO I=1,NNodeGL
         READ(INCOR,*) TEMP,TEMP,CORRG(I)
       END DO
     END IF
     CLOSE(INCOR)
# endif

!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
     ALLOCATE(CORR(0:NTNode)) ; CORR = 0.0_SP
     IF(SERIAL) CORR = CORRG

#  if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,NNode
         CORR(I) = CORRG(NGID(I))
       END DO
       DO I=1,NHN
         CORR(I+NNode) = CORRG(HN_LST(I))
       END DO
     END IF
#  endif

!==============================================================================|
!   COMPUTE FACE CENTER VALUES FOR GRID, DEPTH, AND CORIOLIS PARAMETER         |
!==============================================================================|

# if defined (SPHERICAL)
     IF(SERIAL) XC = XCG; YC = YCG

#  if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,MElem
         XC(I) = XCG(EGID(I))
         YC(I) = YCG(EGID(I))
       ENDDO
       DO I=1,MHE
         XC(I+MElem) = XCG(HE_LST(I))
         YC(I+MElem) = YCG(HE_LST(I))
       END DO
     END IF
#  endif
     COR = YC
     DO I=1,MTElem
       H1(I)  = SUM( H(NV(I,1:3)))/3.0_SP
       COR(I) = 2.*7.292e-5_SP*SIN(COR(I)*2.0_SP*3.14159_SP/360.0_SP)
     END DO

# else
     DO I=1,MTElem
!       XC(I)  = SUM(VX(NV(I,1:3)))/3.0
       XC(I)  = (VX(NV(I,1)) + VX(NV(I,2)) + VX(NV(I,3)))/3.0_SP
       YC(I)  = (VY(NV(I,1)) + VY(NV(I,2)) + VY(NV(I,3)))/3.0_SP
!       YC(I)  = SUM(VY(NV(I,1:3)))/3.0
       H1(I)  = SUM( H(NV(I,1:3)))/3.0_SP
       COR(I) = CORR(NV(I,1)) + CORR(NV(I,2)) + CORR(NV(I,3))
       COR(I) = COR(I)/3.0_SP
!       COR(I) = SUM(CORR(NV(I,1:3)))/3.0
       COR(I) = 2.*7.292e-5_SP*SIN(COR(I)*2.0_SP*3.14159_SP/360.0_SP)
     END DO
#  endif

!==============================================================================|
!   COMPUTE GRAVITY VARIED WITH LATITUDE                                       |
!==============================================================================|

     ALLOCATE(GRAV_N(0:NTNode),GRAV_E(0:MTElem))
# if defined (SPHERICAL)
     DO I=1,NTNode
       GRAV_N(I) = 9.78032+0.005172*SIN(VY(I))**2-0.00006*SIN(2.0*VY(I))**2
     END DO  
     DO I=1,MTElem
       GRAV_E(I) = 9.78032+0.005172*SIN(YC(I))**2-0.00006*SIN(2.0*YC(I))**2
     END DO  
# else
     GRAV_N = GRAV
     GRAV_E = GRAV
# endif            

!==============================================================================|
!   COMPUTE SPONGE LAYER FOR OPEN BOUNDARY DAMPING                             |
!==============================================================================|

!--READ NUMBER OF SPONGE NODES AND ALLOCATE ARRAYS-----------------------------|

     READ(INSPO,*) N_SPONGE
     IF(N_SPONGE > 0 )THEN

     ALLOCATE( N_SPG(N_SPONGE) , R_SPG(N_SPONGE) , C_SPG(N_SPONGE) )

!--READ IN INDICES OF SPONGE NODES --------------------------------------------|

     DO I=1,N_SPONGE
       READ(INSPO,*) N_SPG(I),R_SPG(I),C_SPG(I)
     END DO
     CLOSE(INSPO)



!--SET SPONGE PARAMETERS-------------------------------------------------------|

     CC_SPONGE = 0.0_SP

     DO I=1,MTElem
       DO I1=1,N_SPONGE
         I2=N_SPG(I1)
# if defined (SPHERICAL)
         X1_DP=XC(I)
         Y1_DP=YC(I)
         X2_DP=XG(I2)
         Y2_DP=YG(I2)
         CALL ARC(X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP)
         DTMP=DTMP_DP/R_SPG(I1)
# else
         DTMP=(XC(I)-XG(I2))**2+(YC(I)-YG(I2))**2
         DTMP=SQRT(DTMP)/R_SPG(I1)
# endif
         IF(DTMP <= 1.) THEN
           C_SPONGE=C_SPG(I1)*(1.-DTMP)
           CC_SPONGE(I)=MAX(C_SPONGE,CC_SPONGE(I))
         END IF
       END DO
     END DO

     DEALLOCATE(N_SPG,R_SPG,C_SPG)

   END IF !! N_SPONGE > 0

   IF(MSR)WRITE(IPT,*)'!  # SPONGE LAYER SET BY :',N_SPONGE

!==============================================================================|
!   WRITE TO SMS GRID FILE WHILE GLOBAL VALUES EXIST                           |
!==============================================================================|

   IF(MSR)THEN
     WRITE(IOSMSD,*)'scat2d'
     WRITE(IOSMSD,*)'xyd ',NNodeGL,' dep ',1,' dep '
     DO I=1,NNodeGL
       WRITE(IOSMSD,*) XG(I),YG(I),HG(I)
     END DO
     CLOSE(IOSMSD)
   END IF
   DEALLOCATE(CORR,CORRG)

   RETURN
   END SUBROUTINE PDOMDEC
!==============================================================================|
