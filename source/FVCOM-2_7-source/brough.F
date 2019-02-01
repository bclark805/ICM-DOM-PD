!==============================================================================|
!  Calculate Bottom Drag Coefficient based on Bottom Roughness                 !
!   note:                                                                      !
!   when the log function derived from the constant stress log-viscous         !
!   layer is applied to an estuary, if the value of z0 is close to             !
!   (zz(kbm1)-z(kb)*dt1, drag coefficient "cbc" could become a huge            !
!   number due to near-zero value of alog function. In our application         !
!   we simply cutoff at cbc=0.005. One could adjust this cutoff value          !
!   based on observations or his or her experiences.                           !   
!   CALCULATES:   WUBOT(MElem), WVBOT(MElem) : BOTTOM SHEAR STRESSES           !
!==============================================================================|

   SUBROUTINE BOTTOM_ROUGHNESS

!==============================================================================!
   USE ALL_VARS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif      
   IMPLICIT NONE
   INTEGER :: I,II
   REAL(SP), PARAMETER  :: KAPPA = .40_SP   !!VON KARMAN LENGTH SCALE
   REAL(SP), PARAMETER  :: VK2   = .160_SP  !!KAPPA SQUARED
   REAL(SP)             :: CBCMIN,Z0,ZTEMP,BTPS,RR,U_TAUB,Z0B_GOTM,Z0B_TEMP

#  if defined (TWO_D_MODEL)
!  FOR 2D MODEL ONLY
   REAL(SP), PARAMETER  :: CONST_CD=.0015_SP      !! CD SET CONSTANT TO THIS VALUE
   REAL(SP), PARAMETER  :: ALFA = .166667_SP,   & !! POWER OF WATER DEPTH
                           NN   = 0.02_SP         !! FACTOR TO DIVIDE
   REAL(SP), PARAMETER  :: CFMIN   = .0025_SP,  & !! DEEP WATER VALUE
                           H_BREAK = 1.0_SP,    & !! 
                           THETA   = 10._SP,    & !!
                           LAMB    = 0.3333333333_SP
#  endif
!==============================================================================!
#  if !defined (TWO_D_MODEL)
!
!  SET CONSTANTS
!
   CBCMIN = BFRIC
   Z0     = Z0B

!==============================================================================|
   IF(BROUGH_TYPE == 'orig')THEN !USE ORIGINAL FVCOM FORM FOR BOTTOM FRICTION  |
!==============================================================================|
     DO I=1,MElem
       IF(DT1(I) > 3.0)THEN
        ZTEMP=(ZZ1(I,KBM1)-Z1(I,KB))*DT1(I)/Z0
        CBC(I) = MAX(CBCMIN,VK2/(LOG(ZTEMP))**2)
       ELSE
        ZTEMP=(ZZ1(I,KBM1)-Z1(I,KB))*3.0/Z0
        CBC(I) = MAX(CBCMIN,VK2/(LOG(ZTEMP))**2)
       END IF
     END DO
!==============================================================================|
   ELSE IF(BROUGH_TYPE == 'gotm')THEN !GOTM FORMULATION FOR BOTTOM FRICTION    |
!==============================================================================|

!----Convert Input Z0B to GOTMS H0B
!     H0B = Z0B/.03  
! DAS fixed bug to match gotm's friction.f90
     DO I=1,MElem
     U_TAUB = 0.
     DO II=1,40       
        IF (UMOL <= 0.) THEN
!           Z0B_GOTM= 0.03*Z0B      
           Z0B_GOTM=Z0B   !0.03*H0B 
        ELSE
           Z0B_GOTM=0.1*UMOL/MAX(UMOL,U_TAUB)+Z0B !0.03*H0B
        END IF
          ztemp=(zz1(I,kbm1)-z1(I,kb))*dt1(i)
        RR=KAPPA/(LOG((Z0B_GOTM+ZTEMP)/Z0B_GOTM))
      U_TAUB = RR*SQRT( U(I,KBM1)*U(I,KBM1) + V(I,KBM1)*V(I,KBM1) )
     END DO
     CBC(I) =   RR*RR
     END DO

!==============================================================================|
   ELSE IF(BROUGH_TYPE == 'user_defined')THEN !Use User Defined broud_ud.F     | 
!==============================================================================|
   
     CALL BROUGH_UD

   END IF
#  else
   DO I=1,MElem
!  1.CONSTANT CD
!    CBC(I) = CONST_CD
!  2.  formula 1
     IF(DT1(I) < 3.0_SP)THEN
       CBC(I) = 0.0027_SP
     ELSE
       CBC(I) = GRAV*(DT1(I)**ALFA/NN)**(-2)
     END IF
!  3. formula 2
!    CBC(I) = CFMIN*(1+(H_BREAK/DT1(I))**THETA)**(LAMB/THETA)
!    IF(MOD(I,1000) == 0) print*,CBC(I),DT1(I),CFMIN,H_BREAK,LAMB/THETA
         
   END DO
#  endif

!==============================================================================|
!  CALCULATE SHEAR STRESS ON BOTTOM  --> WUBOT/WVBOT                           |
!==============================================================================|
!   open(333, file = 'strs.out')
!   write(333,*) thour
   DO  I = 1, MElem
     IF(D1(I) > 0.0_SP) THEN
#      if !defined (TWO_D_MODEL)     
       BTPS= CBC(I)*SQRT(U(I,KBM1)**2+V(I,KBM1)**2)
       WUBOT(I) = -BTPS * U(I,KBM1)
       WVBOT(I) = -BTPS * V(I,KBM1)
#      else
       BTPS= CBC(I)*SQRT(UA(I)**2+VA(I)**2)
       WUBOT(I) = -BTPS * UA(I)
       WVBOT(I) = -BTPS * VA(I)
#      endif       
     ELSE
       WUBOT(I) = 0.0_SP
       WVBOT(I) = 0.0_SP
     END IF
!   write(333,*) wubot(i), wvbot(i)
   
   END DO
!   close(333)   !!T.W. test for shear stress output

   RETURN
   END SUBROUTINE BOTTOM_ROUGHNESS
!==============================================================================|
