!------------------------------------------------------------------------------------!
!----This is the Kelp module for Uranium extraction project, T.W. April 2013---------!

!---- This module always uses double precision for calculations ----!

!--- The effect of kelp on flow is modeled as the friction force in the momentum
!--- equations, with additional modification to the turbulence scheme.

!--- Corresponding changes in other source code is commented by T.W.------!
!
!--------Last modified by T.W. in June, 2013-------------!


MODULE MOD_KELP
# if defined (KELP)
 !  USE MOD_TYPES
   USE MOD_PREC
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   IMPLICIT NONE
   SAVE

   INTEGER :: KELP_MAP
   LOGICAL :: C_KELP

   CONTAINS !-----------------------------------------------------------------------|
            ! DATA_RUN_KELP : Input Parameters Which Control KELP MODULE            |
            ! IOFILES_KELP  : Open Input Files for KELP MODULE                      |
	    
!==============================================================================|
!   Input Parameters                                                           |
!==============================================================================|

   SUBROUTINE DATA_RUN_KELP            
!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_INP
   IMPLICIT NONE
   INTEGER  INTVEC(150), ISCAN
   CHARACTER(LEN=120) :: FNAME


!==============================================================================|
!   READ IN VARIABLES AND SET VALUES                                           |
!==============================================================================|

   FNAME = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------|
! Kelp module on/off flag, scan the run file to see if it is activated or not
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"C_KELP",LVAL = C_KELP )
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING C_KELP: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "C_KELP" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF

!-- output diagnostic infor. to run file...
 IF(MSR) THEN
   !  by T.W.
   IF(C_KELP)THEN
     WRITE(IPT,*)'!  # KELP MODULE :  ACTIVE'
   ELSE
     WRITE(IPT,*)'!  # KELP MODULE :  INACTIVE'
   END IF
!  end by T.W.
 END IF

   
   RETURN
   END SUBROUTINE DATA_RUN_KELP    
!------------------------------------------------------------------------------|


!==============================================================================!
!   Open & Read Input Files for KELP Calculation                                !
!==============================================================================!

   SUBROUTINE IOFILES_KELP
   USE ALL_VARS
   USE MOD_UTILS
   USE CONTROL
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif

   IMPLICIT NONE

   INTEGER :: eid,kid, I,J,K
   CHARACTER(LEN=80)  :: ISTR

   INTEGER, ALLOCATABLE :: n_of_kelp_gl(:,:)

   double precision, ALLOCATABLE :: drag_coef_gl(:,:), kelp_area_gl(:,:)

   KELP_MAP = 101

   ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)
!
!-----------------OPEN KELP_MAP FILE (SPATIAL VARYING Kelp bed MAP)-----------!

   CALL FOPEN(KELP_MAP, TRIM(ISTR)//'_kelp_map.dat', "cfr")
   
 
   open(102,file='check_input.out')
   
   ALLOCATE(n_of_kelp_gl(0:MElemGL,kbm1),kelp_area_gl(0:MElemGL,kbm1),drag_coef_gl(0:MElemGL,kbm1))

!initialized arrays to ZERO
   n_of_kelp_gl       = 0
   kelp_area_gl       = 0.
   drag_coef_gl       = 0.

   rewind(KELP_MAP)
   READ(KELP_MAP,*)  !skip the title line

   DO I=1, MElemGL       !read in element-based inputs
      do K=1, KBM1
         READ(KELP_MAP,*) eid, kid, n_of_kelp_gl(i,k), kelp_area_gl(i,k),&
                          drag_coef_gl(i,k)

   !      write(102,*)    'Reading ',EID  !eid, kid, n_of_kelp_gl(i,k), kelp_area_gl(i,k),&
         !                 drag_coef_gl(i,k)
      end do
write(102,*)'reading ',EID
   END DO
   
CLOSE(KELP_MAP)


!!transfer to local domains from global (3 variables read from _kelp_map.dat only.)
!FOR SERIAL MODE FIRST...
   IF(SERIAL) THEN

     n_of_kelp(1:MElem,:)    = n_of_kelp_gl(1:MElem,:)
     kelp_area(1:MElem,:)    = kelp_area_gl(1:MElem,:)
     drag_coef(1:MElem,:)    = drag_coef_gl(1:MElem,:)

   END IF

!FOR PARALLEL MODE
#  if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,MElem

         n_of_kelp(i,:)  = n_of_kelp_gl(EGID(i),:)
         kelp_area(i,:)  = kelp_area_gl(EGID(i),:)
         drag_coef(i,:)  = drag_coef_gl(EGID(i),:)

write(102,*)'broadcasting ',I
       END DO
       
       DO I=1,MHE

         n_of_kelp(i+MElem,:)  = n_of_kelp_gl(HE_LST(I),:)
         kelp_area(i+MElem,:)  = kelp_area_gl(HE_LST(I),:)
         drag_coef(i+MElem,:)  = drag_coef_gl(HE_LST(I),:)

       END DO
     END IF
#  endif

   DEALLOCATE(n_of_kelp_gl, kelp_area_gl, drag_coef_gl)

   RETURN
   END SUBROUTINE IOFILES_KELP
!==============================================================================!


!==============================================================================|
!   CALCULATES:  momentum sink by kelp                                         !
!==============================================================================|

   SUBROUTINE MS_KELP
!==============================================================================!
   USE ALL_VARS
   USE CONTROL

   IMPLICIT NONE
   INTEGER :: i,k
   double precision :: rtemp

!calculate element based momentum sink by kelp, EMS(1:N,kbm1)
! momentum_sink = 0.5 * Ce * A * u^2

!re-initialize array to make sure they are cleared to ZERO

   EMS_X = 0.
   EMS_Y = 0.

!start calculation   
   DO i=1,MElem

      do k=1,kbm1
    
	
	! rtemp =    (-1.0)* float(n_of_kelp(i,k))*0.5*drag_coef(i,k)*kelp_area(i,k)*Z1(i,k)*D*SQRT(U(i,k)**2+V(i,k)**2)
	  ! Previous version used a static depth, this will update the depth with each time step, it is more realistic, especially for
	  ! intertidal areas  ! B Clark, Oct 2017
	  !uses the face centered values for Z1 (Sigma level thickness at element) ad D1 (dynamic depth at element)
		  
   rtemp =   float(n_of_kelp(i,k))*0.5*drag_coef(i,k)*kelp_area(i,k)*DZ1(i,k)*D1(i)*SQRT(U(i,k)**2+V(i,k)**2) 
  if(rtemp < 0.) THEN
	write(*,*)'Rtemp  ='
   write(*,*)Rtemp
   endif

   EMS_X(i,k) = rtemp*U(i,k)
   EMS_Y(i,k) = rtemp*V(i,k)


      end do
   
   END DO
   

   RETURN
   END SUBROUTINE MS_KELP
!==============================================================================|

# endif

END MODULE MOD_KELP


