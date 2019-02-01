!==============================================================================|
!  Module: cell_block by T.W. for CELL/EDGE/layer blockage in FVCOM V2.7.1     |
!  Corresponding changes commented by TW in other source files.                |
!  This module is specially written for bridge or other underwater structures  |
!  that may block water flow. (e.g., Hood Canal Floating Bridge)               |
!  Last updated on April, 2011                                                 |
!==============================================================================|

#  if defined (V_BLOCK)
   MODULE VELOCITY_BLOCK
   
   USE MOD_PREC
   USE CONTROL

#  if defined (MULTIPROCESSOR)
   USE MOD_PAR  
#  endif

   IMPLICIT NONE
   SAVE

             LOGICAL     :: BLOCK              !defined if the module is turned on
   INTEGER,  ALLOCATABLE :: CELL_BLOCK(:,:)    !by TW, to save global array of blocked cells
   real(sp), ALLOCATABLE, DIMENSION(:) :: dep_block(:)  !sigma layer

   INTEGER,  ALLOCATABLE :: kount(:)           !to count total # of layer blockage
   REAL(SP), ALLOCATABLE :: bwt(:)             !by TW, to save vertical layer block weight for external mode
!   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CELL_BLOCK_TEMP    !by TW

   CONTAINS

! this subroutine will be if cell_block is defined in makefile   
   SUBROUTINE ALLOC_VAR_BLOCK
   USE ALL_VARS
!   use mod_inp
   IMPLICIT NONE
   
   ALLOCATE(CELL_BLOCK(0:MTElem,KBM1),kount(0:MTElem),bwt(0:MTElem))    !allocate local variables

   allocate(dep_block(0:NTNode))   !sigma

   CELL_BLOCK = 0     !default - no block
   kount = 0          !default - no cell/layer block
   bwt = 0.0          !default - no cell/layer block

   dep_block = 0.0
   
   RETURN
   END SUBROUTINE ALLOC_VAR_BLOCK


   SUBROUTINE SET_BLOCK_PARAM
   USE ALL_VARS
   USE MOD_PREC
   USE CONTROL
   use mod_inp

   IMPLICIT NONE
   INTEGER ::  INTVEC(150),ISCAN,KTEMP
   INTEGER ::  INBLK, inblk2   !by TW for cell velocity blockage
   INTEGER :: I,J,K, ITEMP
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=80)  :: ISTR
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CELL_BLOCK_TEMP
   real, ALLOCATABLE, DIMENSION(:) :: dep_block_temp    !by TW

   FNAME = "./"//trim(casename)//"_run.dat"

   ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)

!------------------------------------------------------------------------------|
!     "BLOCK"   !!CONTROLS CELL VELOCITY BLOCKAGE CALCULATION by TW
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BLOCK",LVAL = BLOCK)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BLOCKAGE SWITCH: ',ISCAN
     CALL PSTOP
   END IF

! output run file...
 IF(MSR) THEN
   !  by TW
   IF(BLOCK)THEN
     WRITE(IPT,*)'!  # VELOCITY BLOCKAGE :  ACTIVE'
   ELSE
     WRITE(IPT,*)'!  # VELOCITY BLOCKAGE :  INACTIVE'
   END IF
!  end by TW
 END IF

 IF(BLOCK) THEN
     ALLOCATE(CELL_BLOCK_TEMP(0:MElemGL,KBM1))   

!sigma     
     allocate(dep_block_temp(0:NNodeGL))
     
     inblk = 333
     CALL FOPEN(INBLK, TRIM(ISTR)//'_block_cell.dat', "cfr")       

! for sigma layer
     inblk2 = 334
     CALL FOPEN(inblk2, TRIM(ISTR)//'_block_node.dat', "cfr")
! end sigma

       READ(INBLK,*)     !Skip headline
       DO I=1, MElemGL
          READ(INBLK,*) ITEMP, (CELL_BLOCK_TEMP(I,J), J=1,KBM1)
       END DO
     close(inblk)

   !sigma
    read(inblk2,*)
    do i=1, NNodeGL
       read(inblk2,*) itemp, dep_block_temp(i)
    end do

    close(inblk2)

 END IF
   
   IF(BLOCK) THEN 
       IF(SERIAL) THEN
         CELL_BLOCK = CELL_BLOCK_TEMP

   !sigma
         dep_block = dep_block_temp

!----update vertical layer#
          kount = 0
          bwt = 0.0
          do i=1,MElemGL
             do j=1,kbm1
                if(CELL_BLOCK(I,J)==1) kount(i) = kount(i) + 1      !count the blocked vertical layer#
                bwt(i) = bwt(i) + float(CELL_BLOCK(I,J))*DZ1(i,j)
!                write(999,*) bwt(i)
!                if(bwt(i).ge.1.0) stop 'Serial!'
             end do
          end do
!----end update
       END IF

#  if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,MElem
         CELL_BLOCK(I,:) = CELL_BLOCK_TEMP(EGID(I),:)
       END DO

!sigma
       do i=1,NNode
          dep_block(i) = dep_block_temp(ngid(i))
       end do

!----update vertical layer#
          kount = 0
          bwt = 0.0
          do i=1,MElem
             do j=1,kbm1
                if(CELL_BLOCK(I,J)==1) kount(i) = kount(i) + 1      !count the blocked vertical layer#
                bwt(i) = bwt(i) + float(CELL_BLOCK(I,J))*DZ1(i,j)
       !         if(bwt(i).ge.1.0) stop 'Parellel!'
             end do
          end do
!----end update

       DO I=1,MHE
         CELL_BLOCK(I+MElem,:) = CELL_BLOCK_TEMP(HE_LST(I),:)
       END DO

!sigma
       do i=1,nhn
         dep_block(i+NNode) = dep_block_temp(hn_lst(i))
       end do

!----update vertical layer#
          do i=1,MHE
             do j=1,kbm1
                if(CELL_BLOCK(I+MElem,J)==1) kount(i+MElem) = kount(i+MElem) + 1      !count the blocked vertical layer#
                bwt(i+MElem) = bwt(i+MElem) + float(CELL_BLOCK(I+MElem,J))*DZ1(i+MElem,j)
     !        if(bwt(i+MElem).ge.1.0) stop 'Parellel!'
             end do
          end do
!----end update
     END IF
# endif

  DEALLOCATE(CELL_BLOCK_TEMP)
  deallocate(dep_block_temp)

  END IF

  RETURN

   END SUBROUTINE SET_BLOCK_PARAM
!==============================================================================|
  END MODULE VELOCITY_BLOCK
# endif
