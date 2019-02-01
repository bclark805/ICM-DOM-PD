!==============================================================================|
!   Write Binary Data For Iteration Number "IINTT"                             |
!==============================================================================|

   SUBROUTINE OUT_BINARY(IINTT)   

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
#  if defined (NG_OI_ASSIM)
   USE MOD_ASSIM
#  endif
#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif   
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IINTT
   INTEGER :: I,K
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,WWTMP,KMTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: UATMP,VATMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: T1TMP,S1TMP,R1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP

#  if defined (DYE_RELEASE)
   REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: DYETMP
#  endif

   CHARACTER(LEN=4) :: FILENUMBER
   CHARACTER(LEN=120) :: DIR 
   CHARACTER(LEN=120) :: FNAME,FNAME1,FNAME2

#  if defined (WET_DRY)
#  if defined (MULTIPROCESSOR)
   INTEGER, ALLOCATABLE, DIMENSION(:)  ::  ISWETNTMP,ISWETCTMP
#  endif
#  endif

!==============================================================================|
   

!------------------------------------------------------------------------------!
!  OPEN FILE (Name Based on Iteration Number)                                  !
!------------------------------------------------------------------------------!

   IF(MSR)THEN
     WRITE(FILENUMBER,'(I4.4)') IINTT
     FNAME1 = TRIM(CASENAME)//'_assim'//FILENUMBER//'.dat'
     FNAME2 = TRIM(CASENAME)//'_sim'//FILENUMBER//'.dat'
     DIR = TRIM(OUTDIR)//"/medm"
#    if defined (NG_OI_ASSIM)
     IF(SST_ASSIM)       FNAME = FNAME1 
     IF(.NOT. SST_ASSIM) FNAME = FNAME2 
#    else
     FNAME = FNAME2 
#    endif
     OPEN(1,FILE=TRIM(DIR)//"/"//TRIM(FNAME),STATUS='unknown',FORM='unformatted') 
     REWIND(1)
     WRITE(IPT,*)'DUMPING MEDM FILE: ',TRIM(FNAME)
   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Single Processor Case)                                !
!------------------------------------------------------------------------------!

   IF(SERIAL)THEN

     !! ELEMENT BASED VALUES
     WRITE(1) IINTT,MElemGL,NNodeGL,THOUR
#    if !defined (TWO_D_MODEL)
     DO I=1,MElemGL
       WRITE(1) (U(I,K),V(I,K),WW(I,K),KM1(I,K),K=1,KBM1)
     END DO
     
     !! NODE BASED VALUES
#       if !defined (DYE_RELEASE)     
     DO I=1,NNode
       WRITE(1) EL(I),(T1(I,K),S1(I,K),RHO1(I,K),K=1,KBM1)
     END DO
#       else
           DO I=1,NNode
              WRITE(1) EL(I),(T1(I,K),S1(I,K),RHO1(I,K),DYE(I,K),K=1,KBM1)
           END DO
#       endif   
#    else
     DO I=1,MElemGL
       WRITE(1) UA(I),VA(I)
     END DO
     
     !! NODE BASED VALUES
     DO I=1,NNode
       WRITE(1) EL(I)
     END DO
#    endif
   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Multi Processor Case)                                 !
!------------------------------------------------------------------------------!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN

#    if !defined (TWO_D_MODEL)
     !!GATHER AND WRITE ELEMENT-BASED QUANTITIES (U,V,WW,KH)
     ALLOCATE(UTMP(MElemGL,KB))
     ALLOCATE(VTMP(MElemGL,KB))
     ALLOCATE(WWTMP(MElemGL,KB))
     ALLOCATE(KMTMP(MElemGL,KB))
     CALL GATHER(LBOUND(U,1),  UBOUND(U,1),  MElem,MElemGL,KB,MYID,NPROCS,EMAP,U,  UTMP)
     CALL GATHER(LBOUND(V,1),  UBOUND(V,1),  MElem,MElemGL,KB,MYID,NPROCS,EMAP,V,  VTMP)
     CALL GATHER(LBOUND(WW,1), UBOUND(WW,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP,WW, WWTMP)
     CALL GATHER(LBOUND(KM1,1), UBOUND(KM1,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP,KM1, KMTMP)
     IF(MSR)THEN
       WRITE(1) IINTT,MElemGL,NNodeGL,THOUR
       DO I=1,MElemGL
         WRITE(1) (UTMP(I,K),VTMP(I,K),WWTMP(I,K),KMTMP(I,K),K=1,KBM1)
       END DO
     END IF
     DEALLOCATE(UTMP,VTMP,WWTMP,KMTMP)

     !!GATHER AND WRITE NODE-BASED QUANTITIES (EL,T1,S1,RHO1)
     ALLOCATE(ELTMP(NNodeGL))
     ALLOCATE(T1TMP(NNodeGL,KB))
     ALLOCATE(S1TMP(NNodeGL,KB))
     ALLOCATE(R1TMP(NNodeGL,KB))
     CALL GATHER(LBOUND(EL,1),  UBOUND(EL,1),  NNode,NNodeGL, 1,MYID,NPROCS,NMAP,EL,  ELTMP)
     CALL GATHER(LBOUND(T1,1),  UBOUND(T1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,T1,  T1TMP)
     CALL GATHER(LBOUND(S1,1),  UBOUND(S1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,S1,  S1TMP)
     CALL GATHER(LBOUND(RHO1,1),UBOUND(RHO1,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,RHO1,R1TMP)

#    if defined (DYE_RELEASE)
     ALLOCATE(DYETMP(NNodeGL,KB))
     CALL GATHER(LBOUND(DYE,1),UBOUND(DYE,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,DYE,DYETMP)
#    endif     
     IF(MSR)THEN
#    if !defined (DYE_RELEASE)
       DO I=1,NNodeGL
       WRITE(1) ELTMP(I),(T1TMP(I,K),S1TMP(I,K),R1TMP(I,K),K=1,KBM1)
       END DO
#    else
       DO I=1,NNodeGL
         WRITE(1) ELTMP(I),(T1TMP(I,K),S1TMP(I,K),R1TMP(I,K),&
	          DYETMP(I,K),K=1,KBM1)
       END DO
#     endif	 
     END IF
     DEALLOCATE(ELTMP,T1TMP,S1TMP,R1TMP)
#    else
     !!GATHER AND WRITE ELEMENT-BASED QUANTITIES (UA,VA)
     ALLOCATE(UATMP(MElemGL))
     ALLOCATE(VATMP(MElemGL))
     CALL GATHER(LBOUND(UA,1),UBOUND(UA,1),MElem,MElemGL,1,MYID,NPROCS,EMAP,UA,UATMP)
     CALL GATHER(LBOUND(VA,1),UBOUND(VA,1),MElem,MElemGL,1,MYID,NPROCS,EMAP,VA,VATMP)
     IF(MSR)THEN
       WRITE(1) IINTT,MElemGL,NNodeGL,THOUR
       DO I=1,MElemGL
         WRITE(1) UATMP(I),VATMP(I)
       END DO
     END IF
     DEALLOCATE(UATMP,VATMP)

     !!GATHER AND WRITE NODE-BASED QUANTITIES (EL)
     ALLOCATE(ELTMP(NNodeGL))
     CALL GATHER(LBOUND(EL,1),  UBOUND(EL,1),  NNode,NNodeGL, 1,MYID,NPROCS,NMAP,EL,  ELTMP)
     IF(MSR)THEN
       DO I=1,NNodeGL
       WRITE(1) ELTMP(I)
       END DO
     END IF
     DEALLOCATE(ELTMP)
#    endif       
   END IF
#  endif

#  if defined (WET_DRY)
!
!---WRITE NODE VALUES FOR WET_DRY CASE (ISWET_NODE_CURRENTSTEP,ISWET_CELL_CURRENTSTEP)------------------------!
!
   IF(SERIAL)THEN
       WRITE(1) (ISWET_NODE_CURRENTSTEP(I),I=1,NNode)
       WRITE(1) (ISWET_CELL_CURRENTSTEP(I),I=1,MElem)
#    if defined (MULTIPROCESSOR)
   ELSE

!    GATHER QUANTITIES
     ALLOCATE(ISWETNTMP(NNodeGL))
     ALLOCATE(ISWETCTMP(MElemGL))
     CALL IGATHER(LBOUND(ISWET_NODE_CURRENTSTEP,1),UBOUND(ISWET_NODE_CURRENTSTEP,1),NNode,NNodeGL,1,MYID,NPROCS,NMAP, &
                 ISWET_NODE_CURRENTSTEP,ISWETNTMP)
     CALL IGATHER(LBOUND(ISWET_CELL_CURRENTSTEP,1),UBOUND(ISWET_CELL_CURRENTSTEP,1),MElem,MElemGL,1,MYID,NPROCS,EMAP, &
                 ISWET_CELL_CURRENTSTEP,ISWETCTMP)
     IF(MSR)THEN
       WRITE(1) (ISWETNTMP(I),I=1,NNodeGL)
       WRITE(1) (ISWETCTMP(I),I=1,MElemGL)
     END IF
     DEALLOCATE(ISWETNTMP,ISWETCTMP)
#    endif
   END IF
#  endif

   IF(MSR) CLOSE(1)

   RETURN
   END SUBROUTINE OUT_BINARY
!==============================================================================|
