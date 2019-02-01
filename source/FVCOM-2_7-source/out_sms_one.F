!==============================================================================|
!   Write Output Files Used by SMS Post-Processing Software                    |
!==============================================================================|

   SUBROUTINE OUT_SMS_ONE         

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,T1TMP,S1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: UATMP,VATMP,UUWINDTMP,VVWINDTMP,ELTMP 
   INTEGER :: I1 

!==============================================================================|
   
!------------------------------------------------------------------------------!
!  OPEN AND REWIND FILES                                                       !
!------------------------------------------------------------------------------!

   IF(MSR)THEN
     OPEN(IOSMSV,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//'_uvi_uva.xy',STATUS='unknown')
     OPEN(IOSMST,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//'_elts.xy',STATUS='unknown')
     REWIND(IOSMSV)
     REWIND(IOSMST)
   END IF

!------------------------------------------------------------------------------!
!  WRITE TO FILES (SERIAL EXECUTION)                                           !
!------------------------------------------------------------------------------!

!  WRITE NODAL SURFACE GRID COORDINATES, ELEVATION, SALINITY, AND TEMPERATURE-
   IF(SERIAL)THEN
     WRITE(IOSMST,10)
     WRITE(IOSMST,30) NNode,3
     DO I1=1,NNode
       WRITE(IOSMST,'(5E17.8)')VX(I1)+VXMIN,VY(I1)+VYMIN,EL(I1),T1(I1,1),S1(I1,1)
     END DO
     WRITE(IOSMST,*) 'IINT==',IINT
     WRITE(IOSMST,*) 'T1EL==',TIME*86400.0_SP
     WRITE(IOSMST,*) 'T2EL==',T2EL+DELTT
     
 
!  WRITE ELEMENT SURFACE GRID, VELOCITY, VERT-AVGED VELOCITY, AND WIND DATA
   WRITE(IOSMSV,10)
   WRITE(IOSMSV,20) MElem,6
   DO I1=1,MElem
     WRITE(IOSMSV,'(8E17.8)') &
       XC(I1)+VXMIN,YC(I1)+VYMIN,U(I1,1),V(I1,1),UA(I1),VA(I1),UUWIND(I1),VVWIND(I1)
   END DO
   WRITE(IOSMSV,*) 'IINT==',IINT
   END IF

!------------------------------------------------------------------------------!
!  WRITE TO FILES (MULTIPROCESSOR)                                             !
!------------------------------------------------------------------------------!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN

!  WRITE NODAL SURFACE GRID COORDINATES, ELEVATION, SALINITY, AND TEMPERATURE-
   ALLOCATE(ELTMP(NNodeGL))
   ALLOCATE(T1TMP(NNodeGL,KB))
   ALLOCATE(S1TMP(NNodeGL,KB))
   CALL GATHER(LBOUND(EL,1), UBOUND(EL,1), NNode,NNodeGL,1 ,MYID,NPROCS,NMAP,EL, ELTMP)
   CALL GATHER(LBOUND(T1,1), UBOUND(T1,1), NNode,NNodeGL,KB,MYID,NPROCS,NMAP,T1, T1TMP)
   CALL GATHER(LBOUND(S1,1), UBOUND(S1,1), NNode,NNodeGL,KB,MYID,NPROCS,NMAP,S1, S1TMP)
   IF(MSR)THEN
     WRITE(IOSMST,10)
     WRITE(IOSMST,30) NNodeGL,3
     DO I1=1,NNodeGL
       WRITE(IOSMST,'(5E17.8)')XG(I1)+VXMIN,YG(I1)+VYMIN,ELTMP(I1),T1TMP(I1,1),S1TMP(I1,1)
     END DO
     WRITE(IOSMST,*) 'IINT==',IINT
     WRITE(IOSMST,*) 'T1EL==',TIME*86400.0_SP
     WRITE(IOSMST,*) 'T2EL==',T2EL+DELTT
   END IF
   DEALLOCATE(ELTMP,T1TMP,S1TMP)

!  WRITE ELEMENT SURFACE GRID, VELOCITY, VERT-AVGED VELOCITY, AND WIND DATA
   ALLOCATE(UTMP(MElemGL,KB))
   ALLOCATE(VTMP(MElemGL,KB))
   ALLOCATE(UATMP(MElemGL))
   ALLOCATE(VATMP(MElemGL))
   ALLOCATE(UUWINDTMP(MElemGL))
   ALLOCATE(VVWINDTMP(MElemGL))
   CALL GATHER(LBOUND(U,1),     UBOUND(U,1),   MElem,MElemGL,KB,MYID,NPROCS,EMAP,U, UTMP) 
   CALL GATHER(LBOUND(V,1),     UBOUND(V,1),   MElem,MElemGL,KB,MYID,NPROCS,EMAP,V, VTMP) 
   CALL GATHER(LBOUND(UA,1),    UBOUND(UA,1),  MElem,MElemGL,1 ,MYID,NPROCS,EMAP,UA,UATMP) 
   CALL GATHER(LBOUND(VA,1),    UBOUND(VA,1),  MElem,MElemGL,1 ,MYID,NPROCS,EMAP,VA,VATMP) 
   CALL GATHER(LBOUND(UUWIND,1),UBOUND(UUWIND,1),MElem,MElemGL,1,MYID,NPROCS,EMAP,UUWIND,UUWINDTMP) 
   CALL GATHER(LBOUND(VVWIND,1),UBOUND(VVWIND,1),MElem,MElemGL,1,MYID,NPROCS,EMAP,VVWIND,VVWINDTMP) 
   IF(MSR)THEN
       WRITE(IOSMSV,10)
     WRITE(IOSMSV,20) MElemGL,6
     DO I1=1,MElemGL
       WRITE(IOSMSV,'(8E17.8)') &
         XCG(I1)+VXMIN,YCG(I1)+VYMIN,UTMP(I1,1),VTMP(I1,1), &
         UATMP(I1),VATMP(I1),UUWINDTMP(I1),VVWINDTMP(I1)
     END DO
     WRITE(IOSMSV,*) 'IINT==',IINT
   END IF
   DEALLOCATE(UTMP,VTMP,UATMP,VATMP,UUWINDTMP,VVWINDTMP)

   END IF !!PAR
#  endif


!------------------------------------------------------------------------------!
!  CLOSE UP FILES                                                              !
!------------------------------------------------------------------------------!

   IF(MSR)CLOSE(IOSMSV)
   IF(MSR)CLOSE(IOSMST)


!
!--FORMATS---------------------------------------------------------------------!
!

10 FORMAT('scat2d')
20 FORMAT('xyd ',I10,' suva ',I2,' su sv ua va uwind vwind')
30 FORMAT('xyd ',I10,' elts ',I2,' el t s ')

   RETURN
   END SUBROUTINE OUT_SMS_ONE
!==============================================================================|




