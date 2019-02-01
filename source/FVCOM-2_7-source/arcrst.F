!==============================================================================|
!   DUMP DATA FILE FOR RESTART                                                 |
!==============================================================================|

   SUBROUTINE ARCRST            

!------------------------------------------------------------------------------|

   USE ALL_VARS

#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif   
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
#  if defined (ENKF_ASSIM)
   USE MOD_ENKF
#  endif

   IMPLICIT NONE
   REAL(SP), ALLOCATABLE :: RTP(:)
   INTEGER I,K,ME,NPC,N1
   CHARACTER(LEN=6) :: FILENUMBER
   CHARACTER(LEN=120) :: DIR,FNAME
   CHARACTER(LEN=8) :: CH8
!==============================================================================|
   
   ME = MYID ; NPC = NPROCS 
   ALLOCATE(RTP(MElem)) ; RTP = 0.0_SP

   IF(MSR)THEN
#  if defined (ENKF_ASSIM)
     WRITE(FILENUMBER,'(I6.6)') IENS
     DIR = TRIM(OUTDIR)//"/fct/"
     FNAME = 'restart'//FILENUMBER//'.dat'
     OPEN(1,FILE=TRIM(DIR)//TRIM(FNAME),FORM='UNFORMATTED',STATUS='UNKNOWN')
#  elif defined (ETKF_ASSIM)
     WRITE(FILENUMBER,'(I6.6)') IENS
     DIR = TRIM(OUTDIR)//"/fct/"
     FNAME = 'restart'//FILENUMBER//'.dat'
     OPEN(1,FILE=TRIM(DIR)//TRIM(FNAME),FORM='UNFORMATTED',STATUS='UNKNOWN')
#  elif defined (RRK_PRE)        
     write(CH8,'(I8.8)') IINT
     OPEN(1,FILE='re_'//trim(CH8),FORM='UNFORMATTED',STATUS='UNKNOWN')        
# else

     if (IRESTART .NE. 0) then
        WRITE(FILENUMBER,'(I6.6)') IINT/IRESTART
        OPEN(1,FILE='re_'//FILENUMBER//'.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')      
     else
        WRITE(FILENUMBER,'(I6.6)') IINT
        OPEN(1,FILE='re_'//FILENUMBER//'.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')      
     end if

#  endif
     REWIND(1)
     WRITE(1) IINT
   END IF

   IF(SERIAL)THEN
     WRITE(1) ((U(I,K),    K=1,KB),I=0,MElem)
     WRITE(1) ((V(I,K),    K=1,KB),I=0,MElem)
     WRITE(1) ((W(I,K),    K=1,KB),I=0,MElem)
#    if defined (GOTM)
     WRITE(1) ((TKE(I,K),   K=1,KB),I=0,NNode)
     WRITE(1) ((TEPS(I,K),  K=1,KB),I=0,NNode)
#    else
     WRITE(1) ((Q2(I,K),   K=1,KB),I=0,NNode)
     WRITE(1) ((Q2L(I,K),  K=1,KB),I=0,NNode)
     WRITE(1) ((L(I,K)  ,  K=1,KB),I=0,NNode)
#    endif
     WRITE(1) ((S(I,K),    K=1,KB),I=0,MElem)
     WRITE(1) ((T(I,K),    K=1,KB),I=0,MElem)
     WRITE(1) ((RHO(I,K),  K=1,KB),I=0,MElem)
     WRITE(1) ((TMEAN(I,K),K=1,KB),I=0,MElem)
     WRITE(1) ((SMEAN(I,K),K=1,KB),I=0,MElem)
     WRITE(1) ((RMEAN(I,K),K=1,KB),I=0,MElem)

     WRITE(1) ((S1(I,K),    K=1,KB),I=1,NNode)
     WRITE(1) ((T1(I,K),    K=1,KB),I=1,NNode)
     WRITE(1) ((RHO1(I,K),  K=1,KB),I=1,NNode)
     WRITE(1) ((TMEAN1(I,K),K=1,KB),I=1,NNode)
     WRITE(1) ((SMEAN1(I,K),K=1,KB),I=1,NNode)
     WRITE(1) ((RMEAN1(I,K),K=1,KB),I=1,NNode)

     WRITE(1) ((KM(I,K),K=1,KB),I=1,NNode)
     WRITE(1) ((KH(I,K),K=1,KB),I=1,NNode)
     WRITE(1) ((KQ(I,K),K=1,KB),I=1,NNode)

     WRITE(1) (UA(I), I=0,MElem)
     WRITE(1) (VA(I), I=0,MElem)

     WRITE(1) (EL1(I), I=1,MElem)
     WRITE(1) (ET1(I), I=1,MElem)
     WRITE(1) (H1(I),  I=1,MElem)
     WRITE(1) (D1(I),  I=1,MElem)
     WRITE(1) (DT1(I), I=1,MElem)
     WRITE(1) (RTP(I), I=1,MElem)

     WRITE(1) (EL(I), I=1,NNode)
     WRITE(1) (ET(I), I=1,NNode)
     WRITE(1) (H(I),  I=1,NNode)
     WRITE(1) (D(I),  I=1,NNode)
     WRITE(1) (DT(I), I=1,NNode)

#    if defined (EQUI_TIDE)
     WRITE(1) (EL_EQI(I), I=1,NNode)
#    endif
#    if defined (ATMO_TIDE)
     WRITE(1) (EL_ATMO(I), I=1,NNode)
#    endif

#    if defined (WATER_QUALITY)
     DO N1 = 1, NB
       WRITE(1) ((WQM(I,K,N1),K=1,KB),I=1,NNode)
     END DO
#    endif

#    if defined (DYE_RELEASE)
     IF(IINT.GT.IINT_SPE_DYE_B) THEN
     WRITE(1) ((DYE(I,K),K=1,KB),I=1,NNode)
     WRITE(1) ((DYEMEAN(I,K),K=1,KB),I=1,NNode)
     ENDIF
#    endif


   ELSE
#     if defined (MULTIPROCESSOR)
      CALL PWRITE(1,ME,NPC,U,    LBOUND(U,1),    UBOUND(U,1),    MElem,MElemGL,KB,EMAP,0,"U"    )
      CALL PWRITE(1,ME,NPC,V,    LBOUND(V,1),    UBOUND(V,1),    MElem,MElemGL,KB,EMAP,0,"V"    )
      CALL PWRITE(1,ME,NPC,W,    LBOUND(W,1),    UBOUND(W,1),    MElem,MElemGL,KB,EMAP,0,"W"    )
#     if defined (GOTM)
      CALL PWRITE(1,ME,NPC,TKE,  LBOUND(TKE,1),  UBOUND(TKE,1),  NNode,NNodeGL,KB,NMAP,0,"TKE"  )
      CALL PWRITE(1,ME,NPC,TEPS, LBOUND(TEPS,1), UBOUND(TEPS,1), NNode,NNodeGL,KB,NMAP,0,"TEPS" )
#     else
      CALL PWRITE(1,ME,NPC,Q2,   LBOUND(Q2,1),   UBOUND(Q2,1),   NNode,NNodeGL,KB,NMAP,0,"Q2"   )
      CALL PWRITE(1,ME,NPC,Q2L,  LBOUND(Q2L,1),  UBOUND(Q2L,1),  NNode,NNodeGL,KB,NMAP,0,"Q2L"  )
      CALL PWRITE(1,ME,NPC,L,    LBOUND(L,1),    UBOUND(L,1),    NNode,NNodeGL,KB,NMAP,0,"L"  )
#     endif
      CALL PWRITE(1,ME,NPC,S,    LBOUND(S,1),    UBOUND(S,1),    MElem,MElemGL,KB,EMAP,0,"S"    )
      CALL PWRITE(1,ME,NPC,T,    LBOUND(T,1),    UBOUND(T,1),    MElem,MElemGL,KB,EMAP,0,"T"    )
      CALL PWRITE(1,ME,NPC,RHO,  LBOUND(RHO,1),  UBOUND(RHO,1),  MElem,MElemGL,KB,EMAP,0,"RHO"  )
      CALL PWRITE(1,ME,NPC,TMEAN,LBOUND(TMEAN,1),UBOUND(TMEAN,1),MElem,MElemGL,KB,EMAP,0,"TMEAN")
      CALL PWRITE(1,ME,NPC,SMEAN,LBOUND(SMEAN,1),UBOUND(SMEAN,1),MElem,MElemGL,KB,EMAP,0,"SMEAN")
      CALL PWRITE(1,ME,NPC,RMEAN,LBOUND(RMEAN,1),UBOUND(RMEAN,1),MElem,MElemGL,KB,EMAP,0,"RMEAN")

      CALL PWRITE(1,ME,NPC,S1,    LBOUND(S1,1),    UBOUND(S1,1),    NNode,NNodeGL,KB,NMAP,1,"S1"    )
      CALL PWRITE(1,ME,NPC,T1,    LBOUND(T1,1),    UBOUND(T1,1),    NNode,NNodeGL,KB,NMAP,1,"T1"    )
      CALL PWRITE(1,ME,NPC,RHO1,  LBOUND(RHO1,1),  UBOUND(RHO1,1),  NNode,NNodeGL,KB,NMAP,1,"RHO1"  )
      CALL PWRITE(1,ME,NPC,TMEAN1,LBOUND(TMEAN1,1),UBOUND(TMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"TMEAN1")
      CALL PWRITE(1,ME,NPC,SMEAN1,LBOUND(SMEAN1,1),UBOUND(SMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"SMEAN1")
      CALL PWRITE(1,ME,NPC,RMEAN1,LBOUND(RMEAN1,1),UBOUND(RMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"RMEAN1")

      CALL PWRITE(1,ME,NPC,KM,LBOUND(KM,1),UBOUND(KM,1),NNode,NNodeGL,KB,NMAP,1,"KM")
      CALL PWRITE(1,ME,NPC,KH,LBOUND(KH,1),UBOUND(KH,1),NNode,NNodeGL,KB,NMAP,1,"KH")
      CALL PWRITE(1,ME,NPC,KQ,LBOUND(KQ,1),UBOUND(KQ,1),NNode,NNodeGL,KB,NMAP,1,"KQ")

      CALL PWRITE(1,ME,NPC,UA,LBOUND(UA,1),UBOUND(UA,1),MElem,MElemGL,1,EMAP,0,"UA")
      CALL PWRITE(1,ME,NPC,VA,LBOUND(VA,1),UBOUND(VA,1),MElem,MElemGL,1,EMAP,0,"VA")

      CALL PWRITE(1,ME,NPC,EL1,LBOUND(EL1,1),UBOUND(EL1,1),MElem,MElemGL,1,EMAP,1,"EL1")
      CALL PWRITE(1,ME,NPC,ET1,LBOUND(ET1,1),UBOUND(ET1,1),MElem,MElemGL,1,EMAP,1,"ET1")
      CALL PWRITE(1,ME,NPC,H1, LBOUND(H1,1), UBOUND(H1,1), MElem,MElemGL,1,EMAP,1,"H1" )
      CALL PWRITE(1,ME,NPC,D1, LBOUND(D1,1), UBOUND(D1,1), MElem,MElemGL,1,EMAP,1,"D1" )
      CALL PWRITE(1,ME,NPC,DT1,LBOUND(DT1,1),UBOUND(DT1,1),MElem,MElemGL,1,EMAP,1,"DT1")
      CALL PWRITE(1,ME,NPC,RTP,LBOUND(RTP,1),UBOUND(RTP,1),MElem,MElemGL,1,EMAP,1,"RTP")

      CALL PWRITE(1,ME,NPC,EL,LBOUND(EL,1),UBOUND(EL,1),NNode,NNodeGL,1,NMAP,1,"EL")
      CALL PWRITE(1,ME,NPC,ET,LBOUND(ET,1),UBOUND(ET,1),NNode,NNodeGL,1,NMAP,1,"ET")
      CALL PWRITE(1,ME,NPC,H, LBOUND(H,1), UBOUND(H,1), NNode,NNodeGL,1,NMAP,1,"H" )
      CALL PWRITE(1,ME,NPC,D, LBOUND(D,1), UBOUND(D,1), NNode,NNodeGL,1,NMAP,1,"D" )
      CALL PWRITE(1,ME,NPC,DT,LBOUND(DT,1),UBOUND(DT,1),NNode,NNodeGL,1,NMAP,1,"DT")

#     if defined (EQUI_TIDE)
      CALL PWRITE(1,ME,NPC,EL_EQI,LBOUND(EL_EQI,1),UBOUND(EL_EQI,1),NNode,NNodeGL,1,NMAP,1,"EL_EQI")
#     endif
#     if defined (ATMO_TIDE)
      CALL PWRITE(1,ME,NPC,EL_ATMO,LBOUND(EL_ATMO,1),UBOUND(EL_ATMO,1),NNode,NNodeGL,1,NMAP,1,"EL_ATMO")
#     endif

#     if defined (WATER_QUALITY)
      DO N1=1,NB
        CALL PWRITE(1,ME,NPC,WQM(1:NNode,1:KB,N1),LBOUND(WQM(1:NNode,1:KB,N1),1),     &
                      UBOUND(WQM(1:NNode,1:KB,N1),1),NNode,NNodeGL,KB,NMAP,1,"WQM")
      END DO
#     endif
#     if defined (DYE_RELEASE)
     IF(IINT.GT.IINT_SPE_DYE_B) THEN
      CALL PWRITE(1,ME,NPC,DYE,    LBOUND(DYE,1),    UBOUND(DYE,1),    NNode,NNodeGL,KB,NMAP,1,"DYE"    )
      CALL PWRITE(1,ME,NPC,DYEMEAN,LBOUND(DYEMEAN,1),UBOUND(DYEMEAN,1),NNode,NNodeGL,KB,NMAP,1,"DYEMEAN"    )
     ENDIF
#     endif

#     endif
   END IF
   IF(MSR) CLOSE(1)
   DEALLOCATE(RTP)

   RETURN
   END SUBROUTINE ARCRST
!==============================================================================|



