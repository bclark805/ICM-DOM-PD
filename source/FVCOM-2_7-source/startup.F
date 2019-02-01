!==============================================================================|
!   Begin Restart Run From Specified Time                                      |
!==============================================================================|

   SUBROUTINE STARTUP             

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (WET_DRY)
   USE MOD_WD
#  endif   
#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif   
#  if defined (RRK_PRE)
   USE MOD_RRK 
#  endif
#  if defined (ENKF_ASSIM)
   USE MOD_ENKF
#  endif
   USE BCS

   IMPLICIT NONE

   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=8)   :: RRKINP1
   CHARACTER(LEN=4)   :: RRKINP2
   CHARACTER(LEN=4)   :: ENKINP
!==============================================================================|
   



#  if defined (WET_DRY)
   IF(WET_DRY_ON) CALL SET_WD_DATA
#  endif
!
!--Set up Temperature, Salinity, and Turbulence Quantity Fields----------------!
! 
      
   IF((RESTART == 'cold_start').AND.(S_TYPE == 'non-julian'))THEN

     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    COLD_START'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    NON-JULIAN'
     CALL INITIAL_TS
#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL SET_WD_DATA
#    endif
     CALL INITIAL_QQL
#    if defined (WATER_QUALITY)
     CALL INITIAL_WQM
#    endif

#    if defined (DYE_RELEASE)
     CALL INITIAL_DYE
#    endif

   ELSE IF((RESTART=='cold_start').AND.(S_TYPE=='julian'))THEN

     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    COLD_START'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    JULIAN'
     CALL INITIAL_TS
     CALL INITIAL_UVEL
#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL SET_WD_DATA
#    endif
     CALL INITIAL_QQL
#    if defined (WATER_QUALITY)
     CALL INITIAL_WQM
#    endif

#    if defined (DYE_RELEASE)
     CALL INITIAL_DYE
#    endif

   ELSE IF((RESTART=='hot_cold_s').AND.(S_TYPE=='julian'))THEN
          
     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    HOT_COLD_S'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    JULIAN'
     CALL HOT_START_DATA
     CALL INITIAL_TS
#    if defined (WATER_QUALITY)
     CALL INITIAL_WQM
#    endif
#    if defined (WET_DRY)
#    if defined(RRK_PRE)      
     IF(RRK_RUNCONTR == 2)THEN
       WRITE(RRKINP1,'(I8.8)') REF_TIME1
       FNAME=TRIM(OUTDIR)//'/rrktemp/'//'re_'//RRKINP1//'_wd'        
       IF(WET_DRY_ON) CALL WD_READ(FNAME) 
     ELSE IF(RRK_RUNCONTR > 2)THEN
       WRITE(RRKINP2,'(I4.4)') RRK_RUNCONTR-2
       FNAME=TRIM(OUTDIR)//'/rrktemp/'//'eofini'//RRKINP2//'_wd'
       IF(WET_DRY_ON) CALL WD_READ(FNAME)
     ENDIF            
#    elif defined(ENKF_ASSIM)
     IF(ICYC >= 1)THEN
       WRITE(ENKINP,'(I4.4)') IENS
       FNAME=TRIM(OUTDIR)//'/anl/restart'//enkinp//'_wd.dat'
       IF(WET_DRY_ON) CALL WD_READ(FNAME)
     ELSE
       FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
       IF(WET_DRY_ON) CALL WD_READ(FNAME)	    
     ENDIF 
#    else
     FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
     IF(WET_DRY_ON) CALL WD_READ(FNAME)    
#    endif
#    endif

#    if defined (DYE_RELEASE)
     CALL INITIAL_DYE
#    endif

   ELSE IF(RESTART == 'hot_start') THEN
     
     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    HOT_START'
     CALL HOT_START_DATA
     IF(MSR)WRITE(IPT,*)  '!  RESTART DATA          :    READ     '
#    if defined (WET_DRY)
#    if defined(RRK_PRE)      
     IF(RRK_RUNCONTR == 2)THEN
       WRITE(RRKINP1,'(I8.8)') REF_TIME1
       FNAME=TRIM(OUTDIR)//'/rrktemp/'//'re_'//RRKINP1//'_wd'        
       IF(WET_DRY_ON) CALL WD_READ(FNAME) 
     ELSE IF(RRK_RUNCONTR > 2)THEN
       WRITE(RRKINP2,'(I4.4)') RRK_RUNCONTR-2
       FNAME=TRIM(OUTDIR)//'/rrktemp/'//'eofini'//RRKINP2//'_wd'
       IF(WET_DRY_ON) CALL WD_READ(FNAME)
     ENDIF            
#    elif defined(ENKF_ASSIM)
     IF(ICYC > 0)THEN
       WRITE(ENKINP,'(I4.4)') iens
       FNAME=TRIM(OUTDIR)//'/anl/restart'//enkinp//'_wd.dat'
       IF(WET_DRY_ON) CALL WD_READ(FNAME)
     ELSE
       FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
       IF(WET_DRY_ON) CALL WD_READ(FNAME)	    
     ENDIF 
#    else
     FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
     IF(WET_DRY_ON) CALL WD_READ(FNAME)    
#    endif
#    endif

#    if defined (DYE_RELEASE)
     CALL INITIAL_DYE
#    endif

   ELSE
         
     PRINT*,'RESTAR AND S_TYPE DEFINITION NOT CORRECT'
     PRINT*,'RESTAR==',RESTART
     PRINT*,'S_TYPE==',S_TYPE
     CALL PSTOP
         
   END IF
   
   !
!--Set Water Depth-Using Bathymetry and Free Surface Elevation-----------------!
!

   CALL WATER_DEPTH

!
!--Set Values in the Halos-----------------------------------------------------!
! 

   IF(SERIAL)RETURN
# if defined (MULTIPROCESSOR)
   CALL EXCHANGE_ALL
# endif

   RETURN
   END SUBROUTINE STARTUP
!==============================================================================|


!==============================================================================|
!   Exchange All Flow Variables                                                |
!==============================================================================|

   SUBROUTINE EXCHANGE_ALL 

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
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
#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif   

   IMPLICIT NONE

!==============================================================================|
#  if defined (MULTIPROCESSOR)
#  if defined (GOTM)
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,TKE,TEPS,L)
#  else
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,Q2,Q2L,L)
#  endif
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,KM,KQ,KH)
   CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,T,S,RHO)
   CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,TMEAN,SMEAN,RMEAN)
   CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,U,V,W)
   CALL EXCHANGE(EC,MTElem,1 ,MYID,NPROCS,UA,VA)
   CALL EXCHANGE(EC,MTElem,1 ,MYID,NPROCS,EL1,D1,H1)
   CALL EXCHANGE(EC,MTElem,1 ,MYID,NPROCS,ET1,DT1)

   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,TMEAN1,SMEAN1,RMEAN1)
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,S1,T1,RHO1)
   CALL EXCHANGE(NC,NTNode,1 ,MYID,NPROCS,EL,D,H)
   CALL EXCHANGE(NC,NTNode,1 ,MYID,NPROCS,ET,DT)

#  if defined (EQUI_TIDE)
   CALL EXCHANGE(NC,NTNode,1 ,MYID,NPROCS,EL_EQI)
#  endif
#  if defined (ATMO_TIDE)
   CALL EXCHANGE(NC,NTNode,1 ,MYID,NPROCS,EL_ATMO)
#  endif

#  if defined (WATER_QUALITY)
   CALL EXCHANGE_WQM  
   IF(WQM_ON)THEN
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,WMEAN(:,:,1),WMEAN(:,:,2),WMEAN(:,:,3))
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,WMEAN(:,:,4),WMEAN(:,:,5),WMEAN(:,:,6))
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,WMEAN(:,:,7),WMEAN(:,:,8))
   END IF
#  endif

#  if defined (DYE_RELEASE)
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,DYEMEAN)
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,DYE)
#  endif   

#  endif

   CALL RHO_MEAN

   RETURN
   END SUBROUTINE EXCHANGE_ALL
!==============================================================================|
