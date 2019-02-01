!==============================================================================|
!     ARCHIVE THE MODEL RESULTS                                                |
!==============================================================================|

   SUBROUTINE ARCHIVE

!==============================================================================|
#  if defined (NETCDF_IO)
   USE MOD_NCDIO
   USE MOD_NCDAVE
#  endif
#  if defined (WET_DRY)
   USE MOD_WD
#  endif   
#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
#  endif
#  if defined (SEDIMENT)
   USE MOD_SED
#  endif
   USE CONTROL  
#  if defined (RRK_PRE)
   USE MOD_RRK
#  endif
#  if defined (ICE)
   USE MOD_ICE
#  endif

   IMPLICIT NONE
   CHARACTER(LEN=8)   :: CH8
   CHARACTER(LEN=100) :: COPYFILE
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=4)   :: CH4
!==============================================================================|

! comment by T.W., copied from ZY's sediment model source code
!
!--DUMP MEDM FILE EVERY IRECORD  ITERATIONS (NO DUMP IF IRECORD  = 0)----------!
!
   IF(IRECORD /= 0)THEN
     IF(MOD(IINT,IRECORD) == 0)THEN
       CALL OUT_BINARY(IINT/IRECORD)

#  if defined(SEDIMENT)
       IF(SEDIMENT_ON) CALL Out_Binary_Sed(IINT/IRECORD)
#  endif

     END IF
   END IF
!
!----------------------DUMP SMS FILE-------------------------------------------!
!
   IF(IDMPSMS /= 0)THEN
     IF(MOD(IINT,IDMPSMS) == 0)THEN
       IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    SMS FILE'    
       CALL OUT_SMS_ONE
     END IF
   END IF
!
!--DUMP RESTART FILE EVERY IRESTART ITERATIONS (NO DUMP IF IRESTART = 0)-------!
!

#  if !defined (RRK_PRE)
   IF(IRESTART /= 0)THEN
     IF(MOD(IINT,IRESTART) == 0)THEN
       IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    RESTART FILE'
       CALL ARCRST

#      if defined(WET_DRY)
       IF(WET_DRY_ON) THEN
          WRITE(CH8,'(I8.8)') IINT
	  FNAME = 're_'//trim(CH8)//'_wd'
          CALL WD_DUMP(FNAME)
       ENDIF
#      endif       
#      if defined(SEDIMENT)
       IF(SEDIMENT_ON) CALL ARCHIVE_SED
#      endif       

#      if defined(ICE)
       CALL ARCHIVE_ICE
#      endif

     END IF
   END IF
#  else 
   IF(IRESTART /= 0 .OR. REF_TIME1 /= 0)THEN
     IF(IRESTART /= 0 .AND. (MOD(IINT,IRESTART) == 0 .OR. IINT==REF_TIME1))THEN

       IF(MSR) WRITE(IPT,*)  '!  DUMPING               :    RESTART FILE'
       CALL ARCRST

       IF(MSR) THEN
         IF(IINT==REF_TIME1) THEN
           WRITE(CH8,'(I8.8)') IINT
#          if !defined (DOS)
           COPYFILE = "cp -f "//"re_"//trim(CH8)//" "//TRIM(OUTDIR)//"/rrktemp"  
	   IF(MSR) THEN
#            if !defined (CRAY)
	     CALL SYSTEM( TRIM(COPYFILE) )  
#            endif
#            if defined (CRAY)
             CALL CRAY_SYSTEM_CALL(TRIM(COPYFILE)) 
#            endif             
           ENDIF
#          endif     
         ENDIF
       ENDIF

#      if defined(WET_DRY)
       IF(WET_DRY_ON) THEN
         write(CH8,'(I8.8)') IINT
         FNAME = 're_'//trim(CH8)//'_wd'	   
         CALL WD_DUMP(FNAME)

         IF(MSR) THEN
           IF(IINT==REF_TIME1) THEN
             WRITE(CH8,'(I8.8)') IINT
#            if !defined (DOS)
             COPYFILE = "cp -f "//"re_"//trim(CH8)//"_wd "//TRIM(OUTDIR)//"/rrktemp"  
             IF(MSR) THEN
#              if !defined (CRAY)
               CALL SYSTEM( TRIM(COPYFILE) )  
#              endif
#              if defined (CRAY)
               CALL CRAY_SYSTEM_CALL(TRIM(COPYFILE)) 
#              endif             
             ENDIF
#            endif     
           ENDIF
         ENDIF
       ENDIF
	     
#      endif       

     END IF
   END IF
#  endif

!
!--FLOW FIELD AVERAGING AND OUTPUT OF FLOWFIELD AVERAGES-----------------------!
!
   IF(AVGE_ON) CALL OUT_AVGE
#  if defined (NETCDF_IO)
   IF(CDF_OUT_AVE) CALL AVGE_FIELDS(IINT) 
#  endif

!
!--NETCDF OUTPUT---------------------------------------------------------------!
!
#  if defined (NETCDF_IO)
   IF(CDF_INT /= 0 .AND. CDF_OUT)THEN
     IF(MOD(IINT,CDF_INT)==0) CALL OUT_NETCDF
   END IF
#  endif

#  if defined (BALANCE_2D)
   IF(OUT_BALANCE)THEN
     CALL OUT_TIMESERIES_BALANCE
   ENDIF
#  endif

   RETURN
   END SUBROUTINE ARCHIVE
!==============================================================================|
