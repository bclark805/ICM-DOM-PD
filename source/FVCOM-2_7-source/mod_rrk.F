MODULE RRKVAL  !Kalman filter
# if defined (RRK_PRE) || defined(RRK_ASSIM)
   
   USE CONTROL
   IMPLICIT NONE
  
   REAL(DP),ALLOCATABLE :: STTEMP0(:)
   REAL(DP),ALLOCATABLE :: STTEMP1(:)
   REAL(DP),ALLOCATABLE :: STEOF(:)
   REAL(DP),ALLOCATABLE :: SDEOF(:)
   REAL(DP),ALLOCATABLE :: TRANS(:,:)
   
   REAL(DP),ALLOCATABLE :: RRKEL(:)
   REAL(DP),ALLOCATABLE :: RRKU(:,:)
   REAL(DP),ALLOCATABLE :: RRKV(:,:)
   REAL(DP),ALLOCATABLE :: RRKT(:,:)
   REAL(DP),ALLOCATABLE :: RRKS(:,:)

# endif
END MODULE RRKVAL

MODULE MOD_RRK
# if defined (RRK_PRE) || defined(RRK_ASSIM)
   
   USE CONTROL
   IMPLICIT NONE
   SAVE
   
   INTEGER  ::  RRK_RUNCONTR = 1
   INTEGER      REF_TIME1          !!GLOBAL NUMBER OF THE START TIME OF FALSE STATE
   INTEGER      REF_TIME2          !!GLOBAL NUMBER OF THE START TIME OF TRUE STATE
   INTEGER      DELTA_ASS 
   INTEGER      REF_INT           !!GLOBAL NUMBER OF THE READING FILE INTERVALS 
   INTEGER      RRK_NVAR
   INTEGER      RRK_NOBSMAX
   INTEGER      RRK_OPTION        !!OPTION 1 FOR BAROTROPIC CASE; OPTION 2 FOR BAROCLINIC CASE
   INTEGER      RRK_NEOF          !!NUMBER OF THE EOF  
   REAL(SP) ::  RRK_PSIZE         !!PERTURBATION SIZE  
   REAL(SP) ::  RRK_PSCALE        !!PSEUDO MODEL ERROR    
   REAL(SP) ::  RRK_RSCALE        !!SCALE FACTOR APPLIED TO ONE STANDARD DEVIATION FOR R
   INTEGER      RRK_START         !!ASSIMILATION START TIME (MUST BE EQUAL OR GREATER THAN REF_TIME)     
   INTEGER      RRK_END           !!ASSIMILATION END TIME (MUST BE EQUAL OR GREATER THAN REF_TIME)

   LOGICAL  ::  EL_ASSIM          !!OPTION FOR CHOSING ELEVATION AS ASSIMILATION VARIABLES
   LOGICAL  ::  UV_ASSIM          !!OPTION FOR CHOSING CURRENT AS ASSIMILATION VARIABLES  
   LOGICAL  ::  T_ASSIM           !!OPTION FOR CHOSING TEMPERATURE AS ASSIMILATION VARIABLES
   LOGICAL  ::  S_ASSIM           !!OPTION FOR CHOSING SALINITY AS ASSIMILATION VARIABLES

   LOGICAL  ::  EL_OBS            !!OPTION FOR ELEVATION OBSERVATION DATA
   LOGICAL  ::  UV_OBS            !!OPTION FOR CURRENT OBSERVATION DATA  
   LOGICAL  ::  T_OBS             !!OPTION FOR TEMPERATURE OBSERVATION DATA
   LOGICAL  ::  S_OBS             !!OPTION FOR SALINITY OBSERVATION DATA

   REAL(SP)  BC_AMP_ERR(6)
   REAL(SP)  BC_PHA_ERR(6)

   INTEGER      INORRK            !!FILE I/O PIPE NUMBER 

   CONTAINS
   
   SUBROUTINE SET_RRK_PARAM 

!------------------------------------------------------------------------------|
!  READ IN PARAMETERS CONTROLLING WET/DRY TREATMENT                            |
!------------------------------------------------------------------------------|

   USE CONTROL
   USE MOD_INP
   IMPLICIT NONE
   INTEGER  :: I,ISCAN, KTEMP
   CHARACTER(LEN=120) :: FNAME
   REAL(SP) REALVEC(150)

!------------------------------------------------------------------------------|
!   READ IN VARIABLES AND SET VALUES                                           |
!------------------------------------------------------------------------------|

   FNAME = TRIM(INPDIR)//"/"//trim(casename)//"_assim_rrkf.dat"

!------------------------------------------------------------------------------|
!   REF_TIME1  - GLOBAL NUMBER OF THE START TIME OF FALSE STATE                 |
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"REF_TIME1",ISCAL = REF_TIME1) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING REF_TIME1: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE REF_TIME NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!   REF_TIME2  - GLOBAL NUMBER OF THE START TIME OF TRUE STATE
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"REF_TIME2",ISCAL = REF_TIME2) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING REF_TIME2: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF

!------------------------------------------------------------------------------|
!   DELTA_ASS  - GLOBAL NUMBER OF THE DATA ASSIMILATION TIME INTERVAL          |
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"DELTA_ASS",ISCAL = DELTA_ASS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DELTA_ASS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE DELTA_ASS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF
!------------------------------------------------------------------------------|
!   REF_INT  - GLOBAL NUMBER OF THE READING FILE INTERVALS 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"REF_INT",ISCAL = REF_INT) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING REF_INT: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE REF_INT NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF
!------------------------------------------------------------------------------|
!   RRK_NVAR  -  
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_NVAR",ISCAL = RRK_NVAR) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_NVAR: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE RRK_NVAR NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF
!------------------------------------------------------------------------------|
!   RRK_NOBSMAX  -  
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_NOBSMAX",ISCAL = RRK_NOBSMAX) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_NOBSMAX: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE RRK_NOBSMAX NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF
!------------------------------------------------------------------------------|
!   RRK_OPTION  -  OPTION 1 FOR BAROTROPIC CASE; OPTION 2 FOR BAROCLINIC CASE
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_OPTION",ISCAL = RRK_OPTION) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_OPTION: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE RRK_OPTION NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF      
!------------------------------------------------------------------------------|
!   RRK_NEOF  -  THE NUMBER OF EOFs
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_NEOF",ISCAL = RRK_NEOF) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_NEOF: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE RRK_NEOF NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF  
!------------------------------------------------------------------------------|
!   RRK_PSIZE  -  PERTURBATION SIZE
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_PSIZE",FSCAL = RRK_PSIZE) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_PSIZE: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE RRK_PSIZE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF  
!------------------------------------------------------------------------------|
!   RRK_PSCALE  -  PSEUDO MODEL ERROR
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_PSCALE",FSCAL = RRK_PSCALE) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_PSCALE: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE PSEUDO NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF  
!------------------------------------------------------------------------------|
!   RRK_RSCALE  -  SCALE FACTOR APPLIED TO ONE STANDARD DEVIATION FOR R
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_RSCALE",FSCAL = RRK_RSCALE) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_RSCALE: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE RRK_RSCALE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF  
!-----------------------------------------------------------------------------------|
!   RRK_START  -  ASSIMILATION START TIME                                           | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_START",ISCAL = RRK_START) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_START: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE RRK_START NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF  

!-----------------------------------------------------------------------------------|
!   RRK_END  -  ASSIMILATION END TIME                                               | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"RRK_END",ISCAL = RRK_END) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING RRK_END: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE RRK_END NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF     

!------------------------------------------------------------------------------|
!     "BC_AMP_ERR"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BC_AMP_ERR",FVEC = REALVEC,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BC_AMP_ERR: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= 6)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED TIDAL COMSTITUENTS ERROR ARE NOT EQUAL TO 6' 
     END IF
   END IF
  
   BC_AMP_ERR(1:6)= REALVEC(1:6)
!------------------------------------------------------------------------------|
!     "BC_PHA_ERR"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BC_PHA_ERR",FVEC = REALVEC,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BC_PHA_ERR: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= 6)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED TIDAL COMSTITUENTS ERROR ARE NOT EQUAL TO 6' 
     END IF
   END IF
  
   BC_PHA_ERR(1:6)= REALVEC(1:6)
!-----------------------------------------------------------------------------------|
!   EL_ASSIM  -  OPTION FOR CHOSING ELEVATION AS ASSIMILATION VARIABLES             | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"EL_ASSIM",LVAL = EL_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING EL_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE EL_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF   
   
!-----------------------------------------------------------------------------------|
!   EL_OBS  -  ELEVATION OBSERVATION DATA OPTION                                    | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"EL_OBS",LVAL = EL_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING EL_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE EL_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF   

!-----------------------------------------------------------------------------------|
!   UV_ASSIM  -  OPTION FOR CHOSING CURRENT AS ASSIMILATION VARIABLES               | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"UV_ASSIM",LVAL = UV_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING UV_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE UV_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!-----------------------------------------------------------------------------------|
!   UV_OBS  -  CURRENT OBSERVATION DATA OPTION                                      | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"UV_OBS",LVAL = UV_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING UV_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE UV_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!-----------------------------------------------------------------------------------|
!   T_ASSIM  -   OPTION FOR CHOSING TEMPERATURE AS ASSIMILATION VARIABLES           | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"T_ASSIM",LVAL = T_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING T_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE T_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!-----------------------------------------------------------------------------------|
!   T_OBS  -  TEMPERATURE OBERVATION DATA OPTION                                    | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"T_OBS",LVAL = T_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING T_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE T_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!-----------------------------------------------------------------------------------|
!   S_ASSIM  -  OPTION FOR CHOSING SALINITY AS ASSIMILATION VARIABLES               | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"S_ASSIM",LVAL = S_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING S_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE S_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF

!-----------------------------------------------------------------------------------|
!   S_OBS  -  SALINITY OBSERVATION DATA OPTION                                      | 
!-----------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"S_OBS",LVAL = S_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING S_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE S_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF

!------------------------------------------------------------------------------|
!            SCREEN REPORT OF SET VARIABlES                                    !
!------------------------------------------------------------------------------|
    IF(MSR) THEN
    WRITE(IPT,*)
    WRITE(IPT,*)'!            KALMAN FILTER PARAMETERS        '  
    WRITE(IPT,*)'!  # KALMAN FILTER FALSE STATE START TIME   :', REF_TIME1
    WRITE(IPT,*)'!  # KALMAN FILTER TRUE STATE START TIME    :', REF_TIME2
    WRITE(IPT,*)'!  # DATA ASSIMILATION TIME STEP            :', DELTA_ASS
    WRITE(IPT,*)'!  # KALMAN FILTER REF FILE READING INTERVAL:', REF_INT
    WRITE(IPT,*)'!  # KALMAN FILTER PARAMETER                :', RRK_NVAR
    WRITE(IPT,*)'!  # KALMAN FILTER PARAMETER                :', RRK_NOBSMAX
    WRITE(IPT,*)'!  # KALMAN FILTER PARAMETER                :', RRK_OPTION
    WRITE(IPT,*)'!  # KALMAN FILTER PARAMETER                :', RRK_NEOF
    WRITE(IPT,*)'!  # KALMAN FILTER PERTURBATION SIZE        :', RRK_PSIZE
    WRITE(IPT,*)'!  # KALMAN FILTER PSEUDO MODEL ERROR       :', RRK_PSCALE
    WRITE(IPT,*)'!  # SCALE FACTOR FOR R                     :', RRK_RSCALE
    WRITE(IPT,*)'!  # ASSIMILATION START TIME                :', RRK_START
    WRITE(IPT,*)'!  # ASSIMILATION END TIME                  :', RRK_END
    WRITE(IPT,*)'!  # TIDAL AMPLITUDE ERROR RANGE SPECIFIED  :',(BC_AMP_ERR(I),I=1,6)
    WRITE(IPT,*)'!  # TIDAL PHASE ERROR RANGE SPECIFIED      :',(BC_PHA_ERR(I),I=1,6)
    WRITE(IPT,*)'!  # ELEVATION AS ASSIMILATION VARIABLES    :', EL_ASSIM
    WRITE(IPT,*)'!  # ELEVATION OBSERVATION DATA OPTION      :', EL_OBS
    WRITE(IPT,*)'!  # CURRENT AS ASSIMILATION VARIABLES      :', UV_ASSIM
    WRITE(IPT,*)'!  # CURRENT OBSERVATION DATA OPTION        :', UV_OBS
    WRITE(IPT,*)'!  # TEMPERATURE AS ASSIMILATION VARIABLES  :', T_ASSIM
    WRITE(IPT,*)'!  # TEMPERATURE OBSERVATION DATA OPTION    :', T_OBS
    WRITE(IPT,*)'!  # SALINITY AS ASSIMILATION VARIABLES     :', S_ASSIM
    WRITE(IPT,*)'!  # SALINITY OBSERVATION DATA OPTION       :', S_OBS
    ENDIF
   
   RETURN

   END SUBROUTINE SET_RRK_PARAM


!------------------------------------------------------------------------------|
!  READ IN *_sim.dat FILES TO CREAT REF.CDF FILES                              |
!------------------------------------------------------------------------------|   
   
   SUBROUTINE RRK_REF(STATE) 

   USE LIMS
   USE CONTROL
   IMPLICIT NONE

#include "/hosts/salmon01/data00/medm/src/netcdf-3.6.0-p1/src/fortran/netcdf.inc"   
!#include "/usr/local/include/nedcdf.inc" 
   CHARACTER(LEN=160) FILENAME
   CHARACTER(LEN=160) OUTNAME
   CHARACTER(LEN=4)   FCYC
   INTEGER FIRST
   INTEGER STATE
   INTEGER IINDEX
   INTEGER IDVQ1
   INTEGER IDVQ2
   INTEGER IDVQ3
   INTEGER IDVQ4
   INTEGER IDVQ5
   INTEGER DIMS(3)
   INTEGER START(3)
   INTEGER COUNT(3)
   INTEGER I,II,K
   INTEGER RCODE
   INTEGER NSTATE
   INTEGER NTIME
   INTEGER TEMP1
   REAL(SP) :: TEMP2
   
   REAL(SP), ALLOCATABLE :: RRKU(:,:,:)
   REAL(SP), ALLOCATABLE :: RRKV(:,:,:)
   REAL(SP), ALLOCATABLE :: RRKTMP1(:,:)
   REAL(SP), ALLOCATABLE :: RRKEL(:,:)
   REAL(SP), ALLOCATABLE :: RRKT(:,:,:)
   REAL(SP), ALLOCATABLE :: RRKS(:,:,:)
   REAL(SP), ALLOCATABLE :: RRKTMP2(:,:)
 
   IF(STATE==1) THEN
      OUTNAME=TRIM(OUTDIR)//'/rrktemp/'//'ref.cdf'
      NTIME= (NSTEPS-REF_TIME1)/DELTA_ASS/REF_INT 
      WRITE(IPT,*) 'Creating ref.cdf file......'
   ELSE
      OUTNAME=TRIM(OUTDIR)//'/rrktemp/'//'ref2.cdf' 
      NTIME= (NSTEPS-REF_TIME2)/DELTA_ASS/REF_INT
      WRITE(IPT,*) 'Creating ref2.cdf file......'
   ENDIF 
   
!  TEMPORARILY ALLOCATE ARRY TO ARRYS  

   ALLOCATE(RRKU(MElemGL,KBM1,NTIME))              ;RRKU    = ZERO
   ALLOCATE(RRKV(MElemGL,KBM1,NTIME))              ;RRKV    = ZERO
   ALLOCATE(RRKTMP1(MElemGL,KBM1))                 ;RRKTMP1 = ZERO 
   ALLOCATE(RRKEL(NNodeGL,NTIME))                  ;RRKEL   = ZERO
   ALLOCATE(RRKT(NNodeGL,KBM1,NTIME))              ;RRKT    = ZERO
   ALLOCATE(RRKS(NNodeGL,KBM1,NTIME))              ;RRKS    = ZERO 
   ALLOCATE(RRKTMP2(NNodeGL,KBM1))                 ;RRKTMP2 = ZERO  

!  END ALLOCATION

   RCODE = nf_create(OUTNAME,NF_CLOBBER,II)     

   IF(EL_ASSIM) THEN
     RCODE = nf_def_dim(II,'n_el',NNodeGL,DIMS(1))
     RCODE = nf_def_dim(II,'t_el',NTIME,DIMS(2))
     RCODE = nf_def_var(II,'el',NF_DOUBLE,2,DIMS,IDVQ1)
   ENDIF
   
   IF(UV_ASSIM) THEN
     RCODE = nf_def_dim(II,'c_u',MElemGL,DIMS(1))
     RCODE = nf_def_dim(II,'z_u',KBM1,DIMS(2))
     RCODE = nf_def_dim(II,'t_u',NTIME,DIMS(3))
     RCODE = nf_def_var(II,'u',NF_DOUBLE,3,DIMS,IDVQ2)
     RCODE = nf_def_var(II,'v',NF_DOUBLE,3,DIMS,IDVQ3)
   ENDIF
   
   IF(T_ASSIM) THEN
     RCODE = nf_def_dim(II,'n_t',NNodeGL,DIMS(1))
     RCODE = nf_def_dim(II,'z_t',KBM1,DIMS(2))
     RCODE = nf_def_dim(II,'t_t',NTIME,DIMS(3))
     RCODE = nf_def_var(II,'temp',NF_DOUBLE,3,DIMS,IDVQ4)
   ENDIF

   IF(S_ASSIM) THEN
     RCODE = nf_def_dim(II,'n_s',NNodeGL,DIMS(1))
     RCODE = nf_def_dim(II,'z_s',KBM1,DIMS(2))
     RCODE = nf_def_dim(II,'t_s',NTIME,DIMS(3))
     RCODE = nf_def_var(II,'sal',NF_DOUBLE,3,DIMS,IDVQ5)
   ENDIF
   
   RCODE = nf_enddef(II)

   IINDEX = 0
   IF(STATE==1) THEN
      FIRST = REF_TIME1/DELTA_ASS+1
   ELSE
      FIRST = REF_TIME2/DELTA_ASS+1
   ENDIF
      DO NSTATE = FIRST, NSTEPS/DELTA_ASS, REF_INT
         IINDEX = IINDEX +1
         WRITE(FCYC,'(I4.4)') NSTATE  
         FILENAME='../medm_bck/'//TRIM(CASENAME)//'_sim'//FCYC//'.dat'
         OPEN(72,FILE=FILENAME,FORM='UNFORMATTED')

         READ(72) TEMP1, TEMP1, TEMP1, TEMP2
	 DO I=1, MElemGL
 	    READ(72) (RRKU(I,K,IINDEX),RRKV(I,K,IINDEX),RRKTMP1(I,K),RRKTMP1(I,K),K=1,KBM1)
!            READ(72) RRKU(I,1,IINDEX),RRKV(I,1,IINDEX) 
	 ENDDO
         
         DO I=1, NNodeGL
	    READ(72) RRKEL(I,IINDEX), (RRKT(I,K,IINDEX),RRKS(I,K,IINDEX),RRKTMP2(I,K),K=1,KBM1)
!            READ(72) RRKEL(I,IINDEX)
	 ENDDO
!        WRITE(IPT,*) NSTATE,RRKEL(1,IINDEX)
         CLOSE(72)
      ENDDO 

      IF(EL_ASSIM) THEN
         START(1) = 1
         START(2) = 1
         COUNT(1) = NNodeGL
         COUNT(2) = NTIME
         RCODE = nf_put_var_double(II,IDVQ1,RRKEL)
      ENDIF
      
      IF(UV_ASSIM) THEN
         START(1) = 1
         START(2) = 1
         START(3) = 1
         COUNT(1) = MElemGL
         COUNT(2) = KBM1
         COUNT(3) = NTIME
         RCODE = nf_put_var_double(II,IDVQ2,RRKU)
         
	 RCODE = nf_put_var_double(II,IDVQ3,RRKV)
      ENDIF
      
      IF(T_ASSIM) THEN
         START(1) = 1
         START(2) = 1
         START(3) = 1
         COUNT(1) = NNodeGL
         COUNT(2) = KBM1
         COUNT(3) = NTIME
         RCODE = nf_put_var_double(II,IDVQ4,RRKT)
      ENDIF
      
      IF(S_ASSIM) THEN
         START(1) = 1
         START(2) = 1
         START(3) = 1
         COUNT(1) = NNodeGL
         COUNT(2) = KBM1
         COUNT(3) = NTIME
         RCODE = nf_put_var_double(II,IDVQ5,RRKS)
      ENDIF

      RCODE = nf_close(II)

   DEALLOCATE(RRKU,RRKV,RRKTMP1,RRKEL,RRKT,RRKS,RRKTMP2)
   END SUBROUTINE RRK_REF

!------------------------------------------------------------------------------|
!  READ IN ref.cdf FILE TO CREAT eof.cdf, eigenvalue.dat etc FILES             |
!------------------------------------------------------------------------------|   
   
   SUBROUTINE RRK_EOF

   USE LIMS
   USE CONTROL
   IMPLICIT NONE

#include "/hosts/salmon01/data00/medm/src/netcdf-3.6.0-p1/src/fortran/netcdf.inc"
!#include "/usr/local/include/nedcdf.inc" 
   INTEGER I,J,K,II

! WORK ARRAYS FOR LAPACK SUBROUTINE
   INTEGER LWORK4, LDVT
   REAL(DP),ALLOCATABLE :: WORK4(:)
   REAL(DP) :: VT
   REAL(DP),ALLOCATABLE :: RKSF(:,:)
   REAL(DP),ALLOCATABLE :: RKSF1(:,:)
   REAL(DP),ALLOCATABLE :: SEOF(:,:)
   REAL(DP),ALLOCATABLE :: SFSF(:,:)
   REAL(DP),ALLOCATABLE :: SFD(:)
   REAL(DP),ALLOCATABLE :: SFU(:,:)
   REAL(DP),ALLOCATABLE :: STVAR(:)
   REAL(DP),ALLOCATABLE :: STSD(:)
   REAL(DP) :: ELSD, USD, TSD, SSD
   REAL(DP) :: SUM0
   CHARACTER(LEN=80) INAME
   CHARACTER(LEN=80) ONAME
   INTEGER IDVQ1,IDVQ2,IDVQ3,IDVQ4,IDVQ5,IDVQ6,IDVQ7,IDVQ8,IDVQ9
   INTEGER IDVQ10,IDVQ11,IDVQ12,IDVQ13,IDVQ14,IDVQ15
   INTEGER NCINP
   INTEGER VINPID
   INTEGER DIMS(3)
   INTEGER START(3)
   INTEGER COUNT(3)
   INTEGER STATUS
   INTEGER RCODE
   INTEGER I_FIRST
   INTEGER I_LAST
   INTEGER I_STEP
   INTEGER NSTEP
   INTEGER IDUMMY,IDUMMY1
   INTEGER I_SELECT
   REAL(DP),ALLOCATABLE :: RRKU(:,:)
   REAL(DP),ALLOCATABLE :: RRKV(:,:)
   REAL(DP),ALLOCATABLE :: RRKEL(:)
   REAL(DP),ALLOCATABLE :: RRKT(:,:)
   REAL(DP),ALLOCATABLE :: RRKS(:,:)
   REAL(DP),ALLOCATABLE :: RRKU2(:,:,:)
   REAL(DP),ALLOCATABLE :: RRKV2(:,:,:)
   REAL(DP),ALLOCATABLE :: RRKEL2(:,:)
   REAL(DP),ALLOCATABLE :: RRKT2(:,:,:)
   REAL(DP),ALLOCATABLE :: RRKS2(:,:,:)
   REAL(DP),ALLOCATABLE :: STMEAN(:)
   INTEGER SS_DIM
   INTEGER STDIM
        
   WRITE(IPT,*) 'Calculate the EOFs from the control run......'
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
    
   SS_DIM = RRK_NVAR*(NSTEPS-REF_TIME1)/DELTA_ASS
   LWORK4 = 5*SS_DIM
   LDVT   = 1 

! TEMPORARILY ALLOCATE ARRY TO ARRYS

   ALLOCATE(WORK4(LWORK4))             ;WORK4   = ZERO
   ALLOCATE(RKSF(STDIM,SS_DIM))        ;RKSF    = ZERO
   ALLOCATE(RKSF1(STDIM,SS_DIM))       ;RKSF1   = ZERO
   ALLOCATE(SEOF(STDIM,SS_DIM))        ;SEOF    = ZERO
   ALLOCATE(SFSF(SS_DIM,SS_DIM))       ;SFSF    = ZERO
   ALLOCATE(SFD(SS_DIM))               ;SFD     = ZERO
   ALLOCATE(SFU(SS_DIM,SS_DIM))        ;SFU     = ZERO
   ALLOCATE(STVAR(STDIM))              ;STVAR   = ZERO
   ALLOCATE(STSD(STDIM))               ;STSD    = ZERO
   ALLOCATE(RRKU(MElemGL,KBM1))            ;RRKU    = ZERO
   ALLOCATE(RRKV(MElemGL,KBM1))            ;RRKV    = ZERO
   ALLOCATE(RRKELNNodeGL)                ;RRKEL   = ZERO
   ALLOCATE(RRKT(NNodeGL,KBM1))            ;RRKT    = ZERO
   ALLOCATE(RRKS(NNodeGL,KBM1))            ;RRKS    = ZERO
   ALLOCATE(RRKU2(MElemGL,KBM1,SS_DIM))    ;RRKU2   = ZERO
   ALLOCATE(RRKV2(MElemGL,KBM1,SS_DIM))    ;RRKV2   = ZERO
   ALLOCATE(RRKEL2(NNodeGL,SS_DIM))        ;RRKEL2  = ZERO
   ALLOCATE(RRKT2(NNodeGL,KBM1,SS_DIM))    ;RRKT2   = ZERO
   ALLOCATE(RRKS2(NNodeGL,KBM1,SS_DIM))    ;RRKS2   = ZERO
   ALLOCATE(STMEAN(STDIM))             ;STMEAN  = ZERO   

! END ALLOCATION   


    INAME=TRIM(OUTDIR)//'/rrktemp/'//'ref.cdf'
    ONAME=TRIM(OUTDIR)//'/rrktemp/'//'eof.cdf'

    I_FIRST  = 1
    I_LAST   = (NSTEPS-REF_TIME1)/DELTA_ASS
    I_STEP   = 1
    I_SELECT = RRK_OPTION
    
    OPEN(72,FILE=TRIM(OUTDIR)//'/rrktemp/'//'eigenvalue.dat')
    OPEN(73,FILE=TRIM(OUTDIR)//'/rrktemp/'//'avgstd.dat')
    OPEN(74,FILE=TRIM(OUTDIR)//'/rrktemp/'//'stmean.dat')
    OPEN(75,FILE=TRIM(OUTDIR)//'/rrktemp/'//'eof1.dat')
    
    STATUS = nf_open(INAME,NF_NOWRITE,NCINP)
    NSTEP  = 0
 
    DO I = I_FIRST, I_LAST, I_STEP
      NSTEP = NSTEP + 1
      IF (STATUS/=0) THEN
        WRITE(IPT,*) 'Could not open cdf file:', INAME
        CALL PSTOP    
      ENDIF
    
! STORE THE ELEVATION IN RKSF
      
      IDUMMY = 0      

      IF(EL_ASSIM) THEN 
        START(1) = 1
        START(2) = I
        COUNT(1) = NNodeGL
        COUNT(2) = 1
        STATUS = nf_inq_varid(NCINP,'el',VINPID)
        STATUS = nf_get_vara_double(NCINP,VINPID,START,COUNT,RRKEL)
      
        DO J = 1, NNodeGL
          IDUMMY = IDUMMY + 1
	  RKSF(IDUMMY,NSTEP) =  RRKEL(J)
        ENDDO 
      ENDIF

      IF(UV_ASSIM) THEN
        START(1) = 1
        START(2) = 1
        START(3) = I
        COUNT(1) = MElemGL
        COUNT(2) = KBM1
        COUNT(3) = 1
        STATUS = nf_inq_varid(NCINP,'u',VINPID)
        STATUS = nf_get_vara_double(NCINP,VINPID,START,COUNT,RRKU)
        DO K = 1, KBM1
          DO J = 1, MElemGL
	     IDUMMY = IDUMMY + 1
	     RKSF(IDUMMY,NSTEP) = RRKU(J,K)
	  ENDDO
        ENDDO	

        STATUS = nf_inq_varid(NCINP,'v',VINPID)
        STATUS = nf_get_vara_double(NCINP,VINPID,START,COUNT,RRKV)
        DO K = 1, KBM1
          DO J = 1, MElemGL
	     IDUMMY = IDUMMY + 1
	     RKSF(IDUMMY,NSTEP) = RRKV(J,K)
	  ENDDO
        ENDDO     
      ENDIF

      IF(T_ASSIM) THEN
        START(1) = 1
	START(2) = 1
        START(3) = I
        COUNT(1) = NNodeGL
        COUNT(2) = KBM1
        COUNT(3) = 1
        STATUS = nf_inq_varid(NCINP,'temp',VINPID)
        STATUS = nf_get_vara_double(NCINP,VINPID,START,COUNT,RRKT)
        DO K = 1, KBM1
          DO J = 1, NNodeGL
	     IDUMMY = IDUMMY + 1
	     RKSF(IDUMMY,NSTEP) = RRKT(J,K)
	  ENDDO
        ENDDO	
      ENDIF

      IF(S_ASSIM) THEN
        START(1) = 1
	START(2) = 1
        START(3) = I
        COUNT(1) = NNodeGL
        COUNT(2) = KBM1
        COUNT(3) = 1
        STATUS = nf_inq_varid(NCINP,'sal',VINPID)
        STATUS = nf_get_vara_double(NCINP,VINPID,START,COUNT,RRKS)
	
	DO K = 1, KBM1
          DO J = 1, NNodeGL
	     IDUMMY = IDUMMY + 1
	     RKSF(IDUMMY,NSTEP) = RRKS(J,K)
	  ENDDO
        ENDDO	
      ENDIF
      
    ENDDO                ! END OF READING THE TIME SERIES FROM THE CONTROL RUN

    STATUS = nf_close(NCINP)
    
! CALCULATE THE MEAN AND THE STANDARD DEVIATION

    DO I=1, STDIM
      SUM0=0.0_DP
      DO J=1, NSTEP
        SUM0=SUM0+RKSF(I,J)
      ENDDO
     
      STMEAN(I) = SUM0/DBLE(NSTEP)
      SUM0=0.0_DP
      DO J=1, NSTEP
        RKSF1(I,J)=(RKSF(I,J)-STMEAN(I))
	SUM0=SUM0+RKSF1(I,J)**2.0_DP
      ENDDO
           
      STVAR(I)=SUM0
      STSD(I)=DSQRT(SUM0/DBLE(NSTEP-1))
      WRITE(74,'(I10,E15.7)') I, STMEAN(I)
    ENDDO    

    IDUMMY  = 0
    IDUMMY1 = 0
    ELSD    = 0.0_DP
    USD     = 0.0_DP
    TSD     = 0.0_DP
    SSD     = 0.0_DP
    
        

    IF(EL_ASSIM) THEN
      SUM0=0.0_DP
      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        SUM0 = SUM0 + STVAR(IDUMMY)
      ENDDO
      
      ELSD=DSQRT(SUM0/DBLE(NNodeGL)/DBLE(NSTEP-1))

      DO I=1, NNodeGL
        DO J=1, NSTEP
          RKSF1(I,J)=RKSF1(I,J)/ELSD/DSQRT(DBLE(NSTEP-1)) ! NORMALIZED THE ELEVATION BY ITS SPATIALLY AVERAGED S.D.
        ENDDO
      ENDDO
    
      IDUMMY1 = IDUMMY
    ENDIF

    IF(UV_ASSIM) THEN
      SUM0=0.0_DP
      DO K=1, 2*KBM1
        DO I=1, MElemGL
           IDUMMY = IDUMMY + 1
	   SUM0 = SUM0 +STVAR(IDUMMY)
        ENDDO
      ENDDO
      USD=DSQRT(SUM0/DBLE(MElemGL*KBM1)/DBLE(NSTEP-1))
 
      DO I=IDUMMY1+1, IDUMMY1+2*MElemGL*KBM1 
        DO J=1, NSTEP
           IF(I_SELECT == 0) THEN
! NORMALIZE THE CURRENT (RRKU, RRKV) BY ITS SPATIALLY AVERAGED S.D. 
	     RKSF1(I,J)=RKSF1(I,J)/USD/DSQRT(DBLE(NSTEP-1))
	   ELSE
! NORMALIZE THE CURRENT (RRKU, RRKV) BY ITS SPATIALLY AVERAGED S.D. AND IT IS ALSO DEVIDED BY THE NUMBER 
! OF LAYERS TO REPRESENT ONLY THE BAROTROPIC RESPONSE	    
	     RKSF1(I,J)=RKSF1(I,J)/USD/DSQRT(DBLE(KBM1))/DSQRT(DBLE(NSTEP-1))
	   ENDIF         
        ENDDO 
      ENDDO
    
      IDUMMY1 = IDUMMY
    ENDIF  
      
    IF(T_ASSIM) THEN
      SUM0=0.0_DP
      DO K=1, KBM1
        DO I=1, NNodeGL
           IDUMMY = IDUMMY + 1
	   SUM0 = SUM0 +STVAR(IDUMMY)
        ENDDO
      ENDDO
      TSD=DSQRT(SUM0/DBLE(NNodeGL*KBM1)/DBLE(NSTEP-1))
 
      DO I=IDUMMY1+1, IDUMMY1+NNodeGL*KBM1 
        DO J=1, NSTEP
! NORMALIZE TEMPERATURE BY ITS SPATIALLY AVERAGED S.D. 
	  RKSF1(I,J)=RKSF1(I,J)/TSD/DSQRT(DBLE(NSTEP-1))
        ENDDO 
      ENDDO

      IDUMMY1 = IDUMMY
    ENDIF        
      
    IF(S_ASSIM) THEN
      SUM0=0.0_DP
      DO K=1, KBM1
        DO I=1, NNodeGL
           IDUMMY = IDUMMY + 1
	   SUM0 = SUM0 +STVAR(IDUMMY)
        ENDDO
      ENDDO
      SSD=DSQRT(SUM0/DBLE(NNodeGL*KBM1)/DBLE(NSTEP-1))
 
      DO I=IDUMMY1+1, IDUMMY1+NNodeGL*KBM1 
        DO J=1, NSTEP
! NORMALIZE TEMPERATURE BY ITS SPATIALLY AVERAGED S.D. 
	  RKSF1(I,J)=RKSF1(I,J)/SSD/DSQRT(DBLE(NSTEP-1))
        ENDDO 
      ENDDO

      IDUMMY1 = IDUMMY
    ENDIF       

!    WRITE(IPT,*) 'STDIM', STDIM, IDUMMY
!    WRITE(IPT,*) 'S.D.', ELSD, USD
    WRITE(73,*) 'SPATIALLY AVERAGED S.D. OF EL, (U,V), TEMP, SAL'
    WRITE(73,'(4E15.7)') ELSD, USD, TSD, SSD

! CALCULATE THE EOFs BY SVD OF RKSF1' RKSF1 = C LAMBDA C', RKSF1 RKSF1' = E LAMBDA E', E = RKSF1 C LAMBDA ^-1/2

    DO I=1, NSTEP
      DO J=1, NSTEP
        SUM0=0.0_DP
	DO K=1, STDIM
	  SUM0 = SUM0 + RKSF1(K,I)*RKSF1(K,J)
	ENDDO
	SFSF(I,J) =SUM0
      ENDDO
    ENDDO
 
    CALL DGESVD('A','N',NSTEP,NSTEP,SFSF,SS_DIM,SFD,SFU,SS_DIM,VT,LDVT,WORK4,LWORK4,RCODE)
    
    DO I=1, NSTEP
!      WRITE(IPT,*) 'LAMBDA', I, SFD(I)
      WRITE(72,'(I5,E15.7)') I, SFD(I)   
    ENDDO
    DO I=1, NSTEP
      DO J=1, NSTEP 
        SFU(I,J)=SFU(I,J)/DSQRT(SFD(J))
      ENDDO
    ENDDO

    DO I=1, STDIM
      DO J=1, NSTEP
         SUM0=0.0_DP 
         DO K=1, NSTEP
	    SUM0=SUM0+RKSF1(I,K)*SFU(K,J)
	 ENDDO
         SEOF(I,J)=SUM0 
      ENDDO
    ENDDO
    
! STORE THE MEAN, ONE STANDARD DEVIATION, ABD EOFs IN NETCDF FORMAT
  
    RCODE = nf_create(ONAME, NF_CLOBBER, II)

    IF(EL_ASSIM) THEN
      RCODE = nf_def_dim(II,'n_mel',NNodeGL,DIMS(1))
      RCODE = nf_def_var(II,'mean_el',NF_DOUBLE,1,DIMS,IDVQ1)
      RCODE = nf_def_var(II,'sd_el',NF_DOUBLE,1,DIMS,IDVQ2)
    ENDIF
    
    IF(UV_ASSIM) THEN
      RCODE = nf_def_dim(II,'c_mu',MElemGL,DIMS(1))
      RCODE = nf_def_dim(II,'z_mu',KBM1,DIMS(2))
      RCODE = nf_def_var(II,'mean_u',NF_DOUBLE,2,DIMS,IDVQ3)
      RCODE = nf_def_var(II,'sd_u',NF_DOUBLE,2,DIMS,IDVQ4)
      RCODE = nf_def_var(II,'mean_v',NF_DOUBLE,2,DIMS,IDVQ5)
      RCODE = nf_def_var(II,'sd_v',NF_DOUBLE,2,DIMS,IDVQ6)
    ENDIF
    
    IF(T_ASSIM) THEN
      RCODE = nf_def_dim(II,'n_mtemp',NNodeGL,DIMS(1))
      RCODE = nf_def_dim(II,'z_mtemp',KBM1,DIMS(2))
      RCODE = nf_def_var(II,'mean_temp',NF_DOUBLE,2,DIMS,IDVQ10)
      RCODE = nf_def_var(II,'sd_temp',NF_DOUBLE,2,DIMS,IDVQ11)
    ENDIF
    
    IF(S_ASSIM) THEN
      RCODE = nf_def_dim(II,'n_msal',NNodeGL,DIMS(1))
      RCODE = nf_def_dim(II,'z_msal',KBM1,DIMS(2))
      RCODE = nf_def_var(II,'mean_sal',NF_DOUBLE,2,DIMS,IDVQ12)
      RCODE = nf_def_var(II,'sd_sal',NF_DOUBLE,2,DIMS,IDVQ13)
    ENDIF

    IF(EL_ASSIM) THEN
      RCODE = nf_def_dim(II,'n_el',NNodeGL,DIMS(1))
      RCODE = nf_def_dim(II,'t_el',NSTEP,DIMS(2))
      RCODE = nf_def_var(II,'eof_el',NF_DOUBLE,2,DIMS,IDVQ7)
    ENDIF
    
    IF(UV_ASSIM) THEN
      RCODE = nf_def_dim(II,'c_u',MElemGL,DIMS(1))
      RCODE = nf_def_dim(II,'z_u',KBM1,DIMS(2))
      RCODE = nf_def_dim(II,'t_u',NSTEP,DIMS(3))
      RCODE = nf_def_var(II,'eof_u',NF_DOUBLE,3,DIMS,IDVQ8)
      RCODE = nf_def_var(II,'eof_v',NF_DOUBLE,3,DIMS,IDVQ9)
    ENDIF
    
    IF(T_ASSIM) THEN
      RCODE = nf_def_dim(II,'n_t',NNodeGL,DIMS(1))
      RCODE = nf_def_dim(II,'z_t',KBM1,DIMS(2))
      RCODE = nf_def_dim(II,'t_t',NSTEP,DIMS(3))
      RCODE = nf_def_var(II,'eof_temp',NF_DOUBLE,3,DIMS,IDVQ14)
    ENDIF
    
    IF(S_ASSIM) THEN
      RCODE = nf_def_dim(II,'n_s',NNodeGL,DIMS(1))
      RCODE = nf_def_dim(II,'z_s',KBM1,DIMS(2))
      RCODE = nf_def_dim(II,'t_s',NSTEP,DIMS(3))
      RCODE = nf_def_var(II,'eof_sal',NF_DOUBLE,3,DIMS,IDVQ15)
    ENDIF    
    
    RCODE = nf_enddef(II)

! WRITE THE MEAN AND ONE STANDARD DEVIATION
    IDUMMY  = 0
    IDUMMY1 = 0

    IF(EL_ASSIM) THEN
      START(1)=1
      START(2)=1
      START(3)=1
      COUNT(1)=NNodeGL
      COUNT(2)=1
      COUNT(3)=1

      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        RRKEL(I) = STMEAN(IDUMMY)
      ENDDO  
!     WRITE(IPT,*) 'mean_el', RRKEL(1)
      RCODE = nf_put_var_double(II,IDVQ1,RRKEL)
    
      IDUMMY = 0
      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        RRKEL(I) = STSD(IDUMMY)
      ENDDO   
      RCODE = nf_put_var_double(II,IDVQ2,RRKEL)  
   
      IDUMMY1 = IDUMMY
    ENDIF    
 
    IF(UV_ASSIM) THEN
      START(1)=1
      START(2)=1
      COUNT(1)=MElemGL
      COUNT(2)=KBM1
      
      IDUMMY = IDUMMY1
      DO K=1, KBM1
        DO I=1, MElemGL 
          IDUMMY = IDUMMY + 1 
          RRKU(I,K) = STMEAN (IDUMMY)
        ENDDO
      ENDDO
!     WRITE(IPT,*) 'mean_u', RRKU(1,1)     
      RCODE = nf_put_var_double(II,IDVQ3,RRKU)
    
      IDUMMY = IDUMMY1
      DO K=1, KBM1
        DO I=1, MElemGL 
          IDUMMY = IDUMMY + 1 
          RRKU(I,K) = STSD(IDUMMY) 
        ENDDO
      ENDDO
      RCODE = nf_put_var_double(II,IDVQ4,RRKU)

      IDUMMY1 = IDUMMY
      IDUMMY  = IDUMMY1
      DO K=1, KBM1
        DO I=1, MElemGL
          IDUMMY = IDUMMY + 1
          RRKV(I,K) = STMEAN(IDUMMY)
        ENDDO
      ENDDO
      RCODE = nf_put_var_double(II,IDVQ5,RRKV)
    
      IDUMMY  = IDUMMY1
      DO K=1, KBM1
        DO I=1, MElemGL
          IDUMMY = IDUMMY + 1
          RRKV(I,K) = STSD(IDUMMY)
        ENDDO
      ENDDO
      RCODE = nf_put_var_double(II,IDVQ6,RRKV)
  
      IDUMMY1 = IDUMMY
    ENDIF

    IF(T_ASSIM) THEN
      START(1)=1
      START(2)=1
      COUNT(1)=NNodeGL
      COUNT(2)=KBM1
      
      IDUMMY = IDUMMY1
      DO K=1, KBM1
        DO I=1, NNodeGL 
          IDUMMY = IDUMMY + 1 
          RRKT(I,K) = STMEAN (IDUMMY)
        ENDDO
      ENDDO
!     WRITE(IPT,*) 'mean_temp', RRKT(1,1)     
      RCODE = nf_put_var_double(II,IDVQ10,RRKT)
    
      IDUMMY = IDUMMY1
      DO K=1, KBM1
        DO I=1, NNodeGL 
          IDUMMY = IDUMMY + 1 
          RRKT(I,K) = STSD(IDUMMY) 
        ENDDO
      ENDDO
      RCODE = nf_put_var_double(II,IDVQ11,RRKT)

      IDUMMY1 = IDUMMY
    ENDIF

    IF(S_ASSIM) THEN
      START(1)=1
      START(2)=1
      COUNT(1)=NNodeGL
      COUNT(2)=KBM1
      
      IDUMMY = IDUMMY1
      DO K=1, KBM1
        DO I=1, NNodeGL 
          IDUMMY = IDUMMY + 1 
          RRKS(I,K) = STMEAN (IDUMMY)
        ENDDO
      ENDDO
!     WRITE(IPT,*) 'mean_sal', RRKS(1,1)     
      RCODE = nf_put_var_double(II,IDVQ12,RRKS)
    
      IDUMMY = IDUMMY1
      DO K=1, KBM1
        DO I=1, NNodeGL 
          IDUMMY = IDUMMY + 1 
          RRKS(I,K) = STSD(IDUMMY) 
        ENDDO
      ENDDO
      RCODE = nf_put_var_double(II,IDVQ13,RRKS)

      IDUMMY1 = IDUMMY
    ENDIF

! WRITE THE EOFs
    IDUMMY  = 0
    IDUMMY1 = 0  
   
    IF(EL_ASSIM) THEN
      DO K=1, NSTEP
        IDUMMY = 0
        DO I=1, NNodeGL
          IDUMMY = IDUMMY + 1
          RRKEL2(I,K) = SEOF(IDUMMY,K)
        ENDDO
        WRITE(75,*) K,RRKEL2(1,K)
      ENDDO

      START(1) = 1
      START(2) = 1
      COUNT(1) = NNodeGL
      COUNT(2) = NSTEP
      RCODE = nf_put_var_double(II,IDVQ7,RRKEL2)
   
      IDUMMY1 = IDUMMY
    ENDIF
    
    IF(UV_ASSIM) THEN
      DO K=1, NSTEP
        IDUMMY = IDUMMY1
        DO J=1, KBM1
          DO I=1, MElemGL
	    IDUMMY = IDUMMY + 1 
            RRKU2(I,J,K) = SEOF(IDUMMY,K)
          ENDDO
        ENDDO    
      ENDDO
     
      START(1) = 1
      START(2) = 1 
      START(3) = 1
      COUNT(1) = MElemGL
      COUNT(2) = KBM1
      COUNT(3) = NSTEP
      RCODE = nf_put_var_double(II,IDVQ8,RRKU2)
    
      IDUMMY1 = IDUMMY

      DO K=1, NSTEP
        IDUMMY = IDUMMY1
        DO J=1, KBM1
          DO I=1, MElemGL
            IDUMMY = IDUMMY + 1
            RRKV2(I,J,K) = SEOF(IDUMMY,K)
          ENDDO    
        ENDDO
      ENDDO  
      RCODE = nf_put_var_double(II,IDVQ9,RRKV2)
    
      IDUMMY1 = IDUMMY    
    ENDIF

    IF(T_ASSIM) THEN
      DO K=1, NSTEP
        IDUMMY = IDUMMY1
        DO J=1, KBM1
          DO I=1, NNodeGL
	    IDUMMY = IDUMMY + 1 
            RRKT2(I,J,K) = SEOF(IDUMMY,K)
          ENDDO
        ENDDO    
      ENDDO
     
      START(1) = 1
      START(2) = 1 
      START(3) = 1
      COUNT(1) = NNodeGL
      COUNT(2) = KBM1
      COUNT(3) = NSTEP
      RCODE = nf_put_var_double(II,IDVQ14,RRKT2)
    
      IDUMMY1 = IDUMMY
    ENDIF

    IF(S_ASSIM) THEN
      DO K=1, NSTEP
        IDUMMY = IDUMMY1
        DO J=1, KBM1
          DO I=1, NNodeGL
	    IDUMMY = IDUMMY + 1 
            RRKS2(I,J,K) = SEOF(IDUMMY,K)
          ENDDO
        ENDDO    
      ENDDO
     
      START(1) = 1
      START(2) = 1 
      START(3) = 1
      COUNT(1) = NNodeGL
      COUNT(2) = KBM1
      COUNT(3) = NSTEP
      RCODE = nf_put_var_double(II,IDVQ15,RRKS2)
    
      IDUMMY1 = IDUMMY
    ENDIF

    RCODE = nf_close(II)  
      
    CLOSE(72)
    CLOSE(73)
    CLOSE(74)
    CLOSE(75)
      
   DEALLOCATE(WORK4,RKSF,RKSF1,SEOF,SFSF,RRKT,RRKS,RRKT2,RRKS2)
   DEALLOCATE(SFD,SFU,STVAR,STSD,RRKU,RRKV,RRKEL,RRKU2,RRKV2,RRKEL2,STMEAN)
   END SUBROUTINE RRK_EOF

!------------------------------------------------------------------------------|
!  Main program to calculate the reduced rank kalman gain matrix               |
!------------------------------------------------------------------------------|   
   
   SUBROUTINE RRK_RRK(CHOICE)

   USE LIMS
   USE CONTROL
   USE RRKVAL
   IMPLICIT NONE
    
    INTEGER SS_DIM
    INTEGER STDIM
    INTEGER I_EOF
    INTEGER I,J,K
    INTEGER CHOICE
    INTEGER NLOC
    INTEGER RCODE
    INTEGER IT
    INTEGER ILAST
    INTEGER IDUMMY
    INTEGER SS_TOT
    INTEGER,ALLOCATABLE  :: STLOC(:)
    REAL(DP) ERR1
    REAL(DP),ALLOCATABLE :: KAL(:,:)
    REAL(DP),ALLOCATABLE :: HEOF(:,:)
    REAL(DP),ALLOCATABLE :: R(:,:)
    REAL(DP),ALLOCATABLE :: HTR(:,:)
    REAL(DP),ALLOCATABLE :: R2(:,:)
    REAL(DP),ALLOCATABLE :: PHT(:,:)
    REAL(DP),ALLOCATABLE :: RALPHA(:,:)
    REAL(DP),ALLOCATABLE :: BETA2(:,:)
    REAL(DP),ALLOCATABLE :: GAMMA(:,:)
    REAL(DP),ALLOCATABLE :: GAMMA2(:,:)
    REAL(DP),ALLOCATABLE :: TRANS2(:,:)
    
! WORK ARRAYS FOR LAPACK SUBROUTINES
    INTEGER LWORK, LWORK2, LDVT
    INTEGER,ALLOCATABLE :: IPIV(:)
    INTEGER,ALLOCATABLE :: IPIV2(:)
    REAL(DP),ALLOCATABLE :: WORK(:)
    REAL(DP),ALLOCATABLE :: WORK2(:)
    REAL(DP) :: VT
    

! CHARACTER STRINGS
    CHARACTER(LEN=80) KALNAME
    CHARACTER(LEN=80) FILENAME
    CHARACTER(LEN=4)  STRCYC
    CHARACTER(LEN=8)  CH8
    CHARACTER(LEN=80) INAME    
    CHARACTER(LEN=80) INAME2
    
    REAL(DP),ALLOCATABLE :: SEOF(:,:)
    REAL(DP),ALLOCATABLE :: STSD(:)
    REAL(DP) :: ELSD, USD, TSD, SSD
    REAL(DP),ALLOCATABLE :: HU(:,:)
    REAL(DP),ALLOCATABLE :: HUL(:,:)
    REAL(DP),ALLOCATABLE :: EVAL(:)
    REAL(DP) :: RRKSUM
      
    STDIM = 0
    IF(EL_ASSIM) STDIM = STDIM + NNodeGL
    IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
    IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
    IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1

    SS_DIM = RRK_NVAR*RRK_NEOF
    LWORK  = 4*SS_DIM
    LWORK2 = 4*RRK_NOBSMAX
    SS_TOT = (NSTEPS-REF_TIME1)/DELTA_ASS
    LDVT   = 1 
    ILAST  = 10  ! THE NUMBER OF ITERATION IN DOUBLING ALGORITHM

! TEMPORARILY ALLOCATE ARRY TO ARRYS

   ALLOCATE(STLOC(RRK_NOBSMAX))            ;STLOC    = 0
   ALLOCATE(IPIV(SS_DIM))                  ;IPIV     = 0
   ALLOCATE(IPIV2(RRK_NOBSMAX))            ;IPIV2    = 0
   ALLOCATE(KAL(SS_DIM,RRK_NOBSMAX))       ;KAL      = ZERO
   ALLOCATE(HEOF(RRK_NOBSMAX,SS_DIM))      ;HEOF     = ZERO
   ALLOCATE(R(RRK_NOBSMAX,RRK_NOBSMAX))    ;R        = ZERO
   ALLOCATE(HTR(SS_DIM,RRK_NOBSMAX))       ;HTR      = ZERO
   ALLOCATE(R2(RRK_NOBSMAX,RRK_NOBSMAX))   ;R2       = ZERO
   ALLOCATE(PHT(SS_DIM,RRK_NOBSMAX))       ;PHT      = ZERO
   ALLOCATE(RALPHA(SS_DIM,SS_DIM))         ;RALPHA   = ZERO
   ALLOCATE(BETA2(SS_DIM,SS_DIM))          ;BETA2    = ZERO
   ALLOCATE(GAMMA(SS_DIM,SS_DIM))          ;GAMMA    = ZERO
   ALLOCATE(GAMMA2(SS_DIM,SS_DIM))         ;GAMMA2   = ZERO
   ALLOCATE(TRANS2(SS_DIM,SS_DIM))         ;TRANS2   = ZERO 
   ALLOCATE(WORK(LWORK))                   ;WORK     = ZERO
   ALLOCATE(WORK2(LWORK2))                 ;WORK2    = ZERO
   ALLOCATE(SEOF(STDIM,SS_DIM))            ;SEOF     = ZERO
   ALLOCATE(STSD(STDIM))                   ;STSD     = ZERO
   ALLOCATE(HU(RRK_NOBSMAX,SS_TOT))        ;HU       = ZERO      
   ALLOCATE(HUL(RRK_NOBSMAX,SS_TOT))       ;HUL      = ZERO
   ALLOCATE(EVAL(SS_TOT))                  ;EVAL     = ZERO

! END ALLOCATION

   IF(CHOICE == 1) THEN
     CALL RRK_ALLOC_VAR
     WRITE(CH8,'(I8.8)') REF_TIME1
     IF(MSR) WRITE(IPT,*) 'Starting perturb the base state in the eigenvector direction......'
     INAME = TRIM(OUTDIR)//'/rrktemp/'//'re_'//TRIM(CH8)  
     CALL PERTURB(INAME)       

     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,SEOF,STSD,HU,HUL,EVAL) 
     DEALLOCATE(STTEMP0,STTEMP1,STEOF,SDEOF,TRANS,RRKEL,RRKU,RRKV,RRKT,RRKS)
     RETURN
   ENDIF

   IF(CHOICE == 2) THEN
     CALL RRK_ALLOC_VAR
     IF(MSR) WRITE(IPT,*) 'Starting calculate the linearized model matrix in the reduced space......'  
     INAME = TRIM(OUTDIR)//'/rrktemp/'//'basfct'
! CALCULATE THE LINEARIZED MODEL MATRIX IN THE REDUCED SPACE
     CALL MREDUCED(INAME)

     IF(MSR) THEN
! OUTPUT LINEAR TRANSITION MATRIX   
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Amatr.dat' 
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         WRITE(INORRK) (TRANS(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)

 ! ALSO IN ASCII
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Amatr1.txt'
       OPEN(INORRK, FILE=TRIM(FILENAME),FORM='FORMATTED')
       DO I=1, SS_DIM
         DO J=1, SS_DIM
           WRITE(INORRK,*) TRANS(J,I)
         ENDDO
       ENDDO
       CLOSE(INORRK)
     ENDIF  
       
     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,SEOF,STSD,HU,HUL,EVAL) 
     DEALLOCATE(STTEMP0,STTEMP1,STEOF,SDEOF,TRANS,RRKEL,RRKU,RRKV,RRKT,RRKS)
     RETURN
   ENDIF

   IF(CHOICE == 3) THEN
     CALL RRK_ALLOC_VAR
     IF(MSR) WRITE(IPT,*) 'Starting perturb the base in the eigenvector directions with the magnitude of eigenvalues......'     
     INAME = TRIM(OUTDIR)//'rrktemp'//'???'     
! PERTURB THE BASE STATE      
     CALL PERTURBL(INAME)

     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,SEOF,STSD,HU,HUL,EVAL) 
     DEALLOCATE(STTEMP0,STTEMP1,STEOF,SDEOF,TRANS,RRKEL,RRKU,RRKV,RRKT,RRKS)
     RETURN
   ENDIF
   
   IF(CHOICE == 4) THEN
     CALL RRK_ALLOC_VAR
     WRITE(IPT,*) 'Calculate the Kalman gain matrix by doubling algorith......'
     KALNAME = TRIM(OUTDIR)//'/rrktemp/rrK.dat'

! READ THE EIGENVALUES AND SET THE INITIAL ERROR COVARIANCE   
     INAME = TRIM(OUTDIR)//'/rrktemp/eigenvalue.dat' 
     OPEN(INORRK,FILE=INAME,STATUS='OLD')
     
     DO I=1, SS_TOT
       READ(INORRK,*) J, EVAL(I)
     ENDDO
     
     DO I=1, SS_DIM
       DO J=1, SS_DIM
         GAMMA2(I,J) = 0._SP
       ENDDO
       GAMMA2(I,I) = EVAL(I)*RRK_PSCALE
     ENDDO
     CLOSE(INORRK)
   
! READ THE LINEAR TRANSITION MATRIX     
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Amatr.dat' 
     OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
     DO J=1, SS_DIM
       READ(INORRK) (TRANS(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)
     
! READ THE OBSERVATION NUMBER AND LOCATION   
     CALL READOBS(STLOC,NLOC)
   
! SET RALPHA   
     DO I=1, SS_DIM
       DO J=1, SS_DIM
         TRANS2(I,J)=TRANS(J,I)
       ENDDO
!      WRITE(IPT,*) 'Amatr:', I, TRANS(I,I)     
     ENDDO
   
! READ ONE STANDARD DEVIATION (D) OF THE CONTROL RUN

     INAME2 = TRIM(OUTDIR)//'/rrktemp/avgstd.dat'   
     OPEN(INORRK, FILE=INAME2, STATUS='OLD')
     
     READ(INORRK,*)
     READ(INORRK,*) ELSD, USD, TSD, SSD
      
     IDUMMY = 0
     IF(EL_ASSIM) THEN  
       DO I=1, NNodeGL
         IDUMMY = IDUMMY + 1
         STSD(IDUMMY) = ELSD
       ENDDO
     ENDIF
     IF(UV_ASSIM) THEN
       DO I=1, 2*KBM1
         DO J=1, MElemGL
           IDUMMY = IDUMMY + 1
           IF(RRK_OPTION == 0) THEN
	     STSD(IDUMMY)= USD
	   ELSE
	     STSD(IDUMMY) = USD*DSQRT(DBLE(KBM1)) 
	   ENDIF
         ENDDO
       ENDDO
     ENDIF
     IF(T_ASSIM) THEN
       DO I=1, KBM1
         DO J=1, NNodeGL
           IDUMMY = IDUMMY + 1
           STSD(IDUMMY)= TSD
         ENDDO
       ENDDO
     ENDIF
     IF(S_ASSIM) THEN
       DO I=1, KBM1
         DO J=1, NNodeGL
           IDUMMY = IDUMMY + 1
           STSD(IDUMMY)= SSD
         ENDDO
       ENDDO
     ENDIF    
     CLOSE(INORRK)

! CALCULATE OBSERVATION MATRIX (EOF TO OBS TRANSFORMATION) H*D*EOFs
     DO J=1, SS_DIM
       CALL READEOF(J,2)
       DO I=1, NLOC
         
         HEOF(I,J)=STTEMP1(STLOC(I))*STSD(STLOC(I))
!        WRITE(IPT,*) 'Heof:', I,J,HEOF(I,J)        
       ENDDO
     ENDDO
    

     FILENAME = TRIM(OUTDIR)//'/rrktemp/ObsOp.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       WRITE(INORRK) (HEOF(I,J), I=1, NLOC)    
     ENDDO 
     CLOSE(INORRK)     
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Heof.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
     DO J=1, SS_DIM
       DO I=1, NLOC
         WRITE(INORRK,*) I,J,HEOF(I,J)
       ENDDO    
     ENDDO 
     CLOSE(INORRK)

! OUTPUT RALPHA MATRIX, WHICH IS TRANSPOSE OF TRANS
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Alpha.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       WRITE(INORRK) (TRANS2(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)     

! OUTPUT GAMMA MATRIX
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Gamma.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       WRITE(INORRK) (GAMMA2(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK) 
     
! CALCULATE REPRESENTATIVE OBSERVATION ERROR DIRECTLY USING EOFs IN UNRESOLVED SUBSPACE
     DO I=SS_DIM+1, SS_TOT
       CALL READEOF(I,2)
       DO J=1, NLOC
          HU(J,I) = STTEMP1(STLOC(J))*STSD(STLOC(J))       
       ENDDO
     ENDDO
     DO I=1, NLOC
       RRKSUM=0.0_DP
       DO J=SS_DIM+1, SS_TOT
         HUL(I,J)=RRKSUM+HU(I,J)*EVAL(J)
       ENDDO
     ENDDO
     
     DO I=1, NLOC
       DO J=1, NLOC
          RRKSUM=0.0_DP
	  DO K=SS_DIM+1,SS_TOT
             RRKSUM = RRKSUM + HUL(I,K)*HU(J,K)
          ENDDO
          R(I,J) = RRKSUM
       ENDDO
       
! ASSUME THAT THE MEASUREMENT ERROR IS 1% OF THE ONE STANDARD DEVIATION OF THE CONTROL RUN
       IF(RRK_OPTION == 0) THEN
          R(I,I) = R(I,I) + (STSD(STLOC(I))*0.01_DP)**2.0_DP
       ELSE
          R(I,I) = R(I,I) + (STSD(STLOC(I))/DSQRT(DBLE(KBM1))*0.01_DP)**2.0_DP	  
       ENDIF
     ENDDO
   
     DO I=1, NLOC
       DO J=1, I-1
         R(I,J) = 0.5_DP*(R(I,J)+R(J,I))
         R(J,I) = R(I,J)
       ENDDO
     ENDDO
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Robs.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, NLOC
       WRITE(INORRK) (R(I,J), I=1, NLOC)    
     ENDDO 
     CLOSE(INORRK)  
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Robs.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
     DO J=1, NLOC
       DO I=1, NLOC
         WRITE(INORRK,*) I,J, R(I,J)    
       ENDDO	 
     ENDDO 
     CLOSE(INORRK)   
     
! INVERT OBSERVATION ERROR COVARIANCE     
     CALL DGETRF(NLOC,NLOC,R,RRK_NOBSMAX,IPIV2,RCODE)
     IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'     
     call DGETRI(NLOC,R,RRK_NOBSMAX,IPIV2,WORK2,LWORK2,RCODE)
     IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'

! FORM BETA MATRIX = H_T*R**(-1)*H, WHERE H IS A NORMALIZED MATRIX H = H_ORIG*S
     CALL DGEMM('t','n',SS_DIM,NLOC,NLOC,1.0d0,HEOF,RRK_NOBSMAX,R,RRK_NOBSMAX,0.0d0,HTR,SS_DIM) 
     CALL DGEMM('n','n',SS_DIM,SS_DIM,NLOC,1.0d0,HTR,SS_DIM,HEOF,RRK_NOBSMAX,0.0d0,BETA2,SS_DIM)    
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Beta.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)     
     
     WRITE(IPT,*) 'Running the doubling algorithm......'
   
     DO IT =1, ILAST 
       
! READ GAMMA FROM FILE
       FILENAME =  TRIM(OUTDIR)//'/rrktemp/Gamma.dat'    
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)
                    
! READ BETA FROM FILE
       FILENAME =  TRIM(OUTDIR)//'/rrktemp/Beta.dat'    
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (BETA2(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)
      
! COMPUTE EYE + BETA*GAMMA = STORE TEMPORARILY IN RALPHA (RALPHA = EYE + BETA*GAMMA)
         
       DO I=1, SS_DIM
         DO J=1, SS_DIM
           RALPHA(I,J)=0.0_DP
         ENDDO
         RALPHA(I,I)=1.0_DP
       ENDDO     
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,BETA2,SS_DIM,GAMMA,SS_DIM,1.0d0,RALPHA,SS_DIM)
     
! COMPUTE INVERSE OF ABOBE RALPHA (=EYE+BETA*GAMMA) BY LAPACK ROUTINES AND STORE IT IN BETA (BETA=(EYE+BETA*GAMMA)**(-1))
       CALL DGETRF(SS_DIM,SS_DIM,RALPHA,SS_DIM,IPIV,RCODE)
  
       IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'
       CALL DGETRI(SS_DIM,RALPHA,SS_DIM,IPIV,WORK,LWORK,RCODE)
       IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'
       DO I=1, SS_DIM
         DO J=1, SS_DIM
           BETA2(I,J) = RALPHA(I,J)
	 ENDDO
       ENDDO
     
! READ RALPHA FROM FILE
       FILENAME =  TRIM(OUTDIR)//'/rrktemp/Alpha.dat'    
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)

! COMPUTE PRODUCT (EYE+BETA*GAMMA)**(-1)*RALPHA (GAMMA = BETA'*RALPHA)
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,BETA2,SS_DIM,RALPHA,SS_DIM,0.0d0,GAMMA,SS_DIM)  
          
! OUTPUT THIS PRODUCT TO TEMPORARY FILE
       FILENAME = TRIM(OUTDIR)//'/rrktemp/temp.02'
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         WRITE(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)     

! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1) AND STORE IN GAMMA (GAMMA = RALPHA*BETA')
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,BETA2,SS_DIM,0.0d0,GAMMA,SS_DIM) 
   
! READ BACK OLD FILE
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Beta.dat'
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (Beta2(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)  

! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1)*BETA (RALPHA' = GAMMA'*BETA')
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,GAMMA,SS_DIM,BETA2,SS_DIM,0.0d0,RALPHA,SS_DIM)

! READ BACK OLD RALPHA AND PUT IN GAMMA
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Alpha.dat'
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)  

! COMPUTE BETA+RALPHA*(EYE+BETA*GAMMA)**(-1)*BETA*RALPHA_T (BETA = BETA+RALPHA'*GAMMA_T)
       CALL DGEMM('n','t',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,GAMMA,SS_DIM,1.0d0,BETA2,SS_DIM)

! MAKE SURE BETA IS SYMMETRIC
       DO I=1, SS_DIM
          DO J=1, I-1
             BETA2(I,J) = 0.5_DP*(BETA2(I,J)+BETA2(J,I))
             BETA2(J,I) = BETA2(I,J) 
	  ENDDO
       ENDDO

! SAVE THIS NEW BETA TO FILE
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Beta.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
         DO J=1, SS_DIM
         WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)      

! READ BACK OLD GAMMA
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Gamma.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)     

! READ TEMPORARY FILE INTO RALPHA
       FILENAME = TRIM(OUTDIR)//'/rrktemp/temp.02'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)       

! COMPUTE GAMMA*(EYE+BETA*GAMMA)**(-1)*RALPHA (BETA'=GAMMA*RALPHA')
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,GAMMA,SS_DIM,RALPHA,SS_DIM,0.0d0,BETA2,SS_DIM)

! READ BACK OLD RALPHA
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Alpha.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)     

! COMPUTE GAMMA + RALPHA_T*GAMMA*(EYE+BETA*GAMMA)**(-1)*RALPHA (GAMMA = GAMMA+RALPHA_T*BETA')
       CALL DGEMM('t','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,BETA2,SS_DIM,1.0d0,GAMMA,SS_DIM)

! ENSURE GAMMA IS SYMMETRIC
       DO I=1, SS_DIM
          DO J=1, I-1
             GAMMA(I,J) = 0.5_DP*(GAMMA(I,J)+GAMMA(J,I))
             GAMMA(J,I) = GAMMA(I,J) 
	  ENDDO
       ENDDO

! SAVE THIS NEW GAMMA TO FILE
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Gamma.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         WRITE(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)      
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Pfctr.txt'
       OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
         DO J=1, SS_DIM
           DO I=1, SS_DIM
	      WRITE(INORRK,*) I,J,GAMMA(I,J)    
           ENDDO
	 ENDDO 
       CLOSE(INORRK)

! READ TEMPORARY FILE INTO GAMMA   
       FILENAME = TRIM(OUTDIR)//'/rrktemp/temp.02'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)
   
! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1)*RALPHA (BETA = RALPHA*GAMMA)  
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,GAMMA,SS_DIM,0.0d0,BETA2,SS_DIM)
   
! SAVE THIS NEW RALPHA TO FILE
       FILENAME = TRIM(OUTDIR)//'/rrktemp/Alpha.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
         DO J=1, SS_DIM
         WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)

!===================================================================================================
! END OF ITERATION FOR THE DOUBLING ALGORITHM          
     ENDDO   
!===================================================================================================

     WRITE(IPT,*) 'Setting up the Kalman gain matrix in the reduced space......'   

! READ IN PF FROM DOUBLING ALGORITHM 
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Gamma.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)
   
! READ IN OBSERVATION OPERATOR   
     FILENAME = TRIM(OUTDIR)//'/rrktemp/ObsOp.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       READ(INORRK) (HEOF(I,J), I=1, NLOC)    
     ENDDO 
     CLOSE(INORRK)
      
! READ OBSERVATION ERROR COVARIANCE      
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Robs.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, NLOC
       READ(INORRK) (R(I,J), I=1, NLOC)    
     ENDDO 
     CLOSE(INORRK)
   
! COMPUTE P*H^T (PHT = GAMMA * HEOF^T)   
     CALL DGEMM('n','t',SS_DIM,NLOC,SS_DIM,1.0d0,GAMMA,SS_DIM,HEOF,RRK_NOBSMAX,0.0d0,PHT,SS_DIM)

! COMPUTE H*P*H^T (R2 = H*P*HT)  
     CALL DGEMM('n','n',NLOC,NLOC,SS_DIM,1.0d0,HEOF,RRK_NOBSMAX,PHT,SS_DIM,0.0d0,R2,RRK_NOBSMAX)

! OUPUT BOTH FORCAST AND OBSERVATION ERROR STANDARD DEVIATION IN ASCII
     FILENAME =  TRIM(OUTDIR)//'/rrktemp/Fctobserr.txt'  
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
     DO I=1, NLOC
       WRITE(INORRK,*) STLOC(I),DSQRT(R2(I,I)),DSQRT(R(I,I))
     ENDDO
     CLOSE(INORRK)
   
! ADD R + (H P H^t) ---> R   
     DO I = 1, NLOC
        DO J=1, NLOC
           R(J,I) = R(J,I) + R2(J,I)
        ENDDO
     ENDDO
   
! INVERT R (H*P*H^T+R)   
     CALL DGETRF(NLOC,NLOC,R,RRK_NOBSMAX,IPIV2,RCODE)
     IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'
     CALL DGETRI(NLOC,R,RRK_NOBSMAX,IPIV2,WORK2,LWORK2,RCODE)
     IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'
     
! COMPUTE KALMAN GAIN: K = P*HT*(H*P*H^T+R)^(-1)   
     CALL DGEMM('n','n',SS_DIM,NLOC,NLOC,1.0d0,PHT,SS_DIM,R,RRK_NOBSMAX,0.0d0,KAL,SS_DIM)
   
! OUTPUT RESULT   
     OPEN(INORRK,FILE=KALNAME, FORM='UNFORMATTED') 
     DO J=1, NLOC
       WRITE(INORRK) (KAL(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)     
     FILENAME = TRIM(OUTDIR)//'/rrktemp/rrK.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
     DO J=1, NLOC
       DO I=1, SS_DIM
	 WRITE(INORRK,*) I,J,KAL(I,J)    
       ENDDO
     ENDDO 
     CLOSE(INORRK)
     
! OUTPUT KALMAN GAIN MATRIX IN THE FULL SPACE: D*ER*KAL
     DO J=1, SS_DIM
       CALL READEOF(J,2)
       DO I=1, STDIM
         SEOF(I,J)=STTEMP1(I)
       ENDDO
     ENDDO
     
     DO J=1, NLOC
       WRITE(STRCYC,'(I4.4)') J
       DO I=1, STDIM
          RRKSUM = 0.0_DP
          DO K=1, SS_DIM
             RRKSUM = RRKSUM + STSD(I)*SEOF(I,K)*KAL(K,J) 
          ENDDO 
	  STTEMP1(I) = RRKSUM
       ENDDO
       FILENAME = TRIM(OUTDIR)//'/rrktemp/'//'BigK'//STRCYC//'.txt'     
       OPEN(INORRK,FILE=FILENAME,FORM='FORMATTED')
       DO I=1, STDIM
          WRITE(INORRK,*) STTEMP1(I)
       ENDDO
       CLOSE(INORRK)
     ENDDO
     
! COMPUTE ANALYSIS ERROR COVARIANCE: P_AN1 = (I-K*H)*P (JUST AS DIAGNOSIS)
     CALL DGEMM('n','t',SS_DIM,SS_DIM,NLOC,-1.0d0,KAL,SS_DIM,PHT,SS_DIM,1.0d0,GAMMA,SS_DIM)
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Panlr.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
     DO J=1, SS_DIM
       DO I=1, SS_DIM
	 WRITE(INORRK,*) I,J,GAMMA(I,J)    
       ENDDO
     ENDDO 
     CLOSE(INORRK)
     
     CALL DGEMM('n','n',SS_DIM,SS_DIM,NLOC,1.0d0,KAL,SS_DIM,HEOF,RRK_NOBSMAX,0.0d0,RALPHA,SS_DIM)
     
! COMPUTE (I-KH)M AND ITS SINGULAR VALUE     
     DO I=1, SS_DIM
       DO J=1, SS_DIM
         IF(I/=J) THEN
            RALPHA(I,J) = (-1.0_DP)*RALPHA(I,J)
         ELSE
	        RALPHA(I,I) = 1.0_DP - RALPHA(I,I)
         ENDIF
       ENDDO
     ENDDO
     
     FILENAME = TRIM(OUTDIR)//'/rrktemp/I_KH.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
     DO J=1, SS_DIM
       DO I=1, SS_DIM
	 WRITE(INORRK,*) I,J,RALPHA(I,J)    
       ENDDO
     ENDDO 
     CLOSE(INORRK)     
     
     FILENAME = TRIM(OUTDIR)//'/rrktemp/Amatr.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       READ(INORRK) (TRANS(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)     
     CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,TRANS,SS_DIM,0.0d0,RALPHA,SS_DIM)
    
     FILENAME = TRIM(OUTDIR)//'/rrktemp/M_KHM.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
     DO J=1, SS_DIM
       DO I=1, SS_DIM
	 WRITE(INORRK,*) I,J,RALPHA(I,J)    
       ENDDO
     ENDDO 
     CLOSE(INORRK) 
     
     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,SEOF,STSD,HU,HUL,EVAL) 
     DEALLOCATE(STTEMP0,STTEMP1,STEOF,SDEOF,TRANS,RRKU,RRKV,RRKEL,RRKT,RRKS)
   ENDIF
   
   RETURN
   END SUBROUTINE RRK_RRK


! UTILITIES PROGRAMS
 
  SUBROUTINE RRK_ALLOC_VAR

   USE LIMS
   USE CONTROL
   USE RRKVAL
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif   
   IMPLICIT NONE
   
   INTEGER STDIM
   INTEGER SS_DIM
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   
   SS_DIM = RRK_NVAR*RRK_NEOF
   
! ALLOCATE ARRYS

   ALLOCATE(STTEMP0(STDIM))                ;STTEMP0   = ZERO
   ALLOCATE(STTEMP1(STDIM))                ;STTEMP1   = ZERO
   ALLOCATE(STEOF(STDIM))                  ;STEOF     = ZERO
   ALLOCATE(SDEOF(STDIM))                  ;SDEOF     = ZERO
   ALLOCATE(TRANS(SS_DIM,SS_DIM))          ;TRANS     = ZERO
   ALLOCATE(RRKEL(0:NNodeGL))                  ;RRKEL     = ZERO
   ALLOCATE(RRKU(0:MElemGL,KB))                ;RRKU      = ZERO
   ALLOCATE(RRKV(0:MElemGL,KB))                ;RRKV      = ZERO
   ALLOCATE(RRKT(0:NNodeGL,KB))                ;RRKT      = ZERO
   ALLOCATE(RRKS(0:NNodeGL,KB))                ;RRKS      = ZERO
   
! END ALLOCATION
   
   RETURN
  END SUBROUTINE RRK_ALLOC_VAR

!------------------------------------------------------------------------------|
!  PERTURB THE BASE STATE IN THE EIGENVECTOR DIRECTIONS                        |
!------------------------------------------------------------------------------|   
   
  SUBROUTINE PERTURB(IFILE)

#  if defined(WET_DRY)
     USE MOD_WD
#  endif
   USE LIMS
   USE ALL_VARS
   USE RRKVAL
   USE CONTROL
#  if defined (MULTIPROCESSOR)
     USE MOD_PAR
#  endif
   IMPLICIT NONE
    
    INTEGER SS_DIM
    INTEGER STDIM
    INTEGER I_EOF
    INTEGER II,I,J
    INTEGER IDUMMY
    INTEGER I_START
    CHARACTER(LEN=80) IFILE
    CHARACTER(LEN=80) IFILE2
    CHARACTER(LEN=80) TEMP1FILE
    CHARACTER(LEN=80) TEMP2FILE
    CHARACTER(LEN=4)  FEOF
    REAL(DP) ::  ELSD, USD, TSD, SSD
  
    STDIM = 0
    IF(EL_ASSIM) STDIM = STDIM + NNodeGL
    IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
    IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
    IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
    
    SS_DIM = RRK_NVAR*RRK_NEOF

    IFILE2=TRIM(OUTDIR)//'/rrktemp/avgstd.dat'
    OPEN(INORRK,FILE=IFILE2, STATUS='OLD')
   
    READ(INORRK,*)
    READ(INORRK,*) ELSD, USD, TSD, SSD
    CLOSE(INORRK)
   
    IDUMMY  = 0
    IF(EL_ASSIM) THEN
      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        SDEOF(IDUMMY) = ELSD
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN
      DO I=1, 2*KBM1
        DO J=1, MElemGL 
          IDUMMY = IDUMMY + 1
          IF(RRK_OPTION == 0) THEN 
            SDEOF(IDUMMY)=USD
          ELSE
            SDEOF(IDUMMY)=USD*DSQRT(DBLE(KBM1))     
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL 
          IDUMMY = IDUMMY + 1
          SDEOF(IDUMMY)=TSD
        ENDDO
      ENDDO    
    ENDIF
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL 
          IDUMMY = IDUMMY + 1
          SDEOF(IDUMMY)=SSD
        ENDDO
      ENDDO    
    ENDIF
    
    
    CALL READRESTART(IFILE,I_START)
    CALL GR2ST(0)
    
! READ THE EOF AND PERTURB THE BASE STATE IN THIS DIRECTION

   DO I_EOF=1, SS_DIM
       CALL READEOF(I_EOF,1) 
       WRITE(FEOF,'(I4.4)') I_EOF
       TEMP1FILE = TRIM(OUTDIR)//'/rrktemp/'//'eofini'//FEOF
#      if defined(WET_DRY)
       TEMP2FILE = TRIM(OUTDIR)//'/rrktemp/'//'eofini'//FEOF//'_wd' 
#      endif

! PERTURB THE BASE STATE IN THE DIRECTION OF THE I_EOF'TH EOF AND STORE IT IN THE MP1FILE (FOR RESTART)
       DO I=1, STDIM
          STTEMP1(I) = STTEMP0(I) + STEOF(I)*SDEOF(I)*DBLE(RRK_PSIZE)
       ENDDO
       CALL ST2GR

   IF(SERIAL) THEN
     EL =  RRKEL 
!     D  = H + EL
!     ET = EL
!     DT = D
     DO I=1, MElemGL
       EL1(I)=(EL(NVG(I,1)) + EL(NVG(I,2)) + EL(NVG(I,3)) )/3.0_DP
     ENDDO
!     D1  = H1 + EL1
!     ET1 = EL1
!     DT1 = D1

     U  =  RRKU
     V  =  RRKV
!     DO I=1, MElemGL 
!       UA(I)  =  RRKU(I,1) 
!       VA(I)  =  RRKV(I,1)       
!     ENDDO
     T1 =  RRKT
     S1 =  RRKS
   ELSE 
#  if defined (MULTIPROCESSOR)
     DO I=1, NNode
        EL(I)=RRKEL(NGID(I))
     ENDDO
!     D  = H + EL
!     ET = EL
!     DT = D
     DO I=1, MElem
       EL1(I)=(EL(NV(I,1)) + EL(NV(I,2)) + EL(NV(I,3)) )/3.0_DP
     ENDDO
!     D1  = H1 + EL1
!     ET1 = EL1
!     DT1 = D1

     DO I=1, MElem
       DO J=1, KB
         U(I,J)=RRKU(EGID(I),J)
         V(I,J)=RRKV(EGID(I),J)
       ENDDO
     ENDDO
!     DO I=1, MElem
!       UA(I)=RRKU(EGID(I),1)
!       VA(I)=RRKV(EGID(I),1)
!     ENDDO

     DO I=1, NNode
       DO J=1, KB
         T1(I,J)=RRKT(NGID(I),J)
         S1(I,J)=RRKS(NGID(I),J)
       ENDDO
     ENDDO
#  endif
   ENDIF

#    if defined(WET_DRY)

     CALL WET_JUDGE_EL 
  
     OPEN(INORRK,FILE=TEMP2FILE,FORM='FORMATTED')
     CALL WD_DUMP_EL(INORRK,I_START)       

#    endif 

! WRITE THE PERTURBED STATE IN TEMP1FILE
     CALL WRITERESTART(TEMP1FILE,I_START)

   ENDDO
  
! READ THE EOF AND PERTURB THE BASE STATE IN THE NEGATIVE DIRECTION

   DO I_EOF=1, SS_DIM
     CALL READEOF(I_EOF,1) 
     WRITE(FEOF,'(I4.4)') I_EOF
     TEMP1FILE = TRIM(OUTDIR)//'/rrktemp/'//'eofini'//FEOF//'n'
#    if defined(WET_DRY)
     TEMP2FILE = TRIM(OUTDIR)//'/rrktemp/'//'eofini'//FEOF//'n_wd' 
#    endif

! PERTURB THE BASE STATE IN THE DIRECTION OF THE I_EOF'TH EOF AND STORE IT IN THE MP1FILE (FOR RESTART)
    
     DO I=1, STDIM
       STTEMP1(I) = STTEMP0(I) - STEOF(I)*SDEOF(I)*DBLE(RRK_PSIZE)
     ENDDO
     CALL ST2GR

   IF(SERIAL) THEN
     EL =  RRKEL 
!     D  = H + EL
!     ET = EL
!     DT = D
     DO I=1, MElemGL
       EL1(I)=(EL(NVG(I,1)) + EL(NVG(I,2)) + EL(NVG(I,3)) )/3.0_DP
     ENDDO
!     D1  = H1 + EL1
!     ET1 = EL1
!     DT1 = D1

     U  =  RRKU 
     V  =  RRKV       
!     DO I=1, MElemGL 
!       UA(I)  =  RRKU(I,1) 
!       VA(I)  =  RRKV(I,1)       
!     ENDDO
     T1 =  RRKT
     S1 =  RRKS
   ELSE 
#  if defined (MULTIPROCESSOR)
     DO I=1, NNode
        EL(I)=RRKEL(NGID(I))
     ENDDO
!     D  = H + EL
!     ET = EL
!     DT = D
     DO I=1, MElem
       EL1(I)=(EL(NV(I,1)) + EL(NV(I,2)) + EL(NV(I,3)) )/3.0_DP
     ENDDO
!     D1  = H1 + EL1
!     ET1 = EL1
!     DT1 = D1

     DO I=1, MElem
       DO J=1, KB
         U(I,J)=RRKU(EGID(I),J)
         V(I,J)=RRKV(EGID(I),J)
       ENDDO
     ENDDO
!     DO I=1, MElem
!       UA(I)=RRKU(EGID(I),1)
!       VA(I)=RRKV(EGID(I),1)
!     ENDDO
     
     DO I=1, NNode
       DO J=1, KB
         T1(I,J)=RRKT(NGID(I),J)
         S1(I,J)=RRKS(NGID(I),J)
       ENDDO
     ENDDO
#  endif
   ENDIF

#    if defined(WET_DRY)

     CALL WET_JUDGE_EL 
     
     OPEN(INORRK,FILE=TEMP2FILE,FORM='FORMATTED')
     CALL WD_DUMP_EL(INORRK,I_START)

#    endif 

! WRITE THE PERTURBED STATE IN TEMP1FILE
     CALL WRITERESTART(TEMP1FILE,I_START)

  ENDDO    

  RETURN
  END SUBROUTINE PERTURB

!------------------------------------------------------------------------------/
!  PERTURB THE BASE STATE IN THE EIGENVECTOR DIRECTIONS WITH THE MAGNITUDE     /
!  OF EIGENVALUES                                                              /                          
!------------------------------------------------------------------------------/   
   
  SUBROUTINE PERTURBL(IFILE)

#  if defined(WET_DRY)
     USE MOD_WD
#  endif
   USE LIMS
   USE ALL_VARS
   USE RRKVAL
   USE CONTROL
#  if defined (MULTIPROCESSOR)
     USE MOD_PAR
#  endif
   IMPLICIT NONE
    
    INTEGER SS_DIM
    INTEGER STDIM
    INTEGER I_EOF
    INTEGER II,I,J
    INTEGER I_START
    INTEGER IDUMMY
    CHARACTER(LEN=80) IFILE
    CHARACTER(LEN=80) IFILE2
    CHARACTER(LEN=80) TEMP1FILE
    CHARACTER(LEN=80) TEMP2FILE
    CHARACTER(LEN=4)  FEOF
    REAL(DP) ::  ELSD, USD, TSD, SSD
    REAL(DP),ALLOCATABLE :: EVAL(:)
   
    STDIM = 0
    IF(EL_ASSIM) STDIM = STDIM + NNodeGL
    IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
    IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
    IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
    
    SS_DIM = RRK_NVAR*RRK_NEOF

    ALLOCATE(EVAL(SS_DIM))            ; EVAL = ZERO

    IFILE2=TRIM(OUTDIR)//'/rrktemp/avgstd.dat'
    OPEN(INORRK,FILE=IFILE2, STATUS='OLD')
   
    READ(INORRK,*)
    READ(INORRK,*) ELSD, USD, TSD, SSD
    CLOSE(INORRK)
   
    IDUMMY = 0
    IF(EL_ASSIM) THEN
      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        SDEOF(IDUMMY) = ELSD
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN
      DO I=1, 2*KBM1
        DO J=1, MElemGL 
          IDUMMY = IDUMMY + 1
          IF(RRK_OPTION == 0) THEN 
            SDEOF(IDUMMY)=USD
          ELSE
            SDEOF(IDUMMY)=USD*DSQRT(DBLE(KBM1))     
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL 
          IDUMMY = IDUMMY + 1
          SDEOF(IDUMMY)=TSD
        ENDDO
      ENDDO    
    ENDIF
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL 
          IDUMMY = IDUMMY + 1
          SDEOF(IDUMMY)=SSD
        ENDDO
      ENDDO    
    ENDIF

! READ THE EIGENVALUES

    IFILE2=TRIM(OUTDIR)//'/rrktemp/eigenvalue.dat'
    OPEN(INORRK,FILE=IFILE2, STATUS='OLD')
   
    DO I=1, SS_DIM
      READ(INORRK,*) IDUMMY, EVAL(I)
    ENDDO
    CLOSE(INORRK)
   
    CALL READRESTART(IFILE,I_START)
    CALL GR2ST(0)
   
! READ THE EOF AND PERTURB THE BASE STATE IN THIS DIRECTION

   DO I_EOF=1, SS_DIM
     CALL READEOF(I_EOF,1) 
     WRITE(FEOF,'(I4.4)') I_EOF
     TEMP1FILE = TRIM(OUTDIR)//'/rrktemp/'//'eofiniL'//FEOF
#    if defined(WET_DRY)
     TEMP2FILE = TRIM(OUTDIR)//'/rrktemp/'//'eofiniL'//FEOF//'_wd' 
#    endif

! PERTURB THE BASE STATE IN THE DIRECTION OF THE I_EOF'TH EOF AND STORE IT IN THE MP1FILE (FOR RESTART)
    
     DO I=1, STDIM
        STTEMP1(I) = STTEMP0(I) + STEOF(I)*SDEOF(I)*DBLE(RRK_PSIZE)*EVAL(I_EOF)
     ENDDO
     CALL ST2GR

   IF(SERIAL) THEN
     EL =  RRKEL 
!     D  = H + EL
!     ET = EL
!     DT = D
     DO I=1, MElemGL
       EL1(I)=(EL(NVG(I,1)) + EL(NVG(I,2)) + EL(NVG(I,3)) )/3.0_DP
     ENDDO
!     D1  = H1 + EL1
!     ET1 = EL1
!     DT1 = D1

     U  =  RRKU 
     V  =  RRKV       
!     DO I=1, MElemGL 
!       UA(I)  =  RRKU(I,1) 
!       VA(I)  =  RRKV(I,1)       
!     ENDDO
     T1 =  RRKT
     S1 =  RRKS
   ELSE 
#  if defined (MULTIPROCESSOR)
     DO I=1, NNode
        EL(I)=RRKEL(NGID(I))
     ENDDO
!     D  = H + EL
!     ET = EL
!     DT = D
     DO I=1, MElem
       EL1(I)=(EL(NV(I,1)) + EL(NV(I,2)) + EL(NV(I,3)) )/3.0_DP
     ENDDO
!     D1  = H1 + EL1
!     ET1 = EL1
!     DT1 = D1

     DO I=1, MElem
       DO J=1, KB
         U(I,J)=RRKU(EGID(I),J)
         V(I,J)=RRKV(EGID(I),J)
       ENDDO
     ENDDO
!     DO I=1, MElem
!       UA(I)=RRKU(EGID(I),1)
!       VA(I)=RRKV(EGID(I),1)
!     ENDDO

     DO I=1, NNode
       DO J=1, KB
         T1(I,J)=RRKT(NGID(I),J)
         S1(I,J)=RRKS(NGID(I),J)
       ENDDO
     ENDDO
#  endif
   ENDIF 

#    if defined(WET_DRY)

     CALL WET_JUDGE_EL 
     
     OPEN(INORRK,FILE=TEMP2FILE,FORM='FORMATTED')
     CALL WD_DUMP_EL(INORRK,I_START)
     
#    endif 

! WRITE THE PERTURBED STATE IN TEMP1FILE
     CALL WRITERESTART(TEMP1FILE,I_START)

   ENDDO
  
! READ THE EOF AND PERTURB THE BASE STATE IN THE NEGATIVE DIRECTION

   DO I_EOF=1, SS_DIM
     CALL READEOF(I_EOF,1) 
     WRITE(FEOF,'(I4.4)') I_EOF
     TEMP1FILE = TRIM(OUTDIR)//'/rrktemp/'//'eofiniL'//FEOF//'n'
#    if defined(WET_DRY)
     TEMP2FILE = TRIM(OUTDIR)//'/rrktemp/'//'eofiniL'//FEOF//'n_wd' 
#    endif

! PERTURB THE BASE STATE IN THE DIRECTION OF THE I_EOF'TH EOF AND STORE IT IN THE MP1FILE (FOR RESTART)
    
     DO I=1, STDIM
        STTEMP1(I) = STTEMP0(I) - STEOF(I)*SDEOF(I)*DBLE(RRK_PSIZE)*EVAL(I_EOF)
     ENDDO
     CALL ST2GR

   IF(SERIAL) THEN
     EL =  RRKEL 
!     D  = H + EL
!     ET = EL
!     DT = D
     DO I=1, MElemGL
       EL1(I)=(EL(NVG(I,1)) + EL(NVG(I,2)) + EL(NVG(I,3)) )/3.0_DP
     ENDDO
!     D1  = H1 + EL1
!     ET1 = EL1
!     DT1 = D1

     U  =  RRKU 
     V  =  RRKV
!     DO I=1, MElemGL 
!       UA(I)  =  RRKU(I,1) 
!       VA(I)  =  RRKV(I,1)       
!     ENDDO
     T1 =  RRKT
     S1 =  RRKS       
   ELSE 
#  if defined (MULTIPROCESSOR)
     DO I=1, NNode
        EL(I)=RRKEL(NGID(I))
     ENDDO
!     D  = H + EL
!     ET = EL
!     DT = D
     DO I=1, MElem
       EL1(I)=(EL(NV(I,1)) + EL(NV(I,2)) + EL(NV(I,3)) )/3.0_DP
     ENDDO
!     D1  = H1 + EL1
!     ET1 = EL1
!     DT1 = D1

     DO I=1, MElem
       DO J=1, KB
         U(I,J)=RRKU(EGID(I),J)
         V(I,J)=RRKV(EGID(I),J)
       ENDDO
     ENDDO
!     DO I=1, MElem
!       UA(I)=RRKU(EGID(I),1)
!       VA(I)=RRKV(EGID(I),1)
!     ENDDO

     DO I=1, NNode
       DO J=1, KB
         T1(I,J)=RRKT(NGID(I),J)
         S1(I,J)=RRKS(NGID(I),J)
       ENDDO
     ENDDO
#  endif
   ENDIF

#    if defined(WET_DRY)

     CALL WET_JUDGE_EL 
     
     OPEN(INORRK,FILE=TEMP2FILE,FORM='FORMATTED')
     CALL WD_DUMP_EL(INORRK,I_START)
     
#    endif 

! WRITE THE PERTURBED STATE IN TEMP1FILE
     CALL WRITERESTART(TEMP1FILE,I_START)

   ENDDO    
   
   DEALLOCATE(EVAL) 
  RETURN
  END SUBROUTINE PERTURBL

!------------------------------------------------------------------------------|
!  CALCULATE THE LINEARIZED MODEL IN THE REDUCED SPACE                         |
!------------------------------------------------------------------------------|   
   
  SUBROUTINE MREDUCED(IFILE)

   USE MOD_WD
   USE LIMS
   USE CONTROL
   USE RRKVAL
   IMPLICIT NONE
    
   INTEGER SS_DIM
   INTEGER STDIM
   INTEGER I_EOF
   INTEGER I,J
   INTEGER IDUMMY
   INTEGER I_START   
   CHARACTER(LEN=80) IFILE
   CHARACTER(LEN=80) IFILE2
   CHARACTER(LEN=80) TEMPFILE
   CHARACTER(LEN=4)  FEOF
   REAL(DP) ::  ELSD, USD, TSD, SSD
   REAL(DP) ::  SUM0

   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + NNodeGL
   IF(UV_ASSIM) STDIM = STDIM + 2*MElemGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + NNodeGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + NNodeGL*KBM1

   SS_DIM = RRK_NVAR*RRK_NEOF

! READ ONE STANDARD DEVIATION OF THE CONTROL RUN    
   IFILE2=TRIM(OUTDIR)//'/rrktemp/avgstd.dat'
   OPEN(INORRK,FILE=IFILE2, STATUS='OLD')
   
   READ(INORRK,*)
   READ(INORRK,*) ELSD, USD, TSD, SSD
   CLOSE(INORRK)
   
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1, NNodeGL
       IDUMMY = IDUMMY + 1
       SDEOF(IDUMMY) = ELSD
     ENDDO
   ENDIF
   IF(UV_ASSIM) THEN
     DO I =1, 2*KBM1
       DO J=1, MElemGL
         IDUMMY = IDUMMY + 1
         IF(RRK_OPTION == 0) THEN
            SDEOF(IDUMMY) = USD
         ELSE
            SDEOF(IDUMMY) = USD*DSQRT(DBLE(KBM1))
         ENDIF
       ENDDO
     ENDDO
   ENDIF
   IF(T_ASSIM) THEN
     DO I =1, KBM1
       DO J=1, NNodeGL
         IDUMMY = IDUMMY + 1
         SDEOF(IDUMMY) = TSD
       ENDDO
     ENDDO
   ENDIF
   IF(S_ASSIM) THEN
     DO I =1, KBM1
       DO J=1, NNodeGL
         IDUMMY = IDUMMY + 1
         SDEOF(IDUMMY) = SSD
       ENDDO
     ENDDO
   ENDIF   
   
   
   CALL READRESTART(IFILE,I_START)  
   CALL GR2ST(0)
   
! READ THE FORECAST WHICH IS PROPAGATED FROM THE PERTURBED STATE   
   DO I_EOF =1, SS_DIM
     WRITE(FEOF,'(I4.4)') I_EOF
     TEMPFILE=TRIM(OUTDIR)//'/rrktemp/'//'eoffct'//FEOF   

! READ THE EACH PERTURBED STATE AT THE END OF LINEARUZATION TIME STEP   
   CALL READRESTART(TEMPFILE,I_START)   
   CALL GR2ST(1)
   
! PROJECT  EVOLVED PERTURBATION ONTO EOFs
     DO I=1, SS_DIM
       CALL READEOF(I,1)
       SUM0 = 0.0_DP
       DO J=1, STDIM
         SUM0 = SUM0 + STEOF(J)/SDEOF(J)*(STTEMP1(J)-STTEMP0(J))
       ENDDO
       TRANS(I,I_EOF) = SUM0/DBLE(RRK_PSIZE)
           
     ENDDO
! EDN OF LOOP OVER EACH EOF   
   ENDDO
   
  END SUBROUTINE MREDUCED 


!=====================================================================================/
!   READ THE STATE FROM THE RESTART FILE                                              /
!=====================================================================================/
  SUBROUTINE READRESTART(INFILE,I_START)

   USE LIMS
   USE RRKVAL
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE

    INTEGER I,K,N1
    INTEGER I_START
    INTEGER INF
    CHARACTER(LEN=80) INFILE
    REAL(SP), ALLOCATABLE :: RTP(:)
    
    ALLOCATE(RTP(MElemGL))          ; RTP   = 0.0_SP
        
    INF = INORRK
    OPEN(INF,FILE=INFILE, FORM='UNFORMATTED') 

    IF(SERIAL)THEN
     REWIND(INF)
     READ(INF) I_START
     READ(INF) ((U(I,K),K=1,KB),I=0,MElem)
     READ(INF) ((V(I,K),K=1,KB),I=0,MElem)
     READ(INF) ((W(I,K),K=1,KB),I=0,MElem)
#    if defined (GOTM)
     READ(INF) ((TKE(I,K),K=1,KB),I=0,NNode)
     READ(INF) ((TEPS(I,K),K=1,KB),I=0,NNode)
#    else
     READ(INF) ((Q2(I,K),K=1,KB),I=0,NNode)
     READ(INF) ((Q2L(I,K),K=1,KB),I=0,NNode)
     READ(INF) ((L(I,K),K=1,KB),I=0,NNode)
#    endif
     READ(INF) ((S(I,K),K=1,KB),I=0,MElem)
     READ(INF) ((T(I,K),K=1,KB),I=0,MElem)
     READ(INF) ((RHO(I,K),K=1,KB),I=0,MElem)
     READ(INF) ((TMEAN(I,K),K=1,KB),I=0,MElem)
     READ(INF) ((SMEAN(I,K),K=1,KB),I=0,MElem)
     READ(INF) ((RMEAN(I,K),K=1,KB),I=0,MElem)

     READ(INF) ((S1(I,K),K=1,KB),I=1,NNode)
     READ(INF) ((T1(I,K),K=1,KB),I=1,NNode)
     READ(INF) ((RHO1(I,K),K=1,KB),I=1,NNode)
     READ(INF) ((TMEAN1(I,K),K=1,KB),I=1,NNode)
     READ(INF) ((SMEAN1(I,K),K=1,KB),I=1,NNode)
     READ(INF) ((RMEAN1(I,K),K=1,KB),I=1,NNode)
     READ(INF) ((KM(I,K),K=1,KB),I=1,NNode)
     READ(INF) ((KH(I,K),K=1,KB),I=1,NNode)
     READ(INF) ((KQ(I,K),K=1,KB),I=1,NNode)

     READ(INF) (UA(I), I=0,MElem)
     READ(INF) (VA(I), I=0,MElem)
     READ(INF) (EL1(I), I=1,MElem)
     READ(INF) (ET1(I), I=1,MElem)
     READ(INF) (H1(I), I=1,MElem)
     READ(INF) (D1(I), I=1,MElem)
     READ(INF) (DT1(I), I=1,MElem)
     READ(INF) (RTP(I), I=1,MElem)

     READ(INF) (EL(I), I=1,NNode)
     READ(INF) (ET(I), I=1,NNode)
     READ(INF) (H(I), I=1,NNode)
     READ(INF) (D(I), I=1,NNode)
     READ(INF) (DT(I), I=1,NNode)

#    if defined (WATER_QUALITY)
     DO N1=1,NB
       READ(INF) ((WQM(I,K,N1),K=1,KB),I=1,NNode)
     END DO
#    endif

     CLOSE(INF)
   ELSE
#  if defined (MULTIPROCESSOR)
     REWIND(INF)
     READ(INF) I_START
     READ(INF) ((RRKU(I,K),K=1,KB),I=0,MElemGL)
     READ(INF) ((RRKV(I,K),K=1,KB),I=0,MElemGL)
     CALL PREAD(INF,W     ,LBOUND(W,1),    UBOUND(W,1),    MElem,MElemGL,KB,EGID(1),0,"W"     )
#    if defined (GOTM)
     CALL PREAD(INF,TKE    ,LBOUND(TKE,1),   UBOUND(TKE,1),   NNode,MElemGL,KB,EGID(1),0,"TKE"    )
     CALL PREAD(INF,TEPS   ,LBOUND(TEPS,1),  UBOUND(TEPS,1),  NNode,MElemGL,KB,EGID(1),0,"TEPS"   )
#    else
     CALL PREAD(INF,Q2    ,LBOUND(Q2,1),   UBOUND(Q2,1),   NNode,MElemGL,KB,EGID(1),0,"Q2"    )
     CALL PREAD(INF,Q2L   ,LBOUND(Q2L,1),  UBOUND(Q2L,1),  NNode,MElemGL,KB,EGID(1),0,"Q2L"   )
     CALL PREAD(INF,L     ,LBOUND(L,1  ),  UBOUND(L,1),    NNode,MElemGL,KB,EGID(1),0,"L"   )
#    endif
     CALL PREAD(INF,S     ,LBOUND(S,1),    UBOUND(S,1),    MElem,MElemGL,KB,EGID(1),0,"S"     )
     CALL PREAD(INF,T     ,LBOUND(S,1),    UBOUND(T,1),    MElem,MElemGL,KB,EGID(1),0,"T"     )
     CALL PREAD(INF,RHO   ,LBOUND(RHO,1),  UBOUND(RHO,1),  MElem,MElemGL,KB,EGID(1),0,"RHO"   )
     CALL PREAD(INF,TMEAN ,LBOUND(TMEAN,1),UBOUND(TMEAN,1),MElem,MElemGL,KB,EGID(1),0,"TMEAN" )
     CALL PREAD(INF,SMEAN ,LBOUND(SMEAN,1),UBOUND(SMEAN,1),MElem,MElemGL,KB,EGID(1),0,"SMEAN" )
     CALL PREAD(INF,RMEAN ,LBOUND(RMEAN,1),UBOUND(RMEAN,1),MElem,MElemGL,KB,EGID(1),0,"RMEAN" )

     READ(INF) ((RRKS(I,K),K=1,KB),I=1,NNodeGL)
     READ(INF) ((RRKT(I,K),K=1,KB),I=1,NNodeGL)
     CALL PREAD(INF,RHO1  ,LBOUND(RHO1,1),  UBOUND(RHO1,1),  NNode,NNodeGL,KB,NGID,1,"RHO1"   )
     CALL PREAD(INF,TMEAN1,LBOUND(TMEAN1,1),UBOUND(TMEAN1,1),NNode,NNodeGL,KB,NGID,1,"TMEAN1" )
     CALL PREAD(INF,SMEAN1,LBOUND(SMEAN1,1),UBOUND(SMEAN1,1),NNode,NNodeGL,KB,NGID,1,"SMEAN1" )
     CALL PREAD(INF,RMEAN1,LBOUND(RMEAN1,1),UBOUND(RMEAN1,1),NNode,NNodeGL,KB,NGID,1,"RMEAN1" )
  
     CALL PREAD(INF,KM  ,LBOUND(KM,1),UBOUND(KM,1),NNode,MElemGL,KB,EGID(1),1,"KM" )
     CALL PREAD(INF,KH  ,LBOUND(KH,1),UBOUND(KH,1),NNode,MElemGL,KB,EGID(1),1,"KH" )
     CALL PREAD(INF,KQ  ,LBOUND(KQ,1),UBOUND(KQ,1),NNode,MElemGL,KB,EGID(1),1,"KQ" )

     READ(INF) ((RRKU(I,K),K=1,1),I=0,MElemGL)
     READ(INF) ((RRKV(I,K),K=1,1),I=0,MElemGL)
     CALL PREAD(INF,EL1 ,LBOUND(EL1,1),UBOUND(EL1,1),MElem,MElemGL,1 ,EGID(1),1,"EL1" )
     CALL PREAD(INF,ET1 ,LBOUND(ET1,1),UBOUND(ET1,1),MElem,MElemGL,1 ,EGID(1),1,"ET1" )
     CALL PREAD(INF,H1  ,LBOUND(H1,1), UBOUND(H1,1), MElem,MElemGL,1 ,EGID(1),1,"H1"  )
     CALL PREAD(INF,D1  ,LBOUND(D1,1), UBOUND(D1,1), MElem,MElemGL,1 ,EGID(1),1,"D1"  )
     CALL PREAD(INF,DT1 ,LBOUND(DT1,1),UBOUND(DT1,1),MElem,MElemGL,1 ,EGID(1),1,"DT1" )
     CALL PREAD(INF,RTP ,LBOUND(RTP,1),UBOUND(RTP,1),MElem,MElemGL,1 ,EGID(1),1,"RTP" )

     READ(INF) (RRKEL(I), I=1,NNodeGL)
     CALL PREAD(INF,ET  ,LBOUND(ET,1),UBOUND(ET,1),NNode,NNodeGL,1 ,NGID,1,"ET"   )
     CALL PREAD(INF,H   ,LBOUND(H,1), UBOUND(H,1), NNode,NNodeGL,1 ,NGID,1,"H"    )
     CALL PREAD(INF,D   ,LBOUND(D,1), UBOUND(D,1), NNode,NNodeGL,1 ,NGID,1,"D"    )
     CALL PREAD(INF,DT  ,LBOUND(DT,1),UBOUND(DT,1),NNode,NNodeGL,1 ,NGID,1,"DT"   )

#    if defined (WATER_QUALITY)
     DO N1 = 1, NB
       CALL PREAD(INF,WQM(:,:,N1),LBOUND(WQM(:,:,N1),1),UBOUND(WQM(:,:,N1),1), &
                  NNode,NNodeGL,KB,NGID,1,"WQM")
     END DO
#    endif

     CLOSE(INF)
#    endif
   END IF

   IF(SERIAL) THEN
     RRKEL = EL
     RRKU  = U
     RRKV  = V
     DO I=1, MElemGL 
       RRKU(I,1)  = UA(I)
       RRKV(I,1)  = VA(I)
     ENDDO
     RRKT  = T1
     RRKS  = S1
   ENDIF

    DEALLOCATE(RTP)
  RETURN
  END SUBROUTINE READRESTART

!=====================================================================================/
!   WRITE THE STATE IN THE RESTART FILE                                              /
!=====================================================================================/
  SUBROUTINE WRITERESTART(OFILE,I_START)

   USE LIMS
   USE RRKVAL
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif   
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif

   IMPLICIT NONE

    INTEGER I,K,N1,ME,NPC
    INTEGER I_START
    INTEGER INF
    CHARACTER(LEN=80) OFILE
    REAL(SP), ALLOCATABLE :: RTP(:)
    
    ME = MYID ; NPC = NPROCS
    ALLOCATE(RTP(MElemGL)) ; RTP = 0.0_SP
        
    INF = INORRK
    IF(MSR)THEN
     OPEN(INF,FILE=TRIM(OFILE), FORM='UNFORMATTED') 
     REWIND(INF)
     WRITE(INF) I_START
    END IF

    IF(SERIAL)THEN
     WRITE(INF) ((U(I,K),    K=1,KB),I=0,MElem)
     WRITE(INF) ((V(I,K),    K=1,KB),I=0,MElem)
     WRITE(INF) ((W(I,K),    K=1,KB),I=0,MElem)
#    if defined (GOTM)
     WRITE(INF) ((TKE(I,K),   K=1,KB),I=0,NNode)
     WRITE(INF) ((TEPS(I,K),  K=1,KB),I=0,NNode)
#    else
     WRITE(INF) ((Q2(I,K),   K=1,KB),I=0,NNode)
     WRITE(INF) ((Q2L(I,K),  K=1,KB),I=0,NNode)
     WRITE(INF) ((L(I,K)  ,  K=1,KB),I=0,NNode)
#    endif
     WRITE(INF) ((S(I,K),    K=1,KB),I=0,MElem)
     WRITE(INF) ((T(I,K),    K=1,KB),I=0,MElem)
     WRITE(INF) ((RHO(I,K),  K=1,KB),I=0,MElem)
     WRITE(INF) ((TMEAN(I,K),K=1,KB),I=0,MElem)
     WRITE(INF) ((SMEAN(I,K),K=1,KB),I=0,MElem)
     WRITE(INF) ((RMEAN(I,K),K=1,KB),I=0,MElem)

     WRITE(INF) ((S1(I,K),    K=1,KB),I=1,NNode)
     WRITE(INF) ((T1(I,K),    K=1,KB),I=1,NNode)
     WRITE(INF) ((RHO1(I,K),  K=1,KB),I=1,NNode)
     WRITE(INF) ((TMEAN1(I,K),K=1,KB),I=1,NNode)
     WRITE(INF) ((SMEAN1(I,K),K=1,KB),I=1,NNode)
     WRITE(INF) ((RMEAN1(I,K),K=1,KB),I=1,NNode)

     WRITE(INF) ((KM(I,K),K=1,KB),I=1,NNode)
     WRITE(INF) ((KH(I,K),K=1,KB),I=1,NNode)
     WRITE(INF) ((KQ(I,K),K=1,KB),I=1,NNode)

     WRITE(INF) (UA(I), I=0,MElem)
     WRITE(INF) (VA(I), I=0,MElem)

     WRITE(INF) (EL1(I), I=1,MElem)
     WRITE(INF) (ET1(I), I=1,MElem)
     WRITE(INF) (H1(I),  I=1,MElem)
     WRITE(INF) (D1(I),  I=1,MElem)
     WRITE(INF) (DT1(I), I=1,MElem)
     WRITE(INF) (RTP(I), I=1,MElem)

     WRITE(INF) (EL(I), I=1,NNode)
     WRITE(INF) (ET(I), I=1,NNode)
     WRITE(INF) (H(I),  I=1,NNode)
     WRITE(INF) (D(I),  I=1,NNode)
     WRITE(INF) (DT(I), I=1,NNode)

#    if defined (WATER_QUALITY)
     DO N1 = 1, NB
       WRITE(INF) ((WQM(I,K,N1),K=1,KB),I=1,NNode)
     END DO
#    endif

   ELSE
#     if defined (MULTIPROCESSOR)
      CALL PWRITE(INF,ME,NPC,U,    LBOUND(U,1),    UBOUND(U,1),    MElem,MElemGL,KB,EMAP,0,"U"    )
      CALL PWRITE(INF,ME,NPC,V,    LBOUND(V,1),    UBOUND(V,1),    MElem,MElemGL,KB,EMAP,0,"V"    )
      CALL PWRITE(INF,ME,NPC,W,    LBOUND(W,1),    UBOUND(W,1),    MElem,MElemGL,KB,EMAP,0,"W"    )
#     if defined (GOTM)
      CALL PWRITE(INF,ME,NPC,TKE,  LBOUND(TKE,1),  UBOUND(TKE,1),  NNode,MElemGL,KB,EMAP,0,"TKE"  )
      CALL PWRITE(INF,ME,NPC,TEPS, LBOUND(TEPS,1), UBOUND(TEPS,1), NNode,MElemGL,KB,EMAP,0,"TEPS" )
#     else
      CALL PWRITE(INF,ME,NPC,Q2,   LBOUND(Q2,1),   UBOUND(Q2,1),   NNode,MElemGL,KB,EMAP,0,"Q2"   )
      CALL PWRITE(INF,ME,NPC,Q2L,  LBOUND(Q2L,1),  UBOUND(Q2L,1),  NNode,MElemGL,KB,EMAP,0,"Q2L"  )
      CALL PWRITE(INF,ME,NPC,L,    LBOUND(L,1),    UBOUND(L,1),    NNode,MElemGL,KB,EMAP,0,"L"  )
#     endif
      CALL PWRITE(INF,ME,NPC,S,    LBOUND(S,1),    UBOUND(S,1),    MElem,MElemGL,KB,EMAP,0,"S"    )
      CALL PWRITE(INF,ME,NPC,T,    LBOUND(T,1),    UBOUND(T,1),    MElem,MElemGL,KB,EMAP,0,"T"    )
      CALL PWRITE(INF,ME,NPC,RHO,  LBOUND(RHO,1),  UBOUND(RHO,1),  MElem,MElemGL,KB,EMAP,0,"RHO"  )
      CALL PWRITE(INF,ME,NPC,TMEAN,LBOUND(TMEAN,1),UBOUND(TMEAN,1),MElem,MElemGL,KB,EMAP,0,"TMEAN")
      CALL PWRITE(INF,ME,NPC,SMEAN,LBOUND(SMEAN,1),UBOUND(SMEAN,1),MElem,MElemGL,KB,EMAP,0,"SMEAN")
      CALL PWRITE(INF,ME,NPC,RMEAN,LBOUND(RMEAN,1),UBOUND(RMEAN,1),MElem,MElemGL,KB,EMAP,0,"RMEAN")

      CALL PWRITE(INF,ME,NPC,S1,    LBOUND(S1,1),    UBOUND(S1,1),    NNode,NNodeGL,KB,NMAP,1,"S1"    )
      CALL PWRITE(INF,ME,NPC,T1,    LBOUND(T1,1),    UBOUND(T1,1),    NNode,NNodeGL,KB,NMAP,1,"T1"    )
      CALL PWRITE(INF,ME,NPC,RHO1,  LBOUND(RHO1,1),  UBOUND(RHO1,1),  NNode,NNodeGL,KB,NMAP,1,"RHO1"  )
      CALL PWRITE(INF,ME,NPC,TMEAN1,LBOUND(TMEAN1,1),UBOUND(TMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"TMEAN1")
      CALL PWRITE(INF,ME,NPC,SMEAN1,LBOUND(SMEAN1,1),UBOUND(SMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"SMEAN1")
      CALL PWRITE(INF,ME,NPC,RMEAN1,LBOUND(RMEAN1,1),UBOUND(RMEAN1,1),NNode,NNodeGL,KB,NMAP,1,"RMEAN1")

      CALL PWRITE(INF,ME,NPC,KM,LBOUND(KM,1),UBOUND(KM,1),NNode,MElemGL,KB,EMAP,1,"KM")
      CALL PWRITE(INF,ME,NPC,KH,LBOUND(KH,1),UBOUND(KH,1),NNode,MElemGL,KB,EMAP,1,"KH")
      CALL PWRITE(INF,ME,NPC,KQ,LBOUND(KQ,1),UBOUND(KQ,1),NNode,MElemGL,KB,EMAP,1,"KQ")

      CALL PWRITE(INF,ME,NPC,UA,LBOUND(UA,1),UBOUND(UA,1),MElem,MElemGL,1,EMAP,0,"UA")
      CALL PWRITE(INF,ME,NPC,VA,LBOUND(VA,1),UBOUND(VA,1),MElem,MElemGL,1,EMAP,0,"VA")

      CALL PWRITE(INF,ME,NPC,EL1,LBOUND(EL1,1),UBOUND(EL1,1),MElem,MElemGL,1,EMAP,1,"EL1")
      CALL PWRITE(INF,ME,NPC,ET1,LBOUND(ET1,1),UBOUND(ET1,1),MElem,MElemGL,1,EMAP,1,"ET1")
      CALL PWRITE(INF,ME,NPC,H1, LBOUND(H1,1), UBOUND(H1,1), MElem,MElemGL,1,EMAP,1,"H1" )
      CALL PWRITE(INF,ME,NPC,D1, LBOUND(D1,1), UBOUND(D1,1), MElem,MElemGL,1,EMAP,1,"D1" )
      CALL PWRITE(INF,ME,NPC,DT1,LBOUND(DT1,1),UBOUND(DT1,1),MElem,MElemGL,1,EMAP,1,"DT1")
      CALL PWRITE(INF,ME,NPC,RTP,LBOUND(RTP,1),UBOUND(RTP,1),MElem,MElemGL,1,EMAP,1,"RTP")

      CALL PWRITE(INF,ME,NPC,EL,LBOUND(EL,1),UBOUND(EL,1),NNode,NNodeGL,1,NMAP,1,"EL")
      CALL PWRITE(INF,ME,NPC,ET,LBOUND(ET,1),UBOUND(ET,1),NNode,NNodeGL,1,NMAP,1,"ET")
      CALL PWRITE(INF,ME,NPC,H, LBOUND(H,1), UBOUND(H,1), NNode,NNodeGL,1,NMAP,1,"H" )
      CALL PWRITE(INF,ME,NPC,D, LBOUND(D,1), UBOUND(D,1), NNode,NNodeGL,1,NMAP,1,"D" )
      CALL PWRITE(INF,ME,NPC,DT,LBOUND(DT,1),UBOUND(DT,1),NNode,NNodeGL,1,NMAP,1,"DT")

#     if defined (WATER_QUALITY)
      DO N1=1,NB
        CALL PWRITE(INF,ME,NPC,WQM(1:NNode,1:KB,N1),LBOUND(WQM(1:NNode,1:KB,N1),1),     &
                      UBOUND(WQM(1:NNode,1:KB,N1),1),NNode,NNodeGL,1,NMAP,1,"WQM")
      END DO
#     endif

#     endif
   END IF
   IF(MSR) CLOSE(INF)
   DEALLOCATE(RTP)

  RETURN
  END SUBROUTINE WRITERESTART

!=====================================================================================/
!  CONVERT THE STATE FROM THE GRID SPACE TO THE STATE VECTOR SPACE                    /
!=====================================================================================/
  SUBROUTINE GR2ST(OPT)
 
   USE LIMS
   USE ALL_VARS
   USE RRKVAL
   IMPLICIT NONE
   
    INTEGER IDUMMY, I, J
    INTEGER OPT

    IDUMMY=0
    
    IF(EL_ASSIM) THEN
      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        IF (OPT == 0) THEN
          STTEMP0(IDUMMY) = RRKEL(I)
        ELSE
          STTEMP1(IDUMMY) = RRKEL(I)
        ENDIF	 	
      ENDDO
    ENDIF

    IF(UV_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1
          IF (OPT == 0) THEN
	    STTEMP0(IDUMMY) = RRKU(J,I)
	  ELSE
	    STTEMP1(IDUMMY) = RRKU(J,I)
	  ENDIF    
        ENDDO
      ENDDO
      DO I=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1
          IF (OPT == 0) THEN
	    STTEMP0(IDUMMY) = RRKV(J,I)
	  ELSE
	    STTEMP1(IDUMMY) = RRKV(J,I)
	  ENDIF    
        ENDDO
      ENDDO		
    ENDIF  
  
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1
          IF (OPT == 0) THEN
	    STTEMP0(IDUMMY) = RRKT(J,I)
	  ELSE
	    STTEMP1(IDUMMY) = RRKT(J,I)
	  ENDIF    
        ENDDO
      ENDDO
    ENDIF
  
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1
          IF (OPT == 0) THEN
	    STTEMP0(IDUMMY) = RRKS(J,I)
	  ELSE
	    STTEMP1(IDUMMY) = RRKS(J,I)
	  ENDIF    
        ENDDO
      ENDDO
    ENDIF

  RETURN
  END SUBROUTINE GR2ST

!=====================================================================================/
!  CONVERT THE STATE FROM THE STATE VECTOR SPACE TO THE GRID SPACE                    /
!=====================================================================================/
  SUBROUTINE ST2GR
 
   USE LIMS
   USE ALL_VARS
   USE RRKVAL
   IMPLICIT NONE
    
    INTEGER IDUMMY, I, J

    IDUMMY=0
    IF(EL_ASSIM) THEN
      DO I=1, NNodeGL
        IDUMMY = IDUMMY + 1
        RRKEL(I) = STTEMP1(IDUMMY) 
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1
          RRKU(J,I) = STTEMP1(IDUMMY) 
        ENDDO
      ENDDO
      DO I=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1
          RRKV(J,I) = STTEMP1(IDUMMY)
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1
          RRKT(J,I) = STTEMP1(IDUMMY) 
        ENDDO
      ENDDO
    ENDIF
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1
          RRKS(J,I) = STTEMP1(IDUMMY) 
        ENDDO
      ENDDO
    ENDIF
    
  RETURN
  END SUBROUTINE ST2GR

!=====================================================================================/
!  READ THE EOFs FROM eof.cdf                                                         /
!=====================================================================================/
  SUBROUTINE READEOF(IEOF,OPT)
 
   USE LIMS
   USE RRKVAL
   IMPLICIT NONE

#include "/hosts/salmon01/data00/medm/src/netcdf-3.6.0-p1/src/fortran/netcdf.inc"    
!#include "/usr/local/include/nedcdf.inc" 
    INTEGER OPT
    INTEGER I,J,K
    INTEGER RCODE
    INTEGER IEOF
    INTEGER START(3)
    INTEGER COUNT(3)
    INTEGER STATUS
    INTEGER NCID
    INTEGER VARID
    INTEGER IDUMMY
    CHARACTER(LEN=80) FNAME
    REAL(DP),ALLOCATABLE :: UTMP(:,:)
    REAL(DP),ALLOCATABLE :: VTMP(:,:)
    REAL(DP),ALLOCATABLE :: ELTMP(:)
    REAL(DP),ALLOCATABLE :: TTMP(:,:)
    REAL(DP),ALLOCATABLE :: STMP(:,:)

    ALLOCATE(UTMP(MElemGL,KBM1))         ; UTMP  = ZERO
    ALLOCATE(VTMP(MElemGL,KBM1))         ; VTMP  = ZERO
    ALLOCATE(ELTMP(NNodeGL))             ; ELTMP = ZERO
    ALLOCATE(TTMP(NNodeGL,KBM1))         ; TTMP  = ZERO
    ALLOCATE(STMP(NNodeGL,KBM1))         ; STMP  = ZERO    
    
    
    FNAME=TRIM(OUTDIR)//'/rrktemp/'//'eof.cdf'
    STATUS = nf_open(FNAME,NF_NOWRITE,NCID)
    
    IDUMMY = 0

    IF(EL_ASSIM) THEN
      STATUS = nf_inq_varid(NCID,'eof_el',VARID)
    
      START(1) = 1
      COUNT(1) = NNodeGL
      START(2) = IEOF
      COUNT(2) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,ELTMP)
    
      DO J=1, NNodeGL
        IDUMMY = IDUMMY +1
        STEOF(IDUMMY) = ELTMP(J)
      ENDDO
    ENDIF
 
    IF(UV_ASSIM) THEN
      STATUS = nf_inq_varid(NCID,'eof_u',VARID)
     
      START(1) = 1
      COUNT(1) = MElemGL
      START(2) = 1
      COUNT(2) = KBM1
      START(3) = IEOF
      COUNT(3) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,UTMP)
    
      DO K=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1 
          STEOF(IDUMMY) = UTMP(J,K)
        ENDDO
      ENDDO

      STATUS = nf_inq_varid(NCID,'eof_v',VARID)
     
      START(1) = 1
      COUNT(1) = MElemGL
      START(2) = 1
      COUNT(2) = KBM1
      START(3) = IEOF
      COUNT(3) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,VTMP)
    
      DO K=1, KBM1
        DO J=1, MElemGL
          IDUMMY = IDUMMY + 1 
          STEOF(IDUMMY) = VTMP(J,K)
        ENDDO
      ENDDO
    ENDIF
   
! READ THE EOF FOR THE TEMPERATURE 
    IF(T_ASSIM) THEN
      STATUS = nf_inq_varid(NCID,'eof_temp',VARID)
     
      START(1) = 1
      COUNT(1) = NNodeGL
      START(2) = 1
      COUNT(2) = KBM1
      START(3) = IEOF
      COUNT(3) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,TTMP)
    
      DO K=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1 
          STEOF(IDUMMY) = TTMP(J,K)
        ENDDO
      ENDDO   
    ENDIF     

! READ THE EOF FOR THE SALINITY 
    IF(S_ASSIM) THEN
      STATUS = nf_inq_varid(NCID,'eof_sal',VARID)
     
      START(1) = 1
      COUNT(1) = NNodeGL
      START(2) = 1
      COUNT(2) = KBM1
      START(3) = IEOF
      COUNT(3) = 1
      RCODE = nf_get_vara_double(NCID,VARID,START,COUNT,STMP)
    
      DO K=1, KBM1
        DO J=1, NNodeGL
          IDUMMY = IDUMMY + 1 
          STEOF(IDUMMY) = STMP(J,K)
        ENDDO
      ENDDO   
    ENDIF
   
    RCODE = nf_close(NCID)

    IF(OPT == 2) STTEMP1=STEOF     

    DEALLOCATE(UTMP,VTMP,ELTMP,TTMP,STMP)
    
    RETURN
  END SUBROUTINE READEOF

!=====================================================================================/
!  DETERMINE IF NODES/ELEMENTS ARE WET OR DRY                                         /
!=====================================================================================/
  SUBROUTINE WET_JUDGE_EL

   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined(WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   REAL(SP) :: DTMP
   INTEGER  :: ITA_TEMP
   INTEGER  :: I,IL,IA,IB,K1,K2,K3,K4,K5,K6

#  if defined(WET_DRY)

!
!--Determine If Node Points Are Wet/Dry Based on Depth Threshold---------------!
!
   ISWET_NODE_CURRENTSTEP = 1
   DO I = 1, NNode
     DTMP = H(I) + EL(I)
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWET_NODE_CURRENTSTEP(I) = 0
   END DO

!
!--Determine if Cells are Wet/Dry Based on Depth Threshold---------------------!
!
   ISWET_CELL_CURRENTSTEP = 1
   DO I = 1, MElem
     DTMP =  MAX(EL(NV(I,1)),EL(NV(I,2)),EL(NV(I,3)))  + &
             MIN(  H(NV(I,1)),  H(NV(I,2)),  H(NV(I,3)))
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWET_CELL_CURRENTSTEP(I) = 0
   END DO

!
!--A Secondary Condition for Nodal Dryness-(All Elements Around Node Are Dry)--!
!
   DO I = 1, NNode
     IF(SUM(ISWET_CELL_CURRENTSTEP(NBVE(I,1:NTVE(I)))) == 0)  ISWET_NODE_CURRENTSTEP(I) = 0
   END DO

!
!--Adjust Water Surface So It Does Not Go Below Minimum Depth------------------!
!
   EL = MAX(EL,-H + MIN_DEPTH)

!
!--Recompute Element Based Depths----------------------------------------------!
!
   DO I = 1, MElem
     EL1(I) = ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))
   END DO

!
!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     FWET_N_N = FLOAT(ISWET_NODE_CURRENTSTEP)
     FWET_C_C = FLOAT(ISWET_CELL_CURRENTSTEP)
     CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,FWET_C_C)
     CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,FWET_N_N)
     ISWET_NODE_CURRENTSTEP = INT(FWET_N_N+.5)
     ISWET_CELL_CURRENTSTEP = INT(FWET_C_C+.5)
   END IF
#  endif

#  endif

   RETURN

  END SUBROUTINE WET_JUDGE_EL

!==============================================================================|
!   DUMP WET/DRY FLAG DATA FOR RESTART                                         |
!==============================================================================|

   SUBROUTINE WD_DUMP_EL(INF,I_START)

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined(WET_DRY)
   USE MOD_WD
#  endif

   IMPLICIT NONE
   INTEGER, ALLOCATABLE,DIMENSION(:) :: NTEMP1,NTEMP2
   INTEGER I, INF
   INTEGER I_START
!==============================================================================|

#  if defined(WET_DRY)

   IF(MSR)THEN
     REWIND(INF)
     WRITE(INF,*) I_START
     WRITE(INF,*) MElemGL,NNodeGL
   END IF

   IF(SERIAL)THEN
     WRITE(INF,*) (ISWET_CELL_CURRENTSTEP(I), I=1,MElem)
     WRITE(INF,*) (ISWET_NODE_CURRENTSTEP(I), I=1,NNode)
   ELSE
   ALLOCATE(NTEMP1(MElemGL),NTEMP2(NNodeGL))
#  if defined (MULTIPROCESSOR)
   CALL IGATHER(LBOUND(ISWET_CELL_CURRENTSTEP,1),UBOUND(ISWET_CELL_CURRENTSTEP,1),MElem,MElemGL,1,MYID,NPROCS,EMAP,ISWET_CELL_CURRENTSTEP,NTEMP1)
   CALL IGATHER(LBOUND(ISWET_NODE_CURRENTSTEP,1),UBOUND(ISWET_NODE_CURRENTSTEP,1),NNode,NNodeGL,1,MYID,NPROCS,NMAP,ISWET_NODE_CURRENTSTEP,NTEMP2)
   IF(MSR)THEN
     WRITE(INF,*) (NTEMP1(I), I=1,MElemGL)
     WRITE(INF,*) (NTEMP2(I), I=1,NNodeGL)
   END IF
   DEALLOCATE(NTEMP1,NTEMP2)
#  endif
   END IF

   CLOSE(INF)

#  endif

   RETURN
   END SUBROUTINE WD_DUMP_EL

  SUBROUTINE READOBS(STLOC,NLOC)
    
   USE LIMS
   USE CONTROL
   IMPLICIT NONE
   
     INTEGER ::  NUM  = 0 
     INTEGER ::  NLOC
     INTEGER ::  SWITCH = 0
     INTEGER ::  J,K
     INTEGER ::  IDUMMY = 0
     INTEGER ::  TMP
     CHARACTER(LEN=80) FILENAME
     CHARACTER(LEN=24) HEADINFO
     INTEGER STLOC(RRK_NOBSMAX)
     INTEGER LAY(RRK_NOBSMAX)

     FILENAME = TRIM(INPDIR)//"/"//trim(casename)//"_assim_rrkf.dat"
    
     OPEN(73,FILE=TRIM(FILENAME),FORM='FORMATTED')

     NLOC = 0
     
 100 READ(73,'(A24)',END=200) HEADINFO
     IF(SWITCH/=1) THEN
       IF(HEADINFO=='!=== READ IN OBSERVATION') THEN
         SWITCH = 1
         GOTO 100
       ELSE
         GOTO 100
       ENDIF
     ENDIF 
     
     IF(TRIM(HEADINFO)=='!EL') THEN
       IF(EL_OBS) THEN
         READ(73,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(73,*)  (STLOC(K), K=1,NLOC)	
       ENDIF  
       
       IF(EL_ASSIM) THEN
	 IDUMMY = IDUMMY + NNodeGL
       ENDIF
            
     ENDIF
     
     IF(TRIM(HEADINFO)=='!UV') THEN
       IF(UV_OBS) THEN       
         READ(73,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC+NUM>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc+num, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(73,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)
	 READ(73,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MElemGL*(LAY(K)-1)
	 ENDDO   
         IDUMMY = IDUMMY + MElemGL*KBM1
	 
	 NLOC = NLOC + NUM
	 DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K-NUM)+MElemGL*KBM1+MElemGL*(LAY(K-NUM)-1)	   
	 ENDDO
       ENDIF
       
       IF(UV_ASSIM) THEN  
          IDUMMY = IDUMMY + MElemGL*KBM1
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!T') THEN
       IF(T_OBS) THEN
         READ(73,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(73,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)       
         READ(73,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)       
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NNodeGL*(LAY(K)-1)
	 ENDDO   
         
       ENDIF   
       
       IF(T_ASSIM) THEN
         IDUMMY = IDUMMY + NNodeGL*KBM1
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!S') THEN
       IF(S_OBS) THEN
         READ(73,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(73,*)  (STLOC(K),K=NLOC-NUM+1,NLOC)        
         READ(73,*)  (LAY(K),  K=NLOC-NUM+1,NLOC)      
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NNodeGL*(LAY(K)-1)
	 ENDDO   
 
       ENDIF
       
       IF(S_ASSIM) THEN
         IDUMMY = IDUMMY + NNodeGL*KBM1 
       ENDIF 
     ENDIF
     
     GOTO 100
 200 CONTINUE

     DO J=1, NLOC-1
       DO K=2, NLOC
	 IF(STLOC(K)<STLOC(J)) THEN
	   TMP = STLOC(J)
	   STLOC(J) = STLOC(K)
	   STLOC(K) = TMP 
         ENDIF
       ENDDO
     ENDDO

     CLOSE(73)
    
  RETURN 
  END SUBROUTINE READOBS  

# endif
END MODULE MOD_RRK
 
