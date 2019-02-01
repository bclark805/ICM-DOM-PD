!==============================================================================|
!   CONTROL VARIABLES                                                          |
!==============================================================================|

MODULE MOD_HEATFLUX
# if defined (HEAT_FLUX)
   USE CONTROL, ONLY : DEG2RAD
   USE MOD_TYPES
   USE MOD_PREC
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   IMPLICIT NONE
   SAVE

!--File Unit Specifiers -------------------------------------------------------!
   INTEGER  INCWH_AIR,INHFX_AIR
   
!----------------surface forcing: nouniform meteo conditions-------------------!

   TYPE(BC)               :: HFX1_TM        !!TIME MAP FOR TIME DEPENDENT OF CALCU HF DATA
   REAL(SP), ALLOCATABLE  :: T_AIRN(:,:)    !!BULK AIR TEMPERATURE
   REAL(SP), ALLOCATABLE  :: RH_AIRN(:,:)   !!RELATIVE HUMIDITY
   REAL(SP), ALLOCATABLE  :: PA_AIRN(:,:)   !!SURFACE PRESSURE
   REAL(SP), ALLOCATABLE  :: DLW_AIRN(:,:)  !!DOWNWARD LONGWAVE RADIATION
   REAL(SP), ALLOCATABLE  :: DSW_AIRN(:,:)  !!DOWNWARD SHORTWAVE RADIATION
   REAL(SP), ALLOCATABLE  :: EL_AIR(:)
   REAL(SP), ALLOCATABLE  :: EGF_AIR(:)
   REAL(SP), ALLOCATABLE  :: ELF_AIR(:)
   REAL(SP), ALLOCATABLE  :: ELRK_AIR(:)

   REAL(SP), ALLOCATABLE  :: QA_AIRN(:,:)        !!SPECIFIC HUMIDITY
   REAL(SP), ALLOCATABLE  :: CLOUDN(:,:)         !!Total Cloud Cover
   
!----------------boundary conditions: uniform meteo conditions-----------------!

   TYPE(BC)               :: UMF1_TM           !!TIME MAPPING FOR UNIFORM METEOS
   REAL(SP), ALLOCATABLE  :: T_AIRU(:)         !!BULK AIR TEMPERATURE
   REAL(SP), ALLOCATABLE  :: RH_AIRU(:)        !!RELATIVE HUMIDITY
   REAL(SP), ALLOCATABLE  :: PA_AIRU(:)        !!SURFACE PRESSURE
   REAL(SP), ALLOCATABLE  :: DLW_AIRU(:)       !!DOWNWARD LONGWAVE RADIATION
   REAL(SP), ALLOCATABLE  :: DSW_AIRU(:)       !!DOWNWARD SHORTWAVE RADIATION

   REAL(SP) :: ZU                             !!HEIGHT OF WIND SPEED (M) 
   REAL(SP) :: ZT                             !!HEIGHT OF AIR TEMPERATURE (M)
   REAL(SP) :: ZQ                             !!HEIGHT OF RELATIVE HUMIDITY (M)
   
   LOGICAL  :: C_HFX                          !!FLAG TO CONTROL HEAT FLUX CALCULATION
   
   REAL(SP), ALLOCATABLE :: CORRG(:),CORR(:)           !!LATITUDE OF NODES
!==============================================================================|


   CONTAINS !-------------------------------------------------------------------|
            ! DATA_RUN_HFX : Input Parameters Which Control the Calculation of  |
	    !                Heat Flux                                          |
            ! IOFILES_HFX  : Open Input Files for Heat Flux Calculation         !
            ! BCS_FORCE_HFX: Set Up the Surface BCs for Calculating Heat Flux   |
	    ! BCOND_HFX    : Interprates the Surface BCs for Calculating Heat   |
	    !                Flux                                               |
	    ! COARE26Z     : Calculates Sensible Heat Flux and Latent Heat Flux |
            ! PSIT_26      : Computes Temperature Structure Function            |
	    ! PSIU_26      : Computes Velocity Structure Function               |
	    ! QSAT26       : Computes Saturation Specific Humidity              |
	    ! GRV          : Computes g Given Lat in Deg                        |
	    
!==============================================================================|
!   Input Parameters Which Control the Calculation of Heat Flux                |
!==============================================================================|

   SUBROUTINE DATA_RUN_HFX            

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_INP
   IMPLICIT NONE
   INTEGER  INTVEC(150),ISCAN
   CHARACTER(LEN=120) :: FNAME


!==============================================================================|
!   READ IN VARIABLES AND SET VALUES                                           |
!==============================================================================|

   FNAME = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------|
!  HEAT FLUX CALCULATION FLAG
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"C_HFX",LVAL = C_HFX)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING C_HFX: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "C_HFX" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
                                                                                                                            
!------------------------------------------------------------------------------|
!     "ZU"   !!HEIGHT OF WIND SPEED
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ZU",FSCAL = ZU)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ZU: ',ISCAN
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!     "ZT"   !!HEIGHT OF AIR TEMPERATURE
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ZT",FSCAL = ZT)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ZT: ',ISCAN
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!     "ZQ"   !!HEIGHT OF RELATIVE HUMIDITY
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ZQ",FSCAL = ZQ)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ZQ: ',ISCAN
     CALL PSTOP
   END IF

!==============================================================================|
!            SCREEN REPORT OF SET VARIABLES                                    !
!==============================================================================|
   IF(MSR)THEN
   
   IF(C_HFX)THEN
     WRITE(IPT,*)'!  # CALCULATE HEAT FLUX :  ACTIVE'
     WRITE(IPT,*)'!  # ZU (m)              :',ZU   
     WRITE(IPT,*)'!  # ZT (m)              :',ZT   
     WRITE(IPT,*)'!  # ZQ (m)              :',ZQ   
   END IF  

   END IF

   RETURN
   END SUBROUTINE DATA_RUN_HFX    
!------------------------------------------------------------------------------|


!==============================================================================!
!   Open Input Files for Heat Flux Calculation                                 !
!==============================================================================!

   SUBROUTINE IOFILES_HFX
   USE ALL_VARS
   USE MOD_UTILS
   IMPLICIT NONE

   CHARACTER(LEN=80)  :: ISTR

!==============================================================================!
!                Definitions of input files                                    !
!                                                                              !
! incwh_air: casename_mc_air.dat : boundary input meteorlogical forcing:       !
!                                  bulk air temperature at height zt,          !
!                                  relative humidity at height zq,             !
!                                  surface pressure,                           !
!                                  downward longwave radiation,                !
!                                  downward shortwave radiation.               !
! inhfx_air: casename_hfx_air.dat: real-time field of atomospheric input. The  !
!                                  variables are same as above.                !
!==============================================================================!
         
   INCWH_AIR =34
   INHFX_AIR =35
   ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)
!
!-----------------OPEN METEOROLOGICAL FORCING FILES----------------------------!
!
   IF(M_TYPE == 'uniform')THEN
     CALL FOPEN(INCWH_AIR, TRIM(ISTR)//'_mc_air.dat',"cfr")
   ELSE 
     CALL FOPEN(INHFX_AIR, TRIM(ISTR)//'_hfx_air.dat' ,"cur")
   END IF

   RETURN
   END SUBROUTINE IOFILES_HFX
!==============================================================================!

!==============================================================================|
!   Set Up the Following Boundary Conditions for Calculating Heat Flux:        |
!     Meteorological Forcing 						       |
!==============================================================================|

   SUBROUTINE BCS_FORCE_HFX           

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE BCS
   USE MOD_CLOCK
   USE MOD_UTILS

   IMPLICIT NONE
   CHARACTER(LEN=80) :: COMT
   REAL(SP) :: TTIME
   REAL(SP) :: TEMP
   REAL(SP) :: FTEMP1,FTEMP2,FTEMP3,FTEMP4,FTEMP5
   REAL(SP) :: FTEMP6,FTEMP7
   REAL(SP) :: RBUF1,RBUF2,RBUF3,RBUF4,RBUF5
   REAL(SP) :: RBUF6,RBUF7
   REAL(SP), ALLOCATABLE :: RTEMP(:),RTEMP1(:,:),RTEMP2(:,:),RTEMP3(:,:)
   REAL(SP), ALLOCATABLE :: RTEMP4(:,:),RTEMP5(:,:)
   REAL(SP), ALLOCATABLE :: RTEMP6(:,:),RTEMP7(:,:)
   INTEGER   I,J,IOS,NCNT,IERR
   CHARACTER(LEN=13) :: TSTRING

   CHARACTER(LEN=80)  :: ISTR
   ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)
!------------------------------------------------------------------------------|
!--------------READ IN LATITUDE------------------------------------------------!

   ALLOCATE(CORRG(0:NNodeGL))  ; CORRG = 0.0_SP
!#  if defined (SPHERICAL)
!   CORRG=YG
!#  else
   IF (MSR) THEN
   CALL FOPEN(INCOR, TRIM(ISTR)//'_cor.dat',"cfr")
   REWIND(INCOR)
   DO I=1,NNodeGL
     READ(INCOR,*) TEMP,TEMP,CORRG(I)
   END DO
   REWIND(INCOR)
   ENDIF
!#  endif

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

!----------------------------REPORT--------------------------------------------!
   IF(MSR)WRITE(IPT,*)'!'
   IF(MSR)WRITE(IPT,*)'!            SETTING UP PRESCRIBED BOUNDARY CONDITIONS  '
   IF(MSR)WRITE(IPT,*)'!Meteorological Forcing Info for Calculating Heat Flux  '
   IF(MSR)WRITE(IPT,*)'!'

!==============================================================================|
!   Input Meteorological Boundary Conditions for Calculating Heat Flux         |
!==============================================================================|
!    bulk air temperature at height 2m:   degree(C)    "t_air"                 |
!    relative humidity at height 2m:      (%)          "rh_air"                |
!    surface pressure:                    mb           "pa_air"                |
!    downward longwave radiation:         w/m^2        "dlw_air"               |
!    downward shortwave radiation:        w/m^2        "dsw_air"               |
!==============================================================================|

   IF(M_TYPE == 'uniform')THEN

!==============================================================================|
!   UNIFORM METEOLOGICAL CONDITIONS                                            |
!==============================================================================|

     READ(INCWH_AIR,1000) COMT
     IF(MSR)WRITE(IOPRT,*)'Meteorological Forcing Info for Calculating Heat Flux'
     IF(MSR)WRITE(IOPRT,1000) COMT

!
!----Determine Number of Data Times--------------------------------------------!
!
     NCNT = 0
     DO WHILE(.TRUE.)
       READ(INCWH_AIR,*,END=15,IOSTAT=IOS)
       READ(INCWH_AIR,*,END=15,IOSTAT=IOS)
       IF(IOS < 0)EXIT
       NCNT = NCNT + 1
     END DO
 15  CONTINUE
     IF(NCNT == 0)CALL PERROR(6,"NO UNIFORM METEO DATA FOR CALCULATING &
                  HEAT FLUX PROVIDED")

     REWIND(INCWH_AIR) ; READ(INCWH_AIR,*)

!
!----Read in Bulk Air Temperature/Relative Humidity/Surface Pressure/
!    Downward Longwave Radiation/Downward Shortwave Radiation Data at Each Time-!
!

     UMF1_TM%NTIMES = NCNT
     ALLOCATE(UMF1_TM%TIMES(NCNT))
     ALLOCATE(T_AIRU(NCNT),RH_AIRU(NCNT))
     ALLOCATE(PA_AIRU(NCNT))
     ALLOCATE(DLW_AIRU(NCNT),DSW_AIRU(NCNT))

     DO I=1,NCNT
       READ(INCWH_AIR ,*) TTIME
       IF(MSR)WRITE(IOPRT,*) TTIME
       UMF1_TM%TIMES(I) = TTIME

!                         air temperature, relative humidity, pressure, long wave radation downward (?)
!                                                                         short wave radiation (?)
       READ(INCWH_AIR ,*) T_AIRU(I), RH_AIRU(I), PA_AIRU(I), DLW_AIRU(I), DSW_AIRU(I)

       IF(MSR)WRITE(IOPRT,5000) T_AIRU(I), RH_AIRU(I), PA_AIRU(I), &
              DLW_AIRU(I), DSW_AIRU(I)
     END DO

     IF(WINDTYPE /= 'speed') THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'TO CALCULATE HEAT FLUX, WINDTYPE MUST BE "speed"'
       WRITE(IPT,*)'WINDTYPE IS NOT CORRECT, --->',WINDTYPE
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

     CLOSE(INCWH_AIR)

!
!--REPORT RESULTS--------------------------------------------------------------!
!

   IF(MSR)THEN
     WRITE(IPT,*)'!'
     WRITE(IPT,*    )'!  UNIFORM METEO FOR HEAT FLUX :    SET'
      IF(UMF1_TM%NTIMES > 0)THEN
        CALL GETTIME(TSTRING,INT(3600.*UMF1_TM%TIMES(1)))
        WRITE(IPT,102)'!  METEO DATA BEGIN            :  ',TSTRING
        CALL GETTIME(TSTRING,INT(3600.*UMF1_TM%TIMES(UMF1_TM%NTIMES)))
        WRITE(IPT,102)'!  METEO DATA END              :  ',TSTRING
      END IF
    END IF

!==============================================================================|
!   NON-UNIFORM METEOLOGICAL CONDITIONS                                        |
!==============================================================================|

   ELSE IF (M_TYPE == 'non-uniform')THEN

     REWIND(INHFX_AIR)
!
!----Input Number of Data Times for Meteological conditions to calculate-------| 
!----Heat Flux ----------------------------------------------------------------|
!
     NCNT = 0
     IF(MSR)THEN

     DO WHILE(.TRUE.)
       READ(INHFX_AIR,END=10)FTEMP1
       READ(INHFX_AIR)
       NCNT = NCNT + 1
     END DO
 10  CONTINUE
     REWIND(INHFX_AIR)

     IF(NCNT == 0)CALL PERROR(6,"NO DATA PROVIDED FOR CALCULATING HEAT FLUX")
     END IF

     HFX1_TM%NTIMES =  NCNT 

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(HFX1_TM%NTIMES,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#    endif

!
!----Read in Data Times and Global Meteological conditions for Heat Flux-------!
!

     ALLOCATE(HFX1_TM%TIMES(HFX1_TM%NTIMES))
     ALLOCATE(RTEMP1(NNodeGL,HFX1_TM%NTIMES),RTEMP2(NNodeGL,HFX1_TM%NTIMES))
     ALLOCATE(RTEMP3(NNodeGL,HFX1_TM%NTIMES),RTEMP4(NNodeGL,HFX1_TM%NTIMES))
     ALLOCATE(RTEMP5(NNodeGL,HFX1_TM%NTIMES))
!      for ice model
     ALLOCATE(RTEMP6(NNodeGL,HFX1_TM%NTIMES),RTEMP7(NNodeGL,HFX1_TM%NTIMES))

     IF(MSR)THEN
       DO J=1,HFX1_TM%NTIMES
         READ(INHFX_AIR) HFX1_TM%TIMES(J) 
         READ(INHFX_AIR) (RTEMP1(I,J),RTEMP2(I,J),RTEMP3(I,J),RTEMP4(I,J),  &
                      RTEMP5(I,J),RTEMP6(I,J),RTEMP7(I,J),I=1,NNodeGL)
!                      RTEMP5(I,J),I=1,NNodeGL)
       END DO
     END IF
!
!----Broadcast Data------------------------------------------------------------!
!


#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       CALL MPI_BCAST(HFX1_TM%TIMES,HFX1_TM%NTIMES,MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP1,HFX1_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP2,HFX1_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP3,HFX1_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP4,HFX1_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP5,HFX1_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP6,HFX1_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(RTEMP7,HFX1_TM%NTIMES*NNodeGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
     END IF
#    endif

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     ALLOCATE(T_AIRN(NNode,HFX1_TM%NTIMES),RH_AIRN(NNode,HFX1_TM%NTIMES))
     ALLOCATE(PA_AIRN(NNode,HFX1_TM%NTIMES))
     ALLOCATE(DLW_AIRN(NNode,HFX1_TM%NTIMES),DSW_AIRN(NNode,HFX1_TM%NTIMES))
!    for ice
     ALLOCATE(QA_AIRN(NNode,HFX1_TM%NTIMES),CLOUDN(NNode,HFX1_TM%NTIMES))
     IF(SERIAL)THEN
       T_AIRN(1:NNodeGL,:)   = RTEMP1(1:NNodeGL,:)
       RH_AIRN(1:NNodeGL,:)  = RTEMP2(1:NNodeGL,:)
       PA_AIRN(1:NNodeGL,:)  = RTEMP3(1:NNodeGL,:)
       DLW_AIRN(1:NNodeGL,:) = RTEMP4(1:NNodeGL,:)
       DSW_AIRN(1:NNodeGL,:) = RTEMP5(1:NNodeGL,:)
       QA_AIRN(1:NNodeGL,:) =  RTEMP6(1:NNodeGL,:)
       CLOUDN(1:NNodeGL,:) =   RTEMP7(1:NNodeGL,:)
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,NNode
         T_AIRN(I,:)   = RTEMP1(NGID(I),:) 
         RH_AIRN(I,:)  = RTEMP2(NGID(I),:) 
         PA_AIRN(I,:)  = RTEMP3(NGID(I),:) 
         DLW_AIRN(I,:) = RTEMP4(NGID(I),:) 
         DSW_AIRN(I,:) = RTEMP5(NGID(I),:) 
         QA_AIRN(I,:)  = RTEMP6(NGID(I),:) 
         CLOUDN(I,:)   = RTEMP7(NGID(I),:) 
       END DO
     END IF
#   endif

!     DEALLOCATE(RTEMP1,RTEMP2,RTEMP3,RTEMP4,RTEMP5)
     DEALLOCATE(RTEMP1,RTEMP2,RTEMP3,RTEMP4,RTEMP5,RTEMP6,RTEMP7)
!     IF(MSR)WRITE(IPT,101)'!  METEO CONDITIONS FOR HFLUX     :    COMPLETE'
      
!
!--REPORT RESULTS--------------------------------------------------------------!
!

!     IF(MSR)WRITE(IPT,*)'!'
     IF(MSR)WRITE(IPT,*    )'!  NON-UNIFORM METEO FOR HFLUX     :    SET'
     IF(HFX1_TM%NTIMES > 0)THEN
       CALL GETTIME(TSTRING,3600*INT(HFX1_TM%TIMES(1)))
       IF(MSR)WRITE(IPT,102)'!  DATA BEGIN   :  ',TSTRING        
       CALL GETTIME(TSTRING,3600*INT(HFX1_TM%TIMES(HFX1_TM%NTIMES)))
       IF(MSR)WRITE(IPT,102)'!  DATA END     :  ',TSTRING
     END IF

     FTEMP1 = SUM(T_AIRN)/FLOAT(NNode*HFX1_TM%NTIMES)
     FTEMP2 = SUM(RH_AIRN)/FLOAT(NNode*HFX1_TM%NTIMES)
     FTEMP3 = SUM(PA_AIRN)/FLOAT(NNode*HFX1_TM%NTIMES)
     FTEMP4 = SUM(DLW_AIRN)/FLOAT(NNode*HFX1_TM%NTIMES)
     FTEMP5 = SUM(DSW_AIRN)/FLOAT(NNode*HFX1_TM%NTIMES)
     FTEMP6 = SUM(QA_AIRN)/FLOAT(NNode*HFX1_TM%NTIMES)
     FTEMP7 = SUM(CLOUDN)/FLOAT(NNode*HFX1_TM%NTIMES)

     IF(SERIAL)THEN
       RBUF1 = FTEMP1
       RBUF2 = FTEMP2
       RBUF3 = FTEMP3
       RBUF4 = FTEMP4
       RBUF5 = FTEMP5
       RBUF6 = FTEMP6
       RBUF7 = FTEMP7
     END IF
#  if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_REDUCE(FTEMP1,RBUF1,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_REDUCE(FTEMP2,RBUF2,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_REDUCE(FTEMP3,RBUF3,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_REDUCE(FTEMP4,RBUF4,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_REDUCE(FTEMP5,RBUF5,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
!    ggao
     IF(PAR)CALL MPI_REDUCE(FTEMP6,RBUF6,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
     IF(PAR)CALL MPI_REDUCE(FTEMP7,RBUF7,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
#  endif
     IF(MSR)WRITE(IPT,103)'!  AVE BULK AIR TEMPERATURE   :',RBUF1/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,103)'!  AVE RELATIVE HUMIDITY      :',RBUF2/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,103)'!  AVE SURFACE PRESSURE       :',RBUF3/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,103)'!  AVE DOWNWARD LONGWAVE RAD  :',RBUF4/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,103)'!  AVE DOWNWARD SHORTWAVE RAD :',RBUF5/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,103)'!  AVE SPECIFIC HUMIDITY      :',RBUF6/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,103)'!  AVE CLOUD COVER            :',RBUF7/FLOAT(NPROCS)

   ELSE
     WRITE(IPT,*)'==================ERROR=================================='
     WRITE(IPT,*)'M_TYPE NOT CORRECT, --->',M_TYPE
     WRITE(IPT,*)'MUST BE "uniform" or "non-uniform"'
     WRITE(IPT,*)'========================================================='
     CALL PSTOP
   END IF
    
   ALLOCATE(EL_AIR(0:NTNode));   EL_AIR = 0.0
   ALLOCATE(EGF_AIR(0:NTNode));   EGF_AIR = 0.0
   ALLOCATE(ELF_AIR(0:NTNode));   ELF_AIR = 0.0
   ALLOCATE(ELRK_AIR(0:NTNode));   ELRK_AIR = 0.0
   
!
!--Format Statements-----------------------------------------------------------!
!

   101  FORMAT(1X,A26,F10.4)  
   102  FORMAT(1X,A28,A13)  
   103  FORMAT(1X,A26,F10.1)  
   1000 FORMAT(A80)
   5000 FORMAT(8E14.5)

   RETURN
   END SUBROUTINE BCS_FORCE_HFX
!==============================================================================|

!==============================================================================|
!     SET BOUNDARY CONDITIONS FOR CALCULATING HEAT FLUX                        |
!                                                                              |
!==============================================================================|

   SUBROUTINE BCOND_HFX

!==============================================================================|
   USE ALL_VARS
   USE BCS


#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   REAL(SP) :: FACT,UFACT,TX,TY
   REAL(SP) :: WDS,T_AIRF,RH_AIRF,PA_AIRF,DLW_AIRF,DSW_AIRF
!    for ice model
   REAL(SP) :: QA_AIRF,CLOUDF
   REAL(SP) :: HSB,HLB,TAU
   REAL(SP), DIMENSION(0:MTElem) :: WDSN
   INTEGER  I,J,L1,L2,IERR

!==============================================================================|
!   SURFACE FORCING FOR INTERNAL MODE                                          !
!==============================================================================|
   IF(M_TYPE == 'uniform') THEN

     CALL BRACKET(UMF1_TM,THOUR,L1,L2,FACT,UFACT,IERR)
!
!--- Heat Flux for the Internal Mode-----------------------------------------!
!
     IF(WINDTYPE == 'speed')THEN
       TX = UFACT*UWIND(L1) + FACT*UWIND(L2) 
       TY = UFACT*VWIND(L1) + FACT*VWIND(L2) 
       WDS=SQRT(TX*TX+TY*TY)
       T_AIRF  = UFACT*T_AIRU(L1) + FACT*T_AIRU(L2)
       RH_AIRF = UFACT*RH_AIRU(L1) + FACT*RH_AIRU(L2)
       PA_AIRF = UFACT*PA_AIRU(L1) + FACT*PA_AIRU(L2)
       DLW_AIRF = UFACT*DLW_AIRU(L1) + FACT*DLW_AIRU(L2)
       DSW_AIRF = UFACT*DSW_AIRU(L1) + FACT*DSW_AIRU(L2)
       DO I=1,NNode
         CALL COARE26Z(WDS,ZU,T_AIRF,ZT,RH_AIRF,ZQ,PA_AIRF,T1(I,1),  & ! Uncommented for uniform hfx forcing 
!!              DLW_AIRF,DSW_AIRF,TAU,HSB,HLB,CORR(I))					! B Clark
              DLW_AIRF,DSW_AIRF,TAU,HSB,HLB,CORR(I),QA_AIRF)

         WTSURF(I) = DLW_AIRF+DSW_AIRF+HSB+HLB
         SWRAD(I)  = DSW_AIRF
!        transfer data to other meteo forcing
       T_AIR(I) =T_AIRF         !!BULK AIR TEMPERATURE
       RH_AIR(I)=RH_AIRF        !!RELATIVE HUMIDITY
       QA_AIR(I)=QA_AIRF        !!SPECIFIC HUMIDITY
       PA_AIR(I)=PA_AIRF        !!SURFACE PRESSURE
       DLW_AIR(I)=DLW_AIRF      !!DOWNWARD LONGWAVE RADIATION
       DSW_AIR(I)=DSW_AIRF      !!DOWNWARD SHORTWAVE RADIATION

       END DO

     ELSE IF(WINDTYPE == 'stress') THEN
       CALL PSTOP
     END IF

   END IF !! M_TYPE='uniform'
        
       
   IF(M_TYPE == 'non-uniform')THEN

!
!--- Heat flux and short wave radiation  --------------------------------------!
!

     CALL BRACKET(HFX1_TM,THOUR,L1,L2,FACT,UFACT,IERR)
        
     IF(IERR == -1)THEN
       RETURN
     ELSE 

       IF(WINDTYPE == 'speed')THEN
         DO I=1,MElem
           TX = UFACT*DTX(I,L1) + FACT*DTX(I,L2)
           TY = UFACT*DTY(I,L1) + FACT*DTY(I,L2)
           WDSN(I)=SQRT(TX*TX+TY*TY)
         END DO
# if defined (MULTIPROCESSOR)
         IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,WDSN)
# endif

         DO I=1,NNode
           WDS = 0.0_SP
           DO J=1,NTVE(I)
             WDS = WDS + WDSN(NBVE(I,J))
           END DO
           WDS = WDS/FLOAT(NTVE(I))
  
           T_AIRF    = UFACT*T_AIRN(I,L1)   + FACT*T_AIRN(I,L2) 
           RH_AIRF   = UFACT*RH_AIRN(I,L1)  + FACT*RH_AIRN(I,L2) 
           PA_AIRF   = UFACT*PA_AIRN(I,L1)  + FACT*PA_AIRN(I,L2) 
           DLW_AIRF  = UFACT*DLW_AIRN(I,L1) + FACT*DLW_AIRN(I,L2) 
           DSW_AIRF  = UFACT*DSW_AIRN(I,L1) + FACT*DSW_AIRN(I,L2) 
!        for ice model
!----------------------------------------------
!----------------------------------------------
           QA_AIR(I) = UFACT*QA_AIRN(I,L1) + FACT*QA_AIRN(I,L2) 
           CLOUD(I)  = UFACT*CLOUDN(I,L1)   + FACT*CLOUDN(I,L2) 

          PA_AIR(I)=PA_AIRF
          T_AIR(I) =T_AIRF
          RH_AIR(I)=RH_AIRF
          DLW_AIR(I)= DLW_AIRF
          DSW_AIR(I)= DSW_AIRF
        if(msr)then
!         write(101,'(10f10.2)')vx(i),vy(i),T_AIR(I),QA_AIR(I)
        end if

!----------------------------------------------
!----------------------------------------------
!!         PA  --->  mb
!!         Ta--->deg in ----Kelvin  out

           CALL COARE26Z(WDS,ZU,T_AIRF,ZT,RH_AIRF,ZQ,PA_AIRF/100.,T1(I,1),  &
!                DLW_AIRF,DSW_AIRF,TAU,HSB,HLB,CORR(I))
                DLW_AIRF,DSW_AIRF,TAU,HSB,HLB,CORR(I),QA_AIR(I))
!           WTSURF(I) = DLW_AIRF+DSW_AIRF+HSB+HLB
!           SWRAD(I)  = DSW_AIRF
!           QA_AIR(I) = UFACT*QA_AIRN(I,L1) + FACT*QA_AIRN(I,L2)

!         ggao
        if(msr)then
!         write(102,'(10f10.2)')vx(i),vy(i),T_AIR(I),QA_AIR(I)
        end if
!        T_AIRF = T_AIRF -273.15

         END DO
       ELSE IF(WINDTYPE == 'stress') THEN
         CALL PSTOP
       END IF
     END IF

# if defined (MULTIPROCESSOR)

         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,PA_AIR)
         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,T_AIR)
         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,RH_AIR)
         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,QA_AIR)
         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,DLW_AIR)
         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,DSW_AIR)
!         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,SWRAD)
!         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,WTSURF)
         IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,CLOUD)
# endif
         do i=1,NTNode
!         if(msr)then
!         write(102,'(10f10.4)')vx(i),vy(i),T_AIR(I),QA_AIR(I)
!         end if
         end do
 
   END IF !! MTYPE='non-uniform'

   !EL_AIR = -(PA_AIR-SLP0)/(1025._SP*GRAV_N)
   
   RETURN
   END SUBROUTINE BCOND_HFX
   
   
   
   
!==============================================================================|
!     SET BOUNDARY CONDITIONS FOR CALCULATING atmosphere pressure external model|
!                                                                              |
!==============================================================================|

   SUBROUTINE BCOND_PA_AIR

!==============================================================================|
   USE ALL_VARS
   USE BCS

#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   REAL(SP) :: FACT,UFACT,TX,TY
!    for ice model
   REAL(SP) :: pa_air2
   INTEGER  I,J,L1,L2,IERR

   REAL(SP):: Thour_P

   Thour_P = TIME * 24.0_SP



!==============================================================================|
!   SURFACE FORCING FOR INTERNAL MODE                                          !
!==============================================================================|
   IF(M_TYPE == 'uniform') THEN

!
!--- PA for the external Mode-----------------------------------------!
!
     !CALL BRACKET(UMF1_TM,THOUR1,L1,L2,FACT,UFACT,IERR)
     !PA_AIR2= UFACT*PA_AIRU(L1) + FACT*PA_AIRU(L2)        !!SURFACE PRESSURE
     ELF_AIR = 0.0 !-(PA_AIR2-SLP0)/(1025._SP*GRAV_N)

   END IF !! M_TYPE='uniform'
        
       
   IF(M_TYPE == 'non-uniform')THEN

     CALL BRACKET(HFX1_TM,THOUR_P,L1,L2,FACT,UFACT,IERR)
        
     IF(IERR == -1)THEN
       RETURN
     ELSE 

      DO I=1,NNode
          PA_AIR2 = UFACT*PA_AIRN(I,L1)  + FACT*PA_AIRN(I,L2) 
          ELF_AIR(I) = -(PA_AIR2-SLP0)/(1025._SP*GRAV_N(I))
      END DO
     END IF

# if defined (MULTIPROCESSOR)
     IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,ELF_AIR)
# endif
 
   END IF !! MTYPE='non-uniform'

   
   RETURN
   END SUBROUTINE BCOND_PA_AIR
   
   
!==============================================================================|

!==============================================================================|
   SUBROUTINE COARE26Z(UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,TAU,HSB,HLB,LAT,Q)
!   SUBROUTINE COARE26Z(UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,TAU,HSB,HLB,LAT)
!  A=coare26(u,zu,Ta,zt,rh,zq,Pa,Ts,dlw,dsw)
! Simplified non-vectorized version of coare2.6 code
! with cool skin option retained but warm layer and 
! surface wave options removed, and rain set to zero. 
! Assumes input are single scalars and that u is the
! magitude of the difference between wind and surface 
! current vectors, if latter available. Output:
! A = [tau qsen qlat Cd Ch Ce Cdn_10 Chn_10 Cen_10].
 
   REAL(SP) :: UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,TAU,HSB,HLB,LAT,LE,L,JCOOL,L10
   REAL(SP) :: U,US,T,P,RL,RS,QSAT26S,QS,Q,ZI,BETA,VON,FDG,TDK,GRVS,GRAV,RGAS
   REAL(SP) :: CPA,CPV,RHOA,VISA,AL,BE,CPW,RHOW,VISW,TCW,BIGC,WETC,RNS,RNL,DU
   REAL(SP) :: DT,DQ,UG,DTER,DQER,UT,U10,USR,ZO10,CD10,CH10,CT10,ZOT10,CD,CT,CC
   REAL(SP) :: RIBCU,RIBU,NITS,ZETU,PSIU_26S,PSIT_26S,TSR,QSR,TKT,CHARN,ZET
   REAL(SP) :: ZO,ZOQ,ZOT,BF,QOUT,DELS,QCOL,ALQ,XLAMX,CH,CE,CDN_10,CHN_10,CEN_10,RR
   
   INTEGER  :: I
! Set jcool=0 if Ts is surface, =1 if Ts is bulk.
! rcb checked 6/9/04
! set jcool=1 if Ts is bulk, 0 if Ts is true skin jcool=1;
   JCOOL=1.
! rename variables from fairall et al coare3 code
 

! wind speed (m/s) at height zr (m)
   U=UR
!  surface current speed in the wind direction(m/s)
   US=0*UR
!  water temperature (deg C)
   TS=TS
!  BULK AIR TEMPERATURE (C) AT HEIGHT ZT(m)
   T=TA
!  RELATIVE HUMIDITY (%) AT HEIGHT zq(M)
   RH=RH
!  SURFACE PRESSURE (mb)
   P=PA
!  DOWNWARD LONGWAVE RADIATION (W/m2)
   RL=DLW
!  DOWNWARD SHORTWAVE RADIATION (W/m2)
   RS=DSW
!  CONVERT RH TO SPECIFIC HUMIDITY (G/KG)
   CALL QSAT26(TS,P,QSAT26S)
   QS=0.98_SP*QSAT26S/1000._SP
!  SPECIFIC HUMIDITY OF AIR (G/KG)  
   CALL QSAT26(T,P,QSAT26S)
   Q=(0.01_SP*RH)*QSAT26S/1000._SP
!            print *,'qq',Q
!   SET RAIN TO ZERO RAIN=0*U
!   SET RAIN RATE (MM/HR) - KEEP AS OPTION

!  ***********SET LOCAL CONSTANTS *********
!    PBL HEIGHT (M)
   ZI=600._SP
!     LATITUDE (DEG,N=+)- GEORGES BANK
!        LAT=42.
! ************SET CONSTANTS **************
   BETA=1.2_SP
   VON=0.4_SP
   FDG=1.00_SP
   TDK=273.16_SP
   CALL GRV(LAT,GRVS)
   GRAV=GRVS
! ************AIR CONSTANTS **************
   RGAS=287.1_SP
   LE=(2.501_SP-0.00237_SP*TS)*1000000._SP
   CPA=1004.67_SP
   CPV=CPA*(1._SP+0.84_SP*Q)
   RHOA=P*100._SP/(RGAS*(T+TDK)*(1+0.61_SP*Q))
   VISA=1.326_SP*0.00001_SP*(1+6.542_SP*0.001_SP*T+  &
        8.301_SP*0.000001_SP*T*T-4.84_SP*0.000000001_SP*T*T*T)
! ***********COOL SKIN CONSTANTS******************
   AL=2.1_SP*0.00001_SP*((TS+3.2_SP)**0.79_SP)
   BE=0.026_SP
   CPW=4000._SP
   RHOW=1022._SP
   VISW=0.000001_SP
   TCW=0.6_SP
   BIGC=16._SP*GRAV*CPW*((RHOW*VISW)**3)/(TCW*TCW*RHOA*RHOA)
   WETC=0.622_SP*LE*QS/(RGAS*((TS+TDK)**2))
! ***************COMPUTE AUX STUFF *********
   RNS=RS*0.945_SP
   RNL=0.97_SP*(5.67_SP*0.00000001_SP*((TS-0.3_SP*JCOOL+TDK)**4)-RL)


! **************BEGIN BULK LOOP ***********

! **************FIRST GUESS ***************
   DU=U-US
   DT=TS-T-0.0098_SP*ZT
   DQ=QS-Q
   TA=T+TDK
   UG=0.5_SP
   DTER=0.3_SP
   DQER=WETC*DTER
   UT=SQRT(DU*DU+UG*UG) 
   U10=UT*LOG(10._SP/1e-4)/LOG(ZU/1e-4)
   USR=0.035_SP*U10

   ZO10=0.011_SP*USR*USR/GRAV+0.11_SP*VISA/USR
   CD10=(VON/LOG(10._SP/ZO10))**2
   CH10=0.00115_SP
   CT10=Ch10/SQRT(CD10)
   ZOT10=10._SP/EXP(VON/CT10)
   CD=(VON/LOG(ZU/ZO10))**2
   CT=VON/LOG(ZT/ZOT10)
   CC=VON*CT/CD
   RIBCU=-ZU/(ZI*0.004_SP*(BETA**3))
   RIBU=-GRAV*ZU*((DT-DTER*JCOOL)+.61_SP*TA*DQ)/(TA*(UT**2))
!       same as edson
   NITS=6._SP   
   IF(RIBU < 0)THEN
     ZETU=CC*RIBU/(1._SP+RIBU/RIBCU)
   ELSE
     ZETU=CC*RIBU/(1._SP+27._SP/(9*RIBU*CC))
   END IF
   L10=ZU/ZETU
   IF(ZETU > 50)THEN
     NITS=1
   END IF
   CALL PSIU_26(ZU/L10,PSIU_26S)
   USR=UT*VON/(LOG(ZU/ZO10)-PSIU_26S)
   CALL PSIT_26(ZT/L10,PSIT_26S)
   TSR=-(DT-DTER*JCOOL)*VON*FDG/(LOG(ZT/ZOT10)-PSIT_26S)
   CALL PSIT_26(ZQ/L10,PSIT_26S)
   QSR=-(DQ-WETC*DTER*JCOOL)*VON*FDG/(LOG(ZQ/ZOT10)-PSIT_26S)
   TKT=.001_SP
   CHARN=0.011_SP
   IF(UT > 10)THEN 
     CHARN=0.011_SP+(UT-10)/(18-10)*(0.018_SP-0.011_SP)
   END IF
   IF(UT > 18)THEN
     CHARN=0.018_SP
   END IF

!*************** bulk loop ************
!    do I=1,nits
   DO I=1,1
     ZET=VON*GRAV*ZU/TA*(TSR*(1._SP+0.61_SP*Q)+0.61_SP*TA*QSR)/   &
         (USR*USR)/(1._SP+0.61_SP*Q)
     ZO=CHARN*USR*USR/GRAV+0.11_SP*VISA/USR
     RR=ZO*USR/VISA
     L=ZU/ZET
#    if defined (DOUBLE_PRECISION)
     ZOQ=DMIN1(1.15e-4_SP,5.5e-5_SP/(RR**0.6_SP)) 
#    else
     ZOQ=AMIN1(1.15e-4_SP,5.5e-5_SP/(RR**0.6_SP)) 
#    endif     
     ZOT=ZOQ
     CALL PSIU_26(ZU/L,PSIU_26S)
     USR=UT*VON/(LOG(ZU/ZO)-PSIU_26S)
     CALL PSIT_26(ZT/L,PSIT_26S)
     TSR=-(DT-DTER*JCOOL)*VON*FDG/(LOG(ZT/ZOT)-PSIT_26S)
     CALL PSIT_26(ZQ/L,PSIT_26S)
     QSR=-(DQ-WETC*DTER*JCOOL)*VON*FDG/(LOG(ZQ/ZOQ)-PSIT_26S) 
     BF=-GRAV/TA*USR*(TSR+.61_SP*TA*QSR)               
!      print *,usr,tsr,qsr,Bf
     IF(BF > 0)THEN
       UG=BETA*((BF*ZI)**0.333_SP)
     ELSE
       UG=0.2_SP
     END IF
     UT=SQRT(DU*DU+UG*UG) 
     RNL=0.97_SP*(5.67_SP*0.00000001_SP*((TS-DTER*JCOOL+TDK)**4)-RL)
     HSB=-RHOA*CPA*USR*TSR

     HLB=-RHOA*LE*USR*QSR
     QOUT=RNL+HSB+HLB
 
     DELS=RNS*(.065_SP+11*TKT-6.6_SP*0.00001_SP/(TKT*(1-EXP(-TKT/8.0_SP*0.0001_SP))))
     QCOL=QOUT-DELS
     ALQ=AL*QCOL+BE*HLB*CPW/LE
     IF(ALQ > 0)THEN
       XLAMX=6._SP/(1._SP+(BIGC*ALQ/USR**4)**0.75_SP)**0.333_SP
       TKT=XLAMX*VISW/(SQRT(RHOA/RHOW)*USR)
     ELSE
       XLAMX=6.0_SP
#    if defined (DOUBLE_PRECISION)
       TKT=DMIN1(.01_SP,XLAMX*VISW/(SQRT(RHOA/RHOW)*USR))
#    else
       TKT=AMIN1(.01_SP,XLAMX*VISW/(SQRT(RHOA/RHOW)*USR))
#    endif     
     END IF
!           cool skin
     DTER=QCOL*TKT/TCW
     DQER=WETC*DTER
   END DO
!  of  end bulk iter loop

!****** compute fluxes ******************************************
!             wind stress
   TAU=RHOA*USR*USR*DU/UT
!              sensible heat flux
   HSB=RHOA*CPA*USR*TSR
!              latent heat flux
   HLB=RHOA*LE*USR*QSR


!****** compute transfer coeffs relative to ut @ meas. ht ********
#  if defined (DOUBLE_PRECISION)
   CD=TAU/RHOA/UT/DMAX1(.1_SP,DU)
#  else   
   CD=TAU/RHOA/UT/AMAX1(.1_SP,DU)
#  endif   
   CH=-USR*TSR/UT/(DT-DTER*JCOOL)
   CE=-USR*QSR/(DQ-DQER*JCOOL)/UT

!****** compute 10-m neutral coeff relative to ut ****************
   CDN_10=VON*VON/LOG(10._SP/ZO)/LOG(10._SP/ZO)
   CHN_10=VON*VON*FDG/LOG(10._SP/ZO)/LOG(10._SP/ZOT)
   CEN_10=VON*VON*FDG/LOG(10._SP/ZO)/LOG(10._SP/ZOQ)

!******** rain heat flux (save to use if desired) *************
! dwat=2.11e-5*((t+tdk)/tdk)^1.94; %! water vapour diffusivity
! dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa); %!heat diffusivity
! alfac= 1/(1+(wetc*Le*dwat)/(cpa*dtmp)); %! wet bulb factor
! RF= rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600;
!**************************************************************
!----------------------------------------------------------
   RETURN
   END SUBROUTINE COARE26Z
   

!-----------------------------------------------------------------------|
   SUBROUTINE PSIT_26(ZET,PSI)
! computes temperature structure function
   IMPLICIT NONE
   REAL(SP) :: ZET,PSI,X,PSIK,PSIC,F,C
   
   IF(ZET < 0)THEN
     X=(1._SP-15._SP*ZET)**.5_SP
     PSIK=2._SP*LOG((1._SP+X)/2._SP)
     X=(1._SP-34.15_SP*ZET)**.3333_SP
     PSIC=1.5_SP*LOG((1._SP+X+X*X)/3._SP)-SQRT(3._SP)    &
          *ATAN((1._SP+2._SP*X)/SQRT(3._SP))+4._SP*ATAN(1._SP)/SQRT(3._SP)
     F=ZET*ZET/(1._SP+ZET*ZET)
     PSI=(1._SP-F)*PSIK+F*PSIC
   ELSE
#    if defined (DOUBLE_PRECISION)
     C=DMIN1(50._SP,.35_SP*ZET)
#    else     
     C=AMIN1(50._SP,.35_SP*ZET)
#    endif     
     PSI=-((1._SP+2._SP/3._SP*ZET)**1.5_SP+.6667_SP     &
         *(ZET-14.28_SP)/EXP(C)+8.525_SP)
   END IF
   
   RETURN
   END SUBROUTINE PSIT_26
  

!----------------------------------------------------------
   SUBROUTINE PSIU_26(ZET,PSI)
! computes velocity structure function
   IMPLICIT NONE
   REAL(SP) :: ZET,PSI,X,PSIK,PSIC,F,C

   IF(ZET < 0)THEN
     X=(1._SP-15._SP*ZET)**.25_SP
     PSIK=2._SP*LOG((1._SP+X)/2._SP)+LOG((1._SP+X*X)/2._SP)   &
          -2._SP*ATAN(X)+2._SP*ATAN(1._SP)
     X=(1._SP-10.15_SP*ZET)**.3333_SP
     PSIC=1.5_SP*LOG((1._SP+X+X*X)/3._SP)-SQRT(3._SP)   &
          *ATAN((1._SP+2._SP*X)/SQRT(3._SP))+4._SP*ATAN(1._SP)/SQRT(3._SP)
     F=ZET*ZET/(1._SP+ZET*ZET)
     PSI=(1._SP-F)*PSIK+F*PSIC
   ELSE
#    if defined (DOUBLE_PRECISION)
     C=DMIN1(50._SP,.35_SP*ZET)
#    else     
     C=AMIN1(50._SP,.35_SP*ZET)
#    endif     
     PSI=-((1._SP+1.0_SP*ZET)**1.0_SP+.667_SP*(ZET-14.28_SP)/EXP(C)+8.525_SP)
   END IF
   
   RETURN
   END SUBROUTINE PSIU_26


!-----------------------------------------------------------
   SUBROUTINE QSAT26(T,P,QS)
! computes saturation specific humidity
   IMPLICIT NONE
   REAL(SP) :: T,P,QS,ES
   
   ES=6.112_SP*EXP(17.502_SP*T/(T+241.0_SP))     &
      *(1.0007_SP+3.46_SP*0.000001_SP*P)
   QS=ES*622/(P-.378_SP*ES)
   
   RETURN
   END SUBROUTINE QSAT26


!------function g=grv(lat)----------------------------------------------|
   SUBROUTINE GRV(LAT,G)
   IMPLICIT NONE
   REAL(SP) :: PI,GAMMA,C1,C2,C3,C4,PHI,X,G,LAT

   PI=3.1415926_SP
! computes g given lat in deg
   GAMMA=9.7803267715_SP
   C1=0.0052790414_SP
   C2=0.0000232718_SP
   C3=0.0000001262_SP
   C4=0.0000000007_SP
   PHI=LAT*DEG2RAD
   X=SIN(PHI)
   G=GAMMA*(1._SP+C1*X**2+C2*X**4+C3*X**6+C4*X**8)
   
   RETURN
   END SUBROUTINE GRV

# endif

END MODULE MOD_HEATFLUX


