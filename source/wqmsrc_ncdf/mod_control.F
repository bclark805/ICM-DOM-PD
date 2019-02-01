
        !==============================================================================|
        !   CONTROL VARIABLES                                                          |
        !==============================================================================|

        MODULE MOD_CONTROL

           USE MOD_PREC, ONLY : SP!
           IMPLICIT NONE
           SAVE

           LOGICAL SERIAL                  !!TRUE IF SINGLE PROCESSOR
           LOGICAL MSR                     !!TRUE IF MASTER PROCESSOR (MYID==1)
           LOGICAL PAR                     !!TRUE IF MULTIPROCESSOR RUN

           CHARACTER(LEN=80) CASENAME      !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
           CHARACTER(LEN=120) CASETITLE    !!CASE TITLE                                 

           REAL(SP) :: HMAX         !!GLOBAL MAXIMUM DEPTH
           REAL(SP) :: HMIN         !!GLOBAL MINIMUM DEPTH
           REAL(SP) :: UMOL         !!VERTICAL DIFFUSION COEFFICIENT
           REAL(SP) :: HORCON       !!HORIZONTAL DIFFUSION COEFFICIENT
           REAL(SP) :: DTI          !!internal time step

           CHARACTER(LEN=80) HORZMIX       !!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
           CHARACTER(LEN=100) :: FILENUMBER,PREFF,INPDIR,GEOAREA,RIV_FILENUMBER
           CHARACTER(LEN=80) INFLOW_TYPE   !!SPECIFIED RIVER INFLOW TYPE (edge/node) 
           CHARACTER(LEN=80) POINT_ST_TYPE !!(calculated/specified)
           CHARACTER(LEN=80) PNT_SOURCE    !!point_source
           INTEGER :: DAY, in_jday
           
        END MODULE MOD_CONTROL
