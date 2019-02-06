        MODULE string_operation_mod
        CONTAINS

                function readCString()
                      use iso_c_binding
                      character (kind=c_char, len=512)  :: readCString
                      character (kind=c_char, len=511) :: str
                      read(1, *) str
                      readCString = TRIM(str) // C_NULL_CHAR
                      return
                end function readCString

                FUNCTION cascvt(strng)
                ! CASe ConVerTer
                ! function to convert all upper case to lower case of a string
                ! strng: input string
                ! lngt:  length of string strng

                INTEGER lngt,i
                CHARACTER(LEN=*)strng
                CHARACTER(LEN=LEN_TRIM(strng))cascvt
                
                INTEGER jgap

                lngt=LEN_TRIM(strng)

                jgap=ICHAR('a')-ICHAR('A')
                ! this could be replaced by jgap=32, since the ASCII code for "a" is
                ! 97 and for "A" is 65, using the ICHAR function ensures robustness
                ! in case other systems are used (like EBCDIC--hey I'm an old timer)
                IF(lngt.GT.0) THEN
                  ! if string is empty, nothing needs to be changed
                   DO i=1,lngt
                     IF(strng(i:i).LE.'Z') THEN
                       IF(strng(i:i).GE.'A') strng(i:i)=CHAR(ICHAR(strng(i:i))+jgap)
                     ENDIF
                   ENDDO
                ENDIF
                cascvt=strng

                RETURN
                END FUNCTION cascvt

        END MODULE string_operation_mod


        PROGRAM fvcom_medm2txt
           !Program to convert the medm/CASENAME_sim*** files to text format for easy reading in matlab
           !
           !Creation: Wen Long @ PNNL,  May 3, 2012
           !
           USE string_operation_mod

           INTEGER, PARAMETER :: &
               ISHORT = SELECTED_INT_KIND(8),            &    !Short integer      8-digit
               ILONG = SELECTED_INT_KIND(16),            &    !Long integer      16-digit
               SPREC = SELECTED_REAL_KIND(6,30),         &    !singple precision, maximum value is 6^30
               DPREC = SELECTED_REAL_KIND(12,300),       &    !double precision,  maximum value is 12^300
               STRLEN_SHORT = 72,                        &    !Short string length
               STRLEN_LONG  = 1024                            !Long string length

           DOUBLE PRECISION, ALLOCATABLE :: u(:,:),  &  !Eastward velocity (m/s)
                                v(:,:),  &  !Northward velocity (m/s)
                                w(:,:),  &  !Upward velocity (m/s)
                               KM(:,:),  &  !Momentum Diffusivity (m^2/s)
                                 el(:),  &  !Surface elevation (m)
                                T(:,:),  &  !Temperature (degC)
                                S(:,:),  &  !Salinity (psu)
                              RHO(:,:)      !density (kg/m^3)

!           REAL TMP1, TMP2

           INTEGER :: IINTT,   & !internal mode time integration step
                      NGL,     & !number of elements
                      MGL        !number of nodes

           DOUBLE PRECISION :: THOUR            !time in hours 

           CHARACTER(LEN=1024) :: SIMFILE     !file name to read from 
           CHARACTER(LEN=1024) :: SIMFILE_TXT !file name to write text data
           CHARACTER(LEN=1024) :: FORMAT_STR
           CHARACTER(LEN=72)   :: VARNAME     !name of variable to read from file 

           LOGICAL :: BAD_INPUT = .TRUE.
           LOGICAL :: FILE_EXISTS = .FALSE. 

           INTEGER(ILONG)  :: KB,    &   !total number of levels
                              KBM1       !total number of layers (KBM1=KB-1)
           INTEGER(ILONG) :: I,K

           CHARACTER(LEN=1024) :: ARGUMENT_STR
           INTEGER(ISHORT)     :: NARG

           !obtain arguments from command line
           NARG=iargc()                      !Call iargc() function to get number of arguments


           IF(NARG<4)THEN
              WRITE(*,*)'Usage: fvcom_medm2txt INPUTFILE OUTPUTFILE VARNAME NLEVEL               '
              WRITE(*,*)'where                                                                   '
              WRITE(*,*)'     INPUTFILE  --- file name of FVCOM medm output file                 '
              WRITE(*,*)'     OUTPUTFILE --- file name for outputing the chosen variable         '
              WRITE(*,*)'     VARNAME    --- variable name to print into OUTPUTFILE              '
              WRITE(*,*)'                    VARNAME should be one of: U,V,W,KM,EL,T,S,RHO,THOUR '
              WRITE(*,*)'     NLEVEL     --- number of sigma levels (e.g. 11 for Puget Sound)    '
              WRITE(*,*)'                                                                        '
              WRITE(*,*)' e.g.: fvcom_medm2txt psm_sim6607.dat psm_sim6607_u.txt u 11            '
              STOP
           ELSE
               DO I=1,NARG
                  CALL GETARG(I,ARGUMENT_STR)
                  ! WRITE(*,*)'Argument', I, '=', TRIM(ARGUMENT_STR)
                  IF(I==1)SIMFILE=TRIM(ARGUMENT_STR)
                  IF(I==2)SIMFILE_TXT=TRIM(ARGUMENT_STR)
                  IF(I==3)THEN
                     VARNAME=cascvt(TRIM(ARGUMENT_STR))
                     !check to make sure variable name is correct
                     IF(  TRIM(VARNAME)/='u'  .AND.  &
                          TRIM(VARNAME)/='v'  .AND.  &
                          TRIM(VARNAME)/='w'  .AND.  &
                          TRIM(VARNAME)/='km' .AND.  &
                          TRIM(VARNAME)/='el' .AND.  &
                          TRIM(VARNAME)/='t'  .AND.  &
                          TRIM(VARNAME)/='s'  .AND.  &
                          TRIM(VARNAME)/='rho'.AND.  &
                          TRIM(VARNAME)/='thour'  )THEN
                          WRITE(*,*)'Oops, variable name :', TRIM(VARNAME), ' is invalid!'
                          STOP
                    ENDIF
                  ENDIF
                  IF(I==4)READ(ARGUMENT_STR,'(i10)')KB
               ENDDO
          ENDIF

!          WRITE(*,*)'SIMFILE     = ',TRIM(SIMFILE)
!          WRITE(*,*)'SIMFILE_TXT = ',TRIM(SIMFILE_TXT)
!          WRITE(*,*)'VARNAME     = ',TRIM(VARNAME)
!          WRITE(*,*)'KB          = ',KB

          KBM1= KB -1 

!           BAD_INPUT=.TRUE.
!           DO WHILE(BAD_INPUT)
!              WRITE(*,'(A)', ADVANCE='NO')'Please input MEDM SIM  file name:'
!              READ(*,*)SIMFILE
              INQUIRE(FILE=TRIM(SIMFILE), EXIST=FILE_EXISTS)   !test if file exist
              IF(FILE_EXISTS)THEN
                BAD_INPUT=.FALSE.
              ELSE
                WRITE(*,'(A)') 'Oops file does not exist, please try again'
                BAD_INPUT=.TRUE.
                STOP
              ENDIF
!           ENDDO


!           BAD_INPUT=.TRUE.
!           FILE_EXISTS=.FALSE.
!           DO WHILE(BAD_INPUT)
!              WRITE(*,'(A)', ADVANCE='NO')'Please give output txt file name:'
!              READ(*,*)SIMFILE_TXT
              INQUIRE(FILE=TRIM(SIMFILE_TXT), EXIST=FILE_EXISTS)   !test if file exist
              IF(FILE_EXISTS)THEN
                BAD_INPUT=.TRUE.
                WRITE(*,'(A)') 'Oops file ' //TRIM(SIMFILE_TXT)//' already exists, please try again'
                STOP
              ELSE
                BAD_INPUT=.FALSE.
              ENDIF
!           ENDDO

          !read the SIM file data and print as texts
          OPEN(UNIT=11,FILE=TRIM(SIMFILE),STATUS='OLD', FORM='UNFORMATTED')

            READ(11)IINTT, NGL, MGL, THOUR

!            WRITE(*,*)'IINTT,NGL,MGL,THOUR=',IINTT, NGL, MGL, THOUR

             IF(TRIM(VARNAME)=='thour')THEN
                OPEN(UNIT=12,FILE=TRIM(SIMFILE_TXT),STATUS='NEW');
                   WRITE(12,'(F20.10)')THOUR
                CLOSE(12)
                STOP         !Do not go further if only wants THOUR 
             ENDIF

            !ALLOCATE VARIABLES

             IF(.NOT.ALLOCATED(u))ALLOCATE(u(NGL,KBM1))
             IF(.NOT.ALLOCATED(v))ALLOCATE(v(NGL,KBM1))
             IF(.NOT.ALLOCATED(w))ALLOCATE(w(NGL,KBM1))
             IF(.NOT.ALLOCATED(KM))ALLOCATE(KM(NGL,KBM1))

             IF(.NOT.ALLOCATED(el))ALLOCATE(el(MGL))
             IF(.NOT.ALLOCATED(T))ALLOCATE(T(MGL,KBM1))
             IF(.NOT.ALLOCATED(S))ALLOCATE(S(MGL,KBM1))
             IF(.NOT.ALLOCATED(RHO))ALLOCATE(RHO(MGL,KBM1))

              u=0;               v=0;                w=0; 
              KM=0;              el=0; 
              T=0;               S=0;               RHO=0;

            !read u,v,w and KM

             DO I=1,NGL
                READ(11) (u(I,K),v(I,K),w(I,K),KM(I,K),K=1,KBM1)
             ENDDO

            !read elevation, T, S, Density
             DO I=1,MGL
                READ(11) el(I),(T(I,K),S(I,K),RHO(I,K),K=1,KBM1)
             ENDDO

          CLOSE(11)


          !Open output file and dump data
           FORMAT_STR='(:D30.12,1000(:1X,D30.12))'
           OPEN(UNIT=12,FILE=TRIM(SIMFILE_TXT),STATUS='NEW');
             SELECT CASE (TRIM(VARNAME))

                CASE ('u')
                    DO I=1,NGL
                       WRITE(12,FORMAT_STR)(u(I,K),K=1,KBM1)
                    ENDDO

                CASE ('v')
                    DO I=1,NGL
                       WRITE(12,FORMAT_STR)(v(I,K),K=1,KBM1)
                    ENDDO

                CASE ('w')
                    DO I=1,NGL
                       WRITE(12,FORMAT_STR)(w(I,K),K=1,KBM1)
                    ENDDO

                CASE ('km')
                    DO I=1,NGL
                       WRITE(12,FORMAT_STR)(KM(I,K),K=1,KBM1)
                    ENDDO

                CASE ('el')
                    DO I=1,MGL
                       WRITE(12,FORMAT_STR) el(I)
                    ENDDO

                CASE ('t')
                    DO I=1,MGL
                       WRITE(12,FORMAT_STR)(T(I,K),K=1,KBM1)
                    ENDDO

                CASE ('s')
                    DO I=1,MGL
                       WRITE(12,FORMAT_STR)(S(I,K),K=1,KBM1)
                    ENDDO

                CASE ('rho')
                    DO I=1,MGL
                       WRITE(12,FORMAT_STR)(RHO(I,K),K=1,KBM1)
                    ENDDO
                CASE DEFAULT
                    WRITE(*,'(A)')'Wrong choice of variable, please try again'
                    GOTO 100
            END SELECT

100         CLOSE(12)

           !DEALLOCATE variables

            IF(ALLOCATED(u))DEALLOCATE(u);
            IF(ALLOCATED(v))DEALLOCATE(v);
            IF(ALLOCATED(w))DEALLOCATE(w);
            IF(ALLOCATED(KM))DEALLOCATE(KM);

            IF(ALLOCATED(el))DEALLOCATE(el);
            IF(ALLOCATED(T))DEALLOCATE(T);
            IF(ALLOCATED(S))DEALLOCATE(S);
            IF(ALLOCATED(RHO))DEALLOCATE(RHO);

            STOP

        END PROGRAM

