!------------------------------------------------------------------------------------!
!----This is the output_flux module so that the model can output fluxes on TCE edges 
! for selected TCE edges, and later we can aggregate and calculate fluxes through 
! a transect
!
! Wen Long, 04/28/2014

MODULE MOD_OUTPUT_FLUX
# if defined (OUTPUT_FLUX)
 ! USE MOD_TYPES
   USE MOD_PREC
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   IMPLICIT NONE
   SAVE
   
   INTEGER :: NEL_GL    ! number of TCE edges to output in global domain
   INTEGER :: NEL_LOC   ! number of TCE edges in local domain to output
   
   INTEGER,  ALLOCATABLE :: I_edgeline_ele_GL(:)  !global element ID of TCE edges (each TCE edge belongs to an element)
   INTEGER,  ALLOCATABLE :: I_edgeline_ele(:)     !local element ID of TCE edges
   
   INTEGER, ALLOCATABLE  :: I_edgeline_ele_edge_ia_GL(:) !global node number of global TCE edge's ia point
   INTEGER, ALLOCATABLE  :: I_edgeline_ele_edge_ib_GL(:) !global node number of global TCE edge's ib point

   INTEGER, ALLOCATABLE  :: I_edgeline_ele_edge_ia(:) !local node number of local TCE edge's ia point
   INTEGER, ALLOCATABLE  :: I_edgeline_ele_edge_ib(:) !local node number of local TCE edge's ib point
   
   INTEGER, ALLOCATABLE  :: I_edgeline_flux_sign_GL(:)   !sign to apply to global TCE edge's flux
                                                         !after applying this sign to each global TCE edge
                                                         !if positive, flux is positive towards ib point
                                                         !if negative, flux is positive towards ia point
   INTEGER, ALLOCATABLE :: I_edgeline_transect_ID_GL(:)  !transect ID (starts from 1) for the edge line segments
                                                         !those TCE edges with same transect ID will be used in
                                                         !postprocessing for grouping to the same transect

   INTEGER, ALLOCATABLE :: I_edgeline_flux_sign(:)    !local sign to apply to TCE edge's flux

   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_SALT_GL(:,:)        !salt flux through the TCE edge
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_DIFFUSION_SALT_GL(:,:) !diffusive salt flux through the TCE edge 
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_ADVECTION_SALT_GL(:,:)  !advective salt flux through the TCE edge = salt flux - diffusive salt flux
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_ADVECTION_SALTWATER_GL(:,:)   !saltwater flux through the TCE edge

   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_FRSHWATER_GL(:,:)   !freshwater flux through the TCE edge
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_WATER_GL(:,:)       !total (advection salt water + fresh water) volume flux through the TCE edge
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_DZ_GL(:,:)          !thickness of each layer on the edge line (m)   

   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_SALT(:,:)            !salt flux through the local TCE edges
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_DIFFUSION_SALT(:,:)  !diffusive salt flux through the TCE edge
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_ADVECTION_SALT(:,:)  !advective salt flux through the TCE edge = salt flux - diffusive salt flux

   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_ADVECTION_SALTWATER(:,:) !saltwater flux through the TCE edge
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_FRSHWATER(:,:)   !freshwater flux through the local TCE edge
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_WATER(:,:)       !total (salt water + fresh water) volume flux through the local TCE edge
   REAL(SP) , ALLOCATABLE :: EDGELINE_FLUX_DZ(:,:)          !thickness of each layer on the edge line (m)

   
   INTEGER :: INF_TCE_EDGE_LINES
   LOGICAL :: WRITE_TCE_FLUX
   REAL(SP) ::Smax_Flux       !maximum salinity  for water that is deemed as full salt water in flux calculation
   
   INTEGER :: Flux_INT        !interval (internal mode step of of writing out flux)
   
CONTAINS    !-----------------------------------------------------------------------|
            ! DATA_RUN_FLUX  : Input Parameters Which Control MOD_OUTPUT_FLUX       |
            ! IOFILES_FLUX   : Open Input Files for FLUX MODULE                     |


!==============================================================================|
!   Input Parameters                                                           |
!==============================================================================|

   SUBROUTINE DATA_RUN_FLUX
!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_INP
   IMPLICIT NONE
   INTEGER  INTVEC(150), ISCAN
   CHARACTER(LEN=120) :: FNAME


!==============================================================================|
!   READ IN VARIABLES AND SET VALUES                                           |
!==============================================================================|

   FNAME = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------------|
! OUTPUT_FLUX module on/off flag, scan the run file to see if it is activated or not |
!------------------------------------------------------------------------------------|
   
   ISCAN = SCAN_FILE(TRIM(FNAME),"WRITE_TCE_FLUX",LVAL = WRITE_TCE_FLUX )
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING WRITE_TCE_FLUX: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "WRITE_TCE_FLUX" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
   
   !Also find the Smax_Flux value in the casename_run.dat file    
   
   ISCAN = SCAN_FILE(TRIM(FNAME),"Smax_Flux",FSCAL = Smax_flux )
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING Smax_flux: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD FLOATING VARIABLE "Smax_flux" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF

   !Also find the Flux_INT value in the casename_run.dat file    
   
   ISCAN = SCAN_FILE(TRIM(FNAME),"Flux_INT",ISCAL = Flux_INT )
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING Flux_INT: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(IPT,*)'PLEASE ADD FLOATING VARIABLE "Flux_INT" TO INPUT FILE'
     END IF
     CALL PSTOP
   END IF
   
   
!-- output diagnostic infor. to run file...
   IF(MSR) THEN
       IF(WRITE_TCE_FLUX)THEN
           WRITE(IPT,*)'!  # OUTPUT_FLUX MODULE :  ACTIVE'
           WRITE(IPT,*)'!  # Maximum Salinity is :  ', Smax_flux
       ELSE
           WRITE(IPT,*)'!  # OUTPUT_FLUX MODULE :  INACTIVE'
       END IF
   END IF
   
   RETURN
   END SUBROUTINE DATA_RUN_FLUX
!------------------------------------------------------------------------------|


!==============================================================================!
!   Open & Read Input Files for FLUX Calculation                                !
!==============================================================================!

   SUBROUTINE IOFILES_FLUX
   USE ALL_VARS
   USE MOD_UTILS
   USE CONTROL
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif

   IMPLICIT NONE

   INTEGER :: I,J,K,I1,                 &
              tce_edge_tce_id,          &     !global id of the tce edge
              tce_edge_ele_id,          &     !global id of the element that the tce edge is in
              tce_edge_ele_edge_id,     &     !global id of the element edge that the tce edge is intersecting with
              tce_edge_flux_sign              !sign to be applied to flux of the tce edge when outputing it

!   !Outputs of the module
!   REAL(SP), ALLOCATABLE :: edgeline_flux_salt(:,:),               & !size is N_TCE_EL by KBM1
!                            edgeline_flux_freshwater(:,:),         & !        N_TCE_EL by KBM1
!                            edgeline_flux_volume(:,:)                !        N_TCE_EL by KBM1

!   !Inputs of the module
!   INTEGER, ALLOCATABLE :: edgeline_tce_id(:),                   & ! --- the global index of triangle elements to which the TCE edge belongs
!                           edgeline_ele_id(:),                   & ! --- the global index of TCE edges that make up the TCE edge line
!                           edgeline_ele_edge_id(:),              & ! --- the global index of triangle edge with which the TCE edge is intersecting
!                           edgeline_flux_sign(:)                   ! --- the sign applied to the flux when outputing the flux on the edge
                                                                   !         so that the flux positive pointing to RHS of the edge

!   INTEGER :: N_TCE_EL,  &  !number of TCE edges each edge line has
!              NEL           !total number of edge lines


   CHARACTER(LEN=80)  :: ISTR
   INTEGER ::NCNT,ia_GL, ib_GL

   REAL(SP), ALLOCATABLE :: TEMP1(:), TEMP2(:), TEMP3(:), TEMP4(:)

   INF_TCE_EDGE_LINES = 105  !file unit for input definition of edge lines
                             !each edge line consists of the following input information:
                             !
                             !    edgeline_tce_id(1:NEL)      
                             !    edgeline_ele_id(1:NEL)     
                             !    edgeline_ele_edge_id(1:NEL) 
                             !    edgeline_flux_sign(1:NEL)  
       

   ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)
!
!-----------------OPEN INF_TCE_EDGE_LINES FILE (SPATIAL VARYING Kelp bed MAP)-----------!

   !
   !Ideally we should only allow master proc to read the data and then broadcast to others
   !
   
   CALL FOPEN(INF_TCE_EDGE_LINES, TRIM(ISTR)//'_flux_map.dat', "cfr")
   
   IF(MSR)THEN
        open(106,file='check_tceflux_input.out')  !print errors if any while reading input data of edgelines
   ENDIF

   NEL_GL       = 0
   
  !for each edge line, read number of TCE edges in the edge line 
  !    read ele_id, ele_edge_ia, ele_edge_ib and flux_sign
   READ(INF_TCE_EDGE_LINES,*)NEL_GL !read total number of TCE edges

   IF(MSR)THEN
        WRITE(106,*)NEL_GL
   ENDIF
   
   IF(NEL_GL>0)THEN
   
        !ALLOCATE global edge line arrays 
        ALLOCATE(I_edgeline_ele_GL(1:NEL_GL))
        ALLOCATE(I_edgeline_ele_edge_ia_GL(1:NEL_GL))
        ALLOCATE(I_edgeline_ele_edge_ib_GL(1:NEL_GL))
        ALLOCATE(I_edgeline_flux_sign_GL(1:NEL_GL))
        ALLOCATE(I_edgeline_transect_ID_GL(1:NEL_GL))
  
   
        ALLOCATE(EDGELINE_FLUX_SALT_GL(1:NEL_GL,1:KBM1))        !salt flux through the TCE edge
		ALLOCATE(EDGELINE_FLUX_DIFFUSION_SALT_GL(1:NEL_GL,1:KBM1))        !diffusive salt flux through the TCE edge
		ALLOCATE(EDGELINE_FLUX_ADVECTION_SALT_GL(1:NEL_GL,1:KBM1))        !advective salt flux through the TCE edge
        ALLOCATE(EDGELINE_FLUX_ADVECTION_SALTWATER_GL(1:NEL_GL,1:KBM1))   !advective saltwater flux through the TCE edge
        ALLOCATE(EDGELINE_FLUX_FRSHWATER_GL(1:NEL_GL,1:KBM1))   !freshwater flux through the TCE edge
        ALLOCATE(EDGELINE_FLUX_WATER_GL(1:NEL_GL,1:KBM1))       !total (salt water + fresh water) volume flux through the TCE edge
        ALLOCATE(EDGELINE_FLUX_DZ_GL(1:NEL_GL,1:KBM1))          !thickness of each layer of TCE edge

        !initialize the arrays as zeros
        I_edgeline_ele_GL=0;
        I_edgeline_ele_edge_ia_GL=0;
        I_edgeline_ele_edge_ib_GL=0;
        I_edgeline_flux_sign_GL=0;
        I_edgeline_transect_ID_GL=0;
  
        EDGELINE_FLUX_SALT_GL=0.0
		EDGELINE_FLUX_DIFFUSION_SALT_GL=0.0
		EDGELINE_FLUX_ADVECTION_SALT_GL=0.0
        EDGELINE_FLUX_ADVECTION_SALTWATER_GL=0.0
        EDGELINE_FLUX_FRSHWATER_GL=0.0
        EDGELINE_FLUX_WATER_GL=0.0
        EDGELINE_FLUX_DZ_GL=0.0
  
        DO I=1,NEL_GL
 
            READ(INF_TCE_EDGE_LINES,*)I_edgeline_transect_ID_GL(I), & !transect ID of the TCE edge
                                      I_edgeline_ele_GL(I),         & !global element IDs of the TCE edge
                                      I_edgeline_ele_edge_ia_GL(I), & !global node IDs of the TCE edge's ia point
                                      I_edgeline_ele_edge_ib_GL(I), & !global node IDs of the TCE edge's ib point
                                      I_edgeline_flux_sign_GL(I)      !sign to apply before write-out of flux

            IF(MSR)THEN
                WRITE(106,*)I_edgeline_transect_ID_GL(I),            &
                            I_edgeline_ele_GL(I),                    &
                            I_edgeline_ele_edge_ia_GL(I),            &
                            I_edgeline_ele_edge_ib_GL(I),            &
                            I_edgeline_flux_sign_GL(I)
            ENDIF

        ENDDO
 
        CLOSE(INF_TCE_EDGE_LINES)
        
        IF(MSR)THEN
            CLOSE(106)
        ENDIF
!
!!transfer to local domains from global (3 variables read from _flux_map.dat only.)
!!FOR SERIAL MODE FIRST...
!
        IF(SERIAL) THEN

            NEL_LOC=NEL_GL
            ALLOCATE(I_edgeline_ele(1:NEL_LOC))
            ALLOCATE(I_edgeline_ele_edge_ia(1:NEL_LOC))
            ALLOCATE(I_edgeline_ele_edge_ib(1:NEL_LOC))
            ALLOCATE(I_edgeline_flux_sign(1:NEL_LOC))   

            I_edgeline_ele(1:NEL_LOC)=I_edgeline_ele_GL(1:NEL_GL);
            I_edgeline_ele_edge_ia(1:NEL_LOC)=I_edgeline_ele_edge_ia_GL(1:NEL_GL)
            I_edgeline_ele_edge_ib(1:NEL_LOC)=I_edgeline_ele_edge_ib_GL(1:NEL_GL)
            I_edgeline_flux_sign(1:NEL_LOC)=I_edgeline_flux_sign_GL(1:NEL_GL)
     
            !Allocate the local flux arrays
            ALLOCATE(EDGELINE_FLUX_SALT(1:NEL_LOC,1:KBM1))     			 ; EDGELINE_FLUX_SALT=0.0                  !salt flux through the local TCE edges
            ALLOCATE(EDGELINE_FLUX_DIFFUSION_SALT(1:NEL_LOC,1:KBM1))     ; EDGELINE_FLUX_DIFFUSION_SALT=0.0        !salt flux through the local TCE edges
            ALLOCATE(EDGELINE_FLUX_ADVECTION_SALT(1:NEL_LOC,1:KBM1))     ; EDGELINE_FLUX_ADVECTION_SALT=0.0        !salt flux through the local TCE edges
            ALLOCATE(EDGELINE_FLUX_ADVECTION_SALTWATER(1:NEL_LOC,1:KBM1)); EDGELINE_FLUX_ADVECTION_SALTWATER=0.0   !saltwater flux through the local TCE edge
            ALLOCATE(EDGELINE_FLUX_FRSHWATER(1:NEL_LOC,1:KBM1))          ; EDGELINE_FLUX_FRSHWATER=0.0             !freshwater flux through the local TCE edge
            ALLOCATE(EDGELINE_FLUX_WATER(1:NEL_LOC,1:KBM1))              ; EDGELINE_FLUX_WATER=0.0                 !total (salt water + fresh water) volume flux through the local TCE edge
            ALLOCATE(EDGELINE_FLUX_DZ(1:NEL_LOC,1:KBM1))                 ; EDGELINE_FLUX_DZ=0.0                    !thickness of each layer (m)

       END IF

!FOR PARALLEL MODE
#  if defined (MULTIPROCESSOR)
       IF(PAR)THEN
            !individual processors pick you their portion of the global edge line
            !find number of TCE edges for the processor
            NCNT=0
            ALLOCATE(TEMP1(NEL_GL))
            ALLOCATE(TEMP2(NEL_GL))
            ALLOCATE(TEMP3(NEL_GL))
            ALLOCATE(TEMP4(NEL_GL))
            
            DO I=1,NEL_GL
        
                I1 = ELID( I_edgeline_ele_GL(I) )  !get the local element ID of the global element that the global TCE edge is in
         
                IF(I1 /= 0)THEN
                    NCNT = NCNT + 1
                    TEMP1(NCNT) = I1               !record local element number
           
                    ia_GL=I_edgeline_ele_edge_ia_GL(I)       !global ia node number of the TCE edge           
                    TEMP2(NCNT) = NLID(ia_GL)                !set local ia node number by converting from global ia node number 
                    ib_GL=I_edgeline_ele_edge_ib_GL(I)       !global ib node number of the TCE edge
                    TEMP3(NCNT) = NLID(ib_GL)                !set local ib node number by converting from global ib node number
                    TEMP4(NCNT) = I_edgeline_flux_sign_GL(I) !set the flux sign
                END IF
            END DO
       
            NEL_LOC = NCNT                    !number of TCE edge line edges in local domain

            IF(NCNT > 0)THEN
                ALLOCATE(I_edgeline_ele(1:NEL_LOC));             I_edgeline_ele=TEMP1(1:NCNT)
                ALLOCATE(I_edgeline_ele_edge_ia(1:NEL_LOC));     I_edgeline_ele_edge_ia=TEMP2(1:NCNT)
                ALLOCATE(I_edgeline_ele_edge_ib(1:NEL_LOC));     I_edgeline_ele_edge_ib=TEMP3(1:NCNT)
                ALLOCATE(I_edgeline_flux_sign(1:NEL_LOC))  ;     I_edgeline_flux_sign=TEMP4(1:NCNT)
            END IF

            DEALLOCATE(TEMP1,TEMP2,TEMP3,TEMP4) 

            !allocate the local flux arrays
            ALLOCATE(EDGELINE_FLUX_SALT(1:NEL_LOC,1:KBM1))               ; EDGELINE_FLUX_SALT=0.0                !salt flux through the local TCE edges
            ALLOCATE(EDGELINE_FLUX_DIFFUSION_SALT(1:NEL_LOC,1:KBM1))     ; EDGELINE_FLUX_DIFFUSION_SALT=0.0      !diffusive salt flux through the local TCE edges
            ALLOCATE(EDGELINE_FLUX_ADVECTION_SALT(1:NEL_LOC,1:KBM1))     ; EDGELINE_FLUX_ADVECTION_SALT=0.0      !advective salt flux through the local TCE edges
            ALLOCATE(EDGELINE_FLUX_ADVECTION_SALTWATER(1:NEL_LOC,1:KBM1)); EDGELINE_FLUX_ADVECTION_SALTWATER=0.0 !advective saltwater flux through the local TCE edge
            ALLOCATE(EDGELINE_FLUX_FRSHWATER(1:NEL_LOC,1:KBM1))          ; EDGELINE_FLUX_FRSHWATER=0.0           !freshwater flux through the local TCE edge
            ALLOCATE(EDGELINE_FLUX_WATER(1:NEL_LOC,1:KBM1))              ; EDGELINE_FLUX_WATER=0.0               !total (salt water + fresh water) volume flux through the local TCE edge
            ALLOCATE(EDGELINE_FLUX_DZ(1:NEL_LOC,1:KBM1))                 ; EDGELINE_FLUX_DZ=0.0                  !thickness of layers

       END IF
#  endif
    
    !
    !open output file for write by master, in future this should be replaced
    !by netcdf output file
    !and file name should be set in the casename_run.dat file
    !
    
       IF(MSR)THEN !Master opens a file for output 
    
        OPEN(107,FILE='tceflux_output.dat') 
    
       ENDIF
    
   ENDIF  !END of NEL_GL>0 
   
   RETURN
   END SUBROUTINE IOFILES_FLUX
!==============================================================================!

   SUBROUTINE WRITE_EDGELINE_FLUX
!============================================================================================!
!  Gather edgeline flux from local processors and writeout to output file by master processor!
!============================================================================================!
   
   USE ALL_VARS
   USE CONTROL
   USE MOD_PREC
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   
   IMPLICIT NONE
   INTEGER :: I,J,K,II
   
# if defined (MULTIPROCESSOR)  
   !gather values to the root only
   INTEGER :: NELWC    !Water column total number of TCE edges = NEL_LOC*KBM1
   INTEGER :: NEL_LOC1 !NEL_LOC +1
   INTEGER :: I_DISPLACEMENT
   INTEGER :: I_DISPLACEMENT_FLX
   INTEGER :: N_NEL_LOC  !number of local TCE edges in a proc
   INTEGER :: N_FLX_LOC  !number of local TCE edge fluxes in a proc
   INTEGER :: IE_EL_GL   !global element ID of a TCE edge
   INTEGER :: ia_EL_GL   !global node ID of a TCE edge's ia point
   INTEGER :: ib_EL_GL   !global node ID of a TCE edge's ib point
   
   !arrays for gathering edge line TCE edge local element ID and local node IDs for 
   INTEGER, ALLOCATABLE :: COUNTS_EL_ELE(:), DISPLACEMENTS_EL_ELE(:) 
   INTEGER, ALLOCATABLE :: TEMP1_INT_SND(:), TEMP2_INT_SND(:), TEMP3_INT_SND(:)!, TEMP4_INT_SND(:) !temp arrays for sending
   INTEGER, ALLOCATABLE :: TEMP1_INT_RCV(:), TEMP2_INT_RCV(:), TEMP3_INT_RCV(:)!, TEMP4_INT_RCV(:) !temp arrays for receiving
   
   !arrays for gathering flux to master proc
   INTEGER, ALLOCATABLE :: COUNTS_FLX(:), DISPLACEMENTS_FLX(:) !
   REAL(SP), ALLOCATABLE :: TEMP1_SND(:), TEMP2_SND(:), TEMP3_SND(:), TEMP4_SND(:),TEMP5_SND(:), TEMP6_SND(:),TEMP7_SND(:) !temp arrays for sending
   REAL(SP), ALLOCATABLE :: TEMP1_RCV(:), TEMP2_RCV(:), TEMP3_RCV(:), TEMP4_RCV(:),TEMP5_RCV(:), TEMP6_RCV(:),TEMP7_RCV(:) !temp arrays for receiving

   INTEGER :: IERR
   
IF(PAR)THEN       
   
   IF(MSR)THEN
       !number of data points for element numbers
        ALLOCATE(COUNTS_EL_ELE(1:NPROCS))               ;COUNTS_EL_ELE=0; !ARRAY to store number of data values from all processors
        ALLOCATE(DISPLACEMENTS_EL_ELE(1:NPROCS)); DISPLACEMENTS_EL_ELE=0; !ARRAY to store index in receiving buffer 
   
       !number of data points for fluxes (including vertical column)
        ALLOCATE(COUNTS_FLX(1:NPROCS))       ;    COUNTS_FLX=0;           !ARRAY to store number of data values from all processors
        ALLOCATE(DISPLACEMENTS_FLX(1:NPROCS));    DISPLACEMENTS_FLX=0.0;  !ARRAY to store index in receiving buffer 
   ENDIF

   !
   !gather TCE edge's local element and ia, ib node numbers so master proc can output flux in correct sequence
   !
   !We don't have a mapping like bmap for TCE edge lines in the model
   !so we would have to use the element number and ia, ib node numbers of the TCE edge to 
   !match local TCE edge with global TCE edge so that in the end the master proc knows
   !the correspondence between them and to output the sequences in I_edgeline_ele_GL correctly
   !
   !WLong: this information is static, so should be done for just once
   !       rather than passing it every time step (this is left for further refactoring and improvement)
   !
   
   NEL_LOC1=NEL_LOC+1          !at least send one data point (if NEL_LOC=0, then send one dummy value)
   CALL MPI_GATHER(NEL_LOC1,  1,MPI_INTEGER,COUNTS_EL_ELE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR) !gatger NELWC to COUNTS_EL_ELE
    !e.g. NEL_LOC=2,0,0,3 for totally 4 procs, then NEL_LOC1=COUNTS_EL_ELE=(3, 1,1,4) ==> total number of values to send is 9
    !==> DISPLACEMENTS_EL_ELE(1) = 1
    !                        (2) = 1+3 = 4
    !                        (3) = 4+1 = 5
    !                        (4) = 5+1 = 6
    !
   
    !SUM(COUNTS_EL_ELE) = 9
    !
    !so in the receiving buffer TEMP1_INT_RCV(1:3) is for 1st proc  
    !                                        (4:4) is for 2nd proc
    !                                        (5:5) is for 3rd proc
    !                                        (6:9) is for 4th proc
    !
   IF(MSR)THEN
        DISPLACEMENTS_EL_ELE(1)  = 0  !index for values received from first proc
        DO I=2,NPROCS                 !indices for values received from other procs
           DISPLACEMENTS_EL_ELE(I)=COUNTS_EL_ELE(I-1)+DISPLACEMENTS_EL_ELE(I-1)
        ENDDO
        
            !(1) = 1
            !(2) = 1+3 = 4
            !(3) = 4+1 = 5
            !(4) = 5+1 = 6


        ALLOCATE(TEMP1_INT_RCV(1:SUM(COUNTS_EL_ELE)))  !to receive local element number of the TCE edge
        ALLOCATE(TEMP2_INT_RCV(1:SUM(COUNTS_EL_ELE)))  !to receive local ia node number of the TCE edge
        ALLOCATE(TEMP3_INT_RCV(1:SUM(COUNTS_EL_ELE)))  !to receive local ib node number of the TCE edge
                 
   ENDIF
   
   !WRITE(*,*)'MYID=',MYID,'Here Here Here'  


   IF(NEL_LOC1>0)THEN
   
   !prepare to send local element and node information to master proc
   !Since NEL_LOC1 >=1, there is at least one value of zero to send, 
   !truly useful values starts from the second one in the arrays
   
    ALLOCATE(TEMP1_INT_SND(1:NEL_LOC1)); TEMP1_INT_SND=0; !to send local element number
    ALLOCATE(TEMP2_INT_SND(1:NEL_LOC1)); TEMP2_INT_SND=0; !to send local ia node nubmer of the TCE edge
    ALLOCATE(TEMP3_INT_SND(1:NEL_LOC1)); TEMP3_INT_SND=0; !to send local ib node number of the TCE edge
    
   ENDIF
   
   !fill out the sending buffers
   IF(NEL_LOC>0)THEN
                
    DO I=1,NEL_LOC  !loop through all selected local TCE edges    
        !not values are put in sending buffers starting from index 2
        !if NEL_LOC= 2 then 2:3  total 3 points for proc 1
        !            0 then 1:1  total 1 points for proc 2
        !            0 then 1:1  total 1 points for proc 3
        !            3 then 2:4  total 4 points for proc 4
        !--------------------------------------
        !                              9 points for all procs 1 to 4
        !
        TEMP1_INT_SND(I+1)=EGID(I_edgeline_ele(I))          !local element number converted to global element number
        TEMP2_INT_SND(I+1)=NGID(I_edgeline_ele_edge_ia(I))  !local ia node converted to global node number
        TEMP3_INT_SND(I+1)=NGID(I_edgeline_ele_edge_ib(I))  !local ib node converted to global node number
    ENDDO

   ENDIF
   
   
   !Gather element information of local procs
   CALL MPI_GATHERV(             TEMP1_INT_SND,  &  !sending buffer
                                      NEL_LOC1,  &  !number of points to send
                                   MPI_INTEGER,  &  !type of sending array
                                 TEMP1_INT_RCV,  &  !receiving buffer
                                 COUNTS_EL_ELE,  &  !number of points to receive
                          DISPLACEMENTS_EL_ELE,  &  !locations in TEMP2 to receive
                                   MPI_INTEGER,  &  !type of receiving buffer
                                             0,  &  !receiving proc
                                MPI_COMM_WORLD,  &  !communicator
                                          IERR) 
   
   
   !Gather ia node information of local procs
   CALL MPI_GATHERV(             TEMP2_INT_SND,  &  !sending buffer
                                      NEL_LOC1,  &  !number of points to send
                                   MPI_INTEGER,  &  !type of sending array
                                 TEMP2_INT_RCV,  &  !receiving buffer
                                 COUNTS_EL_ELE,  &  !number of points to receive
                          DISPLACEMENTS_EL_ELE,  &  !locations in TEMP2 to receive
                                   MPI_INTEGER,  &  !type of receiving buffer
                                             0,  &  !receiving proc
                                MPI_COMM_WORLD,  &  !communicator
                                          IERR) 
   
   
   !Gather ib node information of local procs
   CALL MPI_GATHERV(             TEMP3_INT_SND,  &  !sending buffer
                                      NEL_LOC1,  &  !number of points to send
                                   MPI_INTEGER,  &  !type of sending array
                                 TEMP3_INT_RCV,  &  !receiving buffer
                                 COUNTS_EL_ELE,  &  !number of points to receive
                          DISPLACEMENTS_EL_ELE,  &  !locations in TEMP2 to receive
                                   MPI_INTEGER,  &  !type of receiving buffer
                                             0,  &  !receiving proc
                                MPI_COMM_WORLD,  &  !communicator
                                          IERR) 
   
   !WRITE(*,*)'Here2 Here2'

   !
   !gather salt flux, saltwater flux, freshwater flux and total volume flux
   !
   
   
   NELWC=(NEL_LOC+1)*KBM1 !at least send a dummy vaue in (by adding 1)
                          !gatgerv may not work empty array in children (when NEL_LOC==0)
   !sendbuf cnt sndtype   rcvbuf rcvcnt rcvtype rcvroot COMM   IERR
   !all processors send their NEL_LOC to master (0)'s COUNTS_FLX array
   
   CALL MPI_GATHER(NELWC,  1,MPI_INTEGER,COUNTS_FLX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR) !gatger NELWC to COUNTS_FLX
   
   !master allocates enough memory to hold the values to gather from everyone
   !e.g. for NEL_LOC= 2,0, 0,3 for four procs, we have NELWC=COUNTS_FLX=(3*KBM1,1*KBM1, 1*KBM1,4*KBM1)
   !    with KBM1=10, we have COUNTS_FLX=(30,10,10,40)
   !==>DISPLACEMENTS_FLX(1) = 1
   !                    (2) = 1+30 = 31
   !                    (3) = 31+10 = 41
   !                    (4) = 41+10 = 51
   !
   !SUM(COUNTS_FLX)=30+10+10+40 = 90
   !
   !==>TEMP1_RCV (1:30) is for  1st proc with 21:30 being truly useful values
   !            (31:40) is for 2nd proc with no useful values
   !            (41:50) is for 3rd proc wtth no useful values
   !            (51:90) is for 4th proc with 61:90 being truly useful values
   
   IF(MSR)THEN
        DISPLACEMENTS_FLX(1)  = 0   !index for values received from first proc 
        DO I=2,NPROCS               !indices for values received from other procs
           DISPLACEMENTS_FLX(I)=COUNTS_FLX(I-1)+DISPLACEMENTS_FLX(I-1)
        ENDDO                      
        
        ALLOCATE(TEMP1_RCV(1:SUM(COUNTS_FLX)))  !to receive salt flux
        ALLOCATE(TEMP2_RCV(1:SUM(COUNTS_FLX)))  !to receive advective salt water flux
        ALLOCATE(TEMP3_RCV(1:SUM(COUNTS_FLX)))  !to receive freshwaer flux
        ALLOCATE(TEMP4_RCV(1:SUM(COUNTS_FLX)))  !to receive total voluem (water) flux
        ALLOCATE(TEMP5_RCV(1:SUM(COUNTS_FLX)))  !to receive layer thickness (m)
		
		ALLOCATE(TEMP6_RCV(1:SUM(COUNTS_FLX)))  !to receive diffusive salt flux
		ALLOCATE(TEMP7_RCV(1:SUM(COUNTS_FLX)))  !to receive advective salt flux
                 
   ENDIF
   
   
   IF(NELWC>0)THEN
   
   !prepare to send flux to master proc
    ALLOCATE(TEMP1_SND(1:NELWC)); TEMP1_SND=0; !to send salt flux
    ALLOCATE(TEMP2_SND(1:NELWC)); TEMP2_SND=0; !to send advective salt water flux
    ALLOCATE(TEMP3_SND(1:NELWC)); TEMP3_SND=0; !to send freshwaer flux
    ALLOCATE(TEMP4_SND(1:NELWC)); TEMP4_SND=0; !to send total voluem (water) flux
    ALLOCATE(TEMP5_SND(1:NELWC)); TEMP5_SND=0; !to send thickness

    ALLOCATE(TEMP6_SND(1:NELWC)); TEMP6_SND=0; !to send diffusive salt flux
	ALLOCATE(TEMP7_SND(1:NELWC)); TEMP7_SND=0; !to send advective salt flux


    !e.g. if NEL_LOC=2,0,0,3 and KBM1=10 ==> NELWC=30,10,10,40
    !and TEMP1_SND(1:10)=0,                                         
    !            (11:20)= EDGELINE_FLUX_SALT(1,1:10)  of proc 1
    !            (21:30)= EDGELINE_FLUX_SALT(2,1:10)  of proc 1
    !
    !            (31:40)=0
    !
    !            (41:50)=0
    !
    !            (51:60)=0
    !            (61:70)= EDGELINE_FLUX_SALT(1,1:10) of proc 4
    !            (71:80)= EDGELINE_FLUX_SALT(2,1:10) of proc 4
    !            (81:90)= EDGELINE_FLUX_SALT(3,1:10) of proc 4
    
    !fill out the sending buffers
            
    DO I=1,NEL_LOC
       DO K=1,KBM1
        !Note I>=1, so the data are put in starting from KBM1+1 in the sending buffers
        !when master proc receives it, it should also start retrieving data from KBM1+1
        !e.g. when KBM1=10, then 11:20 is for I=1
        !                        21:30 is for I=2
        !
        TEMP1_SND(I*KBM1+K)=EDGELINE_FLUX_SALT(I,K)
        TEMP2_SND(I*KBM1+K)=EDGELINE_FLUX_ADVECTION_SALTWATER(I,K)
        TEMP3_SND(I*KBM1+K)=EDGELINE_FLUX_FRSHWATER(I,K)
        TEMP4_SND(I*KBM1+K)=EDGELINE_FLUX_WATER(I,K)
        TEMP5_SND(I*KBM1+K)=EDGELINE_FLUX_DZ(I,K)

		TEMP6_SND(I*KBM1+K)=EDGELINE_FLUX_DIFFUSION_SALT(I,K)
		TEMP7_SND(I*KBM1+K)=EDGELINE_FLUX_ADVECTION_SALT(I,K)		

       ENDDO
    ENDDO
   ENDIF
   
   
 
   !gather total salt flux
   CALL MPI_GATHERV(             TEMP1_SND,  &  !sending buffer
                                     NELWC,  &  !number of points to send
                                     MPI_F,  &  !type of sending array
                                 TEMP1_RCV,  &  !receiving buffer
                                COUNTS_FLX,  &  !number of points to receive
                         DISPLACEMENTS_FLX,  &  !locations in TEMP2 to receive
                                     MPI_F,  &  !type of receiving buffer
                                         0,  &  !receiving proc
                            MPI_COMM_WORLD,  &  !communicator
                                      IERR) 

   !gather diffsive salt flux
   CALL MPI_GATHERV(             TEMP6_SND,  &  !sending buffer
                                     NELWC,  &  !number of points to send
                                     MPI_F,  &  !type of sending array
                                 TEMP6_RCV,  &  !receiving buffer
                                COUNTS_FLX,  &  !number of points to receive
                         DISPLACEMENTS_FLX,  &  !locations in TEMP2 to receive
                                     MPI_F,  &  !type of receiving buffer
                                         0,  &  !receiving proc
                            MPI_COMM_WORLD,  &  !communicator
                                      IERR) 

   !gather advection salt flux
   CALL MPI_GATHERV(             TEMP7_SND,  &  !sending buffer
                                     NELWC,  &  !number of points to send
                                     MPI_F,  &  !type of sending array
                                 TEMP7_RCV,  &  !receiving buffer
                                COUNTS_FLX,  &  !number of points to receive
                         DISPLACEMENTS_FLX,  &  !locations in TEMP2 to receive
                                     MPI_F,  &  !type of receiving buffer
                                         0,  &  !receiving proc
                            MPI_COMM_WORLD,  &  !communicator
                                      IERR) 
   
   !gather advective saltwater volme flux
   CALL MPI_GATHERV(             TEMP2_SND,  &  !sending buffer
                                     NELWC,  &  !number of points to send
                                     MPI_F,  &  !type of sending array
                                 TEMP2_RCV,  &  !receiving buffer
                                COUNTS_FLX,  &  !number of points to receive
                         DISPLACEMENTS_FLX,  &  !locations in TEMP2 to receive
                                     MPI_F,  &  !type of receiving buffer
                                         0,  &  !receiving proc
                            MPI_COMM_WORLD,  &  !communicator
                                      IERR)    
   
   !gather freshwater volume flux
   CALL MPI_GATHERV(             TEMP3_SND,  &  !sending buffer
                                     NELWC,  &  !number of points to send
                                     MPI_F,  &  !type of sending array
                                 TEMP3_RCV,  &  !receiving buffer
                                COUNTS_FLX,  &  !number of points to receive
                         DISPLACEMENTS_FLX,  &  !locations in TEMP2 to receive
                                     MPI_F,  &  !type of receiving buffer
                                         0,  &  !receiving proc
                            MPI_COMM_WORLD,  &  !communicator
                                      IERR)    
   
   !gather total volume flux (=advective salt water + freshater volume flux)
   CALL MPI_GATHERV(             TEMP4_SND,  &  !sending buffer
                                     NELWC,  &  !number of points to send
                                     MPI_F,  &  !type of sending array
                                 TEMP4_RCV,  &  !receiving buffer
                                COUNTS_FLX,  &  !number of points to receive
                         DISPLACEMENTS_FLX,  &  !locations in TEMP2 to receive
                                     MPI_F,  &  !type of receiving buffer
                                         0,  &  !receiving proc
                            MPI_COMM_WORLD,  &  !communicator
                                      IERR)

   !gather thickness
   CALL MPI_GATHERV(             TEMP5_SND,  &  !sending buffer
                                     NELWC,  &  !number of points to send
                                     MPI_F,  &  !type of sending array
                                 TEMP5_RCV,  &  !receiving buffer
                                COUNTS_FLX,  &  !number of points to receive
                         DISPLACEMENTS_FLX,  &  !locations in TEMP2 to receive
                                     MPI_F,  &  !type of receiving buffer
                                         0,  &  !receiving proc
                            MPI_COMM_WORLD,  &  !communicator
                                      IERR)

   DEALLOCATE(TEMP1_INT_SND)
   DEALLOCATE(TEMP2_INT_SND)
   DEALLOCATE(TEMP3_INT_SND)

   DEALLOCATE(TEMP1_SND)
   DEALLOCATE(TEMP2_SND)
   DEALLOCATE(TEMP3_SND)
   DEALLOCATE(TEMP4_SND)
   DEALLOCATE(TEMP5_SND)
   DEALLOCATE(TEMP6_SND)   
   DEALLOCATE(TEMP7_SND)      

   !master proc puts gathered results into the global array
   IF(MSR)THEN
       
        !ALLOCATE global edge line arrays 
        !I_edgeline_ele_GL(1:NEL_GL)
        !I_edgeline_ele_edge_ia_GL(1:NEL_GL)
        !I_edgeline_ele_edge_ib_GL(1:NEL_GL)
        !I_edgeline_flux_sign_GL(1:NEL_GL)
  
        !loop through the selected global TCE edges, find them in the received buffers
        !and the put the values in the following corresponding global output array for edge line fluxes
   
        !EDGELINE_FLUX_SALT_GL(1:NEL_GL,1:KBM1)        !salt flux through the TCE edge
        !EDGELINE_FLUX_SALTWATER_GL(1:NEL_GL,1:KBM1)   !saltwater flux through the TCE edge
        !EDGELINE_FLUX_FRSHWATER_GL(1:NEL_GL,1:KBM1)   !freshwater flux through the TCE edge
        !EDGELINE_FLUX_WATER_GL(1:NEL_GL,1:KBM1)       !total (salt water + fresh water) volume flux through the TCE edge
        !EDGELINE_FLUX_DZ_GL(1:NEL_GL,1:KBM1)          !thickness of layers (m)


        IF(NEL_GL>0)THEN

            DO I=1,NEL_GL  !Loop through global TCE edges
                    
                DO J=1,NPROCS !Loop through all procs
                    
                    N_NEL_LOC=COUNTS_EL_ELE(J)-1            !number of local TCE edges received from proc J
                    I_DISPLACEMENT=DISPLACEMENTS_EL_ELE(J)  !starting indices in receiving buffer
                
                    N_FLX_LOC=COUNTS_FLX(J)-KBM1            !number of local TCE edge flux columns received from proc J
                    I_DISPLACEMENT_FLX=DISPLACEMENTS_FLX(J) !starting indices in receiving buffer for flux
   
                    DO II=1,N_NEL_LOC    !loop all received local TCE edges for proc J
                    
                        !get global element and node numbers of the II'th TCE edge for local TCE edge line in proc J
                    
                        !e.g. NEL_LOC = 2,0,0,3 ==> COUNTS_EL_ELE=(3,1,1,4)==>N_NEL_LOC=(2,0,0,3)
                        !==> I_DISPLACEMENT =(1,4,5,6)
                        !==> corresponding chunks in receiving buffers are (1:3) (4:4) (5:5) (6:9)
                        !                                                     2    0     0     3
                        ! I_DISPLACEMENT(J)+1:I_DISPLACEMENT(J)+N_NEL_LOC =(2:3) ()    ()    (7:9)  

                        IE_EL_GL=TEMP1_INT_RCV(I_DISPLACEMENT+II+1) !global element number of II'th TCE edge 
                        ia_EL_GL=TEMP2_INT_RCV(I_DISPLACEMENT+II+1) !global node ia number of II'th TCE edge  
                        ib_EL_GL=TEMP3_INT_RCV(I_DISPLACEMENT+II+1) !global node ib number of II'th TCE edge 
    
!                            WRITE(*,*)'Here6 '
!                            WRITE(*,*)'Proc J=',J,' of ',NPROCS
!                            WRITE(*,*)'loc edgeline index II=',II,' of ', N_NEL_LOC
!                            WRITE(*,*)'global edgeline element IE_GL of I=',I, ' is ',I_edgeline_ele_GL(I)
!                            WRITE(*,*)'received edgeline element global ID IE_EL_GL of II=',II,' is: ',IE_EL_GL
!                            WRITE(*,*)'global ia ib are: ia=',I_edgeline_ele_edge_ia_GL(I),'ib=',I_edgeline_ele_edge_ib_GL(I)
!                            WRITE(*,*)'local ia ib are : ia=',ia_EL_GL,'ib=',ib_EL_GL
!                            WRITE(*,*)'I_DISPLACEMENT of proc J=',I_DISPLACEMENT
!                            WRITE(*,*)'TEMP1_INT_RCV=',(TEMP1_INT_RCV(K),K=1,SUM(COUNTS_EL_ELE))
!                            WRITE(*,*)'TEMP2_INT_RCV=',(TEMP2_INT_RCV(k),K=1,SUM(COUNTS_EL_ELE))
!                            WRITE(*,*)'TEMP3_INT_RCV=',(TEMP3_INT_RCV(k),K=1,SUM(COUNTS_EL_ELE))
!                            WRITE(*,*)'SUM(COUNTS_FLX)=',SUM(COUNTS_FLX)
!                            WRITE(*,*)'I_DISPLACEMENT_FLX of proc J=',I_DISPLACEMENT_FLX
!                            WRITE(*,*)'TEMP1_RCV=',(TEMP1_RCV(K),K=1,SUM(COUNTS_FLX))
!                            DO K=1,KBM1
!                               WRITE(*,*)'FLX index of II,K=',K,' is ',I_DISPLACEMENT_FLX + II*KBM1+K
!                            ENDDO

                        IF (IE_EL_GL==I_edgeline_ele_GL(I) .AND.         &
                            ia_EL_GL==I_edgeline_ele_edge_ia_GL(I) .AND. &
                            ib_EL_GL==I_edgeline_ele_edge_ib_GL(I))THEN
                        
                           !find a correspondence between global TCE edge I and local TCE edge II
                           !for their global element numbers and global node numbers match
                        
                            DO K=1,KBM1
                                !e.g. if NEL_LOC=2,0,0,3, KBM1=10
                                !==>      I_DISPLACEMENT_FLX=(1,31,41,51)
                                !
                                !         (II-1)*KBM1+K is of range (II-1)*KBM1+1 : (II*KBM1)  ==> 1:10 for II =1, K=1:KBM1
                                !                                                                 11:20 for II =2, K=1,KBM1
                                !                                                                 21:30 for II =3, K=1,KBM1
                                !
                                !       ==>I_DISPLACEMENT_FLX + II*KBM1+K-1 which is of range     11:20 for II=1, K=1:KBM1
                                !                                                                 21:30 for II=2, K=1:KBM1
                                !                                                            
                                !                                                   
                                !for J=2, we have NEW_LOC=0
                                !for J=4, we have NEW_LOC=3
                                !       ==> I_DISPLACEMENT_FLX = 51
                                !       ==> I_DISPLACEMENT_FLX + II*KBM1+K-1 => 51+10+K-1 = 61:70 for II=1 K=1:KBM1
                                !                                                           71:80 for II=2,K=1:KBM1
                                

                                EDGELINE_FLUX_SALT_GL(I,K)                =TEMP1_RCV(I_DISPLACEMENT_FLX + II*KBM1+K) !-1) !
                                EDGELINE_FLUX_ADVECTION_SALTWATER_GL(I,K) =TEMP2_RCV(I_DISPLACEMENT_FLX + II*KBM1+K) !-1)
                                EDGELINE_FLUX_FRSHWATER_GL(I,K)           =TEMP3_RCV(I_DISPLACEMENT_FLX + II*KBM1+K) !-1)
                                EDGELINE_FLUX_WATER_GL(I,K)               =TEMP4_RCV(I_DISPLACEMENT_FLX + II*KBM1+K) !-1)
                                EDGELINE_FLUX_DZ_GL(I,K)                  =TEMP5_RCV(I_DISPLACEMENT_FLX + II*KBM1+K) !-1)
                                EDGELINE_FLUX_DIFFUSION_SALT_GL(I,K)      =TEMP6_RCV(I_DISPLACEMENT_FLX + II*KBM1+K) !-1) !																
                                EDGELINE_FLUX_ADVECTION_SALT_GL(I,K)      =TEMP7_RCV(I_DISPLACEMENT_FLX + II*KBM1+K) !-1) !																								
                                
                            ENDDO

!                            WRITE(*,*)'Here7 '
!                            WRITE(*,*)'Proc J=',J,' of ',NPROCS
!                            WRITE(*,*)'loc edgeline index II=',II,' of ', N_NEL_LOC
!                            WRITE(*,*)'global edgeline element IE_GL of I=',I, ' is ',I_edgeline_ele_GL(I)
!                            WRITE(*,*)'received edgeline element global ID IE_EL_GL of II=',II,' is: ',IE_EL_GL
!
!                            WRITE(*,*)'global ia ib are: ia=',I_edgeline_ele_edge_ia_GL(I),'ib=',I_edgeline_ele_edge_ib_GL(I)
!                            WRITE(*,*)'local ia ib are : ia=',ia_EL_GL,'ib=',ib_EL_GL
!                            WRITE(*,*)'I_DISPLACEMENT of proc J=',I_DISPLACEMENT
!                            WRITE(*,*)'TEMP1_INT_RCV=',(TEMP1_INT_RCV(K),K=1,SUM(COUNTS_EL_ELE))
!                            WRITE(*,*)'TEMP2_INT_RCV=',(TEMP2_INT_RCV(k),K=1,SUM(COUNTS_EL_ELE))
!                            WRITE(*,*)'TEMP3_INT_RCV=',(TEMP3_INT_RCV(k),K=1,SUM(COUNTS_EL_ELE))
!                            WRITE(*,*)'SUM(COUNTS_FLX)=',SUM(COUNTS_FLX)
!                            WRITE(*,*)'I_DISPLACEMENT_FLX of proc J=',I_DISPLACEMENT_FLX
!                            WRITE(*,*)'TEMP1_RCV=',(TEMP1_RCV(K),K=1,SUM(COUNTS_FLX))
!                            DO K=1,KBM1
!                               WRITE(*,*)'FLX index of II,K=',K,' is ',I_DISPLACEMENT_FLX + II*KBM1+K
!                            ENDDO
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
            
        ENDIF

        !Master proc deallocates these gathering variables
        DEALLOCATE(TEMP1_RCV)
        DEALLOCATE(TEMP2_RCV)
        DEALLOCATE(TEMP3_RCV)
        DEALLOCATE(TEMP4_RCV)
        DEALLOCATE(TEMP5_RCV)
        DEALLOCATE(TEMP6_RCV)
        DEALLOCATE(TEMP7_RCV)

        DEALLOCATE(TEMP1_INT_RCV)
        DEALLOCATE(TEMP2_INT_RCV)
        DEALLOCATE(TEMP3_INT_RCV)

        !Master proc deallocates COUTNS_FLX,COUNTS_EL_ELE

        DEALLOCATE(COUNTS_FLX)
        DEALLOCATE(COUNTS_EL_ELE)
        DEALLOCATE(DISPLACEMENTS_EL_ELE)
        DEALLOCATE(DISPLACEMENTS_FLX)

   ENDIF
END IF
# endif

IF (SERIAL)THEN
    !Simply retrieve values into global array
    IF(NEL_GL>0)THEN
        EDGELINE_FLUX_SALT_GL(1:NEL_GL,:)=EDGELINE_FLUX_SALT(1:NEL_LOC,:)
        EDGELINE_FLUX_ADVECTION_SALTWATER_GL(1:NEL_GL,:)=EDGELINE_FLUX_ADVECTION_SALTWATER(1:NEL_LOC,:)
        EDGELINE_FLUX_FRSHWATER_GL(1:NEL_GL,:)=EDGELINE_FLUX_FRSHWATER(1:NEL_LOC,:)
        EDGELINE_FLUX_WATER_GL(1:NEL_GL,:)=EDGELINE_FLUX_WATER(1:NEL_LOC,:)    
        EDGELINE_FLUX_DZ_GL(1:NEL_GL,:)=EDGELINE_FLUX_DZ(1:NEL_LOC,:)
        EDGELINE_FLUX_DIFFUSION_SALT_GL(1:NEL_GL,:)=EDGELINE_FLUX_DIFFUSION_SALT(1:NEL_LOC,:)
        EDGELINE_FLUX_ADVECTION_SALT_GL(1:NEL_GL,:)=EDGELINE_FLUX_ADVECTION_SALT(1:NEL_LOC,:)
    ENDIF
ENDIF
   
   !WRITE(*,*)'Here4 Here4'
  !Master writes the results to a file
   IF(MSR)THEN
        IF(NEL_GL>0)THEN
            DO I=1,NEL_GL                                       
               WRITE(107,'(F12.5,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,4000(F12.2,1X))') THOUR,     &
                                      I_edgeline_transect_ID_GL(I),                            &
                                      I_edgeline_ele_GL(I),                                    &
                                      I_edgeline_ele_edge_ia_GL(I),                            &
                                      I_edgeline_ele_edge_ib_GL(I),                            &
                                      I_edgeline_flux_sign_GL(I),                              &
                                      (EDGELINE_FLUX_SALT_GL(I,K),     K=1,KBM1),              &
                                      (EDGELINE_FLUX_ADVECTION_SALTWATER_GL(I,K),K=1,KBM1),    &
                                      (EDGELINE_FLUX_FRSHWATER_GL(I,K),K=1,KBM1),              &
                                      (EDGELINE_FLUX_WATER_GL(I,K),    K=1,KBM1),              &
                                      (EDGELINE_FLUX_DZ_GL(I,K),K=1,KBM1),                     &
                                      (EDGELINE_FLUX_ADVECTION_SALT_GL(I,K),K=1,KBM1),         &
                                      (EDGELINE_FLUX_DIFFUSION_SALT_GL(I,K),K=1,KBM1)
            ENDDO
        ENDIF
   ENDIF
   RETURN
   END SUBROUTINE WRITE_EDGELINE_FLUX
!==============================================================================|

!============================================================================================!
!  Deallocate variables for edgeline flux calculations                                       !
!============================================================================================!
   
   SUBROUTINE DEALLOC_EDGELINE_FLUX
   USE ALL_VARS
   USE CONTROL
   IMPLICIT NONE
   INTEGER :: i,k
   
   IF(NEL_GL>0)THEN
   
        !DEALLOCATE GLOBAL ARRAY
        DEALLOCATE(I_edgeline_ele_GL)
        DEALLOCATE(I_edgeline_ele_edge_ia_GL)
        DEALLOCATE(I_edgeline_ele_edge_ib_GL)
        DEALLOCATE(I_edgeline_flux_sign_GL)
        DEALLOCATE(I_edgeline_transect_ID_GL)
   
        DEALLOCATE(EDGELINE_FLUX_SALT_GL)  
        DEALLOCATE(EDGELINE_FLUX_DIFFUSION_SALT_GL)  
        DEALLOCATE(EDGELINE_FLUX_ADVECTION_SALT_GL)  
        DEALLOCATE(EDGELINE_FLUX_ADVECTION_SALTWATER_GL)
        DEALLOCATE(EDGELINE_FLUX_FRSHWATER_GL)
        DEALLOCATE(EDGELINE_FLUX_WATER_GL)
        DEALLOCATE(EDGELINE_FLUX_DZ_GL)
        
        IF(NEL_LOC>0)THEN

            !DEALLOCATE LOCAL ARRAY
            DEALLOCATE(I_edgeline_ele)
            DEALLOCATE(I_edgeline_ele_edge_ia)
            DEALLOCATE(I_edgeline_ele_edge_ib)
            DEALLOCATE(I_edgeline_flux_sign)
   
            DEALLOCATE(EDGELINE_FLUX_SALT)
            DEALLOCATE(EDGELINE_FLUX_DIFFUSION_SALT)
            DEALLOCATE(EDGELINE_FLUX_ADVECTION_SALT)
            DEALLOCATE(EDGELINE_FLUX_ADVECTION_SALTWATER)
            DEALLOCATE(EDGELINE_FLUX_FRSHWATER)
            DEALLOCATE(EDGELINE_FLUX_WATER)

        ENDIF

        IF(MSR)THEN !Master closes TCE flux output file
            CLOSE(107)
        ENDIF

   ENDIF

   RETURN
   END SUBROUTINE DEALLOC_EDGELINE_FLUX

!==============================================================================|

# endif

END MODULE MOD_OUTPUT_FLUX
