
        MODULE MOD_BCS
		
	      USE MOD_TYPES , ONLY:   & 
			!GMAP, 		&	! grid mapping type
			!NSIZE, 	&	! size
			!LOC_2_GL, 	&	! conversion from local to global grid for a mapping
			!COMM, 		&	! communication world type
			!NSND, 		&	! sending buffer size
			!NRCV, 		&	! receiving buffer size
			!RCPT, 		&	! pointer in receiving buffer
			!SNDP,      &    ! sending buffer
			!RCVP, 		&	! receiving buffer
			!MLTP, 		&	! muliplicity of each processor
			BC!, 		&	!boundary type
			!NTIMES,	&	!
			!TIMES,		&	!
			!LABEL

           USE MOD_PREC, ONLY: SP
		   
		   USE MOD_CONTROL, ONLY: PAR
		   
           IMPLICIT NONE
           SAVE
        !----------------boundary conditions: fresh water discharge-------------------------!
           INTEGER, ALLOCATABLE  :: INODEQ(:)        !!LOCAL FRESH WATER INFLOW NODES
           INTEGER, ALLOCATABLE  :: ICELLQ(:)        !!LOCAL FRESH WATER INFLOW ELEMENTS 
           INTEGER, ALLOCATABLE  :: RIV_GL2LOC(:)    !!GLOBAL TO LOCAL MAP OF FW NODES
           INTEGER, ALLOCATABLE  :: INOPNT(:)        !!LOCAL NON-POINT SOURCE NODES
		   
           TYPE(BC)              :: QBC_TM           !!TIME MAP FOR RIVER DISCHARGE DATA 
           TYPE(BC)              :: PNT_TM           !!TIME MAP FOR NON POINT SOURCE DATA
		   
		   !WLong moved this to mod_bcmap.F
           !TYPE(BC)              :: NUT_TM           !!TIME MAPPING FOR NUTRIENT OBC  !should be 
		   
		   !WLong these are never used for we are not soving hydrodynamics 
           !REAL(SP), ALLOCATABLE     ::  VQDIST(:,:)     !!DISCHARGE VERTICAL DISTRIBUTION
           !REAL(SP), ALLOCATABLE     ::   DQDIS(:,:)     !!WATER FLUX DISCHARGE DATA
		   
           REAL(SP), ALLOCATABLE     ::   DWDIS(:,:,:)   !!WATER QUALITY DISCHARGE DATA
           REAL(SP), ALLOCATABLE     ::    WDIS(:,:)     !!CURRENT TIME FRESH WATER QUALITY 
		   REAL(SP), ALLOCATABLE     :: WVQDIST(:,:)     !!DISCHARGE VERTICAL DISTRIBUTION for point source


           REAL(SP), ALLOCATABLE     ::  PQDIS(:)        !!Current DISCHARGE at non-point source node
           REAL(SP), ALLOCATABLE     :: PDQDIS(:,:)      !!DISCHARGE at non-point source node
           REAL(SP), ALLOCATABLE     :: PWQDIS(:,:)      !!Current WATER QUALITY at non-point source node
		   
           REAL(SP), ALLOCATABLE     :: PDWQDIS(:,:,:)   !!WATER QUALITY DATA at non-point source node		   
		   
		   !Wlong moved this to mod_bcmap.F
		   !REAL(SP), ALLOCATABLE     :: WQOBC(:,:,:,:)   !!WATER QUALITY DATA AT BOUNDARY	!WLong should not be here
        !===================================================================================!
		   
		   INTEGER, ALLOCATABLE :: N_ICELLQ(:,:)     !!node number of the edge for element edge source
		   
	CONTAINS 
			
		!Subroutine BCS_ALLOC()
		!Subroutine BCS_DEALLOC()
		   
		   SUBROUTINE BCS_ALLOC

				!INODEQ(:)        !!LOCAL FRESH WATER INFLOW NODES		  !allocated in bcs_force 
				!ICELLQ(:)        !!LOCAL FRESH WATER INFLOW ELEMENTS 				!in bcs_force	
				!RIV_GL2LOC(:)    !!GLOBAL TO LOCAL MAP OF FW NODES					!in bcs_force	
				!INOPNT(:)        !!LOCAL NON-POINT SOURCE NODES					!in bcs_force	
				
				!DWDIS(:,:,:)     !!WATER QUALITY DISCHARGE DATA					!in bcs_force	
				
				!WDIS(:,:)        !!CURRENT TIME FRESH WATER QUALITY				!in bcond_wqm.F	
				!WVQDIST(:,:)     !!DISCHARGE VERTICAL DISTRIBUTION for point source!in bcs_force	

				!PQDIS(:)         !!Current DISCHARGE at non-point source node		!in adv_wqm.F	
				!PDQDIS(:,:)      !!DISCHARGE at non-point source node				!in bcs_force.F 
				!PWQDIS(:,:)      !!Current WATER QUALITY at non-point source node  !in adv_wqm.F	

				!PDWQDIS(:,:,:)   !!WATER QUALITY DATA at non-point source node		!in bcs_force.F 

		   		!QBC_TM%TIMES(:)  !in bcs_force.F
				!PNT_TM%TIMES(:)  !in bcs_force.F
		   
				!!VQDIST(:,:)      !!DISCHARGE VERTICAL DISTRIBUTION				!never used		
				!!DQDIS(:,:)       !!WATER FLUX DISCHARGE DATA						!never used		

			!===================================================================================!
		   
		   END SUBROUTINE BCS_ALLOC
		   
		   SUBROUTINE BCS_DEALLOC
		   
				IF(ALLOCATED(INODEQ))DEALLOCATE(INODEQ)				
				IF(ALLOCATED(ICELLQ))DEALLOCATE(ICELLQ)				
#if defined (MULTIPROCESSOR)
				IF(PAR)THEN
					IF(ALLOCATED(RIV_GL2LOC))DEALLOCATE(RIV_GL2LOC)	 	
				ENDIF
#endif				
				IF(ALLOCATED(INOPNT))DEALLOCATE(INOPNT)				
				IF(ALLOCATED(DWDIS))DEALLOCATE(DWDIS)				
				IF(ALLOCATED(WVQDIST))DEALLOCATE(WVQDIST)			
				IF(ALLOCATED(PDQDIS))DEALLOCATE(PDQDIS)				
				IF(ALLOCATED(PDWQDIS))DEALLOCATE(PDWQDIS)			
				
				!IF(ALLOCATED(WDIS))DEALLOCATE(WDIS)					!locally deallocated in bcond_wqm.F
				!IF(ALLOCATED(PQDIS))DEALLOCATE(PQDIS)					!locally deallocated in adv_wqm.F
				!IF(ALLOCATED(PWQDIS))DEALLOCATE(PWQDIS)				!locally deallocated in adv_wqm.F
				
				!!IF(ALLOCATED(VQDIST))DEALLOCATE(VQDIST)				!never used for not solving flow equation
				!!IF(ALLOCATED(DQDIS))DEALLOCATE(DQDIS)				!never used for not solving flow equation

				IF(ALLOCATED(QBC_TM%TIMES)) DEALLOCATE(QBC_TM%TIMES)
				IF(ALLOCATED(PNT_TM%TIMES)) DEALLOCATE(PNT_TM%TIMES)
                
                
                !LB: these were missing:
                IF(ALLOCATED(N_ICELLQ))     DEALLOCATE(N_ICELLQ)
				
		   END SUBROUTINE BCS_DEALLOC
		   
        END MODULE  MOD_BCS
