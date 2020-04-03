MODULE MOD_ZOOP

	USE MOD_PREC, ONLY: SP
	USE MOD_LIMS, ONLY: MTLOC, KBM1

	USE MOD_SIZES, ONLY: MGL
	
	USE MOD_WQM, ONLY:			&!
			    ZFOFN,          &!
                ZOOFN!,         &!
	
 !Zooplankton parameters
				
			INTEGER  :: 	KTBSZ,		&
							KTBLZ
							
			REAL(SP),DIMENSION(366) ::  TVPRSZ, &
									TVPRLZ
										
                                 
             REAL(SP) :: CTSZ,    &  !carbon threshold for grazing (gC/m^3) for microzooplankton
                     CTLZ,    &  !carbon threshold for grazing (gC/m^3) for mesozooplankton
                     KHCSZ,   &  !prey density at which grazing is halved (gC/m^3) for microzooplankton
                     KHCLZ,   &  !prey density at which grazing is halved (gC/m^3) for mesozooplankton
                     MZEROSZ, &  !mortality at zero dissolved oxygen (1/day) for microzooplankton
                     MZEROLZ     !mortality at zero dissolved oxygen (1/day) for mesozooplankton

             REAL(SP) :: UB1SZ, &   !Utilization of algal group 1 by microzooplankton, range [0,1]
                     UB2SZ, &   !Utilization of algal group 2 by microzooplankton, range [0,1]
                     UB3SZ, &   !Utilization of algal group 3 by microzooplankton, range [0,1]
                     UB1LZ, &   !Utilization of algal group 1 by mesozooplankton, range [0,1]
                     UB2LZ, &   !Utilization of algal group 2 by mesozooplankton, range [0,1]
                     UB3LZ, &   !Utilization of algal group 3 by mesozooplankton, range [0,1]                         
                     UDSZ,  &   !Never used, replaced by ULDSZ and URDSZ below 
                     ULDSZ, &   !Utilization of labile dissolved organic carbon by microzooplankton, range [0,1]
                     URDSZ, &   !Utilization of refractory dissolved organic carbon by microzooplankton, range [0,1]
                     ULPSZ, &   !microzooplankton utilization of labile particulate organic carbon, range [0,1]
                     URPSZ, &   !microzooplankton utilization of refractory particulate organic carbon, range [0,1]
                     ULLZ,  &   !mesozooplankton utilization of labile particulate organic carbon, range [0,1]
                     URLZ,  &   !mesozooplankton utilization of refractory particulate organic carbon, range [0,1]
                     USZLZ, &   !mesozooplankton utilization of microzooplankton, range [0,1], (~1.0)
                     TRSZ,  &   !reference temperature for microzooplankton metabolism (degC)
                     TRLZ,  &   !reference temperature for mesozooplankton metabolism (degC)
                     DOCRITSZ, &!dissolved oxygen below which mortality occurs for microzooplankton (gO2/m^3)
                     DOCRITLZ, &!dissolved oxygen below which mortality occurs for mesozooplankton (gO2/m^3)
                     ANCSZ,  &  !microzooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     ANCLZ,  &  !mesozooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     APCSZ,  &  !microzooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
                     APCLZ,  &  !mesozooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
                     AOCRSZ, &  !ratio of oxygen consumed to microzooplankton carbon metabolized (gO2/gC) (~2.67)
                     AOCRLZ, &  !ratio of oxygen consumed to mesozooplankton carbon metabolized (gO2/gC) (~2.67)
                     FRSASZ, &  !fraction of microzooplankton silica recycled to dissolved silica pool, range [0,1]
                     FRSALZ     !fraction of mesozooplankton silica recycled to dissolved silica pool, range [0,1]
                                 
             REAL(SP) ::&  !FDOCSZ,&    !fraction of microzooplankton carbon released to NCDOC1, range [0,1]  ! Have 3 coefficients for to support the NCDOC vars, at this point zooplankton and algae only create non-colored DOC
						  ! FDOC2SZ,&    !fraction of microzooplankton carbon released to NCDOC2, range [0,1]
					      ! FDOC3SZ,&    !fraction of microzooplankton carbon released to NCDOC3, range [0,1]
					 
                     FLPOCSZ,&    !fraction of microzooplankton carbon released to LPOC, range [0,1]
                     FRPOCSZ,&    !fraction of microzooplankton carbon released to RPOC, range [0,1]
					 
                         !  FDON1SZ,&    !fraction of microzooplankton nitrogen released to NCDON1, range [0,1] ! B Clark add in 3 coefficients to support 3 G classes of NCDOM
                          ! FDON2SZ,&    !fraction of microzooplankton nitrogen released to NCDON2, range [0,1]
					     !  FDON3SZ,&    !fraction of microzooplankton nitrogen released to NCDON3, range [0,1]
					 
                     FLPONSZ,&    !fraction of microzooplankton nitrogen released to LPON, range [0,1]
                     FRPONSZ,&    !fraction of microzooplankton nitrogen released to RPON, range [0,1]
					 
                         ! FDOP1SZ,&    !fraction of microzooplankton phosphorus released to NCDOP1 range [0,1]    ! B Clark add in 3 coefficients to support 3 G classes of NCDOM
                         ! FDOP2SZ,&    !fraction of microzooplankton phosphorus released to NCDOP2 range [0,1]  
					      !FDOP3SZ,&    !fraction of microzooplankton phosphorus released to NCDOP3 range [0,1]  
					 
                     FLPOPSZ,&    !fraction of microzooplankton phosphorus released to LPOP range [0,1]    
                     FRPOPSZ,&    !fraction of microzooplankton phosphorus released to RPOP range [0,1]   
					 
                     FNH4SZ, &    !fraction of microzooplankton nitrogen recycled to DIN as NH4, range [0,1]  ! inorganic nutrient production
                     FPO4SZ, &    !fraction of microzooplankton phosphorus recycled to DIP as PO4, range [0,1]
                        
                        ! FDOC1LZ,&    !fraction of mesozooplankton carbon released to NCDOC1, range [0,1]   ! B Clark add in 3 coefficients to support 3 G classes of NCDOM
                        ! FDOC2LZ,&    !fraction of mesozooplankton carbon released to NCDOC2, range [0,1]
					    ! FDOC3LZ,&    !fraction of mesozooplankton carbon released to NCDOC3, range [0,1]
					  
                     FLPOCLZ,&    !fraction of mesozooplankton carbon released to LPOC, range [0,1]
                     FRPOCLZ,&    !fraction of mesozooplankton carbon released to RPOC, range [0,1]
					 
                        ! FDON1LZ,&    !fraction of mesozooplankton nitrogen released to NCDON1, range [0,1] ! B Clark add in 3 coefficients to support 3 G classes of NCDOM
                        ! FDON2LZ,&    !fraction of mesozooplankton nitrogen released to NCDON2, range [0,1]
					    ! FDON3LZ,&    !fraction of mesozooplankton nitrogen released to NCDON3, range [0,1]
					 
                     FLPONLZ,&    !fraction of mesozooplankton nitrogen released to LPON, range [0,1]
                     FRPONLZ,&    !fraction of mesozooplankton nitrogen released to RPON, range [0,1]
					 
                         !FDOP1LZ,&    !fraction of mesozooplankton phosphorus released to NCDOP1 range [0,1]    ! B Clark add in 3 coefficients to support 3 G classes of NCDOM
                        ! FDOP2LZ,&    !fraction of mesozooplankton phosphorus released to NCDOP2 range [0,1]
					     !FDOP3LZ,&    !fraction of mesozooplankton phosphorus released to NCDOP3 range [0,1] 
					 
                     FLPOPLZ,&    !fraction of mesozooplankton phosphorus released to LPOP range [0,1]    
                     FRPOPLZ,&    !fraction of mesozooplankton phosphorus released to RPOP range [0,1]    
					 
                     FNH4LZ, &    !fraction of mesozooplankton nitrogen recycled to DIN as NH4, range [0,1]  ! Inorganics
                     FPO4LZ, &    !fraction of mesozooplankton phosphorus recycled to DIP as PO4, range [0,1]
					 
                     FUREASZ,&    !Never used !!! 					 
                     FUREALZ      !Never used!!!

             REAL(SP),ALLOCATABLE::  B1ASZ(:,:),   & 
                                 B2ASZ(:,:),   &
                                 B3ASZ(:,:),   &
                                 LPOCASZ(:,:), &
                                 RPOCASZ(:,:), &
                                 PRASZ(:,:),   &
                                 B1ALZ(:,:),   &
                                 B2ALZ(:,:),   &
                                 B3ALZ(:,:),   &
                                 SZALZ(:,:),   &
                                 LPOCALZ(:,:), &
                                 RPOCALZ(:,:), &
                                 PRALZ(:,:),   &
                                 CLSZ(:,:),    & 
                                 CLLZ(:,:),    &
                                 RSZ(:,:),     &
                                 RLZ(:,:),     &
                                 RMAXSZ(:,:),  &
                                 RMAXLZ(:,:),  &
                                 BMSZ(:,:),    &
                                 BMLZ(:,:),    &
                                 BMRSZ(:,:),   &
                                 BMRLZ(:,:),   &
                                 MSZ(:,:),     &
                                 MLZ(:,:),     &
                                 PRSZLZ(:,:),  &	!Predation of micro-zoop by meso-zoop
                                 GSZ(:,:),     &
                                 GLZ(:,:),     &
                                 ESZ(:,:),     &
                                 ELZ(:,:),     &
                                 RFSZ(:,:),    &
                                 RFLZ(:,:),    &
                                 PRSZ(:,:),    &
                                 PRLZ(:,:),    &
								 
                                    DOC1ASZ(:,:), & ! Zooplankton grazing of DOC1 over the carbon threshold, or zero
												 
                                    DOC2ASZ(:,:), & ! zooplankton grazing of DOC1 over the carbon threshold, no DOC grazing of refractory DOC
								 
                                 BPRSZ(:,:),   & 
                                 BPRLZ(:,:),   &

                                    DOCSZ(:,:),  &   ! small zooplankton DOC1 mg/l
								 
                                 LPOCSZ(:,:),  &
                                 RPOCSZ(:,:),  &
								 
                                    DOCLZ(:,:),  &  ! large zooplankton DOC1 mg/l
								 
                                 LPOCLZ(:,:),  &
                                 RPOCLZ(:,:),  &
                                 NH4SZ(:,:),   &
								 
                                    DONSZ(:,:),  & ! small zooplankton DOC1 mg/l
								 
                                 LPONSZ(:,:),  &
                                 RPONSZ(:,:),  &
                                 NH4LZ(:,:),   &
								 
                                   DONLZ(:,:),  & ! large zooplankton DOC1 mg/l
								 
                                 LPONLZ(:,:),  &
                                 RPONLZ(:,:),  &
                                 PO4SZ(:,:),   &
								 
                                    DOPSZ(:,:),  & ! small zooplankton DOC1 mg/l
								 
                                 LPOPSZ(:,:),  &
                                 RPOPSZ(:,:),  &
                                 PO4LZ(:,:),   &
								 
                                    DOPLZ(:,:),  &! large zooplankton DOC1 mg/l
								 
                                 LPOPLZ(:,:),  &
                                 RPOPLZ(:,:),  &
								 
                                    !DOC2SZ(:,:),  & ! small zooplankton DOC1 mg/l
								 
                                    !DON2SZ(:,:),  & ! small zooplankton DOC1 mg/l
                                    !DOP2SZ(:,:),  & ! small zooplankton DOC1 mg/l
								 
                                    !DOC2LZ(:,:),  & ! small zooplankton DOC1 mg/l
                                    !DON2LZ(:,:),  & ! small zooplankton DOC1 mg/l
                                    !DOP2LZ(:,:),  & ! small zooplankton DOC1 mg/l
								 
                                 PIB1SZ(:,:),  &
                                 PIB2SZ(:,:),  &
                                 PIB3SZ(:,:),  &
                                 PIB1LZ(:,:),  &
                                 PIB2LZ(:,:),  &
                                 PIB3LZ(:,:),  &
                                 B1SZ(:,:),    &
                                 B2SZ(:,:),    &
                                 B3SZ(:,:),    &
                                 B1LZ(:,:),    &
                                 B2LZ(:,:),    &
                                 B3LZ(:,:),    &
                                 DOSZ(:,:),    &
                                 DOLZ(:,:),    &
                                 SASZ(:,:),    &
                                 SUSZ(:,:),    &
                                 SALZ(:,:),    &
                                 SULZ(:,:)  
	
								 
             REAL(SP),ALLOCATABLE	:: 	ACLSZ(:,:),	&
									ACLLZ(:,:),	&
									ARSZ(:,:),	&
									ARLZ(:,:),   &
									ABMSZ(:,:),	&
									ABMLZ(:,:),	&
									AMSZ(:,:),	&
									AMLZ(:,:),   &
									APRSZLZ(:,:),&
									AGSZ(:,:),	&
									AGLZ(:,:), 	&
									ADOCSZ(:,:),	&
									APOCSZ(:,:),	&
									ADOCLZ(:,:),	&
									APOCLZ(:,:), &
									ANH4SZ(:,:),	&
									ADONSZ(:,:),	&
									APONSZ(:,:),	&
									ANH4LZ(:,:),	&
									ADONLZ(:,:),	&
									APONLZ(:,:),	&
									APO4SZ(:,:),	&
									ADOPSZ(:,:),	&
									APOPSZ(:,:),	&
									APO4LZ(:,:),	&
									ADOPLZ(:,:),	&
									APOPLZ(:,:),	&
									APRSZ(:,:),		&
									APRLZ(:,:),		&
									APISZ(:,:),		&
									APILZ(:,:)

             REAL(SP),ALLOCATABLE:: 	AB1SZ(:,:), 	&
									AB2SZ(:,:),		&
									AB3SZ(:,:),		&
									AB1LZ(:,:),		&
									AB2LZ(:,:),		&
									AB3LZ(:,:),		&
									ADOSZ(:,:),		&
									ADOLZ(:,:),		&
									ASASZ(:,:),		&
									ASUSZ(:,:),		&
									ASALZ(:,:),		&
									ASULZ(:,:)

             REAL(SP),DIMENSION(-50:400) :: FTLZ,	&
										FTSZ,	&
										FTBMSZ,	&
										FTBMLZ,	&
										FTPRSZ,	&
										FTPRLZ	
										
			REAL(SP), ALLOCATABLE  :: 		SZ(:,:), &
										LZ(:,:)
			
			
		CONTAINS
		
		!Subroutines
		!	subroutine 	ZOOP_READ()
		!	subroutine	ZOOP_ALLOC()
		!	subroutine	ZOOP_DEALLOC()
		!
	
	!************************************************************************
	!**             S U B R O U T I N E   Z O OP _ R E A D                 **
	!************************************************************************

	   SUBROUTINE ZOOP_READ
	       
	       USE MOD_FILEINFO, ONLY : 		&!
					!,DIA			&!
					!,CBC 			&!
					!,S1			&!
					!,S2			&!
					!,S3			&!                 
					!,BFI			&!
			    !,BAI           &!
					!,MET			&!
					!,BFO			&!
					!,KEI			&!
					!,ATM			&!
					!,STL			&!
					!,AGR			&!
					!,SVI			&!
					!,SVO			&!
					!,KFL			&!
					ZOO				&!
					,ZFO	!		&!
					!,ALO      		&!
					!,CON			&!
					!,RSO			&!
					!,SNP			&!
					!,PLT			&!
					!,APL 			&!
					!,TFL 			&!
					!,OPL			&!
					!,SFI			&!
					!,SFO			&!
					!,MAP 			&!
					!,ICI 			&!
					!,ICO			&!
					!,MRL			&!
					!,MBL			&!
					!,RSI			&!
					!,UNIT_LINKAGE	&!
					!,UNIT_STN		&!
					!,UNIT_HIS		&!           
					!,CNAME			&!
					!,CONFN

			USE MOD_LIMS, ONLY : MLOC

		!Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
	    USE MOD_CONTROL, ONLY : 		&
				SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
				,MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
				,PAR        !	&           !!TRUE IF MULTIPROCESSOR RUN
				!,CASENAME  	&   		!!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
				!,CASETITLE  	&  			!!CASE TITLE                                 
				!,HMAX       	&  			!!GLOBAL MAXIMUM DEPTH
				!,HMIN       	&  			!!GLOBAL MINIMUM DEPTH
				!,UMOL       	&  			!!VERTICAL DIFFUSION COEFFICIENT
				!,HORCON     	&  			!!HORIZONTAL DIFFUSION COEFFICIENT
				!,DTI        	&  			!!internal time step
				!,HORZMIX    	&   		!!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
				!,FILENUMBER	&			!!
				!,PREFF			&			!!
				!,INPDIR		&			!!
				!,GEOAREA		&			!!
				!,RIV_FILENUMBER	&			!!
		    !,INFLOW_TYPE   	&			!!SPECIFIED RIVER INFLOW TYPE (edge/node) 
		    !,POINT_ST_TYPE 	&			!!(calculated/specified)
		    !,PNT_SOURCE    	&			!!point_source
		    !,DAY				&
				!,in_jday		
#if defined (MULTIPROCESSOR)
			USE MOD_PAR, ONLY:   &	  
				!EL_PID,		&!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
				!ELID,			&!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
				!NLID,			&!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
				!ELID_X,		&!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
				!NLID_X,		&!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
				!OBN_GL2LOC,	&!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
				!EGID,			&!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
				NGID,			&!(:)       	!!GLOBAL ID OF LOCAL NODE 
				!NHE,        	&!			!NUMBER OF HALO ELEMENTS
				!HE_LST,		&!(:)       	!!GLOBAL IDENTITIES OF HALO ELEMENTS
				!HE_OWN,		&!(:)       	!!OWNER OF HALO ELEMENTS
				!NBN ,       	&!          	!!NUMBER OF BOUNDARY NODES
				!MX_MLT,     	&!          	!!MAX MULTIPLICITY OF BOUNDARY NODES
				!BN_LST, 		&!	   		!!GLOBAL IDENTITY OF BOUNDARY NODES
				!BN_LOC, 		&!			!!LOCAL IDENTITY OF BOUNDARY NODES
				!BN_MLT, 		&!			!!MULTIPLICITY OF BOUNDARY NODES
				!BN_NEY, 		&!			!!NODE OWNER LIST
				!NDE_ID, 		&!			!! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
				NHN,  			&!			!!NUMBER OF HALO NODES
				HN_LST !, 		&!			!!LIST OF HALO NODES 
				!HN_OWN, 		&!			!!PRIMARY OWNER OF HALO NODES 
				!EC,			&!
				!NC,			&!
				!BNC, 			&!
				!EMAP,			&!
				!NMAP, 			&!
				!PNE,			&!(:)        !!NUMBER OF ELEMENTS IN EACH PROC
				!PNN,			&!(:)        !!NUMBER OF NODES IN EACH PROC
				!PNHE,			&!(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
				!PNBN,			&!(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
				!PMBM,			&!(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
				!PNHN,			&!(:)       !!NUMBER OF HALO NODES IN EACH PROC
				!FILELOCK,   	&
				!MYFILELOCK, 	&
				!FILELOCKED , 	&
				!NODE_MATCH, 	&
				!NODE_MATCH_ALL, &
				!EXCHANGE, 		&
				!EXCHANGE_ALL, 	&
				!GATHER, 		&
				!IGATHER, 		&
				!INIT_MPI_ENV,	&!			
				!SORT,			&
				!PREAD,			&
				!PWRITE,		&
				!PPRINT, 		&
				!LOCK, 			&
				!UNLOCK, 		&
				!GETLOC 
#endif     
	   IMPLICIT NONE
	   CHARACTER(LEN=72) :: TITLE(6)
	   
	   CHARACTER(LEN=8) :: SPVARSZ, SPVARLZ, PRINTSZ, PRINTLZ, TVARSZ,      &
			       TPRINTSZ, TVARLZ, TPRINTLZ

	! vjp 10/12/04  These variables already declared in module MOD_WQM
	!     REAL(SP)         KTGSZ1, KTGSZ2, KTBSZ, KTGLZ1, KTGLZ2, KTBLZ, 
	!    $             KTPRSZ, KTPRLZ

	   REAL(SP) :: FTB_ZOOP(0:35)
	   
	   INTEGER :: I,J,K
	   
	   !temperature control of micro-zooplankton
	   REAL(SP) :: 	TMSZ, 		&!center temperature (degC) of micro-zoop growth
				KTGSZ1,		&!left hand shape factor of temperature curve for micro-zoop growth (1/deg^2)
				KTGSZ2,		&!right hand shape factor of temperature curve for micro-zoop growth (1/deg^2)
				TRPRSZ,		&!reference temperature of predation rate
				KTPRSZ		 !shape coef of temperature curve for micro-zoop predation rate (1/degC)
				
	   !temperature control of meso-zooplankton
	   REAL(SP) :: 	TMLZ, 		&!center temperature (degC) of meso-zoop growth
				KTGLZ1, 	&!left hand shape factor of temperature curve for meso-zoop growth (1/deg^2)
				KTGLZ2, 	&!right hand shape factor of temperature curve for meso-zoop growth (1/deg^2)
				TRPRLZ,		&!reference temperature of predation rate
				KTPRLZ		 !shape coef of temperature curve for meso-zoop predation rate (1/degC)
				
	   REAL(SP) :: TLOOK	!look up able temperature for temperature control
	   
	   !temporary variables for reading data 
	   REAL(SP),ALLOCATABLE :: 	RTMP21(:,:),		&
							RTMP22(:,:),		&
							RTMP23(:,:),      	&
							RTMP24(:,:),		&
							RTMP25(:,:)

	! TITLE CARDS of zooplanton input file
							
	   OPEN(ZOO,FILE=ZOOFN,STATUS='OLD')
	   READ(ZOO,1010) (TITLE(J),J=1,6)

	! READ SPATIALLY-INVARIANT PARAMETERS FOR GROUP 1

	   READ(ZOO,1030) CTSZ, KHCSZ, DOCRITSZ, MZEROSZ
	   READ(ZOO,1030) ANCSZ, APCSZ, AOCRSZ
	   READ(ZOO,1030) UB1SZ, UB2SZ, UB3SZ, ULDSZ, URDSZ, ULPSZ, URPSZ
	   READ(ZOO,1030) TMSZ, KTGSZ1, KTGSZ2, TRSZ, KTBSZ, TRSZ,KTPRSZ
	   READ(ZOO,1030) FLPOCSZ, FRPOCSZ                      !FDOCSZ, FNCDOC2SZ,
	   READ(ZOO,1030) FNH4SZ, FUREASZ, FLPONSZ,FRPONSZ      ! FNCDON1SZ, FNCDON2SZ
	   READ(ZOO,1030) FPO4SZ, FLPOPSZ, FRPOPSZ              ! FNCDOP1SZ, FNCDOP2SZ,
	   READ(ZOO,1030) FRSASZ

	! READ SPATIALLY-INVARIANT PARAMETERS FOR GROUP 2

	   READ(ZOO,1030) CTLZ, KHCLZ, DOCRITLZ, MZEROLZ
	   READ(ZOO,1030) ANCLZ, APCLZ, AOCRLZ
	   READ(ZOO,1030) UB1LZ, UB2LZ, UB3LZ, USZLZ, ULLZ, URLZ
	   READ(ZOO,1030) TMLZ, KTGLZ1, KTGLZ2, TRLZ, KTBLZ, TRPRLZ,KTPRLZ
	   READ(ZOO,1030) FLPOCLZ, FRPOCLZ                                    !FDOCLZ, FNCDOC2LZ,
	   READ(ZOO,1030) FNH4LZ, FUREALZ, FLPONLZ,FRPONLZ                   !FDONLZ, FNCDON2LZ
	   READ(ZOO,1030) FPO4LZ,  FLPOPLZ, FRPOPLZ                       !FDOPLZ, FNCDOP2LZ,
	   READ(ZOO,1030) FRSALZ

	! ARE REMAINING GROUP 1 PARAMETERS SPATIALLY VARYING?

	   READ(ZOO,1040) SPVARSZ, PRINTSZ
	   IF (SPVARSZ == 'CONSTANT') THEN
	     READ(ZOO,1030) RMAXSZ(1,1), ESZ(1,1), RFSZ(1,1), BMRSZ(1,1), BPRSZ(1,1)
	     DO K=1,KBM1
	       DO I=1,MTLOC
		 RMAXSZ(I,K) = RMAXSZ(1,1)
		 ESZ(I,K) = ESZ(1,1)
		 RFSZ(I,K) = RFSZ(1,1)
		 BMRSZ(I,K) = BMRSZ(1,1)
		 BPRSZ(I,K) = BPRSZ(1,1)
	       ENDDO
	     ENDDO  
	   ELSE
	    ALLOCATE(RTMP21(MGL,KBM1));     RTMP21 = 0.0
	    ALLOCATE(RTMP22(MGL,KBM1));     RTMP22 = 0.0
	    ALLOCATE(RTMP23(MGL,KBM1));     RTMP23 = 0.0
	    ALLOCATE(RTMP24(MGL,KBM1));     RTMP24 = 0.0
	    ALLOCATE(RTMP25(MGL,KBM1));     RTMP25 = 0.0
	     DO K=1,KBM1
	       DO I=1,MGL
		 READ(ZOO,1030) RTMP21(I,K), RTMP22(I,K), RTMP23(I,K),    &
				RTMP24(I,K), RTMP25(I,K)
	       ENDDO
	     ENDDO      
	     IF(SERIAL)THEN
	       RMAXSZ = RTMP21
	       ESZ    = RTMP22
	       RFSZ   = RTMP23
	       BMRSZ  = RTMP24
	       BPRSZ  = RTMP25
	     ENDIF
#if defined (MULTIPROCESSOR)
	     IF(PAR)THEN
	       DO K=1,KBM1
		 DO I=1,MLOC
		   RMAXSZ(I,K) = RTMP21(NGID(I),K)
		   ESZ(I,K)    = RTMP22(NGID(I),K)
		   RFSZ(I,K)   = RTMP23(NGID(I),K)
		   BMRSZ(I,K)  = RTMP24(NGID(I),K)
		   BPRSZ(I,K)  = RTMP25(NGID(I),K)
		 ENDDO
		 DO I=1,NHN
		   RMAXSZ(I+MLOC,K) = RTMP21(HN_LST(I),K)
		   ESZ(I+MLOC,K)    = RTMP22(HN_LST(I),K)
		   RFSZ(I+MLOC,K)   = RTMP23(HN_LST(I),K)
		   BMRSZ(I+MLOC,K)  = RTMP24(HN_LST(I),K)
		   BPRSZ(I+MLOC,K)  = RTMP25(HN_LST(I),K)
		 ENDDO
	       ENDDO
	     ENDIF
#endif
	     DEALLOCATE(RTMP21,RTMP22,RTMP23,RTMP24,RTMP25)
	   ENDIF      
	      
	! ARE REMAINING GROUP 2 PARAMETERS SPATIALLY VARYING?

	   READ(ZOO,1040) SPVARLZ, PRINTLZ
	   IF (SPVARLZ == 'CONSTANT') THEN
	     READ(ZOO,1030) RMAXLZ(1,1), ELZ(1,1), RFLZ(1,1), BMRLZ(1,1), BPRLZ(1,1)
	     DO K=1,KBM1
	       DO I=1,MLOC
		 RMAXLZ(I,K) = RMAXLZ(1,1)
		 ELZ(I,K) = ELZ(1,1)
		 RFLZ(I,K) = RFLZ(1,1)
		 BMRLZ(I,K) = BMRLZ(1,1)
		 BPRLZ(I,K) = BPRLZ(1,1)
	       ENDDO
	     ENDDO  
	   ELSE
	    ALLOCATE(RTMP21(MGL,KBM1));     RTMP21 = 0.0
	    ALLOCATE(RTMP22(MGL,KBM1));     RTMP22 = 0.0
	    ALLOCATE(RTMP23(MGL,KBM1));     RTMP23 = 0.0
	    ALLOCATE(RTMP24(MGL,KBM1));     RTMP24 = 0.0
	    ALLOCATE(RTMP25(MGL,KBM1));     RTMP25 = 0.0
	     DO K=1,KBM1
	       DO I=1,MGL
		 READ(ZOO,1030) RTMP21(I,K), RTMP22(I,K), RTMP23(I,K),       &
				RTMP24(I,K), RTMP25(I,K)
	       ENDDO
	     ENDDO   
	     IF(SERIAL)THEN
	       RMAXLZ = RTMP21
	       ELZ    = RTMP22
	       RFLZ   = RTMP23
	       BMRLZ  = RTMP24
	       BPRLZ  = RTMP25
	     ENDIF
#if defined (MULTIPROCESSOR)
	     IF(PAR)THEN
	       DO K=1,KBM1
		 DO I=1,MLOC
		   RMAXLZ(I,K) = RTMP21(NGID(I),K)
		   ELZ(I,K)    = RTMP22(NGID(I),K)
		   RFLZ(I,K)   = RTMP23(NGID(I),K)
		   BMRLZ(I,K)  = RTMP24(NGID(I),K)
		   BPRLZ(I,K)  = RTMP25(NGID(I),K)
		 ENDDO
		 DO I=1,NHN
		   RMAXLZ(I+MLOC,K) = RTMP21(HN_LST(I),K)
		   ELZ(I+MLOC,K)    = RTMP22(HN_LST(I),K)
		   RFLZ(I+MLOC,K)   = RTMP23(HN_LST(I),K)
		   BMRLZ(I+MLOC,K)  = RTMP24(HN_LST(I),K)
		   BPRLZ(I+MLOC,K)  = RTMP25(HN_LST(I),K)
		 ENDDO
	       ENDDO
	     ENDIF
#endif
	     DEALLOCATE(RTMP21,RTMP22,RTMP23,RTMP24,RTMP25)
	   ENDIF

	! TIME DEPENDENCE OF PREDATION ON GROUP 1

	   READ(ZOO,1040) TVARSZ, TPRINTSZ
	   IF (TVARSZ == 'CONSTANT') THEN
	     READ(ZOO,1060) TVPRSZ(1)
	     DO J=2,366
	       TVPRSZ(J) = TVPRSZ(1)
	     ENDDO
	   ELSE
	     READ(ZOO,1060) (TVPRSZ(J), J=1,366)
	   ENDIF 

	! TIME DEPENDENCE OF PREDATION ON GROUP 2

	   READ(ZOO,1040) TVARLZ, TPRINTLZ
	   IF (TVARLZ == 'CONSTANT') THEN
	     READ(ZOO,1060) TVPRLZ(1)
	     DO J=2,366
	       TVPRLZ(J) = TVPRLZ(1)
	     ENDDO
	   ELSE
	     READ(ZOO,1060) (TVPRLZ(J), J=1,366)
	   ENDIF 

	! CREATE LOOK-UP TABLE FOR TEMPERATURE EFFECTS         

	   DO I=-50,350
	     TLOOK = FLOAT(I)/10.
	     IF (TLOOK < TMSZ) THEN
	       FTSZ(I) = EXP(-KTGSZ1*(TLOOK-TMSZ)**2)
	     ELSE
	       FTSZ(I) = EXP(-KTGSZ2*(TMSZ-TLOOK)**2)
	     ENDIF
	!       IF (T(B) < TMLZ) THEN
	!         FTLZ = EXP(-KTGLZ1*(T(B)-TMLZ)**2)
	!       ELSE
	!         FTLZ = EXP(-KTGLZ2*(TMLZ-T(B))**2)
	!       ENDIF
	     FTBMSZ(I) = EXP(KTBSZ*(TLOOK-TRSZ))
	     FTBMLZ(I) = EXP(KTBLZ*(TLOOK-TRLZ))
	     FTPRSZ(I) = EXP(KTPRSZ*(TLOOK-TRPRSZ))
	     FTPRLZ(I) = EXP(KTPRLZ*(TLOOK-TRPRLZ))
	   ENDDO

	! READ BI-MODAL TEMPERATURE FUNCTION FOR GROUP TWO

	   READ(ZOO,1050) (FTB_ZOOP(I), I=0,35)
	   DO I=-50,0
	     FTLZ(I)=0.
	   ENDDO
	   DO I=0,34
	     DO J=0,10
	       K = 10*I+J
	       FTLZ(K) = FTB_ZOOP(I)+FLOAT(J)*(FTB_ZOOP(I+1)-FTB_ZOOP(I))/10.
	     ENDDO
	   ENDDO

	   CLOSE(38)

	     
	!***** Input FORMAT statements

1010 FORMAT(A72)
1030 FORMAT(//(8X,9F8.0))
1040 FORMAT(//8X,8A8)
1050 FORMAT(//(8X,F8.0))
1060 FORMAT(//(16X,F8.0))


	! OUTPUT WHAT WAS INPUT

	! KURT GLAESEMANN - only open file on MSR
	   IF(MSR)OPEN(ZFO,FILE=ZFOFN,STATUS='UNKNOWN')
	   IF(MSR)WRITE(ZFO,2002) (TITLE(J),J=1,6)

	! WRITE SPATIALLY-INVARIANT PARAMETERS FOR GROUP 1

	   IF(MSR)WRITE(ZFO,2000)
	   IF(MSR)WRITE(ZFO,2010) CTSZ, KHCSZ, DOCRITSZ, MZEROSZ
	   IF(MSR)WRITE(ZFO,2020) ANCSZ, APCSZ, AOCRSZ
	   IF(MSR)WRITE(ZFO,2030) UB1SZ, UB2SZ, UB3SZ, ULDSZ, URDSZ, ULPSZ, URPSZ
	   IF(MSR)WRITE(ZFO,2040) TMSZ, KTGSZ1, KTGSZ2, TRSZ, KTBSZ, TRPRSZ, KTPRSZ
	   IF(MSR)WRITE(ZFO,2045)
	   IF(MSR)WRITE(ZFO,2050)  FLPOCSZ, FRPOCSZ   !FDOCSZ, FNCDOC2SZ,
	   IF(MSR)WRITE(ZFO,2060) FNH4SZ, FUREASZ,  FLPONSZ,FRPONSZ    !FNCDON1SZ, FNCDON2SZ,
	   IF(MSR)WRITE(ZFO,2070) FPO4SZ,  FLPOPSZ, FRPOPSZ    !FNCDOP1SZ, FNCDOP2SZ,
	   IF(MSR)WRITE(ZFO,2080) FRSASZ

	! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 1

	   IF (SPVARSZ == 'CONSTANT') THEN
	     IF(MSR)WRITE(ZFO,2085)
	   ELSE
	     IF(MSR)WRITE(ZFO,2086)
	   ENDIF 
	   IF(MSR)WRITE(ZFO,2090)
	   IF (PRINTSZ /= '     ALL') THEN
	     IF(MSR)WRITE(ZFO,3000) RMAXSZ(1,1), ESZ(1,1), RFSZ(1,1), BMRSZ(1,1), BPRSZ(1,1)
	   ELSE
	     DO K=1,KBM1
	       DO I=1,MLOC
		 IF(MSR)WRITE(ZFO,3010) I, K, RMAXSZ(I,K), ESZ(I,K), RFSZ(I,K), BMRSZ(I,K), &
			    BPRSZ(I,K)
	       ENDDO
	     ENDDO             
	   ENDIF
	      
	! WRITE TEMPORALLY-VARYING PARAMETERS FOR GROUP 1

	   IF (TVARSZ == 'CONSTANT') THEN
	     IF(MSR)WRITE(ZFO,2087)
	   ELSE
	     IF(MSR)WRITE(ZFO,2088)
	   ENDIF 
	   IF(MSR)WRITE(ZFO,2092)
	   IF (TPRINTSZ /= '     ALL') THEN
	     IF(MSR)WRITE(ZFO,3040) (J, TVPRSZ(J), J=1,1)
	   ELSE
	     IF(MSR)WRITE(ZFO,3040) (J, TVPRSZ(J), J=1,366)
	   ENDIF
	      
	! WRITE SPATIALLY-INVARIANT PARAMETERS FOR GROUP 2

	   IF(MSR)WRITE(ZFO,2005)
	   IF(MSR)WRITE(ZFO,2010) CTLZ, KHCLZ, DOCRITLZ, MZEROLZ
	   IF(MSR)WRITE(ZFO,2020) ANCLZ, APCLZ, AOCRLZ
	   IF(MSR)WRITE(ZFO,2035) UB1LZ, UB2LZ, UB3LZ, USZLZ, ULLZ, URLZ
	   IF(MSR)WRITE(ZFO,2040) TMLZ, KTGLZ1, KTGLZ2, TRLZ, KTBLZ, TRPRLZ, KTPRLZ
	   IF(MSR)WRITE(ZFO,2045) 
	   IF(MSR)WRITE(ZFO,2050)  FLPOCLZ, FRPOCLZ          !FDOCLZ, FNCDOC2LZ,
	   IF(MSR)WRITE(ZFO,2060) FNH4LZ, FUREALZ, FLPONLZ, FRPONLZ    ! FDONLZ, FNCDON2LZ,
	   IF(MSR)WRITE(ZFO,2070) FPO4LZ,FLPOPLZ, FRPOPLZ      ! FDOPLZ, FNCDOP2LZ, 
	   IF(MSR)WRITE(ZFO,2080) FRSALZ

	! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 2

	   IF (SPVARLZ == 'CONSTANT') THEN
	     IF(MSR)WRITE(ZFO,2085)
	   ELSE
	     IF(MSR)WRITE(ZFO,2086)
	   ENDIF 
	   IF(MSR)WRITE(ZFO,2090)
	   IF (PRINTLZ /= '     ALL') THEN
	     IF(MSR)WRITE(ZFO,3000) RMAXLZ(1,1), ELZ(1,1), RFLZ(1,1), BMRLZ(1,1), BPRLZ(1,1)
	   ELSE
	     DO K=1,KBM1
	       DO I=1,MLOC
		 IF(MSR)WRITE(ZFO,3010) I, K, RMAXLZ(I,K), ELZ(I,K), RFLZ(I,K), BMRLZ(I,K),  &
				BPRLZ(I,K)
	       ENDDO
	     ENDDO
	   ENDIF

	! WRITE TEMPORALLY-VARYING PARAMETERS FOR GROUP 2

	   IF (TVARLZ == 'CONSTANT') THEN
	     IF(MSR)WRITE(ZFO,2087)
	   ELSE
	     IF(MSR)WRITE(ZFO,2088)
	   ENDIF 
	   IF(MSR)WRITE(ZFO,2092)
	   IF (TPRINTLZ /= '     ALL') THEN
	     IF(MSR)WRITE(ZFO,3040) (J, TVPRLZ(J), J=1,1)
	   ELSE
	     IF(MSR)WRITE(ZFO,3040) (J, TVPRLZ(J), J=1,366)
	   ENDIF
	      
	! WRITE GROUP 2 TEMPERATURE FUNCTION

	   IF(MSR)WRITE(ZFO,3020)
	   DO I=0,350,10
	     IF(MSR)WRITE(ZFO,3030) I/10, FTLZ(I)
	   ENDDO

	! KURT GLAESEMANN - only close file on MSR
	   if(MSR)CLOSE(39)
		

	!***** Output FORMAT statements

2000 FORMAT(/' ZOOPLANKTON GROUP 1')
2002 FORMAT(1X,A72)
2005 FORMAT(/' ZOOPLANKTON GROUP 2')
2010 FORMAT(/' CARBON THRESHOLD = ',F8.2,' GM C/M**3'/                  &
       ' HALF-SATURATION CONCENTRATION FOR CARBON UPTAKE = ',F8.2,      &
       ' GM C/M**3'/                                                    &
       ' DO AT WHICH MORTALITY COMMENCES = ',F8.2,' GM DO/M**3'/        &
       ' MORTALITY UNDER ANOXIC CONDITIONS = ',F8.2,' /DAY')
2020 FORMAT(' NITROGEN TO CARBON RATIO = ',F8.3,' GM N/GM C'/           &
       ' PHOSPHORUS TO CARBON RATIO = ',F8.3, ' GM P/GM C'/             &
       ' RESPIRATION RATIO = ',F8.3,' GM DO/GM C')
2030 FORMAT(' UTILIZATION OF ALGAL GROUP 1   ',F8.2/                    &
	    '                ALGAL GROUP 2   ',F8.2/                    &
	    '                ALGAL GROUP 3   ',F8.2/                    &
	    '                LABILE DOC      ',F8.2/                    &
	    '                REFRACTORY DOC  ',F8.2/                    &
	    '                LABILE PARTICLES',F8.2/                    &
	    '                REFRACTORY PTCLS',F8.2)
2035 FORMAT(' UTILIZATION OF ALGAL GROUP 1   ',F8.2/                    &
	    '                ALGAL GROUP 2   ',F8.2/                    &
	    '                ALGAL GROUP 3   ',F8.2/                    &
	    '                ZOOPL GROUP 1   ',F8.2/                    &
	    '                LABILE PARTICLES',F8.2/                    &
	    '                REFRACTORY PTCLS',F8.2)
2040 FORMAT(' OPTIMUM PRODUCTION SPECIFIED AT ',F8.2,' C.'/             & 
	    ' KT1 = ',F8.3,' KT2 = ',F8.3,' PER DEGREE**2'/             &
	    ' METABOLISM SPECIFIED AT ',F8.2,' C.'/                     &
	    ' TEMPERATURE EFFECT = ',F8.3,' PER DEGREE'/                &
	    ' PREDATION SPECIFIED AT ',F8.2,' C.'/                      &
	    ' TEMPERATURE EFFECT = ',F8.3,' PER DEGREE')
2045 FORMAT(/' RECYCLING BY ZOOPLANKTON'/                               &
       '            DIS INORG  LAB DISS  REF DISS  LAB PART  REF PART')
2050 FORMAT(' CARBON    ',10X,4F10.3)
2060 FORMAT(' NITROGEN  ',6F10.3)
2070 FORMAT(' PHOSPHORUS',5F10.3)
2080 FORMAT(' SILICA    ',F10.3)
2085 FORMAT(/' REMAINING PARAMETERS ARE SPATIALLY-INVARIANT')
2086 FORMAT(/' REMAINING PARAMETERS VARY SPATIALLY')
2087 FORMAT(/' HIGHER-LEVEL PREDATION IS TEMPORALLY-INVARIANT')
2088 FORMAT(/' HIGHER-LEVEL PREDATION VARIES TEMPORALLY')
2090 FORMAT(/'    BOX   RATION   EFNCY  RESPFR   METAB  PRLOSS'/        &
	     '           1/DAY                   1/DAY M**3/GM C/'/     &
	     '                                            DAY')
2092 FORMAT(/'    DAY     TVPR')
3000 FORMAT(8X,5F8.3)
3010 FORMAT(2I8,5F8.3)
3020 FORMAT(/' TEMPERATURE FUNCTION FOR GROUP 2'/                       &
	    '    TEMP      f(T)'/)
3030 FORMAT(I8,F10.4)
3040 FORMAT(I8,F8.3)

		   RETURN
	   END SUBROUTINE ZOOP_READ



       SUBROUTINE ZOOP_ALLOC	!WLong added this for allocating zooplankton related variables

				ALLOCATE(B1ASZ(MTLOC,KBM1));        B1ASZ = 0.0 
				ALLOCATE(B2ASZ(MTLOC,KBM1));        B2ASZ = 0.0
				ALLOCATE(B3ASZ(MTLOC,KBM1));        B3ASZ = 0.0
				ALLOCATE(LPOCASZ(MTLOC,KBM1));      LPOCASZ = 0.0
				ALLOCATE(RPOCASZ(MTLOC,KBM1));      RPOCASZ = 0.0
				ALLOCATE(PRASZ(MTLOC,KBM1));        PRASZ = 0.0
				ALLOCATE(B1ALZ(MTLOC,KBM1));        B1ALZ = 0.0
				ALLOCATE(B2ALZ(MTLOC,KBM1));        B2ALZ = 0.0
				ALLOCATE(B3ALZ(MTLOC,KBM1));        B3ALZ = 0.0
				ALLOCATE(SZALZ(MTLOC,KBM1));        SZALZ = 0.0
				ALLOCATE(LPOCALZ(MTLOC,KBM1));      LPOCALZ = 0.0
				ALLOCATE(RPOCALZ(MTLOC,KBM1));      RPOCALZ = 0.0
				ALLOCATE(PRALZ(MTLOC,KBM1));        PRALZ = 0.0
				ALLOCATE(CLSZ(MTLOC,KBM1));         CLSZ = 0.0
				ALLOCATE(CLLZ(MTLOC,KBM1));         CLLZ = 0.0
				ALLOCATE(RSZ(MTLOC,KBM1));          RSZ = 0.0
				ALLOCATE(RLZ(MTLOC,KBM1));          RLZ = 0.0
				ALLOCATE(RMAXSZ(MTLOC,KBM1));       RMAXSZ = 0.0
				ALLOCATE(RMAXLZ(MTLOC,KBM1));       RMAXLZ = 0.0
				ALLOCATE(BMSZ(MTLOC,KBM1));         BMSZ = 0.0
				ALLOCATE(BMLZ(MTLOC,KBM1));         BMLZ = 0.0
				ALLOCATE(BMRSZ(MTLOC,KBM1));        BMRSZ = 0.0
				ALLOCATE(BMRLZ(MTLOC,KBM1));        BMRLZ = 0.0
				ALLOCATE(MSZ(MTLOC,KBM1));          MSZ = 0.0
				ALLOCATE(MLZ(MTLOC,KBM1));          MLZ = 0.0
				ALLOCATE(PRSZLZ(MTLOC,KBM1));       PRSZLZ = 0.0
				ALLOCATE(GSZ(MTLOC,KBM1));          GSZ = 0.0
				ALLOCATE(GLZ(MTLOC,KBM1));          GLZ = 0.0
				ALLOCATE(ESZ(MTLOC,KBM1));          ESZ = 0.0
				ALLOCATE(ELZ(MTLOC,KBM1));          ELZ = 0.0
				ALLOCATE(RFSZ(MTLOC,KBM1));         RFSZ = 0.0
				ALLOCATE(RFLZ(MTLOC,KBM1));         RFLZ = 0.0
				ALLOCATE(PRSZ(MTLOC,KBM1));         PRSZ = 0.0
				ALLOCATE(PRLZ(MTLOC,KBM1));         PRLZ = 0.0
				
				ALLOCATE(DOC1ASZ(MTLOC,KBM1));      DOC1ASZ = 0.0
				
				ALLOCATE(BPRSZ(MTLOC,KBM1));        BPRSZ = 0.0
				ALLOCATE(BPRLZ(MTLOC,KBM1));        BPRLZ = 0.0
				
				ALLOCATE(DOC2ASZ(MTLOC,KBM1));      DOC2ASZ = 0.0

				ALLOCATE(DOCSZ(MTLOC,KBM1));       DOCSZ = 0.0
				ALLOCATE(DOCLZ(MTLOC,KBM1));      DOCLZ = 0.0
				
				ALLOCATE(LPOCSZ(MTLOC,KBM1));       LPOCSZ = 0.0
				ALLOCATE(RPOCSZ(MTLOC,KBM1));       RPOCSZ = 0.0
				ALLOCATE(LPOCLZ(MTLOC,KBM1));       LPOCLZ = 0.0
				ALLOCATE(RPOCLZ(MTLOC,KBM1));       RPOCLZ = 0.0
				
				ALLOCATE(NH4SZ(MTLOC,KBM1));        NH4SZ = 0.0
				ALLOCATE(NH4LZ(MTLOC,KBM1));        NH4LZ = 0.0
				
				ALLOCATE(LPONSZ(MTLOC,KBM1));       LPONSZ = 0.0
				ALLOCATE(RPONSZ(MTLOC,KBM1));       RPONSZ = 0.0
				ALLOCATE(LPONLZ(MTLOC,KBM1));       LPONLZ = 0.0
				ALLOCATE(RPONLZ(MTLOC,KBM1));       RPONLZ = 0.0
				
				ALLOCATE(DONSZ(MTLOC,KBM1));       DONSZ = 0.0
				ALLOCATE(DONLZ(MTLOC,KBM1));       DONLZ = 0.0
				
				ALLOCATE(PO4SZ(MTLOC,KBM1));        PO4SZ = 0.0
				ALLOCATE(PO4LZ(MTLOC,KBM1));        PO4LZ = 0.0
				
				ALLOCATE(LPOPSZ(MTLOC,KBM1));       LPOPSZ = 0.0
				ALLOCATE(RPOPSZ(MTLOC,KBM1));       RPOPSZ = 0.0
				ALLOCATE(LPOPLZ(MTLOC,KBM1));       LPOPLZ = 0.0
				ALLOCATE(RPOPLZ(MTLOC,KBM1));       RPOPLZ = 0.0				
				
				ALLOCATE(DOPSZ(MTLOC,KBM1));       DOPSZ = 0.0
				ALLOCATE(DOPLZ(MTLOC,KBM1));       DOPLZ = 0.0

				
				!ALLOCATE(NCDOC2SZ(MTLOC,KBM1));       NCDOC2SZ = 0.0  ! GOING TO MAKE ALL DOC INTO ONE POOL FROM ZOOPLANKTON AND THEN SPLIT WITH FRACTIONATIONS
				!ALLOCATE(NCDON2SZ(MTLOC,KBM1));       NCDON2SZ = 0.0
				!ALLOCATE(NCDOP2SZ(MTLOC,KBM1));       NCDOP2SZ = 0.0
				!ALLOCATE(NCDOC2LZ(MTLOC,KBM1));       NCDOC2LZ = 0.0
				!ALLOCATE(NCDON2LZ(MTLOC,KBM1));       NCDON2LZ = 0.0
				!ALLOCATE(NCDOP2LZ(MTLOC,KBM1));       NCDOP2LZ = 0.0
				
				ALLOCATE(PIB1SZ(MTLOC,KBM1));       PIB1SZ = 0.0
				ALLOCATE(PIB2SZ(MTLOC,KBM1));       PIB2SZ = 0.0
				ALLOCATE(PIB3SZ(MTLOC,KBM1));       PIB3SZ = 0.0
				ALLOCATE(PIB1LZ(MTLOC,KBM1));       PIB1LZ = 0.0
				ALLOCATE(PIB2LZ(MTLOC,KBM1));       PIB2LZ = 0.0
				ALLOCATE(PIB3LZ(MTLOC,KBM1));       PIB3LZ = 0.0
		   
				ALLOCATE(B1SZ(MTLOC,KBM1));         B1SZ = 0.0    
				ALLOCATE(B2SZ(MTLOC,KBM1));         B2SZ = 0.0
				ALLOCATE(B3SZ(MTLOC,KBM1));         B3SZ = 0.0
				ALLOCATE(B1LZ(MTLOC,KBM1));         B1LZ = 0.0
				ALLOCATE(B2LZ(MTLOC,KBM1));         B2LZ = 0.0
				ALLOCATE(B3LZ(MTLOC,KBM1));         B3LZ = 0.0
				ALLOCATE(DOSZ(MTLOC,KBM1));         DOSZ = 0.0
				ALLOCATE(DOLZ(MTLOC,KBM1));         DOLZ = 0.0
				ALLOCATE(SASZ(MTLOC,KBM1));         SASZ = 0.0
				ALLOCATE(SUSZ(MTLOC,KBM1));         SUSZ = 0.0
				ALLOCATE(SALZ(MTLOC,KBM1));         SALZ = 0.0
				ALLOCATE(SULZ(MTLOC,KBM1));         SULZ = 0.0
						
				!WLong moved the following from mod_wqm.F
				ALLOCATE(ACLSZ(MTLOC,KBM1));        ACLSZ = 0.0
				ALLOCATE(ACLLZ(MTLOC,KBM1));        ACLLZ = 0.0
				ALLOCATE(ARSZ(MTLOC,KBM1));         ARSZ = 0.0
				ALLOCATE(ARLZ(MTLOC,KBM1));         ARLZ = 0.0
				ALLOCATE(ABMSZ(MTLOC,KBM1));        ABMSZ = 0.0
				ALLOCATE(ABMLZ(MTLOC,KBM1));        ABMLZ = 0.0
				ALLOCATE(AMSZ(MTLOC,KBM1));         AMSZ = 0.0
				ALLOCATE(AMLZ(MTLOC,KBM1));         AMLZ = 0.0
				ALLOCATE(APRSZLZ(MTLOC,KBM1));      APRSZLZ = 0.0
				ALLOCATE(AGSZ(MTLOC,KBM1));         AGSZ = 0.0
				ALLOCATE(AGLZ(MTLOC,KBM1));         AGLZ = 0.0
				ALLOCATE(ADOCSZ(MTLOC,KBM1));       ADOCSZ = 0.0
				ALLOCATE(APOCSZ(MTLOC,KBM1));       APOCSZ = 0.0
				ALLOCATE(ADOCLZ(MTLOC,KBM1));       ADOCLZ = 0.0
				ALLOCATE(APOCLZ(MTLOC,KBM1));       APOCLZ = 0.0
				ALLOCATE(ANH4SZ(MTLOC,KBM1));       ANH4SZ = 0.0
				ALLOCATE(ADONSZ(MTLOC,KBM1));       ADONSZ = 0.0
				ALLOCATE(APONSZ(MTLOC,KBM1));       APONSZ = 0.0
				ALLOCATE(ANH4LZ(MTLOC,KBM1));       ANH4LZ = 0.0
				ALLOCATE(ADONLZ(MTLOC,KBM1));       ADONLZ = 0.0
				ALLOCATE(APONLZ(MTLOC,KBM1));       APONLZ = 0.0
				ALLOCATE(APO4SZ(MTLOC,KBM1));       APO4SZ = 0.0
				ALLOCATE(ADOPSZ(MTLOC,KBM1));       ADOPSZ = 0.0
				ALLOCATE(APOPSZ(MTLOC,KBM1));       APOPSZ = 0.0
				ALLOCATE(APO4LZ(MTLOC,KBM1));       APO4LZ = 0.0
				ALLOCATE(ADOPLZ(MTLOC,KBM1));       ADOPLZ = 0.0
				ALLOCATE(APOPLZ(MTLOC,KBM1));       APOPLZ = 0.0
				ALLOCATE(APRSZ(MTLOC,KBM1));        APRSZ = 0.0
				ALLOCATE(APRLZ(MTLOC,KBM1));        APRLZ = 0.0
				ALLOCATE(APISZ(MTLOC,KBM1));        APISZ = 0.0
				ALLOCATE(APILZ(MTLOC,KBM1));        APILZ = 0.0
				ALLOCATE(AB1SZ(MTLOC,KBM1));        AB1SZ = 0.0  
				ALLOCATE(AB2SZ(MTLOC,KBM1));        AB2SZ = 0.0
				ALLOCATE(AB3SZ(MTLOC,KBM1));        AB3SZ = 0.0
				ALLOCATE(AB1LZ(MTLOC,KBM1));        AB1LZ = 0.0
				ALLOCATE(AB2LZ(MTLOC,KBM1));        AB2LZ = 0.0
				ALLOCATE(AB3LZ(MTLOC,KBM1));        AB3LZ = 0.0
				ALLOCATE(ADOSZ(MTLOC,KBM1));        ADOSZ = 0.0
				ALLOCATE(ADOLZ(MTLOC,KBM1));        ADOLZ = 0.0
				ALLOCATE(ASASZ(MTLOC,KBM1));        ASASZ = 0.0
				ALLOCATE(ASUSZ(MTLOC,KBM1));        ASUSZ = 0.0
				ALLOCATE(ASALZ(MTLOC,KBM1));        ASALZ = 0.0
				ALLOCATE(ASULZ(MTLOC,KBM1));        ASULZ = 0.0

				!WLong moved here from mod_wqm.F
				ALLOCATE(SZ(0:MTLOC,KBM1));    SZ   = 0.0
				ALLOCATE(LZ(0:MTLOC,KBM1));    LZ   = 0.0

           END SUBROUTINE ZOOP_ALLOC


           SUBROUTINE ZOOP_DEALLOC	!WLong added this for de-allocating zooplankton related variables

				IF(ALLOCATED(B1ASZ))		DEALLOCATE(B1ASZ)
				IF(ALLOCATED(B2ASZ))		DEALLOCATE(B2ASZ)
				IF(ALLOCATED(B3ASZ))		DEALLOCATE(B3ASZ)
				IF(ALLOCATED(LPOCASZ))	    DEALLOCATE(LPOCASZ)
				IF(ALLOCATED(RPOCASZ))		DEALLOCATE(RPOCASZ)
				IF(ALLOCATED(PRASZ))		DEALLOCATE(PRASZ)
				IF(ALLOCATED(B1ALZ))		DEALLOCATE(B1ALZ)
				IF(ALLOCATED(B2ALZ))		DEALLOCATE(B2ALZ)
				IF(ALLOCATED(B3ALZ))		DEALLOCATE(B3ALZ)
				IF(ALLOCATED(SZALZ))		DEALLOCATE(SZALZ)
				IF(ALLOCATED(LPOCALZ))		DEALLOCATE(LPOCALZ)
				IF(ALLOCATED(RPOCALZ))		DEALLOCATE(RPOCALZ)
				IF(ALLOCATED(PRALZ))		DEALLOCATE(PRALZ)
				IF(ALLOCATED(CLSZ))			DEALLOCATE(CLSZ)
				IF(ALLOCATED(CLLZ))			DEALLOCATE(CLLZ)
				IF(ALLOCATED(RSZ))			DEALLOCATE(RSZ)
				IF(ALLOCATED(RLZ))			DEALLOCATE(RLZ)
				IF(ALLOCATED(RMAXSZ))		DEALLOCATE(RMAXSZ)
				IF(ALLOCATED(RMAXLZ))		DEALLOCATE(RMAXLZ)
				IF(ALLOCATED(BMSZ))			DEALLOCATE(BMSZ)
				IF(ALLOCATED(BMLZ))			DEALLOCATE(BMLZ)
				IF(ALLOCATED(BMRSZ))		DEALLOCATE(BMRSZ)
				IF(ALLOCATED(BMRLZ))		DEALLOCATE(BMRLZ)
				IF(ALLOCATED(MSZ))			DEALLOCATE(MSZ)
				IF(ALLOCATED(MLZ))			DEALLOCATE(MLZ)
				IF(ALLOCATED(PRSZLZ))		DEALLOCATE(PRSZLZ)
				
				IF(ALLOCATED(GSZ))			DEALLOCATE(GSZ)
				IF(ALLOCATED(GLZ))			DEALLOCATE(GLZ)
				IF(ALLOCATED(ESZ))			DEALLOCATE(ESZ)
				IF(ALLOCATED(ELZ))			DEALLOCATE(ELZ)
				IF(ALLOCATED(RFSZ))			DEALLOCATE(RFSZ)
				IF(ALLOCATED(RFLZ))			DEALLOCATE(RFLZ)
				IF(ALLOCATED(PRSZ))			DEALLOCATE(PRSZ)
				IF(ALLOCATED(PRLZ))			DEALLOCATE(PRLZ)
				IF(ALLOCATED(DOC1ASZ))		DEALLOCATE(DOC1ASZ)
				IF(ALLOCATED(BPRSZ))		DEALLOCATE(BPRSZ)
				IF(ALLOCATED(BPRLZ))		DEALLOCATE(BPRLZ)
				IF(ALLOCATED(DOC2ASZ))		DEALLOCATE(DOC2ASZ)

				IF(ALLOCATED(DOCSZ))		DEALLOCATE(DOCSZ)
				IF(ALLOCATED(LPOCSZ))		DEALLOCATE(LPOCSZ)
				IF(ALLOCATED(RPOCSZ))		DEALLOCATE(RPOCSZ)
				IF(ALLOCATED(DOCLZ))		DEALLOCATE(DOCLZ)
				IF(ALLOCATED(LPOCLZ))		DEALLOCATE(LPOCLZ)
				IF(ALLOCATED(RPOCLZ))		DEALLOCATE(RPOCLZ)
				IF(ALLOCATED(NH4SZ))		DEALLOCATE(NH4SZ)
			!	IF(ALLOCATED(NCDON1SZ))		DEALLOCATE(NCDON1SZ)
				IF(ALLOCATED(LPONSZ))		DEALLOCATE(LPONSZ)
				IF(ALLOCATED(RPONSZ))		DEALLOCATE(RPONSZ)
				IF(ALLOCATED(NH4LZ))		DEALLOCATE(NH4LZ)
				IF(ALLOCATED(DONSZ))		DEALLOCATE(DONSZ)
				IF(ALLOCATED(DONLZ))		DEALLOCATE(DONLZ)
				IF(ALLOCATED(LPONLZ))		DEALLOCATE(LPONLZ)
				IF(ALLOCATED(RPONLZ))		DEALLOCATE(RPONLZ)
				IF(ALLOCATED(PO4SZ))		DEALLOCATE(PO4SZ)
				
				IF(ALLOCATED(LPOPSZ))		DEALLOCATE(LPOPSZ)
				IF(ALLOCATED(RPOPSZ))		DEALLOCATE(RPOPSZ)
				IF(ALLOCATED(PO4LZ))		DEALLOCATE(PO4LZ)
				IF(ALLOCATED(DOPSZ))		DEALLOCATE(DOPSZ)
				IF(ALLOCATED(DOPLZ))		DEALLOCATE(DOPLZ)
				IF(ALLOCATED(LPOPLZ))		DEALLOCATE(LPOPLZ)
				IF(ALLOCATED(RPOPLZ))		DEALLOCATE(RPOPLZ)
			
				
				IF(ALLOCATED(PIB1SZ))		DEALLOCATE(PIB1SZ)
				IF(ALLOCATED(PIB2SZ))		DEALLOCATE(PIB2SZ)
				IF(ALLOCATED(PIB3SZ))		DEALLOCATE(PIB3SZ)
				IF(ALLOCATED(PIB1LZ))		DEALLOCATE(PIB1LZ)
				IF(ALLOCATED(PIB2LZ))		DEALLOCATE(PIB2LZ)
				IF(ALLOCATED(PIB3LZ))		DEALLOCATE(PIB3LZ)
		   
				!WLong moved here from wqm_main.F
			
				IF(ALLOCATED(B1SZ))			DEALLOCATE(B1SZ)
				IF(ALLOCATED(B2SZ))			DEALLOCATE(B2SZ)
				IF(ALLOCATED(B3SZ))			DEALLOCATE(B3SZ)
				IF(ALLOCATED(B1LZ))			DEALLOCATE(B1LZ)
				IF(ALLOCATED(B2LZ))			DEALLOCATE(B2LZ)
				IF(ALLOCATED(B3LZ))			DEALLOCATE(B3LZ)
				IF(ALLOCATED(DOSZ))			DEALLOCATE(DOSZ)
				IF(ALLOCATED(DOLZ))			DEALLOCATE(DOLZ)
				IF(ALLOCATED(SASZ))			DEALLOCATE(SASZ)
				IF(ALLOCATED(SUSZ))			DEALLOCATE(SUSZ)
				IF(ALLOCATED(SALZ))			DEALLOCATE(SALZ)
				IF(ALLOCATED(SULZ))			DEALLOCATE(SULZ)
						
				IF(ALLOCATED(ACLSZ))		DEALLOCATE(ACLSZ)	
				IF(ALLOCATED(ACLLZ))		DEALLOCATE(ACLLZ)
				IF(ALLOCATED(ARSZ))			DEALLOCATE(ARSZ)
				IF(ALLOCATED(ARLZ))			DEALLOCATE(ARLZ)
				IF(ALLOCATED(ABMSZ))		DEALLOCATE(ABMSZ)
				IF(ALLOCATED(ABMLZ))		DEALLOCATE(ABMLZ)
				IF(ALLOCATED(AMSZ))			DEALLOCATE(AMSZ)
				IF(ALLOCATED(AMLZ))			DEALLOCATE(AMLZ)
				IF(ALLOCATED(APRSZLZ))		DEALLOCATE(APRSZLZ)
				IF(ALLOCATED(AGSZ))			DEALLOCATE(AGSZ)
				IF(ALLOCATED(AGLZ))			DEALLOCATE(AGLZ)
				IF(ALLOCATED(ADOCSZ))		DEALLOCATE(ADOCSZ)
				IF(ALLOCATED(APOCSZ))		DEALLOCATE(APOCSZ)
				IF(ALLOCATED(ADOCLZ))		DEALLOCATE(ADOCLZ)
				IF(ALLOCATED(APOCLZ))		DEALLOCATE(APOCLZ)
				IF(ALLOCATED(ANH4SZ))		DEALLOCATE(ANH4SZ)
				IF(ALLOCATED(ADONSZ))		DEALLOCATE(ADONSZ)
				IF(ALLOCATED(APONSZ))		DEALLOCATE(APONSZ)
				IF(ALLOCATED(ANH4LZ))		DEALLOCATE(ANH4LZ)
				IF(ALLOCATED(ADONLZ))		DEALLOCATE(ADONLZ)
				IF(ALLOCATED(APONLZ))		DEALLOCATE(APONLZ)
				IF(ALLOCATED(APO4SZ))		DEALLOCATE(APO4SZ)
				IF(ALLOCATED(ADOPSZ))		DEALLOCATE(ADOPSZ)
				IF(ALLOCATED(APOPSZ))		DEALLOCATE(APOPSZ)		
				IF(ALLOCATED(APO4LZ))		DEALLOCATE(APO4LZ)
				IF(ALLOCATED(ADOPLZ))		DEALLOCATE(ADOPLZ)
				IF(ALLOCATED(APOPLZ))		DEALLOCATE(APOPLZ)
				IF(ALLOCATED(APRSZ))		DEALLOCATE(APRSZ)
				IF(ALLOCATED(APRLZ))		DEALLOCATE(APRLZ)
				IF(ALLOCATED(APISZ))		DEALLOCATE(APISZ)
				IF(ALLOCATED(APILZ))		DEALLOCATE(APILZ)
				IF(ALLOCATED(AB1SZ))		DEALLOCATE(AB1SZ)
				IF(ALLOCATED(AB2SZ))		DEALLOCATE(AB2SZ)
				IF(ALLOCATED(AB3SZ))		DEALLOCATE(AB3SZ)
				IF(ALLOCATED(AB1LZ))		DEALLOCATE(AB1LZ)
				IF(ALLOCATED(AB2LZ))		DEALLOCATE(AB2LZ)
				IF(ALLOCATED(AB3LZ))		DEALLOCATE(AB3LZ)
				IF(ALLOCATED(ADOSZ))		DEALLOCATE(ADOSZ)
				IF(ALLOCATED(ADOLZ))		DEALLOCATE(ADOLZ)
				IF(ALLOCATED(ASASZ))		DEALLOCATE(ASASZ)
				IF(ALLOCATED(ASUSZ))		DEALLOCATE(ASUSZ)
				IF(ALLOCATED(ASALZ))		DEALLOCATE(ASALZ)
				IF(ALLOCATED(ASULZ))		DEALLOCATE(ASULZ)
		   
				!WLong moved this from wqm_main.F
				IF(ALLOCATED(SZ)) 			DEALLOCATE (SZ)
				IF(ALLOCATED(LZ)) 			DEALLOCATE (LZ)
		   
           END SUBROUTINE ZOOP_DEALLOC

END MODULE MOD_ZOOP
