!************************************************************************
!**                                                                    **
!** Deposition Feeder Module                                           **
!**                                                                    **
!**    Wen Long,10/05/2014                                             **
!**                                                                    **
!************************************************************************
!**                                                                    **
!** Inputs:                                                            **
!**                                                                    **
!**            Required inputs for benthic algae model                 **
!**                                                                    **
!** Outpus:                                                            **
!**                                                                    **
!**            Required inputs for benthic algae model                 **
!**                                                                    **
!************************************************************************

    MODULE MOD_DF
		
		USE MOD_WQM, ONLY:		&!
                CTEMP,          &!
                DFIFN,          &!
				DFOFN,			&!
                DLT,            &!
                DOXG,           &!
                MGL,            &!
                MTLOC!,         &!

		
        USE MOD_FILEINFO, ONLY : 		& !
				!,DIA					&!
				!,CBC 					&!
				!,S1					&!
				!,S2					&!
				!,S3					&!                 
				!BAI					&!
				!,MET					&!
				!,BAO	!				&!
				!,KEI					&!
				!,ATM					&!
				!,STL					&!
				!,AGR					&!
				!,SVI					&!
				!,SVO					&!
				!,KFL					&!
				!,ZOO					&!
				!,ZFO					&!
				!,ALO      				&!
				!,CON					&!
				!,RSO					&!
				!,SNP					&!
				!,PLT					&!
				!,APL 					&!
				!,TFL 					&!
				!,OPL					&!
				 DFI	           		&!
				,DFO!					&!
				!,SFI					&!
				!,SFO					&!
                !,MAP 					&!
				!,ICI 					&!
				!,ICO					&!
				!,MRL					&!
				!,MBL					&!
				!,RSI					&!
                !,UNIT_LINKAGE			&!
				!,UNIT_STN				&!
				!,UNIT_HIS				& !           
				!,CNAME					&!
				!,CONFN
				
USE MOD_SIZES, ONLY: MGL

USE MOD_LIMS, ONLY : MLOC, KBM1, MTLOC

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
		
	USE MOD_PREC, ONLY: SP
	
	IMPLICIT NONE
    SAVE
   
		REAL(SP) :: XKMI0,    &
				ING0,     &
				THTAI0,   &
				R,        &
				THTAR,    &
				BETA,     &
				THBETA,   &
				AMCN,     &
				AMCP,     &
				AA1,      &
				AA2,      &
				XKMG1,    &
				XKMG2,    &
				XKBO2,    &
				TDD

   		
 	REAL(SP) :: 	DFDOH,    & !oxygen concentration for hypoxia control on bioturbation
				DFDOQ,	  &	!oxygen concentration for hypoxia control on bioturbation
				DOLOW		!Never used
	
   
        !WLong moved the following from mod_wqm.F

             LOGICAL :: DFEEDER,		&!Flag to have deposition feeder calculated
						HYPOXFX_DF		!Flag for hypoxia effect on deposition feeder
						
			 LOGICAL :: DF_INITIALIZED !flag to check if deposition feeder is initialized
             
             REAL(SP) ::	XKR,					&	!respiration rate of deposition feeder
						XKI0,					&
						XKBETA

			REAL(SP)     :: RDD, 					&
						RMORT,					& !mortality rate of deposition feeder
						XPOC1LIM,				&
						TEMP,					&
						XPOC2LIM,				&
						DFSOD					 !SOD due to deposition feeder respiration
						
             REAL(SP) :: 	DFEED,					&
						DFEEDM1,                &
						ZHTAI0(350),			&
						ZHTAR(350),				&
						ZHTABETA(350)

             REAL(SP),ALLOCATABLE:: DFEEDM1S(:),		&
								DF_GROW(:),			& 
								DF_GROW_POC1(:),	& 
								DF_GROW_POC2(:),	& 
								DF_RESP(:), 		&
								DF_PRED(:), 		&
								DF_MORT(:),			&
								DF_SOD(:), 			&
								ADF_GROW(:),		&
								ADF_RESP(:),		&
								ADF_PRED(:),		&
								ADF_MORT(:)
								
		     REAL(SP), ALLOCATABLE :: DFEEDM1S_GL(:)

			LOGICAL :: DFEEDER_OUTPUT !flag for deposition feeder model output
			
			REAL(SP) ::DLT_DF		!time step (days) for deposition feeder						

   CONTAINS

   !subroutines:
		!subroutine DF_INIT()
		!subroutine DF_ALLOC()
		!subroutine DF_DEALLOC()
		!subroutine DF_READ()
		!subroutine DF_CALC()
		!subroutine DF_INT()
   
!********************************************************************************
!**                    S U B R O U T I N E   DF_INIT                           **
!********************************************************************************
  SUBROUTINE DF_INIT

    INTEGER :: I
	
	IF(.NOT. DFEEDER)THEN 
       ING0=0.0
       R=0.0
       BETA=0.0
    ENDIF

	DF_INITIALIZED = .TRUE. 
	
  END SUBROUTINE DF_INIT

!********************************************************************************
!**                    S U B R O U T I N E   DF_ALLOC                          **
!********************************************************************************
  SUBROUTINE DF_ALLOC


	       ALLOCATE(DFEEDM1S(MTLOC));           DFEEDM1S = 0.0
           ALLOCATE(DF_GROW(MTLOC));            DF_GROW = 0.0
		   ALLOCATE(DF_GROW_POC1(MTLOC));       DF_GROW_POC1 = 0.0
		   ALLOCATE(DF_GROW_POC2(MTLOC));       DF_GROW_POC2 = 0.0
           ALLOCATE(DF_RESP(MTLOC));            DF_RESP = 0.0
           ALLOCATE(DF_PRED(MTLOC));            DF_PRED = 0.0
           ALLOCATE(DF_MORT(MTLOC));            DF_MORT = 0.0
		   ALLOCATE(DF_SOD(MTLOC));             DF_SOD  = 0.0
           ALLOCATE(ADF_GROW(MTLOC));           ADF_GROW = 0.0
           ALLOCATE(ADF_RESP(MTLOC));           ADF_RESP = 0.0
           ALLOCATE(ADF_PRED(MTLOC));           ADF_PRED = 0.0
           ALLOCATE(ADF_MORT(MTLOC));           ADF_MORT = 0.0
		   ALLOCATE(DFEEDM1S_GL(MGL));          DFEEDM1S_GL = 0.0

	 
END SUBROUTINE DF_ALLOC

!********************************************************************************
!**                    S U B R O U T I N E   DF_DEALLOC                        **
!********************************************************************************

SUBROUTINE DF_DEALLOC

     !WLong moved here from wqm_main.F
	  
	
     IF(ALLOCATED(DFEEDM1S)) 		DEALLOCATE (DFEEDM1S)
     IF(ALLOCATED(DF_GROW)) 		DEALLOCATE (DF_GROW)
	 IF(ALLOCATED(DF_GROW_POC1)) 	DEALLOCATE (DF_GROW_POC1)
	 IF(ALLOCATED(DF_GROW_POC2)) 	DEALLOCATE (DF_GROW_POC2)
     IF(ALLOCATED(DF_RESP)) 		DEALLOCATE (DF_RESP)
     IF(ALLOCATED(DF_PRED)) 		DEALLOCATE (DF_PRED)
     IF(ALLOCATED(DF_MORT)) 		DEALLOCATE (DF_MORT)
	 IF(ALLOCATED(DF_SOD))  		DEALLOCATE(DF_SOD)
     IF(ALLOCATED(ADF_GROW)) 		DEALLOCATE (ADF_GROW)
     IF(ALLOCATED(ADF_RESP)) 		DEALLOCATE (ADF_RESP)
     IF(ALLOCATED(ADF_PRED)) 		DEALLOCATE (ADF_PRED)
     IF(ALLOCATED(ADF_MORT)) 		DEALLOCATE (ADF_MORT)
	 
	 IF(ALLOCATED(DFEEDM1S_GL))		DEALLOCATE(DFEEDM1S_GL)
	 
END SUBROUTINE DF_DEALLOC

!********************************************************************************
!**                    S U B R O U T I N E   D F _ R E A D                     **
!********************************************************************************
   SUBROUTINE DF_READ

   IMPLICIT NONE
   SAVE    

!***** Variable declarations

	IF(MSR)THEN
		IF (DFEEDER_OUTPUT)THEN
			OPEN (DFO, FILE=DFOFN)
		ENDIF
	ENDIF
		
	OPEN (DFI,FILE=DFIFN,STATUS='OLD')
	
		!Read the deposition feeder input file
		READ(DFI,1015,ERR=10100) XKMI0,ING0,THTAI0,R,THTAR,BETA,THBETA
		READ(DFI,1015,ERR=10100) AMCN,AMCP,AA1,AA2,XKMG1,XKMG2
		READ(DFI,1015,ERR=10100) XKBO2,TDD,DOLOW,DFDOH,DFDOQ

1015 FORMAT(//8X,8F8.1)

   CLOSE(DFI)

   DF_INITIALIZED = .FALSE.
   
   !***** Error output FORMAT's



   RETURN

!***** Error traps

10100 CONTINUE

	IF(MSR)THEN
		IF (DFEEDER_OUTPUT)THEN
			WRITE(DFO,3010)
			CLOSE(DFO)
		ENDIF
	ENDIF
3010 FORMAT(/' Read error in Deposition Feeder input deck')	
   STOP 'DF_READ'
    
   RETURN
   END SUBROUTINE DF_READ


!********************************************************************************
!**                    S U B R O U T I N E   D F _ C A L C                     **
!**                       Deposition Feeder Calculations                       **
!********************************************************************************

   SUBROUTINE DF_CALC
   USE MOD_SED_DF_EXCHANGE_VARS, ONLY: &
		POC1TM1S_SED_DF, 				&
		POC2TM1S_SED_DF,				&
		M2_SED_DF
		
   IMPLICIT NONE                     
   
   SAVE    
   INTEGER :: I,JT, ITEMP, KWC
   REAL(SP) :: TEMP20, TEMP202
   REAL(SP) :: TEMPDAY  !temperature day of year
   REAL(SP) :: LOGICT
   REAL(SP) :: O20  !oxygen concentration of overylying water (mgO2/L)
   REAL(SP) :: POC1TM1S_DF, POC2TM1S_DF	!(mgC/m^3)
   
   !Look up table for tempreature control of deposition feeder rates
   
   DO JT=1,350
     TEMP         = REAL(JT-1)/10.+0.05
     TEMP20       = TEMP-20.
     TEMP202      = TEMP20/2.
    
     ZHTAI0(JT)   = ING0*THTAI0**TEMP20           ! deposit feeders
     ZHTAR(JT)    = R*THTAR**TEMP20               ! deposit feeders
     ZHTABETA(JT) = BETA*THBETA**TEMP20           ! deposit feeders
   ENDDO

!***** Update Depisition feedr concentration

  
   KWC=KBM1
   DO I=1,MLOC
   
    O20  = AMAX1(DOXG(I,KWC),0.010)
	TEMPDAY=CTEMP(I)		!WLong: sediment temperature (degC) temporary variable
    ITEMP    = 10.*TEMPDAY+1	!calculate temperature days (using temperature to figure out the lookup
							!index in the temperature dependent rates
	
	POC1TM1S_DF=POC1TM1S_SED_DF(I)  !take POC G1 concentration of sediment
	POC2TM1S_DF=POC2TM1S_SED_DF(I)  !take POC G2 concentration of sediment
	
	IF(ITEMP<1)ITEMP=1      !make sure index >=1
	IF(ITEMP>350)ITEMP=350  !mkake sure index <=350
	
     dfeedm1  = dfeedm1s(I)  !deposition feeder befor integrating to next time step
!
! Deposit feeding ingestion rate 
!
     xki0=zhtai0(itemp)
!
! respiration rate
!
     xkr=zhtar(itemp)                           ! deposit feeders
!
! quadratic predation
!
     xkbeta=zhtabeta(itemp)                     ! deposit feeders


!MBM 970109 control switch for hypoxic effects
     if ( HYPOXFX_DF ) then

       LOGICT=1.0/(1.0 + EXP(MAX(1.1*(DFDOh-o20)/(DFDOh-DFDOq),-25.)))

!  Reduce ingestion when o2 is low
!       xki0=xki0*(o20/(o20+xkmi0))
       xki0=xki0*LOGICT

! Mortality due to hypoxia (adds to sediment POM pools)
       rdd=4.6/tdd                   ! ln(1/100) for 99% kill in time tdd
       rmort = rdd*(1.0-logict)

! Reduce predation when o2 low
       xkbeta=xkbeta*o20/(o20+xkbo2)
     endif

!
! Growth rate limitation by available POC1
!
     xpoc1lim=xkmg1/(POC1TM1S_DF+xkmg1)             ! deposit feeders
     xpoc2lim=xkmg2/(POC2TM1S_DF+xkmg2)             ! deposit feeders
!
! calculate deposit feeders biomass (in layer 2 of sediments, layer 1 ignored)
!
     
     DF_GROW_POC1(I) = aa1*xki0/(M2_SED_DF*1e+09)*POC1TM1S_DF*xpoc1lim*dfeedm1   !growth based on grazing POC1
     DF_GROW_POC2(I) = aa2*xki0/(M2_SED_DF*1e+09)*POC2TM1S_DF*xpoc2lim*dfeedm1    !growth based on grazing POC2

     DF_GROW(I) =DF_GROW_POC1(I)+DF_GROW_POC2(I)
				 
				 
     DF_RESP(I) = xkr*dfeedm1 
     DF_PRED(I) = xkbeta*dfeedm1*dfeedm1
     DF_MORT(I) = rmort*dfeedm1
	 DLT_DF = DLT/86400.0		!day
	dfeed = dfeedm1 + DLT_DF*          &
			( 	  DF_GROW(I)		 &  !growth by consuming POC1 POC2
				- DF_RESP(I)         &  !loss due to resipration
				- DF_PRED(I)         &  !loss due to predation
				- DF_MORT(I)		&  !loss due to mortality
			)


   ! dont let it go negative
     dfeed=max(dfeed,0.1)

   !calculation of SOD due to deposition feeder
	 
     DFSOD = XKR*DFEEDM1*2.667E-3       ! DEPOSIT FEEDERS RESP.
	  
	 DF_SOD(I)= DFSOD  !This will be returned to sediment module
	 
	 DFEEDM1S(I)  = DFEED
	 
	ENDDO
	   
   RETURN
   END SUBROUTINE DF_CALC
   
   SUBROUTINE DF_INT()

   INTEGER :: I
   
   RETURN
   END SUBROUTINE DF_INT


   END MODULE MOD_DF
