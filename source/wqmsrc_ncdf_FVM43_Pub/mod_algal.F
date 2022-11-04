!************************************************************************
!** This is a new version of algal input and kinetics.  Created for    **
!** The Tributary Refinements portion of the Chesapeake Bay study.     **
!**                                                                    **
!** Refinements include creation of an algal input routine, a new      **
!** P vs. I curve, a silica fraction for all algal groups, option to   **
!** spatially vary inputs, lookup table for temperature functions,     **
!** switch to uE from langleys, option to specify minimum or product   **
!** function for growth limitations.                                   **
!**                                                                    **
!**                  Algal  Subroutines for CE-QUAL-ICM                **
!**                                                                    **
!**                  Chesapeake Bay Tributary Refinements              **
!**                         February 7, 1996                           **
!**                                                                    **
!**                    Water Quality Modeling Group                    **
!**                    U.S. Army Corps of Engineers                    **
!**                    Waterways Experiment Station                    **
!**                    Vicksburg, Mississippi 39180                    **
!**                                                                    **
!************************************************************************

   MODULE MOD_ALGAL
   
	   USE MOD_PREC, ONLY: SP
	   USE MOD_SIZES, ONLY : MGL
	   USE MOD_LIMS, ONLY: KBM1, MTLOC, MLOC
	   
	   USE MOD_WQM, ONLY :		&!				
                ACHLC1,         &!
                ACHLC2,         &!
                ACHLC3,         &!
                ALOFN,          &!
                ALPHMIN1,       &!
                ALPHMIN2,       &!
                ALPHMIN3,       &!
                ALPHRAT1,       &!
                ALPHRAT2,       &!
                ALPHRAT3,       &!
                ANC1,           &!
                ANC2,           &!
                ANC3,           &!
                ASC1,           &!
                ASC2,           &!
                ASC3,           &!
                ASRAT,          &!
                B1,             &!
                B2,             &!
                ! B2GR,           &!
                B3,             &!
                BCHLC1,         &!
                BCHLC2,         &!
                BCHLC3,         &!
                BM1,            &!
                BM2,            &!
                BMR1,           &!
                BMR2,           &!
                BMR3,           &!
                BPR1,           &!
                BPR2,           &!
                BPR3,           &!
                CCHL1,          &!
                CCHL2,          &!
                CCHL3,          &!
                CCHLC1,         &!
                CCHLC2,         &!
                CCHLC3,         &!
                ! CFIX,           &!
                CHLCMN1,        &!
                CHLCMN2,        &!
                CHLCMN3,        &!
                !FCLD1,          &!now FCD1
                !FCLD2,          &!now FCD2
                !FCLD3,          &!now FCD3
               ! FCLDP,          &!!now FCDP
                FCLP1,          &!
                FCLP2,          &!
                FCLP3,          &!
                FCLPP,          &!
                !FCRD1,          &!now FCD1
                !FCRD2,          &!!now FCD2
                !FCRD3,          &!!now FCD3
                !FCRDP,          &!now FCDP
                FCRP1,          &!
                FCRP2,          &!
                FCRP3,          &!
                FCRPP,          &!
                FDOP,           &!
                FI1,            &!
                FI2,            &!
                FNFIX,          &!
                FNI1,           &!
                FNI2,           &!
                FNI3,           &!
                FNIP,           &!
                !FNLD1,          &!
                !FNLD2,          &!
                !FNLD3,          &!
             !   FNLDP,          &! ! now FNDP
                FNLP1,          &!
                FNLP2,          &!
                FNLP3,          &!
                FNLPP,          &!
                !FNRD1,          &! Now FND1
                !FNRD2,          &! Now FND2
                !FNRD3,          &!Now FND3
               ! FNRDP,          &!now FNDP
                FNRP1,          &!
                FNRP2,          &!
                FNRP3,          &!
                FNRPP,          &!
                FNUP,           &!  !Not used at all
                FPI1,           &!
                FPI2,           &!
                FPI3,           &!
                FPIP,           &!
                !FPLD1,          &! now FPD1
                !FPLD2,          &!now FPD2
                !FPLD3,          &!now FPD3
                !FPLDP,          &! now FPDP
                FPLP1,          &!
                FPLP2,          &!
                FPLP3,          &!
                FPLPP,          &!
                !FPRD1,          &!!now FPD1
                !FPRD2,          &!!now FPD2
                !FPRD3,          &!!now FPD3
               ! FPRDP,          &!! now FPDP
                FPRP1,          &!
                FPRP2,          &!
                FPRP3,          &!
                FPRPP,          &!
                FSAP,           &!
                !FT1,            &!
                !FT2,            &!
                !FT3,            &!
                !FTBM1,          &!
                !FTBM2,          &!
                !FTBM3,          &!
                !FTPR,           &!
                GPP,            &!
                IK1,            &!
                IK2,            &!
                IT,             &!
                JDAY,           &!
                KHN1,           &!
                KHN2,           &!
                KHN3,           &!
                KHNAVG,         &!
                KHNFIX,         &!
                KHP1,           &!
                KHP2,           &!
                KHP3,           &!
                KHPAVG,         &!
                !KHR1,           &!
                !KHR2,           &!
                !KHR3,           &!
                KHS1,           &!
                KHS2,           &!
                KHS3,           &!
                KHST1,          &!
                KHST2,          &!
                KHST3,          &!
                MINPROD,        &!
                NH4,            &!
                NL1,            &!
                NL2,            &!
                NO3,            &!
                NPP,            &!
                P1,             &!
                P1NNF,          &!
                P2,             &!
                PL1,            &!
                PL2,            &!
                PM1,            &!
                PM2,            &!
                PM3,            &!
                PN1,            &!
                PN2,            &!
                PO4,            &!
                PR1,            &!
                PR2,            &!
                PRSP1,          &!
                PRSP2,          &!
                PRSP3,          &!
                Q01,            &!
                Q02,            &!
                Q03,            &!
                SEDIMENT_CALC,  &!
                ! SNFIX,          &!
                STF1,           &!
                STF2,           &!
                STF3,           &!
                T,              &!
                TOTAL_NETPP,    &!
                ! TVPR,           &!
                V2,             &!
                VMAX1,          &!
                VMAX2,          &!
                VMAX3,          &!
                WS1,            &!
                WS1NET,         &!
                WS2,            &!
                WS2NET,         &!
                WS3,            &!
                WS3NET
				
       USE MOD_FILEINFO, ONLY : AGR	,ALO      	!	&!

		

		USE MOD_HYDROVARS, ONLY: DZ	,D	 
	
	!Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
		USE MOD_CONTROL, ONLY : 	&
			SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			,MSR        	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			,PAR        !	&           !!TRUE IF MULTIPROCESSOR RUN



   USE MOD_OWQ, ONLY:	IAVG

	USE MOD_ZOOP, ONLY:			&!
					 B1SZ,    	&!
					 B2SZ,    	&!

					 B1LZ,    	&!
					 B2LZ!,    	&!
					
						 
   IMPLICIT NONE
   
   !LB: moved the definitions below to this location to have acces to TRPR, TMP1, TR1, TMP2, TR2, TMP3, TR3 withing get_f*** functions
   REAL(SP) :: KTG11, KTG12, KTB1, KTG21, KTG22, KTB2, KTG31
   REAL(SP) :: KTG32, KTB3,  KTPR
   REAL(SP) :: TRPR, TMP1, TR1, TMP2, TR2, TMP3, TR3 
   
		!
		!   KURT GLAESEMANN 8 SEPT 2009 MOVE ALL VARIABLES OUT OF GLOBAL (HERE)
		!   AND DOWN TO WHERE THEY BELONG
		!
   
   CONTAINS

   
!************************************************************************
!**             S U B R O U T I N E   A L G _ R E A D                  **
!************************************************************************

   SUBROUTINE ALG_READ
#  if defined (MULTIPROCESSOR)
	USE MOD_PAR, ONLY:   &	  
			!EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			!ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			!NLID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			!ELID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			!NLID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,	&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			!EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
			!NHE,        &   !			!NUMBER OF HALO ELEMENTS
			!HE_LST,		&	!(:)       	!!GLOBAL IDENTITIES OF HALO ELEMENTS
			!HE_OWN,		&	!(:)       	!!OWNER OF HALO ELEMENTS
			!NBN ,       &   !          	!!NUMBER OF BOUNDARY NODES
			!MX_MLT,     &   !          	!!MAX MULTIPLICITY OF BOUNDARY NODES
			!BN_LST, 	&	!	   		!!GLOBAL IDENTITY OF BOUNDARY NODES
			!BN_LOC, 	&	!			!!LOCAL IDENTITY OF BOUNDARY NODES
			!BN_MLT, 	&	!			!!MULTIPLICITY OF BOUNDARY NODES
			!BN_NEY, 	&	!			!!NODE OWNER LIST
			!NDE_ID, 	&	!			!! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
			NHN,  		&	!			!!NUMBER OF HALO NODES
			HN_LST!, 	&	!			!!LIST OF HALO NODES 


#  endif
   use wc_dom, only: sanity_check,sant, FND1, FND2, FPD1, FPD2, FCD1, FCD2   
	USE MOD_WQM, only : FNDP, FCDP, FPDP
	
   IMPLICIT NONE
   real(sp) :: sanity
   INTEGER :: I, J, F, K
   REAL(SP), ALLOCATABLE :: RTMP21(:,:),RTMP22(:,:),RTMP23(:,:)
   CHARACTER(LEN=72) :: ALG_TITLE(6)
   CHARACTER(LEN=8)  :: SPVAR1, SPVAR2, SPVAR3
   CHARACTER(LEN=8)  :: PRINT1, PRINT2, PRINT3
   CHARACTER(LEN=8)  :: TPVAR,  TPRINT, TB2GR,  PRINTB2

! TITLE CARDS

     READ(AGR,1010) (ALG_TITLE(J),J=1,6)

! READ PARAMETERS WHICH ARE IDENTICAL FOR ALL GROUPS

     READ (AGR,1040)  MINPROD
     READ (AGR,*) 
     READ (AGR,*)  TRPR,   KTPR
     READ (AGR,*) 
     READ (AGR,*)  FNIP,   FNUP, FNDP,   FNLPP,   FNRPP  !need to add up to 1, FNUP not used at all, FNLDP,  FNRDP,
	 sanity = FNIP +  FNUP + FNDP + FNLPP + FNRPP
	 sant='Algal Predation Nitrogen'
	 CALL sanity_check(sant,sanity,1.0,0.011)
   READ (AGR,*) 
   
    READ (AGR,*)  FPIP,   FPDP,  FPLPP,  FPRPP !FPLDP,  FPRDP
	  sanity = FPIP +  FPDP + FPLPP + FPRPP
	  sant='Algal phoshorous'
	 CALL sanity_check(sant,sanity,1.0,0.011)
   READ (AGR,*) 
   
     READ (AGR,*)  FDOP,   FCDP, FCLPP,  FCRPP  !FCLDP,  FCRDP, moved to wc_dom into one variable, FCDP 
	 sanity =   FDOP +  FCDP + FCLPP + FCRPP
	 sant='Algal Carbon'
	 CALL sanity_check(sant,sanity,1.0,0.011)
   READ (AGR,*) 

     READ (AGR,*)  FSAP
     READ (AGR,*) 
     READ(AGR,*) ANC1, Q01, ASC1, STF1                    !
     READ (AGR,*) 
     READ(AGR,*) CHLCMN1, ACHLC1, BCHLC1, CCHLC1          !CHL to C ratios etc
     READ (AGR,*) 
     READ(AGR,*) KHN1,KHP1,KHS1,KHST1,KHNFIX,FNFIX   !KH half constants   !,KHR1 B Clark moved to wc_dom
     READ (AGR,*) 
     READ(AGR,*) ALPHMIN1,ALPHRAT1,PRSP1,VMAX1
     READ (AGR,*) 
     READ(AGR,*) TMP1,TR1                !TR is temperature reference, TMP1 is optimal temperature for grwoth
     READ (AGR,*) 
     READ(AGR,*) KTG11,KTG12,KTB1
     READ (AGR,*) 
     READ(AGR,*) FNI1,FNLP1,FNRP1 !WLong: should check if they add up to 1 after reading  
     READ (AGR,*) 
     READ(AGR,*) FPI1,FPLP1,FPRP1  !,FPLD1,FPRD1
     READ (AGR,*) 
     READ(AGR,*)FCLP1,FCRP1   ! algal carbon from respiration
     READ (AGR,*) 
     READ(AGR,*) ANC2, Q02, ASC2, STF2
     READ (AGR,*) 
     READ(AGR,*) CHLCMN2, ACHLC2, BCHLC2, CCHLC2     
     READ (AGR,*)  
     READ(AGR,*) KHN2,KHP2,KHS2,KHST2   !KHR2 B clark moved to wc_dom
     READ (AGR,*) 
     READ(AGR,*) ALPHMIN2,ALPHRAT2,PRSP2,VMAX2
     READ (AGR,*) 
     READ(AGR,*) TMP2,TR2                !TR2 : temperaure reference   !TMP2: optimal temperature for growth
     READ (AGR,*) 
     READ(AGR,*) KTG21,KTG22,KTB2
     READ (AGR,*) 
     READ(AGR,*) FNI2, FNLP2,FNRP2  !,FNLD2,FNRD2,
     READ (AGR,*) 
     READ(AGR,*) FPI2,FPLP2,FPRP2   !,FPLD2,FPRD2
     READ (AGR,*) 
     READ(AGR,*) FCLP2,FCRP2  !FCLD2,FCRD2,
     READ (AGR,*) 
     READ(AGR,*) PM1, BMR1, BPR1
     if(MSR)THEN
     WRITE(*,*) 'PM1, BMR1, BPR1'
     WRITE(*,*) PM1, BMR1, BPR1
     endif
     READ (AGR,*) 
     READ(AGR,*) PM2, BMR2, BPR2  
     if(MSR) THEN  
     WRITE(*,*) 'PM2, BMR2, BPR2'
     WRITE(*,*) PM2, BMR2, BPR2
     ENDIF

!***** Input FORMAT statements

1010 FORMAT(A72)
1030 FORMAT(//(8X,9F8.0))
1032 FORMAT(8X,9F8.0)
1040 FORMAT(//8X,8A8)
1060 FORMAT(//(16X,F8.0))

! OUTPUT WHAT WAS INPUT
     OPEN(ALO,FILE=ALOFN)
     IF(MSR)WRITE(ALO,2002) (ALG_TITLE(J),J=1,6)

! WRITE ALGAL PROPERTIES WHICH ARE NOT SPATIALLY VARYING

     IF(MSR)WRITE(ALO,3000)
     IF(MSR)WRITE(ALO,3010) CHLCMN1, CHLCMN2, CHLCMN3
     IF(MSR)WRITE(ALO,3012) ACHLC1, ACHLC2, ACHLC3
     IF(MSR)WRITE(ALO,3014) BCHLC1, BCHLC2, BCHLC3
     IF(MSR)WRITE(ALO,3016) CCHLC1, CCHLC2, CCHLC3
     IF(MSR)WRITE(ALO,3020) ANC1, ANC2, ANC3
     IF(MSR)WRITE(ALO,3022) Q01,  Q02,  Q03
     IF(MSR)WRITE(ALO,3030) ASC1, ASC2, ASC3
     IF(MSR)WRITE(ALO,3072) STF1, STF2, STF3
     IF(MSR)WRITE(ALO,3040) KHN1, KHN2, KHN3
     IF(MSR)WRITE(ALO,3050) KHP1, KHP2, KHP3
     IF(MSR)WRITE(ALO,3042) KHNFIX
     IF(MSR)WRITE(ALO,3043) FNFIX
     IF(MSR)WRITE(ALO,3060) KHS1, KHS2, KHS3
   !  IF(MSR)WRITE(ALO,3070) KHR1, KHR2, KHR3
     IF(MSR)WRITE(ALO,3074) KHST1,KHST2,KHST3
     IF(MSR)WRITE(ALO,3076) VMAX1, VMAX2, VMAX3
     IF(MSR)WRITE(ALO,3080) ALPHMIN1, ALPHMIN2, ALPHMIN3
     IF(MSR)WRITE(ALO,3090) ALPHRAT1, ALPHRAT2, ALPHRAT3
     IF(MSR)WRITE(ALO,4000) PRSP1, PRSP2, PRSP3
     IF(MSR)WRITE(ALO,4010) TMP1, TMP2, TMP3
     IF(MSR)WRITE(ALO,5020) TR1, TR2, TR3
     IF(MSR)WRITE(ALO,5030) KTG11, KTG21, KTG31
     IF(MSR)WRITE(ALO,5040) KTG12, KTG22, KTG32
     IF(MSR)WRITE(ALO,5050) KTB1, KTB2, KTB3

! WRITE DISTRIBUTION OF BIOMASS UPON MORTALITY

     IF(MSR)WRITE(ALO,4020)
     IF(MSR)WRITE(ALO,4030) FNI1,FNLP1,FNRP1 !,FNLD1,FNRD1   ! FNI1,FNLP1 etc + FND1 must add up to 1!!!!!
     IF(MSR)WRITE(ALO,4040) FPI1,FPLP1,FPRP1 !,FPLD1,FPRD1
     IF(MSR)WRITE(ALO,4050) FCLP1,FCRP1  ! FCLD1,FCRD1,
     IF(MSR)WRITE(ALO,4060)
     IF(MSR)WRITE(ALO,4030) FNI2,FNLP2,FNRP2 !,FNLD2,FNRD2
     IF(MSR)WRITE(ALO,4040) FPI2,FPLP2,FPRP2 !,FPLD2,FPRD2
     IF(MSR)WRITE(ALO,4050) FCLP2,FCRP2 !FCLD2,FCRD2
     IF(MSR)WRITE(ALO,4070)
     IF(MSR)WRITE(ALO,4030) FNI3,FNLP3,FNRP3 !,FNLD3,FNRD3
     IF(MSR)WRITE(ALO,4040) FPI3,FPLP3,FPRP3 !,FPLD3,FPRD3
     IF(MSR)WRITE(ALO,4050) FCLP3,FCRP3 !FCLD3,FCRD3
     IF(MSR)WRITE(ALO,4090)
     IF(MSR)WRITE(ALO,4030) FNIP,FNDP,FNLPP,FNRPP
     IF(MSR)WRITE(ALO,4040) FPIP,FPDP,FPLPP,FPRPP
     IF(MSR)WRITE(ALO,4080) FDOP,FCDP,FCLPP,FCRPP
     IF(MSR)WRITE(ALO,4082) FSAP,(1.-FSAP)
 
     IF (MINPROD == ' MINIMUM' .OR. MINPROD == ' PRODUCT') THEN
       IF(MSR)WRITE(ALO,2010) MINPROD
     ELSE
       IF(MSR)WRITE(ALO,2020) MINPROD
       STOP
     ENDIF

! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 1

     IF(MSR)WRITE(ALO,2000)
     IF (SPVAR1 == 'CONSTANT') THEN
       IF(MSR)WRITE(ALO,2085)
     ELSE
       IF(MSR)WRITE(ALO,2086)
     ENDIF
     IF(MSR)WRITE(ALO,5060)
     IF (PRINT1 /= '     ALL') THEN
       IF(MSR)WRITE(ALO,5070) PM1, BMR1, BPR1
     ENDIF
      
! WRITE SPATIALLY-VARYING PARAMETERS FOR GROUP 2

     IF(MSR)WRITE(ALO,2005)
     IF (SPVAR2 == 'CONSTANT') THEN
       IF(MSR)WRITE(ALO,2085)
     ELSE
       IF(MSR)WRITE(ALO,2086)
     ENDIF
     IF(MSR)WRITE(ALO,5060)
     IF (PRINT2 /= '     ALL') THEN
       IF(MSR)WRITE(ALO,5070) PM2, BMR2, BPR2
     ENDIF


      
!***** Output FORMAT statements

2000 FORMAT(/' ALGAL GROUP 1')
2002 FORMAT(1X,A72)
2005 FORMAT(/' ALGAL GROUP 2')
2006 FORMAT(/' ALGAL GROUP 3')
2010 FORMAT(/A8,' FORMULATION SELECTED FOR LIGHT AND NUTRIENT',' LIMITATIONS')
2020 FORMAT(/' MINPROD INCORRECTLY SPECIFIED AS ',A8)
2040 FORMAT(' PREDATION SPECIFIED AT ',F8.2,' C.'/              &
            ' TEMPERATURE EFFECT = ',F8.3,' PER DEGREE')
2085 FORMAT(/' REMAINING PARAMETERS ARE SPATIALLY-INVARIANT')
2086 FORMAT(/' REMAINING PARAMETERS ARE SPATIALLY-VARYING')
2087 FORMAT(/' HIGHER-LEVEL PREDATION IS TEMPORALLY-INVARIANT')
2088 FORMAT(/' HIGHER-LEVEL PREDATION VARIES TEMPORALLY')
2090 FORMAT(I8,F8.3)
2092 FORMAT(/'    DAY     TVPR')
2093 FORMAT(/'    DAY     B2GR')
2094 FORMAT(/' BASE GROUP 2 GROWTH IS TEMPORALLY-INVARIANT')
2095 FORMAT(/' BASE GROUP 2 VARIES TEMPORALLY')
3000 FORMAT(/' ALGAL PROPERTIES',30X,'GROUP 1   GROUP 2   GROUP 3')
3010 FORMAT(/' MINIMUM GM CHL/GM C    ',20X,3F10.4)
3012 FORMAT(' ACHCL (GM CHL/GM C)    ',20X,3F10.4)
3014 FORMAT(' BCHCL (1/DEGREE C)     ',20X,3F10.4)
3016 FORMAT(' CCHCL (M**2 DAY/E)     ',20X,3F10.4)
3020 FORMAT(' GM N/GM C              ',20X,3F10.4)
3022 FORMAT(' MINIMUM GM P/GM C      ',20X,3F10.4)
3030 FORMAT(' GM SI/GM C             ',20X,3F10.4)
3040 FORMAT(' KHN (GM N/M**3)        ',20X,3F10.4)
3042 FORMAT(' KHNFIX (GM N/M**3)     ',20X,F10.4)
3043 FORMAT(' FRACTION N FIXERS      ',20X,F10.4)
3050 FORMAT(' KHP (GM P/M**3)        ',20X,3F10.4)
3060 FORMAT(' KHS (GM SI/M**3)       ',20X,3F10.4)
3070 FORMAT(' KHR (GM DO/M**3)       ',20X,3F10.4)
3072 FORMAT(' SALT TOXICITY (1/DAY)  ',20X,3F10.4)
3074 FORMAT(' KHST (ppt)             ',20X,3F10.4)
3076 FORMAT(' VMAX (GM P/GM C/DAY)   ',20X,3F10.4)
3080 FORMAT(' ALPHA (GM C/GM CHL/DAY)/(uE/M**2/SEC)      ',3F10.3)
3090 FORMAT(' ALPHA TO P RATIO (1/(uE/M**2/SEC))',9X,3F10.3)
4000 FORMAT(' PHOTORESPIRATION FRACTION     ',13X,3F10.3)
4010 FORMAT(' OPTIMAL TEMPERATURE FOR PRODUCTION (C)',5X,3F10.3) 
4020 FORMAT(/' DISTRIBUTION OF ALGAE UPON MORTALITY'//                       &
     ' GROUP 1 RESPIRATION  DIS INORG  LAB DISS  REF DISS  LAB PART  REF PART')
4030 FORMAT(' NITROGEN            ',5F10.3)
4040 FORMAT(' PHOSPHORUS          ',5F10.3)
4050 FORMAT(' CARBON              ',10X,4F10.3)
4060 FORMAT(/' GROUP 2 RESPIRATION')
4070 FORMAT(/' GROUP 3 RESPIRATION')
4080 FORMAT(' CARBON              ',5F10.3)
4082 FORMAT(' SILICA              ',F10.3,10X,F10.3)
4090 FORMAT(/' PREDATION')
5000 FORMAT(' PHOSPHORUS     ',4F10.3)
5010 FORMAT(' CARBON         ',4F10.3)
5020 FORMAT(' REFERENCE TEMPERATURE FOR RESPIRATION (C)  ',3F10.3)
5030 FORMAT(' EFFECT OF SUBOPTIMAL TEMP ON PRODUCTION    ',3F10.4) 
5040 FORMAT(' EFFECT OF SUPEROPTIMAL TEMP ON PRODUCTION  ',3F10.4) 
5050 FORMAT(' EFFECT OF TEMPERATURE ON RESPIRATION       ',3F10.4)
5060 FORMAT('       CELL     PMAX     METAB     PRDTN',/                   &
            '             C/CHL/DAY   1/DAY     1/DAY'/)
5070 FORMAT(10X,3F10.3)
5080 FORMAT(2I10,3F10.3)

     RETURN
     END SUBROUTINE ALG_READ
      


!************************************************************************
!**                  S U B R O U T I N E   A L G A E                   **
!************************************************************************

   SUBROUTINE ALGAE(DTB1,DTB2,DTB3)
   IMPLICIT NONE
   INTEGER  ::  I, J, F, K, ITEMP
   REAL(SP),DIMENSION(0:MTLOC,KBM1) :: DTB1, DTB2,DTB3

   REAL(SP) :: XL=0.
   REAL(SP) :: FI01, FI02, FI03, GPP1, GPP2, GPP3, TREC
   REAL(SP) :: NETP1, NETP2, NETP3, IK, ITTOP, ITAVG, NUTLIM
   REAL(SP) :: SALTOX, DIN, DSIL, ALPHA, OPTDEPTH



! DETERMINE JULIAN DAY

   J = 1.0 + AMOD(JDAY,365.25)
!     if (jday .gt. 730.) fnfix = 0.9
   total_netPP(:)=0.0 
  DO K=1,KBM1
     DO I=1,MLOC

       P1(I,K) = PM1*GET_FT(T(I,K),TMP1,KTG11,KTG21)            !!LB: replace FT1(ITEMP)-lookup-table with "get_FT1(T)" -NOW get_FT(T,TMP,KTG1,KTG2)
       P2(I,K) = PM2*GET_FT(T(I,K),TMP2,KTG21,KTG22)!*B2GR(J)   ! B2Gr = 1

       BM1(I,K) = BMR1*GET_FTBM(T(I,K),TR1,KTB1)   !ALG1 metabolism 
       BM2(I,K) = BMR2*GET_FTBM(T(I,K),TR2,KTB2)  !ALG2 metabolism
      
	      PR1(I,K) = BPR1*B1(I,K)*B1(I,K)*GET_FTPR(T(I,K),TRPR,KTPR)!  !gC/m^3/day ??  !Laura: should not repeat B1	   
        PR2(I,K) = BPR2*B2(I,K)*B2(I,K)*GET_FTPR(T(I,K),TRPR,KTPR)!*TVPR(J)	 !Laura: should not repeat B2

       DIN = NH4(I,K)+NO3(I,K)

!Wen Long : fix bug in nitrogen limiation (one should consider ammonia preference PN)
!           when PN -->1,       (i.e. KHN becomes small), NL should be NH4/(KHN+NH4) 
!           when PN --> NH4/DIN,(i.e. KHN becomes large), NL should be DIN/(KHN+DIN)

           !ammonia preference  (moved from wqm_kin.F to here)  
           PN1(I,K) = NH4(I,K)*(NO3(I,K)/((KHN1+NH4(I,K))*(KHN1+NO3(I,K)))+KHN1  &
                /((1.E-30+NH4(I,K)+NO3(I,K))*(KHN1+NO3(I,K))))
           PN2(I,K) = NH4(I,K)*(NO3(I,K)/((KHN2+NH4(I,K))*(KHN2+NO3(I,K)))+KHN2  &
                /((1.E-30+NH4(I,K)+NO3(I,K))*(KHN2+NO3(I,K))))

           NL1(I,K) = ( NH4(I,K) + DIN - PN1(I,K)*DIN)   &
                     /(KHN1+NH4(I,K) + DIN - PN1(I,K)*DIN+1E-30)

           NL2(I,K) = (NH4(I,K) + DIN - PN2(I,K)*DIN)   &
                     /(KHN2+NH4(I,K) + DIN - PN2(I,K)*DIN+1E-30)

!Michaelis-Menten here...
#if defined(PHOSPHORUS)                     
        PL1(I,K)= PO4(I,K)/(KHP1+PO4(I,K))
        PL2(I,K)= PO4(I,K)/(KHP2+PO4(I,K))
!        PL3(I,K)= PO4(I,K)/(KHP3+PO4(I,K))
#else
    PL1(I,K)=1.0
    PL2(I,K)=1.0
#endif
!! SILICA 

!   DO K=1,KBM1
!     DO I=1,MLOC
!       DSIL = SIAT(I,K)/(1.+KADSA*SSI(I,K))             !Wen Long: SIAT is total concentration of silca
!       SL1(I,K) = (DSIL+1.0E-30)/(KHS1+DSIL+1.0E-30)
!       SL2(I,K) = (DSIL+1.0E-30)/(KHS2+DSIL+1.0E-30)
!       SL3(I,K) = (DSIL+1.0E-30)/(KHS3+DSIL+1.0E-30)
!     ENDDO
!   ENDDO  

! PHOTOSYNTHESIS VS IRRADIANCE

       ALPHA = ALPHMIN1+ALPHRAT1*P1(I,K)

!RGL changed IK to reflect Pmax/alpha - not modified P
       IK= P1(I,K)/ALPHA        !Wen Long, I believe this is correct
                                          !I agree with RGL that IK calculation should not have 
                                          !nitrogen limitation in it
       IK1(I,K)=IK

       FI1= IAVG(I,K)/SQRT(IK*IK+IAVG(I,K)*IAVG(I,K))

       ALPHA = ALPHMIN2+ALPHRAT2*P2(I,K)

       IK= P2(I,K)/ALPHA      !Wen Long, I believe this is correct 
	   
       IK2(I,K)=IK

       FI2 = IAVG(I,K)/SQRT(IK*IK+IAVG(I,K)*IAVG(I,K))

!     *** Compute carbon to chlorophyll ratio in each cell 
! 14 SEPT 2009 KURT GLAESEMANN - THESE ARE SOMETIME ALL ZERO - MUCH SIMPLER CODE
   IF (ACHLC1 .EQ. 0 .AND. ACHLC2 .EQ. 0 .AND. ACHLC3 .EQ. 0) THEN

       CCHL1(I,K) = 1./CHLCMN1
       CCHL2(I,K) = 1./CHLCMN2
      !  CCHL3(I,K) = 1./CHLCMN3

   ELSE ! 14 SEPT 2009
 
       NUTLIM = AMIN1(NL1(I,K),PL1(I,K))

       CCHL1(I,K) = 1./(CHLCMN1+ACHLC1*EXP(BCHLC1*T(I,K))*NUTLIM/EXP(CCHLC1*IAVG(I,K)))
       NUTLIM = AMIN1(NL2(I,K),PL2(I,K))

       CCHL2(I,K) = 1./(CHLCMN2+ACHLC2*EXP(BCHLC2*T(I,K))*NUTLIM/EXP(CCHLC2*IAVG(I,K)))

   ENDIF

   IF (MINPROD == ' MINIMUM') THEN

!RGL below without P or Si limitation
         P1(I,K)    = P1(I,K)*AMIN1(FI1,NL1(I,K),PL1(I,K))/CCHL1(I,K)
         P2(I,K)    = P2(I,K)*AMIN1(FI2,NL2(I,K),PL2(I,K))/CCHL2(I,K)

   ELSE
!RGL below without P or Si limitation
         P1(I,K)    = P1(I,K)*FI1*AMIN1(   NL1(I,K),PL1(I,K))/CCHL1(I,K)
         P2(I,K)    = P2(I,K)*FI2*AMIN1(   NL2(I,K),PL2(I,K))/CCHL2(I,K)

   ENDIF
! RATE OF CHANGE DUE TO PRODUCTION, RESPIRATION, PREDATION

       NETP1 = (P1(I,K)*(1.-PRSP1)-BM1(I,K))*B1(I,K)  ! predation - metabolism
       NETP2 = (P2(I,K)*(1.-PRSP2)-BM2(I,K))*B2(I,K)

       GPP1 = P1(I,K)*B1(I,K)
       GPP2 = P2(I,K)*B2(I,K)

! PRIMARY PRODUCTION IN GM C/M**2/DAY
        GPP(I,K) = (GPP1+GPP2)*D(I)*DZ(K)
        NPP(I,K) = (NETP1+NETP2)*D(I)*DZ(K)

        DTB1(I,K) = (NETP1-PR1(I,K))/86400.

        DTB2(I,K) = (NETP2-PR2(I,K))/86400.
        total_netPP(I)=total_netPP(I)+NPP(I,K)
 
#if defined(ZOOPLANKTON)      
 
  DTB1(I,K)=DTB1(I,K)-(B1SZ(I,K)-B1LZ(I,K))/86400.
  DTB2(I,K)=DTB2(I,K)-(B2SZ(I,K)-B2LZ(I,K))/86400.
#endif
      
     ENDDO
   ENDDO 

! RATE OF CHANGE DUE TO SETTLING
   DO I=1,MLOC
! 14 SEPT 2009 - KURT GLAESEMANN CHANGED X/Y/86400 TO X/(Y*86400)
     DTB1(I,1) = DTB1(I,1)-WS1*B1(I,1)/(D(I)*DZ(1)*86400.)
     DTB2(I,1) = DTB2(I,1)-WS2*B2(I,1)/(D(I)*DZ(1)*86400.)
     DTB3(I,1) = DTB3(I,1)-WS3*B3(I,1)/(D(I)*DZ(1)*86400.)
   ENDDO          

   DO K=2,KBM1
     DO I=1,MLOC
! 14 SEPT 2009 - KURT GLAESEMANN CHANGED X/Y/86400 TO X/(Y*86400)
       DTB1(I,K) = DTB1(I,K)+(WS1*B1(I,K-1)-WS1*B1(I,K))  &
                  /(D(I)*DZ(K)*86400.)
       DTB2(I,K) = DTB2(I,K)+(WS2*B2(I,K-1)-WS2*B2(I,K))  &
                  /(D(I)*DZ(K)*86400.)                                    
       DTB3(I,K) = DTB3(I,K)+(WS3*B3(I,K-1)-WS3*B3(I,K))  &
                  /(D(I)*DZ(K)*86400.)
     ENDDO
   ENDDO  


!  RESUSPENSION

   IF (SEDIMENT_CALC) THEN 
! 14 SEPT 2009 - KURT GLAESEMANN MAKE THIS SIMPLER LOOP
   
   !WLong and Laura Bianucci comments:
   !     resuspension rate should be (WS1-WS1NET ) in m/d 
   !  and resuspension is counter-acting on the sinking of WS1 
   ! so the net sinking is WS1NET = WS1 - (WS1-WS1NET) , i.e. sinking - resuspension 
   ! 
   
     	DO I=1,MLOC
		  XL=1.0/(D(I)*DZ(KBM1)*86400)
		  DTB1(I,KBM1) = DTB1(I,KBM1)+(WS1-WS1NET)*B1(I,KBM1)*XL
		  DTB2(I,KBM1) = DTB2(I,KBM1)+(WS2-WS2NET)*B2(I,KBM1)*XL
		  DTB3(I,KBM1) = DTB3(I,KBM1)+(WS3-WS3NET)*B3(I,KBM1)*XL
		ENDDO

   ENDIF 


   RETURN
   END SUBROUTINE ALGAE
   
   
   
   FUNCTION GET_FT(TEMPVAL,TMP,KTG1,KTG2)     !m/d  
   REAL(SP) ::     TEMPVAL,TMP,KTG1,KTG2             !temperature (degC)
   REAL(SP) ::    GET_FT       
       IF (TEMPVAL < TMP) THEN                  !shape is different for < TMP1 and > TMP1
         GET_FT = EXP(-KTG1*(TEMPVAL-TMP)**2)  !TMP1 is optimal temperature 
       ELSE                                    !KTG11 is the shape parameter (bell shape) for T < TMP1
         GET_FT = EXP(-KTG2*(TMP-TEMPVAL)**2)  !KTG12 is for T > TMP1
       ENDIF
   RETURN
   END FUNCTION GET_FT
   
   FUNCTION GET_FTBM(TEMPVAL,TR,KTB)     !m/d  
   REAL(SP) ::     TEMPVAL,TR,KTB             !temperature (degC)
   REAL(SP) ::    GET_FTBM      
       GET_FTBM=EXP(KTB*(TEMPVAL-TR))
   RETURN
   END FUNCTION GET_FTBM
   
   FUNCTION GET_FTPR(TEMPVAL,TRPR,KTPR)     !m/d  
   REAL(SP) ::     TEMPVAL,TRPR,KTPR             !temperature (degC)
   REAL(SP) ::    GET_FTPR    
       GET_FTPR=EXP(KTPR*(TEMPVAL-TRPR))
   RETURN
   END FUNCTION GET_FTPR
   
   END MODULE MOD_ALGAL
