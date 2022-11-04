
MODULE MOD_OWQ
      
!********************************************************************************
!**                             OWQ Model Setup (Ocean Water Color)            **
!********************************************************************************
 USE MOD_WQM, ONLY :			&		
                B1,             &!
                B2,             &!
                B3,             &!
                CHLCMN1,        &!
                CHLCMN2,        &!
                CHLCMN3,        &!
                I0,             &!
                JDAY,           &!
               ! LDOC,           &!
                LPOC,           &!
              !  RDOC,           &!
                RPOC,           &!
                SREFLECT,       &!
                SSI,   &!,            
				SPVARKE,PRINTKE

   !USE MOD_LIMS, ONLY: MLOC
   USE MOD_HYDROVARS, ONLY:     &
		  DZ	&		!DELTA-SIGMA VALUE
		 ,D	!&			!CURRENT DEPTH 

        USE MOD_LIMS, ONLY: MTLOC, 		&
						    KBM1,		&
							MLOC
							
		USE MOD_SIZES, ONLY: MGL
		USE MOD_FILEINFO, ONLY: KEI
		
		USE MOD_CONTROL, ONLY: PAR,	&
								SERIAL		
# if defined (MULTIPROCESSOR)

		USE MOD_PAR, ONLY:  NHN,		&
							NGID,		&
							HN_LST
# endif
							
	    USE MOD_PREC, ONLY: SP										
		
		USE MOD_ZOOP, ONLY : 	SZ,		&
								LZ
		
								
#if defined(PHOTODEGRADATION)				 
        USE WC_DOM, only:   WC_CDOC1,WC_CDOC2,WC_CDOC3,  & !get our colored DOC from mod_wc_dom B Clark Sep 2015		
						   ICDOM_abs, PHOTODEG, NWAVELx,   &
						   N_photons, simpson
																		   
#else						   
        USE WC_DOM, only:   WC_CDOC1,WC_CDOC2,WC_CDOC3,NWAVELx!,  & !get our colored DOC from mod_wc_dom B Clark Sep 2015	
		
#endif		
						   
IMPLICIT NONE		


		REAL(SP) :: 	INITKE, 			&
					INTKE,				&
					NXKEI
							
		REAL(SP) :: 	LAT, 				&	
					LON
 

		
			REAL(SP) ::	PHISTAR676,			&
					    BNSTARNTU,			&
					    BNSTARCHL,      	&
					   PSTARINT,			&
					   PSTAR440,			&
					   PSTARCOLOR,			&
					   SCDOM,				&
					   EZEROINT,    		&
					   TCHL2CHL,			&
					   DOFFSET,			&
					   DEPTHKE		



        INTEGER :: NWAVEL, NWAVEL_PAR!number of wavelengths, and number of wavelengths in photosynthesis region

        REAL(SP),ALLOCATABLE,DIMENSION(:,:) :: TOTCHL
												
		REAL(SP), ALLOCATABLE:: TSS(:,:)   ! B CLark June 2017 make TSS allocatable so that we can store it before it gets passed into photochemistry
			REAL(SP):: TOTPOC										
#if defined (PHOTODEGRADATION)

		REAL(SP) ::  ICDOM_bot      ! POC   					
		
		REAL(SP), ALLOCATABLE :: IAVG_pd(:,:,:)   ! avg light in photoedgradation (MTLOC,KBM1,NWAVEL)
		REAL(SP), ALLOCATABLE :: IAVG_in(:)   ! light multiplied by wavelength and converted to photons before integration

	    REAL(SP),ALLOCATABLE ::  ABSCDOM(:,:,:,:) !,   & ! absorbance due to cdom at each wavelevnth (MTLOC,KBM1,NWAVEL,3)
		
		REAL(SP) :: WtoP = 4.3459E20/6.0221409E+23! convert watts to photons   		  

	   
	   INTEGER ::  PAR_MARK  ! where visible light begins in our wave spectra
	   
	   REAL(SP), ALLOCATABLE :: I0_lambda(:)     ! spectrally dependnet surface irradiance (Calcualte from SWRAD and broken down into spectra)
	   REAL(SP), ALLOCATABLE :: KDLAMB_out(:,:,:) ,& ! our spectral dependency of the total Kd (m^-1) (MTLOC,KBM1,NWAVEL)
								KDLAMB_int(:,:)  ! integrated Kd over PAR spectrum (400-700nm)
	   REAL(SP), ALLOCATABLE :: KD_vec(:)
	   REAL(SP),ALLOCATABLE  :: GLAMB(:,:)   ! The absorbance curves for our CDOM 1,2,3 at each wavelenght. read from the input file kei  ! should make these allocatable
	   REAL(SP):: GSTAR440		! absorbance of our CDOM 1,2,3 at 440 nm (reference wavelength)
#else
						
	  REAL(SP), DIMENSION(150) ::	GLAMB  ! not going to have 3 in originial model
	  REAL(SP) :: GSTAR440
					
#endif	
	  
												
        REAL(SP),SAVE :: 	OPZ(9), 	&
						G1(9),		&
						G2(9)
						
	     REAL(SP),ALLOCATABLE, DIMENSION(:) ::  ALAMB, 	&! !water		  ! B Clark: Before these vars were just assigned an arbitiraty dimension to fill, now they will be allocated to the size of the wavelength 
																	! array  
								
								PLAMB,			& ! particles
								PHILAMB,        &! chla 
                                BLAMB,			& !scattering of particles
								EZERO, EZ	  ! irradiance distribution (total = 1
								
		INTEGER, ALLOCATABLE :: 	WAVEL(:)			 !wavelength	
		INTEGER :: LAMBDA_INT !interval of the wavelengths

						
									
        REAL(SP),DIMENSION(10) :: GMUZERO
		
		REAL(SP), ALLOCATABLE::	IAVG(:,:),	& ! Average irradiance for each cell at layer k (MTLOC,KBM1)
											
							
								IATBOT_out(:)  !Irradiance reaching the sediments (MTLOC)
												
		 INTEGER:: I,J,K,IDOM,Iscat
! from subroutine LIGHT ATN		 
   REAL(SP)    :: ITHTAA

   REAL(SP),PARAMETER :: PI=3.141592653589793, &
                      TWOPI = 6.283185307179586
!vjp 9/19/2005 added declarations
   REAL(SP) :: DOY, HOUR, SOLDATE, DECL, TAU, CTHTAA, THTAW, AVREFCOS
   REAL(SP) ::TOTDOC, IATTOP, OPTDEPTH 


#if defined(PHOTODEGRADATION)
    REAL(SP), ALLOCATABLE :: IATBOT(:)   ! need to make left over irradiance!allocatble if using the photodegradation
	
#else
    REAL(SP) :: IATBOT
#endif
   
   REAL(SP) :: CDOM_optdepth ! Added by B Clark for photodegradation  Nov 2015
    
		 
! from subroutine photochemistry  
   REAL(SP) :: PLR,DEPTH,REFSCAT,EZINT

   REAL(SP) :: ABCHL, ABNAPM, ABTOT, SCAT, KDLAMB	 
	
      REAL(SP) SUMI

CONTAINS

!=========================================================================
!
!=========================================================================
   
   SUBROUTINE OWQ_READ
       !USE MOD_LIMS, ONLY: MLOC

	   REAL(SP), ALLOCATABLE :: RTMP2(:,:),RTMP21(:,:)
	   
       INTEGER :: I, K, CDOM_wavel, iCDOM
	   	   
#if defined(PHOTODEGRADATION)       
    !    write(*,*)'Reading the Light Field definition file'
       READ(KEI,1032)
       READ(KEI,1038) INTKE, INITKE
	   


	 ! 
	   	   
       READ(KEI,1037) GSTAR440,PHISTAR676,BNSTARNTU,BNSTARCHL,PSTARINT,   &   ! dont need a spectral GSTAR440 if we have the CDOM normalized absorption curves
                      PSTAR440,PSTARCOLOR,SCDOM  ! just read in the same as before
			  
					  


       READ(KEI,1036) LAT, LON, NWAVEL, TCHL2CHL, DOFFSET, DEPTHKE  ! can just modify the input file so it reads in more wavelengths, pretty easy!!!

	CALL OWQ_ALLOC   
       NXKEI = INTKE/24.


! CONVERT FROM DECIMAL DEGREES TO RADIANS
       LAT = 6.2832*LAT/360.
       LON = 6.2832*LON/360.

       READ(KEI,1032)
	   


	!   ALLOCATE(CDOM1_ABSin(CDOM_wavel,CDOM_counter))  ! allocate the variable to the number of wavelengths x the number of spectra
	IF(PHOTODEG)THEN  
	  DO I =1 , NWAVEL   
	    !write(*,*)'Reading wavelenght ',I
		
	     READ(KEI,*) WAVEL(I),ALAMB(I),GLAMB(I,1),GLAMB(I,2),GLAMB(I,3),PHILAMB(I),PLAMB(I),      &
                    BLAMB(I),EZERO(I)
			
	  END DO   
	
      BLAMB = 0.0099	 ! B Clark updated Dec. 2017 Following ROSE et al., 2017
                                ! For a constant scattering for particles across
                                ! all wavelengths

   ELSE
        
      READ(KEI,*) (WAVEL(I),ALAMB(I),GLAMB(I,1),PHILAMB(I),PLAMB(I),      &
               BLAMB(I),EZERO(I),I=1,NWAVEL)
        
   ENDIF

#endif

	   
#if defined(PHOTODEGRADATION)	   
	 ! find where our visible light begins for when we want to integrate the PAR curve for the light dependcy of algae growth

	   DO I = 1,NWAVEL
	   
			IF (WAVEL(I) == 400)THEN
			
		     	PAR_MARK = I
			
			  EXIT
			
	        ENDIF 
	   ENDDO

	   NWAVEL_PAR=SIZE(WAVEL(PAR_MARK:NWAVEL))
	!    write(*,*)'NWAVEL_PAR = ',NWAVEL_PAR
	   LAMBDA_INT=WAVEL(2)-WAVEL(1)	   
	!   20 continue
	   
#endif

1032 FORMAT(/)
1036 FORMAT(//8X,2F8.1,I8,3F8.1)

1037 FORMAT(//8X,8F8.1)
10371 FORMAT(//8X,10F8.1)  ! B Clark for new input file

1038 FORMAT(//8X,8F8.1)
1082 FORMAT(//8X,2A8//)
1087 FORMAT(8X,4F8.0)

   END SUBROUTINE OWQ_READ
   
!=======================================================================
!			Allocate vars
!=======================================================================
   
 SUBROUTINE OWQ_ALLOC

	  ALLOCATE(TSS(MTLOC,KBM1));   TSS = 0.0

	  ALLOCATE(IAVG(0:MTLOC,KBM1));   IAVG   = 0.0	!WLong moved this from mod_wqm.F
		 
	  ALLOCATE(TOTCHL(MTLOC,KBM1)); TOTCHL = 0.0	! B Clark April 2016, made these allocatble (includine below)
	  ALLOCATE(IATBOT_out(MTLOC));IATBOT_out = 0.0;
	  
      ALLOCATE(ALAMB(NWAVEL)) 			
	  ALLOCATE(PLAMB(NWAVEL))
	  ALLOCATE(PHILAMB(NWAVEL))
	  ALLOCATE(BLAMB(NWAVEL))
      ALLOCATE(WAVEL(NWAVEL))	
      ALLOCATE(EZERO(NWAVEL))
	  ALLOCATE(EZ(NWAVEL))

                               	   
#if defined(PHOTODEGRADATION)
	  
	  ALLOCATE(ABSCDOM(MTLOC,KBM1,NWAVEL,3));ABSCDOM =0.0

	  ALLOCATE(KDLAMB_out(MTLOC,KBM1,NWAVEL)); KDLAMB_out = 0.0
	  ALLOCATE(KD_VEC(NWAVEL)); KD_VEC=0.

	  ALLOCATE(KDLAMB_int(0:MTLOC,KBM1)); KDLAMB_int = 0.0
	  ALLOCATE(IAVG_pd(MTLOC,KBM1,NWAVEL)); IAVG_pd = 0.0
	  ALLOCATE(IAVG_in(NWAVEL)); IAVG_in = 0.0
	  ALLOCATE(I0_lambda(NWAVEL)); I0_lambda=0.0
	 IF(PHOTODEG) THEN 
	  	  ALLOCATE(GLAMB(NWAVEL,3))
        ELSE
                
	    ALLOCATE(GLAMB(NWAVEL,1))
        ENDIF

          ALLOCATE(IATBOT(NWAVEL))
#endif
   END SUBROUTINE OWQ_ALLOC

!=======================================================================
!		 	Deallocate Vars
!=======================================================================
   
   SUBROUTINE OWQ_DEALLOC

		 IF(ALLOCATED(TSS)) 	DEALLOCATE (TSS)

		IF(ALLOCATED(IAVG)) 	DEALLOCATE (IAVG)	!WLong moved this to here from wqm_main.F
		
		
		IF(ALLOCATED(TOTCHL)) DEALLOCATE (TOTCHL)   ! B Clark see above
        IF(ALLOCATED(IATBOT_out)) DEALLOCATE(IATBOT_out)
        IF(ALLOCATED(ALAMB)) DEALLOCATE(ALAMB)
        IF(ALLOCATED(PLAMB)) DEALLOCATE(PLAMB)		
        IF(ALLOCATED(PHILAMB)) DEALLOCATE(PHILAMB)	
	    IF(ALLOCATED(BLAMB)) DEALLOCATE(BLAMB)
		IF(ALLOCATED(WAVEL)) DEALLOCATE(WAVEL)
		IF(ALLOCATED(EZERO)) DEALLOCATE(EZERO)
		IF(ALLOCATED(EZ)) DEALLOCATE(EZ)		
		
		
#if defined(PHOTODEGRADATION)

		 IF(ALLOCATED(ABSCDOM)) DEALLOCATE(ABSCDOM)
	    IF(ALLOCATED(KDLAMB_Out)) DEALLOCATE(KDLAMB_out)
		 IF(ALLOCATED(KD_VEC)) DEALLOCATE(KD_VEC)
       IF(ALLOCATED(KDLAMB_int)) DEALLOCATE(KDLAMB_int)

		 IF(ALLOCATED(IAVG_pd)) DEALLOCATE(IAVG_pd)
		 IF(ALLOCATED(IAVG_in)) DEALLOCATE(IAVG_in)
		 IF(ALLOCATED(I0_lambda)) DEALLOCATE(I0_lambda)
		 IF(ALLOCATED(GLAMB)) DEALLOCATE(GLAMB)
	    IF(ALLOCATED(IATBOT)) DEALLOCATE(IATBOT)	

		 
#endif
		
   END SUBROUTINE OWQ_DEALLOC

!=========================================================================

! Added in to remove the older function that was clunky
   
#if defined(PHOTODEGRADATION)    

 SUBROUTINE PHOTOCHEMISTRY  ! This is modified form the OWQ function to have vector operations on the whole matrices and solve for the CDOM1,2,3 absorbance

							! B Clark, April 2016
   IMPLICIT NONE
 
 !!!!!!!!!!!!===================!!!!!!!!!!!!!!  
  ! B Clark Commented out on June 2017;
  !updated specific absorption spectra from Rose,Neale, Tzortziou et al., (L&O in prep) that do not require the coefficients from Gallegos model
  ! Can now use the specific absorption and the concentrations of each absorbing constituent directly
  
!!!!!!!!!!===================!!!!!!!!!!!!!!!!

! write(*,*)'PHOTOCHEMISTY Loop Activated'
DO J=1,NWAVEL  
  DO K = 1,KBM1
    DO I=1,MLOC

   ! Updated code uses specific absorption curves and concentration of each absorbing compound
   	 
	 ABSCDOM(I,K,J,1) = GLAMB(J,1)*WC_CDOC1(I,K) !  ! get all 3 cdoms and add em up, want to keep the array for calculating the irradiance lost due to each CDOM
  
	 ABSCDOM(I,K,J,2) = GLAMB(J,2)*WC_CDOC2(I,K)

	 ABSCDOM(I,K,J,3) = GLAMB(J,3)*WC_CDOC3(I,K)
	 
								! Older version needed PHISTAR676, with new data we have the actually naperian absorption spectra...			
     ABCHL  = PHILAMB(J)*TOTCHL(I,K)  !*APHI676(I,K)       ! the total CHL is in g chla per m^3 so need to multiply PHILAMB (m^2 mg chla^-1) by 1000

								
     ABNAPM = PLAMB(J) * TSS(I,K) ! B Clark June 2017 !* AP440(I,K)  ! need to figure out the AP440 part, what is going on here?
	 
     ABTOT  = ALAMB(J) + ABCHL + ABNAPM +SUM(ABSCDOM(I,K,J,:))  ! get the total absorbance with all 3 cdoms

     SCAT  = BLAMB(J)*TSS(I,K)  ! B Clark June 2017  !*B440(I,K)

     ! B Clark updated to the Lee et al., model Feb 2020
! 	 ! B Clark put all this in th eloop so now the subrotuine is only called once and it calculates everything
 
   !   180*CTHTAA/pi
     !CTHTAA is the elevation of the sun relative to the horizon in radians
     !convert to degrees below then it is subtracted from 90° in KD_LEE to yield the zenith angle
    CALL KD_LEE(180.*CTHTAA/PI,ABTOT,SCAT,KDLAMB_out(I,K,J))
!       CALL KD_LEE(DECL*180.,ABTOT,SCAT,KDLAMB_out(I,K,J))

! COMPUTE REMAINING SPECTRAL IRRADIANCE

      EZ(J) = I0*EZERO(J) * EXP( -KDLAMB_out(I,K,J) * DEPTH )	 


    ENDDO  ! NWAVEL loop  
  ENDDO
ENDDO
!write(*,*)' now onto the new CDOM absorption parts'
! Still need to get a formula to convert SWRAD to UV-VIS energy
!write(*,*)'abscdom  = ' , ABSCDOM 
! DO J=1,NWAVEL
    I0_lambda= I0 * EZERO  ! This gives the spectral distribution of our shortwave radiation, before the model took par,
! END DO                      !we are using shortwave radiation now
   IAVG = 0.0	 ! set to zero before integration
   KDLAMB_int=0.0
 DO I=1,MLOC
 
   DO J = 1,NWAVEL    ! This will need to be spectrally dependent to get the photobleaching correct

     IATTOP = I0_lambda(J)*(1.-SREFLECT/100.)  ! get the surface irradiance

     OPTDEPTH = KDLAMB_out(I,1,J)*D(I)*DZ(1) !  optical depth for the surface cell from photochem subroutine

     IATBOT(J) = IATTOP*EXP(-OPTDEPTH)    ! ! irradiance left at the bottom of the cell, this gets carried to the next cell

     IAVG_pd(I,1,J) = (IATTOP-IATBOT(J))/OPTDEPTH ! average light level in the cell, assumption is the cell is vertically well mixed

         DO Idom= 1,3		
              ICDOM_abs(I,1,J,Idom) = IATTOP*ABSCDOM(I,1,J,Idom)/ KDLAMB_out(I,1,J)- & 												 
									 	        IATTOP*ABSCDOM(I,1,J,Idom)/ KDLAMB_out(I,1,J)*EXP(-OPTDEPTH)! W / m^2 / day
		  
		  	   if(ICDOM_abs(I,1,J,idom) < 0.0) THEN
				   ICDOM_abs(I,1,J,idom)=0.0
			   ENDIF
			   N_photons(I,1,j,Idom)  = (ICDOM_abs(I,1,j,Idom)*WAVEL(J)*5.03E15)/6.0221409E+23  
        END DO
	 END DO 
	 
	 		!	write(*,*)'ICDOM_abs = ',ICDOM_ABS(:,:,51,:)
	 
! get our PAR light remaining from 400:700 nm
  ! still need to convert the energy flux to number of photons so we can use in P vs I curve (Einstens/m^2)
	!https://www.berthold.com/en/bio/how-do-i-convert-irradiance-photon-flux
!	   IAVG(I,K) = 0.0	 
	 
	!watts/m^2 * lambda(nm) * [Energy / photon / lambda conversion * 86400 s / day[] / avogradros # = moles Photons/m^2/day
	 
	 !  The number of photons Np can be calculated by (Note: the nm value for λ is used):
	   ! h = Planck's constant, c = speed of light in vacuum
	 
     !  Np = E/Ep= E•((λ•10-9)/h•c) = E [W/m2]•λ•10-9[m]•/ (1.988•10-25) [J/s•m/s]
  !       =  E•λ•5.03•1015 [1/(m2•s)] 


	IAVG_in=IAVG_pd(I,1,:)*(WAVEL(:)*4.3459E20)/6.0221409E+23    
	   
	CALL simpson(IAVG_in(PAR_MARK:NWAVEL),1,NWAVEL_PAR,IAVG(I,1),NWAVEL_PAR,LAMBDA_INT)	  
        
	CALL simpson(KDLamb_out(I,1,PAR_MARK:NWAVEL),1,NWAVEL_PAR,KDLAMB_int(I,1),NWAVEL_PAR,LAMBDA_INT)	  

	KDLAMB_int(I,1)=KDLAMB_int(I,1)/(WAVEL(NWAVEL)-WAVEL(PAR_MARK)) ! Get the a

	  
    DO K=2,KBM1 ! loop through the rest of the grid vertically
	
	       DO J = 1,NWAVEL
            IATTOP = IATBOT(J)
            OPTDEPTH = KDLAMB_out(I,K,J)*D(I)*DZ(K)   ! now need to make this spectrally depenedent so that we get the spectral irradiance left over 
			
            IATBOT(J) = IATTOP*EXP(-OPTDEPTH)
			
            IAVG_pd(I,K,J) = (IATTOP-IATBOT(J))/OPTDEPTH 

              DO Idom = 1,3
		   
			! B Clark updated code to get light absorbed by each CDOM according to Hu et al., 2002, L&O
										  
                ICDOM_abs(I,K,J,Idom) = IATTOP*ABSCDOM(I,K,J,Idom)/ KDLAMB_out(I,K,J)- &
                                        IATTOP*ABSCDOM(I,K,J,Idom)/ KDLAMB_out(I,K,J)*EXP(-OPTDEPTH)! W / m^2 / day     
                                                  
				   if(ICDOM_abs(I,K,J,idom) < 0.0) THEN
				     ICDOM_abs(I,K,J,idom)=0.0
				   ENDIF   
				! mol Photons/m^2/s, still need to convert to mols Photons for use with our AQY in mod_wc_dom
				   N_photons(I,K,j,Idom)  = (ICDOM_abs(I,K,j,Idom)*WAVEL(J)*5.03E15 ) / 6.0221409E+23    
			  END DO
       END DO

	! ! updated by B Clark to use a simpson integration over the PAR wavelengths 
	! !    !when the wavebands are not spaced by 1 nm

	   IAVG_in=IAVG_pd(I,K,:)*(WAVEL(:)*4.3459E20)/6.0221409E+23    
	   
		CALL simpson(IAVG_in(PAR_MARK:NWAVEL),1,NWAVEL_PAR,IAVG(I,K),NWAVEL_PAR,LAMBDA_INT)	  
		  
		CALL simpson(KDLamb_out(I,K,PAR_MARK:NWAVEL),1,NWAVEL_PAR,KDLAMB_int(I,K),NWAVEL_PAR,LAMBDA_INT)	  

	   KDLAMB_int(I,K)=KDLAMB_int(I,K)/(WAVEL(NWAVEL)-WAVEL(PAR_MARK)) ! Get the a
	  
	  END DO  ! K

	END DO  ! I


END SUBROUTINE PHOTOCHEMISTRY
   
#endif
   
SUBROUTINE LGHT_ATTN
   

!vjp  10/11/04 initialize arrays OPZ, G1, and G2

   OPZ(1:9) = (/ 0.0, 0.01, 0.025, 0.033, 0.05, 0.1, 0.2, 0.286, 0.4 /)
   G1(1:9) =  (/ 0.425, 0.425, 0.364, 0.3383, 0.2996, 0.2175,           &
                 0.1232, 0.073167, 0.01546 /)
   G2(1:9) =  (/ 0.19, 0.19, 0.1396, 0.1187, 0.0875, 0.02288,           &
                -0.048, -0.08367, -0.12386 /)
                
! CONVERT MODEL JULIAN DAY TO DAY AND HOUR

   DOY = AMOD(JDAY,365.)
   HOUR = 24.*(DOY-AINT(DOY))
   DOY = AINT(DOY)+DOFFSET

! COMPUTE COSIN OF IN-WATER SOLAR ANGLE

! SOLAR DATE
   SOLDATE =  TWOPI*DOY/365.
  ! SOLDATE = 5
 !  write(*,*)'SOLDATE = ',SOLDATE

! DECLINATION
   DECL = TWOPI*(0.39637-(22.9133*COS(SOLDATE))+(4.02543*SIN(SOLDATE))-  &
         (0.3872*COS(2.*SOLDATE))+(0.052*SIN(2*SOLDATE)))/360.

! TIME OF DAY
   TAU = TWOPI*HOUR/24.
   
   
   ! B Clark fix the tau to be at Noon
!   TAU = TWOPI/2.
 !  write(*,*)'ALERT!!!!! TAU is FIXED TO HAVE THE SUN AT NOON!!!!!!'

!write(*,*)'Tau =',TAU
! COSIN OF THETA IN AIR
   CTHTAA=((SIN(LAT))*(SIN(DECL)))-((COS(LAT))*(COS(DECL))*(COS(TAU)))

! INCIDENCE OF THETA IN AIR
   ITHTAA = ACOS(CTHTAA)

! THETA IN WATER
! Wen Long, here 1.33 is the water refraction index of light
   THTAW=ASIN(SIN(ITHTAA)/1.33)

! REFLECTANCE
   IF (ITHTAA > TWOPI/4.) THEN
     SREFLECT = 100.
   ELSE
     SREFLECT = EXP(0.6148+ITHTAA**3)
   ENDIF

! SCATTERING FUNCTION
   AVREFCOS = COS(THTAW)
   GMUZERO(1) = 0.198
   DO I=2,9
     GMUZERO(I) = AVREFCOS*G1(I)-G2(I)
   ENDDO

! EVALUATE ATTENUATION FOR EACH SURFACE CELL

  DO K = 1, KBM1 ! B Clark added in vertical distribution of differential KD values Sep 2015
   DO I=1,MLOC
   
     TOTCHL(I,K) = (B1(I,k)*CHLCMN1 + B2(I,k)*CHLCMN2 + B3(I,k)*CHLCMN3)*TCHL2CHL *1000 ! convert grams of Chl to mg Chl


    !calculate total
    TOTPOC   = B1(I,k)+B2(I,k)+B3(I,k)+LPOC(I,k)+RPOC(I,k)
    
#if defined(ZOOPLANKTON)
     TOTPOC=TOTPOC+SZ(I,k)+LZ(I,k);
#endif
     TSS(I,K)    = MAX(SSI(I,k) + 2.5*TOTPOC, 1.0) 
	 

	END DO 
  END DO
	
#if defined(PHOTODEGRADATION)

  ! B Clark Commented out on June 2017;
  !updated specific absorption spectra from Rose,Neale, Tzortziou et al., (L&O in prep) that do not require the coefficients from Gallegos model
  ! Can now use the specific absorption and the concentrations of each absorbing constituent directly
!  write(*,*)'Calling Photochemistry'
  CALL PHOTOCHEMISTRY 
!  write(*,*)'Done Calling Photochemistry'


#endif
     
   
   RETURN
   END SUBROUTINE LGHT_ATTN

  SUBROUTINE KD_LEE(SZA,a_total,beta_tss,KD_out)
    !This is the Kd light attenuation model based on Lee et al.,
	!SZA is the solar zenith angle in degrees, can take the declination value and multiply by 180
	!a_total is the total absorption form all in water consituents
	!beta_tss is the scattering due to total suspended sediment
    REAL(SP) :: m0=0.0
	REAL(SP) :: m1=4.259
	REAL(SP) :: m2=0.52
    REAL(SP) :: m3=-10.80
    real(sp) ::m4,m5,m6!these are used for the inside of paranthesis

    REAL(SP), INTENT(IN):: SZA
    REAL(SP), INTENT(IN):: a_total
    REAL(SP), INTENT(IN):: beta_tss
    REAL(SP), INTENT(OUT):: KD_out

    !Calculate KD
        m0=1.0+0.005*(90-SZA)
        m4=m3*a_total
        m5=1.-m2*exp(m4)
       !Kd = m0.* a_total + m1 * (1.0 - m2 * exp( m3 .* a_total)) .* beta_tss;

    KD_out=m0 * a_total + m1 * m5 * beta_tss

	RETURN

END SUBROUTINE KD_LEE

   
END MODULE MOD_OWQ
