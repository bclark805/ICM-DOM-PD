MODULE mod_ncdio
!==============================================================================!
!  NetCDF Output for FVCOM using CF Metadata Convention                        !
!                                                                              !
!    see: http://www.cgd.ucar.edu/cms/eaton/cf-metadata/ for info              !
!                                                                              !
!    current time dependent variables set up                                   !
!         el:    surface elevation                                             !
!          u:    x-velocity. In spherical coordinate,lon-velocity              !                         
!          v:    y-velocity. In spherical coordinate,lat-velocity              !                        
!         ww:    z-velocity                                                    !
!         kh:    turbulent diffusivity                                         !
!         km:    turbulent viscosity                                           !
!         t1:    temperature                                                   !
!         s1:    salinity                                                      !
!         ua:    vertically-averaged x-velocity                                !
!                In spherical coordinate,vertically-averaged lon-velocity      !
!         va:    vertically-averaged y-velocity                                !
!                In spherical coordinate,vertically-averaged lat-velocity      !
!          d:    depth at nodes                                                !
!        dye:    dye at nodes                                                  !
!       aice:    ice concentration on nodes
!       vice:    ice thichness on nodes
!      uuice:    ice x-velocity 
!      vvice:    ice y-velocity   
!     uuwind:    wind speed in x direction
!     vvwind:    wind speed in y direction

!       wd:      wet/dry flag (0 or 1)

!    Two-dimensional varialbles used in Sed_diag or NPZD                      !
!        amm:   ammonia in N units (mgN/L) from mod_npzd.F                     !
!        doo:   dissolved oxygen in O2 units (mgO2/m^3 --> mgO2/L by /1000) by using doxyg from mod_npzd.F           !
!		 phy
!		 zoo
!  		 dtn
! 		 dtc
!		 nit



!    One-dimensional varialbles used in Sed_diag                             									    !
!        jcin:  flux of particulate organic phosphorus form POC (mgC/m^2/d --> mgC/m^2/d) by aggregate CPOC in 2nd dimention from mod_sed_diag.F     !
!        jnin:  flux of particulate organic phosphorus form PON (mgC/m^2/d --> mgC/m^2/d) by aggregate CPON in 2nd dimention from mod_sed_diag.F      !
!        jpin:  flux of particulate organic phosphorus form POP (mgC/m^2/d --> mgC/m^2/d) by aggregate CPOP in 2nd dimention from mod_sed_diag.F     !
!        jsin:  flux of particulate organic silica(mgC/m^2/d) from mod_sed_diag.F                   !

!        NH40:  dissolved ammonia (mgN/L)in overlying water column n sigmalayer KBM1                !
!        NO30:  dissolved NO3 		(mgN/L)in overlying water column in sigmalayer KBM1		        !
!        PO40:  dissolved PO4 		(mgN/L)in overlying water column in sigmalayer KBM1             !  
!        SI0 :  dissolved silicate  (mgN/L)in overlying water column in sigmalayer KBM1             !
!        CH40:  dissolved methane   (mgN/L)in overlying water column in sigmalayer KBM1             !  
!                                               
!        SOD:   sediment oxygen demand (mgO2/m^2/d)                            !  
!
!        JNH4:  flux of dissolved ammonia (mgN/m^2/d--> gN/m^2/d by/1000) back to the overlying water column           !  
!        JNO3:  flux of dissolved nitrate (mgN/m^2/d--> gN/m^2/d by/1000) back to the overlying water column           !
!        JDenitT: flux of denitrification (gN/m^2/d --> gN/m^2/d by/1) back to the overlying water column, equal to BENDEN from mod_sed_diag.F             ! 
!        JCH4:    flux of dissolved methane (gO2/m^2/d) back to the overlying water column            !
!        JCH4G:   flux of geseous methane (gO2/m^2/d) back to the overlying water column              !
!        JHS:  	  flux of dissolved HS  (gO2/m^2/d) back to the overlying water column             ! 
!        JPO4:    flux of dissolved PO4 (mgP/m^2/d--> gP/m^2/d by/1000) back to the overlying water column             ! 
!        JSI:  	  flux of dissolved silicate (mgSi/m^2/d--> gSi/m^2/d by/1000) back to the overlying water column       ! 

!        NH4(1):  dissolved ammonia in the 1st layer of sediment (mgN/m^3--> mgN/L(=gN/m^3) by/1000) equals to NH41TM1S (previous time step) in mod_sed_diag.F   ! 
!        NH4(2):  dissolved ammonia in the 2nd layer of sediment (mgN/m^3--> mgN/L by/1000) equals to NH41TM1S (previous time step) in mod_sed_diag.F   ! 
!        NO3(1):  dissolved nitrate in the 1st layer of sediment (mgN/m^3--> mgN/L by/1000)     !
!        NO3(2):  dissolved nitrate in the 2nd layer of sediment (mgN/m^3--> mgN/L by/1000)     !   
!        PO4(1):  dissolved PO4     in the 1st layer of sediment (mgP/m^3--> mgP/L by/1000) equals to PO41TM1S (previous time step) in mod_sed_diag.F   !
!        PO4(2):  dissolved PO4     in the 2nd layer of sediment (mgP/m^3--> mgP/L by/1000) equals to PO42TM1S (previous time step) in mod_sed_diag.F   ! 
!        Si(1):   dissolved silicate in the 1st layer of sediment (mgSi/m^3--> mgSi/L by/1000) equals to SI1TM1S (previous time step) in mod_sed_diag.F   ! 
!        Si(2):   dissolved silicate in the 2nd layer of sediment (mgSi/m^3--> mgSi/L by/1000) equals to SI2TM1S (previous time step) in mod_sed_diag.F   ! 
!        CH4(1):  dissolved methane in the 1st layer of sediment (mgC/m^3--> mgC/L by/1000) equals to CH41TM1S (previous time step) in mod_sed_diag.F   ! 
!        CH4(2):  dissolved methane in the 2nd layer of sediment (mgC/m^3--> mgC/L by/1000) equals to CH42TM1S (previous time step) in mod_sed_diag.F   !  
!        HS(1):   dissolved HS in the 1st layer of sediment (mgO2/m^3--> mgO2/L by/1000) equals to HS1TM1S (previous time step) in mod_sed_diag.F   ! 
!        HS(2):   dissolved HS in the 2nd layer of sediment (mgO2/m^3--> mgO2/L by/1000) equals to HS2TM1S (previous time step) in mod_sed_diag.F   !  
!        POC2G1: particular organic carbon     (POC) for G1 class in 2nd layer of sediment (mgC/m^3 --> gO2/m^3 by/1000.0*2.667)  see CPOC(I,1) in mod_sed_diag.F  ! 
!        POC2G2: particular organic carbon     (POC) for G2 class in 2nd layer of sediment (mgC/m^3 --> gO2/m^3 by/1000.0*2.667)  see CPOC(I,2) in mod_sed_diag.F  ! 
!        POC2G3: particular organic carbon     (POC) for G3 class in 2nd layer of sediment (mgC/m^3 --> gO2/m^3 by/1000.0*2.667)  see CPOC(I,3) in mod_sed_diag.F  ! 
!        PON2G1: particular organic nitrogen   (PON) for G1 class in 2nd layer of sediment (mgN/m^3 --> gN/m^3 by/1000) refer to CPON! 
!        PON2G2: particular organic nitrogen   (PON) for G2 class in 2nd layer of sediment (mgN/m^3 --> gN/m^3 by/1000)  !
!        PON2G3: particular organic nitrogen   (PON) for G3 class in 2nd layer of sediment (mgN/m^3 --> gN/m^3 by/1000)  ! 
!        POP2G1: particular organic phosphorus (POP) for G1 class in 2nd layer of sediment (mgP/m^3 --> gP/m^3 by/1000) refer to CPOP! 
!        POP2G2: particular organic phosphorus (POP) for G2 class in 2nd layer of sediment (mgP/m^3 --> gP/m^3 by/1000)  !   
!        POP2G3: particular organic phosphorus (POP) for G3 class in 2nd layer of sediment (mgP/m^3 --> gP/m^3 by/1000)  !    
!        POS2:    particular biogenic selicate  (POS)              in 2nd layer of sediment (mgSi/m^3 --> gSi/m^3 by/1000) refer to CPOS!  
!        H1:  	  sediment layer thickness (m --> mm by*1000) refer to HSED1 in mod_sed_diag ! 
!        BENSTR: bentic stress (deimensionless) refer to BENSTRTM1S in mod_sed_diag.F       !
!                 Benthic shear stress regulates the convection of volatile substances such as oxygen across the sediment�water interface. To be cont'd !
!        Cont'd   Benthic shear stress is the primary determinant of where particulate matter scours, rolls, saltates, settles, and consolidates. refer to Perterson 1999 [doi:10.1016/S0144-8609(99)00025-4]! 

!    to add additional variables:                                              !
!      1.) add to list above                                                   !
!      2.) add *_vid to variables vid in section "new variable vid"            !
!      3.) go to definition section "new variable definition"                  !
!      4.) add output section "new variable output"                            !
!==============================================================================!


! NOTES:
!
! David: NV was written as a real -> changed to integer
! David: added interface to PUTVAR to handle integers, reals, one and
  ! two dimensional variables. Now we can properly write values in
  ! their native format. 1/19/07

!#if defined (NETCDF_IO)

   USE MOD_INP, ONLY: SCAN_FILE
   
   !USE mod_prec
   
   USE netcdf  !netcdf must be compatiable with the functions here and also
               !compiled using the same compiler
               !netcdf is written in C, gcc or intel c to compile
               !
    USE MOD_CONTROL, ONLY :         &
            SERIAL          &           !!TRUE IF SINGLE PROCESSOR
            ,MSR            &           !!TRUE IF MASTER PROCESSOR (MYID==1)
            ,PAR            &           !!TRUE IF MULTIPROCESSOR RUN
            ,CASENAME,CASETITLE      !    &           !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
			
			
	USE MOD_HYDROVARS, ONLY: UU, VV, EL, &
							 DZ, VX, VY, VXMIN, VYMIN, H,Z, &
							 SIGLEV,SIGLAY
							 
	!USE MOD_HYDRO, ONLY: WWL
	
	 USE MOD_PREC, ONLY: SP
	 
    USE MOD_SIZES, ONLY: NCP, NSSFP, MGL,NGL
	
    USE MOD_LIMS, ONLY: MYID, NPROCS, IINT, NTLOC, MLOC,KB,KBM1,NLOC
	
    USE MOD_TGE, ONLY: NV, TGE_DEALLOC
	
	USE MOD_WQM, ONLY: T, SALT, &
						DOXG, &
						total_NetPP, &
						B1, B2, &
						NH4, NO3, PO4, &
						LPOC, RPOC, &
						LPON, RPON, &
						SSI, RHO1, ShearStress, &
						!nitrogen dt
						ALGDON, &
						ALGPON, &
						ALGNH4,&
						ALGNO3, &
						DENNO3,  &
						NT,  &
						DENITN,  &
						HDRLPON, &
						HDRRPON, &

						!carbon dt
						ALGDOC, &
						ALGPOC, &
						DENIT, &
						HDRLPOC, &
						HDRRPOC, WW_OUT

#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY:   &      
            !EL_PID,    &        !(:)           !!PROCESSOR OWNER OF GLOBAL ELEMENT
            !ELID,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT
            !NLID,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE 
            !ELID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
            !NLID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
            !OBN_GL2LOC,    &    !(:)           !!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
            EGID,        &    !(:)           !!GLOBAL ID OF LOCAL ELEMENT
            NGID,        &    !(:)           !!GLOBAL ID OF LOCAL NODE 
            !NHE,        &   !            !NUMBER OF HALO ELEMENTS
            !HE_LST,        &    !(:)           !!GLOBAL IDENTITIES OF HALO ELEMENTS
            !HE_OWN,        &    !(:)           !!OWNER OF HALO ELEMENTS
            !NBN ,       &   !              !!NUMBER OF BOUNDARY NODES
            !MX_MLT,     &   !              !!MAX MULTIPLICITY OF BOUNDARY NODES
            !BN_LST,     &    !               !!GLOBAL IDENTITY OF BOUNDARY NODES
            !BN_LOC,     &    !            !!LOCAL IDENTITY OF BOUNDARY NODES
            !BN_MLT,     &    !            !!MULTIPLICITY OF BOUNDARY NODES
            !BN_NEY,     &    !            !!NODE OWNER LIST
            !NDE_ID,     &    !            !! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
            !NHN,          &    !            !!NUMBER OF HALO NODES
            !HN_LST,     &    !            !!LIST OF HALO NODES 
            !HN_OWN,     &    !            !!PRIMARY OWNER OF HALO NODES
            EC,            &    !
            NC,            &    !
            BNC,         &    !
            EMAP,        &    !
            NMAP,         &    !
            !PNE,        &    !(:)        !!NUMBER OF ELEMENTS IN EACH PROC
            !PNN,        &    !(:)        !!NUMBER OF NODES IN EACH PROC
            !PNHE,        &    !(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
            !PNBN,        &    !(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
            !PMBM,        &    !(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
            !PNHN,        &    !(:)       !!NUMBER OF HALO NODES IN EACH PROC
            !FILELOCK,   & !
            !MYFILELOCK, & !
            !FILELOCKED , & !
            !NODE_MATCH, & !
             NODE_MATCH_ALL, &  !
            !EXCHANGE, & !
             EXCHANGE_ALL, & !
             GATHER , & !
            !IGATHER, & !
             INIT_MPI_ENV,&!
            !SORT,    & !
            !PREAD,    & !
            !PWRITE,& !
            !PPRINT, & !
            !LOCK, & !
            !UNLOCK, & !
            !GETLOC,  &!
             PAR_DEALLOC
#  endif						
						

						
	
	USE WC_DOM, ONLY: WC_CDOC1, WC_CDOC2, WC_CDOC3, &
             WC_NCDOC1, WC_NCDOC2, WC_NCDOC3, &

            XFLUX_CDOC1, XFLUX_CDOC2, XFLUX_CDOC3, &
            XFLUX_NCDOC1, XFLUX_NCDOC2, XFLUX_NCDOC3, &

					WC_CDON1, WC_CDON2, WC_CDON3, &
             WC_NCDON1, WC_NCDON2, WC_NCDON3, & 
             
					WC_CDOP1, WC_CDOP2, WC_CDOP3, &
				   	WC_NCDOP1, WC_NCDOP2, WC_NCDOP3, &
						COAGN, & 
						MNLDON1, &
						MNLDON2, &
						MNLDON3, &
						COAGC, &
						MNLDOC1, &
						MNLDOC2, &
						MNLDOC3, &
						PHOTOCHEM_DOC_ACCUM
					
					
	USE MOD_SED, ONLY:  SODTM1S, JNH4TM1S, JNO3TM1S, &
						BENDEN, JCH4TM1S, JCH4GTM1S, &
						JHSTM1S, JPO4TM1S, &
						NH41TM1S, NH42TM1S, &
						NO31TM1S, NO32TM1S,  &
						PO41TM1S, PO42TM1S, &
						CH41TM1S, CH42TM1S, &
						HS1TM1S, & HS2TM1S, &

						HSED1, BENSTR, & 
						JPOC, JPON, JPOP, &
						CPOC, CPON, CPOP
						
	USE MOD_SED_DOM, ONLY: CDOC1, CDOC2, CDOC3, &
						   NCDOC1, NCDOC2, NCDOC3, &
						   CDON1, CDON2, CDON3, &
						   NCDON1, NCDON2, NCDON3, &
						   CDOP1, CDOP2, CDOP3, &
						   NCDOP1, NCDOP2, NCDOP3, & 
						   JDOM_OUT
						   
	USE MOD_OWQ, ONLY: IAVG, KDLAMB_int
	
						   
   implicit none
   save

!--Control Variables----------------------------------------------!
   logical:: cdf_out            !!true to activate netcdf input/output
   integer :: nout_vars          !!number of variables to output
   integer :: cdf_int            !!output every cdf_int seconds (used to be intervals)
   integer :: cdf_stk            !!cdf_stk outputs are put in each file
   integer :: cdf_start  !! starting day to output netcdf
   integer :: cdf_end    !! ending day to output netcdf
   !!CDF_STK=0: ALL OUTPUTS IN SINGLE FILE
   integer :: stck_cnt           !!counts number of outputs in each file
   integer :: out_cnt            !!counts number of outputs
   real(sp) :: CDF_COUNTER  ! counter in seconds to keep track of the time for outputting netcdf
   character(len=120) :: cdfname !!netcdf file name
   character (len=80) ,allocatable, dimension(:) :: cdf_vdp
	character(len=256) :: outdir	
	integer :: n,m

!--NetCDF IDs----------------------------------------------------!

   !--NetCDF File 
   integer,private:: nc_ofid 

   !--Dimensions
   integer,private :: node_did,nele_did
   integer,private :: scl_did,siglay_did,siglev_did
   integer,private:: three_did,four_did
   integer,private :: time_did

   !--Grid Variables
   integer,private :: nprocs_vid,partition_vid
   integer,private :: idens_vid
   integer,private :: x_vid,y_vid,lat_vid,lon_vid
   integer,private :: nv_vid,nbe_vid
   integer,private :: aw0_vid,awx_vid,awy_vid
   integer,private :: a1u_vid,a2u_vid
   integer,private :: siglay_vid,siglev_vid,siglay_shift_vid
  
   !--Flow Variables 
   integer ,private:: time_vid
   integer,private:: iint_vid

   integer,private :: s1_vid
   integer,private :: t1_vid
   integer,private :: el_vid
   integer,private :: h_vid

   integer,private :: d_vid
   integer,private :: dz_vid
   !Velocity
   integer, private ::uu_vid
   integer, private ::vv_vid
   integer, private ::ww_vid
   integer, private ::wts_vid

		integer,private :: jcin_vid     
		integer,private :: jnin_vid     
		integer,private :: jpin_vid   
		integer,private :: jsin_vid  
		integer,private :: O20_vid
		integer,private :: NH40_vid    
		integer,private :: NO30_vid
   
		!integer,private :: PO40_vid    
		!integer,private :: SI0_vid 
		!integer,private :: CH40_vid    
   
		integer,private :: SOD_vid 
		integer,private :: JNH4_vid  
		integer,private :: JNO3_vid
		integer,private :: JDenitT_vid
		integer,private :: JCH4_vid
		integer,private :: JCH4G_vid
		integer,private :: JHS_vid  
		integer,private :: JPO4_vid
		integer,private :: JSI_vid
   
		integer,private :: NH41_vid
		integer,private :: NH42_vid
		integer,private :: NO31_vid
		integer,private :: NO32_vid
		integer,private :: PO41_vid
		integer,private :: PO42_vid
		integer,private :: SI1_vid
		integer,private :: SI2_vid
		integer,private :: CH41_vid
		integer,private :: CH42_vid
		integer,private :: HS1_vid
		integer,private :: HS2_vid
		integer,private :: POC2G1_vid
		integer,private :: POC2G2_vid
		integer,private :: POC2G3_vid
		integer,private :: PON2G1_vid
		integer,private :: PON2G2_vid
		integer,private :: PON2G3_vid
		integer,private :: POP2G1_vid
		integer,private :: POP2G2_vid
		integer,private :: POP2G3_vid   
		integer,private :: POS2_vid 
		integer,private :: H1_vid
		integer,private :: BENSTR_vid
   
!#		if defined (SED_DOM)
		integer,private :: SED_CDOC1_vid  ! The Colored G1 DOC  gC/m^3  ~ mgC/L
		integer,private :: SED_CDOC2_vid  ! The colored G2 DOC  gC/m^3  
		integer,private :: SED_CDOC3_vid  ! The colored G3 DOC  gC/m^3

		integer,private :: SED_NCDOC1_vid ! The non-colored G1 DOC gC/m^3
		integer,private :: SED_NCDOC2_vid ! The non-colored G2 DOC gC/m^3
		integer,private :: SED_NCDOC3_vid ! The non-colored G3 DOC gC/m^3

		integer,private :: SED_CDON1_vid  ! The Colored G1 DON  gN/m^3
		integer,private :: SED_CDON2_vid  ! The colored G2 DON  gN/m^3
		integer,private :: SED_CDON3_vid  ! The colored G3 DON  gN/m^3

		integer,private :: SED_NCDON1_vid ! The non-colored G1 DON gN/m^3
		integer,private :: SED_NCDON2_vid ! The non-colored G2 DON gN/m^3
		integer,private :: SED_NCDON3_vid ! The non-colored G3 DON gN/m^3

		integer,private :: SED_CDOP1_vid  ! The Colored G1 DOP  gP/m^3
		integer,private :: SED_CDOP2_vid  ! The colored G2 DOP  gP/m^3
		integer,private :: SED_CDOP3_vid  ! The colored G3 DOP  gP/m^3

		integer,private :: SED_NCDOP1_vid ! The non-colored G1 DOP gP/m^3
		integer,private :: SED_NCDOP2_vid ! The non-colored G2 DOP gP/m^3
		integer,private :: SED_NCDOP3_vid ! The non-colored G3 DOP gP/m^3
!#		endif   

		integer,private :: WC_NCDOC1_vid
		integer,private :: WC_NCDOC2_vid
		integer,private :: WC_NCDOC3_vid
		integer,private :: WC_CDOC1_vid
		integer,private :: WC_CDOC2_vid
    integer,private :: WC_CDOC3_vid
    
    integer,private :: XFLUX_NCDOC1_vid
		integer,private :: XFLUX_NCDOC2_vid
		integer,private :: XFLUX_NCDOC3_vid
		integer,private :: XFLUX_CDOC1_vid
		integer,private :: XFLUX_CDOC2_vid
		integer,private :: XFLUX_CDOC3_vid
									
		integer,private :: WC_NCDON1_vid
		integer,private :: WC_NCDON2_vid
		integer,private :: WC_NCDON3_vid
		integer,private :: WC_CDON1_vid
		integer,private :: WC_CDON2_vid
    integer,private :: WC_CDON3_vid   
    
    integer,private :: XFLUX_NCDON1_vid
		integer,private :: XFLUX_NCDON2_vid
		integer,private :: XFLUX_NCDON3_vid
		integer,private :: XFLUX_CDON1_vid
		integer,private :: XFLUX_CDON2_vid
		integer,private :: XFLUX_CDON3_vid   

		integer,private :: WC_NCDOP1_vid
		integer,private :: WC_NCDOP2_vid
		integer,private :: WC_NCDOP3_vid
		integer,private :: WC_CDOP1_vid
		integer,private :: WC_CDOP2_vid
		integer,private :: WC_CDOP3_vid

!#	if defined (SED_DOM)
		integer, private :: JWCDOC1_vid
		integer, private :: JWCDOC2_vid
		integer, private :: JWCDOC3_vid
		integer, private :: JWNCDOC1_vid
		integer, private :: JWNCDOC2_vid
		integer, private :: JWNCDOC3_vid
		
		integer, private :: JWCDON1_vid
		integer, private :: JWCDON2_vid
		integer, private :: JWCDON3_vid
		integer, private :: JWNCDON1_vid
		integer, private :: JWNCDON2_vid
		integer, private :: JWNCDON3_vid		
		
		integer, private :: JWCDOP1_vid
		integer, private :: JWCDOP2_vid
		integer, private :: JWCDOP3_vid
		integer, private :: JWNCDOP1_vid
		integer, private :: JWNCDOP2_vid
		integer, private :: JWNCDOP3_vid		
!#   endif	

  ! Below is added by B Clark for offline ICM model output for many many variables
		integer, private :: DOXG_vid	
		integer, private :: NPP_vid	
		integer, private :: B1_vid	
		integer, private :: B2_vid	
		integer, private :: NH4_vid	
		integer, private :: NO3_vid	
		integer, private :: PO4_vid	
		integer, private :: PAR_vid	
		integer, private :: LPOC_vid	
		integer, private :: RPOC_vid
		integer, private :: KD_vid	
		integer, private :: SSI_vid	
		integer, private :: RHO1_vid	
		integer, private :: Tau_vid	

		! Nitrogen 
		integer, private :: ALGDON_vid	
		integer, private :: ALGPON_vid	
		integer, private :: ALGNH4_vid
		integer, private :: ALGNO3_vid	
		integer, private :: NT_vid	
		integer, private :: DENNO3_vid
		integer, private :: DENITN_vid	
		integer, private :: HDRLPON_vid	
		integer, private :: HDRRPON_vid	
		integer, private :: COAGN_vid	
		integer, private :: MNLDON1_vid	
		integer, private :: MNLDON2_vid	
		integer, private :: MNLDON3_vid	
		integer, private :: LPON_vid	
		integer, private :: RPON_vid	
		
		!integer, private :: DTWCDON1_vid	
		!integer, private :: DTWCDON2_vid
		!integer, private :: DTWCDON3_vid
		!integer, private :: DTWNCDON1_vid	
		!integer, private :: DTWNCDON2_vid
		!integer, private :: DTWNCDON3_vid
		
		!integer, private :: DTNH4_vid	
		!integer, private :: DTNO3_vid
		!integer, private :: DTLPON_vid
		!integer, private :: DTRPON_vid
		
		! Carbon 
		integer, private :: ALGDOC_vid
		integer, private :: ALGPOC_vid
		integer, private :: DENIT_vid		
		integer, private :: HDRLPOC_vid	
		integer, private :: HDRRPOC_vid
		integer, private :: COAGC_vid	
		integer, private :: MNLDOC1_vid	
		integer, private :: MNLDOC2_vid	
		integer, private :: MNLDOC3_vid	
		
		!integer, private :: DTWCDOC1_vid	
		!integer, private :: DTWCDOC2_vid
		!integer, private :: DTWCDOC3_vid
		!integer, private :: DTWNCDOC1_vid	
		!integer, private :: DTWNCDOC2_vid
		!integer, private :: DTWNCDOC3_vid	
		
!		integer, private :: DTLPOC_vid
!		integer, private :: DTRPOC_vid
!!Photochemistry		
!		integer, private :: PD32_vid
		integer, private :: PD31_vid
		integer, private :: PD30_vid
	    integer, private :: PD3N_vid
		
		integer, private :: PD21_vid
		integer, private :: PD20_vid
		integer, private :: PD2N_vid	
		
		integer, private :: PD10_vid
		integer, private :: PD1N_vid		
		
		integer :: ipt		
		
		

   !--Info Variables
   character(len=256),public :: institution
   character(len=256),public :: netcdf_timestring 

   INTERFACE PUTVAR 
      MODULE PROCEDURE PUTVAR1D_INT
      MODULE PROCEDURE PUTVAR1D_REAL
      MODULE PROCEDURE PUTVAR2D_INT
      MODULE PROCEDURE PUTVAR2D_REAL
   END INTERFACE

   contains !------------------------------------------------------------------!
            ! handle_ncerr        :   deal with netcdf error                   !
            ! set_ncd_io          :   read assimilation parameters from input  !
            ! write_netcdf_setup  :   set up dimensioning and write grid       !
            ! out_netcdf          :   write time dependent data                ! 
            ! putvar              :   collect variable to global and dump      ! 
            ! -----------------------------------------------------------------!

!==============================================================================|
!==============================================================================|

!------------------------------------------------------------------------------|
!  CHECK NETCDF ERROR AND WRITE MESSAGE TO FILE UNIT IPT                       |
!------------------------------------------------------------------------------|
   SUBROUTINE handle_ncerr(status,errmsge,ipt)
   integer, intent(in) :: status,ipt
   character(len=*)    :: errmsge
   if(status /=nf90_noerr)then
     write(ipt,*)trim(errmsge)
     write(ipt,*)trim(nf90_strerror(status))
     call pstop
   end if
   END SUBROUTINE handle_ncerr

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
!  READ IN PARAMETERS CONTROLLING INPUT/OUTPUT FROM RUNTIME PARAM FILE         |
!==============================================================================|
  SUBROUTINE set_ncd_io   
   use mod_prec
  ! use all_vars
   use mod_inp
   use netcdf
   implicit none
!--Local Vars----------------------------|
   real(sp)           :: realvec(150)
   integer            :: intvec(150)
   integer            :: iscan,ingood
   character(len=120) :: fname
   character(len=80), dimension(100) :: charvec
   integer            :: i
!----------------------------------------|
!==============================================================================!										  
!! Set up IPT file 										  
!==============================================================================!
	
   write(*,*)'Setting Up Netcdf FIles'
	OPEN(unit=9832,file = 'netcdf_check.dat', iostat=ingood)								  
      if (ingood /= 0) THEN
	     write(*,*) 'Problem opening the Netcdf Check file'
	     write(*,*) ' error code is ',ingood
	     write(*,*) 'Exiting'
	    stop
	   ELSE
         write(*,*)'opened netcdf check file file'
	ENDIF									  
										  
	ipt=9832

   outdir='outputs'
   out_cnt = 0

   fname = "./"//trim(casename)//"_netcdf.dat"

!------------------------------------------------------------------------------|
!   cdf_out: netcdf activation flag        
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_OUT",LVAL = CDF_OUT)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_out: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  cdf_int: output is performed every cdf_int iterations
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_INT",ISCAL = CDF_INT)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_int: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  cdf_stk: # dumps / file                                
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_STK",ISCAL = CDF_STK)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_stk: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if
   
!------------------------------------------------------------------------------|
!  cdf_start cdf_end                             
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_START",ISCAL = CDF_START)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_start: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if
   
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_END",ISCAL = CDF_END)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_end: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if
      
!------------------------------------------------------------------------------|
!     cdf_vdp: list of variables to write to output file
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_VDP",CVEC = CHARVEC,NSZE = NOUT_VARS)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_vdp: ',iscan
     call pstop
   end if
   if(nout_vars <= 0)then
     write(ipt,*)'incorrect number of netcdf cdf_vdp variables specified'
     write(ipt,*)'in input file',nout_vars
     call pstop
   end if

   allocate(cdf_vdp(nout_vars))

   cdf_vdp(1:nout_vars)= charvec(1:nout_vars)

   
   
   
!------------------------------------------------------------------------------|
!            SCREEN REPORT OF SET VARIABLES                                    !
!------------------------------------------------------------------------------|
   if(msr)then
     write(ipt,*)''
     write(ipt,*)'!        netcdf parameters                  '
     if(cdf_out)then
       write(ipt,*)'!  netcdf i/o            :  active'
       write(ipt,*)'!  output every # its    : ',cdf_int
       write(ipt,*)'!  # dumps / file        : ',cdf_stk
	   write(ipt,*)'!  # Starting Day        : ',cdf_start
	   write(ipt,*)'!  # Ending Day          : ',cdf_end
       write(ipt,*)'!  # variables to write  : ',nout_vars
       do i=1,nout_vars
         write(ipt,999)i,trim(cdf_vdp(i))
       end do
     else
       write(ipt,*)'!  # netcdf i/o          :  not active'
     end if
   end if


   return
   999 format(' !  variable #',i4,'        :',a13)
   END SUBROUTINE set_ncd_io  
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
!  Write NetCDF Header and Static Variables                                    |
!==============================================================================|
   SUBROUTINE write_netcdf_setup(filecnt) 

   !use all_vars

#  if defined (MULTIPROCESSOR)
   use mod_par 
#  endif
   use netcdf
   use mod_types
   use mod_utils
   implicit none
   integer, intent(in)   :: filecnt
   integer, dimension(3) :: dynm3de_lev,dynm3de_lay
   integer, dimension(3) :: dynm3dn_lev,dynm3dn_lay
   integer, dimension(2) :: stat3de_lev,stat3de_lay 
   integer, dimension(2) :: stat3dn_lev,stat3dn_lay 
   integer, dimension(2) :: specdim
   integer, dimension(2) :: dynm2de,dynm2dn
   integer, dimension(1) :: stat2de,stat2dn
   integer, dimension(1) :: stat_lev,stat_lay,dynmtime ,stat_scl
   character(len=100)    :: netcdf_convention
   character(len=100)    :: timestamp ,temp
   integer               :: i,j,ierr,i1,i2,ipt,m
   integer               :: maxnode,maxnodep,maxelem,maxelemp,itmp
   real(sp), allocatable :: tmp(:,:),tvec(:)
   integer, allocatable  :: tmpint(:,:)
   character(len=4)      :: nchar
	character(len=256) :: fvcom_version='2.7.6'
	character(len=256) :: fvcom_website='fvcom_website'					

!==============================================================================|

!==============================================================================|
!  Set up Constants and Initialize Counters                                    |
!==============================================================================|

!--Initialize Stack Count
   stck_cnt = 1

!--NetCDF Convention String
   netcdf_convention = 'CF-1.0'

!--Time Stamp for History
   call get_timestamp(temp)
   timestamp = 'model started at: '//trim(temp)
										  
										  
   m=mloc
   n=nloc

!============================================================|
!  OPEN FILE AND DEFINE VARIABLES                                              |
!==============================================================================|
				
   if(msr)then

!--Define NetCDF Output Filename 
   write(nchar,'(I4)')filecnt
   if(filecnt < 10)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_000'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 100)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_00'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 1000)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_0'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 10000)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_'//trim(adjustl(nchar))//'.nc'
   else
     write(*,*)'error in netcdf module'
     write(*,*)'# history files > 10000'
     stop
   endif

!--Create File 
   write(*,*)'cdfname = ',cdfname
   ierr = nf90_create(path=cdfname,cmode=nf90_clobber,ncid=nc_ofid)
   if(ierr /= nf90_eexist)then
     call handle_ncerr(ierr,"file creation error",ipt)
   else
     write(ipt,*)'file :',cdfname,' already exists'
     write(ipt,*)'exiting'
     stop
   end if

!--Description of File Contents
   ierr = nf90_put_att(nc_ofid,nf90_global,"title"      ,trim(casetitle))
   ierr = nf90_put_att(nc_ofid,nf90_global,"institution",trim(institution))
   ierr = nf90_put_att(nc_ofid,nf90_global,"source"     ,trim(fvcom_version))
   ierr = nf90_put_att(nc_ofid,nf90_global,"history"    ,trim(timestamp))
   ierr = nf90_put_att(nc_ofid,nf90_global,"references" ,trim(fvcom_website))
   ierr = nf90_put_att(nc_ofid,nf90_global,"Conventions",trim(netcdf_convention))
#  if defined (SPHERICAL)
   ierr = nf90_put_att(nc_ofid,nf90_global,"CoordinateSystem","GeoReferenced")
#  endif



!--Define Fixed Model Dimensions 
   ierr = nf90_def_dim(nc_ofid,"scalar" ,1      ,scl_did    )        
   ierr = nf90_def_dim(nc_ofid,"node"   ,mgl    ,node_did   )        
   ierr = nf90_def_dim(nc_ofid,"nele"   ,ngl    ,nele_did   )
   ierr = nf90_def_dim(nc_ofid,"siglay" ,kbm1   ,siglay_did )
   ierr = nf90_def_dim(nc_ofid,"siglev" ,kb     ,siglev_did )
   ierr = nf90_def_dim(nc_ofid,"three"  ,3      ,three_did  )
   ierr = nf90_def_dim(nc_ofid,"four"   ,4      ,four_did   )

!--Define Unlimited Model Dimension
   ierr = nf90_def_dim(nc_ofid,"time"   ,nf90_unlimited,time_did)

!--Set Up Data Dimensioning - Static Vars
   stat_scl     = (/scl_did/)             !!scalar variable               
   stat_lay     = (/siglay_did/)          !!vertical variables at layers
   stat_lev     = (/siglev_did/)          !!vertical variables at levels
   stat2de      = (/nele_did/)            !!2d element vars
   stat2dn      = (/node_did/)            !!2d nodal vars
   stat3de_lay  = (/nele_did,siglay_did/) !!3d element vars at layers
   stat3de_lev  = (/nele_did,siglev_did/) !!3d element vars at levels
   stat3dn_lay  = (/node_did,siglay_did/) !!3d node    vars at layers
   stat3dn_lev  = (/node_did,siglev_did/) !!3d node    vars at levels

!--Set Up Data Dimensioning - Dynamic Vars 
   dynm2de      = (/nele_did,time_did/)            !!2d element vars
   dynm2dn      = (/node_did,time_did/)            !!2d nodal vars
   dynm3de_lay  = (/nele_did,siglay_did,time_did/) !!3d elem vars at layers
   dynm3de_lev  = (/nele_did,siglev_did,time_did/) !!3d elem vars at levels
   dynm3dn_lay  = (/node_did,siglay_did,time_did/) !!3d node vars at layers
   dynm3dn_lev  = (/node_did,siglev_did,time_did/) !!3d node vars at levels
   dynmtime     = (/time_did/)   

!--Define Coordinate Variables and Attributes
   !!====NPROCS: Number of Processors=======================!
   ierr = nf90_def_var(nc_ofid,"nprocs",nf90_int,stat_scl,nprocs_vid)
   ierr = nf90_put_att(nc_ofid,nprocs_vid,"long_name","number of processors")

   !!====PARTITION: Partion Number of Element===============!
   ierr = nf90_def_var(nc_ofid,"partition",nf90_int,stat2de,partition_vid)
   ierr = nf90_put_att(nc_ofid,partition_vid,"long_name","partition")
   ierr = nf90_put_att(nc_ofid,partition_vid,"grid","TWOD_MESH")

   !!====Initial Density (Used for Constructing 3D Domain)==!
   ierr = nf90_def_var(nc_ofid,"Initial_Density",nf90_float,stat3dn_lay,idens_vid)
   ierr = nf90_put_att(nc_ofid,idens_vid,"long_name","Initial Density")
   ierr = nf90_put_att(nc_ofid,idens_vid,"grid","SigmaLayer_Mesh")

   !!====X Grid Coordinate at Nodes (VX) (Meters)===========!
   ierr = nf90_def_var(nc_ofid,"x",nf90_float,stat2dn,x_vid)
   ierr = nf90_put_att(nc_ofid,x_vid,"long_name","nodal x-coordinate")
   ierr = nf90_put_att(nc_ofid,x_vid,"units","meters")
   ierr = nf90_put_att(nc_ofid,x_vid,"grid","TWOD_MESH")

   !!====Y Grid Coordinate at Nodes (VY) (Meters)===========!
   ierr = nf90_def_var(nc_ofid,"y",nf90_float,stat2dn,y_vid)
   ierr = nf90_put_att(nc_ofid,y_vid,"long_name","nodal y-coordinate")
   ierr = nf90_put_att(nc_ofid,y_vid,"units","meters")
   ierr = nf90_put_att(nc_ofid,y_vid,"grid","TWOD_MESH")

   !!====Longitudinal Coordinate at Nodes (LON) (degrees)===!
   ierr = nf90_def_var(nc_ofid,"lon",nf90_float,stat2dn,lon_vid)
   ierr = nf90_put_att(nc_ofid,lon_vid,"long_name","Longitude")
   ierr = nf90_put_att(nc_ofid,lon_vid,"standard_name","longitude")
   ierr = nf90_put_att(nc_ofid,lon_vid,"units","degrees_east")
   ierr = nf90_put_att(nc_ofid,lon_vid,"grid","Bathymetry_Mesh")

   !!====Latitudinal  Coordinate at Nodes (LAT) (degrees)===!
   ierr = nf90_def_var(nc_ofid,"lat",nf90_float,stat2dn,lat_vid)
   ierr = nf90_put_att(nc_ofid,lat_vid,"long_name","Latitude")
   ierr = nf90_put_att(nc_ofid,lat_vid,"standard_name","latitude")
   ierr = nf90_put_att(nc_ofid,lat_vid,"units","degrees_north")
   ierr = nf90_put_att(nc_ofid,lat_vid,"grid","Bathymetry_Mesh")

   !!====Sigma Coordinate for Sigma Layers (ZZ)  (-)========!
   ierr = nf90_def_var(nc_ofid,"siglay",nf90_float,stat3dn_lay,siglay_vid)
   ierr = nf90_put_att(nc_ofid,siglay_vid,"long_name","Sigma Layers")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"standard_name","ocean_sigma/general_coordinate")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"formula_terms","siglay:siglay eta:zeta depth:depth")

   !!====Shifted Sigma Layer Coordinate for Viz ============!
   ierr = nf90_def_var(nc_ofid,"siglay_shift",nf90_float,stat3dn_lay,siglay_shift_vid)
   ierr = nf90_put_att(nc_ofid,siglay_shift_vid,"long_name","Shifted Sigma Layers")

   !!====Sigma Coordinate for Sigma Levels (Z)   (-)========!
   ierr = nf90_def_var(nc_ofid,"siglev",nf90_float,stat3dn_lev,siglev_vid)
   ierr = nf90_put_att(nc_ofid,siglev_vid,"long_name","Sigma Levels")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"standard_name","ocean_sigma/general_coordinate")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"formula_terms","siglev:siglev eta:zeta depth:depth")


!--Define Mesh Relevant Variables and Attributes

   !!====Bathymetry at Nodes (H) (meters)===================!
   ierr = nf90_def_var(nc_ofid,"h",nf90_float,stat2dn,h_vid)
   ierr = nf90_put_att(nc_ofid,h_vid,"long_name","Bathymetry")   
   ierr = nf90_put_att(nc_ofid,h_vid,"units","meters")
   ierr = nf90_put_att(nc_ofid,h_vid,"positive","down")
   ierr = nf90_put_att(nc_ofid,h_vid,"standard_name","depth")
   ierr = nf90_put_att(nc_ofid,h_vid,"grid","fvcom_grid")

   !!====Nodes surrounding each Element (NV)================!
   specdim = (/nele_did,three_did/) 
   ierr = nf90_def_var(nc_ofid,"nv",nf90_int,specdim,nv_vid)
   ierr = nf90_put_att(nc_ofid,nv_vid,"long_name","nodes surrounding element")     

   !!!====Momentum Stencil Interpolation Coefficients========!
   !specdim = (/nele_did,four_did/) 
   !ierr = nf90_def_var(nc_ofid,"a1u",nf90_float,specdim,a1u_vid)
   !ierr = nf90_put_att(nc_ofid,a1u_vid,"long_name","a1u")
   !ierr = nf90_def_var(nc_ofid,"a2u",nf90_float,specdim,a2u_vid)
   !ierr = nf90_put_att(nc_ofid,a2u_vid,"long_name","a2u")

   !!!====Element Based Interpolation Coefficients===========!
   !specdim = (/nele_did,three_did/) 
   !ierr = nf90_def_var(nc_ofid,"aw0",nf90_float,specdim,aw0_vid)
   !ierr = nf90_put_att(nc_ofid,aw0_vid,"long_name","aw0")
   !ierr = nf90_def_var(nc_ofid,"awx",nf90_float,specdim,awx_vid)
   !ierr = nf90_put_att(nc_ofid,awx_vid,"long_name","awx")
   !ierr = nf90_def_var(nc_ofid,"awy",nf90_float,specdim,awy_vid)
   !ierr = nf90_put_att(nc_ofid,awy_vid,"long_name","awy")

!--Define Model Time Variables and Attributes    
   ierr = nf90_def_var(nc_ofid,"time",nf90_float,dynmtime,time_vid)
   ierr = nf90_put_att(nc_ofid,time_vid,"long_name","Time")
   ierr = nf90_put_att(nc_ofid,time_vid,"units",trim(netcdf_timestring))
   ierr = nf90_put_att(nc_ofid,time_vid,"calendar","none")
   ierr = nf90_def_var(nc_ofid,"iint",nf90_int,dynmtime,iint_vid)
   ierr = nf90_put_att(nc_ofid,iint_vid,"long_name","internal mode iteration number")

!--Define Time Dependent Flow Variables (selected by user from input file)
   do i=1,nout_vars


   !write(*,*)'case is ',trim(cdf_vdp(i))
   
     select case(trim(cdf_vdp(i)))
	 
	 case("uu") !!===============U velocity======================================!
     ierr = nf90_def_var(nc_ofid,"uu",nf90_float,dynm3de_lay,uu_vid)
     ierr = nf90_put_att(nc_ofid,uu_vid,"long_name","U Velocity")
     ierr = nf90_put_att(nc_ofid,uu_vid,"standard_name","East-ward velocity")
     ierr = nf90_put_att(nc_ofid,uu_vid,"units","m s^-1")
     ierr = nf90_put_att(nc_ofid,uu_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,uu_vid,"type","data")

	 case("vv") !!===============V velocity======================================!
     ierr = nf90_def_var(nc_ofid,"vv",nf90_float,dynm3de_lay,vv_vid)
     ierr = nf90_put_att(nc_ofid,vv_vid,"long_name","V Velocity")
     ierr = nf90_put_att(nc_ofid,vv_vid,"standard_name","North-ward velocity")
     ierr = nf90_put_att(nc_ofid,vv_vid,"units","m s^-1")
     ierr = nf90_put_att(nc_ofid,vv_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,vv_vid,"type","data")	 
	 

	 case("ww") !!===============W velocity at node======================================!
     ierr = nf90_def_var(nc_ofid,"ww",nf90_float,dynm3de_lay,ww_vid)
     ierr = nf90_put_att(nc_ofid,ww_vid,"long_name","Vert. Velocity")
     ierr = nf90_put_att(nc_ofid,ww_vid,"standard_name","Vertical velocity")
     ierr = nf90_put_att(nc_ofid,ww_vid,"units","m s^-1")
     ierr = nf90_put_att(nc_ofid,ww_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,ww_vid,"type","data")	 	 
	 
	 
     case("t1") !!===============t1======================================!
     ierr = nf90_def_var(nc_ofid,"temp",nf90_float,dynm3dn_lay,t1_vid)
     ierr = nf90_put_att(nc_ofid,t1_vid,"long_name","temperature")
     ierr = nf90_put_att(nc_ofid,t1_vid,"standard_name","sea_water_temperature")
     ierr = nf90_put_att(nc_ofid,t1_vid,"units","degrees_C")
     ierr = nf90_put_att(nc_ofid,t1_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,t1_vid,"type","data")

     case("s1") !!===============s1======================================!
	 write(*,*)'Case = ',trim(cdf_vdp(i))
     ierr = nf90_def_var(nc_ofid,"salinity",nf90_float,dynm3dn_lay,s1_vid)
     ierr = nf90_put_att(nc_ofid,s1_vid,"long_name","salinity")
     ierr = nf90_put_att(nc_ofid,s1_vid,"standard_name","sea_water_salinity")
     ierr = nf90_put_att(nc_ofid,s1_vid,"units","1e-3")
     ierr = nf90_put_att(nc_ofid,s1_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,s1_vid,"type","data")

     case("el") !!===============el======================================!
     ierr = nf90_def_var(nc_ofid,"zeta",nf90_float,dynm2dn,el_vid)
     ierr = nf90_put_att(nc_ofid,el_vid,"long_name","Water Surface Elevation")
     ierr = nf90_put_att(nc_ofid,el_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,el_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,el_vid,"standard_name","sea_surface_elevation")
     ierr = nf90_put_att(nc_ofid,el_vid,"type","data")

     case("d") !!===============d=======================================!
     ierr = nf90_def_var(nc_ofid,"depth",nf90_float,dynm2dn,d_vid)
     ierr = nf90_put_att(nc_ofid,d_vid,"long_name","Water Depth")
     ierr = nf90_put_att(nc_ofid,d_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,d_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,d_vid,"type","data")
	
	 case("dz") !!===============dz======================================!
     ierr = nf90_def_var(nc_ofid,"dz",nf90_float,dynm3dn_lay,dz_vid)
     ierr = nf90_put_att(nc_ofid,dz_vid,"long_name","Delta Sigma Value")
     ierr = nf90_put_att(nc_ofid,dz_vid,"standard_name","DZ")
     ierr = nf90_put_att(nc_ofid,dz_vid,"units","dimensionless")
     ierr = nf90_put_att(nc_ofid,dz_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,dz_vid,"type","data")						 					 

	case("WC_CDOC1") !!===============WC_CDOC1====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDOC1",nf90_float,dynm3dn_lay,WC_CDOC1_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDOC1_vid,"long_name","Water Column Colored DOC1")
    ierr = nf90_put_att(nc_ofid,WC_CDOC1_vid,"standard_name","WC_CDOC1")
    ierr = nf90_put_att(nc_ofid,WC_CDOC1_vid,"units","mgC/L")
    ierr = nf90_put_att(nc_ofid,WC_CDOC1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDOC1_vid,"type","data")

	case("WC_CDOC2") !!===============WC_CDOC2====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDOC2",nf90_float,dynm3dn_lay,WC_CDOC2_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDOC2_vid,"long_name","Water Column Colored DOC2")
    ierr = nf90_put_att(nc_ofid,WC_CDOC2_vid,"standard_name","WC_CDOC2")
    ierr = nf90_put_att(nc_ofid,WC_CDOC2_vid,"units","mgC/L")
    ierr = nf90_put_att(nc_ofid,WC_CDOC2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDOC2_vid,"type","data")

	case("WC_CDOC3") !!===============WC_CDOC3====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDOC3",nf90_float,dynm3dn_lay,WC_CDOC3_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDOC3_vid,"long_name","Water Column Colored DOC3")
    ierr = nf90_put_att(nc_ofid,WC_CDOC3_vid,"standard_name","WC_CDOC3")
    ierr = nf90_put_att(nc_ofid,WC_CDOC3_vid,"units","mgC/L")
    ierr = nf90_put_att(nc_ofid,WC_CDOC3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDOC3_vid,"type","data")
	
	case("WC_NCDOC1") !!===============WC_NCDOC1====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDOC1",nf90_float,dynm3dn_lay,WC_NCDOC1_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDOC1_vid,"long_name","Water Column Non-Colored DOC1")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC1_vid,"standard_name","WC_NCDOC1")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC1_vid,"units","mgC/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC1_vid,"type","data")

	case("WC_NCDOC2") !!===============WC_NCDOC2====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDOC2",nf90_float,dynm3dn_lay,WC_NCDOC2_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDOC2_vid,"long_name","Water Column Non-Colored DOC2")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC2_vid,"standard_name","WC_NCDOC2")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC2_vid,"units","mgC/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC2_vid,"type","data")

	case("WC_NCDOC3") !!===============WC_NCDOC3====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDOC3",nf90_float,dynm3dn_lay,WC_NCDOC3_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDOC3_vid,"long_name","Water Column Non-Colored DOC3")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC3_vid,"standard_name","WC_NCDOC3")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC3_vid,"units","mgC/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDOC3_vid,"type","data")

  case("XFLUX_CDOC1") !!===============XFLUX_CDOC1====================================!
    ierr = nf90_def_var(nc_ofid,"XFLUX_CDOC1",nf90_float,dynm3dn_lay,XFLUX_CDOC1_vid)
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC1_vid,"long_name","Colored DOC1 Volume Flux")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC1_vid,"standard_name","XFLUX_CDOC1")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC1_vid,"units","gC/s")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC1_vid,"type","data")

	case("XFLUX_CDOC2") !!===============XFLUX_CDOC2====================================!
    ierr = nf90_def_var(nc_ofid,"XFLUX_CDOC2",nf90_float,dynm3dn_lay,XFLUX_CDOC2_vid)
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC2_vid,"long_name","Colored DOC2 Volume Flux")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC2_vid,"standard_name","XFLUX_CDOC2")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC2_vid,"units","gC/s")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC2_vid,"type","data")

	case("XFLUX_CDOC3") !!===============XFLUX_CDOC3====================================!
    ierr = nf90_def_var(nc_ofid,"XFLUX_CDOC3",nf90_float,dynm3dn_lay,XFLUX_CDOC3_vid)
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC3_vid,"long_name","Colored DOC3 Volume Flux")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC3_vid,"standard_name","XFLUX_CDOC3")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC3_vid,"units","gC/s")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,XFLUX_CDOC3_vid,"type","data")
	
  case("XFLUX_NCDOC1") !!===============XFLUX_NCDOC1====================================!
    ierr = nf90_def_var(nc_ofid,"XFLUX_NCDOC1",nf90_float,dynm3dn_lay,XFLUX_NCDOC1_vid)
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC1_vid,"long_name","Non-Colored DOC1 Volume Flux")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC1_vid,"standard_name","XFLUX_NCDOC1")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC1_vid,"units","gC/s")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC1_vid,"type","data")

  case("XFLUX_NCDOC2") !!===============XFLUX_NCDOC2====================================!
    ierr = nf90_def_var(nc_ofid,"XFLUX_NCDOC2",nf90_float,dynm3dn_lay,XFLUX_NCDOC2_vid)
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC2_vid,"long_name","Non-Colored DOC2 Volume Flux")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC2_vid,"standard_name","XFLUX_NCDOC2")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC2_vid,"units","gC/s")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC2_vid,"type","data")


  case("XFLUX_NCDOC3") !!===============XFLUX_NCDOC3====================================!
    ierr = nf90_def_var(nc_ofid,"XFLUX_NCDOC3",nf90_float,dynm3dn_lay,XFLUX_NCDOC3_vid)
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC3_vid,"long_name","Non-Colored DOC3 Volume Flux")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC3_vid,"standard_name","XFLUX_NCDOC3")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC3_vid,"units","gC/s")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,XFLUX_NCDOC3_vid,"type","data")
	
	!DON
	case("WC_CDON1") !!===============WC_CDON1====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDON1",nf90_float,dynm3dn_lay,WC_CDON1_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDON1_vid,"long_name","Water Column Colored DON1")
    ierr = nf90_put_att(nc_ofid,WC_CDON1_vid,"standard_name","WC_CDON1")
    ierr = nf90_put_att(nc_ofid,WC_CDON1_vid,"units","mgN/L")
    ierr = nf90_put_att(nc_ofid,WC_CDON1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDON1_vid,"type","data")

	case("WC_CDON2") !!===============WC_CDON2====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDON2",nf90_float,dynm3dn_lay,WC_CDON2_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDON2_vid,"long_name","Water Column Colored DON2")
    ierr = nf90_put_att(nc_ofid,WC_CDON2_vid,"standard_name","WC_CDON2")
    ierr = nf90_put_att(nc_ofid,WC_CDON2_vid,"units","mgN/L")
    ierr = nf90_put_att(nc_ofid,WC_CDON2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDON2_vid,"type","data")

	case("WC_CDON3") !!===============WC_CDON3====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDON3",nf90_float,dynm3dn_lay,WC_CDON3_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDON3_vid,"long_name","Water Column Colored DON3")
    ierr = nf90_put_att(nc_ofid,WC_CDON3_vid,"standard_name","WC_CDON3")
    ierr = nf90_put_att(nc_ofid,WC_CDON3_vid,"units","mgN/L")
    ierr = nf90_put_att(nc_ofid,WC_CDON3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDON3_vid,"type","data")
	
	case("WC_NCDON1") !!===============WC_NCDON1====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDON1",nf90_float,dynm3dn_lay,WC_NCDON1_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDON1_vid,"long_name","Water Column Non-Colored DON1")
    ierr = nf90_put_att(nc_ofid,WC_NCDON1_vid,"standard_name","WC_NCDON1")
    ierr = nf90_put_att(nc_ofid,WC_NCDON1_vid,"units","mgN/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDON1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDON1_vid,"type","data")

	case("WC_NCDON2") !!===============WC_NCDON2====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDON2",nf90_float,dynm3dn_lay,WC_NCDON2_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDON2_vid,"long_name","Water Column Non-Colored DON2")
    ierr = nf90_put_att(nc_ofid,WC_NCDON2_vid,"standard_name","WC_NCDON2")
    ierr = nf90_put_att(nc_ofid,WC_NCDON2_vid,"units","mgN/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDON2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDON2_vid,"type","data")

	case("WC_NCDON3") !!===============WC_NCDON3====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDON3",nf90_float,dynm3dn_lay,WC_NCDON3_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDON3_vid,"long_name","Water Column Non-Colored DON3")
    ierr = nf90_put_att(nc_ofid,WC_NCDON3_vid,"standard_name","WC_NCDON3")
    ierr = nf90_put_att(nc_ofid,WC_NCDON3_vid,"units","mgN/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDON3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDON3_vid,"type","data")
	
	case("WC_CDOP1") !!===============WC_CDOP1====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDOP1",nf90_float,dynm3dn_lay,WC_CDOP1_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDOP1_vid,"long_name","Water Column Colored DOP1")
    ierr = nf90_put_att(nc_ofid,WC_CDOP1_vid,"standard_name","WC_CDOP1")
    ierr = nf90_put_att(nc_ofid,WC_CDOP1_vid,"units","mgP/L")
    ierr = nf90_put_att(nc_ofid,WC_CDOP1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDOP1_vid,"type","data")

	case("WC_CDOP2") !!===============WC_CDOP2====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDOP2",nf90_float,dynm3dn_lay,WC_CDOP2_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDOP2_vid,"long_name","Water Column Colored DOP2")
    ierr = nf90_put_att(nc_ofid,WC_CDOP2_vid,"standard_name","WC_CDOP2")
    ierr = nf90_put_att(nc_ofid,WC_CDOP2_vid,"units","mgP/L")
    ierr = nf90_put_att(nc_ofid,WC_CDOP2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDOP2_vid,"type","data")

	case("WC_CDOP3") !!===============WC_CDOP3====================================!
    ierr = nf90_def_var(nc_ofid,"WC_CDOP3",nf90_float,dynm3dn_lay,WC_CDOP3_vid)
    ierr = nf90_put_att(nc_ofid,WC_CDOP3_vid,"long_name","Water Column Colored DOP3")
    ierr = nf90_put_att(nc_ofid,WC_CDOP3_vid,"standard_name","WC_CDOP3")
    ierr = nf90_put_att(nc_ofid,WC_CDOP3_vid,"units","mgP/L")
    ierr = nf90_put_att(nc_ofid,WC_CDOP3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_CDOP3_vid,"type","data")
	
	case("WC_NCDOP1") !!===============WC_NCDOP1====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDOP1",nf90_float,dynm3dn_lay,WC_NCDOP1_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDOP1_vid,"long_name","Water Column Non-Colored DOP1")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP1_vid,"standard_name","WC_NCDOP1")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP1_vid,"units","mgP/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP1_vid,"type","data")

	case("WC_NCDOP2") !!===============WC_NCDOP2====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDOP2",nf90_float,dynm3dn_lay,WC_NCDOP2_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDOP2_vid,"long_name","Water Column Non-Colored DOP2")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP2_vid,"standard_name","WC_NCDOP2")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP2_vid,"units","mgP/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP2_vid,"type","data")

	case("WC_NCDOP3") !!===============WC_NCDOP3====================================!
    ierr = nf90_def_var(nc_ofid,"WC_NCDOP3",nf90_float,dynm3dn_lay,WC_NCDOP3_vid)
    ierr = nf90_put_att(nc_ofid,WC_NCDOP3_vid,"long_name","Water Column Non-Colored DOP3")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP3_vid,"standard_name","WC_NCDOP3")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP3_vid,"units","mgP/L")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,WC_NCDOP3_vid,"type","data")
	
!# 	if defined (SEDDIAG) 
		
	 case("JCIN") !!===============jcin======================================!
     ierr = nf90_def_var(nc_ofid,"JCIN",nf90_float,dynm2dn,jcin_vid)
     ierr = nf90_put_att(nc_ofid,jcin_vid,"long_name","flux of particulate organic carbon form POC")
     ierr = nf90_put_att(nc_ofid,jcin_vid,"units","gC/m^2/d")
     ierr = nf90_put_att(nc_ofid,jcin_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,jcin_vid,"standard_name","POC flux")
     ierr = nf90_put_att(nc_ofid,jcin_vid,"type","data")

	 case("JNIN") !!===============jnin======================================!
     ierr = nf90_def_var(nc_ofid,"JNIN",nf90_float,dynm2dn,jnin_vid)
     ierr = nf90_put_att(nc_ofid,jnin_vid,"long_name","flux of particulate organic nitrogen form PON")
     ierr = nf90_put_att(nc_ofid,jnin_vid,"units","gN/m^2/d")
     ierr = nf90_put_att(nc_ofid,jnin_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,jnin_vid,"standard_name","PON flux")
     ierr = nf90_put_att(nc_ofid,jnin_vid,"type","data")
	 
	 case("JPIN") !!===============jpin======================================!
     ierr = nf90_def_var(nc_ofid,"JPIN",nf90_float,dynm2dn,jpin_vid)
     ierr = nf90_put_att(nc_ofid,jpin_vid,"long_name","flux of particulate organic phosphorus form POP")
     ierr = nf90_put_att(nc_ofid,jpin_vid,"units","gP/m^2/d")
     ierr = nf90_put_att(nc_ofid,jpin_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,jpin_vid,"standard_name","POP flux")
     ierr = nf90_put_att(nc_ofid,jpin_vid,"type","data")
	 

	 case("SOD") !!===============SOD==============================================!
     ierr = nf90_def_var(nc_ofid,"SOD",nf90_float,dynm2dn,SOD_vid)
     ierr = nf90_put_att(nc_ofid,SOD_vid,"long_name","sediment demand oxygen")
     ierr = nf90_put_att(nc_ofid,SOD_vid,"units","gO2/m^2/d")
     ierr = nf90_put_att(nc_ofid,SOD_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,SOD_vid,"standard_name","sediment demand oxygen")
     ierr = nf90_put_att(nc_ofid,SOD_vid,"type","data")

	 case("JNH4") !!===============JNH4============================================!
     ierr = nf90_def_var(nc_ofid,"JNH4",nf90_float,dynm2dn,JNH4_vid)
     ierr = nf90_put_att(nc_ofid,JNH4_vid,"long_name","flux of dissolved ammonia back to the overlying water column ")
     ierr = nf90_put_att(nc_ofid,JNH4_vid,"units","gN/m^2/d")
     ierr = nf90_put_att(nc_ofid,JNH4_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,JNH4_vid,"standard_name","dissolved NH4 flux")
     ierr = nf90_put_att(nc_ofid,JNH4_vid,"type","data")
	 
	 case("JNO3") !!===============JNO3============================================!
     ierr = nf90_def_var(nc_ofid,"JNO3",nf90_float,dynm2dn,JNO3_vid)
     ierr = nf90_put_att(nc_ofid,JNO3_vid,"long_name","flux of dissolved nitrate back to the overlying water column ")
     ierr = nf90_put_att(nc_ofid,JNO3_vid,"units","gN/m^2/d")
     ierr = nf90_put_att(nc_ofid,JNO3_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,JNO3_vid,"standard_name","dissolved NO3 flux")
     ierr = nf90_put_att(nc_ofid,JNO3_vid,"type","data")
	 
	 case("JDenitT") !!===============JDenitT======================================!
     ierr = nf90_def_var(nc_ofid,"JDenitT",nf90_float,dynm2dn,JDenitT_vid)
     ierr = nf90_put_att(nc_ofid,JDenitT_vid,"long_name","flux of denitrification back to the overlying water column")
     ierr = nf90_put_att(nc_ofid,JDenitT_vid,"units","gN/m^2/d")
     ierr = nf90_put_att(nc_ofid,JDenitT_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,JDenitT_vid,"standard_name","dissolved denitrification flux")
     ierr = nf90_put_att(nc_ofid,JDenitT_vid,"type","data")

	 case("JCH4") !!===============JCH4===========================================!
     ierr = nf90_def_var(nc_ofid,"JCH4",nf90_float,dynm2dn,JCH4_vid)
     ierr = nf90_put_att(nc_ofid,JCH4_vid,"long_name","flux of dissolved methane back to the overlying water column")
     ierr = nf90_put_att(nc_ofid,JCH4_vid,"units","gO2/m^2/d")
     ierr = nf90_put_att(nc_ofid,JCH4_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,JCH4_vid,"standard_name","dissolved methane flux")
     ierr = nf90_put_att(nc_ofid,JCH4_vid,"type","data")
	 
	 case("JCH4G") !!===============JCH4G=========================================!
     ierr = nf90_def_var(nc_ofid,"JCH4G",nf90_float,dynm2dn,JCH4G_vid)
     ierr = nf90_put_att(nc_ofid,JCH4G_vid,"long_name","flux of gaseous methane back to the overlying water column")
     ierr = nf90_put_att(nc_ofid,JCH4G_vid,"units","gO2/m^2/d")
     ierr = nf90_put_att(nc_ofid,JCH4G_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,JCH4G_vid,"standard_name","gaseous methane flux")
     ierr = nf90_put_att(nc_ofid,JCH4G_vid,"type","data")
	 
	 case("JHS") !!===============JHS=============================================!
     ierr = nf90_def_var(nc_ofid,"JHS",nf90_float,dynm2dn,JHS_vid)
     ierr = nf90_put_att(nc_ofid,JHS_vid,"long_name","flux of dissolved HS back to the overlying water column")
     ierr = nf90_put_att(nc_ofid,JHS_vid,"units","gO2/m^2/d")
     ierr = nf90_put_att(nc_ofid,JHS_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,JHS_vid,"standard_name","dissolved HS flux")
     ierr = nf90_put_att(nc_ofid,JHS_vid,"type","data")

	 case("JPO4") !!===============JPO4===========================================!
     ierr = nf90_def_var(nc_ofid,"JPO4",nf90_float,dynm2dn,JPO4_vid)
     ierr = nf90_put_att(nc_ofid,JPO4_vid,"long_name","flux of dissolved PO4 back to the overlying water column")
     ierr = nf90_put_att(nc_ofid,JPO4_vid,"units","gP/m^2/d")
     ierr = nf90_put_att(nc_ofid,JPO4_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,JPO4_vid,"standard_name","dissolved PO4 flux")
     ierr = nf90_put_att(nc_ofid,JPO4_vid,"type","data")


	 case("NH41") !!===============NH41===========================================!
     ierr = nf90_def_var(nc_ofid,"NH41",nf90_float,dynm2dn,NH41_vid)
     ierr = nf90_put_att(nc_ofid,NH41_vid,"long_name","dissolved ammonia in the 1st layer of sediment")
     ierr = nf90_put_att(nc_ofid,NH41_vid,"units","gN/m^3")
     ierr = nf90_put_att(nc_ofid,NH41_vid,"standard_name","1st layer sediment NH4")
     ierr = nf90_put_att(nc_ofid,NH41_vid,"type","data")
	 
	 case("NH42") !!===============NH42===========================================!
     ierr = nf90_def_var(nc_ofid,"NH42",nf90_float,dynm2dn,NH42_vid)
     ierr = nf90_put_att(nc_ofid,NH42_vid,"long_name","dissolved ammonia in the 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,NH42_vid,"units","gN/m^3")
     ierr = nf90_put_att(nc_ofid,NH42_vid,"standard_name","2nd layer sediment NH4")
     ierr = nf90_put_att(nc_ofid,NH42_vid,"type","data")

	 case("NO31") !!===============NO31===========================================!
     ierr = nf90_def_var(nc_ofid,"NO31",nf90_float,dynm2dn,NO31_vid)
     ierr = nf90_put_att(nc_ofid,NO31_vid,"long_name","dissolved nitrate in the 1st layer of sediment")
     ierr = nf90_put_att(nc_ofid,NO31_vid,"units","gN/m^3")
     ierr = nf90_put_att(nc_ofid,NO31_vid,"standard_name","1st layer sediment NO3")
     ierr = nf90_put_att(nc_ofid,NO31_vid,"type","data")
	 
	 case("NO32") !!===============NO32===========================================!
     ierr = nf90_def_var(nc_ofid,"NO32",nf90_float,dynm2dn,NO32_vid)
     ierr = nf90_put_att(nc_ofid,NO32_vid,"long_name","dissolved nitrate in the 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,NO32_vid,"units","gN/m^3")
     ierr = nf90_put_att(nc_ofid,NO32_vid,"standard_name","2nd layer sediment NO3")
     ierr = nf90_put_att(nc_ofid,NO32_vid,"type","data")
	 
 	 case("PO41") !!===============PO41===========================================!
     ierr = nf90_def_var(nc_ofid,"PO41",nf90_float,dynm2dn,PO41_vid)
     ierr = nf90_put_att(nc_ofid,PO41_vid,"long_name","dissolved PO4 in the 1st layer of sediment")
     ierr = nf90_put_att(nc_ofid,PO41_vid,"units","gP/m^3")
     ierr = nf90_put_att(nc_ofid,PO41_vid,"standard_name","1st layer sediment PO4")
     ierr = nf90_put_att(nc_ofid,PO41_vid,"type","data") 
	 
	 case("PO42") !!===============PO42===========================================!
     ierr = nf90_def_var(nc_ofid,"PO42",nf90_float,dynm2dn,PO42_vid)
     ierr = nf90_put_att(nc_ofid,PO42_vid,"long_name","dissolved PO4 in the 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,PO42_vid,"units","gP/m^3")
     ierr = nf90_put_att(nc_ofid,PO42_vid,"standard_name","2nd layer sediment PO4")
     ierr = nf90_put_att(nc_ofid,PO42_vid,"type","data") 	 
	 
	 case("CH41") !!===============CH41===========================================!
     ierr = nf90_def_var(nc_ofid,"CH41",nf90_float,dynm2dn,CH41_vid)
     ierr = nf90_put_att(nc_ofid,CH41_vid,"long_name","dissolved methane in the 1st layer of sediment")
     ierr = nf90_put_att(nc_ofid,CH41_vid,"units","mgO2-equivalents/L")
     ierr = nf90_put_att(nc_ofid,CH41_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,CH41_vid,"standard_name","1st layer sediment CH4")
     ierr = nf90_put_att(nc_ofid,CH41_vid,"type","data")
	 
	 case("CH42") !!===============CH42===========================================!
     ierr = nf90_def_var(nc_ofid,"CH42",nf90_float,dynm2dn,CH42_vid)
     ierr = nf90_put_att(nc_ofid,CH42_vid,"long_name","dissolved methane in the 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,CH42_vid,"units","mgO2-equivalents/L")
     ierr = nf90_put_att(nc_ofid,CH42_vid,"standard_name","2nd layer sediment CH4")
     ierr = nf90_put_att(nc_ofid,CH42_vid,"type","data")
	 
	 case("HS1") !!===============HS1=============================================!
     ierr = nf90_def_var(nc_ofid,"HS1",nf90_float,dynm2dn,HS1_vid)
     ierr = nf90_put_att(nc_ofid,HS1_vid,"long_name","dissolved bisulfide in the 1st layer of sediment")
     ierr = nf90_put_att(nc_ofid,HS1_vid,"units","mgO2-equivalents/L")
     ierr = nf90_put_att(nc_ofid,HS1_vid,"standard_name","1st layer sediment HS")
     ierr = nf90_put_att(nc_ofid,HS1_vid,"type","data")
	 
	 case("HS2") !!===============HS2============================================!
     ierr = nf90_def_var(nc_ofid,"HS2",nf90_float,dynm2dn,HS2_vid)
     ierr = nf90_put_att(nc_ofid,HS2_vid,"long_name","dissolved bisulfide in the 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,HS2_vid,"units","mgO2-equivalents/L")
     ierr = nf90_put_att(nc_ofid,HS2_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,HS2_vid,"standard_name","2nd layer sediment HS")
     ierr = nf90_put_att(nc_ofid,HS2_vid,"type","data")
	 
	 case("POC2G1") !!===============POC2G1======================================!
     ierr = nf90_def_var(nc_ofid,"POC2G1",nf90_float,dynm2dn,POC2G1_vid)
     ierr = nf90_put_att(nc_ofid,POC2G1_vid,"long_name","particular organic carbon (POC) for G1 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,POC2G1_vid,"units","gO2/m^3")
     ierr = nf90_put_att(nc_ofid,POC2G1_vid,"standard_name","G1's POC in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,POC2G1_vid,"type","data")

	 case("POC2G2") !!===============POC2G2======================================!
     ierr = nf90_def_var(nc_ofid,"POC2G2",nf90_float,dynm2dn,POC2G2_vid)
     ierr = nf90_put_att(nc_ofid,POC2G2_vid,"long_name","particular organic carbon (POC) for G2 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,POC2G2_vid,"units","gO2/m^3")
     ierr = nf90_put_att(nc_ofid,POC2G2_vid,"standard_name","G2's POC in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,POC2G2_vid,"type","data")
	 
	 case("POC2G3") !!===============POC2G3======================================!
     ierr = nf90_def_var(nc_ofid,"POC2G3",nf90_float,dynm2dn,POC2G3_vid)
     ierr = nf90_put_att(nc_ofid,POC2G3_vid,"long_name","particular organic carbon (POC) for G3 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,POC2G3_vid,"units","gO2/m^3")
     ierr = nf90_put_att(nc_ofid,POC2G3_vid,"standard_name","G3's POC in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,POC2G3_vid,"type","data")
	 
	 case("PON2G1") !!===============PON2G1======================================!
     ierr = nf90_def_var(nc_ofid,"PON2G1",nf90_float,dynm2dn,PON2G1_vid)
     ierr = nf90_put_att(nc_ofid,PON2G1_vid,"long_name","particular organic nitrogen (PON) for G1 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,PON2G1_vid,"units","gN/m^3")
     ierr = nf90_put_att(nc_ofid,PON2G1_vid,"standard_name","G1's PON in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,PON2G1_vid,"type","data")

	 case("PON2G2") !!===============PON2G2======================================!
     ierr = nf90_def_var(nc_ofid,"PON2G2",nf90_float,dynm2dn,PON2G2_vid)
     ierr = nf90_put_att(nc_ofid,PON2G2_vid,"long_name","particular organic nitrogen (PON) for G2 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,PON2G2_vid,"units","gN/m^3")
     ierr = nf90_put_att(nc_ofid,PON2G2_vid,"standard_name","G2's PON in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,PON2G2_vid,"type","data")
	 
	 case("PON2G3") !!===============PON2G3======================================!
     ierr = nf90_def_var(nc_ofid,"PON2G3",nf90_float,dynm2dn,PON2G3_vid)
     ierr = nf90_put_att(nc_ofid,PON2G3_vid,"long_name","particular organic nitrogen (PON) for G3 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,PON2G3_vid,"units","gN/m^3")
     ierr = nf90_put_att(nc_ofid,PON2G3_vid,"standard_name","G3's PON in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,PON2G3_vid,"type","data")
	 
	 case("POP2G1") !!===============POP2G1======================================!
     ierr = nf90_def_var(nc_ofid,"POP2G1",nf90_float,dynm2dn,POP2G1_vid)
     ierr = nf90_put_att(nc_ofid,POP2G1_vid,"long_name","particular organic phosphorus (POP) for G1 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,POP2G1_vid,"units","gP/m^3")
     ierr = nf90_put_att(nc_ofid,POP2G1_vid,"standard_name","G1's POP in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,POP2G1_vid,"type","data")

	 case("POP2G2") !!===============POP2G2======================================!
     ierr = nf90_def_var(nc_ofid,"POP2G2",nf90_float,dynm2dn,POP2G2_vid)
     ierr = nf90_put_att(nc_ofid,POP2G2_vid,"long_name","particular organic phosphorus (POP) for G2 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,POP2G2_vid,"units","gP/m^3")
     ierr = nf90_put_att(nc_ofid,POP2G2_vid,"standard_name","G2's POP in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,POP2G2_vid,"type","data")
	 
	 case("POP2G3") !!===============POP2G3======================================!
     ierr = nf90_def_var(nc_ofid,"POP2G3",nf90_float,dynm2dn,POP2G3_vid)
     ierr = nf90_put_att(nc_ofid,POP2G3_vid,"long_name","particular organic phosphorus (POP)) for G3 class in 2nd layer of sediment")
     ierr = nf90_put_att(nc_ofid,POP2G3_vid,"units","gP/m^3")
     ierr = nf90_put_att(nc_ofid,POP2G3_vid,"standard_name","G3's POP in 2nd layer sediment")
     ierr = nf90_put_att(nc_ofid,POP2G3_vid,"type","data")

	 case("H1") !!===============H1======================================!
     ierr = nf90_def_var(nc_ofid,"H1",nf90_float,dynm2dn,H1_vid)
     ierr = nf90_put_att(nc_ofid,H1_vid,"long_name","sediment layer thickness")
     ierr = nf90_put_att(nc_ofid,H1_vid,"units","meter")
     ierr = nf90_put_att(nc_ofid,H1_vid,"standard_name","sediment layer thickness")
     ierr = nf90_put_att(nc_ofid,H1_vid,"type","data")
	 
	 case("BENSTR") !!===============BENSTR======================================!
     ierr = nf90_def_var(nc_ofid,"BENSTR",nf90_float,dynm2dn,BENSTR_vid)
     ierr = nf90_put_att(nc_ofid,BENSTR_vid,"long_name","benthic stress")
     ierr = nf90_put_att(nc_ofid,BENSTR_vid,"units","dimensionless")
     ierr = nf90_put_att(nc_ofid,BENSTR_vid,"standard_name","bethic stress")
     ierr = nf90_put_att(nc_ofid,BENSTR_vid,"type","data")
 
!#    if defined (SED_DOM)
	 
	 case("SED_CDOC1") !!===============SED_CDOC1======================================!
     ierr = nf90_def_var(nc_ofid,"SED_CDOC1",nf90_float,dynm2dn,SED_CDOC1_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDOC1_vid,"long_name","sediment colored DOC1")
	 ierr = nf90_put_att(nc_ofid,SED_CDOC1_vid,"standard_name","SED_CDOC1")
     ierr = nf90_put_att(nc_ofid,SED_CDOC1_vid,"units","mgC/L")
     ierr = nf90_put_att(nc_ofid,SED_CDOC1_vid,"type","data")
	 
	 case("SED_CDOC2") !!===============SED_CDOC2======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_CDOC2",nf90_float,dynm2dn,SED_CDOC2_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDOC2_vid,"long_name","sediment colored DOC2")
	 ierr = nf90_put_att(nc_ofid,SED_CDOC2_vid,"standard_name","SED_CDOC2")
     ierr = nf90_put_att(nc_ofid,SED_CDOC2_vid,"units","mgC/L")
     ierr = nf90_put_att(nc_ofid,SED_CDOC2_vid,"type","data")

 	 case("SED_CDOC3") !!===============SED_CDOC3======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_CDOC3",nf90_float,dynm2dn,SED_CDOC3_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDOC3_vid,"long_name","sediment colored DOC3")
	 ierr = nf90_put_att(nc_ofid,SED_CDOC3_vid,"standard_name","SED_CDOC3")
     ierr = nf90_put_att(nc_ofid,SED_CDOC3_vid,"units","mgC/L")
     ierr = nf90_put_att(nc_ofid,SED_CDOC3_vid,"type","data")
	 
	 case("SED_NCDOC1") !!===============SED_NCDOC1======================================!
     ierr = nf90_def_var(nc_ofid,"SED_NCDOC1",nf90_float,dynm2dn,SED_NCDOC1_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDOC1_vid,"long_name","sediment non-colored DOC1")
	 ierr = nf90_put_att(nc_ofid,SED_NCDOC1_vid,"standard_name","SED_NCDOC1")
     ierr = nf90_put_att(nc_ofid,SED_NCDOC1_vid,"units","mgC/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDOC1_vid,"type","data")
	 
	 case("SED_NCDOC2") !!===============SED_NCDOC2======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_NCDOC2",nf90_float,dynm2dn,SED_NCDOC2_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDOC2_vid,"long_name","sediment non-colored DOC2")
	 ierr = nf90_put_att(nc_ofid,SED_NCDOC2_vid,"standard_name","SED_NCDOC2")
     ierr = nf90_put_att(nc_ofid,SED_NCDOC2_vid,"units","mgC/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDOC2_vid,"type","data")

 	 case("SED_NCDOC3") !!===============SED_NCDOC3======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_NCDOC3",nf90_float,dynm2dn,SED_NCDOC3_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDOC3_vid,"long_name","sediment non-colored DOC3")
	 ierr = nf90_put_att(nc_ofid,SED_NCDOC3_vid,"standard_name","SED_NCDOC3")
     ierr = nf90_put_att(nc_ofid,SED_NCDOC3_vid,"units","mgC/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDOC3_vid,"type","data")

	!DON
	 
	 case("SED_CDON1") !!===============SED_CDON1======================================!
     ierr = nf90_def_var(nc_ofid,"SED_CDON1",nf90_float,dynm2dn,SED_CDON1_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDON1_vid,"long_name","sediment colored DON1")
	 ierr = nf90_put_att(nc_ofid,SED_CDON1_vid,"standard_name","SED_CDON1")
     ierr = nf90_put_att(nc_ofid,SED_CDON1_vid,"units","mgN/L")
     ierr = nf90_put_att(nc_ofid,SED_CDON1_vid,"type","data")
	 
	 case("SED_CDON2") !!===============SED_CDON2======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_CDON2",nf90_float,dynm2dn,SED_CDON2_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDON2_vid,"long_name","sediment colored DON2")
	 ierr = nf90_put_att(nc_ofid,SED_CDON2_vid,"standard_name","SED_CDON2")
     ierr = nf90_put_att(nc_ofid,SED_CDON2_vid,"units","mgN/L")
     ierr = nf90_put_att(nc_ofid,SED_CDON2_vid,"type","data")

 	 case("SED_CDON3") !!===============SED_CDON3======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_CDON3",nf90_float,dynm2dn,SED_CDON3_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDON3_vid,"long_name","sediment colored DON3")
	 ierr = nf90_put_att(nc_ofid,SED_CDON3_vid,"standard_name","SED_CDON3")
     ierr = nf90_put_att(nc_ofid,SED_CDON3_vid,"units","mgN/L")
     ierr = nf90_put_att(nc_ofid,SED_CDON3_vid,"type","data")
	 

	 case("SED_NCDON1") !!===============SED_NCDON1======================================!
     ierr = nf90_def_var(nc_ofid,"SED_NCDON1",nf90_float,dynm2dn,SED_NCDON1_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDON1_vid,"long_name","sediment non-colored DON1")
	 ierr = nf90_put_att(nc_ofid,SED_NCDON1_vid,"standard_name","SED_NCDON1")
     ierr = nf90_put_att(nc_ofid,SED_NCDON1_vid,"units","mgN/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDON1_vid,"type","data")
	 
	 case("SED_NCDON2") !!===============SED_NCDON2======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_NCDON2",nf90_float,dynm2dn,SED_NCDON2_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDON2_vid,"long_name","sediment non-colored DON2")
	 ierr = nf90_put_att(nc_ofid,SED_NCDON2_vid,"standard_name","SED_NCDON2")
     ierr = nf90_put_att(nc_ofid,SED_NCDON2_vid,"units","mgN/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDON2_vid,"type","data")

 	 case("SED_NCDON3") !!===============SED_NCDON3======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_NCDON3",nf90_float,dynm2dn,SED_NCDON3_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDON3_vid,"long_name","sediment non-colored DON3")
	 ierr = nf90_put_att(nc_ofid,SED_NCDON3_vid,"standard_name","SED_NCDON3")
     ierr = nf90_put_att(nc_ofid,SED_NCDON3_vid,"units","mgN/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDON3_vid,"type","data")
	 
	 !DOP
	 case("SED_CDOP1") !!===============SED_CDOP1======================================!
     ierr = nf90_def_var(nc_ofid,"SED_CDOP1",nf90_float,dynm2dn,SED_CDOP1_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDOP1_vid,"long_name","sediment colored DOP1")
	 ierr = nf90_put_att(nc_ofid,SED_CDOP1_vid,"standard_name","SED_CDOP1")
     ierr = nf90_put_att(nc_ofid,SED_CDOP1_vid,"units","mgP/L")
     ierr = nf90_put_att(nc_ofid,SED_CDOP1_vid,"type","data")
	 
	 case("SED_CDOP2") !!===============SED_CDOP2======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_CDOP2",nf90_float,dynm2dn,SED_CDOP2_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDOP2_vid,"long_name","sediment colored DOP2")
	 ierr = nf90_put_att(nc_ofid,SED_CDOP2_vid,"standard_name","SED_CDOP2")
     ierr = nf90_put_att(nc_ofid,SED_CDOP2_vid,"units","mgP/L")
     ierr = nf90_put_att(nc_ofid,SED_CDOP2_vid,"type","data")

 	 case("SED_CDOP3") !!===============SED_CDOP3======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_CDOP3",nf90_float,dynm2dn,SED_CDOP3_vid)
     ierr = nf90_put_att(nc_ofid,SED_CDOP3_vid,"long_name","sediment colored DOP3")
	 ierr = nf90_put_att(nc_ofid,SED_CDOP3_vid,"standard_name","SED_CDOP3")
     ierr = nf90_put_att(nc_ofid,SED_CDOP3_vid,"units","mgP/L")
     ierr = nf90_put_att(nc_ofid,SED_CDOP3_vid,"type","data")
	 
	 case("SED_NCDOP1") !!===============SED_NCDOP1======================================!
     ierr = nf90_def_var(nc_ofid,"SED_NCDOP1",nf90_float,dynm2dn,SED_NCDOP1_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDOP1_vid,"long_name","sediment non-colored DOP1")
	 ierr = nf90_put_att(nc_ofid,SED_NCDOP1_vid,"standard_name","SED_NCDOP1")
     ierr = nf90_put_att(nc_ofid,SED_NCDOP1_vid,"units","mgP/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDOP1_vid,"type","data")
	 
	 case("SED_NCDOP2") !!===============SED_NCDOP2======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_NCDOP2",nf90_float,dynm2dn,SED_NCDOP2_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDOP2_vid,"long_name","sediment non-colored DOP2")
	 ierr = nf90_put_att(nc_ofid,SED_NCDOP2_vid,"standard_name","SED_NCDOP2")
     ierr = nf90_put_att(nc_ofid,SED_NCDOP2_vid,"units","mgP/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDOP2_vid,"type","data")

 	 case("SED_NCDOP3") !!===============SED_NCDOP3======================================!
	 ierr = nf90_def_var(nc_ofid,"SED_NCDOP3",nf90_float,dynm2dn,SED_NCDOP3_vid)
     ierr = nf90_put_att(nc_ofid,SED_NCDOP3_vid,"long_name","sediment non-colored DOP3")
	 ierr = nf90_put_att(nc_ofid,SED_NCDOP3_vid,"standard_name","SED_NCDOP3")
     ierr = nf90_put_att(nc_ofid,SED_NCDOP3_vid,"units","mgP/L")
     ierr = nf90_put_att(nc_ofid,SED_NCDOP3_vid,"type","data")
	 
!#    endif	 
	 
	 
!#	 if defined (SED_DOM) 

 	 case("JWCDOC1") !!===============JWCDOC1======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDOC1",nf90_float,dynm2dn,JWCDOC1_vid)
     ierr = nf90_put_att(nc_ofid,JWCDOC1_vid,"long_name","sediment-water flux of colored DOC1")
	 ierr = nf90_put_att(nc_ofid,JWCDOC1_vid,"standard_name","JWCDOC1")
     ierr = nf90_put_att(nc_ofid,JWCDOC1_vid,"units","gC/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDOC1_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDOC1_vid,"positive","up")

 	 case("JWCDOC2") !!===============JWCDOC2======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDOC2",nf90_float,dynm2dn,JWCDOC2_vid)
     ierr = nf90_put_att(nc_ofid,JWCDOC2_vid,"long_name","sediment-water flux of colored DOC2")
	 ierr = nf90_put_att(nc_ofid,JWCDOC2_vid,"standard_name","JWCDOC2")
     ierr = nf90_put_att(nc_ofid,JWCDOC2_vid,"units","gC/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDOC2_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDOC2_vid,"positive","up")
	 
	 
 	 case("JWCDOC3") !!===============JWCDOC3======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDOC3",nf90_float,dynm2dn,JWCDOC3_vid)
     ierr = nf90_put_att(nc_ofid,JWCDOC3_vid,"long_name","sediment-water flux of colored DOC3")
	 ierr = nf90_put_att(nc_ofid,JWCDOC3_vid,"standard_name","JWCDOC3")
     ierr = nf90_put_att(nc_ofid,JWCDOC3_vid,"units","gC/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDOC3_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDOC3_vid,"positive","up")	 
	 
	 case("JWNCDOC1") !!===============JWNCDOC1======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDOC1",nf90_float,dynm2dn,JWNCDOC1_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDOC1_vid,"long_name","sediment-water flux of non-colored DOC1")
	 ierr = nf90_put_att(nc_ofid,JWNCDOC1_vid,"standard_name","JWNCDOC1")
     ierr = nf90_put_att(nc_ofid,JWNCDOC1_vid,"units","gC/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDOC1_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDOC1_vid,"positive","up")


 	 case("JWNCDOC2") !!===============JWNCDOC2======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDOC2",nf90_float,dynm2dn,JWNCDOC2_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDOC2_vid,"long_name","sediment-water flux of non-colored DOC2")
	 ierr = nf90_put_att(nc_ofid,JWNCDOC2_vid,"standard_name","JWNCDOC2")
     ierr = nf90_put_att(nc_ofid,JWNCDOC2_vid,"units","gC/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDOC2_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDOC2_vid,"positive","up")	 
	 
 	 case("JWNCDOC3") !!===============JWNCDOC3======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDOC3",nf90_float,dynm2dn,JWNCDOC3_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDOC3_vid,"long_name","sediment-water flux of non-colored DOC3")
	 ierr = nf90_put_att(nc_ofid,JWNCDOC3_vid,"standard_name","JWNCDOC3")
     ierr = nf90_put_att(nc_ofid,JWNCDOC3_vid,"units","gC/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDOC3_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDOC3_vid,"positive","up")	 
	 
	 !DON
	 
	 
 	 case("JWCDON1") !!===============JWCDON1======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDON1",nf90_float,dynm2dn,JWCDON1_vid)
     ierr = nf90_put_att(nc_ofid,JWCDON1_vid,"long_name","sediment-water flux of colored DON1")
	 ierr = nf90_put_att(nc_ofid,JWCDON1_vid,"standard_name","JWCDON1")
     ierr = nf90_put_att(nc_ofid,JWCDON1_vid,"units","gN/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDON1_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDON1_vid,"positive","up")
	 
 	 case("JWCDON2") !!===============JWCDON2======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDON2",nf90_float,dynm2dn,JWCDON2_vid)
     ierr = nf90_put_att(nc_ofid,JWCDON2_vid,"long_name","sediment-water flux of colored DON2")
	 ierr = nf90_put_att(nc_ofid,JWCDON2_vid,"standard_name","JWCDON2")
     ierr = nf90_put_att(nc_ofid,JWCDON2_vid,"units","gN/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDON2_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDON2_vid,"positive","up")	 
	 
	 
 	 case("JWCDON3") !!===============JWCDON3======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDON3",nf90_float,dynm2dn,JWCDON3_vid)
     ierr = nf90_put_att(nc_ofid,JWCDON3_vid,"long_name","sediment-water flux of colored DON3")
	 ierr = nf90_put_att(nc_ofid,JWCDON3_vid,"standard_name","JWCDON3")
     ierr = nf90_put_att(nc_ofid,JWCDON3_vid,"units","gN/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDON3_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDON3_vid,"positive","up")	 
	 
	 case("JWNCDON1") !!===============JWNCDON1======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDON1",nf90_float,dynm2dn,JWNCDON1_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDON1_vid,"long_name","sediment-water flux of non-colored DON1")
	 ierr = nf90_put_att(nc_ofid,JWNCDON1_vid,"standard_name","JWNCDON1")
     ierr = nf90_put_att(nc_ofid,JWNCDON1_vid,"units","gN/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDON1_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDON1_vid,"positive","up")

 	 case("JWNCDON2") !!===============JWNCDON2======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDON2",nf90_float,dynm2dn,JWNCDON2_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDON2_vid,"long_name","sediment-water flux of non-colored DON2")
	 ierr = nf90_put_att(nc_ofid,JWNCDON2_vid,"standard_name","JWNCDON2")
     ierr = nf90_put_att(nc_ofid,JWNCDON2_vid,"units","gN/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDON2_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDON2_vid,"positive","up")	 
	 
	 
 	 case("JWNCDON3") !!===============JWNCDON3======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDON3",nf90_float,dynm2dn,JWNCDON3_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDON3_vid,"long_name","sediment-water flux of non-colored DON3")
	 ierr = nf90_put_att(nc_ofid,JWNCDON3_vid,"standard_name","JWNCDON3")
     ierr = nf90_put_att(nc_ofid,JWNCDON3_vid,"units","gN/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDON3_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDON3_vid,"positive","up")	 
	 
	  !DOP
 	 case("JWCDOP1") !!===============JWCDOP1======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDOP1",nf90_float,dynm2dn,JWCDOP1_vid)
     ierr = nf90_put_att(nc_ofid,JWCDOP1_vid,"long_name","sediment-water flux of colored DOP1")
	 ierr = nf90_put_att(nc_ofid,JWCDOP1_vid,"standard_name","JWCDOP1")
     ierr = nf90_put_att(nc_ofid,JWCDOP1_vid,"units","gP/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDOP1_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDOP1_vid,"positive","up")	 

 	 case("JWCDOP2") !!===============JWCDOP2======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDOP2",nf90_float,dynm2dn,JWCDOP2_vid)
     ierr = nf90_put_att(nc_ofid,JWCDOP2_vid,"long_name","sediment-water flux of colored DOP2")
	 ierr = nf90_put_att(nc_ofid,JWCDOP2_vid,"standard_name","JWCDOP2")
     ierr = nf90_put_att(nc_ofid,JWCDOP2_vid,"units","gP/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDOP2_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDOP2_vid,"positive","up")	 

 	 case("JWCDOP3") !!===============JWCDOP3======================================!
	 ierr = nf90_def_var(nc_ofid,"JWCDOP3",nf90_float,dynm2dn,JWCDOP3_vid)
     ierr = nf90_put_att(nc_ofid,JWCDOP3_vid,"long_name","sediment-water flux of colored DOP3")
	 ierr = nf90_put_att(nc_ofid,JWCDOP3_vid,"standard_name","JWCDOP3")
     ierr = nf90_put_att(nc_ofid,JWCDOP3_vid,"units","gP/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWCDOP3_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWCDOP3_vid,"positive","up")	 
	 
	 case("JWNCDOP1") !!===============JWNCDOP1======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDOP1",nf90_float,dynm2dn,JWNCDOP1_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDOP1_vid,"long_name","sediment-water flux of non-colored DOP1")
	 ierr = nf90_put_att(nc_ofid,JWNCDOP1_vid,"standard_name","JWNCDOP1")
     ierr = nf90_put_att(nc_ofid,JWNCDOP1_vid,"units","gP/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDOP1_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDOP1_vid,"positive","up")	 


 	 case("JWNCDOP2") !!===============JWNCDON2======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDOP2",nf90_float,dynm2dn,JWNCDOP2_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDOP2_vid,"long_name","sediment-water flux of non-colored DOP2")
	 ierr = nf90_put_att(nc_ofid,JWNCDOP2_vid,"standard_name","JWNCDON2")
     ierr = nf90_put_att(nc_ofid,JWNCDOP2_vid,"units","gN/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDOP2_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDOP2_vid,"positive","up")	 
	 
 	 case("JWNCDOP3") !!===============JWNCDOP3======================================!
	 ierr = nf90_def_var(nc_ofid,"JWNCDOP3",nf90_float,dynm2dn,JWNCDOP3_vid)
     ierr = nf90_put_att(nc_ofid,JWNCDOP3_vid,"long_name","sediment-water flux of non-colored DOP3")
	 ierr = nf90_put_att(nc_ofid,JWNCDOP3_vid,"standard_name","JWNCDOP3")
     ierr = nf90_put_att(nc_ofid,JWNCDOP3_vid,"units","gP/m^2/day")
     ierr = nf90_put_att(nc_ofid,JWNCDOP3_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,JWNCDOP3_vid,"positive","up")	 
	 
!#    endif
	 
!#	endif
	
	case("DOXG") !!===============DOXG====================================!
    ierr = nf90_def_var(nc_ofid,"DOXG",nf90_float,dynm3dn_lay,DOXG_vid)
    ierr = nf90_put_att(nc_ofid,DOXG_vid,"long_name","Dissolved O2 Concentration")
    ierr = nf90_put_att(nc_ofid,DOXG_vid,"standard_name","ALG 1")
    ierr = nf90_put_att(nc_ofid,DOXG_vid,"units","g O2 m^-3")
    ierr = nf90_put_att(nc_ofid,DOXG_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,DOXG_vid,"type","data")						 

	 
	 case("NPP") !!===============NPP======================================!
	 ierr = nf90_def_var(nc_ofid,"NPP",nf90_float,dynm2dn,NPP_vid)
     ierr = nf90_put_att(nc_ofid,NPP_vid,"long_name","Net Primary Production")
	 ierr = nf90_put_att(nc_ofid,NPP_vid,"standard_name","NPP")
     ierr = nf90_put_att(nc_ofid,NPP_vid,"units","gC m^-2 d^-1")
     ierr = nf90_put_att(nc_ofid,NPP_vid,"type","data")
	 ierr = nf90_put_att(nc_ofid,NPP_vid,"positive","up")	 					 
						 
						 
	case("B1") !!===============Algae 1====================================!
    ierr = nf90_def_var(nc_ofid,"B1",nf90_float,dynm3dn_lay,B1_vid)
    ierr = nf90_put_att(nc_ofid,B1_vid,"long_name","Algae 1 Concentration")
    ierr = nf90_put_att(nc_ofid,B1_vid,"standard_name","ALG 1")
    ierr = nf90_put_att(nc_ofid,B1_vid,"units","gC m^-3")
    ierr = nf90_put_att(nc_ofid,B1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,B1_vid,"type","data")						 
						 
						 
	case("B2") !!===============Algae 2====================================!
    ierr = nf90_def_var(nc_ofid,"B2",nf90_float,dynm3dn_lay,B2_vid)
    ierr = nf90_put_att(nc_ofid,B2_vid,"long_name","Algae 2 Concentration")
    ierr = nf90_put_att(nc_ofid,B2_vid,"standard_name","ALG 2")
    ierr = nf90_put_att(nc_ofid,B2_vid,"units","gC m^-3")
    ierr = nf90_put_att(nc_ofid,B2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,B2_vid,"type","data")								 
						 
	case("NH4") !!===============NH4====================================!
    ierr = nf90_def_var(nc_ofid,"NH4",nf90_float,dynm3dn_lay,NH4_vid)
    ierr = nf90_put_att(nc_ofid,NH4_vid,"long_name","Ammonium Concentration")
    ierr = nf90_put_att(nc_ofid,NH4_vid,"standard_name","NH4")
    ierr = nf90_put_att(nc_ofid,NH4_vid,"units","gN m^-3")
    ierr = nf90_put_att(nc_ofid,NH4_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,NH4_vid,"type","data")							 
						 
	case("NO3") !!===============NO3====================================!
    ierr = nf90_def_var(nc_ofid,"NO3",nf90_float,dynm3dn_lay,NO3_vid)
    ierr = nf90_put_att(nc_ofid,NO3_vid,"long_name","Nitrate Concentration")
    ierr = nf90_put_att(nc_ofid,NO3_vid,"standard_name","NO3")
    ierr = nf90_put_att(nc_ofid,NO3_vid,"units","gN m^-3")
    ierr = nf90_put_att(nc_ofid,NO3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,NO3_vid,"type","data")							 
						 
	case("PO4") !!===============PO4====================================!
    ierr = nf90_def_var(nc_ofid,"PO4",nf90_float,dynm3dn_lay,PO4_vid)
    ierr = nf90_put_att(nc_ofid,PO4_vid,"long_name","Phosphate Concentration")
    ierr = nf90_put_att(nc_ofid,PO4_vid,"standard_name","PO4")
    ierr = nf90_put_att(nc_ofid,PO4_vid,"units","gP m^-3")
    ierr = nf90_put_att(nc_ofid,PO4_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PO4_vid,"type","data")							 
						 
	case("PAR") !!===============PAR====================================!
    ierr = nf90_def_var(nc_ofid,"PAR",nf90_float,dynm3dn_lay,PAR_vid)
    ierr = nf90_put_att(nc_ofid,PAR_vid,"long_name","Photosynthetically Active Radiation")
    ierr = nf90_put_att(nc_ofid,PAR_vid,"standard_name","PAR")
    ierr = nf90_put_att(nc_ofid,PAR_vid,"units","mols Photons m^-2 d^-1")
    ierr = nf90_put_att(nc_ofid,PAR_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PAR_vid,"type","data")							 
						 
	case("LPOC") !!===============LPOC====================================!
    ierr = nf90_def_var(nc_ofid,"LPOC",nf90_float,dynm3dn_lay,LPOC_vid)
    ierr = nf90_put_att(nc_ofid,LPOC_vid,"long_name","Labile POC")
    ierr = nf90_put_att(nc_ofid,LPOC_vid,"standard_name","LPOC")
    ierr = nf90_put_att(nc_ofid,LPOC_vid,"units","gC m^-3")
    ierr = nf90_put_att(nc_ofid,LPOC_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,LPOC_vid,"type","data")								 

	case("RPOC") !!===============RPOC====================================!
    ierr = nf90_def_var(nc_ofid,"RPOC",nf90_float,dynm3dn_lay,RPOC_vid)
    ierr = nf90_put_att(nc_ofid,RPOC_vid,"long_name","Refractory POC")
    ierr = nf90_put_att(nc_ofid,RPOC_vid,"standard_name","RPOC")
    ierr = nf90_put_att(nc_ofid,RPOC_vid,"units","gC m^-3")
    ierr = nf90_put_att(nc_ofid,RPOC_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,RPOC_vid,"type","data")	
						
	case("LPON") !!===============LPON====================================!
    ierr = nf90_def_var(nc_ofid,"LPON",nf90_float,dynm3dn_lay,LPON_vid)
    ierr = nf90_put_att(nc_ofid,LPON_vid,"long_name","Labile PON")
    ierr = nf90_put_att(nc_ofid,LPON_vid,"standard_name","LPON")
    ierr = nf90_put_att(nc_ofid,LPON_vid,"units","gN m^-3")
    ierr = nf90_put_att(nc_ofid,LPON_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,LPON_vid,"type","data")								 

	case("RPON") !!===============RPON====================================!
    ierr = nf90_def_var(nc_ofid,"RPON",nf90_float,dynm3dn_lay,RPON_vid)
    ierr = nf90_put_att(nc_ofid,RPON_vid,"long_name","Refractory PON")
    ierr = nf90_put_att(nc_ofid,RPON_vid,"standard_name","RPON")
    ierr = nf90_put_att(nc_ofid,RPON_vid,"units","gN m^-3")
    ierr = nf90_put_att(nc_ofid,RPON_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,RPON_vid,"type","data")							
												
	case("KD") !!===============KD====================================!
    ierr = nf90_def_var(nc_ofid,"KD",nf90_float,dynm3dn_lay,KD_vid)
    ierr = nf90_put_att(nc_ofid,KD_vid,"long_name","Diffuse Attenuation Coef.")
    ierr = nf90_put_att(nc_ofid,KD_vid,"standard_name","KD")
    ierr = nf90_put_att(nc_ofid,KD_vid,"units","m^-1")
    ierr = nf90_put_att(nc_ofid,KD_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,KD_vid,"type","data")							
						
	case("SSI") !!===============SSI====================================!
    ierr = nf90_def_var(nc_ofid,"SSI",nf90_float,dynm3dn_lay,SSI_vid)
    ierr = nf90_put_att(nc_ofid,SSI_vid,"long_name","Suspended Sediment")
    ierr = nf90_put_att(nc_ofid,SSI_vid,"standard_name","SSI")
    ierr = nf90_put_att(nc_ofid,SSI_vid,"units","g m^-3")
    ierr = nf90_put_att(nc_ofid,SSI_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,SSI_vid,"type","data")	
						
	case("RHO1") !!===============RHO1====================================!
    ierr = nf90_def_var(nc_ofid,"RHO1",nf90_float,dynm3dn_lay,RHO1_vid)
    ierr = nf90_put_att(nc_ofid,RHO1_vid,"long_name","Density")
    ierr = nf90_put_att(nc_ofid,RHO1_vid,"standard_name","RHO1")
    ierr = nf90_put_att(nc_ofid,RHO1_vid,"units","kg m^-3")
    ierr = nf90_put_att(nc_ofid,RHO1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,RHO1_vid,"type","data")	
						
	case("TAU") !!===============TAU====================================!
    ierr = nf90_def_var(nc_ofid,"TAU",nf90_float,dynm2dn,TAU_vid)
    ierr = nf90_put_att(nc_ofid,TAU_vid,"long_name","Shear Stress")
    ierr = nf90_put_att(nc_ofid,TAU_vid,"standard_name","Tau")
    ierr = nf90_put_att(nc_ofid,TAU_vid,"units","Pa")
    ierr = nf90_put_att(nc_ofid,TAU_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,TAU_vid,"type","data")	
						
! B Clark added below terms for budget calculations
#if defined (BUDGET)					
!!!!!! NITROGEN						

	case("ALGDON") !!===============Algae DON====================================!
    ierr = nf90_def_var(nc_ofid,"ALGDON",nf90_float,dynm3dn_lay,ALGDON_vid)
    ierr = nf90_put_att(nc_ofid,ALGDON_vid,"long_name","Algae Derived DON Production")
    ierr = nf90_put_att(nc_ofid,ALGDON_vid,"standard_name","ALGDON")
    ierr = nf90_put_att(nc_ofid,ALGDON_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,ALGDON_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,ALGDON_vid,"type","data")	
	
	case("ALGPON") !!===============Algae PON====================================!
    ierr = nf90_def_var(nc_ofid,"ALGPON",nf90_float,dynm3dn_lay,ALGPON_vid)
    ierr = nf90_put_att(nc_ofid,ALGPON_vid,"long_name","Algae Derived PON Production")
    ierr = nf90_put_att(nc_ofid,ALGPON_vid,"standard_name","ALGPON")
    ierr = nf90_put_att(nc_ofid,ALGPON_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,ALGPON_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,ALGPON_vid,"type","data")

	case("ALGNH4") !!===============Algae NH4====================================!
    ierr = nf90_def_var(nc_ofid,"ALGNH4",nf90_float,dynm3dn_lay,ALGNH4_vid)
    ierr = nf90_put_att(nc_ofid,ALGNH4_vid,"long_name","Algae Derived NH4 Uptake")
    ierr = nf90_put_att(nc_ofid,ALGNH4_vid,"standard_name","ALGNH4")
    ierr = nf90_put_att(nc_ofid,ALGNH4_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,ALGNH4_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,ALGNH4_vid,"type","data")
	
	case("ALGNO3") !!===============Algae NO3====================================!
    ierr = nf90_def_var(nc_ofid,"ALGNO3",nf90_float,dynm3dn_lay,ALGNO3_vid)
    ierr = nf90_put_att(nc_ofid,ALGNO3_vid,"long_name","Algae Derived NO3 Uptake")
    ierr = nf90_put_att(nc_ofid,ALGNO3_vid,"standard_name","ALGNO3")
    ierr = nf90_put_att(nc_ofid,ALGNO3_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,ALGNO3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,ALGNO3_vid,"type","data")
		
	case("DENNO3") !!===============Denitrification NO3====================================!
    ierr = nf90_def_var(nc_ofid,"DENNO3",nf90_float,dynm3dn_lay,DENNO3_vid)
    ierr = nf90_put_att(nc_ofid,DENNO3_vid,"long_name","Denitrification NO3 Loss")
    ierr = nf90_put_att(nc_ofid,DENNO3_vid,"standard_name","DENNO3")
    ierr = nf90_put_att(nc_ofid,DENNO3_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,DENNO3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,DENNO3_vid,"type","data")
	
	case("NT") !!===============Denitrification NO3====================================!
    ierr = nf90_def_var(nc_ofid,"NT",nf90_float,dynm3dn_lay,NT_vid)
    ierr = nf90_put_att(nc_ofid,NT_vid,"long_name","Nitrification NH4 Loss")
    ierr = nf90_put_att(nc_ofid,NT_vid,"standard_name","NT")
    ierr = nf90_put_att(nc_ofid,NT_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,NT_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,NT_vid,"type","data")	
	
	case("DENITN") !!===============Denitrification DON Production====================================!
    ierr = nf90_def_var(nc_ofid,"DENITN",nf90_float,dynm3dn_lay,DENITN_vid)
    ierr = nf90_put_att(nc_ofid,DENITN_vid,"long_name","Denitrification DON Loss")
    ierr = nf90_put_att(nc_ofid,DENITN_vid,"standard_name","DENITN")
    ierr = nf90_put_att(nc_ofid,DENITN_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,DENITN_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,DENITN_vid,"type","data")	

	case("HDRLPON") !!===============Hydrolysis of LPON====================================!
    ierr = nf90_def_var(nc_ofid,"HDRLPON",nf90_float,dynm3dn_lay,HDRLPON_vid)
    ierr = nf90_put_att(nc_ofid,HDRLPON_vid,"long_name","Hydrolysis of LPON")
    ierr = nf90_put_att(nc_ofid,HDRLPON_vid,"standard_name","HDRLPON")
    ierr = nf90_put_att(nc_ofid,HDRLPON_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,HDRLPON_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,HDRLPON_vid,"type","data")	
	
	
	case("HDRRPON") !!===============Hydrolysis of RPON====================================!
    ierr = nf90_def_var(nc_ofid,"HDRRPON",nf90_float,dynm3dn_lay,HDRRPON_vid)
    ierr = nf90_put_att(nc_ofid,HDRRPON_vid,"long_name","Hydrolysis of RPON")
    ierr = nf90_put_att(nc_ofid,HDRRPON_vid,"standard_name","HDRRPON")
    ierr = nf90_put_att(nc_ofid,HDRRPON_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,HDRRPON_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,HDRRPON_vid,"type","data")		

	case("COAGN") !!===============Coagulation of DON====================================!
    ierr = nf90_def_var(nc_ofid,"COAGN",nf90_float,dynm3dn_lay,COAGN_vid)
    ierr = nf90_put_att(nc_ofid,COAGN_vid,"long_name","Coagulation of DON")
    ierr = nf90_put_att(nc_ofid,COAGN_vid,"standard_name","COAGN")
    ierr = nf90_put_att(nc_ofid,COAGN_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,COAGN_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,COAGN_vid,"type","data")		

	case("MNLDON1") !!===============Remineralization of DON1====================================!
    ierr = nf90_def_var(nc_ofid,"MNLDON1",nf90_float,dynm3dn_lay,MNLDON1_vid)
    ierr = nf90_put_att(nc_ofid,MNLDON1_vid,"long_name","Remin. of DON1")
    ierr = nf90_put_att(nc_ofid,MNLDON1_vid,"standard_name","MNLDON1")
    ierr = nf90_put_att(nc_ofid,MNLDON1_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,MNLDON1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,MNLDON1_vid,"type","data")	
	
	case("MNLDON2") !!===============Remineralization of DON2====================================!
    ierr = nf90_def_var(nc_ofid,"MNLDON2",nf90_float,dynm3dn_lay,MNLDON2_vid)
    ierr = nf90_put_att(nc_ofid,MNLDON2_vid,"long_name","Remin. of DON2")
    ierr = nf90_put_att(nc_ofid,MNLDON2_vid,"standard_name","MNLDON2")
    ierr = nf90_put_att(nc_ofid,MNLDON2_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,MNLDON2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,MNLDON2_vid,"type","data")		
	
	
	case("MNLDON3") !!===============Remineralization of DON3====================================!
    ierr = nf90_def_var(nc_ofid,"MNLDON3",nf90_float,dynm3dn_lay,MNLDON3_vid)
    ierr = nf90_put_att(nc_ofid,MNLDON3_vid,"long_name","Remin. of DON3")
    ierr = nf90_put_att(nc_ofid,MNLDON3_vid,"standard_name","MNLDON3")
    ierr = nf90_put_att(nc_ofid,MNLDON3_vid,"units","g N m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,MNLDON3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,MNLDON3_vid,"type","data")	
	
!Carbon budget terms	
	
	case("ALGDOC") !!===============Algae DOC====================================!
    ierr = nf90_def_var(nc_ofid,"ALGDOC",nf90_float,dynm3dn_lay,ALGDOC_vid)
    ierr = nf90_put_att(nc_ofid,ALGDOC_vid,"long_name","Algae Derived DOC Production")
    ierr = nf90_put_att(nc_ofid,ALGDOC_vid,"standard_name","ALGDOC")
    ierr = nf90_put_att(nc_ofid,ALGDOC_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,ALGDOC_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,ALGDOC_vid,"type","data")	
	
	case("ALGPOC") !!===============Algae POC====================================!
    ierr = nf90_def_var(nc_ofid,"ALGPOC",nf90_float,dynm3dn_lay,ALGPOC_vid)
    ierr = nf90_put_att(nc_ofid,ALGPOC_vid,"long_name","Algae Derived POC Production")
    ierr = nf90_put_att(nc_ofid,ALGPOC_vid,"standard_name","ALGPOC")
    ierr = nf90_put_att(nc_ofid,ALGPOC_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,ALGPOC_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,ALGPOC_vid,"type","data")	
	
	case("DENIT") !!===============Denitrification DOC Loss====================================!
    ierr = nf90_def_var(nc_ofid,"DENIT",nf90_float,dynm3dn_lay,DENIT_vid)
    ierr = nf90_put_att(nc_ofid,DENIT_vid,"long_name","Denitrification DOC Loss")
    ierr = nf90_put_att(nc_ofid,DENIT_vid,"standard_name","DENIT")
    ierr = nf90_put_att(nc_ofid,DENIT_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,DENIT_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,DENIT_vid,"type","data")	
	
	case("HDRLPOC") !!===============Hydrolysis of LPOC====================================!
    ierr = nf90_def_var(nc_ofid,"HDRLPOC",nf90_float,dynm3dn_lay,HDRLPOC_vid)
    ierr = nf90_put_att(nc_ofid,HDRLPOC_vid,"long_name","Hydrolysis of LPOC")
    ierr = nf90_put_att(nc_ofid,HDRLPOC_vid,"standard_name","HDRLPOC")
    ierr = nf90_put_att(nc_ofid,HDRLPOC_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,HDRLPOC_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,HDRLPOC_vid,"type","data")	
	
	
	case("HDRRPOC") !!===============Hydrolysis of RPOC====================================!
    ierr = nf90_def_var(nc_ofid,"HDRRPOC",nf90_float,dynm3dn_lay,HDRRPOC_vid)
    ierr = nf90_put_att(nc_ofid,HDRRPOC_vid,"long_name","Hydrolysis of RPOC")
    ierr = nf90_put_att(nc_ofid,HDRRPOC_vid,"standard_name","HDRRPOC")
    ierr = nf90_put_att(nc_ofid,HDRRPOC_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,HDRRPOC_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,HDRRPOC_vid,"type","data")	
	
	case("COAGC") !!===============Coagulation of DOC====================================!
    ierr = nf90_def_var(nc_ofid,"COAGC",nf90_float,dynm3dn_lay,COAGC_vid)
    ierr = nf90_put_att(nc_ofid,COAGC_vid,"long_name","Coagulation of DOC")
    ierr = nf90_put_att(nc_ofid,COAGC_vid,"standard_name","COAGC")
    ierr = nf90_put_att(nc_ofid,COAGC_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,COAGC_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,COAGC_vid,"type","data")	
	
	case("MNLDOC1") !!===============Remineralization of DOC1====================================!
    ierr = nf90_def_var(nc_ofid,"MNLDOC1",nf90_float,dynm3dn_lay,MNLDOC1_vid)
    ierr = nf90_put_att(nc_ofid,MNLDOC1_vid,"long_name","Remin. of DOC1")
    ierr = nf90_put_att(nc_ofid,MNLDOC1_vid,"standard_name","MNLDOC1")
    ierr = nf90_put_att(nc_ofid,MNLDOC1_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,MNLDOC1_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,MNLDOC1_vid,"type","data")	
	
	case("MNLDOC2") !!===============Remineralization of DOC2====================================!
    ierr = nf90_def_var(nc_ofid,"MNLDOC2",nf90_float,dynm3dn_lay,MNLDOC2_vid)
    ierr = nf90_put_att(nc_ofid,MNLDOC2_vid,"long_name","Remin. of DOC2")
    ierr = nf90_put_att(nc_ofid,MNLDOC2_vid,"standard_name","MNLDOC2")
    ierr = nf90_put_att(nc_ofid,MNLDOC2_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,MNLDOC2_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,MNLDOC2_vid,"type","data")		
	
	
	case("MNLDOC3") !!===============Remineralization of DOC3====================================!
    ierr = nf90_def_var(nc_ofid,"MNLDOC3",nf90_float,dynm3dn_lay,MNLDOC3_vid)
    ierr = nf90_put_att(nc_ofid,MNLDOC3_vid,"long_name","Remin. of DOC3")
    ierr = nf90_put_att(nc_ofid,MNLDOC3_vid,"standard_name","MNLDOC3")
    ierr = nf90_put_att(nc_ofid,MNLDOC3_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,MNLDOC3_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,MNLDOC3_vid,"type","data")	
	
!photochemistry	
	case("PD31") !!===============Photochem DOC3 to DOC1====================================!
    ierr = nf90_def_var(nc_ofid,"PD31",nf90_float,dynm3dn_lay,PD31_vid)
    ierr = nf90_put_att(nc_ofid,PD31_vid,"long_name","Photochem of 3 to 1")
    ierr = nf90_put_att(nc_ofid,PD31_vid,"standard_name","PD31")
    ierr = nf90_put_att(nc_ofid,PD31_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,PD31_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PD31_vid,"type","data")		

	case("PD30") !!===============Photochem DOC3 to DIC====================================!
    ierr = nf90_def_var(nc_ofid,"PD30",nf90_float,dynm3dn_lay,PD30_vid)
    ierr = nf90_put_att(nc_ofid,PD30_vid,"long_name","Photochem of 3 to DIC")
    ierr = nf90_put_att(nc_ofid,PD30_vid,"standard_name","PD30")
    ierr = nf90_put_att(nc_ofid,PD30_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,PD30_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PD30_vid,"type","data")	
	
	case("PD3N") !!===============Photochem DOC3 to NCDOC====================================!
    ierr = nf90_def_var(nc_ofid,"PD3N",nf90_float,dynm3dn_lay,PD3N_vid)
    ierr = nf90_put_att(nc_ofid,PD3N_vid,"long_name","Photochem of 3 to NCDOC")
    ierr = nf90_put_att(nc_ofid,PD3N_vid,"standard_name","PD3N")
    ierr = nf90_put_att(nc_ofid,PD3N_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,PD3N_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PD3N_vid,"type","data")		

	case("PD21") !!===============Photochem DOC2 to DOC1====================================!
    ierr = nf90_def_var(nc_ofid,"PD21",nf90_float,dynm3dn_lay,PD21_vid)
    ierr = nf90_put_att(nc_ofid,PD21_vid,"long_name","Photochem of 2 to 1")
    ierr = nf90_put_att(nc_ofid,PD21_vid,"standard_name","PD21")
    ierr = nf90_put_att(nc_ofid,PD21_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,PD21_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PD21_vid,"type","data")		

	case("PD20") !!===============Photochem DOC2 to DIC====================================!
    ierr = nf90_def_var(nc_ofid,"PD20",nf90_float,dynm3dn_lay,PD20_vid)
    ierr = nf90_put_att(nc_ofid,PD20_vid,"long_name","Photochem of 2 to DIC")
    ierr = nf90_put_att(nc_ofid,PD20_vid,"standard_name","PD20")
    ierr = nf90_put_att(nc_ofid,PD20_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,PD20_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PD20_vid,"type","data")	
	
	case("PD2N") !!===============Photochem DOC2 to NCDOC====================================!
    ierr = nf90_def_var(nc_ofid,"PD2N",nf90_float,dynm3dn_lay,PD2N_vid)
    ierr = nf90_put_att(nc_ofid,PD2N_vid,"long_name","Photochem of 2 to NCDOC")
    ierr = nf90_put_att(nc_ofid,PD2N_vid,"standard_name","PD2N")
    ierr = nf90_put_att(nc_ofid,PD2N_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,PD2N_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PD2N_vid,"type","data")
	
	
	case("PD10") !!===============Photochem DOC1 to DIC====================================!
    ierr = nf90_def_var(nc_ofid,"PD10",nf90_float,dynm3dn_lay,PD10_vid)
    ierr = nf90_put_att(nc_ofid,PD10_vid,"long_name","Photochem of 1 to DIC")
    ierr = nf90_put_att(nc_ofid,PD10_vid,"standard_name","PD10")
    ierr = nf90_put_att(nc_ofid,PD10_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,PD10_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PD10_vid,"type","data")	
	
	case("PD1N") !!===============Photochem DOC1 to NCDOC====================================!
    ierr = nf90_def_var(nc_ofid,"PD1N",nf90_float,dynm3dn_lay,PD1N_vid)
    ierr = nf90_put_att(nc_ofid,PD1N_vid,"long_name","Photochem of 1 to NCDOC")
    ierr = nf90_put_att(nc_ofid,PD1N_vid,"standard_name","PD1N")
    ierr = nf90_put_att(nc_ofid,PD1N_vid,"units","g C m^-3 d^-1")
    ierr = nf90_put_att(nc_ofid,PD1N_vid,"grid","fvcom_grid")
    ierr = nf90_put_att(nc_ofid,PD1N_vid,"type","data")	
	
	
#endif
	
     case default 
       write(ipt,*)'variable',cdf_vdp(i),' not set up for netcdf output'
       write(ipt,*)'modify module mod_ncdio.F' 
       call pstop
     end select

   end do

!--Exit Define Mode
   
   ierr = nf90_enddef(nc_ofid)
   ierr = nf90_close(nc_ofid)

   end if !(msr)

!==============================================================================|
!  WRITE VARIABLES TO FILE                                                     |
!==============================================================================|
   if(msr)then
     ierr = nf90_open(cdfname,nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error",ipt)
     end if
   end if
   
   !!====Longitude at Nodes (LON) ==========================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,m,mgl,1,1,"n",vx+vxmin,nc_ofid,lon_vid,myid&
        &,nprocs,ipt, stck_cnt)

   !!====Latitude  at Nodes (LAT) ==========================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,m,mgl,1,1,"n",vy+vymin,nc_ofid,lat_vid,myid&
        &,nprocs,ipt, stck_cnt) 

   !!====Number of Processors (NPROCS) =====================!
   if(msr)then 
   ierr = nf90_put_var(nc_ofid,nprocs_vid,nprocs)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing nprocs variable to netcdf",ipt)
   end if
#  if defined (MULTIPROCESSOR)
   ierr = nf90_put_var(nc_ofid,partition_vid,el_pid)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing el_pid variable to netcdf",ipt)
   end if
#  endif
   end if

   !!!==== Density Field==============================!
   !i1 = lbound(rho1,1) ; i2 = ubound(rho1,1)
   !call putvar(i1,i2,m,mgl,kb,kb-1,"n",rho1,nc_ofid,idens_vid,myid&
   !     &,nprocs,ipt, stck_cnt) 

   !!====X Grid Coordinate at Nodes (VX)====================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,m,mgl,1,1,"n",vx+vxmin,nc_ofid,x_vid,myid,nprocs&
        &,ipt, stck_cnt) 

   !!====Y Grid Coordinate at Nodes (VY)====================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,m,mgl,1,1,"n",vy+vymin,nc_ofid,y_vid,myid,nprocs&
        &,ipt, stck_cnt) 

   !!====Bathymetry at Nodes (H)============================!
   i1 = lbound(h,1) ; i2 = ubound(h,1)
   call putvar(i1,i2,m,mgl,1,1,"n",h,nc_ofid,h_vid,myid,nprocs,ipt,&
        & stck_cnt) 

   !!====Nodes surrounding each Element (NV)================!
!   allocate(tmpint(0:ngl,3))
!   if(serial)then
!     tmpint(0:ngl,1:3) = nv(0:ngl,1:3) 
!   end if
!#  if defined (MULTIPROCESSOR)
!   if(par)then
!   do j=1,3
!   do i=1,nloc
!     tmpint(i,j) = egid(nv(i,j))
!   end do
!   end do
!   end if
!#  endif
   !write(*,*)'NV= ',NV
      !!====Nodes surrounding each Element (NV)================!
   allocate(tmp(0:nloc,3))
   if(serial)then
     tmp(0:nloc,1:3) = real(nv(0:nloc,1:3),sp) 
   end if
#  if defined (MULTIPROCESSOR)
   if(par)then
   do j=1,3
   do i=1,nloc
     tmp(i,j) = real(ngid(nv(i,j)),sp)
   end do
   end do
   end if
#  endif
   i1 = lbound(tmp,1) ; i2 = ubound(tmp,1)
   call putvar(i1,i2,n,ngl,3,3,"e",tmp,nc_ofid,nv_vid,myid,nprocs,ipt, stck_cnt)
   deallocate(tmp)
   !i1 = lbound(NV,1) ; i2 = ubound(NV,1)
   !call putvar(i1,i2,n,ngl,3,3,"e",NV,nc_ofid,nv_vid,myid,nprocs&
   !     &,ipt, stck_cnt) 
   !deallocate(tmpint)
!   !!====Momentum Stencil Interpolation Coefficients========!
!   i1 = lbound(a1u,1) ; i2 = ubound(a1u,1)
!   call putvar(i1,i2,n,ngl,4,4,"e",a1u,nc_ofid,a1u_vid,myid,nprocs&
!        &,ipt, stck_cnt) 
!   i1 = lbound(a2u,1) ; i2 = ubound(a2u,1)
!   call putvar(i1,i2,n,ngl,4,4,"e",a2u,nc_ofid,a2u_vid,myid,nprocs&
!        &,ipt, stck_cnt) 

!   !!====Element Based Interpolation Coefficients===========!
!   i1 = lbound(aw0,1) ; i2 = ubound(aw0,1)
!   call putvar(i1,i2,n,ngl,3,3,"e",aw0,nc_ofid,aw0_vid,myid,nprocs&
!        &,ipt, stck_cnt) 
!   i1 = lbound(awx,1) ; i2 = ubound(awx,1)
!   call putvar(i1,i2,n,ngl,3,3,"e",awx,nc_ofid,awx_vid,myid,nprocs&
!        &,ipt, stck_cnt) 
!   i1 = lbound(awy,1) ; i2 = ubound(awy,1)
!   call putvar(i1,i2,n,ngl,3,3,"e",awy,nc_ofid,awy_vid,myid,nprocs&
!        &,ipt, stck_cnt) 

   !!====Sigma Layers (z)==================================!
   i1 = lbound(siglay,1) ; i2 = ubound(siglay,1)
   call putvar(i1,i2,m,mgl,kb-1,kb-1,"n",siglay,nc_ofid,siglay_vid,myid&
        &,nprocs,ipt, stck_cnt) 

   !!====Sigma Layers Shift(zz)==================================!
   !allocate(tmp(0:mt,kbm1))
   !tmp(:,1:kbm1) = z(:,2:kb)
   !i1 = lbound(tmp,1) ; i2 = ubound(tmp,1)
   !call putvar(i1,i2,m,mgl,kb-1,kb-1,"n",dz,nc_ofid,siglay_shift_vid&
   !     &,myid,nprocs,ipt, stck_cnt) 
   !deallocate(tmp)

   !!====Sigma Levels (z)=================================!
   i1 = lbound(siglev,1) ; i2 = ubound(siglev,1)
   call putvar(i1,i2,m,mgl,kb,kb,"n",siglev,nc_ofid,siglev_vid,myid,nprocs&
        &,ipt, stck_cnt) 

   
!==============================================================================|
!  close the file                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)

   return
   end subroutine write_netcdf_setup
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   subroutine out_netcdf 
!==============================================================================|
!   Write Time Dependent NetCDF Data to File                                   |
!==============================================================================|

   !use netcdf
   
		
    !USE MOD_PREC, ONLY: SP
    !USE MOD_SIZES, ONLY: NCP, NSSFP, MGL
    !USE MOD_LIMS, ONLY: MYID, NPROCS, IINT, NTLOC, MLOC
    !USE MOD_TGE, ONLY: NV, TGE_DEALLOC
    
    USE MOD_HYDROVARS
        
    USE MOD_BCMAP, ONLY :    &
         
            BCMAP_DEALLOC

       !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY :         &
            SERIAL          &           !!TRUE IF SINGLE PROCESSOR
            ,MSR            &           !!TRUE IF MASTER PROCESSOR (MYID==1)
            ,PAR            &           !!TRUE IF MULTIPROCESSOR RUN
            ,CASENAME     !    &           !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
         
            
       USE MOD_FILEINFO, ONLY : &
                DIA                &
                ,INIT_FILE_INFO

    !USE MOD_WQM
	

   !USE MOD_SED!, ONLY: JPOC, JPON, JPOP, JPOS, SODTM1S,                   							&
   !!               JNH4TM1S,	 JNO3TM1S,	BENDEN,	   JCH4TM1S, JCH4GTM1S, JHSTM1S, JPO4TM1S, JSITM1S,		&
   !               NH41TM1S,	 NH42TM1S,  NO31TM1S,  NO32TM1S, PO41TM1S,  PO42TM1S, SI1TM1S, SI2TM1S,     &
!				  CH41TM1S,  CH42TM1S,  HS1TM1S,   HS2TM1S,  											&
   !               CPOP,	     CPON,	    CPOC,	   CPOS,	 HSED1,	    BENSTRTM1S			

	!USE WC_DOM, ONLY: MNLDON1, MNLDON1_GL, COAGC, 
	!				  !MNLDOC, &
					  !MNLDOC1,MNLDOC2,MNLDOC3,   &!
					 
					  !MNLDON1,MNLDON2,MNLDON3,   &!

					  !DTWCDOC1,DTWCDOC2,DTWCDOC3,  &
					  !DTWNCDOC1,DTWNCDOC2,DTWNCDOC3,  &
					  !DTWCDON1,DTWCDON2,DTWCDON3,  &
					  !DTWNCDON1,DTWNCDON2,DTWNCDON3, &
					  !DTWCDOP1,DTWCDOP2,DTWCDOP3,  &
					  !DTWNCDOP1,DTWNCDOP2,DTWNCDOP3, &
					  !DTNH4,DTNO3, DTLPON,DTLPOC,DTRPON,DTRPOC, &
					  !WC_CDOC1, WC_CDOC2, WC_CDOC3, & 
					  !WC_NCDOC1, WC_NCDOC2, WC_NCDOC3, &
					  !WC_CDON1, WC_CDON2, WC_CDON3, & 
					  !WC_NCDON1, WC_NCDON2, WC_NCDON3, &
					  !WC_CDOP1, WC_CDOP2, WC_CDOP3, & 
					  !WC_NCDOP1, WC_NCDOP2, WC_NCDOP3
	!				  
	!				  
    !USE MOD_OWQ, ONLY: KDLAMB_int
	
	!USE MOD_KIN
			
!#if defined(SEDIMENT_DOM)			
			!USE MOD_SED_DOM, ONLY : SED_DOM, &        ! B Clark DOM Flag
			!						SED_DOM_OUTPUT,JDOM_OUT, & 
			!						CDOC1, CDOC2, CDOC3,  &
			!						NCDOC1, NCDOC2, NCDOC3, &
			!						CDON1, CDON2, CDON3, &
			!						NCDON1,NCDON2,NCDON3, &
			!						CDOP1, CDOP2, CDOP3,  &
			!						NCDOP1, NCDOP2, NCDOP3
			!						
!#endif
									
!#if defined(PHOTODEGRADATION)
!	USE WC_DOM, ONLY : PHOTOCHEM_DOC_ACCUM, PHOTOCHEM_DOC_ACCUM_GL								
!									
!#endif
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY:   &      
            !EL_PID,    &        !(:)           !!PROCESSOR OWNER OF GLOBAL ELEMENT
            !ELID,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT
            !NLID,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE 
            !ELID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
            !NLID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
            !OBN_GL2LOC,    &    !(:)           !!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
            !EGID,        &    !(:)           !!GLOBAL ID OF LOCAL ELEMENT
            !NGID,        &    !(:)           !!GLOBAL ID OF LOCAL NODE 
            !NHE,        &   !            !NUMBER OF HALO ELEMENTS
            !HE_LST,        &    !(:)           !!GLOBAL IDENTITIES OF HALO ELEMENTS
            !HE_OWN,        &    !(:)           !!OWNER OF HALO ELEMENTS
            !NBN ,       &   !              !!NUMBER OF BOUNDARY NODES
            !MX_MLT,     &   !              !!MAX MULTIPLICITY OF BOUNDARY NODES
            !BN_LST,     &    !               !!GLOBAL IDENTITY OF BOUNDARY NODES
            !BN_LOC,     &    !            !!LOCAL IDENTITY OF BOUNDARY NODES
            !BN_MLT,     &    !            !!MULTIPLICITY OF BOUNDARY NODES
            !BN_NEY,     &    !            !!NODE OWNER LIST
            !NDE_ID,     &    !            !! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
            !NHN,          &    !            !!NUMBER OF HALO NODES
            !HN_LST,     &    !            !!LIST OF HALO NODES 
            !HN_OWN,     &    !            !!PRIMARY OWNER OF HALO NODES
            EC,            &    !
            NC,            &    !
            BNC,         &    !
            EMAP,        &    !
            NMAP,         &    !
            !PNE,        &    !(:)        !!NUMBER OF ELEMENTS IN EACH PROC
            !PNN,        &    !(:)        !!NUMBER OF NODES IN EACH PROC
            !PNHE,        &    !(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
            !PNBN,        &    !(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
            !PMBM,        &    !(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
            !PNHN,        &    !(:)       !!NUMBER OF HALO NODES IN EACH PROC
            !FILELOCK,   & !
            !MYFILELOCK, & !
            !FILELOCKED , & !
            !NODE_MATCH, & !
             NODE_MATCH_ALL, &  !
            !EXCHANGE, & !
             EXCHANGE_ALL, & !
             GATHER , & !
            !IGATHER, & !
             INIT_MPI_ENV,&!
            !SORT,    & !
            !PREAD,    & !
            !PWRITE,& !
            !PPRINT, & !
            !LOCK, & !
            !UNLOCK, & !
            !GETLOC,  &!
             PAR_DEALLOC
#  endif
             
   USE MOD_FILENAMES, ONLY: FNAMEPREFIX, FNAMEEXT!, ifindext
   
#  if defined (MULTIPROCESSOR)   
   USE MOD_BUFFERS, ONLY: RCVBUF, SNDBUF
    USE MOD_PAR, ONLY : NGID
   
#  endif

   USE MOD_DF, ONLY:            &
            DFEEDER,            &!
            DFEEDM1S_GL,        &!
            DFEEDM1S,            &!
            !DF_READ,            &!
            !DF_INIT,            &!
            DF_DEALLOC!,        &!
            !DF_ALLOC           



   implicit none
   integer :: i,ierr,i1,i2,k,icheck,m
   integer :: dims(1)
   real(sp), allocatable :: ftemp(:)
   
!# if defined (SEDDIAG) 
	Real (SP), Allocatable :: 	JCIN(:) ! by feng in order to input Jcin, Jnin, Jpin and Jsin
	Real (SP), Allocatable ::   JNIN(:)
	Real (SP), Allocatable ::   JPIN(:)
	Real (SP), Allocatable ::   JSIN(:)
	
	Real (SP), Allocatable ::   POC2G1(:)
	Real (SP), Allocatable ::   POC2G2(:)
	Real (SP), Allocatable ::   POC2G3(:)
	
	Real (SP), Allocatable ::   PON2G1(:)
	Real (SP), Allocatable ::   PON2G2(:)
	Real (SP), Allocatable ::   PON2G3(:)
	
	Real (SP), Allocatable ::   POP2G1(:)
	Real (SP), Allocatable ::   POP2G2(:)
	Real (SP), Allocatable ::   POP2G3(:)

	Real (SP), Allocatable ::   vartmp(:)
	
	Allocate(JCIN(0:MTLOC)); JCIN = 0.0
	Allocate(JNIN(0:MTLOC)); JNIN = 0.0
	Allocate(JPIN(0:MTLOC)); JPIN = 0.0
	Allocate(JSIN(0:MTLOC)); JSIN = 0.0
	
	Allocate(POC2G1(0:MTLOC)); POC2G1 = 0.0
	Allocate(POC2G2(0:MTLOC)); POC2G2 = 0.0
	Allocate(POC2G3(0:MTLOC)); POC2G3 = 0.0

	Allocate(PON2G1(0:MTLOC)); PON2G1 = 0.0
	Allocate(PON2G2(0:MTLOC)); PON2G2 = 0.0
	Allocate(PON2G3(0:MTLOC)); PON2G3 = 0.0

	Allocate(POP2G1(0:MTLOC)); POP2G1 = 0.0
	Allocate(POP2G2(0:MTLOC)); POP2G2 = 0.0
	Allocate(POP2G3(0:MTLOC)); POP2G3 = 0.0
	Allocate(vartmp(0:MTLOC)); vartmp = 0.0;
	
	!convert JCIN POC, PON, POS etc to the output units
	DO I=1,MTLOC
		JCIN(I) = (JPOC(I, 1)+JPOC(I, 2)+JPOC(I, 3)) / 1000.0 !* 2.667	!mgC/m^2/day --> gO2/m^2/day
		JNIN(I) = (JPON(I, 1)+JPON(I, 2)+JPON(I, 3)) / 1000.0			!mgN/m^2/day --> gO2/m^2/day
		JPIN(I) = (JPOP(I, 1)+JPOP(I, 2)+JPOP(I, 3)) / 1000.0			!mgP/m^2/day --> gO2/m^2/day
		
		POC2G1(I) = CPOC(I, 1)/ 1000.0! * 2.667	!mgC/m^3 --> gO2/m^3	
		POC2G2(I) = CPOC(I, 2)/ 1000.0! * 2.667  !mgC/m^3 --> gO2/m^3
		POC2G3(I) = CPOC(I, 3)/ 1000.0! * 2.667  !mgC/m^3 --> gO2/m^3
		
		PON2G1(I) = CPON(I, 1)/ 1000.0			!mgN/m^3 --> gO2/m^3
		PON2G2(I) = CPON(I, 2)/ 1000.0			!mgN/m^3 --> gO2/m^3
		PON2G3(I) = CPON(I, 3)/ 1000.0			!mgN/m^3 --> gO2/m^3
		
		POP2G1(I) = CPOP(I, 1)/ 1000.0			!mgP/m^3 --> gO2/m^3
		POP2G2(I) = CPOP(I, 2)/ 1000.0			!mgP/m^3 --> gO2/m^3
		POP2G3(I) = CPOP(I, 3)/ 1000.0			!mgP/m^3 --> gO2/m^3
		
	END DO
	
	
	!JSIN(1:MTLOC)=JPOS(1:MTLOC)/1000.0			!mgSi/m^2/day --> gSi/m^2/day
	
	
	
!# endif  
   
	
!==============================================================================|
   

!--Update Counter
   out_cnt = out_cnt + 1
   stck_cnt = stck_cnt + 1 
   
   m=mloc

!--Write Header information if first output of file
   if(cdf_stk == 0)then
     if(out_cnt == 1) call write_netcdf_setup(1)
   else
     icheck = mod(out_cnt-1,cdf_stk)
     if(icheck ==0 .or. out_cnt==1)call write_netcdf_setup((out_cnt-1)/cdf_stk+1)
   endif

   dims(1) = stck_cnt

!--Open File
   if(msr)then
     ierr = nf90_open(cdfname,nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error",ipt)
     end if

!--Dump Time/IINT to File
      ierr    = nf90_put_var(nc_ofid,iint_vid,iint,START=dims)
      ierr    = nf90_put_var(nc_ofid,time_vid,thour*3600.,START=dims)
   end if

!--Write Variables to File
   if(msr) write(ipt,*)'dumping to netcdf file: ',trim(cdfname),stck_cnt

   do i=1,nout_vars

      !IF(MSR)WRITE(*,*)'varibale ', i, '=', trim(cdf_vdp(i))

     select case(trim(cdf_vdp(i)))
	 
	case("uu") !!===============uu======================================!
       i1 = lbound(uu,1) ; i2 = ubound(uu,1)
       call putvar(i1,i2,n,ngl,kb,kb-1,"e",uu,nc_ofid,uu_vid,myid,nprocs&
            &,ipt, stck_cnt) 
	   
	case("vv") !!===============vv======================================!
       i1 = lbound(vv,1) ; i2 = ubound(vv,1)
       call putvar(i1,i2,n,ngl,kb,kb-1,"e",vv,nc_ofid,vv_vid,myid,nprocs&
            &,ipt, stck_cnt) 	   
	   
	case("ww") !!===============ww======================================!
       i1 = lbound(ww_out,1) ; i2 = ubound(ww_out,1)
       call putvar(i1,i2,n,ngl,kb,kb-1,"e",ww_out,nc_ofid,ww_vid,myid,nprocs&
            &,ipt,stck_cnt) 
	   
     case("el") !!===============EL======================================!
       i1 = lbound(el,1) ; i2 = ubound(el,1)
       call putvar(i1,i2,m,mgl,1,1,"n",el,nc_ofid,el_vid,myid,nprocs&
            &,ipt, stck_cnt) 

     case("d") !!===============D=======================================!
       i1 = lbound(d,1) ; i2 = ubound(d,1)
       call putvar(i1,i2,m,mgl,1,1,"n",d,nc_ofid,d_vid,myid,nprocs&
            &,ipt, stck_cnt) 

     case("t1") !!===============T1======================================!
       i1 = lbound(T,1) ; i2 = ubound(T,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",T,nc_ofid,t1_vid,myid&
            &,nprocs,ipt, stck_cnt) 

     case("s1") !!===============S1======================================!
       i1 = lbound(salt,1) ; i2 = ubound(salt,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",salt,nc_ofid,s1_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	 !DOC  	   	   
	   case("WC_CDOC1") !!===============WC_CDOC1===============================!
       i1 = lbound(WC_CDOC1,1) ; i2 = ubound(WC_CDOC1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDOC1,nc_ofid,WC_CDOC1_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_CDOC2") !!===============WC_CDOC2===============================!
       i1 = lbound(WC_CDOC2,1) ; i2 = ubound(WC_CDOC2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDOC2,nc_ofid,WC_CDOC2_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_CDOC3") !!===============WC_CDOC3===============================!
       i1 = lbound(WC_CDOC3,1) ; i2 = ubound(WC_CDOC3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDOC3,nc_ofid,WC_CDOC3_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_NCDOC1") !!===============WC_NCDOC1===============================!
       i1 = lbound(WC_NCDOC1,1) ; i2 = ubound(WC_NCDOC1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDOC1,nc_ofid,WC_NCDOC1_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_NCDOC2") !!===============WC_NCDOC2===============================!
       i1 = lbound(WC_NCDOC2,1) ; i2 = ubound(WC_NCDOC2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDOC2,nc_ofid,WC_NCDOC2_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_NCDOC3") !!===============WC_NCDOC3===============================!
       i1 = lbound(WC_NCDOC3,1) ; i2 = ubound(WC_NCDOC3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDOC3,nc_ofid,WC_NCDOC3_vid,myid&
            &,nprocs,ipt, stck_cnt)     
	   !DON
	   case("WC_CDON1") !!===============WC_CDON1===============================!
       i1 = lbound(WC_CDON1,1) ; i2 = ubound(WC_CDON1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDON1,nc_ofid,WC_CDON1_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_CDON2") !!===============WC_CDON2===============================!
       i1 = lbound(WC_CDON2,1) ; i2 = ubound(WC_CDON2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDON2,nc_ofid,WC_CDON2_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_CDON3") !!===============WC_CDON3===============================!
       i1 = lbound(WC_CDON3,1) ; i2 = ubound(WC_CDON3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDON3,nc_ofid,WC_CDON3_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   
	   case("WC_NCDON1") !!===============WC_NCDON1===============================!
       i1 = lbound(WC_NCDON1,1) ; i2 = ubound(WC_NCDON1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDON1,nc_ofid,WC_NCDON1_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_NCDON2") !!===============WC_NCDON2===============================!
       i1 = lbound(WC_NCDON2,1) ; i2 = ubound(WC_NCDON2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDON2,nc_ofid,WC_NCDON2_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_NCDON3") !!===============WC_NCDON3===============================!
       i1 = lbound(WC_NCDON3,1) ; i2 = ubound(WC_NCDON3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDON3,nc_ofid,WC_NCDON3_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   !DOP
	   case("WC_CDOP1") !!===============WC_CDOP1===============================!
       i1 = lbound(WC_CDOP1,1) ; i2 = ubound(WC_CDOP1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDOP1,nc_ofid,WC_CDOP1_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_CDOP2") !!===============WC_CDOP2===============================!
       i1 = lbound(WC_CDOP2,1) ; i2 = ubound(WC_CDOP2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDOP2,nc_ofid,WC_CDOP2_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_CDOP3") !!===============WC_CDOP3===============================!
       i1 = lbound(WC_CDOP3,1) ; i2 = ubound(WC_CDOP3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_CDON3,nc_ofid,WC_CDOP3_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   	   
	   case("WC_NCDOP1") !!===============WC_NCDOP1===============================!
       i1 = lbound(WC_NCDOP1,1) ; i2 = ubound(WC_NCDOP1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDOP1,nc_ofid,WC_NCDOP1_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_NCDOP2") !!===============WC_NCDOP2===============================!
       i1 = lbound(WC_NCDOP2,1) ; i2 = ubound(WC_NCDOP2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDOP2,nc_ofid,WC_NCDOP2_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	   
	   case("WC_NCDOP3") !!===============WC_NCDOP3===============================!
       i1 = lbound(WC_NCDOP3,1) ; i2 = ubound(WC_NCDOP3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",WC_NCDOP3,nc_ofid,WC_NCDOP3_vid,myid&
            &,nprocs,ipt, stck_cnt) 
  	
		!====== added-by-feng, April-2016  ===============================!
	   case("JCIN") !!===============jcin======================================!
       i1 = lbound(jcin,1) ; i2 = ubound(jcin,1)
       call putvar(i1,i2,m,mgl,1,1,"n",jcin,nc_ofid,jcin_vid,myid,nprocs&
            &,ipt, stck_cnt) 
	   
	   case("JNIN") !!===============jnin======================================!
       i1 = lbound(jnin,1) ; i2 = ubound(jnin,1)
       call putvar(i1,i2,m,mgl,1,1,"n",jnin,nc_ofid,jnin_vid,myid,nprocs&
            &,ipt, stck_cnt) 
	   
	   case("JPIN") !!===============jpin======================================!
       i1 = lbound(jpin,1) ; i2 = ubound(jpin,1)
       call putvar(i1,i2,m,mgl,1,1,"n",jpin,nc_ofid,jpin_vid,myid,nprocs&
            &,ipt, stck_cnt) 

	 case("SOD") !!===============SOD======================================!
      i1 = lbound(SODTM1S,1) ; i2 = ubound(SODTM1S,1)
      call putvar(i1,i2,m,mgl,1,1,"n",SODTM1S,nc_ofid,SOD_vid,myid,nprocs&
            &,ipt, stck_cnt) 

	 case("JNH4") !!===============JNH4======================================! 
	  vartmp = JNH4TM1S/1000  !convert from mgN/m^2/day to gN/m2/day
      i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
      call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,JNH4_vid,myid,nprocs&
            &,ipt, stck_cnt)
	                  		  
	 case("JNO3") !!===============JNO3======================================!      
	  vartmp = JNO3TM1S/1000  !convert from mgN/m^2/day to gN/m2/day 
	  i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
      call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,JNO3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JDenitT") !!===============JDenitT======================================!
      i1 = lbound(BENDEN,1) ; i2 = ubound(BENDEN,1)
       call putvar(i1,i2,m,mgl,1,1,"n",BENDEN,nc_ofid,JDenitT_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("JCH4") !!===============JCH4======================================!
     i1 = lbound(JCH4TM1S,1) ; i2 = ubound(JCH4TM1S,1)
     call putvar(i1,i2,m,mgl,1,1,"n",JCH4TM1S,nc_ofid,JCH4_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JCH4G") !!===============JCH4G======================================!
     i1 = lbound(JCH4GTM1S,1) ; i2 = ubound(JCH4GTM1S,1)
     call putvar(i1,i2,m,mgl,1,1,"n",JCH4GTM1S,nc_ofid,JCH4G_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JHS") !!===============JHS======================================!
     i1 = lbound(JHSTM1S,1) ; i2 = ubound(JHSTM1S,1)
     call putvar(i1,i2,m,mgl,1,1,"n",JHSTM1S,nc_ofid,JHS_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("JPO4") !!===============JPO4======================================!
	 vartmp = JPO4TM1S/1000	!convert from mgNP/m^2/day to gP/m2/day
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,JPO4_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("NH41") !!===============NH41 or in name of NH4(1)======================================!
	 vartmp = NH41TM1S/1000	!convert from mgN/m^3 to gN/m3
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,NH41_vid,myid,nprocs&
            &,ipt, stck_cnt)
  	
	 case("NH42") !!===============NH42 or in name of NH4(2)======================================!
	 vartmp= NH42TM1S/1000	!convert from mgN/m^3 to gN/m3
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,NH42_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("NO31") !!===============NO31 or in name of NO3(1)======================================!
	 vartmp = NO31TM1S/1000	!convert from mgN/m^3 to gN/m3
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,NO31_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("NO32") !!===============NO32 or in name of NO3(2)======================================!
	 vartmp = NO32TM1S/1000	!convert from mgN/m^3 to gN/m3
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,NO32_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("PO41") !!===============PO41 or in name of PO4(1)======================================!
	 vartmp = PO41TM1S/1000	!convert from mgP/m^3 to gP/m3
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,PO41_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("PO42") !!===============PO42 or in name of PO4(2)======================================!
	 vartmp = PO42TM1S/1000 !convert from mgP/m^3 to gP/m3
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,PO42_vid,myid,nprocs&
            &,ipt, stck_cnt) 

	 case("CH41") !!===============CH41 or in name of CH4(1)===================================!
	 vartmp = CH41TM1S	   !convert from gO2/m^3 to mgO2/L 
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,CH41_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("CH42") !!===============CH42 or in name of CH4(2)===================================!
	 vartmp = CH42TM1S		
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,CH42_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("HS1") !!===============HS1 or in name of HS(1)======================================!
	 vartmp = HS1TM1S
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,HS1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("HS2") !!===============HS2 or in name of HS(2)======================================!
	 vartmp= HS2TM1S
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,HS2_vid,myid,nprocs&
            &,ipt, stck_cnt)	 
	 
	 case("POC2G1") !!===============POC2G1 ======================================!
     i1 = lbound(POC2G1,1) ; i2 = ubound(POC2G1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",POC2G1,nc_ofid,POC2G1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("POC2G2") !!===============POC2G2 ======================================!
     i1 = lbound(POC2G2,1) ; i2 = ubound(POC2G2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",POC2G2,nc_ofid,POC2G2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("POC2G3") !!===============POC2G3 ======================================!
     i1 = lbound(POC2G3,1) ; i2 = ubound(POC2G3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",POC2G3,nc_ofid,POC2G3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("PON2G1") !!===============PON2G1 ======================================!
     i1 = lbound(PON2G1,1) ; i2 = ubound(PON2G1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",PON2G1,nc_ofid,PON2G1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	case("PON2G2") !!===============PON2G2 ======================================!
     i1 = lbound(PON2G2,1) ; i2 = ubound(PON2G2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",PON2G2,nc_ofid,PON2G2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("PON2G3") !!===============PON2G3 ======================================!
     i1 = lbound(PON2G3,1) ; i2 = ubound(PON2G3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",PON2G3,nc_ofid,PON2G3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("POP2G1") !!===============POP2G1 ======================================!
     i1 = lbound(POP2G1,1) ; i2 = ubound(POP2G1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",POP2G1,nc_ofid,POP2G1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("POP2G2") !!===============POP2G2 ======================================!
     i1 = lbound(POP2G2,1) ; i2 = ubound(POP2G2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",POP2G2,nc_ofid,POP2G2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("POP2G3") !!===============POP2G3 ======================================!
     i1 = lbound(POP2G3,1) ; i2 = ubound(POP2G3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",POP2G3,nc_ofid,POP2G3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 
	 case("H1") !!===============H1 ======================================!
	 vartmp = HSED1
	 i1 = lbound(vartmp,1) ; i2 = ubound(vartmp,1)
	 !i1=0, i2=MT
     call putvar(i1,i2,m,mgl,1,1,"n",vartmp,nc_ofid,H1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 

     Deallocate(JCIN,JNIN,JPIN,JSIN,POC2G1,POC2G2,POC2G3,PON2G1,PON2G2,PON2G3,POP2G1,POP2G2,POP2G3,vartmp)
	 
!#	if defined (SED_DOM)

	 case("SED_CDOC1") !!===============SED_CDOC1 =====================================!     
	 i1 = lbound(CDOC1,1) ; i2 = ubound(CDOC1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDOC1,nc_ofid,SED_CDOC1_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("SED_CDOC2") !!===============SED_CDOC2 =====================================!     
	 i1 = lbound(CDOC2,1) ; i2 = ubound(CDOC2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDOC2,nc_ofid,SED_CDOC2_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("SED_CDOC3") !!===============SED_CDOC3 =====================================!     
	 i1 = lbound(CDOC3,1) ; i2 = ubound(CDOC3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDOC3,nc_ofid,SED_CDOC3_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("SED_NCDOC1") !!===============SED_NCDOC1 =====================================!     
	 i1 = lbound(CDOC1,1) ; i2 = ubound(NCDOC1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDOC1,nc_ofid,SED_NCDOC1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 	 

	 case("SED_NCDOC2") !!===============SED_NCDOC2 =====================================!     
	 i1 = lbound(NCDOC2,1) ; i2 = ubound(NCDOC2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDOC2,nc_ofid,SED_NCDOC2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("SED_NCDOC3") !!===============SED_NCDOC3 =====================================!     
	 i1 = lbound(NCDOC3,1) ; i2 = ubound(NCDOC3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDOC3,nc_ofid,SED_NCDOC3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 	 	 	 
	 !DON
	 
	 case("SED_CDON1") !!===============SED_CDON1 =====================================!     
	 i1 = lbound(CDON1,1) ; i2 = ubound(CDON1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDON1,nc_ofid,SED_CDON1_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("SED_CDON2") !!===============SED_CDON2 =====================================!     
	 i1 = lbound(CDON2,1) ; i2 = ubound(CDON2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDON2,nc_ofid,SED_CDON2_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("SED_CDON3") !!===============SED_CDON3 =====================================!     
	 i1 = lbound(CDON3,1) ; i2 = ubound(CDON3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDON3,nc_ofid,SED_CDON3_vid,myid,nprocs&
            &,ipt, stck_cnt)


	 case("SED_NCDON1") !!===============SED_NCDON1 =====================================!     
	 i1 = lbound(NCDON1,1) ; i2 = ubound(NCDON1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDON1,nc_ofid,SED_NCDON1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 	 

	 case("SED_NCDON2") !!===============SED_NCDON2 =====================================!     
	 i1 = lbound(NCDON2,1) ; i2 = ubound(NCDON2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDON2,nc_ofid,SED_NCDON2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	  case("SED_NCDON3") !!===============SED_NCDON3 =====================================!     
	 i1 = lbound(NCDON3,1) ; i2 = ubound(NCDON3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDON3,nc_ofid,SED_NCDON3_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 !DOP
	 
	 case("SED_CDOP1") !!===============SED_CDOP1 =====================================!     
	 i1 = lbound(CDOP1,1) ; i2 = ubound(CDOP1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDOP1,nc_ofid,SED_CDOP1_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("SED_CDOP2") !!===============SED_CDOP2 =====================================!     
	 i1 = lbound(CDOP2,1) ; i2 = ubound(CDOP2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDOP2,nc_ofid,SED_CDOP2_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("SED_CDOP3") !!===============SED_CDOP3 =====================================!     
	 i1 = lbound(CDOP3,1) ; i2 = ubound(CDOP3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",CDOP3,nc_ofid,SED_CDOP3_vid,myid,nprocs&
            &,ipt, stck_cnt)


	 case("SED_NCDOP1") !!===============SED_NCDOP1 =====================================!     
	 i1 = lbound(NCDOP1,1) ; i2 = ubound(NCDOP1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDOP1,nc_ofid,SED_NCDOP1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 	 

	 case("SED_NCDOP2") !!===============SED_NCDOP2 =====================================!     
	 i1 = lbound(NCDOP2,1) ; i2 = ubound(NCDOP2,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDOP2,nc_ofid,SED_NCDOP2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("SED_NCDOP3") !!===============SED_NCDOP3 =====================================!     
	 i1 = lbound(NCDOP3,1) ; i2 = ubound(NCDOP3,1)
     call putvar(i1,i2,m,mgl,1,1,"n",NCDOP3,nc_ofid,SED_NCDOP3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWCDOC1") !!===============JWCDOC1 =====================================!     
	 i1 = lbound(JDOM_OUT(:,1),1) ; i2 = ubound(JDOM_OUT(:,1),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,1),nc_ofid,JWCDOC1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWCDOC2") !!===============JWCDOC2 =====================================!     
	 i1 = lbound(JDOM_OUT(:,2),1) ; i2 = ubound(JDOM_OUT(:,2),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,2),nc_ofid,JWCDOC2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWCDOC3") !!===============JWCDOC3 =====================================!     
	 i1 = lbound(JDOM_OUT(:,3),1) ; i2 = ubound(JDOM_OUT(:,3),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,3),nc_ofid,JWCDOC3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 	 	 
	 case("JWNCDOC1") !!===============JWNCDOC1 =====================================!     
	 i1 = lbound(JDOM_OUT(:,4),1) ; i2 = ubound(JDOM_OUT(:,4),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,4),nc_ofid,JWNCDOC1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWNCDOC2") !!===============JWNCDOC2 =====================================!     
	 i1 = lbound(JDOM_OUT(:,5),1) ; i2 = ubound(JDOM_OUT(:,5),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,5),nc_ofid,JWNCDOC2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWNCDOC3") !!===============JWNCDOC3 =====================================!     
	 i1 = lbound(JDOM_OUT(:,6),1) ; i2 = ubound(JDOM_OUT(:,6),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,6),nc_ofid,JWNCDOC3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 !DON
	 case("JWCDON1") !!===============JWCDON1 =====================================!     
	 i1 = lbound(JDOM_OUT(:,7),1) ; i2 = ubound(JDOM_OUT(:,7),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,7),nc_ofid,JWCDON1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWCDON2") !!===============JWCDON2 =====================================!     
	 i1 = lbound(JDOM_OUT(:,8),1) ; i2 = ubound(JDOM_OUT(:,8),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,8),nc_ofid,JWCDON2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWCDON3") !!===============JWCDON3 =====================================!     
	 i1 = lbound(JDOM_OUT(:,9),1) ; i2 = ubound(JDOM_OUT(:,9),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,9),nc_ofid,JWCDON3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 	 
	 case("JWNCDON1") !!===============JWNCDON1 =====================================!     
	 i1 = lbound(JDOM_OUT(:,10),1) ; i2 = ubound(JDOM_OUT(:,10),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,10),nc_ofid,JWNCDON1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWNCDON2") !!===============JWNCDON2 =====================================!     
	 i1 = lbound(JDOM_OUT(:,11),1) ; i2 = ubound(JDOM_OUT(:,11),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,11),nc_ofid,JWNCDON2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWNCDON3") !!===============JWNCDON3 =====================================! 
	 i1 = lbound(JDOM_OUT(:,12),1) ; i2 = ubound(JDOM_OUT(:,12),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,12),nc_ofid,JWNCDON3_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 !DOP
	 case("JWCDOP1") !!===============JWCDOP1 =====================================!     
	 i1 = lbound(JDOM_OUT(:,13),1) ; i2 = ubound(JDOM_OUT(:,13),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,13),nc_ofid,JWCDOP1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWCDOP2") !!===============JWCDOP2 =====================================!     
	 i1 = lbound(JDOM_OUT(:,14),1) ; i2 = ubound(JDOM_OUT(:,14),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,14),nc_ofid,JWCDOP2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWCDOP3") !!===============JWCDOP3 =====================================!     
	 i1 = lbound(JDOM_OUT(:,15),1) ; i2 = ubound(JDOM_OUT(:,15),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,15),nc_ofid,JWCDOP3_vid,myid,nprocs&
            &,ipt, stck_cnt)

	 case("JWNCDOP1") !!===============JWNCDOP1 =====================================!     
	 i1 = lbound(JDOM_OUT(:,16),1) ; i2 = ubound(JDOM_OUT(:,16),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,16),nc_ofid,JWNCDOP1_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWNCDOP2") !!===============JWNCDOP2 =====================================!     
	 i1 = lbound(JDOM_OUT(:,17),1) ; i2 = ubound(JDOM_OUT(:,17),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,17),nc_ofid,JWNCDOP2_vid,myid,nprocs&
            &,ipt, stck_cnt)
	 
	 case("JWNCDOP3") !!===============JWNCDOP3 =====================================!     
	 i1 = lbound(JDOM_OUT(:,18),1) ; i2 = ubound(JDOM_OUT(:,18),1)
     call putvar(i1,i2,m,mgl,1,1,"n",JDOM_OUT(:,18),nc_ofid,JWNCDOP3_vid,myid,nprocs&
            &,ipt, stck_cnt)
!#   endif

	! Below added by B Clark for outputting other water quality Variables	
	 
	   case("DOXG") !!===============DOXG===============================!
       i1 = lbound(DOXG,1) ; i2 = ubound(DOXG,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",DOXG,nc_ofid,DOXG_vid,myid&
            &,nprocs,ipt, stck_cnt) 
	 
	 
	 case("NPP") !!===============NPP=====================================!     
	 i1 = lbound(total_netPP,1) ; i2 = ubound(total_netPP,1)
     call putvar(i1,i2,m,mgl,1,1,"n",total_netPP,nc_ofid,NPP_VID,myid,nprocs&
            &,ipt, stck_cnt)			
				
	   case("B1") !!===============Algae 1===============================!
       i1 = lbound(B1,1) ; i2 = ubound(B1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",B1,nc_ofid,B1_vid,myid&
            &,nprocs,ipt, stck_cnt) 			
			
	   case("B2") !!===============Algae 2===============================!
       i1 = lbound(B2,1) ; i2 = ubound(B2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",B2,nc_ofid,B2_vid,myid&
            &,nprocs,ipt, stck_cnt) 	
			
	   case("NH4") !!===============NH4===============================!
       i1 = lbound(NH4,1) ; i2 = ubound(NH4,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",NH4,nc_ofid,NH4_vid,myid&
            &,nprocs,ipt, stck_cnt) 				
			
	  case("NO3") !!===============NO3===============================!
       i1 = lbound(NO3,1) ; i2 = ubound(NO3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",NO3,nc_ofid,NO3_vid,myid&
            &,nprocs,ipt, stck_cnt) 		
			
	  case("PO4") !!===============PO4===============================!
       i1 = lbound(PO4,1) ; i2 = ubound(PO4,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PO4,nc_ofid,PO4_vid,myid&
            &,nprocs,ipt, stck_cnt)			
			
	  case("PAR") !!===============PAR===============================!
       i1 = lbound(IAVG,1) ; i2 = ubound(IAVG,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",IAVG,nc_ofid,PAR_vid,myid&
            &,nprocs,ipt, stck_cnt)				
			
	  case("LPOC") !!===============LPOC===============================!
       i1 = lbound(LPOC,1) ; i2 = ubound(LPOC,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",LPOC,nc_ofid,LPOC_vid,myid&
            &,nprocs,ipt, stck_cnt)				
			
	  case("RPOC") !!===============RPOC===============================!
       i1 = lbound(RPOC,1) ; i2 = ubound(RPOC,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",RPOC,nc_ofid,RPOC_vid,myid&
            &,nprocs,ipt, stck_cnt)					

	  case("LPON") !!===============LPON===============================!
       i1 = lbound(LPON,1) ; i2 = ubound(LPON,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",LPON,nc_ofid,LPON_vid,myid&
            &,nprocs,ipt, stck_cnt)				
			
	  case("RPON") !!===============RPON===============================!
       i1 = lbound(RPON,1) ; i2 = ubound(RPON,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",RPON,nc_ofid,RPON_vid,myid&
            &,nprocs,ipt, stck_cnt)			
		
	  case("KD") !!===============KD===============================!
       i1 = lbound(KDlamb_int,1) ; i2 = ubound(KDlamb_int,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",KDlamb_int,nc_ofid,KD_vid,myid&
            &,nprocs,ipt, stck_cnt)	
			
       case("SSI") !!===============SSI===============================!
       i1 = lbound(SSI,1) ; i2 = ubound(SSI,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",SSI,nc_ofid,SSI_vid,myid&
            &,nprocs,ipt, stck_cnt)				
			
       case("RHO1") !!===============RHO1===============================!
       i1 = lbound(RHO1,1) ; i2 = ubound(RHO1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",RHO1,nc_ofid,RHO1_vid,myid&
            &,nprocs,ipt, stck_cnt)				
			
       case("TAU") !!===============TAU===============================!
       i1 = lbound(ShearStress,1) ; i2 = ubound(ShearStress,1)
       call putvar(i1,i2,m,mgl,1,1,"n",ShearStress,nc_ofid,TAU_vid,myid&
            &,nprocs,ipt, stck_cnt)				

#if defined (BUDGET)
	   
	   case("ALGDON") !!===============ALgae DON===============================!
       i1 = lbound(ALGDON,1) ; i2 = ubound(ALGDON,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",AlgDON,nc_ofid,ALGDON_vid,myid&
            &,nprocs,ipt, stck_cnt)	

	   case("ALGPON") !!===============ALgae PON===============================!
       i1 = lbound(ALGPON,1) ; i2 = ubound(ALGPON,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",AlgPON,nc_ofid,ALGPON_vid,myid&
            &,nprocs,ipt, stck_cnt)		   
	   
	   case("ALGNH4") !!===============ALgae NH4===============================!
       i1 = lbound(ALGNH4,1) ; i2 = ubound(ALGNH4,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",ALGNH4,nc_ofid,ALGNH4_vid,myid&
            &,nprocs,ipt, stck_cnt)		
	   
	   case("ALGNO3") !!===============ALgae NO3===============================!
       i1 = lbound(ALGNO3,1) ; i2 = ubound(ALGNO3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",ALGNO3,nc_ofid,ALGNO3_vid,myid&
            &,nprocs,ipt, stck_cnt)		   
	   
	   case("DENNO3") !!===============Denit. NO3===============================!
       i1 = lbound(DENNO3,1) ; i2 = ubound(DENNO3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",DENNO3,nc_ofid,DENNO3_vid,myid&
            &,nprocs,ipt, stck_cnt)		   
	   
	   case("NT") !!===============Nitrif. NH4===============================!
       i1 = lbound(NT,1) ; i2 = ubound(NT,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",NT,nc_ofid,NT_vid,myid&
            &,nprocs,ipt, stck_cnt)
	   
	   case("DENITN") !!===============Denit. to DON===============================!
       i1 = lbound(DENITN,1) ; i2 = ubound(DENITN,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",DENITN,nc_ofid,DENITN_vid,myid&
            &,nprocs,ipt, stck_cnt)		   

	   case("HDRLPON") !!===============Hydro. LPON===============================!
       i1 = lbound(HDRLPON,1) ; i2 = ubound(HDRLPON,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",HDRLPON,nc_ofid,HDRLPON_vid,myid&
            &,nprocs,ipt, stck_cnt)	   

	   case("HDRRPON") !!===============Hydro. RPON===============================!
       i1 = lbound(HDRRPON,1) ; i2 = ubound(HDRRPON,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",HDRRPON,nc_ofid,HDRRPON_vid,myid&
            &,nprocs,ipt, stck_cnt)		   
	   
	   case("COAGN") !!===============Coag DON===============================!
       i1 = lbound(COAGN,1) ; i2 = ubound(COAGN,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",COAGN,nc_ofid,COAGN_vid,myid&
            &,nprocs,ipt, stck_cnt)	   
	   
	   case("MNLDON1") !!===============Remin DON1===============================!
       i1 = lbound(MNLDON1,1) ; i2 = ubound(MNLDON1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",MNLDON1,nc_ofid,MNLDON1_vid,myid&
            &,nprocs,ipt, stck_cnt)		   

	   case("MNLDON2") !!===============Remin DON2===============================!
       i1 = lbound(MNLDON2,1) ; i2 = ubound(MNLDON2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",MNLDON2,nc_ofid,MNLDON2_vid,myid&
            &,nprocs,ipt, stck_cnt)	   
	   
	   case("MNLDON3") !!===============Remin DON3===============================!
       i1 = lbound(MNLDON3,1) ; i2 = ubound(MNLDON3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",MNLDON3,nc_ofid,MNLDON3_vid,myid&
            &,nprocs,ipt, stck_cnt)	   
!Carbon
	   case("ALGDOC") !!===============ALgae DOC===============================!
       i1 = lbound(ALGDOC,1) ; i2 = ubound(ALGDOC,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",AlgDOC,nc_ofid,ALGDOC_vid,myid&
            &,nprocs,ipt, stck_cnt)	

	   case("ALGPOC") !!===============ALgae POC===============================!
       i1 = lbound(ALGPOC,1) ; i2 = ubound(ALGPOc,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",AlgPOC,nc_ofid,ALGPOC_vid,myid&
            &,nprocs,ipt, stck_cnt)		   

	   case("HDRLPOC") !!===============Hydro. LPOC===============================!
       i1 = lbound(HDRLPOC,1) ; i2 = ubound(HDRLPOC,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",HDRLPOC,nc_ofid,HDRLPOC_vid,myid&
            &,nprocs,ipt, stck_cnt)	   

	   case("HDRRPOC") !!===============Hydro. RPOC===============================!
       i1 = lbound(HDRRPOC,1) ; i2 = ubound(HDRRPOC,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",HDRRPOC,nc_ofid,HDRRPOC_vid,myid&
            &,nprocs,ipt, stck_cnt)	   

	   case("COAGC") !!===============Coag DOC===============================!
       i1 = lbound(COAGC,1) ; i2 = ubound(COAGC,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",COAGC,nc_ofid,COAGC_vid,myid&
            &,nprocs,ipt, stck_cnt)	
			
	  case("DENIT") !!===============Denit. loss of DOC===============================!
       i1 = lbound(DENIT,1) ; i2 = ubound(DENIT,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",DENIT,nc_ofid,DENIT_vid,myid&
            &,nprocs,ipt, stck_cnt)				

	   case("MNLDOC1") !!===============Remin DOC1===============================!
       i1 = lbound(MNLDOC1,1) ; i2 = ubound(MNLDOC1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",MNLDOC1,nc_ofid,MNLDOC1_vid,myid&
            &,nprocs,ipt, stck_cnt)		   

	   case("MNLDOC2") !!===============Remin DOC2===============================!
       i1 = lbound(MNLDOC2,1) ; i2 = ubound(MNLDOC2,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",MNLDOC2,nc_ofid,MNLDOC2_vid,myid&
            &,nprocs,ipt, stck_cnt)	   
	   
	   case("MNLDOC3") !!===============Remin DOC3===============================!
       i1 = lbound(MNLDOC3,1) ; i2 = ubound(MNLDOC3,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",MNLDOC3,nc_ofid,MNLDOC3_vid,myid&
            &,nprocs,ipt, stck_cnt)	   
	   
	   case("PD31") !!===============Photochem 3 to 1===============================!
       i1 = lbound(PHOTOCHEM_DOC_ACCUM(:,:,2),1) ; i2 = ubound(PHOTOCHEM_DOC_ACCUM(:,:,2),1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PHOTOCHEM_DOC_ACCUM(:,:,2),nc_ofid,PD31_vid,myid&
            &,nprocs,ipt, stck_cnt)		   

	   case("PD30") !!===============Photochem 3 to 0===============================!
       i1 = lbound(PHOTOCHEM_DOC_ACCUM(:,:,3),1) ; i2 = ubound(PHOTOCHEM_DOC_ACCUM(:,:,3),1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PHOTOCHEM_DOC_ACCUM(:,:,3),nc_ofid,PD30_vid,myid&
            &,nprocs,ipt, stck_cnt)	   

	   case("PD3N") !!===============Photochem 3 to N===============================!
       i1 = lbound(PHOTOCHEM_DOC_ACCUM(:,:,4),1) ; i2 = ubound(PHOTOCHEM_DOC_ACCUM(:,:,4),1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PHOTOCHEM_DOC_ACCUM(:,:,4),nc_ofid,PD3N_vid,myid&
            &,nprocs,ipt, stck_cnt)	   

	   case("PD21") !!===============Photochem 2 to 1===============================!
       i1 = lbound(PHOTOCHEM_DOC_ACCUM(:,:,5),1) ; i2 = ubound(PHOTOCHEM_DOC_ACCUM(:,:,5),1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PHOTOCHEM_DOC_ACCUM(:,:,5),nc_ofid,PD21_vid,myid&
            &,nprocs,ipt, stck_cnt)		   

	   case("PD20") !!===============Photochem 2 to 0===============================!
       i1 = lbound(PHOTOCHEM_DOC_ACCUM(:,:,6),1) ; i2 = ubound(PHOTOCHEM_DOC_ACCUM(:,:,6),1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PHOTOCHEM_DOC_ACCUM(:,:,6),nc_ofid,PD20_vid,myid&
            &,nprocs,ipt, stck_cnt)	   

	   case("PD2N") !!===============Photochem 2 to N===============================!
       i1 = lbound(PHOTOCHEM_DOC_ACCUM(:,:,7),1) ; i2 = ubound(PHOTOCHEM_DOC_ACCUM(:,:,7),1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PHOTOCHEM_DOC_ACCUM(:,:,7),nc_ofid,PD2N_vid,myid&
            &,nprocs,ipt, stck_cnt)		   
	   
	   case("PD10") !!===============Photochem 1 to 0===============================!
       i1 = lbound(PHOTOCHEM_DOC_ACCUM(:,:,8),1) ; i2 = ubound(PHOTOCHEM_DOC_ACCUM(:,:,8),1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PHOTOCHEM_DOC_ACCUM(:,:,8),nc_ofid,PD10_vid,myid&
            &,nprocs,ipt, stck_cnt)	   

	   case("PD1N") !!===============Photochem 1 to N===============================!
       i1 = lbound(PHOTOCHEM_DOC_ACCUM(:,:,9),1) ; i2 = ubound(PHOTOCHEM_DOC_ACCUM(:,:,9),1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",PHOTOCHEM_DOC_ACCUM(:,:,9),nc_ofid,PD2N_vid,myid&
            &,nprocs,ipt, stck_cnt)		   
	   
	   
	   
#endif
			
			
     case default
       if(msr)then
         write(ipt,*)'variable',cdf_vdp(i),' not set up for netcdf output'
         write(ipt,*)'modify module MOD_NCDIO.f' 
         call pstop
       end if
     end select

   end do



!==============================================================================|
!  CLOSE THE FILE                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)
   return
 end subroutine out_netcdf
 
 
 !==============================================================================|
 !  Collect Data to Global Array and Write to Netcdf File                       |
 !==============================================================================|
 
 
 ! David added interface to putvar to select the right subroutine for
 ! the data type.                                       
 
 SUBROUTINE PUTVAR1D_REAL(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid&
      &,myid,nprocs,ipt,stk)
   
   !------------------------------------------------------------------------------|
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs&
        &,ipt,stk
   character(len=*),intent(in)   :: map_type
   real(sp), dimension(i1:i2) :: var
   
   real(sp), allocatable, dimension(:,:) :: temp

   allocate(temp(i1:i2,kt))
   temp(i1:i2,1)=var

   CALL PUTVAR2D_REAL(i1,i2,n1,n1gl,kt,k1,map_type,temp,nc_fid,vid&
        &,myid,nprocs,ipt,stk)
   
   deallocate(temp)

 END SUBROUTINE PUTVAR1D_REAL
   

 subroutine PUTVAR2D_REAL(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid&
      &,myid,nprocs,ipt,stk)
!------------------------------------------------------------------------------|

#  if defined (MULTIPROCESSOR)
   use mod_par
#  endif
   use mod_types
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs&
        &,ipt, stk
   character(len=*),intent(in)   :: map_type
   real(sp), dimension(i1:i2,kt) :: var

   real(sp), allocatable, dimension(:,:) :: temp,gtemp
   integer :: ierr,k1m1
   integer, allocatable :: dims(:)
   integer :: i,k  

   k1m1 = k1 
   if(k1m1 == 1)then
     allocate(dims(2))
     dims(1) = 1 
     dims(2) = stk
   else
     allocate(dims(3))
     dims(1) = 1 
     dims(2) = 1 
     dims(3) = stk      
   end if
     

   if(map_type(1:1) /= "e" .and. map_type(1:1) /= "n")then
     write(ipt,*)'map_type input to putvar should be "e" OR "n"'
     call pstop
   end if

   if(nprocs==1)then
     allocate(temp(n1,k1m1))  ; temp(1:n1,1:k1m1) = var(1:n1,1:k1m1)
   end if

#  if defined (MULTIPROCESSOR)
   if(nprocs > 1)then
     allocate(gtemp(n1gl,kt))
     if(map_type(1:1) == "e")then  !collect all data from 1:kt, but only print 1:k1m1 (kt>=k1m1==k1)
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,emap,var,gtemp)
     else 
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,nmap,var,gtemp)
     end if
     if(myid==1)then
        allocate(temp(n1gl,k1m1))  
        temp(1:n1gl,1:k1m1)=gtemp(1:n1gl,1:k1m1)
     endif
   end if
#  endif

!   if(myid /= 1) return
   if(myid==1) then
     !write(*,*)'Hey I am here 1' ! Written by Wen Long  !double precision , parallel, restart has problem here
     ierr = nf90_put_var(nc_fid,vid,temp,START=dims)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"error writing variable to netcdf 1",ipt)
     end if
     deallocate(temp)
   end if  
   deallocate(dims)

   return
 end subroutine PUTVAR2D_REAL
!==============================================================================|

 SUBROUTINE PUTVAR1D_INT(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid&
      &,myid,nprocs,ipt,stk )
   
   !------------------------------------------------------------------------------|
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs&
        &,ipt, stk
   character(len=*),intent(in)   :: map_type
   INTEGER, dimension(i1:i2) :: var
   
   INTEGER, allocatable, dimension(:,:) :: temp
   
   allocate(temp(i1:i2,kt))
   temp(i1:i2,kt)= var
   
   call PUTVAR2D_INT(i1,i2,n1,n1gl,kt,k1,map_type,temp,nc_fid,vid&
        &,myid,nprocs,ipt, stk)
   
   deallocate(temp)
   
 END SUBROUTINE PUTVAR1D_INT
 
 subroutine PUTVAR2D_INT(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid&
      &,myid,nprocs,ipt, stk)
   
   !------------------------------------------------------------------------------|
   
#  if defined (MULTIPROCESSOR)
   use mod_par
#  endif
!   use mod_types
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs&
        &,ipt,stk
   character(len=*),intent(in)   :: map_type
   INTEGER, dimension(i1:i2,kt) :: var

   INTEGER, allocatable, dimension(:,:) :: temp,gtemp
   integer :: ierr,k1m1
   integer, allocatable :: dims(:)
   

   k1m1 = k1 
   if(k1m1 == 1)then
     allocate(dims(2))
     dims(1) = 1 
     dims(2) = stk
   else
     allocate(dims(3))
     dims(1) = 1 
     dims(2) = 1 
     dims(3) = stk 
   end if
     

   if(map_type(1:1) /= "e" .and. map_type(1:1) /= "n")then
     write(ipt,*)'map_type input to putvar should be "e" OR "n"'
     call pstop
   end if

   if(nprocs==1)then
     allocate(temp(n1,k1m1))  ; temp(1:n1,1:k1m1) = var(1:n1,1:k1m1)
   end if

#  if defined (MULTIPROCESSOR)
   if(nprocs > 1)then
     allocate(gtemp(n1gl,kt))
     if(map_type(1:1) == "e")then
       call igather(i1,i2,n1,n1gl,kt,myid,nprocs,emap,var,gtemp)
     else 
       call igather(i1,i2,n1,n1gl,kt,myid,nprocs,nmap,var,gtemp)
     end if
     IF(myid==1)THEN
          allocate(temp(n1gl,k1m1))  ; temp(1:n1gl,1:k1m1) = gtemp(1:n1gl,1:k1m1)
     ENDIF
     deallocate(gtemp)
   end if
#  endif

!   if(myid /= 1) return
   if(myid==1) then
     ierr = nf90_put_var(nc_fid,vid,temp,START=dims)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"error writing variable to netcdf 2",ipt)
     end if
     deallocate(temp)
   endif
   deallocate(dims)
   return
 end subroutine PUTVAR2D_INT
!==============================================================================|
!Subroutines :
!	Subroutine NCD_READ_GRID()
!	Subroutine NCD_READ_SHAPE()
!	Subroutine NCD_READ_OPEN()
!	Subroutine NCD_READ()

SUBROUTINE NCD_READ_GRID(INFILE)
  !---------------------------------------------------------------------
  ! READ DIMENSIONS IN A NETCDF FILES
  !---------------------------------------------------------------------
  USE MOD_NCD, ONLY : NC_FID, GETDIM, GETSVAR, GETDVAR
  USE NETCDF
  USE MOD_LIMS, ONLY : MLOC, NLOC, NTLOC, KBM2, KB, KBM1, MTLOC
  USE MOD_SIZES, ONLY: NGL, MGL
  USE MOD_PREC, ONLY : SP
  	!Wen Long took CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY : 		&
			!SERIAL  		&           !!TRUE IF SINGLE PROCESSOR
			MSR  !      	&           !!TRUE IF MASTER PROCESSOR (MYID==1)
			!,PAR        !	&           !!TRUE IF MULTIPROCESSOR RUN
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

  	
  IMPLICIT NONE
  !----------------------------------------------------------------------------!
  CHARACTER(LEN=1024), INTENT(IN) :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: N_ELEMS,N_NODES,N_SIG_M1,N_SIG
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP
  !----------------------------------------------------------------------------!

  !--Open NetCDF DATA FILE
  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
  if(MSR)WRITE(*,*)'opening netcdf file, INFILE, NC_FID, NGL '//TRIM(INFILE),NC_FID,NGL
  IF(IERR /=NF90_NOERR)THEN
     WRITE(*,*)'ERROR READING ',TRIM(INFILE)
     WRITE(*,*)TRIM(NF90_STRERROR(IERR))
     CALL PSTOP
  ENDIF

  !--Get Model Dimensions
  N_ELEMS   = GETDIM(NC_FID,LEN_TRIM('nele'),'nele')
  N_NODES   = GETDIM(NC_FID,LEN_TRIM('node'),'node')
  N_SIG_M1  = GETDIM(NC_FID,LEN_TRIM('siglay'),'siglay')

  N_SIG     = GETDIM(NC_FID,LEN_TRIM('siglev'),'siglev')
  IF(N_ELEMS /= NGL)THEN
    IF(MSR)THEN 
      WRITE(*,*) 'NGL is inconsistent with the # of ELEMENTS in NetCDF file'
      WRITE(*,*) 'NGL should be equal to',N_ELEMS,NGL
      CALL PSTOP
    ENDIF
  ENDIF    
  IF(N_NODES /= MGL)THEN
    IF(MSR)THEN 
      WRITE(*,*) 'MGL IS inconsistent with the # of NODES in NetCDF file'
      WRITE(*,*) 'MGL should be equal to',N_NODES
      CALL PSTOP
    ENDIF
  ENDIF    
  IF(N_SIG /= KB)THEN
    IF(MSR)THEN 
      WRITE(*,*) 'KB IS inconsistent with the # of SIGMA LAYERS in NetCDF file'
      WRITE(*,*) 'KB should be equal to',N_SIG
      CALL PSTOP
    ENDIF
  ENDIF    
  
  MLOC=N_NODES
  NLOC=N_ELEMS
  MTLOC = MLOC
  NTLOC = NLOC
  
  KB=N_SIG    
  KBM1=N_SIG_M1
  KBM2=KB-2

  !--close file
  IERR = NF90_CLOSE(NC_FID)

  RETURN
END SUBROUTINE NCD_READ_GRID

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE NCD_READ_SHAPE(INFILE)
  !---------------------------------------------------------------------
  ! READ BATHYMETRY, SIGMA LEVELS AND GRID COEFFICIENTS IN A NETCDF FILES
  !---------------------------------------------------------------------
  USE MOD_NCD, ONLY : NC_FID, GETSVAR
  USE NETCDF
  USE MOD_LIMS, ONLY: KBM2
  USE MOD_HYDROVARS, ONLY: &
   		!GRAV	&		!
		!,PI	&		!
		!,PI2	&		!
		!,ZERO	&		!
		!,ONE_THIRD	&	!
		NVG	&		!
		,XG	&		!GLOBAL X-COORD AT NODE 
		,YG	&		!GLOBAL X-COORD AT NODE 
		,HG	&		!GLOBAL DEPTH AT NODE 
		!,XCG	&		!GLOBAL X-COORD AT FACE CENTER 
		!,YCG	&		!GLOBAL X-COORD AT FACE CENTER 
		!,VXMIN	&		!
		!,VYMIN	&		!
		!,VXMAX	&		!
		!,VYMAX	&		!
		!,XC	&		!X-COORD AT FACE CENTER 
		!,YC	&		!Y-COORD AT FACE CENTER
		!,VX	&		!X-COORD AT GRID POINT
		!,VY	&		!Y-COORD AT GRID POINT
		!,ART	&		!AREA OF ELEMENT
		!,ART1	&		!AREA OF NODE-BASE CONTROl VOLUME
		!,ART2	&		!AREA OF ELEMENTS AROUND NODE
		!,NV	&		!NODE NUMBERING FOR ELEMENTS
		!,NBE	&		!INDICES OF ELMNT NEIGHBORS
		!,NTVE	&		!
		!,NTSN	&		!
		!,ISONB	&		!NODE MARKER = 0,1,2 
		!,ISBC	&		!
		!,ISBCE	&		!
		!,IEC	&		!
		!,IENODE &		!
		!,NBSN	&		!
		!,NIEC	&		!
		!,NTRG	&		!
		!,NBVE	&		!
		!,NBVT	&		!
		!,LISBCE_1	&	!LIST OF ELEMENTS WITH ISBCE=1
		!,LISBCE_2	&	!LIST OF ELEMENTS WITH ISBCE=2
		!,LISBCE_3	&	!LIST OF ELEMENTS WITH ISBCE=3
		!,DLTXC	&		!
		!,DLTYC	&		!
		!,DLTXYC	&	!
		!,DLTXE	&		!
		!,DLTYE	&		!
		!,DLTXYE	&	!
		!,SITAC	&		!
		!,SITAE	&		!
		!,XIJC	&		!
		!,YIJC	&		!
		!,XIJE	&		!
		!,YIJE	&		!
		!,EPOR	&		!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
		!,IBCGEO	&	!LOCAL GEOSTROPHIC FRICTION CORRECTION NODES

		,Z	&			!SIGMA COORDINATE VALUE 
		,ZZ	&		!INTRA LEVEL SIGMA VALUE
		,DZ	&		!DELTA-SIGMA VALUE
		,DZZ	&		!DELTA OF INTRA LEVEL SIGMA 
		!,H1	&		!BATHYMETRIC DEPTH 
		!,H	&			!BATHYMETRIC DEPTH 
		!,D	&			!CURRENT DEPTH 
		!,DT	&		!DEPTH AT PREVIOUS TIME STEP
		!,DT1	&		!DEPTH AT PREVIOUS TIME STEP
		!,EL	&		!CURRENT SURFACE ELEVATION
		!,ET	&		!SURFACE ELEVATION AT PREVIOUS TIME STEP
		!,DTFA	&		!ADJUSTED DEPTH FOR MASS CONSERVATION
		!,UU	&		!X-VELOCITY
		!,VV	&		!Y-VELOCITY
		!,UUT	&		!X-VELOCITY FROM PREVIOUS TIMESTEP
		!,VVT	&		!Y-VELOCITY FROM PREVIOUS TIMESTEP
		!,WWT	&		!Z-VELOCITY FROM PREVIOUS TIMESTEP
		!,WTST	&		!Vertical velocity in sigma from PREVIOUS TIMESTEP
		!,UARD_OBCNT	&!tykim
		!,XFLUX_OBCT	&!tykim
		!,DTFAT	&		!tykim
		!,TT_T	&		!tykim
		!,SALTT	&		!tykim
		!,WTS	&		!VERTICAL VELOCITY IN SIGMA SYSTEM
		!,UARD_OBCN	&	! tykim 
		!,XFLUX_OBC	&	! tykim 
		!,WTTS	&		!VERTICAL VELOCITY IN SIGMA SYSTEM 
		!,KH	&		!TURBULENT DIFFUSIVITY
		!,A1U	&		!
		!,A2U	&		!
		!,AWX	&		!
		!,AWY	&		!
		!,AW0	&		!
		!,VISCOFH	&	!
		!,UNC1	&		!
		!,VNC1	&		!
		!,WNC1	&		!
		!,WTSNC1	&		!
		!,UARD_OBCNNC1	&	!
		!,XFLUX_OBCNC1	&	!
		!,DTFANC1	&		!
		!,KHNC1	&		!
		!,TNC1	&		!
		!,SNC1	&		!
		!,ELNC1	&		!
		!,UNC2	&		!
		!,VNC2	&		!
		!,WNC2	&		!
		!,WTSNC2	&	!
		!,UARD_OBCNNC2	&!
		!,XFLUX_OBCNC2	&!
		!,DTFANC2	&	!
		!,KHNC2	&		!
		!,TNC2	&		!
		!,SNC2	&		!
		!,ELNC2	&		!
		,num_hyd_ints  !	&!number of records in each hydrodynamics netcdf file
		!,TIME_MAP	&	!
		!,THOUR1	&	!SIMULATION TIME AT END OF CURRENT EXTERNAL STEP (IEXT) IN HOURS
		!,THOUR	&		!
		!,NCFILE_DIR	&!
		!,NCFILE_PREFIX	&!
		!,NCFILE_SUFFIX	&!
		!,NCFILE_NUMBER	&!
		!,FORMAT_STR	&!
		!,hydro_dir, 	&	! directory name where hydrodynamics results (netcdf) files are stored
		!,hydro_prefix, &	! prefix of file name, e.g. 'psm_'
		!,hydro_suffix	&	! suffix of filename, e.g. '.nc'
		!,hydro_filenumwidth, &	! number of digits in filename following hydro_prefix, e.g. 4 for psm_0002.nc
		!,hydro_filenumstart, &	! starting number of the file name in the digital part of the file name, e.g. 185 for psm_0185.nc
		!,hydro_Nrec	&		! number of records in each of hydrodynamics file
		!,hydro_dlt	&			! time step in hydrodynamics file (in seconds), e.g. 100 for 100sec
		!,t_his_start	&		!
		!,t_his_end	&			!
		!,t_his_dlt	&			!starting time, ending time, and interval of history outputs (days)
		!,Nstation	&			!
		!,NstationNum_GL	&	!maximum number of station is NstationMax!
		!,t_stn_start	&		!
		!,t_stn_end	&			!
		!,t_stn_dlt	&			!starting time, ending time, and interval of station outputs (days)
		!,STNFN	&				!file name for station output
		!,HISFN	&				!file name for history output
		!,HISFN_PREFIX	&		!prefix of history output file
		!,HISFN_EXT	&			!extention name of history output file
		!,HISFN_FINAL	&		! 
		!,HISFN_SPLIT_BYLEVEL	&!True or False for splitting history output in files level by level (default is .FALSE.)
		!,hydro_netcdf	&		!
		!,wqm_history	&		!
		!,wqm_stations	&		!
		!,IFNC	&				!file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise 
		!,NTRECNC	&			!time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file. 
		!,NTHYDRO				!overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read

	!Wen Long took CONTROL out of MOD_HYDROVARS and put the used variables here
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
			
  USE MOD_LIMS, ONLY :  KB, KBM1
  USE MOD_SIZES, ONLY: NGL, MGL
  USE MOD_PREC, ONLY : SP
 
  IMPLICIT NONE
  !----------------------------------------------------------------------------!
  CHARACTER(LEN=1024), INTENT(IN) :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: I,K
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP

  !--OPEN NETCDF DATA FILE
  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
  if(MSR)WRITE(*,*)'opening netcdf file, INFILE, NC_FID'//TRIM(INFILE),NC_FID
  IF(IERR /= NF90_NOERR)THEN
    WRITE(*,*)'ERROR READING ',TRIM(INFILE)
    WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    STOP
  ENDIF

  !--Get Node Coordinates
  !WLong moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F 
  !ALLOCATE(XG(0:MGL),YG(0:MGL)) ; XG = 0.0_SP ; YG = 0.0_SP
  
  ALLOCATE(TEMP(MGL,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('x'),'x',MGL,1,TEMP)
  XG(1:MGL) = TEMP(1:MGL,1)
  DEALLOCATE(TEMP)

  ALLOCATE(TEMP(MGL,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('y'),'y',MGL,1,TEMP)
  YG(1:MGL) = TEMP(1:MGL,1)
  DEALLOCATE(TEMP)

  !--Get Node Numbering
  !WLong moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F 
  !ALLOCATE(NVG(0:NGL,4)); NVG = 0

  ALLOCATE(TEMP(NGL,3))
  CALL GETSVAR(NC_FID,LEN_TRIM('nv'),'nv',NGL,3,TEMP)
  NVG(1:NGL,1:3) = TEMP(1:NGL,1:3)
  DEALLOCATE(TEMP)
  NVG(:,4) = NVG(:,1)

  !--Get Bathymetry
  !!WLong moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F 
  !ALLOCATE(HG(0:MGL))  ; HG = 0.0_SP

  ALLOCATE(TEMP(MGL,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('h'),'h',MGL,1,TEMP)
  HG(1:MGL) = TEMP(1:MGL,1)
  DEALLOCATE(TEMP)

 !WLong moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F 
 ! ALLOCATE(Z(KB))  ; Z   = 0.0_SP    !!SIGMA COORDINATE VALUE 
 ! ALLOCATE(ZZ(KB)) ; ZZ  = 0.0_SP    !!INTRA LEVEL SIGMA VALUE
 ! ALLOCATE(DZ(KB)) ; DZ  = 0.0_SP    !!DELTA-SIGMA VALUE
 ! ALLOCATE(DZZ(KB)); DZZ = 0.0_SP    !!DELTA OF INTRA LEVEL SIGMA 

  !--Get Sigma levels
  ALLOCATE(TEMP(KB,1))
  CALL GETSVAR(NC_FID,LEN_TRIM('siglev'),'siglev',KB,1,TEMP)
  Z(1:KB) = TEMP(1:KB,1)
  ! 	 DO I =1,MGL
!	    DO K = 1,KB
!	    SIGLEV(I,K)=Z(K)
!		END DO
!		DO K=1,KBM1
!		siglay(I,K)=ZZ(K)
!		ENDDO
  !   END DO
!	 write(*,*) 'Z = ',Z
  !write(*,*)'Siglev = ', SIGLEV
  !write(*,*)'Siglay = ',SIGLAY
  DEALLOCATE(TEMP)
!      do K=1,KB
!!tykim00
!!         Z(K)=0.0-(K-1)*(1.0/(KB*1.0-1.0))

!!
!!WLong: this vertical sigma coordinate formulation is hardired
!!       we should have this information in the netcdf file itself
!!       or we should get this **1.5 as a parameter from input control
!!
!        Z(K) = -((K-1)/FLOAT(KB-1))**1.5  
!               
!      enddo
  !--Compute derivative and intra-sigma levels
  DO K=1,KBM1
    ZZ(K)=0.5_SP*(Z(K)+Z(K+1))
    DZ(K)=Z(K)-Z(K+1)
	!SIGLAY(K,:)= 0.5_SP*(SIGLEV(K,:)+SIGLEV(K+1,:))
  ENDDO
  ZZ(KB)=2.0_SP*ZZ(KBM1)-ZZ(KBM2)
  !SIGLAY(KB,:)=2.0_SP*SIGLAY(KBM1,:)-SIGLAY(KBM2,:)
  
  DO K=1,KBM2
    DZZ(K)=ZZ(K)-ZZ(K+1)
  ENDDO
  DZZ(KB-1)=0.0
  DZ(KB)=0.0

  !--Close file
  IERR = NF90_CLOSE(NC_FID)

  RETURN
END SUBROUTINE NCD_READ_SHAPE

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
! 15 SEPT 2009, KURT GLAESEMANN, REMOVE ALL THE CODE FROM THIS ROUTINE AND CALL DUPLICATE CODE IN NCD_READ
SUBROUTINE NCD_READ_OPEN(INFILE,UL,VL,WTSL,UARD_OBCNL,XFLUX_OBCL,DTFAL,KHL,ELL,TL,SL,HO)
  !---------------------------------------------------------------------
  ! READ DATA FROM DAILY NETCDF FILES
  !---------------------------------------------------------------------

  USE MOD_NCD, ONLY : NF90_OPEN,NF90_NOWRITE,NF90_STRERROR,NF90_NOERR,NC_FID
  USE MOD_SIZES, ONLY : NOBTY
  USE MOD_LIMS, ONLY: MTLOC, NTLOC, KBM1, KB
  USE MOD_PREC, ONLY: SP
  
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
  IMPLICIT NONE
  !----------------------------------------------------------------------------!

  REAL(SP), DIMENSION(0:NTLOC,KB),INTENT(OUT)   :: UL,VL!,WWL     !,WWL
  REAL(SP), DIMENSION(0:MTLOC,KB),INTENT(OUT)   :: KHL,WTSL
  REAL(SP), DIMENSION(0:MTLOC,KBM1),INTENT(OUT) :: TL,SL
  REAL(SP), DIMENSION(0:MTLOC),INTENT(OUT)      :: ELL,DTFAL
  REAL(SP), DIMENSION(0:NOBTY+1),INTENT(OUT)    :: UARD_OBCNL
  REAL(SP), DIMENSION(0:NOBTY,KBM1),INTENT(OUT) :: XFLUX_OBCL
  
  INTEGER, INTENT(IN)                       :: HO
  CHARACTER(LEN=1024), INTENT(IN)           :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
!  CHARACTER(LEN=300) TASK

  !--Open NetCDF Datafile 
  IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
! KURT GLAESEMANN added MSR to print out
  if(MSR)WRITE(*,*)'opening netcdf file, INFILE, NC_FID'//TRIM(INFILE),NC_FID

  IF(IERR /=NF90_NOERR)THEN
     WRITE(*,*)'ERROR READING ',TRIM(INFILE)
     WRITE(*,*)TRIM(NF90_STRERROR(IERR))
     STOP
  ENDIF

  CALL NCD_READ(INFILE,UL,VL,WTSL,UARD_OBCNL,XFLUX_OBCL,DTFAL,KHL,ELL,TL,SL,HO)

  RETURN
END SUBROUTINE NCD_READ_OPEN


!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE NCD_READ(INFILE,UL,VL,WTSL,UARD_OBCNL,XFLUX_OBCL,DTFAL,KHL,ELL,TL,SL,HO)
  !---------------------------------------------------------------------
  ! READ DATA FROM DAILY NETCDF FILES
  !---------------------------------------------------------------------
  USE MOD_PREC, ONLY: SP
  USE MOD_NCD!, ONLY : NC_FID, GETDVAR
  USE NETCDF
  
  USE MOD_WQM, ONLY: 		&!
			DTFAL_GL,       &!
            ELL_GL,         &!
            KHL_GL,         &!
            SL_GL,          &!
            TL_GL,          &!
            UARD_OBCN_GL,   &!
            UL_GL,          &!
            VL_GL,          &!
			WW_GL,          &!
            WTSL_GL,        &!
            XFLUX_OBC_GL

  USE MOD_LIMS, ONLY: MTLOC, NTLOC, MLOC, NLOC, KB,KBM1
  USE MOD_SIZES, ONLY: MGL,NGL, NOBTY
  
  
  
  USE MOD_HYDROVARS, ONLY : NUM_HYD_INTS
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
		
  USE MOD_BCMAP, ONLY: 	&
		IOBCN,			&
		I_OBC_GL		

   
# if defined (MULTIPROCESSOR)
		USE MOD_PAR, ONLY:   &	  
			!EL_PID,	&		!(:)       	!!PROCESSOR OWNER OF GLOBAL ELEMENT
			!ELID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT
			 NLID,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE 
			!ELID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
			!NLID_X,	&		!(:)       	!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
			!OBN_GL2LOC,	&	!(:)       	!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
			 EGID,		&	!(:)       	!!GLOBAL ID OF LOCAL ELEMENT
			 NGID,		&	!(:)       	!!GLOBAL ID OF LOCAL NODE 
			 NHE,        &   !			!NUMBER OF HALO ELEMENTS
			 HE_LST,		&	!(:)       	!!GLOBAL IDENTITIES OF HALO ELEMENTS
			!HE_OWN,		&	!(:)       	!!OWNER OF HALO ELEMENTS
			 NBN ,       &   !          	!!NUMBER OF BOUNDARY NODES
			!MX_MLT,     &   !          	!!MAX MULTIPLICITY OF BOUNDARY NODES
			!BN_LST, 	&	!	   		!!GLOBAL IDENTITY OF BOUNDARY NODES
			 BN_LOC, 	&	!			!!LOCAL IDENTITY OF BOUNDARY NODES
			 BN_MLT, 	&	!			!!MULTIPLICITY OF BOUNDARY NODES
			!BN_NEY, 	&	!			!!NODE OWNER LIST
			!NDE_ID, 	&	!			!! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
			 NHN,  		&	!			!!NUMBER OF HALO NODES
			 HN_LST, 	&	!			!!LIST OF HALO NODES 
			!HN_OWN, 	&	!			!!PRIMARY OWNER OF HALO NODES
			!EC,			&	!
			!NC,			&	!
			 BNC!, 		&	!
			!EMAP,		&	!
			!NMAP, 		&	!
			!PNE,		&	!(:)        !!NUMBER OF ELEMENTS IN EACH PROC
			!PNN,		&	!(:)        !!NUMBER OF NODES IN EACH PROC
			!PNHE,		&	!(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
			!PNBN,		&	!(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
			!PMBM,		&	!(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
			!PNHN,		&	!(:)       !!NUMBER OF HALO NODES IN EACH PROC
			!FILELOCK,   &
			!MYFILELOCK, &
			!FILELOCKED , &
			!NODE_MATCH, &
			!NODE_MATCH_ALL, &
			!EXCHANGE, &
			!EXCHANGE_ALL, &
			!GATHER, &
			!IGATHER, &
			!INIT_MPI_ENV,&!			 
			!SORT,	&
			!PREAD,	&
			!PWRITE,&
			!PPRINT, &
			!LOCK, &
			!UNLOCK, &
			!GETLOC 

# endif  
  IMPLICIT NONE

!----------------------------------------------------------------------------!
!
!  REAL(CDF_PREC), DIMENSION(NGL,KBM1)   :: UL_GL,VL_GL                 !u and v
!  REAL(CDF_PREC), DIMENSION(MGL,KB)     :: WTSL_GL,KHL_GL              !wts, kh
!  REAL(CDF_PREC), DIMENSION(MGL,KBM1)   :: SL_GL,TL_GL                 !salinity,temp
!  REAL(CDF_PREC), DIMENSION(MGL)        :: ELL_GL,DTFAL_GL             !zeta,dtfa
!  REAL(CDF_PREC), DIMENSION(NOBTY+1)    :: UARD_OBCN_GL                !uard_obcn
!  REAL(CDF_PREC), DIMENSION(NOBTY,KBM1) :: XFLUX_OBC_GL                !xflux_obc
!

  REAL(SP), DIMENSION(0:NTLOC,KB),INTENT(OUT)   :: UL,VL!, WWL
  REAL(SP), DIMENSION(0:MTLOC,KB),INTENT(OUT)   :: KHL,WTSL
  REAL(SP), DIMENSION(0:MTLOC,KBM1),INTENT(OUT) :: TL,SL
  REAL(SP), DIMENSION(0:MTLOC),INTENT(OUT)      :: ELL,DTFAL
  REAL(SP), DIMENSION(0:NOBTY+1),INTENT(OUT)    :: UARD_OBCNL
  REAL(SP), DIMENSION(0:NOBTY,KBM1),INTENT(OUT) :: XFLUX_OBCL
  INTEGER, INTENT(IN)                       :: HO
  CHARACTER(LEN=1024), INTENT(IN)           :: INFILE
  !----------------------------------------------------------------------------!
  INTEGER            :: IERR
  INTEGER            :: HT
  INTEGER            :: I,K,J

!--Adjustement to read in Netcdf file
!  HT=HO+1
   HT=HO 

  !---------------------------------------------------------------------
  ! Read Data from file INFILE at time level ht
  !---------------------------------------------------------------------

! KURT GLAESEMANN SEPT 24 2009 - REORDER DATA READS TO MATCH ORDER IN FILE - FASTER
  !--U velocity 
  CALL GETDVAR(NC_FID,LEN_TRIM('u'),'u',NGL,KBM1,UL_GL,HT)

  !--V velocity
  CALL GETDVAR(NC_FID,LEN_TRIM('v'),'v',NGL,KBM1,VL_GL,HT)
  
  ! -- W Velocity at element
  CALL GETDVAR(NC_FID,LEN_TRIM('ww'),'ww',NGL,KBM1,WW_GL,HT)

!--Wen Long debugging----
!  IF(MSR)THEN
!     WRITE(*,'(A32,1x,F20.10,1x,F20.10,1x,F20.10)')'Wen Long (ELTMS,U,V) Surface == ',ELTMS,UL_GL(565,1),VL_GL(565,1)
!  ENDIF
!------------------------

  !--WTS velocity
  CALL GETDVAR(NC_FID,LEN_TRIM('wts'),'wts',MGL,KB,WTSL_GL,HT)

  !--UARD_OBCN
  CALL GETDVAR(NC_FID,LEN_TRIM('uard_obcn'),'uard_obcn',NOBTY,1,UARD_OBCN_GL,HT)

  !--XFLUX_OBCN
  CALL GETDVAR(NC_FID,LEN_TRIM('xflux_obc'),'xflux_obc',NOBTY,KBM1,XFLUX_OBC_GL,HT)

  !--free surface elevation
  CALL GETDVAR(NC_FID,LEN_TRIM('dtfa'),'dtfa',MGL,1,DTFAL_GL,HT)

! KURT GLAESEMAN 23 SEPT 2009 - fix dimensions of KH (KB, not KBM1)- here and above and below
  !--KH
  CALL GETDVAR(NC_FID,LEN_TRIM('kh'),'kh',MGL,KB,KHL_GL,HT)

  !--free surface elevation
  CALL GETDVAR(NC_FID,LEN_TRIM('zeta'),'zeta',MGL,1,ELL_GL,HT)

  !--salinity
  CALL GETDVAR(NC_FID,LEN_TRIM('salinity'),'salinity',MGL,KBM1,SL_GL,HT)

  !--temperature
  CALL GETDVAR(NC_FID,LEN_TRIM('temp'),'temp',MGL,KBM1,TL_GL,HT)
  
  ! B Clark add sanity check for bad temperatures to correct to reasonable value;
    DO I = 1,MGL
	  DO K = 1,KBM1
	    IF (TL_GL(I,K) < -5.0) THEN
	        TL_GL(I,K) = -5.0
			write(*,*)'corrected temp at',MGL
			write(*,*)'back to a value of',TL_GL(I,K)
        ENDIF
      END DO
    END DO
! 15 SEPT KURT GLAESEMANN - MOVED ALL COPIES TOGETHER
  IF(SERIAL) THEN
     ELL(1:MGL) = ELL_GL(1:MGL)
     SL(1:MGL,1:KBM1) = SL_GL(1:MGL,1:KBM1)
     TL(1:MGL,1:KBM1) = TL_GL(1:MGL,1:KBM1)
     UL(1:NGL,1:KBM1) = UL_GL(1:NGL,1:KBM1)
     VL(1:NGL,1:KBM1) = VL_GL(1:NGL,1:KBM1)
     WTSL(1:MGL,1:KB) = WTSL_GL(1:MGL,1:KB)
     UARD_OBCNL(1:NOBTY) = UARD_OBCN_GL(1:NOBTY)
     XFLUX_OBCL(1:NOBTY,1:KBM1) = XFLUX_OBC_GL(1:NOBTY,1:KBM1)
     DTFAL(1:MGL) = DTFAL_GL(1:MGL)
     KHL(1:MGL,1:KB) = KHL_GL(1:MGL,1:KB)
  ENDIF

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    DO I=1,MLOC
	  !WRITE(*,*)'ELL_GL Debug here here 05 '
      ELL(I) = ELL_GL(NGID(I))
	  !WRITE(*,*)'ELL_GL Debug here here 5 '
    ENDDO

    DO I=1,NHN
	  !WRITE(*,*)'ELL_GL Debug here here 06 '
      ELL(I+MLOC) = ELL_GL(HN_LST(I))
	  !WRITE(*,*)'ELL_GL Debug here here 6 '
    ENDDO
!
    DO K=1,KBM1
      DO I=1,MLOC
        SL(I,K) = SL_GL(NGID(I),K)
      ENDDO

      DO I=1,NHN
        SL(I+MLOC,K) = SL_GL(HN_LST(I),K)
      ENDDO
    ENDDO
!
    DO K=1,KBM1
      DO I=1,MLOC
        TL(I,K) = TL_GL(NGID(I),K)
      ENDDO

      DO I=1,NHN
        TL(I+MLOC,K) = TL_GL(HN_LST(I),K)
      ENDDO
    ENDDO
!
    DO K=1,KBM1
      DO I=1,NLOC
        UL(I,K) = UL_GL(EGID(I),K)
      ENDDO

      DO I=1,NHE
        UL(I+NLOC,K) = UL_GL(HE_LST(I),K)
      ENDDO
    ENDDO
!
    DO K=1,KBM1
      DO I=1,NLOC
        VL(I,K) = VL_GL(EGID(I),K)
      ENDDO

      DO I=1,NHE
        VL(I+NLOC,K) = VL_GL(HE_LST(I),K)
      ENDDO
    ENDDO
	
	!DO K=1,KBM1
    !  DO I=1,NLOC
    !    WWL(I,K) = WW_GL(EGID(I),K)
    !  ENDDO

    !  DO I=1,NHE
    !    WWL(I+NLOC,K) = WW_GL(HE_LST(I),K)
    !  ENDDO
    !ENDDO
!
    DO K=1,KB
      DO I=1,MLOC
        WTSL(I,K) = WTSL_GL(NGID(I),K)
      ENDDO

      DO I=1,NHN
        WTSL(I+MLOC,K) = WTSL_GL(HN_LST(I),K)
      ENDDO
    ENDDO

       UARD_OBCNL = -99999
       K = 0
       DO I = 1, NOBTY
         J = NLID(I_OBC_GL(I))
         IF ( J .NE. 0) THEN
            K = K + 1
            UARD_OBCNL(K) = UARD_OBCN_GL(I)
         ENDIF
       ENDDO
       IF (K .NE. IOBCN) THEN
         write(*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
         call pstop
       endif
       XFLUX_OBCL = -99999
       K = 0
       DO I = 1, NOBTY
         J = NLID(I_OBC_GL(I))
         IF ( J .NE. 0) THEN
            K = K + 1
            XFLUX_OBCL(K,:) = XFLUX_OBC_GL(I,:)
         ENDIF
       ENDDO
       IF (K .NE. IOBCN) THEN
         write(*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
         call pstop
       endif
!
    DO I=1,MLOC
      DTFAL(I) = DTFAL_GL(NGID(I))
    ENDDO

    DO I=1,NHN
      DTFAL(I+MLOC) = DTFAL_GL(HN_LST(I))
    ENDDO
!
    DO K=1,KB
      DO I=1,MLOC
        KHL(I,K) = KHL_GL(NGID(I),K)
      ENDDO

      DO I=1,NHN
        KHL(I+MLOC,K) = KHL_GL(HN_LST(I),K)
      ENDDO
    ENDDO
  ENDIF
# endif

  IF(HT.eq.num_hyd_ints) THEN
   IERR = NF90_CLOSE(NC_FID)
  ENDIF 

  RETURN
END SUBROUTINE NCD_READ

!#  endif

   END MODULE mod_ncdio
