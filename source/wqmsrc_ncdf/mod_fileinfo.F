        MODULE MOD_FILEINFO
			 USE MOD_SIZES, ONLY :        &	!
			NCP !,            &  !
			!NQFP,           &  !
            !NHQP,           &  !
			!NS1P,           &  !
            !NS2P,           &  !
            !NS3P,           &  !
            !NBCP,           &  !
            !NDP,            &  !
            !NFLP,           &  !
            !NOIP,           &  !
            !NSSFP,			 &  !
            !MGL,            &  !
			!NGL,            &  !
			!OBCGL,          &  !
			!NOBTY  
  
           IMPLICIT NONE
           SAVE
           INTEGER :: & ! HYD, 				!Wen removed HYD
					DIA, 				&
					CBC,				&
					S1,					&
					S2,					&
					S3, 				&             
                    BFI, 				&   !Benthic Flux data file or Sediment Diagenesis module input file
					BAI,                &   !WLong benthic algae input file
					DFI,				&   !WLong deposition feeder input file
					WCL,				&	!WLong sediment overlying water column input file
					MET,				&
					BFO, 				&
					BAO,                &  !WLong Benthic algae output file
					DFO,				&  !WLong Deposition feeder output file
					KEI,				&
                    ATM,				&
					STL,				&
					AGR,				&
					SVI,				&
					SVO,				&
					KFL, 				&
					ZOO,				&
					ZFO,				&
					ALO,				&
                    CON, 				&
					RSO,				&
					SNP,				&
!					GEO,				&  !Wen Long removed GEO
					PLT,				&
					APL,				&
					TFL,				&
					OPL,				&
					SFI,				&
					SFO,				&
                   ! MAP, 				& !Wen Long removed MAP
					ICI,				&
					ICO,				&
					MRL,				&
					MBL,				&
					RSI,				&
                    UNIT_LINKAGE,		&
					UNIT_STN,			&
					UNIT_HIS                  

          !CHARACTER(LEN=20) :: SSNAME(17)  !Wen Long moved to wqm_sed.F
           CHARACTER(LEN=24) :: CNAME(NCP)  !Wen Long this should be moved to wqm_alg.F and also wqm_alg and wqm_kin should be
                                            !combined, otherwise, we would need to separate out zooplankton etc
                                            !
           
           CHARACTER(LEN=18) :: CONFN ='wqm_con.npt'

        !--------------------end of data declarations---------------------------------c

           CONTAINS

			!subroutines:
			! subroutine	INIT_FILE_INFO()
		   
           SUBROUTINE INIT_FILE_INFO()
        
           UNIT_LINKAGE=7   !file handle for LINKAGE input control 
           UNIT_STN    =200 !file handle for station output control
           UNIT_HIS    =100 !file handle for history output control
                            !actual handle will be UNIT_HIS+k, where k is layer number   
           CON = 10 
          ! HYD = 20	!Wen Long removed  HYD
           BFI = 21
		   BAI = 22		!WLong benthic alage input  
		   DFI = 32
            S1 = 23
            S2 = 24 
            S3 = 34
           ATM = 25 
           STL = 26 
           AGR = 27
           SVI = 28 
           MET = 29 
           KEI = 30 
           SVO = 31
           ZOO = 38
           ZFO = 39
           ALO = 50
           DIA = 40 
           BFO = 41 
		   BAO = 76 !WLong bentic algae output
		   DFO = 77 !WLong, Deposition feeder output
           KFL = 47
!           MAP = 32 
!           GEO = 33  
           ICI = 35 
           ICO = 36
           MRL = 37 
           SFI = 58
           SFO = 59
           RSO = 42 
           SNP = 43 
           PLT = 44 
           APL = 45 
           TFL = 46 
           OPL = 48
           MBL = 49 
           RSI = 75
           CBC = 60 
           
           !WLong The following SSNAME moved to wqm_sed.F
           !SSNAME(1)  = 'Sediment Temperature'
           !SSNAME(2)  = 'Sediment POP        '
           !SSNAME(3)  = 'Sediment PON        '  
           !SSNAME(4)  = 'Sediment POC        ' 
           !SSNAME(5)  = 'Sediment PBS        ' 
           !SSNAME(6)  = 'Sediment PO4        '
           !SSNAME(7)  = 'Sediment NH4        '
           !SSNAME(8)  = 'Sediment NO3        '
           !SSNAME(9)  = 'Sediment HS         '
           !SSNAME(10) = 'Sediment CH4        '
           !SSNAME(11) = 'Sediment CH4        '
           !SSNAME(12) = 'Sediment SO4        '
           !SSNAME(13) = 'Sediment DSIL       '
           !SSNAME(14) = 'Benthic Stress      '
           !SSNAME(15) = 'Benthic Algae       '
           !SSNAME(16) = 'Deposit Feeders     '
           !SSNAME(17) = 'Suspension Feeders  '

           CNAME(1)  =  'Temperature'
           CNAME(2)  =  'Salinity'
           CNAME(3)  =  'Inorganic Solids'
           CNAME(4)  =  'Algal Group 1'
           CNAME(5)  =  'Algal Group 2' 
           CNAME(6)  =  'Algal Group 3'
           CNAME(7)  =  'Zooplankton Group 1'
           CNAME(8)  =  'Zooplankton Group 2'
           CNAME(9)  =  'Labile DOC'           
           CNAME(10) =  'Refractory DOC'      
           CNAME(11) =  'Labile POC'                 
           CNAME(12) =  'Refractory POC'      
           CNAME(13) =  'Ammonium'
           CNAME(14) =  'Nitrate-nitrite' 
           CNAME(15) =  'Urea'             
           CNAME(16) =  'Labile DON'                 
           CNAME(17) =  'Refractory DON'      
           CNAME(18) =  'Labile PON'                 
           CNAME(19) =  'Refractory PON'      
           CNAME(20) =  'Total phosphate'
           CNAME(21) =  'Labile DOP'                 
           CNAME(22) =  'Refractory DOP'      
           CNAME(23) =  'Labile POP'                 
           CNAME(24) =  'Refractory POP'      
           CNAME(25) =  'Particulate Inorganic P'
           CNAME(26) =  'COD'
           CNAME(27) =  'Dissolved oxygen'
           CNAME(28) =  'Particulate silica  '
           CNAME(29) =  'Dissolved silica'
           CNAME(30) =  'Internal P Group 1'  !Phosphorus in algae 1
           CNAME(31) =  'Internal P Group 2'  !Phosphorus in algae 2
           CNAME(32) =  'Internal P Group 3'  !Phosphorus in algae 3

           RETURN
           END SUBROUTINE INIT_FILE_INFO

           END MODULE MOD_FILEINFO

