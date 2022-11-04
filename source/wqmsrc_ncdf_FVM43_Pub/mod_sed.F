
!************************************************************************
!**                                                                    **
!**                  Chesapeake Bay Sediment Model                     **
!**                                                                    **
!** Third version (speeded up) received from Fitzpatrick May 8, 1990   **
!**              First modified by C. Cerco May 9, 1990                **
!**           Final modifications D.M. Di Toro, Jan 27, 1992           **
!**                                                                    **
!** Modified (susp feeder inputs; deposit feeder dynamics)             **
!**      -- HydroQual, Inc. (S. Lowe, M. Meyers) Latest: 18 Apr 97     **
!**                                                                    **
!** Modified (updated methanogenesis-sulfate reduction algorithms)     **
!**      -- HydroQual, Inc. (N. Churi, D. Di Toro) Latest: 18 Apr 97   **
!**                                                                    **
!** org flux to sediment ==> diagenesis ==> inorg flux to water column **
!**                                                                    **
!**     3-G model - G1=labile, G2=refractory, G3=slow refractory       **
!**                                                                    **
!** Problems with CH4SAT, DDSO4 fixed by CFC.  DDSO4 11/3/06           **
!**                                                                    **
!************************************************************************
!**                                                                    **
!** Inputs                                                             **
!**                                                                    **
!**            Required inputs for sediment sub-model                  **
!**                                                                    **
!** A.  Passed to sediment subroutine from water quality subroutine    **
!**                                                                    **
!**   1.  Overlying water column segment volume (V1)                   **
!**   2.  Overlying water column segment depth (BL(ISEG,3))            **
!**   3.  Overlying water column depth                                 **
!**   4.  Overlying water column segment temperature and salinity      **
!**   5.  Overlying water column ammonia, nitrate, phosphate, silica   **
!**       and dissolved oxygen concentrations                          **
!**                                                                    **
!** B. Inputs supplied via direct input to the sediment subroutine     **
!**                                                                    **
!**  Variable names        Description                         Units   **
!**                                                                    **
!**                                                                    **
!**     HSEDALL      Depth of sediment layer  (h2)               cm    **
!**      DIFFT       Water column-sediment layer diffusion             **
!**                  coefficient for temperature                       **
!**                  cm^2/sec on input, model converts it to m^/sec    **
!**     SALTSW       Salinity concentration for determining      ppt   **
!**                  whether methane or sulfide SOD formulation        **
!**                  is to be used.  Also determines PO4 sorption.     **
!**     SALTND       Determines whether fresh or saltwater             **
!**                  nitrification/denitrification rates are used      **
!**                                                                    **
!**   Diagenesis stoichiometry                                         **
!**                                                                    **
!**                 Fractions of G1, G2, and G3 contained in ...       **
!**                                                                    **
!**   FRPALG1(3)     Algal group no.1 phosphorus                        **
!**   FRPALG2(3)     Algal group no.2 phosphorus                        **
!**   FRPALG3(3)     Algal group no.3 phosphorus                        **
!**   FRPOP(MTLOC,3) Non-algal particulate organic phosphorus           **
!**   FRNALG1(3)     Algal group no.1 nitrogen                          **
!**   FRNALG2(3)     Algal group no.2 nitrogen                          **
!**   FRNALG3(3)     Algal group no.3 nitrogen                          **
!**   FRPON(MTLOC,3) Non-algal particulate organic nitrogen             **
!**   FRCALG1(3)     Algal group no.1 carbon                            **
!**   FRCALG2(3)     Algal group no.2 carbon                            **
!**   FRCALG3(3)     Algal group no.3 carbon                            **
!**   FRPOC(MTLOC,3) Non-algal particulate organic carbon               **
!**                                                                    **
!**   Diagenesis kinetics                                              **
!**                                                                    **
!**   KPDIAG(3)     Reaction rates for POP G1, G2, and G3        /day  **
!**   DPTHTA(3)     Temperature thetas for POP G1, G2, and G3          **
!**   KNDIAG(3)     Reaction rates for PON G1, G2, and G3        /day  **
!**   DNTHTA(3)     Temperature thetas for PON G1, G2, and G3          **
!**   KCDIAG(3)     Reaction rates for POC G1, G2, and G3        /day  **
!**   DCTHTA(3)     Temperature thetas for POC G1, G2, and G3          **
!**   KSI           Reaction rate for Part. Biogenic Si (PSISED) /day  **
!**   THTASI        Temperature theta for PSISED                       **
!**   THTASISAT     Temperature theta for Si saturation concentration  **
!**                                                                    **
!**   Solids and transport                                             **
!**                                                                    **
!**   VPMIX(MTLOC)   Particulate diffusion coefficient (Dpp)   m**2/day  **
!**   THTADP        Temperature theta for Dpp                           **
!**   VDMIX(MTLOC)   Porewater diffusion coefficient (Ddp)     m**2/day  **
!**   THTADD        Temperature theta for Ddp                           **
!**      M1         Concentration of solids in layer 1       kg/l      **
!**      M2         Concentration of solids in layer 2       kg/l      **
!**                                                                    **
!**   Reaction kinetics                                                **
!**                                                                    **
!**    KAPPNH4F     Nitrification reaction velocity                    **
!**                 for freshwater in layer 1                m/day     **
!**    KAPPNH4S     Nitrification reaction velocity                    **
!**                 for saltwater  in layer 1                m/day     **
!**    PIENH4       Ammonia partition coefficient            L/kg      **
!**    THTANH4      Theta for nitrification reaction velicities        **
!**    KMNH4        Nitrification half saturation constant             **
!**                 for ammonia                              mg N/m3   **
!**    KMNH4O2      Nitrification half saturation constant             **
!**                 for oxygen                               mg O2/L   **
!**    KAPPNO3F     Denitrification reaction velocity                  **
!**                 for freshwater in layer 1                m/day     **
!**    KAPPNO3S     Denitrification reaction velocity                  **
!**                 for saltwater  in layer 1                m/day     **
!**    K2NO3        Denitrification reaction velocity                  **
!**                 in layer 2                               m/day     **
!**    THTANO3      Theta for denitrification                          **
!**    KAPP1HSD       Reaction velocity for dissolved sulfide            **
!**                 oxidation in layer 1                     m/day     **
!**    KAPP1HSP       Reaction velocity for particulate sulfide          **
!**                 oxidation in layer 1                     m/day     **
!**    PIE1HS        Partition coefficient for sulfide                  **
!**                 in layer 1                               L/kg      **
!**    PIE2HS        Partition coefficient for sulfide                  **
!**                 in layer 2                               L/kg      **
!**    THTAH2S      Theta for both dissolved and particulate           **
!**                 sulfide oxidation                                  **
!**    KMHSO2       Sulfide oxidation normalization constant           **
!**                 for oxygen                               mg O2/L   **
!**    CSISAT       Saturation concentration for pore water            **
!**                 silica                                     mg Si/m3  **
!**    DPIE1SI      Incremental partition coefficient for              **
!**                 silica in layer 1                        L/kg      **
!**    PIE2SI       Partition coefficient for silica in                **
!**                 layer 2                                  L/kg      **
!**    KMPSI        Particulate biogenic silica half saturation        **
!**                 constant for dissolution                 mg Si/m3  **
!**    O2CRITSI     Critical dissolved oxygen concentration            **
!**                 for layer 1 incremental silica sorption  mg O2/L   **
!**    JSIDETR      Detrital biogenic silica source to                 **
!**                 sediment                                 mg Si/m2-d**
!**    DPIE1PO4F    Incremental partition coefficient for              **
!**                 phosphate in layer 1 (freshwater)        unitless      **
!**    DPIE1PO4S    Incremental partition coefficient for              **
!**                 phosphate in layer 1 (saltwater)         L/kg      **
!**    PIE2PO4      Partition coefficient for phosphate                **
!**                 in layer 2                               L/kg      **
!**    O2CRITPO4       Critical dissolved oxygen concentration for        **
!**                 layer 1 incremental phosphate sorption   mg O2/L   **
!**    KMO2DP       Particle mixing half saturation constant           **
!**                 for oxygen                               mg O2/L   **
!**    TEMPBEN      Temperature at which benthic stress                **
!**                 accumulation is reset to zero            deg C     **
!**    KBENSTR      Decay constant for benthic stress        /day      **
!**    KLBNTH       Ratio of bio-irrigation to bioturbation            **
!**    DPMIN        Minimum particle diffusion coefficient   m2/day    **
!**    KAPPCH4      methane oxidation reaction velocity      m/day     **
!**    THTACH4      theta for methane oxidation                        **
!**                                                                    **
!** Output                                                             **
!**                                                                    **
!**    The subroutine returns fluxes for                               **
!**                                                                    **
!**     JSOD, JAQSOD, JCH4AQ and JCH4G  [g o2*/m2-day]                **
!**     JNH4, JPO4, JNO3 and JSI  [mg/m2-day]                          **
!**                                                                    **
!**    via array BFLUX in COMMON /BENTHC/                              **
!**                                                                    **
!************************************************************************

MODULE MOD_SED
    
        USE MOD_PREC,     ONLY: SP
        
        USE MOD_SIZES,     ONLY: MGL, NFLP
        
        USE MOD_WQM, ONLY:        &!
                BFOFN,          &!
                WCLPTR,            &!
                NTVWCLF,        &!
                DLT,            &!
				TMSTRT,			&
                SETTLING,       &!
                B,              &!
                B1,             &!
                B2,             &!
                B3,             &!

                LPOP,           &!
                Q1,             &!
                Q2,             &!
                Q3,             &!
                SAV_CALC,       &!
                WS1NET,         &!
                WS2NET,         &!
                WS3NET,         &!
                WSLNET,         &!
                WSRNET,         &!
                WSSNET,         &!
                WSUNET,         &!
                ANC1,           &!
                ANC2,           &!
                ANC3,           &!
                ASC1,           &!
                ASC2,           &!
                ASC3,           &!
                
!                MTVEL,          &!                
                
!                JPOC,         &!
!                JPON,         &!
!                JPOP,         &!                
!                JPOS,         &!

!                CCH4,           &!
!                CHS,            &!
!                CSI,            &!
!                CSO4,           &!                
!                CNH4,           &!
!                CNO3,           &!
!                CPO4,           &!
!                CPOC,           &!
!                CPON,           &!
!                CPOP,           &!
!                CPOS,           &!
                CTEMP,          &!
                
                LPOC,           &!
                LPON,           &!
                PCFWS,          &!
                PNFWS,          &!
                PPFWS,          &!
                PSFWS,          &!
                SIUPB,           &!
                RPOC,           &!
                RPON,           &!
                RPOP,           &!
                SSFWS,          &!
                SSI,            &!
!                BB,             &!
                BENCH4A,        &!
                BENCH4G,        &!
                BENCOD,         &!
                BENDEN,         &!gN/m^2/day
                BENDO,          &!
               !  BENDOC,         &!
                BENNH4,         &!
                BENNO3,         &!
                BENPO4,         &!
                BENSA,          &!
                COD,            &!
                DOXG,           &!
                SIAT,            &!
                JDAY,           &!
                KADPO4,         &!
                KADSA,          &!

                NH4,            &!
                NO3,            &!
                PIP,            &!
                PO4,            &!
                SALT,           &!
!                DIAGN,          &!
                T                 !
        
        USE MOD_FILEINFO, ONLY :         & !
                !,DIA            &!
                !,CBC             &!
                !,S1            &!
                !,S2            &!
                !,S3            & !                 
                BFI                &!
                ,WCL            &!
                !,MET            &!
                ,BFO    !        & !
                !,KEI            &!
                !,ATM            &!
                !,STL            & !
                !,AGR            & !
                !,SVI            & !
                !,SVO            & !
                !,KFL            & !
                !,ZOO            & !
                !,ZFO            & !
                !,ALO              &!
                !,CON            &!
                !,RSO            &!
                !,SNP            &!
                !,PLT            &!
                !,APL             &!
                !,TFL             &!
                !,OPL            &!
                !,SFI            &!
                !,SFO            &!
                !,MAP             &!
                !,ICI             &!
                !,ICO            &!
                !,MRL            &!
                !,MBL            &!
                !,RSI            &!
                !,UNIT_LINKAGE    &!
                !,UNIT_STN        &!
                !,UNIT_HIS        & !           
                !,CNAME            &!
                !,CONFN
                
USE MOD_LIMS, ONLY : MLOC, MTLOC, KBM1

#  if defined (MULTIPROCESSOR)
        USE MOD_PAR, ONLY:   &      
            !EL_PID,    &        !(:)           !!PROCESSOR OWNER OF GLOBAL ELEMENT
            !ELID,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT
            !NLID,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE 
            !ELID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
            !NLID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
            !OBN_GL2LOC,    &    !(:)           !!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
            !EGID,        &    !(:)           !!GLOBAL ID OF LOCAL ELEMENT
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
            NHN,          &    !            !!NUMBER OF HALO NODES
            HN_LST!,     &    !            !!LIST OF HALO NODES 
            !HN_OWN,     &    !            !!PRIMARY OWNER OF HALO NODES
            !EC,            &    !
            !NC,            &    !
            !BNC,         &    !
            !EMAP,        &    !
            !NMAP,         &    !
            !PNE,        &    !(:)        !!NUMBER OF ELEMENTS IN EACH PROC
            !PNN,        &    !(:)        !!NUMBER OF NODES IN EACH PROC
            !PNHE,        &    !(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
            !PNBN,        &    !(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
            !PMBM,        &    !(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
            !PNHN,        &    !(:)       !!NUMBER OF HALO NODES IN EACH PROC
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
            !SORT,    &
            !PREAD,    &
            !PWRITE,&
            !PPRINT, &
            !LOCK, &
            !UNLOCK, &
            !GETLOC
#  endif  

USE MOD_HYDROVARS, ONLY:   &
           !GRAV    &        !
        !,PI    &        !
        !,PI2    &        !
        !,ZERO    &        !
        !,ONE_THIRD    &    !
        !,NVG    &        !
        !,XG    &        !GLOBAL X-COORD AT NODE 
        !,YG    &        !GLOBAL X-COORD AT NODE 
        !,HG    &        !GLOBAL DEPTH AT NODE 
        !,XCG    &        !GLOBAL X-COORD AT FACE CENTER 
        !,YCG    &        !GLOBAL X-COORD AT FACE CENTER 
        !,VXMIN    &        !
        !,VYMIN    &        !
        !,VXMAX    &        !
        !,VYMAX    &        !
        !,XC    &        !X-COORD AT FACE CENTER 
        !,YC    &        !Y-COORD AT FACE CENTER
        !,VX    &        !X-COORD AT GRID POINT
        !,VY    &        !Y-COORD AT GRID POINT
        !,ART    &        !AREA OF ELEMENT
        ART1    &        !AREA OF NODE-BASE CONTROl VOLUME
        !,ART2    &        !AREA OF ELEMENTS AROUND NODE
        !,NV    &        !NODE NUMBERING FOR ELEMENTS
        !,NBE    &        !INDICES OF ELMNT NEIGHBORS
        !,NTVE    &        !
        !,NTSN    &        !
        !,ISONB    &        !NODE MARKER = 0,1,2 
        !,ISBC    &        !
        !,ISBCE    &        !
        !,IEC    &        !
        !,IENODE &        !
        !,NBSN    &        !
        !,NIEC    &        !
        !,NTRG    &        !
        !,NBVE    &        !
        !,NBVT    &        !
        !,LISBCE_1    &    !LIST OF ELEMENTS WITH ISBCE=1
        !,LISBCE_2    &    !LIST OF ELEMENTS WITH ISBCE=2
        !,LISBCE_3    &    !LIST OF ELEMENTS WITH ISBCE=3
        !,DLTXC    &        !
        !,DLTYC    &        !
        !,DLTXYC    &    !
        !,DLTXE    &        !
        !,DLTYE    &        !
        !,DLTXYE    &    !
        !,SITAC    &        !
        !,SITAE    &        !
        !,XIJC    &        !
        !,YIJC    &        !
        !,XIJE    &        !
        !,YIJE    &        !
        !,EPOR    &        !ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
        !,IBCGEO    &    !LOCAL GEOSTROPHIC FRICTION CORRECTION NODES

        !,Z    &            !SIGMA COORDINATE VALUE 
        !,ZZ    &        !INTRA LEVEL SIGMA VALUE
        ,DZ    &        !DELTA-SIGMA VALUE
        !,DZZ    &        !DELTA OF INTRA LEVEL SIGMA 
        !,H1    &        !BATHYMETRIC DEPTH 
        ,H    &            !BATHYMETRIC DEPTH 
        ,D    !&            !CURRENT DEPTH 
        !,DT    &        !DEPTH AT PREVIOUS TIME STEP
        !,DT1    &        !DEPTH AT PREVIOUS TIME STEP
        !,EL    &        !CURRENT SURFACE ELEVATION
        !,ET    &        !SURFACE ELEVATION AT PREVIOUS TIME STEP
        !,DTFA    &        !ADJUSTED DEPTH FOR MASS CONSERVATION
        !,UU    &        !X-VELOCITY
        !,VV    &        !Y-VELOCITY
        !,UUT    &        !X-VELOCITY FROM PREVIOUS TIMESTEP
        !,VVT    &        !Y-VELOCITY FROM PREVIOUS TIMESTEP
        !,WWT    &        !Z-VELOCITY FROM PREVIOUS TIMESTEP
        !,WTST    &        !Vertical velocity in sigma from PREVIOUS TIMESTEP
        !,UARD_OBCNT    &!tykim
        !,XFLUX_OBCT    &!tykim
        !,DTFAT    &        !tykim
        !,TT_T    &        !tykim
        !,SALTT    &        !tykim
        !,WTS    &        !VERTICAL VELOCITY IN SIGMA SYSTEM
        !,UARD_OBCN    &    ! tykim 
        !,XFLUX_OBC    &    ! tykim 
        !,WTTS    &        !VERTICAL VELOCITY IN SIGMA SYSTEM 
        !,KH    &        !TURBULENT DIFFUSIVITY
        !,A1U    &        !
        !,A2U    &        !
        !,AWX    &        !
        !,AWY    &        !
        !,AW0    &        !
        !,VISCOFH    &    !
        !,UNC1    &        !
        !,VNC1    &        !
        !,WNC1    &        !
        !,WTSNC1    &        !
        !,UARD_OBCNNC1    &    !
        !,XFLUX_OBCNC1    &    !
        !,DTFANC1    &        !
        !,KHNC1    &        !
        !,TNC1    &        !
        !,SNC1    &        !
        !,ELNC1    &        !
        !,UNC2    &        !
        !,VNC2    &        !
        !,WNC2    &        !
        !,WTSNC2    &    !
        !,UARD_OBCNNC2    &!
        !,XFLUX_OBCNC2    &!
        !,DTFANC2    &    !
        !,KHNC2    &        !
        !,TNC2    &        !
        !,SNC2    &        !
        !,ELNC2    &        !
        !,num_hyd_ints    &!number of records in each hydrodynamics netcdf file
        !,TIME_MAP    &    !
        !,THOUR1    &    !SIMULATION TIME AT END OF CURRENT EXTERNAL STEP (IEXT) IN HOURS
        !,THOUR    &        !
        !,NCFILE_DIR    &!
        !,NCFILE_PREFIX    &!
        !,NCFILE_SUFFIX    &!
        !,NCFILE_NUMBER    &!
        !,FORMAT_STR    &!
        !,hydro_dir,     &    ! directory name where hydrodynamics results (netcdf) files are stored
        !,hydro_prefix, &    ! prefix of file name, e.g. 'psm_'
        !,hydro_suffix    &    ! suffix of filename, e.g. '.nc'
        !,hydro_filenumwidth, &    ! number of digits in filename following hydro_prefix, e.g. 4 for psm_0002.nc
        !,hydro_filenumstart, &    ! starting number of the file name in the digital part of the file name, e.g. 185 for psm_0185.nc
        !,hydro_Nrec    &        ! number of records in each of hydrodynamics file
        !,hydro_dlt    &            ! time step in hydrodynamics file (in seconds), e.g. 100 for 100sec
        !,t_his_start    &        !
        !,t_his_end    &            !
        !,t_his_dlt    &            !starting time, ending time, and interval of history outputs (days)
        !,Nstation    &            !
        !,NstationNum_GL    &    !maximum number of station is NstationMax!
        !,t_stn_start    &        !
        !,t_stn_end    &            !
        !,t_stn_dlt    &            !starting time, ending time, and interval of station outputs (days)
        !,STNFN    &                !file name for station output
        !,HISFN    &                !file name for history output
        !,HISFN_PREFIX    &        !prefix of history output file
        !,HISFN_EXT    &            !extention name of history output file
        !,HISFN_FINAL    &        ! 
        !,HISFN_SPLIT_BYLEVEL    &!True or False for splitting history output in files level by level (default is .FALSE.)
        !,hydro_netcdf    &        !
        !,wqm_history    &        !
        !,wqm_stations    &        !
        !,IFNC    &                !file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise 
        !,NTRECNC    &            !time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file. 
        !,NTHYDRO                !overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read

    !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY :         &
            SERIAL          &           !!TRUE IF SINGLE PROCESSOR
            ,MSR            &           !!TRUE IF MASTER PROCESSOR (MYID==1)
            ,PAR        !    &           !!TRUE IF MULTIPROCESSOR RUN
            !,CASENAME      &           !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
            !,CASETITLE      &              !!CASE TITLE                                 
            !,HMAX           &              !!GLOBAL MAXIMUM DEPTH
            !,HMIN           &              !!GLOBAL MINIMUM DEPTH
            !,UMOL           &              !!VERTICAL DIFFUSION COEFFICIENT
            !,HORCON         &              !!HORIZONTAL DIFFUSION COEFFICIENT
            !,DTI            &              !!internal time step
            !,HORZMIX        &           !!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
            !,FILENUMBER    &            !!
            !,PREFF            &            !!
            !,INPDIR        &            !!
            !,GEOAREA        &            !!
            !,RIV_FILENUMBER    &            !!
            !,INFLOW_TYPE       &            !!SPECIFIED RIVER INFLOW TYPE (edge/node) 
            !,POINT_ST_TYPE     &            !!(calculated/specified)
            !,PNT_SOURCE        &            !!point_source
            !,DAY                &
            !,in_jday        
        
        USE MOD_SF, ONLY:         & 
                SFEEDER!,         &    !
                !HYPOXFX_SF,        &    !
                !FRSASF,        &    !
                !NSPECIES,        &    !
                !SEDTYPE,        &    !    (:,:)
                !FILT,            &    !                
                !SFA1,            &    !
                !SFA2,            &    !
                !SFA3,            &    !
                !SFA4,            &    !
                !SFA5,          &    !
                !MAXING,        &    !
                !SFCN,            &    !
                !SFCP,          &    !
                !SFRESP,        &    !
                !SFPRED,        &    !
                !SFTMN,         &    !
                !THTAFILT,      &    !
                !THTARESP,         &    !
                !THTAPRED,         &    !
                !XKPO2,         &    !
                !SFTD,          &    !
                !FILTFACT,         &    !
                !RESPFACT,      &    !
                !SFDOH,         &    !
                !SFDOQ,         &    !
                !SFATURB,         &    !
                !SFBTURB,        &    !
                !SF_JLPOC,        &!(:),            
                !SF_JLPON,        &!(:),            
                !SF_JLPOP,        &!(:),            
                !SF_JRPOC,        &!(:),              
                !SF_JRPON,        &!(:),            
                !SF_JRPOP,        &!(:),            
                !JNH4SF,            &!(:),            
                !JPO4SF,            &!(:),            
                !SODSF,            &!(:),            
                !JSASF,            &!(:),            
                !JSUSF,            &!(:),            
                !SF_SSI,            &!(:),            ! mg Solids/m^2/day flux (positive to sediments)
                !SF_SU,            &!(:),        
                !SF_SA,            &!(:),            
                !SF_PIP,            &!(:),            
                !SFEED !,        &!(:,:),        
                !SFEED_GL,        &!(:,:)
                !SFFILTCT,        &!(:,:),        
                !SFRESPCT,        &!(:,:),          
                !SFPREDCT,        &!(:,:),        
                !SFRATEI,        &!(:,:),           
                !SFGCMAX,        &!(:,:),        
                !SFGMIN,        &!(:,:),        
                !SFGCING,        &!(:,:),        
                !SFCFECES,        &!(:,:),          
                !SFNFECES,        &!(:,:),        
                !SFPFECES,        &!(:,:),          
                !SFCPSFEC,        &!(:,:),        
                !SFNPSFEC,        &!(:,:),          
                !SFPPSFEC,        &!(:,:),        
                !SF_SFGC,        &!(:,:),         
                !SF_RESP,        &!(:,:),        
                !SF_PRED,        &!(:,:),        
                !SF_RMORT,        &!(:,:),        
                !ASF_SFGC,        &!(:,:),        
                !ASF_RESP,        &!(:,:),        
                !ASF_PRED,        &!(:,:),        
                !ASF_RMORT,        &!(:,:)        
                !SUSPFEED,        &!        
                !SF_ALLOC,        &!    
                !SF_DEALLOC    &!
                !SF_READ,        &!    
                !SF_INIT
                
    USE MOD_SED_SF_EXCHANGE_VARS, ONLY: JLPOC_SED_SF,&
                                        JLPOC_SED_SF,&
                                        JLPON_SED_SF,&
                                        JLPOP_SED_SF,&
                                        JRPOC_SED_SF,&
                                        JRPON_SED_SF,&
                                        JRPOP_SED_SF,&
                                        JNH4_SED_SF,&
                                        JPO4_SED_SF,&
                                        SOD_SED_SF,&
                                        JSA_SED_SF,&
                                        JSU_SED_SF,&
                                        SSI_SED_SF,&
                                        SU_SED_SF,&
                                        SA_SED_SF,&
                                        PIP_SED_SF
                
 
    USE MOD_BA, ONLY:     BALGAE_CALC,     &
                        BA_CALC,        &  
                        BANH4,            &
                        BANO3,            & 
                        BAPO4,             &
                        BADOC,             &
                        BADO,            &
                        DTPOC1_BA,         &
                        DTPOC2_BA,        &
                        DTPOC3_BA,        &
                        DTPON1_BA,        &
                        DTPON2_BA,        &
                        DTPON3_BA,        &
                        DTPOP1_BA,        &
                        DTPOP2_BA,        &
                        DTPOP3_BA,        &
                        DTAG3CFL_BA,    &
                        DTAG3NFL_BA,    &
                        DTAG3PFL_BA
                    
    USE MOD_DF, ONLY: DFEEDER,              &    !
                      DF_GROW_POC1,        &    !
                      DF_GROW_POC2,        &    !
                      DF_MORT,            &    !
                      DF_RESP,            &    !                      
                      DF_PRED,            &    !
                      DF_SOD,            &    !
                      DFEEDM1S,            &    !
                      AMCN,                &    !
                      AMCP,                &    !
                     !DFEEDM1S_GL,        &    !
                     !DFEEDM1S,            &    !
                     !DF_READ,            &    !
                     !DF_INIT,            &    !
                     !DF_DEALLOC,        &    !
                     !DF_ALLOC,            &    !
                      DF_CALC
					  
#if defined(SEDIMENT_DOM)		  
    USE MOD_SED_DOM, ONLY: SED_DOM, MARSH_SEDIMENTS, &! SEDIMENT DOM AND MARSH SEDIMENTS FLAG	
							SED_DOM_CALC,   HYDRO_FLAG, &
							VMAX_MARSH_NH4, VMAX_MARSH_NH4_OUT, WETLAND_FLAG, PLANT_NH4_UPTAKE, kmnh4_marsh,MARSH_FRAC,MARSH_MASS

   USE MOD_SED_DOM_EXCHANGE, ONLY:     &
	    JDOCX_SHARE1, &
		JDOCX_SHARE2, &
	    JDONX_SHARE1, &
		
	    JDONX_SHARE2, &
	    JDOPX_SHARE1, &
		    JDOPX_SHARE2, &
		SED_TEMP_SHARE,   &
		ZSED_SHARE,   &
		W2_SHARE,   &
		DLTS_SHARE, &
		FX_POC1TM1,  &   !  AdFded by B Clark Nov 2015 for hydrolysis
		FX_POC2TM1,   &
		FX_POC3TM1,  &
		FX_PON1TM1,  &
		FX_PON2TM1,   &
		FX_PON3TM1,    &	
		FX_POP1TM1,  &
		FX_POP2TM1,   &
		FX_POP3TM1, &
		SED_DOM_FLAG, & S_SHARE, KL12_SHARE, HSED1_SHARE

#endif   
    IMPLICIT NONE
    SAVE
	
	!LOGICAL :: SED_DOM = .TRUE.

!
!variables that are used for as flags of model simulation 
!

      !WLong moved this from wqm_modules.F

      LOGICAL :: BENTHIC_OUTPUT,              & !flag for benthic model output
                 STEADY_STATE_SED_G3,        & !flag for steady state simulation for G3  
                 STEADY_STATE_SED_IC,        & !flag for setting initial condition of sediments using steady state solution
                                               !     with first time step forcing from water column to sediments and steady state formulation
                 STEADY_STATE_SED_TS,        & ! quasi-steady state time series calculation of sediment flux 
                 STEADY_STATE_SED_AVG,       &   !flag for calculating steady state using time averaged forcing
                                              !    from water colum to sediments and using steady state formulation
                 AGGREGATE_POM_FLUX           ! If true, aggregate all POM reaching t/JDOChe sediments before distributing into 3 sediment labilities
				 
 !     LOGICAL :: SED_DOM =.FALSE.                     ! B Clark , TRUE if sediment DOM is included in model			
                                               

!WLong                 BALGAE_CALC            !flag for activating benthic algae
!WLong                                        !calculation
				 
!	LOGICAL :: SED_DOM = .TRUE.  ! B Clark added in for easy compilation	 

    
!
! variables that are local to sediment module
!
    REAL(SP) :: DLTS                                   !sediment simulation time step (day)
                                                    !==DLT/86400, where DLT is in seconds
 
    REAL(SP),ALLOCATABLE,DIMENSION(:) :: HSED             !sediment layer thickness (H2 (m))    
    REAL(SP),ALLOCATABLE,DIMENSION(:) :: HSED1             !sediment layer thickness (H1 (m))    
!
! input variables that control the sediment module's simulation
!
    INTEGER      :: NDTSED    !LB number of WC timesteps between sedim calculations (DTsed = NDTSED * DTwater)
    INTEGER      :: SSTATEG3       !==1 for steady state calculation of G3 POC PON POP in sediments
                                !    based on time average of the simulation at the end of the integration
                                !
    
    INTEGER     :: SSTATEIC     !==1 for setting initial condition of sediment using steady state solution
                                !   based on overlying water column concentration and fluxes of the first time step

    INTEGER     :: QUASISS      !==1 for quasi-steady state calculation, meaning that even the forcing flux (JPOC, JPON, JPOP, JPOS)
                                !   and overlying water column is time varying
                                !   the calculation of sediment module is using
                                !   steady state. 

                                !
    
    INTEGER     :: SSTATEAVG    !==1 for calculating sediment model using time averaged overlying water column
                                !   forcing and time averaged POM fluxes to sediments, and steady state formulation
                                !   The results of this can be used for subsequent simulations as initial conditions
                                !   when SSTATEIC /= 1
    
    INTEGER      :: SWITCH_LAB  !==1 to aggregate all POM reaching the sediments before distributing into 3 sediment labilities
                                !  if 0, then LPOMflux goes to G1 in the sediments and RPOM is divided in 65% to G2 and 35% to G3
    
    
    INTEGER        :: ITVWCLF        !==1 for forcing with time varying overlying water column condition input file
    REAL(SP)    :: WCLFDAYS        !length of water column forcing for each WCLFN file (days)
                        
    REAL(SP)     :: HSEDALL     !Depth of sediment layer (h1+h2 ~= h2) (cm)
    REAL(SP)      :: DIFFT        !water column-sediment diffusion rate (m^2/sec)
    
    REAL(SP)      :: SALTSW      !Salinity concentration for determining whether methane or sulfide SOD formulation (ppt)
                                    !if salt water  (> SALTSW) use sulfide, 
                                    !if freshwater (<=SALTSW) use methane
    
    REAL(SP)     :: SALTND        !Salinity for determing whether to use freshwater (<=SALTND) or salt water (>SALTND)
                                !nitrification/denitrification (ppt)
    
    REAL(SP),DIMENSION(3) ::     FRPALG1,    & !Algal group 1 Phosphorus fraction for
                                                !G1,G2,G3 in phytoplankton
                                FRPALG2,    & !Algal group 2 Phosphorus fraction for.
                                                !G1,G2,G3 in phytoplankton
                                FRPALG3,    & !Algal group 3 Phosphorus fraction for.
                                                !G1,G2,G3 in phytoplankton
                                FRNALG1,    & !Algal group 1 Nitrogen fractions for
                                                !G1,G2,G3 in phytoplankton
                                FRNALG2,    & !Algal group 2 Nitrogen fractions for.
                                                !G1,G2,G3 in phytoplankton
                                FRNALG3,    & !Algal group 3 Nitrogen fractions for G1
                                                !G2 G3 in phytoplankton
                                FRCALG1,    & !Algal group 1 Carbon fractions for G1,
                                                !G2, G3 in phytoplankton
                                FRCALG2,    & !Algal group 2 Carbon fractions for G1,G2
                                                !G3 in phytoplankton
                                FRCALG3!,    & !Algal group 3 Carbon fractions for G1,
!                                                !G2, G3 in phytoplankton
!                                 FRCPHB,    & !Fraction of POC generated for G1,G2,G3
!                                                !due to predation on benthic algae 
!                                 FRNPHB,    & !Fraction of PON genereated for G1,G2,G3
!                                               !due to predation on benthic algae
!                                 FRPPHB       !Fraction of POP generated for G1,G2,G3
!                                                !due to predation on benthic algae
                         
    REAL(SP),ALLOCATABLE::  FRPOP(:,:), & !MTLOC x 3, Fraction of non-algal particulate organic phosphorus (dimensionless)
                            FRPON(:,:), & !MTLOC x 3, Fraction of non-algal particulate organic nitrogen (dimensionless)
                            FRPOC(:,:)    !MTLOC x 3, Fraction of non-algal particulate organic carbon (dimensionless)

    REAL(SP)     ::         FRPON1,     &! Fraction of G1 in PON (detritus) flux from water column
                            FRPOP1,        &! Fraction of G1 in POP (detritus) flux from water column
                            FRPOC1         ! Fraction of G1 in POC (detritus) flux from water column

                            
    REAL(SP)     ::         FRACL1,     &! Fraction of total POM flux from water column sent to G1 in sediments
                            FRACL2,     &! Fraction of total POM flux from water column sent to G2 in sediments
                            FRACL3       ! Fraction of total POM flux from water column sent to G2 in sediments
                            
  !WLong moved here from wqm_modules.F
  !sediment reaction kinetics.
  
    REAL(SP)    ::             KPDIAG(3), &  !G1,G2,G3 diagenesis rate for P (1/day)
                            DPTHTA(3), &  !Temperature control theta for P
                            KNDIAG(3), &  !G1,G2,G3 diagenesis rate for N (1/day)
                            DNTHTA(3), &  !Temperature control theta for N
                            KCDIAG(3), &  !G1,G2,G3 diagenesis rate for C (1/day)
                            DCTHTA(3)     !Temperature control theta for C

    REAL(SP) :: KSI,       &  !reaction rate for particulate biogenic Si (PSISED) (1/day)
            THTASI,        &  !temperature control theta for KSI
            THTASISAT         !temperature control theta for Si saturation concentration 
                    
    
     REAL(SP) :: M1,       & !layer 1 solids concentration  (kg/L)
             M2,       & !layer 2 solids concentration  (kg/L)
             THTADP,   & !temperature theta for particle diffusion coefficient
                         !(Dpp)
             THTADD      !temperature tehat for porewater diffusion coefficient
                         !(Ddp)
    
    !Ammonia nitrification kinetics.
    REAL(SP) :: KAPPNH4F, & !Nitrification reaction velocity (m/day) (integrated over layer 1 thickness)
                        !for freshwater in layer 1
            KAPPNH4S, & !Nitrification reaction velocity (m/day) for salt water
                        !in layer 1 (integrated over layer 1 thickness)
            PIENH4,   & !Ammonia partition coefficient (L/kg)
            THTANH4,  & !temperature control theta for Theta for nitrification
                        !reaction velicities.
            KMNH4,    & !Nitrification half saturation constant for ammonia in layer 1
                        ! mg N/m^3 (for total NH4 in layer 1) (DMD book (3.30)
            THTAKMNH4,& !temperature contro theta for nitrification half saturation constant KMNH4
            KMNH4O2     !Nitrification half saturation constant for oxygen
                        ! mgO2/L
    
    
    !NO3 denitrification kinetics
    REAL(SP) :: KAPPNO3F,  &  !Denitrification reaction velocity (m/day) for fresh
                          !water in layer 1 (integrated over layer 1 thickness)
            KAPPNO3S,  &  !Denitrification reaction velocity (m/day) for salt
                          !water in layer 1 (integrated over layer 1 thickness)
            K2NO3,     &  !Denitrification reaction velocity in layer 2 (m/day) (DMD p 105 table 4.1)
                          !(integrated over the thickness of layer 2)
            THTANO3       !temperature control theta for denitrification
    
    !HS sulfide oxidation kinetics
    REAL(SP) :: KAPP1HSD,     & !Reaction velocity for dissolved sulfide oxidation in layer 1
                          !(m/day)
            KAPP1HSP,     & !Reaction velocity for particulate sulfide oxidation
                          !in layer 1 (m/day)
            PIE1HS,      & !Partition coefficient for sulfide in layer 1 (L/kg)
            PIE2HS,      & !Partition coefficient for sulfide in layer 2 (L/kg)
            THTAH2S,    & !temperature control theta for both dissolved and
                          !particulate sulfide
            KMHSO2        !Sulfide oxidation normalization constant
                          !for oxygen (mg O2/L)
    !Silica kinetics
    REAL(SP) :: CSISAT,  & !Saturation concentration for pore water silica (mgSi/m^3) 
                CSISATT20,  & !Saturation concentration for pore water silica (mgSi/m^3) at T=20degC
                DPIE1SI, & !Incremental partition coefficient for silica in layer 1 (unitless)
                PIE2SI,  & !Partition coefficient for silica in layer 2 (L/kg)
                KMPSI      !Particulate biogenic silica half saturation for silica.
                        !for dissolution (mgSi/m^3)
    
    REAL(SP) ::  O2CRITSI, & !Critical dissolved oxygen concentration for layer 1
                         ! incremental silica sorption  mg O2/L
                JSIDETR     !Detrital biogenic silica source to
                         ! sediment mg Si/m2-d

    !PO4 sorption kinetics
    REAL(SP) :: DPIE1PO4F,  &    !Incremental partition coefficient for
                             !phosphate in layer 1 (freshwater) (unitless)
             DPIE1PO4S, &    !Incremental partition coefficient for
                             !phosphate in layer 1 (saltwater)  (unitless)
             PIE2PO4,   &    !Partition coefficient for phosphate
                             !in layer 2                         (L/kg)
             O2CRITPO4,    &    !Critical dissolved oxygen concentration for
                             !layer 1 incremental phosphate sorption  (mg O2/L)
             KMO2DP          !Particle mixing half saturation constant
                             !for oxygen (mgO2/L) (bioturbation)

    !benthic stress and bioturbation etc mixing control
    REAL(SP) ::  TEMPBEN, &      !temperature at which benthic stress is reset to
                             !zero (degC)
             KBENSTR, &      !Decay constant for benthic stress (1/day)
             KLBNTH,  &      !Ratio of bio-irrigation to bioturbation
             DPMIN           !Minimum particle diffusion coefficient (m2/day)

    !methane kinetics
    REAL(SP) :: KAPPCH4, &      !methane oxidation reaction velocity      m/day
            THTACH4, &      !temperature control theta for methane oxidation
            KMCH4O2, &      !methane oxidation saturation concentration for oxygen (mgO2/L)
            KMSO4           !saturation concentration for SO4 in layer 2 for methane oxidation by sulfate (mgO2/L)
                            !

    !WLong moved here from wqm_modules.F
    !Deposit feeder.
            
!WLong     REAL(SP) :: XKMI0,    &
!WLong            ING0,     &
!WLong            THTAI0,   &
!WLong            R,        &
!WLong            THTAR,    &
!WLong            BETA,     &
!WLong            THBETA,   &
!WLong            AMCN,     &
!WLong            AMCP,     &
!WLong            AA1,      &
!WLong            AA2,      &
!WLong            XKMG1,    &
!WLong            XKMG2,    &
!WLong            XKBO2,    &
!WLong            TDD
            
!WLong     REAL(SP) :: DFDOH,    & !oxygen concentration for hypoxia control on bioturbation
!WLong          DFDOQ,      &    !oxygen concentration for hypoxia control on bioturbation
!WLong             DOLOW        !
    
    !WLong moved from subroutine SED_READ
!    CHARACTER(LEN=3)  :: BALC
                           
    !benthic algae growth kinetics    
    !WLong moved from wqm_modules.F
    !REAL(SP) :: &! PMB,         !
             !ANCB,         &!
             !APCB,         &!
             !KTGB1,     &!
             !KTGB2,     &!
             !TMB,       &!
!WLong             KESED,     &!
!Wlong             KEBALG,     &!
!WLong             KHNB,         &!
!WLong             KHPB,         &!
!WLong             KHRB,         &!
!WLong             BMRB,         &!
!WLong             BPRB,         &!
!WLong             KTBB,        &!
!WLong             TRB,          &!1
!WLong             FNIB,         &!
!WLong             FPIB           !

!    REAL(SP) :: ALPHB!,    &  !alpha of P-I curve in light calculation for benthic algae growth
!            CCHLB        !C/CHLA ratio of benthic algae gC/gCHLA
    
    !WLong moved from subroutine SED_READ to here
    CHARACTER(LEN=8)  :: SPVARS,   &
                         PRINTS

    !WLong has checked the units of these settling velocities
    !net settling rate is supposed to be smaller than normal settling rate
    !net settling rate = (normal settling - suspension rate) 
    !so the suspension rate (m/d) = (normal settling - net settling ) (m/d)
    REAL(SP) ::  WSSBNET,  & !net settling rate of inorganic suspended sediments (m/d) 
                                            WSLBNET,  & !net settling rate of LPOM (m/d)
                                            WSRBNET,  & !net settling rate of RPOM (m/d)
                                            WS1BNET,  & !net settling rate of alg 1 (m/d)
                                            WS2BNET,  & !net settling rate of alg 2 (m/d)
                                            WS3BNET,  & !net settling rate of alg 3 (m/d)
                                            WSUBNET     !net settling rate of particulate biogenic silica (m/d)
                              
    !WLong moved from subroutine SED_READ
    CHARACTER(LEN=8)  ::    SPVARB, &
                            PRINTB
                            
    REAL(SP), ALLOCATABLE,DIMENSION (:) ::     VSED,  &    !cm/yr on input, model converts to m/day
                                            VPMIX, &    !m^2/day
                                            VDMIX        !m^2/day
                        
                        
    !WLong moved from subroutine SED_READ
    CHARACTER(LEN=8)  :: SPVARLR,  &
                         PRINTLR
    
    
    
!
!other local variables that are used in the sediment module
!

       !WLong moved the following from wqm_modules.F
       !WLong renamed DP to DPP and DD to DDP here
       REAL(SP) :: DPP,       &   !Particulate diffusion coefficient (m2/day)
                                  !                                   <==VPMIX
               W2,            &   !Burial velocity (m/day)            <==VSED
               DDP,           &   ! porewater diffusion rate (m2/day) <==VDMIX
               H2                 !sediment layer thickness (m).


       REAL(SP) :: DPIE1PO4  !sediment 1st layer incremental partition coefficient for phosphate  (dimensionless)

       REAL(SP)  :: KAPPNH4,   &   !sediment layer 1 nitrification velocity (m/day)
                               !not used, replaced by KAPPNH4F and KAPPNH4S
                    KAPPNO3        !sediment layer 1 denitrification velocity
                               !(m/day), never used, replaced by KAPPNO3F,
                               !KAPPNO3F (integrated over layer 1 thickness)
       !WLong commented KAPP1CH4, it is duplicated with KAPPCH4
       !WLong   KAPP1CH4        !sediment layer 1 CH4 oxidation velocity (m/day)


      INTEGER ::    KWC, & 
                 IERR, &
                 ITEMP

    !WLong moved from wqm_modules.F
    CHARACTER(LEN=20) :: SSNAME(17)  
    
    !WLong moved from wqm_modules.F
     REAL(SP) ::     K0H1D,   &        !dissolved constituent reaction velocity of first  layer (m/d) (Wlong notes (3-32-7)(3-32-8)
                    K0H1P,   &        !particulate constituent reaction velocity of first layer (m/d) (always zero) (WLong notes (3-32-9)
                    K1H1D,   &      !dissolved constituent reaction velocity in first layer (m/d) (WLong notes (3-32-5)(3-32-6))
                    K1H1P,   &        !particulate constituent reaction velocity of first layer (m/d) (WLong notes (3-32-10))
                    K2H2D,   &        
                    K2H2P,   &
                    K3,      &    !dissolution rate of particulate biogenic silicate  PSISED (m/d)
                    J1,      &  !**** need to be careful with J1 here, better rename to JJ1, J1 is used else where
                    J2,      &  !**** need to be careful with J2 here, better rename to JJ2, J2 is used else where
                    KMC1,    &    !half reaction constant for total concentration in layer 1 (CT1) 
                                !(used for NH4, i.e. mgN/m^3, for other constituents, KMC1 is zero)
                    KL12,    &
                    KL12NOM, &    !KL12 after temperature regulation (m/day)
                    KHS_1,   &
                    KL12SO4
    


        !sediment POM concentration (all two layers ~= layer 2)
       REAL(SP) ::     PON1,          & !sediment PON G1 (mgN/m^3)    at new time step
                    PON1TM1,       & !sediment PON G1 (mgN/m^3) at previous time step
                    
                    PON2,          & !sediment PON G2 at new time step
                    PON2TM1,       & !sediment PON G2 (mgN/m^3) at previous time step
                    
                    PON3,          & !sediment PON G3 (mgN/m^3) at new time step
                    PON3TM1,       & !sediment PON G3 (mgN/m^3) at previous time step
                    
                    POC1,          & !sediment POC G1 (mgC/m^3) at new time step
                    POC1TM1,       & !sediment POC G1 (mgC/m^3) at previous time step
                    
                    POC2,          & !sediment POC G2 (mgC/m^3) at new time step
                    POC2TM1,       & !sediment POC G2 (mgC/m^3) at previous time step
                    
                    POC3,          & !sediment POC G3 (mgC/m^3) at new time step
                    POC3TM1,       & !sediment POC G3 (mgC/m^3) at previous time step
                    
                    POP1,          & !sediment POP G1 (mgP/m^3) at new time step
                    POP1TM1,       & !sediment POP G1 (mgP/m^3) at previous time step
                    
                    POP2,          & !sediment POP G2 (mgP/m^3) at new time step
                    POP2TM1,       & !sediment POP G2 (mgP/m^3) at previous time step
                    
                    POP3,          & !sediment POP G3 (mgP/m^3) at new time step
                    POP3TM1             !sediment POP G3 (mgP/m^3) at previous time step

  REAL(SP) ::         PSISED,        & !sediment particulate biogenic Si (mgSi/m^3)
                    PSISEDTM1        !sediment particulate biogenic Si (mgSi/m^3)
                    

  REAL(SP) ::         JNX,           & !depth integrated PON diagenesis rate (mgN/m^2/day)
					JNX_Layer1,           & !depth integrated PON diagenesis rate (mgN/m^2/day) in layer 1 from DON  B Clark APril 2017
                    JCX,           & !depth integrated POC diagenesis rate (mgC/m^2/day)
                    JPX    ,   &           !depth integrated POP diagenesis rate (mgP/m^2/day)
					JPX_LAYER1             !depth integrated POP diagenesis rate (mgP/m^2/day) from DOP remin in layer 1
  
  REAL(SP) ::         XJCNO3,        & !Sulfide/methhane oxidation flux due to NO3 used in denitrifcation  (gO2/m^2/d)
                    XJCNO31,         &  !Sulfide/methhane oxidation flux due to NO3 used in denitrifcation  (gO2/m^2/d) in layer 1
                    XJCNO32,       & !Sulfide/methhane oxidation flux due to NO3 used in denitrifcation  (gO2/m^2/d) in layer 2
                    XJCO2,         & !never used
                    XJC1             ! (gO2/m^2/day) cabon diagenesis flux converted to oxygen equivalents.
                                     !discounting organic matter used for denirification in layer
                                     !1 and layer 2                       
    !WLong moved from wqm_modules.F
    
        REAL(SP) :: NH40,     &        !mgN/m^3
                    NH41,     &        !mgN/m^3

                   NH41TM1,  &        !mgN/m^3        
                    NH4T1,    &        !mgN/m^3        
                    NH42,     &        !mgN/m^3
                    NH4T2,    &        !mgN/m^3
                    NH4T2TM1, &        !mgN/m^3
                    JNH4            !mgN/m^2/day
    
        REAL(SP) :: NO30,     &        !mgN/m^3
                    NO31,     &        !mgN/m^3
                    NO31TM1,  &        !mgN/m^3
                    NO3T1,    &        !mgN/m^3
                    NO32,     &        !mgN/m^3
                    NO3T2,    &        !mgN/m^3
                    NO3T2TM1, &        !mgN/m^3
                    JNO3            !mgN/m^2/day
                
       REAL(SP) ::     PO40,          &!overlaying water PO4 (mgP/m^3)
                    PO41,          &!sediment layer 1 PO4 (mgP/m^3)
                    PO42,          &!sediment layer 2 PO4 (mgP/m^3)
                    PO4T1,         &!sediment layer 1 PO4 (mgP/m^3)
                    PO4T2,         &!sediment layer 2 PO4 (mgP/m^3)
                    PO41TM1,       &!sediment layer 1 PO4 (mgP/m^3)
                    PO4T2TM1,       &!sediment layer 2 PO4 (mgP/m^3)
                    JPO4            !(mgP/m^2/day)
                
 !WLong moved these here from wqm_modules.F

       REAL(SP) ::     SI0,           &!overlying water Si. (mgSi/m^3)
                    SI1,           &!sediment layer 1 Si (mgSi/m^3)
                    SI2,           &!sediment layer 2 Si (mgSi/m^3)
                    SIT1,          &!sediment layer 1 Si (mgSi/m^3)
                    SIT2,          &!sediment layer 2 Si (mgSi/m^3)
                    SI1TM1,        &!sediment layer 1 Si (mgSi/m^3)
                    SIT2TM1,       &!sediment layer 2 Si (mgSi/m^3)
                    JSI             !(mgSi/m^2/day)
                    
        REAL(SP) :: HS0,             &!sediment overlying water HS- concentration (gO2/m^3)
                    HS1,           &!sediment 1st layer HS (gO2/m^3)
                    HS2,           &!sediment 2nd layer HS (gO2/m^3)
                    HST1,          &!sediment 1st layer HS (gO2/m^3)
                    HST2,          &!sediment 2nd layer HS (gO2/m^3)
                    HS1TM1,        &!sediment 1st layer HS (gO2/m^3)
                    HST2TM1,       &!sediment 2nd layer HS (gO2/m^3)
                    JHS                !(gO2/m^2/day)

        REAL(SP) :: CH40,     &    !(gO2/m^3)
                    CH41,     & !(gO2/m^3)
                    CH41TM1,  &    !(gO2/m^3)
                    CH4T1,    & !(gO2/m^3)
                    CH42,     &    !(gO2/m^3)
                    CH42AV,   &    !(gO2/m^3)                
                    CH4T2,    &    !(gO2/m^3)    
                    CH4T2AV,  & !(gO2/m^3)
                    CH4T2TM1, & !(gO2/m^3)
                    JCH4G,    & !(gO2/m^2/day)
                    JCH4AQ,   & !(gO2/m^2/day)
                    JCH4        !(gO2/m^2/day)
    
        REAL(SP) :: SO40,     &!mgO2/L (gO2/m^2/day)
                    SO40MG,   &!overlying SO4 concentration in mgSO4/L (1 mole SO4 = 98g = 32+16*4) where 32 is S, 16 is O
                    SO41,     &!(gO2/m^3)
                    SO41TM1,  &!(gO2/m^3)
                    SO4T1,    &!(gO2/m^3)
                    SO42,     &!(gO2/m^3)
                    SO42AV,   &!(gO2/m^3)
                    SO4T2,    &!(gO2/m^3)
                    SO4T2AV,  &!(gO2/m^3)
                    SO4T2TM1, &!(gO2/m^3)
                    JSO4,     &!(gO2/m^3)
                    SO4T2SAV   !never used
            
    REAL(SP) :: JO2NH4,   & !gO2/m2/day
                CH4T2SAV, & !never used
                JCH4GASO2,& !never used
                JGAS        !never used
                    
    REAL(SP) :: KPOP1,    &      !G1 POP diagenesis rate (1/day)
                KPOP2,    &      !G2 POP diagenesis rate (1/day)
                KPOP3,    &      !G3 POP diagenesis rate (1/day)
                KPON1,    &      !G1 PON diagenesis rate (1/day)
                KPON2,    &      !G2 PON diagenesis rate (1/day)
                KPON3,    &      !G3 PON diagenesis rate (1/day) 
                KPOC1,    &      !G1 POC diagenesis rate (1/day)
                KPOC2,    &      !G2 POC diagenesis rate (1/day)
                KPOC3            !G3 POC diagenesis rate (1/day)


    REAL(SP) :: ISEDMN,   &!total amount of nitrogen in sediments (PON+NH4+NO3) (kgN) (previous time step step)
                ISEDMP,   &!total amount of phosphorus in sediments (POP+PO4) (kgP) (previous time step)
                ISEDMC       !total amount of carbon in sediments (POC) (previous time step)

   REAL(SP) ::     SEDMN,      & !total sediment nitrogen  (kgN) (new time step) 
                SEDMP,      & !new time step total sediment phosphorus (kgP) (new time step)
                SEDMC         !total sediment carbon  (kgC) (new time step)

    
       REAL(SP) ::      WSSINETMMD,     &!Inorganic suspended solids settling velocity (mm/day)
                    WLPOMNETMMD,   &!Labile POM settling velocity (mm/day)
                    WRPOMNETMMD,   &!Refractory POM settling velocity (mm/day)
                    WB1NETMMD,     &!Algae 1 settling rate (mm/day)
                    WB2NETMMD,     &!Algae 2 settling rate (mm/day)
                    WB3NETMMD,     &!Algae 3 settling rate (mm/day)
                    WPSINETMMD      !Particulate biogenic silica settling rate (mm/day) 

    REAL(SP) ::     STP20,         &!sediment temperature subtracted by 20 degC (degC)
                    DF,         &!dissolved fraction of a constituent 
                    PO4AVL,     &!available PO4 for benthic algae growth (=DF*PO4) in bottom layer of water column (gP/L ?)
                    BFOR,         &!overlying water oxygen saturation coefficient for particle mixing (limit to bioturbation) (non-dimensional)
                    W12MIN         !minimum mixing speed of layer 1 and 2 (cm/day?)
                
    REAL(SP)     :: LOGICT
                
	
	! XKPOCn is specified as a reaction rate, 1/d, and gets multiplied by the layer 2 depth H2 to yield a reaction velocity for numerical integration
    REAL(SP)     :: XKPOC1, &!diagenesis reaction velocity (layer integrated ~ kpoc1*H2) (m/d) of POC1 with temperature control
                    XKPOC2, &!diagenesis reaction velocity (layer integrated ~ kpoc2*H2) (m/d) of POC2 with temperature control
                    XKPOC3, &!diagenesis reaction velocity (layer integrated ~ kpoc3*H2) (m/d) of POC3 with temperature control
                    XKPOP1, &!diagenesis reaction velocity (layer integrated ~ kpop1*H2) (m/d) of POP1 with temperature control
                    XKPOP2, &!diagenesis reaction velocity (layer integrated ~ kpop2*H2) (m/d) of POP2 with temperature control
                    XKPOP3, &!diagenesis reaction velocity (layer integrated ~ kpop3*H3) (m/d) of POP3 with temperature control
                    XKPON1, &!diagenesis reaction velocity (layer integrated ~ kpon1*H2) (m/d) of PON1 with temperature control
                    XKPON2, &!diagenesis reaction velocity (layer integrated ~ kpon2*H2) (m/d) of PON2 with temperature control
                    XKPON3   !diagenesis reaction velocity (layer integrated ~ kpon3*H2) (m/d) of PON3 with temperature control
                    

    REAL(SP) ::        TEMP_SED
    
    REAL(SP) ::     SFEED_TMP_1,SFEED_TMP_2,SFEED_TMP_3 !3 types of suspension feeders
                
    REAL(SP) ::        DLTS_H2  !sedimen time step (day) divided by second sediment layer thickness (m), final unit: 1/(m*d)
                
    REAL(SP)     :: ERROR,    &
                    A1,        &
                    AO2N        !oxygen to nitogren ratio for denitrification (WLong (4-45)*(4-46)) (gO2/gN)

    REAL(SP)     :: FP1SO4,    &
                    FP2SO4,    &
                    HS2AV,    &
                    XJ2,    &
                    XJ2CH4,    &        !Source of methane in layer 2 (DMD book eqn(10.34))
                    X1J2

    REAL(SP)     :: PF,        &!
                    PPO4        !particulate iorganic ortho phosphorate 

    REAL(SP)     :: FLUXHS,        &
                    FLUXHSCH4,    &
                    VJCH4G

    !sediment POM and particulate biogenic silica
    REAL(SP),ALLOCATABLE,DIMENSION(:):: PON1TM1S,        &!sediment G1 PON (mgN/m^3)
                                        PON2TM1S,        &!sediment G2 PON (mgN/m^3)
                                        PON3TM1S,        &!sediment G3 PON (mgN/m^3)
                                        POC1TM1S,        &!sediment G1 POC (mgC/m^3)
                                        POC2TM1S,        &!sediment G2 POC (mgC/m^3)
                                        POC3TM1S,        &!sediment G3 POC (mgC/m^3)
                                        POP1TM1S,        &!sediment G1 POP (mgP/m^3)
                                        POP2TM1S,        &!sediment G2 POP (mgP/m^3)
                                        POP3TM1S,        &!sediment G3 POP (mgP/m^3)
                                        PSISEDTM1S           !sediment particulate organic silica  (mgSi/m^3)

    !sediment inorganic chemicals.
    REAL(SP),ALLOCATABLE,DIMENSION(:):: NH41TM1S,  &!sediment disslolved NH4 (1st layer) (mgN/m^3)
                                        NH42TM1S,  &!sediment dissolved NH4 (2nd layer) (mgN/m^3)
                                        NH4T2TM1S, &!sediment total NH4 (2nd layer) (mgN/m^3)
                                        
                                        NO31TM1S,  &!sediment dissolved NO3 (1st layer) (mgN/m^3)
                                        NO32TM1S,  &!sediment dissolved NO3 (2nd layer) (mgN/m^3)
                                        NO3T2TM1S, &!sediment total NO3 (2nd layer) (mgN/m^3)
                                        
                                        
                                        PO41TM1S,  &!sediment dissolved PO4 (1st layer) (mgP/m^3)
                                        PO42TM1S,  &!sediment dissolved PO4 (2nd layer) (mgP/m^3)
                                        PO4T2TM1S, &!sediment total PO4 (2nd layer) (mgP/m^3)
                                        
                                        
                                        SI1TM1S,   &!sediment dissolved silica (1st layer) (mgSi/m^3)
                                        SI2TM1S,   &!sediment dissolved silica (2nd layer) (mgSi/m^3)
                                        SIT2TM1S,  &!sediment total silica (2nd layer) (mgSi/m^3)    
                                        
                                        HS1TM1S,   &!sediment dissolved HS (1st layer) (gO2/m^3)
                                        HS2TM1S,   &!sediment dissolved HS (2nd layer) (gO2/m^3)
                                        HST2TM1S,  &!sediment total HS (2nd layer) (gO2/m^3)
                                        
                                        
                                        CH41TM1S,  &!sediment dissolved CH4 (1st layer) (gO2/m^3)
                                        CH42TM1S,  &!sediment dissolved CH4 (2nd layer) (gO2/m^3)
                                        CH4T2TM1S, &!sediment total CH4 (2nd layer) (gO2/m^3)
                                        
                                        !SO41TM1S,  &!sediment dissolved SO4 (1st layer) (gO2/m^3)
                                        SO4T2TM1S   !sediment total SO4 (2nd layer) (gO2/m^3)
                                        
                                           
    !WLong moved the following here from wqm_modules.F                                        
    !benthic stress calculation
    REAL(SP),ALLOCATABLE,DIMENSION(:):: BENSTRTM1S,        &     !sediment sediment benthic stres
                                        BFORMAXS,         &     !used for calculating benthic stress accumulation
                                        ISWBENS             !Array indicating whether benthic stress accumulation is
                                                            !set to zero
                                                            !If sediment temperature > TEMPBEN, then
                                                            !    yes, it is set to zero   and ISWBENS=1, BFORMAX=0
                                                            !otherwise ISWBENS is 0,
                                                            !      BFOR=BFORMAX                                            
                                        
    REAL(SP), ALLOCATABLE, DIMENSION(:) :: SODTM1S,  & !Sediment oxygen demand in previous time step (mgO2/m^2/d)
                                           JNH4TM1S, & !dissolved NH4 flux to water column (mgN/m^2/d) in previous time step
                                           JNO3TM1S, & !NO3 flux to water column (mgN/m^2/d) in previous time step
                                           JPO4TM1S, & !dissolved PO4 flux to water colmn (mgP/m^2/d)
                                           JCH4TM1S, & !dissolved CH4 flux to water column (gO2/m^2/d) in previous time step
                                           JCH4GTM1S,& !CH4 gas flux to water column (gO2/m^2/d)
                                           JHSTM1S,     & !dissolved HS flux to water column (gO2/m^2/d)
                                           JSITM1S       !dissolved Si flux to water column (gO2/m^2/d)
! #if defined(BUDGET)										   
    REAL(SP),ALLOCATABLE,DIMENSION(:) :: JNH4_OUt, JNO3_OUT ! Integrated variables for calculating budgets ,
															!are integrated in time and space in wqm_main
! #endif 
    !WLong moved here from wqm_modules.F
    !sediment burial fluxes of C, N, P
     REAL(SP),ALLOCATABLE,DIMENSION(:) ::     BURIALC,    & !Burial flux of carbon    (mgC/m^2/day)
                                                            !(positive leaving sediment)
                                            BURIALN,    & !Burial flux of nitrogen    (mgN/m^2/day)
                                                            !(positive leaving sediment)
                                            BURIALP       !Burial flux of phosphorus (mgP/m^2/day)
                                                            !(positive leaving sediments)
    
    !WLong moved here from wqm_modules.F
    REAL(SP), ALLOCATABLE, DIMENSION (:) :: DIAGENC     !carbon diagenesis flux (gC/m^2/day)


    !WLong moved these to wqm_sed.F.
     !look up tables for temperathre control used for interpolation at given
     !temperature
    REAL(SP),DIMENSION(350) :: ZHTADP,          &  !not used
                              ZHTADD,           &  !not used
                              ZHTANH4F,         &  !nitrification velocity after applying temperature control (m/day) for fresh water.
                              ZHTANH4S,         &  !nitrification velocity after applying temperature control (m/day) for salt water
                              ZHTANO3F,         &  !layer 1 denitrification velocity after applying temperature control (m/day) for fresh water
                              ZHTANO3S,         &  !layer 1 denitrification velocity (m/day) after applying temperature control.
                              ZHTAK2NO3,         &  !layer 2 denitrifcationvelocity after applying temperature control (m/day)
                              ZHTA1HSD,         &  !dissolved sulfide oxidationrate (m/day) after temperature control
                              ZHTA1HSP,         &  !particulate sulfide oxidation velocity (m/day) after applying temperature control
                              ZHTASI,           &  !particulate biogenic silica dissolution velocity after applying temperature effect (1/day)
                              ZL12NOM,          &  !temperature control of sediment porewater mixing velocity across layer 1 and layer 2 interface( unitless)
                              ZW12NOM,          &  !sediment particle mixing rate across layer 1, 2 interface after applying temperature control (m/day)
                              ZHTAPON1,         &  !diagenesis G1 PON reaction rate look up table (1/day) after applying temperature control
                              ZHTAPON2,         &  !diagenesis G2 PON reaction rate (1/day) after applying temperature control
                              ZHTAPON3,         &  !diagenesis G3 PON reaction rate (1/day) after applying temperature control
                              ZHTAPOC1,         &  !diagenesis G1 POC reaction rate (1/day) after applying temperature control
                              ZHTAPOC2,         &  !sediment diagenesis G2 POC reaction rate look up table (1/day) after applying temperature control.
                              ZHTAPOC3,         &  !sediment diagenesis G3 POC reaction rate (1/day) after applying temperature control
                              ZHTAPOP1,         &  !sediment diagenesis G1 POP reaction rate (1/day) after applying temperature control
                              ZHTAPOP2,         &  !sediment diagenesis G2 POP reaction rate (1/day) after applying temperature control.
                              ZHTAPOP3,         &  !sediment diagenesis G3 POP reaction rate (1/day) after applying temperature control
                              ZHTACH4              !sediment methane oxidation rate look up table (m/day) after applying temperature effect

    !WLong moved this to wqm_sed.F.
    !POM diagenesis rate temperature theta.
    !used in setting the temperature control look up tables
    !WLong: I'm surprised that these are not in the input controls.
    !       Well these are actually setup as equivalents to DPTHTA etc ****
    !       I have now set this to the corresponding values in SED_READ() subroutine
    !       and removed the equivalence declaration originally in wqm_modules.F - WLong,Feb 7/2014
    REAL(SP) ::  THTAPOP1,  &  !sediment POP diagenesis reaction temperature theta
                           !(G1)
             THTAPOP2,  &  !(G2)
             THTAPOP3,  &  !(G3)
             THTAPON1,  &  !sediment PON diagenesis reaction temperature theta
                           !(G1)
             THTAPON2,  &  !(G2)
             THTAPON3,  &  !(G3)
             THTAPOC1,  &  !sediment POC diagenesis reaction temperature theta
                           !(G1)
             THTAPOC2,  &  !(G2)
             THTAPOC3      !(G3)

    !WLong moved this here from wqm_modules.F
    !time averaged G3 particulate organic flux (mg/m^2/day)
    REAL(SP),ALLOCATABLE,DIMENSION(:) ::AG3CFL,  &     !G3 sediment POC flux
                                                    !(mgC/m^2/day)
                                        AG3NFL,  &     !G3 sediment PON flux
                                                    !(mgN/m^2/day)
                                        AG3PFL         !G3 sediment POP flux
                                                    !(mgP/m^2/day)

    !WLong moved this here from wqm_modules.F
    !time averaged sediment temperature (degC)
    REAL(SP),ALLOCATABLE,DIMENSION(:) :: ASDTMP !time averaged sediment
                                             !temperature  (degC)

    !Benthic algae growth limitation and kinetic fluxes
    !WLong moved this to wqm_sed.F.
    !benthic algae nutrient limitation and kinetic fluxes
    
!    REAL(SP),ALLOCATABLE,DIMENSION(:) :: !FIB,      &  !benthic algae growth light
                                                  !limitation (dimensionless)
!Wlong                                     NLB,      &  !benthic algae growth nitrogen
                                                  !limitation ***need to check
                                                  !NH4 limitation ***
                                                  !(dimensionless)
!WLong                                     PLB,      &  !benthic algae growth
                                                  !phosphorus limitation
                                                  !(dimensionless)
!WLong                                     NPPB,     &  !net primary production of
                                                  !benthic algae (gC/m^2/day ?)
!WLong                                     BBM,      &  !benthic algae biomass
                                                  !(gC/m^2 ?)
!WLong                                     BLITE,    &  !benthic algae light (light
                                                  !at bottom of water column)
                                                  !(W/m^2)
!WLong                                     BMB,      &  !benthic algae basal
                                                  !metabolism  after
                                                  !temperature control (1/day)
!WLong                                     PB,       &  !primary production rate of
                                                  !benthic algae
!Wlong                                     PRB,      &  !predation rate of benthic
                                                  !algae (1/day)
!WLong                                     BANH4,    &  !adjustment to benthic NH4
                                                  !flux gN/m^2/day by benthic algae
!WLong                                     BANO3,    &  !adjustment to benthic NO3
                                                  !flux gN/m^2/day by benthic
                                                  !algae
!WLOng                                     BAPO4,    &  !adjustment to benthic PO4
                                                  !flux gP/m^2/day by benthic
                                                  !algae
!WLong                                     BADOC,    &  !adjustment to
                                                  !benthic DOC flux by benthic
                                                  !algae (gC/m^2/day)
!WLong                                     BADO,     &  !adjustment to
                                                  !benthic DOXG flux by benthic
                                                  !algae (gO2/m^2/day)
!WLong                                     BAPOC,    &  !adjustment to benthic POC
                                                  !flux by benthic algae
                                                  !(gC/m^2/day)
!Wlong                                     BAPON,    &  !adjustment to benthic PON
                                                  !flux by benthic algae
                                                  !(gN/m^2/day)
!WLong                                     BAPOP        !adjustment to benthic POP
                                                  !flux by benthic algae
                                                  !(gP/m^2/day)

  REAL(SP) ::     PIE1,          & !sediment 1st layer partition coefficient (L/kg)
                                    !temporary variable
                PIE2,          & !sediment 2nd layer partition coefficient (L/kg)
                                    !temporary variable
                W12,           & !particle mixing velocity (cm/yr) between Layer 1
                                    !and layer 2
                TEMPD,         & !sediment temperature (degC) temporary variable
                O20,           & !overlying water concentration (mgO2/L)
                CH4SAT,        & !saturation concentration of CH4 in sediment pore
                                    !water (mgO2/L)
                SAL              !overlying water salinity  (ppt)


  REAL(SP) ::     XAPPNH4,    & !nitrification velocity after applying temp control (m/day)
                XAPP1HSD,     & !sediment dissolved sulfide reaction rate after applying temperature effects (m/day)
                XAPP1HSP,     & !sediment particulate sulfide reaction rate after applying temperature effects (m/day)
                XAPP1NO3,   & !sediment denitrification rate after applying temmperature effects (m/day)
                XK2NO3,     & !sediment denitrification rate after applying temperature effects and multiplied by layer thickness H2 (for layer 2) (m/day)
                XKSI,       & !sediment particulate bigenic silica dissolution rate after applying temperature control (m/day), DMD eqn 7.7)
                XAPPCH4,    & !sediment methane oxidation rate after applying temperature effects (m/day)
                TEMP20,     & !temperature (degC) - 20
                TEMP202,    & !TEMP20/2.0
                FD1,        & !dissolved fraction of layer 1 FD1 = 1./(1.+M1*PIE1)
                FP1,        & !particulate fraction of layer 1 FP1=M1*PIE1/(1.+M1*PIE1)
                FD2,        & !dissolved fraction of layer 2 FD2 = 1./(1.+M2*PIE2)
                FP2,        & !particulate fraction of layer 2 FP2 =M2*PIE2/(1.+M2*PIE2)
                SOD,        & !sediment oxygen demand (gO2/m^2/day)
                CSOD,       & !chemical (? or carbonaceous ?)oxygem demand (gO2/m^2/day)
                S,          & !SOD/O20 (m/day) sediment powerwater difffusion velocity (m/day)
                W12NOM,     & !particle mixing velcoity (W12*) (cm/yr), DMD eqn. (13.2)
                HSO4,       & !depth of SO4 reaction (cm)
                DDSO4,      & !diffusion rate of porewater SO4 (m/day)
                CSODHS,     & !chemical oxygen demand due to HS oxidation (gO2/m^2/day)
                CSODCH4,    & !chemical oxygen demand due to CH4 oxidation. (gO2/m^2/day)
                BENSTR,     & !benthic stress
                BENSTRS,    & !benthic stress    (never used)
                BENSTRTM1,  & !benthic stress at previous time step
                SODTM1,     & !last time step SOD
                ISWBEN,     & !flag for re-setting benthic stress to zero
                BFORMAX       !benthic stress limitation temporary variable (dimensionless)
        
   REAL(SP) ::     ZHTANH4,    & !look up table value for nitrification rate after applying temperature control
                                    !not used, replaced by ZHTANH4S, ZHTANH4F
                ZHTANO3       !loop up table value for denitrification rate after
                                    !applying temperature control
                                    !not used, replaced by ZHTANO3F. ZHTANO3S

   !WLong moved here from wqm_modules.F    
       REAL(SP) :: TINTIM !time (days) for steady state simulation of sediments
            
    !global variables fror initial condition and also output
    REAL(SP), ALLOCATABLE ::    CTEMP_GL(:),    &
                            PO4T2TM1S_GL(:),    &
                            NH4T2TM1S_GL(:),    &
                            NO3T2TM1S_GL(:),    &
                             HST2TM1S_GL(:),     &
                            CH4T2TM1S_GL(:),    &
                             CH41TM1S_GL(:),    &
                            SO4T2TM1S_GL(:),     &
                             SIT2TM1S_GL(:),    &
                             BENSTRTM1S_GL(:),    &
                             SODTM1S_GL(:),     &
                                CPOP_GL(:,:),    &
                               CPON_GL(:,:),    &
                               CPOC_GL(:,:),    &
                               CPOS_GL(:)
    
!used for outputing sediment flux variables into history and station files
        
        REAL(SP), ALLOCATABLE ::    JPOC_GL(:,:),    &    
                                    JPON_GL(:,:),    &    
                                    JPOP_GL(:,:),    &    
                                    JPOS_GL(:),    &
                                    !O20_GL(:),        &!given by C2_GL(:,:,27)
                                    !D_GL(:),        &!already defined
                                    !T_GL(:,:),        &!already defined
                                    !NH40_GL(:),    &!given by C2_GL
                                    !NO30_GL(:)),    &!given by C2_GL
                                    !SI0_GL(:)),    &!given by C2_GL
                                    CH40_GL(:),        &
                                    !SAL_GL(:),        &!given by C2_GL
                                    !SOD_GL(:),        &!given by SODTM1S_GL
                                    JNH4_GL(:),        &    
                                    JNO3_GL(:),        &    
                                    BENDEN_GL(:),    &    
                                    JCH4_GL(:),        &    
                                    JCH4G_GL(:),    &    
                                    JHS_GL(:),        &    
                                    JPO4_GL(:),        &    
                                    JSI_GL(:),        &    
                                    NH41_GL(:),        &    
                                    NH42_GL(:),        &    
                                    NO31_GL(:),        &    
                                    NO32_GL(:),        &    
                                    PO41_GL(:),        &    
                                    PO42_GL(:),        &    
                                    SI1_GL(:),        &    
                                    SI2_GL(:),        &    
                                    !CH41_GL(:),    &!given by CH41TM1S_GL
                                    CH42_GL(:),        &    
                                    HS1_GL(:),        &    
                                    HS2_GL(:),        &    
                                    !POC1_GL(:),    &!given by CPOC_GL
                                    !POC2_GL(:),    &    
                                    !POC3_GL(:),    &    
                                    !PON1_GL(:),    &!given by CPON_GL
                                    !PON2_GL(:),    &    
                                    !PON3_GL(:),    &    
                                    !POP1_GL(:),    &!given by CPOP_GL
                                    !POP2_GL(:),    &    
                                    !POP3_GL(:),    &    
                                    !PSISED_GL(:),    &!given by CPOS_GL
                                    HSED1_GL(:)!,    &    
                                    !BENSTR_GL(:)     !given by BENSTRTM1S_GL
                                         
                                 
    !diagnostic flux from P,N,C, Si for sediment module
    REAL(SP), ALLOCATABLE,DIMENSION (:) :: DIAGP,        &!not used
                                             DIAGN,        &!    
                                           DIAGC,        &!not used
                                           DIAGS        !!not used

    REAL(SP), ALLOCATABLE,DIMENSION (:) :: MTVEL    !mass transfer velocity (m/sec) (=SOD/O20)
                                 
    !The following needs to be moved to SED_INIT 
    REAL(SP),ALLOCATABLE                ::     CPOP(:,:),    &!Sediment POP (mgP/m^3)
                                            CPON(:,:),    &!Sediment PON (mgN/m^3)
                                            CPOC(:,:),  &!Sediment POC (mgC/m^3)
                                            CPOS(:)         !Sediment particulate organic silica (mgSi/m^3)
    
    REAL(SP),ALLOCATABLE,DIMENSION(:) ::     CPO4,   &!Sediment inorganic phosphorus (mgP/m^3) (PO4)
                                            CNO3,   &!sediment NO3 concentration (2nd layer) (mgC/m^3)
                                            CNH4,   &!sediment NH4 concentration (mgN/m^3)    
                                            CSI,    &!Silicate concentration (mgSi/m^3)                                            
                                            CCH4,   &!sediment CH4 concentration (gO2/m^3)
                                            CSO4,   &!sediment SO4 concentration (gO2/m^3)
                                            CHS      !Sediment HS concentration (gO2/m^3)


    !Initial conditions of sediment constitutients (moved here from mod_wqminit.F) for individual cell
    REAL(SP)::     CPOPI(3),    &!mgP/m^3
                CPONI(3),     &!mgN/m^3
                CPOCI(3),    &!mgC/m^3
                CPOSI,        &!mgSi/m^3    
                PO41TI,        &!mgP/m^3 Added by WLong
                PO4T2I,        &!mgP/m^3
                NO31TI,        &!mgN/m^3 Added by WLong
                NO3T2I,        &!mgN/m^3
                NH41TI,        &!mgN/m^3 Added by WLong
                NH4T2I,        &!mgN/m^3
                CH4T2I,        &!gO2/m^3
                CH41TI,        &!gO2/m^3
                SO41TI,        &!gO2/m^3 Added by WLong
                SO4T2I,        &!gO2/m^3
                HS1TI,        &!gO2/m^3 Added by WLong
                HST2I,         &!gO2/m^3
                SI1TI,        &!mgSi/m^3 Addec by WLong
                SIT2I,         &!mgSi/m^3
                CTEMPI,        &!degC
                BENSTI,        &!dimensionless
                SODI         !Initial SOD mgO2/m^2/day
                
                                            
    REAL(SP),ALLOCATABLE:: JPOP(:,:), &!POP flux (mgP/m^2/day) positive into sediments
                           JPON(:,:), &!PON flux (mgN/m^2/day) positive into sediments
                           JPOC(:,:), &!POC flux (mgC/m^2/day) positive into sediments
                           JPOS(:)       !POSi settlong flux(mgSi/m^2/day) positive into sediments
! #if defined(BUDGET)						   
	REAL(SP),ALLOCATABLE, DIMENSION(:,:):: JPOC_out,JPON_out, JPOP_out ! Added by B Clark for Budget Calculation
! #endif                           
    REAL(SP),ALLOCATABLE:: JPOPaccum(:,:), &!accumulative POP flux (mgP/m^2/day) positive into sediments
                           JPONaccum(:,:), &!accumulative PON flux (mgN/m^2/day) positive into sediments
                           JPOCaccum(:,:), &!accumulative POC flux (mgC/m^2/day) positive into sediments
                           JPOSaccum(:)     !accumulative POSi settlong flux(mgSi/m^2/day) positive into sediments
                           
                           
    REAL(SP) POC1R  !mgC/gSediments  reference POC concentration in sediments for calculating W12NOM (DMD page 277)
                    !=0.1mgC/gSediments ~ 100000*M2 mgC/m^3 ( if M2 has unit kg/L for particle concentration of layer 2)
                    !WLong notes eqn(3-26) (DMD book page 278, table 13.1, eqn(13.1))

    INTEGER :: IDEBUG_SED
    
    CHARACTER(LEN=72), DIMENSION(NFLP) ::WCLFN  !sediment overlying water column forcing filenames
    
    !
    !variables related to records in WCLFN which provids overlying water column condition as forcing to 
    !sediment flux model 
    !
    !WLong : will need to make this as field arrays for all TCE cells rather then a domain wide uniform value
    !
    
    REAL(SP) :: NXJDAY
    REAL(SP) :: NXWCL
    REAL(SP) :: JDAY_R1                    !JDAY of 1st record kept in memory                    (day)
    REAL(SP) :: JDAY_R2                 !JDAY of 2nd record kept in memory                    (day)
        
    REAL(SP) :: JCIN_R1,    JCIN_R2        !1st and 2nd record of Jcin                         (gO2/m^2/d)
    REAL(SP) :: JNIN_R1,    JNIN_R2        !1st and 2nd record of Jnin                         (gN/m^2/d)
    REAL(SP) :: JPIN_R1,    JPIN_R2        !1st and 2nd record of Jpin                         (gP/m^2/d)
    REAL(SP) :: JSIN_R1,    JSIN_R2        !1st and 2nd record of Jsin                         (gSi/m^2/d)    
        
    REAL(SP) :: O20_R1,        O20_R2        !1st and 2nd record of O20                             (mgO2/L)
    REAL(SP) :: D_R1,        D_R2        !1st and 2nd record of depth                         (m)
    REAL(SP) :: TW_R1,        TW_R2        !1st and 2nd record of overlying water temperature     (degC)
        
    REAL(SP) :: NH30_R1,    NH30_R2        !1st and 2nd record of overlying ammonia             (mgN/L)
    REAL(SP) :: NO30_R1,    NO30_R2        !1st and 2nd record of overlying nitrate             (mgN/L)
    REAL(SP) :: PO40_R1,    PO40_R2        !1st and 2nd record of overlying phosphate             (mgP/L)
    REAL(SP) :: SIAT0_R1,    SIAT0_R2    !1st and 2nd record of overlying total silicate     (mgSi/L)
    REAL(SP) :: CH40_R1,    CH40_R2        !1st and 2nd record of overlying water CH4             (mgO2/L)
    REAL(SP) :: SALT0_R1,    SALT0_R2    !1st and 2nd record of overlying water salinity     (psu)
    
	REAL(SP) :: SODMIN, SODMAX   !min and max value of SOD 
   CONTAINS

   !subroutines:
        !subroutine SED_INIT()
        !subroutine SED_ALLOC()
        !subroutine SED_DEALLOC()
        !subroutine SED_READ()
        !subroutine POM_ACCUMUL
        !subroutine SED_CALC()
        !subroutine SED_DIAGENESIS_G3()
   
        !subroutine SEDTSFNL()
        !subroutine SEDSSFNL()

   !functions:
        !function SEDF()
        !function ZBRENT()
   
   
   
!********************************************************************************
!**                    S U B R O U T I N E   SED_INIT                          **
!********************************************************************************
  SUBROUTINE SED_INIT
 
        
    !WLong moved from wqm_modules.F
    SSNAME(1)  = 'Sediment Temperature'
    SSNAME(2)  = 'Sediment POP        '
    SSNAME(3)  = 'Sediment PON        '
    SSNAME(4)  = 'Sediment POC        '
    SSNAME(5)  = 'Sediment PBS        '
    SSNAME(6)  = 'Sediment PO4        '
    SSNAME(7)  = 'Sediment NH4        '
    SSNAME(8)  = 'Sediment NO3        '
    SSNAME(9)  = 'Sediment HS         '
    SSNAME(10) = 'Sediment CH4        '
    SSNAME(11) = 'Sediment CH4        '
    SSNAME(12) = 'Sediment SO4        '
    SSNAME(13) = 'Sediment DSIL       '
    SSNAME(14) = 'Benthic Stress      '
    SSNAME(15) = 'Benthic Algae       '
    SSNAME(16) = 'Deposit Feeders     '
    SSNAME(17) = 'Suspension Feeders  '

    !allocate variables related to sediment diagenesis
    CALL SED_ALLOC
    
    !initialize time varying overlying water colum input related variables
    ITVWCLF=0        !Assume no input file is needed 
    NTVWCLF=0        !assume no time varying overlying water column condition input file    
    WCLPTR=0        !pointer in WCLFN() array
    WCLFDAYS=0.0    !lenght of file set to zero

    CALL SED_INIT_ICI  !get the input initial conditions first
    
  END SUBROUTINE SED_INIT

!********************************************************************************
!**                    S U B R O U T I N E   SED_ALLOC                         **
!********************************************************************************
  SUBROUTINE SED_ALLOC
  
    !USE MOD_SED_DF_EXCHANGE_VARS, ONLY: &
    !    POC1TM1S_SED_DF,               &
    !    POC2TM1S_SED_DF

     !WLong moved here from wqm_modules.F
     ALLOCATE(FRPOP(MTLOC,3));            FRPOP = 0.0
     ALLOCATE(FRPON(MTLOC,3));            FRPON = 0.0
     ALLOCATE(FRPOC(MTLOC,3));            FRPOC = 0.0

   !   !WLong moved here from wqm_modules.F
   !   ALLOCATE(WSSBNET(MTLOC));            WSSBNET = 0.0
   !   ALLOCATE(WSLBNET(MTLOC));            WSLBNET = 0.0
   !   ALLOCATE(WSRBNET(MTLOC));            WSRBNET = 0.0
   !   ALLOCATE(WS1BNET(MTLOC));            WS1BNET = 0.0
   !   ALLOCATE(WS2BNET(MTLOC));            WS2BNET = 0.0
   !   ALLOCATE(WS3BNET(MTLOC));            WS3BNET = 0.0
   !   ALLOCATE(WSUBNET(MTLOC));            WSUBNET = 0.0
  

     !WLong moved here from wqm_modules.F
     ALLOCATE(VSED(MTLOC));              VSED   = 0.0
     ALLOCATE(VPMIX(MTLOC));             VPMIX  = 0.0
     ALLOCATE(VDMIX(MTLOC));             VDMIX  = 0.0

     !WLong moved here from wqm_modules.F
     ALLOCATE(PON1TM1S(MTLOC));          PON1TM1S = 0.0
     ALLOCATE(PON2TM1S(MTLOC));          PON2TM1S = 0.0
     ALLOCATE(PON3TM1S(MTLOC));          PON3TM1S = 0.0
     ALLOCATE(POC1TM1S(MTLOC));          POC1TM1S =0.0
     ALLOCATE(POC2TM1S(MTLOC));          POC2TM1S = 0.0
     ALLOCATE(POC3TM1S(MTLOC));          POC3TM1S = 0.0
     ALLOCATE(POP1TM1S(MTLOC));          POP1TM1S = 0.0
     ALLOCATE(POP2TM1S(MTLOC));          POP2TM1S = 0.0
     ALLOCATE(POP3TM1S(MTLOC));          POP3TM1S = 0.0
     ALLOCATE(PSISEDTM1S(MTLOC));           PSISEDTM1S  = 0.0
    
     !WLong added POC1TM1S_SED_DF
!     ALLOCATE(POC1TM1S_SED_DF(MTLOC));    POC1TM1S_SED_DF=0.0
!     ALLOCATE(POC2TM1S_SED_DF(MTLOC));  POC2TM1S_SED_DF=0.0
     
     !WLong moved here from wqm_modules.F
     ALLOCATE(BENSTRTM1S(MTLOC));         BENSTRTM1S = 0.0
     ALLOCATE(BFORMAXS(MTLOC));           BFORMAXS = 0.0
     ALLOCATE(ISWBENS(MTLOC));            ISWBENS  = 0.0
    
     ALLOCATE(SODTM1S(MTLOC));            SODTM1S=0.0
     ALLOCATE(JNH4TM1S(MTLOC));           JNH4TM1S=0.0
     ALLOCATE(JNO3TM1S(MTLOC));           JNO3TM1S=0.0
     ALLOCATE(JPO4TM1S(MTLOC));           JPO4TM1S=0.0
     ALLOCATE(JSITM1S(MTLOC));            JSITM1S=0.0
     ALLOCATE(JCH4TM1S(MTLOC));           JCH4TM1S=0.0
     ALLOCATE(JCH4GTM1S(MTLOC));          JCH4GTM1S=0.0
     ALLOCATE(JHSTM1S(MTLOC));            JHSTM1S=0.0
#if defined(BUDGET)    
     ALLOCATE(JNH4_out(MTLOC));           JNH4_out=0.0
     ALLOCATE(JNO3_out(MTLOC));           JNO3_out=0.0
#endif	 
     !WLong moved here from wqm_modules.F
     ALLOCATE(NH41TM1S(MTLOC));           NH41TM1S  = 0.0
     ALLOCATE(NH42TM1S(MTLOC));           NH42TM1S  = 0.0
     ALLOCATE(NH4T2TM1S(MTLOC));          NH4T2TM1S = 0.0
     ALLOCATE(NO31TM1S(MTLOC));           NO31TM1S  = 0.0
     ALLOCATE(NO32TM1S(MTLOC));           NO32TM1S  = 0.0
     ALLOCATE(NO3T2TM1S(MTLOC));          NO3T2TM1S = 0.0
     ALLOCATE(HS1TM1S(MTLOC));            HS1TM1S   = 0.0
     ALLOCATE(HS2TM1S(MTLOC));            HS2TM1S   = 0.0
     ALLOCATE(HST2TM1S(MTLOC));           HST2TM1S  = 0.0
     ALLOCATE(SI1TM1S(MTLOC));            SI1TM1S   = 0.0
     ALLOCATE(SI2TM1S(MTLOC));            SI2TM1S   = 0.0
     ALLOCATE(SIT2TM1S(MTLOC));           SIT2TM1S  = 0.0
     ALLOCATE(PO41TM1S(MTLOC));           PO41TM1S  = 0.0
     ALLOCATE(PO42TM1S(MTLOC));           PO42TM1S  = 0.0
     ALLOCATE(PO4T2TM1S(MTLOC));          PO4T2TM1S = 0.0
    
     !WLong moved here from wqm_modules.F
     ALLOCATE(CH41TM1S(MTLOC));           CH41TM1S     = 0.0
     ALLOCATE(CH42TM1S(MTLOC));           CH42TM1S     = 0.0
     ALLOCATE(CH4T2TM1S(MTLOC));          CH4T2TM1S = 0.0
     ALLOCATE(SO4T2TM1S(MTLOC));          SO4T2TM1S = 0.0
     
     !WLong moved here from wqm_modules.F
     ALLOCATE(BURIALN(MTLOC));           BURIALN = 0.0
     ALLOCATE(BURIALP(MTLOC));           BURIALP = 0.0
     ALLOCATE(BURIALC(MTLOC));           BURIALC = 0.0

     !WLong moved here from wqm_modules.F
     ALLOCATE(DIAGENC(MTLOC));            DIAGENC = 0.0
    
    !WLong moved this here from wqm_modules.F
    ALLOCATE(AG3CFL(MTLOC));             AG3CFL = 0.0
    ALLOCATE(AG3NFL(MTLOC));             AG3NFL = 0.0
    ALLOCATE(AG3PFL(MTLOC));             AG3PFL = 0.0
    
    !WLong moved this here from wqm_modules.F
    ALLOCATE(ASDTMP(MTLOC));             ASDTMP = 0.0
    

    
        ALLOCATE(HSED(MTLOC));                HSED   = 0.0  !Moved from mod_wqm.F
        ALLOCATE(HSED1(MTLOC));               HSED1   = 0.0  !Moved from mod_wqm.F
        
        !used for writing and reading restart file
        ALLOCATE(CTEMP_GL(MGL));             CTEMP_GL     = 0.0
        ALLOCATE(CPOP_GL(MGL,3));            CPOP_GL      = 0.0
        ALLOCATE(CPON_GL(MGL,3));            CPON_GL      = 0.0
        ALLOCATE(CPOC_GL(MGL,3));            CPOC_GL      = 0.0
        ALLOCATE(CPOS_GL(MGL));              CPOS_GL      = 0.0
        ALLOCATE(PO4T2TM1S_GL(MGL));         PO4T2TM1S_GL = 0.0
        ALLOCATE(NH4T2TM1S_GL(MGL));         NH4T2TM1S_GL = 0.0
        ALLOCATE(NO3T2TM1S_GL(MGL));         NO3T2TM1S_GL = 0.0
        ALLOCATE(HST2TM1S_GL(MGL));          HST2TM1S_GL  = 0.0
        ALLOCATE(CH4T2TM1S_GL(MGL));         CH4T2TM1S_GL = 0.0
        ALLOCATE(CH41TM1S_GL(MGL));          CH41TM1S_GL  = 0.0
        ALLOCATE(SO4T2TM1S_GL(MGL));         SO4T2TM1S_GL = 0.0
        ALLOCATE(SIT2TM1S_GL(MGL));          SIT2TM1S_GL  = 0.0
        
        ALLOCATE(BENSTRTM1S_GL(MGL));          BENSTRTM1S_GL  = 0.0
        ALLOCATE(SODTM1S_GL(MGL));          SODTM1S_GL(:) = 0.d0
    
        !used for outputing sediment flux variables into history and station files
        
        ALLOCATE(JPOC_GL(MGL,3));            JPOC_GL        =0.d0
        ALLOCATE(JPON_GL(MGL,3));            JPON_GL        =0.d0
        ALLOCATE(JPOP_GL(MGL,3));            JPOP_GL        =0.d0
        ALLOCATE(JPOS_GL(MGL));                JPOS_GL        =0.d0
        !ALLOCATE(O20_GL(MGL));                O20_GL        =0.d0    !given by C2_GL(:,:,27)
        !ALLOCATE(D_GL(MGL));                D_GL        =0.d0    !already defined
        !ALLOCATE(T_GL(MGL,KBM1));            T_GL        =0.d0    !already defined
        !ALLOCATE(NH40_GL(MGL));            NH40_GL        =0.d0    !given by C2_GL
        !ALLOCATE(NO30_GL(MGL));            NO30_GL        =0.d0    !given by C2_GL
        !ALLOCATE(SI0_GL(MGL));                SI0_GL        =0.d0    !given by C2_GL
         ALLOCATE(CH40_GL(MGL));            CH40_GL        =0.d0
        !ALLOCATE(SAL_GL(MGL));                SAL_GL        =0.d0    !given by C2_GL
        !ALLOCATE(SOD_GL(MGL));                SOD_GL        =0.d0    !given by SODTM1S_GL
        ALLOCATE(JNH4_GL(MGL));                JNH4_GL        =0.d0    !
        ALLOCATE(JNO3_GL(MGL));                JNO3_GL        =0.d0
        ALLOCATE(BENDEN_GL(MGL));            BENDEN_GL    =0.d0
        ALLOCATE(JCH4_GL(MGL));                JCH4_GL        =0.d0
        ALLOCATE(JCH4G_GL(MGL));            JCH4G_GL    =0.d0
        ALLOCATE(JHS_GL(MGL));                JHS_GL        =0.d0
        ALLOCATE(JPO4_GL(MGL));                JPO4_GL        =0.d0
        ALLOCATE(JSI_GL(MGL));                JSI_GL        =0.d0
        ALLOCATE(NH41_GL(MGL));                NH41_GL        =0.d0
        ALLOCATE(NH42_GL(MGL));                NH42_GL        =0.d0

        ALLOCATE(NO31_GL(MGL));                NO31_GL        =0.d0
        ALLOCATE(NO32_GL(MGL));                NO32_GL        =0.d0
        ALLOCATE(PO41_GL(MGL));                PO41_GL        =0.d0
        ALLOCATE(PO42_GL(MGL));                PO42_GL        =0.d0        
        ALLOCATE(SI1_GL(MGL));                SI1_GL        =0.d0
        ALLOCATE(SI2_GL(MGL));                SI2_GL        =0.d0
        !ALLOCATE(CH41_GL(MGL));            CH41_GL        =0.d0    !given by CH41TM1S_GL
        ALLOCATE(CH42_GL(MGL));                CH42_GL        =0.d0
        ALLOCATE(HS1_GL(MGL));                HS1_GL        =0.d0
        ALLOCATE(HS2_GL(MGL));                HS2_GL        =0.d0                                

        ALLOCATE(HSED1_GL(MGL));            HSED1_GL    =0.d0    !
        !ALLOCATE(BENSTR_GL(MGL));            BENSTR_GL    =0.d0    !given by BENSTRTM1S_GL
        
      
        ALLOCATE(MTVEL(MTLOC));              MTVEL  = 0.0    !WL moved from mod_wqm.F        
        ALLOCATE(CPOC(MTLOC,3));             CPOC   = 0.0
        ALLOCATE(CPOP(MTLOC,3));             CPOP   = 0.0
        ALLOCATE(CPON(MTLOC,3));             CPON   = 0.0
        ALLOCATE(CPOS(MTLOC));               CPOS   = 0.0      
        ALLOCATE(CPO4(MTLOC));               CPO4   = 0.0      
        ALLOCATE(CNO3(MTLOC));               CNO3   = 0.0
        ALLOCATE(CNH4(MTLOC));               CNH4   = 0.0     


        ALLOCATE(CCH4(MTLOC));               CCH4   = 0.0      
        ALLOCATE(CSO4(MTLOC));               CSO4   = 0.0     
        ALLOCATE(CHS(MTLOC));                CHS    = 0.0
        ALLOCATE(CSI(MTLOC));                CSI    = 0.0        
        
        ALLOCATE(DIAGP(MTLOC));             DIAGP  = 0.0
        ALLOCATE(DIAGN(MTLOC));             DIAGN  = 0.0
        ALLOCATE(DIAGC(MTLOC));             DIAGC  = 0.0
        ALLOCATE(DIAGS(MTLOC));             DIAGS  = 0.0
   
        ALLOCATE(JPOC(MTLOC,3));           JPOC = 0.0
        ALLOCATE(JPOP(MTLOC,3));           JPOP = 0.0 
        ALLOCATE(JPON(MTLOC,3));           JPON = 0.0
        ALLOCATE(JPOS(MTLOC));             JPOS = 0.0 
#if defined(BUDGET)			
		ALLOCATE(JPOC_out(MTLOC,3)); JPOC_out=0.		
		ALLOCATE(JPOP_out(MTLOC,3)); JPOP_out=0.
		ALLOCATE(JPON_out(MTLOC,3)); JPON_out=0.
#endif        
        ALLOCATE(JPOCaccum(MTLOC,3));           JPOCaccum = 0.0
        ALLOCATE(JPOPaccum(MTLOC,3));           JPOPaccum = 0.0 
        ALLOCATE(JPONaccum(MTLOC,3));           JPONaccum = 0.0
        ALLOCATE(JPOSaccum(MTLOC));             JPOSaccum = 0.0 
        
END SUBROUTINE SED_ALLOC

!********************************************************************************
!**                    S U B R O U T I N E   SED_DEALLOC                       **
!********************************************************************************

  SUBROUTINE SED_DEALLOC
  
     !WLong moved here from wqm_main.F

     IF(ALLOCATED(FRPOP)) DEALLOCATE (FRPOP)
     IF(ALLOCATED(FRPON)) DEALLOCATE (FRPON)
     IF(ALLOCATED(FRPOC)) DEALLOCATE (FRPOC)

   !   !WLong moved here from wqm_main.F
   !   IF(ALLOCATED(WSSBNET)) DEALLOCATE (WSSBNET)
   !   IF(ALLOCATED(WSLBNET)) DEALLOCATE (WSLBNET)
   !   IF(ALLOCATED(WSRBNET)) DEALLOCATE (WSRBNET)
   !   IF(ALLOCATED(WS1BNET)) DEALLOCATE (WS1BNET)
   !   IF(ALLOCATED(WS2BNET)) DEALLOCATE (WS2BNET)
   !   IF(ALLOCATED(WS3BNET)) DEALLOCATE (WS3BNET)
   !   IF(ALLOCATED(WSUBNET)) DEALLOCATE (WSUBNET)
    
     !WLong moved here from wqm_main.F
     IF(ALLOCATED(VSED))  DEALLOCATE (VSED)
     IF(ALLOCATED(VPMIX)) DEALLOCATE (VPMIX)
     IF(ALLOCATED(VDMIX)) DEALLOCATE (VDMIX)
    
     !WLong moved here from wqm_main.F
     IF(ALLOCATED(PON1TM1S)) DEALLOCATE (PON1TM1S)
     IF(ALLOCATED(PON2TM1S)) DEALLOCATE (PON2TM1S)
     IF(ALLOCATED(PON3TM1S)) DEALLOCATE (PON3TM1S)
     IF(ALLOCATED(POC1TM1S)) DEALLOCATE (POC1TM1S)
     IF(ALLOCATED(POC2TM1S)) DEALLOCATE (POC2TM1S)
     IF(ALLOCATED(POC3TM1S)) DEALLOCATE (POC3TM1S)
     IF(ALLOCATED(POP1TM1S)) DEALLOCATE (POP1TM1S)
     IF(ALLOCATED(POP2TM1S)) DEALLOCATE (POP2TM1S)
     IF(ALLOCATED(POP3TM1S)) DEALLOCATE (POP3TM1S)
     IF(ALLOCATED(PSISEDTM1S)) DEALLOCATE (PSISEDTM1S)
  
     !WLong added POC1TM1S_SED_DF POC2TM1S_SED_DF
     
!     IF(ALLOCATED(POC1TM1S_SED_DF))DEALLOCATE(POC1TM1S_SED_DF)
!     IF(ALLOCATED(POC2TM1S_SED_DF))DEALLOCATE(POC2TM1S_SED_DF)
     
     
     !WLong moved these from wqm_main.F
     IF(ALLOCATED(BENSTRTM1S)) DEALLOCATE (BENSTRTM1S)
     IF(ALLOCATED(BFORMAXS)) DEALLOCATE (BFORMAXS)
     IF(ALLOCATED(ISWBENS)) DEALLOCATE (ISWBENS)

     IF(ALLOCATED(SODTM1S)) DEALLOCATE (SODTM1S)
     IF(ALLOCATED(JNH4TM1S)) DEALLOCATE (JNH4TM1S)
     IF(ALLOCATED(JNO3TM1S)) DEALLOCATE (JNO3TM1S)
     IF(ALLOCATED(JPO4TM1S)) DEALLOCATE (JPO4TM1S)
     IF(ALLOCATED(JSITM1S)) DEALLOCATE (JSITM1S)
     IF(ALLOCATED(JCH4TM1S)) DEALLOCATE (JCH4TM1S)
     IF(ALLOCATED(JCH4GTM1S)) DEALLOCATE (JCH4GTM1S)
     IF(ALLOCATED(JHSTM1S)) DEALLOCATE (JHSTM1S)

#if defined(BUDGET)     
	  IF(ALLOCATED(JNH4_out)) DEALLOCATE (JNH4_out)
     IF(ALLOCATED(JNO3_out)) DEALLOCATE (JNO3_out)
#endif     
     !WLong moved these from wqm_main.F
     IF(ALLOCATED(NH41TM1S)) DEALLOCATE (NH41TM1S)
     IF(ALLOCATED(NH42TM1S)) DEALLOCATE (NH42TM1S)
     IF(ALLOCATED(NH4T2TM1S)) DEALLOCATE (NH4T2TM1S)
     IF(ALLOCATED(NO31TM1S)) DEALLOCATE (NO31TM1S)
     IF(ALLOCATED(NO32TM1S)) DEALLOCATE (NO32TM1S)
     IF(ALLOCATED(NO3T2TM1S)) DEALLOCATE (NO3T2TM1S)
     IF(ALLOCATED(HS1TM1S)) DEALLOCATE (HS1TM1S)
     IF(ALLOCATED(HS2TM1S)) DEALLOCATE (HS2TM1S)
     IF(ALLOCATED(HST2TM1S)) DEALLOCATE (HST2TM1S)
     IF(ALLOCATED(SI1TM1S)) DEALLOCATE (SI1TM1S)
     IF(ALLOCATED(SI2TM1S)) DEALLOCATE (SI2TM1S)
     IF(ALLOCATED(SIT2TM1S)) DEALLOCATE (SIT2TM1S)
     IF(ALLOCATED(PO41TM1S)) DEALLOCATE (PO41TM1S)
     IF(ALLOCATED(PO42TM1S)) DEALLOCATE (PO42TM1S)
     IF(ALLOCATED(PO4T2TM1S)) DEALLOCATE (PO4T2TM1S)

     !WLong moved these here  wqm_main.F
     IF(ALLOCATED(CH41TM1S)) DEALLOCATE (CH41TM1S)
     IF(ALLOCATED(CH42TM1S)) DEALLOCATE (CH42TM1S)
     IF(ALLOCATED(CH4T2TM1S)) DEALLOCATE (CH4T2TM1S)
     IF(ALLOCATED(SO4T2TM1S)) DEALLOCATE (SO4T2TM1S)
     !WLong moved these here  wqm_main.F
    
     IF(ALLOCATED(BURIALN)) DEALLOCATE (BURIALN)
     IF(ALLOCATED(BURIALP)) DEALLOCATE (BURIALP)
     IF(ALLOCATED(BURIALC)) DEALLOCATE (BURIALC)

     !WLong moved here from wqm_main.F
     IF(ALLOCATED(DIAGENC)) DEALLOCATE (DIAGENC)
    
     !WLong moved here from wqm_main.F
     IF(ALLOCATED(AG3CFL)) DEALLOCATE (AG3CFL)
     IF(ALLOCATED(AG3NFL)) DEALLOCATE (AG3NFL)
     IF(ALLOCATED(AG3PFL)) DEALLOCATE (AG3PFL)
     IF(ALLOCATED(ASDTMP)) DEALLOCATE (ASDTMP)
     
     IF(ALLOCATED(HSED))DEALLOCATE(HSED)
     IF(ALLOCATED(HSED))DEALLOCATE(HSED1)
     
     !global arrays used as input initial conditions and output
     IF(ALLOCATED(CTEMP_GL))      DEALLOCATE(CTEMP_GL)
     IF(ALLOCATED(CPOP_GL))       DEALLOCATE(CPOP_GL)
     IF(ALLOCATED(CPON_GL))       DEALLOCATE(CPON_GL)
     IF(ALLOCATED(CPOC_GL))       DEALLOCATE(CPOC_GL)
     IF(ALLOCATED(CPOS_GL))       DEALLOCATE(CPOS_GL)
     IF(ALLOCATED(PO4T2TM1S_GL))  DEALLOCATE(PO4T2TM1S_GL)
     IF(ALLOCATED(NH4T2TM1S_GL))  DEALLOCATE(NH4T2TM1S_GL)
     IF(ALLOCATED(NO3T2TM1S_GL))  DEALLOCATE(NO3T2TM1S_GL)
     IF(ALLOCATED(HST2TM1S_GL))   DEALLOCATE(HST2TM1S_GL)
     IF(ALLOCATED(CH4T2TM1S_GL))  DEALLOCATE(CH4T2TM1S_GL)
     IF(ALLOCATED(CH41TM1S_GL))   DEALLOCATE(CH41TM1S_GL)
     IF(ALLOCATED(SO4T2TM1S_GL))  DEALLOCATE(SO4T2TM1S_GL)
     IF(ALLOCATED(SIT2TM1S_GL))   DEALLOCATE(SIT2TM1S_GL)
     IF(ALLOCATED(BENSTRTM1S_GL)) DEALLOCATE(BENSTRTM1S_GL)
     IF(ALLOCATED(SODTM1S_GL))    DEALLOCATE(SODTM1S_GL)
     
     !global arrays used for output to history and station outputs
     IF(ALLOCATED(JPOC_GL))        DEALLOCATE(JPOC_GL)
     IF(ALLOCATED(JPON_GL))        DEALLOCATE(JPON_GL)
     IF(ALLOCATED(JPOP_GL))        DEALLOCATE(JPOP_GL)

#if defined(BUDGET)	
! Added by B Clark for budget calculation
	 IF(ALLOCATED(JPOC_out))        DEALLOCATE(JPOC_out)
     IF(ALLOCATED(JPON_out))        DEALLOCATE(JPON_out)
     IF(ALLOCATED(JPOP_out))        DEALLOCATE(JPOP_out)
#endif	 
	 
     IF(ALLOCATED(JPOS_GL))        DEALLOCATE(JPOS_GL)
     !IF(ALLOCATED(O20_GL))        DEALLOCATE(O20_GL)        !given by C2_GL(:,:,27)
     !IF(ALLOCATED(D_GL))        DEALLOCATE(D_GL)        !already defined
     !IF(ALLOCATED(T_GL))        DEALLOCATE(T_GL)        !already defined
     !IF(ALLOCATED(NH40_GL))    DEALLOCATE(NH40_GL)        !given by C2_GL
     !IF(ALLOCATED(NO30_GL))    DEALLOCATE(NO30_GL)        !given by C2_GL
     !IF(ALLOCATED(SI0_GL))        DEALLOCATE(SI0_GL)         !given by C2_GL
      IF(ALLOCATED(CH40_GL))    DEALLOCATE(CH40_GL)    
     !IF(ALLOCATED(SAL_GL))        DEALLOCATE(SAL_GL)        !given by C2_GL
     !IF(ALLOCATED(SOD_GL))        DEALLOCATE(SOD_GL)        !given by SODTM1S_GL
     IF(ALLOCATED(JNH4_GL))        DEALLOCATE(JNH4_GL)        !
     IF(ALLOCATED(JNO3_GL))        DEALLOCATE(JNO3_GL)
     IF(ALLOCATED(BENDEN_GL))   DEALLOCATE(BENDEN_GL)
     IF(ALLOCATED(JCH4_GL))        DEALLOCATE(JCH4_GL)
     IF(ALLOCATED(JCH4G_GL))    DEALLOCATE(JCH4G_GL)
     IF(ALLOCATED(JHS_GL))        DEALLOCATE(JHS_GL)
     IF(ALLOCATED(JPO4_GL))        DEALLOCATE(JPO4_GL)
     IF(ALLOCATED(JSI_GL))        DEALLOCATE(JSI_GL)
     IF(ALLOCATED(NH41_GL))       DEALLOCATE(NH41_GL)
     IF(ALLOCATED(NH42_GL))        DEALLOCATE(NH42_GL)

     IF(ALLOCATED(NO31_GL))        DEALLOCATE(NO31_GL)
     IF(ALLOCATED(NO32_GL))        DEALLOCATE(NO32_GL)
     IF(ALLOCATED(PO41_GL))        DEALLOCATE(PO41_GL)
     IF(ALLOCATED(PO42_GL))        DEALLOCATE(PO42_GL)
     IF(ALLOCATED(SI1_GL))        DEALLOCATE(SI1_GL)
     IF(ALLOCATED(SI2_GL))        DEALLOCATE(SI2_GL)
     !IF(ALLOCATED(CH41_GL))    DEALLOCATE(CH41_GL)    !given by CH41TM1S_GL
     IF(ALLOCATED(CH42_GL))        DEALLOCATE(CH42_GL)
     IF(ALLOCATED(HS1_GL))        DEALLOCATE(HS1_GL)
     IF(ALLOCATED(HS2_GL))        DEALLOCATE(HS2_GL)

     !IF(ALLOCATED(POC1_GL))    DEALLOCATE(POC1_GL)
     !IF(ALLOCATED(POC2_GL))    DEALLOCATE(POC2_GL)
     !IF(ALLOCATED(POC3_GL))    DEALLOCATE(POC3_GL)

     !IF(ALLOCATED(PON1_GL))    DEALLOCATE(PON1_GL)    !given by CPON_GL
     !IF(ALLOCATED(PON2_GL))    DEALLOCATE(PON2_GL)
     !IF(ALLOCATED(PON3_GL))    DEALLOCATE(PON3_GL)
                                
     !IF(ALLOCATED(POP1_GL))    DEALLOCATE(POP1_GL)    !given by CPOP_GL
     !IF(ALLOCATED(POP2_GL))    DEALLOCATE(POP2_GL)
     !IF(ALLOCATED(POP3_GL))    DEALLOCATE(POP3_GL)
                                
     !IF(ALLOCATED(PSISED_GL))  DEALLOCATE(PSISED_GL)    !given by CPOS_GL
     IF(ALLOCATED(HSED1_GL))    DEALLOCATE(HSED1_GL)
     !IF(ALLOCATED(BENSTR_GL))    DEALLOCATE(BENSTR_GL)    !given by BENSTRTM1S_GL
        
     !used as local calculations 
     IF(ALLOCATED(MTVEL))DEALLOCATE(MTVEL)        !WLong moved to mod_sed.F
     
     IF(ALLOCATED(CPOC))DEALLOCATE(CPOC)
     IF(ALLOCATED(CPON))DEALLOCATE(CPON)             
     IF(ALLOCATED(CPOP))DEALLOCATE(CPOP)
     IF(ALLOCATED(CPOS))DEALLOCATE(CPOS)
     
     IF(ALLOCATED(CPO4))DEALLOCATE(CPO4)
     IF(ALLOCATED(CNO3))DEALLOCATE(CNO3)
     IF(ALLOCATED(CNH4))DEALLOCATE(CNH4)

     IF(ALLOCATED(CCH4))DEALLOCATE(CCH4)
     IF(ALLOCATED(CSO4))DEALLOCATE(CSO4)
     IF(ALLOCATED(CHS))DEALLOCATE(CHS)
     IF(ALLOCATED(CSI))DEALLOCATE(CSI)        
     
     IF(ALLOCATED(DIAGC))DEALLOCATE(DIAGC)           
     IF(ALLOCATED(DIAGN))DEALLOCATE(DIAGN)           
     IF(ALLOCATED(DIAGP))DEALLOCATE(DIAGP)
     IF(ALLOCATED(DIAGS))DEALLOCATE(DIAGS)
     
     IF(ALLOCATED(JPOC))DEALLOCATE(JPOC)           
     IF(ALLOCATED(JPON))DEALLOCATE(JPON)
     IF(ALLOCATED(JPOP))DEALLOCATE(JPOP)
     IF(ALLOCATED(JPOS))DEALLOCATE(JPOS)
                
     IF(ALLOCATED(JPOCaccum))DEALLOCATE(JPOCaccum)           
     IF(ALLOCATED(JPONaccum))DEALLOCATE(JPONaccum)
     IF(ALLOCATED(JPOPaccum))DEALLOCATE(JPOPaccum)
     IF(ALLOCATED(JPOSaccum))DEALLOCATE(JPOSaccum)
     
  END SUBROUTINE SED_DEALLOC

!********************************************************************************
!**                    S U B R O U T I N E   S E D _ R E A D                   **
!********************************************************************************
   SUBROUTINE SED_READ
   
   USE MOD_SED_DF_EXCHANGE_VARS, ONLY:    &
        POC1TM1S_SED_DF    ,              &
        POC2TM1S_SED_DF    ,              &
        M1_SED_DF,                        &
        M2_SED_DF

#  if defined (MULTIPROCESSOR)
        USE MOD_PAR, ONLY:   &      
            !EL_PID,    &        !(:)           !!PROCESSOR OWNER OF GLOBAL ELEMENT
            !ELID,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT
            !NLID,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE 
            !ELID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
            !NLID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
            !OBN_GL2LOC,    &    !(:)           !!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
            !EGID,        &    !(:)           !!GLOBAL ID OF LOCAL ELEMENT
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
            NHN,          &    !            !!NUMBER OF HALO NODES
            HN_LST!,     &    !            !!LIST OF HALO NODES 
            !HN_OWN,     &    !            !!PRIMARY OWNER OF HALO NODES
            !EC,            &    !
            !NC,            &    !
            !BNC,         &    !
            !EMAP,        &    !
            !NMAP,         &    !
            !PNE,        &    !(:)        !!NUMBER OF ELEMENTS IN EACH PROC
            !PNN,        &    !(:)        !!NUMBER OF NODES IN EACH PROC
            !PNHE,        &    !(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
            !PNBN,        &    !(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
            !PMBM,        &    !(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
            !PNHN,        &    !(:)       !!NUMBER OF HALO NODES IN EACH PROC
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
            !SORT,    &
            !PREAD,    &
            !PWRITE,&
            !PPRINT, &
            !LOCK, &
            !UNLOCK, &
            !GETLOC
#  endif  

   IMPLICIT NONE
   SAVE    

!***** Variable declarations

   CHARACTER(LEN=24) :: FRNAME(14)

   INTEGER :: I,J, JG, JT
   REAL(SP), ALLOCATABLE :: RTMP11(:),& 
                        RTMP12(:),    &
                        RTMP13(:),    &
                        RTMP14(:),    &
                        RTMP15(:),    &
                        RTMP16(:),    &
                        RTMP17(:) 

!***** Data declarations

   DATA FRNAME                                  &
        /'Group 1 algal phosphorus',            &     !1    P
         'Group 2 algal phosphorus',            &     !2
         'Group 3 algal phosphorus',            &     !3
         'Detrital org phosphorus ',            &     !4
         'Group 1 algal nitrogen  ',            &     !5    N
         'Group 2 algal nitrogen  ',            &     !6
         'Group 3 algal nitrogen  ',            &     !7
         'Detrital org nitrogen   ',            &     !8
         'Group 1 algal carbon    ',            &     !9    C
         'Group 2 algal carbon    ',            &     !10
         'Group 3 algal carbon    ',            &     !11
         'Benthic algal carbon    ',            &     !12   Benthic algae C
         'Benthic algal nitrogen  ',            &     !13   Benthic algae N
         'Benthic algal phosphorus'/                  !14   Benthic algae P

    IF (BENTHIC_OUTPUT) THEN
        IF(MSR)THEN
            WRITE(*,*)'BFOFN=',BFOFN
            OPEN (BFO,FILE=BFOFN)
        ENDIF
    ENDIF
         
!********************************************************************************
!**                                  Inputs                                    **
!********************************************************************************


999  FORMAT(:///8X,F8.0,I8,I8,I8,I8,I8)      !LB added ,I8 for NDTSED
9999 FORMAT(/)   
1020 FORMAT(://(8X,7F8.1))
1022 FORMAT(://(8X,3F8.1))
1040 FORMAT(://(8X,6F8.1))
1060 FORMAT(//8X,2A8,F8.0)
1070 FORMAT(://8x,I8,I8,F8.0)
1080 FORMAT(://8X,9F8.0)
1090 FORMAT(//(8X,A72))
1099 FORMAT(://(8X,I8,3F8.1))

   READ(BFI,9999,ERR=10100)                      !move to third line with format '(/)', 
                                                !READ empty moves to second line, / moves another line   
   
   READ(BFI,999,ERR=10100)  HSEDALL, SSTATEG3,SSTATEIC ,SSTATEAVG, NDTSED, QUASISS  ! depth of sediments, steady state flag for G3, number of WCtimesteps per sedim timestep (DTsed=NDTSED*DTwc) -added by LB
   
   IF(MSR)WRITE(*,*),'HSEDALL=', HSEDALL,         &!
                     'SSTATEG3=',SSTATEG3,        &!
                     'SSTATEIC=',SSTATEIC,        &!
                     'SSTATEAVG=',SSTATEAVG,      &!
                     'QUASISS=', QUASISS

   
   READ(BFI,1070,ERR=10100)  ITVWCLF,NTVWCLF,WCLFDAYS    !flag for focing with time varying water column information, length in days
                                                        !in the file (all files sholuld have same length)
   
   IF(MSR)WRITE(*,*)'ITVWCLF=',ITVWCLF        !==1 to have overyling water forcing input files
   IF(MSR)WRITE(*,*)'NTVWCLF=',NTVWCLF        !number of files
   IF(MSR)WRITE(*,*)'WCLFDAYS=',WCLFDAYS    !lenght for each file (days)
   
   IF(ITVWCLF==1 .AND. WCLFDAYS<=0.0 )THEN
    WRITE(*,*)'Error, WCLFDAYS should be greater than zero'
    GOTO 10100
   ENDIF    
   IF(ITVWCLF==1 .AND. NTVWCLF <1)THEN    
    WRITE(*,*)'Error, NTVWCLF should be >=1'
    GOTO 10100
   ENDIF
   
   READ(BFI,1080,ERR=10100)  DIFFT              !diffusion rate at sediment water interface
   DIFFT = 0.0001*DIFFT                          !convert from cm^2/sec to m^2/sec
   
   IF(MSR)WRITE(*,*)'DIFFT=', DIFFT
   
   READ(BFI,1080,ERR=10100)  SALTSW, SALTND     !freshwater salt water formulation for SOD  and nitrification/denitrification
   IF(MSR)WRITE(*,*)'SALTSW=', SALTSW,'SALTND=', SALTND

   READ(BFI,1080,ERR=10100)  FRPALG1            !Fractions of P in water column algae group 1 for G1,G2,G3
   READ(BFI,1080,ERR=10100)  FRPALG2             !Fractions of P in water column algae group 2 for G1,G2,G3
   READ(BFI,1080,ERR=10100)  FRPALG3             !Fractions of P in water column algae group 3 for G1,G2,G3
   IF(MSR)WRITE(*,*)'FRPALG1=', FRPALG1
   IF(MSR)WRITE(*,*)'FRPALG2=', FRPALG2   
   IF(MSR)WRITE(*,*)'FRPALG3=', FRPALG3   

   READ(BFI,1080,ERR=10100)  FRNALG1            !Fractions of N in water column algae group 1 for G1,G2,G3
   READ(BFI,1080,ERR=10100)  FRNALG2              !Fractions of N in water column algae group 2 for G1,G2,G3
   READ(BFI,1080,ERR=10100)  FRNALG3            !Fractions of N in water column algae group 3 for G1,G2,G3
   IF(MSR)WRITE(*,*)'FRNALG1=', FRNALG1
   IF(MSR)WRITE(*,*)'FRNALG2=', FRNALG2   
   IF(MSR)WRITE(*,*)'FRNALG3=', FRNALG3   

   READ(BFI,1080,ERR=10100)  FRCALG1            !Fractions of C in water column algae group 1 for G1, G2, G3
   READ(BFI,1080,ERR=10100)  FRCALG2            !Fractions of C in water column algae group 2 for G1, G2, G3
   READ(BFI,1080,ERR=10100)  FRCALG3            !Fractions of C in water column algae group 3 for G1, G2, G3
   IF(MSR)WRITE(*,*)'FRCALG1=', FRCALG1
   IF(MSR)WRITE(*,*)'FRCALG2=', FRCALG2   
   IF(MSR)WRITE(*,*)'FRCALG3=', FRCALG3   
   
   READ(BFI,1080,ERR=10100) (KPDIAG(JG),DPTHTA(JG),JG=1,3)
    
   KPOP1=KPDIAG(1)
   KPOP2=KPDIAG(2)
   KPOP3=KPDIAG(3)

   THTAPOP1=DPTHTA(1)
   THTAPOP2=DPTHTA(2)
   THTAPOP3=DPTHTA(3)
    
   IF(MSR)WRITE(*,*)'KPOP1=', KPOP1
   IF(MSR)WRITE(*,*)'KPOP2=', KPOP2   
   IF(MSR)WRITE(*,*)'KPOP3=', KPOP3       
   
   IF(MSR)WRITE(*,*)'THTAPOP1=', THTAPOP1
   IF(MSR)WRITE(*,*)'THTAPOP2=', THTAPOP2   
   IF(MSR)WRITE(*,*)'THTAPOP3=', THTAPOP3   
   
   READ(BFI,1080,ERR=10100) (KNDIAG(JG),DNTHTA(JG),JG=1,3)

    !WLong added the following and removed the equivalence declaration in wqm_modules.F
    KPON1=KNDIAG(1)
    KPON2=KNDIAG(2)
    KPON3=KNDIAG(3)
      THTAPON1=DNTHTA(1)
    THTAPON2=DNTHTA(2)
    THTAPON3=DNTHTA(3)
    
   IF(MSR)WRITE(*,*)'KPON1=', KPON1
   IF(MSR)WRITE(*,*)'KPON2=', KPON2   
   IF(MSR)WRITE(*,*)'KPON3=', KPON3       
   
   IF(MSR)WRITE(*,*)'THTAPON1=', THTAPON1
   IF(MSR)WRITE(*,*)'THTAPON2=', THTAPON2   
   IF(MSR)WRITE(*,*)'THTAPON3=', THTAPON3       

   READ(BFI,1080,ERR=10100) (KCDIAG(JG),DCTHTA(JG),JG=1,3)

    !WLong added the following and removed the equivalence declaration in wqm_modules.F
    KPOC1=KCDIAG(1)
    KPOC2=KCDIAG(2)
    KPOC3=KCDIAG(3)
    THTAPOC1=DCTHTA(1)
    THTAPOC2=DCTHTA(2)
    THTAPOC3=DCTHTA(3)
    
   IF(MSR)WRITE(*,*)'KPOC1=', KPOC1
   IF(MSR)WRITE(*,*)'KPOC2=', KPOC2   
   IF(MSR)WRITE(*,*)'KPOC3=', KPOC3       
   
   IF(MSR)WRITE(*,*)'THTAPOC1=', THTAPOC1
   IF(MSR)WRITE(*,*)'THTAPOC2=', THTAPOC2   
   IF(MSR)WRITE(*,*)'THTAPOC3=', THTAPOC3       
   
   READ(BFI,1080,ERR=10100)  KSI,THTASI,THTASISAT
   IF(MSR)WRITE(*,*)'KSI=', KSI
   IF(MSR)WRITE(*,*)'THTASI=', THTASI   
   IF(MSR)WRITE(*,*)'THTASISAT=', THTASISAT   

   
   READ(BFI,1080,ERR=10100)  M1,M2,THTADP,THTADD
   IF(MSR)WRITE(*,*)'M1=', M1
   IF(MSR)WRITE(*,*)'M2=', M2   
   IF(MSR)WRITE(*,*)'THTADP=', THTADP
   IF(MSR)WRITE(*,*)'THTADD=', THTADD
   
   
   READ(BFI,1080,ERR=10100)  KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2,THTAKMNH4
   
   IF(MSR)WRITE(*,*)'KAPPNH4F=', KAPPNH4F
   IF(MSR)WRITE(*,*)'KAPPNH4S=', KAPPNH4S   
   IF(MSR)WRITE(*,*)'PIENH4=', PIENH4
   IF(MSR)WRITE(*,*)'THTANH4=', THTANH4
   IF(MSR)WRITE(*,*)'KMNH4=', KMNH4   
   IF(MSR)WRITE(*,*)'KMNH4O2=', KMNH4O2
   IF(MSR)WRITE(*,*)'THTAKMNH4=',THTAKMNH4
   
   READ(BFI,1080,ERR=10100)  KAPPNO3F,KAPPNO3S,K2NO3,THTANO3
   
   IF(MSR)WRITE(*,*)'KAPPNO3F=', KAPPNO3F
   IF(MSR)WRITE(*,*)'KAPPNO3S=', KAPPNO3S   
   IF(MSR)WRITE(*,*)'K2NO3=', K2NO3
   IF(MSR)WRITE(*,*)'THTANO3=', THTANO3
   
   READ(BFI,1080,ERR=10100)  KAPP1HSD,KAPP1HSP,PIE1HS,PIE2HS,THTAH2S,KMHSO2
   
   IF(MSR)WRITE(*,*)'KAPP1HSD=', KAPP1HSD
   IF(MSR)WRITE(*,*)'KAPP1HSP=', KAPP1HSP   
   IF(MSR)WRITE(*,*)'PIE1HS=', PIE1HS
   IF(MSR)WRITE(*,*)'PIE2HS=', PIE2HS
   IF(MSR)WRITE(*,*)'THTAH2S=', THTAH2S   
   IF(MSR)WRITE(*,*)'KMHSO2=', KMHSO2
   
   READ(BFI,1080,ERR=10100)  CSISATT20,DPIE1SI,PIE2SI,KMPSI  !KMPSI is gSi/m^3 on input!!!!
   KMPSI=KMPSI*1000.d0  !convert from gSi/m^3 to mgSi/m^3
   IF(MSR)WRITE(*,*)'CSISATT20=', CSISATT20
   IF(MSR)WRITE(*,*)'DPIE1SI=', DPIE1SI   
   IF(MSR)WRITE(*,*)'PIE2SI=', PIE2SI
   IF(MSR)WRITE(*,*)'KMPSI=', KMPSI
   
   READ(BFI,1080,ERR=10100)  O2CRITSI,JSIDETR
   IF(MSR)WRITE(*,*)'O2CRITSI=', O2CRITSI
   IF(MSR)WRITE(*,*)'JSIDETR=', JSIDETR  
   
   READ(BFI,1080,ERR=10100)  DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRITPO4,KMO2DP
   IF(MSR)WRITE(*,*)'DPIE1PO4F=', DPIE1PO4F
   IF(MSR)WRITE(*,*)'DPIE1PO4S=', DPIE1PO4S   
   IF(MSR)WRITE(*,*)'PIE2PO4=', PIE2PO4
   IF(MSR)WRITE(*,*)'O2CRITPO4=', O2CRITPO4
   IF(MSR)WRITE(*,*)'KMO2DP=', KMO2DP   
   
   READ(BFI,1080,ERR=10100)  TEMPBEN,KBENSTR,KLBNTH,DPMIN
   
   IF(MSR)WRITE(*,*)'TEMPBEN=', TEMPBEN
   IF(MSR)WRITE(*,*)'KBENSTR=', KBENSTR   
   IF(MSR)WRITE(*,*)'KLBNTH=', KLBNTH
   IF(MSR)WRITE(*,*)'DPMIN=', DPMIN
   
   
   READ(BFI,1080,ERR=10100)  KAPPCH4,THTACH4,KMCH4O2,KMSO4

   IF(MSR)WRITE(*,*)'KAPPCH4=', KAPPCH4
   IF(MSR)WRITE(*,*)'THTACH4=', THTACH4   
   IF(MSR)WRITE(*,*)'KMCH4O2=', KMCH4O2
   IF(MSR)WRITE(*,*)'KMSO4=', KMSO4
   
   !
   !read information of water colum time-varying conditions (forcing with water column input file)
   !

   IF(ITVWCLF>=1)THEN
        READ (BFI,1090,ERR=10100) (WCLFN(J),J=1,NTVWCLF)

        IF(MSR)THEN
            DO J=1,NTVWCLF
                WRITE(*,*)'WCLFN(', J, ')=', WCLFN(J)
            ENDDO
        ENDIF
   ENDIF
   
   M1_SED_DF=M1
   M2_SED_DF=M2  !will be used by mod_df


    READ(BFI,*)
    READ(BFI,*) 
     READ(BFI,*,ERR=10100) WSSBNET , WSLBNET,WSRBNET,WS1BNET,WS2BNET,WS3BNET,WSUBNET            !settling rate of particulate biogenic silica (m/d)
                              
    IF(MSR)WRITE(*,*)'WSSBNET=', WSSBNET
    IF(MSR)WRITE(*,*)'WSLBNET=', WSLBNET
    IF(MSR)WRITE(*,*)'WSRBNET=', WSRBNET
    IF(MSR)WRITE(*,*)'WS1BNET=', WS1BNET
    IF(MSR)WRITE(*,*)'WS2BNET=', WS2BNET
    IF(MSR)WRITE(*,*)'WS3BNET=', WS3BNET
    IF(MSR)WRITE(*,*)'WSUBNET=', WSUBNET

! burial and mixing rates
   READ(BFI,1080,ERR=10100)   POC1R  !mgC/gSediment
   IF(MSR)WRITE(*,*)'POC1R=', POC1R
   
   READ(BFI,1060,ERR=10100)  SPVARB,PRINTB
   IF(MSR)WRITE(*,*)'SPVARB=', SPVARB
   IF(MSR)WRITE(*,*)'PRINTB=', PRINTB
   
   
   IF (SPVARB == 'CONSTANT') THEN
     READ(BFI,1022,ERR=10100) VSED(1),VPMIX(1),VDMIX(1)
     DO I=2,MTLOC
       VSED(I)=VSED(1)
       VPMIX(I)=VPMIX(1)
       VDMIX(I)=VDMIX(1)
     ENDDO
     
      IF(MSR)WRITE(*,*)'VSED(1)=',  VSED(1)
     IF(MSR)WRITE(*,*)'VPMIX(1)=', VPMIX(1)
     IF(MSR)WRITE(*,*)'VDMIX(1)=', VDMIX(1)
     

   ELSE
     ALLOCATE(RTMP11(MGL));     RTMP11 = 0.0
     ALLOCATE(RTMP12(MGL));     RTMP12 = 0.0
     ALLOCATE(RTMP13(MGL));     RTMP13 = 0.0
     DO I=1,MGL
       READ(BFI,1022,ERR=10100) RTMP11(I),RTMP12(I),RTMP13(I)
     ENDDO  
     IF(SERIAL)THEN
       VSED  = RTMP11
       VPMIX = RTMP12
       VDMIX = RTMP13
     ENDIF
#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,MLOC
         VSED(I)  = RTMP11(NGID(I))
         VPMIX(I) = RTMP12(NGID(I))
         VDMIX(I) = RTMP13(NGID(I))
       ENDDO
       DO I=1,NHN
         VSED(I+MLOC)  = RTMP11(HN_LST(I))
         VPMIX(I+MLOC) = RTMP12(HN_LST(I))
         VDMIX(I+MLOC) = RTMP13(HN_LST(I))
       ENDDO
     ENDIF
#    endif
     
      IF(MSR)WRITE(*,*)'VSED(1)=', VSED(1)
     IF(MSR)WRITE(*,*)'VPMIX(1)=', VPMIX(1)
     IF(MSR)WRITE(*,*)'VDMIX(1)=', VDMIX(1)

     DEALLOCATE(RTMP11,RTMP12,RTMP13)
   ENDIF !end of SPVARB switch

! splits of refractory water column into G2, G3 sediments 

   READ(BFI,1060,ERR=10100)  SPVARLR,PRINTLR
   
   IF(MSR)WRITE(*,*)'SPVARLR=', SPVARLR
   IF(MSR)WRITE(*,*)'PRINTLR=', PRINTLR
   
   IF (SPVARLR == 'CONSTANT') THEN
     READ(BFI,1040,ERR=10100)  FRPOP(1,2),FRPOP(1,3)    &
                              ,FRPON(1,2),FRPON(1,3)    &
                              ,FRPOC(1,2),FRPOC(1,3)
                              
    IF(MSR)WRITE(*,*)' FRPOP(1,2)=', FRPOP(1,2)
    IF(MSR)WRITE(*,*)' FRPOP(1,3)=', FRPOP(1,3)
   
    IF(MSR)WRITE(*,*)' FRPON(1,2)=', FRPON(1,2)
    IF(MSR)WRITE(*,*)' FRPON(1,3)=', FRPON(1,3)
   
    IF(MSR)WRITE(*,*)' FRPOC(1,2)=', FRPOC(1,2)
    IF(MSR)WRITE(*,*)' FRPOC(1,3)=', FRPOC(1,3)
    
     DO I=2,MTLOC
       FRPOP(I,2)=FRPOP(1,2)
       FRPOP(I,3)=FRPOP(1,3)
       FRPON(I,2)=FRPON(1,2)
       FRPON(I,3)=FRPON(1,3)
       FRPOC(I,2)=FRPOC(1,2)
       FRPOC(I,3)=FRPOC(1,3)
     ENDDO
   ELSE
     ALLOCATE(RTMP11(MGL));     RTMP11 = 0.0
     ALLOCATE(RTMP12(MGL));     RTMP12 = 0.0
     ALLOCATE(RTMP13(MGL));     RTMP13 = 0.0
     ALLOCATE(RTMP14(MGL));     RTMP14 = 0.0
     ALLOCATE(RTMP15(MGL));     RTMP15 = 0.0
     ALLOCATE(RTMP16(MGL));     RTMP16 = 0.0
     DO I=1,MGL
       READ(BFI,1040,ERR=10100) RTMP11(I),RTMP12(I),RTMP13(I),    &
                                RTMP14(I),RTMP15(I),RTMP16(I)
     ENDDO
     
     IF(SERIAL)THEN
       FRPOP(:,2) = RTMP11
       FRPOP(:,3) = RTMP12
       FRPON(:,2) = RTMP13
       FRPON(:,3) = RTMP14
       FRPOC(:,2) = RTMP15
       FRPOC(:,3) = RTMP16
     ENDIF
#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,MLOC
         FRPOP(I,2) = RTMP11(NGID(I))
         FRPOP(I,3) = RTMP12(NGID(I))
         FRPON(I,2) = RTMP13(NGID(I))
         FRPON(I,3) = RTMP14(NGID(I))
         FRPOC(I,2) = RTMP15(NGID(I))
         FRPOC(I,3) = RTMP16(NGID(I))
       ENDDO
       DO I=1,NHN
         FRPOP(I+MLOC,2) = RTMP11(HN_LST(I))
         FRPOP(I+MLOC,3) = RTMP12(HN_LST(I))
         FRPON(I+MLOC,2) = RTMP13(HN_LST(I))
         FRPON(I+MLOC,3) = RTMP14(HN_LST(I))
         FRPOC(I+MLOC,2) = RTMP15(HN_LST(I))
         FRPOC(I+MLOC,3) = RTMP16(HN_LST(I))
       ENDDO
     ENDIF
#    endif
     
     IF(MSR)WRITE(*,*)' FRPOP(1,2)=', FRPOP(1,2)
     IF(MSR)WRITE(*,*)' FRPOP(1,3)=', FRPOP(1,3)
        IF(MSR)WRITE(*,*)' FRPON(1,2)=', FRPON(1,2)
     IF(MSR)WRITE(*,*)' FRPON(1,3)=', FRPON(1,3)
     IF(MSR)WRITE(*,*)' FRPOC(1,2)=', FRPOC(1,2)
     IF(MSR)WRITE(*,*)' FRPOC(1,3)=', FRPOC(1,3)
    
                               
     DEALLOCATE(RTMP11,RTMP12,RTMP13,RTMP14,RTMP15,RTMP16)
   ENDIF !end of SPVARLR switch
 
   
   READ(BFI,1099,ERR=10100) SWITCH_LAB,FRACL1,FRACL2,FRACL3
   
   
   IF(MSR)WRITE(*,*)'SWITCH_LAB=', SWITCH_LAB
   IF(MSR)WRITE(*,*)'FRACTION Lab1=', FRACL1
   IF(MSR)WRITE(*,*)'FRACTION Lab2=', FRACL2
   IF(MSR)WRITE(*,*)'FRACTION Lab3=', FRACL3
   
!***** Define logical variables

   STEADY_STATE_SED_G3 = SSTATEG3 == 1
   STEADY_STATE_SED_IC = SSTATEIC == 1
   STEADY_STATE_SED_TS = QUASISS  == 1
   STEADY_STATE_SED_AVG = SSTATEAVG == 1
   AGGREGATE_POM_FLUX = SWITCH_LAB == 1


!********************************************************************************
!**                                 Outputs                                    **
!********************************************************************************

   IF (BENTHIC_OUTPUT) THEN
     IF(MSR)THEN
        WRITE(BFO,2000)
2000     FORMAT(///34X,'Sediment-water column linkages and sediment ',      &    !start printing on 4th line
                    'depths and volumes'/)                                        !and add one more line
        WRITE(BFO,2020) HSEDALL
2020     FORMAT(/' ACTIVE LAYER DEPTH ',F8.3,' CM')                                !print on 7th line

        IF (STEADY_STATE_SED_G3) THEN
            WRITE(BFO,2022)
2022         FORMAT(/' STEADY-STATE VALUES OF G3 COMPUTED'/)                        !print on 9th line and add one line
        ELSE
            WRITE(BFO,2025)
2025         FORMAT(/' NO STEADY-STATE VALUES OF G3 COMPUTED'/)            
        ENDIF
        WRITE(BFO,2030)
2030     FORMAT(////33X,'S E D I M E N T   I N I T I A L   C O N D I T ',   &
            'I O N S'/)
        
        WRITE(BFO,2040)  SSNAME(1)    !temperature
2040     FORMAT(//25X,'Sediment initial conditions for ',A20/)        

        WRITE(BFO,2050) (CTEMP(I),I=1,MLOC)              !WLong: need to use global array
2050     FORMAT(13X,3(7X,1PE11.4))        

        WRITE(BFO,2060)  SSNAME(2)    !POP
2060     FORMAT(//25X,'Sediment initial conditions for ',A20/               &
            37X,'G1',22X,'G2',22X,'G3'/)        
        WRITE(BFO,2070) ((CPOP(I,JG),JG=1,3),I=1,MLOC)     !WLong: need to use global array
2070     FORMAT(18X,3(2X,1PE11.4))

        WRITE(BFO,2060)  SSNAME(3)    !PON
        WRITE(BFO,2070) ((CPON(I,JG),JG=1,3),I=1,MLOC)    !WLong: need to use global array
        
        WRITE(BFO,2060)  SSNAME(4)    !POC
        WRITE(BFO,2070) ((CPOC(I,JG),JG=1,3),I=1,MLOC)    !WLong: need to use global array
            
        WRITE(BFO,2060)  SSNAME(5)    !PBS    particulate biogenic silica
        WRITE(BFO,2070) (CPOS(I),I=1,MLOC)                !WLong: need to use global array
        
        WRITE(BFO,2040)  SSNAME(6)    !total PO4 in layer 2
        WRITE(BFO,2050) (PO4T2TM1S(I),I=1,MLOC)            !WLong: need to use global array

        WRITE(BFO,2040)  SSNAME(7)    !total NH4 in layer 2
        WRITE(BFO,2050) (NH4T2TM1S(I),I=1,MLOC)            !WLong: need to use global array

        WRITE(BFO,2040)  SSNAME(8)    !total NO3 in layer 2
        WRITE(BFO,2050) (NO3T2TM1S(I),I=1,MLOC)            !WLong: need to use global array

        WRITE(BFO,2040)  SSNAME(9)    !total H2S in layer 2
        WRITE(BFO,2050) (HST2TM1S(I),I=1,MLOC)            !WLong: need to use global array

        WRITE(BFO,2040)  SSNAME(10)    !total DSIL (dissolved silicate) in layer 2
        WRITE(BFO,2050) (SIT2TM1S(I),I=1,MLOC)            !WLong: need to use global array

        WRITE(BFO,2040)  SSNAME(11)    !benthic stress
        WRITE(BFO,2050) (BENSTRTM1S(I),I=1,MLOC)        !WLong: need to use global array

        WRITE(BFO,2040)  SSNAME(13)    !total CH4 in layer 2
        WRITE(BFO,2050) (CH4T2TM1S(I),I=1,MLOC)            !WLong: need to use global array

        WRITE(BFO,2040)  SSNAME(14)    !total SO4 in layer 2
        WRITE(BFO,2050) (SO4T2TM1S(I),I=1,MLOC)            !WLong: need to use global array
        
        WRITE(BFO,2080)  10000.0*DIFFT                    !convert from m^2/s to cm^2/sec
2080     FORMAT(//30X,'Temperature diffusion coefficient ',E10.3,               &
            ' cm**2/sec')

        WRITE(BFO,2090)  SALTSW,SALTND
2090     FORMAT(//31X,'If salinity < ',F10.3,' ppt, methane formed',/           &
                 31X,'If salinity < ',F10.3,' ppt, high nit/denit used')
        WRITE(BFO,2100)
2100     FORMAT(//30X,'Particulate organic matter G-model splits'/              &
            10X,'fraction of....',5X,'recycled to',5X,'G1',5X,'G2',             &
            5X,'G3')        
        WRITE(BFO,2110)  FRNAME(1)    ,FRPALG1
        WRITE(BFO,2110)  FRNAME(2)    ,FRPALG2
        WRITE(BFO,2110)  FRNAME(3)    ,FRPALG3
        WRITE(BFO,2110)  FRNAME(5)    ,FRNALG1
        WRITE(BFO,2110)  FRNAME(6)    ,FRNALG2
        WRITE(BFO,2110)  FRNAME(7)    ,FRNALG3
        WRITE(BFO,2110)  FRNAME(9)    ,FRCALG1
        WRITE(BFO,2110)  FRNAME(10)    ,FRCALG2
        WRITE(BFO,2110)  FRNAME(11)    ,FRCALG3
2110     FORMAT(6X,A24,11X,3F7.2)

        WRITE(BFO,2250) (I    ,FRPOP(I,2),FRPOP(I,3)                            &
                            ,FRPON(I,2),FRPON(I,3)                              &
                            ,FRPOC(I,2),FRPOC(I,3)                              &
                            ,I=1,MLOC                                           &    !WLong: need to use global array
                        )                                                            
2250     FORMAT(//31X,'G2 - G3 splits for Refractory Particulates'/             &
            '  NODE   FRG2P   FRG3P   FRG2N   FRG3N   FRG2C   FRG3C'/           &
            (I5,6F8.3))

        WRITE(BFO,2120)  (KPDIAG(JG),DPTHTA(JG),JG=1,3)                         &
                        ,(KNDIAG(JG),DNTHTA(JG),JG=1,3)                         &
                        ,(KCDIAG(JG),DCTHTA(JG),JG=1,3)                         &
                        ,KSI,THTASI
2120     FORMAT(    //30X,'Diagenesis rates (1/day) | Temp corr factor'/        &
                30X,'Phosphorus'/                                               &
                39X,'G1',E11.3,5X,F7.3/                                         &
                39X,'G2',E11.3,5X,F7.3/                                         &
                39X,'G3',E11.3,5X,F7.3/                                         &
                30X,'Nitrogen'/                                                 &
                39X,'G1',E11.3,5X,F7.3/                                         &
                39X,'G2',E11.3,5X,F7.3/                                         &
                39X,'G3',E11.3,5X,F7.3/                                         &
                30X,'Carbon'/                                                   &
                39X,'G1',E11.3,5X,F7.3/                                         &
                39X,'G2',E11.3,5X,F7.3/                                         &
                39X,'G3',E11.3,5X,F7.3/                                         &
                30X,'Silica'/                                                   &
                41X,E11.3,5X,F7.3)                        

        WRITE(BFO,2170)  M1,M2,THTADP,THTADD
2170     FORMAT(    //35X,'Additional constants'/                               &
                30X,'M1........',F8.2,' kg/L'/                                  &
                30X,'M2........',F8.2,' kg/L'/                                  &
                30X,'THTADP....',F8.3,/                                         &
                30X,'THTADD....',F8.3)        
        WRITE(BFO,2180)  KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2
2180     FORMAT(    30X,'KAPPNH4F..',F8.3,' m/day'/                             &
                30X,'KAPPNH4S..',F8.3,' m/day'/                                 &
                30X,'PIENH4....',F8.3,' L/kg'/                                  &
                30X,'THTANH4...',F8.3,/                                         &
                30X,'KMNH4.....',F8.3,' mg n/m**3'/                             &
                30X,'KMNH4O2...',F8.3,' mg o2/L')        
        WRITE(BFO,2190)  KAPPNO3F,KAPPNO3S,K2NO3,THTANO3
2190     FORMAT(    30X,'KAPPNO3F..',F8.3,' m/day'/                             &
                30X,'KAPPNO3S..',F8.3,' m/day'/                                 &
                30X,'K2NO3.....',F8.3,' /day'/                                  &
                30X,'THTANO3...',F8.3)        
        WRITE(BFO,2200)  KAPP1HSD,KAPP1HSP,PIE1HS,PIE2HS,THTAH2S,KMHSO2
2200     FORMAT(    30X,'KAPP1HSD....',F8.3,' m/day'/                           &
                30X,'KAPP1HSP....',F8.3,' m/day'/                               &
                30X,'PIE1HS......',F8.3,' L/kg'/                                &
                30X,'PIE2HS......',F8.3,' L/kg'/                                &
                30X,'THTAH2S...',F8.3,/                                         &
                30X,'KMHSO2....',F8.3,' mg o2/L')        
        WRITE(BFO,2210)  CSISATT20,DPIE1SI,PIE2SI,KMPSI,O2CRITSI,JSIDETR
2210     FORMAT(    30X,'CSISATT20..',F8.1,' mg si/m**3'/                       &
                30X,'DPIE1SI....',F8.3,' L/kg'/                                 &
                30X,'PIE2SI....',F8.3,' L/kg'/                                  &
                30X,'KMPSI.....',E8.2,' mg si/m**3'/                            &
                30X,'O2CRITSI..',E8.2,' mg O2/L'/                               &
                30X,'JSIDETR...',E8.2,' mg Si/m2-d')        
        WRITE(BFO,2220)  DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRITPO4,KMO2DP
2220     FORMAT(    30X,'DPIE1PO4F..',F8.3,' l/kg'/                             &
                30X,'DPIE1PO4S..',F8.3,' l/kg'/                                 &
                30X,'PIE2PO4..',F8.3,' l/kg'/                                   &
                30X,'O2CRITPO4.',F8.3,' mg o2/l'/                               &
                30X,'KMO2DP....',F8.3,' mg o2/l')        
        WRITE(BFO,2230)  TEMPBEN,KBENSTR,KLBNTH,DPMIN
2230     FORMAT(    30X,'TEMPBEN...',F8.3,' deg c'/                             &
                30X,'KBENSTR...',F8.3,' /day'/                                  &
                30X,'KLBNTH....',F8.3,'---'/                                    &
                30X,'DPMIN.....',F8.3,' m2/d')                   
        WRITE(BFO,2240)  KAPPCH4,THTACH4
2240     FORMAT(    30X,'KAPPCH4...',F8.3,' m/day'/                             &
                30X,'THTACH4...',F8.3)        

        WRITE(BFO,2130) I,WSSBNET,WSLBNET,WSRBNET,                    &
                        WS1BNET,WS2BNET,WS3BNET,                       &
                        WSUBNET    !WLong: need to use global array
2130     FORMAT(//3X,'NODE',6X,'WSSNET',4X,'WSLNET',4X,'WSRNET',4X,             &
            'WS1NET',4X,'WS2NET',4X,'WS3NET',4X,'WSUNET'/                       &
            (I7,7F10.3))        
        WRITE(BFO,2140) (I,VSED(I),I=1,MLOC)    !WLong: need to use global array
        
2140     FORMAT(//31X,'Sedimentation (burial) rates (cm/yr)'/                       &
            10X,8(I5,F6.2))        
        WRITE(BFO,2150) (I,VPMIX(I),I=1,MLOC)    !WLong: need to use global array
2150     FORMAT(//30X,'Sediment solid-phase mixing rates (m**2/day)'/               &
            10X,8(I5,F6.2))        
        WRITE(BFO,2160) (I,VDMIX(I),I=1,MLOC)    !WLong: need to use global array
2160     FORMAT(//30X,'Sediment dissolved-phase mixing rates (m**2/day)'/           &
            10X,8(I5,F12.6))        
        CLOSE(BFO)
     ENDIF
   ENDIF
   
!***** Convert cell heights and burial velocities to sediment units

   !MTLOC    --include the halo nodes
   !MLOC    --does not include halo
   
   DO I=1,MLOC
     HSED(I) = HSEDALL*0.01                !WLong convert from cm to m
     VSED(I) = VSED(I)*2.73791E-5        !cm/yr to m/day = 0.01/365.0 = 2.73791E-5
     
 
     !1cm/yr --> 0.01m/365.25 (m/day)
     
     
     
   ENDDO
   
   !WLong and LB: Should we do exchange here?????
   !Answer: this gather() calls to HSED and HSED1 etc in main program does not
   !need halo node information, hence it is not needed
   
   
   !Open the first time variable water column condition file if need be
   !and the read past the header lines
   IF(ITVWCLF==1)THEN
     WCLPTR = 1
     OPEN (WCL,FILE=WCLFN(WCLPTR),STATUS='OLD')
     READ (WCL,1000)    !read past the 3 header lines

     !READ (WCL,1010)     NXWCL,        &!jday of next record
    !                            JCIN_R2,    &!C OM flux
    !                            JNIN_R2,    &!N OM flux
    !                            JPIN_R2,    &!P OM flux
    !                            JSIN_R2,    &!Si OM flux
    !                            O20_R2,        &!O20
    !                            D_R2,        &!total water depth
    !                            TW_R2,        &!water temperature
    !                            NH30_R2,    &!nitrate
    !                            NO30_R2,    &!ammonia
    !                            PO40_R2,    &!phosphate concentration 
    !                            SIAT0_R2,    &!silicate
    !                            CH40_R2,    &!methane
    !                            SALT0_R2     !salinity

    !    IF(MSR)THEN
    !        WRITE(*,1010) NXWCL,        &!jday of next record
    !                            JCIN_R2,    &!C OM flux
    !                            JNIN_R2,    &!N OM flux
    !                            JPIN_R2,    &!P OM flux
    !                            JSIN_R2,    &!Si OM flux
    !                            O20_R2,        &!O20
    !                            D_R2,        &!total water depth
    !                            TW_R2,        &!water temperature
    !                            NH30_R2,    &!nitrate
    !                            NO30_R2,    &!ammonia
    !                            PO40_R2,    &!phosphate concentration 
    !                            SIAT0_R2,    &!silicate
    !                            CH40_R2,    &!methane
    !                            SALT0_R2     !salinity
    !                            
    !        WRITE(*,*)'here I ''m doing nothing '
    !        READ(*,*)
     !   ENDIF
        
1000 FORMAT(//)        !Note only need two slashes to move to the third line! (Weird here)
1010 FORMAT(14(F12.5,1X))     
   ENDIF
   
   !CALL SED_INIT2

   RETURN

!***** Error traps

10100 CONTINUE

   IF (BENTHIC_OUTPUT) THEN
        WRITE(BFO,3010)
        CLOSE(BFO)        
   ENDIF
3010 FORMAT(/' Read error in sediment input deck')        
   STOP 'SED_READ'
    
   RETURN
   END SUBROUTINE SED_READ
    
!********************************************************************************
!**                             Initializations                                **
!********************************************************************************

   !set the initial condition to sediment if it is read from ICI file
   SUBROUTINE SED_INIT_ICI
   
   USE MOD_WQM, ONLY : C2
   
   IMPLICIT NONE
   SAVE
   INTEGER :: I, JG
   
    !if(MSR)write(*,*)'Init temp sedim = ',CTEMPI,'*******************************************'
         DO I=1,MLOC

        !initialize temperature of sediments from water column instead
            CTEMP(I)=T(I,KBM1)  !take bottom layer of water column instead 
        
!#if defined (test2a) || defined (test2b) ||defined(test2c) ||defined(test2d) ||defined (test2e)       
!            !initial condition of concentrations in sediments
!            CTEMP(I)= CTEMPI
!#endif
            
           !DO JG=1,3
             !CPOP(I,JG)  = 100!CPOPI(JG)
             !CPON(I,JG)  = 1000!CPONI(JG)
             !CPOC(I,JG)  = 10000!CPOCI(JG)
           !ENDDO
           !DO JG=1,3   ! inputs are from sediment flux model output for R-64
			
                    CPOC(:,1)= 486*1000 !988*1000!!mgC/m^3     
                    CPOC(:,2)= 2.4256E3*1000 !mgC/m^3 1.93E31.93E3
                    CPOC(:,3)= 7.1432E3*1000       !mgC/m^3
                    CPON(:,1)=85.967*1000           !mgN/m^3
                    CPON(:,2)=427.111*1000          !mgN/m^3
                    CPON(:,3)=947.25*1000       !mgN/m^3
                    CPOP(:,1)=11.92*1000                  !mgP/m^3
                    CPOP(:,2)=44.198*1000          !mgP/m^3
                    CPOP(:,3)=93.895*1000  
           !ENDDO
					
					!CPOC=CPOC/10
					!CPON=CPON/10
					!CPOP=CPOP/10
	
		   !if(NGID(I) == 352)write(*,*)' CPON 1 inital = ',CPON(352,1)/1000, ' g N m^-3 '
           !if(NGID(I) == 352) write(*,*)' CPON 2 inital = ',CPON(352,2)/1000, ' g N m^-3'
		   !if(NGID(I) == 352)write(*,*)' CPON 3 inital = ',CPON(352,3)/1000, ' g N m^-3 '
		   !if(NGID(I) == 352)write(*,*)' CPON total inital = ',sum(CPON(352,:))/1000, ' g N m^-3 '	
		   
           !if(NGID(I) == 352) write(*,*)' CPOC 1 inital = ',CPOC(352,1)/1000, ' g C m^-3'
           !if(NGID(I) == 352) write(*,*)' CPOC 2 inital = ',CPOC(352,2)/1000, ' g C m^-3'		   
           !if(NGID(I) == 352) write(*,*)' CPOC 3 inital = ',CPOC(352,3)/1000, ' g C m^-3'
		   !if(NGID(I) == 352)write(*,*)' CPOC total inital = ',sum(CPOC(352,:))/1000, ' g N m^-3 '	
		   
           CPOS(I)       = CPOSI
           
           PO41TM1S(I)  = PO41TI
           PO4T2TM1S(I)  = PO4T2I
           
           NH41TM1S(I)     = NH41TI
           NH4T2TM1S(I)  = NH4T2I
           
           NO31TM1S(I)   = NO31TI
           NO3T2TM1S(I)  = NO3T2I
           
           HS1TM1S(I)     = HS1TI
           HST2TM1S(I)   = HST2I
           
           CH41TM1S(I)   = CH41TI
           CH4T2TM1S(I)  = CH4T2I
           
           !SO41TM1S(I)   = SO41TI 
           SO4T2TM1S(I)  = SO4T2I
           
           SI1TM1S(I)       = SI1TI
           SIT2TM1S(I)   = SIT2I
           
           BENSTRTM1S(I)   = BENSTI
           SODTM1S(I)     = SODI
           
    !WLong: may also need to do the following to initialize 2nd layer dissolved concentrations
    !these are not needed for equation solver to march to next time step, but they might be 
    !written to output files at the beginning of the model (before first time step finished)
    !Answer: actually in wqm_main.F these is written to outputs after call to SED_CALC is done
    !       hence it is not a problem.
    !
           !PIE1=PIENH4
           !PIE2=PIENH4
           !FD1=1/(1.0+m1*PIE1)
           !FD2=1/(1.0+m2*PIE2)
           !NH42I=NH4T2I*FD2
           
           !NO32I=  NO3T2I
           
           !IF(SALT(I,KBM1)<SALTND)THEN
            !        DPIE1PO4=DPIE1PO4F
           !ELSE
            !        DPIE1PO4=DPIE1PO4S
           !ENDIF
           !IF(O20>O2CRITPO4)THEN
            !        PIE1=PIE2PO4*DPIE1PO4
           !ELSE
            !        PIE1=PIE2PO4*DPIE1PO4**(O20/O2CRITPO4)
           !ENDIF
           !PIE2=PIE2PO4
           !FD1=1/(1.0+m1*PIE1)
           !FD2=1/(1.0+m2*PIE2)
           !PO42I=PO4T2I*FD2
           
           !IF(O20>O2CRITSI)THEN
           !        PIE1=PIE2SI*DPIE1SI
           !ELSE
           !        PIE1=PIE2SI*DPIE1SI**(O20/O2CRITSI)
           !ENDIF
           !PIE2=PIE2SI
           !FD1=1/(1.0+m1*PIE1)
           !FD2=1/(1.0+m2*PIE2)
           !SI2I= SIT2I*FD2
           !          
           !PIE1=PIE1HS
           !PIE2=PIE2HS
           !FD1=1/(1.0+m1*PIE1)
           !FD2=1/(1.0+m2*PIE2)
           !HS2I= HST2I *FD2

    !Wen Long: we should move the following to mod_ba       
           !IF(BALGAE_CALC)THEN
        !        !biomass of benthic algae
        !        BBM(I)        = BBMI
           !ENDIF
           
    !Wen Long: we should move the following to mod_sf
           !!biomass of suspension feeders
           !IF(SFEEDER)THEN
        !        DO N=1,NSPECIES
        !            SFEED(I,N)      = SFEEDI(N)
        !        ENDDO
           !ENDIF
        !    
           
    !Wen Long: we should move the following to mod_df
           !!biomass of deposition feeder
           !IF(DFEEDER)THEN
        !        DFEEDM1S(I)   = DFEEDI
           !ENDIF
         
         ENDDO

         
END SUBROUTINE SED_INIT_ICI
   
   
SUBROUTINE SED_INIT2
   
   USE MOD_SED_DF_EXCHANGE_VARS, ONLY:  &
        POC1TM1S_SED_DF    ,            &
        POC2TM1S_SED_DF    ,            &
        M1_SED_DF,                      &
        M2_SED_DF

#  if defined (MULTIPROCESSOR)
        USE MOD_PAR, ONLY:   &      
            !EL_PID,    &        !(:)           !!PROCESSOR OWNER OF GLOBAL ELEMENT
            !ELID,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT
            !NLID,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE 
            !ELID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
            !NLID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
            !OBN_GL2LOC,    &    !(:)           !!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
            !EGID,        &    !(:)           !!GLOBAL ID OF LOCAL ELEMENT
            NGID,        &    !(:)           !!GLOBAL ID OF LOCAL NODE 
            !NHE,        &   !            !NUMBER OF HALO ELEMENTS
            !HE_LST,        &    !(:)           !!GLOBAL IDENTITIES OF HALO ELEMENTS
            !HE_OWN,        &    !(:)           !!OWNER OF HALO ELEMENTS
            !NBN ,       &   !              !!NUMBER OF BOUNDARY NODES
            !MX_MLT,     &   !              !!MAX MULTIPLICITY OF BOUNDARY NODES/
            !BN_LST,     &    !               !!GLOBAL IDENTITY OF BOUNDARY NODES
            !BN_LOC,     &    !            !!LOCAL IDENTITY OF BOUNDARY NODES
            !BN_MLT,     &    !            !!MULTIPLICITY OF BOUNDARY NODES
            !BN_NEY,     &    !            !!NODE OWNER LIST
            !NDE_ID,     &    !            !! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
            NHN,          &    !            !!NUMBER OF HALO NODES
            HN_LST!,     &    !            !!LIST OF HALO NODES 
            !HN_OWN,     &    !            !!PRIMARY OWNER OF HALO NODES
            !EC,            &    !
            !NC,            &    !
            !BNC,         &    !
            !EMAP,        &    !
            !NMAP,         &    !
            !PNE,        &    !(:)        !!NUMBER OF ELEMENTS IN EACH PROC
            !PNN,        &    !(:)        !!NUMBER OF NODES IN EACH PROC
            !PNHE,        &    !(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
            !PNBN,        &    !(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
            !PMBM,        &    !(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
            !PNHN,        &    !(:)       !!NUMBER OF HALO NODES IN EACH PROC
            !FILELOCK,   &
            !MYFILELOCK, &
            !FILELOCKED , &
            !NODE_MATCH, &
            !NODE_MATCH_ALL, &
            !EXCHANGE, &
            !EXCHANGE_ALL, &
            !, &
            !IGATHER, &
            !INIT_MPI_ENV,&!            
            !SORT,    &
            !PREAD,    &
            !PWRITE,&
            !PPRINT, &
            !LOCK, &
            !UNLOCK, &
            !GETLOC
#  endif  

   IMPLICIT NONE
   SAVE    
    INTEGER :: I, JG, JT

    
    
    
!***** Set sediment concentrations to initial concentrations
 
    
    
   DO I=1,MLOC
     POP1TM1S(I) = CPOP(I,1)
     POP2TM1S(I) = CPOP(I,2)
     POP3TM1S(I) = CPOP(I,3)
     PON1TM1S(I) = CPON(I,1)
     PON2TM1S(I) = CPON(I,2)
     PON3TM1S(I) = CPON(I,3)
     POC1TM1S(I) = CPOC(I,1)
     POC2TM1S(I) = CPOC(I,2)
     POC3TM1S(I) = CPOC(I,3)
     PSISEDTM1S(I)  = CPOS(I)

     !WLong added POC1TM1S_SED_DF and POC2TM1S_SED_DF
     IF(DFEEDER)THEN
        POC1TM1S_SED_DF(I)=CPOC(I,1)
        POC2TM1S_SED_DF(I)=CPOC(I,2)
        POC2TM1S_SED_DF(I)=CPOC(I,2)
     ENDIF
                
   ENDDO

!***** Initialize mass balance variables (kgC, kgN, kgP)
    ISEDMC = 0.
    ISEDMN = 0.
    ISEDMP = 0.
    

    DO I=1,MLOC
     CPO4(I) = PO4T2TM1S(I)
     CNO3(I) = NO3T2TM1S(I)
     CNH4(I) = NH4T2TM1S(I)

     !initial total C, N, P in sediments (kg)
     ISEDMC   = ISEDMC+(CPOC(I,1)+CPOC(I,2)+CPOC(I,3))             &
                *ART1(I)*HSED(I)/1.E6                               
                    !unit of ISEDMC is mgC/m^3 * m^2 * m /1000/1000 
                    !==> mgC/1000/1000 ==>gC/1000 ==> kgC            
      
     ISEDMN   = ISEDMN+(CPON(I,1)+CPON(I,2)+CPON(I,3)+CNH4(I)      &
               +CNO3(I))*ART1(I)*HSED(I)/1.E6                       
                    !unit of ISEDMN is mgN/m^3 * m^2 * m /1000/1000 
                    !==> mgN/1000/1000 ==>gN/1000 ==> kgN
     ISEDMP   = ISEDMP+(CPOP(I,1)+CPOP(I,2)+CPOP(I,3)+CPO4(I))     &  !WLong, why not including CPO4??
                *ART1(I)*HSED(I)/1.E6                    
                    !unit of ISEDMP is mgP/m^3 * m^2 * m /1000/1000 
                    !==> mgP/1000/1000 ==>gP/1000 ==> kgP
    ENDDO

!***** Set up reaction rates in table look-up form

   DO JT=1,350
     TEMP_SED     = REAL(JT-1)/10.+0.05
     TEMP20       = TEMP_SED-20.
     TEMP202      = TEMP20/2.
     
     ZHTANH4F(JT) = KAPPNH4F*THTANH4**TEMP202
     ZHTANH4S(JT) = KAPPNH4S*THTANH4**TEMP202
     ZHTA1HSD(JT) = KAPP1HSD*THTAH2S**TEMP202
     ZHTA1HSP(JT) = KAPP1HSP*THTAH2S**TEMP202
     ZHTANO3F(JT) = KAPPNO3F*THTANO3**TEMP202
     ZHTANO3S(JT) = KAPPNO3S*THTANO3**TEMP202
     
     ZHTAK2NO3(JT) = K2NO3*THTANO3**TEMP20        !m/d
     
     ZL12NOM(JT)  = THTADD**TEMP20
     ZW12NOM(JT)  = THTADP**TEMP20  !NO UNIT
     
     ZHTAPON1(JT) = KPON1*THTAPON1**TEMP20
     ZHTAPON2(JT) = KPON2*THTAPON2**TEMP20
     ZHTAPON3(JT) = KPON3*THTAPON3**TEMP20
     ZHTAPOC1(JT) = KPOC1*THTAPOC1**TEMP20
     ZHTAPOC2(JT) = KPOC2*THTAPOC2**TEMP20
     
        IDEBUG_SED=0
          !IDEBUG_SED=1
     
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'KPOC2=',KPOC2
        !    WRITE(*,*)'THTAPOC2=',THTAPOC2
        !    WRITE(*,*)'TEMP20=',TEMP20
        !    WRITE(*,*)'JT=',JT
        !    READ(*,*)
        !ENDIF                

     
     ZHTAPOC3(JT) = KPOC3*THTAPOC3**TEMP20
     ZHTAPOP1(JT) = KPOP1*THTAPOP1**TEMP20
     ZHTAPOP2(JT) = KPOP2*THTAPOP2**TEMP20
     ZHTAPOP3(JT) = KPOP3*THTAPOP3**TEMP20
     ZHTASI(JT)   = KSI*THTASI**TEMP20            !
     ZHTACH4(JT)  = KAPPCH4*THTACH4**TEMP202
     
   ENDDO

!***** Turn off settling

   IF (.NOT.SETTLING) THEN

       WSSBNET = 0.
       WSLBNET = 0.
       WSRBNET = 0.
       WS1BNET = 0.
       WS2BNET = 0.
       WS3BNET = 0.
       WSUBNET = 0.

   ENDIF

!***** Initialize accumulators for steady-state computations

   IF (STEADY_STATE_SED_G3) THEN
     TINTIM = 0.
    !WLong, initialize TINTIM_BA for benthic algae
    !WLong moved these to mod_ba.F
    ! IF(BALGAE_CALC)THEN
    !    TINTIM_BA = 0.
    ! ENDIF
     DO I=1,MLOC
       AG3CFL(I) = 0.
       AG3NFL(I) = 0.
       AG3PFL(I) = 0.
       ASDTMP(I) = 0.
     ENDDO
   ENDIF

   RETURN

   END SUBROUTINE SED_INIT2
   
   
     
!********************************************************************************
!**              S U B R O U T I N E   P O M _ A C C U M U L                   **
!**      Accumulate POM from water column, before reaching the sedims          **
!********************************************************************************

   SUBROUTINE POM_ACCUMUL
   
      USE MOD_SAV, ONLY:  NSAVCELL,        &
                          SAVCELL,         &
                          NSAVSPC,         &
                          LEAF,            &
                          STEM,            &
                        SAVFRAC,           &
                        WSSSAV,            &
                        WSLSAV,            &
                        WSRSAV,            &
                        WS1SAV,            &
                        WS2SAV,            &
                        WS3SAV,            &
                        WSUSAV,            &
                        SEDPOCSAV,         &
                        FRPOCSAV,          &
                        SEDPONSAV,         &
                        FRPONSAV,          &
                        SEDPOPSAV,         &
                        FRPOPSAV,          &
                        SEDNH4SAV,         &
                        SEDNO3SAV,         &
                        SEDDOSAV,          &
                        SEDPO4SAV  
      IMPLICIT NONE
   
      INTEGER :: I, N
      REAL(SP) :: SAVEFCT
   !******* Assign base net settling rates

!   IDEBUG=1
!   IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
   

!		
          WSSNET = WSSBNET
          WSLNET = WSLBNET
          WSRNET = WSRBNET
          WS1NET = WS1BNET
          WS2NET = WS2BNET
          WS3NET = WS3BNET
          WSUNET = WSUBNET

!******* Adjust net settling for SAV effect
!commented out by B Clark. If SAV, settling will be accounted for due to increased drag
          !and decreased velocity
   ! IF (SAV_CALC) THEN
   !   DO I=1,NSAVCELL
   !     B=SAVCELL(I)
   !     DO N=1,NSAVSPC(B)
   !       SAVEFCT = (LEAF(B,N)+STEM(B,N)) * SAVFRAC(B,N)
   !       WSSNET = WSSNET + WSSSAV * SAVEFCT  !increase of settling rate (m/d) to suspended solids
   !       WSLNET = WSLNET + WSLSAV * SAVEFCT  !increase of settling rate (m/d) to LPOM
   !       WSRNET = WSRNET + WSRSAV * SAVEFCT  !increase of settling rate (m/d) to RPOM
   !       WS1NET = WS1NET + WS1SAV * SAVEFCT  !increase of settling rate (m/d) to alg 1
   !       WS2NET = WS2NET + WS2SAV * SAVEFCT  !increase of settling rate (m/d) to alg 2
   !       WS3NET = WS3NET + WS3SAV * SAVEFCT  !increase of settling rate (m/d) to alg 3
   !       WSUNET = WSUNET + WSUSAV * SAVEFCT  !increase of settling rate (m/d) to particulate biogenic silicate (unavaiable)
   !     ENDDO
   !   ENDDO
   ! ENDIF

!   IDEBUG=3
!   IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
   
!***** Calculate fluxes

   DO I=1,MLOC
     KWC = KBM1

!******* Settling rate (LHS has unit mm/day) ****

     WSSINETMMD    = 1000.*WSSNET
     WLPOMNETMMD   = 1000.*WSLNET
     WRPOMNETMMD   = 1000.*WSRNET
     WB1NETMMD     = 1000.*WS1NET
     WB2NETMMD     = 1000.*WS2NET
     WB3NETMMD     = 1000.*WS3NET
     WPSINETMMD    = 1000.*WSUNET

!*** Fluxes***(mg/m^2/d) (positive into sediments)***
     
    IF(ITVWCLF/=1)THEN !use settling calculated from overlying water settling
                 
        !G1 flux of POC, Unit : [(mm/d)*(1)*gC/m^3 ] ==> (mgC/m^2/day)                   
        JPOCaccum(I,1) = WB1NETMMD*FRCALG1(1)*B1(I,KWC)                    &
                   +WB2NETMMD*FRCALG2(1)*B2(I,KWC)                         &
                   +WB3NETMMD*FRCALG3(1)*B3(I,KWC)                         &
                   +WLPOMNETMMD*LPOC(I,KWC)                                &!test case has LPOC=0
                   + JPOCaccum(I,1)                         !full channel has LPOC meaningful
                       
            IF(SFEEDER)THEN
                JPOCaccum(I,1) =JPOCaccum(I,1) +JLPOC_SED_SF(I)*FRCALG3(1)                !suspension feeders
            ENDIF


        !G2 flux of POC, Unit : [(mm/d)*(1)*gC/m^3 ] ==> (mgC/m^2/day)                                      
        JPOCaccum(I,2) = WB1NETMMD*FRCALG1(2)*B1(I,KWC)                    &
                   +WB2NETMMD*FRCALG2(2)*B2(I,KWC)                         &
                   +WB3NETMMD*FRCALG3(2)*B3(I,KWC)                         &
                   +WRPOMNETMMD*RPOC(I,KWC)*FRPOC(I,2)/                    &
                    (FRPOC(I,2)+FRPOC(I,3))                                &
                    + JPOCaccum(I,2)
                    
        IF(SFEEDER)THEN                    
            JPOCaccum(I,2)=JPOCaccum(I,2)  +JLPOC_SED_SF(I)*FRCALG3(2)                    !suspension feeders
        ENDIF


        !G3 flux of POC, Unit : [(mm/d)*(1)*gC/m^3 ] ==> (mgC/m^2/day)                                      
        JPOCaccum(I,3) = WB1NETMMD*FRCALG1(3)*B1(I,KWC)+                   &
                        WB2NETMMD*FRCALG2(3)*B2(I,KWC)                     &
                        +WB3NETMMD*FRCALG3(3)*B3(I,KWC)                    &
                        +WRPOMNETMMD*RPOC(I,KWC)*FRPOC(I,3)/               &
                        (FRPOC(I,2)+FRPOC(I,3))                            &
                        + JPOCaccum(I,3)
                    
        IF(SFEEDER)THEN
            JPOCaccum(I,3) = JPOCaccum(I,3) +JLPOC_SED_SF(I)*FRCALG3(3)                   !suspension feeders
                        !WLong and LB: why should LPOP go here? 
                        !this needs to be sorted out
        ENDIF
                   
        !G1 flux of PON, Unit : [(mm/d)*(gN/gC)*(1)*gC/m^3 ] ==> (mgN/m^2/day)
        JPONaccum(I,1) = WB1NETMMD*ANC1*FRNALG1(1)*B1(I,KWC)               &
                    +WB2NETMMD*ANC2*FRNALG2(1)*B2(I,KWC)                   &
                    +WB3NETMMD*ANC3*FRNALG3(1)*B3(I,KWC)                   &
                    +WLPOMNETMMD*LPON(I,KWC)                               &
                    + JPONaccum(I,1)
                    
        IF(SFEEDER)THEN                   
            JPONaccum(I,1)= JPONaccum(I,1)+JLPON_SED_SF(I)*FRNALG3(1)                     !suspension feeders
        ENDIF                   


        !G2 flux of PON, Unit : [(mm/d)*(gN/gC)*(1)*gC/m^3 ] ==> (mgN/m^2/day)                   
        JPONaccum(I,2) = WB1NETMMD*ANC1*FRNALG1(2)*B1(I,KWC)                       &
                    +WB2NETMMD*ANC2*FRNALG2(2)*B2(I,KWC)                    &
                    +WB3NETMMD*ANC3*FRNALG3(2)*B3(I,KWC)                    &
                    +WRPOMNETMMD*RPON(I,KWC)*FRPON(I,2)/                    &
                        (FRPON(I,2)+FRPON(I,3))                             &
                    + JPONaccum(I,2)
    
        IF(SFEEDER)THEN                    
            JPONaccum(I,2) = JPONaccum(I,2) +JLPON_SED_SF(I)*FRNALG3(2)                   !suspension feeders
        ENDIF                   

        !G3 flux of PON, Unit : [(mm/d)*(gN/gC)*(1)*gC/m^3 ] ==> (mgN/m^2/day)                   
        JPONaccum(I,3) = WB1NETMMD*ANC1*FRNALG1(3)*B1(I,KWC)                       &
                    +WB2NETMMD*ANC2*FRNALG2(3)*B2(I,KWC)                    &
                    +WB3NETMMD*ANC3*FRNALG3(3)*B3(I,KWC)                    &
                    +WRPOMNETMMD*RPON(I,KWC)*FRPON(I,3)/                    &
                        (FRPON(I,2)+FRPON(I,3))                             &
                    + JPONaccum(I,3)   
                    
        IF(SFEEDER)THEN                    
            JPONaccum(I,3)=JPONaccum(I,3) + JLPON_SED_SF(I)*FRNALG3(3)                      !suspension feeders
                                !WLong and LB: why should LPOP go here? 
                                !this needs to be sorted out
        ENDIF                   
#if defined(PHOSPHORUS)     
        !G1 flux of POP, Unit : [(mm/d)*(gP/gC)*(1)*gC/m^3 ] ==> mgP/m^2/day)     
        JPOPaccum(I,1) =  WB1NETMMD*Q1(I,KWC)*FRPALG1(1)*B1(I,KWC)               &
                    +WB2NETMMD*Q2(I,KWC)*FRPALG2(1)*B2(I,KWC)               &
                    +WB3NETMMD*Q3(I,KWC)*FRPALG3(1)*B3(I,KWC)               &
                    +WLPOMNETMMD*LPOP(I,KWC)                                &
                    + JPOPaccum(I,1)
                    
        IF(SFEEDER)THEN
            JPOPaccum(I,1)=JPOPaccum(I,1)+JLPOP_SED_SF(I)*FRPALG3(1)                        !suspension feeders
        ENDIF                   
                    
        !G2 flux of POP, Unit : [(mm/d)*(gP/gC)*(1)*gC/m^3 ] ==> (mgP/m^2/day)
        JPOPaccum(I,2) =  WB1NETMMD*Q1(I,KWC)*FRPALG1(2)*B1(I,KWC)               &
                    +WB2NETMMD*Q2(I,KWC)*FRPALG2(2)*B2(I,KWC)               &
                    +WB3NETMMD*Q3(I,KWC)*FRPALG3(2)*B3(I,KWC)               &
                    +WRPOMNETMMD*RPOP(I,KWC)*FRPOP(I,2)/                    &
                        (FRPOP(I,2)+FRPOP(I,3))                             &
                    + JPOPaccum(I,2)
                    
        IF(SFEEDER)THEN
            JPOPaccum(I,2)=JPOPaccum(I,2)+JLPOP_SED_SF(I)*FRPALG3(2)                        !suspension feeders
        ENDIF                   
        
        !G3 flux of POP, Unit : [(mm/d)*(gP/gC)*(1)*gC/m^3 ] ==> (mgP/m^2/day)
        JPOPaccum(I,3) =  WB1NETMMD*Q1(I,KWC)*FRPALG1(3)*B1(I,KWC)               &
                    +WB2NETMMD*Q2(I,KWC)*FRPALG2(3)*B2(I,KWC)               &
                    +WB3NETMMD*Q3(I,KWC)*FRPALG3(3)*B3(I,KWC)               &
                    +WRPOMNETMMD*RPOP(I,KWC)*FRPOP(I,3)/                    &
                        (FRPOP(I,2)+FRPOP(I,3))                             &
                    + JPOPaccum(I,3)
                    
        IF(SFEEDER)THEN                    
            JPOPaccum(I,3)=JPOPaccum(I,3)+JLPOP_SED_SF(I)*FRPALG3(3)                        !suspension feeders
                                !WLong and LB: why should LPOP go here? 
                                !this needs to be sorted out
        ENDIF                   
#endif
        !Particulate Biogeonic Silica flux to sediments (positive to sediments)
        !Unit: (mm/d) * (gSi/gC) * (g/m^3) ==> mgSi/m^2/d
        JPOSaccum(I) =     WB1NETMMD*ASC1*B1(I,KWC)                         &
                    +WB2NETMMD*ASC2*B2(I,KWC)                               &
                    +WB3NETMMD*ASC3*B3(I,KWC)                               &!WLong we need to add contribution due to
                    +WPSINETMMD*SIUPB(I,KWC)                                &!SIUPB (particulate biogenic)
                    + JPOSaccum(I)                                !Note that SIAT is not included for it's not organic
                    
                    !Should add RPOM contribution to Silicate based Si/C ratio in detritus
                    !and JSIDETR is supposed to caculate that instead of being given from input
                    
        IF(SFEEDER)THEN                    
            JPOSaccum(I)=JPOSaccum(I)+ JSU_SED_SF(I)        !+    SU_SED_SF(I)    !suspension feeder
                                                !WLong commented SU_SED_SF for not understanding it
                                                !It seems arleady accounted for in SF
        ENDIF                    

     ELSE   

            !for test2E with POM from a file 
            !set Overlying water colum settling flux according to the values read fromt the file

            !NXWCL,        &!jday of next record
            !JCIN_R1,    &!C OM flux    (gO2/m^2/day)
            !JNIN_R1,    &!N OM flux (gN/m^2/day)
            !JPIN_R1,    &!P OM flux (gP/m^2/day)
            !JSIN_R1,    &!Si OM flux (gSi/m^2/day)             

            !JPOC(I,1)   = JCIN_R1*FRCALG1(1)*1000.0/2.667  !G1 POC flux convert from gO2/m^2/day to mgC/m^2/day 
            !JPOC(I,2)   = JCIN_R1*FRCALG1(2)*1000.0/2.667  !G2 POC flux convert from gO2/m^2/day to mgC/m^2/day
            !JPOC(I,3)   = JCIN_R1*FRCALG1(3)*1000.0/2.667  !G3 POC flux convert from gO2/m^2/day to mgC/m^2/day
     
            !JPON(I,1)   = JNIN_R1*FRNALG1(1)*1000.0        !G1 PON flux converted from gN/m^2/day to mgN/m^2/day
            !JPON(I,2)   = JNIN_R1*FRNALG1(2)*1000.0        !G2 PON flux converted from gN/m^2/day to mgN/m^2/day
            !JPON(I,3)   = JNIN_R1*FRNALG1(3)*1000.0        !G3 PON flux converted from gN/m^2/day to mgN/m^2/day

            !JPOP(I,1)   = JPIN_R1*FRPALG1(1)*1000.0        !G1 PON flux converted from gN/m^2/day to mgN/m^2/day
            !JPOP(I,2)   = JPIN_R1*FRPALG1(2)*1000.0        !G2 PON flux converted from gN/m^2/day to mgN/m^2/day
            !JPOP(I,3)   = JPIN_R1*FRPALG1(3)*1000.0        !G3 PON flux converted from gN/m^2/day to mgN/m^2/day
            
            !JPOS(I)     = JSIN_R1*1000.0                   !Silicate flux from gSi/m^2/day to mgSi/m^2/day

            JPOCaccum(I,1)   = JCIN_R1*FRCALG1(1)*1000.0/2.667  !G1 POC flux converted from gO2/m^2/day to mgC/m^2/day
            JPOCaccum(I,2)   = JCIN_R1*FRCALG1(2)*1000.0/2.667  !G2 POC flux converted from gO2/m^2/day to mgC/m^2/day
            JPOCaccum(I,3)   = JCIN_R1*FRCALG1(3)*1000.0/2.667  !G3 POC flux converted from gO2/m^2/day to mgC/m^2/day

            JPONaccum(I,1)   = JNIN_R1*FRNALG1(1)*1000.0        !G1 PON flux converted from gN/m^2/day to mgN/m^2/day
            JPONaccum(I,2)   = JNIN_R1*FRNALG1(2)*1000.0        !G2 PON flux converted from gN/m^2/day to mgN/m^2/day
            JPONaccum(I,3)   = JNIN_R1*FRNALG1(3)*1000.0        !G3 PON flux converted from gN/m^2/day to mgN/m^2/day

#if defined(PHOSPHORUS)
            JPOPaccum(I,1)   = JPIN_R1*FRPALG1(1)*1000.0        !G1 POP flux converted from gP/m^2/day to mgP/m^2/day
            JPOPaccum(I,2)   = JPIN_R1*FRPALG1(2)*1000.0        !G2 POP flux converted from gP/m^2/day to mgP/m^2/day
            JPOPaccum(I,3)   = JPIN_R1*FRPALG1(3)*1000.0        !G3 POP flux converted from gP/m^2/day to mgP/m^2/day
#endif


            JPOSaccum(I) = JSIN_R1*1000.0                   !Silicate flux from gSi/m^2/day to mgSi/m^2/day

            !       IF(MSR)THEN                          &!
            !           WRITE(*,*)'JDAY_R1=', JDAY,      &!
            !                     'JCIN_R1=',JCIN_R1,    &!
            !                     'JNIN_R1=',JNIN_R1,    &!
            !                     'JPIN_R1=',JPIN_R1,    &!
            !                     'JSIN_R1=',JSIN_R1
            !       ENDIF

     ENDIF
   ENDDO  
   
   END SUBROUTINE POM_ACCUMUL
     
!********************************************************************************
!**                    S U B R O U T I N E   S E D _ C A L C                   **
!**                           Sediment Calculations                            **
!********************************************************************************

   SUBROUTINE SED_CALC(TF_SSTATE)
   
   USE MOD_SED_DF_EXCHANGE_VARS, ONLY: &
        POC1TM1S_SED_DF,                 &
        POC2TM1S_SED_DF	
        
   USE MOD_SAV, ONLY:  NSAVCELL,        & 
                        SAVCELL,         &
                        NSAVSPC,        &
                        LEAF,            &
                        STEM,             &
                      SAVFRAC,            &
                      WSSSAV,             &
                      WSLSAV,             &
                      WSRSAV,             &
                      WS1SAV,            &
                      WS2SAV,            &
                      WS3SAV,             &
                      WSUSAV,             &
                      SEDPOCSAV,        &
                      FRPOCSAV,            &
                      SEDPONSAV,        &
                      FRPONSAV,            &
                      SEDPOPSAV,        &
                      FRPOPSAV,            &
                      SEDNH4SAV,        &
                      SEDNO3SAV,        &
                      SEDDOSAV,            &
                      SEDPO4SAV  
                      
                    !Wen Long: We are not tracking DIC in water column
                    !          otherwise, we would need SEDDICSAV etc for Ocean Acidification
   IMPLICIT NONE                     
! added next line 10/17/05
   SAVE    
   INTEGER :: I, J, JSF, N, IERR, IDEBUG
   REAL(SP) :: SAVEFCT
   REAL(SP) :: SOD!,ZBOUT
   
   REAL(SP) :: DFEED_TMP, DF_SOD_TMP
   
   REAL(SP) :: SAV_SOD_TMP
   
   LOGICAL :: TF_SSTATE   !Type of integration  FALSE - time intergration from current time step to next time step
                             !                        TRUE  - no time integration, calculate steady steate solution

     
   DLTS = NDTSED * DLT/86400.        !time step in days for sediment flux model - sediment timestep is NDTSED times larger than water column time step (DLT) -added by LB
#if defined(SEDIMENT_DOM)  
   DLTS_SHARE=DLTS !B Clark PASS THE TIME STEP OUT TO THE DOM MODULE
#endif  

   IF (STEADY_STATE_SED_G3) TINTIM = TINTIM+DLTS

!***** Initialize sediment nutrient masses
   
   
   SEDMC = 0.
   SEDMN = 0.
   SEDMP = 0.

     
   
   DO I=1,MLOC   !LB moved DO here because of moved calculations of JPOM
     KWC = KBM1
!
! *** Sum particulate fluxes to water column, negative into sediments (postitive into water)
! used for detecting the mass balance from the water column's point of view
! (LHS is g/m^2/day)
                   
    !Sum of particuplate P, N, C, Particulate Biogenic Si, Suspended Solids to water column (g/m^2/day)
     PPFWS(I)=-0.001*(JPOP(I,1)+JPOP(I,2)+JPOP(I,3))  !converted from mgP/m^2/d --> gP/m^2/day
     PNFWS(I)=-0.001*(JPON(I,1)+JPON(I,2)+JPON(I,3))  !converted from mgN/m^2/d --> gN/m^2/day
     PCFWS(I)=-0.001*(JPOC(I,1)+JPOC(I,2)+JPOC(I,3))  !converted from mgC/m^2/d --> gC/m^2/day
     
     PSFWS(I)=-0.001*JPOS(I)                          !converted from mgSi/m^2/d --> gSi/m^2/day
     
     
        IDEBUG_SED=IDEBUG_SED+1

     
     !Unit: (m/day)*(g/m^3) --> gSolids/m^2/day
     SSFWS(I)=-WSSNET*SSI(I,KWC)
     IF(SFEEDER)THEN
        SSFWS(I)=SSFWS(I)-0.001*SSI_SED_SF(I)   !converted from mgSolids/m^2/day -->gSolids/m^2/day
     ENDIF
     
   ENDDO

!  IDEBUG=4
!  IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
   
! ADD IN THE FLUX FROM ROOT MORTALITY OF SAV

   IF (SAV_CALC) THEN
     DO I=1,NSAVCELL
       B=SAVCELL(I)
       DO J=1,3
         JPOC(B,J) = JPOC(B,J)+1000.*SEDPOCSAV(B)*FRPOCSAV(J)
         JPON(B,J) = JPON(B,J)+1000.*SEDPONSAV(B)*FRPONSAV(J)
         JPOP(B,J) = JPOP(B,J)+1000.*SEDPOPSAV(B)*FRPOPSAV(J)
       ENDDO
     ENDDO
   ENDIF

!    *** Accumulate fluxes for steady-state computation for G3

   IF (STEADY_STATE_SED_G3) THEN
     DO I=1,MLOC
       AG3CFL(I) = AG3CFL(I)+JPOC(I,3)*DLTS
       AG3NFL(I) = AG3NFL(I)+JPON(I,3)*DLTS
       AG3PFL(I) = AG3PFL(I)+JPOP(I,3)*DLTS
     ENDDO
   ENDIF
   
!   IDEBUG=5
!   IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
   

!***** Assign previous timestep concentrations to particulate organics
  !IF(MSR)write(*,*)'LBnote a: CPOC_GL(38,2)=',CPOC_GL(38,2), '  CPOC(38,3)=',CPOC(38,2)
   DO I=1,MLOC
     CPOP(I,1) = POP1TM1S(I)
     CPOP(I,2) = POP2TM1S(I)
     CPOP(I,3) = POP3TM1S(I)
     CPON(I,1) = PON1TM1S(I)
     CPON(I,2) = PON2TM1S(I)
     CPON(I,3) = PON3TM1S(I)
     CPOC(I,1) = POC1TM1S(I)
     CPOC(I,2) = POC2TM1S(I)
     CPOC(I,3) = POC3TM1S(I)
     CPOS(I)   = PSISEDTM1S(I)
     
        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'CPOP(I,1)=',CPOP(I,1)
        !    WRITE(*,*)'CPOP(I,2)=',CPOP(I,2)
        !    WRITE(*,*)'CPOP(I,3)=',CPOP(I,3)
        !    WRITE(*,*)'CPON(I,1)=',CPON(I,1)
        !    WRITE(*,*)'CPON(I,2)=',CPON(I,2)
        !    WRITE(*,*)'CPON(I,3)=',CPON(I,3)
        !    WRITE(*,*)'CPOC(I,1)=',CPOC(I,1)
        !    WRITE(*,*)'CPOC(I,2)=',CPOC(I,2)
        !    WRITE(*,*)'CPOC(I,3)=',CPOC(I,3)
        !    WRITE(*,*)'CPOS(I)=',CPOS(I)
        !    READ(*,*)
        !ENDIF                
        
   ENDDO
   !IF(MSR)write(*,*)'LBnote b: CPOC_GL(38,2)=',CPOC_GL(38,2), '  CPOC(38,2)=',CPOC(38,2)

!***** Update sediment concentrations

   !calculation deposition feeder and flux terms
   !due to deposition feeder
   IF(DFEEDER)THEN
      CALL DF_CALC
   ENDIF
  
   
!   IDEBUG=6
!   IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
   
   ! B Clark shared variables to pass out to the DOM calc subroutine 
#if defined(SEDIMENT_DOM)
      W2_SHARE = VSED;
	  SED_TEMP_SHARE = CTEMP
	  ZSED_SHARE = HSED
!	  write(*,*) 'sed Depth', ZSED_SHARE ! B Clark Debug
	  
IF(HYDRO_FLAG) THEN  ! Added by B Clark for hydrolysis of POC --> DOC  Nov 2015
   DO I = 1,MLOC 

	 TEMPD = CTEMP(I)

     XKPOC1 = GET_ZHTAPOC1(TEMPD)!*H2 ! keep in 1/d for DOM module
     XKPOC2 = GET_ZHTAPOC2(TEMPD)!*H2
     XKPOC3 = GET_ZHTAPOC3(TEMPD)!*H2
     XKPOP1 = GET_ZHTAPOP1(TEMPD)!*H2
     XKPOP2 = GET_ZHTAPOP2(TEMPD)!*H2
     XKPOP3 = GET_ZHTAPOP3(TEMPD)!*H2
     XKPON1 = GET_ZHTAPON1(TEMPD)!*H2
     XKPON2 = GET_ZHTAPON2(TEMPD)!*H2
     XKPON3 = GET_ZHTAPON3(TEMPD)!*H2
	 
	 FX_POC1TM1(I)= XKPOC1 *POC1TM1S(I)/1000
	 FX_POC2TM1(I)= XKPOC2*POC2TM1S(I)/1000
	 FX_POC3TM1(I)= XKPOC3*POC3TM1S(I)/1000	
	 FX_PON1TM1(I)= XKPON1*PON1TM1S(I)/1000
	 FX_PON2TM1(I)= XKPON2*PON2TM1S(I)/1000
	 FX_PON3TM1(I)= XKPON3*PON3TM1S(I)/1000
	 FX_POP1TM1(I)= XKPOP1*POP1TM1S(I)/1000
	 FX_POP2TM1(I)= XKPOP2*POP2TM1S(I)/1000
	 FX_POP3TM1(I)= XKPOP3*POP3TM1S(I)/1000
	 
   END DO
	 
END IF  

  ! Diffusion rate for application the DOM solver for layer m/d ---> 1/d
   IF(SED_DOM_FLAG) THEN
      
      CALL SED_DOM_CALC   ! B CLARK SED_DOM UPDATE
	   
   ENDIF
#endif   
   
   !IF(MSR)write(*,*)'LBnote i: POC2TM1S(38)=',POC2TM1S(38),'  POC2=',POC2
   DO I=1,MLOC

    !get previous time step concentrations of POM and inorganics
!write(*,*)'POC1TM1S = ',POC1TM1S(I)
     POC1TM1  = POC1TM1S(I)        
     POC2TM1  = POC2TM1S(I)
     POC3TM1  = POC3TM1S(I)
     
     PON1TM1  = PON1TM1S(I)
     PON2TM1  = PON2TM1S(I)
     PON3TM1  = PON3TM1S(I)
     
     POP1TM1  = POP1TM1S(I)
     POP2TM1  = POP2TM1S(I)
     POP3TM1  = POP3TM1S(I)
     
     PSISEDTM1   = PSISEDTM1S(I)
     
     BENSTRTM1  = BENSTRTM1S(I)
     
     SODTM1        = SODTM1S(I)
     
     NH41TM1  = NH41TM1S(I)
     NH4T2TM1 = NH4T2TM1S(I)
     
     NO31TM1  = NO31TM1S(I)
     NO3T2TM1 = NO3T2TM1S(I)
     
     PO41TM1  = PO41TM1S(I)
     PO4T2TM1 = PO4T2TM1S(I)
     
     SI1TM1   = SI1TM1S(I)
     SIT2TM1  = SIT2TM1S(I)
     
     
     HS1TM1   = HS1TM1S(I)
     HST2TM1  = HST2TM1S(I)
     
     
     CH41TM1  = CH41TM1S(I)            ! CH4 of layer 1
     CH4T2TM1 = CH4T2TM1S(I)           ! total CH4 of layer 2
     
     SO4T2TM1 = SO4T2TM1S(I)           ! total SO4 of layer 2
     
     !WLong dfeedm1  = dfeedm1s(I)  !deposition feeder befor integrating to next time step
                                    !moved to mod_df.F
     
    ! ACCOUNT FOR SAV NUTRIENT UPTAKE
    ! NOTE THIS IS OVER ALL CELLS, NOT JUST SAV CELLS SO SEDNH4SAV
    ! MUST BE ZEROED OUT FOR ALL CELLS

     IF (SAV_CALC) THEN
       NH4T2TM1 = NH4T2TM1 - 1000.*SEDNH4SAV(I)*DLTS/HSED(I)
       NO3T2TM1 = NO3T2TM1 - 1000.*SEDNO3SAV(I)*DLTS/HSED(I)  !Wen Long: Added SEDNO3SAV here      
       PO4T2TM1 = PO4T2TM1 - 1000.*SEDPO4SAV(I)*DLTS/HSED(I)
     ENDIF
        
     
        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'POC1TM1=',POC1TM1
        !    WRITE(*,*)'POC2TM1=',POC2TM1
        !    WRITE(*,*)'POC3TM1=',POC3TM1
        !    WRITE(*,*)'PON1TM1=',PON1TM1
        !    WRITE(*,*)'PON2TM1=',PON2TM1
        !    WRITE(*,*)'PON3TM1=',PON3TM1
        !    WRITE(*,*)'POP1TM1=',POP1TM1
        !    WRITE(*,*)'POP2TM1=',POP2TM1
        !    WRITE(*,*)'POP3TM1=',POP3TM1
        !    
        !    WRITE(*,*)'PSISEDTM1=',PSISEDTM1
        !    WRITE(*,*)'BENSTRTM1=',BENSTRTM1
        !    
        !    WRITE(*,*)'NH41TM1=',NH41TM1
        !    WRITE(*,*)'NH4T2TM1=',NH4T2TM1
        !    
        !    WRITE(*,*)'NO31TM1=',NO31TM1
        !    WRITE(*,*)'NO3T2TM1=',NO3T2TM1
        !    
        !    WRITE(*,*)'PO41TM1=',PO41TM1
        !    WRITE(*,*)'PO4T2TM1=',PO4T2TM1
        !    
        !    WRITE(*,*)'SI1TM1=',SI1TM1
        !    WRITE(*,*)'SIT2TM1=',SIT2TM1
        !    
        !    WRITE(*,*)'HS1TM1=',HS1TM1
        !    WRITE(*,*)'HST2TM1=',HST2TM1
        !    
        !    WRITE(*,*)'CH41TM1=',CH41TM1
        !    WRITE(*,*)'CH4T2TM1=',CH4T2TM1
        !    
        !    WRITE(*,*)'SO4T2TM1=',SO4T2TM1
        !    
        !    READ(*,*)
        !ENDIF                     
     
!    IDEBUG=7
!    IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
        
     BFORMAX  = BFORMAXS(I)
     ISWBEN   = ISWBENS(I)
     H2       = HSED(I)

    !******* Sedimentation, mixing rates, and sediment temperature

     W2    = VSED(I)        !m/d
     DPP    = VPMIX(I)        !Particle mixing coefficient m**2/day
     !WLong renamed DD to DDP so that is is not duplicated with DD elsewhere out of sediment module
     !DD    = VDMIX(I)        
     DDP    = VDMIX(I)        !dissolved phase diffusion (mixing) m**2/day
     TEMPD = CTEMP(I)        !temperature of sediments
     STP20 = TEMPD-20.

 
     CSISAT= CSISATT20*THTASISAT**STP20    !WLong calculation saturation concentration based on temperature control

        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'BFORMAX=',BFORMAX
        !    WRITE(*,*)'ISWBEN=',ISWBEN
        !    WRITE(*,*)'H2=',H2
        !    WRITE(*,*)'W2=',W2
        !    WRITE(*,*)'DDP=',DDP
        !    WRITE(*,*)'TEMPD=',TEMPD
        !    WRITE(*,*)'STP20=',STP20
        !    WRITE(*,*)'CSISAT=',CSISAT            
        !    READ(*,*)
        !ENDIF                     
          
!     
! Overlying water concentrations
! Converted to mg/m^3
!     
     
     KWC  = KBM1
     
     !PO40 (mgP/m^3)
     
     !steady state estimate (isotherm) of dissolved PO4 in overlying water column
     DF   = 1./(1.+KADPO4*SSI(I,KWC))
     PO4AVL = DF*PO4(I,KWC)                !avaiable PO4 for algae growth (dissolved)
     PO40 = PO4AVL*1000.                !convert from gP/m^3 to mgP/m^3

     
        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'DF=',DF
        !    WRITE(*,*)'PO4AVL=',PO4AVL
        !    WRITE(*,*)'PO40=',PO40
        !    READ(*,*)
        !ENDIF                
     
     !NH40 (mgN/m^3)
     NH40 = NH4(I,KWC)*1000.            !convert from gN/m^3 to mgN/m^3
     NO30 = NO3(I,KWC)*1000.            !convert from gN/m^3 to mgN/m^3
     
    ! if(I.eq.38) then
    !    write(*,*)'LBnote: NO3 water at node 38 (mgN/m^3):'
    !      write(*,*)'K=',1,'  NO3(38,K)=',NO3(38,1)*1000.
    !      write(*,*)'K=',5,'  NO3(38,K)=',NO3(38,5)*1000.
    !      write(*,*)'K=',10,'  NO3(38,K)=',NO3(38,10)*1000.
    ! endif
     
     !SI0 (mgSi/m^3)
     DF   = 1./(1.+KADSA*SSI(I,KWC))
     SI0  = DF*SIAT(I,KWC)*1000.        !convert from gSi/m^3 to mgSi/m^3
     
         !IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'DF=',DF
        !    WRITE(*,*)'SI0=',SI0/1000.0
        !    READ(*,*)
        !ENDIF                
     
     !O20 (mgO2/L) = (gO2/m^3)
     O20  = AMAX1(DOXG(I,KWC),0.010)    !mgO2/L
     
     !HS0 (gO2/m^3)
     HS0  = COD(I,KWC)                    !mgO2/L (in oxygen equivalents)
     
     !salinity (ppt)
     SAL  = SALT(I,KWC)                    !overlying water salinity (ppt)
     

     !SO40 (gO2/m^3)
     
    ! lets flag DOXG if it starts near zero as it may cause root finder problems
    ! with suspension and deposit feeders on.
    !

    !       Regression function to get SO4 concentration from SAL
    !       [SO4] = 20 mgO2/L        for        [Cl] < 6 mgO2/L
    !             = (10/3)[Cl]       for        [Cl] > 6 mgO2/L
    !       1 ppt = 607.445 mg/L Cl

     IF (SAL > 0.0099) THEN
       SO40MG = 20.0 + (27./190.)*607.445*SAL        !mgSO4/L
     ELSE
       SO40MG = 20.0                                !mgSO4/L
     ENDIF
     !**** units: so4 in o2 equivalents
     !     SO4 (mg so4/L)* 1 mmol SO4 /98 mg SO4 * 2 mmol O2/ 1 mmol SO4
     !     * 32 mg O2 / mmol O2= 0.65306122
     !
     !==> 1 mgSO4/L ~ 0.65306122 mgO2/L   
 
     SO40=SO40MG*0.65306122 !mgO2/L, i.e gO2/m^3

    !CH40 (gO2/m^3)
     IF(ITVWCLF/=1)THEN
        CH40 = 0.
     ELSE
         CH40 = CH40_R1
     ENDIF
    
    !Methane saturation (gO2/m^3) in oxygen equivalents
     
    !WLong   CH4SAT = 99.*(1.+(ZD(I,KWC)+D(I)*DZ(KWC)+HSED(I))/10.)*0.9759**STP20
!     CH4SAT = 99.*(1.+(D(I)+HSED(I))/10.)*0.9759**STP20  
     CH4SAT = 100.*(1.+(D(I)+HSED(I))/10.)*1.024**(-STP20) !WLong, 0.9759 should be 1.024 according to QAPP
                                                           !Also exponent is 20-T rather than T-20
                                                           !And 1/1.024 = 0.9765625
                                                           !     1/0.9759=1.0246951531919255
                                                           !which means they are close!
        IDEBUG_SED=IDEBUG_SED+1

        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'NH40=',NH40
        !    WRITE(*,*)'NO30=',NO30
        !    WRITE(*,*)'DF=',DF
        !    WRITE(*,*)'SI0=',SI0
        !    WRITE(*,*)'O20=',O20
        !    WRITE(*,*)'HS0=',HS0
        !    WRITE(*,*)'SO40MG=',SO40MG
        !    WRITE(*,*)'SO40=',SO40
        !    WRITE(*,*)'CH4SAT=',CH4SAT
        !    WRITE(*,*)'CH40=',CH40
        !    READ(*,*)
        !ENDIF                     
     
!******* Evaluate the temperature dependent coefficients

     ITEMP    = 10.*TEMPD+1

!******* Salinity dependence of nitrification and denitrification

     IF (SAL < SALTND) THEN
       XAPPNH4  = GET_ZHTANH4F(TEMPD)
       XAPP1NO3 = GET_ZHTANO3F(TEMPD)
     ELSE
       XAPPNH4  = GET_ZHTANH4S(TEMPD)
       XAPP1NO3 = GET_ZHTANO3S(TEMPD)
     ENDIF
     
     XAPP1HSD   = GET_ZHTA1HSD(TEMPD)
     XAPP1HSP   = GET_ZHTA1HSP(TEMPD)
     
     XK2NO3   = GET_ZHTAK2NO3(TEMPD)!*H2        !m/d     !WLong (3-36-9)    !WLong, no need to multiply by H
     XKSI     = GET_ZHTASI(TEMPD)*H2            !m/d     !DMD book eqn(7.7)
     
     XAPPCH4  = GET_ZHTACH4(TEMPD)

!    KL12NOM  = DDP/H2*GET_ZL12NOM(TEMPD)         !OLD code uses H2 (this depends on definition of DD)
     
     KL12NOM  = DDP/(H2/2)*GET_ZL12NOM(TEMPD)   !normalized porewater diffusion rate across layer1,2 (m/day)
                                                !We see new DDP is half of old DDP to keep KL12NOM the same
                                                !DMD(13.6)

#if defined(SEDIMENT_DOM)	 
	 
	 KL12_SHARE(I)=KL12NOM !m/d
	 
#endif	 
	 
   !   S_SHARE = KL12NOM   
         !
         !POC1R = 100,000 mgC/m^3 ~= 0.1 mgC/gSediment if M2 =1kg/L 
         !***Note it is based on M2, not M1 ****
         
    !    W12NOM   = DPP/H2*GET_ZW12NOM(TEMPD)*POC1/1.0E5     !Old code uses H2
         
         IF(TF_SSTATE)THEN
         
            W12NOM   = DPP/(H2/2)*GET_ZW12NOM(TEMPD)   !WLong, Greg Pellitier does not have POC1 contribution
                                                       !to W12NOM
         ELSE
         
            W12NOM   = DPP/(H2/2)*GET_ZW12NOM(TEMPD)*POC1TM1/(POC1R*1.0E6*M2) !New code uses H2/2 
                                                  !we see new DPP is half of old DPP to keep  WL12NOM same
			
       
		   
                !Here "*1.0E6*M2 converts  POC1R from mgC/gSediment to mgC/m^3
         
                !normalized particle diffusion rate across layer1,2 !WLong changed POC1 to POC1TM1
                !WLong: POC1R = (1.0E5*M2) (mgC/m^3) if POC1R= 0.1mgC/gSediment 
                !WLong notes eqn(3-26) (DMD book page 278, table 13.1, eqn(13.1))
            !
            !M2 has unit kg/L = 1000kg/m^3 =1000,000 g/m^3 
            !
            !POC1R = 0.1mgC/gSediment *M2 ~ 0.1mgC *1000,000 gSediment/m^3  = 1.0E5 *M2 (mgC/m^3)
            !
         
            !
            !In excel program: w12base = Dp * ThtaDp ^ (Tw - 20) / (H2 / 2) * (POC2(1) / 1000# / m2) / POC1R
            !where H2 is repaced by H2/2
            !POC2(1) is POC in layer 2 for G1 class in gO2/m^3Sediment unit (~ 1 mgO2/Lsediment)
            !m2 is sediment concentration in kg/L
            !POC1R is reference POC2(1) for bioturbation , ~= 0.2667 mgO2/gSediment 
            !POC2(1)/1000/m2 = POC2(1)/(1000*m2) converts to mgO2/gSediment 
            !==> POC2(1)/1000/m2 /POC1R finally has unit of (mgO2/gSediment)/(mgO2/gSediment) ==> unitness
            !
            !WLong Derivation:
            !POC1R in DMD eqn (13.1) is 0.1mgC/g sediment ~ 0.1gC/kg 
            !m2 has unit  kg/L            gC      kg
            !==> POC1R = 0.1gC/kg ~0.1 ----- *m2---- = 0.1*m2 gC/L  ~ 0.1*m2 * 1000gC/m^3 = 100*m2 (gC/m^3) =100000*m2 (mgC/m^3)
            !                             kg      L 
            !            ~ 0.1*10E6*m2 (mgC/m^3)    !i.e. POC1R*1.0E6*M2 used in W12NOM calculation above
         
         
         ENDIF
         
         
            IDEBUG_SED=IDEBUG_SED+1


        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'XAPPNH4=',XAPPNH4
        !    WRITE(*,*)'XAPP1NO3=',XAPP1NO3
        !    WRITE(*,*)'XAPP1HSD=',XAPP1HSD
        !    WRITE(*,*)'XAPP1HSP=',XAPP1HSP
        !    WRITE(*,*)'XK2NO3=',XK2NO3
        !    WRITE(*,*)'XKSI=',XKSI
        !    WRITE(*,*)'XAPPCH4=',XAPPCH4
        !    WRITE(*,*)'KL12NOM=',KL12NOM
        !    WRITE(*,*)'W12NOM=',W12NOM
        !    
        !    READ(*,*)
        !ENDIF             
!---------------     
     !IF (ISWBEN == 0) THEN
     !  IF (TEMPD >= TEMPBEN) THEN    !If time exceeding benthic stress time, then reset benthic pressure to zero
     !    ISWBEN  = 1
     !    BFORMAX = 0.
     !  ENDIF
     !  BFOR = KMO2DP/(KMO2DP+O20)    !oxygen regulation for benthic community's bioturbation
     !ELSE
     !  IF (TEMPD < TEMPBEN) THEN
     !    ISWBEN = 0                    !do not set benthic stress to zero
     !  ENDIF
     !  BFORMAX = AMAX1(KMO2DP/(KMO2DP+O20),BFORMAX)
     !  BFOR    = BFORMAX
     !ENDIF
        
    !Wen Long commented the upper method of resetting BFORMAX to zero
!------------------
        
        
     BFOR= KMO2DP/(KMO2DP+O20)
    
     IF(.NOT.TF_SSTATE)THEN  
     
        !predict benthic pressure at new time setp based on previous value    
        BENSTR = (BENSTRTM1+DLTS*BFOR)/(1.+KBENSTR*DLTS)  !WLong (3-44)
                
     ELSE !steady state
     
        BENSTR = BFOR/KBENSTR                              !WLong (3-42)
        
        BENSTR = 0.d0  !Greg set BENSTR to zero for steady state
        
     ENDIF
     
     BFORMAX=AMAX1(BENSTR,BFORMAX)
     
!## -- add minimum mixing term and bio-irrigation formulation
!##
!##     W12    = W12NOM*(1.-KBENSTR*BENSTR)
!##     KL12   = KL12NOM
!##     w12min = Dpmin/(h2/2) is minimum particle mixing
     
!     W12MIN = DPMIN/H2 !Old Code
     W12MIN = DPMIN/(H2/2.0) !New Code (WLong)
                             !We see that new DPMIN is half of old DPMIN to keep WL12MIN the same
     W12    = W12NOM*(1.-KBENSTR*BFORMAX)!+W12MIN    !WLong removed minimum mixing 
        
                            !Greg used BFORMAX above for W12 calculation

!## -- klbnth is ratio of bio-irrigation to bio-particle mixing
     KL12   = KL12NOM + KLBNTH*W12        !dissolved phase mixing (KL12) affected by bio-particle mixing W12

!******* Lookup reaction rates (m/d)

     ITEMP  = 10.*TEMPD+1

     XKPOC1 = GET_ZHTAPOC1(TEMPD)*H2
     XKPOC2 = GET_ZHTAPOC2(TEMPD)*H2
     XKPOC3 = GET_ZHTAPOC3(TEMPD)*H2
     XKPOP1 = GET_ZHTAPOP1(TEMPD)*H2
     XKPOP2 = GET_ZHTAPOP2(TEMPD)*H2
     XKPOP3 = GET_ZHTAPOP3(TEMPD)*H2
     XKPON1 = GET_ZHTAPON1(TEMPD)*H2
     XKPON2 = GET_ZHTAPON2(TEMPD)*H2
     XKPON3 = GET_ZHTAPON3(TEMPD)*H2
     
     !IF(MSR)WRITE(*,*)'LBnote: GET_ZHTAPOC2(TEMPD)=',GET_ZHTAPOC2(TEMPD),'  H2 (m)=',H2,' TEMPD=',TEMPD
     
        IDEBUG_SED=IDEBUG_SED+1

        IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'BFOR=',BFOR
        !    WRITE(*,*)'BENSTR=',BENSTR
        !    WRITE(*,*)'W12MIN=',W12MIN
        !    WRITE(*,*)'W12=',W12
            !WRITE(*,*)'KL12=',KL12
        !    WRITE(*,*)'ITEMP=',ITEMP
        !    WRITE(*,*)'XKPOC1=',XKPOC1
        !    
        !    WRITE(*,*)'XKPOC2=',XKPOC2
        !    WRITE(*,*)'XKPOC3=',XKPOC3
        !    WRITE(*,*)'XKPOP1=',XKPOP1
        !    WRITE(*,*)'XKPOP2=',XKPOP2
        !    WRITE(*,*)'XKPOP3=',XKPOP3
        !    WRITE(*,*)'XKPON1=',XKPON1
        !    WRITE(*,*)'XKPON2=',XKPON2
        !    WRITE(*,*)'XKPON3=',XKPON3
        !    
            !READ(*,*)
        ENDIF                  
!******* Calculate sediment concentrations

     DLTS_H2=DLTS/H2   !multipler dt/H2 in front of sediment diagensis equations (see WLong notes eqn(2-12))

     frpon1=1.-(frpon(I,2)+frpon(I,3))
     frpop1=1.-(frpop(I,2)+frpop(I,3))
     frpoc1=1.-(frpoc(I,2)+frpoc(I,3))

!     IDEBUG=8
!     IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
   
    !sediment silicate
     FD2  = 1./(1.+M2*PIE2SI)                        !WLong (4-6-7)(3-7)
    K3   = XKSI*(CSISAT-FD2*SIT2TM1)/(PSISEDTM1+KMPSI)    !Dissolution rate (m/day) (WLong eqn (2-16))

	!write(*,*)'TF_SSTATE = ',TF_SSTATE
	
     IF(.NOT.TF_SSTATE)THEN  !non-steady state calculation, integration one time step DLTS
     
     
            !non-steady state solution of sediment diagenesis (mgN/m^3, mgC/m^3, mgP/m^3)
            !based on assuming H1 <<H2
            !See WLong notes eqn (2-13-1) for over all solution of POC, PON, POP
            !
            !WLong notes (2-15) for POCi
	 
			! B Clark we take the loss term that is being passed to diagenesis calculations, XKPOC1* POC1 and transfer to colored/non-colored portions of DOC in mod_sed_doc.F
     !write(*,*)'POC1TM1 = ',POC1TM1    
	 POC1 = (JPOC(I,1)*DLTS_H2+POC1TM1)                       & 
                    /(1.+(XKPOC1+W2)*DLTS_H2)
					
            POC2 = (JPOC(I,2)*DLTS_H2+POC2TM1)                       &
                    /(1.+(XKPOC2+W2)*DLTS_H2)
                    
            POC3 = (JPOC(I,3)*DLTS_H2+POC3TM1)                          &
                    /(1.+(XKPOC3+W2)*DLTS_H2)
			
					
        
            IDEBUG_SED=IDEBUG_SED+1
                !IF(MSR)THEN        
                !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
                !    WRITE(*,*)'POC1=',POC1
                !    WRITE(*,*)'POC2=',POC2
                !    WRITE(*,*)'JPOC(I,2)=',JPOC(I,2)
                !    WRITE(*,*)'POC2TM1=',POC2TM1
                !    WRITE(*,*)'XKPOC2=',XKPOC2
                !    WRITE(*,*)'DLTS_H2=',DLTS_H2
                !    WRITE(*,*)'POC3=',POC3
                !    READ(*,*)
                !ENDIF                
            
            !WLong notes (2-18) for PONi
            PON1 =  (JPON(I,1)*DLTS_H2+PON1TM1)                      &
                    /(1.+(XKPON1+W2)*DLTS_H2)

            PON2 = (JPON(I,2)*DLTS_H2+PON2TM1)                       &
                    /(1.+(XKPON2+W2)*DLTS_H2)
        
            PON3 = (JPON(I,3)*DLTS_H2+PON3TM1)                          &
                    /(1.+(XKPON3+W2)*DLTS_H2)
     
#if defined(PHOSPHORUS)    
            !WLong notes (2-18-2) for POPi
            POP1 = (JPOP(I,1)*DLTS_H2+POP1TM1)                       &
                    /(1.+(XKPOP1+W2)*DLTS_H2)

            POP2 = (JPOP(I,2)*DLTS_H2+POP2TM1)                          &
                    /(1.+(XKPOP2+W2)*DLTS_H2)

            POP3 = (JPOP(I,3)*DLTS_H2+POP3TM1)                          &
                    /(1.+(XKPOP3+W2)*DLTS_H2)
#endif
            !## -- modification for detrital Si input to sediment
            !##  PSISED  = (JPOS(I)*DLTS/H2+PSISEDTM1)/(1.+(K3+W2)*DLTS/H2)
     
            PSISED  = ((JPOS(I)+JSIDETR)*DLTS_H2+PSISEDTM1)/(1.+(K3+W2)*DLTS_H2)    !WLong notes (2-16) (DMD book 13.31, 13.33)
            !Maybe we should put JSIDTER in JPOS directly based on Si/C ration in RPOC and LPOC
            
            !Bellow are Greg Pelletier's calculation: (same as James L Martin)
     
            PSISED  = (JPOS(I)*DLTS_H2+PSISEDTM1)/(1.+(XKSI+W2)*DLTS_H2)            !
     

        IDEBUG_SED=IDEBUG_SED+1

        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'POC1=',POC1
        !    WRITE(*,*)'POC2=',POC2
        !    WRITE(*,*)'POC3=',POC3
        !    WRITE(*,*)'PON1=',PON1
        !    WRITE(*,*)'PON2=',PON2
        !    WRITE(*,*)'PON3=',PON3
        !    WRITE(*,*)'POP1=',POP1
        !    WRITE(*,*)'POP2=',POP2
        !    WRITE(*,*)'POP3=',POP3
        !    WRITE(*,*)'FD2 silicate=',FD2
        !    WRITE(*,*)'PSISED=',PSISED
        !    WRITE(*,*)'K3=',K3
        !    READ(*,*)
        !ENDIF             
     
        !WLong JSIDETR here should be included in JPOS calculation
        !see DMD p391 JdetrSi (16A.21)
     
        
        !Modification to PON, POC, POP due to deposition feeder fluxes
        IF(DFEEDER)THEN
     
            POC1=POC1+(                                            &
                -DF_GROW_POC1(I)                                  &   ! deposit feeders
                                            *dlts/h2            &
                +frpoc1*(DF_MORT(I)+DF_PRED(I))                 &    !WLong frpoc1 is questionable here
                                            *dlts/h2               &     !should any mortality and predation loss
            )/(1.+(XKPOC1+W2)*DLTS_H2)                                !contribute to POC1 at all?
     
            POC2=POC2+(                                            &
                -DF_GROW_POC2(I)                                &    ! deposit feeders
                                       *dlts/h2                    &
                +frpoc(I,2)*(DF_MORT(I)+DF_PRED(I))              &
                                       *dlts/h2                 &
                )/(1.+(XKPOC2+W2)*DLTS_H2)
     
            PON1=PON1+(                                            &
                -DF_GROW_POC1(I)                                &    ! deposit feeders
                                       *dlts/h2/amcn            &
                +frpon1*(DF_MORT(I)+DF_PRED(I))                 &     !WLong frpoN1 is questionable here
                                       *dlts/h2/amcn            &     !should any mortality and predation loss
                )/(1.+(XKPON1+W2)*DLTS_H2)                             !contribute to PON1 at all? Yes
        
            PON2=PON2+(                                            &
                -DF_GROW_POC2(I)                                &    ! deposit feeders
                                       *dlts/h2/amcn            &
                +frpon(I,2)*(DF_MORT(I)+DF_PRED(I))              &
                                        *dlts/h2/amcn            &
                )/(1.+(XKPON2+W2)*DLTS_H2)

            POP1=POP1+(                                            &
                -DF_GROW_POC1(I)                                &    ! deposit feeders
                                       *dlts/h2/amcp            &
                +frpop1*(DF_MORT(I)+DF_PRED(I))                 &     !WLong frpoN1 is questionable here
                                        *dlts/h2/amcp            &    !should any mortality and predation loss
                )/(1.+(XKPOP1+W2)*DLTS_H2)                             !contribute to POP1 at all? Yes

            POP2 = POP2+(                                        &
                -DF_GROW_POC2(I)                                &    ! deposit feeders
                                       *dlts/h2/amcp            &
                +frpop(I,2)*(DF_MORT(I)+DF_PRED(I))              &
                                        *dlts/h2/amcp            &
                )/(1.+(XKPOP2+W2)*DLTS_H2)

        ENDIF

     ELSE    !STEADY STATE CALCULATION
       !IF(MSR)write(*,*)'LBnote ii: XKPOC2=',XKPOC2,'  W2=',W2
            POC1 = JPOC(I,1)/(XKPOC1+W2)    !WLong (2-23-1)
            POC2 = JPOC(I,2)/(XKPOC2+W2)
            POC3 = JPOC(I,3)/(XKPOC3+W2)
       !IF(MSR)write(*,*)'LBnote iii: JPOC(I,2)=',JPOC(I,2),'  POC2=',POC2
            PON1 = JPON(I,1)/(XKPON1+W2)    !WLong (2-23-1)
            PON2 = JPON(I,2)/(XKPON2+W2)
            PON3 = JPON(I,3)/(XKPON3+W2)
#if defined(PHOSPHORUS)     
            POP1 = JPOP(I,1)/(XKPOP1+W2)    !WLong (2-23-1)
            POP2 = JPOP(I,2)/(XKPOP2+W2)
            POP3 = JPOP(I,3)/(XKPOP3+W2)
#endif     
            !## -- modification for detrital Si input to sediment 
            PSISED  = (JPOS(I)+JSIDETR)/(K3+W2)    !WLong notes (2-16) (DMD book 13.31, 13.33)
            !Bellow are Greg Pelletier's calculation: (same as James L Martin)     
            PSISED  = JPOS(I)/(XKSI+W2)        !WLong (2-23-1)
            
            IF(DFEEDER)THEN
    !write(*,*)'DFEEDER = ',DFEEDER
                POC1 = (JPOC(I,1)-DF_GROW_POC1(I)                            &!
                                 +frpoc1*(DF_MORT(I)+DF_PRED(I))            &!
                        )/(XKPOC1+W2)    !WLong (2-23-2)
                
                
                POC2 = (JPOC(I,2)-DF_GROW_POC2(I)                            &!
                                 +frpoc(I,2)*(DF_MORT(I)+DF_PRED(I))        &!
                        )/(XKPOC2+W2)    !WLong (2-23-2)

                PON1 = (JPON(I,1)-DF_GROW_POC1(I)/amcn                         &!
                                 +frpon1*(DF_MORT(I)+DF_PRED(I))/amcn         &!
                        )/(XKPON1+W2)    !WLong (2-23-2)
                
                
                PON2 = (JPON(I,2)-DF_GROW_POC2(I)/amcn                         &!
                                 +frpon(I,2)*(DF_MORT(I)+DF_PRED(I))/amcn     &!
                        )/(XKPON2+W2)    !WLong (2-23-2)

#if defined(PHOSPHORUS)
                POP1 = (JPOP(I,1)-DF_GROW_POC1(I)/amcp                        &!
                                 +frpop1*(DF_MORT(I)+DF_PRED(I))/amcp        &!
                        )/(XKPOP1+W2)    !WLong (2-23-2)
                
                
                POP2 = (JPOP(I,2)-DF_GROW_POC2(I)/amcp                         &!
                                 +frpop(I,2)*(DF_MORT(I)+DF_PRED(I))/amcp    &!
                        )/(XKPOP2+W2)    !WLong (2-23-2)
#endif
                !WLong: if have W2 and POM G3 class as no relaction, how could it reach steady state?
                !yes, it can. Should we also calculate POC3 and POP3? Yes
                
                !WLong: should we also have DFEEDER contribution to silicate in sediments?
                !    -maybe
                
            ENDIF

     ENDIF
     
     
     !
     !Make sure PON, POC, POP greater than zero (this should not happen if
     !an implicit scheme is used and the reactions are of first order to PON, POC, POP
     !for the deposition feeder model
     !
     
     IF(PON1 < 0.0)PON1=0.0
     IF(PON2 < 0.0)PON2=0.0
     IF(POC1 < 0.0)POC1=0.0 
     IF(POC2 < 0.0)POC2=0.0
     IF(POP1 < 0.0)POP1=0.0
     IF(POP2 < 0.0)POP2=0.0
     IF(PSISED < 0.0)PSISED=0.0  !WLong make sure PSISED >=0
     
!******* Assign diagenesis values for sediment model
!basically the decomposition of PON, POC, POP in sediments
!aided by respiration of deposition feeder
     
     !JCX = XKPOC1*POC1+XKPOC2*POC2+XKPOC3*POC3    !WLong notes (2-7)               (m/d * mgC/m^3  ==> mgC/m^2/day)
     !JNX = XKPON1*PON1+XKPON2*PON2+XKPON3*PON3    !WLong notes (2-8)'s JN term (mgN/m^2/day)
     !JPX = XKPOP1*POP1+XKPOP2*POP2+XKPOP3*POP3  !WLong notes (2-9)'s JP term (mgP/m^2/day) 


#if defined(SEDIMENT_DOM)	 
IF(HYDRO_FLAG .AND. SED_DOM) THEN

    ! write(*,*)'JDOC = ',JDOCX_SHARE(I) 
	  JCX=JDOCX_SHARE2(I)*1000    ! B Clark convert from g/m^2/day to mg/m^2/day
	  JNX=JDONX_SHARE2(I)*1000
	  JPX=JDOPX_SHARE2(I)*1000
						  
	! Layer 1 input of DON remineralization to Layer 1 NH4 					  
      JNX_LAYER1 =  JDONX_SHARE1(I)*HSED1(I)*1000  ! g N m^-3 d * depth of layer 1 (m) * 1000 mg g^-1
		
!	 write(*,*)'JCX before layer 1',JCX									 										 ! ----> mg N m^-2 d^-1     
	  JCX = JCX+JDOCX_SHARE1(I)*HSED1(I)*1000  !JCX only goes to SOD
	  JPX_LAYER1 =  JDOPX_SHARE1(I)*HSED1(I)*1000  ! g P m^-3 d * depth of layer 1 (m) * 1000 mg g^-1

ELSE
#endif
     JCX = XKPOC1*POC1+XKPOC2*POC2+XKPOC3*POC3    !WLong notes (2-7)               (m/d * mgC/m^3  ==> mgC/m^2/day)
     JNX = XKPON1*PON1+XKPON2*PON2+XKPON3*PON3    !WLong notes (2-8)'s JN term (mgN/m^2/day)
	 !write(*,*)'JNX = ',JNX
     JPX = XKPOP1*POP1+XKPOP2*POP2+XKPOP3*POP3  !WLong notes (2-9)'s JP term (mgP/m^2/day) 
	!write(*,*)'JCX in 1 =',JCX 
#if defined(SEDIMENT_DOM)	 
	IF(SED_DOM)THEN  ! B Clark Sediment DOM addition
	  
	  JCX=JCX+JDOCX_SHARE2(I)*1000    ! B Clark convert from g/m^2/day to mg/m^2/day
	  JNX=JNX+JDONX_SHARE2(I)*1000
	  JPX=JPX+JDOPX_SHARE2(I)*1000
!	write(*,*)'JCX in 2  = ',JCX  
						  
	! Layer 1 input of DON remineralization to Layer 1 NH4 					  
      JNX_LAYER1 =  JDONX_SHARE1(I)*HSED1(I)*1000  ! g N m^-3 d * depth of layer 1 (m) * 1000 mg g^-1
	  JPX_LAYER1 =  JDOPX_SHARE1(I)*HSED1(I)*1000  ! g P m^-3 d * depth of layer 1 (m) * 1000 mg g^-1										 ! ----> mg N m^-2 d^-1     
  	  JCX = JCX+JDOCX_SHARE1(I)*HSED1(I)*1000
	  
	ENDIF
	
ENDIF

	IF(WETLAND_FLAG) THEN
		!	write(*,*)'JCX in 2  = ',JCX  
	  	 !***** Take temperature integration step
		!  B CLark get root mass normalized NH4 uptake rate
		!	  write(*,*) 'Vmax marsh before = ', VMAX_MARSH_NH4

		VMAX_MARSH_NH4_out =   MARSH_MASS(I)*MARSH_FRAC(I) * (VMAX_MARSH_NH4) * 1.08**(CTEMP(I)-20)  

		!! add in scaling factor for the sediment layer 1 depth so 
	!	write(*,*) 'Vmax marsh out = ', VMAX_MARSH_NH4_out
		!write(*,*)' Marsh Mass = ',MARSH_MASS(I)

		! Add in our Michaelis Menton marsh plant uptake rate to the layer two flux
		!From Bradley and Morris 1990 marsh plant uptake in ~20 PSU salinity
		 ! Flux into layer two from remin (mg N/m^2/d) - wetland plant uptake (mg N /m^2/d) * layer two thickness (m)  = mg N /m^2 /d
!		write(*,*)' JNX before = ',JNX
		
                Plant_NH4_UPTAKE(I) = VMAX_MARSH_NH4_out * NH42TM1S(I)/(KMNH4_marsh + NH42TM1S(I))
		
	!	IF(Plant_NH4_UPTAKE(I) > 0.0) THEN
	!	write(*,*)'Plant_NH4_UPTAKE = ',Plant_NH4_UPTAKE(I)
		
	!		write(*,*)',JNX before',JNX
		
			JNX = JNX - PLANT_NH4_UPTAKE(I)
		       ! JNX = -10.0
	!		write(*,*)'JNX after =',JNX
			
	!	ENDIF
   
	!	write(*,*)' JNX after= ',JNX
     ENDIF



#endif

        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'JCX=',JCX
        !    WRITE(*,*)'JNX=',JNX
        !    WRITE(*,*)'JPX=',JPX
        !    READ(*,*)
        !ENDIF             

    IF(DFEEDER)THEN
        JCX=JCX+DF_RESP(I)                                    !deposition feeder contribution
        JNX=JNX+DF_RESP(I)*(1./AMCN)                        !deposition feeder contribution converted from C to N
        JPX=JPX+DF_RESP(I)*(1./AMCP)                        !deposition feeder contribution converted from C to P
    ENDIF
    
    !WLong: The sediment flux are solved in the following procedures 
    ! (1) take initial guess of SOD and H1
    ! (2) solve WLong notes (3-21) (3-22) for CT1 and CT2 at t+dt time for NH4,NO3, HS and CH4
    !       (a) solve for NH4 at t+dt in layer 1 and layer 2
    !       (b) compute oxygen demand consumed by nitrification, called NSOD
    !       (c) solve for NO3 at t+dt in layer 1 and layer 2
    !       (d) compute CH4 (freshwater) or HS( salt water) oxidation
    !             <i> for salt water, compute HS reaction terms and compute SOD due to HS oxidation
    !             <ii> for fresh water, compute CH4 flux by establishing the chemical specific conditions
    !                    [1] compare computed and saturation concentrations of CH4 and correct it
    !                     [2] calculate CSOD due to methane oxidation
    !             <1>compute CSOD  due to HS or CH4
    !             <2>compute flux of JCH4 or JHS
    !             <3>compute total SOD = NSOD + CSOD
    !       (e) refine estimate of SOD, make new estimate of SOD
    !  (3)go to (2) if no covergence of SOD and H1
    !  (4)solve Si and PO4 using steady state equations
    
!
!****** Evaluate the NH4, NO3, SO4, HS, CH4 and SOD equations *****
!

    !initial guess of SOD (including SOD due to deposition feeder)
     IF(DFEEDER)THEN
        DF_SOD_TMP=DF_SOD(I)
     ELSE
        DF_SOD_TMP=0.0
     ENDIF
     
     IF(SAV_CALC)THEN
        SAV_SOD_TMP=-SEDDOSAV(I) !SAV_SOD(I)  !SEDDOSAV
     ELSE
        SAV_SOD_TMP=0.0
     ENDIF
     
!       IDEBUG=10
!     IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

    ! REAL(SP) ::     DFSOD_SED !SOD due to deposition feeder    
    !REAL(SP) ::     SAVSOD_SED !SOD due to SAV root respiration
    
     !
     !solve NH4, NO3, SO4, HS, CH4 using ZBRENT function which
     !calls SEDF, where SEDF solves NH4, NO3, SO4, HS, CH4
     !and finally find out what SOD is until ZBRENT converges
     !
     IF(TF_SSTATE)THEN
        SODTM1=  2.667*1.0E-3*JCX  &!mgC/m^2/day converted to mgO2/m^2/day then to gO2/m^2/day
                +1.714*1.0E-3*JNX   !mgN/m^2/day converted to mgO2/m^2/day then to gO2/m^2/day
     ELSE
        !SODTM1=SODTM1  !no change
        
     ENDIF
     
     SODMIN = 1.E-15  
     SODMAX = 1E15   
    
  
     SOD = ZBRENT(IERR,DF_SOD_TMP,SAV_SOD_TMP,TF_SSTATE, SODTM1) !returns IERR and SOD
     
     !if(MSR)write(*,*)'LBnote: I=',I,'  NO30=',NO30
     
     !ZBRENT returns IERR and value ZBRENT as SOD
     !if IERR = 0   good, find soultion and SOD converges
     !   IERR = 1    did not find solution of SOD within [SODMIN, SODMAX] range
     !   IERR = 2    did not find the solution within the IMAX interations
     

        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'SOD=',SOD
        !    WRITE(*,*)'DF_SOD_TMP=',DF_SOD_TMP
        !    WRITE(*,*)'IERR=',IERR
        !    READ(*,*)
        !ENDIF             
     
!      IDEBUG=11
!     IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

!          
! Wen Long --we should add oxygenb consumption due to ROOT and Tuber respiration which contributed
!            to CO2 generation by losing ROOT and Tuber mass with a fraction of FDOSAV 
!
       !
       !Wen Long: We should have SEDDOSAV included here so that it afffects SOD.
       !          Or this is tricky as DOXG in sediment top layer is assumed to be
       !          a linear distribution. Basically SEDDOSAV should affect
       !          the oxygen flux (SOD) and hence the thickness of the top layer
       !
       !The presense of SAV roots in sediments should change the oxygen content
       ! in sediment. Strickly speaking roots are distrbuted over both aerobic and anaerobic layers
       ! of sediments and two layered structure of sediment chemistry analysis based on DOXG
       ! is not completely valid instantaneously. Over a longer period, the two-layer assumption is still okay
       ! 
       !Without going into detailed effects of SEDDOSAV (root and tuber respiration) on DOXG distribution
       !in sediments, we can simply aggregate vertically and make it affect SOD directly
       !SOD finally act as a sink term to water column, so if we add the following section here
       !we will NOT need the SEDDOSAV terms in DTDO calculation in wqm_kin.F's OXYGEN() subroutines
!
!WL IF(SAV_CALC)THEN          
!WL      SOD=SOD+SEDDOSAV(I)    !need to make sure SEDDOSAV is greater than zero for respiration
!WL ENDIF
!
!WL Alternatively, we can add SEDDOSAV to BENDO calcuation in wqm_sed.F, which finally gets in DTDO in wqm_kin.F
!WL and in that case we are not considering SEDDOSAV as part of SOD which is logically questionable
!WL Ideally we would want to have the leaf production (production of DOXG by photosynethesis) of oxygen NOT counted in SOD,
!WL instead count that as water column DDOXGO source, and have the tuber and root part respiration 
!WL consumption of DOXG counted in SOD
!
    
!       IF (IERR.NE.0.AND.BENTHIC_OUTPUT) WRITE(BFO,3000) IERR,I

        !print error and diagnostic information if estimate of SOD does not work
     
     IF(IERR ==1 .OR. IERR == 2) THEN
     
        IF(DFEEDER)THEN
            DF_SOD_TMP = DF_SOD(I)
            DFEED_TMP = DFEEDM1S(I)
        ELSE
            DF_SOD_TMP = 0.0
            DFEED_TMP = 0.0
        ENDIF
        
        IF(SAV_CALC)THEN
            SAV_SOD_TMP=SEDDOSAV(I) !SAV_SOD(I)
        ELSE
            SAV_SOD_TMP= 0.0
        ENDIF
        
		!this is to use SODMIN and SODMAX to estimate the bounds of error
        SODMIN = 1.E-15
        SODMAX = 1.E15
        
        IF(SFEEDER)THEN !three species of suspension feeder and their contribution to sediment diagenesis
        !    SFEED_TMP_1= SFEED(I,1)    !this needs to be treated as SFEED_SED_SF in mod_sed_sf_exchange_vars.F
        !    SFEED_TMP_2= SFEED(I,2)
        !    SFEED_TMP_3= SFEED(I,3)
            SFEED_TMP_1= 0.    !disabled for now
            SFEED_TMP_1= 0.
            SFEED_TMP_1= 0.
        ELSE
            SFEED_TMP_1= 0.
            SFEED_TMP_2= 0.
            SFEED_TMP_3= 0.
        ENDIF
        
        IF(MSR)WRITE(6,9000)   JDAY,IERR,I,SAL,SO40MG,DFEED_TMP             &
                          ,SFEED_TMP_1,SFEED_TMP_2,SFEED_TMP_3                &!print 3 types of sfeeders
                          ,SODMIN,SODMAX
                          
        IF(MSR)WRITE(6,9911) CSODHS, CSODCH4, CSOD
9911    FORMAT(/1X,' CSODHS, CSODCH4, CSOD'/3E10.3)
        IF(MSR)WRITE(6,9910)  CH41,CH42,HST1,HS1,HS2
9910    FORMAT(/1X,' CH41   CH42   HST1   HS1   HS2'/5E10.3)

        IF(IERR == 2)  THEN     !SOD=ZBRENT(SOD) did not converge 
        
         IF(MSR)WRITE(6,9900)  JDAY,CTEMP(I),POP1,POP2,POP3
         IF(MSR)WRITE(6,9901)  PON1,PON2,PON3,POC1,POC2,POC3
         IF(MSR)WRITE(6,9902)  PO4T2,HST2,SIT2,PSISED
         IF(MSR)WRITE(6,9903)                                    &
              (JPOP(I,1)+JPOP(I,2)+JPOP(I,3))                    &
             ,(JPON(I,1)+JPON(I,2)+JPON(I,3))                    &
             ,(JPOC(I,1)+JPOC(I,2)+JPOC(I,3))
         IF(MSR)WRITE(6,9904)  O20,CSOD,DF_SOD_TMP,SOD,S          &
                       ,H2,HSED(I),VSED(I)
         IF(MSR)WRITE(6,9905)  JPX,JNX,JCX,JO2NH4,XJC1
         IF(MSR)WRITE(6,9906)  JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BENSTR
         IF(MSR)WRITE(6,9907)  PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4T2
         IF(MSR)WRITE(6,9908)  NO30,NO31,NO32,NO3T2,HS1,HS2,HST2
         IF(MSR)WRITE(6,9909)  SI0,SI1,SI2,SIT2
#if defined(MULTIPROCESSOR)
	
         WRITE(*,*)'Oops, zbrent did not converge for box number I=', NGID(I) !WLong we should report global box NGID(I)!for parallel runs
		 
#else  
		  WRITE(*,*)'Oops, zbrent did not converge for box number I=', I
         STOP  
#endif   
		 
		 
       ELSE !(IERR==1) !did not find solution within SODMIN and SODMAX
       
          !If IERR /=0 and /=2 then redo the calculation with minimum SOD and find out the size of error
          !                    also redo the calculation with maximum SOD and find out the size of error
       
         IF(DFEEDER)THEN
            DF_SOD_TMP=DF_SOD(I)
         ELSE
            DF_SOD_TMP=0
         ENDIF
         
         IF(SAV_CALC)THEN
            SAV_SOD_TMP=SEDDOSAV(I) !SAV_SOD(I)
         ELSE
            SAV_SOD_TMP= 0.0
         ENDIF
         !calculate NO3, NH4, SO4, HS etc using SEDF function with a)minimum SOD, b)maximum SOD
         
         !
         !calculate ERROR using minimum SOD, basically report the lower bound of error
         !
         
         ERROR=SEDF(SODMIN,DF_SOD_TMP,SAV_SOD_TMP,TF_SSTATE) 
         
         IF(MSR)WRITE(6,9889)  JDAY,SODMIN,ERROR
         IF(MSR)WRITE(6,9900)  JDAY,CTEMP(I),POP1,POP2,POP3
         IF(MSR)WRITE(6,9901)  PON1,PON2,PON3,POC1,POC2,POC3
         IF(MSR)WRITE(6,9902)  PO4T2,HST2,SIT2,PSISED
         IF(MSR)WRITE(6,9903)  (JPOP(I,1)+JPOP(I,2)+JPOP(I,3))              &
                            ,(JPON(I,1)+JPON(I,2)+JPON(I,3))                &
                            ,(JPOC(I,1)+JPOC(I,2)+JPOC(I,3))
             
         IF(MSR)WRITE(6,9904)  O20,CSOD,DF_SOD_TMP,SOD,S                       &
                                ,H2,HSED(I),VSED(I)
                       
         IF(MSR)WRITE(6,9911)  CSODHS, CSODCH4, CSOD
         IF(MSR)WRITE(6,9905)  JPX,JNX,JCX,JO2NH4,XJC1
         IF(MSR)WRITE(6,9906)  JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BENSTR
         IF(MSR)WRITE(6,9907)  PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4T2
         IF(MSR)WRITE(6,9908)  NO30,NO31,NO32,NO3T2,HS1,HS2,HST2
         IF(MSR)WRITE(6,9909)  SI0,SI1,SI2,SIT2
         
         
         !
         !calculate ERROR using maximum SOD, basically report the upper bound of error
         !
         IF(DFEEDER)THEN
            DF_SOD_TMP=DF_SOD(I)
         ELSE
            DF_SOD_TMP=0
         ENDIF
         
          IF(SAV_CALC)THEN
            SAV_SOD_TMP=SEDDOSAV(I) !SAV_SOD(I)
         ELSE
            SAV_SOD_TMP= 0.0
         ENDIF
         
         !calculate ERROR using maximum SOD
         ERROR=SEDF(SODMAX,DF_SOD_TMP,SAV_SOD_TMP,TF_SSTATE)
         
         IF(MSR)WRITE(6,9889)  JDAY,SODMAX,ERROR
         IF(MSR)WRITE(6,9900)  JDAY,CTEMP(I),POP1,POP2,POP3
         IF(MSR)WRITE(6,9901)  PON1,PON2,PON3,POC1,POC2,POC3
         IF(MSR)WRITE(6,9902)  PO4T2,HST2,SIT2,PSISED
         
         IF(MSR)WRITE(6,9903)  (JPOP(I,1)+JPOP(I,2)+JPOP(I,3))                 &
                              ,(JPON(I,1)+JPON(I,2)+JPON(I,3))                 & 
                              ,(JPOC(I,1)+JPOC(I,2)+JPOC(I,3))
                              
         IF(MSR)WRITE(6,9904)  O20,CSOD,DF_SOD_TMP,SOD,S                        &
                              ,H2,HSED(I),VSED(I)
         IF(MSR)WRITE(6,9911)  CSODHS, CSODCH4, CSOD
         IF(MSR)WRITE(6,9905)  JPX,JNX,JCX,JO2NH4,XJC1
         IF(MSR)WRITE(6,9906)  JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BENSTR
         IF(MSR)WRITE(6,9907)  PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4T2
         IF(MSR)WRITE(6,9908)  NO30,NO31,NO32,NO3T2,HS1,HS2,HST2
         IF(MSR)WRITE(6,9909)  SI0,SI1,SI2,SIT2
         
         !stop the calculation

		 
#if defined(MULTIPROCESSOR)        
         WRITE(*,*) 'Oops, did not find SOD soultion with SODMIN<SOD<SODMAX, SODMIN=', &
		 
                SODMIN, 'SODMAX=',SODMAX,'for box I=',NGID(I) !WLong we should report global box NGID(I)
#else                                                     
				
        WRITE(*,*) 'Oops, did not find SOD soultion with SODMIN<SOD<SODMAX, SODMIN=', &
		 
                SODMIN, 'SODMAX=',SODMAX,'for box I=',I
#endif				
         STOP
         
       ENDIF
	   
	   

9889 FORMAT(/5X,'ZBRENT DIAGNOSTICS AT TIME =',F8.3,               &
            ' FOR SOD =',F8.4,' ERROR =',E12.3/)
9000 FORMAT(/                                                      &
            5X,'ZBRENT FAILURE AT TIME =',F8.3,' WITH IERR=',I2/   &
            5X,'IN SEDIMENT SEGMENT IR=',I5/                       &
            5X,'WITH SALT, SO40MG     =',2E10.3/                   &
            5X,'DFEED=',F10.3,' SFEED=',3F11.3/                    &
            5X,'(SODMIN,SODMAX=',F6.4,F6.1,')'/                    &
            5X,'PROGRAM TERMINATION FOLLOWS DIAGNOSTIC DUMPS')
9900 FORMAT(/1X,' TIME,CTEMP,POP1,POP2,POP3'/8E10.3)
9901 FORMAT(/1X,' PON1,PON2,PON3,POC1,POC2,POC3'/ 8E10.3)
9902 FORMAT(/1X,' PO4T2,HST2,SIT2,PSISED'/8E10.3)
9903 FORMAT(/1X,' JPOP,JPON,JPOC'/8E10.3)
9904 FORMAT(/1X,' O20,CSOD,DFSOD,SOD,S,H2',',HSED,VSED'/10E10.3)
9905 FORMAT(/1X,' JP,JN,JC,JO2NH4,XJC1'/8E10.3)
9906 FORMAT(/1X,' JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BENSTR'/ 8E10.3)
9907 FORMAT(/1X,' PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4T2'/8E10.3)
9908 FORMAT(/1X,' NO30,NO31,NO32,NO3T2,HS1,HS2,HST2'/8E10.3)
9909 FORMAT(/1X,' SI0,SI1,SI2,SIT2'/8E10.3)


     ENDIF !Done with error reporting for node I
	 
	   !  


       !BENCOD(I)=CSOD ! Added by B Clark July 2018, to get the benthic contribution of  H2S and CH4 to COD for the water column
	 
       !IF(NGID(I)==352)THEN ! added by B Clark July 2018

       !    write(*,*)'BENCOD in mod_sed = ',BENCOD(I)
       !ENDIF


!******* accumulate temperature for steady state calculation of POC, PON, POP for G3

      !IDEBUG=12
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

     IF (STEADY_STATE_SED_G3) THEN
       ASDTMP(I) = ASDTMP(I)+TEMPD*DLTS
     ENDIF
     
!Calculate Silica and PO4 for they do not affect SOD

!calculate Si in sediments

     K0H1D = 0.
     K0H1P = 0.
     KMC1  = 0.
     K1H1D = S
     K1H1P = 0.
     K2H2D = 0.
     K2H2P = 0.
     J1    = S*SI0

    ! Oxygen dependency of pie1 (L/kg) (DMD book page 156, 391, 571)

     IF (O20 < O2CRITSI) THEN
       PIE1 = PIE2SI*DPIE1SI**(O20/O2CRITSI)        !WLong notes (4-6-6), DMD  (7.19)
     ELSE    
       PIE1 = PIE2SI*DPIE1SI                        !WLong notes (4-6-5 ), DMD (7.18)
     ENDIF
     
     PIE2 = PIE2SI

        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'K0H1D=',K0H1D
        !    WRITE(*,*)'K0H1P=',K0H1P
        !    WRITE(*,*)'KMC1=',KMC1
        !    WRITE(*,*)'K1H1D=',K1H1D
        !    WRITE(*,*)'K1H1P=',K1H1P
        !    WRITE(*,*)'K2H2D=',K2H2D
        !    WRITE(*,*)'K2H2P=',K2H2P
        !    WRITE(*,*)'J1 Si=',J1
        !    WRITE(*,*)'PIE1=',PIE1
        !    WRITE(*,*)'PIE2=',PIE2
        !    READ(*,*)
        !ENDIF             
     

    !******* Silica dissolution kinetics

     FD2 = 1./(1.+M2*PIE2)
     
     K3  = XKSI*PSISED/(PSISED+KMPSI)*FD2     !WLong notes eqn (4-6-3) (3-36) (3-36-3)
     
     !particle fraction in overlying water column 
     PF  = KADSA*SSI(I,KWC)/(1.+KADSA*SSI(I,KWC))  
     
     J2  = XKSI*PSISED/(PSISED+KMPSI)*CSISAT  ! &        !WLong (4-6-4)
          
        !+ WSSINETMMD*PF*SIAT(I,KWC)   !source of Silicate due to suspended sediments settling (mgSi/m^2/day)

        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'FD2 Si=',FD2
        !    WRITE(*,*)'K3=',K3
        !    WRITE(*,*)'PF=',PF
        !    WRITE(*,*)'J2 Si=',J2
        !    READ(*,*)
        !ENDIF             
     
     IF(SFEEDER)THEN
          J2=J2 + SA_SED_SF(I)                           ! Suspension Feeders
     ENDIF           
           
           !
           !Here PF*SIAT(I,KWC) is the Silicate concentration (gSi/m^3) that 
           !is adsorbed on suspended solids SSI (gSolids/m^)
           !PF is ratio of particulate silicate concentration vs dissolved silicate concentration 
           !based on adsorption kinetics (steady state, isotherm)
           !
      !IDEBUG=13
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

     
      CALL SEDTSFNL (SI1,SI2,SIT1,SIT2,SI1TM1,SIT2TM1,TF_SSTATE)
     
      !IDEBUG=14
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

     JSI = S*(SI1-SI0)
     
        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'SI1=',SI1
        !    WRITE(*,*)'SI2=',SI2
        !    WRITE(*,*)'SIT1=',SIT1
        !    
        !    WRITE(*,*)'SIT2=',SIT2
        !    
        !    WRITE(*,*)'SI1TM1=',SI1TM1
        !    
        !    WRITE(*,*)'SIT2TM1=',SIT2TM1
        !    
        !    READ(*,*)
        !ENDIF             
     
!calculate PO4 in sediments
        
     K0H1D = 0.
     K0H1P = 0.
     KMC1  = 0.
     K1H1D = S
     K1H1P = 0.
     K2H2D = 0.
     K2H2P = 0.
     
    
     J1    = S*PO40 
#if defined(SEDIMENT_DOM)	 
	 J1 = J1 + JPX_LAYER1! add in DOP layer 1 remineralizatoin flux
#endif
     K3    = 0.
     
     PF    = KADPO4*SSI(I,KWC)/(1.+KADPO4*SSI(I,KWC))
     
     PPO4  = PF*PO4(I,KWC)                        !PO4 adsorbed to SSI in overlying water (gP/m^3) 
                                                !and that will get into sediments through SSI settling
     
     J2    = JPX !                                &!diagenesis of PO4 from POP decay (mgP/m^2/day)
            !+ WSSINETMMD*PPO4                     &!settling of dissolved PO4 adsorbed on to Suspended solids (mgP/m^2/day)
            !+ WSSINETMMD*PIP(I,KWC)                 !settling of particulate inorganic phosphorus (mgP/m^2/day)
            
     IF(SFEEDER)THEN
            J2=J2+ PIP_SED_SF(I) ! Suspension Feeders
     ENDIF            

    !Salinity dependence of DPIE1PO4 which is  (delta pi PO4,1 ) in eqn (16A.18) of DMD page 191

     IF (SAL < SALTSW) THEN
       DPIE1PO4=DPIE1PO4F        !use fresh water incremental partitioning coefficient     
     ELSE
       DPIE1PO4=DPIE1PO4S
     ENDIF

    ! Oxygen dependency of pie1

     IF (O20 < O2CRITPO4) THEN
       PIE1 = PIE2PO4*DPIE1PO4**(O20/O2CRITPO4)      !DMD (16A.19)
     ELSE
       PIE1 = PIE2PO4*DPIE1PO4                        !DMD (16A.18)
     ENDIF
     
     !write(*,*)'1/(1+m1*PIE1), PIE1,m1'
     !write(*,*)1/(1+m1*PIE1), PIE1,m1
     !READ(*,*)
     
     PIE2 = PIE2PO4 !(L/kg), partitioning coef of PO4
     
     !!IDEBUG=15
     !!IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
    !    IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'PIE1=',PIE1
        !    WRITE(*,*)'PIE2=',PIE2
        !    
        !    WRITE(*,*)'K0H1D=',K0H1D
        !    WRITE(*,*)'K0H1P=',K0H1P
        !    WRITE(*,*)'KMC1=',KMC1
        !    WRITE(*,*)'KL12=',KL12
        !    WRITE(*,*)'K1H1D=',K1H1D
        !    WRITE(*,*)'K1H1P=',K1H1P
        !    WRITE(*,*)'K2H2D=',K2H2D
        !    WRITE(*,*)'K2H2P=',K2H2P
        !    WRITE(*,*)'J1 PO4=',J1
        !    WRITE(*,*)'K3 PO4=',K3
        !    WRITE(*,*)'PF PO4=',PF
        !    WRITE(*,*)'PPO4=',PPO4
        !    WRITE(*,*)'J2 PO4=',J2
        !    WRITE(*,*)'S=',S
        !    !READ(*,*)
        !ENDIF             
     
     CALL SEDTSFNL (PO41,PO42,PO4T1,PO4T2,PO41TM1,PO4T2TM1,TF_SSTATE)
     
     
     !IDEBUG=16
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

     JPO4 = S*(PO41-PO40)

        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'PIE1=',PIE1
        !    WRITE(*,*)'PIE2=',PIE2
        !    WRITE(*,*)'PO41=',PO41
        !    WRITE(*,*)'PO42=',PO42
        !    WRITE(*,*)'PO4T1=',PO4T1
        !    WRITE(*,*)'PO4T2=',PO4T2
        !    WRITE(*,*)'PO41TM1=',PO41TM1
        !    WRITE(*,*)'PO4T2TM1=',PO4T2TM1
        !    WRITE(*,*)'JPO4=',JPO4
        !    READ(*,*)
        !ENDIF             

!******* Assign flux-flux results to mod_wqm arrays

     BENDO(I)   = - SOD                          !gO2/m^2/day
	 
     IF(SFEEDER)THEN
        BENDO(I)   = BENDO(I) - SOD_SED_SF(I)        !suspension feeders
     ENDIF
                 
                  
     MTVEL(I)   = SOD/O20                        !mass transfer velocity (m/day)
     
     !calculate layer 1 thickness (m)
     HSED1(I)=KL12*H2/S                    !Estimated thickness of layer 1 (H1) (Di Toro Appendix B) Book page 576
                                            ![m/d ]* [m] /[m/d] ==> [m]
     
!	 write(*,*)' S before DOM = ',S  ! B CLARK S_SHARE GETS PASSED TO THE MOD_SED_DOM 
#if defined(SEDIMENT_DOM)	 
   ! S_SHARE(I) = S/HSED1(I)

   !   IF(D(I)*DZ(KBM1)<=0.055)THEN
   !    S=1.e-10;
   !   ENDIF
	 S_SHARE(I) = S/((HSED1(I)+D(I)*DZ(KBM1))*0.5)   
	 
	 HSED1_SHARE(I)=HSED1(I)
	 
	 !S_SHARE=0.1
#endif    
	 
	 !IDEBUG=161
     !IF(MSR)THEN
     
        !WRITE(*,*)'DEBUG=',IDEBUG
        !WRITE(*,*)'H1=',HSED1*1000,'mm'
        !WRITE(*,*)'H2=',H2*1000,'mm'
        !WRITE(*,*)'KL12=',KL12
        !WRITE(*,*)'S=',S
        !WRITE(*,*)'O20=',O20
        !WRITE(*,*)'SOD=',SOD
        !READ(*,*)
     !ENDIF

	 !
     !BENNH4 is in gN/m^2/day
     !where JNH4 is in mgN/m^2/day 
	 !
	 
     BENNH4(I)  = JNH4/1000.                
     IF(SFEEDER)THEN
        BENNH4(I)= BENNH4(I)  + JNH4_SED_SF(I)/1000.             ! suspension feeders     
     ENDIF
     
     !IDEBUG=162
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

     BENNO3(I)  = JNO3/1000.   !nitrate flux calculated basedon KL01*(NO31-NO30) i.e. diffusion to water column converted from mgN/m^2/day to gN/m^2/d and saved in BENNO3
								!then BENNO3 is added to bottom concentration solver of water column as source term
     BENPO4(I)  = JPO4/1000.                
     IF(SFEEDER)THEN
         BENPO4(I) =BENPO4(I)  + JPO4_SED_SF(I)/1000.             ! suspension feeders
     ENDIF
     
  

     BURIALN(I) = (PON1+PON2+PON3+NO3T2+NH4T2)*W2 !mgN/m^3*m/d ==> mgN/m^2/d
     !IDEBUG=167
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
     
     BURIALP(I) = (POP1+POP2+POP3+PO4T2)*W2
     !IDEBUG=168
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
     
     BURIALC(I) = (POC1+POC2+POC3)*W2
     
     !IDEBUG=169
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG
    
    ! Diagenesis of carbon forms (gC/m^2/day)

     DIAGENC(I) = JCX/1000.        !gC/m^2/day 
     
     !IDEBUG=170
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

     !IDEBUG=17
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

        IDEBUG_SED=IDEBUG_SED+1
        !IF(MSR)THEN        
        !    WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
        !    WRITE(*,*)'BENDEN=',BENDEN(I)
        !    WRITE(*,*)'BURIALN=',BURIALN(I)
        !    WRITE(*,*)'BURIALC=',BURIALC(I)
        !    WRITE(*,*)'BURIALP=',BURIALP(I)
        !    WRITE(*,*)'DIAGENC=',DIAGENC(I)
        !    READ(*,*)
        !ENDIF             
     
     IF (BALGAE_CALC) THEN
     
        CALL BA_CALC

!******* Phosphorus available for benthic algal growth

!WLong: adjust BENNH4 BENNO3, BENPO4, BENDOC, BENDO etc flux based on benthic algae
       BENNH4(I) = BENNH4(I) + BANH4(I)  !(gN/m2/day)
       BENNO3(I) = BENNO3(I) + BANO3(I) 
       
       BENPO4(I) = BENPO4(I) + BAPO4(I)     !(gP/m^2/day)
       
      !  BENDOC(I) = BENDOC(I) + BADOC(I)  !(gC/m^2/day)
       BENDO(I)  = BENDO(I)  + BADO(I)   !(gO2/m^2/day)

! COMPUTE EFFECTS OF ALGAL ACTIVITY ON ORGANIC PARTICULATES (MG/M**3)


!WLong : acumulate POCi using DTPOCi_BA returned from MOD_BA
       POC1 = POC1 + DTPOC1_BA(I)    !(mgC/m^3)
       POC2    = POC2 + DTPOC2_BA(I)    !(mgC/m^3)
       POC3 = POC3 + DTPOC3_BA(I)    !(mgC/m^3)
       
!Wlong: accumulate PONi using DTPONi_BA returned from MOD_BA       
       PON1 = PON1 + DTPON1_BA(I)      !mgN/m^3
       PON2 = PON2 + DTPON2_BA(I)     !mgN/m^3
       PON3    = PON3 + DTPON3_BA(I)    !mgN/m^3
       
!WLong: accumulate POPi using DTPOPi_BA returned from MOD_BA       
       POP1 = POP1 + DTPOP1_BA(I)
       POP2 = POP2 + DTPOP2_BA(I)
       POP3 = POP3 + DTPOP3_BA(I)
       
!******* Accumulate fluxes for steady-state computation for G3

       IF (STEADY_STATE_SED_G3) THEN
        
         AG3CFL(I) = AG3CFL(I)+DTAG3CFL_BA(I) !1000.*PRB(I)*FRCPHB(3)*BBM(I)*DLTS         !(mgC/m^2/day)
         AG3NFL(I) = AG3NFL(I)+DTAG3NFL_BA(I) !1000.*PRB(I)*FRNPHB(3)*ANCB*BBM(I)*DLTS    !(mgN/m^2/day
         AG3PFL(I) = AG3PFL(I)+DTAG3PFL_BA(I) !1000.*PRB(I)*FRPPHB(3)*APCB*BBM(I)*DLTS  !(mgP/m^2/day)
       
       ENDIF

     ENDIF
     
     
     
! TEMPORARY FIX UP TO EXAMINE EFFECT OF SAV ON SEDIMENTS
66666   CONTINUE
!        NH4T2 = NH4T2TM1
!        PO4T2 = PO4T2TM1
!        NO3T2 = NO3T2TM1
!        HST2  = HST2TM1
! END TEMPORARY FIX UP

!******* Total sediment nutrient mass

!Wen Long     SEDMN = SEDMN+(PON1+PON2+PON3+NH4T2+NO3T2)*SFA(I)*H2/1.E6
!Wen Long     SEDMP = SEDMP+(POP1+POP2+POP3+PO4T2)*SFA(I)*H2/1.E6
!Wen Long     SEDMC = SEDMC+(POC1+POC2+POC3)*SFA(I)*H2/1.E6

     SEDMN = SEDMN+(PON1+PON2+PON3+NH4T2+NO3T2)*ART1(I)*H2/1.E6  !(kgN) Wen Long replaced SFA by ART1 
     SEDMP = SEDMP+(POP1+POP2+POP3+PO4T2      )*ART1(I)*H2/1.E6  !(kgP) Wen Long replaced SFA by ART1
     SEDMC = SEDMC+(POC1+POC2+POC3            )*ART1(I)*H2/1.E6  !(kgC) Wen Long replaced SFA by ART1

     !                (mgC/m^3) *(m^2)*m /1000/1000 = mgC/1000/1000 = gC/1000 = kgC
     
   
   !cycle the new clauclations to n-1 time step for next step marching
     
     
     BENSTRTM1S(I)= BENSTR
     
     SODTM1S(I)   =SOD
     JNH4TM1S(I)  =JNH4
     JNO3TM1S(I)  =JNO3
     JPO4TM1S(I)  =JPO4
     JCH4TM1S(I)  =JCH4
     JCH4GTM1S(I) =JCH4G
     JHSTM1S(I)   =JHS
     JSITM1S(I)   =JSI
                                            
     PON1TM1S(I)  = PON1
     PON2TM1S(I)  = PON2
     PON3TM1S(I)  = PON3
     !IF(MSR)write(*,*)'LBnote iv: POC2TM1S(38)=',POC2TM1S(38),'  POC2=',POC2
     POC1TM1S(I)  = POC1
     POC2TM1S(I)  = POC2
     POC3TM1S(I)  = POC3
     !IF(MSR)write(*,*)'LBnote v: POC2TM1S(38)=',POC2TM1S(38),'  POC2=',POC2
     
     POP1TM1S(I)  = POP1
     POP2TM1S(I)  = POP2
     POP3TM1S(I)  = POP3
     
     PSISEDTM1S(I)   = PSISED    !sediment particulate organic silicate

     NH41TM1S(I)  = NH41        !dissolved ammonia in layer 1
     NH42TM1S(I)  = NH42        !dissolved ammonia in layer 2
     NH4T2TM1S(I) = NH4T2        !total ammonia in layer 2
     
     NO31TM1S(I)  = NO31        !dissolved nitrate in layer 1
     NO32TM1S(I)  = NO32        !dissolved nitrate in layer 2
     NO3T2TM1S(I) = NO3T2        !total nitrate in layer 2
     
     PO41TM1S(I)  = PO41        !dissolved PO4 in layer 1
     PO42TM1S(I)  = PO42        !dissolved PO4 in layer 2
     PO4T2TM1S(I) = PO4T2        !total PO4 in layer 2

     HS1TM1S(I)   = HS1            !dissolved HS in layer 1
     HS2TM1S(I)   = HS2            !dissolved HS in layer 2
     HST2TM1S(I)  = HST2        !total HS in layer 2
     
     SI1TM1S(I)   = SI1         !dissolved Si in layer 1
     SI2TM1S(I)   = SI2         !dissolved Si in layer 2
     SIT2TM1S(I)  = SIT2        !total Si in layer 2
     
     !IDEBUG=18
     !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

     !WLong added POC1TM1S_SED_DF
     IF(DFEEDER)THEN
        POC1TM1S_SED_DF(I) = POC1
        POC2TM1S_SED_DF(I) = POC2
     ENDIF        
     
     BFORMAXS(I)  = BFORMAX
     ISWBENS(I)   = ISWBEN
     
!WLong     DFEEDM1S(I)  = DFEED   !moved to mod_df.F
     CH41TM1S(I)  = CH41                ! dissolved CH4 in layer 1
     CH42TM1S(I)  = CH42                ! dissolved CH4 in layer 2
     CH4T2TM1S(I) = CH4T2               ! total CH4 in sediment layer 2
     SO4T2TM1S(I) = SO4T2               ! SO4
     
     DIAGN(I)     = JNX                 ! CFC

     IDEBUG_SED=IDEBUG_SED+1
     !IF(MSR)THEN        
    !        WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'CH41=',CH41
    !        WRITE(*,*)'CH4T2=',CH4T2
    !        WRITE(*,*)'SO4T2=',SO4T2
    !        WRITE(*,*)'JNX=',JNX
    !        READ(*,*)
     !ENDIF             
   ENDDO

!***** Assign concentrations to plot variables !
   !IF(MSR)write(*,*)'LBnote vi: POC2TM1S(38)=',POC2TM1S(38),'  CPOC(38,2)=',CPOC(38,2)
   DO I=1,MLOC
     CPOC(I,1) = POC1TM1S(I)
     CPOC(I,2) = POC2TM1S(I)
     CPOC(I,3) = POC3TM1S(I)
     CPON(I,1) = PON1TM1S(I)
     CPON(I,2) = PON2TM1S(I)
     CPON(I,3) = PON3TM1S(I)
     CPOP(I,1) = POP1TM1S(I)
     CPOP(I,2) = POP2TM1S(I)
     CPOP(I,3) = POP3TM1S(I)
     
     CNH4(I)   = NH4T2TM1S(I)
     CNO3(I)   = NO3T2TM1S(I)
     CPO4(I)   = PO4T2TM1S(I)

     CPOS(I)   = PSISEDTM1S(I)
     CCH4(I)   = CH4T2TM1S(I)
     CSO4(I)   = SO4T2TM1S(I)
     CHS(I)    = HST2TM1S(I)  
     CSI(I)    = SIT2TM1S(I)  
   ENDDO
    !IF(MSR)write(*,*)'LBnote c: CPOC_GL(38,2)=',CPOC_GL(38,2), '  CPOC(38,2)=',CPOC(38,2)
    
   !IDEBUG=19
   !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

!***** Take temperature integration step

   DO I=1,MLOC
   
     IDEBUG_SED=IDEBUG_SED+1
     !IF(MSR)THEN        
    !        WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'CTEMP(I)=',CTEMP(I)
    !        WRITE(*,*)'HSED(I)=',HSED(I)
    !        WRITE(*,*)'T(I,KBM1)=',T(I,KBM1)
    !        READ(*,*)
     !ENDIF                              
     
!---Wen Long commented out calculation of sediment temperature through vertical diffusion from water column ----     
     
!     IF(DLT*DIFFT/HSED(I)/HSED(I)>0.5)THEN
!        IF(MSR)THEN
!            WRITE(*,*)'Oops time step does not meet stability criterion for solving sediment temperature' 
!        ENDIF
        
!        WRITE(*,*)'statbility criterion: DLT*DIFFT/HSED(I)/HSED(I)=',DLT*DIFFT/HSED(I)/HSED(I)
!        STOP 'time step too large for sediment temperature calculation'
!
!     ENDIF     
!     CTEMP(I) = CTEMP(I)+DLT*DIFFT/HSED(I)/HSED(I)          &  !WLong: discretization scheme of this
!                 *(T(I,KBM1)-CTEMP(I))                         !       needs to be checked
!                                                               !Note here DIFFT has unit m^2/s so DLT is in sec

    
     
     CTEMP(I)=T(I,KBM1)     !Use water column temperature directly
     
     IDEBUG_SED=IDEBUG_SED+1
     !IF(MSR)THEN        
    !        WRITE(*,*)'JDAY=',JDAY,'I=',I,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'CTEMP(I)=',CTEMP(I)
    !        WRITE(*,*)'HSED(I)=',HSED(I)
    !        WRITE(*,*)'T(I,KBM1)=',T(I,KBM1)
    !        WRITE(*,*)'DLT=',DLT
    !        WRITE(*,*)'statbility criterion: DLT*DIFFT/HSED(I)/HSED(I)=',DLT*DIFFT/HSED(I)/HSED(I)
    !        READ(*,*)
     !ENDIF
   ENDDO

   !IDEBUG=20
   !IF(MSR)WRITE(*,*)'DEBUG=',IDEBUG

   RETURN
   END SUBROUTINE SED_CALC
   
   SUBROUTINE SED_DIAGENESIS_G3()

   INTEGER :: I

!***** Compute time-average values 

!   DO I=1,MLOC
!     AG3CFL(I) = AG3CFL(I)/TINTIM
!     AG3NFL(I) = AG3NFL(I)/TINTIM
!     AG3PFL(I) = AG3PFL(I)/TINTIM
!     ASDTMP(I) = ASDTMP(I)/TINTIM
!   ENDDO

!***** Compute G3 organic concentrations

   DO I=1,MLOC
     CPOC(I,3) = AG3CFL(I)/(KCDIAG(3)*DCTHTA(3)**(ASDTMP(I)-20.)   &
                  *HSED(I)+VSED(I))
     CPON(I,3) = AG3NFL(I)/(KNDIAG(3)*DNTHTA(3)**(ASDTMP(I)-20.)   &
                  *HSED(I)+VSED(I))
     CPOP(I,3) = AG3PFL(I)/(KPDIAG(3)*DPTHTA(3)**(ASDTMP(I)-20.)   &
                  *HSED(I)+VSED(I))
   ENDDO
   
   !IF(MSR)write(*,*)'LBnote d: CPOC_GL(38,2)=',CPOC_GL(38,2), '  CPOC(38,2)=',CPOC(38,2)
   
   RETURN
   
   END SUBROUTINE SED_DIAGENESIS_G3


!********************************************************************************
!**                          F U N C T I O N   S E D F                         **
!********************************************************************************

   FUNCTION SEDF(SOD1, DF_SOD_TMP,SAV_SOD_TMP,TF_SSTATE) RESULT(SEDFOUT)

   !
   !Function that calculates the imbalance of equation
   !
   ! SOD=CSOD(SOD1)+NSOD(SOD1) + DFSOD(SOD1) + SAVSOD(SOD1)
   !
   !where SOD1 is given SOD, and CSOD is calculated using CSOD_HS , CSOD_CH4
   !NSOD is calculated as NSOD_NH4
   !
   !     CSOD = CSOD_HS  for salt > SALTSW (salt water system)
   !or
   !     CSOD = CSOD_CH4 for salt < SALTSW (fresh water system)
   !
   !     NSOD = oxygen consumed during nitrification of NH4 to NO3 (DMD book eqn (1.25) reaction)
   ! 
   !Imbalance = SOD-SOD1
   !
   !Aim of the program is to find SOD=SOD1 so that Imbalance is zero, this is done through zbrent() function
   !
   
   !
   !First of all,give any SOD1 value (gO2/m^2/day), calculate NO3, NH4, SO4, CH4, HS based on XJC1
   !
   !Then caluclate CSOD and NSOD, and put together CSOD+NSOD+DFSOD+SAVSOD
   !
   !Finally calculate the Imbalance SEDF=SOD-SOD1=(CSOD+NSOD+DFSOD+SAVSOD) - SOD1
   !
!   USE MOD_PAR, ONLY: NGID

   IMPLICIT NONE
                          
   SAVE    
   LOGICAL  ::  TF_SSTATE
   INTEGER :: I
   REAL(SP) ::     SOD1,             &    !guessed SOD (gO2/m^2/day)
                DF_SOD_TMP,     &    !SOD due to deposption feeder
                SAV_SOD_TMP,    &    !SOD due to SAV root respiration
                SEDFOUT                !returned error SOD-SOD1 for checking convergence (gO2/m^3/d)
                
   REAL(SP) ::     DCH4T2,            &    !rate of total  CH4 generation in layer 2 (gO2/m^3 * m/d  = gO2/m^2/d)
                DHST2                !rate of HS generation in layer 2 (gO2/m^3/d * m/d = gO2/m^2/d)

   REAL(SP) ::  SAVSOD_SED,        &
                DFSOD_SED
   
! KURT GLAESEMANN Changed H to HH to avoid name conflict with module
   REAL(SP) AD(4,4), BX(4), G(2), HH(2,2)
   
   REAL(SP) DBLSO41,     &
            DBLSO42,     &
            RA0,         &
            RA1,         &
            RA2,         &
            R1,         &
            R2,         &
            DISC,         &
            SN1

   REAL(SP) CSODMAX,    &    !maximum value of CSOD (gO2*/m^2/d)
            SECH_ARG        !argument of hyperbolic secant function (no unit)
            
!***** Compute the NH4, NO3, and SOD fluxes

   !surface diffusion velocity (m/day)

   S = SOD1/O20        !(mgO2/m^2/day) / (mgO2/m^3) ==> (m/day)

!***** Ammonia flux

   !layer integrated reaction rates (m/d)
   K0H1P = 0. 
   K1H1P = 0.
   
   K2H2D = 0.
   K2H2P = 0.
   
   IF (KMNH4 /= 0.) THEN                !KMNH4 is K_{NH4} in WLong notes (4-11)
     K0H1D = XAPPNH4**2/S                &
            *KMNH4*THTAKMNH4**(STP20)    &    !WLong fixed based on Excel version
            *(O20/(2*KMNH4O2+O20))          !O20/(KMNH4O2+O20) is F_{O2} in WLong notes (4-10)
                                                !WLong (3-32-8)
     K1H1D = S                                    !WLong (3-32-6)
   ELSE

     K0H1D = 0.                                      !WLong (3-32-7)
     K1H1D = XAPPNH4**2/S*(O20/(2*KMNH4O2+O20))+S !WLong (3-32-5)
   ENDIF
   
   J1   = S*NH40                                !(m/d) * (mgN/m^3) ==> mgN/m^2/day WLong (4-6) 
   
#if defined(SEDIMENT_DOM) 
   ! Potentially at in NH4 flux from DON remineralization in layer 1
  ! write(*,*)' J1 before = ',J1
   
   J1 = J1 + JNX_LAYER1
   
 !  write(*,*)' J1 after = ',J1
   
#endif  
   
   K3   = 0.                                    !WLong (4-8)
   J2   = JNX                                    !WLong (4-9)
   PIE1 = PIENH4
   PIE2 = PIENH4
   KMC1 = KMNH4*THTAKMNH4**(STP20) !Michael=Mention half constant for layer 1 reaction, DMD eqn (3.30)

     IDEBUG_SED=IDEBUG_SED+1
     !IF(MSR)THEN        
    !        WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'S=',S
    !        WRITE(*,*)'SOD1=',SOD1
    !        WRITE(*,*)'K0H1D=',K0H1D
    !        WRITE(*,*)'K1H1D=',K1H1D
    !        WRITE(*,*)'KMNH4=',KMNH4
    !        WRITE(*,*)'J1=',J1
    !        WRITE(*,*)'K3=',K3
    !        WRITE(*,*)'J2=',J2
    !        WRITE(*,*)'PIE1=',PIE1
    !        WRITE(*,*)'PIE2=',PIE2
    !        WRITE(*,*)'KMC1=',KMC1
    !        READ(*,*)
     !ENDIF    
      
   CALL SEDTSFNL (NH41,NH42,NH4T1,NH4T2,NH41TM1,NH4T2TM1,TF_SSTATE)
   
   IF(TF_SSTATE)THEN !for steady state make sure NH41TM1=NH41 just calculated
      NH41TM1=NH41      !WLong: this is important for steady state
   ENDIF
   
  ! write(*,*)'S before NH4 = ',S
   
   JNH4 = S*(NH41-NH40)            !(m/d)*(mgN/m^3) ==> mgN/m^2/d
   
   IDEBUG_SED=IDEBUG_SED+1
   !IF(MSR)THEN        
!            WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
!            WRITE(*,*)'NH41=',NH41
!            WRITE(*,*)'NH42=',NH42
!            WRITE(*,*)'NH4T1=',NH4T1
!            WRITE(*,*)'NH4T2=',NH4T2
!            WRITE(*,*)'NH41TM1=',NH41TM1
!            WRITE(*,*)'NH4T2TM1=',NH4T2TM1
!            WRITE(*,*)'NH40=',NH40
!            WRITE(*,*)'JNH4=',JNH4
!            READ(*,*)
   !ENDIF    

!***** Oxygen consumed by nitrification (NSOD)

  
   A1 = 0.0045714                            !WLong notes (4-14-1) see DMD book eqn (1.25) and p206,p 189
   
  !calculate NSOD (gO2/m^2/d)
   
   IF (KMNH4 /= 0.) THEN                    
     JO2NH4 = A1*K0H1D*NH41/(KMNH4*THTAKMNH4**(STP20)+NH41TM1)
                            !WLong (3-32-14) (4-14)
                            !Note oxygen regulation is already in K0H1D, DMD book eqn (9.12)(9.14)
   ELSE                                        
     JO2NH4 = A1*(K1H1D-S)*NH41  !WLong notes (4-14-2) (3-32-15) 
   ENDIF                        
                            !(mgO2/mgN)* (m/d) *(mgN/m^3)/1000 ==> mgO2/m^2/d/1000 ==> gO2/m^2/d
   
   IDEBUG_SED=IDEBUG_SED+1
   !IF(MSR)THEN        
!            WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
!            WRITE(*,*)'JO2NH4=',JO2NH4
!            WRITE(*,*)'A1=',A1        
!            WRITE(*,*)'KMNH4=',KMNH4                    
!            WRITE(*,*)'K1H1D=',K1H1D                                
!            READ(*,*)
   !ENDIF    

!***** Denitrification 

   K0H1D = 0.                !WLong (3-32-7)
   K0H1P = 0.                !WLong (3-32-9)
   KMC1  = 0.
   K1H1D = XAPP1NO3**2/S+S  !WLong (3-32-5)
   K1H1P = 0.                !WLong (3-32-10, should not be zero but equivalent to zero since FP1 = 0 for NO3)
                            !    see WLong(3-32-2), it is multiplied by FP1 anyway
   K2H2D = XK2NO3            !WLong (3-36-7) 
   K2H2P = 0.                !WLong (3-36-8)
   IF (KMNH4 /= 0.) THEN    
     J1 = S*NO30+XAPPNH4**2/S                                                &
            *KMNH4*THTAKMNH4**(STP20)/(KMNH4*THTAKMNH4**(STP20)+NH41TM1)        &
            *(O20/(2*KMNH4O2+O20))*NH41                        !WLong (4-40-1) with KMNH4 !DMD(3.34)
          
   ELSE
     J1 = S*NO30+XAPPNH4**2/S*(O20/(2*KMNH4O2+O20))*NH41    !WLong (4-40-1) with no KMNH4
   ENDIF
   K3   = 0.
   J2   = 0.
   PIE1 = 0.
   PIE2 = 0.

   IDEBUG_SED=IDEBUG_SED+1
   !IF(MSR)THEN        
!             
!            WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
!            WRITE(*,*)'S=',S
!            WRITE(*,*)'SOD1=',SOD1
!            WRITE(*,*)'K0H1D=',K0H1D
!            WRITE(*,*)'K1H1D=',K1H1D
!            WRITE(*,*)'KMNH4=',KMNH4
!            WRITE(*,*)'J1=',J1
!            WRITE(*,*)'K3=',K3
!            WRITE(*,*)'J2=',J2
!            WRITE(*,*)'PIE1=',PIE1
!            WRITE(*,*)'PIE2=',PIE2
!            WRITE(*,*)'KMC1=',KMC1
!            READ(*,*)
   !ENDIF    
   
   CALL SEDTSFNL(NO31,NO32,NO3T1,NO3T2,NO31TM1,NO3T2TM1,TF_SSTATE)
   JNO3 = S*(NO31-NO30)            !WLong (4-44)  (S is the diffusion coefficient (m/d) unit)
   !IF(MSR)WRITE(*,*)'LBnote: NO30=',NO30  !,'  S=',S
   !IF(MSR)WRITE(*,*)'        Jno3=',JNO3,' NO31=',NO31

   IDEBUG_SED=IDEBUG_SED+1
   !IF(MSR)THEN        
!            WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
!            WRITE(*,*)'NO31=',NO31
!            WRITE(*,*)'NO32=',NO32
!            WRITE(*,*)'NO3T1=',NO3T1
!            WRITE(*,*)'NO3T2=',NO3T2
!            WRITE(*,*)'NO31TM1=',NO31TM1
!            WRITE(*,*)'NO3T2TM1=',NO3T2TM1
!            WRITE(*,*)'NO30=',NO30
!            WRITE(*,*)'JNO3=',JNO3
!            READ(*,*)
   !ENDIF

!***** Sulfide/methane oxidation

   AO2N      = 2.857142857142857       !AO2N= 2.85714 is alpha_O2,C and alpha_C,N in (4-56) of Wen Long notes
                                    !alpha_O2C = 2.667gO2/gC (4-46) of WLong notes
                                    !alpha_C,N = 1.071 gC/gN (4-45) of WLong notes
                        
                                    != 2.0*16/12.0*10/8.0*12/14.0 
   
                                    !
                                    !==> AO2N=2.667*1.071 = 2.85714 gO2/gN
                                    !
                        
   XJCNO31  = AO2N*XAPP1NO3**2/S*NO31/1000.0  !layer 1 denitrification flux (gO2/m^2/d) (WLong (4-56) RHS 2nd term)
   XJCNO32  = AO2N*XK2NO3*NO32/1000.0          !layer 2 denitrification flux (gO2/m^2/d) (WLong (4-57) RHS 2nd term)
                                              !mgO2/mgN  * m/d * mgN/m^3 /1000= mgO2/m^2/d /1000 = gO2/m^2/d
                
   
   
!***** Add the aerobic and first anaerobic layer to keep mass balance

   XJCNO3 = XJCNO31+XJCNO32        !JO2NO3T in Greg's excel version

!***** Convert carbon diagenesis flux to O2 units

   XJC1 = AMAX1(2.667E-3*JCX-XJCNO3,1.0e-10)  !WLong, note here XJC1 is in gO2/m^2/day as 2.667 converts from C to O2
                                               !and 1.0E-3 converts from mgC to gC
                                               !WLong notes eqn(4-57) + eqn (4-56) for both layer 1 and layer 2
   
   !XJC1 is JH2S,1+JH2S,2 (i.e. WLong (4-56) , (4-57) together ) 
   
    IDEBUG_SED=IDEBUG_SED+1
    !IF(MSR)THEN        

    !        WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'AO2N=',AO2N
    !        WRITE(*,*)'XJCNO31=',XJCNO31
    !        WRITE(*,*)'XJCNO3=',XJCNO3
    !        WRITE(*,*)'XK2NO3=',XK2NO3
    !        WRITE(*,*)'XJC1=',XJC1
    !        READ(*,*)

    !ENDIF    
    
#if(1) /* James L. Martin's method of calculating HS and CH4 : solve HS if SAL<SALSW using DMD chapter 10 */
   
    IF(SAL>=SALTSW)THEN
    
            !dissolved and particle fractions for HS
   
            PIE1=PIE1HS    
            PIE2=PIE2HS
   
            KMC1=0.0
            K0H1D=0.                                        !WLong (3-32-7)
            K0H1P=0.                                        !WLong (3-32-9)
   
            !   K1H1D=XAPP1HSD**2/S*(O20/KMHSO2) + S    
            !   K1H1P=XAPP1HSP**2/S*(O20/KMHSO2)
   
            K1H1D=XAPP1HSD**2/S*(O20/(KMHSO2)) + S        !WLong notes (4-58)
                                                        !WLong (3-32-5)
            K1H1P=XAPP1HSP**2/S*(O20/(KMHSO2))            !WLong (3-32-10)
   
            K2H2D=0.        !WLong (3-36-8)
            K2H2P=0.        !WLong (3-36-8-1)
   
            J1=0.        !WLong (4-52) J1 term with J1=S*fd0*HST0 + JH2ST,1, where HST0=0, JH2ST,1 is put in J2
                        !instead
   
            K3=0.0        !WLong (3-36-6)
   
            J2=XJC1        !gO2/m2/day, carbon diagenesis in layer 2 in oxygen units after discounting diaggenesis
                        !WLong (4-57) + (4-56) with JH2ST,1 and JH2ST,2 are put in JH2ST,2 together and ignoring 
                !layer 1 
   
            CALL SEDTSFNL(HS1,HS2,HST1,HST2,HS1TM1,HST2TM1,TF_SSTATE) 
   
            JHS        =    S*(HS1-HS0)         !gO2/m^2/day (flux of dissolved HS to overlying water column)
                                    !WLong (4-61)
    
            FD1=1./(1.+M1*PIE1)
            FP1=M1*PIE1/(1.+M1*PIE1)
            FD2=1./(1.+M2*PIE2)
            FP2=M2*PIE2/(1.+M2*PIE2)  
   
   
            !CSOD due to HS oxidation in layer 1 limited by avaiable O2 in overlying water
            ! (m/d)*(gO2/m^3)*(1) ==> unit is gO2/m^2/d
   
            
!WLong fix slight error, should use XAPP1HSD here
!            CSODHS= (  XAPP1HSP**2/S*FD1*HST1                 &!             
            
            CSODHS= (  XAPP1HSD**2/S*FD1*HST1                 & !oxidation of dissolved HS in layer 1            
                    + XAPP1HSP**2/S*FP1*HST1                & !oxidation of particulate HS in layer 1
                    )*(O20/KMHSO2)                               !limitation by avaiable O20
                                                              !WLong (4-62)

            !No methane generation, no CSOD due to methane oxidation to CO2
            
            CSODCH4    =0.d0
            JCH4AQ    =0.d0
            JCH4G    =0.d0
            JCH4     =0.0d0      !LB: KURT GLAESSEMANN found that JCH4 was not being zeroed here!  28 April 2015
            
            CSOD    =CSODCH4+CSODHS
            
            !
            !dummy CH4 solution as no CH4 will be generated??WLong: should we retrain them to be
            !same as previous time step? consider SAL<SALTSW is altenatingly satisfied.
            !WLong: may need to set CH40 to zero as well?
            !
            
            CH4T2=0
            CH4T1=0
            CH41=0
            CH42=0
   
    ELSE    !methane generation with salt is < saltsw
    
    
            CSODMAX=MIN(SQRT(2.d0*KL12*CH4SAT*XJC1),XJC1)  !WLong 4-65)
            
            SECH_ARG=XAPPCH4/S                               !WLong (4-70), DMD book (10.39)
            
            CSODCH4=0.d0
            
                
            IF(.FALSE.)THEN    !--the analytical solution
                IF(SECH_ARG<400) THEN
                    CSODCH4=CSODMAX*(1.0-(2.0/(EXP(SECH_ARG)+EXP(-SECH_ARG))))  !CSODMAX(1-sech(sech_arg))
                            !where hyperbolic secant is : sech(x)=2/[exp(x)+exp(-x)]
                ELSE                    
                    CSODCH4=CSODMAX
                ENDIF

            ELSE            !--the numerical solution
                        
                CH41 = (CSODMAX + s * CH40) / (XAPPCH4**2/s + s)    !analytical solution of CH4(1)
                CSODCH4 = (XAPPCH4**2/s)*CH41
            ENDIF
            
            JCH4AQ=CSODMAX-CSODCH4            !DMD book (10.35a)            
            JCH4G =XJC1-JCH4AQ-CSODCH4        !DMD book (10.33)
          
			
            CH41 = (CSODMAX + s * CH40) / (XAPPCH4**2/s + s)    !analytical solution of CH4(1)            
            
            JCH4      = JCH4AQ !S*(CH41-CH40)
                        
            !
            !dummy HS solution as no HS will be generated??WLong: should we retrain them to be
            !same value as previous time step? Consider SAL>SALTSW condition is alternatingly satisfied.
            !May need to set HS0 to zero as well?
            !
            
            HST1=0
            HST2=0
            HS1=0
            HS2=0
            
            CSODHS=0.d0
            JHS=0.d0

            CSOD=CSODHS+CSODCH4

    ENDIF


    !do not solve for SO4 at all (make all SO4 terms zero)
    !Need to set SO40 to zero as well ?
    SO4T1=0
    SO4T2=0
    SO41=0
    SO42=0
    
#else /*solve HS, SO4 and CH4 together following DMD book chapter 11*/
    
    CALL GET_SO4_HS_CH4_CSOD_FLUXES(TF_SSTATE)

#endif 
   
!**** Volumetric methane and total gas flux (L/m2-d)
    VJCH4G=22.4/64.*JCH4G                  !JCH4G is mgO2/m^2/day
!   JGAS=VJN2GAS+VJCH4G                   !vjn2gas not computed    !DMD book (10.54), pp 210
    

    !
    !calculate final SOD=CSOD+NSOD=(CSOD_CH4+CSOD_HS)+NSOD  
    !where     CSOD_CH4 is SOD due to CH4 
    !        CSOD_HS is SOD due to HS
    
    SOD  = CSOD+JO2NH4                    !WLong (4-9-14)
    
    !SOD=(SODold+CSOD+JO2NH4)/2.0        !Greg's calculation of SOD is an average with previous time step

    !include SOD due to deposition feeder
    !IF(DFEEDER)THEN
        DFSOD_SED = DF_SOD_TMP            ! deposit feeders
        SOD=SOD+DFSOD_SED                ! gO2/m^2/day
    !ENDIF
    
    
    !include SOD due to SAV root respiratioon
    !IF(SAVC_ALC)THEN
      SAVSOD_SED=SAV_SOD_TMP            
      SOD=SOD+SAVSOD_SED                !gO2/m^2/day
    !ENDIF

   SEDFOUT = SOD-SOD1                    !calculate error in SOD and SOD1

   RETURN
   END FUNCTION SEDF
  
  
   
!********************************************************************************
!**                        F U N C T I O N   Z B R E N T                       **
!********************************************************************************

   FUNCTION ZBRENT(IERR, DF_SOD_TMP,SAV_SOD_TMP,TF_SSTATE, SODold) RESULT(ZBOUT)  
   !
   !Trying to solve the equation : SEDF[SOD]== SOD-(CSOD+NSOD)[SOD] =0 
   !by guessing SOD and then call SEDF, here [] means argument to function 
   !and and eventually converge to a value within the range of SODMIN and SODMAX
   !
   !In order to solve for f(SOD)=0, i.e. SEDF[SOD]=0, with f(SOD)=SOD-CSOD(SOD)+NSOD(SOD)
   !   where SOD is unknown, CSOD and NSOD are functions of SOD as well, this program ZBRENT
   !takes guess of what SOD is by using SODMIN and SODMAX initially, then narrow down the
   !search successively by moving B point, [SODMIN, A, B, SODMAX],
   !until A, B gets closer and closer with f(A),f(B) ~=0
   !
   !zbrent finds the root of a function without knowing the derivative
   !http://en.wikipedia.org/wiki/Brent%27s_method 
   
   !
   !ZBRENT returns IERR and value of SOD as ZBOUT
   !
   !if IERR = 0   good, find soultion and SOD converges and the solution is in ZBOUT=B
   !   IERR = 1    did not find solution of SOD within [SODMIN, SODMAX] range
   !   IERR = 2    did not find the solution within the IMAX interations, last point B is saved in ZBOUT
   !
     
   IMPLICIT NONE
   LOGICAL :: TF_SSTATE 
   INTEGER,PARAMETER      :: IMAX=500
   REAL(SP),PARAMETER   :: EPS=3.E-10, TOL=1.E-8
   
   INTEGER                :: IERR &    !error flag, 1 : not able to find SOD solution within SODMIN and SODMAX
                                    !             2 : good, able to find solution and the value is returned 
                                    !                by ZBRENT
                        , I
   REAL(SP)             ::     ZBOUT,      &  
!                            SODMIN,       &!minimum SOD (gO2/m^2/day) !moved to module definition
!                            SODMAX,       &!maximum SOD (gO2/m^2/day)
                            TOL1

   REAL(SP)             :: A, B, C, D, E, P, Q, R, XM
   REAL(SP)             :: FA, FB, FC
   REAL(SP)                :: DF_SOD_TMP,SAV_SOD_TMP
   REAL(SP)                :: SODold        !previous time step SOD
 
  
   !SODMIN = 1.E-4  !WLong and LB: instead of hardwire SODMIN and SODMAX here
   !SODMAX = 100.   !it is hardwire before the call to zbrent 
   
   
   
   ZBOUT  = 0.0
!***** Initialize upper and lower limits for solution

   IERR = 0
   
   !check if SODMIN and SODMAX are bracketing the soultion 
   !i.e. if there would be a positive SOD value, such that SODMIN<SOD<SODMAX
   !and  the zbrent function SEDF(SOD, DF_SOD_TMP) = 0
   !
   
   A    = SODMIN
   B    = SODMAX
   
   FA   = SEDF(A,DF_SOD_TMP,SAV_SOD_TMP,TF_SSTATE)  !sediment flux calculation (SOD) using minimum SOD
   FB   = SEDF(B,DF_SOD_TMP,SAV_SOD_TMP,TF_SSTATE)  !sediment flux calculation (SOD) using maximum SOD
 
    !***** Root must bracket ZBRENT

   IF (FB*FA > 0.) THEN !if f(A)*f(B)>0, then no solution within [A B] bracket
     IERR = 1 
     RETURN      
   ENDIF
   
   !WLong initialized C and E to prevent a undefined error using ifort 
   !The error message without initializing C and E was:
   !   Run-Time Check Failure. The variable 'mod_sed_mp_zbrent_$C' is being used without being defined
   !
   !C=0.5*(A+B)         !? WLong how to initialize this one?
   
   C    = B  !WLong followed JLM code and initialize C as B instead of (A+B)/2
                
   E    = B-A                !? WLong how to initialize this one?
            !WLong following JLM, no need to initiate E here

   !start with FC = FB, i.e. searching from SODMAX end towards SODMIN
   FC = FB
   
   !iterate IMAX times until found solution or failed
   DO I=1,IMAX
   
   
!     IF(MSR)WRITE(*,*) 'ZBRENT=3, IMAX=',IMAX,'I=',I
     IF (FB*FC > 0.) THEN    !If f(B)f(C) >0 then start lookfing for A side instead
       C  = A
       FC = FA            
       D  = B-A                !D keeps the search window width = B-A
       E  = D                !E also keeps a copy of D
     ENDIF
     
!
!makre sure |f(A)| > |f(B)| 
!
     
     
!     IF(MSR)WRITE(*,*) 'ZBRENT=4, IMAX=',IMAX,'I=',I
     IF (ABS(FC) < ABS(FB)) THEN    !if f(C) is closer to zero than f(B) then 
                                    !replace search widnow [A B] by new window [B C] 
                                    !and again start guessing root from the C side
     
                                    !                    
                                    !                  ~
                                    !                /  \
                                     !before ------A(C)-------B-----
                                    !               /     \
                                    !         .- ~ *          `_
                                    !              f(C)        \
                                    !                        * f(B)
                                    !            
                                    !               |<---------->|
                                    !                    E
                                    !
     
                                    !                  ~
                                    !                /  \
                                     !after  -------(B)----(C)A-----
                                    !               _^     \
                                    !         .- ~ *          `_
                                    !              f(B)        \
                                    !                        * f(A),f(C)
                                     !
                                    !               |<------->|
                                    !                  XM*2        !new window width
       A  = B
        
!       IF(MSR)WRITE(*,*) 'ZBRENT=71, IMAX=',IMAX,'I=',I
     
       B  = C
!       IF(MSR)WRITE(*,*) 'ZBRENT=72, IMAX=',IMAX,'I=',I
     
       C  = A
!       IF(MSR)WRITE(*,*) 'ZBRENT=73, IMAX=',IMAX,'I=',I
     
       FA = FB
!       IF(MSR)WRITE(*,*) 'ZBRENT=74, IMAX=',IMAX,'I=',I
     
       FB = FC
!       IF(MSR)WRITE(*,*) 'ZBRENT=75, IMAX=',IMAX,'I=',I
     
       FC = FA
!       IF(MSR)WRITE(*,*) 'ZBRENT=76, IMAX=',IMAX,'I=',I
     
     ENDIF

     !
     !now we have |f(A)| greater than |f(B)|, which is either of the following two cases
     !
     !case 1:
     !                    
                                    !                  ~
                                    !                /  \
                                     !       ------A(C)-------B--------------->x
                                    !                /      \
                                    !               /       `\* f(B)
                                    !            ~ *              \
                                    !          /                `
                                    !            
                                    !               |<---------->|
                                    !                    E
                                    !
                                    !             |<--------->|
                                    !                  -XM*2
                                    !
     
     !case 2:
     
                                    !                  ~
                                    !                /  \
                                     !      -------(B)----(C)A---------------->x
                                    !               _^     \
                                    !         .- ~ *          `_
                                    !              f(B)        \
                                    !                        * f(A),f(C)
                                     !
                                    !               |<------->|
                                    !                  XM*2        !new window width
                                    !            
     
     
!     IF(MSR)WRITE(*,*) 'ZBRENT=5, IMAX=',IMAX,'I=',I
     TOL1 = 2.*EPS*ABS(B)+0.5*TOL    !update tolerance value based on size of B
     XM   = 0.5*(C-B)                !XM records half width of new window [B C]

     IF (ABS(XM) <= TOL1 .OR. FB == 0.) THEN        !if widnow is very samall or FB is already zero
                                                    !then B must be the solution to f(x)=0
       ZBOUT  = B
       
       !ZBOUT = (ZBOUT+SODold)/2.d0
       
!       IF(MSR)WRITE(*,*) 'ZBRENT=6, IMAX=',IMAX,'I=',I
       RETURN
     ENDIF
     
     
     IF (ABS(E) >= TOL1 .AND. ABS(FA) > ABS(FB)) THEN    !If window XM is still too large
                                                        !then check if previoius widnow E is greater 
                                                        !than tolreance and if f(A) is farther away 
                                                        !from zero than f(B)
                                                        !
                                                        !if yes, then calculate the slope of ratio FB/FA
                                                        !which should have |f(B)/f(A)|<1
     
       S = FB/FA
       
       IF (A == C) THEN                    !if A and C are the same then set P to scaled down size of window
         P = 2.*XM*S                    !by S, and have Q as ratio shortened
         Q = 1.-S
       ELSE
         Q = FA/FC                        !if A and C are not the same, e.g. C is with B, and f(C)*f(B)
                                        !
                                        !            ----A-------(C)B-----
                                        !
                                        !then have Q=f(A)/f(C), R = f(B)/f(C)
                                        !
                                        !
            
         R = FB/FC
         P = S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
         Q = (Q-1.)*(R-1.)*(S-1.)
       ENDIF
       
       !recalculate P and Q
       
       IF (P > 0.) Q = -Q
       P = ABS(P)
       
       IF (2.*P < MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
         E = D
         D = P/Q
       ELSE
         D = XM
         E = D
       ENDIF
     ELSE
       D = XM
       E = D
     ENDIF
     
!     IF(MSR)WRITE(*,*) 'ZBRENT=8, IMAX=',IMAX,'I=',I
     
     !set A and f(A) to B point
     
     A  = B
     FA = FB
     
     !
     !move B by half window size XM
     !and E still keeps a copy of windows size D
     !
     
     IF (ABS(D) > TOL1) THEN
       B = B+D
     ELSE
       B = B+SIGN(TOL1,XM)
     ENDIF
!     IF(MSR)WRITE(*,*) 'ZBRENT=9, IMAX=',IMAX,'I=',I
     
     !recalculate f(B)
     
     FB = SEDF(B,DF_SOD_TMP,SAV_SOD_TMP,TF_SSTATE)
     
!     IF(MSR)WRITE(*,*) 'ZBRENT=10, IMAX=',IMAX,'I=',I
     
   ENDDO
   
   !Oops, do loop of IMAX interation finished and still not converging
   !save the current results in ZBOUT
   IERR   = 2
   ZBOUT  = B
   !ZBOUT = (ZBOUT+SODold)/2.d0   !WLong, use gp method to average with previous step
   
!   IF(MSR)WRITE(*,*) 'ZBRENT=11, IERR=',IERR,'B=',B 

   RETURN
   END FUNCTION ZBRENT

!********************************************************************************
!**                    S U B R O U T I N E   S E D T S F N L                   **
!********************************************************************************
   !
   !Time series solution of sediment diagenesis
   !
   
   SUBROUTINE SEDTSFNL(C1S,C2S,CT1S,CT2S,C1TM1S,CT2TM1S,STEADYSTATE) 
   !inputs: C1TM1S         !total concentration at last time step
   !        CT2TM1S      !total concentration at last time step
   !
   !outputs: C1S        !dissolved concentration in layer 1
   !         C2S        !dissolved concentration in layer 2
   !         CT1S        !total concentration in layer 1
   !         CT2S        !total concentration in layer 2

   IMPLICIT NONE
   REAL(SP) :: C1TM1S,CT2TM1S
   LOGICAL :: STEADYSTATE
   
   REAL(SP) :: C1S,C2S,CT1S,CT2S
   REAL(SP) :: A11, A12, A21, A22, B_1, B_2
   REAL(SP) :: F12, F21, XK0, XK1, XK2, DELTA

!***** Initialize constants
     
   FD1 = 1./(1.+M1*PIE1)
   FP1 = M1*PIE1/(1.+M1*PIE1)
   FD2 = 1./(1.+M2*PIE2)
   FP2 = M2*PIE2/(1.+M2*PIE2)
   
   F12 = W12*FP1+KL12*FD1        !WLong (3-32-1-0)
   
   F21 = W12*FP2+KL12*FD2        !WLong (3-36-2)

!***** Evaluate the MM term at time level t-1

   IF (KMC1 /= 0.) THEN

     XK0 = (K0H1D*FD1+K0H1P*FP1)/(KMC1+C1TM1S)            !WLong (3-32-4)
                                                        !Note C1TM1S is dissolved fraction in layer 1
   ELSE
     XK0 = 0.                                            !WLong (3-32-3)
   ENDIF
   
   XK1 = XK0+K1H1D*FD1+K1H1P*FP1
   
   XK2 = K2H2D*FD2+K2H2P*FP2
   
  !Excel
  !If steadystate Then
  !  a22 = -fd2 * KL12 - fp2 * w12 - w2
  !  b2 = -Jn
  !   
  !Else  !unsteady state
  !  a22 = -fd2 * KL12 - fp2 * w12 - w2 -H2 / tc
  !  b2 = -Jn - H2 / tc * NH3Tp(2)
  !End If
  !
  
   A11 = -F12-XK1-W2                    !eqn(3-32)
   A21 = F12+W2                            !eqn(3-35)
   A12 = F21                            !eqn(3-33) 
   B_1 = -J1                            !eqn(3-34)
   
   IF(STEADYSTATE)THEN
        A22 = -F21-XK2-W2-K3            !eqn(3-36)
        B_2 = -J2                        !eqn(3-37)
   ELSE
        A22 = -F21-XK2-W2-K3-H2/DLTS    !eqn(3-36)
        B_2 = -J2-H2/DLTS*CT2TM1S        !eqn(3-37)
   ENDIF
  
!***** Solve the 2x2 set of linear equations

   DELTA = A11*A22-A12*A21
   IF (DELTA == 0.) THEN
     PRINT *,'DELTA is singular: A11,A12,A21,A22'
     write(*,*) f12, xk1, w2
     write(*,*) w12, fp1, kl12, fd1
     write(*,*) xk0, k1h1d, fd1, k1h1p, fp1
     PRINT *,A11,A12,A21,A22
     STOP
   ENDIF


!***** Assign results

   CT1S = (B_1*A22-B_2*A12)/DELTA    !Total concentration in layer 1
   CT2S = (B_2*A11-B_1*A21)/DELTA    !Total concentration in layer 2
   
   C1S  = FD1*CT1S                    !dissolved concentration in layer 1
   C2S  = FD2*CT2S                    !dissolved concentration in layer 2
   
 !  
 !!Do debugging here for PO4
 !  
 !  IF(PIE1==400.0) THEN !must be PO4 now
 !  
!        WRITE(*,*)'FD1', FD1
!        WRITE(*,*)'FD2', FD2
!        WRITE(*,*)'FP1', FP1
!        WRITE(*,*)'FP2', FP2
!        
!        WRITE(*,*)'B1, B2, DELTA, A11, A21', B1, B2, DELTA, A11, A21
!        
!        WRITE(*,*)'A12, A22',A12, A22
!        
!        WRITE(*,*)'F21',F21
!        
!        WRITE(*,*)'-F21-XK2-W2-K3, A22',-F21-XK2-W2-K3,A22
!        WRITE(*,*)'F21 XK2 W2 K3', F21,XK2, W2, K3  !XK2 =0 
!        WRITE(*,*)'XK1',XK1,'XK0+K1H1D*FD1+K1H1P*FP1',XK0+K1H1D*FD1+K1H1P*FP1
!        WRITE(*,*)'XK0, K1H1D, FD1, K1H1P, FP1', XK0,K1H1D,FD1,K1H1P,FP1
!        
!        
!        WRITE(*,*)'CT1S, CT2S', CT1S, CT2S
!        
!  ENDIF
!  

! 
!IF(MSR)THEN
!     WRITE(*,*)'KL12=', KL12
!        WRITE(*,*)'W12=', W12
!        WRITE(*,*)'W2=', W2
!        WRITE(*,*)'FD1=',FD1
!        WRITE(*,*)'FP1=',FP1
!        
!        WRITE(*,*)'FD2=',FD2
!        WRITE(*,*)'FP2=',FP2
!        
!        WRITE(*,*)'A11=',A11
!        WRITE(*,*)'A12=',A12
!        WRITE(*,*)'A21=',A21
!        WRITE(*,*)'A22=',A22
!        WRITE(*,*)'B_1=',B_1
!        WRITE(*,*)'B_2=',B_2
!!        READ(*,*)
!ENDIF
!   
   
   
   RETURN
   END SUBROUTINE SEDTSFNL

!********************************************************************************
!**                    S U B R O U T I N E   S E D S S F N L                   **
!********************************************************************************
   !steady state solution of sediment diagenesis
   SUBROUTINE SEDSSFNL(C1,C2,C2AV,CT1,CT2,CT2AV,C1TM1,CT2TM1,ITYPE) !TM1 ~ t-1
   IMPLICIT NONE
   INTEGER ITYPE
   REAL(SP) :: C1,C2,C2AV,CT1,CT2,CT2AV,C1TM1,CT2TM1

!  This subroutine translates between SEDTSFNL and SEDSSFNL

   IF(ITYPE.eq.1)THEN
      CALL SEDTSFNL (C1,C2,CT1,CT2,C1TM1,CT2TM1,.TRUE.)      !call time series solution with steadystate flag
      C2AV  = C2
      CT2AV = CT2
   ELSE
      STOP 'Oops, you are not supposed to call this function with ITYPE/=1'
   ENDIF

   RETURN
   END SUBROUTINE SEDSSFNL

   SUBROUTINE GET_SO4_HS_CH4_CSOD_FLUXES(TF_SSTATE)    
   !**** *******************************************************************************
   !**** New code for methane formation.  CH4 starts forming                           *
   !**** once all sulfate is used up. (Chapter 11 of DMD book instead of Chapter 10)   *
   !**** *******************************************************************************

   !
   !WLong: I moved this here (was coded in WASP and CE-QUAL-ICM but not
   !    verified if it works
   
   LOGICAL :: TF_SSTATE !Type of integration FALSE - time intergration from current step to next step
                        !                    TRUE  - steady state solution
!**** Sulfide and sulfate in O2 equivalents

   
   REAL(SP) ::     DCH4T2,            &    !rate of total  CH4 generation in layer 2 (gO2/m^3 * m/d  = gO2/m^2/d)
                DHST2                !rate of HS generation in layer 2 (gO2/m^3/d * m/d = gO2/m^2/d)

!   REAL(SP) :: SAVSOD_SED,        &
!               DFSOD_SED
   
! KURT GLAESEMANN Changed H to HH to avoid name conflict with module
   REAL(SP) AD(4,4), BX(4), G(2), HH(2,2)
   
   REAL(SP) DBLSO41,     &
            DBLSO42,     &
            RA0,         &
            RA1,         &
            RA2,         &
            R1,         &
            R2,         &
            DISC,         &
            SN1
            
   KMC1=0.0
   K0H1D=0.                                        !WLong (3-32-7)
   K0H1P=0.                                        !WLong (3-32-9)
   
   
   K1H1D=XAPP1HSD**2/S*(O20/(KMHSO2)) + S        !WLong notes 4-58 with O20 downstairs
                                                !WLong (3-32-5)
   K1H1P=XAPP1HSP**2/S*(O20/(KMHSO2))            !WLong (3-32-10)
   
   K2H2D=0.        !WLong (3-36-8)
   K2H2P=0.        !WLong (3-36-8-1)
   
   J1=0.        !WLong (4-52) J1 term with J1=S*fd0*HST0 + JH2ST,1, where HST0=0, JH2ST,1 is put in J2
                !instead
   
   K3=0.0        !WLong (3-36-6)
   
   J2=XJC1        !gO2/m2/day, carbon diagenesis in layer 2 in oxygen units after discounting diaggenesis
                !WLong (4-57) + (4-56) with JH2ST,1 and JH2ST,2 are put in JH2ST,2 together and ignoring 
                !layer 1 
   
   PIE1=PIE1HS    
   PIE2=PIE2HS

!**** Set KL12 using HSO4    

    ITEMP = 10.*TEMPD+1
   
!   DDSO4 = GET_ZL12NOM(TEMPD)*H2                !KL12=DD/H2 ==> DD=KL12*H2                
!   DDSO4 = GET_ZL12NOM(TEMPD)*DDP                 !Fixed by CFC as per JF 11/11/06 !WLong this CFC fix is WRONG
   
    DDSO4=KL12NOM*H2/2.0        !WLong fixed this
                                !  KL12NOM  = DDP/(H2/2)*GET_ZL12NOM(TEMPD)            !(m/d) (normalized KL12)
                                !==> DDP*GET_ZL12NOM(TEMPD)=(H2/2)*KL12NOM            !
                                !i.e. left hand is diffusion rate with temperature control, 
                                !right hand is (H2/2)*KL12NOM
    
#if(0)  /* Wen Long: not very clear on how the quadratic equation for SO42 is derived, based on WASP cand ICM code */
            
    
   !*
   !* Not sure how the quadratic equation for SO42 is derived
   !* overall it is a combination of HS and SO4 equations, with HS elimited to arrive two linear equations
   !* for SO41 and SO42, and then according to DMD book pp.232 eqn(11.27c) a third equation is needed for solving 
   !* the third unknown HSO4
   !*
   !* However the following statement calculates HSO4 directly based on DMD book pp.229 eqn (11.17)
   !* since HSO4 already calculated, then it is hard to understand why we need the third equation (11.27c) anymore
   !* and why we need to derive a quadratic equation for SO42 
   !*
   
!   HSO4  =SQRT(2.*DDSO4*SO40*H2/XJC1/)                !Equation (4-142) of WLong notes, where we seem to be missing "alpha_SO4,C"
                                                    !wihich is 1/2 *********
                                                    !Well, it depends on how DDSO4 is defined
                                                    !
   HSO4  =SQRT(2.*DDSO4*SO40*H2/(XJC1*(1./2.0)))    !WLong included 1/2.0 here which is alpha_SO4,C in eqn (4-142) 
                                                    !where XJC1 already included alpha_O2,C = 2.667gO2/gC
                                                    !i.e. WLong notes (4-8-18)
   
   
!**** No deeper than H2
   IF(HSO4 > H2) HSO4=H2
   
   KL12SO4=KL12*H2/HSO4        !KL12=DD/H2 ==> DD= KL12*H2 
                            !KL12SO4 = DD/HSO4
                            !==> KL12SO4 = KL12*H2/HSO4  (see WLong notes 4-143-1) (4-8-17)
 
    IDEBUG_SED=IDEBUG_SED+1
    !IF(MSR)THEN
    !        WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'K0H1D=',K0H1D
    !        WRITE(*,*)'K0H1P=',K0H1P
    !        WRITE(*,*)'K1H1D=',K1H1D
    !        WRITE(*,*)'K2H2D=',K2H2D
    !        WRITE(*,*)'K1H1P=',K1H1P
    !        WRITE(*,*)'K2H2P=',K2H2P
    !        WRITE(*,*)'J1=',J1
    !        WRITE(*,*)'K3=',K3
    !        WRITE(*,*)'J2=',J2
    !        WRITE(*,*)'PIE1=',PIE1
    !        WRITE(*,*)'PIE2=',PIE2
    !        WRITE(*,*)'DDSO4=',DDSO4
    !        WRITE(*,*)'HSO4=',HSO4
    !        WRITE(*,*)'KL12SO4=',KL12SO4
    !        READ(*,*)
    !ENDIF

!**** Fractions and overall decay reaction velocity
   FD1=1./(1.+M1*PIE1)
   FP1=M1*PIE1/(1.+M1*PIE1)
   FD2=1./(1.+M2*PIE2)
   FP2=M2*PIE2/(1.+M2*PIE2)
   FP1SO4=FP1            !WLong: ??? seems FP1SO4 and FP2SO4 should be zero *WLong notes (4-128))4-129)
   FP2SO4=FP2            !WLong: ??? seems FP1SO4 and FP2SO4 should be zero *WLong notes (4-128)(4-129)
   
   !
   !DMD book page 244 to 246
   !WLong notes (4-133), (4-144)
   !
   
   !equations for HS and SO4
   !
   !i.e.  AD(1,1)*SO41+ AD(1,2)*SO42+AD(1,3)*HST1 +BX(1) = 0  !layer 1 of SO4 eqn, see WLong (4-133)
   !                                                          !where HST1 is from (4-135)
   !      AD(2,1)*SO41+ AD(2,2)*SO42              +BX(2) = 0  !layer 2 of SO4 eqn, see WLong (4-134)
   !
   !       AD(3,3)*HST1+ AD(3,4)*HST2 +BX(3) = 0               !layer 1 of HS eqn, see WLong (4-48)
   !      AD(4,3)*HST1+ AD(4,4)*HST2 +BX(4) = 0               !layer 2 of HS eqn, see WLong (4-49)
   !
 
    !
   
   !
   !equations for SO4
   !

   AD(1,1) = -S-KL12SO4                                        !WLong eqn (4-8-7) assuming H1~=0, W2~=0
   AD(1,2) = KL12SO4                                                !WLong eqn (4-8-8) assuming H1~=0
   AD(2,1) = KL12SO4                                                !WLong eqn (4-8-10) H1~= constant (i.e. H1 ~=0) and W2=0

   IF(.NOT.TF_SSTATE)THEN
        AD(2,2) = -(DLTS*KL12SO4+H2)/DLTS            !WLong eqn (4-8-11) assuming H1~=0, W2=0 and H2,SO4 ~= H2
        !WLong used HSO4 rather than H2 according to WLong notes eqn (4-8-11) 
        !   AD(2,2) = -(DLTS*KL12SO4+HSO4)/DLTS        !WLong eqn (4-8-11) assuming H1~=0, W2=0 and H2,SO4 kept
   ELSE    
        AD(2,2) = -KL12SO4            !WLong eqn (4-8-11) assuming H1~=0, W2=0 and H2,SO4 ~= H2
   ENDIF
                                
   
   BX(1) = S*SO40                                    !Wlong (4-8-15) assuming H1~=0
   
   IF(.NOT.TF_SSTATE)THEN   
        BX(2) = H2*SO4T2TM1/DLTS                    !Wlong (4-8-12) assuming JSO4,2 = 0 and H2,SO4 ~= H2
        !BX(2) = HSO4*SO4T2TM1/DLTS                    !Wlong (4-8-12) assuming JSO4,2 = 0 and H2,SO4 ~= H2   
   ELSE
        BX(2) = 0.d0                                !Wlong (4-8-12) assuming JSO4,2 = 0        
        !WLong: should consider JSO4,2
   ENDIF
   
   KHS_1=     FP1*XAPP1HSP**2/S*(O20/(KMHSO2))    &    !WLong eqn (4-135), (4,48),(4-58)
            +FD1*XAPP1HSD**2/S*(O20/(KMHSO2))
            
   AD(1,3) = KHS_1    !Basically diagenesis of SO4 (oxidation of HS) in layer 1, i.e. JSO4,1 at t+dt time
   

   !
   !equations for HS 
   !
   
   
   !(after assuming H1~=0, JH2ST,2=0, JH2ST,1=0 (***WLong: ok to assume JH2ST,1=0, but 
   !should not assume JH2ST,2=0???), i.e. BX(4) should  have -J2 in it. However, here it is added to RA1
   !term in the solution for SO42)
   !

   AD(3,3) = -(W2)-(FP1)*(W12)-(FD1)*(S)         &
             -(FD1)*(KL12SO4)-(KHS_1)                            
!             -(FD1)*(KL12)-(KHS_1)                            !WLong, *****should be KL12 *****   
                                                            !WLong notes eqn (4-8-19)
             
   AD(3,4) = (FP2)*(W12)+(FD2)*(KL12SO4)        
!   AD(3,4) = (FP2)*(W12)+(FD2)*(KL12)                        !WLong, *****should be KL12 ****   
                                                            !WLong notes eqn (4-8-20)
   
   AD(4,3) = (W2)+(FP1)*(W12)+(FD1)*(KL12SO4)            
!   AD(4,3) = (W2)+(FP1)*(W12)+(FD1)*(KL12)                    !WLong, **** should be KL12*****
                                                            !WLong eqn (4-8-22)
   
   IF(.NOT.TF_SSTATE)THEN
        AD(4,4) = -((DLTS)*(FP2)*(W12)                        &!!WLong eqn (4-8-23)
                    +(DLTS)*(FD2)*(KL12SO4)                      &!
    !              +(DLTS)*(FD2)*(KL12)                          &!!WLong, ****should be KL12 ****          
                    +(DLTS)*(W2)                               &!
                    +(H2)                                      &!
                    )/(DLTS)
   ELSE
        AD(4,4) = -((FP2)*(W12)                       &!WLong eqn (4-8-23)
                    +(FD2)*(KL12SO4)                  &!
    !               +(FD2)*(KL12)                      &!WLong, ****should be KL12 ****          
                    +(W2)                               &!
                    )
   ENDIF

   BX(3) = (HS0)*(S)                                !WLong eqn (4-8-21)
   
   IF(.NOT.TF_SSTATE)THEN
        BX(4) =  (H2)*(HST2TM1)/(DLTS)    !&!            !WLong eqn (4-8-24)
        !           +J2                                !**** should not really ignore JH2ST,2****
        
   ELSE
        BX(4) =  0.d0         !&!!WLong eqn (4-8-24)
        !           +J2        !**** should not really ignore JH2ST,2****
   ENDIF
   
!   
!maybe the J2 is implemented only when SO42 is solved with J2, see RA1 term below
!  

   !
   !re-arrange the equations for easier solution
   !
   !    SO41*HH(1,1)+HH(1,2)*SO42 + G(1) = 0
   !    SO41*HH(2,1)+HH(2,2)*SO42 + G(2) =0
   !           
   !find solution of HST1 in terms of AD(3,3), AD(3,4), AD(4,3),AD(4,4),BX(3),BX(4)
   !then substitute that into equations of SO4
   !
   !See WLong notes eqn(4-8-29) and (4-8-27)
           
   G(1) = ( BX(1)*AD(3,3)*AD(4,4)             &
          - BX(3)*AD(1,3)*AD(4,4)            &
          - BX(1)*AD(3,4)*AD(4,3)             &
          + BX(4)*AD(1,3)*AD(3,4)             &
          + BX(2)*AD(1,3)*AD(3,4)            &
          )/(AD(1,3)*AD(3,4))

   !
   !i.e. AD(1,3)*AD(3,4)*G(1) =  AD(3,3)*AD(4,4)                    *BX(1)
   !                            -AD(3,4)*AD(4,3)                    *BX(1)
   !                            +AD(1,3)*AD(3,4)                    *BX(2) 
   !                            -AD(1,3)*AD(4,4)                    *BX(3)
   !                            +AD(1,3)*AD(3,4)                    *BX(4)
   !

   
   G(2) = (   BX(1)*AD(3,3)*AD(4,4)         &
            - BX(3)*AD(1,3)*AD(4,4)         &
            - BX(1)*AD(3,4)*AD(4,3)         &
            + BX(4)*AD(1,3)*AD(3,4)            &
          )/(AD(1,3)*AD(3,4))

   !
   !i.e. AD(1,3)*AD(3,4)*G(2) =    AD(3,3)*AD(4,4)                    * BX(1) 
   !                             - AD(3,4)*AD(4,3)                    * BX(1)
   !                             - AD(1,3)*AD(4,4)                    * BX(3)
   !                              +AD(1,3)*AD(3,4)                    * BX(4)

   
   
   HH(1,1) = (     AD(1,1)*AD(3,3)*AD(4,4)    &
                -AD(1,1)*AD(3,4)*AD(4,3)    &
                +AD(1,3)*AD(2,1)*AD(3,4)    &
             )/(AD(1,3)*AD(3,4))
   
   HH(1,2) = ( AD(1,2)*AD(3,3)*AD(4,4)        &
              -AD(1,2)*AD(3,4)*AD(4,3)        &
              +AD(1,3)*AD(2,2)*AD(3,4)        &
              )/(AD(1,3)*AD(3,4))
   
   HH(2,1) = ( AD(1,1)*AD(3,3)*AD(4,4)        &
              -AD(1,1)*AD(3,4)*AD(4,3)        &
              )/(AD(1,3)*AD(3,4))
   
   HH(2,2) = ( AD(1,2)*AD(3,3)*AD(4,4)        &
              -AD(1,2)*AD(3,4)*AD(4,3)        &
              )/(AD(1,3)*AD(3,4))

   
   !WLong: the above are still only linear equations to solve, NOT sure where does the 
   !nonlinear (quadratic) equation for SO41 come from. And it is related to KMSO4
   !which is saturation concentration of SO4 for generation of CH4.
   !ie. presense of SO4 in layer 2 limits production of CH4
   
   !
   !    SO41*HH(1,1)+HH(1,2)*SO42 + G(1) = 0
   !    SO41*HH(2,1)+HH(2,2)*SO42 + G(2) =0
   !
   
   
   !WLong: I guess this is same as DMD book page 245 except it is written as a quadratic equatin of 
   !SO42 instead of HSO4. (Need to re-derive this from (4-152) of WLong notes)
   !
   !The derivation of RA0, RA1, RA2 should be similar to the lines of DMD book page 248, except
   !the book only listed them for solving steady state CH4 equations, the derivation was done
   !using MAXIMA software
   !
   !W. Long did simiar derivations for SO42 by hand, see notes (4-8-63) to (4-8-84), except
   !that AD(1,3), i.e. A13 (4-8-16) is assumed zero (JSO4,1=0, see comments above (4-8-50) in WLong notes ). 
   !The deriviation by Wen Long had J2 (JSO4,2) included (4-8-55), (4-8-12), however here BX(2) has J2 ignored
   !Also BX(4) has J2 ignored as well! Following W.Long's hand derivation we see that in (4-8-82), i.e. RA1
   !has J2 in it. Below, RA1 does have J2 included, so alrough BX(2) and BX(4) above did not contain J2
   !somehow it is included below in RA1 term.
   !
   
   !
   !Find the solution for SO42 based on quadratic equation
   !
   !    RA2 *  x^2 + RA1 * x + RA0 = 0
   !
   !where x is SO42 concentration (mgO2/L) in oxygen equivalents
   !
   !     -RA1 +/- sqrt( RA1^2 - 4.0*RA0* RA2)       DISC
   ! x= --------------------------------------- = ---------
   !           2 * RA2                               RA2
   !
   !                -[RA1 -/+ sqrt(RA1^2 - 4.0 * RA0 * RA2)]
   !where DISC = -------------------------------------------
   !                            2
   !
   !The solution method is based on Vieta's formula: (http://en.wikipedia.org/wiki/Quadratic_equation)
   !
   !  x1+x2= - RA1/RA2    (i.e. x1+x2 = -b/a for equation a*x^2+b*x+c=0)
   !  x1*x2=   RA0/RA2    (i.e. x1*x2 =  c/a for equation a*x^2+b*x+c=0)
   !
   
   RA2 = HH(1,1)*HH(2,2)-HH(1,2)*HH(2,1)
   
   RA0 = (  HH(1,1)*G(2)                                    &
           -HH(2,1)*G(1)                                    &
          )*(KMSO4)        !Here KMSO4 is saturation concentration of SO4 in oxygen equivalents
                            !(mgO2/L) for methane generation, DMD book does not have this term
                            !
   
   RA1 = ( HH(1,1)*G(2)                                     &
         - HH(2,1)*G(1)                                     &
         )                                                    &   
         + (  HH(1,1)*HH(2,2)                                &
             -HH(1,2)*HH(2,1)                                &
           )*(KMSO4)                                     &         
         + HH(1,1)*J2          !WLong, why add J2 here ???? This shouuld be in the JHST,2 in BX(4) above
         
   !     
   !i.e.   (H11*H22-H12*H21)               *x^2 
   !     + [ (H11*G2-H21*G1) + 
   !        +(H11*H22-H12*H21)KMSO4 
   !        + H11*J2 
   !       ]                               *x 
   !     + (H11*G2-H21*G1)*KMSO4 = 0
   !

   !WLong: then solve this equation: a2*x^2+a1*x+a0=0
  
   !Sign infront of sqrt
   SN1 = 1.                          !solution of a2*x^2+a1*x+a0=0

   IF (RA1 <= 0.0) SN1 = -1.         !see Num Rec p178
   
   DISC = -( RA1                                        &
            +SN1*DSQRT(RA1**2-4.0D0*RA2*RA0)            &
            )/2.0D0
   
   IF (ABS(DISC) /= 0.0D0) THEN  
   
     R1 = DISC / RA2            !first root
                                !==> DISC= R1*RA2
     R2 = RA0 / DISC            !==> basically x2 = c/(a*x1)
    
     
   ELSE   ! vjp 11/17/2005 added logic for case (a2*a0) = 0 
     IF (ABS(RA2) == 0.0D0) THEN   !  a2 = 0, a1/=0     
     
        !RA1*x + RA0 = 0 ==> x = - RA0/RA1
     
       R1 = -RA0/RA1                !first root
       R2 = R1                        !2nd root
       
     ELSE                           !  a0 = 0, a2 /= 0
     
        !(RA2*x+RA1)x = 0
        !==> x = 0 or x= -RA1/RA2
       R1 = -RA1/RA2                !first root
       R2 = 0.0D0                    !2nd root
     ENDIF
   ENDIF

   DBLSO42 = R1

   IF (DBLSO42 < 0.) DBLSO42 = R2
   

   !solve for SO41 from 
   !
   !  SO41*HH(1,1)+HH(1,2)*SO42 + G(1) = 0
   !
   
   DBLSO41 = -(HH(1,2)*DBLSO42+G(1))/HH(1,1)
   
!   HST1 = -( AD(1,2)*DBLSO42                &        !James Martin code
!             +AD(1,1)*DBLSO41+BX(1)        &
!             )/AD(1,3)

   HST1 = -( AD(1,2)*DBLSO42                &
            +AD(1,1)*DBLSO41                &
            +BX(1)                            &
            )/AD(1,3)
   !
   !i.e.  AD(1,1)*SO41+ AD(1,2)*SO42+AD(1,3)*HST1 +BX(1) = 0  !layer 1 of SO4 eqn, see WLong (4-133)
   !                                                          !where HST1 is from (4-135)
   !      AD(2,1)*SO41+ AD(2,2)*SO42+             +BX(2) = 0  !layer 2 of SO4 eqn, see WLong (4-134)
   

   HST2 =     (    AD(1,2)*AD(3,3)*DBLSO42        &    !James Martin code
                 +AD(1,1)*AD(3,3)*DBLSO41        &
               +BX(1)*AD(3,3)                &
               -AD(1,3)*BX(3)                &
             )/(AD(1,3)*AD(3,4))
   
!   HST2 = ( (     AD(1,2)*DBLSO42            &
!                +AD(1,1)*DBLSO41            &
!                +BX(1))*AD(3,3)                &
!           -AD(1,3)*BX(3)                    &
!           )/(AD(1,3)*AD(3,4))

   !
   !i.e. HST2 = (-HST1 * AD(3,3)/AD(3,4) - BX(3)/AD(3,4)            
   !
   !i.e.          -HST1*AD(3,3)-BX(3) 
   !     HST2 = ---------------------
   !                   AD(3,4)
   !

   !i.e. 
   !
   !    AD(3,3)*HST1 + AD(3,4)*HST2  +BX(3) = 0 (layer 1 of HS eqn)  !see WLong notes
   !    AD(4,3)*HST1 + AD(4,4)*HST2  +BX(4) = 0 (layer 2 of HS eqn)
   
   
   !
   !Check here and make sure HST1 and HST2 do not go negative ***WLONG***
   !    

!   IF(HST1.LT.0.0)HST1=0
!   IF(HST2.LT.0.0)HST2=0
   
#else /* Wen Long included the following code for calculating HSO4, SO42,SO41, HS1, HS2, Chapter 11 of DMD */
   
    !
    !The solution procedures follows DMD book page 231 to 232 equations (11.27a) ~(11.27c)
    !
    !Detailed derivations are based on Wen Long notes (4-8-1) to (4-8-82)
    !
    !First solve 3 unknowns: SO41, SO42, HSO4
    !
    !Then solve 2 unknowns: HST1, HST2
    !
   
   !
   !equations for solving SO41, SO42, HSO4
   !
   
   !
   !A11*SO41 + A12 * SO42 + B1PRIME=0                            (4-8-14-1) with A13 ~=0
   !A21*SO41 + A22 * SO42 + B2=0                                    (4-8-14-2)
   !
   !
   !A11     = -S-KL12SO4-W2-{H1/DLTS} - {(DH1DT+DH1DTM)}        (4-8-7)
   !A12     = KL12SO4 + {DH1DTP    }                                (4-8-8)
   !B1PRIME    = S*SO40+JSO41+{H1*SO41TM1/DLTS}                    (4-8-9)
   !
   !A21     = KL12SO4+W2 + {DH1DTM}                                (4-8-10)
   !A22        =-KL12SO4-W2-{(DH1DTP-DH1DT)} - {HSO4/DLTS}            (4-8-11)
   !B2        = JSO42 + {HSO4*SO42TM1/DLTS}                        (4-8-12)
   !
   !where {} terms drop if assuming steady state
   !
   !Assume steady steate and W2=0 and JSO4,1=0 (i.e. JH2ST,1=0)
   !A11 = -S - DDSO4/HSO4                                (4-8-50)
   !A12 = DDSO4/HSO4                                    (4-8-51)
   !B1  =BPRIME = S * SO40                                (4-8-52)
   !        
   !A21= DDSO4/HSO4                                        (4-8-53)
   !A22= -DDSO4/HSO4                                    (4-8-54)
   !B2 = -alpha_SO4C * alpha_O2C * JC2 * HSO4/H2        (4-8-55)
   !

   !
   !The third equation: 
   !
   ! SO42=SO41/2.0                     (4-8-65)
   !
   !
   !putting (4-8-14-1) (4-8-14-2) (4-8-65) together, we can eliminate SO42,SO41 and obtain 
   !the following quadratic equation for HSO4:
   !
   !a_HSO4* HSO4**2 + b_HSO4 *HSO4 + c_HSO4 = 0        (4-8-71)
   !
   !where 
   !
   !a_HSO4 = RA2 = 2*S*alpha_SO4C*alpha_O2C*JC2        (4-8-72)
   !b_HSO4 = RA1 = DDSO4*alpha_SO4C*alpha_O2C*JC2    (4-8-73)
   !c_HSO4 = RA0 = -S*H2*DDSO4*SO40                    (4-8-74)
   !
   !
   !==> HSO4 =[ -b_HSO4 +/- sqrt(b_HSO4^2-4*a_HSO4*c_HSO4) ]/(2*a_HSO4)  (4-8-75)
   !
   !

   !
   !equations for solving HST1, HST2
   !
   ! A33 * HST1 + A34*HST2 + B3 = 0                                                    (4-8-3)
   ! A43 * HST1 + A4$*HST2 + B4 = 0                                                    (4-8-4)
   !
   ! A33 = -S*FD1-W12*FP1-KL12*FD1-W2-JH2ST1 - {H1/DLTS} - {(DH1DTM + DH1DT)}        (4-8-19)
   ! A34 =  W12*FP2 + KL12*FD2+ {DH1DTP}                                            (4-8-20)
   ! B3  = S*HS0 + JH2ST1 + {H1*HST1TM1/DLTS}                                        (4-8-21)
   !
   ! A43 = W12*FP1+KL12*FD1+W2+{DH1DTM}                                                (4-8-22)
   ! A44 = -W2-KL12*FD2-W12*FP2-{H2/DLTS} - {(DH1DTP-DH1DT)}                        (4-8-23)
   ! B4  = JH2ST2 + {H2*HST2TM1/DLTS}                                                (4-8-24)
   !
   !Assuming steady state,the {} terms in above coefficients are dropped (zero)
   !

   !
   !JH2ST1 = - JSO41 = alpha_O2C*JPOC1-alpha_O2C*alpha_CN*Kappa_NO31**2*THETANO3**(T-20)/S (4-56) (4-8-13)
   !
   !JH2ST2 = - JSO42 = alpha_O2C*JPOC2-alpha_O2C*alpha_CN*Kappa_NO32**2*THETANO3**(T-20)/S (4-57)
   !
   !The above JH1ST1, JH2ST2 calculation assumes carbon generated by diagenesis (after discounting denitrification) 
   !will be all oxidized SO4 and generat H2S in layer 2 and will be eventually oxidized by O2 in layer1 
   !
   !For convenience and under the assumption that H1 << H2, the  JH2ST1 and JH2ST2 are combined into JH2ST2
   !such that
   !
   !JH2ST1=-JSO41 = 0
   !JH2ST2=-JSO42 = (HSO4/H2)* {alpha_O2C*(JPOC1+JPOC2) - [ alpha_O2C*alpha_CN*Kappa_NO31**2*NO31/S
   !                                    +alpha_O2C*alpha_CN*Kappa_NO32**2*NO32/S]*THETANO3**(T-20)}     
   !                                                                                (4-8-55)
   !if HSO4< H2, there will be a portion left for methane genereation, i.e. the portion
   !oxidized is assumed to be HSO4/H2
   !
   
   
!   HSO4  =SQRT(2.*DDSO4*SO40*H2/XJC1/)        !Equation (4-142) of WLong notes, where we seem to be missing "alpha_SO4,C"
                                            !wihich is 1/2 *********
                                            !Well, it depends on how DDSO4 is defined
                                            !
!   HSO4  =SQRT(2.*DDSO4*SO40*H2/(XJC1*(1./2.0)))    !WLong included 1/2.0 here which is alpha_SO4,C in eqn (4-142) 
                                                    !where XJC1 already included alpha_O2,C = 2.667gO2/gC
                                                    !i.e. WLong notes (4-8-18)
   
    IDEBUG_SED=IDEBUG_SED+1
    !IF(MSR)THEN
    !        WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'K0H1D=',K0H1D
    !        WRITE(*,*)'K0H1P=',K0H1P
    !        WRITE(*,*)'K1H1D=',K1H1D
    !        WRITE(*,*)'K2H2D=',K2H2D
    !        WRITE(*,*)'K1H1P=',K1H1P
    !        WRITE(*,*)'K2H2P=',K2H2P
    !        WRITE(*,*)'J1=',J1
    !        WRITE(*,*)'K3=',K3
    !        WRITE(*,*)'J2=',J2
    !        WRITE(*,*)'PIE1=',PIE1
    !        WRITE(*,*)'PIE2=',PIE2
    !        WRITE(*,*)'DDSO4=',DDSO4
    !        WRITE(*,*)'HSO4=',HSO4
    !        WRITE(*,*)'KL12SO4=',KL12SO4
    !        READ(*,*)
    !ENDIF

!**** Fractions and overall decay reaction velocity
    FD1=1./(1.+M1*PIE1)
    FP1=M1*PIE1/(1.+M1*PIE1)
    FD2=1./(1.+M2*PIE2)
    FP2=M2*PIE2/(1.+M2*PIE2)
    FP1SO4=FP1
    FP2SO4=FP2
   
   
   !
   !DMD book page 244 to 246
   !WLong notes (4-133), (4-144)
   !
    
   !
   !equation for HSO4 (WLong (4-8-71), assuming W2=0, steady steate, JSO41=0
   !
   !     RA2*HSO4**2 + RA1*HSO4 + RA0 = 0        (quadratic equation for HSO4)
   !
   !then solve for HSO4
   
   RA2=2*0.5*S*J2            !here J2 is in gO2/m^2/d and has carbon diagenesis discounted by denitrification
                            !what's left is to be oxidized by sulfate or for generating methane
                            !0.5 is alphas_SO4_C in oxygen equivalents (gO2/gO2)
                            !
   RA1=DDSO4*0.5*J2
   
   RA0=-S*H2*DDSO4*SO40
   
   !IDEBUG_SED=IDEBUG_SED+1
   !IF(MSR)THEN        
!            WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
!            WRITE(*,*)'RA2=',RA2
!            WRITE(*,*)'RA1=',RA1
!            WRITE(*,*)'RA0=',RA0
!            !READ(*,*)
   !ENDIF       
   
      !Sign in front of sqrt
   SN1 = 1.                          !solution of a2*x^2+a1*x+a0=0

   IF (RA1 <= 0.0) SN1 = -1.         !see Num Rec p178
   
   DISC = -( RA1                                        &
            +SN1*DSQRT(RA1**2-4.0D0*RA2*RA0)            &
            )/2.0D0
   
   IF (ABS(DISC) /= 0.0D0) THEN   ! B Clark changed from DABS to ABS for single precesion compilation
   
     R1 = DISC / RA2            !first root
                                !==> DISC= R1*RA2
     R2 = RA0 / DISC            !==> basically x2 = c/(a*x1)
    
     
   ELSE   !logic for case (a2*a0) = 0
     IF (ABS(RA2) == 0.0D0) THEN   !  a2 = 0, a1/=0       ! B Clark changed from DABS to ABS for single precesion compilation
     
        !RA1*x + RA0 = 0 ==> x = - RA0/RA1
     
       R1 = -RA0/RA1                !first root
       R2 = R1                        !2nd root
       
     ELSE                           !  a0 = 0, a2 /= 0
     
        !(RA2*x+RA1)x = 0
        !==> x = 0 or x= -RA1/RA2
       R1 = -RA1/RA2                !first root
       R2 = 0.0D0                    !2nd root
     ENDIF
   ENDIF

    !take the positive root as HSO4
   
      HSO4 = R1
      IF (HSO4 < 0.0) HSO4 = R2
      
      IF(HSO4<0.0)THEN
        WRITE(*,*)'oops, HSO4 is less than zero!!'
        HSO4=0.0000001  !make a very small number to avoid singularity
      ENDIF
      
      IF(HSO4>H2) HSO4=H2
    
   !
   !equations for SO4 (steady state, JSO41=0, H1~=0)
   !
      
   KL12SO4=DDSO4/HSO4        !KL12=DD/H2 ==> DD= KL12*H2 
                            !KL12SO4 = DD/HSO4
                            !==> KL12SO4 = KL12*H2/HSO4  (see WLong notes 4-143-1) (4-8-17)
 
   !A11*SO41 + A12 * SO42 + B1PRIME=0                    (4-8-14-1) with A13 ~=0
   !A21*SO41 + A22 * SO42 + B2=0                            (4-8-14-2)

   !A11 = -S - DDSO4/HSO4                                (4-8-50)
   !A12 = DDSO4/HSO4                                    (4-8-51)
   !B1  = S * SO40                                        (4-8-52)
   !        
   !A21= DDSO4/HSO4                                        (4-8-53)
   !A22= -DDSO4/HSO4                                    (4-8-54)
   !B2 = -alpha_SO4C * alpha_O2C * JC2 * HSO4/H2        (4-8-55)
   !

   
   AD(1,1)     = -(S)-(KL12SO4)                                    !(4-8-7) assuming H1~=0, W2~=0
   AD(1,2)     = (KL12SO4)                                                !(4-8-8) assuming H1~=0
   BX(1)     = (S)*(SO40)                                        !(4-8-15) assuming H1~=0
   
   AD(2,1)     = (KL12SO4)                                                !(4-8-10) H1~= constant (i.e. H1 ~=0) and W2=0

   AD(2,2)     = -(KL12SO4) !-(HSO4))/(DLTS)                    !(4-8-11) assuming H1~=0, W2=0
   
   
!  BX(2)     = (H2)*(SO4T2TM1)/(DLTS)                        !(4-8-12) assuming JSO4,2 = 0
   BX(2)     = -J2*0.5 !+(HSO4)*(SO4T2TM1)/(DLTS)            !(4-8-12) assuming steady state

            !here J2 is JH2ST in gO2/m^2/d and 0.5 is alpha_SO4C in gO2/gO2

   
   !solve linear equation for SO42 and SO41 (based on (3-30)(3-31) of WLong notes)
   
    DBLSO41=   (-AD(2,2)*BX(1)+AD(1,2)*BX(2))/                &
            (AD(1,1)*AD(2,2)-AD(1,2)*AD(2,1))
    
    
   
    DBLSO42=   (-AD(1,1)*BX(2)+AD(2,1)*BX(1))/                &
            (AD(1,1)*AD(2,2)-AD(1,2)*AD(2,1))
            

    IF(DBLSO41<0.0)DBLSO41=0.0        
    IF(DBLSO42<0.0)DBLSO42=0.0

   !check if SO42= SO41/2.0 is ture
    
    IF(DABS(DBLSO41-DBLSO42*2.0)>0.00001)THEN
        WRITE(*,*)'Oops, SO41 /= SO42*2 '
        WRITE(*,*)'SO41=', DBLSO41, 'SO42*2=',DBLSO42*2
    ENDIF
        
   !equation for HST1, HST2 (total HS in layer 1 and layer 2 (gO2/m^3))
    

   ! A33 * HST1 + A34*HST2 + B3 = 0                                                    (4-8-3)
   ! A43 * HST1 + A4$*HST2 + B4 = 0                                                    (4-8-4)
   !
   ! A33 = -S*FD1-W12*FP1-KL12*FD1-W2-JH2ST1 - {H1/DLTS} - {(DH1DTM + DH1DT)}        (4-8-19)
   ! A34 =  W12*FP2 + KL12*FD2+ {DH1DTP}                                            (4-8-20)
   ! B3  = S*HS0 + JH2ST1 + {H1*HST1TM1/DLTS}                                        (4-8-21)
   !
   ! A43 = W12*FP1+KL12*FD1+W2+{DH1DTM}                                                (4-8-22)
   ! A44 = -W2-KL12*FD2-W12*FP2-{H2/DLTS} - {(DH1DTP-DH1DT)}                        (4-8-23)
   ! B4  = JH2ST2 + {H2*HST2TM1/DLTS}                                                (4-8-24)
   !
    
   KHS_1=     FP1*XAPP1HSP**2/S*(O20/(KMHSO2))                &        !WLong eqn (4-135), (4,48),(4-58)
            +FD1*XAPP1HSD**2/S*(O20/(KMHSO2))
    
  !Assuming H1~=0, JH2ST,1=0.0

   AD(3,3) = -(W2)-(FP1)*(W12)-(FD1)*(S)         &
!             -(FD1)*(KL12SO4)-(KHS_1)                        
             -(FD1)*(KL12)-(KHS_1)                                !(4-8-19)
             
!  AD(3,4) = (FP2)*(W12)+(FD2)*(KL12SO4)        
   AD(3,4) = (FP2)*(W12)+(FD2)*(KL12)                        !(4-8-20)
   
!  AD(4,3) = (W2)+(FP1)*(W12)+(FD1)*(KL12SO4)            
   AD(4,3) = (W2)+(FP1)*(W12)+(FD1)*(KL12)                !(4-8-22)
   
   IF(.NOT.TF_SSTATE)THEN
        AD(4,4) = -((DLTS)*(FP2)*(W12)                   &!(4-8-23)
!                   +(DLTS)*(FD2)*(KL12SO4)                 &!
                    +(DLTS)*(FD2)*(KL12)                 &!
                    +(DLTS)*(W2)                          &!
                    +(H2)                                 &!
                    )/(DLTS)
   
   ELSE
           AD(4,4) = -((FP2)*(W12)                          &!        !(4-8-23)
!                   +(FD2)*(KL12SO4)                     &!
                    +(FD2)*(KL12)                         &!
                    +(W2)                                  &!
                    )
   ENDIF

   BX(3) =  (HS0)*(S)                                    !(4-8-21)
   
   IF(.NOT.TF_SSTATE)THEN
        BX(4) =  (H2)*(HST2TM1)/(DLTS)+J2*0.5            !(4-8-24)
        !here J2 is JH2ST,2 in gO2/m^2/d, 0.5 is to convert to SO4 unit also in gO2 (4-57)
   ELSE
        BX(4) =  +J2*0.5                                !(4-8-24)
   ENDIF
   
   !
   !solve for HST1, HST2 according to (3-30)(3-31) of WLong's notes
   !
   
    HST1=    (-AD(4,4)*BX(3)+AD(3,4)*BX(4))/                &
            (AD(3,3)*AD(4,4)-AD(3,4)*AD(4,3))
    
    HST2=   (-AD(3,3)*BX(4)+AD(4,3)*BX(3))/                &
            (AD(3,3)*AD(4,4)-AD(3,4)*AD(4,3))
            
   !check to make sure both HST1 HST2 >0
   IF(HST1< 0.0)HST1=0
   IF(HST2< 0.0)HST2=0   
   
#endif    /*end solving SO4 and HSO4 and HS */
   
   HS1=FD1*HST1
   HS2=FD2*HST2
   HS2AV=FD2*HST2
   
   IDEBUG_SED=IDEBUG_SED+1
   !IF(MSR)THEN        
!            WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
!            WRITE(*,*)'HST1=',HST1
!            WRITE(*,*)'HST2=',HST2            
!            WRITE(*,*)'DBLSO42=',DBLSO42
!            WRITE(*,*)'DBLSO41=',DBLSO41
!            WRITE(*,*)'HS1=',HS1
!            WRITE(*,*)'HS2=',HS2
!            WRITE(*,*)'HS2AV=',HS2AV
!            WRITE(*,*)'HSO4=',HSO4
!            READ(*,*)
   !ENDIF    

   SO41        =    DBLSO41
   SO42        =    DBLSO42
   SO42AV    =    SO42
   SO4T2     =    SO42

   JHS        =    S*(HS1-HS0) !gO2/m^2/day (flux of dissolved HS to overlying water column)
                            !WLong (4-61)

  !CSOD due to HS oxidation in layer 1 limited by avaiable O2 in overlying water
  ! (m/d)*(gO2/m^3)*(1) ==> unit is gO2/m^2/d
   
   CSODHS= (  XAPP1HSP**2/S*FD1*HST1                 & !oxidation of dissolved HS in layer 1
            + XAPP1HSP**2/S*FP1*HST1                & !oxidation of particulate HS in layer 1
           )*(O20/KMHSO2)                               !limitation by avaiable O20
                                    !WLong (4-62)

    IDEBUG_SED=IDEBUG_SED+1
    !IF(MSR)THEN
    !        WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'SO41=',SO41            
    !        WRITE(*,*)'SO42=',SO42
    !        WRITE(*,*)'SO42AV=',SO42AV            
    !        WRITE(*,*)'SO4T2=',SO4T2
    !        WRITE(*,*)'JHS=',JHS
    !        WRITE(*,*)'HS1=',HS1
    !        WRITE(*,*)'HS0=',HS0
    !        WRITE(*,*)'CSODHS=',CSODHS
    !        !READ(*,*)
    !ENDIF    

!
!solve for methane equation
!
     
    XJ2        =    J2*KMSO4/(SO42+KMSO4)  !WLong (4-70) (4-157-2)
                                       !Where J2 = XJC1 = JCX-XJCNO3
   
    !gO2/m^2/day !generation of CH4 in layer 2 regulated by abundance level of SO4 in layer 2
    !If there is plenty of SO4 then there will be no production of CH4
    !as SO4 will be used to oxidize the organic matter CHO and generate CO2 directly
   
    !
    !WLong: the use of SO4 in layer 2 for oxidizing carbon should be limited by HSO4 rather than
    !SO42, hence I'm not so convinced here on the use of KMSO4, unless SO42 is SO42*HSO4/H2
    !
   
   XJ2CH4    =    XJ2    !WLong (4-70) flux of methane generation (gO2/m^2/d)
   
   X1J2        =    J2*DBLSO42/(SO42+KMSO4)     !never used = J2-XJ2
                                            !e.g. the amount of J2 oxidized by SO42

     !WLong, I think we should use J2-J2*HSO4/H2 as the diagenesis for JCH4 (4-70) (i.e. 4-157))
   
   
    IDEBUG_SED=IDEBUG_SED+1
    !IF(MSR)THEN
    !        WRITE(*,*)'JDAY=',JDAY,'IDEBUG_SED=',IDEBUG_SED
    !        WRITE(*,*)'XJ2=',XJ2
    !        WRITE(*,*)'XJ2CH4=',XJ2CH4
    !        READ(*,*)
    !ENDIF    
      
!**** Methane solved assuming steady state

   KMC1=0.0
   PIE1=0.0
   PIE2=0.0
   
   K0H1D=0.        !WLong (3-32-7)
   K0H1P=0.     !WLong (3-32-9)
   
   K1H1P=0.     !WLong (3-32-10, which should be non-zero but equivalent to zero with fp1 =0 in (3-32-11)
   K1H1D=XAPPCH4**2/S*(O20/(KMCH4O2+O20))+S   !WLong (4-162-1)(3-32-11 with fp1=0) (3-32-5)
   
   K2H2D=0.        !WLong (3-36-8)
   K2H2P=0.     !WLong (3-36-8-1)
   
   J1=S*CH40    !WLong (4-164)

   K3=0.0       !WLong (3-36-6)
   
   J2=XJ2        !gO2/m^2/day  !carbon diagenesis flux converted from mgC/m^3/day to gO2/m^3/day
                !WLong: should be X1J2 instead?
                !WLong (4-167) (4-70) (4-71)
   
   !call steady state calculations
   CALL SEDSSFNL(CH41,CH42,CH42AV,CH4T1,CH4T2,CH4T2AV,CH41TM1,CH4T2TM1,1) 

   IF(CH42 > CH4SAT) THEN
     CH42=CH4SAT
     CH41 = (CH40*S**2+CH42*KL12*S)/                                  &
            (S**2+KL12*S+XAPPCH4**2*(O20/(KMCH4O2+O20)))        !WLong (4-169) with W2~=0
            
   ENDIF

!**** Calculate changes in CH4 and HS stored in the sediment
   
   IF(.NOT.TF_SSTATE)THEN
        DCH4T2 = (CH4T2 - CH4T2TM1)*H2/DLTS
        DHST2  = (HST2 - HST2TM1)*H2/DLTS  !never used elsewhere??? WLONG (used in James Martin's code)
        
   ELSE
        DCH4T2=0.d0
        DHST2=0.d0
   ENDIF

!**** Calculate CSOD (gO2/M2/day)
   
   CSODCH4 = XAPPCH4**2/S*(O20/(KMCH4O2+O20))*CH41  !WLong (4-170)
   
   CSOD    = CSODCH4+CSODHS


!**** Calculate Fluxes    
   JCH4      = S*(CH41-CH40)
   JCH4AQ    = S*(CH41-CH40)        !same as JCH4
   
!   FLUXHS    = S*FD1*HS1               !WLong (this should be HST1)
   FLUXHS    = S*FD1*HST1  
   FLUXHSCH4 = JCH4AQ + FLUXHS            !Never used

!**** If not flux or SOD or stored then it must escape as gas flux
   JCH4G = 0.
   IF (CH42 == CH4SAT) THEN
     JCH4G = XJ2CH4 - DCH4T2 - CSODCH4 - JCH4AQ
     
     !
     !ie. JCH4G+JCH4AQ + CSODCH4 = XJ2CH4  -DCH4T2 
     !i.e. DMD equation (10.34), pp 204 where XJ2CH4 is JCH4 there
     !
     !==> JCH4G=XJ2CH4-DCH4T2-JCH4AQ-CSODCH4
     !
     
     !Another calculation
     JCH4G = XJC1 - DCH4T2 - DHST2 - CSOD - FLUXHSCH4    !WLong used Jame Martin's method 
     !
     
   ENDIF
   RETURN
   END SUBROUTINE GET_SO4_HS_CH4_CSOD_FLUXES
   
   FUNCTION GET_ZHTANH4F(TEMPVAL)    !m/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTANH4F 
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTANH4F = KAPPNH4F*THTANH4**TEMP202
   RETURN
   END FUNCTION GET_ZHTANH4F
   
   FUNCTION GET_ZHTANH4S(TEMPVAL)   !m/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTANH4S 
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTANH4S = KAPPNH4S*THTANH4**TEMP202
   RETURN
   END FUNCTION GET_ZHTANH4S
   
   FUNCTION GET_ZHTA1HSD(TEMPVAL)     !m/d  
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTA1HSD 
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTA1HSD = KAPP1HSD*THTAH2S**TEMP202
   RETURN
   END FUNCTION GET_ZHTA1HSD
   
   FUNCTION GET_ZHTA1HSP(TEMPVAL)   !m/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTA1HSP 
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTA1HSP = KAPP1HSP*THTAH2S**TEMP202
   RETURN
   END FUNCTION GET_ZHTA1HSP
   
   FUNCTION GET_ZHTANO3F(TEMPVAL)   !m/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTANO3F 
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTANO3F = KAPPNO3F*THTANO3**TEMP202
   RETURN
   END FUNCTION GET_ZHTANO3F
   
   FUNCTION GET_ZHTANO3S(TEMPVAL)   !m/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTANO3S
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTANO3S = KAPPNO3S*THTANO3**TEMP202  
   RETURN
   END FUNCTION GET_ZHTANO3S
   
   FUNCTION GET_ZHTAK2NO3(TEMPVAL)  !DMD p. 105 eqn (4.9) (m/d) (anerobic layer denitrification velocity)
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAK2NO3
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAK2NO3 = K2NO3*THTANO3**TEMP20        !m/d
   RETURN
   END FUNCTION GET_ZHTAK2NO3
   
   FUNCTION GET_ZL12NOM(TEMPVAL)       !temperature contol on diffusion rate of layer 1, 2 (unitless)
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZL12NOM
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZL12NOM = THTADD**TEMP20 !unitless
   RETURN
   END FUNCTION GET_ZL12NOM
    
   FUNCTION GET_ZW12NOM(TEMPVAL)    !temperature control on particle mixing rate of layer 1,2 (unitless)
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZW12NOM
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZW12NOM = THTADP**TEMP20 !unitless    
   RETURN
   END FUNCTION GET_ZW12NOM

   FUNCTION GET_ZHTAPON1(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPON1
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPON1 = KPON1*THTAPON1**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPON1
   
   FUNCTION GET_ZHTAPON2(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPON2
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPON2 = KPON2*THTAPON2**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPON2
   
   FUNCTION GET_ZHTAPON3(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPON3
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPON3 = KPON3*THTAPON3**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPON3
  
   FUNCTION GET_ZHTAPOC1(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPOC1
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPOC1 = KPOC1*THTAPOC1**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPOC1
   
   FUNCTION GET_ZHTAPOC2(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPOC2
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPOC2 = KPOC2*THTAPOC2**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPOC2
   
   FUNCTION GET_ZHTAPOC3(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPOC3
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPOC3 = KPOC3*THTAPOC3**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPOC3
   
   FUNCTION GET_ZHTAPOP1(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPOP1
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPOP1 = KPOP1*THTAPOP1**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPOP1
   
   FUNCTION GET_ZHTAPOP2(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPOP2
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPOP2 = KPOP2*THTAPOP2**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPOP2
   
   FUNCTION GET_ZHTAPOP3(TEMPVAL)   !1/d
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTAPOP3
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTAPOP3 = KPOP3*THTAPOP3**TEMP20
   RETURN
   END FUNCTION GET_ZHTAPOP3

   FUNCTION GET_ZHTASI(TEMPVAL)       !Silicate dissolution rate (1/day), DMD book page P 151 
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTASI
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTASI = KSI*THTASI**TEMP20
   RETURN
   END FUNCTION GET_ZHTASI
   
   FUNCTION GET_ZHTACH4(TEMPVAL)       !m/d !WLong (4-70), DMD (40.39)
   REAL(SP) ::     TEMPVAL             !temperature (degC)
   REAL(SP) ::    GET_ZHTACH4
        TEMP20=TEMPVAL-20.d0
        TEMP202=TEMP20/2.d0
        GET_ZHTACH4 = KAPPCH4*THTACH4**TEMP202    !WLong (4-70)
   RETURN
   END FUNCTION GET_ZHTACH4
   

    SUBROUTINE SEDTEST_UPDATE_WATERCOLUMN(NXJDAY)
        !
        !read the overlying water column data information to make sure overlying water column concentration of current and next time
        !steps are given correctly
        !
        IMPLICIT NONE
        REAL(SP) :: NXJDAY
        !
        !Input file looks like the following (e.g. for test2c copied from SedFlux_ver17b12_test2c.xlsm
        !
!======================        
!Sediment overyling water column condition forcing file
!        Segment ID    Date+time    Jcin        Jnin        Jpin        Jsin        O20            Depth        Tw        NH30        NO30        PO40        SI0            CH40        SALw
!                                  (gO2/m^2/d)    (gN/m^2/d)    (gP/m^2/d)    (gSi/m^2/d)    (mgO2/L)    (m)        (deg C)        (mgN/L)        (mgN/L)        (mgP/L)        (mgSi/L)    (mgO2/L)    (ppt)        
!        1    1/1/00 12:00 AM        0.093450    0.051800    0.000750    0.097125    12.0        4.30    0.0            0.100        0.300        0.040        0.750        1.3    20.0
!        1    1/1/00 6:00 AM        0.092830    0.051066    0.000743    0.095748    12.1        4.32    0.3            0.100        0.300        0.040        0.750        1.3    20.0
!        1    1/1/00 12:00 PM        0.092670    0.050807    0.000740    0.095264    12.1        4.34    0.3            0.100        0.300        0.040        0.751        1.3    20.0
!        1    1/1/00 6:00 PM        0.092371    0.050443    0.000736    0.094580    12.2        4.35    0.4            0.100        0.300        0.040        0.751        1.3    20.0
!        1    1/2/00 12:00 AM        0.092201    0.050215    0.000732    0.094153    12.2        4.35    0.4            0.100        0.301        0.040        0.751        1.3    20.0
!        1    1/2/00 6:00 AM        0.092018    0.049986    0.000729    0.093723    12.2        4.35    0.4            0.100        0.301        0.040        0.751        1.3    20.0
!        1    1/2/00 12:00 PM        0.091872    0.049799    0.000726    0.093373    12.2        4.36    0.4            0.100        0.301        0.040        0.751        1.3    20.0
!        1    1/2/00 6:00 PM        0.091734    0.049624    0.000723    0.093045    12.2        4.36    0.4            0.100        0.301        0.040        0.751        1.3    20.0
!        1    1/3/00 12:00 AM        0.091610    0.049466    0.000721    0.092749    12.3        4.36    0.4            0.100        0.301        0.040        0.751        1.3    20.0
!        1    1/3/00 6:00 AM        0.091494    0.049317    0.000718    0.092469    12.3        4.36    0.4            0.100        0.301        0.040        0.751        1.2    20.0
!        1    1/3/00 12:00 PM        0.091386    0.049176    0.000716    0.092205    12.3        4.36    0.4            0.100        0.301        0.040        0.751        1.2    20.0
!....
!        1    12/30/03 6:00 AM    0.139803    0.059920    0.001033    0.112350    12.2        4.43    0.6            0.073        0.273        0.033        0.647        1.1    20.0
!        1    12/30/03 12:00 PM    0.139610    0.059949    0.001032    0.112405    12.2        4.43    0.6            0.073        0.272        0.033        0.647        1.1    20.0
!        1    12/30/03 6:00 PM    0.139419    0.059978    0.001031    0.112460    12.2        4.43    0.6            0.072        0.272        0.033        0.646        1.1    20.0
!======================        

133   CONTINUE
      DO WHILE(JDAY >= NXWCL)

        !cycle back the record
        JCIN_R1    =    JCIN_R2
        JNIN_R1    =    JNIN_R2
        JPIN_R1    =    JPIN_R2    
        JSIN_R1 =   JSIN_R2
        O20_R1    =    O20_R2        
        D_R1    =    D_R2        
        TW_R1    =    TW_R2 ! + 4.5 ! B Clark increase temperature       
        NH30_R1    =    NH30_R2        
        NO30_R1    =    NO30_R2        
        PO40_R1    =    PO40_R2        
        SIAT0_R1=    SIAT0_R2        
        CH40_R1    =    CH40_R2        
        SALT0_R1=   SALT0_R2 !-10.   

        !read in a new record, if end of file, jump to 134 and open next file
        
        READ (WCL,*,END=134) NXWCL,        &!jday of next record
                                JCIN_R2,    &!C OM flux
                                JNIN_R2,    &!N OM flux
                                JPIN_R2,    &!P OM flux
                                JSIN_R2,    &!Si OM flux
                                O20_R2,        &!O20
                                D_R2,        &!total water depth
                                TW_R2,        &!water temperature
                                NH30_R2,    &!ammonia
                                NO30_R2,    &!nitrate
                                PO40_R2,    &!phosphate concentration 
                                SIAT0_R2,    &!silicate
                                CH40_R2,    &!methane
                                SALT0_R2     !salinity
D_R2 = 1.0
        IF(MSR)THEN
            WRITE(*,1010) NXWCL,        &!jday of next record
                                JCIN_R2,    &!C OM flux
                                JNIN_R2,    &!N OM flux
                                JPIN_R2,    &!P OM flux
                                JSIN_R2,    &!Si OM flux
                                O20_R2,        &!O20
                                D_R2,        &!total water depth
                                TW_R2,        &!water temperature
                                NH30_R2,    &!nitrate
                                NO30_R2,    &!ammonia
                                PO40_R2,    &!phosphate concentration 
                                SIAT0_R2,    &!silicate
                                CH40_R2,    &!methane
                                SALT0_R2     !salinity
                                
            !WRITE(*,*)'here I ''m doing nothing '
            !READ(*,*)
        ENDIF

        NXWCL = (WCLPTR-1)*WCLFDAYS+NXWCL  !global record time in all forcing of WCL
        
      ENDDO
      
	  JCIN_R2 = JCIN_R2!0.04*12*2.667 !JCIN_R2*2   ! B Clark use a constant flux
								! 40 mmol/m^2/day (Brady et al., 2013) * 12 mg C/ mmol C *2.667 mg O2/ mg C
      JNIN_R2 = JNIN_R2!*7.3/6.625!0.04*12/6.625!JNIN_R2*2  !  &!N OM flux
      JPIN_R2 = JPIN_R2!0.04*12/41!JPIN_R2*2  !  &!P OM flux
      JSIN_R2 = JSIN_R2!0.04*12/45.3!JSIN_R2*2    !&!Si OM flux
	  
	  
      GO TO 135

    !
    !Open next data file and read one record if necesary
    !

134   CONTINUE

      WCLPTR = WCLPTR+1
      WRITE (*,*)     'Opening water column condition ',   &
                    'file ',WCLFN(WCLPTR),' at ',        &
                    'day ',JDAY
      CLOSE (WCL)
      
      OPEN (WCL,FILE=WCLFN(WCLPTR),STATUS='OLD')
      
        READ (WCL,1000)         !Read past the header lines
        READ (WCL,*,END=134) NXWCL,        &!jday of next record
                                JCIN_R2,    &!C OM flux
                                JNIN_R2,    &!N OM flux
                                JPIN_R2,    &!P OM flux
                                JSIN_R2,    &!Si OM flux
                                O20_R2,        &!O20
                                D_R2,        &!total water depth
                                TW_R2,        &!water temperature
                                NH30_R2,    &!nitrate
                                NO30_R2,    &!ammonia
                                PO40_R2,    &!phosphate concentration 
                                SIAT0_R2,    &!silicate
                                CH40_R2,    &!methane
                                SALT0_R2     !salinity
        IF(MSR)THEN
            WRITE(*,1010) NXWCL,        &!jday of next record
                                JCIN_R2,    &!C OM flux
                                JNIN_R2,    &!N OM flux
                                JPIN_R2,    &!P OM flux
                                JSIN_R2,    &!Si OM flux
                                O20_R2,        &!O20
                                D_R2,        &!total water depth
                                TW_R2,        &!water temperature
                                NH30_R2,    &!nitrate
                                NO30_R2,    &!ammonia
                                PO40_R2,    &!phosphate concentration 
                                SIAT0_R2,    &!silicate
                                CH40_R2,    &!methane
                                SALT0_R2     !salinity
                                
            !WRITE(*,*)'here I ''m doing nothing 2'    
            !READ(*,*)
        ENDIF

      NXWCL = (WCLPTR-1)*WCLFDAYS+NXWCL !global record time in all forcing of WCL
      
      IF (JDAY >= NXWCL) GOTO 133       !keep reading if JDAY still larger than NXWCL
      
      
135   CONTINUE

      NXJDAY = MAX(NXJDAY,NXWCL)        !global record time in all forcing of WCL
      WRITE(*,*)'NXJDAY=',NXJDAY
1000 FORMAT(//)
!!!---1010 FORMAT(13F8.0,:/(:8X,14F8.0))
1010 FORMAT(14(F12.5,1X))

    END SUBROUTINE SEDTEST_UPDATE_WATERCOLUMN
   
END MODULE MOD_SED



