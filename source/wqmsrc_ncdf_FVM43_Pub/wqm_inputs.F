!Subroutine INPUTS()

    SUBROUTINE INPUTS() 
	
	USE MOD_NCDIO, ONLY: SET_NCD_IO!, NCD_READ_SHAPE
        
       USE MOD_FILEINFO, ONLY : &
                 DIA            &
                ,CBC             &
                ,S1                &
                ,S2                &
                ,S3                &
                ,BFI            &
                ,BAI            &
                ,MET            &
                ,BFO            & 
                ,KEI            &
                ,ATM            &
                ,STL            & 
                ,AGR            & 
                ,SVI            & 
                !,SVO            & 
                !,KFL            & 
                !,ZOO            & 
                !,ZFO            & 
                !,ALO              &
                ,CON            &
                !,RSO            &
                ,SNP            &
                ,PLT            &
                ,APL             &
                ,TFL             &
                ,OPL            &
                ,SFI            &
                ,SFO            &
                !,MAP             &
                ,ICI             &
                !,ICO            &
                ,MRL            &
                ,MBL            &
                ,RSI            &
                ,UNIT_LINKAGE    &
                !,UNIT_STN        &
                !,UNIT_HIS        &            
                ,CNAME            &
                ,CONFN

    USE MOD_LIMS, ONLY:  MTLOC, NTLOC, NPROCS, MYID, NSTATIONMAX, MLOC, KB, KBM1
    USE MOD_SIZES, ONLY: MGL,NGL,NOBTY,NCP
    USE MOD_PREC, ONLY : SP
        
    USE MOD_TGE, ONLY :                 &!
                ISBCE,                    &!
                ISONB,                    &!
                NV,                        &
                TRIANGLE_GRID_EDGE,        &
                TGE_ALLOC
                           
    USE MOD_BCMAP, ONLY: BCMAP
    
    USE MOD_HYDROVARS, ONLY:  &
 
        ART    &        !AREA OF ELEMENT
        ,ART1    &        !AREA OF NODE-BASE CONTROl VOLUME
     
        ,DZ    &        !DELTA-SIGMA VALUE
        !,DZZ    &        !DELTA OF INTRA LEVEL SIGMA 
        !,H1    &        !BATHYMETRIC DEPTH 
        ,H    &            !BATHYMETRIC DEPTH 
        ,D    &            !CURRENT DEPTH 
        ,DT    &        !DEPTH AT PREVIOUS TIME STEP
        ,DT1    &        !DEPTH AT PREVIOUS TIME STEP
        ,EL    &        !CURRENT SURFACE ELEVATION
        ,ET    &        !SURFACE ELEVATION AT PREVIOUS TIME STEP
        ,DTFA    &        !ADJUSTED DEPTH FOR MASS CONSERVATION
        ,UU    &        !X-VELOCITY
        ,VV    &        !Y-VELOCITY
        ,UUT    &        !X-VELOCITY FROM PREVIOUS TIMESTEP
        ,VVT    &        !Y-VELOCITY FROM PREVIOUS TIMESTEP
        !,WWT    &        !Z-VELOCITY FROM PREVIOUS TIMESTEP
        ,WTST    &        !Vertical velocity in sigma from PREVIOUS TIMESTEP
        ,UARD_OBCNT    &!tykim
        ,XFLUX_OBCT    &!tykim
        ,DTFAT    &        !tykim
        !,TT_T    &        !tykim
        !,SALTT    &        !tykim
        ,WTS    &        !VERTICAL VELOCITY IN SIGMA SYSTEM
        ,UARD_OBCN    &    ! tykim 
        ,XFLUX_OBC    &    ! tykim 
        !,WTTS    &        !VERTICAL VELOCITY IN SIGMA SYSTEM 
        ,KH    &        !TURBULENT DIFFUSIVITY
        !,A1U    &        !
        !,A2U    &        !
        !,AWX    &        !
        !,AWY    &        !
        !,AW0    &        !
        ,VISCOFH    &    !
        ,UNC1    &        !
        ,VNC1    &        !
        ,WNC1    &        !
        ,WTSNC1    &        !
        ,UARD_OBCNNC1    &    !
        ,XFLUX_OBCNC1    &    !
        ,DTFANC1    &        !
        ,KHNC1    &        !
        ,TNC1    &        !
        ,SNC1    &        !
        ,ELNC1    &        !
        ,UNC2    &        !
        ,VNC2    &        !
        ,WNC2    &        !
        ,WTSNC2    &    !
        ,UARD_OBCNNC2    &!
        ,XFLUX_OBCNC2    &!
        ,DTFANC2    &    !
        ,KHNC2    &        !
        ,TNC2    &        !
        ,SNC2    &        !
        ,ELNC2    &        !
        !,num_hyd_ints    &!number of records in each hydrodynamics netcdf file
        !,TIME_MAP    &    !
        !,THOUR1    &    !SIMULATION TIME AT END OF CURRENT EXTERNAL STEP (IEXT) IN HOURS
        !,THOUR    &        !
        !,NCFILE_DIR    &!
        ,NCFILE_PREFIX    &!
        ,NCFILE_SUFFIX    &!
        ,NCFILE_NUMBER    &!
        ,FORMAT_STR    &!
        ,hydro_dir       &    ! directory name where hydrodynamics results (netcdf) files are stored
        ,hydro_prefix  &    ! prefix of file name, e.g. 'psm_'
        ,hydro_suffix    &    ! suffix of filename, e.g. '.nc'
        ,hydro_filenumwidth &    ! number of digits in filename following hydro_prefix, e.g. 4 for psm_0002.nc
        ,hydro_filenumstart &    ! starting number of the file name in the digital part of the file name, e.g. 185 for psm_0185.nc
        ,hydro_Nrec    &        ! number of records in each of hydrodynamics file
        !,hydro_dlt    &            ! time step in hydrodynamics file (in seconds), e.g. 100 for 100sec
        !,t_his_start    &        !
        !,t_his_end    &            !
        ,t_his_dlt    &            !starting time, ending time, and interval of history outputs (days)
        ,Nstation    &            !
        ,NstationNum_GL    &    !maximum number of station is NstationMax!
        !,t_stn_start    &        !
        !,t_stn_end    &            !
        ,t_stn_dlt    &            !starting time, ending time, and interval of station outputs (days)
         ,hydro_netcdf    &        !
         ,wqm_history    &        !
         ,wqm_stations    &        !
         ,IFNC    &                !file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise 
         ,NTRECNC    &            !time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file. 
         ,NTHYDRO     &            !overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read
         ,HYDRO_GEOM_ALLOC      &!subroutine that allocates geometry info
         ,HYDRO_ALLOC
         
    !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY :         &
            SERIAL          &           !!TRUE IF SINGLE PROCESSOR
            ,MSR            &           !!TRUE IF MASTER PROCESSOR (MYID==1)
            ,PAR &!            &           !!TRUE IF MULTIPROCESSOR RUN

            ,DTI  !          &              !!internal time step

   
   USE MOD_WQM, ONLY :             &!
                AHMDLT,         &!
                FILGTH,         &!
                JDAY,           &!
                NAPL,           &!
                NDLT,           &!
                NSNP,           &!
                TMEND,          &!
                TMSTRT,         &!
                NDIA,           &!
                NKFL,           &!
                NMBL,           &!
                NOINT,          &!
                NOPL,           &!
                NPLT,           &!
                NRSO,           &!
                NTFL,           &!
                REDS1C,         &!
                REDS1N,         &!
                REDS1P,         &!
                REDS2C,         &!
                TH,             &!
                XYDF,           &!
                ZDFBCK,         &!
                ZDFMUL,         &!        
                ATMFN,          &!
                KEIFN,          &!
                METFN,          &!
                NHYDF,          &!
                NTVDF,          &!
                REDCBC,         &!
                REDCBN,         &!
                REDCBP,         &!
                REDS2N,         &!
                REDS2P,         &!
                REDS3C,         &!
                REDS3N,         &!
                REDS3P,         &!
                S1FN,           &!
                S2FN,           &!
                S3FN,           &!
                SVIFN,          &!
                ZOOFN,          &!
                ALOFN,          &!
                ATMOS_LOADS,    &!
                BAIFN,          &!
                BAOFN,          &!
                BENTHIC_FLUXES, &!
                BFIFN,          &!
                BFOFN,          &!
                BOUNDARY_CONC,  &!
                CBCFN,          &!
                DFIFN,          &!
                DFOFN,          &!
                SAV_CALC,       &!
                SAV_LOADS,      &!
                SETTLING,       &!
                SOURCE_ONE,     &!
                SOURCE_THR,     &!
                SOURCE_TWO,     &!
                SVOFN,          &!
                ZFOFN,          &!
                AC,             &!
                ATMPTR,         &!
                BAOPTR,         &!
                BFIPTR,         &!
                CBCPTR,         &!
                
                HYDPTR,         &!
                IJDAY,          &!
                
                KEIPTR,         &!
                KINETIC_FLUXES, &!
                LIGHT_EXTINCTION,       &!
                METPTR,         &!
                NAC,            &!

                S1PTR,          &!
                S2PTR,          &!
                S3PTR,          &!
                SAVPTR,         &!
                SEDIMENT_CALC,    &!
                BBMI,           &!

                C1,             &!
                C2,             &!
                

                CTEMP,          &!

                DFEEDI,         &!
                EPI,            &!
                
                AANOX,          &!
                ANDC,           &!
                AOCR,           &!
                AONT,           &!
                AREAR,          &!
                BREAR,          &!
                CREAR,          &!
                
                KADPO4,         &!
                KADSA,          &!
                KHCOAG,         &!
                KHNDN,          &!
                KHNNT,          &!
                KHOCOD,         &!
                KHODOC,         &!
                KHONT,          &!
              !  KLDC,           &  ! B Clark moved to wc_dom mod Sep 2015
                KTCOD,          &!
                KTHDR,          &!
                KTMNL,          &!
                KTNT1,          &!
                KTNT2,          &!
                KTSUA,          &!
                TMNT,           &!
                TRCOD,          &!
                TRHDR,          &!
                TRMNL,          &!
                TRSUA,          &!

              !  KLDN,           &! ! B Clark moved to wc_dom mod Sep 2015
                KLPC,           &!
             !   KRDC,           &!  ! B Clark moved to wc_dom mod Sep 2015
              !  KRDN,           &!! B Clark moved to wc_dom mod Sep 2015
                KRPC,           &!

              !  KLDP,           &! ! B Clark moved to wc_dom mod Sep 2015
                KLPN,           &!
                KLPP,           &!
               ! KRDP,           &! ! B Clark moved to wc_dom mod Sep 2015
                KRPN,            &
                
                KCOD,           &!
                KDCALG,         &!
                KLCALG,         &!
                KRPP,           &!
                KSUA,           &!
                
                KDNALG,         &!
                KDPALG,         &!
                KLNALG,         &!
                KLPALG,         &!
                KRCOAG2in,         &!
                KRCOAG3in,         &!
                NCB,            &!
                NTM,            &!
                S1LB,           &!
                S1LN,           &!
                S2LB,           &!
                S2LN,           &!
                S3LB,           &!
                S3LN,           &!
                WS1,            &!
                WS2,            &!
                WSL,            &!
                WSR,            &!
                WSS,            &!
                
                Ksed,  &!
                S_vmax, &!
                SedOriginal, &!
                SedHyp, &!

                WSSHI,   &!
                WSSLIM,  &!
                FLOW,           &!
                KHSO,           &!
                KSDOC,          &!
                KSNH4,          &!
                KSNO3,          &!
                KSO,            &!
                KSPO4,          &!
                KSSA,           &!
                MTCNO3,         &!
                SEDNO3,         &!
                TRSDOC,         &!
                TRSNH4,         &!
                TRSNO3,         &!
                TRSO,           &!
                TRSPO4,         &!
                TRSSA,          &!
                WS3,            &!
                WSU,            &!

                AVERAGE_PLOTS,          &!
                CONSERVE_MASS,          &!
                DIAGNOSTICS,            &!

                JDAYMBL,                &!
                KFLDP,                  &!
                MASS_BALANCE,           &!
                MBLDP,                  &!
                NIT,                    &!
                PIP_CALC,               &!
                QUALITY_DIAG,           &!
                SAV_PLOTS,              &!
                SEDIMENT_DIAG,          &!
                SOLIDS_CALC,            &!
                STEP_BOUNDARY,          &!
                TRANSPORT_FLUXES,       &!
                XY_DIFFUSION,           &!
                Z_DIFFUSION,            &!
                DLT,            &!
                ELTMS,          &!
                ELTMSKFL,       &!
                ELTMSPLT,       &!
                ELTMSTFL,       &!
                FNDLT,          &!
                NHMR,           &!
                NWQMR,          &!
                NXATM,          &!
                NXBFI,          &!
                NXCBC,          &!
                NXMET,          &!
                NXS1,           &!
                NXS2,           &!
                NXS3,           &!
                NXSAV,          &!

                ANC1,           &!
                ANC2,           &!
                ANC3,           &!

                ASC1,           &!
                ASC2,           &!
                ASC3,           &!
                CMASS,          &!
                DTC,            &!
                INFLOW,         &!
                IWCMC,          &!
                IWCMN,          &!
                IWCMP,          &!
                IWCMS,          &!
                JCS1MAX,        &!
                JCS2MAX,        &!
                JCS3MAX,        &!
                
                NHQF,           &!
                WQM_ALLOC
                
   USE MOD_OWQ, ONLY:   INITKE,            &    !
                        INTKE,            &    !
                        NXKEI,            &    !
                        LAT,             &    !
                        LON,            &    !
                        GSTAR440,        &    !
                        PHISTAR676,        &    !
                        BNSTARNTU,        &    !
                        BNSTARCHL,      &    !
                        PSTARINT,        &    !
                        PSTAR440,        &    !
                        PSTARCOLOR,        &    !
                        SCDOM,            &    !
                        EZEROINT,        &    !
                        TCHL2CHL,        &    !
                        DOFFSET,        &    !
                        DEPTHKE,        &    !
                        NWAVEL,            &    !
!                        OPZ,             &    !
!                        G1,                &    !
!                        G2,                &    !
                        ALAMB,             &    !
                        GLAMB,            &    !
                        PLAMB,            &    !
                        PHILAMB,        &    !
                        BLAMB,            &    !
                        WAVEL,            &    !
                        EZERO,            &    !
!                        GMUZERO,        &    !
!                        IAVG,            &    !
!                        IATBOT,            &    !
                        OWQ_ALLOC,        &    !
                        !OWQ_DEALLOC,    &    !
                        OWQ_READ!,        &    !
                        !LGHT_ATTN
   
   USE MOD_SF, ONLY : SFEEDER,    &
                    NSPECIES,    &!

                    HYPOXFX_SF,    &!

                    SFEED,        &!
                    SFEED_GL,    &!(:,:)

                    SF_ALLOC,   &!
                    SF_READ,    &!
                    SF_INIT
                    
    USE MOD_DF, ONLY:            &
            DFEEDER,            &!
            DFEEDER_OUTPUT,        &
            HYPOXFX_DF,            &!
            DFEEDM1S_GL,        &!
            DFEEDM1S,            &!
            DF_READ,            &!
            DF_INIT,            &!
            !DF_DEALLOC,        &!
            DF_ALLOC    
			
#if defined(SEDIMENT_DOM)
USE MOD_SED_DOM, ONLY: 	SED_DOM_ALLOCATE_VARS,   &
		                   SED_DOM_INITIALIZE,   &
						   SED_DOM_INPUT
	USE MOD_SED_DOM_EXCHANGE, ONLY : SED_DOM_FLAG		
	
# endif	
   USE WC_DOM, ONLY : WC_DOM_ALLOCATE, WC_DOM_INPUT
						   
   USE MOD_WQMINIT, ONLY :        &        !
            RSODP,                 &        !
            DLTDP,               &        !
            SNPDP,                &        !
            TFLDP,                &        !
            PLTDP,                &        !
            APLDP,               &        !
            OPLDP,                &        !
            DIADP,                &        !
            !COURB,                &        !
            !COURFS,            &        !
            !COURBS,            &        !
            F,                    &        !            
            !SB,                &        !
            S1LNMAX,            &        !
            S2LNMAX,            &        !
            S3LNMAX,            &        !
            !DIFFFS,            &        !
            !DIFFBS,              &        !
            NXSNP,                 &        !
            NXPLT,                &        !
            NXAPL,                &        !
            NXTFL,                &        !
            NXKFL,                &        !
            NXTVD,               &        !
            NXOPL,                &        !
            NXMBL,                &        !
            NXDIA,                &        !

            LEAFI,                &        !
            STEMI,                &        !
            ROOTI,                 &        !
            TUBERI,                &        !
            !MAXDLT,            &        !
            MXDLT,                &        !
!            INTKE,                &        !    !WLong moved this to mod_owq.F
!            INITKE,                &        !    !WLong moved this to mod_owq.F
            !TM1,                &        !
            !TM2,                &        !
            !TM3,                &        !
            MINSTEP,            &        !
!            NXKEI,                &        !
            SNPC,                &        !
            RSIC,                 &        !
            RSOC,                &        !
            RSOD,                &        !            
            BCC,                &        !
            S1C,                   &        !
            S2C,                &        !
            S3C,                  &        !
            MDC,                &!*3,    !    
            PLTC,                &        !
            FLC,                   &        !
            MBLC,                &        !
            BFC,                &        !
            VBC,                &        !
            QPLTC,                &        !
            XYDFC,                &        !
            XYDFU,              &        !
            ZDFC,                &        !
            ICOC,                &        !
            ATMC,                &        !
            SAVLC,                &        !
            SEDC,                  &        !
            AUTOC,                &        !
            SPLTC,                &        !
            TFLC,                &        !
            DIAC,                &        !
            STLC,                  &        !
            APLTC,                &        !
            KFLC,                &        !
            OPLC,                &        !
            BFOC,                &        !
            BAOC,               &       ! WLong benthic algae output flag
            DFOC,                &        ! WLong, deposition feeder output flag
            SAVMC,                &        !
            SAVPLTC,            &        !
            DFLC,                &        !
            SFLC,                &        !
            SFLOXC,                &        ! WLong: hypoxia effects ON/OFF on suspension feeder
            DFLOXC,                &       ! WLong: hypoxia effects ON/OFF on deposition feeder
                                 
            !EXT1,                 &!*1,   !    
            !EXT2,                &!*2,    !    
            !EXT3,                &!*3,    !    
            !EXT4,                &!*4,    !    
            !EXT_R,                &!*1,    !
            !EXT_R2,            &!*2    !
            
            SLC,                &        !
            !HYDC,                &        !
            BNDTC,                &        !
            CONSC,                &        !
            ICIC,                &        !
                                
            !TIMEA,                &        !
            !TIMEB,                &        !
            SPVARM,                &        !
            PRINTM,                &        !
                                
            TITLE,                &!(6),    !
            OLDTITLE,            &!(6),    !
            !FILENAME,            &        !
                                    
            !FILENAME2,            &        !
                
            MAPFN,                &        !
            GEOFN,                &        !
            ICIFN,                &        !
            RSIFN,              &        !
            AGRFN,                &        !
            STLFN,                &        !
            MRLFN,                &        !
            KFLFN,                &        !
            ICOFN,                &        !
            
            SNPFN,                &        !
            PLTFN,                &        !
            APLFN,                 &        !
            DIAFN,                 &        !
            TFLFN,              &        !
            RSOFN,                 &        !
            OPLFN,                 &        !
            MBLFN,                 &        !
            SFIFN,                 &        !
            SFOFN,                &        !

            RESTART_OUT,        &        !
            SNAPSHOTS,            &        !
            END_RUN,            &        !
            MODIFY_ICONC,       &          !
            VOLUME_BALANCE,     &        !
            QUICKEST,            &        !
            UPWIND,                &        !
            ICOND_OUT,          &        !
            UNI_ICON_IN,           &        !
            BIN_ICON_IN,        &        !
            RES_ICON_IN,        &        !
            AUTO_STEPPING,      &        !
            !STOP_RUN,           &        !
            PLOTS,                &        !
            OXYGEN_PLOTS,       &        !
            ZOO_CALC,           &        !
            !NEW_VOLUMES,          &        !
            RESTART_IN,            &        !
                                
            TEMPERATURE_CALC,     &        !
            ALGAE_CALC,            &        !
            CARBON_CALC,           &        !
            NITROGEN_CALC,         &        !
            PHOSPHORUS_CALC,       &        !
            COD_CALC,            &        !
            OXYGEN_CALC,        &        !
            SILICA_CALC,        &        !
                                
            LEFT_FLOWB,         &!(NHQP),        !        
            RIGHT_FLOWB,        &!(NHQP),        !    
            !LEFTM1_BOUNDARY,    &!(NHQP),         !
            !RIGHTP1_BOUNDARY,    &!(NHQP)          !

            !DEN1,                &        !
            !DEN2,                &        !
            !DEN3,               &        !
            !TP1,                 &        !
            !TP2,                &        !
            !TP3,                &        !
            !T2,                &        !
            !SF2,                &        !
                                        
            !T1,                &        !
            !T3,                &        !
            !SF1,                &        !
                    
            !COUR,                &        !
            
            !GRAD,                &    !(NHQP,3),        !
            !TERM,                &    !(NHQP,3),        !

            !GRAD1,                &        !
            !GRAD2,                &        !
            !GRAD3,                &        !
            !TERM1,                &        !
            !TERM2,                &        !
            !TERM3,                &        !

            IFLOWP,                &    !(NBCP),        !
            !IBT,                &    !(NHQP),          !
            !COURBSV,            &    !(NHQP),        !
            !COURVSV,            &    !(NHQP)            !

            ! C1MIN,              &                !
            ! C1MAX,                &                !

            DOVDAYS,            &    !(:,:,:)    !
            OINT,                &    !(NOIP)        !
            ACC,                &    !(NCP)        !
             
            CIC,                &    !(NCP),        !
            !COUT,                &    !(NCP),        !
            SFEEDI,                &    !(10)        !
            
            DLTD,                &!                
            SNPD,                &!
            
            SNPF,                 &!
            DLTVAL,                &!
            DLTMAX,               &!
            DLTFTN,             &!
            PLTD,                &!
            PLTF,                 &!
            APLTD,                &!
            APLF,                &!
            TFLD,                 &!
            TFLF,                &!
            KFLD,                &! 
            KFLF,                 &!
            OPLD,                &!
            OPLF,                &!
            DIAD,                 &!
            DIAF,                 &!
            MBLD,                 &!
            MBLF,                &!

            !DOCLIT ,             &!
            !LPOCLIT,              &!
            !POCLIT,              &!
            !PBSLIT,              &!
            !DSILLIT,            &!
            
            WQMINIT_ALLOC        !Subroutine to initialize MOD_WQM data
            
            
   USE MOD_ALGAL, ONLY :     & 
        ALG_READ!,            &
        !ALGAE    

   USE MOD_BA, ONLY : BBM,                 &
                      BBM_GL,            &
                      BALGAE_CALC,        &
                      BA_OUTPUT,        &
                      BALC
        
   USE MOD_SED, ONLY :         &!
           BENTHIC_OUTPUT,     &!
    !    STEADY_STATE_SED_G3,   &!
           STEADY_STATE_SED_IC, &  !added by LB!
    !    DLTS,                &!
    !    SSTATEG3,               &!
    !    HSEDALL,               &!
    !    DIFFT,                &!
    !    SALTSW,             &!
    !    SALTND,                &!    
    !    FRPALG1,                &! 
    !      FRPALG2,                &! 
    !      FRPALG3,                &! 
    !    FRNALG1,                &! 
    !   FRNALG2,               &! 
    !   FRNALG3,               &!
    !    FRCALG1,                &!
    !   FRCALG2,               &!
    !   FRCALG3,               &!
    !   FRCPHB,               &!
    !   FRNPHB,               &!
    !   FRPPHB,               &!
    !    KPDIAG,             &!
    !   DPTHTA,             &!  
    !   KNDIAG,             &!  
    !   DNTHTA,             &!  
    !   KCDIAG,             &!  
    !   DCTHTA,             &!  
    !    KSI,                &!  
    !   THTASI,             &!  
    !    M1,                   &! 
    !   M2,                   &! 
    !   THTADP,               &! 
    !   THTADD,               &!
    !    KAPPNH4F,             &! 
    !   KAPPNH4S,             &! 
    !   PIENH4,               &! 
    !   THTANH4,              &! 
    !   KMNH4,                &! 
    !   KMNH4O2,              &! 
    !    KAPPNO3F,              &!  
    !   KAPPNO3S,              &!  
    !   K2NO3,                 &!  
    !   THTANO3,               &! 
    !    KAPP1HSD,             &! 
    !   KAPP1HSP,            &! 
    !   PIE1HS,             &! 
    !   PIE2HS,             &! 
    !   THTAH2S,               &! 
    !   KMHSO2,             &!  
    !    CSISAT,              &! 
    !   DPIE1SI,             &! 
    !   PIE2SI,              &! 
    !   KMPSI,               &!  
    !    O2CRITSI,             &! 
    !   JSIDETR               &!
    !    DPIE1PO4F,          &!
    !   DPIE1PO4S,            &!
    !   PIE2PO4,               &!
    !   O2CRITPO4,            &!
    !   KMO2DP,                &!
    !    TEMPBEN,             &!
    !   KBENSTR,             &!
    !   KLBNTH,              &!
    !   DPMIN,               &!
    !    KAPPCH4,             &!
    !   THTACH4,             &!
    !   KMCH4O2,            &!
    !   KMSO4,               &! 
    !    XKMI0,                &!
    !    ING0,                 &!
    !    THTAI0,               &!
    !    R,                    &!
    !    THTAR,                &!
    !    BETA,                 &!
    !    THBETA,               &!
    !    AMCN,                 &!
    !    AMCP,                 &!
    !    AA1,                  &!
    !    AA2,                  &!
    !    XKMG1,                &!
    !    XKMG2,                &!
    !    XKBO2,                &!
    !    TDD,                  &!
    !    DFDOH,                &!
    !   DFDOQ,                  &!
    !    DOLOW,                &!
    !    BALC,                &!
    !    PMB,                 &!
    !    ANCB,                 &!
    !    APCB,                 &!
    !    KTGB1,                 &!
    !    KTGB2,                 &!
    !    TMB,                   &!
    !    KESED,                 &!
    !    KEBALG,             &!
    !    KHNB,                 &!
    !    KHPB,                 &!
    !    KHRB,                 &!
    !    BMRB,                 &!
    !    BPRB,                 &!
    !    KTBB,                &!
    !    TRB,                  &!
    !    BALGMIN,             &!
    !    FNIB,                 &!
    !    FPIB,                  &!
    !    ALPHB,                &! 
    !   CCHLB,                 &!  
    !    SPVARS,               &!
    !    PRINTS,                &!
    !    WSSBNET,              &!
    !    WSLBNET,              &!
    !    WSRBNET,              &!
    !    WS1BNET,              &!
    !    WS2BNET,              &!
    !    WS3BNET,              &!
    !    WSUBNET,            &!
    !    SPVARB,             &!
    !    PRINTB,                &!
    !    VSED,                  &!
    !    VPMIX,                &!
    !    VDMIX,                &!
    !    SPVARLR,              &!
    !    PRINTLR,            &!
    !    FRPOP,                 &!
    !    FRPON,                 &! 
    !   FRPOC,                &! 
    !   DPP,                &!   
    !   W2,                 &!   
    !   DDP,                &!   
    !   H2,                   &!   
    !   DPIE1PO4,           &!   
       !KAPPNH4,               &!   
       !KAPPNO3,            &!      
        !KWC,                 &! 
        !IERR,                &!
        !ITEMP,                &!
        !SSNAME,            &!
        !K0H1D,               &!
        !K0H1P,               &!
        !K1H1D,               &!
        !K1H1P,               &!
        !K2H2D,               &!
        !K2H2P,               &!
        !K3,                  &!
        !J1,                  &! 
        !J2,                  &!  
        !KMC1,                &!
        !KL12,                &!
        !KL12NOM,             &!
        !KHS_1,               &!
        !KL12SO4,             &!
        !NH40,                 &!
        !NH41,                 &!
        !NH42,                 &!
        !NH4T1,                &!
        !NH4T2,                &!
        !NH41TM1,              &!
        !NH4T2TM1,             &!
        !JNH4,                &!    
        !NO30,                 &!
        !NO31,                 &!
        !NO32,                 &!
        !NO3T1,                &!
        !NO3T2,                &!
        !NO31TM1,              &!
        !NO3T2TM1,             &!
        !JNO3,                &!
        !JSI,                  &!
        !JPO4,                &!
        !JCH4AQ,               &!
        !JO2NH4,               &!
        !JCH4G,                &!
        !CH40,                 &!
        !CH41,                 &!
        !CH42,                 &!
        !CH42AV,               &!
        !CH4T1,                &!
        !CH4T2,                &!    
        !CH4T2AV,              &!
        !CH41TM1,              &!
        !CH4T2TM1,             &!
        !JCH4,                 &!
        !CH4T2SAV,             &!    
        !JCH4GASO2,            &!     
        !JGAS,                &!      
        !SO40,                 &!
        !SO41,                 &!
         !SO42,                 &!
        !SO42AV,               &!
        !SO4T1,                &!
         !SO4T2,                &!
        !SO4T2AV,              &!
        !SO40MG,               &!
         !SO41TM1,              &!
        !SO4T2TM1,             &!
         !JSO4,                 &!
        !SO4T2SAV,          &!    
        !KPOP1,                &!
        !KPOP2,                &!
        !KPOP3,                &!     
        !KPON1,                &!      
        !KPON2,                &!      
        !KPON3,                &!      
        !KPOC1,                &!      
        !KPOC2,                &!
        !KPOC3,             &!     
        !HS0,                  &!
         !JHS,                &!
        !ISEDMN,               &!
        !ISEDMP,               &!
        !ISEDMC,            &!
        !NH4AVL,               &!    
        !NO3AVL,               &!
        !KETOTL,            &!
        !WSSINETMMD,             &!
        !WLPOMNETMMD,            &!
        !WRPOMNETMMD,            &!
        !WB1NETMMD,            &!
        !WB2NETMMD,            &!   
        !WB3NETMMD,            &!   
        !WPSINETMMD,            &!   
        !STP20,             &!
        !DF,                 &!   
        !PO4AVL,             &!   
        !BFOR,                 &! 
        !W12MIN,            &!             
        !LOGICT,             &! 
        !IK,                &!
        !XKPOP1,             &!
        !XKPOP2,             &!
        !XKPOP3,             &!
        !XKPON1,             &!
        !XKPON2,             &!
        !XKPON3,              &!
        !XKPOC1,             &!
         !XKPOC2,             &!
         !XKPOC3,            &!
        !RDD,                 &!
         !RMORT,                &!
         !XPOC1LIM,            &!
         !TEMP,                &!
         !XPOC2LIM,            &!
         !DLTS_H2,                &!
        !FRPON1,             &!
        !FRPOP1,             &!
        !FRPOC1,             &!

        !ERROR,                &!
         !A1,                &!
         !A2,                &!
        !XJCNO31,            &!
         !FP1SO4,            &!
         !FP2SO4,            &!
         !HS2AV,                &!
         !XJ2,                &!
         !XJ2CH4,            &!
         !X1J2,                &!
        !PF,                &!
         !PPO4,                &!

         !PRNB,                &!
         !FRDOB,                &!
        !FLUXHS,            &!
         !FLUXHSCH4,            &!
         !VJCH4G,            &!
        !PON1TM1S,           &!    
        !PON2TM1S,           &!     
        !PON3TM1S,           &!      
        !POC1TM1S,           &!      
        !POC2TM1S,           &!      
        !POC3TM1S,           &!      
        !POP1TM1S,           &!      
        !POP2TM1S,           &!      
        !POP3TM1S,           &!      
        !PSISEDTM1S,            &!      
        BENSTRTM1S,           &!
        SODTM1S,            &!
    !    BFORMAXS,           &!      
    !    ISWBENS             &!    
        NH41TM1S,              &!  
        NH4T2TM1S,         &!  
        NO31TM1S,          &! 
        NO3T2TM1S,         &!  
        HS1TM1S,           &! 
        HST2TM1S,              &! 
        SI1TM1S,           &! 
        SIT2TM1S,              &!  
        PO41TM1S,          &!  
        PO4T2TM1S,             &!  
        CH4T2TM1S,         &!  
        CH41TM1S,          &!  
        SO4T2TM1S,            &!     
    !    BURIALC,            &!     
    !    BURIALN,            &!     
    !    BURIALP,               &!
    !    DIAGENC,              &!   
    !    ZHTADP,             &!  
    !    ZHTADD,             &!  
    !   ZHTANH4F,           &! 
    !    ZHTANH4S,           &!  
    !    ZHTANO3F,           &!  
    !    ZHTANO3S,           &!  
    !  ZHTAK2NO3,               &! 
    !   ZHTA1HSD,             &! 
    !    ZHTA1HSP,             &!  
    !    ZHTASI,             &!  
    !    ZL12NOM,            &!  
    !    ZW12NOM,            &!  
    !    ZHTAPON1,           &!  
    !    ZHTAPON2,           &!  
    !    ZHTAPON3,           &! 
    !    ZHTAPOC1,           &!
    !    ZHTAPOC2,           &!
    !   ZHTAPOC3,           &!
    !    ZHTAPOP1,           &!
    !    ZHTAPOP2,           &!
    !    ZHTAPOP3,           &!
    !    ZHTACH4,            &!
    !    THTAPOP1,              &! 
    !    THTAPOP2,              &! 
    !    THTAPOP3,              &! 
    !    THTAPON1,              &! 
    !    THTAPON2,              &! 
    !   THTAPON3,              &! 
    !   THTAPOC1,              &! 
    !   THTAPOC2,              &! 
    !   THTAPOC3,            &!
    !    AG3CFL,              &! 
    !   AG3NFL,              &! 
    !   AG3PFL               &!  
    !    ASDTMP,             &!
    !    FIB,                  &!  
    !   NLB,                  &! 
    !    PLB,                  &!  
    !   NPPB,                 &! 
    !    BBM,                  &!  
    !   BLITE,                &! 
    !   BMB,                  &! 
    !   PB,                   &! 
    !   PRB,                  &! 
    !   BANH4,                &! 
    !   BANO3,                &! 
    !   BAPO4,                &! 
    !   BADOC,                &! 
    !   BADO,                 &! 
    !   BAPOC,                &! 
    !   BAPON,                &! 
    !   BAPOP,                &! 
    !    HS1,                   &! 
    !   HS2,                   &!
    !   HST1,                  &!
    !   HST2,                  &!
    !   HS1TM1,                &!
    !   HST2TM1                &!
    !    SI0,                   &!
    !    SI1,                   &!
    !    SI2,                   &!
    !    SIT1,                  &! 
    !    SIT2,                  &!     
    !    SI1TM1,                &! 
    !    SIT2TM1,               &! 
    !    PO40,                  &! 
    !   PO41,                  &! 
    !   PO42,                  &! 
    !   PO4T1,                 &! 
    !   PO4T2,                 &! 
    !   PO41TM1,               &! 
    !   PO4T2TM1,              &! 
    !    PON1,                  &! 
    !    PON1TM1,               &!
    !    PON2,                  &! 
    !    PON2TM1,               &!
    !    PON3,                  &! 
    !    PON3TM1,               &!
    !    POC1,                  &! 
    !    POC1TM1,               &! 
    !    POC2,                  &! 
    !    POC2TM1,               &!
    !    POC3,                  &! 
    !    POC3TM1,               &!
    !    POP1,                  &! 
    !    POP1TM1,               &!
    !    POP2,                  &! 
    !    POP2TM1,               &!
    !    POP3,                  &!
    !    POP3TM1,               &!
    !    JNX,                   &!
    !    JCX,                   &!
    !    JPX                    &!
    !    PSISED,                &!
    !    PSISEDTM1                 &!
    !    XJCNO3,                &!
    !    XJCO2,                 &! 
    !    XJC1,                  &!                          
    !    PIE1,                  &! 
    !    PIE2,                  &! 
    !    W12,                   &! 
    !    TEMPD,                 &! 
    !    O20,                   &! 
    !    CH4SAT,                &! 
    !    SAL,                 &!
    !    XAPPNH4,            &! 
    !    XAPP1HSD,             &!
    !    XAPP1HSP,             &! 
    !    XAPP1NO3,           &! 
    !    XK2NO3,             &! 
    !    XKSI,               &! 
    !    XAPPCH4,            &! 
    !    TEMP20,             &! 
    !    TEMP202,            &! 
    !    FD1,                &! 
    !    FP1,                &! 
    !    FD2,                &! 
    !    FP2,                &! 
    !    SOD,                &! 
    !    CSOD,               &! 
    !    S,                  &! 
    !    W12NOM,             &!                      
    !    HSO4,               &! 
    !    DDSO4,              &! 
    !    CSODHS,             &! 
    !    CSODCH4,            &! 
    !    BENSTR,             &! 
    !    BENSTRS,            &! 
    !    BENSTRTM1,            &! 
    !    ISWBEN,             &! 
    !    BFORMAX,            &! 
    !    ZHTANH4,            &! 
    !   ZHTANO3,            &! 
    !    TINTIM,                &! 
        CTEMP_GL,            &!
        CPOP_GL,            &!
        CPON_GL,            &!
        CPOC_GL,            &!
        CPOS_GL,            &!
        PO4T2TM1S_GL,        &!
        NH4T2TM1S_GL,        &!
        NO3T2TM1S_GL,        &!
        HST2TM1S_GL,        &!
        CH4T2TM1S_GL,        &!
        CH41TM1S_GL,        &!
        SO4T2TM1S_GL,        &!
        SIT2TM1S_GL,        &!
        BENSTRTM1S_GL,        &!
        SODTM1S_GL,         &!
                CPOC,       &!
                CPON,       &!
                CPOP,       &!
                CPOS,       &!
                CPOPI,        &!Initial conditions
                CPONI,         &!
                CPOCI,        &!
                CPOSI,        &!    
                PO41TI,        &!WLong 
                PO4T2I,        &!
                NO31TI,        &!WLong
                NO3T2I,        &!
                NH41TI,        &!WLong
                NH4T2I,        &!
                CH4T2I,        &!
                CH41TI,        &!
                SO41TI,        &!WLong
                SO4T2I,        &!
                HS1TI,        &!WLong
                HST2I,         &!
                SI1TI,        &!WLong
                SIT2I,         &!
                CTEMPI,        &!
                BENSTI,        &!
                SODI,       &!WLong
        SED_INIT,            &!
        SED_INIT2,            &!
    !    SED_ALLOC,            &!
    !    SED_DEALLOC,        &!
        SED_READ,            &!
         SED_CALC!,            &!!LB
    !    SED_DIAGENESIS_G3,            &!
    !    SEDTSFNL,            &!
    !    SEDSSFNL,            &!
       !    SEDF,                &!
    !    ZBRENT        
   
   USE MOD_SAV, ONLY : NSAVM,   &             !variables used 
                      ALEAF,    &
                      AROOT,    &   
                      ASTEM,    &                                      
                      ATUBER,   &                                      
                      AEP,      &
                      APLEAF,   &
                      ABMLEAF,  & 
                      ABMTUBER, &
                      APEP,     &
                      ABMEP,    &
                      APREP,    &
                      ASLSH,    &
                      ANLSAV,   &
                      APLSAV,   &
                      ANLEPI,   &
                      APLEPI,   &
                      AFNSED,   &
                      AFPSED,   &
                      AFHS,     &                          
                      AEPATN,   &                                     
                      AWATATN,  &
                      AFISH,    &
                      AFIEP,    &                                                                                                                
                      ANPPSAV,  &
                      ANPPEPI,  &
                      ADOSAVW,  &  !average DO  flux due to SAV gO2/m^2/day
                      ADOCSAVW, &  !average DOC flux due to SAV
                      APOCSAVW, &  !average POC flux due to SAV
                      ANH4SAVW, &  !average NH4 flux due to SAV
                      ANO3SAVW, &  !average NO3 flux due to SAV
                      ADONSAVW, &  !average DON flux due to SAV              
                      APONSAVW, &  !average PON flux due to SAV
                      APO4SAVW, &  !average PO4 flux due to SAV
                      ADOPSAVW, &  !average DOP flux due to SAV
                      APOPSAVW, &  !average POL flux due to SAV
                      ADOEPIW,  &  !
                      ADOCEPIW, &  !
                      APOCEPIW, &
                      ANH4EPIW, &
                      ANO3EPIW, &
                      ADONEPIW, &
                      APONEPIW, &                                      
                      APO4EPIW, &
                      ADOPEPIW, &
                      APOPEPIW, & 
                     ASEDDOSAV, &  !DO  flux due to SAV,  positive leaving sediments 
                    ASEDPOCSAV, &  !POC flux due to SAV, positive into sediments
                    ASEDPONSAV, &  !PON flux due to SAV, positive into sediments
                    ASEDPOPSAV, &  !POP flux due to SAV, positive into sediments
                    ASEDNH4SAV, &  !NH4 flux due to SAV, positive leaving sediemnts
                    ASEDPO4SAV, &  !PO4 flux due to SAV, positive leaving sediments   
                !!subroutines used from MOD_SAV
                SAV_ALLOC,        &    !subroutine to allocate SAV variables
                !SAV_DEALLOC,    &    !subroutine to deallocate SAV variables
            SAV_LOADS_ALLOC,    &    !
            !SAV_LOADS_DEALLOC,    &    !subroutine to deallocate SAV loads
                SAV_READ!,        &    !
                !SAV_COMP            !subroutine that does SAV calculation
          
#  if defined (MULTIPROCESSOR)
            USE MOD_PAR, ONLY:   &      
            EL_PID,    &        !(:)           !!PROCESSOR OWNER OF GLOBAL ELEMENT
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
            HN_LST,     &    !            !!LIST OF HALO NODES 
            !HN_OWN,     &    !            !!PRIMARY OWNER OF HALO NODES
            EC,            &    !
            NC,            &    !
            !BNC,         &    !
            !EMAP,        &    !
            !NMAP,         &    !
            !PNE,        &    !(:)        !!NUMBER OF ELEMENTS IN EACH PROC
            !PNN,        &    !(:)        !!NUMBER OF NODES IN EACH PROC
            !PNHE,        &    !(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
            !PNBN,        &    !(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
            !PMBM,        &    !(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
            !PNHN,        &    !(:)       !!NUMBER OF HALO NODES IN EACH PROC
            !FILELOCK,   & !
            !MYFILELOCK, & !
            !FILELOCKED , & !
            ! NODE_MATCH, & !
            ! NODE_MATCH_ALL, & !
             EXCHANGE!, & !
            !EXCHANGE_ALL, & !
            !GATHER, & !
            !IGATHER, & !
            !INIT_MPI_ENV,&!             
            !SORT,    & !
            !PREAD,    & !
            !PWRITE,& !
            !PPRINT, & !
            !LOCK, & !
            !UNLOCK, & !
            !GETLOC !

#  endif

   USE MOD_ZOOP, ONLY : &            
                     ANCSZ,      &!microzooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     ANCLZ,      &!mesozooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     APCSZ,      &!microzooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
                     APCLZ,      &!mesozooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
                     !AOCRSZ,     &!ratio of oxygen consumed to microzooplankton carbon metabolized (gO2/gC) (~2.67)
                     !AOCRLZ,     &!ratio of oxygen consumed to mesozooplankton carbon metabolized (gO2/gC) (~2.67)
                     !FRSASZ,     &!fraction of microzooplankton silica recycled to dissolved silica pool, range [0,1]
                     !FRSALZ,   &!fraction of mesozooplankton silica recycled to dissolved silica pool, range [0,1]
                     !FLDOCSZ,    &!fraction of microzooplankton carbon released to LDOC, range [0,1]
                     !FRDOCSZ,    &!fraction of microzooplankton carbon released to RDOC, range [0,1]
                     !FLPOCSZ,    &!fraction of microzooplankton carbon released to LPOC, range [0,1]
                     !FRPOCSZ,    &!fraction of microzooplankton carbon released to RPOC, range [0,1]
                     !FLDONSZ,    &!fraction of microzooplankton nitrogen released to LDON, range [0,1]
                     !FRDONSZ,    &!fraction of microzooplankton nitrogen released to RDON, range [0,1]
                     !FLPONSZ,    &!fraction of microzooplankton nitrogen released to LPON, range [0,1]
                     !FRPONSZ,    &!fraction of microzooplankton nitrogen released to RPON, range [0,1]
                     !FLDOPSZ,    &!fraction of microzooplankton phosphorus released to LDOP range [0,1]    
                     !FRDOPSZ,    &!fraction of microzooplankton phosphorus released to RDOP range [0,1]    
                     !FLPOPSZ,    &!fraction of microzooplankton phosphorus released to LPOP range [0,1]    
                     !FRPOPSZ,    &!fraction of microzooplankton phosphorus released to RPOP range [0,1]    
                     !FNH4SZ,     &!fraction of microzooplankton nitrogen recycled to DIN as NH4, range [0,1]
                     !FPO4SZ,     &!fraction of microzooplankton phosphorus recycled to DIP as PO4, range [0,1]
                     !FLDOCLZ,    &!fraction of mesozooplankton carbon released to LDOC, range [0,1]
                     !FRDOCLZ,    &!fraction of mesozooplankton carbon released to RDOC, range [0,1]
                     !FLPOCLZ,    &!fraction of mesozooplankton carbon released to LPOC, range [0,1]
                     !FRPOCLZ,    &!fraction of mesozooplankton carbon released to RPOC, range [0,1]
                     !FLDONLZ,    &!fraction of mesozooplankton nitrogen released to LDON, range [0,1]
                     !FRDONLZ,    &!fraction of mesozooplankton nitrogen released to RDON, range [0,1]
                     !FLPONLZ,    &!fraction of mesozooplankton nitrogen released to LPON, range [0,1]
                     !FRPONLZ,    &!fraction of mesozooplankton nitrogen released to RPON, range [0,1]
                     !FLDOPLZ,    &!fraction of mesozooplankton phosphorus released to LDOP range [0,1]    
                     !FRDOPLZ,    &!fraction of mesozooplankton phosphorus released to RDOP range [0,1]    
                     !FLPOPLZ,    &!fraction of mesozooplankton phosphorus released to LPOP range [0,1]    
                     !FRPOPLZ,    &!fraction of mesozooplankton phosphorus released to RPOP range [0,1]    
                     !FNH4LZ,    &!fraction of mesozooplankton nitrogen recycled to DIN as NH4, range [0,1]
                     !FPO4LZ,    &!fraction of mesozooplankton phosphorus recycled to DIP as PO4, range [0,1]
                     !FUREASZ,    &!Never used !!!                                              
                     !FUREALZ,    &!Never used!!!
                     !B1ASZ,       &!
                     !B2ASZ,       &!
                     !B3ASZ,       &!
                     !LPOCASZ,    &!
                     !RPOCASZ,    &!
                     !PRASZ,      &!
                     !B1ALZ,       &!
                     !B2ALZ,       &!
                     !B3ALZ,       &!
                     !SZALZ,       &!
                     !LPOCALZ,     &!
                     !RPOCALZ,     &!
                     !PRALZ,       &!
                     !CLSZ,        &!
                     !CLLZ,        &!
                     !RSZ,         &!
                     !RLZ,         &!
                     !RMAXSZ,      &!
                     !RMAXLZ,      &!
                     !BMSZ,        &!
                     !BMLZ,        &!
                     !BMRSZ,       &!
                     !BMRLZ,       &!
                     !MSZ,         &!
                     !MLZ,         &!
                     !PRSZLZ,      &!
                     !GSZ,         &!
                     !GLZ,         &!
                     !ESZ,         &!
                     !ELZ,         &!
                     !RFSZ,        &!
                     !RFLZ,        &!
                     !PRSZ,        &!
                     !PRLZ,        &!
                     !LDOCASZ,     &!
                     !BPRSZ,       &!
                     !BPRLZ,       &!
                     !RDOCASZ,     &!
                     !LDOCSZ,      &!
                     !LPOCSZ,      &!
                     !RPOCSZ,      &!
                     !LDOCLZ,      &!
                     !LPOCLZ,      &!
                     !RPOCLZ,      &!
                     !NH4SZ,       &!
                     !LDONSZ,      &!
                     !LPONSZ,      &!
                     !RPONSZ,      &!
                     !NH4LZ,       &!
                     !LDONLZ,      &!
                     !LPONLZ,      &!
                     !RPONLZ,      &!
                     !PO4SZ,       &!
                     !LDOPSZ,      &!
                     !LPOPSZ,      &!
                     !RPOPSZ,      &!
                     !PO4LZ,       &!
                     !LDOPLZ,      &!
                     !LPOPLZ,      &!
                     !RPOPLZ,      &!
                     !RDOCSZ,      &!
                     !RDONSZ,      &!
                     !RDOPSZ,      &!
                     !RDOCLZ,      &!
                     !RDONLZ,      &!
                     !RDOPLZ,      &!
                     !PIB1SZ,      &!
                     !PIB2SZ,      &!
                     !PIB3SZ,      &!
                     !PIB1LZ,      &!
                     !PIB2LZ,      &!
                     !PIB3LZ,      &!
                     !B1SZ,        &!
                     !B2SZ,        &!
                     !B3SZ,        &!
                     !B1LZ,        &!
                     !B2LZ,        &!
                     !B3LZ,        &!
                     !DOSZ,        &!
                     !DOLZ,        &!
                     !SASZ,        &!
                     !SUSZ,        &!
                     !SALZ,        &!
                     !SULZ,        &!
                     
                     ACLSZ,        &!
                     ACLLZ,        &!
                     ARSZ,        &!
                     ARLZ,       &!
                     ABMSZ,        &!
                     ABMLZ,        &!
                     AMSZ,        &!
                     AMLZ,       &!
                     APRSZLZ,    &!
                     AGSZ,        &!
                     AGLZ,         &!
                     ADOCSZ,    &!        
                     APOCSZ,    &!
                      ADOCLZ,    &!
                     APOCLZ,     &!
                     ANH4SZ,    &!
                     ADONSZ,    &!
                     APONSZ,    &!
                     ANH4LZ,    &!
                     ADONLZ,    &!
                     APONLZ,    &!
                     APO4SZ,    &!
                     ADOPSZ,    &!
                     APOPSZ,    &!
                     APO4LZ,    &!
                     ADOPLZ,    &!
                     APOPLZ,    &!
                     APRSZ,        &!
                     APRLZ,        &!
                     APISZ,        &!
                     APILZ,        &!
                     AB1SZ,     &!
                     AB2SZ,        &!
                     AB3SZ,        &!
                     AB1LZ,        &!
                     AB2LZ,        &!
                     AB3LZ,        &!
                     
                     ADOSZ,        &!
                     ADOLZ,        &!
                     ASASZ,        &!
                     ASUSZ,        &!
                     ASALZ,        &!
                     ASULZ,        &!

                     !FTLZ,        &!
                     !FTSZ,        &!
                     !FTBMSZ,    &!
                     !FTBMLZ,    &!
                     !FTPRSZ,    &!
                     !FTPRLZ,    &!
                     !SZ,        &!
                     !LZ,        &!
                     ZOOP_ALLOC, &!Subroutine                                                                   
                     ZOOP_READ
             
   IMPLICIT NONE
   INTEGER :: J, JF, JG,     &
                JCON,        &    !Index for constituents
                K, L, N, NC_ID2    
   INTEGER :: I, II, ISF, JP, JNP
   INTEGER :: JS1, JS2, JS3
   INTEGER :: YEARSTRT
   CHARACTER(LEN=1024) :: NCFILE
   CHARACTER(LEN=4) :: YEAR
   REAL(SP), ALLOCATABLE :: FTEMP(:)      !,FTEMP2(:)
   REAL(SP), ALLOCATABLE :: RTMP(:,:,:)
   REAL(SP), ALLOCATABLE :: RTMP2(:,:),RTMP21(:,:),RTMP22(:,:),     &
                        RTMP23(:,:),RTMP24(:,:),RTMP25(:,:),RTMP26(:,:)

    !WLong moved this to mod_sed.F
    !REAL(SP), ALLOCATABLE ::    CTEMP_GL(:),    &
    !                        PO4T2TM1S_GL(:),    &
    !                        NH4T2TM1S_GL(:),    &
    !                        NO3T2TM1S_GL(:),    &
    !                         HST2TM1S_GL(:),     &
    !                        CH4T2TM1S_GL(:),    &
    !                         CH41TM1S_GL(:),    &
    !                        SO4T2TM1S_GL(:),     &
    !                         SIT2TM1S_GL(:),    &
    !                         BENSTRTM1S_GL(:),    &
    !                            CPOP_GL(:,:),    &
    !                           CPON_GL(:,:),    &
    !                           CPOC_GL(:,:),    &
    !                             CPOS_GL(:),    &
    !
    !!WLong moved BBM_GL to mod_ba.F
    !                              BBM_GL(:)

!WLong moved this to mod_sf.F
!   REAL(SP), ALLOCATABLE :: SEDTYPE_GL(:,:)  
!   REAL(SP), ALLOCATABLE :: SFEED_GL(:,:)
!WLong moved this to mod_df.F                        
!   REAL(SP), ALLOCATABLE :: DFEEDM1S_GL(:)

   CHARACTER(LEN=120)  :: ERROR_MSG

!************************************************************************
!**                              Inputs                                **
!************************************************************************

!WLong: added UNIT_LINKAGE and wqm_linage.in for input control of linkage
!       to FVCOM
!***** Open Linkage control file (Linkage to FVCOM)

   OPEN(UNIT_LINKAGE,FILE='wqm_linkage.in',STATUS='OLD')

     NstationNum_GL=0
     Nstation=0

     READ(UNIT_LINKAGE, NML=hydro_netcdf, ERR=10)              !read $hydro namelist
     READ(UNIT_LINKAGE, NML=wqm_history,ERR=20)                !read history  output namelist
     READ(UNIT_LINKAGE, NML=wqm_stations,ERR=30)               !read stations output namelist

!---WLong, below is what wqm_linkage.in should look like-----
!$hydro_netcdf
!     hydro_dir='/home/long075/fvcom/outputs/'                 !directory name where the fvcom
!     hydro_prefix='psm_'
!     hydro_suffix='.nc'                                       !extentionname
!     hydro_filenumstart= 1                                    !number in first netcdf input file
!     hydro_Nrec   = 72                                        !number of records in each netcdf file
!     hydro_filenumwidth = 4                                   !number of digits in file numbe following hydro_prefix
!$end
!$wqm_history
!
!    HISFN = 'output/wqm_history.out'
!    t_his_start = 1
!    t_his_end  =  365
!    t_his_dlt  =  0.25
!    HISFN_SPLIT_BYLEVEL = .FALSE.
!
!$end
!$wqm_stations
!    STNFN       = 'output/wqm_station.out'                    !output file name for stations
!    Nstation    = 12                                          !total number of stations
!    NStationLoc = 1 2 3 4 5 6 7 8 9 10 11 12                  !node numbers in fvcom model grid for outputing station timeseries
!    t_stn_start = 1
!    t_stn_end   = 365
!    t_stn_dlt   = 0.25
!$end
!-------------------------------------------------------------

     GOTO 91                                                   !if nothing wrong, goto 91 directly

10   ERROR_MSG='Error when reading hydro, pls check wqm_linkage.in' 
     GOTO 90

20   ERROR_MSG='ERROR when reading wqm_history, pls check wqm_linkage.in'
     GOTO 90

30   ERROR_MSG='ERROR when reading wqm_stations, pls check wqm_linkage.in'
     GOTO 90

90   CLOSE(UNIT_LINKAGE)
     CALL ERROR_CHECK(ERROR_MSG, -1)

91   CLOSE(UNIT_LINKAGE)
     IF(MSR)WRITE(*,*)'Finished reading wqm_linkage.in ...'  

!find out the FVCOM netcdf output file name (first file) 

     NCFILE_PREFIX=TRIM(hydro_prefix)
     NCFILE_SUFFIX=TRIM(hydro_suffix)
     WRITE(FORMAT_STR,'(A4,I1,A1)')'(I0.',hydro_filenumwidth,')'
     NCFILE_NUMBER=''         !empty string for NCFILE_NUMBER
                              !Populate it with the number hydro_filenumstart
     WRITE(NCFILE_NUMBER(1:hydro_filenumwidth),TRIM(FORMAT_STR))hydro_filenumstart
     NCFILE=TRIM(hydro_dir)//TRIM(NCFILE_PREFIX)//TRIM(NCFILE_NUMBER)//TRIM(NCFILE_SUFFIX)

!***** Open control file

   OPEN (CON,FILE=CONFN,STATUS='OLD')

!***** Title cards

   READ (CON,1000)
   READ (CON,1010) (TITLE(J),J=1,6)

!***** Grid setup -- TYKIM -- NOBTY (open boundary numbers)

   READ (CON,1020)  MGL,NGL,KB,NOBTY

!***** Time control -- TYKIM -- DTI (FVCOM internal time steps)

   READ (CON,1031)  TMSTRT,TMEND,YEARSTRT,DTI
   READ (CON,1040)  NDLT
   READ (CON,1030) (DLTD(J),J=1,NDLT)
   READ (CON,1030) (DLTVAL(J),J=1,NDLT)
   READ (CON,1030) (DLTMAX(J),J=1,NDLT)
   READ (CON,1030) (DLTFTN(J),J=1,NDLT)

!TYKIM
! FILGTH IS LENGTH (IN DAYS) OF PERIOD SPANNED BY MET AND OTHER INPUT
!RGL FILGTH is only specified once in the hydro part of the con file??
! AHMDLT IS FVCOM output interval
   READ (CON,1030)  AHMDLT, FILGTH
! RGl moved JDAY setting earlier to allow for correct netcdf file find
   JDAY = TMSTRT
!--SELECT NETCDF FILE TO READ FOR GRID CONSTRUCTION
   IF(MSR)WRITE(YEAR,'(I4.4)') YEARSTRT


!   NCFILE="/home/long075/work/input-netcdf/"//TRIM(CASENAME)//"_0001.nc"


!--DETERMINE NUMBER OF ELEMENTS AND NODES IN THE MODEL
   
   !read global dimensions
   CALL NCD_READ_GRID(NCFILE)
   
   !allocate global variables for geometry related hydrodynamics
   CALL HYDRO_GEOM_ALLOC 
   
   IF(MSR)WRITE(*,*)  '!  # OF NODES            :',MGL
   IF(MSR)WRITE(*,*)  '!  # OF ELEMENTS         :',NGL
   IF(MSR)WRITE(*,*)  '!  # OF SIGMA LEVELS     :',KB     !KBM1
   IF(MSR)WRITE(*,*)
   IF(MSR)WRITE(*,*)  '!      MESH READING      :    FINISHED'

   !Read global grid coordiantes, sigma levels etc 
   CALL NCD_READ_SHAPE(NCFILE)
   
!
!--DECOMPOSE DOMAIN BY ELEMENTS USING METIS
!
#  if defined (MULTIPROCESSOR) 
   ALLOCATE(EL_PID(NGL))  ;  EL_PID = 1
   IF(PAR)CALL DOMDEC(NGL,NPROCS,EL_PID,MSR)
#  endif 

!
!--GENERATE GLOBAL<==>LOCAL ELEMENT/NODE MAPPING
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL GENMAP
#  endif

   !Allocate variables related to local hydrodynamics 
    CALL HYDRO_ALLOC
   
   !Allocate water quality related variables for local domain
    
  ! CALL WQM_ALLOC !Moved after CON file is read
#if defined(ZOOPLANKTON)
   CALL ZOOP_ALLOC 
#endif
   CALL WQMINIT_ALLOC

!
!--MAP OPEN BOUNDARY CONDITION NODES TO LOCAL DOMAIN
!
   CALL BCMAP

!
!--INPUT AND SETUP BOUNDARY FORCING (HEAT/RIVERS/WIND/etc)
!
   CALL BCS_FORCE

!
!--SHIFT GRID/CORIOLIS/BATHYMETRY TO LOCAL DOMAIN
!
   
   CALL TGE_ALLOC      !allocate NV (element to node connectivity)     
   CALL PDOMDEC         !calculate local bathymetry, grid locations etc (using NV)

!
!--SET UP GRID METRICS (FLUX EDGES/CONTROL VOLUMES/ETC)
!
   CALL TRIANGLE_GRID_EDGE      !Set up fluxes and control Volumes
!   write(*,*)'called tri grid edge subrotuine'
  
   CALL CELL_AREA               !Calculate Element and Control Volume Areas

!  write(*,*)'called cell area'
!--SET ISBCE AND ISONB CORRECTLY IN HALO CELLS/NODES
#  if defined (MULTIPROCESSOR)
   ALLOCATE(FTEMP(0:NTLOC)) ; FTEMP = ISBCE
   IF(PAR)CALL EXCHANGE(EC,NTLOC,1,MYID,NPROCS,FTEMP)
   ISBCE = FTEMP
   DEALLOCATE(FTEMP)
   ALLOCATE(FTEMP(0:MTLOC)) ; FTEMP = ISONB
   IF(PAR)CALL EXCHANGE(NC,MTLOC,1,MYID,NPROCS,FTEMP)
   ISONB = FTEMP
   DEALLOCATE(FTEMP)
#  endif

!
!--EXCHANGE SHAPE FACTOR INFORMATION
!
#  if defined (MULTIPROCESSOR)
!   IF(PAR)CALL EXCHANGE(EC,NTLOC,4,MYID,NPROCS,A1U,A2U)
!   IF(PAR)CALL EXCHANGE(EC,NTLOC,3,MYID,NPROCS,AWX,AWY,AW0)
!   IF(PAR)CALL EXCHANGE(EC,NTLOC,1,MYID,NPROCS,ALPHA) 
    IF(PAR)CALL EXCHANGE(EC,NTLOC,1,MYID,NPROCS,ART)
#  endif

!***** Output control

   
   READ (CON,1050)  SNPC,NSNP
   READ (CON,1030) (SNPD(J),J=1,NSNP)
   READ (CON,1030) (SNPF(J),J=1,NSNP)
   READ (CON,1055)  PLTC,QPLTC,SPLTC,SAVPLTC,NPLT
   READ (CON,1030) (PLTD(J),J=1,NPLT)
   READ (CON,1030) (PLTF(J),J=1,NPLT)
   READ (CON,1050)  APLTC,NAPL
   READ (CON,1030) (APLTD(J),J=1,NAPL)
   READ (CON,1030) (APLF(J),J=1,NAPL)
   READ (CON,1050)  TFLC,NTFL
   READ (CON,1030) (TFLD(J),J=1,NTFL)
   READ (CON,1030) (TFLF(J),J=1,NTFL)
   READ (CON,1052)  KFLC,NKFL
   READ (CON,1030) (KFLD(J),J=1,NKFL)
   READ (CON,1030) (KFLF(J),J=1,NKFL)
   READ (CON,1052)  OPLC,NOPL,NOINT
   READ (CON,1030) (OINT(J),J=1,NOINT)
   READ (CON,1030) (OPLD(J),J=1,NOPL)
   READ (CON,1030) (OPLF(J),J=1,NOPL)
   READ (CON,1050)  MBLC,NMBL
   READ (CON,1030) (MBLD(J),J=1,NMBL)
   READ (CON,1030) (MBLF(J),J=1,NMBL)
   READ (CON,1050)  DIAC,NDIA
   READ (CON,1030) (DIAD(J),J=1,NDIA)
   READ (CON,1030) (DIAF(J),J=1,NDIA)
   READ (CON,1050)  RSOC,NRSO,RSIC
   READ (CON,1030) (RSOD(J),J=1,NRSO)
   
   
!***** Hydrodynamic solution scheme

   READ (CON,1060)  SLC, CONSC, TH, MINSTEP

  
!***** Controls

   READ (CON,1071)  SEDC, AUTOC, VBC, BFOC, STLC, ICIC, ICOC, SAVMC, BAOC !WLong added BAOC

   
  
! suspension and deposit feeders and benthic algae

   READ (CON,1072)  SFLC, DFLC, SFLOXC,DFLOXC, BALC, DFOC  !WLong: added DFLOXC, and DFOC 

  
!***** Dead sea case

   READ (CON,1070)  FLC,XYDFC,ZDFC

  
!***** Dispersion

   READ (CON,1030)  XYDF,ZDFMUL,ZDFBCK
   
   

!***** Dispersion type (closure (ON) or constant (OFF))

   READ (CON,1070)  XYDFU

   
!***** Constituent control cards

   READ (CON,1070)  BCC,    S1C,    S2C,    S3C,  MDC,    BFC,  ATMC,   SAVLC
   
   
  
   READ (CON,1080)  REDS1C, REDS1N, REDS1P, REDS2C, REDS2N, REDS2P,     &
                    REDS3C, REDS3N, REDS3P
                    
        
   READ (CON,1080)  REDCBC, REDCBN, REDCBP
   
   
  
   READ (CON,1065)  BNDTC
   
 
   READ (CON,1070) (ACC(JCON),JCON=1,NCP)

   
  
!***** Input filenames

   READ (CON,1020)  NHYDF,NTVDF
   
 
   READ (CON,1090)  MAPFN
   
  
   READ (CON,1090)  GEOFN
   
   
  
   READ (CON,1090)  ICIFN
   

  
   READ (CON,1090)  RSIFN
   
   
   READ (CON,1090)  AGRFN
   
   
  
   READ (CON,1090)  ZOOFN
   
   
   READ (CON,1090)  SFIFN
   
  
   READ (CON,1090)  STLFN
   
   
  
   READ (CON,1090)  MRLFN                            !MNOEL 2/20/93
  
  
   READ (CON,1090)  KEIFN
      ! write(*,*)'KEI file = ',KEIFN
   
   READ (CON,1090) (METFN(J),J=1,NTVDF)
   IF(MSR)THEN
    WRITE(*,*)(METFN(J),J=1,NTVDF)
   ENDIF
   READ (CON,1090) (S1FN(J),J=1,NTVDF)
   READ (CON,1090) (S2FN(J),J=1,NTVDF)
   READ (CON,1090) (S3FN(J),J=1,NTVDF)
   READ (CON,1090) (ATMFN(J),J=1,NTVDF)
   READ (CON,1090) (SVIFN(J),J=1,NTVDF)
   READ (CON,1090) (BFIFN(J),J=1,NTVDF)
   READ (CON,1090) (CBCFN(J),J=1,NTVDF)      !JLM 2/23/01
   READ (CON,1090) BAIFN                      !WLong Benthic algae input file
   READ (CON,1090) DFIFN                      !WLong Deposition feeder input file name
!***** Output filenames

   READ (CON,1090)  ICOFN                   !MNOEL  2-5-93
   READ (CON,1090)  SNPFN
   READ (CON,1090)  RSOFN
   READ (CON,1090)  PLTFN
   READ (CON,1090)  APLFN
   READ (CON,1090)  DIAFN
   READ (CON,1090)  TFLFN
   READ (CON,1090)  KFLFN
   READ (CON,1090)  OPLFN
   READ (CON,1090)  MBLFN
   READ (CON,1090)  ALOFN
   READ (CON,1090)  ZFOFN
   READ (CON,1090)  BFOFN
   READ (CON,1090)  SVOFN
   READ (CON,1090)  SFOFN
   READ (CON,1090)  BAOFN  !WLong read benthic algae output file name
   READ (CON,1090)  DFOFN  !WLong read Deposition feeder output file name
   CLOSE (CON)
 ! write(*,*)'Read the Con FIle'  ! B Clark Debug
!***** Initialize I/O logical control variables

   BOUNDARY_CONC    = BCC == ' ON'
   MODIFY_ICONC     = MDC == ' ON'
   SOURCE_ONE       = S1C == ' ON'           !MNOEL   1-25-93
   BENTHIC_FLUXES   = BFC == ' ON'
   SOURCE_TWO       = S2C == ' ON'           !MNOEL   1-25-93
   SOURCE_THR       = S3C == ' ON'           !MNOEL   1-25-93
   ATMOS_LOADS      = ATMC == ' ON'
   SAV_LOADS        = SAVLC == ' ON'
   SAV_CALC         = SAVMC == ' ON'

   !Wen Long, make sure SAV_LOADS and SAV_CALC are not on at the same time
   !we can have only one of them or neither of them
   IF(SAV_LOADS .AND. SAV_CALC)THEN
       ERROR_MSG='Oops, SAV_LOADS and SAV_CALC can''t be on at same time!!'
       ERROR_MSG=TRIM(ERROR_MSG)//'Please check SAVMC and SAVLC in input file'
       CALL ERROR_CHECK(ERROR_MSG,-1)
   ENDIF
   
   SETTLING         = STLC == ' ON'
   ICOND_OUT        = ICOC == ' ON'              !MNOEL 2-5-93
   SEDIMENT_CALC    = SEDC == ' ON'
   
   !PRINTE error if BOTH SEDC and BFC are ' ON'
   
   IF(MSR)THEN
      IF(BENTHIC_FLUXES .AND. SEDIMENT_CALC)THEN !fix the format of this error message later
        WRITE(*,*)'Oops, should not have BFC and SEDC ON at same time'
        CALL PSTOP
      ENDIF
   ENDIF
   
   BALGAE_CALC      = BALC == ' ON'                 !WLong benthic algae calculation
   KINETIC_FLUXES   = KFLC == ' ON'
   LIGHT_EXTINCTION = ((ACC(4) == ' ON').OR.(ACC(5) == ' ON')      &
                      .OR.(ACC(6) == ' ON'))
   ZOO_CALC         = ((ACC(7) == ' ON').OR.(ACC(8) == ' ON'))
   UNI_ICON_IN      = ICIC == ' UNIFORM'              !MNOEL 2-5-93
   BIN_ICON_IN      = ICIC == '  BINARY'              !MNOEL 2-5-93
   RES_ICON_IN      = ICIC == ' RESTART'
! deposit and suspension feeder control

   DFEEDER          = DFLC == ' ON'
   SFEEDER          = SFLC == ' ON'
   HYPOXFX_SF          = SFLOXC == ' ON'
   HYPOXFX_DF          = DFLOXC == ' ON'

!***** Active Constituents

   NAC = 0
   DO JCON=1,NCP
     IF (ACC(JCON) == ' ON') THEN
       NAC     = NAC+1
       AC(NAC) = JCON
     ENDIF
   ENDDO
  !  write(*,*) ' NAC = ', NAC
  !  write(*,*)'Active consituents are'
  !  write(*,*),AC

     VOLUME_BALANCE   = VBC == ' ON'
     FLOW             = FLC == ' ON'
     PLOTS            = PLTC == ' ON'
     MASS_BALANCE     = MBLC == ' ON'
     OXYGEN_PLOTS     = OPLC == ' ON'
     SNAPSHOTS        = SNPC == ' ON'
     BENTHIC_OUTPUT   = BFOC == ' ON'
     DFEEDER_OUTPUT   = DFOC == ' ON' !Wen Long, deposition feeder output
     BA_OUTPUT        = BAOC == ' ON' !Wen Long, benthic algae output
     TRANSPORT_FLUXES = TFLC == ' ON'.OR.MBLC == ' ON'
     RESTART_OUT      = RSOC == ' ON'
     RESTART_IN       = RSIC == ' ON'
     DIAGNOSTICS      = DIAC == ' ON'
     Z_DIFFUSION      = ZDFC == ' ON'
     AVERAGE_PLOTS    = APLTC == ' ON'
     QUALITY_DIAG     = QPLTC == ' ON'
     AUTO_STEPPING    = AUTOC == ' ON'
     XY_DIFFUSION     = XYDFC == ' ON'
     UPWIND           = SLC == '  UPWIND'
     QUICKEST         = SLC == 'QUICKEST'
     STEP_BOUNDARY    = BNDTC == '    STEP'
     CONSERVE_MASS    = CONSC == '    MASS'
     SEDIMENT_DIAG    = SPLTC == ' ON'.AND.(SEDIMENT_CALC .OR. BENTHIC_FLUXES)
     SAV_PLOTS        = SAVPLTC == ' ON'.AND.(SAV_LOADS.OR.SAV_CALC)
     END_RUN          = .FALSE.

     IF (ACC(1) == ' ON')  TEMPERATURE_CALC = .TRUE.
     IF (ACC(3) == ' ON')  SOLIDS_CALC      = .TRUE.
     IF ((ACC(4) == ' ON').OR.(ACC(5) == ' ON').OR.(ACC(6) == ' ON'))  &
                           ALGAE_CALC   = .TRUE.
     IF (ACC(9) == ' ON')  CARBON_CALC      = .TRUE.
     IF (ACC(13) == ' ON') NITROGEN_CALC    = .TRUE.
     IF (ACC(20) == ' ON') PHOSPHORUS_CALC  = .TRUE.
     IF (ACC(26) == ' ON') COD_CALC         = .TRUE.
     IF (ACC(27) == ' ON') OXYGEN_CALC      = .TRUE.
     IF (ACC(29) == ' ON') SILICA_CALC      = .TRUE.
     IF (ACC(25) == ' ON') PIP_CALC         = .TRUE.
if(MSR) write(*,*)'About to allocate the vars'
   !Have to allocate after LIGHT_EXTINCTION is set to True or False
    CALL WQM_ALLOC

 IF(MSR) write(*,*) 'vars are Allocated'
!***** Open remaining input files

   HYDPTR = 1
   METPTR = 1
   CBCPTR = 1
   S1PTR  = 1
   S2PTR  = 1
   S3PTR  = 1
   BFIPTR = 1
   BAOPTR = 1
   KEIPTR = 1
   ATMPTR = 1
   SAVPTR = 1

   OPEN (MET,FILE=METFN(METPTR),STATUS='OLD')

  IF(MSR) write(*,*)'Opened the met file'

   IF (SETTLING)         OPEN (STL,FILE=STLFN,STATUS='OLD')
!   IF (LIGHT_EXTINCTION) OPEN (KEI,FILE=KEIFN,STATUS='OLD')    !WLong moved to downstairs
   IF (ALGAE_CALC) OPEN (AGR,FILE=AGRFN,STATUS='OLD')            !Wen Long LIGHT_EXTINCTION -> ALGAE_CALC
   IF (SEDIMENT_CALC.OR.BENTHIC_FLUXES)                      &
                         OPEN (BFI,FILE=BFIFN(BFIPTR),STATUS='OLD')
   
   IF (BALGAE_CALC)                                           &
                         OPEN (BAI,FILE=BAIFN,STATUS='OLD')
   
   IF (BOUNDARY_CONC)    OPEN (CBC,FILE=CBCFN(CBCPTR),STATUS='OLD') 
   IF (SOURCE_ONE)       OPEN (S1, FILE=S1FN(S1PTR),  STATUS='OLD')
   IF (SOURCE_TWO)       OPEN (S2, FILE=S2FN(S2PTR),  STATUS='OLD')
   IF (SOURCE_THR)       OPEN (S3, FILE=S3FN(S3PTR),  STATUS='OLD')
   IF (ATMOS_LOADS)      OPEN (ATM,FILE=ATMFN(ATMPTR),STATUS='OLD')

     IF(Nstation > NstationMax)THEN
       ERROR_MSG=''
       WRITE(ERROR_MSG,'(I5)')NstationMAX
       ERROR_MSG='Nstation exceeding maximum value allowed:'//TRIM(ERROR_MSG)
       CALL ERROR_CHECK(ERROR_MSG,-1)
     ENDIF

     IF(Nstation <=0 .OR. MAXVAL(NstationNum_GL(1:Nstation),1) > MGL  &
                     .OR. MINVAL(NstationNum_GL(1:Nstation),1) < 1 )THEN
       ERROR_MSG='Error in wqm_stations section, please check wqm_linkage.in'
       CALL ERROR_CHECK(ERROR_MSG,-1)
     ENDIF

     IF(t_stn_dlt<0.0)THEN
       ERROR_MSG='Error wqm_stations section, t_stn_dlt must be > 0. Pls check wqm_linkage.in'
       CALL ERROR_CHECK(ERROR_MSG,-1)
     ENDIF

     IF(t_his_dlt<=0.0)THEN
       ERROR_MSG='Error in wqm_history section, t_his_dlt must be > 0. Pls check wqm_linkage.in'
       CALL ERROR_CHECK(ERROR_MSG,-1)
     ENDIF


! WLong: 
! Read two time steps in FVCOM netcdf output files. Note that this may have to be moved
! to after the control file reading is finished. For model restart, one would 
! need to find out what time record to read based on restarting time, as opposed to 
! reading the first two records of all files combined.
!

! Tykim FVCOM netcdf outfile is daily base - RGL -test change to hourly/semihourly for coarse grid  
! coarsegrid = read every timestep, open new file every 2 simulation hours
!  DAY = DAY + 1
   IJDAY=AINT(JDAY)
   IF(MSR)WRITE(*,*)'IJDAY, JDAY', IJDAY,JDAY

!WLong: Note here R_HOUR was never initiallized, and we should have set R_HOUR to 0 at beginning 
!       Also netcdf files do not necessary only have length of one hour.  Hence calling the
!       file numbers by R_HOUR is misleading. We should change these to IFNC etc for file number
!       counting and NTRECNC for time record counting within each netcdf file 


   IFNC = 0       !initialize indices for file
   NTRECNC = 0    !indices for time record in each file
   NTHYDRO = 0    !indices for overall hydrodyanmics 

! The first record

   IFNC = 1       !read from the first file (WL: for restart, this has to be different)
   NTRECNC =1     !read from the first record for the current file
   NTHYDRO = 1    !the first record over all

    
   NCFILE_NUMBER=''         !empty string for NCFILE_NUMBER
                            !Populate it with the number with hydro_filenumstart

   WRITE(NCFILE_NUMBER(1:hydro_filenumwidth),TRIM(FORMAT_STR))(IFNC+hydro_filenumstart-1)
   NCFILE=TRIM(hydro_dir)//TRIM(NCFILE_PREFIX)//TRIM(NCFILE_NUMBER)//TRIM(NCFILE_SUFFIX)

   CALL NCD_READ_OPEN(NCFILE,UNC1,VNC1,WTSNC1,UARD_OBCNNC1,XFLUX_OBCNC1,DTFANC1,&
        KHNC1,ELNC1,TNC1,SNC1,NTRECNC)  !reading NTRECNC'th record in NCFILE

   UUT=UNC1
   VVT=VNC1
   WTST=WTSNC1
   UARD_OBCNT=UARD_OBCNNC1
   XFLUX_OBCT=XFLUX_OBCNC1
   DTFAT=DTFANC1
   ET=ELNC1
   DT=H+ET
  DO I=1,NTLOC
    DT1(I)=DT(NV(I,1))+DT(NV(I,2))+DT(NV(I,3))
    DT1(I)=DT1(I)/3.0
  ENDDO
 
!The seocnd record 

    NTRECNC=NTRECNC+1  !point to next record
    NTHYDRO=NTHYDRO+1
   
    IF(NTRECNC .gt. hydro_Nrec) THEN

       NTRECNC =1                   !reset time index for file to be 1
       IFNC = IFNC +1               !A new file 

      NCFILE_NUMBER=''
      WRITE(NCFILE_NUMBER(1:hydro_filenumwidth),TRIM(FORMAT_STR))(IFNC+hydro_filenumstart-1)
      NCFILE=TRIM(hydro_dir)//TRIM(NCFILE_PREFIX)//TRIM(NCFILE_NUMBER)//TRIM(NCFILE_SUFFIX)

      IF(MSR) WRITE(*,*) TRIM(NCFILE)

     CALL NCD_READ_OPEN(NCFILE,UNC2,VNC2,WTSNC2,UARD_OBCNNC2,XFLUX_OBCNC2,DTFANC2,&
          KHNC2,ELNC2,TNC2,SNC2,NTRECNC)

   ELSE

     CALL NCD_READ(NCFILE,UNC2,VNC2,WTSNC2,UARD_OBCNNC2,XFLUX_OBCNC2,DTFANC2,KHNC2,ELNC2,TNC2,SNC2,NTRECNC)

   ENDIF

    !***** Suspension feeders
    IF(SFEEDER) THEN
        CALL SF_ALLOC
    ENDIF
     
 
    IF (SFEEDER) THEN
     
       CALL SF_READ
     
        ! initialize arrays that are passed in case susp. feeder model not used
        CALL SF_INIT
        
    ENDIF

    IF(DFEEDER)THEN
        CALL DF_ALLOC
        CALL DF_READ
        CALL DF_INIT
    ENDIF


#if defined(SEDIMENT_DOM)   

	    CALL SED_DOM_ALLOCATE_VARS
		  CALL SED_DOM_INITIALIZE
      CALL SED_DOM_INPUT
      IF(MSR)THEN
      WRITE(*,*)'READ DOM INPUT'
      ENDIF

#endif		 
 
IF(MSR)WRITE(*,*)'Allocated all WQM vars and Got hydrodyamic data'
	CALL WC_DOM_ALLOCATE
	CALL WC_DOM_INPUT
	
IF(MSR)WRITE(*,*)'Allocated all DOM vars'

!************************************************************************
!**                   Parameters for SAV Submodel                      **
!************************************************************************

   IF (SAV_CALC) THEN
       !read SVIFN file and find maximum number of SAV species - NSAVM
   
       OPEN (SVI,FILE=SVIFN,STATUS='OLD')   
   
       READ(SVI,'(///)')    !jump 4th lines (first skip 3 lines and then read nothing from 4th line)   
       READ(SVI,'(A72)') (TITLE(J),J=1,6)  !read six lines
       READ(SVI, '(//8X,I8)') NSAVM     !go down two lines and read the third line
       CLOSE(SVI)
       
       !WLong, it's probably beter to ut NSAVM in CON file and  then skip the upper SVI reading thing
       
       !Use the maximum number of SAV species to allocate SAV arrays
       CALL SAV_ALLOC
       !read SAV parameters
       OPEN (SVI,FILE=SVIFN,STATUS='OLD')          
       CALL SAV_READ       
       CLOSE(SVI)
   ENDIF
   
   !Wen Long: Open SVI for tvds.F to read if SAV_LOADS is ON
   !          Note that SAV_LOADS and SAV_CALC can't be ON at the same time!
   IF (SAV_LOADS) THEN
      CALL SAV_LOADS_ALLOC
      OPEN (SVI,FILE=SVIFN(SAVPTR),STATUS='OLD')
   ENDIF
   
!******* Initial conditions 

!******* Uniform constituent initial concentrations

     IF (UNI_ICON_IN) THEN                              !MNOEL 2-5-93 
       OPEN (ICI,FILE=ICIFN,STATUS='OLD')

! KRG 22-JUNE-2010, the above now in parallel, KURT GLAESEMANN
       ALLOCATE(RTMP(0:MGL,KBM1,NCP));   RTMP = 0.0  !WLong, note index starts from zero
       DO I=1,MGL
         DO K=1,KBM1

           read(ICI,*)(RTMP(I,K,JCON),JCON=1,NCP)

         ENDDO
       ENDDO
       
       IF(SERIAL)C1=RTMP
       
#      if defined (MULTIPROCESSOR)
       IF(PAR) THEN
         DO JCON = 1,NCP
           DO K=1,KBM1
             DO I=1,MLOC
               C1(I,K,JCON) = RTMP(NGID(I),K,JCON)
              !  C2(I,K,JCON) = RTMP(NGID(I),K,JCON)
             ENDDO
             DO I=1,NHN
               C1(I+MLOC,K,JCON) = RTMP(HN_LST(I),K,JCON)
              !  C2(I+MLOC,K,JCON) = RTMP(HN_LST(I),K,JCON)
             ENDDO
           ENDDO
         ENDDO
       ENDIF
#      endif
       DEALLOCATE(RTMP)
       C2 = C1
      ! !  C1MIN = 1.E10
      ! !  C1MAX = 0.
! END KRG



!
!WLong Note below we are still using constant IC for Sediment flux model
!      suspention feeder model, deposition feeder model and sav model
!
       
       !LB: We read this case only if STEADY_STATE_SED_IC=FALSE (ie, if not initializing sediments with steady state)

       !IF (.NOT. STEADY_STATE_SED_IC .AND. SEDIMENT_CALC ) THEN  !--> unfortunately, this option doesn't work because STEADY_STATE_SED_IC is F at this time. 
                                                                  !    It gets updated to T if SSTATEIC=1 when calling SED_READ (line 4138 of this subroutine)
    !   IF (SEDIMENT_CALC ) THEN
             
    !     READ (ICI,1030) CTEMPI
         
    !     READ (ICI,1030) (CPOPI(JG),JG=1,3)
    !     READ (ICI,1030) (CPONI(JG),JG=1,3)
    !     READ (ICI,1030) (CPOCI(JG),JG=1,3)
    !     READ (ICI,1030) CPOSI
         
    !     READ (ICI,1030) PO41TI,PO4T2I    !WLong changed input file reading sequence here
    !     READ (ICI,1030) NH41TI,NH4T2I
    !     READ (ICI,1030) NO31TI,NO3T2I
    !     READ (ICI,1030) CH41TI,CH4T2I         
    !     READ (ICI,1030) HS1TI ,HST2I
    !     READ (ICI,1030) SO41TI,SO4T2I
    !     READ (ICI,1030) SI1TI,SIT2I
    !     READ (ICI,1030) BENSTI, SODI  !Wlong input initial SOD          
    !     IF(BALGAE_CALC)THEN
    !        READ (ICI,1030) BBMI       !benthic algae initial condition
    !     ENDIF
         
    !   ENDIF

    !   IF (DFEEDER) READ (ICI,1030) DFEEDI
    !   IF (SFEEDER) READ (ICI,1030) (SFEEDI(N), N=1,NSPECIES) 
    !   IF (SAV_CALC) READ (ICI,1030) LEAFI, STEMI, ROOTI, TUBERI, EPI

       CLOSE(ICI)   

!RGL changed binary initial condition file to restart file
!******* RESTART constituent initial concentrations

     ELSE IF (RES_ICON_IN) THEN                              !MNOEL 2-5-93
       OPEN (RSI,FILE=RSIFN,STATUS='OLD',FORM='UNFORMATTED') 
       READ (RSI) OLDTITLE

       ALLOCATE(RTMP(0:MGL,KBM1,NCP));   RTMP = 0.0            !WLong: note index from zero
       READ (RSI) (((RTMP(I,K,JCON),I=1,MGL),K=1,KBM1),JCON=1,NCP)

       IF(SERIAL) C1 = RTMP
#      if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         DO JCON = 1,NCP
            DO K=1,KBM1
                DO I=1,MLOC
                    C1(I,K,JCON) = RTMP(NGID(I),K,JCON)
                ENDDO
                DO I=1,NHN
                    C1(I+MLOC,K,JCON) = RTMP(HN_LST(I),K,JCON)
                ENDDO
            ENDDO
            ENDDO
       ENDIF
#      endif             
       DEALLOCATE(RTMP)
    
       DO JCON=1,NAC
         DO K=1,KBM1
           DO I=1,MLOC
             C1(I,K,AC(JCON)) = MAX(C1(I,K,AC(JCON)),0.0)     
             C2(I,K,AC(JCON)) = C1(I,K,AC(JCON))
       ENDDO  
         ENDDO
       ENDDO    

       IF (SEDIMENT_CALC) THEN

     
         READ (RSI) (CTEMP_GL(I),I=1,MGL),                                 &
                ((CPOP_GL(I,J),I=1,MGL),J=1,3),                         &
                ((CPON_GL(I,J),I=1,MGL),J=1,3),                         &
                ((CPOC_GL(I,J),I=1,MGL),J=1,3),                         &
                (CPOS_GL(I),I=1,MGL),                                   &
                (PO4T2TM1S_GL(I),I=1,MGL), (NH4T2TM1S_GL(I),I=1,MGL),   &
                (NO3T2TM1S_GL(I),I=1,MGL), (HST2TM1S_GL(I),I=1,MGL),    &
                (CH4T2TM1S_GL(I),I=1,MGL), (CH41TM1S_GL(I),I=1,MGL),    &
                (SO4T2TM1S_GL(I),I=1,MGL), (SIT2TM1S_GL(I),I=1,MGL),    &
                (BENSTRTM1S_GL(I),I=1,MGL), (SODTM1S_GL(I),I=1,MGL)
                
        IF(BALGAE_CALC)THEN
            READ(RSI)(BBM_GL(I),I=1,MGL)
        ENDIF
          
        IF(SERIAL)THEN 
        
            CTEMP     = CTEMP_GL    
            CPOP      = CPOP_GL     
            CPON      = CPON_GL     
            CPOC      = CPOC_GL     
            CPOS      = CPOS_GL     
            
            PO4T2TM1S = PO4T2TM1S_GL
            NH4T2TM1S = NH4T2TM1S_GL
            NO3T2TM1S = NO3T2TM1S_GL
            HST2TM1S  = HST2TM1S_GL 
            CH4T2TM1S = CH4T2TM1S_GL
            CH41TM1S  = CH41TM1S_GL 
            SO4T2TM1S = SO4T2TM1S_GL
            SIT2TM1S  = SIT2TM1S_GL 
            BENSTRTM1S  = BENSTRTM1S_GL 
            SODTM1S     = SODTM1S_GL
            IF(BALGAE_CALC)THEN
                BBM       = BBM_GL      
            ENDIF
        ENDIF  

# if defined (MULTIPROCESSOR)
        IF(PAR)THEN
           DO I=1,MLOC
                CTEMP(I)      = CTEMP_GL(NGID(I))    
                CPOP(I,:)     = CPOP_GL(NGID(I),:)     
                CPON(I,:)     = CPON_GL(NGID(I),:)     
                CPOC(I,:)     = CPOC_GL(NGID(I),:)     
                CPOS(I)       = CPOS_GL(NGID(I))     
                PO4T2TM1S(I)  = PO4T2TM1S_GL(NGID(I))
                NH4T2TM1S(I)  = NH4T2TM1S_GL(NGID(I))
                NO3T2TM1S(I)  = NO3T2TM1S_GL(NGID(I))
                HST2TM1S(I)   = HST2TM1S_GL(NGID(I))
                CH4T2TM1S(I)  = CH4T2TM1S_GL(NGID(I))
                CH41TM1S(I)   = CH41TM1S_GL(NGID(I)) 
                SO4T2TM1S(I)  = SO4T2TM1S_GL(NGID(I))
                SIT2TM1S(I)   = SIT2TM1S_GL(NGID(I)) 
                BENSTRTM1S(I) = BENSTRTM1S_GL(NGID(I)) 
                SODTM1S(I)    = SODTM1S_GL(NGID(I))
                IF(BALGAE_CALC)THEN
                    BBM(I)       = BBM_GL(NGID(I))      
                ENDIF
           ENDDO
           DO I=1,NHN
                CTEMP(I+MLOC)     = CTEMP_GL(HN_LST(I))    
                CPOP(I+MLOC,:)    = CPOP_GL(HN_LST(I),:)     
                CPON(I+MLOC,:)    = CPON_GL(HN_LST(I),:)     
                CPOC(I+MLOC,:)    = CPOC_GL(HN_LST(I),:)     
                CPOS(I+MLOC)      = CPOS_GL(HN_LST(I))     
                PO4T2TM1S(I+MLOC) = PO4T2TM1S_GL(HN_LST(I))
                NH4T2TM1S(I+MLOC) = NH4T2TM1S_GL(HN_LST(I))
                NO3T2TM1S(I+MLOC) = NO3T2TM1S_GL(HN_LST(I))
                HST2TM1S(I+MLOC)  = HST2TM1S_GL(HN_LST(I)) 
                CH4T2TM1S(I+MLOC) = CH4T2TM1S_GL(HN_LST(I))
                CH41TM1S(I+MLOC)  = CH41TM1S_GL(HN_LST(I)) 
                SO4T2TM1S(I+MLOC) = SO4T2TM1S_GL(HN_LST(I))
                SIT2TM1S(I+MLOC)  = SIT2TM1S_GL(HN_LST(I)) 
                BENSTRTM1S(I+MLOC)= BENSTRTM1S_GL(HN_LST(I)) 
                SODTM1S(I+MLOC)  = SODTM1S_GL(HN_LST(I))
                IF(BALGAE_CALC)THEN
                    BBM(I+MLOC)       = BBM_GL(HN_LST(I))
                ENDIF
           ENDDO
         ENDIF
# endif

          
     ENDIF

     IF (DFEEDER)THEN
!WLong moved this to mod_df.F     
!Wlong   ALLOCATE(DFEEDM1S_GL(MGL));      DFEEDM1S_GL = 0.0
         READ (RSI) (DFEEDM1S_GL(I),I=1,MGL)
         IF(SERIAL) DFEEDM1S = DFEEDM1S_GL
#        if defined (MULTIPROCESSOR)
         IF(PAR)THEN
           DO I=1,MLOC
             DFEEDM1S(I) = DFEEDM1S_GL(NGID(I))
           ENDDO
           DO I=1,NHN
             DFEEDM1S(I+MLOC) = DFEEDM1S_GL(HN_LST(I))
           ENDDO
         ENDIF
#        endif
!WLong         DEALLOCATE(DFEEDM1S_GL)
       ENDIF
       
       IF (SFEEDER)THEN
!         ALLOCATE(SFEED_GL(MGL,NSSFP));    SFEED_GL = 0.0  !Wlong, I think dimension NCP is wrong here
!                                                               !Now takaen care of in mod_sf sf_alloc()
         READ (RSI) ((SFEED_GL(I,K),I=1,MGL),K=1,NSPECIES)
         IF(SERIAL) SFEED = SFEED_GL
#        if defined (MULTIPROCESSOR)
         IF(PAR)THEN
           DO K = 1,NSPECIES
             DO I = 1,MLOC
               SFEED(I,K) = SFEED_GL(NGID(I),K)
             ENDDO
             DO I=1,NHN
               SFEED(I+MLOC,K) = SFEED_GL(HN_LST(I),K)
             ENDDO
           ENDDO
         ENDIF
#        endif
!WLong         DEALLOCATE(SFEED_GL)  !now taken care of in mod_sf.F
       ENDIF
       CLOSE (RSI)
       
!       CLOSE (ICI)   !Wen Logn, removed CLOSE(ICI) HERE

     ELSE
       WRITE(*,*) 'initial conditions file specified incorrectly'
       STOP
     ENDIF

!***** Mineralization rates

     OPEN (MRL,FILE=MRLFN,STATUS='OLD')                        !MNOEL 2-20-93
     READ(MRL,*)
     READ(MRL,*)
     READ(MRL,*)
     READ(MRL,*)
      
!***** Spatially-invariant kinetics parameters

     READ (MRL,*)  KHONT,  KHNNT,  KHOCOD, KHODOC, KHNDN
     if(MSR)THEN
     write (*,*)  'KHONT,  KHNNT,  KHOCOD, KHODOC, KHNDN'
     write (*,*)  KHONT,  KHNNT,  KHOCOD, KHODOC, KHNDN
    endif
     READ(MRL,*)
     READ (MRL,*)  AOCR,   AONT,   KHCOAG
     if(MSR)THEN
     write(*,*)  'AOCR,   AONT,   KHCOAG'
     write(*,*)  AOCR,   AONT,   KHCOAG
    endif
     READ(MRL,*)
     READ (MRL,*)  TRCOD,  TRMNL,  TRHDR,  TRSUA
     if(MSR)THEN
     write (*,*)  'TRCOD,  TRMNL,  TRHDR,  TRSUA'
     write (*,*)  TRCOD,  TRMNL,  TRHDR,  TRSUA
    endif
     READ(MRL,*)
     READ (MRL,*)  KTCOD,  KTMNL,  KTHDR,  KTSUA
     if(MSR)THEN
     write (*,*)  'KTCOD,  KTMNL,  KTHDR,  KTSUA'
     write (*,*)  KTCOD,  KTMNL,  KTHDR,  KTSUA
    endif
     READ(MRL,*)
     READ (MRL,*)  KTNT1,  KTNT2,  TMNT
     if(MSR)THEN
     write (*,*)  'KTNT1,  KTNT2,  TMNT'
     write (*,*)  KTNT1,  KTNT2,  TMNT
    endif
     READ(MRL,*)
     READ (MRL,*)  KADPO4, KADSA
     if(MSR)THEN
     write (*,*)  'KADPO4, KADSA'
     write (*,*)  KADPO4, KADSA
    endif
     READ(MRL,*)
     READ (MRL,*)  AANOX,  ANDC
     if(MSR)THEN
     write (*,*)  'AANOX,  ANDC'
     write (*,*)  AANOX,  ANDC
    endif
     READ(MRL,*)                             !MNOEL  2-20-93
     READ (MRL,*)  AREAR,  BREAR,  CREAR
     if(MSR)THEN
     write(*,*)  'AREAR,  BREAR,  CREAR'
     write(*,*)  AREAR,  BREAR,  CREAR
    endif
     READ(MRL,*)
     READ(MRL,*)KLPC, KRPC
     if(MSR)THEN
     write(*,*) 'KLPC, KRPC'
     write(*,*) KLPC, KRPC
    endif
     READ(MRL,*)
     READ(MRL,*)KLPN, KRPN
     if(MSR)THEN
     write(*,*)'KLPN, KRPN'
     write(*,*)KLPN, KRPN
    endif
     READ(MRL,*)
     READ(MRL,*)KLPP, KRPP
     if(MSR)THEN
     write(*,*)'KLPP, KRPP'
     write(*,*)KLPP, KRPP
    endif
     READ(MRL,*)
     READ(MRL,*)KSUA
     if(MSR)THEN
     write(*,*)'KSUA'
     write(*,*) KSUA
    endif
     READ(MRL,*)
     READ(MRL,*)KCOD 
     if(MSR)THEN
     write(*,*)'KCOD' 
     write(*,*) KCOD 
    endif
     READ(MRL,*)
     READ(MRL,*)KDCALG, KLCALG  
     if(MSR)THEN
     write(*,*) 'KDCALG, KLCALG' 
     write(*,*) KDCALG, KLCALG 
    endif  
     READ(MRL,*)
     READ(MRL,*)KRCOAG2in, KRCOAG3in 
     if(MSR)THEN
     write(*,*)'KRCOAG2in, KRCOAG3in'
     write(*,*)KRCOAG2in, KRCOAG3in  
    endif
     READ(MRL,*)
     READ(MRL,*)KDNALG, KLNALG 
     if(MSR)THEN
     write(*,*)'KDNALG, KLNALG' 
     write(*,*)KDNALG, KLNALG 
     endif
     READ(MRL,*)
     READ(MRL,*)KDPALG, KLPALG   
     if(MSR)THEN
     write(*,*)'KDPALG, KLPALG'
     write(*,*)KDPALG, KLPALG   
    endif

     READ(MRL,*)
     READ(MRL,*)NTM      
     if(MSR)THEN
     write(*,*) 'NTM'    
     write(*,*) NTM    
     endif
     CLOSE (MRL)

!***** Boundary concentrations

     IF (BOUNDARY_CONC) THEN
       READ (CBC,1100)
       READ (CBC,1020) (NCB(JCON),JCON=1,NCP)
       READ (CBC,1100)
     ENDIF

!***** Sources One
!vjp  OPEN(UNIT=130,file='NPS96.OUT',status='unknown') 
     IF (SOURCE_ONE) THEN                         !MNOEL   1-25-93
       READ (S1,1100)
       READ (S1,1020) (S1LN(JCON),JCON=1,NCP)
       DO JCON=1,NCP
         READ (S1,1020) (S1LB(JP,JCON),JP=1,S1LN(JCON))
       ENDDO
       READ (S1,1100)
     ENDIF

!***** Sources Two

     IF (SOURCE_TWO) THEN                         !MNOEL   1-25-93
       READ (S2,1100)
       READ (S2,1020) (S2LN(JCON),JCON=1,NCP)
       DO JCON=1,NCP
         READ (S2,1020) (S2LB(JNP,JCON),JNP=1,S2LN(JCON))
       ENDDO
       READ (S2,1100)
     ENDIF

!***** Sources Three

     IF (SOURCE_THR) THEN                         
       READ (S3,1100)
       READ (S3,1020) (S3LN(JCON),JCON=1,NCP)
       DO JCON=1,NCP
         READ (S3,1020) (S3LB(JNP,JCON),JNP=1,S3LN(JCON))
       ENDDO
       READ (S3,1100)
     ENDIF

!***** Light extinction and algal growth

     IF (LIGHT_EXTINCTION) THEN
        OPEN (KEI,FILE=KEIFN,STATUS='OLD')     
         !   CALL OWQ_ALLOC
            CALL OWQ_READ
        CLOSE(KEI)
     ENDIF

!***** Algae

     IF (ALGAE_CALC) THEN    !WLong corrected this IF , LIGHT_EXTINCTION  -> ALGAE_CALC
       CALL ALG_READ
     ENDIF

!***** Zooplankton 
#if defined(ZOOPLANKTON)
     IF (ZOO_CALC) CALL ZOOP_READ
#endif
!***** Settling rates

     IF (SETTLING) THEN
      READ(STL,*)
      READ(STL,*)
      READ(STL,*)
      READ(STL,*)

         READ(STL,*) WSS,WSL,WSR,WS1,WS2,WS3,WSU,WSSHI,WSSLIM,SedOriginal,SedHyp,Ksed,S_vmax

         if(MSR)THEN
             write(*,*)'WSS,WSL,WSR,WS1,WS2,WS3,WSU,WSSHI,WSSLIM,SedOriginal,SedHyp,Ksed,S_vmax'
             write(*,*)WSS,WSL,WSR,WS1,WS2,WS3,WSU,WSSHI,WSSLIM,SedOriginal,SedHyp,Ksed,S_vmax
         ENDIF

       CLOSE (STL)
     ENDIF

!***** Atmospheric loads

     IF (ATMOS_LOADS) READ (ATM,1000)

     
!***** Benthic fluxes
     IF (BENTHIC_FLUXES) THEN
!       WRITE(*,*)'reading data from BFI file'
       READ (BFI,11040)                                       !skip 3 title lines 
                                                              !(the format skips two lines, the read itself skips one line)
       READ (BFI,1080) KSDOC,KSNH4,KSNO3,KSPO4,KSO,KSSA       !skip 2 lines and read
!       WRITE(*,11051) KSDOC, KSNH4, KSNO3, KSPO4, KSO, KSSA  
       READ (BFI,1080) TRSDOC,TRSNH4,TRSNO3,TRSPO4,TRSO,TRSSA !skip 2 lines and read
!       WRITE(*,11051) TRSDOC,TRSNH4,TRSNO3,TRSPO4,TRSO,TRSSA
       READ (BFI,1080) MTCNO3, SEDNO3, KHSO                   !skip 2 lines and read
!       WRITE(*,11051) MTCNO3, SEDNO3, KHSO
       READ (BFI,1100)                                        !skip 2 lines
                                                              !note here 1100 format skips one line
                                                              !but the read itself skips one line
     ENDIF
!***** Meteorologic data

     READ (MET,1000)

!***** Input FORMAT statements
1000 FORMAT(///)
1001 FORMAT(I6,10(I4,F10.1))
1002 FORMAT(7F10.3)
1003 FORMAT(8F10.3)
1005 FORMAT(//F10.1)
1010 FORMAT(A72)
1020 FORMAT(//(8X,9I8))
1025 FORMAT(8X,9I8)
1027 FORMAT(//(11X,8I8))
1030 FORMAT(//(8X,9F8.0))
1031 FORMAT(//(8X,2F8.0,I8,F8.2))
1032 FORMAT(/)
666 FORMAT(////)
1033 FORMAT(//:(8X,9F8.0))
1035 FORMAT(:///10(10F8.0:/))
!1036 FORMAT(//8X,2F8.1,I8,3F8.1)
!1037 FORMAT(//8X,8F8.1)
!1038 FORMAT(//8X,8F8.1)
1040 FORMAT(//8X,I8,8F8.0)
1050 FORMAT(//13X,A3,I8,5X,A3)
1052 FORMAT(//13X,A3,9I8)
1055 FORMAT(//8X,4(5X,A3),I8)
1060 FORMAT(//8X,2A8,2F8.0)
1065 FORMAT(://(8X,9A8))
1070 FORMAT(//(8X,9(5X,A3)))
1071 FORMAT(//(8X,5(5X,A3),A8,5X,A3,5X,A3,5X,A3))
1072 FORMAT(//(8X,6(5X,A3)))
1080 FORMAT(://8X,9F8.0)
!1082 FORMAT(//8X,2A8//)
1085 FORMAT(8X,9F8.0)
1087 FORMAT(8X,4F8.0)
1090 FORMAT(//(8X,A72))
1100 FORMAT(/)
1110 FORMAT(:////////(8X,5I8))
1120 FORMAT(////(8X,I8))
1130 FORMAT(8X,9I8)
1140 FORMAT(://(2I8))
1150 FORMAT(5X,3F15.0,F18.0,F12.0,I10)
1160 format(://(13X,F13.0))
1170 FORMAT(://(2I10))
1190 FORMAT(I8,9F8.2)
1200 FORMAT(/(10F12.0))
1220 FORMAT(//(:8X,6F8.0))
1230 FORMAT(8X,1F8.0)
11103 FORMAT(10F8.0)
11040 FORMAT(//)
11051 FORMAT(://8X,9F8.4)  !go down 2 lines and write, but do not
                           !continue going down when writing is finished

!************************************************************************
!**                 Initialize Computational Variables                 **
!************************************************************************


!***** Time variables

     NIT      = 0
!     JDAY     = TMSTRT
     JDAYMBL  = TMSTRT
     DLTDP    = 1
     TFLDP    = 1
     KFLDP    = 1
     RSODP    = 1
     SNPDP    = 1
     PLTDP    = 1
     APLDP    = 1
     OPLDP    = 1
     MBLDP    = 1
     DIADP    = 1
     NWQMR    = 0
     DLT      = DLTVAL(DLTDP)
     ELTMSPLT = JDAY*86400.
     ELTMSTFL = JDAY*86400.
     ELTMSKFL = JDAY*86400.
      
     NHMR  = 0
     NXTVD = JDAY
     ELTMS = JDAY*86400.
     MXDLT = DLTMAX(DLTDP)
     FNDLT = DLTFTN(DLTDP)
     NXSNP = SNPD(SNPDP)
     NXPLT = PLTD(PLTDP)
     NXTFL = TFLD(TFLDP)
     NXKFL = KFLD(KFLDP)
     NXOPL = OPLD(OPLDP)
     NXMBL = MBLD(MBLDP)
     NXDIA = DIAD(DIADP)
     NXAPL = APLTD(APLDP)

!***** Output control variables

     SNPD(NSNP+1)  = TMEND+1.
     PLTD(NPLT+1)  = TMEND+1.
     OPLD(NOPL+1)  = TMEND+1.
     MBLD(NMBL+1)  = TMEND+1.
     TFLD(NTFL+1)  = TMEND+1.
     KFLD(NKFL+1)  = TMEND+1.
     RSOD(NRSO+1)  = TMEND+1.
     DLTD(NDLT+1)  = TMEND+1.
     DIAD(NDIA+1)  = TMEND+1.
     APLTD(NAPL+1) = TMEND+1.

     NXCBC = 0.
! KURT GLAESEMANN
     NXMET = 0.
     NXS1  = 0.
     NXS2  = 0.
     NXS3  = 0.
     NXBFI = 0.
     NXATM = 0.
     NXSAV = 0.
! END KURT


     !Accumulation of stuff for diagonostics output
     ACLSZ    = 0.0
     ACLLZ    = 0.0
     ARSZ     = 0.0
     ARLZ     = 0.0
     ABMSZ    = 0.0
     ABMLZ    = 0.0
     AMSZ     = 0.0
     AMLZ     = 0.0
     APRSZLZ  = 0.0
     AGSZ     = 0.0
     AGLZ     = 0.0
     ADOCSZ   = 0.0
     APOCSZ   = 0.0
     ADOCLZ   = 0.0
     APOCLZ   = 0.0
     ANH4SZ   = 0.0
     ADONSZ   = 0.0
     APONSZ   = 0.0
     ANH4LZ   = 0.0
     ADONLZ   = 0.0
     APONLZ   = 0.0
     APO4SZ   = 0.0
     ADOPSZ   = 0.0
     APOPSZ   = 0.0
     APO4LZ   = 0.0
     ADOPLZ   = 0.0
     APOPLZ   = 0.0
     APRSZ    = 0.0
     APRLZ    = 0.0

     AB1SZ    = 0.0
     AB2SZ    = 0.0
     AB3SZ    = 0.0
     AB1LZ    = 0.0
     AB2LZ    = 0.0
     AB3LZ    = 0.0
     ADOSZ    = 0.0
     ADOLZ    = 0.0
     ASASZ    = 0.0
     ASUSZ    = 0.0
     ASALZ    = 0.0
     ASULZ    = 0.0



     IF(SAV_CALC)THEN
     !initialize time average quantities for SAV module
     !Wen Long, this should be moved to SAV_INIT subroutine in wqm_sav.F
        ALEAF  = 0.0
        AROOT  = 0.0    
        ASTEM  = 0.0
        ATUBER = 0.0
        AEP    = 0.0

        APLEAF   = 0.0
        ABMLEAF  = 0.0
        ABMTUBER = 0.0
     
        APEP     = 0.0
        ABMEP    = 0.0
        APREP    = 0.0
        ASLSH    = 0.0

        ANLSAV = 0.0
        APLSAV = 0.0
        ANLEPI = 0.0
        APLEPI = 0.0
        AFNSED = 0.0
        AFPSED = 0.0
        AFHS   = 0.0
     
        AEPATN = 0.0
        AWATATN= 0.0
        AFISH  = 0.0
        AFIEP  = 0.0
        ANPPSAV= 0.0
        ANPPEPI = 0.0
     
        ADOCSAVW    = 0.0
        APOCSAVW    = 0.0
        ADOCEPIW    = 0.0
        APOCEPIW    = 0.0

        ANH4SAVW   = 0.0
        ANO3SAVW   = 0.0
        ADONSAVW   = 0.0
        APONSAVW   = 0.0
        ANH4EPIW    = 0.0
        ANO3EPIW    = 0.0
        ADONEPIW    = 0.0
        APONEPIW    = 0.0

        APO4SAVW   = 0.0
        ADOPSAVW   = 0.0
        APOPSAVW   = 0.0
        APO4EPIW    = 0.0
        ADOPEPIW    = 0.0
        APOPEPIW    = 0.0

        ADOSAVW     = 0.0
        ADOEPIW     = 0.0

        ASEDDOSAV  = 0.0     
        ASEDPOCSAV   = 0.0
        ASEDPONSAV   = 0.0
        ASEDNH4SAV = 0.0          
        ASEDPOPSAV   = 0.0
        ASEDPO4SAV = 0.0
 
     ENDIF
!***** Initial water column mass

     IF (MASS_BALANCE) THEN
       DO JCON=1,NAC
         CMASS(AC(JCON)) = 0.0
         DO K=1,KBM1
            DO I=1,MTLOC
             CMASS(AC(JCON)) = CMASS(AC(JCON))+                       &
               C1(I,K,AC(JCON))*ART1(I)*D(I)*DZ(K)/1000.0
            ENDDO
         ENDDO
       ENDDO
       IWCMN = ANC1*CMASS(4)+ANC2*CMASS(5)+ANC3*CMASS(6)+CMASS(12)    &
               +ANCSZ*CMASS(7)+ANCLZ*CMASS(8)                         &
               +CMASS(13)+CMASS(14)+CMASS(15)+CMASS(16)               &
               +CMASS(17)+CMASS(18)+CMASS(19)
       IWCMP = CMASS(20)+CMASS(21)+CMASS(22)+CMASS(23)+CMASS(24)      &
               +CMASS(25)+APCSZ*CMASS(7)+APCLZ*CMASS(8)               &
               +CMASS(30)+CMASS(31)+CMASS(32)
       IWCMC = CMASS(4)+CMASS(5)+CMASS(6)+CMASS(7)+CMASS(8)           &
               +CMASS(9)+CMASS(10)+CMASS(11)+CMASS(12)
       IWCMS = ASC1*CMASS(4)*ASC2*CMASS(5)+ASC3*CMASS(6)+CMASS(28)    &
               +CMASS(29)
     ENDIF

!***** Change in concentrations

     DTC = 0.0
    
!
!    *** DO days
!
     DOVDAYS = 0.0

!***** Maximum number of Sources One through Three

     DO JCON=1,NCP
       IF (S1LN(JCON) > S1LNMAX) THEN
         S1LNMAX = S1LN(JCON)
         JCS1MAX = JCON
       ENDIF
       IF (S2LN(JCON) > S2LNMAX) THEN
         S2LNMAX = S2LN(JCON)
         JCS2MAX = JCON
       ENDIF
       IF (S3LN(JCON) > S3LNMAX) THEN
         S3LNMAX = S3LN(JCON)
         JCS3MAX = JCON
       ENDIF
     ENDDO

!***** Horizontal diffusion

     VISCOFH = 0.0
     IF (XY_DIFFUSION) THEN
       DO K=1,KBM1
         DO I=1,MLOC
           VISCOFH(I,K) = XYDF
         ENDDO  
       ENDDO
     ENDIF


!************************************************************************
!**                   Parameters for Sediment Submodels                **
!************************************************************************

    IF(SEDIMENT_CALC)THEN
        CALL SED_INIT()          !  Setup sediment diagenesis module (including SED_ALLOC and SED_INIT_ICI)
        CALL SED_READ            !
        CALL SED_INIT2()         !  initialize the values from input read by INPUTS() above

     ENDIF	 
	 
!Set up netcdf files --------------------------------------------------

! B Clark added for netcdf output Sep 2018
	!  write(*,*)'About to call set NCDIO'
       CALL SET_NCD_IO
	!  write(*,*)'Called set NCDIO'

     RETURN

END SUBROUTINE INPUTS

     SUBROUTINE ERROR_CHECK(HEADER, FLAG)
    !SUBROUTINE to print an error message and stop the code
    !           or print a warning message only
     IMPLICIT NONE
     CHARACTER(LEN=120) HEADER
     INTEGER FLAG
     IF (FLAG .EQ. -1) THEN
       WRITE(*,"(a120)")HEADER !print error message in report.txt
       STOP
     ELSE
       WRITE(*,"('***** WARNING *****')") !print warning message on the screen
       WRITE(*,"(a120)")HEADER
     ENDIF
     RETURN
     END SUBROUTINE ERROR_CHECK
