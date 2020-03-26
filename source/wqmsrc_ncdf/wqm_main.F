!************************************************************************
!** Modified by JLM 2/22/01 to allow multiple boundary files           **
!**                                                                    **
!** added temporal diffusion variation to agree with hw's ch3d         **
!**  modified so that flux area is not dependent upon VFN              **
!**  bwb 1-12-93                                                       **
!**  corrected ratio calculation for non step boundary conditions      **
!**  bwb 1-12-93                                                       **
!**  modified error calculations so that mass balance can be run       **
!**  when sediment model is OFF  bwb 1-12-93                           **
!**  modified hardwired diffusion for Indian River to correspond to    **
!**  Harry Wang's CH3D diffusion scheme.  bwb 1-19-93                  **
!**  inserted Carl's new autostepping stability criteria bwb 1-20-93   **
!**  references to SALINITY_RUN have been removed   MNOEL 1-15-93      **
!**  references to HydroQual and CH3d have been removed now            **
!**  referenced as ASCII and binary    MNOEL 1-25-93                   **
!**  references to point and non-point sources have been replaced with **
!**  Source One and Source Two   MNOEL  1-25-93                        **
!**  Saving of final output for binary initial conditions added        **
!**  MNOEL 2-5-93                                                      **
!**  Removal of Ocean boundary conditions faces from DATA statements   **
!**  read as ASCII from wqm_map file   MNOEL  2-10-93                  **
!**  additional parameters made spacially variable-- note preprocessor **
!**  needed to produce ASCII input file (parameters specified by       **C
!**  Cerco)   MNOEL 2-20-93                                            **
!**  Modified Kinetic flux output to be aggregated temporally and to   **
!**  write for every box in a binary output file                       **
!**  Benthos inputs from HydroQual's Oct 97 transmission               **
!**                                                                    **
!**                                                                    **
!**                                                                    **
!**                                                                    **
!**                           CE-QUAL-IC                               **
!**                                                                    **
!**                    An Integrated Compartment                       **
!**                       Water Quality Model                          **
!**                                                                    **
!**                           Version 1.1                              **
!**                          July 28, 1992                             **
!**                                                                    **
!**                           Developed by                             **
!**                                                                    **
!**             Carl F. Cerco      : Water quality scheme              **
!**             Raymond S. Chapman : Numerical solution scheme         **
!**             Thomas M. Cole     : Computer algorithms & coding      **
!**             Hydroqual          : Sediment compartment              **
!**                                                                    **
!**                    Water Quality Modeling Group                    **
!**                    U.S. Army Corps of Engineers                    **
!**                    Waterways Experiment Station                    **
!**                    Vicksburg, Mississippi 39180                    **
!**                                                                    **
!************************************************************************

   PROGRAM  PARWQM   
   
    USE MOD_PREC, ONLY: SP
    USE MOD_SIZES, ONLY: NCP, NSSFP, NGL
    USE MOD_LIMS, ONLY: MYID, NPROCS, IINT, NTLOC,NLOC, MLOC, KB
    
    USE MOD_TGE, ONLY: NV, TGE_DEALLOC
    
    USE MOD_HYDROVARS, ONLY:  &
           !GRAV    &        !
          PI    &        !
        !,PI2    &        !
        !,ZERO    &        !
        !,ONE_THIRD    &    !
!         ,NVG    &        !
!         ,XG    &        !GLOBAL X-COORD AT NODE 
!         ,YG    &        !GLOBAL X-COORD AT NODE 
!         ,HG    &        !GLOBAL DEPTH AT NODE 
!         ,XCG    &        !GLOBAL X-COORD AT FACE CENTER 
!         ,YCG    &        !GLOBAL X-COORD AT FACE CENTER 
        !,VXMIN    &        !
        !,VYMIN    &        !
        !,VXMAX    &        !
        !,VYMAX    &        !
!         ,XC    &        !X-COORD AT FACE CENTER 
!         ,YC    &        !Y-COORD AT FACE CENTER
!         ,VX    &        !X-COORD AT GRID POINT
!         ,VY    &        !Y-COORD AT GRID POINT
!         ,ART    &        !AREA OF ELEMENT
         ,ART1    &        !AREA OF NODE-BASE CONTROl VOLUME
!         ,ART2    &        !AREA OF ELEMENTS AROUND NODE
!         ,NV    &        !NODE NUMBERING FOR ELEMENTS
!         ,NBE    &        !INDICES OF ELMNT NEIGHBORS
!         ,NTVE    &        !
!         ,NTSN    &        !
!         ,ISONB    &        !NODE MARKER = 0,1,2 
!         ,ISBC    &        !
!         ,ISBCE    &        !
!         ,IEC    &        !
!         ,IENODE &        !
!         ,NBSN    &        !
!         ,NIEC    &        !
!         ,NTRG    &        !
!         ,NBVE    &        !
!         ,NBVT    &        !
!         ,LISBCE_1    &    !LIST OF ELEMENTS WITH ISBCE=1
!         ,LISBCE_2    &    !LIST OF ELEMENTS WITH ISBCE=2
!         ,LISBCE_3    &    !LIST OF ELEMENTS WITH ISBCE=3
!         ,DLTXC    &        !
!         ,DLTYC    &        !
!         ,DLTXYC    &    !
!         ,DLTXE    &        !
!         ,DLTYE    &        !
!         ,DLTXYE    &    !
!         ,SITAC    &        !
!         ,SITAE    &        !
!         ,XIJC    &        !
!         ,YIJC    &        !
!         ,XIJE    &        !
!         ,YIJE    &        !
!         ,EPOR    &        !ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
        !,IBCGEO    &    !LOCAL GEOSTROPHIC FRICTION CORRECTION NODES

         ,Z    &            !SIGMA COORDINATE VALUE 
         ,ZZ    &        !INTRA LEVEL SIGMA VALUE
         ,DZ    &        !DELTA-SIGMA VALUE
         ,DZZ    &        !DELTA OF INTRA LEVEL SIGMA 
!         ,H1    &        !BATHYMETRIC DEPTH 
		 !,siglev &
		 !,siglay &
		 
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
         ,UARD_OBCNT    &!tykim previous time step
         ,XFLUX_OBCT    &!tykim previous time step
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
         ,THOUR    &        !
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
         ,t_his_start    &        !
         ,t_his_end    &            !
         ,t_his_dlt    &            !starting time, ending time, and interval of history outputs (days)
         ,Nstation    &            !
         ,NstationNum_GL    &    !maximum number of station is NstationMax!
         ,t_stn_start    &        !
         ,t_stn_end    &            !
         ,t_stn_dlt    &            !starting time, ending time, and interval of station outputs (days)
         ,STNFN    &                !file name for station output
         ,HISFN    &                !file name for history output
         ,HISFN_PREFIX    &        !prefix of history output file
         ,HISFN_EXT    &            !extention name of history output file
         ,HISFN_FINAL    &        ! 
         ,HISFN_SPLIT_BYLEVEL & !True or False for splitting history output in files level by level (default is .FALSE.)
        !,hydro_netcdf    &        !
        !,wqm_history    &        !
        !,wqm_stations    &        !
        !,IFNC            &        !file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise 
        !,NTRECNC        &        !time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file. 
        !,NTHYDRO        &        !overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read
         ,HYDRO_GEOM_DEALLOC &    !
         ,HYDRO_DEALLOC
		 
    USE MOD_HYDRO, ONLY:			&!
#if defined (MULTIPROCESSOR)
		 BROADCAST_HYDRO,	&		!BROADCAST_HYDRO()
		 BROADCAST_HYDRO_REDUCED, &	!BROADCAST_HYDRO_REDUCED()
		 SCATTER, &					!SCATTER()
#endif		 
        HYDRO

    USE MOD_BCMAP, ONLY :    &
            !I_OBC_N,        &!
            !I_OBC_GL,        &!
            BCMAP_DEALLOC

       !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY :         &
            SERIAL          &           !!TRUE IF SINGLE PROCESSOR
            ,MSR            &           !!TRUE IF MASTER PROCESSOR (MYID==1)
            ,PAR            &           !!TRUE IF MULTIPROCESSOR RUN
            ,CASENAME      !    &           !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
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
            
       USE MOD_FILEINFO, ONLY : &
                DIA                &
                !,CBC             &
                !,S1            &
                !,S2            &
                !,S3            &                  
                !,BFI            &
                !,BAI           &
                !,MET            &
                ,BFO            & 
                !,KEI            &
                !,ATM            &
                !,STL            & 
                !,AGR            & 
                !,SVI            & 
                !,SVO            & 
                ,KFL            & 
                !,ZOO            & 
                !,ZFO            & 
                !,ALO              &
                !,CON            &
                ,RSO            &
                ,SNP            &
                ,PLT            &
                ,APL             &
                ,TFL             &
                ,OPL            &
                !,SFI            &
                ,SFO            &
                !,MAP             &
                !,ICI             &
                ,ICO            &
                !,MRL            &
                ,MBL            &
                !,RSI            &
                !,UNIT_LINKAGE    &
                ,UNIT_STN        &
                ,UNIT_HIS        &            
                ,CNAME            &
                !,CONFN
                ,INIT_FILE_INFO

    
    USE MOD_SF, ONLY:             &     !
                SFEEDER,         &    !
                !HYPOXFX_SF,        &    !
                !FRSASF,        &    !
                NSPECIES,        &    !
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
                !JNH4SF,        &!(:),            
                !JPO4SF,        &!(:),            
                !SODSF,            &!(:),            
                !JSASF,            &!(:),            
                !JSUSF,            &!(:),            
                !SF_SSI,        &!(:),            
                !SF_SU,            &!(:),        
                !SF_SA,            &!(:),            
                !SF_PIP,        &!(:),            
                SFEED,            &!(:,:),
                SFEED_GL,        &!(:,:),
                !SEDTYPE_GL,        &!(:,:)
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
                SUSPFEED,        &!        
                !SF_ALLOC        &!
                SF_DEALLOC!,    &!
                !SF_READ,        &!    
                !SF_INIT
                
    USE MOD_WQM, ONLY:             &!
                JDAY,           &!    
                DLT,            &!
                AHMDLT,         &!                
                DLTAV,          &!
                ELTMS,          &!
                ELTMS1,         &!
                ELTMJD,         &!
                XFLUX, &!

                ! XFLUX_OUT,  &
                WQ_DT,          &!    

                SEDIMENT_CALC,  &!

                T,              &!
                C2,                &!
                SALT,            &!
                SSI,            &!
                B1,                &!
                B2,                &!
                B3,                &!
              !  LDOC,            &!
             !   RDOC,            &!
                LPOC,            &!
                RPOC,            &!
                NH4,            &!
                NO3,            &!
                UREA,            &!
            !    LDON,            &!
             !   RDON,            &!
                LPON,            &!
                RPON,            &!
                PO4,            &!
             !   LDOP,            &!
             !   RDOP,            &!
                LPOP,            &!
                RPOP,            &!
                PIP,            &!
                COD,            &!
                DOXG,            &!
                SIUPB,            &!
                SIAT,            &!
                PIB1,            &!
                PIB2,            &!
                PIB3,            &!
                Q1,             &!
                Q2,             &!
                Q3,             &!
                BENTHIC_FLUXES, &!
                DSSR,           &!
                DTC,            &!
                FD,             &!
                FLOW,           &!
                I0,             &!
                IT,             &!
                LIGHT_EXTINCTION, &!
                RELTMS,         &!
                SAV_CALC,       &!
                SAVENT,         &!
                SEDENT,         &!
                TODS,           &!
                TTSS,           &!
                Z_DIFFUSION,    &!
                AC,             &!
                ANC1,           &!
                ANC2,           &!
                ANC3,           &!
                ASC1,           &!
                ASC2,           &!
                ASC3,           &!                
                C1,             &!
                C2F,            &!
                CMASS,          &!
                DTM,            &!

                FLUXS,          &!
                FNDLT,          &!
                KBM1,           &!
                MASS_BALANCE,           &!
                NAC,            &!
                NIT,            &!
                SFENT,          &!
                SOLIDS_CALC,            &!
                TMEND,          &!
                TMSTRT,         &!
                
                WCMC,           &!
                WCMN,           &!
                WCMP,           &!
                WCMS,           &!
                
                XYDF,           &!
                XY_DIFFUSION,   &!
                

                ATMFLXN,                &!
                ATMFLXNB,               &!
                ATMFLXP,                &!
                ATMFLXPB,       &!
                ATMLDON,        &!
                ATMLDOP,        &!
                ATMNH4,         &!
                ATMNO3,         &!
                ATMPO4,         &!
                ATMRDON,        &!
                ATMRDOP,        &!

                BENFLXDN,       &!
                BENFLXPN,       &!
                BENFLXPNB,      &!

                DLWCMC,         &!
                DLWCMN,         &!
                DLWCMP,         &!

                IWCMC,          &!
                IWCMN,          &!
                IWCMP,          &!

                PRECIP,         &!
                SFATOT,         &!
                
                ANDC,           &!
                BENFLXDNB,              &!
                BENFLXDP,               &!
                BENFLXDPB,              &!
                BENFLXPC,               &!
                BENFLXPCB,              &!
                BENFLXPP,               &!
                BENFLXPPB,              &!
                
                BENNH4,         &!
                BENNO3,         &!
                BENSA,          &!                
                BENPO4,         &!
                
                BNDN,           &!
                DENIT,          &!
                DLALGC,         &!
                DLWCKMC,                &!
                DLWCKMCB,               &!
                DLWCKMN,                &!
                DLWCKMNB,               &!
                DTMP1,          &!

                FLXTLDON,               &!
                FLXTLPON,               &!
                FLXTNH4,                &!
                FLXTNO3,                &!
                FLXTRDON,               &!
                INFLOW,         &!
                JCB,            &!
            !    MNLLDOC,                &!  ! B Clark moved into 3 pools in mod_wc_dom Sep 2015
           !     MNLRDOC,                &!  ! B Clark moved into 3 pools in mod_wc_dom Sep 2015
                PCFWS,          &!
                PPFWS,          &!
                RESPC,          &!

                BNDC,           &!
                BNDFLXC,                &!
                BNDFLXN,                &!
                BNDFLXP,                &!
                BNDP,           &!
                DLWCSSMN,               &!
                DLWCSSMP,               &!


                FLXT1,          &!
                FLXT2,          &!
                FLXT3,          &!
                FLXTLDOC,               &!
                FLXTLDOP,               &!
                FLXTLPOC,               &!
                FLXTLPOP,               &!
                FLXTLZ,         &!
                FLXTPIB1,               &!
                FLXTPIB2,               &!
                FLXTPIB3,               &!
                FLXTPO4,                &!
                FLXTRDOC,               &!
                FLXTRDOP,               &!
                FLXTRPOC,               &!
                FLXTRPON,               &!
                FLXTRPOP,               &!
                FLXTSZ,         &!
                S1FLXN,         &!
                S1FLXP,         &!
                S2FLXN,         &!
                S3FLXN,         &!

                BENDEN,         &!
                BURIALFLXC,             &!
                BURIALFLXCB,            &!
                BURIALFLXN,             &!
                BURIALFLXNB,            &!
                BURIALFLXP,             &!
                BURIALFLXPB,            &!
                DFTMP,          &!
                DLSEDKC,                &!
                DLSEDKCB,               &!
                DLSEDKN,                &!
                DLSEDKNB,               &!
                DLSEDMC,                &!
                DLSEDMN,                &!
                DLSEDMP,                &!
                DLSEDSSMC,              &!
                DLSEDSSMN,              &!
                DLSEDSSMP,              &!
                DLWCSSMC,               &!
                DTMP2,          &!
                ERRSEDN,                &!
                ERRWCC,         &!
                ERRWCN,         &!
                ERRWCP,         &!

                S1FLXC,         &!
                S2FLXC,         &!
                S2FLXP,         &!
                S3FLXC,         &!
                S3FLXP,         &!

                C2_GL,          &!
                CCHL1,          &!
                CCHL1_GL,               &!
                D_GL,           &!
                EL_GL,          &!
                ERRSEDC,                &!
                ERRSEDP,                &!

                H_GL,           &!
                MGL,            &!
                MTLOC,          &!
                P1,             &!
                S_GL,           &!
                T_GL,           &!
                
                B1SZ_GL,                &!
                BM1,            &!
                BM1_GL,         &!
                BM2,            &!
                BM2_GL,         &!

                FI1,            &!
                FI1_GL,         &!
                FI2,            &!
                FI2_GL,         &!
                NL1,            &!
                NL1_GL,         &!
                NL2,            &!
                NL2_GL,         &!
                P1_GL,          &!
                P2,             &!
                P2_GL,          &!
                PL1,            &!
                PL1_GL,         &!
                PL2,            &!
                PL2_GL,         &!

                B1LZ_GL,        &!
                B2LZ_GL,        &!
                B2SZ_GL,        &!
                DCOD,           &!
                DCOD_GL,        &!
                DDOC,           &!
                DDOC_GL,        &!

                IAVG_GL,        &!
				KDLAMB_int_GL,   & ! B Clark May 2018,
				SSI_GL,  & ! B Clark June 2018
                IK1,            &!
                IK1_GL,         &!
                IK2,            &!
                IK2_GL,         &!
                NITRIF,         &!
                NITRIF_GL,              &!
                PR1,            &!
                PR1_GL,         &!
                PR2,            &!
                PR2_GL,         &!

                BFLUX,          &!
                BFLUX_GL,       &!
                COD_GL,         &!
                COURMX,         &!
                DIAGNOSTICS,    &!
                DOLZ_GL,        &!
                DOPR,           &!
                DOPR_GL,        &!
                DORALG,         &!
                DORALG_GL,      &!
                DOSZ_GL,        &!
                FTCOD,          &!
                FTCOD_GL,       &!
                REAERDO,        &!
                REAERDO_GL,     &!
                RESP,           &!
                RESP_GL,        &!
                PNFWS,          &!                
                CTEMP,             &!
                DIFFMX,            &!

                B,              &!
                BB,             &!                
                IJDAY,          &!
                LENGTH,         &!
                SAV_LOADS,      &!
                WQM_DEALLOC, &  !
				UL_GL,VL_GL, &
				SL_GL,TL_GL, &
				WTSL_GL,	&
				KHL_GL, &
				total_netPP, &
				ALGDON,    &!Algae DON production
				ALGDON_OUT,    &! OUTobal Algae DON production

				
				ALGPON,    &!Algae PON production from predation
				ALGPON_OUT,    &! OUTobal Algae PON production
				
				ALGNH4,  & ! Algal delta NH4
				ALGNH4_OUT ,   &!
				
				ALGNO3 ,   &!  Algal delta NO3
				ALGNO3_OUT,  &  
				
				DENNO3,   & ! Denitrification loss of NO3
				DENNO3_OUT,   &!
				
				NT,   & ! Nitrification loss rate of NH4 (production rate of NO3 from nitrification)
				NT_OUT,   &
				
				DENITN,   & ! Denitrificaion loss rate of DON (remin of DON from denitrification)
				DENITN_OUT,  & !
				
			    DENIT,          & ! loss of DOC due to denitrification
				DENIT_OUT,   & !
				
				LPON,            &!
				
                RPON,            &!
		
				
				HDRLPON,  & ! Hydrolysis of LPON to DON
				HDRLPON_OUT,  & !
				
				HDRRPON,  & ! Hydrolysis of RPON to DON
				HDRRPON_OUT,  & ! 	
				
				ALGDOC,    & ! Algae DOC production
				ALGDOC_OUT,    & ! Global Algae DOC production,
			
				
				ALGPOC,    & !Algae POC production from predation
				ALGPOC_OUT,    & ! OUTobal Algae POC production
				
				HDRLPOC,  & ! Hydrolysis of LPOC to DOC
				HDRLPOC_OUT,  &! 
				
				HDRRPOC,  & ! Hydrolysis of RPOC to DOC
				HDRRPOC_OUT,  &
				DTLPOC_out, DTRPOC_out,  &
				DTLPON_out, DTRPON_out, &
				DTNH4_out, DTNO3_out

                
    USE MOD_OWQ, ONLY:                  & 
                        !INITKE,        &!
                        INTKE,            &!
                        NXKEI,            &!
                        !LAT,             &!
                        !LONG,            &!
                        !GSTAR440,        &!
                        !PHISTAR676,    &!
                        !BNSTARNTU,        &!
                        !BNSTARCHL,     &!
                        !PSTARINT,        &!
                        !PSTAR440,        &!
                        !PSTARCOLOR,    &!
                        !SCDOM,            &!
                        !EZEROINT,        &!
                        !TCHL2CHL,        &!
                        !DOFFSET,        &!
                        !DEPTHKE,        &!
                        !NWAVEL,        &!
                        !COLOR,         &!
                        !TURB,            &!
                        !OPZ,             &!
                        !G1,            &!
                        !G2,            &!
                        !ALAMB,         &!
                        !GLAMB,            &!
                        !PLAMB,            &!
                        !PHILAMB,       &!
                        !BLAMB,            &!
                        !WAVEL,            &!
                        !EZERO,            &!
                        !GMUZERO,        &!
                        IAVG,            &!
                        KESS,            &!
                        !IATBOT,        &!            
                        !OWQ_ALLOC,        &!
                        OWQ_DEALLOC,    &!
                        !OWQ_READ,        &!
                        LGHT_ATTN
                        
    USE MOD_WQMINIT, ONLY :            &        !
            RSOD,                &        !    
            RSODP,                 &        !
            RSOFN,                 &        !
            DLTD,                &        !
            DLTDP,               &        !
            !SNPDP,                &        !
            !TFLDP,                &        !
            !PLTDP,                &        !
            !APLDP,               &        !
            !OPLDP,                &        !
            DIAD,                 &        !            
            DIADP,                &        !
            DIAF,                 &        !
            
            !COURB,                &        !
            COURFS,            &        !
            COURBS,            &        !
            F,                    &        !            
            SB,                    &        !
            !S1LNMAX,            &        !
            !S2LNMAX,            &        !
            !S3LNMAX,            &        !
            DIFFFS,                &        !
            DIFFBS,              &        !
                
            !NXSNP,             &        !
            !NXPLT,                &        !
            !NXAPL,                &        !
            !NXTFL,                &        !
            !NXKFL,                &        !
            NXTVD,               &        !
            !NXOPL,                &        !
            !NXMBL,                &        !
            NXDIA,                &        !

            !LEAFI,                &        !
            !STEMI,                &        !
            !ROOTI,             &        !
            !TUBERI,            &        !
            
            !MAXDLT,            &        !
            MXDLT,                &        !
            !INTKE,                &        !
            !INITKE,            &        !
            TM1,                &        !
            TM2,                &        !
            !TM3,                &        !
            !MINSTEP,            &        !
            !NXKEI,                &        !
                        
            !SNPC,                &        !
            !RSIC,                 &        !
            !RSOC,                &        !
            !BCC,                &        !
            !S1C,               &        !
            !S2C,                &        !
            !MDC,                &!*3,    !    
            !PLTC,                &        !
            !FLC,               &        !
            !MBLC,                &        !
            !BFC,                &        !
            !VBC,                &        !
            !QPLTC,                &        !
            !XYDFC,                &        !
            !XYDFU,              &        !
            !ZDFC,                &        !
            !ICOC,                &        !
            !ATMC,                &        !
            !SAVLC,                &        !
            !SEDC,              &        !
            !AUTOC,                &        !
            !SPLTC,                &        !
            !TFLC,                &        !
            !DIAC,                &        !
            !STLC,              &        !
            !APLTC,                &        !
            !KFLC,                &        !
            !OPLC,                &        !
            !BFOC,                &        !
            !DFOC,                &        !
            !S3C,                  &        !
            !SAVMC,                &        !
            !SAVPLTC,            &        !
            !DFLC,                &        !
            !SFLC,                &        !
            !SFLOXC,            &        !
            !DFLOXC,            &        !
            EXT1,                 &!*1,   !    
            EXT2,                &!*2,    !    
            EXT3,                &!*3,    !    
            EXT4,                &!*4,    !    
            !EXT_R,                &!*1,    !
            !EXT_R2,            &!*2    !
             
            !SLC,                &        !
            !HYDC,                &        !
            !BNDTC,                &        !
            !CONSC,                &        !
            !ICIC,                &        !
                                
            !TIMEA,                &        !
            !TIMEB,                &        !
            !SPVARM,            &        !
            !PRINTM,            &        !
                                
            TITLE,                &!(6),    !
            !OLDTITLE            &!(6),    !
            FILENAME,            &        !
                                    
            !FILENAME2,            &        !
                
            !MAPFN,                &        !
            !GEOFN,                &        !
            !ICIFN,                &        !
            !RSIFN,              &        !
            !AGRFN,                &        !
            !STLFN,                &        !
            !MRLFN,                &        !
            !KFLFN,                &        !
            !ICOFN,                &        !
            
            !SNPFN,                &        !
            !PLTFN,                &        !
            !APLFN,             &        !
            !DIAFN,             &        !
            !TFLFN,              &        !

            !OPLFN,             &        !
            !MBLFN,             &        !
            !SFIFN,             &        !
            !SFOFN,                &        !

            RESTART_OUT,        &        !
            !SNAPSHOTS,            &        !
            END_RUN,            &        !
            !MODIFY_ICONC,       &          !
            VOLUME_BALANCE,     &        !
            !QUICKEST,            &        !
            !UPWIND,            &        !
            ICOND_OUT,          &        !
            !UNI_ICON_IN,       &        !
            !BIN_ICON_IN,        &        !
            !RES_ICON_IN,        &        !
            AUTO_STEPPING,      &        !
            !STOP_RUN,           &        !
            !PLOTS,                &        !
            !OXYGEN_PLOTS,       &        !
            
            NEW_VOLUMES,          &        !
            !RESTART_IN,        &        !
                                
            !TEMPERATURE_CALC,     &        !
            ALGAE_CALC,        &        !
            ZOO_CALC,           &        !
            CARBON_CALC,       &        !
            NITROGEN_CALC,     &        !
            PHOSPHORUS_CALC,   &        !
            COD_CALC,            &        !
            OXYGEN_CALC,        &        !
            SILICA_CALC,        &        !
                                
            LEFT_FLOWB,         &!(NHQP),        !        
            !RIGHT_FLOWB,        &!(NHQP),        !    
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

            IFLOWP,            &    !(NBCP),        !
            !IBT,                &    !(NHQP),          !
            !COURBSV,            &    !(NHQP),        !
            !COURVSV,            &    !(NHQP)            !

            C1MIN,              &                !
            C1MAX,                &                !

            DOVDAYS,            &    !(:,:,:)    !
            !OINT,                &    !(NOIP)        !
            ACC,                &    !(NCP)        !
             
            !CIC,                &    !(NCP),        !
            !COUT,                &    !(NCP),        !

            !SFEEDI,            &    !(10)        !
            
            !SNPD,                &!
            !SNPF,                 &!
            
            DLTVAL,                &!
            DLTMAX,               &!
            DLTFTN,             &!
            !PLTD,                &!
            !PLTF,                 &!
            !APLTD,                &!
            !APLF,                &!
            !TFLD,                 &!
            !TFLF,                &!
            !KFLD,                &! 
            !KFLF,                 &!
            !OPLD,                &!
            !OPLF,                &!
            

            !MBLD,                 &!
            !MBLF,                &!

            !DOCLIT ,             &!
            !LPOCLIT,              &!
            !POCLIT,              &!
            !PBSLIT,              &!
            !DSILLIT,            &!
            !WQMINIT_ALLOC,        &!Subroutine to initialize MOD_WQM data
             WQMINIT_DEALLOC
             
   USE MOD_ALGAL, ONLY :         &
            !ALG_READ,            &
            ALGAE    
   
   USE MOD_ZOOP, ONLY :            &
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
                     CLSZ,        & 
                     CLLZ,        &
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
                      B1SZ,        &!
                      B2SZ,        &!
                     !B3SZ,        &!
                     B1LZ,        &!
                     B2LZ,        &!
                     !B3LZ,        &!
                     DOSZ,        &!
                     DOLZ,        &!
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
                     SZ,        &!
                     LZ,        &!
                     ZOOP_DEALLOC
                     
    USE MOD_KIN, ONLY :                &
                        TEMPER,        &
                        SOLIDS,        &
                        ZOOPL,        &
                        PARTICULATE_CARBON,        &
                        NITROG,        &
                        PHOSPH,        &
                        CODMND,        &
                        OXYGEN,        &
                        SILICA,        &
                        BEN_FLUX

   USE MOD_SED, ONLY :         &!
       !    BENTHIC_OUTPUT,     &!
        STEADY_STATE_SED_G3,&!
        STEADY_STATE_SED_IC,&!
        STEADY_STATE_SED_TS,&!
        DLTS,                &!
        NDTSED,              &!LB
    !    SSTATEG3,           &!
    !    HSEDALL,               &!
    !    DIFFT,                &!
    !    SALTSW,             &!
    !    SALTND,                &!    
    !    FRPALG1,            &! 
    !      FRPALG2,            &! 
    !      FRPALG3,            &! 
    !    FRNALG1,            &! 
    !   FRNALG2,               &! 
    !   FRNALG3,               &!
    !    FRCALG1,            &!
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
        IERR,                &!
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
        ISEDMN,               &!
        ISEDMP,               &!
        ISEDMC,                &!
        SEDMN,              &!total sediment organic nitrogen  (kgN)
        SEDMP,              &!total sediment organic phosphorus (kgP)
        SEDMC,                &!total sediment organic carbon  (kgC)
        !NH4AVL,               &!    
        !NO3AVL,               &!
        !KETOTL,            &!
        !WSSINETMMD,         &!
        !WLPOMNETMMD,        &!
        !WRPOMNETMMD,        &!
        !WB1NETMMD,            &!
        !WB2NETMMD,            &!   
        !WB3NETMMD,            &!   
        !WPSINETMMD,        &!   
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
    !    BFORMAXS,           &!      
    !    ISWBENS             &!    
        SODTM1S,            &
        JNH4TM1S,     & 
        JNO3TM1S,     & 
		
		JNH4_out, &!Added By B Clark Feb 2019 for Budget calculation
		JNO3_out, &
		
        JPO4TM1S,     & 
        JCH4TM1S,     & 
        JCH4GTM1S,    & 
        JHSTM1S,    &
        JSITM1S,    &
        NH41TM1S,              &!  
        NH42TM1S,              &!  
        NH4T2TM1S,             &!  
        NO31TM1S,          &! 
        NO32TM1S,          &! 
        NO3T2TM1S,             &!  
        HS1TM1S,           &! 
        HS2TM1S,           &! 
        HST2TM1S,              &! 
        SI1TM1S,           &! 
        SI2TM1S,           &! 
        SIT2TM1S,              &!  
        PO41TM1S,          &!  
        PO42TM1S,          &!  
        PO4T2TM1S,             &!  
        CH41TM1S,              &!  
        CH42TM1S,              &!  
        CH4T2TM1S,             &!  
        SO4T2TM1S,            &!     
        BURIALC,            &!     
        BURIALN,            &!     
        BURIALP,               &!
        DIAGENC,              &!   
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
    !   BBM,                  &!  
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
    !    PSISEDTM1             &!
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
        BENDEN_GL,            &!
        JPOC_GL,            &
        JPON_GL,            &
        JPOP_GL,            &
        JPOS_GL,            &
        CH40_GL,            &
        JNH4_GL,            &
        JNO3_GL,            &
        JCH4_GL,            &
        JCH4G_GL,            &
        JHS_GL,                &
        JPO4_GL,            &
        JSI_GL,                &
        NH41_GL,            &
        NH42_GL,            &
        NO31_GL,            &
        NO32_GL,            &
        PO41_GL,            &
        PO42_GL,            &
        SI1_GL,                &
        SI2_GL,                &
        CH42_GL,            &
        HS1_GL,                &
        HS2_GL,                &
        HSED1_GL,            &
        HSED1,                &
        
       JPOC,         &!
      JPON,         &!
      JPOP,         &!  
	  
	  !B Clark Feb 2019
	  JPOC_out, JPON_out, &
	  
      JPOS,         &!
      JPOCaccum,         &!
      JPONaccum,         &!
      JPOPaccum,         &!                
      JPOSaccum,         &!
!        CCH4,           &!
!       CHS,            &!
!       CSI,            &!
!       CSO4,           &!                
        CPOC,               &!
        CPON,               &!
        CPOP,               &!
        CPOS,               &!
        CPO4,                 &!
        CNH4,                 &!
        CNO3,                 &!                
        CHS,                 &!
        CSI,                 &!                
        DIAGN,                 &!                    
        !SED_INIT,            &!
        SED_INIT2,            &!
        SED_INIT_ICI,         &!
    !    SED_ALLOC,            &!
        SED_DEALLOC,        &!
    !    SED_READ,            &!
        POM_ACCUMUL,         &!LB
        SED_CALC,            &!
        SED_DIAGENESIS_G3,            &!
    !    SEDTSFNL,            &!
    !    SEDSSFNL,            &!
       !    SEDF,                &!
    !    ZBRENT        
        AGGREGATE_POM_FLUX,   &!LB  Option to accumulate all labilities from water column and then distribute them to labilities in sediments
        FRACL1,               &!LB
        FRACL2,               &!LB
        FRACL3                 !LB
        
        
        
   USE MOD_SED_DF_EXCHANGE_VARS, ONLY:     SED_DF_EXCHANGE_ALLOC,    &!
                                        SED_DF_EXCHANGE_DEALLOC
   
   USE MOD_SED_SF_EXCHANGE_VARS, ONLY:     SED_SF_EXCHANGE_ALLOC,    &!
                                        SED_SF_EXCHANGE_DEALLOC
   
   USE MOD_SED_SAV_EXCHANGE_VARS, ONLY: SED_SAV_EXCHANGE_ALLOC,    &!
                                        SED_SAV_EXCHANGE_DEALLOC
   
   USE MOD_BA, ONLY : BALGAE_CALC, BBM,BBM_GL, BA_DEALLOC
   
   USE MOD_OBCS, ONLY:  &
        !NEXT_OBC ,       & !!INTERIOR NEIGHBOR OF OPEN BOUNDARY NODE
        !NEXT_OBC2,       & !!INTERIOR NEIGHBOR OF NEXT_OBC
        !!TYPE_OBC ,       & !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
        !TYPE_OBC_GL,     & !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
        !!IBCN       ,     & !!NUMBER OF EACH TYPE OF OBN IN LOCAL  DOM
        !!IBCN_GL    ,     & !!NUMBER OF EACH TYPE OF OBN IN GLOBAL DOM
        !OBC_LST    ,     & !!MAPPING OF OPEN BOUNDARY ARRAYS TO EACH TYPE
        !NADJN_OBC  ,     & !!NUMBER OF ADJACENT OPEN BOUNDARY NODES TO OBN
        !ADJN_OBC   ,     & !!ADJACENT OBN's of OBN
        !NADJC_OBC  ,     & !!NUMBER OF ADJACENT OPEN BOUNDARY CELLS TO OBN
        !ADJC_OBC   ,     & !!ADJACENT OPEN BOUNDARY CELLS
        !!NFLUXF_OBC ,     & !!NUMBER OF FLUX SEGMENTS TO OBN
        !!FLUXF_OBC  ,     & !!FLUX FRACTION ON EACH SIDE OF OBN
        !NXOBC      ,     & !!INWARD POINTING X-NORMAL OF OBN
        !NYOBC      ,     & !!INWARD POINTING Y-NORMAL OF OBN
        !DLTN_OBC   ,     & !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY
        !NUT_OBC    ,     & !!NUTRIENT AT OPEN BOUNDARY
        !NUT_OBC_GL ,     & !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY

        !SEPARATE_OBC,     &!
        SETUP_OBC ,        &!
        OBCS_DEALLOC    
   
   USE MOD_BCS, ONLY : &
        ! INODEQ            &  !!LOCAL FRESH WATER INFLOW NODES
        !,ICELLQ         &  !!LOCAL FRESH WATER INFLOW ELEMENTS 
        !,RIV_GL2LOC      &  !!GLOBAL TO LOCAL MAP OF FW NODES
        !,INOPNT          &  !!LOCAL NON-POINT SOURCE NODES
        !,QBC_TM         &  !!TIME MAP FOR RIVER DISCHARGE DATA 
        !,PNT_TM        &  !!TIME MAP FOR NON POINT SOURCE DATA
        !,VQDIST        &  !!DISCHARGE VERTICAL DISTRIBUTION
        !,WVQDIST       &  !!DISCHARGE VERTICAL DISTRIBUTION for point source
        !,DQDIS         &  !!WATER FLUX DISCHARGE DATA
        !,DWDIS           !&  !!WATER QUALITY DISCHARGE DATA
        !,WDIS          &  !!FRESH WATER QUALITY AT CURRENT TIME
        !,PWQDIS        &  !!Current WATER QUALITY at non-point source node
        !,PDWQDIS         &  !!WATER QUALITY DATA at non-point source node
        !,PQDIS           &  !!Current DISCHARGE at non-point source node
        !,PDQDIS         &  !!DISCHARGE at non-point source node
        BCS_DEALLOC

  USE MOD_SAV, ONLY    :        &
            !SPNAM,           &    !species name of the SAV
            !EPINAM,          &   !name of the epiphyte
            !NSAVM,           &      !global maximum number of SAV species in all of tracer control elements (TCE's) - polygon
                                !surrounding a given node
            
            !KHSTOX,         &   !Sulfide toxicity coefficient in calculating HS^- limitation on SAV growth (m^3/mole-S2)
            !KESAV,      &   !Effect on light attenuation by SAV leaf and stem, KESAV (gC/m^3)  
            !KEEPI,      &   !Shading effect on light on SAV shoots (leaf) due to Epiphytes ((m^2-leaf-surface/gDW-epiphyte)
            !KHEP,       &  !half saturation desnity of epiphyte (gC-Epiphyte/gC-leaf) for epiphyte growth on SAV leaf
            !NLIMEPI,    &  !Never used ???                      
             
            !KHNLEAF,    &  !Half saturation constant for N uptake by leaf (shoots) (gN/m^3)
            !KHNROOT,    &  !Half saturation constant for N uptake by root (gN/m^3)
            !KHPLEAF,    &  !Half saturation constant for P uptake by leaf (gP/m^3)
            !KHPROOT,    &  !Half saturation constant for P uptake by root (gP/m^3)
            !KHNEPI,     &  !Half saturation constant for N uptake by epiphytes (gN/m^3)
            !KHPEPI,     &  !Half saturation constant for P uptake by epiphytes (gN/m^3)
            !KHNPSAV,    &  !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for SAV
            !KHNPEPI,    &  !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for epiphytes

            !PMSAV,      &   !maximum production rate of SAV, gC/gDW/day  where DW is dry weight
            !PRSPSAVB,   &   !fraction of leaf production rate that is counted as photochemical respiration at salt << saltmax
                            !at salt >> saltmax, PRSPSAV approaches 1, and there will be no production of SAV leaf hence
                            !no photosynthesis for SAV growth
            !BMSAV,      &   !Basal metabolism rate of SAV (1/day)
            !BMTBRREF,   &   !Reference basal metabolism rate of tube [1/day] BMTUBE=BMTBRREF*f(T), i.e. function of temperature
            !FDOSR,      &   !Fraction of DOXG production discounted due to route of Leaf production to root
            !SALMAX,     &   !Salinity at which PRSPSAV is approaching 1 (maximum) from PRSPSAVB, [psu]
            !ALPHSAV,    &   !Initial slope of P-I curve for SAV [ (gC/gDW)/(E/m^2) ]
            !ALAC,       &   !m^2-leaf-area/gC-leaf (leaf area to leaf carbon ratio) 
                            !Note in 2002 version, it is the inverse ACLA: gC-leaf/m^2-leaf-area
            !SLSAV,      &   !???? Never used !Sloughing rate of SAV ?  1/day
            !ACDWSAV,    &   !Carbon to Dry Weight ratio of SAV (gC/gDW)
            !ANDWSAV,    &   !Nitrogen to Dry Weight ratio of SAV (gN/gDW)
            !APDWSAV,    &   !Phosphorus to Dry Weight ratio of Sav (gP/gDW)
            !ANCSAV,     &   !Nitrogen to Carbon ratio of SAV (gN/gC)             
            !APCSAV,     &   !P to C ratio of SAV (gP/gC)
            !HCAN,       &   !Canopy height (m), HCAN = ACAN + BCAN *( leaf+stem )
            !ACAN,       &   !basic canopy height (m)
            !BCAN,       &     !coefficient for calculating leaf (gC/m^2) and stem (gC/m^2) 
                             !contribution to canopy height,  [m^3/gC]

            !FNISAV,       &   !Fraction of NH4 Nitrogen generated due to SAV (leaf, stem etc) respiration (mortality)
            !FNLDSAV,      &   !Fraction of Nitrogen routed to LDON on SAV mortality
            !FNRDSAV,      &   !Fraction of Nitrogen routed to RDON on SAV mortality
            !FNLPSAV,      &   !Fraction of Nitrogen routed to LPON on SAV mortality
            !FNRPSAV,      &   !Fraction of Nitrogen routed to RPON on SAV mortality
            !FPISAV,       &   !Fraction of Phosphorus routed to PO4 in water column due to SAV (leaf,stem) respiration
            !FPLDSAV,      &   !Fraction of P routed to LDOP in water column due to SAV mortality
            !FPRDSAV,      &   !Fraction of P routed to RDOP in water column due to SAV mortality
            !FPLPSAV,      &   !Fraction of P routed to LPOP in water column due to SAV mortality
            !FPRPSAV,      &   !Fraction of P routed to RPOP in water column due to SAV mortality
            !FCLDSAV,      &   !Fraction of Carbon routed to LDOC in water column due to SAV mortality
            !FCRDSAV,      &   !Fraction of Carbon routed to RDOC in water column due to SAV mortality
            !FCLPSAV,      &   !Fraction of Carbon routed to LPOC in water column due to SAV mortality
            !FCRPSAV,      &   !Fraction of Carbon routed to RPOC in water column due to SAV mortality
            !FDOSAV,     &   !Fraction of Carbon routed to CO2 (oxygen consumption) in water column due to SAV mortality
                            !Note FDOSAV + FCLDSAV+ FCRDSAV + FCLPSAV + FCRPSAV = 1
                            !Wen Long: we really need to rename FDOSAV to FCO2SAV and make it equivalent to FDOSAV!

            !WSSSAV,       &   !Increase of fixed solids (ISS) settling by SAV leaf and stem, (m/day)/(gC/m^2)==>m^3/gC/day
            !WSLSAV,       &   !Increase of LPOM settling by SAV leaf and stem, m^3/gC/day
            !WSRSAV,       &   !Increase of RPOM settling by SAV leaf and stem, m^3/gC/day
            !WS1SAV,       &   !Increase of Algae group 1 settling by SAV leaf and stem, m^3/gC/day
            !WS2SAV,       &   !                        2
            !WS3SAV,       &   !                        3
            !WSUSAV,       &   !Increase of particulate biogenic silica settling rate by SAV leaf and stem, m^3/gC/day
                     
            !NDDI,         &     !number of degree days for entries of look up table of FPLEAF etc for Senescence of leaf etc
            !NSAVCELL,     &     !total number of grid cells (ICM boxes) that have SAV --> need to be changed to nodes of FVCOM
            !NSAVSPCM,     &   !maximum number of species in a cell
              
            !PMEPI,      &     !maximum photosynthesis production rate (gC/gCHLA/day) of epiphytes
            !BMEPI,      &     !basal metabolism loss rate of epiphyte (1/day)
            !PREPI,      &     !predation rate on epiphytes (gC-leaf/gC-epiphyte/day)
            !PRSPEPI,    &     !Photochemical respiration rate of epiphyte --> BMEP

            !ALPHEPI,    &     !slope of PI curve for epiphyte photothenesis (gC/gDW)/(E/m^2)
            !CCHLEPI,    &     !Carbon to Chlorophyll ratio of epiphyte (gC/gCHLA)
            !ANCEPI,     &     !N to C ratio of epiphytes (gN/gC)
            !APCEPI,     &     !P to C ratio of epiphytes (gP/gC)
            !ADWCEPI,    &      !ratio of total epiphyte dry weight to viable carbon (gDW/gC)
             
            !FCLDEPI,    &   !Fraction of C routed to LDOC due to basal metabolism of epiphyte
            !FCRDEPI,    &   !Fraction of C routed to RDOC due to basal metabolism of epiphyte
                            !WLong: missing FCLPEPI and FCRPEPI to routt to LPOC and RPOC ??? (seems not coded)
                            !FCLPEPI and FCRPEPI should be added (Wen Long)
                            !also missing FDOEPI compared to FDOSAV for SAV
                            !code assumes that no carbon will be routed to LPOP and RPOP due to metabolism
                            !and also FDOEPI=1-FCLDEPI-FCRDEPI, and FDOEPI is only implicitly used in calculating 
                            !DOEPIW
             
            !FNIEPI,     &   !Fraction of N routed to NH4  due to epiphyte basal metabolism
            !FNLDEPI,    &   !Fraction of N routed to LDON due to epiphyte basal metabolism
            !FNRDEPI,    &   !Fraction of N routed to RDON due to epiphyte basal metabolism
            !FNLPEPI,    &   !Fraction of N routed to LPON due to epiphyte basal metabolism
            !FNRPEPI,    &   !Fraction of N routed to RPON due to epiphyte basal metabolism
            !FPIEPI,     &   !Fraction of P routed to  PO4 due to basal metabolism of epiphyte
            !FPLDEPI,    &   !Fraction of P routed to LDOP due to basal metabolism of epiphyte
            !FPRDEPI,    &   !Fraction of P routed to RPOP due to basal metabolism of epiphyte
            !FPLPEPI,    &   !Fraction of P routed to LPOP due to basal metabolism of epiphyte
            !FPRPEPI,    &   !Fraction of P routed to RPOP due to basal metabolism of epiphyte             

            !FCLDPEP,    &   !Fraction of C routed to LDOC due to predation and sloughing on epiphytes
            !FCRDPEP,    &   !Fraction of C routed to RDOC due to predation and sloughing on epiphytes
            !FCLPPEP,    &   !Fraction of C routed to LPOC due to predation and sloughing on epiphytes
            !FCRPPEP,    &   !Fraction of C routed to RPOC due to predation and sloughing on epiphytes
            !FNIPEP,     &   !Fraction of Nitrogen routed to NH4  due to predation on epiphytes and Slough
            !FNLDPEP,    &   !Fraction of Nitrogen routed to LDON due to predation on epiphytes and Slough
            !FNRDPEP,    &   !Fraction of Nitrogen routed to RDON due to predation on epiphytes and Slough
            !FNLPPEP,    &   !Fraction of Nitrogen routed to LPON due to predation on epiphytes and Slough
            !FNRPPEP,    &   !Fraction of Nitrogen routed to RPON due to predation on epiphytes and Slough
            !FPIPEP,     &   !Fraction of P routed to PO4  due to predation on epiphytes and slough
            !FPLDPEP,    &   !Fraction of P routed to LDOP due to predation on epiphytes and slough
            !FPRDPEP,    &   !Fraction of P routed to RDOP due to predation on epiphytes and slough
            !FPLPPEP,    &   !Fraction of P routed to LPOP due to predation on epiphytes and slough
            !FPRPPEP,    &   !Fraction of P routed to RPOP due to predation on epiphytes and slough            

            !FPSR,          &    !Never used!!!! This is supposed to be the fraction of leaf production routed to root
                            !(replaced by FPROOT)

     
            !FPLEAF,     &     !fraction of production assigned to leaf growth
            !FPROOT,     &     !fraction of production assigned to root growth
            !FPSTEM,     &     !fraction of production assigned to stem growth
            !FPTUBER,    &     !fraction of production assigned to tuber growth
                            !WLong: FPLEAF +FPROOT+FPSTEM+FPTUBER = 1
                            !all growth is based on photosynthesis by leaf
                            !(so stem and tuber have no photosynthesis?? I thought they look green )
            !TRTBRLF,    &   !transfer rate from tuber to leaf (gC-Leaf/day/gC-Tuber)
                               
            !!Look up table of temperature control in a year on growth and mortality rates of SAV and epiphytes
            !FTPSAV,     &     !temperature control on SAV growth
            !FTRSAV,     &     !temperature control on SAV resipration 
            !FTPEP,      &     !temperature control on epiphyte growth rate
            !FTREP,      &     !temperature control on epiphyte resipration (basal metabolism)
            !FTPREP,    &   !temperature control on predation on epiphyte
                               
            !!fluxes to sediments by SAV mortality
            !FRPOCSAV,     &    !
            !FRPONSAV,     &    !
            !FRPOPSAV,  &    !3G partitioning of POC, PON, POP sources to sediments from SAV (non-dimensional)
  
            !!State variables of SAV model 
            !LEAF,       &     !Leaf biomass (gC/m^2) 
            !ROOT,       &     !Root biomass (gC/m^2)
            !STEM,       &     !Stem biomass (gC/m^2)
            !TUBER,      &     !Tuber biomass (gC/m^2)
            !EP,        &   !Epiphytes biomass (gC-epiphytes/gm^2-leaf)
                            !Note in 2002 version of ICM, EP's unit was (gC-epiphytes/gC-leaf)
                                 
            !SAVAREA,    &   !area of SAV coverage for each surface cell (m^2) --> needs to be changed for FVCOM grid                                 
            !SAVCELL,    &   !grid cell number for all cells that have SAV in it
                            !size 1xNSAVCELL --> need to change to element number or node number!!
            
            !NSAVSPC,    &    !(:) !number of SAV species in each grid cell
                            !size 1xNSAVCELL

            !SAVDPH,      &      !Depth of SAV (mean depth of the plant) for all species in each grid cell
                            !size, NSAVCELL x NSAVSPC(B) where B is cell number, B=1,...,NSAVCELL
            !SAVFRAC,    &      !coverage for each SAV cell and each vertical deph increment
                            !i.e. coverage is function of grid cell and SAV species

                            !Growth and mortality rates of SAV
            !PLEAF,        &    !Leaf growth rate (1/day), P (SAV shoots growth) in C-M 2001 paper's SH equation
            !BMLEAF,     &      !respiration rate of leaf of SAV (1/day) with temperature control 
            !BMSTEM,     &      !resipration rate of stem of SAV (1/day) with temperature control
            !BMROOT,     &      !respiration rate of root of SAV (1/day) with temperature control
            !BMTUBER,   &    !respiration rate of tuber of SAV (1/day) with temperature control                                      
     
            !Growth and mortality rates of SAV                                      
            !PEP,       &      !epiphytes growth rate (1/day)  = PMEPI*light_limitation*deinsity_limitation/CCHLEPI
            !BMEP,      &      !epiphytes basal metabolism rate (1/day) (also including photochemical respiration)
            !PREP,      &   !predation rate on epiphytes (1/day) = PREPI*EP
                            !where PREPI has unit gC-shoot/gC-epiphyte/day

            !SLSH,        &    !(:)      !epiphyte sloughing rate (1/day)
     
            !!Nutrient limitations on SAV and epiphytes growth
            !NLSAV,        &      !Nitrogen limiation on SAV growth (non-dimensional)
            !PLSAV,        &      !Phosphorus limitation on SAV growth (non-dimensional)
            !FNSEDSAV,     &      !fraction of nitrogen uptake by SAV due to nitrogen in sediment pore water
                            !Wen Long: I think we should further split into NH4 and NO3 
                            !          uptake by SAV from sediments
            !FPSEDSAV,     &      !fraction of phosphorus uptake by SAV due to P in sediment pore water
            !NLEPI,        &      !Nitrogen limitation on epiphytes growth (non-dimensional)
            !PLEPI,        &      !Phosphorus limitation on epiphytes growth (non-dimensional)
            !FHS,        &   !fresh water SAV species growth limitation due to sulfide HS^- 
                            !(in sediment anerobic layer) toxicity (dimensionless)

            !EPATN,      &     !Light attenuation by epiphyte
                            !e^(-Kep*Acla*Adwcep*EP) in Cerco-Moore 2001 paper for Ish eqn
                            !e^(-Kep*Adwcep*EP) in new version as EP changed to gC-Epi/m^2-leaf unit
            !SAVATN,     &     !Light attenuation by SAV (shoots)
                            ![1-e^(-Ksh*SH+KESS*HCAN)]/(Ksh*SH+KESS*HCAN) in C-M 2001 Iwc eqn
            !WATATN,     &     !water column attenuation of light(dimensionless) at canopy top
                            != exp(-KESS*ZTOCNPY) where ZTOCNPY is distance from surface
                            !                           to canopy top
            !FISH,       &     !light limitation function for LEAF (shoots) f(I)=Ish/sqrt(Ish^2+Ik^2)
            !FIEP,       &     !light limitation function for epiphytes, f(I)=Iep/sqrt(Iep^2+Ik^2)
            !NPPSAV,     &     !net leaf production per day (growth-basal metabolism) (gC/m^2/day)
            !NPPEPI,    &   !net epiphyte production per day (growth - basal metabolism) (gC-epi/gC-leaf/day)

            !!fluxes to water column
            !DOSAVW,     &     !DOXG source (flux) due to SAV (source-sink) (gO2/m^2/day ) 
                            !DOSAVW --> DTDO (gO2/m^3/sec)--> DTC (gO2/m^3/sec)
                            !-->DTM(gO2/m^3/sec)-->DTM*DLT-->C1(:,:,27)(gO2/m^3) !27'th constituents is DOXG
            !LDOCSAVW,     &     !LDOC flux due to SAV (gC/m^2/day)
            
                            !WL         !***ALSO used by TVDS.F ****
                            !LDOCSAVW --> DTDOC (gC/m^3/sec) -> DTC --> DTM --> DTM*DLT --> C1(:,:,9)
            !RDOCSAVW,     &      !RDOC flux due to SAV (gC/m^2/day)                                 
            !LPOCSAVW,     &      !LPOC flux due to SAV (gC/m^2/day) 
                            !WL           !***ALSO used by TVDS.F ****
            !RPOCSAVW,     &      !RPOC flux due to SAV (gC/m^2/day)
                            !WL            !***ALSO used by TVDS.F ****
            !NH4SAVW,      &      !NH4 flux to water column (respiration source - uptake sink by SAV growth) (gN/m^2/day)
                            !NH4SAVW (gN/m^2/day) --> DTNH4  (gN/m^3/sec) --> DTC --> DTM --> C1(:,:,13)
            !NO3SAVW,      &      !NO3 flux to water column (due to uptake by SAV, positive increasing NO3 in water) (gN/m^2/day)
            !LDONSAVW,     &      !LDON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf, basal metabolism of stem) (gN/m^2/day)
            !RDONSAVW,     &      !RDON flux                                     
            !LPONSAVW,     &      !LPON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem (gN/m^2/day)
            !RPONSAVW,     &      !RPON flux 
            !PO4SAVW,      &      !PO4 flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem and uptake by leaf) (gP/m^2/day)
            !LDOPSAVW,     &      !LDOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day)
            !RDOPSAVW,     &      !RDOP flux (gP/m^2/day)
            !LPOPSAVW,     &      !LPOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day) 
            !RPOPSAVW,     &      !RPOP flux (gP/m^2/day)
            !  DOEPIW,     &      !DOXG   flux (positive into water column) due to epiphytes (gO2/m^2/day)
            !LDOCEPIW,     &      !LDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
            !RDOCEPIW,     &      !RDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
            !LPOCEPIW,     &      !LPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
            !RPOCEPIW,     &      !RPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
            !NH4EPIW,     &      !NH4  flux (positive into water column) due to epiphytes (gN/m^2/day) (basal meabolism + predation + sloughing - uptake by photosynthesis)
            !NO3EPIW,     &      !NO3  flux
            !LDONEPIW,     &      !LDON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
            !RDONEPIW,     &      !RDON flux (positive into water column) (gN/m^2/day)
            !LPONEPIW,     &      !LPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
            !RPONEPIW,     &      !RPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
            !PO4EPIW,      &      !PO4  flux (positive into water column) due to epiphytes mortality (basal metabolism, salt toxicity, sloughing) and uptake (sink) (growth) (gP/m^2/day)                                    
            !LDOPEPIW,     &      !LDOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
            !RDOPEPIW,     &      !RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
            !LPOPEPIW,     &      !LPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
            !RPOPEPIW,     &    !RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)

            !!fluxes to sediments                                                                         
            !SEDPOCSAV,     &      !POC flux to sediments (positive into sediments) due to SAV (gC/m^2/day)
                            !SEDPOCSAV (gC/m^2/day) --> JPOC (mgC/m^2/day) -> POC1, POC2, POC3 (source of particulate organic carbon in sediments)
            !SEDPONSAV,     &      !PON flux to sediments (positive into sediments) due to SAV (gN/m^2/day)
                            !SEDPONSAV (gN/m^2/day) --> JPON (mgN/m^2/day) -> PON1, PON2, PON3 (source of particulate organic nitrogen in sediments)
            !SEDPOPSAV,     &      !POP flux to sedimetns (positive into sediments) due to SAV (gP/m^2/day)
                            !SEDPOPSAV (gP/m^2/day) --> JPOP (mgP/m^2/day) -> POP1, POP2, POP3 (source of particulate organic phosphorus in sediments)
            !SEDDOSAV,      &      !DOXG flux to water column (positive into water) due to SAV growth in sediments (source) and mortality of SAV in sediments (tuber and root basal metabolism)(sink)
                            !(gO2/m^2/day) (WLong: does this have to be to water ?? I thought the sinks should be in sediments directly for root metabolism)
            !SEDNH4SAV,     &      !NH4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
                            !(gN/m^2/day), SEDNH4SAV --> NH4T2TM1 (mg/N/m^3) decerase in sediments 
            !SEDNO3SAV,     &      !NO3 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
                            !(gN/m^2/day), SEDNO3SAV --> NO3T2TM1 (mg/N/m^3) decerase in sediments 
            !SEDPO4SAV,     &   !PO4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
            
            !!misc
            !DGRDAYS,       &      !Degree day calculated from temperature of water

            !!time average values 
                        
            !ALEAF,       &     !Average leaf biomass (gC/m^2)
            !AROOT,       &     !Average root biomass (gC/m^2)
            !ASTEM,       &     !Average stem biomass (gC/m^2)       
            !ATUBER,       &     !Average tuber biomass (gC/m^2)       
            !AEP,        &   !Average epiphyte biomass (gC-epi/gC-leaf)

            !APLEAF,       &   !Averaged leaf growth rate (1/day)
            !ABMLEAF,       &   !Averaged basal metabolism rate of leaf (1/day)
            !ABMTUBER,   &   !Average basal metabolism rate of tuber (1/day)

            !APEP,        &   !Average growth rate  for epiphytes (1/day)
            !ABMEP,       &   !Average basal metabolism rate for epiphytes (1/day)
            !APREP,        &    !Average predation rate on epiphytes (1/day) = APREP*EP

            !ASLSH,        &    !Average epiphute sloughing rate (1/day)
                
            !ANLSAV,        &   !Averaged N limitation on SAV 
            !APLSAV,        &   !Averaged P limitation on SAV
            !ANLEPI,        &   !Averaged N limitation on epiphytes
            !APLEPI,        &   !Averaged P limitation on epiphytes
            !AFNSED,        &   !Averaged fraction of uptake of N from sediment
            !AFPSED,        &   !Averaged fraction of uptake of P from sediment
            !AFHS,       &     !Averaged SAV growth limiation due to sulfide HS^-            
    
        !attenuation and light limitation
        !AEPATN,           &   !Averaged attenuation by epiphyte
        !AWATATN,           &   !Averaged water column attenuation of light at canopy top
        !AFISH,           &   !Averaged light limitation for LEAF (shoots)
        !AFIEP,        &    !Averaged light limitation for Epiphytes      

        !net production                                          
        !ANPPSAV,           &   !Averaged net leaf production per day (gC/m^2/day)
        !ANPPEPI,        &    !Averaged net epiphyte production per day (gC-epi/gC-leaf/day)

        !time averaged flux to water column
         !  ADOSAVW,      &      !average DOXG  flux due to SAV gO2/m^2/day
         ! ADOCSAVW,      &      !average DOC flux due to SAV
         ! APOCSAVW,      &      !average POC flux due to SAV
         ! ANH4SAVW,      &      !average NH4 flux due to SAV
         ! ANO3SAVW,      &      !average NO3 flux due to SAV
         ! ADONSAVW,      &      !average DON flux due to SAV              
         ! APONSAVW,      &      !average PON flux due to SAV
         ! APO4SAVW,      &      !average PO4 flux due to SAV
         ! ADOPSAVW,      &      !average DOP flux due to SAV
         ! APOPSAVW,      &      !average POL flux due to SAV
         !  ADOEPIW,      &      !DOXG  flux (positive into water column ) due to epiphytes (gO2/m2/day)
         ! ANH4EPIW,      &      !NH4 flux (positive into water column) due to epiphytes (gN/m2/day)
         ! ANO3EPIW,      &      !NO3 flux (positive into water column) due to epiphytes (gN/m2/day)
         ! APO4EPIW,      &      !PO4 flux (positive into water column) due to epiphytes (gP/m2/day)

         ! ADOCEPIW,      &      ! DOC flux to water column due to epiphytes (gC/m^2/day)
         ! APOCEPIW,      &      ! POC flux to water column due to epiphytes (gC/m^2/day)
         ! ADONEPIW,      &      ! DON flux to water column due to epiphytes (gN/m^2/day)
         ! APONEPIW,      &      ! PON flux to water column due to epiphytes (gN/m^2/day)                  
         ! ADOPEPIW,      &      ! DOP flux to water column due to epiphytes (gP/m^2/day)
         ! APOPEPIW,        &      ! POP flux to water column due to epiphytes (gP/m^2/day)

         ! ASEDDOSAV,       &   !DOXG  flux due to SAV,  positive leaving sediments
         ! ASEDPOCSAV,      &      !POC flux due to SAV, positive into sediments
         ! ASEDPONSAV,      &      !PON flux due to SAV, positive into sediments
         ! ASEDPOPSAV,      &      !POP flux due to SAV, positive into sediments
         ! ASEDNH4SAV,      &      !NH4 flux due to SAV, positive leaving sediemnts
         ! ASEDPO4SAV,        &    !PO4 flux due to SAV, positive leaving sediments
  
              !SAV_ALLOC,        &    !subroutine to allocate SAV variables
            SAV_DEALLOC,    &    !subroutine to deallocate SAV variables
        !SAV_LOADS_ALLOC,    &    !
        SAV_LOADS_DEALLOC,    &    !subroutine to deallocate SAV loads
            !SAV_READ,        &    !
            SAV_COMP            !subroutine that does SAV calculation
          

#  if defined (MULTIPROCESSOR)

   USE MOD_PAR, ONLY:   &      
            !EL_PID,    &        !(:)           !!PROCESSOR OWNER OF GLOBAL ELEMENT
            ELID,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT
            NLID,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE 
            ELID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
            NLID_X,    &        !(:)           !!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
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
#  endif

   USE MOD_DF, ONLY:            &
            DFEEDER,            &!
            DFEEDM1S_GL,        &!
            DFEEDM1S,            &!
            !DF_READ,            &!
            !DF_INIT,            &!
            DF_DEALLOC!,        &!
            !DF_ALLOC       
#if defined(SEDIMENT_DOM)			
 USE MOD_SED_DOM, ONLY:               &
			SED_DOM,SED_DOM_DEALLOC,   &
			SED_DOM_OUTPUT
			
 USE MOD_SED_DOM_EXCHANGE, ONLY:SED_DOM_FLAG,  SED_DOM_SHARE_ALLOCATE,  &
								SED_DOM_SHARE_DEALLOCATE, &
                                JNH4_GL_SHARE,    &  !NH4 flux in sediments (mg N / m^2/day)
                                    JNO3_GL_SHARE,    &  !NO3 flux in sediments (mg N / m^2/day)
                                    NO31_GL_SHARE,&        !sediment layer 1 NO3 conc
                                NO32_GL_SHARE,&        !sediment layer 2 NO3 conc
                                NH41_GL_SHARE,&        !sediment layer 1 NH4 conc
                                NH42_GL_SHARE       !sediment layer 2 NH4 conc!,  &
								!J_LDOC_SHARE, J_RDOC_SHARE, &    ! THE  FLUX OF DOM FROM THE SEDIMENTS.  POSITIVE if the water column is less than seds,
							    !J_LDON_SHARE, J_RDON_SHARE, &    ! i.e. water column is gaining DOM from sediments
								!J_LDOP_SHARE, J_RDOP_SHARE       ! g /m^2/day
#endif								 
 USE WC_DOM, ONLY : WC_DOM_ALLOCATE,WC_DOM_DEALLOCATE,      &
                        ZOOP_DOM, DOC, DON, DOP, & 

                        WC_CDOC1, WC_CDOC2, WC_CDOC3,   &
                        WC_NCDOC1, WC_NCDOC2, WC_NCDOC3,    &
                        XFLUX_CDOC1, XFLUX_CDOC2, XFLUX_CDOC3,   &
                        XFLUX_NCDOC1, XFLUX_NCDOC2, XFLUX_NCDOC3,    &   

						WC_CDON1, WC_CDON2, WC_CDON3, &
						WC_NCDON1, WC_NCDON2, WC_NCDON3,  &
						WC_CDOP1, WC_CDOP2, WC_CDOP3,  &
						WC_NCDOP1, WC_NCDOP2, WC_NCDOP3, MNLDOC,   &
						DTWCDOC1,DTWNCDOC1,DTWCDOC2,DTWNCDOC2,   &
						DTWCDOC3,DTWNCDOC3,   & 
						DTWCDON1,DTWNCDON1,DTWCDON2,DTWNCDON2,   &
						DTWCDON3,DTWNCDON3,   & 
						DTWCDOP1,DTWNCDOP1,DTWCDOP2,DTWNCDOP2,   &
						DTWCDOP3,DTWNCDOP3  ,  &
						DTWCDOC1_out,DTWNCDOC1_out,DTWCDOC2_out,DTWNCDOC2_out,   &
						DTWCDOC3_out,DTWNCDOC3_out,   & 
						DTWCDON1_out,DTWNCDON1_out,DTWCDON2_out,DTWNCDON2_out,   &
						DTWCDON3_out,DTWNCDON3_out,   & 
						DTWCDOP1_out,DTWNCDOP1_out,DTWCDOP2_out,DTWNCDOP2_out,   &
						DTWCDOP3_out,DTWNCDOP3_out  ,  &
   					  COAGN_OUT, COAGN, COAGC_OUT, COAGC,   &!
					  MNLDOC, &
					  MNLDOC1_OUT,MNLDOC2_OUT,MNLDOC3_OUT,   &!
					  MNLDOC1,MNLDOC2,MNLDOC3,   &!
					  MNLDON1_OUT,MNLDON2_OUT,MNLDON3_OUT,   &!
					  MNLDON1,MNLDON2,MNLDON3, &
					  PHOTOCHEM_DOC_ACCUM, PHOTOCHEM_DOC_ACCUM_OUT, &
					  MyVolumeTM1

						
USE MOD_NCDIO, ONLY : CDF_OUT, &
					  CDF_COUNTER, &
					  OUT_NETCDF, & 
					  CDF_INT, &
            CDF_START,CDF_END
            
!B Clark March 2020
#if defined(XFLUX_FLAG)            
 USE MOD_OUTPUT_FLUX, ONLY: EDGELINE_FLUX,  &
                            EDGELINE_FLUX_IN, &
                            WRITE_EDGELINE_FLUX, &
                            IOFILES_FLUX, &
                            DEALLOC_EDGELINE_FLUX

#endif                                     
						
   IMPLICIT NONE
# if defined (MULTIPROCESSOR)
   INCLUDE "mpif.h"  
# endif
   INTEGER :: J, JF, JG,         &
                    JCON,        &    !Index for constituents
                    K, L, N
   INTEGER :: I, II, III, ISF, JP, JNP
   INTEGER :: JS1, JS2, JS3, PRO2

   INTEGER :: MASTER_IINT             !Wen Long test synchronization of while loop
   REAL(SP) :: MASTER_ELTMS            !Wen Long test synchronization of while loop

!  KURT GLAESEMANN remove bug - ADJUSTL requires source and target to be same
!  CHARACTER(LEN=100) :: TEMPSTR
   CHARACTER(LEN=80 ) :: TEMPSTR

   INTEGER :: I1,I2
   REAL(SP)    :: ITT,TTMP1,TTMP2
   
   REAL(SP) :: MyVolume
   
!RGl added WQ_DT, above
!WLong moved this to mod_sed.F
! Restart printout variables below - KURT GLAESEMANN
   !REAL(SP),ALLOCATABLE :: CTEMP_GL(:)
   !REAL(SP),ALLOCATABLE :: CPOP_GL(:,:)
   !REAL(SP),ALLOCATABLE :: CPON_GL(:,:)
   !REAL(SP),ALLOCATABLE :: CPOC_GL(:,:)
   !REAL(SP),ALLOCATABLE :: CPOS_GL(:)
   !REAL(SP),ALLOCATABLE :: PO4T2TM1S_GL(:)
   !REAL(SP),ALLOCATABLE :: NH4T2TM1S_GL(:)
   !REAL(SP),ALLOCATABLE :: NO3T2TM1S_GL(:)
   !REAL(SP),ALLOCATABLE :: HST2TM1S_GL(:)
   !REAL(SP),ALLOCATABLE :: CH4T2TM1S_GL(:)
   !REAL(SP),ALLOCATABLE :: CH41TM1S_GL(:)
   !REAL(SP),ALLOCATABLE :: SO4T2TM1S_GL(:)
   !REAL(SP),ALLOCATABLE :: SIT2TM1S_GL(:)
   !REAL(SP),ALLOCATABLE :: BENSTRTM1S_GL(:)
   
   !WLong moved this to mod_ba.F
   !REAL(SP),ALLOCATABLE :: BBM_GL(:)

!WLong moved to mod_df.F
!   REAL(SP),ALLOCATABLE :: DFEEDM1S_GL(:)

!WLong moved SFEED_GL to mod_sf.F
!   REAL(SP),ALLOCATABLE :: SFEED_GL(:,:)

!WLong Added logicals for writing station and history outputs
   LOGICAL :: OPEN_HIS,PRINT_HIS
   LOGICAL :: OPEN_STN,PRINT_STN
   INTEGER(4) :: CTR_HIS, CTR_STN    !record counter for station and history outputs

   INTEGER(4) :: NTSTN_OLD, NTSTN
   INTEGER(4) :: NTHIS_OLD, NTHIS
 
!   CHARACTER(1024):: fnameprefix, fnameext
   CHARACTER(1024):: HISFN_FILENUMBER
            
!Wen Long time the simulation here

   REAL(8) :: WTIME_MOD=0.d0
   REAL(8) :: WTIME_SYS=0.d0
   REAL(8) :: TOTAL_TIME=0.d0
   REAL(4) :: ETIME,DTIME    !ETIME and DTIME are intrinsic functions to calculate time
   REAL(4) :: ELAPSED(2)
   
!Wen Long and Laura B tried to profile program using ELAPSED_PROFILE(2)
   INTEGER :: NPROFILE=0
   REAL(4) :: ELAPSED_PROFILE(2)
   REAL(8) :: WTIME_MOD_PROFILE=0.d0
   REAL(8) :: WTIME_SYS_PROFILE=0.d0
   REAL(8) :: TOTAL_TIME_PROFILE=0.d0

   CHARACTER(LEN=1024) :: TITLE_STN

   CHARACTER(LEN=1024)::COMMAND_STR,SLEEPTIME_STR
   REAL(SP) :: SLEEPTIME

   LOGICAL :: END_RUN_CHILD =.FALSE.
   INTEGER :: ISNDTAG =0 , IRCVTAG =0

!************************************************************************
!**                       Program Initialization                       **
!************************************************************************

!==============================================================================!
!   SETUP PARALLEL ENVIRONMENT                                                 !
!==============================================================================!

   SERIAL = .TRUE. 
   PAR    = .FALSE. 
   MSR    = .TRUE.
   MYID   = 1
   NPROCS = 1
#  if defined (MULTIPROCESSOR)
   CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR)
#  endif

!--IMPORT CASENAME FROM COMMAND LINE  
   IF(MSR)THEN
     CALL GETARG(1,TEMPSTR)

     IF(LEN_TRIM(TEMPSTR) == 0)THEN
       WRITE(*,*) 'PLEASE PROVIDE CASENAME ON COMMAND LINE'
       WRITE(*,*) 'STOPPING...'
       CALL PSTOP
     ENDIF
     CASENAME = ADJUSTL(TEMPSTR)
   ENDIF        

   IF(MSR)THEN
      TOTAL_TIME=TOTAL_TIME+DTIME(ELAPSED)  !CALL DTIME() and return user time in ELAPSED(1), system time in ELAPSED(2)
      WTIME_MOD=WTIME_MOD+ELAPSED(1)
      WTIME_SYS=WTIME_SYS+ELAPSED(2)
   ENDIF


#  if defined (MULTIPROCESSOR)
   CALL MPI_BCAST(CASENAME,80,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
#  endif

   CALL INIT_FILE_INFO()    !  Setup MOD_WQM file information
   
   CALL INPUTS()            !  Initialze 

# if defined(XFLUX_FLAG)   
   CALL IOFILES_FLUX
# endif   
  
   IF(SEDIMENT_CALC)THEN
   
!        CALL SED_INIT()          !  Setup sediment diagenesis module (including SED_ALLOC)
        !CALL SED_INIT2()         !  initialize the values from input read by INPUTS() above
        IF(DFEEDER)THEN
            CALL SED_DF_EXCHANGE_ALLOC
        ENDIF
#if defined(SEDIMENT_DOM)       
		IF(SED_DOM)THEN
		
      CALL SED_DOM_SHARE_ALLOCATE  ! B CLARK SED_DOM UPDATE
			SED_DOM_FLAG = SED_DOM  ! logical to pass to mod_wc_dom
	    ENDIF
#endif			
		
        CALL SED_SF_EXCHANGE_ALLOC
        
        IF(SAV_CALC)THEN
            CALL SED_SAV_EXCHANGE_ALLOC
        ENDIF
   ENDIF
      
   CALL SETUP_OBC
   
  !CALL SET_BNDRY              !WLong: Was in FVCOM, but now removed, 
                            !which should be used for calculating N_ICELLQ(:,1:2)
                            ! maybe we should add it back so that WQM can take edge type of nutrient flux
                            !   

   CTR_STN = 0              !counter for station output records
   CTR_HIS = 0              !counter for history output records
   PRINT_HIS = .FALSE.
   PRINT_STN = .FALSE.
   OPEN_HIS  = .FALSE.
   OPEN_STN  = .FALSE. 

   !counting numbe of steps past t_his_start and t_stn_start
   !this may need to be different for RESTART runs, basically we should read these numbers from
   !restart file there 
   NTHIS_OLD = FLOOR((JDAY-t_his_start)/t_his_dlt)-1  !number of steps for history output
   NTSTN_OLD = FLOOR((JDAY-t_stn_start)/t_stn_dlt)-1  !number of steps for station output
   NTHIS=0
   NTSTN=0 
   
!below RGl added WQ_DT
   WQ_DT =DLT
   T(:,:)    = C2(:,:,1)
   SALT(:,:) = C2(:,:,2)
   SSI(:,:)  = C2(:,:,3)
   B1(:,:)   = C2(:,:,4)
   B2(:,:)   = C2(:,:,5)
   B3(:,:)   = C2(:,:,6)
   SZ(:,:)   = C2(:,:,7)
   LZ(:,:)   = C2(:,:,8)
   
   WC_CDOC1(:,:) = C2(:,:,9)
   WC_NCDOC1(:,:) = C2(:,:,10)
  
   !LDOC(:,:) = C2(:,:,9)
   !RDOC(:,:) = C2(:,:,10)
   LPOC(:,:) = C2(:,:,11)
   RPOC(:,:) = C2(:,:,12)
   NH4(:,:)  = C2(:,:,13)
   NO3(:,:)  = C2(:,:,14)
   UREA(:,:) = C2(:,:,15)
   
   WC_CDON1(:,:) = C2(:,:,16)
   WC_NCDON1(:,:) = C2(:,:,17)
   !LDON(:,:) = C2(:,:,16)
   !RDON(:,:) = C2(:,:,17)
   LPON(:,:) = C2(:,:,18)
   RPON(:,:) = C2(:,:,19)
   PO4(:,:)  = C2(:,:,20)
   
   WC_CDOP1(:,:) = C2(:,:,21)
   WC_NCDOP1(:,:) = C2(:,:,22)
   !LDOP(:,:) = C2(:,:,21)
   !RDOP(:,:) = C2(:,:,22)
   LPOP(:,:) = C2(:,:,23)
   RPOP(:,:) = C2(:,:,24)
   PIP(:,:)  = C2(:,:,25)
   COD(:,:)  = C2(:,:,26)
   DOXG(:,:)   = C2(:,:,27)
   SIUPB(:,:)   = C2(:,:,28) 
   SIAT(:,:)   = C2(:,:,29)
   PIB1(:,:) = C2(:,:,30)
   PIB2(:,:) = C2(:,:,31)
   PIB3(:,:) = C2(:,:,32)
   
   WC_CDOC2(:,:) = C2(:,:,33)
   WC_NCDOC2(:,:) = C2(:,:,34)

   WC_CDOC3(:,:) = C2(:,:,35)
   WC_NCDOC3(:,:) = C2(:,:,36)


   WC_CDON2(:,:) = C2(:,:,37)
   WC_NCDON2(:,:) = C2(:,:,38)
   WC_CDON3(:,:) = C2(:,:,39)
   WC_NCDON3(:,:) = C2(:,:,40)
   WC_CDOP2(:,:) = C2(:,:,41)
   WC_NCDOP2(:,:) = C2(:,:,42)
   WC_CDOP3(:,:) = C2(:,:,43)
   WC_NCDOP3(:,:) = C2(:,:,44)
   
   IINT=0

!************************************************************************
!**                          Begin Simulation                          **
!************************************************************************
   D   = H+ELNC1
   DO I = 1,MLOC
    DO K=1,KBM1
     MyVolumeTM1(I,K)=(D(I)*DZ(K)*ART1(I)) 
    !  write(*,*)'Initial VOlume = ',MyVolumeTM1(I,K)
    !  write(*,*)'Initial DZ = ',DZ(K)
    !  write(*,*) 'Initial D =',D(I)
    !  write(*,*)'Initial ART1 = ',ART1(i)
    ENDDO
  ENDDO


   DO WHILE(.NOT.END_RUN)
	 NPROFILE=0
     IINT=IINT+1
     THOUR = ELTMS/3600.0
     
     IF(MSR)WRITE(*,*)'SIMULATION TIME: ',THOUR,' HOURS',JDAY,'JDAYS'
     
	 
     
!*******  time-varying kinetics data
! RGL need to modify this for hourly read-in of met data
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 
     IF (JDAY >= NXTVD) CALL TVDS (NXTVD)
	 
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," TVDS "
	 ENDIF
	 
!******* Update flows
!     IF (NINT(ELTMS) >= NXHYD)  THEN
!RGL added below

# if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       !Wen Long
       !Here we should make sure ELTMS is the same for all processors 
       !as we are going to have a broadcasting from master process
       !and only master process is reading the hydrodynamics.
       !It will be disastrous if time-sensitive information (hydrodynamic update) is broadcast
       !to child processes who are not yet up to the same page in time (while loop)
       !This may not occur very often since the master process is in general doing more work and 
       !slower than child processors, however, in case of meeting the odds, 
       !a strict check point here is definitely necessary.
       !
       !In theory, we could do an MPI_BARRIER() here, but that will slow down the program
       !and also MPI_BARRIER() is not well suited  for being in while-loops (where number of looping
       !interations is unknown apriori owing to variable timestep)

       !Check loop index IINT

        IF(MSR)THEN                !Master calculates its own IINT
            MASTER_IINT=IINT       !Then broadcast it to all processes
            MASTER_ELTMS=ELTMS
        ENDIF
        CALL MPI_BCAST(MASTER_IINT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(MASTER_ELTMS,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)   !Notice MPI_REAL8

        IF(MASTER_IINT /= IINT)THEN !child checks for consistency
          WRITE(*,*)'MYID=',MYID,'. Oops, inconsistent time step here with master process'
        ENDIF

        IF(MASTER_ELTMS /= ELTMS )THEN !child checks for consistency
          WRITE(*,*)'MYID=',MYID,'. Oops, inconsistent time, ELTMS=', ELTMS
        ENDIF

       !Check time ELTMS  (seconds elapsed from begining of time)

     ENDIF
# endif

   IF (MOD(ELTMS,AHMDLT).eq.0)  THEN
      !Save previous time step results
       UNC1  = UNC2
       VNC1  = VNC2
       WNC1  = WNC2
       WTSNC1 = WTSNC2
       UARD_OBCNNC1=UARD_OBCNNC2
       XFLUX_OBCNC1=XFLUX_OBCNC2
       DTFANC1=DTFANC2
       KHNC1 = KHNC2
       ELNC1 = ELNC2
       TNC1  = TNC2
       SNC1  = SNC2
	   
	   
		IF(MSR)THEN
			TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
			NPROFILE=NPROFILE+1
		ENDIF
	 
     ! ELNC2_TM1 = ELNC2  ! B Clark added this here so that we can get the old surface elevation to know the difference between the two
      
      !update the master process' hydrodynamics
       IF(MSR)CALL HYDRO()
	   IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," HYDRO"
		ENDIF
	 
      !update the child process' hydrodynamics by broadcasting from master
# if defined (MULTIPROCESSOR)  
		
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 
	!
    !Method A): scatter the hydrodynamics from master proc to other procs
    !           This is supposed to be faster as the message size is smaller. However 
    !

	!!Wen Long: this does not work well yet
    !!  now BROAD_CAST_HYDRO has been downsized to only deal with UARD_OBCNNC2 and XFLUX_OBCNC2 and ELL_GL, DTFAL_GL
    !  CALL BROADCAST_HYDRO_REDUCED(0,UARD_OBCNNC2,XFLUX_OBCNC2,DTFANC2,ELNC2)
    !!  scatter the rest of variables instead of broadcasting (scattering is faster with smaller message size)
    !  CALL SCATTER(0,NLOfC,NTLOC,NGL,KBM1,MYID,NPROCS,EMAP,EC,ELID,ELID_X,  UL_GL,UNC2)   !velocity UL_GL->UNC2(0:NTLOC,KB), VL_GL->VNC2(0:NTLOC,KB)
	!  CALL SCATTER(0,NLOC,NTLOC,NGL,KBM1,MYID,NPROCS,EMAP,EC,ELID,ELID_X,  VL_GL,VNC2)   !velocity UL_GL->UNC2(0:NTLOC,KB), VL_GL->VNC2(0:NTLOC,KB)	  
    !  CALL SCATTER(0,MLOC,MTLOC,MGL,KB,  MYID,NPROCS,NMAP,NC,NLID,NLID_X,WTSL_GL,WTSNC2) !velocity WTSL_GL->WTSNC2(0:MTLOC,KB) 
    !!                                                                                    !      and KHL_GL->KHNC2(0:MTLOC,KB)
	!  CALL SCATTER(0,MLOC,MTLOC,MGL,KB,  MYID,NPROCS,NMAP,NC,NLID,NLID_X, KHL_GL,KHNC2)  !velocity WTSL_GL->WTSNC2(0:MTLOC,KB) 	  
    !  CALL SCATTER(0,MLOC,MTLOC,MGL,KBM1,MYID,NPROCS,NMAP,NC,NLID,NLID_X,  SL_GL,SNC2)   !SL_GL->SNC2(0:MTLOC,KB)  TL_GL->TNC2(0:MTLOC,KB)
	!  CALL SCATTER(0,MLOC,MTLOC,MGL,KBM1,MYID,NPROCS,NMAP,NC,NLID,NLID_X,  TL_GL,TNC2)   !SL_GL->SNC2(0:MTLOC,KB)  TL_GL->TNC2(0:MTLOC,KB)	  
    !!                                                                                                ! Wen Long: seems SNC2 and TNC2 should have 
     
    !Method B) broadcast hydrodynamics from master proc to other procs
       CALL BROADCAST_HYDRO(0,UNC2,VNC2,WTSNC2,UARD_OBCNNC2,XFLUX_OBCNC2,DTFANC2,KHNC2,ELNC2,TNC2,SNC2)
  
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A8)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," BCHYDRO"
	 ENDIF
# endif

    ENDIF

!--Time interpolation of physical fields
     UU  = UNC1
     VV  = VNC1
     WTS = WTSNC1
     UARD_OBCN=UARD_OBCNNC1
     XFLUX_OBC=XFLUX_OBCNC1
     DTFA=DTFANC1  
     KH  = KHNC1
     T   = TNC1
     SALT= SNC1
     EL  = ELNC1
     D   = H+EL

!RGL uncommented below
   C2(:,:,1)  = T(:,:)
   C2(:,:,2)  = SALT(:,:)

!******* Dead sea case
    !  IF(.NOT. FLOW)THEN
    !    UU  = 0.0
    !    VV  = 0.0
    !    WTS = 0.0
    !  ENDIF

     IF(.NOT. Z_DIFFUSION) KH = 0.0  

     ELTMS1 = REAL(ELTMS,8)

!******* Time-variable solar radiation

!--- TIME OF DAY (SECONDS)
     RELTMS = ELTMS
     TODS = MOD(RELTMS,86400.)  ! calcualtes remainder of model time and seconds in aday to get time of day in seconds

! KURT GLAESEMANN added MSR to prints
!     if(MSR)WRITE(*,*)'TODS', TODSrite
!     if(MSR)WRITE(*,*)'TTSS', TTSS
     IF (TODS < TTSS) THEN
       I0 = 0.
     ELSE
!--- DAYS SINCE SUNRISE (DAYS)
       DSSR = (TODS-TTSS)/86400.
	   ! FD is fixed in tvds so that light is always being calculated, but it is turned on or off
	   !based on the time of day and the surface reflectance in mod_owq
	   !this is based off of the lat of the model

! KURT GLAESEMANN - change 3.14.. to "PI", for more accuracy

       I0 = PI*IT*SIN(PI*DSSR/FD)/2./FD
       I0 = MAX(I0,0.)
! KURT GLAESEMANN add MSR to print out

     ENDIF

!--- HARDWIRE LIGHT HERE
!        I0 = 1.0        
 !    write(*,*)'light extinction ' ,light_extinction  
!******* Light Attenuation
    IF(LIGHT_EXTINCTION)THEN !Wen Long added this flag of LIGHT_EXTINCTION
        IF (SAV_CALC .OR. ALGAE_CALC) THEN
            IF (JDAY >= NXKEI) THEN
                NXKEI = NXKEI + INTKE/24. 

                  IF (I0 > 0. .OR. IAVG(1,1) > 0.) CALL LGHT_ATTN    !calculate light attenuation 
            ENDIF
        ENDIF        
    ENDIF


!******* Determine cell quota for this time step
!WLong and LB reactivated calculation of Q1,Q2,Q3
     DO K=1,KBM1
       DO I=1,MLOC
         Q1(I,K) = PIB1(I,K)/(B1(I,K)+1.0E-9)  !Phosphorus to Carbon Ration in Alg 1 (gP/gC)
         Q2(I,K) = PIB2(I,K)/(B2(I,K)+1.0E-9)  !Phosphorus to Carbon Ration in Alg 2 (gP/gC)
         Q3(I,K) = PIB3(I,K)/(B3(I,K)+1.0E-9)    !Phosphorus to Carbon Ration in Alg 1 (gP/gC)
         
         IF(SEDIMENT_CALC.AND.ACC(30)=='OFF')THEN
                Q1(I,K)=0.0243
         ENDIF
         IF(SEDIMENT_CALC.AND.ACC(31)=='OFF')THEN
                Q2(I,K)=0.0243
         ENDIF
         IF(SEDIMENT_CALC.AND.ACC(32)=='OFF')THEN
                Q3(I,K)=0.0243
         ENDIF
         
       ENDDO  
     ENDDO

	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 
!******* Water column kinetics
!     IF (TEMPERATURE_CALC) CALL TEMPER(DTC(0,1,1))
!tykim-this makes change in temperature zero - need to comment out if want to calc temp
        DTC(:,:,1)=0.0
!In next section first subscript changed to 0 per KURT 
		
	 
       IF (SOLIDS_CALC)      CALL SOLIDS(DTC(0,1,3),FLUXS(0,1,1))  !WLong and LB, we should move this into the IF block below as well??

       IF(MOD(ELTMS,WQ_DT).eq.0) THEN
! !RGL added sub cycling of ecosystem model, above
! !WL  IF(MSR) WRITE(*,*)'calling ecosystem model'

!     !    IF (ZOO_CALC)    THEN     
! 	  !     CALL ZOOPL(DTC(0,1,7),DTC(0,1,8))
! 		!   CALL ZOOP_DOM
! 	  !  ENDIF
	   
       IF (ALGAE_CALC)         &
        CALL ALGAE(DTC(0,1,4),DTC(0,1,5),DTC(0,1,6),FLUXS(0,1,2),      &
                  FLUXS(0,1,3),FLUXS(0,1,4))
       IF (CARBON_CALC)     THEN 
	       CALL PARTICULATE_CARBON(DTC(0,1,11),DTC(0,1,12) , FLUXS(0,1,5))   ! DTC(0,1,9) DTC(0,1,10)  ! old DOC
	       CALL DOC
		     DTC(:,:,9) = DTWCDOC1
			   DTC(:,:,10) = DTWNCDOC1
			   DTC(:,:,33) = DTWCDOC2
			   DTC(:,:,34) = DTWNCDOC2
			   DTC(:,:,35) = DTWCDOC3
			   DTC(:,:,36) = DTWNCDOC3
	     ENDIF 
       IF (NITROGEN_CALC)     THEN
	        CALL NITROG(DTC(0,1,13) ,DTC(0,1,14), &  !,DTC(0,1,16)  DTC(0,1,17), old ldon and rdon
                                      DTC(0,1,18),DTC(0,1,19),  &
                                       FLUXS(0,1,6))
	        CALL DON
          DTC(:,:,16) = DTWCDON1    
			    DTC(:,:,17) = DTWNCDON1   
			    DTC(:,:,37) = DTWCDON2   
			    DTC(:,:,38) = DTWNCDON2 
			    DTC(:,:,39) = DTWCDON3
			    DTC(:,:,40) = DTWNCDON3
			 ENDIF
		
!       IF (PHOSPHORUS_CALC)    THEN  
!           CALL PHOSPH(DTC(0,1,20),DTC(0,1,23),&  !DTC(0,1,21),DTC(0,1,22), these are old ldop and rdop
!                   DTC(0,1,24),DTC(0,1,25),DTC(0,1,30),DTC(0,1,31),    &
!                   DTC(0,1,32),FLUXS(0,1,7),FLUXS(0,1,8),FLUXS(0,1,9))
		
!	       CALL DOP 
!		   
!		    DTC(:,:,21) = DTWCDOP1
!				DTC(:,:,22) = DTWNCDOP1
!				DTC(:,:,41) = DTWCDOP2
!				DTC(:,:,42) = DTWNCDOP2
!				DTC(:,:,43) = DTWCDOP3
!				DTC(:,:,44) = DTWNCDOP3
				
!	     ENDIF
	   
	   
         IF (COD_CALC)         CALL CODMND(DTC(0,1,26))
	   
         IF (OXYGEN_CALC)      CALL OXYGEN(DTC(0,1,27))
        IF (SILICA_CALC)      CALL SILICA(DTC(0,1,28),DTC(0,1,29),FLUXS(0,1,10))

      !  SFENT = -9999
      !  IF (SFEEDER) THEN
      !    DO I=1,NSPECIES
      !      CALL SUSPFEED(I,DTC(0,1,4), DTC(0,1,5), DTC(0,1,6), DTC(0,1,30), &
      !                    DTC(0,1,31),DTC(0,1,32),DTC(0,1,11),DTC(0,1,23), &
      !                    DTC(0,1,18),DTC(0,1,12),DTC(0,1,24),DTC(0,1,19), &
      !                    DTC(0,1,3), DTC(0,1,28),DTC(0,1,29),DTC(0,1,20))
      !    ENDDO
      !  ENDIF
	

!     *** Convert water column change in concentration to change in mass
!
      ! DO JCON=1,NAC
      !  II=AC(JCON)
	!	write(*,*)'II = ',II
       ! DO K=1,KBM1
        !  DO I=1,MLOC
         ! DTM(I,K,II) = DTC(I,K,II)       !*V1(I,K)
            DTM= DTC    !*V1(I,K)

         ! ENDDO  
       ! ENDDO
      ! ENDDO  

      ENDIF

	 
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," WQKIN"
	 ENDIF
	 
	  
! END of RGL sub-cycling of water quality equations
!***** Horizontal diffusion
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 
	
     IF (XY_DIFFUSION .AND. XYDF < 0.0) CALL VISCOF_H
	IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," VISCOFH"
	 ENDIF


	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF

     CALL EXCHANGE_WQM(C2)
	 
	  IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," EXCHG1"
	 ENDIF
	 
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 
	 
	 CALL DENS2 ! Added by B Clark to calculate density in the model, June 2018
	 
     CALL ADV_WQM
	 
	  IF(MSR)THEN
	
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," WQADV"
	 ENDIF
	 
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 

	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF

     CALL VDIF_WQM(C2F)
	  IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," WQVDIF"
	 ENDIF
	 
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 
!    CALL EXCHANGE_WQM                              !Interprocessor Exchange   !
     CALL BCOND_WQM                                 !Boundary Conditions       !
	  IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," WQBCOND"
	 ENDIF
	 
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 
! KURT GLAESEMANN - have to exchange again, after boundary conditions
! because boundary condition code does NOT know about the buffer regions
     CALL EXCHANGE_WQM(C2F)                         !Interprocessor Exchange   !
	  IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," EXCHG2"
	 ENDIF
	
	 
     C2 = C2F
     ET  = EL
     DT  = D
 
#if(1) /* Laura and Wen moved sediment calculation after WQ update */
!******* Sediment kinetics

     SEDENT = -9999

     IF (SEDIMENT_CALC) THEN

        !set initial condition using steady state calculation 
        IF((STEADY_STATE_SED_IC).AND.(JDAY<TMSTRT+DLT/86400.0))THEN
        
            CALL SED_INIT_ICI       !WL: make sure the sediment CTEMP get the initial condition of temperature from
                                    !just updated water column (end of first time step)
                                    !(this is becaluse excel version uses water column temperature AND DO NOT REALLY
                                    !HAVE INITIAL TEMPERATURE, EXCEL VERSION DOESN NOT HAVE A TEMPERATURE EQUATION)


            CALL POM_ACCUMUL
            
            IF(AGGREGATE_POM_FLUX) THEN !LB: Option to accumulate all labilities from water column and then distribute them to labilities in sediments
                DO I=1,MLOC
                !Accumulate all labilities into JPOM(I,1)
                  JPOCaccum(I,1) = JPOCaccum(I,1) + JPOCaccum(I,2) + JPOCaccum(I,3)
                  JPONaccum(I,1) = JPONaccum(I,1) + JPONaccum(I,2) + JPONaccum(I,3)
                  JPOPaccum(I,1) = JPOPaccum(I,1) + JPOPaccum(I,2) + JPOPaccum(I,3)
                  
                !Divide JPOM(I,1) into labilities according to FRACL1,FRACL2,FRACL3 (given in BFI input)
                  JPOC(I,1) = JPOCaccum(I,1) * FRACL1
                  JPON(I,1) = JPONaccum(I,1) * FRACL1
                  JPOP(I,1) = JPOPaccum(I,1) * FRACL1

                  JPOC(I,2) = JPOCaccum(I,1) * FRACL2
                  JPON(I,2) = JPONaccum(I,1) * FRACL2
                  JPOP(I,2) = JPOPaccum(I,1) * FRACL2
                    
                  JPOC(I,3) = JPOCaccum(I,1) * FRACL3
                  JPON(I,3) = JPONaccum(I,1) * FRACL3
                  JPOP(I,3) = JPOPaccum(I,1) * FRACL3
                    
                  JPOS(I) = JPOSaccum(I)
                ENDDO
            
            ELSE
                DO I=1,MLOC
                  DO K=1,3
                    JPOC(I,K) = JPOCaccum(I,K)
                    JPON(I,K) = JPONaccum(I,K)
                    JPOP(I,K) = JPOPaccum(I,K)
                  ENDDO
                JPOS(I) = JPOSaccum(I)
                ENDDO
            ENDIF
        
            CALL SED_CALC(.TRUE.)  !True means use steady state solution as initial condition in sediment module
            CALL SED_INIT2         !set the initial conditions using the calcualated value            
!            CALL WQM_OUTPUT(NTHIS_OLD, NTSTN_OLD, NTHIS,NTSTN,CTR_HIS,CTR_STN)            
            
            !set accumulative fluxes back to zero
            JPOCaccum = 0.0
            JPONaccum = 0.0
            JPOPaccum = 0.0
            JPOSaccum = 0.0
        ENDIF
        
        
	IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	ENDIF
	
	
        !LB: Here we call POM_ACCUMUL to accumulate organic matter in JPOM variables during NDTSED time steps.
        !    Then, we divide JPOM/NDTSED to calculate average JPOM in that period. 
        !    The resulting JPOM is used in SED_CALC, which is only called every NDTSED timesteps
        
        CALL POM_ACCUMUL
         
	
	 
        !if(MSR)write(*,*)'LBnote: NIT=',NIT,'  NDTSED=',NDTSED
        !IF(MSR)WRITE(*,*)'LBnote: MOD(NIT,NDTSED)=',MOD(NIT,NDTSED)
        IF(MOD(NIT,NDTSED) .EQ. 0 .AND. NIT.NE.0) THEN   !LB: if NIT is multiple of NDTSED, then average JPOM and call sed_calc
            !if(MSR) write(*,*)'LBnote call sed_calc'
        
            IF(AGGREGATE_POM_FLUX) THEN !LB: Option to accumulate all labilities from water column and then distribute them to labilities in sediments
                DO I=1,MLOC
                !Accumulate all labilities into JPOM(I,1)
                  JPOCaccum(I,1) = JPOCaccum(I,1) + JPOCaccum(I,2) + JPOCaccum(I,3)
                  JPONaccum(I,1) = JPONaccum(I,1) + JPONaccum(I,2) + JPONaccum(I,3)
                  JPOCaccum(I,1) = JPOPaccum(I,1) + JPOPaccum(I,2) + JPOPaccum(I,3)

                 
                !Divide JPOM(I,1) into labilities according to FRACL1,FRACL2,FRACL3 (given in BFI input)
                  JPOC(I,1) = JPOCaccum(I,1) * FRACL1 / REAL(NDTSED)
                  JPON(I,1) = JPONaccum(I,1) * FRACL1 / REAL(NDTSED)
                  JPOP(I,1) = JPOPaccum(I,1) * FRACL1 / REAL(NDTSED)

                  JPOC(I,2) = JPOCaccum(I,1) * FRACL2 / REAL(NDTSED)
                  JPON(I,2) = JPONaccum(I,1) * FRACL2 / REAL(NDTSED)
                  JPOP(I,2) = JPOPaccum(I,1) * FRACL2 / REAL(NDTSED)
                    
                  JPOC(I,3) = JPOCaccum(I,1) * FRACL3 / REAL(NDTSED)
                  JPON(I,3) = JPONaccum(I,1) * FRACL3 / REAL(NDTSED)
                  JPOP(I,3) = JPOPaccum(I,1) * FRACL3 / REAL(NDTSED)
                    
                  JPOS(I) = JPOSaccum(I) / REAL(NDTSED)
                ENDDO
            
            ELSE
                DO I=1,MLOC
                  DO K=1,3
                    JPOC(I,K) = JPOCaccum(I,K)/ REAL(NDTSED)
                    JPON(I,K) = JPONaccum(I,K)/ REAL(NDTSED)
                    JPOP(I,K) = JPOPaccum(I,K)/ REAL(NDTSED)
                  ENDDO
                  JPOS(I) = JPOSaccum(I)/ REAL(NDTSED)
                ENDDO
            ENDIF                                                                               
 !       write(*,*)'About to call SED model'
            IF(STEADY_STATE_SED_TS)THEN    
               CALL SED_CALC(.TRUE.)  !quasi-steady state calculation 
            ELSE
               CALL SED_CALC(.FALSE.) !FALSE means use time-varying calculation
            ENDIF
            
            !set accumulative fluxes back to zero
            JPOCaccum = 0.0
            JPONaccum = 0.0
            JPOPaccum = 0.0
            JPOSaccum = 0.0
            
        ENDIF
    
        IF(MSR)THEN
                TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
                WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," SEDCALC"
                ENDIF   !SEDCALC 
        ELSE IF (BENTHIC_FLUXES) THEN
                CALL BEN_FLUX
        ENDIF

#endif

        IF(MSR)THEN
                TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
                NPROFILE=NPROFILE+1
        ENDIF
	
     
!******* Update time related variables

     NIT     = NIT+1
     ELTMS = ELTMS+DLT
     DLTAV   = ELTMS/NIT
     JDAY    = ELTMS/86400.
     ELTMJD  = JDAY-TMSTRT

     IF(JDAY.GE.TMEND) END_RUN=.TRUE. 
                    
#if defined(MULTIPROCESSOR)
     
     IF(PAR)THEN
        !Check and make sure all procs end the run at the same time
        IF(MSR)THEN  !Master receive END_RUN from each child
           DO I=2,NPROCS
              IRCVTAG=I+1000
              CALL MPI_RECV(END_RUN_CHILD,1,MPI_LOGICAL,I-1,IRCVTAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,IERR)
              IF(END_RUN)WRITE(*,*)'Child ID=',I,'CHILD END_RUN=',END_RUN_CHILD
           ENDDO
        ELSE  !Child procs send to master
           ISNDTAG=MYID+1000
           CALL MPI_SEND(END_RUN,1,MPI_LOGICAL,0,ISNDTAG,MPI_COMM_WORLD,IERR)
        ENDIF
     ENDIF
     
#endif

     IF(MSR)THEN
         IF(END_RUN)WRITE(*,*)'MYID =',MYID,'END_RUN=',END_RUN   !master print its own END_RUN
     ENDIF

     IF (JDAY >= DLTD(DLTDP+1)) THEN
       DLTDP = DLTDP+1
       IF (.NOT.AUTO_STEPPING) DLT = DLTVAL(DLTDP)
       MXDLT = DLTMAX(DLTDP)
       FNDLT = DLTFTN(DLTDP)
     ENDIF

! KURT GLAESEMANN, use local number of elements
!    DO I=1,NGL
     DO I=1,NTLOC
       DT1(I)=DT(NV(I,1))+DT(NV(I,2))+DT(NV(I,3))
       DT1(I)=DT1(I)/3.0
     ENDDO
 	 
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," KILLTREES1"
	 ENDIF	 



	 ! Add in Budget variables here and will output at the very end by changing
	 ! the output frequency to the length of the model run
  
#if defined (BUDGET)	  
	  !this is all in grams m^-3 timestep^-1, so need to multiply by the layer and nodal volume
	  !!get the volume
	  DO I = 1,MLOC
	
	  !Integrate sediment nutrogen fluxes in time and space
	  !in addition to the settling flux of POM from the overlying water column
	  
	   JNH4_out(I)=JNH4_out(I)+JNH4TM1S(I)*ART1(I)/86400*DLT/1E6   !kg N per nodal area
		 JNO3_out(I)=JNO3_out(I)+JNO3TM1S(I)*ART1(I)/86400*DLT/1E6		 !kg N per nodal area
		 JPOC_out(I,:)=JPOC_out(I,:)+JPOC(I,:)*ART1(I)/86400*DLT/1E6     !kg C per node
		 JPON_out(I,:)=JPON_out(I,:)+JPON(I,:)*ART1(I)/86400*DLT/1E6	 !kg N per node
		 !JPOP_out(I,:)=JPOP_out(I,:)+JPOP(I,:)*ART1(I)/86400*DLT/1E6	 !kg P per node
		 
	  DO K=1,KBM1
					 
	  MyVolume=(D(I)*DZ(K)*ART1(I))
								
	  ALGDON_OUT(I,K)=ALGDON_OUT(I,K)+ALGDON(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  ALGPON_OUT(I,K)=ALGPON_OUT(I,K)+ALGPON(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  ALGNH4_OUT(I,K)=ALGNH4_OUT(I,K)+ALGNH4(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  ALGNO3_OUT(I,K)=ALGNO3_OUT(I,K)+ALGNO3(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  DENNO3_OUT(I,K)=DENNO3_OUT(I,K)+DENNO3(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  NT_OUT(I,K)=NT_OUT(I,K)+NT(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  DENITN_OUT(I,K)=DENITN_OUT(I,K)+DENITN(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  HDRLPON_OUT(I,K)=HDRLPON_OUT(I,K)+HDRLPON(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  HDRRPON_OUT(I,K)=HDRRPON_OUT(I,K)+HDRRPON(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  COAGN_OUT(I,K)=COAGN_OUT(I,K)+COAGN(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  MNLDON1_OUT(I,K)=MNLDON1_OUT(I,K)+MNLDON1(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  MNLDON2_OUT(I,K)=MNLDON2_OUT(I,K)+MNLDON2(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  MNLDON3_OUT(I,K)=MNLDON3_OUT(I,K)+MNLDON3(I,K)/86400*DLT*MyVolume/1E3 !kg N per node per depth
	  	  
	  ALGDOC_OUT(I,K)=ALGDOC_OUT(I,K)+ALGDOC(I,K)/86400*DLT*MyVolume/1E3 !kg C per node per depth
	  ALGPOC_OUT(I,K)=ALGPOC_OUT(I,K)+ALGPOC(I,K)/86400*DLT*MyVolume/1E3 !kg C per node per depth
	  DENIT_OUT(I,K)=DENIT_OUT(I,K)+DENIT(I,K)/86400*DLT*MyVolume/1E3 !kg C per node per depth
	  HDRLPOC_OUT(I,K)=HDRLPOC_OUT(I,K)+HDRLPOC(I,K)/86400*DLT*MyVolume/1E3 !kg C per node per depth
	  HDRRPOC_OUT(I,K)=HDRRPOC_OUT(I,K)+HDRRPOC(I,K)/86400*DLT*MyVolume/1E3 !kg C per node per depth
	  COAGC_OUT(I,K)=COAGC_OUT(I,K)+COAGC(I,K)/86400*DLT*MyVolume/1E3 !kg C per node per depth
	  
	  MNLDOC1_OUT(I,K)=MNLDOC1_OUT(I,K)+MNLDOC1(I,K)/86400*DLT*MyVolume/1E3 !kg C per node per depth
	  MNLDOC2_OUT(I,K)=MNLDOC2_OUT(I,K)+MNLDOC2(I,K)/86400*DLT*MyVolume/1E3 !kg C per node per depth
	  MNLDOC3_OUT(I,K)=MNLDOC3_OUT(I,K)+MNLDOC3(I,K)/86400*DLT*MyVolume /1E3 !kg C per node per depth	 
	  
	  
	  PHOTOCHEM_DOC_ACCUM_OUT(I,K,:)=PHOTOCHEM_DOC_ACCUM_OUT(I,K,:)+PHOTOCHEM_DOC_ACCUM(I,K,:)/86400*DLT*MyVolume/1E3!kg C per node per depth	
	  	  
	  !! Below will be for calculating the change in concentration over time, with advection
	  !! C2 is the updated concentration after advection, the old one is from the previous time step
	  
	  !Delta DOC
	       DTWCDOC1_out(I,K)=DTWCDOC1_out(I,K)+(C2(I,K,9)*MyVolume-WC_CDOC1(I,K)*MyVolumeTM1(I,K))/1E3!kg C per node per depth	
	       DTWCDOC2_out(I,K)=DTWCDOC2_out(I,K)+(C2(I,K,33)*MyVolume-WC_CDOC2(I,K)*MyVolumeTM1(I,K))/1E3!kg C per node per depth	
		     DTWCDOC3_out(I,K)=DTWCDOC3_out(I,K)+(C2(I,K,35)*MyVolume-WC_CDOC3(I,K)*MyVolumeTM1(I,K))/1E3!kg C per node per depth	
		   
	       DTWNCDOC1_out(I,K)=DTWNCDOC1_out(I,K)+(C2(I,K,10)*MyVolume-WC_NCDOC1(I,K)*MyVolumeTM1(I,K))/1E3!kg C per node per depth	
	       DTWNCDOC2_out(I,K)=DTWNCDOC2_out(I,K)+(C2(I,K,34)*MyVolume-WC_NCDOC2(I,K)*MyVolumeTM1(I,K))/1E3!kg C per node per depth	
		     DTWNCDOC3_out(I,K)=DTWNCDOC3_out(I,K)+(C2(I,K,36)*MyVolume-WC_NCDOC3(I,K)*MyVolumeTM1(I,K))/1E3!kg C per node per depth	
	  !Delta DON
		     DTWCDON1_out(I,K)=DTWCDON1_out(I,K)+(C2(I,K,16)*MyVolume-WC_CDON1(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
	       DTWCDON2_out(I,K)=DTWCDON2_out(I,K)+(C2(I,K,37)*MyVolume-WC_CDON2(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
		     DTWCDON3_out(I,K)=DTWCDON3_out(I,K)+(C2(I,K,39)*MyVolume-WC_CDON3(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
		   
	       DTWNCDON1_out(I,K)=DTWNCDON1_out(I,K)+(C2(I,K,17)*MyVolume-WC_NCDON1(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
	       DTWNCDON2_out(I,K)=DTWNCDON2_out(I,K)+(C2(I,K,38)*MyVolume-WC_NCDON2(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
		     DTWNCDON3_out(I,K)=DTWNCDON3_out(I,K)+(C2(I,K,40)*MyVolume-WC_NCDON3(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
		   
	   !Delta POC
		     DTLPOC_out(I,K)=DTLPOC_out(I,K)+(C2(I,K,11)*MyVolume-LPOC(I,K)*MyVolumeTM1(I,K))/1E3!kg C per node per depth	
		     DTRPOC_out(I,K)=DTRPOC_out(I,K)+(C2(I,K,12)*MyVolume-RPOC(I,K)*MyVolumeTM1(I,K))/1E3!kg C per node per depth	

	   !Delta PON
		     DTLPON_out(I,K)=DTLPON_out(I,K)+(C2(I,K,18)*MyVolume-LPON(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
		     DTRPON_out(I,K)=DTRPON_out(I,K)+(C2(I,K,19)*MyVolume-RPON(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
		   
	   !Delta NH4
		     DTNH4_out(I,K)=DTNH4_out(I,K)+(C2(I,K,13)*MyVolume-NH4(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
		   
	   !Delta NO3
		     DTNO3_out(I,K)=DTNO3_out(I,K)+(C2(I,K,14)*MyVolume-NO3(I,K)*MyVolumeTM1(I,K))/1E3!kg N per node per depth	
		   
		     MyVolumeTM1(I,K)=MyVolume
	     enddo
	  enddo	
	  
#endif	  
	    
	  
!RGL commented below
! T(:,:)    = C2(:,:,1)
!   SALT(:,:) = C2(:,:,2)
   SSI(:,:)  = C2(:,:,3)
   B1(:,:)   = C2(:,:,4)
   B2(:,:)   = C2(:,:,5)
   B3(:,:)   = C2(:,:,6)
   SZ(:,:)   = C2(:,:,7)
   LZ(:,:)   = C2(:,:,8)
   
   WC_CDOC1(:,:) = C2(:,:,9)
   WC_NCDOC1(:,:) = C2(:,:,10)  


   !LDOC(:,:) = C2(:,:,9)
   !RDOC(:,:) = C2(:,:,10)
   LPOC(:,:) = C2(:,:,11)
   RPOC(:,:) = C2(:,:,12)
   NH4(:,:)  = C2(:,:,13)
   NO3(:,:)  = C2(:,:,14)
   UREA(:,:) = C2(:,:,15)
   
   WC_CDON1(:,:) = C2(:,:,16)
   WC_NCDON1(:,:) = C2(:,:,17)
   !LDON(:,:) = C2(:,:,16)
   !RDON(:,:) = C2(:,:,17)
   LPON(:,:) = C2(:,:,18)
   RPON(:,:) = C2(:,:,19)
   PO4(:,:)  = C2(:,:,20)
   
   WC_CDOP1(:,:) = C2(:,:,21)
   WC_NCDOP1(:,:) = C2(:,:,22)
   !LDOP(:,:) = C2(:,:,21)
   !RDOP(:,:) = C2(:,:,22)
   LPOP(:,:) = C2(:,:,23)
   RPOP(:,:) = C2(:,:,24)
   PIP(:,:)  = C2(:,:,25)
   COD(:,:)  = C2(:,:,26)
   DOXG(:,:)   = C2(:,:,27)
   SIUPB(:,:)   = C2(:,:,28) 
   SIAT(:,:)   = C2(:,:,29)
   PIB1(:,:) = C2(:,:,30)
   PIB2(:,:) = C2(:,:,31)
   PIB3(:,:) = C2(:,:,32)
   
   WC_CDOC2(:,:) = C2(:,:,33)
   WC_NCDOC2(:,:) = C2(:,:,34)

   WC_CDOC3(:,:) = C2(:,:,35)
   WC_NCDOC3(:,:) = C2(:,:,36)


   WC_CDON2(:,:) = C2(:,:,37)
   WC_NCDON2(:,:) = C2(:,:,38)
   WC_CDON3(:,:) = C2(:,:,39)
   WC_NCDON3(:,:) = C2(:,:,40)
   WC_CDOP2(:,:) = C2(:,:,41)
   WC_NCDOP2(:,:) = C2(:,:,42)
   WC_CDOP3(:,:) = C2(:,:,43)
   WC_NCDOP3(:,:) = C2(:,:,44)
     
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
   
#if defined(XFLUX_FLAG)   
! write(*,*)' Calcuating the XFLUX out'
!Added by B Clark to time integrate xflux for each WQ var
   EDGELINE_FLUX_IN=EDGELINE_FLUX_in+EDGELINE_FLUX
! write(*,*)'Done Calcuating the XFLUX out'
#endif

	 !set up loop to output the netcdf at the correct interval
	 IF(CDF_OUT) THEN
	     IF(JDAY >= CDF_START .AND. JDAY <= CDF_END)THEN
	   CDF_COUNTER=CDF_COUNTER+DLT
	      IF(CDF_COUNTER .EQ. CDF_INT) THEN
		  !write(*,*)'NPP = ',total_netPP(1)
             CALL OUT_NETCDF
            !  write(*,*)'Wrote NETCDF'
#if defined(XFLUX_FLAG)             
             CALL WRITE_EDGELINE_FLUX
             EDGELINE_FLUX_IN=0.
#endif             
	!		 write(*,*)'Wrote NETCDF'
			 CDF_COUNTER = 0.
		  ENDIF
	    ENDIF
	  ENDIF	  
	  
    CALL  WQM_OUTPUT(NTHIS_OLD, NTSTN_OLD, NTHIS,NTSTN,CTR_HIS,CTR_STN)
	
	!IF(SED_DOM) THEN! B Clark output sediment DOM variables
	
	!	CALL SED_DOM_OUTPUT
	!	
	!END IF
	! 
	! 
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," OUTPUT"
	 ENDIF	 
	
!************************************************************************
!**                           Output Results                           **
!************************************************************************

!******* Snapshots
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF



	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," KILLTREES2"
	 ENDIF	 
	   
!       
!******* Timestep, volume balance, and mass balance diagnostics
!

! KURT GLAESEMANN - add MSR to printout to DIA
	 
 
	  IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	
       IF (DIAGNOSTICS  .AND. MSR) THEN
         IF (JDAY >= NXDIA .OR. JDAY >= DIAD(DIADP+1)) THEN
           IF (JDAY >= DIAD(DIADP+1)) THEN
             DIADP = DIADP+1
             NXDIA = DIAD(DIADP)
           ENDIF
           NXDIA = NXDIA+DIAF(DIADP)
           WRITE(DIA,4000) INT(JDAY),(JDAY-INT(JDAY))*24.,NIT
           IF (AUTO_STEPPING) THEN
             WRITE(DIA,4010) INT(DLT),INT(DLTAV),COURMX,COURFS,COURBS
             WRITE(DIA,4012) INT(DLT),INT(DLTAV),DIFFMX,DIFFFS,DIFFBS
           ENDIF
           IF (VOLUME_BALANCE .AND. NEW_VOLUMES) THEN

           ENDIF
           IF (MASS_BALANCE) THEN
             WRITE (DIA,4050) 'Constituent total mass',(CNAME(AC(JCON)),   &
                                CMASS(AC(JCON)),JCON=1,NAC)
             WRITE (DIA,4060) ERRWCN,ERRWCP,ERRWCC,ERRSEDN,ERRSEDP,      &
                              ERRSEDC
           ENDIF
         ENDIF
       ENDIF !END DIAGNOSTICS

!       
!******* Restarts
!

! KURT GLAESEMANN - restarts
       IF (RESTART_OUT) THEN
            IF (JDAY >= RSOD(RSODP) .OR. END_RUN) THEN
!                IF(.NOT. ALLOCATED(CTEMP_GL)) THEN
!                    ALLOCATE(CTEMP_GL(1:MGL))
!                    IF (SEDIMENT_CALC) THEN
!                        ALLOCATE(CPOP_GL(1:MGL,1:3))
!                        ALLOCATE(CPON_GL(1:MGL,1:3))
!                        ALLOCATE(CPOC_GL(1:MGL,1:3))
!                        ALLOCATE(CPOS_GL(1:MGL))
!                        ALLOCATE(PO4T2TM1S_GL(1:MGL))
!                        ALLOCATE(NH4T2TM1S_GL(1:MGL))
!                        ALLOCATE(NO3T2TM1S_GL(1:MGL))
!                        ALLOCATE(HST2TM1S_GL(1:MGL))
!                        ALLOCATE(CH4T2TM1S_GL(1:MGL))
!                        ALLOCATE(CH41TM1S_GL(1:MGL))
!                        ALLOCATE(SO4T2TM1S_GL(1:MGL))
!                        ALLOCATE(SIT2TM1S_GL(1:MGL))
!                        ALLOCATE(BENSTRTM1S_GL(1:MGL))
!                        ALLOCATE(BBM_GL(1:MGL))
!                    ENDIF

!!WLong moved this to mod_df.F           
!!                   IF(DFEEDER) THEN
!!                         ALLOCATE(DFEEDM1S_GL(1:MGL))
!!                   ENDIF

!!WLong moved this to SF_MOD.F           
!!                   IF(SFEEDER) THEN
!!                        ALLOCATE(SFEED_GL(1:MGL,1:NSSFP))
!!                   ENDIF
!           
!                ENDIF
# if defined (MULTIPROCESSOR)
                if(PAR) then
                
                    DO J=1,NCP
                        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C1(0,1,J),C2_GL(1,1,J))
                    ENDDO
                    IF (SEDIMENT_CALC) THEN
					
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CTEMP    ,CTEMP_GL)
                        call gather(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPOP     ,CPOP_GL)
                        call gather(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPON     ,CPON_GL)
                        call gather(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPOC     ,CPOC_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CPOS     ,CPOS_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,PO4T2TM1S,PO4T2TM1S_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NH4T2TM1S,NH4T2TM1S_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NO3T2TM1S,NO3T2TM1S_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,HST2TM1S ,HST2TM1S_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CH4T2TM1S,CH4T2TM1S_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CH41TM1S ,CH41TM1S_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SO4T2TM1S,SO4T2TM1S_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SIT2TM1S ,SIT2TM1S_GL)
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,BENSTRTM1S ,BENSTRTM1S_GL)
                        !call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SODTM1S, SODTM1S_GL)
                        IF(BALGAE_CALC)THEN
                            call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,BBM      ,BBM_GL)
                        ENDIF
                    ENDIF
                    IF(DFEEDER) THEN
                        call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,DFEEDM1S ,DFEEDM1S_GL)
                    ENDIF
                    IF(SFEEDER) THEN
                        call gather(1,MTLOC,MLOC,MGL,NSSFP,MYID,NPROCS,NMAP,SFEED,SFEED_GL)
                    ENDIF
                    
                else
#endif
                    DO J=1,NCP
                        C2_GL(1:MGL,1:KBM1,J)  = C2(1:MGL,1:KBM1,J)
                    ENDDO
                    IF (SEDIMENT_CALC) THEN
                        CTEMP_GL(1:MGL)     = CTEMP(1:MGL)
                        CPOP_GL(1:MGL,1:3)  = CPOP(1:MGL,1:3)
                        CPON_GL(1:MGL,1:3)  = CPON(1:MGL,1:3)
                        CPOC_GL(1:MGL,1:3)  = CPOC(1:MGL,1:3)
                        CPOS_GL(1:MGL)      = CPOS(1:MGL)
                        PO4T2TM1S_GL(1:MGL) = PO4T2TM1S(1:MGL)
                        NH4T2TM1S_GL(1:MGL) = NH4T2TM1S(1:MGL)
                        NO3T2TM1S_GL(1:MGL) = NO3T2TM1S(1:MGL)
                        HST2TM1S_GL(1:MGL)  = HST2TM1S(1:MGL)
                        CH4T2TM1S_GL(1:MGL) = CH4T2TM1S(1:MGL)
                        CH41TM1S_GL(1:MGL)  = CH41TM1S(1:MGL)
                        SO4T2TM1S_GL(1:MGL) = SO4T2TM1S(1:MGL)
                        SIT2TM1S_GL(1:MGL)  = SIT2TM1S(1:MGL)
                        BENSTRTM1S_GL(1:MGL)  = BENSTRTM1S(1:MGL)
                        !SODTM1S_GL(1:MGL)     = SODTM1S(1:MGL)
                        !changed below for consistency RGL/Kurt
                        IF(BALGAE_CALC)THEN
                            BBM_GL(1:MGL) = BBM(1:MGL)
                        ENDIF
                            
                    ENDIF
                    IF(DFEEDER) THEN
                        DFEEDM1S_GL(1:MGL) = DFEEDM1S(1:MGL)
                    ENDIF
                    IF(SFEEDER) THEN
                        SFEED_GL(1:MGL,1:NSPECIES) = SFEED(1:MGL,1:NSPECIES)
                    ENDIF
# if defined (MULTIPROCESSOR)
                endif
#endif
                RSODP  = RSODP+1
                IF(MSR) THEN
                    LENGTH = INDEX(RSOFN,'.')
                    IJDAY  = AINT(JDAY)
                    IF (IJDAY < 10) THEN
                        WRITE (EXT1,'(I1)') IJDAY
                        FILENAME = RSOFN(1:LENGTH)//EXT1
                    ELSE IF (IJDAY < 100) THEN
                        WRITE (EXT2,'(I2)') IJDAY
                        FILENAME = RSOFN(1:LENGTH)//EXT2
                    ELSE IF (IJDAY < 1000) THEN
                        WRITE (EXT3,'(I3)') IJDAY
                        FILENAME = RSOFN(1:LENGTH)//EXT3
                    ELSE
                        WRITE (EXT4,'(I4)') IJDAY
                        FILENAME = RSOFN(1:LENGTH)//EXT4
                    ENDIF

                    OPEN  (RSO,FILE=FILENAME, FORM='UNFORMATTED')
                    WRITE (RSO) TITLE
                    WRITE (RSO) (((C2_GL(I,K,J),I=1,MGL ),K=1,KBM1),J=1,NCP)
        
                    IF (SEDIMENT_CALC)THEN
                    
                            WRITE (RSO) (CTEMP_GL(B),B=1,MGL ),                                &
                                ((CPOP_GL(B,J),B=1,MGL ),J=1,3),                              &
                                ((CPON_GL(B,J),B=1,MGL ),J=1,3),                              &
                                ((CPOC_GL(B,J),B=1,MGL ),J=1,3),                              &
                                (CPOS_GL(B),B=1,MGL ),                                        &
                                (PO4T2TM1S_GL(B),B=1,MGL ), (NH4T2TM1S_GL(B),B=1,MGL ),       &
                                (NO3T2TM1S_GL(B),B=1,MGL ), (HST2TM1S_GL(B),B=1,MGL ),        &
                                (CH4T2TM1S_GL(B),B=1,MGL ), (CH41TM1S_GL(B),B=1,MGL ),        &
                                (SO4T2TM1S_GL(B),B=1,MGL ), (SIT2TM1S_GL(B),B=1,MGL ),        &
                                (BENSTRTM1S_GL(B),B=1,MGL), (SODTM1S_GL(B),B=1,MGL)
                        
                                IF(BALGAE_CALC)THEN
                                    WRITE(RSO)(BBM_GL(B),B=1,MGL)
                                ENDIF
                    ENDIF
                    
                    IF (DFEEDER) WRITE (RSO) (DFEEDM1S_GL(B),B=1,MGL )
                    IF (SFEEDER) WRITE (RSO) ((SFEED_GL(B,K),B=1,MGL ),K=1,NSPECIES)

                    CLOSE (RSO)
          ENDIF
       ENDIF
     ENDIF  !END RESTART_OUT
	 

	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=TOTAL_TIME_PROFILE+DTIME(ELAPSED_PROFILE)
		WRITE(*,'(A24,I5.0,A3,EN15.3,A20)')"Total time used in part ",NPROFILE," = ",TOTAL_TIME_PROFILE," RESTART"
	 ENDIF	 	 

!Wen Long Debugging
!     WRITE(*,*)'MYID=',MYID,'END_RUN=','END_RUN'

   ENDDO   !END DO WHILE(.NOT.END_RUN)


!************************************************************************
!**                          End Simulation                            **
!************************************************************************

!WLong, close station and history outputs if open
   IF(MSR)THEN
      IF(OPEN_HIS) THEN
         SELECT CASE (HISFN_SPLIT_BYLEVEL)
            CASE (.TRUE.)
              DO K=1,KBM1
                 CLOSE(UNIT_HIS+K)
              ENDDO
            CASE (.FALSE.)
              CLOSE(UNIT_HIS)
         ENDSELECT
      ENDIF
      IF(OPEN_STN) CLOSE(UNIT_STN)
   ENDIF

!Wen Long debug
# if defined (MULTIPROCESSOR)
   IF(PAR)THEN
!      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      !WRITE(*,*)'MYID =', MYID, 'Debug 1'
   ENDIF
# endif   


!***** Write results in binary form for use as initial conditions

   IF (ICOND_OUT) THEN

!***** Integrate sediments POM to steady state                               !MNOEL 2-5-93

     IF (STEADY_STATE_SED_G3) CALL SED_DIAGENESIS_G3  

   ENDIF


!Wen Long debug
# if defined (MULTIPROCESSOR)
   IF(PAR)THEN
!      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      !WRITE(*,*)'MYID =', MYID, 'Debug 2'
   ENDIF
# endif


!***** Snapshot FORMAT statements

 3000 FORMAT('1')
 3010 FORMAT(1X,A72)
 3020 FORMAT(//                                                        &
             1X,'Time related variables'/'+',22('_')//                 &
             3X,'Julian day',T25,'=',F8.2/                             &
             3X,'Elapsed time',T25,'=',I8,' days ',F6.2,' hours'/      &
             3X,'Number of iterations',T25,'=',I8/                     &
             3X,'Timestep',T25,'=',F8.0,' sec'/                        &
             3X,'Average time step',T25,'=',F8.0)
 3025 FORMAT(//                                                          &
             1X,'After steady-state sediment computation'/'+',22('_')//  &
             3X,'Julian day',T25,'=',F8.2/                               &
             3X,'Elapsed time',T25,'=',I8,' days ',F6.2,' hours'/        &
             3X,'Number of iterations',T25,'=',I8/                       &
             3X,'Timestep',T25,'=',F8.0,' sec'/                          &
             3X,'Average time step',T25,'=',F8.0)  
 3030 FORMAT(/1X,A24/10(/10F8.4))
 3032 FORMAT(/1X,A20/10(/10F8.3))
 3034 FORMAT(/1X,A20/10(/10(1PE8.2E1)))
 3035 FORMAT(/1X,'G',I1,1X,A20/10(/10(1PE8.2E1)))
 3036 FORMAT(/1X,'N',I1,1X,A20/10(/10(1PE8.2E1)))

!***** Diagnostic FORMAT statements

 4000 FORMAT(/                                                        &
             1X,79('*')//                                             &
             1X,'Date         =',I8,'  days ',F5.2,' hours'/          &
             1X,'# iterations =',I8)
 4010 FORMAT(1X,'Timestep',T15,'=',I8,2X,'Average timestep',T43,      &
               '=',I8/                                                &
             1X,'Courant #',T15,'=',F8.3/                             &
             1X,'at face',T15,'=',I8,2X,'and box',T43,'=',I8/)
 4012 FORMAT(1X,'Timestep',T15,'=',I8,2X,'Average timestep',T43,      &
               '=',I8/                                                &
             1X,'Diffusion #',T15,'=',F8.3/                           &
             1X,'at face',T15,'=',I8,2X,'and box',T43,'=',I8/)
 4020 FORMAT(/                                                        &
             1X,A14,' performed at ',I6,' days and ',F5.2,' hours')
! 4030 FORMAT(5X,'V2(',I4,') =',1PE15.6E1,' m*3',3X,'HMV(',I4,') =',   &
!               1PE15.6E1,' m*3',3X,'Difference =',1PE15.6E1,' m*3')
 4040 FORMAT(/                                                          &
             3X,'Water quality model volume      =',1PE20.12E2,' m*3'/  &
             3X,'Difference in total grid volume =',1PE20.3E2,' %'/     &
             3X,'Total difference in volume      =',1PE20.3E2,' m*3')
 4050 FORMAT(/                                                          &
             1X,A22//                                                   &
             (:3X,A24,T31,'=',1PE16.8E2,' kg'))
 4060 FORMAT(/                                                          &
             1X,'Water column mass balance error'//                     &
             3X,'Nitrogen   = ',1PE10.3E2,' %'/                         &
             3X,'Phosphorus = ',1PE10.3E2,' %'/                         &
             3X,'Carbon     = ',1PE10.3E2,' %'//                        &
             1X,'Sediment mass balance error'//                         &
             3X,'Nitrogen   = ',1PE10.3E2,' %'/                         &
             3X,'Phosphorus = ',1PE10.3E2,' %'/                         &
             3X,'Carbon     = ',1PE10.3E2,' %'/)
 4070 FORMAT(9F8.2)
 4080 FORMAT(8X,A24)

!***** Mass balance FORMAT's

 5000 FORMAT(10(1X,1PE10.2E2))
!KURT below
 5001 CONTINUE

 
!!---Wen Long deallocate buffers
#if defined (MULTIPROCESSOR)
!   IF(ALLOCATED(SNDBUF))DEALLOCATE(SNDBUF)
!   IF(ALLOCATED(RCVBUF))DEALLOCATE(RCVBUF)
#endif
   
     CALL  WQM_DEALLOC !WLong repaced the above section by WQM_DEALLOC

# if defined(XFLUX_FLAG)     
     CALL DEALLOC_EDGELINE_FLUX
# endif
    
     CALL ZOOP_DEALLOC !WLong replaced the above by ZOOP_DEALLOC()
     IF(LIGHT_EXTINCTION)CALL OWQ_DEALLOC    !WLong replace the above with OWQ_DEALLOC
 
!====The following are now moved to SAV_DEALLOC() function in wqm_sav.F
     IF(SAV_CALC) THEN
        CALL SAV_DEALLOC
     ENDIF
     !Wen Long: deallocate SAV fluxes to water column if SAV_LOADS is on 
     IF(SAV_LOADS) THEN
        CALL SAV_LOADS_DEALLOC
     ENDIF
     
!=========Above are moved to SAV_DEALLOC() in wqm_sav.F ====
     
     IF(SFEEDER)THEN
        CALL SF_DEALLOC
     ENDIF
     
     CALL WQMINIT_DEALLOC !WLong  replaced the above by WQMINIT_DEALLOC()
     

     CALL HYDRO_GEOM_DEALLOC !WLong replaced the above by HYDRO_GEOM_DEALLOC()

       CALL OBCS_DEALLOC
  
    !WLong moved the above to HYDRO_DEALLOC() in mod_hydrovars.F
     CALL HYDRO_DEALLOC
     
    !WLong move the above to BCS_DEALLOC in mod_bcs.F
     CALL BCS_DEALLOC
     !WLong moved the above to BCMAP_DEALLOC() in mod_bcmap.F
     CALL BCMAP_DEALLOC
     
     CALL TGE_DEALLOC !replaced the above by TGE_DEALLOC
          
     !WLong moved to mod_tge.F and remove NBE usage in wqm_main.F here
     !should move to mod_hydrovars.F
     

#if defined (MULTIPROCESSOR)
     if(PAR) THEN
     
        DO I = 1, NPROCS
            IF(ALLOCATED(EMAP(I)%LOC_2_GL)) DEALLOCATE(EMAP(I)%LOC_2_GL)
            IF(ALLOCATED(NMAP(I)%LOC_2_GL)) DEALLOCATE(NMAP(I)%LOC_2_GL)
            IF(ALLOCATED(EC(I)%RCVP)) DEALLOCATE(EC(I)%RCVP)
            IF(ALLOCATED(EC(I)%SNDP)) DEALLOCATE(EC(I)%SNDP)
            IF(ALLOCATED(BNC(I)%RCVP)) DEALLOCATE(BNC(I)%RCVP)
            IF(ALLOCATED(BNC(I)%SNDP)) DEALLOCATE(BNC(I)%SNDP)
            IF(ALLOCATED(BNC(I)%MLTP)) DEALLOCATE(BNC(I)%MLTP)
            IF(ALLOCATED(NC(I)%RCVP)) DEALLOCATE(NC(I)%RCVP)
            IF(ALLOCATED(NC(I)%SNDP)) DEALLOCATE(NC(I)%SNDP)
        ENDDO
        
     ENDIF
     
     CALL PAR_DEALLOC
#endif
     
    IF(SEDIMENT_CALC)THEN 

       !Wen Long deallocating sediment diagenesis related variables 
        CALL SED_DEALLOC
    
       !Wen Long deallocating benthic algae related variables
        IF(BALGAE_CALC)    CALL BA_DEALLOC
       
        IF(SAV_CALC)     CALL SED_SAV_EXCHANGE_DEALLOC
        
        !WLong deallocating deposition feeder related variables
        IF(DFEEDER)THEN
            CALL DF_DEALLOC
            CALL SED_DF_EXCHANGE_DEALLOC
        ENDIF
#if defined(SEDIMENT_DOM)		
		IF (SED_DOM)THEN
		
			CALL SED_DOM_DEALLOC
			CALL SED_DOM_SHARE_DEALLOCATE
			
		ENDIF
#endif
        IF(SEDIMENT_CALC)THEN
            CALL SED_SF_EXCHANGE_DEALLOC
        ENDIF

    ENDIF       
    
    CALL WC_DOM_DEALLOCATE
	
!Wen Long, report model time 
    IF(MSR)THEN
        TOTAL_TIME=TOTAL_TIME+DTIME(ELAPSED)    !CALL DTIME() and return user time in ELAPSED(1), system time in ELAPSED(2)
        WTIME_MOD =WTIME_MOD+ELAPSED(1)
        WTIME_SYS =WTIME_SYS+ELAPSED(2)
        WRITE(*,'(A22,EN15.3,A8)')"Time used by model   =",WTIME_MOD, " seconds"
        WRITE(*,'(A22,EN15.3,A8)')"Time used by system  =",WTIME_SYS, " seconds"
        WRITE(*,'(A22,EN15.3,A8)')"Total time used      =",TOTAL_TIME," seconds"
    ENDIF

   !CALL MPI_FINALIZE(IERR)
   STOP     !Wen Long and Laura found this does not stop the run properly
   

   END PROGRAM PARWQM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   SUBROUTINE EXCHANGE_WQM(CCC)
!==============================================================================!
!     PERFORM DATA EXCHANGE FOR WATER QUALITY VARIABLES                        |
!==============================================================================!
   USE MOD_PREC,  ONLY:  SP
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
            NBN ,       &   !              !!NUMBER OF BOUNDARY NODES
            !MX_MLT,     &   !              !!MAX MULTIPLICITY OF BOUNDARY NODES
            !BN_LST,     &    !               !!GLOBAL IDENTITY OF BOUNDARY NODES
            BN_LOC,     &    !            !!LOCAL IDENTITY OF BOUNDARY NODES
            BN_MLT,     &    !            !!MULTIPLICITY OF BOUNDARY NODES
            !BN_NEY,     &    !            !!NODE OWNER LIST
            !NDE_ID,     &    !            !! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
            !NHN,          &    !            !!NUMBER OF HALO NODES
            !HN_LST!,     &    !            !!LIST OF HALO NODES 
            !HN_OWN,     &    !            !!PRIMARY OWNER OF HALO NODES
            !EC,        &    !
            NC,            &    !
            BNC,         &    !
            !EMAP,        &    !
            !NMAP,         &    !
            !PNE,        &    !(:)        !!NUMBER OF ELEMENTS IN EACH PROC
            !PNN,        &    !(:)        !!NUMBER OF NODES IN EACH PROC
            !PNHE,        &    !(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
            !PNBN,        &    !(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
            !PMBM,        &    !(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
            !PNHN,        &    !(:)       !!NUMBER OF HALO NODES IN EACH PROC
            !FILELOCK,  & !
            !MYFILELOCK, & !
            !FILELOCKED, & !
            !NODE_MATCH, & !
             NODE_MATCH_ALL, & !
            !EXCHANGE, & !
             EXCHANGE_ALL !, & !
            !GATHER,  & !
            !IGATHER, & !
            !INIT_MPI_ENV,&!             
            !SORT,     & !
            !PREAD,     & !
            !PWRITE, & !
            !PPRINT, & !
            !LOCK,   & !
            !UNLOCK, & !
            !GETLOC
#  endif     
             
    USE MOD_SIZES, ONLY : NCP
    
    USE MOD_LIMS, ONLY :  &
          NLOC,         &       !!NUMBER OF ELEMENTS 
         ! MLOC,         &       !!NUMBER OF NODES
          !NISBCE_1,     &      !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 1
          !NISBCE_2,     &      !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 2
          !NISBCE_3,     &      !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 3
          !KB,           &      !!NUMBER OF SIGMA LEVELS
          KBM1,          &      !!NUMBER OF SIGMA LEVELS-1
          !KBM2 ,        &      !!NUMBER OF SIGMA LEVELS-2
          MYID ,         &      !!UNIQUE PROCESSOR ID (1 => NPROCS)
          NPROCS,        &      !!NUMBER OF PROCESSORS
          !NE   ,        &      !!NUMBER OF UNIQUE EDGES
          !NCV  ,        &      !!NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)
          !IINT  ,       &      !!TYKIM added for nudging
          !NCV_I ,       &      !!NUMBER OF INTERNAL CONTROL VOLUMES (LOCAL ONLY)
          !NTLOC ,       &      !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
          MTLOC! ,       &      !!TOTAL OF LOCAL INTERNAL + HALO NODES
          !NCT   ,       &
          !MX_NBR_ELEM , &      !!MAX NUMBER OF ELEMENTS SURROUNDING A NODE
          !NUMQBC_GL ,   &
          !NUMPNT_GL ,   &
          !NUMQBC,       &
          !NUMPNT,       &
          !NstationMax         !Maximum number of stations 

       !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY :         &
    !        SERIAL          &           !!TRUE IF SINGLE PROCESSOR
            MSR            &           !!TRUE IF MASTER PROCESSOR (MYID==1)
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

   IMPLICIT NONE
# if defined (MULTIPROCESSOR)
!   include "mpif.h"
# endif   
   
   REAL(SP), INTENT(INOUT) :: CCC(0:MTLOC,KBM1,NCP)
   INTEGER ::     NJ,            &!
                JCON,        &!Index for constituents
                NJ1,NJ2,NJ3,EXTRA

                
                
        !IF(MSR)write(*,*)'LBnote 1.1: size(C2,1)=',SIZE(CCC,1)
        !IF(MSR)write(*,*)'LBnote 1.1: size(C2,2)=',SIZE(CCC,2)
        !IF(MSR)write(*,*)'LBnote 1.1: size(C2,3)=',SIZE(CCC,3)
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
   
   

! KURT GLAESEMANN 16 SEPT 2009 ONLY LOOP OVER USED VARIABLES
!    DO NJ = 1,NCP
!    DO JCON=1,NAC
!      NJ = AC(JCON)
!    KURT GLAESEMANN - fix dimensions
!      CALL EXCHANGE(NC,MTLOC,KB,MYID,NPROCS,C2F(:,:,NJ))
!      CALL EXCHANGE(NC,MTLOC,KBM1,MYID,NPROCS,C2F(:,:,NJ))
!    ENDDO
! KURT GLAESEMANN 17 SEPT 2009 DO MULTIPLE VARIABLES AT ONCE
! ask Kurt if this just saves time -- why is he using 3 as extra
!   
#if(0)
     EXTRA = MOD(NAC,3)
!     if(MSR)write(*,*)'LBnote: extra=',EXTRA
     IF (EXTRA .eq. 1) THEN
       NJ1= AC(1)
       CALL EXCHANGE(NC,MTLOC,KBM1,MYID,NPROCS,CCC(:,:,NJ1))
     ELSE IF (EXTRA .eq. 2) THEN
       NJ1= AC(1)
       NJ2= AC(2)
       CALL EXCHANGE(NC,MTLOC,KBM1,MYID,NPROCS,CCC(:,:,NJ1),CCC(:,:,NJ2))
     ENDIF
        !IF(MSR)write(*,*)'LBnote 1.2: size(C2,1)=',SIZE(CCC,1)
        !IF(MSR)write(*,*)'LBnote 1.2: size(C2,2)=',SIZE(CCC,2)
        !IF(MSR)write(*,*)'LBnote 1.2: size(C2,3)=',SIZE(CCC,3)
     DO JCON=EXTRA+1,NAC-2,3
       NJ1= AC(JCON+0)
       NJ2= AC(JCON+1)
       NJ3= AC(JCON+2)
       CALL EXCHANGE(NC,MTLOC,KBM1,MYID,NPROCS,CCC(:,:,NJ1),CCC(:,:,NJ2),CCC(:,:,NJ3))
     ENDDO
        !IF(MSR)write(*,*)'LBnote 1.3: size(C2,1)=',SIZE(CCC,1)
        !IF(MSR)write(*,*)'LBnote 1.3: size(C2,2)=',SIZE(CCC,2)
        !IF(MSR)write(*,*)'LBnote 1.3: size(C2,3)=',SIZE(CCC,3)
#else
!  KURT GLAESEMANN 22 SEPT 2009 - DO ALL AT ONCE
     CALL EXCHANGE_ALL(NC,MTLOC,KBM1,NCP,MYID,NPROCS,CCC)
        !IF(MSR)write(*,*)'LBnote 1.4: size(C2,1)=',SIZE(CCC,1)
        !IF(MSR)write(*,*)'LBnote 1.4: size(C2,2)=',SIZE(CCC,2)
        !IF(MSR)write(*,*)'LBnote 1.4: size(C2,3)=',SIZE(CCC,3)
#endif
     CALL NODE_MATCH_ALL(1,NBN,BN_MLT,BN_LOC,BNC,MTLOC,KBM1,NCP,MYID,NPROCS,CCC)
   ENDIF
   
#  endif
   
        !IF(MSR)write(*,*)'LBnote 1.5: size(C2,1)=',SIZE(CCC,1)
        !IF(MSR)write(*,*)'LBnote 1.5: size(C2,2)=',SIZE(CCC,2)
        !IF(MSR)write(*,*)'LBnote 1.5: size(C2,3)=',SIZE(CCC,3)
   
   RETURN
   END SUBROUTINE EXCHANGE_WQM
!==============================================================================!


   SUBROUTINE WQM_OUTPUT(NTHIS_OLD, NTSTN_OLD, NTHIS,NTSTN,CTR_HIS,CTR_STN)

!subroutine to write history and station outputs
   
   
    USE MOD_PREC, ONLY: SP
    USE MOD_SIZES, ONLY: NCP, NSSFP
    USE MOD_LIMS, ONLY: MYID, NPROCS, IINT, NTLOC, MLOC
    
    USE MOD_TGE, ONLY: NV, TGE_DEALLOC
    
    USE MOD_HYDROVARS, ONLY:  &
           !GRAV    &        !
          PI    &        !
		  
        !,PI2    &        !
        !,ZERO    &        !
        !,ONE_THIRD    &    !
!         ,NVG    &        !
!         ,XG    &        !GLOBAL X-COORD AT NODE 
!         ,YG    &        !GLOBAL X-COORD AT NODE 
!         ,HG    &        !GLOBAL DEPTH AT NODE 
!         ,XCG    &        !GLOBAL X-COORD AT FACE CENTER 
!         ,YCG    &        !GLOBAL X-COORD AT FACE CENTER 
        !,VXMIN    &        !
        !,VYMIN    &        !
        !,VXMAX    &        !
        !,VYMAX    &        !
!         ,XC    &        !X-COORD AT FACE CENTER 
!         ,YC    &        !Y-COORD AT FACE CENTER
!         ,VX    &        !X-COORD AT GRID POINT
!         ,VY    &        !Y-COORD AT GRID POINT
!         ,ART    &        !AREA OF ELEMENT
         ,ART1   &!,ART1_GL    &        !AREA OF NODE-BASE CONTROl VOLUME
!         ,ART2    &        !AREA OF ELEMENTS AROUND NODE
!         ,NV    &        !NODE NUMBERING FOR ELEMENTS
!         ,NBE    &        !INDICES OF ELMNT NEIGHBORS
!         ,NTVE    &        !
!         ,NTSN    &        !
!         ,ISONB    &        !NODE MARKER = 0,1,2 
!         ,ISBC    &        !
!         ,ISBCE    &        !
!         ,IEC    &        !
!         ,IENODE &        !
!         ,NBSN    &        !
!         ,NIEC    &        !
!         ,NTRG    &        !
!         ,NBVE    &        !
!         ,NBVT    &        !
!         ,LISBCE_1    &    !LIST OF ELEMENTS WITH ISBCE=1
!         ,LISBCE_2    &    !LIST OF ELEMENTS WITH ISBCE=2
!         ,LISBCE_3    &    !LIST OF ELEMENTS WITH ISBCE=3
!         ,DLTXC    &        !
!         ,DLTYC    &        !
!         ,DLTXYC    &    !
!         ,DLTXE    &        !
!         ,DLTYE    &        !
!         ,DLTXYE    &    !
!         ,SITAC    &        !
!         ,SITAE    &        !
!         ,XIJC    &        !
!         ,YIJC    &        !
!         ,XIJE    &        !
!         ,YIJE    &        !
!         ,EPOR    &        !ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
        !,IBCGEO    &    !LOCAL GEOSTROPHIC FRICTION CORRECTION NODES

         ,Z    &            !SIGMA COORDINATE VALUE 
         ,ZZ    &        !INTRA LEVEL SIGMA VALUE
         ,DZ    &        !DELTA-SIGMA VALUE
         ,DZZ    &        !DELTA OF INTRA LEVEL SIGMA 
!         ,H1    &        !BATHYMETRIC DEPTH 
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
         ,UARD_OBCNT    &!tykim previous time step
         ,XFLUX_OBCT    &!tykim previous time step
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
         ,THOUR    &        !
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
         ,t_his_start    &        !
         ,t_his_end    &            !
         ,t_his_dlt    &            !starting time, ending time, and interval of history outputs (days)
         ,Nstation    &            !
         ,NstationNum_GL    &    !maximum number of station is NstationMax!
         ,t_stn_start    &        !
         ,t_stn_end    &            !
         ,t_stn_dlt    &            !starting time, ending time, and interval of station outputs (days)
         ,STNFN    &                !file name for station output
         ,HISFN    &                !file name for history output
         ,HISFN_PREFIX    &        !prefix of history output file
         ,HISFN_EXT    &            !extention name of history output file
         ,HISFN_FINAL    &        ! 
         ,HISFN_SPLIT_BYLEVEL & !True or False for splitting history output in files level by level (default is .FALSE.)
        !,hydro_netcdf    &        !
        !,wqm_history    &        !
        !,wqm_stations    &        !
        !,IFNC            &        !file number index for hydrodynamics netcdf files, set to hydro_filenumstart initially for cold start, set otherwise 
        !,NTRECNC        &        !time record index for a particular hydrodynamics netcdf file, reset to 1 upon opening new file. 
        !,NTHYDRO        &        !overall time record index for all netcdf files, increment by 1 each time a hydrodynamics record is read
         ,HYDRO_GEOM_DEALLOC &    !
         ,HYDRO_DEALLOC
        
    USE MOD_BCMAP, ONLY :    &
            !I_OBC_N,        &!
            !I_OBC_GL,        &!
            BCMAP_DEALLOC

       !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY :         &
            SERIAL          &           !!TRUE IF SINGLE PROCESSOR
            ,MSR            &           !!TRUE IF MASTER PROCESSOR (MYID==1)
            ,PAR            &           !!TRUE IF MULTIPROCESSOR RUN
            ,CASENAME      !    &           !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
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
            
       USE MOD_FILEINFO, ONLY : &
                DIA                &
                !,CBC             &
                !,S1            &
                !,S2            &
                !,S3            &                  
                !,BFI            &
                !,BAI           &
                !,MET            &
                ,BFO            & 
                !,KEI            &
                !,ATM            &
                !,STL            & 
                !,AGR            & 
                !,SVI            & 
                !,SVO            & 
                ,KFL            & 
                !,ZOO            & 
                !,ZFO            & 
                !,ALO              &
                !,CON            &
                ,RSO            &
                ,SNP            &
                ,PLT            &
                ,APL             &
                ,TFL             &
                ,OPL            &
                !,SFI            &
                ,SFO            &
                !,MAP             &
                !,ICI             &
                ,ICO            &
                !,MRL            &
                ,MBL            &
                !,RSI            &
                !,UNIT_LINKAGE    &
                ,UNIT_STN        &
                ,UNIT_HIS        &            
                ,CNAME            &
                !,CONFN
                ,INIT_FILE_INFO

    
    USE MOD_SF, ONLY:             &     !
                SFEEDER,         &    !
                !HYPOXFX_SF,        &    !
                !FRSASF,        &    !
                NSPECIES,        &    !
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
                !JNH4SF,        &!(:),            
                !JPO4SF,        &!(:),            
                !SODSF,            &!(:),            
                !JSASF,            &!(:),            
                !JSUSF,            &!(:),            
                !SF_SSI,        &!(:),            
                !SF_SU,            &!(:),        
                !SF_SA,            &!(:),            
                !SF_PIP,        &!(:),            
                SFEED,            &!(:,:),
                SFEED_GL,        &!(:,:),
                !SEDTYPE_GL,        &!(:,:)
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
                SUSPFEED,        &!        
                !SF_ALLOC        &!
                SF_DEALLOC!,    &!
                !SF_READ,        &!    
                !SF_INIT
                
    USE MOD_WQM, ONLY:             &!
                JDAY,           &!    
                DLT,            &!
                AHMDLT,         &!                
                DLTAV,          &!
                ELTMS,          &!
                ELTMS1,         &!
                ELTMJD,         &!

                WQ_DT,          &!    

                SEDIMENT_CALC,  &!
				total_netPP,   &
				total_netPP_gl,  & !

                T,              &!
                C2,                &!
                SALT,            &!
                SSI,            &!
                B1,                &!
                B2,                &!
                B3,                &!
          !      LDOC,            &!
           !     RDOC,            &!
                LPOC,            &!
                RPOC,            &!
                NH4,            &!
                NO3,            &!
                UREA,            &!
            !    LDON,            &!
           !     RDON,            &!

                PO4,            &!
             !   LDOP,            &!
            !    RDOP,            &!
                LPOP,            &!
                RPOP,            &!
                PIP,            &!
                COD,            &!
                DOXG,            &!
                SIUPB,            &!
                SIAT,            &!
                PIB1,            &!
                PIB2,            &!
                PIB3,            &!
                Q1,             &!
                Q2,             &!
                Q3,             &!
                BENTHIC_FLUXES, &!
                DSSR,           &!
                DTC,            &!
                FD,             &!
                FLOW,           &!
                I0,             &!
                IT,             &!
                LIGHT_EXTINCTION, &!
                RELTMS,         &!
                SAV_CALC,       &!
                SAVENT,         &!
                SEDENT,         &!
                TODS,           &!
                TTSS,           &!
                Z_DIFFUSION,    &!
                AC,             &!
                ANC1,           &!
                ANC2,           &!
                ANC3,           &!
                ASC1,           &!
                ASC2,           &!
                ASC3,           &!                
                C1,             &!
                C2F,            &!
                CMASS,          &!
                DTM,            &!

                FLUXS,          &!
                FNDLT,          &!
                KBM1,           &!
                MASS_BALANCE,           &!
                NAC,            &!
                NIT,            &!
                SFENT,          &!
                SOLIDS_CALC,            &!
                TMEND,          &!
                TMSTRT,         &!
                
                WCMC,           &!
                WCMN,           &!
                WCMP,           &!
                WCMS,           &!
                
                XYDF,           &!
                XY_DIFFUSION,   &!
                

                ATMFLXN,                &!
                ATMFLXNB,               &!
                ATMFLXP,                &!
                ATMFLXPB,       &!
                ATMLDON,        &!
                ATMLDOP,        &!
                ATMNH4,         &!
                ATMNO3,         &!
                ATMPO4,         &!
                ATMRDON,        &!
                ATMRDOP,        &!

                BENFLXDN,       &!
                BENFLXPN,       &!
                BENFLXPNB,      &!

                DLWCMC,         &!
                DLWCMN,         &!
                DLWCMP,         &!

                IWCMC,          &!
                IWCMN,          &!
                IWCMP,          &!

                PRECIP,         &!
                SFATOT,         &!
                
                ANDC,           &!
                BENFLXDNB,              &!
                BENFLXDP,               &!
                BENFLXDPB,              &!
                BENFLXPC,               &!
                BENFLXPCB,              &!
                BENFLXPP,               &!
                BENFLXPPB,              &!
                
                BENNH4,         &!
                BENNO3,         &!
                BENSA,          &!                
                BENPO4,         &!
                
                BNDN,           &!
           
                DLALGC,         &!
                DLWCKMC,                &!
                DLWCKMCB,               &!
                DLWCKMN,                &!
                DLWCKMNB,               &!
                DTMP1,          &!

                FLXTLDON,               &!
                FLXTLPON,               &!
                FLXTNH4,                &!
                FLXTNO3,                &!
                FLXTRDON,               &!
                INFLOW,         &!
                JCB,            &!
                !MNLLDOC,                &!
                !MNLRDOC,                &!
                PCFWS,          &!
                PPFWS,          &!
                RESPC,          &!

                BNDC,           &!
                BNDFLXC,                &!
                BNDFLXN,                &!
                BNDFLXP,                &!
                BNDP,           &!
                DLWCSSMN,               &!
                DLWCSSMP,               &!


                FLXT1,          &!
                FLXT2,          &!
                FLXT3,          &!
                FLXTLDOC,               &!
                FLXTLDOP,               &!
                FLXTLPOC,               &!
                FLXTLPOP,               &!
                FLXTLZ,         &!
                FLXTPIB1,               &!
                FLXTPIB2,               &!
                FLXTPIB3,               &!
                FLXTPO4,                &!
                FLXTRDOC,               &!
                FLXTRDOP,               &!
                FLXTRPOC,               &!
                FLXTRPON,               &!
                FLXTRPOP,               &!
                FLXTSZ,         &!
                S1FLXN,         &!
                S1FLXP,         &!
                S2FLXN,         &!
                S3FLXN,         &!

                BENDEN,         &!
                BURIALFLXC,             &!
                BURIALFLXCB,            &!
                BURIALFLXN,             &!
                BURIALFLXNB,            &!
                BURIALFLXP,             &!
                BURIALFLXPB,            &!
                DFTMP,          &!
                DLSEDKC,                &!
                DLSEDKCB,               &!
                DLSEDKN,                &!
                DLSEDKNB,               &!
                DLSEDMC,                &!
                DLSEDMN,                &!
                DLSEDMP,                &!
                DLSEDSSMC,              &!
                DLSEDSSMN,              &!
                DLSEDSSMP,              &!
                DLWCSSMC,               &!
                DTMP2,          &!
                ERRSEDN,                &!
                ERRWCC,         &!
                ERRWCN,         &!
                ERRWCP,         &!

                S1FLXC,         &!
                S2FLXC,         &!
                S2FLXP,         &!
                S3FLXC,         &!
                S3FLXP,         &!

                C2_GL,          &!
                CCHL1,          &!
                CCHL1_GL,               &!
                D_GL,           &!
                EL_GL,          &!
                ERRSEDC,                &!
                ERRSEDP,                &!

                H_GL,           &!
                MGL,            &!
                MTLOC,          &!
                P1,             &!
                S_GL,           &!
                T_GL,           &!
                
                B1SZ_GL,                &!
                BM1,            &!
                BM1_GL,         &!
                BM2,            &!
                BM2_GL,         &!

                FI1,            &!
                FI1_GL,         &!
                FI2,            &!
                FI2_GL,         &!
                NL1,            &!
                NL1_GL,         &!
                NL2,            &!
                NL2_GL,         &!
                P1_GL,          &!
                P2,             &!
                P2_GL,          &!
                PL1,            &!
                PL1_GL,         &!
                PL2,            &!
                PL2_GL,         &!

                B1LZ_GL,        &!
                B2LZ_GL,        &!
                B2SZ_GL,        &!
                DCOD,           &!
                DCOD_GL,        &!
                DDOC,           &!
                DDOC_GL,        &!

                IAVG_GL,        &!
				KDLAMB_int_GL,   &!B Clark may 2018
				SSI_GL,  & ! B Clark June 2018
				
       RHO1_GL,RHO1,ShearStress,ShearStress_GL,  &  ! B Clark June 2018
                IK1,            &!
                IK1_GL,         &!
                IK2,            &!
                IK2_GL,         &!
                NITRIF,         &!
                NITRIF_GL,              &!
                PR1,            &!
                PR1_GL,         &!
                PR2,            &!
                PR2_GL,         &!

                BFLUX,          &!
                BFLUX_GL,       &!
                COD_GL,         &!
                COURMX,         &!
                DIAGNOSTICS,    &!
                DOLZ_GL,        &!
                DOPR,           &!
                DOPR_GL,        &!
                DORALG,         &!
                DORALG_GL,      &!
                DOSZ_GL,        &!
                FTCOD,          &!
                FTCOD_GL,       &!
                REAERDO,        &!
                REAERDO_GL,     &!
                RESP,           &!
                RESP_GL,        &!
                PNFWS,          &!                
                CTEMP,             &!
                DIFFMX,            &!

                B,              &!
                BB,             &!                
                IJDAY,          &!
                LENGTH,         &!
                SAV_LOADS,      &!
                WQM_DEALLOC,     &!
				
				! B Clark added below for outputting budget terms, July 2018
				
				ALGDON_OUT,    &!Algae DON production
				ALGDON_GL,    &! Global Algae DON production
				
				ALGPredDON, &
				ALGPredDON_GL, &
				ALGPPDON,  &
				ALGPPDON_GL, &!
				
				ALGPON_OUT,    &!Algae PON production from predation
				ALGPON_GL,    &! Global Algae PON production
				
				ALGNH4_OUT,  & ! Algal delta NH4
				ALGNH4_GL ,   &!
				
				ALGNO3_OUT ,   &!  Algal delta NO3
				ALGNO3_GL,  &  
				
				DENNO3_OUT,   & ! Denitrification loss of NO3
				DENNO3_GL,   &!
				
				NT_OUT,   & ! Nitrification loss rate of NH4 (production rate of NO3 from nitrification)
				NT_GL,   &
				
				DENITN_OUT,   & ! Denitrificaion loss rate of DON (remin of DON from denitrification)
				DENITN_GL,  & !
				
			    DENIT_OUT,          & ! loss of DOC due to denitrification
				DENIT_GL,   & !
				
				LPON,            &!
				
                RPON,            &!
		
				
				HDRLPON_OUT,  & ! Hydrolysis of LPON to DON
				HDRLPON_GL,  & !
				
				HDRRPON_OUT,  & ! Hydrolysis of RPON to DON
				HDRRPON_GL,  & ! 	
				
				ALGDOC_OUT ,   & ! Algae DOC production
				
				ALGDOC_GL,    & ! Global Algae DOC production,

			
				
				ALGPOC_OUT,    & !Algae POC production from predation
				ALGPOC_GL,    & ! Global Algae POC production
				
				HDRLPOC_OUT,  & ! Hydrolysis of LPOC to DOC
				HDRLPOC_GL,  &! 
				
				HDRRPOC_OUT,  & ! Hydrolysis of RPOC to DOC
				HDRRPOC_GL, &!
				DTLPOC_out, &!
				DTRPOC_out, &!
				DTRPON_out, &!
				DTLPON_out, &!
				DTNH4_out, &!
				DTNO3_out
				
				

	USE WC_DOM, ONLY: COAGN_OUT, COAGN_GL, COAGC_OUT, COAGC_GL,   &!
					  MNLDOC, &
					  MNLDOC1_OUT,MNLDOC2_OUT,MNLDOC3_OUT,   &!
					  MNLDOC1_GL,MNLDOC2_GL,MNLDOC3_GL,   &!
					  MNLDON1_OUT,MNLDON2_OUT,MNLDON3_OUT,   &!
					  MNLDON1_GL,MNLDON2_GL,MNLDON3_GL, &
					  DTWCDOC1_GL,DTWCDOC2_GL,DTWCDOC3_GL,  &
					  DTWNCDOC1_GL,DTWNCDOC2_GL,DTWNCDOC3_GL,  &
					  DTWCDON1_GL,DTWCDON2_GL,DTWCDON3_GL,  &
					  DTWNCDON1_GL,DTWNCDON2_GL,DTWNCDON3_GL, &
					  DTNH4_GL,DTNO3_GL,DTLPON_GL,DTLPOC_GL,DTRPON_GL,DTRPOC_GL, &
					  PHOTOCHEM_DOC_ACCUM_OUT, &
					  	DTWCDOC1_out,DTWNCDOC1_out,DTWCDOC2_out,DTWNCDOC2_out,   &
						DTWCDOC3_out,DTWNCDOC3_out,   & 
						DTWCDON1_out,DTWNCDON1_out,DTWCDON2_out,DTWNCDON2_out,   &
						DTWCDON3_out,DTWNCDON3_out,   & 
						DTWCDOP1_out,DTWNCDOP1_out,DTWCDOP2_out,DTWNCDOP2_out
					  
    USE MOD_OWQ, ONLY:                  & 
                        !INITKE,        &!
                        INTKE,            &!
                        NXKEI,            &!
                        !LAT,             &!
                        !LONG,            &!
                        !GSTAR440,        &!
                        !PHISTAR676,    &!
                        !BNSTARNTU,        &!
                        !BNSTARCHL,     &!
                        !PSTARINT,        &!
                        !PSTAR440,        &!
                        !PSTARCOLOR,    &!
                        !SCDOM,            &!
                        !EZEROINT,        &!
                        !TCHL2CHL,        &!
                        !DOFFSET,        &!
                        !DEPTHKE,        &!
                        !NWAVEL,        &!
                        !COLOR,         &!
                        !TURB,            &!
                        !OPZ,             &!
                        !G1,            &!
                        !G2,            &!
                        !ALAMB,         &!
                        !GLAMB,            &!
                        !PLAMB,            &!
                        !PHILAMB,       &!
                        !BLAMB,            &!
                        !WAVEL,            &!
                        !EZERO,            &!
                        !GMUZERO,        &!
                        IAVG,            &!
                        KESS,            &!
                        !IATBOT,        &!            
                        !OWQ_ALLOC,        &!
                        OWQ_DEALLOC,    &!
                        !OWQ_READ,        &!
                        LGHT_ATTN, &
						KDLAMB_int, & ! Added by B Clark May 2018 to output integrated kd for PAR
						SSI   ! B Clark June 2018
                        
    USE MOD_WQMINIT, ONLY :            &        !
            RSOD,                &        !    
            RSODP,                 &        !
            RSOFN,                 &        !
            DLTD,                &        !
            DLTDP,               &        !
            !SNPDP,                &        !
            !TFLDP,                &        !
            !PLTDP,                &        !
            !APLDP,               &        !
            !OPLDP,                &        !
            DIAD,                 &        !            
            DIADP,                &        !
            DIAF,                 &        !
            
            !COURB,                &        !
            COURFS,            &        !
            COURBS,            &        !
            F,                    &        !            
            SB,                    &        !
            !S1LNMAX,            &        !
            !S2LNMAX,            &        !
            !S3LNMAX,            &        !
            DIFFFS,                &        !
            DIFFBS,              &        !
                
            !NXSNP,             &        !
            !NXPLT,                &        !
            !NXAPL,                &        !
            !NXTFL,                &        !
            !NXKFL,                &        !
            NXTVD,               &        !
            !NXOPL,                &        !
            !NXMBL,                &        !
            NXDIA,                &        !

            !LEAFI,                &        !
            !STEMI,                &        !
            !ROOTI,             &        !
            !TUBERI,            &        !
            
            !MAXDLT,            &        !
            MXDLT,                &        !
            !INTKE,                &        !
            !INITKE,            &        !
            TM1,                &        !
            TM2,                &        !
            !TM3,                &        !
            !MINSTEP,            &        !
            !NXKEI,                &        !
                        
            !SNPC,                &        !
            !RSIC,                 &        !
            !RSOC,                &        !
            !BCC,                &        !
            !S1C,               &        !
            !S2C,                &        !
            !MDC,                &!*3,    !    
            !PLTC,                &        !
            !FLC,               &        !
            !MBLC,                &        !
            !BFC,                &        !
            !VBC,                &        !
            !QPLTC,                &        !
            !XYDFC,                &        !
            !XYDFU,              &        !
            !ZDFC,                &        !
            !ICOC,                &        !
            !ATMC,                &        !
            !SAVLC,                &        !
            !SEDC,              &        !
            !AUTOC,                &        !
            !SPLTC,                &        !
            !TFLC,                &        !
            !DIAC,                &        !
            !STLC,              &        !
            !APLTC,                &        !
            !KFLC,                &        !
            !OPLC,                &        !
            !BFOC,                &        !
            !DFOC,                &        !
            !S3C,                  &        !
            !SAVMC,                &        !
            !SAVPLTC,            &        !
            !DFLC,                &        !
            !SFLC,                &        !
            !SFLOXC,            &        !
            !DFLOXC,            &        !
            EXT1,                 &!*1,   !    
            EXT2,                &!*2,    !    
            EXT3,                &!*3,    !    
            EXT4,                &!*4,    !    
            !EXT_R,                &!*1,    !
            !EXT_R2,            &!*2    !
             
            !SLC,                &        !
            !HYDC,                &        !
            !BNDTC,                &        !
            !CONSC,                &        !
            !ICIC,                &        !
                                
            !TIMEA,                &        !
            !TIMEB,                &        !
            !SPVARM,            &        !
            !PRINTM,            &        !
                                
            TITLE,                &!(6),    !
            !OLDTITLE            &!(6),    !
            FILENAME,            &        !
                                    
            !FILENAME2,            &        !
                
            !MAPFN,                &        !
            !GEOFN,                &        !
            !ICIFN,                &        !
            !RSIFN,              &        !
            !AGRFN,                &        !
            !STLFN,                &        !
            !MRLFN,                &        !
            !KFLFN,                &        !
            !ICOFN,                &        !
            
            !SNPFN,                &        !
            !PLTFN,                &        !
            !APLFN,             &        !
            !DIAFN,             &        !
            !TFLFN,              &        !

            !OPLFN,             &        !
            !MBLFN,             &        !
            !SFIFN,             &        !
            !SFOFN,                &        !

            RESTART_OUT,        &        !
            !SNAPSHOTS,            &        !
            END_RUN,            &        !
            !MODIFY_ICONC,       &          !
            VOLUME_BALANCE,     &        !
            !QUICKEST,            &        !
            !UPWIND,            &        !
            ICOND_OUT,          &        !
            !UNI_ICON_IN,       &        !
            !BIN_ICON_IN,        &        !
            !RES_ICON_IN,        &        !
            AUTO_STEPPING,      &        !
            !STOP_RUN,           &        !
            !PLOTS,                &        !
            !OXYGEN_PLOTS,       &        !
            
            NEW_VOLUMES,          &        !
            !RESTART_IN,        &        !
                                
            !TEMPERATURE_CALC,     &        !
            ALGAE_CALC,        &        !
            ZOO_CALC,           &        !
            CARBON_CALC,       &        !
            NITROGEN_CALC,     &        !
            PHOSPHORUS_CALC,   &        !
            COD_CALC,            &        !
            OXYGEN_CALC,        &        !
            SILICA_CALC,        &        !
                                
            LEFT_FLOWB,         &!(NHQP),        !        
            !RIGHT_FLOWB,        &!(NHQP),        !    
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

            IFLOWP,            &    !(NBCP),        !
            !IBT,                &    !(NHQP),          !
            !COURBSV,            &    !(NHQP),        !
            !COURVSV,            &    !(NHQP)            !

            C1MIN,              &                !
            C1MAX,                &                !

            DOVDAYS,            &    !(:,:,:)    !
            !OINT,                &    !(NOIP)        !
            ACC,                &    !(NCP)        !
             
            !CIC,                &    !(NCP),        !
            !COUT,                &    !(NCP),        !

            !SFEEDI,            &    !(10)        !
            
            !SNPD,                &!
            !SNPF,                 &!
            
            DLTVAL,                &!
            DLTMAX,               &!
            DLTFTN,             &!
            !PLTD,                &!
            !PLTF,                 &!
            !APLTD,                &!
            !APLF,                &!
            !TFLD,                 &!
            !TFLF,                &!
            !KFLD,                &! 
            !KFLF,                 &!
            !OPLD,                &!
            !OPLF,                &!
            

            !MBLD,                 &!
            !MBLF,                &!

            !DOCLIT ,             &!
            !LPOCLIT,              &!
            !POCLIT,              &!
            !PBSLIT,              &!
            !DSILLIT,            &!
            !WQMINIT_ALLOC,        &!Subroutine to initialize MOD_WQM data
             WQMINIT_DEALLOC
             
   USE MOD_ALGAL, ONLY :         &
            !ALG_READ,            &
            ALGAE    
   
   USE MOD_ZOOP, ONLY :            &
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
                     CLSZ,        & 
                     CLLZ,        &
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
                      B1SZ,        &!
                      B2SZ,        &!
                     !B3SZ,        &!
                     B1LZ,        &!
                     B2LZ,        &!
                     !B3LZ,        &!
                     DOSZ,        &!
                     DOLZ,        &!
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
                     SZ,        &!
                     LZ,        &!
                     ZOOP_DEALLOC
                     
    USE MOD_KIN, ONLY :                &
                        TEMPER,        &
                        SOLIDS,        &
                        ZOOPL,        &
                        PARTICULATE_CARBON,        &
                        NITROG,        &
                        PHOSPH,        &
                        CODMND,        &
                        OXYGEN,        &
                        SILICA,        &
                        BEN_FLUX

   USE MOD_SED, ONLY :         &!
       !    BENTHIC_OUTPUT,     &!
        STEADY_STATE_SED_G3,&!
        STEADY_STATE_SED_IC,&!
        DLTS,                &!
    !    SSTATEG3,           &!
    !    HSEDALL,               &!
    !    DIFFT,                &!
    !    SALTSW,             &!
    !    SALTND,                &!    
    !    FRPALG1,            &! 
    !      FRPALG2,            &! 
    !      FRPALG3,            &! 
    !    FRNALG1,            &! 
    !   FRNALG2,               &! 
    !   FRNALG3,               &!
    !    FRCALG1,            &!
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
        IERR,                &!
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
        ISEDMN,               &!
        ISEDMP,               &!
        ISEDMC,                &!
        SEDMN,              &!total sediment organic nitrogen  (kgN)
        SEDMP,              &!total sediment organic phosphorus (kgP)
        SEDMC,                &!total sediment organic carbon  (kgC)
        !NH4AVL,               &!    
        !NO3AVL,               &!
        !KETOTL,            &!
        !WSSINETMMD,         &!
        !WLPOMNETMMD,        &!
        !WRPOMNETMMD,        &!
        !WB1NETMMD,            &!
        !WB2NETMMD,            &!   
        !WB3NETMMD,            &!   
        !WPSINETMMD,        &!   
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
    !    BFORMAXS,           &!      
    !    ISWBENS             &!    
        SODTM1S,            &
        JNH4TM1S,     & 
        JNO3TM1S,     & 
        JPO4TM1S,     & 
        JCH4TM1S,     & 
        JCH4GTM1S,    & 
        JHSTM1S,    &
        JSITM1S,    &
        NH41TM1S,              &!  
        NH42TM1S,              &!  
        NH4T2TM1S,             &!  
        NO31TM1S,          &! 
        NO32TM1S,          &! 
        NO3T2TM1S,             &!  
        HS1TM1S,           &! 
        HS2TM1S,           &! 
        HST2TM1S,              &! 
        SI1TM1S,           &! 
        SI2TM1S,           &! 
        SIT2TM1S,              &!  
        PO41TM1S,          &!  
        PO42TM1S,          &!  
        PO4T2TM1S,             &!  
        CH41TM1S,              &!  
        CH42TM1S,              &!  
        CH4T2TM1S,             &!  
        SO4T2TM1S,            &!     
        BURIALC,            &!     
        BURIALN,            &!     
        BURIALP,               &!
        DIAGENC,              &!   
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
    !   BBM,                  &!  
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
    !    PSISEDTM1             &!
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
        BENDEN_GL,            &!
        JPOC_GL,            &
        JPON_GL,            &
        JPOP_GL,            &
        JPOS_GL,            &
        CH40_GL,            &
        JNH4_GL,            &
        JNO3_GL,            &
        JCH4_GL,            &
        JCH4G_GL,            &
        JHS_GL,                &
        JPO4_GL,            &
        JSI_GL,                &
        NH41_GL,            &
        NH42_GL,            &
        NO31_GL,            &
        NO32_GL,            &
        PO41_GL,            &
        PO42_GL,            &
        SI1_GL,                &
        SI2_GL,                &
        CH42_GL,            &
        HS1_GL,                &
        HS2_GL,                &
        HSED1_GL,            &
        HSED1,                &
        
       JPOC,         &!
      JPON,         &!
	  
	  JPOC_out, &!
	  JPON_out, &!
	  JNH4_out, JNO3_out, &
      JPOP,         &!                
      JPOS,         &!
!        CCH4,           &!
!       CHS,            &!
!       CSI,            &!
!       CSO4,           &!                
        CPOC,               &!
        CPON,               &!
        CPOP,               &!
        CPOS,               &!
        CPO4,                 &!
        CNH4,                 &!
        CNO3,                 &!                
        CHS,                 &!
        CSI,                 &!                
        DIAGN,                 &!                    
        !SED_INIT,            &!
        SED_INIT2,            &!
    !    SED_ALLOC,            &!
        SED_DEALLOC,        &!
    !    SED_READ,            &!
        SED_CALC,            &!
        SED_DIAGENESIS_G3!,            &!
    !    SEDTSFNL,            &!
    !    SEDSSFNL,            &!
       !    SEDF,                &!
    !    ZBRENT        

   USE MOD_SED_DF_EXCHANGE_VARS, ONLY:     SED_DF_EXCHANGE_ALLOC,    &!
                                        SED_DF_EXCHANGE_DEALLOC
   
   USE MOD_SED_SF_EXCHANGE_VARS, ONLY:     SED_SF_EXCHANGE_ALLOC,    &!
                                        SED_SF_EXCHANGE_DEALLOC
   
   USE MOD_SED_SAV_EXCHANGE_VARS, ONLY: SED_SAV_EXCHANGE_ALLOC,    &!
                                        SED_SAV_EXCHANGE_DEALLOC
   
   USE MOD_BA, ONLY : BALGAE_CALC, BBM,BBM_GL, BA_DEALLOC
   
   USE MOD_OBCS, ONLY:  &
        !NEXT_OBC ,       & !!INTERIOR NEIGHBOR OF OPEN BOUNDARY NODE
        !NEXT_OBC2,       & !!INTERIOR NEIGHBOR OF NEXT_OBC
        !!TYPE_OBC ,       & !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
        !TYPE_OBC_GL,     & !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
        !!IBCN       ,     & !!NUMBER OF EACH TYPE OF OBN IN LOCAL  DOM
        !!IBCN_GL    ,     & !!NUMBER OF EACH TYPE OF OBN IN GLOBAL DOM
        !OBC_LST    ,     & !!MAPPING OF OPEN BOUNDARY ARRAYS TO EACH TYPE
        !NADJN_OBC  ,     & !!NUMBER OF ADJACENT OPEN BOUNDARY NODES TO OBN
        !ADJN_OBC   ,     & !!ADJACENT OBN's of OBN
        !NADJC_OBC  ,     & !!NUMBER OF ADJACENT OPEN BOUNDARY CELLS TO OBN
        !ADJC_OBC   ,     & !!ADJACENT OPEN BOUNDARY CELLS
        !!NFLUXF_OBC ,     & !!NUMBER OF FLUX SEGMENTS TO OBN
        !!FLUXF_OBC  ,     & !!FLUX FRACTION ON EACH SIDE OF OBN
        !NXOBC      ,     & !!INWARD POINTING X-NORMAL OF OBN
        !NYOBC      ,     & !!INWARD POINTING Y-NORMAL OF OBN
        !DLTN_OBC   ,     & !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY
        !NUT_OBC    ,     & !!NUTRIENT AT OPEN BOUNDARY
        !NUT_OBC_GL ,     & !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY

        !SEPARATE_OBC,     &!
        SETUP_OBC ,        &!
        OBCS_DEALLOC    
   
   USE MOD_BCS, ONLY : &
        ! INODEQ            &  !!LOCAL FRESH WATER INFLOW NODES
        !,ICELLQ         &  !!LOCAL FRESH WATER INFLOW ELEMENTS 
        !,RIV_GL2LOC      &  !!GLOBAL TO LOCAL MAP OF FW NODES
        !,INOPNT          &  !!LOCAL NON-POINT SOURCE NODES
        !,QBC_TM         &  !!TIME MAP FOR RIVER DISCHARGE DATA 
        !,PNT_TM        &  !!TIME MAP FOR NON POINT SOURCE DATA
        !,VQDIST        &  !!DISCHARGE VERTICAL DISTRIBUTION
        !,WVQDIST       &  !!DISCHARGE VERTICAL DISTRIBUTION for point source
        !,DQDIS         &  !!WATER FLUX DISCHARGE DATA
        !,DWDIS           !&  !!WATER QUALITY DISCHARGE DATA
        !,WDIS          &  !!FRESH WATER QUALITY AT CURRENT TIME
        !,PWQDIS        &  !!Current WATER QUALITY at non-point source node
        !,PDWQDIS         &  !!WATER QUALITY DATA at non-point source node
        !,PQDIS           &  !!Current DISCHARGE at non-point source node
        !,PDQDIS         &  !!DISCHARGE at non-point source node
        BCS_DEALLOC

  USE MOD_SAV, ONLY    :        &
            !SPNAM,           &    !species name of the SAV
            !EPINAM,          &   !name of the epiphyte
            !NSAVM,           &      !global maximum number of SAV species in all of tracer control elements (TCE's) - polygon
                                !surrounding a given node
            
            !KHSTOX,         &   !Sulfide toxicity coefficient in calculating HS^- limitation on SAV growth (m^3/mole-S2)
            !KESAV,      &   !Effect on light attenuation by SAV leaf and stem, KESAV (gC/m^3)  
            !KEEPI,      &   !Shading effect on light on SAV shoots (leaf) due to Epiphytes ((m^2-leaf-surface/gDW-epiphyte)
            !KHEP,       &  !half saturation desnity of epiphyte (gC-Epiphyte/gC-leaf) for epiphyte growth on SAV leaf
            !NLIMEPI,    &  !Never used ???                      
             
            !KHNLEAF,    &  !Half saturation constant for N uptake by leaf (shoots) (gN/m^3)
            !KHNROOT,    &  !Half saturation constant for N uptake by root (gN/m^3)
            !KHPLEAF,    &  !Half saturation constant for P uptake by leaf (gP/m^3)
            !KHPROOT,    &  !Half saturation constant for P uptake by root (gP/m^3)
            !KHNEPI,     &  !Half saturation constant for N uptake by epiphytes (gN/m^3)
            !KHPEPI,     &  !Half saturation constant for P uptake by epiphytes (gN/m^3)
            !KHNPSAV,    &  !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for SAV
            !KHNPEPI,    &  !Threshold for NH4 in calculating NH4 preference over NO3 in N uptake(gN/m^3) for epiphytes

            !PMSAV,      &   !maximum production rate of SAV, gC/gDW/day  where DW is dry weight
            !PRSPSAVB,   &   !fraction of leaf production rate that is counted as photochemical respiration at salt << saltmax
                            !at salt >> saltmax, PRSPSAV approaches 1, and there will be no production of SAV leaf hence
                            !no photosynthesis for SAV growth
            !BMSAV,      &   !Basal metabolism rate of SAV (1/day)
            !BMTBRREF,   &   !Reference basal metabolism rate of tube [1/day] BMTUBE=BMTBRREF*f(T), i.e. function of temperature
            !FDOSR,      &   !Fraction of DOXG production discounted due to route of Leaf production to root
            !SALMAX,     &   !Salinity at which PRSPSAV is approaching 1 (maximum) from PRSPSAVB, [psu]
            !ALPHSAV,    &   !Initial slope of P-I curve for SAV [ (gC/gDW)/(E/m^2) ]
            !ALAC,       &   !m^2-leaf-area/gC-leaf (leaf area to leaf carbon ratio) 
                            !Note in 2002 version, it is the inverse ACLA: gC-leaf/m^2-leaf-area
            !SLSAV,      &   !???? Never used !Sloughing rate of SAV ?  1/day
            !ACDWSAV,    &   !Carbon to Dry Weight ratio of SAV (gC/gDW)
            !ANDWSAV,    &   !Nitrogen to Dry Weight ratio of SAV (gN/gDW)
            !APDWSAV,    &   !Phosphorus to Dry Weight ratio of Sav (gP/gDW)
            !ANCSAV,     &   !Nitrogen to Carbon ratio of SAV (gN/gC)             
            !APCSAV,     &   !P to C ratio of SAV (gP/gC)
            !HCAN,       &   !Canopy height (m), HCAN = ACAN + BCAN *( leaf+stem )
            !ACAN,       &   !basic canopy height (m)
            !BCAN,       &     !coefficient for calculating leaf (gC/m^2) and stem (gC/m^2) 
                             !contribution to canopy height,  [m^3/gC]

            !FNISAV,       &   !Fraction of NH4 Nitrogen generated due to SAV (leaf, stem etc) respiration (mortality)
            !FNLDSAV,      &   !Fraction of Nitrogen routed to LDON on SAV mortality
            !FNRDSAV,      &   !Fraction of Nitrogen routed to RDON on SAV mortality
            !FNLPSAV,      &   !Fraction of Nitrogen routed to LPON on SAV mortality
            !FNRPSAV,      &   !Fraction of Nitrogen routed to RPON on SAV mortality
            !FPISAV,       &   !Fraction of Phosphorus routed to PO4 in water column due to SAV (leaf,stem) respiration
            !FPLDSAV,      &   !Fraction of P routed to LDOP in water column due to SAV mortality
            !FPRDSAV,      &   !Fraction of P routed to RDOP in water column due to SAV mortality
            !FPLPSAV,      &   !Fraction of P routed to LPOP in water column due to SAV mortality
            !FPRPSAV,      &   !Fraction of P routed to RPOP in water column due to SAV mortality
            !FCLDSAV,      &   !Fraction of Carbon routed to LDOC in water column due to SAV mortality
            !FCRDSAV,      &   !Fraction of Carbon routed to RDOC in water column due to SAV mortality
            !FCLPSAV,      &   !Fraction of Carbon routed to LPOC in water column due to SAV mortality
            !FCRPSAV,      &   !Fraction of Carbon routed to RPOC in water column due to SAV mortality
            !FDOSAV,     &   !Fraction of Carbon routed to CO2 (oxygen consumption) in water column due to SAV mortality
                            !Note FDOSAV + FCLDSAV+ FCRDSAV + FCLPSAV + FCRPSAV = 1
                            !Wen Long: we really need to rename FDOSAV to FCO2SAV and make it equivalent to FDOSAV!

            !WSSSAV,       &   !Increase of fixed solids (ISS) settling by SAV leaf and stem, (m/day)/(gC/m^2)==>m^3/gC/day
            !WSLSAV,       &   !Increase of LPOM settling by SAV leaf and stem, m^3/gC/day
            !WSRSAV,       &   !Increase of RPOM settling by SAV leaf and stem, m^3/gC/day
            !WS1SAV,       &   !Increase of Algae group 1 settling by SAV leaf and stem, m^3/gC/day
            !WS2SAV,       &   !                        2
            !WS3SAV,       &   !                        3
            !WSUSAV,       &   !Increase of particulate biogenic silica settling rate by SAV leaf and stem, m^3/gC/day
                     
            !NDDI,         &     !number of degree days for entries of look up table of FPLEAF etc for Senescence of leaf etc
            !NSAVCELL,     &     !total number of grid cells (ICM boxes) that have SAV --> need to be changed to nodes of FVCOM
            !NSAVSPCM,     &   !maximum number of species in a cell
              
            !PMEPI,      &     !maximum photosynthesis production rate (gC/gCHLA/day) of epiphytes
            !BMEPI,      &     !basal metabolism loss rate of epiphyte (1/day)
            !PREPI,      &     !predation rate on epiphytes (gC-leaf/gC-epiphyte/day)
            !PRSPEPI,    &     !Photochemical respiration rate of epiphyte --> BMEP

            !ALPHEPI,    &     !slope of PI curve for epiphyte photothenesis (gC/gDW)/(E/m^2)
            !CCHLEPI,    &     !Carbon to Chlorophyll ratio of epiphyte (gC/gCHLA)
            !ANCEPI,     &     !N to C ratio of epiphytes (gN/gC)
            !APCEPI,     &     !P to C ratio of epiphytes (gP/gC)
            !ADWCEPI,    &      !ratio of total epiphyte dry weight to viable carbon (gDW/gC)
             
            !FCLDEPI,    &   !Fraction of C routed to LDOC due to basal metabolism of epiphyte
            !FCRDEPI,    &   !Fraction of C routed to RDOC due to basal metabolism of epiphyte
                            !WLong: missing FCLPEPI and FCRPEPI to routt to LPOC and RPOC ??? (seems not coded)
                            !FCLPEPI and FCRPEPI should be added (Wen Long)
                            !also missing FDOEPI compared to FDOSAV for SAV
                            !code assumes that no carbon will be routed to LPOP and RPOP due to metabolism
                            !and also FDOEPI=1-FCLDEPI-FCRDEPI, and FDOEPI is only implicitly used in calculating 
                            !DOEPIW
             
            !FNIEPI,     &   !Fraction of N routed to NH4  due to epiphyte basal metabolism
            !FNLDEPI,    &   !Fraction of N routed to LDON due to epiphyte basal metabolism
            !FNRDEPI,    &   !Fraction of N routed to RDON due to epiphyte basal metabolism
            !FNLPEPI,    &   !Fraction of N routed to LPON due to epiphyte basal metabolism
            !FNRPEPI,    &   !Fraction of N routed to RPON due to epiphyte basal metabolism
            !FPIEPI,     &   !Fraction of P routed to  PO4 due to basal metabolism of epiphyte
            !FPLDEPI,    &   !Fraction of P routed to LDOP due to basal metabolism of epiphyte
            !FPRDEPI,    &   !Fraction of P routed to RPOP due to basal metabolism of epiphyte
            !FPLPEPI,    &   !Fraction of P routed to LPOP due to basal metabolism of epiphyte
            !FPRPEPI,    &   !Fraction of P routed to RPOP due to basal metabolism of epiphyte             

            !FCLDPEP,    &   !Fraction of C routed to LDOC due to predation and sloughing on epiphytes
            !FCRDPEP,    &   !Fraction of C routed to RDOC due to predation and sloughing on epiphytes
            !FCLPPEP,    &   !Fraction of C routed to LPOC due to predation and sloughing on epiphytes
            !FCRPPEP,    &   !Fraction of C routed to RPOC due to predation and sloughing on epiphytes
            !FNIPEP,     &   !Fraction of Nitrogen routed to NH4  due to predation on epiphytes and Slough
            !FNLDPEP,    &   !Fraction of Nitrogen routed to LDON due to predation on epiphytes and Slough
            !FNRDPEP,    &   !Fraction of Nitrogen routed to RDON due to predation on epiphytes and Slough
            !FNLPPEP,    &   !Fraction of Nitrogen routed to LPON due to predation on epiphytes and Slough
            !FNRPPEP,    &   !Fraction of Nitrogen routed to RPON due to predation on epiphytes and Slough
            !FPIPEP,     &   !Fraction of P routed to PO4  due to predation on epiphytes and slough
            !FPLDPEP,    &   !Fraction of P routed to LDOP due to predation on epiphytes and slough
            !FPRDPEP,    &   !Fraction of P routed to RDOP due to predation on epiphytes and slough
            !FPLPPEP,    &   !Fraction of P routed to LPOP due to predation on epiphytes and slough
            !FPRPPEP,    &   !Fraction of P routed to RPOP due to predation on epiphytes and slough            

            !FPSR,          &    !Never used!!!! This is supposed to be the fraction of leaf production routed to root
                            !(replaced by FPROOT)

     
            !FPLEAF,     &     !fraction of production assigned to leaf growth
            !FPROOT,     &     !fraction of production assigned to root growth
            !FPSTEM,     &     !fraction of production assigned to stem growth
            !FPTUBER,    &     !fraction of production assigned to tuber growth
                            !WLong: FPLEAF +FPROOT+FPSTEM+FPTUBER = 1
                            !all growth is based on photosynthesis by leaf
                            !(so stem and tuber have no photosynthesis?? I thought they look green )
            !TRTBRLF,    &   !transfer rate from tuber to leaf (gC-Leaf/day/gC-Tuber)
                               
            !!Look up table of temperature control in a year on growth and mortality rates of SAV and epiphytes
            !FTPSAV,     &     !temperature control on SAV growth
            !FTRSAV,     &     !temperature control on SAV resipration 
            !FTPEP,      &     !temperature control on epiphyte growth rate
            !FTREP,      &     !temperature control on epiphyte resipration (basal metabolism)
            !FTPREP,    &   !temperature control on predation on epiphyte
                               
            !!fluxes to sediments by SAV mortality
            !FRPOCSAV,     &    !
            !FRPONSAV,     &    !
            !FRPOPSAV,  &    !3G partitioning of POC, PON, POP sources to sediments from SAV (non-dimensional)
  
            !!State variables of SAV model 
            !LEAF,       &     !Leaf biomass (gC/m^2) 
            !ROOT,       &     !Root biomass (gC/m^2)
            !STEM,       &     !Stem biomass (gC/m^2)
            !TUBER,      &     !Tuber biomass (gC/m^2)
            !EP,        &   !Epiphytes biomass (gC-epiphytes/gm^2-leaf)
                            !Note in 2002 version of ICM, EP's unit was (gC-epiphytes/gC-leaf)
                                 
            !SAVAREA,    &   !area of SAV coverage for each surface cell (m^2) --> needs to be changed for FVCOM grid                                 
            !SAVCELL,    &   !grid cell number for all cells that have SAV in it
                            !size 1xNSAVCELL --> need to change to element number or node number!!
            
            !NSAVSPC,    &    !(:) !number of SAV species in each grid cell
                            !size 1xNSAVCELL

            !SAVDPH,      &      !Depth of SAV (mean depth of the plant) for all species in each grid cell
                            !size, NSAVCELL x NSAVSPC(B) where B is cell number, B=1,...,NSAVCELL
            !SAVFRAC,    &      !coverage for each SAV cell and each vertical deph increment
                            !i.e. coverage is function of grid cell and SAV species

                            !Growth and mortality rates of SAV
            !PLEAF,        &    !Leaf growth rate (1/day), P (SAV shoots growth) in C-M 2001 paper's SH equation
            !BMLEAF,     &      !respiration rate of leaf of SAV (1/day) with temperature control 
            !BMSTEM,     &      !resipration rate of stem of SAV (1/day) with temperature control
            !BMROOT,     &      !respiration rate of root of SAV (1/day) with temperature control
            !BMTUBER,   &    !respiration rate of tuber of SAV (1/day) with temperature control                                      
     
            !Growth and mortality rates of SAV                                      
            !PEP,       &      !epiphytes growth rate (1/day)  = PMEPI*light_limitation*deinsity_limitation/CCHLEPI
            !BMEP,      &      !epiphytes basal metabolism rate (1/day) (also including photochemical respiration)
            !PREP,      &   !predation rate on epiphytes (1/day) = PREPI*EP
                            !where PREPI has unit gC-shoot/gC-epiphyte/day

            !SLSH,        &    !(:)      !epiphyte sloughing rate (1/day)
     
            !!Nutrient limitations on SAV and epiphytes growth
            !NLSAV,        &      !Nitrogen limiation on SAV growth (non-dimensional)
            !PLSAV,        &      !Phosphorus limitation on SAV growth (non-dimensional)
            !FNSEDSAV,     &      !fraction of nitrogen uptake by SAV due to nitrogen in sediment pore water
                            !Wen Long: I think we should further split into NH4 and NO3 
                            !          uptake by SAV from sediments
            !FPSEDSAV,     &      !fraction of phosphorus uptake by SAV due to P in sediment pore water
            !NLEPI,        &      !Nitrogen limitation on epiphytes growth (non-dimensional)
            !PLEPI,        &      !Phosphorus limitation on epiphytes growth (non-dimensional)
            !FHS,        &   !fresh water SAV species growth limitation due to sulfide HS^- 
                            !(in sediment anerobic layer) toxicity (dimensionless)

            !EPATN,      &     !Light attenuation by epiphyte
                            !e^(-Kep*Acla*Adwcep*EP) in Cerco-Moore 2001 paper for Ish eqn
                            !e^(-Kep*Adwcep*EP) in new version as EP changed to gC-Epi/m^2-leaf unit
            !SAVATN,     &     !Light attenuation by SAV (shoots)
                            ![1-e^(-Ksh*SH+KESS*HCAN)]/(Ksh*SH+KESS*HCAN) in C-M 2001 Iwc eqn
            !WATATN,     &     !water column attenuation of light(dimensionless) at canopy top
                            != exp(-KESS*ZTOCNPY) where ZTOCNPY is distance from surface
                            !                           to canopy top
            !FISH,       &     !light limitation function for LEAF (shoots) f(I)=Ish/sqrt(Ish^2+Ik^2)
            !FIEP,       &     !light limitation function for epiphytes, f(I)=Iep/sqrt(Iep^2+Ik^2)
            !NPPSAV,     &     !net leaf production per day (growth-basal metabolism) (gC/m^2/day)
            !NPPEPI,    &   !net epiphyte production per day (growth - basal metabolism) (gC-epi/gC-leaf/day)

            !!fluxes to water column
            !DOSAVW,     &     !DOXG source (flux) due to SAV (source-sink) (gO2/m^2/day ) 
                            !DOSAVW --> DTDO (gO2/m^3/sec)--> DTC (gO2/m^3/sec)
                            !-->DTM(gO2/m^3/sec)-->DTM*DLT-->C1(:,:,27)(gO2/m^3) !27'th constituents is DOXG
            !LDOCSAVW,     &     !LDOC flux due to SAV (gC/m^2/day)
            
                            !WL         !***ALSO used by TVDS.F ****
                            !LDOCSAVW --> DTDOC (gC/m^3/sec) -> DTC --> DTM --> DTM*DLT --> C1(:,:,9)
            !RDOCSAVW,     &      !RDOC flux due to SAV (gC/m^2/day)                                 
            !LPOCSAVW,     &      !LPOC flux due to SAV (gC/m^2/day) 
                            !WL           !***ALSO used by TVDS.F ****
            !RPOCSAVW,     &      !RPOC flux due to SAV (gC/m^2/day)
                            !WL            !***ALSO used by TVDS.F ****
            !NH4SAVW,      &      !NH4 flux to water column (respiration source - uptake sink by SAV growth) (gN/m^2/day)
                            !NH4SAVW (gN/m^2/day) --> DTNH4  (gN/m^3/sec) --> DTC --> DTM --> C1(:,:,13)
            !NO3SAVW,      &      !NO3 flux to water column (due to uptake by SAV, positive increasing NO3 in water) (gN/m^2/day)
            !LDONSAVW,     &      !LDON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf, basal metabolism of stem) (gN/m^2/day)
            !RDONSAVW,     &      !RDON flux                                     
            !LPONSAVW,     &      !LPON flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem (gN/m^2/day)
            !RPONSAVW,     &      !RPON flux 
            !PO4SAVW,      &      !PO4 flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem and uptake by leaf) (gP/m^2/day)
            !LDOPSAVW,     &      !LDOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day)
            !RDOPSAVW,     &      !RDOP flux (gP/m^2/day)
            !LPOPSAVW,     &      !LPOP flux to water column (due to photochemical respiration of leaf, basal metabolism of leaf and stem) (gP/m^2/day) 
            !RPOPSAVW,     &      !RPOP flux (gP/m^2/day)
            !  DOEPIW,     &      !DOXG   flux (positive into water column) due to epiphytes (gO2/m^2/day)
            !LDOCEPIW,     &      !LDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
            !RDOCEPIW,     &      !RDOC flux (positive into water column) due to epiphytes (gC/m^2/day)
            !LPOCEPIW,     &      !LPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
            !RPOCEPIW,     &      !RPOC flux (positive into water column) due to epiphytes (gC/m^2/day)
            !NH4EPIW,     &      !NH4  flux (positive into water column) due to epiphytes (gN/m^2/day) (basal meabolism + predation + sloughing - uptake by photosynthesis)
            !NO3EPIW,     &      !NO3  flux
            !LDONEPIW,     &      !LDON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
            !RDONEPIW,     &      !RDON flux (positive into water column) (gN/m^2/day)
            !LPONEPIW,     &      !LPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
            !RPONEPIW,     &      !RPON flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gN/m^2/day)
            !PO4EPIW,      &      !PO4  flux (positive into water column) due to epiphytes mortality (basal metabolism, salt toxicity, sloughing) and uptake (sink) (growth) (gP/m^2/day)                                    
            !LDOPEPIW,     &      !LDOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
            !RDOPEPIW,     &      !RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
            !LPOPEPIW,     &      !LPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)
            !RPOPEPIW,     &    !RPOP flux (positive into water column) due to epiphytes mortality (basal metabolism, predation, and sloughing) (gP/m^2/day)

            !!fluxes to sediments                                                                         
            !SEDPOCSAV,     &      !POC flux to sediments (positive into sediments) due to SAV (gC/m^2/day)
                            !SEDPOCSAV (gC/m^2/day) --> JPOC (mgC/m^2/day) -> POC1, POC2, POC3 (source of particulate organic carbon in sediments)
            !SEDPONSAV,     &      !PON flux to sediments (positive into sediments) due to SAV (gN/m^2/day)
                            !SEDPONSAV (gN/m^2/day) --> JPON (mgN/m^2/day) -> PON1, PON2, PON3 (source of particulate organic nitrogen in sediments)
            !SEDPOPSAV,     &      !POP flux to sedimetns (positive into sediments) due to SAV (gP/m^2/day)
                            !SEDPOPSAV (gP/m^2/day) --> JPOP (mgP/m^2/day) -> POP1, POP2, POP3 (source of particulate organic phosphorus in sediments)
            !SEDDOSAV,      &      !DOXG flux to water column (positive into water) due to SAV growth in sediments (source) and mortality of SAV in sediments (tuber and root basal metabolism)(sink)
                            !(gO2/m^2/day) (WLong: does this have to be to water ?? I thought the sinks should be in sediments directly for root metabolism)
            !SEDNH4SAV,     &      !NH4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
                            !(gN/m^2/day), SEDNH4SAV --> NH4T2TM1 (mg/N/m^3) decerase in sediments 
            !SEDNO3SAV,     &      !NO3 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
                            !(gN/m^2/day), SEDNO3SAV --> NO3T2TM1 (mg/N/m^3) decerase in sediments 
            !SEDPO4SAV,     &   !PO4 flux to sediment (positive leaving sediment) due to SAV photosynthetic uptake of NH4 in sediments
            
            !!misc
            !DGRDAYS,       &      !Degree day calculated from temperature of water

            !!time average values 
                        
            !ALEAF,       &     !Average leaf biomass (gC/m^2)
            !AROOT,       &     !Average root biomass (gC/m^2)
            !ASTEM,       &     !Average stem biomass (gC/m^2)       
            !ATUBER,       &     !Average tuber biomass (gC/m^2)       
            !AEP,        &   !Average epiphyte biomass (gC-epi/gC-leaf)

            !APLEAF,       &   !Averaged leaf growth rate (1/day)
            !ABMLEAF,       &   !Averaged basal metabolism rate of leaf (1/day)
            !ABMTUBER,   &   !Average basal metabolism rate of tuber (1/day)

            !APEP,        &   !Average growth rate  for epiphytes (1/day)
            !ABMEP,       &   !Average basal metabolism rate for epiphytes (1/day)
            !APREP,        &    !Average predation rate on epiphytes (1/day) = APREP*EP

            !ASLSH,        &    !Average epiphute sloughing rate (1/day)
                
            !ANLSAV,        &   !Averaged N limitation on SAV 
            !APLSAV,        &   !Averaged P limitation on SAV
            !ANLEPI,        &   !Averaged N limitation on epiphytes
            !APLEPI,        &   !Averaged P limitation on epiphytes
            !AFNSED,        &   !Averaged fraction of uptake of N from sediment
            !AFPSED,        &   !Averaged fraction of uptake of P from sediment
            !AFHS,       &     !Averaged SAV growth limiation due to sulfide HS^-            
    
        !attenuation and light limitation
        !AEPATN,           &   !Averaged attenuation by epiphyte
        !AWATATN,           &   !Averaged water column attenuation of light at canopy top
        !AFISH,           &   !Averaged light limitation for LEAF (shoots)
        !AFIEP,        &    !Averaged light limitation for Epiphytes      

        !net production                                          
        !ANPPSAV,           &   !Averaged net leaf production per day (gC/m^2/day)
        !ANPPEPI,        &    !Averaged net epiphyte production per day (gC-epi/gC-leaf/day)

        !time averaged flux to water column
         !  ADOSAVW,      &      !average DOXG  flux due to SAV gO2/m^2/day
         ! ADOCSAVW,      &      !average DOC flux due to SAV
         ! APOCSAVW,      &      !average POC flux due to SAV
         ! ANH4SAVW,      &      !average NH4 flux due to SAV
         ! ANO3SAVW,      &      !average NO3 flux due to SAV
         ! ADONSAVW,      &      !average DON flux due to SAV              
         ! APONSAVW,      &      !average PON flux due to SAV
         ! APO4SAVW,      &      !average PO4 flux due to SAV
         ! ADOPSAVW,      &      !average DOP flux due to SAV
         ! APOPSAVW,      &      !average POL flux due to SAV
         !  ADOEPIW,      &      !DOXG  flux (positive into water column ) due to epiphytes (gO2/m2/day)
         ! ANH4EPIW,      &      !NH4 flux (positive into water column) due to epiphytes (gN/m2/day)
         ! ANO3EPIW,      &      !NO3 flux (positive into water column) due to epiphytes (gN/m2/day)
         ! APO4EPIW,      &      !PO4 flux (positive into water column) due to epiphytes (gP/m2/day)

         ! ADOCEPIW,      &      ! DOC flux to water column due to epiphytes (gC/m^2/day)
         ! APOCEPIW,      &      ! POC flux to water column due to epiphytes (gC/m^2/day)
         ! ADONEPIW,      &      ! DON flux to water column due to epiphytes (gN/m^2/day)
         ! APONEPIW,      &      ! PON flux to water column due to epiphytes (gN/m^2/day)                  
         ! ADOPEPIW,      &      ! DOP flux to water column due to epiphytes (gP/m^2/day)
         ! APOPEPIW,        &      ! POP flux to water column due to epiphytes (gP/m^2/day)

         ! ASEDDOSAV,       &   !DOXG  flux due to SAV,  positive leaving sediments
         ! ASEDPOCSAV,      &      !POC flux due to SAV, positive into sediments
         ! ASEDPONSAV,      &      !PON flux due to SAV, positive into sediments
         ! ASEDPOPSAV,      &      !POP flux due to SAV, positive into sediments
         ! ASEDNH4SAV,      &      !NH4 flux due to SAV, positive leaving sediemnts
         ! ASEDPO4SAV,        &    !PO4 flux due to SAV, positive leaving sediments
  
              !SAV_ALLOC,        &    !subroutine to allocate SAV variables
            SAV_DEALLOC,    &    !subroutine to deallocate SAV variables
        !SAV_LOADS_ALLOC,    &    !
        SAV_LOADS_DEALLOC,    &    !subroutine to deallocate SAV loads
            !SAV_READ,        &    !
            SAV_COMP            !subroutine that does SAV calculation
			
#if defined(SEDIMENT_DOM)			
			USE MOD_SED_DOM, ONLY : SED_DOM, &        ! B Clark DOM Flag
									SED_DOM_OUTPUT,JDOM_OUT_GL
									
#endif
									
#if defined(PHOTODEGRADATION)
	USE WC_DOM, ONLY : PHOTOCHEM_DOC_ACCUM, PHOTOCHEM_DOC_ACCUM_GL, PHOTOCHEM_DOC_ACCUM_OUT								
									
#endif
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


# if defined (SEDIMENT_DOM)

  USE MOD_SED_DOM_EXCHANGE, ONLY:  SED_DOM_FLAG,   &

                            JNH4_GL_SHARE,    &  !NH4 flux in sediments (mg N / m^2/day)
                            JNO3_GL_SHARE,    &  !NO3 flux in sediments (mg N / m^2/day)
                            NO31_GL_SHARE,&        !sediment layer 1 NO3 conc
                            NO32_GL_SHARE,&        !sediment layer 2 NO3 conc
                            NH41_GL_SHARE,&        !sediment layer 1 NH4 conc
                            NH42_GL_SHARE       !sediment layer 2 NH4 conc!,  &
# endif

   IMPLICIT NONE
# if defined (MULTIPROCESSOR)

   INCLUDE "mpif.h"

# endif

   INTEGER:: I,K

   CHARACTER(LEN=1024) :: TITLE_STN
   LOGICAL :: OPEN_HIS,PRINT_HIS
   LOGICAL :: OPEN_STN,PRINT_STN
   INTEGER(4) :: CTR_HIS, CTR_STN    !record counter for station and history outputs

   INTEGER(4) :: NTSTN_OLD, NTSTN
   INTEGER(4) :: NTHIS_OLD, NTHIS
 
!   CHARACTER(1024):: fnameprefix, fnameext
   CHARACTER(1024):: HISFN_FILENUMBER
   
   REAL(SP) :: JDAY_TMP            
   !REAL(SP) :: Jcin,Jnin,Jpin
   
   
!WLong: Output - field history and station timeseries
!                this should be moved into write_history.F and wrt_station.F later to
!                minimize the size of the main program and make it modular 


!WLong: Now starting and ending times of history and station outputs are controlled
!       in wqm_linkage.in file 

     NTHIS = FLOOR((JDAY-t_his_start)/t_his_dlt) 
     NTSTN = FLOOR((JDAY-t_stn_start)/t_stn_dlt)

     IF(((NTHIS>NTHIS_OLD).AND.(JDAY >= t_his_start)  &
                          .AND.(JDAY <= t_his_end))   &
       .OR.                                         &
        ((NTSTN>NTSTN_OLD).AND.(JDAY >= t_stn_start)  &
                          .AND.(JDAY <= t_stn_end))   &
       )THEN
	 

# if defined (MULTIPROCESSOR)
     
     
       if(PAR) then
        ! We are doing all K for some variables, K=1 for others
        call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,CCHL1(1,1),CCHL1_GL(1,1))
		
        call gather(0,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,total_netPP(0),total_netPP_GL(1))
		
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,27),C2_GL(1,1,27)) !dissolved oxygen

        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,11 ),C2_GL(1,1,11 )) ! LPOC
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,12 ),C2_GL(1,1,12 ))  ! RPOC
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,4 ),C2_GL(1,1,4 ))
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,5), C2_GL(1,1,5))
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,13),C2_GL(1,1,13))
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,14),C2_GL(1,1,14))
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,20),C2_GL(1,1,20))
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,T(0,1),T_GL(1,1))    !Wen Long, correction, pass all layers
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,SALT(0,1),S_GL(1,1)) !Wen Long, correction, pass all layers
        call gather(0,MTLOC,MLOC,MGL,1,   MYID,NPROCS,NMAP,D(0),D_GL(1))        !Wen Long, add D_GL for total depth
        call gather(0,MTLOC,MLOC,MGL,1,   MYID,NPROCS,NMAP,EL(0),EL_GL(1))      !Wen Long, add EL_GL for surface elevation 
        call gather(0,MTLOC,MLOC,MGL,1,   MYID,NPROCS,NMAP,H(0),H_GL(1))        !Wen Long, add H_GL for still water depth
     
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,9 ),C2_GL(1,1,9 )) ! CDOC1
        call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,10 ),C2_GL(1,1,10 )) !NCDOC1
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,33),C2_GL(1,1,33)) ! CDOC2
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,34),C2_GL(1,1,34)) !NCDOC2
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,35),C2_GL(1,1,35)) !CDOC3
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,36),C2_GL(1,1,36)) !NCDOC3
		
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,16),C2_GL(1,1,16)) !CDON1
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,17),C2_GL(1,1,17)) !NCDON1
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,37),C2_GL(1,1,37)) !CDON2
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,38),C2_GL(1,1,38)) !NCDON2
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,39),C2_GL(1,1,39)) !CDON3
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,40),C2_GL(1,1,40)) !NCDON3
																			
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,18),C2_GL(1,1,18)) !LPON
		call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,19),C2_GL(1,1,19)) !RPON
																			
		!call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,21),C2_GL(1,1,21)) !CDOP1
		!call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,22),C2_GL(1,1,22)) !NCDOP1
		!call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,41),C2_GL(1,1,41)) !CDOP2
		!call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,42),C2_GL(1,1,42)) !NCDOP2
		!call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,43),C2_GL(1,1,43)) !CDOP3
		!call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C2(0,1,44),C2_GL(1,1,44)) !NCDOP3
		
		
		
! Wen Long added the following for debugging algae

        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,P1(1,1)  ,   P1_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,P2(1,1)  ,   P2_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,BM1(1,1) ,  BM1_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,BM2(1,1) ,  BM2_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,NL1(1,1) ,  NL1_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,NL2(1,1) ,  NL2_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PL1(1,1) ,  PL1_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PL2(1,1) ,  PL2_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,FI1(1,1) ,  FI1_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,FI2(1,1) ,  FI2_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,B1SZ(1,1), B1SZ_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,B2SZ(1,1), B2SZ_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,B1LZ(1,1), B1LZ_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,B2LZ(1,1), B2LZ_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PR1(1,1) ,  PR1_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PR2(1,1) ,  PR2_GL(1,1))
        IF(LIGHT_EXTINCTION)THEN
                CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,IAVG(1,1), IAVG_GL(1,1))
                CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,KDLAMB_int(1,1), KDLAMB_int_GL(1,1))				
        ENDIF

		CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,ShearStress(1) ,  ShearStress_GL(1))		
		CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,RHO1(1,1) ,  RHO1_GL(1,1))
		
		CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,SSI(1,1),   SSI_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,IK1(1,1),   IK1_GL(1,1))
        CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,IK2(1,1),   IK2_GL(1,1))

! Wen Long added the following for debugging DOXG

       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DDOC(1,1)  ,    DDOC_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DCOD(1,1)  ,    DCOD_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,NITRIF(1,1),  NITRIF_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DOSZ(1,1)  ,    DOSZ_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DOLZ(1,1)  ,    DOLZ_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DOPR(1,1)  ,    DOPR_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DORALG(1,1),  DORALG_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,RESP(1,1)  ,    RESP_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,FTCOD(1,1) ,   FTCOD_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,COD(1,1)   ,     COD_GL(1,1))
       CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,REAERDO(1,1),REAERDO_GL(1,1))

!Wen Long added the following for debugging benthic flux including benthic SOD 
       CALL GATHER(1,MTLOC,MLOC,MGL,9   ,MYID,NPROCS,NMAP,BFLUX(1,1),BFLUX_GL(1,1))
	   
	   

	   
	   ! B Clark add in photochemical fluxes if we are using the module
#if defined(PHOTODEGRADATION) && defined(BUDGET)	
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,1),          PHOTOCHEM_DOC_ACCUM_GL(1,1,1))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,2),          PHOTOCHEM_DOC_ACCUM_GL(1,1,2))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,3),          PHOTOCHEM_DOC_ACCUM_GL(1,1,3))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,4),          PHOTOCHEM_DOC_ACCUM_GL(1,1,4))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,5),          PHOTOCHEM_DOC_ACCUM_GL(1,1,5))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,6),          PHOTOCHEM_DOC_ACCUM_GL(1,1,6))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,7),          PHOTOCHEM_DOC_ACCUM_GL(1,1,7))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,8),          PHOTOCHEM_DOC_ACCUM_GL(1,1,8))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,9),          PHOTOCHEM_DOC_ACCUM_GL(1,1,9))																													  
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,10),          PHOTOCHEM_DOC_ACCUM_GL(1,1,10))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,11),          PHOTOCHEM_DOC_ACCUM_GL(1,1,11))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,12),          PHOTOCHEM_DOC_ACCUM_GL(1,1,12))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,13),          PHOTOCHEM_DOC_ACCUM_GL(1,1,13))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,14),          PHOTOCHEM_DOC_ACCUM_GL(1,1,14))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,15),          PHOTOCHEM_DOC_ACCUM_GL(1,1,15))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,16),          PHOTOCHEM_DOC_ACCUM_GL(1,1,16))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,17),          PHOTOCHEM_DOC_ACCUM_GL(1,1,17))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM_OUT(1,1,18),          PHOTOCHEM_DOC_ACCUM_GL(1,1,18))		
#elif defined(PHOTODEGRADATION)
	   !This has the instantanous flux output
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,1),          PHOTOCHEM_DOC_ACCUM_GL(1,1,1))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,2),          PHOTOCHEM_DOC_ACCUM_GL(1,1,2))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,3),          PHOTOCHEM_DOC_ACCUM_GL(1,1,3))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,4),          PHOTOCHEM_DOC_ACCUM_GL(1,1,4))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,5),          PHOTOCHEM_DOC_ACCUM_GL(1,1,5))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,6),          PHOTOCHEM_DOC_ACCUM_GL(1,1,6))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,7),          PHOTOCHEM_DOC_ACCUM_GL(1,1,7))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,8),          PHOTOCHEM_DOC_ACCUM_GL(1,1,8))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,9),          PHOTOCHEM_DOC_ACCUM_GL(1,1,9))																													  
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,10),          PHOTOCHEM_DOC_ACCUM_GL(1,1,10))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,11),          PHOTOCHEM_DOC_ACCUM_GL(1,1,11))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,12),          PHOTOCHEM_DOC_ACCUM_GL(1,1,12))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,13),          PHOTOCHEM_DOC_ACCUM_GL(1,1,13))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,14),          PHOTOCHEM_DOC_ACCUM_GL(1,1,14))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,15),          PHOTOCHEM_DOC_ACCUM_GL(1,1,15))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,16),          PHOTOCHEM_DOC_ACCUM_GL(1,1,16))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,17),          PHOTOCHEM_DOC_ACCUM_GL(1,1,17))
	   CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,PHOTOCHEM_DOC_ACCUM(1,1,18),          PHOTOCHEM_DOC_ACCUM_GL(1,1,18))																													  
#endif
	   
																													  
! B CLark added below for budget terms, July 2018
		!!! Nitrogen	
#if defined (BUDGET)	   
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGDON_OUT(1,1), ALGDON_GL(1,1))																											  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGPON_OUT(1,1), ALGPON_GL(1,1))	
			
			!CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGPPDON(1,1), ALGPPDON_GL(1,1))	
			!CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGPredDON(1,1), ALGPredDON_GL(1,1))	
																					  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGNH4_OUT(1,1), ALGNH4_GL(1,1))	
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGNO3_OUT(1,1), ALGNO3_GL(1,1))	
																					  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DENNO3_OUT(1,1), DENNO3_GL(1,1))																						  		  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,NT_OUT(1,1), NT_GL(1,1))																						  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DENITN_OUT(1,1), DENITN_GL(1,1))		
																					  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,HDRLPON_OUT(1,1), HDRLPON_GL(1,1))																				  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,HDRRPON_OUT(1,1), HDRRPON_GL(1,1))
																						
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,COAGN_OUT(1,1), COAGN_GL(1,1))																							
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,COAGC_OUT(1,1), COAGC_GL(1,1))
																					
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,MNLDON1_OUT(1,1), MNLDON1_GL(1,1))																							
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,MNLDON2_OUT(1,1), MNLDON2_GL(1,1))																					
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,MNLDON3_OUT(1,1), MNLDON3_GL(1,1))				
																					
	!!!! CARBON																				
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,MNLDOC1_OUT(1,1), MNLDOC1_GL(1,1))																							
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,MNLDOC2_OUT(1,1), MNLDOC2_GL(1,1))																					
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,MNLDOC3_OUT(1,1), MNLDOC3_GL(1,1))
																						
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGDOC_OUT(1,1), ALGDOC_GL(1,1))																											  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGPOC_OUT(1,1), ALGPOC_GL(1,1))	
			
			!CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGPPDOC(1,1), ALGPPDOC_GL(1,1))	
			!CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGPredDOC(1,1), ALGPredDOC_GL(1,1))	
			
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DENIT_OUT(1,1), DENIT_GL(1,1))
			
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,HDRLPOC_OUT(1,1), HDRLPOC_GL(1,1))																				  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,HDRRPOC_OUT(1,1), HDRRPOC_GL(1,1))	
																						
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDOC1_out(1,1),DTWCDOC1_GL(1,1)) ! CDOC1
        call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDOC2_out(1,1),DTWCDOC2_GL(1,1)) !NCDOC1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDOC3_out(1,1),DTWCDOC3_GL(1,1)) ! CDOC2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDOC1_out(1,1),DTWNCDOC1_GL(1,1)) !NCDOC2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDOC2_out(1,1),DTWNCDOC2_GL(1,1)) !CDOC3
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDOC3_out(1,1),DTWNCDOC3_GL(1,1)) !NCDOC3
		
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDON1_out(1,1),DTWCDON1_GL(1,1))!CDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDON2_out(1,1),DTWCDON2_GL(1,1)) !NCDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDON3_out(1,1),DTWCDON3_GL(1,1))!CDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDON1_out(1,1),DTWNCDON1_GL(1,1)) !NCDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDON2_out(1,1),DTWNCDON2_GL(1,1)) !CDON3
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDON3_out(1,1),DTWNCDON3_GL(1,1))!NCDON3
		!																			
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTNH4_out,DTNH4_GL(1,1))!CDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTNO3_out,DTNO3_GL(1,1)) !NCDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPOC_out,DTLPOC_GL(1,1))!CDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTRPOC_out,DTRPOC_GL(1,1)) !NCDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPON_out,DTLPON_GL(1,1)) !CDON3
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPON_out,DTRPON_GL(1,1))!NCDON3
		!																				
#else
	    call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDOC1(1,1),DTWCDOC1_GL(1,1)) ! CDOC1
        call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDOC2(1,1),DTWCDOC2_GL(1,1)) !NCDOC1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDOC3(1,1),DTWCDOC3_GL(1,1)) ! CDOC2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDOC1(1,1),DTWNCDOC1_GL(1,1)) !NCDOC2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDOC2(1,1),DTWNCDOC2_GL(1,1)) !CDOC3
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDOC3(1,1),DTWNCDOC3_GL(1,1)) !NCDOC3
		
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDON1(1,1),DTWCDON1_GL(1,1))!CDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDON2(1,1),DTWCDON2_GL(1,1)) !NCDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWCDON3(1,1),DTWCDON3_GL(1,1))!CDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDON1(1,1),DTWNCDON1_GL(1,1)) !NCDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDON2(1,1),DTWNCDON2_GL(1,1)) !CDON3
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTWNCDON3(1,1),DTWNCDON3_GL(1,1))!NCDON3
		!																			
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTNH4,DTNH4_GL(1,1))!CDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTNO3,DTNO3_GL(1,1)) !NCDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPOC,DTLPOC_GL(1,1))!CDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTRPOC,DTRPOC_GL(1,1)) !NCDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPON,DTLPON_GL(1,1)) !CDON3
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPON,DTRPON_GL(1,1))!NCDON3	
		
#endif
		
!WLong: Note above, C2 and C2_GL all have first indices starting from zero
!    I wonder why would we pass C2(0,,) to C2_GL(1,,) instead of C2_GL(0,*,*) 
!    but SALT starts from 0, and S_GL starts from 1
!    There are might be problems here to sort out 


!WLong: We ought to gather all other components of C2 to C2_GL as well
!       Right now, above statements only deals with afew state variables :27,9,10,11,12,4,5,13,14,20)
!       To fully parallelize ICM, the rest of them must also be gathered and 
!       output to RSO and ICO for restart purpose
!       As well as output to profiles and plots etc
    
       
       IF(SEDIMENT_CALC)THEN !LBcomment
       
        !   
        !gather sediment diagenesis related variables
        !
       
           !all of the following are used for outputing sediment flux variables into history and station files
        !Make sure all of them are collected
    
            !(JPOC_GL(MGL,3));
            !(JPON_GL(MGL,3));
            !(JPOP_GL(MGL,3));
            !(JPOS_GL(MGL));    
#if defined(BUDGET)
	        CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,JPOC_out(1,1)  ,    JPOC_GL(1,1))
            CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,JPON_out(1,1)  ,    JPON_GL(1,1))   
			CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JNH4_out(1)  ,        JNH4_GL(1))               
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JNO3_out(1)  ,        JNO3_GL(1))   
	   
#else
            CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,JPOC(1,1)  ,    JPOC_GL(1,1))
            CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,JPON(1,1)  ,    JPON_GL(1,1))
            CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,JPOP(1,1)  ,    JPOP_GL(1,1))
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JPOS(1)  ,        JPOS_GL(1))
			CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JNH4TM1S(1)  ,        JNH4_GL(1))
            !(JNO3_GL(MGL));                
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JNO3TM1S(1)  ,        JNO3_GL(1))   
#endif       
            !!(O20_GL(MGL));            !given by C2_GL(:,:,27)
            !!(D_GL(MGL));                !already defined
            !(CTEMP_GL(MGL));             !
        
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CTEMP(1)  ,        CTEMP_GL(1))

            
            CH40_GL=0.d0                 !temporarily set to zero
            !!(SAL_GL(MGL));            !given by S_GL
            !(SODTM1S_GL(MGL));      
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SODTM1S(1)  ,        SODTM1S_GL(1))
            !!(SOD_GL(MGL));            !given by SODTM1S_GL
            !(JNH4_GL(MGL));                     
            !(BENDEN_GL(MGL));        
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,BENDEN(1)  ,        BENDEN_GL(1))        
            !(JCH4_GL(MGL));                
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JCH4TM1S(1)  ,        JCH4_GL(1))        
            !(JCH4G_GL(MGL));        
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JCH4GTM1S(1)  ,        JCH4G_GL(1))        
            !(JHS_GL(MGL));                
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JHSTM1S(1)  ,        JHS_GL(1))        
            !(JPO4_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JPO4TM1S(1)  ,        JPO4_GL(1))        
            !(JSI_GL(MGL));                
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JSITM1S(1)  ,        JSI_GL(1))        
            !(NH41_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NH41TM1S(1)  ,        NH41_GL(1))        
            !(NH42_GL(MGL));        
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NH42TM1S(1)  ,        NH42_GL(1))        
            !(NO31_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NO31TM1S(1)  ,        NO31_GL(1))        
            !(NO32_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NO32TM1S(1)  ,        NO32_GL(1))    
            !(PO41_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,PO41TM1S(1)  ,        PO41_GL(1))        
            !(PO42_GL(MGL));        
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,PO42TM1S(1)  ,        PO42_GL(1))        
            !(SI1_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SI1TM1S(1)  ,        SI1_GL(1))        
            !(SI2_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SI2TM1S(1)  ,        SI2_GL(1))        
            !(CH41TM1S_GL(MGL));
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CH41TM1S(1)  ,        CH41TM1S_GL(1))        
            !!(CH41_GL(MGL));            !given by CH41TM1S_GL
            !(CH42_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CH42TM1S(1)  ,        CH42_GL(1))    
            !(HS1_GL(MGL));            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,HS1TM1S(1)  ,        HS1_GL(1))    
            !(HS2_GL(MGL));                                        
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,HS2TM1S(1)  ,        HS2_GL(1))
            !(CPOP_GL(MGL,3));            
            CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPOP(1,1)  ,    CPOP_GL(1,1))
            !(CPON_GL(MGL,3));            
            CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPON(1,1)  ,    CPON_GL(1,1))
            !(CPOC_GL(MGL,3));            
            CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPOC(1,1)  ,    CPOC_GL(1,1))
            !(CPOS_GL(MGL));          
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CPOS(1)  ,        CPOS_GL(1))

            !!(POC1_GL(MGL));            !given by CPOC_GL
            !!(POC2_GL(MGL));            
            !!(POC3_GL(MGL));            
        
            !!(PON1_GL(MGL));            !given by CPON_GL
            !!(PON2_GL(MGL));            
            !!(PON3_GL(MGL));            

            !!(POP1_GL(MGL));            !given by CPOP_GL
            !!(POP2_GL(MGL));            
            !!(POP3_GL(MGL));            
            
            !!(PSISED_GL(MGL));            !given by CPOS_GL
        
            !(HSED1_GL(MGL));            !
            
            
            
            !-----
            !===================================================================================|
   !SUBROUTINE GATHER(N1,N2,NT,NTG,KT,MYID,NPROCS,GM,A,AG) 
   !INTEGER,   INTENT(IN)           :: N1,N2,NT,NTG,KT,MYID,NPROCS
   !TYPE(GMAP), INTENT(IN)           :: GM(NPROCS)
   !REAL(SP),    INTENT(IN)           :: A(N1:N2,KT)  !source array bounds
   !REAL(SP),    INTENT(OUT)          :: AG(NTG,KT)
  !            -----
            
            !first 1 (N1) is starting index of local array HSED1
            !MTLOC is end index of loal array HSED1
            !NT is not used yet
            !1:MGL is the dimension of the global array HSED1_GL
            !second 1 (KT) (before MYID) is the vertical dimension (last dimension) of local array and global array
            !HSED1(1) is the starting memory address of local array to communicate
            !HSED1_GL is the starting memory address of global array to communicate
            
            
            !trying to map A(N1:N2,KT) to global array AG(1:NTG, KT)
            
            !instead of passing A and AG, it only passes the referene (starting address of the array)
            !in theory, 
            !----
            
            
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,HSED1(1)  ,  HSED1_GL(1))        
            !(BENSTRTM1S_GL(MGL));      
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,BENSTRTM1S(1)  ,        BENSTRTM1S_GL(1))        
            !!(BENSTR_GL(MGL));            !given by BENSTRTM1S_GL        
            
            !
            !these are the outputs of the excel version (Greg Pelletier) 
            !                                
            !        Print Quantity                        Variable            FortranUnit            ExceUnit
                !!',',(JPOC(I,1)+JPOC(I,2)+JPOC(I,3))/1000.0*2.667,    &    
                                                        !'Jcin'                mgC/m^2/d            gO2/m^2/d
                !!',',(JPON(I,1)+JPON(I,2)+JPON(I,3))/1000.0,        &
                                                        !'Jnin'                *mgN/m^2/d            gN/m^2/d
                !!',',(JPOP(I,1)+JPOP(I,2)+JPOP(I,3))/1000.0,        &
                                                        !'Jpin'                *mgP/m^2/d            gP/m^2/d
                !!',',JPOS(I)/1000.0,&                    !'Jsin                 *mgSi/m^2/d            gSi/m^2/d
                !!',',O20,            &                    !'O20'                mgO2/L                mgO2/L
                !!',',D(I),            &                    !'Depth'            m                    m
                !!',',T(I,KBM1),    &                    !'Tw'                degC                degC
                !!',',NH40/1000.0,    &                    !'NH30'                *mgN/m^3            mgN/L
                !!',',NO30/1000.0,    &                    !'NO30'                *mgN/m^3            mgN/L
                !!',',PO40/1000.0,    &                    !'PO40'                *mgP/m^3            mgP/L
                !!',',SI0/1000.0,    &                    !'SI0'                *mgSi/m^3             mgSi/L        (dissiolved SIAT)
                !!',',CH40,            &                    !'CH40'                gO2/m^3                mgO2/L
                !!',',SAL,            &                    !'SALw'              ppt                    ppt            !i.e. SALT(I,KWC), ppt
                !!',',SOD,            &                    !'SOD'                gO2/m^2/day            gO2/m^2/d
                !!',',JNH4/1000.0,    &                    !'Jnh4                *mgN/m^2/day        gN/m^2/d
                !!',',JNO3/1000.0,    &                    !'Jno3                *mgN/m^2/day        gN/m^2/d
                !!',',BENDEN(1),    &                    !'JDenitT'            gN/m^2/day            gN/m^2/d
                !!',',JCH4,            &                    !'Jch4'                gO2/m^2/day            gO2/m^2/d
                !!',',JCH4G,        &                    !'Jch4g'            gO2/m^2/day            gO2/m^2/d
                !!',',JHS,            &                    !'Jhs'                gO2/m^2/day            gO2/m^2/d
                !!',',JPO4/1000.0,    &                    !'Jpo4'                *mgP/m^2/day        gP/m^2/d
                !!',',JSI/1000.0,    &                    !'Jsi'                *mgSi/m^2/day        gSi/m^2/d
                !!',',NH41/1000.0,    &                    !'NH3(1)'            *mgN/m^3             mgN/L        (dissolved)    
                !!',',NH42/1000.0,    &                    !'NH3(2)'            *mgN/m^3             mgN/L        (dissolved)
                !!',',NO31/1000.0,    &                    !'NO3(1)'            *mgN/m^3             mgN/L        (dissolved)
                !!',',NO32/1000.0,    &                    !'NO3(2)'            *mgN/m^3             mgN/L        (dissolved)
                !!',',PO41/1000.0,    &                    !'PO4(1)'            *mgP/m^3             mgP/L        (dissolved)
                !!',',PO42/1000.0,    &                    !'PO4(2)'            *mgP/m^3             mgP/L        (dissolved)
                !!',',SI1/1000.0,    &                    !'Si(1)'            *mgSi/m^3             mgSi/L        (dissolved)
                !!',',SI2/1000.0,    &                    !'Si(2)'            *mgSi/m^3             mgSi/L        (dissolved)
                !!',',CH41,            &                    !'CH4(1)'            gO2/m^3             mgO2/L        (dissolved)
                !!',',CH42,            &                    !'CH4(2)'            gO2/m^3             mgO2/L        (dissolved)
                !!',',HS1,            &                    !'HS(1)'            gO2/m^3             mgO2/L        (dissolved)
                !!',',HS2,            &                    !'HS(2)'            gO2/m^3             mgO2/L        (dissolved)
                !!',',POC1/1000.0*2.667,&                !'POC2(1)'            *mgC/m^3             gO2/m^3        (G1 of POC in layer 2)
                !!',',POC2/1000.0*2.667,&                !'POC2(2)'             *mgC/m^3             gO2/m^3        (G2 of POC in layer 2)
                !!',',POC3/1000.0*2.667,&                !'POC2(3)'             *mgC/m^3             gO2/m^3        (G3 of POC in layer 2)
                !!',',PON1/1000.0,    &                    !'PON2(1)'            *mgN/m^3             gN/m^3        (G1 of PON in layer 2)
                !!',',PON2/1000.0,    &                    !'PON2(2)'             *mgN/m^3             gN/m^3        (G2 of PON in layer 2)
                !!',',PON3/1000.0,    &                    !'PON2(3)'             *mgN/m^3             gN/m^3        (G3 of PON in layer 2)
                !!',',POP1/1000.0,    &                    !'POP2(1)'            *mgP/m^3             gP/m^3        (G1 of POP in layer 2)
                !!',',POP2/1000.0,    &                    !'POP2(2)'            *mgP/m^3             gP/m^3        (G2 of POP in layer 2)
                !!',',POP3/1000.0,    &                    !'POP2(3)'             *mgP/m^3             gP/m^3        (G3 of POP in layer 2)
                !!',',PSISED/1000.0,&                    !'POS2'                *mgSi/m^3             gSi/m^3        (particulate biogenic Silicate in layer2)
                !!',',HSED1(I)*1000.0,&                    !'H1'                m                    mm
                !!',',BENSTR                             !'BEN_STR'             dimensionless        dimensionless



# if defined(SEDIMENT_DOM)

                IF(SED_DOM) THEN
                    JNH4_GL_SHARE = JNH4_GL!NH4 flux in sediments (mg N / m^2/day)
                    JNO3_GL_SHARE = JNO3_GL!NO3 flux in sediments (mg N / m^2/day)
                    NO31_GL_SHARE = NO31_GL   !sediment layer 1 NO3 conc
                    NO32_GL_SHARE = NO32_GL      !sediment layer 2 NO3 conc
                    NH41_GL_SHARE = NH41_GL      !sediment layer 1 NH4 conc
                    NH42_GL_SHARE = NH42_GL	
					CALL SED_DOM_OUTPUT

                ENDIF
# endif


        ENDIF !!end LBcomment
       else
# endif
        CCHL1_GL(1:MGL,1:KBM1) = CCHL1(1:MGL,1:KBM1)
       total_netPP_GL(1:MGL)= total_netPP(1:MGL)
        C2_GL(1:MGL,1:KBM1,27) = C2(1:MGL,1:KBM1,27)
	  C2_GL(1:MGL,1:KBM1,26) = C2(1:MGL,1:KBM1,26)
      !  C2_GL(1:MGL,1:KBM1,9 ) = C2(1:MGL,1:KBM1,9 )
     !   C2_GL(1:MGL,1:KBM1,10 ) = C2(1:MGL,1:KBM1,10 )
        C2_GL(1:MGL,1:KBM1,11 ) = C2(1:MGL,1:KBM1,11 )
        C2_GL(1:MGL,1:KBM1,12 ) = C2(1:MGL,1:KBM1,12 )
        C2_GL(1:MGL,1:KBM1,4 ) = C2(1:MGL,1:KBM1,4 )
        C2_GL(1:MGL,1:KBM1,5)  = C2(1:MGL,1:KBM1,5)
        C2_GL(1:MGL,1:KBM1,13) = C2(1:MGL,1:KBM1,13)
        C2_GL(1:MGL,1:KBM1,14) = C2(1:MGL,1:KBM1,14)
        C2_GL(1:MGL,1:KBM1,20) = C2(1:MGL,1:KBM1,20)
   
		C2_GL(:,:,9) = C2(:,:,9)
		C2_GL(:,:,10) = C2(:,:,10)
		C2_GL(:,:,33) = C2(:,:,33)
		C2_GL(:,:,34) = C2(:,:,34)
		C2_GL(:,:,35) = C2(:,:,35)
		C2_GL(:,:,36) = C2(:,:,36)
		
		C2_GL(:,:,16) = C2(:,:,16)
		C2_GL(:,:,17) = C2(:,:,17)
		C2_GL(:,:,37) = C2(:,:,37)
		C2_GL(:,:,38) = C2(:,:,38)
		C2_GL(:,:,39) = C2(:,:,39)
		C2_GL(:,:,40) = C2(:,:,40)
		
		C2_GL(:,:,21) = C2(:,:,21)
		C2_GL(:,:,22) = C2(:,:,22)
		C2_GL(:,:,41) = C2(:,:,41)
		C2_GL(:,:,42) = C2(:,:,42)
		C2_GL(:,:,43) = C2(:,:,43)
		C2_GL(:,:,44) = C2(:,:,44)
	
        T_GL(1:MGL,1:KBM1) = T(1:MGL,1:KBM1)            ! Wen Long, correction to assign values for all layers
        S_GL(1:MGL,1:KBM1) = SALT(1:MGL,1:KBM1)         ! Wen Long, correction to assign values for all layers

        D_GL(1:MGL)= D(1:MGL)                           ! Wen Long, add D_GL, EL_GL, H_GL
        EL_GL(1:MGL)= EL(1:MGL)
        H_GL(1:MGL) = H(1:MGL)

! Wen Long added following for debugging algae

        P1_GL   = P1
        P2_GL   = P2
        BM1_GL  = BM1
        BM2_GL  = BM2
        NL1_GL  = NL1
        NL2_GL  = NL2
        PL1_GL  = PL1
        PL2_GL  = PL2
        FI1_GL  = FI1
        FI2_GL  = FI2
        B1SZ_GL = B1SZ
        B2SZ_GL = B2SZ
        B1LZ_GL = B1LZ
        B2LZ_GL = B2LZ
        PR1_GL  = PR1
        PR2_GL  = PR2
        IF(LIGHT_EXTINCTION)THEN
            IAVG_GL = IAVG
			KDLAMB_int_GL=KDLAMB_int
        ENDIF
        IK1_GL  = IK1
        IK2_GL  = IK2

! Wen Long added the following for debugging DOXG

        DDOC_GL   = DDOC
        DCOD_GL   = DCOD
        NITRIF_GL = NITRIF
        DOSZ_GL   = DOSZ
        DOLZ_GL   = DOLZ
        DOPR_GL   = DOPR
        DORALG_GL = DORALG
        RESP_GL   = RESP
        FTCOD_GL  = FTCOD
        COD_GL    = COD
        REAERDO_GL = REAERDO

! Wen Long added the following for debuguing benthic fluxes
        BFLUX_GL = BFLUX

        !collect sediment outputs to global arrays
        IF(SEDIMENT_CALC)THEN
            JPOC_GL=JPOC
            JPON_GL=JPON
            JPOP_GL=JPOP
            JPOS_GL=JPOS
            CTEMP_GL=CTEMP
            CH40_GL=0.d0                    !temporarily set to zero for no CH4 in water column model !0.d0 instead of 0.0 to make it equal to serial run- 
            SODTM1S_GL=SODTM1S
            JNH4_GL=JNH4TM1S
            JNO3_GL=JNO3TM1S
            BENDEN_GL=BENDEN
            JCH4_GL=JCH4TM1S
            JCH4G_GL=JCH4GTM1S
            JHS_GL=JHSTM1S
            JPO4_GL=JPO4TM1S
            JSI_GL=JSITM1S
            NH41_GL=NH41TM1S
            NH42_GL=NH42TM1S
            NO31_GL=NO31TM1S
            NO32_GL=NO32TM1S
            PO41_GL=PO41TM1S
            PO42_GL=PO42TM1S
            SI1_GL=SI1TM1S
            SI2_GL=SI2TM1S
            CH41TM1S_GL=CH41TM1S
            CH42_GL=CH42TM1S
            HS1_GL=HS1TM1S
            HS2_GL=HS2TM1S
            CPOP_GL=CPOP
            CPON_GL=CPON
            CPOC_GL=CPOC
            CPOS_GL=CPOS
            HSED1_GL=HSED1
            BENSTRTM1S_GL=BENSTRTM1S

!# if defined(SEDIMENT_DOM)

!           IF(SED_DOM) THEN

!                JNH4_GL_SHARE = JNH4_GL!NH4 flux in sediments (mg N / m^2/day)
!                JNO3_GL_SHARE = JNO3_GL!NO3 flux in sediments (mg N / m^2/day)
!                NO31_GL_SHARE = NO31_GL   !sediment layer 1 NO3 conc
!                NO32_GL_SHARE = NO32_GL      !sediment layer 2 NO3 conc
!                NH41_GL_SHARE = NH41_GL      !sediment layer 1 NH4 conc
!                NH42_GL_SHARE = NH42_GL

!		    !	CALL SED_DOM_OUTPUT

!				
!           ENDIF
!	   
!# endif
        ENDIF

# if defined (MULTIPROCESSOR)
       endif
# endif

	   
!WLong: Instead of outputing different layers to different files, here, I write all history records
!       into the same file HISFN, where HISFN is input in wqm_linkage.in

       PRINT_HIS = .FALSE. 

!      IF(AMOD(JDAY-t_his_start,t_his_dlt).eq.0)THEN
       IF(NTHIS > NTHIS_OLD) THEN
          CTR_HIS = CTR_HIS + 1
          PRINT_HIS = .TRUE. 
          NTHIS_OLD = NTHIS
       ENDIF

       IF (MSR) THEN
         WRITE(*,*)'CTR_HIS=', CTR_HIS
         IF(CTR_HIS.eq.1) THEN  !Open file if the first record

             SELECT CASE (HISFN_SPLIT_BYLEVEL)

                 CASE(.TRUE.)  !open one file for each vertical level
                   !find file names based on layer number

                    HISFN_PREFIX=TRIM(fnameprefix(HISFN))  !obtain prefix (removing extension)
                    HISFN_EXT   =TRIM(fnameext(HISFN))    !obtain extension if any

                    DO K=1,KBM1
                       WRITE(HISFN_FILENUMBER,'(I5.5)')K
                       HISFN_FINAL=TRIM(HISFN_PREFIX)//'_'//TRIM(HISFN_FILENUMBER)//TRIM(HISFN_EXT)
                       WRITE(*,*)'Opening history output file:',TRIM(HISFN_FINAL)
                       OPEN(UNIT=UNIT_HIS+K,FILE=TRIM(HISFN_FINAL),STATUS='UNKNOWN')
                    ENDDO

                 CASE(.FALSE.) !all levels in the same file (open just one file)

                    OPEN(UNIT=UNIT_HIS,FILE=TRIM(HISFN),STATUS='UNKNOWN')

              ENDSELECT

              OPEN_HIS = .TRUE.     !flag that history file(s) is(are) open

         ENDIF

        !make sure we output a prior step if steady state is called for setting up diagenesis initial condition
         IF(STEADY_STATE_SED_IC.AND.(JDAY<TMSTRT+DLT/86400.0))THEN
            JDAY_TMP=JDAY-DLT    
         ELSE
            JDAY_TMP=JDAY
         ENDIF
                  
         IF(PRINT_HIS)THEN   !output history of fields

           DO K=1,KBM1
             WRITE(*,*)'Writing to history outputs, JDAY= ',JDAY_TMP,' K= ',K

             SELECT CASE (HISFN_SPLIT_BYLEVEL)
                  CASE(.TRUE.)
                  
                  
                      WRITE(UNIT_HIS+K,'(:F8.4,1x,I1,1x,I8/(:8(E12.6,1x),E12.6))')JDAY_TMP,4,MGL,  &
                                     (H_GL(I),I=1,MGL),            &  !Wen Long output still water depth                    1
                                     (EL_GL(I),I=1,MGL),           &  !Wen Long output surface elevation					2
                                     (D_GL(I),I=1,MGL),            &  !Wen Long output total depth							3	
                                     (-D_GL(I)*ZZ(K),I=1,MGL),     &  !Wen Long,depth of this layer (from surface)			4									 
! B Clark deleted state variable output above, now all state variables can be output in netcdf format
!!Below are the volume and time integrated budget terms, these are the only values output, preferably at t0 and tfinal
!use netcdf for state variable time series, these will be for the budget
! Changed Feb, 2019									 
	!=============== Sediment water column Flux ============================
                        (JNH4_GL(I),I=1,MGL),          &   !														5
                        (JNO3_GL(I),I=1,MGL),          &					!										6										
						
						(JPOC_GL(I,1),I=1,MGL),        &      !'Jcin1'                 								7		      
						(JPOC_GL(I,2), I=1,MGL),       &      !'Jcin2'                  							8		      
						(JPOC_GL(I,3),I=1,MGL),        &      !'Jcin3'               							    9
						
					    (JPON_GL(I,1),I=1,MGL),        &      !'Jnin1'                      						10
						(JPON_GL(I,2),I=1,MGL),        &      !'Jnin2'                       						11
						(JPON_GL(I,3),I=1,MGL),        &      !'Jnin3'             									12
							
						(JDOM_OUT_GL(I,1),I=1,MGL),    &      ! JCDOC1  											13
					    (JDOM_OUT_GL(I,4),I=1,MGL),    &  !	JNCDOC1													14	
						 
						(JDOM_OUT_GL(I,2),I=1,MGL),    & !    JCDOC2  												15
					    (JDOM_OUT_GL(I,5),I=1,MGL),    & ! JNCDOC2													16
						 
						(JDOM_OUT_GL(I,3) ,I=1,MGL),   &  !    JCDOC3  												17
					    (JDOM_OUT_GL(I,6),I=1,MGL),    & ! 	JNCDOC3													18

						(JDOM_OUT_GL(I,7),I=1,MGL),    &   !    JCDON1  											19
					    (JDOM_OUT_GL(I,10),I=1,MGL),   &  !	JNCDON1													20
						 
						(JDOM_OUT_GL(I,8) ,I=1,MGL),   & !      JCDON2  											21	
					    (JDOM_OUT_GL(I,11),I=1,MGL),   & ! 	JNCDON2													22
						
						(JDOM_OUT_GL(I,9),I=1,MGL),    &!      JCDON3  												23
					    (JDOM_OUT_GL(I,12),I=1,MGL),   & ! 	JNCDON3													24
						 
! ==================== Water column changes ============================
				!!!!! NITROGEN		 
						  (ALGDON_GL(I,K),I=1,MGL),   &   ! Algae DON production (exudation and predation) 			25
						  (ALGPON_GL(I,K),I=1,MGL),   &   ! Algae PON production (exudation and predation) 			26
              (ALGNH4_GL(I,K),I=1,MGL),   &   ! Algae NH4 balance (uptake/excretion/predation) 			27
              
						  (ALGNO3_GL(I,K),I=1,MGL),   &   ! Algae NO3 Uptake 										28
						  (DENNO3_GL(I,K),I=1,MGL),   &   ! Denitrification loss of NO3								29		
              (NT_GL(I,K),I=1,MGL),     & ! Loss of Nh4 due to nitrification, production of NO3 		30
            
						  (DENITN_GL(I,K),I=1,MGL),  & !remin of DON due to denitrification  (makes NH4)       		31
						  (HDRLPON_GL(I,K),I=1,MGL), & ! hydrolytic loss of LPON to DON                     		32
              (HDRRPON_GL(I,K),I=1,MGL), & ! hydrolytic loss of RPON to DON                        		33
              
              (COAGN_GL(I,K),I=1,MGL),  & ! coagulation of DON3 to RPON 								34
						  (MNLDON1_GL(I,K),I=1,MGL), & ! remin of DON 1 											35
              (MNLDON2_GL(I,K),I=1,MGL), & ! remin of DON 2 											36			
              
						  (MNLDON3_GL(I,K),I=1,MGL), & ! remin of DON 3 											37
						  (C2_GL(I,K,18),I=1,MGL),   & ! LPON 														38
              (C2_GL(I,K,19),I=1,MGL),   & ! RPON 														39
              
						  (DTWCDON1_GL(I,K),I=1,MGL),&		!														40
						  (DTWCDON2_GL(I,K),I=1,MGL),&!																41
              (DTWCDON3_GL(I,K),I=1,MGL), &	!														    42
              
						  (DTWNCDON1_GL(I,K),I=1,MGL),&!															43
						  (DTWNCDON2_GL(I,K),I=1,MGL),&!															44	
              (DTWNCDON3_GL(I,K),I=1,MGL),&	! 															45
              
						  (DTNH4_GL(I,K),I=1,MGL),  &	!															46
              (DTNO3_GL(I,K),I=1,MGL),  & 	!								!							47
              
						  (DTLPON_GL(I,K),I=1,MGL),    &!															48
						  (DTRPON_GL(I,K),I=1,MGL),  & !															49
 
				!!!Carbon	
						  (ALGDOC_GL(I,K),I=1,MGL),   &   ! Algae DOC production (exudation and predation) 			50
						  (ALGPOC_GL(I,K),I=1,MGL),   &   ! Algae POC production (exudation and predation) 			51	
              (DENIT_GL(I,K),I=1,MGL),   & ! denitrification remin of DOC 								52
              
						  (HDRLPOC_GL(I,K),I=1,MGL),  & ! Hydrolysis of LPOC to DOC 								53
						  (HDRRPOC_GL(I,K),I=1,MGL),  & ! Hydrolysis of RPOC to DOC 								54	  
              (COAGC_GL(I,K),I=1,MGL),   & ! coagulation of DOC3 to RPOC 								55	
              					  
						  (MNLDOC1_GL(I,K),I=1,MGL), & ! remin of DOC 1 											56
						  (MNLDOC2_GL(I,K),I=1,MGL), & ! remin of DOC 2 											57				
              (MNLDOC3_GL(I,K),I=1,MGL), & ! remin of DOC 3 											58
              
						  (DTWCDOC1_GL(I,K), I=1,MGL), & !															59
						  (DTWCDOC2_GL(I,K), I=1,MGL), &!															60
              (DTWCDOC3_GL(I,K), I=1,MGL), &!															61
            
						  (DTWNCDOC1_GL(I,K), I=1,MGL), &!															62
						  (DTWNCDOC2_GL(I,K), I=1,MGL), &!															63
						  (DTWNCDOC3_GL(I,K), I=1,MGL), & !															64	
						  
						  (DTLPOC_GL(I,K),I=1,MGL) ,   & !														65
						  (DTRPOC_GL(I,K),I=1,MGL),   & 								!							66
		 
! ==================== Photochemistry ==================================						 
						 
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,1),I=1,MGL),&! Phototransformation from CDOC3 --> CDOC2 		67
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,2),I=1,MGL),&! Phototransformation from CDOC3 --> CDOC1 		68	
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,3),I=1,MGL),&! Phototransformation from CDOC3 --> DIC 		69
              (PHOTOCHEM_DOC_ACCUM_GL(I,K,4),I=1,MGL),&! Phototransformation from CDOC3 --> NCDOC 		70
              				 
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,5),I=1,MGL),&! Phototransformation from CDOC2 --> CDOC1 		71
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,6),I=1,MGL),&! Phototransformation from CDOC2 --> DIC			72	
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,7),I=1,MGL),&! Phototransformation from CDOC2 --> CNDOC 		73
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,8),I=1,MGL),&! Phototransformation from CDOC1 --> DIC 		74		
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,9),I=1,MGL), & ! Phototransformation from CDOC1 --> NCDOC 	75
						  !transformation of nitrogen
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,10),I=1,MGL),&! Phototransformation from CDOC3 --> CDOC2 		76
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,11),I=1,MGL),&! Phototransformation from CDOC3 --> CDOC1 		77
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,12),I=1,MGL),&! Phototransformation from CDOC3 --> DIC 		78
              (PHOTOCHEM_DOC_ACCUM_GL(I,K,13),I=1,MGL),&! Phototransformation from CDOC3 --> NCDOC 		79
              		 
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,14),I=1,MGL),&! Phototransformation from CDOC2 --> CDOC1 		80
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,15),I=1,MGL),&! Phototransformation from CDOC2 --> DIC		81	
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,16),I=1,MGL),&! Phototransformation from CDOC2 --> CNDOC 		82
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,17),I=1,MGL),&! Phototransformation from CDOC1 --> DIC 		83
						  (PHOTOCHEM_DOC_ACCUM_GL(I,K,18),I=1,MGL) ! Phototransformation from CDOC1 --> NCDOC		84
				  
						  
            !write sediment diagenesis model outputs   !LBcomment
         !IF(SEDIMENT_CALC .AND. K.eq.KBM1)THEN 
               
            !these are corresponding to the outputs of the excel version (Greg Pelletier)  

                !WRITE(UNIT_HIS+K,'(:8(E13.6e3,1x),E13.6e3)')        &!LB used to "E" rather than "F" for sedims
                !    !Print Quantity                                        Variable            FortranUnit            ExceUnit
                !                    !(( JPOC_GL(I,1)             &
                !                    !+ JPOC_GL(I,2)              &
                !                    !+ JPOC_GL(I,3))             &
                !                    ! /1000.0*2.667,I=1,MGL),    &        !'Jcin'                mgC/m^2/d            gO2/m^2/d
                !                    !((JPON_GL(I,1)              &
                !                    !+ JPON_GL(I,2)              &
                !                    !+ JPON_GL(I,3))     cd         &
                !                    ! /1000.0,I=1,MGL),          &        !'Jnin'                *mgN/m^2/d            gN/m^2/d
                !                    !((JPOP_GL(I,1)              &
                !                    !+ JPOP_GL(I,2)              &
                !                    !+ JPOP_GL(I,3))             &
                !                    ! /1000.0,I=1,MGL),          &    !'Jpin'                *mgP/m^2/d            gP/m^2/d
                !            (JPOC_GL(I,1)/1000.0*2.667,I=1,MGL), &       !'Jcin1'                mgC/m^2/d            gO2/m^2/d
                !            (JPOC_GL(I,2)/1000.0*2.667,I=1,MGL), &       !'Jcin2'                mgC/m^2/d            gO2/m^2/d
                !            (JPOC_GL(I,3)/1000.0*2.667,I=1,MGL), &       !'Jcin3'                mgC/m^2/d            gO2/m^2/d
                !            (JPON_GL(I,1)/1000.0,I=1,MGL),       &       !'Jnin1'                *mgN/m^2/d            gN/m^2/d
                !            (JPON_GL(I,2)/1000.0,I=1,MGL),       &       !'Jnin2'                *mgN/m^2/d            gN/m^2/d
                !            (JPON_GL(I,3)/1000.0,I=1,MGL),       &       !'Jnin3'                *mgN/m^2/d            gN/m^2/d
                !            (JPOP_GL(I,1)/1000.0,I=1,MGL),       &       !'Jpin1'                *mgP/m^2/d            gP/m^2/d
                !            (JPOP_GL(I,2)/1000.0,I=1,MGL),       &       !'Jpin2'                *mgP/m^2/d            gP/m^2/d
                !            (JPOP_GL(I,3)/1000.0,I=1,MGL),       &       !'Jpin3'                *mgP/m^2/d            gP/m^2/d
                !            (JPOS_GL(I)/1000.0,I=1,MGL),        &        !'Jsin                 *mgSi/m^2/d            gSi/m^2/d
                !            (C2_GL(I,KBM1,27),I=1,MGL),         &        !'O20'                mgO2/L                mgO2/L
                !            (D_GL(I),I=1,MGL),                  &        !'Depth'            m                    m
                !            (T_GL(I,KBM1),I=1,MGL),             &        !'Tw'                degC                degC
                !            (C2_GL(I,KBM1,13),I=1,MGL),         &        !'NH30'                *gN/m^3            mgN/L
                !            (C2_GL(I,KBM1,14),I=1,MGL),         &        !'NO30'                *gN/m^3            mgN/L
                !            (C2_GL(I,KBM1,20),I=1,MGL),         &        !'PO40'                *gP/m^3            mgP/L
                !            (C2_GL(I,KBM1,29),I=1,MGL),         &        !'SI0'                *gSi/m^3             mgSi/L        (dissiolved SIAT)
                !            (CH40_GL(I),I=1,MGL),               &        !'CH40'                gO2/m^3                mgO2/L
                !            (S_GL(I,KBM1),I=1,MGL),             &        !'SALw'              ppt            ppt            !i.e. SALT(I,KWC), ppt
                !            (SODTM1S_GL(I),I=1,MGL),            &        !'SOD'                gO2/m^2/day            gO2/m^2/d
                !            (JNH4_GL(I)/1000.0,I=1,MGL),        &        !'Jnh4                *mgN/m^2/day        gN/m^2/d
                !            (JNO3_GL(I)/1000.0,I=1,MGL),        &        !'Jno3                *mgN/m^2/day        gN/m^2/d
                !            (BENDEN_GL(I),I=1,MGL),             &        !'JDenitT'            gN/m^2/day            gN/m^2/d
                !            (JCH4_GL(I),I=1,MGL),               &        !'Jch4'                gO2/m^2/day            gO2/m^2/d
                !            (JCH4G_GL(I),I=1,MGL),              &        !'Jch4g'            gO2/m^2/day            gO2/m^2/d
                !            (JHS_GL(I),I=1,MGL),                &        !'Jhs'                gO2/m^2/day            gO2/m^2/d
                !            (JPO4_GL(I)/1000.0,I=1,MGL),        &        !'Jpo4'                *mgP/m^2/day        gP/m^2/d
                !            (JSI_GL(I)/1000.0,I=1,MGL),         &        !'Jsi'                *mgSi/m^2/day        gSi/m^2/d
                !            (NH41_GL(I)/1000.0,I=1,MGL),        &        !'NH3(1)'            *mgN/m^3             mgN/L        (dissolved)    
                !            (NH42_GL(I)/1000.0,I=1,MGL),        &        !'NH3(2)'            *mgN/m^3             mgN/L        (dissolved)
                !            (NO31_GL(I)/1000.0,I=1,MGL),        &        !'NO3(1)'            *mgN/m^3             mgN/L        (dissolved)
                !            (NO32_GL(I)/1000.0,I=1,MGL),        &        !'NO3(2)'            *mgN/m^3             mgN/L        (dissolved)
                !            (PO41_GL(I)/1000.0,I=1,MGL),        &        !'PO4(1)'            *mgP/m^3             mgP/L        (dissolved)
                !            (PO42_GL(I)/1000.0,I=1,MGL),        &        !'PO4(2)'            *mgP/m^3             mgP/L        (dissolved)
                !            (SI1_GL(I)/1000.0,I=1,MGL),         &        !'Si(1)'            *mgSi/m^3             mgSi/L        (dissolved)
                !            (SI2_GL(I)/1000.0,I=1,MGL),         &        !'Si(2)'            *mgSi/m^3             mgSi/L        (dissolved)
                !            (CH41TM1S_GL(I),I=1,MGL),           &        !'CH4(1)'            gO2/m^3             mgO2/L        (dissolved)
                !            (CH42_GL(I),I=1,MGL),               &        !'CH4(2)'            gO2/m^3             mgO2/L        (dissolved)
                !            (HS1_GL(I),I=1,MGL),                &        !'HS(1)'            gO2/m^3             mgO2/L        (dissolved)
                !            (HS2_GL(I),I=1,MGL),                &        !'HS(2)'            gO2/m^3             mgO2/L        (dissolved)
                !            (CPOC_GL(I,1)/1000.0*2.667,I=1,MGL), &        !'POC2(1)'            *mgC/m^3             gO2/m^3        (G1 of POC in layer 2)
                !            (CPOC_GL(I,2)/1000.0*2.667,I=1,MGL), &        !'POC2(2)'             *mgC/m^3             gO2/m^3        (G2 of POC in layer 2)
                !            (CPOC_GL(I,3)/1000.0*2.667,I=1,MGL), &        !'POC2(3)'             *mgC/m^3             gO2/m^3        (G3 of POC in layer 2)
                !            (CPON_GL(I,1)/1000.0,I=1,MGL),      &        !'PON2(1)'            *mgN/m^3             gN/m^3        (G1 of PON in layer 2)
                !            (CPON_GL(I,2)/1000.0,I=1,MGL),      &        !'PON2(2)'             *mgN/m^3             gN/m^3        (G2 of PON in layer 2)
                !            (CPON_GL(I,3)/1000.0,I=1,MGL),      &        !'PON2(3)'             *mgN/m^3             gN/m^3        (G3 of PON in layer 2)
                !            (CPOP_GL(I,1)/1000.0,I=1,MGL),      &        !'POP2(1)'            *mgP/m^3             gP/m^3        (G1 of POP in layer 2)
                !            (CPOP_GL(I,2)/1000.0,I=1,MGL),      &        !'POP2(2)'            *mgP/m^3             gP/m^3        (G2 of POP in layer 2)
                !            (CPOP_GL(I,3)/1000.0,I=1,MGL),      &        !'POP2(3)'             *mgP/m^3             gP/m^3        (G3 of POP in layer 2)
                !            (CPOS_GL(I)/1000.0,I=1,MGL),        &        !'POS2'                *mgSi/m^3             gSi/m^3        (particulate biogenic Silicate in layer2)
                !            (HSED1_GL(I)*1000.0,I=1,MGL),       &        !'H1'                m            mm
                !            (BENSTRTM1S_GL(I),I=1,MGL)                  !'BEN_STR'             dimensionless        dimensionless
                        !ENDIF !end LBcomment

                       FLUSH(UNIT_HIS+K)
                   CASE(.FALSE.)
                   
                       WRITE(UNIT_HIS,'(:F8.4,1x,I1,1x,I8/(:8(F12.6,1x),F12.6))')JDAY_TMP,4,MGL,  &
                                     (H_GL(I),I=1,MGL),            &  !Wen Long output still water depth
                                     (EL_GL(I),I=1,MGL),           &  !Wen Long output surface elevation
                                     (D_GL(I),I=1,MGL),            &  !Wen Long output total depth
                                     (-D_GL(I)*ZZ(K),I=1,MGL),     &  !Wen Long,depth of this layer (from surface)
                                     (C2_GL(I,K,27),I=1,MGL),      &  !DOXG
                                     !(C2_GL(I,K,9),I=1,MGL),       &  !LDOC
                                     (C2_GL(I,K,4),I=1,MGL),       &  !ALG1
                                     (C2_GL(I,K,5),I=1,MGL),       &  !ALG2
                                     (C2_GL(I,K,13),I=1,MGL),      &  !NH4
                                     (C2_GL(I,K,14),I=1,MGL),      &  !NO3
                                     (C2_GL(I,K,20),I=1,MGL),      &  !PO4
!                                    (total_netPP_GL(I),I=1,MGL),  &  !
                                     (T_GL(I,K),I=1,MGL),          &  !Temperature
                                     (S_GL(I,K),I=1,MGL),          &  !Salinity
                                     (IAVG_GL(I,K),I=1,MGL),       &  !PAR (E/m^2/day)
                                     !(C2_GL(I,K,10),I=1,MGL),      &  !NCDOC1
                                     (C2_GL(I,K,11),I=1,MGL),      &  !LPOC
                                     (C2_GL(I,K,12),I=1,MGL),      &   !RPOC
	
									 (C2_GL(I,K,9),I=1,MGL),       &  !CDOC1
									 (C2_GL(I,K,10),I=1,MGL),      &  !NCDOC1
									 (C2_GL(I,K,33),I=1,MGL),       &  !CDOC2
									 (C2_GL(I,K,34),I=1,MGL),      &  !NCDOC2
									 (C2_GL(I,K,35),I=1,MGL),       &  !CDOC3
									 (C2_GL(I,K,36),I=1,MGL),      &  !NCDOC3
									
									 (C2_GL(I,K,16),I=1,MGL),       &  !CDON1
									 (C2_GL(I,K,17),I=1,MGL),      &  !NCDON1
									 (C2_GL(I,K,37),I=1,MGL),       &  !CDON2
									 (C2_GL(I,K,38),I=1,MGL),      &  !NCDON2
									 (C2_GL(I,K,39),I=1,MGL),       &  !CDON3
									 (C2_GL(I,K,40),I=1,MGL),      &  !NCDON3
									
									 (C2_GL(I,K,21),I=1,MGL),       &  !CDOP1
									 (C2_GL(I,K,22),I=1,MGL),      &  !NCDOP1
									 (C2_GL(I,K,41),I=1,MGL),       &  !CDOP2
									 (C2_GL(I,K,42),I=1,MGL),      &  !NCDOP2
									 (C2_GL(I,K,43),I=1,MGL),       &  !CDOP3
									 (C2_GL(I,K,44),I=1,MGL) !,  &      !NCDOP3
									 
                                     
                      IF(SEDIMENT_CALC)THEN  !
                                  !
            !these are the outputs of the excel version (Greg Pelletier) 
            !                                
            !        Print Quantity                        Variable            FortranUnit            ExceUnit
                !!',',(JPOC(I,1)+JPOC(I,2)+JPOC(I,3))/1000.0*2.667,    &    
                                                        !'Jcin'                mgC/m^2/d            gO2/m^2/d
                !!',',(JPON(I,1)+JPON(I,2)+JPON(I,3))/1000.0,        &
                                                        !'Jnin'                *mgN/m^2/d            gN/m^2/d
                !!',',(JPOP(I,1)+JPOP(I,2)+JPOP(I,3))/1000.0,        &
                                                        !'Jpin'                *mgP/m^2/d            gP/m^2/d
                !!',',JPOS(I)/1000.0,&                    !'Jsin                 *mgSi/m^2/d            gSi/m^2/d
                !!',',O20,            &                    !'O20'                mgO2/L                mgO2/L
                !!',',D(I),            &                    !'Depth'            m                    m
                !!',',T(I,KBM1),    &                    !'Tw'                degC                degC
                !!',',NH40/1000.0,    &                    !'NH30'                *mgN/m^3            mgN/L
                !!',',NO30/1000.0,    &                    !'NO30'                *mgN/m^3            mgN/L
                !!',',PO40/1000.0,    &                    !'PO40'                *mgP/m^3            mgP/L
                !!',',SI0/1000.0,    &                    !'SI0'                *mgSi/m^3             mgSi/L        (dissiolved SIAT)
                !!',',CH40,            &                    !'CH40'                gO2/m^3                mgO2/L
                !!',',SAL,            &                    !'SALw'              ppt                    ppt            !i.e. SALT(I,KWC), ppt
                !!',',SOD,            &                    !'SOD'                gO2/m^2/day            gO2/m^2/d
                !!',',JNH4/1000.0,    &                    !'Jnh4                *mgN/m^2/day        gN/m^2/d
                !!',',JNO3/1000.0,    &                    !'Jno3                *mgN/m^2/day        gN/m^2/d
                !!',',BENDEN(1),    &                    !'JDenitT'            gN/m^2/day            gN/m^2/d
                !!',',JCH4,            &                    !'Jch4'                gO2/m^2/day            gO2/m^2/d
                !!',',JCH4G,        &                    !'Jch4g'            gO2/m^2/day            gO2/m^2/d
                !!',',JHS,            &                    !'Jhs'                gO2/m^2/day            gO2/m^2/d
                !!',',JPO4/1000.0,    &                    !'Jpo4'                *mgP/m^2/day        gP/m^2/d
                !!',',JSI/1000.0,    &                    !'Jsi'                *mgSi/m^2/day        gSi/m^2/d
                !!',',NH41/1000.0,    &                    !'NH3(1)'            *mgN/m^3             mgN/L        (dissolved)    
                !!',',NH42/1000.0,    &                    !'NH3(2)'            *mgN/m^3             mgN/L        (dissolved)
                !!',',NO31/1000.0,    &                    !'NO3(1)'            *mgN/m^3             mgN/L        (dissolved)
                !!',',NO32/1000.0,    &                    !'NO3(2)'            *mgN/m^3             mgN/L        (dissolved)
                !!',',PO41/1000.0,    &                    !'PO4(1)'            *mgP/m^3             mgP/L        (dissolved)
                !!',',PO42/1000.0,    &                    !'PO4(2)'            *mgP/m^3             mgP/L        (dissolved)
                !!',',SI1/1000.0,    &                    !'Si(1)'            *mgSi/m^3             mgSi/L        (dissolved)
                !!',',SI2/1000.0,    &                    !'Si(2)'            *mgSi/m^3             mgSi/L        (dissolved)
                !!',',CH41,            &                    !'CH4(1)'            gO2/m^3             mgO2/L        (dissolved)
                !!',',CH42,            &                    !'CH4(2)'            gO2/m^3             mgO2/L        (dissolved)
                !!',',HS1,            &                    !'HS(1)'            gO2/m^3             mgO2/L        (dissolved)
                !!',',HS2,            &                    !'HS(2)'            gO2/m^3             mgO2/L        (dissolved)
                !!',',POC1/1000.0*2.667,&                !'POC2(1)'            *mgC/m^3             gO2/m^3        (G1 of POC in layer 2)
                !!',',POC2/1000.0*2.667,&                !'POC2(2)'             *mgC/m^3             gO2/m^3        (G2 of POC in layer 2)
                !!',',POC3/1000.0*2.667,&                !'POC2(3)'             *mgC/m^3             gO2/m^3        (G3 of POC in layer 2)
                !!',',PON1/1000.0,    &                    !'PON2(1)'            *mgN/m^3             gN/m^3        (G1 of PON in layer 2)
                !!',',PON2/1000.0,    &                    !'PON2(2)'             *mgN/m^3             gN/m^3        (G2 of PON in layer 2)
                !!',',PON3/1000.0,    &                    !'PON2(3)'             *mgN/m^3             gN/m^3        (G3 of PON in layer 2)
                !!',',POP1/1000.0,    &                    !'POP2(1)'            *mgP/m^3             gP/m^3        (G1 of POP in layer 2)
                !!',',POP2/1000.0,    &                    !'POP2(2)'            *mgP/m^3             gP/m^3        (G2 of POP in layer 2)
                !!',',POP3/1000.0,    &                    !'POP2(3)'             *mgP/m^3             gP/m^3        (G3 of POP in layer 2)
                !!',',PSISED/1000.0,&                    !'POS2'                *mgSi/m^3             gSi/m^3        (particulate biogenic Silicate in layer2)
                !!',',HSED1(I)*1000.0,&                    !'H1'                m                    mm
                !!',',BENSTR                             !'BEN_STR'             dimensionless        dimensionless
                                     
    
                      ENDIF
                      
                      FLUSH(UNIT_HIS)
             ENDSELECT
           ENDDO
         ENDIF !END PRINT_HIS SWITCH
       ENDIF  !END MSR SWITCH

       PRINT_STN = .FALSE. 
       IF(NTSTN>NTSTN_OLD)THEN
           CTR_STN   = CTR_STN+1                      !Counter of steps in station output file
           PRINT_STN = .TRUE. 
           NTSTN_OLD  = NTSTN
       ENDIF

       IF(MSR) THEN

        WRITE(*,*)'CTR_STN=', CTR_STN
        IF(CTR_STN.eq.1) THEN   !Open file for the first time 

!WLong: Instead of opening one file for each station, here I simply use a single file for all 
!       stations

            WRITE(*,*)'Opening station output file:',TRIM(STNFN)
            OPEN(UNIT=UNIT_STN, FILE=TRIM(STNFN),STATUS='UNKNOWN')
            WRITE(UNIT_STN,*)"Nstation,Nlayer"
            WRITE(UNIT_STN,*)Nstation, KBM1       !number of stations and number of layers

            TITLE_STN='Variables="StationID","Node","Layer","depth(m)","DO","NO3","NH4","Alg1","Alg2","LDOC","RDOC","LPOC","RPOC","PO4"'
!Wen Long debug algae
            TITLE_STN=TRIM(TITLE_STN)//',"T","S","P1","P2","BM1","BM2","NL1","NL2","PL1","PL2","FI1","FI2","B1SZ","B2SZ","B1LZ","B2LZ","PR1","PR2","IAVG"'
            TITLE_STN=TRIM(TITLE_STN)//',"DDOC","DCOD","NITRIF","DOSZ","DOLZ","DOPR","DORALG","RESP","FTCOD","COD"'
            TITLE_STN=TRIM(TITLE_STN)//',"REAERDO","BFLUXDO","BFLUXNO3","BFLUXNH4","IK1","IK2"'
            IF(SEDIMENT_CALC)THEN
               TITLE_STN=TRIM(TITLE_STN)//',"Jcin1","Jcin2","Jcin3","Jnin1","Jnin2","Jnin3","Jpin1","Jpin2","Jpin3"'
               TITLE_STN=TRIM(TITLE_STN)//',"Jsin","O20","Depth","Tw","NH30","NO30","PO40","SI0","CH40"'
               TITLE_STN=TRIM(TITLE_STN)//',"SALw","SOD","Jnh4","Jno3","JDenitT"'
               TITLE_STN=TRIM(TITLE_STN)//',"Jch4"'
               TITLE_STN=TRIM(TITLE_STN)//',"Jch4g","Jhs","Jpo4","Jsi","NH31","NH32","NO31","NO32"'
               TITLE_STN=TRIM(TITLE_STN)//',"PO41","PO42","Si1","Si2","CH41","CH42","HS1","HS2","POC21","POC22","POC23""'
               TITLE_STN=TRIM(TITLE_STN)//',"PON21","PON22","PON23,"POP21","POP22","POP23","POS2","H1","BEN_STR"'
               !TITLE_STN=TRIM(TITLE_STN)//',"FRPOC1","FRPOC2","FRPON1","FRPON2","FRPOP1","FRPOP2"'
            ENDIF
            WRITE(UNIT_STN,'(A)')TRIM(TITLE_STN)

            OPEN_STN = .TRUE.
        ENDIF !END CTR_STN switch
		

        IF(PRINT_STN) THEN

          WRITE(*,*)'Writing to station outputs, JDAY=',JDAY_TMP
          WRITE(UNIT_STN,*)Nstation,KBM1,JDAY_TMP
          DO I=1,Nstation
             DO K=1,KBM1
                WRITE(UNIT_STN,'(:I8,1x,I8,1x,I8/(:8(F12.6,1X),F12.6))')     &
                                                   I,NstationNum_GL(I),K,    &  !ID, global station node number, layer 
                                       -D_GL(NstationNum_GL(I))*ZZ(K),       &  !1 depth of k'th layer, positive from surface downward
                                        C2_GL(NstationNum_GL(I),K,27),       &  !2 DO
                                        C2_GL(NstationNum_GL(I),K,14),       &  !3 NO3
                                        C2_GL(NstationNum_GL(I),K,13),       &  !4 NH4
                                        C2_GL(NstationNum_GL(I),K,4),        &  !5 ALG1
                                        C2_GL(NstationNum_GL(I),K,5),        &  !6 ALG2
                                        C2_GL(NstationNum_GL(I),K,9),        &  !7 LDOC
                                        C2_GL(NstationNum_GL(I),K,10),       &  !8 RDOC
                                        C2_GL(NstationNum_GL(I),K,11),       &  !9 LPOC
                                        C2_GL(NstationNum_GL(I),K,12),       &  !10 RPOC
                                        C2_GL(NstationNum_GL(I),K,20)        &  !11 PO4
                                       ,T_GL(NstationNum_GL(I),K),    &         !12 temperature 
                                        S_GL(NstationNum_GL(I),K),    &         !13 salt
                                        P1_GL(NstationNum_GL(I),K),   &         !14 Wen Long debug algae   ALG1 gross growth rate (1/day)
                                        P2_GL(NstationNum_GL(I),K),   &         !15 ALG2 growth rate (1/day)
                                        BM1_GL(NstationNum_GL(I),K),  &         !16 ALG1 basal metabolism rate (1/day)
                                        BM2_GL(NstationNum_GL(I),K),  &         !17 ALG2                       (1/day)
                                        NL1_GL(NstationNum_GL(I),K),  &         !18 Nitrogen limitation of ALG 1 (nondimensional)
                                        NL2_GL(NstationNum_GL(I),K),  &         !19                            2
                                        PL1_GL(NstationNum_GL(I),K),  &         !20 Phosphorus                 1
                                        PL2_GL(NstationNum_GL(I),K),  &         !21         !                            2
                                        FI1_GL(NstationNum_GL(I),K),  &         !22            ! Light                      1
                                        FI2_GL(NstationNum_GL(I),K),  &         !23           !                            2
                                        B1SZ_GL(NstationNum_GL(I),K), &         !24            ! microzoop consumption      1  (mgC/L/day)
                                        B2SZ_GL(NstationNum_GL(I),K), &         !25            !                            2
                                        B1LZ_GL(NstationNum_GL(I),K), &         !26            ! macrozoop                  1
                                        B2LZ_GL(NstationNum_GL(I),K), &         !27            !                            2
                                        PR1_GL(NstationNum_GL(I),K),  &         !28            ! predation rate             1  (mgC/L/day)
                                        PR2_GL(NstationNum_GL(I),K),  &         !29            !                            2 
                                        IAVG_GL(NstationNum_GL(I),K)  &         !30            ! light intensity               (E/m^2/day)
                                      ,-DDOC_GL(NstationNum_GL(I),K)  &         !31 oxygen loss due to DOC resmin      (gO2/m^3/day)
                                      ,-DCOD_GL(NstationNum_GL(I),K)  &         !32 oxygen loss due to COD oxidation  (gO2/m^3/day)
                                      ,-NITRIF_GL(NstationNum_GL(I),K)&         !33 oxygen loss due to nitrification (gO2/m^3/day)
                                      ,-DOSZ_GL(NstationNum_GL(I),K)  &         !34 oxygen loss due to micro zoop respiration (gO2/m^3/day)
                                      ,-DOLZ_GL(NstationNum_GL(I),K)  &         !35 oxygen loss due to macro zoop respiration (gO2/m^3/day)
                                      ,-DOPR_GL(NstationNum_GL(I),K)  &         !36 oxygen loss due to predation on alg (gO2/m^3/day)
                                      ,+DORALG_GL(NstationNum_GL(I),K)&         !37 oxygen production due to photosynthesis -photorespiration -basal metabolism
                                      ,-RESP_GL(NstationNum_GL(I),K)  &         !38 sum of oxygen loss due to all respirations (gO2/m^3/day)
                                      , FTCOD_GL(NstationNum_GL(I),K) &         !39 COD decay rate after applying temperature control (1/day)
                                      , COD_GL(NstationNum_GL(I),K)   &         !40 COD cocnentration (mgO2-equiv/L)
                                      , REAERDO_GL(NstationNum_GL(I),K) &       !41 areation rate of DOXG (gO2/m^3/day)
                                      , BFLUX_GL(NstationNum_GL(I),8)/(D_GL(NstationNum_GL(I))*DZ(KBM1)) & !42 Benthic DOXG flux (gO2/m^3/day),
                                                                                                           !   D*DZ layer thickness
                                      , BFLUX_GL(NstationNum_GL(I),3)/(D_GL(NstationNum_GL(I))*DZ(KBM1)) & !43 Benthic NO3 flux (gN/m^3/day)
                                      , BFLUX_GL(NstationNum_GL(I),2)/(D_GL(NstationNum_GL(I))*DZ(KBM1)) & !44 Benthic NH4 flux (gN/m^3/day)
                                      , IK1_GL(NstationNum_GL(I),K) &                                      !45
                                      , IK2_GL(NstationNum_GL(I),K)                                        !46
             ENDDO
             
             !write sediment diagenesis model outputs   !LBcomment
			 
             IF(SEDIMENT_CALC)THEN 
			 
			 !write(*,*)'made to sed dom loop' ! B Clark debug
			!	IF(SED_DOM) THEN
			!		write(*,*)'inside sed dom out loop' ! B clark debug
			!		
			!		CALL SED_DOM_OUTPUT	
			!	ENDIF
            
             !!pre-calculation of Jcin,Jnin,Jpin IN FORTRAN UNITS
             !Jcin = JPOC_GL(NstationNum_GL(I),1)                &
             !       +JPOC_GL(NstationNum_GL(I),2)               &
             !       +JPOC_GL(NstationNum_GL(I),3)
             !Jnin = JPON_GL(NstationNum_GL(I),1)                &
             !        +JPON_GL(NstationNum_GL(I),2)              &
             !        +JPON_GL(NstationNum_GL(I),3)
             !Jpin = JPOP_GL(NstationNum_GL(I),1)                &
             !        +JPOP_GL(NstationNum_GL(I),2)              &
             !        +JPOP_GL(NstationNum_GL(I),3)
                     
                     
                !these are corresponding to the outputs of the excel version (Greg Pelletier) 
                !write(*,*)'LBnote: variables with formatting problems:'  !The following variables are quite large!
                !write(*,*)CPOC_GL(NstationNum_GL(1),3)/1000.0*2.667,CPON_GL(NstationNum_GL(1),3)/1000.0,CPOP_GL(NstationNum_GL(1),3)/1000.0

                !these are corresponding to the outputs of the excel version (Greg Pelletier)

                !WRITE(UNIT_STN,'(:8(F12.6,1X),F12.6)')         &!LB commented and tried to fix formatting issue switching to "E" rather than using "F"
                WRITE(UNIT_STN,'(:8(E13.6e3,1X),E13.6e3)')          &
                    !Print Quantity                                        Variable            FortranUnit            ExceUnit
                    !Jcin/1000.0*2.667,                          &   !'Jcin'                mgC/m^2/d            gO2/m^2/d
                    !Jnin/1000.0,                                &  !'Jnin'                *mgN/m^2/d            gN/m^2/d
                    !Jpin/1000.0,                                &  !'Jpin'                *mgP/m^2/d            gP/m^2/d
                    JPOC_GL(NstationNum_GL(I),1)/1000.0*2.667,    &      !'Jcin1'               mgC/m^2/d            gO2/m^2/d
                    JPOC_GL(NstationNum_GL(I),2)/1000.0*2.667,    &      !'Jcin2'               mgC/m^2/d            gO2/m^2/d
                    JPOC_GL(NstationNum_GL(I),3)/1000.0*2.667,    &      !'Jcin3'               mgC/m^2/d            gO2/m^2/d
                    JPON_GL(NstationNum_GL(I),1)/1000.0,         &      !'Jnin1'               *mgN/m^2/d            gN/m^2/d
                    JPON_GL(NstationNum_GL(I),2)/1000.0,         &      !'Jnin2'               *mgN/m^2/d            gN/m^2/d
                    JPON_GL(NstationNum_GL(I),3)/1000.0,         &      !'Jnin3'               *mgN/m^2/d            gN/m^2/d
                    JPOP_GL(NstationNum_GL(I),1)/1000.0,         &      !'Jpin1'               *mgP/m^2/d            gP/m^2/d
                    JPOP_GL(NstationNum_GL(I),2)/1000.0,         &      !'Jpin2'               *mgP/m^2/d            gP/m^2/d
                    JPOP_GL(NstationNum_GL(I),3)/1000.0,         &      !'Jpin3'               *mgP/m^2/d            gP/m^2/d
                    JPOS_GL(NstationNum_GL(I))/1000.0,           &      !'Jsin                 *mgSi/m^2/d            gSi/m^2/d
                    C2_GL(NstationNum_GL(I),KBM1,27),            &       !'O20'                mgO2/L                mgO2/L
                    D_GL(NstationNum_GL(I)),                    &        !'Depth'            m                    m
                    T_GL(NstationNum_GL(I),KBM1),                &       !'Tw'                degC                degC
                    C2_GL(NstationNum_GL(I),KBM1,13),            &        !'NH30'                *gN/m^3            mgN/L
                    C2_GL(NstationNum_GL(I),KBM1,14),            &        !'NO30'                *gN/m^3            mgN/L
                    C2_GL(NstationNum_GL(I),KBM1,20),            &        !'PO40'                *gP/m^3            mgP/L
                    C2_GL(NstationNum_GL(I),KBM1,29),            &        !'SI0'                *gSi/m^3             mgSi/L        (dissiolved SIAT)
                    CH40_GL(NstationNum_GL(I)),                 &        !'CH40'             gO2/m^3                mgO2/L
                    S_GL(NstationNum_GL(I),KBM1),                &        !'SALw'              ppt                    ppt            !i.e. SALT(I,KWC), ppt
                    SODTM1S_GL(NstationNum_GL(I)),                &        !'SOD'                gO2/m^2/day            gO2/m^2/d
                    JNH4_GL(NstationNum_GL(I))/1000.0,            &        !'Jnh4                *mgN/m^2/day        gN/m^2/d
                    JNO3_GL(NstationNum_GL(I))/1000.0,            &        !'Jno3                *mgN/m^2/day        gN/m^2/d
                    BENDEN_GL(NstationNum_GL(I)),                &        !'JDenitT'            gN/m^2/day            gN/m^2/d
                    JCH4_GL(NstationNum_GL(I)),                    &        !'Jch4'                gO2/m^2/day            gO2/m^2/d
                    JCH4G_GL(NstationNum_GL(I)),                &        !'Jch4g'            gO2/m^2/day            gO2/m^2/d
                    JHS_GL(NstationNum_GL(I)),                    &        !'Jhs'                gO2/m^2/day            gO2/m^2/d
                    JPO4_GL(NstationNum_GL(I))/1000.0,            &        !'Jpo4'                *mgP/m^2/day        gP/m^2/d
                    JSI_GL(NstationNum_GL(I))/1000.0,            &        !'Jsi'                *mgSi/m^2/day        gSi/m^2/d
                    NH41_GL(NstationNum_GL(I))/1000.0,            &        !'NH3(1)'            *mgN/m^3             mgN/L        (dissolved)
                    NH42_GL(NstationNum_GL(I))/1000.0,            &        !'NH3(2)'            *mgN/m^3             mgN/L        (dissolved)
                    NO31_GL(NstationNum_GL(I))/1000.0,            &        !'NO3(1)'            *mgN/m^3             mgN/L        (dissolved)
                    NO32_GL(NstationNum_GL(I))/1000.0,            &        !'NO3(2)'            *mgN/m^3             mgN/L        (dissolved)
                    PO41_GL(NstationNum_GL(I))/1000.0,            &        !'PO4(1)'            *mgP/m^3             mgP/L        (dissolved)
                    PO42_GL(NstationNum_GL(I))/1000.0,            &        !'PO4(2)'            *mgP/m^3             mgP/L        (dissolved)
                    SI1_GL(NstationNum_GL(I))/1000.0,            &        !'Si(1)'            *mgSi/m^3             mgSi/L        (dissolved)
                    SI2_GL(NstationNum_GL(I))/1000.0,            &        !'Si(2)'            *mgSi/m^3             mgSi/L        (dissolved)
                    CH41TM1S_GL(NstationNum_GL(I)),                    &        !'CH4(1)'            gO2/m^3             mgO2/L        (dissolved)
                    CH42_GL(NstationNum_GL(I)),                    &        !'CH4(2)'            gO2/m^3             mgO2/L        (dissolved)
                    HS1_GL(NstationNum_GL(I)),                    &        !'HS(1)'            gO2/m^3             mgO2/L        (dissolved)
                    HS2_GL(NstationNum_GL(I)),                    &        !'HS(2)'            gO2/m^3             mgO2/L        (dissolved)
                    CPOC_GL(NstationNum_GL(I),1)/1000.0*2.667,    &        !'POC2(1)'            *mgC/m^3             gO2/m^3        (G1 of POC in layer 2)
                    CPOC_GL(NstationNum_GL(I),2)/1000.0*2.667,    &        !'POC2(2)'             *mgC/m^3             gO2/m^3        (G2 of POC in layer 2)
                    CPOC_GL(NstationNum_GL(I),3)/1000.0*2.667,    &        !'POC2(3)'             *mgC/m^3             gO2/m^3        (G3 of POC in layer 2)
                    CPON_GL(NstationNum_GL(I),1)/1000.0,        &        !'PON2(1)'            *mgN/m^3             gN/m^3        (G1 of PON in layer 2)
                    CPON_GL(NstationNum_GL(I),2)/1000.0,        &        !'PON2(2)'             *mgN/m^3             gN/m^3        (G2 of PON in layer 2)
                    CPON_GL(NstationNum_GL(I),3)/1000.0,        &        !'PON2(3)'             *mgN/m^3             gN/m^3        (G3 of PON in layer 2)
                    CPOP_GL(NstationNum_GL(I),1)/1000.0,        &        !'POP2(1)'            *mgP/m^3             gP/m^3        (G1 of POP in layer 2)
                    CPOP_GL(NstationNum_GL(I),2)/1000.0,        &        !'POP2(2)'            *mgP/m^3             gP/m^3        (G2 of POP in layer 2)
                    CPOP_GL(NstationNum_GL(I),3)/1000.0,        &        !'POP2(3)'             *mgP/m^3             gP/m^3        (G3 of POP in layer 2)
                    CPOS_GL(NstationNum_GL(I))/1000.0,          &        !'POS2'                *mgSi/m^3             gSi/m^3        (particulate biogenic Silicate in layer2)
                    HSED1_GL(NstationNum_GL(I))*1000.0,         &        !'H1'                m                    mm
                    BENSTRTM1S_GL(NstationNum_GL(I)) !,         &           !'BEN_STR'             dimensionless        dimensionless
                    !JPOC_GL(NstationNum_GL(I),1)/Jcin,          & !FRPOC1                                           dimensionless
                    !JPOC_GL(NstationNum_GL(I),2)/Jcin,          & !FRPOC2                                           dimensionless
                    !JPON_GL(NstationNum_GL(I),1)/Jnin,          & !FRPON1                                           dimensionless
                    !JPON_GL(NstationNum_GL(I),2)/Jnin,          & !FRPON2                                           dimensionless
                    !JPOP_GL(NstationNum_GL(I),1)/Jpin,          & !FRPOP1                                           dimensionless
                    !JPOP_GL(NstationNum_GL(I),2)/Jpin             !FRPOP2                                           dimensionless

            ENDIF  !end LBcomment
          ENDDO
         !Push the write through onto disk
          FLUSH(UNIT_STN)
        ENDIF  !END PRINT_STN SWITCH

       ENDIF  !END MSR SWITCH 

!rgl check if below is necessary
!        if(MSR) stop
! KURT GLAESEMANN DEBUG with 999 stop

     ENDIF !END AMOD(JDAY  ... ) SWITCH
   END SUBROUTINE WQM_OUTPUT
