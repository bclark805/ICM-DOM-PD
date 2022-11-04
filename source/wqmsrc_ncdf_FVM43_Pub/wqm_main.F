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

                ATMLDON,        &!
                ATMLDOP,        &!
                ATMNH4,         &!
                ATMNO3,         &!
                ATMPO4,         &!
                ATMRDON,        &!
                ATMRDOP,        &!


                DLWCMC,         &!
                DLWCMN,         &!
                DLWCMP,         &!

                IWCMC,          &!
                IWCMN,          &!
                IWCMP,          &!

                PRECIP,         &!
                SFATOT,         &!
                
                ANDC,           &!

                
                BENNH4,         &!
                BENNO3,         &!
                BENSA,          &!                
                BENPO4,         &!
                
                BNDN,           &!
                DENIT,          &!
         
                DTMP1,          &!

                INFLOW,         &!
                JCB,            &!
            !    MNLLDOC,                &!  ! B Clark moved into 3 pools in mod_wc_dom Sep 2015
           !     MNLRDOC,                &!  ! B Clark moved into 3 pools in mod_wc_dom Sep 2015
                PCFWS,          &!
                PPFWS,          &!
                ! RESPC,          &!

                BNDC,           &!

                BNDP,           &!
                DLWCSSMN,               &!
                DLWCSSMP,               &!

                S1FLXN,         &!
                S1FLXP,         &!
                S2FLXN,         &!
                S3FLXN,         &!

                BENDEN,         &!

                DFTMP,          &!

     
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
                ! DCOD,           &!
                ! ! DCOD_GL,        &!
                ! DDOC,           &!
                ! DDOC_GL,        &!

                IAVG_GL,        &!
				KDLAMB_int_GL,   & ! B Clark May 2018,
				SSI_GL,  & ! B Clark June 2018
                IK1,            &!
                IK1_GL,         &!
                IK2,            &!
                IK2_GL,         &!
                ! NITRIF,         &!
                ! NITRIF_GL,              &!
                PR1,            &!
                PR1_GL,         &!
                PR2,            &!
                PR2_GL,         &!


                COD_GL,         &!
                COURMX,         &!
                DIAGNOSTICS,    &!
                ! DOLZ_GL,        &!
                ! DOPR,           &!
                ! DOPR_GL,        &!
                ! DORALG,         &!
                ! DORALG_GL,      &!
                ! DOSZ_GL,        &!
                FTCOD,          &!
                FTCOD_GL,       &!
                REAERDO,        &!
                REAERDO_GL,     &!
                RESP,           &!
    
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

                
    USE MOD_OWQ, ONLY: INTKE,            &!
                        NXKEI,            &!
   
                        IAVG,            &!

                        OWQ_DEALLOC,    &!

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

            ! C1MIN,              &                !
            ! C1MAX,                &                !

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
#if defined(ZOOPLANKTON)   
   USE MOD_ZOOP, ONLY :            &
                     ANCSZ,      &!microzooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     ANCLZ,      &!mesozooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     APCSZ,      &!microzooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
                     APCLZ,      &!mesozooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
 
                     CLSZ,        & 
                     CLLZ,        &
 
                      B1SZ,        &!
                      B2SZ,        &!
                     !B3SZ,        &!
                     B1LZ,        &!
                     B2LZ,        &!
                     !B3LZ,        &!
                     DOSZ,        &!
                     DOLZ,        &!
 
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
 
                     SZ,        &!
                     LZ,        &!
                     ZOOP_DEALLOC
#endif                    
    USE MOD_KIN, ONLY :                &
                        TEMPER,        &
                        SOLIDS,        &
#if defined(ZOOPLANKTON)
                         ZOOPL,        &
#endif
                        PARTICULATE_CARBON,        &
                        NITROG,        &
#if defined(PHOSPHORUS)
                         PHOSPH,        &
#endif
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
        IERR,                &!
        ISEDMN,               &!
        ISEDMP,               &!
        ISEDMC,                &!
        SEDMN,              &!total sediment organic nitrogen  (kgN)
        SEDMP,              &!total sediment organic phosphorus (kgP)
        SEDMC,                &!total sediment organic carbon  (kgC)     
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

        SETUP_OBC ,        &!
        OBCS_DEALLOC    
   
   USE MOD_BCS, ONLY : &
 
        BCS_DEALLOC

  USE MOD_SAV, ONLY    :        &
  
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

            NGID,        &    !(:)           !!GLOBAL ID OF LOCAL NODE 

            EC,            &    !
            NC,            &    !
            BNC,         &    !
            EMAP,        &    !
            NMAP,         &    !

             NODE_MATCH_ALL, &  !

             EXCHANGE_ALL, & !
             GATHER , & !

             INIT_MPI_ENV,&!

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
                        DOC, DON, & 
#if defined(PHOSPHORUS)
                       DOP, &  
#endif
#if defined(ZOOPLANKTON)
                        ZOOP_DOM,  &
#endif
                        WC_CDOC1, WC_CDOC2, WC_CDOC3,   &
                        WC_NCDOC1, WC_NCDOC2, WC_NCDOC3,    &

						WC_CDON1, WC_CDON2, WC_CDON3, &
						WC_NCDON1, WC_NCDON2, WC_NCDON3,  &
						WC_CDOP1, WC_CDOP2, WC_CDOP3,  &
						WC_NCDOP1, WC_NCDOP2, WC_NCDOP3, MNLDOC,   &
						DTWCDOC1,DTWNCDOC1,DTWCDOC2,DTWNCDOC2,   &
						DTWCDOC3,DTWNCDOC3,   & 
						DTWCDON1,DTWNCDON1,DTWCDON2,DTWNCDON2,   &
						DTWCDON3,DTWNCDON3,   & 
						DTWCDOC1_out,DTWNCDOC1_out,DTWCDOC2_out,DTWNCDOC2_out,   &
						DTWCDOC3_out,DTWNCDOC3_out,   & 
						DTWCDON1_out,DTWNCDON1_out,DTWCDON2_out,DTWNCDON2_out,   &
                  DTWCDON3_out,DTWNCDON3_out,   & 
#if defined(PHOSPHORUS)				
		DTWCDOP1,DTWNCDOP1,DTWCDOP2,DTWNCDOP2,   &
						DTWCDOP3,DTWNCDOP3,   &                  
#endif
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
		
        ! CALL SED_SF_EXCHANGE_ALLOC
        
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
 !  C2(:,:,3)=1.0
   SSI(:,:)  = C2(:,:,3)
   B1(:,:)   = C2(:,:,4)
   B2(:,:)   = C2(:,:,5)
   B3(:,:)   = C2(:,:,6)

#if defined(ZOOPLANKTON)
      SZ(:,:)   = C2(:,:,7)
      LZ(:,:)   = C2(:,:,8)
#endif
   
   WC_CDOC1(:,:) = C2(:,:,9)
   WC_NCDOC1(:,:) = C2(:,:,10)
  
   !LDOC(:,:) = C2(:,:,9)
   !RDOC(:,:) = C2(:,:,10)
   LPOC(:,:) = C2(:,:,11)
   RPOC(:,:) = C2(:,:,12)
   NH4(:,:)  = C2(:,:,13)
   NO3(:,:)  = C2(:,:,14)
!    UREA(:,:) = C2(:,:,15)
   
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
     
   CALL OUT_NETCDF

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
#if defined(PHOSPHORUS)    
     DO K=1,KBM1
       DO I=1,MLOC
         Q1(I,K) = PIB1(I,K)/(B1(I,K)+1.0E-9)  !Phosphorus to Carbon Ration in Alg 1 (gP/gC)
         Q2(I,K) = PIB2(I,K)/(B2(I,K)+1.0E-9)  !Phosphorus to Carbon Ration in Alg 2 (gP/gC)
     !    Q3(I,K) = PIB3(I,K)/(B3(I,K)+1.0E-9)    !Phosphorus to Carbon Ration in Alg 1 (gP/gC)
         
         IF(SEDIMENT_CALC.AND.ACC(30)=='OFF')THEN
                Q1(I,K)=0.0243
         ENDIF
         IF(SEDIMENT_CALC.AND.ACC(31)=='OFF')THEN
                Q2(I,K)=0.0243
         ENDIF
        !  IF(SEDIMENT_CALC.AND.ACC(32)=='OFF')THEN
        !         Q3(I,K)=0.0243
        !  ENDIF
         
       ENDDO  
     ENDDO
#endif
	 IF(MSR)THEN
		TOTAL_TIME_PROFILE=DTIME(ELAPSED_PROFILE)
		NPROFILE=NPROFILE+1
	 ENDIF
	 
!******* Water column kinetics
!     IF (TEMPERATURE_CALC) CALL TEMPER(DTC(0,1,1))
!tykim-this makes change in temperature zero - need to comment out if want to calc temp
        DTC(:,:,1)=0.0
!In next section first subscript changed to 0 per KURT 
		
	 
       IF (SOLIDS_CALC)      CALL SOLIDS(DTC(0,1,3)) !WLong and LB, we should move this into the IF block below as well??
   !     DTC(:,:,3)=0.0    
	IF(MOD(ELTMS,WQ_DT).eq.0) THEN
! !RGL added sub cycling of ecosystem model, above
! !WL  IF(MSR) WRITE(*,*)'calling ecosystem model'

#if defined(ZOOPLANKTON)
       IF (ZOO_CALC)    THEN     
	      CALL ZOOPL(DTC(0,1,7),DTC(0,1,8))
		  CALL ZOOP_DOM
     ENDIF
#endif     
	   
       IF (ALGAE_CALC)         &
        CALL ALGAE(DTC(0,1,4),DTC(0,1,5),DTC(0,1,6))
       IF (CARBON_CALC)     THEN 
	       CALL PARTICULATE_CARBON(DTC(0,1,11),DTC(0,1,12) )   ! DTC(0,1,9) DTC(0,1,10)  ! old DOC
	       CALL DOC
		     DTC(:,:,9) = DTWCDOC1
			   DTC(:,:,10) = DTWNCDOC1
			   DTC(:,:,33) = DTWCDOC2
			   DTC(:,:,34) = DTWNCDOC2
			   DTC(:,:,35) = DTWCDOC3
			   DTC(:,:,36) = DTWNCDOC3
	     ENDIF 
       IF (NITROGEN_CALC)     THEN
	        CALL NITROG(DTC(0,1,13) ,DTC(0,1,14), DTC(0,1,18),DTC(0,1,19))
	        CALL DON
          DTC(:,:,16) = DTWCDON1    
			    DTC(:,:,17) = DTWNCDON1   
			    DTC(:,:,37) = DTWCDON2   
			    DTC(:,:,38) = DTWNCDON2 
			    DTC(:,:,39) = DTWCDON3
			    DTC(:,:,40) = DTWNCDON3
			 ENDIF
#if defined(PHOSPHORUS)		
      IF (PHOSPHORUS_CALC)    THEN  
          CALL PHOSPH(DTC(0,1,20),DTC(0,1,23),&  !DTC(0,1,21),DTC(0,1,22), these are old ldop and rdop
                  DTC(0,1,24),DTC(0,1,25),DTC(0,1,30),DTC(0,1,31),    &
                  DTC(0,1,32))
		
	       CALL DOP 
		   
		    DTC(:,:,21) = DTWCDOP1
				DTC(:,:,22) = DTWNCDOP1
				DTC(:,:,41) = DTWCDOP2
				DTC(:,:,42) = DTWNCDOP2
				DTC(:,:,43) = DTWCDOP3
				DTC(:,:,44) = DTWNCDOP3
				
	     ENDIF
#endif
         IF (COD_CALC)         CALL CODMND(DTC(0,1,26))
	   
         IF (OXYGEN_CALC)      CALL OXYGEN(DTC(0,1,27))
        IF (SILICA_CALC)      CALL SILICA(DTC(0,1,28),DTC(0,1,29))

         DTM= DTC    
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
     !        write(*,*)  'Sed Calc is true, doing steady state' 
             ELSE
               CALL SED_CALC(.FALSE.) !FALSE means use time-varying calculation
      !       write(*,*)  'Sed Calc is FALSE, doing time series calc' 
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

     IF(JDAY.GE.TMEND) THEN
       END_RUN=.TRUE. 
                    
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
     ENDIF

     IF(MSR)THEN
         IF(END_RUN)WRITE(*,*)'MYID =',MYID,'END_RUN=',END_RUN   !master print its own END_RUN
     ENDIF

   !   IF (JDAY >= DLTD(DLTDP+1)) THEN
   !     DLTDP = DLTDP+1
   !     IF (.NOT.AUTO_STEPPING) DLT = DLTVAL(DLTDP)
   !     MXDLT = DLTMAX(DLTDP)
   !     FNDLT = DLTFTN(DLTDP)
   !   ENDIF

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
	  
#else
DTNH4_OUT=DTC(0,1,13) 
DTNO3_OUT=DTC(0,1,14)
DTLPOC_OUT=DTC(0,1,11)
DTRPOC_OUT=DTC(0,1,12)
DTLPON_OUT=DTC(0,1,18)
DTRPON_OUT=DTC(0,1,19)
#endif
	    
	  
!RGL commented below
! T(:,:)    = C2(:,:,1)
!   SALT(:,:) = C2(:,:,2)
  ! C2(:,:,3)=1.0
   SSI(:,:)  = C2(:,:,3)
   B1(:,:)   = C2(:,:,4)
   B2(:,:)   = C2(:,:,5)
   B3(:,:)   = C2(:,:,6)
#if  defined(ZOOPLANKTON)
   SZ(:,:)   = C2(:,:,7)
   LZ(:,:)   = C2(:,:,8)
#endif
   
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
#if defined(PHOSPHORUS)   
   PIP(:,:)  = C2(:,:,25)
#endif   
   COD(:,:)  = C2(:,:,26)
   DOXG(:,:)   = C2(:,:,27)
   SIUPB(:,:)   = C2(:,:,28) 
   SIAT(:,:)   = C2(:,:,29)
#if defined(PHOSPHORUS)     
   PIB1(:,:) = C2(:,:,30)
   PIB2(:,:) = C2(:,:,31)
   PIB3(:,:) = C2(:,:,32)
#endif   
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
   
!#if defined(XFLUX_FLAG)   
! write(*,*)' Calcuating the XFLUX out'
!Added by B Clark to time integrate xflux for each WQ var
   EDGELINE_FLUX_IN=EDGELINE_FLUX_in+EDGELINE_FLUX
! write(*,*)'Done Calcuating the XFLUX out'
!#endif

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
	  
   !  CALL  WQM_OUTPUT(NTHIS_OLD, NTSTN_OLD, NTHIS,NTSTN,CTR_HIS,CTR_STN)

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
	
!       
!******* Restarts
!

! KURT GLAESEMANN - restarts
!        IF (RESTART_OUT) THEN
!             IF (JDAY >= RSOD(RSODP) .OR. END_RUN) THEN

! # if defined (MULTIPROCESSOR)
!                 if(PAR) then
                
!                     DO J=1,NCP
!                         call gather(0,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,C1(0,1,J),C2_GL(1,1,J))
!                     ENDDO
!                     IF (SEDIMENT_CALC) THEN
					
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CTEMP    ,CTEMP_GL)
!                         call gather(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPOP     ,CPOP_GL)
!                         call gather(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPON     ,CPON_GL)
!                         call gather(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,CPOC     ,CPOC_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CPOS     ,CPOS_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,PO4T2TM1S,PO4T2TM1S_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NH4T2TM1S,NH4T2TM1S_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,NO3T2TM1S,NO3T2TM1S_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,HST2TM1S ,HST2TM1S_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CH4T2TM1S,CH4T2TM1S_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,CH41TM1S ,CH41TM1S_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SO4T2TM1S,SO4T2TM1S_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SIT2TM1S ,SIT2TM1S_GL)
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,BENSTRTM1S ,BENSTRTM1S_GL)
!                         !call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,SODTM1S, SODTM1S_GL)
!                         IF(BALGAE_CALC)THEN
!                             call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,BBM      ,BBM_GL)
!                         ENDIF
!                     ENDIF
!                     IF(DFEEDER) THEN
!                         call gather(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,DFEEDM1S ,DFEEDM1S_GL)
!                     ENDIF
!                     IF(SFEEDER) THEN
!                         call gather(1,MTLOC,MLOC,MGL,NSSFP,MYID,NPROCS,NMAP,SFEED,SFEED_GL)
!                     ENDIF
                    
!                 else
! #endif
!                     DO J=1,NCP
!                         C2_GL(1:MGL,1:KBM1,J)  = C2(1:MGL,1:KBM1,J)
!                     ENDDO
!                     IF (SEDIMENT_CALC) THEN
!                         CTEMP_GL(1:MGL)     = CTEMP(1:MGL)
!                         CPOP_GL(1:MGL,1:3)  = CPOP(1:MGL,1:3)
!                         CPON_GL(1:MGL,1:3)  = CPON(1:MGL,1:3)
!                         CPOC_GL(1:MGL,1:3)  = CPOC(1:MGL,1:3)
!                         CPOS_GL(1:MGL)      = CPOS(1:MGL)
!                         PO4T2TM1S_GL(1:MGL) = PO4T2TM1S(1:MGL)
!                         NH4T2TM1S_GL(1:MGL) = NH4T2TM1S(1:MGL)
!                         NO3T2TM1S_GL(1:MGL) = NO3T2TM1S(1:MGL)
!                         HST2TM1S_GL(1:MGL)  = HST2TM1S(1:MGL)
!                         CH4T2TM1S_GL(1:MGL) = CH4T2TM1S(1:MGL)
!                         CH41TM1S_GL(1:MGL)  = CH41TM1S(1:MGL)
!                         SO4T2TM1S_GL(1:MGL) = SO4T2TM1S(1:MGL)
!                         SIT2TM1S_GL(1:MGL)  = SIT2TM1S(1:MGL)
!                         BENSTRTM1S_GL(1:MGL)  = BENSTRTM1S(1:MGL)
!                         !SODTM1S_GL(1:MGL)     = SODTM1S(1:MGL)
!                         !changed below for consistency RGL/Kurt
!                         IF(BALGAE_CALC)THEN
!                             BBM_GL(1:MGL) = BBM(1:MGL)
!                         ENDIF
                            
!                     ENDIF
!                     IF(DFEEDER) THEN
!                         DFEEDM1S_GL(1:MGL) = DFEEDM1S(1:MGL)
!                     ENDIF
!                     IF(SFEEDER) THEN
!                         SFEED_GL(1:MGL,1:NSPECIES) = SFEED(1:MGL,1:NSPECIES)
!                     ENDIF
! # if defined (MULTIPROCESSOR)
!                 endif
! #endif
!                 RSODP  = RSODP+1
!                 IF(MSR) THEN
!                     LENGTH = INDEX(RSOFN,'.')
!                     IJDAY  = AINT(JDAY)
!                     IF (IJDAY < 10) THEN
!                         WRITE (EXT1,'(I1)') IJDAY
!                         FILENAME = RSOFN(1:LENGTH)//EXT1
!                     ELSE IF (IJDAY < 100) THEN
!                         WRITE (EXT2,'(I2)') IJDAY
!                         FILENAME = RSOFN(1:LENGTH)//EXT2
!                     ELSE IF (IJDAY < 1000) THEN
!                         WRITE (EXT3,'(I3)') IJDAY
!                         FILENAME = RSOFN(1:LENGTH)//EXT3
!                     ELSE
!                         WRITE (EXT4,'(I4)') IJDAY
!                         FILENAME = RSOFN(1:LENGTH)//EXT4
!                     ENDIF

!                     OPEN  (RSO,FILE=FILENAME, FORM='UNFORMATTED')
!                     WRITE (RSO) TITLE
!                     WRITE (RSO) (((C2_GL(I,K,J),I=1,MGL ),K=1,KBM1),J=1,NCP)
        
!                     IF (SEDIMENT_CALC)THEN
                    
!                             WRITE (RSO) (CTEMP_GL(B),B=1,MGL ),                                &
!                                 ((CPOP_GL(B,J),B=1,MGL ),J=1,3),                              &
!                                 ((CPON_GL(B,J),B=1,MGL ),J=1,3),                              &
!                                 ((CPOC_GL(B,J),B=1,MGL ),J=1,3),                              &
!                                 (CPOS_GL(B),B=1,MGL ),                                        &
!                                 (PO4T2TM1S_GL(B),B=1,MGL ), (NH4T2TM1S_GL(B),B=1,MGL ),       &
!                                 (NO3T2TM1S_GL(B),B=1,MGL ), (HST2TM1S_GL(B),B=1,MGL ),        &
!                                 (CH4T2TM1S_GL(B),B=1,MGL ), (CH41TM1S_GL(B),B=1,MGL ),        &
!                                 (SO4T2TM1S_GL(B),B=1,MGL ), (SIT2TM1S_GL(B),B=1,MGL ),        &
!                                 (BENSTRTM1S_GL(B),B=1,MGL), (SODTM1S_GL(B),B=1,MGL)
                        
!                                 IF(BALGAE_CALC)THEN
!                                     WRITE(RSO)(BBM_GL(B),B=1,MGL)
!                                 ENDIF
!                     ENDIF
                    
!                     IF (DFEEDER) WRITE (RSO) (DFEEDM1S_GL(B),B=1,MGL )
!                     IF (SFEEDER) WRITE (RSO) ((SFEED_GL(B,K),B=1,MGL ),K=1,NSPECIES)

!                     CLOSE (RSO)
!           ENDIF
!        ENDIF
!      ENDIF  !END RESTART_OUT
	 

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

#if defined(ZOOPLANKTON)   
     CALL ZOOP_DEALLOC !WLong replaced the above by ZOOP_DEALLOC()
#endif

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
            NBN ,       &   !              !!NUMBER OF BOUNDARY NODES
            BN_LOC,     &    !            !!LOCAL IDENTITY OF BOUNDARY NODES
            BN_MLT,     &    !            !!MULTIPLICITY OF BOUNDARY NODES
            NC,            &    !
            BNC,         &    !
             NODE_MATCH_ALL, & !
             EXCHANGE_ALL !, & !
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

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
   

     CALL EXCHANGE_ALL(NC,MTLOC,KBM1,NCP,MYID,NPROCS,CCC)

     CALL NODE_MATCH_ALL(1,NBN,BN_MLT,BN_LOC,BNC,MTLOC,KBM1,NCP,MYID,NPROCS,CCC)
   ENDIF
   
#  endif
   

   
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
          PI    &        !
         ,ART1   &
         ,Z    &            !SIGMA COORDINATE VALUE 
         ,ZZ    &        !INTRA LEVEL SIGMA VALUE
         ,DZ    &        !DELTA-SIGMA VALUE
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
         ,WTST    &        !Vertical velocity in sigma from PREVIOUS TIMESTEP
         ,UARD_OBCNT    &!tykim previous time step
         ,XFLUX_OBCT    &!tykim previous time step
         ,DTFAT    &        !tykim
         ,WTS    &        !VERTICAL VELOCITY IN SIGMA SYSTEM
         ,UARD_OBCN    &    ! tykim 
         ,XFLUX_OBC    &    ! tykim 
         ,KH    &        !TURBULENT DIFFUSIVITY
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
         ,THOUR    &        !
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
         ,HYDRO_GEOM_DEALLOC &    !
         ,HYDRO_DEALLOC
        
    USE MOD_BCMAP, ONLY :    &
            BCMAP_DEALLOC

       !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
    USE MOD_CONTROL, ONLY :         &
            SERIAL          &           !!TRUE IF SINGLE PROCESSOR
            ,MSR            &           !!TRUE IF MASTER PROCESSOR (MYID==1)
            ,PAR            &           !!TRUE IF MULTIPROCESSOR RUN
            ,CASENAME      !    &           !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
            
       USE MOD_FILEINFO, ONLY : &
                DIA                &
               ,BFO            & 
                ,KFL            & 
               ,RSO            &
                ,SNP            &
                ,PLT            &
                ,APL             &
                ,TFL             &
                ,SFO            &
                ,ICO            &
                ,MBL            &
                ,UNIT_STN        &
                ,UNIT_HIS        &            
                ,CNAME            &
                ,INIT_FILE_INFO

    
    USE MOD_SF, ONLY:             &     !
                SFEEDER,         &    !
                NSPECIES,        &    !         
                SFEED,            &!(:,:),
                SFEED_GL,        &!(:,:),   
                SUSPFEED,        &!        
                SF_DEALLOC!,    &!

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
                LPOC,            &!
                RPOC,            &!
                NH4,            &!
                NO3,            &!
                UREA,            &!
                PO4,            &!
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

                ATMLDON,        &!
                ATMLDOP,        &!
                ATMNH4,         &!
                ATMNO3,         &!
                ATMPO4,         &!
                ATMRDON,        &!
                ATMRDOP,        &!

                DLWCMC,         &!
                DLWCMN,         &!
                DLWCMP,         &!
                IWCMC,          &!
                IWCMN,          &!
                IWCMP,          &!
                PRECIP,         &!
                SFATOT,         &!               
                ANDC,           &!
       
                BENNH4,         &!
                BENNO3,         &!
                BENSA,          &!                
                BENPO4,         &!               
                BNDN,           &!           
                ! DLALGC,         &!
  
                DTMP1,          &!
 
                INFLOW,         &!
                JCB,            &!
                PCFWS,          &!
                PPFWS,          &!
                ! RESPC,          &!
                BNDC,           &!
 
                BNDP,           &!
                DLWCSSMN,               &!
                DLWCSSMP,               &!

                S1FLXN,         &!
                S1FLXP,         &!
                S2FLXN,         &!
                S3FLXN,         &!
                BENDEN,         &!

                DFTMP,          &!
 
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
                ! DCOD,           &!
                ! DCOD_GL,        &!
                ! DDOC,           &!
                ! DDOC_GL,        &!
                IAVG_GL,        &!
				KDLAMB_int_GL,   &!B Clark may 2018
				SSI_GL,  & ! B Clark June 2018				
       RHO1_GL,RHO1,ShearStress,ShearStress_GL,  &  ! B Clark June 2018
                IK1,            &!
                IK1_GL,         &!
                IK2,            &!
                IK2_GL,         &!
                ! NITRIF,         &!
                ! NITRIF_GL,              &!
                PR1,            &!
                PR1_GL,         &!
                PR2,            &!
                PR2_GL,         &!

                COD_GL,         &!
                COURMX,         &!
                DIAGNOSTICS,    &!
                ! DOLZ_GL,        &!
                ! DOPR,           &!
                ! DOPR_GL,        &!
                ! DORALG,         &!
                ! DORALG_GL,      &!
                ! DOSZ_GL,        &!
                FTCOD,          &!
                FTCOD_GL,       &!
                REAERDO,        &!
                REAERDO_GL,     &!
                RESP,           &!

                PNFWS,          &!                
                CTEMP,             &!
                DIFFMX,            &!
                B,              &!
                BB,             &!                
                IJDAY,          &!
                LENGTH,         &!
                SAV_LOADS,      &!
                WQM_DEALLOC,     &!
				ALGDON_OUT,    &!Algae DON production
				ALGDON_GL,    &! Global Algae DON production							
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
					  DTWCDOC1,DTWCDOC2,DTWCDOC3,  &
					  DTWNCDOC1,DTWNCDOC2,DTWNCDOC3,  &                      
					  DTWCDON1_GL,DTWCDON2_GL,DTWCDON3_GL,  &
                      DTWNCDON1_GL,DTWNCDON2_GL,DTWNCDON3_GL, &
                      DTWCDON1,DTWCDON2,DTWCDON3,  &
					  DTWNCDON1,DTWNCDON2,DTWNCDON3,  & 
					  DTNH4_GL,DTNO3_GL,DTLPON_GL,DTLPOC_GL,DTRPON_GL,DTRPOC_GL, &
					  PHOTOCHEM_DOC_ACCUM_OUT, &
					  	DTWCDOC1_out,DTWNCDOC1_out,DTWCDOC2_out,DTWNCDOC2_out,   &
						DTWCDOC3_out,DTWNCDOC3_out,   & 
						DTWCDON1_out,DTWNCDON1_out,DTWCDON2_out,DTWNCDON2_out,   &
						DTWCDON3_out,DTWNCDON3_out  
				  
    USE MOD_OWQ, ONLY:                  & 
                        INTKE,            &!
                        NXKEI,            &!
                        IAVG,            &!
                        OWQ_DEALLOC,    &!
                        LGHT_ATTN, &
						KDLAMB_int, & ! Added by B Clark May 2018 to output integrated kd for PAR
						SSI   ! B Clark June 2018
                        
    USE MOD_WQMINIT, ONLY :            &        !
            RSOD,                &        !    
            RSODP,                 &        !
            RSOFN,                 &        !
            DLTD,                &        !
            DLTDP,               &        !
            DIAD,                 &        !            
            DIADP,                &        !
            DIAF,                 &        !
            COURFS,            &        !
            COURBS,            &        !
            F,                    &        !            
            SB,                    &        !
            DIFFFS,                &        !
            DIFFBS,              &        !
            NXTVD,               &        !
            NXDIA,                &        !
            TM1,                &        !
            TM2,                &        !
            EXT1,                 &!*1,   !    
            EXT2,                &!*2,    !    
            EXT3,                &!*3,    !    
            EXT4,                &!*4,    !    
            FILENAME,            &        !
            RESTART_OUT,        &        !
            END_RUN,            &        !
            VOLUME_BALANCE,     &        !
            ICOND_OUT,          &        !
            AUTO_STEPPING,      &                   
            NEW_VOLUMES,          &        !
            ALGAE_CALC,        &        !
            ZOO_CALC,           &        !
            CARBON_CALC,       &        !
            NITROGEN_CALC,     &        !
            PHOSPHORUS_CALC,   &        !
            COD_CALC,            &        !
            OXYGEN_CALC,        &        !
            SILICA_CALC,        &        !                                
            LEFT_FLOWB,         &!(NHQP),        !        
            IFLOWP,            &    !(NBCP),        !
            ! C1MIN,              &                !
            ! C1MAX,                &                !
            DOVDAYS,            &    !(:,:,:)    !
            ACC,                &    !(NCP)        !
            DLTVAL,                &!
            DLTMAX,               &!
            DLTFTN,             &!
             WQMINIT_DEALLOC
             
   USE MOD_ALGAL, ONLY :         &
            ALGAE    
#if  defined(ZOOPLANKTON) 
   USE MOD_ZOOP, ONLY :            &
                     ANCSZ,      &!microzooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     ANCLZ,      &!mesozooplankton nitrogen to carbon ratio (gN/gC)  (~ 0.2)
                     APCSZ,      &!microzooplankton phosphorus to carbon ratio (gP/gC) (~0.02)
                     APCLZ,      &!mesozooplankton phosphorus to carbon ratio (gP/gC) (~0.02) 
                     CLSZ,        & 
                     CLLZ,        &
                      B1SZ,        &!
                      B2SZ,        &!
                     B1LZ,        &!
                     B2LZ,        &!
                     DOSZ,        &!
                     DOLZ,        &!

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
                     SZ,        &!
                     LZ,        &!
                     ZOOP_DEALLOC
#endif                     
    USE MOD_KIN, ONLY :                &
                        TEMPER,        &
                        SOLIDS,        &
#if defined(ZOOPLANKTON)
                         ZOOPL,        &
#endif
                        PARTICULATE_CARBON,        &
                        NITROG,        &
#if defined(PHOSPHORUS)
                         PHOSPH,        &
#endif
                        CODMND,        &
                        OXYGEN,        &
                        SILICA,        &
                        BEN_FLUX

   USE MOD_SED, ONLY :         &!
        STEADY_STATE_SED_G3,&!
        STEADY_STATE_SED_IC,&!
        DLTS,                &!
        IERR,                &!
    
        ISEDMN,               &!
        ISEDMP,               &!
        ISEDMC,                &!
        SEDMN,              &!total sediment organic nitrogen  (kgN)
        SEDMP,              &!total sediment organic phosphorus (kgP)
        SEDMC,                &!total sediment organic carbon  (kgC)
     
        BENSTRTM1S,           &!                             
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
        SED_INIT2,            &!
        SED_DEALLOC,        &!
        SED_CALC,            &!
        SED_DIAGENESIS_G3!,            &!    

   USE MOD_SED_DF_EXCHANGE_VARS, ONLY:     SED_DF_EXCHANGE_ALLOC,    &!
                                        SED_DF_EXCHANGE_DEALLOC
   
   USE MOD_SED_SF_EXCHANGE_VARS, ONLY:     SED_SF_EXCHANGE_ALLOC,    &!
                                        SED_SF_EXCHANGE_DEALLOC
   
   USE MOD_SED_SAV_EXCHANGE_VARS, ONLY: SED_SAV_EXCHANGE_ALLOC,    &!
                                        SED_SAV_EXCHANGE_DEALLOC
   
   USE MOD_BA, ONLY : BALGAE_CALC, BBM,BBM_GL, BA_DEALLOC
   
   USE MOD_OBCS, ONLY:  &

        SETUP_OBC ,        &!
        OBCS_DEALLOC    
   
   USE MOD_BCS, ONLY : &
  
        BCS_DEALLOC

#if defined(SEDIMENT_DOM)			
			USE MOD_SED_DOM, ONLY : SED_DOM, &        ! B Clark DOM Flag
									SED_DOM_OUTPUT,JDOM_OUT_GL
									
#endif
									
#if defined(PHOTODEGRADATION)
	USE WC_DOM, ONLY : PHOTOCHEM_DOC_ACCUM, PHOTOCHEM_DOC_ACCUM_GL, PHOTOCHEM_DOC_ACCUM_OUT								
									
#endif
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY:   &              EC,            &    !
            NC,            &    !
            BNC,         &    !
            EMAP,        &    !
            NMAP,         &    !
             NODE_MATCH_ALL, &  !
             EXCHANGE_ALL, & !
             GATHER , & !
             INIT_MPI_ENV,&!
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
            DF_DEALLOC!,        &!

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


!WLong: Now starting and ending times of history and station outputs are controlled
!       in wqm_linkage.in file 

     NTHIS = FLOOR((JDAY-t_his_start)/t_his_dlt) 
     NTSTN = FLOOR((JDAY-t_stn_start)/t_stn_dlt)

IF((NTHIS>NTHIS_OLD).AND.(JDAY >= t_his_start) .AND.(JDAY <= t_his_end)) THEN
	 
#if defined (MULTIPROCESSOR)

    if(PAR) then
 
! Wen Long added the following for debugging DOXG

#if defined(BUDGET)	
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
   
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGDON_OUT(1,1), ALGDON_GL(1,1))																											  
			CALL GATHER(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,ALGPON_OUT(1,1), ALGPON_GL(1,1))	
																						  
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
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTNH4_OUT,DTNH4_GL(1,1))!CDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTNO3_OUT,DTNO3_GL(1,1)) !NCDON1
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPOC_OUT,DTLPOC_GL(1,1))!CDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTRPOC_OUT,DTRPOC_GL(1,1)) !NCDON2
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPON_OUT,DTLPON_GL(1,1)) !CDON3
		call gather(1,MTLOC,MLOC,MGL,KBM1,MYID,NPROCS,NMAP,DTLPON_OUT,DTRPON_GL(1,1))!NCDON3	

    IF(SEDIMENT_CALC)THEN !LBcomment
	         CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,JPOC_out(1,1)  ,    JPOC_GL(1,1))
            CALL GATHER(1,MTLOC,MLOC,MGL,3,MYID,NPROCS,NMAP,JPON_out(1,1)  ,    JPON_GL(1,1))   
			   CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JNH4_out(1)  ,        JNH4_GL(1))               
            CALL GATHER(1,MTLOC,MLOC,MGL,1,MYID,NPROCS,NMAP,JNO3_out(1)  ,        JNO3_GL(1))   
    ENDIF  
#endif  
   !sed calc

    endif! par
#endif

#if defined(BUDGET)
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

         ENDIF !ctr_his

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

                      FLUSH(UNIT_HIS)
             ENDSELECT
           ENDDO
         ENDIF !END PRINT_HIS SWITCH
       ENDIF  !END MSR SWITCH
#endif
    
    ENDIF !END AMOD(JDAY  ... ) SWITCH
END SUBROUTINE WQM_OUTPUT
