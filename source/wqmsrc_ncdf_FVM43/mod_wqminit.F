        MODULE MOD_WQMINIT
        
         USE MOD_PREC, ONLY: SP
         USE MOD_SIZES, ONLY :        &    !
            NCP,            &  !
            !NQFP,           &  !
            NHQP ,           &  !
            !NS1P,           &  !
            !NS2P,           &  !
            !NS3P,           &  !
            NBCP,           &  !
            NDP,            &  !
            !NFLP,           &  !
            NOIP ,           &  !
            NSSFP!,             &  !
            !MGL,            &  !
            !NGL,            &  !
            !OBCGL,          &  !
            !NOBTY  

        USE MOD_LIMS, ONLY:     &!
                KBM1,           &!
                MTLOC
                
        IMPLICIT NONE
        
            SAVE
        !***** Variable declarations
            
            INTEGER ::     RSODP,         &
                        DLTDP,       &
                        SNPDP,        &
                        TFLDP,        &
                        PLTDP,        &
                        APLDP,       &
                        OPLDP,        &
                        DIADP        
            INTEGER :: COURB,        &
                        COURFS,        &
                        COURBS,        &
                        F,            &            
                        SB,            &
                        S1LNMAX,    &
                        S2LNMAX,    &
                        S3LNMAX,    &
                        DIFFFS,        &
                        DIFFBS      
             
            REAL(SP)    ::     NXSNP,         &
                        NXPLT,        &
                        NXAPL,        &
                        NXTFL,        &
                        NXKFL,        &
                        NXTVD,       &
                        NXOPL,        &
                        NXMBL,        &
                        NXDIA
                        
!            REAL(SP)::  NO3T2I,        &!    WLong moved this to mod_sed.F
!                        NH4T2I
                        
            REAL(SP)::    LEAFI,        &
                
                        STEMI,        &
                        ROOTI,         &
                        TUBERI    

             REAL(SP)    :: MAXDLT,        &    !
                        MXDLT!,        &    !
!                        INTKE!,        &    !Wen Long moved to mod_owq.F
!                        INITKE            !Wen Long moved to mod_owq.F

             REAL(SP) :: TM1,        &    !
                        TM2,        &    !
                        TM3,        &    !
                        MINSTEP!,    &    !
!                        NXKEI            !moved to mod_owq.F
                        
             CHARACTER(LEN=3) :: SNPC,        &
                                 RSIC,         &
                                 RSOC,        &
                                 BCC,        &
                                 S1C,       &
                                 S2C,        &
                                 MDC*3,        &
                                 PLTC,        &
                                 FLC,       &
                                 MBLC,        &
                                 BFC,        &
                                 VBC,        &
                                 QPLTC,        &
                                 XYDFC,        &
                                 XYDFU,      &
                                 ZDFC,        &
                                 ICOC,        &
                                 ATMC,        &
                                 SAVLC,        &
                                 SEDC,      &
                                 AUTOC,        &
                                 SPLTC,        &
                                 TFLC,        &
                                 DIAC,        &
                                 STLC,      &
                                 APLTC,        &
                                 KFLC,        &
                                 OPLC,        &
                                 BFOC,        &  
                                 BAOC,      & !WLong benthic algae ouytput flag
                                 DFOC,        & !WLong, deposition feeder output flag
                                 S3C,       &
                                 SAVMC,        &
                                 SAVPLTC,    &
                                 DFLC,        &
                                 SFLC,        &
                                 SFLOXC,    & !WLong, hypoxia effects on suspension feeder
                                 DFLOXC          !WLong, hypoxia effects on deposition feeder
                                 
             CHARACTER ::        EXT1*1,       &
                                EXT2*2,        &
                                EXT3*3,        &
                                EXT4*4,        &
                                EXT_R*1,    &
                                EXT_R2*2    
             
             CHARACTER(LEN=8)::    SLC,        &
                                HYDC,        &
                                BNDTC,        &
                                CONSC,        &
                                ICIC   
                                
             CHARACTER(LEN=8)::    TIMEA,        &
                                TIMEB,        &
                                SPVARM,        &
                                PRINTM
                                
             CHARACTER(LEN=72)::TITLE(6),    &
                                OLDTITLE(6),&
                                FILENAME
                                
             CHARACTER(LEN=100) :: FILENAME2
             
             CHARACTER(LEN=72)::MAPFN,        &
                                GEOFN,        &
                                ICIFN,        &
                                RSIFN,      &
                                AGRFN,        &
                                STLFN,        &
                                MRLFN,        &
                                KFLFN,        &
                                ICOFN
            
            CHARACTER(LEN=72):: SNPFN,        &
                                PLTFN,        & 
                                APLFN,         &
                                DIAFN,         &
                                TFLFN,      &
                                RSOFN,         &
                                OPLFN,         &
                                MBLFN,         &
                                SFIFN,         &
                                SFOFN

            LOGICAL ::          RESTART_OUT,          &
                                SNAPSHOTS,            &
                                END_RUN,            &
                                MODIFY_ICONC,       &  
                                VOLUME_BALANCE,     &
                                QUICKEST,            &
                                UPWIND,                &
                                ICOND_OUT,          &
                                UNI_ICON_IN,           &
                                BIN_ICON_IN,        &
                                RES_ICON_IN,        &
                                AUTO_STEPPING,       &
                                STOP_RUN,           &
                                PLOTS,                &
                                OXYGEN_PLOTS,        &
                                ZOO_CALC,           &
                                NEW_VOLUMES,          &
                                RESTART_IN
                        
             LOGICAL :: TEMPERATURE_CALC,             &
                        ALGAE_CALC,                 &
                        CARBON_CALC,                   &
                        NITROGEN_CALC,                 &
                        PHOSPHORUS_CALC,             &
                        COD_CALC,                    &
                        OXYGEN_CALC,                &
                        SILICA_CALC
                        
             LOGICAL ::    LEFT_FLOWB(NHQP),            &
                        RIGHT_FLOWB(NHQP),          &
                        LEFTM1_BOUNDARY(NHQP),         &
                        RIGHTP1_BOUNDARY(NHQP)      

        !***** Dimension declarations

             REAL(SP),DIMENSION(NHQP,2) ::     DEN1,        &
                                        DEN2,        &
                                        DEN3,       &
                                        TP1,         &
                                        TP2,        &
                                        TP3,        &
                                        T2,            &
                                        SF2
                                        
             REAL(SP),DIMENSION(NHQP) ::     T1,            &
                                        T3,            &
                                        SF1
                                        
             REAL(SP),DIMENSION(NHQP) ::     TERM1,        &
                                        TERM2,         &
                                        TERM3,      &
                                        GRAD1,         &
                                        GRAD2,        &
                                        GRAD3,      &
                                        COUR
            
             REAL(SP) ::     GRAD(NHQP,3),        &
                        TERM(NHQP,3),        &
                        IFLOWP(NBCP),        &
                        IBT(NHQP),          &
                        COURBSV(NHQP),        &
                        COURVSV(NHQP)

            EQUIVALENCE (GRAD1,GRAD(1,1))
            EQUIVALENCE    (GRAD2,GRAD(1,2))
            EQUIVALENCE (GRAD3,GRAD(1,3))
            EQUIVALENCE (TERM1,TERM(1,1))
            EQUIVALENCE    (TERM2,TERM(1,2))
            EQUIVALENCE (TERM3,TERM(1,3))

        !     REAL(SP) ::     CR1(0:NLP,NSBP), &
        !                    CR2(0:NLP,NSBP), & 
        !                    GAMFM1(NLP,NSBP),&
        !                   GAMF(NLP,NSBP),  &
        !                    WDILB(NLP,NSBP), & 
        !                    WDIB(NLP,NSBP),  & 
        !                   WAILB(NLP,NSBP), &
        !                    WAIB(NLP,NSBP)
        !     REAL(SP) :: DIFFS(NQFP)
        !     REAL(SP),DIMENSION(NLP,NSBP) :: AT,    BT,     CT
        !     REAL(SP)),DIMENSION(NLP,NCP) :: DT,     GAMMA,   BTA
        !     REAL(SP) :: VT(NBP,NCP)
        !     REAL(SP),DIMENSION(NBP,NCP) :: C1MIN,   C1MAX
            
             REAL(SP),ALLOCATABLE,DIMENSION(:,:,:) ::     C1MIN,  &
                                                    C1MAX

        !    REAL(SP),DIMENSION(NBP) :: V2VB, DLBV
             REAL(SP),ALLOCATABLE ::     DOVDAYS(:,:,:)
             REAL(SP) ::                 OINT(NOIP)
             CHARACTER(LEN=3) ::     ACC(NCP)
             
             REAL(SP) :: CIC(NCP),     &
                        COUT(NCP)!
                        
!            REAL(SP)::     CPOPI(3),    &        !WLong moved to mod_sed.F
!                        CPONI(3),     &
!                        CPOCI(3)
                     
            REAL(SP) ::  SFEEDI(NSSFP) !WLong changed dimension size from 10 to NSSFP
                     
            REAL(SP) :: DLTD(NDP),        &
                     SNPD(NDP),        &
                     RSOD(NDP),     &
                     SNPF(NDP),     &
                     DLTVAL(NDP),    &
                     DLTMAX(NDP),   &
                     DLTFTN(NDP),     &
                     PLTD(NDP),        &
                     PLTF(NDP),     &
                     APLTD(NDP),    &
                     APLF(NDP),        &
                     TFLD(NDP),     &
                     TFLF(NDP),        &
                     KFLD(NDP),        & 
                     KFLF(NDP),     &
                     OPLD(NDP),        &
                     OPLF(NDP),        &
                     DIAD(NDP),     &
                     DIAF(NDP),     &
                     MBLD(NDP),     &
                     MBLF(NDP)
             
        !***** Data declarations
             REAL(SP),PARAMETER ::     DOCLIT  = 0.0,     &
                                LPOCLIT = 0.0,  &
                                POCLIT = 0.0,      &
                                PBSLIT  = 0.0,  &
                                DSILLIT = 0.0

        CONTAINS
        
           !subroutines:
            !subroutine     WQMINIT_ALLOC()
        
            SUBROUTINE WQMINIT_ALLOC
                IMPLICIT NONE
                
                ALLOCATE(C1MIN(MTLOC,KBM1,NCP)); C1MIN = 0.0   
                ALLOCATE(C1MAX(MTLOC,KBM1,NCP)); C1MAX = 0.0
                ALLOCATE(DOVDAYS(MTLOC,KBM1,NOIP)); DOVDAYS = 0.0
                PHOSPHORUS_CALC = .FALSE.
                RETURN
            END SUBROUTINE WQMINIT_ALLOC

            SUBROUTINE WQMINIT_DEALLOC
                IMPLICIT NONE
                IF(ALLOCATED(C1MIN))DEALLOCATE(C1MIN)
                IF(ALLOCATED(C1MAX))DEALLOCATE(C1MAX)
                IF(ALLOCATED(DOVDAYS))DEALLOCATE(DOVDAYS)
                RETURN
            END SUBROUTINE WQMINIT_DEALLOC
            
        END MODULE MOD_WQMINIT

