
    MODULE MOD_SED_SAV_EXCHANGE_VARS

           USE MOD_PREC, ONLY: SP
           USE MOD_LIMS, ONLY: MTLOC
           IMPLICIT NONE
           SAVE
         
           REAL(SP), ALLOCATABLE :: NH4T2TM1S_SHARE(:)  !NH4T2TM1S to be shared between SAV and sediment module
           REAL(SP), ALLOCATABLE :: PO4T2TM1S_SHARE(:)  !PO4T2TM1S to be shared between SAV and sediment mnodule  
           REAL(SP), ALLOCATABLE :: HST2TM1S_SHARE(:)   !
           REAL(SP) ::              M2_SHARE            !kg/L
           REAL(SP) ::              PIE2HS_SHARE        !partitioning coef. of H2S in sediment layer 2 (L/kg)

    CONTAINS
    
        SUBROUTINE     SED_SAV_EXCHANGE_ALLOC
        
            ALLOCATE(NH4T2TM1S_SHARE(MTLOC));     NH4T2TM1S_SHARE    =0.0
            ALLOCATE(PO4T2TM1S_SHARE(MTLOC));     PO4T2TM1S_SHARE    =0.0
            ALLOCATE(HST2TM1S_SHARE(MTLOC));     HST2TM1S_SHARE    =0.0
            
        END SUBROUTINE SED_SAV_EXCHANGE_ALLOC

        SUBROUTINE     SED_SAV_EXCHANGE_DEALLOC
        
            IF(ALLOCATED(NH4T2TM1S_SHARE))    DEALLOCATE(NH4T2TM1S_SHARE)
            IF(ALLOCATED(PO4T2TM1S_SHARE))    DEALLOCATE(PO4T2TM1S_SHARE)
            IF(ALLOCATED(HST2TM1S_SHARE))    DEALLOCATE(HST2TM1S_SHARE)
            
        END SUBROUTINE SED_SAV_EXCHANGE_DEALLOC
    END

