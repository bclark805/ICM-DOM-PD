
MODULE MOD_SED_DF_EXCHANGE_VARS

    USE MOD_PREC, ONLY: SP
    USE MOD_LIMS, ONLY: MTLOC
    IMPLICIT NONE
    SAVE
         
    REAL(SP) :: M1_SED_DF
    REAL(SP) :: M2_SED_DF
    REAL(SP), ALLOCATABLE ::    POC1TM1S_SED_DF(:), &
                                POC2TM1S_SED_DF(:)    !POC1TM1S and POC2TM1S shared between DF and sediment module
                                

                                    
    CONTAINS
                    
        SUBROUTINE SED_DF_EXCHANGE_ALLOC
        
            ALLOCATE(POC1TM1S_SED_DF(MTLOC));    POC1TM1S_SED_DF=0.0
            ALLOCATE(POC2TM1S_SED_DF(MTLOC));   POC2TM1S_SED_DF=0.0
            
        END SUBROUTINE SED_DF_EXCHANGE_ALLOC

        SUBROUTINE SED_DF_EXCHANGE_DEALLOC
        
            IF(ALLOCATED(POC1TM1S_SED_DF))DEALLOCATE(POC1TM1S_SED_DF)
            IF(ALLOCATED(POC2TM1S_SED_DF))DEALLOCATE(POC2TM1S_SED_DF)

        END SUBROUTINE SED_DF_EXCHANGE_DEALLOC
		
END

