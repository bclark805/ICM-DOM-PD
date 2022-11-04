
MODULE MOD_SED_SF_EXCHANGE_VARS

    USE MOD_PREC, ONLY: SP
    USE MOD_LIMS, ONLY: MTLOC
    IMPLICIT NONE
    SAVE
         
        REAL(SP), ALLOCATABLE ::JLPOC_SED_SF(:), &
                                JLPON_SED_SF(:), &
                                JLPOP_SED_SF(:), &
                                JRPOC_SED_SF(:), &                        
                                JRPON_SED_SF(:), &                
                                JRPOP_SED_SF(:), &                
                                JNH4_SED_SF(:),  &               
                                JPO4_SED_SF(:),  &               
                                SOD_SED_SF(:),   &              
                                JSA_SED_SF(:),   &              
                                JSU_SED_SF(:),   &              
                                SSI_SED_SF(:),   &
                                SU_SED_SF(:),    &     
                                SA_SED_SF(:),    &             
                                PIP_SED_SF(:)
        
        CONTAINS
        
        SUBROUTINE SED_SF_EXCHANGE_ALLOC
        
                ALLOCATE(JLPOC_SED_SF(MTLOC));  JLPOC_SED_SF=0.0
                ALLOCATE(JLPON_SED_SF(MTLOC));  JLPON_SED_SF=0.0
                ALLOCATE(JLPOP_SED_SF(MTLOC));  JLPOP_SED_SF=0.0
                ALLOCATE(JRPOC_SED_SF(MTLOC));  JRPOC_SED_SF=0.0
                ALLOCATE(JRPON_SED_SF(MTLOC));  JRPON_SED_SF=0.0
                ALLOCATE(JRPOP_SED_SF(MTLOC));  JRPOP_SED_SF=0.0

                ALLOCATE(JNH4_SED_SF(MTLOC));  JNH4_SED_SF=0.0
                ALLOCATE(JPO4_SED_SF(MTLOC));  JPO4_SED_SF=0.0
                ALLOCATE(SOD_SED_SF(MTLOC));  SOD_SED_SF=0.0
                ALLOCATE(JSA_SED_SF(MTLOC));  JSA_SED_SF=0.0
                ALLOCATE(JSU_SED_SF(MTLOC));  JSU_SED_SF=0.0
                ALLOCATE(SSI_SED_SF(MTLOC));  SSI_SED_SF=0.0

                ALLOCATE(SU_SED_SF(MTLOC));  SU_SED_SF=0.0
                ALLOCATE(SA_SED_SF(MTLOC));  SA_SED_SF=0.0
                ALLOCATE(PIP_SED_SF(MTLOC)); PIP_SED_SF=0.0


        
        END SUBROUTINE SED_SF_EXCHANGE_ALLOC

        SUBROUTINE SED_SF_EXCHANGE_DEALLOC
        
               IF(ALLOCATED(JLPOC_SED_SF))DEALLOCATE(JLPOC_SED_SF)
               IF(ALLOCATED(JLPON_SED_SF))DEALLOCATE(JLPON_SED_SF)
               IF(ALLOCATED(JLPOP_SED_SF))DEALLOCATE(JLPOP_SED_SF)
               IF(ALLOCATED(JRPOC_SED_SF))DEALLOCATE(JRPOC_SED_SF)
               IF(ALLOCATED(JRPON_SED_SF))DEALLOCATE(JRPON_SED_SF)
               IF(ALLOCATED(JRPOP_SED_SF))DEALLOCATE(JRPOP_SED_SF)

               IF(ALLOCATED(JNH4_SED_SF))DEALLOCATE(JNH4_SED_SF)
               IF(ALLOCATED(JPO4_SED_SF))DEALLOCATE(JPO4_SED_SF)
               IF(ALLOCATED(SOD_SED_SF))DEALLOCATE(SOD_SED_SF)
               IF(ALLOCATED(JSA_SED_SF))DEALLOCATE(JSA_SED_SF)
               IF(ALLOCATED(JSU_SED_SF))DEALLOCATE(JSU_SED_SF)
               IF(ALLOCATED(SSI_SED_SF))DEALLOCATE(SSI_SED_SF)
               IF(ALLOCATED(SU_SED_SF))DEALLOCATE(SU_SED_SF)
               IF(ALLOCATED(SA_SED_SF))DEALLOCATE(SA_SED_SF)
               IF(ALLOCATED(PIP_SED_SF))DEALLOCATE(PIP_SED_SF)

        END SUBROUTINE SED_SF_EXCHANGE_DEALLOC
END

