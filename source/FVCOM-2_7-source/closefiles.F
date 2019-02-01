 
!==============================================================================!
!   CLOSE UP FILES                                                             !
!==============================================================================!

   SUBROUTINE CLOSEFILES
   USE ALL_VARS
   IMPLICIT NONE
   REAL(SP) ::  COMT(80)

   CLOSE(IOPRT)
   CLOSE(IOPLT)
   CLOSE(IOTSR)
   
   IF(M_TYPE=='non-uniform')THEN
     CLOSE(INHFX) ; CLOSE(INWND)
     IF(EVP_FLAG) CLOSE(INEVP) 
   END IF
     
   IF(S_TYPE == 'julian') CLOSE(INJUL)
     

   RETURN
   END SUBROUTINE CLOSEFILES
!==============================================================================!




