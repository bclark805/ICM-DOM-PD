!==============================================================================|
!   Set Boundary Conditions on Temperature and Salinity                        |
!    NCON2 = 1:  SET CONDITIONS SPECIFIC TO TEMPERATURE                        |
!    NCON2 = 2:  SET CONDITIONS SPECIFIC TO SALINITY                           |
!==============================================================================|

   SUBROUTINE BCOND_TS(NCON2)     

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
#  if defined (TS_OBC)
   USE MOD_TSOBC
#  endif

   IMPLICIT NONE
   REAL(SP) :: S2D,S2D_NEXT,S2D_OBC,T2D,T2D_NEXT,T2D_OBC,XFLUX2D,TMP,RAMP_TS

#  if defined (TS_OBC)
   REAL(SP), ALLOCATABLE :: TOBC_TMP(:,:),SOBC_TMP(:,:)
   REAL(SP) ::FACT,UFACT
   INTEGER  ::L1,L2,IERR
#  endif 

   INTEGER  :: I,J,K,J1,J11,J22,NCON2
   REAL(SP), ALLOCATABLE :: TTMP(:,:),STMP(:,:)

   REAL(SP) ::TMAX,TMIN,SMAX,SMIN

!------------------------------------------------------------------------------|


!
!--SET CONDITIONS FOR FRESH WATER INFLOW---------------------------------------|
!
   IF(POINT_ST_TYPE == 'specified') THEN
     IF(NUMQBC > 0) THEN
       IF(INFLOW_TYPE == 'node') THEN
         DO I=1,NUMQBC
           J11=INODEQ(I)
           DO K=1,KBM1
             TF1(J11,K)=TDIS(I)
             SF1(J11,K)=SDIS(I)
           END DO
         END DO
       ELSE IF(INFLOW_TYPE == 'edge') THEN
         DO I=1,NUMQBC
           J11=N_ICELLQ(I,1)
           J22=N_ICELLQ(I,2)
           DO K=1,KBM1
             TF1(J11,K)=TDIS(I)
             SF1(J11,K)=SDIS(I)
             TF1(J22,K)=TDIS(I)
             SF1(J22,K)=SDIS(I)
           END DO
         END DO
       END IF
     END IF
   END IF

       
   IF(IOBCN > 0) THEN
#    if defined (TS_OBC)
     ALLOCATE(TOBC_TMP(1:IOBCN,1:KBM1))
     ALLOCATE(SOBC_TMP(1:IOBCN,1:KBM1))
     CALL BRACKET(TSO_TM,THOUR,L1,L2,FACT,UFACT,IERR)
     IF(IERR==-1)THEN
       TOBC_TMP = 0.0_SP
       SOBC_TMP  = 0.0_SP
     ELSE
       TOBC_TMP(1:IOBCN,1:KBM1) = UFACT*TEMPOBC(1:IOBCN,1:KBM1,L1)  +&
                                 FACT*TEMPOBC(1:IOBCN,1:KBM1,L2)
       SOBC_TMP(1:IOBCN,1:KBM1)  = UFACT*SALTOBC(1:IOBCN,1:KBM1,L1)  +&
                                 FACT*SALTOBC(1:IOBCN,1:KBM1,L2)
     END IF
#    endif
!
!  SET TEMPERATURE CONDITIONS ON OUTER BOUNDARY
!
   RAMP_TS = TANH(FLOAT(IINT)/FLOAT(IRAMP+1))
   IF(NCON2 == 1) THEN
   ALLOCATE(TTMP(IOBCN,KBM1));  TTMP = 0.0_SP
     DO I=1,IOBCN
       J=I_OBC_N(I)
       J1=NEXT_OBC(I)
       T2D=0.0_SP
       T2D_NEXT=0.0_SP
       XFLUX2D=0.0_SP
       DO K=1,KBM1
         T2D=T2D+T1(J,K)*DZ(J,K)
         T2D_NEXT=T2D_NEXT+TF1(J1,K)*DZ(J1,K)
         XFLUX2D=XFLUX2D+XFLUX_OBC(I,K)           !*DZ(J,K)
       END DO
  
       IF(UARD_OBCN(I) > 0.0_SP) THEN
         TMP=XFLUX2D+T2D*UARD_OBCN(I)
         T2D_OBC=(T2D*DT(J)-TMP*DTI/ART1(J))/D(J)

         CALL BCOND_T_PERTURBATION(T2D_NEXT,T2D,TTMP,I,J,J1)
	 
         DO K=1,KBM1
           TF1(J,K)=T2D_OBC+TTMP(I,K)
!           TF1(J,K)=T2D_OBC+(TF1(J1,K)-T2D_NEXT)
         END DO

         DO K=1,KBM1
           TMAX = MAXVAL(T1(NBSN(J,1:NTSN(J)),K))
           TMIN = MINVAL(T1(NBSN(J,1:NTSN(J)),K))
         
           IF(K == 1)THEN
            TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K+1)))
            TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K+1)))
           ELSE IF(K == KBM1)THEN
            TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)))
            TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  & 
                   (DZ(J,K)+DZ(J,K-1)))
           ELSE
            TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)),                             &
                   (T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/           &
                   (DZ(J,K)+DZ(J,K+1)))
            TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)),                             &
                   (T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/           &
                   (DZ(J,K)+DZ(J,K+1)))
           END IF
 
           IF(TMIN-TF1(J,K) > 0.0_SP)TF1(J,K) = TMIN
           IF(TF1(J,K)-TMAX > 0.0_SP)TF1(J,K) = TMAX

         END DO

        ELSE
         DO K=1,KBM1

#          if defined (TS_OBC)
          IF(IERR.NE.-1)THEN
                TF1(J,K) = T1(J,K) - ALPHA_SERIES_OBC*RAMP_TS*(T1(J,K)-TOBC_TMP(I,K))   !T at OBC is equal to T specified by the imnput file 
              ELSE
                TF1(J,K) = T1(J,K)
	      ENDIF 	
#          else
              TF1(J,K) = T1(J,K) - ALPHA_OBC*RAMP_TS*(T1(J,K)-TEMP_OBC(I))
#          endif

         END DO
       END IF
     END DO
     DEALLOCATE(TTMP)


!
!  SET SALINITY CONDITIONS ON OUTER BOUNDARY
!
   ELSE IF(NCON2 == 2) THEN
   ALLOCATE(STMP(IOBCN,KBM1));  STMP = 0.0_SP
     DO I=1,IOBCN
       J=I_OBC_N(I)
       J1=NEXT_OBC(I)
       S2D=0.0_SP
       S2D_NEXT=0.0_SP
       XFLUX2D=0.0_SP
       DO K=1,KBM1
         S2D=S2D+S1(J,K)*DZ(J,K)
         S2D_NEXT=S2D_NEXT+SF1(J1,K)*DZ(J1,K)
         XFLUX2D=XFLUX2D+XFLUX_OBC(I,K)             !*DZ(J,K)
       END DO
 
       IF(UARD_OBCN(I) > 0.0_SP) THEN
         TMP=XFLUX2D+S2D*UARD_OBCN(I)
         S2D_OBC=(S2D*DT(J)-TMP*DTI/ART1(J))/D(J)

         CALL BCOND_S_PERTURBATION(S2D_NEXT,S2D,STMP,I,J,J1)
	 
         DO K=1,KBM1
           SF1(J,K)=S2D_OBC+STMP(I,K)  
!           SF1(J,K)=S2D_OBC+(SF1(J1,K)-S2D_NEXT)  
          END DO

         DO K=1,KBM1
           SMAX = MAXVAL(S1(NBSN(J,1:NTSN(J)),K))
           SMIN = MINVAL(S1(NBSN(J,1:NTSN(J)),K))

           IF(K == 1)THEN
            SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K+1)))
            SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K+1)))
           ELSE IF(K == KBM1)THEN
            SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)))
            SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)))
           ELSE
            SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)),                             &
                   (S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/           &
                   (DZ(J,K)+DZ(J,K+1)))
            SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)),                             &
                   (S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/           &
                   (DZ(J,K)+DZ(J,K+1)))
           END IF

           IF(SMIN-SF1(J,K) > 0.0_SP) SF1(J,K) = SMIN
           IF(SF1(J,K)-SMAX > 0.0_SP) SF1(J,K) = SMAX

         END DO
        ELSE
         DO K=1,KBM1

#          if defined (TS_OBC)
              IF(IERR.NE.-1)THEN
                SF1(J,K) = S1(J,K) - ALPHA_SERIES_OBC*RAMP_TS*(S1(J,K)-SOBC_TMP(I,K))
              ELSE
                SF1(J,K) = S1(J,K)
              ENDIF
#          else
              SF1(J,K) = S1(J,K) - ALPHA_OBC*RAMP_TS*(S1(J,K)-SALT_OBC(I))
#          endif

         END DO
       END IF
     END DO
     DEALLOCATE(STMP)
   ELSE
     PRINT*, 'NCON2 NOT IN THE LIST'
     PRINT*, 'MUST BE 1 OR 2'
     CALL PSTOP
   END IF

#  if defined (TS_OBC)
     DEALLOCATE(TOBC_TMP,SOBC_TMP)
#  endif

   END IF

!
!--SET BOUNDARY CONDITIONS-----------------------------------------------------|
!
   DO K=1,KBM1
     T(0,K)=0.0_SP
     S(0,K)=0.0_SP
   END DO

   RETURN
   END SUBROUTINE BCOND_TS
!==============================================================================|
