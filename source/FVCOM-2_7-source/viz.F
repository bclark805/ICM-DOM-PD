#    if defined (PV3)
!==============================================================================!
      SUBROUTINE VIZ_INIT(IOPT)
!==============================================================================!
      USE ALL_VARS
      INTEGER, INTENT(INOUT) :: IOPT
      INTEGER, PARAMETER :: NKEYS=26,NCUTS=3
      INTEGER :: FKEYS(NKEYS),IKEYS(NKEYS)
      REAL*4 :: FLIMS(2,NKEYS)
      CHARACTER(LEN=16) TKEYS(NKEYS),TCUTS(NCUTS)
      CHARACTER(LEN=32) INFILE
      CHARACTER(LEN=80) DATA(5), TTL
      REAL*4 :: VSCALE
      INTEGER:: MIRR,KNBLOCK,NCNT
      INTEGER:: ISTAT,IERR
!
      DATA TTL/'      FVCOM PV3 CLIENT                '/
      VSCALE = SUM(DLTXYC(1:NE))/FLOAT(NE) !!VECTOR LENGTH SCALE
#     if defined (MULTIPROCESSOR)
      CALL MPI_BCAST(VSCALE,1,MPI_REAL,0,MPI_COMM_WORLD,IERR)
#     endif
       NCNT = 0
! FREE SURFACE (EL) -e
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "FREE SURFACE EL"
       IKEYS(NCNT) = 101
       FLIMS(1,NCNT) = MINVAL(EL(1:NNode))
       FLIMS(2,NCNT) = MAXVAL(EL(1:NNode))
! UA -u
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "UA"
       TKEYS(NCNT) = "UA"
       IKEYS(NCNT) = 117
       FLIMS(1,NCNT) = MINVAL(UA(1:MElem))
       FLIMS(2,NCNT) = MAXVAL(UA(1:MElem))
! VA -v
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "VA"
       IKEYS(NCNT) = 118
       FLIMS(1,NCNT) = MINVAL(VA(1:MElem))
       FLIMS(2,NCNT) = MAXVAL(VA(1:MElem))
! DEPTH -b
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "BATHYMETRIC DEPTH"
       IKEYS(NCNT) = 98
       FLIMS(1,NCNT) = MINVAL(H(1:NNode))
       FLIMS(2,NCNT) = MAXVAL(H(1:NNode))
! SALINITY S
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "SALINITY"
       IKEYS(NCNT) = 83
       FLIMS(1,NCNT) = MINVAL(S1(1:NNode,1:KBM1))
       FLIMS(2,NCNT) = MAXVAL(S1(1:NNode,1:KBM1))
! TEMPERATURE T
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "TEMPERATURE"
       IKEYS(NCNT) = 84
       FLIMS(1,NCNT) = MINVAL(T1(1:NNode,1:KBM1))
       FLIMS(2,NCNT) = MAXVAL(T1(1:NNode,1:KBM1))
! X-VELOCITY U
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "U"
       IKEYS(NCNT) = 85
       FLIMS(1,NCNT) = MINVAL(U(1:MElem,1:KBM1))
       FLIMS(2,NCNT) = MAXVAL(U(1:MElem,1:KBM1))
! Y-VELOCITY V
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "V"
       IKEYS(NCNT) = 86
       FLIMS(1,NCNT) = MINVAL(V(1:MElem,1:KBM1))
       FLIMS(2,NCNT) = MAXVAL(V(1:MElem,1:KBM1))
! TURB LENGTH SCALE L
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "TURBULENT LENGTH"
       IKEYS(NCNT) = 76
       FLIMS(1,NCNT) = MINVAL(L(1:MElem,1:KBM1))
       FLIMS(2,NCNT) = MAXVAL(L(1:MElem,1:KBM1))
! TURB LENGTH SCALE Q
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "TURBULENT KE"
       IKEYS(NCNT) = 81
       FLIMS(1,NCNT) = MINVAL(Q2(1:MElem,1:KBM1))
       FLIMS(2,NCNT) = MAXVAL(Q2(1:MElem,1:KBM1))
! WTSURF h
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "SURFACE HEATING"
       IKEYS(NCNT) = 104
       FLIMS(1,NCNT) = 0.
       FLIMS(2,NCNT) = 30.
! SHORTWAVE RADIATION r
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "SWRAD"
       IKEYS(NCNT) = 114
       FLIMS(1,NCNT) = 0.
       FLIMS(2,NCNT) = 30.
! BOTTOM FRICTION f
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "BOTTOM STRESS"
       IKEYS(NCNT) = 102
       FLIMS(1,NCNT) = 0.
       FLIMS(2,NCNT) = 10.
! PROCESSOR ID p
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "PROC ID"
       IKEYS(NCNT) = 112
       FLIMS(1,NCNT) = 1.
       FLIMS(2,NCNT) = FLOAT(NPROCS)
! VERTICAL VELOCITY W
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "VERTICAL VELOCITY"
       IKEYS(NCNT) = 87
       FLIMS(1,NCNT) = MINVAL(WW(1:MElem,1:KBM1))
       FLIMS(2,NCNT) = MAXVAL(WW(1:MElem,1:KBM1))
! WET/DRY               
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "WET/DRY"
       IKEYS(NCNT) = 119 
       FLIMS(1,NCNT) = 0
       FLIMS(2,NCNT) = 1
! OMEGA VELOCITY
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "OMEGA"
       IKEYS(NCNT) = 79 
       FLIMS(1,NCNT) = 0
       FLIMS(2,NCNT) = 1
! BOTTOM DRAG COEFFICIENT 
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "CBC"
       IKEYS(NCNT) = 99 
       FLIMS(1,NCNT) = 0
       FLIMS(2,NCNT) = .005 
!      DEPTH (D) -d
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "DEPTH"
       IKEYS(NCNT) = 100
       FLIMS(1,NCNT) = MINVAL(D(1:NNode))
       FLIMS(2,NCNT) = MAXVAL(D(1:NNode))

!      SAND  -s
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "sand"
       IKEYS(NCNT) = 115
       FLIMS(1,NCNT) = 0. 
       FLIMS(2,NCNT) = 10. 

!      MUD  -m
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "mud"
       IKEYS(NCNT) = 109 
       FLIMS(1,NCNT) = 0.
       FLIMS(2,NCNT) = 10.

!      SEDIMENT Thickness Change -t
       NCNT = NCNT + 1
       FKEYS(NCNT) = 1
       TKEYS(NCNT) = "SED Thickness Change"
       IKEYS(NCNT) = 116
       FLIMS(1,NCNT) = -1.
       FLIMS(2,NCNT) = 1.
       





!====================VECTORS===================================
! UA/VA FLOW VECTORS A
       NCNT = NCNT + 1
       FKEYS(NCNT) = 2
       TKEYS(NCNT) = "UA/VA VECS"
       IKEYS(NCNT) = 65
       FLIMS(1,NCNT) = VSCALE
       FLIMS(2,NCNT) = 1.
! UWIND FLOW VECTORS B
       NCNT = NCNT + 1
       FKEYS(NCNT) = 2
       TKEYS(NCNT) = "UWIND VECS"
       IKEYS(NCNT) = 66
       FLIMS(1,NCNT) = VSCALE
       FLIMS(2,NCNT) = 1.
! UWIND FLOW VECTORS C
       NCNT = NCNT + 1
       FKEYS(NCNT) = 2
       TKEYS(NCNT) = "UWIND2 VECS"
       IKEYS(NCNT) = 67
       FLIMS(1,NCNT) = VSCALE
       FLIMS(2,NCNT) = 1.
! 3D VELOCITY FLOW VECTS D
       NCNT = NCNT + 1
       FKEYS(NCNT) = 2
       TKEYS(NCNT) = "3D VELOCITY VECTORS"
       IKEYS(NCNT) = 68
       FLIMS(1,NCNT) = VSCALE
       FLIMS(2,NCNT) = 1.
      DATA TCUTS /'MERIDIONAL','ZONAL','BATHYMETRIC'/
! 1 CONSTANTS
!
! IOPT = 0 STEADY GRID AND DATA
! IOPT = 1 STEADY GRID AND UNSTEADY DATA
! IOPT = 2 UNSTEADY GRID AND DATA
!
! MIRR = 0 NO MIRRORING
! MIRR = 1 MIRRORING ABOUT X=0.0 PLANE
! MIRR = 2 MIRRORING ABOUT Y=0.0 PLANE
! MIRR = 3 MIRRORING ABOUT K=0.0 PLANE
! KNBLOCK NUMBER OF STRUCTURED BLOCKS
      ISTAT = 0
      MIRR = 0
      KNBLOCK = 0
!
! INITIALIZE AND REPORT
!
      CALL PV_INIT(TTL,MYID,' ',' ',IOPT,NCUTS, &
       TCUTS,NKEYS,IKEYS,TKEYS,FKEYS,FLIMS,MIRR,0,KNBLOCK,ISTAT)
      IF(ISTAT /=0) WRITE(*,*) 'ISTAT FROM PV_INIT: ',ISTAT
      WRITE(*,*)'=========VISUALIZATION SERVER INITIALIZED======='
      RETURN
      END SUBROUTINE VIZ_INIT
!===============================================================================
!===============================================================================
         SUBROUTINE pVStruc(knode, kequiv, kcel1, kcel2, kcel3, kcel4, &
                          knptet, kptet, knblck, blocks, khedra, &
                          ksurf, knsurf, hint)
!===============================================================================
      USE ALL_VARS
      INTEGER :: KNODE,KEQUIV,KCEL1,KCEL2,KCEL3,KCEL4,KNPTET,KPTET
      INTEGER :: KNBLCK,BLOCKS(3,*),KHEDRA,KSURF,KNSURF,HINT
!===============================================================================
      KNODE = NNode * KBM1
      KEQUIV = 0
      KCEL1 = 0
      KCEL2 = 0
      KCEL3 = MElem * (KBM1-1)
      KCEL4 = 0
      KNPTET = 0
      KPTET = 0
      KNBLCK = 0
      KHEDRA = 0
      KSURF = MElem * 6
      KNSURF = 2
      RETURN
      END SUBROUTINE PVSTRUC
!===============================================================================
!===============================================================================
      SUBROUTINE PVCELL(CEL1,CEL2,CEL3,CEL4,NPTET,PTET)
!-------------------------------------------------------------------------------
      USE ALL_VARS
      INTEGER, INTENT(OUT) :: CEL1(4,1),CEL2(5,1),CEL3(6,*),CEL4(8,1)
      INTEGER, INTENT(OUT) :: NPTET(8,1),PTET(1)
      INTEGER NNCNT,NCC,K
!
! where: cel1 - node numbers for tetrahedral cells
! (filled to kcel31)
! cel2 - node numbers for pyramid cells
! (filled to kcel32)
! cel3 - node numbers for prism cells
! (filled to kcel33)
! cel4 - node numbers for hexahedral cells
! (filled to kcel34)
! nptet - poly-tetrahedra header information:
! fill first 4 entries for each strip.
! the first entry points to the last position
! in ptet for the end of the strip,
! the next entries are filled with the
! first 3 nodes of the strip.
! ptet - the rest of the poly-tetra, 1 node per
! additional cell
!
! Note: if KCELn is zero the corresponding CELn must NOT be
! filled!
!-------------------------------------------------------------------------------
      NCC = 0
      DO K=1,KBM1-1
      DO NNCNT=1,MElem
        NCC = NCC + 1
        CEL3(1,NCC) = NV(NNCNT,1)+(K-1)*NNode
        CEL3(2,NCC) = NV(NNCNT,1)+(K*NNode)
        CEL3(3,NCC) = NV(NNCNT,2)+(K*NNode)
        CEL3(4,NCC) = NV(NNCNT,2)+(K-1)*NNode
        CEL3(5,NCC) = NV(NNCNT,3)+(K*NNode)
        CEL3(6,NCC) = NV(NNCNT,3)+(K-1)*NNode
        END DO
      END DO
      RETURN
      END SUBROUTINE PVCELL
!===============================================================================
!===============================================================================
!
      SUBROUTINE PVSURFACE(NSURF,SCON,SCEL,TSURF)
!
!-------------------------------------------------------------------------------
      USE ALL_VARS
      INTEGER, INTENT(OUT) :: NSURF(3,*),SCON(1),SCEL(4,*)
      CHARACTER(LEN=20), INTENT(OUT) :: TSURF(*)
      INTEGER :: NNCNT,NCC
!
!
! where: nsurf - pointers to the end of each surface group and
! flag to indicate if surface is on
! scon - connection data for multi-client cases
! scel - node numbers for surfaces
! tsurf - surface group title (optional)
!-------------------------------------------------------------------------------
!
!
! FREE SURFACE
!
      NCC = 0
      DO NNCNT=1,MElem
        NCC = NCC + 1
        SCEL(1,NCC) = NV(NNCNT,1)
        SCEL(2,NCC) = NV(NNCNT,2)
        SCEL(3,NCC) = NV(NNCNT,3)
        SCEL(4,NCC) = 0
      END DO
      NSURF(1,1) = NCC
      NSURF(2,1) = 16
      TSURF(1) = 'Free Surface'
!
! BOTTOM
!
      DO NNCNT=1,MElem
        NCC = NCC + 1
        SCEL(1,NCC) = NV(NNCNT,1) + NNode*(KBM1-1)
        SCEL(2,NCC) = NV(NNCNT,2) + NNode*(KBM1-1)
        SCEL(3,NCC) = NV(NNCNT,3) + NNode*(KBM1-1)
        SCEL(4,NCC) = 0
      END DO
      NSURF(1,2) = NCC
      NSURF(2,2) = 0
      TSURF(2) = 'Bottom'
!
      RETURN
      END SUBROUTINE PVSURFACE
!===============================================================================
!===============================================================================
      SUBROUTINE PVGRID(XYZ)
!-------------------------------------------------------------------------------
      USE ALL_VARS
      REAL*4,INTENT(OUT) :: XYZ(3,*)
      INTEGER :: K,KK,KOUT
      REAL*4 :: LSCALE,DPTHSCALE
!-------------------------------------------------------------------------------
      LSCALE = MAX(ABS(VXMAX-VXMIN),ABS(VYMAX-VYMIN))
      DPTHSCALE = .2*LSCALE/MAXVAL(H)
      KK = 0
      DO KOUT = 1,KBM1
      DO K=1,NNode
        KK = KK + 1
        XYZ(1,KK) = VX(K)
        XYZ(2,KK) = VY(K)
        XYZ(3,KK) = Z(K,KOUT)*H(K)*DPTHSCALE
      END DO
      END DO
      RETURN
      END SUBROUTINE PVGRID
!===============================================================================
!===============================================================================
      SUBROUTINE PVSCAL(KEY,VV)
!-------------------------------------------------------------------------------
      USE ALL_VARS
#     if defined (NG_OI_ASSIM)
      USE MOD_ASSIM
#     endif
#     if defined (WET_DRY)
      USE MOD_WD
#     endif
#     if defined (SEDIMENT)
      USE MOD_SED
#     endif
      REAL*4, INTENT(OUT) :: VV(*)
      INTEGER, INTENT(IN) :: KEY
      INTEGER :: NNCNT,KK,KOUT
      REAL*4, ALLOCATABLE :: VART(:,:)
!-------------------------------------------------------------------------------
      KK = 0
      SELECT CASE (KEY)
! EL
      CASE(1)
      DO KOUT =1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = EL(NNCNT)
      END DO
      END DO
! UA
      CASE(2)
      DO KOUT =1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = SUM(UA(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
      END DO
      END DO
! VA
      CASE(3)
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = SUM(VA(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
      END DO
      END DO
! H
      CASE(4)
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = H(NNCNT)
      END DO
      END DO
! SALINITY
      CASE(5)
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = S1(NNCNT,KOUT)
      END DO
      END DO
! POTENTIAL TEMPERATURE
      CASE(6)
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = T1(NNCNT,KOUT)
      END DO
      END DO
! X VELOCITY
      CASE(7)
      ALLOCATE(VART(NNode,KBM1))
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        VART(NNCNT,KOUT) = &
           SUM( U(NBVE(NNCNT,1:NTVE(NNCNT)),KOUT))/FLOAT(NTVE(NNCNT))
      END DO
      END DO
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = VART(NNCNT,KOUT)
      END DO
      END DO
      DEALLOCATE(VART)
! Y VELOCITY
      CASE(8)
      ALLOCATE(VART(NNode,KBM1))
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        VART(NNCNT,KOUT) = &
          SUM( V(NBVE(NNCNT,1:NTVE(NNCNT)),KOUT))/FLOAT(NTVE(NNCNT))
      END DO
      END DO
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = VART(NNCNT,KOUT)
      END DO
      END DO
      DEALLOCATE(VART)
! TURB LENGTH SCALE
      CASE(9)
      ALLOCATE(VART(NNode,KBM1))
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        VART(NNCNT,KOUT) = &
           SUM( L(NBVE(NNCNT,1:NTVE(NNCNT)),KOUT))/FLOAT(NTVE(NNCNT))
      END DO
      END DO
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = VART(NNCNT,KOUT)
      END DO
      END DO
      DEALLOCATE(VART)
! TURB KINETIC ENERGY
      CASE(10)
      ALLOCATE(VART(NNode,KBM1))
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        VART(NNCNT,KOUT) = &
         SUM( Q2(NBVE(NNCNT,1:NTVE(NNCNT)),KOUT))/FLOAT(NTVE(NNCNT))
      END DO
      END DO
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = VART(NNCNT,KOUT)
        END DO
        END DO
        DEALLOCATE(VART)
! WTSURF
        CASE(11)
        write(*,*)'setting wtsurf'
        DO KOUT=1,KBM1
        DO NNCNT=1,NNode
          KK = KK + 1
          VV(KK) = WTSURF(NNCNT)
        END DO
        END DO
! SWRAD
        CASE(12)
        DO KOUT=1,KBM1
        DO NNCNT=1,NNode
          KK = KK + 1
          VV(KK) = SWRAD(NNCNT)
        END DO
        END DO
! BOTTOM FRICTION
        CASE(13)
        DO KOUT =1,KBM1
        DO NNCNT=1,NNode
          KK = KK + 1
          VV(KK) = &
            SUM(WUBOT(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
        END DO
        END DO
! PROC ID
        CASE(14)
        DO KOUT=1,KBM1
        DO NNCNT=1,NNode
          KK = KK + 1
          VV(KK) = FLOAT(MYID)
!          VV(KK) = DA_SST(NNCNT) 
!          VV(KK) = &
!            SUM(DA_CUR(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
        END DO
        END DO
! VERTICAL VELOCITY
        CASE(15)
        ALLOCATE(VART(NNode,KBM1))
        DO KOUT=1,KBM1
        DO NNCNT=1,NNode
          VART(NNCNT,KOUT) = &
           SUM( WW(NBVE(NNCNT,1:NTVE(NNCNT)),KOUT))/FLOAT(NTVE(NNCNT))
        END DO
        END DO
        DO KOUT=1,KBM1
        DO NNCNT=1,NNode
          KK = KK + 1
          VV(KK) = VART(NNCNT,KOUT)
        END DO
        END DO
        DEALLOCATE(VART)
! WET/DRY FLAG       
        CASE(16)
        DO KOUT =1,KBM1
        DO NNCNT=1,NNode
          KK = KK + 1
#         if defined (WET_DRY)
          VV(KK) = ISWET_NODE_CURRENTSTEP(NNCNT)
#         else
          VV(KK) = 1.0
#         endif
        END DO
        END DO

!     OMEGA 
      CASE(17)
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = WTS(NNCNT,KOUT)
      END DO
      END DO

!     CBC 
      CASE(18)
      DO KOUT =1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = SUM(CBC(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
      END DO
      END DO
!     DEPTH
      CASE(19)
      DO KOUT =1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = D(NNCNT)
      END DO
      END DO
!     SAND SEDIMENT 
      CASE(20)
#     if defined(SEDIMENT)
      IF(SEDIMENT_ON)THEN
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = SED(1)%CONC(NNCNT,KOUT)
      END DO
      END DO
      END IF
!     MUD SEDIMENT 
      CASE(21)
      IF(SEDIMENT_ON .AND. NSED > 1)THEN
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = SED(2)%CONC(NNCNT,KOUT)
      END DO
      END DO
      END IF

!     SED THICKNESS CHANGE
      CASE(22)
      DO KOUT=1,KBM1
      DO NNCNT=1,NNode
        KK = KK + 1
        VV(KK) = BOTTOM(NNCNT,DTHCK)
      END DO
      END DO
#     endif





        END SELECT

        RETURN
        END SUBROUTINE PVSCAL
!===============================================================================
!===============================================================================
        SUBROUTINE PVVECT(KEY,VV)
        USE ALL_VARS
        REAL*4, INTENT(OUT) :: VV(3,*)
        INTEGER, INTENT(IN) :: KEY
        INTEGER :: NNCNT,KK,KOUT
        SELECT CASE (KEY)
        CASE(23)
        KK = 0
        DO KOUT=1,KBM1
        DO NNCNT=1,NNode
          KK = KK +1
          VV(1,KK) = &
            SUM(UA(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
         VV(2,KK) = &
             SUM(VA(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
         VV(3,KK) = 0.
       END DO
       END DO
       CASE(24)
       KK = 0
       DO KOUT=1,KBM1
       DO NNCNT=1,NNode
         KK = KK +1
         VV(1,KK) = &
            SUM(WUSURF(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
         VV(2,KK) = &
            SUM(WVSURF(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
         VV(3,KK) = 0.
       END DO
       END DO
       CASE(25)
       KK = 0
       DO KOUT=1,KBM1
       DO NNCNT=1,NNode
         KK = KK +1
         VV(1,KK) = &
            SUM(WUSURF2(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
         VV(2,KK) = &
            SUM(WVSURF2(NBVE(NNCNT,1:NTVE(NNCNT))))/FLOAT(NTVE(NNCNT))
         VV(3,KK) = 0.
       END DO
       END DO
       CASE(26)
       KK = 0
       DO KOUT=1,KBM1
       DO NNCNT=1,NNode
         KK = KK +1
         VV(1,KK) = &
           SUM(U(NBVE(NNCNT,1:NTVE(NNCNT)),KOUT))/FLOAT(NTVE(NNCNT))
         VV(2,KK) = &
           SUM(V(NBVE(NNCNT,1:NTVE(NNCNT)),KOUT))/FLOAT(NTVE(NNCNT))
         VV(3,KK) = &
           SUM(W(NBVE(NNCNT,1:NTVE(NNCNT)),KOUT))/FLOAT(NTVE(NNCNT))
       END DO
       END DO
       END SELECT
       RETURN
       END
!===============================================================================
!===============================================================================
!===============================================================================
      SUBROUTINE PVZPRIME(IDCUT,XYZ,NNODE,ZP,ZPRIME,XPC,YPC,HALFW)
!-------------------------------------------------------------------------------
      USE ALL_VARS
      INTEGER, INTENT(IN) :: IDCUT,NNODE
      REAL*4,INTENT(IN) :: XYZ(3,*)
      REAL*4,INTENT(OUT) :: ZP(*)
      REAL*4, INTENT(OUT) :: ZPRIME,XPC,YPC,HALFW
      REAL*4 XMAX,XMIN,YMAX,YMIN,ZMIN,ZMAX,LSCALE
      INTEGER I
!-------------------------------------------------------------------------------
      XMAX = MAXVAL(VX(1:NNode))
      YMAX = MAXVAL(VY(1:NNode))
      ZMAX = MAXVAL(H(1:NNode))
      XMIN = MINVAL(VX(1:NNode))
      YMIN = MINVAL(VY(1:NNode))
      ZMIN = MINVAL(H(1:NNode))
      XMID = .5*(XMAX+XMIN)
      YMID = .5*(YMAX+YMIN)
      ZMID = .5*(ZMAX+ZMIN)
      LXYSCALE = MAX(ABS(XMAX-XMIN),ABS(YMAX-YMIN))
      LXZSCALE = MAX(ABS(XMAX-XMIN),ABS(ZMAX-ZMIN))
      LYZSCALE = MAX(ABS(ZMAX-ZMIN),ABS(YMAX-YMIN))
      SELECT CASE(IDCUT)
! CONSTANT X CUT
      CASE(1)
      DO I=1,NNODE
        ZP(I) = XYZ(1,I)
      END DO
      ZPRIME = XMID
      XPC = YMID
      YPC = ZMID
      HALFW = .6*LYZSCALE
! CONSTANT Y CUT
      CASE(2)
      DO I=1,NNODE
        ZP(I) = XYZ(2,I)
      END DO
      ZPRIME = YMID
      XPC = XMID
      YPC = ZMID
      HALFW = .6*LXZSCALE
! CONSTANT Z CUT
      CASE(3)
      DO I=1,NNODE
        ZP(I) = XYZ(3,I)
      END DO
      ZPRIME = -.0001
      XPC = XMID
      YPC = YMID
      HALFW = .6*LXYSCALE
      END SELECT
      RETURN
      END SUBROUTINE PVZPRIME
!===============================================================================
!===============================================================================
      SUBROUTINE PVXYPRIME(ZPRIME,KN,XYZ,MElem,XYP)
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KN(*)
      REAL*4,INTENT(IN) :: XYZ(3,*)
      INTEGER, INTENT(IN) :: MElem
      REAL*4,INTENT(IN) :: ZPRIME
      REAL*4, INTENT(OUT) :: XYP(2,*)
      INTEGER I
!-------------------------------------------------------------------------------
      DO I=1,MElem
        XYP(1,I) = XYZ(1,KN(I))
        XYP(2,I) = XYZ(3,KN(I))
      END DO
      RETURN
      END SUBROUTINE PVXYPRIME
#     else
      SUBROUTINE DUM4
      END SUBROUTINE DUM4
#     endif
