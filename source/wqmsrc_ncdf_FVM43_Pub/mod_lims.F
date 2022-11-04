        !==============================================================================|
        !   GLOBAL LIMITS AND ARRAY SIZING PARAMETERS                                  !
        !==============================================================================|

        MODULE MOD_LIMS

          IMPLICIT NONE
          SAVE

          INTEGER NLOC               !!NUMBER OF ELEMENTS 
          INTEGER MLOC               !!NUMBER OF NODES
          INTEGER NISBCE_1           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 1
          INTEGER NISBCE_2           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 2
          INTEGER NISBCE_3           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 3

          INTEGER KB                 !!NUMBER OF SIGMA LEVELS
          INTEGER KBM1               !!NUMBER OF SIGMA LEVELS-1
          INTEGER KBM2               !!NUMBER OF SIGMA LEVELS-2
          INTEGER MYID               !!UNIQUE PROCESSOR ID (1 => NPROCS)
          INTEGER NPROCS             !!NUMBER OF PROCESSORS
          INTEGER NE                 !!NUMBER OF UNIQUE EDGES
          INTEGER NCV                !!NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)
           
          INTEGER IINT               !!TYKIM added for nudging
          INTEGER NCV_I              !!NUMBER OF INTERNAL CONTROL VOLUMES (LOCAL ONLY)
          INTEGER NTLOC              !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
          INTEGER MTLOC              !!TOTAL OF LOCAL INTERNAL + HALO NODES
          INTEGER NCT				 !!(NTLOC) *3
          INTEGER MX_NBR_ELEM        !!MAX NUMBER OF ELEMENTS SURROUNDING A NODE

          INTEGER NUMQBC_GL,NUMPNT_GL
          INTEGER NUMQBC,NUMPNT
          INTEGER NstationMax
          PARAMETER (NstationMax=200)

                                     !Maximum number of stations 
                                     !note this is predifined here because this is going into
                                     !defintion of NstationNum_GL for NAMELIST /wqm_stations/ 
                                     !fortran90 does not support dynamic arrays in namelist yet
                                     
        END MODULE MOD_LIMS

