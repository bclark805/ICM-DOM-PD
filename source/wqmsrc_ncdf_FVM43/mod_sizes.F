
        MODULE MOD_SIZES
        !***** Parameter declarations
        !***** ST. JOHNS RIVER COARSE GRID:6 LAYER MODEL
           INTEGER,PARAMETER :: NCP=44,                               & ! changed to include new CDOM formulations
                                NQFP=32686,                           &
                                NHQP=22086,                           &
                                NS1P=378,                             &
                                NS2P=378,                             &
                                NS3P=378,                             &
                                NBCP=492,                             &
                                NDP=500,                              &
                                NFLP=20,                              &
                                NOIP=10,                              &
                                NSSFP=10 
           INTEGER :: MGL,NGL,OBCGL,NOBTY  
        END MODULE MOD_SIZES
