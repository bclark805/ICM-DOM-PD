MODULE MOD_BUFFERS

#if defined MULTIPROCESSOR
	USE MOD_PREC, ONLY : CDF_PREC,SP,DP, MPI_DP,MPI_F
#else
	USE MOD_PREC, ONLY : CDF_PREC,SP,DP
#endif
	

  REAL(SP), ALLOCATABLE, DIMENSION(:) :: RCVBUF,SNDBUF          !floating point buffers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: RCVBUF_INT, SNDBUF_INT  !integer buffers

!  CONTAINS

  !intialize buffer


  !delete buffer


  !expand buffer


  !shrink buffer 
 

END MODULE MOD_BUFFERS

