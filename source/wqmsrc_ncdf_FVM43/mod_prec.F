!===============================================================================!
! DEFINE FLOATING POINT PRECISION USING KIND                                    !
!===============================================================================!
MODULE MOD_PREC
   IMPLICIT NONE
#  if defined (MULTIPROCESSOR)
   include "mpif.h"
#  endif

!--Single Precision Coding------------------------------------------------------!
#  if !defined (DOUBLE_PRECISION)
   INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,30)
#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_F = MPI_REAL
#  endif
   
#  else 
!--Double Precision Coding------------------------------------------------------!
   INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(12,300)
#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_F = MPI_DOUBLE_PRECISION
#  endif
#  endif


  
   INTEGER, PARAMETER :: DP     = SELECTED_REAL_KIND(12,300)
#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_DP = MPI_DOUBLE_PRECISION
#  endif
   
! KURT GLAESEMANN SEPT 22 2009 - USE CORRECT PRECISION WHEN READ NETCDF
  INTEGER, PARAMETER :: CDF_PREC = 4

! KURT GLAESEMANN 14 April 2015 - USE CORRECT PRECISION WHEN USING MPI ON CDF DATA
#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_CDF = MPI_REAL
#  endif

END MODULE MOD_PREC

