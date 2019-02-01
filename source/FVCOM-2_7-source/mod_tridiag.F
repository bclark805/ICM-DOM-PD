!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mtridiagonal --- solver for tri-diagonal matrices \label{sec:tridiagonal}
!
! !INTERFACE:
   MODULE mtridiagonal_scal
   use mod_prec
!
! !DESCRIPTION: 
!
!  Solves a linear system of equations with a tridiagonal matrix
!  using Gaussian elimination. 
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_tridiagonal_scal,tridiagonal_scal
!
! !PUBLIC DATA MEMBERS:
   real(sp), dimension(:), allocatable     :: au_m,bu_m,cu_m,du_m
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  $Log: mod_tridiag.F,v $
!  Revision 1.1  2006/06/20 00:14:52  gcowles
!  *** empty log message ***
!
!  Revision 1.4  2003/03/28 09:20:36  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:06:33  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:54:16  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
!  private data members
   real(sp), private, dimension(:),allocatable  ::  ru_m,qu_m
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory
!
! !INTERFACE:
   subroutine init_tridiagonal_scal(N)
!
! !DESCRIPTION:
!  This routines allocates memory necessary to perform the Gaussian 
! elimination.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!
!
!-----------------------------------------------------------------------
!BOC
   if(allocated(au_m))then
     deallocate(au_m)
     deallocate(bu_m)
     deallocate(cu_m)
     deallocate(du_m)
     deallocate(ru_m)
     deallocate(qu_m)
   endif

   allocate(au_m(0:N)) ; au_m = 0.
   allocate(bu_m(0:N)) ; bu_m = 0.
   allocate(cu_m(0:N)) ; cu_m = 0.
   allocate(du_m(0:N)) ; du_m = 0.
   allocate(ru_m(0:N)) ; ru_m = 0.
   allocate(qu_m(0:N)) ; qu_m = 0.

   return
   end subroutine init_tridiagonal_scal
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Simplified Gaussian elimination 
!
! !INTERFACE:
   subroutine tridiagonal_scal(N,fi,lt,value)
!
! !DESCRIPTION:
! A linear equation with tridiagonal matrix is solved here. The main
! diagonal is stored on {\tt bu}, the upper diagonal on {\tt au}, and the
! lower diagonal on {\tt cu}, the right hand side is stored on {\tt du}. 
! The method used here is the simplified Gauss elimination, also called 
! \emph{Thomas algorithm}.  
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N,fi,lt
!
! !OUTPUT PARAMETERS:
   real(sp)                                    :: value(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  $Log: mod_tridiag.F,v $
!  Revision 1.1  2006/06/20 00:14:52  gcowles
!  *** empty log message ***
!
!  Revision 1.4  2003/03/28 09:20:36  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:06:33  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:54:16  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
   ru_m(lt)=au_m(lt)/bu_m(lt)
   qu_m(lt)=du_m(lt)/bu_m(lt)

   do i=lt-1,fi+1,-1
      ru_m(i)=au_m(i)/(bu_m(i)-cu_m(i)*ru_m(i+1))
      qu_m(i)=(du_m(i)-cu_m(i)*qu_m(i+1))/(bu_m(i)-cu_m(i)*ru_m(i+1))
   end do

   qu_m(fi)=(du_m(fi)-cu_m(fi)*qu_m(fi+1))/(bu_m(fi)-cu_m(fi)*ru_m(fi+1))

   value(fi)=qu_m(fi)
   do i=fi+1,lt
      value(i)=qu_m(i)-ru_m(i)*value(i-1)
   end do


   return
   end subroutine tridiagonal_scal
!EOC

!-----------------------------------------------------------------------

   end module mtridiagonal_scal

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
