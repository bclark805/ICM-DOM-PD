!=======================================================================
! FVCOM Sediment Initialization Routine (User Defined)
!   Configuration:  default    
!=======================================================================
  Subroutine Init_Sed  
# if defined (SEDIMENT)
  Use Mod_Prec 
  Use Mod_Sed
  Use Lims, only: NNode,kbm1
  implicit none 
  integer :: i,k,ised
  real(sp) :: bed_thickness

  !--------------------------------------------------
  !Initialize Bed Properties
  !--------------------------------------------------
  bed_thickness = 10.00 !meters 

  Do k=1,Nbed
    Do i=1,NNode
       bed(i,k,iaged) = 0.00
       bed(i,k,ithck) = bed_thickness/float(Nbed)
       bed(i,k,iporo) = 0.40
     End Do
   End do

  !--------------------------------------------------
  !Initialize Bed_Frac properties
  !--------------------------------------------------

  Do k=1,Nbed
    Do i=1,NNode
       Do ised=1,Nsed
         sed(ised)%frac(i,k)=1.0/float(nsed)   
       End Do
     End Do
  End Do

  !--------------------------------------------------
  !Initialize Sediment Concentrations
  !--------------------------------------------------
  Do k=1,Kbm1
    Do i=1,NNode
      Do ised=1,Nsed
        sed(ised)%conc(i,k) = 1.0
        sed(ised)%cnew(i,k) = sed(ised)%conc(i,k) 
      End Do
    End do
  End do

  Return
# endif
  End Subroutine Init_Sed
!=======================================================================
