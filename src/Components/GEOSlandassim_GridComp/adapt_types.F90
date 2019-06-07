
module adapt_types

  ! definition of types for adaptive filtering
  !
  ! reichle, 19 Jul 2007
  !
  ! -------------------------------------------------------------------
  
  implicit none

  save
  
  ! everything is private by default unless made public
  
  private
  
  public :: adapt_misc_param_type
  
  ! ----------------------------------------------------------------------
    
  type :: adapt_misc_param_type
     
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

     real    :: gamma_P
     real    :: gamma_R
     
     real    :: delta_P
     real    :: delta_R
     
     real    :: beta_P
     real    :: beta_R
     
     real    :: min_alpha_P
     real    :: max_alpha_P

     real    :: min_alpha_R
     real    :: max_alpha_R
          
  end type adapt_misc_param_type

  ! -------------------------------------------------------------------
  
end module adapt_types

! ======================== EOF ===========================================
