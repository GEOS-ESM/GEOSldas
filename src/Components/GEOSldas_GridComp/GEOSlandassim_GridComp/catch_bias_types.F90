
module catch_bias_types

  ! definition of bias types for Catchment land surface model
  !
  ! reichle,        17 Oct 2005
  ! reichle+draper, 26 Mar 2013 - prep for separation of model and obs bias 
  ! reichle+draper, 28 Aug 2013 - added "obs_bias" type
  !
  ! -------------------------------------------------------------------
  
  use catch_constants,                  ONLY:     &
       N_snow => CATCH_N_SNOW,                    &
       N_gt   => CATCH_N_GT
  
  use catch_types,                      ONLY:     &
       cat_progn_type
  
  implicit none

  save
  
  ! everything is private by default unless made public
  
  private
  
  public :: cat_progn_int_type
  public :: cat_bias_param_type

  public :: obs_bias_type

  ! ----------------------------------------------------------------------
  !
  ! *INTEGER* Catchment model prognostic variables
  !
  ! THESE MUST MATCH THE "real" cat_progn TYPE IN catch_types.f90

  type :: cat_progn_int_type

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

     integer :: tc1     ! surface/canopy temperature 
     integer :: tc2
     integer :: tc4
            
     integer :: qa1     ! specific humidity in canopy air
     integer :: qa2
     integer :: qa4
            
     integer :: capac   ! canopy interception water
            
     integer :: catdef  ! catchment deficit
     integer :: rzexc   ! root zone excess
     integer :: srfexc  ! surface excess
            
     integer, dimension(N_gt)   :: ght     ! ground heat content
            
     integer, dimension(N_snow) :: wesn    ! snow water equivalent
     integer, dimension(N_snow) :: htsn    ! snow heat content
     integer, dimension(N_snow) :: sndz    ! snow depth
  
  end type cat_progn_int_type

  ! -----------------------------------------
  
  type :: cat_bias_param_type

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

     type(cat_progn_type)     :: tconst
     type(cat_progn_type)     :: trelax
     
     type(cat_progn_int_type) :: Nparam
     
  end type cat_bias_param_type
  
  ! -----------------------------------------
  
  type :: obs_bias_type
     
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     
     real                     :: bias	 ! observation bias 
     integer, dimension(2)    :: tcount  ! count for time since last obs  [seconds]
                                         ! at start of assim cycle, for each tile: 
					 ! tcount(1) = time since most recent obs.
					 ! tcount(1) = time since 2nd most recent obs.
          
  end type obs_bias_type
  
  ! -------------------------------------------------------------------
    
end module catch_bias_types

! ======================== EOF ===========================================
