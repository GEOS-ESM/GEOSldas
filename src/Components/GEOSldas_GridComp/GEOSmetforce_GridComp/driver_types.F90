
module driver_types
  
  ! definition of types and associated operators for Catchment Model driver
  !
  ! IMPORTANT:
  ! When adding a field to any of the derived types, must also update
  ! the associated assignment and operator definitions.
  ! THERE IS NO WARNING/ERROR IF OPERATOR IS NOT DEFINED FOR ALL FIELDS!
  !
  ! reichle, 10 May 2005
  ! reichle, 10 Jun 2005 - converted met_force_type to ALMA
  ! reichle, 10 Sep 2007 - added modis_alb_param_type
  ! reichle+qliu,  8 Oct 2008 (and earlier) - added fields "RefH" and "SWnet"
  !                        for DAS/LDASsa integration
  ! reichle, 23 Feb 2009 - added fields ParDrct, ParDffs for MERRA
  ! reichle,  5 Mar 2009 - deleted ParDrct, ParDffs after testing found no impact
  ! reichle, 30 Sep 2009 - changed "out_avg_type" to "out_select_type"
  ! reichle,  8 Dec 2011 - added "veg_param_type" and "bal_diagn_type"
  ! reichle, 20 Dec 2011 - reinstated met_force fields "PARdrct" and "PARdffs"
  !                        for MERRA-Land output specs
  ! reichle,  5 Apr 2013 - removed modis_alb_param_type fields "sc_albvr", "sc_albnr" 
  ! reichle, 23 Jul 2013 - renamed "modis_alb_param" --> "alb_param"
  !
  ! --------------------------------------------------------------------------
  
  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: met_force_type, veg_param_type, bal_diagn_type
  public :: out_dtstep_type
  public :: out_select_type, out_select_sub_type
  public :: out_choice_type, out_choice_time_type
  public :: alb_param_type
  public :: assignment (=), operator (/), operator (+), operator (*)
  
  ! ---------------------------------------------------------------------
  !
  ! meteorological forcing variables
  !
  ! The Catchment model driver requires forcing fields of the types
  ! and units in the structure met_force_type.
  ! Make sure to convert the native types and units of each forcing
  ! data set into the types and units of met_force_type right after
  ! the native data have been read.  See for example get_Berg_netcdf().
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: met_force_type
     real :: Tair                 ! air temperature at RefH                 [K]
     real :: Qair                 ! specific humidity at RefH               [kg/kg]
     real :: Psurf                ! surface pressure                        [Pa]
     real :: Rainf_C              ! convective rainfall                     [kg/m2/s]
     real :: Rainf                ! total rainfall                          [kg/m2/s]
     real :: Snowf                ! total snowfall                          [kg/m2/s]
     real :: LWdown               ! downward longwave radiation             [W/m2]
     real :: SWdown               ! downward shortwave radiation            [W/m2]
     real :: SWnet                ! downward net shortwave radiation        [W/m2]
     real :: PARdrct              ! Photosynth. Active Radiation (direct)   [W/m2]
     real :: PARdffs              ! Photosynth. Active Radiation (diffuse)  [W/m2]
     real :: Wind                 ! wind speed at RefH                      [m/s]
     real :: RefH                 ! reference height for Tair, Qair, Wind   [m]
  end type met_force_type

  ! ---------------------------------------------------------------------
  !
  ! vegetation variables
  !
  ! The Catchment model requires seasonally varying greenness and leaf-area-index
    
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: veg_param_type
     real :: grn                  ! vegetation greenness fraction           [-]
     real :: lai                  ! leaf-area-index                         [m2/m2]
  end type veg_param_type

  ! ---------------------------------------------------------------------
  !
  ! water and energy balance diagnostic variables

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: bal_diagn_type
     real :: etotl                ! total energy store (in all model progn)  [J/m2]
     real :: echng                ! energy change per unit time (model only) [W/m2]
     real :: eincr                ! energy analysis increment per unit time  [W/m2]
     real :: wtotl                ! total water store (in all model progn)   [kg/m2]
     real :: wchng                ! water change per unit time (model only)  [kg/m2/s]
     real :: wincr                ! water analysis increment per unit time   [kg/m2/s]
  end type bal_diagn_type
  
  ! ---------------------------------------------------------------
  !
  ! type output time steps

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: out_dtstep_type
     integer :: rstrt
     integer :: inst
     integer :: xhourly
  end type out_dtstep_type
  
  ! ---------------------------------------------------------------
  !
  ! type for reading in output choices from namelist file

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type :: out_select_sub_type
     logical :: inst
     logical :: xhourly
     logical :: daily
     logical :: pentad
     logical :: monthly
  end type out_select_sub_type
  
  type :: out_select_type
     type(out_select_sub_type) :: tile
     type(out_select_sub_type) :: grid
  end type out_select_type

  ! ---------------------------------------------------------------
  !
  ! type for output choices *after* processing in read_driver_inputs()

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  type :: out_choice_space_type
     logical :: tile
     logical :: grid
     logical :: any
  end type out_choice_space_type
  
  type :: out_choice_type
     type(out_choice_space_type) :: inst
     type(out_choice_space_type) :: xhourly
     type(out_choice_space_type) :: daily
     type(out_choice_space_type) :: pentad
     type(out_choice_space_type) :: monthly
     type(out_choice_space_type) :: any
  end type out_choice_type
      
  type :: out_choice_time_type
     logical :: rstrt
     logical :: inst
     logical :: xhourly
     logical :: daily
     logical :: pentad
     logical :: monthly
     logical :: any_non_rstrt
  end type out_choice_time_type
  
  ! ---------------------------------------------------------------
  !
  ! albedo scaling factors

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: alb_param_type     
     real :: sc_albvf    ! Scaling factor for diffuse visible  or whitesky 0.3-0.7 
     real :: sc_albnf    ! Scaling factor for diffuse infrared or whitesky 0.7-5.0
  end type alb_param_type

  ! --------------------------------------------------------------
  
  interface assignment (=)
     module procedure scalar2met_force
     module procedure scalar2veg_param
     module procedure scalar2bal_diagn
     module procedure scalar2alb_param
  end interface
  
  interface operator (/)
     module procedure met_force_div_scalar
     module procedure veg_param_div_scalar
     module procedure bal_diagn_div_scalar
  end interface
  
  interface operator (*)
     module procedure alb_param_times_scalar ! Need both definitions to 
     module procedure scalar_times_alb_param ! define a*b and b*a
  end interface
  
  interface operator (+)
     module procedure add_met_force
     module procedure add_veg_param
     module procedure add_bal_diagn
     module procedure add_alb_param
  end interface
  
contains

  ! --------------------------------------------------

  subroutine scalar2met_force( met_force, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(met_force_type), intent(out) :: met_force
    
    met_force%Tair     = scalar
    met_force%Qair     = scalar
    met_force%Psurf    = scalar
    met_force%Rainf_C  = scalar
    met_force%Rainf    = scalar
    met_force%Snowf    = scalar
    met_force%LWdown   = scalar
    met_force%SWdown   = scalar
    met_force%SWnet    = scalar
    met_force%PARdrct  = scalar
    met_force%PARdffs  = scalar
    met_force%Wind     = scalar
    met_force%RefH     = scalar

  end subroutine scalar2met_force
  
  ! ---------------------------------------------------
  
  function met_force_div_scalar( met_force, scalar )
    
    implicit none
    
    type(met_force_type)             :: met_force_div_scalar
    type(met_force_type), intent(in) :: met_force
    
    real, intent(in) :: scalar
    
    met_force_div_scalar%Tair     =     met_force%Tair     / scalar
    met_force_div_scalar%Qair     =     met_force%Qair     / scalar
    met_force_div_scalar%Psurf    =     met_force%Psurf    / scalar
    met_force_div_scalar%Rainf_C  =     met_force%Rainf_C  / scalar
    met_force_div_scalar%Rainf    =     met_force%Rainf    / scalar
    met_force_div_scalar%Snowf    =     met_force%Snowf    / scalar
    met_force_div_scalar%LWdown   =     met_force%LWdown   / scalar
    met_force_div_scalar%SWdown   =     met_force%SWdown   / scalar
    met_force_div_scalar%SWnet    =     met_force%SWnet    / scalar
    met_force_div_scalar%PARdrct  =     met_force%PARdrct  / scalar
    met_force_div_scalar%PARdffs  =     met_force%PARdffs  / scalar
    met_force_div_scalar%Wind     =     met_force%Wind     / scalar
    met_force_div_scalar%RefH     =     met_force%RefH     / scalar
        
  end function met_force_div_scalar

  ! -----------------------------------------------------------

  function add_met_force( met_force_1, met_force_2 )
    
    implicit none

    type(met_force_type)             :: add_met_force
    type(met_force_type), intent(in) :: met_force_1, met_force_2
    
    add_met_force%Tair     = met_force_1%Tair     + met_force_2%Tair    
    add_met_force%Qair     = met_force_1%Qair     + met_force_2%Qair     
    add_met_force%Psurf    = met_force_1%Psurf    + met_force_2%Psurf     
    add_met_force%Rainf_C  = met_force_1%Rainf_C  + met_force_2%Rainf_C  
    add_met_force%Rainf    = met_force_1%Rainf    + met_force_2%Rainf   
    add_met_force%Snowf    = met_force_1%Snowf    + met_force_2%Snowf    
    add_met_force%LWdown   = met_force_1%LWdown   + met_force_2%LWdown   
    add_met_force%SWdown   = met_force_1%SWdown   + met_force_2%SWdown   
    add_met_force%SWnet    = met_force_1%SWnet    + met_force_2%SWnet   
    add_met_force%PARdrct  = met_force_1%PARdrct  + met_force_2%PARdrct   
    add_met_force%PARdffs  = met_force_1%PARdffs  + met_force_2%PARdffs   
    add_met_force%Wind     = met_force_1%Wind     + met_force_2%Wind    
    add_met_force%RefH     = met_force_1%RefH     + met_force_2%RefH    
    
  end function add_met_force

  ! ---------------------------------------------------

  subroutine scalar2veg_param( veg_param, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(veg_param_type), intent(out) :: veg_param
    
    veg_param%grn      = scalar
    veg_param%lai      = scalar

  end subroutine scalar2veg_param
  
  ! -----------------------------------------------------------
  
  function veg_param_div_scalar( veg_param, scalar )
    
    implicit none
    
    type(veg_param_type)             :: veg_param_div_scalar
    type(veg_param_type), intent(in) :: veg_param
    
    real, intent(in) :: scalar
    
    veg_param_div_scalar%grn      =     veg_param%grn      / scalar
    veg_param_div_scalar%lai      =     veg_param%lai      / scalar
    
  end function veg_param_div_scalar
  
  ! -----------------------------------------------------------

  function add_veg_param( veg_param_1, veg_param_2 )
    
    implicit none

    type(veg_param_type)             :: add_veg_param
    type(veg_param_type), intent(in) :: veg_param_1, veg_param_2
    
    add_veg_param%grn      = veg_param_1%grn      + veg_param_2%grn    
    add_veg_param%lai      = veg_param_1%lai      + veg_param_2%lai     
    
  end function add_veg_param

  ! ---------------------------------------------------

  subroutine scalar2bal_diagn( bal_diagn, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(bal_diagn_type), intent(out) :: bal_diagn
    
    bal_diagn%etotl    = scalar
    bal_diagn%echng    = scalar
    bal_diagn%eincr    = scalar
    bal_diagn%wtotl    = scalar
    bal_diagn%wchng    = scalar
    bal_diagn%wincr    = scalar

  end subroutine scalar2bal_diagn
  
  ! -----------------------------------------------------------
  
  function bal_diagn_div_scalar( bal_diagn, scalar )
    
    implicit none
    
    type(bal_diagn_type)             :: bal_diagn_div_scalar
    type(bal_diagn_type), intent(in) :: bal_diagn
    
    real, intent(in) :: scalar
    
    bal_diagn_div_scalar%etotl    =     bal_diagn%etotl    / scalar
    bal_diagn_div_scalar%echng    =     bal_diagn%echng    / scalar
    bal_diagn_div_scalar%eincr    =     bal_diagn%eincr    / scalar
    bal_diagn_div_scalar%wtotl    =     bal_diagn%wtotl    / scalar
    bal_diagn_div_scalar%wchng    =     bal_diagn%wchng    / scalar
    bal_diagn_div_scalar%wincr    =     bal_diagn%wincr    / scalar
        
  end function bal_diagn_div_scalar

  ! -----------------------------------------------------------

  function add_bal_diagn( bal_diagn_1, bal_diagn_2 )
    
    implicit none

    type(bal_diagn_type)             :: add_bal_diagn
    type(bal_diagn_type), intent(in) :: bal_diagn_1, bal_diagn_2
    
    add_bal_diagn%etotl    = bal_diagn_1%etotl    + bal_diagn_2%etotl  
    add_bal_diagn%echng    = bal_diagn_1%echng    + bal_diagn_2%echng  
    add_bal_diagn%eincr    = bal_diagn_1%eincr    + bal_diagn_2%eincr   
    add_bal_diagn%wtotl    = bal_diagn_1%wtotl    + bal_diagn_2%wtotl  
    add_bal_diagn%wchng    = bal_diagn_1%wchng    + bal_diagn_2%wchng  
    add_bal_diagn%wincr    = bal_diagn_1%wincr    + bal_diagn_2%wincr  
    
  end function add_bal_diagn

  ! -----------------------------------------------------------

  subroutine scalar2alb_param( alb_param, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(alb_param_type), intent(out) :: alb_param
    
    alb_param%sc_albvf = scalar
    alb_param%sc_albnf = scalar
    
  end subroutine scalar2alb_param

  ! ---------------------------------------------------
  
  function alb_param_times_scalar( alb_param, scalar )
    
    implicit none
    
    type(alb_param_type)             :: alb_param_times_scalar
    type(alb_param_type), intent(in) :: alb_param
    
    real, intent(in) :: scalar
    
    alb_param_times_scalar%sc_albvf = alb_param%sc_albvf * scalar
    alb_param_times_scalar%sc_albnf = alb_param%sc_albnf * scalar
    
  end function alb_param_times_scalar

  ! ---------------------------------------------------
  
  function scalar_times_alb_param( scalar, alb_param )
    
    implicit none
    
    type(alb_param_type)             :: scalar_times_alb_param
    type(alb_param_type), intent(in) :: alb_param
    
    real, intent(in) :: scalar
    
    scalar_times_alb_param%sc_albvf = alb_param%sc_albvf * scalar
    scalar_times_alb_param%sc_albnf = alb_param%sc_albnf * scalar
    
  end function scalar_times_alb_param

  ! -----------------------------------------------------------

  function add_alb_param( alb_param_1, alb_param_2 )
    
    implicit none

    type(alb_param_type)             :: add_alb_param
    type(alb_param_type), intent(in) :: alb_param_1, alb_param_2
    
    add_alb_param%sc_albvf = alb_param_1%sc_albvf + alb_param_2%sc_albvf
    add_alb_param%sc_albnf = alb_param_1%sc_albnf + alb_param_2%sc_albnf 
    
  end function add_alb_param

  ! -----------------------------------------------------------
  
end module driver_types


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if 0

program test
  
  use driver_types

  type(alb_param_type), dimension(2,3) :: map1, map2
  type(alb_param_type), dimension(2)   :: map3

  integer :: n, m

  n=1
  m=3

  map1(n,m) =  1.111
  map2(n,m) =  2.222
  map3(n) =  0.

  write (*,*) map1(n,m)
  write (*,*) map2(n,m)
  write (*,*) map3(n)
  write (*,*)

  map3(n) = map1(n,m) + map2(n,m)

  write (*,*) map1(n,m)
  write (*,*) map2(n,m)
  write (*,*) map
  write (*,*)

  map3(n)%sc_albnf = map2(n,m)%sc_albnf * 40.
  
  write (*,*) map1(n,m)
  write (*,*) map2(n,m)
  write (*,*) map3(n)
  write (*,*)

  map3(n) = map2(n,m) * 40. + map1(n,m) * 3.

  write (*,*) map1(n,m)
  write (*,*) map2(n,m)
  write (*,*) map3(n)
  write (*,*)

end program test
  
#endif

! =========== EOF =======================================================

