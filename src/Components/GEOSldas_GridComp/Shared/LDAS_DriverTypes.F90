
module LDAS_DriverTypes
  
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
!
!GOSWIN
!
     real :: DUDP001              ! below all units are                     [kg/m2/s]
     real :: DUDP002
     real :: DUDP003
     real :: DUDP004
     real :: DUDP005
     real :: DUSV001
     real :: DUSV002
     real :: DUSV003
     real :: DUSV004
     real :: DUSV005
     real :: DUWT001
     real :: DUWT002
     real :: DUWT003
     real :: DUWT004
     real :: DUWT005
     real :: DUSD001
     real :: DUSD002
     real :: DUSD003
     real :: DUSD004
     real :: DUSD005
     real :: BCDP001
     real :: BCDP002
     real :: BCSV001
     real :: BCSV002
     real :: BCWT001
     real :: BCWT002
     real :: BCSD001
     real :: BCSD002
     real :: OCDP001
     real :: OCDP002
     real :: OCSV001
     real :: OCSV002
     real :: OCWT001
     real :: OCWT002
     real :: OCSD001
     real :: OCSD002
     real :: SUDP003
     real :: SUSV003
     real :: SUWT003
     real :: SUSD003
     real :: SSDP001
     real :: SSDP002
     real :: SSDP003
     real :: SSDP004
     real :: SSDP005
     real :: SSSV001
     real :: SSSV002
     real :: SSSV003
     real :: SSSV004
     real :: SSSV005
     real :: SSWT001
     real :: SSWT002
     real :: SSWT003
     real :: SSWT004
     real :: SSWT005
     real :: SSSD001
     real :: SSSD002
     real :: SSSD003
     real :: SSSD004
     real :: SSSD005

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

   met_force%DUDP001     = scalar
   met_force%DUDP002     = scalar
   met_force%DUDP003     = scalar
   met_force%DUDP004     = scalar
   met_force%DUDP005     = scalar
   met_force%DUSV001     = scalar
   met_force%DUSV002     = scalar
   met_force%DUSV003     = scalar
   met_force%DUSV004     = scalar
   met_force%DUSV005     = scalar
   met_force%DUWT001     = scalar
   met_force%DUWT002     = scalar
   met_force%DUWT003     = scalar
   met_force%DUWT004     = scalar
   met_force%DUWT005     = scalar
   met_force%DUSD001     = scalar
   met_force%DUSD002     = scalar
   met_force%DUSD003     = scalar
   met_force%DUSD004     = scalar
   met_force%DUSD005     = scalar
   met_force%BCDP001     = scalar
   met_force%BCDP002     = scalar
   met_force%BCSV001     = scalar
   met_force%BCSV002     = scalar
   met_force%BCWT001     = scalar
   met_force%BCWT002     = scalar
   met_force%BCSD001     = scalar
   met_force%BCSD002     = scalar
   met_force%OCDP001     = scalar
   met_force%OCDP002     = scalar
   met_force%OCSV001     = scalar
   met_force%OCSV002     = scalar
   met_force%OCWT001     = scalar
   met_force%OCWT002     = scalar
   met_force%OCSD001     = scalar
   met_force%OCSD002     = scalar
   met_force%SUDP003     = scalar
   met_force%SUSV003     = scalar
   met_force%SUWT003     = scalar
   met_force%SUSD003     = scalar
   met_force%SSDP001     = scalar
   met_force%SSDP002     = scalar
   met_force%SSDP003     = scalar
   met_force%SSDP004     = scalar
   met_force%SSDP005     = scalar
   met_force%SSSV001     = scalar
   met_force%SSSV002     = scalar
   met_force%SSSV003     = scalar
   met_force%SSSV004     = scalar
   met_force%SSSV005     = scalar
   met_force%SSWT001     = scalar
   met_force%SSWT002     = scalar
   met_force%SSWT003     = scalar
   met_force%SSWT004     = scalar
   met_force%SSWT005     = scalar
   met_force%SSSD001     = scalar
   met_force%SSSD002     = scalar
   met_force%SSSD003     = scalar
   met_force%SSSD004     = scalar
   met_force%SSSD005     = scalar

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

   met_force_div_scalar%DUDP001     =    met_force%DUDP001   / scalar
   met_force_div_scalar%DUDP002     =    met_force%DUDP002   / scalar
   met_force_div_scalar%DUDP003     =    met_force%DUDP003   / scalar
   met_force_div_scalar%DUDP004     =    met_force%DUDP004   / scalar
   met_force_div_scalar%DUDP005     =    met_force%DUDP005   / scalar
   met_force_div_scalar%DUSV001     =    met_force%DUSV001   / scalar
   met_force_div_scalar%DUSV002     =    met_force%DUSV002   / scalar
   met_force_div_scalar%DUSV003     =    met_force%DUSV003   / scalar
   met_force_div_scalar%DUSV004     =    met_force%DUSV004   / scalar
   met_force_div_scalar%DUSV005     =    met_force%DUSV005   / scalar
   met_force_div_scalar%DUWT001     =    met_force%DUWT001   / scalar
   met_force_div_scalar%DUWT002     =    met_force%DUWT002   / scalar
   met_force_div_scalar%DUWT003     =    met_force%DUWT003   / scalar
   met_force_div_scalar%DUWT004     =    met_force%DUWT004   / scalar
   met_force_div_scalar%DUWT005     =    met_force%DUWT005   / scalar
   met_force_div_scalar%DUSD001     =    met_force%DUSD001   / scalar
   met_force_div_scalar%DUSD002     =    met_force%DUSD002   / scalar
   met_force_div_scalar%DUSD003     =    met_force%DUSD003   / scalar
   met_force_div_scalar%DUSD004     =    met_force%DUSD004   / scalar
   met_force_div_scalar%DUSD005     =    met_force%DUSD005   / scalar
   met_force_div_scalar%BCDP001     =    met_force%BCDP001   / scalar
   met_force_div_scalar%BCDP002     =    met_force%BCDP002   / scalar
   met_force_div_scalar%BCSV001     =    met_force%BCSV001   / scalar
   met_force_div_scalar%BCSV002     =    met_force%BCSV002   / scalar
   met_force_div_scalar%BCWT001     =    met_force%BCWT001   / scalar
   met_force_div_scalar%BCWT002     =    met_force%BCWT002   / scalar
   met_force_div_scalar%BCSD001     =    met_force%BCSD001   / scalar
   met_force_div_scalar%BCSD002     =    met_force%BCSD002   / scalar
   met_force_div_scalar%OCDP001     =    met_force%OCDP001   / scalar
   met_force_div_scalar%OCDP002     =    met_force%OCDP002   / scalar
   met_force_div_scalar%OCSV001     =    met_force%OCSV001   / scalar
   met_force_div_scalar%OCSV002     =    met_force%OCSV002   / scalar
   met_force_div_scalar%OCWT001     =    met_force%OCWT001   / scalar
   met_force_div_scalar%OCWT002     =    met_force%OCWT002   / scalar
   met_force_div_scalar%OCSD001     =    met_force%OCSD001   / scalar
   met_force_div_scalar%OCSD002     =    met_force%OCSD002   / scalar
   met_force_div_scalar%SUDP003     =    met_force%SUDP003   / scalar
   met_force_div_scalar%SUSV003     =    met_force%SUSV003   / scalar
   met_force_div_scalar%SUWT003     =    met_force%SUWT003   / scalar
   met_force_div_scalar%SUSD003     =    met_force%SUSD003   / scalar
   met_force_div_scalar%SSDP001     =    met_force%SSDP001   / scalar
   met_force_div_scalar%SSDP002     =    met_force%SSDP002   / scalar
   met_force_div_scalar%SSDP003     =    met_force%SSDP003   / scalar
   met_force_div_scalar%SSDP004     =    met_force%SSDP004   / scalar
   met_force_div_scalar%SSDP005     =    met_force%SSDP005   / scalar
   met_force_div_scalar%SSSV001     =    met_force%SSSV001   / scalar
   met_force_div_scalar%SSSV002     =    met_force%SSSV002   / scalar
   met_force_div_scalar%SSSV003     =    met_force%SSSV003   / scalar
   met_force_div_scalar%SSSV004     =    met_force%SSSV004   / scalar
   met_force_div_scalar%SSSV005     =    met_force%SSSV005   / scalar
   met_force_div_scalar%SSWT001     =    met_force%SSWT001   / scalar
   met_force_div_scalar%SSWT002     =    met_force%SSWT002   / scalar
   met_force_div_scalar%SSWT003     =    met_force%SSWT003   / scalar
   met_force_div_scalar%SSWT004     =    met_force%SSWT004   / scalar
   met_force_div_scalar%SSWT005     =    met_force%SSWT005   / scalar
   met_force_div_scalar%SSSD001     =    met_force%SSSD001   / scalar
   met_force_div_scalar%SSSD002     =    met_force%SSSD002   / scalar
   met_force_div_scalar%SSSD003     =    met_force%SSSD003   / scalar
   met_force_div_scalar%SSSD004     =    met_force%SSSD004   / scalar
   met_force_div_scalar%SSSD005     =    met_force%SSSD005   / scalar
        
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

   add_met_force%DUDP001     = met_force_1%DUDP001     + met_force_2%DUDP001
   add_met_force%DUDP002     = met_force_1%DUDP002     + met_force_2%DUDP002
   add_met_force%DUDP003     = met_force_1%DUDP003     + met_force_2%DUDP003
   add_met_force%DUDP004     = met_force_1%DUDP004     + met_force_2%DUDP004
   add_met_force%DUDP005     = met_force_1%DUDP005     + met_force_2%DUDP005
   add_met_force%DUSV001     = met_force_1%DUSV001     + met_force_2%DUSV001
   add_met_force%DUSV002     = met_force_1%DUSV002     + met_force_2%DUSV002
   add_met_force%DUSV003     = met_force_1%DUSV003     + met_force_2%DUSV003
   add_met_force%DUSV004     = met_force_1%DUSV004     + met_force_2%DUSV004
   add_met_force%DUSV005     = met_force_1%DUSV005     + met_force_2%DUSV005
   add_met_force%DUWT001     = met_force_1%DUWT001     + met_force_2%DUWT001
   add_met_force%DUWT002     = met_force_1%DUWT002     + met_force_2%DUWT002
   add_met_force%DUWT003     = met_force_1%DUWT003     + met_force_2%DUWT003
   add_met_force%DUWT004     = met_force_1%DUWT004     + met_force_2%DUWT004
   add_met_force%DUWT005     = met_force_1%DUWT005     + met_force_2%DUWT005
   add_met_force%DUSD001     = met_force_1%DUSD001     + met_force_2%DUSD001
   add_met_force%DUSD002     = met_force_1%DUSD002     + met_force_2%DUSD002
   add_met_force%DUSD003     = met_force_1%DUSD003     + met_force_2%DUSD003
   add_met_force%DUSD004     = met_force_1%DUSD004     + met_force_2%DUSD004
   add_met_force%DUSD005     = met_force_1%DUSD005     + met_force_2%DUSD005
   add_met_force%BCDP001     = met_force_1%BCDP001     + met_force_2%BCDP001
   add_met_force%BCDP002     = met_force_1%BCDP002     + met_force_2%BCDP002
   add_met_force%BCSV001     = met_force_1%BCSV001     + met_force_2%BCSV001
   add_met_force%BCSV002     = met_force_1%BCSV002     + met_force_2%BCSV002
   add_met_force%BCWT001     = met_force_1%BCWT001     + met_force_2%BCWT001
   add_met_force%BCWT002     = met_force_1%BCWT002     + met_force_2%BCWT002
   add_met_force%BCSD001     = met_force_1%BCSD001     + met_force_2%BCSD001
   add_met_force%BCSD002     = met_force_1%BCSD002     + met_force_2%BCSD002
   add_met_force%OCDP001     = met_force_1%OCDP001     + met_force_2%OCDP001
   add_met_force%OCDP002     = met_force_1%OCDP002     + met_force_2%OCDP002
   add_met_force%OCSV001     = met_force_1%OCSV001     + met_force_2%OCSV001
   add_met_force%OCSV002     = met_force_1%OCSV002     + met_force_2%OCSV002
   add_met_force%OCWT001     = met_force_1%OCWT001     + met_force_2%OCWT001
   add_met_force%OCWT002     = met_force_1%OCWT002     + met_force_2%OCWT002
   add_met_force%OCSD001     = met_force_1%OCSD001     + met_force_2%OCSD001
   add_met_force%OCSD002     = met_force_1%OCSD002     + met_force_2%OCSD002
   add_met_force%SUDP003     = met_force_1%SUDP003     + met_force_2%SUDP003
   add_met_force%SUSV003     = met_force_1%SUSV003     + met_force_2%SUSV003
   add_met_force%SUWT003     = met_force_1%SUWT003     + met_force_2%SUWT003
   add_met_force%SUSD003     = met_force_1%SUSD003     + met_force_2%SUSD003
   add_met_force%SSDP001     = met_force_1%SSDP001     + met_force_2%SSDP001
   add_met_force%SSDP002     = met_force_1%SSDP002     + met_force_2%SSDP002
   add_met_force%SSDP003     = met_force_1%SSDP003     + met_force_2%SSDP003
   add_met_force%SSDP004     = met_force_1%SSDP004     + met_force_2%SSDP004
   add_met_force%SSDP005     = met_force_1%SSDP005     + met_force_2%SSDP005
   add_met_force%SSSV001     = met_force_1%SSSV001     + met_force_2%SSSV001
   add_met_force%SSSV002     = met_force_1%SSSV002     + met_force_2%SSSV002
   add_met_force%SSSV003     = met_force_1%SSSV003     + met_force_2%SSSV003
   add_met_force%SSSV004     = met_force_1%SSSV004     + met_force_2%SSSV004
   add_met_force%SSSV005     = met_force_1%SSSV005     + met_force_2%SSSV005
   add_met_force%SSWT001     = met_force_1%SSWT001     + met_force_2%SSWT001
   add_met_force%SSWT002     = met_force_1%SSWT002     + met_force_2%SSWT002
   add_met_force%SSWT003     = met_force_1%SSWT003     + met_force_2%SSWT003
   add_met_force%SSWT004     = met_force_1%SSWT004     + met_force_2%SSWT004
   add_met_force%SSWT005     = met_force_1%SSWT005     + met_force_2%SSWT005
   add_met_force%SSSD001     = met_force_1%SSSD001     + met_force_2%SSSD001
   add_met_force%SSSD002     = met_force_1%SSSD002     + met_force_2%SSSD002
   add_met_force%SSSD003     = met_force_1%SSSD003     + met_force_2%SSSD003
   add_met_force%SSSD004     = met_force_1%SSSD004     + met_force_2%SSSD004
   add_met_force%SSSD005     = met_force_1%SSSD005     + met_force_2%SSSD005
    
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
  
end module LDAS_DriverTypes


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

