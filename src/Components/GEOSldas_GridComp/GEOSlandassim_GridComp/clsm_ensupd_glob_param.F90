
module clsm_ensupd_glob_param
  
  ! global parameters for CLSM ens driver update
  !
  ! must re-compile if any of these change
  !
  ! reichle, 18 Jul 2005
  ! reichle, 14 Apr 2006 - "update_type" now in nml
  
  use MAPL_ConstantsMod,  ONLY: MAPL_TICE
  
  implicit none

  private

  public :: N_obs_species_nml
  public :: N_state_per_cat
  public :: unitnum_obslog
  public :: scale_temp   
  public :: scale_qa     
  public :: scale_capac  
  public :: scale_catdef 
  public :: scale_rzexc  
  public :: scale_srfexc 
  public :: scale_ght1   
  public :: scale_ght2   
  public :: scale_ght3   
  public :: scale_ght4   
  public :: scale_ght5   
  public :: scale_ght6   
  public :: scale_wesn   
  public :: scale_htsn   
  public :: scale_sndz   
  public :: FT_ANA_FT_THRESHOLD   
  public :: FT_ANA_LOWERBOUND_ASNOW
  public :: FT_ANA_LOWERBOUND_TEFF
  public :: FT_ANA_UPPERBOUND_TEFF

  public :: echo_clsm_ensupd_glob_param

  ! -----------------------------------------------------------------------
  !
  ! total number of all obs species defined in namelist file
  ! (regardless of whether "assim" flag is true or false)
  
  integer, parameter :: N_obs_species_nml = 52
  
  ! ----------------------------------------------------------------------
  !
  ! state and measurement dimensions
  
  integer, parameter :: N_state_per_cat = 25

  ! ----------------------------------------------------------------------
  !
  ! unit number for obslog file
  
  integer, parameter :: unitnum_obslog = 4567
  
  ! ----------------------------------------------------------------------
  !  
  ! scaling parameters for states
  !
  ! Catchment prognostic variables have very different scales.
  ! The parameters below attempt to have all elements of the 
  ! State vector vary on a comparable scale (roughly 0-1)

  real, parameter    :: scale_temp              = 273.16
                                                
  real, parameter    :: scale_qa                = 1.e-2 
                                                
  real, parameter    :: scale_capac             = .5     
                                                
  real, parameter    :: scale_catdef            = 500.
  real, parameter    :: scale_rzexc             = 50.
  real, parameter    :: scale_srfexc            = 5.
                                                
  ! for non-frozen conditions, ght(i) ~ 2.e6*dzgt(i)*tp(i) [deg C]
  !                                             
  ! assuming tp ~ 10 deg C, e.g. scale_ght1 = 2.e6 (with dzgt(1)~0.1 m)
                                                
  real, parameter    :: scale_ght1              = 2.e6
  real, parameter    :: scale_ght2              = 4.e6
  real, parameter    :: scale_ght3              = 1.e7
  real, parameter    :: scale_ght4              = 2.e7
  real, parameter    :: scale_ght5              = 3.e7
  real, parameter    :: scale_ght6              = 2.e8
                                                
  real, parameter    :: scale_wesn              = 1.     ! needs work
  real, parameter    :: scale_htsn              = 1.     ! needs work
  real, parameter    :: scale_sndz              = 1.     ! needs work
 
  ! ----------------------------------------------------------------
  !
  ! parameter for freeze/thaw (FT) analysis 

  real, parameter    :: FT_ANA_FT_THRESHOLD     =  0.5
  
  real, parameter    :: FT_ANA_LOWERBOUND_ASNOW =  0.15

  real, parameter    :: FT_ANA_LOWERBOUND_TEFF  = -1.0 + MAPL_TICE  ! [Kelvin]
  real, parameter    :: FT_ANA_UPPERBOUND_TEFF  = +1.0 + MAPL_TICE  ! [Kelvin]

  ! ----------------------------------------------------------------

contains
  
  subroutine echo_clsm_ensupd_glob_param(unitnum)

    implicit none
    
    integer, intent(in) :: unitnum

    ! echo all global parameters 
    
    write (unitnum,*)
    write (unitnum,*) '--------------------------------------------------------'
    write (unitnum,*) 'echo_clsm_ensupd_glob_param():'
    write (unitnum,*)
    write (unitnum,*) 'N_obs_species_nml       = ', N_obs_species_nml
    write (unitnum,*)                          
    write (unitnum,*) 'N_state_per_cat         = ', N_state_per_cat
    write (unitnum,*)                          
    write (unitnum,*) 'unitnum_obslog          = ', unitnum_obslog
    write (unitnum,*)                          
    write (unitnum,*) 'scale_temp              = ', scale_temp
    write (unitnum,*)                                                              
    write (unitnum,*) 'scale_qa                = ', scale_qa
    write (unitnum,*)                                                              
    write (unitnum,*) 'scale_capac             = ', scale_capac
    write (unitnum,*)                                                              
    write (unitnum,*) 'scale_catdef            = ', scale_catdef
    write (unitnum,*) 'scale_rzexc             = ', scale_rzexc
    write (unitnum,*) 'scale_srfexc            = ', scale_srfexc
    write (unitnum,*)                                                              
    write (unitnum,*) 'scale_ght1              = ', scale_ght1
    write (unitnum,*) 'scale_ght2              = ', scale_ght2
    write (unitnum,*) 'scale_ght3              = ', scale_ght3
    write (unitnum,*) 'scale_ght4              = ', scale_ght4
    write (unitnum,*) 'scale_ght5              = ', scale_ght5
    write (unitnum,*) 'scale_ght6              = ', scale_ght6
    write (unitnum,*)                                                         
    write (unitnum,*) 'scale_wesn              = ', scale_wesn
    write (unitnum,*) 'scale_htsn              = ', scale_htsn
    write (unitnum,*) 'scale_sndz              = ', scale_sndz
    write (unitnum,*)                                                         
                                               
    write (unitnum,*) 'FT_ANA_FT_THRESHOLD     = ', FT_ANA_FT_THRESHOLD     
    write (unitnum,*) 'FT_ANA_LOWERBOUND_ASNOW = ', FT_ANA_LOWERBOUND_ASNOW
    write (unitnum,*) 'FT_ANA_LOWERBOUND_TEFF  = ', FT_ANA_LOWERBOUND_TEFF  
    write (unitnum,*) 'FT_ANA_UPPERBOUND_TEFF  = ', FT_ANA_UPPERBOUND_TEFF  

    write (unitnum,*) 'end echo_clsm_ensupd_glob_param()'
    write (unitnum,*) '--------------------------------------------------------'
    write (unitnum,*)
        
  end subroutine echo_clsm_ensupd_glob_param
  
end module clsm_ensupd_glob_param


!======== EOF ==============================================================
