
module enkf_types  
  
  ! definition of types for the EnKF
  !
  ! reichle, 19 Jul 2005
  ! 
  ! added new field to obs_param_type for adaptive filtering 
  ! - reichle, 15 Dec 2006
  !
  ! revised "obs_param_type" for SMOS angles and downscaling ("FOV")
  ! - reichle, 13 Jun 2011
  !
  ! added "%bias_*" fields for obs bias estimation
  ! - reichle+draper, 28 Aug 2013
  !
  ! added "%coarsen_pert" to obs_param_type
  ! - reichle, 6 Dec 2013
  !
  ! reichle, 31 Jan 2014: added "%time" to "obs_type"
  !
  ! reichle, 8 Jun 2017: added "%flistpath" and "%flistname" to "obs_param_type"
  !
  ! -------------------------------------------------------------------
  
  implicit none

  save
  
  ! everything is private by default unless made public
  
  private
  
  public :: obs_type
  public :: obs_param_type
  public :: write_obs_param

  public :: N_obs_ang_max
  
  ! -------------------------------------------------------------------------
  !
  ! N_obs_ang_max = max # obs angles permitted per obs type in in nml file
  
  integer, parameter :: N_obs_ang_max = 7

  ! -----------------------------------------------------------------------
  
  ! obs_type is basic element of vector "Observations" (length N_obs), 
  ! which contains all observations of all types that are available
  ! at a given update time
  
  ! added innov information for adaptive filtering - reichle, 14 Dec 2006
  ! added OminusA information for adaptive filtering - reichle, 1 Feb 2007
  ! added "varname" field to "obs_param_type" - reichle 14 Jun 2011
  ! major revisions to "obs_type" fields - reichle 16 Jun 2011
  ! added "units" field to "obs_param_type" - reichle 22 Nov 2011

  type :: obs_type

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     logical :: assim     ! .T. if obs is assimilated (ie., used in state update), 
                          ! .F. if only "get_innov" or only used for obs bias update
     integer :: species   ! identifier for type of observation
     integer :: tilenum   ! number of tile within (full) domain
     integer :: DUMMYGAP  ! fill gap so that MPI STRUCT works with "-align" compiler flag
     real*8  :: time      ! time of obs (J2000 seconds w/ 'TT12' epoch; see date_time_util.F90)
     real    :: lon       ! longitude of obs
     real    :: lat       ! latitude of obs
     real    :: obs       ! observed value
     real    :: obsvar    ! obs error var
     real    :: fcst      ! "forecast": value of obs pred before EnKF update (ens mean)
     real    :: fcstvar   ! forecast error var (in obs space), a.k.a. HPHt
     real    :: ana       ! "analysis": value of obs pred after EnKF update (ens mean)
     real    :: anavar    ! analysis error var (in obs space), a.k.a. HAHt
     
  end type obs_type
  
  ! ----------------------------------------------------------------------
  !
  ! vector obs_param contains information about each species of observations 
  
  type :: obs_param_type

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

     character(40)                 :: descr     ! description
     integer                       :: species   ! identifier for type of measurement

     integer                       :: orbit     ! type of (half-)orbit
                                                !                     0 = n/a  [eg., in situ obs]
                                                !                     1 = ascending
                                                !                     2 = descending
                                                !                     3 = ascending or descending
                                                !                     4 = geostationary

     integer                       :: pol       ! polarization
                                                ! 0 = n/a  [eg., multi-pol. retrieval]
                                                ! 1 = horizontal
                                                ! 2 = vertical 
                                                ! 3 = ...
                                                ! [add 3rd/4th Stokes, HH, HV, VH, VV]

     integer                       :: N_ang     ! # angles in species

     real, &
          dimension(N_obs_ang_max) :: ang       ! vector of angles
     
     real                          :: freq      ! frequency [Hz]

     real                          :: FOV       ! field-of-view *radius* 
                                                ! if FOV==0. equate obs footprint w/ tile
                                                ! for details see LDASsa_DEFAULT_inputs ensupd.nml
     character(40)                 :: FOV_units ! FOV units ('km' or 'deg') 

     logical                       :: assim     ! assimilate yes/no? (see also "obs_type")
     logical                       :: scale     ! scale yes/no?
     logical                       :: getinnov  ! compute innovs? (.T. if assim==.T.)

     integer                       :: RTM_ID    ! ID of radiative transfer model 

     integer                       :: bias_Npar ! number of bias states tracked per day
     integer                       :: bias_trel ! e-folding time scale of obs bias memory [s]
     integer                       :: bias_tcut ! cutoff time for confident obs bias est [s]

     real                          :: nodata    ! no-data-value

     character(40)                 :: varname   ! equivalent model variable name (Obs_pred)
     character(40)                 :: units     ! units (eg., 'K' or 'm3/m3')

     character(200)                :: path      ! path to measurements file 
     character(80)                 :: name      ! name identifier for measurements 
     character(200)                :: scalepath ! path to file with scaling parameters
     character(80)                 :: scalename ! filename for scaling parameters
     character(200)                :: flistpath ! path to file with list of obs file names
     character(80)                 :: flistname ! name of file with list of obs file names

     real                          :: errstd    ! default obs error std
                                  
     real                          :: std_normal_max  ! see pert_param_type
     logical                       :: zeromean        ! see pert_param_type
     logical                       :: coarsen_pert    ! see pert_param_type ("%coarsen")
     real                          :: xcorr           ! see pert_param_type
     real                          :: ycorr           ! see pert_param_type
     
     integer                       :: adapt     ! identifier for adaptive filtering
     
  end type obs_param_type
  
  ! ----------------------------------------------------------------------

contains

  subroutine write_obs_param(unitnumber, N_obs_param, obs_param)

    implicit none
    
    integer,                                      intent(in) :: unitnumber
    integer,                                      intent(in) :: N_obs_param
    
    type(obs_param_type), dimension(N_obs_param), intent(in) :: obs_param
    
    ! local variables

    integer :: i
    
    ! --------------------------------------------------------------------

    write (unitnumber,*) N_obs_param

    do i=1,N_obs_param
       
       write (unitnumber, '(42A)') "'" // trim(obs_param(i)%descr)     // "'"
       write (unitnumber,       *) obs_param(i)%species   
       write (unitnumber,       *) obs_param(i)%orbit     
       write (unitnumber,       *) obs_param(i)%pol       
       write (unitnumber,       *) obs_param(i)%N_ang     
       
       write (unitnumber,       *) obs_param(i)%ang(1:obs_param(i)%N_ang)       
       
       write (unitnumber,       *) obs_param(i)%freq      
       write (unitnumber,       *) obs_param(i)%FOV
       write (unitnumber, '(42A)') "'" // trim(obs_param(i)%FOV_units) // "'"       
       write (unitnumber,       *) obs_param(i)%assim     
       write (unitnumber,       *) obs_param(i)%scale     
       write (unitnumber,       *) obs_param(i)%getinnov  
       write (unitnumber,       *) obs_param(i)%RTM_ID
       write (unitnumber,       *) obs_param(i)%bias_Npar
       write (unitnumber,       *) obs_param(i)%bias_trel
       write (unitnumber,       *) obs_param(i)%bias_tcut
       write (unitnumber,       *) obs_param(i)%nodata    
       write (unitnumber, '(42A)') "'" // trim(obs_param(i)%varname)   // "'"
       write (unitnumber, '(42A)') "'" // trim(obs_param(i)%units)     // "'"
       write (unitnumber,'(202A)') "'" // trim(obs_param(i)%path)      // "'"    
       write (unitnumber, '(82A)') "'" // trim(obs_param(i)%name)      // "'"      
       write (unitnumber,'(202A)') "'" // trim(obs_param(i)%scalepath) // "'" 
       write (unitnumber, '(82A)') "'" // trim(obs_param(i)%scalename) // "'" 
       write (unitnumber,'(202A)') "'" // trim(obs_param(i)%flistpath) // "'" 
       write (unitnumber, '(82A)') "'" // trim(obs_param(i)%flistname) // "'" 
       write (unitnumber,       *) obs_param(i)%errstd    
       write (unitnumber,       *) obs_param(i)%std_normal_max
       write (unitnumber,       *) obs_param(i)%zeromean
       write (unitnumber,       *) obs_param(i)%coarsen_pert
       write (unitnumber,       *) obs_param(i)%xcorr         
       write (unitnumber,       *) obs_param(i)%ycorr         
       write (unitnumber,       *) obs_param(i)%adapt     
       
    end do

  end subroutine write_obs_param
  
end module enkf_types

! ================== EOF ===============================================
