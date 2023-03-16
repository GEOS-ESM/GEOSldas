#include "MAPL_Generic.h"

!BOP
module GEOS_EnsGridCompMod

  ! !USES
  !! This grid comp behaves like a coupler. The set service, initialization are compliant with MAPL grid comp concept. 
  use ESMF
  use MAPL_Mod
  use catch_constants, only: DZGT => CATCH_DZGT
  use catch_types,     only: cat_progn_type
  use catch_types,     only: cat_param_type

  use, intrinsic :: ieee_arithmetic
  
  implicit none

  private

  public :: SetServices
  public :: catch_progn
  public :: catch_param

  ! !DESCRIPTION: This GridComp collects ensemble members and then averages the variables from Catchment.
  !               For select variables, the ensemble standard deviation is also computed.

  !EOP
  integer            :: NUM_ENSEMBLE
  integer            :: collect_land_counter
  integer            :: collect_force_counter
  integer, parameter :: NUM_SUBTILES=4
  real               :: enavg_nodata_threshold

  type(cat_progn_type),dimension(:,:), allocatable :: catch_progn
  type(cat_param_type),dimension(:  ), allocatable :: catch_param


contains

  !BOP

  ! !IROTUINE: SetServices -- Set ESMF services for this component

  ! !INTERFACE:

  subroutine SetServices(gc, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc ! gridded component
    integer, optional                  :: rc ! return code

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name


    ! Get my name and setup traceback handle
    Iam = 'SetServices'
    call ESMF_GridCompGet(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::" // Iam

    ! Register services for this component
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_INITIALIZE,                                                &
         Initialize,                                                            &
         rc=status                                                              &
         )
    VERIFY_(status)

    ! phase one: collect forcing ensemble 
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         Collect_force_ens,                                                     &
         rc=status                                                              &
         )
    VERIFY_(status)
 
   ! phase two : collect ensemble out from land
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         Collect_land_ens,                                                      &
         rc=status                                                              &
         )
    VERIFY_(status)

    !phase 3 : get cat_param
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         GET_CATCH_PARAM ,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)

    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_FINALIZE,                                                  &
         Finalize,                                                              &
         rc=status                                                              &
         )
    VERIFY_(status)

    ! Set the state variable specs
    !BOS

    ! !IMPORT STATE:

    ! this grid comp will take the other's export as import

    ! !EXPORT STATE:
    ! relay LAI to landassim
    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'LAI'                                     ,&
       LONG_NAME  = 'leaf_area_index'                         ,&
       UNITS      = '1'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)

!! exports for ens average (and, for a few variables, the ens std)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'canopy_temperature'        ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TC'                        ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'canopy_specific_humidity'  ,&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QC'                        ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'interception_reservoir_capac',&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CAPAC'                     ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'catchment_deficit'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CATDEF'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'root_zone_excess'          ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'RZEXC'                     ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'surface_excess'            ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'SRFEXC'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_1' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT1'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_2' ,&
    UNITS              = 'J_m-2'                     ,&
    SHORT_NAME         = 'GHTCNT2'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_3' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT3'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_4' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT4'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_5' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT5'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_6' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT6'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

!  call MAPL_AddExportSpec(GC                  ,&
!    LONG_NAME          = 'mean_catchment_temp_incl_snw',&
!    UNITS              = 'K'                         ,&
!    SHORT_NAME         = 'TSURF'                     ,&
!    DIMS               = MAPL_DimsTileOnly           ,&
!    VLOCATION          = MAPL_VLocationNone          ,&
!    RESTART            = RESTART_IN_FILE             ,&
!                                           RC=STATUS  )
!  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'snow_mass_layer_1'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN1'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'snow_mass_layer_2'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN2'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'snow_mass_layer_3'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN3'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'heat_content_snow_layer_1' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN1'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'heat_content_snow_layer_2' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN2'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'heat_content_snow_layer_3' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN3'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'snow_depth_layer_1'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN1'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'snow_depth_layer_2'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN2'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'snow_depth_layer_3'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN3'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

!  call MAPL_AddExportSpec(GC                  ,&
!    LONG_NAME          = 'surface_heat_exchange_coefficient',&
!    UNITS              = 'kg m-2 s-1'                ,&
!    SHORT_NAME         = 'CH'                        ,&
!    DIMS               = MAPL_DimsTileTile           ,&
!    NUM_SUBTILES       = NUM_SUBTILES                ,&
!    VLOCATION          = MAPL_VLocationNone          ,&
!    RESTART            = RESTART_IN_FILE             ,&
!                                           RC=STATUS  )
!  VERIFY_(STATUS)

!  call MAPL_AddExportSpec(GC                  ,&
!    LONG_NAME          = 'surface_momentum_exchange_coefficient',&
!    UNITS              = 'kg m-2 s-1'                ,&
!    SHORT_NAME         = 'CM'                        ,&
!    DIMS               = MAPL_DimsTileTile           ,&
!    NUM_SUBTILES       = NUM_SUBTILES                ,&
!    VLOCATION          = MAPL_VLocationNone          ,&
!    RESTART            = RESTART_IN_FILE             ,&
!                                           RC=STATUS  )
!  VERIFY_(STATUS)
!
!  call MAPL_AddExportSpec(GC                  ,&
!    LONG_NAME          = 'surface_moisture_exchange_coffiecient',&
!    UNITS              = 'kg m-2 s-1'                ,&
!    SHORT_NAME         = 'CQ'                        ,&
!    DIMS               = MAPL_DimsTileTile           ,&
!    NUM_SUBTILES       = NUM_SUBTILES                ,&
!    VLOCATION          = MAPL_VLocationNone          ,&
!    RESTART            = RESTART_IN_FILE             ,&
!                                           RC=STATUS  )
!  VERIFY_(STATUS)

!  call MAPL_AddExportSpec(GC                  ,&
!    LONG_NAME          = 'subtile_fractions'         ,&
!    UNITS              = '1'                         ,&
!    SHORT_NAME         = 'FR'                        ,&
!    DIMS               = MAPL_DimsTileTile           ,&
!    NUM_SUBTILES       = NUM_SUBTILES                ,&
!    VLOCATION          = MAPL_VLocationNone          ,&
!    RESTART            = RESTART_IN_FILE             ,&
!                                           RC=STATUS  )
!  VERIFY_(STATUS)

!  call MAPL_AddExportSpec(GC,                  &
!    SHORT_NAME         = 'WW',                        &
!    LONG_NAME          = 'vertical_velocity_scale_squared', &
!    UNITS              = 'm+2 s-2',                   &
!    DIMS               = MAPL_DimsTileTile,           &
!    NUM_SUBTILES       = NUM_SUBTILES                ,&
!    VLOCATION          = MAPL_VLocationNone,          &
!    RESTART            = RESTART_IN_FILE             ,&
!                                           RC=STATUS  )
!  VERIFY_(STATUS)

!  call MAPL_AddExportSpec(GC,                  &
!    SHORT_NAME         = 'DCH',                        &
!    LONG_NAME          = 'ch difference, optional in louissurface', &
!    UNITS              = '1',                   &
!    DIMS               = MAPL_DimsTileTile,           &
!    NUM_SUBTILES       = NUM_SUBTILES                ,&
!    VLOCATION          = MAPL_VLocationNone,          &
!    RESTART            = .false.                     ,&
!                                           RC=STATUS  )
!  VERIFY_(STATUS)
!
!  call MAPL_AddExportSpec(GC,                  &
!    SHORT_NAME         = 'DCQ',                        &
!    LONG_NAME          = 'cq difference, optional in louissurface', &
!    UNITS              = '1',                   &
!    DIMS               = MAPL_DimsTileTile,           &
!    NUM_SUBTILES       = NUM_SUBTILES                ,&
!    VLOCATION          = MAPL_VLocationNone,          &
!    RESTART            = .false.                     ,&
!                                           RC=STATUS  )
!  VERIFY_(STATUS)


!  !EXPORT STATE:

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'evaporation'               ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'EVAPOUT'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                     &
    LONG_NAME          = 'sublimation'               ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'SUBLIM'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'upward_sensible_heat_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SHOUT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'runoff_flux'               ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNOFF'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'interception_loss_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPINT'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baresoil_evap_energy_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSOI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'transpiration_energy_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPVEG'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_ice_evaporation_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPICE'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil moisture in Upper 10cm'     ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WAT10CM'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'totoal soil moisture'      ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WATSOI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil frozen water content' ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'ICESOI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowpack_evaporation_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSNO'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baseflow_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'BASEFLOW'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'overland_runoff_including_throughflow'  ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNSURF'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowmelt_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'SMELT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_outgoing_longwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'HLWUP'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    LONG_NAME          = 'surface_net_downward_longwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'LWNDSRF'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
    VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    LONG_NAME          = 'surface_net_downward_shortwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SWNDSRF'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
    VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'total_latent_energy_flux'  ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'HLATN'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'rainwater_infiltration_flux',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'QINFIL'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  !call MAPL_AddExportSpec(GC,                    &
  !  LONG_NAME          = 'areal_fraction_saturated_zone',&
  !  UNITS              = '1'                         ,&
  !  SHORT_NAME         = 'AR1'                       ,&
  !  DIMS               = MAPL_DimsTileOnly           ,&
  !  VLOCATION          = MAPL_VLocationNone          ,&
  !                                         RC=STATUS  )
  !VERIFY_(STATUS)

  !call MAPL_AddExportSpec(GC,                    &
  !  LONG_NAME          = 'areal_fraction_transpiration_zone',&
  !  UNITS              = '1'                         ,&
  !  SHORT_NAME         = 'AR2'                       ,&
  !  DIMS               = MAPL_DimsTileOnly           ,&
  !  VLOCATION          = MAPL_VLocationNone          ,&
  !                                         RC=STATUS  )
  !VERIFY_(STATUS)

  !call MAPL_AddExportSpec(GC,                    &
  !  LONG_NAME          = 'root_zone_equilibrium_moisture',&
  !  UNITS              = 'kg m-2'                    ,&
  !  SHORT_NAME         = 'RZEQ'                      ,&
  !  DIMS               = MAPL_DimsTileOnly           ,&
  !  VLOCATION          = MAPL_VLocationNone          ,&
  !                                         RC=STATUS  )
  !VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'ground_energy_flux'        ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'GHFLX'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'ave_catchment_temp_incl_snw',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSURF'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'ave_catchment_temp_incl_snw_ensstd',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSURF_ENSSTD'             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'temperature_top_snow_layer',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSNOW'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'temperature_unsaturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPUNST'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'temperature_saturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSAT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'temperature_wilted_zone'   ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPWLT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_land_snowcover',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ASNOW'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'downward_heat_flux_into_snow',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SHSNOW'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'averaged_snow_temperature' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'AVETSNOW'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_saturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRSAT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_unsaturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRUST'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_wilting_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRWLT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_mass'                 ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'SNOWMASS'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_depth'                ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNOWDP'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_soil_wetness'      ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET1'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'root_zone_soil_wetness'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET2'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'ave_prof_soil__moisture'   ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET3'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_surface_layer'       ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCSF'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_surface_layer_ensstd' ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCSF_ENSSTD'               ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_root_zone'           ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCRZ'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_root_zone_ensstd'    ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCRZ_ENSSTD'               ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_ave_prof'            ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCPR'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_ave_prof_ensstd'     ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCPR_ENSSTD'               ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_1' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TSOIL1TILE'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_1_ensstd' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TSOIL1TILE_ENSSTD'         ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_2' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TSOIL2TILE'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_3' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TSOIL3TILE'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_4' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TSOIL4TILE'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_5' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TSOIL5TILE'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_6' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TSOIL6TILE'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_emissivity'        ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'EMIS'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_albedo_visible_beam',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBVR'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_albedo_visible_diffuse',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBVF'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_albedo_near_infrared_beam',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBNR'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_albedo_near_infrared_diffuse',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBNF'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'change_surface_skin_temperature',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'DELTS'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'change_surface_specific_humidity',&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'DELQS'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  !call MAPL_AddExportSpec(GC,                    &
  !  LONG_NAME          = 'change_evaporation'        ,&
  !  UNITS              = 'kg m-2 s-1'                ,&
  !  SHORT_NAME         = 'DELEVAP'                   ,&
  !  DIMS               = MAPL_DimsTileOnly           ,&
  !  VLOCATION          = MAPL_VLocationNone          ,&
  !                                         RC=STATUS  )
  !VERIFY_(STATUS)

  !call MAPL_AddExportSpec(GC,                    &
  !  LONG_NAME          = 'change_upward_sensible_energy_flux',&
  !  UNITS              = 'W m-2'                     ,&
  !  SHORT_NAME         = 'DELSH'                     ,&
  !  DIMS               = MAPL_DimsTileOnly           ,&
  !  VLOCATION          = MAPL_VLocationNone          ,&
  !                                         RC=STATUS  )
  !VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_skin_temperature'  ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TST'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'land_surface_skin_temperature'  ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'LST'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_specific_humidity' ,&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QST'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'turbulence_surface_skin_temperature',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TH'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'turbulence_surface_skin_specific_hum',&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QH'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_heat_exchange_coefficient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CHT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_momentum_exchange_coefficient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CMT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_moisture_exchange_coefficient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CQT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'neutral_drag_coefficient'  ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_bulk_richardson_number',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'RIT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_roughness'         ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'Z0'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT2M',                     &
        LONG_NAME          = 'temperature 2m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ2M',                     &
        LONG_NAME          = 'humidity 2m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU2M',                    &
        LONG_NAME          = 'zonal 2m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV2M',                    &
        LONG_NAME          = 'meridional 2m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT10M',                     &
        LONG_NAME          = 'temperature 10m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ10M',                     &
        LONG_NAME          = 'humidity 10m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU10M',                    &
        LONG_NAME          = 'zonal 10m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV10M',                    &
        LONG_NAME          = 'meridional 10m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU50M',                    &
        LONG_NAME          = 'zonal 50m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV50M',                    &
        LONG_NAME          = 'meridional 50m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_roughness_for_heat',&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'Z0H'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'zero_plane_displacement_height',&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'D0'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GUST',                      &
    LONG_NAME          = 'gustiness',                 &
    UNITS              = 'm s-1',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'VENT',                      &
    LONG_NAME          = 'surface_ventilation_velocity',&
    UNITS              = 'm s-1',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'ACCUM',                             &
        LONG_NAME          = 'net_ice_accumulation_rate',         &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'EVLAND',                    &
    LONG_NAME          = 'Evaporation_land',          &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'PRLAND',                    &
    LONG_NAME          = 'Total_precipitation_land',  &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SNOLAND',                   &
    LONG_NAME          = 'snowfall_land',             &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DRPARLAND',                 &
    LONG_NAME          = 'surface_downwelling_par_beam_flux', &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DFPARLAND',                 &
    LONG_NAME          = 'surface_downwelling_par_diffuse_flux', &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LHSNOW',                    &
    LONG_NAME          = 'Latent_heat_flux_snow',     &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWNETSNOW',                    &
    LONG_NAME          = 'Net_shortwave_snow',        &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWUPSNOW',                    &
    LONG_NAME          = 'Net_longwave_snow',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWDNSNOW',                    &
    LONG_NAME          = 'Net_longwave_snow',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TCSORIG',                   &
    LONG_NAME          = 'Input_tc_for_snow',         &
    UNITS              = 'K',                         &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TPSN1IN',                   &
    LONG_NAME          = 'Input_temp_of_top_snow_lev',&
    UNITS              = 'K',                         &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TPSN1OUT',                  &
    LONG_NAME          = 'Output_temp_of_top_snow_lev',&
    UNITS              = 'K',                         &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHSNOW',                    &
    LONG_NAME          = 'Ground_heating_snow',       &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LHLAND',                    &
    LONG_NAME          = 'Latent_heat_flux_land',     &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SHLAND',                    &
    LONG_NAME          = 'Sensible_heat_flux_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWLAND',                    &
    LONG_NAME          = 'Net_shortwave_land',        &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWDOWNLAND',                &
    LONG_NAME          = 'Incident_shortwave_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWLAND',                    &
    LONG_NAME          = 'Net_longwave_land',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHLAND',                    &
    LONG_NAME          = 'Ground_heating_land',       &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHTSKIN',                   &
    LONG_NAME          = 'Ground_heating_skin_temp',  &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SMLAND',                    &
    LONG_NAME          = 'Snowmelt_flux_land',        &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TWLAND',                    &
    LONG_NAME          = 'Avail_water_storage_land',  &
    UNITS              = 'kg m-2',                    &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TELAND',                    &
    LONG_NAME          = 'Total_energy_storage_land', &
    UNITS              = 'J m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TSLAND',                    &
    LONG_NAME          = 'Total_snow_storage_land',   &
    UNITS              = 'kg m-2',                    &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DWLAND',                    &
    LONG_NAME          = 'rate_of_change_of_total_land_water',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DHLAND',                    &
    LONG_NAME          = 'rate_of_change_of_total_land_energy',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPLAND',                    &
    LONG_NAME          = 'rate_of_spurious_land_energy_source',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPWATR',                    &
    LONG_NAME          = 'rate_of_spurious_land_water_source',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPSNOW',                    &
    LONG_NAME          = 'rate_of_spurious_snow_energy',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'depth_to_water_table_from_surface_in_peat',&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'PEATCLSM_WATERLEVEL'       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RC=STATUS  ) 

  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'change_in_free_surface_water_reservoir_on_peat',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'PEATCLSM_FSWCHANGE'        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RC=STATUS  ) 
  
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_exposed_leaf-area_index',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNLAI'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_leaf-area_index'  ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNTLAI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_exposed_stem-area_index',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNSAI'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_carbon'           ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNTOTC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_vegetation_carbon',&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNVEGC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
 call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_net_primary_production' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNNPP'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_gross_primary_production',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNGPP'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_soil_respiration' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNSR'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_net_ecosystem_exchange' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNNEE'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'abstract_C_pool_to_meet_excess_MR_demand' ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNXSMR'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_added_to_maintain_positive_C' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNADD'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_carbon_loss_to_fire'    ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNLOSS'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_fractional_area_burn_rate' ,&
    UNITS              = 's-1'                       ,&
    SHORT_NAME         = 'CNBURN'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_root_C'           ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNROOT'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_fine_root_carbon'       ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNFROOTC'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'fire season length'        ,&
    UNITS              = 'days'                      ,&
    SHORT_NAME         = 'CNFSEL'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'absorbed_PAR'              ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'PARABS'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'incident_PAR'              ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'PARINC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

 call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'saturated_stomatal_conductance' ,&
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'SCSAT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'unstressed_stomatal_conductance' ,&
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'SCUNS'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'transpiration coefficient' ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BTRANT'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'solar induced fluorescence',&
    UNITS              = 'umol m-2 sm s-1'           ,&
    SHORT_NAME         = 'SIF'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)


  ! !EXPORT FORCING STATE:

   call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "TA",                                                 &
         LONG_NAME  = "perturbed_surface_air_temperature",                      &
         UNITS      = "K",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "QA",                                                 &
         LONG_NAME  = "perturbed_surface_air_specific_humidity",                &
         UNITS      = "kg kg-1",                                                &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PS",                                                 &
         LONG_NAME  = "surface_pressure",                             &
         UNITS      = "Pa",                                                     &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "UU",                                                 &
         LONG_NAME  = "perturbed_surface_wind_speed",                           &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "UWINDLMTILE",                                        &
         LONG_NAME  = "perturbed_levellm_uwind",                                &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "VWINDLMTILE",                                        &
         LONG_NAME  = "perturbed_levellm_vwind",                                &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)


    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PCU",                                                &
         LONG_NAME  = "perturbed_liquid_water_convective_precipitation",        &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PLS",                                                &
         LONG_NAME  = "perturbed_liquid_water_large_scale_precipitation",       &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "RainfSnowf",                                                  &
         LONG_NAME  = "rainf+snowf",                                         &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "SNO",                                                &
         LONG_NAME  = "perturbed_snowfall",                                     &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRPAR",                                              &
         LONG_NAME  = "surface_downwelling_par_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFPAR",                                              &
         LONG_NAME  = "surface_downwelling_par_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

   call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRNIR",                                              &
         LONG_NAME  = "surface_downwelling_nir_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFNIR",                                              &
         LONG_NAME  = "surface_downwelling_nir_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRUVR",                                              &
         LONG_NAME  = "surface_downwelling_uvr_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFUVR",                                              &
         LONG_NAME  = "surface_downwelling_uvr_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "LWDNSRF",                                            &
         LONG_NAME  = "perturbed_surface_downwelling_longwave_flux",            &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DZ",                                                   &
         LONG_NAME  = "reference_height_for_Tair_Qair_Wind",                    &
         UNITS      = "m",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    !EOS

    ! Set profiling timers
    call MAPL_TimerAdd(gc, name="Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Collect_force", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Collect_land", rc=status)
    VERIFY_(status)

    ! Call SetServices for children
    call MAPL_GenericSetServices(gc, rc=status)
    VERIFY_(status)
    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


  !BOP

  ! !IROTUINE: Initialize -- initialize method for LDAS GC

  ! !INTERFACE:

  subroutine Initialize(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(MAPL_LocStream) :: locstream
    integer :: land_nt_local

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Initialize"

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

   ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Initialize")

    call MAPL_GetResource ( MAPL, NUM_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)

    ! nodata handling for ensemble average
    _ASSERT( (MAPL_UNDEF>.9999*1.e15), Iam // ': nodata handling for ensemble average requires MAPL_UNDEF to be a very large number') 
    enavg_nodata_threshold   = MAPL_UNDEF/(NUM_ENSEMBLE+1)   

    collect_land_counter     = 0
    collect_force_counter    = 0

    ! Get number of land tiles
    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)
    call MAPL_LocStreamGet(locstream, NT_LOCAL=land_nt_local,rc=status)
    VERIFY_(status)

    allocate(catch_progn(land_nt_local, NUM_ENSEMBLE))
    allocate(catch_param(land_nt_local))

    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    call MAPL_TimerOff(MAPL, "Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")  

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


  !BOP

  ! !IROTUINE: collecting and averaging

  subroutine Collect_force_ens(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! ESMF variables
    type(ESMF_VM) :: vm

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null() ! MAPL obj

    ! Pointers to imports
    real, pointer :: TApert(:)=>null()
    real, pointer :: QApert(:)=>null()
    real, pointer :: PSpert(:)=>null()
    real, pointer :: UUpert(:)=>null()
    real, pointer :: PCUpert(:)=>null()
    real, pointer :: PLSpert(:)=>null()
    real, pointer :: SNOpert(:)=>null()
    real, pointer :: DRPARpert(:)=>null()
    real, pointer :: DFPARpert(:)=>null()
    real, pointer :: DRNIRpert(:)=>null()
    real, pointer :: DFNIRpert(:)=>null()
    real, pointer :: DRUVRpert(:)=>null()
    real, pointer :: DFUVRpert(:)=>null()
    real, pointer :: LWDNSRFpert(:)=>null()
    real, pointer :: DZpert(:)=>null()

    real, pointer :: TA_enavg(:)=>null()
    real, pointer :: QA_enavg(:)=>null()
    real, pointer :: PS_enavg(:)=>null()
    real, pointer :: UU_enavg(:)=>null()
    real, pointer :: PCU_enavg(:)=>null()
    real, pointer :: PLS_enavg(:)=>null()
    real, pointer :: SNO_enavg(:)=>null()
    real, pointer :: DRPAR_enavg(:)=>null()
    real, pointer :: DFPAR_enavg(:)=>null()
    real, pointer :: DRNIR_enavg(:)=>null()
    real, pointer :: DFNIR_enavg(:)=>null()
    real, pointer :: DRUVR_enavg(:)=>null()
    real, pointer :: DFUVR_enavg(:)=>null()
    real, pointer :: LWDNSRF_enavg(:)=>null()
    real, pointer :: DZ_enavg(:)=>null()
    real, pointer :: RainfSnowf(:)=>null()


   ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_collect_force_ens"

    if(NUM_ENSEMBLE ==1) then
    !   RETURN_(ESMF_SUCCESS)
    endif

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Collect_force")


    call MAPL_GetPointer(import, TApert, 'TApert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, QApert, 'QApert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PSpert, 'PSpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, UUpert, 'UUpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PCUpert, 'PCUpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PLSpert, 'PLSpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SNOpert, 'SNOpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRPARpert, 'DRPARpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFPARpert, 'DFPARpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRNIRpert, 'DRNIRpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFNIRpert, 'DFNIRpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRUVRpert, 'DRUVRpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFUVRpert, 'DFUVRpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LWDNSRFpert, 'LWDNSRFpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DZpert, 'DZpert', rc=status)
    VERIFY_(status)

   ! Pointers to exports (allocate memory)
    call MAPL_GetPointer(export, TA_enavg, 'TA', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QA_enavg, 'QA', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PS_enavg, 'PS', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, UU_enavg, 'UU', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PCU_enavg, 'PCU',alloc=.true.,  rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PLS_enavg, 'PLS', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNO_enavg, 'SNO', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRPAR_enavg, 'DRPAR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFPAR_enavg, 'DFPAR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRNIR_enavg, 'DRNIR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFNIR_enavg, 'DFNIR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRUVR_enavg, 'DRUVR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFUVR_enavg, 'DFUVR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LWDNSRF_enavg, 'LWDNSRF', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DZ_enavg, 'DZ', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RainfSnowf, 'RainfSnowf', alloc=.true., rc=status)
    VERIFY_(status)

    if(collect_force_counter ==0) then
       if(associated(TA_enavg))  TA_enavg = 0.0
       if(associated(QA_enavg))  QA_enavg = 0.0
       if(associated(PS_enavg))  PS_enavg = 0.0
       if(associated(UU_enavg))  UU_enavg = 0.0
       if(associated(PCU_enavg))  PCU_enavg = 0.0
       if(associated(PLS_enavg))  PLS_enavg = 0.0
       if(associated(SNO_enavg))  SNO_enavg = 0.0
       if(associated(DRPAR_enavg))  DRPAR_enavg = 0.0
       if(associated(DFPAR_enavg))  DFPAR_enavg = 0.0
       if(associated(DRNIR_enavg))  DRNIR_enavg = 0.0
       if(associated(DFNIR_enavg))  DFNIR_enavg = 0.0
       if(associated(DRUVR_enavg))  DRUVR_enavg = 0.0
       if(associated(DFUVR_enavg))  DFUVR_enavg = 0.0
       if(associated(LWDNSRF_enavg))  LWDNSRF_enavg = 0.0
       if(associated(DZ_enavg))  DZ_enavg = 0.0
    endif

    if(associated(TA_enavg))   &
       TA_enavg = TA_enavg + TApert
    if(associated(QA_enavg))   &
       QA_enavg = QA_enavg + QApert
    if(associated(PS_enavg))   &
       PS_enavg = PS_enavg + PSpert
    if(associated(UU_enavg))   &
       UU_enavg = UU_enavg + UUpert
    if(associated(PCU_enavg))   &
       PCU_enavg = PCU_enavg + PCUpert
    if(associated(PLS_enavg))   &
       PLS_enavg = PLS_enavg + PLSpert
    if(associated(SNO_enavg))   &
       SNO_enavg = SNO_enavg + SNOpert
    if(associated(DRPAR_enavg))   &
       DRPAR_enavg = DRPAR_enavg + DRPARpert
    if(associated(DFPAR_enavg))   &
       DFPAR_enavg = DFPAR_enavg + DFPARpert
    if(associated(DRNIR_enavg))   &
       DRNIR_enavg = DRNIR_enavg + DRNIRpert
    if(associated(DFNIR_enavg))   &
       DFNIR_enavg = DFNIR_enavg + DFNIRpert
    if(associated(DRUVR_enavg))   &
       DRUVR_enavg = DRUVR_enavg + DRUVRpert
    if(associated(DFUVR_enavg))   &
       DFUVR_enavg = DFUVR_enavg + DFUVRpert
    if(associated(LWDNSRF_enavg))   &
       LWDNSRF_enavg = LWDNSRF_enavg + LWDNSRFpert
    if(associated(DZ_enavg))   &
       DZ_enavg = DZ_enavg + DZpert

    collect_force_counter = collect_force_counter + 1
    if(collect_force_counter == NUM_ENSEMBLE) then
      collect_force_counter = 0

      if(associated(TA_enavg))  TA_enavg =TA_enavg /NUM_ENSEMBLE
      if(associated(QA_enavg))  QA_enavg =QA_enavg /NUM_ENSEMBLE
      if(associated(PS_enavg))  PS_enavg =PS_enavg /NUM_ENSEMBLE
      if(associated(UU_enavg))  UU_enavg =UU_enavg /NUM_ENSEMBLE
      if(associated(PCU_enavg))  PCU_enavg =PCU_enavg /NUM_ENSEMBLE
      if(associated(PLS_enavg))  PLS_enavg =PLS_enavg /NUM_ENSEMBLE
      if(associated(SNO_enavg))  SNO_enavg =SNO_enavg /NUM_ENSEMBLE
      if(associated(DRPAR_enavg))  DRPAR_enavg =DRPAR_enavg /NUM_ENSEMBLE
      if(associated(DFPAR_enavg))  DFPAR_enavg =DFPAR_enavg /NUM_ENSEMBLE
      if(associated(DRNIR_enavg))  DRNIR_enavg =DRNIR_enavg /NUM_ENSEMBLE
      if(associated(DFNIR_enavg))  DFNIR_enavg =DFNIR_enavg /NUM_ENSEMBLE
      if(associated(DRUVR_enavg))  DRUVR_enavg =DRUVR_enavg /NUM_ENSEMBLE
      if(associated(DFUVR_enavg))  DFUVR_enavg =DFUVR_enavg /NUM_ENSEMBLE
      if(associated(LWDNSRF_enavg))  LWDNSRF_enavg =LWDNSRF_enavg /NUM_ENSEMBLE
      if(associated(DZ_enavg))   DZ_enavg      =DZ_enavg /NUM_ENSEMBLE
      if(associated(RainfSnowf)) RainfSnowf    =PLS_enavg + PCU_enavg + SNO_enavg

    endif
    ! Turn timers off
    call MAPL_TimerOff(MAPL, "Collect_force")
    call MAPL_TimerOff(MAPL, "TOTAL")
       
   RETURN_(ESMF_SUCCESS)

  end subroutine Collect_force_ens

  subroutine Collect_land_ens(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! ESMF variables
    type(ESMF_VM) :: vm

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null() ! MAPL obj

    real, pointer :: poros(:) =>null()
    real, pointer :: cond(:) =>null()
    real, pointer :: psis(:) =>null()
    real, pointer :: bee(:) =>null()
    real, pointer :: wpwet(:) =>null()
    real, pointer :: gnu(:) =>null()
    real, pointer :: vgwmax(:) =>null()
    real, pointer :: bf1(:) =>null()
    real, pointer :: bf2(:) =>null()
    real, pointer :: bf3(:) =>null()
    real, pointer :: cdcr1(:) =>null()
    real, pointer :: cdcr2(:) =>null()
    real, pointer :: ars1(:) =>null()
    real, pointer :: ars2(:) =>null()
    real, pointer :: ars3(:) =>null()
    real, pointer :: ara1(:) =>null()
    real, pointer :: ara2(:) =>null()
    real, pointer :: ara3(:) =>null()
    real, pointer :: ara4(:) =>null()
    real, pointer :: arw1(:) =>null()
    real, pointer :: arw2(:) =>null()
    real, pointer :: arw3(:) =>null()
    real, pointer :: arw4(:) =>null()
    real, pointer :: tsa1(:) =>null()
    real, pointer :: tsa2(:) =>null()
    real, pointer :: tsb1(:) =>null()
    real, pointer :: tsb2(:) =>null()
    real, pointer :: atau(:) =>null()
    real, pointer :: btau(:) =>null()
    real, pointer :: ity(:) =>null()

    real, pointer :: in_lai(:) =>null()
    real, pointer :: out_lai(:) =>null()


    real, dimension(:,:),pointer :: TC,TC_enavg
    real, dimension(:,:),pointer :: QC,QC_enavg
    real, dimension(:),pointer :: CAPAC,CAPAC_enavg
    real, dimension(:),pointer :: CATDEF,CATDEF_enavg
    real, dimension(:),pointer :: RZEXC,RZEXC_enavg
    real, dimension(:),pointer :: SRFEXC,SRFEXC_enavg
    real, dimension(:),pointer :: GHTCNT1,GHTCNT1_enavg
    real, dimension(:),pointer :: GHTCNT2,GHTCNT2_enavg
    real, dimension(:),pointer :: GHTCNT3,GHTCNT3_enavg
    real, dimension(:),pointer :: GHTCNT4,GHTCNT4_enavg
    real, dimension(:),pointer :: GHTCNT5,GHTCNT5_enavg
    real, dimension(:),pointer :: GHTCNT6,GHTCNT6_enavg
    real, dimension(:),pointer :: WESNN1,WESNN1_enavg
    real, dimension(:),pointer :: WESNN2,WESNN2_enavg
    real, dimension(:),pointer :: WESNN3,WESNN3_enavg
    real, dimension(:),pointer :: HTSNNN1,HTSNNN1_enavg
    real, dimension(:),pointer :: HTSNNN2,HTSNNN2_enavg
    real, dimension(:),pointer :: HTSNNN3,HTSNNN3_enavg
    real, dimension(:),pointer :: SNDZN1,SNDZN1_enavg
    real, dimension(:),pointer :: SNDZN2,SNDZN2_enavg
    real, dimension(:),pointer :: SNDZN3,SNDZN3_enavg
    real, dimension(:),pointer :: EVAPOUT,EVAPOUT_enavg
    real, dimension(:),pointer :: SUBLIM,SUBLIM_enavg
    real, dimension(:),pointer :: SHOUT,SHOUT_enavg
    real, dimension(:),pointer :: RUNOFF,RUNOFF_enavg
    real, dimension(:),pointer :: EVPINT,EVPINT_enavg
    real, dimension(:),pointer :: EVPSOI,EVPSOI_enavg
    real, dimension(:),pointer :: EVPVEG,EVPVEG_enavg
    real, dimension(:),pointer :: EVPICE,EVPICE_enavg
    real, dimension(:),pointer :: WAT10CM,WAT10CM_enavg
    real, dimension(:),pointer :: WATSOI,WATSOI_enavg
    real, dimension(:),pointer :: ICESOI,ICESOI_enavg
    real, dimension(:),pointer :: EVPSNO,EVPSNO_enavg
    real, dimension(:),pointer :: BASEFLOW,BASEFLOW_enavg
    real, dimension(:),pointer :: RUNSURF,RUNSURF_enavg
    real, dimension(:),pointer :: SMELT,SMELT_enavg
    real, dimension(:),pointer :: HLWUP,HLWUP_enavg
    real, dimension(:),pointer :: LWNDSRF,LWNDSRF_enavg
    real, dimension(:),pointer :: SWNDSRF,SWNDSRF_enavg
    real, dimension(:),pointer :: HLATN,HLATN_enavg
    real, dimension(:),pointer :: QINFIL,QINFIL_enavg
    real, dimension(:),pointer :: GHFLX,GHFLX_enavg
    real, dimension(:),pointer :: TPSURF,TPSURF_enavg,TPSURF_enstd
    real, dimension(:),pointer :: TPSNOW,TPSNOW_enavg
    real, dimension(:),pointer :: TPUNST,TPUNST_enavg
    real, dimension(:),pointer :: TPSAT,TPSAT_enavg
    real, dimension(:),pointer :: TPWLT,TPWLT_enavg
    !real, dimension(:),pointer :: ASNOW,ASNOW_enavg
    real, dimension(:),pointer :: ASNOW_enavg
    real, dimension(:),pointer :: SHSNOW,SHSNOW_enavg
    real, dimension(:),pointer :: AVETSNOW,AVETSNOW_enavg
    real, dimension(:),pointer :: FRSAT,FRSAT_enavg
    real, dimension(:),pointer :: FRUST,FRUST_enavg
    real, dimension(:),pointer :: FRWLT,FRWLT_enavg
    real, dimension(:),pointer :: SNOWMASS,SNOWMASS_enavg
    real, dimension(:),pointer :: SNOWDP,SNOWDP_enavg
    real, dimension(:),pointer :: WET1,WET1_enavg
    real, dimension(:),pointer :: WET2,WET2_enavg
    real, dimension(:),pointer :: WET3,WET3_enavg
    real, dimension(:),pointer :: WCSF,WCSF_enavg,WCSF_enstd
    real, dimension(:),pointer :: WCRZ,WCRZ_enavg,WCRZ_enstd
    real, dimension(:),pointer :: WCPR,WCPR_enavg,WCPR_enstd
    real, dimension(:),pointer :: TP1,TP1_enavg,TP1_enstd
    real, dimension(:),pointer :: TP2,TP2_enavg
    real, dimension(:),pointer :: TP3,TP3_enavg
    real, dimension(:),pointer :: TP4,TP4_enavg
    real, dimension(:),pointer :: TP5,TP5_enavg
    real, dimension(:),pointer :: TP6,TP6_enavg
    real, dimension(:),pointer :: EMIS,EMIS_enavg
    real, dimension(:),pointer :: ALBVR,ALBVR_enavg
    real, dimension(:),pointer :: ALBVF,ALBVF_enavg
    real, dimension(:),pointer :: ALBNR,ALBNR_enavg
    real, dimension(:),pointer :: ALBNF,ALBNF_enavg
    real, dimension(:),pointer :: DELTS,DELTS_enavg
    real, dimension(:),pointer :: DELQS,DELQS_enavg
    real, dimension(:),pointer :: TST,TST_enavg
    real, dimension(:),pointer :: LST,LST_enavg
    real, dimension(:),pointer :: QST,QST_enavg
    real, dimension(:),pointer :: TH,TH_enavg
    real, dimension(:),pointer :: QH,QH_enavg
    real, dimension(:),pointer :: CHT,CHT_enavg
    real, dimension(:),pointer :: CMT,CMT_enavg
    real, dimension(:),pointer :: CQT,CQT_enavg
    real, dimension(:),pointer :: CNT,CNT_enavg
    real, dimension(:),pointer :: RIT,RIT_enavg
    real, dimension(:),pointer :: Z0,Z0_enavg
    real, dimension(:),pointer :: MOT2M,MOT2M_enavg
    real, dimension(:),pointer :: MOQ2M,MOQ2M_enavg
    real, dimension(:),pointer :: MOU2M,MOU2M_enavg
    real, dimension(:),pointer :: MOV2M,MOV2M_enavg
    real, dimension(:),pointer :: MOT10M,MOT10M_enavg
    real, dimension(:),pointer :: MOQ10M,MOQ10M_enavg
    real, dimension(:),pointer :: MOU10M,MOU10M_enavg
    real, dimension(:),pointer :: MOV10M,MOV10M_enavg
    real, dimension(:),pointer :: MOU50M,MOU50M_enavg
    real, dimension(:),pointer :: MOV50M,MOV50M_enavg
    real, dimension(:),pointer :: Z0H,Z0H_enavg
    real, dimension(:),pointer :: D0,D0_enavg
    real, dimension(:),pointer :: GUST,GUST_enavg
    real, dimension(:),pointer :: VENT,VENT_enavg
    real, dimension(:),pointer :: ACCUM,ACCUM_enavg
    real, dimension(:),pointer :: EVLAND,EVLAND_enavg
    real, dimension(:),pointer :: PRLAND,PRLAND_enavg
    real, dimension(:),pointer :: SNOLAND,SNOLAND_enavg
    real, dimension(:),pointer :: DRPARLAND,DRPARLAND_enavg
    real, dimension(:),pointer :: DFPARLAND,DFPARLAND_enavg
    real, dimension(:),pointer :: LHSNOW,LHSNOW_enavg
    real, dimension(:),pointer :: SWNETSNOW,SWNETSNOW_enavg
    real, dimension(:),pointer :: LWUPSNOW,LWUPSNOW_enavg
    real, dimension(:),pointer :: LWDNSNOW,LWDNSNOW_enavg
    real, dimension(:),pointer :: TCSORIG,TCSORIG_enavg
    real, dimension(:),pointer :: TPSN1IN,TPSN1IN_enavg
    real, dimension(:),pointer :: TPSN1OUT,TPSN1OUT_enavg
    real, dimension(:),pointer :: GHSNOW,GHSNOW_enavg
    real, dimension(:),pointer :: LHLAND,LHLAND_enavg
    real, dimension(:),pointer :: SHLAND,SHLAND_enavg
    real, dimension(:),pointer :: SWLAND,SWLAND_enavg
    real, dimension(:),pointer :: SWDOWNLAND,SWDOWNLAND_enavg
    real, dimension(:),pointer :: LWLAND,LWLAND_enavg
    real, dimension(:),pointer :: GHLAND,GHLAND_enavg
    real, dimension(:),pointer :: GHTSKIN,GHTSKIN_enavg
    real, dimension(:),pointer :: SMLAND,SMLAND_enavg
    real, dimension(:),pointer :: TWLAND,TWLAND_enavg
    real, dimension(:),pointer :: TELAND,TELAND_enavg
    real, dimension(:),pointer :: TSLAND,TSLAND_enavg
    real, dimension(:),pointer :: DWLAND,DWLAND_enavg
    real, dimension(:),pointer :: DHLAND,DHLAND_enavg
    real, dimension(:),pointer :: SPLAND,SPLAND_enavg
    real, dimension(:),pointer :: SPWATR,SPWATR_enavg
    real, dimension(:),pointer :: SPSNOW,SPSNOW_enavg
    real, dimension(:),pointer :: PEATCLSM_WATERLEVEL,PEATCLSM_WATERLEVEL_enavg
    real, dimension(:),pointer :: PEATCLSM_FSWCHANGE, PEATCLSM_FSWCHANGE_enavg

    real, dimension(:), pointer :: CNLAI,    CNLAI_enavg 
    real, dimension(:), pointer :: CNTLAI,  CNTLAI_enavg
    real, dimension(:), pointer :: CNSAI,    CNSAI_enavg
    real, dimension(:), pointer :: CNTOTC,  CNTOTC_enavg 
    real, dimension(:), pointer :: CNVEGC,  CNVEGC_enavg
    real, dimension(:), pointer :: CNROOT,  CNROOT_enavg
    real, dimension(:), pointer :: CNFROOTC,  CNFROOTC_enavg
    real, dimension(:), pointer :: CNNPP,    CNNPP_enavg
    real, dimension(:), pointer :: CNGPP,    CNGPP_enavg
    real, dimension(:), pointer :: CNSR,      CNSR_enavg
    real, dimension(:), pointer :: CNNEE,    CNNEE_enavg
    real, dimension(:), pointer :: CNXSMR,  CNXSMR_enavg
    real, dimension(:), pointer :: CNADD,    CNADD_enavg
    real, dimension(:), pointer :: PARABS,  PARABS_enavg
    real, dimension(:), pointer :: PARINC,  PARINC_enavg
    real, dimension(:), pointer :: SCSAT,    SCSAT_enavg
    real, dimension(:), pointer :: SCUNS,    SCUNS_enavg
    real, dimension(:), pointer :: BTRANT,  BTRANT_enavg
    real, dimension(:), pointer :: SIF,        SIF_enavg
    real, dimension(:), pointer :: CNLOSS,  CNLOSS_enavg
    real, dimension(:), pointer :: CNBURN,  CNBURN_enavg
    real, dimension(:), pointer :: CNFSEL,  CNFSEL_enavg

    real ::  Nm1, NdivNm1    

    ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_ens_averaging"
   
    if(NUM_ENSEMBLE ==1) then
       !RETURN_(ESMF_SUCCESS)
    endif

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Collect_land")

    call MAPL_GetPointer(import, poros,  'POROS' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, cond,  'COND' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, psis,  'PSIS' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, bee,  'BEE' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, wpwet,  'WPWET' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, gnu,  'GNU' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, vgwmax,  'VGWMAX' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, bf1,  'BF1' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, bf2,  'BF2' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, bf3,  'BF3' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, cdcr1,  'CDCR1' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, cdcr2,  'CDCR2' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ars1,  'ARS1' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ars2,  'ARS2' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ars3,  'ARS3' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ara1,  'ARA1' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ara2,  'ARA2' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ara3,  'ARA3' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ara4,  'ARA4' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, arw1,  'ARW1' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, arw2,  'ARW2' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, arw3,  'ARW3' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, arw4,  'ARW4' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, tsa1,  'TSA1' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, tsa2,  'TSA2' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, tsb1,  'TSB1' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, tsb2,  'TSB2' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, atau,  'ATAU' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, btau,  'BTAU' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ity,  'ITY' , rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, in_lai,  'LAI' , rc=status)
    VERIFY_(status)

    call MAPL_GetPointer(import, TC,  'TC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, QC,  'QC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, CAPAC,  'CAPAC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, CATDEF,  'CATDEF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, RZEXC,  'RZEXC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SRFEXC,  'SRFEXC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHTCNT1,  'GHTCNT1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHTCNT2,  'GHTCNT2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHTCNT3,  'GHTCNT3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHTCNT4,  'GHTCNT4' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHTCNT5,  'GHTCNT5' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHTCNT6,  'GHTCNT6' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WESNN1,  'WESNN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WESNN2,  'WESNN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WESNN3,  'WESNN3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, HTSNNN1,  'HTSNNN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, HTSNNN2,  'HTSNNN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, HTSNNN3,  'HTSNNN3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SNDZN1,  'SNDZN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SNDZN2,  'SNDZN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SNDZN3,  'SNDZN3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, EVAPOUT,  'EVAPOUT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SUBLIM,  'SUBLIM' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SHOUT,  'SHOUT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, RUNOFF,  'RUNOFF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, EVPINT,  'EVPINT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, EVPSOI,  'EVPSOI' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, EVPVEG,  'EVPVEG' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, EVPICE,  'EVPICE' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WAT10CM,  'WAT10CM' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WATSOI,  'WATSOI' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ICESOI,  'ICESOI' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, EVPSNO,  'EVPSNO' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, BASEFLOW,  'BASEFLOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, RUNSURF,  'RUNSURF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SMELT,  'SMELT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, HLWUP,  'HLWUP' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LWNDSRF,  'LWNDSRF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SWNDSRF,  'SWNDSRF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, HLATN,  'HLATN' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, QINFIL,  'QINFIL' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHFLX,  'GHFLX' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TPSURF,  'TPSURF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TPSNOW,  'TPSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TPUNST,  'TPUNST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TPSAT,  'TPSAT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TPWLT,  'TPWLT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SHSNOW,  'SHSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, AVETSNOW,  'AVETSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, FRSAT,  'FRSAT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, FRUST,  'FRUST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, FRWLT,  'FRWLT' ,rc=status)
    VERIFY_(status)
    ! for offline model , there is no 'ASNOW', recompute 
    ! call MAPL_GetPointer(import, ASNOW,  'ASNOW' ,rc=status)
    ! VERIFY_(status)
    call MAPL_GetPointer(import, SNOWMASS,  'SNOWMASS' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SNOWDP,  'SNOWDP' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WET1,  'WET1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WET2,  'WET2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WET3,  'WET3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WCSF,  'WCSF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WCRZ,  'WCRZ' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, WCPR,  'WCPR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TP1,  'TP1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TP2,  'TP2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TP3,  'TP3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TP4,  'TP4' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TP5,  'TP5' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TP6,  'TP6' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, EMIS,  'EMIS' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ALBVR,  'ALBVR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ALBVF,  'ALBVF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ALBNR,  'ALBNR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ALBNF,  'ALBNF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DELTS,  'DELTS' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DELQS,  'DELQS' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TST,  'TST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LST,  'LST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, QST,  'QST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TH,  'TH' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, QH,  'QH' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, CHT,  'CHT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, CMT,  'CMT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, CQT,  'CQT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, CNT,  'CNT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, RIT,  'RIT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, Z0,  'Z0' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOT2M,  'MOT2M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOQ2M,  'MOQ2M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOU2M,  'MOU2M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOV2M,  'MOV2M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOT10M,  'MOT10M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOQ10M,  'MOQ10M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOU10M,  'MOU10M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOV10M,  'MOV10M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOU50M,  'MOU50M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, MOV50M,  'MOV50M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, Z0H,  'Z0H' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, D0,  'D0' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GUST,  'GUST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, VENT,  'VENT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ACCUM,  'ACCUM' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, EVLAND,  'EVLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PRLAND,  'PRLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SNOLAND,  'SNOLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRPARLAND,  'DRPARLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFPARLAND,  'DFPARLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LHSNOW,  'LHSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SWNETSNOW,  'SWNETSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LWUPSNOW,  'LWUPSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LWDNSNOW,  'LWDNSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TCSORIG,  'TCSORIG' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TPSN1IN,  'TPSN1IN' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TPSN1OUT,  'TPSN1OUT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHSNOW,  'GHSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LHLAND,  'LHLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SHLAND,  'SHLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SWLAND,  'SWLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SWDOWNLAND,  'SWDOWNLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LWLAND,  'LWLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHLAND,  'GHLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, GHTSKIN,  'GHTSKIN' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SMLAND,  'SMLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TWLAND,  'TWLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TELAND,  'TELAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, TSLAND,  'TSLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DWLAND,  'DWLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DHLAND,  'DHLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SPLAND,  'SPLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SPWATR,  'SPWATR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SPSNOW,  'SPSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PEATCLSM_WATERLEVEL,'PEATCLSM_WATERLEVEL' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PEATCLSM_FSWCHANGE, 'PEATCLSM_FSWCHANGE'  ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ITY,  'ITY' ,rc=status)
    VERIFY_(status)

    ! CatchCN-specific variables (not available in standard Catch)
    
    call MAPL_GetPointer(import, CNLAI  ,   'CNLAI' ,   notFoundOK=.true.,  _RC)   
    call MAPL_GetPointer(import, CNTLAI ,   'CNTLAI',   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNSAI  ,   'CNSAI' ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNTOTC ,   'CNTOTC',   notFoundOK=.true.,  _RC) 
    call MAPL_GetPointer(import, CNVEGC ,   'CNVEGC',   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNROOT ,   'CNROOT',   notFoundOK=.true.,  _RC)    ! CatchCNCLM45 only
    call MAPL_GetPointer(import, CNFROOTC , 'CNFROOTC', notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNNPP  ,   'CNNPP' ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNGPP  ,   'CNGPP' ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNSR   ,   'CNSR'  ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNNEE  ,   'CNNEE' ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNXSMR ,   'CNXSMR',   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNADD  ,   'CNADD' ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, PARABS ,   'PARABS',   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, PARINC ,   'PARINC',   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, SCSAT  ,   'SCSAT' ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, SCUNS  ,   'SCUNS' ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, BTRANT ,   'BTRANT',   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, SIF    ,   'SIF'   ,   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNLOSS ,   'CNLOSS',   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNBURN ,   'CNBURN',   notFoundOK=.true.,  _RC)
    call MAPL_GetPointer(import, CNFSEL ,   'CNFSEL',   notFoundOK=.true.,  _RC)



    call MAPL_GetPointer(export, TC_enavg,  'TC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QC_enavg,  'QC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CAPAC_enavg,  'CAPAC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CATDEF_enavg,  'CATDEF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RZEXC_enavg,  'RZEXC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SRFEXC_enavg,  'SRFEXC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT1_enavg,  'GHTCNT1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT2_enavg,  'GHTCNT2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT3_enavg,  'GHTCNT3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT4_enavg,  'GHTCNT4' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT5_enavg,  'GHTCNT5' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT6_enavg,  'GHTCNT6' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN1_enavg,  'WESNN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN2_enavg,  'WESNN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN3_enavg,  'WESNN3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN1_enavg,  'HTSNNN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN2_enavg,  'HTSNNN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN3_enavg,  'HTSNNN3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN1_enavg,  'SNDZN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN2_enavg,  'SNDZN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN3_enavg,  'SNDZN3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, EVAPOUT_enavg,  'EVAPOUT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SUBLIM_enavg,  'SUBLIM' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SHOUT_enavg,  'SHOUT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RUNOFF_enavg,  'RUNOFF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, EVPINT_enavg,  'EVPINT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, EVPSOI_enavg,  'EVPSOI' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, EVPVEG_enavg,  'EVPVEG' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, EVPICE_enavg,  'EVPICE' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WAT10CM_enavg,  'WAT10CM' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WATSOI_enavg,  'WATSOI' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, ICESOI_enavg,  'ICESOI' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, EVPSNO_enavg,  'EVPSNO' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, BASEFLOW_enavg,  'BASEFLOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RUNSURF_enavg,  'RUNSURF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SMELT_enavg,  'SMELT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HLWUP_enavg,  'HLWUP' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LWNDSRF_enavg,  'LWNDSRF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SWNDSRF_enavg,  'SWNDSRF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HLATN_enavg,  'HLATN' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QINFIL_enavg,  'QINFIL' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHFLX_enavg,  'GHFLX' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TPSURF_enavg,  'TPSURF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TPSURF_enstd,  'TPSURF_ENSSTD' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TPSNOW_enavg,  'TPSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TPUNST_enavg,  'TPUNST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TPSAT_enavg,  'TPSAT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TPWLT_enavg,  'TPWLT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, ASNOW_enavg,  'ASNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SHSNOW_enavg,  'SHSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, AVETSNOW_enavg,  'AVETSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, FRSAT_enavg,  'FRSAT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, FRUST_enavg,  'FRUST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, FRWLT_enavg,  'FRWLT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNOWMASS_enavg,  'SNOWMASS' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNOWDP_enavg,  'SNOWDP' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WET1_enavg,  'WET1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WET2_enavg,  'WET2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WET3_enavg,  'WET3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WCSF_enavg,  'WCSF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WCSF_enstd,  'WCSF_ENSSTD' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WCRZ_enavg,  'WCRZ' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WCRZ_enstd,  'WCRZ_ENSSTD' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WCPR_enavg,  'WCPR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WCPR_enstd,  'WCPR_ENSSTD' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TP1_enavg,  'TSOIL1TILE' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TP1_enstd,  'TSOIL1TILE_ENSSTD' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TP2_enavg,  'TSOIL2TILE' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TP3_enavg,  'TSOIL3TILE' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TP4_enavg,  'TSOIL4TILE' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TP5_enavg,  'TSOIL5TILE' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TP6_enavg,  'TSOIL6TILE' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, EMIS_enavg,  'EMIS' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, ALBVR_enavg,  'ALBVR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, ALBVF_enavg,  'ALBVF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, ALBNR_enavg,  'ALBNR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, ALBNF_enavg,  'ALBNF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DELTS_enavg,  'DELTS' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DELQS_enavg,  'DELQS' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TST_enavg,  'TST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LST_enavg,  'LST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QST_enavg,  'QST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TH_enavg,  'TH' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QH_enavg,  'QH' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CHT_enavg,  'CHT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CMT_enavg,  'CMT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CQT_enavg,  'CQT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CNT_enavg,  'CNT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RIT_enavg,  'RIT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, Z0_enavg,  'Z0' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOT2M_enavg,  'MOT2M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOQ2M_enavg,  'MOQ2M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOU2M_enavg,  'MOU2M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOV2M_enavg,  'MOV2M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOT10M_enavg,  'MOT10M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOQ10M_enavg,  'MOQ10M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOU10M_enavg,  'MOU10M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOV10M_enavg,  'MOV10M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOU50M_enavg,  'MOU50M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, MOV50M_enavg,  'MOV50M' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, Z0H_enavg,  'Z0H' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, D0_enavg,  'D0' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GUST_enavg,  'GUST' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, VENT_enavg,  'VENT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, ACCUM_enavg,  'ACCUM' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, EVLAND_enavg,  'EVLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PRLAND_enavg,  'PRLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNOLAND_enavg,  'SNOLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRPARLAND_enavg,  'DRPARLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFPARLAND_enavg,  'DFPARLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LHSNOW_enavg,  'LHSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SWNETSNOW_enavg,  'SWNETSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LWUPSNOW_enavg,  'LWUPSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LWDNSNOW_enavg,  'LWDNSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TCSORIG_enavg,  'TCSORIG' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TPSN1IN_enavg,  'TPSN1IN' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TPSN1OUT_enavg,  'TPSN1OUT' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHSNOW_enavg,  'GHSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LHLAND_enavg,  'LHLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SHLAND_enavg,  'SHLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SWLAND_enavg,  'SWLAND',alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SWDOWNLAND_enavg,  'SWDOWNLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LWLAND_enavg,  'LWLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHLAND_enavg,  'GHLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTSKIN_enavg,  'GHTSKIN' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SMLAND_enavg,  'SMLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TWLAND_enavg,  'TWLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TELAND_enavg,  'TELAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TSLAND_enavg,  'TSLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DWLAND_enavg,  'DWLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DHLAND_enavg,  'DHLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SPLAND_enavg,  'SPLAND' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SPWATR_enavg,  'SPWATR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SPSNOW_enavg,  'SPSNOW' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PEATCLSM_WATERLEVEL_enavg,'PEATCLSM_WATERLEVEL' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PEATCLSM_FSWCHANGE_enavg, 'PEATCLSM_FSWCHANGE'  ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, out_lai,  'LAI' , alloc=.true., rc=status)
    VERIFY_(status)

    call MAPL_GetPointer(export,   CNLAI_enavg  ,   'CNLAI' , _RC)   
    call MAPL_GetPointer(export,  CNTLAI_enavg  ,   'CNTLAI', _RC)
    call MAPL_GetPointer(export,   CNSAI_enavg  ,   'CNSAI' , _RC)
    call MAPL_GetPointer(export,  CNTOTC_enavg  ,   'CNTOTC', _RC) 
    call MAPL_GetPointer(export,  CNVEGC_enavg  ,   'CNVEGC', _RC)
    call MAPL_GetPointer(export,  CNROOT_enavg  ,   'CNROOT', _RC)
    call MAPL_GetPointer(export,  CNFROOTC_enavg,   'CNFROOTC', _RC)
    call MAPL_GetPointer(export,   CNNPP_enavg  ,   'CNNPP' , _RC)
    call MAPL_GetPointer(export,   CNGPP_enavg  ,   'CNGPP' , _RC)
    call MAPL_GetPointer(export,    CNSR_enavg  ,   'CNSR'  , _RC)
    call MAPL_GetPointer(export,   CNNEE_enavg  ,   'CNNEE' , _RC)
    call MAPL_GetPointer(export,   CNXSMR_enavg ,   'CNXSMR', _RC)
    call MAPL_GetPointer(export,   CNADD_enavg  ,   'CNADD' , _RC)
    call MAPL_GetPointer(export,  PARABS_enavg  ,   'PARABS', _RC)
    call MAPL_GetPointer(export,  PARINC_enavg  ,   'PARINC', _RC)
    call MAPL_GetPointer(export,   SCSAT_enavg  ,   'SCSAT' , _RC)
    call MAPL_GetPointer(export,   SCUNS_enavg  ,   'SCUNS' , _RC)
    call MAPL_GetPointer(export,  BTRANT_enavg  ,   'BTRANT', _RC)
    call MAPL_GetPointer(export,     SIF_enavg  ,   'SIF'   , _RC)
    call MAPL_GetPointer(export,  CNLOSS_enavg  ,   'CNLOSS', _RC)
    call MAPL_GetPointer(export,  CNBURN_enavg  ,   'CNBURN', _RC)
    call MAPL_GetPointer(export,  CNFSEL_enavg  ,   'CNFSEL', _RC)

    out_lai = in_lai
    if (collect_land_counter == 0) then

        if(associated(TC_enavg)) TC_enavg =  0.0
        if(associated(QC_enavg)) QC_enavg =  0.0
        if(associated(CAPAC_enavg)) CAPAC_enavg =  0.0
        if(associated(CATDEF_enavg)) CATDEF_enavg =  0.0
        if(associated(RZEXC_enavg)) RZEXC_enavg =  0.0
        if(associated(SRFEXC_enavg)) SRFEXC_enavg =  0.0
        if(associated(GHTCNT1_enavg)) GHTCNT1_enavg =  0.0
        if(associated(GHTCNT2_enavg)) GHTCNT2_enavg =  0.0
        if(associated(GHTCNT3_enavg)) GHTCNT3_enavg =  0.0
        if(associated(GHTCNT4_enavg)) GHTCNT4_enavg =  0.0
        if(associated(GHTCNT5_enavg)) GHTCNT5_enavg =  0.0
        if(associated(GHTCNT6_enavg)) GHTCNT6_enavg =  0.0
        if(associated(WESNN1_enavg)) WESNN1_enavg =  0.0
        if(associated(WESNN2_enavg)) WESNN2_enavg =  0.0
        if(associated(WESNN3_enavg)) WESNN3_enavg =  0.0
        if(associated(HTSNNN1_enavg)) HTSNNN1_enavg =  0.0
        if(associated(HTSNNN2_enavg)) HTSNNN2_enavg =  0.0
        if(associated(HTSNNN3_enavg)) HTSNNN3_enavg =  0.0
        if(associated(SNDZN1_enavg)) SNDZN1_enavg =  0.0
        if(associated(SNDZN2_enavg)) SNDZN2_enavg =  0.0
        if(associated(SNDZN3_enavg)) SNDZN3_enavg =  0.0
        if(associated(EVAPOUT_enavg)) EVAPOUT_enavg =  0.0
        if(associated(SUBLIM_enavg)) SUBLIM_enavg =  0.0
        if(associated(SHOUT_enavg)) SHOUT_enavg =  0.0
        if(associated(RUNOFF_enavg)) RUNOFF_enavg =  0.0
        if(associated(EVPINT_enavg)) EVPINT_enavg =  0.0
        if(associated(EVPSOI_enavg)) EVPSOI_enavg =  0.0
        if(associated(EVPVEG_enavg)) EVPVEG_enavg =  0.0
        if(associated(EVPICE_enavg)) EVPICE_enavg =  0.0
        if(associated(WAT10CM_enavg)) WAT10CM_enavg =  0.0
        if(associated(WATSOI_enavg)) WATSOI_enavg =  0.0
        if(associated(ICESOI_enavg)) ICESOI_enavg =  0.0
        if(associated(EVPSNO_enavg)) EVPSNO_enavg =  0.0
        if(associated(BASEFLOW_enavg)) BASEFLOW_enavg =  0.0
        if(associated(RUNSURF_enavg)) RUNSURF_enavg =  0.0
        if(associated(SMELT_enavg)) SMELT_enavg =  0.0
        if(associated(HLWUP_enavg)) HLWUP_enavg =  0.0
        if(associated(LWNDSRF_enavg)) LWNDSRF_enavg =  0.0
        if(associated(SWNDSRF_enavg)) SWNDSRF_enavg =  0.0
        if(associated(HLATN_enavg)) HLATN_enavg =  0.0
        if(associated(QINFIL_enavg)) QINFIL_enavg =  0.0
        if(associated(GHFLX_enavg)) GHFLX_enavg =  0.0
        if(associated(TPSURF_enavg)) TPSURF_enavg =  0.0
        if(associated(TPSURF_enstd)) TPSURF_enstd =  0.0
        if(associated(TPSNOW_enavg)) TPSNOW_enavg =  0.0
        if(associated(TPUNST_enavg)) TPUNST_enavg =  0.0
        if(associated(TPSAT_enavg)) TPSAT_enavg =  0.0
        if(associated(TPWLT_enavg)) TPWLT_enavg =  0.0
        if(associated(ASNOW_enavg)) ASNOW_enavg =  0.0
        if(associated(SHSNOW_enavg)) SHSNOW_enavg =  0.0
        if(associated(AVETSNOW_enavg)) AVETSNOW_enavg =  0.0
        if(associated(FRSAT_enavg)) FRSAT_enavg =  0.0
        if(associated(FRUST_enavg)) FRUST_enavg =  0.0
        if(associated(FRWLT_enavg)) FRWLT_enavg =  0.0
        if(associated(SNOWMASS_enavg)) SNOWMASS_enavg =  0.0
        if(associated(SNOWDP_enavg)) SNOWDP_enavg =  0.0
        if(associated(WET1_enavg)) WET1_enavg =  0.0
        if(associated(WET2_enavg)) WET2_enavg =  0.0
        if(associated(WET3_enavg)) WET3_enavg =  0.0
        if(associated(WCSF_enavg)) WCSF_enavg =  0.0
        if(associated(WCSF_enstd)) WCSF_enstd =  0.0
        if(associated(WCRZ_enavg)) WCRZ_enavg =  0.0
        if(associated(WCRZ_enstd)) WCRZ_enstd =  0.0
        if(associated(WCPR_enavg)) WCPR_enavg =  0.0
        if(associated(WCPR_enstd)) WCPR_enstd =  0.0
        if(associated(TP1_enavg)) TP1_enavg =  0.0
        if(associated(TP1_enstd)) TP1_enstd =  0.0
        if(associated(TP2_enavg)) TP2_enavg =  0.0
        if(associated(TP3_enavg)) TP3_enavg =  0.0
        if(associated(TP4_enavg)) TP4_enavg =  0.0
        if(associated(TP5_enavg)) TP5_enavg =  0.0
        if(associated(TP6_enavg)) TP6_enavg =  0.0
        if(associated(EMIS_enavg)) EMIS_enavg =  0.0
        if(associated(ALBVR_enavg)) ALBVR_enavg =  0.0
        if(associated(ALBVF_enavg)) ALBVF_enavg =  0.0
        if(associated(ALBNR_enavg)) ALBNR_enavg =  0.0
        if(associated(ALBNF_enavg)) ALBNF_enavg =  0.0
        if(associated(DELTS_enavg)) DELTS_enavg =  0.0
        if(associated(DELQS_enavg)) DELQS_enavg =  0.0
        if(associated(TST_enavg)) TST_enavg =  0.0
        if(associated(LST_enavg)) LST_enavg =  0.0
        if(associated(QST_enavg)) QST_enavg =  0.0
        if(associated(TH_enavg)) TH_enavg =  0.0
        if(associated(QH_enavg)) QH_enavg =  0.0
        if(associated(CHT_enavg)) CHT_enavg =  0.0
        if(associated(CMT_enavg)) CMT_enavg =  0.0
        if(associated(CQT_enavg)) CQT_enavg =  0.0
        if(associated(CNT_enavg)) CNT_enavg =  0.0
        if(associated(RIT_enavg)) RIT_enavg =  0.0
        if(associated(Z0_enavg)) Z0_enavg =  0.0
        if(associated(MOT2M_enavg)) MOT2M_enavg =  0.0
        if(associated(MOQ2M_enavg)) MOQ2M_enavg =  0.0
        if(associated(MOU2M_enavg)) MOU2M_enavg =  0.0
        if(associated(MOV2M_enavg)) MOV2M_enavg =  0.0
        if(associated(MOT10M_enavg)) MOT10M_enavg =  0.0
        if(associated(MOQ10M_enavg)) MOQ10M_enavg =  0.0
        if(associated(MOU10M_enavg)) MOU10M_enavg =  0.0
        if(associated(MOV10M_enavg)) MOV10M_enavg =  0.0
        if(associated(MOU50M_enavg)) MOU50M_enavg =  0.0
        if(associated(MOV50M_enavg)) MOV50M_enavg =  0.0
        if(associated(Z0H_enavg)) Z0H_enavg =  0.0
        if(associated(D0_enavg)) D0_enavg =  0.0
        if(associated(GUST_enavg)) GUST_enavg =  0.0
        if(associated(VENT_enavg)) VENT_enavg =  0.0
        if(associated(ACCUM_enavg)) ACCUM_enavg =  0.0
        if(associated(EVLAND_enavg)) EVLAND_enavg =  0.0
        if(associated(PRLAND_enavg)) PRLAND_enavg =  0.0
        if(associated(SNOLAND_enavg)) SNOLAND_enavg =  0.0
        if(associated(DRPARLAND_enavg)) DRPARLAND_enavg =  0.0
        if(associated(DFPARLAND_enavg)) DFPARLAND_enavg =  0.0
        if(associated(LHSNOW_enavg)) LHSNOW_enavg =  0.0
        if(associated(SWNETSNOW_enavg)) SWNETSNOW_enavg =  0.0
        if(associated(LWUPSNOW_enavg)) LWUPSNOW_enavg =  0.0
        if(associated(LWDNSNOW_enavg)) LWDNSNOW_enavg =  0.0
        if(associated(TCSORIG_enavg)) TCSORIG_enavg =  0.0
        if(associated(TPSN1IN_enavg)) TPSN1IN_enavg =  0.0
        if(associated(TPSN1OUT_enavg)) TPSN1OUT_enavg =  0.0
        if(associated(GHSNOW_enavg)) GHSNOW_enavg =  0.0
        if(associated(LHLAND_enavg)) LHLAND_enavg =  0.0
        if(associated(SHLAND_enavg)) SHLAND_enavg =  0.0
        if(associated(SWLAND_enavg)) SWLAND_enavg =  0.0
        if(associated(SWDOWNLAND_enavg)) SWDOWNLAND_enavg =  0.0
        if(associated(LWLAND_enavg)) LWLAND_enavg =  0.0
        if(associated(GHLAND_enavg)) GHLAND_enavg =  0.0
        if(associated(GHTSKIN_enavg)) GHTSKIN_enavg =  0.0
        if(associated(SMLAND_enavg)) SMLAND_enavg =  0.0
        if(associated(TWLAND_enavg)) TWLAND_enavg =  0.0
        if(associated(TELAND_enavg)) TELAND_enavg =  0.0
        if(associated(TSLAND_enavg)) TSLAND_enavg =  0.0
        if(associated(DWLAND_enavg)) DWLAND_enavg =  0.0
        if(associated(DHLAND_enavg)) DHLAND_enavg =  0.0
        if(associated(SPLAND_enavg)) SPLAND_enavg =  0.0
        if(associated(SPWATR_enavg)) SPWATR_enavg =  0.0
        if(associated(SPSNOW_enavg)) SPSNOW_enavg =  0.0
        if(associated(PEATCLSM_WATERLEVEL_enavg)) PEATCLSM_WATERLEVEL_enavg =  0.0
        if(associated(PEATCLSM_FSWCHANGE_enavg))  PEATCLSM_FSWCHANGE_enavg  =  0.0

        if(associated(    CNLAI_enavg))      CNLAI_enavg = 0.0 
        if(associated(   CNTLAI_enavg))     CNTLAI_enavg = 0.0
        if(associated(    CNSAI_enavg))      CNSAI_enavg = 0.0
        if(associated(   CNTOTC_enavg))     CNTOTC_enavg = 0.0
        if(associated(   CNVEGC_enavg))     CNVEGC_enavg = 0.0
        if(associated(   CNROOT_enavg))     CNROOT_enavg = 0.0
        if(associated(   CNFROOTC_enavg))     CNFROOTC_enavg = 0.0
        if(associated(    CNNPP_enavg))      CNNPP_enavg = 0.0
        if(associated(    CNGPP_enavg))      CNGPP_enavg = 0.0
        if(associated(     CNSR_enavg))       CNSR_enavg = 0.0
        if(associated(    CNNEE_enavg))      CNNEE_enavg = 0.0
        if(associated(    CNXSMR_enavg))    CNXSMR_enavg = 0.0
        if(associated(    CNADD_enavg))      CNADD_enavg = 0.0
        if(associated(   PARABS_enavg))     PARABS_enavg = 0.0
        if(associated(   PARINC_enavg))     PARINC_enavg = 0.0
        if(associated(    SCSAT_enavg))      SCSAT_enavg = 0.0
        if(associated(    SCUNS_enavg))      SCUNS_enavg = 0.0
        if(associated(   BTRANT_enavg))     BTRANT_enavg = 0.0
        if(associated(      SIF_enavg))        SIF_enavg = 0.0
        if(associated(   CNLOSS_enavg))     CNLOSS_enavg = 0.0
        if(associated(   CNBURN_enavg))     CNBURN_enavg = 0.0
        if(associated(   CNFSEL_enavg))     CNFSEL_enavg = 0.0
    endif

    if(associated(TC_enavg) .and. associated(TC))   & 
        TC_enavg = TC_enavg + TC
    if(associated(QC_enavg) .and. associated(QC))   & 
        QC_enavg = QC_enavg + QC
    if(associated(CAPAC_enavg) .and. associated(CAPAC))   & 
        CAPAC_enavg = CAPAC_enavg + CAPAC
    if(associated(CATDEF_enavg) .and. associated(CATDEF))   & 
        CATDEF_enavg = CATDEF_enavg + CATDEF
    if(associated(RZEXC_enavg) .and. associated(RZEXC))   & 
        RZEXC_enavg = RZEXC_enavg + RZEXC
    if(associated(SRFEXC_enavg) .and. associated(SRFEXC))   & 
        SRFEXC_enavg = SRFEXC_enavg + SRFEXC
    if(associated(GHTCNT1_enavg) .and. associated(GHTCNT1))   & 
        GHTCNT1_enavg = GHTCNT1_enavg + GHTCNT1
    if(associated(GHTCNT2_enavg) .and. associated(GHTCNT2))   & 
        GHTCNT2_enavg = GHTCNT2_enavg + GHTCNT2
    if(associated(GHTCNT3_enavg) .and. associated(GHTCNT3))   & 
        GHTCNT3_enavg = GHTCNT3_enavg + GHTCNT3
    if(associated(GHTCNT4_enavg) .and. associated(GHTCNT4))   & 
        GHTCNT4_enavg = GHTCNT4_enavg + GHTCNT4
    if(associated(GHTCNT5_enavg) .and. associated(GHTCNT5))   & 
        GHTCNT5_enavg = GHTCNT5_enavg + GHTCNT5
    if(associated(GHTCNT6_enavg) .and. associated(GHTCNT6))   & 
        GHTCNT6_enavg = GHTCNT6_enavg + GHTCNT6
    if(associated(WESNN1_enavg) .and. associated(WESNN1))   & 
        WESNN1_enavg = WESNN1_enavg + WESNN1
    if(associated(WESNN2_enavg) .and. associated(WESNN2))   & 
        WESNN2_enavg = WESNN2_enavg + WESNN2
    if(associated(WESNN3_enavg) .and. associated(WESNN3))   & 
        WESNN3_enavg = WESNN3_enavg + WESNN3
    if(associated(HTSNNN1_enavg) .and. associated(HTSNNN1))   & 
        HTSNNN1_enavg = HTSNNN1_enavg + HTSNNN1
    if(associated(HTSNNN2_enavg) .and. associated(HTSNNN2))   & 
        HTSNNN2_enavg = HTSNNN2_enavg + HTSNNN2
    if(associated(HTSNNN3_enavg) .and. associated(HTSNNN3))   & 
        HTSNNN3_enavg = HTSNNN3_enavg + HTSNNN3
    if(associated(SNDZN1_enavg) .and. associated(SNDZN1))   & 
        SNDZN1_enavg = SNDZN1_enavg + SNDZN1
    if(associated(SNDZN2_enavg) .and. associated(SNDZN2))   & 
        SNDZN2_enavg = SNDZN2_enavg + SNDZN2
    if(associated(SNDZN3_enavg) .and. associated(SNDZN3))   & 
        SNDZN3_enavg = SNDZN3_enavg + SNDZN3
    if(associated(EVAPOUT_enavg) .and. associated(EVAPOUT))   & 
        EVAPOUT_enavg = EVAPOUT_enavg + EVAPOUT
    if(associated(SUBLIM_enavg) .and. associated(SUBLIM))   & 
        SUBLIM_enavg = SUBLIM_enavg + SUBLIM
    if(associated(SHOUT_enavg) .and. associated(SHOUT))   & 
        SHOUT_enavg = SHOUT_enavg + SHOUT
    if(associated(RUNOFF_enavg) .and. associated(RUNOFF))   & 
        RUNOFF_enavg = RUNOFF_enavg + RUNOFF
    if(associated(EVPINT_enavg) .and. associated(EVPINT))   & 
        EVPINT_enavg = EVPINT_enavg + EVPINT
    if(associated(EVPSOI_enavg) .and. associated(EVPSOI))   & 
        EVPSOI_enavg = EVPSOI_enavg + EVPSOI
    if(associated(EVPVEG_enavg) .and. associated(EVPVEG))   & 
        EVPVEG_enavg = EVPVEG_enavg + EVPVEG
    if(associated(EVPICE_enavg) .and. associated(EVPICE))   & 
        EVPICE_enavg = EVPICE_enavg + EVPICE
    if(associated(WAT10CM_enavg) .and. associated(WAT10CM))   & 
        WAT10CM_enavg = WAT10CM_enavg + WAT10CM
    if(associated(WATSOI_enavg) .and. associated(WATSOI))   & 
        WATSOI_enavg = WATSOI_enavg + WATSOI
    if(associated(ICESOI_enavg) .and. associated(ICESOI))   & 
        ICESOI_enavg = ICESOI_enavg + ICESOI
    if(associated(EVPSNO_enavg) .and. associated(EVPSNO))   & 
        EVPSNO_enavg = EVPSNO_enavg + EVPSNO
    if(associated(BASEFLOW_enavg) .and. associated(BASEFLOW))   & 
        BASEFLOW_enavg = BASEFLOW_enavg + BASEFLOW
    if(associated(RUNSURF_enavg) .and. associated(RUNSURF))   & 
        RUNSURF_enavg = RUNSURF_enavg + RUNSURF
    if(associated(SMELT_enavg) .and. associated(SMELT))   & 
        SMELT_enavg = SMELT_enavg + SMELT
    if(associated(HLWUP_enavg) .and. associated(HLWUP))   & 
        HLWUP_enavg = HLWUP_enavg + HLWUP
    if(associated(LWNDSRF_enavg) .and. associated(LWNDSRF))   & 
        LWNDSRF_enavg = LWNDSRF_enavg + LWNDSRF
    if(associated(SWNDSRF_enavg) .and. associated(SWNDSRF))   & 
        SWNDSRF_enavg = SWNDSRF_enavg + SWNDSRF
    if(associated(HLATN_enavg) .and. associated(HLATN))   & 
        HLATN_enavg = HLATN_enavg + HLATN
    if(associated(QINFIL_enavg) .and. associated(QINFIL))   & 
        QINFIL_enavg = QINFIL_enavg + QINFIL
    if(associated(GHFLX_enavg) .and. associated(GHFLX))   & 
        GHFLX_enavg = GHFLX_enavg + GHFLX
    if(associated(TPSURF_enavg) .and. associated(TPSURF))   & 
        TPSURF_enavg = TPSURF_enavg + TPSURF
    if(associated(TPSURF_enstd) .and. associated(TPSURF))   &
        TPSURF_enstd = TPSURF_enstd + TPSURF*TPSURF
    if(associated(TPSNOW_enavg) .and. associated(TPSNOW))   & 
        TPSNOW_enavg = TPSNOW_enavg + TPSNOW
    if(associated(TPUNST_enavg) .and. associated(TPUNST))   & 
        TPUNST_enavg = TPUNST_enavg + TPUNST
    if(associated(TPSAT_enavg) .and. associated(TPSAT))   & 
        TPSAT_enavg = TPSAT_enavg + TPSAT
    if(associated(TPWLT_enavg) .and. associated(TPWLT))   & 
        TPWLT_enavg = TPWLT_enavg + TPWLT
    !if(associated(ASNOW_enavg) .and. associated(ASNOW))   & 
    !    ASNOW_enavg = ASNOW_enavg + ASNOW
    if(associated(SHSNOW_enavg) .and. associated(SHSNOW))   & 
        SHSNOW_enavg = SHSNOW_enavg + SHSNOW
    if(associated(AVETSNOW_enavg) .and. associated(AVETSNOW))   & 
        AVETSNOW_enavg = AVETSNOW_enavg + AVETSNOW
    if(associated(FRSAT_enavg) .and. associated(FRSAT))   & 
        FRSAT_enavg = FRSAT_enavg + FRSAT
    if(associated(FRUST_enavg) .and. associated(FRUST))   & 
        FRUST_enavg = FRUST_enavg + FRUST
    if(associated(FRWLT_enavg) .and. associated(FRWLT))   & 
        FRWLT_enavg = FRWLT_enavg + FRWLT
    if(associated(SNOWMASS_enavg) .and. associated(SNOWMASS))   & 
        SNOWMASS_enavg = SNOWMASS_enavg + SNOWMASS
    if(associated(SNOWDP_enavg) .and. associated(SNOWDP))   & 
        SNOWDP_enavg = SNOWDP_enavg + SNOWDP
    if(associated(WET1_enavg) .and. associated(WET1))   & 
        WET1_enavg = WET1_enavg + WET1
    if(associated(WET2_enavg) .and. associated(WET2))   & 
        WET2_enavg = WET2_enavg + WET2
    if(associated(WET3_enavg) .and. associated(WET3))   & 
        WET3_enavg = WET3_enavg + WET3
    if(associated(WCSF_enavg) .and. associated(WCSF))   & 
        WCSF_enavg = WCSF_enavg + WCSF
    if(associated(WCSF_enstd) .and. associated(WCSF))   &
        WCSF_enstd = WCSF_enstd + WCSF*WCSF
    if(associated(WCRZ_enavg) .and. associated(WCRZ))   & 
        WCRZ_enavg = WCRZ_enavg + WCRZ
    if(associated(WCRZ_enstd) .and. associated(WCRZ))   &
        WCRZ_enstd = WCRZ_enstd + WCRZ*WCRZ
    if(associated(WCPR_enavg) .and. associated(WCPR))   & 
        WCPR_enavg = WCPR_enavg + WCPR
    if(associated(WCPR_enstd) .and. associated(WCPR))   &
        WCPR_enstd = WCPR_enstd + WCPR*WCPR
    if(associated(TP1_enavg) .and. associated(TP1))   & 
        TP1_enavg = TP1_enavg + TP1
    if(associated(TP1_enstd) .and. associated(TP1))   &
        TP1_enstd = TP1_enstd + TP1*TP1
    if(associated(TP2_enavg) .and. associated(TP2))   & 
        TP2_enavg = TP2_enavg + TP2
    if(associated(TP3_enavg) .and. associated(TP3))   & 
        TP3_enavg = TP3_enavg + TP3
    if(associated(TP4_enavg) .and. associated(TP4))   & 
        TP4_enavg = TP4_enavg + TP4
    if(associated(TP5_enavg) .and. associated(TP5))   & 
        TP5_enavg = TP5_enavg + TP5
    if(associated(TP6_enavg) .and. associated(TP6))   & 
        TP6_enavg = TP6_enavg + TP6
    if(associated(EMIS_enavg) .and. associated(EMIS))   & 
        EMIS_enavg = EMIS_enavg + EMIS
    if(associated(ALBVR_enavg) .and. associated(ALBVR))   & 
        ALBVR_enavg = ALBVR_enavg + ALBVR
    if(associated(ALBVF_enavg) .and. associated(ALBVF))   & 
        ALBVF_enavg = ALBVF_enavg + ALBVF
    if(associated(ALBNR_enavg) .and. associated(ALBNR))   & 
        ALBNR_enavg = ALBNR_enavg + ALBNR
    if(associated(ALBNF_enavg) .and. associated(ALBNF))   & 
        ALBNF_enavg = ALBNF_enavg + ALBNF
    if(associated(DELTS_enavg) .and. associated(DELTS))   & 
        DELTS_enavg = DELTS_enavg + DELTS
    if(associated(DELQS_enavg) .and. associated(DELQS))   & 
        DELQS_enavg = DELQS_enavg + DELQS
    if(associated(TST_enavg) .and. associated(TST))   & 
        TST_enavg = TST_enavg + TST
    if(associated(LST_enavg) .and. associated(LST))   & 
        LST_enavg = LST_enavg + LST
    if(associated(QST_enavg) .and. associated(QST))   & 
        QST_enavg = QST_enavg + QST
    if(associated(TH_enavg) .and. associated(TH))   & 
        TH_enavg = TH_enavg + TH
    if(associated(QH_enavg) .and. associated(QH))   & 
        QH_enavg = QH_enavg + QH
    if(associated(CHT_enavg) .and. associated(CHT))   & 
        CHT_enavg = CHT_enavg + CHT
    if(associated(CMT_enavg) .and. associated(CMT))   & 
        CMT_enavg = CMT_enavg + CMT
    if(associated(CQT_enavg) .and. associated(CQT))   & 
        CQT_enavg = CQT_enavg + CQT
    if(associated(CNT_enavg) .and. associated(CNT))   & 
        CNT_enavg = CNT_enavg + CNT
    if(associated(RIT_enavg) .and. associated(RIT))   & 
        RIT_enavg = RIT_enavg + RIT
    if(associated(Z0_enavg) .and. associated(Z0))   & 
        Z0_enavg = Z0_enavg + Z0
    if(associated(MOT2M_enavg) .and. associated(MOT2M))   & 
        MOT2M_enavg = MOT2M_enavg + MOT2M
    if(associated(MOQ2M_enavg) .and. associated(MOQ2M))   & 
        MOQ2M_enavg = MOQ2M_enavg + MOQ2M
    if(associated(MOU2M_enavg) .and. associated(MOU2M))   & 
        MOU2M_enavg = MOU2M_enavg + MOU2M
    if(associated(MOV2M_enavg) .and. associated(MOV2M))   & 
        MOV2M_enavg = MOV2M_enavg + MOV2M
    if(associated(MOT10M_enavg) .and. associated(MOT10M))   & 
        MOT10M_enavg = MOT10M_enavg + MOT10M
    if(associated(MOQ10M_enavg) .and. associated(MOQ10M))   & 
        MOQ10M_enavg = MOQ10M_enavg + MOQ10M
    if(associated(MOU10M_enavg) .and. associated(MOU10M))   & 
        MOU10M_enavg = MOU10M_enavg + MOU10M
    if(associated(MOV10M_enavg) .and. associated(MOV10M))   & 
        MOV10M_enavg = MOV10M_enavg + MOV10M
    if(associated(MOU50M_enavg) .and. associated(MOU50M))   & 
        MOU50M_enavg = MOU50M_enavg + MOU50M
    if(associated(MOV50M_enavg) .and. associated(MOV50M))   & 
        MOV50M_enavg = MOV50M_enavg + MOV50M
    if(associated(Z0H_enavg) .and. associated(Z0H))   & 
        Z0H_enavg = Z0H_enavg + Z0H
    if(associated(D0_enavg) .and. associated(D0))   & 
        D0_enavg = D0_enavg + D0
    if(associated(GUST_enavg) .and. associated(GUST))   & 
        GUST_enavg = GUST_enavg + GUST
    if(associated(VENT_enavg) .and. associated(VENT))   & 
        VENT_enavg = VENT_enavg + VENT
    if(associated(ACCUM_enavg) .and. associated(ACCUM))   & 
        ACCUM_enavg = ACCUM_enavg + ACCUM
    if(associated(EVLAND_enavg) .and. associated(EVLAND))   & 
        EVLAND_enavg = EVLAND_enavg + EVLAND
    if(associated(PRLAND_enavg) .and. associated(PRLAND))   & 
        PRLAND_enavg = PRLAND_enavg + PRLAND
    if(associated(SNOLAND_enavg) .and. associated(SNOLAND))   & 
        SNOLAND_enavg = SNOLAND_enavg + SNOLAND
    if(associated(DRPARLAND_enavg) .and. associated(DRPARLAND))   & 
        DRPARLAND_enavg = DRPARLAND_enavg + DRPARLAND
    if(associated(DFPARLAND_enavg) .and. associated(DFPARLAND))   & 
        DFPARLAND_enavg = DFPARLAND_enavg + DFPARLAND
    if(associated(LHSNOW_enavg) .and. associated(LHSNOW))   & 
        LHSNOW_enavg = LHSNOW_enavg + LHSNOW
    if(associated(SWNETSNOW_enavg) .and. associated(SWNETSNOW))   & 
        SWNETSNOW_enavg = SWNETSNOW_enavg + SWNETSNOW
    if(associated(LWUPSNOW_enavg) .and. associated(LWUPSNOW))   & 
        LWUPSNOW_enavg = LWUPSNOW_enavg + LWUPSNOW
    if(associated(LWDNSNOW_enavg) .and. associated(LWDNSNOW))   & 
        LWDNSNOW_enavg = LWDNSNOW_enavg + LWDNSNOW
    if(associated(TCSORIG_enavg) .and. associated(TCSORIG))   & 
        TCSORIG_enavg = TCSORIG_enavg + TCSORIG
    if(associated(TPSN1IN_enavg) .and. associated(TPSN1IN))   & 
        TPSN1IN_enavg = TPSN1IN_enavg + TPSN1IN
    if(associated(TPSN1OUT_enavg) .and. associated(TPSN1OUT))   & 
        TPSN1OUT_enavg = TPSN1OUT_enavg + TPSN1OUT
    if(associated(GHSNOW_enavg) .and. associated(GHSNOW))   & 
        GHSNOW_enavg = GHSNOW_enavg + GHSNOW
    if(associated(LHLAND_enavg) .and. associated(LHLAND))   & 
        LHLAND_enavg = LHLAND_enavg + LHLAND
    if(associated(SHLAND_enavg) .and. associated(SHLAND))   & 
        SHLAND_enavg = SHLAND_enavg + SHLAND
    if(associated(SWLAND_enavg) .and. associated(SWLAND))   & 
        SWLAND_enavg = SWLAND_enavg + SWLAND
    if(associated(SWDOWNLAND_enavg) .and. associated(SWDOWNLAND))   & 
        SWDOWNLAND_enavg = SWDOWNLAND_enavg + SWDOWNLAND
    if(associated(LWLAND_enavg) .and. associated(LWLAND))   & 
        LWLAND_enavg = LWLAND_enavg + LWLAND
    if(associated(GHLAND_enavg) .and. associated(GHLAND))   & 
        GHLAND_enavg = GHLAND_enavg + GHLAND
    if(associated(GHTSKIN_enavg) .and. associated(GHTSKIN))   & 
        GHTSKIN_enavg = GHTSKIN_enavg + GHTSKIN
    if(associated(SMLAND_enavg) .and. associated(SMLAND))   & 
        SMLAND_enavg = SMLAND_enavg + SMLAND
    if(associated(TWLAND_enavg) .and. associated(TWLAND))   & 
        TWLAND_enavg = TWLAND_enavg + TWLAND
    if(associated(TELAND_enavg) .and. associated(TELAND))   & 
        TELAND_enavg = TELAND_enavg + TELAND
    if(associated(TSLAND_enavg) .and. associated(TSLAND))   & 
        TSLAND_enavg = TSLAND_enavg + TSLAND
    if(associated(DWLAND_enavg) .and. associated(DWLAND))   & 
        DWLAND_enavg = DWLAND_enavg + DWLAND
    if(associated(DHLAND_enavg) .and. associated(DHLAND))   & 
        DHLAND_enavg = DHLAND_enavg + DHLAND
    if(associated(SPLAND_enavg) .and. associated(SPLAND))   & 
        SPLAND_enavg = SPLAND_enavg + SPLAND
    if(associated(SPWATR_enavg) .and. associated(SPWATR))   & 
        SPWATR_enavg = SPWATR_enavg + SPWATR
    if(associated(SPSNOW_enavg) .and. associated(SPSNOW))   & 
        SPSNOW_enavg = SPSNOW_enavg + SPSNOW
    if(associated(PEATCLSM_WATERLEVEL_enavg) .and. associated(PEATCLSM_WATERLEVEL))   & 
        PEATCLSM_WATERLEVEL_enavg = PEATCLSM_WATERLEVEL_enavg + PEATCLSM_WATERLEVEL
    if(associated(PEATCLSM_FSWCHANGE_enavg)  .and. associated(PEATCLSM_FSWCHANGE))    & 
        PEATCLSM_FSWCHANGE_enavg  = PEATCLSM_FSWCHANGE_enavg  + PEATCLSM_FSWCHANGE

    if(associated(  CNLAI_enavg) .and. associated( CNLAI))  CNLAI_enavg =  CNLAI_enavg +  CNLAI
    if(associated( CNTLAI_enavg) .and. associated(CNTLAI)) CNTLAI_enavg = CNTLAI_enavg + CNTLAI
    if(associated(  CNSAI_enavg) .and. associated( CNSAI))  CNSAI_enavg =  CNSAI_enavg +  CNSAI
    if(associated( CNTOTC_enavg) .and. associated(CNTOTC)) CNTOTC_enavg = CNTOTC_enavg + CNTOTC
    if(associated( CNVEGC_enavg) .and. associated(CNVEGC)) CNVEGC_enavg = CNVEGC_enavg + CNVEGC
    if(associated( CNROOT_enavg) .and. associated(CNROOT)) CNROOT_enavg = CNROOT_enavg + CNROOT
    if(associated( CNFROOTC_enavg) .and. associated(CNFROOTC)) CNFROOTC_enavg = CNFROOTC_enavg + CNFROOTC
    if(associated(  CNNPP_enavg) .and. associated( CNNPP))  CNNPP_enavg =  CNNPP_enavg +  CNNPP
    if(associated(  CNGPP_enavg) .and. associated( CNGPP))  CNGPP_enavg =  CNGPP_enavg +  CNGPP
    if(associated(   CNSR_enavg) .and. associated(  CNSR))   CNSR_enavg =   CNSR_enavg +   CNSR
    if(associated(  CNNEE_enavg) .and. associated( CNNEE))  CNNEE_enavg =  CNNEE_enavg +  CNNEE
    if(associated(  CNXSMR_enavg).and. associated( CNXSMR))CNXSMR_enavg =  CNXSMR_enavg+ CNXSMR
    if(associated(  CNADD_enavg) .and. associated( CNADD))  CNADD_enavg =  CNADD_enavg +  CNADD
    if(associated( PARABS_enavg) .and. associated(PARABS)) PARABS_enavg = PARABS_enavg + PARABS
    if(associated( PARINC_enavg) .and. associated(PARINC)) PARINC_enavg = PARINC_enavg + PARINC
    if(associated(  SCSAT_enavg) .and. associated( SCSAT))  SCSAT_enavg =  SCSAT_enavg +  SCSAT
    if(associated(  SCUNS_enavg) .and. associated( SCUNS))  SCUNS_enavg =  SCUNS_enavg +  SCUNS
    if(associated( BTRANT_enavg) .and. associated(BTRANT)) BTRANT_enavg = BTRANT_enavg + BTRANT
    if(associated(    SIF_enavg) .and. associated(   SIF))    SIF_enavg =    SIF_enavg +    SIF
    if(associated( CNLOSS_enavg) .and. associated(CNLOSS)) CNLOSS_enavg = CNLOSS_enavg + CNLOSS
    if(associated( CNBURN_enavg) .and. associated(CNBURN)) CNBURN_enavg = CNBURN_enavg + CNBURN
    if(associated( CNFSEL_enavg) .and. associated(CNFSEL)) CNFSEL_enavg = CNFSEL_enavg + CNFSEL

    ! This counter is relative to ens_id
    collect_land_counter = collect_land_counter + 1
    !collect catch_progn

    catch_progn(:,collect_land_counter)%tc1 = TC(:,1)
    catch_progn(:,collect_land_counter)%tc2 = TC(:,2)
    catch_progn(:,collect_land_counter)%tc4 = TC(:,3)

    catch_progn(:,collect_land_counter)%qa1 = QC(:,1)
    catch_progn(:,collect_land_counter)%qa2 = QC(:,2)
    catch_progn(:,collect_land_counter)%qa4 = QC(:,3)

    catch_progn(:,collect_land_counter)%capac  = CAPAC(:)
    catch_progn(:,collect_land_counter)%catdef = catdef(:)
    catch_progn(:,collect_land_counter)%rzexc  = rzexc(:)
    catch_progn(:,collect_land_counter)%srfexc = srfexc(:)

    catch_progn(:,collect_land_counter)%ght(1) = GHTCNT1(:)
    catch_progn(:,collect_land_counter)%ght(2) = GHTCNT2(:)
    catch_progn(:,collect_land_counter)%ght(3) = GHTCNT3(:)
    catch_progn(:,collect_land_counter)%ght(4) = GHTCNT4(:)
    catch_progn(:,collect_land_counter)%ght(5) = GHTCNT5(:)
    catch_progn(:,collect_land_counter)%ght(6) = GHTCNT6(:)

    catch_progn(:,collect_land_counter)%wesn(1) = WESNN1(:)
    catch_progn(:,collect_land_counter)%wesn(2) = WESNN2(:)
    catch_progn(:,collect_land_counter)%wesn(3) = WESNN3(:)

    catch_progn(:,collect_land_counter)%htsn(1) = HTSNNN1(:)
    catch_progn(:,collect_land_counter)%htsn(2) = HTSNNN2(:)
    catch_progn(:,collect_land_counter)%htsn(3) = HTSNNN3(:)

    catch_progn(:,collect_land_counter)%sndz(1) = SNDZN1(:)
    catch_progn(:,collect_land_counter)%sndz(2) = SNDZN2(:)
    catch_progn(:,collect_land_counter)%sndz(3) = SNDZN3(:)


    if(collect_land_counter == NUM_ENSEMBLE) then

        Nm1 = real(NUM_ENSEMBLE-1) 
        if (NUM_ENSEMBLE>1) NdivNm1 = real(NUM_ENSEMBLE)/Nm1

        collect_land_counter = 0
        if(associated(TC_enavg)) TC_enavg = TC_enavg/NUM_ENSEMBLE
        if(associated(QC_enavg)) QC_enavg = QC_enavg/NUM_ENSEMBLE
        if(associated(CAPAC_enavg)) CAPAC_enavg = CAPAC_enavg/NUM_ENSEMBLE
        if(associated(CATDEF_enavg)) CATDEF_enavg = CATDEF_enavg/NUM_ENSEMBLE
        if(associated(RZEXC_enavg)) RZEXC_enavg = RZEXC_enavg/NUM_ENSEMBLE
        if(associated(SRFEXC_enavg)) SRFEXC_enavg = SRFEXC_enavg/NUM_ENSEMBLE
        if(associated(GHTCNT1_enavg)) GHTCNT1_enavg = GHTCNT1_enavg/NUM_ENSEMBLE
        if(associated(GHTCNT2_enavg)) GHTCNT2_enavg = GHTCNT2_enavg/NUM_ENSEMBLE
        if(associated(GHTCNT3_enavg)) GHTCNT3_enavg = GHTCNT3_enavg/NUM_ENSEMBLE
        if(associated(GHTCNT4_enavg)) GHTCNT4_enavg = GHTCNT4_enavg/NUM_ENSEMBLE
        if(associated(GHTCNT5_enavg)) GHTCNT5_enavg = GHTCNT5_enavg/NUM_ENSEMBLE
        if(associated(GHTCNT6_enavg)) GHTCNT6_enavg = GHTCNT6_enavg/NUM_ENSEMBLE
        if(associated(WESNN1_enavg)) WESNN1_enavg = WESNN1_enavg/NUM_ENSEMBLE
        if(associated(WESNN2_enavg)) WESNN2_enavg = WESNN2_enavg/NUM_ENSEMBLE
        if(associated(WESNN3_enavg)) WESNN3_enavg = WESNN3_enavg/NUM_ENSEMBLE
        if(associated(HTSNNN1_enavg)) HTSNNN1_enavg = HTSNNN1_enavg/NUM_ENSEMBLE
        if(associated(HTSNNN2_enavg)) HTSNNN2_enavg = HTSNNN2_enavg/NUM_ENSEMBLE
        if(associated(HTSNNN3_enavg)) HTSNNN3_enavg = HTSNNN3_enavg/NUM_ENSEMBLE
        if(associated(SNDZN1_enavg)) SNDZN1_enavg = SNDZN1_enavg/NUM_ENSEMBLE
        if(associated(SNDZN2_enavg)) SNDZN2_enavg = SNDZN2_enavg/NUM_ENSEMBLE
        if(associated(SNDZN3_enavg)) SNDZN3_enavg = SNDZN3_enavg/NUM_ENSEMBLE
        if(associated(EVAPOUT_enavg)) EVAPOUT_enavg = EVAPOUT_enavg/NUM_ENSEMBLE
        if(associated(SUBLIM_enavg)) SUBLIM_enavg = SUBLIM_enavg/NUM_ENSEMBLE
        if(associated(SHOUT_enavg)) SHOUT_enavg = SHOUT_enavg/NUM_ENSEMBLE
        if(associated(RUNOFF_enavg)) RUNOFF_enavg = RUNOFF_enavg/NUM_ENSEMBLE
        if(associated(EVPINT_enavg)) EVPINT_enavg = EVPINT_enavg/NUM_ENSEMBLE
        if(associated(EVPSOI_enavg)) EVPSOI_enavg = EVPSOI_enavg/NUM_ENSEMBLE
        if(associated(EVPVEG_enavg)) EVPVEG_enavg = EVPVEG_enavg/NUM_ENSEMBLE
        if(associated(EVPICE_enavg)) EVPICE_enavg = EVPICE_enavg/NUM_ENSEMBLE
        if(associated(WAT10CM_enavg)) WAT10CM_enavg = WAT10CM_enavg/NUM_ENSEMBLE
        if(associated(WATSOI_enavg)) WATSOI_enavg = WATSOI_enavg/NUM_ENSEMBLE
        if(associated(ICESOI_enavg)) ICESOI_enavg = ICESOI_enavg/NUM_ENSEMBLE
        if(associated(EVPSNO_enavg)) EVPSNO_enavg = EVPSNO_enavg/NUM_ENSEMBLE
        if(associated(BASEFLOW_enavg)) BASEFLOW_enavg = BASEFLOW_enavg/NUM_ENSEMBLE
        if(associated(RUNSURF_enavg)) RUNSURF_enavg = RUNSURF_enavg/NUM_ENSEMBLE
        if(associated(SMELT_enavg)) SMELT_enavg = SMELT_enavg/NUM_ENSEMBLE
        if(associated(HLWUP_enavg)) HLWUP_enavg = HLWUP_enavg/NUM_ENSEMBLE
        if(associated(LWNDSRF_enavg)) LWNDSRF_enavg = LWNDSRF_enavg/NUM_ENSEMBLE
        if(associated(SWNDSRF_enavg)) SWNDSRF_enavg = SWNDSRF_enavg/NUM_ENSEMBLE
        if(associated(HLATN_enavg)) HLATN_enavg = HLATN_enavg/NUM_ENSEMBLE
        if(associated(QINFIL_enavg)) QINFIL_enavg = QINFIL_enavg/NUM_ENSEMBLE
        !if(associated(AR1_enavg)) AR1_enavg = AR1_enavg/NUM_ENSEMBLE
        !if(associated(AR2_enavg)) AR2_enavg = AR2_enavg/NUM_ENSEMBLE
        !if(associated(RZEQ_enavg)) RZEQ_enavg = RZEQ_enavg/NUM_ENSEMBLE
        if(associated(GHFLX_enavg)) GHFLX_enavg = GHFLX_enavg/NUM_ENSEMBLE
        if(associated(TPSURF_enavg)) TPSURF_enavg = TPSURF_enavg/NUM_ENSEMBLE
        if((NUM_ENSEMBLE>1) .and. associated(TPSURF_enstd) .and. associated(TPSURF_enavg)) then
           TPSURF_enstd = sqrt( TPSURF_enstd/Nm1 - NdivNm1*(TPSURF_enavg**2) )
        else if (associated(TPSURF_enstd)) then
           TPSURF_enstd = MAPL_UNDEF
        end if
        if(associated(TPSNOW_enavg)) TPSNOW_enavg = TPSNOW_enavg/NUM_ENSEMBLE
        if(associated(TPUNST_enavg)) TPUNST_enavg = TPUNST_enavg/NUM_ENSEMBLE
        if(associated(TPSAT_enavg)) TPSAT_enavg = TPSAT_enavg/NUM_ENSEMBLE
        if(associated(TPWLT_enavg)) TPWLT_enavg = TPWLT_enavg/NUM_ENSEMBLE
        if(associated(SHSNOW_enavg)) SHSNOW_enavg = SHSNOW_enavg/NUM_ENSEMBLE
        if(associated(AVETSNOW_enavg)) AVETSNOW_enavg = AVETSNOW_enavg/NUM_ENSEMBLE
        if(associated(FRSAT_enavg)) FRSAT_enavg = FRSAT_enavg/NUM_ENSEMBLE
        if(associated(FRUST_enavg)) FRUST_enavg = FRUST_enavg/NUM_ENSEMBLE
        if(associated(FRWLT_enavg)) FRWLT_enavg = FRWLT_enavg/NUM_ENSEMBLE
        if(associated(ASNOW_enavg)) ASNOW_enavg = max(min(1.0-(FRSAT_enavg+FRUST_enavg+FRWLT_enavg),1.0),0.0)
        if(associated(SNOWMASS_enavg)) SNOWMASS_enavg = SNOWMASS_enavg/NUM_ENSEMBLE
        if(associated(SNOWDP_enavg)) SNOWDP_enavg = SNOWDP_enavg/NUM_ENSEMBLE
        if(associated(WET1_enavg)) WET1_enavg = WET1_enavg/NUM_ENSEMBLE
        if(associated(WET2_enavg)) WET2_enavg = WET2_enavg/NUM_ENSEMBLE
        if(associated(WET3_enavg)) WET3_enavg = WET3_enavg/NUM_ENSEMBLE
        if(associated(WCSF_enavg)) WCSF_enavg = WCSF_enavg/NUM_ENSEMBLE
        if((NUM_ENSEMBLE>1) .and. associated(WCSF_enstd) .and. associated(WCSF_enavg)) then
           WCSF_enstd = sqrt( WCSF_enstd/Nm1 - NdivNm1*(WCSF_enavg**2) )
        else if (associated(WCSF_enstd)) then
           WCSF_enstd = MAPL_UNDEF
        end if
        if(associated(WCRZ_enavg)) WCRZ_enavg = WCRZ_enavg/NUM_ENSEMBLE
        if((NUM_ENSEMBLE>1) .and. associated(WCRZ_enstd) .and. associated(WCRZ_enavg)) then
           WCRZ_enstd = sqrt( WCRZ_enstd/Nm1 - NdivNm1*(WCRZ_enavg**2) )
        else if (associated(WCRZ_enstd)) then
           WCRZ_enstd = MAPL_UNDEF
        end if
        if(associated(WCPR_enavg)) WCPR_enavg = WCPR_enavg/NUM_ENSEMBLE
        if((NUM_ENSEMBLE>1) .and. associated(WCPR_enstd) .and. associated(WCPR_enavg)) then
           WCPR_enstd = sqrt( WCPR_enstd/Nm1 - NdivNm1*(WCPR_enavg**2) ) 
        else if (associated(WCPR_enstd)) then
           WCPR_enstd = MAPL_UNDEF
        end if
        if(associated(TP1_enavg)) TP1_enavg = TP1_enavg/NUM_ENSEMBLE                  ! units K 
        if((NUM_ENSEMBLE>1) .and. associated(TP1_enstd) .and. associated(TP1_enavg)) then
           TP1_enstd = sqrt( TP1_enstd/Nm1 - NdivNm1*(TP1_enavg**2) ) 
        else if (associated(TP1_enstd)) then
           TP1_enstd = MAPL_UNDEF
        end if
        if(associated(TP2_enavg)) TP2_enavg = TP2_enavg/NUM_ENSEMBLE                  ! units now K, rreichle & borescan, 6 Nov 2020
        if(associated(TP3_enavg)) TP3_enavg = TP3_enavg/NUM_ENSEMBLE                  ! units now K, rreichle & borescan, 6 Nov 2020
        if(associated(TP4_enavg)) TP4_enavg = TP4_enavg/NUM_ENSEMBLE                  ! units now K, rreichle & borescan, 6 Nov 2020
        if(associated(TP5_enavg)) TP5_enavg = TP5_enavg/NUM_ENSEMBLE                  ! units now K, rreichle & borescan, 6 Nov 2020
        if(associated(TP6_enavg)) TP6_enavg = TP6_enavg/NUM_ENSEMBLE                  ! units now K, rreichle & borescan, 6 Nov 2020
        if(associated(EMIS_enavg)) EMIS_enavg = EMIS_enavg/NUM_ENSEMBLE
        if(associated(ALBVR_enavg)) ALBVR_enavg = ALBVR_enavg/NUM_ENSEMBLE
        if(associated(ALBVF_enavg)) ALBVF_enavg = ALBVF_enavg/NUM_ENSEMBLE
        if(associated(ALBNR_enavg)) ALBNR_enavg = ALBNR_enavg/NUM_ENSEMBLE
        if(associated(ALBNF_enavg)) ALBNF_enavg = ALBNF_enavg/NUM_ENSEMBLE
        if(associated(DELTS_enavg)) DELTS_enavg = DELTS_enavg/NUM_ENSEMBLE
        if(associated(DELQS_enavg)) DELQS_enavg = DELQS_enavg/NUM_ENSEMBLE
        !if(associated(DELEVAP_enavg)) DELEVAP_enavg = DELEVAP_enavg/NUM_ENSEMBLE
        !if(associated(DELSH_enavg)) DELSH_enavg = DELSH_enavg/NUM_ENSEMBLE
        if(associated(TST_enavg)) TST_enavg = TST_enavg/NUM_ENSEMBLE
        if(associated(LST_enavg)) LST_enavg = LST_enavg/NUM_ENSEMBLE
        if(associated(QST_enavg)) QST_enavg = QST_enavg/NUM_ENSEMBLE
        if(associated(TH_enavg)) TH_enavg = TH_enavg/NUM_ENSEMBLE
        if(associated(QH_enavg)) QH_enavg = QH_enavg/NUM_ENSEMBLE
        if(associated(CHT_enavg)) CHT_enavg = CHT_enavg/NUM_ENSEMBLE
        if(associated(CMT_enavg)) CMT_enavg = CMT_enavg/NUM_ENSEMBLE
        if(associated(CQT_enavg)) CQT_enavg = CQT_enavg/NUM_ENSEMBLE
        if(associated(CNT_enavg)) CNT_enavg = CNT_enavg/NUM_ENSEMBLE
        if(associated(RIT_enavg)) RIT_enavg = RIT_enavg/NUM_ENSEMBLE
        if(associated(Z0_enavg)) Z0_enavg = Z0_enavg/NUM_ENSEMBLE
        if(associated(MOT2M_enavg)) MOT2M_enavg = MOT2M_enavg/NUM_ENSEMBLE
        if(associated(MOQ2M_enavg)) MOQ2M_enavg = MOQ2M_enavg/NUM_ENSEMBLE
        if(associated(MOU2M_enavg)) MOU2M_enavg = MOU2M_enavg/NUM_ENSEMBLE
        if(associated(MOV2M_enavg)) MOV2M_enavg = MOV2M_enavg/NUM_ENSEMBLE
        if(associated(MOT10M_enavg)) MOT10M_enavg = MOT10M_enavg/NUM_ENSEMBLE
        if(associated(MOQ10M_enavg)) MOQ10M_enavg = MOQ10M_enavg/NUM_ENSEMBLE
        if(associated(MOU10M_enavg)) MOU10M_enavg = MOU10M_enavg/NUM_ENSEMBLE
        if(associated(MOV10M_enavg)) MOV10M_enavg = MOV10M_enavg/NUM_ENSEMBLE
        if(associated(MOU50M_enavg)) MOU50M_enavg = MOU50M_enavg/NUM_ENSEMBLE
        if(associated(MOV50M_enavg)) MOV50M_enavg = MOV50M_enavg/NUM_ENSEMBLE
        if(associated(Z0H_enavg)) Z0H_enavg = Z0H_enavg/NUM_ENSEMBLE
        if(associated(D0_enavg)) D0_enavg = D0_enavg/NUM_ENSEMBLE
        if(associated(GUST_enavg)) GUST_enavg = GUST_enavg/NUM_ENSEMBLE
        if(associated(VENT_enavg)) VENT_enavg = VENT_enavg/NUM_ENSEMBLE
        if(associated(ACCUM_enavg)) ACCUM_enavg = ACCUM_enavg/NUM_ENSEMBLE
        if(associated(EVLAND_enavg)) EVLAND_enavg = EVLAND_enavg/NUM_ENSEMBLE
        if(associated(PRLAND_enavg)) PRLAND_enavg = PRLAND_enavg/NUM_ENSEMBLE
        if(associated(SNOLAND_enavg)) SNOLAND_enavg = SNOLAND_enavg/NUM_ENSEMBLE
        if(associated(DRPARLAND_enavg)) DRPARLAND_enavg = DRPARLAND_enavg/NUM_ENSEMBLE
        if(associated(DFPARLAND_enavg)) DFPARLAND_enavg = DFPARLAND_enavg/NUM_ENSEMBLE
        if(associated(LHSNOW_enavg)) LHSNOW_enavg = LHSNOW_enavg/NUM_ENSEMBLE
        if(associated(SWNETSNOW_enavg)) SWNETSNOW_enavg = SWNETSNOW_enavg/NUM_ENSEMBLE
        if(associated(LWUPSNOW_enavg)) LWUPSNOW_enavg = LWUPSNOW_enavg/NUM_ENSEMBLE
        if(associated(LWDNSNOW_enavg)) LWDNSNOW_enavg = LWDNSNOW_enavg/NUM_ENSEMBLE
        if(associated(TCSORIG_enavg)) TCSORIG_enavg = TCSORIG_enavg/NUM_ENSEMBLE
        if(associated(TPSN1IN_enavg)) TPSN1IN_enavg = TPSN1IN_enavg/NUM_ENSEMBLE
        if(associated(TPSN1OUT_enavg)) TPSN1OUT_enavg = TPSN1OUT_enavg/NUM_ENSEMBLE
        if(associated(GHSNOW_enavg)) GHSNOW_enavg = GHSNOW_enavg/NUM_ENSEMBLE
        if(associated(LHLAND_enavg)) LHLAND_enavg = LHLAND_enavg/NUM_ENSEMBLE
        if(associated(SHLAND_enavg)) SHLAND_enavg = SHLAND_enavg/NUM_ENSEMBLE
        if(associated(SWLAND_enavg)) SWLAND_enavg = SWLAND_enavg/NUM_ENSEMBLE
        if(associated(SWDOWNLAND_enavg)) SWDOWNLAND_enavg = SWDOWNLAND_enavg/NUM_ENSEMBLE
        if(associated(LWLAND_enavg)) LWLAND_enavg = LWLAND_enavg/NUM_ENSEMBLE
        if(associated(GHLAND_enavg)) GHLAND_enavg = GHLAND_enavg/NUM_ENSEMBLE
        if(associated(GHTSKIN_enavg)) GHTSKIN_enavg = GHTSKIN_enavg/NUM_ENSEMBLE
        if(associated(SMLAND_enavg)) SMLAND_enavg = SMLAND_enavg/NUM_ENSEMBLE
        if(associated(TWLAND_enavg)) TWLAND_enavg = TWLAND_enavg/NUM_ENSEMBLE
        if(associated(TELAND_enavg)) TELAND_enavg = TELAND_enavg/NUM_ENSEMBLE
        if(associated(TSLAND_enavg)) TSLAND_enavg = TSLAND_enavg/NUM_ENSEMBLE
        if(associated(DWLAND_enavg)) DWLAND_enavg = DWLAND_enavg/NUM_ENSEMBLE
        if(associated(DHLAND_enavg)) DHLAND_enavg = DHLAND_enavg/NUM_ENSEMBLE
        if(associated(SPLAND_enavg)) SPLAND_enavg = SPLAND_enavg/NUM_ENSEMBLE
        if(associated(SPWATR_enavg)) SPWATR_enavg = SPWATR_enavg/NUM_ENSEMBLE
        if(associated(SPSNOW_enavg)) SPSNOW_enavg = SPSNOW_enavg/NUM_ENSEMBLE
        if(associated(PEATCLSM_WATERLEVEL_enavg)) PEATCLSM_WATERLEVEL_enavg = PEATCLSM_WATERLEVEL_enavg/NUM_ENSEMBLE
        if(associated(PEATCLSM_FSWCHANGE_enavg))  PEATCLSM_FSWCHANGE_enavg  = PEATCLSM_FSWCHANGE_enavg /NUM_ENSEMBLE

        if(associated(   CNLAI_enavg))     CNLAI_enavg =    CNLAI_enavg/NUM_ENSEMBLE 
        if(associated(  CNTLAI_enavg))    CNTLAI_enavg =   CNTLAI_enavg/NUM_ENSEMBLE
        if(associated(   CNSAI_enavg))     CNSAI_enavg =    CNSAI_enavg/NUM_ENSEMBLE
        if(associated(  CNTOTC_enavg))    CNTOTC_enavg =   CNTOTC_enavg/NUM_ENSEMBLE
        if(associated(  CNVEGC_enavg))    CNVEGC_enavg =   CNVEGC_enavg/NUM_ENSEMBLE
        if(associated(  CNROOT_enavg))    CNROOT_enavg =   CNROOT_enavg/NUM_ENSEMBLE
        if(associated(  CNFROOTC_enavg))    CNFROOTC_enavg =   CNFROOTC_enavg/NUM_ENSEMBLE
        if(associated(   CNNPP_enavg))     CNNPP_enavg =    CNNPP_enavg/NUM_ENSEMBLE
        if(associated(   CNGPP_enavg))     CNGPP_enavg =    CNGPP_enavg/NUM_ENSEMBLE
        if(associated(    CNSR_enavg))      CNSR_enavg =     CNSR_enavg/NUM_ENSEMBLE
        if(associated(   CNNEE_enavg))     CNNEE_enavg =    CNNEE_enavg/NUM_ENSEMBLE
        if(associated(   CNXSMR_enavg))   CNXSMR_enavg =   CNXSMR_enavg/NUM_ENSEMBLE
        if(associated(   CNADD_enavg))     CNADD_enavg =    CNADD_enavg/NUM_ENSEMBLE
        if(associated(  PARABS_enavg))    PARABS_enavg =   PARABS_enavg/NUM_ENSEMBLE
        if(associated(  PARINC_enavg))    PARINC_enavg =   PARINC_enavg/NUM_ENSEMBLE
        if(associated(   SCSAT_enavg))     SCSAT_enavg =    SCSAT_enavg/NUM_ENSEMBLE
        if(associated(   SCUNS_enavg))     SCUNS_enavg =    SCUNS_enavg/NUM_ENSEMBLE
        if(associated(  BTRANT_enavg))    BTRANT_enavg =   BTRANT_enavg/NUM_ENSEMBLE
        if(associated(     SIF_enavg))       SIF_enavg =      SIF_enavg/NUM_ENSEMBLE
        if(associated(  CNLOSS_enavg))    CNLOSS_enavg =   CNLOSS_enavg/NUM_ENSEMBLE
        if(associated(  CNBURN_enavg))    CNBURN_enavg =   CNBURN_enavg/NUM_ENSEMBLE
        if(associated(  CNFSEL_enavg))    CNFSEL_enavg =   CNFSEL_enavg/NUM_ENSEMBLE

        ! Deal with no-data-values
        !
        ! Surface temperature components may be nodata in some but not all ensemble members.
        !   (Nodata values are assigned in GEOS_CatchGridComp.F90 when the associated 
        !    area fraction is zero.)
        !
        ! For now, the ensemble average is set to nodata if any member has a nodata value.
        !
        ! Alternatively, only ensemble members with good values could be averaged, or the 
        !   averaging could use the associated area fraction as averaging weights.
        !
        ! The simple detection implemented here relies on MAPL_UNDEF being many orders of 
        !   magnitude larger than any valid values, which works fine for Earth surface 
        !   temperatures as long as MAPL_UNDEF is 1.e15.
        !
        ! - reichle, 29 May 2020
                
        if(associated(TPSNOW_enavg))  where (TPSNOW_enavg   > enavg_nodata_threshold)  TPSNOW_enavg              = MAPL_UNDEF
        if(associated(TPSAT_enavg ))  where (TPSAT_enavg    > enavg_nodata_threshold)  TPSAT_enavg               = MAPL_UNDEF
        if(associated(TPWLT_enavg ))  where (TPWLT_enavg    > enavg_nodata_threshold)  TPWLT_enavg               = MAPL_UNDEF
        if(associated(TPUNST_enavg))  where (TPUNST_enavg   > enavg_nodata_threshold)  TPUNST_enavg              = MAPL_UNDEF

        ! restore exact no-data-values for PEATCLSM diagnostics in mineral tiles

        if(associated(PEATCLSM_WATERLEVEL_enavg)) where (PEATCLSM_WATERLEVEL_enavg > enavg_nodata_threshold)  PEATCLSM_WATERLEVEL_enavg = MAPL_UNDEF
        if(associated(PEATCLSM_FSWCHANGE_enavg))  where (PEATCLSM_FSWCHANGE_enavg  > enavg_nodata_threshold)  PEATCLSM_FSWCHANGE_enavg  = MAPL_UNDEF

     end if     ! collect_land_counter==NUM_ENSEMBLE

    ! Turn timers off
    call MAPL_TimerOff(MAPL, "Collect_land")
    call MAPL_TimerOff(MAPL, "TOTAL")
    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Collect_land_ens


  subroutine GET_CATCH_PARAM( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),intent(inout) :: GC     !Gridded component
    type(ESMF_State),   intent(inout) :: IMPORT !Import state
    type(ESMF_State),   intent(inout) :: EXPORT !Export state
    type(ESMF_Clock),   intent(inout) :: CLOCK  !The clock
    integer,optional,   intent(out  ) :: RC     !Error code:

!EOP
! ErrLog Variables

    character(len=ESMF_MAXSTR) :: IAm
    integer :: STATUS
    character(len=ESMF_MAXSTR) :: COMP_NAME
!

! Locals
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    logical :: firsttime = .true.

    real, pointer :: poros(:) =>null()
    real, pointer :: cond(:) =>null()
    real, pointer :: psis(:) =>null()
    real, pointer :: bee(:) =>null()
    real, pointer :: wpwet(:) =>null()
    real, pointer :: gnu(:) =>null()
    real, pointer :: vgwmax(:) =>null()
    real, pointer :: bf1(:) =>null()
    real, pointer :: bf2(:) =>null()
    real, pointer :: bf3(:) =>null()
    real, pointer :: cdcr1(:) =>null()
    real, pointer :: cdcr2(:) =>null()
    real, pointer :: ars1(:) =>null()
    real, pointer :: ars2(:) =>null()
    real, pointer :: ars3(:) =>null()
    real, pointer :: ara1(:) =>null()
    real, pointer :: ara2(:) =>null()
    real, pointer :: ara3(:) =>null()
    real, pointer :: ara4(:) =>null()
    real, pointer :: arw1(:) =>null()
    real, pointer :: arw2(:) =>null()
    real, pointer :: arw3(:) =>null()
    real, pointer :: arw4(:) =>null()
    real, pointer :: tsa1(:) =>null()
    real, pointer :: tsa2(:) =>null()
    real, pointer :: tsb1(:) =>null()
    real, pointer :: tsb2(:) =>null()
    real, pointer :: atau(:) =>null()
    real, pointer :: btau(:) =>null()
    real, pointer :: ity(:) =>null()
    real, pointer :: z2ch(:) =>null()

    real :: SURFLAY, x
    integer :: i

    if (firsttime) then
       firsttime = .false.
       call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_GetPointer(import, poros, 'POROS', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, cond, 'COND', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, psis, 'PSIS', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, bee, 'BEE', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, wpwet, 'WPWET', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, gnu, 'GNU', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, vgwmax, 'VGWMAX', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, bf1, 'BF1', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, bf2, 'BF2', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, bf3, 'BF3', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, cdcr1, 'CDCR1', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, cdcr2, 'CDCR2', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, ars1, 'ARS1', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, ars2, 'ARS2', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, ars3, 'ARS3', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, ara1, 'ARA1', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, ara2, 'ARA2', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, ara3, 'ARA3', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, ara4, 'ARA4', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, arw1, 'ARW1', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, arw2, 'ARW2', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, arw3, 'ARW3', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, arw4, 'ARW4', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, tsa1, 'TSA1', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, tsa2, 'TSA2', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, tsb1, 'TSB1', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, tsb2, 'TSB2', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, atau, 'ATAU', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, btau, 'BTAU', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, ity, 'ITY', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(import, z2ch, 'Z2CH', rc=status)
       VERIFY_(status)

       catch_param(:)%dzgt(1) = dzgt(1)
       catch_param(:)%dzgt(2) = dzgt(2)
       catch_param(:)%dzgt(3) = dzgt(3)
       catch_param(:)%dzgt(4) = dzgt(4)
       catch_param(:)%dzgt(5) = dzgt(5)
       catch_param(:)%dzgt(6) = dzgt(6)
       catch_param(:)%poros = poros
       catch_param(:)%cond  = cond
       catch_param(:)%psis  = psis
       catch_param(:)%bee   = bee
       catch_param(:)%wpwet = wpwet
       catch_param(:)%gnu   = gnu
       catch_param(:)%vgwmax= vgwmax
       catch_param(:)%bf1   = bf1
       catch_param(:)%bf2   = bf2
       catch_param(:)%bf3   = bf3
       catch_param(:)%cdcr1 = cdcr1
       catch_param(:)%cdcr2 = cdcr2
       catch_param(:)%ars1 = ars1
       catch_param(:)%ars2 = ars2
       catch_param(:)%ars3 = ars3
       catch_param(:)%ara1 = ara1
       catch_param(:)%ara2 = ara2
       catch_param(:)%ara3 = ara3
       catch_param(:)%ara4 = ara4
       catch_param(:)%arw1 = arw1
       catch_param(:)%arw2 = arw2
       catch_param(:)%arw3 = arw3
       catch_param(:)%arw4 = arw4
       catch_param(:)%tsa1 = tsa1
       catch_param(:)%tsa2 = tsa2
       catch_param(:)%tsb1 = tsb1
       catch_param(:)%tsb2 = tsb2
       catch_param(:)%atau = atau
       catch_param(:)%btau = btau
       catch_param(:)%vegcls  = nint(ity)
       catch_param(:)%veghght = z2ch

       call MAPL_GetResource(MAPL, SURFLAY, Label="SURFLAY:", DEFAULT=50.0, rc=status)

       catch_param(:)%dzsf = SURFLAY
       catch_param(:)%dzpr = (cdcr2/(1.-wpwet)) / poros
       catch_param(:)%dzrz = vgwmax/poros

       !assign NaN to other fields
       x = ieee_value(x,ieee_quiet_nan)
       catch_param(:)%soilcls30  = transfer(x,i)
       catch_param(:)%soilcls100 = transfer(x,i)
       catch_param(:)%gravel30   = x
       catch_param(:)%orgC30     = x
       catch_param(:)%orgC       = x
       catch_param(:)%sand30     = x
       catch_param(:)%clay30     = x
       catch_param(:)%sand       = x
       catch_param(:)%clay       = x
       catch_param(:)%wpwet30    = x
       catch_param(:)%poros30    = x
       catch_param(:)%dpth       = x
    endif
    RETURN_(ESMF_SUCCESS)
end subroutine GET_CATCH_PARAM


  subroutine Finalize(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! This Finalize routine cleans up the Ldas GridComp

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! Local variables
    type(MAPL_MetaComp), pointer :: MAPL=>null() ! MAPL obj

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Finalize"

   ! Call Finalize for every child
    call MAPL_GenericFinalize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Finalize

end module GEOS_EnsGridCompMod
