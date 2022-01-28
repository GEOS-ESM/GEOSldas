#include "MAPL_Generic.h"
!BOP
! !MODULE: GEOS_LandPertGridCompMod - Module to compute perturbations
module GEOS_LandPertGridCompMod

  ! !USES

  use ESMF
  use ESMF_CFIOMod, only: ESMF_CFIOStrTemplate
  use MAPL_Mod
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer, c_ptr
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_PertTypes, only: allocate_pert_param
  use LDAS_PertTypes, only: T_LANDPERT_STATE
  use LDAS_PertTypes, only: LANDPERT_WRAP

  use nr_ran2_gasdev, only: NRANDSEED, init_randseed
  use LDAS_ConvertMod, only: esmf2ldas
  use LDAS_ensdrv_Globals, only: nodata_generic, nodata_tol_generic
  use LDAS_DriverTypes, only: met_force_type
  use LDAS_DateTimeMod, only: date_time_type, date_time_print
  use RepairForcingMod, only: repair_forcing
  use LDAS_TileCoordType, only: grid_def_type
  use LDAS_TileCoordType, only: tile_coord_type
  use LDAS_TileCoordType, only: T_TILECOORD_STATE
  use LDAS_TileCoordType, only: TILECOORD_WRAP
  use land_pert_routines, only: get_pert, propagate_pert
  use land_pert_routines, only: get_init_pert_rseed
  use LDAS_PertRoutinesMod, only: apply_pert
  use LDAS_PertRoutinesMod, only: get_force_pert_param
  use LDAS_PertRoutinesMod, only: get_progn_pert_param
  use LDAS_PertRoutinesMod, only: read_ens_prop_inputs
  use LDAS_PertRoutinesMod, only: echo_pert_param
  use LDAS_PertRoutinesMod, only: get_pert_grid
  use LDAS_PertRoutinesMod, only: interpolate_pert_to_timestep
  use LDAS_PertRoutinesMod, only: check_pert_dtstep
  use LDAS_PertRoutinesMod, only: GEOSldas_FORCE_PERT_DTSTEP
  use LDAS_PertRoutinesMod, only: GEOSldas_PROGN_PERT_DTSTEP
  use LDAS_PertRoutinesMod, only: GEOSldas_NUM_ENSEMBLE
  use LDAS_PertRoutinesMod, only: GEOSldas_FIRST_ENS_ID
  use LDAS_TileCoordRoutines, only: grid2tile, tile2grid_simple, tile_mask_grid
  use force_and_cat_progn_pert_types, only: N_FORCE_PERT_MAX
  use force_and_cat_progn_pert_types, only: N_PROGN_PERT_MAX
 
  implicit none

  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: SetServices

  ! !DESCRIPTION:
  !EOP

  integer, parameter :: NUM_SUBTILES = 4
  ! the two global varaibles store the mean of pert. They are shared across ensemble
  real,allocatable :: fpert_enavg(:,:,:)
  real,allocatable :: ppert_enavg(:,:,:)
  logical :: phase2_initialized

  integer, public :: N_force_pert, N_progn_pert
  type(pert_param_type),dimension(:), pointer,public :: progn_pert_param =>null()
  type(pert_param_type),dimension(:), pointer,public :: force_pert_param =>null()

  integer,dimension(:,:),pointer,public :: pert_iseed=>null()
  integer :: lat1, lat2, lon1, lon2
  integer :: FIRST_ENS_ID
  logical :: COLDSTART
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
    character(len=ESMF_MAXSTR) :: GridName
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    ! Local variables
    type(T_LANDPERT_STATE), pointer :: internal
    type(LANDPERT_WRAP) :: wrap
    integer :: ens_id

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
    ! -phase-1 : phase2_initilization
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         Phase2_Initialize,                                                     &
         rc=status                                                              &
         )
    VERIFY_(status)
    ! -phase-2 :generate ntrmdt without adjusting
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         GenerateRaw_ntrmdt,                                                    &
         rc=status                                                              &
         )
    VERIFY_(status)
    ! -phase-3 : force-pert-
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         ApplyForcePert,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)
    ! -phase-4 : progn-pert-
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         ApplyPrognPert,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)
    ! -phase-5 : update rseed from assim
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         Update_pert_rseed,                                                     &
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

    ! Allocate an instance of the internal state and put it in wrapper
    ! Then, save the pointer to the wrapper internal state in the GridComp
    allocate(internal, stat=status)
    VERIFY_(status)
    wrap%ptr => internal
    internal%isCubedSphere = .false.

    call ESMF_UserCompSetInternalState(gc, 'Landpert_state', wrap, status)
    VERIFY_(status)

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    call MAPL_GetResource ( MAPL, internal%NUM_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    GEOSldas_NUM_ENSEMBLE = internal%NUM_ENSEMBLE

    call MAPL_GetResource ( MAPL, GridName, Label="GEOSldas.GRIDNAME:", DEFAULT="EASE", RC=STATUS)
    VERIFY_(STATUS)
    if (index(GridName,"-CF") /=0) internal%isCubedSphere = .true.

    call MAPL_GetResource(MAPL, GEOSldas_FIRST_ENS_ID, 'FIRST_ENS_ID:',DEFAULT=0, rc=status)
    VERIFY_(status)

    FIRST_ENS_ID = GEOSldas_FIRST_ENS_ID
    ens_id = FIRST_ENS_ID
    if ( internal%NUM_ENSEMBLE > 1) then
       !landpertxxxx
       read(comp_name(9:12),*) ens_id
    endif
    internal%ens_id= ens_id

    ! Set the state variable specs
    !IMPORT STATE:
    ! ForcePert

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Tair",                                                   &
         LONG_NAME  = "air_temperature_at_RefH",                                &
         UNITS      = "K",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Qair",                                                   &
         LONG_NAME  = "specific_humidity_at_RefH",                              &
         UNITS      = "kg kg-1",                                                &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Psurf",                                                  &
         LONG_NAME  = "surface_pressure",                                       &
         UNITS      = "Pa",                                                     &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Rainf_C",                                                &
         LONG_NAME  = "convective_rainfall",                                    &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Rainf",                                                  &
         LONG_NAME  = "liquid_water_large_scale_precipitation",                 &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Snowf",                                                  &
         LONG_NAME  = "total_snowfall",                                         &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "LWdown",                                                 &
         LONG_NAME  = "downward_longwave_radiation",                            &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "SWdown",                                                 &
         LONG_NAME  = "downward_shortwave_radiation",                           &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PARdrct",                                                &
         LONG_NAME  = "photosynth_active_radiation_direct",                     &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PARdffs",                                                &
         LONG_NAME  = "photosynth_active_radiation_diffuse",                    &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Wind",                                                   &
         LONG_NAME  = "wind_speed_at_RefH",                                     &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "RefH",                                                   &
         LONG_NAME  = "reference_height_for_Tair_Qair_Wind",                    &
         UNITS      = "m",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    ! PrognPert, the connected to land's exports which are catchment's Internal

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'TCPert',                                                 &
         LONG_NAME  = 'canopy_temperature',                                     &
         UNITS      = 'K',                                                      &
         DIMS       = MAPL_DimsTileTile,                                        &
         NUM_SUBTILES = NUM_SUBTILES,                                           &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'CATDEFPert',                                             &
         LONG_NAME  = 'catchment_deficit',                                      &
         UNITS      = 'kg m-2',                                                 &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'RZEXCPert',                                              &
         LONG_NAME  = 'root_zone_excess',                                       &
         UNITS      = 'kg m-2',                                                 &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'SRFEXCPert',                                             &
         LONG_NAME  = 'surface_excess',                                         &
         UNITS      = 'kg m-2',                                                 &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'GHTCNT1Pert',                                            &
         LONG_NAME  = 'soil_heat_content_layer_1',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'GHTCNT2Pert',                                            &
         LONG_NAME  = 'soil_heat_content_layer_2',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'GHTCNT3Pert',                                            &
         LONG_NAME  = 'soil_heat_content_layer_3',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'GHTCNT4Pert',                                            &
         LONG_NAME  = 'soil_heat_content_layer_4',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'GHTCNT5Pert',                                            &
         LONG_NAME  = 'soil_heat_content_layer_5',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'GHTCNT6Pert',                                            &
         LONG_NAME  = 'soil_heat_content_layer_6',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)


    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'WESNN1Pert',                                             &
         LONG_NAME  = 'snow_mass_layer_1',                                      &
         UNITS      = 'kg m-2',                                                 &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'WESNN2Pert',                                             &
         LONG_NAME  = 'snow_mass_layer_2',                                      &
         UNITS      = 'kg m-2',                                                 &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'WESNN3Pert',                                             &
         LONG_NAME  = 'snow_mass_layer_3',                                      &
         UNITS      = 'kg m-2',                                                 &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'HTSNNN1Pert',                                            &
         LONG_NAME  = 'heat_content_snow_layer_1',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'HTSNNN2Pert',                                            &
         LONG_NAME  = 'heat_content_snow_layer_2',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'HTSNNN3Pert',                                            &
         LONG_NAME  = 'heat_content_snow_layer_3',                              &
         UNITS      = 'J m-2',                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'SNDZN1Pert',                                             &
         LONG_NAME  = 'snow_delth_layer_1',                                     &
         UNITS      = 'm',                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'SNDZN2Pert',                                             &
         LONG_NAME  = 'snow_delth_layer_2',                                     &
         UNITS      = 'm',                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = 'SNDZN3Pert',                                             &
         LONG_NAME  = 'snow_delth_layer_3',                                     &
         UNITS      = 'm',                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VLocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(STATUS)

    ! !EXPORT STATE:

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "TApert",                                                 &
         LONG_NAME  = "perturbed_surface_air_temperature",                      &
         UNITS      = "K",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "QApert",                                                 &
         LONG_NAME  = "perturbed_surface_air_specific_humidity",                &
         UNITS      = "kg kg-1",                                                &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PSpert",                                                 &
         LONG_NAME  = "Perturbed_surface_pressure",                             &
         UNITS      = "Pa",                                                     &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "UUpert",                                                 &
         LONG_NAME  = "perturbed_surface_wind_speed",                           &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "UWINDLMTILEpert",                                        &
         LONG_NAME  = "perturbed_levellm_uwind",                                &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "VWINDLMTILEpert",                                        &
         LONG_NAME  = "perturbed_levellm_vwind",                                &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PCUpert",                                                &
         LONG_NAME  = "perturbed_liquid_water_convective_precipitation",        &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PLSpert",                                                &
         LONG_NAME  = "perturbed_liquid_water_large_scale_precipitation",       &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "SNOpert",                                                &
         LONG_NAME  = "perturbed_snowfall",                                     &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRPARpert",                                              &
         LONG_NAME  = "surface_downwelling_par_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFPARpert",                                              &
         LONG_NAME  = "surface_downwelling_par_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRNIRpert",                                              &
         LONG_NAME  = "surface_downwelling_nir_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFNIRpert",                                              &
         LONG_NAME  = "surface_downwelling_nir_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRUVRpert",                                              &
         LONG_NAME  = "surface_downwelling_uvr_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFUVRpert",                                              &
         LONG_NAME  = "surface_downwelling_uvr_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "LWDNSRFpert",                                            &
         LONG_NAME  = "perturbed_surface_downwelling_longwave_flux",            &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DZpert",                                                   &
         LONG_NAME  = "reference_height_for_Tair_Qair_Wind",                    &
         UNITS      = "m",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    ! !INTERNAL STATE:
    if ( .not. internal%isCubedSphere ) then

       call MAPL_AddInternalSpec(                                                  &
            gc,                                                                    &
            SHORT_NAME = "pert_rseed",                                             &
            LONG_NAME  = "Perturbations_rseed",                                    &
            UNITS      = "1",                                                      &
            PRECISION  = ESMF_KIND_R8,                                             &
            FRIENDLYTO = trim(COMP_NAME),                                          &
            DIMS       = MAPL_DimsNone,                                            &
            UNGRIDDED_DIMS = (/NRANDSEED/),                                        &
            VLOCATION  = MAPL_VlocationNone,                                       &
            rc         = status                                                    &
            )
       VERIFY_(status)

       call MAPL_AddInternalSpec(                                                  &
            gc,                                                                    &
            SHORT_NAME = "fpert_ntrmdt",                                           &
            LONG_NAME  = "force_pert_intermediate",                                &
            UNITS      = "1",                                                      &
            DIMS       = MAPL_DimsHorzOnly,                                        &
            UNGRIDDED_DIMS = (/N_FORCE_PERT_MAX/),                                 &
            VLOCATION  = MAPL_VlocationNone,                                       &
            rc         = status                                                    &
            )
       VERIFY_(status)

       call MAPL_AddInternalSpec(                                                  &
            gc,                                                                    &
            SHORT_NAME = "ppert_ntrmdt",                                           &
            LONG_NAME  = "progn_pert_intermediate",                                &
            UNITS      = "1",                                                      &
            DIMS       = MAPL_DimsHorzOnly,                                        &
            UNGRIDDED_DIMS = (/N_PROGN_PERT_MAX/),                                 &
            VLOCATION  = MAPL_VlocationNone,                                       &
            rc         = status                                                    &
            )
       VERIFY_(status)

    endif

    !EOS

    ! Set profiling timers
    call MAPL_TimerAdd(gc, name="Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="phase2_Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="GenerateRaw", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Run_ApplyForcePert", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="-GetPert", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="-ApplyPert", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="-MetForcing2Catch", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="-LocStreamTransform", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Run_ApplyPrognPert", rc=status)
    VERIFY_(status)

    ! Call SetServices for children
    call MAPL_GenericSetServices(gc, rc=status)
    VERIFY_(status)

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


  !BOP

  ! !IROUTINE: Initialize -- initialize method for LDAS GC

  ! !INTERFACE:

  subroutine Initialize(gc, import, export, clock, rc)

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
    character(len=ESMF_MAXSTR) :: rst_fname, rst_fname_tmp

    ! ESMF variables
    type(ESMF_VM) :: vm
    type(ESMF_TimeInterval) :: ModelTimeStep
    type(ESMF_Time) :: CurrentTime, StopTime
    type(ESMF_Alarm) :: ForcePertAlarm, PrognPertAlarm
    type(ESMF_TimeInterval) :: ForcePert_DT, PrognPert_DT
    type(ESMF_State) :: MINTERNAL

    ! LDAS variables
    type(date_time_type) :: stop_time, current_time

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(MAPL_LocStream) :: locstream
    character(len=300) :: out_path
    character(len=ESMF_MAXSTR) :: exp_id

    ! Internal private state variables
    type(T_LANDPERT_STATE), pointer :: internal=>null()
    type(LANDPERT_WRAP) :: wrap
 
    type(TILECOORD_WRAP) :: tcwrap
    type(tile_coord_type), pointer :: tile_coord(:)=>null()

    ! MAPL internal pointers
    real, pointer :: pert_ptr(:)=>null()
    real, pointer :: fpert_ntrmdt(:,:,:)=>null()
    real, pointer :: ppert_ntrmdt(:,:,:)=>null()
    real(kind=ESMF_KIND_R8), pointer :: pert_rseed_r8(:)=>null()

    ! Misc variables
    integer :: imjm(7), imjm_global(7) ! we need just the first 2
    integer :: model_dtstep
    integer :: land_nt_local,m,n, i1, in, j1, jn
    logical :: IAmRoot, f_exist
    integer :: ipert,n_lon,n_lat, n_lon_g, n_lat_g
    integer, allocatable :: pert_rseed(:)
    real :: dlon, dlat,locallat,locallon
    type(ESMF_Grid) :: Grid
    character(len=ESMF_MAXSTR) :: id_string
    integer :: ens_id_width
    ! Begin...

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Initialize"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    call MAPL_TimerOn(MAPL, "TOTAL")

    call MAPL_TimerOn(MAPL, "Initialize")

    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    ! Get pointer to the private internal state
    call ESMF_UserCompGetInternalState(gc, 'Landpert_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr

   ! Get component's internal tile_coord variable
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tile_coord => tcwrap%ptr%tile_coord

    ! Are we perturbing variables?
    call MAPL_GetResource(MAPL, internal%PERTURBATIONS, 'PERTURBATIONS:', default=0, rc=status)
    VERIFY_(status)

    if (internal%PERTURBATIONS == 0) then ! no perturbations
       allocate(progn_pert_param(0))
       allocate(force_pert_param(0))
       call MAPL_TimerOff(MAPL, "Initialize")
       call MAPL_TimerOff(MAPL, "TOTAL")
       RETURN_(ESMF_SUCCESS)
    end if

    ! MPI
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    IAmRoot = MAPL_Am_I_Root(vm)


    call MAPL_GetResource ( MAPL, out_path, Label="OUT_PATH:", DEFAULT="./", RC=STATUS)
    call MAPL_GetResource ( MAPL, exp_id, Label="EXP_ID:", DEFAULT="NO_ID", RC=STATUS)
   
    call MAPL_GetResource(MAPL, internal%ForcePert%dtstep, 'FORCE_PERT_DTSTEP:',DEFAULT=10800, rc=status)
     VERIFY_(status)
    call ESMF_TimeIntervalSet(ForcePert_DT, s=internal%ForcePert%dtstep, rc=status)
    VERIFY_(status)
    ! -PrognPert-
    call MAPL_GetResource(MAPL, internal%PrognPert%dtstep, 'PROGN_PERT_DTSTEP:',DEFAULT=10800, rc=status)
    VERIFY_(status)
    call ESMF_TimeIntervalSet(PrognPert_DT, s=internal%PrognPert%dtstep, rc=status)
    VERIFY_(status)
     
    

    GEOSldas_FORCE_PERT_DTSTEP = internal%ForcePert%dtstep
    GEOSldas_PROGN_PERT_DTSTEP = internal%PrognPert%dtstep

    call get_pert_grid(tcwrap%ptr%grid_g,internal%pgrid_g)
    call get_pert_grid(tcwrap%ptr%grid_f,internal%pgrid_f)
    call get_pert_grid(tcwrap%ptr%grid_l,internal%pgrid_l)
 
    n_lon = internal%pgrid_l%n_lon
    n_lat = internal%pgrid_l%n_lat

    call MAPL_GetResource( MAPL, ens_id_width,"ENS_ID_WIDTH:", default=4, RC=STATUS)
    VERIFY_(status)

   ! Pointers to mapl internals
    if ( internal%isCubedSphere ) then
       n_lon_g = internal%pgrid_g%n_lon
       n_lat_g = internal%pgrid_g%n_lat
       allocate(internal%fpert_ntrmdt(n_lon_g, n_lat_g, N_FORCE_PERT_MAX), source=0.0)
       allocate(internal%ppert_ntrmdt(n_lon_g, n_lat_g, N_PROGN_PERT_MAX), source=0.0)
       allocate(internal%pert_rseed_r8(NRANDSEED), source=0.0d0)
       
       fpert_ntrmdt => internal%fpert_ntrmdt      
       ppert_ntrmdt => internal%ppert_ntrmdt      
       pert_rseed_r8 => internal%pert_rseed_r8      

       call MAPL_GetResource(MAPL, rst_fname_tmp, 'LANDPERT_INTERNAL_RESTART_FILE:',DEFAULT='NONE', rc=status)
       VERIFY_(status)

       id_string=""
       if (internal%NUM_ENSEMBLE > 1) then
         n = len(trim(COMP_NAME))
         id_string = COMP_NAME(n-ens_id_width+1:n)
       endif

       call ESMF_CFIOStrTemplate(rst_fname, trim(adjustl(rst_fname_tmp)),'GRADS', xid = trim(id_string), stat=status)

       if (index(rst_fname, 'NONE') == 0 ) then
          f_exist = .false.
          if ( IAmRoot) then
            inquire(file=rst_fname, exist=f_exist)
            if (f_exist) call read_pert_rst(trim(rst_fname), fpert_ntrmdt, ppert_ntrmdt, pert_rseed_r8) 
          endif
          call MAPL_CommsBcast(vm, data=f_exist,  N=1, ROOT=0,rc=status) 
          if (f_exist) then
             n = n_lat_g*n_lon_g*N_FORCE_PERT_MAX

             block
               type(c_ptr) :: cptr
               cptr = c_loc(fpert_ntrmdt(1,1,1))
               call c_f_pointer(cptr, pert_ptr, [n])
             end block

             call MAPL_CommsBcast(vm, data=pert_ptr,  N=n, ROOT=0,rc=status)
             VERIFY_(status)
             pert_ptr=>null()

             n = n_lat_g*n_lon_g*N_PROGN_PERT_MAX
             block
               type(c_ptr) :: cptr
               cptr = c_loc(ppert_ntrmdt(1,1,1))
               call c_f_pointer(cptr, pert_ptr, [n])
             end block
             call MAPL_CommsBcast(vm, data=pert_ptr,  N=n, ROOT=0,rc=status)
             VERIFY_(status)
             pert_ptr=>null()

             call MAPL_CommsBcast(vm, data=pert_rseed_r8, N=NRANDSEED, ROOT=0,rc=status)
             VERIFY_(status)
          endif
       endif 
       lon1 = internal%pgrid_l%i_offg + 1
       lon2 = internal%pgrid_l%i_offg + n_lon
       lat1 = internal%pgrid_l%j_offg + 1
       lat2 = internal%pgrid_l%j_offg + n_lat
    else
       call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=MINTERNAL, rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, fpert_ntrmdt, 'fpert_ntrmdt', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, ppert_ntrmdt, 'ppert_ntrmdt', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, pert_rseed_r8, 'pert_rseed', rc=status)
       VERIFY_(status)
      ! Get grid info from the gridcomp
       call ESMF_GridCompGet(gc, grid=Grid, rc=status)
       VERIFY_(status)
       call ESMF_GRID_INTERIOR  (GRID, I1,IN,J1,JN)

       lon1 = internal%pgrid_l%i_offg + 1 ! global index, starting from 1
       lon1 = lon1 - i1 + 1               ! relative to local
       lon2 = lon1 + n_lon - 1
       lat1 = internal%pgrid_l%j_offg + 1 ! global index, starting from 1
       lat1 = lat1 - j1 +1                ! relative to local
       lat2 = lat1 + n_lat - 1
    endif

    ! Convert pert_rseed_r8 to integer
    allocate(pert_rseed(size(pert_rseed_r8)), source=0, stat=status)
    VERIFY_(status)
    pert_rseed = nint(pert_rseed_r8)
    if( .not. associated(pert_iseed)) then
       allocate(pert_iseed(size(pert_rseed_r8),internal%NUM_ENSEMBLE), source=0, stat=status)
    endif
    VERIFY_(status)
    ! Check if we need to coldstart -
    ! If the MAPL internal state variables are zero, a restart file was
    ! not available to read - in that case we cold-start
    COLDSTART = .false.
    if (all(pert_rseed==0)) COLDSTART = .true.

    ! Get number of land tiles
    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)
    call MAPL_LocStreamGet(locstream, NT_LOCAL=land_nt_local,rc=status)
    VERIFY_(status)

   
    allocate(internal%i_indgs(land_nt_local),stat=status)
    VERIFY_(status) 
    allocate(internal%j_indgs(land_nt_local),stat=status)
    VERIFY_(status)
    internal%i_indgs(:)=tile_coord(:)%hash_i_indg
    internal%j_indgs(:)=tile_coord(:)%hash_j_indg

    ! Get pert options from *default* namelist files
    ! WARNING: get_force/progn_pert_param() calls allocate memory

    call get_force_pert_param(internal%pgrid_l, internal%ForcePert%npert, internal%ForcePert%param)
    _ASSERT(internal%ForcePert%npert==size(internal%ForcePert%param), "ForcePert: param size does not match npert")

    internal%ForcePert%fft_npert = internal%ForcePert%npert
    call MAPL_CommsBcast(vm, data=internal%ForcePert%fft_npert, N=1, ROOT=0,rc=status)
    if (size(internal%ForcePert%param) == 0 .and. internal%ForcePert%fft_npert >0 ) then
       allocate(internal%ForcePert%param(internal%ForcePert%fft_npert))
    endif
    do n = 1, internal%ForcePert%fft_npert
       call MAPL_CommsBcast(vm, data=internal%ForcePert%param(n)%xcorr,    N=1, ROOT=0,rc=status)
       call MAPL_CommsBcast(vm, data=internal%ForcePert%param(n)%ycorr,    N=1, ROOT=0,rc=status)
       call MAPL_CommsBcast(vm, data=internal%ForcePert%param(n)%tcorr,    N=1, ROOT=0,rc=status)
       call MAPL_CommsBcast(vm, data=internal%ForcePert%param(n)%coarsen,  N=1, ROOT=0,rc=status)
    enddo

    call get_progn_pert_param(internal%pgrid_l, internal%PrognPert%npert, internal%PrognPert%param)
    _ASSERT(internal%PrognPert%npert==size(internal%PrognPert%param), "PrognPert: param size does not match npert")

    internal%PrognPert%fft_npert = internal%PrognPert%npert
    call MAPL_CommsBcast(vm, data=internal%PrognPert%fft_npert, N=1, ROOT=0,rc=status)
    if (size(internal%PrognPert%param) == 0 .and. internal%PrognPert%fft_npert > 0) then
       allocate(internal%PrognPert%param(internal%PrognPert%fft_npert))
    endif

    do n = 1, internal%PrognPert%fft_npert
       call MAPL_CommsBcast(vm, data=internal%PrognPert%param(n)%xcorr,    N=1, ROOT=0,rc=status)
       call MAPL_CommsBcast(vm, data=internal%PrognPert%param(n)%ycorr,    N=1, ROOT=0,rc=status)
       call MAPL_CommsBcast(vm, data=internal%PrognPert%param(n)%tcorr,    N=1, ROOT=0,rc=status)
       call MAPL_CommsBcast(vm, data=internal%PrognPert%param(n)%coarsen,  N=1, ROOT=0,rc=status)
    enddo

    N_force_pert = internal%ForcePert%npert
    N_progn_pert = internal%PrognPert%npert

    ! params are the same across ensemble
    if (.not. associated (progn_pert_param)) then
        progn_pert_param=>internal%PrognPert%param
    endif

    if (.not. associated (force_pert_param)) then
        force_pert_param=>internal%ForcePert%param
    endif

    !allocate(internal%fpert_ntrmdt(n_lon,n_lat,internal%ForcePert%npert), source=0., stat=status)
    !VERIFY_(status)
    !allocate(internal%ppert_ntrmdt(n_lon,n_lat,internal%PrognPert%npert), source=0., stat=status)
    !VERIFY_(status)

    !fpert_ntrmdt=>internal%fpert_ntrmdt
    !ppert_ntrmdt=>internal%ppert_ntrmdt

    ! allocate the global vaiable
    if( .not. allocated(fpert_enavg) ) then
       allocate(fpert_enavg(n_lon,n_lat,internal%ForcePert%npert), source=0., stat=status)
       VERIFY_(status)
    endif
    if( .not. allocated(ppert_enavg) ) then
       allocate(ppert_enavg(n_lon,n_lat,internal%PrognPert%npert), source=0., stat=status)
       VERIFY_(status)
    endif

    if (IAmRoot .and. internal%ens_id == FIRST_ENS_ID) then
       call echo_pert_param( internal%ForcePert%npert, internal%ForcePert%param, 1, 1 )
       call echo_pert_param( internal%PrognPert%npert, internal%PrognPert%param, 1, 1 )
    endif

    ! Allocate and initialize  memory for pvt internal state variables
    allocate(internal%ForcePert%DataPrv(land_nt_local, internal%ForcePert%npert), stat=status)
    VERIFY_(status)
    allocate(internal%ForcePert%DataNxt(land_nt_local, internal%ForcePert%npert), stat=status)
    VERIFY_(status)
    allocate(internal%PrognPert%DataPrv(land_nt_local, internal%PrognPert%npert), stat=status)
    VERIFY_(status)
    allocate(internal%PrognPert%DataNxt(land_nt_local, internal%PrognPert%npert), stat=status)
    VERIFY_(status)
    internal%ForcePert%DataPrv = MAPL_UNDEF
    internal%ForcePert%DataNxt = MAPL_UNDEF
    internal%PrognPert%DataPrv = MAPL_UNDEF
    internal%PrognPert%DataNxt = MAPL_UNDEF

    ! Coldstart
    if (COLDSTART) then
       if (IAmRoot) print *, trim(Iam)//'::WARNING: Cold-starting '// trim(COMP_NAME) // ' GridComp'
       ! -pert_rseed-
       call get_init_pert_rseed(internal%ens_id, pert_rseed(1))
       call init_randseed(pert_rseed)
       ! -ForcePert-
       call propagate_pert(                                                     &
            internal%ForcePert%fft_npert,                                       &
            1,                                                                  &
            internal%pgrid_l, internal%pgrid_f,                                 &
            ! arbitrary dtstep
            -1.0,                                                               &
            pert_rseed,                                                         &
            internal%ForcePert%param,                                           &
            fpert_ntrmdt(lon1:lon2,lat1:lat2,                                   &
                         1:internal%ForcePert%npert),                           &
            ! initialize
            .true.                                                              &
            )

       ! -prognostics-
       call propagate_pert(                                                     &
            internal%PrognPert%fft_npert,                                       &
            1,                                                                  &
            internal%pgrid_l, internal%pgrid_f,                                 &
            ! arbitrary dtstep
            -1.0,                                                               &
            pert_rseed,                                                         &
            internal%PrognPert%param,                                           &
            ppert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%PrognPert%npert),       &
            ! initialize
            .true.                                                              &
            )

   
    end if

    if(internal%ens_id == FIRST_ENS_ID ) fpert_enavg(:,:,:)=0. 

    do m = 1,internal%ForcePert%npert
       call tile_mask_grid(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), fpert_ntrmdt(lon1:lon2,lat1:lat2,m))
       if(internal%ForcePert%param(m)%zeromean .and. internal%NUM_ENSEMBLE >2) then
          fpert_enavg(:,:,m)=fpert_enavg(:,:,m)+fpert_ntrmdt(lon1:lon2,lat1:lat2,m)
          if( internal%ens_id-FIRST_ENS_ID == internal%NUM_ENSEMBLE-1) then
             fpert_enavg(:,:,m) = -fpert_enavg(:,:,m)/real(internal%NUM_ENSEMBLE)
          endif
       endif
    enddo

    if(internal%ens_id == FIRST_ENS_ID) ppert_enavg(:,:,:)=0. 

    do m = 1,internal%PrognPert%npert  
       call tile_mask_grid(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), ppert_ntrmdt(lon1:lon2,lat1:lat2,m))
       if(internal%PrognPert%param(m)%zeromean .and. internal%NUM_ENSEMBLE >2) then
          ppert_enavg(:,:,m)=ppert_enavg(:,:,m)+ppert_ntrmdt(lon1:lon2,lat1:lat2,m)
          if( internal%ens_id - FIRST_ENS_ID == internal%NUM_ENSEMBLE-1) then
             ppert_enavg(:,:,m) = -ppert_enavg(:,:,m)/real(internal%NUM_ENSEMBLE)
          endif
       endif
    enddo

    ! Check force/progn pert dtsteps against model dtstep
    ! -Get-model-times-
    call ESMF_ClockGet(                                                         &
         clock,                                                                 &
         currTime=CurrentTime,                                                  &
         timeStep=ModelTimeStep,                                                &
         stopTime=StopTime                                                      &
         )
    VERIFY_(status)
    ! -model-dtstep-in-seconds-
    call ESMF_TimeIntervalGet(ModelTimeStep, s=model_dtstep)
    VERIFY_(status)
    ! -model-times-in-LDAS-datetime-format-
    call esmf2ldas(StopTime, stop_time, rc=status)
    VERIFY_(status)
    call esmf2ldas(CurrentTime, current_time, rc=status)
    VERIFY_(status)

    if( internal%ens_id == FIRST_ENS_ID .and. IAmRoot) then
       ! write out the input file
       call read_ens_prop_inputs(write_nml = .true. , work_path = trim(out_path), &
            exp_id = trim(exp_id), date_time = current_time)
    endif


    ! -Now-check-pert-dtstep-
    call check_pert_dtstep(                                                     &
         model_dtstep,                                                          &
         current_time, stop_time,                                               &
         internal%PrognPert%npert, internal%ForcePert%npert,                    &
         internal%PrognPert%dtstep, internal%ForcePert%dtstep                   &
         )

    ! Create (non-sticky) alarms for force and progn perturbations
    ! -ForcePert-
    ForcePertAlarm = ESMF_AlarmCreate(                                          &
         clock,                                                                 &
         name='ForcePert',                                                      &
         ringTime=CurrentTime,                                                  &
         ringInterval=ForcePert_DT,                                             &
         ringTimeStepCount=1,                                                   &
         sticky=.false.,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)
    ! -PrognPert-
    PrognPertAlarm = ESMF_AlarmCreate(                                          &
         clock,                                                                 &
         name='PrognPert',                                                      &
         ringTime=CurrentTime,                                                  &
         ringInterval=PrognPert_DT,                                             &
         ringTimeStepCount=1,                                                   &
         sticky=.false.,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)

    ! Perturbation times
    ! -force-
    internal%ForcePert%TimePrv = CurrentTime
    internal%ForcePert%TimeNxt = CurrentTime
    ! -progn-
    internal%PrognPert%TimePrv = CurrentTime
    internal%PrognPert%TimeNxt = CurrentTime

    ! Update the r4 version of pert_rseed
    pert_rseed_r8 = real(pert_rseed,kind=ESMF_KIND_R8)
    pert_iseed(:,internal%ens_id + 1 - FIRST_ENS_ID ) = pert_rseed
    ! Clean up

    if (allocated(pert_rseed)) then ! integer version of MINTERNAL state
       deallocate(pert_rseed, stat=status)
       VERIFY_(status)
    end if

    Phase2_initialized = .false.

    ! Turn timers off
    call MAPL_TimerOff(MAPL, "Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


  subroutine Phase2_Initialize(gc, import, export, clock, rc)

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
    type(ESMF_State) :: MINTERNAL

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(MAPL_LocStream) :: locstream

    ! Internal private state variables
    type(T_LANDPERT_STATE), pointer :: internal=>null()
    type(LANDPERT_WRAP) :: wrap
 
    type(TILECOORD_WRAP) :: tcwrap
    type(tile_coord_type), pointer :: tile_coord(:)=>null()

    ! MAPL internal pointers
    real, pointer :: fpert_ntrmdt(:,:,:)=>null()
    real, pointer :: ppert_ntrmdt(:,:,:)=>null()
    real(kind=ESMF_KIND_R8), pointer :: pert_rseed_r8(:)=>null()

    ! Misc variables
    real, allocatable :: fpert_grid(:,:,:), ppert_grid(:,:,:)
    integer,allocatable :: pert_rseed(:)

    integer :: land_nt_local,n,ipert,i,j,n_lon,n_lat
    logical :: IAmRoot
    real :: locallat,locallon

    ! Begin...
    ! phase2_initialized is a global variables shared by all ensemble member
    if( phase2_initialized) then
       RETURN_(ESMF_SUCCESS)
    endif

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::phase2_Initialize"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    call MAPL_TimerOn(MAPL, "TOTAL")

    call MAPL_TimerOn(MAPL, "phase2_Initialize")

    ! Get pointer to the private internal state
    call ESMF_UserCompGetInternalState(gc, 'Landpert_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr
    n_lon = internal%pgrid_l%n_lon
    n_lat = internal%pgrid_l%n_lat

   ! Get component's internal tile_coord variable
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tile_coord => tcwrap%ptr%tile_coord

    if (internal%PERTURBATIONS == 0) then ! no perturbations
       call MAPL_TimerOff(MAPL, "phase2_Initialize")
       call MAPL_TimerOff(MAPL, "TOTAL")
       RETURN_(ESMF_SUCCESS)
    end if

    ! MPI
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    IAmRoot = MAPL_Am_I_Root(vm)

    !if (IAmRoot) print *, trim(Iam)//':: run'

    ! Get number of land tiles
    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)
    call MAPL_LocStreamGet(locstream, NT_LOCAL=land_nt_local,rc=status)
    VERIFY_(status)

    ! Pointers to mapl internals
    if( internal%isCubedSphere) then
       fpert_ntrmdt => internal%fpert_ntrmdt
       ppert_ntrmdt => internal%ppert_ntrmdt
       pert_rseed_r8 => internal%pert_rseed_r8
    else
       call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=MINTERNAL, rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, fpert_ntrmdt, 'fpert_ntrmdt', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, ppert_ntrmdt, 'ppert_ntrmdt', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, pert_rseed_r8, 'pert_rseed', rc=status)
       VERIFY_(status)
    endif


    ! Convert pert_rseed_r8 to integer
    allocate(pert_rseed(size(pert_rseed_r8)), source=0, stat=status)
    VERIFY_(status)
    pert_rseed = nint(pert_rseed_r8)

    ! Allocate perturbation arrays on grid
    allocate(fpert_grid(n_lon,n_lat, internal%ForcePert%npert), source=MAPL_UNDEF, stat=status)
    VERIFY_(status)
    allocate(ppert_grid(n_lon,n_lat, internal%PrognPert%npert), source=MAPL_UNDEF, stat=status)
    VERIFY_(status)

    ! Get pertubations on the underlying grid and convert grid data to tile data
    !
    ! -ForcePert-
    !
    ! adjust mean after cold start  (if fpert_ntrmdt is from restart file,
    !   mean was adjusted in ApplyForcePert before restart was written)
    if (COLDSTART)                                                              &
         fpert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%ForcePert%npert) =         &
         fpert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%ForcePert%npert) +         &
         fpert_enavg(:,:,:)

    call get_pert(                                                              &
         internal%ForcePert%npert,                                              &
         internal%ForcePert%fft_npert,                                          &
         1,                                                                     &
         internal%pgrid_l, internal%pgrid_f,                                    &
         real(internal%ForcePert%dtstep),                                       &
         internal%ForcePert%param,                                              &
         pert_rseed,                                                            &
         fpert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%ForcePert%npert),          &
         fpert_grid,                                                            &
         initialize_rseed=.false.,                                              &
         initialize_ntrmdt=.false.,                                             &
         ! propagate_pert is NOT called
         diagnose_pert_only=.true.                                              &
         )

    do ipert=1,internal%ForcePert%npert
       call grid2tile( internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), &
                fpert_grid(:,:,ipert), internal%ForcePert%DataPrv(:,ipert))
       call tile_mask_grid(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), fpert_ntrmdt(lon1:lon2,lat1:lat2,ipert))
    end do
    internal%ForcePert%DataNxt = internal%ForcePert%DataPrv

    ! -PrognPert-
    !
    ! adjust mean after cold start  (if ppert_ntrmdt is from restart file,
    !   mean was adjusted in ApplyPrognPert before restart was written)
    if (COLDSTART)                                                              &
         ppert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%PrognPert%npert) =         &
         ppert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%PrognPert%npert) +         &
         ppert_enavg(:,:,:)
    
    call get_pert(                                                              &
         internal%PrognPert%npert,                                              &
         internal%PrognPert%fft_npert,                                          &
         1,                                                                     &
         internal%pgrid_l, internal%pgrid_f,                                    &
         real(internal%PrognPert%dtstep),                                       &
         internal%PrognPert%param,                                              &
         pert_rseed,                                                            &
         ppert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%PrognPert%npert),          &
         ppert_grid,                                                            &
         initialize_rseed=.false.,                                              &
         initialize_ntrmdt=.false.,                                             &
        ! propagate_pert is NOT called
         diagnose_pert_only=.true.                                              &
         )

    do ipert=1,internal%PrognPert%npert
       call grid2tile( internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), ppert_grid(:,:,ipert), &
                internal%PrognPert%DataPrv(:,ipert))
       call tile_mask_grid(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), ppert_ntrmdt(lon1:lon2,lat1:lat2,ipert))
    end do
    internal%PrognPert%DataNxt = internal%PrognPert%DataPrv

    ! Update the r8 version of pert_rseed
    pert_rseed_r8 = real(pert_rseed,kind=ESMF_kind_r8)
    pert_iseed(:,internal%ens_id+1-FIRST_ENS_ID) = pert_rseed

    ! Clean up
    if (allocated(fpert_grid)) then
       deallocate(fpert_grid, stat=status)
       VERIFY_(status)
    end if
    if (allocated(ppert_grid)) then
       deallocate(ppert_grid, stat=status)
       VERIFY_(status)
    end if
    if (allocated(pert_rseed)) then ! integer version of MINTERNAL state
       deallocate(pert_rseed, stat=status)
       VERIFY_(status)
    end if

    ! Turn timers off
    call MAPL_TimerOff(MAPL, "phase2_Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")

    if(internal%ens_id - FIRST_ENS_ID == internal%NUM_ENSEMBLE -1) phase2_initialized = .true. 
    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Phase2_Initialize

  subroutine GenerateRaw_ntrmdt(gc, import, export, clock, rc)

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
    character(len=ESMF_MAXSTR) :: chk_fname
    character(len=ESMF_MAXSTR) :: id_string
    character(len=14)          :: datestamp

    ! ESMF variables
    type(ESMF_Alarm) :: ForcePertAlarm, PrognPertAlarm
    type(ESMF_State) :: MINTERNAL

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    ! Internal private state variables
    type(T_LANDPERT_STATE), pointer :: internal=>null()
    type(LANDPERT_WRAP) :: wrap
    type(TILECOORD_WRAP) :: tcwrap
    type(MAPL_LocStream) :: locstream
 
    ! MAPL internal pointers
    real, pointer :: fpert_ntrmdt(:,:,:)=>null()
    real, pointer :: ppert_ntrmdt(:,:,:)=>null()
    real(kind=ESMF_KIND_R8), pointer :: pert_rseed_r8(:)=>null()

    ! Misc variables
    type(ESMF_VM) :: vm
    logical :: IAmRoot
    integer, allocatable :: pert_rseed(:)
    integer :: m,n_lon,n_lat, land_nt_local, ens_id_width

    integer :: nfpert, nppert, n_tile
    type(tile_coord_type), pointer :: tile_coord_f(:)=>null()    
    type (ESMF_Grid)    :: tilegrid
    integer, pointer    :: mask(:)
    real, allocatable, dimension(:,:) :: tile_data_f, tile_data_p, tile_data_f_all, tile_data_p_all

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // ":: GenerateRaw_ntrmdt"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    call MAPL_TimerOn(MAPL, "GenerateRaw")
    ! Get pointer to the private internal state
    call ESMF_UserCompGetInternalState(gc, 'Landpert_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr

    n_lon=internal%pgrid_l%n_lon
    n_lat=internal%pgrid_l%n_lat

    if (internal%PERTURBATIONS == 0) then ! no perturbations
       call MAPL_TimerOff(MAPL, "GenerateRaw")
       RETURN_(ESMF_SUCCESS)
    end if

    ! Get locstream
    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)

    ! Get number of land tiles
    call MAPL_LocStreamGet(locstream, NT_LOCAL=land_nt_local,rc=status)
    VERIFY_(status)

    ! MPI
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    IAmRoot = MAPL_Am_I_Root(vm)
    ! Get alarm
    call ESMF_ClockGetAlarm(clock, 'ForcePert', ForcePertAlarm, rc=status)
    VERIFY_(status)
    call ESMF_ClockGetAlarm(clock, 'PrognPert', PrognPertAlarm, rc=status)
    VERIFY_(status)
    call MAPL_GetResource( MAPL, ens_id_width,"ENS_ID_WIDTH:", default=4, RC=STATUS)
    VERIFY_(status)
    ! Pointers to mapl internals

    if( internal%isCubedSphere) then
       fpert_ntrmdt => internal%fpert_ntrmdt
       ppert_ntrmdt => internal%ppert_ntrmdt
       pert_rseed_r8 => internal%pert_rseed_r8
    else
       call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=MINTERNAL, rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, fpert_ntrmdt, 'fpert_ntrmdt', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, ppert_ntrmdt, 'ppert_ntrmdt', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, pert_rseed_r8, 'pert_rseed', rc=status)
       VERIFY_(status)
    endif
    ! Convert pert_rseed_r8 to integer
    allocate(pert_rseed(size(pert_rseed_r8)), source=0, stat=status)
    VERIFY_(status)
    pert_rseed = nint(pert_rseed_r8)


    if (MAPL_RecordAlarmIsRinging(MAPL, rc=status) .and. internal%isCubedSphere) then

       call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TileMaskGet(tilegrid,  mask, rc=status)
       VERIFY_(STATUS)
       call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
       VERIFY_(status)

       nfpert = internal%ForcePert%npert
       nppert = internal%PrognPert%npert
       tile_coord_f => tcwrap%ptr%tile_coord_f
       n_tile =  size(tile_coord_f,1)
       ! 1) grid2tile
       allocate(tile_data_f(land_nt_local,nfpert))
       allocate(tile_data_p(land_nt_local,nppert))
       do m = 1, nfpert
         call grid2tile(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), &
                fpert_ntrmdt(lon1:lon2,lat1:lat2,m), tile_data_f(:,m)) 
       enddo
       do m = 1, nppert
         call grid2tile(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), &
                ppert_ntrmdt(lon1:lon2,lat1:lat2,m), tile_data_p(:,m)) 
       enddo
       ! 2) gather tiledata
       if (IAmRoot) then
          allocate(tile_data_f_all(n_tile,nfpert), stat=status)
          VERIFY_(STATUS)
          allocate(tile_data_p_all(n_tile,nppert), stat=status)
          VERIFY_(STATUS)
       else
          allocate(tile_data_f_all(0,nfpert), stat=status)
          VERIFY_(STATUS)
          allocate(tile_data_p_all(0,nppert), stat=status)
          VERIFY_(STATUS)
       end if

       do m = 1, nfpert
         call ArrayGather(tile_data_f(:,m), tile_data_f_all(:,m), tilegrid, mask=mask, rc=status)
         VERIFY_(STATUS)
       enddo
       do m = 1, nppert
         call ArrayGather(tile_data_p(:,m), tile_data_p_all(:,m), tilegrid, mask=mask, rc=status)
         VERIFY_(STATUS)
       enddo
       if (IamRoot) then
       ! 3) tile2grid. simple reverser of grid2tile without weighted averaging/no-data-handling
          do m = 1, nfpert
             call tile2grid_simple( N_tile, tile_coord_f%hash_i_indg, tile_coord_f%hash_j_indg, internal%pgrid_g, tile_data_f_all(:,m), internal%fpert_ntrmdt(:,:,m))
          enddo
          do m = 1, nppert
             call tile2grid_simple( N_tile, tile_coord_f%hash_i_indg, tile_coord_f%hash_j_indg, internal%pgrid_g, tile_data_p_all(:,m), internal%ppert_ntrmdt(:,:,m))
          enddo
       
       ! 4) writing
          call MAPL_DateStampGet(clock, datestamp, rc=status)
          VERIFY_(STATUS)

          id_string=''
          if (internal%NUM_ENSEMBLE > 1) then
            m = len(trim(COMP_NAME))
            id_string = COMP_NAME(m-ens_id_width+1:m)
          endif

          chk_fname = 'landpert'//trim(id_string)//'_internal_checkpoint.'//datestamp//'.nc4'

          call write_pert_checkpoint(trim(chk_fname),internal%fpert_ntrmdt, internal%ppert_ntrmdt, internal%pert_rseed_r8)
       endif
       deallocate(tile_data_f, tile_data_p, tile_data_f_all, tile_data_p_all)
    endif

    if (ESMF_AlarmIsRinging(ForcePertAlarm)) then

    ! -ForcePert-
       call propagate_pert(                                                   &
          internal%ForcePert%fft_npert,                                       &
          1,                                                                  &
          internal%pgrid_l, internal%pgrid_f,                                 &
          real(internal%ForcePert%dtstep),                                    &
          pert_rseed,                                                         &
          internal%ForcePert%param,                                           &
          fpert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%ForcePert%npert),       &
          .false.                                                             &
        )

       if(internal%ens_id == FIRST_ENS_ID ) fpert_enavg(:,:,:)=0.    

       do m = 1,internal%ForcePert%npert
          call tile_mask_grid(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), fpert_ntrmdt(lon1:lon2,lat1:lat2,m))
          if(internal%ForcePert%param(m)%zeromean .and. internal%NUM_ENSEMBLE >2) then
             fpert_enavg(:,:,m)=fpert_enavg(:,:,m)+fpert_ntrmdt(lon1:lon2,lat1:lat2,m)
             if( internal%ens_id - FIRST_ENS_ID == internal%NUM_ENSEMBLE-1) then              
                fpert_enavg(:,:,m) = -fpert_enavg(:,:,m)/real(internal%NUM_ENSEMBLE)
             endif
          endif
       enddo

    endif

    if (ESMF_AlarmIsRinging(PrognPertAlarm)) then

    ! -prognostics-
        call propagate_pert(                                                  &
           internal%PrognPert%fft_npert,                                      &
           1,                                                                 &
           internal%pgrid_l, internal%pgrid_f,                                &
           real(internal%PrognPert%dtstep),                                   &
           pert_rseed,                                                        &
           internal%PrognPert%param,                                          &
           ppert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%PrognPert%npert),      &
           .false.                                                            &
          )

       if(internal%ens_id == FIRST_ENS_ID) ppert_enavg(:,:,:)=0.       

       do m = 1,internal%PrognPert%npert
          call tile_mask_grid(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), ppert_ntrmdt(lon1:lon2,lat1:lat2,m))
          if(internal%PrognPert%param(m)%zeromean .and. internal%NUM_ENSEMBLE >2) then
              ppert_enavg(:,:,m)=ppert_enavg(:,:,m)+ppert_ntrmdt(lon1:lon2,lat1:lat2,m)
              if( internal%ens_id - FIRST_ENS_ID == internal%NUM_ENSEMBLE -1) then                      
                 ppert_enavg(:,:,m) = -ppert_enavg(:,:,m)/real(internal%NUM_ENSEMBLE)
              endif
          endif
       enddo


    endif
    ! Update the r4 version of pert_rseed
    pert_rseed_r8 = real(pert_rseed,kind=ESMF_KIND_R8)
    pert_iseed(:,internal%ens_id+1 - FIRST_ENS_ID) = pert_rseed            
 
    call MAPL_TimerOff(MAPL, "GenerateRaw")
    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine GenerateRaw_ntrmdt

  ! IROTUINE: ApplyForcePert -- Compute and apply perts to MetForcing vars

  ! INTERFACE:

  subroutine ApplyForcePert(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! Compute and apply perturbations to Prognostic variables

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(MAPL_LocStream) :: locstream

    ! ESMF variables
    type(ESMF_VM) :: vm
    type(ESMF_Alarm) :: ForcePertAlarm
    type(ESMF_State) :: MINTERNAL
    type(ESMF_Time) :: ModelTimeCur, ModelTimeNxt, tmpTime
    type(ESMF_TimeInterval) :: ModelTimeStep, ntrvl

    ! Internal private state variables
    type(T_LANDPERT_STATE), pointer :: internal=>null()
    type(LANDPERT_WRAP) :: wrap

    ! MAPL internal pointers
    real, pointer :: fpert_ntrmdt(:,:,:)=>null()
    real(kind=ESMF_KIND_R8), pointer :: pert_rseed_r8(:)=>null()

    ! LDAS variables
    type(date_time_type) :: model_time_nxt
    type(date_time_type) :: fpert_time_prv

    ! Pointers to imports
    real, pointer :: Tair(:)=>null()
    real, pointer :: Qair(:)=>null()
    real, pointer :: Psurf(:)=>null()
    real, pointer :: Rainf_C(:)=>null()
    real, pointer :: Rainf(:)=>null()
    real, pointer :: Snowf(:)=>null()
    real, pointer :: LWdown(:)=>null()
    real, pointer :: SWdown(:)=>null()
    real, pointer :: PARdrct(:)=>null()
    real, pointer :: PARdffs(:)=>null()
    real, pointer :: Wind(:)=>null()
    real, pointer :: RefH(:)=>null()

    ! Perturbed variables
    type(met_force_type), allocatable :: mfPert(:)

    ! Pointers to exports
    real, pointer :: TApert(:)=>null()
    real, pointer :: QApert(:)=>null()
    real, pointer :: PSpert(:)=>null()
    real, pointer :: UUpert(:)=>null()
    real, pointer :: UWINDLMTILEpert(:)=>null()
    real, pointer :: VWINDLMTILEpert(:)=>null()
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

    ! Misc variables
    real, allocatable :: FORCEPERT(:,:)
    real, allocatable :: fpert_grid(:,:,:)
    type(TILECOORD_WRAP) :: tcwrap
    type(tile_coord_type), pointer :: tile_coord(:)=>null()

    integer :: n_lon,n_lat
    integer :: ipert, itile
    integer, allocatable :: pert_rseed(:)
    logical :: IAmRoot
    integer :: land_nt_local
    type(pert_param_type), pointer :: PertParam=>null() ! pert param
    real :: tmpRealArrDim1(1)
    real, allocatable :: tmpreal(:)

    ! Begin...

    ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_ApplyForcePert"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Run_ApplyForcePert")

    ! MPI
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    IAmRoot = MAPL_Am_I_Root(vm)

    ! Get component's internal private state
    call ESMF_UserCompGetInternalState(gc, 'Landpert_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr
    n_lon = internal%pgrid_l%n_lon
    n_lat = internal%pgrid_l%n_lat

   ! Get component's internal tile_coord variable
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tile_coord => tcwrap%ptr%tile_coord

    ! Get locstream
    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)
    
    ! Get number of land tiles
    call MAPL_LocStreamGet(locstream, NT_LOCAL=land_nt_local,rc=status)
    VERIFY_(status)

    ! TODO: this is really, really kludgy and needs to be cleaned
    if (internal%PERTURBATIONS /=0 ) then

       ! Compute FORCEPERT

       ! Get current time
       call ESMF_ClockGet(clock, currTime=ModelTimeCur, rc=status)
       VERIFY_(status)

       ! Compute time stamp of Next model step - convert to LDAS datetime
       call ESMF_ClockGet(clock, timeStep=ModelTimeStep, rc=status)
       VERIFY_(status)
       ModelTimeNxt = ModelTimeCur + ModelTimeStep
       call esmf2ldas(ModelTimeNxt, model_time_nxt, rc=status)
       VERIFY_(status)

       !if(IamRoot) print *, trim(Iam)//'::model_time_nxt: ', date_time_print(model_time_nxt)

       ! Pointers to mapl internals
       if( internal%isCubedSphere) then
          fpert_ntrmdt  => internal%fpert_ntrmdt
          pert_rseed_r8 => internal%pert_rseed_r8
       else
          call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=MINTERNAL, rc=status)
          VERIFY_(status)
          call MAPL_GetPointer(MINTERNAL, fpert_ntrmdt, 'fpert_ntrmdt', rc=status)
          VERIFY_(status)
          call MAPL_GetPointer(MINTERNAL, pert_rseed_r8, 'pert_rseed', rc=status)
          VERIFY_(status)
       endif

       ! Convert pert_rseed_r8 to integer
       allocate(pert_rseed(size(pert_rseed_r8)), source=0, stat=status)
       VERIFY_(status)
       pert_rseed = nint(pert_rseed_r8)

       ! Get alarm
       call ESMF_ClockGetAlarm(clock, 'ForcePert', ForcePertAlarm, rc=status)
       VERIFY_(status)

       ! Allocate and initialize perturbation arrays on grid
       allocate(fpert_grid(n_lon,n_lat, internal%ForcePert%npert), stat=status)
       VERIFY_(status)
       fpert_grid = MAPL_UNDEF

       ! Get forcing perturbations on tiles if alarm is ringing
       if (ESMF_AlarmIsRinging(ForcePertAlarm)) then

          ! -update-times-
          tmpTime = internal%ForcePert%TimeNxt
          internal%ForcePert%TimePrv = tmpTime
          call ESMF_TimeIntervalSet(ntrvl, s=internal%ForcePert%dtstep, rc=status)
          VERIFY_(status)
          internal%ForcePert%TimeNxt = tmpTime + ntrvl

          ! -nxt-pert-data-becomes-prv-
          internal%ForcePert%DataPrv = internal%ForcePert%DataNxt

          ! -get-nxt-forcing-perturbations-on-grid-
          call MAPL_TimerOn(MAPL, '-GetPert')

          ! adjust mean
          fpert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%ForcePert%npert) =                           &
                fpert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%ForcePert%npert)+fpert_enavg(:,:,:)

          call get_pert(                                                           &
               internal%ForcePert%npert,                                           &
               internal%ForcePert%fft_npert,                                       &
               1,                                                                  &
               internal%pgrid_l, internal%pgrid_f,                                 &
               real(internal%ForcePert%dtstep),                                    &
               internal%ForcePert%param,                                           &
               pert_rseed,                                                         &
               fpert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%ForcePert%npert),       &
               fpert_grid,                                                         &
               initialize_rseed=.false.,                                           &
               initialize_ntrmdt=.false.,                                          &
               ! Weiyuan notes: propagate_pert is called in GenerateRaw, not here
               diagnose_pert_only=.true.                                           &
               )

          call MAPL_TimerOff(MAPL, '-GetPert')

          ! -convert-nxt-gridded-perturbations-to-tile-
          call MAPL_TimerOn(MAPL, '-LocStreamTransform')
          do ipert=1,internal%ForcePert%npert
             call grid2tile( internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), fpert_grid(:,:,ipert), &
                  internal%ForcePert%DataNxt(:,ipert))
             call tile_mask_grid(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), fpert_ntrmdt(lon1:lon2,lat1:lat2,ipert))
          end do
       
          call MAPL_TimerOff(MAPL, '-LocStreamTransform')

       end if

       ! Allocate and initialize memory
       allocate(FORCEPERT(land_nt_local, internal%ForcePert%npert), stat=status)
       VERIFY_(status)
       FORCEPERT = MAPL_UNDEF

       ! Interpolate perts on tiles to the end of the model integration time step
       call esmf2ldas(internal%ForcePert%TimePrv, fpert_time_prv, rc=status)
       VERIFY_(status)

       !if(IamRoot) print *, trim(Iam)//'::fpert_time_prv: ', date_time_print(fpert_time_prv)

       call interpolate_pert_to_timestep(                                          &
            model_time_nxt,                                                        &
            fpert_time_prv,                                                        &
            real(internal%ForcePert%dtstep),                                       &
            internal%ForcePert%DataPrv,                                            &
            internal%ForcePert%DataNxt,                                            &
            FORCEPERT(:,1:internal%ForcePert%npert)                                &
            )

    end if

    ! Pointers to imports
    call MAPL_GetPointer(import, Tair, 'Tair', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, Qair, 'Qair', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, Psurf, 'Psurf', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, Rainf_C, 'Rainf_C', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, Rainf, 'Rainf', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, Snowf, 'Snowf', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LWdown, 'LWdown', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SWdown, 'SWdown', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PARdrct, 'PARdrct', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PARdffs, 'PARdffs', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, Wind, 'Wind', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, RefH, 'RefH', rc=status)
    VERIFY_(status)

    ! Allocate memory for perturbed arrays
    allocate(mfPert(land_nt_local), stat=status)
    VERIFY_(status)

    call MAPL_TimerOn(MAPL, '-ApplyPert')

    ! Compute exports
    mfPert%Tair = Tair
    mfPert%Qair = Qair
    mfPert%Rainf_C = Rainf_C
    mfPert%Rainf = Rainf
    mfPert%Snowf = Snowf
    mfPert%LWdown = LWdown
    mfPert%SWdown = SWdown
    mfPert%PARdrct = PARdrct
    mfPert%PARdffs = PARdffs
    mfPert%Wind = Wind

    if (internal%PERTURBATIONS /=0) then

       ! Apply FORCEPERT to MetForcing variables

       do ipert=1,internal%ForcePert%npert
          PertParam => internal%ForcePert%param(ipert) ! shorthand
          select case (trim(PertParam%descr))
          case('pcp')
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%Rainf)
             call repair_forcing(land_nt_local, mfPert, fieldname='Rainf')
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%Rainf_C)
             call repair_forcing(land_nt_local, mfPert, fieldname='Rainf_C')
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%Snowf)
             call repair_forcing(land_nt_local, mfPert, fieldname='Snowf')
          case('sw')
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%SWdown)
             call repair_forcing(land_nt_local, mfPert, fieldname='SWdown')
             ! reichle, 20 Dec 2011 - add perts to "PARdrct" and "PARdffs"
             ! wjiang+reichle, 22 Apr 2021 - "PARdrct" and "PARdffs" now 
             !   backfilled in get_forcing(), arrive here with only "good" values
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%PARdrct)
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%PARdffs)
             ! must repair "PARdrct" and "PARdffs" together
             call repair_forcing(land_nt_local, mfPert, fieldname='PAR')
          case('lw')
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%LWdown)
             call repair_forcing(land_nt_local, mfPert, fieldname='LWdown')
          case('tair')
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%Tair)
             call repair_forcing(land_nt_local, mfPert, fieldname='Tair')
          case('qair')
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%Qair)
             call repair_forcing(land_nt_local, mfPert, fieldname='Qair')
          case('wind')
             call apply_pert(PertParam, FORCEPERT(:,ipert), mfPert%Wind)
             call repair_forcing(land_nt_local, mfPert, fieldname='Wind')
          case default
             RETURN_(ESMF_FAILURE)
          end select
       end do

    end if

    call MAPL_TimerOff(MAPL, '-ApplyPert')

    call MAPL_TimerOn(MAPL, '-MetForcing2Catch')

    ! Pointers to exports (allocate memory)
    call MAPL_GetPointer(export, TApert, 'TApert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QApert, 'QApert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PSpert, 'PSpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, UUpert, 'UUpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, UWINDLMTILEpert, 'UWINDLMTILEpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, VWINDLMTILEpert, 'VWINDLMTILEpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PCUpert, 'PCUpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PLSpert, 'PLSpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNOpert, 'SNOpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRPARpert, 'DRPARpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFPARpert, 'DFPARpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRNIRpert, 'DRNIRpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFNIRpert, 'DFNIRpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRUVRpert, 'DRUVRpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFUVRpert, 'DFUVRpert', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LWDNSRFpert, 'LWDNSRFpert',alloc=.true.,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DZpert, 'DZpert',alloc=.true.,rc=status)
    VERIFY_(status)

    ! Set exports
    TApert = mfPert%Tair
    QApert = mfPert%Qair
    UUpert = mfPert%Wind
    UWINDLMTILEpert = mfPert%Wind
    VWINDLMTILEpert = 0.
    PCUpert = mfPert%Rainf_C
    PLSpert = mfPert%Rainf - mfPert%Rainf_C
    SNOpert = mfPert%Snowf
    ! no pert for psurf
    PSpert = Psurf
    DZpert = RefH
    ! -par-
    ! wjiang+reichle, 22 Apr 2021 - "PARdrct" and "PARdffs" now 
    !   backfilled in get_forcing(), arrive here with only "good" values
    DRPARpert = mfPert%PARdrct
    DFPARpert = mfPert%PARdffs
    ! -nir-and-uvr-
    ! S-V=I+U where S=SWdown, V=DRPAR+DFPAR, I=DRNIR+DFNIR, U=DRUVR+DFUVR
    ! => U=0.5*S-V, I=0.5*S
    allocate(tmpreal(land_nt_local), stat=status)
    VERIFY_(status)
    tmpreal = 0.5*mfPert%SWdown ! I = DRNIR+DFNIR
    DRNIRpert = 0.5*tmpreal
    DFNIRpert = 0.5*tmpreal
   ! tmpreal = tmpreal - (DRPARpert + DFPARpert) ! U = DRUVR+DFUVR
    DRUVRpert = 0.5*tmpreal-DRPARpert
    DFUVRpert = 0.5*tmpreal-DFPARpert
    if (allocated(tmpreal)) deallocate(tmpreal)
    LWDNSRFpert = mfPert%LWdown

    call MAPL_TimerOff(MAPL, '-MetForcing2Catch')

    ! Update the r8 version of pert_rseed
    if (internal%PERTURBATIONS /=0 ) then
       pert_rseed_r8 = real(pert_rseed,kind=ESMF_kind_r8)
       pert_iseed(:,internal%ens_id+1-FIRST_ENS_ID) = pert_rseed           
    endif

    ! Clean up
    if (allocated(mfPert)) then
       deallocate(mfPert, stat=status)
       VERIFY_(status)
    end if
    if (allocated(FORCEPERT)) then
       deallocate(FORCEPERT, stat=status)
       VERIFY_(status)
    end if
    if (allocated(fpert_grid)) then
       deallocate(fpert_grid, stat=status)
       VERIFY_(status)
    end if
    if (allocated(pert_rseed)) then
       deallocate(pert_rseed, stat=status)
       VERIFY_(status)
    end if

    ! Turn timers off
    call MAPL_TimerOff(MAPL, "Run_ApplyForcePert")
    call MAPL_TimerOff(MAPL, "TOTAL")

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine ApplyForcePert


  !BOP

  ! !IROTUINE: ApplyPrognPert -- Compute and apply perts to Prognostic vars

  ! !INTERFACE:

  subroutine ApplyPrognPert(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! Apply perturbations to CATCH GridComp's prognostic variables

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(MAPL_LocStream) :: locstream

    ! ESMF variables
    type(ESMF_VM) :: vm
    type(ESMF_Alarm) :: PrognPertAlarm
    type(ESMF_State) :: MINTERNAL
    type(ESMF_Time) :: ModelTimeCur, ModelTimeNxt, tmpTime
    type(ESMF_TimeInterval) :: ModelTimeStep, ntrvl

    ! Internal private state variables
    type(T_LANDPERT_STATE), pointer :: internal=>null()
    type(LANDPERT_WRAP) :: wrap

    ! MAPL internal pointers
    real, pointer :: ppert_ntrmdt(:,:,:)=>null()
    real(kind=ESMF_KIND_R8), pointer :: pert_rseed_r8(:)=>null()

    ! LDAS variables
    type(date_time_type) :: model_time_nxt
    type(date_time_type) :: ppert_time_prv

    ! Pointers to imports
    real, pointer :: tcPert(:,:)=>null()
    real, pointer :: catdefPert(:)=>null()
    real, pointer :: rzexcPert(:)=>null()
    real, pointer :: srfexcPert(:)=>null()
    real, pointer :: ghtcnt1Pert(:)=>null()
    real, pointer :: ghtcnt2Pert(:)=>null()
    real, pointer :: ghtcnt3Pert(:)=>null()
    real, pointer :: ghtcnt4Pert(:)=>null()
    real, pointer :: ghtcnt5Pert(:)=>null()
    real, pointer :: ghtcnt6Pert(:)=>null()
    real, pointer :: wesnn1Pert(:)=>null()
    real, pointer :: wesnn2Pert(:)=>null()
    real, pointer :: wesnn3Pert(:)=>null()
    real, pointer :: htsnnn1Pert(:)=>null()
    real, pointer :: htsnnn2Pert(:)=>null()
    real, pointer :: htsnnn3Pert(:)=>null()
    real, pointer :: sndzn1Pert(:)=>null()
    real, pointer :: sndzn2Pert(:)=>null()
    real, pointer :: sndzn3Pert(:)=>null()

    ! Misc variables
    real, allocatable :: PROGNPERT(:,:)
    real, allocatable :: ppert_grid(:,:,:)
    type(TILECOORD_WRAP) :: tcwrap
    type(tile_coord_type), pointer :: tile_coord(:)=>null()

    integer :: n_lon,n_lat
    integer :: ipert,ntiles
    integer, allocatable :: pert_rseed(:)
    logical :: IAmRoot
    integer :: land_nt_local
    type(pert_param_type), pointer :: PertParam=>null() ! pert param
    integer :: model_dtstep
    real :: dtmh
    integer :: m
    ! Begin...

    ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_ApplyPrognPert"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)


    ! MPI
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    IAmRoot = MAPL_Am_I_Root(vm)
    
    ! Get component's internal private state
    call ESMF_UserCompGetInternalState(gc, 'Landpert_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr
    ! if no perturbation, do nothing
    if(internal%PERTURBATIONS == 0) then
       RETURN_(ESMF_SUCCESS)
    endif
    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Run_ApplyPrognPert")

    n_lon = internal%pgrid_l%n_lon
    n_lat = internal%pgrid_l%n_lat

   ! Get component's internal tile_coord variable
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tile_coord => tcwrap%ptr%tile_coord
 
    ! Pointers to imports
    call MAPL_GetPointer(import, tcPert, 'TCPert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, catdefPert, 'CATDEFPert',  rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, rzexcPert, 'RZEXCPert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, srfexcPert, 'SRFEXCPert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ghtcnt1Pert, 'GHTCNT1Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ghtcnt2Pert, 'GHTCNT2Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ghtcnt3Pert, 'GHTCNT3Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ghtcnt4Pert, 'GHTCNT4Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ghtcnt5Pert, 'GHTCNT5Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, ghtcnt6Pert, 'GHTCNT6Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, wesnn1Pert, 'WESNN1Pert',  rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, wesnn2Pert, 'WESNN2Pert',  rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, wesnn3Pert, 'WESNN3Pert',  rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, htsnnn1Pert, 'HTSNNN1Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, htsnnn2Pert, 'HTSNNN2Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, htsnnn3Pert, 'HTSNNN3Pert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, sndzn1Pert, 'SNDZN1Pert',  rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, sndzn2Pert, 'SNDZN2Pert',  rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, sndzn3Pert, 'SNDZN3Pert',  rc=status)
    VERIFY_(status)

    call MAPL_TimerOn(MAPL, '-ApplyPert')

    ! Compute PROGNPERT

    ! Get current time
    call ESMF_ClockGet(clock, currTime=ModelTimeCur, rc=status)
    VERIFY_(status)

    ! Compute time stamp of Next model step - convert to LDAS datetime
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep, rc=status)
    VERIFY_(status)
    ModelTimeNxt = ModelTimeCur + ModelTimeStep
    call esmf2ldas(ModelTimeNxt, model_time_nxt, rc=status)
    VERIFY_(status)


    ! Pointers to mapl internals
    if( internal%isCubedSphere) then
       ppert_ntrmdt  => internal%ppert_ntrmdt
       pert_rseed_r8 => internal%pert_rseed_r8
    else
       call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=MINTERNAL, rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, ppert_ntrmdt, 'ppert_ntrmdt', rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, pert_rseed_r8, 'pert_rseed', rc=status)
       VERIFY_(status)
    endif

    ! Convert pert_rseed_r8 to integer
    allocate(pert_rseed(size(pert_rseed_r8)), source=0, stat=status)
    VERIFY_(status)
    pert_rseed = nint(pert_rseed_r8)

    ! Get alarm
    call ESMF_ClockGetAlarm(clock, 'PrognPert', PrognPertAlarm, rc=status)
    VERIFY_(status)

    ! Get locstream
    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)

    ! Get number of land tiles
    call MAPL_LocStreamGet(locstream, NT_LOCAL=land_nt_local,rc=status)
    VERIFY_(status)

    ! Allocate and initialize perturbation arrays on grid
    allocate(ppert_grid(n_lon, n_lat, internal%PrognPert%npert), stat=status)
    VERIFY_(status)
    ppert_grid = MAPL_UNDEF

    ! Get forcing perturbations on tiles if alarm is ringing
    if (ESMF_AlarmIsRinging(PrognPertAlarm)) then

       ! -update-times-
       tmpTime = internal%PrognPert%TimeNxt
       internal%PrognPert%TimePrv = tmpTime
       call ESMF_TimeIntervalSet(ntrvl, s=internal%PrognPert%dtstep, rc=status)
       VERIFY_(status)
       internal%PrognPert%TimeNxt = tmpTime + ntrvl

       ! -nxt-pert-data-becomes-prv-
       internal%PrognPert%DataPrv = internal%PrognPert%DataNxt

       ! -get-nxt-forcing-perturbations-on-grid-
       call MAPL_TimerOn(MAPL, '-GetPert')

       ! adjust mean
       ppert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%PrognPert%npert) =                           &
          ppert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%PrognPert%npert)+ppert_enavg(:,:,:)

       call get_pert(                                                           &
            internal%PrognPert%npert,                                           &
            internal%PrognPert%fft_npert,                                       &
            1,                                                                  &
            internal%pgrid_l, internal%pgrid_f,                                 &
            real(internal%PrognPert%dtstep),                                    &
            internal%PrognPert%param,                                           &
            pert_rseed,                                                         &
            ppert_ntrmdt(lon1:lon2,lat1:lat2,1:internal%PrognPert%npert),       &
            ppert_grid,                                                         &
            initialize_rseed=.false.,                                           &
            initialize_ntrmdt=.false.,                                          &
            ! Weiyuan notes: propagate_pert is called in GenerateRaw, not here
            diagnose_pert_only=.true.                                           &
            )
       call MAPL_TimerOff(MAPL, '-GetPert')

       ! -convert-nxt-gridded-perturbations-to-tile-
       call MAPL_TimerOn(MAPL, '-LocStreamTransform')
       do ipert=1,internal%PrognPert%npert
          call grid2tile( internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), ppert_grid(:,:,ipert), &
                internal%PrognPert%DataNxt(:,ipert))
          call tile_mask_grid(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), ppert_ntrmdt(lon1:lon2,lat1:lat2,ipert))
       end do
       call MAPL_TimerOff(MAPL, '-LocStreamTransform')

    end if

    ! Allocate and initialize memory
    allocate(PROGNPERT(land_nt_local, internal%PrognPert%npert), stat=status)
    VERIFY_(status)
    PROGNPERT = MAPL_UNDEF

    ! Interpolate perts on tiles to the end of the model integration time step
    call esmf2ldas(internal%PrognPert%TimePrv, ppert_time_prv, rc=status)
    VERIFY_(status)

    !if(IamRoot) print *, trim(Iam)//'::ppert_time_prv: ', date_time_print(ppert_time_prv)

    call interpolate_pert_to_timestep(                                          &
         model_time_nxt,                                                        &
         ppert_time_prv,                                                        &
         real(internal%PrognPert%dtstep),                                       &
         internal%PrognPert%DataPrv,                                            &
         internal%PrognPert%DataNxt,                                            &
         PROGNPERT(:,1:internal%PrognPert%npert)                                &
         )
    ! Compute export (perturbed arrays)

    ! -model-dtstep-in-hours-
    call ESMF_TimeIntervalGet(ModelTimeStep, s=model_dtstep)
    VERIFY_(status)
    dtmh = real(model_dtstep)/3600.


    if (internal%PERTURBATIONS /=0) then

       ! Apply PROGNPERT to Prognostic variables
       do ipert=1,internal%PrognPert%npert
          PertParam => internal%PrognPert%param(ipert) ! shorthand
          select case (trim(PertParam%descr))
          case ('catdef')
             call apply_pert(PertParam, PROGNPERT(:,ipert), catdefPert, dtmh)
          case ('rzexc')
             call apply_pert(PertParam, PROGNPERT(:,ipert), rzexcPert, dtmh)
          case ('srfexc')
             call apply_pert(PertParam, PROGNPERT(:,ipert), srfexcPert, dtmh)
          case ('snow')
             _ASSERT(PertParam%typ==1, 'ONLY multiplicative snow perturbations implemented')
             call apply_pert(PertParam, PROGNPERT(:,ipert), wesnn1Pert, dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), wesnn2Pert, dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), wesnn3Pert, dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), htsnnn1Pert, dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), htsnnn2Pert, dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), htsnnn3Pert, dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), sndzn1Pert, dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), sndzn2Pert, dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), sndzn3Pert, dtmh)
          case ('tc')
             call apply_pert(PertParam, PROGNPERT(:,ipert), tcPert(:,1), dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), tcPert(:,2), dtmh)
             call apply_pert(PertParam, PROGNPERT(:,ipert), tcPert(:,4), dtmh)
          case ('ght1')
             call apply_pert(PertParam, PROGNPERT(:,ipert), ghtcnt1Pert, dtmh)
          case ('ght2')
             call apply_pert(PertParam, PROGNPERT(:,ipert), ghtcnt2Pert, dtmh)
          case ('ght3')
             call apply_pert(PertParam, PROGNPERT(:,ipert), ghtcnt3Pert, dtmh)
          case ('ght4')
             call apply_pert(PertParam, PROGNPERT(:,ipert), ghtcnt4Pert, dtmh)
          case ('ght5')
             call apply_pert(PertParam, PROGNPERT(:,ipert), ghtcnt5Pert, dtmh)
          case ('ght6')
             call apply_pert(PertParam, PROGNPERT(:,ipert), ghtcnt6Pert, dtmh)
          case default
             RETURN_(ESMF_FAILURE)
          end select
       end do

!  Removing call to check_cat_progns (wrapper for check_catch_progn) in prep for SMAP L4_SM Version 5.
!  Call was inserted for compatibility of GEOSldas with LDASsa tag used for SMAP L4_SM Version 4 product (Tv4034).
!  Earlier testing without the call (Tv4033) did not result in crashes of catchment() and yielded slightly drier 
!   soil moisture in deserts.
!   - reichle, 17 Jan 2020
!
!       call check_cat_progns(land_nt_local, cat_param, tcPert(:,1), tcPert(:,2), tcPert(:,4),               &
!!          qa1,qa2,qa4, capac                       &
!          catdefPert,                                  &
!          rzexcPert, srfexcPert,                                  &
!         ghtcnt1Pert,ghtcnt2Pert,ghtcnt3Pert,ghtcnt4Pert,ghtcnt5Pert,ghtcnt6Pert, &
!         wesnn1Pert,wesnn2Pert,wesnn3Pert, &
!         htsnnn1Pert,htsnnn2Pert,htsnnn3Pert, &
!         sndzn1Pert, sndzn2Pert,sndzn3Pert)
    end if

    call MAPL_TimerOff(MAPL, '-ApplyPert')

    ! Update the r8 version of pert_rseed
    pert_rseed_r8 = real(pert_rseed,kind=ESMF_kind_r8)
    pert_iseed(:,internal%ens_id+1-FIRST_ENS_ID) = pert_rseed              

    ! Clean up
    if (allocated(PROGNPERT)) then
       deallocate(PROGNPERT, stat=status)
       VERIFY_(status)
    end if
    if (allocated(ppert_grid)) then
       deallocate(ppert_grid, stat=status)
       VERIFY_(status)
    end if
    if (allocated(pert_rseed)) then
       deallocate(pert_rseed, stat=status)
       VERIFY_(status)
    end if

    ! Turn timers off
    call MAPL_TimerOff(MAPL, "Run_ApplyPrognPert")
    call MAPL_TimerOff(MAPL, "TOTAL")

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine ApplyPrognPert

  subroutine Update_pert_rseed(gc,import,export,clock,rc)
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

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()

    ! Internal private state variables
    type(T_LANDPERT_STATE), pointer :: internal=>null()
    type(LANDPERT_WRAP) :: wrap
    type(ESMF_State) :: MINTERNAL
    real(kind=ESMF_KIND_R8), pointer :: pert_rseed_r8(:)=>null()

    ! Begin...

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // ":: Update_pert_rseed"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Get component's private internal state
    call ESMF_UserCompGetInternalState(gc, 'Landpert_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr

    if( internal%isCubedSphere) then
       pert_rseed_r8 => internal%pert_rseed_r8
    else
       call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=MINTERNAL, rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(MINTERNAL, pert_rseed_r8, 'pert_rseed', rc=status)
       VERIFY_(status)
    endif

    pert_rseed_r8(:) = real(pert_iseed(:,internal%ens_id+1-FIRST_ENS_ID),kind=ESMF_KIND_R8)
    
  ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Update_pert_rseed

  !BOP

  ! !IROTUINE: Finalize -- finalize method for LDAS GC

  ! !INTERFACE:

  subroutine Finalize(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! Clean up the private internal state

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name
    character(len=ESMF_MAXSTR) :: chk_fname
    character(len=ESMF_MAXSTR) :: id_string

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()

    ! Internal private state variables
    type(T_LANDPERT_STATE), pointer :: internal=>null()
    type(LANDPERT_WRAP) :: wrap
    type(MAPL_LocStream) :: locstream
    type(TILECOORD_WRAP) :: tcwrap
    integer :: m,n_lon,n_lat, land_nt_local, ens_id_width

    integer :: nfpert, nppert, n_tile
    type(tile_coord_type), pointer :: tile_coord_f(:)=>null()
    type (ESMF_Grid)    :: tilegrid
    integer, pointer    :: mask(:)
    real, allocatable, dimension(:,:) :: tile_data_f, tile_data_p, tile_data_f_all, tile_data_p_all
    ! Begin...


    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Finalize"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    call MAPL_GetResource( MAPL, ens_id_width,"ENS_ID_WIDTH:", default=4, RC=STATUS)
    VERIFY_(status)
    ! Get component's private internal state
    call ESMF_UserCompGetInternalState(gc, 'Landpert_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr

    if ( internal%isCubedSphere .and. internal%PERTURBATIONS /= 0) then
       call MAPL_Get(MAPL, LocStream=locstream,rc=status)
       VERIFY_(status)
       call MAPL_LocStreamGet(locstream, NT_LOCAL=land_nt_local,rc=status)
       VERIFY_(status)
       call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TileMaskGet(tilegrid,  mask, rc=status)
       VERIFY_(STATUS)
       call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
       VERIFY_(status)

       nfpert = internal%ForcePert%npert
       nppert = internal%PrognPert%npert
       tile_coord_f => tcwrap%ptr%tile_coord_f
       n_tile =  size(tile_coord_f,1)
       ! 1) grid2tile
       allocate(tile_data_f(land_nt_local,nfpert))
       allocate(tile_data_p(land_nt_local,nppert))
       do m = 1, nfpert
         call grid2tile(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), &
                internal%fpert_ntrmdt(lon1:lon2,lat1:lat2,m), tile_data_f(:,m)) 
       enddo
       do m = 1, nppert
         call grid2tile(internal%pgrid_l, land_nt_local, internal%i_indgs(:),internal%j_indgs(:), &
                internal%ppert_ntrmdt(lon1:lon2,lat1:lat2,m), tile_data_p(:,m)) 
       enddo
       ! 2) gather tiledata
       if (MAPL_am_I_Root()) then
          allocate(tile_data_f_all(n_tile,nfpert), stat=status)
          VERIFY_(STATUS)
          allocate(tile_data_p_all(n_tile,nppert), stat=status)
          VERIFY_(STATUS)
       else
          allocate(tile_data_f_all(0,nfpert), stat=status)
          VERIFY_(STATUS)
          allocate(tile_data_p_all(0,nppert), stat=status)
          VERIFY_(STATUS)
       end if

       do m = 1, nfpert
         call ArrayGather(tile_data_f(:,m), tile_data_f_all(:,m), tilegrid, mask=mask, rc=status)
         VERIFY_(STATUS)
       enddo
       do m = 1, nppert
         call ArrayGather(tile_data_p(:,m), tile_data_p_all(:,m), tilegrid, mask=mask, rc=status)
         VERIFY_(STATUS)
       enddo

       if (MAPL_am_I_Root()) then
       ! 3) tile2grid 
            ! this step is simply a reverse of grid2tile without any weighted   
          do m = 1, nfpert
             call tile2grid_simple( N_tile, tile_coord_f%hash_i_indg, tile_coord_f%hash_j_indg, internal%pgrid_g, tile_data_f_all(:,m), internal%fpert_ntrmdt(:,:,m))
          enddo
          do m = 1, nppert
             call tile2grid_simple( N_tile, tile_coord_f%hash_i_indg, tile_coord_f%hash_j_indg, internal%pgrid_g, tile_data_p_all(:,m), internal%ppert_ntrmdt(:,:,m))
          enddo

        ! 4) writing
          id_string=''
          if (internal%NUM_ENSEMBLE > 1) then
            m = len(trim(COMP_NAME))
            id_string = COMP_NAME(m-ens_id_width+1:m)
          endif

          chk_fname = 'landpert'//trim(id_string)//'_internal_checkpoint'
          call write_pert_checkpoint(trim(chk_fname),internal%fpert_ntrmdt, internal%ppert_ntrmdt, internal%pert_rseed_r8)
       endif
       deallocate(tile_data_f, tile_data_p, tile_data_f_all, tile_data_p_all)
    endif

    ! Clean up private internal state
    if (associated(internal%ForcePert%param)) then
       deallocate(internal%ForcePert%param, stat=status)
       VERIFY_(status)
    end if
    if (associated(internal%PrognPert%param)) then
       deallocate(internal%PrognPert%param, stat=status)
       VERIFY_(status)
    end if
    if (allocated(internal%ForcePert%DataPrv)) then
       deallocate(internal%ForcePert%DataPrv, stat=status)
       VERIFY_(status)
    end if
    if (allocated(internal%ForcePert%DataNxt)) then
       deallocate(internal%ForcePert%DataNxt, stat=status)
       VERIFY_(status)
    end if
    if (allocated(internal%PrognPert%DataPrv)) then
       deallocate(internal%PrognPert%DataPrv, stat=status)
       VERIFY_(status)
    end if
    if (allocated(internal%PrognPert%DataNxt)) then
       deallocate(internal%PrognPert%DataNxt, stat=status)
       VERIFY_(status)
    end if

    if (allocated(internal%fpert_ntrmdt)) then
       deallocate(internal%fpert_ntrmdt, stat=status)
       VERIFY_(status)
    end if

    if (allocated(internal%ppert_ntrmdt)) then
       deallocate(internal%ppert_ntrmdt, stat=status)
       VERIFY_(status)
    end if

    if (allocated(internal%pert_rseed_r8)) then
       deallocate(internal%pert_rseed_r8, stat=status)
       VERIFY_(status)
    end if

    if (allocated(fpert_enavg)) then
       deallocate(fpert_enavg, stat=status)
       VERIFY_(status)
    end if
    if (allocated(ppert_enavg)) then
       deallocate(ppert_enavg, stat=status)
       VERIFY_(status)
    end if

    ! Call Finalize for every child
    call MAPL_GenericFinalize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Finalize

  subroutine read_pert_rst(rst_fname,fpert, ppert,pert_rseed_r8)
     use netcdf
     character(len=*),intent(in) :: rst_fname
     real,intent(inout) :: fpert(:,:,:)
     real,intent(inout) :: ppert(:,:,:)
     real(kind=ESMF_KIND_R8),intent(inout) :: pert_rseed_r8(:)
     integer :: ncid, varid

     call check( nf90_open(rst_fname, NF90_NOWRITE, ncid) )

  ! Get the varid of the data variable, based on its name.
     call check( nf90_inq_varid(ncid, "fpert_ntrmdt", varid) )
     call check( nf90_get_var(ncid, varid, fpert) )
     call check( nf90_inq_varid(ncid, "ppert_ntrmdt", varid) )
     call check( nf90_get_var(ncid, varid, ppert) )
     call check( nf90_inq_varid(ncid, "pert_rseed", varid) )
     call check( nf90_get_var(ncid, varid, pert_rseed_r8) )

  ! Close the file, freeing all resources.
     call check( nf90_close(ncid) )

     contains
        subroutine check(status)
           integer, intent ( in) :: status
    
           if(status /= nf90_noerr) then 
              print *, trim(nf90_strerror(status))
              stop 1
           end if
       end subroutine check  
  end subroutine

  subroutine write_pert_checkpoint(chk_fname, fpert,ppert, pert_rseed_r8)
     use netcdf
     character(len=*),intent(in) :: chk_fname
     real,intent(inout) :: fpert(:,:,:)
     real,intent(inout) :: ppert(:,:,:)
     real(kind=ESMF_KIND_R8),intent(inout) :: pert_rseed_r8(:)
     character(len=*), parameter  :: SHORT_NAME = "SHORT_NAME"
     character(len=*), parameter  :: LONG_NAME  = "LONG_NAME"
     character(len=*), parameter  :: UNITS      = "UNITS"
     character(len=*), parameter  :: f_SHORT = "fpert_ntrmdt"
     character(len=*), parameter  :: p_SHORT = "pert_ntrmdt"
     character(len=*), parameter  :: s_SHORT = "pert_rseed"
     character(len=*), parameter  :: f_long = "force_pert_intermediate"
     character(len=*), parameter  :: p_long = "progn_pert_intermediate"
     character(len=*), parameter  :: s_long  = "Perturbations_rseed"
     character(len=*), parameter  :: units_ = "1"
     integer :: n_lon, n_lat, nseeds, n_f_max, n_p_max
     integer :: ncid, p_varid,f_varid, s_varid
     integer :: dimids(3),lat_dimid, lon_dimid, seed_dimid, n_f_dimid, n_p_dimid

     n_lon   = size(fpert,1)
     n_lat   = size(fpert,2)
     n_f_max = size(fpert,3)
     n_p_max = size(ppert,3)
     nseeds  = size(pert_rseed_r8)

 ! Create the file. 
     call check( nf90_create(trim(chk_fname), nf90_clobber + NF90_NETCDF4, ncid) )

! Define the dimensions.
     call check( nf90_def_dim(ncid, "latitude",   n_lat, lat_dimid) )
     call check( nf90_def_dim(ncid, "longtitude", n_lon, lon_dimid) )
     call check( nf90_def_dim(ncid, "N_FORCE_MAX", n_f_max, n_f_dimid) )
     call check( nf90_def_dim(ncid, "N_PROGN_MAX", n_p_max, n_p_dimid) )
     call check( nf90_def_dim(ncid, "NRANDSEED",  nseeds, seed_dimid) )

     dimids = (/ lon_dimid, lat_dimid,n_f_dimid /)
     call check( nf90_def_var(ncid, 'fpert_ntrmdt', NF90_REAL,   dimids, f_varid) )
     dimids = (/ lon_dimid, lat_dimid,n_p_dimid /)
     call check( nf90_def_var(ncid, 'ppert_ntrmdt', NF90_REAL,   dimids, p_varid) )
     call check( nf90_def_var(ncid, 'pert_rseed',   NF90_DOUBLE, seed_dimid, s_varid) )

  !   call check( nf90_def_var_deflate(ncid, f_varid, 1, 1, 2))
  !   call check( nf90_def_var_deflate(ncid, p_varid, 1, 1, 2))
  ! Assign attribute
     call check( nf90_put_att(ncid, f_varid, UNITS, units_) )
     call check( nf90_put_att(ncid, p_varid, UNITS, units_) )
     call check( nf90_put_att(ncid, s_varid, UNITS, units_) )

     call check( nf90_put_att(ncid, f_varid, SHORT_NAME, f_short) )
     call check( nf90_put_att(ncid, p_varid, SHORT_NAME, p_short) )
     call check( nf90_put_att(ncid, s_varid, SHORT_NAME, s_short) )

     call check( nf90_put_att(ncid, f_varid, LONG_NAME, f_long) )
     call check( nf90_put_att(ncid, p_varid, LONG_NAME, p_long) )
     call check( nf90_put_att(ncid, s_varid, LONG_NAME, s_long) )

  ! End define mode.
     call check( nf90_enddef(ncid) )

  ! write varaible
     call check( nf90_put_var(ncid, p_varid, ppert) )
     call check( nf90_put_var(ncid, f_varid, fpert) )
     call check( nf90_put_var(ncid, s_varid, pert_rseed_r8) )

  ! Close the file.
     call check( nf90_close(ncid) )

     contains
        subroutine check(status)
           integer, intent ( in) :: status
    
           if(status /= nf90_noerr) then 
              print *, trim(nf90_strerror(status))
              stop 1
           end if
       end subroutine check  
  end subroutine

end module GEOS_LandPertGridCompMod
