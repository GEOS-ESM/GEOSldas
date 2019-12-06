#include "MAPL_Generic.h"

!BOP
! !MODULE: GEOS_LddataatmGridCompMod - Data atmosphere GridComp for Catchment.
module GEOS_MetforceGridCompMod

  ! !USES

  use ESMF
  use MAPL_Mod

  use LDAS_ensdrv_Globals, only: nodata_generic, nodata_tol_generic
  use LDAS_ensdrv_Globals, only: logunit,master_logit,logit
  use LDAS_DateTimeMod, only: date_time_type, date_time_print
  use LDAS_TileCoordType, only: tile_coord_type
  use LDAS_TileCoordType, only: T_TILECOORD_STATE 
  use LDAS_TileCoordType, only: TILECOORD_WRAP
  use LDAS_ForceMod, only: LDAS_GetForcing => get_forcing
  use LDAS_ForceMod, only: LDAS_move_new_force_to_old
  use LDAS_ForceMod, only: FileOpenedHash,GEOS_closefile
  use LDAS_DriverTypes, only: met_force_type, assignment(=)
  use LDAS_ConvertMod, only: esmf2ldas
  use LDAS_InterpMod, only: LDAS_TInterpForcing=>metforcing_tinterp
  !use force_and_cat_progn_pert_types, only: N_FORCE_PERT_MAX

  use StieglitzSnow, only : NUM_DUDP, NUM_DUSV, NUM_DUWT, NUM_DUSD, &
                            NUM_BCDP, NUM_BCSV, NUM_BCWT, NUM_BCSD, &
                            NUM_OCDP, NUM_OCSV, NUM_OCWT, NUM_OCSD, &
                            NUM_SUDP, NUM_SUSV, NUM_SUWT, NUM_SUSD, &
                            NUM_SSDP, NUM_SSSV, NUM_SSWT, NUM_SSSD
  implicit none

  private

  real, parameter      :: daylen = 86400.

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: SetServices

  ! !DESCRIPTION: This GridComp read MetForcing files

  !EOP
  include 'mpif.h'

  ! MetForcing type
  type T_MET_FORCING
     integer :: hinterp ! 1 => Bilin interp
     ! Start/End points of forcing internal
     type(ESMF_Time) :: TimePrv
     type(ESMF_Time) :: TimeNxt
     ! Length of forcing internal
     type(ESMF_TimeInterval) :: ntrvl
     ! File path/tag
     character(len=ESMF_MAXSTR) :: Path
     character(len=ESMF_MAXSTR) :: Tag
     ! Average zenith angle over daylight path of forcing interval
     real, allocatable :: zenav(:)
     ! Met forcing data
     type(met_force_type), pointer, contiguous :: DataPrv(:)
     type(met_force_type), pointer, contiguous :: DataNxt(:)
  end type T_MET_FORCING

  ! Internal state and its wrapper
  type T_DATAATM_STATE
     private
     type(T_MET_FORCING) :: mf
  end type T_DATAATM_STATE
  type DATAATM_WRAP
     type(T_DATAATM_STATE), pointer :: ptr=>null()
  end type DATAATM_WRAP

  !! Wrapper to the tile_coord variable
  !type T_TILECOORD_STATE
  !   type(tile_coord_type), pointer, contiguous :: tile_coord(:)=>null()
  !end type T_TILECOORD_STATE
  !type TILECOORD_WRAP
  !   type(T_TILECOORD_STATE), pointer :: ptr=>null()
  !end type TILECOORD_WRAP

contains

  !BOP

  ! !IROTUINE: SetServices -- Set ESMF services for this component

  ! !INTERFACE:

  subroutine SetServices(gc, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc ! gridded component
    integer, optional                  :: rc ! return code

    ! !DESCRIPTION:
    ! da..da...da....

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! Local variables
    type(T_DATAATM_STATE), pointer :: internal
    type(DATAATM_WRAP) :: wrap

    ! Begin...

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
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         Run,                                                                   &
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
    ! Then, save the pointer to the wrapped internal state in the GridComp
    allocate(internal, stat=status)
    VERIFY_(status)
    wrap%ptr => internal
    call ESMF_UserCompSetInternalState(gc, 'Dataatm_state', wrap, status)
    VERIFY_(status)

    ! Set the state variable specs
    !BOS

    ! !IMPORT STATE:



    ! !EXPORT STATE:

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Tair",                                                   &
         LONG_NAME  = "air_temperature_at_RefH",                                &
         UNITS      = "K",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Qair",                                                   &
         LONG_NAME  = "specific_humidity_at_RefH",                              &
         UNITS      = "kg kg-1",                                                &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Psurf",                                                  &
         LONG_NAME  = "surface_pressure",                                       &
         UNITS      = "Pa",                                                     &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Rainf_C",                                                &
         LONG_NAME  = "convective_rainfall",                                    &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Rainf",                                                  &
         LONG_NAME  = "total_liquid_water_precipitation",                       &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Snowf",                                                  &
         LONG_NAME  = "total_snowfall",                                         &
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
         SHORT_NAME = "LWdown",                                                 &
         LONG_NAME  = "downward_longwave_radiation",                            &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "SWdown",                                                 &
         LONG_NAME  = "downward_shortwave_radiation",                           &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "SWnet",                                                  &
         LONG_NAME  = "downward_net_shortwave_radiation",                       &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PARdrct",                                                &
         LONG_NAME  = "photosynth_active_radiation_direct",                     &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PARdffs",                                                &
         LONG_NAME  = "photosynth_active_radiation_diffuse",                    &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "Wind",                                                   &
         LONG_NAME  = "wind_speed_at_RefH",                                     &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "RefH",                                                   &
         LONG_NAME  = "reference_height_for_Tair_Qair_Wind",                    &
         UNITS      = "m",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)
!
! extra export for GOWIN
!
    call MAPL_AddExportSpec(GC,                          &
         LONG_NAME          = 'dust_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         LONG_NAME          = 'dust_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                  &
         LONG_NAME          = 'dust_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         LONG_NAME          = 'dust_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                 &
         SHORT_NAME         = 'DUSD',                       &
         DIMS               = MAPL_DimsTileOnly,            &
         UNGRIDDED_DIMS     = (/NUM_DUSD/),                 &
         VLOCATION          = MAPL_VLocationNone,           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'black_carbon_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'black_carbon_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'black_carbon_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'black_carbon_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'organic_carbon_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'organic_carbon_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'organic_carbon_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'organic_carbon_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'sulfate_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'sulfate_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'sulfate_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'sulfate_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'sea_salt_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'sea_salt_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'sea_salt_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                         &
         LONG_NAME          = 'sea_salt_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    !EOS

    ! Set profiling timers
    call MAPL_TimerAdd(gc, name="Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Run_GetForcing", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Run_RepairForcing", rc=status)
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

    ! !DESCRIPTION:
    ! da...da..da.

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! ESMF variables
    type(ESMF_Time) :: CurrentTime
    type(ESMF_Alarm) :: MetForcingAlarm
    type(ESMF_TimeInterval) :: Forcing_DT

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()

    ! LDAS variables
    type(date_time_type) :: force_time_prv

    ! MetForcing variable
    type(T_MET_FORCING) :: mf

    ! Internal private state variables
    type(T_DATAATM_STATE), pointer :: internal=>null()
    type(DATAATM_WRAP) :: wrap
    type(TILECOORD_WRAP) :: tcwrap
    type(tile_coord_type), pointer :: tile_coord(:)=>null()

    ! Misc variables
    integer :: land_nt_local
    integer :: ForceDtStep
    type(met_force_type) :: mf_nodata
    logical :: MERRA_file_specs,GEOS_Forcing

    integer :: AEROSOL_DEPOSITION
    type(MAPL_LocStream) :: locstream

    ! Begin...

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Initialize"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Initialize")

    ! Get current time
    call ESMF_ClockGet(clock, currTime=CurrentTime, rc=status)
    VERIFY_(status)

    ! Get component's internal private state
    call ESMF_UserCompGetInternalState(gc, 'Dataatm_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr

     ! Get component's internal tile_coord variable
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tile_coord => tcwrap%ptr%tile_coord

    ! Number of land tiles (on local PE)
    call MAPL_Get(MAPL, LocStream=locstream)
    VERIFY_(status)
    call MAPL_LocStreamGet(                                                     &
         locstream,                                                             &
         NT_LOCAL=land_nt_local,                                                &
         rc=status                                                              &
         )
    VERIFY_(status)

    call MAPL_GetResource ( MAPL, AEROSOL_DEPOSITION, Label="AEROSOL_DEPOSITION:", &
         DEFAULT=0, RC=STATUS)

    ! Get MetForcing values and put them in Ldas' internal state
    ! Get resources needed to call LDAS_ForceMod::get_forcing()
    ! - hinterp=1 => Bilinear Interpolation -
    call MAPL_GetResource(                                                      &
         MAPL,                                                                  &
         mf%hinterp,                                                            &
         'MET_HINTERP:',                                                        &
         default=1,                                                             &
         rc=status                                                              &
         )
    VERIFY_(status)
    ! -previous/next-times-
    mf%TimePrv = CurrentTime
    mf%TimeNxt = CurrentTime
    ! -time-interval-
    call MAPL_GetResource(                                                      &
         MAPL,                                                                  &
         ForceDtStep,                                                           &
         'FORCE_DTSTEP:',                                                       &
         default=3600,                                                          &
         rc=status                                                              &
         )
    VERIFY_(status)
    call ESMF_TimeIntervalSet(Forcing_DT, s=ForceDtStep, rc=status)
    VERIFY_(status)
    mf%ntrvl = Forcing_DT
    ! -path-
    call MAPL_GetResource(MAPL, mf%Path, 'MET_PATH:', rc=status)
    VERIFY_(status)
    ! -tag-
    call MAPL_GetResource(MAPL, mf%Tag, 'MET_TAG:', rc=status)
    VERIFY_(status)
    ! -allocate-memory-for-metforcing-data-
    mf_nodata = nodata_generic
    allocate(mf%DataPrv(land_nt_local), source=mf_nodata, stat=status)
    VERIFY_(status)
    allocate(mf%DataNxt(land_nt_local), source=mf_nodata, stat=status)
    VERIFY_(status)
    ! -allocate-memory-for-avg-zenith-angle
    allocate(mf%zenav(land_nt_local), source=nodata_generic, stat=status)
    VERIFY_(status)
    ! Put MetForcing in Ldas' pvt internal state
    internal%mf = mf

    ! Create alarm for MetForcing
    ! -create-nonsticky-alarm-
    MetForcingAlarm = ESMF_AlarmCreate(                                         &
         clock,                                                                 &
         name='MetForcing',                                                     &
         ringTime=CurrentTime,                                                  &
         ringInterval=Forcing_DT,                                               &
         ringTimeStepCount=1,                                                   &
         sticky=.false.,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)

    ! Get "prv" forcing
    ! -convert-mf%TimePrv-to-LDAS-datetime-
    call esmf2ldas(mf%TimePrv, force_time_prv, rc=status)
    VERIFY_(status)
    ! -now-get-the-initial-forcings-
    call LDAS_GetForcing(                                                       &
         force_time_prv,                                                        &
         ForceDtStep,                                                           &
         internal%mf%Path,                                                      &
         internal%mf%Tag,                                                       &
         land_nt_local,                                                         &
         tile_coord,                                                            &
         internal%mf%hinterp,                                                   &
         MERRA_file_specs,                                                      &
         GEOS_Forcing,                                                          &
         internal%mf%DataNxt,                                                   &
         AEROSOL_DEPOSITION,                                                    &
         .true.                                                                 &
         )
    VERIFY_(status)
    call LDAS_move_new_force_to_old(internal%mf%DataNxt,internal%mf%DataPrv,   &
           MERRA_file_specs,GEOS_Forcing, AEROSOL_DEPOSITION)

    ! DataPrv is not well defined here
    ! print *, 'prv%tair max/min: ', maxval(internal%mf%DataPrv%Tair), minval(internal%mf%DataPrv%Tair)
    ! print *, 'nxt%tair max/min: ', maxval(internal%mf%DataNxt%Tair), minval(internal%mf%DataNxt%Tair)

    ! Turn timer off
    call MAPL_TimerOff(MAPL, "Initialize")

    ! Call Initialize for every child
    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    ! End
    call MAPL_TimerOff(MAPL, "TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


  !BOP

  ! !IROTUINE: Run_GetForcing - a wrapper around Rolf's get_forcing()

  ! !INTERFACE:

  subroutine Run(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! Reads met_forcing files.

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! ESMF variables
    type(ESMF_VM) :: vm
    type(ESMF_Time) :: ModelTimeCur, tmpTime, ModelTimeNxt
    type(ESMF_Alarm) :: MetForcingAlarm
    type(ESMF_TimeInterval) :: ModelTimeStep

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null() ! MAPL obj
    type(MAPL_LocStream) :: locstream
    type(MAPL_SunOrbit) :: orbit

    ! LDAS variables
    type(date_time_type) :: force_time_prv, force_time_nxt, model_time_nxt

    ! Private internal state variables
    type(T_DATAATM_STATE), pointer :: internal=>null()
    type(DATAATM_WRAP) :: wrap
    type(TILECOORD_WRAP) :: tcwrap ! LDAS' tile_coord variable
    type(tile_coord_type), pointer :: tile_coord(:)

    ! Misc variables
    integer :: land_nt_local ! number of LAND tiles in local PE
    integer :: comm
    logical :: IAmRoot
    integer :: fdtstep
    integer :: YEAR, DAY_OF_YEAR, SEC_OF_DAY,n
    real, pointer :: LandTileLats(:)
    real, pointer :: LandTileLons(:)
    real, allocatable :: zth(:), slr(:), zth_tmp(:)
    type(met_force_type), allocatable :: mfDataNtp(:)
    type(met_force_type), pointer :: DataTmp(:)=>null()
    real, allocatable :: tmpreal(:)
    type(met_force_type) :: mf_nodata

    logical :: MERRA_file_specs,GEOS_Forcing
    integer :: AEROSOL_DEPOSITION
    ! Export pointers
    real, pointer :: Tair(:)=>null()
    real, pointer :: Qair(:)=>null()
    real, pointer :: Psurf(:)=>null()
    real, pointer :: Rainf_C(:)=>null()
    real, pointer :: Rainf(:)=>null()
    real, pointer :: Snowf(:)=>null()
    real, pointer :: RainfSnowf(:)=>null()
    real, pointer :: LWdown(:)=>null()
    real, pointer :: SWdown(:)=>null()
    real, pointer :: SWnet(:)=>null()
    real, pointer :: PARdrct(:)=>null()
    real, pointer :: PARdffs(:)=>null()
    real, pointer :: Wind(:)=>null()
    real, pointer :: RefH(:)=>null()

    real,pointer :: DUDP(:,:)=>null()
    real,pointer :: DUSV(:,:)=>null()
    real,pointer :: DUWT(:,:)=>null()
    real,pointer :: DUSD(:,:)=>null()
    real,pointer :: BCDP(:,:)=>null()
    real,pointer :: BCSV(:,:)=>null()
    real,pointer :: BCWT(:,:)=>null()
    real,pointer :: BCSD(:,:)=>null()
    real,pointer :: OCDP(:,:)=>null()
    real,pointer :: OCSV(:,:)=>null()
    real,pointer :: OCWT(:,:)=>null()
    real,pointer :: OCSD(:,:)=>null()
    real,pointer :: SUDP(:,:)=>null()
    real,pointer :: SUSV(:,:)=>null()
    real,pointer :: SUWT(:,:)=>null()
    real,pointer :: SUSD(:,:)=>null()
    real,pointer :: SSDP(:,:)=>null()
    real,pointer :: SSSV(:,:)=>null()
    real,pointer :: SSWT(:,:)=>null()
    real,pointer :: SSSD(:,:)=>null()
    ! Begin...

    ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_GetForcing"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Run_GetForcing")

    ! MPI stuff - communicator, root etc.
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    call ESMF_VmGet(vm, mpicommunicator=comm, rc=status)
    VERIFY_(status)
    IAmRoot = MAPL_Am_I_Root(vm)

    ! Get current time
    call ESMF_ClockGet(clock, currTime=ModelTimeCur, rc=status)
    VERIFY_(status)

    ! Get component's internal private state
    call ESMF_UserCompGetInternalState(gc, 'Dataatm_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr
 
   call MAPL_GetResource ( MAPL, AEROSOL_DEPOSITION, Label="AEROSOL_DEPOSITION:", &
         DEFAULT=1, RC=STATUS)

    ! Get number of tiles, tile lats/lons from LocStream
    call MAPL_Get(MAPL, LocStream=locstream)
    VERIFY_(status)
    call MAPL_LocStreamGet(                                                     &
         locstream,                                                             &
         NT_LOCAL=land_nt_local,                                                &
         TILELATS=LandTileLats,                                                 &
         TILELONS=LandTileLons,                                                 &
         rc=status                                                              &
         )
    VERIFY_(status)

    ! Sun's orbit
    call MAPL_Get(MAPL, orbit=orbit)

    ! Allocate memory for zenith angle
    allocate(zth(land_nt_local), source=nodata_generic, stat=status)
    VERIFY_(status)
    allocate(slr(land_nt_local), source=nodata_generic, stat=status)
    VERIFY_(status)
    allocate(zth_tmp(land_nt_local), source=nodata_generic, stat=status)
    VERIFY_(status)

    ! Convert forcing time interval to seconds
    call ESMF_TimeIntervalGet(internal%mf%ntrvl, s=fdtstep, rc=status)

    ! MetForcing alarm
    call ESMF_ClockGetAlarm(clock, 'MetForcing', MetForcingAlarm, rc=status)
    VERIFY_(status)

    ! Get component's internal tile_coord variable
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tile_coord => tcwrap%ptr%tile_coord

    ! Time stamp of next model step
    ! -get-model-time-step-
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep)
    VERIFY_(status)
    ! -time-stamp-
    ModelTimeNxt = ModelTimeCur + ModelTimeStep

    ! Get forcing data if MetForcing alarm is ringing
    if (ESMF_AlarmIsRinging(MetForcingAlarm)) then

       ! -update-forcing-times-
       tmpTime = internal%mf%TimeNxt
       internal%mf%TimePrv = tmpTime
       internal%mf%TimeNxt = tmpTime + internal%mf%ntrvl

       ! -update-forcing-data-
       ! -swap-DataPrv-and-DataNxt-
       DataTmp => internal%mf%DataPrv
       internal%mf%DataPrv => internal%mf%DataNxt
       internal%mf%DataNxt => DataTmp
       nullify(DataTmp)

       ! -convert-mf%TimeNxt-to-LDAS-datetime-
       call esmf2ldas(internal%mf%TimeNxt, force_time_nxt, rc=status)
       VERIFY_(status)

       call LDAS_GetForcing(                                                    &
            force_time_nxt,                                                     &
            fdtstep,                                                            &
            internal%mf%Path,                                                   &
            internal%mf%Tag,                                                    &
            land_nt_local,                                                      &
            tile_coord,                                                         &
            internal%mf%hinterp,                                                &
            MERRA_file_specs,                                                   &
            GEOS_Forcing,                                                       &
            internal%mf%DataNxt,                                                &
            AEROSOL_DEPOSITION,                                                 &
            .false.                                                             &
            )
       call LDAS_move_new_force_to_old(internal%mf%DataNxt,internal%mf%DataPrv, &
           MERRA_file_specs,GEOS_Forcing,AEROSOL_DEPOSITION)

       !if(master_logit) write(logunit,*) trim(Iam)//'::force_time_nxt: ', date_time_print(force_time_nxt)

       ! -compute-average-zenith-angle-over-daylight-part-of-forcing-interval-
       call MAPL_SunGetInsolation(                                              &
            LandTileLons,                                                       &
            LandTileLats,                                                       &
            orbit,                                                              &
            zth_tmp,                                                                &
            slr,                                                                &
            currTime=internal%mf%TimePrv,                                       &
            INTV=internal%mf%ntrvl,                                             &
            ZTHB=internal%mf%zenav,                                             &
            STEPSIZE=150.0,                                                     &
            rc=status                                                           &
            )
       VERIFY_(STATUS) 

      ! call ESMF_TimeGet(internal%mf%TimePrv, YY=YEAR, S=SEC_OF_DAY, &
      !          dayOfYear=DAY_OF_YEAR, RC=STATUS)
      ! VERIFY_(STATUS) 

      ! call zenith(DAY_OF_YEAR,SEC_OF_DAY,fdtstep,ModelTimeStep,land_nt_local,tile_coord%com_lon,                                                    &
      !         tile_coord%com_lat,internal%mf%zenav)


       ! -checks-on-computed-zenith-angles-
       if (any(internal%mf%zenav<0)) then
          RETURN_(ESMF_FAILURE)
       end if

    end if

    !if(master_logit) write(logunit,*) trim(Iam)//'::zenav max/min: ', maxval(internal%mf%zenav), minval(internal%mf%zenav)
    !if(logit) write(logunit,*) trim(Iam)//'::zenav max/min: ', maxval(internal%mf%zenav), minval(internal%mf%zenav)

    ! Compute zenith angle at the next time step
    call MAPL_SunGetInsolation(                                                 &
         LandTileLons,                                                          &
         LandTileLats,                                                          &
         orbit,                                                                 &
         zth_tmp,                                                               &
         slr,                                                                   &
         INTV=ModelTimeStep,                                                    &
         ZTHB=zth,                                                              &
         currTime=ModelTimeCur,                                                 &
         STEPSIZE=150.0,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)

    !call ESMF_TimeGet(ModelTimeNxt, YY=YEAR, S=SEC_OF_DAY, &
    !          dayOfYear=DAY_OF_YEAR, RC=STATUS)
    !VERIFY_(STATUS)
    !do n=1, land_nt_local
    !  call solar(tile_coord(n)%com_lon,tile_coord(n)%com_lat, DAY_OF_YEAR,SEC_OF_DAY,zth(n),slr(n))
    !enddo

    if (any(zth<0.)) then
       RETURN_(ESMF_FAILURE)
    end if

    !if(master_logit) write(logunit,*)  trim(Iam)//'::zth max/min: ', maxval(zth), minval(zth)

    ! -convert-mf%TimePrv-to-LDAS-datetime-
    call esmf2ldas(internal%mf%TimePrv, force_time_prv, rc=status)
    VERIFY_(status)

    ! -convert-ModelTimeNxt-to-LDAS-datetime-
    call esmf2ldas(ModelTimeNxt, model_time_nxt, rc=status)

    !if(master_logit) write(logunit,*) trim(Iam)//'::force_time_prv: ', date_time_print(force_time_prv)

    !if(master_logit) write(logunit,*) trim(Iam)//'::model_time_nxt: ', date_time_print(model_time_nxt)

    ! Allocate memory for interpolated MetForcing data
    mf_nodata = nodata_generic
    allocate(mfDataNtp(land_nt_local), source=mf_nodata, stat=status)
    VERIFY_(status)

    ! Interpolate MetForcing data to the end of model integration time step
    call LDAS_TInterpForcing(                                                   &
         tile_coord%com_lon,                                                    &
         tile_coord%com_lat,                                                    &
         zth,                                                                   &
         internal%mf%zenav,                                                     &
         force_time_prv,                                                        &
         model_time_nxt,                                                        &
         fdtstep,                                                               &
         internal%mf%DataPrv,                                                   &
         internal%mf%DataNxt,                                                   &
         mfDataNtp,                                                             &
         AEROSOL_DEPOSITION,                                                    &
         rc=status                                                              &
         )
    VERIFY_(status)
    !if(master_logit) write(logunit,*) trim(Iam)//'::mf_ntp%tair max/min: ', maxval(mfDataNtp%Tair), minval(mfDataNtp%Tair)

    ! Pointers to exports (allocate memory)
    call MAPL_GetPointer(export, Tair, 'Tair', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, Qair, 'Qair', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, Psurf, 'Psurf', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, Rainf_C, 'Rainf_C', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, Rainf, 'Rainf', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, Snowf, 'Snowf', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RainfSnowf, 'RainfSnowf', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LWdown, 'LWdown', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SWdown, 'SWdown', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SWnet, 'SWnet', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PARdrct, 'PARdrct', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PARdffs, 'PARdffs', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, Wind, 'Wind', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RefH, 'RefH', alloc=.true., rc=status)
    VERIFY_(status)

    if (AEROSOL_DEPOSITION /=0 ) then
       call MAPL_GetPointer(export, DUDP, 'DUDP' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, DUSV, 'DUSV' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, DUWT, 'DUWT' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, DUSD, 'DUSD' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, BCDP, 'BCDP' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, BCSV, 'BCSV' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, BCWT, 'BCWT' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, BCSD, 'BCSD' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, OCDP, 'OCDP' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, OCSV, 'OCSV' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, OCWT, 'OCWT' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, OCSD, 'OCSD' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, SUDP, 'SUDP' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, SUSV, 'SUSV' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, SUWT, 'SUWT' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, SUSD, 'SUSD' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, SSDP, 'SSDP' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, SSSV, 'SSSV' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, SSWT, 'SSWT' , alloc=.true., rc=status)
       VERIFY_(status)
       call MAPL_GetPointer(export, SSSD, 'SSSD' , alloc=.true., rc=status)
       VERIFY_(status)
    endif ! AEROSOL_DEPOSITION /=0


    ! Set exports
    Tair = mfDataNtp%Tair
    Qair = mfDataNtp%Qair
    Psurf = mfDataNtp%Psurf
    Rainf_C = mfDataNtp%Rainf_C
    Rainf = mfDataNtp%Rainf
    Snowf = mfDataNtp%Snowf
    ! *daylen convert [kg/m2/s] into [kg/m2/day]
    !RainfSnowf= (Rainf+Snowf)*daylen 
    RainfSnowf= Rainf+Snowf 
    LWdown = mfDataNtp%LWdown
    SWdown = mfDataNtp%SWdown
    SWnet = mfDataNtp%SWnet
    PARdrct = mfDataNtp%PARdrct
    PARdffs = mfDataNtp%PARdffs
    Wind = mfDataNtp%Wind
    RefH = mfDataNtp%RefH

    if(AEROSOL_DEPOSITION /=0) then
      DUDP(:, 1) = mfDataNtp%DUDP001
      DUDP(:, 2) = mfDataNtp%DUDP002
      DUDP(:, 3) = mfDataNtp%DUDP003
      DUDP(:, 4) = mfDataNtp%DUDP004
      DUDP(:, 5) = mfDataNtp%DUDP005
      DUSV(:, 1) = mfDataNtp%DUSV001
      DUSV(:, 2) = mfDataNtp%DUSV002
      DUSV(:, 3) = mfDataNtp%DUSV003
      DUSV(:, 4) = mfDataNtp%DUSV004
      DUSV(:, 5) = mfDataNtp%DUSV005
      DUWT(:, 1) = mfDataNtp%DUWT001
      DUWT(:, 2) = mfDataNtp%DUWT002
      DUWT(:, 3) = mfDataNtp%DUWT003
      DUWT(:, 4) = mfDataNtp%DUWT004
      DUWT(:, 5) = mfDataNtp%DUWT005
      DUSD(:, 1) = mfDataNtp%DUSD001
      DUSD(:, 2) = mfDataNtp%DUSD002
      DUSD(:, 3) = mfDataNtp%DUSD003
      DUSD(:, 4) = mfDataNtp%DUSD004
      DUSD(:, 5) = mfDataNtp%DUSD005
      BCDP(:, 1) = mfDataNtp%BCDP001
      BCDP(:, 2) = mfDataNtp%BCDP002
      BCSV(:, 1) = mfDataNtp%BCSV001
      BCSV(:, 2) = mfDataNtp%BCSV002
      BCWT(:, 1) = mfDataNtp%BCWT001
      BCWT(:, 2) = mfDataNtp%BCWT002
      BCSD(:, 1) = mfDataNtp%BCSD001
      BCSD(:, 2) = mfDataNtp%BCSD002
      OCDP(:, 1) = mfDataNtp%OCDP001
      OCDP(:, 2) = mfDataNtp%OCDP002
      OCSV(:, 1) = mfDataNtp%OCSV001
      OCSV(:, 2) = mfDataNtp%OCSV002
      OCWT(:, 1) = mfDataNtp%OCWT001
      OCWT(:, 2) = mfDataNtp%OCWT002
      OCSD(:, 1) = mfDataNtp%OCSD001
      OCSD(:, 2) = mfDataNtp%OCSD002
      SUDP(:, 1) = mfDataNtp%SUDP003
      SUSV(:, 1) = mfDataNtp%SUSV003
      SUWT(:, 1) = mfDataNtp%SUWT003
      SUSD(:, 1) = mfDataNtp%SUSD003
      SSDP(:, 1) = mfDataNtp%SSDP001
      SSDP(:, 2) = mfDataNtp%SSDP002
      SSDP(:, 3) = mfDataNtp%SSDP003
      SSDP(:, 4) = mfDataNtp%SSDP004
      SSDP(:, 5) = mfDataNtp%SSDP005
      SSSV(:, 1) = mfDataNtp%SSSV001
      SSSV(:, 2) = mfDataNtp%SSSV002
      SSSV(:, 3) = mfDataNtp%SSSV003
      SSSV(:, 4) = mfDataNtp%SSSV004
      SSSV(:, 5) = mfDataNtp%SSSV005
      SSWT(:, 1) = mfDataNtp%SSWT001
      SSWT(:, 2) = mfDataNtp%SSWT002
      SSWT(:, 3) = mfDataNtp%SSWT003
      SSWT(:, 4) = mfDataNtp%SSWT004
      SSWT(:, 5) = mfDataNtp%SSWT005
      SSSD(:, 1) = mfDataNtp%SSSD001
      SSSD(:, 2) = mfDataNtp%SSSD002
      SSSD(:, 3) = mfDataNtp%SSSD003
      SSSD(:, 4) = mfDataNtp%SSSD004
      SSSD(:, 5) = mfDataNtp%SSSD005
    endif ! AEROSOL_DEPOSITION /=0
    ! Clean up
    if (allocated(mfDataNtp)) then
       deallocate(mfDataNtp, stat=status)
       VERIFY_(status)
    end if
    if (allocated(zth)) then
       deallocate(zth, stat=status)
       VERIFY_(status)
    end if
    if (allocated(slr)) then
       deallocate(slr, stat=status)
       VERIFY_(status)
    end if

    ! Turn timers off
    call MAPL_TimerOff(MAPL, "Run_GetForcing")
    call MAPL_TimerOff(MAPL, "TOTAL")

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Run


  !BOP

  ! !IROTUINE: Finalize -- Finalize method for LDAS GridComp

  ! !INTERFACE:

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
    type(T_DATAATM_STATE), pointer :: internal
    type(DATAATM_WRAP) :: wrap
    type(ESMF_Alarm) :: MetForcing
    !external :: GEOS_closefile
    ! Begin...

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Initialize"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Get component's internal private state
    call ESMF_UserCompGetInternalState(gc, 'Dataatm_state', wrap, status)
    VERIFY_(status)
    internal => wrap%ptr

    call FileOpenedHash%free(GEOS_closefile,.true.)

    ! Clean-up private internal state
    if (allocated(internal%mf%zenav)) then
       deallocate(internal%mf%zenav)
    end if
    if (associated(internal%mf%DataPrv)) then
       deallocate(internal%mf%DataPrv)
    end if
    if (associated(internal%mf%DataNxt)) then
       deallocate(internal%mf%DataNxt)
    end if

    ! Destroy MetForcingAlarm
    ! -get-the-alarm
    call ESMF_ClockGetAlarm(clock, 'MetForcing', MetForcing, rc=status)
    VERIFY_(status)
    ! -destroy-it-
    call ESMF_AlarmDestroy(MetForcing, rc=status)
    VERIFY_(status)

    ! Call Finalize for every child
    call MAPL_GenericFinalize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Finalize

end module GEOS_MetforceGridCompMod
