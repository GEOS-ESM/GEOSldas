#include "MAPL_Generic.h"

!=============================================================================
module GEOS_LandAssimGridCompMod

!BOP
! !DESCRIPTION:
!
!   {\tt Obs} is a gridded component to
!   {\tt Obs} has no children.

!
! !USES:

  use ESMF
  use MAPL_Mod
  use ESMF_CFIOMOD, only:  ESMF_CFIOstrTemplate
  !USE GEOS_MOD

  use LDAS_TileCoordType, only: tile_coord_type
  use LDAS_TileCoordType, only: grid_def_type
  use LDAS_TileCoordType, only: T_TILECOORD_STATE
  use LDAS_TileCoordType, only: TILECOORD_WRAP

  use enkf_types,      only: obs_type,obs_param_type
  use nr_ran2_gasdev,                   ONLY:     &
       NRANDSEED, init_randseed
  use land_pert_routines, only: get_init_pert_rseed
  use LDAS_ensdrv_mpi, only: mpicomm,numprocs,myid
  use LDAS_ensdrv_mpi, only: master_proc
  use LDAS_ensdrv_mpi, only: MPI_obs_param_type 
  
  use LDAS_DateTimeMod,ONLY: date_time_type 
  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename 
  use LDAS_ensdrv_init_routines, only : GEOS_read_catparam
  use LDAS_ensdrv_Globals, only: logunit

  use LDAS_ConvertMod, ONLY: esmf2ldas
  use LDAS_DriverTypes,                     ONLY:     &
       met_force_type

  use GEOS_LandPertGridCompMod, only: N_force_pert, N_progn_pert
  use GEOS_LandPertGridCompMod, only: progn_pert_param
  use GEOS_LandPertGridCompMod, only: force_pert_param
!  use GEOS_LandPertGridCompMod, only: pert_rseed=>pert_iseed

  use lsm_routines, only: DZGT
  use GEOS_EnsGridCompMod,      only: cat_progn=>catch_progn
  use GEOS_EnsGridCompMod,      only: cat_param=>catch_param
  use mwRTM_types, only: mwRTM_param_type
  use catch_bias_types, only: obs_bias_type
  use catch_bias_types, only: cat_bias_param_type
  use catch_types, only: cat_progn_type
  use catch_types, only: cat_param_type
  use catch_types, only: assignment(=), operator (+), operator (/)
  use clsm_bias_routines, only: initialize_obs_bias
  use clsm_bias_routines, only: read_cat_bias_inputs 

  use clsm_ensupd_upd_routines, only: read_ens_upd_inputs
  use clsm_ensupd_upd_routines, only: finalize_obslog
  use clsm_ensupd_glob_param, only: echo_clsm_ensupd_glob_param
  use clsm_ensupd_enkf_update, only: get_enkf_increments 
  use clsm_ensupd_enkf_update, only: apply_enkf_increments 
  use clsm_ensupd_enkf_update, only: output_incr_etc
  use clsm_ensupd_enkf_update, only: write_smapL4SMaup 
  use clsm_ensdrv_out_routines, only: init_log, GEOS_output_smapL4SMlmc 
  use, intrinsic :: ieee_arithmetic    

implicit none

include 'mpif.h'

private

! !PUBLIC MEMBER FUNCTIONS:

public :: SetServices
!
!EOP
!
integer, parameter :: NUM_SUBTILES = 4
integer           :: NUM_ENSEMBLE
integer           :: FIRST_ENS_ID

type(met_force_type), allocatable :: mfPert_ensavg(:)

type(obs_param_type),pointer :: obs_param(:)=>null()
logical :: need_mwRTM_param
integer :: update_type, dtstep_assim
logical :: centered_update
real    :: xcompact, ycompact
integer :: N_obs_param
logical :: out_obslog
logical :: out_ObsFcstAna
logical :: out_smapL4SMaup
integer :: N_obsbias_max

integer,dimension(:),pointer :: N_catl_vec,low_ind
integer :: N_catf
!reordered tile_coord_rf and mapping l2rf
integer,dimension(:),pointer :: l2rf, rf2l,rf2g, rf2f
type(tile_coord_type), dimension(:), pointer :: tile_coord_rf => null()
integer, allocatable :: Pert_rseed(:,:)
real(kind=ESMF_KIND_R8), allocatable :: pert_rseed_r8(:,:)

contains

!BOP
! !IROUTINE: SetServices -- Sets ESMF services for component
! !INTERFACE:

subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),intent(INOUT) :: GC
    integer, optional,  intent(  OUT) :: RC

! !DESCRIPTION:

!EOP
!
! ErrLog Variables

    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: COMP_NAME
    integer                    :: STATUS

! Local Variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(ESMF_Config) :: CF
! Begin...
! --------

! Get my name and set-up traceback handle
! ------------------------------------------------------------------------------

    Iam='SetServices'
    call ESMF_GridCompGet ( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//trim(Iam)

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

  ! Register services for this component
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_INITIALIZE,                                                &
         Initialize,                                                            &
         rc=status                                                              &
         )
    VERIFY_(status)

    !phase 1: assimilation run
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         RUN,                                                                   &
         rc=status                                                              &
         )
    VERIFY_(status)

    !phase 2: feed back to change catch_progn  
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         UPDATE_ASSIM,                                                          &
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



! Set the state variable specs.
! -----------------------------
!BOS
!
! IMPORT STATE:
!
! ---------------------------------

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'soil_porosity'             ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'POROS'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'sfc_sat_hydraulic_conduct' ,&
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'COND'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'saturated_matric_potential',&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'PSIS'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'clapp_hornberger_b'        ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BEE'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'wetness_at_wilting_point'  ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WPWET'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'vertical_transmissivity'   ,&
    UNITS              = 'm-1'                       ,&
    SHORT_NAME         = 'GNU'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'max_rootzone_water_content',&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'VGWMAX'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'topo_baseflow_param_1'     ,&
    UNITS              = 'kg m-4'                    ,&
    SHORT_NAME         = 'BF1'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'topo_baseflow_param_2'     ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'BF2'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'topo_baseflow_param_3'     ,&
    UNITS              = 'log(m)'                    ,&
    SHORT_NAME         = 'BF3'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'moisture_threshold'        ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CDCR1'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'max_water_content'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CDCR2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'wetness_param_1'           ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARS1'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'wetness_param_2'           ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARS2'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'wetness_param_3'           ,&
    UNITS              = 'm+4 kg-2'                  ,&
    SHORT_NAME         = 'ARS3'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'shape_param_1'             ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARA1'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'shape_param_2'             ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ARA2'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'shape_param_3'             ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARA3'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'shape_param_4'             ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ARA4'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'min_theta_param_1'         ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARW1'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'min_theta_param_2'         ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARW2'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'min_theta_param_3'          ,&
    UNITS              = 'm+4 kg-2'                  ,&
    SHORT_NAME         = 'ARW3'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'min_theta_param_4'         ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ARW4'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_1'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TSA1'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_2'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TSA2'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_3'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TSB1'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_4'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TSB2'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_5'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ATAU'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_6'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BTAU'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC                     ,&
     SHORT_NAME = 'ITY'                              ,&
     LONG_NAME  = 'vegetation_type'                  ,&
     UNITS      = '1'                                ,&
     DIMS       = MAPL_DimsTileOnly                  ,&
     VLOCATION  = MAPL_VLocationNone                 ,&
     RC=STATUS  )
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC                     ,&
     SHORT_NAME = 'Z2CH'                             ,&
     LONG_NAME  = 'vegetation_height'                ,&
     UNITS      = 'm'                                ,&
     DIMS       = MAPL_DimsTileOnly                  ,&
     VLOCATION  = MAPL_VLocationNone                 ,&
     RC=STATUS  )
  VERIFY_(STATUS)

!
! Export for incr
!

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_canopy_temperature_saturated_zone'  ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TCFSAT_INCR'                        ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_canopy_temperature_transition_zone'   ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TCFTRN_INCR'                        ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_canopy_temperature_wilting_zone'    ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TCFWLT_INCR'                        ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_canopy_specific_humidity_saturated_zone'  ,&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QCFSAT_INCR'                        ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_canopy_specific_humidity_transition_zone'  ,&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QCFTRN_INCR'                        ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_canopy_specific_humidity_wilting_zone'  ,&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QCFWLT_INCR'                        ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_interception_reservoir_capac',&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CAPAC_INCR'                     ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_catchment_deficit'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CATDEF_INCR'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_root_zone_excess'          ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'RZEXC_INCR'                     ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_surface_excess'            ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'SRFEXC_INCR'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_soil_heat_content_layer_1' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT1_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_soil_heat_content_layer_2' ,&
    UNITS              = 'J_m-2'                     ,&
    SHORT_NAME         = 'GHTCNT2_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_soil_heat_content_layer_3' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT3_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_soil_heat_content_layer_4' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT4_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_soil_heat_content_layer_5' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT5_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_soil_heat_content_layer_6' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT6_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_snow_mass_layer_1'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN1_INCR'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_snow_mass_layer_2'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN2_INCR'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_snow_mass_layer_3'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN3_INCR'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_heat_content_snow_layer_1' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN1_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_heat_content_snow_layer_2' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN2_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_heat_content_snow_layer_3' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN3_INCR'                   ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_snow_depth_layer_1'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN1_INCR'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_snow_depth_layer_2'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN2_INCR'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC                  ,&
    LONG_NAME          = 'increment_snow_depth_layer_3'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN3_INCR'                    ,&
!    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

!
! INTERNAL STATE
!
   call MAPL_AddInternalSpec(GC                ,&
      LONG_NAME   = 'L-band Microwave RTM: Vegetation class. Type is Unsigned32'   ,&
      UNITS       = '1'                        ,&
      SHORT_NAME  = 'MWRTM_VEGCLS'             ,&
      DIMS        = MAPL_DimsTileOnly          ,&
      VLOCATION   = MAPL_VLocationNone         ,&
      FRIENDLYTO  = trim(COMP_NAME)             ,&
      DEFAULT     = MAPL_UNDEF                , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                  ,&
      LONG_NAME   = 'L-band Microwave RTM: Soil class. Type is Unsigned32'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_SOILCLS'             ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                  ,&
      LONG_NAME   = 'L-band Microwave RTM: Sand fraction'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_SAND'                ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                  ,&
      LONG_NAME   = 'L-band Microwave RTM: Clay fraction'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_CLAY'                ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: Porosity'   ,&
      UNITS       = 'm3 m-3'                    ,&
      SHORT_NAME  = 'MWRTM_POROS'               ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: Wang dielectric model transition soil moisture'   ,&
      UNITS       = 'm3 m-3'                    ,&
      SHORT_NAME  = 'MWRTM_WANGWT'              ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: Wang dielectric model wilting point soil moisture'   ,&
      UNITS       = 'm3 m-3'                    ,&
      SHORT_NAME  = 'MWRTM_WANGWP'              ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: Minimum microwave roughness parameter'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_RGHHMIN'             ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: Maximum microwave roughness parameter'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_RGHHMAX'             ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                  ,&
      LONG_NAME   = 'L-band Microwave RTM: Soil moisture value below which maximum microwave roughness parameter is used'   ,&
      UNITS       = 'm3 m-3'                    ,&
      SHORT_NAME  = 'MWRTM_RGHWMIN'             ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                 , &
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                  ,&
      LONG_NAME   = 'L-band Microwave RTM: Soil moisture value above which minimum microwave roughness parameter is used'   ,&
      UNITS       = 'm3 m-3'                    ,&
      SHORT_NAME  = 'MWRTM_RGHWMAX'             ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                  ,&
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: H-pol. Exponent for rough reflectivity parameterization'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_RGHNRH'              ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                  ,&
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: V-pol. Exponent for rough reflectivity parameterization'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_RGHNRV'              ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                  ,&
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: Polarization mixing parameter'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_RGHPOLMIX'           ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                  ,&
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: Scattering albedo'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_OMEGA'               ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                  ,&
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: H-pol. Vegetation b parameter'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_BH'                  ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                  ,&
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: V-pol. Vegetation b parameter'   ,&
      UNITS       = '1'                         ,&
      SHORT_NAME  = 'MWRTM_BV'                  ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                  ,&
      RC=STATUS)

   call MAPL_AddInternalSpec(GC                 ,&
      LONG_NAME   = 'L-band Microwave RTM: Parameter to transform leaf area index into vegetation water content'   ,&
      UNITS       = 'kg m-2'                    ,&
      SHORT_NAME  = 'MWRTM_LEWT'                ,&
      DIMS        = MAPL_DimsTileOnly           ,&
      VLOCATION   = MAPL_VLocationNone          ,&
      DEFAULT     = MAPL_UNDEF                  ,&
      RC=STATUS)

    call MAPL_TimerAdd(GC, name="Initialize"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GenericSetServices ( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROTUINE: Initialize -- initialize method for LandAssim GC

! !INTERFACE:
subroutine Initialize(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

  ! ESMF variables
    type(ESMF_Time)  :: CurrentTime
    type(ESMF_Alarm) :: LandAssimAlarm
    type(ESMF_TimeInterval) :: LandAssim_DT
    integer :: LandAssimDTstep
    type(ESMF_TimeInterval) :: ModelTimeStep

    ! locals
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(MAPL_LocStream) :: locstream

    character(len=300) :: out_path,fname
    character(len=ESMF_MAXSTR) :: exp_id, GridName
    integer :: model_dtstep
    type(date_time_type) :: start_time

   ! LDAS' tile_coord variable
    type(T_TILECOORD_STATE), pointer :: tcinternal
    type(TILECOORD_WRAP) :: tcwrap

    type(tile_coord_type), dimension(:), pointer :: tile_coord_f => null()
    type(tile_coord_type), dimension(:), pointer :: tile_coord_l => null()

    integer :: land_nt_local,i,mpierr, ens
    ! mapping f to re-orderd f so it is continous for mpi_gather
    ! rf -- ordered by processors. Within the processor, ordered by MAPL grid
    integer,allocatable :: f2rf(:) ! mapping re-orderd rf to f for the LDASsa output
    type(grid_def_type) :: tile_grid_g
    type(grid_def_type) :: tile_grid_f
    character(len=300) :: seed_fname
    character(len=300) :: fname_tpl
    character(len=14) :: datestamp
    character(len=4)  :: id_string
    integer :: nymd, nhms
!! from LDASsa

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
    
    call MAPL_GetResource ( MAPL, out_path, Label="OUT_PATH:", DEFAULT="./", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, exp_id, Label="EXP_ID:", DEFAULT="exp_id", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, NUM_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, FIRST_ENS_ID, Label="FIRST_ENS_ID:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    call init_log( myid, numprocs, master_proc )

    ! Get current time
    call ESMF_ClockGet(clock, currTime=CurrentTime, rc=status)
    VERIFY_(status)
    call esmf2ldas(CurrentTime, start_time, rc=status)
    VERIFY_(status)

    call ESMF_ClockGet(clock, timeStep=ModelTimeStep,rc=status)
    VERIFY_(status)
    call ESMF_TimeIntervalGet(ModelTimeStep, s=model_dtstep,rc=status)
    VERIFY_(status)

    ! Create alarm for Land assimilation
    ! -create-nonsticky-alarm-
       ! -time-interval-
    call MAPL_GetResource(                                                      &
         MAPL,                                                                  &
         LandAssimDtStep,                                                           &
         'LANDASSIM_DTSTEP:',                                                       &
         default=10800,                                                          &
         rc=status                                                              &
         )
    VERIFY_(status)

    call ESMF_TimeIntervalSet(LandAssim_DT, s=LandAssimDtStep, rc=status)
    VERIFY_(status)

    LandAssimAlarm = ESMF_AlarmCreate(                                          &
         clock,                                                                 &
         name='LandAssim',                                                      &
         ringTime=CurrentTime+Landassim_DT-ModelTimeStep,                       &
         ringInterval=Landassim_DT,                                             &
         ringTimeStepCount=1,                                                   &
         sticky=.false.,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)

    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tcinternal   =>tcwrap%ptr
    tile_coord_l =>tcinternal%tile_coord
  ! Get number of land tiles
    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)
    call MAPL_LocStreamGet(locstream, NT_LOCAL=land_nt_local,rc=status)
    VERIFY_(status)

    allocate(Pert_rseed(NRANDSEED, NUM_ENSEMBLE), source = 0)
    allocate(Pert_rseed_r8(NRANDSEED, NUM_ENSEMBLE), source = 0.0d0)

    if (master_proc) then

       call MAPL_GetResource ( MAPL, fname_tpl, Label="LANDASSIM_OBSPERTRSEED_RESTART_FILE:", DEFAULT="../intput/restart/landassim_obspertrseed%s_rst", RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_DateStampGet( clock, datestamp, rc=status)
       VERIFY_(STATUS)
       read(datestamp(1:8),*)   nymd
       read(datestamp(10:13),*) nhms
       nhms = nhms*100
       do ens = 0, NUM_ENSEMBLE-1
         write(id_string,'(I4.4)') ens + FIRST_ENS_ID
         seed_fname = ""
         call ESMF_CFIOStrTemplate(seed_fname,fname_tpl,'GRADS', xid=id_string,nymd=nymd,nhms=nhms,stat=status)
         call read_pert_rseed(seed_fname,Pert_rseed_r8(:,ens+1))
  
         Pert_rseed(:,ens+1) = nint(Pert_rseed_r8(:,ens+1))
         if (all(Pert_rseed(:,ens+1) == 0)) then
            call get_init_pert_rseed(ens, pert_rseed(1,ens+1))
            call init_randseed(pert_rseed(:,ens+1))
         endif
       enddo
    endif
    call MPI_Bcast(pert_rseed, NRANDSEED*NUM_ENSEMBLE, MPI_INTEGER, 0, mpicomm, mpierr)


    allocate(N_catl_vec(numprocs))
    allocate(low_ind(numprocs))
    allocate(l2rf(land_nt_local))

    call MPI_AllGATHER(land_nt_local,1,MPI_INTEGER,N_catl_vec,1,MPI_INTEGER,mpicomm,mpierr)

    low_ind(1) = 1
    do i = 2, numprocs
       low_ind(i) = low_ind(i-1) + N_catl_vec(i-1)
    enddo
    N_catf = sum(N_catl_vec)
    allocate(rf2f(N_catf))
    allocate(f2rf(N_catf))
    
    call MPI_AllGATHERV(tcinternal%l2f,  land_nt_local, MPI_INTEGER,              &
                     rf2f, N_catl_vec,    low_ind-1,  MPI_INTEGER,  &
                     mpicomm,mpierr)

    allocate(tile_coord_rf(N_catf))
    tile_coord_rf(:) = tcwrap%ptr%tile_coord_f(rf2f(:))
    allocate(rf2g(N_catf))
    rf2g(:) = tile_coord_rf(:)%tile_id

    do i=1,N_catf
      f2rf(rf2f(i))= i
      tile_coord_rf(i)%f_num = i
    enddo

    do i=1, land_nt_local
       l2rf(i) = low_ind(myid+1) + i - 1
    end do

    tcwrap%ptr%tile_coord%f_num = l2rf

  ! invert mapping from local to full grid (get f2l from l2f)

    allocate(rf2l(N_catf))

    rf2l = -9999

    do i=1,land_nt_local
      rf2l( l2rf(i) ) = i
    end do
    
    if (master_proc) then
      call read_ens_upd_inputs(                 &
       trim(out_path),                          &
       trim(exp_id),                            &
       start_time,                              &
       model_dtstep,                            &
       N_catf,             tile_coord_rf,       &
       N_progn_pert,       progn_pert_param,    &
       N_force_pert,       force_pert_param,    &
       need_mwRTM_param,                        &
       update_type,                             &
       dtstep_assim,                            &
       centered_update,                         &
       xcompact, ycompact,                      &
       N_obs_param,                             &
       obs_param,                               &
       out_obslog,                              &
       out_ObsFcstAna,                          &
!       out_incr,                                &
!       out_incr_format,                         &
       out_smapL4SMaup,                         &
       N_obsbias_max                            &
       )
      call MAPL_GetResource ( MAPL, GridName, Label="GEOSldas.GRIDNAME:", DEFAULT="EASE", RC=STATUS)
      VERIFY_(STATUS)
      if (index(GridName,"-CF") /=0) out_smapL4SMaup = .false. ! no out_smap for now if it is cs frid
    endif

    call MPI_BCAST(need_mwRTM_param,      1, MPI_LOGICAL,        0,MPICOMM,mpierr)
    call MPI_BCAST(update_type,           1, MPI_INTEGER,        0,MPICOMM,mpierr)
    call MPI_BCAST(dtstep_assim,          1, MPI_INTEGER,        0,MPICOMM,mpierr)
    call MPI_BCAST(centered_update,       1, MPI_LOGICAL,        0,MPICOMM,mpierr)
    call MPI_BCAST(xcompact,              1, MPI_REAL,           0,MPICOMM,mpierr)
    call MPI_BCAST(ycompact,              1, MPI_REAL,           0,MPICOMM,mpierr)
    call MPI_BCAST(N_obs_param,           1, MPI_INTEGER,        0,MPICOMM,mpierr)
    call MPI_BCAST(out_obslog,            1, MPI_LOGICAL,        0,MPICOMM,mpierr)
    call MPI_BCAST(out_ObsFcstAna,        1, MPI_LOGICAL,        0,MPICOMM,mpierr)
!    call MPI_BCAST(out_incr,              1, MPI_LOGICAL,        0,MPICOMM,mpierr)
!    call MPI_BCAST(out_incr_format,       1, MPI_INTEGER,        0,MPICOMM,mpierr)
    call MPI_BCAST(out_smapL4SMaup,       1, MPI_LOGICAL,        0,MPICOMM,mpierr)
    call MPI_BCAST(N_obsbias_max,         1, MPI_INTEGER,        0,MPICOMM,mpierr)

    if (.not. master_proc)  allocate(obs_param(N_obs_param))

    call MPI_BCAST(obs_param,   N_obs_param, MPI_OBS_PARAM_TYPE, 0,MPICOMM,mpierr)

    if (master_proc) call echo_clsm_ensupd_glob_param(logunit) 

    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    call MAPL_TimerOff(MAPL, "Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")

    RETURN_(ESMF_SUCCESS)

end subroutine Initialize

! !IROUTINE: RUN 
! !INTERFACE:
subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

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
!  time
!  
    type(ESMF_Time)  :: ModelTimeCur, ModelTimeNxt
    type(ESMF_Alarm) :: LandAssimAlarm
    type(ESMF_TimeInterval) :: ModelTimeStep


! Locals
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(TILECOORD_WRAP)           :: tcwrap
    type(tile_coord_type), pointer :: tile_coord_l(:)=>null()
    type(T_TILECOORD_STATE), pointer :: tcinternal
    type(mwRTM_param_type),dimension(:),allocatable ::  mwRTM_param

    type(ESMF_State) :: INTERNAL
    type(date_time_type) :: start_time
    type(date_time_type) :: date_time_new
    character(len=14)    :: datestamp

    integer :: N_catl, N_catg,N_obsl_max, n_e, i

    character(len=300) :: out_path
    character(len=ESMF_MAXSTR) :: exp_id
    character(40)  ::  exp_domain
    integer :: model_dtstep
    type(met_force_type), dimension(:), allocatable :: met_force

    integer :: N_adapt_R
    type(MAPL_LocStream) :: locstream
    integer,dimension(:),allocatable :: obs_pert_adapt_param
    real,dimension(:,:), allocatable :: Pert_adapt_R
    real,dimension(:,:), allocatable :: Obs_pert
    type(obs_bias_type),       dimension(:,:,:), allocatable :: obs_bias 

    type(cat_progn_type), dimension(:,:), allocatable :: cat_progn_incr
    type(cat_progn_type), dimension(:), allocatable :: cat_progn_incr_ensavg
    type(obs_type),       dimension(:), pointer :: Observations_l => null()
    logical  :: fresh_incr
    integer  :: N_obsf,N_obsl
    integer  :: secs_in_day
! -----------------------------------------------------
! INTERNAL Pointers
! -----------------------------------------------------

    real, dimension(:), pointer :: VEGCLS
    real, dimension(:), pointer :: SOILCLS
    real, dimension(:), pointer :: SAND
    real, dimension(:), pointer :: CLAY
    real, dimension(:), pointer :: mw_POROS
    real, dimension(:), pointer :: WANGWT
    real, dimension(:), pointer :: WANGWP
    real, dimension(:), pointer :: RGHHMIN
    real, dimension(:), pointer :: RGHHMAX
    real, dimension(:), pointer :: RGHWMAX
    real, dimension(:), pointer :: RGHWMIN
    real, dimension(:), pointer :: RGHNRH
    real, dimension(:), pointer :: RGHNRV
    real, dimension(:), pointer :: RGHPOLMIX
    real, dimension(:), pointer :: OMEGA
    real, dimension(:), pointer :: BH
    real, dimension(:), pointer :: BV
    real, dimension(:), pointer :: LEWT

!! import ensemble forcing

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
    real, pointer :: SWLAND(:)=>null()
    real, pointer :: LAI(:)=>null()

!! export incr progn

    real, dimension(:),pointer :: TC1_incr=>null()
    real, dimension(:),pointer :: TC2_incr=>null()
    real, dimension(:),pointer :: TC4_incr=>null()
    real, dimension(:),pointer :: QC1_incr=>null()
    real, dimension(:),pointer :: QC2_incr=>null()
    real, dimension(:),pointer :: QC4_incr=>null()
    real, dimension(:),pointer :: CAPAC_incr=>null()
    real, dimension(:),pointer :: CATDEF_incr=>null()
    real, dimension(:),pointer :: RZEXC_incr=>null()
    real, dimension(:),pointer :: SRFEXC_incr=>null()
    real, dimension(:),pointer :: GHTCNT1_incr=>null()
    real, dimension(:),pointer :: GHTCNT2_incr=>null()
    real, dimension(:),pointer :: GHTCNT3_incr=>null()
    real, dimension(:),pointer :: GHTCNT4_incr=>null()
    real, dimension(:),pointer :: GHTCNT5_incr=>null()
    real, dimension(:),pointer :: GHTCNT6_incr=>null()
    real, dimension(:),pointer :: WESNN1_incr=>null()
    real, dimension(:),pointer :: WESNN2_incr=>null()
    real, dimension(:),pointer :: WESNN3_incr=>null()
    real, dimension(:),pointer :: HTSNNN1_incr=>null()
    real, dimension(:),pointer :: HTSNNN2_incr=>null()
    real, dimension(:),pointer :: HTSNNN3_incr=>null()
    real, dimension(:),pointer :: SNDZN1_incr=>null()
    real, dimension(:),pointer :: SNDZN2_incr=>null()
    real, dimension(:),pointer :: SNDZN3_incr=>null()


    logical :: spin
    logical, save              :: firsttime=.true.
    type(cat_bias_param_type)  :: cat_bias_param
    integer                    :: N_catbias
    character(len=300) :: seed_fname
    character(len=300) :: fname_tpl
    character(len=4) :: id_string
    integer:: ens, nymd, nhms

#ifdef DBG_LANDASSIM_INPUTS
        ! vars for debugging purposes
        type(ESMF_Grid)                 :: TILEGRID
        integer, pointer                :: mask(:)
        integer                         :: nt,ens_id
        integer, save                   :: unit_i=0
        integer                         :: unit
        integer                         :: NT_GLOBAL,mpierr,i
        real,allocatable                :: metTair(:),metTair_l(:)
        integer,allocatable             :: ids(:)
#endif


    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//"::RUN"

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)
! Start timers
! ------------
    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN")
    call ESMF_ClockGetAlarm(clock, 'LandAssim', LandAssimAlarm, rc=status)
    VERIFY_(status)

    call MAPL_GetResource ( MAPL, out_path, Label="OUT_PATH:", DEFAULT="./", RC=STATUS)
    call MAPL_GetResource ( MAPL, exp_id, Label="EXP_ID:", DEFAULT="exp_id", RC=STATUS)

   ! Get component's internal variable
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tcinternal => tcwrap%ptr
    tile_coord_l => tcwrap%ptr%tile_coord

    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, rc=status)
    VERIFY_(status)

    ! Get current time
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep,rc=status)
    VERIFY_(status)
    call ESMF_ClockGet(clock, currTime=ModelTimeCur, rc=status)
    VERIFY_(status)
    call esmf2ldas(ModelTimeCur+ModelTimeStep, date_time_new, rc=status)
    VERIFY_(status)

    call esmf2ldas(ModelTimeCur, start_time, rc=status)
    VERIFY_(status)

    ! Get number of land tiles
    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)
    call MAPL_LocStreamGet(locstream, NT_LOCAL=N_catl,rc=status)
    VERIFY_(status)

! Pointers to internals
!----------------------
    if (need_mwRTM_param) then
       call MAPL_GetPointer(INTERNAL, SAND     , 'MWRTM_SAND'     ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL,SOILCLS   , 'MWRTM_SOILCLS'  ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, VEGCLS   , 'MWRTM_VEGCLS'   ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, CLAY     , 'MWRTM_CLAY'     ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, mw_POROS , 'MWRTM_POROS' ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, WANGWT   , 'MWRTM_WANGWT'   ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, WANGWP   , 'MWRTM_WANGWP'   ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, RGHHMIN  , 'MWRTM_RGHHMIN'  ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, RGHHMAX  , 'MWRTM_RGHHMAX'  ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, RGHWMIN  , 'MWRTM_RGHWMIN'  ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, RGHWMAX  , 'MWRTM_RGHWMAX'  ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, RGHNRH   , 'MWRTM_RGHNRH'   ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, RGHNRV   , 'MWRTM_RGHNRV'   ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, RGHPOLMIX, 'MWRTM_RGHPOLMIX',    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, OMEGA  , 'MWRTM_OMEGA'      ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, BH     , 'MWRTM_BH'         ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, BV     , 'MWRTM_BV'         ,    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, LEWT   , 'MWRTM_LEWT'       ,    RC=STATUS)
       VERIFY_(STATUS)
 
       allocate(mwRTM_param(N_catl))
       mwRTM_param(:)%sand      = SAND(:)
       mwRTM_param(:)%vegcls    = nint(VEGCLS(:))
       mwRTM_param(:)%soilcls   = nint(SOILCLS(:))
       mwRTM_param(:)%clay      = CLAY(:)
       mwRTM_param(:)%poros     = mw_POROS(:)
       mwRTM_param(:)%wang_wt   = WANGWT(:)
       mwRTM_param(:)%wang_wp   = WANGWP(:)
       mwRTM_param(:)%rgh_hmin  = RGHHMIN(:)
       mwRTM_param(:)%rgh_hmax  = RGHHMAX(:)
       mwRTM_param(:)%rgh_wmin  = RGHWMIN(:)
       mwRTM_param(:)%rgh_wmax  = RGHWMAX(:)
       mwRTM_param(:)%rgh_Nrh   = RGHNRH(:)
       mwRTM_param(:)%rgh_Nrv   = RGHNRV(:)
       mwRTM_param(:)%rgh_polmix= RGHPOLMIX(:)
       mwRTM_param(:)%omega     = OMEGA(:)
       mwRTM_param(:)%bh        = BH(:)
       mwRTM_param(:)%bv        = bv(:)
       mwRTM_param(:)%lewt      = LEWT(:)
    endif

    if (firsttime) then
       firsttime = .false.
       if (need_mwRTM_param) &
          call GEOS_output_smapL4SMlmc( GC, start_time, trim(out_path), trim(exp_id), &
           N_catl, tile_coord_l, cat_param, mwRTM_param )
       if (master_proc) then 
           ! for out put
          call read_cat_bias_inputs(  trim(out_path), trim(exp_id), start_time, update_type, &
          cat_bias_param, N_catbias)
       endif
    endif
 
    ! The time is one model time step behind Current time, so record the checkpoint here 
    if (MAPL_RecordAlarmIsRinging(MAPL)) then
       if (master_proc) then
          Pert_rseed_r8 = Pert_rseed
          call MAPL_GetResource ( MAPL, fname_tpl, Label="LANDASSIM_OBSPERTRSEED_CHECKPOINT_FILE:", DEFAULT="landassim_obspertrseed%s_checkpoint", RC=STATUS)
          VERIFY_(STATUS)
          fname_tpl = trim(fname_tpl) //".%y4%m2%d2_%h2%n2z.nc4"
          call MAPL_DateStampGet( clock, datestamp, rc=status)
          VERIFY_(STATUS)
          read(datestamp(1:8),*)   nymd
          read(datestamp(10:13),*) nhms
          nhms = nhms*100
          do ens = 0, NUM_ENSEMBLE-1
             write(id_string,'(I4.4)') ens + FIRST_ENS_ID
             seed_fname = ""
             call ESMF_CFIOStrTemplate(seed_fname,fname_tpl,'GRADS', xid=id_string,nymd=nymd,nhms=nhms,stat=status)
             VERIFY_(STATUS)
             call write_pert_rseed(trim(seed_fname), Pert_rseed_r8(:,ens+1))
          enddo
        endif 
    endif


    if ( .not. ESMF_AlarmIsRinging(LandAssimAlarm)) then
       call MAPL_TimerOff ( MAPL, "RUN"  )
       call MAPL_TimerOff ( MAPL, "TOTAL" )
       RETURN_(ESMF_SUCCESS)
    endif

    N_obsl_max = N_catl*N_obs_param
   
!! get import from ens to get ensemble average forcing

    call MAPL_GetPointer(import, TA_enavg, 'TA', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, QA_enavg, 'QA', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PS_enavg, 'PS', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, UU_enavg, 'UU', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PCU_enavg, 'PCU', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PLS_enavg, 'PLS', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SNO_enavg, 'SNO', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRPAR_enavg, 'DRPAR', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFPAR_enavg, 'DFPAR', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRNIR_enavg, 'DRNIR', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFNIR_enavg, 'DFNIR', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRUVR_enavg, 'DRUVR', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFUVR_enavg, 'DFUVR', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LWDNSRF_enavg, 'LWDNSRF', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DZ_enavg, 'DZ', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SWLAND, 'SWLAND', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LAI, 'LAI', rc=status)
    VERIFY_(status)
!
! export for incr
!
    call MAPL_GetPointer(export, TC1_incr,  'TCFSAT_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TC2_incr,  'TCFTRN_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, TC4_incr,  'TCFWLT_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QC1_incr,  'QCFSAT_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QC2_incr,  'QCFTRN_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QC4_incr,  'QCFWLT_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CAPAC_incr,  'CAPAC_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CATDEF_incr,  'CATDEF_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RZEXC_incr,  'RZEXC_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SRFEXC_incr,  'SRFEXC_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT1_incr,  'GHTCNT1_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT2_incr,  'GHTCNT2_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT3_incr,  'GHTCNT3_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT4_incr,  'GHTCNT4_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT5_incr,  'GHTCNT5_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT6_incr,  'GHTCNT6_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN1_incr,  'WESNN1_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN2_incr,  'WESNN2_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN3_incr,  'WESNN3_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN1_incr,  'HTSNNN1_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN2_incr,  'HTSNNN2_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN3_incr,  'HTSNNN3_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN1_incr,  'SNDZN1_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN2_incr,  'SNDZN2_INCR' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN3_incr,  'SNDZN3_INCR' ,rc=status)
    VERIFY_(status)

    allocate(met_force(N_catl))
    met_force(:)%Tair    = TA_enavg(:)
    met_force(:)%Qair    = QA_enavg(:)    
    met_force(:)%Psurf   = PS_enavg(:)    
    met_force(:)%Rainf_c = PCU_enavg(:)    
    met_force(:)%Rainf   = PCU_enavg(:) + PLS_enavg(:)   
    met_force(:)%Snowf   = SNO_enavg(:)
    met_force(:)%LWdown  = LWDNSRF_enavg(:)
    met_force(:)%SWdown  = DRPAR_enavg(:)+DFPAR_enavg(:)+DRNIR_enavg(:) + &
                           DFNIR_enavg(:)+DRUVR_enavg(:)+DFUVR_enavg(:)
    met_force(:)%SWnet   = SWLAND(:)
    met_force(:)%PARdrct = DRPAR_enavg(:)
    met_force(:)%PARdffs = DFPAR_enavg(:)
    met_force(:)%wind    = UU_enavg(:)
    met_force(:)%RefH    = DZ_enavg(:)

!   Weiyuan note: dummy adapt for now
    N_adapt_R = 0
    !allocate(obs_pert_adapt_param(N_obs_param))
    !allocate(Pert_adapt_R(N_adapt_R,NUM_ENSEMBLE))
    !allocate(Obs_pert(N_obsl_max,NUM_ENSEMBLE))

    if (N_obsbias_max>0) then
        allocate(obs_bias(N_catl,N_obs_param,N_obsbias_max))
        call initialize_obs_bias( N_catf, N_obs_param, N_obsbias_max, trim(out_path), &
                      trim(exp_id), start_time, N_catl, numprocs, N_catl_vec, low_ind, obs_bias)
    end if

    allocate(cat_progn_incr(N_catl,NUM_ENSEMBLE))
    allocate(cat_progn_incr_ensavg(N_catl))
    allocate(Observations_l(N_obsl_max))

    !WY note: temportary


#ifdef DBG_LANDASSIM_INPUTS

    if (firsttime) then
        firsttime = .false.
        call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS)
        VERIFY_(STATUS)

        call MAPL_TileMaskGet(tilegrid,  mask, rc=status)
        VERIFY_(STATUS)

        allocate(metTair(N_catf),metTair_l(N_catl))
        allocate(ids(N_catf))

        metTair_l(:) = met_force(:)%Tair
        ids(:) = tile_coord_rf(:)%tile_id

        call MPI_AllGATHERV(metTair_l,         N_catl,     MPI_REAL,  &
                        metTair,   N_catl_vec, low_ind-1,  MPI_REAL,  &
                        mpicomm,mpierr)
        

        if(myid ==0) then
          open(unit=10,file='metTair.txt',action="write",status="replace")
          do i = 1, N_catf
            write(10,*) ids(i), metTair(i)
          enddo
          close(10)
        endif

        unit = GETFILE( "landassim_force_inputs.bin", form="unformatted", RC=STATUS )
        VERIFY_(STATUS)
! Inputs
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%Tair, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%Qair, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%Psurf, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%Rainf_c, mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, met_force(:)%Rainf,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%Snowf,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%LWdown, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%SWdown, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%SWnet,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%PARdrct,mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%PARdffs, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%wind, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, met_force(:)%RefH, mask=mask, rc=status); VERIFY_(STATUS)

        unit = GETFILE( "landassim_catprogn_inputs.bin", form="unformatted", RC=STATUS )
        VERIFY_(STATUS)
       
        ens_id = 1

        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%tc1,     mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%tc2,     mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%tc4,     mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%qa1,     mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%qa2,     mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%qa4,     mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%capac,   mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%catdef,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%rzexc,   mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%srfexc,  mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%ght(1),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%ght(2),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%ght(3),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%ght(4),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%ght(5),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%ght(6),  mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%wesn(1), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%wesn(2), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%wesn(3), mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%htsn(1), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%htsn(2), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%htsn(3), mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%sndz(1), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%sndz(2), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, cat_progn(:,ens_id)%sndz(3), mask=mask, rc=status); VERIFY_(STATUS)


        unit = GETFILE( "landassim_mwrtm_inputs.bin", form="unformatted", RC=STATUS )
        VERIFY_(STATUS)
       
        call MAPL_VarWrite(unit, tilegrid,real(mwRTM_param(:)%vegcls),    mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,real(mwRTM_param(:)%soilcls),   mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%sand,      mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%clay,      mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%poros,     mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%wang_wt,   mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%wang_wp,   mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_hmin,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_hmax,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_wmin,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_wmax,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_Nrh,   mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_Nrv,   mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_polmix,mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%omega,     mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%bh,        mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%bv,        mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%lewt,      mask=mask, rc=status); VERIFY_(STATUS)

        !unit = GETFILE( "landassim_catparam_inputs.bin", form="unformatted", RC=STATUS )
        !VERIFY_(STATUS)
       
     endif



#endif

       call get_enkf_increments(                                            &
          date_time_new,                                                    &
          NUM_ENSEMBLE, N_catl, N_catf, N_obsl_max,                         &
          trim(out_path), trim(exp_id), exp_domain,                         &
          met_force, lai, cat_param, mwRTM_param,                           &
          tile_coord_l, tile_coord_rf, tcinternal%grid_f,                   &
          tcinternal%grid_f, tcinternal%grid_l, tcinternal%grid_g,          &
          N_catl_vec, low_ind, l2rf, rf2l,                                  &
          N_force_pert, N_progn_pert, force_pert_param, progn_pert_param,   &
          update_type,                                                      &
          dtstep_assim, centered_update,                                    &
          xcompact, ycompact,                                               &
          N_obs_param, obs_param, N_obsbias_max,                            &
          out_obslog, out_smapL4SMaup,                                      &
          cat_progn,                                                        &
          Pert_rseed, obs_bias,                                             &
          cat_progn_incr, fresh_incr,                                       &
          N_obsf, N_obsl, Observations_l,                                   &
        ! below  are dummy for now
          N_adapt_R, obs_pert_adapt_param, Pert_adapt_R)
          !Obs_pert )

        ! forced to apply
       spin = .false.

       if (.not. spin) then
          if (fresh_incr) then
          ! apply EnKF increments 
          ! (without call to subroutine recompute_diagnostics())
             call apply_enkf_increments( N_catl, NUM_ENSEMBLE, update_type, cat_param, &
                   cat_progn_incr, cat_progn )

           end if ! fresh_incr

           ! if requested, write incr and/or ObsFcstAna files whenever it was
           ! time for assimilation, even if there were no observations
           ! - reichle, 29 Aug 2014           

           secs_in_day = &
                date_time_new%hour*3600 + date_time_new%min*60 + date_time_new%sec

           if (centered_update)  secs_in_day = secs_in_day + dtstep_assim/2
       
           ! WY note : Here N_catg is not the global land tile number
           ! but a maximum global_id this simulation covers. 
           ! Need to find the number 
           N_catg = maxval(rf2g)

           if (mod(secs_in_day, dtstep_assim)==0) then

               call output_incr_etc( out_ObsFcstAna,             &
                date_time_new, trim(out_path), trim(exp_id),            &
                N_obsl, N_obs_param, NUM_ENSEMBLE,                           &
                N_catl, tile_coord_l,                                        &
                N_catf, tile_coord_rf, tcinternal%grid_f, tcinternal%grid_g, &
                N_catl_vec, low_ind, rf2l, N_catg, rf2g,                     &
                obs_param,                                                   &
                met_force, lai,                     &
                cat_param, cat_progn, cat_progn_incr, mwRTM_param,           &
                Observations_l, rf2f=rf2f )


              do i = 1, N_catl
                 cat_progn_incr_ensavg(i) = 0.0
                 do n_e=1, NUM_ENSEMBLE
                    cat_progn_incr_ensavg(i) = cat_progn_incr_ensavg(i) &
                       + cat_progn_incr(i,n_e)
                 end do
                 cat_progn_incr_ensavg(i) = cat_progn_incr_ensavg(i)/real(NUM_ENSEMBLE)
              enddo

              if(associated(TC1_incr))  TC1_incr(:) = cat_progn_incr_ensavg(:)%tc1
              if(associated(TC2_incr))  TC2_incr(:) = cat_progn_incr_ensavg(:)%tc2
              if(associated(TC4_incr))  TC4_incr(:) = cat_progn_incr_ensavg(:)%tc4
              if(associated(QC1_incr))  QC1_incr(:) = cat_progn_incr_ensavg(:)%qa1
              if(associated(QC2_incr))  QC2_incr(:) = cat_progn_incr_ensavg(:)%qa2
              if(associated(QC4_incr))  QC4_incr(:) = cat_progn_incr_ensavg(:)%qa4

              if(associated(CAPAC_incr))   CAPAC_incr(:)  = cat_progn_incr_ensavg(:)%capac
              if(associated(CATDEF_incr))  CATDEF_incr(:) = cat_progn_incr_ensavg(:)%catdef
              if(associated(RZEXC_incr))   RZEXC_incr(:)  = cat_progn_incr_ensavg(:)%rzexc
              if(associated(SRFEXC_incr))  SRFEXC_incr(:) = cat_progn_incr_ensavg(:)%srfexc

              if(associated(GHTCNT1_incr)) GHTCNT1_incr(:) = cat_progn_incr_ensavg(:)%ght(1)
              if(associated(GHTCNT2_incr)) GHTCNT2_incr(:) = cat_progn_incr_ensavg(:)%ght(2)
              if(associated(GHTCNT3_incr)) GHTCNT3_incr(:) = cat_progn_incr_ensavg(:)%ght(3)
              if(associated(GHTCNT4_incr)) GHTCNT4_incr(:) = cat_progn_incr_ensavg(:)%ght(4)
              if(associated(GHTCNT5_incr)) GHTCNT5_incr(:) = cat_progn_incr_ensavg(:)%ght(5)
              if(associated(GHTCNT6_incr)) GHTCNT6_incr(:) = cat_progn_incr_ensavg(:)%ght(6)

              if(associated(WESNN1_incr))  WESNN1_incr(:) = cat_progn_incr_ensavg(:)%wesn(1)
              if(associated(WESNN2_incr))  WESNN2_incr(:) = cat_progn_incr_ensavg(:)%wesn(2)
              if(associated(WESNN3_incr))  WESNN3_incr(:) = cat_progn_incr_ensavg(:)%wesn(3)

              if(associated(HTSNNN1_incr)) HTSNNN1_incr(:) = cat_progn_incr_ensavg(:)%htsn(1)
              if(associated(HTSNNN2_incr)) HTSNNN2_incr(:) = cat_progn_incr_ensavg(:)%htsn(2)
              if(associated(HTSNNN3_incr)) HTSNNN3_incr(:) = cat_progn_incr_ensavg(:)%htsn(3)

              if(associated(SNDZN1_incr))  SNDZN1_incr(:)  = cat_progn_incr_ensavg(:)%sndz(1)
              if(associated(SNDZN2_incr))  SNDZN2_incr(:)  = cat_progn_incr_ensavg(:)%sndz(2)
              if(associated(SNDZN3_incr))  SNDZN3_incr(:)  = cat_progn_incr_ensavg(:)%sndz(3)



           ! write analysis fields into SMAP L4_SM aup file 
           ! whenever it was time for assimilation (regardless 
           ! of whether obs were actually assimilated and fresh
           ! increments were computed)

              if (out_smapL4SMaup)                                                  &
                   call write_smapL4SMaup( 'analysis', date_time_new, trim(out_path),     &
                   trim(exp_id), NUM_ENSEMBLE, N_catl, N_catf, N_obsl, tile_coord_rf,     &
                   tcinternal%grid_g, N_catl_vec, low_ind,                          &
                   N_obs_param, obs_param, Observations_l, cat_param, cat_progn   )

           end if


          fresh_incr = .false.

     endif !spin



!--------------------
! Pointers to inputs
!--------------------
    deallocate(cat_progn_incr)
    deallocate(cat_progn_incr_ensavg)
    deallocate(Observations_l)

    call MAPL_TimerOff ( MAPL, "RUN"  )
    call MAPL_TimerOff ( MAPL, "TOTAL" )

    RETURN_(ESMF_SUCCESS)

end subroutine RUN

 ! !IROTUINE: collecting and averaging

subroutine UPDATE_ASSIM(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    ! this export is from land grid come
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam='UPDATE_ASSIM'
    character(len=ESMF_MAXSTR) :: comp_name

    ! ESMF variables
    type(ESMF_Alarm) :: LandAssimAlarm
    type(ESMF_VM) :: vm

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null() ! MAPL obj

    real, dimension(:,:),pointer :: TC
    real, dimension(:,:),pointer :: QC
    real, dimension(:),pointer :: CAPAC
    real, dimension(:),pointer :: CATDEF
    real, dimension(:),pointer :: RZEXC
    real, dimension(:),pointer :: SRFEXC
    real, dimension(:),pointer :: GHTCNT1
    real, dimension(:),pointer :: GHTCNT2
    real, dimension(:),pointer :: GHTCNT3
    real, dimension(:),pointer :: GHTCNT4
    real, dimension(:),pointer :: GHTCNT5
    real, dimension(:),pointer :: GHTCNT6
    real, dimension(:),pointer :: WESNN1
    real, dimension(:),pointer :: WESNN2
    real, dimension(:),pointer :: WESNN3
    real, dimension(:),pointer :: HTSNNN1
    real, dimension(:),pointer :: HTSNNN2
    real, dimension(:),pointer :: HTSNNN3
    real, dimension(:),pointer :: SNDZN1
    real, dimension(:),pointer :: SNDZN2
    real, dimension(:),pointer :: SNDZN3

    integer,save :: ens_id = 0

  !BOP

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//"::RUN"

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_ClockGetAlarm(clock, 'LandAssim', LandAssimAlarm, rc=status)
    VERIFY_(status)
    if ( .not. ESMF_AlarmIsRinging(LandAssimAlarm)) then
       RETURN_(ESMF_SUCCESS)
    endif

    call MAPL_GetPointer(export, TC,  'TC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QC,  'QC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CAPAC,  'CAPAC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, CATDEF,  'CATDEF' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RZEXC,  'RZEXC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SRFEXC,  'SRFEXC' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT1,  'GHTCNT1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT2,  'GHTCNT2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT3,  'GHTCNT3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT4,  'GHTCNT4' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT5,  'GHTCNT5' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, GHTCNT6,  'GHTCNT6' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN1,  'WESNN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN2,  'WESNN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, WESNN3,  'WESNN3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN1,  'HTSNNN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN2,  'HTSNNN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, HTSNNN3,  'HTSNNN3' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN1,  'SNDZN1' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN2,  'SNDZN2' ,rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNDZN3,  'SNDZN3' ,rc=status)
    VERIFY_(status)

    ! This counter is relative to ens_id
    ens_id = ens_id + 1
    !distrbute catch_progn

    TC(:,1)     = cat_progn(:,ens_id)%tc1
    TC(:,2)     = cat_progn(:,ens_id)%tc2
    TC(:,3)     = cat_progn(:,ens_id)%tc4

    QC(:,1)     = cat_progn(:,ens_id)%qa1
    QC(:,2)     = cat_progn(:,ens_id)%qa2
    QC(:,3)     = cat_progn(:,ens_id)%qa4

    CAPAC(:)    = cat_progn(:,ens_id)%capac
    CATDEF(:)   = cat_progn(:,ens_id)%catdef
    RZEXC(:)    = cat_progn(:,ens_id)%rzexc
    SRFEXC(:)   = cat_progn(:,ens_id)%srfexc
    GHTCNT1(:)  = cat_progn(:,ens_id)%ght(1)
    GHTCNT2(:)  = cat_progn(:,ens_id)%ght(2)
    GHTCNT3(:)  = cat_progn(:,ens_id)%ght(3)
    GHTCNT4(:)  = cat_progn(:,ens_id)%ght(4)
    GHTCNT5(:)  = cat_progn(:,ens_id)%ght(5)
    GHTCNT6(:)  = cat_progn(:,ens_id)%ght(6)

    WESNN1(:)   = cat_progn(:,ens_id)%wesn(1)
    WESNN2(:)   = cat_progn(:,ens_id)%wesn(2)
    WESNN3(:)   = cat_progn(:,ens_id)%wesn(3)

    HTSNNN1(:)  = cat_progn(:,ens_id)%htsn(1)
    HTSNNN2(:)  = cat_progn(:,ens_id)%htsn(2)
    HTSNNN3(:)  = cat_progn(:,ens_id)%htsn(3)

    SNDZN1(:)   = cat_progn(:,ens_id)%sndz(1)
    SNDZN2(:)   = cat_progn(:,ens_id)%sndz(2)
    SNDZN3(:)   = cat_progn(:,ens_id)%sndz(3)

    if(ens_id == NUM_ENSEMBLE ) ens_id = 0

    ! End
    RETURN_(ESMF_SUCCESS)

end subroutine UPDATE_ASSIM

subroutine read_pert_rseed(seed_fname,pert_rseed_r8)
     use netcdf
     character(len=*),intent(in) :: seed_fname
     real(kind=ESMF_KIND_R8),intent(inout) :: pert_rseed_r8(:)

     integer :: ncid, s_varid, en_dim, n_ens, id_varid, i, pos
     logical :: file_exist
     
     inquire (file = trim(seed_fname), exist=file_exist)
     if ( .not. file_exist) then
        pert_rseed_r8 = 0
        return
     endif

     call check( nf90_open(seed_fname, NF90_NOWRITE, ncid) )
  ! Get the varid of the data variable, based on its name.
     call check( nf90_inq_varid(ncid, "pert_rseed", s_varid) )
     call check( nf90_get_var(ncid, s_varid, pert_rseed_r8) )

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
end subroutine read_pert_rseed

subroutine write_pert_rseed(chk_fname, pert_rseed_r8)
     use netcdf
     character(len=*),intent(in) :: chk_fname
     real(kind=ESMF_KIND_R8),intent(in) :: pert_rseed_r8(:)
     character(len=*), parameter  :: SHORT_NAME = "SHORT_NAME"
     character(len=*), parameter  :: LONG_NAME  = "LONG_NAME"
     character(len=*), parameter  :: UNITS      = "UNITS"
     character(len=*), parameter  :: s_SHORT = "obspert_rseed"
     character(len=*), parameter  :: s_long  = "Observation_Perturbations_rseed"
     character(len=*), parameter  :: units_ = "1"
     integer :: nseeds
     integer :: ncid, s_varid
     integer :: seed_dimid

     nseeds   = size(pert_rseed_r8)

 ! Create the file. 
     call check( nf90_create(trim(chk_fname), nf90_clobber + NF90_NETCDF4, ncid) )
! Define the dimensions.
     call check( nf90_def_dim(ncid, "NRANDSEED",   nseeds, seed_dimid) )
     call check( nf90_def_var(ncid, 'pert_rseed',   NF90_DOUBLE, [seed_dimid], s_varid) )

  ! Assign attribute
     call check( nf90_put_att(ncid, s_varid, UNITS, units_) )
     call check( nf90_put_att(ncid, s_varid, SHORT_NAME, s_short) )
     call check( nf90_put_att(ncid, s_varid, LONG_NAME, s_long) )

  ! End define mode.
     call check( nf90_enddef(ncid) )

  ! write varaible
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
end subroutine write_pert_rseed

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

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    character(len=300) :: seed_fname
    character(len=300) :: fname_tpl
    character(len=300) :: out_path
    character(len=ESMF_MAXSTR) :: exp_id
    character(len=4) :: id_string
    character(len=14):: datestamp
    integer :: ens, nymd, nhms
    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Finalize"

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, out_path, Label="OUT_PATH:", DEFAULT="./", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, exp_id, Label="EXP_ID:", DEFAULT="exp_id", RC=STATUS)
    VERIFY_(STATUS)

    if (master_proc) then
       call finalize_obslog()
       Pert_rseed_r8 = Pert_rseed
       call MAPL_GetResource ( MAPL, fname_tpl, Label="LANDASSIM_OBSPERTRSEED_CHECKPOINT_FILE:", DEFAULT="landassim_obspertrseed%s_checkpoint", RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_DateStampGet( clock, datestamp, rc=status)
       VERIFY_(STATUS)

       read(datestamp(1:8),*)   nymd
       read(datestamp(10:13),*) nhms
       nhms = nhms*100
       do ens = 0, NUM_ENSEMBLE-1
          write(id_string,'(I4.4)') ens + FIRST_ENS_ID
          seed_fname = ""
          call ESMF_CFIOStrTemplate(seed_fname,fname_tpl,'GRADS', xid=id_string,nymd=nymd,nhms=nhms,stat=status)
          VERIFY_(STATUS)
          call write_pert_rseed(trim(seed_fname), Pert_rseed_r8(:,ens+1))
       enddo
    endif

    ! Call Finalize for every child
    call MAPL_GenericFinalize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    ! End
    RETURN_(ESMF_SUCCESS)

end subroutine Finalize

end module GEOS_LandAssimGridCompMod
