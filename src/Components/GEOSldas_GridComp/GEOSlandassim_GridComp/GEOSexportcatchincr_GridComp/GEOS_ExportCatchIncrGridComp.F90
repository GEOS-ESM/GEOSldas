#include "MAPL_Generic.h"

!=============================================================================
module GEOS_ExportCatchIncrGridCompMod

  !BOP
  ! !DESCRIPTION:
  !
  !   This is a gridded component to export analysis increments for Catchment.
  !   It has no children.
  
  !
  ! !USES:
  
  use ESMF
  use MAPL_Mod
  use ESMF_CFIOMOD,              only: ESMF_CFIOstrTemplate
  
  implicit none

  include 'mpif.h'
  
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices
  !
  !EOP
  !
  integer, parameter :: NUM_SUBTILES = 4

contains
  
  ! ******************************************************************************

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
    
    character(len=ESMF_MAXSTR)   :: Iam
    character(len=ESMF_MAXSTR)   :: COMP_NAME
    integer                      :: STATUS
    
    ! Local Variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(ESMF_Config)            :: CF

    ! Begin...
    ! --------
    
    ! Get my name and set-up traceback handle
    ! ------------------------------------------------------------------------------
    
    Iam='SetServices'
    call ESMF_GridCompGet ( GC, NAME=COMP_NAME, RC=STATUS )
    _VERIFY(STATUS)
    Iam=trim(COMP_NAME)//trim(Iam)
    
    ! Register services for this component
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_INITIALIZE,                                                &
         Initialize,                                                            &
         rc=status                                                              &
         )
    _VERIFY(status)
    
    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_FINALIZE,                                                  &
         Finalize,                                                              &
         rc=status                                                              &
         )
    _VERIFY(status)
        
  ! Exports for Catchment prognostics increments

  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_canopy_temperature_saturated_zone'         ,&
       UNITS              = 'K'                                                   ,&
       SHORT_NAME         = 'TCFSAT_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_canopy_temperature_transition_zone'        ,&
       UNITS              = 'K'                                                   ,&
       SHORT_NAME         = 'TCFTRN_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_canopy_temperature_wilting_zone'           ,&
       UNITS              = 'K'                                                   ,&
       SHORT_NAME         = 'TCFWLT_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_canopy_specific_humidity_saturated_zone'   ,&
       UNITS              = 'kg kg-1'                                             ,&
       SHORT_NAME         = 'QCFSAT_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_canopy_specific_humidity_transition_zone'  ,&
       UNITS              = 'kg kg-1'                                             ,&
       SHORT_NAME         = 'QCFTRN_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)

  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_canopy_specific_humidity_wilting_zone'     ,&
       UNITS              = 'kg kg-1'                                             ,&
       SHORT_NAME         = 'QCFWLT_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_interception_reservoir_capac'              ,&
       UNITS              = 'kg m-2'                                              ,&
       SHORT_NAME         = 'CAPAC_INCR'                                          ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_catchment_deficit'                         ,&
       UNITS              = 'kg m-2'                                              ,&
       SHORT_NAME         = 'CATDEF_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_root_zone_excess'                          ,&
       UNITS              = 'kg m-2'                                              ,&
       SHORT_NAME         = 'RZEXC_INCR'                                          ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
    
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_surface_excess'                            ,&
       UNITS              = 'kg m-2'                                              ,&
       SHORT_NAME         = 'SRFEXC_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_soil_heat_content_layer_1'                 ,&
       UNITS              = 'J m-2'                                               ,&
       SHORT_NAME         = 'GHTCNT1_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_soil_heat_content_layer_2'                 ,&
       UNITS              = 'J_m-2'                                               ,&
       SHORT_NAME         = 'GHTCNT2_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_soil_heat_content_layer_3'                 ,&
       UNITS              = 'J m-2'                                               ,&
       SHORT_NAME         = 'GHTCNT3_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_soil_heat_content_layer_4'                 ,&
       UNITS              = 'J m-2'                                               ,&
       SHORT_NAME         = 'GHTCNT4_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)

  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_soil_heat_content_layer_5'                 ,&
       UNITS              = 'J m-2'                                               ,&
       SHORT_NAME         = 'GHTCNT5_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_soil_heat_content_layer_6'                 ,&
       UNITS              = 'J m-2'                                               ,&
       SHORT_NAME         = 'GHTCNT6_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_snow_mass_layer_1'                         ,&
       UNITS              = 'kg m-2'                                              ,&
       SHORT_NAME         = 'WESNN1_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_snow_mass_layer_2'                         ,&
       UNITS              = 'kg m-2'                                              ,&
       SHORT_NAME         = 'WESNN2_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_snow_mass_layer_3'                         ,&
       UNITS              = 'kg m-2'                                              ,&
       SHORT_NAME         = 'WESNN3_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_heat_content_snow_layer_1'                 ,&
       UNITS              = 'J m-2'                                               ,&
       SHORT_NAME         = 'HTSNNN1_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_heat_content_snow_layer_2'                 ,&
       UNITS              = 'J m-2'                                               ,&
       SHORT_NAME         = 'HTSNNN2_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_heat_content_snow_layer_3'                 ,&
       UNITS              = 'J m-2'                                               ,&
       SHORT_NAME         = 'HTSNNN3_INCR'                                        ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_snow_depth_layer_1'                        ,&
       UNITS              = 'm'                                                   ,&
       SHORT_NAME         = 'SNDZN1_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_snow_depth_layer_2'                        ,&
       UNITS              = 'm'                                                   ,&
       SHORT_NAME         = 'SNDZN2_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)

  call MAPL_AddExportSpec(GC                                                      ,&
       LONG_NAME          = 'increment_snow_depth_layer_3'                        ,&
       UNITS              = 'm'                                                   ,&
       SHORT_NAME         = 'SNDZN3_INCR'                                         ,&
       DIMS               = MAPL_DimsTileOnly                                     ,&
       VLOCATION          = MAPL_VLocationNone                                    ,&
       RC=STATUS  )
  _VERIFY(STATUS)
  
    call MAPL_TimerAdd(GC, name="Initialize"    ,RC=STATUS)
    _VERIFY(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    _VERIFY(STATUS)

    call MAPL_GenericSetServices ( GC, RC=STATUS )
    _VERIFY(STATUS)
    
    RETURN_(ESMF_SUCCESS)
    
  end subroutine SetServices
  
  !BOP
  ! !INTERFACE:
  subroutine Initialize(gc, import, export, clock, rc)
    ! !ARGUMENTS:
    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code
    
    integer :: status
    character(len=ESMF_MAXSTR)   :: Iam
    character(len=ESMF_MAXSTR)   :: comp_name
    type(MAPL_MetaComp), pointer :: MAPL=>null()   
 
    ! Begin...
    
    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    _VERIFY(status)
    Iam = trim(comp_name) // "::Initialize"
    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    _VERIFY(status)
    
    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Initialize")
    
    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    _VERIFY(status)
       
    call MAPL_TimerOff(MAPL, "Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")
    _RETURN(ESMF_SUCCESS)
    
  end subroutine Initialize
  
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
    integer                      :: status
    character(len=ESMF_MAXSTR)   :: Iam
    character(len=ESMF_MAXSTR)   :: comp_name

    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    _VERIFY(status)
    Iam = trim(comp_name) // "::Finalize"
    
    call MAPL_GenericFinalize(gc, import, export, clock, rc=status)
    _VERIFY(status)
    
    RETURN_(ESMF_SUCCESS)
    
  end subroutine Finalize

end module GEOS_ExportCatchIncrGridCompMod 
