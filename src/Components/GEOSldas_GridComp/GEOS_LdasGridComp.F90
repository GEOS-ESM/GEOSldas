#include "MAPL_Generic.h"

!BOP
! !MODULE: GEOS_LdasGridCompMod - Module to combine children GridComps
module GEOS_LdasGridCompMod

  ! !USES

  use ESMF
  use MAPL_Mod

  use GEOS_MetforceGridCompMod, only: MetforceSetServices => SetServices
  use GEOS_LandGridCompMod, only: LandSetServices => SetServices
  use GEOS_LandPertGridCompMod, only: LandPertSetServices => SetServices
  use GEOS_EnsGridCompMod, only: EnsSetServices => SetServices
  use GEOS_LandAssimGridCompMod, only: LandAssimSetServices => SetServices

  use LDAS_EASE_conv, only: ease_inverse
  use LDAS_TileCoordType, only: tile_coord_type , T_TILECOORD_STATE, TILECOORD_WRAP
  use LDAS_TileCoordType, only: grid_def_type, io_grid_def_type
  use LDAS_TileCoordRoutines, only: get_tile_grid, get_ij_ind_from_latlon
  use LDAS_ConvertMod, only: esmf2ldas
  use LDAS_PertRoutinesMod, only: get_pert_grid
  use LDAS_ensdrv_functions,ONLY:  get_io_filename 
  use LDAS_DateTimeMod,ONLY: date_time_type
  use LDAS_ensdrv_mpi, only: MPI_tile_coord_type, MPI_grid_def_type
  use LDAS_ensdrv_mpi, only: init_MPI_types,mpicomm,numprocs,myid 
  use LDAS_ensdrv_mpi, only: root_proc
  use LDAS_ensdrv_init_routines, only: io_domain_files
  use LDAS_ensdrv_Globals, only: logunit,logit,root_logit,echo_clsm_ensdrv_glob_param
  use lsm_routines,  only: lsmroutines_echo_constants  
  use StieglitzSnow, only: StieglitzSnow_echo_constants
  use SurfParams,    only: SurfParams_init

  implicit none

  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  ! !DESCRIPTION: This gridded component (GC) combines the GridComps:
  !     METFORCE, LAND, LANDPERT, ENSAVG, and LANDASSIM
  !  into a new composite LDAS GricComp.
  !  Include later: LAKE, LANDICE(?), SALTWATER(?)

  !EOP

  include 'mpif.h'

  ! All children
  integer,allocatable :: LAND(:)
  integer,allocatable :: LANDPERT(:)
  integer             :: METFORCE, ENSAVG, LANDASSIM

  ! other global variables
  integer :: NUM_ENSEMBLE
  logical :: land_assim
  logical :: mwRTM
  
contains

  !BOP

  ! !IROTUINE: SetServices -- Set ESMF services for this component

  ! !INTERFACE:

  subroutine SetServices(gc, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc ! gridded component
    integer, optional                  :: rc ! return code

    ! ensemble set up:

    integer :: i
    integer,allocatable :: ens_id(:)
    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(ESMF_GridComp), pointer :: gcs(:)=>null() ! Children gridcomps
    character(len=ESMF_MAXSTR), pointer :: gcnames(:)=>null() ! Children's names
    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name
    character(len=ESMF_MAXSTR) :: id_string,childname, fmt_str
    character(len=ESMF_MAXSTR) :: LAND_ASSIM_STR, mwRTM_file
    integer                    :: ens_id_width
    ! Local variables
    type(T_TILECOORD_STATE), pointer :: tcinternal
    type(TILECOORD_WRAP) :: tcwrap

    type(ESMF_Config) :: CF
    integer :: LSM_CHOICE
    integer :: FIRST_ENS_ID
    
    ! Begin...

    ! Get my name and setup traceback handle
    Iam = 'SetServices'
    call ESMF_GridCompGet(gc, name=comp_name,CONFIG=CF, rc=status)
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


    ! Save the tile_coord variable as an internal state of the GridComp
    ! memory for tcwrap%tile_coord is allocated in Initialize()
    allocate(tcinternal, stat=status)
    VERIFY_(status)
    tcwrap%ptr => tcinternal
    call ESMF_UserCompSetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
   
    !create ensemble children
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)
    call MAPL_GetResource ( MAPL, NUM_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1,       RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, ens_id_width, Label="ENS_ID_WIDTH:",      DEFAULT=0,       RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, FIRST_ENS_ID, Label="FIRST_ENS_ID:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, LAND_ASSIM_STR, Label="LAND_ASSIM:", DEFAULT="NO", RC=STATUS)
    VERIFY_(STATUS)
    LAND_ASSIM_STR =  ESMF_UtilStringUpperCase(LAND_ASSIM_STR, rc=STATUS)
    VERIFY_(STATUS)
    land_assim = (trim(LAND_ASSIM_STR) /= 'NO')

    call MAPL_GetResource ( MAPL, mwRTM_file, Label="LANDASSIM_INTERNAL_RESTART_FILE:", DEFAULT='', RC=STATUS)
    VERIFY_(STATUS)
    mwRTM = ( len_trim(mwRTM_file) /= 0 )

    call MAPL_GetResource ( MAPL, LSM_CHOICE, Label="LSM_CHOICE:", DEFAULT=1, RC=STATUS)
    if (LSM_CHOICE /=1 ) then
      _ASSERT( .not. (mwRTM .or. land_assim), "CatchCN is Not Ready for assimilation or mwRTM")
    endif

    METFORCE = MAPL_AddChild(gc, name='METFORCE', ss=MetforceSetServices, rc=status)
    VERIFY_(status)

    allocate(ens_id(NUM_ENSEMBLE),LAND(NUM_ENSEMBLE),LANDPERT(NUM_ENSEMBLE))
    write (fmt_str, "(A2,I1,A1,I1,A1)") "(I", ens_id_width,".",ens_id_width,")"   ! BUG? only works for ens_id_width<10)    (reichle, 11 Jun 2020)
    do i=1,NUM_ENSEMBLE
       ens_id(i) = i-1 + FIRST_ENS_ID ! id start form FIRST_ENS_ID
       if(NUM_ENSEMBLE == 1 ) then
          id_string=''
       else
          write(id_string, fmt_str) ens_id(i)
       endif

       id_string=trim(id_string)

       childname='LANDPERT'//trim(id_string)
       LANDPERT(i) = MAPL_AddChild(gc, name=childname, ss=LandPertSetServices, rc=status)
       VERIFY_(status)

       childname='LAND'//trim(id_string)
       LAND(i) = MAPL_AddChild(gc, name=childname, ss=LandSetServices, rc=status)
       VERIFY_(status)
    enddo

    ENSAVG    = MAPL_AddChild(gc, name='ENSAVG', ss=EnsSetServices, rc=status)
    VERIFY_(status)
    
    if(land_assim .or. mwRTM ) then
       LANDASSIM = MAPL_AddChild(gc, name='LANDASSIM', ss=LandAssimSetServices, rc=status)
       VERIFY_(status)
    endif

    ! Connections
    do i=1,NUM_ENSEMBLE
    ! -METFORCE-feeds-LANDPERT's-imports-
       call MAPL_AddConnectivity(                                                  &
            gc,                                                                    &
            SHORT_NAME = ['Tair   ', 'Qair   ', 'Psurf  ', 'Rainf_C', 'Rainf  ',   &
                          'Snowf  ', 'LWdown ', 'SWdown ', 'SWnet  ', 'PARdrct',   &
                          'PARdffs', 'Wind   ', 'RefH   '],                        &
            SRC_ID = METFORCE,                                                     &
            DST_ID = LANDPERT(i),                                                  &
            rc = status                                                            &
            )
       VERIFY_(status)
    ! -LANDPERT-feeds-LAND's-imports-
       call MAPL_AddConnectivity(                                                  &
            gc,                                                                    &
            SRC_NAME = ['TApert         ', 'QApert         ', 'UUpert         ',   &
                        'UWINDLMTILEpert', 'VWINDLMTILEpert', 'PCUpert        ',   &
                        'PLSpert        ', 'SNOpert        ', 'DRPARpert      ',   &
                        'DFPARpert      ', 'DRNIRpert      ', 'DFNIRpert      ',   &
                        'DRUVRpert      ', 'DFUVRpert      ', 'LWDNSRFpert    '],  &
            SRC_ID = LANDPERT(i),                                                  &
            DST_NAME = ['TA         ', 'QA         ', 'UU         ', 'UWINDLMTILE',&
                        'VWINDLMTILE', 'PCU        ', 'PLS        ', 'SNO        ',&
                        'DRPAR      ', 'DFPAR      ', 'DRNIR      ', 'DFNIR      ',&
                        'DRUVR      ', 'DFUVR      ', 'LWDNSRF    '],              &
            DST_ID = LAND(i),                                                      &
            rc = status                                                            &
            )
          VERIFY_(status)
    ! -METFORCE-feeds-LAND's-imports-
       call MAPL_AddConnectivity(                                                  &
            gc,                                                                    &
            SRC_NAME = ['Psurf', 'RefH ',                                          &
                        'DUDP ', 'DUSV ', 'DUWT ', 'DUSD ', 'BCDP ', 'BCSV ',      &
                        'BCWT ', 'BCSD ', 'OCDP ', 'OCSV ', 'OCWT ', 'OCSD ',      &
                        'SUDP ', 'SUSV ', 'SUWT ', 'SUSD ', 'SSDP ', 'SSSV ' ],    &
            SRC_ID = METFORCE,                                                     &
            DST_NAME = ['PS  ', 'DZ  ',                                            &
                        'DUDP', 'DUSV', 'DUWT', 'DUSD', 'BCDP', 'BCSV',            &
                        'BCWT', 'BCSD', 'OCDP', 'OCSV', 'OCWT', 'OCSD',            &
                        'SUDP', 'SUSV', 'SUWT', 'SUSD', 'SSDP', 'SSSV' ],          &
            DST_ID = LAND(i),                                                      &
            rc = status                                                            &
            )
       VERIFY_(status)
    ! -LAND-feeds-LANDPERT's-imports-
       call MAPL_AddConnectivity(                                                                  &
            gc,                                                                                    &
            SRC_NAME =  ['TC     ','CATDEF ','RZEXC  ','SRFEXC ','WESNN1 ','WESNN2 ','WESNN3 ',    &
               'GHTCNT1','GHTCNT2','GHTCNT3','GHTCNT4','GHTCNT5','GHTCNT6',                        &
               'HTSNNN1','HTSNNN2','HTSNNN3','SNDZN1 ','SNDZN2 ','SNDZN3 '],                       &
            SRC_ID = LAND(i),                                                                      &
            DST_NAME =     ['TCPert     ','CATDEFPert ','RZEXCPert  ','SRFEXCPert ','WESNN1Pert ', &
              'WESNN2Pert ','WESNN3Pert ','GHTCNT1Pert','GHTCNT2Pert',                             &
              'GHTCNT3Pert','GHTCNT4Pert','GHTCNT5Pert','GHTCNT6Pert',                             &
              'HTSNNN1Pert','HTSNNN2Pert','HTSNNN3Pert','SNDZN1Pert ',                             &
              'SNDZN2Pert ','SNDZN3Pert '],                                                        &
            DST_ID = LANDPERT(i),                                                                  &
            rc = status                                                                            &
            )
       VERIFY_(status)
    enddo

    if(land_assim .or. mwRTM) then
       ! -LAND-feeds-LANDASSIM's-imports-
       ! Catchment model parameters from first LAND ens member, assumes no parameter perturbations!
       call MAPL_AddConnectivity(                                                                  &
            gc,                                                                                    &
            SHORT_NAME = ['POROS ', 'COND  ','PSIS  ','BEE   ','WPWET ','GNU   ','VGWMAX',         &
                          'BF1   ', 'BF2   ','BF3   ','CDCR1 ','CDCR2 ','ARS1  ',                  &
                          'ARS2  ', 'ARS3  ','ARA1  ','ARA2  ','ARA3  ','ARA4  ',                  &
                          'ARW1  ', 'ARW2  ','ARW3  ','ARW4  ','TSA1  ','TSA2  ','TSB1  ',         &
                          'TSB2  ', 'ATAU  ','BTAU  ','ITY   ','Z2CH  ' ],                         &
            SRC_ID = LAND(1),                                                                      &  ! Note (1) !
            DST_ID = LANDASSIM,                                                                    &
            rc = status                                                                            &
            )
       VERIFY_(status)
    endif

    call MAPL_TimerAdd(gc, name="Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="-LocStreamCreate", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name='Run', rc=status)
    VERIFY_(status)

    ! Terminate all imports
    call MAPL_TerminateImport(gc, ALL=.true., rc=status)
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
    !use MAPL_LatLonToCubeRegridderMod
    !use MAPL_CubeToLatLonRegridderMod
    !use MAPL_CubeToCubeRegridderMod
    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! The Initialize routine creates the grid, locstream for all surface
    ! (surf\_locstream). It then splits this 'Surface' locstream based on mask
    ! (land, lake etc.) and attaches the sub-locstream to the corresponding
    ! child's GridComp.

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name
    character(len=ESMF_MAXSTR) :: gridname

    ! ESMF variables
    type(ESMF_Grid) :: agrid ! atospheric grid
    type(ESMF_GridComp), pointer :: gcs(:)=>null() ! Children gridcomps
    character(len=ESMF_MAXSTR), pointer :: gcnames(:)=>null() ! Children's names
    type(ESMF_DELayout) :: layout
    type(ESMF_DistGrid) :: distgrid

    ! MAPL variables
    type(MAPL_LocStream) :: surf_locstream
    type(MAPL_LocStream) :: land_locstream
    type(MAPL_MetaComp), pointer :: MAPL=>null() ! GC's MAPL obj
    type(MAPL_MetaComp), pointer :: CHILD_MAPL=>null() ! Child's MAPL obj

    ! LDAS' tile_coord variable
    type(T_TILECOORD_STATE), pointer :: tcinternal
    type(TILECOORD_WRAP) :: tcwrap

    ! Misc variables
    character(len=ESMF_MAXSTR) :: tilingfile
    character(len=300)         :: out_path,decomf
    character(len=ESMF_MAXSTR) :: exp_id
    character(len=ESMF_MAXSTR) :: LDAS_logit
    character(len=ESMF_MAXSTR) :: LAND_PARAMS 
    character(len=ESMF_MAXSTR) :: grid_type 

    integer :: total_nt,land_nt_local,i,j
    real, pointer :: LandTileLats(:)
    real, pointer :: LandTileLons(:)
    integer, pointer :: local_id(:)
    real(ESMF_KIND_R8), pointer     :: centerX(:,:)
    real(ESMF_KIND_R8), pointer     :: centerY(:,:)

    logical :: isEASEv1
    logical :: isEASEv2
    integer :: I1,IN,J1,JN
    real :: lat,lon
    type(ESMF_VM) :: vm
    integer :: mpierr
    logical :: IamRoot

    type(tile_coord_type), dimension(:), pointer :: tile_coord_f => null()
    type(tile_coord_type), dimension(:), pointer :: tile_coord_l => null()

    integer,dimension(:),pointer :: f2g
    integer :: N_catf

    type(grid_def_type) :: tile_grid_g
    type(grid_def_type) :: tile_grid_f
    type(grid_def_type) :: tile_grid_l
    type(date_time_type):: start_time
    type(ESMF_Time)     :: CurrentTime
    !type(CubedSphereGridFactory) :: cubed_sphere_factory
    !type (CubeToLatLonRegridder) :: cube_to_latlon_prototype
    !type (LatLonToCubeRegridder) :: latlon_to_cube_prototype
    !type (CubeToCubeRegridder) :: cube_to_cube_prototype
    real :: DT, DT_Solar
    type(ESMF_Alarm) :: SolarAlarm
    type(ESMF_TimeInterval) :: Solar_DT

    ! Begin...

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Initialize"
    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! MPI
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    IAmRoot = MAPL_Am_I_Root(vm)
    call ESMF_VmGet(vm, mpicommunicator=mpicomm, rc=status)
    VERIFY_(status)
    call MPI_COMM_RANK(mpicomm, myid,mpierr)
    call MPI_COMM_SIZE(mpicomm, numprocs, mpierr )
    root_proc = IAmRoot
    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Initialize")

    call ESMF_ClockGet(clock, currTime=CurrentTime, rc=status)
    VERIFY_(status)
    call esmf2ldas(CurrentTime, start_time, rc=status)
    VERIFY_(status)

    call MAPL_GetResource(MAPL,LDAS_logit,'LDAS_logit:',default = "NO",rc = status)
    VERIFY_(status)

    logit = (trim(LDAS_logit) /= 'NO')
    root_logit = (IamRoot .and. logit)

    ! Init catchment constants, currently different in GCM and GEOSldas
    call MAPL_GetResource(MAPL, LAND_PARAMS,Label="LAND_PARAMS:",DEFAULT="Icarus",RC=STATUS)
    VERIFY_(STATUS)
    call SurfParams_init(LAND_PARAMS)

    
    call MAPL_GetResource(MAPL, grid_type,Label="GEOSldas.GRID_TYPE:",RC=STATUS)
    VERIFY_(STATUS)

   ! if (trim(grid_type) == "Cubed-Sphere") then
   !    call grid_manager%add_prototype("Cubed-Sphere", cubed_sphere_factory)
   !    associate (method => REGRID_METHOD_BILINEAR, mgr => regridder_manager)
   !       call mgr%add_prototype('Cubed-Sphere', 'LatLon', method, cube_to_latlon_prototype)
   !       call mgr%add_prototype('LatLon', 'Cubed-Sphere', method, latlon_to_cube_prototype)
   !       call mgr%add_prototype('Cubed-Sphere', 'Cubed-Sphere', method, cube_to_cube_prototype)
   !    end associate
   ! endif

   ! Create atmospheric (single level atm grid covers all of surface) grid
    call MAPL_GridCreate(gc, rc=status)
    VERIFY_(status)

    ! Get grid info from the gridcomp
    call ESMF_GridCompGet(gc, grid=agrid, rc=status)
    VERIFY_(status)

    ! Get distgrid info from grid
    call ESMF_GridGet(agrid, distgrid=distgrid, rc=status)
    VERIFY_(status)

    ! Get DElayout info from distgrid
    call ESMF_DistGridGet(distgrid, delayout=layout, rc=status)
    VERIFY_(status)

    ! get grid name
    call ESMF_GridGet(agrid, name=gridname, rc=status)
    VERIFY_(STATUS)
    isEASEv1 =.false.
    isEASEv2 =.false.

    if (index(gridname,'EASEv2') /=0) then
       isEASEv2 = .true.
    else if (index(gridname,'EASE') /=0) then
       isEASEv1 = .true.
    endif
    if( isEASEv1) then
    ! To be implemented
    endif
    if( isEASEv2) then
! Retrieve the coordinates so we can set them
       call ESMF_GridGetCoord(agrid, coordDim=1, localDE=0, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=centerX, rc=status)
       VERIFY_(STATUS)

       call ESMF_GridGetCoord(agrid, coordDim=2, localDE=0, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=centerY, rc=status)
       VERIFY_(STATUS)

       call ESMF_GRID_INTERIOR(agrid,I1,IN,J1,JN)

       do I = 1,size(centerX,1)
          call ease_inverse(gridname,1.0*(I+I1-2),0.0,lat,lon)
          centerX(I,:) = lon
       enddo

       do J = 1,size(centerY,2)
          call ease_inverse(gridname,0.0,1.0*(J+J1-2),lat,lon)
          centerY(:,J) = lat
       enddo

    endif

    ! Get tile file location
    call MAPL_GetResource(                                                      &
         MAPL,                                                                  &
         tilingfile,                                                            &
         'TILING_FILE:',                                                        &
         default = "tile.data",                                                 &
         rc = status                                                            &
         )
    VERIFY_(status)


    !print*," Create LocStream for all surface (land, lake, landice, saltwater)"
    call MAPL_TimerOn(MAPL, "-LocStreamCreate")
    call MAPL_LocStreamCreate(                                                  &
         surf_locstream,                                                        &
         layout = layout,                                                       &
         filename = tilingfile,                                                 &
         name = "Surface",                                                      &
         grid = agrid,                                                          &
         rc = status                                                            &
         )
    VERIFY_(status)

    call MAPL_TimerOff(MAPL, "-LocStreamCreate")

    ! Get children and their im/ex states from MAPL obj
    call MAPL_Get(MAPL, GCS=gcs, GCNAMES=gcnames, rc=status)
    VERIFY_(status)

    ! Create LAND's locstreams as subset of Surface locstream
    ! and add it to the children's MAPL objects

    call MAPL_TimerOn(MAPL, "-LocStreamCreate")
    call MAPL_LocStreamCreate(                                                  &
         land_locstream,                                                        &
         surf_locstream,                                                        &
         name=gcnames(LAND(1)),                                                    &
         mask=[MAPL_LAND],                                                      &
         rc=status                                                              &
         )
    VERIFY_(status)
    call MAPL_TimerOff(MAPL, "-LocStreamCreate")
    ! Convert LAND's LocStream to LDAS' tile_coord and save it in the GridComp
    ! -get-tile-information-from-land's-locstream-
    call MAPL_LocStreamGet(                                                     &
         land_locstream,                                                        &
         NT_LOCAL=land_nt_local,                                                &
         TILELATS=LandTileLats,                                                 &
         TILELONS=LandTileLons,                                                 &
         LOCAL_ID=local_id    ,                                                 &
         rc=status                                                              &
         )
    VERIFY_(status)

    ! -get-component's-internal-state-
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tcinternal => tcwrap%ptr
    ! -allocate-memory-for-tile-coord-
    allocate(tcinternal%tile_coord(land_nt_local), stat=status)
    VERIFY_(status)
    allocate(tcinternal%l2f(land_nt_local))
    VERIFY_(status)
    

    call MAPL_GetResource ( MAPL, out_path, Label="OUT_PATH:", DEFAULT="./", RC=STATUS)
    call MAPL_GetResource ( MAPL, exp_id, Label="EXP_ID:", DEFAULT="./", RC=STATUS)

    call init_MPI_types()
    
    call MPI_Reduce(land_nt_local,total_nt,1,MPI_INT,MPI_SUM,0,mpicomm,mpierr);

    decomf = get_io_filename(trim(out_path), trim(exp_id), 'ldas_domdecomp', date_time=start_time, &
             dir_name='rc_out', file_ext='.txt')

    if (IamRoot) then
       call io_domain_files('r',trim(out_path), trim(exp_id),N_catf,f2g,tile_coord_f,tile_grid_g,tile_grid_f,RC)
       ! WY notes: f2g == tile_coord_f%tile_id
       deallocate(f2g)
       print*, "Number of tiles: ", N_catf
       if(N_catf /= total_nt) then
         print*, "total_nt = ", total_nt
         stop "tiles number not equal"
       endif
       open(10,file= trim(decomf), action='write')
       write(10,*) N_catf
       close(10)
       call io_grid_def_type('w', logunit, tile_grid_f, 'tile_grid_f')

       block 
          type(grid_def_type) :: latlon_tmp_g
          integer :: perturbations

          call MAPL_GetResource(MAPL, perturbations, 'PERTURBATIONS:', default=0, rc=status)
          if(trim(grid_type) == "Cubed-Sphere" ) then

            ASSERT_(index(tile_grid_g%gridtype, 'c3') /=0)
            !1) save original index
            tile_coord_f%cs_i_indg = tile_coord_f%i_indg 
            tile_coord_f%cs_j_indg = tile_coord_f%j_indg 
            
            !2) generate a lat-lon grid for landpert and land assim ( 4*N_lonX3*N_lon)
            call get_pert_grid(tile_grid_g, latlon_tmp_g)
            tile_grid_g = latlon_tmp_g
            !3) change the index
            !   need to chang min_lon, max_lon, min_lat , max_lat? 
            do i = 1, N_catf
               call get_ij_ind_from_latlon(latlon_tmp_g,tile_coord_f(i)%com_lat,tile_coord_f(i)%com_lon, &
                 tile_coord_f(i)%i_indg,tile_coord_f(i)%j_indg)
            enddo
            !3) re-generate tile_grid_f in Lat-Lon
            call get_tile_grid(N_catf, tile_coord_f, tile_grid_g, tile_grid_f)
            
          endif
       end block 

    endif
    
    call MPI_BCAST(N_catf,1,MPI_INTEGER,0,mpicomm,mpierr)
    if (.not. IamRoot) allocate(tile_coord_f(N_catf))

    call MPI_BCAST(tile_coord_f,N_catf,    MPI_tile_coord_type,0,mpicomm, mpierr)
    call MPI_BCAST(tile_grid_g, 1,         MPI_grid_def_type,  0,mpicomm, mpierr)
    call MPI_BCAST(tile_grid_f, 1,         MPI_grid_def_type,  0,mpicomm, mpierr)

    block
      integer, allocatable :: f2tile_id(:), tile_id2f(:)
      integer :: max_id
      allocate(f2tile_id(N_catf))
      f2tile_id = tile_coord_f%tile_id

      max_id = maxval(f2tile_id)
      allocate(tile_id2f(max_id),source = 0)
      do i = 1, N_catf
         tile_id2f(f2tile_id(i)) = i
      enddo
      tcinternal%l2f = tile_id2f(local_id)
      tcinternal%tile_coord = tile_coord_f(tcinternal%l2f)
      deallocate(f2tile_id, tile_id2f)
    end block

    do i = 0, numprocs-1
      if( i == myid) then
         open(10,file= trim(decomf), action='write',position='append')
         do j = 1, land_nt_local
            write(10,*) local_id(j), myid
         enddo
         close(10)
      endif
      call MPI_Barrier(mpicomm,mpierr)
    enddo   

    allocate(tcinternal%tile_coord_f,source = tile_coord_f)
    
    call get_tile_grid(land_nt_local,tcinternal%tile_coord,tile_grid_g,tile_grid_l)
   
    ! re-arrange tile_coord_f

    tcinternal%grid_g = tile_grid_g
    tcinternal%grid_f = tile_grid_f
    tcinternal%grid_l = tile_grid_l

    call MAPL_GetObjectFromGC(gcs(METFORCE), CHILD_MAPL, rc=status)
    VERIFY_(status) ! CHILD = METFORCE
    call MAPL_Set(CHILD_MAPL, LocStream=land_locstream, rc=status)
    VERIFY_(status)

    call MAPL_GetObjectFromGC(gcs(ENSAVG), CHILD_MAPL, rc=status)
    VERIFY_(status) ! CHILD = ens_avg
    call MAPL_Set(CHILD_MAPL, LocStream=land_locstream, rc=status)
    VERIFY_(status)
    call ESMF_UserCompSetInternalState(gcs(METFORCE), 'TILE_COORD', tcwrap, status)
    VERIFY_(status)

    do i = 1,NUM_ENSEMBLE
       call MAPL_GetObjectFromGC(gcs(LAND(i)), CHILD_MAPL, rc=status)
       VERIFY_(status)
       call MAPL_Set(CHILD_MAPL, LocStream=land_locstream, rc=status)
       VERIFY_(status)
       call MAPL_GetObjectFromGC(gcs(LANDPERT(i)), CHILD_MAPL, rc=status)
       VERIFY_(status) ! CHILD = LANDPERT
       call MAPL_Set(CHILD_MAPL, LocStream=land_locstream, rc=status)
       VERIFY_(status)
       ! Add LAND's tile_coord to children's GridComps
       call ESMF_UserCompSetInternalState(gcs(LAND(i)), 'TILE_COORD', tcwrap, status)
       VERIFY_(status)
       call ESMF_UserCompSetInternalState(gcs(LANDPERT(i)), 'TILE_COORD', tcwrap, status)
       VERIFY_(status)
    enddo

    if (land_assim .or. mwRTM) then
       call MAPL_GetObjectFromGC(gcs(LANDASSIM), CHILD_MAPL, rc=status)
       VERIFY_(status) 
       call MAPL_Set(CHILD_MAPL, LocStream=land_locstream, rc=status)
       VERIFY_(status)

       call ESMF_UserCompSetInternalState(gcs(LANDASSIM), 'TILE_COORD', tcwrap, status)
       VERIFY_(status)
    endif

    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    VERIFY_(status)

   ! solar alarm is created in solar gridcomp. Since GEOSldas doesnot have that gridcomp, it is created here
   ! -create-nonsticky-alarm-
    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DT_Solar, Label="SOLAR_DT:", DEFAULT=DT, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_TimeIntervalSet(SOLAR_DT, s=NINT(DT_Solar), rc=status)
    VERIFY_(status)

    SolarAlarm = ESMF_AlarmCreate(                                         &
         clock,                                                                 &
         name='SOLAR_Alarm',                                                     &
         ringTime=CurrentTime,                                                  &
         ringInterval=SOLAR_DT,                                               &
         sticky=.false.,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)


    if ( IamRoot) call echo_clsm_ensdrv_glob_param()
    if ( IamRoot) call lsmroutines_echo_constants(logunit)
    if ( IamRoot) call StieglitzSnow_echo_constants(logunit)


    ! Turn timer off
    call MAPL_TimerOff(MAPL, "Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


  !BOP

  ! !IROTUINE: Run -- Run method for the composite Ldas GridComp

  ! !INTERFACE:

  subroutine Run(gc, import, export, clock, rc)

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! Calls children's Run methods.

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! ESMF variables
    type(ESMF_VM) :: vm
    type(ESMF_GridComp), pointer :: gcs(:)
    type(ESMF_State), pointer :: gim(:)
    type(ESMF_State), pointer :: gex(:)
    character(len=ESMF_MAXSTR), pointer :: gcnames(:)
    type(ESMF_Time) :: ModelTimeCur

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL

    ! Misc variables
    integer :: igc,i
    logical :: IAmRoot
    integer :: mpierr
    integer :: LSM_CHOICE

    ! Begin...

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run"

    ! Get MAPL obj
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Run")

    ! Get information about children
    call MAPL_Get(MAPL, GCS=gcs, GIM=gim, GEX=gex, GCNAMES=gcnames, rc=status)
    VERIFY_(STATUS)

    ! MPI
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    IAmRoot = MAPL_Am_I_Root(vm)
    !call ESMF_VmGet(vm, mpicommunicator=mpicomm, rc=status)
    !VERIFY_(status)

    call MAPL_GetResource ( MAPL, LSM_CHOICE, Label="LSM_CHOICE:", DEFAULT=1, RC=STATUS)

    ! Get current time
    call ESMF_ClockGet(clock, currTime=ModelTimeCur, rc=status)
    VERIFY_(status)
    if (IAmRoot) then
       call ESMF_TimePrint(ModelTimeCur, options='string', rc=status)
       VERIFY_(status)
    end if
    
    !phase2 initialization ( executed once)
    !adjust mean of perturbed forcing or Progn
    do i  = 1,NUM_ENSEMBLE
       igc = LANDPERT(i)
       call MAPL_TimerOn(MAPL, gcnames(igc))
       call ESMF_GridCompRun(gcs(igc), importState=gim(igc), exportState=gex(igc), clock=clock, phase=1, userRC=status)
       VERIFY_(status)
       call MAPL_TimerOff(MAPL, gcnames(igc))
    enddo

    ! Run children GridComps (in order)
    ! Generate raw perturbed force and progn
    do i  = 1,NUM_ENSEMBLE
       igc = LANDPERT(i)
       call MAPL_TimerOn(MAPL, gcnames(igc))
       call ESMF_GridCompRun(gcs(igc), importState=gim(igc), exportState=gex(igc), clock=clock, phase=2, userRC=status)
       VERIFY_(status)
       call MAPL_TimerOff(MAPL, gcnames(igc))
    enddo


    igc = METFORCE
    call MAPL_TimerOn(MAPL, gcnames(igc))
    call ESMF_GridCompRun(gcs(igc), importState=gim(igc), exportState=gex(igc), clock=clock, userRC=status)
    VERIFY_(status)
    call MAPL_TimerOff(MAPL, gcnames(igc))

    do i  = 1,NUM_ENSEMBLE

       !ApplyForcePert
       igc = LANDPERT(i)
       call MAPL_TimerOn(MAPL, gcnames(igc))
       call ESMF_GridCompRun(gcs(igc), importState=gim(igc), exportState=gex(igc), clock=clock, phase=3, userRC=status)
       VERIFY_(status)
       call MAPL_TimerOff(MAPL, gcnames(igc))

       ! Use landpert's output as the input to calculate the ensemble average forcing
       ! W.J note: So far it is only for the Catchment model. 
       ! To make CatchmentCN work with assim, the export from landgrid and catchmentCN grid need to be modified.  
       if ( LSM_CHOICE == 1 ) then
          call ESMF_GridCompRun(gcs(ENSAVG), importState=gex(igc), exportState=gex(ENSAVG), clock=clock,phase=1, userRC=status)
          VERIFY_(status)
       endif

       ! Run the land model
       igc = LAND(i)
       call MAPL_TimerOn(MAPL, gcnames(igc))
       call ESMF_GridCompRun(gcs(igc), importState=gim(igc), exportState=gex(igc), clock=clock, phase=1, userRC=status)
       VERIFY_(status)
       call ESMF_GridCompRun(gcs(igc), importState=gim(igc), exportState=gex(igc), clock=clock, phase=2, userRC=status)
       VERIFY_(status)
       call MAPL_TimerOff(MAPL, gcnames(igc))

       ! ApplyPrognPert - moved: now before calculating ensemble average that is picked up by land analysis and HISTORY; reichle 28 May 2020 
       igc = LANDPERT(i)
       call MAPL_TimerOn(MAPL, gcnames(igc))
       call ESMF_GridCompRun(gcs(igc), importState=gim(igc), exportState=gex(igc), clock=clock, phase=4, userRC=status)
       VERIFY_(status)
       call MAPL_TimerOff(MAPL, gcnames(igc))

       ! Use LAND's output as the input to calculate the ensemble average
       igc = LAND(i)
       if (LSM_CHOICE == 1) then
          ! collect cat_param 
          call ESMF_GridCompRun(gcs(ENSAVG), importState=gex(igc), exportState=gex(ENSAVG), clock=clock,phase=3, userRC=status)
          VERIFY_(status)
          call ESMF_GridCompRun(gcs(ENSAVG), importState=gex(igc), exportState=gex(ENSAVG), clock=clock,phase=2, userRC=status)
          VERIFY_(status)

          if( mwRTM ) then
             ! Calculate ensemble-average L-band Tb using LAND's output (add up and normalize after last member has been added)
             call ESMF_GridCompRun(gcs(LANDASSIM), importState=gex(igc), exportState=gex(LANDASSIM), clock=clock,phase=3, userRC=status)
             VERIFY_(status)
          endif
       endif

    enddo

    ! Run land analysis
    if (land_assim) then 
       igc = LANDASSIM
       call MAPL_TimerOn(MAPL, gcnames(igc))
       ! Get EnKF increments and apply to "cat_progn" (imported from ENSAVG via "use" statement!); otherwise import state is export from ENSAVG
       call ESMF_GridCompRun(gcs(igc), importState=gex(ENSAVG), exportState=gex(igc), clock=clock, phase=1, userRC=status)
       VERIFY_(status)

       do i = 1, NUM_ENSEMBLE
          ! Extract updated exports from "cat_progn" 
          call ESMF_GridCompRun(gcs(igc), importState=gim(igc), exportState=gex(LAND(i)), clock=clock, phase=2, userRC=status)
          VERIFY_(status)
          
       enddo
       call MAPL_TimerOff(MAPL, gcnames(igc))
    endif

    ! Turn timers off
    call MAPL_TimerOff(MAPL, "Run")
    call MAPL_TimerOff(MAPL, "TOTAL")

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Run


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
    ! Clean-up.

    !EOP

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! Internal state variables
    type(T_TILECOORD_STATE), pointer :: tcinternal
    type(TILECOORD_WRAP) :: tcwrap

    ! Begin...

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Finalize"

    ! Get component's internal state
    call ESMF_UserCompGetInternalState(gc, 'TILE_COORD', tcwrap, status)
    VERIFY_(status)
    tcinternal => tcwrap%ptr

    ! ! Clean up internal state's tile_coord variable
    ! if (associated(tcinternal%tile_coord)) deallocate(tcinternal%tile_coord)

    ! Call Finalize for every child
    call MAPL_GenericFinalize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Finalize

end module GEOS_LdasGridCompMod
