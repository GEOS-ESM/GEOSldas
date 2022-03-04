
! ********************************************************************************

!   IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT  

!    When a derived type is defined elsewhere in the code and then replicated 
!    as an MPI STRUCTURE here, make sure the original type is defined in such 
!    a way that it is compatible with the "-align" compiler flag.
!
!    See page 106 of "MPI: A Message-Passing Interface Standard, Version 3.0", 
!    Message Passing Interface Forum, September 21, 2012:
!
!    Structures combining different basic datatypes should be defined so that
!     there will be no gaps based on alignment rules.  If such a datatype is used
!     to create an array of structures, users should *also* avoid an alignment-gap
!     at the end of the structure. 
!
!     [...]
!
!     Example: Instead of
!    
!       TYPE, BIND(C) :: my_data
!       REAL, DIMENSION(3) :: x
!       ! there may be a gap of the size of one REAL
!       ! if the alignment of a DOUBLE PRECISION is
!       ! two times the size of a REAL
!       DOUBLE PRECISION :: p
!       END TYPE
!
!     one should define
!
!       TYPE, BIND(C) :: my_data
!       REAL, DIMENSION(3) :: x
!       REAL :: gap1
!       DOUBLE PRECISION :: p
!       END TYPE
!       
!     and also include gap1 in the matching MPI derived datatype.

!   IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT  

! ********************************************************************************

module LDAS_ensdrv_mpi

  use catch_constants, only: N_gt => CATCH_N_GT
  
  use catch_types,     only: N_cat_progn, N_cat_diagS, N_cat_diagF

  use enkf_types,      only: N_obs_ang_max
  
  implicit none
  
  include 'mpif.h'
  
  public :: init_MPI_types
  ! initialize to non-MPI values
 
  integer, public  :: myid=0, numprocs=1, mpicomm
  integer, public  :: mpierr, mpistatus(MPI_STATUS_SIZE)
  
  logical, public  :: root_proc=.true.
  
  integer, public  :: MPI_tile_coord_type, MPI_grid_def_type
  integer, public  :: MPI_cat_param_type,  MPI_cat_progn_type
  integer, public  :: MPI_cat_diagS_type,  MPI_cat_diagF_type
  integer, public  :: MPI_met_force_type,  MPI_veg_param_type
  integer, public  :: MPI_bal_diagn_type,  MPI_alb_param_type  
  integer, public  :: MPI_date_time_type
  integer, public  :: MPI_mwRTM_param_type,MPI_obs_type,        MPI_obs_param_type
  integer, public  :: MPI_cat_progn_int_type,                   MPI_cat_bias_param_type
  integer, public  :: MPI_obs_bias_type

  integer, private :: N_real, N_int
  
  
contains
  
  ! *****************************************************************************
  
  subroutine init_MPI_types()
    
    integer                                                   :: icount
    integer,                        allocatable, dimension(:) :: iblock, itype
    integer(KIND=MPI_ADDRESS_KIND), allocatable, dimension(:) :: idisp

    ! ---------------------------------------------------------------------
    !
    ! WARNING: do not confuse date_time_type with the f90 intrisic
    !          function "date_and_time()"
    !
    !  type :: date_time_type              
    !   integer :: year               ! 4-digit year
    !   integer :: month              ! month in year
    !   integer :: day                ! day in month
    !   integer :: hour               ! hour of day
    !   integer :: min                ! minute of hour
    !   integer :: sec                ! seconds of minute
    !   integer :: pentad             ! pentad of year
    !   integer :: dofyr              ! day of year
    !  end type date_time_type

    icount = 1

    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_INTEGER
    
    iblock(1) = 8
    
    idisp(1)  = 0
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_date_time_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_date_time_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)


    ! --------------------------------------------------------------------------------
    !
    !  type :: tile_coord_type
    !     
    !     integer :: tile_id      ! unique tile ID
    !     integer :: f_num        ! full domain ID
    !     integer :: typ          ! (0=MAPL_Ocean, 100=MAPL_Land, 19=MAPL_Lake, 20=MAPL_LandIce) 
    !     integer :: pfaf         ! Pfafstetter number (for land tiles, NOT unique)
    !     real    :: com_lon      ! center-of-mass longitude
    !     real    :: com_lat      ! center-of-mass latitude
    !     real    :: min_lon      ! minimum longitude (bounding box for tile)
    !     real    :: max_lon      ! maximum longitude (bounding box for tile)
    !     real    :: min_lat      ! minimum latitude (bounding box for tile)
    !     real    :: max_lat      ! maximum latitude (bounding box for tile)
    !     integer :: i_indg       ! i index (w.r.t. *global* grid that cuts tiles) 
    !     integer :: j_indg       ! j index (w.r.t. *global* grid that cuts tiles)
    !     ! For cubed-sphere tile spaces, hash_[x]_indg refers to a lat-lon "hash" grid that will 
    !     !   be created at runtime to support efficient mapping for perturbations and the EnKF analysis.
    !     ! For EASE and LatLon tile spaces, hash_[x]_indg is identical to [x]_indg
    !     integer :: hash_i_indg  ! i index (w.r.t. *global* "hash" grid for perts and EnKF) 
    !     integer :: hash_j_indg  ! j index (w.r.t. *global* "hash" grid for perts and EnKF)
    !     real    :: frac_cell    ! area fraction of grid cell covered by tile
    !     real    :: frac_pfaf    ! fraction of Pfafstetter catchment for land tiles 
    !     real    :: area         ! area [km^2]
    !     real    :: elev         ! elevation above sea level [m]

    icount = 4

    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_INTEGER
    itype(2)  = MPI_REAL
    itype(3)  = MPI_INTEGER
    itype(4)  = MPI_REAL
    
    iblock(1) = 4
    iblock(2) = 6
    iblock(3) = 4 ! i_indg, j_indg, hash_i_indg and hash_j_indg
    iblock(4) = 4
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    idisp(3)  = idisp(2) + iblock(2)*4
    idisp(4)  = idisp(3) + iblock(3)*4
    
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_tile_coord_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_tile_coord_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    ! --------------------------------------------------------------------------------
    !
    !  type :: grid_def_type              
    ! 
    ! character(40) :: gridtype ! grid type, eg. "LatLon", "EASE_M36", "EASE_M09", ...
    ! integer       :: ind_base ! 0=zero-based indices (EASE), 1=one-based indices (LatLon)
    ! integer       :: i_dir    ! direction of indices (+1=W-to-E, -1=E-to-W)
    ! integer       :: j_dir    ! direction of indices (+1=S-to-N, -1=N-to-S)
    ! integer       :: N_lon    ! number of longitude nodes
    ! integer       :: N_lat    ! number of latitude nodes
    ! integer       :: i_offg   ! minimum lon index (offset from global grid)
    ! integer       :: j_offg   ! minimum lat index (offset from global grid)
    ! real          :: ll_lon   ! lower left  longitude of grid cell edge [deg]
    ! real          :: ll_lat   ! lower left  latitude of grid cell edge [deg]
    ! real          :: ur_lon   ! upper right longitude of grid cell edge [deg]
    ! real          :: ur_lat   ! upper right latitude of grid cell edge [deg]
    ! real          :: dlon     ! longitude grid spacing [deg]
    ! real          :: dlat     ! latitude grid spacing [deg]
    
    icount = 3

    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_CHARACTER
    itype(2)  = MPI_INTEGER
    itype(3)  = MPI_REAL

    iblock(1) = 40
    iblock(2) = 7
    iblock(3) = 6
        
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*1    ! each character is 1 byte
    idisp(3)  = idisp(2) + iblock(2)*4    ! each integer   is 4 bytes
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_grid_def_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_grid_def_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    ! ---------------------------------------------------------------------------------
    !
    ! type cat_param_type
    ! 
    ! real :: dpth  ! depth to bedrock from data file (dpth/=dzpr in general!)
    ! real :: dzsf  ! "surface layer"                     formerly zdep1
    ! real :: dzrz  ! "root zone layer"                   formerly zdep2
    ! real :: dzpr  ! "profile layer" (unsaturated zone)  formerly zdep3
    ! real, dimension(N_gt) :: dzgt  
    ! real :: poros   ! porosity
    ! real :: cond    ! saturated hydraulic conductivity
    ! real :: psis    ! Clapp-Hornberger parameter
    ! real :: bee     ! Clapp-Hornberger parameter
    ! real :: wpwet   ! wilting poing wetness
    ! real :: gnu     ! vertical decay factor for transmissivity
    ! real :: vgwmax      ! max amount of water available to vegetation 
    ! integer :: vegcls         ! vegetation class
    ! integer :: soilcls30      ! soil_class_top (0- 30cm)
    ! integer :: soilcls100     ! soil_class_com (0-100cm)
    ! real :: bf1
    ! real :: bf2
    ! real :: bf3
    ! real :: cdcr1
    ! real :: cdcr2
    ! real :: ars1 
    ! real :: ars2
    ! real :: ars3
    ! real :: ara1
    ! real :: ara2
    ! real :: ara3
    ! real :: ara4
    ! real :: arw1
    ! real :: arw2
    ! real :: arw3
    ! real :: arw4
    ! real :: tsa1
    ! real :: tsa2
    ! real :: tsb1
    ! real :: tsb2
    ! real :: atau
    ! real :: btau
    ! real :: gravel30
    ! real :: orgC30  
    ! real :: orgC    
    ! real :: sand30  
    ! real :: clay30  
    ! real :: sand    
    ! real :: clay    
    ! real :: wpwet30 
    ! real :: poros30 
    ! real :: veghght

    icount = 3
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_REAL
    itype(2)  = MPI_INTEGER
    itype(3)  = MPI_REAL

    iblock(1) = 4+N_gt+7
    iblock(2) = 3
    iblock(3) = 32
        
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    idisp(3)  = idisp(2) + iblock(2)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_cat_param_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_cat_param_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)
    
    ! ---------------------------------------------------------------------------------
    !
    ! integer, parameter :: N_cat_progn = 10 + N_gt + 3*N_snow
    !
    ! type :: cat_progn_type
    !     
    ! real :: tc1     ! surface/canopy temperature 
    ! real :: tc2
    ! real :: tc4
    ! real :: qa1     ! specific humidity in canopy air
    ! real :: qa2
    ! real :: qa4
    ! real :: capac   ! canopy interception water
    ! real :: catdef  ! catchment deficit
    ! real :: rzexc   ! root zone excess
    ! real :: srfexc  ! surface excess
    ! real, dimension(N_gt)   :: ght     ! ground heat content
    ! real, dimension(N_snow) :: wesn    ! snow water equivalent
    ! real, dimension(N_snow) :: htsn    ! snow heat content
    ! real, dimension(N_snow) :: sndz    ! snow depth
    
    N_real = N_cat_progn
    
    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    ! split MPI type into two blocks of real numbers
    ! (having just one with N_cat_progn MPI_REAL entries did not work)
    
    itype(1)  = MPI_REAL
    itype(2)  = MPI_REAL
    
    iblock(1) = 1
    iblock(2) = N_real-1
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_cat_progn_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_cat_progn_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)
    
    ! ---------------------------------------------------------------------------------
    !
    ! integer, parameter :: N_cat_diagS = 7 + N_gt + N_snow
    !
    ! type :: cat_diagS_type
    ! 
    ! real :: ar1      ! area fraction of saturated zone
    ! real :: ar2      ! area fraction of unsaturated and unstressed zone
    ! real :: asnow    ! area fraction of snow
    ! real :: sfmc     ! surface moisture content
    ! real :: rzmc     ! root zone moisture content
    ! real :: prmc     ! profile moisture content
    ! real :: tsurf    ! mean surface temperature over entire catchment 
    ! real, dimension(N_gt)   :: tp     ! temperature of soil layers
    ! real, dimension(N_snow) :: tpsn   ! temperature of snow layers
    
    N_real = N_cat_diagS

    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    ! split MPI type into two blocks of real numbers
    ! (having just one with N_cat_diagS MPI_REAL entries did not work)

    itype(1)  = MPI_REAL
    itype(2)  = MPI_REAL
    
    iblock(1) = 1
    iblock(2) = N_real-1
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_cat_diagS_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_cat_diagS_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)
        
    ! --------------------------------------------------------------------------------
    !
    ! integer, parameter :: N_cat_diagF = 22
    !
    ! type :: cat_diagF_type
    ! 
    ! real :: shflux   ! sensible heat flux
    ! real :: lhflux   ! total latent heat flux
    ! real :: ghflux   ! ground heat flux to top soil layer
    ! real :: evap     ! total evaporation 
    ! real :: eint     ! interception loss
    ! real :: esoi     ! evaporation from bare soil
    ! real :: eveg     ! transpiration 
    ! real :: esno     ! evaporation from snow
    ! real :: runoff   ! total runoff
    ! real :: runsrf   ! surface runoff
    ! real :: bflow    ! baseflow
    ! real :: snmelt   ! snow melt
    ! real :: lwup     ! outgoing/upward longwave radiation
    ! real :: swup     ! outgoing/upward shortwave radiation
    ! real :: qinfil   ! infiltration
    ! real :: hsnacc   ! accounting term for energy related to snowfall etc.
    ! real :: evacc    ! accounting term for evaporation   (see catchment()) 
    ! real :: shacc    ! accounting term for sensible heat (see catchment()) 
    ! real :: lhacc    ! accounting term for latent heat   (see catchment()) 
    ! real :: eacc_0   ! accounting term for oscillations  (see catchment()) 
    ! real :: t2m      ! air temperature at 2m above the displacement height
    ! real :: q2m      ! specific humidity at 2m above the displacement height
    
    N_real = N_cat_diagF

    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    ! split MPI type into two blocks of real numbers
    ! (having just one with N_cat_diagF MPI_REAL entries did not work)

    itype(1)  = MPI_REAL
    itype(2)  = MPI_REAL
    
    iblock(1) = 1
    iblock(2) = N_real-1
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_cat_diagF_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_cat_diagF_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)
        
    ! --------------------------------------------------------------------------------
    !
    ! type met_force_type
    !
    ! real :: Tair                 ! air temperature at RefH                 [K]
    ! real :: Qair                 ! specific humidity at RefH               [kg/kg]
    ! real :: Psurf                ! surface pressure                        [Pa]
    ! real :: Rainf_C              ! convective rainfall                     [kg/m2/s]
    ! real :: Rainf                ! total rainfall                          [kg/m2/s]
    ! real :: Snowf                ! total snowfall                          [kg/m2/s]
    ! real :: LWdown               ! downward longwave radiation             [W/m2]
    ! real :: SWdown               ! downward shortwave radiation            [W/m2]
    ! real :: PARdrct              ! Photosynth. Active Radiation (direct)   [W/m2]
    ! real :: PARdffs              ! Photosynth. Active Radiation (diffuse)  [W/m2]
    ! real :: Wind                 ! wind speed at RefH                      [m/s]
    ! real :: RefH                 ! reference height for Tair, Qair, Wind   [m]
    
    N_real = 12

    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    ! split MPI type into two blocks of real numbers
    ! (having just one with N_met_force MPI_REAL entries did not work)
    
    itype(1)  = MPI_REAL
    itype(2)  = MPI_REAL
    
    iblock(1) = 1
    iblock(2) = N_real-1
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_met_force_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_met_force_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    ! --------------------------------------------------------------------------------
    !    
    ! type veg_param_type
    ! 
    ! real :: grn                  ! vegetation greenness fraction           [-]
    ! real :: lai                  ! leaf-area-index                         [m2/m2]
    
    N_real = 2   

    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    ! split MPI type into two blocks of real numbers
    ! (having just one with N_veg_param MPI_REAL entries did not work)
    
    itype(1)  = MPI_REAL
    itype(2)  = MPI_REAL
    
    iblock(1) = 1
    iblock(2) = N_real-1
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_veg_param_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_veg_param_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    ! --------------------------------------------------------------------------------
    !    
    ! type bal_diagn_type
    ! 
    ! real :: etotl                ! total energy store (in all model progn)  [J/m2]
    ! real :: echng                ! energy change per unit time (model only) [W/m2]
    ! real :: eincr                ! energy analysis increment per unit time  [W/m2]
    ! real :: wtotl                ! total water store (in all model progn)   [kg/m2]
    ! real :: wchng                ! water change per unit time (model only)  [kg/m2/s]
    ! real :: wincr                ! water analysis increment per unit time   [kg/m2/s]
    
    N_real = 6

    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    ! split MPI type into two blocks of real numbers
    ! (having just one with N_bal_diagn MPI_REAL entries did not work)
    
    itype(1)  = MPI_REAL
    itype(2)  = MPI_REAL
    
    iblock(1) = 1
    iblock(2) = N_real-1
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_bal_diagn_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_bal_diagn_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    ! --------------------------------------------------------------------------------
    !    
    ! type alb_param_type
    ! 
    ! reichle,  5 Apr 2013 - removed alb_param_type fields "sc_albvr" and "sc_albnr" 
    !
    ! !real :: sc_albvr    ! Scaling factor for direct  visible  or blacksky 0.3-0.7 
    ! !real :: sc_albnr    ! Scaling factor for direct  infrared or blacksky 0.7-5.0
    ! real :: sc_albvf    ! Scaling factor for diffuse visible  or whitesky 0.3-0.7 
    ! real :: sc_albnf    ! Scaling factor for diffuse infrared or whitesky 0.7-5.0
    
    N_real = 2   

    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    ! split MPI type into two blocks of real numbers
    ! (having just one with N_alb_param MPI_REAL entries did not work)
    
    itype(1)  = MPI_REAL
    itype(2)  = MPI_REAL
    
    iblock(1) = 1
    iblock(2) = N_real-1
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_alb_param_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_alb_param_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)
    
    ! --------------------------------------------------------------------------------
    !
    ! type mwRTM_param_type
    ! 
    ! integer :: vegcls     
    ! integer :: soilcls    
    ! real    :: sand       
    ! real    :: clay       
    ! real    :: poros
    ! real    :: wang_wt    
    ! real    :: wang_wp    
    ! real    :: rgh_hmin   
    ! real    :: rgh_hmax   
    ! real    :: rgh_wmin   
    ! real    :: rgh_wmax   
    ! real    :: rgh_Nrh    
    ! real    :: rgh_Nrv    
    ! real    :: rgh_polmix 
    ! real    :: omega      
    ! real    :: bh         
    ! real    :: bv         
    ! real    :: lewt       
    ! real    :: vegopacity  

    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_INTEGER
    itype(2)  = MPI_REAL

    iblock(1) = 2
    iblock(2) = 17
        
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_mwRTM_param_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_mwRTM_param_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    ! ---------------------------------------------------------------------------------
    !
    ! type obs_type
    ! 
    !  logical :: assim     
    !  integer :: species   
    !  integer :: tilenum   
    !  integer :: DUMMYGAP  ! fill gap so that MPI STRUCT works with "-align" compiler flag
    !  real*8  :: time
    !  real    :: lon       
    !  real    :: lat       
    !  real    :: obs       
    !  real    :: obsvar    
    !  real    :: fcst      
    !  real    :: fcstvar   
    !  real    :: ana       
    !  real    :: anavar    

    icount = 4
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_LOGICAL
    itype(2)  = MPI_INTEGER
    itype(3)  = MPI_REAL8
    itype(4)  = MPI_REAL

    iblock(1) = 1
    iblock(2) = 3
    iblock(3) = 1
    iblock(4) = 8
        
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4          
    idisp(3)  = idisp(2) + iblock(2)*4
    idisp(4)  = idisp(3) + iblock(3)*8    ! real*8
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_obs_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_obs_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    ! ---------------------------------------------------------------------------------
    !
    ! type obs_param_type
    ! 
    !  character(40)                 :: descr                 ! block  #1 (character)
    !  integer                       :: species               ! block  #2 (integer)
    !  integer                       :: orbit   
    !  integer                       :: pol     
    !  integer                       :: N_ang   
    !  real, &
    !       dimension(N_obs_ang_max) :: ang                   ! block  #3 (real)
    !  real                          :: freq    
    !  real                          :: FOV     
    !  character(40)                 :: FOV_units             ! block  #4 (character)
    !  logical                       :: assim                 ! block  #5 (logical)
    !  logical                       :: scale   
    !  logical                       :: getinnov
    !  integer                       :: RTM_ID                ! block  #6 (integer)  
    !  integer                       :: bias_Npar
    !  integer                       :: bias_trel             
    !  integer                       :: bias_tcut
    !  real                          :: nodata                ! block  #7 (real)
    !  character(40)                 :: varname               ! block  #8 (character)
    !  character(40)                 :: units
    !  character(200)                :: path    
    !  character(80)                 :: name    
    !  character(200)                :: scalepath
    !  character(80)                 :: scalename
    !  character(200)                :: flistpath
    !  character(80)                 :: flistname
    !  real                          :: errstd                ! block  #9 (real)
    !  real                          :: std_normal_max
    !  logical                       :: zeromean              ! block #10 (logical)
    !  logical                       :: coarsen_pert       
    !  real                          :: xcorr                 ! block #11 (real)
    !  real                          :: ycorr         
    !  integer                       :: adapt                 ! block #12 (integer)
     
    icount = 12
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))

    itype( 1)  = MPI_CHARACTER    
    itype( 2)  = MPI_INTEGER
    itype( 3)  = MPI_REAL
    itype( 4)  = MPI_CHARACTER
    itype( 5)  = MPI_LOGICAL
    itype( 6)  = MPI_INTEGER
    itype( 7)  = MPI_REAL
    itype( 8)  = MPI_CHARACTER
    itype( 9)  = MPI_REAL
    itype(10)  = MPI_LOGICAL
    itype(11)  = MPI_REAL
    itype(12)  = MPI_INTEGER

    iblock( 1) = 40
    iblock( 2) = 4
    iblock( 3) = N_obs_ang_max+2
    iblock( 4) = 40
    iblock( 5) = 3  
    iblock( 6) = 4
    iblock( 7) = 1                  
    iblock( 8) = 40+40+200+80+200+80+200+80      
    iblock( 9) = 2                  
    iblock(10) = 2                  
    iblock(11) = 2                  
    iblock(12) = 1                  
    
    idisp( 1)  = 0
    idisp( 2)  = idisp( 1) + iblock( 1)         ! CHARACTER*1 !!!
    idisp( 3)  = idisp( 2) + iblock( 2)*4
    idisp( 4)  = idisp( 3) + iblock( 3)*4
    idisp( 5)  = idisp( 4) + iblock( 4)         ! CHARACTER*1 !!!
    idisp( 6)  = idisp( 5) + iblock( 5)*4
    idisp( 7)  = idisp( 6) + iblock( 6)*4
    idisp( 8)  = idisp( 7) + iblock( 7)*4
    idisp( 9)  = idisp( 8) + iblock( 8)         ! CHARACTER*1 !!!
    idisp(10)  = idisp( 9) + iblock( 9)*4
    idisp(11)  = idisp(10) + iblock(10)*4
    idisp(12)  = idisp(11) + iblock(11)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_obs_param_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_obs_param_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)




    !--------------------------------------------------------------------------------
    !
    ! type :: cat_progn_int_type
    !
    !   integer :: tc1     ! surface/canopy temperature 
    !   integer :: tc2
    !   integer :: tc4
    !   integer :: qa1     ! specific humidity in canopy air
    !   integer :: qa2
    !   integer :: qa4
    !   integer :: capac   ! canopy interception water
    !   integer :: catdef  ! catchment deficit
    !   integer :: rzexc   ! root zone excess
    !   integer :: srfexc  ! surface excess
    !   integer, dimension(N_gt)   :: ght     ! ground heat content
    !   integer, dimension(N_snow) :: wesn    ! snow water equivalent
    !   integer, dimension(N_snow) :: htsn    ! snow heat content
    !   integer, dimension(N_snow) :: sndz    ! snow depth
    
    N_int  = N_cat_progn
    
    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_INTEGER
    itype(2)  = MPI_INTEGER
    
    iblock(1) = 1
    iblock(2) = N_int-1
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_cat_progn_int_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_cat_progn_int_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    ! ---------------------------------------------------------------------------------
    !
    ! type :: cat_bias_param_type
    !
    !   type(cat_progn_type)     :: tconst
    !   type(cat_progn_type)     :: trelax
    !   type(cat_progn_int_type) :: Nparam
    
    N_real = 2*N_cat_progn
    N_int  =   N_cat_progn
    
    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_REAL
    itype(2)  = MPI_INTEGER
    
    iblock(1) = N_real
    iblock(2) = N_int
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_cat_bias_param_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_cat_bias_param_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)

    
    ! ---------------------------------------------------------------------------------
    !
    ! type :: obs_bias_type
    !
    !   real       	      :: bias
    !   integer, dimension(2) :: tcount 
    ! 

    icount = 2
    
    allocate(iblock(icount))
    allocate(idisp( icount))
    allocate(itype( icount))
    
    itype(1)  = MPI_REAL
    itype(2)  = MPI_INTEGER
    
    iblock(1) = 1
    iblock(2) = 2
    
    idisp(1)  = 0
    idisp(2)  = idisp(1) + iblock(1)*4
    
    call MPI_TYPE_CREATE_STRUCT( icount, iblock, idisp, itype, &
         MPI_obs_bias_type, mpierr )
    
    call MPI_TYPE_COMMIT(MPI_obs_bias_type, mpierr)
    
    deallocate(iblock)
    deallocate(idisp)
    deallocate(itype)


    ! ---------------------------------------------------------------------------------
    
  end subroutine init_MPI_types
  
  ! *****************************************************************************
  
!  subroutine mpi_call_out(my_message)
!    
!    character(200) :: my_message
!    
!    character(  8) :: date_string
!    character( 10) :: time_string
!
!    character(  2) :: tmpstr2
!
!    call date_and_time(date_string, time_string)
!    
!    write (tmpstr2,'(i2.2)') myid
!
!    write (*,*) trim(my_message), ": myid ", tmpstr2, ', ', date_string, time_string
!    
!  end subroutine mpi_call_out
  
  ! *****************************************************************************
  
end module LDAS_ensdrv_mpi
