!
! type definitions for module to generate land surface perturbations
!
! reichle,  14 Apr 2006 - split land_pert.F90 into 2 files to avoid 
!                         having more than one module per file
!
! ------------------------------------------------------------

module LDAS_PertTypes
  
  ! reichle, 26 May 2005
  
  use ESMF
  use LDAS_TileCoordType, only: grid_def_type

  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: pert_param_type
  public :: allocate_pert_param
  public :: deallocate_pert_param

  public :: T_LANDPERT_STATE
  public :: LANDPERT_WRAP  

  ! --------------------------------------------------------------------
  !
  ! parameters for each kind of perturbation (precip, radiation, 
  ! soil moisture, etc)
  
  type :: pert_param_type

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          calls to MPI_BCAST in clsm_ensdrv_main.F90 are also updated.
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     character(40)     :: descr    ! 'pcp' (precip), 'sw' (shortwave), etc

     ! add or multiply perturbation?
     !
     ! additive:                     typ = 0
     ! multiplicative and lognormal: typ = 1

     integer           :: typ      
     
     ! max allowed normalized perturbation (relative to N(0,1))
     
     real              :: std_normal_max  
     
     ! if .true. enforce zeromean across ensemble 
     ! (implies mean=1 for multiplicative perturbations)
     ! (not applicable if only one ensemble member is done at a time)
     
     logical           :: zeromean        ! enforce zero mean across ensemble

     ! Allow perturbations to be computed on coarsened grid?
     ! Coarse grid spacing automatically determined as a function of model 
     ! grid spacing and spatial correlation scales (see random_fields.F90)

     logical           :: coarsen
     
     ! Mean and std are allowed to vary in space (dimension(N_x,N_y)).
     
     real, dimension(:,:), pointer :: mean      ! mean
     real, dimension(:,:), pointer :: std       ! standard deviation
     
     ! Cross-correlations between different kinds of perturbations
     ! (eg. between precip and shortwave perturbations) are allowed to vary
     ! in space (dimension(N_pert_kind,N_x,N_y)).
     
     real, dimension(:,:,:), pointer :: ccorr
     
     ! Spatial and temporal correlation scales must be constant in space.
     ! For non-zero cross-correlations they must also be the same for
     ! all kinds for perturbations (eg. if precip and radiation
     ! perturbations are cross-correlated, their xcorr, ycorr and tcorr
     ! must be the same).
     
     real             :: xcorr  ! correlation length along latitudes   [deg]
     real             :: ycorr  ! correlation length along longitudes  [deg]
     real             :: tcorr  ! temporal correlation length          [s]
     
  end type pert_param_type
  
  ! **********************************************************************

  type T_PERT
    ! private
     integer :: npert ! number of perturbations
     integer :: dtstep
     type(ESMF_Time) :: TimePrv, TimeNxt
     type(pert_param_type), pointer :: param(:)=>null()
     real, allocatable :: DataPrv(:,:), DataNxt(:,:)
   end type T_PERT

  ! Internal state and its wrapper
   type T_LANDPERT_STATE
     !private
      integer :: PERTURBATIONS ! 1: perturb variables; 0: no perturbation
      integer :: ens_id
      integer :: NUM_ENSEMBLE
      logical :: isCubedSphere
     ! pert grids - local and full
      type(grid_def_type) :: pgrid_l, pgrid_f,pgrid_g
      integer,allocatable :: i_indgs(:)
      integer,allocatable :: j_indgs(:)
     ! if it is cubed-sphere grid, swith to internal start
      real,allocatable    :: fpert_ntrmdt(:,:,:)
      real,allocatable    :: ppert_ntrmdt(:,:,:)
      real(kind=ESMF_KIND_R8), allocatable :: pert_rseed_r8(:)
     ! force/progn perturbations
      type(T_PERT) :: ForcePert, PrognPert
   end type T_LANDPERT_STATE
   type LANDPERT_WRAP
     type(T_LANDPERT_STATE), pointer :: ptr=>null()
   end type LANDPERT_WRAP
  
contains  
  
  subroutine allocate_pert_param(N_pert, N_x, N_y, pert_param)
    
    implicit none
    
    integer, intent(in) :: N_pert, N_x, N_y
    
    type(pert_param_type), dimension(:), pointer :: pert_param
    
    ! local variables
    
    integer :: k
    
    ! --------------------------------------------------------
    
    nullify(pert_param)
        
    allocate(pert_param(N_pert))
    
    do k=1,N_pert
       
       allocate(pert_param(k)%mean(N_x,N_y))
       allocate(pert_param(k)%std(N_x,N_y))
       allocate(pert_param(k)%ccorr(N_pert,N_x,N_y))
       
    end do
    
  end subroutine allocate_pert_param
  
  ! **********************************************************************

  subroutine deallocate_pert_param(N_pert, pert_param)
    
    implicit none
    
    integer, intent(in) :: N_pert
    
    type(pert_param_type), dimension(:), pointer :: pert_param
    
    ! local variables
    
    integer :: k
    
    ! --------------------------------------------------------
    
    do k=1,N_pert
       
       deallocate(pert_param(k)%mean)
       deallocate(pert_param(k)%std)
       deallocate(pert_param(k)%ccorr)
       
    end do
    
    deallocate(pert_param)
        
  end subroutine deallocate_pert_param

  ! **********************************************************************
  
end module LDAS_PertTypes


! =============== EOF =================================================
