
! type definitions for tile coordinates and domain 
!
! reichle, 26 Jan 2005
! reichle, 14 Apr 2006 - split tile_coord.F90 into 2 files to avoid 
!                         having more than one module per file
! reichle,  2 Aug 2020 - removed tile_typ_* (use MAPL_Ocean, MAPL_Land, etc instead) 
!
! ========================================================================
!
! type definition for tile coordinates

module LDAS_TileCoordType
  
  ! TODO: Replace ldas_abort with MAPL_ABORT
  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use LDAS_ensdrv_Globals,              ONLY:     &
       nodata_generic

  use iso_fortran_env

  implicit none
  
  private

  public :: tile_coord_type
  public :: grid_def_type

  public :: init_grid_def_type
  public :: io_grid_def_type
  public :: io_tile_coord_type
  public :: T_TILECOORD_STATE
  public :: TILECOORD_WRAP
  
  public :: operator (==) 
  
  ! ------------------------------------------------------------

  type :: tile_coord_type
     
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     integer :: tile_id      ! unique tile ID
     integer :: f_num        ! full domain ID
     integer :: typ          ! (0=MAPL_Ocean, 100=MAPL_Land, 19=MAPL_Lake, 20=MAPL_LandIce) 
     integer :: pfaf         ! Pfafstetter number (for land tiles, NOT unique)
     real    :: com_lon      ! center-of-mass longitude
     real    :: com_lat      ! center-of-mass latitude
     real    :: min_lon      ! minimum longitude (bounding box for tile)
     real    :: max_lon      ! maximum longitude (bounding box for tile)
     real    :: min_lat      ! minimum latitude (bounding box for tile)
     real    :: max_lat      ! maximum latitude (bounding box for tile)
     integer :: i_indg       ! i index (w.r.t. *global* grid that cuts tiles) 
     integer :: j_indg       ! j index (w.r.t. *global* grid that cuts tiles)
     ! For cubed-sphere tile spaces, pert_[x]_indg refers to a lat-lon perturbation grid that will 
     !   be created at runtime to support efficient mapping for perturbations and the EnKF analysis.
     ! For EASE and LatLon tile spaces, pert_[x]_indg is identical to [x]_indg
     integer :: pert_i_indg  ! i index (w.r.t. *global* pert_grid for perts and EnKF) 
     integer :: pert_j_indg  ! j index (w.r.t. *global* pert_grid for perts and EnKF)
     real    :: frac_cell    ! area fraction of grid cell covered by tile
     real    :: frac_pfaf    ! fraction of Pfafstetter catchment for land tiles 
     real    :: area         ! area [km^2]
     real    :: elev         ! elevation above sea level [m]

     
  end type tile_coord_type

  
  ! ------------------------------------------------------------
  ! 
  ! definition of *rectangular* (regular) grid
  !
  ! Possible grid types (structure field "gridtype"):
  !
  !   - "LatLon"     : regular lat/lon grid    (constant dlon, dlat)
  !   - "EASEv1_Mxx" : cylindrical EASEv1 grid (constant dlon, variable dlat)
  !   - "EASEv2_Mxx" : cylindrical EASEv2 grid (constant dlon, variable dlat)
  !
  ! Grid orientation (convention):
  !
  !  "LatLon"     : 1-based indexing, SouthWest to NorthEast
  !  "EASEv1_Mxx" : 0-based indexing, NorthWest to SouthEast
  !  "EASEv2_Mxx" : 0-based indexing, NorthWest to SouthEast
  !
  ! Grids are defined by the following fields:
  !    
  !  ---------------------------------------------------------
  !  |           ||     "LatLon"       |    "EASEv1_Mxx"     |
  !  |           ||                    |    "EASEv2_Mxx"     |
  !  ---------------------------------------------------------
  !  | indexing  ||  ind_base          |  ind_base           |
  !  |           ||  i_dir,  j_dir     |  i_dir,  j_dir      |             
  !  ---------------------------------------------------------
  !  | extent    ||  N_lon,  N_lat     |  N_lon,  N_lat      |             
  !  ---------------------------------------------------------
  !  | position  ||  ll_lon, ll_lat    |  i_offg, j_offg     |               
  !  ---------------------------------------------------------
  !  | spacing   ||  dlon,   dlat      |  gridtype ('Mxx')   |
  !  ---------------------------------------------------------
  !
  ! All other fields are derived from the above parameters 
  !  (except "descr", which is just a descriptor)
  !
  ! Lon, lat convention is -180<=lon<=180, -90<=lat<=90
  !
  ! ll_lon/ll_lat denote the coordinates of the lower left hand
  !  corner of the lower left grid cell (=southwestern corner of domain)
  ! ur_lon/ur_lat denote the coordinates of the upper right hand
  !  corner of the upper right grid cell (=northeastern corner of domain)
  
  type :: grid_def_type              

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     character(40) :: gridtype ! type of grid, eg. "LatLon", "EASEv1_M36", "EASEv2_M09", ...
     integer       :: ind_base ! 0=zero-based indices (EASE), 1=one-based indices (LatLon)
     integer       :: i_dir    ! direction of indices (+1=W-to-E, -1=E-to-W)
     integer       :: j_dir    ! direction of indices (+1=S-to-N, -1=N-to-S)
     integer       :: N_lon    ! number of longitude nodes
     integer       :: N_lat    ! number of latitude nodes
     integer       :: i_offg   ! minimum lon index (offset from *global* grid)
     integer       :: j_offg   ! minimum lat index (offset from *global* grid)
     !                         !  "LatLon"    : i_offg -> westernmost longitude
     !                         !                j_offg -> southernmost latitude
     !                         !  "EASEv1_Mxx": i_offg -> westernmost longitude
     !                         !                j_offg -> northernmost latitude
     !                         !  "EASEv2_Mxx": i_offg -> westernmost longitude
     !                         !                j_offg -> northernmost latitude
     real          :: ll_lon   ! lower left  longitude of grid cell edge [deg]
     real          :: ll_lat   ! lower left  latitude of grid cell edge [deg]
     real          :: ur_lon   ! upper right longitude of grid cell edge [deg]
     real          :: ur_lat   ! upper right latitude of grid cell edge [deg]
     real          :: dlon     ! longitude grid spacing [deg]
     real          :: dlat     ! latitude grid spacing [deg]
     !                         !  "LatLon"    : constant
     !                         !  "EASEv1_Mxx": *average* dlat over grid extent
     !                         !  "EASEv2_Mxx": *average* dlat over grid extent

     !CLSM_ENSDRV_MPI is NOT updated. will not be saved and bcasted
     real(kind=REAL64), dimension(:,:,:), pointer :: LonEdge =>null()
     real(kind=REAL64), dimension(:,:,:), pointer :: LatEdge =>null()
          
  end type grid_def_type

  ! ------------------------------------------------------------
  
  ! Wrapper for tile_coord structure
  
  type T_TILECOORD_STATE
     type(tile_coord_type), pointer, contiguous :: tile_coord(:)=>null()
     type(tile_coord_type), pointer, contiguous :: tile_coord_f(:)=>null()
     integer, pointer    :: l2f(:)=>null()
     type(grid_def_type) :: tgrid_g         ! tile_grid_g
     type(grid_def_type) :: pgrid_g         ! pert_grid_g
     type(grid_def_type) :: pgrid_f         ! pert_grid_f
     type(grid_def_type) :: pgrid_l         ! pert_grid_l
  end type T_TILECOORD_STATE

  type TILECOORD_WRAP
     type(T_TILECOORD_STATE), pointer :: ptr=>null()
  end type TILECOORD_WRAP

  ! --------------------------------------------------------------
  
  interface operator (==)
     module procedure eq_grid_def_type
  end interface
  
  ! *******************************************************************

contains

  subroutine init_grid_def_type( grid )
    
    implicit none
    
    type(grid_def_type), intent(out) :: grid
    
    grid%gridtype = 'undefined'
    grid%ind_base = nint(nodata_generic)
    grid%i_dir    = nint(nodata_generic)
    grid%j_dir    = nint(nodata_generic)
    grid%N_lon    = nint(nodata_generic)
    grid%N_lat    = nint(nodata_generic)
    grid%i_offg   = nint(nodata_generic)
    grid%j_offg   = nint(nodata_generic)
    grid%ll_lon   = nodata_generic
    grid%ll_lat   = nodata_generic
    grid%ur_lon   = nodata_generic
    grid%ur_lat   = nodata_generic
    grid%dlon     = nodata_generic
    grid%dlat     = nodata_generic
    
  end subroutine init_grid_def_type

  ! *******************************************************************

  function eq_grid_def_type( grid_1, grid_2 )
    
    ! reichle, 24 July 2010

    logical :: eq_grid_def_type

    type(grid_def_type), intent(in) :: grid_1, grid_2
    
    if ( index(grid_1%gridtype,trim(grid_2%gridtype))>0        .and. &
         grid_1%ind_base == grid_2%ind_base                    .and. &
         grid_1%i_dir    == grid_2%i_dir                       .and. &
         grid_1%j_dir    == grid_2%j_dir                       .and. &
         grid_1%N_lon    == grid_2%N_lon                       .and. &
         grid_1%N_lat    == grid_2%N_lat                       .and. &
         grid_1%i_offg   == grid_2%i_offg                      .and. &
         grid_1%j_offg   == grid_2%j_offg                      .and. &
         abs(grid_1%ll_lon-grid_2%ll_lon)<1e-4                 .and. &
         abs(grid_1%ll_lat-grid_2%ll_lat)<1e-4                 .and. &
         abs(grid_1%ur_lon-grid_2%ur_lon)<1e-4                 .and. &
         abs(grid_1%ur_lat-grid_2%ur_lat)<1e-4                 .and. &
         abs(grid_1%dlon  -grid_2%dlon  )<1e-4                 .and. &
         abs(grid_1%dlat  -grid_2%dlat  )<1e-4                       &
         )                                                     then
       
       eq_grid_def_type = .true.
       
    else
       
       eq_grid_def_type = .false.
       
    end if
    
  end function eq_grid_def_type
  

  ! *******************************************************************

  subroutine io_grid_def_type( action, unitnum, grid, varname  )
    
    ! reichle, 24 July 2010
    
    implicit none
    
    character, intent(in) :: action
    
    integer, intent(in) :: unitnum
    
    type(grid_def_type), intent(inout) :: grid
    
    character(*), intent(in), optional :: varname
    
    ! local variables
    
    character(40) :: vname, tmpstr40, tmpstr40b

    character(len=*), parameter :: Iam = 'io_grid_def_type'
    character(len=400) :: err_msg
    
    ! -------------------------------------------------------------
    
    if (present(varname)) then
       vname = varname
    else
       vname = 'grid'
    end if
       
    inquire(unit=unitnum, form=tmpstr40)
    
    select case (action)
       
    case ('w','W')

       if ( index(tmpstr40,'UNFORMATTED') /=0 .or.          &
            index(tmpstr40,'unformatted') /=0     ) then
          
          ! unformatted output
          
          write (unitnum)                    & 
               grid%gridtype, grid%ind_base,    &     
               grid%i_dir,    grid%j_dir,       &        
               grid%N_lon,    grid%N_lat,       &        
               grid%i_offg,   grid%j_offg,      &       
               grid%ll_lon,   grid%ll_lat,      &       
               grid%ur_lon,   grid%ur_lat,      &       
               grid%dlon,     grid%dlat         
          
       else
          
          ! write formatted file for easy reading into Matlab
          
          write(unitnum,*) trim(vname)//".gridtype = '", trim(grid%gridtype),"' ;"
          write(unitnum,*) trim(vname)//".ind_base =  ",      grid%ind_base, "  ;"
          write(unitnum,*) trim(vname)//".i_dir    =  ",      grid%i_dir,    "  ;"  
          write(unitnum,*) trim(vname)//".j_dir    =  ",      grid%j_dir,    "  ;" 
          write(unitnum,*) trim(vname)//".N_lon    =  ",      grid%N_lon,    "  ;"  
          write(unitnum,*) trim(vname)//".N_lat    =  ",      grid%N_lat,    "  ;" 
          write(unitnum,*) trim(vname)//".i_offg   =  ",      grid%i_offg,   "  ;"
          write(unitnum,*) trim(vname)//".j_offg   =  ",      grid%j_offg,   "  ;"
          write(unitnum,*) trim(vname)//".ll_lon   =  ",      grid%ll_lon,   "  ;"
          write(unitnum,*) trim(vname)//".ll_lat   =  ",      grid%ll_lat,   "  ;"
          write(unitnum,*) trim(vname)//".ur_lon   =  ",      grid%ur_lon,   "  ;"
          write(unitnum,*) trim(vname)//".ur_lat   =  ",      grid%ur_lat,   "  ;"
          write(unitnum,*) trim(vname)//".dlon     =  ",      grid%dlon,     "  ;"  
          write(unitnum,*) trim(vname)//".dlat     =  ",      grid%dlat,     "  ;"  
          write(unitnum,*)
       
       end if
       
    case ('r','R')
       
       if ( index(tmpstr40,'UNFORMATTED') /=0 .or.          &
            index(tmpstr40,'unformatted') /=0     ) then
          
          ! unformatted output
          
          read (unitnum)                     & 
               grid%gridtype, grid%ind_base,    &     
               grid%i_dir,    grid%j_dir,       &        
               grid%N_lon,    grid%N_lat,       &        
               grid%i_offg,   grid%j_offg,      &       
               grid%ll_lon,   grid%ll_lat,      &       
               grid%ur_lon,   grid%ur_lat,      &       
               grid%dlon,     grid%dlat         
          
       else
                    
          ! read formatted file
          
          read(unitnum,*) tmpstr40, tmpstr40b, grid%gridtype  
          read(unitnum,*) tmpstr40, tmpstr40b, grid%ind_base  
          read(unitnum,*) tmpstr40, tmpstr40b, grid%i_dir       
          read(unitnum,*) tmpstr40, tmpstr40b, grid%j_dir      
          read(unitnum,*) tmpstr40, tmpstr40b, grid%N_lon       
          read(unitnum,*) tmpstr40, tmpstr40b, grid%N_lat      
          read(unitnum,*) tmpstr40, tmpstr40b, grid%i_offg    
          read(unitnum,*) tmpstr40, tmpstr40b, grid%j_offg    
          read(unitnum,*) tmpstr40, tmpstr40b, grid%ll_lon    
          read(unitnum,*) tmpstr40, tmpstr40b, grid%ll_lat    
          read(unitnum,*) tmpstr40, tmpstr40b, grid%ur_lon    
          read(unitnum,*) tmpstr40, tmpstr40b, grid%ur_lat    
          read(unitnum,*) tmpstr40, tmpstr40b, grid%dlon        
          read(unitnum,*) tmpstr40, tmpstr40b, grid%dlat        
          read(unitnum,*)

       end if
       
    case default
       
       err_msg = 'unknown action ' // action
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
  end subroutine io_grid_def_type
  
  ! *******************************************************************

  subroutine io_tile_coord_type( action, unitnum, N_tile, tile_coord  )
    
    ! reichle, 24 July 2010
    ! reichle,  2 May  2013 - changed N_tile to intent(in)
    ! reichle,  7 Jan  2014 - changed to binary (unformatted) I/O
    ! wjiang, reichle, 18 Aug 2020 - Added initialization of %pert_[x]_indg during read.
    !                                Note that %pert_[x]_indg is NOT written out.

    implicit none
    
    character, intent(in) :: action
    
    integer,   intent(in) :: unitnum

    integer,                                  intent(in)    :: N_tile
    
    type(tile_coord_type), dimension(N_tile), intent(inout) :: tile_coord
    
    ! local

    integer       :: n, istat, N_tile_tmp

    character(len=*), parameter :: Iam = 'io_tile_coord_type'
    character(len=400) :: err_msg
    integer, allocatable :: tmp_int(:)
    real , allocatable :: tmp_real(:)
 
    ! -------------------------------------------------------------
    
    select case (action)
       
    case ('w','W')
       
       write (unitnum) N_tile
       
       write (unitnum) (tile_coord(n)%tile_id,   n=1,N_tile)
       write (unitnum) (tile_coord(n)%typ,       n=1,N_tile)
       write (unitnum) (tile_coord(n)%pfaf,      n=1,N_tile)
       write (unitnum) (tile_coord(n)%com_lon,   n=1,N_tile)
       write (unitnum) (tile_coord(n)%com_lat,   n=1,N_tile)
       write (unitnum) (tile_coord(n)%min_lon,   n=1,N_tile)
       write (unitnum) (tile_coord(n)%max_lon,   n=1,N_tile)
       write (unitnum) (tile_coord(n)%min_lat,   n=1,N_tile)
       write (unitnum) (tile_coord(n)%max_lat,   n=1,N_tile)
       write (unitnum) (tile_coord(n)%i_indg,    n=1,N_tile)
       write (unitnum) (tile_coord(n)%j_indg,    n=1,N_tile)
       write (unitnum) (tile_coord(n)%frac_cell, n=1,N_tile)
       write (unitnum) (tile_coord(n)%frac_pfaf, n=1,N_tile)
       write (unitnum) (tile_coord(n)%area,      n=1,N_tile)
       write (unitnum) (tile_coord(n)%elev,      n=1,N_tile)  
       
    case ('r','R')
       
       read (unitnum, iostat=istat) N_tile_tmp
       
       if (istat>0) then
          err_msg = 'ERROR reading tile_coord_file. (1)'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       if (N_tile/=N_tile_tmp) then
          err_msg = 'inconsistent N_tile'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       err_msg = 'ERROR reading tile_coord_file. (2)'
       allocate(tmp_int(N_tile))
       allocate(tmp_real(N_tile))

       read (unitnum, iostat=istat) tmp_int; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%tile_id = tmp_int(:)
       read (unitnum, iostat=istat) tmp_int; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%typ = tmp_int(:)
       read (unitnum, iostat=istat) tmp_int; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%pfaf = tmp_int(:)

       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%com_lon = tmp_real(:)
       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%com_lat= tmp_real(:)
       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%min_lon= tmp_real(:)
       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%max_lon= tmp_real(:)
       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%min_lat= tmp_real(:)
       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%max_lat= tmp_real(:)

       read (unitnum, iostat=istat) tmp_int; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%i_indg = tmp_int(:)
       read (unitnum, iostat=istat) tmp_int; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%j_indg = tmp_int(:)

       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%frac_cell = tmp_real(:)
       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%frac_pfaf = tmp_real(:)
       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%area= tmp_real(:)
       read (unitnum, iostat=istat) tmp_real; if (istat>0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       tile_coord(:)%elev= tmp_real(:)

       ! Initialize [x]_indg to pert_[x]_indg.  For cs tile spaces, pert_[x]_indg will be redefined
       tile_coord(:)%pert_i_indg = tile_coord(:)%i_indg
       tile_coord(:)%pert_j_indg = tile_coord(:)%j_indg
       tile_coord(:)%f_num = -9999 ! not assigned values yet
       deallocate(tmp_int, tmp_real)
    case default
       
       err_msg = 'unknown action ' // action
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
  end subroutine io_tile_coord_type
  
  ! *******************************************************************
  
end module LDAS_TileCoordType

! =====  EOF ==============================================================

