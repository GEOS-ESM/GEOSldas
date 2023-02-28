
! this file contains types and subroutines for tile coordinates and domain 
!
! reichle, 26 Jan 2005
! reichle,  14 Apr 2006 - split tile_coord.F90 into 2 files to avoid 
!                         having more than one module per file
! reichle,   5 Apr 2013 - added EASEv2 grid, minimal change to max lat/lon for EASE (v1)  
!
! ========================================================================

module LDAS_TileCoordRoutines

  use LDAS_TileCoordType,               ONLY:     &
       tile_coord_type,                           &
       grid_def_type,                             &
       init_grid_def_type,                        &
       io_grid_def_type,                          &
       io_tile_coord_type

  use LDAS_ensdrv_Globals,              ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_generic,                            &
       nodata_tol_generic
  
  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_RADIUS                                    ! Earth radius
  
  use EASE_conv,                        ONLY:     &
       ease_convert,                              &
       ease_inverse,                              &
       ease_extent

  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename
  
  implicit none

  private
  
  public :: get_tile_grid
  public :: get_number_of_tiles_in_cell_ij
  public :: get_tile_num_in_cell_ij
  public :: get_tile_num_in_ellipse
  public :: get_tile_num_from_latlon
  public :: get_ij_ind_from_latlon
  public :: tile2grid_simple
  public :: grid2tile
  public :: tile_mask_grid
  public :: LDAS_create_grid_g
  public :: io_domain_files

  character(10) :: tmpstring10
  character(40) :: tmpstring40

  interface grid2tile
    module procedure grid2tile_real, grid2tile_real8
  end interface grid2tile
   
contains
  
  ! **********************************************************************  
  
  subroutine io_domain_files( action, work_path, exp_id,               &
       N_cat_domain, d2g, tile_coord, tile_grid_g, tile_grid_d, rc )
    
    ! reichle, 23 July 2010
    ! reichle,  7 Jan  2014 - changed tile_coord and tile_grids I/O to binary
    ! reichle,  3 Aug  2020 - moved here from LDAS_ensdrv_init_routines.F90
    
    implicit none
    
    character,                          intent(in)    :: action
    
    character(*),                       intent(in)    :: work_path
    character(*),                       intent(in)    :: exp_id

    integer,                             intent(inout) :: N_cat_domain
    
    integer,               dimension(:), pointer       :: d2g

    type(tile_coord_type), dimension(:), pointer       :: tile_coord  ! inout

    type(grid_def_type),                 intent(inout) :: tile_grid_g
    type(grid_def_type),                 intent(inout) :: tile_grid_d

    integer,                             intent(  out) :: rc

    ! local 

    integer, parameter  :: unitnumber = 10
    
    integer             :: n, istat
    
    logical             :: writing
    
    character(300)      :: fname
    character( 40)      :: file_tag, dir_name, file_ext, tmp_action, tmp_status

    character(len=*), parameter :: Iam = 'io_domain_files'
    character(len=400) :: err_msg

    ! -----------------------------------------------------------    
    !
    ! read or write? 
    
    select case (action)
       
    case ('w','W')
       
       tmp_action = 'write'
       tmp_status = 'unknown'

       writing    = .true.

    case ('r','R')
       
       tmp_action = 'read'
       tmp_status = 'old'
       
       writing    = .false.
       
    case default
       
       err_msg = 'io_domain_files: unknown action ' // action
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
    ! -----------------------------------------------------------    
    !
    !  *.ldas_domain*txt file
    
    file_tag = 'ldas_domain'
    dir_name = 'rc_out'
    file_ext = '.txt'
 
    fname = get_io_filename( trim(work_path), trim(exp_id), file_tag, &
         dir_name=dir_name, option=1, file_ext=file_ext )
 
    if (logit) write (logunit,'(400A)') '  ' // trim(tmp_action) // ' ' // trim(fname)

    open(unitnumber, file=fname, form='formatted', action=trim(tmp_action), &
         iostat=istat, status=trim(tmp_status))
    
    if (istat/=0) then
       
       if (logit) write (logunit,*) 'cannot open for ', trim(tmp_action)
       if (logit) write (logunit,*) 

       rc = 1
       
       return
       
    else
       
       if (writing) then
          
          write (unitnumber,*) N_cat_domain
          
          do n=1,N_cat_domain
             write (unitnumber,*) tile_coord(n)%tile_id, d2g(n)
          end do
          
       else
          
          read  (unitnumber,*) N_cat_domain

          allocate(tile_coord(N_cat_domain))
          allocate(d2g(       N_cat_domain))

          do n=1,N_cat_domain
             read  (unitnumber,*) tile_coord(n)%tile_id, d2g(n)
          end do
          
       end if
       
       close (unitnumber,status='keep')
       
       if (logit) write (logunit,*) 'done with ', trim(tmp_action)
       if (logit) write (logunit,*) 
       
    end if
    
    ! -----------------------------------------------------------    
    !
    ! *.ldas_tile_coord.txt file
    
    file_tag = 'ldas_tilecoord'
    dir_name = 'rc_out'
    file_ext = '.bin'
    
    fname = get_io_filename( work_path, exp_id, file_tag, &
         dir_name=dir_name, option=1, file_ext=file_ext )
    
    if (logit) write (logunit,'(400A)') ' ' // trim(tmp_action) // ' ' // trim(fname)

    open(unitnumber, file=fname, form='unformatted', action=trim(tmp_action), &
         iostat=istat, status=trim(tmp_status))
    
    if (istat/=0) then
       
       if (logit) write (logunit,*) 'cannot open for ', trim(tmp_action)
       if (logit) write (logunit,*) 

       rc = 2
       
       return
       
    else

       call io_tile_coord_type( action, unitnumber, N_cat_domain, tile_coord )
       
       close (unitnumber,status='keep')
       
       if (logit) write (logunit,*) 'done with ', trim(tmp_action)
       if (logit) write (logunit,*) 
       
    end if

    ! -----------------------------------------------------------    
    !
    ! *.ldas_tilegrids.txt file
    
    file_tag = 'ldas_tilegrids'
    dir_name = 'rc_out'
    file_ext = '.bin'
    
    fname = get_io_filename( work_path, exp_id, file_tag, &
         dir_name=dir_name, option=1, file_ext=file_ext )

    if (logit) write (logunit,'(400A)') ' ' // trim(tmp_action) // ' ' // trim(fname)
    
    open(unitnumber, file=fname, form='unformatted', action=trim(tmp_action), &
         iostat=istat, status=trim(tmp_status))
    
    if (istat/=0) then
       
       if (logit) write (logunit,*) 'cannot open for ', trim(tmp_action)
       if (logit) write (logunit,*)        

       rc = 3
       
       return
       
    else
       
       ! read/write 'tile_grid_g'
       call io_grid_def_type( action, unitnumber, tile_grid_g )
       
       ! read/write 'tile_grid_d'
       call io_grid_def_type( action, unitnumber, tile_grid_d )
       
       close (unitnumber,status='keep')
       
       if (logit) write (logunit,*) 'done with ', trim(tmp_action)
       if (logit) write (logunit,*) 

    end if
    
    ! ------------------------------------------------------------------------
    
    rc = 0    ! successful read or write of all files
    
  end subroutine io_domain_files

  ! **********************************************************************

  subroutine LDAS_create_grid_g( gridname, n_lon, n_lat,       &
       tile_grid, i_indg_offset, j_indg_offset, cell_area)
    
    ! inputs:
    !  gridname, n_lon, n_lat
    !
    ! inouts:
    !  tile_grid   : parameters of tile definition grid  
    ! 
    ! outputs:
    !  offsets
    !  cell_area  [m^2]   (optional, for EASE  grids only)
    
    implicit none
    
    character(*),                  intent(in)    :: gridname
    integer,                       intent(in)    :: n_lon, n_lat
    type(grid_def_type),           intent(inout) :: tile_grid
    integer,                       intent(out)   :: i_indg_offset, j_indg_offset
    real,                optional, intent(out)   :: cell_area
    
    ! locals

    real    :: ease_cell_area
    logical :: date_line_on_center, pole_on_center
    logical :: ease_grid, c3_grid, latlon_grid
    logical :: file_exist
    integer :: k, rows, cols
    character(len=*),  parameter :: Iam = 'create global ldas_grid '
    character(len=400)           :: err_msg
    
    ! initialize all fields to no-data values
    
    i_indg_offset = 0
    j_indg_offset = 0
    
    call init_grid_def_type(tile_grid)
    
    tile_grid%N_lon = N_lon
    tile_grid%N_lat = N_lat
    
    tile_grid%i_offg = 0  ! tile_grid refers to *global* grid
    tile_grid%j_offg = 0  ! tile_grid refers to *global* grid
    
    date_line_on_center = .false.
    pole_on_center      = .false.
    ease_grid           = .false.
    c3_grid             = .false.
    latlon_grid         = .true.
    
    if (index(gridname,"DC") /=0) then
       date_line_on_center = .true.
    endif
    
    if (index(gridname,"PC") /=0) then
       pole_on_center      = .true.
    endif
    
    if( index(gridname,"FV") /=0 ) then
       pole_on_center      = .true.
    endif
    
    if (index(gridname,"EASEv") /=0) then
       ease_grid      = .true.
       latlon_grid    = .false.
    endif
    
    if (index(gridname,"CF") /=0) then
       c3_grid        = .true.
       latlon_grid    = .false.
    endif
    
    ! special cases , inconsistent of naming
    ! find out whether date line is on edge or through center of grid cell
    if( index(gridname,"FV_380x180") /=0) then
       pole_on_center      = .false.
    endif
    
    ! Weiyuan Note, we should fix the tile file and the naming 
    if( index(gridname,"PE_720x360_DE") /=0) then
       i_indg_offset = 110
       j_indg_offset = 230
    endif
    if( index(gridname,"PE_2880x1440_DE") /=0) then
       i_indg_offset = 440
       j_indg_offset = 920
    endif
    
    ! ----------------
    
    if (ease_grid) then
       
       ! gridname may be EASEv2-M36 or EASEv2_M36 (to be cleaned up later)

       k = index(gridname, 'EASEv')
       
       if (k == 0) then   
          err_msg = 'unknown EASE grid tile defs, gridname = ' // trim( gridname)
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       tile_grid%gridtype = trim(gridname(k:))
       
       tile_grid%ind_base = 0
       
       ! global cylindrical EASE grid 
       
       tile_grid%i_dir  = +1
       tile_grid%j_dir  = -1
       
       call ease_extent (                   &
            gridname, cols, rows,           &
            cell_area = ease_cell_area,     &             ! [m^2]
            ll_lon    = tile_grid%ll_lon,   &
            ll_lat    = tile_grid%ll_lat,   &
            ur_lon    = tile_grid%ur_lon,   &
            ur_lat    = tile_grid%ur_lat    &
            )

       tile_grid%dlon   = 360./real(tile_grid%N_lon) 
       tile_grid%dlat   = (tile_grid%ur_lat-tile_grid%ll_lat)/real(tile_grid%N_lat) ! *avg* dlat!

       if(present(cell_area)) then
          cell_area=ease_cell_area
       endif
       
    endif ! EASE grid
    
    if (latlon_grid) then !  regular LatLon grid
       
       tile_grid%gridtype = 'LatLon'
       
       tile_grid%ind_base = 1
       
       tile_grid%dlon = 360./real(tile_grid%N_lon)
       
       tile_grid%i_dir = +1
       tile_grid%j_dir = +1
       
       if (pole_on_center) then
          
          tile_grid%dlat   = 180./real(tile_grid%N_lat-1)
          
          tile_grid%ll_lat = -90. - tile_grid%dlat/2.
          tile_grid%ur_lat =  90. + tile_grid%dlat/2.
          
       else
          
          tile_grid%dlat = 180./real(tile_grid%N_lat)
          
          tile_grid%ll_lat = -90. 
          tile_grid%ur_lat =  90. 
          
       end if
       
       if (date_line_on_center) then
          
          tile_grid%ll_lon = -180. + tile_grid%dlon/2.  ! fixed sign (sqz 2023)
          tile_grid%ur_lon =  180. - tile_grid%dlon/2.  ! fixed 20 sep 2010, reichle
          
       else
          
          tile_grid%ll_lon = -180. 
          tile_grid%ur_lon =  180. 
          
       end if
       
    end if ! lat lon grid
    
    if( c3_grid) then
       
       tile_grid%gridtype = 'c3'
       tile_grid%ind_base =    1
       tile_grid%i_dir    =   +1
       tile_grid%j_dir    =   +1
       tile_grid%ll_lon   = -180. 
       tile_grid%ur_lon   =  180. 
       tile_grid%ll_lat   =  -90. 
       tile_grid%ur_lat   =   90.
       
       ! dlon and dlat are approximate!
       tile_grid%dlon     = 360./real(4*tile_grid%N_lon)
       tile_grid%dlat     = tile_grid%dlon
       
    endif
    
  end subroutine LDAS_create_grid_g
  
  ! *******************************************************************
  
  subroutine get_tile_grid( N_tile, tc_i_indg, tc_j_indg, tc_minlon, tc_minlat, tc_maxlon, tc_maxlat, &
       tile_grid_g, tile_grid )     
    
    ! get matching tile_grid for given tile_coord and (global) tile_grid_g
    !
    ! make sure to pass in consistent tile_coord (tc) and tile_grid_g inputs:
    !
    !   iff tile_grid_g is the grid that is associated with the tile space definition,
    !   then tc_[x]_indg must be tile_coord%[x]_indg
    !
    !   iff tile_grid_g is the grid that supports efficient mapping for perts and the EnKF analysis,
    !   then tc_[x]_indg must be tile_coord%hash_[x]_indg
    !
    ! reichle, 20 June 2012 -- moved from within domain_setup() 
    !                           for re-use in get_obs_pred() 
    !
    ! reichle+wjiang, 19 Aug 2020 -- changed interface to generically accommodate use of
    !                                tile_coord%[x]_indg or tile_coord%hash_[x]_indg
    !
    ! -------------------------------------------------------------------
    
    integer,               intent(in)                    :: N_tile

    integer,               intent(in), dimension(N_tile) :: tc_i_indg, tc_j_indg
    real,                  intent(in), dimension(N_tile) :: tc_minlon, tc_minlat
    real,                  intent(in), dimension(N_tile) :: tc_maxlon, tc_maxlat
    
    type(grid_def_type),   intent(in)                    :: tile_grid_g
    type(grid_def_type),   intent(out)                   :: tile_grid
    
    ! local variables
    
    integer :: n
    
    real    :: this_minlon, this_minlat, this_maxlon, this_maxlat
    
    real    :: min_min_lon, min_min_lat, max_max_lon, max_max_lat
    
    integer :: ind_i_min, ind_i_max, ind_j_min, ind_j_max
    
    integer :: this_i_indg, this_j_indg
    integer , allocatable :: i_indg_(:), j_indg_(:)
    logical :: c3_grid 
    character(len=*), parameter :: Iam = 'get_tile_grid'
    character(len=400) :: err_msg

    ! -------------------------------------------------
    
    min_min_lon =  180.                  ! initialize
    min_min_lat =   90.
    max_max_lon = -180.
    max_max_lat =  -90.
    
    ind_i_min = tile_grid_g%N_lon+1      ! initialize
    ind_j_min = tile_grid_g%N_lat+1
    ind_i_max = -1
    ind_j_max = -1

    ! THIS COMMENT SEEMS OUTDATED (reichle, 2 Aug 2020)
    ! for c3 grid, only get the ll_,ur_ lat and lon, the index is meaningless;
    ! it will be used in creating the lat_lon pert_grid
    
    if(index(tile_grid_g%gridtype,"c3") /=0) then

       ! for cube-sphere grids, do NOT zoom in 
       
       tile_grid=tile_grid_g
       
       return
       
    endif   

    if (N_tile == 0) then
      tile_grid=tile_grid_g
      tile_grid%n_lon = 0
      tile_grid%n_lat = 0
      return
    endif

    do n=1,N_tile
       
       this_minlon  = tc_minlon(n)
       this_minlat  = tc_minlat(n)
       this_maxlon  = tc_maxlon(n)
       this_maxlat  = tc_maxlat(n)
       
       this_i_indg  = tc_i_indg(n)
       this_j_indg  = tc_j_indg(n)
       
       min_min_lon = min( min_min_lon, this_minlon)
       min_min_lat = min( min_min_lat, this_minlat)
       max_max_lon = max( max_max_lon, this_maxlon)
       max_max_lat = max( max_max_lat, this_maxlat)
       
       ind_i_min   = min( ind_i_min,   this_i_indg)
       ind_j_min   = min( ind_j_min,   this_j_indg)
       ind_i_max   = max( ind_i_max,   this_i_indg)
       ind_j_max   = max( ind_j_max,   this_j_indg)
       
    end do
    
    ! assemble tile_grid (revised 20 Sep 2010, reichle)
       
    tile_grid%N_lon  = ind_i_max - ind_i_min + 1
    tile_grid%N_lat  = ind_j_max - ind_j_min + 1
    
    tile_grid%i_offg = ind_i_min - tile_grid_g%ind_base
    tile_grid%j_offg = ind_j_min - tile_grid_g%ind_base
       
    tile_grid%dlon   = tile_grid_g%dlon
       
    if     (index(tile_grid_g%gridtype,'LatLon')/=0) then
       
       tile_grid%dlat   = tile_grid_g%dlat
       
       tile_grid%ll_lon = tile_grid_g%ll_lon + real(tile_grid%i_offg)*tile_grid%dlon
       tile_grid%ll_lat = tile_grid_g%ll_lat + real(tile_grid%j_offg)*tile_grid%dlat
       
       tile_grid%ur_lon = tile_grid%ll_lon + real(tile_grid%N_lon)*tile_grid%dlon
       tile_grid%ur_lat = tile_grid%ll_lat + real(tile_grid%N_lat)*tile_grid%dlat
       
    elseif ( index(tile_grid_g%gridtype,'EASEv')  /=0 ) then
       
       ! *average* dlat over the domain
       
       tile_grid%dlat   = (max_max_lat - min_min_lat)/real(tile_grid%N_lat)
       
       tile_grid%ll_lon = min_min_lon
       tile_grid%ll_lat = min_min_lat
       
       tile_grid%ur_lon = max_max_lon
       tile_grid%ur_lat = max_max_lat
       
    else
       
       err_msg = 'not yet implemented for ' // tile_grid_g%gridtype
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    tile_grid%gridtype = tile_grid_g%gridtype
    
    tile_grid%ind_base = tile_grid_g%ind_base   ! should this be inherited???
    
    tile_grid%i_dir    = tile_grid_g%i_dir
    tile_grid%j_dir    = tile_grid_g%j_dir
    
  end subroutine get_tile_grid

  ! **********************************************************************

  subroutine get_number_of_tiles_in_cell_ij( N_tile, tc_i_indg, tc_j_indg, tile_grid, &
       N_tile_in_cell_ij)
    
    ! find out how many tiles are in a given tile definition grid cell
    ! reichle, 22 Jul 2005
    !
    ! make sure to pass in consistent tile_coord (tc) and tile_grid_g inputs:
    !
    !   iff tile_grid_g is the grid that is associated with the tile space definition,
    !   then tc_[x]_indg must be tile_coord%[x]_indg
    !
    !   iff tile_grid_g is the grid that supports efficient mapping for perts and the EnKF analysis,
    !   then tc_[x]_indg must be tile_coord%hash_[x]_indg
    !
    ! wjiang(?) -- split off from LDASsa legacy subroutine get_tile_num_in_cell_ij()
    !
    ! reichle+wjiang, 19 Aug 2020 -- changed interface to generically accommodate use of
    !                                tile_coord%[x]_indg or tile_coord%hash_[x]_indg
    !
    ! ----------------------------------------------------------
    
    implicit none
    
    integer,             intent(in)                    :: N_tile
    
    integer,             intent(in), dimension(N_tile) :: tc_i_indg, tc_j_indg
    
    type(grid_def_type), intent(in)                    :: tile_grid
        
    integer, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(inout) :: &
         N_tile_in_cell_ij
    
    ! locals 
    
    integer :: i, j, k, n, off_i, off_j
    

    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    ! (re-)initialize
    
    N_tile_in_cell_ij = 0
    
    do n=1,N_tile
       
       i = tc_i_indg(n) - off_i
       j = tc_j_indg(n) - off_j
       
       N_tile_in_cell_ij(i,j) = N_tile_in_cell_ij(i,j) + 1
       
    enddo

    if (logit) write (logunit,*) &
        'Maximum number of tiles in tile def grid cell = ', maxval(N_tile_in_cell_ij)
    if (logit) write (logunit,*)

  end subroutine get_number_of_tiles_in_cell_ij

  ! **********************************************************************

  subroutine get_tile_num_in_cell_ij( N_tile, tc_i_indg, tc_j_indg, tile_grid, &
       max_N_tile_in_cells, tile_num_in_cell_ij )
    
    ! find out tile_num in given cells
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domain of interest.  
    ! With these offsets tile2grid() can be used to map from a
    ! subgrid of "tile_grid_g" to tile space
    !
    ! reichle, 22 Jul 2005
    !
    ! make sure to pass in consistent tile_coord (tc) and tile_grid_g inputs:
    !
    !   iff tile_grid_g is the grid that is associated with the tile space definition,
    !   then tc_[x]_indg must be tile_coord%[x]_indg
    !
    !   iff tile_grid_g is the grid that supports efficient mapping for perts and the EnKF analysis,
    !   then tc_[x]_indg must be tile_coord%hash_[x]_indg
    !
    ! wjiang(?) -- split off from LDASsa legacy subroutine get_tile_num_in_cell_ij()
    !
    ! reichle+wjiang, 19 Aug 2020 -- changed interface to generically accommodate use of
    !                                tile_coord%[x]_indg or tile_coord%hash_[x]_indg
    !
    ! ----------------------------------------------------------
    
    implicit none
    
    integer,             intent(in)                    :: N_tile
    
    integer,             intent(in), dimension(N_tile) :: tc_i_indg, tc_j_indg
    
    type(grid_def_type), intent(in)                    :: tile_grid
        
    integer,             intent(in)                    ::  max_N_tile_in_cells
    
    ! the pointer is an output arguments that is allocated here
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij

    integer, dimension(tile_grid%N_lon,tile_grid%N_lat) :: &
         N_tile_in_cell_ij

    ! locals 
    
    integer, parameter :: nodata = -9999.
    
    integer :: i, j, k, n, off_i, off_j
    
    ! -----------------------------------------------------------------
    !
    ! allocate and initialize pointers if present
    
       
    allocate(tile_num_in_cell_ij(tile_grid%N_lon,tile_grid%N_lat,    &
          max_N_tile_in_cells))
       
    tile_num_in_cell_ij = nodata
       
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    ! (re-)initialize
    
    N_tile_in_cell_ij = 0
    
    do n=1,N_tile
       
       i = tc_i_indg(n) - off_i
       j = tc_j_indg(n) - off_j
       
       N_tile_in_cell_ij(i,j) = N_tile_in_cell_ij(i,j) + 1
       
       k = N_tile_in_cell_ij(i,j)
          
       tile_num_in_cell_ij(i,j,k) = n 
          
    end do
    
  end subroutine get_tile_num_in_cell_ij
  
  ! *******************************************************************
  
  subroutine get_tile_num_from_latlon(N_catd, tile_coord,                      &
       tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij, N_latlon, lat, lon,  &
       tile_num, max_dist_x, max_dist_y )

    ! bug fix re. "check that lat/lon is inside tile_grid"
    ! - reichle, 2005/11/17
    !
    ! added optional input "max_dist" that permits returning a valid tile_num
    ! even if lat/lon is not within the bounding box of the chosen tile
    ! - reichle, 2008/03/28
    !
    ! distance units are in [deg]; added functionality for latitude-dependent
    !  meridional max distance (max_dist_x)
    ! - reichle, 2014/12/23
    
    implicit none
    
    integer,                                             intent(in)  :: N_catd, N_latlon
    
    type(tile_coord_type), dimension(:),                 pointer     :: tile_coord ! input
    
    type(grid_def_type),                                 intent(in)  :: tile_grid    
    
    integer, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in)  :: N_tile_in_cell_ij
    
    integer, dimension(:,:,:),                           pointer     :: tile_num_in_cell_ij ! input    
    
    real,    dimension(N_latlon),                        intent(in)  :: lat, lon
    
    integer, dimension(N_latlon),                        intent(out) :: tile_num
    
    real,    dimension(N_latlon),              optional, intent(in)  :: max_dist_x  ! vector [deg]
    real,                                      optional, intent(in)  :: max_dist_y  ! scalar [deg]
    

    ! local variables
    
    integer, parameter           :: nodata_tilenum = -9999
    
    integer                      :: n, k, i_ind, j_ind, this_tile_num

    logical                      :: outside_bbox, too_far_away
    
    real                         :: tmp_dist,   min_dist
    real                         :: tmp_dist_y, min_dist_y, max_dist_tmp_y
    real                         :: tmp_dist_x, min_dist_x
    real, dimension(N_latlon)    ::                         max_dist_tmp_x 
    
    character(len=*),  parameter :: Iam = 'get_tile_num_from_latlon'
    character(len=400)           :: err_msg

    ! -----------------------------------------------------------
    
    ! initialize
    
    max_dist_tmp_x = 0.
    max_dist_tmp_y = 0.

    if (present(max_dist_x))  max_dist_tmp_x = max_dist_x
    if (present(max_dist_y))  max_dist_tmp_y = max_dist_y
    
    tile_num = nodata_tilenum          ! initialize to negative value
    
    ! loop through observations
    
    do n=1,N_latlon
       
       ! Make sure lat/lon is *inside* tile_grid, otherwise do nothing.
       ! Do *not* allow obs outside tile_grid because obs perturbations 
       ! routines are not set up to generate perturbations for such obs, 
       ! and obs_pred would be questionable anyway.
       ! - reichle+csdraper, 29 Jan 2016

       if ( tile_grid%ll_lat   < (lat(n)          )  .and.  &
            tile_grid%ll_lon   < (lon(n)          )  .and.  &
            (lat(n)          ) < tile_grid%ur_lat    .and.  &
            (lon(n)          ) < tile_grid%ur_lon              ) then  

          ! min_dist = distance betw lat/lon in question and center-of-mass of
          !            matching tile 
          
          min_dist   = 1.e10    ! initialize   (bug fix, csdraper+reichle, 30 Jun 2015)
          min_dist_x = 1.e10    ! initialize 
          min_dist_y = 1.e10    ! initialize 
          
          ! determine grid cell that contains lat/lon 
          
          call get_ij_ind_from_latlon( tile_grid, lat(n), lon(n), i_ind, j_ind )
          
          ! make sure that i/j_ind is still within bounds 
          ! (works in conjunction with if statement above re. ll/ur_lat/lon)
          
          i_ind = min( max(i_ind, 1), tile_grid%N_lon )
          j_ind = min( max(j_ind, 1), tile_grid%N_lat )
          
          ! map from i_ind, j_ind to tile_num
          
          if   ( index(tile_grid%gridtype, 'EASEv')  /=0 ) then
             
             ! ASSUMPTION: tiles match EASE grid cells exactly
             !             (unless "outside" the domain, eg. water surface)
             
             if     (N_tile_in_cell_ij(i_ind,j_ind)==1) then
                
                tile_num(n)=tile_num_in_cell_ij(i_ind,j_ind,1)

                min_dist_x = abs(lon(n) - tile_coord(tile_num(n))%com_lon)
                min_dist_y = abs(lat(n) - tile_coord(tile_num(n))%com_lat) 
                
             elseif (N_tile_in_cell_ij(i_ind,j_ind)==0) then
                
                ! Do nothing.  If given EASE grid cell is not land, 
                ! tile_num will not change from its initialized value.
                
             else
                
                err_msg = 'something wrong for EASE grid'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                
             end if
             
          elseif (index(tile_grid%gridtype, 'LatLon')/=0) then
             
             ! Loop through all tiles within grid cell that contain lat/lon and
             ! find minimum distance (Minkowski norm).
             ! If there are no land tiles in given tile definition grid cell, 
             ! tile_num will not change from its initialized value.
             
             do k=1,N_tile_in_cell_ij(i_ind,j_ind) 
                
                this_tile_num = tile_num_in_cell_ij(i_ind,j_ind,k)

                tmp_dist_x = abs(lon(n) - tile_coord(this_tile_num)%com_lon)
                tmp_dist_y = abs(lat(n) - tile_coord(this_tile_num)%com_lat) 

                tmp_dist   = tmp_dist_x + tmp_dist_y    ! Minkowski norm
                
                if (tmp_dist<min_dist) then
                   
                   min_dist    = tmp_dist
                   min_dist_x  = tmp_dist_x
                   min_dist_y  = tmp_dist_y
                   
                   tile_num(n) = this_tile_num
                   
                end if
                
             end do

          else
             
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown grid type')
             
          end if
                    

          ! check whether lat/lon is inside bounding box of given tile
          ! (if domain consists of a single tile and tile_grid is coarse,
          !  it is likely that lat/lon is outside bounding box of tile
          !  even if the lat/lon is inside the tile definition grid cell that
          !  contains the given tile) - reichle, 18 Aug 2005
          
          ! return a valid tile_num even if lat/lon is not within the 
          ! bounding box of the chosen tile, as long as lat/lon is no further
          ! from the com_lat/com_lon of the tile than max_dist
          ! - reichle, 2008/03/28
          
          ! this might be important for assimilation of 36km EASE grid obs 
          ! into the 9 km EASE grid tile space (because the lat/lon of
          ! the 36 km EASE grid obs is exactly on the corners of four
          ! 9 km EASE grid cells)
          ! - reichle, 2013/06/21          

          if (tile_num(n)>0) then

             outside_bbox  = (                                          &
                  lon(n) < tile_coord(tile_num(n))%min_lon    .or.      &
                  lon(n) > tile_coord(tile_num(n))%max_lon    .or.      &
                  lat(n) < tile_coord(tile_num(n))%min_lat    .or.      &
                  lat(n) > tile_coord(tile_num(n))%max_lat           )      
             
             too_far_away = (                                           &
                  min_dist_x > max_dist_tmp_x(n)              .or.      &
                  min_dist_y > max_dist_tmp_y                        )
             
             ! keep tile_num unless obs is outside the bounding box *and* too far away
             
             if (outside_bbox .and. too_far_away)  tile_num(n) = nodata_tilenum
             
          end if
          
       end if
       
    end do
        
  end subroutine get_tile_num_from_latlon
  
  ! *******************************************************************
  
  subroutine get_ij_ind_from_latlon( tile_grid, lat, lon, i_ind, j_ind )

    ! NOTE order of input arguments ("lat", "lon", "lon_ind", "lat_ind")

    ! find i/j_ind of grid cell that contains lat/lon
    !
    ! that is, for given lat/lon, find out corresponding i/j_ind into an array
    ! of size N_lon-by-N_lat whose lat/lon coordinates are defined by tile_grid
    
    ! NOTE: ALL arrays in LDASsa are declared with base 1
    
    ! typical use: have lat/lon of an obs, want to know i/j_ind within the 
    ! grid (array) that defines the tile space ("tile_grid_d")

    ! major revision by reichle, 11 May 2011
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    real, intent(in) :: lat, lon
    integer, intent(out) :: i_ind, j_ind
    
    real    :: lats(1),    lons(1)
    integer :: i_inds(1), j_inds(1)

    ! local variables
    
    real :: r, s, i_indg, j_indg

    character(len=*), parameter :: Iam = 'get_ij_ind_from_latlon'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------
    !
    ! select grid type
    
    if     (index(tile_grid%gridtype, 'EASEv')/=0) then

       ! EASE grid lat/lon to index provides *global*, *0-based* index!

       call ease_convert(tile_grid%gridtype, lat, lon, r, s)
       
       i_indg = nint(r)    ! i_ind or lon_ind
       j_indg = nint(s)    ! j_ind or lat_ind

       ! convert to index into array defined by tile_grid_d
       
       i_ind = i_indg - tile_grid%i_offg - (tile_grid%ind_base - 1)
       j_ind = j_indg - tile_grid%j_offg - (tile_grid%ind_base - 1)
       
    elseif (index(tile_grid%gridtype, 'LatLon')/=0) then
       
       ! ll_lon and ll_lat refer to lower left  corner of grid cell
       ! ur_lon and ur_lat refer to upper right corner of grid cell
       ! (as opposed to the grid point in the center of the grid cell)
       !
       ! ALL arrays in LDASsa are declared with base 1 (--> use "ceiling")

       if (tile_grid%i_dir==1) then
          
          i_ind = ceiling( (lon - tile_grid%ll_lon)/tile_grid%dlon )
          
       else
          
          i_ind = ceiling( (tile_grid%ur_lon - lon)/tile_grid%dlon )

       end if


       if (tile_grid%j_dir==1) then
          
          j_ind = ceiling( (lat - tile_grid%ll_lat)/tile_grid%dlat )
          
       else
          
          j_ind = ceiling( (tile_grid%ur_lat - lat)/tile_grid%dlat )
          
       end if

    else

       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown grid type')

    end if
    
  end subroutine get_ij_ind_from_latlon
  
  ! *******************************************************************

  subroutine get_tile_num_in_ellipse( lon, lat, r_x, r_y,                 &
       N_tile, tile_coord, tile_grid,                                     &
       N_tile_in_cell_ij, tile_num_in_cell_ij,                            &
       N_tile_in_ellipse,  tile_num_in_ellipse,  norm_square_distance   )
            
    ! reichle, 22 Feb 2015

    ! return tile numbers of all tiles within ellipse with center at lat/lon
    !  and axes r_x, r_y
    !
    ! also return normalized square distance of tile center-of-mass lat/lon from lat/lon

    ! TO DO: what about dateline? poles?  [ignore for now]
    
    ! NOTE: get_tile_num_in_ellipse() returns zero tiles (N_tile_in_ellipse=0)
    !        if ellipse straddles date line and EASE[v2] grid is used
    !       - reichle, 17 Apr 2017
    
    
    implicit none
    
    real,                                     intent(in)    :: lon, lat, r_x, r_y

    integer,                                  intent(in)    :: N_tile
    
    type(tile_coord_type), dimension(:),      pointer       :: tile_coord ! input
    
    type(grid_def_type),                      intent(in)    :: tile_grid
    
    integer,               dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in) :: &
         N_tile_in_cell_ij
        
    integer,               dimension(:,:,:),  pointer       :: tile_num_in_cell_ij ! input
    
    integer,                                  intent(out)   :: N_tile_in_ellipse
        
    integer,               dimension(N_tile), intent(out)   :: tile_num_in_ellipse

    real,                  dimension(N_tile), intent(out)   :: norm_square_distance
    
    ! local variables
    
    real    :: minlat, maxlat, minlon, maxlon, tmp_dx, tmp_dy, tmp_dist
    real    :: r_x_square, r_y_square
    
    integer :: i_ind_ll, i_ind_ur, i_min, i_max
    integer :: j_ind_ll, j_ind_ur, j_min, j_max

    integer :: ii, jj, kk, mm, nn
    
    ! -------------------------------------------------------------
    !
    ! initialize
    
    N_tile_in_ellipse   = 0
    
    tile_num_in_ellipse = nodata_generic
    
    r_x_square          = r_x**2
    r_y_square          = r_y**2
    
    ! identify rectangle +/- r_x and +/- r_y around lat/lon
    
    tmp_dx = tile_grid%dlon  ! broaden rectangle to avoid missing tiles
    tmp_dy = tile_grid%dlon  ! use *dlon* because for EASE grid *dlat* is ill-defined   
    
    minlon = lon - r_x - tmp_dx
    maxlon = lon + r_x + tmp_dx
    minlat = lat - r_y - tmp_dy
    maxlat = lat + r_y + tmp_dy
    
    ! find i,j indices of two opposite corners (lower left and upper right)
    
    call get_ij_ind_from_latlon( tile_grid, minlat, minlon, i_ind_ll, j_ind_ll )
    call get_ij_ind_from_latlon( tile_grid, maxlat, maxlon, i_ind_ur, j_ind_ur )
    
    ! restrict i_ind_* and j_ind_* to fall within tile_grid
    
    i_ind_ll = min( tile_grid%N_lon, max( 1, i_ind_ll ))
    j_ind_ll = min( tile_grid%N_lat, max( 1, j_ind_ll ))
    i_ind_ur = min( tile_grid%N_lon, max( 1, i_ind_ur ))
    j_ind_ur = min( tile_grid%N_lat, max( 1, j_ind_ur ))

    ! find smaller and larger of the two indices
    ! (index may run north-to-south/west-to-east depending on i_dir/j_dir)
    
    i_min = min( i_ind_ll, i_ind_ur )
    i_max = max( i_ind_ll, i_ind_ur )

    j_min = min( j_ind_ll, j_ind_ur )
    j_max = max( j_ind_ll, j_ind_ur )
    
    ! loop through all grid cells "inside" rectangle and identify tiles w/in ellipse

    norm_square_distance = nodata_generic

    mm = 0
    
    do ii=i_min,i_max
       do jj=j_min,j_max
          
          do kk=1,N_tile_in_cell_ij(ii,jj)
             
             nn = tile_num_in_cell_ij(ii,jj,kk)
             
             tmp_dist = (                                            &
                  (tile_coord(nn)%com_lon - lon)**2 / r_x_square +    &
                  (tile_coord(nn)%com_lat - lat)**2 / r_y_square   )
             
             if (tmp_dist <= 1.) then
                
                ! record tile number and distance
                
                mm = mm+1
                
                tile_num_in_ellipse(mm) = nn
                
                norm_square_distance(mm) = tmp_dist
                
             end if
                          
          end do
          
       end do
    end do
    
    N_tile_in_ellipse = mm
    
  end subroutine get_tile_num_in_ellipse

  ! **********************************************************************
  
  subroutine tile2grid_simple( N_tile, i_indgs, j_indgs, tile_grid, tile_data, grid_data)

    ! Interpolate from tile space to grid space without interpolation/weighted/no-data-handling.
    ! Simply assign the tile value to the grid cell (last assignment prevails)

    implicit none
    
    integer, intent(in) :: N_tile

    integer, intent(in) :: i_indgs(:)     ! dimension(N_tile)
    integer, intent(in) :: j_indgs(:)     ! dimension(N_tile)
     
    type(grid_def_type), intent(in) :: tile_grid
    
    real, dimension(N_tile), intent(in) :: tile_data
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(out) :: grid_data
        
    ! local variables
    
    integer :: n, i, j, off_i, off_j
    real, parameter :: no_data= -9999.   
    character(len=*), parameter :: Iam = 'tile2grid_simple'
    character(len=400) :: err_msg

    ! ------------------------------------
    
    if (size(i_indgs)/=N_tile .or. size(j_indgs)/=N_tile) then
       err_msg = '[i,j]_indg and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    !
    ! adjust for 0-based indexing (eg., EASE grids)
    !
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    ! loop through tile space

    grid_data = no_data
    do n=1,N_tile
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j
       grid_data(i,j) = tile_data(n)
    end do
    
  end subroutine tile2grid_simple

  ! *******************************************************************
  !
  ! The subroutine tile2grid_full() is the original LDASsa tile2grid() 
  !  subroutine with weighted averaging and no-data-handling.
  !
  ! This subroutine should no longer be used and is provided here for
  !   reference only.  Use MAPL LocationStream instead.
  !
  ! - reichle, 13 July 2020

  subroutine tile2grid_full( N_tile, tile_coord, tile_grid, tile_data,    &
       grid_data, no_data_value, no_data_tol, echo )
    
    ! map from tile space to tile definition grid
    !
    ! NOTE: tile_coord must match tile_data
    !
    ! optional inputs: 
    !   no_data_value :
    !   no_data_tol   : tolerance when checking tile_data 
    !                     against no_data_value
    !   echo          : echo no_data_value and tolerance
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domaim of interest.  
    ! With these offsets tile2grid() can be used to map from tile space to a 
    ! subgrid of "tile_grid_g"
    
    ! reichle, 28 Jan 2005
    ! reichle, 16 May 2005 - added offset for "domain" grid
    ! reichle, 21 May 2010 - off_i, off_j now part of grid_def_type

    implicit none
    
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid
    
    real, dimension(N_tile), intent(in) :: tile_data
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(out) :: grid_data
    
    real, intent(in), optional :: no_data_value, no_data_tol
    
    logical, intent(in), optional :: echo
        
    ! local variables
    
    integer :: n, i, j, off_i, off_j
    
    real :: w, no_data, tol
    
    real, parameter :: no_data_default     = -9999.
    real, parameter :: no_data_tol_default = 1e-4
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat) :: wgrid
    
    character(len=*), parameter :: Iam = 'tile2grid'
    character(len=400) :: err_msg

    ! ------------------------------------
    
    if (size(tile_coord)/=N_tile) then
       err_msg = 'tile_coord and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    if (present(no_data_value)) then
       no_data = no_data_value
    else
       no_data = no_data_default
    end if
    
    if (present(no_data_tol)) then
       tol = no_data_tol
    else
       tol = no_data_tol_default
    end if
    
    if (present(echo)) then
       if (echo .and. logit) then
          write (logunit,*) 'tile2grid: using no-data-value = ', no_data , &
               ' with tolerance = ', tol
       end if
    end if
    
    ! ------------------------------------------------------
    !
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    ! initialize
    
    grid_data = 0.
    wgrid     = 0.
    
    ! loop through tile space
    
    do n=1,N_tile
       
       i = tile_coord(n)%i_indg - off_i
       j = tile_coord(n)%j_indg - off_j
       
       w = tile_coord(n)%frac_cell
       
       if (abs(tile_data(n)-no_data)>tol) then
          
          grid_data(i,j) = grid_data(i,j) + w*tile_data(n)
          
          wgrid(i,j) = wgrid(i,j) + w
          
       end if
       
    end do
    
    ! normalize and set no-data-value
    
    do i=1,tile_grid%N_lon
       do j=1,tile_grid%N_lat
          
          if (wgrid(i,j)>0.) then
             
             grid_data(i,j) = grid_data(i,j)/wgrid(i,j)
             
          else
             
             grid_data(i,j) = no_data
             
          end if
          
       end do
    end do
    
  end subroutine tile2grid_full
  
  ! *******************************************************************

  subroutine tile_mask_grid( tile_grid, N_tile, i_indgs,j_indgs, grid_data)
    
    ! set grid cell to no value if there is no tile in it
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    integer, intent(in) :: N_tile
    
    integer, intent(in) :: i_indgs(:)
    integer, intent(in) :: j_indgs(:)
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(inout) :: grid_data
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat) :: grid
        
    ! local variables
    real, parameter :: no_data= -9999.
    integer :: n, i, j, off_i, off_j

    character(len=*), parameter :: Iam = 'tile_mask_grid'
    character(len=400) :: err_msg
    
    ! ------------------------------------
    
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
       
    grid = no_data
    do n=1,N_tile
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j
       grid(i,j) = grid_data(i,j)
    end do
    grid_data = grid

  end subroutine tile_mask_grid

  ! **********************************************************************

  subroutine grid2tile_real( tile_grid, N_tile, i_indgs,j_indgs, grid_data, tile_data)
    
    ! map from grid to tile space
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domain of interest.  
    ! With these offsets grid2tile() can be used to map from a
    ! subgrid of "tile_grid_g" to tile space
    !
    ! reichle, 28 Jan 2005
    ! reichle, 16 May 2005 - added offset for "domain" grid
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    integer, intent(in) :: N_tile
    
    integer, intent(in) :: i_indgs(:)
    integer, intent(in) :: j_indgs(:)
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in) :: grid_data
    
    real, dimension(N_tile), intent(out) :: tile_data
        
    ! local variables
    
    integer :: n, i, j, off_i, off_j

    character(len=*), parameter :: Iam = 'grid2tile_real'
    character(len=400) :: err_msg
    
    ! ------------------------------------
    
    if (size(i_indgs)/=N_tile .or. size(j_indgs)/=N_tile) then
       err_msg = '[i,j]_indg and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
        
    do n=1,N_tile
       
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j

       tile_data(n) = grid_data(i,j)
       
    end do
    
  end subroutine grid2tile_real

  ! **********************************************************************
  
  subroutine grid2tile_real8( tile_grid, N_tile, i_indgs,j_indgs, grid_data_8, tile_data_8)
    
    ! same as grid2tile_real but for real*8
    !
    ! reichle,  3 Feb 2014
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    integer, intent(in) :: N_tile
    integer, intent(in) :: i_indgs(:)
    integer, intent(in) :: j_indgs(:)
    
    real*8, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in) :: grid_data_8
    
    real*8, dimension(N_tile), intent(out) :: tile_data_8
    
    ! local variables
    
    integer :: n, i, j, off_i, off_j
    
    character(len=*), parameter :: Iam = 'grid2tile_real8'
    character(len=400) :: err_msg

    ! ------------------------------------
    
    if (size(i_indgs)/=N_tile .or. size(j_indgs)/=N_tile) then
       err_msg = '[i,j]_indg and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
        
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
        
    do n=1,N_tile
       
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j
       
       tile_data_8(n) = grid_data_8(i,j)
       
    end do
    
  end subroutine grid2tile_real8

  ! **********************************************************************

end module LDAS_TileCoordRoutines


! ================================= EOF =======================================
