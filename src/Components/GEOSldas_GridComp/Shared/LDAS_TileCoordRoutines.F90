
! this file contains types and subroutines for tile coordinates and domain 
!
! reichle, 26 Jan 2005
! reichle,  14 Apr 2006 - split tile_coord.F90 into 2 files to avoid 
!                         having more than one module per file
! reichle,   5 Apr 2013 - added EASEv2 grid, minimal change to max lat/lon for EASE (v1)  
!
! ========================================================================

module LDAS_TileCoordRoutines

  use MAPL_BaseMod, only: MAPL_GetHorzIJIndex
  use LDAS_TileCoordType,                 ONLY:     &
       tile_coord_type,                           &
       grid_def_type,                             &
       init_grid_def_type,                        &
       io_grid_def_type,                          &
       io_tile_coord_type,                        &
       tile_typ_land,                             &
       N_cont_max

  use LDAS_ensdrv_Globals,           ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_generic,                            &
       nodata_tol_generic
  
  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_RADIUS                                    ! Earth radius
  
  use LDAS_EASE_conv,                      ONLY:     &
       easev1_convert,easev2_convert

  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

 use LDAS_ensdrv_functions,            ONLY:     &
       is_in_list,                                &
       is_in_domain 
  
  implicit none

  private
  
  public :: read_til_file
  public :: get_tile_grid
  public :: get_land_tile_info
  public :: get_number_of_tiles_in_cell_ij
  public :: get_tile_num_in_cell_ij
  public :: get_tile_num_in_ellipse
  public :: get_tile_num_from_latlon
  public :: get_ij_ind_from_latlon
  public :: tile2grid
  public :: tile_mask_grid
  public :: grid2tile, grid2tile_real8
  public :: is_cat_in_box
  public :: reorder_tiles
  public :: get_is_land
  public :: LDAS_create_grid_g
  public :: LDAS_read_land_tile
  character(10) :: tmpstring10
  character(40) :: tmpstring40

  interface grid2tile
    module procedure grid2tile_new, grid2tile_old
  end interface  
contains

  subroutine LDAS_create_grid_g( gridname,n_lon,n_lat, tile_grid,i_indg_offset,j_indg_offset,cell_area)
    
    ! inputs:
    !  grid name, n_lon, n_lat
    ! outputs:
    !
    !  tile_grid   : parameters of tile definition grid
    !  offsets
    implicit none
    
    character(*),intent(in) :: gridname
    integer,intent(in) :: n_lon,n_lat
    type(grid_def_type), intent(inout) :: tile_grid
    integer,intent(out) :: i_indg_offset, j_indg_offset
    real,optional,intent(out)    :: cell_area
    
    ! locals

    real    :: ease_cell_area
    logical :: date_line_on_center, pole_on_center
    logical :: ease_grid,c3_grid,latlon_grid
    logical :: file_exist

    character(len=*), parameter :: Iam = 'create global ldas_grid '
    character(len=400) :: err_msg


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

     if (index(gridname,"EASE") /=0) then
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

          tile_grid%ind_base = 0
             
        ! global cylindrical EASE grid 
             
          tile_grid%ll_lon = -180.
          tile_grid%ur_lon =  180.
             
          tile_grid%i_dir  = +1
          tile_grid%j_dir  = -1
          ! It is sloopy that the name may be EASEv2-M36 or EASEv2_M36
          if (index(gridname, 'EASEv2_M36')/=0 .or. index(gridname, 'EASEv2-M36')/=0) then  ! version *2*
                
             tile_grid%gridtype = 'EASEv2_M36'                          
                
             tile_grid%ll_lat =  -85.04456
             tile_grid%ur_lat =   85.04456

             ease_cell_area   = 1298.320938704616
              
          elseif (index(gridname, 'EASEv2_M09')/=0 .or. index(gridname, 'EASEv2-M09')/=0) then  ! version *2*

             tile_grid%gridtype = 'EASEv2_M09'  

             tile_grid%ll_lat =  -85.04456
             tile_grid%ur_lat =   85.04456
               
             ease_cell_area   =   81.145058669038477

          elseif (index(gridname, 'EASE_M36')/=0 .or. index(gridname, 'EASE-M36')/=0) then
                
             tile_grid%gridtype = 'EASE_M36'                          
                
             tile_grid%ll_lat =  -86.62256 ! minimal change, reichle, 5 Apr 2013
             tile_grid%ur_lat =   86.62256 ! minimal change, reichle, 5 Apr 2013

             ease_cell_area   = 1296.029001087600 
                
          elseif (index(gridname, 'EASE_M09')/=0 .or. index(gridname, 'EASE-M09')/=0) then

             tile_grid%gridtype = 'EASE_M09'                          

             tile_grid%ll_lat =  -86.62256 ! minimal change, reichle, 5 Apr 2013
             tile_grid%ur_lat =   86.62256 ! minimal change, reichle, 5 Apr 2013
                
             ease_cell_area   =   81.001812568020028

          elseif (index(gridname, 'EASE_M25')/=0 .or. index(gridname, 'EASE-M25')/=0 ) then

             tile_grid%gridtype = 'EASE_M25'                          

             tile_grid%ll_lat =  -86.7167 ! need to double-check (reichle, 11 May 2011)
             tile_grid%ur_lat =   86.7167 ! need to double-check (reichle, 11 May 2011)
                
             ease_cell_area   =   628.38080962
                
          else
               
              err_msg = 'unknown EASE grid tile defs, grid name = ' &
                     // trim( gridname)
              call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                
          end if
             
          tile_grid%dlon   = 360./real(tile_grid%N_lon) 
          tile_grid%dlat   = &
                (tile_grid%ur_lat-tile_grid%ll_lat)/real(tile_grid%N_lat) ! *avg* dlat!
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
             
             tile_grid%ll_lon = -180. - tile_grid%dlon/2.
             tile_grid%ur_lon =  180. - tile_grid%dlon/2.  ! fixed 20 sep 2010, reichle
             
         else
             
             tile_grid%ll_lon = -180. 
             tile_grid%ur_lon =  180. 
          
         end if
          
       end if ! lat lon grid

       if( c3_grid) then

          tile_grid%gridtype ='c3'
          tile_grid%ind_base = 1
          tile_grid%i_dir = +1
          tile_grid%j_dir = +1
          tile_grid%ll_lon = -180. 
          tile_grid%ur_lon =  180. 
          tile_grid%ll_lat = -90. 
          tile_grid%ur_lat =  90. 
          ! dlon and dlat are approximate
          tile_grid%dlon = 360./real(4*tile_grid%N_lon)
          tile_grid%dlat = tile_grid%dlon
          
       endif

  end subroutine LDAS_create_grid_g

  subroutine LDAS_read_land_tile( tile_file,catch_file, tile_grid_g, tile_coord_land,f2g )
    ! inputs:
    !
    !  tile_file : full path + name 
    !
    ! outputs:
    !
    !  tile_grid   : parameters of tile definition grid
    !
    !  tile_coord : coordinates of tiles (see tile_coord_type),
    !               implemented as pointer which is allocated in 
    !               this subroutine
    !               NOTE: number of tiles can be diagnosed 
    !                     with size(tile_coord)
    ! optional:
    !    if the tile file type 1100 ,which is land excluded
    !    f2g  : the full domain id to the global id
    ! "tile_id" is no longer read from *.til file and is now set in this 
    ! subroutine to match order of tiles in *.til file
    ! - reichle, 22 Aug 2013
    !
    ! -------------------------------------------------------------

    implicit none

    character(*), intent(in) :: tile_file
    character(*), intent(in) :: catch_file
    type(grid_def_type), intent(inout):: tile_grid_g
    type(tile_coord_type), dimension(:), pointer :: tile_coord_land ! out
    integer, dimension(:), optional,pointer :: f2g ! out

    ! locals
    type(tile_coord_type),dimension(:),allocatable :: tile_coord
    integer, dimension(:), allocatable :: f2g_tmp  ! out

    real    :: ease_cell_area
    integer :: i, N_tile,N_grid,tmpint1, tmpint2, tmpint3, tmpint4
    integer :: i_indg_offset, j_indg_offset, col_order
    integer :: N_tile_land,n_lon,n_lat

    logical ::  ease_grid
    integer :: typ,k,fid
    character(200) :: tmpline,gridname
    character(300) :: fname

    character(len=*), parameter :: Iam = 'LDAS_read_tile_file'
    character(len=400) :: err_msg

    ! ---------------------------------------------------------------

    i_indg_offset = 0
    j_indg_offset = 0

   ! call LDAS_create_grid(tile_file,tile_grid_g,i_indg_offset,j_indg_offset)

   ! read file header 

    if (logit) write (logunit,'(400A)') 'LDAS_read_tile_file(): reading from' // trim(tile_file)


    open (10, file=trim(tile_file), form='formatted', action='read')

    read (10,*) N_tile
    read (10,*) N_grid          ! some number (?)
    read (10,*) gridname         ! some string describing tile definition grid (?)
    read (10,*) n_lon
    read (10,*) n_lat
    if(N_grid == 2) then
       read (10,*)          ! some string describing ocean grid                   (?)
       read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
       read (10,*)          ! # ocean grid cells in latitude direction (N_j_ocn)  (?)
    endif

    ease_grid = .false.
    col_order = 0
    
    call LDAS_create_grid_g(gridname,n_lon,n_lat,tile_grid_g,i_indg_offset,j_indg_offset,ease_cell_area)

    if (index(tile_grid_g%gridtype,'EASE')/=0)  ease_grid = .true.  ! 'EASE' and 'EASEv2'
    if (index(tile_grid_g%gridtype,'SiB2')/=0)  col_order=1  ! Weiyuan note: grid name for SiB2??
    allocate(tile_coord(N_tile))
    allocate(f2g_tmp(N_tile))
    i= 0
    fid = 0

    ! WJ notes: i and k are the same---global ids
    !           fid --- num in simulation domain
    do k=1,N_tile
       
        read(10,'(A)') tmpline 
        read(tmpline,*) typ
        if(typ == 100 .or. typ ==1100) then ! if it is land
           i=i+1
           tile_coord(i)%tile_id   = k

           if (typ == 100) then
              fid=fid+1
              f2g_tmp(fid) = k
           endif

           if(ease_grid .or. N_grid ==1) then

             ! EASE grid til file has fewer columns 
             ! (excludes "tile_id", "frac_pfaf", and "area")

             read (tmpline,*)                          &
                  tile_coord(i)%typ,              &   !  1
                  tile_coord(i)%pfaf,             &   !  2
                  tile_coord(i)%com_lon,          &   !  3
                  tile_coord(i)%com_lat,          &   !  4
                  tile_coord(i)%i_indg,           &   !  5
                  tile_coord(i)%j_indg,           &   !  6
                  tile_coord(i)%frac_cell             !  7

             tile_coord(i)%frac_pfaf = nodata_generic
             tile_coord(i)%area      = ease_cell_area*tile_coord(i)%frac_cell

           else ! not ease grid

             if (col_order==1) then
                ! old "SiB2_V2" file format
                read (tmpline,*)                       &
                     tile_coord(i)%typ,           &   !  1  
                     tile_coord(i)%pfaf,          &   !  2  *
                     tile_coord(i)%com_lon,       &   !  3
                     tile_coord(i)%com_lat,       &   !  4
                     tile_coord(i)%i_indg,        &   !  5
                     tile_coord(i)%j_indg,        &   !  6
                     tile_coord(i)%frac_cell,     &   !  7
                     tmpint1,                     &   !  8
                     tmpint2,                     &   !  9  *
                     tmpint3,                     &   ! 10
                     tile_coord(i)%frac_pfaf,     &   ! 11
                     tmpint4,                     &   ! 12  (previously "tile_id")
                     tile_coord(i)%area               ! 13

             else

                read (tmpline,*)                       &
                     tile_coord(i)%typ,           &   !  1
                     tile_coord(i)%area,          &   !  2  *
                     tile_coord(i)%com_lon,       &   !  3
                     tile_coord(i)%com_lat,       &   !  4
                     tile_coord(i)%i_indg,        &   !  5
                     tile_coord(i)%j_indg,        &   !  6
                     tile_coord(i)%frac_cell,     &   !  7
                     tmpint1,                     &   !  8
                     tile_coord(i)%pfaf,          &   !  9  *
                     tmpint2,                     &   ! 10
                     tile_coord(i)%frac_pfaf,     &   ! 11
                     tmpint3                          ! 12  * (previously "tile_id")

                ! change units of area to [km^2]  - 23 Sep 2010: fixed units, reichle

                tile_coord(i)%area = tile_coord(i)%area*MAPL_RADIUS*MAPL_RADIUS/1000./1000.

             end if ! col_order 1

           end if  ! (ease_grid)

          ! fix i_indg and j_indg such that they refer to a global grid
          ! (see above)

           tile_coord(i)%i_indg = tile_coord(i)%i_indg + i_indg_offset
           tile_coord(i)%j_indg = tile_coord(i)%j_indg + j_indg_offset
        else
          exit ! land comes first in the til file
        endif
    end do
    close(10)

    N_tile_land=i
    allocate(tile_coord_land(N_tile_land))
    tile_coord_land=tile_coord(1:N_tile_land)

    if(present(f2g)) then
       allocate(f2g(fid))
       f2g = f2g_tmp(1:fid)
    endif

    call read_catchment_def( catch_file, N_tile_land, tile_coord_land )

    ! ----------------------------------------------------------------------
    !
    ! if still needed read gridded elevation (check only first tile!)

    if ( abs(tile_coord_land(1)%elev-nodata_generic)<nodata_tol_generic ) then

       i=index(catch_file,'/clsm/')
       fname = catch_file(1:i)//'topo_DYN_ave_*.data'
       call system('ls '//trim(fname) // ' >topo_DYN_ave.file')
       open(10,file='topo_DYN_ave.file', action='read')
       fname= ''
       read(10,'(A)') fname
       !close(10,status='DELETE')
       close(10)
       call read_grid_elev( trim(fname), tile_grid_g, N_tile_land, tile_coord_land )
    end if

    
    if ( abs(tile_coord_land(1)%elev-nodata_generic)<nodata_tol_generic ) then

       if (logit) write (logunit,*) 'WARNING: tile elevation NOT avaialable'

    end if
    ! ----------------------------------------------------------------------
    !
    ! fix dateline bug that existed up to and including MERRA version of
    !  *.til and catchment.def files

    call fix_dateline_bug_in_tilecoord( N_tile_land, tile_grid_g, tile_coord_land ) 
    deallocate(tile_coord)
    deallocate(f2g_tmp)

  end subroutine LDAS_read_land_tile

  subroutine read_til_file( tile_coord_path, tile_coord_file, tile_grid, tile_coord )
    
    ! read tile coordinates from GEOS5 *.til file 
    ! optionally assemble parameters of tile definition grid
    !
    ! IMPORTANT: overlaying grid must be *global*
    !            (with matching i_indg and j_indg entries in til file)
    !
    ! inputs:
    !
    !  tile_coord_file : full path to *.til file
    !  tile_coord_file : name of *.til file
    !
    ! outputs:
    !
    !  tile_grid   : parameters of tile definition grid
    !
    !  tile_coord : coordinates of tiles (see tile_coord_type),
    !               implemented as pointer which is allocated in 
    !               this subroutine
    !               NOTE: number of tiles can be diagnosed 
    !                     with size(tile_coord)
    !
    !
    ! "tile_id" is no longer read from *.til file and is now set in this 
    ! subroutine to match order of tiles in *.til file
    ! - reichle, 22 Aug 2013
    !
    ! -------------------------------------------------------------
    
    implicit none
    
    character(*), intent(in) :: tile_coord_path
    character(*), intent(in) :: tile_coord_file
    
    type(grid_def_type), intent(out), optional :: tile_grid
    
    type(tile_coord_type), dimension(:), pointer, optional :: tile_coord ! out
    
    ! locals

    real    :: ease_cell_area
    
    integer :: i, N_tile, tmpint1, tmpint2, tmpint3, tmpint4
    integer :: i_indg_offset, j_indg_offset, col_order

    logical :: date_line_on_center, pole_on_center, ease_grid

    character(300) :: fname

    character(len=*), parameter :: Iam = 'read_til_file'
    character(len=400) :: err_msg

    ! ---------------------------------------------------------------

    fname = trim(tile_coord_path) // '/' // trim(tile_coord_file)

    i_indg_offset = 0
    j_indg_offset = 0
    
    ! read file header 
    
    if (logit) write (logunit,'(400A)') 'read_til_file(): reading from' // trim(fname)
    if (logit) write (logunit,*)
    
    open (10, file=trim(fname), form='formatted', action='read') 
    
    read (10,*) N_tile
    read (10,*)          ! some number (?)
    read (10,*)          ! some string describing tile definition grid (?)
    read (10,*) tmpint1
    read (10,*) tmpint2
    read (10,*)          ! some string describing ocean grid                   (?)
    read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
    read (10,*)          ! # ocean grid cells in latitude direction (N_j_ocn)  (?)
    
    if (logit) write (logunit,*) 'file contains coordinates for ', N_tile, ' tiles' 
    if (logit) write (logunit,*)
        
    ! ---------------------------------------------------------------
    
    ! find out order of columns in *.til file
    
    col_order = 0      ! default: "Fortuna-2" and newer CVS tags
    
    if (index(tile_coord_path, 'SiB2_V2')/=0)  col_order=1  ! prior to "Fortuna-2"
    
    ! find out grid type (default: assume LatLon grid)

    ease_grid = .false.
    
    if (index(tile_coord_file, 'EASE')/=0)  ease_grid = .true.  ! 'EASE' and 'EASEv2'
    
    ! ---------------------------------------------------------------
    !
    ! determine tile_grid 
    
    if (present(tile_grid)) then       
       
       ! initialize all fields to no-data values

       call init_grid_def_type(tile_grid)
       
       tile_grid%N_lon = tmpint1
       tile_grid%N_lat = tmpint2
       
       tile_grid%i_offg = 0  ! tile_grid refers to *global* grid
       tile_grid%j_offg = 0  ! tile_grid refers to *global* grid

       ! ----------------
       !
       ! find out whether date line is on edge or through center of grid cell

       select case (trim(tile_coord_file))
          
       case('FV_144x91_DC_360x180_DE.til',   &  ! 2   deg by 2.5  deg
            'FV_144x91_DC_576x540_DE.til',   &  ! 2   deg by 2.5  deg
            'FV_288x181_DC_360x180_DE.til',  &  ! 1   deg by 1.25 deg
            'FV_540x361_DC_360x180_DE.til',  &  ! 1/2 deg by 2/3  deg
            'FV_576x361_DC_360x180_DE.til',  &  ! 1/2 deg by 5/8  deg
            'FV_1152x721_DC_360x180_DE.til'  &  ! 1/4 deg by 5/16 deg
            )

          ! Fortuna and before

          date_line_on_center = .true.
          pole_on_center      = .true.

       case('DC0144xPC0091_DE0360xPE0180-Pfafstetter.til', & ! 2   deg by 2.5  deg 
            'DC0144xPC0091_DE1440xPE0720-Pfafstetter.til', & ! 2   deg by 2.5  deg 
            'DC0144xPC0091_DE2880xPE1440-Pfafstetter.til', & ! 2   deg by 2.5  deg 
            'DC0288xPC0181_DE0360xPE0180-Pfafstetter.til', & ! 1   deg by 1.25 deg
            'DC0288xPC0181_DE1440xPE0720-Pfafstetter.til', & ! 1   deg by 1.25 deg
            'DC0288xPC0181_DE2880xPE1440-Pfafstetter.til', & ! 1   deg by 1.25 deg
            'DC0540xPC0361_DE0360xPE0180-Pfafstetter.til', & ! 1/2 deg by 2/3  deg
            'DC0576xPC0361_DE0360xPE0180-Pfafstetter.til', & ! 1/2 deg by 5/8  deg
            'DC0576xPC0361_DE1440xPE0720-Pfafstetter.til', & ! 1/2 deg by 5/8  deg
            'DC0576xPC0361_DE2880xPE1440-Pfafstetter.til', & ! 1/2 deg by 5/8  deg
            'DC1152xPC0721_DE0360xPE0180-Pfafstetter.til', & ! 1/4 deg by 5/16 deg
            'DC1152xPC0721_DE1440xPE0720-Pfafstetter.til', & ! 1/4 deg by 5/16 deg
            'DC1152xPC0721_DE2880xPE1440-Pfafstetter.til'  & ! 1/4 deg by 5/16 deg
            )

          ! Ganymed

          date_line_on_center = .true.
          pole_on_center      = .true.

       case('tile.data') ! GEOS5 GCM/DAS convention for tile coord file in work directory
          
          date_line_on_center = .true.
          pole_on_center      = .true.
          
       case('FV_144x91_DE_360x180_DE.til')     ! GEOS-5 2 deg by 2.5 deg DE grid
          
          date_line_on_center = .false.
          pole_on_center      = .true.
          
       case('FV_360x180_DE_360x180_DE.til', &  ! GSWP2 1 deg by 1 deg (orig GSWP-2 resolution)
            'PE_360x180_DE_288x270_DE.til'  &  ! GSWP2 1 deg by 1 deg w/ irregular tiles
            )
          
          date_line_on_center = .false.
          pole_on_center      = .false.
          
       case('PE_720x360_DE_115x48_US.til')     ! CONUS 0.5-degree grid
          
          date_line_on_center = .false.
          pole_on_center      = .false.
          
          ! fix i_indg and j_indg (in the *.til file, i_indg and j_indg are NOT 
          ! relative to the global grid, contrary to the convention in the LDASsa
          ! driver)

          ! fix i_indg and j_indg:
          ! in the *.til file, i_indg and j_indg are relative to the a grid with
          !    ll_lon = -125 
          !    ll_lat =   25
          ! contrary to the convention in the LDASsa driver that needs i_indg and
          ! j_indg relative to a global grid
          
          i_indg_offset = 110
          j_indg_offset = 230
          
       case('PE_2880x1440_DE_464x224_NLDAS.til')  ! NLDAS 1/8-degree grid
          
          date_line_on_center = .false.
          pole_on_center      = .false.
          
          ! fix i_indg and j_indg:
          ! in the *.til file, i_indg and j_indg are relative to the a grid with
          !    ll_lon = -125 
          !    ll_lat =   25
          ! contrary to the convention in the LDASsa driver that needs i_indg and
          ! j_indg relative to a global grid
          
          i_indg_offset = 440
          j_indg_offset = 920
          
       case('SMAP_EASEv2_M36_964x406.til',   &  ! SMAP 'M36' cylindrical EASEv2 grid 
            'SMAP_EASEv2_M09_3856x1624.til', &  ! SMAP 'M09' cylindrical EASEv2 grid 
            'SMAP_EASE_M36_963x408.til' ,    &  ! SMAP 'M36' cylindrical EASE   grid
            'SMAP_EASE_M09_3852x1632.til',   &  ! SMAP 'M09' cylindrical EASE   grid
            'EASE_M25_1383x586.til'          &  !      'M25' cylindrical EASE   grid
            )
          
          date_line_on_center = .false.
          pole_on_center      = .false.

       case default
          
          err_msg = 'unknown tile definitions, filename='//trim(tile_coord_file)
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end select
       
       ! ----------------
       
       if (ease_grid) then

          tile_grid%ind_base = 0
          
          if   (                                                       &
               (index(tile_coord_file, 'EASE_M')  /=0) .or.            &
               (index(tile_coord_file, 'EASEv2_M')/=0)       )  then   

             
             ! global cylindrical EASE grid 
             
             tile_grid%ll_lon = -180.
             tile_grid%ur_lon =  180.
             
             tile_grid%i_dir  = +1
             tile_grid%j_dir  = -1
             
             if     (index(tile_coord_file, 'EASEv2_M36')/=0) then  ! version *2*
                
                tile_grid%gridtype = 'EASEv2_M36'                          
                
                tile_grid%ll_lat =  -85.04456
                tile_grid%ur_lat =   85.04456

                ease_cell_area   = 1298.320938704616
                
             elseif (index(tile_coord_file, 'EASEv2_M09')/=0) then  ! version *2*

                tile_grid%gridtype = 'EASEv2_M09'  

                tile_grid%ll_lat =  -85.04456
                tile_grid%ur_lat =   85.04456
                
                ease_cell_area   =   81.145058669038477

             elseif (index(tile_coord_file, 'EASE_M36')/=0) then
                
                tile_grid%gridtype = 'EASE_M36'                          
                
                tile_grid%ll_lat =  -86.62256 ! minimal change, reichle, 5 Apr 2013
                tile_grid%ur_lat =   86.62256 ! minimal change, reichle, 5 Apr 2013

                ease_cell_area   = 1296.029001087600 
                
             elseif (index(tile_coord_file, 'EASE_M09')/=0) then

                tile_grid%gridtype = 'EASE_M09'                          

                tile_grid%ll_lat =  -86.62256 ! minimal change, reichle, 5 Apr 2013
                tile_grid%ur_lat =   86.62256 ! minimal change, reichle, 5 Apr 2013
                
                ease_cell_area   =   81.001812568020028

             elseif (index(tile_coord_file, 'EASE_M25')/=0) then

                tile_grid%gridtype = 'EASE_M25'                          

                tile_grid%ll_lat =  -86.7167 ! need to double-check (reichle, 11 May 2011)
                tile_grid%ur_lat =   86.7167 ! need to double-check (reichle, 11 May 2011)
                
                ease_cell_area   =   628.38080962
                
             else
                
                err_msg = 'unknown EASE grid tile defs, file= ' &
                     // trim(tile_coord_file)
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                
             end if
             
             tile_grid%dlon   = 360./real(tile_grid%N_lon) 
             tile_grid%dlat   = &
                  (tile_grid%ur_lat-tile_grid%ll_lat)/real(tile_grid%N_lat) ! *avg* dlat!
             
          else
             
             err_msg = 'unknown EASE grid tile defs, file= ' &
                  //trim(tile_coord_file)
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

          end if
          
       else ! assume regular LatLon grid

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
             
             tile_grid%ll_lon = -180. - tile_grid%dlon/2.
             tile_grid%ur_lon =  180. - tile_grid%dlon/2.  ! fixed 20 sep 2010, reichle
             
          else
             
             tile_grid%ll_lon = -180. 
             tile_grid%ur_lon =  180. 
          
          end if
          
       end if

       if (logit) then
          
          write (logunit,*) 'tile definition grid:'
          if (date_line_on_center) then
             write (logunit,*) '  date line on center'
          else
             write (logunit,*) '  date line on edge'
          end if
          if (pole_on_center) then
             write (logunit,*) '  pole on center'
          else
             write (logunit,*) '  pole on edge'
          end if
          
          tmpstring40 = 'tile_grid_g'

          call io_grid_def_type( 'w', logunit, tile_grid, tmpstring40 )
          
       end if
       
    end if
    
    ! -----------------------------------------------------------------
    
    if (present(tile_coord)) then        ! read body of file

       allocate(tile_coord(N_tile))
       
       do i=1,N_tile
          
          ! set "tile_id" to match order in which tiles are read from *.til file

          tile_coord(i)%tile_id   = i
          
          if (ease_grid) then
                
             ! EASE grid til file has fewer columns 
             ! (excludes "tile_id", "frac_pfaf", and "area")
             
             read (10,*)                          &   
                  tile_coord(i)%typ,              &   !  1
                  tile_coord(i)%pfaf,             &   !  2
                  tile_coord(i)%com_lon,          &   !  3
                  tile_coord(i)%com_lat,          &   !  4
                  tile_coord(i)%i_indg,           &   !  5
                  tile_coord(i)%j_indg,           &   !  6
                  tile_coord(i)%frac_cell             !  7
                
             tile_coord(i)%frac_pfaf = nodata_generic
             tile_coord(i)%area      = ease_cell_area*tile_coord(i)%frac_cell
                
          else

             if (col_order==1) then
             
                ! old "SiB2_V2" file format
                
                read (10,*)                       &   
                     tile_coord(i)%typ,           &   !  1  
                     tile_coord(i)%pfaf,          &   !  2  *
                     tile_coord(i)%com_lon,       &   !  3
                     tile_coord(i)%com_lat,       &   !  4
                     tile_coord(i)%i_indg,        &   !  5
                     tile_coord(i)%j_indg,        &   !  6
                     tile_coord(i)%frac_cell,     &   !  7
                     tmpint1,                     &   !  8
                     tmpint2,                     &   !  9  *
                     tmpint3,                     &   ! 10
                     tile_coord(i)%frac_pfaf,     &   ! 11
                     tmpint4,                     &   ! 12  (previously "tile_id")
                     tile_coord(i)%area	              ! 13
             
             else
             
                read (10,*)                       &   
                     tile_coord(i)%typ,           &   !  1
                     tile_coord(i)%area,          &   !  2  *
                     tile_coord(i)%com_lon,       &   !  3
                     tile_coord(i)%com_lat,       &   !  4
                     tile_coord(i)%i_indg,        &   !  5
                     tile_coord(i)%j_indg,        &   !  6
                     tile_coord(i)%frac_cell,     &   !  7
                     tmpint1,                     &   !  8
                     tile_coord(i)%pfaf,          &   !  9  *
                     tmpint2,                     &   ! 10
                     tile_coord(i)%frac_pfaf,     &   ! 11
                     tmpint3                          ! 12  * (previously "tile_id")
                
                ! change units of area to [km^2]  - 23 Sep 2010: fixed units, reichle
                
                tile_coord(i)%area = tile_coord(i)%area*MAPL_RADIUS*MAPL_RADIUS/1000./1000.
                
             end if
             
          end if  ! (ease_grid)
                    
          ! fix i_indg and j_indg such that they refer to a global grid
          ! (see above)
          
          tile_coord(i)%i_indg = tile_coord(i)%i_indg + i_indg_offset
          tile_coord(i)%j_indg = tile_coord(i)%j_indg + j_indg_offset

       end do
       
    end if

    ! -----------------------------------------------------------------

    close(10, status='keep')
    
  end subroutine read_til_file

  ! **********************************************************************
  
  subroutine get_tile_grid( N_tile, tile_coord, tile_grid_g, tile_grid )     
    
    ! get matching tile_grid for given tile_coord and (global) tile_grid_g
    !
    ! reichle, 20 June 2012 -- moved from within domain_setup() 
    !                           for re-use in get_obs_pred() 
    !
    ! -------------------------------------------------------------------
    
    integer,               intent(in)            :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! in
    
    type(grid_def_type),   intent(in)            :: tile_grid_g
    type(grid_def_type),   intent(out)           :: tile_grid
    
    ! local variables
    
    integer :: n
    
    real    :: this_minlon, this_minlat, this_maxlon, this_maxlat
    
    real    :: min_min_lon, min_min_lat, max_max_lon, max_max_lat
    
    integer :: ind_i_min, ind_i_max, ind_j_min, ind_j_max
    
    integer :: this_i_indg, this_j_indg
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

    ! for c3 grid, only get the ll_,ur_ lat and lon, the index is meaning less
    ! it will be used in creating the lat_lon pert_grid
    if(index(tile_grid_g%gridtype,"c3") /=0) then
 
      ! do n=1,N_tile
      ! 
      !    this_minlon  = tile_coord(n)%com_lon
      !    this_minlat  = tile_coord(n)%com_lat
      !    this_maxlon  = tile_coord(n)%com_lon
      !    this_maxlat  = tile_coord(n)%com_lat
      !    min_min_lon = min( min_min_lon, this_minlon)
      !    min_min_lat = min( min_min_lat, this_minlat)
      !    max_max_lon = max( max_max_lon, this_maxlon)
      !    max_max_lat = max( max_max_lat, this_maxlat)
      ! enddo
       tile_grid=tile_grid_g
      ! tile_grid%ll_lon= min_min_lon
      ! tile_grid%ur_lon= max_max_lon
      ! tile_grid%ll_lat= min_min_lat
      ! tile_grid%ur_lat= max_max_lat
       return
    endif   

    do n=1,N_tile
       
       this_minlon  = tile_coord(n)%min_lon
       this_minlat  = tile_coord(n)%min_lat
       this_maxlon  = tile_coord(n)%max_lon
       this_maxlat  = tile_coord(n)%max_lat
       
       this_i_indg  = tile_coord(n)%i_indg
       this_j_indg  = tile_coord(n)%j_indg
       
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
       
    elseif ( (index(tile_grid_g%gridtype,'EASE_M')  /=0)     .or.          &
             (index(tile_grid_g%gridtype,'EASEv2_M')/=0)          )  then
       
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
  subroutine get_number_of_tiles_in_cell_ij( N_tile, tile_coord, tile_grid, &
       N_tile_in_cell_ij)
    
    ! find out how many tiles are in a given tile definition grid cell
    ! reichle, 22 Jul 2005
    ! split for the old get_tile_num_in_cell_ij
    ! ----------------------------------------------------------
    
    implicit none
    
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! input
    
    type(grid_def_type), intent(in) :: tile_grid
        
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
       
       i = tile_coord(n)%i_indg - off_i
       j = tile_coord(n)%j_indg - off_j
       
       N_tile_in_cell_ij(i,j) = N_tile_in_cell_ij(i,j) + 1
       
    enddo

    if (logit) write (logunit,*) &
        'Maximum number of tiles in tile def grid cell = ', maxval(N_tile_in_cell_ij)
    if (logit) write (logunit,*)

  end subroutine get_number_of_tiles_in_cell_ij

  subroutine get_tile_num_in_cell_ij( N_tile, tile_coord, tile_grid, &
       max_N_tile_in_cells, tile_num_in_cell_ij )
    
    ! find out tile_num in given cells
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domaim of interest.  
    ! With these offsets tile2grid() can be used to map from a
    ! subgrid of "tile_grid_g" to tile space
    !
    ! reichle, 22 Jul 2005
    !
    ! ----------------------------------------------------------
    
    implicit none
    
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! input
    
    type(grid_def_type), intent(in) :: tile_grid
        
    integer, intent(in) ::  max_N_tile_in_cells
    
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
       
       i = tile_coord(n)%i_indg - off_i
       j = tile_coord(n)%j_indg - off_j
       
       N_tile_in_cell_ij(i,j) = N_tile_in_cell_ij(i,j) + 1
       
       k = N_tile_in_cell_ij(i,j)
          
       tile_num_in_cell_ij(i,j,k) = n 
          
    end do
    
  end subroutine get_tile_num_in_cell_ij

  ! *******************************************************************
  
  subroutine extract_land_tiles( N_tile_global, tile_coord_global, &
       N_tile_land, tile_coord_land )
    
    ! extract land tiles from tile_coord_global
    !
    ! When called without optional arguments only counts number of tiles.
    ! When called with optional arguments allocates and fills pointer.
    !
    ! reichle, 28 Jan 2005
    ! reichle, 22 Jul 2005
    !
    ! ----------------------------------------------------------
    
    implicit none
    
    integer, intent(in) :: N_tile_global
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord_global ! input
    
    integer, intent(inout) :: N_tile_land
    
    ! the pointer is an output arguments that is allocated here
    
    type(tile_coord_type), dimension(:), pointer, optional :: tile_coord_land
    
    ! locals 
    
    integer :: n
    
    ! -----------------------------------------------------------------
    !
    ! allocate and initialize pointers if present
    
    if (present(tile_coord_land)) allocate(tile_coord_land(N_tile_land))
    
    ! (re-)initialize
    
    N_tile_land           = 0
    
    do n=1,N_tile_global
       
       ! count number of land tiles
       
       if (tile_coord_global(n)%typ == tile_typ_land) then
          
          N_tile_land = N_tile_land + 1
          
          if (present(tile_coord_land)) &          
               tile_coord_land(N_tile_land) = tile_coord_global(n)
          
       end if
       
    end do
    
    if (logit) write (logunit,*) 'Number of land tiles = ', N_tile_land 
    if (logit) write (logunit,*)
    
  end subroutine extract_land_tiles
  
  ! *******************************************************************

  subroutine get_land_tile_info( tile_coord_path, tile_coord_file,     &
       catchment_def_path, catchment_def_file, res_ftag,               &
       tile_grid, N_tile, tile_coord )
    
    ! get land tile coordinates, parameters of tile definition grid, and 
    ! grid-to-tile mapping from GEOS5 *.til file
    !
    ! note use of optional inout arguments
    !
    ! reichle, 28 Jan 2005
    !
    ! ----------------------------------------------------------
    
    implicit none
    
    character(*), intent(in) :: tile_coord_path, catchment_def_path
    character(*), intent(in) :: tile_coord_file
    character(*), intent(in) :: catchment_def_file, res_ftag

    type(grid_def_type), intent(out) :: tile_grid
    
    integer, intent(out) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! out
    
    ! locals

    character(300) :: fname

    integer :: N_tile_tmp
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord_tmp

    logical :: file_exists

    ! --------------------------------------------------------------
    
    nullify(tile_coord_tmp)
    
    ! get parameters of tile definition grid and global tiles
        
    call read_til_file( tile_coord_path, tile_coord_file,  &
         tile_grid=tile_grid, tile_coord=tile_coord_tmp )
    
    N_tile_tmp = size(tile_coord_tmp)
    
    ! first call counts land tiles in tile_coord_file
    
    call extract_land_tiles( N_tile_tmp, tile_coord_tmp, N_tile )
    
    ! second call allocates and fills tile_coord with land tiles
    
    call extract_land_tiles( N_tile_tmp, tile_coord_tmp, N_tile, tile_coord )
    
    deallocate(tile_coord_tmp)
    
    ! ----------------------------------------------------------------------
    !
    ! add lat-lon bounding box and (possibly elevation) for each tile 
    ! from "catchment.def" file
    ! 

    fname = trim(catchment_def_path) // '/' // trim(catchment_def_file)
    
    inquire(file=fname, exist=file_exists)
    
    if (.not. file_exists) then
    
       fname = trim(catchment_def_path) // '/clsm/' // trim(catchment_def_file)
       
       inquire(file=fname, exist=file_exists)
       
       if (.not. file_exists) then
          
          if (logit) write(logunit,*) 'cannot find catchment def file'
          if (logit) write(logunit,*) 'catchment_def_path = ', trim(catchment_def_path)
          if (logit) write(logunit,*) 'catchment_def_file = ', trim(catchment_def_file)
          
       end if
       
    end if
    
    call read_catchment_def( fname, N_tile, tile_coord )
    
    ! ----------------------------------------------------------------------
    !
    ! if still needed read gridded elevation (check only first tile!)
    
    if ( abs(tile_coord(1)%elev-nodata_generic)<nodata_tol_generic ) then
       
       fname = trim(tile_coord_path) 
       
       fname = trim(fname) // '/' // 'topo_DYN_ave_' // trim(res_ftag) // '.data'
       
       call read_grid_elev( fname, tile_grid, N_tile, tile_coord )    
       
    end if
    
    ! check whether elevation data is now available

    if ( abs(tile_coord(1)%elev-nodata_generic)<nodata_tol_generic ) then
    
       if (logit) write (logunit,*) 'WARNING: tile elevation NOT avaialable'
       
    end if
    
    ! ----------------------------------------------------------------------
    !
    ! fix dateline bug that existed up to and including MERRA version of
    !  *.til and catchment.def files
    
    call fix_dateline_bug_in_tilecoord( N_tile, tile_grid, tile_coord )
   
  end subroutine get_land_tile_info

  ! *******************************************************************
  
  subroutine read_grid_elev( fname, tile_grid, N_tile, tile_coord )

    ! read gridded elevation file (for GEOS-5 discretizations; NOT available
    ! for EASE grids, where elevation information is in catchment.def file)
    
    ! reichle,  8 Dec 2011: bug fix -- bin elev data is stored in single record

    implicit none
    
    character(*),      intent(in) :: fname
    
    type(grid_def_type), intent(in) :: tile_grid
        
    integer,             intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! inout
    
    ! local variables
    
    integer :: istat, i, j

    real, dimension(tile_grid%N_lon,tile_grid%N_lat) :: grid_elev
    
    ! ---------------------------------------------------------------
    
    if (logit) write (logunit,'(400A)') 'read_grid_elev(): reading from' // trim(fname)
    if (logit) write (logunit,*)
    
    ! open, read, and close file
    
    open(10, file=fname, form='unformatted', status='old', &
         convert='little_endian', action='read', iostat=istat)
 
    if (istat/=0) then
       
       if (logit) write (logunit,*) 'WARNING: cannot open file, returning'

       grid_elev = nodata_generic

    else
       
       !  binary elevation data is stored in single Fortran record
       
       read (10) (( grid_elev(i,j), i=1,tile_grid%N_lon), j=1,tile_grid%N_lat)
       
       close (10,status='keep')

       if (logit) write (logunit,*) 'done reading file'
       if (logit) write (logunit,*)    
       
    end if

    ! ---------------------------
    
    ! map elevation to tiles
    
    do i=1,N_tile
       
       tile_coord(i)%elev = grid_elev( tile_coord(i)%i_indg, tile_coord(i)%j_indg )
       
    end do
    
  end subroutine read_grid_elev
  
  ! *******************************************************************

  subroutine fix_dateline_bug_in_tilecoord( N_tile, tile_grid, tile_coord )

    ! bug in com_lon and minlon/maxlon for tiles straddling the dateline
    ! existed through (and including) MERRA tag
    !
    ! for now do not to allow any tiles that straddle the dateline
    !
    ! reichle,  5 Feb 2008

    implicit none
    
    integer, intent(in) :: N_tile
    
    type(grid_def_type), intent(in) :: tile_grid
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! inout
    
    ! locals

    real, parameter :: latlon_tol = 1e-4
    
    integer :: k
    
    real :: this_minlon, this_maxlon
    
    ! --------------------------------------------
    !
    ! check whether there could be tiles crossing the dateline
    
    if ( (tile_grid%ll_lon<-180.) .and.  &
         (.not. (tile_grid%ll_lon-nodata_generic)<nodata_tol_generic) ) then
       
       do k=1,N_tile
          
          ! min/max longitude of tile definition grid cell
          
          this_minlon = &
               tile_grid%ll_lon + real(tile_coord(k)%i_indg-1)*tile_grid%dlon
          
          this_maxlon = this_minlon + tile_grid%dlon
          
          ! fix min/max lon
          
          if (abs(tile_coord(k)%max_lon-tile_coord(k)%min_lon)>tile_grid%dlon+latlon_tol) then
             
             if (logit) write (logunit,*) &
                  'resetting min/maxlon in tile_id=', tile_coord(k)%tile_id
             
             ! "push" tile to the east of the dateline
             
             tile_coord(k)%min_lon = -180.
             tile_coord(k)%max_lon = -180. + 0.5*tile_grid%dlon
             
          end if
          
          
          ! check that com_lon is within tile definition grid cell
          
          if ( .not. (  &
               (tile_coord(k)%com_lon >= this_minlon) .and.                 &
               (tile_coord(k)%com_lon <= this_maxlon)        ) ) then
             if (logit) write (logunit,*) &
                  'resetting com_lon in tile_id=', tile_coord(k)%tile_id
             
             tile_coord(k)%com_lon = &
                  0.5*(tile_coord(k)%min_lon + tile_coord(k)%max_lon)
             
          end if
          
       end do
    end if
    
  end subroutine fix_dateline_bug_in_tilecoord

  ! **************************************************************

  subroutine reorder_tiles( reorder, pfaf_system, N_tile, tile_coord, d2g, N_tiles_cont )
    
    ! Re-order tile_coord (and d2g) for more optimal domain decomposition and 
    ! to address the date-line issue (subroutine get_obs_pred() requires that 
    ! local domains must not have tiles on both sides of the dateline).
    !
    ! If input argument "reorder" is ".false." and "pfaf_system=0", assume that tiles have 
    ! already been reordered.  Check for obvious violations and only return "N_tiles_cont".
    ! 
    ! Typically done only by the master process (because the re-ordering requires
    ! a second copy of the full domain tile coord structure).
    !
    ! reichle, 26 June 2012
    ! reichle, 14 June 2013 - L1bas bug fix
    ! reichle, 31 Oct  2017 - added functionality for SRTM-based tiles
    !
    ! ---------------------------------------------------------------
        
    implicit none
    
    logical,                                      intent(in)    :: reorder
    
    integer,                                      intent(in)    :: pfaf_system, N_tile
    
    type(tile_coord_type), dimension(N_tile),     intent(inout) :: tile_coord  
    
    integer,               dimension(N_tile),     intent(inout) :: d2g

    integer,               dimension(N_cont_max), intent(out)   :: N_tiles_cont
    
    ! local variables
    
    ! (max) number of Level 1 Pfafstetter basins per continent    
    
    integer, parameter :: N_L1bas_per_cont = 10   
    
    integer                                          :: i, j, j_prev, ii, jj, kk

    integer, dimension(N_cont_max)                   :: cont_beg, cont_end, jj_cont
    
    integer, dimension(N_cont_max*N_L1bas_per_cont)  :: N_tiles_L1bas, j_L1bas
    integer, dimension(N_cont_max*N_L1bas_per_cont)  :: jstart_L1bas, jend_L1bas
    
    integer, dimension(N_tile)                       :: cont, L1bas, d2g_tmp
    
    type(tile_coord_type), dimension(:), pointer     :: tc_tmp => null()    
    
    character(len=*), parameter :: Iam = 'reorder_tiles'
    character(len=400) :: err_msg
    character(len=9)   :: tmpstr9
    character(len=3)   :: tmpstr3

    ! -----------------------------------------------------------------------------

    select case (pfaf_system)
       
    case (0)    
       
       ! Level 6 Pfafstetter numbers based on Hydro-1k
       !
       ! - for boundary condition files through ~2016
       ! - tile_coord%pfaf is actual Level-6 Pfafstetter number

    ! -----------------------------------------------------------------------------
    !
    ! Step 1: Count tiles on each continent and re-assign "problem" tiles 
    !         to continents that are on the same side of the dateline.
    !         Also process Level 1 basin information for later use
    
    N_tiles_cont  = 0
    N_tiles_L1bas = 0
    
    j_prev = 0
    
    do i=1,N_tile
       
       ! definition of continents
       
       !  NAMER  = 1
       !  SAMER  = 2
       !  EUROPE = 3
       !  AFRICA = 4
       !  ASIA   = 5
       !  AUST   = 6
       
       cont( i) =  (tile_coord(i)%pfaf/1000000)      + 1;
       
       ! definition of Level 1 basins
       !
       ! bug fix - reichle, 14 Jun 2013 

       L1bas(i) = ((tile_coord(i)%pfaf/ 100000) - (cont(i)-1)*10) + 1;

       ! -------------
       
       ! reassign some islands to make bounding boxes more reasonable
       
       ! if SAMER and west of 20W reassign to AFRICA (S Atlantic)
       
       if (cont(i)==2 .and. tile_coord(i)%com_lon>-20.) then
          
          cont( i) = 4
          L1bas(i) = 6 
          
       end if
       
       ! if AUST and west of dateline reassign to NAMER (Hawaii)
       
       if ( cont(i)==6                 .and. &
            tile_coord(i)%com_lon< 0.  .and. &
            tile_coord(i)%com_lat>15.)       &
            cont(i) = 1
       
       ! NOTE: could introduce additional "continents", eg., for
       !       the Pacific islands
       
       ! -------------
       
       ! reassignments to address dateline issue
       
       ! if NAMER and west of dateline reassign to ASIA
       
       if (cont(i)==1 .and. tile_coord(i)%com_lon>0.)  cont(i) = 5
       
       ! if SAMER and west of dateline reassign to AUST
       
       if (cont(i)==2 .and. tile_coord(i)%com_lon>0.)  cont(i) = 6
       
       ! if ASIA and east of dateline reassign to NAMER
       
       if (cont(i)==5 .and. tile_coord(i)%com_lon<0.) then
          
          cont(i)  = 1
          L1bas(i) = 2
          
       end if
       
       ! if AUST and east of dateline reassign to SAMER
       
       if (cont(i)==6 .and. tile_coord(i)%com_lon<0.) then
          
          cont(i)  = 2
          L1bas(i) = 2
          
       end if
       
       ! count number of tiles on each continent
       
       N_tiles_cont(cont(i)) = N_tiles_cont(cont(i)) + 1
       
       ! count number of tiles in each basin
       
       j = (cont(i)-1)*N_L1bas_per_cont + L1bas(i)
       
       N_tiles_L1bas(j) = N_tiles_L1bas(j) + 1

       ! check for obvious violations of the assumption that tiles have already
       ! been reordered
       
       if (.not. reorder) then
          
          if (j<j_prev) then
             err_msg = 'tiles cannot have been reordered'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          j_prev = j
          
       end if
       
    end do
    
    ! ---------------------------------------------------------------
    !
    ! Step 2: Re-arrange tile_coord in order of revised continent and
    !         Level 1 basin IDs
    
    if (reorder) then
       
       jstart_L1bas(1)                         = 1
       jend_L1bas(N_cont_max*N_L1bas_per_cont) = N_tile
       
       do j=1,(N_cont_max*N_L1bas_per_cont-1)
          
          jend_L1bas(j)     = jstart_L1bas(j) + N_tiles_L1bas(j) - 1
          jstart_L1bas(j+1) = jend_L1bas(j) + 1
          
       end do
       
       j_L1bas = jstart_L1bas
       
       allocate(tc_tmp(N_tile))
       
       do i=1,N_tile
          
          j = (cont(i)-1)*N_L1bas_per_cont + L1bas(i)
          
          tc_tmp( j_L1bas(j)) = tile_coord(i)
          
          d2g_tmp(j_L1bas(j)) = d2g(       i)
          
          j_L1bas(j) = j_L1bas(j) + 1
          
       end do
       
       tile_coord = tc_tmp
       
       d2g        = d2g_tmp       
       
       deallocate(tc_tmp)
       
    end if

       ! ------------------------------------------------------------

       
    case (1)   

       ! Level 12 Pfafstetter catchments based on SRTM
       !
       ! - for boundary condition files from ~2017
       ! - tile_coord%pfaf is a consecutive number (not the actual Level-12 Pfafstetter number)
       
       ! -----------------------------------------------------------------------------
       !
       ! Step 1: Count tiles on each continent and re-assign "problem" tiles 
       !         to continents that are on the same side of the dateline.
       
       N_tiles_cont  = 0
       
       do i=1,N_tile
          
          ! definition of continents (retain order of continents from pfaf_system=0)
          
          ! Asia   :      1- 75368  -->  continent = 5
          ! Africa :  75369-140751  -->  continent = 4
          ! NA     : 140752-189105  -->  continent = 1
          ! Europe : 189106-229074  -->  continent = 3
          ! SA     : 229075-267083  -->  continent = 2
          ! AU     : 267084-291284  -->  continent = 6
          
          select case (tile_coord(i)%pfaf)
             
          case (      1 :  75368);  cont(i) = 5
          case (  75369 : 140751);  cont(i) = 4
          case ( 140752 : 189105);  cont(i) = 1
          case ( 189106 : 229074);  cont(i) = 3
          case ( 229075 : 267083);  cont(i) = 2
          case ( 267084 : 291284);  cont(i) = 6
             
          case default          
             
             write(tmpstr9,'(i9)') tile_coord(i)%pfaf
             write(tmpstr3,'(i3)') pfaf_system
             
             err_msg = 'unknown tile_coord(i)%pfaf for i=' // tmpstr9 // ', pfaf_system=' // tmpstr3 
             
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             
          end select
          
          
          ! -------------
             
          ! reassignments to address dateline issue
          
          ! if ASIA and east of dateline reassign to NAMER
          
          if (cont(i)==5 .and. tile_coord(i)%com_lon<0.)  cont(i)  = 1
          
          ! if AUST and east of dateline reassign to SAMER
          
          if (cont(i)==6 .and. tile_coord(i)%com_lon<0.)  cont(i)  = 2

          ! -------------
                    
          ! count number of tiles on each continent
          
          N_tiles_cont(cont(i)) = N_tiles_cont(cont(i)) + 1
          
       end do
       
       ! ---------------------------------------------------------------
       !
       ! Step 2: ALWAYS reorder tile_coord such that "continents"
       !                are arranged in contiguous chunks
       
       ! find index range for each continent
       !
       ! NOTE: if N_cont(kk)=0, then cont_beg(kk) and cont_end(kk) are not used 
       !                             (and cont_beg(kk)=cont_end(kk)+1)

       cont_beg(1) = 1;
       cont_end(1) = N_tiles_cont(1);
       
       do kk=2,N_cont_max
          
          cont_end(kk) = cont_end(kk-1) + N_tiles_cont(kk);
          
          cont_beg(kk) = cont_end(kk)   - N_tiles_cont(kk) + 1;
          
       end do
       
       ! re-arrange order by continent (order within continent is unchanged)

       allocate(tc_tmp(N_tile))

       jj_cont = cont_beg
       
       do ii=1,N_tile
          
          jj = jj_cont( cont(ii) )
          
          tc_tmp( jj) = tile_coord(ii)
          
          d2g_tmp(jj) = d2g(       ii)
          
          jj_cont( cont(ii) ) = jj_cont( cont(ii) ) + 1
          
       end do
       
       tile_coord = tc_tmp
       
       d2g        = d2g_tmp       
       
       deallocate(tc_tmp)
              
       ! ------------------------------------------------------------

    case default
                    
       write(tmpstr3,'(i3)') pfaf_system
       
       err_msg = 'unknown pfaf_system=' // tmpstr3 
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
   
  end subroutine reorder_tiles
  
  ! *******************************************************************
  
  subroutine read_catchment_def( catchment_def_file, N_tile, tile_coord )
    
    ! reichle, 17 May 2011: read elevation data if available
    
    ! format of catchment.def file
    !
    ! Header line: N_tile
    !
    ! Columns: tile_id, Pfaf, min_lon, max_lon, min_lat, max_lat, [elev]
    !
    ! Elevation [m] is ONLY available for EASE grid tile definitions

    implicit none
    
    character(*), intent(in) :: catchment_def_file
    
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! inout
    
    ! locals
    
    integer :: i, istat, tmpint1, sweep
    
    integer, dimension(N_tile) :: tmp_tileid, tmp_pfaf
    
    character(len=*), parameter :: Iam = 'read_catchment_def'
    character(len=400) :: err_msg
    
    ! ---------------------------------------------------------------
    
    ! read file header 
    
    if (logit) write (logunit,'(400A)') &
         'read_catchment_def(): reading from' // trim(catchment_def_file)
    if (logit) write (logunit,*)
    
    ! sweep=1: Try reading 7 columns.  If this fails, try again.
    ! sweep=2: Read only 6 columns.
    
    do sweep=1,2
       
       if (logit) write (logunit,*) 'starting sweep ', sweep

       open (10, file=trim(catchment_def_file), form='formatted', action='read') 
       
       read (10,*) tmpint1
       
       if (logit) write (logunit,*) 'file contains coordinates for ', tmpint1, ' tiles' 
       if (logit) write (logunit,*)
       
       if (N_tile/=tmpint1) then
          print*,"need :", N_tile,"but have: ",tmpint1
          err_msg = 'tile_coord_file and catchment_def_file mismatch. (1)'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       do i=1,N_tile
          
          if (sweep==1) then
             
             ! read 7 columns, avoid using exact format specification
             
             read (10,*, iostat=istat) tmp_tileid(i), tmp_pfaf(i), &
                  tile_coord(i)%min_lon,    &
                  tile_coord(i)%max_lon,    &
                  tile_coord(i)%min_lat,    &
                  tile_coord(i)%max_lat,    &
                  tile_coord(i)%elev
             
          else

             ! read 6 columns, avoid using exact format specification
             
             read (10,*, iostat=istat) tmp_tileid(i), tmp_pfaf(i), &
                  tile_coord(i)%min_lon,    &
                  tile_coord(i)%max_lon,    &
                  tile_coord(i)%min_lat,    &
                  tile_coord(i)%max_lat
             
             tile_coord(i)%elev = nodata_generic
             
          end if
          
          if (istat/=0) then   ! read error
             
             if (sweep==1) then
                
                close(10,status='keep')
                
                if (logit) write (logunit,*) 'sweep 1 failed, trying sweep 2'
                
                exit  ! exit sweep 1, try sweep 2
                
             else
                
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'sweep 2 failed')
                
             end if
             
          end if
          
          if (i==N_tile) then  ! reached end of tile loop w/o read error
             
             close(10,status='keep')
             
             if (logit) write (logunit,*) 'sweep ', sweep, 'successfully completed'
             
             return
             
          end if
          
       end do   ! loop through tiles
       
    end do      ! loop through sweeps
    
    if ( any(tile_coord(1:N_tile)%tile_id/=tmp_tileid) .or.       &
         any(tile_coord(1:N_tile)%pfaf   /=tmp_pfaf)            ) then

       err_msg = 'tile_coord_file and catchment_def_file mismatch. (2)'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if
    
    ! -----------------------------------------------------------------
    
  end subroutine read_catchment_def
  
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
          
          if   (                                                         &
               (index(tile_grid%gridtype, 'EASE_M')  /=0) .or.           &
               (index(tile_grid%gridtype, 'EASEv2_M')/=0)       )  then
             
             ! ASSUMPTION: tiles match EASE or EASEv2 grid cells exactly
             !             (unless "outside" the domain, eg. water surface)
             
             if     (N_tile_in_cell_ij(i_ind,j_ind)==1) then
                
                tile_num(n)=tile_num_in_cell_ij(i_ind,j_ind,1)

                min_dist_x = abs(lon(n) - tile_coord(tile_num(n))%com_lon)
                min_dist_y = abs(lat(n) - tile_coord(tile_num(n))%com_lat) 
                
             elseif (N_tile_in_cell_ij(i_ind,j_ind)==0) then
                
                ! Do nothing.  If given EASE or EASEv2 grid cell is not land, 
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

  logical function is_cat_in_box(                                &
       this_minlon, this_minlat, this_maxlon, this_maxlat,       &
       minlon, minlat, maxlon, maxlat        )
    
    ! determine whether catchment is within bounding box - reichle, 7 May 2003
    
    implicit none
    
    real :: this_minlon, this_minlat, this_maxlon, this_maxlat
    real :: minlon, minlat, maxlon, maxlat
    
    if ( (this_minlon >= minlon) .and.        &
         (this_maxlon <= maxlon) .and.        &
         (this_minlat >= minlat) .and.        & 
         (this_maxlat <= maxlat)       )    then
       is_cat_in_box = .true. 
    else
       is_cat_in_box = .false.
    end if
    
  end function is_cat_in_box
  
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
    
    if     (index(tile_grid%gridtype, 'EASE')/=0) then

       ! EASE grid lat/lon to index provides *global*, *0-based* index!
       
       if     (index(tile_grid%gridtype, 'EASE_M')  /=0) then
          
          call easeV1_convert( tile_grid%gridtype(6:8),  lat, lon, r, s)
          
       elseif (index(tile_grid%gridtype, 'EASEv2_M')/=0) then
          
          call easeV2_convert( tile_grid%gridtype(8:10), lat, lon, r, s)
          
       else
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown EASE grid type')
          
       end if
       
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
  
  
  ! *******************************************************************
  
  subroutine tile2grid( N_tile, tile_coord, tile_grid, tile_data,    &
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
    
  end subroutine tile2grid
  
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

  subroutine grid2tile_new( tile_grid, N_tile, i_indgs,j_indgs, grid_data, tile_data)
    
    ! map from tile definition grid to tile space
    !
    ! NOTE: tile_coord must match tile_data
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domaim of interest.  
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

    character(len=*), parameter :: Iam = 'grid2tile_new'
    character(len=400) :: err_msg
    
    ! ------------------------------------
    
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
        
    do n=1,N_tile
       
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j

       tile_data(n) = grid_data(i,j)
       
    end do
    
  end subroutine grid2tile_new
  
  subroutine grid2tile_old( tile_grid, N_tile, tile_coord, grid_data, tile_data)
    
    ! map from tile definition grid to tile space
    !
    ! NOTE: tile_coord must match tile_data
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domaim of interest.  
    ! With these offsets grid2tile() can be used to map from a
    ! subgrid of "tile_grid_g" to tile space
    !
    ! reichle, 28 Jan 2005
    ! reichle, 16 May 2005 - added offset for "domain" grid
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in) :: grid_data
    
    real, dimension(N_tile), intent(out) :: tile_data
        
    ! local variables
    
    integer :: n, i, j, off_i, off_j

    character(len=*), parameter :: Iam = 'grid2tile_Old'
    character(len=400) :: err_msg
    
    ! ------------------------------------
    
    if (size(tile_coord)/=N_tile) then
       err_msg = 'tile_coord and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
       
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
        
    do n=1,N_tile
       
       i = tile_coord(n)%i_indg - off_i
       j = tile_coord(n)%j_indg - off_j
       
       tile_data(n) = grid_data(i,j)
       
    end do
    
  end subroutine grid2tile_old
  
  ! *******************************************************************
  
  subroutine grid2tile_real8( tile_grid, N_tile, tile_coord, grid_data_8, tile_data_8)
    
    ! map from tile definition grid to tile space for real*8
    !
    ! reichle,  3 Feb 2014
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    real*8, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in) :: grid_data_8
    
    real*8, dimension(N_tile), intent(out) :: tile_data_8
    
    ! local variables
    
    integer :: n, i, j, off_i, off_j
    
    character(len=*), parameter :: Iam = 'grid2tile_real8'
    character(len=400) :: err_msg

    ! ------------------------------------
    
    if (size(tile_coord)/=N_tile) then
       err_msg = 'tile_coord and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
        
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
        
    do n=1,N_tile
       
       i = tile_coord(n)%i_indg - off_i
       j = tile_coord(n)%j_indg - off_j
       
       tile_data_8(n) = grid_data_8(i,j)
       
    end do
    
  end subroutine grid2tile_real8
  
  ! *************************************************************************

  subroutine get_is_land(tile_coord, tile_grid, is_land, rc)

    ! return the logical array is_land with
    !    is_land(i,j) = .true.     if the cell (i,j) has land
    !    is_land(i,j) = .false.    otherwise

    implicit none

    ! input
    type(tile_coord_type), dimension(:), pointer, intent(in) :: tile_coord
    type(grid_def_type), intent(in) :: tile_grid

    ! output
    logical, dimension(tile_grid%N_lon, tile_grid%N_lat), intent(out) :: is_land
    integer, intent(out), optional :: rc

    ! local
    real :: no_land_tol, dummy, w
    integer :: N_tile, n, i, j, off_i, off_j

    ! status report
    if (present(rc)) rc = 1

    no_land_tol = 10.0*epsilon(dummy)

    ! adjust for 0-based indexing (e.g EASE grids)
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    
    N_tile = size(tile_coord)
    is_land = .false.
    do n=1,N_tile
       i = tile_coord(n)%i_indg - off_i
       j = tile_coord(n)%j_indg - off_j

       w = tile_coord(n)%frac_cell
       if (w>no_land_tol) is_land(i,j) = .true.
    end do
       
    ! status report - all is well
    if (present(rc)) rc = 0

  end subroutine get_is_land

  subroutine read_black_or_whitelist(N_cat, fname, blacklist, N_black)

    ! read numbers/IDs of blacklisted catchments 
    !
    ! format of blacklist file: ASCII list of "Pfafstetter+3" numbers
    !
    ! N_black = number of blacklisted catchments
    !
    ! reichle, 2 May 2003
    !
    ! --------------------------------------------------------------

    implicit none

    ! N_cat = max number of catchments allowed in list 
    !         (use N_cat_global when calling this subroutine)

    integer,        intent(in)  :: N_cat
    character(*), intent(in)  :: fname

    integer,        intent(out) :: N_black

    integer, dimension(N_cat), intent(out) :: blacklist

    ! locals

    integer :: istat, tmpint

    logical :: file_exists

    character(len=*), parameter :: Iam = 'read_black_or_whitelist'
    character(len=400) :: err_msg

    ! -----------------------------------------------------------

    N_black = 0

    inquire( file=fname, exist=file_exists)

    if (file_exists) then

       open(10, file=fname, form='formatted', action='read', &
            status='old', iostat=istat)

       if (istat==0) then

          if (logit) write (logunit,*) &
               'reading black- or whitelist from ', trim(fname)
          if (logit) write (logunit,*)

          do
             read(10,*,iostat=istat) tmpint

             if (istat==-1) then
                if (logit) write (logunit,*) ' found ', N_black, ' catchments on list'
                exit
             else if (istat/=0) then
                err_msg = 'read error other than end-of-file'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             else
                N_black = N_black+1
                blacklist(N_black) = tmpint
             end if

             if (N_black>N_cat) then

                write (tmpstring10,*) N_cat
                write (tmpstring40,*) N_black

                err_msg = 'N_black=' // trim(tmpstring40) &
                     // ' > N_cat=' // trim(tmpstring10)
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

             end if
          end do

          close(10,status='keep')

       else

          if (logit) write (logunit,*) &
               'could not open black- or whitelist file ', trim(fname)

       end if

    else

       if (logit) write (logunit,*) &
            'black- or whitelist file does not exist: ', trim(fname)

    end if

    if (logit) write (logunit,*)

  end subroutine read_black_or_whitelist



end module LDAS_TileCoordRoutines


