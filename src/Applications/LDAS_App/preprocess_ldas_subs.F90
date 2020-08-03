
! Subroutines needed for LDAS pre-processing. Moved from LDAS_TileCoordRoutines.F90.
! - wjiang, reichle, 2 Aug 2020

module preprocess_ldas_subs

  use MAPL_BaseMod,                     ONLY:     &
       MAPL_Land
  
  use LDAS_TileCoordType,               ONLY:     &	
       tile_coord_type,                           &	
       grid_def_type
  
  use LDAS_TileCoordRoutines,           ONLY:     &
       LDAS_create_grid_g
  
  use LDAS_ensdrv_Globals,              ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_generic,                            &
       nodata_tol_generic
  
  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_RADIUS                                    ! Earth radius
  
  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR
  
  implicit none

  private
  
  public :: LDAS_read_til_file
  public :: MAPL_Land_ExcludeFromDomain

  ! Tile type for land that is to be excluded from the simulation domain.
  ! (GEOSldas  allows for non-global simulations and repeated "zooming"
  !  of the domain while MAPL generally assumes a complete (global) tile
  !  space.  The *_ExcludeFromDomain tile type makes it possible to work
  !  with complete (global) tile files (ie, make use of MAPL functionality)
  !  and also maintain GEOSldas functionality.
  
  integer, parameter :: MAPL_Land_ExcludeFromDomain = 1100
  
contains
  
  subroutine LDAS_read_til_file( tile_file, catch_file, tile_grid_g, tile_coord_land, f2g )
    
    ! read land tile information from *.til file
    !
    ! This is the GEOSldas version of the LDASsa subroutine read_til_file().
    !
    ! inputs:
    !  tile_file       : *.til tile definition file (full path + name)
    !  catch_file      : catchment.def file         (full path + name) 
    !
    ! outputs:
    !  tile_grid_g     : parameters of tile definition grid
    !  tile_coord_land : coordinates of tiles (see tile_coord_type),
    !                    implemented as pointer which is allocated in 
    !                    this subroutine
    !                    NOTE: number of *land* tiles can be diagnosed with size(tile_coord)
    ! optional:
    !    f2g           : the full domain id to the global id
    !
    ! "tile_id" is no longer read from *.til file and is now set in this 
    ! subroutine to match order of tiles in *.til file
    ! - reichle, 22 Aug 2013
    !
    ! improved documentation of bug in some EASE *.til files (header says N_grid=1 but has two grid defs)
    ! and minor clean-up
    ! - reichle,  2 Aug 2020
    !
    ! -------------------------------------------------------------

    implicit none

    character(*),                                          intent(in)   :: tile_file
    character(*),                                          intent(in)   :: catch_file
    type(grid_def_type),                                   intent(inout):: tile_grid_g
    type(tile_coord_type), dimension(:),           pointer              :: tile_coord_land ! out
    integer,               dimension(:), optional, pointer              :: f2g ! out

    ! locals
    type(tile_coord_type), dimension(:), allocatable :: tile_coord
    integer,               dimension(:), allocatable :: f2g_tmp  ! out

    real    :: ease_cell_area
    integer :: i, N_tile, N_grid,tmpint1, tmpint2, tmpint3, tmpint4
    integer :: i_indg_offset, j_indg_offset, col_order
    integer :: N_tile_land, n_lon, n_lat
    logical :: ease_grid
    integer :: typ,k,fid
    
    character(200) :: tmpline,gridname
    character(300) :: fname

    character(len=*), parameter :: Iam = 'LDAS_read_til_file'

    ! ---------------------------------------------------------------

    i_indg_offset = 0
    j_indg_offset = 0
    
    ! read *.til file header 
    
    if (logit) write (logunit,'(400A)') trim(Iam), '(): reading from ' // trim(tile_file)
    
    open (10, file=trim(tile_file), form='formatted', action='read')
    
    read (10,*) N_tile      ! number of all tiles in *.til file, incl non-land types
    read (10,*) N_grid          
    read (10,*) gridname        
    read (10,*) n_lon
    read (10,*) n_lat

    ! NOTE:
    !  There is a bug in at least some EASE *.til files through at least Icarus-NLv4.
    !  Affected files state "N_grid=1" in line 2 of the header, but the header still includes
    !  three additional lines for a second grid.
    !  LDAS pre-processing corrects for this bug through subroutine correctEase() in
    !  preprocess_LDAS.F90, which creates a second, corrected version of the *.til file during
    !  ldas_setup.  Here, this corrected *.til file is read!
    
    if(N_grid==2) then
       read (10,*)          ! some string describing ocean grid                   (?)
       read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
       read (10,*)          ! # ocean grid cells in latitude direction (N_j_ocn)  (?)
    endif

    ease_grid = .false.
    col_order = 0
    
    call LDAS_create_grid_g( gridname, n_lon, n_lat,                        &
         tile_grid_g, i_indg_offset, j_indg_offset, ease_cell_area )
    
    if (index(tile_grid_g%gridtype,'EASE')/=0)  ease_grid = .true.  ! 'EASE' and 'EASEv2'
    if (index(tile_grid_g%gridtype,'SiB2')/=0)  col_order=1         ! old bcs
    
    allocate(tile_coord(N_tile))
    allocate(f2g_tmp(N_tile))
    
    i   = 0
    fid = 0
    
    ! WJ notes: i and k are the same---global ids
    !           fid --- num in simulation domain
    
    do k=1,N_tile
       
       read(10,'(A)')  tmpline 
       read(tmpline,*) typ

       ! tile type "MAPL_Land_ExcludeFromDomain" identifies land tiles to exclude
       !  when non-global domain is created

       if (typ==MAPL_Land .or. typ==MAPL_Land_ExcludeFromDomain) then     ! all land

          i=i+1
          tile_coord(i)%tile_id = k

          ! now keep only tiles that are not excluded by way of MAPL_Land_ExcludeFromDomain
          
          if (typ==MAPL_Land) then
             fid=fid+1
             f2g_tmp(fid) = k
          end if
             
          ! Not sure ".or. N_grid==1" will always work in the following conditional.
          ! Some Tripolar grid *.til files may have N_grid=1.
          ! - reichle, 2 Aug 2020
          
          if (ease_grid .or. N_grid==1) then  
             
             ! EASE grid til file has fewer columns 
             ! (excludes "tile_id", "frac_pfaf", and "area")
             
             read (tmpline,*)                     &
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

                read (tmpline,*)                  &
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
                
                read (tmpline,*)                  &
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

          ! exit if not land
          
          if (logit) then
             write (logunit,*) 'WARNING: Encountered first non-land tile in *.til file.'
             write (logunit,*) '         Stop reading *.til file under the assumption that'
             write (logunit,*) '           land tiles are first in *.til file.'
             write (logunit,*) '         This is NOT a safe assumption beyond Icarus-NLv[x] tile spaces!!'
          end if
          
          exit ! assuming land comes first in the til file
          
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
    ! if elevation info is still needed, read *gridded* elevation data (check only first tile!)
    
    ! gridded elevation file is NOT available for EASE grids, where elevation information
    !  is in catchment.def file
    
    if ( abs(tile_coord_land(1)%elev-nodata_generic)<nodata_tol_generic ) then
       
       i=index(catch_file,'/clsm/')
       fname = catch_file(1:i)//'topo_DYN_ave_*.data'
       call Execute_command_line('ls '//trim(fname) // ' >topo_DYN_ave.file')
       open(10,file='topo_DYN_ave.file', action='read')
       fname= ''
       read(10,'(A)') fname
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
    
  end subroutine LDAS_read_til_file

  ! *************************************************************************************
  
  subroutine read_grid_elev( fname, tile_grid, N_tile, tile_coord )

    ! read gridded elevation file (for GEOS-5 discretizations; NOT available
    ! for EASE grids, where elevation information is in catchment.def file)
    
    ! reichle,  8 Dec 2011: bug fix -- bin elev data is stored in single record

    implicit none
    
    character(*),        intent(in) :: fname
    
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
    
    integer,             intent(in) :: N_tile
    
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

  ! **********************************************************************
  
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
    
    integer,      intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! inout
    
    ! locals
    
    integer :: i, istat, tmpint1, sweep
    
    integer, dimension(N_tile)  :: tmp_tileid, tmp_pfaf
    
    character(len=*), parameter :: Iam = 'read_catchment_def'
    character(len=400)          :: err_msg
    
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
    
    if ( any(tile_coord(1:N_tile)%tile_id/=tmp_tileid) .or.          &
         any(tile_coord(1:N_tile)%pfaf   /=tmp_pfaf)         ) then

       err_msg = 'tile_coord_file and catchment_def_file mismatch. (2)'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if
    
    ! -----------------------------------------------------------------
    
  end subroutine read_catchment_def
  
end module preprocess_ldas_subs

! ==================================== EOF =====================================
