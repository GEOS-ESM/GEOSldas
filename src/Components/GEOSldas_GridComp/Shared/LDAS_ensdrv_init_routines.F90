#include "MAPL_Generic.h"

module LDAS_ensdrv_init_routines
  
  ! collection of LDASsa initialization subroutines, including some I/O routines
  !  that are also used for output
  !
  ! (originally in clsm_ensdrv_drv_routines.F90)
  !
  ! reichle, 22 Aug 2014

  use ESMF
  use GEOS_MOD

  use LDAS_ensdrv_Globals,           ONLY:     &
       log_master_only,                           &
       logunit,                                   &
       logit,                                     &
       nodata_generic
       
  use MAPL_ConstantsMod,                ONLY:     &
       Tzero  => MAPL_TICE

  use MAPL_BaseMod,                     ONLY:     &
       NTYPS  => MAPL_NumVegTypes

  use LDAS_TileCoordType,                 ONLY:     &
       tile_coord_type,                           &
       grid_def_type,                             &
       operator (==),                              &
       N_cont_max,                                &
       io_tile_coord_type,                        &
       io_grid_def_type
  
  use LDAS_DateTimeMod,                   ONLY:     &
       date_time_type,                            &
       get_dofyr_pentad

  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename,                           &
       is_in_list,                                &
       is_in_domain,                              &
       open_land_param_file,                      &
       word_count
       
!  use clsm_ensdrv_drv_routines,         ONLY:     &
!       check_cat_progn

  use LDAS_TileCoordRoutines,              ONLY:     & 
       is_cat_in_box,                             &
       !reorder_tiles,                             &
       get_tile_grid,                             &
       read_til_file
  
  use LDAS_ExceptionsMod,                  ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use catch_types,                      ONLY:     &
       cat_param_type,                            &
       cat_progn_type,                            &
       N_gt,                                      &
       N_snow

  use LDAS_ensdrv_mpi,                  ONLY:     &
       master_proc

  use ESMF
  use MAPL_Mod

  implicit none

  private
 
  public :: add_domain_to_path
  public :: domain_setup
  public :: read_cat_param
  public :: clsm_ensdrv_get_command_line
  public :: io_domain_files

  !integer ,parameter :: N_gt=6, N_snow=3

  character(10), private :: tmpstring10
  character(40), private :: tmpstring40

contains
  
  ! ********************************************************************


  character(200) function add_domain_to_path( pathname, exp_domain )
    
    ! make sure "exp_domain" is always added to "pathname" in the same way
    ! (so that comparison of strings does not depend on extra slashes)
    ! - reichle, 2 Apr 2014
    
    implicit none
    
    character(200) :: pathname
    character(40)  :: exp_domain
    
    add_domain_to_path = trim(pathname)  // '/' // trim(exp_domain) // '/'
    
  end function add_domain_to_path
    
  ! ****************************************************************
  
  subroutine domain_setup(                                             &
       N_cat_global, tile_coord_global,                                &
       tile_grid_g,                                                    &
       black_path, black_file, white_path, white_file,                 &
       work_path, exp_domain, exp_id,                                  &
       minlon, minlat, maxlon, maxlat,                                 &
       N_cat_domain, d2g, tile_coord, tile_grid_d,                     &
       N_catd_cont )
    
    ! Set up modeling domain and determine index vectors mapping from the
    ! domain to global catchment space.
    ! Determine actual bounding box for domain.
    ! Also return tile_coord for domain and tile_grid_d for domain.
    !
    ! -----------------------
    !
    ! The domain is set up using (if present) a "blacklist" of catchments 
    ! to be excluded, a "whitelist" (if present) of catchments to be included,
    ! and the bounding box of a rectangular "zoomed" area (as specified
    ! in driver_inputs). 
    !
    ! order of precedence:
    !  exclude blacklisted catchments and catchments outside continent
    !  include whitelisted catchments or catchments within rectangular domain
    !
    ! input: 
    !
    ! input/output:
    !  tile_grid_g :        def of global tile definition grid
    !  minlon, maxlon, etc: coordinates of bounding box of domain units as 
    !                        in tile_coord file, that is longitude -180:180, 
    !                        latitude -90:90
    !
    ! output:
    !  N_cat_domain = number of catchments in zoomed/whitelisted domain
    !                   (for which model integration is conducted) 
    !  d2g          = index from domain to global tiles
    !  tile_coord_d = tile_coord vector for domain 
    !  tile_grid_d  = def of smallest subgrid of global tile_grid_g that contains
    !                  all catchments (or tiles) in the domain (tile_grid_d%i_offg, 
    !                  tile_grid_d%j_offg are offsets in indices between tile_grid_g
    !                  and tile_grid_d)
    !  N_catd_cont  = number of catchments of (full) domain on each continent 
    !
    !
    ! - reichle, May  7, 2003
    ! - reichle, Nov  7, 2003 - computation of bounding box of actual domain
    ! - reichle, Jul 20, 2004 - fixed initialization of min_min_lon etc 
    ! - reichle, May 11, 2005 - minor output path changes for redesign
    ! - reichle, May 16, 2005 - add output of tile_grid_d
    ! - reichle, Aug 18, 2005 - reinstated minlon, maxlon, minlat, maxlat
    ! - reichle, Jul 23, 2010 - major overhaul
    !
    ! ----------------------------------------------------------
    
    implicit none
        
    integer, intent(in) :: N_cat_global
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord_global ! input
    
    type(grid_def_type),   intent(in)            :: tile_grid_g
    
    character(*),        intent(in)            :: black_path, white_path
    character(*),        intent(in)            :: black_file, white_file
    
    character(*),        intent(in)            :: work_path
    
    character(*),         intent(in)            :: exp_domain, exp_id
    
    real,                  intent(in)            :: minlon, minlat  ! from nml inputs
    real,                  intent(in)            :: maxlon, maxlat  ! from nml inputs
    
    integer,               intent(out)           :: N_cat_domain
    
    integer,               dimension(:), pointer :: d2g        ! output
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! output
    
    type(grid_def_type),   intent(out)           :: tile_grid_d

    integer, dimension(N_cont_max), intent(out)  :: N_catd_cont

    ! locals 
    
    integer :: n, this_tileid, this_catpfaf, N_black, N_white, indomain, rc
    
    integer, dimension(N_cat_global) :: blacklist, whitelist, tmp_d2g
    
    real :: this_minlon, this_minlat, this_maxlon, this_maxlat
    
    logical :: this_cat_black, this_cat_white, this_cat_in_box
    
    integer :: this_i_indg, this_j_indg
    
    type(grid_def_type) :: tmp_grid_def
    logical :: c3_grid 
    character(300) :: fname

    character(len=*), parameter :: Iam = 'domain_setup'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------
    
    if (logit) write (logunit,*) 'Setting up domain: '
    if (logit) write (logunit,*)
    
    ! ------------------------------------------------------------
    !
    ! try reading *domain.txt, *tilecoord.txt, and *tilegrids.txt files 

    call io_domain_files( 'r', work_path, exp_id, &
         N_cat_domain, d2g, tile_coord, tmp_grid_def, tile_grid_d, rc )
    
    if (rc==0) then        ! read was successful
       
       ! minimal consistency check 
       
       if  (.not. tile_grid_g==tmp_grid_def) then
          err_msg = 'existing domain files inconsistent with ' // &
               'global tile_grid_g from tile_coord_file'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       write (logunit,*) 'Domain successfully defined from existing files above.'
       write (logunit,*)

       ! assume that tile_coord and d2g have been reordered previously,
       ! just get number of tiles on each continent 
       
       !call reorder_tiles(.false., N_cat_domain, tile_coord, d2g, N_catd_cont)        

    else           
       
       print*, "Creating domain..., reading white and black lists if there have ones..." 
       ! ------------------------------------------------------------
       !
       ! load blacklist: catchments listed in this file will be excluded
       
       fname = trim(black_path) // '/' // trim(black_file)

       call read_black_or_whitelist(N_cat_global, fname, blacklist, N_black) 
       
       ! load whitelist: catchments listed in this file will be included
       ! (unless excluded via blacklist)
       
       fname = trim(white_path) // '/' // trim(white_file)
       
       call read_black_or_whitelist(N_cat_global, fname, whitelist, N_white) 
       ! -----------------
       !
       ! find and count catchments that are in the domain
       
       c3_grid = .false.
       if(index(tile_grid_g%gridtype,"c3")/=0) c3_grid = .true.

       indomain    = 0     ! initialize

       do n=1,N_cat_global
          
          this_tileid  = tile_coord_global(n)%tile_id
          
          if( .not. c3_grid) then
             this_minlon  = tile_coord_global(n)%min_lon
             this_minlat  = tile_coord_global(n)%min_lat
             this_maxlon  = tile_coord_global(n)%max_lon
             this_maxlat  = tile_coord_global(n)%max_lat
          else ! c3 grid can straddle the lat-lon
             this_minlon  = tile_coord_global(n)%com_lon
             this_minlat  = tile_coord_global(n)%com_lat
             this_maxlon  = tile_coord_global(n)%com_lon
             this_maxlat  = tile_coord_global(n)%com_lat
          endif
          
          
          this_cat_black = is_in_list(N_black,blacklist(1:N_black),this_tileid)
          this_cat_white = is_in_list(N_white,whitelist(1:N_white),this_tileid)
          
          this_cat_in_box =                                                     &
               is_cat_in_box(this_minlon,this_minlat,this_maxlon,this_maxlat,   &
               minlon, minlat, maxlon, maxlat        )
           
          if (is_in_domain(                                                 &
               this_cat_black, this_cat_white, this_cat_in_box ))  then
             
             indomain = indomain + 1
             tmp_d2g(indomain) = n

          end if
          
       end do
       
       N_cat_domain    = indomain
       
       if (N_cat_domain .eq. 0) then
          err_msg = 'No catchments found in domain'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       else
          if (logit) then
             write (logunit,*) 'Number of catchments in domain = ', N_cat_domain 
             write (logunit,*)
          end if
       end if
       
       ! -------------------------------------------------------------------
       !
       ! assemble d2g, tile_coord, tile_grid_d
       
       allocate(d2g(       N_cat_domain))
       allocate(tile_coord(N_cat_domain))
       
       d2g(1:N_cat_domain) = tmp_d2g(1:N_cat_domain)
       
       tile_coord          = tile_coord_global(d2g)
              
       ! reorder tile_coord and d2g to fix dateline issue and for domain decomposition
       ! No need to reorder in GEOSldas because it devides domain among processors first and then assign tiles
       ! In old LDAS, it dedvides the tile vector among processors.

        !call reorder_tiles(.true., N_cat_domain, tile_coord, d2g, N_catd_cont) 
       
       ! finalize extent of actual domain:
       !  determine smallest subgrid of tile_grid_d that contains all
       !  catchments/tiles in domain
       
       call get_tile_grid( N_cat_domain, tile_coord, tile_grid_g, tile_grid_d ) 

       ! output domain files
       
       tmp_grid_def = tile_grid_g  ! cannot use intent(in) tile_grid_g w/ io_domain_files
      
       call io_domain_files( 'w', work_path, exp_id, &
            N_cat_domain, d2g, tile_coord, tmp_grid_def, tile_grid_d, rc )
       
    end if   ! domain/tilecoord/tilegrids files exist
    
    ! output extent of domain and tile_grid_d to logunit
    
    if (logit) write (logunit,*) 'Actual extent of domain grid:'
    if (logit) write (logunit,*) 'min lon = ', tile_grid_d%ll_lon
    if (logit) write (logunit,*) 'max lon = ', tile_grid_d%ur_lon
    if (logit) write (logunit,*) 'min lat = ', tile_grid_d%ll_lat
    if (logit) write (logunit,*) 'max lat = ', tile_grid_d%ur_lat
    if (logit) write (logunit,*) 
    
    tmpstring40 = 'tile_grid_d'
    
    if (logit) call io_grid_def_type('w', logunit, tile_grid_d, tmpstring40)
    print*,"Finish domain setup..." 
  end subroutine domain_setup

  ! **********************************************************************  
  
  subroutine io_domain_files( action, work_path, exp_id,               &
       N_cat_domain, d2g, tile_coord, tile_grid_g, tile_grid_d, rc )
    
    ! reichle, 23 July 2010
    ! reichle,  7 Jan  2014 - changed tile_coord and tile_grids I/O to binary

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


    !fname = trim(work_path)//'/'//trim(exp_id)//trim(file_tag)//trim(file_ext)
 
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

   ! fname = trim(work_path)//'/'//trim(exp_id)//trim(file_tag)//trim(file_ext)
    
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

   ! fname = trim(work_path)//'/'//trim(exp_id)//trim(file_tag)//trim(file_ext)

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

  subroutine read_cat_param(                                                    &
       N_catg, N_catf, f2g, tile_coord_f, dzsf, veg_path, soil_path, top_path,  &
       cp )

    ! Reads soil properties and topographic parameters from global files 
    ! and extracts data for the (full) domain.
    !    
    ! Additional parameters are derived from the ones that have been read 
    ! from files.
    !
    !  cp = cat_param_f
    !
    ! reichle, 12 May 2003 
    ! reichle,  6 Jun 2005 - adapted to read "SiB2_V2" parameters 
    ! reichle,  5 Apr 2013 - removed alb_param_type fields "sc_albvr", "sc_albnr" 
    ! reichle, 25 Jul 2013 - removed LAI, GRN, and albedo inputs, renamed subroutine
    !                        from "read_land_parameters()" to "read_cat_param()"
    ! reichle, 16 Nov 2015 - read static (JPL) veg height from boundary condition file
    ! 
    ! -------------------------------------------------------------------

    implicit none

    integer,                                    intent(in)  :: N_catg, N_catf

    type(tile_coord_type), dimension(:),   pointer :: tile_coord_f ! intent(in)

    real,                                       intent(in)  :: dzsf

    integer,               dimension(N_catf),   intent(in)  :: f2g

    character(*),                             intent(in)  :: veg_path
    character(*),                             intent(in)  :: soil_path
    character(*),                             intent(in)  :: top_path

    type(cat_param_type),  dimension(N_catf),   intent(out) :: cp

    ! local variables

    integer, parameter :: N_search_dir_max =  5
    integer, parameter :: N_col_max        = 18  ! "v15" soil_param.dat had 22 columns 

    character( 80) :: fname
    character(999) :: tmpstr999

    character(100), dimension(N_search_dir_max) :: search_dir

    integer :: n, k, m, dummy_int, dummy_int2, istat, N_search_dir, N_col

    integer, dimension(N_catg)           :: tmpint, tmpint2, tmptileid

    real,    dimension(N_catg,N_col_max) :: tmpreal

    real    :: dummy_real, dummy_real2, z_in_m, term1, term2

    logical :: dummy_logical

    character(len=*), parameter          :: Iam = 'read_cat_param'
    character(len=400)                   :: err_msg

    real,    dimension(NTYPS)            :: VGZ2

    ! legacy vegetation height look-up table (for backward compatibility)
    !
    DATA VGZ2 /35.0, 20.0, 17.0, 0.6, 0.5, 0.6/      ! Dorman and Sellers (1989)

    ! ---------------------------------------------------------------------

    if (logit) write (logunit,*) 'reading Catchment model parameters'
    if (logit) write (logunit,*)

    ! -----------------------------

    ! Vegetation class

    if (logit) write (logunit,*) 'Reading vegetation class and, if available, height'

    fname = '/mosaic_veg_typs_fracs'

    N_search_dir = 2         ! specify sub-dirs of veg_path to search for file "fname"

    search_dir(1) = 'clsm'
    search_dir(2) = 'VEGETATION-GSWP2'

    ! find out how many columns are in the (formatted) file

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, trim(veg_path), search_dir)

    read(10,'(a)') tmpstr999    ! read first line

    close(10, status='keep')

    ! count words in first line (delimited by space)

    N_col = word_count( tmpstr999 )

    ! read parameters

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, trim(veg_path), search_dir)

    tmptileid = 0

    tmpreal   = nodata_generic

    select case (N_col)

    case (6)

       ! legacy vegetation height from look-up table

       if (logit) write (logunit,*) 'Using vegetation height look-up table'

       do n=1,N_catg

          read (10,*) tmptileid(n), dummy_int, tmpint(n)

       end do

       if ( (any(tmpint<1)) .or. (any(tmpint>NTYPS)) ) then

          err_msg = 'veg type (class) exceeds allowed min/max'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end if

       do n=1,N_catg

          tmpreal(n,1) = VGZ2( tmpint(n) )

       end do


    case (7,8)

       ! vegetation height from boundary condition file

       if (logit) write (logunit,*) 'reading vegetation height from file'

       do n=1,N_catg

          ! 7-th column contains veg height in m
          ! 8-th column contains ASCAT z0 values (IGNORED for now, reichle, 31 Oct 2017)

          read (10,*) tmptileid(n), dummy_int, tmpint(n),         &
               dummy_int2, dummy_real, dummy_real2, tmpreal(n,1)

       end do

    case default

       err_msg = 'unknown number of columns in ' // trim(fname)
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end select

    close (10,status='keep')

    if (logit) write (logunit,*) 'done reading'
    if (logit) write (logunit,*)

    do k=1,N_catf

       ! this check works only for "SiB2_V2" and newer versions

       if (tile_coord_f(k)%tile_id/=tmptileid(f2g(k))) then
          err_msg = 'something wrong with veg parameters'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       cp(k)%vegcls  = tmpint( f2g(k)  )
       cp(k)%veghght = tmpreal(f2g(k),1)

    end do

    ! -----------------------------------

    ! Soil parameters, surface layer time scales, and topographical parameters

    N_search_dir = 2         ! specify sub-dirs of path to search for file "fname"

    search_dir(1) = 'clsm'
    search_dir(2) = '.'

    ! ---------------------
    !
    ! Soil parameters

    if (logit) write (logunit,*) 'Reading soil parameters'

    fname = '/soil_param.dat'

    ! find out how many columns are in the (formatted) file

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, soil_path, search_dir)

    read(10,'(a)') tmpstr999    ! read first line

    close(10, status='keep')

    ! count words in first line (delimited by space)

    N_col = word_count( tmpstr999 )

    ! read parameters

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, soil_path, search_dir)

    tmptileid = 0

    do n=1,N_catg

       ! "SiB2_V2" version  

       read (10,*) tmptileid(n), dummy_int, tmpint(n), tmpint2(n), &
            (tmpreal(n,m), m=1,N_col-4)

    end do

    close (10,status='keep')

    if (logit) write (logunit,*) 'done reading'
    if (logit) write (logunit,*)

    do k=1,N_catf

       if (tile_coord_f(k)%tile_id/=tmptileid(f2g(k))) then
          err_msg = 'something wrong with soil parameters'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       cp(k)%bee   = tmpreal(f2g(k),1)
       cp(k)%psis  = tmpreal(f2g(k),2)
       cp(k)%poros = tmpreal(f2g(k),3)
       cp(k)%cond  = tmpreal(f2g(k),4)
       cp(k)%wpwet = tmpreal(f2g(k),5)
       cp(k)%dpth  = tmpreal(f2g(k),6)

       cp(k)%soilcls30  = tmpint( f2g(k))
       cp(k)%soilcls100 = tmpint2(f2g(k))

    end do

    ! additional soil parameters from recent versions of "soil_param.dat"
    ! (eg. for use in calibration of the microwave radiative transfer model)
    ! - reichle,  1 Apr 2015

    select case (N_col)

    case (19)

       ! starting with "v16" (De Lannoy et al., 2014, doi:10.1002/2014MS000330),
       !  soil_param.dat has 19 columns

       do k=1,N_catf

          cp(k)%gravel30 = tmpreal(f2g(k), 7)
          cp(k)%orgC30   = tmpreal(f2g(k), 8)
          cp(k)%orgC     = tmpreal(f2g(k), 9)
          cp(k)%sand30   = tmpreal(f2g(k),10)
          cp(k)%clay30   = tmpreal(f2g(k),11)
          cp(k)%sand     = tmpreal(f2g(k),12)
          cp(k)%clay     = tmpreal(f2g(k),13)
          cp(k)%wpwet30  = tmpreal(f2g(k),14)
          cp(k)%poros30  = tmpreal(f2g(k),15)

       end do

    case default

       do k=1,N_catf

          cp(k)%gravel30 = nodata_generic
          cp(k)%orgC30   = nodata_generic
          cp(k)%orgC     = nodata_generic
          cp(k)%sand30   = nodata_generic
          cp(k)%clay30   = nodata_generic
          cp(k)%sand     = nodata_generic
          cp(k)%clay     = nodata_generic
          cp(k)%wpwet30  = nodata_generic
          cp(k)%poros30  = nodata_generic

       end do

    end select

    ! ------------------------------------

    ! Surface layer timescales

    if (logit) write (logunit,*) 'Reading surface layer timescales atau/btau'

    fname = '/tau_param.dat'

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, soil_path, search_dir)

    tmptileid = 0

    do n=1,N_catg

       read (10,*) tmptileid(n), dummy_int, (tmpreal(n,m), m=1,4)

    end do

    close (10,status='keep')

    if (logit) write (logunit,*) 'done reading'
    if (logit) write (logunit,*)

    do k=1,N_catf

       ! this check works only for "SiB2_V2" version

       if (tile_coord_f(k)%tile_id/=tmptileid(f2g(k))) then
          err_msg = 'something wrong with tau parameters'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       ! select atau and btau depending on surface layer depth 

       if     (abs(dzsf-20.)<1e-4 ) then     ! use atau2, btau2

          cp(k)%atau = tmpreal(f2g(k),1)
          cp(k)%btau = tmpreal(f2g(k),2)

       elseif (abs(dzsf-50.)<1e-4 ) then     ! use atau5, btau5

          cp(k)%atau = tmpreal(f2g(k),3)
          cp(k)%btau = tmpreal(f2g(k),4)

       else

          err_msg = 'unknown value for dzsf'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end if

    end do

    ! make sure atau and btau are not unphysical

    if (any(cp%atau<=0)) then
       err_msg = 'unphysical atau value(s)'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    if (any(cp%btau<=0)) then
       err_msg = 'unphysical btau value(s)'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if


    ! ------------------------------------

    ! Topographical parameters

    if (logit) write (logunit,*) 'Reading topo parameters (ar)'

    fname = '/ar.new'

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, top_path, search_dir)

    tmptileid = 0

    do n=1,N_catg

       read (10,*) tmptileid(n), dummy_int, (tmpreal(n,m), m=1,12)

    end do

    close (10,status='keep')

    if (logit) write (logunit,*) 'done reading'
    if (logit) write (logunit,*)

    do k=1,N_catf

       ! this check works only for "SiB2_V2" version

       if (tile_coord_f(k)%tile_id/=tmptileid(f2g(k))) then
          err_msg = 'something wrong with ar parameters'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       cp(k)%gnu   = tmpreal(f2g(k),1)
       cp(k)%ars1  = tmpreal(f2g(k),2)
       cp(k)%ars2  = tmpreal(f2g(k),3)
       cp(k)%ars3  = tmpreal(f2g(k),4)
       cp(k)%ara1  = tmpreal(f2g(k),5)
       cp(k)%ara2  = tmpreal(f2g(k),6)
       cp(k)%ara3  = tmpreal(f2g(k),7)
       cp(k)%ara4  = tmpreal(f2g(k),8)
       cp(k)%arw1  = tmpreal(f2g(k),9)
       cp(k)%arw2  = tmpreal(f2g(k),10)
       cp(k)%arw3  = tmpreal(f2g(k),11)
       cp(k)%arw4  = tmpreal(f2g(k),12)

    end do

    ! --------------------

    if (logit) write (logunit,*) 'Reading topo parameters (bf)'

    fname = '/bf.dat'

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, top_path, search_dir)

    tmptileid = 0

    do n=1,N_catg

       read (10,*) tmptileid(n), dummy_int, (tmpreal(n,m), m=1,4)

    end do

    close (10,status='keep')

    if (logit) write (logunit,*) 'done reading'
    if (logit) write (logunit,*)

    do k=1,N_catf

       ! this check works only for "SiB2_V2" version

       if (tile_coord_f(k)%tile_id/=tmptileid(f2g(k))) then
          err_msg = 'something wrong with bf parameters'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       ! ---------

       if (cp(k)%gnu/=tmpreal(f2g(k),1)) then
          err_msg = 'land(): something wrong with gnu'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       cp(k)%bf1  = tmpreal(f2g(k),2)
       cp(k)%bf2  = tmpreal(f2g(k),3)
       cp(k)%bf3  = tmpreal(f2g(k),4)

    end do

    ! --------------------

    if (logit) write (logunit,*) 'Reading topo parameters (ts)'

    fname = '/ts.dat'

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, top_path, search_dir)

    tmptileid = 0

    do n=1,N_catg

       read (10,*) tmptileid(n), dummy_int, (tmpreal(n,m), m=1,5)

    end do

    close (10,status='keep')

    if (logit) write (logunit,*) 'done reading'
    if (logit) write (logunit,*)
    do k=1,N_catf

       ! this check works only for "SiB2_V2" version

       if (tile_coord_f(k)%tile_id/=tmptileid(f2g(k))) then
          err_msg = 'something wrong with ts parameters'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       ! -------

       if (cp(k)%gnu/=tmpreal(f2g(k),1)) then
          err_msg = 'land(): something wrong with gnu'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       cp(k)%tsa1  = tmpreal(f2g(k),2)
       cp(k)%tsa2  = tmpreal(f2g(k),3)
       cp(k)%tsb1  = tmpreal(f2g(k),4)
       cp(k)%tsb2  = tmpreal(f2g(k),5)

    end do

    ! ---------------------------------------------------------------------

    if (logit) write (logunit,*) 'computing derived land surface parameters...'
    if (logit) write (logunit,*)

    do k=1,N_catf

       ! Three soil depths for soil moisture model: 
       !
       ! dzsf: surface layer
       ! dzrz: root zone        -> water capacity of the root zone
       ! dzpr: unsaturated zone -> approx depth-to-bedrock
       !
       ! NOTE: Units of dz** are [mm] while excess/deficits from catchment()
       !       are in SI units (ie kg/m^2) or loosely speaking, in mm of water.
       !       In other words, density of water (1000 kg/m^3) is built 
       !       into dz** (reichle, 5 Feb 04).

       cp(k)%dzsf = dzsf

       cp(k)%dzrz = 1000.

       ! changed re-setting of dzrz back to earlier value because
       ! Sarith parameters are in fact consistent that the earlier version
       ! reichle, 12 Sep 2007
       !
       ! cp(k)%dzpr = max(1500., cp(k)%dpth)
       !
       ! previously, root zone depth ranged from .75m to 1m, which
       ! is inconsistent with subroutine catchment(), where root
       ! zone depth is hard-wired to 1m, and with the time scale
       ! parameters, that have been derived for 1m root zone depth
       ! (THE LATTER IS IN FACT *NOT* TRUE - reichle, 12 Sep 2007)
       ! - reichle, 30 May 2003

       cp(k)%dzpr = max(1000., cp(k)%dpth)

       if (cp(k)%dzrz > 0.75*cp(k)%dzpr)     cp(k)%dzrz = 0.75*cp(k)%dzpr

       ! soil storages

       cp(k)%vgwmax = cp(k)%poros*cp(k)%dzrz

       z_in_m = cp(k)%dzpr/1000.

       term1 = -1.+((cp(k)%psis-z_in_m)/cp(k)%psis)**((cp(k)%bee-1.)/cp(k)%bee)

       term2 = cp(k)%psis*cp(k)%bee/(cp(k)%bee-1)

       cp(k)%cdcr1 = 1000.*cp(k)%poros*(z_in_m-(-term2*term1))

       cp(k)%cdcr2 = (1.-cp(k)%wpwet)*cp(k)%poros*cp(k)%dzpr

       ! soil depths for ground temperature model

       if (N_gt/=6) then

          write (tmpstring10,*) N_gt

          err_msg = 'using N_gt = ' // trim(tmpstring10) // &
               'but only 6 layer depths are specified.'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end if

       cp(k)%dzgt(1) =  0.0988
       cp(k)%dzgt(2) =  0.1952
       cp(k)%dzgt(3) =  0.3859
       cp(k)%dzgt(4) =  0.7626
       cp(k)%dzgt(5) =  1.5071
       cp(k)%dzgt(6) = 10.0000

    end do

  end subroutine read_cat_param

  ! *************************************************************************

  subroutine read_VEG_Height(                                                    &
       N_catg, veg_path, V_HEIGHT )

    ! addapted from read_cat_param

    implicit none
    
    integer,                                    intent(in)  :: N_catg
    character(*),                             intent(in)  :: veg_path
    real, dimension(N_catg),intent(inout) :: V_HEIGHT
       
    character( 80) :: fname
    character(999) :: tmpstr999
    
    character(100), dimension(2) :: search_dir

    integer :: n, k, m, dummy_int, dummy_int2, istat, N_search_dir, N_col
    
    integer, dimension(N_catg)           :: tmpint, tmpint2, tmptileid
    
    real,    dimension(N_catg,7) :: tmpreal
    
    real    :: dummy_real, dummy_real2, z_in_m, term1, term2

    logical :: dummy_logical

    character(len=*), parameter          :: Iam = 'read_Veg_Hight'
    character(len=400)                   :: err_msg

    real,    dimension(NTYPS)            :: VGZ2

    ! legacy vegetation height look-up table (for backward compatibility)
    !
    DATA VGZ2 /35.0, 20.0, 17.0, 0.6, 0.5, 0.6/      ! Dorman and Sellers (1989)
    
    ! ---------------------------------------------------------------------
    
    if (logit) write (logunit,*) 'reading VEG_HEIGHT'
    if (logit) write (logunit,*)
    
    ! -----------------------------
    
    ! Vegetation class

    if (logit) write (logunit,*) 'Reading vegetation class and, if available, height'
    
    fname = '/mosaic_veg_typs_fracs'

    N_search_dir = 2         ! specify sub-dirs of veg_path to search for file "fname"
    
    search_dir(1) = 'clsm'
    search_dir(2) = 'VEGETATION-GSWP2'

    ! find out how many columns are in the (formatted) file
    
    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, veg_path, search_dir)
        
    read(10,'(a)') tmpstr999    ! read first line
    ! count words in first line (delimited by space)
    N_col = word_count( tmpstr999 )

    ! get line number or the real N_catg
    n =1
    do while (.true.) 
       read(10,*,iostat= istat) tmpstr999
       if(IS_IOSTAT_END(istat)) exit
       n=n+1  
    enddo

    if(n /= N_catg) stop " Please don't add vegheight to REGIONAOL veg restart"

    close(10, status='keep')
    
    
    ! read parameters

    istat = open_land_param_file(                                                 &
         10, .true., dummy_logical, N_search_dir, fname, veg_path, search_dir)

    tmptileid = 0

    tmpreal   = nodata_generic
    
    select case (N_col)
       
    case (6)

       ! legacy vegetation height from look-up table
       
       if (logit) write (logunit,*) 'Using vegetation height look-up table'

       do n=1,N_catg
          
          read (10,*) tmptileid(n), dummy_int, tmpint(n)
          
       end do

       if ( (any(tmpint<1)) .or. (any(tmpint>NTYPS)) ) then
          
          err_msg = 'veg type (class) exceeds allowed min/max'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
       do n=1,N_catg

          tmpreal(n,1) = VGZ2( tmpint(n) )

       end do

    case (7)

       ! vegetation height from boundary condition file
       
       if (logit) write (logunit,*) 'reading vegetation height from file'

       do n=1,N_catg
          
          ! 7-th column contains veg height in m
          
          read (10,*) tmptileid(n), dummy_int, tmpint(n),         &
               dummy_int2, dummy_real, dummy_real2, tmpreal(n,1)
          
       end do
       
    case default
       
       err_msg = 'unknown number of columns in ' // trim(fname)
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
    close (10,status='keep')
    
    if (logit) write (logunit,*) 'done reading'
    if (logit) write (logunit,*) 
    
    V_HEIGHT = tmpreal(:,1)
      
  end subroutine read_VEG_Height
  
  ! **********************************************************************
  
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
  
  ! ****************************************************************
  
  subroutine clsm_ensdrv_get_command_line(   &
       driver_inputs_path,      &
       driver_inputs_file,      &
       start_time,              &
       end_time,                &
       resolution,              &
       exp_domain,              &
       exp_id,                  &
       work_path,               &
       run_path,                &
       restart_path,            &
       restart_domain,          &
       restart_id,              &
       tile_coord_path,         &
       tile_coord_file,         &
       catchment_def_path,      &
       catchment_def_file,      &
       met_tag,                 &
       met_path,                &
       force_dtstep,            &
       restart,                 &
       spin,                    &
       ens_prop_inputs_path,    &
       ens_prop_inputs_file,    &
       N_ens,                   &
       first_ens_id             &
       )
    
    ! get inputs from command line 
    !
    ! if present, command line arguments overwrite inputs from
    !  driver_inputs or ens_prop_inputs namelist files
    !
    ! command line should look something like
    !
    ! a.out -start_year 1979 -restart true -driver_inputs_file fname.nml 
    !
    ! NOTE: Arguments that are used for assimilation ("ensupd") must be 
    !       listed here explicitly (and will be ignored).  Otherwise it
    !       would not be possible to stop for unknown arguments.
    !
    ! reichle, 29 Aug 02
    ! 
    ! modified for namelist input file path and name
    ! - reichle, 6 May 03
    ! converted to optional arguments, added arguments for EnKF inputs
    ! - reichle, 29 Mar 04
    !
    ! reichle,  2 Aug 2005 - consistency with clsm_ensdrv_get_command_line()
    ! reichle,  6 Mar 2008 - added force_dtstep for DAS/MERRA integration
    !
    ! ----------------------------------------------------------------
    
    implicit none
    
    character(200), intent(inout), optional :: driver_inputs_path
    character(200), intent(inout), optional :: work_path, run_path
    character(200), intent(inout), optional :: restart_path
    
    character(40), intent(inout),  optional :: driver_inputs_file    
    
    type(date_time_type), intent(inout), optional :: start_time
    type(date_time_type), intent(inout), optional :: end_time
    
    character(40), intent(inout), optional :: exp_domain, exp_id, resolution
    character(40), intent(inout), optional :: restart_domain, restart_id

    character(200), intent(inout), optional :: tile_coord_path
    character(200), intent(inout), optional :: catchment_def_path
    character(80),  intent(inout), optional :: tile_coord_file
    character(40),  intent(inout), optional :: catchment_def_file
    
    character(200), intent(inout), optional :: met_path
    character(80),  intent(inout), optional :: met_tag
    
    integer,        intent(inout), optional :: force_dtstep

    logical,        intent(inout), optional :: restart, spin 
    
    character(200), intent(inout), optional :: ens_prop_inputs_path
    character(40),  intent(inout), optional :: ens_prop_inputs_file    
    
    integer, intent(inout), optional :: N_ens, first_ens_id

    ! -----------------------------------------------------------------
    
    integer :: N_args, iargc, i
    
    character(40) :: arg

    logical :: outlog
    
    character(len=*), parameter :: Iam = 'clsm_ensdrv_get_command_line'
    character(len=400) :: err_msg

    !external getarg, iargc
    
    ! -----------------------------------------------------------------

    ! make sure log file has already been opened

    inquire( logunit, opened=outlog )

    ! do not write log output for non-master processes as requested
    
    if ((log_master_only) .and. (.not. master_proc)) outlog = .false.
    
    N_args = iargc()
    
    i=0
    
    do while ( i < N_args )
       
       i = i+1
       
       call getarg(i,arg)
       
       if     ( trim(arg) == '-driver_inputs_path' ) then
          i = i+1
          if (present(driver_inputs_path))  call getarg(i,driver_inputs_path)
          
       elseif ( trim(arg) == '-driver_inputs_file' ) then
          i = i+1
          if (present(driver_inputs_file))  call getarg(i,driver_inputs_file)
          
       elseif ( trim(arg) == '-start_year' ) then
          i = i+1
          call getarg(i,arg)
          if (present(start_time))  read (arg,*) start_time%year
          
       elseif ( trim(arg) == '-start_month' ) then
          i = i+1
          call getarg(i,arg)
          if (present(start_time))  read (arg,*) start_time%month

       elseif ( trim(arg) == '-start_day' ) then
          i = i+1
          call getarg(i,arg)
          if (present(start_time))  read (arg,*) start_time%day

       elseif ( trim(arg) == '-start_hour' ) then
          i = i+1
          call getarg(i,arg)
          if (present(start_time))  read (arg,*) start_time%hour

       elseif ( trim(arg) == '-start_min' ) then
          i = i+1
          call getarg(i,arg)
          if (present(start_time))  read (arg,*) start_time%min

       elseif ( trim(arg) == '-start_sec' ) then
          i = i+1
          call getarg(i,arg)
          if (present(start_time))  read (arg,*) start_time%sec

       elseif ( trim(arg) == '-end_year' ) then
          i = i+1
          call getarg(i,arg)
          if (present(end_time))  read (arg,*) end_time%year
          
       elseif ( trim(arg) == '-end_month' ) then
          i = i+1
          call getarg(i,arg)
          if (present(end_time))  read (arg,*) end_time%month

       elseif ( trim(arg) == '-end_day' ) then
          i = i+1
          call getarg(i,arg)
          if (present(end_time))  read (arg,*) end_time%day

       elseif ( trim(arg) == '-end_hour' ) then
          i = i+1
          call getarg(i,arg)
          if (present(end_time))  read (arg,*) end_time%hour

       elseif ( trim(arg) == '-end_min' ) then
          i = i+1
          call getarg(i,arg)
          if (present(end_time))  read (arg,*) end_time%min

       elseif ( trim(arg) == '-end_sec' ) then
          i = i+1
          call getarg(i,arg)
          if (present(end_time))  read (arg,*) end_time%sec

       elseif ( trim(arg) == '-resolution' ) then
          i = i+1
          if (present(resolution))  call getarg(i,resolution)

       elseif ( trim(arg) == '-exp_domain' ) then
          i = i+1
          if (present(exp_domain))  call getarg(i,exp_domain)

       elseif ( trim(arg) == '-exp_id' ) then
          i = i+1
          if (present(exp_id))  call getarg(i,exp_id)

       elseif ( trim(arg) == '-work_path' ) then
          i = i+1
          if (present(work_path))  call getarg(i,work_path)
          
       elseif ( trim(arg) == '-run_path' ) then
          i = i+1
          if (present(run_path))  call getarg(i,run_path)
          
       elseif ( trim(arg) == '-restart_path' ) then
          i = i+1
          if (present(restart_path))  call getarg(i,restart_path)
          
       elseif ( trim(arg) == '-restart_domain' ) then
          i = i+1
          if (present(restart_domain))  call getarg(i,restart_domain)

       elseif ( trim(arg) == '-restart_id' ) then
          i = i+1
          if (present(restart_id))  call getarg(i,restart_id)
          
       elseif ( trim(arg) == '-tile_coord_path' ) then
          i = i+1
          if (present(tile_coord_path))     call getarg(i,tile_coord_path)
          
       elseif ( trim(arg) == '-tile_coord_file' ) then
          i = i+1
          if (present(tile_coord_file))     call getarg(i,tile_coord_file)
          
       elseif ( trim(arg) == '-catchment_def_path' ) then
          i = i+1
          if (present(catchment_def_path))  call getarg(i,catchment_def_path)
          
       elseif ( trim(arg) == '-catchment_def_file' ) then
          i = i+1
          if (present(catchment_def_file))  call getarg(i,catchment_def_file)

       elseif ( trim(arg) == '-met_path' ) then
          i = i+1
          if (present(met_path))  call getarg(i,met_path)
          
       elseif ( trim(arg) == '-met_tag' ) then
          i = i+1
          if (present(met_tag))  call getarg(i,met_tag)
          
       elseif ( trim(arg) == '-force_dtstep' ) then
          i = i+1
          call getarg(i,arg)
          if (present(force_dtstep))  read (arg,*) force_dtstep
          
       elseif ( trim(arg) == '-restart' ) then
          i = i+1
          call getarg(i,arg)
          if (present(restart))  read (arg,*) restart

       elseif ( trim(arg) == '-spin' ) then
          i = i+1
          call getarg(i,arg)
          if (present(spin))  read (arg,*) spin

       elseif ( trim(arg) == '-ens_prop_inputs_path' ) then
          i = i+1
          if (present(ens_prop_inputs_path))  &
               call getarg(i,ens_prop_inputs_path)
          
       elseif ( trim(arg) == '-ens_prop_inputs_file' ) then
          i = i+1
          if (present(ens_prop_inputs_file)) &
               call getarg(i,ens_prop_inputs_file)
          
       elseif ( trim(arg) == '-N_ens' ) then
          i = i+1
          call getarg(i,arg)
          if (present(N_ens))  read (arg,*) N_ens

       elseif ( trim(arg) == '-first_ens_id' ) then
          i = i+1
          call getarg(i,arg)
          if (present(first_ens_id))  read (arg,*) first_ens_id
          
          
          ! ignore arguments for assimilation, bias, adaptive estimation
          
       elseif ( trim(arg) == '-ens_upd_inputs_path' ) then
          i = i+1
          if (outlog) &
               write (logunit,*) 'clsm_ensdrv_get_command_line(): IGNORING argument = ', &
               trim(arg)
       
       elseif ( trim(arg) == '-ens_upd_inputs_file' ) then
          i = i+1
          if (outlog) &
               write (logunit,*) 'clsm_ensdrv_get_command_line(): IGNORING argument = ', &
               trim(arg)
          
       elseif ( trim(arg) == '-cat_bias_inputs_path' ) then
          i = i+1
          if (outlog) &
               write (logunit,*) 'clsm_ensdrv_get_command_line(): IGNORING argument = ', &
               trim(arg)
          
       elseif ( trim(arg) == '-cat_bias_inputs_file' ) then
          i = i+1
          if (outlog) &
               write (logunit,*) 'clsm_ensdrv_get_command_line(): IGNORING argument = ', &
               trim(arg)

       elseif ( trim(arg) == '-adapt_inputs_path' ) then
          i = i+1
          if (outlog) &
               write (logunit,*) 'clsm_ensdrv_get_command_line(): IGNORING argument = ', &
               trim(arg)
          
       elseif ( trim(arg) == '-adapt_inputs_file' ) then
          i = i+1
          if (outlog) &
               write (logunit,*) 'clsm_ensdrv_get_command_line(): IGNORING argument = ', &
               trim(arg)
          
          ! stop for any other arguments
          
       else
          
          err_msg = 'unknown argument = ' // trim(arg)
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       endif
       
    end do
    
  end subroutine clsm_ensdrv_get_command_line

  ! ***********************************************************************

end module LDAS_ensdrv_init_routines

! *********** EOF ******************************************************
