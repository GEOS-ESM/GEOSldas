#include "MAPL_Generic.h"

module LDAS_ensdrv_init_routines
  
  ! collection of LDASsa initialization subroutines, including some I/O routines
  !  that are also used for output
  !
  ! (originally in clsm_ensdrv_drv_routines.F90)
  !
  ! reichle, 22 Aug 2014
  ! reichle, 22 May 2020 - cleanup

  use ESMF
  use GEOS_MOD

  use LDAS_ensdrv_Globals,              ONLY:     &
       logunit,                                   &
       logit,                                     &
       nodata_generic
       
  use MAPL_BaseMod,                     ONLY:     &
       NTYPS  => MAPL_NumVegTypes

  use LDAS_TileCoordType,               ONLY:     &
       tile_coord_type,                           &
       grid_def_type,                             &
       operator (==),                             &
       io_tile_coord_type,                        &
       io_grid_def_type
  
  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename,                           &
       is_in_list,                                &
       is_in_domain,                              &
       open_land_param_file,                      &
       word_count
       
  use LDAS_TileCoordRoutines,           ONLY:     & 
       is_cat_in_box,                             &
       get_tile_grid
  
  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use catch_types,                      ONLY:     &
       cat_param_type,                            &
       cat_progn_type,                            &
       N_gt,                                      &
       N_snow

  use ESMF
  use MAPL_Mod

  implicit none

  private
 
  public :: domain_setup
  public :: read_cat_param
  public :: io_domain_files

  character(10), private :: tmpstring10
  character(40), private :: tmpstring40

contains
  
  ! ********************************************************************
  
  subroutine domain_setup(                                             &
       N_cat_global, tile_coord_global,                                &
       tile_grid_g,                                                    &
       exclude_path, exclude_file, include_path, include_file,         &
       work_path, exp_domain, exp_id,                                  &
       minlon, minlat, maxlon, maxlat,                                 &
       N_cat_domain, d2g, tile_coord, tile_grid_d )
    
    ! Set up modeling domain and determine index vectors mapping from the
    ! domain to global catchment space.
    ! Determine actual bounding box for domain.
    ! Also return tile_coord for domain and tile_grid_d for domain.
    !
    ! -----------------------
    !
    ! The domain is set up using (if present) an "ExcludeList" of catchments 
    ! to be excluded, an "IncludeList" (if present) of catchments to be included,
    ! and the bounding box of a rectangular "zoomed" area (as specified
    ! in the "exeinp" file used in ldas_setup). 
    !
    ! order of precedence:
    !  1. exclude catchments in ExcludeList
    !  2. include catchments in IncludeList or catchments within rectangular domain
    ! (i.e., catchments in ExcludeList are *always* excluded)
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
    !  N_cat_domain = number of catchments in zoomed domain
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
    
    character(*),          intent(in)            :: exclude_path, include_path
    character(*),          intent(in)            :: exclude_file, include_file
    
    character(*),          intent(in)            :: work_path
    
    character(*),          intent(in)            :: exp_domain, exp_id
    
    real,                  intent(in)            :: minlon, minlat  ! from nml inputs
    real,                  intent(in)            :: maxlon, maxlat  ! from nml inputs
    
    integer,               intent(out)           :: N_cat_domain
    
    integer,               dimension(:), pointer :: d2g        ! output
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! output
    
    type(grid_def_type),   intent(out)           :: tile_grid_d

    ! locals 
    
    integer :: n, this_tileid, this_catpfaf, N_exclude, N_include, indomain, rc
    
    integer, dimension(N_cat_global) :: ExcludeList, IncludeList, tmp_d2g
    
    real :: this_minlon, this_minlat, this_maxlon, this_maxlat
    
    logical :: this_cat_exclude, this_cat_include, this_cat_in_box
    
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

    else           
       
       print*, "Creating domain..., reading IncludeList and ExludeList if present..." 
       ! ------------------------------------------------------------
       !
       ! load ExcludeList: catchments listed in this file will *always* be excluded
       
       fname = trim(exclude_path) // '/' // trim(exclude_file)

       call read_exclude_or_includelist(N_cat_global, fname, ExcludeList, N_exclude) 
       
       ! load IncludeList: catchments listed in this file will be included
       ! (unless excluded via ExcludeList)
       
       fname = trim(include_path) // '/' // trim(include_file)
       
       call read_exclude_or_includelist(N_cat_global, fname, IncludeList, N_include) 
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
          
          
          this_cat_exclude = is_in_list( N_exclude, ExcludeList(1:N_exclude), this_tileid )
          this_cat_include = is_in_list( N_include, IncludeList(1:N_include), this_tileid )
          
          this_cat_in_box =                                                     &
               is_cat_in_box(this_minlon,this_minlat,this_maxlon,this_maxlat,   &
               minlon, minlat, maxlon, maxlat        )
           
          if (is_in_domain(                                                     &
               this_cat_exclude, this_cat_include, this_cat_in_box ))  then
             
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
  
  subroutine read_exclude_or_includelist(N_cat, fname, MyList, N_list) 
    
    ! read numbers/IDs of catchments in MyList (ExcludeList or IncludeList)
    !
    ! format of MyList file: ASCII list of tile IDs
    !
    ! N_list = number of catchments in MyList
    !
    ! reichle, 2 May 2003
    !
    ! --------------------------------------------------------------
    
    implicit none
    
    ! N_cat = max number of catchments allowed in list 
    !         (use N_cat_global when calling this subroutine)
    
    integer,      intent(in)  :: N_cat     
    character(*), intent(in)  :: fname
    
    integer,      intent(out) :: N_list   
    
    integer, dimension(N_cat), intent(out) :: MyList
    
    ! locals
    
    integer :: istat, tmpint
    
    logical :: file_exists
    
    character(len=*), parameter :: Iam = 'read_exclude_or_includelist'
    character(len=400) :: err_msg

    ! -----------------------------------------------------------
    
    N_list = 0
        
    inquire( file=fname, exist=file_exists)
    
    if (file_exists) then
       
       open(10, file=fname, form='formatted', action='read', &
            status='old', iostat=istat)
       
       if (istat==0) then
          
          if (logit) write (logunit,*) &
               'reading ExcludeList or IncludeList from ', trim(fname)
          if (logit) write (logunit,*) 

          do
             read(10,*,iostat=istat) tmpint
             
             if (istat==-1) then
                if (logit) write (logunit,*) ' found ', N_list, ' catchments on list'
                exit
             else if (istat/=0) then
                err_msg = 'read error other than end-of-file'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             else
                N_list = N_list+1
                MyList(N_list) = tmpint
             end if
             
             if (N_list>N_cat) then
                
                write (tmpstring10,*) N_cat
                write (tmpstring40,*) N_list
                
                err_msg = 'N_list=' // trim(tmpstring40) &
                     // ' > N_cat=' // trim(tmpstring10)
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

             end if
          end do
          
          close(10,status='keep')
          
       else
          
          if (logit) write (logunit,*) &
               'could not open ExcludeList or IncludeList file ', trim(fname)
          
       end if
       
    else
       
       if (logit) write (logunit,*) &
            'ExcludeList or IncludeList file does not exist: ', trim(fname)
       
    end if
    
    if (logit) write (logunit,*) 
       
  end subroutine read_exclude_or_includelist
  
  ! ***********************************************************************

end module LDAS_ensdrv_init_routines

! *********** EOF ******************************************************
