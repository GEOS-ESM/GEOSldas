
module preprocess_ldas_routines

  ! collection of subroutines and functions needed for GEOSldas pre-processing
  !
  ! The *.F90 module i was created as follows:
  !
  !  1.) git mv preprocess_ldas.F90 preprocess_ldas_routines.F90   (for best possible git diff)
  !  2.) removed main program from file and put into new file preprocess_ldas.F90
  !  3.) moved additional helper subroutines and functions to here:
  !      - LDAS_read_til_file()                  [from LDAS_TileCoordRoutines.F90]
  !      - read_grid_elev()                      [from LDAS_TileCoordRoutines.F90]
  !      - fix_dateline_bug_in_tilecoord()       [from LDAS_TileCoordRoutines.F90]
  !      - read_catchment_def()                  [from LDAS_TileCoordRoutines.F90]
  !      - is_cat_in_box()                       [from LDAS_TileCoordRoutines.F90]
  !      - domain_setup()                        [from LDAS_ensdrv_init_routines.F90]
  !      - read_exclude_or_includelist()         [from LDAS_ensdrv_init_routines.F90]
  !      - read_cat_param()                      [from LDAS_ensdrv_init_routines.F90]
  !      - is_in_list()                          [from LDAS_ensdrv_functions.F90]
  !      - is_in_domain()                        [from LDAS_ensdrv_functions.F90]
  !      - word_count()                          [from LDAS_ensdrv_functions.F90]
  !      - open_land_param_file()                [from LDAS_ensdrv_functions.F90]
  
  use netcdf
  
  use MAPL
  
  use MAPL_BaseMod,                    ONLY:   &
       NTYPS  => MAPL_NumVegTypes,             &
       MAPL_Land
  
  use MAPL_ConstantsMod,               ONLY:   &
       MAPL_RADIUS                                    ! Earth radius
  
  use LDAS_ensdrv_Globals,             ONLY:   &
       logit,                                  &
       logunit,                                &
       nodata_generic,                         &
       nodata_tol_generic
  
  use LDAS_TileCoordType,              ONLY:   &
       tile_coord_type,                        &
       grid_def_type,                          &
       operator (==),                          &
       io_grid_def_type
  
  use LDAS_TileCoordRoutines,          ONLY:   &
       LDAS_create_grid_g,                     &
       get_minExtent_grid,                     &
       io_domain_files
  
  use nr_ran2_gasdev,                  ONLY:   &
       NRANDSEED
  
  use LDAS_DateTimeMod,                ONLY:   &
       date_time_type
  
  use force_and_cat_progn_pert_types,  ONLY:   &
       N_progn_pert_max,                       &
       N_force_pert_max
  
  use catch_types,                     ONLY:   &
       cat_param_type

  use catch_constants,                 ONLY:   &
       N_gt   => CATCH_N_GT 

  use LDAS_ensdrv_functions,           ONLY:   &
       get_io_filename
  
  use LDAS_ExceptionsMod,              ONLY:   &
       ldas_abort,                             &
       LDAS_GENERIC_ERROR
  
  use gFTL_StringVector
  
  use pFIO
  
  implicit none

  private

  public :: createf2g
  public :: createLocalTilefile
  public :: createLocalBC
  public :: createLocalCatchRestart
  public :: createLocalVegRestart
  public :: createLocalmwRTMRestart
  public :: correctEase
  public :: optimize_latlon
  public :: convert_pert_rst

  character(10), private :: tmpstring10
  character(40), private :: tmpstring40
  
  ! Tile type for land that is to be excluded from the simulation domain.
  ! (GEOSldas  allows for non-global simulations and repeated "zooming"
  !  of the domain while MAPL generally assumes a complete (global) tile
  !  space.  The *_ExcludeFromDomain tile type makes it possible to work
  !  with complete (global) tile files (ie, make use of MAPL functionality)
  !  and also maintain GEOSldas functionality.
  
  integer, parameter :: MAPL_Land_ExcludeFromDomain = 1100
  
contains
  
  ! ********************************************************************

  subroutine createf2g(orig_tile,domain_def,out_path,catch_def_file,exp_id,ymdhm, SURFLAY)
    
    implicit none
    character(*) :: orig_tile
    character(*) :: domain_def
    character(*) :: out_path
    character(*) :: catch_def_file
    character(*) :: exp_id
    character(*) :: ymdhm
    character(*) :: SURFLAY
    
    real :: minlon,maxlon,minlat,maxlat
    character(len=512):: exclude_file,include_file
    character(len=512):: bcs_path
    logical :: file_exist
    logical :: d_exist,c_exist
    
    integer :: n
    
    type(grid_def_type) :: tile_grid_g,tile_grid_d
    type(tile_coord_type), dimension(:), pointer :: tile_coord_g => null()
    type(tile_coord_type), dimension(:), pointer :: tile_coord_d => null()
    integer, dimension(:), pointer     :: f2g => null()
    integer, dimension(:), pointer     :: d2g => null()
    integer, dimension(:), pointer     :: d2f => null()
    integer :: N_catg, N_catd,n1,n2,N_catf
    
    type(cat_param_type),  dimension(:), allocatable :: cp
    real :: dzsf
    
    namelist / domain_inputs /                              &
         minlon, maxlon,minlat,maxlat,                     &
         exclude_file,include_file
    
    inquire(file=trim(orig_tile),exist=file_exist)
    if( .not. file_exist) stop ("original tile file not exist")
    
    inquire(file=trim(domain_def),exist=d_exist)
    if( .not. d_exist) then
       print*,"no domain definition file"
    endif
    
    inquire(file=trim(catch_def_file),exist=c_exist)
    if( .not. c_exist) then
       print*,"no catchment definition file:" , catch_def_file
    endif
    
    
    if(d_exist) then
       open (10, file=trim(domain_def), delim='apostrophe', action='read', status='old')
       read (10, nml= domain_inputs)
       close(10)
    else
       minlon = -180.
       maxlon = 180.
       minlat = -90.
       maxlat = 90.
       exclude_file = ' '
       include_file = ' '
    endif
    
    call LDAS_read_til_file(orig_tile,catch_def_file,tile_grid_g,tile_coord_g,f2g)
    
    N_catg=size(tile_coord_g)
    
    ! include and exclude files are absolute
    
    call domain_setup(                                               &
         N_catg, tile_coord_g,                                        &
         tile_grid_g,                                                 &
         ' ', exclude_file, ' ', include_file,                         &
         trim(out_path), 'exp_domain ', trim(exp_id),                 &
         minlon, minlat, maxlon, maxlat,                              &
         N_catd, d2g, tile_coord_d,                                   &
         tile_grid_d )
    
    allocate(cp(N_catd))
    
    read(SURFLAY,*) dzsf
    print*, "SURFLAY: ", dzsf
    n1 = index(catch_def_file,'/clsm/')
    bcs_path(1:n1-1) = catch_def_file(1:n1-1)
    call read_cat_param( N_catg, N_catd, d2g, tile_coord_d, dzsf, bcs_path(1:n1-1), bcs_path(1:n1-1),bcs_path(1:n1-1),  &
         cp )
    call write_cat_param(cp,N_catd)
    
    allocate(d2f(N_catd))
    d2f = 0
    N_catf = size(f2g)
    if( N_catf /= N_catg) then
       n = 1
       do n1 = 1,N_catd
          do n2 = n, N_catf
             if (d2g(n1) == f2g(n2)) then
                d2f(n1) = n2
                n = n2+1
                exit
             endif
          enddo
       enddo
       if(any(d2f == 0)) stop " Domain includes those excluded tiles"
       print*," f2g now is d2f "
    else
       d2f = d2g
    endif
    open(40,file='f2g.txt',form='formatted',action='write')
    write(40,*)N_catf
    write(40,*)N_catd
    do n=1,N_catd
       write(40,*)d2f(n)
    enddo
    do n=1,N_catd
       write(40,*)d2g(n)
    enddo
    close(40)
    if (associated(f2g)) deallocate(f2g)
    if (associated(d2g)) deallocate(d2g)
    if (associated(d2f)) deallocate(d2f)
    
  contains

    ! ********************************************************************
    
    logical function is_in_list(N_list, list, this_one)
      
      ! checks whether "this_one" is element of list
      
      ! reichle, 2 May 2003
      
      implicit none
      
      integer :: N_list, this_one
      integer, dimension(N_list) :: list
      
      integer :: n
      
      ! ------------------------------------
      
      is_in_list = .false.
      
      do n=1,N_list
         
         if (list(n)==this_one) then
            is_in_list = .true.
            exit
         end if
      end do
      
    end function is_in_list
    
    ! ******************************************************************
    
    logical function is_in_domain(                                               &
         this_cat_exclude, this_cat_include, this_cat_in_box )
      
      ! determine whether catchment is in domain
      !
      ! The domain is set up using (if present) an "ExcludeList" of catchments 
      ! to be excluded, an "IncludeList" (if present) of catchments to be included,
      ! and the bounding box of a rectangular "zoomed" area (as specified
      ! in the "exeinp" file used in ldas_setup). 
      !
      ! order of precedence:
      !  1. exclude catchments on ExcludeList
      !  2. include catchments on IncludeList or catchments within rectangular domain
      ! (i.e., catchments in ExcludeList are *always* excluded)    
      !
      ! reichle, 7 May 2003
      ! reichle, 9 May 2005 - redesign (no more continents)
      !
      ! ----------------------------------------------------------------
      
      implicit none
      
      logical :: this_cat_include, this_cat_exclude, this_cat_in_box
      
      is_in_domain = .false.
      
      ! if catchment is NOT in ExcludeList 
      
      if (.not. this_cat_exclude)  then     
         
         ! if catchment is within bounding box OR in IncludeList
         
         if ((this_cat_in_box) .or. (this_cat_include)) then
            
            is_in_domain = .true.
            
         end if
      end if
      
    end function is_in_domain
    
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
      character(512) :: fname
      
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
         
         tile_grid_d = get_minExtent_grid( N_cat_domain, tile_coord%i_indg, tile_coord%j_indg, &
              tile_coord%min_lon, tile_coord%min_lat, tile_coord%max_lon, tile_coord%max_lat,  &
              tile_grid_g) 
         
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

      print*, "Done with " // trim(Iam)

    end subroutine domain_setup
    
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
    
    ! ****************************************************************
    
    integer function word_count( mystring )
      
      ! count number of words in "mystring" (delimited by space)
      !
      ! - reichle, 31 Mar 2015
      
      implicit none
      
      character(len=*) :: mystring
      
      integer :: N_words, N_string, ii
      
      logical :: current_is_space, next_is_space
      
      N_words = 0
      
      current_is_space = (mystring(1:1)==' ')
      
      if (.not. current_is_space)  N_words = N_words + 1
      
      do ii=2,len(mystring)
         
         next_is_space = (mystring(ii:ii)==' ')
         
         if (current_is_space .and. .not. next_is_space)  N_words = N_words + 1
         
         current_is_space = next_is_space
         
      end do
      
      word_count = N_words
      
    end function word_count
    
    ! ***********************************************************************
  
    integer function open_land_param_file( unitnumber, formatted_file, is_big_endian, &
         N_search_dir, fname, pathname, search_dir, ignore_stop )
      
      ! reichle, 13 Dec 2010
      ! reichle, 21 Oct 2011 - added optional output "istat" 
      ! reichle, 11 Dec 2013 - moved from "clsm_ensdrv_drv_routines.F90"
      !                         and converted to function
      
      ! try reading land or mwRTM parameter files from various sub-dirs for
      ! compatibility with old and new parameter directory structures
      
      ! fname      = file name (without path) of parameter file
      ! pathname   = path to parameter file
      ! search_dir = vector (length N_search_dir) of subdirectories to search 
      !               for file fname
      
      ! ignore_stop = optional input, if present and .true., skip call to "stop_it()" 
      
      implicit none
      
      integer                                  :: unitnumber, N_search_dir
      
      logical                                  :: formatted_file
      
      logical                                  :: is_big_endian
      
      character(*)                           :: fname
      
      character(*)                           :: pathname
      
      character(*), dimension(:)  :: search_dir
      
      logical,        optional                 :: ignore_stop
      
      ! local variables
      
      character(len=512) :: filename
      
      integer :: i, istat
      
      logical :: ignore_stop_tmp
      
      character(len=*), parameter :: Iam = 'open_land_param_file'
      character(len=400) :: err_msg
      
      ! ------------------------------------------------------------------
      !
      ! try opening file

      do i=1,N_search_dir
         
         filename = trim(pathname) // '/' // trim(search_dir(i)) // '/' // trim(fname)
         
         if (formatted_file) then
            
            open(unitnumber, file=filename, form='formatted',           &
                 action='read', status='old', iostat=istat)
            
         else
            
            if (is_big_endian) then
               
               open(unitnumber, file=filename, form='unformatted',      &
                    convert='big_endian',                               &
                    action='read', status='old', iostat=istat)
               
            else
               
               open(unitnumber, file=filename, form='unformatted',      &
                    convert='little_endian',                            &
                    action='read', status='old', iostat=istat)
               
            end if
            
         end if
         
         if (istat==0) exit  ! exit loop when first successful
         
      end do
      
      ! report back opened filename or stop (unless requested otherwise)
      
      if (istat==0) then
         
         if (logit) write (logunit,'(400A)') 'Reading from: ' // trim(filename)     
         
      else
         
         if (logit) then
            
            write (logunit,*) 'Cannot find file ', trim(fname), ' in: '
            
            do i=1,N_search_dir
               write (logunit,*) trim(pathname) // '/' // trim(search_dir(i))
            end do
            
         end if
         
         ! figure out whether to stop 
         
         ignore_stop_tmp = .false.  ! default: stop if file not opened successfully
         
         if (present(ignore_stop))  ignore_stop_tmp = ignore_stop
         
         if (.not. ignore_stop_tmp)  then
            err_msg = 'ERROR opening file'
            call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
         end if
         
      end if
      
      if (logit) write (logunit,*) 
      
      open_land_param_file = istat
      
    end function open_land_param_file
    
    ! *****************************************************************************************
    
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
      ! reichle, 14 Jul 2020 - work around for new "peat fraction" column in Icarus-NLv4 (ignore for now)
      ! 
      ! -------------------------------------------------------------------
      
      implicit none
      
      integer,                                  intent(in)  :: N_catg, N_catf
      
      type(tile_coord_type), dimension(:),      pointer     :: tile_coord_f ! intent(in)
      
      real,                                     intent(in)  :: dzsf
      
      integer,               dimension(N_catf), intent(in)  :: f2g
      
      character(*),                             intent(in)  :: veg_path
      character(*),                             intent(in)  :: soil_path
      character(*),                             intent(in)  :: top_path
      
      type(cat_param_type),  dimension(N_catf), intent(out) :: cp
      
      ! local variables
      
      integer, parameter :: N_search_dir_max =  5
      integer, parameter :: N_col_real_max   = 18  ! "v15" soil_param.dat had 22 columns (incl. first 4 columns with integers)
      
      character( 80) :: fname
      character(999) :: tmpstr999
      
      character(100), dimension(N_search_dir_max) :: search_dir
      
      integer :: n, k, m, dummy_int, dummy_int2, istat, N_search_dir, N_col
      
      integer, dimension(N_catg)                  :: tmpint, tmpint2, tmptileid
      
      real,    dimension(N_catg,N_col_real_max)   :: tmpreal
      
      real    :: dummy_real, dummy_real2, z_in_m, term1, term2
      
      logical :: dummy_logical
      
      character(len=*), parameter                 :: Iam = 'read_cat_param'
      character(len=400)                          :: err_msg
      
      real,    dimension(NTYPS)                   :: VGZ2
      
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
      
      tmpreal   = nodata_generic
      
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
         
      case (19,20)
         
         ! starting with "v16" (De Lannoy et al., 2014, doi:10.1002/2014MS000330),
         !  soil_param.dat has 19 columns
         
         ! "Icarus-NLv4" has 20 columns (new, last column is peat fraction, ignore for now)
         
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
      
      tmpreal   = nodata_generic
      
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
      
      tmpreal   = nodata_generic
      
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
      
      tmpreal   = nodata_generic
      
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
      
      tmpreal   = nodata_generic
      
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
    
    ! **********************************************************************************
  
    subroutine write_cat_param(cat_param, N_catd)
      type(cat_param_type), intent(in) :: cat_param(:)
      integer,intent(in) :: N_catd
      character(len=512):: fname
      type(date_time_type) :: start_time
      
      integer :: k,n
      
      read(ymdhm(1:4),*) start_time%year  ! 4-digit year
      read(ymdhm(5:6),*) start_time%month ! month in year
      read(ymdhm(7:8),*) start_time%day   ! day in month
      read(ymdhm(9:10),*) start_time%hour  ! hour of day
      read(ymdhm(11:12),*) start_time%min  
      start_time%sec = 0
      start_time%pentad  = -9999        ! pentad of year
      start_time%dofyr   = -9999   
      
      fname = get_io_filename(trim(out_path), trim(exp_id),'ldas_catparam', date_time=start_time, &
           dir_name='rc_out', file_ext='.bin')
      
      open(10, file=trim(fname), form='unformatted', status='unknown', action='write')
      
      print*, 'Writing catparam file : ' // trim(fname)
      
      write (10) (cat_param(n)%dpth,       n=1,N_catd)
      
      write (10) (cat_param(n)%dzsf,       n=1,N_catd)
      write (10) (cat_param(n)%dzrz,       n=1,N_catd)
      write (10) (cat_param(n)%dzpr,       n=1,N_catd)
      
      do k=1,N_gt
         write (10) (cat_param(n)%dzgt(k), n=1,N_catd)
      end do
      
      write (10) (cat_param(n)%poros,      n=1,N_catd)
      write (10) (cat_param(n)%cond,       n=1,N_catd)
      write (10) (cat_param(n)%psis,       n=1,N_catd)
      write (10) (cat_param(n)%bee,        n=1,N_catd)
      
      write (10) (cat_param(n)%wpwet,      n=1,N_catd)
      
      write (10) (cat_param(n)%gnu,        n=1,N_catd)
      
      write (10) (cat_param(n)%vgwmax,     n=1,N_catd)
      
      write (10) (real(cat_param(n)%vegcls),     n=1,N_catd)
      write (10) (real(cat_param(n)%soilcls30),  n=1,N_catd)
      write (10) (real(cat_param(n)%soilcls100), n=1,N_catd)
      
      write (10) (cat_param(n)%bf1,        n=1,N_catd)
      write (10) (cat_param(n)%bf2,        n=1,N_catd)
      write (10) (cat_param(n)%bf3,        n=1,N_catd)
      write (10) (cat_param(n)%cdcr1,      n=1,N_catd)
      write (10) (cat_param(n)%cdcr2,      n=1,N_catd)
      write (10) (cat_param(n)%ars1,       n=1,N_catd)
      write (10) (cat_param(n)%ars2,       n=1,N_catd)
      write (10) (cat_param(n)%ars3,       n=1,N_catd)
      write (10) (cat_param(n)%ara1,       n=1,N_catd)
      write (10) (cat_param(n)%ara2,       n=1,N_catd)
      write (10) (cat_param(n)%ara3,       n=1,N_catd)
      write (10) (cat_param(n)%ara4,       n=1,N_catd)
      write (10) (cat_param(n)%arw1,       n=1,N_catd)
      write (10) (cat_param(n)%arw2,       n=1,N_catd)
      write (10) (cat_param(n)%arw3,       n=1,N_catd)
      write (10) (cat_param(n)%arw4,       n=1,N_catd)
      write (10) (cat_param(n)%tsa1,       n=1,N_catd)
      write (10) (cat_param(n)%tsa2,       n=1,N_catd)
      write (10) (cat_param(n)%tsb1,       n=1,N_catd)
      write (10) (cat_param(n)%tsb2,       n=1,N_catd)
      write (10) (cat_param(n)%atau,       n=1,N_catd)
      write (10) (cat_param(n)%btau,       n=1,N_catd)
      
      write (10) (cat_param(n)%gravel30,   n=1,N_catd)
      write (10) (cat_param(n)%orgC30  ,   n=1,N_catd)
      write (10) (cat_param(n)%orgC    ,   n=1,N_catd)
      write (10) (cat_param(n)%sand30  ,   n=1,N_catd)
      write (10) (cat_param(n)%clay30  ,   n=1,N_catd)
      write (10) (cat_param(n)%sand    ,   n=1,N_catd)
      write (10) (cat_param(n)%clay    ,   n=1,N_catd)
      write (10) (cat_param(n)%wpwet30 ,   n=1,N_catd)
      write (10) (cat_param(n)%poros30 ,   n=1,N_catd)
      
      write (10) (cat_param(n)%veghght ,   n=1,N_catd)
      
      close (10,status='keep')
      
    end subroutine write_cat_param
    
  end subroutine createf2g
  
  ! ********************************************************************
  
  subroutine readsize(N_catg,N_catf)
    
    implicit none
    integer,intent(out) :: N_catg
    integer,intent(out) :: N_catf
    
    logical :: file_exist
    
    inquire(file=trim('f2g.txt'),exist=file_exist)
    if(file_exist) then
       open(40,file='f2g.txt',form='formatted',action='read',status='old')
       read(40,*)N_catg
       read(40,*)N_catf
       close(40)
    else
       print*, " wrong, no f2g.txt"
    endif
  end subroutine readsize
  
  ! ********************************************************************
  
  subroutine readf2g(N_catf,f2g)
    
    implicit none
    integer,intent(in) :: N_catf
    integer,dimension(N_catf),intent(inout) :: f2g
    
    integer :: N_catg
    logical :: file_exist
    integer :: local_size,n
    
    inquire(file=trim('f2g.txt'),exist=file_exist)
    if(file_exist) then
       open(40,file='f2g.txt',form='formatted',action='read',status='old')
       read(40,*)N_catg
       read(40,*)local_size
       
       if(local_size /= N_catf) print*, "wrong f2g.txt"
       
       if(N_catg == N_catf) then
          close(40)
          return
       endif
       
       do n=1,N_catf
          read(40,*)f2g(n)
       enddo
       close(40)
       ! call MAPL_sort(this%f2g)
    else
       print*, " wrong, no f2g.txt"
    endif
    
  end subroutine readf2g
  
  ! ********************************************************************
  
  subroutine createLocalTilefile(orig_tile,new_tile)
    
    implicit none
    character(*), intent(in) :: orig_tile
    character(*), intent(in) :: new_tile
    
    character(len=256) :: line
    character(len=3)   :: MAPL_Land_STRING
    character(len=4)   :: MAPL_Land_ExcludeFromDomain_STRING
    character(len=400) :: err_msg
    
    logical :: file_exist
    
    integer, dimension(:),allocatable :: f2g 
    integer :: N_catg, N_catf,n,stat, ty
    integer :: N_tile,N_grid,g_id
    
    character(len=*), parameter :: Iam = 'createLocalTilefile'
    
    ! string handling below relies on MAPL_Land and MAPL_Land_ExcludeFromDomain
    !  falling into a certain range
    
    ! verify that MAPL_Land has three digits
    
    if (MAPL_Land<100 .or. MAPL_Land>999) then   
       err_msg = 'string handling implemented only for 100<=MAPL_Land<=999'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    ! verify that MAPL_Land_ExcludeFromDomain has four digits
    
    if (MAPL_Land_ExcludeFromDomain<1000 .or. MAPL_Land_ExcludeFromDomain>9999) then   
       err_msg = 'string handling implemented only for 1000<=MAPL_Land_ExcludeFromDomain<=9999'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    ! convert integers to appropriate-length strings
    
    write (MAPL_Land_STRING,                  '(i3)') MAPL_Land   
    write (MAPL_Land_ExcludeFromDomain_STRING,'(i4)') MAPL_Land_ExcludeFromDomain
    
    inquire(file=trim(orig_tile),exist=file_exist)
    if( .not. file_exist) stop ("original tile file does not exist")
    
    ! Set default local tile file name
    call readsize(N_catg,N_catf)
    if(N_catg == N_catf) then
       print*, "It is global domain..."
       return
    endif
    allocate(f2g(N_catf))
    call readf2g(N_catf,f2g)   
    
    open(40,file=trim(orig_tile),action="read")
    open(50,file=trim(new_tile),action="write")
    
    ! copy the header back into the output tile file
    ! (also corrects bug in EASE *.til files that have "N_grid=1" in line 2 but
    !  still contain three additional lines for second grid definition)
    do n=1,5
       read(40,'(A)') line
       if(n==1) then
          read(line,*) N_tile
       endif
       if(n==2) then
          read(line,*) N_grid
       endif
       write(50,'(A)') trim(line)
    enddo
    if (N_grid==2) then
       do n=1,3
          read(40,'(A)') line
          write(50,'(A)') trim(line)
       enddo
    endif
    
    g_id = 0
    do while(.true.)
       ! read one line of *.til file
       read(40,'(A)',IOSTAT=stat) line
       if(IS_IOSTAT_END(stat)) exit
       ! extract first "integer" in "line" and put into "ty"
       read(line,*) ty
       if( ty == MAPL_Land ) then
          ! find index where MAPL_Land ("100") starts in "line"
          n=index(line,MAPL_Land_STRING)
          ! make sure that a space is available in front of MAPL_Land ("100")
          if (n<=1) then   
             err_msg = 'string handling requires at least one blank space in first column of *.til file'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          ! here g_id is (consecutive) id of the global *land* tiles
          g_id=g_id+1
          if(.not. any( f2g(:) == g_id)) then
             ! if tile is not in local domain, replace " 100" in "line" with "1100"
             line(n-1:n+2)=MAPL_Land_ExcludeFromDomain_STRING
          endif
       endif
       ! write "line" into the output tile file
       write(50,'(A)') trim(line)
    enddo
    close(40)
    close(50)
    
  end subroutine createLocalTilefile
  
  ! ********************************************************************
  
  subroutine createLocalBC(orig_BC, new_BC)
    
    implicit none
    character(*),intent(in) :: orig_BC
    character(*),intent(in) :: new_BC
    
    real,dimension(14) :: tmprealvec14
    real,allocatable ::   tmpvec(:)  
    integer :: istat, N_catg,N_catf
    integer,dimension(:),allocatable :: f2g
    
    call readsize(N_catg,N_catf)
    if(N_catg==N_catf) return
    allocate(f2g(N_catf))
    call readf2g(N_catf,f2g)  
    
    allocate(tmpvec(N_catg))
    open(10,file=trim(orig_BC),form='unformatted',action='read',status='old',iostat=istat)
    open(20,file=trim(new_BC),form='unformatted',action='write')
    
    do while(.true.)
       read(10,iostat=istat) tmprealvec14
       if(IS_IOSTAT_END(istat)) exit
       read(10) tmpvec
       write(20) tmprealvec14
       write(20) tmpvec(f2g)
    enddo
    close(10)
    close(20)
    deallocate(tmpvec)
  end subroutine createLocalBC

  ! ********************************************************************
  
  subroutine createLocalCatchRestart(orig_catch, new_catch)
    
    implicit none
    character(*),intent(in):: orig_catch
    character(*),intent(in):: new_catch
    integer,parameter :: subtile=4
    integer :: istat, filetype, rc,i, j, ndims
    real,allocatable :: tmp1(:)
    real,allocatable :: tmp2(:,:)
    type(Netcdf4_FileFormatter) :: InFmt,OutFmt
    type(FileMetadata)        :: OutCfg
    type(FileMetadata)        :: InCfg
    integer                   :: dim1,dim2
    type(StringVariableMap), pointer :: variables
    type(Variable), pointer :: var
    type(StringVariableMapIterator) :: var_iter
    type(StringVector), pointer :: var_dimensions
    character(len=:), pointer :: vname,dname
    integer ::n, N_catg,N_catf
    integer,dimension(:),allocatable :: f2g
    
    call readsize(N_catg,N_catf)
    if(N_catg == N_catf) return
    allocate(f2g(N_catf))
    call readf2g(N_catf,f2g) 
    
    allocate(tmp1(N_catg))
    allocate(tmp2(N_catg,subtile))
    
    ! check file type
    
    call MAPL_NCIOGetFileType(orig_catch, filetype,rc=rc)
    
    if (filetype /= 0) then
       
       print*, "Catchment restart is binary"
       
       ! binary 
       
       open(10,file=trim(orig_catch),form='unformatted',action='read',status='old',iostat=istat)
       open(20,file=trim(new_catch),form='unformatted',action='write')
       
       do n=1,30
          read(10) tmp1
          write(20) tmp1(f2g)
       enddo
       
       do n=1,2
          read(10) tmp2
          write(20) tmp2(f2g,:)
       enddo
       
       do n=1,20
          read(10) tmp1
          write(20) tmp1(f2g)
       enddo
       ! note : the offline restart does not have the last five variables
       do n=1,4
          read(10,iostat=istat) tmp2
          if(.not. IS_IOSTAT_END(istat)) write(20) tmp2(f2g,:)
       enddo
       ! 57 WW
       read(10,iostat=istat) tmp2
       if(.not. IS_IOSTAT_END(istat)) write(20) tmp2(f2g,:)
       
       close(10)
       close(20)
    else
       
       ! filetype = 0 : nc4 output file will also be nc4
       
       call InFmt%open(trim(orig_catch), pFIO_READ,rc=rc)
       InCfg  = InFmt%read(rc=rc)
       OutCfg = InCfg
       
       call OutCfg%modify_dimension('tile', size(f2g), rc=rc)
       
       call OutFmt%create(trim(new_catch),rc=rc)
       call OutFmt%write(OutCfg,rc=rc)
       
       variables => InCfg%get_variables()
       var_iter = variables%begin()
       do while (var_iter /= variables%end())
          
          vname => var_iter%key()
          var => var_iter%value()
          var_dimensions => var%get_dimensions()
          
          ndims = var_dimensions%size()
          
          if (trim(vname) =='time') then
             call var_iter%next()
             cycle
          endif
          
          if (ndims == 1) then
             call MAPL_VarRead (InFmt,vname,tmp1)
             call MAPL_VarWrite(OutFmt,vname,tmp1(f2g))
          else if (ndims == 2) then
             
             dname => var%get_ith_dimension(2)
             dim1=InCfg%get_dimension(dname)
             do j=1,dim1
                call MAPL_VarRead ( InFmt,vname,tmp1 ,offset1=j)
                call MAPL_VarWrite(OutFmt,vname,tmp1(f2g),offset1=j)
             enddo
             
          else if (ndims == 3) then
             
             dname => var%get_ith_dimension(2)
             dim1=InCfg%get_dimension(dname)
             dname => var%get_ith_dimension(3)
             dim2=InCfg%get_dimension(dname)
             do i=1,dim2
                do j=1,dim1
                   call MAPL_VarRead ( InFmt,vname,tmp1 ,offset1=j,offset2=i)
                   call MAPL_VarWrite(OutFmt,vname,tmp1(f2g) ,offset1=j,offset2=i)
                enddo
             enddo
             
          end if
          call var_iter%next()
       enddo
       call inFmt%close(rc=rc)
       call OutFmt%close(rc=rc)
    end if ! file type nc4
    print*, "done create local catchment restart"
  end subroutine createLocalCatchRestart

  ! ********************************************************************
  
  subroutine createLocalmwRTMRestart(orig_mwrtm, new_mwrtm)
    
    implicit none
    character(*),intent(in):: orig_mwrtm
    character(*),intent(in):: new_mwrtm
    integer,parameter :: subtile=4
    integer :: rc
    real,allocatable :: tmp1(:)
    type(Netcdf4_FileFormatter) :: InFmt,OutFmt
    type(FileMetadata)        :: OutCfg
    type(FileMetadata)        :: InCfg
    
    type(StringVariableMap), pointer :: variables
    type(StringVariableMapIterator) :: var_iter
    character(len=:), pointer :: vname
    integer :: N_catg,N_catf
    integer,dimension(:),allocatable :: f2g
    
    call readsize(N_catg,N_catf)
    if(N_catg == N_catf) return
    allocate(f2g(N_catf))
    call readf2g(N_catf,f2g) 
    
    allocate(tmp1(N_catg))
    
    ! nc4 in and out file will also be nc4
    call InFmt%open(trim(orig_mwrtm), pFIO_READ,rc=rc)
    InCfg = InFmt%read(rc=rc)
    OutCfg = InCfg
    
    call OutCfg%modify_dimension('tile', size(f2g), rc=rc)
    
    call OutFmt%create(trim(new_mwrtm),rc=rc)
    call OutFmt%write(OutCfg,rc=rc)
    
    variables => InCfg%get_variables()
    var_iter = variables%begin()
    do while (var_iter /= variables%end())
       vname => var_iter%key()
       call MAPL_VarRead (InFmt,vname,tmp1)
       call MAPL_VarWrite(OutFmt,vname,tmp1(f2g))
       call var_iter%next()
    enddo
    
    call inFmt%close(rc=rc)
    call OutFmt%close(rc=rc)
    
    deallocate(f2g,tmp1)
    
  end subroutine createLocalmwRTMRestart
  
  ! ********************************************************************
  
  subroutine createLocalVegRestart(orig_veg, new_veg)
    
    implicit none
    character(*),intent(in):: orig_veg
    character(*),intent(in):: new_veg
    integer :: istat
    real,allocatable :: rity(:)
    real,allocatable :: z2(:)
    real,allocatable :: ascatz0(:)
    real,allocatable :: tmp(:)
    
    integer :: N_catg,N_catf
    integer,dimension(:),allocatable :: f2g
    integer :: filetype
    type(Netcdf4_FileFormatter) :: InFmt,OutFmt
    type(FileMetadata)        :: OutCfg
    type(FileMetadata)        :: InCfg
    
    type(StringVariableMap), pointer :: variables
    type(StringVariableMapIterator) :: var_iter
    character(len=:), pointer :: vname
    integer :: rc
    
    call readsize(N_catg,N_catf)
    if(N_catg == N_catf) return
    allocate(f2g(N_catf))
    call readf2g(N_catf,f2g)  
    
    allocate(rity(N_catg))
    allocate(z2(N_catg))
    allocate(ascatz0(N_catg))
    
    call MAPL_NCIOGetFileType(orig_veg, filetype,rc=rc)
    
    if (filetype /=0) then
       open(10,file=trim(orig_veg),form='unformatted',action='read',status='old',iostat=istat)
       open(20,file=trim(new_veg),form='unformatted',action='write')
       read(10) rity 
       read(10) z2 
       read(10) ascatz0 
       write(20) rity(f2g)
       write(20) z2(f2g) 
       write(20) ascatz0(f2g) 
       
       close(10)
       close(20)
    else
       ! nc4 in and out file will also be nc4
       call InFmt%open(trim(orig_veg), pFIO_READ,rc=rc)
       InCfg = InFmt%read(rc=rc)
       OutCfg = InCfg
       
       call OutCfg%modify_dimension('tile', size(f2g), rc=rc)
       
       call OutFmt%create(trim(new_veg),rc=rc)
       call OutFmt%write(OutCfg,rc=rc)
       
       variables => InCfg%get_variables()
       var_iter = variables%begin()
       allocate(tmp(N_catg))
       do while (var_iter /= variables%end())
          vname => var_iter%key()
          call MAPL_VarRead (InFmt,vname,tmp)
          call MAPL_VarWrite(OutFmt,vname,tmp(f2g))
          call var_iter%next()
       enddo
       
       call inFmt%close(rc=rc)
       call OutFmt%close(rc=rc)
       deallocate(tmp)
    endif
    deallocate(f2g)
    
  end subroutine createLocalVegRestart
  
  ! ********************************************************************
  
  subroutine correctEase(orig_ease,new_ease)
    
    ! This subroutine corrects for a bug that is present in some EASE *.til files
    ! through at least Icarus-NLv4.
    ! Affected files state "N_grid=1" in line 2 of the header, but the header still includes
    ! three additional lines for a second grid, which throws off the canonical *.til reader
    ! (subroutine LDAS_read_til_file()).
    !
    ! This subroutine creates a second, corrected version of the *.til file that can be
    ! read with the canonical reader during ldas_setup.
    !
    ! - reichle,  2 Aug 2020
    
    implicit none
    character(*),intent(in) :: orig_ease
    character(*),intent(in) :: new_ease
    logical :: file_exist,is_oldEASE
    integer :: i, N_tile, N_grid
    character(len=256) :: tmpline
    
    inquire(file=trim(orig_ease),exist=file_exist)
    if( .not. file_exist) stop (" no ease_tile_file")
    
    open(55,file=trim(orig_ease),action='read')
    read(55,*) N_tile
    read(55,*) N_grid
    read(55,*)
    read(55,*)
    read(55,*)
    read(55,'(A)') tmpline
    close(55)
    
    is_oldEASE= .false.
    if(N_grid==1 .and. index(tmpline,'OCEAN')/=0) is_oldEASE=.true.
    
    if( is_oldEASE) then
       open(55,file=trim(orig_ease),action='read')
       open(56,file=trim(new_ease),action='write')
       do i =1,5
          read(55,'(A)')tmpline
          write(56,'(A)')trim(tmpline)
       enddo
       read(55,*)
       read(55,*)
       read(55,*)
       do i=1,N_tile
          read(55,'(A)')tmpline
          write(56,'(A)')trim(tmpline)
       enddo
       close(56)
       close(55)
    end if
  end subroutine correctEase
  
  ! ********************************************************************
  !
  ! subroutine to optimize the domain setup (domain decomposition and processor layout)
  !
  ! The domain is cut into N_proc stripes, where
  !   - N_proc is the number of processors used for the model simulation (i.e., excl. OSERVER tasks)
  !   - a stripe cuts across the entire north-south extent of the domain (lat-lon, EASE) or
  !      across an entire face (cube-sphere) -- see below
  !   - each stripe must be at least two grid cells thick
  !   - each stripe should contain roughly the same number of land *tiles*
  !
  ! For lat-lon and EASE grid tile spaces:
  !   - cut into N_proc stripes of size N_lat-by-IMS(k),  k=1:N_proc, IMS(k)>=2,
  !     where IMS(k)=no. of grid cells in longitude direction ("thickness" of stripe k) is
  !     written into file IMS.rc
  ! For cube-sphere grid tile spaces:
  !   - cut into N_proc stripes of size N_face-by-JMS(k), k=1:N_proc, JMS(k)>=2,
  !     where JMS(k)=no. of grid cells ("thickness" of stripe k) is written
  !     into file JMS.rc
  !
  !         cube-sphere           lat-lon/EASE
  ! ---------------------------------------------------
  ! NX:     1                     N_proc
  ! NY:     N_proc                1
  !         JMS.rc                IMS.rc
  
  subroutine optimize_latlon(fname_tilefile, N_proc_string)
    
    implicit none
    
    character(*), intent(in) :: fname_tilefile  ! file name (with path) of tile file (*.til)
    character(*), intent(in) :: N_proc_string   ! *string* w/ no. of processors (or tasks), excl. OSERVER tasks
        
    ! local variables
    integer :: N_proc
    integer :: N_tile,N_lon,N_lat,N_grid
    integer,allocatable :: landPosition(:)
    integer,allocatable :: IMS(:),JMS(:)
    integer,allocatable :: local_land(:)
    integer :: total_land
    integer :: n,typ,tmpint
    real ::  tmpreal
    integer :: avg_land,n0,local
    integer :: i,s,e,j,k,n1,n2
    logical :: file_exist
    character(len=256):: tmpLine
    character(len=128):: gridname
    real :: rate,rates(60),maxf(60)
    integer :: IMGLOB, JMGLOB
    integer :: face(6),face_land(6)
    logical :: forward

    ! -----------------------------
    
    read (N_proc_string,*) N_proc   ! input is string for historical reasons...
    
    ! get tile info

    inquire(file=trim(fname_tilefile),exist=file_exist)
    if( .not. file_exist) stop ( "tile file does not exist")
    
    open (10, file=trim(fname_tilefile), form='formatted', action='read')
    read (10,*) N_tile
    read (10,*) N_grid         ! some number (?)
    read (10,*) gridname       ! some string describing tile definition grid (?)
    read (10,*) N_lon
    read (10,*) n_lat
    
    if (index(gridname,"CF") /=0) then    ! cube-sphere tile space
       
       IMGLOB = N_lon                     ! e.g.,  180 for c180
       JMGLOB = N_lat                     ! e.g., 1080 for c180  (6*180=1080)
       if(JMGLOB/6 /= IMGLOB) stop " wrong im, jm"
       
       allocate(landPosition(JMGLOB))
       landPosition = 0
       total_land   = 0
       
       if(N_grid==2) then
          read (10,*)          ! some string describing ocean grid                   (?)
          read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
          read (10,*)
       endif
       
       do n = 1,N_tile
          read (10,*)  &
               typ,           &   !  1
               tmpreal,       &   !  2  *
               tmpreal,       &   !  3
               tmpreal,       &   !  4
               i ,            &   !  5
               j                  !  6
          !tmpreal,       &   !  7
          !tmpint,        &   !  8
          !tmpreal,       &   !  9  *
          !tmpint,        &   ! 10
          !tmpreal,       &   ! 11
          !tmpint             ! 12  * (previously "tile_id")
          if(typ==MAPL_Land) then
             total_land=total_land+1
             landPosition(j) = landPosition(j)+1
          endif
          
          ! assume all land tiles are at the beginning
          ! UNSAFE ASSUMPTION! - reichle, 2 Aug 2020
          
          if (typ/=MAPL_Land .and. typ/=MAPL_Land_ExcludeFromDomain) then   ! exit if not land
             
             if (logit) then
                write (logunit,*) 'WARNING: Encountered first non-land tile in *.til file.'
                write (logunit,*) '         Stop reading *.til file under the assumption that'
                write (logunit,*) '           land tiles are first in *.til file.'
                write (logunit,*) '         This is NOT a safe assumption beyond Icarus-NLv[x] tile spaces!!'
             end if
             
             exit ! assuming land comes first in the til file
             
          end if
          
       enddo
       close(10)
       
       if(mod(N_proc,6) /=0) then
          print*,"WARNING: ntasks should be adjusted to multiple of 6 for cubed-sphere grid :",N_proc
          N_proc = N_proc-mod(N_proc,6)
       endif
       
       print*, "total tiles", total_land
       
       if(sum(landPosition) /= total_land) print*, "wrong counting of land"
       
       do k=1,6
          n1 = (k-1)*IMGLOB+1
          n2 = k*IMGLOB
          face_land(k) = sum(landPosition(n1:n2)) 
          face(k) = nint(1.0*face_land(k)/total_land * N_proc)
          if ( face(k) == 0) face(k) = 1
       enddo
       
       ! now make sure sum(face) == N_proc
       k=sum(face)-N_proc
       
       if (k < 0) then
          do i=1, -k
             n=minloc(face,DIM=1)
             face(n) = face(n)+1
          enddo
       else
          do i = 1,k
             n=maxloc(face,DIM=1)
             face(n) = face(n)-1
          enddo
       endif
       
       if (sum(face) /= N_proc) stop " wrong proc face"
       
       ! 2) each process should have average land tiles
       
       ALLOCATE(JMS(N_proc))
       allocate(local_land(N_Proc))
       JMS = 0
       local_land = 0
       
       local  = 0
       n0     = 0
       j = 0
       do k=1,6
          n1 = (k-1)*IMGLOB+1
          n2 = k*IMGLOB
          
          do i=1,60
             rates(i) = -0.3 + i*0.01
          enddo
          
          maxf=rms_cs(rates)
          i=minloc(maxf,DIM=1)
          rate = rates(i)
          avg_land = ceiling(1.0*face_land(k)/face(k))
          avg_land = avg_land - nint(rate*avg_land)
          
          tmpint = 0
          j = j+face(k) 
          forward = .true.
          do n = n1,n2
             tmpint=tmpint+landPosition(n)
             if((local+1) == j .and. n < n2) cycle
             if(n==n2) then
                local = local + 1
                local_land(local)=tmpint
                JMS(local)=n-n0
                tmpint=0
                n0=n
                cycle
             endif
             if(tmpint .ge. avg_land) then
                local = local + 1
                if (n-n0 == 1) forward =.true.
                if (forward) then
                   local_land(local)=tmpint
                   JMS(local)=n-n0
                   tmpint=0
                   n0=n
                   forward = .false.
                else
                   local_land(local)=tmpint - landPosition(n)
                   JMS(local)=n-1-n0
                   tmpint=landPosition(n)
                   n0=n-1
                   forward = .true.
                endif
            endif
          enddo
          local = j
       enddo
       if( sum(JMS) /= JMGLOB) then
          print*, sum(JMS), JMGLOB
          stop ("wrong cs-domain distribution in the first place")
       endif 
       ! adjust JMS.rc to make sure each processor has at least 2 grid cells in j dimension 
       ! (i.e., each proc's subdomain must include at least 2 latitude stripes;
       !  stripes of grid cells may or may not contain land tiles)
       j = 1
       do k = 1,6
          n1 = j
          n2 = j+face(k)-1
          do i = n1,n2
             if(JMS(i) == 0) then
                n = maxloc(JMS(n1:n2),DIM=1)
                JMS(i) = 1
                JMS(n+n1-1) = JMS(n+n1-1)-1
             endif
             if(JMS(i) == 1) then
                n = maxloc(JMS(n1:n2),DIM=1)
                JMS(i) = 2
                JMS(n+n1-1) = JMS(n+n1-1)-1
             endif
          enddo
          j=j+face(k)
       enddo
       
       print*,"land_distribute: ",local_land
       print*, "JMS.rc", JMS
       if( sum(JMS) /= JMGLOB) then
          print*, sum(JMS), JMGLOB
          stop ("wrong cs-domain distribution")
       endif
       tmpint = 0
       k = 0
       do n = 1, N_proc
          tmpint= tmpint+JMS(n)
          if( tmpint == IMGLOB) then
             k=k+1
             tmpint = 0
          endif
       enddo
       
       if( k /=6 ) stop ("one or more processes may accross the face")
       
       open(10,file="optimized_distribution",action='write')
       write(10,'(A)')    "GEOSldas.GRIDNAME:  " // trim(gridname)
       write(10,'(A)')    "GEOSldas.GRID_TYPE:  Cubed-Sphere"
       write(10,'(A)')    "GEOSldas.NF:  6"
       write(10,'(A,I6)') "GEOSldas.IM_WORLD: ", IMGLOB
       write(10,'(A)')    "GEOSldas.LM:   1"
       write(10,'(A,I5)') "NY: ",N_proc
       write(10,'(A)')    "NX:   1"
       write(10,'(A)')    "GEOSldas.JMS_FILE:    JMS.rc"
       close(10)
       
       open(10,file="JMS.rc",action='write')
       write(10,'(I5,I5)') N_proc, maxval(face)
       do n=1,N_proc
          write(10,'(I8)') JMS(n)
       enddo
       close(10)
       
    else
       
       allocate(IMS(N_Proc))
       allocate(local_land(N_Proc))
       IMS=0
       local_land = 0
       
       ! NOTE:
       !  There is a bug in at least some EASE *.til files through at least Icarus-NLv4.
       !  Affected files state "N_grid=1" in line 2 of the header, but the header still includes
       !  three additional lines for a second grid.
       !
       !  The "else" block below corrects for this bug.
       !
       !  Elsewhere, LDAS pre-processing corrects for this bug through subroutine correctEase() in
       !  preprocess_LDAS.F90, which creates a second, corrected version of the *.til file during
       !  ldas_setup.
       !
       ! -reichle, 2 Aug 2020
       
       if(N_grid==2) then
          read (10,*)          ! some string describing ocean grid                   (?)
          read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
          read (10,*)   
          read(10,'(A)') tmpLine
       else
          read(10,'(A)') tmpLine
          if (index(tmpLine,"OCEAN") /=0) then
             read (10,*)          ! some string describing ocean grid                   (?)
             read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
             read (10,*)   
             read(10,'(A)') tmpLine
          endif
       endif
       
       if (index(gridname,'EASE') /=0) then
          s=0
          e=N_lon-1
       else
          s=1
          e=N_lon
       endif
       allocate(landPosition(s:e))
       
       landPosition=0
       total_land= 0
       
       ! 1) read through tile file, put the land tile into the N_lon of bucket
       
       read (tmpLine,*)    &
            typ,           &   !  1
            tmpreal,       &   !  2  *
            tmpreal,       &   !  3
            tmpreal,       &   !  4
            i                  !  5
       if(typ==MAPL_Land) then
          total_land=total_land+1
          landPosition(i) = landPosition(i)+1
       endif
       
       do n = 2,N_tile
          read (10,*)         &
               typ,           &   !  1
               tmpreal,       &   !  2  *
               tmpreal,       &   !  3
               tmpreal,       &   !  4
               i                  !  5
          !tmpint,        &   !  6
          !tmpreal,       &   !  7
          !tmpint,        &   !  8
          !tmpreal,       &   !  9  *
          !tmpint,        &   ! 10
          !tmpreal,       &   ! 11
          !tmpint             ! 12  * (previously "tile_id")
          if(typ==MAPL_Land) then
             total_land=total_land+1
             landPosition(i) = landPosition(i)+1
          endif
          
          ! assume all land tiles are at the beginning
          ! UNSAFE ASSUMPTION! - reichle, 2 Aug 2020
          
          if (typ/=MAPL_Land .and. typ/=MAPL_Land_ExcludeFromDomain) then   ! exit if not land
             
             if (logit) then
                write (logunit,*) 'WARNING: Encountered first non-land tile in *.til file.'
                write (logunit,*) '         Stop reading *.til file under the assumption that'
                write (logunit,*) '           land tiles are first in *.til file.'
                write (logunit,*) '         This is NOT a safe assumption beyond Icarus-NLv[x] tile spaces!!'
             end if
             
             exit ! assuming land comes first in the til file
             
          end if
          
       enddo
       
       close(10)
       
       if(sum(landPosition) /= total_land) print*, "wrong counting of land"
       
       do n=1,60
          rates(n) = -0.3 + (n-1)*0.01
       enddo
       
       maxf=rms(rates)
       n=minloc(maxf,DIM=1)
       rate = rates(n)
       
       ! 2) each process should have average land tiles
       
       avg_land = ceiling(1.0*total_land/N_proc)
       print*,"avg_land",avg_land
       
       ! rate is used to readjust the avg_land
       ! in case that the last processors don't have any land tiles,
       ! we can increase ther rates
       
       avg_land = avg_land - nint(rate*avg_land)
       print*,"re adjust the avg_land",avg_land
       tmpint = 0
       local = 1
       n0 = s-1
       forward = .true.
       do n=s,e
          tmpint=tmpint+landPosition(n)
          if(local == N_proc .and. n < e) cycle ! all lefteover goes to the last process
          if( n==e ) then
             local_land(local)=tmpint
             IMS(local)=n-n0
             exit
          endif

          if( tmpint .ge. avg_land ) then
             if (forward .or. n-n0 == 1 ) then
                local_land(local)=tmpint
                IMS(local)=n-n0
                tmpint=0
                n0=n
                forward = .false.
             else
                local_land(local) = tmpint - landPosition(n)
                IMS(local)=(n-1)-n0
                tmpint= landPosition(n)
                n0 = n-1
                forward = .true.
             endif
             local = local + 1
          endif
       enddo
       print*,"rms rate: ", rms(rate)
       
       print*,"land_distribute: ",local_land
       
       if( sum(local_land) /= total_land) stop ("wrong distribution")
       if( sum(IMS) /= N_lon) stop ("wrong domain distribution")

       ! redistribute IMS and try to make it >=2 (may be impossible for large N_Proc)
       do i = 1, N_proc
          if(IMS(i) == 0) then
            n = maxloc(IMS,DIM=1)
            IMS(i) = 1
            IMS(n) = IMS(n)-1
          endif
          if(IMS(i) == 1) then
            n = maxloc(IMS,DIM=1)
            IMS(i) = 2
            IMS(n) = IMS(n)-1
          endif
       enddo
       if( any(IMS <=1) ) stop ("Each processor must have at least 2 longitude stripes. Request fewer processors.")  

       open(10,file="optimized_distribution",action='write')
       write(10,'(A)')    "GEOSldas.GRID_TYPE:  LatLon"
       write(10,'(A)')    "GEOSldas.GRIDNAME:   "//trim(gridname)
       write(10,'(A)')    "GEOSldas.LM:         1"
       write(10,'(A)')    "GEOSldas.POLE:       PE"
       write(10,'(A)')    "GEOSldas.DATELINE:   DE"
       write(10,'(A,I6)') "GEOSldas.IM_WORLD:   ", N_lon
       write(10,'(A,I6)') "GEOSldas.JM_WORLD:   ", N_lat
       
       write(10,'(A,I5)') "NX:                  ",N_proc
       write(10,'(A)')    "NY:                  1"
       
       write(10,'(A)')    "GEOSldas.IMS_FILE:   IMS.rc"
       close(10)
       
       open(10,file="IMS.rc",action='write')
       write(10,'(I5)') N_proc
       do n=1,N_proc
          write(10,'(I8)') IMS(n)
       enddo
       close(10)
       
    endif
    
  contains 
    
    ! ***************************************************************************
    
    elemental function rms(rates) result (f)
      real :: f
      real,intent(in) :: rates
      integer :: tmpint,local
      integer :: n0,proc,n
      integer :: avg_land
      integer,allocatable :: local_land(:)
      logical :: forward
 
      allocate (local_land(N_proc))
      local_land = 0
      avg_land = ceiling(1.0*total_land/N_proc)
      avg_land = avg_land -nint(rates*avg_land)

      forward = .true.      
      tmpint = 0
      local = 1
      n0 = s-1
      do n=s,e
         tmpint=tmpint+landPosition(n)
         if(local == N_proc .and. n < e) cycle ! all lefteover goes to the last process
         if( n==e ) then
             local_land(local)=tmpint
             exit
          endif

          if( tmpint .ge. avg_land ) then
             if (forward .or. n-n0 == 1 ) then
                local_land(local)=tmpint
                tmpint=0
                n0=n
                forward = .false.
             else
                local_land(local) = tmpint - landPosition(n)
                tmpint= landPosition(n)
                n0 = n-1
                forward = .true.
             endif
             local = local + 1
          endif
      enddo
      f = 0.0
      do proc = 1, N_proc
         f =max(f,1.0*abs(local_land(proc)-avg_land))
      enddo
      deallocate(local_land)
    end function rms
    
    ! ***************************************************************************
    
    elemental function rms_cs(rates) result (f)
      real :: f
      real,intent(in) :: rates
      integer :: tmpint,local
      integer :: proc,n
      integer :: avg_land
      integer,allocatable :: local_land(:)
      integer :: n1,n2,n0
      logical :: forward

      allocate (local_land(face(k)))
      local_land = 0
      avg_land = ceiling(1.0*face_land(k)/face(k))
      avg_land = avg_land -nint(rates*avg_land)
      if (avg_land <=0) then
         f = face_land(k)
         return
      endif
      
      tmpint = 0
      local = 1
      
      n1 = (k-1)*IMGLOB+1
      n2 = k*IMGLOB
      tmpint = 0
      forward = .true.
      n0 = n1-1
      do n = n1,n2
         tmpint=tmpint+landPosition(n)
         if(local == face(k) .and. n < n2) cycle ! all lefteover goes to the last process
         if(n==n2) then
            local_land(local)=tmpint
            local = local + 1
            cycle
         endif
         if(tmpint .ge. avg_land) then
            if (n -n0 == 1) forward = .true. ! if only one step, should not got backward
            if (forward) then
               local_land(local)=tmpint
               tmpint=0
               n0 = n
               forward = .false.
            else
               local_land(local) = tmpint - landPosition(n)
               tmpint = landPosition(n)
               n0 = n-1
               forward = .true.
            endif
            local = local + 1
         endif
      enddo
      
      f = 0.0
      do proc = 1, face(k)
         ! punish for no land tiles
         f =max(f,1.0*abs(local_land(proc)-avg_land))
      enddo
      deallocate(local_land)
    end function rms_cs
    
  end subroutine optimize_latlon
  
  ! ********************************************************************
  
  subroutine convert_pert_rst(pfile_name,pfile_nc4,in_path,exp_id)
    
    implicit none
    character(*),intent(in) :: pfile_name
    character(*),intent(in) :: in_path
    character(*),intent(in) :: exp_id
    character(*),intent(in) :: pfile_nc4
    
    integer :: N_catf,N_lon,N_lat,N_lonf,N_latf
    integer :: N_force_pert,N_progn_pert
    integer,pointer :: f2g(:)
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord_f => null()
    
    type(grid_def_type) :: pert_grid_g
    type(grid_def_type) :: pert_grid_f 
    integer :: RC,istat
    integer,allocatable :: Pert_rseed(:)
    real,allocatable :: Force_pert_ntrmdt_f(:,:,:)
    real,allocatable :: Progn_pert_ntrmdt_f(:,:,:)
    
    call io_domain_files('r',in_path, trim(exp_id),N_catf,f2g,tile_coord_f,pert_grid_g,pert_grid_f,RC) 

    N_lon = pert_grid_g%N_lon
    N_lat = pert_grid_g%N_lat
    N_lonf= pert_grid_f%N_lon
    N_latf= pert_grid_f%N_lat

    call i_pert_ldas(RC)

    call o_pert_GEOSldas(rc)

  contains 
    
    ! ***************************************************************************
    
    subroutine i_pert_ldas(rc)
      integer,intent(inout),optional :: rc
      
      integer :: nrandseed_tmp
      type(grid_def_type) :: pert_grid_f_tmp 
      character(len=*), parameter :: Iam = 'io_pert_rstrt'
      integer :: k
      real, allocatable :: real_tmp(:)
      
      open(10, file=pfile_name, convert='big_endian',form='unformatted', status='old', &
           action='read', iostat=istat)
      
      ! one additional header line (as of 21 May 2010)!!!
      
      call io_grid_def_type( 'r', 10, pert_grid_f_tmp )
      
      read (10) nrandseed_tmp, N_force_pert, N_progn_pert
      
      ! check whether entries in file match passed arguments
      ! (check does *not* include *_pert_param!)
      
      if ( (nrandseed_tmp          /= NRANDSEED) ) then !          .or.               &
         !     (N_force_pert_tmp       /= N_force_pert)         .or.               &
         !     (N_progn_pert_tmp       /= N_progn_pert) ) then
         stop 'pert.rstrt file not compatible (1)'
      end if
      
      allocate(Pert_rseed(NRANDSEED))
      allocate(Force_pert_ntrmdt_f(N_lonf,N_latf, N_Force_pert)) 
      allocate(Progn_pert_ntrmdt_f(N_lonf,N_latf, N_Progn_pert))
      
      if ( index(pert_grid_f%gridtype,'LatLon')/=0  .or.         &
           index(pert_grid_f%gridtype,'LATLON')/=0  .or.         &
           index(pert_grid_f%gridtype,'latlon')/=0       ) then
         
         if ( (pert_grid_f_tmp%N_lon  /= pert_grid_f%N_lon)    .or.            &
              (pert_grid_f_tmp%N_lat  /= pert_grid_f%N_lat)    .or.            &
              (abs(pert_grid_f_tmp%ll_lon - pert_grid_f%ll_lon) > 1e-4) .or.   &
              (abs(pert_grid_f_tmp%ll_lat - pert_grid_f%ll_lat) > 1e-4) .or.   &
              (abs(pert_grid_f_tmp%dlon   - pert_grid_f%dlon)   > 1e-4) .or.   &
              (abs(pert_grid_f_tmp%dlat   - pert_grid_f%dlat)   > 1e-4)      ) then
            stop 'pert.rstrt file not compatible (2)'
         end if
         
      else
         
         if ( index(pert_grid_f_tmp%gridtype,pert_grid_f%gridtype)==0  .or.   &
              (pert_grid_f_tmp%N_lon  /= pert_grid_f%N_lon)            .or.   &
              (pert_grid_f_tmp%N_lat  /= pert_grid_f%N_lat)            .or.   &
              (pert_grid_f_tmp%i_offg /= pert_grid_f%i_offg)           .or.   &
              (pert_grid_f_tmp%j_offg /= pert_grid_f%j_offg)                ) then
            stop 'pert.rstrt file not compatible (3)'
         end if
         
      end if
      
      ! reading
      read (10) Pert_rseed(:)
      allocate(real_tmp(N_lonf*N_latf))
      do k=1,N_force_pert
         !read (10) ((Force_pert_ntrmdt_f(i,j,k), i=1,N_lonf),j=1,N_latf)
         read (10) real_tmp(:)
         Force_pert_ntrmdt_f(:,:,k) = reshape(real_tmp,[N_lonf, N_latf])
      end do
      
      do k=1,N_progn_pert
         !read (10) ((Progn_pert_ntrmdt_f(i,j,k), i=1,N_lonf),j=1,N_latf)
         read (10) real_tmp(:)
         Progn_pert_ntrmdt_f(:,:,k) = reshape(real_tmp,[N_lonf, N_latf])
      end do
      
      close(10)
      deallocate(real_tmp)
      rc = 0
    end subroutine i_pert_ldas
    
    ! ********************************************************************
    
    subroutine o_pert_GEOSldas(rc)
      integer,intent(inout) :: rc
      integer :: NCFOutID, STATUS
      integer :: seeddim,latdim, londim, Nforce,NProgn
      integer :: dims(3), seedid,forceid,prognid        
      integer :: xstart, ystart
      integer :: shuffle, deflate, deflate_level
      real    :: fill_value
      
      fill_value = -9999. !1.0e+15
      shuffle = 1
      deflate = 1
      deflate_level = 2
      
      !1) create file
      status = NF90_CREATE (trim(pfile_nc4), NF90_NOCLOBBER + NF90_HDF5, NCFOutID)
      
      !2) define dims
      ! status = NF_DEF_DIM(NCFOutID, 'nprogn' , N_progn_pert_max, Nprogn)
      status = NF90_DEF_DIM(NCFOutID, 'nseed' , NRANDSEED, seeddim)
      status = NF90_DEF_DIM(NCFOutID, 'lat' , N_lat, latdim)
      status = NF90_DEF_DIM(NCFOutID, 'lon' , N_lon, londim)
      status = NF90_DEF_DIM(NCFOutID, 'nforce' , N_force_pert_max, Nforce)
      status = NF90_DEF_DIM(NCFOutID, 'nprogn' , N_progn_pert_max, Nprogn)
      
      ! 3) define vars
      status = NF90_DEF_VAR(NCFOutID,'pert_rseed',NF90_DOUBLE,seeddim,seedid)
      status = NF90_PUT_ATT(NCFOutID, seedid, 'long_name','perturbations_rseed')
      status = NF90_PUT_ATT(NCFOutID, seedid, 'units', '1')
      
      dims(1)= londim
      dims(2)= latdim
      dims(3)= Nforce
      
      status = NF90_DEF_VAR(NCFOutID,'fpert_ntrmdt',NF90_REAL, dims, forceid)
      !status = nf90_def_var_deflate(NCFOutID, forceid, shuffle, deflate, deflate_level)
      status = NF90_PUT_ATT(NCFOutID, forceid, 'long_name', 'force_pert_intermediate')
      status = NF90_PUT_ATT(NCFOutID, forceid, 'units', '1')
      status = nf90_put_att(NCFOutID, forceid, '_FillValue', fill_value)
      dims(1)= londim
      dims(2)= latdim
      dims(3)= Nprogn
      
      status = NF90_DEF_VAR(NCFOutID, 'ppert_ntrmdt', NF90_REAL, dims, prognid)
      !status = nf90_def_var_deflate(NCFOutID, prognid, shuffle, deflate, deflate_level)
      status = NF90_PUT_ATT(NCFOutID, prognid, 'long_name', 'progn_pert_intermediate')
      status = NF90_PUT_ATT(NCFOutID, prognid, 'units', '1')
      status = nf90_put_att(NCFOutID, prognid, '_FillValue', fill_value)
      
      
      status = nf90_enddef(NCFOutID)
      ! 4) writing
      
      status= NF90_PUT_VAR(NCFOutID,seedid ,real(Pert_rseed,kind=8))
      
      xstart = 1 + pert_grid_f%i_offg
      ystart = 1 + pert_grid_f%j_offg
      
      ! will change to MAPL default 1.0e+15
      !do i = 1, N_lonf
      !   do j = 1, N_latf
      !      do k = 1, N_force_pert
      !         if (Force_pert_ntrmdt_f(i,j,k) < -9998) Force_pert_ntrmdt_f(i,j,k)=fill_value
      !      enddo
      !   enddo
      !enddo
      
      status= NF90_PUT_VAR(NCFOutID, forceid, Force_pert_ntrmdt_f, start=[xstart,ystart,1], &
           count=[N_lonf,N_latf,N_force_pert])
      
      ! will change to MAPL default 1.0e+15
      !do i = 1, N_lonf
      !   do j = 1, N_latf
      !      do k = 1, N_progn_pert
      !         if (Progn_pert_ntrmdt_f(i,j,k) < -9998) Progn_pert_ntrmdt_f(i,j,k)=fill_value
      !      enddo
      !   enddo
      !enddo
      
      status= NF90_PUT_VAR(NCFOutID, prognid, Progn_pert_ntrmdt_f, start=[xstart,ystart,1], &
           count=[N_lonf,N_latf,N_progn_pert])
      
      STATUS   = NF90_CLOSE (NCFOutID)
      
      deallocate(Force_pert_ntrmdt_f, Progn_pert_ntrmdt_f)
      
      rc = status
    end subroutine o_pert_GEOSldas
        
  end subroutine convert_pert_rst

  ! **************************************************************************************************
  
  subroutine LDAS_read_til_file( tile_file, catch_file, tile_grid_g, tile_coord_land, f2g )
    
    ! read land tile information from *.til file
    !
    ! This subroutine:
    !  - is the GEOSldas version of the LDASsa subroutine read_til_file() and
    !  - was known as LDAS_read_land_tile() when in LDAS_TileCoordRoutines.F90.
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
    
    character(256) :: tmpline
    character(128) :: gridname
    character(512) :: fname

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
    
    if (index(tile_grid_g%gridtype,'EASE')/=0)  ease_grid = .true.  ! 'EASEv1' or 'EASEv2'
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

             ! compute area of tile in [km^2]  (units convention in tile_coord structure)

             tile_coord(i)%area      = ease_cell_area*tile_coord(i)%frac_cell/1000./1000.  ! [km^2]
             
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
    ! pert_[x]_indg is not written into the tile_coord file and not needed in preprocessing
    tile_coord_land%pert_i_indg = nint(nodata_generic)
    tile_coord_land%pert_j_indg = nint(nodata_generic)
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
    
  contains
    
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
      
    end subroutine read_catchment_def

    ! -----------------------------------------------------------------
    
  end subroutine LDAS_read_til_file

  ! ************************************************************************************
  
end module preprocess_ldas_routines

! ====================== EOF =======================================================
  
