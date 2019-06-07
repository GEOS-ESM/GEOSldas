
module LDAS_ensdrv_vegalb_routines
  
  ! collection of LDASsa subroutines for vegetation and albedo parameters
  !
  ! (originally in clsm_ensdrv_drv_routines.F90)
  !
  ! reichle, 22 Aug 2014

  use LDAS_ensdrv_Globals,           ONLY:     &
       logunit,                                &
       logit
  
  use LDAS_DriverTypes,                 ONLY:     &
       veg_param_type,                            &
       alb_param_type
  
!  use clsm_ensdrv_mpi,                  ONLY:     &
!       MPI_INTEGER,                               &
!       MPI_COMM_WORLD,                            &
!       mpierr,                                    &
!       MPI_DATE_TIME_TYPE,                        &
!       numprocs,                                  &
!       master_proc

  use LDAS_DateTimeMod,                 ONLY:     &
       date_time_type,                            &
       datetime_lt_refdatetime,                   &
       datetime_eq_refdatetime,                   &
       datetime2_minus_datetime1,                 &
       augment_date_time,                         &
       get_dofyr_pentad

  use LDAS_ensdrv_functions,            ONLY:     &
       open_land_param_file
         
  use clsm_ensdrv_drv_routines,         ONLY:     &
       f2l_real
  
  use LDAS_ExceptionsMod,                  ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR
  
  implicit none
  
  private
  
  public :: get_veg_and_alb_times
  public :: get_VEG
  public :: get_ALB
  
contains
  
  ! ********************************************************************

  subroutine open_file_VEG( field_name, unitnum, file_format_VEG, veg_path, &
       res_ftag )
    
    implicit none

    integer,        intent(in) :: unitnum, file_format_VEG
    
    character(  3), intent(in) :: field_name
    character(200), intent(in) :: veg_path
    character( 40), intent(in) :: res_ftag
    
    ! local variables
 
    integer, parameter :: N_search_dir_max = 5

    integer            :: N_search_dir, istat

    logical            :: is_big_endian

    character(100), dimension(N_search_dir_max) :: search_dir
    character( 80)                              :: fname
    character( 20)                              :: ftag
    
    character(len=*), parameter :: Iam = 'open_file_VEG'

    ! -------------------------------------------------------------------------
    
    select case (field_name)
       
    case ('GRN');  ftag = 'green'
    case ('LAI');  ftag = 'lai'
       
    case default

       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown field_name')
       
    end select
    
    ! -----------------------------------

    select case (file_format_VEG)
       
    case (0) 
       
       N_search_dir = 2    ! specify sub-dirs of veg_path to search for file "fname"
       
       search_dir(1) = 'clsm'
       search_dir(2) = 'VEGETATION-GSWP2/LAI_GRN_CLIMATOLOGY'

       ! 'green.dat'
       ! 'lai.dat'
       
       fname = '/' // trim(ftag) // '.dat'

       is_big_endian = .true.

    case (1)
       
       N_search_dir = 1    ! specify sub-dirs of veg_path to search for file "fname"
       
       search_dir(1) = './'
       
       ! 'green_clim_180x1080.data'     - MERRA-2 on cube-sphere grid
       ! 'green_clim_540x361_DC.data'   - MERRA DC grid with MERRA-2 tiling

       ! 'lai_clim_180x1080.data'       - MERRA-2 on cube-sphere grid
       ! 'lai_clim_540x361_DC.data'     - MERRA DC grid with MERRA-2 tiling
       
       fname = '/' // trim(ftag) // '_clim_' // trim(res_ftag) // '.data'

       is_big_endian = .false.

    case default
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown file_format_VEG')
       
    end select
    
    istat = open_land_param_file(                                                       &
         unitnum, .false., is_big_endian, N_search_dir, fname, veg_path, search_dir)
    
  end subroutine open_file_VEG
    
  ! ********************************************************************
  
  subroutine open_file_ALB( field_name, unitnum, file_format_ALB, alb_path, &
       res_ftag )

    implicit none
    
    integer,        intent(in) :: unitnum, file_format_ALB
    
    character(  5), intent(in) :: field_name
    character(200), intent(in) :: alb_path
    character( 40), intent(in) :: res_ftag

    ! local variables
    
    integer, parameter :: N_search_dir_max = 5
    
    integer            :: N_search_dir, istat
    
    logical            :: is_big_endian

    character(100), dimension(N_search_dir_max) :: search_dir
    character( 80)                              :: fname
    character( 20)                              :: ftag    
    
    character(len=*), parameter :: Iam = 'open_file_ALB'

    ! -------------------------------------------------------------------------

    select case (file_format_ALB)
       
    case (0) 
       
       N_search_dir = 2    ! specify sub-dirs of alb_path to search for file "fname"
       
       search_dir(1) = 'clsm'
       search_dir(2) = 'MODIS_alb'

       select case (field_name)
          
       case ('ALBnf');  fname = '/modis_scale_factor.albnf.clim' 
       case ('ALBvf');  fname = '/modis_scale_factor.albvf.clim' 
          
       case default

          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown field_name')
          
       end select

       is_big_endian = .true.
        
    case (1)

       N_search_dir = 1    ! specify sub-dirs of alb_path to search for file "fname"
       
       search_dir(1) = './'

       ! 'nirdf_180x1080.dat'     - MERRA-2 on cube-sphere grid
       ! 'nirdf_540x361_DC.dat'   - MERRA DC grid with MERRA-2 tiling

       ! 'visdf_180x1080.dat'     - MERRA-2 on cube-sphere grid
       ! 'visdf_540x361_DC.dat'   - MERRA DC grid with MERRA-2 tiling

       ! NOTE: files named "AlbMap.*.dat" contain albedos, not the albedo *scaling*
       !       factors needed here
              
       select case (field_name)
          
       case ('ALBnf');  ftag = 'nirdf'
       case ('ALBvf');  ftag = 'visdf'

       case default

          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown field_name')
          
       end select
       
       fname = '/' // trim(ftag) // '_' // trim(res_ftag) // '.dat'
       
       is_big_endian = .false.       
       
    case default
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown file_format_ALB')
       
    end select
    
    istat = open_land_param_file(                                                      &
         unitnum, .false., is_big_endian, N_search_dir, fname, alb_path, search_dir)
    
  end subroutine open_file_ALB
    
  ! ********************************************************************
  
  subroutine get_veg_and_alb_times( N_catg, N_catf, res_ftag,                 &
       veg_path, alb_path, file_format_VEG, file_format_ALB, this_date_time,  &
       N_GRN,     N_LAI,   N_ALB,                                             &
       mid_GRN, mid_LAI, mid_ALB  )
    
    ! Read and MPI broadcast either 
    !  (i) number of data times (if "mid_*" arguments are NOT present)
    ! (ii) timestamps for veg and albedo files (otherwise)
    !
    ! reichle, 25 Jul 2013
    ! 
    ! -------------------------------------------------------------------
    
    implicit none
    
    integer,                                intent(in)    :: N_catg, N_catf
    
    character( 40),                         intent(in)    :: res_ftag
    
    character(200),                         intent(in)    :: veg_path, alb_path
    
    integer,                                intent(in)    :: file_format_VEG
    integer,                                intent(in)    :: file_format_ALB
    
    type(date_time_type),                   intent(in)    :: this_date_time

    integer,                                intent(inout) :: N_GRN
    integer,                                intent(inout) :: N_LAI
    integer,                                intent(inout) :: N_ALB
    
    type(date_time_type), dimension(N_GRN), intent(out), optional :: mid_GRN
    type(date_time_type), dimension(N_GRN), intent(out), optional :: mid_LAI
    type(date_time_type), dimension(N_GRN), intent(out), optional :: mid_ALB
    
    ! local variables
    
    integer :: unitnum
    character(len=*), parameter :: Iam = 'get_veg_and_alb_times'

    ! ------------------------------------------------------------------------

    ! ensure proper usage
    
    if ( (.not. present(mid_GRN)) .and.         &
         (.not. present(mid_LAI)) .and.         &
         (.not. present(mid_ALB))       ) then
       
       if (logit) write (logunit,*) 'reading number of data times for LAI, GRN, ALB'
       
    elseif (present(mid_GRN) .and. present(mid_LAI) .and. present(mid_ALB) ) then
       
       if (logit) write (logunit,*) 'reading midpoint times for LAI, GRN, ALB'
       
    else
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown usage')       
       
    end if
    
    unitnum = 10
    
    ! ----------------------------------

    if (master_proc) then
       
       ! greenness (GRN)
       
       call open_file_VEG( 'GRN', unitnum, file_format_VEG, veg_path, res_ftag)
       
       if (present(mid_GRN)) then
          
          call read_veg_or_alb_clim( unitnum, file_format_VEG, N_catg, N_catf, N_GRN, &
               this_date_time, mid_GRN ) 
          
       else
          
          call read_veg_or_alb_clim( unitnum, file_format_VEG, N_catg, N_catf, N_GRN ) 
          
       end if
       
       close (unitnum,status='keep')
       
       if (logit) write (logunit,*) 'done reading GRN info'
       
       ! -----------------------------------------
       !
       ! leaf area index (LAI)
       
       call open_file_VEG( 'LAI', unitnum, file_format_VEG, veg_path, res_ftag)
       
       if (present(mid_LAI)) then
          
          call read_veg_or_alb_clim( unitnum, file_format_VEG, N_catg, N_catf, N_LAI, &
               this_date_time, mid_LAI ) 
          
       else
          
          call read_veg_or_alb_clim( unitnum, file_format_VEG, N_catg, N_catf, N_LAI ) 
          
       end if
       
       close (unitnum,status='keep')
       
       if (logit) write (logunit,*) 'done reading LAI info'
       
       ! -----------------------------------------
       
       ! albedo scaling parameters (ALB)
       
       ! assume that N_times matches between ALBnf and ALBvf files
       
       call open_file_ALB( 'ALBnf', unitnum, file_format_ALB, alb_path, res_ftag )
       
       if (present(mid_ALB)) then
          
          call read_veg_or_alb_clim( unitnum, file_format_ALB, N_catg, N_catf, N_ALB, &
               this_date_time, mid_ALB ) 
          
       else
       
          call read_veg_or_alb_clim( unitnum, file_format_ALB, N_catg, N_catf, N_ALB ) 
          
       end if
       
       close (unitnum,status='keep')
       
       if (logit) write (logunit,*) 'done reading ALB info'

    end if   ! master_proc

    ! -----------------------------------------------------------------------
    !
    ! MPI broadcast (simplified "if" construct, see "proper usage" block above)

#ifdef LDAS_MPI
    
    if (.not. present(mid_GRN)) then
         
       call MPI_BCAST(N_GRN, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,mpierr)
       call MPI_BCAST(N_LAI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,mpierr)
       call MPI_BCAST(N_ALB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,mpierr)
       
    else
       
       call MPI_BCAST(mid_GRN, N_GRN, MPI_date_time_type, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(mid_LAI, N_LAI, MPI_date_time_type, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(mid_ALB, N_ALB, MPI_date_time_type, 0, MPI_COMM_WORLD, mpierr)
       
    end if
    
#endif
    
  end subroutine get_veg_and_alb_times
  
  ! ********************************************************************
  
  subroutine get_VEG( field_name, N_catg, N_catf, N_catl, f2g, N_catl_vec, low_ind,  &
       res_ftag, veg_path, file_format_VEG, this_date_time,                          &
       N_VEG, mid_VEG, veg_time_new, veg_time_old, veg_param_new, veg_param_old ) 
    
    ! Read either greenness (GRN) *or* leaf area index (LAI) data and put 
    !  into veg_param
    !
    ! field_name = 'GRN': read GRN 
    ! field_name = 'LAI': read LAI 
    !
    ! veg_time_new: first  available data time *after*  this_date_time
    ! veg_time_old: latest available data time *before* this_date_time
    ! 
    ! veg_param_new: data for veg_time_new
    ! veg_param_old: data for veg_time_old (optional, use for initialization)
    !
    ! data are read for the full domain and MPI-scattered to the local processor
    !
    ! reichle, 26 Jul 2013
    ! 
    ! -------------------------------------------------------------------
    
    implicit none
    
    character(3),                              intent(in)  :: field_name

    integer,                                   intent(in)  :: N_catg, N_catf, N_catl
    
    integer,              dimension(N_catf),   intent(in)  :: f2g

    integer,              dimension(numprocs), intent(in)  :: N_catl_vec, low_ind  
    
    character( 40),                            intent(in)  :: res_ftag    
    character(200),                            intent(in)  :: veg_path
    
    integer,                                   intent(in)  :: file_format_VEG
    
    type(date_time_type),                      intent(in)  :: this_date_time
    
    integer,                                   intent(in)  :: N_VEG
    
    type(date_time_type), dimension(N_VEG),    intent(in)  :: mid_VEG

    type(date_time_type),                      intent(out) :: veg_time_new
    type(date_time_type),                      intent(out) :: veg_time_old
    
    type(veg_param_type), dimension(N_catl),   intent(out)           :: veg_param_new
    type(veg_param_type), dimension(N_catl),   intent(out), optional :: veg_param_old
    
    ! local variables
    
    integer                                                :: unitnum, N_VEG_tmp
    
    type(date_time_type), dimension(N_VEG)                 :: mid_VEG_tmp

    real,                 dimension(N_catl)                :: data_new
    real,                 dimension(N_catl)                :: data_old

    real,                 dimension(:),        allocatable :: data_new_f
    real,                 dimension(:),        allocatable :: data_old_f
        
    ! ------------------------------------------------------------------------

    if (master_proc) then
       
       ! prepare
       
       unitnum     = 10
       
       N_VEG_tmp   = N_VEG
       
       mid_VEG_tmp = mid_VEG
       
       allocate(data_new_f(N_catf))
       
       if (present(veg_param_old))  allocate(data_old_f(N_catf))
       
       ! read full domain data
       
       call open_file_VEG( field_name, unitnum, file_format_VEG, veg_path, res_ftag)
              
       if (present(veg_param_old)) then
          
          call read_veg_or_alb_clim( unitnum, file_format_VEG, N_catg, N_catf,   &
               N_VEG_tmp, this_date_time, mid_VEG_tmp,                           &
               f2g, veg_time_new, veg_time_old, data_new_f, data_old_f )
                    
       else
          
          call read_veg_or_alb_clim( unitnum, file_format_VEG, N_catg, N_catf,   &
               N_VEG_tmp, this_date_time, mid_VEG_tmp,                           &
               f2g, veg_time_new, veg_time_old, data_new_f )
                         
       end if
       
       close (unitnum,status='keep')
    
       if (logit) write (logunit,*) 'done reading ' // field_name
       
    end if
    
    ! map from full to local domain
    
    call f2l_real( N_catf, N_catl, N_catl_vec, low_ind, data_new_f, data_new)

    select case (field_name)
       
    case ('GRN');  veg_param_new%grn = data_new
    case ('LAI');  veg_param_new%lai = data_new
       
    end select
    
    if (master_proc)  deallocate(data_new_f)

    if (present(veg_param_old)) then
       
       call f2l_real( N_catf, N_catl, N_catl_vec, low_ind, data_old_f, data_old)

       select case (field_name)
          
       case ('GRN');  veg_param_old%grn = data_old
       case ('LAI');  veg_param_old%lai = data_old
       
       end select

       if (master_proc)  deallocate(data_old_f)
       
    end if
    
#ifdef LDAS_MPI
    
    call MPI_BCAST(veg_time_new, 1, MPI_date_time_type, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(veg_time_old, 1, MPI_date_time_type, 0, MPI_COMM_WORLD, mpierr)
    
#endif
    
  end subroutine get_VEG
    
  ! ********************************************************************
  
  subroutine get_ALB( N_catg, N_catf, N_catl, f2g, N_catl_vec, low_ind,          &
       res_ftag, alb_path, file_format_ALB, this_date_time,                      &
       N_ALB, mid_ALB, alb_time_new, alb_time_old, alb_param_new, alb_param_old ) 
    
    ! Read albedo scaling parameters
    !
    ! alb_time_new: first  available data time *after*  this_date_time
    ! alb_time_old: latest available data time *before* this_date_time
    ! 
    ! alb_param_new: data for alb_time_new
    ! alb_param_old: data for alb_time_old (optional, use for initialization)
    !
    ! data are read for the full domain and MPI-scattered to the local processor
    !
    ! reichle, 26 Jul 2013
    ! 
    ! -------------------------------------------------------------------
    
    implicit none
    
    integer,                                   intent(in)  :: N_catg, N_catf, N_catl
    
    integer,              dimension(N_catf),   intent(in)  :: f2g
    
    integer,              dimension(numprocs), intent(in)  :: N_catl_vec, low_ind  
    
    character( 40),                            intent(in)  :: res_ftag
    character(200),                            intent(in)  :: alb_path
    
    integer,                                   intent(in)  :: file_format_ALB
    
    type(date_time_type),                      intent(in)  :: this_date_time
    
    integer,                                   intent(in)  :: N_ALB
    
    type(date_time_type), dimension(N_ALB),    intent(in)  :: mid_ALB
    
    type(date_time_type),                      intent(out) :: alb_time_new
    type(date_time_type),                      intent(out) :: alb_time_old
    
    type(alb_param_type), dimension(N_catl),   intent(out)           :: alb_param_new
    type(alb_param_type), dimension(N_catl),   intent(out), optional :: alb_param_old

    ! local variables
    
    integer                                                :: unitnum, N_ALB_tmp, ff

    integer,                                   parameter   :: N_fields = 2

    character(5),         dimension(N_fields)              :: field_names 
    
    type(date_time_type), dimension(N_ALB)                 :: mid_ALB_tmp

    real,                 dimension(N_catl)                :: data_new
    real,                 dimension(N_catl)                :: data_old
    
    real,                 dimension(:),        allocatable :: data_new_f
    real,                 dimension(:),        allocatable :: data_old_f
        
    ! ------------------------------------------------------------------------
    
    ! prepare

    field_names = (/ 'ALBnf', 'ALBvf' /)
    
    if (master_proc) then
       
       unitnum     = 10
       
       N_ALB_tmp   = N_ALB
       
       mid_ALB_tmp = mid_ALB

       allocate(data_new_f(N_catf))
       
       if (present(alb_param_old))  allocate(data_old_f(N_catf))
       
    end if

    ! read data

    do ff=1,N_fields
       
       if (master_proc) then
          
          ! read full domain data
          
          call open_file_ALB( field_names(ff), unitnum, file_format_ALB, alb_path,  &
               res_ftag)
          
          if (present(alb_param_old)) then
             
             call read_veg_or_alb_clim( unitnum, file_format_ALB, N_catg, N_catf,   &
                  N_ALB_tmp, this_date_time, mid_ALB_tmp,                           &
                  f2g, alb_time_new, alb_time_old, data_new_f, data_old_f )
             
          else
             
             call read_veg_or_alb_clim( unitnum, file_format_ALB, N_catg, N_catf,   &
                  N_ALB_tmp, this_date_time, mid_ALB_tmp,                           &
                  f2g, alb_time_new, alb_time_old, data_new_f )               
             
          end if
       
          close (unitnum,status='keep')
          
          if (logit) write (logunit,*) 'done reading ' // field_names(ff)
          
       end if
       
       ! map from full to local domain
    
       call f2l_real( N_catf, N_catl, N_catl_vec, low_ind, data_new_f, data_new)

       select case (field_names(ff))
          
       case ('ALBnf');  alb_param_new%sc_albnf = data_new
       case ('ALBvf');  alb_param_new%sc_albvf = data_new
          
       end select
       
       if (present(alb_param_old)) then
          
          call f2l_real( N_catf, N_catl, N_catl_vec, low_ind, data_old_f, data_old)
          
          select case (field_names(ff))
             
          case ('ALBnf');  alb_param_old%sc_albnf = data_old
          case ('ALBvf');  alb_param_old%sc_albvf = data_old
             
          end select
          
       end if
       
    end do  ! ff=1,N_fields

    if ( master_proc                              )  deallocate(data_new_f)
    if ( master_proc .and. present(alb_param_old) )  deallocate(data_old_f)

#ifdef LDAS_MPI
    
    call MPI_BCAST(alb_time_new, 1, MPI_date_time_type, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(alb_time_old, 1, MPI_date_time_type, 0, MPI_COMM_WORLD, mpierr)
    
#endif

  end subroutine get_ALB
    
  ! ********************************************************************

  subroutine read_veg_or_alb_clim( unitnum, file_format, N_catg, N_catf, N_times,  &
       this_date_time, mid_date_time, f2g, new_date_time, old_date_time, data_new, &
       data_old )
    
    ! Read climatological vegetation (LAI, greenness) or albedo scaling 
    ! parameters from file.
    !
    ! Climatological science data are provided as n-day averages for the 
    ! global tile space.
    !
    ! This subroutine accomodates the following file formats:
    !
    !  file_format=0: legacy format (monthly data, flat binaries, no date/time info)
    !  file_format=1: compatible with MAPL_readforcing()
    !
    ! The subroutine can be called in the following ways:
    !  usage=1: Obtain only N_times
    !  usage=2: Obtain mid-point date/time of data averaging intervals for all N_times
    !  usage=3: a) Read data for interval with mid-point after this_date_time   
    !           b) Read data for intervals w/ mid-points after and before 
    !               this_date_time (use for initialization)  
    !
    ! reichle, 25 Jul 2013
    ! 
    ! -------------------------------------------------------------------
    
    implicit none
    
    integer,                                  intent(in)    :: unitnum
    integer,                                  intent(in)    :: file_format

    integer,                                  intent(in)    :: N_catg, N_catf
    
    integer,                                  intent(inout) :: N_times
    
    type(date_time_type),                     intent(in),    optional :: this_date_time
    
    type(date_time_type), dimension(N_times), intent(inout), optional :: mid_date_time
    
    integer,              dimension(N_catf),  intent(in),    optional :: f2g
    
    type(date_time_type),                     intent(out),   optional :: new_date_time
    
    type(date_time_type),                     intent(out),   optional :: old_date_time
    
    real,                 dimension(N_catf),  intent(out),   optional :: data_new
    
    real,                 dimension(N_catf),  intent(out),   optional :: data_old
    
    
    ! local variables
    
    integer, parameter      :: max_times = 75  ! 73 pentads plus 2 for wrap-around
    
    integer                 :: usage, ii, jj, istat, dim1, dim2
    integer                 :: ind_new, tmp_ind_new, tmp_ind_old
    integer                 :: prev_year, curr_year, next_year

    real, dimension(14)     :: tmprealvec14

    real, dimension(N_catg) :: tmpvec

    type(date_time_type)    :: end_date_time

    type(date_time_type), dimension(:), allocatable :: start_date_time

    character(len=*), parameter :: Iam = 'read_veg_or_alb_clim'
    character(len=400) :: err_msg
    
    ! ------------------------------------------------
    !
    ! determine what is needed
    
    if (                                                        &
         (.not. present(this_date_time))  .and.                 &
         (.not. present(mid_date_time ))  .and.                 &
         (.not. present(f2g           ))  .and.                 &
         (.not. present(new_date_time ))  .and.                 &
         (.not. present(old_date_time ))  .and.                 &
         (.not. present(data_new      ))  .and.                 &
         (.not. present(data_old      ))               ) then 
       
       !  usage=1: Obtain only N_times
       
       usage = 1
       
    elseif (                                                    &
         (      present(this_date_time))  .and.                 &
         (      present(mid_date_time ))  .and.                 &
         (.not. present(f2g           ))  .and.                 &
         (.not. present(new_date_time ))  .and.                 &
         (.not. present(old_date_time ))  .and.                 &
         (.not. present(data_new      ))  .and.                 &
         (.not. present(data_old      ))               ) then 
       
       !  usage=2: Obtain mid-point date/time of data averaging intervals 
       !            for all N_times
       
       usage = 2
       
       allocate(start_date_time(N_times+1))
       
    elseif (                                                    &
         (      present(this_date_time))  .and.                 &
         (      present(mid_date_time ))  .and.                 &
         (      present(f2g           ))  .and.                 &
         (      present(new_date_time ))  .and.                 &
         (      present(old_date_time ))  .and.                 &
         (      present(data_new      ))               ) then 
       
       ! usage=3: 
       ! 
       ! a) Read data for interval with mid-point after this_date_time   
       !               (if "data_old" is NOT present)
       ! b) Read data for intervals w/ mid-points after *and* before this_date_time 
       !               (if "data_old" *is* present --> use for initialization)
       !
       ! in this usage, "mid_date_time" is intent(in)
       
       usage = 3
       
       ! Determine ind_new such that:
       !
       !   mid_date_time(ind_new-1) < this_date_time <= mid_date_time(ind_new)
       !
       ! Note:   2 <= ind_new <= N_times  (by construction)
       !         
       !         ind_old = ind_new-1      (by definition)
       
       ind_new = 2   
       
       do while (ind_new<=N_times) 
          
          ! test whether this_date_time is before mid_date_time(ind_new)
          
          if (datetime_lt_refdatetime( this_date_time, mid_date_time(ind_new)))  exit
          
          ind_new = ind_new+1
          
       end do
       
       new_date_time = mid_date_time(ind_new)
       
       old_date_time = mid_date_time(ind_new-1)
              
    else
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown usage (B)')
       
    end if

    ! ----------------------------------------------------------------------
    !
    ! read data from file
   
    select case (file_format)

    case (0) 
       
       ! file format: flat binaries, monthly data, no date/time info,
       !              wrap-around NOT stored in file 
       !              (legacy file format, e.g. MERRA-Land, Fortuna)
       
       if     (usage==1) then
          
          N_times = 14  ! incl 2 for wrap-around NOT stored in file
          
       elseif (usage==2) then

          ! get N_times+1 (!) start date/times for averaging intervals
          
          !  1 = Dec 1, 0z (year 0)
          
          start_date_time(1)%year  = 0
          start_date_time(1)%month = 12
          start_date_time(1)%day   = 1
          start_date_time(1)%hour  = 0
          start_date_time(1)%min   = 0
          start_date_time(1)%sec   = 0

          !  2 = Jan 1, 0z (year 1)
          !  3 = Feb 1, 0z (year 1)
          !  ...
          ! 13 = Dec 1, 0z (year 1)
          
          do ii=2,(N_times-1)    
             
             start_date_time(ii)%year  = 1
             start_date_time(ii)%month = ii-1
             start_date_time(ii)%day   = 1
             start_date_time(ii)%hour  = 0
             start_date_time(ii)%min   = 0
             start_date_time(ii)%sec   = 0
             
          end do

          ! 14 = Jan 1, 0z (year 2)

          start_date_time(N_times  )%year  = 2
          start_date_time(N_times  )%month = 1
          start_date_time(N_times  )%day   = 1
          start_date_time(N_times  )%hour  = 0
          start_date_time(N_times  )%min   = 0
          start_date_time(N_times  )%sec   = 0

          ! 15 = Feb 1, 0z (year 2)

          start_date_time(N_times+1)%year  = 2
          start_date_time(N_times+1)%month = 2
          start_date_time(N_times+1)%day   = 1
          start_date_time(N_times+1)%hour  = 0
          start_date_time(N_times+1)%min   = 0
          start_date_time(N_times+1)%sec   = 0
          
       elseif (usage==3) then
          
          ! translate ind_new into tmp_ind_new (wrap-around months NOT stored in file)

          if (ind_new==2 .or. ind_new==14) then
             
             tmp_ind_new = 1       ! Jan
             tmp_ind_old = 12      ! Dec

          else
             
             tmp_ind_new = ind_new - 1
             tmp_ind_old = ind_new - 2
             
          end if

          ! disable "tmp_ind_old" if not needed

          if (.not. present(data_old))  tmp_ind_old = -9999
          
          ! read through file and extract months of interest
          
          do ii=1,max(tmp_ind_old,tmp_ind_new)
             
             if     (ii==tmp_ind_new) then
                
                read (unitnum) (tmpvec(jj), jj=1,N_catg)
                
                data_new(1:N_catf) = tmpvec(f2g(1:N_catf)) 
                
             elseif (ii==tmp_ind_old) then
                
                ! per definition of tmp_ind_old above and loop boundaries,
                ! "data_old" must be present in this case
                
                read (unitnum) (tmpvec(jj), jj=1,N_catg)
                
                data_old(1:N_catf) = tmpvec(f2g(1:N_catf)) 
                
             else
                
                read (unitnum)  ! SKIP science data record
                
             end if
             
          end do
          
       end if ! usage
              
       ! -------------------------------------------
       
    case (1)

       ! file format: compatible with MAPL_readforcing()
       !              flat binaries, n-day averages, date/time info
       !              (e.g., MERRA-2)

       if     (usage==1) then
          
          ! determine number of data records in file
          
          ii = 0
          
          do while (ii<max_times+1)    ! limit number of read attempts
             
             read (unitnum, iostat=istat)       ! SKIP date/time record
             
             if (istat/=0)  exit      
             
             read (unitnum, iostat=istat)       ! SKIP science data record
             
             if (istat/=0) then
                err_msg = 'odd number of records in "MAPL_readforcing" file'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if

             ii = ii+1
             
          end do
          
          if (ii>=max_times) then
             err_msg = 'number or data times in file exceeds max allowed'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if

          N_times = ii
   
       elseif (usage==2) then
                    
          do ii=1,N_times
             
             ! read date/time record
          
             read (unitnum) (tmprealvec14(jj), jj=1,14)
             
             ! start date/time of averaging interval
             
             start_date_time(ii)%year  = nint(tmprealvec14( 1))
             start_date_time(ii)%month = nint(tmprealvec14( 2))
             start_date_time(ii)%day   = nint(tmprealvec14( 3))
             start_date_time(ii)%hour  = nint(tmprealvec14( 4)) 
             start_date_time(ii)%min   = nint(tmprealvec14( 5)) 
             start_date_time(ii)%sec   = nint(tmprealvec14( 6)) 
             
             ! sanity check
             
             if (ii>1) then
                
                ! start of current interval must match end of previous interval
                
                if (.not. datetime_eq_refdatetime(                                    &
                     start_date_time(ii), end_date_time )                             &
                     ) then
                   err_msg = 'intervals do not line up'
                   call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                end if
                
             end if
             
             ! end date/time of averaging interval
             
             end_date_time%year        = nint(tmprealvec14( 7))
             end_date_time%month       = nint(tmprealvec14( 8))
             end_date_time%day         = nint(tmprealvec14( 9))
             end_date_time%hour        = nint(tmprealvec14(10)) 
             end_date_time%min         = nint(tmprealvec14(11)) 
             end_date_time%sec         = nint(tmprealvec14(12)) 
             
             ! spatial dimensions
             
             dim1                      = nint(tmprealvec14(13)) 
             dim2                      = nint(tmprealvec14(14)) 
             
             ! sanity check
             
             ! dim1 must match N_catg, dim2 must be 1
             
             if ((dim1/=N_catg) .or. (dim2/=1)) then
                err_msg = 'dimensions do not match'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if
             
             
             read (unitnum)        ! SKIP science data record
             
          end do
          
          ! fill in last element of start_date_time

          start_date_time(N_times+1) = end_date_time
          
          ! additional sanity checks
          
          ! check wrap-around: last three start_date_time entries must match 
          !  first three, resp. (except for year)
          
          if ( (start_date_time(N_times-1)%month/=start_date_time(1)%month) .or.  &
               (start_date_time(N_times-1)%day  /=start_date_time(1)%day  ) .or.  &
               (start_date_time(N_times-1)%hour /=start_date_time(1)%hour ) .or.  &
               (start_date_time(N_times-1)%min  /=start_date_time(1)%min  ) .or.  &
               (start_date_time(N_times-1)%sec  /=start_date_time(1)%sec  )     ) then
             err_msg = 'something wrong with wrap-around (A)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          if ( (start_date_time(N_times  )%month/=start_date_time(2)%month) .or.  &
               (start_date_time(N_times  )%day  /=start_date_time(2)%day  ) .or.  &
               (start_date_time(N_times  )%hour /=start_date_time(2)%hour ) .or.  &
               (start_date_time(N_times  )%min  /=start_date_time(2)%min  ) .or.  &
               (start_date_time(N_times  )%sec  /=start_date_time(2)%sec  )     ) then
             err_msg = 'something wrong with wrap-around (B)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          if ( (start_date_time(N_times+1)%month/=start_date_time(3)%month) .or.  &
               (start_date_time(N_times+1)%day  /=start_date_time(3)%day  ) .or.  &
               (start_date_time(N_times+1)%hour /=start_date_time(3)%hour ) .or.  &
               (start_date_time(N_times+1)%min  /=start_date_time(3)%min  ) .or.  &
               (start_date_time(N_times+1)%sec  /=start_date_time(3)%sec  )     ) then
             err_msg = 'something wrong with wrap-around (C)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          ! check years
          
          prev_year = start_date_time(        1)%year
          curr_year = start_date_time(        2)%year
          next_year = start_date_time(  N_times)%year
          
          if ((prev_year+1/=curr_year) .or. (curr_year+1/=next_year)) then
             err_msg = 'error with years in file (A)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          do ii=3,N_times-1
             
             if (start_date_time(ii)%year/=curr_year) then
                err_msg = 'error with years in file (B)'
             end if

          end do
          
          if (start_date_time(N_times+1)%year/=next_year) then
             err_msg = 'error with years in file (C)'
          end if
          
          
       elseif (usage==3) then
          
          do ii=1,ind_new    

             read (unitnum)     ! SKIP date/time info record

             if (ii<ind_new-1) then
                
                read (unitnum)  ! SKIP science data record
                
             else
                
                read (unitnum) (tmpvec(jj), jj=1,N_catg)
                
             end if
             
             if (ii==ind_new-1) then   ! (ii==ind_old)
                
                if (present(data_old))                              &
                     data_old(1:N_catf) = tmpvec(f2g(1:N_catf)) 
                
             else                      ! (ii==ind_new)
                
                data_new(1:N_catf) = tmpvec(f2g(1:N_catf)) 
                
             end if
             
          end do
          
       end if ! usage
       
       ! ------------------------------------------
       
    case default
       
       err_msg = 'unknown file format in subroutine'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
    ! ----------------------------------------------------------------------
    !
    ! get mid-point date/time

    if (usage==2) then
       
       call get_veg_or_alb_clim_mid_date_time( this_date_time%year, N_times, &
            start_date_time, mid_date_time)
       
    end if
    
  end subroutine read_veg_or_alb_clim
  
  ! ********************************************************************
  
  subroutine get_veg_or_alb_clim_mid_date_time( current_year, N_times, &
       start_date_time, mid_date_time)
    
    implicit none
    
    ! compute mid-point date/time of averaging intervals for vegetation or albedo
    ! climatologies

    ! current year overwrites year provided with the data and thereby takes
    ! leap years into account
    
    ! reichle, 23 July 2013
    
    integer, intent(in)  :: current_year, N_times
    
    type(date_time_type), dimension(N_times+1), intent(in)  :: start_date_time 
    
    type(date_time_type), dimension(N_times),   intent(out) :: mid_date_time 
    
    ! local variables

    integer :: ii, dt

    type(date_time_type), dimension(N_times+1)  :: start_date_time_tmp
    
    ! -----------------------------------------------------
    
    start_date_time_tmp = start_date_time

    ! prepare
    
    do ii=1,N_times+1
       
       ! match years to previous, current, and next year
       
       if     (ii==1)                   then
          
          start_date_time_tmp(ii)%year = current_year-1 
          
       elseif (ii>1 .and. ii<N_times)   then
          
          start_date_time_tmp(ii)%year = current_year   

       elseif (ii>=N_times)             then
                 
          start_date_time_tmp(ii)%year = current_year+1  

       end if
       
       ! recompute day-of-year
       
       call get_dofyr_pentad( start_date_time_tmp(ii) )  
       
    end do
   
    ! compute mid-point date/time
    
    do ii=1,N_times
       
       ! get length of interval (ii:ii+1) in seconds
       
       dt = datetime2_minus_datetime1(                            &
            start_date_time_tmp(ii), start_date_time_tmp(ii+1) )
       
       ! initialize and add dt/2
       
       mid_date_time(ii) = start_date_time_tmp(ii)
       
       call augment_date_time( dt/2, mid_date_time(ii) ) 
       
    end do
    
  end subroutine get_veg_or_alb_clim_mid_date_time

  ! ********************************************************************

end module clsm_ensdrv_vegalb_routines

! *********** EOF ******************************************************
