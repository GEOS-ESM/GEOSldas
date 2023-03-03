! 
! this file contains subroutines for reading and processing observations 
! for the GEOS5 land EnKF update algorithm
!
! reichle, 27 Jan 2005
! reichle, 10 Jan 2011 - replaced "UVA" with "LPRM"
! reichle,  1 Jul 2015 - clarified definition of obs time stamp 
!                          (J2000 seconds w/ 'TT12' epoch)

! added work_path to inputs of many subroutines so that "tmpfname"
! (needed several times for reading AMSR-E hdf files) is distinct for each job
! reichle, 27 Aug 2005

! *********************************************************************

module clsm_ensupd_read_obs
  
  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_TICE

  use io_hdf5,                          ONLY:     &
       hdf5read 

  use EASE_conv,                        ONLY:     &
       ease_convert,                              &
       ease_extent
  
  use LDAS_ensdrv_globals,              ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_tolfrac_generic
  
  use clsm_ensupd_glob_param,           ONLY:     &
       unitnum_obslog

  use LDAS_DateTimeMod,                 ONLY:     &
       date_time_type,                            &
       augment_date_time,                         &
       get_dofyr_pentad,                          &
       datetime_le_refdatetime,                   &
       datetime_lt_refdatetime,                   &
       date_time2string,                          &
       datetime_to_J2000seconds

  use enkf_types,                       ONLY:     &
       obs_type,                                  &
       obs_param_type

  use LDAS_TilecoordType,               ONLY:     &
       tile_coord_type,                           &
       grid_def_type

  use clsm_ensdrv_drv_routines,         ONLY:     &
       f2l_real,                                  &
       f2l_real8,                                 &
       f2l_logical
  
  use LDAS_TilecoordRoutines,           ONLY:     &
       get_tile_num_from_latlon
  
  use LDAS_ensdrv_mpi,                  ONLY:     &
       root_proc,                                 &
       numprocs,                                  &
       mpicomm,                                   &
       MPI_obs_type,                              &
       mpierr

  use LDAS_exceptionsMod,               ONLY:     &
       ldas_abort,                                &
       ldas_warn,                                 &
       LDAS_GENERIC_ERROR,                        &
       LDAS_GENERIC_WARNING

  use clsm_ensupd_upd_routines,         ONLY:     &
       dist_km2deg
  
  implicit none

  include 'mpif.h'  
  
  private
  
  public :: collect_obs  

contains
  
  ! *****************************************************************
  
  subroutine read_ae_l2_sm_hdf( &
       N_files, fnames, N_data, lon, lat, ae_l2_sm, ease_col, ease_row )
    
    ! read soil moisture data from one or more AMSR-E Land hdf files
    !
    ! return ONLY valid data points (ie. excluding no-data-values)
    ! that also pass initial QC (based on "Surface Type" QC flag and on 
    ! Heterogeneity_Index )
    !
    ! reichle, 20 Sep 2005 - added "Surface Type" QC flag
    ! reichle, 17 Nov 2005 - replace hdp with call to hdf library
    ! reichle,  8 Feb 2006 - optionally read EASE row- and column index
    !                      - added Heterogeneity_Index QC
    ! reichle, 10 Jan 2011 - revised "Surface Type" QC flag

    implicit none
    
    integer, intent(in) :: N_files
    
    character(*), dimension(N_files), intent(in) :: fnames
    
    integer, intent(out) :: N_data
    
    real, dimension(:), pointer :: lon, lat, ae_l2_sm  ! output
    
    integer, dimension(:), pointer, optional :: ease_col, ease_row ! output
    
    ! local parameters
    
    integer, parameter :: N_fields = 7
    
    character(19), dimension(N_fields), parameter :: field_names = (/ &
         'Longitude          ',        &      ! 1
         'Latitude           ',        &      ! 2
         'Surface_Type       ',        &      ! 3
         'Soil_Moisture      ',        &      ! 4
         'Row_Index          ',        &      ! 5
         'Column_Index       ',        &      ! 6
         'Heterogeneity_Index'  /)            ! 7
    
    real, parameter :: scale_fac_Soil_Moisture       = 1000.0;
    
    integer, parameter :: nodata = -9999     ! NOTE: integer
    
    ! Initial QC:
    !
    ! "Surface_Type": Only "low vegetation" or "moderate vegetation" (and no precip, 
    ! frozen ground, etc) passes.  
    !
    ! "Heterogeneity_Index": Discard all pixels with 
    !       Heterogeneity_Index>max_Heterogeneity_Index 
    ! because they are likely mixed land/water pixels.  See matlab code
    ! "detect_coast_in_AMSRE_Land.m"
    ! in land01:/home/reichle/NSIPP/AMSR/AMSR_E_Land/matlab/
    ! In this subroutine, Heterogeneity_Index is used on its raw form and 
    ! never scaled into units of Kelvin.
    !
    ! Further info:
    ! http://nsidc.org/data/docs/daac/ae_land_l2b_soil_moisture.gd.html
    ! http://nsidc.org/data/amsre/versions.html
    ! -reichle, 20 Sep 2005
    ! -reichle,  8 Feb 2006
    
    integer, parameter :: max_Heterogeneity_Index      = 500 ! = 5 Kelvin
    
    ! declarations of hdf functions 
    
    integer :: hopen, vfstart, vsfatch, vsqfnelt, vsfseek, vsfsfld, vsfread
    integer :: vsfdtch, vfend, hclose
    
    
    ! declarations of hdf-related parameters and variables
    
    integer, dimension(N_files) :: file_id, vdata_id 
    
    integer :: status, n_read, n_records, record_pos
    
    integer, parameter :: num_dds_block = 0  ! only important for writing hdf
    
    integer, parameter :: vdata_ref = 7      ! works for AMSRE_L2_Land
    
    integer, parameter :: DFACC_READ     = 1 ! from hdf.inc
    integer, parameter :: FULL_INTERLACE = 0 ! from hdf.inc
        
    ! local variables

    logical :: must_stop
    
    integer, dimension(N_files) :: N_data_tmp
    
    integer :: i, j, k, k_off
    
    integer, dimension(:), allocatable :: surface_type_qc_flag
    integer, dimension(:), allocatable :: Heterogeneity_Index
    
    integer*2, dimension(:), allocatable :: tmpint2vec
    real,      dimension(:), allocatable :: tmprealvec
    
    character(len=*), parameter :: Iam = 'read_ae_l2_sm_hdf'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------
    
    ! determine number of data to be read from each file
    
    do j=1,N_files
       
       ! open and "start" hdf file 
       
       file_id(j) = hopen( fnames(j), DFACC_READ, num_dds_block )  
       
       status = vfstart(file_id(j))
         
       ! select vdata block that contains fields of interest
       
       vdata_id(j) = vsfatch(file_id(j), vdata_ref, 'r') 
       
       ! determine number of records in vdata
       
       status = vsqfnelt(vdata_id(j), n_records)
       
       N_data_tmp(j) = n_records
       
    end do
    
    ! allocate pointers (must be deallocated outside this subroutine)
    
    must_stop = .false.
    
    if ( associated(lon) .or. associated(lat) .or. associated(ae_l2_sm) ) then
       must_stop = .true.
    end if
    
    if ( present(ease_col) ) then         
       if (associated(ease_col))  must_stop = .true.
    end if
    
    if ( present(ease_row) ) then
       if (associated(ease_row))  must_stop = .true.
    end if
    
    if (must_stop) then      
       err_msg = 'output pointers must not be associated/allocated on input.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    N_data = sum(N_data_tmp)
    
    allocate(lon(N_data))
    allocate(lat(N_data))
    if (present(ease_col))  allocate(ease_col(N_data))
    if (present(ease_row))  allocate(ease_row(N_data))
    allocate(ae_l2_sm(N_data))
    
    allocate(surface_type_qc_flag(N_data))
    allocate(Heterogeneity_Index(N_data))

    ! read hdf data into arrays, concatenate data from N_files files
    
    k_off = 0
    
    do j=1,N_files
       
       allocate(tmprealvec(N_data_tmp(j)))
       allocate(tmpint2vec(N_data_tmp(j)))
       
       do i=1,N_fields
          
          ! go to start of record (zero-based count)
          
          record_pos = vsfseek(vdata_id(j), 0)
          
          ! pick the field to be read
          
          status = vsfsfld(vdata_id(j), field_names(i))
          
          ! read data
          
          select case (i)
             
          case (1)
             
             n_read = vsfread( vdata_id(j), tmprealvec, &
                  N_data_tmp(j), FULL_INTERLACE)
             
             lon(k_off+1:k_off+N_data_tmp(j)) = tmprealvec
             
          case (2)
             
             n_read = vsfread( vdata_id(j), tmprealvec, &
                  N_data_tmp(j), FULL_INTERLACE)
             
             lat(k_off+1:k_off+N_data_tmp(j)) = tmprealvec
             
          case (3)
             
             n_read = vsfread( vdata_id(j) ,tmpint2vec, &
                  N_data_tmp(j), FULL_INTERLACE)

             surface_type_qc_flag(k_off+1:k_off+N_data_tmp(j)) = tmpint2vec
             
             ! overwrite surface_type_qc_flag with common pass/fail
             ! (negative number -1 means fail)
             
             do k=1,N_data_tmp(j)

                ! keep only data with
                !
                !  surface_type_qc_flag=128  ("moderate veg", and no other bits set)
                !  surface_type_qc_flag=256  ("low veg",      and no other bits set)
                !
                ! http://nsidc.org/data/docs/daac/ae_land_l2b_soil_moisture.gd.html
                
                if (.not. (                                           &
                     (surface_type_qc_flag(k+k_off)==128)  .or.       &
                     (surface_type_qc_flag(k+k_off)==256)       )  )  &
                     surface_type_qc_flag(k+k_off) = -1
                
             end do
             
          case (4)
             
             n_read = vsfread(vdata_id(j), tmpint2vec, &
                  N_data_tmp(j), FULL_INTERLACE)
             
             do k=1,N_data_tmp(j)
                if (tmpint2vec(k)/=nodata) then
                   ae_l2_sm(k+k_off) = real(tmpint2vec(k)) / &
                        scale_fac_Soil_Moisture
                else
                   ae_l2_sm(k+k_off) = real(tmpint2vec(k))
                end if
             end do
             
          case (5)
             
             if (present(ease_row)) then
                
                n_read = vsfread( vdata_id(j) ,tmpint2vec, &
                     N_data_tmp(j), FULL_INTERLACE)
                
                ease_row(k_off+1:k_off+N_data_tmp(j)) = tmpint2vec
                
             end if
             
          case (6)
             
             if (present(ease_col)) then
                
                n_read = vsfread( vdata_id(j) ,tmpint2vec, &
                     N_data_tmp(j), FULL_INTERLACE)
                
                ease_col(k_off+1:k_off+N_data_tmp(j)) = tmpint2vec
                
             end if
             
          case (7)
             
             n_read = vsfread( vdata_id(j) ,tmpint2vec, &
                  N_data_tmp(j), FULL_INTERLACE)
             
             Heterogeneity_Index(k_off+1:k_off+N_data_tmp(j)) = tmpint2vec
             
          case default
             
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown case')

          end select
          
          if (n_read/=N_data_tmp(j)) then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'ERROR reading hdf')
          end if
          
       end do
       
       ! clean up

       deallocate(tmprealvec)
       deallocate(tmpint2vec)
       
       ! close hdf files
       
       status = vsfdtch(vdata_id(j))
  
       status = vfend(file_id(j))
       
       status = hclose(file_id(j))
       
       ! prepare next j
       
       k_off = k_off + N_data_tmp(j)
       
    end do
    
    ! -------------------------------------
    !
    ! eliminate no-data-values and data that fail initial QC
    
    j = 0
    
    do i=1,N_data
       
       if ( (ae_l2_sm(i)>0.)             .and.         &  ! any neg is nodata
            (surface_type_qc_flag(i)>=0) .and.         &  ! any neg is fail
            (Heterogeneity_Index(i)<=max_Heterogeneity_Index) &
            ) then     
          
          j=j+1
          
          ae_l2_sm(j) = ae_l2_sm(i)
          lon(j)      = lon(i)
          lat(j)      = lat(i)
          if (present(ease_col))  ease_col(j) = ease_col(i)
          if (present(ease_row))  ease_row(j) = ease_row(i)
          
       end if
       
    end do
    
    N_data = j
    
    deallocate(surface_type_qc_flag)
    deallocate(Heterogeneity_Index)

  end subroutine read_ae_l2_sm_hdf
    
  ! *****************************************************************
  
  subroutine read_obs_ae_l2_sm(                                  &
       work_path, exp_id,                                        &
       date_time, dtstep_assim, N_catd, tile_coord,              &
       tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
       this_obs_param,                                           &
       found_obs, ae_l2_sm, std_ae_l2_sm )
    
    ! Read observations of surface soil moisture from AMSR-E L2 files
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation.
    !
    ! If there are N > 1 observations in a given tile,
    ! a "super-observation" is computed by averaging the N observations
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! reichle, 25 Jul 2005
    !
    ! --------------------------------------------------------------------
        
    implicit none
    
    ! inputs:
    
    character(*),  intent(in) :: work_path
    character(*),  intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: dtstep_assim, N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: ae_l2_sm
    real,    intent(out), dimension(N_catd) :: std_ae_l2_sm
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals
    
    ! AMSR_E_L2_Land files are available for approximately 50min periods
    ! covering one ascending or descending swath.  The filename indicates
    ! the start time of the swath
    ! "ae_time_offset" is used to find the mean time of the the interval
    ! which is approximately the time of the equator overpass.
    ! This time is assigned to all observations of the swath.
    
    integer, parameter :: ae_time_offset = 1500   ! 25 minutes in seconds
    
    character(4)   :: DDHH
    character(6)   :: YYYYMM
    character(8)   :: date_string
    character(10)  :: time_string
    character(300) :: tmpfname, tmpfname2
    character(400) :: cmd

    type(date_time_type) :: date_time_tmp
    
    integer :: i, ind, N_tmp, N_files

    character(300), dimension(:), allocatable :: fnames
    
    real,    dimension(:), pointer :: tmp_obs, tmp_lat, tmp_lon
    integer, dimension(:), pointer :: tmp_tile_num
    
    integer, dimension(N_catd) :: N_obs_in_tile    

    character(len=*), parameter :: Iam = 'read_obs_ae_l2_sm'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------------
    
    nullify( tmp_obs, tmp_lat, tmp_lon, tmp_tile_num )    
    
    ! ---------------
    
    ! initialize
    
    found_obs = .false.

    ! find files that are within half-open interval 
    ! [date_time-dtstep_assim/2,date_time+dtstep_assim/2)
    
    date_time_tmp = date_time    
    
    call augment_date_time( -(dtstep_assim/2 + ae_time_offset), date_time_tmp )
    
    ! get tmp file name and remove file if it exists
    
    call date_and_time(date_string, time_string)  ! f90 intrinsic function
    
    tmpfname = trim(work_path) // '/' // 'tmp.' // trim(exp_id) &
         // '.' // date_string // time_string

    cmd = '/bin/rm -f ' // tmpfname 
    
    call Execute_command_line(trim(cmd))
    
    ! identify all files within current assimilation interval
    ! (list all files within hourly intervals)
    
    do i=1,(dtstep_assim/3600)
       
       write (YYYYMM,'(i6.6)') date_time_tmp%year*100 + date_time_tmp%month
       write (DDHH,  '(i4.4)') date_time_tmp%day *100 + date_time_tmp%hour
       
       cmd = 'ls ' // trim(this_obs_param%path) // '/' // YYYYMM(1:4) // &
            '/M' // YYYYMM(5:6) // '/' // trim(this_obs_param%name) &
            // '*' // YYYYMM // DDHH // '*'
       
       if     (trim(this_obs_param%descr)=='ae_l2_sm_a') then
          
          cmd = trim(cmd) // '_A.hdf'
          
       elseif (trim(this_obs_param%descr)=='ae_l2_sm_d') then
          
          cmd = trim(cmd) // '_D.hdf'
          
       else
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown descr')

       end if
       
       cmd = trim(cmd) // ' >> ' // trim(tmpfname)
       
       call Execute_command_line(trim(cmd))
       
       call augment_date_time( 3600, date_time_tmp )
       
    end do
    
    ! find out how many need to be read
    
    tmpfname2 = trim(tmpfname) // '.wc'
    
    cmd = 'wc -w ' // trim(tmpfname) // ' > ' // trim(tmpfname2)
    
    call Execute_command_line(trim(cmd))
    
    open(10, file=tmpfname2, form='formatted', action='read')
    
    read(10,*) N_files
    
    close(10,status='delete')
    
    ! load file names into "fnames"
    
    open(10, file=tmpfname,  form='formatted', action='read')
    
    if (N_files>0) then
       
       allocate(fnames(N_files))
       
       do i=1,N_files
          read(10,'(a)') fnames(i)
       end do
       
    end if
    
    close(10,status='delete')
    
    ! read observations:
    !
    ! 1.) read N_tmp observations and their lat/lon info from file
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon
    ! 3.) compute super-obs for each tile from all obs w/in that tile
    !
    ! ----------------------------------------------------------------
    !
    ! 1.) read N_tmp observations and their lat/lon info from file

    if (N_files>0) then

       call read_ae_l2_sm_hdf( &
            N_files, fnames, &
            N_tmp, tmp_lon, tmp_lat, tmp_obs )
              
       if (logit) then
          
          write (logunit,*) 'read_obs_ae_l2_sm: read ', N_tmp,  &
               ' at date_time = ', date_time, ' from '
          do i=1,N_files
             write (logunit,*) trim(fnames(i))
          end do
          write (logunit,*) '----------'
          
       end if

       deallocate(fnames)

    else
       
       N_tmp = 0
       
    end if

    ! ------------------------------------------------------------------

    ! note QC and no-data-value block in subroutine read_ae_l2_sm_hdf()
    
    ! ------------------------------------------------------------------
    !
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon

    if (N_tmp>0) then
       
       allocate(tmp_tile_num(N_tmp))
       
       call get_tile_num_for_obs(N_catd, tile_coord,                 &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
            N_tmp, tmp_lat, tmp_lon,                                 &
            this_obs_param,                                          &
            tmp_tile_num )
       
       
       ! ----------------------------------------------------------------
       !
       ! 3.) compute super-obs for each tile from all obs w/in that tile
       !     (also eliminate observations that are not in domain)
       
       ae_l2_sm = 0.
       N_obs_in_tile  = 0
       
       do i=1,N_tmp
          
          ind = tmp_tile_num(i)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
          
          if (ind>0) then         ! this step eliminates obs outside domain
             
             ae_l2_sm(ind) = ae_l2_sm(ind) + tmp_obs(i)
             
             N_obs_in_tile(ind) = N_obs_in_tile(ind) + 1
             
          end if
          
       end do
       
       ! normalize
       
       do i=1,N_catd
          
          if (N_obs_in_tile(i)>1) then
             
             ae_l2_sm(i) = ae_l2_sm(i)/real(N_obs_in_tile(i))
          
          elseif (N_obs_in_tile(i)==0) then
             
             ae_l2_sm(i) = this_obs_param%nodata
             
          end if
          
       end do
       
       ! clean up
       
       if (associated(tmp_tile_num)) deallocate(tmp_tile_num)
       
       ! --------------------------------
       
       ! set observation error standard deviation
       
       do i=1,N_catd
          std_ae_l2_sm(i) = this_obs_param%errstd
       end do
       
       ! --------------------------------
       
       if (any(N_obs_in_tile>0)) then
          
          found_obs = .true.
          
       else 
          
          found_obs = .false.
          
       end if
       
    end if
    
    ! clean up
    
    if (associated(tmp_obs))      deallocate(tmp_obs)
    if (associated(tmp_lon))      deallocate(tmp_lon)
    if (associated(tmp_lat))      deallocate(tmp_lat)

  end subroutine read_obs_ae_l2_sm
  
  ! ***************************************************************************

  subroutine read_ae_sm_LPRM_bin( &
       this_obs_param, N_files, fnames, N_data, lon, lat, ae_sm_LPRM, ease_col, ease_row )
    
    ! read soil moisture data from one or more AMSR-E LPRM bin files
    !
    ! return ONLY valid data points (ie. excluding no-data-values)
    ! 
    ! no QC in addition to what was done in matlab-preprocessing
    !
    ! reichle, 20 Feb 2009
    
    implicit none
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    integer, intent(in) :: N_files
    
    character(*), dimension(N_files), intent(in) :: fnames
    
    integer, intent(out) :: N_data
    
    real, dimension(:), pointer :: lon, lat, ae_sm_LPRM  ! output
    
    integer, dimension(:), pointer, optional :: ease_col, ease_row ! output
    
    ! local variables
    
    logical :: must_stop
    
    integer, dimension(N_files) :: N_data_tmp
    
    integer :: i, j, k_off
    
    integer, dimension(:), allocatable :: tmpintvec
    real,    dimension(:), allocatable :: tmprealvec

    character(len=*), parameter :: Iam = 'read_ae_sm_LPRM_bin'
    character(len=400) :: err_msg
    
    ! -------------------------------------------------------------
    
    ! make sure pointers are not allocated or associated
    
    must_stop = .false.
    
    if ( associated(lon) .or. associated(lat) .or. associated(ae_sm_LPRM) ) then
       must_stop = .true.
    end if
    
    if ( present(ease_col) ) then         
       if (associated(ease_col))  must_stop = .true.
    end if
    
    if ( present(ease_row) ) then
       if (associated(ease_row))  must_stop = .true.
    end if
    
    if (must_stop) then
       err_msg = 'output pointers must not be associated/allocated on input'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! determine number of data to be read from each file
    
    N_data = 0
    
    do j=1,N_files
       
       ! open file 
       
       open( 10, file=trim(fnames(j)), form='unformatted',convert='big_endian', action='read' )
       
       read( 10) N_data_tmp(j)
       
       close(10, status='keep')
       
    end do
    
    ! allocate pointers (must be deallocated outside this subroutine!)
    
    N_data = sum(N_data_tmp)
    
    allocate(lon(N_data))
    allocate(lat(N_data))
    if (present(ease_col))  allocate(ease_col(N_data))
    if (present(ease_row))  allocate(ease_row(N_data))
    allocate(ae_sm_LPRM(N_data))
    
    ! read data into arrays, concatenate data from N_files files
    
    ! format of AMSR_sm_LPRM_EASE_bin files:
    !
    ! record  1 -- N_data            int*4
    ! record  2 -- lon(  1:N_data)   real*4    
    ! record  3 -- lat(  1:N_data)   real*4    
    ! record  4 -- sm_C( 1:N_data)   real*4    
    ! record  5 -- sm_X( 1:N_data)   real*4    
    ! record  6 -- od_C( 1:N_data)   real*4    
    ! record  7 -- od_X( 1:N_data)   real*4    
    ! record  8 -- res_C(1:N_data)   real*4    
    ! record  9 -- res_X(1:N_data)   real*4    
    ! record 10 -- ts_C( 1:N_data)   real*4       
    ! record 11 -- ts_X( 1:N_data)   real*4        
    ! record 12 -- rfi_C(1:N_data)   int*4    
    ! record 13 -- rfi_X(1:N_data)   int*4    
    ! record 14 -- ind_i(1:N_data)   int*4    zero-based (!) EASE row index
    ! record 15 -- ind_j(1:N_data)   int*4    zero-based (!) EASE col index
    ! record 16 -- time( 1:N_data)   real*4   minutes since beginning of half-orbit
        
    k_off = 0
    
    do j=1,N_files
       
       allocate(tmprealvec(N_data_tmp(j)))
       
       if (present(ease_col))  allocate(tmpintvec(N_data_tmp(j)))
       
       open (10, file=trim(fnames(j)),  form='unformatted',convert='big_endian', action='read' )
       
       ! re-read N_data
       
       read (10) N_data_tmp(j)       
       
       ! read data as needed 
       
       read (10) tmprealvec; lon(k_off+1:k_off+N_data_tmp(j)) = tmprealvec    
       read (10) tmprealvec; lat(k_off+1:k_off+N_data_tmp(j)) = tmprealvec    
       
       if      (this_obs_param%descr(13:14)=='_C') then
          
          read (10) tmprealvec; ae_sm_LPRM(k_off+1:k_off+N_data_tmp(j)) = tmprealvec
          read (10) ! skip sm_X
          
       else if (this_obs_param%descr(13:14)=='_X') then
          
          read (10) ! skip sm_C
          read (10) tmprealvec; ae_sm_LPRM(k_off+1:k_off+N_data_tmp(j)) = tmprealvec
          
       else          
          
          err_msg = 'unknown descr, ' // this_obs_param%descr
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
       if (present(ease_col) .and. present(ease_row)) then
          
          read (10) ! skip od_C
          read (10) ! skip od_X
          read (10) ! skip res_C
          read (10) ! skip res_X
          read (10) ! skip ts_C
          read (10) ! skip ts_X
          read (10) ! skip rfi_C
          read (10) ! skip rfi_X
          
          read (10) tmpintvec; ease_col(k_off+1:k_off+N_data_tmp(j)) = tmpintvec
          read (10) tmpintvec; ease_row(k_off+1:k_off+N_data_tmp(j)) = tmpintvec
          
          ! read (10)  ! skip time
          
       end if
       
       ! clean up
       
       close(10, status='keep')
       
       deallocate(tmprealvec)
       if (allocated(tmpintvec))  deallocate(tmpintvec)
       
       ! prepare next j
       
       k_off = k_off + N_data_tmp(j)
       
    end do
    
    ! -------------------------------------
    !
    ! eliminate no-data-values 
    
    j = 0
    
    do i=1,N_data
       
       if (ae_sm_LPRM(i)>0.) then                     ! any neg is nodata
          
          j=j+1
          
          ae_sm_LPRM(j) = ae_sm_LPRM(i)
          lon(j)       = lon(i)
          lat(j)       = lat(i)
          if (present(ease_col))  ease_col(j) = ease_col(i)
          if (present(ease_row))  ease_row(j) = ease_row(i)
          
       end if
       
    end do
    
    N_data = j
    
  end subroutine read_ae_sm_LPRM_bin
  
  ! *****************************************************************
  
  subroutine read_obs_ae_sm_LPRM(                                 &
       work_path, exp_id,                                        &
       date_time, dtstep_assim, N_catd, tile_coord,              &
       tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
       this_obs_param,                                           &
       found_obs, ae_sm_LPRM, std_ae_sm_LPRM )
    
    ! Read observations of surface soil moisture from AMSR-E sm LPRM files
    ! (Richard de Jeu, Vrije Universiteit Amsterdam)
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation.
    !
    ! If there are N > 1 observations in a given tile,
    ! a "super-observation" is computed by averaging the N observations
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! reichle, 20 Feb 2009
    !
    ! --------------------------------------------------------------------
        
    implicit none
    
    ! inputs:
    
    character(*), intent(in) :: work_path
    character(*),  intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: dtstep_assim, N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: ae_sm_LPRM
    real,    intent(out), dimension(N_catd) :: std_ae_sm_LPRM
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals
    
    ! AMSR_E_L2_Land files are available for approximately 50min periods
    ! covering one ascending or descending swath.  The filename indicates
    ! the start time of the swath
    ! "ae_time_offset" is used to find the mean time of the the interval
    ! which is approximately the time of the equator overpass.
    ! This time is assigned to all observations of the swath.
    
    integer, parameter :: ae_time_offset = 1500   ! 25 minutes in seconds
    
    character(4)   :: DDHH
    character(6)   :: YYYYMM
    character(8)   :: date_string
    character(10)  :: time_string
    character(300) :: tmpfname, tmpfname2
    character(400) :: cmd

    type(date_time_type) :: date_time_tmp
    
    integer :: i, ind, N_tmp, N_files

    character(300), dimension(:), allocatable :: fnames
    
    real,    dimension(:), pointer :: tmp_obs, tmp_lat, tmp_lon
    integer, dimension(:), pointer :: tmp_tile_num
    
    integer, dimension(N_catd) :: N_obs_in_tile    

    character(len=*), parameter :: Iam = 'read_obs_ae_sm_LPRM'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------------
    
    nullify( tmp_obs, tmp_lat, tmp_lon, tmp_tile_num )    
    
    ! ---------------
    
    ! initialize
    
    found_obs = .false.
    
    ! find files that are within half-open interval 
    ! [date_time-dtstep_assim/2,date_time+dtstep_assim/2)
    
    date_time_tmp = date_time    
    
    call augment_date_time( -(dtstep_assim/2 + ae_time_offset), date_time_tmp )
    
    ! get tmp file name and remove file if it exists
    
    call date_and_time(date_string, time_string)  ! f90 intrinsic function
    
    tmpfname = trim(work_path) // '/' // 'tmp.' // trim(exp_id) &
         // '.' // date_string // time_string

    cmd = '/bin/rm -f ' // tmpfname 
    
    call Execute_command_line(trim(cmd))
    
    ! identify all files within current assimilation interval
    ! (list all files within hourly intervals)
    
    do i=1,(dtstep_assim/3600)
       
       write (YYYYMM,'(i6.6)') date_time_tmp%year*100 + date_time_tmp%month
       write (DDHH,  '(i4.4)') date_time_tmp%day *100 + date_time_tmp%hour
       
       cmd = 'ls ' // trim(this_obs_param%path) // '/' // YYYYMM(1:4) // &
            '.' // YYYYMM(5:6) // '/' // trim(this_obs_param%name) &
            // YYYYMM // DDHH(1:2) // '.' // DDHH(3:4) // '??' 
       
       if     (this_obs_param%descr(1:12)=='ae_sm_LPRM_a') then
          
          cmd = trim(cmd) // '_A.bin'
          
       elseif (this_obs_param%descr(1:12)=='ae_sm_LPRM_d') then
          
          cmd = trim(cmd) // '_D.bin'
          
       else
          
          !write (logunit,*) 'read_obs_ae_sm_LPRM(): unknown descr, STOPPING.'
          !stop
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown descr')
          
       end if
       
       cmd = trim(cmd) // ' >> ' // trim(tmpfname)
       
       call Execute_command_line(trim(cmd))
       
       call augment_date_time( 3600, date_time_tmp )
       
    end do
    
    ! find out how many need to be read
    
    tmpfname2 = trim(tmpfname) // '.wc'
    
    cmd = 'wc -w ' // trim(tmpfname) // ' > ' // trim(tmpfname2)
    
    call Execute_command_line(trim(cmd))
    
    open(10, file=tmpfname2, form='formatted', action='read')
    
    read(10,*) N_files
    
    close(10,status='delete')
    
    ! load file names into "fnames"
    
    open(10, file=tmpfname,  form='formatted', action='read')
    
    if (N_files>0) then
       
       allocate(fnames(N_files))
       
       do i=1,N_files
          read(10,'(a)') fnames(i)
       end do
       
    end if
    
    close(10,status='delete')
    
    ! read observations:
    !
    ! 1.) read N_tmp observations and their lat/lon info from file
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon
    ! 3.) compute super-obs for each tile from all obs w/in that tile
    !
    ! ----------------------------------------------------------------
    !
    ! 1.) read N_tmp observations and their lat/lon info from file

    if (N_files>0) then
              
       call read_ae_sm_LPRM_bin( &
            this_obs_param, N_files, fnames, &
            N_tmp, tmp_lon, tmp_lat, tmp_obs )
       
       if (logit) then
          
          write (logunit,*) 'read_obs_ae_sm_LPRM: read ', N_tmp,  &
               ' at date_time = ', date_time, ' from '
          do i=1,N_files
             write (logunit,*) trim(fnames(i))
          end do
          write (logunit,*) '----------'

       end if

       deallocate(fnames)
       
    else
       
       N_tmp = 0
       
    end if
    
    ! ------------------------------------------------------------------

    ! QC is done in matlab pre-processing to bin files
    
    ! ------------------------------------------------------------------
    !
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon

    if (N_tmp>0) then
       
       allocate(tmp_tile_num(N_tmp))
       
       call get_tile_num_for_obs(N_catd, tile_coord,                 &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
            N_tmp, tmp_lat, tmp_lon,                                 &
            this_obs_param,                                          &
            tmp_tile_num )

       
       ! ----------------------------------------------------------------
       !
       ! 3.) compute super-obs for each tile from all obs w/in that tile
       !     (also eliminate observations that are not in domain)
       
       ae_sm_LPRM = 0.
       N_obs_in_tile  = 0
       
       do i=1,N_tmp
          
          ind = tmp_tile_num(i)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
          
          if (ind>0) then         ! this step eliminates obs outside domain
             
             ae_sm_LPRM(ind) = ae_sm_LPRM(ind) + tmp_obs(i)
             
             N_obs_in_tile(ind) = N_obs_in_tile(ind) + 1
             
          end if
          
       end do
       
       ! normalize
       
       do i=1,N_catd
          
          if (N_obs_in_tile(i)>1) then
             
             ae_sm_LPRM(i) = ae_sm_LPRM(i)/real(N_obs_in_tile(i))
          
          elseif (N_obs_in_tile(i)==0) then
             
             ae_sm_LPRM(i) = this_obs_param%nodata
             
          end if
          
       end do
       
       ! clean up
       
       if (associated(tmp_tile_num)) deallocate(tmp_tile_num)
       
       ! --------------------------------
       
       ! set observation error standard deviation
       
       do i=1,N_catd
          std_ae_sm_LPRM(i) = this_obs_param%errstd
       end do
       
       ! --------------------------------
       
       if (any(N_obs_in_tile>0)) then
          
          found_obs = .true.
          
       else 
          
          found_obs = .false.
          
       end if
       
    end if
    
    ! clean up
    
    if (associated(tmp_obs))      deallocate(tmp_obs)
    if (associated(tmp_lon))      deallocate(tmp_lon)
    if (associated(tmp_lat))      deallocate(tmp_lat)

  end subroutine read_obs_ae_sm_LPRM

  ! *****************************************************************

  subroutine read_obs_sm_ASCAT(                                  &
       work_path, exp_id,                                        &
       date_time, dtstep_assim, N_catd, tile_coord,              &
       tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
       this_obs_param,                                           &
       found_obs, sm_ASCAT, std_sm_ASCAT )
    
    !---------------------------------------------------------------------
    ! 
    ! Routine to read in ASCAT surface degree of saturation obs. 
    ! Output is found_obs, sm_ASCAT, std_sm_ASCAT 
    !
    ! Reads in obs provided by Wolfgang Wagner (TUW), converted to
    ! once hourly binary files (and projected onto EASE grid). 
    ! Updating to the EUMETSAT BUFR (DGG) files will require a new
    ! reader, and changes to file-name / time-stamping
    !
    !  Draper, May 2011.  
    !  Based on read_obs_ae_sm_LPRM
    !
    ! --------------------------------------------------------------------
        
    implicit none
    
    ! inputs:
    
    character(*), intent(in) :: work_path
    character(*),  intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: dtstep_assim, N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: sm_ASCAT
    real,    intent(out), dimension(N_catd) :: std_sm_ASCAT
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! For the files generated from the TUW timseries, each 
    ! file is time-stamped with the time of the observations, rounded
    ! down to the nearest obs 
    ! E(minutes past hour)=30 -> use this as the offset

    ! Will need to be updated if using EUMETSAT BUFR files

    integer, parameter :: ae_time_offset = 1800   ! 30 minutes in seconds
    
    character(4)   :: DDHH
    character(6)   :: YYYYMM
    character(8)   :: date_string
    character(10)  :: time_string
    character(300) :: tmpfname, tmpfname2
    character(400) :: cmd

    type(date_time_type) :: date_time_tmp
    
    integer :: i, ind, N_tmp, N_files

    character(300), dimension(:), allocatable :: fnames
    
    real,    dimension(:), pointer :: tmp_obs, tmp_lat, tmp_lon
    integer, dimension(:), pointer :: tmp_tile_num
    
    integer, dimension(N_catd) :: N_obs_in_tile    

    real, parameter :: tol = 1e-2

    character(len=*), parameter :: Iam = 'read_obs_sm_ASCAT'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------------
    
    nullify( tmp_obs, tmp_lat, tmp_lon, tmp_tile_num )    
    
    ! ---------------
    
    ! initialize
    
    found_obs = .false.

    ! find files that are within half-open interval 
    ! [date_time-dtstep_assim/2,date_time+dtstep_assim/2)
    
    date_time_tmp = date_time    
    
    call augment_date_time( -(dtstep_assim/2 + ae_time_offset), date_time_tmp )
    
    ! get tmp file name and remove file if it exists
    
    call date_and_time(date_string, time_string)  ! f90 intrinsic function
    
    tmpfname = trim(work_path) // '/' // 'tmp.' // trim(exp_id) &
         // '.' // date_string // time_string

    cmd = '/bin/rm -f ' // tmpfname 
    
    call Execute_command_line(trim(cmd))
    
    ! identify all files within current assimilation interval
    ! (list all files within hourly intervals)
    
    do i=1,(dtstep_assim/3600)
       
       write (YYYYMM,'(i6.6)') date_time_tmp%year*100 + date_time_tmp%month
       write (DDHH,  '(i4.4)') date_time_tmp%day *100 + date_time_tmp%hour
       
       ! TUW files time stamped with hour only. Update for EUMETSAT BUFR
       cmd = 'ls ' // trim(this_obs_param%path) // '/' // YYYYMM(1:4) // &
            '.' // YYYYMM(5:6) // '/' // trim(this_obs_param%name) &
            // YYYYMM // DDHH(1:2) // '.' // DDHH(3:4) !// '??' 
       
       if     (this_obs_param%descr(1:10)=='ASCAT_SM_A') then
          
          cmd = trim(cmd) // '_A.bin'
          
       elseif (this_obs_param%descr(1:10)=='ASCAT_SM_D') then
          
          cmd = trim(cmd) // '_D.bin'
          
       else
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown descr')
          
       end if

       cmd = trim(cmd) // ' >> ' // trim(tmpfname)
       
       call Execute_command_line(trim(cmd))
       
       
       call augment_date_time( 3600, date_time_tmp )
       
    end do
    
    ! find out how many need to be read
    
    tmpfname2 = trim(tmpfname) // '.wc'
    
    cmd = 'wc -w ' // trim(tmpfname) // ' > ' // trim(tmpfname2)
    
    call Execute_command_line(trim(cmd))
    
    open(10, file=tmpfname2, form='formatted', action='read')
    
    read(10,*) N_files
    
    close(10,status='delete')
    
    ! load file names into "fnames"
    
    open(10, file=tmpfname,  form='formatted', action='read')
    
    if (N_files>0) then
       
       allocate(fnames(N_files))
       
       do i=1,N_files
          read(10,'(a)') fnames(i)
       end do
       
    end if
    
    close(10,status='delete')
    
    ! read observations:
    !
    ! 1.) read N_tmp observations and their lat/lon info from file
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon
    ! 3.) compute super-obs for each tile from all obs w/in that tile
    !
    ! ----------------------------------------------------------------
    !
    ! 1.) read N_tmp observations and their lat/lon info from file

    if (N_files>0) then
              
       call read_sm_ASCAT_bin( &
            this_obs_param, N_files, fnames, &
            N_tmp, tmp_lon, tmp_lat, tmp_obs )
       
       if (logit) then
          
          write (logunit,*) 'read_obs_sm_ASCAT: read ', N_tmp,  &
               ' at date_time = ', date_time, ' from '
          do i=1,N_files
             write (logunit,*) trim(fnames(i))
          end do
          write (logunit,*) '----------'

       end if

       deallocate(fnames)
       
    else
       
       N_tmp = 0
       
    end if
    
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! SOME QC SHOULD BE DONE HERE!!!
    !
    ! MAKE SURE no-data-values ARE DEALT WITH
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! ----------------------------------------------------------------
    !
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon

    if (N_tmp>0) then
       
       allocate(tmp_tile_num(N_tmp))
       
       call get_tile_num_for_obs(N_catd, tile_coord,                 &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
            N_tmp, tmp_lat, tmp_lon,                                 &
            this_obs_param,                                          &
            tmp_tile_num )
 
       
       ! ----------------------------------------------------------------
       !
       ! 3.) compute super-obs for each tile from all obs w/in that tile
       !     (also eliminate observations that are not in domain)
       
       sm_ASCAT = 0.
       N_obs_in_tile  = 0
       
       do i=1,N_tmp
          
          ind = tmp_tile_num(i)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
          
          if (ind>0) then         ! this step eliminates obs outside domain
             
             sm_ASCAT(ind) = sm_ASCAT(ind) + tmp_obs(i)
             
             N_obs_in_tile(ind) = N_obs_in_tile(ind) + 1
             
          end if
          
       end do
       
       ! normalize
       
       do i=1,N_catd
          
          if (N_obs_in_tile(i)>1) then
             
             sm_ASCAT(i) = sm_ASCAT(i)/real(N_obs_in_tile(i))
          
          elseif (N_obs_in_tile(i)==0) then
             
             sm_ASCAT(i) = this_obs_param%nodata
             
          end if
          
       end do
       
       ! clean up
       
       if (associated(tmp_tile_num)) deallocate(tmp_tile_num)
       
       ! --------------------------------
       
       ! set observation error standard deviation

         do i=1,N_catd
	      std_sm_ASCAT(i) = this_obs_param%errstd
         enddo
       ! --------------------------------
       
       if (any(N_obs_in_tile>0)) then
        
          found_obs = .true.
          
       else 
          
          found_obs = .false.
          
       end if
       
    end if
    
    ! clean up
    
    if (associated(tmp_obs))      deallocate(tmp_obs)
    if (associated(tmp_lon))      deallocate(tmp_lon)
    if (associated(tmp_lat))      deallocate(tmp_lat)

  end subroutine read_obs_sm_ASCAT
  
  ! ****************************************************************************

  subroutine read_obs_sm_ASCAT_EUMET(                                  &
       work_path, exp_id,                                        &
       date_time, dtstep_assim, N_catd, tile_coord,              &
       tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
       this_obs_param,                                           &
       found_obs, sm_ASCAT, std_sm_ASCAT )
    
    !---------------------------------------------------------------------
    ! 
    ! Routine to read in ASCAT surface degree of saturation obs. 
    ! Output is found_obs, sm_ASCAT, std_sm_ASCAT 
    !
    ! Read in the EUMETSAT level 2 soil mositure product 25 km (SMO), PPF software version 5.0 
    ! the data correspond to re-sampled (spatially averaged) sigma0 values, on a 25 km
    ! orbit swath grid. The input data files are in BUFR file format.
    !
    !  Q. Liu, Nov. 2019.  
    ! based on read_obs_sm_ASCAT
    ! --------------------------------------------------------------------
        
    implicit none
    
    ! inputs:
    
    character(*), intent(in) :: work_path
    character(*),  intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: dtstep_assim, N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: sm_ASCAT   ! wetness range 0-1
    real,    intent(out), dimension(N_catd) :: std_sm_ASCAT
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! Each obs file contains about 1 hour 40 minutes observations
    ! file name indicates the start time of the swaths.  
    ! "ae_time_offset" is used to find the mean time of the the interval
    ! which is approximately the time of the equator overpass.
    ! This time is assigned to all observations of the swath.

    ! Will need to be updated if using EUMETSAT BUFR files

    integer, parameter :: ae_time_offset = 3600   ! 60 minutes in seconds
    
    character(4)   :: DDHH
    character(6)   :: YYYYMM
    character(8)   :: date_string
    character(10)  :: time_string
    character(300) :: tmpfname, tmpfname2
    character(400) :: cmd

    type(date_time_type) :: date_time_tmp
    type(date_time_type) :: date_time_low
    type(date_time_type) :: date_time_upp
    
    integer :: i, ind, N_tmp, N_files

    character(300), dimension(:), allocatable :: fnames
 
    real(8) :: tmp_data, tmp_vdata(4), tmp_time(6) 
    integer, parameter :: lnbufr = 50
    integer, parameter :: max_rec = 200000 
    integer :: idate,iret,kk
    integer :: ireadmg,ireadsb
    character(8)  :: subset
    real,  dimension(:),     allocatable  :: tmp1_lon, tmp1_lat, tmp1_obs

    real,    dimension(:), pointer :: tmp_obs, tmp_lat, tmp_lon
    integer, dimension(:), pointer :: tmp_tile_num

    integer, dimension(N_catd) :: N_obs_in_tile    

    real, parameter :: tol = 1e-2

    character(len=*), parameter :: Iam = 'read_obs_sm_ASCAT_EUMET'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------------
    
    nullify( tmp_obs, tmp_lat, tmp_lon, tmp_tile_num )    
    
    ! ---------------
    
    ! initialize
    
    found_obs = .false.

    ! find files that are within half-open interval 
    ! [date_time-dtstep_assim/2,date_time+dtstep_assim/2)

    date_time_low = date_time    
    call augment_date_time( -(dtstep_assim/2), date_time_low)
    date_time_upp = date_time    
    call augment_date_time(  (dtstep_assim/2), date_time_upp)

    ! for ASCAT file name stamp
    date_time_tmp = date_time 
    call augment_date_time( -(dtstep_assim/2 + ae_time_offset), date_time_tmp )
    
    ! get tmp file name and remove file if it exists
    
    call date_and_time(date_string, time_string)  ! f90 intrinsic function
    
    tmpfname = trim(work_path) // '/' // 'tmp.' // trim(exp_id) &
         // '.' // date_string // time_string

    cmd = '/bin/rm -f ' // tmpfname 
    
    call Execute_command_line(trim(cmd))
    
    ! identify all files within current assimilation interval
    ! (list all files within hourly intervals)
   
    ! Every EUMETSTA BUFR contains data over ~2 hr sensing period. it's necessary to
    ! search additional files for obs. 
    do i=1,(dtstep_assim/3600)+2
       
       write (YYYYMM,'(i6.6)') date_time_tmp%year*100 + date_time_tmp%month
       write (DDHH,  '(i4.4)') date_time_tmp%day *100 + date_time_tmp%hour
       
       ! EUMETSAT BUFR
       cmd = 'ls ' // trim(this_obs_param%path) // '/Y' // YYYYMM(1:4) // &
            '/M' // YYYYMM(5:6) // '/' // trim(this_obs_param%name) // '*-'&
            // YYYYMM // DDHH // '*Z-*.bfr' 
       
       cmd = trim(cmd) // ' >> ' // trim(tmpfname)
       
       call Execute_command_line(trim(cmd))
       
       call augment_date_time( 3600, date_time_tmp )
       
    end do
    
    ! find out how many need to be read
    
    tmpfname2 = trim(tmpfname) // '.wc'
    
    cmd = 'wc -w ' // trim(tmpfname) // ' > ' // trim(tmpfname2)
    
    call Execute_command_line(trim(cmd))
    
    open(10, file=tmpfname2, form='formatted', action='read')
    
    read(10,*) N_files
    
    close(10,status='delete')
    
    ! load file names into "fnames"
    
    open(10, file=tmpfname,  form='formatted', action='read')
    
    if (N_files>0) then
       
       allocate(fnames(N_files))
       
       do i=1,N_files
          read(10,'(a)') fnames(i)
          write(logunit,*) trim(fnames(i))
       end do
       
    end if
    
    close(10,status='delete')
    
    ! read observations:
    !
    ! 1.) read N_tmp observations and their lat/lon info from file
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon
    ! 3.) compute super-obs for each tile from all obs w/in that tile
    !
    ! ----------------------------------------------------------------
    !
    ! 1.) read N_tmp observations and their lat/lon info from file

    ! read and process data if files are found
    allocate(tmp1_lon(max_rec))
    allocate(tmp1_lat(max_rec))
    allocate(tmp1_obs(max_rec))

    if (N_files>0) then
       
       ! file loop
       N_tmp = 0
       do kk = 1,N_files

          ! open on bufr file
          call closbf(lnbufr)
          open(lnbufr, file=trim(fnames(kk)), action='read',form='unformatted')
          call openbf(lnbufr,'SEC3', lnbufr)
          call MTINFO( trim(this_obs_param%path) // '/BUFR_mastertable/', 51, 52)
          call datelen(10)
         
          msg_report: do while(ireadmg(lnbufr,subset,idate) ==0)  
            loop_report: do while(ireadsb(lnbufr) == 0)
            ! extract sensing time information 
                 call ufbint(lnbufr,tmp_time,6,1,iret,'YEAR MNTH DAYS HOUR MINU SECO')
                 date_time_tmp.year = int(tmp_time(1))
                 date_time_tmp.month = int(tmp_time(2)) 
                 date_time_tmp.day = int(tmp_time(3))
                 date_time_tmp.hour = int(tmp_time(4))
                 date_time_tmp.min = int(tmp_time(5))
                 date_time_tmp.sec = int(tmp_time(6))
                 ! skip if record outside of current assim window
                 if ( datetime_lt_refdatetime( date_time_low, date_time_tmp ) .and.        &
                      datetime_le_refdatetime( date_time_tmp, date_time_upp )) cycle loop_report
                 
                 ! skip if record contain no valid soil moisture value
                 call ufbint(lnbufr,tmp_data,1,1,iret,'SSOM')
                 if(tmp_data > 100. .or. tmp_data < 0.) cycle loop_report

                 ! EUMETSAT file contains data of both ascending and descending orbits. 
                 ! DOMO - “Direction of motion of moving observing platform” is used to seperate Asc and Desc
                 ! because the file doesn't contain any explicit orbit indicator variable.
                 ! according to Pamela Schoebel-Pattiselanno, EUMETSAT User Services Helpdesk 
                 ! "When the value (of DOMO) is between 180 and 270 degrees, it is the descending part 
                 ! of the orbit … when it is between 270 and 360 degrees, it is the ascending part"
                 call ufbint(lnbufr,tmp_data,1,1,iret,'DOMO')
                 if (index(this_obs_param%descr,'_A') /=0 .and. (tmp_data < 270 .or. tmp_data > 360)) cycle loop_report
                 if (index(this_obs_param%descr,'_D') /=0 .and. (tmp_data < 180 .or. tmp_data >= 270)) cycle loop_report
                  
                 ! skip if processing flag is set               
                 call ufbint(lnbufr,tmp_data,1,1,iret,'SMPF')
                 if(int(tmp_data) /= 0) cycle loop_report

                 ! skip if correction flag is set               
                 call ufbint(lnbufr,tmp_data,1,1,iret,'SMCF')
                 ! if (.not. (int(tmp_data) == 0 .or. int(tmp_data) == 4)) cycle loop_report
                 if(int(tmp_data) /= 0) cycle loop_report

                 ! skip if land fraction is missing or < 0.9
                 call ufbint(lnbufr,tmp_data,1,1,iret,'ALFR')
                 if(tmp_data >1 .or. tmp_data < 0.9 ) cycle loop_report

                 ! additioanal QC varibles from file               
                 ! skip if topographic complexity > 10%
                 call ufbint(lnbufr,tmp_data,1,1,iret,'TPCX') ! topo complexity
                 if(tmp_data > 10.) cycle loop_report
                 
                 ! skip if inudatation and wetland faction > 10%
                 call ufbint(lnbufr,tmp_data,1,1,iret,'IWFR') ! Inundation And Wetland Fraction
                 if(tmp_data > 10.) cycle loop_report
                 !call ufbint(lnbufr,tmp_data,1,1,iret,'SNOC') ! snow cover
                 !call ufbint(lnbufr,tmp_data,1,1,iret,'FLSF') ! frozen land fraction

                 N_tmp = N_tmp + 1
                 call ufbint(lnbufr,tmp_vdata,4,1,iret,'CLATH CLONH SSOM EESSM')
                 tmp1_lat(N_tmp) = tmp_vdata(1) 
                 tmp1_lon(N_tmp) = tmp_vdata(2) 
                 tmp1_obs(N_tmp) = tmp_vdata(3)/100.  ! change value from 0-100 to 0-1
                 !tmp_obserr(N_tmp) = tmp_vdata(4)

            end do loop_report

          end do msg_report
          call closbf(lnbufr)
          close(lnbufr)

       end do ! end file loop

       if (logit) then
          
          write (logunit,*) 'read_obs_sm_ASCAT_EUMET: read ', N_tmp,  &
               ' at date_time = ', date_time, ' from '
          do i=1,N_files
             write (logunit,*) trim(fnames(i))
          end do
          write (logunit,*) '----------'
          write (logunit,*) 'max(obs)=',maxval(tmp1_obs(1:N_tmp)), 'min(obs)=',minval(tmp1_obs(1:N_tmp)), &
               ' avg(obs)=',sum(tmp1_obs(1:N_tmp))/N_tmp
       end if

       deallocate(fnames)
    else
       N_tmp = 0
       
    end if

    allocate(tmp_lon(N_tmp))
    allocate(tmp_lat(N_tmp))
    allocate(tmp_obs(N_tmp))

    tmp_lon = tmp1_lon(1:N_tmp)
    tmp_lat = tmp1_lat(1:N_tmp)
    tmp_obs = tmp1_obs(1:N_tmp)
   
    deallocate(tmp1_lon)
    deallocate(tmp1_lat)
    deallocate(tmp1_obs) 

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! SOME QC SHOULD BE DONE HERE!!!
    !
    ! MAKE SURE no-data-values ARE DEALT WITH
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! ----------------------------------------------------------------
    !
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon

    if (N_tmp>0) then
       
       allocate(tmp_tile_num(N_tmp))
       
       call get_tile_num_for_obs(N_catd, tile_coord,                 &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
            N_tmp, tmp_lat, tmp_lon,                                 &
            this_obs_param,                                          &
            tmp_tile_num )
 
       
       ! ----------------------------------------------------------------
       !
       ! 3.) compute super-obs for each tile from all obs w/in that tile
       !     (also eliminate observations that are not in domain)
       
       sm_ASCAT = 0.
       N_obs_in_tile  = 0
       
       do i=1,N_tmp
          
          ind = tmp_tile_num(i)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
          
          if (ind>0) then         ! this step eliminates obs outside domain
             
             sm_ASCAT(ind) = sm_ASCAT(ind) + tmp_obs(i)
             
             N_obs_in_tile(ind) = N_obs_in_tile(ind) + 1
             
          end if
          
       end do
       
       ! normalize
       
       do i=1,N_catd
          
          if (N_obs_in_tile(i)>1) then
             
             sm_ASCAT(i) = sm_ASCAT(i)/real(N_obs_in_tile(i))
          
          elseif (N_obs_in_tile(i)==0) then
             
             sm_ASCAT(i) = this_obs_param%nodata
             
          end if
          
       end do
       
       ! clean up
       
       if (associated(tmp_tile_num)) deallocate(tmp_tile_num)
       
       ! --------------------------------
       
       ! set observation error standard deviation

         do i=1,N_catd
	      std_sm_ASCAT(i) = this_obs_param%errstd
         enddo
       ! --------------------------------
       
       if (any(N_obs_in_tile>0)) then
        
          found_obs = .true.
          
       else 
          
          found_obs = .false.
          
       end if
       
    end if
    
    ! clean up
    
    if (associated(tmp_obs))      deallocate(tmp_obs)
    if (associated(tmp_lon))      deallocate(tmp_lon)
    if (associated(tmp_lat))      deallocate(tmp_lat)

  end subroutine read_obs_sm_ASCAT_EUMET

  ! ***************************************************************************

  subroutine read_sm_ASCAT_bin( &
       this_obs_param, N_files, fnames, N_data, lon, lat, sm_ASCAT, ease_col, ease_row )
    
    ! read soil moisture data from one or more ASCAT bin files
    !
    ! return ONLY valid data points (ie. excluding no-data-values)
    ! 
    ! no QC in addition to what was done in matlab-preprocessing
    !
    ! DRAPER, May 2011
    ! based on read_ae_sm_LPRM_bin
    !
    ! DRAPER, July 2012
    ! updated, for inclusion of SDS  error in ASCAT file 
    ! (error info currently not saved)

    implicit none
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    integer, intent(in) :: N_files
    
    character(*), dimension(N_files), intent(in) :: fnames
    
    integer, intent(out) :: N_data
    
    real, dimension(:), pointer :: lon, lat, sm_ASCAT  ! output
    
    integer, dimension(:), pointer, optional :: ease_col, ease_row ! output
    
    ! local variables
    
    logical :: must_stop
    
    integer, dimension(N_files) :: N_data_tmp
    
    integer :: i, j, k_off
    
    integer, dimension(:), allocatable :: tmpintvec
    real,    dimension(:), allocatable :: tmprealvec
    
    character(len=*), parameter :: Iam = 'read_sm_ASCAT_bin'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------
    
    ! make sure pointers are not allocated or associated
    
    must_stop = .false.
    
    if ( associated(lon) .or. associated(lat) .or. associated(sm_ASCAT) ) then
       must_stop = .true.
    end if
    
    if ( present(ease_col) ) then         
       if (associated(ease_col))  must_stop = .true.
    end if
    
    if ( present(ease_row) ) then
       if (associated(ease_row))  must_stop = .true.
    end if
    
    if (must_stop) then
       err_msg = 'output pointers must not be ' // &
            'associated or allocated on input.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    
    !write (logunit,*) 'read_ae_sm_LPRM_bin(): output pointers must not be ' // &
    !     'associated or allocated on input. STOPPING.'
    !stop
    !
    !end if
    
    ! determine number of data to be read from each file
    
    N_data = 0
    
    do j=1,N_files
       
       ! open file 
       
       open( 10, file=trim(fnames(j)), form='unformatted',convert='big_endian', action='read' )
       
       read( 10) N_data_tmp(j)
       
       close(10, status='keep')
       
    end do
    
    ! allocate pointers (must be deallocated outside this subroutine!)
    
    N_data = sum(N_data_tmp)
    
    allocate(lon(     N_data))
    allocate(lat(     N_data))
    allocate(sm_ASCAT(N_data))
    
    if (present(ease_col))  allocate(ease_col(N_data))
    if (present(ease_row))  allocate(ease_row(N_data))
    
    ! read data into arrays, concatenate data from N_files files
    
    ! format of AMSR_sm_LPRM_EASE_bin files:
    !
    ! record  1 -- N_data            int*4
    ! record  2 -- lon(  1:N_data)   real*4    
    ! record  3 -- lat(  1:N_data)   real*4    
    ! record  4 -- sds(  1:N_data)   real*4    
    ! record  5 -- sds_noise (1:N_data) real*4 
    ! record  6 -- ind_i(1:N_data)   int*4    zero-based (!) EASE row index
    ! record  7 -- ind_j(1:N_data)   int*4    zero-based (!) EASE col index
        
    k_off = 0
    
    do j=1,N_files
       
       allocate(tmprealvec(N_data_tmp(j)))
       
       if (present(ease_col))  allocate(tmpintvec(N_data_tmp(j)))
       
       open (10, file=trim(fnames(j)),  form='unformatted', convert='big_endian',  action='read' )
       
       ! re-read N_data
       
       read (10) N_data_tmp(j)       
       
       ! read data as needed 
       
       read (10) tmprealvec; lon(k_off+1:k_off+N_data_tmp(j)) = tmprealvec    
       read (10) tmprealvec; lat(k_off+1:k_off+N_data_tmp(j)) = tmprealvec    
       
       read (10) tmprealvec; sm_ASCAT(k_off+1:k_off+N_data_tmp(j)) = tmprealvec
       
       ! skip record with error 
       
       read (10) !  tmprealvec; er_ASCAT(k_off+1:k_off+N_data_tmp(j)) = tmprealvec
       
       if (present(ease_col) .and. present(ease_row)) then
          
          read (10) tmpintvec; ease_col(k_off+1:k_off+N_data_tmp(j)) = tmpintvec
          read (10) tmpintvec; ease_row(k_off+1:k_off+N_data_tmp(j)) = tmpintvec
          
       end if
       
       ! clean up
       
       close(10, status='keep')
       
       deallocate(tmprealvec)
       if (allocated(tmpintvec))  deallocate(tmpintvec)
       
       ! prepare next j
       
       k_off = k_off + N_data_tmp(j)
       
    end do
    
    ! -------------------------------------
    !
    ! eliminate no-data-values 
    
    j = 0
    
    do i=1,N_data
       
       if (sm_ASCAT(i)>0.) then                     ! any neg is nodata
          
          j=j+1
          
          sm_ASCAT(j) = sm_ASCAT(i)
          lon(j)       = lon(i)
          lat(j)       = lat(i)
          if (present(ease_col))  ease_col(j) = ease_col(i)
          if (present(ease_row))  ease_row(j) = ease_row(i)
          
       end if
       
    end do
    
    N_data = j
    
  end subroutine read_sm_ASCAT_bin


  ! *****************************************************************
  
  subroutine read_obs_LaRC_Tskin(                                &
       work_path, exp_id,                                        &
       date_time, dtstep_assim, N_catd, tile_coord,              &
       tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
       this_obs_param,                                           &
       found_obs, ts_LARC, std_ts_LARC )

    !---------------------------------------------------------------------
    ! 
    !  Subroutine to read in Tskin from LaRC. 
    !  Draper, June 2012.  
    !
    ! --------------------------------------------------------------------
        
    implicit none
    
    ! inputs:
    
    character(*), intent(in) :: work_path
    character(*),  intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: dtstep_assim, N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord ! input
    
    type(grid_def_type), intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij

    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: ts_LARC
    real,    intent(out), dimension(N_catd) :: std_ts_LARC
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    real, parameter :: min_Tskin = 200.
    real, parameter :: max_Tskin = 370.
    
    real, parameter :: tol = 1.e-2
    
    ! Offsets taken 
    
    integer        :: ts_time_offset   ! in seconds
    
    character(6)   :: DDHHMM
    character(6)   :: YYYYMM
    
    type(date_time_type) :: date_time_tmp
    
    integer :: i, ind, N_tmp, N_files
    
    ! 24 = max number of files in one cycle
    ! (assumes daily assim cycle, and hourly files) 
    
    character(300), dimension(24)              :: fnames
    
    real,           dimension(:),      pointer :: tmp_obs, tmp_lat, tmp_lon
    
    integer,        dimension(N_catd)          :: N_obs_in_tile
    
    integer,        dimension(:),      pointer :: tmp_tile_num
    
    character(300)  :: tmp_fname
        
    logical         :: ex
    
    integer         :: MM 
    
    character(len=*), parameter :: Iam = 'read_obs_LaRC_Tskin'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------------
    
    nullify( tmp_obs, tmp_lat, tmp_lon , tmp_tile_num)    
    
    ! ---------------
    
    ! initialize
        
    ! time stampes are at start of scan. MET-09 scans much faster than others
    
    select case (trim(this_obs_param%descr))
       
    case ('LaRC_tskin-GOESE','LaRC_tskin-GOESW', 'LaRC_tskin-FY2E-') 
       
       ts_time_offset=26*60/2 

    case ('LaRC_tskin-MTST2')
       
       ts_time_offset=28*60/2 
       
    case ('LaRC_tskin-MET09') 
       
       ts_time_offset=12*60/2
       
    case default     
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown obs_param%descr')
       
    end select
    
    found_obs = .false.
    
    ! find files that are within half-open interval 
    ! [date_time-dtstep_assim/2,date_time+dtstep_assim/2)

    select case (trim(this_obs_param%descr))

    case ('LaRC_tskin-GOESW')
       MM=00
    case ('LaRC_tskin-GOESE') 
       MM=45
    case ('LaRC_tskin-MET09') 
       MM=00
    case ('LaRC_tskin-FY2E-')
       MM=00
    case ('LaRC_tskin-MTST2')
       MM=30
    case default     
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown obs_param%descr')
    end select
    
    date_time_tmp = date_time    
    
    call augment_date_time( -(dtstep_assim/2 + ts_time_offset), date_time_tmp )
    
    ! identify all files within current assimilation interval
    ! inquire for the file once per hour - will always be at same minutes past hour
    
    N_files=0 
    
    do i=1,(dtstep_assim/3600)
       
       write (YYYYMM, '(i6.6)') date_time_tmp%year*100  + date_time_tmp%month
       write (DDHHMM, '(i6.6)') date_time_tmp%day*10000 + date_time_tmp%hour*100 + MM
       
       tmp_fname=trim(this_obs_param%path) // '/' // YYYYMM(1:4) // &
            '/' // YYYYMM(5:6) // '/' // DDHHMM(1:2) // '/' //  &
	    trim(this_obs_param%name) &
            // YYYYMM // DDHHMM(1:2) // '_' // DDHHMM(3:6) // 'z.nc4'
       
       inquire(file=trim(tmp_fname),exist=ex) 
       
       if (ex) then 
          
          N_files         = N_files + 1
          
          fnames(N_files) = tmp_fname
          
       end if
       
       call augment_date_time( 3600, date_time_tmp )
       
    end do
    
    ! read observations:
    !
    ! 1.) read N_tmp observations and their lat/lon info from file
    ! 2.) get tile number for each obs from lat/lon
    ! 3.) compute super-obs for each tile from all obs w/in that tile
    !     (also eliminate observations that are not in domain)
    ! ----------------------------------------------------------------
    !
    ! 1.) read N_tmp observations and their lat/lon info from file
    
    if (N_files>0) then
       
       call read_LaRC_Tskin_nc4(                          &
            this_obs_param, N_files, fnames(1:N_files),   &
            N_tmp, tmp_lon, tmp_lat, tmp_obs )
       
       if (logit) then
          
          write (logunit,*) 'read_obs_LaRC_Tskin: read ', N_tmp,  &
               ' at date_time = ', date_time
          do i=1,N_files
             write (logunit,*) trim(fnames(i))
          end do
          write (logunit,*) '----------'
          
       end if

    else
       
       N_tmp = 0
       
    end if

    ! 2.) get tile number for each obs from lat/lon

    if (N_tmp>0) then
       
       allocate(tmp_tile_num(N_tmp))
       
       call get_tile_num_for_obs(N_catd, tile_coord,                 &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
            N_tmp, tmp_lat, tmp_lon,                                 &
            this_obs_param,                                          &
            tmp_tile_num )
     

       ! ----------------------------------------------------------------
       !
       ! 3.) compute super-obs for each tile from all obs w/in that tile
       !     (also eliminate observations that are not in domain)
       
       ts_LARC = 0.

       N_obs_in_tile  = 0
       
       do i=1,N_tmp
          
          ind = tmp_tile_num(i)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
          
          if (ind>0) then         ! this step eliminates obs outside domain
             
             ! make sure obs is within allowed range
             
             if ( (min_Tskin<tmp_obs(i)) .and. (tmp_obs(i)<max_Tskin) ) then
                
                ts_LARC(ind) = ts_LARC(ind) + tmp_obs(i)
                
                N_obs_in_tile(ind) = N_obs_in_tile(ind) + 1
                
             end if

          end if
          
       end do
       
       ! normalize
       
       do i=1,N_catd
          
          if (N_obs_in_tile(i)>1) then
             
             ts_LARC(i) = ts_LARC(i)/real(N_obs_in_tile(i))
             
          elseif (N_obs_in_tile(i)==0) then
             
             ts_LARC(i) = this_obs_param%nodata
             
          end if
          
       end do
       
       if (associated(tmp_tile_num)) deallocate(tmp_tile_num)
       
       do i=1,N_catd
          
          std_ts_LARC(i) = this_obs_param%errstd
          
          ! CSD - temporary - only works over Americas
          ! if keep, base on SZA
          
          if ( (date_time_tmp%hour .GT. 2 ) .AND. (date_time_tmp%hour .LT. 14.5 )) then 

             std_ts_LARC(i) = 1.3 ! 1.5
          
          else 
          
             std_ts_LARC(i) = 2.1 ! 2.2
          
          end if
          
          ! CSD  - end temporary. 
          ! --------------------------------
          
       end do
       
       if (any(N_obs_in_tile>0)) then
          
          found_obs = .true.
          
       else 
          
          found_obs = .false.
          
       end if
       
    end if
    
    ! clean up
    
    if (associated(tmp_obs))      deallocate(tmp_obs)
    if (associated(tmp_lon))      deallocate(tmp_lon)
    if (associated(tmp_lat))      deallocate(tmp_lat)
    
  end subroutine read_obs_LaRC_Tskin

  ! ***************************************************************************

  subroutine read_LaRC_Tskin_nc4( &
       this_obs_param, N_files, fnames, N_data, lon, lat, ts_LaRC )
    
    ! read Tskin from LaRC nc4 files 
    !
    ! Apply QC:
    ! Screen longitude to retain obs only from closest GEOsat
    ! Viewing zenith angle (VZA) 
    ! Cloud fraction
    ! Sun zenith angle (SZA) 
    ! 
    ! returns number of data, and pointers to the lon, lat, and data
    ! return ONLY valid data points (ie. excluding no-data-values)
    !
    ! Currently tailored to read in LaRC files that have been reprocessed to replace
    ! integer data with float data (to enable GFIO to be used)
    ! 
    ! If LaRC fixes the files, will need to change: 
    ! if ( nvars .LT 5 ) - this is included as some LaRC files are missing HIRESTSKIN field
    ! Replace lower case variable names with appropriate case
    ! Replace nodata value (hardwired, as cannot read missing_value from file) 
    ! Replace time setting with actual minutes (currently rounded down to nearest hour,
    ! due to bug/assumption in GFIO read var routines). 

    ! DRAPER, June 2011
    
    implicit none
    
    type(obs_param_type), intent(in) :: this_obs_param

    integer, intent(in) :: N_files
    
    character(*), dimension(N_files), intent(in) :: fnames
    
    integer, intent(out) :: N_data
    
    real, dimension(:), pointer :: lon, lat, ts_LaRC  ! output
    
    ! local variables
    
    integer :: i, j, fid, rc
    integer :: x, x_min, x_max, y, y_min, y_max
    integer :: YYYYMMDD, HHMMSS
    integer :: g_nlon, g_nlat, km, lm, nvars, ngatts
    
    real :: nodata_LARC
    
    real, parameter :: tol=0.01

    ! QC
    
    real, parameter :: max_fcld      = 20. ! max total cloud fraction (%)
    real, parameter :: max_vza       = 60. ! max viewing zenith angle
    real, parameter :: min_excl_sza  = 82. ! min of sza exclusion interval
    real, parameter :: max_excl_sza  = 90. ! max of sza exclusion interval

    logical               :: tskin_ok, fcld_ok, sza_ok, vza_ok
    
    real,    dimension(2) :: lon_range, lat_range

    real                  :: dlat, dlon

    logical :: first 
    
    real,    dimension(:,:), allocatable :: tmp_data, ave_data, sum_data
    
    integer, dimension(:,:), allocatable :: cnt_data
    
    real,    dimension(:,:), allocatable :: tmp_vza, tmp_sza, tmp_fcld

    ! ----------------------
    !
    ! variables to fix viewing zenith angle bug in some of the 2012 GOES-East files 
    !  reprocessed by Ben Scarino
    !
    ! the "tmp_goesEast_vza_*" parameters point to a dummy "vza.nc4" file that contains 
    !  the correct vza values (located within the path of the reprocessed GOES-East files)
    
    character( 40), parameter :: tmp_goesEast_vza_dirname  = 'goesEast201205_201304_RM2'
    character( 40), parameter :: tmp_goesEast_vza_fname    = 'int2float/vza.nc4'
        
    integer,        parameter :: tmp_goesEast_vza_YYYYMMDD = 20120629   
    integer,        parameter :: tmp_goesEast_vza_HHMMSS   =   230000   
    
    character(300)            :: tmp_fname
    
    integer                   :: ind, tmp_fid

    character(len=*), parameter :: Iam = 'read_LaRC_Tskin_nc4'
    character(len=400) :: err_msg
    
    ! -------------------------------------------------------------
    
    ! make sure pointers are not allocated or associated
    
    if ( associated(lon) .or. associated(lat) .or. associated(ts_LaRC) ) then
       err_msg = 'output pointers must not be ' //   &
            'associated or allocated on input.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
       
    !-----------------------------------------
    ! read data into global arrays, perform QC
    !-----------------------------------------

    first = .true.
    
    do j=1,N_files
       
       if (logit) write(logunit,'(400A)') 'reading LaRC Tskin obs from ', trim(fnames(j))

       call Gfio_Open ( fnames(j), 1, fid, rc )
       
       if (rc<0) then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'Error opening gfio file')
       end if
       
       ! temporary catch for files missing HIRESTSKIN (LaRC will fix this) 
       
       call GFIO_DimInquire (fid,g_nlon,g_nlat,km,lm,nvars,ngatts,rc)
       
       dlon=360./real(g_nlon  )
       dlat=180./real(g_nlat-1) 
       
       if (rc<0) then
          err_msg = 'DimInquire error, reading gfio file'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       !if (nvars .LT. 24) then  (for original LaRC files) 
       
       if (nvars .LT. 5) then ! reprocessed files have only 5 variables
          
          write (logunit,*) 'CSD_LaRC file missing a variable ', fnames(j)
          
       else 
          
          if (first) then
             
             ! Attribute cannot be read from LaRC (or reprocessed)  files
             ! no idea why this is not working. ncdump shows attribute is in file
             ! call GFIO_GetRealAtt ( fid, 'missing_value', 1, nodata_LARC, rc )
             !if (rc<0)  call stop_it('Get attribute error, reading nodata from gfio file')
             
             nodata_LARC=-9.99e33 ! For reprocessed files only!!! LaRC use different value
             
             write(logunit,*) &
                  'No-data-value manually set for reprocessed files: ', nodata_LARC
             
             ! get dimensions of grid and allocate temporary arrays
             
             allocate(tmp_data(g_nlon, g_nlat)) 
             allocate(ave_data(g_nlon, g_nlat)) 
             allocate(sum_data(g_nlon, g_nlat)) 
             allocate(cnt_data(g_nlon, g_nlat)) 
             allocate(tmp_vza( g_nlon, g_nlat)) 
             allocate(tmp_sza( g_nlon, g_nlat)) 
             allocate(tmp_fcld(g_nlon, g_nlat)) 
             
             sum_data=0
             cnt_data=0
             
             ! get dimensions of subgrid containing data for this disk

             lat_range=(/-52.0,52.0/)  ! maximum range is +/-52.0 for VZA<60
             
             ! specify boundary for min and max lon (select disk with lowest VZA) 
             
             select case (trim(this_obs_param%descr))
                
             case ('LaRC_tskin-GOESW')
                lon_range=(/-175.   ,-105.   /)   
             case ('LaRC_tskin-GOESE') 
                lon_range=(/-105.   , -36.875/) 
             case ('LaRC_tskin-MET09') 
                lon_range=(/ -36.875,  54.   /) 
             case ('LaRC_tskin-FY2E-')
                ! bounds will not change if not using FY2, as FY2 eastern hemisphere is missing
                lon_range=(/  54.   ,  90.   /)
             case ('LaRC_tskin-MTST2')
                lon_range=(/  90.   , 180.   /)   ! avoid crossing the dateline
             case default     
                err_msg = 'unknown obs_param%descr'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end select

             
             ! find index of smallest lon >  min_lon and of largest lon <= max_lon
             
             x_min=floor( (minval(lon_range)-(-180.))/dlon ) + 2
             x_max=floor( (maxval(lon_range)-(-180.))/dlon ) + 1
             
             ! find index of smallest lat >= min_lat and of largest lat <= max_lat
             
             y_min=floor( (lat_range(1)-(-90.))/dlat ) + 1
             y_max=floor( (lat_range(2)-(-90.))/dlat ) + 1
             
             first = .false. 
             
          endif ! first
          
          ! cannot use actual time, as obs files are rounded to nearest hour
          
          ! ORIGINAL LaRC FILES 
          ! -define time as "minutes since YYYYMMDD, HHMM", and have time=0
          ! read(fnames(j)(len_trim(fnames(j))-8:len_trim(fnames(j))-5),'(i4)') HHMMSS
          ! HHMMSS=HHMMSS*100 ! add seconds
          
          ! lat4d.sh REPROCESSED FILES  
          ! -define time as "hours since YYYYMMDD, HH", and have time=MM/60
          ! The GFIO routine getbegdatetime calculates the time increment in a file by 
          ! reading in first two values
          ! If there is only one time value in the file, the second read (line 240) fails, 
          ! and incSecs=1
          ! This results in TimeIndex in GFIO_GetVar= (MM in seconds)  (rather than 1) 
          ! (in summary, getbegdatetime assumes that if there are not multiple times in 
          ! a file, time must equal 0
          ! which it does not for reprocessed LaRC files)  
          ! Get around this by specifying MM=0 regardless of actual time 
          
          read(fnames(j)(len_trim(fnames(j)) -8:len_trim(fnames(j))-7),'(i2)') HHMMSS
          
          HHMMSS=HHMMSS*10000 ! add seconds, assume minutes are zero
          
          read(fnames(j)(len_trim(fnames(j))-17:len_trim(fnames(j))-9),'(i8)') YYYYMMDD
          

          ! read HIRESTSKIN

          call GFIO_GetVar( fid,'hirestskin',                          &
               YYYYMMDD, HHMMSS, g_nlon, g_nlat,                       &
               0, 1, tmp_data(:,:), rc )
          if (rc<0) then
             err_msg = 'GetVar error, reading hirestskin from gfio file'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          ! Viewing zenith angle (VZA) 

          call GFIO_GetVar( fid,'vza',                                 &
               YYYYMMDD, HHMMSS, g_nlon, g_nlat,                       &
               0, 1, tmp_vza, rc )
          if (rc<0) then
             err_msg = 'GetVar error, reading vza from gfio file'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          ! GOES-East VZA for end Sep, start Oct 2012 are incorrect in the files 
          !  reprocessed by Ben Scarino, replace VZA with that from a hard-coded file
          
          ind = index(fnames(j), trim(tmp_goesEast_vza_dirname))
          
          if (ind/=0) then

             ! extract path to replacement vza file from "fnames(j)"
             ! (=fnames(j) up to and including "tmp_goesEast_vza_dirname")
             
             tmp_fname = fnames(j)(1:ind+len_trim(tmp_goesEast_vza_dirname)-1)
             
             ! append "tmp_goesEast_vza_fname"
             
             tmp_fname = trim(tmp_fname) // '/' // trim(tmp_goesEast_vza_fname)
             
             call Gfio_Open (tmp_fname, 1, tmp_fid, rc )
             
             if (rc<0) then
                err_msg = 'Error opening gfio file'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if
             
             ! read replacement vza

             call GFIO_GetVar( tmp_fid,'vza',                                         &
                  tmp_goesEast_vza_YYYYMMDD, tmp_goesEast_vza_HHMMSS, g_nlon, g_nlat, &
                  0, 1, tmp_vza, rc )
             if (rc<0) then
                err_msg = 'GetVar error, reading vza from gfio file (2)'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if
             
             call GFIO_Close ( tmp_fid, rc )
             
          end if

          
          ! Sun zenith angle (SZA) 

          call GFIO_GetVar( fid,'sza',                                 &
               YYYYMMDD, HHMMSS, g_nlon, g_nlat,                       &
               0, 1, tmp_sza, rc )
          if (rc<0) then
             err_msg = 'GetVar error, reading sza from gfio file'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          

          ! Cloud fraction (FCLD), level 1 is total cloud

          call GFIO_GetVar( fid,'fcld',                                &
               YYYYMMDD, HHMMSS, g_nlon, g_nlat,                       &
               1, 1, tmp_fcld, rc )
          if (rc<0) then
             err_msg = 'GetVar error, reading fcld from gfio file'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if

          ! calculate mean of all data that passes QC
          
          do x=x_min, x_max
             
             do y=y_min, y_max
	    	
                ! Tskin must not be no-data-value
                
                tskin_ok = (abs(tmp_data(x,y) - nodata_LARC) .GE. tol)
                
                ! fcld must not be no-data-value; fcld<=max_cld

                fcld_ok =                                                               &
                     (tmp_fcld(x,y) .LE. max_fcld)              .AND.                   &
                     (abs(tmp_fcld(x,y)-nodata_LARC) .GT. tol)

                ! sza must not be no-data-value; sza<=min_excl_sza; sza>=max_excl_sza
                
                sza_ok =                                                                &
                     (                                                                  &
                     (tmp_sza(x,y) .LE. min_excl_sza)  .OR.                             &
                     (tmp_sza(x,y) .GE. max_excl_sza)                                   &
                     )                                          .AND.                   &
                     (abs(tmp_sza(x,y)-nodata_LARC) .GT. tol)           

                ! vza must not be no-data-value; vza<=vza_max
                
                vza_ok =                                                                &
                     (tmp_vza( x,y) .LE. max_vza)               .AND.                   &
                     (abs(tmp_vza( x,y)-nodata_LARC) .GT. tol)
                
                
                if (tskin_ok .and. fcld_ok .and. sza_ok .and. vza_ok ) then
                   
                   sum_data(x,y) =sum_data(x,y) + tmp_data(x,y)
                   
                   cnt_data(x,y) = cnt_data(x,y)+1
                   
                end if
                
             end do
          end do
          
       end if  ! nvars 
       
       call GFIO_Close ( fid, rc )
       
    end do  ! N_files
    
    
    if ( .not. first) then 
       
       ! only calc averages if found file with HIRESTSKIN - temporary for incomplete files
       
       ! calculate average over appropriate lat/lon range, and count locations with data
       
       N_data=0 
       
       do x=x_min, x_max
          do y=y_min, y_max
             if(cnt_data(x,y)>0) then 
                ave_data(x,y)=sum_data(x,y)/cnt_data(x,y) 	
                N_data=N_data+1
             else
                ave_data(x,y)=nodata_LARC
             endif
          enddo
       enddo
       
       ! allocate pointers for return vectors (must be deallocated outside this subroutine!)
       
       allocate(lon(    N_data))
       allocate(lat(    N_data))
       allocate(ts_LaRC(N_data))
       
       ! pass return data into vectors 
       
       i=1
       do x=x_min, x_max
          do y=y_min, y_max
             if ( abs(ave_data(x,y)-nodata_LARC) .GT. tol ) then 
                ts_LaRC(i)=ave_data(x,y) 
                lon(i)=(x-1)*dlon -180.
                lat(i)=(y-1)*dlat - 90. 
                i=i+1
             endif
          enddo
       enddo
       
       ! clean up 
       
       deallocate(tmp_data) 
       deallocate(sum_data) 
       deallocate(cnt_data) 
       deallocate(ave_data) 
       deallocate(tmp_vza) 
       deallocate(tmp_sza) 
       deallocate(tmp_fcld) 
       
    else
       
       N_data=0
       
       ! OK not to allocate pointers for data?

    endif
    
  end subroutine read_LaRC_Tskin_nc4

  ! *****************************************************************
    
  subroutine read_obs_RedArkOSSE_sm(                             &
       date_time, N_catd, tile_coord, this_obs_param,            &
       found_obs, RedArkOSSE_sm, std_RedArkOSSE_sm )
    
    ! Read observations of surface soil moisture from Wade Crow's
    ! synthetic RedArk OSSE 36km soil moisture files.
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation.
    !
    ! If there are N > 1 observations in a given tile,
    !  a "super-observation" is computed by averaging the N observations.
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! reichle, 25 May 2006
    ! reichle, 26 Sep 2006 - added iostat to open statement
    !
    ! --------------------------------------------------------------------
    
    implicit none
    
    ! inputs:
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: RedArkOSSE_sm
    real,    intent(out), dimension(N_catd) :: std_RedArkOSSE_sm
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals

    ! RedArkOSSE synthetic soil moisture obs at 36km are available once a day.
    !
    ! Mapping of the 36km retrievals to catchment/tile space is precomputed
    ! and stored in p2t_nearest.dat and t2p_nearest.dat
    
    integer, parameter :: N_p2t = 1120  ! # of 36 km pixels
    integer, parameter :: N_t2p =   69  ! # of tiles that need "duplicated" obs
    
    integer, parameter :: HH_obs = 21  ! obs hour of day (UTC) = 3pm CST
    
    ! initial QC parameters

    real,    parameter :: obs_min  = 0.0  ! min allowed obs
    real,    parameter :: obs_max  = 0.45 ! max allowed obs
    real,    parameter :: opac_max = 0.3  ! max allowed vegetation opacity 
    
    integer, parameter :: qc_failed_obs = -888.
    
    integer, dimension(N_p2t)   :: p2t
    integer, dimension(N_t2p,2) :: t2p
    
    character(3)   :: DDD
    character(4)   :: YYYY
    character(300) :: tmpfname
    
    integer :: i, j, ind, N_tmp, tmp_tile_id, istat
    
    real :: tmp_real, tmp_opac
    
    real,    dimension(:), allocatable :: tmp_obs
    integer, dimension(:), allocatable :: tmp_tile_num
    
    integer, dimension(N_catd) :: N_obs_in_tile    

    character(len=*), parameter :: Iam = 'read_obs_RedArkOSSE_sm'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------------
    
    ! initialize
    
    found_obs = .false.
    
    ! obs are available only once per day
    ! (hard-coded b/c time-of-day is not clear from filename)
    
    if (date_time%hour == HH_obs) then
       
       ! read observations:
       !
       ! 1.) read observations, p2t, and t2p from files
       ! 2.) for each observation determine tile_num ("p2t")
       ! 3.) duplicate observations for tiles that are not covered yet ("t2p")
       ! 4.) compute super-obs for each tile from all obs w/in that tile
       !
       ! ----------------------------------------------------------------
       
       ! 1.) a) read obs
       
       write (YYYY,'(i4.4)') date_time%year
       write (DDD, '(i3.3)') date_time%dofyr
       
       tmpfname =  trim(this_obs_param%path) // '/' // YYYY // &
            '/' // trim(this_obs_param%name) // YYYY // '.' // DDD
       
       open(10, file=tmpfname, form='formatted', action='read', iostat=istat)
       
       if (istat==0) then
          
          N_tmp = N_p2t + N_t2p
          
          allocate(tmp_obs(N_tmp))
          allocate(tmp_tile_num(N_tmp))
          
          do i=1,N_p2t
             
             read(10,*) tmp_real, tmp_opac
             
             tmp_obs(i) = tmp_real/100.     ! unit conversion 
             
             ! initial QC
             
             if ( (tmp_opac   > opac_max)     .or.            &
                  (tmp_obs(i) > obs_max)      .or.            &
                  (tmp_obs(i) < obs_min)    )        then
                
                tmp_obs(i) = qc_failed_obs
                
             end if
             
          end do
          
          close(10,status='keep')
          
          ! 1.) b) read p2t
          
          tmpfname =  trim(this_obs_param%path) // '/p2t_nearest.dat'
          
          open(10, file=tmpfname, form='formatted', action='read')
          
          do i=1,N_p2t
             
             read(10,*) p2t(i), tmp_tile_id  ! col 1: tile_num, col 2: tile_id
             
             ! check tile_id for consistency
             
             if (p2t(i)>0) then           ! if statement added 1 May 2007, reichle
                if (tmp_tile_id/=tile_coord(p2t(i))%tile_id) then
                   
                   !write (logunit,*) 'read_obs_RedArkOSSE(): something is wrong (p2t)'
                   !!!write (logunit,*) i, p2t(i), tmp_tile_id, tile_coord(p2t(i))%tile_id
                   !stop

                   err_msg = 'something is wrong (p2t)'
                   call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                   
                end if
             end if
             
          end do
          
          close(10,status='keep')
          
          ! 1.) c) read t2p
          
          tmpfname =  trim(this_obs_param%path) // '/t2p_nearest.dat'
          
          open(10, file=tmpfname, form='formatted', action='read')
          
          do i=1,N_t2p
             
             read(10,*) t2p(i,1:2), tmp_tile_id  ! pixel_num, tile_num, tile_id
             
             if (tmp_tile_id/=tile_coord(t2p(i,2))%tile_id) then
                err_msg = 'something is wrong (t2p)'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if
                
             !write (logunit,*) 'read_obs_RedArkOSSE(): something is wrong (t2p)'
             !stop
             !
             !end if
             
          end do
          
          close(10,status='keep')
          
          ! -----------------------
          !
          ! 2.) for each observation determine tile_num ("p2t")
          
          do i=1,N_p2t
             
             tmp_tile_num(i) = p2t(i)
             
          end do
          
          ! -----------------------
          !
          ! 3.) duplicate observations for tiles not covered yet ("t2p")
          
          do i=1,N_t2p
             
             j = i + N_p2t
             
             tmp_obs(     j) = tmp_obs(t2p(i,1))
             
             tmp_tile_num(j) =         t2p(i,2)
             
          end do
          
          
          ! ------------------------
          !
          ! 4.) compute super-obs for each tile from all obs w/in that tile
          !     (also eliminate observations that are not in domain)
          
          RedArkOSSE_sm = 0.
          N_obs_in_tile = 0
          
          do i=1,N_tmp
             
             ind = tmp_tile_num(i)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
             
             if (ind>0) then         ! this step eliminates obs outside domain
                
                if (tmp_obs(i)>0.) then   ! this step eliminates no-data
                   
                   RedArkOSSE_sm(ind) = RedArkOSSE_sm(ind) + tmp_obs(i)
                   
                   N_obs_in_tile(ind) = N_obs_in_tile(ind) + 1
                   
                end if
                
             end if
             
          end do
          
          ! normalize
          
          do i=1,N_catd
             
             if (N_obs_in_tile(i)>1) then
                
                RedArkOSSE_sm(i) = RedArkOSSE_sm(i)/real(N_obs_in_tile(i))
                
             elseif (N_obs_in_tile(i)==0) then
                
                RedArkOSSE_sm(i) = this_obs_param%nodata
                
             end if
             
          end do
          
          ! --------------------------------
          
          ! set observation error standard deviation
          
          do i=1,N_catd
             std_RedArkOSSE_sm(i) = this_obs_param%errstd
          end do
          
          ! --------------------------------
          
          if (any(N_obs_in_tile>0)) then
             
             found_obs = .true.
             
          else 
             
             found_obs = .false.
             
          end if
          
       end if   ! istat==0
       
       ! clean up
       
       if (allocated(tmp_tile_num)) deallocate(tmp_tile_num)
       if (allocated(tmp_obs))      deallocate(tmp_obs)
       
    end if
    
  end subroutine read_obs_RedArkOSSE_sm
  
  ! *******************************************************************
  
  subroutine read_obs_RedArkOSSE_CLSMsynthSM( date_time, N_catd,       &
       this_obs_param,                                                 &
       found_obs, RedArkOSSE_CLSMsynthSM, std_RedArkOSSE_CLSMsynthSM )
    
    ! Read synthetic observations of surface soil moisture from CLSM
    ! RedArkOSSE integration (generated in matlab from innov output with 
    ! get_RedArk_CLSM_synth_retrievals.m)
    !
    ! synthetic RedArk OSSE 36km soil moisture files.
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation.
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! reichle, 16 Feb 2006
    ! reichle,  5 Apr 2007 - use obs std from default nml input 
    !                        (instead of reading from file)
    !
    ! --------------------------------------------------------------------
    
    implicit none
    
    ! inputs:
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_catd
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: RedArkOSSE_CLSMsynthSM
    real,    intent(out), dimension(N_catd) :: std_RedArkOSSE_CLSMsynthSM
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals

    ! RedArkOSSE CLSM synthetic soil moisture obs are available once a day.
    
    integer, parameter :: HH_obs = 21  ! obs hour of day (UTC) = 3pm CST
    
    real,    parameter :: nodata = -9999.
    
    ! initial QC parameters
    
    real,    parameter :: obs_min  = 0.0  ! min allowed obs
    real,    parameter :: obs_max  = 0.5  ! max allowed obs
    
    character(2)   :: MM, DD
    character(4)   :: YYYY, HHMM
    character(300) :: tmpfname
    
    integer :: i, N_tmp, istat, tmp_tilenum
    
    real :: tmp_obs, tmp_obs_std
    
    ! -------------------------------------------------------------------
    
    ! initialize
    
    found_obs = .false.
    
    ! obs are available only once per day
    ! (hard-coded b/c time-of-day is not very clear in RedArkOSSE
    
    if (date_time%hour == HH_obs) then
       
       write (YYYY,'(i4.4)') date_time%year
       write (MM,  '(i2.2)') date_time%month
       write (DD,  '(i2.2)') date_time%day
       write (HHMM,'(i4.4)') 100*HH_obs
       
       tmpfname =  trim(this_obs_param%path) // '/Y' // YYYY // '/M' // MM // &
            '/' // trim(this_obs_param%name) // YYYY // MM // DD // '_' // HHMM
       
       open(10, file=tmpfname, form='formatted', action='read', iostat=istat)
       
       if (istat==0) then
          
          read(10,*) N_tmp
          read(10,*) ! tmp_obs_std

          ! do NOT use obs std from file 
          ! (s.t. "wrong" obs std can be specified conveniently in nml file)
          ! reichle, 5 Apr 2007
          
          tmp_obs_std = this_obs_param%errstd
          
          RedArkOSSE_CLSMsynthSM(    1:N_catd) = nodata
          std_RedArkOSSE_CLSMsynthSM(1:N_catd) = tmp_obs_std
          
          do i=1,N_tmp
             
             read(10,*) tmp_tilenum, tmp_obs
             
             ! initial QC
             
             tmp_obs = min( obs_max, tmp_obs )
             tmp_obs = max( obs_min, tmp_obs )
             
             RedArkOSSE_CLSMsynthSM( tmp_tilenum ) = tmp_obs
             
          end do
          
          close(10,status='keep')
          
          found_obs = .true.
          
       end if   ! istat==0
       
    end if
    
  end subroutine read_obs_RedArkOSSE_CLSMsynthSM

  ! *****************************************************************
  
  subroutine read_obs_VivianaOK_CLSMsynthSM( date_time, N_catd,      &
       this_obs_param,                                               &
       found_obs, VivianaOK_CLSMsynthSM, std_VivianaOK_CLSMsynthSM )
    
    ! Read synthetic observations of surface soil moisture from CLSM
    ! integration over Viviana's OK domain 
    !
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation.
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! vmaggion + reichle,  4 Aug 2008:
    !    adapted from read_obs_RedArkOSSE_CLSMsynthSM
    !    use obs std from default nml input (instead of reading from file)
    !
    ! --------------------------------------------------------------------
    
    implicit none
    
    ! inputs:
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_catd
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: VivianaOK_CLSMsynthSM
    real,    intent(out), dimension(N_catd) :: std_VivianaOK_CLSMsynthSM
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals

    ! VivianaOK CLSM synthetic soil moisture obs are available once a day.
    
    integer, parameter :: HH_obs = 12  ! obs hour of day (UTC) = 6am CST
    
    real,    parameter :: nodata = -9999.
    
    ! initial QC parameters
    
    real,    parameter :: obs_min  = 0.0  ! min allowed obs
    real,    parameter :: obs_max  = 0.5  ! max allowed obs
    
    character(2)   :: MM, DD, HH
    character(4)   :: YYYY
    character(300) :: tmpfname
    
    integer :: i, istat
    
    real :: tmp_obs, tmp_obs_std
    
    ! -------------------------------------------------------------------
    
    ! initialize
    
    found_obs = .false.
    
    ! obs are available only once per day
    ! (hard-coded b/c time-of-day is not very clear in VivianaOK
    
    if (date_time%hour == HH_obs) then
       
       write (YYYY,'(i4.4)') date_time%year
       write (MM,  '(i2.2)') date_time%month
       write (DD,  '(i2.2)') date_time%day
       write (HH,  '(i2.2)') HH_obs
       
       tmpfname =  trim(this_obs_param%path) // '/' // &
            trim(this_obs_param%name) // YYYY // MM // DD // '_' // HH // 'z.txt'
       
       open(10, file=tmpfname, form='formatted', action='read', iostat=istat)
       
       if (istat==0) then
          
          ! do NOT use obs std from file 
          ! (s.t. "wrong" obs std can be specified conveniently in nml file)
          
          tmp_obs_std = this_obs_param%errstd
          
          VivianaOK_CLSMsynthSM(    1:N_catd) = nodata
          std_VivianaOK_CLSMsynthSM(1:N_catd) = tmp_obs_std
          
          do i=1,N_catd
             
             read(10,*) tmp_obs
             
             ! initial QC
             
             tmp_obs = min( obs_max, tmp_obs )
             tmp_obs = max( obs_min, tmp_obs )
             
             VivianaOK_CLSMsynthSM(i) = tmp_obs
             
             if (abs(tmp_obs-nodata)>1e-4)  found_obs = .true.
             
          end do
          
          close(10,status='keep')
          
       end if   ! istat==0
       
    end if
    
  end subroutine read_obs_VivianaOK_CLSMsynthSM

  ! *****************************************************************
  
  subroutine read_obs_RedArkOSSE_truth(                             &
       date_time, N_catd, this_obs_param,                           &
       found_obs, RedArkOSSE_truth, std_RedArkOSSE_truth )
    
    ! Read "observations" of true surface soil moisture from Wade Crow's
    ! truth soil moisture in CLSM catchment space. 
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation.
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! reichle, 18 Sep 2006
    !
    ! --------------------------------------------------------------------
    
    implicit none
    
    ! inputs:
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_catd
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: RedArkOSSE_truth
    real,    intent(out), dimension(N_catd) :: std_RedArkOSSE_truth
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals

    ! RedArkOSSE soil moisture truth available once a day.
    !
    ! Truth is stored directly in Catchment space
    
    integer, parameter :: HH_obs = 21  ! obs hour of day (UTC) = 3pm CST
    
    character(2), parameter   :: HH = '15'    ! might be CST???
    
    character(3)   :: DDD
    character(4)   :: YYYY
    character(300) :: tmpfname
    
    integer :: i, istat
    
    real :: tmp_real
    
    ! -------------------------------------------------------------------
    
    ! initialize
    
    found_obs = .false.
    
    ! obs are available only once per day
    ! (hard-coded b/c time-of-day is not clear from filename)
    
    if (date_time%hour == HH_obs) then
       
       ! read observations from file
    
       write (YYYY,'(i4.4)') date_time%year
       write (DDD, '(i3.3)') date_time%dofyr
       
       ! file name
       
       tmpfname =  trim(this_obs_param%path) // '/' // YYYY // &
            '/' // trim(this_obs_param%name) // YYYY // '.' // DDD // '.' // HH
       
       ! open file and read obs if available
       
       open(10, file=tmpfname, form='formatted', action='read', iostat=istat)
       
       if (istat==0) then
          
          do i=1,N_catd
             
             read(10,*) tmp_real
             
             RedArkOSSE_truth(i) = tmp_real/100.     ! unit conversion 
             
          end do
          
          close(10,status='keep')
          
          ! set observation error standard deviation
          
          do i=1,N_catd
             std_RedArkOSSE_truth(i) = this_obs_param%errstd
          end do
          
          found_obs = .true.
          
       end if
    end if

  end subroutine read_obs_RedArkOSSE_truth

  ! *****************************************************************
  
  subroutine read_obs_isccp_tskin_gswp2_v1(                      &
       date_time, N_catd, tile_coord,                            &
       tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
       this_obs_param,                                           &
       found_obs, isccp_tskin_gswp2_v1, std_isccp_tskin_gswp2_v1 )
    
    ! Read observations of land skin temperature from ISCCP data
    ! produced by Sarith on GSWP-2 grid
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation.
    !
    ! If there are N > 1 observations in a given tile,
    ! a "super-observation" is computed by averaging the N observations
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! reichle, 26 Sep 2005
    !
    ! --------------------------------------------------------------------
        
    implicit none
    
    ! inputs:
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: isccp_tskin_gswp2_v1
    real,    intent(out), dimension(N_catd) :: std_isccp_tskin_gswp2_v1
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals
    
    integer, parameter :: N_gswp2_compressed = 15238
    
    ! land_i_gswp2 and land_j_gswp2 as stored in 
    ! ISCCP_Tskin_GSWP2_grid_V1 files (by Sarith) follow the GSWP2 convention
    ! for grid orientation, that is counting from north-to-south
    ! and from west-to-east

    real, parameter :: minlon_gswp2 = -180.5
    real, parameter :: maxlat_gswp2 =   90.5
    
    real, parameter :: dx_gswp2 = 1.
    real, parameter :: dy_gswp2 = 1.

    ! parameters for initial quality control and no-data-value treatment
    
    real,    parameter :: tskin_min = 200.
    real,    parameter :: tskin_max = 400.
    
    ! ISCCP_Tskin_GSWP2_grid_V1 files are available every 3h
    
    character(2)   :: HH
    character(4)   :: YYYY, MMDD
    character(300) :: fname
    
    integer :: i, j, ind, istat, N_tmp
    
    real    :: tsclr
    
    integer :: land_i_gswp2, land_j_gswp2
    
    integer, dimension(N_gswp2_compressed) :: tmp_tile_num
    
    real,    dimension(N_gswp2_compressed) :: tmp_obs, tmp_lat, tmp_lon
    
    integer, dimension(N_catd) :: N_obs_in_tile    
    
    ! -------------------------------------------------------------------
    
    ! initialize
    
    found_obs = .false.

    ! assemble file name
    
    write (YYYY,'(i4.4)') date_time%year
    write (MMDD,'(i4.4)') date_time%month*100 + date_time%day
    write (HH,  '(i2.2)') date_time%hour
    
    fname = trim(this_obs_param%path) // '/' // '/Y' // YYYY // &
         '/M' // MMDD(1:2) // '/' // trim(this_obs_param%name) &
         // YYYY // MMDD // '_' // HH // 'z.bin'
    
    open(10, file=fname, form='unformatted', convert='big_endian', &
         action='read', status='old', iostat=istat)
    
    ! read observations:
    !
    ! 1.) read N_tmp observations and their lat/lon info from file
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon
    ! 3.) compute super-obs for each tile from all obs w/in that tile
    !
    ! ----------------------------------------------------------------
    !
    ! 1.) read N_tmp observations and their lat/lon info from file
    
    if (istat==0) then
       
       j = 0
       
       do i=1,N_gswp2_compressed
          
          read (10) land_i_gswp2, land_j_gswp2, tsclr
          
          ! eliminate no-data-values (specify range of acceptable tskin)
          
          if ( (tsclr>tskin_min) .and.         &  
               (tsclr<tskin_max)       ) then     
             
             j=j+1
             
             ! land_i_gswp2 and land_j_gswp2 as stored in 
             ! ISCCP_Tskin_GSWP2_grid_V1 files (by Sarith) follow 
             ! the GSWP2 convention for grid orientation, that is counting 
             ! from north-to-south and from west-to-east
     
             tmp_obs(j) = tsclr
             tmp_lon(j) = minlon_gswp2 + land_i_gswp2*dx_gswp2 
             tmp_lat(j) = maxlat_gswp2 - land_j_gswp2*dy_gswp2 
             
          end if
          
       end do
       
       N_tmp = j
       
       close (10, status='keep')
       
       if (logit) write (logunit,*) 'read_obs_isccp_tskin_gswp2_v1: read ', N_tmp,  &
            ' at date_time = ', date_time, ' from ', trim(fname)
       
    else
       
       N_tmp = 0
       
    end if
        
    ! ----------------------------------------------------------------
    !
    ! 2.) for each observation
    !     a) determine grid cell that contains lat/lon
    !     b) determine tile within grid cell that contains lat/lon
    
    if (N_tmp>0) then
       
       call get_tile_num_for_obs(N_catd, tile_coord,                 &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
            N_tmp, tmp_lat(1:N_tmp), tmp_lon(1:N_tmp),               &
            this_obs_param,                                          &
            tmp_tile_num(1:N_tmp) )
       
       ! ----------------------------------------------------------------
       !
       ! 3.) compute super-obs for each tile from all obs w/in that tile
       !     (also eliminate observations that are not in domain)
       
       isccp_tskin_gswp2_v1 = 0.
       N_obs_in_tile  = 0
       
       do i=1,N_tmp
          
          ind = tmp_tile_num(i)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
          
          if (ind>0) then         ! this step eliminates obs outside domain
             
             isccp_tskin_gswp2_v1(ind) = isccp_tskin_gswp2_v1(ind) + tmp_obs(i)
             
             N_obs_in_tile(ind) = N_obs_in_tile(ind) + 1
             
          end if
          
       end do
       
       ! normalize
       
       do i=1,N_catd
          
          if (N_obs_in_tile(i)>1) then
             
             isccp_tskin_gswp2_v1(i) = &
                  isccp_tskin_gswp2_v1(i)/real(N_obs_in_tile(i))
             
          elseif (N_obs_in_tile(i)==0) then
             
             isccp_tskin_gswp2_v1(i) = this_obs_param%nodata
             
          end if
          
       end do
       
       
       ! --------------------------------
       
       ! set observation error standard deviation
       
       do i=1,N_catd
          std_isccp_tskin_gswp2_v1(i) = this_obs_param%errstd
       end do
       
       ! --------------------------------
       
       if (any(N_obs_in_tile>0)) then
          
          found_obs = .true.
          
       else 
          
          found_obs = .false.
          
       end if
       
    end if
    
  end subroutine read_obs_isccp_tskin_gswp2_v1
  
  ! *****************************************************************
  
  subroutine read_obs_isccp_tskin_ceop3n4(                       &
       date_time, N_catd, tile_coord,                            &
       this_obs_param,                                           &
       found_obs, isccp_tskin_gswp2_v1, std_isccp_tskin_gswp2_v1 )
    
    ! *** ONLY for "CEOP3n4_by_tile_FV_144x91" domain ***
    !
    ! Read observations of land skin temperature from ISCCP data
    ! produced by Sarith on GSWP-2 grid 
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation (even if only "nodata" values in file...).
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! reichle, 27 Jan 2009
    !
    ! --------------------------------------------------------------------
        
    implicit none
    
    ! inputs:
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: isccp_tskin_gswp2_v1
    real,    intent(out), dimension(N_catd) :: std_isccp_tskin_gswp2_v1
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals
    
    integer, parameter :: N_gswp2_compressed = 15238
    
    ! parameters for initial quality control and no-data-value treatment
    
    real,    parameter :: tskin_min = 200.
    real,    parameter :: tskin_max = 400.

    integer, parameter :: N_tiles = 41

    ! the following mapping is from 
    ! /land/l_data/CEOP/EOP3n4/coord/map_CEOP3n4_to_ISCCP_GSWP2.txt
    ! produced with 
    ! land01:/home/reichle/GMAO/station_data/CEOP/EOP3n4/matlab/map_GEOS5_to_ISCCP_GSWP2.m

    ! ISCCP_Tskin_GSWP2_grid_V1 files are available every 3h
    
    character(2)   :: HH
    character(4)   :: YYYY, MMDD
    character(300) :: fname
    
    integer :: i, istat
    
    real, dimension(N_gswp2_compressed)    :: tsclr
    
    real    :: tmp_obs
    
    integer :: land_i_gswp2, land_j_gswp2

    character(len=*), parameter :: Iam = 'read_obs_isccp_tskin_ceop3n4'
    character(len=400) :: err_msg
    
    integer, dimension(2,N_tiles) :: GEOS5_to_ISCCP

    ! ------------------------------------------
         
    
    GEOS5_to_ISCCP = reshape(  &
         (/                    &   
         64402,     96,        & 
         68663,   1687,        &
         68677,   1686,        &
         68771,   1792,        &
         68773,   1630,        &
         68774,   1792,        &
         68775,   1842,        &
         68811,   1791,        &
         68813,   1790,        &
         68814,   1840,        &
         68816,   1686,        &
         68819,   1685,        &
         68836,   1740,        &
         68837,   1629,        &
         68841,   1739,        &
         68842,   1739,        &
         68844,   1684,        &
         68845,   1738,        &
         68849,   1628,        &
         69075,   2146,        &
         69256,   1839,        &
         69301,   1738,        &
         80530,   5920,        &
         80668,   6477,        &
         81200,   7547,        &
         91079,   11144,       &
         91761,   12506,       &
         91762,   12505,       &
         91822,   12569,       &
         92290,   11258,       &
         97082,   13894,       &
         97435,   13833,       &
         99867,   13526,       &
         100925,  11630,       &
         100952,  11522,       &
         101511,  12243,       &
         101689,  12024,       &
         101896,  11578,       &
         101899,  11524,       &
         101901,  11523,       &
         106846,  14836        &
         /),                   &
         shape(GEOS5_to_ISCCP) &
         )
             
    
    ! -------------------------------------------------------------------

    if (N_catd/=N_tiles) then
       err_msg = 'error 1 -- use only for CEOP3n4_by_tile_144x91 domain'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    !write (logunit,*) 'error 1 -- use only for CEOP3n4_by_tile_144x91 domain'
    !stop
    !end if
    
    ! initialize
    
    found_obs = .false.

    isccp_tskin_gswp2_v1     = this_obs_param%nodata
    std_isccp_tskin_gswp2_v1 = this_obs_param%nodata
    
    ! assemble file name
    
    write (YYYY,'(i4.4)') date_time%year
    write (MMDD,'(i4.4)') date_time%month*100 + date_time%day
    write (HH,  '(i2.2)') date_time%hour
    
    fname = trim(this_obs_param%path) // '/' // '/Y' // YYYY // &
         '/M' // MMDD(1:2) // '/' // trim(this_obs_param%name) &
         // YYYY // MMDD // '_' // HH // 'z.bin'
    
    open(10, file=fname, form='unformatted', convert='big_endian', &
         action='read', status='old', iostat=istat)
    
    if (istat==0) then
       
       ! set found_obs=true (even if only nodata values in file)
       
       found_obs = .true.
       
       ! read observations from ISCCP file
       
       do i=1,N_gswp2_compressed
          
          read (10) land_i_gswp2, land_j_gswp2, tsclr(i)
          
       end do
       
       ! extract obs for CEOP3n4_by_tile_FV_144x91 domain
       
       do i=1,N_tiles
          
          ! double-check tile ids
          
          if (GEOS5_to_ISCCP(1,i)/=tile_coord(i)%tile_id) then          
             
             !write (logunit,*) 'error 2 -- only for CEOP3n4_by_tile_144x91 domain'
             !stop

             err_msg = 'error 2 -- only for CEOP3n4_by_tile_144x91 domain'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

          else
             
             tmp_obs = tsclr(GEOS5_to_ISCCP(2,i))
             
             ! basic QC
             
             if ( (tmp_obs>tskin_min) .and. (tmp_obs<tskin_max) ) then     
                
                isccp_tskin_gswp2_v1(i) = tmp_obs
                
             else
                
                isccp_tskin_gswp2_v1(i) = this_obs_param%nodata
                
             end if
             
          end if
          
       end do
       
       close (10, status='keep')
       
       ! --------------------------------
       
       ! set observation error standard deviation
       
       do i=1,N_catd
          std_isccp_tskin_gswp2_v1(i) = this_obs_param%errstd
       end do
       
       
    end if
    
    if (logit) write (logunit,*) 'done with read_obs_isccp_tskin_ceop3n4', &
         ' at date_time = ', date_time, ', read from ', trim(fname)
    
    !write(997,'(41f9.2)') isccp_tskin_gswp2_v1(1:41)
    
  end subroutine read_obs_isccp_tskin_ceop3n4
  
  ! *****************************************************************
  ! *****************************************************************
  
  subroutine read_obs_isccp_ts_ceop3n4_hdASC(          &
       date_time, N_catd, tile_coord,                            &
       this_obs_param,                                           &
       found_obs, isccp_tskin_halfdeg, std_isccp_tskin_halfdeg )
    
    ! *** ONLY for "CEOP3n4_by_tile_FV_144x91" domain ***
    !
    ! read obs from half-deg ASCII ("hdASC") files
    !
    ! Read observations of land skin temperature from ISCCP data
    ! produced by Sarith on half-deg grid and written into ASCII files
    ! for CEOP3n4 locations only
    ! Set flag "found_obs" to true if observations are available 
    !  for assimilation (even if only "nodata" values in file...).
    !
    ! inputs to this subroutine:
    !  date_time = current model date and time 
    !  N_catd    = number of catchments in domain
    !
    ! reichle, 1 Apr 2010
    !
    ! --------------------------------------------------------------------
        
    implicit none
    
    ! inputs:
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: isccp_tskin_halfdeg
    real,    intent(out), dimension(N_catd) :: std_isccp_tskin_halfdeg
    logical, intent(out)                    :: found_obs
    
    ! ---------------
    
    ! locals
    
    integer, parameter :: N_halfdeg = 38
    
    ! parameters for initial quality control and no-data-value treatment
    
    real,    parameter :: tskin_min = 200.
    real,    parameter :: tskin_max = 400.

    integer, parameter :: N_tiles = 41

    ! the following mapping is produced with 
    ! land01:/home/reichle/GMAO/station_data/CEOP/EOP3n4/matlab/map_GEOS5_to_ISCCP_halfdeg.m

    integer, dimension(2,N_tiles) :: GEOS5_to_ISCCP = reshape( (/ &
         64402,        -999,    &
         68663,          12,    &
         68677,          11,    &
         68771,          13,    &
         68773,          18,    &
         68774,          16,    &
         68775,          17,    &
         68811,          21,    &
         68813,          24,    &
         68814,          25,    &
         68816,          15,    &
         68819,          19,    &
         68836,          20,    &
         68837,          22,    &
         68841,          27,    &
         68842,          23,    &
         68844,          26,    &
         68845,          29,    &
         68849,          28,    &
         69075,          32,    &
         69256,          31,    &
         69301,          30,    &
         80530,          37,    &
         80668,          38,    &
         81200,        -999,    &
         91079,           8,    &
         91761,          34,    &
         91762,          35,    &
         91822,          36,    &
         92290,           4,    &
         97082,        -999,    &
         97435,        -999,    &
         99867,          33,    &
         100925,         14,    &
         100952,         10,    &
         101511,          2,    &
         101689,          3,    &
         101896,          6,    &
         101899,          5,    &
         101901,          9,    &
         106846,          1 /), (/2,N_tiles/))

    ! ISCCP files are available every 3h
    
    character(2)   :: HH
    character(4)   :: YYYY, MMDD
    character(300) :: fname
    
    integer :: i, istat
    
    real, dimension(N_halfdeg)    :: tsclr
    
    real    :: tmp_obs
    
    integer :: i_ind, j_ind

    character(len=*), parameter :: Iam = 'read_obs_isccp_ts_ceop3n4_hdASC'
    character(len=400) :: err_msg
    
    ! -------------------------------------------------------------------

    if (N_catd/=N_tiles) then
       err_msg = 'error 1 -- use only for CEOP3n4_by_tile_144x91 domain'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    !write (logunit,*) 'error 1 -- use only for CEOP3n4_by_tile_144x91 domain'
    !stop
    !end if
    
    ! initialize
    
    found_obs = .false.

    isccp_tskin_halfdeg     = this_obs_param%nodata
    std_isccp_tskin_halfdeg = this_obs_param%nodata
    
    ! assemble file name
    
    write (YYYY,'(i4.4)') date_time%year
    write (MMDD,'(i4.4)') date_time%month*100 + date_time%day
    write (HH,  '(i2.2)') date_time%hour
    
    fname = trim(this_obs_param%path) // '/' // '/Y' // YYYY // &
         '/M' // MMDD(1:2) // '/' // trim(this_obs_param%name) &
         // YYYY // MMDD // '_' // HH // 'z.dat'
    
    open(10, file=fname, form='formatted', action='read', status='old', iostat=istat)
    
    if (istat==0) then
       
       ! set found_obs=true (even if only nodata values in file)
       
       found_obs = .true.
       
       ! read observations from ISCCP file
       
       do i=1,N_halfdeg
          
          read (10,*) j_ind, i_ind, tsclr(i)
          
       end do
       
       ! extract obs for CEOP3n4_by_tile_FV_144x91 domain
       
       do i=1,N_tiles
          
          ! double-check tile ids
          
          if (GEOS5_to_ISCCP(1,i)/=tile_coord(i)%tile_id) then          
             
             !write (logunit,*) 'error 2 -- only for CEOP3n4_by_tile_144x91 domain'
             !stop

             err_msg = 'error 2 -- only for CEOP3n4_by_tile_144x91 domain'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             
          else

             if (GEOS5_to_ISCCP(2,i)>0) then
                
                tmp_obs = tsclr(GEOS5_to_ISCCP(2,i))
                
             else
                
                tmp_obs = this_obs_param%nodata
                
             end if
             
             ! basic QC
             
             if ( (tmp_obs>tskin_min) .and. (tmp_obs<tskin_max) ) then     
                
                isccp_tskin_halfdeg(i) = tmp_obs
                
             else
                
                isccp_tskin_halfdeg(i) = this_obs_param%nodata
                
             end if
             
          end if
          
       end do
       
       close (10, status='keep')
       
       ! --------------------------------
       
       ! set observation error standard deviation
       
       do i=1,N_catd
          std_isccp_tskin_halfdeg(i) = this_obs_param%errstd
       end do
       
    end if
    
    if (logit) write (logunit,*) 'done with read_obs_isccp_ts_ceop3n4_hdASC', &
         ' at date_time = ', date_time, ', read from ', trim(fname)
    
    !write(997,'(41f9.2)') isccp_tskin_halfdeg(1:41)
    
  end subroutine read_obs_isccp_ts_ceop3n4_hdASC
  
  ! *****************************************************************

  subroutine read_obs_SMOS( date_time, N_catd, this_obs_param,      &
       dtstep_assim, tile_coord, tile_grid_d,                       &
       N_tile_in_cell_ij, tile_num_in_cell_ij, write_obslog,        &
       found_obs, SMOS_data, std_SMOS_data, SMOS_lon, SMOS_lat   )

    ! reader for preprocessed SMOS Tb and soil moisture (SM) files
    !
    ! - SMOS Tbh/v and SM files contain half-orbit data
    ! - file stamp is UTC to the nearest 30 minutes
    ! - data are already projected onto the SMAP M36 EASE[v2] grid.
    ! - Tb can be "regular" or "fitted" (see LDASsa_DEFAULT_inputs_ensupd.nml)
    !
    ! "this_obs_param%descr" is used to determine which kind of file is read.
    !
    ! "this_obs_param%descr" must contain:
    !
    ! reg_Tb: - 'SMOS_reg_TbX', with X=h or v
    !         - '_A' or '_D', for ascending/descending
    !
    ! fit_Tb: - 'SMOS_fit_TbX', with X=h or v
    !         - '_A' or '_D', for ascending/descending
    !
    ! SM:     - 'SMOS_SM'
    !         - '_A' or '_D', for ascending/descending
    !
    !
    ! GDL, 22Oct10
    ! GDL, 21Mar11 - updated
    ! GDL, 28mar11 - switched fields and angle loops
    !
    ! 26 May 2011, reichle - included into LDASsa
    !  2 Jun 2011, reichle - merged Tb and SM readers
    ! 23 Aug 2013, reichle - added output of SMOS_lon, SMOS_lat
    ! 14 Nov 2013, reichle - added capability to read "fitted" SMOS Tb
    ! 24 Dec 2013, reichle - added output of "obslog" files
    !
    ! --------------------------------------------------------------

    implicit none

    ! inputs:

    type(date_time_type), intent(in) :: date_time
    type(obs_param_type), intent(in) :: this_obs_param
    
    integer,              intent(in) :: dtstep_assim, N_catd

    type(grid_def_type),  intent(in) :: tile_grid_d

    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input

    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input

    logical,              intent(in) :: write_obslog

    ! outputs:

    logical, intent(out)                    :: found_obs

    real,    intent(out), dimension(N_catd) :: SMOS_data, std_SMOS_data
    real,    intent(out), dimension(N_catd) :: SMOS_lon, SMOS_lat

    ! --------------------------------------

    ! local variables
        
    real,      parameter :: SM_min       =   0.   ! min allowed SM
    real,      parameter :: SM_max       =   0.6  ! max allowed SM

    real,      parameter :: SM_DQX_min   =   0.   ! min allowed SM DQX
    real,      parameter :: SM_DQX_max   =   0.2  ! max allowed SM DQX

    real,      parameter :: SM_std_max   =   0.2  ! max SM std-dev  in one M36 grid cell

    integer,   parameter :: SM_cnt_min   =   2    ! min # of SM obs in one M36 grid cell

    real,      parameter :: Tb_min       = 100.0  ! min allowed Tb
    real,      parameter :: Tb_max       = 320.0  ! max allowed Tb
    
    real,      parameter :: reg_Tb_std_max = 7.0  ! max Tb std-dev  in one M36 grid cell
                                                  ! (across contributing 15 km DGG cells)

    integer,   parameter :: reg_Tb_cnt_min = 2    ! min # of Tb obs in one M36 grid cell
                                                  ! (across contributing 15 km DGG cells)

    real,      parameter :: fit_Tb_std_max = 5.0  ! max std-dev between reg. and fitted Tb
                                                  ! (across all angles, a.k.a. fit error)
    
    integer,   parameter :: fit_Tb_cnt_min = 20   ! min # of reg. Tb obs contributing to fit
                                                  ! (across all angles)

    integer,   parameter :: dtstep_file  = 1800   ! time step of pre-proc. SMOS files
    
    integer,   parameter :: unitnum_off  = 10
    
    ! temporarily shift lat/lon of obs for computation of nearest tile to
    ! avoid ambiguous assignment of M09 model tile within M36 obs grid cell
    ! (center of M36 grid cell is equidistant from at least two M09 model 
    !  tiles) -- reichle, 23 Aug 2013
    
    real,      parameter :: tmp_shift_lon = 0.01
    real,      parameter :: tmp_shift_lat = 0.005
    
    ! --------------------------------------------------------

    integer              :: ii, i, j, k, n, istat, secs_in_day
    integer              :: ind_tile, ind_angle, ind_start, ind_end
    integer              :: N_files, N_files_max, N_obs, N_ang
    integer              :: Tb_cnt_min

    logical              :: file_exists, SM_files, reg_Tb_files, hpol, keep_data

    real                 :: M36_col_ind_tile, M36_row_ind_tile  
    real                 :: M36_col_ind_obs,  M36_row_ind_obs
    real                 :: tmpreal, Tb_std_max
    
    type(date_time_type) :: date_time_low, date_time_upp
    
    character(  2)       :: MM, DD, HH, MI, orbit_tag
    character(  4)       :: YYYY
    character( 12)       :: tmpstr12
    character( 16)       :: YYYYMMDD_HHMMSSz
    character( 20)       :: ftag
    character( 80)       :: tmpstr80
    character(300)       :: tmpfname
    
    integer,        dimension(N_catd)           :: N_obs_in_tile

    character(300), dimension(:),   allocatable :: fnames

    integer,        dimension(:),   allocatable :: unitnum    
    integer,        dimension(:,:), allocatable :: start_time, end_time
    integer,        dimension(:),   allocatable :: Asc_flag, N_obs_tmp, N_ang_tmp

    real,           dimension(:),   allocatable :: tmp_lon, tmp_lat, tmp_obs
    real,           dimension(:),   allocatable :: tmp_ang, tmp_std, tmp_DQX
    
    integer,        dimension(:),   allocatable :: tmp_cnt, tmp_tile_num
    integer,        dimension(:),   allocatable :: tmp_file_ind
    
    character(len=*), parameter :: Iam = 'read_obs_SMOS'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------------
    
    ! initialize
    
    found_obs = .false.

    ! read soil moisture or brightness temperature files? 

    if     (index(this_obs_param%descr,'_SM') /= 0) then

       SM_files = .true.                                     ! read SMOS SM files

    elseif (index(this_obs_param%descr,'_Tb') /= 0) then

       SM_files = .false.                                    ! read SMOS Tb files

    else

       err_msg = 'cannot interpret %descr'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if
    
    ! read ascending or descending files? 
    
    if     (index(this_obs_param%descr,'_A') /=0 ) then
       
       orbit_tag = '_A'
       
       if (this_obs_param%orbit/=1) then
          err_msg = 'inconsistent %descr and %orbit'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

    elseif (index(this_obs_param%descr,'_D') /=0 ) then
       
       orbit_tag = '_D'
       
       if (this_obs_param%orbit/=2) then
          err_msg = 'inconsistent %descr and %orbit'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
    else
       
       err_msg = 'unknown %descr or %orbit'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if

    ! read "regular" or "fitted" Tb files?
    
    if (.not. SM_files) then
       
       if     (index(this_obs_param%descr,'_reg') /= 0) then

          reg_Tb_files = .true.                              ! read regular SMOS Tb files
          
          Tb_std_max   = reg_Tb_std_max
          Tb_cnt_min   = reg_Tb_cnt_min
          
       elseif (index(this_obs_param%descr,'_fit') /= 0) then

          reg_Tb_files = .false.                             ! read fitted SMOS Tb files

          Tb_std_max   = fit_Tb_std_max
          Tb_cnt_min   = fit_Tb_cnt_min          

       else
          
          err_msg = 'cannot interpret %descr'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if

    end if
    
    ! ---------------------------
    !
    ! Define search interval for obs
    !
    !   [ date_time - dtstep_assim/2 + one_second, date_time + dtstep_assim/2 ]
    
    ! lower boundary
    
    date_time_low = date_time
    
    call augment_date_time( -dtstep_assim/2+1, date_time_low )

    ! find date_time that is exactly on next half-hour ("dtstep_file")
    
    secs_in_day = date_time_low%hour*3600+date_time_low%min*60+date_time_low%sec
    
    call augment_date_time( -secs_in_day, date_time_low )   ! go to previous 0z
    
    secs_in_day = (secs_in_day/dtstep_file + 1)*dtstep_file 
    
    call augment_date_time(  secs_in_day, date_time_low )   ! go to "next" half-hour
    
    ! upper boundary
    
    date_time_upp = date_time

    call augment_date_time(  dtstep_assim/2,   date_time_upp )

    ! -----------------------------------------------------------------
    !
    ! identify files that exist
    
    N_files_max = (dtstep_assim/dtstep_file) + 1 
    
    allocate(fnames(N_files_max))
    
    k = 0
    
    do while (datetime_le_refdatetime(date_time_low,date_time_upp))
       
       ! assemble file name
       
       if (SM_files) then
          
          ftag = 'SMOS_SM'
          
       else
          
          if (reg_Tb_files) then
             
             ftag = 'SMOS_reg_Tb'
             
          else
             
             ftag = 'SMOS_fit_Tb'
             
          end if
          
       end if
       
       write (YYYY,'(i4.4)') date_time_low%year
       write (MM,  '(i2.2)') date_time_low%month
       write (DD,  '(i2.2)') date_time_low%day 
       write (HH,  '(i2.2)') date_time_low%hour 
       write (MI,  '(i2.2)') date_time_low%min
       
       tmpfname = trim(this_obs_param%path) // '/' // YYYY // MM // '/' // &
            trim(ftag) // '_' // YYYY // MM // DD // '_' // HH // MI // orbit_tag //'.bin' 
       
       inquire(file=tmpfname, exist=file_exists)
       
       if (file_exists) then

          k = k+1
          
          if (k>N_files_max) then
             err_msg = 'too many files found'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          fnames(k) = tmpfname
          
       end if
       
       call augment_date_time( dtstep_file, date_time_low )
       
    end do        ! end do while time step loop
    
    N_files = k

    ! ---------------------------------------------------------------
    !
    ! read data if files were found

    if (N_files>0) then

       allocate(unitnum(    N_files   ))

       allocate(start_time( N_files, 5))
       allocate(end_time(   N_files, 5))
       
       allocate(Asc_flag(   N_files   ))
       allocate(N_obs_tmp(  N_files   ))
       allocate(N_ang_tmp(  N_files   ))
       
       ! open files, read and interpret headers
       
       do k=1,N_files
          
          unitnum(k) = unitnum_off + k
          
          open(unitnum(k), file=trim(fnames(k)),form='unformatted', convert='big_endian',status='old', &
               access='SEQUENTIAL', iostat=istat)
          
          if (istat/=0) then
             
             err_msg = 'could not open file'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             
          else
             
             if (logit) write (logunit,'(400A)') 'reading file ' // trim(fnames(k))
             
          end if
          
          read(unitnum(k)) Asc_flag(  k  )  ! could also read version no.
          read(unitnum(k)) start_time(k,:)
          read(unitnum(k)) end_time(  k,:)
          read(unitnum(k)) N_obs_tmp( k  ),   N_ang_tmp(k)
          
          if (logit) write (logunit,*) '  Asc_flag   = ', Asc_flag(   k  )
          if (logit) write (logunit,*) '  start_time = ', start_time( k,:)
          if (logit) write (logunit,*) '  end_time   = ', end_time(   k,:)
          if (logit) write (logunit,*) '  N_obs_tmp  = ', N_obs_tmp(  k  )
          if (logit) write (logunit,*) '  N_ang_tmp  = ', N_ang_tmp(  k  )
          
       end do
       
       ! make sure N_ang is same in all files
       
       N_ang = N_ang_tmp(1)
       
       do k=2,N_files
          
          if ( N_obs_tmp(k)>0 .and. (N_ang/=N_ang_tmp(k)) ) then
             err_msg = 'angles differ between files'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
       end do
       
       ! allocate data variables
       
       N_obs = sum(N_obs_tmp(1:N_files))
       
       if (N_obs>0) then
          
          allocate(tmp_ang(N_ang))
          
          allocate(tmp_lon(N_obs))
          allocate(tmp_lat(N_obs))
          
          allocate(tmp_obs(N_obs))

          allocate(tmp_std(N_obs))
          allocate(tmp_cnt(N_obs))

          allocate(tmp_file_ind(N_obs))
          
          if (SM_files)  allocate(tmp_DQX(N_obs))
          
          ! loop through files and read data
          
          ind_end   = 0
          
          do k=1,N_files
             
             if (N_obs_tmp(k)>0) then
                
                ind_start = ind_end + 1
                
                ind_end   = ind_start + N_obs_tmp(k) - 1
                
                ! record to which file each obs belongs
                
                tmp_file_ind(ind_start:ind_end) = k

                ! continue with reading file (scalars in header were read above)
                
                read(unitnum(k)) tmp_ang           ! assume same angles in all files
                read(unitnum(k)) tmp_lon(ind_start:ind_end)
                read(unitnum(k)) tmp_lat(ind_start:ind_end)   
                
             
                if (SM_files) then
                   
                   ! read SMOS SM file
                   
                   read(unitnum(k)) tmp_obs(ind_start:ind_end)      !  1 SM
                   read(unitnum(k))                                 !  2 ST
                   read(unitnum(k))                                 !  3 tau
                   read(unitnum(k))                                 !  4 Tbh
                   read(unitnum(k))                                 !  5 Tbv
                   read(unitnum(k)) tmp_DQX(ind_start:ind_end)      !  6 SM RSTD
                   read(unitnum(k))                                 !  7 ST RSTD
                   read(unitnum(k))                                 !  8 tau RSTD
                   read(unitnum(k)) tmp_std(ind_start:ind_end)      !  9 std-dev SM
                   read(unitnum(k)) tmp_cnt(ind_start:ind_end)      ! 10 count SM
                   
                else
                   
                   ! read SMOS Tb file
                   
                   ! first time obs are read: figure out angle and polarization of interest
                   
                   if (ind_start==1) then
                      
                      ! find the index for the angle of interest
                      ! NOTE: after processing of namelist inputs, each species 
                      !       has a unique angle (see subroutine read_ens_upd_inputs())
                      
                      ind_angle = -9999
                      
                      do i=1,N_ang
                         
                         if (abs(tmp_ang(i)-this_obs_param%ang(1))<0.01)  ind_angle = i
                         
                      end do
                      
                      if (ind_angle<0) then
                         err_msg = 'Problem with incidence angle'
                         call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                      end if
                      
                      ! need h-pol or v-pol?
                      
                      if     (this_obs_param%pol==1) then
                         
                         hpol = .true.
                         
                      elseif (this_obs_param%pol==2) then
                         
                         hpol = .false.
                         
                      else
                         
                         err_msg = 'Problem with polarization'
                         call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                         
                      end if
                      
                   end if
                   
                   ! for each field, loop over all angles, keep data at angle of interest
                   
                   do i=1,N_ang           ! (1) Tbh
                      
                      if ( (i==ind_angle) .and. (      hpol) ) then
                         read(unitnum(k)) tmp_obs(ind_start:ind_end)
                      else
                         read(unitnum(k))
                      end if
                      
                   end do
                   
                   do i=1,N_ang           ! (2) Tbv
                      
                      if ( (i==ind_angle) .and. (.not. hpol) ) then
                         read(unitnum(k)) tmp_obs(ind_start:ind_end)
                      else
                         read(unitnum(k))
                      end if
                      
                   end do
                   
                   do i=1,N_ang           ! (3) std-dev Tbh
                      
                      if ( (i==ind_angle) .and. (      hpol) ) then
                         read(unitnum(k)) tmp_std(ind_start:ind_end)
                      else
                         read(unitnum(k))
                      end if
                      
                   end do
                   
                   do i=1,N_ang           ! (4) std-dev Tbv
                      
                      if ( (i==ind_angle) .and. (.not. hpol) ) then
                         read(unitnum(k)) tmp_std(ind_start:ind_end)
                      else
                         read(unitnum(k))
                      end if
                      
                   end do
                   
                   do i=1,N_ang           ! (5) count Tbh
                      
                      if ( (i==ind_angle) .and. (      hpol) ) then
                         read(unitnum(k)) tmp_cnt(ind_start:ind_end)
                      else
                         read(unitnum(k))
                      end if
                      
                   end do
                   
                   do i=1,N_ang           ! (6) count Tbv
                      
                      if ( (i==ind_angle) .and. (.not. hpol) ) then
                         read(unitnum(k)) tmp_cnt(ind_start:ind_end)
                      else
                         read(unitnum(k))
                      end if
                      
                   end do
                   
                   ! additional fields in file that are not currently read:
                   !
                   ! (7-8) RA Tbh-Tbv
                   ! (9-16) repeat the above for T3 and T4
                   
                end if  ! if SM_files

             end if     ! if N_obs_tmp(k)>0
             
          end do        ! loop through files
          
          ! -------------------------------------------------
          !
          ! eliminate no-data-values and data that fail initial QC
          ! and keep track how many obs survived from each file

          N_obs_tmp = 0  ! re-use N_obs_tmp 
          
          j=0
          
          do n=1,N_obs
             
             if (SM_files) then
                
                keep_data =                                 &
                     (tmp_obs( n) >  SM_min)     .and.      &  ! incl: any neg is nodata
                     (tmp_obs( n) <  SM_max)     .and.      &
                     (tmp_DQX( n) >  SM_DQX_min) .and.      &
                     (tmp_DQX( n) <  SM_DQX_max) .and.      &
                     (tmp_std( n) <  SM_std_max) .and.      &
                     (tmp_cnt( n) >= SM_cnt_min)         
                
             else
                
                keep_data =                                 &
                     (tmp_obs( n) >  Tb_min)     .and.      &  ! incl: any neg is nodata
                     (tmp_obs( n) <  Tb_max)     .and.      &
                     (tmp_std( n) <  Tb_std_max) .and.      &
                     (tmp_cnt( n) >= Tb_cnt_min)         
                
             end if
             
             if (keep_data) then
                
                j=j+1
                
                tmp_obs(j) = tmp_obs(n)
                tmp_lon(j) = tmp_lon(n)
                tmp_lat(j) = tmp_lat(n)

                N_obs_tmp(tmp_file_ind(n)) = N_obs_tmp(tmp_file_ind(n)) + 1
                
             end if
             
          end do
          
          N_obs = j  ! Note: This is NOT the final number of valid obs in "SMOS_data"!
          
          if (SM_files)  deallocate(tmp_DQX)
          
          deallocate(tmp_std)
          deallocate(tmp_cnt)
          
          deallocate(tmp_file_ind)          
          
          ! -------------------------------------------------
          !
          ! map obs to tiles
          !
          ! for each observation
          !     a) determine grid cell that contains lat/lon
          !     b) determine tile within grid cell that contains lat/lon
          
          if (N_obs>0) then
             
             allocate(tmp_tile_num(N_obs))
             
             ! temporarily shift lat/lon of obs for computation of nearest tile to
             ! avoid ambiguous assignment of M09 model tile within M36 obs grid cell
             ! (center of M36 grid cell is equidistant from at least two M09 model 
             !  tiles) -- reichle, 23 Aug 2013
             
             call get_tile_num_for_obs(N_catd, tile_coord,                 &
                  tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
                  N_obs, tmp_lat, tmp_lon,                                 &
                  this_obs_param,                                          &
                  tmp_tile_num, &
                  tmp_shift_lat, tmp_shift_lon )
             
             ! make sure center-of-mass of tile that administers obs 
             !  is within EASEv2 M36 obs grid cell, discard obs otherwise
             !  (by setting tmp_tile_num to negative value)
             ! - reichle, 31 Jan 2014
             !
             ! It is not 100 percent clear why this piece code had been added
             !  at the time.
             ! Chances are that it had to do with land-water mask issues and the
             !  distortion of EASE grid cells at high latitudes, and/or Gabrielle's
             !  use of an M36 innovations integration to derive Tb scaling files
             !  for the M09 (SMAP) system.
             ! The problem with the piece of code is that it throws out far too many 
             !  obs (at all latitudes) if the tile space is coarse, e.g., that of the 
             !  1/2 deg Lat/Lon grids of MERRA or MERRA-2.
             ! The piece of code may no longer be needed because of improvements in
             !  the obs readers and in get_obs_pred(), but the impact of removing the
             !  code on the SMAP L4_SM system is not clear.  At this time, just prior
             !  to finalizing the L4_SM "validated release", keep the code for the 
             !  EASEv2 M09 and M36 tile spaces that are relevant for the L4_SM 
             !  system, but drop it for all other tile spaces.
             ! - reichle, 3 Feb 2016
             
             if   (                                                               &
                  (index(tile_grid_d%gridtype, 'EASEv2_M09')  /=0) .or.           &
                  (index(tile_grid_d%gridtype, 'EASEv2_M36')  /=0)       )  then
                
                do ii=1,N_obs
                
                   if (tmp_tile_num(ii)>0) then
                      
                      call ease_convert('EASEv2_M36',             &
                           tile_coord(tmp_tile_num(ii))%com_lat,  &
                           tile_coord(tmp_tile_num(ii))%com_lon,  &
                           M36_col_ind_tile, M36_row_ind_tile  )
                      
                      call ease_convert('EASEv2_M36',             &
                           tmp_lat(ii),                           &
                           tmp_lon(ii),                           &
                           M36_col_ind_obs,  M36_row_ind_obs  )
                      
                      if ( (nint(M36_col_ind_tile)/=nint(M36_col_ind_obs))  .or.     &
                           (nint(M36_row_ind_tile)/=nint(M36_row_ind_obs))        )  &
                           tmp_tile_num(ii) = -9999
                      
                   end if
                   
                end do
             
             end if

             ! compute super-obs for each tile from all obs w/in that tile
             !     (also eliminate observations that are not in domain)
             
             SMOS_data = 0.
             SMOS_lon  = 0.
             SMOS_lat  = 0.

             N_obs_in_tile  = 0
             
             do i=1,N_obs
                
                ind_tile = tmp_tile_num(i)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
                
                if (ind_tile>0) then         ! this step eliminates obs outside domain
                   
                   SMOS_data(ind_tile) = SMOS_data(ind_tile) + tmp_obs(i)
                   SMOS_lon( ind_tile) = SMOS_lon( ind_tile) + tmp_lon(i)
                   SMOS_lat( ind_tile) = SMOS_lat( ind_tile) + tmp_lat(i)
                   
                   N_obs_in_tile(ind_tile) = N_obs_in_tile(ind_tile) + 1
                   
                end if
                
             end do
             
             ! normalize
             
             do i=1,N_catd
                
                if (N_obs_in_tile(i)>1) then
                   
                   tmpreal = real(N_obs_in_tile(i))

                   SMOS_data(i) = SMOS_data(i)/tmpreal
                   SMOS_lon( i) = SMOS_lon( i)/tmpreal
                   SMOS_lat( i) = SMOS_lat( i)/tmpreal
                   
                elseif (N_obs_in_tile(i)==0) then
                   
                   SMOS_data(i) = this_obs_param%nodata
                   SMOS_lon( i) = this_obs_param%nodata
                   SMOS_lat( i) = this_obs_param%nodata
                   
                end if
                
             end do
             
             ! clean up
             
             deallocate(tmp_tile_num)
             
             ! --------------------------------
             
             ! set observation error standard deviation
             
             do i=1,N_catd
                std_SMOS_data(i) = this_obs_param%errstd
             end do
             
             ! --------------------------------
             
             if (any(N_obs_in_tile>0)) then
             
                found_obs = .true.
                
             else 
                
                found_obs = .false.
                
             end if
             
          end if
          
          deallocate(tmp_ang)
       
          deallocate(tmp_lon)
          deallocate(tmp_lat)
          
          deallocate(tmp_obs)
          
       end if
       
       do k=1,N_files
          
          close(unitnum(k),status='keep')
          
       end do
       
       deallocate(unitnum)

       deallocate(start_time)
       deallocate(end_time)         
       deallocate(Asc_flag)  
       deallocate(N_ang_tmp) 
       
    end if    ! if N_files>0

    ! -------------------------------------------------
    !
    ! write "obslog" file
    
    if (write_obslog) then
       
       YYYYMMDD_HHMMSSz = date_time2string(date_time)
       
       tmpstr80 = 'read_obs_SMOS()'  ! name of this subroutine
       
       do k=1,N_files
          
          write (tmpstr12,'(i12)') N_obs_tmp(k)   ! convert integer to string
          
          call add_to_obslog( YYYYMMDD_HHMMSSz, this_obs_param%descr, tmpstr80, &
               tmpstr12, fnames(k) )
          
       end do
       
    end if

    ! clean up

    if (N_files>0) deallocate(N_obs_tmp) 
    
    deallocate(fnames)
    
    if (logit) write (logunit,*) 'read_obs_SMOS(): done.'
    
  end subroutine read_obs_SMOS

  ! *****************************************************************
  
  subroutine read_obs_SMAP_FT( date_time, N_catd, this_obs_param,            &
       dtstep_assim, tile_coord, tile_grid_d,                                &
       N_tile_in_cell_ij, tile_num_in_cell_ij, write_obslog,                 &
       found_obs, SMAP_data, std_SMAP_data, SMAP_lon, SMAP_lat, SMAP_time   )
    
    ! read freeze/thaw (FT) data within the assimilation window from one or more 
    !  SMAP L2_SM_AP half-orbit h5 files  (or L3_FT_A files -- TO BE IMPLEMENTED)
    !
    ! this subroutine reads each species independently of the others; 
    ! TO DO: what if more than one flavor of FT data is read from SMAP?
    !
    ! reichle, 14 Nov 2014
    
    implicit none
    
    ! inputs:
    
    type(date_time_type), intent(in) :: date_time
    
    integer,              intent(in) :: N_catd, dtstep_assim

    type(obs_param_type), intent(in) :: this_obs_param

    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type),  intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input

    logical,              intent(in) :: write_obslog

    ! outputs:
    
    logical, intent(out)                    :: found_obs

    real,    intent(out), dimension(N_catd) :: SMAP_data, std_SMAP_data
    real,    intent(out), dimension(N_catd) :: SMAP_lon,  SMAP_lat
    real*8,  intent(out), dimension(N_catd) :: SMAP_time                 ! J2000 seconds

    ! --------------------------------------

    ! local variables
    
    character(len=*),  parameter :: Iam = 'read_obs_SMAP_FT'

    character(4),      parameter :: J2000_epoch_id = 'TT12'    ! see date_time_util.F90

    ! SMAP data are stored in Yyyyy/Mmm/Ddd directories
    !
    ! Within each directory, there must be an ASCII text file that lists the h5 file 
    !  names of all files in the directory that are suitable for assimilation.
    ! These ASCII files must be named as follows:
    
    character(80), parameter :: fname_of_fname_list_L2_SM_AP_A = 'SMAP_L2_SM_AP_A_list.txt'
    character(80), parameter :: fname_of_fname_list_L2_SM_AP_D = 'SMAP_L2_SM_AP_D_list.txt'
    character(80), parameter :: fname_of_fname_list_L3_FT_A    = 'SMAP_L3_FT_A_list.txt'
    
    logical,   parameter :: tmp_debug    = .false.
    
    real,      parameter :: FT_min       = 0.  ! min allowed FT
    real,      parameter :: FT_max       = 1.  ! max allowed FT

    integer,   parameter :: dt_halforbit     = 50*60 ! seconds
    
    integer,   parameter :: N_halforbits_max = 15    ! max number of half-orbits per day
    
    integer,   parameter :: dtstep_assim_max = 10800 ! max allowed dtstep_assim [seconds]
    
    ! get max number of files to be read  (use 2*dt_halforbit just to be safe)
    
    integer,   parameter :: N_fnames_max =                                             &
         ceiling(real(N_halforbits_max)/86400.*real(dtstep_assim_max+2*dt_halforbit))
    
    ! --------------------------------------------

    type(hdf5read)       :: h5r
    
    logical              :: L2AP_files, keep_data, file_exists
    
    type(date_time_type) :: date_time_low,       date_time_upp
    type(date_time_type) :: date_time_low_fname, date_time_tmp

    integer              :: ii, jj, kk, nn
    integer              :: N_fnames, N_fnames_tmp, N_obs_tmp
    integer              :: dset_rank
    integer              :: ind_tile, ind_start, ind_end

    real                 :: M09_col_ind_tile, M09_row_ind_tile  
    real                 :: M09_col_ind_obs,  M09_row_ind_obs
    real                 :: tmpreal
    
    real*8               :: J2000_seconds_low, J2000_seconds_upp

    character( 12)       :: tmpstr12
    character( 15)       :: SMAP_date_time
    character( 16)       :: YYYYMMDD_HHMMSSz 
    character( 80)       :: fname_of_fname_list, tmpstr80
    character(300)       :: fname_tmp, tmp_err_msg

    character(100)       :: dset_name_lon,  dset_name_lat
    character(100)       :: dset_name_time, dset_name_ft, dset_name_ft_qual_flag

    character(100), dimension(2*N_halforbits_max)  :: fname_list  ! max 2 days of files

    integer,        dimension(7)                   :: dset_size
    integer,        dimension(N_fnames_max)        :: N_obs_kept
    integer,        dimension(N_catd)              :: N_obs_in_tile

    real,           dimension(:),     allocatable  :: tmp_lon, tmp_lat
    
    real,           dimension(:),     allocatable  :: tmp_ft

    real*8,         dimension(:),     allocatable  :: tmp_time

    integer,        dimension(:),     allocatable  :: tmp_ft_qual_flag

    integer,        dimension(:),     allocatable  :: tmp_tile_num

    character(len=400)           :: err_msg

    ! -------------------------------------------------------------------
    
    ! check inputs
    
    ! the subroutine makes sense only if dtstep_assim <= 3 hours
    !
    ! (this avoids that more than 2 different Yyyyy/Mmm/Ddd directories are needed
    !  and that the time mismatch between the observed Tb and the model forecast Tb 
    !  becomes excessive;  in future, add time interpolation of forecast Tb)
    
    if (dtstep_assim > dtstep_assim_max) then
       err_msg = 'dtstep_assim must not exceed 3 hours'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end  if
    
    ! initialize
    
    found_obs = .false.
    
    ! read L2_SM_AP or L3_FT_A files? 
    
    if     (index(this_obs_param%descr,'L2AP') /= 0) then
       
       L2AP_files = .true.                                  ! read SMAP L2_SM_AP files
       
    elseif (index(this_obs_param%descr,'L3FT') /= 0) then
       
       L2AP_files = .false.                                 ! read SMAP L3_FT_A files

    else
       
       err_msg = 'cannot interpret %descr'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if
 
        
    ! determine the name of the file that contain the relevant list of file names
    
    if     (L2AP_files .and. this_obs_param%orbit==1) then
       
       fname_of_fname_list = trim(fname_of_fname_list_L2_SM_AP_A)
       
    elseif (L2AP_files .and. this_obs_param%orbit==2) then
       
       fname_of_fname_list = trim(fname_of_fname_list_L2_SM_AP_D)
       
    else
       
       fname_of_fname_list = trim(fname_of_fname_list_L3_FT_A)
       
    end if
    
    ! ---------------------------
    !
    ! define h5 data set names
    
    if (L2AP_files) then
       
       ! L2_SM_AP
       
       dset_name_lon           = '/Soil_Moisture_Retrieval_Data/longitude'
       dset_name_lat           = '/Soil_Moisture_Retrieval_Data/latitude'
       
       dset_name_time          = '/Soil_Moisture_Retrieval_Data/spacecraft_overpass_time_seconds'

       dset_name_ft            = '/Soil_Moisture_Retrieval_Data/freeze_thaw_fraction'
       
       dset_name_ft_qual_flag  = '/Soil_Moisture_Retrieval_Data/NOT_YET_IMPLEMENTED'

    else
       
       ! L3FT

       err_msg = 'L3FT NOT YET IMPLEMENTED'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    ! ---------------------------
    !
    ! Define search interval for obs
    !
    !   [ date_time - dtstep_assim/2 + one_second, date_time + dtstep_assim/2 ]
    
    ! lower boundary
    
    date_time_low = date_time
    
    call augment_date_time( -dtstep_assim/2+1, date_time_low )
    
    J2000_seconds_low = datetime_to_J2000seconds( date_time_low, J2000_epoch_id )

    ! upper boundary
    
    date_time_upp = date_time

    call augment_date_time(  dtstep_assim/2,   date_time_upp )

    J2000_seconds_upp = datetime_to_J2000seconds( date_time_upp, J2000_epoch_id )

    ! -----------------------------------------------------------------
    !
    ! identify names of files that could contain data within the
    !  assimilation window [date_time_low,date_time_upp]

    date_time_low_fname = date_time_low
    
    call augment_date_time( -dt_halforbit, date_time_low_fname )
        
    ! read file with list of SMAP file names for first day
    
    call read_obs_SMAP_fnames( date_time_low_fname, this_obs_param,                &
         fname_of_fname_list, N_halforbits_max,                                    &
         N_fnames, fname_list(1:N_halforbits_max) )
    
    ! if needed, read file with list of SMAP file names for second day and add
    !  file names into "fname_list"
    
    if (date_time_low_fname%day /= date_time_upp%day) then
       
       call read_obs_SMAP_fnames( date_time_upp, this_obs_param,                   &
            fname_of_fname_list, N_halforbits_max,                                 &
            N_fnames_tmp, fname_list((N_fnames+1):(N_fnames+N_halforbits_max)) )
       
       N_fnames = N_fnames + N_fnames_tmp
       
    end if

    ! ------------------------------------------------------------------
    
    ! extract names of files that could contain data within the assimilation
    !  window
    !
    ! sample file names:
    !
    !   Yyyyy/Mmm/Ddd/SMAP_L2_SM_AP_03073_D_20010730T193828_D04003_000.h5
    !                                       ||||||||||||||| 
    ! counter:   1         2         3         4         5         6
    !   1234567890123456789012345678901234567890123456789012345678901234567890

    if (L2AP_files) then
       
       ind_start = 37
       ind_end   = 51
       
    else
       
       ind_start = -9999999    ! TO DO: IMPLEMENT F3FT
       ind_end   = -9999999
       
    end if

    kk = 0
    
    do ii=1,N_fnames

       SMAP_date_time = fname_list(ii)(ind_start:ind_end)
       
       date_time_tmp = SMAPdatetime_to_DateTimeType( SMAP_date_time )

       ! check whether:  date_time_low_fname < date_time_tmp <= date_time_upp

       if ( datetime_lt_refdatetime( date_time_low_fname, date_time_tmp ) .and.        &
            datetime_le_refdatetime( date_time_tmp,       date_time_upp )       ) then 
          
          kk = kk+1

          ! there can be no more than N_fnames_max files that have data falling into the 
          ! assimilation window (see also dtstep_assim_max)
          
          if (kk>N_fnames_max) then
             err_msg = 'too many files  match assimilation window'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          fname_list(kk) = fname_list(ii)
          
       end if
       
    end do  ! ii=1,N_fnames
    
    N_fnames = kk
    
    ! the "N_fnames" files of interest are now the first "N_fnames" entries
    !  in "fname_list"
    
    ! ---------------------------------------------------------
    !
    ! read and process data if files were found
    
    if (N_fnames==0) then
       
       ! no data files found
       
       SMAP_data     =      this_obs_param%nodata
       SMAP_lon      =      this_obs_param%nodata
       SMAP_lat      =      this_obs_param%nodata
       SMAP_time     = real(this_obs_param%nodata,kind(0.0D0))
       std_SMAP_data =      this_obs_param%nodata
       
    else
       
       ! initialize outputs
       
       SMAP_data = 0.
       SMAP_lon  = 0.
       SMAP_lat  = 0.
       SMAP_time = 0.0D0

       N_obs_in_tile  = 0  ! for normalization after mapping to tile and super-obs
       
       ! loop through files
              
       do kk=1,N_fnames
          
          ! open file
          
          fname_tmp = trim(this_obs_param%path) // '/' // fname_list(kk)
          
          if (logit) write(logunit,'(400A)') 'reading file: ' // trim(fname_tmp)
          
          inquire(file=fname_tmp, exist=file_exists)
          
          if (.not. file_exists) then
             
             err_msg = 'file does NOT exist'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)             

          end if

          call h5r%openFile(fname_tmp)
          
          ! ------------------------------

          ! read h5 datasets

          tmp_err_msg = trim(Iam) // ': inconsistent dataset lengths'

          ! LONGITUDE: query dataset, record size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_lon)

          call h5r%queryDataset(dset_name_lon, dset_rank, dset_size)
          
          N_obs_tmp = dset_size(1)
          
          allocate(tmp_lon(N_obs_tmp))
          
          call h5r%readDataset(tmp_lon)
          
          ! LATITUDE: query dataset, check size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_lat)

          call h5r%queryDataset(dset_name_lat, dset_rank, dset_size)
          
          if (N_obs_tmp/=dset_size(1)) then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
          end if
          
          allocate(tmp_lat(N_obs_tmp))

          call h5r%readDataset(tmp_lat)
                     
          ! TIME: query dataset, check size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_time)
          
          call h5r%queryDataset(dset_name_time, dset_rank, dset_size)
          
          if (N_obs_tmp/=dset_size(1)) then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
          end if
          
          allocate(tmp_time(N_obs_tmp))

          call h5r%readDataset(tmp_time)

          ! FT: query dataset, check size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_ft)
          
          call h5r%queryDataset(dset_name_ft, dset_rank, dset_size)
          
          if (N_obs_tmp/=dset_size(1)) then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
          end if
          
          allocate(tmp_ft(N_obs_tmp))
          
          call h5r%readDataset(tmp_ft)
          
          !! FT_QUAL_FLAG: query dataset, check size, allocate space, read data
          ! 
          !if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_ft_qual_flag)
          !
          !call h5r%queryDataset(dset_name_ft_qual_flag, dset_rank, dset_size)
          ! 
          !if (N_obs_tmp/=dset_size(1)) then
          !   call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
          !end if
          !
          !allocate(tmp_ft_qual_flag(N_obs_tmp))
          !
          !call h5r%readDataset(tmp_ft_qual_flag)

          
          ! close file
          
          call h5r%closeFile

          if (logit) write(logunit,*) 'done reading file'

          ! --------------------------------------------------------
          !
          ! eliminate obs outside desired time window, no-data-values and data 
          !  that fail initial QC
          ! keep track of how many obs survived from current file

          ! #####################################
          ! TO DO: 
          ! - use quality flag once available
          ! #####################################
          
          jj = 0
          
          if (L2AP_files) then
             
             do nn=1,N_obs_tmp
                
                !(mod(tmp_tb_qual_flag_1(nn),2)==0)    .and. & ! lowest bit must be 0
                
                keep_data =                                    &
                     (tmp_time(nn) >  J2000_seconds_low) .and. &
                     (tmp_time(nn) <= J2000_seconds_upp) .and. &
                     (tmp_ft(nn)   >  FT_min)            .and. & ! elim neg nodata
                     (tmp_ft(nn)   <  FT_max)                    ! elim unphysically large value
                
                if (keep_data) then
                   
                   jj=jj+1
                   
                   tmp_lon( jj) = tmp_lon( nn)
                   tmp_lat( jj) = tmp_lat( nn)
                   tmp_ft(  jj) = tmp_ft(  nn)
                   tmp_time(jj) = tmp_time(nn)
                   
                end if
                
             end do

          else
             
             ! L3FT
             
             err_msg = 'L3FT NOT YET IMPLEMENTED'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)             
             
          end if  ! (L2AP_files)
          
          N_obs_kept(kk) = jj
          
          ! -------------------------------------------------
          !
          ! map obs to tiles
          !
          ! for each observation
          !     a) determine grid cell that contains lat/lon
          !     b) determine tile within grid cell that contains lat/lon
          
          if (N_obs_kept(kk)>0) then
              
             allocate(tmp_tile_num(N_obs_kept(kk)))
             
             call get_tile_num_for_obs(N_catd, tile_coord,                 &
                  tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
                  N_obs_kept(kk),                                          &
                  tmp_lat(1:N_obs_kept(kk)),                               &
                  tmp_lon(1:N_obs_kept(kk)),                               &
                  this_obs_param,                                          &
                  tmp_tile_num )
             
             ! make sure center-of-mass of tile that administers obs 
             !  is within EASEv2 M09 obs grid cell, discard obs otherwise
             !  (by setting tmp_tile_num to negative value)
             !
             ! this step eliminates SMAP obs that fall outside of the GEOS-5 
             ! land mask at M09 resolution
             !
             ! not sure what this does if model is run on tile space other than 
             ! EASEv2_M09 - reichle, 11 Nov 2014
             !
             ! It is not 100 percent clear why this piece code had been added
             !  at the time.
             ! Chances are that it had to do with land-water mask issues and the
             !  distortion of EASE grid cells at high latitudes, and/or Gabrielle's
             !  use of an M36 innovations integration to derive Tb scaling files
             !  for the M09 (SMAP) system.
             ! The problem with the piece of code is that it throws out far too many 
             !  obs (at all latitudes) if the tile space is coarse, e.g., that of the 
             !  1/2 deg Lat/Lon grids of MERRA or MERRA-2.
             ! The piece of code may no longer be needed because of improvements in
             !  the obs readers and in get_obs_pred(), but the impact of removing the
             !  code on the SMAP L4_SM system is not clear.  At this time, just prior
             !  to finalizing the L4_SM "validated release", keep the code for the 
             !  EASEv2 M09 and M36 tile spaces that are relevant for the L4_SM 
             !  system, but drop it for all other tile spaces.
             ! - reichle, 3 Feb 2016
             
             if   (                                                               &
                  (index(tile_grid_d%gridtype, 'EASEv2_M09')  /=0) .or.           &
                  (index(tile_grid_d%gridtype, 'EASEv2_M36')  /=0)       )  then
                
                do ii=1,N_obs_kept(kk)
                   
                   if (tmp_tile_num(ii)>0) then
                      
                      call ease_convert('EASEv2_M09',                  &
                           tile_coord(tmp_tile_num(ii))%com_lat,  &
                           tile_coord(tmp_tile_num(ii))%com_lon,  &
                           M09_col_ind_tile, M09_row_ind_tile  )
                      
                      call ease_convert('EASEv2_M09',                  &
                           tmp_lat(ii),                           &
                           tmp_lon(ii),                           &
                           M09_col_ind_obs,  M09_row_ind_obs  )
                      
                      if ( (nint(M09_col_ind_tile)/=nint(M09_col_ind_obs))  .or.     &
                           (nint(M09_row_ind_tile)/=nint(M09_row_ind_obs))        )  &
                           tmp_tile_num(ii) = -9999
                      
                   end if
                   
                end do

             end if

             ! compute super-obs for each tile from all obs w/in that tile
             !     (also eliminate observations that are not in domain)
             
             do ii=1,N_obs_kept(kk)
                
                ind_tile = tmp_tile_num(ii)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
                
                if (ind_tile>0) then          ! this step eliminates obs outside domain
                
                   SMAP_data(ind_tile) = SMAP_data(ind_tile) + tmp_ft(  ii)
                   SMAP_lon( ind_tile) = SMAP_lon( ind_tile) + tmp_lon( ii)
                   SMAP_lat( ind_tile) = SMAP_lat( ind_tile) + tmp_lat( ii)
                   SMAP_time(ind_tile) = SMAP_time(ind_tile) + tmp_time(ii)
                   
                   N_obs_in_tile(ind_tile) = N_obs_in_tile(ind_tile) + 1
                   
                end if
                
             end do
             
             deallocate(tmp_tile_num)
             
          end if
          
          ! clean up
          
          if (allocated(tmp_lon         )) deallocate(tmp_lon         )
          if (allocated(tmp_lat         )) deallocate(tmp_lat         )
          if (allocated(tmp_time        )) deallocate(tmp_time        )
          if (allocated(tmp_ft          )) deallocate(tmp_ft          )
          if (allocated(tmp_ft_qual_flag)) deallocate(tmp_ft_qual_flag)
          
       end do  ! kk=1,N_fnames
          
       ! normalize
       
       do  ii=1,N_catd
          
          if (N_obs_in_tile(ii)>1) then
             
             tmpreal = real(N_obs_in_tile(ii))
             
             SMAP_data(ii) = SMAP_data(ii)/     tmpreal
             SMAP_lon( ii) = SMAP_lon( ii)/     tmpreal
             SMAP_lat( ii) = SMAP_lat( ii)/     tmpreal
             SMAP_time(ii) = SMAP_time(ii)/real(tmpreal,kind(0.0D0))
             
          elseif (N_obs_in_tile(ii)==0) then
             
             SMAP_data(ii) =      this_obs_param%nodata
             SMAP_lon( ii) =      this_obs_param%nodata
             SMAP_lat( ii) =      this_obs_param%nodata
             SMAP_time(ii) = real(this_obs_param%nodata,kind(0.0D0))
             
          end if
          
       end do
       
       ! --------------------------------
       
       ! set observation error standard deviation
       
       do ii=1,N_catd
          std_SMAP_data(ii) = this_obs_param%errstd
       end do
       
       ! --------------------------------
       
       if (any(N_obs_in_tile>0)) then
          
          found_obs = .true.
          
       else 
          
          found_obs = .false.
          
       end if
       
    end if  ! N_fnames==0
     
    ! -------------------------------------------------
    !
    ! write "obslog" file

    if (write_obslog) then

       YYYYMMDD_HHMMSSz = date_time2string(date_time)

       tmpstr80 = 'read_obs_SMAP_FT()'  ! name of this subroutine
              
       do kk=1,N_fnames
          
          fname_tmp = trim(this_obs_param%path) // fname_list(kk)
          
          write (tmpstr12,'(i12)') N_obs_kept(kk)    ! convert integer to string
          
          call add_to_obslog( YYYYMMDD_HHMMSSz, this_obs_param%descr, tmpstr80, &
               tmpstr12, fname_tmp )
          
       end do
       
    end if
    
    ! clean up
    
    if (logit) write (logunit,*) 'read_obs_SMAP_FT(): done.'
    
  end subroutine read_obs_SMAP_FT
  
  ! *****************************************************************

  ! *****************************************************************
  
  subroutine read_obs_SMAP_halforbit_Tb( date_time, N_catd, this_obs_param,  &
       dtstep_assim, tile_coord, tile_grid_d,                                &
       N_tile_in_cell_ij, tile_num_in_cell_ij, write_obslog,                 &
       found_obs, SMAP_data, std_SMAP_data, SMAP_lon, SMAP_lat, SMAP_time   )
    
    ! read brightness temperature data within the assimilation window from one or more 
    !  SMAP half-orbit h5 files (L1C, L1CE (Enhanced), or L2AP)
    !
    ! this subroutine reads each species independently of the others; 
    !  see subroutine turn_off_assim_SMAP_L1CTb() for avoiding the assimilation
    !  of redundant L1C_Tb observations when corresponding L2AP_Tb obs are assimilated
    !
    ! this subroutine is *not* meant to work for SMAP L3 files, but it could
    !  perhaps be extended to read soil moisture data from half-orbit h5 files
    !
    ! reichle, 17 Jan 2014
    ! reichle, 31 Jan 2014: added output of "SMAP_time"
    ! reichle, 31 May 2016: added stats check for L1C fore-minus-aft Tb differences
    ! reichle, 26 Dec 2017: added functionality for L1CE (Enhanced) files, incl. thinning
    ! reichle, 23 Jan 2018: removed stats check for L1C fore-minus-aft Tb differences;
    !                       use avg fore/aft timestamp so that fore and aft Tbs for same 
    !                         location are never used in different assimilation windows
    ! reichle, 22 Apr 2020: resurrected check for L1C fore-minus-aft Tb differences
    !                         after antenna-scan-angle (ASA) issues continued and the
    !                         SMAP Project declined to address these issues in L1 ops

    implicit none
        
    ! inputs:

    type(date_time_type), intent(in) :: date_time

    integer,              intent(in) :: N_catd, dtstep_assim

    type(obs_param_type), intent(in) :: this_obs_param

    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type),  intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input

    logical,              intent(in) :: write_obslog

    ! outputs:
    
    logical, intent(out)                    :: found_obs

    real,    intent(out), dimension(N_catd) :: SMAP_data, std_SMAP_data
    real,    intent(out), dimension(N_catd) :: SMAP_lon,  SMAP_lat
    real*8,  intent(out), dimension(N_catd) :: SMAP_time                 ! J2000 seconds

    ! --------------------------------------

    ! local variables
    
    character(4),      parameter :: J2000_epoch_id = 'TT12'    ! see date_time_util.F90
    
    ! SMAP data are stored in Yyyyy/Mmm/Ddd directories
    !
    ! Within each directory, there must be an ASCII text file that lists the h5 file 
    !  names of all files in the directory that are suitable for assimilation.
    ! These ASCII files must be named as follows:

    character(80), parameter :: fname_of_fname_list_L1C_TB_A   = 'SMAP_L1C_TB_A_list.txt'    
    character(80), parameter :: fname_of_fname_list_L1C_TB_D   = 'SMAP_L1C_TB_D_list.txt'    

    character(80), parameter :: fname_of_fname_list_L1C_TB_E_A = 'SMAP_L1C_TB_E_A_list.txt'    
    character(80), parameter :: fname_of_fname_list_L1C_TB_E_D = 'SMAP_L1C_TB_E_D_list.txt'    
    
    character(80), parameter :: fname_of_fname_list_L2_SM_AP_A = 'SMAP_L2_SM_AP_A_list.txt'
    character(80), parameter :: fname_of_fname_list_L2_SM_AP_D = 'SMAP_L2_SM_AP_D_list.txt'
    
    logical,   parameter :: tmp_debug    = .false.
    
    real,      parameter :: Tb_min       = 100.0  ! min allowed Tb
    real,      parameter :: Tb_max       = 320.0  ! max allowed Tb

    real,      parameter :: max_std_tb_fore_minus_aft = 20.  ! max std-dev L1C[E] fore-minus-aft Tb diffs

    integer,   parameter :: L1CE_spacing = 3  ! thinning of L1C_TB_E in units of 9-km indices ("3" => 27 km)
    
    ! temporarily shift lat/lon of obs for computation of nearest tile to
    ! avoid ambiguous assignment of M09 model tile within M36 obs grid cell
    ! (center of M36 grid cell is equidistant from at least two M09 model 
    !  tiles) -- reichle, 23 Aug 2013
    
    real,      parameter :: tmp_shift_lon    = 0.01
    real,      parameter :: tmp_shift_lat    = 0.005

    integer,   parameter :: dt_halforbit     = 50*60 ! seconds
    
    integer,   parameter :: N_halforbits_max = 15    ! max number of half-orbits per day

    integer,   parameter :: dtstep_assim_max = 10800 ! max allowed dtstep_assim [seconds]
    
    ! get max number of files to be read  (use 2*dt_halforbit just to be safe)
    
    integer,   parameter :: N_fnames_max =                                             &
         ceiling(real(N_halforbits_max)/86400.*real(dtstep_assim_max+2*dt_halforbit))

    ! --------------------------------------------

    type(hdf5read)       :: h5r

    logical              :: L1C_files, L1CE_files, hpol
    logical              :: L1CE_thinning, keep_data_1, keep_data_2, tmp_keep, file_exists

    type(date_time_type) :: date_time_low,       date_time_upp
    type(date_time_type) :: date_time_low_fname, date_time_tmp

    integer              :: ii, jj, kk, nn, mm
    integer              :: N_fnames, N_fnames_tmp, N_obs_tmp
    integer              :: dset_rank
    integer              :: ind_tile, ind_start, ind_end

    real                 :: M36_col_ind_tile, M36_row_ind_tile  
    real                 :: M36_col_ind_obs,  M36_row_ind_obs
    real                 :: tmpreal
    real                 :: tmp_tb_diff, tmp_tb_diff_Sum, tmp_tb_diff_SumOfSq, tmp_var

    real*8               :: J2000_seconds_low, J2000_seconds_upp

    character(  2)       :: orbit_tag
    character( 12)       :: tmpstr12
    character( 15)       :: SMAP_date_time
    character( 16)       :: YYYYMMDD_HHMMSSz 
    character( 80)       :: fname_of_fname_list, tmpstr80
    character(300)       :: fname_tmp, tmp_err_msg

    character(100)       :: dset_name_lon,    dset_name_lat
    character(100)       :: dset_name_col,    dset_name_row
    character(100)       :: dset_name_time_1, dset_name_tb_1, dset_name_tb_qual_flag_1
    character(100)       :: dset_name_time_2, dset_name_tb_2, dset_name_tb_qual_flag_2

    character(100), dimension(2*N_halforbits_max)  :: fname_list  ! max 2 days of files

    integer,        dimension(7)                   :: dset_size
    integer,        dimension(N_fnames_max)        :: N_obs_kept
    integer,        dimension(N_catd)              :: N_obs_in_tile

    real,           dimension(:),     allocatable  :: tmp_lon, tmp_lat

    integer,        dimension(:),     allocatable  :: tmp_col, tmp_row
    
    real,           dimension(:),     allocatable  :: tmp_tb_1
    real,           dimension(:),     allocatable  :: tmp_tb_2
    
    real*8,         dimension(:),     allocatable  :: tmp_time_1
    real*8,         dimension(:),     allocatable  :: tmp_time_2

    integer,        dimension(:),     allocatable  :: tmp_tb_qual_flag_1
    integer,        dimension(:),     allocatable  :: tmp_tb_qual_flag_2
    
    integer,        dimension(:),     allocatable  :: tmp_tile_num

    character(len=*), parameter :: Iam = 'read_obs_SMAP_halforbit_Tb'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------------
    
    ! check inputs
    
    ! the subroutine makes sense only if dtstep_assim <= 3 hours
    !
    ! (this avoids that more than 2 different Yyyyy/Mmm/Ddd directories are needed
    !  and that the time mismatch between the observed Tb and the model forecast Tb 
    !  becomes excessive;  in future, add time interpolation of forecast Tb)
    
    if (dtstep_assim > dtstep_assim_max) then
       err_msg = 'dtstep_assim must not exceed 3 hours'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    ! initialize
    
    found_obs = .false.
    
    ! read Tbs from L1C_TB, L1C_TB_E, or L2_SM_AP files? 

    L1CE_thinning = .false.   ! initialize

    if     (index(this_obs_param%descr,'L1C') /= 0) then
       
       if (index(this_obs_param%descr,'_E')  /= 0) then       
          
          L1CE_files = .true.                                 ! read SMAP L1C_TB_E (Enhanced) files
          L1C_files  = .false.
                    
          if (index(this_obs_param%descr,'_E27')  /= 0)  L1CE_thinning = .true.       
          
       else
          
          L1C_files  = .true.                                 ! read SMAP L1C_TB (standard) files
          L1CE_files = .false.       
          
       end if
       
    elseif (index(this_obs_param%descr,'L2AP') /= 0) then
       
       L1C_files  = .false.                                ! read SMAP L2_SM_AP files
       L1CE_files = .false.

    else

       err_msg = 'cannot interpret %descr'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if
 
        
    ! read ascending or descending files? 
    
    tmp_err_msg = 'inconsistent %descr and %orbit'
    
    if     (index(this_obs_param%descr,'_A') /=0 ) then
       
       orbit_tag = '_A'
       
       if (this_obs_param%orbit/=1) then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
       end if
       
    elseif (index(this_obs_param%descr,'_D') /=0 ) then
       
       orbit_tag = '_D'
       
       if (this_obs_param%orbit/=2) then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
       end if
       
    else
       
       err_msg = 'unknown %descr or %orbit'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if

    ! determine the name of the file that contains the relevant list of file names
    
    if     (this_obs_param%orbit==1) then       ! ascending
       
       if     (L1C_files)  then
          
          fname_of_fname_list = trim(fname_of_fname_list_L1C_TB_A)
          
       elseif (L1CE_files) then
          
          fname_of_fname_list = trim(fname_of_fname_list_L1C_TB_E_A)
          
       else
          
          fname_of_fname_list = trim(fname_of_fname_list_L2_SM_AP_A)
          
       end if
       
    elseif (this_obs_param%orbit==2) then       ! descending
       
       if     (L1C_files)  then
          
          fname_of_fname_list = trim(fname_of_fname_list_L1C_TB_D)
          
       elseif (L1CE_files) then
          
          fname_of_fname_list = trim(fname_of_fname_list_L1C_TB_E_D)
          
       else
          
          fname_of_fname_list = trim(fname_of_fname_list_L2_SM_AP_D)
          
       end if
       
    else
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown obs type')

    end if
        
    ! need h-pol or v-pol?
    
    if     (this_obs_param%pol==1) then
       
       hpol = .true.
       
    elseif (this_obs_param%pol==2) then
       
       hpol = .false.
       
    else
       
       err_msg = 'Problem with polarization'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if


    ! ---------------------------
    !
    ! define h5 data set names
    
    if (L1C_files .or. L1CE_files) then
       
       ! L1C_TB or L1_C_TB_E
       
       dset_name_lon               = '/Global_Projection/cell_lon'
       dset_name_lat               = '/Global_Projection/cell_lat'

       dset_name_col               = '/Global_Projection/cell_column'
       dset_name_row               = '/Global_Projection/cell_row'

       dset_name_time_1            = '/Global_Projection/cell_tb_time_seconds_fore'
       dset_name_time_2            = '/Global_Projection/cell_tb_time_seconds_aft'

       if (hpol) then
          
          dset_name_tb_1           = '/Global_Projection/cell_tb_h_fore'
          dset_name_tb_2           = '/Global_Projection/cell_tb_h_aft'
          
          dset_name_tb_qual_flag_1 = '/Global_Projection/cell_tb_qual_flag_h_fore'
          dset_name_tb_qual_flag_2 = '/Global_Projection/cell_tb_qual_flag_h_aft'
          
       else
          
          dset_name_tb_1           = '/Global_Projection/cell_tb_v_fore'
          dset_name_tb_2           = '/Global_Projection/cell_tb_v_aft'
          
          dset_name_tb_qual_flag_1 = '/Global_Projection/cell_tb_qual_flag_v_fore'
          dset_name_tb_qual_flag_2 = '/Global_Projection/cell_tb_qual_flag_v_aft'
          
       end if
       
    else  
       
       ! L2_SM_AP
       
       dset_name_lon               = '/Soil_Moisture_Retrieval_Data/longitude'
       dset_name_lat               = '/Soil_Moisture_Retrieval_Data/latitude'
       
       dset_name_time_1            = '/Soil_Moisture_Retrieval_Data/spacecraft_overpass_time_seconds'
       dset_name_time_2            = ''  ! *not* used

       if (hpol) then
          
          dset_name_tb_1           = '/Soil_Moisture_Retrieval_Data/tb_h_disaggregated'
          dset_name_tb_2           = ''  ! *not* used
          
          dset_name_tb_qual_flag_1 = '/Soil_Moisture_Retrieval_Data/tb_h_disaggregated_qual_flag'
          dset_name_tb_qual_flag_2 = ''  ! *not* used
          
       else
          
          dset_name_tb_1           = '/Soil_Moisture_Retrieval_Data/tb_v_disaggregated'
          dset_name_tb_2           = ''  ! *not* used
          
          dset_name_tb_qual_flag_1 = '/Soil_Moisture_Retrieval_Data/tb_v_disaggregated_qual_flag'
          dset_name_tb_qual_flag_2 = ''  ! *not* used
          
       end if
       
    end if
    
    ! ---------------------------
    !
    ! Define search interval for obs
    !
    !   [ date_time - dtstep_assim/2 + one_second, date_time + dtstep_assim/2 ]

    ! lower boundary
    
    date_time_low = date_time
    
    call augment_date_time( -dtstep_assim/2+1, date_time_low )

    J2000_seconds_low = datetime_to_J2000seconds( date_time_low, J2000_epoch_id )

    ! upper boundary
    
    date_time_upp = date_time

    call augment_date_time(  dtstep_assim/2,   date_time_upp )

    J2000_seconds_upp = datetime_to_J2000seconds( date_time_upp, J2000_epoch_id )

    ! -----------------------------------------------------------------
    !
    ! identify names of files that could contain data within the
    !  assimilation window [date_time_low,date_time_upp]

    date_time_low_fname = date_time_low
    
    call augment_date_time( -dt_halforbit, date_time_low_fname )
        
    ! read file with list of SMAP file names for first day
    
    call read_obs_SMAP_fnames( date_time_low_fname, this_obs_param,                &
         fname_of_fname_list, N_halforbits_max,                                    &
         N_fnames, fname_list(1:N_halforbits_max) )
    
    ! if needed, read file with list of SMAP file names for second day and add
    !  file names into "fname_list"
    
    if (date_time_low_fname%day /= date_time_upp%day) then
       
       call read_obs_SMAP_fnames( date_time_upp, this_obs_param,                   &
            fname_of_fname_list, N_halforbits_max,                                 &
            N_fnames_tmp, fname_list((N_fnames+1):(N_fnames+N_halforbits_max)) )
       
       N_fnames = N_fnames + N_fnames_tmp
       
    end if

    ! ------------------------------------------------------------------
    
    ! extract names of files that could contain data within the assimilation
    !  window
    !
    ! sample file names:
    !                                     ||||||||||||||| 
    !   Yyyyy/Mmm/Ddd/SMAP_L1C_TB_03027_D_20010727T160914_D04003_000.h5
    !
    !                                       ||||||||||||||| 
    !   Yyyyy/Mmm/Ddd/SMAP_L1C_TB_E_03027_D_20010727T160914_D04003_000.h5
    !   Yyyyy/Mmm/Ddd/SMAP_L2_SM_AP_03073_D_20010730T193828_D04003_000.h5
    !                                       ||||||||||||||| 
    ! counter:   1         2         3         4         5         6
    !   1234567890123456789012345678901234567890123456789012345678901234567890

    if (L1C_files) then
       
       ind_start = 35
       ind_end   = 49

    else  ! L1CE or L2AP files
       
       ind_start = 37
       ind_end   = 51
       
    end if

    kk = 0
    
    do ii=1,N_fnames

       SMAP_date_time = fname_list(ii)(ind_start:ind_end)
       
       date_time_tmp = SMAPdatetime_to_DateTimeType( SMAP_date_time )

       ! check whether:  date_time_low_fname < date_time_tmp <= date_time_upp

       if ( datetime_lt_refdatetime( date_time_low_fname, date_time_tmp ) .and.        &
            datetime_le_refdatetime( date_time_tmp,       date_time_upp )       ) then 
          
          kk = kk+1

          ! there can be no more than N_fnames_max files that have data falling into the 
          ! assimilation window (see also dtstep_assim_max)

          if (kk>N_fnames_max) then
             err_msg = 'too many files  match assimilation window'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if

          fname_list(kk) = fname_list(ii)
          
       end if
       
    end do  ! ii=1,N_fnames
    
    N_fnames = kk

    ! the "N_fnames" files of interest are now the first "N_fnames" entries
    !  in "fname_list"
    
    ! ---------------------------------------------------------
    !
    ! read and process data if files were found

    if (N_fnames==0) then
       
       ! no data files found
       
       SMAP_data     =      this_obs_param%nodata
       SMAP_lon      =      this_obs_param%nodata
       SMAP_lat      =      this_obs_param%nodata
       SMAP_time     = real(this_obs_param%nodata,kind(0.0D0))
       std_SMAP_data =      this_obs_param%nodata
       
    else
       
       ! initialize outputs
       
       SMAP_data = 0.
       SMAP_lon  = 0.
       SMAP_lat  = 0.
       SMAP_time = 0.0D0

       N_obs_in_tile  = 0  ! for normalization after mapping to tile and super-obs
       
       ! loop through files
              
       do kk=1,N_fnames
          
          ! open file
          
          fname_tmp = trim(this_obs_param%path) // '/' // fname_list(kk)
          
          if (logit) write(logunit,'(400A)') 'reading file: ' // trim(fname_tmp)
          
          inquire(file=fname_tmp, exist=file_exists)
          
          if (.not. file_exists) then
             
             err_msg = 'file does NOT exist'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)             

          end if

          call h5r%openFile(fname_tmp)
          
          ! ------------------------------

          ! read h5 datasets

          tmp_err_msg = 'read_obs_SMAP_halforbit_Tb(): inconsistent dataset lengths'

          ! LONGITUDE: query dataset, record size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_lon)

          call h5r%queryDataset(dset_name_lon, dset_rank, dset_size)
          
          N_obs_tmp = dset_size(1)
          
          allocate(tmp_lon(N_obs_tmp))
          
          call h5r%readDataset(tmp_lon)
          
          ! LATITUDE: query dataset, check size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_lat)

          call h5r%queryDataset(dset_name_lat, dset_rank, dset_size)
          
          if (N_obs_tmp/=dset_size(1)) then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
          end if
          
          allocate(tmp_lat(N_obs_tmp))

          call h5r%readDataset(tmp_lat)

          if (L1CE_thinning) then  ! need to read column and row indices

             ! COLUMN: query dataset, record size, allocate space, read data
             
             if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_col)
             
             call h5r%queryDataset(dset_name_col, dset_rank, dset_size)
             
             if (N_obs_tmp/=dset_size(1)) then
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
             end if
             
             allocate(tmp_col(N_obs_tmp))
             
             call h5r%readDataset(tmp_col)
             
             ! ROW: query dataset, check size, allocate space, read data
             
             if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_row)
             
             call h5r%queryDataset(dset_name_row, dset_rank, dset_size)
             
             if (N_obs_tmp/=dset_size(1)) then
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
             end if
             
             allocate(tmp_row(N_obs_tmp))
             
             call h5r%readDataset(tmp_row)

          end if
                     
          ! TIME_1: query dataset, check size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_time_1)
          
          call h5r%queryDataset(dset_name_time_1, dset_rank, dset_size)
          
          if (N_obs_tmp/=dset_size(1)) then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
          end if
          
          allocate(tmp_time_1(N_obs_tmp))

          call h5r%readDataset(tmp_time_1)

          ! TB_1: query dataset, check size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_tb_1)
          
          call h5r%queryDataset(dset_name_tb_1, dset_rank, dset_size)
          
          if (N_obs_tmp/=dset_size(1)) then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
          end if
          
          allocate(tmp_tb_1(N_obs_tmp))

          call h5r%readDataset(tmp_tb_1)

          ! TB_QUAL_FLAG_1: query dataset, check size, allocate space, read data
          
          if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_tb_qual_flag_1)

          call h5r%queryDataset(dset_name_tb_qual_flag_1, dset_rank, dset_size)
          
          if (N_obs_tmp/=dset_size(1)) then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
          end if
          
          allocate(tmp_tb_qual_flag_1(N_obs_tmp))

          call h5r%readDataset(tmp_tb_qual_flag_1)

          ! for L1C_TB or L1C_TB_E files also read "aft"

          if (L1C_files .or. L1CE_files) then
             
             ! *_1 = "fore"
             ! *_2 = "aft"

             ! TIME_2: query dataset, check size, allocate space, read data
             
             if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_time_2)

             call h5r%queryDataset(dset_name_time_2, dset_rank, dset_size)
             
             if (N_obs_tmp/=dset_size(1)) then
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
             end if
             
             allocate(tmp_time_2(N_obs_tmp))

             call h5r%readDataset(tmp_time_2)

             ! TB_2: query dataset, check size, allocate space, read data
             
             if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_tb_2)
             
             call h5r%queryDataset(dset_name_tb_2, dset_rank, dset_size)
             
             if (N_obs_tmp/=dset_size(1)) then
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
             end if

             allocate(tmp_tb_2(N_obs_tmp))

             call h5r%readDataset(tmp_tb_2)
             
             ! TB_QUAL_FLAG_2: query dataset, check size, allocate space, read data
             
             if (tmp_debug .and. logit) write(logunit,*) trim(dset_name_tb_qual_flag_2)
             
             call h5r%queryDataset(dset_name_tb_qual_flag_2, dset_rank, dset_size)
             
             if (N_obs_tmp/=dset_size(1)) then
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, tmp_err_msg)
             end if
             
             allocate(tmp_tb_qual_flag_2(N_obs_tmp))
             
             call h5r%readDataset(tmp_tb_qual_flag_2)

          end if
          
          ! close file
          
          call h5r%closeFile

          if (logit) write(logunit,*) 'done reading file'

          ! --------------------------------------------------------
          !
          ! eliminate obs outside desired time window, no-data-values and data 
          !  that fail initial QC
          ! keep track of how many obs survived from current file

          ! #####################################
          ! TO DO: 
          ! - refine use of quality flag?
          ! - how to QC for proximity to water body? 
          ! - anything else missing that is done for SMOS 
          !        (either in preproc or in obs reader)??
          ! #####################################
          
          jj = 0

          if (L1C_files .or. L1CE_files) then

             ! initialize stats check for fore-minus-aft Tb diffs

             mm = 0
             
             tmp_tb_diff_Sum     = 0.
             tmp_tb_diff_SumOfSq = 0.

             do nn=1,N_obs_tmp

                ! QC

                keep_data_1 =                                    &
                     (mod(tmp_tb_qual_flag_1(nn),2)==0)    .and. & ! lowest bit must be 0
                     (tmp_tb_1(nn)   >  Tb_min)            .and. & ! elim neg nodata
                     (tmp_tb_1(nn)   <  Tb_max)                    ! elim huge pos nodata
                
                keep_data_2 =                                    &
                     (mod(tmp_tb_qual_flag_2(nn),2)==0)    .and. & ! lowest bit must be 0
                     (tmp_tb_2(nn)   >  Tb_min)            .and. & ! elim neg nodata
                     (tmp_tb_2(nn)   <  Tb_max)                    ! elim huge pos nodata

                ! thinning of L1C_TB_E obs
                
                if (L1CE_thinning) then
                   
                   tmp_keep =                                          &
                        ( mod(tmp_col(nn),L1CE_spacing) == 1 )  .and.  &
                        ( mod(tmp_row(nn),L1CE_spacing) == 1 )        
                   
                   keep_data_1 = keep_data_1 .and. tmp_keep
                   keep_data_2 = keep_data_2 .and. tmp_keep
                   
                end if

                ! compute fore and aft average, put into "tb_1", "tmp_time_1"
                
                if     (keep_data_1 .and. keep_data_2) then
                   
                   ! Compute stats for fore-minus-aft Tb differences.
                   ! Excessive diffs are found in bad L1C_TB files, which occur 
                   ! occasionally due to bad ANT_AZ files in L1B processing.
                   ! Includes ALL SURFACES!!!
                   ! - reichle, 22 Apr 2020 (resurrected)
                   ! - reichle, 16 Oct 2020 (bug fix: do stats first, then avg)
                   
                   mm=mm+1
                   
                   tmp_tb_diff         = tmp_tb_1(nn) - tmp_tb_2(nn) 
                   
                   tmp_tb_diff_Sum     = tmp_tb_diff_Sum     + tmp_tb_diff
                   tmp_tb_diff_SumOfSq = tmp_tb_diff_SumOfSq + tmp_tb_diff**2
                   
                   ! put average of "fore" and "aft" into "tb_1", "tmp_time_1"
                   
                   tmp_tb_1(  nn) = 0.5  *( tmp_tb_1(  nn) + tmp_tb_2(  nn) )
                   tmp_time_1(nn) = 0.5D0*( tmp_time_1(nn) + tmp_time_2(nn) ) 
                   
                elseif (keep_data_2)                   then

                   ! put "aft" data into "tb_1", "tmp_time_1"
                   
                   tmp_tb_1(  nn) = tmp_tb_2(  nn)
                   tmp_time_1(nn) = tmp_time_2(nn)

                else

                   ! nothing to do here
                   ! - if only keep_data_1 is true  (tmp_tb_1 and tmp_time_1 are already as needed)
                   ! - if both keep_data_1 and keep_data_2 are false  (next if block ignores data)
                   
                end if

                ! apply QC and thinning, ensure that time stamp is within assimilation window
                
                if ( (keep_data_1 .or. keep_data_2)        .and.           & 
                     (tmp_time_1(nn) >  J2000_seconds_low) .and.           &
                     (tmp_time_1(nn) <= J2000_seconds_upp)                 &
                     )                                             then
                   
                   jj=jj+1
                   
                   tmp_lon(   jj) = tmp_lon(nn)
                   tmp_lat(   jj) = tmp_lat(nn)
                   tmp_tb_1(  jj) = tmp_tb_1(  nn)
                   tmp_time_1(jj) = tmp_time_1(nn)
                   
                end if
                
             end do
             
             ! finalize stats check for fore-minus-aft differences (ALL SURFACES!!!)
             ! - reichle, 22 Apr 2020 (resurrected)
             
             if (mm>1) then
                
                tmp_var = ( tmp_tb_diff_SumOfSq - (tmp_tb_diff_Sum**2)/real(mm) )/(real(mm-1))
                
                if ( tmp_var > max_std_tb_fore_minus_aft**2 ) then
                   
                   write(err_msg, '(e12.5)') sqrt(tmp_var) 
                   
                   err_msg =                                                                                      &
                        'Ignoring ALL obs in halforbit file b/c of excessive std-dev in fore-minus-aft Tbs. ' //  &
                        'std-dev( tb_fore - tb_aft ) = ' // trim(err_msg)
                   
                   call ldas_warn(LDAS_GENERIC_WARNING, Iam, err_msg)
                   
                   jj = 0  ! results in N_obs_kept=0 below
                   
                end if

             end if

          else  ! L2_SM_AP
             
             do nn=1,N_obs_tmp
                
                keep_data_1 =                                    &
                     (mod(tmp_tb_qual_flag_1(nn),2)==0)    .and. & ! lowest bit must be 0
                     (tmp_time_1(nn) >  J2000_seconds_low) .and. &
                     (tmp_time_1(nn) <= J2000_seconds_upp) .and. &
                     (tmp_tb_1(nn)   >  Tb_min)            .and. & ! elim neg nodata
                     (tmp_tb_1(nn)   <  Tb_max)                    ! elim huge pos nodata
                
                if (keep_data_1) then
                   
                   jj=jj+1
                   
                   tmp_lon(   jj) = tmp_lon(   nn)
                   tmp_lat(   jj) = tmp_lat(   nn)
                   tmp_tb_1(  jj) = tmp_tb_1(  nn)
                   tmp_time_1(jj) = tmp_time_1(nn)
                   
                end if
                
             end do
             
          end if  ! (L1C_files .or. L1CE_files)
          
          N_obs_kept(kk) = jj
             
          ! -------------------------------------------------
          !
          ! map obs to tiles
          !
          ! for each observation
          !     a) determine grid cell that contains lat/lon
          !     b) determine tile within grid cell that contains lat/lon
          
          if (N_obs_kept(kk)>0) then
             
             allocate(tmp_tile_num(N_obs_kept(kk)))
             
             ! shift M36 obs lat/lon for proper assignment of M09 tile?
             
             if ( L1C_files .and. (index(tile_grid_d%gridtype, 'EASEv2_M09')  /=0 .or. index(tile_grid_d%gridtype, 'EASEv2-M09')  /=0 )) then
                
                ! temporarily shift lat/lon of obs for computation of nearest tile to
                ! avoid ambiguous assignment of M09 model tile within M36 obs grid cell
                ! (center of M36 grid cell is equidistant from at least two M09 model 
                !  tiles) -- reichle, 23 Aug 2013
                
                call get_tile_num_for_obs(N_catd, tile_coord,                 &
                     tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
                     N_obs_kept(kk),                                          &
                     tmp_lat(1:N_obs_kept(kk)),                               &
                     tmp_lon(1:N_obs_kept(kk)),                               &
                     this_obs_param,                                          &
                     tmp_tile_num,                                            &
                     tmp_shift_lat, tmp_shift_lon )
                
             else
                
                call get_tile_num_for_obs(N_catd, tile_coord,                 &
                     tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,     &
                     N_obs_kept(kk),                                          &
                     tmp_lat(1:N_obs_kept(kk)),                               &
                     tmp_lon(1:N_obs_kept(kk)),                               &
                     this_obs_param,                                          &
                     tmp_tile_num)
                
             end if
             
             ! make sure center-of-mass of tile that administers obs 
             !  is within EASEv2 M36 obs grid cell, discard obs otherwise
             !  (by setting tmp_tile_num to negative value)
             !
             ! this step eliminates SMAP obs that fall outside of the GEOS-5 
             ! land mask at M09 resolution
             !
             ! not sure what this does if model is run on tile space other than 
             ! EASEv2_M09 - reichle, 11 Nov 2014
             !
             ! Bug fix - index used within loop should be "ii" (not "kk").
             ! For global runs the bug should have had only minimal impact.
             ! - reichle, 23 Apr 2015
             !
             ! It is not 100 percent clear why this piece code had been added
             !  at the time.
             ! Chances are that it had to do with land-water mask issues and the
             !  distortion of EASE grid cells at high latitudes, and/or Gabrielle's
             !  use of an M36 innovations integration to derive Tb scaling files
             !  for the M09 (SMAP) system.
             ! The problem with the piece of code is that it throws out far too many 
             !  obs (at all latitudes) if the tile space is coarse, e.g., that of the 
             !  1/2 deg Lat/Lon grids of MERRA or MERRA-2.
             ! The piece of code may no longer be needed because of improvements in
             !  the obs readers and in get_obs_pred(), but the impact of removing the
             !  code on the SMAP L4_SM system is not clear.  At this time, just prior
             !  to finalizing the L4_SM "validated release", keep the code for the 
             !  EASEv2 M09 and M36 tile spaces that are relevant for the L4_SM 
             !  system, but drop it for all other tile spaces.
             ! - reichle, 3 Feb 2016
             
             if   (                                                               &
                  (index(tile_grid_d%gridtype, 'EASEv2_M09')  /=0) .or.           &
                  (index(tile_grid_d%gridtype, 'EASEv2_M36')  /=0)       )  then
                
                do ii=1,N_obs_kept(kk)
                   
                   if (tmp_tile_num(ii)>0) then
                      
                      call ease_convert('EASEv2_M36',                  &
                           tile_coord(tmp_tile_num(ii))%com_lat,  &
                           tile_coord(tmp_tile_num(ii))%com_lon,  &
                           M36_col_ind_tile, M36_row_ind_tile  )
                      
                      call ease_convert('EASEv2_M36',                  &
                           tmp_lat(ii),                           &
                           tmp_lon(ii),                           &
                           M36_col_ind_obs,  M36_row_ind_obs  )
                      
                      if ( (nint(M36_col_ind_tile)/=nint(M36_col_ind_obs))  .or.     &
                           (nint(M36_row_ind_tile)/=nint(M36_row_ind_obs))        )  &
                           tmp_tile_num(ii) = -9999
                      
                   end if
                   
                end do
                
             end if

             ! compute super-obs for each tile from all obs w/in that tile
             !     (also eliminate observations that are not in domain)
             
             do ii=1,N_obs_kept(kk)
                
                ind_tile = tmp_tile_num(ii)   ! 1<=tmp_tile_num<=N_catd (unless nodata)
                
                if (ind_tile>0) then          ! this step eliminates obs outside domain
                
                   SMAP_data(ind_tile) = SMAP_data(ind_tile) + tmp_tb_1(  ii)
                   SMAP_lon( ind_tile) = SMAP_lon( ind_tile) + tmp_lon(   ii)
                   SMAP_lat( ind_tile) = SMAP_lat( ind_tile) + tmp_lat(   ii)
                   SMAP_time(ind_tile) = SMAP_time(ind_tile) + tmp_time_1(ii)
                   
                   N_obs_in_tile(ind_tile) = N_obs_in_tile(ind_tile) + 1
                   
                end if
             
             end do

             deallocate(tmp_tile_num)
             
          end if
          
          ! clean up
          
          if (allocated(tmp_lon           )) deallocate(tmp_lon           )
          if (allocated(tmp_lat           )) deallocate(tmp_lat           )
          if (allocated(tmp_col           )) deallocate(tmp_col           )
          if (allocated(tmp_row           )) deallocate(tmp_row           )
          if (allocated(tmp_time_1        )) deallocate(tmp_time_1        )
          if (allocated(tmp_time_2        )) deallocate(tmp_time_2        )
          if (allocated(tmp_tb_1          )) deallocate(tmp_tb_1          )
          if (allocated(tmp_tb_2          )) deallocate(tmp_tb_2          )
          if (allocated(tmp_tb_qual_flag_1)) deallocate(tmp_tb_qual_flag_1)
          if (allocated(tmp_tb_qual_flag_2)) deallocate(tmp_tb_qual_flag_2)
          
       end do  ! kk=1,N_fnames
          
       ! normalize
       
       do  ii=1,N_catd
          
          if (N_obs_in_tile(ii)>1) then
             
             tmpreal = real(N_obs_in_tile(ii))
             
             SMAP_data(ii) = SMAP_data(ii)/     tmpreal
             SMAP_lon( ii) = SMAP_lon( ii)/     tmpreal
             SMAP_lat( ii) = SMAP_lat( ii)/     tmpreal
             SMAP_time(ii) = SMAP_time(ii)/real(tmpreal,kind(0.0D0))
             
          elseif (N_obs_in_tile(ii)==0) then
             
             SMAP_data(ii) =      this_obs_param%nodata
             SMAP_lon( ii) =      this_obs_param%nodata
             SMAP_lat( ii) =      this_obs_param%nodata
             SMAP_time(ii) = real(this_obs_param%nodata,kind(0.0D0))
             
          end if
          
       end do
       
       ! --------------------------------
       
       ! set observation error standard deviation
       
       do ii=1,N_catd
          std_SMAP_data(ii) = this_obs_param%errstd
       end do
       
       ! --------------------------------
       
       if (any(N_obs_in_tile>0)) then
          
          found_obs = .true.
          
       else 
          
          found_obs = .false.
          
       end if
       
    end if  ! N_fnames==0
    
    ! -------------------------------------------------
    !
    ! write "obslog" file

    if (write_obslog) then

       YYYYMMDD_HHMMSSz = date_time2string(date_time)

       tmpstr80 = 'read_obs_SMAP_halforbit_Tb()'  ! name of this subroutine
              
       do kk=1,N_fnames
          
          fname_tmp = trim(this_obs_param%path) // fname_list(kk)
          
          write (tmpstr12,'(i12)') N_obs_kept(kk)    ! convert integer to string
          
          call add_to_obslog( YYYYMMDD_HHMMSSz, this_obs_param%descr, tmpstr80, &
               tmpstr12, fname_tmp )
          
       end do

    end if

    ! clean up
    
    if (logit) write (logunit,*) 'read_obs_SMAP_halforbit_Tb(): done.'
    
  end subroutine read_obs_SMAP_halforbit_Tb
  
  ! *****************************************************************

  subroutine turn_off_assim_SMAP_L1CTb(N_obs_param, obs_param, N_obsl, Observations_l)
    
    ! this subroutine turns off the assimilation of *individual* L1C_Tb obs 
    !  if corresponding L2AP_Tb obs are assimilated 
    
    ! rationale: SMAP L1C_Tb obs are used to derive the disaggregated L2AP_Tb obs
    !  and should not be assimilated along with L2AP_Tb
    
    ! L1C_Tb  obs are on the EASEv2 36 km grid ("M36")
    ! L2AP_Tb obs are on the EASEv2  9 km grid ("M09")

    ! reichle, 17 Jan 2014
    ! reichle,  5 May 2015 - added *ascending* L2AP
 
    implicit none
    
    integer,                                      intent(in)    :: N_obs_param, N_obsl
    
    type(obs_param_type), dimension(N_obs_param), intent(in)    :: obs_param  
    
    type(obs_type),       dimension(N_obsl),      intent(inout) :: Observations_l  
        
    ! local variables 

    logical            :: turnoff_L1C_Tbh_A, turnoff_L1C_Tbv_A 
    logical            :: turnoff_L1C_Tbh_D, turnoff_L1C_Tbv_D 

    integer            :: n, ii, j, N_obsf, N_cols, N_rows

    integer            :: species_L1C_Tbh_A, species_L2AP_Tbh_A
    integer            :: species_L1C_Tbh_D, species_L2AP_Tbh_D
    integer            :: species_L1C_Tbv_A, species_L2AP_Tbv_A
    integer            :: species_L1C_Tbv_D, species_L2AP_Tbv_D

    integer            :: ind_L1C_Tbh_A, ind_L2AP_Tbh_A
    integer            :: ind_L1C_Tbh_D, ind_L2AP_Tbh_D
    integer            :: ind_L1C_Tbv_A, ind_L2AP_Tbv_A
    integer            :: ind_L1C_Tbv_D, ind_L2AP_Tbv_D

    real               :: col, row
    
    integer,        dimension(numprocs)              :: N_obsl_vec, tmp_low_ind
    
    type(obs_type), dimension(:),        allocatable :: Observations_f
    
    logical,        dimension(:,:),      allocatable :: mask_h_A, mask_v_A
    logical,        dimension(:,:),      allocatable :: mask_h_D, mask_v_D

    character(len=*), parameter :: Iam = 'turn_off_assim_SMAP_L1CTb'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------------------------    
    
    ! determine species and index numbers of conflicting obs species
    
    species_L1C_Tbh_A  = -9999
    species_L1C_Tbh_D  = -9999
    species_L1C_Tbv_A  = -9999
    species_L1C_Tbv_D  = -9999
    species_L2AP_Tbh_A = -9999
    species_L2AP_Tbh_D = -9999
    species_L2AP_Tbv_A = -9999
    species_L2AP_Tbv_D = -9999

    ind_L1C_Tbh_A      = -9999
    ind_L1C_Tbh_D      = -9999
    ind_L1C_Tbv_A      = -9999
    ind_L1C_Tbv_D      = -9999
    ind_L2AP_Tbh_A     = -9999
    ind_L2AP_Tbh_D     = -9999
    ind_L2AP_Tbv_A     = -9999
    ind_L2AP_Tbv_D     = -9999
    
    do j=1,N_obs_param

       ! reichle, 1 Feb 2018: Only worry about assimilated species.
       !                      NOTE: At most one of L1C_Tb, L1C_Tb E09, or L1C_Tb E27 can be 
       !                            assimilated at a time, checked in read_ens_upd_inputs().

       if (obs_param(j)%assim) then     
          
          select case (trim(obs_param(j)%descr))
             
          case('SMAP_L1C_Tbh_A','SMAP_L1C_Tbh_E09_A','SMAP_L1C_Tbh_E27_A'); species_L1C_Tbh_A =obs_param(j)%species; ind_L1C_Tbh_A =j 
          case('SMAP_L1C_Tbh_D','SMAP_L1C_Tbh_E09_D','SMAP_L1C_Tbh_E27_D'); species_L1C_Tbh_D =obs_param(j)%species; ind_L1C_Tbh_D =j 
          case('SMAP_L1C_Tbv_A','SMAP_L1C_Tbv_E09_A','SMAP_L1C_Tbv_E27_A'); species_L1C_Tbv_A =obs_param(j)%species; ind_L1C_Tbv_A =j 
          case('SMAP_L1C_Tbv_D','SMAP_L1C_Tbv_E09_D','SMAP_L1C_Tbv_E27_D'); species_L1C_Tbv_D =obs_param(j)%species; ind_L1C_Tbv_D =j 
          case('SMAP_L2AP_Tbh_A');                                          species_L2AP_Tbh_A=obs_param(j)%species; ind_L2AP_Tbh_A=j 
          case('SMAP_L2AP_Tbh_D');                                          species_L2AP_Tbh_D=obs_param(j)%species; ind_L2AP_Tbh_D=j 
          case('SMAP_L2AP_Tbv_A');                                          species_L2AP_Tbv_A=obs_param(j)%species; ind_L2AP_Tbv_A=j 
          case('SMAP_L2AP_Tbv_D');                                          species_L2AP_Tbv_D=obs_param(j)%species; ind_L2AP_Tbv_D=j 

          end select
          
       end if

    end do
    
    ! determine whether there is possibly a redundancy

    turnoff_L1C_Tbh_A = .false.
    turnoff_L1C_Tbh_D = .false.
    turnoff_L1C_Tbv_A = .false.
    turnoff_L1C_Tbv_D = .false.
    
    if (ind_L1C_Tbh_A>0 .and. ind_L2AP_Tbh_A>0) turnoff_L1C_Tbh_A = .true.
    if (ind_L1C_Tbh_D>0 .and. ind_L2AP_Tbh_D>0) turnoff_L1C_Tbh_D = .true.
    if (ind_L1C_Tbv_A>0 .and. ind_L2AP_Tbv_A>0) turnoff_L1C_Tbv_A = .true.
    if (ind_L1C_Tbv_D>0 .and. ind_L2AP_Tbv_D>0) turnoff_L1C_Tbv_D = .true.

    ! -----------------------------------------------------------------------

    ! proceed only if there is work to be done
    
    if ( turnoff_L1C_Tbh_A .or. turnoff_L1C_Tbv_A .or.                &
         turnoff_L1C_Tbh_D .or. turnoff_L1C_Tbv_D        ) then
       
       ! gather local obs 
       
#ifdef LDAS_MPI
       
       call MPI_Gather(                 &
            N_obsl,     1, MPI_integer, &
            N_obsl_vec, 1, MPI_integer, & 
            0, mpicomm, mpierr ) 
       
#else
       
       N_obsl_vec(1) = N_obsl
       
#endif
       
       if (root_proc) then
          
          N_obsf = sum(N_obsl_vec)
          
          allocate(Observations_f(N_obsf))
          
          tmp_low_ind(1) = 1
          
          do n=1,numprocs-1
             
             tmp_low_ind(n+1) = tmp_low_ind(n) + N_obsl_vec(n)
             
          end do
          
       end if
       
#ifdef LDAS_MPI
       
       call MPI_GATHERV(                                                    &
            Observations_l,  N_obsl,                      MPI_obs_type,     &
            Observations_f,  N_obsl_vec,  tmp_low_ind-1,  MPI_obs_type,     &
            0, mpicomm, mpierr ) 
       
#else             
       
       Observations_f = Observations_l
       
#endif
       
       ! ---------------------------------------------------------    
       !
       ! assemble 36 km EASEv2 mask of L2AP_Tb obs 

       call ease_extent( 'EASEv2_M36', N_cols, N_rows )
       
       allocate( mask_h_A(N_cols,N_rows) )
       allocate( mask_h_D(N_cols,N_rows) )
       allocate( mask_v_A(N_cols,N_rows) )
       allocate( mask_v_D(N_cols,N_rows) )
       
       mask_h_A = .false.   ! initialize
       mask_h_D = .false.   ! initialize
       mask_v_A = .false.   ! initialize
       mask_v_D = .false.   ! initialize
    
       if (root_proc) then
          
          ! mask for H-pol ascending
          
          if (turnoff_L1C_Tbh_A) then
             
             do ii=1,N_obsf
                
                if  (Observations_f(ii)%species==species_L2AP_Tbh_A) then  
                   
                   call ease_convert('EASEv2_M36', Observations_f(ii)%lat, Observations_f(ii)%lon, &
                        col, row)

                   ! set mask=.true. for the M36 grid cell that contains the L2AP_Tb obs;
                   ! note conversion to one-based indices
                   
                   mask_h_A(nint(col)+1,nint(row)+1) = .true.
                   
                end if
                
             end do

          end if

          ! mask for H-pol descending
          
          if (turnoff_L1C_Tbh_D) then
             
             do ii=1,N_obsf
                
                if  (Observations_f(ii)%species==species_L2AP_Tbh_D) then  
                   
                   call ease_convert('EASEv2_M36', Observations_f(ii)%lat, Observations_f(ii)%lon, &
                        col, row)

                   ! set mask=.true. for the M36 grid cell that contains the L2AP_Tb obs;
                   ! note conversion to one-based indices
                   
                   mask_h_D(nint(col)+1,nint(row)+1) = .true.
                   
                end if
                
             end do

          end if
          
          ! mask for V-pol ascending
          
          if (turnoff_L1C_Tbv_A) then
             
             do ii=1,N_obsf
                
                if  (Observations_f(ii)%species==species_L2AP_Tbv_A) then  
                   
                   call ease_convert('EASEv2_M36', Observations_f(ii)%lat, Observations_f(ii)%lon, &
                        col, row)

                   ! set mask=.true. for the M36 grid cell that contains the L2AP_Tb obs;
                   ! note conversion to one-based indices
                   
                   mask_v_A(nint(col)+1,nint(row)+1) = .true.
                   
                end if
                
             end do
             
          end if

          ! mask for V-pol descending
          
          if (turnoff_L1C_Tbv_D) then
             
             do ii=1,N_obsf
                
                if  (Observations_f(ii)%species==species_L2AP_Tbv_D) then  
                   
                   call ease_convert('EASEv2_M36', Observations_f(ii)%lat, Observations_f(ii)%lon, &
                        col, row)

                   ! set mask=.true. for the M36 grid cell that contains the L2AP_Tb obs;
                   ! note conversion to one-based indices
                   
                   mask_v_D(nint(col)+1,nint(row)+1) = .true.
                   
                end if
                
             end do
             
          end if

          deallocate(Observations_f)
          
       end if  ! (root_proc)
       
       ! MPI broadcast masks
       
       call MPI_Bcast(mask_h_A, N_cols*N_rows, MPI_LOGICAL, 0, mpicomm, mpierr)
       call MPI_Bcast(mask_h_D, N_cols*N_rows, MPI_LOGICAL, 0, mpicomm, mpierr)
       call MPI_Bcast(mask_v_A, N_cols*N_rows, MPI_LOGICAL, 0, mpicomm, mpierr)
       call MPI_Bcast(mask_v_D, N_cols*N_rows, MPI_LOGICAL, 0, mpicomm, mpierr)
       
       ! ---------------------------------------------------------    
       !
       ! apply H-pol masks
       
       if (turnoff_L1C_Tbh_A) then
          
          do ii=1,N_obsl
             
             if (Observations_l(ii)%species==species_L1C_Tbh_A) then  
                
                call ease_convert('EASEv2_M36',                                      &
                     Observations_l(ii)%lat, Observations_l(ii)%lon, col, row)
                
                ! note conversion to one-based indices
                
                if (mask_h_A(nint(col)+1,nint(row)+1))  Observations_l(ii)%assim = .false.

             end if

          end do

       end if

       if (turnoff_L1C_Tbh_D) then
          
          do ii=1,N_obsl
             
             if (Observations_l(ii)%species==species_L1C_Tbh_D) then  
                
                call ease_convert('EASEv2_M36',                                      &
                     Observations_l(ii)%lat, Observations_l(ii)%lon, col, row)
                
                ! note conversion to one-based indices
                
                if (mask_h_D(nint(col)+1,nint(row)+1))  Observations_l(ii)%assim = .false.

             end if

          end do

       end if
       
       ! apply V-pol masks
       
       if (turnoff_L1C_Tbv_A) then
          
          do ii=1,N_obsl
             
             if (Observations_l(ii)%species==species_L1C_Tbv_A) then  
                
                call ease_convert('EASEv2_M36',                                      &
                     Observations_l(ii)%lat, Observations_l(ii)%lon, col, row)
                
                ! note conversion to one-based indices
                
                if (mask_v_A(nint(col)+1,nint(row)+1))  Observations_l(ii)%assim = .false.

             end if
       
          end do
          
       end if
       
       if (turnoff_L1C_Tbv_D) then
          
          do ii=1,N_obsl
             
             if (Observations_l(ii)%species==species_L1C_Tbv_D) then  
                
                call ease_convert('EASEv2_M36',                                      &
                     Observations_l(ii)%lat, Observations_l(ii)%lon, col, row)
                
                ! note conversion to one-based indices
                
                if (mask_v_D(nint(col)+1,nint(row)+1))  Observations_l(ii)%assim = .false.

             end if
       
          end do
          
       end if
       
       ! clean up
       
       deallocate( mask_h_A )
       deallocate( mask_h_D )
       deallocate( mask_v_A )
       deallocate( mask_v_D )
       
    end if  ! (turnoff_L1C_Tbh_A .or. turnoff_L1C_Tbv_A .or. turnoff_L1C_Tbh_D .or. turnoff_L1C_Tbv_D)
    
  end subroutine turn_off_assim_SMAP_L1CTb
  
  ! *****************************************************************
  
  subroutine read_obs_SMAP_fnames( date_time, this_obs_param,         &
       fname_of_fname_list, N_max, N_fnames, fname_list )
    
    ! read the file within a SMAP Yyyyy/Mmm/Ddd directory that lists
    !  the SMAP h5 file names; preface file names with "Yyyyy/Mmm/Ddd"
    !
    ! reichle,  3 Jan 2014
    ! reichle,  8 Jun 2017: Use "%flistpath" and "%flistname" from "obs_param_type". 
    !
    ! ---------------------------------------------------------------------------------
    
    implicit none
    
    type(date_time_type),             intent(in)  :: date_time
    
    type(obs_param_type),             intent(in)  :: this_obs_param
    
    character( *),                   intent(in)  :: fname_of_fname_list

    integer,                          intent(in)  :: N_max

    integer,                          intent(out) :: N_fnames

    character(100), dimension(N_max), intent(out) :: fname_list
    
    ! local variables
    
    character(300)       :: fname
    character(200)       :: fpath_tmp
    character( 80)       :: fname_tmp
    character( 80)       :: tmpstr80

    character( 14)       :: YYYYMMDDdir
    character(  4)       :: YYYY
    character(  2)       :: MM, DD

    integer              :: ii, istat

    character(len=*), parameter :: Iam = 'read_obs_SMAP_fnames'
    character(len=400) :: err_msg

    ! ---------------------------------------------------------------------
    
    write (YYYY,'(i4.4)') date_time%year
    write (MM  ,'(i2.2)') date_time%month
    write (DD  ,'(i2.2)') date_time%day
    
    YYYYMMDDdir = 'Y' // YYYY // '/M' // MM // '/D' // DD // '/'

    ! initialize default values

    fpath_tmp = this_obs_param%path
    fname_tmp = fname_of_fname_list

    ! use inputs from ensupd nml file if not empty 

    if (len_trim(this_obs_param%flistpath)>0)  fpath_tmp = this_obs_param%flistpath
    if (len_trim(this_obs_param%flistname)>0)  fname_tmp = this_obs_param%flistname
    
    fname = trim(fpath_tmp) // '/' // YYYYMMDDdir // trim(fname_tmp)
    
    open( 10, file=fname, form='formatted', action='read', iostat=istat)
    
    if (istat/=0) then
       err_msg = 'cannot open file ' // trim(fname)
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
        
    if (logit) write(logunit,'(400A)') 'reading file ' // trim(fname)
    
    ! read list of file names

    istat = 0
    ii    = 0
    
    do while (istat==0)
       
       read(10,*,iostat=istat) tmpstr80
       
       if (istat==0) then
          
          ii = ii+1

          if (ii>N_max) then
             err_msg = 'too many files in list'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if

          ! preface file names with "Yyyyy/Mmm/Ddd"

          fname_list(ii) = YYYYMMDDdir // trim(tmpstr80)
          
       else
          
          exit
             
       end if
       
    end do
    
    close(10, status='keep')

    N_fnames = ii

  end subroutine read_obs_SMAP_fnames

  ! *****************************************************************

  type(date_time_type) function SMAPdatetime_to_DateTimeType( SMAP_date_time )
    
    ! convert SMAP date/time strings to LDASsa date_time_type

    ! reichle,  6 Jan 2014
    
    implicit none
    
    ! input arguments
    
    character(len=*) :: SMAP_date_time
    
    ! local variables
    
    type(date_time_type) :: date_time
        
    character(4) :: YYYY
    character(2) :: MM, DD, HH, MI, SS

    integer      :: YYYY_is, YYYY_ie
    integer      :: MM_is,   MM_ie
    integer      :: DD_is,   DD_ie
    integer      :: HH_is,   HH_ie
    integer      :: MI_is,   MI_ie
    integer      :: SS_is,   SS_ie
    
    character(len=*), parameter :: Iam = 'SMAPdatetime_to_DateTimeType'
    character(len=400) :: err_msg

    ! ---------------------------------------------------
    
    if ( (len_trim(SMAP_date_time)==15) .and.            &
         (SMAP_date_time( 9: 9)=='T')          ) then
       
       ! format: "20010727T160914"  (yyyymmddThhmmss)

       YYYY_is =  1;       YYYY_ie =  4;
       MM_is   =  5;       MM_ie   =  6;
       DD_is   =  7;       DD_ie   =  8;
       HH_is   = 10;       HH_ie   = 11;
       MI_is   = 12;       MI_ie   = 13;
       SS_is   = 14;       SS_ie   = 15;
       
    else if (                                            &
         (len_trim(SMAP_date_time)==24) .and.            &
         (SMAP_date_time( 5: 5)=='-')   .and.            &
         (SMAP_date_time( 8: 8)=='-')   .and.            &
         (SMAP_date_time(11:11)=='T')   .and.            &
         (SMAP_date_time(14:14)==':')   .and.            &
         (SMAP_date_time(17:17)==':')   .and.            &
         (SMAP_date_time(20:20)=='.')   .and.            &
         (SMAP_date_time(24:24)=='Z')          ) then
       
       ! format: "2001-07-27T16:09:14.567Z"
       
       YYYY_is =  1;       YYYY_ie =  4;
       MM_is   =  6;       MM_ie   =  7;
       DD_is   =  9;       DD_ie   = 10;
       HH_is   = 12;       HH_ie   = 13;
       MI_is   = 15;       MI_ie   = 16;
       SS_is   = 18;       SS_ie   = 19;
       
    else       

       err_msg = 'invalid SMAPdatetime string'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
    end if
    
    YYYY = SMAP_date_time(YYYY_is:YYYY_ie)
    MM   = SMAP_date_time(MM_is  :MM_ie  )
    DD   = SMAP_date_time(DD_is  :DD_ie  )
    HH   = SMAP_date_time(HH_is  :HH_ie  )
    MI   = SMAP_date_time(MI_is  :MI_ie  )
    SS   = SMAP_date_time(SS_is  :SS_ie  )
    
    read (YYYY,*) date_time%year
    read (MM  ,*) date_time%month
    read (DD  ,*) date_time%day
    read (HH  ,*) date_time%hour
    read (MI  ,*) date_time%min
    read (SS  ,*) date_time%sec
    
    call get_dofyr_pentad( date_time )

    SMAPdatetime_to_DateTimeType = date_time
    
  end function SMAPdatetime_to_DateTimeType
    
  ! *****************************************************************
  
  subroutine read_obs(                                               &
       work_path, exp_id,                                            &
       date_time, dtstep_assim, N_catd, tile_coord,                  &
       tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,          &
       this_obs_param, write_obslog,                                 &
       found_obs, scaled_obs,                                        &
       tmp_obs, tmp_std_obs, tmp_lon, tmp_lat, tmp_time, tmp_assim )
    
    ! read observations and optionally scale observations to model clim
    !
    ! intended to be called by root_proc

    ! 10 Jun 2011 - removed model-based QC for MPI re-structuring (now done
    !               in connection with get_obs_pred())
    ! 22 Nov 2011 - minor clean-up, renamed scale_obs_*() subroutines
    ! 23 Aug 2013 - added possibility of using lat/lon of obs from reader  
    !                (rather than using lat/lon from model tile_coord)
    ! 31 Jan 2014 - added output of time stamp ("tmp_time")
    ! 14 Jul 2014 - added summary diagnostic "scaled_obs" 
    !                (indicates whether any of the founds obs were scaled)
    !  6 Jun 2016 - added functionality to keep obs that cannot be scaled
    !                (but will not be assimilated; SMAP/SMOS Tb only for now)

    implicit none

    ! inputs:
    
    character(*), intent(in) :: work_path
    character(*),  intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: dtstep_assim, N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid_d
    
    integer, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat), intent(in) :: &
         N_tile_in_cell_ij
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij   ! input
    
    type(obs_param_type), intent(in) :: this_obs_param

    logical,              intent(in) :: write_obslog

    ! outputs:
    
    real,    intent(out), dimension(N_catd) :: tmp_obs
    real,    intent(out), dimension(N_catd) :: tmp_std_obs
    real,    intent(out), dimension(N_catd) :: tmp_lon
    real,    intent(out), dimension(N_catd) :: tmp_lat
    real*8,  intent(out), dimension(N_catd) :: tmp_time
    logical, intent(out), dimension(N_catd) :: tmp_assim

    logical, intent(out)                    :: found_obs, scaled_obs

    ! obs time stamp in LDASsa *must* be in J2000 seconds with 'TT12' epoch
    ! (see enkf_types.F90)

    character(4),     parameter :: J2000_epoch_id = 'TT12'  ! see date_time_util.F90

    character(len=*), parameter :: Iam = 'read_obs'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------

    scaled_obs = .false.  ! initialize

    tmp_assim  = .true.   ! initialize

    ! initialize lat/lon/time info (may later be overwritten by individual reader)
    
    tmp_lon  = tile_coord%com_lon
    tmp_lat  = tile_coord%com_lat

    ! obs time stamp in LDASsa *must* be in J2000 seconds with 'TT12' epoch
    ! (see enkf_types.F90)
    ! if needed, must be converted by obs reader (e.g., ASCAT/EUMETSAT)

    tmp_time = datetime_to_J2000seconds(date_time, J2000_epoch_id )
    
    ! -----------------------------
    
    ! choose appropriate reader
    
    select case (trim(this_obs_param%descr))
       
    case ('ae_l2_sm_a', 'ae_l2_sm_d')
       
       call read_obs_ae_l2_sm(                                        &
            work_path, exp_id,                                        &
            date_time, dtstep_assim, N_catd, tile_coord,              &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
            this_obs_param,                                           &
            found_obs, tmp_obs, tmp_std_obs )
       
       ! scale observations to model climatology
       
       if (this_obs_param%scale .and. found_obs) then

          scaled_obs = .true.
          
          call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
               tmp_obs, tmp_std_obs )
          
       end if

       
    case ('ae_sm_LPRM_a_C', 'ae_sm_LPRM_d_C', 'ae_sm_LPRM_a_X', 'ae_sm_LPRM_d_X' )
       
       call read_obs_ae_sm_LPRM(                                      &
            work_path, exp_id,                                        &
            date_time, dtstep_assim, N_catd, tile_coord,              &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
            this_obs_param,                                           &
            found_obs, tmp_obs, tmp_std_obs )
       
       ! scale observations to model climatology
       
       if (this_obs_param%scale .and. found_obs) then
          
          scaled_obs = .true.
          
          call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
               tmp_obs, tmp_std_obs )

       end if

       
    case ('ASCAT_SM_A', 'ASCAT_SM_D' )
      
        call read_obs_sm_ASCAT(                                        &
             work_path, exp_id,                                        &
             date_time, dtstep_assim, N_catd, tile_coord,              &
             tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
             this_obs_param,                                           &
             found_obs, tmp_obs, tmp_std_obs )
        
        ! scale observations to model climatology
        
        if (this_obs_param%scale .and. found_obs) then

           scaled_obs = .true.

           call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
                tmp_obs, tmp_std_obs )
           
        end if
        
    case ('ASCAT_META_SM_A', 'ASCAT_META_SM_D','ASCAT_METB_SM_A', 'ASCAT_METB_SM_D' )
      
        call read_obs_sm_ASCAT_EUMET(                                        &
             work_path, exp_id,                                        &
             date_time, dtstep_assim, N_catd, tile_coord,              &
             tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
             this_obs_param,                                           &
             found_obs, tmp_obs, tmp_std_obs )
        ! Temporay hack to output if ascending/decending instead of tmp_std_obs
        select case (trim(this_obs_param%descr))

        case ('ASCAT_META_SM_A', 'ASCAT_METB_SM_A')
          
            tmp_std_obs(1:N_catd) = 1

        case ('ASCAT_META_SM_D', 'ASCAT_METB_SM_D')

           tmp_std_obs(1:N_catd) = 2

        end select
 
        ! scale observations to model climatology
        
        if (this_obs_param%scale .and. found_obs) then

           scaled_obs = .true.

           call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
                tmp_obs, tmp_std_obs )
           
        end if        

    case ('isccp_tskin_gswp2_v1')
       
       call read_obs_isccp_tskin_gswp2_v1(                            &
            date_time, N_catd, tile_coord,                            &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
            this_obs_param,                                           &
            found_obs, tmp_obs, tmp_std_obs )
       
       if (this_obs_param%scale .and. found_obs) then
          
          scaled_obs = .true.
          
          call scale_obs_tskin_zscore( N_catd, tile_coord,          &
               date_time, this_obs_param,  tmp_obs, tmp_std_obs )
          
       end if
       
    case ('isccp_tskin_ceop3n4')
              
       call read_obs_isccp_tskin_ceop3n4(                             &
            date_time, N_catd, tile_coord,                            &
            this_obs_param,                                           &
            found_obs, tmp_obs, tmp_std_obs )
       
       if (this_obs_param%scale .and. found_obs) then
          
          scaled_obs = .true.

          call scale_obs_tskin_zscore( N_catd, tile_coord,          &
               date_time, this_obs_param,  tmp_obs, tmp_std_obs )
       
       end if
       
    case ('isccp_tskin_ceop3n4_hdASC')
       
       call read_obs_isccp_ts_ceop3n4_hdASC(                          &
            date_time, N_catd, tile_coord,                            &
            this_obs_param,                                           &
            found_obs, tmp_obs, tmp_std_obs )
       
       if (this_obs_param%scale .and. found_obs) then
          
          scaled_obs = .true.
          
          call scale_obs_tskin_zscore( N_catd, tile_coord,          &
               date_time, this_obs_param,  tmp_obs, tmp_std_obs )

       end if
       
    case ('RedArkOSSE_sm')
       
       call read_obs_RedArkOSSE_sm(                                   &
            date_time, N_catd, tile_coord, this_obs_param,            &
            found_obs, tmp_obs, tmp_std_obs )
       
       ! scale observations to model climatology (use AMSR-E subroutine)
       
       if (this_obs_param%scale .and. found_obs) then
          
          scaled_obs = .true.
          
          call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
               tmp_obs, tmp_std_obs )

       end if
       
    case ('RedArkOSSE_CLSMsynthSM')
       
       call read_obs_RedArkOSSE_CLSMsynthSM(                          &
            date_time, N_catd, this_obs_param,                        &
            found_obs, tmp_obs, tmp_std_obs )
       
       ! scale observations to model climatology (use AMSR-E subroutine)
       
       if (this_obs_param%scale .and. found_obs) then
          
          scaled_obs = .true.
          
          call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
               tmp_obs, tmp_std_obs )

       end if
       
    case ('RedArkOSSE_truth_50mm','RedArkOSSE_truth_400mm')
       
       call read_obs_RedArkOSSE_truth(                                   &
            date_time, N_catd, this_obs_param,                           &
            found_obs, tmp_obs, tmp_std_obs )

       ! scale observations to model climatology (use AMSR-E subroutine)
       
       if (this_obs_param%scale .and. found_obs) then
      
          scaled_obs = .true.
          
          call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
               tmp_obs, tmp_std_obs )

       end if

       ! assimilation NOT implemented for RedArkOSSE_truth obs
       
       if (this_obs_param%assim) then
          err_msg = 'assimilation NOT implemented for RedArkOSSE_truth obs'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       
    case ('VivianaOK_CLSMsynthSM')
       
       call read_obs_VivianaOK_CLSMsynthSM(                           &
            date_time, N_catd, this_obs_param,                        &
            found_obs, tmp_obs, tmp_std_obs )
       
       ! scale observations to model climatology (use AMSR-E subroutine)
       
       if (this_obs_param%scale .and. found_obs) then
          
          scaled_obs = .true.
          
          call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
               tmp_obs, tmp_std_obs )
          
       end if

    case('SMOS_SM_A','SMOS_SM_D')
       
       call read_obs_SMOS(                                            &
            date_time, N_catd, this_obs_param,                        &
            dtstep_assim, tile_coord, tile_grid_d,                    &
            N_tile_in_cell_ij, tile_num_in_cell_ij, write_obslog,     &
            found_obs, tmp_obs, tmp_std_obs, tmp_lon, tmp_lat )
       
       ! scale observations to model climatology
       
       if (this_obs_param%scale .and. found_obs) then

          scaled_obs = .true.

          call scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,   &
               tmp_obs, tmp_std_obs )
          
       end if

    case('SMOS_reg_Tbh_A','SMOS_reg_Tbh_D','SMOS_reg_Tbv_A','SMOS_reg_Tbv_D', &
         'SMOS_fit_Tbh_A','SMOS_fit_Tbh_D','SMOS_fit_Tbv_A','SMOS_fit_Tbv_D')
       
       call read_obs_SMOS(                                            &
            date_time, N_catd, this_obs_param,                        &
            dtstep_assim, tile_coord, tile_grid_d,                    &
            N_tile_in_cell_ij, tile_num_in_cell_ij, write_obslog,     &
            found_obs, tmp_obs, tmp_std_obs, tmp_lon, tmp_lat )
       
       ! scale observations to model climatology XXXXXXX
       
       if (this_obs_param%scale .and. found_obs) then

          scaled_obs = .true.

          call scale_obs_Tb_zscore( N_catd, tile_coord, date_time,  &
               this_obs_param, tmp_obs, tmp_std_obs, tmp_assim )
          
       end if

    case('SMAP_L1C_Tbh_A',     'SMAP_L1C_Tbv_A',        &
         'SMAP_L1C_Tbh_D',     'SMAP_L1C_Tbv_D',        &
         'SMAP_L1C_Tbh_E09_A', 'SMAP_L1C_Tbv_E09_A',    &
         'SMAP_L1C_Tbh_E09_D', 'SMAP_L1C_Tbv_E09_D',    &
         'SMAP_L1C_Tbh_E27_A', 'SMAP_L1C_Tbv_E27_A',    &
         'SMAP_L1C_Tbh_E27_D', 'SMAP_L1C_Tbv_E27_D',    &
         'SMAP_L2AP_Tbh_A',    'SMAP_L2AP_Tbv_A',       &
         'SMAP_L2AP_Tbh_D',    'SMAP_L2AP_Tbv_D' )
       
       call read_obs_SMAP_halforbit_Tb(                                     &
            date_time, N_catd, this_obs_param,                              &
            dtstep_assim, tile_coord, tile_grid_d,                          &
            N_tile_in_cell_ij, tile_num_in_cell_ij, write_obslog,           &
            found_obs, tmp_obs, tmp_std_obs, tmp_lon, tmp_lat, tmp_time )
       
       ! scale observations to model climatology XXXXXXX
       
       if (this_obs_param%scale .and. found_obs) then

          scaled_obs = .true.
          
          call scale_obs_Tb_zscore( N_catd, tile_coord, date_time,  &
               this_obs_param, tmp_obs, tmp_std_obs, tmp_assim )

       end if
       
    case('LaRC_tskin-GOESW', 'LaRC_tskin-GOESE', 'LaRC_tskin-MET09',  & 
         'LaRC_tskin-FY2E-', 'LaRC_tskin-MTST2')
              
       call read_obs_LaRC_Tskin(                                      &
            work_path, exp_id,                                        &
            date_time, dtstep_assim, N_catd, tile_coord,              &
            tile_grid_d, N_tile_in_cell_ij, tile_num_in_cell_ij,      &
            this_obs_param,                                           &
            found_obs, tmp_obs, tmp_std_obs )

       ! NOT IMPLEMENTED: scale observations to model climatology

    case default
       
       err_msg = 'unknown obs_param%descr'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
    
  end subroutine read_obs
  
    
  ! *****************************************************************
  
  subroutine scale_obs_sfmc_cdf( N_catd, tile_coord, this_obs_param,  &
       tmp_obs, tmp_std_obs )
    
    ! scale sfmc obs to model climatology via cdf matching
    ! 
    ! use matlab functions "get_cdf_match_AMSR.m" and "get_model_and_obs_stats.m" 
    ! to create input scaling files
    !
    ! IMPORTANT: Make sure that model and obs data are in the SAME UNITS prior
    !            to generating the input scaling files with the matlab routines.  
    !            Otherwise, the (linear) rescaling outside of the observed 
    !            range (between obs_min and obs_max) will fail!  For example,
    !            ASCAT soil moisture retrievals must first be converted into 
    !            volumetric units that range roughly from 0 to 0.5 m3/m3 
    !            (for some porosity) before they can be rescaled to the 
    !            Catchment model's sfmc (which is in m3/m3).
    
    ! THIS SUBROUTINE COULD USE WORK...
    ! - map from scaling parameter domain to model domain
    ! reichle, 14 Oct 2005
    ! reichle, 26 Sep 2006
    !
    ! bug fix: scaling via linear interpolation outside of [obs_min,obs_max]
    ! reichle, 27 Jan 2006    
    !
    ! reichle, 22 Nov 2011 - renamed subroutine, minor clean-up, added comments
    ! reichle,  9 Nov 2012 - edited to enable use of stats file that do not 
    !                         perfectly match current domain (as in tile_coord)
    
    implicit none
    
    integer, intent(in) :: N_catd
        
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
        
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! inout
    
    real,    intent(inout), dimension(N_catd) :: tmp_obs
    real,    intent(inout), dimension(N_catd) :: tmp_std_obs
    
    ! local variables
    
    real, parameter :: no_data_stats = -9999.
    
    real, parameter :: tol = 1e-2
    
    character(300) :: fname
    
    integer :: i, j, N_sclprm, N_poly, ind
    
    real :: edge_min, edge_max, edge_dx
    
    real :: tmpreal, x, x0, x1, y0, y1
    
    integer, dimension(:),   allocatable :: tmp_tile_id
    
    real,    dimension(:),   allocatable :: std_obs, std_mod, min_obs, max_obs
    
    real,    dimension(:,:), allocatable :: fit_coeff

    character(len=*), parameter :: Iam = 'scale_obs_sfmc_cdf'
    character(len=400) :: err_msg
    
    ! ------------------------------------------------------------------
    
    ! read scaling parameters from file
    
    fname = trim(this_obs_param%scalepath) // '/' // &
         trim(this_obs_param%scalename)
    
    if (logit) write (logunit,*)        'scaling obs species ', this_obs_param%species, ':'
    if (logit) write (logunit,'(400A)') '  reading ', trim(fname)
    
    open(10, file=fname, form='formatted', action='read')
    
    read(10,*) N_sclprm, N_poly, edge_min, edge_max, edge_dx 
    
    ! minimal consistency check
    
    if (N_catd>N_sclprm) then
       err_msg = 'N_sclprm too small'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    allocate( tmp_tile_id( N_sclprm           ))
    allocate( std_obs(     N_sclprm           ))
    allocate( std_mod(     N_sclprm           ))
    allocate( min_obs(     N_sclprm           ))
    allocate( max_obs(     N_sclprm           ))
    allocate( fit_coeff(   N_sclprm, N_poly+1 ))
    
    do i=1,N_sclprm
       
       read (10,*) tmp_tile_id(i), std_obs(i), std_mod(i),  &
            min_obs(i), max_obs(i), fit_coeff(i,:)
       
    end do
    
    close(10,status='keep')
    
    ! --------------------------------------------------------------
    
    ! scale observations
    
    do i=1,N_catd
       
       ! check for no-data-values in observation
       !  (any negative number could be no-data-value for observations)

       if (tmp_obs(i)>=0.)  then                    
          
          ! find ind for current tile id in scaling parameters
          
          if (tmp_tile_id(i)==tile_coord(i)%tile_id) then
             
             ind = i  ! educated guess for global domain
             
          else
             
             do ind=1,N_sclprm
                
                if (tmp_tile_id(ind)==tile_coord(i)%tile_id)  exit
                
             end do
             
             if (ind>N_sclprm) then
                err_msg = 'tile_id not found in scaling parameter file'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if
             
          end if
          
          
          ! make sure obs falls within range valid for scaling and
          ! check for no-data-value in fit parameters
          
          if ( (tmp_obs(i)>=edge_min)                       .and.          &
               (tmp_obs(i)<=edge_max)                       .and.          &
               (abs(fit_coeff(ind,1)-no_data_stats) > tol)  .and.          &
               (abs(fit_coeff(ind,2)-no_data_stats) > tol)        ) then 
             
             ! evaluate polynomial fit
             !
             ! fit_coeff(ind,1) is coeff for highest order
             ! fit_coeff(ind,N_poly+1) is constant term
             
             ! evaluate polynomial at tmp_obs if min_obs<=tmp_obs<=max_obs,
             ! otherwise at min_obs or max_obs (accordingly)
             
             x = min( max( tmp_obs(i), min_obs(ind) ), max_obs(ind) )
             
             tmpreal = fit_coeff(ind,1) 
             
             do j=1,N_poly
                
                tmpreal = tmpreal*x + fit_coeff(ind,j+1)
                
             end do
             
             if (tmp_obs(i)<min_obs(ind)) then
                
                ! linear interpolation between min(edges) and min_obs
                ! (NOTE: model and obs data must be in SAME units, 
                !  reichle, 22 Nov 2011)
                
                y1 = tmpreal
                y0 = edge_min
                x1 = min_obs(ind)
                x0 = edge_min
                
                tmp_obs(i) = (y1-y0)/(x1-x0)*( tmp_obs(i) - x0 ) + y0
                
             elseif (tmp_obs(i)>max_obs(ind)) then
                
                ! linear interpolation between max_obs and max(edges)
                ! (NOTE: model and obs data must be in SAME units, 
                !  reichle, 22 Nov 2011)
                
                y1 = edge_max
                y0 = tmpreal
                x1 = edge_max
                x0 = max_obs(ind)
                
                tmp_obs(i) = (y1-y0)/(x1-x0)*( tmp_obs(i) - x0 ) + y0
                
             else
                
                ! accept polynomial fit as is
                
                tmp_obs(i) = tmpreal
             
             end if
             
             ! scale observation error std
          
             tmp_std_obs(i) = std_mod(ind)/std_obs(ind)*tmp_std_obs(i)
             
          else
          
             tmp_obs(i) = this_obs_param%nodata
             
          end if
          
       end if

       ! qc check after scaling
       
       if ((tmp_obs(i)>edge_max) .or. (tmp_obs(i)<edge_min)) &
            tmp_obs(i) = this_obs_param%nodata
       
       
    end do
    
    deallocate( tmp_tile_id )
    deallocate( std_obs     )
    deallocate( std_mod     )
    deallocate( min_obs     )
    deallocate( max_obs     )
    deallocate( fit_coeff   )    

  end subroutine scale_obs_sfmc_cdf
  
  ! *****************************************************************
  
  subroutine scale_obs_tskin_zscore( N_catd, tile_coord,   &
       date_time, this_obs_param,                          &
       tmp_obs, tmp_std_obs )

    ! scale tskin obs to model climatology via standard-normal-deviate (zscore)
    ! scaling
    ! 
    ! use matlab functions "get_cdf_match_AMSR.m" and "get_model_and_obs_stats.m" 
    ! to create input scaling files
    !
    ! IMPORTANT: Make sure that model and obs data are in the SAME UNITS prior
    !            to generating the input scaling files with the matlab routines.  
    
    ! reichle, 14 Oct 2005
    !
    ! reichle, 22 Nov 2011 - renamed subroutine, minor clean-up, added comments
    
    implicit none
    
    integer, intent(in) :: N_catd
        
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(date_time_type), intent(in) :: date_time
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! inout
    
    real,    intent(inout), dimension(N_catd) :: tmp_obs
    real,    intent(inout), dimension(N_catd) :: tmp_std_obs
    
    ! ----------------------------------------------------------
    
    ! local variables
    
    real, parameter :: no_data_stats = -9999.
    
    real, parameter :: tol = 1e-2
    
    character(3), dimension(12) :: month_string = (/ &
         'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',   &
         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /)
    
    ! -------------------
    
    character(300) :: fname
    
    character(2) :: tmpchar2
    
    integer :: i, ind, N_sclprm
    
    real :: tmpreal
    
    integer, dimension(:), allocatable :: sclprm_tile_id
    
    real, dimension(:), allocatable :: sclprm_lon,      sclprm_lat 
    real, dimension(:), allocatable :: sclprm_mean_obs, sclprm_std_obs
    real, dimension(:), allocatable :: sclprm_mean_mod, sclprm_std_mod
    
    character(len=*), parameter :: Iam = ' scale_obs_tskin_zscore'
    character(len=400) :: err_msg
    
    ! ------------------------------------------------------------------
    
    write (tmpchar2, '(i2.2)') date_time%hour
    
    ! read scaling parameters from file
    
    fname = trim(this_obs_param%scalepath) // '/' // &
         trim(this_obs_param%scalename)    //        &
         month_string(date_time%month)     // '_' // &
         tmpchar2 // 'z.bin'
    
    if (logit) write (logunit,*)        'scaling obs species ', this_obs_param%species, ':'
    if (logit) write (logunit,'(400A)') '  reading ', trim(fname)
    
    open(10, file=fname, form='unformatted',convert='big_endian', action='read')
    
    read(10) N_sclprm
    
    allocate(sclprm_tile_id(N_sclprm))    
    allocate(sclprm_lon(N_sclprm))     
    allocate(sclprm_lat(N_sclprm))          
    allocate(sclprm_mean_obs(N_sclprm))     
    allocate(sclprm_std_obs(N_sclprm))      
    allocate(sclprm_mean_mod(N_sclprm))     
    allocate(sclprm_std_mod(N_sclprm))      
    
    read(10) sclprm_tile_id
    read(10) sclprm_lon
    read(10) sclprm_lat 
    read(10) sclprm_mean_obs
    read(10) sclprm_std_obs
    read(10) sclprm_mean_mod
    read(10) sclprm_std_mod 
    
    close(10,status='keep')
    
    ! minimal consistency check
    
    if (N_catd>N_sclprm) then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'something is wrong')
    end if
    
    ! --------------------------------------------------------------
    
    ! scale observations (at this point all obs are of type 
    ! isccp_tskin_gswp2_v1 because of the way the subroutine is called
    ! from subroutine read_obs())
    
    do i=1,N_catd
       
       ! check for no-data-values in observation (any neg Tskin is no_data)
       
       if (tmp_obs(i)>=0.) then
          
          ! find ind for current tile id in scaling parameters
          
          do ind=1,N_sclprm
             
             if (sclprm_tile_id(ind)==tile_coord(i)%tile_id)  exit
             
          end do
          
          ! sanity check (against accidental use of wrong tile space)
          
          if ( abs(tile_coord(i)%com_lat-sclprm_lat(ind))>tol  .or.             &
               abs(tile_coord(i)%com_lon-sclprm_lon(ind))>tol        ) then
             err_msg = 'something wrong'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          ! check for no-data-values in observation and fit parameters
          ! (any negative number could be no-data-value for observations)
          
          if ( sclprm_mean_obs(ind)>0.                       .and.          &
               sclprm_mean_mod(ind)>0.                       .and.          &
               sclprm_std_obs(ind)>=0.                       .and.          &
               sclprm_std_mod(ind)>=0.                             ) then
             
             ! scale via standard normal deviates
             
             tmpreal = sclprm_std_mod(ind)/sclprm_std_obs(ind) 
             
             tmp_obs(i) = sclprm_mean_mod(ind)                       &
                  + tmpreal*(tmp_obs(i)-sclprm_mean_obs(ind)) 
                          
             ! scale observation error std
             
             tmp_std_obs(i) = tmpreal*tmp_std_obs(i)
             
          else
             
             tmp_obs(i) = this_obs_param%nodata
             
          end if
          
       end if
       
    end do
    
    deallocate(sclprm_tile_id)
    deallocate(sclprm_lon)     
    deallocate(sclprm_lat)          
    deallocate(sclprm_mean_obs)     
    deallocate(sclprm_std_obs)      
    deallocate(sclprm_mean_mod)     
    deallocate(sclprm_std_mod)      
    
  end subroutine scale_obs_tskin_zscore

  ! ********************************************************************************
  
  subroutine scale_obs_Tb_zscore( N_catd, tile_coord, date_time, this_obs_param,   &
       tmp_obs, tmp_std_obs, tmp_assim )
    
    ! Scale Tb obs to model climatology via standard-normal-deviate (zscore)
    ! scaling or mean-only scaling.
    !
    ! The type of scaling is determined by the first 5 characters of the 
    ! nml input "this_obs_param%scalename" as follows:
    !
    !  this_obs_param%scalename = 'ScZS_' : Scale using ZScore (mean and std-dev)
    !  this_obs_param%scalename = 'ScMO_' : Scale using Mean Only
    !
    ! These first 5 characters are NOT part of the file name for the scaling file.
    !
    ! Use matlab functions "get_model_and_obs_clim_stats_oct11.m" 
    ! and "scaling_prep_multi_year.m" by Gabrielle De Lannoy to create 
    ! global input scaling files for each pentad.
    !    
    ! Note that the tiles in the input scaling parameter files do NOT need
    ! to be in the same order  as in the domain tile_coord vector (although
    ! the fastest execution is obtained when the scaling parameter files
    ! and the domain tile_coord vector are for exactly the same tiles in 
    ! the same order.
    !
    ! The vector "tmp_obs" is of length N_catd.  The observations in "tmp_obs"
    ! are mapped into tile space ONLY for book-keeping (each observation has
    ! been assigned to a model tile).  The observations are not necessarily
    ! representative of the assigned tile and could be for a larger area.
    ! The scaling parameter file therefore is specific to a model tile space.
    ! Global scaling parameter files *can* be used with subdomains. 
    !
    ! Valid statistics are not necessarily available for each observed 
    ! location and day-of-year.
    !
    ! Gabrielle De Lannoy (GDL), 10 Sep 2012 - first draft 
    ! reichle,                   12 Oct 2012 - revised and merged into CVS trunk
    ! reichle,                    6 Jun 2016 - keep obs that can *not* be scaled, 
    !                                           but set tmp_assim=.false.

    implicit none
    
    integer, intent(in) :: N_catd
        
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(date_time_type), intent(in) :: date_time
    
    type(obs_param_type), intent(in) :: this_obs_param
    
    ! inout
    
    real,    intent(inout), dimension(N_catd) :: tmp_obs
    real,    intent(inout), dimension(N_catd) :: tmp_std_obs
    logical, intent(inout), dimension(N_catd) :: tmp_assim
    
    ! ----------------------------------------------------------
    
    ! local variables
    
    real,    parameter :: no_data_stats = -9999.
    
    real,    parameter :: tol = 1e-2
            
    ! -------------------
        
    logical        :: hpol, scale_mean_only
    
    character(300) :: fname
    
    character( 80) :: tmpstring80
    character(  2) :: tmpstring2, orbit_flag
    
    integer        :: i, ind, istat, ind_angle
    
    integer        :: asc_flag, N_data_min, N_sclprm, N_ang
    
    real           :: tmpreal
    
    integer, dimension(:), allocatable :: sclprm_tile_id
    
    real,    dimension(:), allocatable :: sclprm_ang
    real,    dimension(:), allocatable :: sclprm_lon,      sclprm_lat 
    real,    dimension(:), allocatable :: sclprm_mean_obs, sclprm_std_obs
    real,    dimension(:), allocatable :: sclprm_mean_mod, sclprm_std_mod
    
    character(len=*), parameter :: Iam = 'scale_obs_Tb_zscore'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------------
    
    ! determine whether to scale mean and std-dev or mean only

    tmpstring80 = this_obs_param%scalename
    
    if     (tmpstring80(1:5)=='ScZS_') then
       
       scale_mean_only = .false.
       
    elseif (tmpstring80(1:5)=='ScMO_') then
       
       scale_mean_only = .true.
       
    else
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown scaling method')
       
    end if

    ! assemble the name of the file with scaling parameters
    !
    ! - first 5 chars of this_obs_param%scalename are NOT part of the file name
    ! - different scaling files for each orbit and pentad
    
    if     (this_obs_param%orbit==1) then
       
       orbit_flag = '_A'
       
    elseif (this_obs_param%orbit==2) then
       
       orbit_flag = '_D'

    else
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown orbit')
       
    end if
    
    write (tmpstring2, '(i2.2)') date_time%pentad
    
    fname =                                                 &
         trim(this_obs_param%scalepath)       // '/'  //    &
         trim(this_obs_param%scalename(6:80)) //            &
         orbit_flag // '_p' // tmpstring2 // '.bin'
    
    if (logit) write (logunit,*)        'scaling obs species ', this_obs_param%species, ':'
    if (logit) write (logunit,'(400A)') '  reading ', trim(fname)
    
    open(10, file=fname, form='unformatted',convert='big_endian',status='old', &
         access='SEQUENTIAL', iostat=istat)
    
    if (istat/=0) then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'could not open file')
    end if
    

    ! read file header 
    !
    ! file format of scaling files mirrors that of pre-processed SMOS obs files

    read(10) asc_flag, N_data_min
    read(10)                 ! start time of interval for stats computation (not used)
    read(10)                 ! end   time of interval for stats computation (not used)
    read(10) N_sclprm, N_ang

    ! minimal consistency checks

    if ( (this_obs_param%orbit==1 .and. asc_flag/=1) .or.       &
         (this_obs_param%orbit==2 .and. asc_flag/=0)      ) then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'orbit mismatch')
    end if
    
    if (logit) write (logunit,*) '  asc_flag   = ', asc_flag
    if (logit) write (logunit,*) '  N_data_min = ', N_data_min
    if (logit) write (logunit,*) '  N_sclprm   = ', N_sclprm
    if (logit) write (logunit,*) '  N_ang      = ', N_ang
    
    allocate(sclprm_ang(     N_ang   ))
    
    allocate(sclprm_lon(     N_sclprm))
    allocate(sclprm_lat(     N_sclprm))
    allocate(sclprm_tile_id( N_sclprm))
    
    allocate(sclprm_mean_obs(N_sclprm))     
    allocate(sclprm_std_obs( N_sclprm))      
    allocate(sclprm_mean_mod(N_sclprm))     
    allocate(sclprm_std_mod( N_sclprm)) 

    ! read angle and location information
    
    read(10) sclprm_ang
    
    read(10) sclprm_lon      !only valid values where obs were available
    read(10) sclprm_lat      !only valid values where obs were available
    read(10) sclprm_tile_id
    
    ! find the index for the angle of interest
    ! NOTE: after processing of namelist inputs, each species 
    !       has a unique angle (see subroutine read_ens_upd_inputs())
    
    ind_angle = -9999
    
    do i=1,N_ang
       
       if (abs(sclprm_ang(i)-this_obs_param%ang(1))<0.01)  ind_angle = i
       
    end do
    
    if (ind_angle<0) then
       err_msg = 'problem with incidence angle'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    ! need h-pol or v-pol?
    
    if     (this_obs_param%pol==1) then
       
       hpol = .true.
       
    elseif (this_obs_param%pol==2) then
       
       hpol = .false.
       
    else
       
       err_msg = 'unknown polarization'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    ! read scaling parameters
    !
    ! for each field, loop over all angles, read stats for angle of interest
    !
    ! blocks (1- 5) after header: Tbh stats
    ! blocks (6-10) after header: Tbv stats
    
    if (.not. hpol) then  ! in case of V-pol, skip through H-pol entries
       
       do i=1,N_ang ! block 1 - mean_obs Tbh           
          read(10) 
       end do
       
       do i=1,N_ang ! block 2 - std_obs Tbh           
          read(10) 
       end do
       
       do i=1,N_ang ! block 3 - mean_mod Tbh           
          read(10) 
       end do
       
       do i=1,N_ang ! block 4 - std_mod Tbh           
          read(10) 
       end do
       
       do i=1,N_ang ! block 5 - N_data Tbh          
          read(10)  
       end do
       
    end if
    
    ! from each block, read stats for angle of interest
    
    do i=1,N_ang           ! block 1 (h-pol) or 6 (v-pol) - mean_obs 
       
       if (i==ind_angle) then
          read(10) sclprm_mean_obs
       else
          read(10)
       end if
       
    end do
    
    do i=1,N_ang           ! block 2 (h-pol) or 7 (v-pol) - std_obs
       
       if (i==ind_angle) then
          read(10) sclprm_std_obs
       else
          read(10)
       end if
       
    end do
    
    do i=1,N_ang           ! block 3 (h-pol) or 8 (v-pol) - mean_mod
       
       if (i==ind_angle) then
          read(10) sclprm_mean_mod
       else
          read(10)
       end if
       
    end do
    
    do i=1,N_ang           ! block 4 (h-pol) or 9 (v-pol) - std_mod
       
       if (i==ind_angle) then
          read(10) sclprm_std_mod
       else
          read(10)
       end if
       
    end do
    
    !do i=1,N_ang           ! block 5 (h-pol) or 10 (v-pol) - N_data
    !   
    !   if (i==ind_angle) then
    !      read(10) sclprm_Ndata
    !   else
    !      read(10)
    !   end if
    !   
    !end do
    
    close(10,status='keep')
    
       
    ! --------------------------------------------------------------
    
    ! scale observations (at this point all obs are of type Tb because 
    ! of the way the subroutine is called from subroutine read_obs())
    
    do i=1,N_catd
       
       ! check for no-data-values in observation (any neg Tb is no_data)
       
       if (tmp_obs(i)>=0.) then
          
          ! find ind for current tile id in scaling parameters

          ind = -9999  ! initialize to negative integer
          
          if (N_sclprm==N_catd) then
             
             ! try an educated guess (to avoid an additional loop)
             
             if (sclprm_tile_id(i)==tile_coord(i)%tile_id)  ind = i  
             
          end if
          
          if (ind<0) then  ! the educated guess failed, try again
             
             do ind=1,N_sclprm
                
                if (sclprm_tile_id(ind)==tile_coord(i)%tile_id)  exit
                
             end do

             if (ind>N_sclprm) then
                err_msg = 'tile_id not found in scaling parameter file'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if

          end if
          
          ! check for no-data-values in observation and fit parameters
          ! (any negative number could be no-data-value for observations)
          
          if ( sclprm_mean_obs(ind)>0.   .and.          &
               sclprm_mean_mod(ind)>0.   .and.          &
               sclprm_std_obs( ind)>0.   .and.          &
               sclprm_std_mod( ind)>0.          ) then
             
             
             ! sanity check (against accidental use of wrong tile space)
             
             if ( abs(tile_coord(i)%com_lat-sclprm_lat(ind))>tol .or.            &
                  abs(tile_coord(i)%com_lon-sclprm_lon(ind))>tol       ) then
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'something wrong')
             end if

             ! scale
             
             if (scale_mean_only) then
                   
                ! adjust mean only
                
                tmp_obs(i) =                                                   &
                     tmp_obs(i) - sclprm_mean_obs(ind) + sclprm_mean_mod(ind) 
                
             else
                
                ! scale via standard normal deviates
                
                tmpreal = sclprm_std_mod(ind)/sclprm_std_obs(ind) 
                
                tmp_obs(i) = sclprm_mean_mod(ind)                              &
                     + tmpreal*(tmp_obs(i)-sclprm_mean_obs(ind)) 
                
                ! scale observation error std
                
                tmp_std_obs(i) = tmpreal*tmp_std_obs(i)
                
             end if
             
          else
             
             ! keep unscaled obs, provide only "do not assimilate" info
             ! reichle, 6 Jun 2016

             tmp_assim(i) = .false.
             
             !tmp_obs(i) = this_obs_param%nodata
             
          end if
          
       end if
       
    end do
    
    deallocate(sclprm_ang)
    
    deallocate(sclprm_tile_id)
    deallocate(sclprm_lon)     
    deallocate(sclprm_lat)          
    
    deallocate(sclprm_mean_obs)     
    deallocate(sclprm_std_obs)      
    deallocate(sclprm_mean_mod)     
    deallocate(sclprm_std_mod) 
    
  end subroutine scale_obs_Tb_zscore
  
  ! *****************************************************************
  
  subroutine collect_obs(                                                             &
       work_path, exp_id, date_time, dtstep_assim,                                    &
       N_catl, tile_coord_l,                                                          &
       N_catf, tile_coord_f, tile_grid_f, N_tile_in_cell_ij_f, tile_num_in_cell_ij_f, &
       N_catl_vec, low_ind, l2f,                                                      &
       N_obs_param, obs_param, N_obsl_max, write_obslog,                              &
       N_obsl, Observations_l, found_obs_f )
    
    ! check for observations that must be assimilated, 
    ! collect into measurement vector
    !
    ! a total of "N_obsl" observations are returned in "Observations_l"
    !
    ! 25 Jul 2005 - rewritten
    ! 10 Jun 2011 - re-structured for MPI, incl. removal of model-based QC 
    !               (now done in connection with get_obs_pred())
    ! 24 Dec 2013 - added "obslog" output files (not yet implemented for all readers!)
    ! 31 Dec 2014 - added time stamp to Observations
    ! 21 Mar 2014 - sort Observations_l by tilenum and then species to avoid lay-out 
    !                dependency for MPI parallel execution

    implicit none
    
    character(*),                               intent(in)  :: work_path
    character(*),                                intent(in)  :: exp_id
    
    type(date_time_type),                         intent(in)  :: date_time
    
    integer,                                      intent(in)  :: dtstep_assim
    integer,                                      intent(in)  :: N_catl, N_catf

    integer, dimension(numprocs),                 intent(in)  :: N_catl_vec, low_ind

    integer, dimension(N_catl),                   intent(in)  :: l2f

    ! tile_coord_f of catchments in domain (length N_catf)
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord_l, tile_coord_f  ! input
    
    type(grid_def_type),                          intent(in)  :: tile_grid_f
    
    ! N_tile_in_cell_ij and tile_num_in_cell_ij are on the "full" domain
    !  and guaranteed to be allocated ONLY for the root_proc
    !  (but may be allocated on all processors depending on obs_param%FOV)

    integer, dimension(:,:),   pointer :: N_tile_in_cell_ij_f   ! input 
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij_f ! input 
    
    integer,                                      intent(in)  :: N_obs_param
    
    type(obs_param_type), dimension(N_obs_param), intent(in)  :: obs_param  
    
    integer,                                      intent(in)  :: N_obsl_max
    
    logical,                                      intent(in)  :: write_obslog
    
    integer,                                      intent(out) :: N_obsl

    type(obs_type), dimension(N_obsl_max),        intent(out) :: Observations_l
    
    logical,                                      intent(out) :: found_obs_f

    ! ----------------------------
    
    ! locals
    
    logical :: found_obs, scaled_obs, any_scaled_obs

    integer :: obs_count, species
    integer :: ii, ind_start, ind_end, N_tmp, this_tilenum, this_tilenum_new
    
    real,   dimension(N_catl) :: tmp_obs,   tmp_std_obs,   tmp_lon,   tmp_lat
    real,   dimension(N_catf) :: tmp_obs_f, tmp_std_obs_f, tmp_lon_f, tmp_lat_f

    real*8, dimension(N_catl) :: tmp_time
    real*8, dimension(N_catf) :: tmp_time_f

    logical, dimension(N_catl) :: tmp_assim
    logical, dimension(N_catf) :: tmp_assim_f

    integer, dimension(numprocs) :: N_obsl_vec
    
    integer, dimension(N_obs_param) :: tmp_species
    
    integer, dimension(:), allocatable :: indx, tilenums
    
    character( 12) :: tmpstr12
    character(300) :: tmpstr300

    character(len=*), parameter :: Iam = 'collect_obs'
    character(len=400) :: err_msg
    
    ! ---------------------------------------------------------------

    if (logit) write (logunit,*) 'collecting observations...' 
    
    ! obs_count serves as counter for total number of observations,
    ! excluding no-data-values
    
    obs_count = 0

    any_scaled_obs = .false.      ! initialize
    
    ! ----------------------------------------------------------
    
    do species = 1,N_obs_param
       
       ! make sure species number here is consistent with
       !  definitions in nml input file
       
       if (obs_param(species)%species .ne. species) then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'something wrong')
       end if
                 
       if (root_proc) then
          
          ! subroutine read_obs() reads all observations in obs files 
          ! (typically global) and returns a vector in (full domain) 
          ! tile space with the values of the observations (at most one
          ! observation per tile and per species)
          
          call read_obs(                                                               &
               work_path, exp_id,                                                      &
               date_time, dtstep_assim, N_catf, tile_coord_f,                          &
               tile_grid_f, N_tile_in_cell_ij_f, tile_num_in_cell_ij_f,                &
               obs_param(species), write_obslog,                                       &
               found_obs, scaled_obs,                                                  &
               tmp_obs_f, tmp_std_obs_f, tmp_lon_f, tmp_lat_f, tmp_time_f, tmp_assim_f)
          
          if (scaled_obs)  any_scaled_obs = .true.

       end if
       
       ! put "tmp_obs" (in "tile" space) into "compressed" vector "Observations_l"
       !
       ! each observation is "managed" ("administered") by the processor that 
       ! contains the tile to which the observation was assigned
       !
       ! the assignment of each observation to tile space does NOT imply that
       ! only the observations from one (local) processor impact the increments
       ! for that processor --> see "halo" below          
       
#ifdef LDAS_MPI
       
       call MPI_BCAST(found_obs, 1, MPI_LOGICAL, 0,mpicomm, mpierr)
       
#endif       
       
       if (found_obs)  then                                            
          
          ! map from full to local domain
          
          call f2l_real(   N_catf,N_catl,N_catl_vec,low_ind, tmp_obs_f,     tmp_obs)
          call f2l_real(   N_catf,N_catl,N_catl_vec,low_ind, tmp_std_obs_f, tmp_std_obs)
          call f2l_real(   N_catf,N_catl,N_catl_vec,low_ind, tmp_lon_f,     tmp_lon)
          call f2l_real(   N_catf,N_catl,N_catl_vec,low_ind, tmp_lat_f,     tmp_lat)

          call f2l_real8(  N_catf,N_catl,N_catl_vec,low_ind, tmp_time_f,    tmp_time)
          
          call f2l_logical(N_catf,N_catl,N_catl_vec,low_ind, tmp_assim_f,   tmp_assim)
          

          ! NOTE: "Observations" here are l(ocal) obs only
          
          call put_into_Observations( obs_param(species), N_obsl_max, N_catl, l2f,  &
               tmp_obs, tmp_std_obs, tmp_lon, tmp_lat, tmp_time, tmp_assim,         &
               obs_count, Observations_l )
          
       end if
       
    end do

#ifdef LDAS_MPI
       
    call MPI_BCAST(any_scaled_obs, 1, MPI_LOGICAL, 0,mpicomm, mpierr)

#endif
    
    N_obsl = obs_count

    ! -----------------------------------------------------------------
    !
    ! sort relevant elements of Observations_l by tilenum and then by species 
    ! to avoid lay-out dependency for MPI parallel execution
    ! - reichle, 21 March 2014
    
    if (N_obsl>1) then                        ! sort only if 2 or more obs
       
       allocate(indx(    N_obsl))
       
       allocate(tilenums(N_obsl))
       
       tilenums = Observations_l(1:N_obsl)%tilenum
       
       ! get index vector,   NOTE:  
       !                       nr_indexx() does not change input "arr"
       !                       nr_indexx() only works with *real* input "arr"
       
       call nr_indexx( N_obsl, real(tilenums), indx(1:N_obsl) )
       
       ! apply sort by tilenum
       
       Observations_l(1:N_obsl) = Observations_l(indx(1:N_obsl))
       
       ! now make sure that within each tilenum, Observations_l are sorted by species
       
       ind_start    = 1
       
       this_tilenum = Observations_l(ind_start)%tilenum
       
       do ii=2,N_obsl
          
          this_tilenum_new = Observations_l(ii)%tilenum
          
          if ( (this_tilenum_new/=this_tilenum) .or. (ii==N_obsl) ) then
             
             if ( (this_tilenum_new/=this_tilenum) .and. (ii<=N_obsl) ) then
                
                ind_end = ii-1
                
             else
                
                ind_end = N_obsl
                
             end if
             
             ! Observations_l(ind_start:ind_end) are the complete subset 
             !  of (local) obs with the same tilenum
             
             N_tmp   = ind_end - ind_start + 1
             
             if (N_tmp>1) then
                
                tmp_species(1:N_tmp) = Observations_l(ind_start:ind_end)%species
                
                ! get index vector for sorting by species (see NOTES above!)
                
                call nr_indexx( N_tmp, real(tmp_species(1:N_tmp)), indx(1:N_tmp) )
                
                ! apply sort by species
                
                indx(1:N_tmp) = indx(1:N_tmp) + ind_start - 1    ! add offset
                
                Observations_l(ind_start:ind_end) = Observations_l(indx(1:N_tmp))

             end if

             ! re-initialize
             
             ind_start    = ii
             
             this_tilenum = this_tilenum_new
             
          end if
          
       end do
       
       ! clean up
       
       deallocate(indx)
       deallocate(tilenums)

    end if

    ! -----------------------------------------------------------------
    !
    ! check whether number of obs exceeds max number allowed
    
    if (N_obsl>N_obsl_max) then
       err_msg = 'N_obsl > N_obsl_max, too many observations'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! make sure L1C_Tb obs are not assimilated at the same time/place where
    !  corresponding disaggregated L2AP_Tb observations are assimilated
    
    call turn_off_assim_SMAP_L1CTb(N_obs_param, obs_param, N_obsl, Observations_l)
    
    ! gather information about number of observations assigned to each processor
    
#ifdef LDAS_MPI
    
    call MPI_AllGather(                     &
         N_obsl,            1, MPI_integer, &
         N_obsl_vec,        1, MPI_integer, & 
         mpicomm, mpierr ) 
    
#else
    
    N_obsl_vec(1) = N_obsl
    
#endif
    
    if (any(N_obsl_vec>0)) then
       
       found_obs_f = .true.    ! found obs somewhere in full domain
       
    else
       
       found_obs_f = .false.    
       
    end if

    ! -----------------------------------------------------------------
    
    ! TO DO: scale Observations (note: no scaling for moisture contents)
    
    ! may be needed for multi-variate assimilation, reichle, 10 Jun 2011

    !!do i=1,N_obs
    !!   
    !!   select case (Observations(i)%species)
    !!      
    !!   case (5)
    !!      
    !!      Observations(i)%obs      = Observations(i)%obs     / scale_temp
    !!      Observations(i)%obsvar   = Observations(i)%obsvar  /(scale_temp**2)
    !!      
    !!   end select
    !!   
    !!end do
    
    ! -----------------------------------------------------------------

    ! determine total number of observations to be assimilated
    !
    ! NOTE: This number may be different from the total number of obs
    !       recorded in the "obslog" file because: 
    !  (i)  the number of obs recorded in the "obslog" file is
    !       before obs from separate files may have been aggregated
    !       to "super-obs", and
    !  (ii) older obs readers may not contribute to the "obslog" file.
    
    write (tmpstr12,'(i12)') sum(N_obsl_vec)   ! convert integer to string
    
    tmpstr300 =  'collect_obs(): read N_obsf = ' //  tmpstr12 //                 &
         ' [after obs-based QC, super-obbing'
    
    if (.not. any_scaled_obs) then
       
       tmpstr300 = trim(tmpstr300) // ']'
       
    else
       
       tmpstr300 = trim(tmpstr300) // ', scaling]'
       
    end if
    
    if (logit) write (logunit,'(400A)') trim(tmpstr300)
    
    ! -------------------------------------
    
  end subroutine collect_obs
  
  
  ! **************************************************************

  subroutine put_into_Observations( this_obs_param, N_obs_max, N_catd, l2f,   &
       tmp_obs, tmp_std_obs, tmp_lon, tmp_lat, tmp_time, tmp_assim,           &
       obs_count, Observations )
    
    ! Put one type of observations into the general "Observations" vector:
    !  throw out no-data-values and keep track of which components have 
    !  been filled ("obs_count")
    ! All observations (except no-data-values) are included at this stage
    !  regardless of whether they will be assimilated or whether only
    !  innovations will be computed.
    !
    ! added "l2f" index vector that maps "local" tilenum to tilenum 
    ! within "full" domain
    ! - reichle, 17 Oct 2011
    !
    ! reichle, 31 Jan 2014: added "tmp_time"
    ! reichle,  6 Jun 2016: added flag "tmp_assim" to facilitate retaining unscaled obs
    
    implicit none
    
    type(obs_param_type), intent(in)  :: this_obs_param
    
    integer, intent(in) :: N_obs_max, N_catd 

    integer, dimension(N_catd), intent(in) :: l2f

    real,    dimension(N_catd), intent(in) :: tmp_obs, tmp_std_obs, tmp_lon, tmp_lat

    real*8,  dimension(N_catd), intent(in) :: tmp_time

    logical, dimension(N_catd), intent(in) :: tmp_assim
    
    integer,                              intent(inout) :: obs_count
    
    type(obs_type), dimension(N_obs_max), intent(inout) :: Observations  
    
    ! ------------------------------------------------------

    integer :: i
    
    real :: nodatavalue, tol
    
    ! ------------------------------------------------------
    
    nodatavalue = this_obs_param%nodata
    
    tol = abs(nodatavalue*nodata_tolfrac_generic)
    
    do i=1,N_catd
       
       if (abs(tmp_obs(i)-nodatavalue) > tol) then  ! check for no-data-value
          
          obs_count = obs_count+1             ! augment observation counter 
          
          Observations(obs_count)%obs  = tmp_obs(i)
          
          ! check if std has been set already; if not, use default value
          ! (no-data-value for std is any negative value)
          
          if (tmp_std_obs(i) > .0) then
             
             Observations(obs_count)%obsvar = tmp_std_obs(i)**2
             
          else
             
             Observations(obs_count)%obsvar = this_obs_param%errstd**2
             
          end if
          
          Observations(obs_count)%tilenum = l2f(i)

          Observations(obs_count)%time    = tmp_time(i)

          Observations(obs_count)%lat     = tmp_lat(i)
          Observations(obs_count)%lon     = tmp_lon(i)
          
          Observations(obs_count)%species = this_obs_param%species 
          
          Observations(obs_count)%assim   = this_obs_param%assim .and. tmp_assim(i)
          
          Observations(obs_count)%fcst    = this_obs_param%nodata
          Observations(obs_count)%fcstvar = this_obs_param%nodata
          
          Observations(obs_count)%ana     = this_obs_param%nodata
          Observations(obs_count)%anavar  = this_obs_param%nodata
          
       end if
    end do
    
  end subroutine put_into_Observations
  
  
  ! *****************************************************************
  
  subroutine add_to_obslog(                                                            &
       date_time_string, obs_param_descr, subroutine_name, num_obs_string, file_name )
    
    implicit none
    
    character( *), intent(in) :: date_time_string    ! format: YYYYMMDD_HHMMSSz
    character( *), intent(in) :: obs_param_descr
    character( *), intent(in) :: subroutine_name
    character( *), intent(in) :: num_obs_string
    character(*), intent(in) :: file_name
    
    ! obslog file format (comma-separated values; CSV): 
    !
    !  analysis time, obs species descriptor, subroutine name, obs_count, file name
    
    write (unitnum_obslog,'(500A)')           &
         date_time_string      // ', ' //     &
         trim(obs_param_descr) // ', ' //     &
         trim(subroutine_name) // ', ' //     &
         num_obs_string        // ', ' //     &
         trim(file_name)             
    
  end subroutine add_to_obslog
  
  ! *****************************************************************
  
  subroutine get_tile_num_for_obs(N_catd, tile_coord,                          &
       tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij, N_latlon, lat, lon,  &
       this_obs_param,                                                         &
       tile_num,                                                               &
       shift_lat, shift_lon )
    
    ! find one tile for each obs that "administers" the obs 
    !
    ! designed to work only for a vector of obs from a *single* species
    ! (to be called from within obs readers)
    !
    ! - reichle, 2015/02/20
    
    implicit none
    
    integer,                                             intent(in)  :: N_catd, N_latlon
    
    type(tile_coord_type), dimension(:),                 pointer     :: tile_coord ! input
    
    type(grid_def_type),                                 intent(in)  :: tile_grid    
    
    integer, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in)  :: N_tile_in_cell_ij
    
    integer, dimension(:,:,:),                           pointer     :: tile_num_in_cell_ij ! input
    
    real,    dimension(N_latlon),                        intent(in)  :: lat, lon

    type(obs_param_type),                                intent(in)  :: this_obs_param
    
    integer, dimension(N_latlon),                        intent(out) :: tile_num
    
    real,                                     optional,  intent(in)  :: shift_lat
    real,                                     optional,  intent(in)  :: shift_lon
    
    ! ------------------------
    
    ! local variables
    
    real,    dimension(N_latlon) :: max_dist_x  ! vector [deg]
    real                         :: max_dist_y  ! scalar [deg]
    
    character(len=*), parameter  :: Iam = 'get_tile_num_for_obs'
    
    real,    dimension(N_latlon) :: tmp_lon, tmp_lat
    
    ! -----------------------------------------------------------------------------
    !
    ! get "max_dist" in deg lat/lon from field-of-view (FOV) 
    !
    ! "max_dist" = Maximum distance allowed between obs lat/lon and tile com_lat/com_lon
    ! when searching for a tile to which the obs will be assigned.
    !
    ! NOTE: Subroutine get_tile_num_from_latlon() computes distances in Minkowski norm.
    
    if     ( trim(this_obs_param%FOV_units)=='deg' ) then
       
       max_dist_y = this_obs_param%FOV
       max_dist_x = this_obs_param%FOV
       
    elseif ( trim(this_obs_param%FOV_units)=='km'  ) then

       ! convert from [km] (FOV) to [deg] (max_dist_*)
       
       call dist_km2deg( this_obs_param%FOV, N_latlon, lat, max_dist_x, max_dist_y )
       
    else
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown FOV_units')
       
    end if
    
    if (max_dist_y<0. .or. any(max_dist_x<0.))  &
         call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'encountered negative max_dist')

    
    ! SMAP 36 km EASE grid readers require special accommodation:
    !
    ! temporarily shift lat/lon of obs for computation of nearest tile to
    ! avoid ambiguous assignment of M09 model tile within M36 obs grid cell
    ! (center of M36 grid cell is equidistant from at least two M09 model 
    !  tiles) -- reichle, 23 Aug 2013
    
    tmp_lat = lat
    tmp_lon = lon
    
    if (present(shift_lat))  tmp_lat = tmp_lat + shift_lat
    if (present(shift_lon))  tmp_lon = tmp_lon + shift_lon

    
    ! find tile numbers

    call get_tile_num_from_latlon(N_catd, tile_coord,             &
         tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij,       &
         N_latlon, tmp_lat, tmp_lon,                              &
         tile_num, max_dist_x, max_dist_y )

    
  end subroutine get_tile_num_for_obs

  ! *****************************************************************

end module clsm_ensupd_read_obs

#if 0

! test programs

program test

  use clsm_ensupd_read_obs
  
  implicit none
  
  integer, parameter :: N_files = 3

  character(200), parameter :: fpath = &
       '/land/l_data/AMSR/data/AMSR_E_L2_Land_V001/2002/M09/'
  character(37), dimension(N_files), parameter :: fname = (/ &
       'AMSR_E_L2_Land_B01_200209050129_A.hdf', &
       'AMSR_E_L2_Land_B01_200209091008_D.hdf', &
       'AMSR_E_L2_Land_B01_200209092002_D.hdf' /)
  
  integer :: N_data
  
  real, dimension(:), pointer :: lon, lat, ae_l2_sm
  
  character(300), dimension(N_files) :: infiles
  
  integer :: i
  
  do i=1,N_files
     
     infiles(i) = trim(fpath) // trim(fname(i))
     
  end do
  
  call read_ae_l2_sm_hdf(N_files, infiles, N_data, lon, lat, ae_l2_sm )
  
  write (*,*) 'N_data = ', N_data
  
  do i=1,3
     write (*,*) i, lon(i), lat(i), ae_l2_sm(i)
  end do
  
  do i=N_data-2,N_data
     write (*,*) i, lon(i), lat(i), ae_l2_sm(i)
  end do
  
end program test

#endif

! *******  EOF *************************************************************

  
