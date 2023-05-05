! this file contains a collection of subroutines that are needed to
! run the Ensemble Kalman filter with the catchment model off-line driver
!
! reichle, 10 Apr 01
!

module clsm_ensupd_upd_routines
  
  use nr_ran2_gasdev,                   ONLY:     &
       NRANDSEED
  
  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_TICE,                                 &
       MAPL_RADIUS,                               &
       MAPL_PI
  
  use LDAS_ensdrv_Globals,              ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_generic,                            &
       nodata_tol_generic
  
  use clsm_ensupd_glob_param,           ONLY:     &
       N_obs_species_nml,                         &
       unitnum_obslog,                            &
       scale_temp,                                &
       scale_catdef,                              &
       scale_rzexc,                               &
       scale_srfexc,                              &
       scale_ght1,                                &
       FT_ANA_FT_THRESHOLD,                       &
       FT_ANA_LOWERBOUND_ASNOW,                   &
       FT_ANA_LOWERBOUND_TEFF,                    &
       FT_ANA_UPPERBOUND_TEFF 
  
  use my_matrix_functions,              ONLY:     &
       row_variance,                              &
       unique_rows_3col
  
  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename,                           &
       is_in_rectangle

  use LDAS_DateTimeMod,                 ONLY:     &
       date_time_type
  
  use catch_types,                      ONLY:     &
       cat_param_type,                            &
       cat_progn_type,                            &
       catprogn2wesn,                             &
       catprogn2htsn,                             &
       catprogn2ghtcnt,                           &
       assignment (=),                            &
       operator (+),                              &
       operator (/)
  
  use enkf_types,                       ONLY:     &
       obs_type,                                  &
       obs_param_type,                            &
       write_obs_param,                           &
       N_obs_ang_max

  use LDAS_DriverTypes,                 ONLY:     &
       met_force_type

  use mwRTM_types,                      ONLY:     &
       mwRTM_param_type

  use LDAS_PertTypes,                   ONLY:     &
       pert_param_type,                           &
       allocate_pert_param,                       &
       deallocate_pert_param
  
  use LDAS_TileCoordType,               ONLY:     &
       tile_coord_type,                           &
       grid_def_type

  use LDAS_TilecoordRoutines,           ONLY:     &
       get_tile_num_in_ellipse,                   &
       get_number_of_tiles_in_cell_ij,            &
       get_tile_num_in_cell_ij,                   &
       get_minExtent_grid,                        &
       get_ij_ind_from_latlon

  use land_pert_routines,               ONLY:     &
       get_pert

  use mwRTM_routines,                   ONLY:     &
       mwRTM_get_Tb,                              &
       catch2mwRTM_vars
  
  use catchment_model,                  ONLY:     &
       catch_calc_tsurf,                          &
       catch_calc_tsurf_excl_snow
  
  use lsm_routines,                     ONLY:     &
       catch_calc_soil_moist,                     &
       catch_calc_tp,                             &
       catch_calc_ght,                            &
       catch_calc_FT

  use catch_constants,                  ONLY:     &
       N_snow => CATCH_N_SNOW,                    &
       N_gt   => CATCH_N_GT,                      &
       PEATCLSM_POROS_THRESHOLD

  use StieglitzSnow,                    ONLY:     &
       StieglitzSnow_calc_asnow

  use LDAS_ensdrv_mpi,                  ONLY:     &
       numprocs,                                  &
       myid,                                      &
       mpicomm,                                   &
       MPI_obs_type,                              &
       mpistatus,                                 &
       mpierr

  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR
  
  use enkf_general,                     ONLY:     &
       enkf_increments
  
  ! -----------------------------------------------------------------------
  
  implicit none

  include 'mpif.h'

  private

  public :: read_ens_upd_inputs
  public :: finalize_obslog
  public :: get_cat_progn_ens_avg
  public :: get_obs_pred
  public :: get_halo_obs
  public :: get_obs_pert
  public :: cat_enkf_increments
  public :: get_ind_obs_assim
  public :: get_ind_obs_lat_lon_box
  public :: get_halo_around_tile
  public :: TileNnzObs
  public :: dist_km2deg

  ! threshold below which FOV is considered zero (regardless of units)
  
  real, parameter, public :: FOV_threshold = 1e-4  
  
  type, public :: halo_type
     real :: minlon, minlat, maxlon, maxlat
  end type halo_type
  
contains
  
  ! ********************************************************************
 
  subroutine read_ens_upd_inputs(               &    
       work_path,                               &
       exp_id,                                  &
       date_time,                               &
       N_catf,             tile_coord_f,        &
       N_progn_pert,       progn_pert_param,    &
       N_force_pert,       force_pert_param,    &
       need_mwRTM_param,                        &
       update_type,                             &
       xcompact, ycompact,                      &
       fcsterr_inflation_fac,                   &
       N_obs_param,                             &
       obs_param,                               &
       out_obslog,                              &
       out_ObsFcstAna,                          &
       out_smapL4SMaup,                         &
       N_obsbias_max                            &
       )
    
    ! read EnKF inputs from namelist file
    !
    ! runtime options are read in three steps:
    !
    ! 1.) read options from default namelist file called 
    !      ens_upd_inputs.nml in working directory (must be present)
    !
    ! 2.) overwrite options from special namelist file (if present)
    !      specified at the command line using -ens_upd_inputs_path 
    !      and -ens_upd_inputs_file
    !
    ! reichle, 19 Jul 2005 
    ! reichle, 14 Apr 2006 - added "update_type" to namelist and outputs
    !                      - removed reading from "stored"/"saved" nml file
    ! reichle, 27 Mar 2014 - added "obslog" 

    implicit none
    
    character(*),         intent(in)    :: work_path
    character(*),         intent(in)    :: exp_id

    type(date_time_type), intent(in)    :: date_time
    
    integer,              intent(in)    :: N_catf, N_progn_pert, N_force_pert
    
    type(tile_coord_type), dimension(N_catf),       intent(in) :: tile_coord_f

    type(pert_param_type), dimension(N_progn_pert), intent(in) :: progn_pert_param
    type(pert_param_type), dimension(N_force_pert), intent(in) :: force_pert_param
    
    logical,              intent(inout) :: need_mwRTM_param

    integer,              intent(out)   :: update_type
    
    real,                 intent(out)   :: xcompact, ycompact
    real,                 intent(out)   :: fcsterr_inflation_fac

    integer,              intent(out)   :: N_obs_param
    
    type(obs_param_type), dimension(:), pointer :: obs_param     ! output
    
    logical,              intent(out)   :: out_obslog
    logical,              intent(out)   :: out_ObsFcstAna
    logical,              intent(out)   :: out_smapL4SMaup

    integer,              intent(out)   :: N_obsbias_max
    
    ! ------------------------
    
    ! locals

    ! frequency range for determining "need_mwRTM_param" 
    
    real, parameter :: min_freq =  1.e9     ! GHz   ! include L-band (SMOS, SMAP)
    real, parameter :: max_freq = 10.e9     ! GHz   ! include X-band (TRMM, AMSR-E)

    ! tolerance for checking xcorr against FOV (xcorr + tol >= "FOV") 
    
    real, parameter :: tol = 1.e-5  ! units of deg lat/lon
    
    character(300)  :: fname

    character(200)  :: ens_upd_inputs_path
    character( 40)  :: ens_upd_inputs_file, dir_name, file_tag, file_ext
    
    integer :: i, j, k, N_tmp, k_hD, k_hA, k_vD, k_vA

    real    :: r_y
    
    real, dimension(1) :: tmp_lat, r_x
    
    logical :: smap_species, smos_species

    type(obs_param_type), dimension(N_obs_species_nml) :: obs_param_nml

    character(len=*), parameter :: Iam = 'read_ens_upd_inputs'
    character(len=400) :: err_msg
    character(len=  6) :: tmpstring6
    logical :: file_exists

    ! -----------------------------------------------------------------
    
    namelist /ens_upd_inputs/      &
         update_type,              &
         out_obslog,               &
         out_ObsFcstAna,           &
         out_smapL4SMaup,          &
         xcompact, ycompact,       &
         fcsterr_inflation_fac,    &
         obs_param_nml
        
    ! ------------------------------------------------------------------
    !
    ! Set default file name for EnKF inputs namelist file
    
    ens_upd_inputs_path = '.'                                       ! set default 
    ens_upd_inputs_file = 'LDASsa_DEFAULT_inputs_ensupd.nml'
    
    ! Read data from default ens_upd_inputs namelist file 
    
    fname = trim(ens_upd_inputs_path) // '/' // trim(ens_upd_inputs_file)
    
    open (10, file=fname, delim='apostrophe', action='read', status='old')
    
    if (logit) write (logunit,*)
    if (logit) write (logunit,'(400A)') 'reading *default* EnKF inputs from ' // trim(fname)
    if (logit) write (logunit,*)

    read (10, nml=ens_upd_inputs)
    
    close(10,status='keep')
    
    
    ! Get name and path for special EnKF inputs file from
    ! command line (if present) 
    
    ens_upd_inputs_path = '.'
    ens_upd_inputs_file = 'LDASsa_SPECIAL_inputs_ensupd.nml'
       
    ! Read data from special EnKF inputs namelist file 
       
    fname = trim(ens_upd_inputs_path)//'/'//trim(ens_upd_inputs_file)
    inquire(file=fname, exist=file_exists)

    if(file_exists) then
 
       open (10, file=fname, delim='apostrophe', action='read', status='old')
          
       if (logit) write (logunit,*)
       if (logit) write (logunit,'(400A)') 'reading *special* EnKF inputs from ' // trim(fname)
       if (logit) write (logunit,*)

       read (10, nml=ens_upd_inputs)

       close(10,status='keep')
       
    end if
    
    ! ---------------------------------
    
    ! overwrite EnKF inputs with command line options, if any
    !
    ! none implemented so far (reichle, 19 Jul 2005)
    
    ! -----------------------------------------------------------------
    !
    ! consistency checks etc
    
    if (update_type==0) then
       err_msg = 'executable was built for assimilation but update_type=0'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    N_obsbias_max = 0  ! initialize

    do i=1,N_obs_species_nml
       
       ! make sure "getinnov" is .true. if innovations are needed for state 
       !  or obs bias updates
       
       if ( obs_param_nml(i)%assim .or. obs_param_nml(i)%bias_Npar>0 )       &
            obs_param_nml(i)%getinnov = .true.
       
       ! determine maximum "bias_Npar"
       
       N_obsbias_max = max( N_obsbias_max, obs_param_nml(i)%bias_Npar )
       
    end do

    ! -----------------------------------------------------------------
    !
    ! Extract only species of interest (i.e., %getinnov=.true.) from nml inputs:
    !
    ! NOTE: multi-angular obs (eg, SMOS Tb h-pol ascending) are defined 
    !       as *one* species obs_param_nml in namelist file and are
    !       split here into a new set of species, each with a unique incidence angle

    ! first loop: count number of species of interest (those with %getinnov=.true.)

    j = 0
    
    do i=1,N_obs_species_nml
       
       if (obs_param_nml(i)%getinnov) then
          
          N_tmp = max( obs_param_nml(i)%N_ang, 1 )   ! some species have N_ang=0
          
          j = j + N_tmp
          
       end if
       
    end do
    
    N_obs_param = j
    
    allocate(obs_param(N_obs_param))
    
    ! second loop: extract species of interest
    
    j = 0
    
    do i=1,N_obs_species_nml
       
       if (obs_param_nml(i)%getinnov) then
          
          ! check for consistency between "varname" and "RTMid"
          
          if ( (trim(obs_param_nml(i)%varname)=='Tb') .and.          &
               (     obs_param_nml(i)%RTM_ID  ==  0 )        ) then
             
             write (tmpstring6,*) i
             
             err_msg = 'inconsistent obs_param_nml%varname and obs_param_nml%RTM_ID '
             err_msg = trim(err_msg) // 'in ensupd nml inputs for species=' // tmpstring6
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             
          end if
          
          ! extract species into obs_param
          
          N_tmp = max( obs_param_nml(i)%N_ang, 1 )
          
          do k=1,N_tmp
             
             j = j + 1
             
             obs_param(j)         = obs_param_nml(i)         ! inherit everything
             
             obs_param(j)%species = j                        ! provide unique species ID
             
             obs_param(j)%N_ang   = 1                        ! overwrite N_ang
             
             obs_param(j)%ang(1)  = obs_param_nml(i)%ang(k)  ! overwrite ang(1)
             
             obs_param(j)%ang(2:N_obs_ang_max) = nodata_generic  ! fill rest with nodata

          end do
          
       end if
       
    end do
    
    ! -----------------------------------------------------------------
    !
    ! echo variables of ens_upd_inputs
    
    if (logit) write (logunit,*) 'EnKF inputs are:'
    if (logit) write (logunit,*)
    if (logit) write (logunit, nml=ens_upd_inputs) 
    if (logit) write (logunit,*)

    ! -----------------------------------------------------------------
    !
    ! more consistency checks (only done on species of interest)
    
    do i=1,N_obs_param
       
       ! check xcorr, ycorr (spatial correlation scale of obs error) 
       !  against some measure of FOV
       
       ! get FOV in units of [deg] 
       
       if     ( trim(obs_param(i)%FOV_units)=='deg' ) then
          
          r_x(1) = obs_param(i)%FOV
          r_y    = obs_param(i)%FOV
          
       elseif ( trim(obs_param(i)%FOV_units)=='km'  ) then
          
          ! compute FOV in units of [deg] at (area-weighted) average abs latitude of tiles 
          
          tmp_lat(1) = sum( abs(tile_coord_f%com_lat) * tile_coord_f%area )/sum(tile_coord_f%area)
          
          call dist_km2deg( obs_param(i)%FOV, 1, tmp_lat, r_x, r_y )
          
       else
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown FOV_units')
          
       end if
       
       ! enforce (xcorr >= "FOV") and (ycorr >= "FOV")  (with some tolerance)
       
       if ( ( obs_param(i)%xcorr < (r_x(1) - tol) )  .or.            &
            ( obs_param(i)%ycorr < (r_y    - tol) )         ) then
          
          if (logit) write (logunit,*) 'i                  = ', i
          if (logit) write (logunit,*) 'obs_param(i)%xcorr = ', obs_param(i)%xcorr
          if (logit) write (logunit,*) 'obs_param(i)%ycorr = ', obs_param(i)%ycorr
          if (logit) write (logunit,*) 'r_x(1)             = ', r_x(1)
          if (logit) write (logunit,*) 'r_y                = ', r_y
          
          err_msg = 'found xcorr<avg(FOV) or ycorr<avg(FOV) for ' //   &
               'obs_param%descr = ' // trim(obs_param(i)%descr)
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    end do
        
    ! check correlation lengths of perturbations against compact
    !  support length scale
    
    call check_compact_support(                                      &
         update_type,                                                &
         N_catf,             tile_coord_f,                           &
         N_progn_pert,       progn_pert_param,                       &
         N_force_pert,       force_pert_param,                       &
         N_obs_param,        obs_param,                              &
         xcompact, ycompact )
    

    ! check whether a microwave radiative transfer model (mwRTM)
    !  is needed for innovations (obs-minus-forecast residuals
    !  or assimilation
    !
    ! only check if "need_mwRTM_param=.false." on input
    
    if (.not. need_mwRTM_param) then 
       
       do i=1,N_obs_param
          
          if ( (trim(obs_param(i)%varname) ==     'Tb' )  .and.            &
               (     obs_param(i)%getinnov             )  .and.            &
               (     obs_param(i)%freq     >= min_freq )  .and.            &
               (     obs_param(i)%freq     <= max_freq )         )   then
             
             need_mwRTM_param = .true.
             
          end if
          
       end do
       
    end if

    ! when L4SMaup files are written, ensure that species of interest do not
    ! simultaneously include "SMAP_L*_Tb*" and "SMOS_fit_Tb*" obs
    
    if (out_smapL4SMaup) then
       
       smap_species = .false.
       smos_species = .false.
       
       do i=1,N_obs_param
          
          select case (trim(obs_param(i)%descr))
             
          case('SMAP_L2AP_Tbh_D',    'SMAP_L2AP_Tbv_D',    &
               'SMAP_L2AP_Tbh_A',    'SMAP_L2AP_Tbv_A',    &
               'SMAP_L1C_Tbh_D', 'SMAP_L1C_Tbv_D',     &
               'SMAP_L1C_Tbh_A',     'SMAP_L1C_Tbv_A',     &
               'SMAP_L1C_Tbh_E09_D', 'SMAP_L1C_Tbv_E09_D', &
               'SMAP_L1C_Tbh_E09_A', 'SMAP_L1C_Tbv_E09_A', & 
               'SMAP_L1C_Tbh_E27_D', 'SMAP_L1C_Tbv_E27_D', &
               'SMAP_L1C_Tbh_E27_A', 'SMAP_L1C_Tbv_E27_A'  & 
               )
             
             smap_species = .true.
             
             
          case('SMOS_fit_Tbh_D', 'SMOS_fit_Tbv_D',     &
               'SMOS_fit_Tbh_A', 'SMOS_fit_Tbv_A'      &
               )
             
             smos_species = .true.
             
          case default
             
             ! do nothing
             
          end select
          
       end do
       
       ! stop if SMAP and SMOS species are present simultaneously
       
       if (smap_species .and. smos_species) then
          
          err_msg = 'out_smapL4SMaup=.true. is *not* compatible with ' // &
               'simultaneously using "SMAP_L*_Tb*" and "SMOS_fit_Tb*"'  // &
               'obs species. Use "out_ObsFcstAna" or remove either' // &
               'SMAP or SMOS obs species.'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    end if
    
    ! stop if two or more SMAP L1C Tb sub-species of same polarization and orbit
    !  direction are assimilated (it is ok to have more than one as passive obs)
    
    k_hD=0
    k_hA=0
    k_vD=0
    k_vA=0

    ! count number of assimilated L1C Tb sub-species (per polarization and orbit dir)    

    do i=1,N_obs_param         
       
       if (obs_param(i)%assim) then
          
          select case (trim(obs_param(i)%descr))
             
          case('SMAP_L1C_Tbh_D','SMAP_L1C_Tbh_E09_D','SMAP_L1C_Tbh_E27_D'); k_hD=k_hD+1 
          case('SMAP_L1C_Tbh_A','SMAP_L1C_Tbh_E09_A','SMAP_L1C_Tbh_E27_A'); k_hA=k_hA+1 
          case('SMAP_L1C_Tbv_D','SMAP_L1C_Tbv_E09_D','SMAP_L1C_Tbv_E27_D'); k_vD=k_vD+1 
          case('SMAP_L1C_Tbv_A','SMAP_L1C_Tbv_E09_A','SMAP_L1C_Tbv_E27_A'); k_vA=k_vA+1 
             
          end select
          
       end if
          
    end do
    
    if (k_hD>1 .or. k_hA>1 .or. k_vD>1 .or. k_vA>1) then
       
       err_msg =                                                            &
            'for given polarization and orbit dir must not assimilate ' //  &
            'more than one of L1C_Tb, L1C_Tb E09, L1C_Tb E27'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    

    
    ! -------------------------------------------------------------
    !
    ! save ens update inputs into *ens_upd_inputs.nml file

    dir_name = 'rc_out'
    file_tag = 'ldas_ensupd_inputs'
    file_ext = '.nml'
    
    fname = get_io_filename( work_path, exp_id, file_tag, date_time=date_time, &
         dir_name=dir_name, file_ext=file_ext )
    
    if (logit) write (logunit,'(400A)') 'writing ens upd inputs to ' // trim(fname)
    if (logit) write (logunit,*)
    
    open (10, file=fname, status='new', action='write', delim='apostrophe' )
    
    write(10, nml=ens_upd_inputs)
    
    close(10, status='keep')

    ! -------------------------------------------------------------
    !
    ! save obs_param into *obsparam.txt file

    dir_name = 'rc_out'
    file_tag = 'ldas_obsparam'
    file_ext = '.txt'
    
    fname = get_io_filename( work_path, exp_id, file_tag, date_time=date_time, &
         dir_name=dir_name, file_ext=file_ext )
    
    if (logit) write (logunit,'(400A)') 'writing obs parameters to ' // trim(fname)
    if (logit) write (logunit,*)

    open (10, file=fname, status='new', action='write', delim='apostrophe' )
    
    call write_obs_param( 10, N_obs_param, obs_param)
    
    close(10, status='keep')
    
    ! -------------------------------------------------------------------------
    !
    ! if requested, open obslog file and write header
    
    if (out_obslog)  call init_obslog( work_path, exp_id, date_time )
    
    ! -------------------------------------------------------------
    
  end subroutine read_ens_upd_inputs
  
  ! ********************************************************************

  subroutine init_obslog( work_path, exp_id, date_time )
    
    ! open obslog file and write header
    
    implicit none
    
    character(*),       intent(in)    :: work_path
    character(*),       intent(in)    :: exp_id

    type(date_time_type), intent(in)    :: date_time
    
    ! local variables
    
    character(300) :: fname
    character( 40) :: dir_name, file_tag, file_ext
    
    integer        :: istat
    
    logical        :: is_open
    
    character(len=*), parameter :: Iam = 'init_obslog'
    character(len=400) :: err_msg

    ! ----------------------------------------------------------------------------
    
    dir_name = 'rc_out'
    file_tag = 'ldas_obslog' 
    file_ext = '.txt'
    
    fname = get_io_filename( work_path, exp_id, file_tag, date_time=date_time, &
         dir_name=dir_name, file_ext=file_ext )
    
    ! make sure "unitnum_obslog" is not already in use
    
    inquire( unit=unitnum_obslog, opened=is_open )
    
    if (is_open) then
       err_msg = '"unitnum_obslog" is taken, edit src code ' // &
            'to automatically detect suitable "unitnum_obslog"'
    end if

    ! open "obslog" file and write header lines
    
    open( unitnum_obslog, file=fname, form='formatted', action='write', iostat=istat)
    
    if (istat/=0) then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'cannot open "obslog" file')
    end if
    
    if (logit) write (logunit,'(400A)') 'writing header of obslog file: ' // trim(fname)
    if (logit) write (logunit,*)

    ! header line 1 -- objective and guidance
    
    write (unitnum_obslog,'(200A)')                                               & 
         'Obs log for EnKF analysis, includes obs counts and file names '      // &
         'from which obs have been read. Obs count is *after* initial,  '      // &
         'obs-based QC and *before* model-based QC. '
    
    ! header line 2 -- warnings
    
    write (unitnum_obslog,'(200A)') & 
         'IMPORTANT:  Provides log ONLY for select obs species/readers!  '     // &
         'SMAP L1C_Tb obs have not yet been cross-masked against L2AP_Tb obs!'
    
    ! header line 3 -- file format
    
    write (unitnum_obslog,'(200A)') & 
         'Format (comma-separated values; last line: "EOF"): '                 // & 
         'EnKF analysis time [YYYYMMDD_HHMMSSz], obs species descriptor, '     // &
         'subroutine name, obs_count, file name'
    
    ! obslog file remains open (similar to "logunit")
        
  end subroutine init_obslog

  ! ********************************************************************
  
  subroutine finalize_obslog()
    
    ! finalize and close "obslog" file
    
    implicit none
    
    write (unitnum_obslog,'(3A)') 'EOF'
    
    close (unitnum_obslog, status='keep')
    
    if (logit) write(logunit,*) 'done writing obslog file'
    
  end subroutine finalize_obslog

  ! ********************************************************************
  
  subroutine get_cat_progn_ens_avg(N_catd, N_ens, cat_progn, cat_progn_ensavg)
    
    implicit none
    
    integer, intent(in) :: N_catd, N_ens
    
    type(cat_progn_type), dimension(N_catd,N_ens), intent(in) :: cat_progn
    
    type(cat_progn_type), dimension(N_catd), intent(out) :: cat_progn_ensavg
    
    ! locals
    
    integer :: i, n_e
    
    ! -------------------------------------
    
    do i=1,N_catd
       
       cat_progn_ensavg(i) = 0.
       
       do n_e=1,N_ens
          
          cat_progn_ensavg(i) = cat_progn_ensavg(i) + cat_progn(i,n_e)
          
       end do
       
       cat_progn_ensavg(i) = cat_progn_ensavg(i)/real(N_ens)
       
    end do
    
  end subroutine get_cat_progn_ens_avg
  
  ! *********************************************************************
  
  ! subroutine recompute_diagnostic( )
  !
  ! moved to clsm_ensdrv_drv_routines.F90
  ! -reichle+csdraper, 30 Oct 2013
  !
  ! end subroutine recompute_diagnostic

  ! ********************************************************************
  
  subroutine assemble_obs_cov(N_obs, N_obs_param, obs_param, Observations, Obs_cov)
    
    ! assemble measurements error covariance
    
    ! reichle, 27 Jul 2005
    
    implicit none
    
    integer, intent(in) :: N_obs, N_obs_param
    
    type(obs_param_type), dimension(N_obs_param), intent(in) :: obs_param

    type(obs_type), dimension(N_obs), intent(in) :: Observations  
    
    real, intent(out), dimension(N_obs,N_obs) :: Obs_cov
    
    ! -------------------------------------------------------------
    
    ! locals
    
    integer :: i, j, i_species, j_species   !! inum, jnum
    
    real :: fac, xcorr_tmp, ycorr_tmp
    
    ! -------------------------------------------------------------

    if (N_obs==0) return
    
    ! assemble measurement error covariance 
    
    ! initialize
    
    Obs_cov = 0.
    
    ! diagonal elements

    do i=1,N_obs
       
       Obs_cov(i,i) = Observations(i)%obsvar
       
    end do
    
    ! off-diagonal elements
    
    do i=1,N_obs
       do j=(i+1),N_obs
          
          i_species = Observations(i)%species
          j_species = Observations(j)%species
          
          ! have non-zero correlation only between observations of same type
          
          if (i_species == j_species) then
             
             xcorr_tmp = obs_param(i_species)%xcorr
             ycorr_tmp = obs_param(i_species)%ycorr
             
             ! check for zero correlation distance 
             
             if (xcorr_tmp>0. .and. ycorr_tmp>0.) then
                
                ! compute correlation between observation locations
                
                !!inum = Observations(i)%tilenum
                !!jnum = Observations(j)%tilenum 
                
                ! compute Gaussian correlation
                
                !!fac =  & 
                !!  ((tile_coord(inum)%com_lon-tile_coord(jnum)%com_lon)**2 &
                !!  /xcorr_tmp**2 )                                         &
                !!  +                                                       &
                !!  ((tile_coord(inum)%com_lat-tile_coord(jnum)%com_lat)**2 & 
                !!  /ycorr_tmp**2 )               
                
                fac =  & 
                     ((Observations(i)%lon-Observations(j)%lon)**2 &
                     /xcorr_tmp**2 )                                         &
                     +                                                       &
                     ((Observations(i)%lat-Observations(j)%lat)**2 & 
                     /ycorr_tmp**2 )               

                fac = exp(-.5*fac)

                ! bug fix
                ! GDL+reichle, 17 Oct 2014
                !Obs_cov(i,j) = Observations(i)%obsvar * fac
                Obs_cov(i,j) = sqrt(Observations(i)%obsvar * Observations(j)%obsvar) * fac
                
                Obs_cov(j,i) = Obs_cov(i,j)
                
             end if
          end if
          
       end do
    end do
        
  end subroutine assemble_obs_cov
  
  ! *********************************************************************
  
  subroutine get_obs_pred(                                       &
       beforeEnKFupdate,                                         &
       N_obs_param, N_ens,                                       &
       N_catl, tile_coord_l,                                     &
       N_catf, tile_coord_f, f2l,                                &
       N_catl_vec, low_ind, tile_grid_g,                         &
       obs_param,                                                &
       met_force, lai, cat_param, cat_progn, mwRTM_param,        &
       N_obsl, Observations_l, Obs_pred_l, obsbias_ok,           &
       fcsterr_inflation_fac )
    
    ! Compute ensemble of measurement predictions from ensemble 
    !  of tile-space Catchment prognostics.
    ! Return only those Observations that pass model-based QC.
    !
    ! Overview:
    ! 1.) Determine which diagnostics are needed.
    ! 2.) Compute diagnostics for local domain, apply model-based QC.
    ! 3.) Get diagnostics from processors within halo. 
    ! 4.) Aggregate from tile space to obs space.
    ! 5.) Deal with no-data-values, compute ens mean and var of Obs_pred.
    !
    !
    ! 27 Jul 2005
    !  7 Jul 2006 - reichle: use temporary dimension(1) vectors for Absoft 
    ! 13 Jun 2011 - reichle: re-structured for MPI
    !                        added field-of-view (FOV) for observations
    !                        moved model-based QC here (from read_obs())
    !  1 Dec 2011 - reichle: added QC for Tb vs. model *soil* temp (RFI-motivated)
    ! 18 Jun 2012 - reichle: rewritten for better memory management w/ MPI
    ! 26 Mar 2014 - reichle: apply all model-based QC only before EnKF update
    ! 25 Sep 2020 - wjiang+reichle: accommodate processors that have no tiles
    
    implicit none
    
    logical,                intent(in)                             :: beforeEnKFupdate
    
    integer,                intent(in)                             :: N_obs_param, N_ens
    integer,                intent(in)                             :: N_catl, N_catf
    
    type(tile_coord_type),  dimension(:),     pointer :: tile_coord_l          ! input
    type(tile_coord_type),  dimension(:),     pointer :: tile_coord_f          ! input
    
    integer,                intent(in),    dimension(numprocs)     :: N_catl_vec, low_ind

    type(grid_def_type),    intent(in)                             :: tile_grid_g

    type(obs_param_type),   intent(in),    dimension(N_obs_param)  :: obs_param

    integer,                intent(in),    dimension(N_catf)       :: f2l
    type(met_force_type),   intent(in),    dimension(N_catl)       :: met_force
    real,                   intent(in),    dimension(N_catl)       :: lai
    type(cat_param_type),   intent(in),    dimension(N_catl)       :: cat_param
    type(cat_progn_type),   intent(in),    dimension(N_catl,N_ens) :: cat_progn
    type(mwRTM_param_type), intent(in),    dimension(N_catl)       :: mwRTM_param

    integer,                intent(inout)                          :: N_obsl   ! InOut !!!

    type(obs_type),         dimension(:),     pointer :: Observations_l        ! InOut
    
    real,                   dimension(:,:),   pointer :: Obs_pred_l            ! output
    
    logical,                intent(in),    dimension(N_obsl), optional :: obsbias_ok       

    real,                   intent(in),                       optional :: fcsterr_inflation_fac      
    
    ! --------------------------------------------------------------------------------
    !
    ! locals
    
    real,    parameter                      :: Tbobs_minus_stemp_max = 5.   ! [K]
    
    real,    parameter                      :: fac_search_FOV_km     = 2.   ! [-]
    
    real,    parameter                      :: EASE_max_water_frac   = 0.05 ! [-]

    integer                                 :: N_catlH, n_e, i, j, k, N_tmp, ii, jj
    integer                                 :: N_fields, N_Tbspecies, N_TbuniqFreqAngRTMid
    integer                                 :: this_species, this_tilenum, this_pol
    integer                                 :: this_Tbspecies, this_TbuniqFreqAngRTMid, RTM_id
    integer                                 :: istart, iend

    real                                    :: this_lon, this_FOV, r_y
    real, dimension(1)                      :: this_lat, r_x
    
    real                                    :: freq, inc_angle

    real, dimension(numprocs)               :: xhalo, yhalo, tmplatvec, tmprx
    
    real                                    :: tmpreal, tmp_stemp, tmp_fraccell

    logical                                 :: tmpRFI, tmpWater, use_distance_weights

    real, dimension(1)                      :: tmpmean, tmpvar
        
    logical                                 :: get_sfmc_l,   get_sfmc_lH 
    logical                                 :: get_rzmc_l,   get_rzmc_lH
    logical                                 :: get_tsurf_l,  get_tsurf_lH 
    logical                                 :: get_tp_l
    logical                                 :: get_Tb_l,     get_Tb_lH
    logical                                 :: get_FT_l,     get_FT_lH
    
    type(grid_def_type)                     :: tile_grid_lH        
    
    integer, dimension(N_obs_param)         :: ind_obsparam2Tbspecies
    integer, dimension(N_obs_param)         :: ind_Tbspecies2TbuniqFreqAngRTMid
    
    real,    dimension(N_obs_param,3)       :: Tb_freq_ang_RTMid
    
    ! dimension "N_catl"
    
    real,    dimension(N_catl)              :: srfexc,   rzexc,    catdef,   prmc_l
    
    real,    dimension(N_catl)              :: ar1_l,    ar2_l,    ar4_l 

    real,    dimension(N_catl)              :: asnow, tsurf_excl_snow
    
    real,    dimension(N_gt,N_catl)         :: tp_l  ! NOTE dims: N_gt-by-N_catl
                                                     !  for consistency w/ calc_tp

    real,    dimension(N_catl)              :: Tb_h_vec, Tb_v_vec
    
    real,    dimension(N_catl)              :: precip, SWE, smoist

    ! dimension "N_catl-[by-N_XXXX-]by-N_ens"  
    ! (need to be communicated between processors)
    
    real,    dimension(N_catl,N_ens)        :: sfmc_l,   rzmc_l
    real,    dimension(N_catl,N_ens)        :: tsurf_l,  stemp_l
    real,    dimension(N_catl,N_ens)        :: FT_l
    
    real,    dimension(:,:,:), allocatable  :: Tb_h_l, Tb_v_l

    real,    dimension(:,:,:), allocatable  :: tile_data_l
    
    ! dimension "N_catlH"

    type(tile_coord_type), dimension(:), pointer :: tile_coord_lH => null()
    
    integer, dimension(:),     allocatable  :: ind_tmp

    real,    dimension(:),     allocatable  :: tmp_ndst2, tmp_wFOV, tmp_weights, tmp_data
                                                             
    ! dimension "N_catlH-[by-N_XXXX-]by-N_ens"  
    ! (need to be communicated between processors)
    
    real,    dimension(:,:),   allocatable  :: sfmc_lH,   rzmc_lH
    real,    dimension(:,:),   allocatable  :: tsurf_lH,  stemp_lH
    real,    dimension(:,:),   allocatable  :: FT_lH
    
    real,    dimension(:,:,:), allocatable  :: Tb_h_lH, Tb_v_lH
    
    real,    dimension(:,:,:), pointer      :: tile_data_lH => null()
    
    ! dimension "N_catlH-by-OTHER"
    
    integer, dimension(:,:),   pointer      :: N_tile_in_cell_ij_lH   => null()
    integer, dimension(:,:,:), pointer      :: tile_num_in_cell_ij_lH => null()

    ! dimension "N_obsl" (as in N_obsl upon input)
    
    logical, dimension(N_obsl) :: obsbias_ok_tmp

    real                       :: inflation_factor

    character(len=*), parameter :: Iam = 'get_obs_pred'
    character(len=400) :: err_msg
    character(len= 10) :: tmpstring10
    
    ! --------------------------------------------------------------
    !
    ! allocate and initialize
    
    allocate(Obs_pred_l(N_obsl,N_ens))

    if (N_catl == 0) return   ! return if processor has no tiles

    if (N_obsl > 0) Obs_pred_l = nodata_generic

    ! deal with optional arguments
    
    if (present(obsbias_ok)) then
       
       obsbias_ok_tmp = obsbias_ok
       
    else
       
       obsbias_ok_tmp = .false.

    end if
        
    if (present(fcsterr_inflation_fac) .and. beforeEnKFupdate) then

       ! ONLY inflate *before* EnKF update!!!

       inflation_factor = fcsterr_inflation_fac
       
    else
       
       inflation_factor = -9999.
       
    end if

    ! --------------------------------------------------------------
    !
    ! determine which diagnostics are needed (based on obs_param because
    ! observations on local proc may not reflect all obs)

    ! get_*_l : may include additional fields needed to compute observed fields
    
    get_sfmc_l   = .false. 
    get_rzmc_l   = .false. 
    get_tsurf_l  = .false.
    get_tp_l     = .false.
    get_FT_l     = .false.
    get_Tb_l     = .false.
    
    ! get_*_lH : directly match observed fields
    
    get_sfmc_lH  = .false. 
    get_rzmc_lH  = .false. 
    get_tsurf_lH = .false.
    get_FT_lH    = .false.
    get_Tb_lH    = .false.

    ! loop through obs_param b/c obs on local proc may not reflect all obs
    
    ind_obsparam2Tbspecies = -999
    
    j = 0
    
    do i=1,N_obs_param  
       
       select case (trim(obs_param(i)%varname))
          
       case ('sfmc')
          
          get_sfmc_l   = .true.
          get_sfmc_lH  = .true.
          get_tsurf_l  = .true.    ! needed for model-based QC
          
       case ('rzmc')
          
          get_rzmc_l   = .true.
          get_rzmc_lH  = .true.
          get_tsurf_l  = .true.    ! needed for model-based QC
          
       case ('tsurf')
          
          get_tsurf_l  = .true.
          get_tsurf_lH = .true.
          get_tp_l     = .true.    ! needed for model-based QC
          get_sfmc_l   = .true.    ! needed to get ar1, ar2, and ar4

       case ('FT')
          
          get_FT_l     = .true.
          get_FT_lH    = .true.
          get_sfmc_l   = .true.    ! needed to get ar1, ar2, and ar4
          get_tp_l     = .true.    ! needed as input to calc_FT
          
       case ('Tb')
          
          j=j+1        ! count number of Tb species
          
          ind_obsparam2Tbspecies(i) = j
          
          Tb_freq_ang_RTMid(j,1) =      obs_param(i)%freq
          Tb_freq_ang_RTMid(j,2) =      obs_param(i)%ang(1)
          Tb_freq_ang_RTMid(j,3) = real(obs_param(i)%RTM_ID)
          
          get_sfmc_l   = .true. 
          get_tsurf_l  = .true.
          get_tp_l     = .true.          
          get_Tb_l     = .true.
          
          get_Tb_lH    = .true.
          
       case default
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown obs_param%varname')
          
       end select
       
    end do
    
    N_Tbspecies = j

    ! determine unique combinations of Tb frequency, angle, and RTM_ID
    !
    ! Step 1:
    ! determine unique combinations of Tb frequency and angle
    ! (obs_param usually has separate species for H- and V-pol and
    !  for ascending and descending orbits, but the mwRTM model
    !  always provides both polarizations and does not depend on the 
    !  orbit direction --> avoid computing and communicating redundant
    !  information)

    call unique_rows_3col(                                                        &
         N_Tbspecies, Tb_freq_ang_RTMid(1:N_Tbspecies,:),                         &
         N_TbuniqFreqAngRTMid, ind_Tbspecies2TbuniqFreqAngRTMid(1:N_Tbspecies) )
        
    if (get_Tb_l)  allocate(Tb_h_l(N_catl,N_TbuniqFreqAngRTMid,N_ens))
    if (get_Tb_l)  allocate(Tb_v_l(N_catl,N_TbuniqFreqAngRTMid,N_ens))
    
    ! -------------------------

    ! determine xhalo, yhalo in units of [deg] based on some measure of FOV
    
    xhalo = 0.   ! initialize
    yhalo = 0.   ! initialize
    
    ! for FOV_units in 'km', all processors need to know the xhalo of each processor, 
    !  which in turn depends on latitude
    
    tmplatvec = 0.

    do jj=1,numprocs
 
       if (N_catl_vec(jj) <= 0) cycle    ! nothing to do for this processor
       
       istart = low_ind(jj)
       iend   = istart + N_catl_vec(jj) - 1
       
       ! largest abs(lat) will have largest FOV
       
       tmplatvec(jj) = maxval( abs( tile_coord_f(istart:iend)%com_lat ))  
       
    end do
    
    ! find maximum FOV in units of [deg] across all obs params 
    
    do ii=1,N_obs_param
       
       if     ( trim(obs_param(ii)%FOV_units)=='deg' ) then
          
          xhalo = max( xhalo, obs_param(ii)%FOV )
          yhalo = max( yhalo, obs_param(ii)%FOV )
          
       elseif ( trim(obs_param(ii)%FOV_units)=='km'  ) then
          
          ! convert from [km] (FOV) to [deg] 
          
          call dist_km2deg( obs_param(ii)%FOV, numprocs, tmplatvec, tmprx, r_y )

          ! for now, ignore what happens to xhalo for processors without tiles (fixed below)
          
          xhalo = max( xhalo, tmprx )
          yhalo = max( yhalo, r_y   )
          
       else
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown FOV_units (i)')
          
       end if
       
    end do

    where (N_catl_vec<=0)  xhalo = 0.    ! set xhalo=0. for processors without tiles

    ! FOV is *radius*, leave some room

    xhalo = 2.5 * xhalo
    yhalo = 2.5 * yhalo

    ! ------------------------------------------------------------------
    !
    ! compute required diagnostics for locally managed tiles

    precip = met_force%Rainf + met_force%Snowf     ! for model-based QC

    do n_e=1,N_ens
       
       ! compute observed fields in tile space (for local domain)
       
       ! total SWE (for model-based QC)
       
       do k=1,N_catl
          
          SWE(k) = sum( cat_progn(k,n_e)%wesn(1:N_snow) )
          
       end do
       
       if (get_sfmc_l .or. get_rzmc_l) then
          
          srfexc = cat_progn(:,n_e)%srfexc
          rzexc  = cat_progn(:,n_e)%rzexc
          catdef = cat_progn(:,n_e)%catdef

          ! updated to new interface - reichle, 3 Apr 2012
          
          call catch_calc_soil_moist(                                &
               N_catl,           cat_param%dzsf,   cat_param%vgwmax, &
               cat_param%cdcr1,  cat_param%cdcr2,  cat_param%psis,   &
               cat_param%bee,    cat_param%poros,  cat_param%wpwet,  &
               cat_param%ars1,   cat_param%ars2,   cat_param%ars3,   &
               cat_param%ara1,   cat_param%ara2,   cat_param%ara3,   &
               cat_param%ara4,   cat_param%arw1,   cat_param%arw2,   &
               cat_param%arw3,   cat_param%arw4,                     &
               cat_param%bf1,    cat_param%bf2,                      &
               srfexc, rzexc, catdef,                                &
               ar1_l,  ar2_l,  ar4_l,                                &
               sfmc_l(:,n_e), rzmc_l(:,n_e), prmc_l )
          
       end if
       
       if (get_tsurf_l) then
          
          ! updated to new interface, 
          ! need ar1, ar2, ar4 from call to catch_calc_soil_moist() above
          ! - reichle, 3 Apr 2012

          call catch_calc_tsurf( N_catl,                                         &
               cat_progn(:,n_e)%tc1, cat_progn(:,n_e)%tc2, cat_progn(:,n_e)%tc4, &
               catprogn2wesn(N_catl,cat_progn(:,n_e)),                           &
               catprogn2htsn(N_catl,cat_progn(:,n_e)),                           &
               ar1_l,  ar2_l,  ar4_l,                                            &
               tsurf_l(:,n_e) )
                    
       end if
       
       if (get_tp_l) then
          
          ! NOTE: "tp" is returned in CELSIUS [for consistency w/ catchment.F90]

          ! updated to new interface - reichle, 3 Apr 2012
          
          call catch_calc_tp( N_catl, cat_param%poros,                  &
               catprogn2ghtcnt(N_catl,cat_progn(:,n_e)), tp_l )
                    
       end if
                    
       if (get_FT_l) then
          
          call StieglitzSnow_calc_asnow( N_snow, N_catl,                         &
               catprogn2wesn(N_catl,cat_progn(:,n_e)),                           &
               asnow )
          
          call catch_calc_tsurf_excl_snow( N_catl,                               &
               cat_progn(:,n_e)%tc1, cat_progn(:,n_e)%tc2, cat_progn(:,n_e)%tc4, &
               ar1_l, ar2_l, ar4_l, tsurf_excl_snow )
          
          ! catch_calc_FT() expects "tp" in CELSIUS
          
          call catch_calc_FT( N_catl, asnow, tp_l(1,:), tsurf_excl_snow, FT_l(:,n_e))
          
       end if
       
       if (get_Tb_l) then
          
          ! convert Catchment model variables into inputs suitable for the mwRTM 
          
          call catch2mwRTM_vars( N_catl, cat_param%vegcls, cat_param%poros,    &
               mwRTM_param%poros, sfmc_l(:,n_e), tsurf_l(:,n_e), tp_l(1,:),    &
               smoist, stemp_l(:,n_e) )
          
          ! calculate brightness temperatures
          
          do j=1,N_TbuniqFreqAngRTMid
             
             freq       = Tb_freq_ang_RTMid(j,1)
             inc_angle  = Tb_freq_ang_RTMid(j,2)
             RTM_id     = Tb_freq_ang_RTMid(j,3)

             ! Select a specific configuration of the RTM via the field 
             ! "RTM_ID" in the "obs_param" type. 
             !
             ! %RTM_ID = ID of radiative transfer model to use for Tb forward modeling
             !           (subroutine get_obs_pred()) 
             !           0 = none
             !           1 = tau-omega model as in De Lannoy et al. 2013 (doi:10.1175/JHM-D-12-092.1)
             !           2 = same as 1 but without Pellarin atmospheric corrections
             !           3 = ...
             
             select case (RTM_id)
                
             case (1)
                
                ! bug fix: previously, mwRTM_get_Tb() was called without specifying the
                !          sub-array of "stemp_l" that corresponds to ensemble member n_e
                !          - reichle, 11 Dec 2013  
                
                call mwRTM_get_Tb(                                              &
                     N_catl, freq, inc_angle, mwRTM_param, tile_coord_l%elev,   &
                     lai, smoist, stemp_l(:,n_e), SWE, met_force%Tair, .true.,  &
                     'wang',Tb_h_vec, Tb_v_vec )
                
             case (2)
                
                call mwRTM_get_Tb(                                              &
                     N_catl, freq, inc_angle, mwRTM_param, tile_coord_l%elev,   &
                     lai, smoist, stemp_l(:,n_e), SWE, met_force%Tair, .false., &
                     'wang',Tb_h_vec, Tb_v_vec )

             case (3)
                
                call mwRTM_get_Tb(                                              &
                     N_catl, freq, inc_angle, mwRTM_param, tile_coord_l%elev,   &
                     lai, smoist, stemp_l(:,n_e), SWE, met_force%Tair, .false., &
                     'mironov',Tb_h_vec, Tb_v_vec )

                
             case default
                
                write (tmpstring10,*) RTM_ID
                
                err_msg = 'unknown or inconsistent RTM_ID=' // tmpstring10
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                
             end select

             Tb_h_l(:,j,n_e) = Tb_h_vec
             Tb_v_l(:,j,n_e) = Tb_v_vec
             
          end do
          
       end if

       ! ----------------------------------------------------------------
       !
       ! model-based QC in *l*ocal tile space:
       !
       ! - overwrite obs_pred fields w/ no-data-values as needed
       ! - whether QC is needed depends on get_*_lH (NOT get_*_l)!
       !
       ! perform model-based QC only before the EnKF update
       ! (not needed afterwards)
       
       if (beforeEnKFupdate) then

          if (get_sfmc_lH)  &
               call qc_model_based_for_sat_sfmc( N_catl, precip, SWE, tsurf_l(:,n_e), &
               sfmc_l(:,n_e) )
          
          if (get_tsurf_lH) &
               call qc_model_based_for_sat_tsurf(N_catl, precip, SWE, tp_l(1,:),      &
               tsurf_l(:,n_e) )
          
          if (get_Tb_lH) then
             
             do j=1,N_TbuniqFreqAngRTMid
                
                call qc_model_based_for_Tb( N_catl, precip, Tb_h_l(:,j,n_e) )
                call qc_model_based_for_Tb( N_catl, precip, Tb_v_l(:,j,n_e) )
                
             end do
             
          end if

       end if

    end do  ! loop through ens members
    
    ! ----------------------------------------------------------------
    !
    ! gather observed fields into local halo domain (for all processors)
    
    ! determine N_catlH and tile_coord_lH  

    N_fields = 0  ! set to zero temporarily, not yet needed
    ! move up the allocation. The input should be allocated in debug mode although it is not used
    ! allocate and assemble tile_data_l
    allocate(tile_data_l(0,0,0))  ! for debugging to pass  
    call get_tiles_in_halo( N_catl, N_fields, N_ens, tile_data_l, tile_coord_l,  &
         N_catf, tile_coord_f, N_catl_vec, low_ind, xhalo, yhalo,                &
         N_catlH, tile_coord_lH=tile_coord_lH )
    
    if (get_sfmc_lH)   allocate(sfmc_lH( N_catlH,                     N_ens))
    if (get_rzmc_lH)   allocate(rzmc_lH( N_catlH,                     N_ens))
    if (get_tsurf_lH)  allocate(tsurf_lH(N_catlH,                     N_ens))
    if (get_FT_lH)     allocate(FT_lH(   N_catlH,                     N_ens))
    if (get_Tb_lH)     allocate(stemp_lH(N_catlH,                     N_ens))
    if (get_Tb_lH)     allocate(Tb_h_lH( N_catlH,N_TbuniqFreqAngRTMid,N_ens))
    if (get_Tb_lH)     allocate(Tb_v_lH( N_catlH,N_TbuniqFreqAngRTMid,N_ens))
    
#ifdef LDAS_MPI
    
    ! count number of fields that need to be communicated (N_fields), allocate as needed
    
    call get_obs_pred_comm_helper( N_catl, N_ens, N_TbuniqFreqAngRTMid,          &
         get_sfmc_lH, get_rzmc_lH, get_tsurf_lH, get_FT_lH, get_Tb_lH, N_fields)
    
    ! allocate and assemble tile_data_l
    
    if (allocated(tile_data_l))  deallocate(tile_data_l)
    allocate(tile_data_l(N_catl,N_fields,N_ens))
    call get_obs_pred_comm_helper( N_catl, N_ens, N_TbuniqFreqAngRTMid,          &
         get_sfmc_lH, get_rzmc_lH, get_tsurf_lH, get_FT_lH, get_Tb_lH, N_fields, &
         option=1, tile_data=tile_data_l,                                        &
         sfmc=sfmc_l, rzmc=rzmc_l, tsurf=tsurf_l, FT=FT_l, stemp=stemp_l,        &
         Tb_h=Tb_h_l, Tb_v=Tb_v_l )
    
    ! communicate tile_data_l as needed and get tile_data_lH
    
    call get_tiles_in_halo( N_catl, N_fields, N_ens, tile_data_l, tile_coord_l,  &
         N_catf, tile_coord_f, N_catl_vec, low_ind, xhalo, yhalo,                &
         N_catlH, tile_data_lH=tile_data_lH )    
    
    ! read out sfmc, rzmc, etc. from tile_data_lH    
    
    call get_obs_pred_comm_helper( N_catlH, N_ens, N_TbuniqFreqAngRTMid,         &
         get_sfmc_lH, get_rzmc_lH, get_tsurf_lH, get_FT_lH, get_Tb_lH, N_fields, &
         option=2, tile_data=tile_data_lH,                                       &
         sfmc=sfmc_lH, rzmc=rzmc_lH, tsurf=tsurf_lH, FT=FT_l, stemp=stemp_lH,    &
         Tb_h=Tb_h_lH, Tb_v=Tb_v_lH )
    
    ! clean up
    
    if (associated(tile_data_lH))  deallocate(tile_data_lH)
    
#else             
    
    if (get_sfmc_lH)   sfmc_lH  = sfmc_l
    if (get_rzmc_lH)   rzmc_lH  = rzmc_l
    if (get_tsurf_lH)  tsurf_lH = tsurf_l
    if (get_FT_lH)     FT_lH    = FT_l
    if (get_Tb_lH)     stemp_lH = stemp_l
    if (get_Tb_lH)     Tb_h_lH  = Tb_h_l
    if (get_Tb_lH)     Tb_v_lH  = Tb_v_l
    
#endif
    if (allocated(tile_data_l))  deallocate(tile_data_l)
    ! ----------------------------------------------------------------
    !
    ! Get additional grid/tile information that is needed to map from tile
    ! to obs space

    if ( any(obs_param(1:N_obs_param)%FOV>FOV_threshold) )  then
       
       ! determine tile_grid_lH from tile_coord_lH
       
       tile_grid_lH = get_minExtent_grid( N_catlH, tile_coord_lH%pert_i_indg, tile_coord_lH%pert_j_indg,&
            tile_coord_lH%min_lon, tile_coord_lH%min_lat, tile_coord_lH%max_lon, tile_coord_lH%max_lat, &
            tile_grid_g) 
       
       allocate(N_tile_in_cell_ij_lH(tile_grid_lH%N_lon,tile_grid_lH%N_lat))
       
       ! first call: count how many tiles are in each tile_grid_lH cell

       call get_number_of_tiles_in_cell_ij( N_catlH,                                   &
            tile_coord_lH%pert_i_indg, tile_coord_lH%pert_j_indg,                      &
            tile_grid_lH, N_tile_in_cell_ij_lH )
       
       ! second call: find out which tiles are in each tile_grid_lH cell
       !              [tile numbers in "tile_num_in_cell_ij_lH" are relative
       !               to local halo ("lH") domain]
       
       call get_tile_num_in_cell_ij( N_catlH,                                          &
            tile_coord_lH%pert_i_indg, tile_coord_lH%pert_j_indg,                      &
            tile_grid_lH, maxval(N_tile_in_cell_ij_lH), tile_num_in_cell_ij_lH )
       
    end if
    
    ! -----------------------
    
    allocate(ind_tmp(    N_catlH))
    allocate(tmp_ndst2(  N_catlH))
    allocate(tmp_wFOV(   N_catlH))
    allocate(tmp_weights(N_catlH))
    allocate(tmp_data(   N_catlH))

    do i=1,N_obsl
       
       this_species       = Observations_l(i)%species
       
       this_Tbspecies     = ind_obsparam2Tbspecies(Observations_l(i)%species)
       
       if (this_Tbspecies>0) then
          
          this_TbuniqFreqAngRTMid = ind_Tbspecies2TbuniqFreqAngRTMid(this_Tbspecies)

       else
          
          this_TbuniqFreqAngRTMid = -999
          
       end if

       this_tilenum       = Observations_l(i)%tilenum  ! tilenum w.r.t. "full" domain
       
       this_lon           = Observations_l(i)%lon
       this_lat(1)        = Observations_l(i)%lat

       this_FOV           = obs_param(this_species)%FOV
       
       this_pol           = obs_param(this_species)%pol
       
       ! ----------------------------------------
       !
       ! map from full domain to this Observation

       use_distance_weights = .false.   ! initialize
       
       if (this_FOV < FOV_threshold) then
          
          ! equate obs footprint with nearest tile
          
          N_tmp      = 1
          
          ind_tmp(1) = f2l(this_tilenum) ! requires that "lH" starts w/ "local" tiles
          
       else
          
          ! find all tiles w/in given distance from obs, see
          !  LDASsa_DEFAULT_inputs_ensupd.nml for details!
          
          ! get appropriate distances in units of [deg] 
          
          if     ( trim(obs_param(this_species)%FOV_units)=='deg' ) then

             ! IMPORTANT: search distance is FOV when FOV_units='deg' !!!
             
             r_x(1) = this_FOV
             r_y    = this_FOV
             
          elseif ( trim(obs_param(this_species)%FOV_units)=='km'  ) then
             
             ! convert from [km] (FOV) to [deg] 
             
             ! IMPORTANT: search distance is fac_search_FOV_km*FOV when FOV_units='km' !!!
             
             call dist_km2deg( fac_search_FOV_km*this_FOV, 1, this_lat, r_x, r_y )
             
             use_distance_weights = .true.
             
          else
             
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown FOV_units (ii)')
             
          end if
          
          call get_tile_num_in_ellipse( this_lon, this_lat(1), r_x(1), r_y,  &
               N_catlH, tile_coord_lH, tile_grid_lH,                         &
               N_tile_in_cell_ij_lH, tile_num_in_cell_ij_lH,                 &
               N_tmp,  ind_tmp, tmp_ndst2                           )
          
          ! N_tmp could be zero (if ellipse straddles dateline)
          ! - reichle, 17 Apr 2017

       end if
       
       
       if (N_tmp>0) then            ! map from tiles to obs space
          
          ! compute weights based on tile area and (if applicable) based on distance from obs
          
          if (.not. use_distance_weights) then
             
             ! ignore tmp_weights from get_tile_num_in_ellipse(),
             ! weights based only on tile area
             
             tmp_wFOV(   1:N_tmp) = 1.
             
             tmp_weights(1:N_tmp) = tile_coord_lH(ind_tmp(1:N_tmp))%area
             
          else
             
             ! use distance weights along with tile area weights
             
             ! convert normalized square distance "ndst2" from get_tile_num_in_ellipse()
             !  into distance-based weights
             
             ! normalized distance from get_tile_num_in_ellipse() is relative to r_x and r_y,
             !  first scale back so that distance is w.r.t. FOV
             
             tmp_ndst2(  1:N_tmp) = (fac_search_FOV_km**2) * tmp_ndst2(1:N_tmp) 
             
             tmp_wFOV(   1:N_tmp) = exp( -0.5*tmp_ndst2(1:N_tmp) )
             
             ! further adjust weights based on tile area
             
             tmp_weights(1:N_tmp) = tmp_wFOV(1:N_tmp) * tile_coord_lH(ind_tmp(1:N_tmp))%area
             
          end if
          
          do n_e=1,N_ens
             
             ! -----------------------------------------
             !
             ! fill Obs_pred with observed field, aggregated
             ! from tiles as appropriate for this Observation
             
             select case (trim(obs_param(this_species)%varname))
                
             case ('sfmc')
                
                tmp_data(1:N_tmp)    = sfmc_lH(  ind_tmp(1:N_tmp), n_e )
                
             case ('rzmc') 
                
                tmp_data(1:N_tmp)    = rzmc_lH(  ind_tmp(1:N_tmp), n_e )
                
             case ('tsurf')   
                
                tmp_data(1:N_tmp)    = tsurf_lH( ind_tmp(1:N_tmp), n_e ) 
                
             case ('FT')   
                
                tmp_data(1:N_tmp)    = FT_lH(    ind_tmp(1:N_tmp), n_e ) 
                
             case('Tb')
                
                ! start with QC based on model *soil* temperature, motivated by RFI
                ! (requires model estimates *and* observation together;
                !  only performed before the EnKF update b/c not needed afterwards)
                
                if (beforeEnKFupdate) then
                   
                   tmp_data(1:N_tmp) = stemp_lH( ind_tmp(1:N_tmp), n_e ) 
                   
                   call tile2obs_helper(                                           &
                        N_tmp, tmp_weights(1:N_tmp), tmp_data(1:N_tmp), tmp_stemp)
                   
                   ! if Tb is too warm, RFI is likely
                   
                   tmpRFI = ((Observations_l(i)%obs-tmp_stemp) > Tbobs_minus_stemp_max)
                   
                else
                   
                   tmpRFI = .false.
                   
                end if
                
                ! for EASE grids *ONLY*: screen for non-land surfaces (e.g., lakes)
                ! - reichle, 28 Mar 2015
                
                if (index(tile_grid_g%gridtype, 'EASEv')  /=0) then
                   
                   ! ASSUMPTIONS: 
                   !  - at most one land tile per grid cell
                   !  - grid cells have the same area (or at least have nearly
                   !     identical areas in the surrounding region; that is, 
                   !     a regular lat/lon grid with just one tile per grid
                   !     cell would be ok if that property could be asserted here)
                   
                   ! compute FOV-weighted average "frac_cell"
                   
                   tmp_data(1:N_tmp) = tile_coord_lH(ind_tmp(1:N_tmp))%frac_cell
                   
                   call tile2obs_helper(                                               &
                        N_tmp, tmp_wFOV(1:N_tmp), tmp_data(1:N_tmp), tmp_fraccell)
                   
                   ! check whether there is too much non-land in FOV 
                   !  (typically water, but could be land-ice)
                   
                   tmpWater = ( 1. - tmp_fraccell > EASE_max_water_frac )
                   
                else
                   
                   tmpWater = .false.
                   
                end if
                
                
                ! compute Obs_pred
                
                if (tmpRFI .or. tmpWater) then
                   
                   ! apply model-based QC (suspected RFI, too much non-land in FOV)
                   
                   tmp_data(1:N_tmp) = nodata_generic  ! results in Obs_pred(i,n_e) = nodata
                   
                else
                   
                   select case (this_pol)
                      
                   case (1)   ! H-pol
                      
                      tmp_data(1:N_tmp) = Tb_h_lH( ind_tmp(1:N_tmp), this_TbuniqFreqAngRTMid, n_e) 
                      
                   case (2)   ! V-pol
                      
                      tmp_data(1:N_tmp) = Tb_v_lH( ind_tmp(1:N_tmp), this_TbuniqFreqAngRTMid, n_e) 
                      
                   case default
                      
                      err_msg = 'unknown obs_param%pol for varname=Tb'
                      call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                      
                   end select
                   
                end if
                
             case default
                
                err_msg = 'unknown obs_param%varname'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                
             end select
             
             ! average from tiles to obs
             
             call tile2obs_helper(N_tmp, tmp_weights(1:N_tmp), tmp_data(1:N_tmp), tmpreal)
             
             Obs_pred_l(i,n_e) = tmpreal
             
          end do     ! loop through ens members
          
       end if        ! N_tmp>0

       ! ----------------------------------------------------
       
       ! additional, model-based QC (only done before EnKF update): 
       !  potentially remove obs if innovation is too large

       ! IMPORTANT: 
       ! The calibration of mwRTM parameters relies on the "getinnov" capability 
       ! of LDASsa in conjunction with a (poor) prior guess of mwRTM_param.
       ! If model-based QC for Tb is added that in any way uses mwRTM_param, 
       ! observations might be discarded when they should not be during runs that 
       ! prepare for the mwRTM parameter calibration
       ! - reichle, 12 Dec 2013
       
       if (beforeEnKFupdate) then          
          
          if (trim(obs_param(this_species)%varname)=='tsurf') then 
             
             ! potentially eliminate obs (except if "bias_Npar>0" and "obsbias_ok==FALSE")
             
             if ( obs_param(this_species)%bias_Npar>0  .and.  (.not. obsbias_ok_tmp(i)) ) then
                
                ! do nothing (ie, keep obs), obs bias estimate is spinning up
                
             else
                
                ! compute ensemble mean innovation
                
                tmpreal = Observations_l(i)%obs - sum(Obs_pred_l(i,1:N_ens))/real(N_ens)
                
                ! check whether innovation is considered too large
                !
                ! rough estimate of innovations variance:  R + HPH^t ~ 2*R 
                !
                ! threshold for acceptance:  innovation**2 < (fac**2) * (2*R)  
                !
                ! changed threshold from fac=2 to fac=5 times innovations std-dev because
                ! cat_bias estimation might never get started otherwise;
                ! typical values of sqrt(obsvar) are between 1.3K (night) and 2.1K (day), 
                ! which corresponds to thresholds of between 9K (night) and 15K (day)
                ! for innovations values;  
                ! egregious outliers of Tskin retrievals that result from sensing cloud
                ! tops would presumably still be filtered out.
                ! - reichle, 30 Jun 2015

                if ( tmpreal**2 > (5.**2)*2.*Observations_l(i)%obsvar ) then
                   
                   Obs_pred_l(i,1:N_ens) = nodata_generic  ! eliminate obs
                   
                end if
                
             end if
             
          end if
          
       end if
       
    end do        ! loop through Observations
    
    ! ----------------------------------------------------------------
    ! 
    ! clean up
    
    if (associated(N_tile_in_cell_ij_lH))    deallocate(N_tile_in_cell_ij_lH)
    if (associated(tile_num_in_cell_ij_lH))  deallocate(tile_num_in_cell_ij_lH)
    
    if (allocated(ind_tmp))                  deallocate(ind_tmp)
    if (allocated(tmp_ndst2))                deallocate(tmp_ndst2)
    if (allocated(tmp_weights))              deallocate(tmp_weights)
    if (allocated(tmp_data))                 deallocate(tmp_data)
    
    if (associated(tile_coord_lH))           deallocate(tile_coord_lH)
    
    if (get_Tb_l)           deallocate(Tb_h_l)
    if (get_Tb_l)           deallocate(Tb_v_l)
    
    if (get_sfmc_lH)        deallocate(sfmc_lH) 
    if (get_rzmc_lH)        deallocate(rzmc_lH) 
    if (get_tsurf_lH)       deallocate(tsurf_lH)
    if (get_FT_lH)          deallocate(FT_lH)
    if (get_Tb_lH)          deallocate(stemp_lH)           
    if (get_Tb_lH)          deallocate(Tb_h_lH) 
    if (get_Tb_lH)          deallocate(Tb_v_lH) 
    
    ! ----------------------------------------------------------------
    ! 
    ! deal with no-data-values, compute ens mean and var of Obs_pred
    
    if (beforeEnKFupdate) then
       
       ! when used for "forecast" delete obs if Obs_pred is no-data-value
       
       j = 0
       
       do i=1,N_obsl
          
          if (all(abs(Obs_pred_l(i,1:N_ens)-nodata_generic)>nodata_tol_generic))  then
             
             ! keep this obs
             
             j = j + 1
             
             Observations_l(j        ) = Observations_l(i        )
             
             Obs_pred_l(    j,1:N_ens) = Obs_pred_l(    i,1:N_ens)
             
             ! fill in fcst and fcstvar
             
             if (N_ens>1) then
                
                call row_variance( 1, N_ens, Obs_pred_l(j,1:N_ens), tmpvar, tmpmean )
                
                ! inflate fcstvar
                
                if (inflation_factor > 0.)  tmpvar(1) = tmpvar(1) * inflation_factor**2
                             
             else
                
                tmpmean(1) = Obs_pred_l(j,1)
                
                tmpvar(1)  = nodata_generic
                
             end if
             
             Observations_l(j)%fcst    = tmpmean(1)
             Observations_l(j)%fcstvar = tmpvar(1)
             
          end if
          
       end do
       
       N_obsl = j
       
    else
       
       ! keep *all* obs even if Obs_pred turns out to be nodata
       !
       ! Obs_pred can still be no-data (e.g., if model soil temperature is too cold,
       !  mwRTM_get_Tb() returns a no-data-value)

       do i=1,N_obsl
          
          ! fill in ana and anavar
          
          if (N_ens>1) then
                          
             if (any(abs(Obs_pred_l(i,1:N_ens)-nodata_generic)<nodata_tol_generic )) then
                
                tmpmean(1) = nodata_generic
                tmpvar(1)  = nodata_generic
                
             else
                
                call row_variance( 1, N_ens, Obs_pred_l(i,1:N_ens), tmpvar, tmpmean )
                
                ! no need to inflate analysis Obs_pred because state increments already included
                !  impact of inflation

             end if
             
          else
             
             tmpmean(1) = Obs_pred_l(i,1)
             
             tmpvar(1)  = nodata_generic
             
          end if
          
          Observations_l(i)%ana    = tmpmean(1)
          Observations_l(i)%anavar = tmpvar(1)
          
       end do

    end if

    ! ----------------------------------------------------------------
    
  end subroutine get_obs_pred

  ! *****************************************************************

  subroutine get_obs_pred_comm_helper(                                           &
       N_cat, N_ens, N_Tb, get_sfmc, get_rzmc, get_tsurf, get_FT, get_Tb,        &
       N_fields, option, tile_data, sfmc, rzmc, tsurf, FT, stemp, Tb_h, Tb_v )
    
    ! bundle/unbundle individual fields into/from single array for more 
    ! efficient communication across processors
    !
    ! option:
    !   1         = assemble tile_data from sfmc, rzmc, etc. 
    !   2         = read sfmc, rzmc, etc. out from tile_data
    !   otherwise = count N_fields only
    !
    ! reichle, 19 June 2012
    
    ! ----------------------------------------------------------------------

    implicit none

    integer, intent(in)    :: N_cat, N_ens, N_Tb
    
    logical, intent(in)    :: get_sfmc, get_rzmc, get_tsurf, get_FT, get_Tb
        
    integer, intent(inout) :: N_fields
    
    integer,                               intent(in),    optional :: option 
    
    real, dimension(N_cat,N_fields,N_ens), intent(inout), optional :: tile_data
    
    real, dimension(N_cat,         N_ens), intent(inout), optional :: sfmc, rzmc
    real, dimension(N_cat,         N_ens), intent(inout), optional :: tsurf, FT, stemp
    real, dimension(N_cat,N_Tb,    N_ens), intent(inout), optional :: Tb_h, Tb_v
    
    ! -----------------------------------

    ! local variables
    
    integer :: k, ks, opt

    character(len=*), parameter :: Iam = 'get_obs_pred_comm_helper'
    character(len=400) :: err_msg
    
    ! -------------------------------------------------------------------------
    
    ! check optional inputs
    
    if (present(option)) then
       
       opt = option
       
    else
       
       opt = 0

    end if
    
    if ((opt==1) .or. (opt==2)) then
       
       if ( ((get_sfmc ) .and. (.not. present(sfmc )))    .or.            &
            ((get_rzmc ) .and. (.not. present(rzmc )))    .or.            &
            ((get_tsurf) .and. (.not. present(tsurf)))    .or.            &
            ((get_FT)    .and. (.not. present(FT   )))    .or.            &
            ((get_Tb)    .and. (.not. present(stemp)))    .or.            &
            ((get_Tb)    .and. (.not. present(Tb_h )))    .or.            &
            ((get_Tb)    .and. (.not. present(Tb_v )))             ) then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error 1')
       end if
       
       if ( (get_sfmc .or. get_rzmc .or. get_tsurf .or. get_FT .or. get_Tb) .and.   &
            (.not. present(tile_data))                                              &
            )  then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error 2')
       end if
       
    end if
    
    ! --------------------------
    
    k = 0
    
    if (get_sfmc)   then       
       
       k = k+1
       
       if (opt==1)  tile_data(1:N_cat,k,1:N_ens) = sfmc
       
       if (opt==2)  sfmc = tile_data(1:N_cat,k,1:N_ens)
       
    end if
    
    if (get_rzmc)   then
       
       k = k+1
              
       if (opt==1)  tile_data(1:N_cat,k,1:N_ens) = rzmc
       
       if (opt==2)  rzmc = tile_data(1:N_cat,k,1:N_ens)
                     
    end if
    
    if (get_tsurf)  then
       
       k = k+1
       
       if (opt==1)  tile_data(1:N_cat,k,1:N_ens) = tsurf
       
       if (opt==2)  tsurf = tile_data(1:N_cat,k,1:N_ens)
                            
    end if

    if (get_FT)     then
       
       k = k+1
       
       if (opt==1)  tile_data(1:N_cat,k,1:N_ens) = FT
       
       if (opt==2)  FT    = tile_data(1:N_cat,k,1:N_ens)
                            
    end if
        
    if (get_Tb)     then

       k = k+1
       
       if (opt==1)  tile_data(1:N_cat,k,1:N_ens) = stemp
       
       if (opt==2)  stemp = tile_data(1:N_cat,k,1:N_ens)
       
    end if
    
    if (get_Tb)     then
       
       ks = k+1     ! "k_start"
       k  = k+N_Tb  ! "k_end"

       if (opt==1)  tile_data(1:N_cat,ks:k,1:N_ens) = Tb_h
       
       if (opt==2)  Tb_h = tile_data(1:N_cat,ks:k,1:N_ens)
       
    end if
    
    if (get_Tb)     then

       ks = k+1     ! "k_start"
       k  = k+N_Tb  ! "k_end"
              
       if (opt==1)  tile_data(1:N_cat,ks:k,1:N_ens) = Tb_v
       
       if (opt==2)  Tb_v = tile_data(1:N_cat,ks:k,1:N_ens)
              
    end if
    
    N_fields = k
    
  end subroutine get_obs_pred_comm_helper

  ! *****************************************************************

  subroutine tile2obs_helper( N_tile, tile_weights, tile_data, obs_data )
    
    ! reichle, 14 Jun 2011
    !
    ! map from tile space to single obs 
    !
    ! inputs:
    !
    ! N_tile       : number of tiles that contribute to obs
    ! tile_data    : must be zoomed into the N_tile tiles to be averaged into obs_data
    ! tile_weights : weights (normalization done in this subroutine)
    
    implicit none
    
    integer,                    intent(in)  :: N_tile
    
    real,    dimension(N_tile), intent(in)  :: tile_weights, tile_data
    
    real,                       intent(out) :: obs_data
    
    ! local variables
    
    real,    dimension(N_tile)  :: weights   ! normalized weights
    
    ! ----------------------------------------
    !
    ! return no-data-value if any of the inputs is a no-data-value

    if (any(abs(tile_data-nodata_generic)<nodata_tol_generic)) then
       
       obs_data = nodata_generic
                  
    else
       
       weights  = tile_weights/sum(tile_weights)   ! normalized weights
       
       obs_data = sum( weights*tile_data )
       
    end if
    
  end subroutine tile2obs_helper
  
  ! *****************************************************************

  subroutine qc_model_based_for_sat_sfmc( N_cat, precip, SWE, tsurf, &
       sfmc )
       
    ! simple model-based quality control for satellite sfmc observations
    !
    ! set sfmc "Obs_pred" to no-data when model indicates difficult
    ! retrieval conditions
    !
    ! reichle, 22 Sep 2006
    ! reichle, 14 Jun 2011 - moved from clsm_ensupd_read_obs.F90 and edited
    ! reichle, 23 Nov 2011 - changed tsurf_threshold and precip_threshold 
    !                        b/c QC now done for individual ensemble members
    !                        rather than the ensemble mean
    !
    ! --------------------------------------------------------------
    
    implicit none
    
    integer,                   intent(in)    :: N_cat
    
    real,    dimension(N_cat), intent(in)    :: precip   ! Rainf+Snowf [kg/m2/s]
    real,    dimension(N_cat), intent(in)    :: SWE      ! total SWE   [kg/m2]  
    real,    dimension(N_cat), intent(in)    :: tsurf    !             [K]
    
    real,    dimension(N_cat), intent(inout) :: sfmc
    
    ! local variables
    
    !real, parameter :: precip_threshold = 5./86400.                   ! [kg/m2/s]
    real, parameter :: precip_threshold = 10./86400.                  ! [kg/m2/s]
     
    real, parameter :: SWE_threshold    = 1.e-4                       ! [kg/m2]
    
    !real, parameter :: tsurf_threshold  = MAPL_TICE + 1.              ! [K]
    real, parameter :: tsurf_threshold  = MAPL_TICE + 0.2             ! [K]
        
    integer :: i
    
    ! ---------------------------------------    

    do i=1,N_cat
       
       ! delete obs 
       ! - if there is snow on the ground 
       ! - if it is raining/snowing
       ! - if surface temperature is around or below freezing

       if ( (precip(i)  > precip_threshold)        .or.     &
            (SWE(i)     > SWE_threshold)           .or.     &
            (tsurf(i)   < tsurf_threshold)               )  &
            sfmc(i) = nodata_generic
       
    end do
    
  end subroutine qc_model_based_for_sat_sfmc
  
  
  ! *****************************************************************

  subroutine qc_model_based_for_sat_tsurf( N_cat, precip, SWE, tp1,    &
       tsurf )
    
    ! simple model-based quality control for satellite tsurf observations
    !
    ! set tsurf "Obs_pred" to no-data when model indicates difficult
    ! retrieval conditions
    !
    ! reichle,  2 Nov 2004
    ! reichle, 14 Jun 2011 - moved from clsm_ensupd_read_obs.F90 and edited
    ! reichle, 23 Nov 2011 - changed precip_threshold b/c QC now done for individual 
    !                        ensemble members rather than the ensemble mean
    ! reichle, 14 Feb 2013 - added hard-coded option to eliminate frozen conditions in QC
    !                         based on tsurf and top layer soil temperature (tp1) 
    ! reichle, 26 Mar 2014 - set "avoid_frozen=.false." following advice from Clara
    !
    ! --------------------------------------------------------------
    
    implicit none
    
    integer,                   intent(in)    :: N_cat
    
    real,    dimension(N_cat), intent(in)    :: precip   ! Rainf+Snowf      [kg/m2/s]
    real,    dimension(N_cat), intent(in)    :: SWE      ! total SWE        [kg/m2]  
    real,    dimension(N_cat), intent(in)    :: tp1      ! soil temperature [C]
    
    real,    dimension(N_cat), intent(inout) :: tsurf    
    
    ! local variables
    
    real,    parameter :: precip_threshold = 10./86400.                  ! [kg/m2/s]
    
    real,    parameter :: SWE_threshold    = 1.e-4                       ! [kg/m2]
    
    logical, parameter :: avoid_frozen     = .false.
    
    real,    parameter :: temp_threshold   = MAPL_TICE + 2.              ! [K]
    
    integer :: i
    
    logical, dimension(N_cat) :: frozen
    
    real,    dimension(N_cat) :: tp1_in_Kelvin

    ! ---------------------------------------    
    
    frozen = .false.
    
    if (avoid_frozen) then
       
       tp1_in_Kelvin = tp1 + MAPL_TICE

       do i=1,N_cat
          
          if ( (tsurf(          i) < temp_threshold)         .or.     &
               (tp1_in_Kelvin(  i) < temp_threshold)                  &
               )                                                      &
               frozen(i) = .true. 
          
       end do
       
    end if
    
    
    do i=1,N_cat
       
       ! delete obs 
       ! - if there is snow on the ground 
       ! - if it is raining/snowing
       ! - if "avoid_frozen" and frozen
       
       if ( (precip(i)  > precip_threshold)        .or.     &
            (SWE(i)     > SWE_threshold)           .or.     &
            (frozen(i))                                     &
            )                                               &
            tsurf(i) = nodata_generic
       
    end do
    
  end subroutine qc_model_based_for_sat_tsurf
  
  
  ! *****************************************************************
  
  subroutine qc_model_based_for_Tb( N_cat, precip, Tb )
    
    ! simple model-based quality control for Tb observations
    !
    ! set Tb "Obs_pred" to no-data when model indicates difficult conditions
    !
    !
    ! GDL,     15 Nov 2010
    ! reichle, 27 May 2011 - included in LDASsa
    ! reichle, 14 Jun 2011 - moved from clsm_ensupd_read_obs.F90 and edited
    ! reichle, 23 Nov 2011 - changed precip_threshold b/c QC now done for individual 
    !                        ensemble members rather than the ensemble mean
    !
    ! --------------------------------------------------------------
    
    implicit none
    
    integer,                   intent(in)    :: N_cat
    
    real,    dimension(N_cat), intent(in)    :: precip   ! Rainf+Snowf      [kg/m2/s]
    
    real,    dimension(N_cat), intent(inout) :: Tb       ! brightness temp. [K]    
    
    ! local variables
    
    ! relatively large threshold for precip indirectly screens for standing water
    
    !real, parameter          :: precip_threshold = 25./86400.        ! [kg/m2/s]
    real, parameter          :: precip_threshold = 50./86400.        ! [kg/m2/s]
    
    integer :: i
    
    ! ---------------------------------------    

    do i=1,N_cat
       
       ! delete obs 
       ! - if there is heavy rain or snow
       
       ! NOTE: subroutine mwRTM_get_Tb already returns no-data-values
       ! - if there is snow on the ground 
       ! - if surface temperature is around or below freezing
       
       if ( (precip(i)   > precip_threshold)         )  &
            Tb(i) = nodata_generic
       
       ! IMPORTANT: 
       ! The calibration of mwRTM parameters relies on the "getinnov" capability 
       ! of LDASsa in conjunction with a (poor) prior guess of mwRTM_param.
       ! If model-based QC for Tb is added that in any way uses mwRTM_param, 
       ! observations might be discarded when they should not be during runs that 
       ! prepare for the mwRTM parameter calibration
       ! - reichle, 12 Dec 2013

    end do
    
  end subroutine qc_model_based_for_Tb
  
  ! *********************************************************************

  subroutine get_halo_obs( N_ens, N_catl, N_obsl, Observations_l, Obs_pred_l,  &
       tile_coord_l, xcompact, ycompact,                                       &
       N_obslH, Observations_lH, Obs_pred_lH )
    
    ! collect observations from other local domains (processors) that are 
    ! within the halo of the current local domain (processor)
    !
    ! ONLY collect obs with flag assim==.true.
    !
    ! Current implementation:
    !
    ! 1. Determine min/max lat/lon for locally managed Observations_l and
    !    min/max lat/lon for local halo
    ! 2. MPI *all*gather this information to all processors
    ! 3. Determine which processors need to communicate
    ! 4. Pairwise MPI_SENDRECV as needed
    !
    !
    ! A shorter and simpler but somewhat wasteful implementation (in terms of memory 
    ! and communications) could be as follows:
    !
    ! 1. MPI *all*gather local Observations_l to all processors
    ! 2. Determine index of (full domain) obs that are in halo of given local domain
    ! 3. MPI *all*gather local Obs_pred and extract as needed
    !
    ! 
    ! reichle,  2 Aug 2011
    ! reichle, 30 Sep 2011
    !
    ! ----------------------------------------------------------------------    
    
    implicit none
    
    integer,                                        intent(in)  :: N_ens, N_catl, N_obsl
    
    type(obs_type),        dimension(N_obsl),       intent(in)  :: Observations_l
    
    real,                  dimension(N_obsl,N_ens), intent(in)  :: Obs_pred_l
    
    type(tile_coord_type), dimension(:),            pointer     :: tile_coord_l    ! in
    
    real,                                           intent(in)  :: xcompact, ycompact

    integer,                                        intent(out) :: N_obslH
    
    type(obs_type),        dimension(:),            pointer     :: Observations_lH ! out
     
    real,                  dimension(:,:),          pointer     :: Obs_pred_lH     ! out
    
    ! local variables
    
    integer :: i, j, m, n, N_obslH_max

    integer :: N_recv, N_recvi, N_sendi, tag1, tag2

    real    :: lon, lat

    real    :: obsl_minlon, obsl_maxlon, obsl_minlat, obsl_maxlat
    real    :: halo_minlon, halo_maxlon, halo_minlat, halo_maxlat

    real    :: minlon, maxlon, minlat, maxlat
    real    :: ll_lon, ur_lon, ll_lat, ur_lat

    integer, dimension(numprocs) :: N_obsl_vec

    real,    dimension(numprocs) :: obsl_minlon_vec, obsl_maxlon_vec
    real,    dimension(numprocs) :: obsl_minlat_vec, obsl_maxlat_vec

    real,    dimension(numprocs) :: halo_minlon_vec, halo_maxlon_vec
    real,    dimension(numprocs) :: halo_minlat_vec, halo_maxlat_vec

    logical, dimension(numprocs,numprocs) :: need_obsl
    
    type(obs_type), dimension(:),   pointer     :: Observations_l_recv => null()
    
    real,           dimension(:,:), allocatable :: Obs_pred_l_recv
        
    ! -------------------------------------------------
    !
    ! determine and communicate min/max lat/lon of Observations_l 
    !  managed by given processor
    !
    ! WARNING: this will most likely create problems if a given local domain
    !  crosses the dateline!
    !
    !
    ! NOTE: All Observations lat/lon values must be "good" 
    !       (ie, *not* no-data-values)
    
    if (N_obsl>0) then
       
       ! make sure the rectangle does not have zero area
       ! (which might happen if there is only one obs)
       
       obsl_minlon = minval(Observations_l%lon) - 0.01
       obsl_maxlon = maxval(Observations_l%lon) + 0.01
       obsl_minlat = minval(Observations_l%lat) - 0.01
       obsl_maxlat = maxval(Observations_l%lat) + 0.01
       
    else
       
       obsl_minlon = nodata_generic
       obsl_maxlon = nodata_generic
       obsl_minlat = nodata_generic
       obsl_maxlat = nodata_generic
       
    end if
    
#ifdef LDAS_MPI
    
    call MPI_AllGather(                          &
         N_obsl,                 1, MPI_integer, &
         N_obsl_vec,             1, MPI_integer, & 
         mpicomm, mpierr ) 
    
    call MPI_AllGather(                          &
         obsl_minlon,            1, MPI_real,    &
         obsl_minlon_vec,        1, MPI_real,    & 
         mpicomm, mpierr ) 

    call MPI_AllGather(                          &
         obsl_maxlon,            1, MPI_real,    &
         obsl_maxlon_vec,        1, MPI_real,    & 
         mpicomm, mpierr ) 

    call MPI_AllGather(                          &
         obsl_minlat,            1, MPI_real,    &
         obsl_minlat_vec,        1, MPI_real,    & 
         mpicomm, mpierr ) 

    call MPI_AllGather(                          &
         obsl_maxlat,            1, MPI_real,    &
         obsl_maxlat_vec,        1, MPI_real,    & 
         mpicomm, mpierr ) 

#else
    
    N_obsl_vec(1)      = N_obsl

    obsl_minlon_vec(1) = obsl_minlon
    obsl_maxlon_vec(1) = obsl_maxlon
    obsl_minlat_vec(1) = obsl_minlat
    obsl_maxlat_vec(1) = obsl_maxlat
    
#endif
        
    ! ------------------------------------------------------------
    !
    ! determine and communicate min/max lat/lon of halo
    
    halo_minlon = minval(tile_coord_l%com_lon) - 1.25*xcompact
    halo_maxlon = maxval(tile_coord_l%com_lon) + 1.25*xcompact
    halo_minlat = minval(tile_coord_l%com_lat) - 1.25*ycompact
    halo_maxlat = maxval(tile_coord_l%com_lat) + 1.25*ycompact

    ! simple approach to dateline issue (cut halo back to at most -180:180, -90:90)
    ! - reichle, 28 May 2013
    
    halo_minlon = max(halo_minlon,-180.)
    halo_maxlon = min(halo_maxlon, 180.)
    halo_minlat = max(halo_minlat, -90.)
    halo_maxlat = min(halo_maxlat,  90.)
        

#ifdef LDAS_MPI
    
    call MPI_AllGather(                       &
         halo_minlon,            1, MPI_real, &
         halo_minlon_vec,        1, MPI_real, & 
         mpicomm, mpierr ) 

    call MPI_AllGather(                       &
         halo_maxlon,            1, MPI_real, &
         halo_maxlon_vec,        1, MPI_real, & 
         mpicomm, mpierr ) 

    call MPI_AllGather(                       &
         halo_minlat,            1, MPI_real, &
         halo_minlat_vec,        1, MPI_real, & 
         mpicomm, mpierr ) 

    call MPI_AllGather(                       &
         halo_maxlat,            1, MPI_real, &
         halo_maxlat_vec,        1, MPI_real, & 
         mpicomm, mpierr ) 

#else

    halo_minlon_vec(1) = halo_minlon
    halo_maxlon_vec(1) = halo_maxlon
    halo_minlat_vec(1) = halo_minlat
    halo_maxlat_vec(1) = halo_maxlat
    
#endif
    
    ! ------------------------------------------------------------
    !
    ! determine which processors need to communicate (directional!)
    !
    ! need_obsl(i,j) = .true. 
    !    ==> processor i needs Observations_l from processor j
    
    need_obsl = .false.
    
    do i=1,numprocs            ! i = "tile-space" 
       
       ll_lon = halo_minlon_vec(i)
       ur_lon = halo_maxlon_vec(i)
       ll_lat = halo_minlat_vec(i)
       ur_lat = halo_maxlat_vec(i)

       do j=1,numprocs         ! j = "obs-space"
          
          if ( (N_obsl_vec(j)>0) .and. (i/=j) ) then
             
             minlon = obsl_minlon_vec(j)
             maxlon = obsl_maxlon_vec(j)
             minlat = obsl_minlat_vec(j)
             maxlat = obsl_maxlat_vec(j)
             
             ! processor i needs Observations_l from processor j 
             ! if bounding box around Observations_l from j overlaps
             ! with bounding box plus halo of processor i
             
             if ( (min(ur_lon,maxlon) - max(ll_lon,minlon))>0.   .and. &
                  (min(ur_lat,maxlat) - max(ll_lat,minlat))>0. )       &
                  need_obsl(i,j) = .true.
             
          else
             
             need_obsl(i,i) = .true.
             
          end if
          
       end do
    end do

    ! determine maximum number of obs within halo (incl locally managed obs)
    
    N_obslH_max = 0            ! note: need_obsl(j,j)=.true.
    
    do j=1,numprocs          
       
       if (need_obsl(myid+1,j))  N_obslH_max = N_obslH_max + N_obsl_vec(j)
       
    end do
    
    ! allocate Observations_lH, Obs_pred_lH

    allocate(Observations_lH(N_obslH_max      ))
    
    allocate(Obs_pred_lH(    N_obslH_max,N_ens))

    ! ------------------------------------------------------------
    !
    ! initialize Observations_lH, Obs_pred_lH with local Observations
    ! (use only those with flag "assim"==.true.)
    
    m = 0
    
    do n=1,N_obsl
       
       if (Observations_l(n)%assim) then
          
          m = m+1
          
          Observations_lH(m  ) = Observations_l(n)
          
          Obs_pred_lH(    m,:) = Obs_pred_l(n,:)
          
       end if
       
    end do

    N_obslH = m

#ifdef LDAS_MPI
    
    ! pairwise communication
    !
    ! use MPI_SENDRECV to avoid "hung" or "deadlocked" communications
    
    do i=1,numprocs
       
       do j=i+1,numprocs   ! loop only through upper triangle of (i,j) matrix
          
          ! determine how many Observations should be 
          !   - received by i (same as sent by j), and 
          !   - sent by i (same as received by j)

          if (need_obsl(i,j)) then

             N_recvi = N_obsl_vec(j)

          else
             
             N_recvi = 0

          end if
          
          if (need_obsl(j,i)) then
             
             N_sendi = N_obsl_vec(i)
             
          else
             
             N_sendi = 0
   
          end if

          tag1 = 1
          tag2 = 2

          if (N_recvi>0 .or. N_sendi>0) then
             
             if     (myid+1==i) then
                
                ! allocate Observations_l_recv, Obs_pred_recv
                
                allocate(Observations_l_recv(N_recvi      ))
                allocate(Obs_pred_l_recv(    N_recvi,N_ens))   
   
                ! obtain Observations_l from processor j
                ! send Observations_l to processor j
                
                call MPI_SENDRECV( &
                     Observations_l,      N_sendi, MPI_obs_type, j-1, tag1, &
                     Observations_l_recv, N_recvi, MPI_obs_type, j-1, tag2, &
                     mpicomm, mpistatus, mpierr ) 

                call MPI_SENDRECV( &
                     Obs_pred_l,      N_sendi*N_ens, MPI_real, j-1, tag1, &
                     Obs_pred_l_recv, N_recvi*N_ens, MPI_real, j-1, tag2, &
                     mpicomm, mpistatus, mpierr ) 

             elseif (myid+1==j) then
                
                ! allocate Observations_l_recv, Obs_pred_recv
                
                allocate(Observations_l_recv(N_sendi      ))
                allocate(Obs_pred_l_recv(    N_sendi,N_ens))   
                
                call MPI_SENDRECV( &
                     Observations_l,      N_recvi, MPI_obs_type, i-1, tag2, &
                     Observations_l_recv, N_sendi, MPI_obs_type, i-1, tag1, &
                     mpicomm, mpistatus, mpierr ) 
                
                call MPI_SENDRECV( &
                     Obs_pred_l,      N_recvi*N_ens, MPI_real, i-1, tag2, &
                     Obs_pred_l_recv, N_sendi*N_ens, MPI_real, i-1, tag1, &
                     mpicomm, mpistatus, mpierr ) 

             end if

             ! put received obs into Observations_lH

             if (myid+1==i .or. myid+1==j) then
                
                if (myid+1==i) then
                   
                   N_recv = N_recvi

                   ll_lon = halo_minlon_vec(i)
                   ur_lon = halo_maxlon_vec(i)
                   ll_lat = halo_minlat_vec(i)
                   ur_lat = halo_maxlat_vec(i)
                   
                else   ! my_id+1==j
                   
                   N_recv = N_sendi
                   
                   ll_lon = halo_minlon_vec(j)
                   ur_lon = halo_maxlon_vec(j)
                   ll_lat = halo_minlat_vec(j)
                   ur_lat = halo_maxlat_vec(j)
                   
                end if

                m = N_obslH
                
                do n=1,N_recv

                   lon = Observations_l_recv(n)%lon
                   lat = Observations_l_recv(n)%lat
                   
                   if ( Observations_l_recv(n)%assim                        .and.       &
                        is_in_rectangle(lon,lat,ll_lon,ll_lat,ur_lon,ur_lat)     ) then
                      
                      m = m+1
                      
                      Observations_lH(m  ) = Observations_l_recv(n)
                      
                      Obs_pred_lH(    m,:) = Obs_pred_l_recv(    n,:)
                      
                   end if
                   
                end do
                
                N_obslH = m

             end if
             
          end if
          
          call MPI_BARRIER( mpicomm, mpierr )

          if (associated(Observations_l_recv))  deallocate(Observations_l_recv)

          if (allocated(Obs_pred_l_recv))       deallocate(Obs_pred_l_recv)
          
       end do
    end do

#else
    
    ! previous block kicks in only for numprocs>1 (not for OpenMP and sequential)
    
#endif
    
  end subroutine get_halo_obs
  
  ! *********************************************************************
  
  subroutine get_tiles_in_halo( N_catl, N_fields, N_ens, tile_data_l, tile_coord_l,  &
       N_catf, tile_coord_f, N_catl_vec, low_ind, xhalo, yhalo,                      &
       N_catlH, tile_coord_lH, tile_data_lH )
    
    ! collect (bundled) tile_data from other local domains (processors) that are 
    ! within the halo of the current local domain (processor) for the purpose
    ! of calculating Obs_pred
    !
    ! Current implementation:
    !
    ! 1. Determine min/max lat/lon for all local domains
    ! 2. Determine which processors need to communicate
    ! 3. Pairwise MPI_SENDRECV as needed
    !
    ! reichle, 19 Jun 2012
    ! reichle, 27 Mar 2014 - renamed to better distinguish from "get_halo_around_tile()"
    !
    ! ----------------------------------------------------------------------    
    
    implicit none
    
    integer,                                    intent(in)  :: N_catl, N_fields
    integer,                                    intent(in)  :: N_ens,  N_catf
    
    real,    dimension(N_catl,N_fields,N_ens),  intent(in)  :: tile_data_l
    
    type(tile_coord_type), dimension(:),        pointer     :: tile_coord_l  ! in
    type(tile_coord_type), dimension(:),        pointer     :: tile_coord_f  ! in
    
    integer,               dimension(numprocs), intent(in)  :: N_catl_vec
    integer,               dimension(numprocs), intent(in)  :: low_ind
    
    real,                  dimension(numprocs), intent(in)  :: xhalo, yhalo
    
    integer,                                    intent(out) :: N_catlH
    
    real,                  dimension(:,:,:), pointer, optional :: tile_data_lH  ! out
    
    type(tile_coord_type), dimension(:),     pointer, optional :: tile_coord_lH ! out
    
    ! local variables
    
    integer :: i, j, istart, iend, istart_f, iend_f
    
    integer :: N_recv, N_recvi, N_sendi, tag1, tag2
    
    real    :: minlon, maxlon, minlat, maxlat
    real    :: ll_lon, ur_lon, ll_lat, ur_lat
    
    real,    dimension(numprocs)           :: catl_minlon_vec, catl_maxlon_vec
    real,    dimension(numprocs)           :: catl_minlat_vec, catl_maxlat_vec
    
    real,    dimension(numprocs)           :: halo_minlon_vec, halo_maxlon_vec
    real,    dimension(numprocs)           :: halo_minlat_vec, halo_maxlat_vec
    
    logical, dimension(numprocs,numprocs)  :: need_catl
    
    real,    dimension(:,:,:), allocatable :: tile_data_l_recv
    
    ! -------------------------------------------------
    !
    ! each processor determines min/max lat/lon of tiles managed by all 
    ! processors (without and with halo) 
    !
    ! WARNING: this will most likely create problems if a given local domain
    !  crosses the dateline!
    
    do i=1,numprocs
      
       if (N_catl_vec(i) <= 0) cycle    ! nothing to do for this processor
 
       istart = low_ind(i)
       iend   = istart + N_catl_vec(i) - 1

       ! use center-of-mass of tiles (rather than min/max_lon, min/max_lat)
       
       catl_minlon_vec(i) = minval( tile_coord_f(istart:iend)%com_lon )
       catl_maxlon_vec(i) = maxval( tile_coord_f(istart:iend)%com_lon )
       catl_minlat_vec(i) = minval( tile_coord_f(istart:iend)%com_lat )
       catl_maxlat_vec(i) = maxval( tile_coord_f(istart:iend)%com_lat )
       
       halo_minlon_vec(i) = catl_minlon_vec(i) - xhalo(i)
       halo_maxlon_vec(i) = catl_maxlon_vec(i) + xhalo(i)
       halo_minlat_vec(i) = catl_minlat_vec(i) - yhalo(i)
       halo_maxlat_vec(i) = catl_maxlat_vec(i) + yhalo(i)
              
    end do
    
    ! ------------------------------------------------------------
    !
    ! determine which processors need to communicate (symmetric!)
    !
    ! need_catl(i,j) = .true. 
    !    ==> processors i and j need to exchange tile_data_l
    
    need_catl = .false.

    do i=1,numprocs            

       if (N_catl_vec(i) >0) need_catl(i,i) = .true.

    end do
    
    do i=1,numprocs            

       ! all tiles within the following rectangle are needed by proc i

       if ( N_catl_vec(i) <= 0) cycle    ! nothing to do for this processor

       ll_lon = halo_minlon_vec(i)
       ur_lon = halo_maxlon_vec(i)
       ll_lat = halo_minlat_vec(i)
       ur_lat = halo_maxlat_vec(i)
       
       do j=i+1,numprocs
     
          if (N_catl_vec(j) <= 0) cycle  ! nothing to do for this processor
             
          minlon = catl_minlon_vec(j)
          maxlon = catl_maxlon_vec(j)
          minlat = catl_minlat_vec(j)
          maxlat = catl_maxlat_vec(j)
          
          ! processor i needs tile_data_l from processor j 
          ! if bounding box around tile_data_l(j) overlaps
          ! with bounding box plus halo of processor i
          
          if ( (min(ur_lon,maxlon) - max(ll_lon,minlon))>0.   .and. &
               (min(ur_lat,maxlat) - max(ll_lat,minlat))>0. )       &
               need_catl(i,j) = .true.
          
          need_catl(j,i) = need_catl(i,j)
       end do
    end do
    
    ! determine number of tiles within halo (incl local tiles)
    
    N_catlH = 0            ! note: need_catl(j,j)=.true.
    
    do j=1,numprocs          
       
       if (need_catl(myid+1,j))  N_catlH = N_catlH + N_catl_vec(j)
       
    end do
    
    ! allocate tile_coord_lH, tile_data_lH
    
    if (present(tile_coord_lH))  allocate(tile_coord_lH(N_catlH               ))
    if (present(tile_data_lH))   allocate(tile_data_lH( N_catlH,N_fields,N_ens))
    
    ! ------------------------------------------------------------
    !
    ! initialize tile_coord_lH, tile_data_lH with local tile_coord_l, 
    ! tile_data_l
    
    N_catlH = N_catl ! will grow as data from other processors are appended
    
    if (present(tile_coord_lH)) tile_coord_lH(1:N_catlH) = tile_coord_l(1:N_catlH)
    
    if (present(tile_data_lH))  tile_data_lH(1:N_catlH,1:N_fields,1:N_ens) = tile_data_l
    
#ifdef LDAS_MPI
    
    ! pairwise communication (symmetric!)
    !
    ! use MPI_SENDRECV to avoid "hung" or "deadlocked" communications
    
    do i=1,numprocs
       
       do j=i+1,numprocs   ! loop only through upper triangle of (i,j) matrix
          
          ! determine how many elements should be 
          !   - received by i (same as sent by j), and 
          !   - sent by i (same as received by j)
          
          if (need_catl(i,j)) then
             
             N_recvi = N_catl_vec(j)
             N_sendi = N_catl_vec(i)
                          
          else
             
             N_recvi = 0
             N_sendi = 0

          end if
          
          tag1 = 1
          tag2 = 2
          
          if (N_recvi>0 .or. N_sendi>0) then
             
             if     (myid+1==i) then
                
                if (present(tile_data_lH)) then
                                      
                   ! allocate tile_data_l_recv
                   
                   allocate(tile_data_l_recv(N_recvi,N_fields,N_ens))   
                   
                   ! obtain tile_data_l from processor j
                   ! send tile_data_l to processor j
                   
                   call MPI_SENDRECV(                                              &
                        tile_data_l,     N_sendi*N_fields*N_ens,MPI_real,j-1,tag1, &
                        tile_data_l_recv,N_recvi*N_fields*N_ens,MPI_real,j-1,tag2, &
                        mpicomm, mpistatus, mpierr ) 

                end if
                
             elseif (myid+1==j) then
                
                if (present(tile_data_lH)) then
                   
                   ! allocate tile_data_l_recv
                   
                   allocate(tile_data_l_recv(N_sendi,N_fields,N_ens)) ! N_sendi=N_recvj
                   
                   ! obtain tile_data_l from processor i
                   ! send tile_data_l to processor i
                   
                   call MPI_SENDRECV(                                              &
                        tile_data_l,     N_recvi*N_fields*N_ens,MPI_real,i-1,tag2, &
                        tile_data_l_recv,N_sendi*N_fields*N_ens,MPI_real,i-1,tag1, &
                        mpicomm, mpistatus, mpierr ) 
                   
                end if
                
             end if

             ! put received data into tile_data_lH
             
             if (myid+1==i .or. myid+1==j) then
                
                if (myid+1==i) then
                   
                   N_recv   = N_recvi
                   istart_f = low_ind(j)

                else   ! my_id+1==j
                   
                   N_recv   = N_sendi     ! N_sendi=N_recvj
                   istart_f = low_ind(i) 
                   
                end if

                iend_f   = istart_f + N_recv - 1

                istart   = N_catlH + 1
                
                iend     = istart + N_recv - 1

                ! append tile_data_l_recv to tile_data_lH 
                ! (similarly for tile_coord_lH)

                if (present(tile_coord_lH))                                          &
                     tile_coord_lH(istart:iend) = tile_coord_f(istart_f:iend_f)
                
                if (present(tile_data_lH))                                           &
                     tile_data_lH(istart:iend,1:N_fields,1:N_ens) = tile_data_l_recv
                
                N_catlH = iend
                
             end if
             
          end if
          
          if (allocated(tile_data_l_recv)) deallocate(tile_data_l_recv)
          
       end do
    end do
    
#else
    
    ! previous block kicks in only for numprocs>1 (not for OpenMP and sequential)
    
#endif
    
  end subroutine get_tiles_in_halo
  
  ! *********************************************************************
  
  subroutine get_obs_pert( N_ens, N_obs, N_obs_param,                  &
       pert_grid_f,                                                    &
       obs_param, Observations,                                        &
       Pert_rseed,                                                     &
       Obs_pert )

    ! reichle, 27 Jul 2005
    ! reichle,  3 Oct 2011 --- rewritten for obs w/in halo
    ! reichle, 14 Feb 2013 --- bug fix in call to get_pert()
    
    implicit none
    
    integer, intent(in) :: N_ens, N_obs, N_obs_param
    
    type(obs_param_type), dimension(N_obs_param), intent(in) :: obs_param
    
    type(obs_type), dimension(N_obs), intent(in) :: Observations
        
    type(grid_def_type), intent(in) :: pert_grid_f                
    
    integer, dimension(NRANDSEED,N_ens), intent(inout) :: Pert_rseed
    
    real, dimension(N_obs,N_ens), intent(out) :: Obs_pert

    ! --------------------------------
    
    ! local variables

    real, parameter :: dtstep = 0.
    
    type(pert_param_type), dimension(:), pointer :: obs_pert_param => null()

    real    :: this_lon, this_lat, delta_lon, delta_lat
    
    real    :: obs_minlon, obs_minlat, obs_maxlon, obs_maxlat

    integer :: i, j, k, lon_ind, lat_ind, N_assim_species, max_species_id

    integer :: ind_minlon, ind_minlat, ind_maxlon, ind_maxlat

    integer :: ind_i_min,  ind_j_min,  ind_i_max,  ind_j_max

    type(grid_def_type)                   :: pert_grid_lH

    integer, dimension(N_obs_param)       :: ind_assim_species
    
    integer, dimension(:),    allocatable :: ind_species2obsparam

    real, dimension(:,:,:,:), allocatable :: Obs_pert_ntrmdt, Obs_pert_grid
    
    character(len=*), parameter :: Iam = 'get_obs_pert'
    character(len=400) :: err_msg

    ! -----------------------------------------------------------------
    
    nullify(obs_pert_param)
    
    ! determine pert_grid_lH 
    !  - pert_grid_lH is the local grid for which perturbations are needed
    !  - pert_grid_lH is larger than pert_grid_l by the "halo"
    
    ! inherit select properties from pert_grid_f
    
    pert_grid_lH%gridtype = pert_grid_f%gridtype
    pert_grid_lH%ind_base = pert_grid_f%ind_base
    pert_grid_lH%i_dir    = pert_grid_f%i_dir   
    pert_grid_lH%j_dir    = pert_grid_f%j_dir   
        
    if (N_obs>0) then
       
       ! determine grid extent from Observations within local halo
       !
       ! inherit simple approach to dateline issue from get_halo_obs(), that is:
       !  obs_minlon>=-180, obs_maxlon<=180, obs_minlat>=-90, obs_maxlat<=90
       
       obs_minlon = minval(Observations%lon)
       obs_maxlon = maxval(Observations%lon)
       obs_minlat = minval(Observations%lat)
       obs_maxlat = maxval(Observations%lat)
    
       ! get i/j_ind of corner grid cells w.r.t. pert_grid_f
       
       call get_ij_ind_from_latlon( pert_grid_f, obs_minlat, obs_minlon, &
            ind_minlon, ind_minlat)
       
       call get_ij_ind_from_latlon( pert_grid_f, obs_maxlat, obs_maxlon, &
            ind_maxlon, ind_maxlat)
       
    else
       
       ! create a dummy 1-by-1 grid
       
       ind_minlon = 1
       ind_maxlon = 1
       ind_minlat = 1
       ind_maxlat = 1
       
    end if
    
    ! convert indices associated with lat/lon to min/max indices 
    ! (note that ind_minlon>=ind_maxlon if j_dir==-1)
    ! (note that ind_minlon, ind_minlat, etc will be needed again below)
        
    ind_i_min = min( ind_minlon, ind_maxlon)
    ind_i_max = max( ind_minlon, ind_maxlon)
    
    ind_j_min = min( ind_minlat, ind_maxlat)
    ind_j_max = max( ind_minlat, ind_maxlat)
    
    pert_grid_lH%N_lon  = ind_i_max - ind_i_min + 1
    pert_grid_lH%N_lat  = ind_j_max - ind_j_min + 1
   
    pert_grid_lH%i_offg = ind_i_min - 1 + pert_grid_f%i_offg
    pert_grid_lH%j_offg = ind_j_min - 1 + pert_grid_f%j_offg
    
    if     (index(pert_grid_lH%gridtype,'LatLon')/=0) then
       
       pert_grid_lH%dlon   = pert_grid_f%dlon
       pert_grid_lH%dlat   = pert_grid_f%dlat
              
       delta_lon = real(ind_minlon-1)*pert_grid_lH%dlon
       delta_lat = real(ind_minlat-1)*pert_grid_lH%dlat
       
       pert_grid_lH%ll_lon = pert_grid_f%ll_lon  + delta_lon
       pert_grid_lH%ll_lat = pert_grid_f%ll_lat  + delta_lat

       delta_lon = real(pert_grid_lH%N_lon)*pert_grid_lH%dlon
       delta_lat = real(pert_grid_lH%N_lat)*pert_grid_lH%dlat
       
       pert_grid_lH%ur_lon = pert_grid_lH%ll_lon + delta_lon
       pert_grid_lH%ur_lat = pert_grid_lH%ll_lat + delta_lat
       
    elseif ( index(pert_grid_lH%gridtype,'EASEv')  /=0 ) then
       
       pert_grid_lH%dlon   = pert_grid_f%dlon
       
       pert_grid_lH%dlat   = nodata_generic      ! not needed here
       
       pert_grid_lH%ll_lon = nodata_generic      ! not needed here
       pert_grid_lH%ll_lat = nodata_generic      ! not needed here
       pert_grid_lH%ur_lon = nodata_generic      ! not needed here
       pert_grid_lH%ur_lat = nodata_generic      ! not needed here
              
    else
       
       err_msg = 'not yet implemented for ' // pert_grid_lH%gridtype
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    ! ---------------------------------------------------------------------
    !
    ! Find out which species are assimilated in general (even there are
    ! none of a given assimilated species among the Observations at this time 
    ! and for this processor, or even none at all of a given species for this
    ! time step).
    ! This is necessary because this subroutine is called by all
    ! processors, and Pert_rseed can only stay consistent across processors
    ! if all processors generate random fields for the same number of 
    ! species.  (Alternatively, could implement a query across the Observations 
    ! from all processors, but this would require communications.)

    ind_assim_species = -9999

    max_species_id    = -9999
    
    j = 0
    
    do i=1,N_obs_param
       
       if (obs_param(i)%assim) then
          
          j = j+1
          
          ind_assim_species(j) = i
          
          max_species_id = max( max_species_id, obs_param(i)%species)

       end if
       
    end do

    N_assim_species = j

    ! get mapping from assimilated species counter (j) back to "species" id
    
    allocate(ind_species2obsparam(max_species_id))

    ind_species2obsparam = -9999
    
    do j=1,N_assim_species

       i = ind_assim_species(j)

       ind_species2obsparam(obs_param(i)%species) = j       
       
    end do
    
    ! --------------------------------------------------------------------
    !
    ! get obs_pert_param for use with get_pert()

    call allocate_pert_param( N_assim_species, &
         pert_grid_lH%N_lon, pert_grid_lH%N_lat, obs_pert_param)
    
    do j=1,N_assim_species
       
       ! Use land_pert module to generate spatially correlated and 
       ! standard-normally distributed perturbations.
       !
       ! For observation perturbations, always use:
       !
       !   tcorr = 0.  (never temporally correlated)
       !   typ   = 0   (always additive)
       !   ccorr = 0.  (never cross-correlated)
       !
       
       i = ind_assim_species(j)

       obs_pert_param(j)%descr           = obs_param(i)%descr
       obs_pert_param(j)%typ             = 0     
       obs_pert_param(j)%std_normal_max  = obs_param(i)%std_normal_max
       obs_pert_param(j)%zeromean        = obs_param(i)%zeromean
       obs_pert_param(j)%coarsen         = obs_param(i)%coarsen_pert 

       obs_pert_param(j)%mean(:,:)       = 0.   ! need N(0,1) perturbations
       obs_pert_param(j)%std( :,:)       = 1.   ! need N(0,1) perturbations
       
       obs_pert_param(j)%ccorr(:,:,:)    = 0.
       obs_pert_param(j)%ccorr(j,:,:)    = 1.
       
       obs_pert_param(j)%xcorr           = obs_param(i)%xcorr
       obs_pert_param(j)%ycorr           = obs_param(i)%ycorr         
       obs_pert_param(j)%tcorr           = 0
       
    end do
    
    ! --------------------------------------------------------------------
    !
    ! get gridded perturbations
    
    allocate(Obs_pert_ntrmdt(pert_grid_lH%N_lon, pert_grid_lH%N_lat, N_assim_species, N_ens))
    allocate(Obs_pert_grid(  pert_grid_lH%N_lon, pert_grid_lH%N_lat, N_assim_species, N_ens))

    call get_pert(                                                        &
         N_assim_species, N_assim_species, N_ens,                         &
         pert_grid_lH, pert_grid_f,                                       &  ! switched order (reichle, 17 Jul 2020)
         dtstep,                                                          &
         obs_pert_param,                                                  &
         Pert_rseed,                                                      &
         Obs_pert_ntrmdt,                                                 &
         Obs_pert_grid )

    ! clean up

    deallocate(Obs_pert_ntrmdt)

    call deallocate_pert_param(N_assim_species, obs_pert_param)
    
    ! --------------------------------------------------------------------
    !
    ! map gridded perturbations to observations one at a time,
    ! scale to desired error standard deviation
    
    do k=1,N_obs
       
       ! map from grid to obs
       
       this_lon = Observations(k)%lon
       this_lat = Observations(k)%lat
       
       call get_ij_ind_from_latlon( pert_grid_lH, this_lat, this_lon, lon_ind, lat_ind )
       
       j = ind_species2obsparam(Observations(k)%species)
       
       Obs_pert(k,1:N_ens) = Obs_pert_grid( lon_ind, lat_ind, j, 1:N_ens )
       
       Obs_pert(k,1:N_ens) = Obs_pert(k,1:N_ens) * sqrt(Observations(k)%obsvar)
       
    end do

    deallocate(Obs_pert_grid)

    deallocate(ind_species2obsparam)
    
    ! -----------------------------------------------------------------
    !
    ! enforce physical constraints 
    
    ! Skip this step.  It has only been used for soil moisture content,
    ! where a perturbation may result in an observation outside of the
    ! physical range, but because of the physical checks on the analysis
    ! this step should not be necessary.  
    ! reichle,  3 Oct 2011

    !call check_obs_pert( N_ens, N_catd, N_obs, cat_param, Observations, &
    !     Obs_pert )
    
  end subroutine get_obs_pert
  
  ! *********************************************************************

  subroutine cat_enkf_increments(                               &
       N_ens, N_obs, N_catd, N_obs_param,                       &
       update_type, obs_param,                                  &
       tile_coord, l2f,                                         &
       Observations, Obs_pred, Obs_pert,                        &
       cat_param,                                               &
       xcompact, ycompact, fcsterr_inflation_fac,               &
       cat_progn, cat_progn_incr )
    
    ! get increments for Catchment prognostic variables
    !
    ! reichle, 27 Jan 2005 - eliminated update of Mod_err
    ! reichle, 27 Jul 2005
    ! reichle, 18 Oct 2005 - return increments (instead of updated cat_progn)
    ! reichle, 17 Oct 2011 - added "l2f" for revised (MPI) analysis
    ! reichle, 20 Feb 2022 - modified update_type 10 for PEATCLSM
    !
    ! --------------------------------------------------------------
    
    ! IMPORTANT:
    ! on input, cat_progn must contain cat_progn_minus(1:N_catd,1:N_ens)
    ! on output, cat_progn_incr contains INCREMENTS
    
    ! type of update is selected by "update_type"
    
    ! *********************************************************************
    ! **************** WARNING  WARNING  WARNING  WARNING  ****************
    ! ********************************************************************* 
    !
    ! "update_types" 1-5 below have NOT been fully tested after extensive
    ! revisions re. obs handling and MPI implementation.
    ! - reichle, 17 Oct 2011

    ! -------------------------------------------------------------------
    
    implicit none
    
    ! inputs
    
    integer, intent(in) :: N_ens, N_obs, N_catd, N_obs_param, update_type

    type(obs_param_type), dimension(N_obs_param), intent(in) :: obs_param               

    type(tile_coord_type), dimension(:), pointer  :: tile_coord     ! input

    integer, dimension(N_catd), intent(in) :: l2f
    
    type(obs_type), intent(in), dimension(N_obs)  :: Observations   ! input
    
    real, intent(in), dimension(N_obs,N_ens)  :: Obs_pred
    real, intent(in), dimension(N_obs,N_ens)  :: Obs_pert
    
    type(cat_param_type), dimension(N_catd), intent(in) :: cat_param
    
    real, intent(in) :: xcompact, ycompact, fcsterr_inflation_fac
    
    type(cat_progn_type), intent(in),  dimension(N_catd,N_ens) :: cat_progn
    
    ! output
    
    type(cat_progn_type), intent(out), dimension(N_catd,N_ens) :: &
         cat_progn_incr
    
    ! -----------------------------
    !
    ! locals

    ! thresholds for identifying snow-free and non-frozen tiles
    ! NOT YET USED
    ! - reichle, 14 Oct 2014

    real, parameter :: SWE_threshold = +HUGE(1.) ! = 1.e-4           ! [kg/m2]
    
    real, parameter :: tp1_threshold = -HUGE(1.) ! = 0.2             ! [CELSIUS]

    integer :: n, n_e, kk, ii
    
    integer :: N_state_max, N_state, N_selected_obs, N_select_varnames, N_select_species
    
    real    :: halo_minlon, halo_maxlon, halo_minlat, halo_maxlat
    real    :: tmp_minlon,  tmp_maxlon,  tmp_minlat,  tmp_maxlat

    real    :: tmp_obs, deltaT

    real    :: fice_minus, tp1_minus, ght1_minus
    real    :: fice_plus,  tp1_plus,  ght1_plus
    
    integer,           dimension(N_obs)   :: ind_obs
    
    real, allocatable, dimension(:,:)     :: State_incr
    real, allocatable, dimension(:,:)     :: Obs_cov      ! measurement error covariance
    
    real, allocatable, dimension(:)       :: State_lon, State_lat

    integer,       dimension(N_obs_param) :: select_species  ! alloc max possible length

    character(40), dimension(N_obs_param) :: select_varnames ! alloc max possible length
    
    integer, dimension(:), allocatable    :: select_tilenum
    
    integer, dimension(:,:),   pointer    :: N_tile_in_cell_ij      => null()
    integer, dimension(:,:,:), pointer    :: tile_num_in_cell_ij    => null()

    character(len=*),  parameter          :: Iam = 'cat_enkf_increments'
    character(len=400)                    :: err_msg

    real, dimension(     N_catd)          :: r_x, tmp_dlon
    real                                  :: r_y, tmp_dlat
    
    real, dimension(     N_catd)          :: srfexc,   rzexc, catdef
    real, dimension(     N_catd)          :: ar1,      ar2,   ar4
    real, dimension(     N_catd)          :: FT_state
    
    real, dimension(     N_catd,N_ens)    :: sfmc,     rzmc,  prmc
    real, dimension(     N_catd,N_ens)    :: tsurf,    tsurf_excl_snow
    real, dimension(     N_catd,N_ens)    :: SWE,      asnow
    real, dimension(     N_catd,N_ens)    :: FT_Teff
    
    real, dimension(N_gt,N_catd,N_ens)    :: tp,       fice

    real, dimension(     N_catd)          :: tsurf_ensavg
    real, dimension(     N_catd)          :: SWE_ensavg
    real, dimension(     N_catd)          :: tp1_ensavg

    type(obs_param_type)                  :: this_obs_param
        
    ! -----------------------------------------------------------------------

    if (logit) write (logunit,*) &
         'cat_enkf_increments(): getting assimilation increments...' 
    
    ! initialize - needed regardless of (local) N_obs
    
    do kk=1,N_catd
       do n_e=1,N_ens
          cat_progn_incr(kk,n_e) = 0.
       end do
    end do
    
    ! avoid unnecessary work or subroutine calls

    if (N_obs<=0) return   ! nothing left to do
    
    ! more initializations

    N_select_varnames = 0
    N_select_species  = 0

    select_varnames   = ''
    select_species    = -8888  ! intentionally differs from init in get_select_species()
    
    ! ----------------------------------------------------------------------
    !
    ! IMPORTANT: do *NOT* add MPI calls to the remainder of this subroutine
    !            or to subroutines called from there because not all 
    !            processors continue past this point
    !
    ! In future, perhaps make sure that all processors can safely proceed, 
    ! which is not clear right now. reichle, 26 Sep 2013
    ! 
    ! ----------------------------------------------------------------------
    
    ! compute soil temperature and snow diagnostics for
    ! - screening tiles for which increments should not be computed (typically 3d updates)
    ! - FT analysis
    
    ! total SWE
    
    do kk=1,N_catd
       do n_e=1,N_ens
          
          SWE(kk,n_e) = sum( cat_progn(kk,n_e)%wesn(1:N_snow) )
          
       end do
    end do

    ! soil moisture, tsurf, and soil temperature diagnostics

    do n_e=1,N_ens

       ! soil moisture

       srfexc = cat_progn(:,n_e)%srfexc
       rzexc  = cat_progn(:,n_e)%rzexc
       catdef = cat_progn(:,n_e)%catdef

       call catch_calc_soil_moist(                                            &
            N_catd,           cat_param%dzsf,   cat_param%vgwmax,             &
            cat_param%cdcr1,  cat_param%cdcr2,  cat_param%psis,               &
            cat_param%bee,    cat_param%poros,  cat_param%wpwet,              &
            cat_param%ars1,   cat_param%ars2,   cat_param%ars3,               &
            cat_param%ara1,   cat_param%ara2,   cat_param%ara3,               &
            cat_param%ara4,   cat_param%arw1,   cat_param%arw2,               &
            cat_param%arw3,   cat_param%arw4,                                 &
            cat_param%bf1,    cat_param%bf2,                                  &
            srfexc, rzexc, catdef, ar1, ar2, ar4,                             &
            sfmc(:,n_e), rzmc(:,n_e), prmc(:,n_e) )
       
       ! tsurf
       
       call catch_calc_tsurf( N_catd,                                         &
            cat_progn(:,n_e)%tc1, cat_progn(:,n_e)%tc2, cat_progn(:,n_e)%tc4, &
            catprogn2wesn(N_catd,cat_progn(:,n_e)),                           &
            catprogn2htsn(N_catd,cat_progn(:,n_e)),                           &
            ar1,  ar2,  ar4, tsurf(:,n_e) )

       ! tsurf excluding snow
       
       call catch_calc_tsurf_excl_snow( N_catd,                               &
            cat_progn(:,n_e)%tc1, cat_progn(:,n_e)%tc2, cat_progn(:,n_e)%tc4, &
            ar1,  ar2,  ar4, tsurf_excl_snow(:,n_e) )
       
       ! soil temperature 
       !
       ! NOTE: "tp" is returned in CELSIUS
       
       call catch_calc_tp( N_catd, cat_param%poros,                           &
            catprogn2ghtcnt(N_catd,cat_progn(:,n_e)),                         &
            tp(:,:,n_e), fice(:,:,n_e) )
       
       ! snow cover fraction
       
       call StieglitzSnow_calc_asnow( N_snow, N_catd,                         &
            catprogn2wesn(N_catd,cat_progn(:,n_e)), asnow(:,n_e) )       
       
    end do

    ! compute ensemble average of select variables

    SWE_ensavg   = 0.
    tsurf_ensavg = 0.
    tp1_ensavg   = 0.
    
    do n_e=1,N_ens
       
       SWE_ensavg   = SWE_ensavg   + SWE(    :,n_e)
       tsurf_ensavg = tsurf_ensavg + tsurf(  :,n_e)
       tp1_ensavg   = tp1_ensavg   + tp(   1,:,n_e)
       
    end do
    
    SWE_ensavg   = SWE_ensavg   /real(N_ens)
    tsurf_ensavg = tsurf_ensavg /real(N_ens)
    tp1_ensavg   = tp1_ensavg   /real(N_ens)
    
    ! ---------------------------------------------------------------------

    select_update_type: select case (update_type)
       
    case (1) select_update_type   ! 1d soil moisture analysis; sfmc obs

       ! this 1d update requires that obs are on same tile space as model
       
       if (logit) write (logunit,*) 'get 1d soil moisture increments; sfmc obs'

       N_select_varnames  = 1
       
       select_varnames(1) = 'sfmc'
       
       call get_select_species(                                           &
            N_select_varnames, select_varnames(1:N_select_varnames),      &
            N_obs_param, obs_param, N_select_species, select_species )
       
       N_state = 3
       
       allocate( State_incr(N_state,N_ens) )
       allocate( select_tilenum(1))
       
       allocate(Obs_cov(N_obs,N_obs))
       
       call assemble_obs_cov( N_obs, N_obs_param, obs_param, &
            Observations(1:N_obs), Obs_cov )

       do n=1,N_catd
          
          ! find observations for catchment n
          
          select_tilenum(1) = l2f(n)
          
          call get_ind_obs(                                           &
               N_obs,            Observations,                        &
               1,                select_tilenum,                      &
               N_select_species, select_species(1:N_select_species),  &
               N_selected_obs,   ind_obs )
          
          if (N_selected_obs > 0) then
             
             ! assemble State_minus
             ! (on input, cat_progn contains cat_progn_minus)
             
             State_incr(1,:) = cat_progn(n,:)%srfexc/scale_srfexc
             State_incr(2,:) = cat_progn(n,:)%rzexc /scale_rzexc
             State_incr(3,:) = cat_progn(n,:)%catdef/scale_catdef
             
             ! EnKF update
             
             call enkf_increments(                                        &
                  N_state, N_selected_obs, N_ens,                         &
                  Observations(ind_obs(1:N_selected_obs)),                &
                  Obs_pred(ind_obs(1:N_selected_obs),:),                  & 
                  Obs_pert(ind_obs(1:N_selected_obs),:),                  & 
                  Obs_cov(                                                &
                  ind_obs(1:N_selected_obs), ind_obs(1:N_selected_obs)),  &
                  State_incr,                                             &
                  fcsterr_inflation_fac=fcsterr_inflation_fac )
             
             ! assemble cat_progn increments
             
             cat_progn_incr(n,:)%srfexc = State_incr(1,:)*scale_srfexc
             cat_progn_incr(n,:)%rzexc  = State_incr(2,:)*scale_rzexc
             cat_progn_incr(n,:)%catdef = State_incr(3,:)*scale_catdef
             
          end if
          
       end do
       
    case (2) select_update_type   ! 3d soil moisture analysis; sfmc obs
       
       ! update each tile separately using all observations within 
       ! the customized halo around each tile
       
       if (logit) write (logunit,*) 'get 3d soil moisture increments; sfmc obs'

       N_select_varnames  = 1
       
       select_varnames(1) = 'sfmc'

       call get_select_species(                                           &
            N_select_varnames, select_varnames(1:N_select_varnames),      &
            N_obs_param, obs_param, N_select_species, select_species )
              
       N_state = 3
       
       allocate( State_incr(N_state,N_ens)) 
       allocate( State_lon( N_state      ))
       allocate( State_lat( N_state      ))
       
       do kk=1,N_catd
          
          ! find observations within halo around tile kk
          
          halo_minlon = tile_coord(kk)%com_lon - xcompact
          halo_maxlon = tile_coord(kk)%com_lon + xcompact
          halo_minlat = tile_coord(kk)%com_lat - ycompact
          halo_maxlat = tile_coord(kk)%com_lat + ycompact
          
          ! simple approach to dateline issue (cut halo back to at most -180:180, -90:90)
          ! - reichle, 28 May 2013
          
          halo_minlon = max(halo_minlon,-180.)
          halo_maxlon = min(halo_maxlon, 180.)
          halo_minlat = max(halo_minlat, -90.)
          halo_maxlat = min(halo_maxlat,  90.)
          
          call get_ind_obs_lat_lon_box(                               &
               N_obs,            Observations,                        &
               halo_minlon, halo_maxlon, halo_minlat, halo_maxlat,    &
               N_select_species, select_species(1:N_select_species),  &
               N_selected_obs,   ind_obs )
          
          if (N_selected_obs>0) then
             
             ! assemble State_minus
             ! (on input, cat_progn contains cat_progn_minus)
             
             State_incr(1,:) = cat_progn( kk,:)%srfexc/scale_srfexc
             State_incr(2,:) = cat_progn( kk,:)%rzexc /scale_rzexc
             State_incr(3,:) = cat_progn( kk,:)%catdef/scale_catdef

             State_lon(   :) = tile_coord(kk  )%com_lon
             State_lat(   :) = tile_coord(kk  )%com_lat
             
             allocate(Obs_cov(N_selected_obs,N_selected_obs))

             call assemble_obs_cov( N_selected_obs, N_obs_param, obs_param, &
                  Observations(ind_obs(1:N_selected_obs)), Obs_cov )
             
             ! EnKF update
             
             call enkf_increments(                                        &
                  N_state, N_selected_obs, N_ens,                         &
                  Observations(ind_obs(1:N_selected_obs)),                &
                  Obs_pred(ind_obs(1:N_selected_obs),:),                  &
                  Obs_pert(ind_obs(1:N_selected_obs),:),                  &
                  Obs_cov,                                                &
                  State_incr, State_lon, State_lat, xcompact, ycompact,   &
                  fcsterr_inflation_fac )             
             
             deallocate(Obs_cov)
             
             ! assemble cat_progn increments
             
             cat_progn_incr(kk,:)%srfexc = State_incr(1,:)*scale_srfexc
             cat_progn_incr(kk,:)%rzexc  = State_incr(2,:)*scale_rzexc
             cat_progn_incr(kk,:)%catdef = State_incr(3,:)*scale_catdef
             
          end if
          
       end do
       
       ! ----------------------------------
              
    case (3) select_update_type   ! 1d Tskin analysis; tskin obs
       
       ! update_type = 3: 1d Tskin (incr NOT applied, use w/ cat bias corr)

       ! this 1d update requires that obs are on same tile space as model
       
       if (logit) write (logunit,*) 'get 1d Tskin increments; tskin obs'
              
       N_select_varnames  = 1
       
       select_varnames(1) = 'tskin'

       call get_select_species(                                           &
            N_select_varnames, select_varnames(1:N_select_varnames),      &
            N_obs_param, obs_param, N_select_species, select_species )
       
       N_state = 3
       
       allocate( State_incr(N_state,N_ens) )
       allocate( select_tilenum(1))

       allocate(Obs_cov(N_obs,N_obs))
       
       call assemble_obs_cov( N_obs, N_obs_param, obs_param, &
            Observations(1:N_obs), Obs_cov )       

       do n=1,N_catd
          
          ! find observations for catchment n
          
          select_tilenum(1) = l2f(n)
          
          call get_ind_obs(                                           &
               N_obs,            Observations,                        &
               1,                select_tilenum,                      &
               N_select_species, select_species(1:N_select_species),  &
               N_selected_obs,   ind_obs )

          if (N_selected_obs > 0) then
             
             ! assemble State_minus
             ! (on input, cat_progn contains cat_progn_minus)
             
             State_incr(1,:) = cat_progn(n,:)%tc1/scale_temp
             State_incr(2,:) = cat_progn(n,:)%tc2/scale_temp
             State_incr(3,:) = cat_progn(n,:)%tc4/scale_temp
             
             ! EnKF update
             
             call enkf_increments(                                        &
                  N_state, N_selected_obs, N_ens,                         &
                  Observations(ind_obs(1:N_selected_obs)),                &
                  Obs_pred(ind_obs(1:N_selected_obs),:),                  & 
                  Obs_pert(ind_obs(1:N_selected_obs),:),                  & 
                  Obs_cov(                                                &
                  ind_obs(1:N_selected_obs), ind_obs(1:N_selected_obs)),  &
                  State_incr,                                             &
                  fcsterr_inflation_fac=fcsterr_inflation_fac )
             
             ! assemble cat_progn increments
                          
             cat_progn_incr(n,:)%tc1 = State_incr(1,:)*scale_temp
             cat_progn_incr(n,:)%tc2 = State_incr(2,:)*scale_temp
             cat_progn_incr(n,:)%tc4 = State_incr(3,:)*scale_temp
             
          end if
          
       end do
       
       ! ----------------------------------

    case (4,5) select_update_type   ! 1d Tskin/ght(1) analysis; tskin obs 
       
       ! update_type = 4: 1d Tskin/ght(1) (incr applied, use w/ or w/o cat bias corr)
       ! update_type = 5: 1d Tskin/ght(1) (incr NOT applied, use w/ cat bias corr)

       ! this 1d update requires that obs are on same tile space as model
       
       if (logit) write (logunit,*) 'get 1d Tskin/ght(1) increments; tskin obs'
       
       N_select_varnames  = 1
       
       select_varnames(1) = 'tskin'

       call get_select_species(                                           &
            N_select_varnames, select_varnames(1:N_select_varnames),      &
            N_obs_param, obs_param, N_select_species, select_species )
       
       N_state = 4
       
       allocate( State_incr(N_state,N_ens) )
       allocate( select_tilenum(1))
       
       allocate(Obs_cov(N_obs,N_obs))
       
       call assemble_obs_cov( N_obs, N_obs_param, obs_param, &
            Observations(1:N_obs), Obs_cov )

       do n=1,N_catd
          
          ! find observations for catchment n
          
          select_tilenum(1) = l2f(n)
          
          call get_ind_obs(                                           &
               N_obs,            Observations,                        &
               1,                select_tilenum,                      &
               N_select_species, select_species(1:N_select_species),  &
               N_selected_obs,   ind_obs )

          if (N_selected_obs > 0) then
             
             ! assemble State_minus
             ! (on input, cat_progn contains cat_progn_minus)
             
             State_incr(1,:) = cat_progn(n,:)%tc1/scale_temp
             State_incr(2,:) = cat_progn(n,:)%tc2/scale_temp
             State_incr(3,:) = cat_progn(n,:)%tc4/scale_temp
             State_incr(4,:) = cat_progn(n,:)%ght(1)/scale_ght1
             
             ! EnKF update
             
             call enkf_increments(                                        &
                  N_state, N_selected_obs, N_ens,                         &
                  Observations(ind_obs(1:N_selected_obs)),                &
                  Obs_pred(ind_obs(1:N_selected_obs),:),                  & 
                  Obs_pert(ind_obs(1:N_selected_obs),:),                  & 
                  Obs_cov(                                                &
                  ind_obs(1:N_selected_obs), ind_obs(1:N_selected_obs)),  &
                  State_incr,                                             &
                  fcsterr_inflation_fac=fcsterr_inflation_fac )
             
             ! assemble cat_progn increments
                          
             cat_progn_incr(n,:)%tc1    = State_incr(1,:)*scale_temp
             cat_progn_incr(n,:)%tc2    = State_incr(2,:)*scale_temp
             cat_progn_incr(n,:)%tc4    = State_incr(3,:)*scale_temp
             cat_progn_incr(n,:)%ght(1) = State_incr(4,:)*scale_ght1
             
          end if
          
       end do

       ! ----------------------------------

    case (6) select_update_type   ! 1d soil moisture/Tskin/ght(1) analysis; Tb obs

       ! this 1d update requires that obs are on same tile space as model
       
       if (logit) write (logunit,*) 'get 1d soil moisture/Tskin/ght(1) increments; Tb obs'

       N_select_varnames  = 1
       
       select_varnames(1) = 'Tb'

       call get_select_species(                                           &
            N_select_varnames, select_varnames(1:N_select_varnames),      &
            N_obs_param, obs_param, N_select_species, select_species )

       N_state = 7  
       
       allocate( State_incr(N_state,N_ens) )
       allocate( select_tilenum(1))

       allocate(Obs_cov(N_obs,N_obs))
       
       call assemble_obs_cov( N_obs, N_obs_param, obs_param, &
            Observations(1:N_obs), Obs_cov )

       do n=1,N_catd

          ! find observations for catchment n

          select_tilenum(1) = l2f(n)
          
          call get_ind_obs(                                           &
               N_obs,            Observations,                        &
               1,                select_tilenum,                      &
               N_select_species, select_species(1:N_select_species),  &
               N_selected_obs,   ind_obs )

          if (N_selected_obs > 0) then
             
             ! assemble State_minus
             ! (on input, cat_progn contains cat_progn_minus)

             State_incr(1,:) = cat_progn(n,:)%srfexc/scale_srfexc
             State_incr(2,:) = cat_progn(n,:)%rzexc /scale_rzexc
             State_incr(3,:) = cat_progn(n,:)%catdef/scale_catdef

             State_incr(4,:) = cat_progn(n,:)%tc1/scale_temp
             State_incr(5,:) = cat_progn(n,:)%tc2/scale_temp
             State_incr(6,:) = cat_progn(n,:)%tc4/scale_temp
             State_incr(7,:) = cat_progn(n,:)%ght(1)/scale_ght1

             ! EnKF update
             
             call enkf_increments(                                        &
                  N_state, N_selected_obs, N_ens,                         &
                  Observations(ind_obs(1:N_selected_obs)),                &
                  Obs_pred(ind_obs(1:N_selected_obs),:),                  &
                  Obs_pert(ind_obs(1:N_selected_obs),:),                  &
                  Obs_cov(                                                &
                  ind_obs(1:N_selected_obs), ind_obs(1:N_selected_obs)),  &
                  State_incr,                                             &
                  fcsterr_inflation_fac=fcsterr_inflation_fac )

             ! assemble cat_progn increments

             cat_progn_incr(n,:)%srfexc = State_incr(1,:)*scale_srfexc
             cat_progn_incr(n,:)%rzexc  = State_incr(2,:)*scale_rzexc
             cat_progn_incr(n,:)%catdef = State_incr(3,:)*scale_catdef

             cat_progn_incr(n,:)%tc1    = State_incr(4,:)*scale_temp
             cat_progn_incr(n,:)%tc2    = State_incr(5,:)*scale_temp
             cat_progn_incr(n,:)%tc4    = State_incr(6,:)*scale_temp
             cat_progn_incr(n,:)%ght(1) = State_incr(7,:)*scale_ght1

          end if

       end do
       
       ! ----------------------------------
       
    case (7) select_update_type   ! 3d Tskin/ght(1) analysis; tskin obs 
       
       ! update each tile separately using all observations within 
       ! the customized halo around each tile
       !
       ! replaces previous approach ("3d update over each grid cell")
       ! - reichle, 26 March 2014

       if (logit) write (logunit,*) 'get 3d Tskin/ght(1) increments; tskin obs'

       N_select_varnames  = 1
       
       select_varnames(1) = 'tskin'
       
       call get_select_species(                                           &
            N_select_varnames, select_varnames(1:N_select_varnames),      &
            N_obs_param, obs_param, N_select_species, select_species )
       
       N_state = 4
       
       allocate( State_incr(N_state,N_ens)) 
       allocate( State_lon( N_state      ))
       allocate( State_lat( N_state      ))
       
       do kk=1,N_catd
          
          ! find observations within halo around tile kk
          
          halo_minlon = tile_coord(kk)%com_lon - xcompact
          halo_maxlon = tile_coord(kk)%com_lon + xcompact
          halo_minlat = tile_coord(kk)%com_lat - ycompact
          halo_maxlat = tile_coord(kk)%com_lat + ycompact
          
          ! simple approach to dateline issue (cut halo back to at most -180:180, -90:90)
          ! - reichle, 28 May 2013
          
          halo_minlon = max(halo_minlon,-180.)
          halo_maxlon = min(halo_maxlon, 180.)
          halo_minlat = max(halo_minlat, -90.)
          halo_maxlat = min(halo_maxlat,  90.)
          
          call get_ind_obs_lat_lon_box(                               &
               N_obs,            Observations,                        &
               halo_minlon, halo_maxlon, halo_minlat, halo_maxlat,    &
               N_select_species, select_species(1:N_select_species),  &
               N_selected_obs,   ind_obs )
          
          if (N_selected_obs>0) then
                          
             ! assemble State_minus
             ! (on input, cat_progn contains cat_progn_minus)
             
             State_incr(1,:) = cat_progn( kk,:)%tc1   /scale_temp
             State_incr(2,:) = cat_progn( kk,:)%tc2   /scale_temp
             State_incr(3,:) = cat_progn( kk,:)%tc4   /scale_temp
             State_incr(4,:) = cat_progn( kk,:)%ght(1)/scale_ght1

             State_lon(   :) = tile_coord(kk  )%com_lon
             State_lat(   :) = tile_coord(kk  )%com_lat
             
             allocate(Obs_cov(N_selected_obs,N_selected_obs))
             
             call assemble_obs_cov( N_selected_obs, N_obs_param, obs_param, &
                  Observations(ind_obs(1:N_selected_obs)), Obs_cov )
             
             ! EnKF update
             
             call enkf_increments(                                        &
                  N_state, N_selected_obs, N_ens,                         &
                  Observations(ind_obs(1:N_selected_obs)),                &
                  Obs_pred(ind_obs(1:N_selected_obs),:),                  &
                  Obs_pert(ind_obs(1:N_selected_obs),:),                  &
                  Obs_cov,                                                &
                  State_incr, State_lon, State_lat, xcompact, ycompact,   &
                  fcsterr_inflation_fac )             
             
             deallocate(Obs_cov)

             ! assemble cat_progn increments
             
             cat_progn_incr(kk,:)%tc1    = State_incr(1,:)*scale_temp
             cat_progn_incr(kk,:)%tc2    = State_incr(2,:)*scale_temp
             cat_progn_incr(kk,:)%tc4    = State_incr(3,:)*scale_temp
             cat_progn_incr(kk,:)%ght(1) = State_incr(4,:)*scale_ght1
             
          end if
          
       end do

       ! ----------------------------------
       
    case (8,10) select_update_type   ! 3d soil moisture/Tskin/ght(1) analysis; Tb obs

       ! update each tile separately using all observations within customized halo around each tile
       !
       ! state vector includes different subsets of Catchment model soil moisture prognostics:
       !
       !  update_type | subset of tiles | state vector
       !  ===================================================================================================
       !       8      | all             | srfexc, rzexc, catdef, tc1, tc2, tc4, ght1
       !  ---------------------------------------------------------------------------------------------------
       !      10      | PEATCLSM tiles  | srfexc, rzexc, catdef, tc1, tc2, tc4, ght1
       !              | otherwise       | srfexc, rzexc,         tc1, tc2, tc4, ght1  (incl. NLv4 peat tiles)
       !  ---------------------------------------------------------------------------------------------------
       !
       ! reichle, 27 Nov 2017
       ! reichle, 20 Feb 2022 - modified for PEATCLSM
       
       if (logit) write (logunit,*) 'get 3d soil moisture/Tskin/ght(1) increments; Tb obs'

       N_select_varnames  = 1
       
       select_varnames(1) = 'Tb'
       
       call get_select_species(                                           &
            N_select_varnames, select_varnames(1:N_select_varnames),      &
            N_obs_param, obs_param, N_select_species, select_species )
       
       N_state_max = 7
       
       allocate( State_incr(N_state_max,N_ens)) 
       allocate( State_lon( N_state_max      ))
       allocate( State_lat( N_state_max      ))
       
       do kk=1,N_catd       
          
          ! compute increments only snow-free and non-frozen tiles
          
          if ( (SWE_ensavg(kk) < SWE_threshold)     .and.            &
               (tp1_ensavg(kk) > tp1_threshold)             ) then  
             
             ! find observations within halo around tile kk
             
             halo_minlon = tile_coord(kk)%com_lon - xcompact
             halo_maxlon = tile_coord(kk)%com_lon + xcompact
             halo_minlat = tile_coord(kk)%com_lat - ycompact
             halo_maxlat = tile_coord(kk)%com_lat + ycompact
             
             ! simple approach to dateline issue (cut halo back to at most -180:180, -90:90)
             ! - reichle, 28 May 2013
             
             halo_minlon = max(halo_minlon,-180.)
             halo_maxlon = min(halo_maxlon, 180.)
             halo_minlat = max(halo_minlat, -90.)
             halo_maxlat = min(halo_maxlat,  90.)
             
             call get_ind_obs_lat_lon_box(                               &
                  N_obs,            Observations,                        &
                  halo_minlon, halo_maxlon, halo_minlat, halo_maxlat,    &
                  N_select_species, select_species(1:N_select_species),  &
                  N_selected_obs,   ind_obs )
             
             if (N_selected_obs>0) then
                
                if ( (update_type== 8                                                    ) .or. &
                     (update_type==10 .and. cat_param(kk)%poros>=PEATCLSM_POROS_THRESHOLD)      &
                     ) then
                   
                   N_state = 7   ! srfexc, rzexc, catdef, tc1, tc2, tc4, ght1
                   
                else
                   
                   N_state = 6   ! srfexc, rzexc,         tc1, tc2, tc4, ght1
                   
                end if
                
                ! assemble State_minus
                ! (on input, cat_progn contains cat_progn_minus)
                
                if ( N_state==7 ) then
                   
                   State_incr(1,:) = cat_progn( kk,:)%srfexc/scale_srfexc
                   State_incr(2,:) = cat_progn( kk,:)%rzexc /scale_rzexc
                   State_incr(3,:) = cat_progn( kk,:)%catdef/scale_catdef   ! catdef in State
                   
                   State_incr(4,:) = cat_progn( kk,:)%tc1   /scale_temp
                   State_incr(5,:) = cat_progn( kk,:)%tc2   /scale_temp
                   State_incr(6,:) = cat_progn( kk,:)%tc4   /scale_temp
                   State_incr(7,:) = cat_progn( kk,:)%ght(1)/scale_ght1
                   
                else
                   
                   State_incr(1,:) = cat_progn( kk,:)%srfexc/scale_srfexc
                   State_incr(2,:) = cat_progn( kk,:)%rzexc /scale_rzexc
                   
                   State_incr(3,:) = cat_progn( kk,:)%tc1   /scale_temp
                   State_incr(4,:) = cat_progn( kk,:)%tc2   /scale_temp
                   State_incr(5,:) = cat_progn( kk,:)%tc4   /scale_temp
                   State_incr(6,:) = cat_progn( kk,:)%ght(1)/scale_ght1

                end if
                
                State_lon(   :) = tile_coord(kk  )%com_lon
                State_lat(   :) = tile_coord(kk  )%com_lat
                
                allocate(Obs_cov(N_selected_obs,N_selected_obs))
                
                call assemble_obs_cov( N_selected_obs, N_obs_param, obs_param, &
                     Observations(ind_obs(1:N_selected_obs)), Obs_cov )
                
                ! EnKF update
                
                call enkf_increments(                                        &
                     N_state, N_selected_obs, N_ens,                         &
                     Observations(ind_obs(1:N_selected_obs)),                &
                     Obs_pred(ind_obs(1:N_selected_obs),:),                  &
                     Obs_pert(ind_obs(1:N_selected_obs),:),                  &
                     Obs_cov,                                                &
                     State_incr(1:N_state,:),                                &
                     State_lon( 1:N_state  ),                                &
                     State_lat( 1:N_state  ),                                &
                     xcompact, ycompact,                                     &
                     fcsterr_inflation_fac )             
                
                deallocate(Obs_cov)
                
                ! assemble cat_progn increments
                
                if (N_state==7) then
                   
                   cat_progn_incr(kk,:)%srfexc = State_incr(1,:)*scale_srfexc
                   cat_progn_incr(kk,:)%rzexc  = State_incr(2,:)*scale_rzexc
                   cat_progn_incr(kk,:)%catdef = State_incr(3,:)*scale_catdef   ! catdef in State
                   
                   cat_progn_incr(kk,:)%tc1    = State_incr(4,:)*scale_temp
                   cat_progn_incr(kk,:)%tc2    = State_incr(5,:)*scale_temp
                   cat_progn_incr(kk,:)%tc4    = State_incr(6,:)*scale_temp
                   cat_progn_incr(kk,:)%ght(1) = State_incr(7,:)*scale_ght1
                
                else
                   
                   cat_progn_incr(kk,:)%srfexc = State_incr(1,:)*scale_srfexc
                   cat_progn_incr(kk,:)%rzexc  = State_incr(2,:)*scale_rzexc
                   
                   cat_progn_incr(kk,:)%tc1    = State_incr(3,:)*scale_temp
                   cat_progn_incr(kk,:)%tc2    = State_incr(4,:)*scale_temp
                   cat_progn_incr(kk,:)%tc4    = State_incr(5,:)*scale_temp
                   cat_progn_incr(kk,:)%ght(1) = State_incr(6,:)*scale_ght1
                   
                end if
                
             end if
             
          end if     ! thresholds
          
       end do

       ! ----------------------------------       

    case (9) select_update_type   ! 1d Tskin/ght(1) analysis; FT obs
       
       if (logit) write (logunit,*) 'get 1d Tskin/ght(1) increments; FT obs'
       
       ! rule-based FT analysis similar to Farhadi et al., JHM, 2014

       ! "1d" update using obs that may or may not be in the model tile space.
       ! This approach differs from the early "1d" updates that assume obs are
       ! provided in the model tile space.
       ! The "1d" update here implies that tiles are not updated using obs 
       ! beyond the obs FOV, even if model errors are spatially correlated
       ! beyond the obs FOV.
       ! - reichle, 20 Oct 2014

       ! determine "effective" temperature for landscape-average FT 
       
       do n_e=1,N_ens
          
          ! tp must be in CELSIUS;  FT_Teff is returned in Kelvin

          call catch_calc_FT( N_catd, asnow(:,n_e), tp(1,:,n_e),          &
               tsurf_excl_snow(:,n_e), FT_state(:), FT_Teff(:,n_e ) )
          
       end do

       ! identify species ID numbers of interest
       
       N_select_varnames  = 1
       
       select_varnames(1) = 'FT'
       
       call get_select_species(                                           &
            N_select_varnames, select_varnames(1:N_select_varnames),      &
            N_obs_param, obs_param, N_select_species, select_species )
       
       

       ! determine appropriate lat/lon distances (in units of [deg])
       !   within which to look for obs (see comment above)
       
       tmp_dlon = FOV_threshold - TINY(tmp_dlon)    ! initialize
       tmp_dlat = FOV_threshold - TINY(tmp_dlat)    ! initialize
       
       ! find maximum FOV in units of [deg] across all obs params 
       
       do ii=1,N_select_species
          
          this_obs_param = obs_param(select_species(ii))
          
          if     ( trim(this_obs_param%FOV_units)=='deg' ) then
             
             tmp_dlon = max( tmp_dlon, this_obs_param%FOV )
             tmp_dlat = max( tmp_dlat, this_obs_param%FOV )
             
          elseif ( trim(this_obs_param%FOV_units)=='km'  ) then
             
             ! convert from [km] (FOV) to [deg] 
             
             call dist_km2deg(                                                 &
                  this_obs_param%FOV, N_catd, tile_coord%com_lat, r_x, r_y )
             
             tmp_dlon = max( tmp_dlon, r_x )
             tmp_dlat = max( tmp_dlat, r_y )
             
          else
             
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown FOV_units')
             
          end if
          
       end do
              

       ! loop through tiles and compute increments
       
       do kk=1,N_catd
          
          ! find observations near tile kk
          
          ! TO DO: MAKE DLAT DEPENDENT ON LAT FOR SMAP EASE GRID OBS??
          !        (For obs on the EASE grid the current approach might not work well.)
          
          tmp_minlon = tile_coord(kk)%com_lon - tmp_dlon(kk)
          tmp_maxlon = tile_coord(kk)%com_lon + tmp_dlon(kk)
          tmp_minlat = tile_coord(kk)%com_lat - tmp_dlat
          tmp_maxlat = tile_coord(kk)%com_lat + tmp_dlat
          
          ! simple approach to dateline issue (cut back to at most -180:180, -90:90)
          
          tmp_minlon = max(tmp_minlon,-180.)
          tmp_maxlon = min(tmp_maxlon, 180.)
          tmp_minlat = max(tmp_minlat, -90.)
          tmp_maxlat = min(tmp_maxlat,  90.)
          
          call get_ind_obs_lat_lon_box(                               &
               N_obs,            Observations,                        &
               tmp_minlon, tmp_maxlon, tmp_minlat, tmp_maxlat,        &
               N_select_species, select_species(1:N_select_species),  &
               N_selected_obs,   ind_obs )
          
          if (N_selected_obs > 0) then
             
             ! compute average obs value
                
             tmp_obs = sum(Observations(ind_obs(1:N_selected_obs))%obs)
             
             tmp_obs = tmp_obs/real(N_selected_obs)

             ! compute increments
             
             do n_e=1,N_ens
                
                ! TO DO: RESTRICT ANALYSIS TO SITUATIONS WHERE OBS AND MODEL 
                !        CONTRADICT EACH OTHER

                ! TO DO: COMPUTE deltaT BASED ON ENS AVG TEFF???
                !        USE ENS AVG OF ASNOW VS. THRESHOLD???

                ! determine target (landscape-average) temperature increment
                
                if     (tmp_obs <= FT_ANA_FT_THRESHOLD) then                    
                   
                   ! obs thawed                                 -->  deltaT >= 0. 
                   
                   deltaT = max( (FT_ANA_LOWERBOUND_TEFF - FT_Teff(kk,n_e)), 0.) 
                   
                elseif (asnow(kk,n_e) <= FT_ANA_LOWERBOUND_ASNOW) then
                   
                   ! obs frozen and model snow below threshold  -->  deltaT <= 0.
                   
                   deltaT = min( (FT_ANA_UPPERBOUND_TEFF - FT_Teff(kk,n_e)), 0.) 
                   
                end if
                
                ! set temperature increment for each component temperature
                ! (do nothing if deltaT=0. because cat_progn_incr was initialized to 0.)
                
                if (abs(deltaT)>0.) then

                   ! TO DO: SHOULD PHASE CHANGE BE PREVENTED FOR TC1, TC2, TC4 AS WELL?
                   !        SHOULD PHASE CHANGE BE PREVENTED BASED ON LAND COVER/CSOIL?

                   ! surface temperature increments
                   
                   cat_progn_incr(kk,n_e)%tc1 = deltaT
                   cat_progn_incr(kk,n_e)%tc2 = deltaT
                   cat_progn_incr(kk,n_e)%tc4 = deltaT
                   
                   ! soil temperature increment
                   
                   ght1_minus = cat_progn(kk,n_e)%ght(1) ! model forecast 
                   
                   fice_minus = fice(1,kk,n_e)           ! model forecast
                   
                   tp1_minus  = tp(1,kk,n_e)             ! model forecast [CELSIUS]       

                   fice_plus  = fice_minus               ! ice fraction does not change

                   tp1_plus   = tp1_minus + deltaT       ! tentative tp1 analysis [CELSIUS]
                   
                   ! avoid phase change of soil temp
                   
                   if ((tp1_minus*tp1_plus) < 0.)  tp1_plus = 0.
                   
                   ! compute ght1_plus from tp1_plus and fice_plus
                   
                   call catch_calc_ght( cat_param(kk)%dzgt(1),                    &
                        cat_param(kk)%poros, tp1_plus, fice_plus, ght1_plus )
                   
                   cat_progn_incr(kk,n_e)%ght(1) = ght1_plus - ght1_minus
                   
                end if

             end do
             
          end if
          
       end do

       ! ----------------------------------
       
    case default
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown update_type')
       
    end select select_update_type
    
    ! clean up
    
    if (allocated( State_incr ))         deallocate( State_incr )
    if (allocated( State_lon ))          deallocate( State_lon )
    if (allocated( State_lat ))          deallocate( State_lat )
    if (allocated( select_tilenum ))     deallocate( select_tilenum )

    if (allocated( Obs_cov        ))     deallocate( Obs_cov )
    
    if (associated(N_tile_in_cell_ij  )) deallocate( N_tile_in_cell_ij   )
    if (associated(tile_num_in_cell_ij)) deallocate( tile_num_in_cell_ij )
        
    ! NO checks of prognostics after update, this is now done after 
    ! increments have been applied.
    ! - reichle, 18 Oct 2005

  end subroutine cat_enkf_increments
  
  ! **********************************************************************

  subroutine get_select_species(                                       &
       N_select_varnames, select_varnames, N_obs_param, obs_param,     &
       N_select_species, select_species )
    
    ! find out obs species ID numbers ("obs_param%species") that correspond 
    !  to a set of obs species variables names ("obs_param%varname")
    ! - reichle, 16 Oct 2014

    implicit none
    
    integer,                                            intent(in)  :: N_select_varnames
    integer,                                            intent(in)  :: N_obs_param
    
    character(len=*),     dimension(N_select_varnames), intent(in)  :: select_varnames
    
    type(obs_param_type), dimension(N_obs_param),       intent(in)  :: obs_param        
    
    integer,                                            intent(out) :: N_select_species
    integer,              dimension(N_obs_param),       intent(out) :: select_species
    
    ! local variables
    
    integer :: ii, kk

    ! -----------------------------------------------------------------

    ! initialize

    select_species = -7777
    
    kk = 0

    if (N_select_varnames > 0) then
       
       do ii=1,N_obs_param
          
          if (any(trim(obs_param(ii)%varname)==select_varnames)) then
             
             kk = kk+1
             
             select_species(kk) = obs_param(ii)%species
             
          end if
          
       end do

    end if

    N_select_species = kk
    
  end subroutine get_select_species
  

  ! **********************************************************************
  
  subroutine get_ind_obs(                 &
       N_obs,            Observations,    &
       N_select_tilenum, select_tilenum,  &
       N_select_species, select_species,  &
       N_selected_obs,   ind_obs )
    
    ! find the (index vector for) observations matching a given selection
    ! of catchments or observation types
    !
    ! the vector select_tilenum of length N_select_tilenum contains the
    ! catchment numbers to be selected.
    ! if N_select_tilenum is zero, all catchments are selected
    !
    ! the vector select_species of length N_select_species contains the
    ! species IDs to be selected.
    ! if N_select_species is zero, all species are selected
    ! 
    ! the indices (relative to "Observations") of matching observations
    ! are returned in the first "N_selected_obs" components of the 
    ! vector "ind_obs" (which must be allocated to (maximum) length "N_obs"
    ! by the calling routine)
    !
    ! simplified using "any()" - reichle, 16 Oct 2014
        
    implicit none
    
    integer,        intent(in) :: N_obs, N_select_tilenum, N_select_species
    
    type(obs_type), intent(in),  dimension(N_obs)            :: Observations
    
    integer,        intent(in),  dimension(N_select_tilenum) :: select_tilenum
    integer,        intent(in),  dimension(N_select_species) :: select_species
    
    integer,        intent(out)                              :: N_selected_obs
    
    integer,        intent(out), dimension(N_obs)            :: ind_obs
    
    ! ---------------------------
    
    ! locals
    
    integer :: i, j, k, m
    
    logical :: selected_obs
    
    ! --------------------------------------------------------------
    
    if (N_select_species==0 .and. N_select_tilenum==0) then
       
       ! all observations are selected
       
       do i=1,N_obs
          ind_obs(i) = i
       end do
       
       N_selected_obs = N_obs
              
    else if (N_select_species==0) then    ! select given catchments
       
       k = 0                              ! counter for selected obs
       
       do i=1,N_obs

          if (any(Observations(i)%tilenum == select_tilenum)) then
             
             k          = k+1
             ind_obs(k) = i
             
          end if
          
       end do
       
       N_selected_obs = k

    else if (N_select_tilenum==0) then    ! select given species
       
       k = 0                              ! counter for selected obs
       
       do i=1,N_obs

          if (any(Observations(i)%species == select_species)) then

             k          = k+1
             ind_obs(k) = i

          end if
          
       end do
       
       N_selected_obs = k
       
    else                                  ! select given species and catchments
       
       k = 0                              ! counter for selected obs
       
       do i=1,N_obs

          if ( any(Observations(i)%tilenum == select_tilenum) .and.         &
               any(Observations(i)%species == select_species)       ) then
             
             k          = k+1
             ind_obs(k) = i
             
          end if
          
       end do
       
       N_selected_obs = k
       
    end if
    
  end subroutine get_ind_obs
  
  
  ! **********************************************************************
  
  subroutine get_ind_obs_assim( N_obs, assim_flag, N_obs_assim, ind_obs_assim )
    
    ! loop through Observations%assim and construct an index vector that maps to
    ! those obs that are meant to be assimilated (%assim==true)
    !
    ! reichle, 27 March 2014

    implicit none
    
    integer, intent(in)                    :: N_obs
    
    logical, intent(in),  dimension(N_obs) :: assim_flag  ! typically "Observations%assim"
    
    integer, intent(out)                   :: N_obs_assim
    integer, intent(out), dimension(N_obs) :: ind_obs_assim
    
    ! ---------------------------
    
    ! locals
    
    integer :: ii, jj
    
    ! --------------------------------------------------------------

    ind_obs_assim = -9999
    
    jj=0
    
    do ii=1,N_obs
       
       if (assim_flag(ii)) then
          
          jj=jj+1
          
          ind_obs_assim(jj) = ii
          
       end if
       
    end do
    
    N_obs_assim = jj

  end subroutine get_ind_obs_assim


  ! ********************************************************************

  function get_halo_around_tile(tile, xcompact, ycompact, skin) result(halo)

    ! determine halo around a tile for the purpose of the EnKF analysis
    ! - pchakrab, 25 March 2014

    ! input/output
    
    type(tile_coord_type), intent(in)           :: tile
    real,                  intent(in)           :: xcompact, ycompact
    real,                  intent(in), optional :: skin
    type(halo_type)                             :: halo ! output
    
    ! local
    real :: this_lon, this_lat, tmp_skin
    
    tmp_skin = 1.0
    if (present(skin)) tmp_skin = skin

    this_lon = tile%com_lon
    this_lat = tile%com_lat

    ! simple approach to dateline issue - 
    ! cut halo back to at most -180:180, -90:90
    halo%minlon = max(this_lon-tmp_skin*xcompact,-180.)
    halo%maxlon = min(this_lon+tmp_skin*xcompact, 180.)
    halo%minlat = max(this_lat-tmp_skin*ycompact, -90.)
    halo%maxlat = min(this_lat+tmp_skin*ycompact,  90.)
    
  end function get_halo_around_tile

  ! **********************************************************************
  
  function TileNnzObs(Observations, halo, select_species) result (nnz)

    ! determine whether or not any Observations (possibly restricted
    ! to a select list of species) fall within a given halo
    !
    ! - pchakrab, 25 March 2014

    ! "nnz" = non-zero
    
    implicit none
    
    ! input/output
    type(obs_type),  intent(in), dimension(:) :: Observations 
    type(halo_type), intent(in)               :: halo
    integer,         intent(in), dimension(:) :: select_species
    logical                                   :: nnz             ! output
    
    ! locals
    integer :: N_obs,   N_select_species
    integer :: iObs,    jSpecies
    real    :: lon_obs, lat_obs

    N_obs            = size(Observations)
    N_select_species = size(select_species) ! size=0 for un-allocated array
    
    nnz = .false.
    
    if (N_select_species==0) then      ! use all species
       do iObs=1,N_obs
          ! center-of-mass coordinates for the given observation
          lon_obs = Observations(iObs)%lon
          lat_obs = Observations(iObs)%lat
          if ( halo%minlon<=lon_obs .and. lon_obs<=halo%maxlon .and. &
               halo%minlat<=lat_obs .and. lat_obs<=halo%maxlat ) then
             nnz = .true.
             exit                   ! out of iObs loop
          end if
       end do                          ! end loop over observations
    else                               ! pick out selected species
       !
       ! pchakrab: THIS SECTION NEEDS TO BE TESTED
       !
       do iObs=1,N_obs
          lon_obs = Observations(iObs)%lon
          lat_obs = Observations(iObs)%lat
          do jSpecies=1,N_select_species 
             if ( halo%minlon<=lon_obs .and. lon_obs<=halo%maxlon .and. &
                  halo%minlat<=lat_obs .and. lat_obs<=halo%maxlat .and. &
                  (Observations(iObs)%species == select_species(jSpecies)) ) then  
                nnz = .true.
                exit                   ! exit to next observation (next iObs)
             end if
          end do                       ! end loop over select_species
          if (nnz) exit                ! out of iObs loop
       end do
    end if
    
  end function TileNnzObs


  ! **********************************************************************

  subroutine get_ind_obs_lat_lon_box(         &
       N_obs,            Observations,        &
       min_lon, max_lon, min_lat, max_lat,    &
       N_select_species, select_species,      &
       N_selected_obs,   ind_obs )
    
    ! find the (index vector for) observations within a given lat/lon box
    ! and for given observation types
    !
    ! min_lon, max_lon, min_lat, max_lat describe the extent of the lat/lon box
    !
    ! the vector select_species of length N_select_species contains the
    ! species IDs to be selected.
    ! if N_select_species is zero, all species are selected
    ! 
    ! the indices (relative to "Observations") of matching observations
    ! are returned in the first "N_selected_obs" components of the 
    ! vector "ind_obs" (which must be allocated to (maximum) length "N_obs"
    ! by the calling routine)
    !
    ! NOTE: 1.) definitions:
    !            update region = geographic region for which the state
    !                             vector is updated
    !            lat-lon-box   = geographic region of observations that 
    !                             influence the states of the udpate region
    !                             (encloses update region, bigger than
    !                              update region by xcompact/ycompact)
    !       2.) update regions must NOT cross the date line
    !       3.) observations that fall into the lat-lon-box but are
    !            on the opposite side of the date line from the
    !            update region are NOT used (of course they will be used
    !            to update their "native" update region)
    !    
    ! uses nr_indexx.f
    !
    ! reichle,  26 Jul 2002 
    ! reichle,   1 Aug 2005 - use Observations%lat/lon instead of tile_coord
    ! reichle,  21 Jun 2013 - added sort to avoid lay-out dependency for parallel runs
    ! pchakrab, 25 Mar 2014 - sort is not needed following fix in read_obs.F90
    ! reichle,  16 Oct 2014 - simplified using "any()"

    implicit none
    
    integer, intent(in) :: N_obs, N_select_species
    
    type(obs_type), intent(in), dimension(N_obs) :: Observations 
    
    real :: min_lon, max_lon, min_lat, max_lat
    
    integer, intent(in), dimension(N_select_species)  :: select_species

    integer, intent(out) :: N_selected_obs
    
    integer, intent(out), dimension(N_obs) :: ind_obs
    
    ! ---------------------------
    
    ! locals
    
    integer :: i, j, k

    real :: lon_obs, lat_obs

    ! --------------------------------------------------------------
               
    k = 0                              ! counter for selected obs
    
    if (N_select_species==0) then      ! use all species
       
       do i=1,N_obs
          
          ! determine center-of-mass coordinates for the given observation
          
          lon_obs = Observations(i)%lon
          lat_obs = Observations(i)%lat
          
          if ( min_lon <= lon_obs  .and.          &
               lon_obs <= max_lon  .and.          &
               min_lat <= lat_obs  .and.          &
               lat_obs <= max_lat         ) then
             
             k          = k+1
             ind_obs(k) = i
             
          end if
          
       end do                          ! end loop through observations
       
       N_selected_obs = k
       
    else                               ! pick out selected species
       
       do i=1,N_obs
          
          ! determine center-of-mass coordinates for the given observation
          
          lon_obs = Observations(i)%lon
          lat_obs = Observations(i)%lat
          
          if ( any(Observations(i)%species == select_species)  .and.         &
               min_lon <= lon_obs                              .and.         &
               lon_obs <= max_lon                              .and.         &
               min_lat <= lat_obs                              .and.         &
               lat_obs <= max_lat                                    ) then
             
             k          = k+1
             ind_obs(k) = i
             
          end if
          
       end do
       
       N_selected_obs = k
       
    end if
        
  end subroutine get_ind_obs_lat_lon_box

  ! *********************************************************************

  subroutine check_compact_support(                                &
       update_type,                                                &
       N_catf,             tile_coord_f,                           &
       N_progn_pert,       progn_pert_param,                       &
       N_force_pert,       force_pert_param,                       &
       N_obs_param,        obs_param,                              &
       xcompact, ycompact )
    
    ! check whether any of the correlation scales exceeds or is comparable to
    !  the compact support length scale
    ! also check whether the compact support length scale may be 
    !  too large for the given correlation scales
    
    implicit none
    
    integer, intent(in) :: update_type

    integer, intent(in) :: N_catf, N_progn_pert, N_force_pert, N_obs_param

    type(tile_coord_type), dimension(N_catf),       intent(in) :: tile_coord_f

    type(pert_param_type), dimension(N_progn_pert), intent(in) :: progn_pert_param
    type(pert_param_type), dimension(N_force_pert), intent(in) :: force_pert_param
        
    type(obs_param_type),  dimension(N_obs_param),  intent(in) :: obs_param
    
    real, intent(inout) :: xcompact, ycompact
    
    ! ---------------------------------
    
    ! locals
    
    integer       :: ii
    
    character(40) :: error_type

    real          :: max_dist_x, max_dist_y, r_y
    
    real, dimension(1) :: r_x, tmp_lat

    character(len=*), parameter :: Iam = 'check_compact_support'

    ! -------------------------------------------------------------
    
    select case (update_type)

    case (1,3,4,5,6,9)   ! "1d" updates

       ! Make xcompact and ycompact just large enough so that 
       ! the EnKF analysis correctly identifies the tiles 
       ! that should be included in the 1d analysis.
       !
       ! In the obs readers, each obs is assigned to a single tile
       ! based on subroutine get_tile_num_from_latlon().  The assignment
       ! is such that the the obs lat/lon and the tile center-of-mass lat/lon
       ! must always be within the same grid cell of the tile_grid [or pert_grid??].
       ! Furthermore, the obs lat/lon 
       !  - must be within the min/max lat/lon boundaries of the tile
       !  OR
       !  - must be within max_dist (defined by obs_param%FOV)
       
       ! Based on the above, compute the maximum distance between the lat/lon 
       ! associated with an obs and the center-of-mass lat/lon for the tile
       ! to which the obs is assigned
       
       max_dist_x = maxval( tile_coord_f%max_lon - tile_coord_f%min_lon )
       max_dist_y = maxval( tile_coord_f%max_lat - tile_coord_f%min_lat )
       
       do ii=1,N_obs_param
          
          if (obs_param(ii)%assim) then
             
             if     ( trim(obs_param(ii)%FOV_units)=='deg' ) then
                
                r_x(1) = obs_param(ii)%FOV                
                r_y    = obs_param(ii)%FOV
                
             elseif ( trim(obs_param(ii)%FOV_units)=='km'  ) then
                
                ! convert from [km] (FOV) to [deg] (max_dist_*)
                
                ! largest FOV in [deg] will be at largest abs(lat)
                
                tmp_lat(1) = maxval(abs(tile_coord_f%com_lat))

                call dist_km2deg( obs_param(ii)%FOV, 1, tmp_lat, r_x, r_y )
                
             else
                
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown FOV_units')
                
             end if
             
             max_dist_x = max( max_dist_x, r_x(1) )
             max_dist_y = max( max_dist_y, r_y    )
             
          end if
          
       end do
       
       ! reset xcompact and ycompact accordingly
       
       xcompact = max( xcompact, max_dist_x )
       ycompact = max( ycompact, max_dist_y )

       if (logit) write (logunit,*)
       if (logit) write (logunit,'(A,ES10.3)')                                         &
            Iam // '(): reset for 1d update_type: xcompact = ', xcompact
       if (logit) write (logunit,'(A,ES10.3)')                                         &
            Iam // '(): reset for 1d update_type: ycompact = ', ycompact
       if (logit) write (logunit,*)
       
    case (2,7,8,10)  ! "3d" updates, check consistency of xcompact, ycompact
       
       ! check xcompact/ycompact against corr scales of model error
       
       do ii=1,N_progn_pert
          
          error_type = progn_pert_param(ii)%descr
          
          call check_compact( xcompact, ycompact,                      &
               progn_pert_param(ii)%xcorr, progn_pert_param(ii)%ycorr, &
               error_type )
          
       end do
              
       ! check xcompact/ycompact against corr scales of forcing perturbations
       
       do ii=1,N_force_pert
          
          error_type = force_pert_param(ii)%descr
          
          call check_compact( xcompact, ycompact,                      &
               force_pert_param(ii)%xcorr, force_pert_param(ii)%ycorr, &
               error_type )
          
       end do
       
       ! check xcompact/ycompact against corr scales of observation errors
       ! (only if obs are assimilated)
       
       do ii=1,N_obs_param
          
          if (obs_param(ii)%assim) then
             
             error_type = obs_param(ii)%descr
             
             call check_compact( xcompact, ycompact,                   &
                  obs_param(ii)%xcorr, obs_param(ii)%ycorr,            &
                  error_type )
             
          end if
          
       end do
       
    case default
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown update_type')
       
    end select
    
  end subroutine check_compact_support
  
  ! ******************************************************************
  
  subroutine check_compact( xcompact, ycompact, xcorr, ycorr, error_type )
    
    ! Helper routine for check_compact_support()
    !
    ! If xcompact/ycompact are less than the largest correlation length times
    !  a multiple, the execution of the program is halted.
    !
    ! If xcompact/ycompact exceed the smallest correlation length times
    !  a multiple, a warning statement is printed.
        
    implicit none
    
    real, intent(in) :: xcompact, ycompact
    real, intent(in) :: xcorr,    ycorr
    
    character(*), intent(in) :: error_type
    
    real, parameter :: min_compact_div_corr = 2.
    real, parameter :: max_compact_div_corr = 5.
        
    character(len=*), parameter :: Iam = 'check_compact'
    character(len=400) :: err_msg

    ! -------------------------------------------------------
    
    ! check whether compact support length scale might be too large
    
    if ( (xcompact > xcorr*max_compact_div_corr) .or. &
         (ycompact > ycorr*max_compact_div_corr)          ) then
       
       if (logit) then
          write (logunit,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write (logunit,*)
          write (logunit,*) 'WARNING: xcompact/ycompact may be too large for'
          write (logunit,*) '         error corr scale of ', trim(error_type)
          write (logunit,*) 
          write (logunit,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write (logunit,*)
       end if
       
    end if
    
    ! check whether compact support length scale is too small
    
    if ( (xcompact < xcorr*min_compact_div_corr) .or.      &
         (ycompact < ycorr*min_compact_div_corr)         ) then
       err_msg = 'xcompact/ycompact too small for ' // &
            'error corr scale of ' // trim(error_type)
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

  end subroutine check_compact

  ! ******************************************************************
  
  subroutine dist_km2deg( dist_km, N_lat, lat, dist_x_deg, dist_y_deg )
    
    ! Convert (horizontal) distance from units of [km] to units of [deg]
    !
    ! Output meridional distance in units of [deg] depends on latitude
    
    implicit none
    
    real,                      intent(in)  :: dist_km     ! [km]
    
    integer,                   intent(in)  :: N_lat
    
    real,    dimension(N_lat), intent(in)  :: lat         ! [deg] (-90:90)
    
    real,    dimension(N_lat), intent(out) :: dist_x_deg  ! [deg] vector (depends on latitude)
    real,                      intent(out) :: dist_y_deg  ! [deg] scalar

    ! local variables
    
    character(len=*),  parameter :: Iam     = 'dist_km2deg'
    
    ! -------------------------------------------------------
    
    ! NOTE: MAPL_radius (Earth radius) is in [m] and dist_km is in [km]
    
    dist_y_deg = dist_km * (180./MAPL_PI) / (MAPL_RADIUS/1000.)
    
    ! NOTE: cos() needs argument in [rad], lat is in [deg] (-90:90)
    
    dist_x_deg = dist_y_deg / cos( MAPL_PI/180. * lat )
    
    if (any(dist_x_deg<0.) .or. dist_y_deg<0.)  &
         call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'encountered negative distance' )
    
  end subroutine dist_km2deg
  
  ! ******************************************************************

! abandoned during changes for obs halo
! to reinstate, need to MPI_Gather "cat_param_f"
!
! reichle, 30 Sep 2011
  
!  subroutine check_obs_pert( N_ens, N_catd, N_obs, &
!       cat_param, Observations, Obs_pert )
!    
!    ! check synthetic observation error for physical constraints, re-set obs 
!    !  error if constraints are violated
!    !
!    ! reichle, 26 Feb 2002
!    ! reichle,  2 Aug 2005
!    ! reichle, 10 Jan 2011 - added checks for AMSR-E/LPRM soil moisture retrievals
!    !
!    ! -------------------------------------------------------------------
!    
!    implicit none
!    
!    integer, intent(in) :: N_ens, N_obs, N_catd
!    
!    type(cat_param_type), dimension(N_catd), intent(in) :: cat_param
!    
!    type(obs_type), dimension(N_obs), intent(in) :: Observations   
!    
!    real, dimension(N_obs,N_ens), intent(inout)  :: Obs_pert
!    
!    ! ----------------------------------------------------------------
!    
!    ! local variables
!    
!    integer :: i, n
!    
!    real :: min_pert, max_pert
!    
!    ! ---------------------------------------------------------------
!    
!    do i=1,N_obs
!       
!       select case (Observations(i)%species)
!          
!       case (1,2,4,7,8,9,10,11,12)
!
!          ! ae_l2_sm_a, ae_l2_sm_d, RedArkOSSE_sm, RedArkOSSE_CLSMsynthSM,
!          ! VivianaOK_CLSMsynthSM, ae_sm_LPRM_a_C, ae_sm_LPRM_d_C,
!          ! ae_sm_LPRM_a_X, ae_sm_LPRM_d_X
!          
!          min_pert = -Observations(i)%obs
!          max_pert = cat_param(Observations(i)%tilenum)%poros - &
!               Observations(i)%obs
!          
!          do n=1,N_ens
!             
!             Obs_pert(i,n) = max( Obs_pert(i,n), min_pert ) 
!             Obs_pert(i,n) = min( Obs_pert(i,n), max_pert )
!             
!          end do
!          
!       case (3)                     ! isccp_tskin_gswp2_v1
!          
!          ! no constraints
!          
!       case default
!          
!          call stop_it('check_obs_pert(): unkown obs species.')
!
!          !write (logunit,*) 'check_obs_pert(): unkown obs species. STOPPING.'
!          !stop
!          
!       end select
!       
!    end do
!    
!  end subroutine check_obs_pert
  
  ! **********************************************************************
  ! **********************************************************************
  ! **********************************************************************
  
end module clsm_ensupd_upd_routines


! **** EOF ******************************************************

