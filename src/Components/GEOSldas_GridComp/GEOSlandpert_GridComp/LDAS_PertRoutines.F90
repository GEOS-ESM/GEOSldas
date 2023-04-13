!
! this file contains a collection of subroutines that are needed to
! run the Ensemble Kalman filter with the catchment model off-line driver
!
! reichle, 10 May 2005
! reichle,  6 Dec 2013 - introduced "progn_pert_type"
!                        (no longer use "cat_progn_type" for progn perts)
! reichle, 21 Nov 2014 - re-interpreted progn_pert as perturbation flux forcing
!                      - added "qair" and "wind" perts in apply_force_pert()
!                      - renamed force_pert_type fields for consistency w/ met_force_type
!                          %tmp2m --> %tair  (but note lower-case!)
!                          %dpt2m --> %qair  (but note lower-case!)
!                          %wnd   --> %wind  (but note lower-case!)
#include "MAPL_Generic.h"

module LDAS_PertRoutinesMod

  use ESMF
  use MAPL_Mod

  use LDAS_ensdrv_Globals,              ONLY:     &
       logunit,                                   &
       root_logit,                                &
       nodata_generic,                            &
       nodata_tolfrac_generic,                    &
       nodata_tol_generic
  use LDAS_ensdrv_mpi, only: mpicomm,numprocs,myid
  use MAPL_ConstantsMod,                ONLY:     &
       Tzero => MAPL_TICE,                        &
       alhe  => MAPL_ALHL,                        &
       alhs  => MAPL_ALHS

  use LDAS_TileCoordType,               ONLY:     &
       tile_coord_type,                           &
       grid_def_type,                             &
       io_grid_def_type

  use LDAS_TileCoordRoutines,           ONLY:     &
        LDAS_create_grid_g

  use LDAS_PertTypes,                   ONLY:     &
       pert_param_type,                           &
       allocate_pert_param

  ! CHANGED: cat_param/progn_types are only used by apply_*_pert
  ! and these routine (apply_*_pert) are no longer part of this module
  ! use catch_types,                      ONLY:     &
  !      cat_param_type,                            &
  !      cat_progn_type
  !
  !use LDAS_DriverTypes,                 ONLY:     &
  !     met_force_type

  use force_and_cat_progn_pert_types,   ONLY:     &
       N_force_pert_max,                          &
       force_pert_real_type,                      &
       force_pert_logi_type,                      &
       force_pert_char_type,                      &
       force_pert_ccor_type,                      &
       N_progn_pert_max,                          &
       progn_pert_real_type,                      &
       progn_pert_logi_type,                      &
       progn_pert_char_type,                      &
       progn_pert_ccor_type,                      &
       struct2vec_force_pert,                     &
       struct2mat_force_pert_ccor,                &
       struct2vec_progn_pert,                     &
       struct2mat_progn_pert_ccor,                &
       assignment (=)

  use LDAS_DateTimeMod,                 ONLY:     &
       date_time_type,                            &
       datetime2_minus_datetime1

  use nr_ran2_gasdev,                   ONLY:     &
       NRANDSEED

  use land_pert_routines,               ONLY:     &
       get_pert,                                  &
       get_sqrt_corr_matrix,                      &
       get_init_Pert_rseed

  ! CHANGED: We are going to use MAPL_LocStreamXform instead
  ! use tile_coord_routines,              ONLY:     &
  !      grid2tile,                                 &
  !      get_is_land

  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename

  ! TODO: maybe copy over check_cat_progn to some file here
  ! use clsm_ensdrv_drv_routines,         ONLY:     &
  !      check_cat_progn

  use RepairForcingMod,                 ONLY:     &
       repair_forcing

  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use netcdf

  implicit none

  include 'mpif.h'

  ! everything is private by default unless made public

  private

  public :: read_ens_prop_inputs
  public :: interpolate_pert_to_timestep
  public :: get_pert_grid
  public :: get_progn_pert_param
  public :: get_force_pert_param
  public :: echo_pert_param
  ! WY note :: io_pert_rstrt() was adapted. read from LDASsa and write to a nc4 file as MAPL internal
  public :: io_pert_rstrt
  public :: check_pert_dtstep
  ! ADDED
  public :: apply_pert
  ! the parameters below will be overwritten by RC file 
  integer,public :: GEOSldas_NUM_ENSEMBLE = -1
  integer,public :: GEOSldas_FIRST_ENS_ID = -1
  integer,public :: GEOSldas_FORCE_PERT_DTSTEP = -1
  integer,public :: GEOSldas_PROGN_PERT_DTSTEP = -1

contains

  ! *******************************************************************

  subroutine read_ens_prop_inputs(   &
       write_nml,                    &
       work_path,                    &
       exp_id,                       &
       date_time,                    &
       kw_echo,                      &
       kw_N_ens,                     &
       kw_ens_id,                    &
       kw_progn_pert_dtstep,         &
       kw_force_pert_dtstep,         &
       kw_descr_progn_pert,          &
       kw_typ_progn_pert,            &
       kw_std_progn_pert,            &
       kw_stdfromfile_progn_pert,    &
       kw_stdfilename_progn_pert,    &
       kw_zeromean_progn_pert,       &
       kw_coarsen_progn_pert,        &
       kw_std_normal_max_progn_pert, &
       kw_xcorr_progn_pert,          &
       kw_ycorr_progn_pert,          &
       kw_tcorr_progn_pert,          &
       kw_ccorr_progn_pert,          &
       kw_descr_force_pert,          &
       kw_typ_force_pert,            &
       kw_std_force_pert,            &
       kw_stdfromfile_force_pert,    &
       kw_stdfilename_force_pert,    &
       kw_zeromean_force_pert,       &
       kw_coarsen_force_pert,        &
       kw_std_normal_max_force_pert, &
       kw_xcorr_force_pert,          &
       kw_ycorr_force_pert,          &
       kw_tcorr_force_pert,          &
       kw_ccorr_force_pert           &
       )

    ! read ensemble propagation inputs from namelist file
    !
    ! runtime options are read in three steps:
    !
    ! 1.) read options from default namelist file called
    !      ens_prop_inputs.nml in working directory (must be present)
    !
    ! pchakrab: Removing the 2nd and 3rd options. Options will be read
    !           ONLY from the 'default' namelist file

    ! 2.) overwrite options from special namelist file (if present)
    !      specified at the command line using -ens_prop_inputs_path
    !      and -ens_prop_inputs_file
    !
    ! reichle, 29 Mar 2004
    ! reichle, 31 Aug 2004 - added tskin_isccp
    ! reichle, 31 May 2005 - redesign for CLSM ens driver

    implicit none

    logical,              intent(in),  optional :: write_nml

    character(*),       intent(in),  optional :: work_path
    character(*),        intent(in),  optional :: exp_id

    type(date_time_type), intent(in),  optional :: date_time

    logical,              intent(in),  optional :: kw_echo

    integer,              intent(out), optional :: kw_N_ens

    integer, dimension(:), pointer,    optional :: kw_ens_id     ! output

    integer,              intent(out), optional :: kw_progn_pert_dtstep
    integer,              intent(out), optional :: kw_force_pert_dtstep

    character(*),       intent(out), optional :: kw_stdfilename_progn_pert
    character(*),       intent(out), optional :: kw_stdfilename_force_pert

    type(progn_pert_char_type), intent(out), optional :: &
         kw_descr_progn_pert

    type(progn_pert_logi_type), intent(out), optional :: &
         kw_zeromean_progn_pert,                         &
         kw_coarsen_progn_pert,                          &
         kw_stdfromfile_progn_pert

    type(progn_pert_real_type), intent(out), optional :: &
         kw_std_normal_max_progn_pert,                   &
         kw_std_progn_pert,                              &
         kw_xcorr_progn_pert,                            &
         kw_ycorr_progn_pert,                            &
         kw_tcorr_progn_pert,                            &
         kw_typ_progn_pert

    type(progn_pert_ccor_type), intent(out), optional :: &
         kw_ccorr_progn_pert

    type(force_pert_char_type), intent(out), optional :: &
         kw_descr_force_pert

    type(force_pert_logi_type), intent(out), optional :: &
         kw_zeromean_force_pert,                         &
         kw_coarsen_force_pert,                          &
         kw_stdfromfile_force_pert

    type(force_pert_real_type), intent(out), optional :: &
         kw_std_normal_max_force_pert,                   &
         kw_std_force_pert,                              &
         kw_xcorr_force_pert,                            &
         kw_ycorr_force_pert,                            &
         kw_tcorr_force_pert,                            &
         kw_typ_force_pert

    type(force_pert_ccor_type), intent(out), optional :: &
         kw_ccorr_force_pert

    ! ------------------------

    ! locals

    character(300)  :: fname

    character(200)  :: ens_prop_inputs_path
    character( 40)  :: ens_prop_inputs_file, dir_name, file_tag, file_ext

    integer :: i, N_ens, first_ens_id, progn_pert_dtstep, force_pert_dtstep

    character(300)  :: stdfilename_progn_pert, stdfilename_force_pert

    type(progn_pert_char_type) :: descr_progn_pert

    type(progn_pert_logi_type) :: zeromean_progn_pert
    type(progn_pert_logi_type) :: coarsen_progn_pert
    type(progn_pert_logi_type) :: stdfromfile_progn_pert

    type(progn_pert_real_type) :: std_normal_max_progn_pert
    type(progn_pert_real_type) :: std_progn_pert
    type(progn_pert_real_type) :: xcorr_progn_pert
    type(progn_pert_real_type) :: ycorr_progn_pert
    type(progn_pert_real_type) :: tcorr_progn_pert
    type(progn_pert_real_type) :: typ_progn_pert

    type(progn_pert_ccor_type) :: ccorr_progn_pert

    type(force_pert_char_type) :: descr_force_pert

    type(force_pert_logi_type) :: zeromean_force_pert
    type(force_pert_logi_type) :: coarsen_force_pert
    type(force_pert_logi_type) :: stdfromfile_force_pert

    type(force_pert_real_type) :: std_normal_max_force_pert
    type(force_pert_real_type) :: std_force_pert
    type(force_pert_real_type) :: xcorr_force_pert
    type(force_pert_real_type) :: ycorr_force_pert
    type(force_pert_real_type) :: tcorr_force_pert
    type(force_pert_real_type) :: typ_force_pert

    type(force_pert_ccor_type) :: ccorr_force_pert

    ! Errlong variables
    integer :: rc, status
    character(len=*), parameter :: Iam = 'read_ens_prop_inputs'

    ! MPI variables
    logical :: root_proc,f_exist

    ! -----------------------------------------------------------------

    namelist /ens_prop_inputs/       &
         N_ens,                      &
         first_ens_id,               &
         progn_pert_dtstep,          &
         force_pert_dtstep,          &
         descr_progn_pert,           &
         typ_progn_pert,             &
         std_progn_pert,             &
         zeromean_progn_pert,        &
         coarsen_progn_pert,         &
         stdfromfile_progn_pert,     &
         stdfilename_progn_pert,     &
         std_normal_max_progn_pert,  &
         xcorr_progn_pert,           &
         ycorr_progn_pert,           &
         tcorr_progn_pert,           &
         ccorr_progn_pert,           &
         descr_force_pert,           &
         typ_force_pert,             &
         std_force_pert,             &
         zeromean_force_pert,        &
         coarsen_force_pert,         &
         stdfromfile_force_pert,     &
         stdfilename_force_pert,     &
         std_normal_max_force_pert,  &
         xcorr_force_pert,           &
         ycorr_force_pert,           &
         tcorr_force_pert,           &
         ccorr_force_pert


    root_proc = (myid ==0)

    ! ---------------------------------------------------------------------
    !
    ! initialize selected name list inputs
    ! (useful if not all fields of a structure are set explicitly
    !  in namelist file)

    ccorr_progn_pert = nodata_generic
    ccorr_force_pert = nodata_generic

    ! ------------------------------------------------------
    !
    ! Set default file name for ens prop inputs namelist file

    ens_prop_inputs_path = './'                                       ! set default
    ens_prop_inputs_file = 'LDASsa_DEFAULT_inputs_ensprop.nml'

    ! Read data from default ens_prop_inputs namelist file

    fname = trim(ens_prop_inputs_path) // '/' // trim(ens_prop_inputs_file)
    open (10, file=fname, delim='apostrophe', action='read', status='old')

    if (present(kw_echo)) then
       if (kw_echo) then

          if(root_logit) then
             write (logunit,*)
             write (logunit,'(400A)') 'reading *default* ens prop inputs from ' // trim(fname)
             write (logunit,*)
          endif
       end if
    end if

    read (10, nml=ens_prop_inputs)
    
    close(10,status='keep')


    fname = './LDASsa_SPECIAL_inputs_ensprop.nml'
    inquire(file=fname,exist=f_exist)

    if (f_exist) then
       open (10, file=fname, delim='apostrophe', action='read', status='old')
       if (present(kw_echo)) then
          if (kw_echo) then

             if(root_logit) then
                write (logunit,*)
                write (logunit,'(400A)') 'reading *SPECIAL* ens prop inputs from ' // trim(fname)
                write (logunit,*)
             endif
          end if
       end if
       read (10, nml=ens_prop_inputs)
       close(10,status='keep')
    endif
    if(     GEOSldas_NUM_ENSEMBLE      == -1 .or. GEOSldas_FIRST_ENS_ID      == -1         &
       .or. GEOSldas_FORCE_PERT_DTSTEP == -1 .or. GEOSldas_PROGN_PERT_DTSTEP == -1 ) then
       stop " GEOSldas_NUM_ENSEMBLE etc. should be initialized"
    endif
    N_ens             = GEOSldas_NUM_ENSEMBLE
    first_ens_id      = GEOSldas_FIRST_ENS_ID
    force_pert_dtstep = GEOSldas_FORCE_PERT_DTSTEP 
    progn_pert_dtstep = GEOSldas_PROGN_PERT_DTSTEP 

    ! echo variables of ens_prop_inputs

    if (present(kw_echo) .and. root_logit) then
       if (kw_echo) then

          write (logunit,*) 'ens_prop inputs are:'
          write (logunit,*)
          write (logunit, nml=ens_prop_inputs)
          write (logunit,*)

       end if
    end if

    ! -------------------------------------------------------------

    if (present(kw_N_ens)) then

       kw_N_ens            = N_ens

       if (present(kw_ens_id)) then

          allocate(kw_ens_id(N_ens))
          do i=1,N_ens
             kw_ens_id(i) = first_ens_id + i - 1
          end do
          if(root_logit) then
             write (logunit,*)
             write (logunit,*) 'ens_id = ', (kw_ens_id(i), i=1,N_ens)
             write (logunit,*)
          endif
       end if
    end if


    ! perturbations time steps

    if (present(kw_progn_pert_dtstep)) &
         kw_progn_pert_dtstep  = progn_pert_dtstep

    if (present(kw_force_pert_dtstep)) &
         kw_force_pert_dtstep  = force_pert_dtstep


    ! other perturbations parameters

    if (present(kw_descr_progn_pert)) &
         kw_descr_progn_pert = descr_progn_pert

    if (present(kw_zeromean_progn_pert)) &
         kw_zeromean_progn_pert = zeromean_progn_pert

    if (present(kw_coarsen_progn_pert)) &
         kw_coarsen_progn_pert = coarsen_progn_pert

    if (present(kw_stdfromfile_progn_pert)) &
         kw_stdfromfile_progn_pert = stdfromfile_progn_pert

    if (present(kw_stdfilename_progn_pert)) &
         kw_stdfilename_progn_pert = stdfilename_progn_pert

    if (present(kw_std_normal_max_progn_pert)) &
         kw_std_normal_max_progn_pert = std_normal_max_progn_pert

    if (present(kw_std_progn_pert  )) kw_std_progn_pert   = std_progn_pert
    if (present(kw_xcorr_progn_pert)) kw_xcorr_progn_pert = xcorr_progn_pert
    if (present(kw_ycorr_progn_pert)) kw_ycorr_progn_pert = ycorr_progn_pert
    if (present(kw_tcorr_progn_pert)) kw_tcorr_progn_pert = tcorr_progn_pert
    if (present(kw_typ_progn_pert  )) kw_typ_progn_pert   = typ_progn_pert

    if (present(kw_ccorr_progn_pert)) kw_ccorr_progn_pert = ccorr_progn_pert

    if (present(kw_descr_force_pert)) &
         kw_descr_force_pert = descr_force_pert

    if (present(kw_zeromean_force_pert)) &
         kw_zeromean_force_pert = zeromean_force_pert

    if (present(kw_coarsen_force_pert)) &
         kw_coarsen_force_pert = coarsen_force_pert

    if (present(kw_stdfromfile_force_pert)) &
         kw_stdfromfile_force_pert = stdfromfile_force_pert

    if (present(kw_stdfilename_force_pert)) &
         kw_stdfilename_force_pert = stdfilename_force_pert

    if (present(kw_std_normal_max_force_pert)) &
         kw_std_normal_max_force_pert = std_normal_max_force_pert

    if (present(kw_std_force_pert  )) kw_std_force_pert   = std_force_pert
    if (present(kw_xcorr_force_pert)) kw_xcorr_force_pert = xcorr_force_pert
    if (present(kw_ycorr_force_pert)) kw_ycorr_force_pert = ycorr_force_pert
    if (present(kw_tcorr_force_pert)) kw_tcorr_force_pert = tcorr_force_pert
    if (present(kw_typ_force_pert  )) kw_typ_force_pert   = typ_force_pert

    if (present(kw_ccorr_force_pert)) kw_ccorr_force_pert = ccorr_force_pert

    ! ------------------------------------------------------------------
    !
    ! save driver inputs into *ens_prop_inputs.nml file

    if (present(write_nml)) then

       if (write_nml) then

          dir_name = 'rc_out'
          file_tag = 'ldas_ensprop_inputs'
          file_ext = '.nml'

          fname = get_io_filename( work_path, exp_id, file_tag, date_time=date_time, &
               dir_name=dir_name, file_ext=file_ext )

          if(root_logit) write (logunit,'(400A)') 'writing ens prop inputs to ' // trim(fname)
          if(root_logit) write (logunit,*)

          open(10, file=fname, status='unknown', action='write', &
               delim='apostrophe')

          write(10, nml=ens_prop_inputs)

          close(10, status='keep')

       end if

    end if

  end subroutine read_ens_prop_inputs

  ! *********************************************************************

  subroutine check_pert_dtstep( model_dtstep, start_time, end_time,       &
       N_progn_pert, N_force_pert, progn_pert_dtstep, force_pert_dtstep )

    ! reichle, 28 May 2013

    ! all time steps are in *seconds*

    implicit none

    integer,               intent(in) :: model_dtstep

    type(date_time_type),  intent(in) :: start_time, end_time

    integer,               intent(in) :: N_progn_pert, N_force_pert

    integer,               intent(in) :: progn_pert_dtstep, force_pert_dtstep

    ! local
    character(len=*), parameter :: Iam = 'check_pert_dtstep'
    character(len=400) :: err_msg

    ! -----------------------------------------------------------------------

    if (N_progn_pert>0) then

       if      (progn_pert_dtstep<=0)    then
          err_msg = 'progn_pert time step must be greater than 0'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       else if (progn_pert_dtstep>86400) then
          err_msg = 'progn_pert time step > 1 day not allowed'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       if      (mod(progn_pert_dtstep,model_dtstep)/=0) then
          err_msg = 'progn_pert and model time steps incompatible'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       else if (mod(86400,progn_pert_dtstep)/=0)        then
          err_msg = 'day not evenly divided by progn_pert time step'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       ! The following checks will eliminate some flexibility in restarting
       ! runs "inbetween" prognostics perturbations intervals in order to maintain
       ! reproducibility of longer runs regardless of the number and times of restarts.
       ! The longest time step in the system dictates the minimum restart interval.
       !
       ! Example: If progn_pert_dtstep=10800 seconds, ie, 3 hours, then runs
       !          can only be restarted (and end) at 0z, 3z, 6z, ...

       if  (mod(start_time%hour*3600+start_time%min*60+start_time%sec,          &
            progn_pert_dtstep)/=0) then
          err_msg = 'progn_pert time step clashes with start_time'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       if  (mod(end_time%hour*3600  +end_time%min*60  +end_time%sec,            &
            progn_pert_dtstep)/=0) then
          err_msg = 'Error: progn_pert time step clashes with end_time'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
    end if

    if (N_force_pert>0) then

       if      (force_pert_dtstep<=0)    then
          err_msg = 'force_pert time step must be greater than 0'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       else if (force_pert_dtstep>86400) then
          err_msg = 'force_pert time step > 1 day not allowed'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       if      (mod(force_pert_dtstep,model_dtstep)/=0) then
          err_msg = 'force_pert and model time steps incompatible'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       else if (mod(86400,force_pert_dtstep)/=0)        then
          err_msg = 'day not evenly divided by force_pert time step'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       ! The following checks will eliminate some flexibility in restarting
       ! runs "inbetween" forcing perturbations intervals in order to maintain
       ! reproducibility of longer runs regardless of the number and times of restarts.
       ! The longest time step in the system dictates the minimum restart interval.
       !
       ! Example: If force_pert_dtstep=10800 seconds, ie, 3 hours, then runs
       !          can only be restarted (and end) at 0z, 3z, 6z, ...

       if  (mod(start_time%hour*3600+start_time%min*60+start_time%sec,          &
            force_pert_dtstep)/=0) then
          err_msg = 'force_pert time step clashes with start_time'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       if  (mod(end_time%hour*3600  +end_time%min*60  +end_time%sec,            &
            force_pert_dtstep)/=0) then
          err_msg = 'force_pert time step clashes with end_time'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
    end if

  end subroutine check_pert_dtstep

  ! *********************************************************************
  function get_pert_grid( tile_grid) result (pert_grid)

    ! reichle, 20 May 2010
    ! jiang,   03/10/2017
    implicit none

    type(grid_def_type), intent(in)  :: tile_grid
    type(grid_def_type)              :: pert_grid
    character(len=30) :: latlon_gridname
    character(len=6) ::  lattmp,lontmp

    type(grid_def_type)  :: latlon_grid_tmp

    integer :: n_x,i_off,j_off,n_lon,n_lat

    ! in future implement such that a coarser grid could be used
    ! 
    ! NOTE: must then also modify grid2tile that is used to diagnose
    !       perturbations in tile space from gridded perturbations fields
    !       (see calls to "grid2tile" in clsm_ensdrv_pert_routines.F90,  
    !        clsm_ensupd_upd_routines.F90, and clsm_adapt_routines.F90)
    
    if(index(tile_grid%gridtype,"c3") ==0) then

       ! If *not* cube-sphere tile space, then for perturbations use the grid that 
       ! defines the tile space (a.k.a. "tile_grid").  E.g., if in EASE grid tile space, 
       ! the pert grid is the EASE grid.

       pert_grid = tile_grid

    else ! cubed-sphere grid
    
       ! For cubed-sphere tile space, use a global lat_lon pert grid with a resolution
       ! similar to that of the grid that defines the tile space.
       
       N_x=tile_grid%n_lon
       
       ! NOTE: The pert grid specification is hard-wired here.
       ! If perturbation stddev is heterogeneous input from a file,  
       ! then the input grid must match this hard-wired grid.  (sqz 2/2023)  
       
       n_lon=4*N_x
       n_lat=3*N_x
       write(lattmp,'(I6.6)') n_lat
       write(lontmp,'(I6.6)') n_lon
       latlon_gridname = "DE"//lontmp//"x"//"PE"//lattmp

       call LDAS_create_grid_g(latlon_gridname,n_lon,n_lat, latlon_grid_tmp,i_off,j_off)

       pert_grid = latlon_grid_tmp  

    endif

  end function get_pert_grid

  ! *********************************************************************

  ! CHANGED: No longer using apply_force/progn_pert()
  ! The GridComp directly calls apply_pert()

  ! *********************************************************************

  subroutine apply_pert( pert_param, Pert, F, dt, rc )

    implicit none

    ! apply the perturbation Pert to a 1d field F
    !
    ! If the optional argument "dt" is present, Pert is interpreted
    ! as a perturbation flux forcing on the field F, and the perturbations
    ! are computed as follows:
    !
    !  F = F + Pert*dt    for additive perturbations
    !  F = F * Pert**dt   for multiplicative perturbations
    !
    ! Note that the units of "dt" must be consistent with those of "Pert".
    !
    ! If "dt" is not present, perturbations are computed as follows:
    !
    !  F = F + Pert       for additive perturbations
    !  F = F * Pert       for multiplicative perturbations
    !
    ! Note that pert_param is a scalar of type pert_param_type, so this
    ! subroutine typically appears within a nested loop from 1 through N_pert.
    !
    ! reichle,  1 Jun 2005
    ! reichle, 21 Nov 2014 - added optional interpretation of Pert as flux
    !
    ! -------------------------------------------------------------------

    type(pert_param_type),   intent(in)              :: pert_param

    real, dimension(:),      intent(in)              :: Pert

    real, dimension(:),      intent(inout)           :: F

    real,                    intent(in),    optional :: dt
    
    integer,                 intent(out),   optional :: rc

    ! local variables

    character(len=*), parameter :: Iam = 'apply_pert'

    ! -----------------------------------------------------------

    _ASSERT(size(Pert)==size(F), "sizes of Pert and perturbed field do not match")

    select case (pert_param%typ)

    case (0)

       if (present(dt)) then

          F = F + Pert*dt

       else

          F = F + Pert

       end if

    case (1)

       if (present(dt)) then

          F = F * Pert**dt

       else

          F = F * Pert

       end if

    case default
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown pert_param%typ')

    end select

  end subroutine apply_pert

  ! ******************************************************************

  subroutine get_progn_pert_param( pert_grid_l, N_progn_pert, &
       progn_pert_param )

    ! get parameters of perturbations to prognostic variables
    !
    ! reichle, 27 Nov 2001
    ! reichle, 31 May 2005 - redesign

    implicit none

    ! ---------------------------------------------------------------

    type(grid_def_type), intent(in) :: pert_grid_l

    integer, intent(out) :: N_progn_pert

    type(pert_param_type), dimension(:), pointer :: progn_pert_param ! output

    ! ---------------------------------------------------------------

    ! local variables

    real, dimension(N_progn_pert_max,pert_grid_l%N_lon,pert_grid_l%N_lat) :: &
         std_progn_pert

    character(40), dimension(N_progn_pert_max) :: descr_progn_pert

    logical, dimension(N_progn_pert_max ):: zeromean_progn_pert
    logical, dimension(N_progn_pert_max ):: coarsen_progn_pert

    real, dimension(N_progn_pert_max) :: std_normal_max_progn_pert
    real, dimension(N_progn_pert_max) :: xcorr_progn_pert, ycorr_progn_pert
    real, dimension(N_progn_pert_max) :: tcorr_progn_pert, typ_progn_pert

    real, dimension(N_progn_pert_max,N_progn_pert_max) :: ccorr_progn_pert

    integer, dimension(N_progn_pert_max)  :: progn_pert_select

    ! ---------------------------------------------------------------

    call get_progn_pert_inputs( pert_grid_l,                               &
         descr_progn_pert, zeromean_progn_pert, coarsen_progn_pert,        &
         std_normal_max_progn_pert, std_progn_pert,                        &
         xcorr_progn_pert, ycorr_progn_pert,                               &
         tcorr_progn_pert, typ_progn_pert,   ccorr_progn_pert        )

    call get_pert_select( N_progn_pert_max, pert_grid_l, std_progn_pert,   &
         N_progn_pert, progn_pert_select )


    if (N_progn_pert>0) then

       call allocate_pert_param(N_progn_pert,    &
            pert_grid_l%N_lon,pert_grid_l%N_lat, &
            progn_pert_param)

       call assemble_pert_param( N_progn_pert_max, N_progn_pert, pert_grid_l, &
            descr_progn_pert, zeromean_progn_pert, coarsen_progn_pert,        &
            std_normal_max_progn_pert,                                        &
            std_progn_pert, xcorr_progn_pert, ycorr_progn_pert,               &
            tcorr_progn_pert, typ_progn_pert, ccorr_progn_pert,               &
            progn_pert_select, progn_pert_param )

    end if

  end subroutine get_progn_pert_param

  ! **********************************************************************

  subroutine get_force_pert_param( pert_grid_l, N_force_pert, force_pert_param )

    ! get parameters of forcing perturbations
    !
    ! reichle, 27 Nov 2001
    ! reichle, 19 Jul 2005

    implicit none

    ! ---------------------------------------------------------------

    type(grid_def_type), intent(in) :: pert_grid_l

    integer, intent(out) :: N_force_pert

    type(pert_param_type), dimension(:), pointer :: force_pert_param ! output

    ! ---------------------------------------------------------------

    ! local variables

    real, dimension(N_force_pert_max,pert_grid_l%N_lon,pert_grid_l%N_lat) :: &
         std_force_pert

    character(40), dimension(N_force_pert_max) :: descr_force_pert

    logical, dimension(N_force_pert_max ):: zeromean_force_pert
    logical, dimension(N_force_pert_max ):: coarsen_force_pert

    real, dimension(N_force_pert_max) :: std_normal_max_force_pert
    real, dimension(N_force_pert_max) :: xcorr_force_pert, ycorr_force_pert
    real, dimension(N_force_pert_max) :: tcorr_force_pert, typ_force_pert

    real, dimension(N_force_pert_max,N_force_pert_max) :: ccorr_force_pert

    integer, dimension(N_force_pert_max)  :: force_pert_select

    ! ---------------------------------------------------------------

    call get_force_pert_inputs( pert_grid_l,                               &
         descr_force_pert, zeromean_force_pert, coarsen_force_pert,        &
         std_normal_max_force_pert, std_force_pert,                        &
         xcorr_force_pert, ycorr_force_pert,                               &
         tcorr_force_pert, typ_force_pert,   ccorr_force_pert        )

    call get_pert_select( N_force_pert_max, pert_grid_l, std_force_pert,   &
         N_force_pert, force_pert_select )

    if (N_force_pert>0) then

       call allocate_pert_param(N_force_pert,    &
            pert_grid_l%N_lon,pert_grid_l%N_lat, &
            force_pert_param)

       call assemble_pert_param( N_force_pert_max, N_force_pert, pert_grid_l, &
            descr_force_pert, zeromean_force_pert, coarsen_force_pert,        &
            std_normal_max_force_pert,                                        &
            std_force_pert, xcorr_force_pert, ycorr_force_pert,               &
            tcorr_force_pert, typ_force_pert, ccorr_force_pert,               &
            force_pert_select, force_pert_param )

    end if

  end subroutine get_force_pert_param

  ! **********************************************************************

  subroutine get_force_pert_inputs( pert_grid_l,                         &
       descr_force_pert, zeromean_force_pert, coarsen_force_pert,        &
       std_normal_max_force_pert, std_force_pert,                        &
       xcorr_force_pert, ycorr_force_pert,                               &
       tcorr_force_pert, typ_force_pert,   ccorr_force_pert       )

    ! get inputs for forcing perturbations for ALL forcing
    ! variables (including zero standard deviations) on a grid (typically,
    ! the subgrid of the tile definition grid that covers the domain, or
    ! a coarser version of that)
    !
    ! in subroutine() get_pert_select all forcing types with a nonzero
    ! standard deviation in at least one grid cell will be included
    ! in the forcing perturbations
    !
    ! this structure should allow for easy implementation of heterogeneous
    ! forcing perturbations std's in the future (in this case the std's would
    ! probably be read from a file or modified from a nominal value
    ! according to the properties of a given catchment)
    !
    ! parameters are obtained from namelist file as structures with
    ! fields corresponding to the kinds of forcing perturbations, then
    ! moved into straight (multidimensional) arrays for further processing
    ! outside of this subroutine
    !
    ! reichle, 27 Nov 2001
    ! reichle, 24 Mar 2004 - revised for enkf inputs namelist
    ! reichle, 27 May 2005 - redesign
    ! reichle+pchakrab, 17 May 2013 - parallelized perturbations
    !                                 EXCEPT I/O of distributed %std, %ccorr
    !
    ! ---------------------------------------------------------------

    implicit none

    ! ---------------------------------------------------------------

    type(grid_def_type), intent(in) :: pert_grid_l

    character(40), dimension(N_force_pert_max), intent(out) :: &
         descr_force_pert

    logical, dimension(N_force_pert_max), intent(out) :: zeromean_force_pert
    logical, dimension(N_force_pert_max), intent(out) :: coarsen_force_pert

    real,                                                                 &
         dimension(N_force_pert_max,pert_grid_l%N_lon,pert_grid_l%N_lat), &
         intent(out) :: std_force_pert

    real, dimension(N_force_pert_max), intent(out) :: &
         std_normal_max_force_pert, &
         xcorr_force_pert, &
         ycorr_force_pert, &
         tcorr_force_pert, &
         typ_force_pert

    real, dimension(N_force_pert_max,N_force_pert_max), intent(out) :: &
         ccorr_force_pert

    ! ---------------------------------------------------------------
    !
    ! local variables

    integer :: ii, jj

    type(force_pert_real_type) :: tmp_force_pert_real

    type(force_pert_logi_type) :: tmp_force_pert_logical

    type(force_pert_char_type) :: tmp_force_pert_character

    real,    dimension(N_force_pert_max) :: tmp_force_pert_vec

    logical, dimension(N_force_pert_max) :: stdfromfile_force_pert

    character(300)                       :: stdfilename_force_pert

    type(force_pert_ccor_type) :: tmp_force_pert_ccorr

    integer :: ccorr_size

    ! pchakrab: variables for reading NetCDF-4 file
    integer :: ivar, iproc ! counters
    integer :: nc4_stat ! return code of nc4 function calls
    integer :: nc4_id, nc4_grpid, nc4_varid ! various IDs
    character( 10) :: nc4_varname
    character(300) :: nc4_file
    real :: nc4_fillval

    integer :: xstart, xcount, ystart, ycount ! for computing local indices

    ! Errlong variables
    integer :: rc, status
    character(len=*), parameter :: Iam = 'get_force_pert_inputs'

    ! MPI variables
    integer :: mpierr
    logical :: root_proc

    ! -----------------------------------------------------------------

    root_proc = (myid==0)

    ! ---------
    !
    ! DESCR

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_descr_force_pert=tmp_force_pert_character)

       ! move data from structure into regular vector

       call struct2vec_force_pert(tmp_force_pert_character,descr_force_pert)
    endif

    ! for now, broadcast each element of the array descr_force_pert individually
    ! TODO: send the array in one MPI_Bcast call
    do ii=1,N_force_pert_max
       call MPI_Bcast(descr_force_pert(ii), 40, MPI_CHARACTER, 0, mpicomm, mpierr)
    end do

    ! ----------
    !
    ! ZEROMEAN

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_zeromean_force_pert=tmp_force_pert_logical)

       ! move data from structure into regular vector

       call struct2vec_force_pert(tmp_force_pert_logical,zeromean_force_pert)
    end if

    ! broadcast zeromean_force_pert
    call MPI_Bcast(zeromean_force_pert, N_force_pert_max, MPI_LOGICAL, 0, mpicomm, mpierr)

    ! ----------
    !
    ! COARSEN

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_coarsen_force_pert=tmp_force_pert_logical)

       ! move data from structure into regular vector

       call struct2vec_force_pert(tmp_force_pert_logical,coarsen_force_pert)
    end if

    ! broadcast coarsen_force_pert
    call MPI_Bcast(coarsen_force_pert, N_force_pert_max, MPI_LOGICAL, 0, mpicomm, mpierr)

    ! -----------------------------------------------------------------
    !
    ! STD

    ! obtain (default) homogeneous std of forcing perturbations

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_std_force_pert=tmp_force_pert_real)

       ! move data from structure into regular vector
       call struct2vec_force_pert(tmp_force_pert_real, tmp_force_pert_vec)
    endif

    ! broadcast tmp_force_pert_vec
    call MPI_Bcast(tmp_force_pert_vec, N_force_pert_max, MPI_REAL, 0, mpicomm, mpierr)


    ! initialize std_force_pert to homogeneous value

    do ivar=1,N_force_pert_max

       std_force_pert(ivar,:,:) = tmp_force_pert_vec(ivar)

    end do


    ! find out whether std_force_pert should be read from file

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_stdfromfile_force_pert=tmp_force_pert_logical)

       ! move data from structure into regular vector

       call struct2vec_force_pert(tmp_force_pert_logical,stdfromfile_force_pert)
    end if

    ! broadcast stdfromfile_force_pert
    call MPI_Bcast(stdfromfile_force_pert, N_force_pert_max, MPI_LOGICAL, 0, mpicomm, mpierr)

    ! read std_force_pert from file as needed

    if (any(stdfromfile_force_pert)) then

       ! find out name (incl full path) of file with std value

       if (root_proc)                                               &
            call read_ens_prop_inputs(                                &
            kw_stdfilename_force_pert = stdfilename_force_pert        &
            )

       call MPI_BCAST(stdfilename_force_pert,300,MPI_CHARACTER,0,mpicomm,mpierr)

       nc4_file = stdfilename_force_pert
       
       ! NOTE: the input file is in netcdf format, with a group 'std_force_pert',
       ! and the grid in the netcdf file must be the *global* pert grid 
       ! (see subroutine get_pert_grid()) 
       
       ! --compute-local-shape-first-
       ! ASSUMPTION: data in file are on the *global* pert grid
       xstart = pert_grid_l%i_offg + 1
       xcount = pert_grid_l%N_lon
       ystart = pert_grid_l%j_offg + 1
       ycount = pert_grid_l%N_lat

       ! pchakrab - 05/13/2014
       ! added reading herterogeneous std_force/progn_pert from file
       ! NOTE: With the current version of Baselibs (3.3.3), we cannot
       ! read 'in parallel'. So, we let each proc read the file in turn
       ! and read *its* data. Once, we transition to Baselibs v4, we will
       ! open the file in parallel and read the relevant part of the data

       do iproc=0,numprocs-1
          if (myid==iproc) then
             ! open file
             nc4_stat = nf90_open(path=nc4_file, mode=NF90_NOWRITE, ncid=nc4_id)
             if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
             nc4_stat = nf90_inq_ncid(nc4_id, 'std_force_pert', nc4_grpid)
             if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
             do ivar=1,N_force_pert_max
                if (stdfromfile_force_pert(ivar)) then
                   nc4_varname = trim(descr_force_pert(ivar))
                   ! id the the ivar-th variable
                   nc4_stat = nf90_inq_varid(nc4_grpid, nc4_varname, nc4_varid)
                   if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
                   ! read variable for iproc into std_force_pert(ivar,:,:)
                   ! --read-iproc's-slice-of-data-
                   nc4_stat = nf90_get_var( &
                        nc4_grpid, &
                        nc4_varid, &
                        values=std_force_pert(ivar,:,:), &
                        start=(/xstart,ystart/), &
                        count=(/xcount,ycount/) &
                        )
                   if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
                   ! get _FillValue for nc4_varname
                   nc4_stat = nf90_get_att(nc4_grpid, nc4_varid, '_FillValue', nc4_fillval)
                   if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
                   ! replace _FillValue by zero
                   where (abs(std_force_pert(ivar,:,:)-nc4_fillval)<nodata_tolfrac_generic) &
                        std_force_pert(ivar,:,:) = 0.
                end if
             end do
          end if
          ! close file
          nc4_stat = nf90_close(nc4_id)
          if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
       end do
       call MPI_Barrier(mpicomm, mpierr)
    end if

    ! --------------------------------------------------------------
    !
    ! STD_NORMAL_MAX, XCORR, YCORR, TCORR, TYP
    !
    ! NOTE: xcorr_force_pert, ycorr_force_pert, and tcorr_force_pert must be
    !       homogeneous (ie same for all catchments, unlike std_force_pert)
    !       typ_force_pert must also be homogeneous

    ! PC: instead of reading one param (and broadcasting it) at a time, it
    ! will be better to read them all and broadcast at one go

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_std_normal_max_force_pert=tmp_force_pert_real)

       ! move data from structure into regular vector

       call struct2vec_force_pert( tmp_force_pert_real, &
            std_normal_max_force_pert )
    end if

    ! broadcast std_normal_max_force_pert
    call MPI_Bcast(std_normal_max_force_pert, N_force_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_xcorr_force_pert=tmp_force_pert_real)

       ! move data from structure into regular vector

       call struct2vec_force_pert( tmp_force_pert_real, xcorr_force_pert )
    end if

    ! broadcast xcorr_force_pert
    call MPI_Bcast(xcorr_force_pert, N_force_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_ycorr_force_pert=tmp_force_pert_real)

       ! move data from structure into regular vector

       call struct2vec_force_pert( tmp_force_pert_real, ycorr_force_pert )
    end if

    ! broadcast ycorr_force_pert
    call MPI_Bcast(ycorr_force_pert, N_force_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_tcorr_force_pert=tmp_force_pert_real)

       ! move data from structure into regular vector

       call struct2vec_force_pert( tmp_force_pert_real, tcorr_force_pert )
    end if

    ! broadcast tcorr_force_pert
    call MPI_Bcast(tcorr_force_pert, N_force_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_typ_force_pert=tmp_force_pert_real)

       ! move data from structure into regular vector

       call struct2vec_force_pert( tmp_force_pert_real, typ_force_pert )
    end if

    ! broadcast type_force_pert
    call MPI_Bcast(typ_force_pert, N_force_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------
    !
    ! only the minimum number of elements of the correlation matrix is
    ! specified in the namelist file, the rest is initialized to nodata_generic
    ! (see subroutine read_ens_prop_inputs)
    ! now fill in the rest of the information (diagonal=1 and symmetry)

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_ccorr_force_pert=tmp_force_pert_ccorr)

       ! move data from structure into regular matrix

       call struct2mat_force_pert_ccor(tmp_force_pert_ccorr, ccorr_force_pert )
    endif

    ! reichle+pchakrab, 17 May 2013 - parallelized perturbations
    !                                 EXCEPT I/O of distributed %std, %ccorr

    ! broadcast ccorr_force_pert
    ccorr_size = N_force_pert_max*N_force_pert_max
    call MPI_Bcast(ccorr_force_pert, ccorr_size, MPI_REAL, 0, mpicomm, mpierr)

    ! set diagonal to 1. and make matrix symmetric

    do ii=1,N_force_pert_max
       do jj=1,N_force_pert_max

          if (ii==jj) then

             ccorr_force_pert(ii,ii) = 1.

          elseif (abs(ccorr_force_pert(ii,jj)-nodata_generic)<nodata_tol_generic) then

             ccorr_force_pert(ii,jj) = ccorr_force_pert(jj,ii)

          end if

       end do
    end do

  end subroutine get_force_pert_inputs

  ! *********************************************************************

  subroutine get_progn_pert_inputs( pert_grid_l,                         &
       descr_progn_pert, zeromean_progn_pert, coarsen_progn_pert,        &
       std_normal_max_progn_pert, std_progn_pert,                        &
       xcorr_progn_pert, ycorr_progn_pert,                               &
       tcorr_progn_pert, typ_progn_pert,   ccorr_progn_pert       )

    ! get inputs for perturbations to prognostic variables for ALL prognostic
    ! variables (including zero standard deviations) on a grid (typically,
    ! the subgrid of the tile definition grid that covers the domain, or
    ! a coarser version of that)
    !
    ! in subroutine() get_err_select all states with a nonzero
    ! standard deviation in at least one grid cell will be included
    ! in the perturbations to prognostic variables
    !
    ! this structure should allow for easy implementation of heterogeneous
    ! perturbations std's in the future (in this case the std's would
    ! probably be read from a file or modified from a nominal value
    ! according to the properties of a given catchment)
    !
    ! parameters are obtained from namelist file as structures with
    ! fields corresponding to the kinds of perturbations, then
    ! moved into straight (multidimensional) arrays for further processing
    ! outside of this subroutine
    !
    ! reichle, 27 Nov 2001
    ! reichle, 24 Mar 2004 - revised for enkf inputs namelist
    ! reichle, 27 May 2005 - redesign
    ! reichle+pchakrab, 17 May 2013 - parallelized perturbations
    !                                 EXCEPT I/O of distributed %std, %ccorr
    !
    ! ---------------------------------------------------------------

    implicit none

    type(grid_def_type), intent(in) :: pert_grid_l

    character(40), dimension(N_progn_pert_max), intent(out) :: &
         descr_progn_pert

    logical, dimension(N_progn_pert_max), intent(out) :: zeromean_progn_pert
    logical, dimension(N_progn_pert_max), intent(out) :: coarsen_progn_pert

    real,                                                                 &
         dimension(N_progn_pert_max,pert_grid_l%N_lon,pert_grid_l%N_lat), &
         intent(out) :: std_progn_pert

    real, dimension(N_progn_pert_max), intent(out) :: &
         std_normal_max_progn_pert, &
         xcorr_progn_pert, &
         ycorr_progn_pert, &
         tcorr_progn_pert, &
         typ_progn_pert

    real, dimension(N_progn_pert_max,N_progn_pert_max), intent(out) :: &
         ccorr_progn_pert

    ! ---------------------------------------------------------------
    !
    ! local variables

    integer :: ii, jj

    type(progn_pert_real_type) :: tmp_progn_pert_real

    type(progn_pert_logi_type) :: tmp_progn_pert_logical

    type(progn_pert_char_type) :: tmp_progn_pert_character

    real,    dimension(N_progn_pert_max) :: tmp_progn_pert_vec

    logical, dimension(N_progn_pert_max) :: stdfromfile_progn_pert

    character(300)                       :: stdfilename_progn_pert

    type(progn_pert_ccor_type) :: tmp_progn_pert_ccorr

    integer :: ccorr_size

    ! pchakrab: variables for reading NetCDF-4 file
    integer :: ivar, iproc ! counters
    integer :: nc4_stat ! return code of nc4 function calls
    integer :: nc4_id, nc4_grpid, nc4_varid ! various IDs
    character( 10) :: nc4_varname
    character(300) :: nc4_file
    real :: nc4_fillval

    integer :: xstart, xcount, ystart, ycount ! for computing local indices

    ! Errlong variables
    integer :: rc, status
    character(len=*), parameter :: Iam = 'get_progn_pert_inputs'

    ! MPI variables
    integer :: mpierr
    logical :: root_proc

    ! -----------------------------------------------------------------

    root_proc = (myid==0)

    ! -------
    !
    ! DESCR

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_descr_progn_pert=tmp_progn_pert_character)

       ! move data from structure into regular vector

       call struct2vec_progn_pert(tmp_progn_pert_character,descr_progn_pert)
    end if

    ! for now, broadcast each element of the array descr_progn_pert individually
    ! TODO: send the array in one MPI_Bcast call
    do ii=1,N_progn_pert_max
       call MPI_Bcast(descr_progn_pert(ii), 40, MPI_CHARACTER, 0, mpicomm, mpierr)
    end do

    ! ----------
    !
    ! ZEROMEAN

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_zeromean_progn_pert=tmp_progn_pert_logical)

       ! move data from structure into regular vector

       call struct2vec_progn_pert(tmp_progn_pert_logical,zeromean_progn_pert)
    end if

    ! broadcast zeromean_progn_pert
    call MPI_Bcast(zeromean_progn_pert, N_progn_pert_max, MPI_LOGICAL, 0, mpicomm, mpierr)

    ! ----------
    !
    ! COARSEN

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_coarsen_progn_pert=tmp_progn_pert_logical)

       ! move data from structure into regular vector

       call struct2vec_progn_pert(tmp_progn_pert_logical,coarsen_progn_pert)
    end if

    ! broadcast coarsen_progn_pert
    call MPI_Bcast(coarsen_progn_pert, N_progn_pert_max, MPI_LOGICAL, 0, mpicomm, mpierr)

    ! -----------------------------------------------------------------
    !
    ! STD
    !
    ! obtain (default) homogeneous std of forcing perturbations

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_std_progn_pert=tmp_progn_pert_real)

       ! move data from structure into regular vector
       call struct2vec_progn_pert(tmp_progn_pert_real, tmp_progn_pert_vec)
    end if

    ! broadcast tmp_progn_pert_vec
    call MPI_Bcast(tmp_progn_pert_vec, N_progn_pert_max, MPI_REAL, 0, mpicomm, mpierr)


    ! initialize std_progn_pert to homogeneous value

    do ivar=1,N_progn_pert_max

       std_progn_pert(ivar,:,:) = tmp_progn_pert_vec(ivar)

    end do


    ! find out whether std_progn_pert should be read from file

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_stdfromfile_progn_pert=tmp_progn_pert_logical)

       ! move data from structure into regular vector

       call struct2vec_progn_pert(tmp_progn_pert_logical,stdfromfile_progn_pert)
    end if

    ! broadcast stdfromfile_progn_pert
    call MPI_Bcast(stdfromfile_progn_pert, N_progn_pert_max, MPI_LOGICAL, 0, mpicomm, mpierr)


    ! read std_progn_pert from file as needed

    if (any(stdfromfile_progn_pert)) then

       ! find out name (incl full path) of file with std value

       if (root_proc)                                               &
            call read_ens_prop_inputs(                                &
            kw_stdfilename_progn_pert = stdfilename_progn_pert        &
            )

       call MPI_BCAST(stdfilename_progn_pert,300,MPI_CHARACTER,0,mpicomm,mpierr)

       nc4_file = stdfilename_progn_pert
       
       ! NOTE: the input file is in netcdf format, with a group 'std_force_pert',
       ! and the grid in the netcdf file must be the *global* pert grid 
       ! (see subroutine get_pert_grid()) 

       ! --compute-local-shape-first-
       ! ASSUMPTION: data in file are on the *global* pert grid 
       xstart = pert_grid_l%i_offg + 1
       xcount = pert_grid_l%N_lon
       ystart = pert_grid_l%j_offg + 1
       ycount = pert_grid_l%N_lat

       ! pchakrab - 05/13/2014
       ! added reading herterogeneous std_force/progn_pert from file
       ! NOTE: With the current version of Baselibs (3.3.3), we cannot
       ! read 'in parallel'. So, we let each proc read the file in turn
       ! and read *its* data. Once, we transition to Baselibs v4, we will
       ! open the file in parallel and read the relevant part of the data

       do iproc=0,numprocs-1
          if (myid==iproc) then
             ! open file
             nc4_stat = nf90_open(path=nc4_file, mode=NF90_NOWRITE, ncid=nc4_id)
             if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
             nc4_stat = nf90_inq_ncid(nc4_id, 'std_progn_pert', nc4_grpid)
             if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
             do ivar=1,N_progn_pert_max
                if (stdfromfile_progn_pert(ivar)) then
                   nc4_varname = trim(descr_progn_pert(ivar))
                   ! id the the ivar-th variable
                   nc4_stat = nf90_inq_varid(nc4_grpid, nc4_varname, nc4_varid)
                   if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
                   ! read variable for iproc into std_progn_pert(ivar,:,:)
                   ! --read-iproc's-slice-of-data-
                   nc4_stat = nf90_get_var( &
                        nc4_grpid, &
                        nc4_varid, &
                        values=std_progn_pert(ivar,:,:), &
                        start=(/xstart,ystart/), &
                        count=(/xcount,ycount/) &
                        )
                   if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
                   ! get _FillValue for nc4_varname
                   nc4_stat = nf90_get_att(nc4_grpid, nc4_varid, '_FillValue', nc4_fillval)
                   if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
                   ! replace _FillValue by zero
                   where (abs(std_progn_pert(ivar,:,:)-nc4_fillval)<nodata_tolfrac_generic) &
                        std_progn_pert(ivar,:,:) = 0.
                end if
             end do
             ! close file
             nc4_stat = nf90_close(nc4_id)
             if (nc4_stat /= nf90_noerr) call handle_nc4_stat(nc4_stat)
          end if
       end do
       call MPI_Barrier(mpicomm, mpierr)
    end if

    ! --------------------------------------------------------------
    !
    ! STD_NORMAL_MAX, XCORR, YCORR, TCORR, TYP
    !
    ! NOTE: xcorr_progn_pert, ycorr_progn_pert, and tcorr_progn_pert must be
    !       homogeneous (ie same for all catchments, unlike std_progn_pert)
    !       typ_progn_pert must also be homogeneous

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_std_normal_max_progn_pert=tmp_progn_pert_real)

       ! move data from structure into regular vector

       call struct2vec_progn_pert( tmp_progn_pert_real, &
            std_normal_max_progn_pert )
    end if

    ! broadcast std_normal_max_progn_pert
    call MPI_Bcast(std_normal_max_progn_pert, N_progn_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_xcorr_progn_pert=tmp_progn_pert_real)

       ! move data from structure into regular vector

       call struct2vec_progn_pert( tmp_progn_pert_real, xcorr_progn_pert )
    end if

    ! broadcast xcorr_progn_pert
    call MPI_Bcast(xcorr_progn_pert, N_progn_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_ycorr_progn_pert=tmp_progn_pert_real)

       ! move data from structure into regular vector

       call struct2vec_progn_pert( tmp_progn_pert_real, ycorr_progn_pert )
    end if

    ! broadcast ycorr_progn_pert
    call MPI_Bcast(ycorr_progn_pert, N_progn_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_tcorr_progn_pert=tmp_progn_pert_real)

       ! move data from structure into regular vector

       call struct2vec_progn_pert( tmp_progn_pert_real, tcorr_progn_pert )
    end if

    ! broadcast tcorr_progn_pert
    call MPI_Bcast(tcorr_progn_pert, N_progn_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_typ_progn_pert=tmp_progn_pert_real)

       ! move data from structure into regular vector

       call struct2vec_progn_pert( tmp_progn_pert_real, typ_progn_pert )
    endif

    ! broadcast type_progn_pert
    call MPI_Bcast(typ_progn_pert, N_progn_pert_max, MPI_REAL, 0, mpicomm, mpierr)

    ! ----------
    !
    ! only the minimum number of elements of the correlation matrix is
    ! specified in the namelist file, the rest is initialize to nodata_generic
    ! (see subroutine read_ens_prop_inputs)
    ! now fill in the rest of the information (diagonal=1 and symmetry)

    if (root_proc) then
       call read_ens_prop_inputs(kw_echo=.false., &
            kw_ccorr_progn_pert=tmp_progn_pert_ccorr)

       ! move data from structure into regular matrix

       call struct2mat_progn_pert_ccor(tmp_progn_pert_ccorr, ccorr_progn_pert )
    end if

    ! reichle+pchakrab, 17 May 2013 - parallelized perturbations
    !                                 EXCEPT I/O of distributed %std, %ccorr

    ! broadcast ccorr_progn_pert
    ccorr_size = N_progn_pert_max*N_progn_pert_max
    call MPI_Bcast(ccorr_progn_pert, ccorr_size, MPI_REAL, 0, mpicomm, mpierr)

    ! set diagonal to 1. and make matrix symmetric

    do ii=1,N_progn_pert_max
       do jj=1,N_progn_pert_max

          if (ii==jj) then

             ccorr_progn_pert(ii,ii) = 1.

          elseif (abs(ccorr_progn_pert(ii,jj)-nodata_generic)<nodata_tol_generic) then

             ccorr_progn_pert(ii,jj) = ccorr_progn_pert(jj,ii)

          end if

       end do
    end do

  end subroutine get_progn_pert_inputs

  ! *********************************************************************

  subroutine assemble_pert_param(                                &
       N_pert_max, N_pert, pert_grid_l,                          &
       descr_pert, zeromean_pert, coarsen_pert,                  &
       std_normal_max_pert,                                      &
       std_pert, xcorr_pert, ycorr_pert, tcorr_pert, typ_pert,   &
       ccorr_pert,                                               &
       pert_select, pert_param )

    ! assemble perturbation parameters for the prognostics or forcings
    !  that have been selected to be subject to perturbations
    !
    ! reichle, 30 Nov 2001
    ! reichle, 27 May 2005 - redesign
    !
    ! -----------------------------------------------------------------

    implicit none

    ! -----------------------------------------------------------------

    integer, intent(in) :: N_pert_max, N_pert

    type(grid_def_type) :: pert_grid_l

    character(40), dimension(N_pert_max), intent(in) :: descr_pert

    logical, dimension(N_pert_max), intent(in) :: zeromean_pert
    logical, dimension(N_pert_max), intent(in) :: coarsen_pert

    real,                                                             &
         dimension(N_pert_max,pert_grid_l%N_lon,pert_grid_l%N_lat),   &
         intent(in) :: std_pert

    real, dimension(N_pert_max), intent(in) :: std_normal_max_pert
    real, dimension(N_pert_max), intent(in) :: xcorr_pert
    real, dimension(N_pert_max), intent(in) :: ycorr_pert
    real, dimension(N_pert_max), intent(in) :: tcorr_pert
    real, dimension(N_pert_max), intent(in) :: typ_pert

    real, dimension(N_pert_max,N_pert_max), intent(in) :: ccorr_pert

    integer, dimension(N_pert_max), intent(in) :: pert_select

    type(pert_param_type), dimension(:), pointer :: pert_param ! output

    ! --------------------------

    ! local variables

    integer :: i,j,k,l,m,n

    real :: tmpreal

    real, dimension(N_pert,N_pert) :: tmpmat1, tmpmat2

    character(len=400) :: err_msg
    character(len=*), parameter :: Iam = 'assemble_pert_param'

    ! -----------------------------------------------------------------

    ! extract info into pert_param of length N_pert
    ! (only those states listed in pert_select are affected by perturbations)

    m = 0

    do k=1,N_pert_max

       if (pert_select(k)>0) then

          m = m+1

          pert_param(m)%descr = trim(descr_pert(k))

          pert_param(m)%zeromean = zeromean_pert(k)
          pert_param(m)%coarsen  = coarsen_pert( k)

          pert_param(m)%std_normal_max = std_normal_max_pert(k)

          pert_param(m)%std   = std_pert(k,:,:)

          pert_param(m)%xcorr = xcorr_pert(k)
          pert_param(m)%ycorr = ycorr_pert(k)
          pert_param(m)%tcorr = tcorr_pert(k)
          pert_param(m)%typ   = nint(typ_pert(k))


          n = 0

          do l=1,N_pert_max

             if (pert_select(l)>0) then

                n = n+1

                pert_param(m)%ccorr(n,:,:) = ccorr_pert(k,l)

             end if

          end do


       end if

    end do

    ! -------------------------------------------------------------
    !
    ! set mean and (if needed) modify standard deviation according to 'typ'
    ! (additive or multiplicative perturbations)

    do m=1,N_pert

       select case (pert_param(m)%typ)

       case (0)         ! additive (mean=0, std as above)

          pert_param(m)%mean  = 0.

       case (1)         ! multiplicative and lognormal (mean=1)

          do i=1,pert_grid_l%N_lon
             do j=1,pert_grid_l%N_lat

                tmpreal = pert_param(m)%std(i,j)

                tmpreal = log( 1. + tmpreal**2)

                pert_param(m)%mean(i,j) = - .5*tmpreal
                pert_param(m)%std(i,j)  = sqrt(tmpreal)

             end do
          end do

       case default

          err_msg = 'unknown typ of perturbation'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end select

    end do

    ! -------------------------------------------------------------
    !
    ! compute sqrt of correlation matrix for each grid point

    do i=1,pert_grid_l%N_lon
       do j=1,pert_grid_l%N_lat

          ! extract local correlation matrix for grid point (i,j)

          do m=1,N_pert
             do n=1,N_pert

                tmpmat1(m,n) = pert_param(m)%ccorr(n,i,j)

             end do
          end do

          ! compute sqrt of local correlation matrix

          call get_sqrt_corr_matrix( N_pert, tmpmat1, tmpmat2 )

          ! overwrite cross-correlations in forcepert_param with square
          ! root of cross-correlation matrix

          do m=1,N_pert
             do n=1,N_pert

                pert_param(m)%ccorr(n,i,j) = tmpmat2(m,n)

             end do
          end do

       end do
    end do

  end subroutine assemble_pert_param


  ! *********************************************************************

  subroutine get_pert_select( N_pert_max, pert_grid_l, std_pert, &
       N_pert, pert_select )

    ! determine which components of the prognostics and forcing perturbations
    ! are turned on (as identified through nonzero standard deviation)
    !
    ! input:
    !  N_pert_max  : number of prognostic or forcing perturbation variables
    !
    ! output:
    !  pert_select : vector of length N_pert_max, components are
    !                - zero if ALL grid points have zero perturbation std
    !                - one if at least one grid point has nonzero std
    !
    ! reichle, 30 Nov 2001
    ! reichle, 27 May 2005 - redesign
    !
    ! ---------------------------------------------------------------

    implicit none

    integer, intent(in) :: N_pert_max

    type(grid_def_type), intent(in) :: pert_grid_l

    real,                                                             &
         dimension(N_pert_max,pert_grid_l%N_lon,pert_grid_l%N_lat),   &
         intent(in) :: std_pert

    integer, intent(out) :: N_pert

    integer, intent(out), dimension(N_pert_max) :: pert_select

    ! ---------------------------------------------------------------
    !
    ! local variables

    integer :: i,j,k, ierr

    ! ---------------------------------------------------------------
    !
    ! record which components of the prognostics/forcings are affected
    !  with non-zero error standard deviation

    pert_select(1:N_pert_max) = 0

    do k=1,N_pert_max
       do i=1,pert_grid_l%N_lon
          do j=1,pert_grid_l%N_lat

             if (std_pert(k,i,j)>0.) pert_select(k) = 1

          end do
       end do
    end do

    call MPI_Allreduce(MPI_IN_PLACE, pert_select, N_pert_max , MPI_INTEGER, MPI_MAX, mpicomm, ierr )

    N_pert = sum( pert_select )

  end subroutine get_pert_select

  ! *********************************************************************

  subroutine interpolate_pert_to_timestep(                 &
       date_time, pert_time_old, pert_dtstep_real,         &
       Pert_old, Pert_new, Pert_ntp )

    ! Linearly interpolates perturbations to model time step
    !
    ! "_old" = at old time
    ! "_new" = at new time
    ! "_ntp" = at current ("interpolated") time
    !
    ! reichle, 26 May 2005
    !
    ! pchakrab, 24 Feb 2015 - using assumed shape arrays

    implicit none

    type(date_time_type), intent(in) :: date_time, pert_time_old

    real, intent(in) :: pert_dtstep_real

    real, dimension(:,:), intent(in)  :: Pert_old, Pert_new

    real, dimension(:,:), intent(out) :: Pert_ntp

    ! ----------------

    ! local variables

    real :: w

    ! ------------------------------------------------------------
    !
    ! weight for interpolation

    w = real(datetime2_minus_datetime1( pert_time_old, date_time ))

    w = w/pert_dtstep_real

    Pert_ntp = (1.-w)*Pert_old + w*Pert_new

  end subroutine interpolate_pert_to_timestep

  ! ***********************************************************************

  subroutine echo_pert_param( N_pert, pert_param, ind_i, ind_j )

    ! echo pert_param for grid point (ind_i,ind_j)

    implicit none

    integer, intent(in) :: N_pert, ind_i, ind_j

    type(pert_param_type), dimension(:), pointer :: pert_param

    ! locals

    integer :: m, n

    ! -------------------------------------------------------------

    if (root_logit) then
       write (logunit,*) 'echo_pert_param():'

       do m=1,N_pert

          write (logunit,*) 'pert_param(',m,')%descr=',          pert_param(m)%descr
          write (logunit,*) 'pert_param(',m,')%typ=',            pert_param(m)%typ
          write (logunit,*) 'pert_param(',m,')%zeromean=',       pert_param(m)%zeromean
          write (logunit,*) 'pert_param(',m,')%coarsen=',        pert_param(m)%coarsen

          write (logunit,*) 'pert_param(',m,')%std_normal_max=',                      &
            pert_param(m)%std_normal_max
          write (logunit,*) 'pert_param(',m,')%xcorr=',          pert_param(m)%xcorr
          write (logunit,*) 'pert_param(',m,')%ycorr=',          pert_param(m)%ycorr
          write (logunit,*) 'pert_param(',m,')%tcorr=',          pert_param(m)%tcorr

          write (logunit,*) 'pert_param(',m,')%mean(',ind_i,',',ind_j,')=',           &
            pert_param(m)%mean(ind_i,ind_j)
          write (logunit,*) 'pert_param(',m,')%std(',ind_i,',',ind_j,')=',            &
            pert_param(m)%std(ind_i,ind_j)

          do n=1,N_pert

             write (logunit,*) 'pert_param(',m,')%ccorr(',n,',',ind_i,',',ind_j,')=', &
               pert_param(m)%ccorr(n,ind_i,ind_j)

          end do
       end do
     endif ! root_logit
  end subroutine echo_pert_param

 !*************************************************************

  ! WY noted :: -l is changed to _g. Only read part was kept  
  subroutine io_pert_rstrt( action, work_path, exp_id, ens_id, &
       date_time, tile_coord_g, pert_grid_g, pert_grid_f,      &
       N_force_pert, N_progn_pert, Pert_rseed,                 &
       Force_pert_ntrmdt_g, Progn_pert_ntrmdt_g, rc )

    ! read or write perturbations re-start file.
    !
    ! reichle, 21 Jun 2005
    ! reichle, 16 Oct 2008 - added optional output "rc" (success/failure of read)
    ! Weiyuan, 30 Apr 2018 addapt to read into MAPL rst. 
    implicit none

    character,      intent(in) :: action     ! read or write

    character(200), intent(in) :: work_path
    character(40),  intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time

    integer, intent(in) :: ens_id, N_force_pert, N_progn_pert

    type(tile_coord_type), dimension(:), pointer :: tile_coord_g  ! input

    type(grid_def_type), intent(in) :: pert_grid_g, pert_grid_f

    integer, dimension(NRANDSEED), intent(inout) :: Pert_rseed

    real, dimension(pert_grid_g%N_lon,pert_grid_g%N_lat,N_force_pert), &
         intent(inout) :: Force_pert_ntrmdt_g

    real, dimension(pert_grid_g%N_lon,pert_grid_g%N_lat,N_force_pert), &
         intent(inout) :: Progn_pert_ntrmdt_g

    integer, intent(out), optional :: rc

    ! local variables

    integer        :: n, k, i, j, istat

    integer        :: nrandseed_tmp, N_force_pert_tmp, N_progn_pert_tmp

    type(grid_def_type) :: pert_grid_f_tmp

    character(300) :: filename

    character(40)  :: file_tag = 'pert_ldas_rst', dir_name='rs', file_ext='.bin'

    ! full array for reading/writing pert ntrmdt
    ! while reading, the data is read into Pert_ntrmdt_f and dispersed
    ! while writing, the data is assembled into Pert_ntrmdt_f and written
    real, dimension(:,:), pointer :: Pert_ntrmdt_f => null()

    character(len=*), parameter :: Iam = 'io_pert_rstrt'
    character(len=400) :: err_msg
    logical :: file_exists
    integer :: itmp,jtmp, xstart,xend, ystart,yend
    
    ! --------------------------------------------------------------------

    if (present(rc)) rc = 9999

    select case (action)

    case ('r','R')

       filename = get_io_filename( work_path, exp_id,               &
            file_tag, date_time=date_time,                          &
            dir_name=dir_name, ens_id=ens_id, file_ext=file_ext )

!!$       if (root_proc) then
       inquire(file=filename,exist=file_exists)
       if(.not. file_exists) then
          write (6,'(400A)') &
               'Warning :: Pert restart file NOT found: ' // trim(filename)
          if (present(rc))  rc = 1
          return
       endif

       open(10, file=filename, convert='big_endian',form='unformatted', status='old', &
            action='read', iostat=istat)
!!$       end if

!!$#ifdef LDAS_MPI
!!$       ! bcast the status of open (all procs need to return if open fails)
!!$       call MPI_Bcast(istat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
!!$#endif

!!$       if (root_proc) then

       write (6,'(400A)') &
               'Reading pert restart file ' // trim(filename)

          ! one additional header line (as of 21 May 2010)!!!

       call io_grid_def_type( 'r', 10, pert_grid_f_tmp )

       read (10) nrandseed_tmp, N_force_pert_tmp, N_progn_pert_tmp

          ! check whether entries in file match passed arguments
          ! (check does *not* include *_pert_param!)

       if ( (nrandseed_tmp          /= NRANDSEED)            .or.               &
               (N_force_pert_tmp       /= N_force_pert)         .or.               &
               (N_progn_pert_tmp       /= N_progn_pert) ) then
             err_msg = 'pert.rstrt file not compatible (1)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

       if ( index(pert_grid_f%gridtype,'LatLon')/=0  .or.         &
               index(pert_grid_f%gridtype,'LATLON')/=0  .or.         &
               index(pert_grid_f%gridtype,'latlon')/=0       ) then

             if ( (pert_grid_f_tmp%N_lon  /= pert_grid_f%N_lon)    .or.            &
                  (pert_grid_f_tmp%N_lat  /= pert_grid_f%N_lat)    .or.            &
                  (abs(pert_grid_f_tmp%ll_lon - pert_grid_f%ll_lon) > 1e-4) .or.   &
                  (abs(pert_grid_f_tmp%ll_lat - pert_grid_f%ll_lat) > 1e-4) .or.   &
                  (abs(pert_grid_f_tmp%dlon   - pert_grid_f%dlon)   > 1e-4) .or.   &
                  (abs(pert_grid_f_tmp%dlat   - pert_grid_f%dlat)   > 1e-4)      ) then
                err_msg = 'pert.rstrt file not compatible (2)'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if

       else

             if ( index(pert_grid_f_tmp%gridtype,pert_grid_f%gridtype)==0  .or.   &
                  (pert_grid_f_tmp%N_lon  /= pert_grid_f%N_lon)            .or.   &
                  (pert_grid_f_tmp%N_lat  /= pert_grid_f%N_lat)            .or.   &
                  (pert_grid_f_tmp%i_offg /= pert_grid_f%i_offg)           .or.   &
                  (pert_grid_f_tmp%j_offg /= pert_grid_f%j_offg)                ) then
                err_msg = 'pert.rstrt file not compatible (3)'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if

       end if


       read (10) (Pert_rseed(n), n=1,NRANDSEED)

          ! allocate memory for full pert ntrmdt array (to be dispersed)
       allocate (Pert_ntrmdt_f(pert_grid_f%N_lon, pert_grid_f%N_lat))

!!$       endif

       itmp = pert_grid_f%i_offg
       xstart = itmp + 1
       xend = itmp + pert_grid_f%N_lon
       jtmp = pert_grid_f%j_offg
       ystart = jtmp + 1
       yend = jtmp + pert_grid_f%N_lat

       do k=1,N_force_pert

!!$          if (root_proc) then
             read (10) ((Pert_ntrmdt_f(i,j), i=1,pert_grid_f%N_lon), &
                  j=1,pert_grid_f%N_lat)

             Force_pert_ntrmdt_g(xstart:xend, ystart:yend,k) = Pert_ntrmdt_f(:,:)
!!$          end if
!!$          call scatter_arr2d_grid(pert_grid_f, pert_grid_l, &
!!$               Pert_ntrmdt_f, Force_pert_ntrmdt_l(k,:,:))
!!$#ifdef LDAS_MPI
!!$          call MPI_Barrier(MPI_COMM_WORLD, mpierr)
!!$#endif

       end do

       do k=1,N_progn_pert

!!$          if (root_proc) then
             read (10) ((Pert_ntrmdt_f(i,j), i=1,pert_grid_f%N_lon), &
                  j=1,pert_grid_f%N_lat)
             Progn_pert_ntrmdt_g(xstart:xend, ystart:yend,k) = Pert_ntrmdt_f(:,:)
!!$          end if
!!$          call scatter_arr2d_grid(pert_grid_f, pert_grid_l, &
!!$               Pert_ntrmdt_f, Progn_pert_ntrmdt_l(k,:,:))
!!$#ifdef LDAS_MPI
!!$          call MPI_Barrier(MPI_COMM_WORLD, mpierr)
!!$#endif
       end do

       deallocate(Pert_ntrmdt_f)

       if (present(rc)) rc = 0


    case ('w','W')

       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'write part not needed any more')

    case default

       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown action')

    end select
    
!!    if (root_proc)  close (10,status='keep')
    close (10,status='keep')

   end subroutine io_pert_rstrt

  ! ******************************************************************

  ! handle return code of nf90_* calls
  ! stop on error
  subroutine handle_nc4_stat(status)
    ! input
    integer, intent (in) :: status
    ! local
    character(len=*), parameter :: Iam = 'handle_nc4_stat'
    character(len=400) :: err_msg

    if(status /= nf90_noerr) then
       err_msg = 'Stopped [' // trim(nf90_strerror(status)) // ']'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
  end subroutine handle_nc4_stat

  ! **********************************************************************

end module LDAS_PertRoutinesMod


! **** EOF ******************************************************
