! this file contains a collection of bias subroutines for the 
! the Ensemble Kalman filter of the Catchment model off-line driver
!
! reichle,        17 Oct 2005
! draper+reichle, 28 Aug 2013 - added obs bias routines
! draper+reichle, 19 Sep 2013 - revised obs bias routines

module clsm_bias_routines

  use catch_constants,                  ONLY:     &
       N_snow => CATCH_N_SNOW,                    &
       N_gt   => CATCH_N_GT
  
  use LDAS_ensdrv_globals,           ONLY:     &
       nodata_generic,                            &
       logit,                                     &
       logunit
  
  use  LDAS_DateTimeMod,                   ONLY:     &
       date_time_type
  
  use catch_types,                      ONLY:     &
       cat_progn_type,                            &
       assignment (=)

  use enkf_types,                       ONLY:     &
       obs_type,                                  &
       obs_param_type

  use catch_bias_types,                 ONLY:     &
       cat_bias_param_type,                       &
       obs_bias_type

  use LDAS_ensdrv_mpi,                  ONLY:     &
       master_proc,                               &
       MPI_obs_bias_type,                         &
       mpicomm,                            &
       MPIERR

  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename

  use LDAS_ensdrv_init_routines,        ONLY:     &
       clsm_ensdrv_get_command_line

  use clsm_ensupd_upd_routines,         ONLY:     &
       get_cat_progn_ens_avg

  use LDAS_exceptionsMod,                  ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR
    
  implicit none

  private
    
  public :: io_rstrt_cat_bias
  public :: read_cat_bias_inputs
  public :: init_cat_bias
  public :: cat_bias_corr
  public :: cat_bias_calcs_update

  public :: io_rstrt_obs_bias
  public :: init_obs_bias
  public :: obs_bias_upd_tcount
  public :: obs_bias_upd_bias_and_Obs
  public :: obs_bias_corr_obs
  public :: output_obs_bias
  public :: initialize_obs_bias
  
  integer, parameter :: max_tcount = 86400*365  ! ~1 year 
  
contains
    
  ! -------------------------------------------------------------------
  !
  ! CAT BIAS ROUTINES 
  !
  ! -------------------------------------------------------------------
  
  subroutine clsm_cat_bias_get_command_line(                    &
       cat_bias_inputs_path, cat_bias_inputs_file               &
       )
    
    ! get some inputs from command line 
    !
    ! if present, command line arguments overwrite inputs from
    ! catbias_inputs namelist files
    !
    ! command line should look something like
    !
    ! a.out -cat_bias_inputs_file fname.nml 
    !
    ! NOTE: This subroutine does NOT stop for unknown arguments!
    !       (If that is desired, all arguments used by 
    !        clsm_ensdrv_get_command_line() must be listed here
    !        explicitly and be ignored.)
    !
    ! reichle, 18 Oct 2005
    ! 
    ! ----------------------------------------------------------------
    
    implicit none
    
    character(*), intent(inout), optional :: cat_bias_inputs_path
    character(*),  intent(inout), optional :: cat_bias_inputs_file    
    
    ! -----------------------------------------------------------------
    
    integer :: N_args, iargc, i
    
    character(40) :: arg
    
    !external getarg, iargc
    
    ! -----------------------------------------------------------------
    
    N_args = iargc()
    
    i=0
    
    do while ( i < N_args )
       
       i = i+1
       
       call getarg(i,arg)
       
       if     ( trim(arg) == '-cat_bias_inputs_path' ) then
          i = i+1
          if (present(cat_bias_inputs_path))  &
               call getarg(i,cat_bias_inputs_path)
          
       elseif ( trim(arg) == '-cat_bias_inputs_file' ) then
          i = i+1
          if (present(cat_bias_inputs_file))  &
               call getarg(i,cat_bias_inputs_file)
          
       else
               
          i=i+1
          if (logit) write (logunit,*)                                    &
               'clsm_cat_bias_get_command_line(): IGNORING argument = ',  &
               trim(arg)
          
       endif
       
    end do
    
  end subroutine clsm_cat_bias_get_command_line
  
  ! ********************************************************************

  subroutine io_rstrt_cat_bias( action, work_path, exp_id, date_time, &
       model_dtstep, N_cat, N_catbias, cat_bias )
    
    ! read or write cat bias re-start file.
    !
    ! bias restart file contains all time-invariant and time-varying bias
    ! parameters
    !
    ! reichle, 18 Oct 2005
    ! reichle+draper, 27 Mar 2013 - revised cat_bias structure
    
    implicit none
    
    character,            intent(in) :: action     ! read or write
    
    character(*),       intent(in) :: work_path
    character(*),       intent(in) :: exp_id
    
    type(date_time_type), intent(in) :: date_time
    
    integer,              intent(in) :: model_dtstep, N_cat, N_catbias
    
    type(cat_progn_type), dimension(N_cat,N_catbias), intent(inout) :: cat_bias
    
    ! local variables
    
    integer        :: j, k, n, model_dtstep_tmp, N_cat_tmp, N_catbias_tmp
    
    character(300) :: filename
    
    character(40)  :: file_tag='catbias_ldas_rst', dir_name='rs', file_ext='.bin'
    
    character(len=*), parameter :: Iam = 'io_rstrt_cat_bias'

    ! --------------------------------------------------------------------
    
    select case (action)
       
    case ('r','R')
       
       filename = get_io_filename( work_path, exp_id,          &
            file_tag, date_time=date_time,                     &
            dir_name=dir_name, ens_id=-1, file_ext=file_ext )
       
       if (logit) write (logunit,*) 'Reading bias restart file ', trim(filename)
       
       open(10, file=filename, form='unformatted', status='old', &
            action='read')

       read (10) model_dtstep_tmp, N_cat_tmp, N_catbias_tmp
       
       if ( model_dtstep_tmp /= model_dtstep ) then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'inconsistent model_dtstep')
       end if
       
       if ( N_cat_tmp        /= N_cat        ) then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'inconsistent num of tiles')
       end if
       
       if ( N_catbias_tmp    /= N_catbias    ) then
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'inconsistent N_catbias')
       end if

       do j=1,N_catbias
          
          read (10) (cat_bias(n,j)%tc1,    n=1,N_cat)
          read (10) (cat_bias(n,j)%tc2,    n=1,N_cat)
          read (10) (cat_bias(n,j)%tc4,    n=1,N_cat)
          
          read (10) (cat_bias(n,j)%qa1,    n=1,N_cat)
          read (10) (cat_bias(n,j)%qa2,    n=1,N_cat)
          read (10) (cat_bias(n,j)%qa4,    n=1,N_cat)
          
          read (10) (cat_bias(n,j)%capac,  n=1,N_cat)
          
          read (10) (cat_bias(n,j)%catdef, n=1,N_cat)
          read (10) (cat_bias(n,j)%rzexc,  n=1,N_cat)
          read (10) (cat_bias(n,j)%srfexc, n=1,N_cat)
          
          do k=1,N_gt
             read (10) (cat_bias(n,j)%ght(k),  n=1,N_cat)
          end do
          
          do k=1,N_snow
             read (10) (cat_bias(n,j)%wesn(k), n=1,N_cat)
          end do
          do k=1,N_snow
             read (10) (cat_bias(n,j)%htsn(k), n=1,N_cat)
          end do
          do k=1,N_snow
             read (10) (cat_bias(n,j)%sndz(k), n=1,N_cat)
          end do
          
       end do
       
       
    case ('w','W')
       
       filename = get_io_filename( work_path, exp_id,         & 
            file_tag, date_time=date_time,                    &
            dir_name=dir_name, ens_id=-1, file_ext=file_ext )
       
       if (logit) write (logunit,*) 'Writing bias restart file ', trim(filename)
       
       open(10, file=filename, form='unformatted', status='unknown', &
            action='write')
       
       ! write header
       
       write (10) model_dtstep, N_cat, N_catbias
       
       ! write bias estimate
       
       do j=1,N_catbias
          
          write (10) (cat_bias(n,j)%tc1,    n=1,N_cat)
          write (10) (cat_bias(n,j)%tc2,    n=1,N_cat)
          write (10) (cat_bias(n,j)%tc4,    n=1,N_cat)
          
          write (10) (cat_bias(n,j)%qa1,    n=1,N_cat)
          write (10) (cat_bias(n,j)%qa2,    n=1,N_cat)
          write (10) (cat_bias(n,j)%qa4,    n=1,N_cat)
          
          write (10) (cat_bias(n,j)%capac,  n=1,N_cat)
          
          write (10) (cat_bias(n,j)%catdef, n=1,N_cat)
          write (10) (cat_bias(n,j)%rzexc,  n=1,N_cat)
          write (10) (cat_bias(n,j)%srfexc, n=1,N_cat)
          
          do k=1,N_gt
             write (10) (cat_bias(n,j)%ght(k),  n=1,N_cat)
          end do
          
          do k=1,N_snow
             write (10) (cat_bias(n,j)%wesn(k), n=1,N_cat)
          end do
          do k=1,N_snow
             write (10) (cat_bias(n,j)%htsn(k), n=1,N_cat)
          end do
          do k=1,N_snow
             write (10) (cat_bias(n,j)%sndz(k), n=1,N_cat)
          end do
          
       end do
       
       
    case default
              
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown action')

    end select
    
    close (10,status='keep')
    
  end subroutine io_rstrt_cat_bias
  
  ! *************************************************************

  subroutine init_cat_bias(                                           &
       work_path, exp_id, date_time, model_dtstep, N_cat, N_catbias,  &
       cat_bias )
    
    ! reichle, 18 Oct 2005
    
    implicit none

    character(*), intent(in) :: work_path
    character(*), intent(in) :: exp_id
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: model_dtstep, N_cat, N_catbias
    
    type(cat_progn_type), dimension(N_cat,N_catbias), intent(out) :: cat_bias

    ! local variables
     
    character(300) :: fname
    
    integer        :: j, n
    logical        :: fexists

    character(40)  :: file_tag = 'catbias_ldas_rst' 
    character(40)  :: dir_name = 'rs'
    character(40)  :: file_ext = '.bin'
    
    ! -------------------------------------------------------------
    
    ! try reading from restart file
    
    fname = get_io_filename( work_path, exp_id,             &
         file_tag, date_time=date_time,                     &
         dir_name=dir_name, ens_id=-1, file_ext=file_ext )
    
    inquire(file=fname, exist=fexists)
    
    if (fexists) then
       
       ! read bias restart file
       
       call io_rstrt_cat_bias(                                                  &
            'r', work_path, exp_id, date_time, model_dtstep, N_cat, N_catbias,  &
            cat_bias)
       
    else
       
       if (logit) then
          write (logunit,*) 'init_cat_bias(): restart file not found ', trim(fname)
          write (logunit,*) '                 initializing cat_bias to zero'
       end if
       
       do n=1,N_cat
          
          do j=1,N_catbias
             
             cat_bias(n,j) = 0.
             
          end do
          
       end do
       
    end if
    
  end subroutine init_cat_bias
  
  ! ********************************************************************
  
  subroutine check_cat_bias_inputs( update_type, cat_bias_param )
    
    ! Check cat bias param inputs against update_type:
    ! Make sure that the increments that are needed for bias estimation
    ! are computed with the selected "update_type".
    !
    ! reichle, 19 Oct 2005
    ! reichle, 11 Jan 2006
    
    implicit none
    
    integer, intent(in) :: update_type
    
    type(cat_bias_param_type), intent(in) :: cat_bias_param

    ! locals

    character(10) :: update_type_string

    character(len=*), parameter :: Iam = 'check_cat_bias_inputs'
    character(len=400) :: err_msg
        
    ! ------------------------

    write(update_type_string,'(i10)') update_type
    
    select_update_type: select case (update_type)
       
    case (1,2)
       
       if (                                                        & 
            (cat_bias_param%Nparam%tc1>0)        .or.              &
            (cat_bias_param%Nparam%tc2>0)        .or.              &
            (cat_bias_param%Nparam%tc4>0)        .or.              &   
            !(cat_bias_param%Nparam%srfexc>0)     .or.              &
            !(cat_bias_param%Nparam%rzexc >0)     .or.              &
            !(cat_bias_param%Nparam%catdef>0)     .or.              &
            (cat_bias_param%Nparam%qa1  >0)      .or.              &
            (cat_bias_param%Nparam%qa2  >0)      .or.              &
            (cat_bias_param%Nparam%qa4  >0)      .or.              &
            (cat_bias_param%Nparam%capac>0)      .or.              &
            any(cat_bias_param%Nparam%ght(:) >0) .or.              &
            any(cat_bias_param%Nparam%wesn(:)>0) .or.              &
            any(cat_bias_param%Nparam%htsn(:)>0) .or.              &
            any(cat_bias_param%Nparam%sndz(:)>0)         ) then
          
          err_msg = 'no increments computed for requested bias corr' // &
               'variables (update_type = ' // trim(update_type_string) // ')'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       else
          
          err_msg = 'untested bias corr requested for srfexc/rzexc/catdef' // &
               ' - are you sure? If yes, disable this abort and start over.'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end if
       
    case (3)
       
       if (                                                        &
            !(cat_bias_param%Nparam%tc1>0)        .or.              &
            !(cat_bias_param%Nparam%tc2>0)        .or.              &
            !(cat_bias_param%Nparam%tc4>0)        .or.              &   
            (cat_bias_param%Nparam%srfexc>0)     .or.              &
            (cat_bias_param%Nparam%rzexc >0)     .or.              &
            (cat_bias_param%Nparam%catdef>0)     .or.              &
            (cat_bias_param%Nparam%qa1  >0)      .or.              &
            (cat_bias_param%Nparam%qa2  >0)      .or.              &
            (cat_bias_param%Nparam%qa4  >0)      .or.              &
            (cat_bias_param%Nparam%capac>0)      .or.              &
            any(cat_bias_param%Nparam%ght(:) >0) .or.              &
            any(cat_bias_param%Nparam%wesn(:)>0) .or.              &
            any(cat_bias_param%Nparam%htsn(:)>0) .or.              &
            any(cat_bias_param%Nparam%sndz(:)>0)         ) then
          
          err_msg = 'no increments computed  for requested bias corr ' // &
               'variables (update_type = ' // trim(update_type_string)  // ')'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    case (4,5,7,9)
       
       if (                                                        & 
            !(cat_bias_param%Nparam%tc1>0)        .or.              &
            !(cat_bias_param%Nparam%tc2>0)        .or.              &
            !(cat_bias_param%Nparam%tc4>0)        .or.              &   
            (cat_bias_param%Nparam%srfexc>0)     .or.              &
            (cat_bias_param%Nparam%rzexc >0)     .or.              &
            (cat_bias_param%Nparam%catdef>0)     .or.              &
            (cat_bias_param%Nparam%qa1  >0)      .or.              &
            (cat_bias_param%Nparam%qa2  >0)      .or.              &
            (cat_bias_param%Nparam%qa4  >0)      .or.              &
            (cat_bias_param%Nparam%capac>0)      .or.              &
            !any(cat_bias_param%Nparam%ght(:) >0) .or.              &
            any(cat_bias_param%Nparam%ght(2:N_gt) >0) .or.         &
            any(cat_bias_param%Nparam%wesn(:)>0) .or.              &
            any(cat_bias_param%Nparam%htsn(:)>0) .or.              &
            any(cat_bias_param%Nparam%sndz(:)>0)         ) then
                    
          err_msg = 'no increments computed for requested bias corr ' // &
               'variables (update_type = ' // trim(update_type_string)  // ')'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    case (6,8)   
       
       if (                                                        &
            !(cat_bias_param%Nparam%tc1>0)        .or.              &
            !(cat_bias_param%Nparam%tc2>0)        .or.              &
            !(cat_bias_param%Nparam%tc4>0)        .or.              &   
            !(cat_bias_param%Nparam%srfexc>0)     .or.              &
            !(cat_bias_param%Nparam%rzexc >0)     .or.              &
            !(cat_bias_param%Nparam%catdef>0)     .or.              &
            (cat_bias_param%Nparam%qa1  >0)      .or.              &
            (cat_bias_param%Nparam%qa2  >0)      .or.              &
            (cat_bias_param%Nparam%qa4  >0)      .or.              &
            (cat_bias_param%Nparam%capac>0)      .or.              &
            !any(cat_bias_param%Nparam%ght(:) >0) .or.              &
            any(cat_bias_param%Nparam%ght(2:N_gt) >0) .or.         &
            any(cat_bias_param%Nparam%wesn(:)>0) .or.              &
            any(cat_bias_param%Nparam%htsn(:)>0) .or.              &
            any(cat_bias_param%Nparam%sndz(:)>0)         ) then
                         
          err_msg = 'no increments computed for requested bias corr ' // &
               'variables (update_type = ' // trim(update_type_string)  // ')'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    case (10)   
       
       if (                                                        &
            !(cat_bias_param%Nparam%tc1>0)        .or.              &
            !(cat_bias_param%Nparam%tc2>0)        .or.              &
            !(cat_bias_param%Nparam%tc4>0)        .or.              &   
            !(cat_bias_param%Nparam%srfexc>0)     .or.              &
            !(cat_bias_param%Nparam%rzexc >0)     .or.              &
            (cat_bias_param%Nparam%catdef>0)     .or.              &
            (cat_bias_param%Nparam%qa1  >0)      .or.              &
            (cat_bias_param%Nparam%qa2  >0)      .or.              &
            (cat_bias_param%Nparam%qa4  >0)      .or.              &
            (cat_bias_param%Nparam%capac>0)      .or.              &
            !any(cat_bias_param%Nparam%ght(:) >0) .or.              &
            any(cat_bias_param%Nparam%ght(2:N_gt) >0) .or.         &
            any(cat_bias_param%Nparam%wesn(:)>0) .or.              &
            any(cat_bias_param%Nparam%htsn(:)>0) .or.              &
            any(cat_bias_param%Nparam%sndz(:)>0)         ) then
                         
          err_msg = 'no increments computed for requested bias corr ' // &
               'variables (update_type = ' // trim(update_type_string)  // ')'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    case default     
              
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown update_type')
       
    end select select_update_type
    
  end subroutine check_cat_bias_inputs
    
  ! ********************************************************************
  
  subroutine read_cat_bias_inputs( work_path, exp_id, date_time,       &
       update_type, cat_bias_param, N_catbias )
    
    implicit none
    
    character(*),            intent(in)  :: work_path
    character(*),             intent(in)  :: exp_id
    
    type(date_time_type),      intent(in)  :: date_time
    
    integer,                   intent(in)  :: update_type

    type(cat_bias_param_type), intent(out) :: cat_bias_param
    
    integer,                   intent(out) :: N_catbias
    
    ! -----------
    
    character(300) :: fname
    character(200) :: cat_bias_inputs_path
    character( 40) :: cat_bias_inputs_file, dir_name, file_tag, file_ext

    logical :: file_exist
    
    ! -------------------------------------------------------------------
    
    namelist /cat_bias_inputs/         &
         cat_bias_param
    
    ! -------------------------------------------------------------------
    
    ! read default cat bias inputs file
    
    cat_bias_inputs_path = './'                                       ! set default 
    !call clsm_ensdrv_get_command_line(run_path=cat_bias_inputs_path)
    cat_bias_inputs_file = 'LDASsa_DEFAULT_inputs_catbias.nml'
      
    fname = trim(cat_bias_inputs_path) // '/' // trim(cat_bias_inputs_file)
    
    open (10, file=fname, delim='apostrophe', action='read', status='old')
    
    if (logit) write (logunit,*)
    if (logit) write (logunit,'(400A)') 'reading *default* cat bias inputs from ' // trim(fname)
    if (logit) write (logunit,*)
    
    read (10, nml=cat_bias_inputs)
    
    close(10,status='keep')
    
       
    ! Get name and path for special cat bias inputs file from
    ! command line (if present) 
    
   ! cat_bias_inputs_path = ''
   ! cat_bias_inputs_file = ''
    
   ! call clsm_cat_bias_get_command_line(                            &
   !      cat_bias_inputs_path=cat_bias_inputs_path,                 &
   !      cat_bias_inputs_file=cat_bias_inputs_file )
    
    cat_bias_inputs_file = 'LDASsa_SPECIAL_inputs_catbias.nml'
   ! if ( trim(cat_bias_inputs_path) /= ''  .and.                &
   !      trim(cat_bias_inputs_file) /= ''          ) then
       
       ! Read data from special cat bias inputs namelist file 
       
    fname = trim(cat_bias_inputs_path)//'/'//trim(cat_bias_inputs_file)

    inquire(file=fname, exist=file_exist)

    if (file_exist) then        

       open (10, file=fname, delim='apostrophe', action='read', status='old')
       
       if (logit) write (logunit,*)
       if (logit) write (logunit,'(400A)') 'reading *special* cat bias inputs from ' // trim(fname)
       if (logit) write (logunit,*)
       
       read (10, nml=cat_bias_inputs)
       
       close(10,status='keep')
       
    end if
    
    ! diagnose N_catbias = max(cat_bias_param%Nparam)
    
    N_catbias = 0
    
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%tc1      ) 
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%tc2      )
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%tc4      )
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%srfexc   )
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%rzexc    ) 
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%catdef   )
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%qa1      )
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%qa2      )
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%qa4      )
    N_catbias = max( N_catbias,        cat_bias_param%Nparam%capac    )
    N_catbias = max( N_catbias, maxval(cat_bias_param%Nparam%ght    ) )
    N_catbias = max( N_catbias, maxval(cat_bias_param%Nparam%wesn   ) )
    N_catbias = max( N_catbias, maxval(cat_bias_param%Nparam%htsn   ) )
    N_catbias = max( N_catbias, maxval(cat_bias_param%Nparam%sndz   ) )


    ! check for consistency with update_type
    
    if (N_catbias>0)  call check_cat_bias_inputs( update_type, cat_bias_param )
    
   
    ! -------------------------------------------------------------
    !
    ! echo variables of namelist cat_bias_inputs
    
    if (logit) write (logunit,*) 'cat bias inputs are:'
    if (logit) write (logunit,*)
    if (logit) write (logunit, nml=cat_bias_inputs) 
    if (logit) write (logunit,*)
    
    ! -------------------------------------------------------------
    !
    ! save cat bias inputs into *catbias_inputs.nml file
    
    dir_name = 'rc_out'
    file_tag = 'ldas_catbias_inputs'
    file_ext = '.nml'
    
    fname = get_io_filename( work_path, exp_id, file_tag, date_time=date_time, &
         dir_name=dir_name, file_ext=file_ext )
    
    if (logit) write (logunit,'(400A)') 'writing cat bias inputs to ' // trim(fname)
    
    open (10, file=fname, status='unknown', action='write', delim='apostrophe')
    
    write(10, nml=cat_bias_inputs)
    
    close(10, status='keep')
    
  end subroutine read_cat_bias_inputs
  
  ! ********************************************************************
  
  !subroutine cat_bias_calcs_corr( )

  ! This subroutine was a wrapper for subroutine cat_bias_corr() with 
  ! an added call to subroutine recompute_diagnostics().
  ! The call to recompute_diagnostics() has been moved, making the wrapper
  ! obsolete
  ! -reichle+csdraper, 30 Oct 2013
      
  !end subroutine cat_bias_calcs_corr
  
  ! ********************************************************************
  
  subroutine cat_bias_calcs_update( date_time, model_dtstep_real,            &
       N_cat, N_ens, N_catbias, cat_progn_incr, fresh_incr, cat_bias_param,  &
       cat_bias )
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
     
    real,    intent(in) :: model_dtstep_real
    integer, intent(in) :: N_cat, N_ens, N_catbias
    
    type(cat_progn_type), dimension(N_cat,N_ens),     intent(in)    :: cat_progn_incr
    
    logical, intent(in) :: fresh_incr
    
    type(cat_bias_param_type),                        intent(in)    :: cat_bias_param
    
    type(cat_progn_type), dimension(N_cat,N_catbias), intent(inout) :: cat_bias
    
    ! local variables
    
    type(cat_progn_type), dimension(N_cat) :: cat_progn_incr_vec
    
    ! --------------------------------------------------------------
    
    if (fresh_incr) then
       
       ! update bias parameters from ensemble average increments
       
       call get_cat_progn_ens_avg(N_cat, N_ens, cat_progn_incr,       &
            cat_progn_incr_vec)
       
       call cat_bias_update(date_time, model_dtstep_real, N_cat,      &
            N_catbias, cat_progn_incr_vec, cat_bias_param, cat_bias )
       
    end if
    
  end subroutine cat_bias_calcs_update
  
  ! ********************************************************************
  
  subroutine cat_bias_corr( date_time, dtstep, N_cat, N_ens, &
       N_catbias, cat_bias_param, cat_bias, cat_progn )
    
    ! apply bias correction to cat_progn, relax bias parameters)
    !
    ! 
    ! reichle, 19 Oct 2005
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    real,    intent(in) :: dtstep
    integer, intent(in) :: N_cat, N_ens, N_catbias   

    type(cat_bias_param_type), intent(in) :: cat_bias_param
    
    type(cat_progn_type), dimension(N_cat,N_catbias), intent(inout) :: cat_bias
    type(cat_progn_type), dimension(N_cat,N_ens),     intent(inout) :: cat_progn
    
    ! local variables
    
    integer :: n
    
    real, dimension(N_catbias) :: const_param

    character(len=*), parameter :: Iam = 'cat_bias_corr'
    character(len=400) :: err_msg
    
    ! -----------------------------------------------------------
    
    call bias_get_const_param( date_time, N_catbias, const_param )
    
    do n=1,N_cat
       
       ! tc1, tc2, tc4
       
       if (cat_bias_param%Nparam%tc1>0)                                                    &
            call bias_corr_helper( date_time, dtstep, N_ens, cat_bias_param%Nparam%tc1,    &
            const_param, cat_bias_param%trelax%tc1, cat_progn(n,:)%tc1,                    &
            cat_bias(n,1:cat_bias_param%Nparam%tc1)%tc1 )
       
       if (cat_bias_param%Nparam%tc2>0)                                                    &
            call bias_corr_helper( date_time, dtstep, N_ens, cat_bias_param%Nparam%tc2,    &
            const_param, cat_bias_param%trelax%tc2, cat_progn(n,:)%tc2,                    &
            cat_bias(n,1:cat_bias_param%Nparam%tc2)%tc2 )
       
       if (cat_bias_param%Nparam%tc4>0)                                                    &
            call bias_corr_helper( date_time, dtstep, N_ens, cat_bias_param%Nparam%tc4,    &
            const_param, cat_bias_param%trelax%tc4, cat_progn(n,:)%tc4,                    &
            cat_bias(n,1:cat_bias_param%Nparam%tc4)%tc4 )
       
       if (cat_bias_param%Nparam%ght(1)>0)                                                 &
            call bias_corr_helper( date_time, dtstep, N_ens, cat_bias_param%Nparam%ght(1), &
            const_param, cat_bias_param%trelax%ght(1), cat_progn(n,:)%ght(1),              &
            cat_bias(n,1:cat_bias_param%Nparam%ght(1))%ght(1) )
       
       
       if ( (cat_bias_param%Nparam%qa1   >0)     .or.              &
            (cat_bias_param%Nparam%qa2   >0)     .or.              &
            (cat_bias_param%Nparam%qa4   >0)     .or.              &
            (cat_bias_param%Nparam%capac >0)     .or.              &
            (cat_bias_param%Nparam%srfexc>0)     .or.              &
            (cat_bias_param%Nparam%rzexc >0)     .or.              &
            (cat_bias_param%Nparam%catdef>0)     .or.              &
            any(cat_bias_param%Nparam%ght(2:N_gt) >0) .or.              &
            any(cat_bias_param%Nparam%wesn(:)>0) .or.              &
            any(cat_bias_param%Nparam%htsn(:)>0) .or.              &
            any(cat_bias_param%Nparam%sndz(:)>0)         ) then
          
          err_msg = 'must add fields to source code'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    end do
    
  end subroutine cat_bias_corr
  
  ! ********************************************************************
  
  subroutine bias_corr_helper( date_time, dtstep, N_ens, Nparam, & 
       const_param, trelax, cat_progn_field, bias_param )
    
    implicit none
    
    ! diagnose/apply bias flux and relax bias parameters for a single 
    ! tile and a single field of cat_progn
    !
    ! reichle, 17 Oct 2005
    !
    ! reichle, 18 Aug 2008 -- added new "time-of-day" bias option
    
    type(date_time_type), intent(in) :: date_time
    
    real,    intent(in) :: dtstep, trelax
    
    integer, intent(in) :: N_ens, Nparam
    
    real, dimension(Nparam), intent(in) :: const_param
    
    real, dimension(N_ens), intent(inout) :: cat_progn_field
    
    real, dimension(Nparam), intent(inout) :: bias_param
    
    ! local variables
    
    real :: tmpflux, relax_fac
    
    integer :: i, ind_start, ind_end
    
    real, dimension(Nparam) :: const_param_tmp
    
    ! -----------------------------------------------------------------    
    
    if (Nparam > 0) then
       
       call bias_options_helper( date_time, Nparam, const_param, &
            ind_start, ind_end, const_param_tmp )
       
       ! diagnose bias flux from bias parameters
       
       tmpflux = 0.       
       
       do i=ind_start,ind_end
          
          tmpflux = tmpflux + bias_param(i)*const_param_tmp(i)
          
       end do
       
       ! apply bias flux
       
       cat_progn_field(:) = cat_progn_field(:) - tmpflux*dtstep
       
       ! relax bias parameters
       
       relax_fac = (1. - dtstep/trelax)
       
       do i=ind_start,ind_end
          
          bias_param(i) = bias_param(i) * relax_fac
          
       end do
       
    end if
    
  end subroutine bias_corr_helper
  
  ! ********************************************************************
       
  subroutine cat_bias_update( date_time, dtstep,             &
       N_cat, N_catbias, cat_progn_incr, cat_bias_param,     &
       cat_bias)
    
    ! update bias parameters from ens avg assimilation increments
    !
    ! reichle, 19 Oct 2005
        
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    real,    intent(in) :: dtstep
    integer, intent(in) :: N_cat, N_catbias
    
    type(cat_progn_type), dimension(N_cat),           intent(in)    :: cat_progn_incr
    
    type(cat_bias_param_type),                        intent(in)    :: cat_bias_param

    type(cat_progn_type), dimension(N_cat,N_catbias), intent(inout) :: cat_bias
        
    ! local variables
    
    integer :: n
    
    real, dimension(N_catbias) :: const_param

    character(len=*), parameter :: Iam = 'cat_bias_update'
    character(len=400) :: err_msg
    
    ! -------------------------------------------------------------------- 

    call bias_get_const_param( date_time, N_catbias, const_param )
    
    do n=1,N_cat
       
       ! tc1, tc2, tc4
       
       if (cat_bias_param%Nparam%tc1>0)                                                   &
            call bias_update_helper(date_time, dtstep, cat_bias_param%Nparam%tc1,         &
            const_param, cat_bias_param%tconst%tc1, cat_progn_incr(n)%tc1,                &
            cat_bias(n,1:cat_bias_param%Nparam%tc1)%tc1 )

       if (cat_bias_param%Nparam%tc2>0)                                                   &
            call bias_update_helper(date_time, dtstep, cat_bias_param%Nparam%tc2,         &
            const_param, cat_bias_param%tconst%tc2, cat_progn_incr(n)%tc2,                &
            cat_bias(n,1:cat_bias_param%Nparam%tc2)%tc2 )
       
       if (cat_bias_param%Nparam%tc4>0)                                                   &
            call bias_update_helper(date_time, dtstep, cat_bias_param%Nparam%tc4,         &
            const_param, cat_bias_param%tconst%tc4, cat_progn_incr(n)%tc4,                &
            cat_bias(n,1:cat_bias_param%Nparam%tc4)%tc4 )
       
       if (cat_bias_param%Nparam%ght(1)>0)                                                &
            call bias_update_helper(date_time, dtstep, cat_bias_param%Nparam%ght(1),      &
            const_param, cat_bias_param%tconst%ght(1), cat_progn_incr(n)%ght(1),          &
            cat_bias(n,1:cat_bias_param%Nparam%ght(1))%ght(1) )
       
       if ( (cat_bias_param%Nparam%qa1   >0)     .or.              &
            (cat_bias_param%Nparam%qa2   >0)     .or.              &
            (cat_bias_param%Nparam%qa4   >0)     .or.              &
            (cat_bias_param%Nparam%capac >0)     .or.              &
            (cat_bias_param%Nparam%srfexc>0)     .or.              &
            (cat_bias_param%Nparam%rzexc >0)     .or.              &
            (cat_bias_param%Nparam%catdef>0)     .or.              &
            any(cat_bias_param%Nparam%ght(2:N_gt) >0) .or.              &
            any(cat_bias_param%Nparam%wesn(:)>0) .or.              &
            any(cat_bias_param%Nparam%htsn(:)>0) .or.              &
            any(cat_bias_param%Nparam%sndz(:)>0)         ) then
          
          err_msg = 'must add fields to source code'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    end do


  end subroutine cat_bias_update
  
  ! ********************************************************************

  subroutine bias_update_helper( date_time, dtstep, Nparam, const_param, tconst, &
       increment, bias_param )
    
    implicit none
    
    ! update bias parameters for a single tile and a single field of cat_progn
    !
    ! reichle, 17 Oct 2005
    !
    ! reichle, 18 Aug 2008 -- added new "time-of-day" bias option
    
    type(date_time_type), intent(in) :: date_time
    
    real, intent(in) :: dtstep, tconst, increment
    
    integer, intent(in) :: Nparam
    
    real, dimension(Nparam), intent(in) :: const_param
    
    real, dimension(Nparam), intent(inout) :: bias_param
    
    ! local variables

    real, dimension(Nparam) :: const_param_tmp
    
    real :: tmp_real
    
    integer :: i, ind_start, ind_end
    
    ! -------------------------------------------------------------
    
    if (Nparam > 0) then
       
       call bias_options_helper( date_time, Nparam, const_param, &
            ind_start, ind_end, const_param_tmp )
       
       ! update bias parameters
       
       tmp_real = - tconst * increment / dtstep
       
       do i=ind_start,ind_end 
          
          bias_param(i) = bias_param(i) + tmp_real*const_param_tmp(i)
          
       end do
       
    end if
    
  end subroutine bias_update_helper
  
  ! ***************************************************************************
  
  subroutine bias_options_helper( date_time, Nparam, const_param_in, &
       ind_start, ind_end, const_param_out )
    
    implicit none
    
    type(date_time_type),                    intent(in)  :: date_time
    integer,                                 intent(in)  :: Nparam     
    real,                 dimension(Nparam), intent(in)  :: const_param_in
    
    integer,                                 intent(out) :: ind_start, ind_end    
    real,                 dimension(Nparam), intent(out) :: const_param_out
    
    ! local variables
    
    !integer :: bias_time_of_day_index
    
    ! -------------------------------------------------------------
    
    select case (Nparam)
       
    case (1,3,5)          ! constant or sine/cosine (semi-)diurnal bias corr
       
       ind_start = 1
       ind_end   = Nparam
       
       const_param_out = const_param_in
       
    case (2,4,8)          ! separate "time-of-day" bias correction
       
       ind_start = bias_time_of_day_index( date_time, Nparam )
       ind_end   = ind_start
       
       const_param_out = 1.
       
    case default
       
       write (*,*) 'bias_options_helper(): unknown Nparam = ', Nparam
       
    end select
    
  end subroutine bias_options_helper
  
  ! *************************************************************************
  
  subroutine bias_get_const_param( date_time, N_const, const_param )
    
    ! compute time-of-day dependent sine/cosine variables for bias calculations
    !
    ! reichle, 18 Aug 2008
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    integer,              intent(in) :: N_const
    
    real, dimension(N_const), intent(out) :: const_param
    
    ! local variables
    
    real, parameter :: MY_PI = 3.14159265
    
    real, parameter :: omega1 = 2.*MY_PI/86400.
    real, parameter :: omega2 = 2.*MY_PI/43200.
    
    real :: secs_in_day, om1_t, om2_t
    
    ! ------------------------------------------------------------------

    const_param = nodata_generic
    
    secs_in_day = real(date_time%hour)*3600. + real(date_time%min)*60. &
         + real(date_time%sec)
    
    om1_t = secs_in_day*omega1
    om2_t = secs_in_day*omega2
    

    if (N_const>=1)  const_param(1) = 1. 
    
    if (N_const>=3) then
       
       const_param(2) = cos(om1_t)
       const_param(3) = sin(om1_t)
       
    end if

    if (N_const>=5) then
       
       const_param(4) = cos(om2_t)
       const_param(5) = sin(om2_t)
       
    end if

  end subroutine bias_get_const_param
  
  ! -------------------------------------------------------------------
  !
  ! OBS BIAS ROUTINES 
  !
  ! -------------------------------------------------------------------
  
  subroutine initialize_obs_bias( N_catf, N_obs_param, N_obsbias_max, work_path, & 
       exp_id, date_time, N_catl, numprocs, N_catl_vec, low_ind, obs_bias) 
    
    ! initialize obs_bias for use with observations bias corrections, 
    ! and handle mpi 
    ! 
    ! draper, Aug 29 2013. 
    
    implicit none		
    
    integer,              intent(in)                         :: N_catf, N_obs_param
    integer,              intent(in)                         :: N_obsbias_max
    
    character(*),       intent(in)                         :: work_path
    character(*),        intent(in)                         :: exp_id
    
    type(date_time_type), intent(in)                         :: date_time
    integer,              intent(in)                         :: N_catl, numprocs 
    integer,              intent(in),    dimension(numprocs) :: N_catl_vec, low_ind 
    
    type(obs_bias_type),  intent(inout), dimension(N_catl,N_obs_param,N_obsbias_max) :: & 
         obs_bias
    
    !local variables
    
    type(obs_bias_type),  dimension(:,:,:), allocatable :: obs_bias_f 

    integer :: i,j

    ! ------------------------------------------------------------------
    
    if (master_proc) then 
       
       allocate(obs_bias_f(N_catf,N_obs_param,N_obsbias_max))
       
       call init_obs_bias(                                                     &
	    work_path, exp_id, date_time, N_catf, N_obs_param, N_obsbias_max,  &
	    obs_bias_f )
       
    end if
	
#ifdef LDAS_MPI
    
    do i=1,N_obs_param
       do j=1,N_obsbias_max

          call MPI_SCATTERV(                                                &
               obs_bias_f(:,i,j), N_catl_vec, low_ind-1, MPI_obs_bias_type, &
               obs_bias(  :,i,j), N_catl,                MPI_obs_bias_type, &
	       0, mpicomm, mpierr )
          
       end do
    end do
    
#else
       
    obs_bias = obs_bias_f
    
#endif
    
    if (master_proc)  deallocate(obs_bias_f)
    
  end subroutine initialize_obs_bias
  
  ! *************************************************************************

  subroutine init_obs_bias(                                              &
       work_path, exp_id, date_time, N_cat, N_obs_param, N_obsbias_max,  &
       obs_bias )
    
    ! draper, 4 April 2013
    ! based on init_cat_bias
    
    implicit none
    
    character(*),       intent(in) :: work_path
    character(*),       intent(in) :: exp_id
    
    type(date_time_type), intent(in) :: date_time
    
    integer,              intent(in) :: N_cat, N_obs_param,N_obsbias_max
    
    type(obs_bias_type),  intent(out), dimension(N_cat,N_obs_param,N_obsbias_max) :: &
         obs_bias
    
    ! local variables
     
    character(300) :: fname
    
    logical        :: fexists
    
    character(40)  :: file_tag = 'obsbias_ldas_rst' 
    character(40)  :: dir_name = 'rs'
    character(40)  :: file_ext = '.bin'
    
    ! -------------------------------------------------------------
    
    ! try reading from restart file
    
    fname = get_io_filename( work_path, exp_id,             &
         file_tag, date_time=date_time,                     &
         dir_name=dir_name, ens_id=-1, file_ext=file_ext )

    inquire(file=fname, exist=fexists)
    
    if (fexists) then
       
       ! read bias restart file
       
       call io_rstrt_obs_bias('r', work_path, exp_id, date_time, & 
            N_cat, N_obs_param, N_obsbias_max, obs_bias )
       
    else
       
       if (logit) then
          write (logunit,*) 'init_obs_bias(): restart file not found ', trim(fname)
       end if
       
       obs_bias(:,:,:)%bias      = 0.0
       
       obs_bias(:,:,:)%tcount(1) = max_tcount 
       obs_bias(:,:,:)%tcount(2) = max_tcount
       
    end if
    
  end subroutine init_obs_bias
  
  ! ********************************************************************
  
  subroutine io_rstrt_obs_bias( action, work_path, exp_id, date_time, &
       N_cat, N_obs_param, N_obsbias_max, obs_bias )
    
    ! read or write obs bias re-start file.
    !
    ! bias restart file contains all time-invariant and time-varying bias
    ! parameters
    !
    ! reichle, 18 Oct 2005
    ! reichle+draper, 27 Mar 2013 - revised cat_bias structure
    
    implicit none
    
    character,            intent(in) :: action     ! read or write
    
    character(*),       intent(in) :: work_path
    character(*),       intent(in) :: exp_id
    
    type(date_time_type), intent(in) :: date_time
    
    integer,              intent(in) ::  N_cat, N_obs_param, N_obsbias_max
    
    type(obs_bias_type),  intent(inout), dimension(N_cat,N_obs_param,N_obsbias_max) :: &
         obs_bias
    
    ! local variables
    
    integer        :: i, j, k, n, N_cat_tmp, N_obs_param_tmp, N_obsbias_max_tmp
    
    character(300) :: filename
    
    character(40)  :: file_tag='obsbias_ldas_rst', dir_name='rs', file_ext='.bin'

    character(len=*), parameter :: Iam = 'io_rstrt_obs_bias'
    character(len=400) :: err_msg
    
    ! --------------------------------------------------------------------
    
    select case (action)
       
    case ('r','R')
       
       filename = get_io_filename( work_path, exp_id,          &
            file_tag, date_time=date_time,                     &
            dir_name=dir_name, ens_id=-1, file_ext=file_ext )
       
       if (logit) write (logunit,*) 'Reading obs bias restart file ', trim(filename)
       
       open(10, file=filename, form='unformatted', status='old',                  &
            action='read')
       
       read (10) N_cat_tmp, N_obs_param_tmp, N_obsbias_max_tmp
       
       if ( N_cat_tmp         /= N_cat        ) then
          err_msg = 'inconsistent number of tiles'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       if ( N_obs_param_tmp   /= N_obs_param    ) then
          err_msg = 'inconsistent N_obs_param'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       if ( N_obsbias_max_tmp /= N_obsbias_max    ) then
          err_msg = 'inconsistent N_catbias'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       do j=1,N_obsbias_max
          do i=1,N_obs_param
             
             read (10) (obs_bias(n,i,j)%bias, n=1,N_cat)
             
             do k=1,2
                read (10) (obs_bias(n,i,j)%tcount(k), n=1,N_cat)
             end do

          end do
       end do
       
    case ('w','W')
       
       filename = get_io_filename( work_path, exp_id,          & 
            file_tag, date_time=date_time,                     &
            dir_name=dir_name, ens_id=-1, file_ext=file_ext )
       
       if (logit) write (logunit,*) 'Writing obs bias restart file ', trim(filename)
       
       open(10, file=filename, form='unformatted', status='unknown',              &
            action='write')
       
       ! write header
       
       write (10)  N_cat, N_obs_param, N_obsbias_max
       
       do j=1,N_obsbias_max
          do i=1,N_obs_param
             
             write (10) (obs_bias(n,i,j)%bias, n=1,N_cat)

             do k=1,2
                write (10) (obs_bias(n,i,j)%tcount(k), n=1,N_cat)
             end do

          end do
       end do
       
    case default
              
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown action')
       
    end select
    
    close (10,status='keep')
    
  end subroutine io_rstrt_obs_bias
  
  ! ********************************************************************
  
  subroutine output_obs_bias(N_obs_param, N_obsbias_max, N_catl, N_catf,      &
       numprocs, N_catl_vec, low_ind, work_path, exp_id, date_time, obs_bias)
    
    ! output the obs_bias restarts
    ! 
    ! draper, Aug 29 2013
    
    implicit none 
    
    integer,              intent(in) :: N_obs_param, N_obsbias_max
    integer,              intent(in) :: N_catl, N_catf, numprocs
    
    integer,              intent(in), dimension(numprocs) :: N_catl_vec, low_ind 
    
    character(*),       intent(in) :: work_path
    character(*),       intent(in) :: exp_id
    type(date_time_type), intent(in) :: date_time
    
    type(obs_bias_type),  intent(in), dimension(N_catl,N_obs_param,N_obsbias_max) :: &
         obs_bias
    
    ! local variables
    
    type(obs_bias_type),       dimension(:,:,:), allocatable :: obs_bias_f 
    
    integer :: i,j

    if (master_proc) allocate(obs_bias_f(N_catf,N_obs_param, N_obsbias_max))
    
#ifdef LDAS_MPI
    
    do i=1,N_obs_param	
       do j=1,N_obsbias_max
          
          call MPI_GATHERV(                                                    &
               obs_bias(:,i,j), N_catl, MPI_obs_bias_type,                     &
               obs_bias_f(  :,i,j), N_catl_vec, low_ind-1, MPI_obs_bias_type,  &
               0, mpicomm, mpierr )
          
       end do
    end do
    
    call MPI_BARRIER( mpicomm, mpierr )
    
#else                 
    
    obs_bias_f = obs_bias
    
#endif
    
    if (master_proc) then
       
       call io_rstrt_obs_bias(                         &
            'w', work_path, exp_id, date_time, N_catf, &
            N_obs_param, N_obsbias_max, obs_bias_f )
       
       deallocate(obs_bias_f)		 
       
    end if
    
  end subroutine output_obs_bias
  
  ! ********************************************************************
  
  subroutine obs_bias_corr_obs(date_time, N_catl, N_catf, N_obsl, N_obs_param, & 
       N_obsbias_max, f2l, obs_param, obs_bias, Observations, obsbias_ok) 

    ! correct observations to remove the obs_bias (obs = obs minus obs_bias)
    ! set obsbias_ok flag to indicate whether have confidence in obs_bias estimate

    ! draper, Apr 2013. 

    implicit none

    type(date_time_type), intent(in)                            :: date_time  

    integer,              intent(in)                            :: N_catl, N_catf, N_obsl
    integer,              intent(in)                            :: N_obs_param
    integer,              intent(in)                            :: N_obsbias_max

    integer,              intent(in),    dimension(N_catf)      :: f2l

    type(obs_param_type), intent(in),    dimension(N_obs_param) :: obs_param  

    type(obs_bias_type),  intent(in),    dimension(N_catl, N_obs_param, N_obsbias_max) :: &
         obs_bias 

    type(obs_type),       intent(inout), dimension(N_obsl)       :: Observations

    logical,              intent(inout), dimension(N_obsl)       :: obsbias_ok

    ! local variables

    integer                         :: i, ind_catl, ind_spec, tcount2

    integer, dimension(N_obs_param) :: indv_time  

    ! ----------------------------------------------------------

    ! get species-dependent time of day index

    do i=1,N_obs_param

       indv_time(i)=bias_time_of_day_index( date_time,obs_param(i)%bias_Npar )

    end do
    
    do i=1,N_obsl

       ind_spec = Observations(i)%species

       if (obs_param(ind_spec)%bias_Npar > 0) then

          ind_catl = f2l(Observations(i)%tilenum)

          ! correct Observation(i) for bias

          Observations(i)%obs = (                                         &
               Observations(i)%obs -                                      &
               obs_bias(ind_catl, ind_spec, indv_time(ind_spec))%bias )

          if (Observations(i)%assim) then
             
             ! determine of obs bias correction is good enough for use of
             ! obs in state update

             tcount2 = obs_bias(ind_catl,ind_spec,indv_time(ind_spec))%tcount(2)
             
             if ( tcount2 < obs_param(ind_spec)%bias_tcut ) then
                
		! obsbias is good, keep assim flag true, set obsbias_ok
                
		obsbias_ok(i) = .TRUE.
  
             else
                
                ! obsbias is not good: switch assim flag to false, leave obsbias_ok as is
                
                Observations(i)%assim = .FALSE.
                
             end if
             
          end if
          
       end if
       
    end do
    
  end subroutine  obs_bias_corr_obs
  
  ! ********************************************************************
  
  subroutine obs_bias_upd_bias_and_Obs(                            &
       date_time, N_catl, N_catf, N_obs,                           &
       N_ens,  N_obs_param, N_obsbias_max, f2l, obs_param,	   & 	 
       Obs_pred, obs_bias, Observations )
    
    ! calculate the bias increment, and use to update obs_bias and the Observations
    !
    ! delta_b = b+ - b_= lambda(y - <Hx-> - b-) 
    !
    ! lambda  = 1 - exp( -Delta_time_since_last_ob/trelax ) 
    !
    ! draper+reichle, Sep 2013. 
    
    implicit none
    
    ! -----------------------------------------------------------
    
    type(date_time_type), intent(in)                             :: date_time  
    
    integer,              intent(in)                             :: N_catl, N_catf
    integer,              intent(in)                             :: N_obs, N_ens 
    integer,              intent(in)                             :: N_obs_param
    integer,              intent(in)                             :: N_obsbias_max
    
    integer,              intent(in),    dimension(N_catf)       :: f2l
    
    type(obs_param_type), intent(in),    dimension(N_obs_param)  :: obs_param   
    
    real,                 intent(in),    dimension(N_obs, N_ens) :: Obs_pred
    
    type(obs_bias_type),  intent(inout), dimension(N_catl, N_obs_param, N_obsbias_max) :: &
         obs_bias

    type(obs_type),       intent(inout), dimension(N_obs)        :: Observations

    ! local variables 
    
    real                            :: lambda, bias_incr

    integer                         :: i, ind_catl, ind_spec, tcount1, trel

    integer, dimension(N_obs_param) :: indv_time 
    
    ! ---------------------------------------------------------------------

    ! get species-dependent time of day index

    do i=1,N_obs_param
       
       indv_time(i) = bias_time_of_day_index( date_time, obs_param(i)%bias_Npar )
       
    end do

    ! update obs_bias
    
    do i=1,N_obs
       
       ind_spec = Observations(i)%species
       
       if (obs_param(ind_spec)%bias_Npar > 0) then
          
	  ind_catl = f2l(Observations(i)%tilenum)
    
          tcount1  = obs_bias(ind_catl,ind_spec,indv_time(ind_spec))%tcount(1)
          
          trel     = obs_param(ind_spec)%bias_trel
          
          lambda   = (1. - exp(-real(tcount1)/real(trel)))
          
          ! get the bias increment
          
          bias_incr =                                                               & 
               lambda*(Observations(i)%obs - sum(Obs_pred(i,1:N_ens))/real(N_ens))
          
          ! update the bias with the bias increment 

          obs_bias(ind_catl,ind_spec,indv_time(ind_spec))%bias =                    & 
               obs_bias(ind_catl,ind_spec,indv_time(ind_spec))%bias + bias_incr
          
          ! update the obs  with the bias increment
                    
          Observations(i)%obs = Observations(i)%obs - bias_incr
          
       end if

    end do
    
  end subroutine obs_bias_upd_bias_and_Obs
  
  ! ********************************************************************
  
  subroutine obs_bias_upd_tcount(date_time, dtstep, N_catf, N_catl, N_Obs,   &
       N_obs_param, N_obsbias_max, f2l, obs_param, Observations, obs_bias  )  
    ! 
    ! update counter recording the time since the last two observations
    ! (times are relevant at start of next assim cycle) 
    !
    ! draper, apr 2013 
    
    type(date_time_type), intent(in)                         :: date_time

    integer,              intent(in)                         :: dtstep    
    integer,              intent(in)                         :: N_catf, N_catl, N_obs
    integer,              intent(in)                         :: N_obs_param, N_obsbias_max

    integer,              intent(in), dimension(N_catf)      :: f2l
    
    type(obs_param_type), intent(in), dimension(N_obs_param) :: obs_param  
    type(obs_type),       intent(in), dimension(N_obs)       :: Observations
    
    type(obs_bias_type),  intent(inout), dimension(N_catl, N_obs_param, N_obsbias_max) :: &
         obs_bias
    
    ! local variables
    
    integer                         :: ind_catl, ind_spec, i, j, k, t
    integer, dimension(N_obs_param) :: indv_time 

    ! ------------------------------------------------------------------
    !
    ! add tstep to tcount 

    do i=1,N_catl
       do j=1,N_obs_param
          do k=1,obs_param(j)%bias_Npar
             do t=1,2
                
                obs_bias(i,j,k)%tcount(t) =                              &
                     min( obs_bias(i,j,k)%tcount(t)+dtstep, max_tcount )
                
             end do
          end do
       end do
    end do

    ! shift bias tcounts if had an observation 

    do i=1, N_obs_param
       indv_time(i)=bias_time_of_day_index( date_time,obs_param(i)%bias_Npar )
    end do

    do i=1, N_obs

       ind_spec=Observations(i)%species
       
       if (obs_param(ind_spec)%bias_Npar > 0) then
          
          ind_catl = f2l(Observations(i)%tilenum)
          
          ! tcount(2) = tcount(1)
          
          obs_bias(ind_catl,ind_spec,indv_time(ind_spec))%tcount(2) =     & 
               obs_bias(ind_catl,ind_spec,indv_time(ind_spec))%tcount(1)
          
          ! reinitialize tcount(1)
          
          obs_bias(ind_catl,ind_spec,indv_time(ind_spec))%tcount(1) = dtstep

       end if

    end do
    
  end subroutine obs_bias_upd_tcount
  
  ! -------------------------------------------------------------------
  !
  ! SHARED (OBS/CAT BIAS) ROUTINES 
  !
  !-------------------------------------------------------------------
  
  integer function bias_time_of_day_index( date_time, Nparam )
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    integer,              intent(in) :: Nparam             
    
    ! local variables
    
    integer :: tmpint, dtstep_tmp
    
    ! -----------------------------------------------
    
    dtstep_tmp = 86400/Nparam
    
    ! seconds-in-day 
    
    tmpint = date_time%hour*3600 + date_time%min*60 + date_time%sec 
    
    ! add half of dtstep_assim (to center interval around nominal time)
    
    tmpint = tmpint + dtstep_tmp/2
    
    ! modulus calculation (yields tmpint ranging from 1 to N+1)
    
    tmpint = tmpint/dtstep_tmp + 1
    
    ! fix "N+1" with another modulus calculation
    
    bias_time_of_day_index = mod( tmpint-1, 86400/dtstep_tmp) + 1
    
  end function bias_time_of_day_index
  
  ! ******************************************************************************
  
end module clsm_bias_routines

! ----------------------------------------------------------------    
! ----------------------------------------------------------------    

! test programs

#if 0

program test_bias_time_helper
  
  ! ifort date_time_util.o tmp.F90
  
  use date_time_util
  use leap_year
  
  implicit none

  type(date_time_type) :: date_time
  
  integer :: dtstep_assim(5), time_ind_bias, dtstep, n, k, tmpintvec(5), ttt
  
  real, dimension(8) :: const_param
  
  date_time%year  = 2008
  date_time%month =    8
  date_time%day   =   18
  date_time%hour  =   20
  date_time%min   =    0
  date_time%sec   =    0
  
  dtstep_assim(1) = 86400
  dtstep_assim(2) = 43200
  dtstep_assim(3) = 21600
  dtstep_assim(4) = 10800
  dtstep_assim(5) =  3600
  
  dtstep = 1200
  
  do ttt=1,100
     
     call augment_date_time(dtstep, date_time)
     
     do n=1,5
        
        call bias_time_helper( date_time, dtstep_assim(n),          &
             const_param, time_ind_bias )
        
        tmpintvec(n)=time_ind_bias
        
     end do
     
     write (*,'(7i5)') date_time%hour, date_time%min, tmpintvec(1:5)
     write (999,'(7i5)') date_time%hour, date_time%min, tmpintvec(1:5)
     
  end do
  
end program test_bias_time_helper

#endif 

! ====================== EOF ==============================================
