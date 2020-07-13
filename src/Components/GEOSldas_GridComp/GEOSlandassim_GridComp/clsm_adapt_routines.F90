! this file contains a collection of subroutines for adaptive features in the 
! the Ensemble Kalman filter of the Catchment model off-line driver
!
! reichle, 14 Dec 2006
! reichle,  1 Feb 2007 - major re-design for separate tuning of P and R
! reichle, 11 Apr 2007 - major re-design for run-time selection of adapt_type
! reichle, 27 Jun 2007 - new adapt_update_7, changed adapt_min/max
! reichle, 28 Jun 2007 - changed useless update_adapt_7 into update_adapt_8 
! reichle,  2 Jul 2007 - added spatial avg of innov stats (update_adapt_9)
! reichle,  3 Jul 2007 - added empirical factor to increase alpha_P (update_adapt_10)
! reichle, 19 Jul 2007 - clean up, add "adapt_misc_param", keep only update_adapt_10
! reichle, 24 Aug 2007 - added update_adapt_12
! reichle, 16 Jun 2011 - updated for new "obs_type" fields - COULD HAVE ADDED BUGS!!!
! reichle, 21 Nov 2014 - renamed force_pert_type fields for consistency w/ met_force_type
!                          %tmp2m --> %tair  (but note lower-case!)
!                          %dpt2m --> %qair  (but note lower-case!)
!                          %wnd   --> %wind  (but note lower-case!)

module clsm_adapt_routines

  use LDAS_ensdrv_globals,           ONLY:     &
       logit,                                     &
       logunit

  use LDAS_DateTimeMod,                   ONLY:     &
       date_time_type
  
  use catch_types,                      ONLY:     &
       cat_progn_type
  
  use force_and_cat_progn_pert_types,   ONLY:     &
       force_pert_real_type

  use LDAS_pertTypes,                  ONLY:     &
       pert_param_type

  use enkf_types,                       ONLY:     &
       obs_type,                                  &
       obs_param_type

  use LDAS_tilecoordtype,                 ONLY:     &
       tile_coord_type,                           &
       grid_def_type
  
  use adapt_types,                      ONLY:     &
       adapt_misc_param_type
  
  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename

  use LDAS_TilecoordRoutines,              ONLY:     &
       grid2tile

  use LDAS_exceptionsMod,                  ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR
  
  implicit none
  
  private
  
  public :: get_adapt_param
  public :: apply_adapt_P
  public :: apply_adapt_R
  public :: io_adapt_5
  public :: update_adapt_10
  public :: update_adapt_12

contains
  
  ! ***********************************************************************
  
  subroutine read_adapt_inputs( work_path, exp_id, date_time,              &
       adapt_type, adapt_misc_param, adapt_progn_pert, adapt_force_pert )
    
    ! read adapt inputs from nml file
    
    implicit none

    character(200),       intent(in) :: work_path
    character(40),        intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time
   
    integer,                     intent(out) :: adapt_type
    
    type(adapt_misc_param_type), intent(out) :: adapt_misc_param
    
    type(cat_progn_type),        intent(out) :: adapt_progn_pert
    type(force_pert_real_type),  intent(out) :: adapt_force_pert
    
    ! local variables

    character(200) :: adapt_inputs_path
    character( 40) :: adapt_inputs_file, dir_name, file_tag, file_ext
    
    character(300) :: fname

    logical        :: file_exists

    ! -----------------------------------------
    
    namelist / adapt_inputs / &
         adapt_type, adapt_misc_param, adapt_progn_pert, adapt_force_pert
    
    ! ---------------------------------------------------------------------
    !
    ! Set default file name for driver inputs namelist file
    
    adapt_inputs_path = './'                                       ! set default 
    adapt_inputs_file = 'LDASsa_DEFAULT_inputs_adapt.nml'
    
    ! Read data from default adapt_inputs namelist file 
    
    fname = trim(adapt_inputs_path) // '/' // trim(adapt_inputs_file)
    
    open (10, file=fname, delim='apostrophe', action='read', status='old')

    if (logit) write (logunit,*)
    if (logit) write (logunit,'(400A)') 'reading *default* adapt inputs from ', trim(fname)
    if (logit) write (logunit,*)
    
    read (10, nml=adapt_inputs)
    
    close(10,status='keep')
    
    
    ! Read from special adapt inputs file (if present) 
    
    adapt_inputs_file = 'LDASsa_SPECIAL_inputs_adapt.nml'
    
    ! Read data from special adapt_inputs namelist file 
    
    fname = trim(adapt_inputs_path) // '/' // trim(adapt_inputs_file)
    
    inquire(file=fname, exist=file_exists)

    if (file_exists) then        

       open (10, file=fname, delim='apostrophe', action='read', status='old')
       
       if (logit) write (logunit,*)
       if (logit) write (logunit,'(400A)') 'reading *special* adapt inputs from ', trim(fname)
       if (logit) write (logunit,*)
       
       read (10, nml=adapt_inputs)
       
       close(10,status='keep')

    end if
       
    ! echo variables of adapt_inputs
    
    if (logit) write (logunit,*) 'adapt inputs are:'
    if (logit) write (logunit,*)
    if (logit) write (logunit, nml=adapt_inputs) 
    if (logit) write (logunit,*)

    ! -------------------------------------------------------------
    !
    ! save adapt inputs into *adapt_inputs.nml file
    
    dir_name = 'rc_out'
    file_tag = 'ldas_adapt_inputs'
    file_ext = '.nml'
    
    fname = get_io_filename( work_path, exp_id, file_tag, date_time=date_time, &
         dir_name=dir_name, file_ext=file_ext )
    
    if (logit) write (logunit,*) 'writing adapt inputs to ', trim(fname)
    
    open (10, file=fname, status='unknown', action='write', delim='apostrophe')
    
    write(10, nml=adapt_inputs)
    
    close(10, status='keep')

  end subroutine read_adapt_inputs
  
  ! -----------------------------------------------------------------
  
  subroutine get_adapt_progn_pert_param( N_progn_pert, &
       progn_pert_param, adapt_progn_pert, progn_pert_adapt_param )
    
    ! transform adapt_progn_pert from nml into progn_pert_adapt_param
    
    implicit none
    
    integer, intent(in) :: N_progn_pert
    
    type(pert_param_type), dimension(N_progn_pert), intent(in) :: &
         progn_pert_param
    
    type(cat_progn_type), intent(in) :: adapt_progn_pert
    
    integer, dimension(N_progn_pert), intent(out) :: progn_pert_adapt_param
    
    ! local variables
    
    integer :: m

    character(len=*), parameter :: Iam = 'get_adapt_progn_pert_param'
    character(len=400) :: err_msg

    ! ------------------------------------------------------

    do m=1,N_progn_pert
       
       select case (trim(progn_pert_param(m)%descr))
          
       case ('tc1'   ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%tc1    
       case ('tc2'   ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%tc2    
       case ('tc4'   ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%tc4    
       case ('qa1'   ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%qa1    
       case ('qa2'   ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%qa2    
       case ('qa4'   ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%qa4    
       case ('capac' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%capac  
       case ('catdef') 
          progn_pert_adapt_param(m) = adapt_progn_pert%catdef 
       case ('rzexc' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%rzexc  
       case ('srfexc') 
          progn_pert_adapt_param(m) = adapt_progn_pert%srfexc 
       case ('ght1'  ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%ght(1) 
       case ('ght2'  ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%ght(2) 
       case ('ght3'  ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%ght(3) 
       case ('ght4'  ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%ght(4) 
       case ('ght5'  ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%ght(5) 
       case ('ght6'  ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%ght(6) 
       case ('wesn1' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%wesn(1)
       case ('wesn2' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%wesn(2)
       case ('wesn3' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%wesn(3)
       case ('htsn1' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%htsn(1)
       case ('htsn2' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%htsn(2)
       case ('htsn3' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%htsn(3)
       case ('sndz1' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%sndz(1)
       case ('sndz2' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%sndz(2)
       case ('sndz3' ) 
          progn_pert_adapt_param(m) = adapt_progn_pert%sndz(3) 
          
       case default
          
          err_msg = 'unknown progn_pert_param%descr = ' &
               // trim(progn_pert_param(m)%descr)
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end select
       
    end do
    
  end subroutine get_adapt_progn_pert_param
  
  ! -----------------------------------------------------------------------
  
  subroutine get_adapt_force_pert_param( N_force_pert, &
       force_pert_param, adapt_force_pert, force_pert_adapt_param )
    
    ! transform adapt_force_pert from nml into force_pert_adapt_param
    
    implicit none
    
    integer, intent(in) :: N_force_pert
    
    type(pert_param_type), dimension(N_force_pert), intent(in) :: &
         force_pert_param
    
    type(force_pert_real_type), intent(in) :: adapt_force_pert
    
    integer, dimension(N_force_pert), intent(out) :: force_pert_adapt_param
    
    ! local variables
    
    integer :: m

    character(len=*), parameter :: Iam = 'get_adapt_force_pert_param'
    character(len=400) :: err_msg
    
    ! ------------------------------------------------------

    do m=1,N_force_pert
       
       select case (trim(force_pert_param(m)%descr))
          
       case ('pcp'    ) 
          force_pert_adapt_param(m) = adapt_force_pert%pcp    
       case ('sw'     ) 
          force_pert_adapt_param(m) = adapt_force_pert%sw    
       case ('lw'     ) 
          force_pert_adapt_param(m) = adapt_force_pert%lw    
       case ('tair'   ) 
          force_pert_adapt_param(m) = adapt_force_pert%tair
       case ('qair'   ) 
          force_pert_adapt_param(m) = adapt_force_pert%qair    
       case ('wind'   ) 
          force_pert_adapt_param(m) = adapt_force_pert%wind    
          
       case default
          
          err_msg = 'unknown force_pert_param%descr = ' &
               // trim(force_pert_param(m)%descr)
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end select
       
    end do
    
  end subroutine get_adapt_force_pert_param
  
  ! -----------------------------------------------------------------------
  
  subroutine get_adapt_obs_pert_param( N_obs_param, obs_param, obs_pert_adapt_param )
    
    ! map between first dimension of Pert_adapt_R and obs_param
    
    implicit none

    integer, intent(in) :: N_obs_param

    type(obs_param_type), dimension(N_obs_param), intent(in) :: obs_param
    
    integer, dimension(N_obs_param), intent(out) :: obs_pert_adapt_param
    
    ! local variables
    
    integer :: i, n
    
    ! -------------------------------------------------------
    
    i = 0
    
    do n=1,N_obs_param
       
       if (obs_param(n)%assim)  then
          
          i = i+1
          
          obs_pert_adapt_param(n) = i
          
       else
          
          obs_pert_adapt_param(n) = -9999
          
       end if
       
    end do
    
  end subroutine get_adapt_obs_pert_param
  
  ! -----------------------------------------------------------------------
  
  subroutine get_adapt_param( work_path, exp_id, date_time,                   &
       N_progn_pert, N_force_pert, N_obs_param,                               &
       progn_pert_param, force_pert_param, obs_param,                         &
       adapt_type, adapt_misc_param, N_adapt_P, N_adapt_R,                    &
       progn_pert_adapt_param, force_pert_adapt_param, obs_pert_adapt_param )
    
    implicit none
    
    character(200),       intent(in) :: work_path
    character(40),        intent(in) :: exp_id

    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_progn_pert, N_force_pert, N_obs_param    

    type(pert_param_type), dimension(N_progn_pert), intent(in) :: progn_pert_param
    type(pert_param_type), dimension(N_force_pert), intent(in) :: force_pert_param
    
    type(obs_param_type), dimension(N_obs_param), intent(in) :: obs_param
    
    integer, intent(out) :: adapt_type, N_adapt_P, N_adapt_R

    type(adapt_misc_param_type), intent(out) :: adapt_misc_param

    integer, dimension(N_progn_pert), intent(out) :: progn_pert_adapt_param
    integer, dimension(N_force_pert), intent(out) :: force_pert_adapt_param
    
    integer, dimension(N_obs_param), intent(out) :: obs_pert_adapt_param
    
    ! local variables
    
    integer :: n
    
    type(cat_progn_type)       :: adapt_progn_pert
    type(force_pert_real_type) :: adapt_force_pert
    
    ! -------------------------------------------------------
    
    ! get adapt params from nml
    
    call read_adapt_inputs( work_path, exp_id, date_time,                    &
         adapt_type, adapt_misc_param, adapt_progn_pert, adapt_force_pert )
    
    ! initialize N_adapt_P and N_adapt_R

    N_adapt_P = 0
    N_adapt_R = 0
    
    if (adapt_type>0) then
       
       ! get N_adapt_P and N_adapt_R from obs_param
       
       do n=1,N_obs_param
          
          N_adapt_P = max( N_adapt_P, obs_param(n)%adapt )
          
          if (obs_param(n)%assim)  N_adapt_R = N_adapt_R + 1
          
       end do
       
       ! IMPORTANT: overwrite N_adapt_R for some cases of adapt_type
       
       select case (adapt_type)
          
       case (3,12)
          
          N_adapt_R = 0
          
       end select
       
       ! set up progn_pert_adapt_param
       
       call get_adapt_progn_pert_param( N_progn_pert, &
            progn_pert_param, adapt_progn_pert, progn_pert_adapt_param )
       
       ! set up force_pert_adapt_param
       
       call get_adapt_force_pert_param( N_force_pert, &
            force_pert_param, adapt_force_pert, force_pert_adapt_param )
       
       ! set up obs_pert_adapt_param
       
       call get_adapt_obs_pert_param( N_obs_param, obs_param, obs_pert_adapt_param )

    end if
    
  end subroutine get_adapt_param
  
  ! -----------------------------------------------------------------------
   
  subroutine io_adapt_5( adapt_tag, action, work_path, exp_id, date_time, &
       N_adapt, N_catd, N_obs_param, obs_param,                           &
       Pert_adapt_MA_XmXxOmB, Pert_adapt_MA_X, Pert_adapt )
    
    ! read/initialize or write Pert_adapt
    !
    ! adapt_tag = 'P' (adapting state err cov...) 
    ! adapt_tag = 'R' (adapting obs err cov...)
    !
    ! action_tag = 'r' or 'R' (read)
    ! action_tag = 'w' or 'W' (write)
    !
    ! reichle, 15 Dec 2006
    ! reichle,  1 Feb 2007 - modified for separate tuning of P and R
    
    implicit none
    
    character,      intent(in) :: adapt_tag  
    
    character,      intent(in) :: action     ! read ('r') or write ('w')
    
    character(200), intent(in) :: work_path
    character(40),  intent(in) :: exp_id
    
    type(date_time_type), intent(in) :: date_time
    
    integer, intent(in) :: N_adapt, N_catd, N_obs_param
    
    type(obs_param_type), dimension(N_obs_param), intent(in) :: obs_param
    
    real, dimension(N_adapt,N_catd), intent(inout) :: Pert_adapt_MA_XmXxOmB
    real, dimension(N_adapt,N_catd), intent(inout) :: Pert_adapt_MA_X
    real, dimension(N_adapt,N_catd), intent(inout) :: Pert_adapt
    
    ! local variables
    
    character(40), parameter :: file_tag = 'adapt'
    
    character( 40) :: tmp_file_tag, dir_name='rs', file_ext='.bin'
    character(300) :: fname
    
    integer :: n, k, istat
    
    logical :: tmp_write_action     

    character(len=*), parameter :: Iam = 'io_adapt_5'
    character(len=400) :: err_msg
    
    ! --------------------------------------------------------------------
    
    tmp_write_action = .false.
    
    if (adapt_tag=='P') then
       tmp_file_tag = 'adaptP_ldas_rst'
    elseif (adapt_tag=='R') then
       tmp_file_tag = 'adaptR_ldas_rst'
    else
       
       err_msg = 'unknown adapt_tag=' // adapt_tag
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    fname = get_io_filename( work_path, exp_id, tmp_file_tag, date_time=date_time, &
         dir_name=dir_name, ens_id=-1, file_ext=file_ext)
    
    if ( action=='r' .or. action=='R') then
       
       open(10, file=fname, form='unformatted', status='old', &
            action='read', iostat=istat)
       
       if (istat==0) then
          
          if (logit) write (logunit,*) 'Reading Pert_adapt file ', trim(fname)
          
          do k=1,N_adapt
             
             read (10) (Pert_adapt(k,n), n=1,N_catd)
             
          end do
          
          do k=1,N_adapt
             
             read (10) (Pert_adapt_MA_XmXxOmB(k,n), n=1,N_catd)
             
          end do
          
          do k=1,N_adapt
             
             read (10) (Pert_adapt_MA_X(k,n), n=1,N_catd)
             
          end do
          
          close(10, status='keep')

       else
          
          ! initialize and set "write" flag
          
          if (logit) write (logunit,*) 'Initializing Pert_adapt=1.'
          
          Pert_adapt = 1. 
          
          do n=1,N_obs_param
             
             if (obs_param(n)%assim) then
                
                k=obs_param(n)%adapt
                
                Pert_adapt_MA_XmXxOmB(k,:) = obs_param(n)%errstd**2
                Pert_adapt_MA_X(      k,:) = obs_param(n)%errstd**2
                
             end if
             
          end do
          
          tmp_write_action = .true.       
          
       end if
       
    end if
    
    if ( action=='w' .or. action=='W' .or. tmp_write_action ) then
       
       open(10,file=fname,form='unformatted',status='unknown',action='write')
       
       if (logit) write (logunit,*) 'Writing Pert_adapt file ', trim(fname)
       
       do k=1,N_adapt
          
          write (10) (Pert_adapt(k,n), n=1,N_catd)
          
       end do
       
       do k=1,N_adapt
          
          write (10) (Pert_adapt_MA_XmXxOmB(k,n), n=1,N_catd)
          
       end do
       
       do k=1,N_adapt
          
          write (10) (Pert_adapt_MA_X(k,n), n=1,N_catd)
          
       end do
       
       close (10,status='keep')
       
    end if
    
    if ( action/='w' .and. action/='W' .and.  &
         action/='r' .and. action/='R'           )    then
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown action')

    end if
    
  end subroutine io_adapt_5
  
  ! -----------------------------------------------------------------------

  subroutine update_adapt_10( N_obs_param, obs_param, N_obs, N_catd, Observations,  &
       obs_pert_adapt_param, N_adapt_P, N_adapt_R, adapt_misc_param,   &
       Pert_adapt_MA_AmBxOmB, Pert_adapt_MA_HPHt, Pert_adapt_P,        &
       Pert_adapt_MA_OmAxOmB, Pert_adapt_MA_R,    Pert_adapt_R      )
    
    ! reichle, 8 Feb 2007 - use tuning of Desroziers et al, QJR Met Soc, 2005
    !                       based on temporal avg of AmBxOmB divided by 
    !                       temporal avg of HPHt etc
    ! reichle, 27 Jun 2007 - correct design flaw in update_adapt_5 
    !                        (temporarily named update_adapt_7)
    ! reichle, 28 Jun 2007 - renamed update_adapt_8, introduced damping parameter delta
    ! reichle,  3 Jul 2007 - same as A0008 but with additional  
    !                         reduction of tmp_real by empirical factor
    ! reichle,  6 Jul 2007 - corrected lower bound of tmp_real to 1./(1.+delta)

    implicit none

    integer, intent(in) :: N_obs_param, N_obs, N_catd, N_adapt_P, N_adapt_R
    
    type(obs_param_type), dimension(N_obs_param), intent(in) :: &
         obs_param
    
    type(adapt_misc_param_type), intent(in) :: adapt_misc_param

    type(obs_type), dimension(N_obs), intent(in) :: Observations
    
    integer, dimension(N_obs_param), intent(in) :: obs_pert_adapt_param
    
    real, dimension(N_adapt_P,N_catd), intent(inout) :: Pert_adapt_MA_AmBxOmB
    real, dimension(N_adapt_P,N_catd), intent(inout) :: Pert_adapt_MA_HPHt
    real, dimension(N_adapt_P,N_catd), intent(inout) :: Pert_adapt_P
    real, dimension(N_adapt_R,N_catd), intent(inout) :: Pert_adapt_MA_OmAxOmB
    real, dimension(N_adapt_R,N_catd), intent(inout) :: Pert_adapt_MA_R
    real, dimension(N_adapt_R,N_catd), intent(inout) :: Pert_adapt_R
    
    ! local variables
    
    integer :: n, k, m
    
    real :: HPHt, AminusB, OminusF, OminusA, tmp_real

    !!!!real, parameter :: beta  = 1.06
    !!!!real, parameter :: delta = 0.005
    
    !!!!real, parameter :: w_P           = adapt_innov_weight_P
    !!!!real, parameter :: one_minus_w_P = 1 - w_P
    
    !!!!real, parameter :: w_R           = adapt_innov_weight_R
    !!!!real, parameter :: one_minus_w_R = 1 - w_R

    real :: w_P, one_minus_w_P, w_R, one_minus_w_R 
    
    ! -----------------------------------------------------------------
    
    w_P           = adapt_misc_param%gamma_P
    one_minus_w_P = 1. - w_P
    
    w_R           = adapt_misc_param%gamma_R
    one_minus_w_R = 1. - w_R
    
    
    do n=1,N_obs
       
       if (Observations(n)%assim) then
          
          m = Observations(n)%tilenum
          
          ! ---------------------
          !
          ! update Pert_adapt_P:
          
          k = obs_param( Observations(n)%species )%adapt
          
          ! compute analysis minus background (in obs space)
          
          AminusB = Observations(n)%ana - Observations(n)%fcst

          OminusF = Observations(n)%obs - Observations(n)%fcst

          OminusA = Observations(n)%obs - Observations(n)%ana
          
          ! compute expected HPHt
          
          HPHt = Observations(n)%fcstvar
          
          ! moving average
          
          Pert_adapt_MA_AmBxOmB(k,m) =                    &
               one_minus_w_P*Pert_adapt_MA_AmBxOmB(k,m) + &
               w_P*(AminusB*OminusF)
          
          Pert_adapt_MA_HPHt(k,m)    =                    &
               one_minus_w_P*Pert_adapt_MA_HPHt(k,m)    + &
               w_P*HPHt
          
          ! update Pert_adapt_P (omit whenever moving avg of HPHt<=0.)
          
          if (Pert_adapt_MA_HPHt(k,m)>0.) then
             
             tmp_real = Pert_adapt_MA_AmBxOmB(k,m)/Pert_adapt_MA_HPHt(k,m)
             
             ! --------------------------------------------------------------
             !
             ! multiply with empirical factor 
             ! A0008 tends to produce alphas that are too small in all cases
             
             tmp_real = tmp_real*adapt_misc_param%beta_P   !!!!!!!!!!!!
             
             ! --------------------------------------------------------------
             
             tmp_real = max(tmp_real,1./(1.+adapt_misc_param%delta_P))
             tmp_real = min(tmp_real,1.+adapt_misc_param%delta_P)
             
             Pert_adapt_P(k,m) = Pert_adapt_P(k,m)*tmp_real
             
          end if
          
          ! ensure adapt_min_P <= Pert_adapt <= adapt_max_P
          
          Pert_adapt_P(k,m) = &
               max( min( Pert_adapt_P(k,m), adapt_misc_param%max_alpha_P), &
               adapt_misc_param%min_alpha_P )
          
          ! ---------------------
          !
          ! update Pert_adapt_R:
          
          k = obs_pert_adapt_param( Observations(n)%species )
          
          ! moving average
          
          Pert_adapt_MA_OmAxOmB(k,m) =                              &
               one_minus_w_R*Pert_adapt_MA_OmAxOmB(k,m) +           &
               w_R*(OminusA*OminusF)
          
          Pert_adapt_MA_R(k,m) =                                    &
               one_minus_w_R*Pert_adapt_MA_R(k,m)       +           &
               w_R*(Observations(n)%obsvar)
          
          ! update Pert_adapt_R (omit whenever moving avg of R<=0.)
          
          if (Pert_adapt_MA_R(k,m)>0.) then
             
             tmp_real = Pert_adapt_MA_OmAxOmB(k,m)/Pert_adapt_MA_R(k,m)
          
             tmp_real = max(tmp_real,1./(1.+adapt_misc_param%delta_R))
             tmp_real = min(tmp_real,1.+adapt_misc_param%delta_R)
             
             Pert_adapt_R(k,m) = Pert_adapt_R(k,m)*tmp_real
   
          end if
          
          ! ensure adapt_min_R <= Pert_adapt <= adapt_max_R
          
          Pert_adapt_R(k,m) = &
               max( min( Pert_adapt_R(k,m), adapt_misc_param%max_alpha_R), &
               adapt_misc_param%min_alpha_R )
          
       end if
       
    end do
        
  end subroutine update_adapt_10

  ! ---------------------------------------------------------------------------

  subroutine update_adapt_12( N_obs_param, obs_param, N_obs, N_catd, Observations,  &
       N_adapt_P, adapt_misc_param,                                    &
       Pert_adapt_MA_AmBxOmB, Pert_adapt_MA_HPHt, Pert_adapt_P )
    
    ! reichle,24 Aug 2007 - same as update_adapt_10 but for Pert_adapt_P *only*

    implicit none

    integer, intent(in) :: N_obs_param, N_obs, N_catd, N_adapt_P
            
    type(obs_param_type), dimension(N_obs_param), intent(in) :: &
         obs_param
    
    type(adapt_misc_param_type), intent(in) :: adapt_misc_param

    type(obs_type), dimension(N_obs), intent(in) :: Observations
            
    real, dimension(N_adapt_P,N_catd), intent(inout) :: Pert_adapt_MA_AmBxOmB
    real, dimension(N_adapt_P,N_catd), intent(inout) :: Pert_adapt_MA_HPHt
    real, dimension(N_adapt_P,N_catd), intent(inout) :: Pert_adapt_P
        
    ! local variables
    
    integer :: n, k, m
    
    real :: HPHt, AminusB, OminusF, OminusA, tmp_real

    real :: w_P, one_minus_w_P
    
    ! -----------------------------------------------------------------
    
    w_P           = adapt_misc_param%gamma_P
    one_minus_w_P = 1. - w_P
    
    do n=1,N_obs
       
       if (Observations(n)%assim) then
          
          m = Observations(n)%tilenum
          
          ! ---------------------
          !
          ! update Pert_adapt_P:
          
          k = obs_param( Observations(n)%species )%adapt
          
          ! compute analysis minus background (in obs space)
          
          AminusB = Observations(n)%ana - Observations(n)%fcst

          OminusF = Observations(n)%obs - Observations(n)%fcst

          OminusA = Observations(n)%obs - Observations(n)%ana          

          ! compute expected HPHt
          
          HPHt = Observations(n)%fcstvar
          
          ! moving average
          
          Pert_adapt_MA_AmBxOmB(k,m) =                    &
               one_minus_w_P*Pert_adapt_MA_AmBxOmB(k,m) + &
               w_P*(AminusB*OminusF)
          
          Pert_adapt_MA_HPHt(k,m)    =                    &
               one_minus_w_P*Pert_adapt_MA_HPHt(k,m)    + &
               w_P*HPHt
          
          ! update Pert_adapt_P (omit whenever moving avg of HPHt<=0.)
          
          if (Pert_adapt_MA_HPHt(k,m)>0.) then
             
             tmp_real = Pert_adapt_MA_AmBxOmB(k,m)/Pert_adapt_MA_HPHt(k,m)
             
             ! --------------------------------------------------------------
             !
             ! multiply with empirical factor 
             ! A0008 tends to produce alphas that are too small in all cases
             
             tmp_real = tmp_real*adapt_misc_param%beta_P   !!!!!!!!!!!!
             
             ! --------------------------------------------------------------
             
             tmp_real = max(tmp_real,1./(1.+adapt_misc_param%delta_P))
             tmp_real = min(tmp_real,1.+adapt_misc_param%delta_P)
             
             Pert_adapt_P(k,m) = Pert_adapt_P(k,m)*tmp_real
             
          end if
          
          ! ensure adapt_min_P <= Pert_adapt <= adapt_max_P
          
          Pert_adapt_P(k,m) = &
               max( min( Pert_adapt_P(k,m), adapt_misc_param%max_alpha_P), &
               adapt_misc_param%min_alpha_P )
          
          
       end if
       
    end do
    
  end subroutine update_adapt_12

  ! ---------------------------------------------------------------------------
  
  subroutine apply_adapt_P( N_pert, pert_param, pert_adapt_param,                  &
       N_adapt_P, N_catd, Pert_adapt_P, pert_grid, tile_coord, N_ens, Pert_tile )
    
    ! re-scale "Force_pert" or "Progn_pert" based on Pert_adapt_P
    ! (pert_param() never changes, ie. "Force_pert" and "Progn_pert" 
    !  are computed with fixed nml inputs for perturbations parameters)
    !
    ! reichle, 15 Dec 2006
    !
    ! reichle, 19 Dec 2006 - empirical rule:
    !                        apply adapt using *square* of Pert_adapt, because
    !                        doubling of std_pert does *not* lead to
    !                        doubling of sqrt(HPHt)
    ! reichle, 28 Jun 2007 - deleted empirical rule (RedArk_adapt A0007 did not work)
    
    implicit none
    
    integer, intent(in) :: N_pert, N_adapt_P, N_catd, N_ens

    type(pert_param_type), dimension(N_pert), intent(in) :: pert_param
    
    integer, dimension(N_pert), intent(in) :: pert_adapt_param
    
    real, dimension(N_adapt_P,N_catd) :: Pert_adapt_P
    
    type(grid_def_type) :: pert_grid
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord
    
    real, dimension(N_pert,N_catd,N_ens) :: Pert_tile
    
    ! local variables
    
    integer :: n, m, k, j
    
    real, dimension(N_catd) :: mu, sg
    
    real :: s_square, tmp_x, tmp_mu, tmp_sg, tmp_real

    character(len=*), parameter :: Iam = 'apply_adapt_P'
    character(len=400) :: err_msg
    
    ! -----------------------------------------------
    
    do n=1,N_pert
       
       j = pert_adapt_param(n)
       
       if (j>0) then
          
          ! get mean and std in tile space 
          
          call grid2tile( pert_grid, N_catd, tile_coord%i_indg,tile_coord%j_indg, & !tile_coord,       &
               pert_param(n)%mean, mu )

          call grid2tile( pert_grid, N_catd, tile_coord%i_indg,tile_coord%j_indg, & !tile_coord,       &
               pert_param(n)%std,  sg )
          
          select case (pert_param(n)%typ)
             
          case (0)         ! additive
             
             do k=1,N_catd
                do m=1,N_ens
                   
                   !Pert_tile(n,k,m) = &
                   !     Pert_adapt_P(j,k)*(Pert_tile(n,k,m)-mu(k)) + mu(k)
                   
                   Pert_tile(n,k,m) = &
                        sqrt(Pert_adapt_P(j,k))*(Pert_tile(n,k,m)-mu(k)) + mu(k)
                   
                end do
             end do
             
          case (1)         ! multiplicative and lognormal (mean=1)
             
             do k=1,N_catd
                do m=1,N_ens
                   
                   ! compute new lognormal parameters mean and std (tmp_mu, 
                   ! tmp_sg) instead of scaling Pert_tile
                   ! tmp_x    = original std_normal perturbation that
                   !             was used to compute Pert_tile
                   ! s_square = original variance of Pert_tile
                   
                   s_square = exp(-2.*mu(k)) - 1.
                   
                   tmp_x    = (log(Pert_tile(n,k,m)) - mu(k))/sg(k)
                   
                   tmp_real = log(1+Pert_adapt_P(j,k)*s_square)
                   !tmp_real = log(1+(Pert_adapt_P(j,k)**2)*s_square)
                   
                   tmp_mu   = -.5*tmp_real
                   
                   tmp_sg   = sqrt(tmp_real);
                   
                   Pert_tile(n,k,m) = exp( tmp_mu + tmp_sg*tmp_x );
                   
                end do
             end do
             
          case default
             
             err_msg = 'unknown typ of perturbation'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

          end select
          
       end if
       
    end do
    
  end subroutine apply_adapt_P

  ! ---------------------------------------------------------------------
  
  subroutine apply_adapt_R( N_obs, N_obs_param, obs_pert_adapt_param,  &
       N_adapt_R, N_catd, Pert_adapt_R, Observations )
    
    ! reichle,  1 Feb 2007
    
    implicit none
    
    integer, intent(in) :: N_obs, N_obs_param, N_adapt_R, N_catd
    
    integer, dimension(N_obs_param), intent(in) :: obs_pert_adapt_param
    
    real, dimension(N_adapt_R,N_catd) :: Pert_adapt_R

    type(obs_type), dimension(N_obs), intent(inout) :: Observations

    ! local variables
    
    integer :: n, k, j
    
    ! -----------------------------------------------
    
    do n=1,N_obs
       
       j = obs_pert_adapt_param(Observations(n)%species)
       k = Observations(n)%tilenum
       
       ! THE FOLLOWING LINE WAS CHANGED WHEN OBS_TYPE FIELDS WERE REVISED
       !
       ! NOT SURE WHETHER THE REVISED LINE IS OK
       !
       ! REICHLE - 16 JUNE 2011

       ! Observations(n)%errstd = sqrt(Pert_adapt_R(j,k))*Observations(n)%errstd 
       
       Observations(n)%obsvar = Pert_adapt_R(j,k)*Observations(n)%obsvar
       
    end do
    
  end subroutine apply_adapt_R
  
  ! -----------------------------------------------------------------------
  
end module clsm_adapt_routines

! ====================== EOF ==============================================
