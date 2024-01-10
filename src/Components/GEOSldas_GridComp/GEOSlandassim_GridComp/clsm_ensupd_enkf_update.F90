
!***********************************************************************
!
!  land EnKF update for off-line CLSM ensemble driver
!
!  Rolf Reichle, 18 Jul 2005
!
!***********************************************************************

module clsm_ensupd_enkf_update

  use MAPL_SortMod,                     ONLY:     &
       MAPL_Sort

  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_TICE

  USE CATCH_CONSTANTS,                  ONLY :    &
       N_gt           => CATCH_N_GT

  use catchment_model,                  ONLY:     &
       catch_calc_tsurf

  use lsm_routines,                     ONLY:     &
       catch_calc_soil_moist,                     &
       catch_calc_tp
  
  use LDAS_ensdrv_globals,              ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_generic,                            &
       nodata_tolfrac_generic

  use LDAS_DateTimeMod,                 ONLY:     &
       date_time_type,                            &
       date_time2string

  use enkf_types,                       ONLY:     &
       obs_param_type,                            &
       obs_type

  use LDAS_DriverTypes,                 ONLY:     &
       met_force_type
       
  use catch_types,                      ONLY:     &
       cat_param_type,                            &
       cat_progn_type,                            &
       catprogn2wesn,                             &
       catprogn2htsn,                             &
       catprogn2ghtcnt,                           &
       assignment (=),                            &
       operator (+),                              &
       operator (/)

  use mwRTM_types,                      ONLY:     &
       mwRTM_param_type
  
  use LDAS_PertTypes,                   ONLY:     &
       pert_param_type
  
  use LDAS_TilecoordType,               ONLY:     &
       tile_coord_type,                           &
       grid_def_type
  
  use catch_bias_types,                 ONLY:     &
       obs_bias_type
  
  use LDAS_TilecoordRoutines,           ONLY:     &
       get_number_of_tiles_in_cell_ij,            &
       get_tile_num_in_cell_ij,                   &
       grid2tile

  use nr_ran2_gasdev,                   ONLY:     &
       NRANDSEED

  use ease_conv,                        ONLY:     &
       ease_convert

  use my_matrix_functions,              ONLY:     &
       row_std

  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename

  use clsm_ensdrv_drv_routines,         ONLY:     &
       check_cat_progn,                           &
       l2f_real

  use clsm_ensupd_upd_routines,         ONLY:     &
       get_obs_pred,                              &
       get_halo_obs,                              &
       get_obs_pert,                              &
       cat_enkf_increments,                       &
       get_ind_obs_assim,                         &
       get_ind_obs_lat_lon_box,                   &
       halo_type,                                 &
       FOV_threshold,                             &
       get_halo_around_tile,                      &
       TileNnzObs

  use clsm_ensupd_read_obs,             ONLY:     &
       collect_obs

  use clsm_bias_routines,               ONLY:     &
       obs_bias_upd_tcount,                       &
       obs_bias_corr_obs,                         &
       obs_bias_upd_bias_and_Obs

  use clsm_adapt_routines,              ONLY:     &
       apply_adapt_R

  use LDAS_ensdrv_mpi,                  ONLY:     &
       MPI_cat_param_type,                        &
       MPI_cat_progn_type,                        &
       root_proc,                                 &
       numprocs,                                  &
       myid,                                      &
       mpierr,                                    &
       mpicomm,                                   &
       MPI_obs_type,                              &
       mpistatus

  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  implicit none

  include 'mpif.h'

  private

  public :: get_enkf_increments
  public :: apply_enkf_increments
  public :: output_ObsFcstAna_wrapper
  public :: write_smapL4SMaup 

contains

  subroutine get_enkf_increments(                                        &
       date_time,                                                        &
       N_ens, N_catl, N_catf, N_obsl_max,                                &
       work_path, exp_id,                                                &
       met_force, lai, cat_param, mwRTM_param,                           &
       tile_coord_l, tile_coord_f,                                       &
       tile_grid_g, pert_grid_f, pert_grid_g,                            &
       N_catl_vec, low_ind, l2f, f2l,                                    &
       update_type,                                                      &
       dtstep_assim,                                                     &
       xcompact, ycompact, fcsterr_inflation_fac,                        &
       N_obs_param, obs_param, N_obsbias_max,                            &
       out_obslog, out_smapL4SMaup,                                      &
       cat_progn,                                                        &
       Pert_rseed, obs_bias,                                             &
       cat_progn_incr, fresh_incr,                                       &
       N_obsf, N_obsl, Observations_l,                                   &
       N_adapt_R, obs_pert_adapt_param, Pert_adapt_R,                    &
       Obs_pert )

    ! -------------------------------------------------------------

    ! return increments instead of updated cat_progn
    ! reichle, 18 Oct 2005

    implicit none

    ! ----------------------------------------------------------------
    !
    ! inputs (via argument list, command line, or namelist file, TBD)

    type(date_time_type), intent(in) :: date_time    ! current date, time

    integer, intent(in) :: N_ens          ! number of ensemble members
    integer, intent(in) :: N_catl         ! # tiles in *l*ocal procs subdomain
    integer, intent(in) :: N_catf         ! # tiles in *f*ull domain
    integer, intent(in) :: N_obsl_max     ! max number of observations allowed

    character(*), intent(in) :: work_path
    character(*), intent(in) :: exp_id


    ! Meteorological forcings, Catchment model and microwave RTM parameters

    type(met_force_type),   dimension(N_catl), intent(in) :: met_force
    real,                   dimension(N_catl), intent(in) :: lai
    type(cat_param_type),   dimension(N_catl), intent(in) :: cat_param

    type(mwRTM_param_type), dimension(N_catl), intent(in) :: mwRTM_param

    ! grid and tile coordinate variables

    type(tile_coord_type), dimension(:), pointer :: tile_coord_l, tile_coord_f  ! input

    type(grid_def_type), intent(in) :: tile_grid_g, pert_grid_f, pert_grid_g

    integer, intent(in), dimension(numprocs) :: N_catl_vec, low_ind

    integer, intent(in), dimension(N_catl)   :: l2f

    integer, intent(in), dimension(N_catf)   :: f2l

    integer, intent(in) :: update_type, dtstep_assim

    real,    intent(in) :: xcompact, ycompact, fcsterr_inflation_fac

    integer, intent(in) :: N_obs_param

    integer, intent(in) :: N_obsbias_max  ! max number of obs bias parameters (bias_Npar)

    type(obs_param_type), dimension(N_obs_param), intent(in) :: obs_param

    logical, intent(in) :: out_obslog
    logical, intent(in) :: out_smapL4SMaup

    type(cat_progn_type), dimension(N_catl,N_ens), intent(in) :: cat_progn

    ! intent(inout) variables:

    integer, dimension(NRANDSEED,N_ens), intent(inout) :: Pert_rseed

    type(obs_bias_type), dimension(N_catl,N_obs_param,N_obsbias_max), intent(inout) :: &
         obs_bias

    ! intent(out) variables:

    type(cat_progn_type),dimension(N_catl,N_ens),intent(out) :: cat_progn_incr

    logical, intent(inout) :: fresh_incr

    type(obs_type), dimension(:), pointer :: Observations_l  ! output

    integer, intent(out) :: N_obsf, N_obsl

    ! must always have some variables for adaptive filtering
    ! (GNUMakefile setup does not permit C-preprocessor directives)

    integer, intent(in) :: N_adapt_R

    integer, dimension(N_obs_param), intent(in) :: obs_pert_adapt_param

    real, dimension(N_adapt_R,N_catl), intent(in) :: Pert_adapt_R

    real, dimension(N_obsl_max,N_ens), intent(out), optional :: Obs_pert

    ! ----------------------------------------------

    ! local variables

    logical :: found_obs_f, assimflag

    type(obs_type), dimension(:), pointer :: Observations_lH => null()  ! obs w/in halo

    ! matrix of measurement predictions

    real, dimension(:,:), pointer :: Obs_pred_l  => null()
    real, dimension(:,:), pointer :: Obs_pred_lH => null()

    ! random realization of measurement error

    real, allocatable, dimension(:,:) :: Obs_pert_tmp

    ! variables for echoing wall clock time

    character(8)  :: date_string
    character(10) :: time_string

    ! additional grid/tile information that is needed for mapping observations
    ! to tiles

    integer, dimension(:,:),   pointer     :: N_tile_in_cell_ij_f
    integer, dimension(:,:,:), pointer     :: tile_num_in_cell_ij_f

    integer :: i, n, n_e

    ! obs bias variables

    logical, dimension(:), allocatable     :: obsbias_ok

    ! pchakrab: variables for analysis load balance (AnaLoadBal)

    type varLenIntArr                                    ! used to store indices on each processor
       integer, dimension(:), allocatable :: ind
    end type varLenIntArr
    integer                            :: N_select_species  ! input to get_ind_obs_lat_lon_box()
    integer, dimension(:), allocatable :: select_species    ! input to get_ind_obs_lat_lon_box() and TileNnzObs()
    type(halo_type)                    :: halo
    integer                            :: N_selected_obs
    integer, dimension(numprocs)       :: tmp_low_ind    ! tmp_low_ind-1 is the displs vector for Gatherv/Scatterv

    ! tiles related
    integer                            :: nTiles_l,   nTilesl_vec(numprocs),  nTiles_f
    integer                            :: nTiles_ana, nTilesAna_vec(numprocs)
    integer, dimension(:), allocatable :: indTiles_l, indTiles_f, indTiles_ana
    type(varLenIntArr)                 :: indTilesAna_vec(numprocs)
    type(tile_coord_type), dimension(:), pointer    :: tile_coord_ana  ! input to cat_enkf_increment() is a pointer
    type(cat_param_type), dimension(:), allocatable :: cat_param_f,         cat_param_ana
    type(cat_progn_type), allocatable               :: cat_progn_f(:),      cat_progn_ana(:,:)
    type(cat_progn_type), allocatable               :: tmp_cat_progn_ana(:)
    type(cat_progn_type), allocatable               :: cat_progn_incr_f(:), cat_progn_incr_ana(:,:)
    type(cat_progn_type), allocatable               :: recvBuf(:)

    ! Obs related
    integer                                         :: nObs_ana
    integer                                         :: nObsAna_vec(numprocs)
    integer                                         :: N_obsf_assim, N_obsl_assim
    integer                                         :: N_obsl_assim_vec(numprocs)
    integer, dimension(:), allocatable              :: indObs_ana
    integer, dimension(:), allocatable, target      :: ind_obsl_assim
    integer, dimension(:), pointer                  :: ptr2indx => null()
    type(varLenIntArr)                              :: indObsAna_vec(numprocs)
    integer, dimension(:), allocatable              :: tmp_ind_obs
    type(obs_type), dimension(:), allocatable       :: Obs_f_assim, Obs_ana  ! collect obs before distributing for ana
    real, allocatable                               :: Obs_pred_f_assim(:), Obs_pred_ana(:,:)

    ! odds and ends
    real               :: t_start, t_end, tmax, tmin            ! for timing routines
    integer            :: iTile, iproc, iEns, ctr               ! counters
    integer            :: quotient, remainder
    integer            :: dest, src, sendct, recvct, sendtag, recvtag

    character(12)      :: tmpstr12

    character(len=*), parameter :: Iam = 'get_enkf_increments'
    character(len=400) :: err_msg

    ! **********************************************************************
    !
    !                     END OF DECLARATIONS
    !
    ! ***********************************************************************

    ! nullify all pointers
    ! (good practice; necessary on halem when -omp is used)

    nullify(N_tile_in_cell_ij_f, tile_num_in_cell_ij_f)
    nullify(Observations_lH,Obs_pred_l,Obs_pred_lH)

    ! *************************************************************************

    ! initialize

    fresh_incr = .false.

    N_obsl     = 0
    N_obsf     = 0

    do n=1,N_catl
       do n_e=1,N_ens
          cat_progn_incr(n,n_e) = 0.
       end do
    end do
    
    ! check if update is needed at all

    if (update_type==0) then

       if (logit) write (logunit,*) 'no EnKF increments b/c update_type=0'

       return    ! nothing else to be done here

    end if

    ! echo analysis time and wall clock time

    call date_and_time(date_string, time_string)  ! f90 intrinsic function

    if (logit) write (logunit,*) 'get_enkf_increments(): enter at ', &
         date_string, ', ', time_string
    if (logit) write (logunit,*) 'get_enkf_increments(): enter at anal time ', &
         date_time2string(date_time)

    if (.true.) then  ! replace obsolete check for analysis time with "if true" to keep indents
       
       ! proceed with update

       ! -----------------------------------------------------------------
       !
       ! Get additional grid/tile information that is needed to map obs
       ! from lat/lon to tiles.  This needs to be done:
       ! - by root process (because of call to read_obs() in collect_obs())
       ! - by all processes if FOV>~0 ("tile_num_in_circle" needed in get_obs_pred())
       if ( (root_proc)                                        .or.              &
            (any(obs_param(1:N_obs_param)%FOV>FOV_threshold))        )  then

          allocate(N_tile_in_cell_ij_f(pert_grid_f%N_lon,pert_grid_f%N_lat))

          ! first call: count how many tiles are in each pert_grid_f cell
          call get_number_of_tiles_in_cell_ij( N_catf,                           &
               tile_coord_f%pert_i_indg, tile_coord_f%pert_j_indg,               &
               pert_grid_f, N_tile_in_cell_ij_f )
          ! second call: find out which tiles are in each pert_grid_f cell

          call get_tile_num_in_cell_ij( N_catf,                                  &
               tile_coord_f%pert_i_indg, tile_coord_f%pert_j_indg,               &
               pert_grid_f, maxval(N_tile_in_cell_ij_f), tile_num_in_cell_ij_f )
       else
          allocate(N_tile_in_cell_ij_f(0,0)) !for debugging
       end if

       ! *********************************************************************
       !
       ! collect observations
       !
       ! *********************************************************************

       ! write original (unscaled) SMAP Tb observations into SMAP L4_SM aup file
       !
       ! NOTE: this requires a call to collect_obs() within output_smapL4SMaup()

       if (out_smapL4SMaup)                                                        &
            call output_smapL4SMaup( date_time, work_path, exp_id, dtstep_assim,   &
            N_ens, N_catl, N_catf, N_obsl_max,                                     &
            tile_coord_f, tile_grid_g, pert_grid_f,                                &
            N_catl_vec, low_ind, l2f, N_tile_in_cell_ij_f, tile_num_in_cell_ij_f,  &
            N_obs_param, obs_param, Observations_l, cat_param, cat_progn  )

       ! check for observations, found_obs_f=.true. if obs availalbe

       call collect_obs(                                              &
            work_path, exp_id, date_time, dtstep_assim,               &
            N_catl,                                                   &
            N_catf, tile_coord_f, pert_grid_f,                        &
            N_tile_in_cell_ij_f, tile_num_in_cell_ij_f,               &
            N_catl_vec, low_ind, l2f,                                 &
            N_obs_param, obs_param, N_obsl_max, out_obslog,           &
            N_obsl, Observations_l, found_obs_f )


       ! --------------------------------------------------------------------

       if (found_obs_f) then

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !
          ! *** Obs bias algorithm ***
          !
          ! Draper, Aug 2013
          !
          ! Obs bias update equations:
          ! i)    b^+ = b^- + l (y - Hx^- - b^- ) for all obs
          ! ii)   x^+ = x^- + K| (y| - b^+ -H| x^- )
          !              y|,H|,K| for obs only where bias estimate is good
          !
          ! B1. Adjust observation to remove the prior bias  (y'=(y - b^-))
          !     Determine whether have good obs_bias estimate based on time
          !     since last bias update
          !     -> set obsbias_ok and Observations%assim flags on each obs
          !     (get_obs_pred call is here - has some post model-based QC)
          ! B2. Calculate the obs_bias incr (b+ - b-)
          !    a Update the obs_bias with obs_bias incr to get posterior bias
          !    b Update the obs with with obs_bias incr (y - b+ = y' - (b^+ - b^-)
          !
          ! Note: get_halo_obs() screens for obs%assim flag, so only obs
          !       that have assim==.T. will be used in state update
          !
          ! CAUTION: observed value (Observations%obs) is now bias-corrected!!!
          !          (and written out as such)

          allocate(obsbias_ok(N_obsl))

          if (N_obsl>0) obsbias_ok = .false.          ! initialize

	  if ( (N_obsl>0) .and. (N_obsbias_max>0) ) then

      ! B1. Adjust observations to remove the obs bias
      !     and set obsbias_ok flag

             call obs_bias_corr_obs(date_time, N_catl, N_catf,           &
                  N_obsl , N_obs_param, N_obsbias_max, f2l, obs_param,   &
                  obs_bias, Observations_l, obsbias_ok)

          end if

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++

          ! adapt obs error covariance R

          ! not yet implemented for MPI

          if (N_adapt_R>0)                                                       &
               call apply_adapt_R( N_obsl, N_obs_param, obs_pert_adapt_param,    &
               N_adapt_R, N_catl, Pert_adapt_R, Observations_l )

          ! ******************************************************************
          !
          ! compute innovations (O-F)
          !
          ! ******************************************************************

          call date_and_time(date_string, time_string)  ! f90 intrinsic function

          if (logit) write (logunit,'(400A)') 'computing innovations starting at ' // &
               date_string // ', ' // time_string
          
          ! compute model forecast of observations
          ! (ensemble mean "obs_pred" is also stored in Observations_l%fcst)

          ! incl. obs bias estimation

          call get_obs_pred(                                      &
               .true.,                                            & ! -> before EnKF update
               N_obs_param, N_ens,                                &
               N_catl, tile_coord_l,                              &
               N_catf, tile_coord_f, f2l,                         &
               N_catl_vec, low_ind, pert_grid_g,                  &
               obs_param,                                         &
               met_force, lai, cat_param, cat_progn, mwRTM_param, &
               N_obsl, Observations_l, Obs_pred_l, obsbias_ok,    &
               fcsterr_inflation_fac )

          if (allocated(obsbias_ok)) deallocate(obsbias_ok)


          ! IF NEEDED, INCLUDE WITHHOLDING SUBROUTINE HERE.
          ! SUCH A SUBROUTINE SHOULD CHANGE Observations(i)%assim TO FALSE
          ! IF THE OBSERVATION IS TO BE WITHHELD
          !
          ! call withhold_obs()


          ! count observations across all processors that are left after
          !  model-based QC (done within get_obs_pred())

#ifdef LDAS_MPI

          call MPI_AllReduce(                             &
               N_obsl, N_obsf, 1, MPI_integer, MPI_SUM,   &
               mpicomm, mpierr )

#else
          N_obsf = N_obsl
#endif

          ! check whether any "assim" flag is set in obs_param

	  ! CSD - if want to skip cat_enkf_incr and apply_incr blocks
	  ! in instances where obs bias correction has discared all obs
	  ! (first few cycles of avail. obs), would need to test here
	  ! for whether obs on any processors are being assimilated.

          assimflag = .false.

          do i=1,N_obs_param

             if (obs_param(i)%assim)  assimflag = .true.

          end do

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++
          !
          ! Obs bias
          !
          ! B2. Update obs_bias and Observations with the obs bias increment

	  if ( (N_obsl>0) .and. (N_obsbias_max>0) )                              &
               call obs_bias_upd_bias_and_Obs(date_time, N_catl, N_catf,           &
               N_obsl, N_ens,  N_obs_param, N_obsbias_max, f2l, obs_param,         &
               Obs_pred_l(1:N_obsl,1:N_ens), obs_bias, Observations_l(1:N_obsl) )

   ! ++++++++++++++++++++++++++++++++++++++++++++++++++++

       end if ! found_obs_f==.true.


       ! ******************************************************************
       !
       ! assimilate (EnKF state update)
       !
       ! ******************************************************************

       if ( (N_obsf>0) .and. assimflag ) then

#ifdef LDAS_MPI

          ! AnaLoadBal (Analysis Load Balance)
          !
          ! pchakrab: 03/25/2014
          !
          ! Issue: Load imbalance of cat_enkf_increment. For single-sensor
          ! assimilation (e.g., SMAP) the assimilated obs end up on only a few 
          ! processors and the majority of processors sit idle.
          !
          ! Solution: Collect and distribute tiles with a non-zero (nnz) number of 
          ! assimilated obs evenly among processors. This is not the ideal solution since
          ! the number of obs still varies considerably from  processor to
          ! processor, but is a major improvement. A 24-hr run (with 24
          ! ensembles and simulated SMAP data using the SMAP S/W Delivery 5 CVS tag) 
          ! on 128 Sandybridge processors went from 4h:40m (original code) 
          ! to 1h:10m (load balanced code). On 64 processors the same run took 1h:40m.
          !
          ! Outline of algorithm: The load balancing is done in 5 steps.
          !
          ! Prereq
          ! ------
          ! Step 1: Combine Observations_l (obs w/ %assim==.true.) into
          ! Obs_f_assim (available on all processors). Obs_f_assim is needed 
          ! to identify the (local) tiles with nnz obs.
          ! NOTE: N_obsf:       number of all obs
          !       N_obsf_assim: number of obs w/ assim==.true.
          ! 
          ! Decomposition (2 steps)
          ! -----------------------
          ! Step 2: Each processor identifies the local tiles with nnz obs,
          ! indTiles_l, collects them on root, indTiles_f, and distributes
          ! them evenly among all procs, indTiles_ana. The corresponding
          ! numbers are nTiles_l, nTiles_f, nTiles_ana. Root needs
          ! nTilesAna_vec, indTilesAna_vec (list of nTiles_ana,
          ! indTiles_ana on each proc) to distribute cat_param, cat_progn.
          !
          ! IMPORTANT: Regardless of update_type, obs from *all* species are
          !            considered (ie, N_select_species=0).  This could result in
          !            poor load balancing if different species are present, 
          !            which will need to be addressed in future. 
          !            - wjiang+reichle, 13 Oct 2020
          !
          ! Step 3: Each processor computes nObs_ana and indObs_ana,
          ! the number and indices of obs affecting tiles in indTiles_ana.
          ! Root needs nObsAna_vec, indObsAna_vec to distribute Obs_pred_l.
          !
          ! Distribute input
          ! ----------------
          ! Step 4: On each proc, create tile_coord_ana, Obs_ana etc. the
          ! load-balanced versions of tile_coord_l, Observations_l etc.
          ! (input to cat_enkf_increments).
          !
          !      call to cat_enkf_increments w/ load-balanced data
          !
          ! Collect output
          ! --------------
          ! Step 5: Output of cat_enkf_increment() is cat_progn_incr_ana
          ! which is then massaged into its 'local' counterpart

          call date_and_time(date_string, time_string)  ! f90 intrinsic function
          
          if (logit) write(logunit,'(400A)')                                       &
               'Dynamic load balancing and analysis (AnaLoadBal) starting at ' //  &
               date_string // ', ' // time_string

          !-AnaLoadBal-Prereq-starts-here-
          ! Step 1a: identify obs w/ obs%assim==.true.
          allocate(ind_obsl_assim(N_obsl), source=-99)
          call get_ind_obs_assim(N_obsl, Observations_l%assim, N_obsl_assim, ind_obsl_assim)
          ! its easier to write ptr2indx than ind_obsl_assim(1:N_obsl_assim)
          ptr2indx => ind_obsl_assim(1:N_obsl_assim)

          ! Step 1b: Observations_l(obs%assim=.true.) -> Obs_f_assim (on all processors)
          ! NOTE: For MPI_Allgatherv, we need N_obsl_assim_vec and tmp_low_ind on all procs
          call MPI_AllGather(N_obsl_assim,1,MPI_INTEGER,                                &
               N_obsl_assim_vec,1,MPI_INTEGER,mpicomm,mpierr )
          N_obsf_assim = sum(N_obsl_assim_vec)
          if (logit) then             
             write (tmpstr12,'(i12)') N_obsf         ! convert integer to string             
             write(logunit,'(400A)') 'Total number of obs:             ' // &
                  tmpstr12, ' [after model-based QC]'
             write (tmpstr12,'(i12)') N_obsf_assim   ! convert integer to string             
             write(logunit,'(400A)') 'Total number of assimilated obs: ' // &
                  tmpstr12, ' [after "assim_flag==true"]'
          end if
          allocate(Obs_f_assim(N_obsf_assim))
          tmp_low_ind(1) = 1
          do n=1,numprocs-1
             tmp_low_ind(n+1) = tmp_low_ind(n) + N_obsl_assim_vec(n)
          end do
          call MPI_Allgatherv(                                                          &
               Observations_l(ptr2indx), N_obsl_assim,                    MPI_obs_type, &
               Obs_f_assim,              N_obsl_assim_vec, tmp_low_ind-1, MPI_obs_type, &
               mpicomm, mpierr)
          !-AnaLoadBal-Prereq-ends-here-

          !-AnaLodaBal-Decomposition-starts-here
          ! Step 2a: compute nTiles_l, indTiles_l, nTilesl_vec (on root), nTiles_f (on root)
          ! NOTE: loop over tile_coord_l, if tile has nnz obs, store the 'full' index
          call cpu_time(t_start)
          allocate(indTiles_l(N_catl), source=-99)
          N_select_species=0                           ! include *all* obs species 
          allocate(select_species(N_select_species))   ! allocate() needed for gcc10
          nTiles_l = 0
          do iTile=1,N_catl
             halo = get_halo_around_tile(tile_coord_l(iTile), xcompact, ycompact)
             if (TileNnzObs(Obs_f_assim, halo, select_species)) then
                nTiles_l = nTiles_l + 1 ! num of tiles w/ nnz obs
                indTiles_l(nTiles_l) = l2f(iTile) ! 'full' index of tile w/ nnz obs
             end if
          end do
          call MPI_Gather(nTiles_l,1,MPI_INTEGER,                                    &
               nTilesl_vec,1,MPI_INTEGER,0,mpicomm,mpierr)
          if (root_proc) nTiles_f = sum(nTilesl_vec)
          call MPI_Bcast(nTiles_f,1,MPI_INTEGER,0,mpicomm,mpierr)
          if (logit) then
             write (tmpstr12,'(i12)') nTiles_f         ! convert integer to string
             write(logunit,'(400A)')                                                 &
                  'AnaLoadBal: Total number of tiles in EnKF analysis: ' // tmpstr12
          end if
             
          ! Step 2b: indTiles_l -> indTiles_f (on root)
          if (root_proc) then
            allocate(indTiles_f(nTiles_f), source=-99)
          else
            allocate(indTiles_f(0)) ! for debugging mode
          endif
  
          if (root_proc) then
             tmp_low_ind(1) = 1
             do iproc=1,numprocs-1
                tmp_low_ind(iproc+1) = tmp_low_ind(iproc) + nTilesl_vec(iproc)
             end do
          end if
          call MPI_Gatherv(                                                          &
               indTiles_l(1:nTiles_l), nTiles_l,                   MPI_INTEGER,      &
               indTiles_f,             nTilesl_vec, tmp_low_ind-1, MPI_INTEGER,      &
               0, mpicomm, mpierr)
          if (allocated(indTiles_l)) deallocate(indTiles_l)

          ! Step 2c: compute nTiles_ana, indTiles_f -> indTiles_ana
          quotient = nTiles_f/numprocs
          remainder = mod(nTiles_f,numprocs)

          do iproc=1,numprocs
             nTilesAna_vec(iproc) = quotient
             if (iproc<=remainder) nTilesAna_vec(iproc) = nTilesAna_vec(iproc) + 1
          end do

          nTiles_ana = nTilesAna_vec(myid+1) ! shorthand
          allocate(indTiles_ana(nTiles_ana), source=-99)
          tmp_low_ind(1) = 1
          do iproc=1,numprocs-1
             tmp_low_ind(iproc+1) = tmp_low_ind(iproc) + nTilesAna_vec(iproc)
          end do

          call MPI_Scatterv(                                            &
               indTiles_f,   nTilesAna_vec, tmp_low_ind-1, MPI_INTEGER, &
               indTiles_ana, nTiles_ana,                   MPI_INTEGER, &
               0, mpicomm, mpierr)
          if (allocated(indTiles_f)) deallocate(indTiles_f)

          ! Step 2d: indTiles_ana -> indTilesAna_vec (on root)
          ! root needs indTiles_ana from each proc to distribute cat_param, cat_progn etc.
          if (root_proc) then
             do iproc=1,numprocs
                allocate(indTilesAna_vec(iproc)%ind(nTilesAna_vec(iproc)))
             end do
          end if
          if (root_proc) then
             indTilesAna_vec(1)%ind = indTiles_ana ! copy contribution from root
             do src=1,numprocs-1
                recvct = nTilesAna_vec(src+1)
                recvtag = src
                call MPI_Recv(indTilesAna_vec(src+1)%ind,recvct,MPI_INTEGER,         &
                     src,recvtag,mpicomm,mpistatus,mpierr)
             end do
          else
             sendtag = myid
             sendct = nTiles_ana
             call MPI_Send(indTiles_ana,sendct,MPI_INTEGER,0,sendtag,mpicomm,mpierr)
          end if
          call cpu_time(t_end)

          ! Step 2: timing info
          call MPI_Reduce(t_end-t_start,tmax,1,MPI_REAL,MPI_MAX,0,mpicomm,mpierr)
          call MPI_Reduce(t_end-t_start,tmin,1,MPI_REAL,MPI_MIN,0,mpicomm,mpierr)
          if (root_proc .and. logit) write (logunit,'(2A,ES10.3,A,ES10.3)')           &
               'AnaLoadBal: Step 2 time taken (create indTiles_ana): ',               &
               '  max =', tmax, ',  min =', tmin

          ! Step 3a: for each proc create nObs_ana, indObs_ana and Obs_ana
          !          [we still have Obs_f_assim in each proc]
          call cpu_time(t_start)
          allocate(indObs_ana(N_obsf_assim), source=-99)
          allocate(tmp_ind_obs(N_obsf_assim), source=-99)
          nObs_ana = 0
          do ctr=1,nTiles_ana
             iTile = indTiles_ana(ctr) ! 'full' index
             halo = get_halo_around_tile(tile_coord_f(iTile), xcompact, ycompact)
             tmp_ind_obs = -1
             call get_ind_obs_lat_lon_box(                            &
                  N_obsf_assim,     Obs_f_assim,                      &
                  halo%minlon, halo%maxlon, halo%minlat, halo%maxlat, &
                  N_select_species, select_species,                   &   ! incl. *all* obs species (N_select_species=0)
                  N_selected_obs,   tmp_ind_obs )
             ! add N_selected_obs indices to indObs_ana. CAREFUL not to duplicate indices
             if (N_selected_obs>0) &
                  call addUniqueInts(tmp_ind_obs(1:N_selected_obs),indObs_ana,nObs_ana)
          end do
          if (allocated(tmp_ind_obs)) deallocate(tmp_ind_obs)
          if (allocated(select_species)) deallocate(select_species)
          ! sort obs indices (for layout independence)
          if (nObs_ana>1) call MAPL_Sort(indObs_ana(1:nObs_ana))

          ! Step 3b: nObs_ana -> nObsAna_vec (on root)
          call MPI_Gather(nObs_ana,1,MPI_INTEGER,                            &
               nObsAna_vec,1,MPI_INTEGER,0,mpicomm,mpierr)
          if (root_proc .and. logit) write (logunit,'(2A,I7,A,I7)')          &
               'AnaLoadBal: nObs_ana statistics:   ',                        &
               'max =', maxval(nObsAna_vec), ',  min =', minval(nObsAna_vec)
          
          ! Step 3c: indObs_ana -> indObsAna_vec (on root)
          ! root needs indObs_ana from each proc (to distribute Obs_pred_l)
          if (root_proc) then
             do iproc=1,numprocs
                allocate(indObsAna_vec(iproc)%ind(nObsAna_vec(iproc)))
             end do
          end if
          if (root_proc) then
             indObsAna_vec(1)%ind = indObs_ana(1:nObs_ana) ! copy contribution from root
             do src=1,numprocs-1
                recvct = nObsAna_vec(src+1)
                recvtag = src
                call MPI_Recv(indObsAna_vec(src+1)%ind,recvct,MPI_INTEGER,   &
                     src,recvtag,mpicomm,mpistatus,mpierr)
             end do
          else
             sendtag = myid
             sendct = nObs_ana
             call MPI_Send(indObs_ana,sendct,MPI_INTEGER,                    &
                  0,sendtag,mpicomm,mpierr)
          end if
          call cpu_time(t_end)

          ! Step 3: timing info
          call MPI_Reduce(t_end-t_start,tmax,1,MPI_REAL,MPI_MAX,0,mpicomm,mpierr)
          call MPI_Reduce(t_end-t_start,tmin,1,MPI_REAL,MPI_MIN,0,mpicomm,mpierr)
          if (root_proc .and. logit) write (logunit,'(2A,ES10.3,A,ES10.3)')            &
               'AnaLoadBal: Step 3 time taken (create indObs_ana): ',                  &
               '  max =', tmax, ',  min =', tmin
          !-AnaLodaBal-decomposition-ends-here

          !-AnaLoadBal-Input-Distribution-starts-here
          ! Step 4a: Obs_ana
          call cpu_time(t_start)
          allocate(Obs_ana(nObs_ana))
          Obs_ana = Obs_f_assim(indObs_ana(1:nObs_ana))
          if (allocated(Obs_f_assim)) deallocate(Obs_f_assim)

          ! step 4b: tile_coord_ana
          allocate(tile_coord_ana(nTiles_ana))
          tile_coord_ana = tile_coord_f(indTiles_ana)

          ! Step 4c: cat_param(N_catl) -> cat_param_f (on root) -> cat_param_ana
          if (root_proc) then
             allocate(cat_param_f(N_catf))
          else
             allocate(cat_param_f(0)) !for debugging mode
          endif
          call MPI_Gatherv(                                             &
               cat_param,   N_catl,                MPI_cat_param_type,  &
               cat_param_f, N_catl_vec, low_ind-1, MPI_cat_param_type,  &
               0, mpicomm, mpierr )
          allocate(cat_param_ana(nTiles_ana))
          if (root_proc) then
             cat_param_ana = cat_param_f(indTilesAna_vec(1)%ind)
             do dest=1,numprocs-1
                sendtag = dest
                sendct = nTilesAna_vec(dest+1)
                call MPI_Send(cat_param_f(indTilesAna_vec(dest+1)%ind), &
                     sendct,MPI_cat_param_type,                         &
                     dest,sendtag,mpicomm,mpierr)
             end do
          else
             ! source = 0
             recvtag = myid
             recvct = nTiles_ana
             call MPI_Recv(cat_param_ana,recvct,MPI_cat_param_type, &
                  0,recvtag,mpicomm,mpistatus,mpierr)
          end if
          if (allocated(cat_param_f)) deallocate(cat_param_f)

          ! Step 4d: cat_progn -> cat_progn_f (on root) -> cat_progn_ana
          ! one ensemble at a time
          if (root_proc) then
             allocate(cat_progn_f(N_catf))
          else
             allocate(cat_progn_f(0)) ! for debugging mode
          endif

          allocate(cat_progn_ana(nTiles_ana,N_ens))
          allocate(tmp_cat_progn_ana(nTiles_ana)) ! CSD-BUGFIX

          do iEns=1,N_ens
             ! cat_progn_ana(:,iEns) -> cat_progn_f (on root)
             call MPI_Gatherv(                                                      &
                  cat_progn(:,iEns),  N_catl,                  MPI_cat_progn_type,  &
                  cat_progn_f,        N_catl_vec,  low_ind-1,  MPI_cat_progn_type,  &
                  0, mpicomm, mpierr )
             if (root_proc) then
                cat_progn_ana(:,iEns) = cat_progn_f(indTilesAna_vec(1)%ind)
                do dest=1, numprocs-1
                   sendtag = dest
                   sendct = nTilesAna_vec(dest+1) ! send count
                   call MPI_Send(cat_progn_f(indTilesAna_vec(dest+1)%ind),          &
                        sendct,MPI_cat_progn_type,                                  &
                        dest,sendtag,mpicomm,mpierr)
                end do
             else
                ! source = 0 (root)
                recvtag = myid
                recvct = nTiles_ana

                ! CSD-BUGFIX (adopted by reichle, 7 Oct 2015)
                ! MPI crashed here if given a zero index first dimension, combined with a non-zero length 
                ! second dimension, and recvct=0, as it inteprets the array as having length 1.
                ! Reading into a zero length vector is OK though 
                !
                !call MPI_Recv(cat_progn_ana(:,iEns),recvct,MPI_cat_progn_type,      &
                !     0,recvtag,mpicomm,mpistatus,mpierr)
                !
                ! Solution: Use MPI_Recv with 1d array, then copy into 2d array.

                call MPI_Recv(tmp_cat_progn_ana,recvct,MPI_cat_progn_type,      &
                     0,recvtag,mpicomm,mpistatus,mpierr)
                
                cat_progn_ana(:,iEns)=tmp_cat_progn_ana
                
             end if
          end do

          if (allocated( tmp_cat_progn_ana))  deallocate(tmp_cat_progn_ana)
          if (allocated(cat_progn_f)) deallocate(cat_progn_f)

          ! Step 4e: Obs_pred_l (obs%assim=.true.) -> Obs_pred_f_assim (on root) -> Obs_pred_ana
          ! one ensemble at a time
          if (root_proc) then
             allocate(Obs_pred_f_assim(N_obsf_assim))
          else
             allocate(Obs_pred_f_assim(0)) ! for debugging mode
          endif
          allocate(Obs_pred_ana(nObs_ana,N_ens), source=0.)
          if (root_proc) then
             tmp_low_ind(1) = 1
             do iproc=1,numprocs-1
                tmp_low_ind(iproc+1) = tmp_low_ind(iproc) + N_obsl_assim_vec(iproc)
             end do
          end if
          do iEns=1,N_ens
             ! Obs_pred_l(:,iEns) -> Obs_pred_f_assim (on root) [only for obs%assim=.true.]
             call MPI_Gatherv(                                                           &
                  Obs_pred_l(ptr2indx,iEns), N_obsl_assim,                    MPI_REAL,  &
                  Obs_pred_f_assim,          N_obsl_assim_vec, tmp_low_ind-1, MPI_REAL,  &
                  0, mpicomm, mpierr )
             ! Obs_pred_f_assim (on root) -> Obs_pred_ana
             if (root_proc) then
                ! copy Obs_pred_ana for root
                Obs_pred_ana(:,iEns) = Obs_pred_f_assim(indObsAna_vec(1)%ind) 
                ! communicate
                do dest=1, numprocs-1
                   sendtag = dest
                   sendct = nObsAna_vec(dest+1) ! send count
                   call MPI_Send(Obs_pred_f_assim(indObsAna_vec(dest+1)%ind(1:sendct)),  &
                        sendct,MPI_REAL,                                                 &
                        dest,sendtag,mpicomm,mpierr)
                end do
             else
                ! source = 0 (root)
                recvtag = myid
                recvct = nObs_ana
                call MPI_Recv(Obs_pred_ana(:,iEns),recvct,MPI_REAL, &
                     0,recvtag,mpicomm,mpistatus,mpierr)
             end if
          end do
          if (allocated(Obs_pred_f_assim)) deallocate(Obs_pred_f_assim)
          if (allocated(ind_obsl_assim))   deallocate(ind_obsl_assim)
          nullify(ptr2indx)
          call cpu_time(t_end)

          ! Step 4: timing info
          call MPI_Reduce(t_end-t_start,tmax,1,MPI_REAL,MPI_MAX,0,mpicomm,mpierr)
          call MPI_Reduce(t_end-t_start,tmin,1,MPI_REAL,MPI_MIN,0,mpicomm,mpierr)
          if (root_proc .and. logit) write (logunit,'(2A,ES10.3,A,ES10.3)')            &
               'AnaLoadBal: Step 4 time taken (distribute inputs): ',                  &
               '  max =', tmax, ',  min =', tmin
          !-AnaLoadBal-Input-Distribution-ends-here

#else

          ! collect observations (within halo) that should be assimilated
          ! (ie, ONLY collect obs with flag assim==.true.)
          !
          ! NOTE: make sure to pass into get_halo_obs() only the portion
          !       of Observations_l and Obs_pred_l that are "good"
          !       [allocation of these arrays in get_obs_pred() is larger
          !        than eventual size]
          call get_halo_obs( N_ens, N_obsl,                                  &
               Observations_l(1:N_obsl), Obs_pred_l(1:N_obsl,1:N_ens),       &
               tile_coord_l, xcompact, ycompact,                             &
               N_obslH, Observations_lH, Obs_pred_lH )

#endif

          ! get observations perturbations for all ensemble members

#ifdef LDAS_MPI

          ! MPI: input full domain pert_grid_f because *_pert_ntrmdt
          !        corresponds to full domain
          !      input local nObs_ana and Obs_ana because *_pert_tile_*
          !        is diagnosed for the local domain only
          !      ALL processors MUST call this subroutine (even if nObs_ana==0)
          !        to keep Pert_rseed consistent across processors 

          allocate(Obs_pert_tmp(nObs_ana,N_ens))
          call get_obs_pert( N_ens, nObs_ana, N_obs_param,                 &
               pert_grid_f,                                                &
               obs_param, Obs_ana,                                         &
               Pert_rseed,                                                 &
               Obs_pert_tmp )

#else

          allocate(Obs_pert_tmp(N_obslH,N_ens))
          call get_obs_pert( N_ens, N_obslH, N_obs_param,                  &
               pert_grid_f,                                                &
               obs_param, Observations_lH,                                 &
               Pert_rseed,                                                 &
               Obs_pert_tmp )
#endif

          ! fill optional outputs if needed

          if (present(Obs_pert)) then

             ! Obs_pert must be returned for adaptive filter
             ! but this has not yet been implemented for MPI
             ! version

             !! Obs_pert(1:N_obslH,1:N_ens) = Obs_pert_tmp

             err_msg = 'ERROR with Obs_pert output - optional ' // &
                  'output variable Obs_pert not yet implemented'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

          end if

          ! get assimilation increments

#ifdef LDAS_MPI
          allocate(cat_progn_incr_ana(nTiles_ana,N_ens))

          call cpu_time(t_start)
          call cat_enkf_increments(                       &
               N_ens, nObs_ana, nTiles_ana, N_obs_param,  &
               update_type, obs_param,                    &
               tile_coord_ana, indTiles_ana,              & ! indTiles_ana is essentially ana2f
               Obs_ana,                                   & ! size: nObs_ana
               Obs_pred_ana,                              & ! size: (nObs_ana,N_ens)
               Obs_pert_tmp,                              &
               cat_param_ana,                             &
               xcompact, ycompact, fcsterr_inflation_fac, &
               cat_progn_ana, cat_progn_incr_ana)
          call cpu_time(t_end)

          
          ! we need this for correct timing info
          call MPI_Barrier(mpicomm, mpierr)

          ! cat_enkf_incr timinig info
          call MPI_Reduce(t_end-t_start,tmax,1,MPI_REAL,MPI_MAX,0,mpicomm,mpierr)
          call MPI_Reduce(t_end-t_start,tmin,1,MPI_REAL,MPI_MIN,0,mpicomm,mpierr)
          if (root_proc .and. logit) write (logunit,'(2A,ES10.3,A,ES10.3)')            &
               'Time taken by cat_enkf_increments:  ',                                 &
               '  max =', tmax, ',  min =', tmin
#else
          ! NOTE: make sure to pass into cat_enkf_increments() only the
          !       "valid" sub-arrays of Observations_lH and Obs_pred_lH
          !       [allocation of these arrays in get_halo_obs() is larger
          !        than eventual size]

          call cat_enkf_increments(                                     &
               N_ens, N_obslH, N_catl, N_obs_param,                     &
               update_type, obs_param,                                  &
               tile_coord_l, l2f,                                       &
               Observations_lH(1:N_obslH),                              &
               Obs_pred_lH(1:N_obslH,1:N_ens),                          &
               Obs_pert_tmp,                                            &
               cat_param,                                               &
               xcompact, ycompact, fcsterr_inflation_fac,               &
               cat_progn, cat_progn_incr )
#endif          

#ifdef LDAS_MPI
          !-AnaLoadBal-Output-Collection-starts-here-
          ! Step 5: do the reverse for cat_progn_incr
          ! cat_progn_incr_ana -> cat_progn_incr_f -> cat_progn_incr
          ! WE PROBABLY SHOULD DO AWAY WITH recvBuf
          call cpu_time(t_start)
          if (root_proc) then
             allocate(cat_progn_incr_f(N_catf))
             allocate(recvBuf(maxval(nTilesAna_vec))) ! temp storage of incoming data
          else
             allocate(cat_progn_incr_f(0)) ! for debugging
          end if
          do iEns=1,N_ens
             ! cat_progn_incr_ana -> cat_progn_incr_f
             if (root_proc) then
                do iTile=1,N_catf ! cannot do cat_progn_incr_f = 0.
                   cat_progn_incr_f(iTile) = 0. ! initialize
                end do
                ! copy contribution from root
                cat_progn_incr_f(indTilesAna_vec(1)%ind(1:nTilesAna_vec(1))) = &
                     cat_progn_incr_ana(:,iEns)
                ! communicate
                do src=1,numprocs-1
                   recvct = nTilesAna_vec(src+1)
                   recvtag = src
                   call MPI_Recv(recvBuf,recvct,MPI_cat_progn_type,                   &
                        src,recvtag,mpicomm,mpistatus,mpierr)
                   ! unpack cat_progn_incr_tmp into cat_progn_incr
                   cat_progn_incr_f(indTilesAna_vec(src+1)%ind(1:recvct))=recvBuf(1:recvct)
                end do
             else
                sendtag = myid
                sendct = nTiles_ana
                call MPI_Send(cat_progn_incr_ana(:,iEns),sendct,MPI_cat_progn_type,   &
                     0,sendtag,mpicomm,mpierr)
             end if
             ! cat_progn_incr_f -> cat_progn_incr
             call MPI_Scatterv(                                                       &
                  cat_progn_incr_f,       N_catl_vec, low_ind-1, MPI_cat_progn_type,  &
                  cat_progn_incr(:,iEns), N_catl,                MPI_cat_progn_type,  &
                  0, mpicomm, mpierr)
          end do
          if (allocated(recvBuf)) deallocate(recvBuf)
          if (allocated(cat_progn_incr_f)) deallocate(cat_progn_incr_f)
          call cpu_time(t_end)

          ! Step 5: timing info
          call MPI_Reduce(t_end-t_start,tmax,1,MPI_REAL,MPI_MAX,0,mpicomm,mpierr)
          call MPI_Reduce(t_end-t_start,tmin,1,MPI_REAL,MPI_MIN,0,mpicomm,mpierr)
          if (root_proc .and. logit) write (logunit,'(2A,ES10.3,A,ES10.3)')           &
               'AnaLoadBal: Step 5 time taken (collect increments): ',                &
               '  max =', tmax, ',  min =', tmin
          !-AnaLoadBal-Output-Collection-ends-here-
#endif

          ! set flag that fresh increments are available

          fresh_incr = .true.

          ! cleanup

          if (allocated( Obs_pert_tmp))       deallocate(Obs_pert_tmp)

#ifdef LDAS_MPI
          if (allocated( cat_progn_incr_ana)) deallocate(cat_progn_incr_ana)
          if (allocated( Obs_pred_ana))       deallocate(Obs_pred_ana)
          do iproc=1,numprocs
             if (allocated(indObsAna_vec(iproc)%ind))     &
                  deallocate(indObsAna_vec(iproc)%ind)
          end do
          if (allocated( Obs_ana))            deallocate(Obs_ana)
          if (allocated( indObs_ana))         deallocate(indObs_ana)
          if (allocated( cat_progn_ana))      deallocate(cat_progn_ana)
          if (allocated( cat_param_ana))      deallocate(cat_param_ana)
          do iproc=1,numprocs
             if (allocated(indTilesAna_vec(iproc)%ind))   &
                  deallocate(indTilesAna_vec(iproc)%ind)
          end do
          if (associated(tile_coord_ana))     deallocate(tile_coord_ana)
          if (allocated( indTiles_ana))       deallocate(indTiles_ana)
#else
          if (associated(Obs_pred_lH))        deallocate(Obs_pred_lH)

#endif

       end if   ! (N_obsf>0) .and. assimflag

       ! --------------------------------------------------------------
       !
       ! Obs bias
       !
       ! B4. Update the bias tcount for latest obs

       if (N_obsbias_max > 0)                                           &
            call obs_bias_upd_tcount(date_time,dtstep_assim,            &
            N_catf, N_catl, N_obsl, N_obs_param, N_obsbias_max, f2l,    &
            obs_param, Observations_l(1:N_obsl),                        &
            obs_bias)

       ! --------------------------------------------------------------

       ! clean up

       if (associated(Obs_pred_l ))            deallocate(Obs_pred_l)

       if (associated(N_tile_in_cell_ij_f))    deallocate(N_tile_in_cell_ij_f)
       if (associated(tile_num_in_cell_ij_f))  deallocate(tile_num_in_cell_ij_f)

       ! write assimilated (scaled) SMAP Tb observations and select forecast states
       ! into SMAP L4_SM aup file

       if (out_smapL4SMaup)                                                          &
            call write_smapL4SMaup( 'obs_fcst', date_time, exp_id, N_ens, &
            N_catl, N_catf, N_obsl, tile_coord_f, tile_grid_g, N_catl_vec, low_ind,  &
            N_obs_param, obs_param, Observations_l, cat_param, cat_progn       )

    end if  ! end if (.true.)
       
    call date_and_time(date_string, time_string)

    if (logit) write (logunit,*) 'get_enkf_increments(): exit at ', &
         date_string, ', ', time_string 

  end subroutine get_enkf_increments

  ! ********************************************************************

  subroutine addUniqueInts(src, dest, n_dest)

    ! add integers from the src array to the dest arr
    ! provided that they don't already exist in dest
    ! upon input, n_dest is the number of existing
    ! entries. both n_dest and dest are updated

    ! input/output
    integer, intent(in) :: src(:)
    integer, intent(inout) :: n_dest
    integer, intent(inout) :: dest(:)

    ! local
    integer :: iInt, nInts

    nInts = size(src)

    if (n_dest==0) then
       dest(1:nInts) = src
       n_dest = nInts
    else
       do iInt=1,nInts
          if(any(dest(1:n_dest)==src(iInt))) then
             ! src(iInt) exists in dset(1:n_dset)
             cycle
          else
             n_dest = n_dest + 1
             dest(n_dest) = src(iInt)
          end if
       end do
    end if

  end subroutine addUniqueInts

  ! ********************************************************************

  subroutine apply_enkf_increments( N_catd, N_ens, update_type, &
       cat_param, cat_progn_incr, cat_progn )

    implicit none

    ! reichle, 19 Oct 2005

    ! removed call to subroutine recompute_diagnostics()
    ! -reichle+csdraper, 30 Oct 2013

    integer, intent(in) :: N_catd, N_ens, update_type

    type(cat_param_type), dimension(N_catd),       intent(in)    :: cat_param

    type(cat_progn_type), dimension(N_catd,N_ens), intent(in)    :: cat_progn_incr

    type(cat_progn_type), dimension(N_catd,N_ens), intent(inout) :: cat_progn

    ! -----------------

    integer :: n, n_e

    logical :: cat_progn_has_changed, check_snow

    character(len=*), parameter :: Iam = 'apply_enkf_increments'

    ! ----------------------------------------------------------------
    !
    ! apply increments

    cat_progn_has_changed = .true.     ! conservative initialization
    
    check_snow            = .true.     ! conservative initialization
    
    select_update_type: select case (update_type)
       
    case (1,2) select_update_type  ! soil moisture update
       
       if (logit) write (logunit,*) &
            'apply_enkf_increments(): applying soil moisture increments'
       
       do n=1,N_catd
          do n_e=1,N_ens
             
             cat_progn(n,n_e)%srfexc = &
                  cat_progn(n,n_e)%srfexc + cat_progn_incr(n,n_e)%srfexc
             cat_progn(n,n_e)%rzexc = &
                  cat_progn(n,n_e)%rzexc  + cat_progn_incr(n,n_e)%rzexc
             cat_progn(n,n_e)%catdef = &
                  cat_progn(n,n_e)%catdef + cat_progn_incr(n,n_e)%catdef

          end do
       end do

       cat_progn_has_changed = .true.

    case (3) select_update_type    ! Tskin update

       if (logit) write (logunit,*) &
            'apply_enkf_increments(): NOT applying Tskin increments'

       cat_progn_has_changed = .false.

    case (4,7) select_update_type  ! Tskin/ght1 update

       if (logit) write (logunit,*) &
            'apply_enkf_increments(): applying Tskin/ght1 increments'

       do n=1,N_catd
          do n_e=1,N_ens

             cat_progn(n,n_e)%tc1 = &
                  cat_progn(n,n_e)%tc1 + cat_progn_incr(n,n_e)%tc1
             cat_progn(n,n_e)%tc2 = &
                  cat_progn(n,n_e)%tc2 + cat_progn_incr(n,n_e)%tc2
             cat_progn(n,n_e)%tc4 = &
                  cat_progn(n,n_e)%tc4 + cat_progn_incr(n,n_e)%tc4
	     cat_progn(n,n_e)%ght(1) = &
                  cat_progn(n,n_e)%ght(1) + cat_progn_incr(n,n_e)%ght(1)

          end do
       end do

       cat_progn_has_changed = .true.

    case (5) select_update_type    ! Tskin/ght1 update

       if (logit) write (logunit,*) &
            'apply_enkf_increments(): NOT applying Tskin/ght1 increments'

       cat_progn_has_changed = .false.

    case (6,8,9,10) select_update_type    ! soil moisture and temperature update

       ! for update_type 10, catdef increments may be zero by design       
       
       if (logit) write (logunit,*) &
            'apply_enkf_increments(): applying soil moisture and Tskin/ght1 increments'
       
       do n=1,N_catd
          do n_e=1,N_ens
             
             cat_progn(n,n_e)%srfexc = &
                  cat_progn(n,n_e)%srfexc + cat_progn_incr(n,n_e)%srfexc
             cat_progn(n,n_e)%rzexc = &
                  cat_progn(n,n_e)%rzexc  + cat_progn_incr(n,n_e)%rzexc
             cat_progn(n,n_e)%catdef = &
                  cat_progn(n,n_e)%catdef + cat_progn_incr(n,n_e)%catdef
             
             cat_progn(n,n_e)%tc1 = &
                  cat_progn(n,n_e)%tc1    + cat_progn_incr(n,n_e)%tc1
             cat_progn(n,n_e)%tc2 = &
                  cat_progn(n,n_e)%tc2    + cat_progn_incr(n,n_e)%tc2
             cat_progn(n,n_e)%tc4 = &
                  cat_progn(n,n_e)%tc4    + cat_progn_incr(n,n_e)%tc4

             cat_progn(n,n_e)%ght(1) = &
                  cat_progn(n,n_e)%ght(1) + cat_progn_incr(n,n_e)%ght(1)

          end do
       end do

       cat_progn_has_changed = .true.

       check_snow            = .false.  ! turn off for now to maintain 0-diff w/ SMAP Tb DA test case

    case default

       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown update_type')

    end select select_update_type

    ! ------------------------------------------------------------------
    !
    ! check (and possibly fix) cat_progn

    if (cat_progn_has_changed) then

       do n_e=1,N_ens

          call check_cat_progn( check_snow, N_catd, cat_param, cat_progn(:,n_e) )

       end do

    end if

  end subroutine apply_enkf_increments

  ! ********************************************************************

  subroutine output_ObsFcstAna(date_time, exp_id, &
       N_obsl, Observations_l, N_obs_param, rf2f)

    ! obs space output: observations, obs space forecast, obs space analysis, and
    ! associated error variances
    !
    ! - reichle, 16 Jun 2011

    implicit none
    
    type(date_time_type),                         intent(in) :: date_time
    
    character(*),                                 intent(in) :: exp_id

    integer,                                      intent(in) :: N_obsl, N_obs_param


    type(obs_type),       dimension(N_obsl),      intent(in) :: Observations_l

    integer,              dimension(:), optional, intent(in) :: rf2f

    ! ---------------------

    ! locals

    character(40),  parameter                 :: file_tag = 'ldas_ObsFcstAna'
    character(40),  parameter                 :: dir_name = 'ana'

    type(obs_type), dimension(:), allocatable :: Observations_f, Observations_tmp

    integer                                   :: n, N_obsf
    integer,        dimension(:), allocatable :: rf_tilenums, tilenums

    integer,        dimension(numprocs)       :: N_obsl_vec, tmp_low_ind

    character(300)                            :: fname

#ifdef LDAS_MPI

    integer                                   :: this_species, ind_tmp, j

    integer,        dimension(N_obs_param)    :: N_per_species
    integer,        dimension(N_obs_param)    :: species_low_ind
    integer,        dimension(N_obs_param)    :: ind_within_species

#endif

    ! --------------------------------------------------------------------

    if (logit) write (logunit,*) 'writing ' // trim(file_tag) //' file'

    ! ---------------------------------------------
    !
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
    else
      allocate(Observations_f(0))
    end if

#ifdef LDAS_MPI

    call MPI_GATHERV(                                                    &
         Observations_l,  N_obsl,                      MPI_obs_type,     &
         Observations_f,  N_obsl_vec,  tmp_low_ind-1,  MPI_obs_type,     &
         0, mpicomm, mpierr )

#else

    Observations_f = Observations_l

#endif

    ! --------------------------------------------------------------
    !
    ! write to file

    if (root_proc) then

#ifdef LDAS_MPI

       ! sort observations according to species so that the ObsFcstAna
       ! file looks the same regardless of the MPI configuration
       !
       ! NOTE: within each species, observations should already be sorted according
       !       to tilenum because of the loop through species in collect_obs()

       ! first loop: count how many obs there are per species

       N_per_species(1:N_obs_param) = 0

       do n=1,N_obsf

          this_species = Observations_f(n)%species

          N_per_species(this_species) = N_per_species(this_species) + 1

       end do

       ! get starting index of each species in sorted vector Observations_f

       species_low_ind(1) = 1

       do j=1,N_obs_param-1

          species_low_ind(j+1) = species_low_ind(j) + N_per_species(j)

       end do

       ! second loop: re-order Observations_f into Observations_tmp

       allocate(Observations_tmp(N_obsf))

       ind_within_species(1:N_obs_param) = 0

       do n=1,N_obsf

          this_species = Observations_f(n)%species

          ind_tmp = species_low_ind(this_species) + ind_within_species(this_species)

          Observations_tmp(ind_tmp) = Observations_f(n)

          ! NOTE: ind_within_species is "zero-based"

          ind_within_species(this_species) = ind_within_species(this_species) + 1

       end do

       ! copy back into Observations_f

       do n=1,N_obsf

          Observations_f(n) = Observations_tmp(n)

       end do


       deallocate(Observations_tmp)

#endif ! LDAS_MPI
       
       ! reorder tilenum, so it is consistent with the order in tile_coord.bin file
       if(present(rf2f)) then
          allocate(rf_tilenums(N_obsf), tilenums(N_obsf))
          rf_tilenums = Observations_f(:)%tilenum
          tilenums = rf2f(rf_tilenums)
          Observations_f(:)%tilenum =tilenums
          deallocate(rf_tilenums, tilenums)
       endif
       
       ! write to file

       fname = get_io_filename( './', exp_id, file_tag, date_time=date_time, &
            dir_name=dir_name, ens_id=-1, no_subdirs=.true. )
         
       open( 10, file=fname, form='unformatted', action='write')

       ! write header

       write (10) N_obsf, date_time%year, date_time%month, &
            date_time%day, date_time%hour, date_time%min, date_time%sec, &
            date_time%dofyr, date_time%pentad

       ! write data

       ! Assuming a linear model and uncorrelated obs/model errors,
       !
       ! the expected var of OminusF is   HPHt + R,   and
       ! the expected var of OminusA is   R - HAHt,   where
       !
       !   P = prior state error covariance
       !   A = posterior state error covariance
       !   H = observation operator

       write (10) (Observations_f(n)%assim,   n=1,N_obsf)
       write (10) (Observations_f(n)%species, n=1,N_obsf)

       write (10) (Observations_f(n)%tilenum, n=1,N_obsf)

       write (10) (Observations_f(n)%lon,     n=1,N_obsf)
       write (10) (Observations_f(n)%lat,     n=1,N_obsf)

       write (10) (Observations_f(n)%obs,     n=1,N_obsf)
       write (10) (Observations_f(n)%obsvar,  n=1,N_obsf)     ! R

       write (10) (Observations_f(n)%fcst,    n=1,N_obsf)
       write (10) (Observations_f(n)%fcstvar, n=1,N_obsf)     ! HPHt

       write (10) (Observations_f(n)%ana,     n=1,N_obsf)
       write (10) (Observations_f(n)%anavar,  n=1,N_obsf)     ! HAHt

       close(10,status='keep')

    end if
    if (allocated(Observations_f)) deallocate(Observations_f)

  end subroutine output_ObsFcstAna

  ! **********************************************************************

  subroutine output_ObsFcstAna_wrapper( out_ObsFcstAna,                      &
       date_time, exp_id,                                                    &
       N_obsl, N_obs_param, N_ens,                                           &
       N_catl, tile_coord_l,                                                 &
       N_catf, tile_coord_f, pert_grid_g,                                    &
       N_catl_vec, low_ind, f2l,                                             &
       obs_param,                                                            &
       met_force, lai, cat_param, cat_progn, mwRTM_param,                    &
       Observations_l, rf2f )

    implicit none

    ! reichle, 5 Jun 2006

    ! changed intent of Observations for adaptive filtering
    ! - reichle, 15 Dec 2006

    ! major revisions for new obs handling and MPI

    logical,                intent(in) :: out_ObsFcstAna


    type(date_time_type),   intent(in) :: date_time

    character(len=*),       intent(in) :: exp_id

    integer,                intent(in) :: N_obsl, N_obs_param, N_ens, N_catl, N_catf

    type(tile_coord_type),  dimension(:),     pointer :: tile_coord_l  ! input
    type(tile_coord_type),  dimension(:),     pointer :: tile_coord_f  ! input

    type(grid_def_type),                              intent(in) :: pert_grid_g

    integer,                dimension(numprocs),      intent(in) :: N_catl_vec, low_ind

    integer,                dimension(N_catf),        intent(in) :: f2l

    type(obs_param_type),   dimension(N_obs_param),   intent(in) :: &
         obs_param

    type(met_force_type),   dimension(N_catl),        intent(in)    :: met_force

    real,                   dimension(N_catl),        intent(in)    :: lai

    type(cat_param_type),   dimension(N_catl),        intent(in)    :: cat_param
    type(cat_progn_type),   dimension(N_catl,N_ens),  intent(in)    :: cat_progn

    type(mwRTM_param_type), dimension(N_catl),        intent(in)    :: mwRTM_param


    type(obs_type),         dimension(:),     pointer :: Observations_l ! inout

    integer,                dimension(N_catf), optional, intent(in) :: rf2f ! re-ordered to LDASsa 

    ! local variables

    real,    dimension(:,:),   pointer :: Obs_pred_l            => null()

    integer :: N_obsl_tmp


    character(len=*), parameter :: Iam = 'output_ObsFcstAna_wrapper'

    ! --------------------------------------------------------------

    nullify(Obs_pred_l)

    ! output "O-A" (obs - analysis) whenever innovations are output

    if (out_ObsFcstAna) then

       ! compute model forecast of observations

       N_obsl_tmp = N_obsl ! cannot pass N_obsl into get_obs_pred() b/c of intent(in)

       call get_obs_pred(                                      &
            .false.,                                           & ! -> after EnKF update
            N_obs_param, N_ens,                                &
            N_catl, tile_coord_l,                              &
            N_catf, tile_coord_f, f2l,                         &
            N_catl_vec, low_ind, pert_grid_g,                  &
            obs_param,                                         &
            met_force, lai, cat_param, cat_progn, mwRTM_param, &
            N_obsl_tmp, Observations_l, Obs_pred_l )

       ! clean up

       if (associated(Obs_pred_l))             deallocate(Obs_pred_l)

       ! write out model, observations, and "OminusA" information

       call output_ObsFcstAna( date_time, exp_id, N_obsl, &
            Observations_l(1:N_obsl), N_obs_param, rf2f=rf2f )

    end if

  end subroutine output_ObsFcstAna_wrapper

  ! **********************************************************************

  subroutine output_smapL4SMaup( date_time, work_path, exp_id, dtstep_assim,    &
       N_ens, N_catl, N_catf, N_obsl_max,                               &
       tile_coord_f, tile_grid_g, pert_grid_f,                                  &
       N_catl_vec, low_ind, l2f, N_tile_in_cell_ij_f, tile_num_in_cell_ij_f,    &
       N_obs_param, obs_param, Observations_l, cat_param, cat_progn )

    ! wrapper for output of original (unscaled) SMAP Tb observations into
    ! SMAP L4_SM "aup" (analysis update)
    !
    ! see subroutine write_smapL4SMaup() for details
    !
    ! reichle,  2 May 2013
    ! reichle, 11 Dec 2013 - added 'SMOS_fit_Tb*' obs species
    !                        (needed for L4_SM_SMOS prototype product)
    !
    ! -------------------------------------------------------------------

    implicit none

    type(date_time_type),  intent(in) :: date_time

    character(*),          intent(in) :: work_path
    character(*),          intent(in) :: exp_id

    integer,               intent(in) :: dtstep_assim
    integer,               intent(in) :: N_ens,  N_catl,     N_catf
    integer,               intent(in) :: N_obsl_max, N_obs_param

    type(tile_coord_type), dimension(:),     pointer :: tile_coord_f  ! input

    type(grid_def_type),                             intent(in) :: tile_grid_g
    type(grid_def_type),                             intent(in) :: pert_grid_f

    integer,               dimension(numprocs),      intent(in) :: N_catl_vec
    integer,               dimension(numprocs),      intent(in) :: low_ind

    integer,               dimension(N_catl),        intent(in) :: l2f

    ! N_tile_in_cell_ij and tile_num_in_cell_ij are on the "full" domain
    !  and guaranteed to be allocated ONLY for the root_proc
    !  (but may be allocated on all processors depending on obs_param%FOV)

    integer, dimension(:,:),   pointer :: N_tile_in_cell_ij_f   ! input
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij_f ! input

    type(obs_param_type),  dimension(N_obs_param),   intent(in) :: obs_param

    type(obs_type),        dimension(:),     pointer :: Observations_l ! input

    type(cat_param_type),  dimension(N_catl),        intent(in) :: cat_param
    type(cat_progn_type),  dimension(N_catl,N_ens),  intent(in) :: cat_progn

    ! ------------------------------------------

    ! local variables

    logical :: found_obs_f

    integer :: j, n, N_obs_param_tmp, N_obsl_tmp

    type(obs_param_type),  dimension(N_obs_param) :: obs_param_tmp

    ! ---------------------------------------------------------------------

    ! create custom obs_param_tmp that contains only relevant SMAP obs
    ! and disables scaling

    j = 0

    do n=1,N_obs_param

       select case (trim(obs_param(n)%descr))

       case('SMAP_L2AP_Tbh_D','SMAP_L2AP_Tbv_D',    &
            'SMAP_L2AP_Tbh_A','SMAP_L2AP_Tbv_A',    &
            'SMAP_L1C_Tbh_D', 'SMAP_L1C_Tbv_D',     &
            'SMAP_L1C_Tbh_A', 'SMAP_L1C_Tbv_A',     &
            'SMAP_L1C_Tbh_E09_D', 'SMAP_L1C_Tbv_E09_D', &
            'SMAP_L1C_Tbh_E09_A', 'SMAP_L1C_Tbv_E09_A', &
            'SMAP_L1C_Tbh_E27_D', 'SMAP_L1C_Tbv_E27_D', &
            'SMAP_L1C_Tbh_E27_A', 'SMAP_L1C_Tbv_E27_A', &
            'SMOS_fit_Tbh_D', 'SMOS_fit_Tbv_D',     &
            'SMOS_fit_Tbh_A', 'SMOS_fit_Tbv_A'      &
            )

          j = j + 1

          obs_param_tmp(j) = obs_param(n)

          obs_param_tmp(j)%scale = .false.  ! disable scaling

       case default

          ! do nothing

       end select

    end do

    N_obs_param_tmp = j

    ! collect relevant observations  (note: out_obslog=.false. in this case)

    call collect_obs(                                                            &
         work_path, exp_id, date_time, dtstep_assim,                             &
         N_catl,                                                                 &
         N_catf, tile_coord_f, pert_grid_f,                                      &
         N_tile_in_cell_ij_f, tile_num_in_cell_ij_f,                             &
         N_catl_vec, low_ind, l2f,                                               &
         N_obs_param_tmp, obs_param_tmp(1:N_obs_param_tmp), N_obsl_max, .false., &
         N_obsl_tmp, Observations_l, found_obs_f )

    ! write appropriate fields (according to 'option') into file

    call write_smapL4SMaup( 'orig_obs', date_time, exp_id, N_ens,     &
         N_catl, N_catf, N_obsl_tmp, tile_coord_f, tile_grid_g,                  &
         N_catl_vec, low_ind,                                                    &
         N_obs_param_tmp, obs_param_tmp(1:N_obs_param_tmp), Observations_l,      &
         cat_param, cat_progn       )

  end subroutine output_smapL4SMaup

  ! **********************************************************************

  subroutine write_smapL4SMaup( option, date_time, exp_id, N_ens,    &
       N_catl, N_catf, N_obsl, tile_coord_f, tile_grid_g, N_catl_vec, low_ind,  &
       N_obs_param, obs_param, Observations_l, cat_param, cat_progn       )

    ! output of custom collection for SMAP L4_SM "aup" (analysis update)
    !
    ! can be used with "SMAP_L*_Tb*" or "SMOS_fit_Tb*" obs species, but *not* with
    ! both simultaneously, which is checked in subroutine read_ens_upd_inputs()
    !
    ! "aup" files are written in three stages (controlled by "option"):
    !
    ! option = 'orig_obs' : write original obs (before scaling) into output file
    !
    !                         - tb_h_obs_time_sec                ! real*8
    !                         - tb_v_obs_time_sec                ! real*8
    !                         - tb_h_resolution_flag             ! integer
    !                         - tb_v_resolution_flag             ! integer
    !                         - tb_h_orbit_flag                  ! integer
    !                         - tb_v_orbit_flag                  ! integer
    !                         - tb_h_obs
    !                         - tb_v_obs
    !
    ! option = 'obs_fcst' : append assimilated (scaled) obs and select fcst
    !                       fields into output file
    !
    !                         - tb_h_obs_assim
    !                         - tb_v_obs_assim
    !                         - tb_h_obs_errstd
    !                         - tb_v_obs_errstd
    !
    !                         - tb_h_forecast
    !                         - tb_v_forecast
    !                         - tb_h_forecast_ensstd
    !                         - tb_v_forecast_ensstd
    !
    !                         - sm_surface_forecast
    !                         - sm_rootzone_forecast
    !                         - sm_profile_forecast
    !                         - surface_temp_forecast
    !                         - soil_temp_layer1_forecast
    
    ! option = 'analysis' : append select analysis fields into output file
    !
    !                         - sm_surface_analysis
    !                         - sm_rootzone_analysis
    !                         - sm_profile_analysis
    !                         - surface_temp_analysis
    !                         - soil_temp_layer1_analysis
    !
    !                         - sm_surface_analysis_ensstd
    !                         - sm_rootzone_analysis_ensstd
    !                         - sm_profile_analysis_ensstd
    !                         - surface_temp_analysis_ensstd
    !                         - soil_temp_layer1_analysis_ensstd
    !
    !
    ! reichle, 26 Apr 2013
    ! reichle, 11 Dec 2013 - added 'SMOS_fit_Tb*' obs species
    !                        (needed for L4_SM_SMOS prototype product)
    ! reichle,  3 Feb 2014 - added output: "tb_[h/v]_obs_time_sec", "tb_[h/v]_orbit_flag"
    ! reichle, 20 Nov 2014 - changed units of soil moisture output to [m3/m3]
    ! reichle,  6 Jun 2016 - added Tb forecast for obs that cannot be scaled (and are not
    !                         assimilated)
    !
    ! ------------------------------------------------------------------------------

    implicit none

    character(*),          intent(in) :: option

    type(date_time_type),  intent(in) :: date_time

    character(*),          intent(in) :: exp_id

    integer,               intent(in) :: N_ens, N_catl, N_catf
    integer,               intent(in) :: N_obsl, N_obs_param

    type(tile_coord_type), dimension(:),     pointer :: tile_coord_f  ! input

    type(grid_def_type),                             intent(in) :: tile_grid_g

    integer,               dimension(numprocs),      intent(in) :: N_catl_vec
    integer,               dimension(numprocs),      intent(in) :: low_ind

    type(obs_param_type),  dimension(N_obs_param),   intent(in) :: obs_param

    type(obs_type),        dimension(:),     pointer :: Observations_l ! input

    type(cat_param_type),  dimension(N_catl),        intent(in) :: cat_param
    type(cat_progn_type),  dimension(N_catl,N_ens),  intent(in) :: cat_progn

    ! --------------------------------------------------------------

    ! local variables

    integer,       parameter :: unitnum  = 10

    character(40), parameter :: file_tag = 'ldas_tile_inst_smapL4SMaup'
    character(40), parameter :: dir_name = 'ana'

    logical                  :: use_real8, is_orbit, is_obsassim

    integer                  :: j, k, m, n, n_e, kmax
    integer                  :: N_obsf, this_species

    real                     :: this_lat, this_lon, col_ind, row_ind
    real                     :: nodatavalue, tol, existing_orbflag

    character(300)           :: fname
    character(40)            :: position

    integer,        dimension(numprocs)         :: N_obsl_vec, tmp_low_ind

    type(obs_type), dimension(:),   allocatable :: Observations_f

    integer,        dimension(:),   allocatable :: col_beg_9km, col_end_9km
    integer,        dimension(:),   allocatable :: row_beg_9km, row_end_9km

    real,           dimension(:),   allocatable :: Obs_f_resflag
    real,           dimension(:),   allocatable :: Obs_f_orbflag
    real,           dimension(:),   allocatable :: Obs_f_tmpdata
    real*8,         dimension(:),   allocatable :: Obs_f_tmpdata_8

    real,           dimension(:,:), allocatable :: data_h_9km_grid,   data_v_9km_grid
    real,           dimension(:),   allocatable :: data_h_9km_tile,   data_v_9km_tile

    real*8,         dimension(:,:), allocatable :: data_h_9km_grid_8, data_v_9km_grid_8
    real*8,         dimension(:),   allocatable :: data_h_9km_tile_8, data_v_9km_tile_8

    integer,        dimension(:,:), allocatable :: ndata_h_9km_grid,  ndata_v_9km_grid

    real,           dimension(:),   allocatable :: tile_data_f

    real, dimension(      N_catl)        :: srfexc, rzexc, catdef
    real, dimension(      N_catl)        :: ar1,    ar2,   ar4

    real, dimension(N_gt, N_catl)        :: tp

    real, dimension(      N_catl, N_ens) :: sfmc,   rzmc,  prmc, tsurf, tp1

    real, dimension(      N_catl, 5)     :: tile_mean_l, tile_std_l

    character(len=*), parameter :: Iam = 'write_smapL4SMaup'
    character(len=400) :: err_msg
    character(len=10)  :: gridname_tmp

    ! --------------------------------------------------------------
    !
    ! smapL4SMaup output only works for 9 km EASE grids

    if ( index(tile_grid_g%gridtype, 'M09') == 0  ) then
       err_msg = 'out_smapL4SMaup requires tile-space for 9 km EASEv1 or EASEv2 grid'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! --------------------------------------------------------------
    !
    ! assemble file name and open file

    if (root_proc) then

       fname = get_io_filename( './', exp_id, file_tag,                      &
            date_time=date_time, dir_name=dir_name, ens_id=-1, no_subdirs=.true.)

       if     (option=='orig_obs')                         then

          position='rewind'          ! open file at the beginning

       elseif (option=='obs_fcst' .or. option=='analysis') then

          position='append'          ! append to file

       else

          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown option')

       end if

       open(unitnum, file=fname, form='unformatted', action='write',  &
            position=position)

    end if

    ! --------------------------------------------------------------
    !
    ! observation-space output

    if (option=='orig_obs' .or. option=='obs_fcst') then

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
       else
         allocate(Observations_f(0))
       end if

#ifdef LDAS_MPI

       call MPI_GATHERV(                                                    &
            Observations_l,  N_obsl,                      MPI_obs_type,     &
            Observations_f,  N_obsl_vec,  tmp_low_ind-1,  MPI_obs_type,     &
            0, mpicomm, mpierr )

#else

       Observations_f = Observations_l

#endif

       ! -----------------------------------------------------

       if (root_proc) then

          ! determine mapping from Observations vector onto global 9 km EASE grid

          allocate(col_beg_9km(  N_obsf))
          allocate(col_end_9km(  N_obsf))
          allocate(row_beg_9km(  N_obsf))
          allocate(row_end_9km(  N_obsf))

          allocate(Obs_f_resflag(N_obsf))
          allocate(Obs_f_orbflag(N_obsf))

          col_beg_9km   = -8888
          col_end_9km   = -8888
          row_beg_9km   = -8888
          row_end_9km   = -8888

          Obs_f_resflag = nodata_generic
          Obs_f_orbflag = nodata_generic

          do n=1,N_obsf
             
             this_species = Observations_f(n)%species
             
             this_lat     = Observations_f(n)%lat
             this_lon     = Observations_f(n)%lon
             
             ! obtain i,j index range of this obs on global 9 km EASE grid
             
             select case (trim(obs_param(this_species)%descr))
                
             case('SMAP_L2AP_Tbh_D',    'SMAP_L2AP_Tbv_D',      &
                  'SMAP_L2AP_Tbh_A',    'SMAP_L2AP_Tbv_A',      &
                  'SMAP_L1C_Tbh_E09_D', 'SMAP_L1C_Tbv_E09_D',   &
                  'SMAP_L1C_Tbh_E09_A', 'SMAP_L1C_Tbv_E09_A'    &
                  )
               
                if (index(tile_grid_g%gridtype, 'M09') /=0) then 
                   call  ease_convert(trim(tile_grid_g%gridtype), this_lat, this_lon, col_ind, row_ind)                
                endif

                ! col_ind and row_ind are zero-based, need one-based index here
                
                col_beg_9km(n) = nint(col_ind)+1
                col_end_9km(n) = col_beg_9km(n)
                
                row_beg_9km(n) = nint(row_ind)+1
                row_end_9km(n) = row_beg_9km(n)
                
                Obs_f_resflag(n) = 2
                Obs_f_orbflag(n) = obs_param(this_species)%orbit
                
             case('SMAP_L1C_Tbh_E27_D', 'SMAP_L1C_Tbv_E27_D',   &
                  'SMAP_L1C_Tbh_E27_A', 'SMAP_L1C_Tbv_E27_A'    &
                  )

                if (index(tile_grid_g%gridtype, 'M09') /=0) then
                   call  ease_convert(trim(tile_grid_g%gridtype), this_lat, this_lon, col_ind, row_ind)             
                endif                
                
                ! col_ind and row_ind are zero-based, need one-based index here
                ! L1C E27 spacing is one every three in each direction (~27-km spacing)
                
                col_beg_9km(n) = max( (nint(col_ind)-1)+1,    1)
                col_end_9km(n) = min( col_beg_9km(n)+2,    3856)
                
                row_beg_9km(n) = max( (nint(row_ind)-1)+1,    1)
                row_end_9km(n) = min( row_beg_9km(n)+2,    1624)
                
                Obs_f_resflag(n) = 3
                Obs_f_orbflag(n) = obs_param(this_species)%orbit
                
             case('SMAP_L1C_Tbh_D', 'SMAP_L1C_Tbv_D',   &
                  'SMAP_L1C_Tbh_A', 'SMAP_L1C_Tbv_A',   &
                  'SMOS_fit_Tbh_D', 'SMOS_fit_Tbv_D',   &
                  'SMOS_fit_Tbh_A', 'SMOS_fit_Tbv_A'    &
                  )
                
                if (index(tile_grid_g%gridtype, 'M09') /=0) then
                   ! subindex (1:7) to get the string EASEvx_
                   gridname_tmp = tile_grid_g%gridtype(1:7)//'M36'
                   call  ease_convert(gridname_tmp, this_lat, this_lon, col_ind, row_ind)                
                endif

                ! col_ind and row_ind are zero-based, need one-based index here
                
                col_beg_9km(n) =  nint(col_ind)   *4 + 1
                col_end_9km(n) = (nint(col_ind)+1)*4
                
                row_beg_9km(n) =  nint(row_ind)   *4 + 1
                row_end_9km(n) = (nint(row_ind)+1)*4
                
                Obs_f_resflag(n) = 1
                Obs_f_orbflag(n) = obs_param(this_species)%orbit
                
             case default

                ! do nothing
                
             end select
             
          end do     ! loop through Observations_f

          ! ----------------------------------------

          ! map Observations%[xx] fields onto global 9km EASE grid,

          allocate(ndata_h_9km_grid(tile_grid_g%N_lon,tile_grid_g%N_lat))
          allocate(ndata_v_9km_grid(tile_grid_g%N_lon,tile_grid_g%N_lat))

          if (option=='orig_obs') then

             kmax = 4

             ! for "orig_obs" and k==1, need real*8
             !  (switch to regular real within "k=1,kmax" loop when "k==2")

             use_real8 = .true.

             allocate(Obs_f_tmpdata_8(N_obsf))

             allocate(data_h_9km_grid_8(tile_grid_g%N_lon,tile_grid_g%N_lat))
             allocate(data_v_9km_grid_8(tile_grid_g%N_lon,tile_grid_g%N_lat))

             allocate(data_h_9km_tile_8(N_catf))
             allocate(data_v_9km_tile_8(N_catf))

          else  ! "obs_fcst"

             kmax = 4

             use_real8 = .false.

             allocate(Obs_f_tmpdata(N_obsf))

             allocate(data_h_9km_grid(tile_grid_g%N_lon,tile_grid_g%N_lat))
             allocate(data_v_9km_grid(tile_grid_g%N_lon,tile_grid_g%N_lat))

             allocate(data_h_9km_tile(N_catf))
             allocate(data_v_9km_tile(N_catf))

          end if

          ! loop through individual output fields

          do k=1,kmax

             if (option=='orig_obs') then

                if (k==2) then

                   ! switch allocatable arrays to default type for real
                   !  (real*8 were needed only for k==1)

                   use_real8 = .false.

                   deallocate(Obs_f_tmpdata_8)

                   deallocate(data_h_9km_grid_8)
                   deallocate(data_v_9km_grid_8)

                   deallocate(data_h_9km_tile_8)
                   deallocate(data_v_9km_tile_8)

                   allocate(Obs_f_tmpdata(N_obsf))

                   allocate(data_h_9km_grid(tile_grid_g%N_lon,tile_grid_g%N_lat))
                   allocate(data_v_9km_grid(tile_grid_g%N_lon,tile_grid_g%N_lat))

                   allocate(data_h_9km_tile(N_catf))
                   allocate(data_v_9km_tile(N_catf))

                end if

                ! "orig_obs"

                is_obsassim = .false.

                select case (k)

                case(1); Obs_f_tmpdata_8 = Observations_f(1:N_obsf)%time; is_orbit = .false.
                case(2); Obs_f_tmpdata   = Obs_f_resflag;                 is_orbit = .false.
                case(3); Obs_f_tmpdata   = Obs_f_orbflag;                 is_orbit = .true.
                case(4); Obs_f_tmpdata   = Observations_f(1:N_obsf)%obs;  is_orbit = .false.

                end select

             else  ! "obs_fcst"

                is_orbit = .false.

                select case (k)
                   
                case(1); Obs_f_tmpdata = Observations_f(1:N_obsf)%obs;     is_obsassim = .true.
                case(2); Obs_f_tmpdata = Observations_f(1:N_obsf)%obsvar;  is_obsassim = .true.  ! see sqrt below
                case(3); Obs_f_tmpdata = Observations_f(1:N_obsf)%fcst;    is_obsassim = .false.
                case(4); Obs_f_tmpdata = Observations_f(1:N_obsf)%fcstvar; is_obsassim = .false. ! see sqrt below
                   
                end select

                ! convert *var into *std (but preserve no-data-values)

                if (k==2 .or. k==4) then

                   do m=1,N_obsf

                      nodatavalue = obs_param(Observations_f(m)%species)%nodata

                      tol = abs(nodatavalue*nodata_tolfrac_generic)

                      if (abs(Obs_f_tmpdata(m)-nodatavalue)>tol) &
                           Obs_f_tmpdata(m) = sqrt(Obs_f_tmpdata(m))

                   end do

                end if

             end if

             ! map "Obs_f_tmpdata" onto global 9 km EASE grid

             ndata_h_9km_grid = 0
             ndata_v_9km_grid = 0

             if (use_real8) then

                data_h_9km_grid_8 = 0.0D0
                data_v_9km_grid_8 = 0.0D0

             else

                data_h_9km_grid   = 0.
                data_v_9km_grid   = 0.

             end if

             ! orbit flag requires special case

             if (is_orbit) then

                data_h_9km_grid   = nodata_generic
                data_v_9km_grid   = nodata_generic

                tol = abs(nodata_generic*nodata_tolfrac_generic)

             end if

             do n=1,N_obsf

                ! write "obs_assim" (and obsvar) only if assimilated

                if ( (is_obsassim .and. Observations_f(n)%assim) .or. (.not. is_obsassim) ) then
                   
                   this_species = Observations_f(n)%species

                   select case (trim(obs_param(this_species)%descr))

                   case('SMAP_L1C_Tbh_D', 'SMAP_L1C_Tbh_A',                             &
                        'SMAP_L1C_Tbh_E09_D', 'SMAP_L1C_Tbh_E09_A',                     &
                        'SMAP_L1C_Tbh_E27_D', 'SMAP_L1C_Tbh_E27_A',                     &
                        'SMAP_L2AP_Tbh_D',    'SMAP_L2AP_Tbh_A',                        &
                        'SMOS_fit_Tbh_D', 'SMOS_fit_Tbh_A'                     )

                      ! H-pol species

                      if (is_orbit) then   ! orbit flag

                         ! determine existing orbit flag (scalar)

                         existing_orbflag = data_h_9km_grid(col_beg_9km(n),row_beg_9km(n))

                         ! if the obs covers more than one 9 km grid cell,
                         !  make sure (existing) orbit flags are the same for all
                         !  9 km grid cells in question

                         if ( (col_beg_9km(n) < col_end_9km(n)) .or.                    &
                              (row_beg_9km(n) < row_end_9km(n))      ) then

                            if (any(                                                    &
                                 abs(data_h_9km_grid(                                   &
                                 col_beg_9km(n):col_end_9km(n),                         &
                                 row_beg_9km(n):row_end_9km(n) )                        &
                                 - existing_orbflag) > 0.01       )) then
                               err_msg = 'orbit flags differ for 9 km grid cells (H-pol)'
                               call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                            end if

                         end if

                         ! set orbit flag of 9 km grid cells in question

                         if     (abs(existing_orbflag - nodata_generic)  <tol) then

                            ! orbit flag was nodata, set for first time

                            data_h_9km_grid(                        &
                                 col_beg_9km(n):col_end_9km(n),     &
                                 row_beg_9km(n):row_end_9km(n)      &
                                 )                                  &
                                 = Obs_f_tmpdata(n)

                         elseif (abs(existing_orbflag - Obs_f_tmpdata(n))>tol) then

                            ! existing orbit flag differs from current, set to 0
                            ! (indicating that obs in aup output were averaged
                            !  across different orbit directions)

                            data_h_9km_grid(                        &
                                 col_beg_9km(n):col_end_9km(n),     &
                                 row_beg_9km(n):row_end_9km(n)      &
                                 )                                  &
                                 = 0.

                         end if

                      else  ! all fields *except* orbit flag

                         ! ndata_grid(ind) = ndata_grid(ind) + 1

                         ndata_h_9km_grid(                          &
                              col_beg_9km(n):col_end_9km(n),        &
                              row_beg_9km(n):row_end_9km(n)         &
                              )                                     &
                              =                                     &
                              ndata_h_9km_grid(                     &
                              col_beg_9km(n):col_end_9km(n),        &
                              row_beg_9km(n):row_end_9km(n)         &
                              )                                     &
                              + 1

                         ! data_grid(ind)  = data_grid(ind)  + Obs_f_tmpdata(n)

                         if (use_real8) then

                            data_h_9km_grid_8(                         &
                                 col_beg_9km(n):col_end_9km(n),        &
                                 row_beg_9km(n):row_end_9km(n)         &
                                 )                                     &
                                 =                                     &
                                 data_h_9km_grid_8(                    &
                                 col_beg_9km(n):col_end_9km(n),        &
                                 row_beg_9km(n):row_end_9km(n)         &
                                 )                                     &
                                 + Obs_f_tmpdata_8(n)

                         else

                            data_h_9km_grid(                           &
                                 col_beg_9km(n):col_end_9km(n),        &
                                 row_beg_9km(n):row_end_9km(n)         &
                                 )                                     &
                                 =                                     &
                                 data_h_9km_grid(                      &
                                 col_beg_9km(n):col_end_9km(n),        &
                                 row_beg_9km(n):row_end_9km(n)         &
                                 )                                     &
                                 + Obs_f_tmpdata(n)

                         end if

                      end if

                   case('SMAP_L1C_Tbv_D', 'SMAP_L1C_Tbv_A',                             &
                        'SMAP_L1C_Tbv_E09_D', 'SMAP_L1C_Tbv_E09_A',                     &
                        'SMAP_L1C_Tbv_E27_D', 'SMAP_L1C_Tbv_E27_A',                     &
                        'SMAP_L2AP_Tbv_D',    'SMAP_L2AP_Tbv_A',                        &
                        'SMOS_fit_Tbv_D', 'SMOS_fit_Tbv_A'                     )

                      ! V-pol species

                      if (is_orbit) then   ! orbit flag

                         ! determine existing orbit flag (scalar)

                         existing_orbflag = data_v_9km_grid(col_beg_9km(n),row_beg_9km(n))

                         ! if the obs covers more than one 9 km grid cell,
                         !  make sure (existing) orbit flags are the same for all
                         !  9 km grid cells in question

                         if ( (col_beg_9km(n) < col_end_9km(n)) .or.                    &
                              (row_beg_9km(n) < row_end_9km(n))      ) then

                            if (any(                                                    &
                                 abs(data_v_9km_grid(                                   &
                                 col_beg_9km(n):col_end_9km(n),                         &
                                 row_beg_9km(n):row_end_9km(n) )                        &
                                 - existing_orbflag) > 0.01       )) then
                               err_msg = 'orbit flags differ for 9 km grid cells (V-pol)'
                               call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                            end if
                         end if

                         ! set orbit flag of 9 km grid cells in question

                         if     (abs(existing_orbflag - nodata_generic)  <tol) then

                            ! orbit flag was nodata, set for first time

                            data_v_9km_grid(                        &
                                 col_beg_9km(n):col_end_9km(n),     &
                                 row_beg_9km(n):row_end_9km(n)      &
                                 )                                  &
                                 = Obs_f_tmpdata(n)

                         elseif (abs(existing_orbflag - Obs_f_tmpdata(n))>tol) then

                            ! existing orbit flag differs from current, set to 0
                            ! (indicating that obs in aup output were averaged
                            !  across different orbit directions)

                            data_v_9km_grid(                        &
                                 col_beg_9km(n):col_end_9km(n),     &
                                 row_beg_9km(n):row_end_9km(n)      &
                                 )                                  &
                                 = 0.

                         end if

                      else  ! all fields *except* orbit flag

                         ! ndata_grid(ind) = ndata_grid(ind) + 1

                         ndata_v_9km_grid(                          &
                              col_beg_9km(n):col_end_9km(n),        &
                              row_beg_9km(n):row_end_9km(n)         &
                              )                                     &
                              =                                     &
                              ndata_v_9km_grid(                     &
                              col_beg_9km(n):col_end_9km(n),        &
                              row_beg_9km(n):row_end_9km(n)         &
                              )                                     &
                              + 1

                         ! data_grid(ind)  = data_grid(ind)  + Obs_f_tmpdata(n)

                         if (use_real8) then

                            data_v_9km_grid_8(                         &
                                 col_beg_9km(n):col_end_9km(n),        &
                                 row_beg_9km(n):row_end_9km(n)         &
                                 )                                     &
                                 =                                     &
                                 data_v_9km_grid_8(                    &
                                 col_beg_9km(n):col_end_9km(n),        &
                                 row_beg_9km(n):row_end_9km(n)         &
                                 )                                     &
                                 + Obs_f_tmpdata_8(n)

                         else

                            data_v_9km_grid(                           &
                                 col_beg_9km(n):col_end_9km(n),        &
                                 row_beg_9km(n):row_end_9km(n)         &
                                 )                                     &
                                 =                                     &
                                 data_v_9km_grid(                      &
                                 col_beg_9km(n):col_end_9km(n),        &
                                 row_beg_9km(n):row_end_9km(n)         &
                                 )                                     &
                                 + Obs_f_tmpdata(n)

                         end if

                      end if

                   case default

                      ! do nothing

                   end select
                   
                end if  ! write "obs_assim" (and obsvar) only if assimilated
                
             end do     ! n=1,N_obsf

             ! normalize and set to nodata as needed

             if (.not. is_orbit) then

                ! H-pol species

                if (use_real8) then

                   where     (ndata_h_9km_grid>1)

                      data_h_9km_grid_8 = data_h_9km_grid_8/real(ndata_h_9km_grid,kind(0.0D0))

                   elsewhere (ndata_h_9km_grid==0)

                      data_h_9km_grid_8 = real(nodata_generic,kind(0.0D0))

                   end where

                else

                   where     (ndata_h_9km_grid>1)

                      data_h_9km_grid   = data_h_9km_grid  /real(ndata_h_9km_grid)

                   elsewhere (ndata_h_9km_grid==0)

                      data_h_9km_grid   =      nodata_generic

                   end where

                end if

                ! V-pol species

                if (use_real8) then

                   where     (ndata_v_9km_grid>1)

                      data_v_9km_grid_8 = data_v_9km_grid_8/real(ndata_v_9km_grid,kind(0.0D0))

                   elsewhere (ndata_v_9km_grid==0)

                      data_v_9km_grid_8 = real(nodata_generic,kind(0.0D0))

                   end where

                else

                   where     (ndata_v_9km_grid>1)

                      data_v_9km_grid   = data_v_9km_grid  /real(ndata_v_9km_grid)

                   elsewhere (ndata_v_9km_grid==0)

                      data_v_9km_grid   =      nodata_generic

                   end where

                end if

             end if  ! (.not. is_orbit)

             ! done mapping to 9 km grid

             ! -----------------------------------
             !
             ! map to tile space

             if (use_real8) then

                data_h_9km_tile_8 = real(nodata_generic,kind(0.0D0)) ! init (not in grid2tile!)
                data_v_9km_tile_8 = real(nodata_generic,kind(0.0D0)) ! init (not in grid2tile!)

                call grid2tile( tile_grid_g, N_catf, tile_coord_f%i_indg, tile_coord_f%j_indg, data_h_9km_grid_8, &
                     data_h_9km_tile_8 )

                call grid2tile( tile_grid_g, N_catf, tile_coord_f%i_indg, tile_coord_f%j_indg, data_v_9km_grid_8, &
                     data_v_9km_tile_8 )

                ! write into file

                write (unitnum) (     data_h_9km_tile_8(n) , n=1,N_catf)
                write (unitnum) (     data_v_9km_tile_8(n) , n=1,N_catf)

             else

                data_h_9km_tile = nodata_generic  ! initialize (not done in grid2tile!)
                data_v_9km_tile = nodata_generic  ! initialize (not done in grid2tile!)

                call grid2tile( tile_grid_g, N_catf, tile_coord_f%i_indg,tile_coord_f%j_indg, data_h_9km_grid, &
                     data_h_9km_tile )

                call grid2tile( tile_grid_g, N_catf, tile_coord_f%i_indg,tile_coord_f%j_indg, data_v_9km_grid, &
                     data_v_9km_tile )

                ! write into file

                if (option=='orig_obs' .and. (k==2 .or. k==3)) then

                   ! convert to integer before writing to file

                   write (unitnum) (nint(data_h_9km_tile(n)), n=1,N_catf)
                   write (unitnum) (nint(data_v_9km_tile(n)), n=1,N_catf)

                else

                   write (unitnum) (     data_h_9km_tile(n) , n=1,N_catf)
                   write (unitnum) (     data_v_9km_tile(n) , n=1,N_catf)

                end if

             end if

          end do ! loop through k=1,kmax

          ! clean up


          deallocate(col_beg_9km)
          deallocate(col_end_9km)
          deallocate(row_beg_9km)
          deallocate(row_end_9km)

          deallocate(Obs_f_resflag)
          deallocate(Obs_f_orbflag)

          deallocate(ndata_h_9km_grid)
          deallocate(ndata_v_9km_grid)

          deallocate(Obs_f_tmpdata)

          deallocate(data_h_9km_grid)
          deallocate(data_v_9km_grid)

          deallocate(data_h_9km_tile)
          deallocate(data_v_9km_tile)

       end if  ! root_proc

       if(allocated(Observations_f)) deallocate(Observations_f)

    end if     ! (option=='orig_obs' .or. option=='obs_fcst')

    ! --------------------------------------------------------------
    !
    ! assemble state-space data for writing as needed

    if (option=='obs_fcst' .or. option=='analysis') then

       ! diagnose variables of interest for the ensemble

       do n_e=1,N_ens

          srfexc = cat_progn(:,n_e)%srfexc
          rzexc  = cat_progn(:,n_e)%rzexc
          catdef = cat_progn(:,n_e)%catdef

          call catch_calc_soil_moist(                                             &
               N_catl,           cat_param%dzsf,   cat_param%vgwmax,              &
               cat_param%cdcr1,  cat_param%cdcr2,  cat_param%psis,                &
               cat_param%bee,    cat_param%poros,  cat_param%wpwet,               &
               cat_param%ars1,   cat_param%ars2,   cat_param%ars3,                &
               cat_param%ara1,   cat_param%ara2,   cat_param%ara3,                &
               cat_param%ara4,   cat_param%arw1,   cat_param%arw2,                &
               cat_param%arw3,   cat_param%arw4,                                  &
               cat_param%bf1,    cat_param%bf2,                                   &
               srfexc,           rzexc,            catdef,                        &
               ar1,  ar2,  ar4, sfmc(:,n_e), rzmc(:,n_e), prmc(:,n_e) )

          call catch_calc_tsurf( N_catl,                                          &
               cat_progn(:,n_e)%tc1, cat_progn(:,n_e)%tc2, cat_progn(:,n_e)%tc4,  &
               catprogn2wesn(N_catl,cat_progn(:,n_e)),                            &
               catprogn2htsn(N_catl,cat_progn(:,n_e)),                            &
               ar1,  ar2,  ar4,                                                   &
               tsurf(:,n_e) )

          ! NOTE: "tp" is returned in CELSIUS [for consistency w/ catchment.F90]

          call catch_calc_tp( N_catl, cat_param%poros,                            &
               catprogn2ghtcnt(N_catl,cat_progn(:,n_e)), tp )

          tp1(:,n_e) = tp(1,:) + MAPL_TICE

       end do

       ! compute ensemble mean values and std-dev

       call row_std( N_catl, N_ens, sfmc,  tile_std_l(:,1), tile_mean_l(:,1) )
       call row_std( N_catl, N_ens, rzmc,  tile_std_l(:,2), tile_mean_l(:,2) )
       call row_std( N_catl, N_ens, prmc,  tile_std_l(:,3), tile_mean_l(:,3) )
       call row_std( N_catl, N_ens, tsurf, tile_std_l(:,4), tile_mean_l(:,4) )
       call row_std( N_catl, N_ens, tp1,   tile_std_l(:,5), tile_mean_l(:,5) )

       ! make sure mean *mc values are between 0. and porosity

       do j=1,N_catl

          tile_mean_l(j,1) = max( min( tile_mean_l(j,1), cat_param(j)%poros), 0.)
          tile_mean_l(j,2) = max( min( tile_mean_l(j,2), cat_param(j)%poros), 0.)
          tile_mean_l(j,3) = max( min( tile_mean_l(j,3), cat_param(j)%poros), 0.)

       end do

       ! make sure std-dev values are non-negative

       tile_std_l = max( tile_std_l, 0.)

       ! write out (append) ensemble mean values
       !
       ! 1: sm_surface_[forecast/analysis]           [m3 m-3]
       ! 2: sm_rootzone_[forecast/analysis]          [m3 m-3]
       ! 3: sm_profile_[forecast/analysis]           [m3 m-3]
       ! 4: surface_temp_[forecast/analysis]         [K]
       ! 5: soil_temp_layer1_[forecast/analysis]     [K]

       allocate(tile_data_f(N_catf))

       do k=1,5   ! write output one field at a time

          ! gatherv tile data from local to full domain and map to grid

          call l2f_real( N_catl, N_catf, N_catl_vec, low_ind,  &
               tile_mean_l(:,k), tile_data_f)

          if (root_proc) write(unitnum) (tile_data_f(n), n=1,N_catf)

       end do

       if (option=='analysis') then

          ! write out (append) ensstd fields for analysis

          ! 1: sm_surface_analysis_ensstd           [m3 m-3]
          ! 2: sm_rootzone_analysis_ensstd          [m3 m-3]
          ! 3: sm_profile_analysis_ensstd           [m3 m-3]
          ! 4: surface_temp_analysis_ensstd         [K]
          ! 5: soil_temp_layer1_analysis_ensstd     [K]

          do k=1,5   ! write output one field at a time

             ! gatherv tile data from local to full domain and map to grid

             call l2f_real( N_catl, N_catf, N_catl_vec, low_ind,  &
                  tile_std_l(:,k), tile_data_f)

             if (root_proc) write(unitnum) (tile_data_f(n), n=1,N_catf)

          end do

       end if  ! (option=='analysis')

       deallocate(tile_data_f)

    end if     ! (option=='obs_fcst' .or. option=='analysis')

    ! --------------------------------------------------------------
    !
    ! close output file

    if (root_proc)  close(unitnum,status='keep')

  end subroutine write_smapL4SMaup

  ! -----------------------------------------------------------------

end module clsm_ensupd_enkf_update

! ============ EOF =========================================================
