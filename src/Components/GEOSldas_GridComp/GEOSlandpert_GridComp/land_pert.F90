!
! modules to generate land surface perturbations:
!  
!   module land_pert_types
!   module land_pert_routines
!
! can be used to perturb forcing fields such as precipitation and radiation
! or model prognostic variables such as soil moisture or soil temperature
!
! MUST initialize random seed and Pert_ntrmdt by calling 
!  get_pert() with initialize=.true. at the start of the driver program
!  (otherwise set initialize=.false.)
!
! compile line for test program:
! cpp -P -C nr_ran2_gasdev.f90 nr_ran2_gasdev.cpp.f90; cpp -P -C random_fields.f90 random_fields.cpp.f90; cpp -P -C land_pert.f90 land_pert.cpp.f90; f90 nr_ran2_gasdev.cpp.f90 random_fields.cpp.f90 land_pert.cpp.f90
!
! reichle, 24 Jan 2005
! reichle, 11 Feb 2005
! reichle, 26 May 2005
! reichle,  7 Jun 2005 - more init options
! reichle, 14 Apr 2006 - split land_pert.F90 into 2 files to avoid 
!                         having more than one module per file
! reichle,  8 Aug 2008 - added "logunit" for use within LDASsa
! reichle,  1 Oct 2009 - added "stop_it" for use within LDASsa
!
! ------------------------------------------------------------

module land_pert_routines

  use ESMF

  use LDAS_PertTypes,                   ONLY:     & 
       pert_param_type,                           &
       allocate_pert_param

  use nr_ran2_gasdev,                   ONLY:     &
       NRANDSEED,                                 &
       init_randseed

  use Random_FieldsMod
  use StringRandom_fieldsMapMod

  use nr_jacobi,                        ONLY:     &
       jacobi

  use LDAS_TileCoordType,               ONLY:     &
       grid_def_type

  use LDAS_ExceptionsMod,               ONLY:     & 
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use LDAS_ensdrv_Globals, only: root_logit, logunit

  use MAPL    
  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: get_pert
  public :: propagate_pert
  public :: assemble_forcepert_param
  public :: get_sqrt_corr_matrix
  public :: get_init_Pert_rseed
  public :: clear_rf

  ! **********************************************************************

  type(StringRandom_fieldsMap) :: random_fieldsMap

contains
  
  ! **********************************************************************
  
  subroutine get_pert(                                &
       N_pert, fft_npert,  N_ens,                     &
       pert_grid_l, pert_grid_f,                      &
       dtstep,                                        &
       pert_param,                                    &
       Pert_rseed,                                    &
       Pert_ntrmdt,                                   &
       Pert,                                          &
       initialize_rseed,                              &
       initialize_ntrmdt,                             &
       diagnose_pert_only       )

    ! get perturbations
    !
    ! reichle, 22 Jun 2005
    ! reichle, 17 Jul 2020 - switched order of input arguments pert_grid_f and pert_grid_l
    !                         for consistency with other subroutines
    !
    ! This subroutine unifies subroutines GEOSldas_get_pert() and LDASsa_get_pert().
    ! GEOSldas_get_pert() was used for forcing and prognostics perturbations and did
    ! not have an array dimension for ensemble members.
    ! LDASsa_get_pert() was used for observations perturbations and included an array
    ! dimension for ensemble members (and a call to subroutine adjust_mean()).
    ! Otherwise the two subroutines were identical.
    ! - wjiang+reichle, 9 Apr 2021

    
    implicit none

    ! N_pert is the number of *perturbation* fields and is not 
    ! necessarily equal to the number of forcing fields.
    ! E.g. generate N_pert=3 perturbations for precip, 
    ! shortwave radiation, and longwave radiation. These can then be
    ! applied to N_force forcing fields, possibly various precip fields 
    ! (incl large-scale & convective precip and snow) and to radiation 
    ! fields.

    integer, intent(in) :: N_pert     ! # different perturbations
    integer, intent(in) :: fft_npert  ! # different perturbations for fft; equals n_pert if proc has tiles.

    integer, intent(in) :: N_ens  ! # ensemble members

    type(grid_def_type), intent(in) :: pert_grid_l
    type(grid_def_type), intent(in) :: pert_grid_f

    real, intent(in) :: dtstep        ! perturbation time step in seconds

    ! Parameter structure for perturbations (see type definition for details).

    type(pert_param_type), dimension(:), pointer :: pert_param

    ! Pert_ntrmdt are intermediate perturbation fields
    ! that need to be remembered between calls to this subroutine.
    ! In essence, they store N_pert mutually uncorrelated
    ! perturbation fields of standard-normal distribution.
    ! Pert_rseed is the random seed for the generation of
    ! Pert_ntrmdt and is treated similarly to a prognostic variable.
    ! Each ensemble member has its own random seed.

    integer, dimension(NRANDSEED,N_ens), intent(inout) :: Pert_rseed

    real, dimension(pert_grid_l%N_lon, pert_grid_l%N_lat, N_pert, N_ens), intent(inout) :: Pert_ntrmdt

    ! Pert are N_pert cross-correlated perturbation
    ! fields that are rotated and scaled versions of Pert_ntrmdt
    ! so that Pert has the the mean values, standard deviations and
    ! cross-correlations specified in pert_param.  
    ! The distribution is lognormal for multiplicative perturbations.  
    ! Pert should be used as follows for field F (eg. large-scale
    ! precip, convective precip, lw radiation, ...)
    !
    ! Fprime = F+Pert   for additive perturbations
    ! Fprime = F*Pert   for multiplicative perturbations
    !
    ! Note that this subroutine does NOT ensure physically meaningful
    ! perturbed fields.  This is best done outside this subroutine
    ! after the perturbations have been applied.    

    real, dimension(pert_grid_l%N_lon, pert_grid_l%N_lat, N_pert, N_ens), intent(  out) :: Pert

    ! If initialize_rseed==.true., set initial random seed vector.
    ! If initialize_rseed==.true., the first row of Pert_rseed must be
    ! filled with a different negative integer for each ensemble member.
    ! See sample subroutine get_init_Pert_rseed().
    !
    ! If initialize_ntrmdt==.true., generate initial Pert_ntrmdt (must be 
    ! allocated!!!).  
    !
    ! Note that when get_pert() is used for generating independent 
    ! perturbations to forcings and prognostic variables, only one
    ! common Pert_rseed should be used, so one of the "ntrmdt" fields
    ! must be initialized without initializing "rseed" again.

    ! If initialize_*==.false. or absent, Pert_rseed and 
    ! Pert_ntrmdt must be what was obtained as output from the last call 
    ! to get_pert().  There is only one exception: if there are no temporal 
    ! correlations, it is not necessary to remember Pert_ntrmdt.

    ! If diagnose_pert_only==.true., Pert_ntrmdt must be available (for
    ! example from a restart file) and Pert will then be "diagnosed" from
    ! Pert_ntrmdt (and initialize_* must be .false.).

    logical, intent(in), optional :: initialize_rseed
    logical, intent(in), optional :: initialize_ntrmdt

    logical, intent(in), optional :: diagnose_pert_only

    ! --------------------------------------------------
    !
    ! local variables

    integer :: i, j, m, mm, n

    real, dimension(pert_grid_l%N_lon, pert_grid_l%N_lat)  :: tmp_grid

    real :: tmpreal

    logical :: init_rseed, init_ntrmdt, diagn_only

    character(len=*), parameter :: Iam = 'get_pert'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------------
    !
    ! initialize random seed if necessary

    init_rseed  = .false.
    init_ntrmdt = .false.

    if (present(initialize_rseed))   init_rseed  = initialize_rseed
    if (present(initialize_ntrmdt))  init_ntrmdt = initialize_ntrmdt

    if (init_rseed) then

       do n=1,N_ens

          call init_randseed(Pert_rseed(:,n))

       end do

    end if

    ! ------------------------------------------------------------------

    diagn_only = .false.

    if (present(diagnose_pert_only))  diagn_only = diagnose_pert_only

    if ( diagn_only .and. (init_rseed .or. init_ntrmdt) ) then
       err_msg = 'contradictory optional inputs'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! ------------------------------------------------------------------
    !
    ! Pert_ntrmdt are standard-normal with desired 
    ! temporal and spatial correlation structure.
    ! Cross-correlations between different fields and scaling to desired
    ! mean and variance is NOT included in Pert_ntrmdt.
    !
    ! Propagate perturbation fields:
    ! On input, Pert_ntrmdt must contain fields from last time step.
    ! (If init_ntrmdt=.true. Pert_ntrmdt is initialized to a 
    ! standard-normal field with the desired spatial correlation structure)

    if (.not. diagn_only)                            &
         call propagate_pert(                        &
         fft_npert, N_ens, pert_grid_l, pert_grid_f, &
         dtstep,                                     &
         Pert_rseed,                                 &
         pert_param,                                 &
         Pert_ntrmdt,                                &
         init_ntrmdt                           )

    ! compute diagnostic "Pert" 

    ! ensure that ensemble mean perturbation is zero
    ! (must have N_ens>2 for this to make sense).
    !
    ! NOTE: since the sample mean model error varies spatially,
    !       this adjustment slightly changes the *spatial* mean and 
    !       covariance of the model error fields
    !       likely, the benefits of the adjustments for small
    !       ensemble sizes outweigh the disadvantages of altering
    !       the statistical properties

    do m=1,N_pert

       if ( (pert_param(m)%zeromean) .and. (N_ens>2)) then

          do i=1,pert_grid_l%N_lon

             call adjust_mean(pert_grid_l%N_lat, N_ens, Pert_ntrmdt(i,:,m,:) )

          end do

       end if

    end do

    ! compute rotated fields to get desired cross-correlations between
    ! different fields, then scale to desired mean and std

    do m=1,N_pert

       do n=1,N_ens

          ! rotate to get desired multivariate correlations

          do j=1,pert_grid_l%N_lat
             do i=1,pert_grid_l%N_lon

                tmp_grid(i,j) = 0.

                do mm=1,N_pert

                   tmp_grid(i,j) = tmp_grid(i,j) + &
                        pert_param(m)%ccorr(mm,i,j) * Pert_ntrmdt(i,j,mm,n)

                end do

             end do
          end do

          ! scale back freak outliers

          call truncate_std_normal( pert_grid_l%N_lon, pert_grid_l%N_lat, &
               pert_param(m)%std_normal_max, tmp_grid )

          ! scale

          do j=1,pert_grid_l%N_lat
             do i=1,pert_grid_l%N_lon

                tmpreal = pert_param(m)%mean(i,j) + &
                     pert_param(m)%std(i,j) * tmp_grid(i,j)

                select case (pert_param(m)%typ)

                case (0)        ! additive

                   Pert(i,j,m,n) = tmpreal

                case (1)        ! multiplicative and lognormal

                   Pert(i,j,m,n) = exp(tmpreal)

                case default

                   err_msg = 'encountered unknown typ_pert'
                   call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

                end select

             end do
          end do

       end do  ! end loop through ensemble members (1:N_ens)
    end do     ! end loop through different fields (1:N_pert)

  end subroutine get_pert

  ! *******************************************************************************
  
  subroutine get_init_Pert_rseed( ens_id, init_Pert_rseed )
    
    ! get initial random seed "init_Pert_rseed" for initializing 
    ! Pert_rseed within get_pert()
    !
    ! A different random seed is necessary for each ensemble member.
    !
    ! This subroutine is meant as a sample for how the initial Pert_rseed 
    ! can be set.  
    ! In this example, ens_id is meant to be a nonnegative small integer.
    
    implicit none
    
    integer, intent(in)  :: ens_id
    
    integer, intent(out) :: init_Pert_rseed
    
    ! --------------------------------------------------
    ! 
    ! local parameter values for initial random seed
    
    integer, parameter :: RSEED_CONST0 =     -1
    integer, parameter :: RSEED_CONST2 =    -10    ! must be negative
    
    ! --------------------------------------------------
    !
    ! local variables
    
    character(len=*), parameter :: Iam = 'get_init_Pert_rseed'
    character(len=400) :: err_msg
    
    ! -------------------------------------------------------------
    
    init_Pert_rseed = RSEED_CONST0 + ens_id*RSEED_CONST2

    ! make sure init_Pert_rseed is negative and no two numbers are the same
    if (init_Pert_rseed>=0) then
       
       err_msg = 'found nonnegative component of init_Pert_rseed'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
       
  end subroutine get_init_Pert_rseed
  
  ! *******************************************************************************
  
  subroutine propagate_pert(                            &
       N_pert, N_ens, pert_grid_l, pert_grid_f, dtstep, &
       Pert_rseed,                                      &
       pert_param,                                      &
       Pert_ntrmdt,                                     &
       initialize              )

    ! generate zero-mean, unit-variance (!!) time series 
    !  of N_pert 2d perturbation fields
    !
    ! can also be used just to get a set of 2d random fields (set dtstep 
    !  to arbitrary number and "initialize" to .true.)
    !
    ! on input, Pert_ntrmdt must contain the corresponding 
    !  perturbations from the previous time step 
    !
    ! accounts for temporal correlation with AR(1) approach
    !  (if pert_param%tcorr==0 then error is white in time) 
    !
    ! adapted from off-line EnKF, subroutine propagate_err() 
    !  from NCAT_59124_tskin in enkf_catchment.f90
    !
    ! reichle, 14 Feb 2005

    ! This subroutine unifies subroutines GEOSldas_propagate_pert() and LDASsa_propagate_pert().
    ! GEOSldas_propagate_pert() was used for forcing and prognostics perturbations and did
    ! not have an array dimension for ensemble members.
    ! LDASsa_propagate_pert() was used for observations perturbations and included an array
    ! dimension for ensemble members and a call to subroutine adjust_mean().
    ! Otherwise the two subroutines were identical.
    ! - wjiang+reichle, 9 Apr 2021

    implicit none

    ! ---------------------------    

    integer, intent(in) :: N_ens, N_pert

    type(grid_def_type), intent(in) :: pert_grid_l, pert_grid_f   ! local/full grids

    type(pert_param_type), dimension(N_pert), intent(in) :: pert_param

    real, intent(in) :: dtstep  ! time step of generation of error fields [s]

    integer, dimension(NRANDSEED,N_ens), intent(inout) :: Pert_rseed

    real, dimension(pert_grid_l%N_lon, pert_grid_l%N_lat, N_pert, N_ens), &
         intent(inout) :: Pert_ntrmdt

    logical, intent(in) :: initialize   ! switch

    ! ---------------------------    

    ! locals

    ! Depending on the input logical switch pert_param%coarsen, perturbations may 
    ! be computed on a coarsened grid with spacing automatically determined by 
    !  - the original perturbation grid spacing, 
    !  - the spatial correlation scale, and 
    !  - the following parameter:

    real, parameter :: coarsen_param = 0.8

    ! Example: For the global SMAP EASEv2 M09 grid with dlon~dlat~0.1 and
    !          xcorr=ycorr=0.5, coarsen_param=0.8 results in xstride=ystride=4, 
    !          which implies that perturbations are effectively computed on 
    !          the EASEv2 M36 grid.
    !
    ! ----------------------

    integer :: i, j, m, n, Nx, Ny, Nx_fft, Ny_fft, xStride, yStride, imax, jmax

    real    :: cc, dd, xCorr, yCorr, tCorr, tmpReal, rdlon, rdlat

    real, dimension(pert_grid_f%N_lon,pert_grid_f%N_lat), target :: rfield, rfield2

    real, dimension(:,:), pointer :: ptr2rfield, ptr2rfield2

    logical :: white_in_time, white_in_space, stored_field

    integer :: tmpInt, xstart, xend, ystart, yend

    type(random_fields), pointer :: rf

    type(ESMF_VM) :: vm
    integer :: mpicomm, status  

    call ESMF_VmGetCurrent(vm, rc=status)
    call ESMF_VmGet(vm, mpicommunicator=mpicomm, rc=status)

    do m=1,N_pert

       ! shorthand

       xCorr = pert_param(m)%xcorr
       yCorr = pert_param(m)%ycorr
       tCorr = pert_param(m)%tcorr

       ! get parameters for temporal correlation

       if ((.not. initialize) .and. (tCorr>0.0)) then

          white_in_time = .false.
          cc = exp( - dtstep / tCorr )
          dd = sqrt( 1 - cc**2 )
       else
          cc = 0.
          dd = 1.
          white_in_time = .true.
       end if

       ! find out whether there are spatial correlations
       if ( (xCorr>0.0) .or. (yCorr>0.0) ) then
          white_in_space = .false.
       else
          white_in_space = .true.
       end if


       ! get grid parameters for generation of new random fields on 
       ! possibly coarsened grid

       call calc_fft_grid(pert_param(m), pert_grid_f, Nx, Ny, Nx_fft, Ny_fft, xStride, yStride, rdlon, rdlat)

       ptr2rfield  => rfield( 1:pert_grid_f%N_lon:xStride,1:pert_grid_f%N_lat:yStride)
       ptr2rfield2 => rfield2(1:pert_grid_f%N_lon:xStride,1:pert_grid_f%N_lat:yStride)

       ! generate new random fields and propagate AR(1)
       !
       ! Note that rfg2d always produces a pair of random fields!
       !
       ! Use logical variable "stored_field" to figure out whether a second
       ! standard-normal random field is available for next ensemble member.
       ! (in other words, this subroutine is most efficient if N_ens=even 
       ! number, and it is least efficient if N_ens=1)

       stored_field = .false.

       ! initialize instance rf of class random_fields
       ! this needs to be done for each pert field
#ifdef MKL_AVAILABLE      
       ! W.J Note: hardcoded comm = mpicomm to activate parallel fft
       rf => find_rf(Nx, Ny, Nx_fft, Ny_fft, comm=mpicomm )
#else
       rf => find_rf(Nx, Ny, Nx_fft, Ny_fft)
#endif

       do n=1,N_ens

          ! generate a random field

          if (white_in_space) then

             call rf%generate_white_field(Pert_rseed(:,n), ptr2rfield)

          else       ! spatially correlated random fields

             ! NOTE: rfg2d_fft() relies on CXML math library (22 Feb 05)
             ! rfg2d_fft() now relies on Intel MKL (19 Jun 13)

             if (.not. stored_field) then
                call rf%rfg2d_fft(Pert_rseed(:,n), ptr2rfield, ptr2rfield2, xCorr, yCorr, rdlon, rdlat)
                stored_field = .true.
             else
                rfield = rfield2
                stored_field = .false.
             end if

          end if

          !! -----------------------------------------------------------
          !!
          !! adjust std of fields to match exactly 1.0
          !!
          !! WARNING: before doing this should check that 
          !!          N_x*dx>>xcorr .and. N_y*dy>>ycorr .and. N_x*N_y>>1 
          !!
          !! Cannot use adjust_std if the field is small relative to 
          !! its spatial correlation scales or if it contains only a few
          !! grid cells (even for white noise in space)
          !!
          !! reichle, 25 Jan 2005
          !!
          !! call adjust_std( loc_grid%N_x, loc_grid%N_y, rfield )
          !!
          !! -----------------------------------------------------------

          ! copy to fine grid

          ! [At a later time, perhaps insert bilinear interpolation here. 
          !  For that, will need rfield grid to extend *beyond* pert_grid_f!]

          do i=1,pert_grid_f%N_lon,xStride

             do j=1,pert_grid_f%N_lat,yStride

                tmpReal = rfield(i,j)

                imax = min( i+xStride-1, pert_grid_f%N_lon )
                jmax = min( j+yStride-1, pert_grid_f%N_lat )

                rfield(i:imax,j:jmax) = tmpReal

             end do
          end do

          ! restrict rfield to local pert grid
          if (pert_grid_l%n_lon /=0) then
             tmpInt = pert_grid_l%i_offg - pert_grid_f%i_offg
             xstart = tmpInt + 1
             xend   = tmpInt + pert_grid_l%N_lon

             tmpInt = pert_grid_l%j_offg - pert_grid_f%j_offg
             ystart = tmpInt + 1
             yend   = tmpInt + pert_grid_l%N_lat

             ! propagate AR(1) 

             if (white_in_time) then
                Pert_ntrmdt(:,:,m,n) = rfield(xstart:xend, ystart:yend)
             else
                Pert_ntrmdt(:,:,m,n) = cc*Pert_ntrmdt(:,:,m,n) + dd*rfield(xstart:xend, ystart:yend)
             end if
           endif
       end do ! n=1,N_ens

       ! finalize rf
       ! The rf map will be destroy in the finalize of GEOSLandperp_Gridcomp
       !call rf%finalize

    end do ! m=1,N_pert

  end subroutine propagate_pert

  subroutine calc_fft_grid(pert_param, pert_grid_f, Nx, Ny, N_x_fft, N_y_fft, xStride, yStride, rdlon, rdlat) 
    type(pert_param_type), intent(in) :: pert_param
    type(grid_def_type), intent(in)   :: pert_grid_f
    integer, intent(out) :: Nx, Ny, N_x_fft, N_y_fft, xStride, yStride
    real, intent(out)    :: rdlon, rdlat

    integer :: Nx_fft, Ny_fft
    real, parameter :: mult_of_xcorr = 2.
    real, parameter :: mult_of_ycorr = 2.
    real, parameter :: coarsen_param = 0.8
    real :: xCorr, yCorr

    xCorr = pert_param%xcorr
    yCorr = pert_param%ycorr

    xStride = 1
    yStride = 1
    if (pert_param%coarsen) then
       xStride = max( 1, floor(coarsen_param * xCorr / pert_grid_f%dlon) )
       yStride = max( 1, floor(coarsen_param * yCorr / pert_grid_f%dlat) )
    endif
    rdlon = real(xStride)*pert_grid_f%dlon
    rdlat = real(yStride)*pert_grid_f%dlat

    ! NOTE: number of grid cells of coarsened grid might not evenly divide 
    ! that of pert_grid_f

    Nx = pert_grid_f%N_lon / xStride
    Ny = pert_grid_f%N_lat / yStride

    if (mod(pert_grid_f%N_lon,xStride)>0) Nx = Nx + 1
    if (mod(pert_grid_f%N_lat,yStride)>0) Ny = Ny + 1
  
    ! add minimum required correlation lengths 
    Nx_fft = Nx + ceiling(mult_of_xcorr*xCorr/rdlon)
    Ny_fft = Ny + ceiling(mult_of_ycorr*yCorr/rdlat)
       
    ! ensure N_x_fft, N_y_fft are powers of two
    N_x_fft = 2**ceiling(log(real(Nx_fft))/log(2.))
    N_y_fft = 2**ceiling(log(real(Ny_fft))/log(2.))

  end subroutine
  ! ******************************************************************
  
  subroutine truncate_std_normal( N_x, N_y, std_normal_max, grid_data )
    
    ! truncate a realization of standard normal variables
    ! (scale back freak outliers)
    
    implicit none
    
    integer, intent(in) :: N_x, N_y
    
    real, intent(in) :: std_normal_max
    
    real, dimension(N_x,N_y), intent(inout) :: grid_data
    
    ! local variables
    
    integer :: i,j
    
    ! --------------------------------------------------------
    
    do i=1,N_x
       do j=1,N_y
          
          ! want: -std_normal_max < cat_data < std_normal_max
          
          grid_data(i,j) = &
               sign( min(abs(grid_data(i,j)),std_normal_max), grid_data(i,j) )
          
       end do
    end do
    
  end subroutine truncate_std_normal
  
  ! **************************************************************
  
  subroutine adjust_mean( N_row, N_col, A, M )
    
    ! adjust N_row by N_col matrix A such that 
    ! mean over columns for each row is given by the
    ! corresponding element in vector M of length N_row
    ! 
    ! vector of mean values M is optional input, if not present 
    ! zero mean is assumed
    
    implicit none
    
    integer, intent(in) :: N_row, N_col
    
    real, intent(inout), dimension(N_row,N_col)  :: A
    
    real, intent(in), optional, dimension(N_row) :: M
    
    ! ----------------------------
    
    ! locals
    
    integer i
    
    real, dimension(N_row) :: correction
    
    ! ------------------------------------------------------------
    
    if (present(M)) then
       correction = M - sum(A,2)/real(N_col) 
    else
       correction = - sum(A,2)/real(N_col) 
    end if
    
    do i=1,N_col
       A(:,i) = A(:,i) + correction
    end do
    
  end subroutine adjust_mean
  
  ! *************************************************************************
  
#if 0
  ! This subroutine is not used.
  ! Note that my_matrix_functions.F90 contains another (commented-out) version.
  ! wjiang + reichle, 25 Nov 2020
  
  subroutine adjust_std( N_row, N_col, A, std )
    
    ! adjust N_row by N_col matrix A such that (sample) standard deviation
    !  of all elements is exactly equal to std
    ! 
    ! std is optional input, if not present std=1 is assumed
    
    implicit none
    
    integer, intent(in) :: N_row, N_col
    
    real, intent(inout), dimension(N_row,N_col)  :: A
    
    real, intent(in), optional :: std
    
    ! ----------------------------
    
    ! locals
    
    integer :: i, j
    
    real :: correction, sample_std
    
    ! ------------------------------------------------------------
    
    ! compute sample std
    
    call matrix_std( N_row, N_col, A, sample_std )
    
    if (present(std)) then
       correction = std/sample_std
    else
       correction = 1./sample_std
    end if
    
    do i=1,N_row
       do j=1,N_col
          A(i,j) = correction*A(i,j)
       end do
    end do
    
  end subroutine adjust_std
#endif
  
  ! ***************************************************************************
  
#if 0
  ! This subroutine is not used.
  ! Note that my_matrix_functions.F90 contains other (commented-out) versions.
  ! wjiang + reichle, 25 Nov 2020
  
  subroutine matrix_std( N_row, N_col, A, std )
    
    ! compute std of all elements of N_row by N_col matrix A
    
    implicit none
    
    integer, intent(in) :: N_row, N_col
    
    real, intent(inout), dimension(N_row,N_col)  :: A
    
    real, intent(out) :: std
    
    ! ----------------------------
    
    ! locals
    
    integer :: i, j
    
    real :: x2, m, N_real, N_real_minus_one
    
    ! ------------------------------------------------------------
    
    N_real = real(N_row)*real(N_col)
    
    N_real_minus_one = N_real - 1.
    
    ! compute sample std
    
    x2 = 0.0
    m  = 0.0
    
    do i=1,N_row
       do j=1,N_col
          m  = m  + A(i,j)
          x2 = x2 + A(i,j)*A(i,j)
       end do
    end do
    
    std = sqrt( ( x2 - m**2/N_real )/N_real_minus_one )
    
  end subroutine matrix_std
#endif

  ! *****************************************************************
  ! *****************************************************************
  
  subroutine assemble_forcepert_param( N_x, N_y,               & 
       N_forcepert, forcepert_param )
    
    ! *sample* subroutine that demonstrates how pert_param can be 
    ! assembled for forcing perturbations
    !
    ! THIS SUBROUTINE IS NOT USED IN LDASsa - reichle, 8/8/2008
    !
    ! return N_force_pert, allocate and assemble structure forcepert_param
    !
    ! make sure order of fields is compatible with your driver
    !
    ! forcing field 1 = precip
    ! forcing field 2 = shortwave
    ! forcing field 3 = longwave
    ! forcing field 4 = air temperature
    !
    ! reichle, 11 Feb 2005
    ! reichle,  8 Jun 2005
    
    implicit none
    
    integer, intent(in) :: N_x, N_y
    
    integer, intent(out) :: N_forcepert 
    
    type(pert_param_type), dimension(:), pointer :: forcepert_param ! out
    
    ! -----------------------------------------------------------------
    
    ! forcing perturbation parameters
    
    ! # of forcing variables that are perturbed
    ! (currently pcp, sw, lw, tair)
    
    integer, parameter :: N_tmp  = 4

    integer, parameter :: ind_pcp  = 1
    integer, parameter :: ind_sw   = 2
    integer, parameter :: ind_lw   = 3
    integer, parameter :: ind_tair = 4
    
    character(40), parameter :: descr_pcp  = 'pcp'
    character(40), parameter :: descr_sw   = 'sw'
    character(40), parameter :: descr_lw   = 'lw'
    character(40), parameter :: descr_tair = 'tair'
    
    
    ! limit on range of random numbers:
    !  specify max absolute value allowed to be drawn from a standard
    !  normal distribution
    
    real, parameter :: std_normal_max =  2.5
    
    ! decide whether to ensure zero mean for synthetically generated errors
    ! (IMPORTANT: this will only have an effect for N_ens>2!!!)
    
    logical, parameter :: zeromean = .true.
    
    ! Allow perturbations to be computed on coarsened grid?
    ! Coarse grid spacing automatically determined as a function of model 
    ! grid spacing and spatial correlation scales (see random_fields.F90)
    
    logical, parameter :: coarsen = .false.

    ! temporal correlation scale
    
    real, parameter :: tcorr = 10800    ! 86400               ! [s]
    
    ! horizontal correlation scales 
    
    real, parameter :: xcorr = 0.
    real, parameter :: ycorr = 0.
    
    ! perturbations are either
    !
    !   typ=0: additive, mean=0
    !   typ=1: multiplicative and lognormal, mean=1
    
    integer, parameter :: typ_pcp  = 1             
    real,    parameter :: std_pcp  = .3         
    
    integer, parameter :: typ_sw    = 1
    real,    parameter :: std_sw    = .15        
    
    integer, parameter :: typ_lw    = 0
    real,    parameter :: std_lw    = 15.       

    integer, parameter :: typ_tair  = 0
    real,    parameter :: std_tair  = 1.       
    
    ! correlation coefficients -1 <= rho <= 1     
    !
    ! (these numbers are made up and are not tested with any data!)
    
    real, parameter :: rho_pcp_sw   = -.5
    real, parameter :: rho_pcp_lw   =  .5
    real, parameter :: rho_pcp_tair = -.3
    real, parameter :: rho_sw_lw    = -.5
    real, parameter :: rho_sw_tair  =  .5
    real, parameter :: rho_lw_tair  =  .4
    
    ! ---------------------------------------------------------------------
    !
    ! local variables
    
    integer   :: i, j, k, l
    real      :: tmpreal
    real, dimension(N_tmp,N_tmp) :: tmpmat1, tmpmat2
        
    ! ---------------------------------------------------------------------
    !
    ! allocate forcepert_param (must not be associated at this time)

    if (associated(forcepert_param)) then
       write (*,*) 'assemble_forcepert_param(): this needs work...'
       write (*,*) 'stopping'
       stop
    end if

    N_forcepert = N_tmp
    
    call allocate_pert_param(N_forcepert, N_x, N_y, forcepert_param)
    
    ! ----------------------------------------
    !
    ! copy inputs into structure
    
    ! precip perturbations
    
    forcepert_param(ind_pcp)%descr          = descr_pcp
    forcepert_param(ind_pcp)%typ            = typ_pcp
    forcepert_param(ind_pcp)%std_normal_max = std_normal_max
    forcepert_param(ind_pcp)%zeromean       = zeromean
    forcepert_param(ind_pcp)%coarsen        = coarsen
    forcepert_param(ind_pcp)%tcorr          = tcorr
    forcepert_param(ind_pcp)%xcorr          = xcorr
    forcepert_param(ind_pcp)%ycorr          = ycorr
    
    forcepert_param(ind_pcp)%std            = std_pcp
    
    forcepert_param(ind_pcp)%ccorr(1,:,:)   = 1.
    forcepert_param(ind_pcp)%ccorr(2,:,:)   = rho_pcp_sw
    forcepert_param(ind_pcp)%ccorr(3,:,:)   = rho_pcp_lw
    forcepert_param(ind_pcp)%ccorr(4,:,:)   = rho_pcp_tair    
    
    ! shortwave perturbations
    
    forcepert_param(ind_sw)%descr          = descr_sw
    forcepert_param(ind_sw)%typ            = typ_sw
    forcepert_param(ind_sw)%std_normal_max = std_normal_max
    forcepert_param(ind_sw)%zeromean       = zeromean
    forcepert_param(ind_sw)%coarsen        = coarsen
    forcepert_param(ind_sw)%tcorr          = tcorr
    forcepert_param(ind_sw)%xcorr          = xcorr
    forcepert_param(ind_sw)%ycorr          = ycorr
    
    forcepert_param(ind_sw)%std            = std_sw
    
    forcepert_param(ind_sw)%ccorr(1,:,:)   = rho_pcp_sw
    forcepert_param(ind_sw)%ccorr(2,:,:)   = 1.
    forcepert_param(ind_sw)%ccorr(3,:,:)   = rho_sw_lw
    forcepert_param(ind_sw)%ccorr(4,:,:)   = rho_sw_tair
    
    ! longwave perturbations
    
    forcepert_param(ind_lw)%descr          = descr_lw
    forcepert_param(ind_lw)%typ            = typ_lw
    forcepert_param(ind_lw)%std_normal_max = std_normal_max
    forcepert_param(ind_lw)%zeromean       = zeromean
    forcepert_param(ind_lw)%coarsen        = coarsen
    forcepert_param(ind_lw)%tcorr          = tcorr
    forcepert_param(ind_lw)%xcorr          = xcorr
    forcepert_param(ind_lw)%ycorr          = ycorr
    
    forcepert_param(ind_lw)%std            = std_lw
    
    forcepert_param(ind_lw)%ccorr(1,:,:)   = rho_pcp_lw
    forcepert_param(ind_lw)%ccorr(2,:,:)   = rho_sw_lw
    forcepert_param(ind_lw)%ccorr(3,:,:)   = 1.
    forcepert_param(ind_lw)%ccorr(4,:,:)   = rho_lw_tair
    
    ! air temperature perturbations
    
    forcepert_param(ind_tair)%descr          = descr_tair
    forcepert_param(ind_tair)%typ            = typ_tair
    forcepert_param(ind_tair)%std_normal_max = std_normal_max
    forcepert_param(ind_tair)%zeromean       = zeromean
    forcepert_param(ind_tair)%coarsen        = coarsen
    forcepert_param(ind_tair)%tcorr          = tcorr
    forcepert_param(ind_tair)%xcorr          = xcorr
    forcepert_param(ind_tair)%ycorr          = ycorr
    
    forcepert_param(ind_tair)%std            = std_tair
    
    forcepert_param(ind_tair)%ccorr(1,:,:)   = rho_pcp_tair
    forcepert_param(ind_tair)%ccorr(2,:,:)   = rho_sw_tair
    forcepert_param(ind_tair)%ccorr(3,:,:)   = rho_lw_tair
    forcepert_param(ind_tair)%ccorr(4,:,:)   = 1.
    
    ! -------------------------------------------------------------
    !
    ! set mean and (if needed) modify standard deviation according to 'typ'
    ! (additive or multiplicative perturbations)    
    
    do k=1,N_forcepert
       
       select case (forcepert_param(k)%typ)
          
       case (0)         ! additive (mean=0, std as above)
          
          forcepert_param(k)%mean  = 0.
          
       case (1)         ! multiplicative and lognormal (mean=1)
          
          do i=1,N_x
             do j=1,N_y
                
                tmpreal = forcepert_param(k)%std(i,j) 
                
                tmpreal = log( 1. + tmpreal**2)
                
                forcepert_param(k)%mean(i,j) = - .5*tmpreal
                forcepert_param(k)%std(i,j)  = sqrt(tmpreal)
                
             end do
          end do
          
       case default
          
          write (*,*) 'assemble_forcepert_param(): encountered unknown'
          write (*,*) 'type of error, stopping...'
          stop
          
       end select
       
    end do
    
    
    
      ! echo part of forcepert_param (mean, std, and ccorr for i=1, j=1 only):
  if(root_logit) then
     do i=1,N_forcepert
     
        write (logunit,*) 'forcepert_param(',i,')%descr=', &
          forcepert_param(i)%descr
        write (logunit,*) 'forcepert_param(',i,')%typ=', &
          forcepert_param(i)%typ
        write (logunit,*) 'forcepert_param(',i,')%zeromean=', &
          forcepert_param(i)%zeromean
        write (logunit,*) 'forcepert_param(',i,')%coarsen=', &
          forcepert_param(i)%coarsen
        write (logunit,*) 'forcepert_param(',i,')%std_normal_max=', &
          forcepert_param(i)%std_normal_max
        write (logunit,*) 'forcepert_param(',i,')%xcorr=', &
          forcepert_param(i)%xcorr
        write (logunit,*) 'forcepert_param(',i,')%ycorr=', &
          forcepert_param(i)%ycorr
        write (logunit,*) 'forcepert_param(',i,')%tcorr=', &
          forcepert_param(i)%tcorr
     
        write (logunit,*) 'forcepert_param(',i,')%mean(1,1)=', &
          forcepert_param(i)%mean(1,1)
        write (logunit,*) 'forcepert_param(',i,')%std(1,1)=', &
          forcepert_param(i)%std(1,1)
     
        do j=1,N_forcepert
        
           write (logunit,*) 'forcepert_param(',i,')%ccorr(',j,',1,1)=', &
             forcepert_param(i)%ccorr(j,1,1)
        end do
     end do
 endif ! root_logit
    
    
    

    ! -------------------------------------------------------------
    !
    ! compute sqrt of correlation matrix for each grid point
    
    do i=1,N_x
       do j=1,N_y

          ! extract local correlation matrix for grid point (i,j)
          
          do k=1,N_forcepert
             do l=1,N_forcepert
                
                tmpmat1(k,l) = forcepert_param(k)%ccorr(l,i,j)
                
             end do
          end do
          
          ! compute sqrt of local correlation matrix
          
          call get_sqrt_corr_matrix( N_forcepert, tmpmat1, tmpmat2 )
          
          ! overwrite cross-correlations in forcepert_param with square 
          ! root of cross-correlation matrix
          
          do k=1,N_forcepert
             do l=1,N_forcepert
                
                forcepert_param(k)%ccorr(l,i,j) = tmpmat2(k,l)
                
             end do
          end do
          
       end do
    end do
    
  end subroutine assemble_forcepert_param

  ! ************************************************************************
  
#if 0
  subroutine get_sqrt_corr_matrix( N, rho12, rho13, rho23, A )
    
    ! get sqrt of correlation matrix
    !
    ! correlation matrix has diagonal=variance=1
    !
    ! corr_matrix = [ 1 rho12 rho13; rho12 1 rho23; rho13 rho23 1 ]
    !
    ! A is lower tri-angular
    !
    ! A = sqrt_corr_matrix = [ a11 0 0; a12 a22 0; a13 a23 a33 ]
    !
    ! A*transpose(A) = corr_matrix 
    !
    ! so far only implemented for 3-by-3 correlation matrices with all
    !  nonnegative eigenvalues (no check for eigenvalues!)
    !
    ! reichle, 25 Jan 2005
    
    implicit none
    
    integer, intent(in) :: N
    real,    intent(in) :: rho12, rho13, rho23
    
    real, intent(out), dimension(N,N) :: A
    
    ! ------------------------------------------------------------------
    
    if (N/=3) then
       
       write (*,*) 'get_sqrt_corr_matrix() not implemented for N<=3'
       
    else
       
       A(3,3) = 1. - rho12*rho12  ! temporary 
       
       A(1,1) = 1.
       
       A(1,2) = 0.
       A(2,1) = rho12
       
       A(2,2) = sqrt(A(3,3))
       
       A(1,3) = 0.
       A(3,1) = rho13
       A(2,3) = 0.
       A(3,2) = (rho23-rho12*rho13)/A(2,2)
       
       A(3,3) =                                                          &
            sqrt(                                                        &
            (A(3,3) - rho13*rho13 - rho23*rho23 + 2*rho12*rho13*rho23)   &
            /                                                            &
            A(3,3) )
       
    end if
    
  end subroutine get_sqrt_corr_matrix
#endif  
    
  ! ************************************************************************

  subroutine get_sqrt_corr_matrix( N, A, S )
    
    ! get sqrt S of real, symmetric (correlation) matrix A
    !
    ! NOTE: there is no check that A is indeed symmetric!
    !
    ! A = S*transpose(S)
    !
    ! reichle, 7 Jun 2005
    
    implicit none
    
    integer, intent(in)                  :: N
    real,    intent(in),  dimension(N,N) :: A
    real,    intent(out), dimension(N,N) :: S

    ! local
    
    integer :: j

    real, dimension(N) :: D

    character(len=*), parameter :: Iam = 'get_sqrt_corr_matrix'
    character(len=400) :: err_msg
    
    ! ------------------------------------------------------------------
    
    call jacobi(A,N,D,S)
        
    do j=1,N
       
       if (D(j)<0.) then
          err_msg = 'negative eigenvalue found, invalid corr matrix'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       S(:,j) = S(:,j)*sqrt(D(j))
       
    end do
    
  end subroutine get_sqrt_corr_matrix
  
  ! ************************************************************************

  function find_rf(Nx, Ny, Nx_fft, Ny_fft, comm) result (rf)
    type(random_fields), pointer :: rf 
    integer, intent(in) :: Nx, Ny, Nx_fft, Ny_fft
    integer, optional, intent(in) :: comm

    type(StringRandom_fieldsMapIterator) :: iter
    Character(len=:), allocatable :: id_string
    type(random_fields) :: rf_tmp

    id_string = i_to_string(Nx)//":"//i_to_string(Ny)//":"//i_to_string(Nx_fft)//":"//i_to_string(Ny_fft)
    iter = random_fieldsMap%find(id_string)
    if (iter == random_fieldsMap%end() ) then
      rf_tmp = random_fields(Nx, Ny, Nx_fft, Ny_fft, comm=comm)
      call random_fieldsMap%insert(id_string, rf_tmp)
      iter = random_fieldsMap%find(id_string) 
    endif
    rf => iter%value()

  end function

  subroutine clear_rf()
    type(StringRandom_fieldsMapIterator) :: iter     
    type(random_fields), pointer :: rf_ptr 
    iter = random_fieldsMap%begin()
    do while (iter /= random_fieldsMap%end())
      rf_ptr => iter%value()
      call rf_ptr%finalize()
      ! remove the files
      call random_fieldsMap%erase(iter)
      iter = random_fieldsMap%begin()
    enddo
  end subroutine

end module land_pert_routines

! ***************************************************************************
! ***************************************************************************

#if 0

program test
  
  ! test land_pert_types
  
  use land_pert_types
  
  implicit none
  
  integer :: i, j, k, l, Nx=4, Ny=2, Nf=3
  
  type(pert_param_type), dimension(:), pointer :: pp
  
  character(40) :: tmpstr
  
  ! --------------------------------------------------------------------
  
  nullify(pp)
    
  call allocate_pert_param(Nf, Nx, Ny, pp)

  ! assemble

  do k=1,Nf
     
     write (tmpstr,'(i4.4)') k
     
     tmpstr = 'descr' // trim(tmpstr)
     
     pp(k)%descr = tmpstr
     
     do i=1,Nx
        do j=1,Ny
           
           pp(k)%mean(i,j) = k**2
           
           do l=1,Nf
              pp(k)%ccorr(l,i,j) = real(i*j)/real(k*l)
           end do
           
        end do
     end do
     
  end do
  
  ! write 

  do k=1,Nf
     
     write (*,*) pp(k)%descr
     do i=1,Nx
        write (*,*) pp(k)%mean(i,:)
     end do
     write (*,*) 'ccorr', k
     do i=1,Nx
        do j=1,Ny
           write (*,*) (pp(k)%ccorr(l,i,j), l=1,Nf)
        end do
     end do
     
  end do
     
end program test

#endif

! -------------------------------------------------------------

#if 0

program test
  
  use land_pert_routines
  use land_pert_types
  use nr_jacobi
  
  implicit none
  
  integer, parameter :: N = 3
  
  real, dimension(N,N) :: A, V, S
  
  real, dimension(N)   :: D
  
  integer :: i,j
  
  ! ------------------------------------------------
  
  A(1,1) = 1.
  A(1,2) = -.8
  A(1,3) =  .5
  A(2,2) = 1.
  A(2,3) = -.5
  A(3,3) = 1.
  A(2,1) = A(1,2)
  A(3,2) = A(2,3)
  A(3,1) = A(1,3)
  
  
  do i=1,N
     write (*,*) A(i,:)     
  end do
  
  call jacobi(A,N,D,V)
    
  do i=1,N
     write (*,*) A(i,:)     
  end do
  
  do i=1,N
     write (*,*) V(i,:)     
  end do
  
  write (*,*) D
  
  
  call get_sqrt_corr_matrix(N,A,S)
  
  do i=1,N
     write (*,*) S(i,:)     
  end do
  
  call get_sqrt_corr_matrix(1,A(1,3),S(1,1))
  
  write (*,*) S(1,1)     
  
  call get_sqrt_corr_matrix(2,A(1:2,1:2),S(1:2,1:2))
  
  do i=1,2
     write (*,*) S(i,1:2)     
  end do
  
end program test
  
#endif

! ------------------------------------------------------------------

#if 0

program test
  
  ! driver routine for testing subroutine precip_rad_perturb()

  ! ifort clsm_ensdrv_glob_param.o random_fields.F90 land_pert_types.o nr_ran2_gasdev.o nr_jacobi.o nr_fft.o land_pert.F90

  ! reichle, 24 Jan 2005
  ! reichle, 16 Feb 2005
  ! reichle,  3 Dec 2013 - updated to new interface of get_pert()

  ! reichle,  9 Apr 2021 - adapted to new order of array dimensions in pert fields after unification of get_pert()
  !                        and propagate_pert().  MAY NEED MORE WORK!!!
  
  ! ------------------------------------------------------------
  
  use land_pert_routines
  use land_pert_types
  use nr_ran2_gasdev
  use tile_coord_types, ONLY: grid_def_type 

  
  implicit none
  
  integer, parameter :: N_ens =  30
  integer, parameter :: N_x   =  4  
  integer, parameter :: N_y   =  8  
  
  integer, parameter :: N_domain = 1
  
  integer :: N_t, i, j, tt, n, N_forcepert
  
  real :: dx, dy, dtstep
  
  integer, dimension(N_ens,N_domain) :: init_Pert_rseed
  integer, dimension(N_ens)          :: ens_id
  integer, dimension(N_domain)       :: domain_id
  
  integer, dimension(:,:), allocatable  :: Forcepert_rseed
  
  real, dimension(:,:,:,:), allocatable :: Forcepert
  real, dimension(:,:,:,:), allocatable :: Forcepert_ntrmdt
  
  type(pert_param_type), dimension(:), pointer :: forcepert_param

  type(grid_def_type) :: pert_grid
  
  ! -----------------------------------------
  
  nullify(forcepert_param)
  
  write (*,*) NRANDSEED
  
  dx = 2.5
  dy = 2.
  
  dtstep = 1800.
  
  N_t = 1000
  
  ! open files for output
  
  open(991, file='tmp_precip.dat', form='formatted', status='unknown', &
       action='write')
  open(992, file='tmp_swdn.dat',   form='formatted', status='unknown', &
       action='write')
  open(993, file='tmp_lwdn.dat',   form='formatted', status='unknown', &
       action='write')
  open(994, file='tmp_tair.dat',   form='formatted', status='unknown', &
       action='write')
  
  ! initialize
  
  call assemble_forcepert_param(N_x, N_y, N_forcepert, forcepert_param)
  
  ! echo part of forcepert_param (mean, std, and ccorr for i=1, j=1 only):
  
  do i=1,N_forcepert
     
     write (*,*) 'forcepert_param(',i,')%descr=', &
          forcepert_param(i)%descr
     write (*,*) 'forcepert_param(',i,')%typ=', &
          forcepert_param(i)%typ
     write (*,*) 'forcepert_param(',i,')%zeromean=', &
          forcepert_param(i)%zeromean
     write (*,*) 'forcepert_param(',i,')%coarsen=', &
          forcepert_param(i)%coarsen
     write (*,*) 'forcepert_param(',i,')%std_normal_max=', &
          forcepert_param(i)%std_normal_max
     write (*,*) 'forcepert_param(',i,')%xcorr=', &
          forcepert_param(i)%xcorr
     write (*,*) 'forcepert_param(',i,')%ycorr=', &
          forcepert_param(i)%ycorr
     write (*,*) 'forcepert_param(',i,')%tcorr=', &
          forcepert_param(i)%tcorr
     
     write (*,*) 'forcepert_param(',i,')%mean(1,1)=', &
          forcepert_param(i)%mean(1,1)
     write (*,*) 'forcepert_param(',i,')%std(1,1)=', &
          forcepert_param(i)%std(1,1)
     
     do j=1,N_forcepert
        
        write (*,*) 'forcepert_param(',i,')%ccorr(',j,',1,1)=', &
             forcepert_param(i)%ccorr(j,1,1)
        
     end do
     
  end do
  
  ! -------------------------------------------------------------------

  allocate(Forcepert_rseed(NRANDSEED,N_ens))
  
  allocate(Forcepert_ntrmdt(N_x,N_y,N_forcepert,N_ens))
  allocate(Forcepert(       N_x,N_y,N_forcepert,N_ens))
  
  Forcepert_ntrmdt = 0.    ! initialize just in case (should not be needed)
  Forcepert        = 0.    ! initialize just in case (should not be needed)
  
  do n=1,N_ens
     ens_id(n) = n-1
  end do
  
  do n=1,N_domain
     domain_id(n) = n
  end do

  ! get different negative integer for each ensemble member and
  ! each domain
  
  !call get_init_Pert_rseed( N_ens, N_domain, ens_id, domain_id, &
  !     init_Pert_rseed )
  call get_init_Pert_rseed( N_ens, ens_id, domain_id, init_Pert_rseed )


  ! initialize first row of Forcepert_rseed (for first domain)
  
  Forcepert_rseed(1,:) = init_Pert_rseed( :, 1)
  
  ! initial call to get_pert  (for first domain)
  !
  ! this initializes Forcepert_ntrmdt and the rest of Forcepert_rseed 
  !
  ! NOTE: after initial call to get_pert use restart files 
  !       for Forcepert_rseed and Forcepert_ntrmdt to continue
  !       the perturbation time series whenever the land model
  !       integration is interrupted

  pert_grid%N_lon = N_x
  pert_grid%N_lat = N_y
  pert_grid%dlon  = dx
  pert_grid%dlat  = dy
    
  call get_pert(                                           &
       N_forcepert, N_ens, pert_grid,                       &
       pert_grid, dtstep,                                     &
       forcepert_param,                                    &
       Forcepert_rseed,                                    &
       Forcepert_ntrmdt,                                   &
       Forcepert,                                          &                
       initialize_rseed=.true.,                            &
       initialize_ntrmdt=.true.            )
  
  
  write (*,*) 'Forcepert_rseed=', Forcepert_rseed  
  
  ! output 
  
  do i=1,N_forcepert
     
     write (990+i, '(3312(e13.5))') Forcepert(:,:,i,:)
     
  end do
  
  write (*,*) Forcepert(:,1,1,1)
  
  ! loop through time
  
  do tt=1,N_t
     
     call get_pert(                                           &
          N_forcepert, N_ens, pert_grid,                      &
          pert_grid, dtstep,                                     &
          forcepert_param,                                    &
          Forcepert_rseed,                                    &
          Forcepert_ntrmdt,                                   &
          Forcepert                )
     
     ! output 
     
     do i=1,N_forcepert
        
        write (990+i, '(3312(e13.5))') Forcepert(:,:,i,:)
        
     end do
     
     write (*,*) Forcepert(:,1,1,1)
     
  end do
  
end program test

#endif

! =============== EOF =================================================
