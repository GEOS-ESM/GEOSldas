

module clsm_ensdrv_drv_routines
  
  ! collection of subroutines for enkf_driver written in f90
  ! reichle, 10 May 2005
  !
  ! reichle, 13 Aug 2008 - moved forcing subroutines into 
  !                         clsm_ensdrv_force_routines.F90
  ! reichle, 28 Oct 2008 - added soilcls30 and soilcls100
  !                      - optimized restart-to-exp-domain mapping in initialize_model()
  ! reichle,  5 Apr 2013 - revised treatment of output collections

  use LDAS_ensdrv_globals,           ONLY:     &
       logunit,                                   &
       logit,                                     &
       nodata_generic,                            &
       nodata_tol_generic

  use catch_constants,                  ONLY:     &
       N_snow => CATCH_N_SNOW,                    &
       N_gt   => CATCH_N_GT

  use catch_incr,                        ONLY:     &
       check_catch_progn
  
  use MAPL_ConstantsMod,                ONLY:     &
       stefan_boltzmann => MAPL_STFBOL,           &
       alhe             => MAPL_ALHL,             &
       alhs             => MAPL_ALHS,             &
       alhm             => MAPL_ALHF,             &
       Tzero            => MAPL_TICE       
  
  use LDAS_DriverTypes,                     ONLY:     &
       met_force_type,                            &
       veg_param_type,                            &
       alb_param_type,                            &
       bal_diagn_type,                            &
       assignment (=),                            &
       operator (+),                              &  
       operator (*)

  use catch_types,                      ONLY:     &
       cat_param_type,                            &
       cat_progn_type,                            &
       cat_diagS_type,                            &
       cat_diagF_type,                            &
       catprogn2wesn,                             &
       catprogn2htsn,                             &
       catprogn2ghtcnt

  use LDAS_DateTimeMod,                   ONLY:     &
       date_time_type,                            &
       datetime2_minus_datetime1

  use LDAS_ensdrv_mpi,                  ONLY:     &
       mpicomm,                            &
       mpierr,                                    &
       numprocs,                                  &
       master_proc

  use catchment_model,                  ONLY:     &
       catch_calc_tsurf,                          &
       catch_calc_etotl

  use lsm_routines,                     ONLY:     &
       catch_calc_soil_moist,                     &
       catch_calc_tp,                             &
       catch_calc_wtotl

  use StieglitzSnow,                    ONLY:     &
       StieglitzSnow_calc_asnow,                  &
       StieglitzSnow_calc_tpsnow

  use LDAS_ExceptionsMod,                  ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  implicit none

  include 'mpif.h'

  private
  
  public :: spin_stuff
  public :: check_cat_progn
  public :: recompute_diagS
!!  public :: interpolate_to_timestep
!!  public :: zenith
  public :: remove_snow
  public :: balance_calcs
  public :: l2f_real
  public :: f2l_real
  public :: f2l_real8
  public :: f2l_logical

  character(10), private :: tmpstring10
  character(40), private :: tmpstring40

contains
  
  ! ********************************************************************
    
  subroutine spin_stuff( start_time, end_time, N_ens, N_force_pert, N_progn_pert,  &
       spin_loop, restart  )
    
    implicit none
    
    type(date_time_type), intent(in) :: start_time, end_time
    
    integer, intent(in) :: N_ens, N_force_pert, N_progn_pert
    
    integer, intent(inout) :: spin_loop
    
    logical, intent(inout) :: restart
    
    ! local
    character(len=*), parameter :: Iam = 'spin_stuff'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------------
        
    ! consistency checks
    
    ! make sure end time is exactly n years after start time
    
    if ( start_time%month /= end_time%month .or. &
         start_time%day   /= end_time%day   .or. &
         start_time%hour  /= end_time%hour  .or. &
         start_time%min   /= end_time%min   .or. &
         start_time%sec   /= end_time%sec   ) then
       err_msg = 'spin up only for full years'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! do not allow N_ens>1 during spin-up
    
    if ((N_ens>1) .or. (N_force_pert>1) .or. (N_progn_pert>1)) then
       err_msg = 'spin-up only for N_ens=1 and N_force_pert=N_progn_pert=0'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
           
    ! use restart files
    
    if (spin_loop==0) then
       
       restart = .false.
       
    else
       
       restart = .true.
       
    end if
    
    ! augment counter for spin years
    
    spin_loop = spin_loop + 1
    
    ! echo spin_loop
    
    if (logit) write (logunit,*)
    if (logit) write (logunit,*) 'at beginning of spin_loop ', spin_loop
    
  end subroutine spin_stuff
    
  ! *********************************************************************
  
  subroutine check_cat_progn( N_cat, cat_param, cat_progn )
    
    ! wrapper for subroutine check_catch_progn() which has been 
    ! moved to "catch_iau.F90" in GEOScatch_GridComp - reichle, 3 Apr 2012
    
    implicit none
    
    integer, intent(in) :: N_cat
    
    type(cat_param_type), dimension(N_cat), intent(in)    :: cat_param
    
    type(cat_progn_type), dimension(N_cat), intent(inout) :: cat_progn
    
    ! ----------------------------------------------------------------
    
    ! local variables
    
    integer :: k 
    
    real, dimension(N_snow,N_cat) :: wesn
    real, dimension(N_snow,N_cat) :: htsn
    real, dimension(N_snow,N_cat) :: sndz

    real, dimension(N_gt,  N_cat) :: ghtcnt
    
    ! ----------------------------------------------------------------
    
    ! copy select cat_progn fields into 2-d arrays
    
    do k=1,N_gt
       
       GHTCNT(k,:) = cat_progn%ght(k)
       
    end do
    
    do k=1,N_snow
       
       WESN(k,:)   = cat_progn%wesn(k)
       HTSN(k,:)   = cat_progn%htsn(k)
       SNDZ(k,:)   = cat_progn%sndz(k)
       
    end do

    ! check for consistency and unphysical values
    
    call check_catch_progn( N_cat, cat_param%vegcls, cat_param%dzsf,         &
         cat_param%vgwmax,  cat_param%cdcr1, cat_param%cdcr2,                &
         cat_param%psis, cat_param%bee, cat_param%poros, cat_param%wpwet,    &
         cat_param%ars1, cat_param%ars2, cat_param%ars3,                     &
         cat_param%ara1, cat_param%ara2, cat_param%ara3, cat_param%ara4,     &
         cat_param%arw1, cat_param%arw2, cat_param%arw3, cat_param%arw4,     &
         cat_progn%tc1, cat_progn%tc2, cat_progn%tc4,                        &
         cat_progn%qa1, cat_progn%qa2, cat_progn%qa4,                        &
         cat_progn%capac, cat_progn%catdef,                                  &
         cat_progn%rzexc, cat_progn%srfexc,                                  &
         ghtcnt, wesn, htsn, sndz            )
    
    ! copy 2-d arrays back into cat_progn fields
    
    do k=1,N_gt
       
       cat_progn%ght(k)  = ghtcnt(k,:)
       
    end do
    
    do k=1,N_snow
       
       cat_progn%wesn(k) = WESN(k,:)
       cat_progn%htsn(k) = HTSN(k,:)
       cat_progn%sndz(k) = SNDZ(k,:)
       
    end do

  end subroutine check_cat_progn

  ! *********************************************************************
  
  subroutine recompute_diagS( N_catd, cat_param, cat_progn, cat_diagS )
    
    ! replace cat_diagS with updated diagnostics
    !
    ! typically call after prognostics perturbations, EnKF update, and/or cat 
    !  bias correction
    ! 
    ! IMPORTANT: cat_progn is intent(inout) because srfexc, rzexc, and catdef 
    !  fields *might* be changed! Such changes can be prevented by first calling 
    !  subroutine check_cat_progn().
    !
    ! reichle, 20 Oct 2004
    ! reichle, 14 Feb 2013 - added recomputation of soil temperature
    ! reichle, 30 Oct 2013 - moved from clsm_ensupd_upd_routines.F90
    !                      - removed dependency on "update_type"
    !                      - added recompute of snow diagnostics
    
    integer, intent(in) :: N_catd
    
    ! note: cat_progn must be "inout" because call to calc_soil_moist
    !       might reset inconsistent sets of srfexc, rzexc, catdef
    
    type(cat_param_type), dimension(N_catd), intent(in)    :: cat_param    
    type(cat_progn_type), dimension(N_catd), intent(inout) :: cat_progn
    type(cat_diagS_type), dimension(N_catd), intent(inout) :: cat_diagS
    
    ! local variables
    
    integer :: i
    
    real, dimension(     N_catd) :: ar4, fices
    
    real, dimension(N_gt,N_catd) :: tp
    
    ! ------------------------------------------------------------------
    
    ! update soil moisture diagnostics
    
    ! note that the call to calc_soil_moist resets srfexc, rzexc, and
    ! catdef if they are not consistent!
    
    ! updated to new interface - reichle, 3 Apr 2012
    
    call catch_calc_soil_moist(                                     &
         N_catd,                                                    &
         cat_param%vegcls, cat_param%dzsf,  cat_param%vgwmax,       &
         cat_param%cdcr1,  cat_param%cdcr2, cat_param%psis,         &
         cat_param%bee,    cat_param%poros, cat_param%wpwet,        &
         cat_param%ars1,   cat_param%ars2,  cat_param%ars3,         &
         cat_param%ara1,   cat_param%ara2,  cat_param%ara3,         &
         cat_param%ara4,   cat_param%arw1,  cat_param%arw2,         &
         cat_param%arw3,   cat_param%arw4,                          &
         cat_progn%srfexc, cat_progn%rzexc, cat_progn%catdef,       &
         cat_diagS%ar1,    cat_diagS%ar2,   ar4,                    &
         cat_diagS%sfmc,   cat_diagS%rzmc,  cat_diagS%prmc)            
    

    ! update snow cover fraction and snow temperatures
    
    call StieglitzSnow_calc_asnow(                                  &
         N_snow, N_catd, catprogn2wesn(N_catd,cat_progn), cat_diagS%asnow )

    do i=1,N_snow
              
       call StieglitzSnow_calc_tpsnow( N_catd,                      &
            cat_progn(1:N_catd)%htsn(i),                            &
            cat_progn(1:N_catd)%wesn(i),                            &
            cat_diagS(1:N_catd)%tpsn(i),                            &
            fices )
       
       cat_diagS%tpsn(i) = cat_diagS%tpsn(i) + Tzero   ! convert to Kelvin
       
    end do
    
    ! update surface temperature
    
    ! updated to new interface, 
    ! need ar1, ar2, ar4 from call to catch_calc_soil_moist() above
    ! - reichle, 3 Apr 2012
    
    call catch_calc_tsurf( N_catd,                                  &
         cat_progn%tc1, cat_progn%tc2, cat_progn%tc4,               &
         catprogn2wesn(N_catd,cat_progn),                           &
         catprogn2htsn(N_catd,cat_progn),                           &
         cat_diagS%ar1, cat_diagS%ar2, ar4,                         &
         cat_diagS%tsurf )
 
    ! update soil temperature
    
    ! NOTE: "tp" is returned in CELSIUS [for consistency w/ catchment.F90]
    
    call catch_calc_tp( N_catd, cat_param%poros,                    &
         catprogn2ghtcnt(N_catd,cat_progn), tp )
    
    do i=1,N_gt
       
       cat_diagS(:)%tp(i) = tp(i,:)
       
    end do
    
  end subroutine recompute_diagS
  
  ! ********************************************************************
   
  ! ******************************************************************
  
!!   subroutine interpolate_to_timestep(                             &
!!        N_catd, vegcls, lat, lon, zenav, date_time_new,            &
!!        force_time_old, force_dtstep,                              &
!!        grn_time_old,   grn_time_new,                              &
!!        lai_time_old,   lai_time_new,                              &
!!       alb_time_old,   alb_time_new,                              &
!!        mf_old,         mf_new,                                    &
!!        veg_param_old,  veg_param_new,                             &
!!        alb_param_old,  alb_param_new,                             &
!!        mf_ntp, sunang_ntp, veg_param_ntp, alb_param_ntp )
!!     
!!     ! Interpolates the forcing, vegetation and albedo data to current timestep.
!!     !
!!     !  date_time_new = date_time at end of model integration time step
!!     !
!!     ! "mf"   = "met_force"
!!     !
!!     ! "mf_old" = at old forcing time 
!!     ! "mf_new" = at new forcing time
!!     ! "mf_ntp" = at current ("interpolated") time
!!     !
!!     ! NOTE: time avg radiative fluxes for the interval between "old" 
!!     !       and "new" time must be stored in mf_old
!!     !
!!     ! reichle, 14 May 2003
!!     ! reichle, 11 Sep 2007 - albedo changes
!!     ! reichle, 23 Feb 2009 - add ParDrct, ParDffs from MERRA forcing 
!!     ! reichle,  6 Mar 2009 - added fractional day-of-year (fdofyr) for monthly interp
!!     !                      - deleted ParDrct, ParDffs after testing found no impact
!!     ! reichle, 20 Dec 2011 - reinstated "PARdrct", "PARdffs" for MERRA-Land file specs
!!     !                      - cleanup
!!     ! reichle, 26 Jul 2013 - revised GRN, LAI and albedo scaling parameter inputs
!! 
!!     implicit none
!! 
!!     integer,                                 intent(in)  :: N_catd
!!     
!!     integer,              dimension(N_catd), intent(in)  :: vegcls
!!     
!!     real,                 dimension(N_catd), intent(in)  :: lat, lon, zenav
!! 
!!     type(date_time_type),                    intent(in)  :: date_time_new
!!     type(date_time_type),                    intent(in)  :: force_time_old
!! 
!!     integer,                                 intent(in)  :: force_dtstep
!!         
!!     type(date_time_type),                    intent(in)  :: grn_time_old, grn_time_new
!!     type(date_time_type),                    intent(in)  :: lai_time_old, lai_time_new
!!     type(date_time_type),                    intent(in)  :: alb_time_old, alb_time_new
!!     
!!     type(met_force_type), dimension(N_catd), intent(in)  :: mf_old
!!     type(met_force_type), dimension(N_catd), intent(in)  :: mf_new
!!     
!!     type(veg_param_type), dimension(N_catd), intent(in)  :: veg_param_old
!!     type(veg_param_type), dimension(N_catd), intent(in)  :: veg_param_new
!!     
!!     type(alb_param_type), dimension(N_catd), intent(in)  :: alb_param_old
!!     type(alb_param_type), dimension(N_catd), intent(in)  :: alb_param_new
!!         
!!     type(met_force_type), dimension(N_catd), intent(out) :: mf_ntp
!!     
!!     real,                 dimension(N_catd), intent(out) :: sunang_ntp
!!     
!!     type(veg_param_type), dimension(N_catd), intent(out) :: veg_param_ntp
!!     
!!     type(alb_param_type), dimension(N_catd), intent(out) :: alb_param_ntp
!!         
!!     ! ----------------
!!     
!!     ! local variables
!!     
!!     real, parameter :: min_grn = 0.0001 ! per GEOS_CatchGridComp.F90 (Ganymed-4_0)
!!     real, parameter :: min_lai = 0.0001 ! per GEOS_CatchGridComp.F90 (Ganymed-4_0)
!!     real, parameter :: min_zth = 0.01   ! per testing Feb 2009 (see below)
!!     
!!     integer :: n, secs_since_old, secs_in_day  
!!     
!!     real    :: zth, slr, w, w_old, w_new, tmpreal
!! 
!!     character(len=*), parameter :: Iam = 'interpolate_to_timestep'
!!         
!!     ! ------------------------------------------------------------
!!     !
!!     ! met forcing interpolation
!!     !
!!     ! get secs_in_day from hh:mm:ss
!!     
!!     secs_in_day = date_time_new%hour*3600 + date_time_new%min*60 &
!!          + date_time_new%sec
!!     
!!     ! weight for forcing "states" interpolation 
!!     ! (temperature, humidity, pressure, wind) 
!!     
!!     secs_since_old = datetime2_minus_datetime1( force_time_old, date_time_new )
!!     
!!     ! use integer division such that w changes from 0. to 1.
!!     ! halfway through the current forcing interval, that is,  
!!     !
!!     !   w = 0.  if                  secs_since_old <  force_dtstep/2
!!     !   w = 0.5 if                  secs_since_old == force_dtstep/2
!!     !   w = 1.  if force_dtstep/2 < secs_since_old <= force_dtstep
!!     !
!!     ! For example, using 15 min model time steps and hourly forcing, 
!!     ! the time interpolation weights are as follows:
!!     !
!!     !   secs_since_old:         900        1800        2700        3600 
!!     !   w:                        0.          0.5         1.          1.
!!     !
!!     ! Note that w=0.5 for secs_since_old==force_dtstep/2 (at the mid-point).
!! 
!!     if (secs_since_old==force_dtstep/2) then
!!        
!!        w = 0.5
!!        
!!     else
!!        
!!        w = real( (secs_since_old-1)/(force_dtstep/2) )
!!        
!!     end if
!!     
!!     ! ---------------------
!!     
!!     do n=1,N_catd
!!        
!!        ! initialize
!!        
!!        mf_ntp(n) = nodata_generic
!!        
!!        ! STATES
!!        !
!!        ! temperature, humidity, pressure and wind 
!!        
!!        mf_ntp(n)%Tair  = (1.-w)*mf_old(n)%Tair  + w*mf_new(n)%Tair
!!        mf_ntp(n)%Qair  = (1.-w)*mf_old(n)%Qair  + w*mf_new(n)%Qair
!!        mf_ntp(n)%Psurf = (1.-w)*mf_old(n)%Psurf + w*mf_new(n)%Psurf
!!        mf_ntp(n)%RefH  = (1.-w)*mf_old(n)%RefH  + w*mf_new(n)%RefH
!!        
!!        ! Wind
!! 
!!        ! LDASsa CVS tags between reichle-LDASsa_m2-10 and reichle-LDASsa_m2-13_p2 
!!        !  worked with *inst*lfo* and *tavg*lfo* G5DAS forcing (as opposed to just 
!!        !  *tavg*lfo* from MERRA) and G5DAS Wind was read from *inst*lfo* files. 
!!        !  But for G5DAS forcing (as for MERRA), Wind was treated as a time-average 
!!        !  field, which implied that, e.g., the instantaneous G5DAS Wind at 0z was 
!!        !  used to force the land at 0:20z, 0:40z, and 1z.
!!        ! In these tags, Wind was treated as a time-average field whenever
!!        !  force_dtstep<=3601, which included MERRA, G5DAS, RedArk, and CONUS
!!        !  forcing.
!!        ! To make things more consistent for G5DAS winds, the "if" statement now
!!        !  checks whether Wind at date_time_new is available (which is true G5DAS
!!        !  and false for MERRA).  The revised "if" statement does not change how
!!        !  Wind is interpolated in MERRA (because Wind is unavailable at date_time_new)
!!        !  and in GLDAS, Princeton, and GWSP (because force_dtstep>3601).
!!        !  The revised "if" statement does change how CONUS and RedArk Wind data
!!        !  are interpolated (now treated as instantaneous, previously treated as 
!!        !  time-average).
!!        ! In summary, as of this change, Wind is treated as instantaneous for all
!!        !  forcing data *except* MERRA.
!!        !
!!        ! - reichle, 31 Jan 2014
!!        
!!        ! if (force_dtstep_real > 3601.) then       
!!        if (abs(mf_new(n)%Wind-nodata_generic)>nodata_tol_generic) then
!!           
!!           ! treat Wind as instantaneous fields (all forcing data sets *except* MERRA)
!!           
!!           mf_ntp(n)%Wind  = (1.-w)*mf_old(n)%Wind  + w*mf_new(n)%Wind
!!           
!!        else
!!           
!!           ! treat Wind as time-average fields (MERRA)
!!           
!!           mf_ntp(n)%Wind  = mf_old(n)%Wind
!!           
!!        end if
!!        
!!        ! FLUXES
!!        
!!        ! precipitation
!!        
!!        mf_ntp(n)%Rainf_C = mf_old(n)%Rainf_C       
!!        mf_ntp(n)%Rainf   = mf_old(n)%Rainf              
!!        mf_ntp(n)%Snowf   = mf_old(n)%Snowf
!!        
!!        ! incoming radiation
!!        
!!        mf_ntp(n)%LWdown = mf_old(n)%LWdown
!!        
!!        call solar(lon(n),lat(n),date_time_new%dofyr,secs_in_day,zth,slr)
!!        
!!        ! changed min sun-angle from 0.01 to 0.0001 for consistency with CatchGridComp
!!        ! reichle, 23 Feb 2009
!!        ! changed min sun-angle back to 0.01 after testing
!!        ! reichle, 27 Feb 2009
!!        
!!        sunang_ntp(n) = max(zth, min_zth)  
!!               
!!        ! changed minimum SWdown to 0. from 0.00001 - reichle, 28 Aug 2008
!! 
!!        if (zth > 0.) then
!!           
!!           if (zenav(n) <= 0.) then
!!              call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'Problem with solar')
!!           end if
!! 
!!           tmpreal = zth/zenav(n)
!! 
!!           mf_ntp(n)%SWdown = mf_old(n)%SWdown*tmpreal
!!                     
!!           ! mf_ntp%SWnet only used if mf_old%SWnet is not no-data value; 
!!           ! protect multiplication for any no-data-value because it could
!!           ! fail (floating point excess) if huge number is used as nodata value
!!           
!!           if(abs(mf_old(n)%SWnet  -nodata_generic)>nodata_tol_generic)  &
!!                mf_ntp(n)%SWnet   = mf_old(n)%SWnet*tmpreal   
!!           
!!           if(abs(mf_old(n)%PARdrct-nodata_generic)>nodata_tol_generic) then
!!              
!!              ! assume that PARdffs is available whenever PARdrct is
!!           
!!              mf_ntp(n)%PARdrct = mf_old(n)%PARdrct*tmpreal  
!!              mf_ntp(n)%PARdffs = mf_old(n)%PARdffs*tmpreal  
!!              
!!           end if
!!           
!!        elseif ((zth <= 0.) .and. (zenav(n) <= 0.)) then
!!           
!!           mf_ntp(n)%SWdown  = max(0., mf_old(n)%SWdown)
!!           mf_ntp(n)%SWnet   = max(0., mf_old(n)%SWnet)    ! no-data handling done below
!!           mf_ntp(n)%PARdrct = max(0., mf_old(n)%PARdrct)  ! no-data handling done below
!!           mf_ntp(n)%PARdffs = max(0., mf_old(n)%PARdffs)  ! no-data handling done below
!!           
!!        else
!!           
!!           mf_ntp(n)%SWdown  = 0.
!!           mf_ntp(n)%SWnet   = 0.   ! no-data handling done below
!!           mf_ntp(n)%PARdrct = 0.   ! no-data handling done below
!!           mf_ntp(n)%PARdffs = 0.   ! no-data handling done below
!!           
!!        end if
!!        
!!        ! cap shortwave radiation at (cosine of) sun angle times solar constant
!!        ! reichle, 14 Aug 2002
!!        
!!        mf_ntp(n)%SWdown = min( mf_ntp(n)%SWdown, 1360.*sunang_ntp(n) )
!!        
!!        ! cap SWnet at SWdown
!!        
!!        mf_ntp(n)%SWnet  = min( mf_ntp(n)%SWnet, mf_ntp(n)%SWdown )
!!               
!!        ! reinstate no-data-values
!!        
!!        if(abs(mf_old(n)%SWnet-nodata_generic)<nodata_tol_generic)  &
!!             mf_ntp(n)%SWnet = nodata_generic
!!        
!!        if ( (abs(mf_old(n)%PARdrct-nodata_generic)<nodata_tol_generic) ) then
!!           
!!           ! assume that PARdffs is no-data whenever PARdrct is
!!           
!!           mf_ntp(n)%PARdrct = nodata_generic
!!           mf_ntp(n)%PARdffs = nodata_generic
!!           
!!        end if
!!        
!!     end do
!!     
!!     ! --------------------------------------------------------------
!!     
!!     ! Slowly varying vegetation (GRN, LAI) and albedo scaling parameters 
!!     
!!     ! vegetation parameters
!!     
!!     ! greenness (GRN)
!!     
!!     call get_time_interpolation_weights(                               &
!!          grn_time_old, grn_time_new, date_time_new, w_old, w_new )
!!     
!!     veg_param_ntp%grn =   w_old * veg_param_old%grn + w_new * veg_param_new%grn
!!     
!!     ! leaf area index (LAI)
!!     
!!     call get_time_interpolation_weights(                               &
!!          lai_time_old, lai_time_new, date_time_new, w_old, w_new )
!!     
!!     veg_param_ntp%lai =   w_old * veg_param_old%lai + w_new * veg_param_new%lai
!!     
!!     ! Avoid values that are too close to zero
!!     ! 
!!     ! Note: In earlier versions of LDASsa, this has been done to the original
!!     !       (monthly) data.  Here it is applied to the interpolated data.    
!!     
!!     veg_param_ntp%grn = max( veg_param_ntp%grn, min_grn )
!!     veg_param_ntp%lai = max( veg_param_ntp%lai, min_lai )
!!     
!!     ! ---------------------------
!! 
!!     ! albedo scaling parameter
!! 
!!     call get_time_interpolation_weights(                               &
!!          alb_time_old, alb_time_new, date_time_new, w_old, w_new )
!! 
!!     do n=1,N_catd
!!        
!!        alb_param_ntp(n) = w_old * alb_param_old(n)  + w_new * alb_param_new(n)
!!        
!!     end do
!!     
!!   end subroutine interpolate_to_timestep
  
  ! ***********************************************************************
  
  subroutine get_time_interpolation_weights(                         &
       date_time_old, date_time_new, date_time_ntp, w_old, w_new )
    
    ! compute linear time interpolation weights w_old and w_new such that:
    ! 
    !   f(date_time_ntp) = w_old * f(date_time_old) + w_new * f(date_time_new) 
    !
    ! reichle, 26 Jul 2013
    
    implicit none
    
    type(date_time_type), intent(in)  :: date_time_old, date_time_new, date_time_ntp
    
    real,                 intent(out) :: w_old, w_new
    
    ! local variables
    
    integer :: dt_old_ntp, dt_ntp_new
    
    ! -------------------------
    
    dt_old_ntp = datetime2_minus_datetime1( date_time_old, date_time_ntp )
    dt_ntp_new = datetime2_minus_datetime1( date_time_ntp, date_time_new  )
    
    w_new = real(dt_old_ntp)/real(dt_old_ntp+dt_ntp_new)
    
    w_old = 1. - w_new
        
  end subroutine get_time_interpolation_weights
    
  ! ***********************************************************************

!!   subroutine zenith(                                                   &
!!        force_time, force_dtstep, model_dtstep, N_catd, lon, lat, zenav )
!!     
!!     ! calculate average zenith angle over the period from force_time
!!     ! to force_time+force_dtstep with a time discretization of model_dtstep
!!     !
!!     ! reichle, 25 May 2005
!!     
!!     implicit none
!!     
!!     type(date_time_type), intent(in) :: force_time
!!     
!!     integer, intent(in) :: N_catd, force_dtstep, model_dtstep
!!     
!!     real, dimension(N_catd), intent(in) :: lat, lon
!!     
!!     real, dimension(N_catd), intent(out) :: zenav
!!     
!!     ! local variables
!!     
!!     integer :: i, n, secs_in_day, N_steps
!!     
!!     real    :: zth, dummy
!!     
!!     ! -----------------------------------------------------------------------
!!     
!!     zenav = 0.
!!     
!!     ! compute seconds in day at force_time
!!     
!!     secs_in_day = force_time%hour*3600+force_time%min*60+force_time%sec &
!!          + model_dtstep
!!     
!!     N_steps = force_dtstep/model_dtstep
!!     
!!     do i=1,N_steps
!!        
!!        do n=1,N_catd
!!           
!!           call solar(lon(n),lat(n),force_time%dofyr,secs_in_day,zth,dummy)
!!           
!!           if (zth <= 0.) zth=0.
!!           
!!           zenav(n) = zenav(n) + zth
!!           
!!        end do
!!        
!!        secs_in_day = secs_in_day + model_dtstep
!!        
!!     end do
!!     
!!     ! normalize
!!     
!!     zenav = zenav/real(N_steps)
!!     
!!   end subroutine zenith
    
  ! *************************************************************
  
  subroutine remove_snow( N_catd, N_ens, date_time, latitude, cat_progn ) 
    
    ! Remove all snow on August 16 (or Feb 16 for the Southern Hemisphere)
    ! to prevent massive accumulations in cold regions during spin-up.
    !
    ! - reichle, 28 May 2003
    
    implicit none
    
    integer, intent(in) :: N_catd, N_ens
    
    type(date_time_type), intent(in) :: date_time
    
    real, dimension(N_catd), intent(in) :: latitude
    
    type(cat_progn_type), dimension(N_catd,N_ens), intent(inout) :: cat_progn
    
    ! local variables
    
    integer :: n, i 
    
    ! --------------------------------------------------------------------
    
    ! Northern Hemisphere, remove snow on Aug 16
    
    if ( date_time%month == 8   .and.       &
         date_time%day   == 16  .and.       &
         date_time%hour  == 0   .and.       &
         date_time%min   == 0   .and.       &
         date_time%sec   == 0 )       then          
       
       if (logit) write (logunit,*) &
            'spin-up mode: removing snow in Northern Hemisphere ' &
            //  ' on date_time = ', date_time
       
       do n=1,N_catd
          
          if (latitude(n)>0.) then
             
             do i=1,N_snow
                
                cat_progn(n,:)%wesn(i)=0.
                cat_progn(n,:)%htsn(i)=0.
                cat_progn(n,:)%sndz(i)=0.
                
             end do
             
          end if
       end do
    end if
    
    ! Southern Hemisphere, remove snow on Feb 16

    if ( date_time%month == 2   .and.       &
         date_time%day   == 16  .and.       &
         date_time%hour  == 0   .and.       &
         date_time%min   == 0   .and.       &
         date_time%sec   == 0 )       then          
       
       if (logit) write (logunit,*) &
            'spin-up mode: removing snow in Southern Hemisphere ' &
            //  ' on date_time = ', date_time
       
       do n=1,N_catd
          
          if (latitude(n)<0.) then
             
             do i=1,N_snow
                
                cat_progn(n,:)%wesn(i)=0.
                cat_progn(n,:)%htsn(i)=0.
                cat_progn(n,:)%sndz(i)=0.
                
             end do

          end if
       end do
    end if
    
  end subroutine remove_snow
  
  ! *************************************************************
  
  subroutine balance_calcs(                                &
       option, N_catd, dtstep,                             &
       cat_param, cat_progn, cat_diagF, met_force,         &
       bal_diagn)
    
    ! mass and energy balance calculations before and after forecasting the
    ! catchment prognostic states and performing the EnKF update
    !
    ! upon input, "bal_diagn" must contain total water and energy stores
    ! ("wtotl" and "etotl") from previous model time step
    !
    ! does not rely on cat_diagF *before* model propagation (option='ini')
    !
    ! reichle, 28 May 2003
    !
    ! reichle+qliu, 12 Aug 2008: amended mass balance calculations (hsnacc,Snowf*ALHM)
    ! reichle,       9 Dec 2011: revised using new "bal_diagn" type
    ! reichle,       4 Jan 2012: corrected "wtotl" and moved "wtotl" and "etotl" calcs
    !                            to catch_diagn_routines.F90
    ! reichle,       3 Apr 2012: revised after moving Catchment diagnostic subroutines
    !                             into GEOScatch_GridComp/catchment.F90 
    
    implicit none
    
    character(3), intent(in) :: option
    
    integer, intent(in) :: N_catd
    
    real, intent(in) :: dtstep
    
    type(cat_param_type), dimension(N_catd), intent(in)    :: cat_param
    type(cat_progn_type), dimension(N_catd), intent(in)    :: cat_progn
    type(cat_diagF_type), dimension(N_catd), intent(in)    :: cat_diagF
    type(met_force_type), dimension(N_catd), intent(in)    :: met_force
    
    type(bal_diagn_type), dimension(N_catd), intent(inout) :: bal_diagn
    
    ! local variables
    
    real, dimension(N_catd) :: wtotl, etotl

    integer :: n
    
    real    :: wflx, eflx
    
    character(len=*), parameter :: Iam = 'balance_calcs'
    character(len=400) :: err_msg

    ! -----------------------------------------------------------------------
    
    ! compute total water and energy stores

    call catch_calc_wtotl( N_catd, cat_param%cdcr2, cat_param%wpwet,           &
         cat_progn%srfexc, cat_progn%rzexc, cat_progn%catdef, cat_progn%capac, &
         catprogn2wesn(N_catd,cat_progn), wtotl )
    
    call catch_calc_etotl( N_catd,          cat_param%vegcls, cat_param%dzsf,  &
         cat_param%vgwmax, cat_param%cdcr1, cat_param%cdcr2,                   &
         cat_param%psis,   cat_param%bee,   cat_param%poros,  cat_param%wpwet, &
         cat_param%ars1,   cat_param%ars2,  cat_param%ars3,                    &
         cat_param%ara1,   cat_param%ara2,  cat_param%ara3,   cat_param%ara4,  &
         cat_param%arw1,   cat_param%arw2,  cat_param%arw3,   cat_param%arw4,  &
         cat_progn%srfexc, cat_progn%rzexc, cat_progn%catdef,                  &
         cat_progn%tc1,    cat_progn%tc2,   cat_progn%tc4,                     &
         catprogn2wesn(  N_catd,cat_progn),                                    &
         catprogn2htsn(  N_catd,cat_progn),                                    &
         catprogn2ghtcnt(N_catd,cat_progn),                                    &
         etotl )
        
    ! ------------------------------------------------------------------
    
    select case (option)
       
    case ('ini')                          ! before model propagation
       
       ! initialize water and energy balances 
       
       bal_diagn%wtotl = wtotl
       
       bal_diagn%etotl = etotl
       
       ! set fluxes and increments (in flux units) to no-data-value
       
       bal_diagn%wchng = nodata_generic
       bal_diagn%wincr = nodata_generic
       bal_diagn%echng = nodata_generic
       bal_diagn%eincr = nodata_generic
       
    case ('fin')                            ! after model propagation
       
       do n=1,N_catd
          
          ! water balance
          
          bal_diagn(n)%wchng = (wtotl(n) - bal_diagn(n)%wtotl)/dtstep
          
          bal_diagn(n)%wtotl =  wtotl(n)
          
          wflx =                                                         &
               + met_force(n)%Rainf                                      &
               + met_force(n)%Snowf                                      &
               - cat_diagF(n)%evap                                       &
               - cat_diagF(n)%runoff 
          
          bal_diagn(n)%wincr = bal_diagn(n)%wchng - wflx
          
          ! ------------------------------------------------------------
          
          bal_diagn(n)%echng = (etotl(n) - bal_diagn(n)%etotl)/dtstep
          
          bal_diagn(n)%etotl =  etotl(n)
          
          eflx =                                                         &
               + met_force(n)%SWdown                                     &
               - cat_diagF(n)%swup                                       &
               + met_force(n)%LWdown                                     &
               - cat_diagF(n)%lwup                                       &
               - cat_diagF(n)%lhflux                                     &
               - cat_diagF(n)%shflux                                     &
	       - cat_diagF(n)%hsnacc                                     &
               + cat_diagF(n)%eacc_0                                     &
	       - met_force(n)%Snowf*ALHM    
          
          bal_diagn(n)%eincr = bal_diagn(n)%echng - eflx
          
       end do
       
    case default
       
       err_msg = 'Invalid option for mass balance calculations'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
  end subroutine balance_calcs
  
  ! *************************************************************
  
  subroutine l2f_real( N_l, N_f, N_l_vec, low_ind, data_l, data_f)
    
    ! wrapper for MPI_GATHERV applied to MPI_REAL vector

    ! reichle, 23 Dec 2011
    
    implicit none
  
    integer,                      intent(in)  :: N_l, N_f
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    real,    dimension(N_l),      intent(in)  :: data_l
    real,    dimension(N_f),      intent(out) :: data_f
    
#ifdef LDAS_MPI
    call MPI_GATHERV( &
         data_l, N_l,                MPI_REAL,          &
         data_f, N_l_vec, low_ind-1, MPI_REAL,          &
         0, mpicomm, mpierr )

    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_f = data_l
#endif              
    
  end subroutine l2f_real

  ! *************************************************************
  
  subroutine l2f_real8( N_l, N_f, N_l_vec, low_ind, data_l, data_f)
    
    ! wrapper for MPI_GATHERV applied to MPI_REAL8 vector
    
    ! reichle, 31 Jan 2014
    
    implicit none
    
    integer,                      intent(in)  :: N_l, N_f
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    real*8,  dimension(N_l),      intent(in)  :: data_l
    real*8,  dimension(N_f),      intent(out) :: data_f
    
#ifdef LDAS_MPI
    call MPI_GATHERV( &
         data_l, N_l,                MPI_REAL8,         &
         data_f, N_l_vec, low_ind-1, MPI_REAL8,         &
         0, mpicomm, mpierr )
    
    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_f = data_l
#endif              
    
  end subroutine l2f_real8
  
  ! *************************************************************
  
  subroutine f2l_real( N_f, N_l, N_l_vec, low_ind, data_f, data_l)
    
    ! wrapper for MPI_SCATTERV applied to MPI_REAL vector

    ! reichle, 25 Jul 2013
    
    implicit none
  
    integer,                      intent(in)  :: N_f, N_l
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    real,    dimension(N_f),      intent(in)  :: data_f
    real,    dimension(N_l),      intent(out) :: data_l
    
#ifdef LDAS_MPI
    call MPI_SCATTERV( &
         data_f, N_l_vec, low_ind-1, MPI_REAL,          &
         data_l, N_l,                MPI_REAL,          &
         0, mpicomm, mpierr )

    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_l = data_f
#endif              
    
  end subroutine f2l_real
  
  ! *************************************************************
  
  subroutine f2l_real8( N_f, N_l, N_l_vec, low_ind, data_f, data_l)
    
    ! wrapper for MPI_SCATTERV applied to MPI_REAL8 vector

    ! reichle, 31 Jan 2014
    
    implicit none
  
    integer,                      intent(in)  :: N_f, N_l
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    real*8,  dimension(N_f),      intent(in)  :: data_f
    real*8,  dimension(N_l),      intent(out) :: data_l
    
#ifdef LDAS_MPI
    call MPI_SCATTERV( &
         data_f, N_l_vec, low_ind-1, MPI_REAL8,          &
         data_l, N_l,                MPI_REAL8,          &
         0, mpicomm, mpierr )

    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_l = data_f
#endif              
    
  end subroutine f2l_real8
  
  ! *************************************************************
  
  subroutine f2l_logical( N_f, N_l, N_l_vec, low_ind, data_f, data_l)
    
    ! wrapper for MPI_SCATTERV applied to MPI_LOGICAL vector

    ! reichle,  6 Jun 2016
    
    implicit none
  
    integer,                      intent(in)  :: N_f, N_l
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    logical, dimension(N_f),      intent(in)  :: data_f
    logical, dimension(N_l),      intent(out) :: data_l
    
#ifdef LDAS_MPI
    call MPI_SCATTERV( &
         data_f, N_l_vec, low_ind-1, MPI_LOGICAL,        &
         data_l, N_l,                MPI_LOGICAL,        &
         0, mpicomm, mpierr )

    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_l = data_f
#endif              
    
  end subroutine f2l_logical
  
  ! *************************************************************

end module CLSM_ensdrv_drv_routines

! *********** EOF **************************************************
