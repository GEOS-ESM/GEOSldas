#include "MAPL_Generic.h"

module LDAS_ForceMod

  ! collection of *forcing* subroutines for enkf_driver
  ! (originally these routines were in clsm_ensdrv_drv_routines.F90)

  ! reichle, 13 Aug 2008
  use, intrinsic :: iso_c_binding
  use ESMF
  use MAPL_Mod
  use MAPL_ShmemMod

  use LDAS_ensdrv_Globals,              ONLY:     &
       logunit,                                   &
       logit,                                     &
       root_logit,                                &
       nodata_generic,                            &
       nodata_tol_generic,                        &
       nodata_tolfrac_generic

  use MAPL_ConstantsMod,                ONLY:     &
       stefan_boltzmann => MAPL_STFBOL,           &
       Tzero => MAPL_TICE

  use MAPL_SatVaporMod,                 ONLY:     &
       MAPL_EQsat

  use LDAS_DriverTypes,                 ONLY:     &
       met_force_type

  use LDAS_TileCoordType,               ONLY:     &
       tile_coord_type

  use LDAS_DateTimeMod,                 ONLY:     &
       date_time_type,                            &
       augment_date_time,                         &
       datetime_lt_refdatetime,                   &
       datetime_le_refdatetime,                   &
       is_leap_year,                              &
       get_dofyr_pentad

  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use RepairForcingMod,                 ONLY:     &
       repair_forcing
  use LDAS_HashTable,                   only:     &
       Hash_Table
  implicit none

  include 'netcdf.inc'

  ! everything is private by default unless made public

  private

  public :: get_forcing
  public :: LDAS_move_new_force_to_old
  public :: GEOS_closefile

  type(Hash_Table), public :: FileOpenedHash

  real, parameter :: DEFAULT_REFH   = 10.   ! m
  
  ! Jun 2021: Revised length of unformatted string from 10 to 100 (must be >= 12 for integer).
  character(100), private :: tmpstring100   

  real, contiguous, pointer :: ptrShForce(:,:)=>null()

  type local_grid
     integer :: N_lon = 0
     integer :: N_lat = 0
     integer :: N_cat = 0
     integer, allocatable :: i1(:),i2(:),j1(:),j2(:)
     real,    allocatable :: x1(:),x2(:),y1(:),y2(:)
  end type local_grid

  type(local_grid), target :: local_info

  ! for cubed sphere forcing checking, initialized by GEOS_MetforceGridComp
  integer, public :: im_world_cs = 0

contains

  ! ********************************************************************

  subroutine get_forcing( date_time, force_dtstep, met_path, met_tag,        &
       N_catd, tile_coord, MET_HINTERP, AEROSOL_DEPOSITION,                  &
       MERRA_file_specs, bkwd_looking_fluxes, met_force_obs_tile_new,        & 
       init )
    
    ! Read and check meteorological forcing data for the domain.
    !
    ! forcing readers must provide ALL of the surface meteorological fields in the
    !   met_force_type structure EXCEPT:
    !   - PARdrct, PARdffs : if not available, will be backfilled as fraction of SWdown
    !   - aerosol forcing  : currently defunct
    !
    ! time convention:
    ! - forcing states (such as Tair) are snapshots at date_time
    ! - forcing fluxes (such as SWdn) are time avg over *subsequent* (*forward-looking*) 
    !    forcing interval (date_time:date_time+force_dtstep)    
    !
    ! The above time convention was inherited from older versions of the 
    ! off-line driver and creates problems with "operational" forcing
    ! data from GEOS.  For "operational" integrations, the forward-looking
    ! forcing fluxes are not available for "met_force_obs_tile_new".  
    !
    ! For datasets that provide fluxes over the forcing interval that *precedes*
    ! date_time_new, the output parameter "bkwd_looking_fluxes" must be set to .true.,
    ! so that subroutine LDAS_move_new_force_to_old() can time-shift the fields
    ! accordingly.
    !
    ! When LDASsa is integrated within the coupled GEOS5 DAS, initial (time-avg)
    ! "tavg1_2d_*_Nx" files are not available.  Use optional "init" flag to 
    ! deal with this situation.
    !
    ! reichle,      28 March 2006
    ! reichle,      13 March 2008 - added optional "init" flag
    ! qliu+reichle, 12 Aug   2008 - new field RefH (reference height) in met_force_type
    ! reichle,      25 Sep   2009 - removed unneeded inputs 
    ! reichle,      23 Feb   2016 - new and more efficient work-around to make GEOS-5 
    !                                forcing work with LDASsa time convention for forcing data
    ! borescan,     01 Feb   2021 - added ERA5_LIS forcing 
    ! reichle,      12 Apr   2021 - removed obsolete optional input "alb_from_SWnet"
    !                             - replaced "GEOS_forcing" switch with "bkwd_looking_fluxes"
    !                             - added checks for supported options of MET_HINTERP and
    !                                AEROSOL_DEPOSITION
    ! reichle,      22 Apr   2021 - clean up:
    !                               - calls to check_forcing_nodata() and repair_forcing()
    !                               - handling nodata-values in PARdrct, PARdffs
    ! reichle,      23 Apr   2021 - clean up:
    !                               - removed SWnet from met_force_type
    
    implicit none

    ! intent in:
    
    type(date_time_type), intent(in) :: date_time
    
    integer,              intent(in) :: force_dtstep

    character(*),         intent(in) :: met_path
    character(*),         intent(in) :: met_tag

    integer,              intent(in) :: N_catd, MET_HINTERP, AEROSOL_DEPOSITION

    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input

    ! intent out:
    
    logical,              intent(out) :: MERRA_file_specs
    logical,              intent(out) :: bkwd_looking_fluxes

    type(met_force_type), dimension(N_catd), intent(out) :: &
         met_force_obs_tile_new

    ! optional:
    
    logical, intent(in), optional  :: init

    ! -------------------
    
    ! local variables
    
    real    :: nodata_forcing, tol

    logical :: PAR_available                         ! indicate whether reader provides PARdrct, PARdffs
    
    logical :: supported_option_MET_HINTERP          ! for consistency check of resource parameter settings
    logical :: supported_option_AEROSOL_DEPOSITION   ! for consistency check of resource parameter settings

    logical :: unlimited_Qair, unlimited_LWdown      ! options for repair_forcing()
    
    type(date_time_type) :: date_time_tmp
    
    character(len=*), parameter :: Iam = 'get_forcing'
    character(len=400) :: err_msg
    
    ! --------------------------------------------------------------
    !
    ! shift forcing date if so indicated by met_tag (for twin experiments, 
    ! see function shift_forcing_date for details) - reichle, 6 Apr 2007
    
    date_time_tmp = shift_forcing_date(met_tag, date_time)
        
    ! set reference height to default value (if appropriate, will be overwriten
    !  within specific subroutine)
    
    met_force_obs_tile_new%RefH  = DEFAULT_REFH

    ! Note that "nodata_forcing" is set to the native nodata-value
    !  in the individual get_*() subroutines and used to communicate with
    !  check_forcing_nodata.  AFTER the call to check_forcing_nodata all forcing 
    !  fields must NOT be nodata values.
    !
    ! reichle+qliu,  8 Oct 2008
    ! reichle, 23 Feb 2009 -- same goes for ParDrct, ParDffs
    ! reichle,  5 Mar 2009 -- deleted ParDrct, ParDffs after testing found no impact
    ! reichle, 22 Jul 2010 -- fixed treatment of SWnet nodata values
    ! reichle, 20 Dec 2011 -- reinstated PARdrct and PARdffs for MERRA-Land file specs
    
    ! ---------------------------------------------------------------------------------
    !
    ! initialize 
    
    MERRA_file_specs                    = .false.

    bkwd_looking_fluxes                 = .false.

    PAR_available                       = .false.  ! default; so far, only GEOS forcing provides PAR
    
    unlimited_Qair                      = .false.  ! default for call to repair_forcing
    unlimited_LWdown                    = .false.  ! default for call to repair_forcing
    
    ! every forcing data reader must support the default settings:
    ! 
    !   MET_HINTER         = 0 : nearest neighbor
    !   AEROSOL_DEPOSITION = 0 : *no* aerosol deposition
    !
    ! initialize "supported_option_*" to .true. for defaults and .false. otherwise; then
    !  set to .true. in individual forcing readers for any non-default options that are
    !  supported for the selected met_tag 
    
    supported_option_MET_HINTERP        = (MET_HINTERP        == 0)

    supported_option_AEROSOL_DEPOSITION = (AEROSOL_DEPOSITION == 0)

    ! ---------------------------------------------------------------------------------
    !
    ! get forcing in tile space

    if     (index(met_tag, 'Berg_netcdf')/=0) then
       
       call get_Berg_netcdf(       date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
    elseif (index(met_tag, 'GLDAS_2x2_5_netcdf')/=0) then
       
       call get_GLDAS_2x2_5_netcdf(date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
    elseif (index(met_tag, 'Viviana_OK')/=0) then
       
       ! vmaggion & reichle, 17 July 2008
       !
       ! use 2x2.5 deg GLDAS for all forcing fields except precip
       
       call get_GLDAS_2x2_5_netcdf(date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
       if (index(met_tag, 'Viviana_OK_nopert')/=0) then
          
          call get_Viviana_OK_precip(10, date_time_tmp, met_path, met_tag,  &
               N_catd, tile_coord, met_force_obs_tile_new)
          
       end if

    elseif (index(met_tag, 'GSWP2_1x1_netcdf')/=0) then
       
       call get_GSWP2_1x1_netcdf(  date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)

    elseif (index(met_tag, 'RedArk_ASCII')/=0) then
       
       call get_RedArk_ASCII(      date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
    elseif (index(met_tag, 'RedArk_GOLD')/=0) then
       
       call get_RedArk_GOLD(       date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
    elseif (index(met_tag, 'RedArk_Princeton')/=0) then
       
       call get_RedArk_Princeton(  date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
    elseif (index(met_tag, 'Princeton_netcdf')/=0) then ! tyamada+reichle, 17 Jul 2007   
       
       call get_Princeton_netcdf(  date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
    elseif (index(met_tag, 'conus_0.5d_netcdf')/=0) then ! sarith+reichle, 17 Jul 2007   
       
       call get_conus_netcdf(      date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
    elseif (index(met_tag, 'ERA5_LIS')/=0) then 
       
       call get_ERA5_LIS(          date_time_tmp, met_path, N_catd, tile_coord, &
            met_force_obs_tile_new, nodata_forcing)
       
       ! Subroutine get_ERA5_LIS() provided backward-looking fluxes.
       ! The time convention stated above is restored through a later call to subroutine
       ! LDAS_move_new_force_to_old().
       
       bkwd_looking_fluxes            = .true.        
       
       ! model-based dataset; call repair_forcing() below without certain limitations
       
       unlimited_Qair                 = .true.
       unlimited_LWdown               = .true.

    elseif (index(met_tag(1:7), 'GEOSs2s')/=0) then

       call get_GEOSs2s( date_time_tmp, met_path, met_tag, N_catd, tile_coord, &
            MET_HINTERP, met_force_obs_tile_new, nodata_forcing, PAR_available)

    else ! assume forcing from GEOS5 GCM ("DAS" or "MERRA") output
       
       if(root_logit) write (logunit,*) 'get_forcing(): assuming GEOS-5 forcing data set'
       
       call get_GEOS( date_time_tmp, force_dtstep, met_path, met_tag,                 &
            N_catd, tile_coord, MET_HINTERP, AEROSOL_DEPOSITION,                      &
            supported_option_MET_HINTERP,                                             &
            supported_option_AEROSOL_DEPOSITION,                                      &
            met_force_obs_tile_new, nodata_forcing, PAR_available, MERRA_file_specs,  &
            init )

       ! subroutine get_GEOS() provided backward-looking fluxes. 
       ! The time convention is restored through a later call to subroutine
       ! LDAS_move_new_force_to_old().
       !
       ! Subroutine get_GEOS() reads forcing fluxes from "previous"
       ! interval, not from "subsequent" interval, because in operational
       ! applications the "subsequent" fluxes for "met_force_new" are not
       ! available.  Note that only "old" fluxes are needed in time interpolation
       ! (module LDAS_InterpMod).
       
       bkwd_looking_fluxes            = .true.
       
       ! model-based dataset; call repair_forcing() below without certain limitations
       !
       ! call repair_forcing with switch "unlimited_Qair=.true." 
       ! (default is to limit Qair so that it does not exceed Qair_sat)
       ! reichle+qliu,  8 Oct 2008
       ! 
       ! likewise for "unlimited_LWdown=.true."
       ! reichle, 11 Feb 2009       
       
       unlimited_Qair                 = .true.
       unlimited_LWdown               = .true.

    end if

    ! ------------------
    !
    ! backfill PAR:
    !
    ! PAR might be missing for one of the following reasons:
    ! - the reader did not provide PAR
    ! - the reader provides PAR but owing to a land mask difference between
    !     the forcing data and the LDAS simulation, some tiles may not have PAR 
    ! - the reader provides PAR but there are gaps in the original dataset
    !
    ! if PAR was not available at all from the forcing reader, assign "nodata_forcing"

    if (.not. PAR_available) then
       met_force_obs_tile_new%PARdrct = nodata_forcing
       met_force_obs_tile_new%PARdffs = nodata_forcing
    end if

    ! the nodata-value for all forcing fields is now "nodata_forcing"
    
    ! where PAR is not available for any of the above reasons:
    ! - assume half of SWdown is photosynthetically active
    ! - assume half of PAR is direct, half diffuse
    
    tol = abs(nodata_forcing*nodata_tolfrac_generic)
    
    where (                                                                       &
         (abs(met_force_obs_tile_new%PARdrct-nodata_forcing) < tol)  .and.        &
         (abs(met_force_obs_tile_new%SWdown -nodata_forcing) > tol)            )
       met_force_obs_tile_new%PARdrct = 0.5*0.5*met_force_obs_tile_new%SWdown
       met_force_obs_tile_new%PARdffs = met_force_obs_tile_new%PARdrct
    end where
    
    ! ------------------
    !
    ! check for nodata values and fill with neighboring data if sensible;
    !   otherwise, check_forcing_nodata() will abort
    
    call check_forcing_nodata( N_catd, tile_coord, nodata_forcing,                &
         met_force_obs_tile_new )
    
    ! ------------------
    !
    ! reset unphysical or inconsistent forcing values
    
    call repair_forcing( N_catd, met_force_obs_tile_new,                          &
         echo=.true., tile_coord=tile_coord,                                      &
         fieldname='all',                                                         &
         unlimited_Qair=unlimited_Qair, unlimited_LWdown=unlimited_LWdown )
    
    ! ------------------
    !
    ! stop if a supported_option_* switch remains .false.
    
    if (.not. supported_option_MET_HINTERP)        then
       err_msg = 'selected MET_HINTERP option not supported for met_tag '         &
            //  trim(met_tag)
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    if (.not. supported_option_AEROSOL_DEPOSITION) then
       err_msg = 'selected AEROSOL_DEPOSITION option not supported for met_tag '  &
            //  trim(met_tag)
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
  end subroutine get_forcing

   !**************************************************************************************

   subroutine LDAS_move_new_force_to_old( MERRA_file_specs, AEROSOL_DEPOSITION, &
        new_force, old_force )

     ! move *flux*-type forcing data from "new" to "old";
     ! call this subroutine ONLY if the forcing reader provides backward-looking fluxes
     ! - reichle, 14 Apr 2021
     
     implicit none
     
     logical, intent(in) :: MERRA_file_specs
     integer, intent(in) :: AEROSOL_DEPOSITION
     
     type(met_force_type), dimension(:), intent(inout) :: new_force
     type(met_force_type), dimension(:), intent(inout) :: old_force
     
     old_force%Rainf_C = new_force%Rainf_C
     old_force%Rainf   = new_force%Rainf
     old_force%Snowf   = new_force%Snowf
     old_force%LWdown  = new_force%LWdown
     old_force%SWdown  = new_force%SWdown
     old_force%PARdrct = new_force%PARdrct
     old_force%PARdffs = new_force%PARdffs

     new_force%Rainf_C = nodata_generic
     new_force%Rainf   = nodata_generic
     new_force%Snowf   = nodata_generic
     new_force%LWdown  = nodata_generic
     new_force%SWdown  = nodata_generic
     new_force%PARdrct = nodata_generic
     new_force%PARdffs = nodata_generic
     
     ! [moved here from below, reichle, 28 Jan 2021]
     ! treat Wind as flux when forcing with MERRA
     if (MERRA_file_specs) then
        old_force%Wind  = new_force%Wind
        new_force%Wind  = nodata_generic
     endif
     
     ! not sure what exactly the following fields are and
     ! whether it makes sense to include them here
     ! - reichle, 28 Jan 2021
     !
     if( AEROSOL_DEPOSITION /=0) then
        
        old_force%DUDP001 = new_force%DUDP001
        old_force%DUDP002 = new_force%DUDP002
        old_force%DUDP003 = new_force%DUDP003
        old_force%DUDP004 = new_force%DUDP004
        old_force%DUDP005 = new_force%DUDP005
        old_force%DUSV001 = new_force%DUSV001
        old_force%DUSV002 = new_force%DUSV002
        old_force%DUSV003 = new_force%DUSV003
        old_force%DUSV004 = new_force%DUSV004
        old_force%DUSV005 = new_force%DUSV005
        old_force%DUWT001 = new_force%DUWT001
        old_force%DUWT002 = new_force%DUWT002
        old_force%DUWT003 = new_force%DUWT003
        old_force%DUWT004 = new_force%DUWT004
        old_force%DUWT005 = new_force%DUWT005
        old_force%DUSD001 = new_force%DUSD001
        old_force%DUSD002 = new_force%DUSD002
        old_force%DUSD003 = new_force%DUSD003
        old_force%DUSD004 = new_force%DUSD004
        old_force%DUSD005 = new_force%DUSD005
        old_force%BCDP001 = new_force%BCDP001
        old_force%BCDP002 = new_force%BCDP002
        old_force%BCSV001 = new_force%BCSV001
        old_force%BCSV002 = new_force%BCSV002
        old_force%BCWT001 = new_force%BCWT001
        old_force%BCWT002 = new_force%BCWT002
        old_force%BCSD001 = new_force%BCSD001
        old_force%BCSD002 = new_force%BCSD002
        old_force%OCDP001 = new_force%OCDP001
        old_force%OCDP002 = new_force%OCDP002
        old_force%OCSV001 = new_force%OCSV001
        old_force%OCSV002 = new_force%OCSV002
        old_force%OCWT001 = new_force%OCWT001
        old_force%OCWT002 = new_force%OCWT002
        old_force%OCSD001 = new_force%OCSD001
        old_force%OCSD002 = new_force%OCSD002
        old_force%SUDP003 = new_force%SUDP003
        old_force%SUSV003 = new_force%SUSV003
        old_force%SUWT003 = new_force%SUWT003
        old_force%SUSD003 = new_force%SUSD003
        old_force%SSDP001 = new_force%SSDP001
        old_force%SSDP002 = new_force%SSDP002
        old_force%SSDP003 = new_force%SSDP003
        old_force%SSDP004 = new_force%SSDP004
        old_force%SSDP005 = new_force%SSDP005
        old_force%SSSV001 = new_force%SSSV001
        old_force%SSSV002 = new_force%SSSV002
        old_force%SSSV003 = new_force%SSSV003
        old_force%SSSV004 = new_force%SSSV004
        old_force%SSSV005 = new_force%SSSV005
        old_force%SSWT001 = new_force%SSWT001
        old_force%SSWT002 = new_force%SSWT002
        old_force%SSWT003 = new_force%SSWT003
        old_force%SSWT004 = new_force%SSWT004
        old_force%SSWT005 = new_force%SSWT005
        old_force%SSSD001 = new_force%SSSD001
        old_force%SSSD002 = new_force%SSSD002
        old_force%SSSD003 = new_force%SSSD003
        old_force%SSSD004 = new_force%SSSD004
        old_force%SSSD005 = new_force%SSSD005


        new_force%DUDP001 = nodata_generic
        new_force%DUDP002 = nodata_generic
        new_force%DUDP003 = nodata_generic
        new_force%DUDP004 = nodata_generic
        new_force%DUDP005 = nodata_generic
        new_force%DUSV001 = nodata_generic
        new_force%DUSV002 = nodata_generic
        new_force%DUSV003 = nodata_generic
        new_force%DUSV004 = nodata_generic
        new_force%DUSV005 = nodata_generic
        new_force%DUWT001 = nodata_generic
        new_force%DUWT002 = nodata_generic
        new_force%DUWT003 = nodata_generic
        new_force%DUWT004 = nodata_generic
        new_force%DUWT005 = nodata_generic
        new_force%DUSD001 = nodata_generic
        new_force%DUSD002 = nodata_generic
        new_force%DUSD003 = nodata_generic
        new_force%DUSD004 = nodata_generic
        new_force%DUSD005 = nodata_generic
        new_force%BCDP001 = nodata_generic
        new_force%BCDP002 = nodata_generic
        new_force%BCSV001 = nodata_generic
        new_force%BCSV002 = nodata_generic
        new_force%BCWT001 = nodata_generic
        new_force%BCWT002 = nodata_generic
        new_force%BCSD001 = nodata_generic
        new_force%BCSD002 = nodata_generic
        new_force%OCDP001 = nodata_generic
        new_force%OCDP002 = nodata_generic
        new_force%OCSV001 = nodata_generic
        new_force%OCSV002 = nodata_generic
        new_force%OCWT001 = nodata_generic
        new_force%OCWT002 = nodata_generic
        new_force%OCSD001 = nodata_generic
        new_force%OCSD002 = nodata_generic
        new_force%SUDP003 = nodata_generic
        new_force%SUSV003 = nodata_generic
        new_force%SUWT003 = nodata_generic
        new_force%SUSD003 = nodata_generic
        new_force%SSDP001 = nodata_generic
        new_force%SSDP002 = nodata_generic
        new_force%SSDP003 = nodata_generic
        new_force%SSDP004 = nodata_generic
        new_force%SSDP005 = nodata_generic
        new_force%SSSV001 = nodata_generic
        new_force%SSSV002 = nodata_generic
        new_force%SSSV003 = nodata_generic
        new_force%SSSV004 = nodata_generic
        new_force%SSSV005 = nodata_generic
        new_force%SSWT001 = nodata_generic
        new_force%SSWT002 = nodata_generic
        new_force%SSWT003 = nodata_generic
        new_force%SSWT004 = nodata_generic
        new_force%SSWT005 = nodata_generic
        new_force%SSSD001 = nodata_generic
        new_force%SSSD002 = nodata_generic
        new_force%SSSD003 = nodata_generic
        new_force%SSSD004 = nodata_generic
        new_force%SSSD005 = nodata_generic

     endif ! AEROSOL_DEPOSITION /=0

  end subroutine LDAS_move_new_force_to_old 
  ! ****************************************************************  
  
  subroutine get_Berg_netcdf(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )

    ! read Berg_NetCDF files and extract forcings in tile space
    ! (uses nearest neighbor interpolation)
    
    ! reichle, 25 May 2005
    ! reichle, 23 Feb 2009 -- revised treatment of Rainf_C

    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    character(*),       intent(in) :: met_path
    
    integer,              intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new
    
    real, intent(out) :: nodata_forcing
  
    ! Berg grid and netcdf parameters 
        
    integer, parameter :: berg_grid_N_lon  = 720
    integer, parameter :: berg_grid_N_lat  = 360
    real,    parameter :: berg_grid_ll_lon = -180.
    real,    parameter :: berg_grid_ll_lat = -90.
    real,    parameter :: berg_grid_dlon   = .5
    real,    parameter :: berg_grid_dlat   = .5
    
    integer, parameter :: N_berg_compressed = 67420
    
    ! Berg forcing time step in hours

    integer, parameter :: dt_berg_in_hours = 6 
    
    integer, parameter :: nciv_land_i = 3
    integer, parameter :: nciv_land_j = 4
    integer, parameter :: nciv_data   = 7
    
    integer, parameter :: N_berg_vars = 8
    
    real,    parameter :: nodata_berg = 1.e20
    
    character(40), dimension(N_berg_vars), parameter :: berg_dir = (/ &
         'PREC-CONV/',    &
         'PREC-TOTL/',    &
         'PRES-SRF/ ',    &
         'RAD-LW/   ',    &
         'RAD-SW/   ',    &
         'TEMP-AIR/ ',    &
         'TEMP-DEW/ ',    &
         'WIND/     '       /)
    
    character(40), dimension(N_berg_vars), parameter :: berg_name = (/ &
         'RainfSnowf_C_ecmwf',       &
         'RainfSnowf_ecmwf  ',       &
         'PSurf_ecmwf       ',       &
         'LWdown_ecmwf      ',       &
         'SWdown_ecmwf      ',       &
         'Tair_ecmwf        ',       &
         'Tdew_ecmwf        ',       &
         'Wind_ecmwf        '      /)

    ! local variables
    
    !!real, parameter :: Tzero = 273.16
    
    real :: tol 
        
    real, dimension(berg_grid_N_lon,berg_grid_N_lat) :: tmp_grid
    
    integer, dimension(N_berg_compressed)   :: land_i_berg, land_j_berg
    integer, dimension(N_catd)              :: i_ind, j_ind
    
    real,    dimension(N_berg_compressed)   :: tmp_vec
    
    real,    dimension(N_catd,N_berg_vars) :: force_array
    
    integer, dimension(2) :: start, icount
    
    integer :: k, n, hours_in_month, berg_var, ierr, ncid
    
    real    :: this_lon, this_lat, dt_berg_in_seconds
    
    character(4) :: YYYY
    character(2) :: MM
    
    character(300) :: fname

    character(len=*), parameter :: Iam = 'get_Berg_netcdf'
    character(len=400) :: err_msg
    
    ! --------------------------------------------------------------------
    
    dt_berg_in_seconds = real(3600*dt_berg_in_hours)
    
    nodata_forcing = nodata_berg
    
    tol = abs(nodata_forcing*nodata_tolfrac_generic)

    ! assemble year and month strings
    
    write (YYYY,'(i4.4)') date_time%year
    write (MM,  '(i2.2)') date_time%month
    
    ! find out which data are needed
    
    ! compressed space dimension (always read global vector)
    
    start(1)  = 1
    icount(1) = N_berg_compressed
    
    ! time dimension (first entry in Berg_NetCDF file is at 0Z)

    if ( (date_time%min/=0) .or. (date_time%sec/=0) .or.           &
         (mod(date_time%hour,dt_berg_in_hours)/=0)        ) then
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'timing error')
              
    end if
    
    hours_in_month = (date_time%day-1)*24 + date_time%hour
    
    start(2)  = hours_in_month / dt_berg_in_hours + 1 
    icount(2) = 1
    
    ! ----------------------------------------------    
    !
    ! compute indices for nearest neighbor interpolation from Berg grid 
    ! to tile space
    !
    ! (NOTE: this should at some point be replaced with a regridding
    !  subroutine that interpolates from the
    !  native forcing grid to the GCM atmospheric grid that is used
    !  to cut catchments into tiles - then "standard" grid2tile
    !  using tile_coord%atm_i and tile_coord%atm_j applies. 
    !  reichle, 26 May 2005)
    
    do k=1,N_catd
       
       ! ll_lon and ll_lat refer to lower left corner of grid cell
       ! (as opposed to the grid point in the center of the grid cell)
       
       this_lon = tile_coord(k)%com_lon
       this_lat = tile_coord(k)%com_lat
       
       i_ind(k) = ceiling( (this_lon - berg_grid_ll_lon)/berg_grid_dlon )
       j_ind(k) = ceiling( (this_lat - berg_grid_ll_lat)/berg_grid_dlat )
       
       ! NOTE: For a "date line on center" grid and (180-dlon/2) < lon < 180 
       ! we now have i_ind=(grid%N_lon+1) 
       ! This would need to be fixed.

       if (i_ind(k)>berg_grid_N_lon)  i_ind(k)=1
              
    end do
    
    ! ------------------------------------------------------
    !
    ! read compression parameters (same for all data variables and time steps)
    
    berg_var = 1
    
    fname = trim(met_path) // trim(berg_dir(berg_var)) // '/' // YYYY       &
         // '/' // trim(berg_name(berg_var)) // '.' // YYYY // MM // '.nc'
        
    if(root_logit) write(logunit,*) 'get netcdf compression params from ' // trim(fname)
    
    ierr = NF_OPEN(fname,NF_NOWRITE,ncid)
    
    if (ierr/=0) then
       err_msg = 'error opening netcdf file'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    ierr = NF_GET_VARA_INT(ncid, nciv_land_i, start, icount, land_i_berg)
    ierr = NF_GET_VARA_INT(ncid, nciv_land_j, start, icount, land_j_berg)
    
    ierr = NF_CLOSE(ncid)

    ! ------------------------------------------------------
    !
    ! get forcing data
    
    do berg_var = 1,N_berg_vars
       
       ! open file, read compressed data, and put on global grid
       
       fname = trim(met_path) // trim(berg_dir(berg_var)) // '/' // YYYY     &
            // '/' // trim(berg_name(berg_var)) // '.' // YYYY // MM // '.nc'
       
       if(root_logit) write (logunit,*) 'opening ' // trim(fname)
       
       ierr = NF_OPEN(fname,NF_NOWRITE,ncid)
       
       if (ierr/=0)  then
          err_msg = 'error opening netcdf file'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       ierr = NF_GET_VARA_REAL(ncid, nciv_data, start, icount, tmp_vec )
       
       ierr = NF_CLOSE(ncid)
       
       tmp_grid = nodata_forcing

       do n=1,N_berg_compressed
          
          tmp_grid(land_i_berg(n), land_j_berg(n) ) = tmp_vec(n)
          
       end do
       
       ! interpolate to tile space

       ! (NOTE: This should at some point be replaced with a regridding
       !  subroutine that interpolates from the
       !  native forcing grid to the GCM atmospheric grid that is used
       !  to cut catchments into tiles - then "standard" grid2tile
       !  using tile_coord%atm_i and tile_coord%atm_j applies. 
       !  reichle, 26 May 2005)
              
       do k=1,N_catd
          
          force_array(k,berg_var) = tmp_grid(i_ind(k), j_ind(k))
          
       end do
       
    end do
    
    ! convert variables and units of force_array to match met_force_type, 
    ! put into structure
    
    ! from Berg files:
    !
    !  force_array(:,1) = PREC-CONV = Rainf_C+Snowf_C   kg/m2 (6h total)  ??? 
    !  force_array(:,2) = PREC-TOTL = RainfSnowf        kg/m2 (6h total)	
    !  force_array(:,3) = PRES-SRF  = PSurf             Pa			
    !  force_array(:,4) = RAD-LW    = LWdown            W/m2			
    !  force_array(:,5) = RAD-SW    = SWdown            W/m2			
    !  force_array(:,6) = TEMP-AIR  = Tair              K			
    !  force_array(:,7) = TEMP-DEW  = Tdew              K			
    !  force_array(:,8) = WIND      = Wind              m/s
    
    met_force_new%Psurf   = force_array(:,3)
    met_force_new%LWdown  = force_array(:,4)
    met_force_new%SWdown  = force_array(:,5)
    met_force_new%Tair    = force_array(:,6)
    met_force_new%Wind    = force_array(:,8)
    
    ! get specific humidity from dew point temperature
    
    do k=1,N_catd
       
       if ( abs(force_array(k,3)-nodata_berg)<tol .or.           &
            abs(force_array(k,7)-nodata_berg)<tol      ) then
          
          met_force_new(k)%Qair = nodata_forcing
          
       else
          
          met_force_new(k)%Qair = MAPL_EQsat(force_array(k,7),force_array(k,3))
          
       end if
       
    end do
    
    ! rain and snow:
    ! convert from Berg 6h totals [kg/m2] (or [mm]) into precipitation 
    ! rates [kg/m2/s] 
    
    met_force_new%Rainf   = 0.
    met_force_new%Rainf_C = 0.
    met_force_new%Snowf   = 0.
    
    do k=1,N_catd
       
       ! set convective precip to zero for nodata-values
       
       if (abs(force_array(k,1)-nodata_berg)<tol) force_array(k,1) = 0.
       
       if (abs(force_array(k,2)-nodata_berg)<tol) then
          
          met_force_new(k)%Snowf   = nodata_forcing
          met_force_new(k)%Rainf   = nodata_forcing
          met_force_new(k)%Rainf_C = nodata_forcing
          
       else
          
          if ( met_force_new(k)%Tair < Tzero ) then
             met_force_new(k)%Snowf   = force_array(k,2)/dt_berg_in_seconds
          else
             met_force_new(k)%Rainf   = force_array(k,2)/dt_berg_in_seconds
             ! leave Rainf_C=0 until sure that "PREC-CONV = Rainf_C+Snowf_C"
             !!met_force_new(k)%Rainf_C = force_array(k,1)/dt_berg_in_seconds
          end if
          
       end if


    end do
    
  end subroutine get_Berg_netcdf

  ! ****************************************************************  
  
  subroutine get_RedArk_ASCII(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )
    
    ! read RedArk ASCII forcing files from Hatim Sharif and Wade Crow
    ! extract forcings in tile space (use nearest neighbor interpolation)
    
    ! reichle, 13 Apr 2006
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    character(*),       intent(in) :: met_path
    
    integer,              intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new
    
    real, intent(out) :: nodata_forcing
    
    ! RedArk parameters 
    
    integer, parameter :: N_redark_basins  = 314
    
    character(40), parameter :: fname_b2t = 'RedArkBasins_2_FV288x181Tiles.dat'
    
    ! RedArk forcing time step in hours
    
    integer, parameter :: dt_redark_in_hours = 1 
    
    integer, parameter :: N_redark_vars = 6
    
    ! not checking for "bad" forcing data from Hatim/Wade, but need
    ! to provide output "nodata_forcing"
    
    real,    parameter :: nodata_redark = 1.e20
    
    real,    parameter :: RedArk_Psurf_in_Pa = 101300.
    
    ! local variables

    type(date_time_type) :: date_time_CST
    
    !!real, parameter :: Tzero = 273.16
    
    integer, dimension(N_catd)                        :: b_ind
    
    real,    dimension(N_redark_basins,N_redark_vars) :: force_array
    
    integer :: i, j, k, tmp_tile_id
    
    real    :: dt_redark_in_seconds, Tdew
    
    character(4) :: YYYY
    character(3) :: DDD
    character(2) :: HH

    character(300) :: fname
    
    character(len=*), parameter :: Iam = 'get_RedArk_ASCII'
    character(len=400) :: err_msg

    ! --------------------------------------------------------------------
    
    nodata_forcing = nodata_redark
    
    dt_redark_in_seconds = real(3600*dt_redark_in_hours)
    
    ! convert GMT (date_time) into CST (Central Standard Time)
    
    date_time_CST = date_time
    
    call augment_date_time( -6*3600, date_time_CST )
    
    ! assemble year, day-of-year, and hour strings
    
    write (YYYY,'(i4.4)') date_time_CST%year
    write (DDD, '(i3.3)') date_time_CST%dofyr
    write (HH,  '(i2.2)') date_time_CST%hour

    ! ----------------------------------------------    
    !
    ! read indices for nearest neighbor interpolation from RedArk basins 
    ! to tile space
    
    fname = trim(met_path) // '/' // fname_b2t
    
    open(10, file=fname, form='formatted', action='read', status='old')
    
    do i=1,N_catd
       
       read (10,*) tmp_tile_id, b_ind(i)
       
       if (tmp_tile_id /= tile_coord(i)%tile_id) then
          err_msg = 'RedArk basin2tile error'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

    end do
    
    close(10, status='keep')
    
    ! ------------------------------------------------------
    !
    ! read forcing data
    
    fname = trim(met_path) // '/' // YYYY // '/'    &
         // 'red_ark_forc' // '.' // YYYY // '.' // DDD // '.' // HH
    
    if(root_logit) write(logunit,*) 'opening ' // trim(fname)
    
    open(10, file=fname, form='formatted', action='read', status='old')
    
    do i=1,N_redark_basins
       
       read (10,*) (force_array(i,j), j=1,N_redark_vars)
       
    end do
    
    close(10, status='keep')
     
    ! ------------------------------------------------------
    !
    ! convert variables and units of force_array to match met_force_type, 
    ! map from basin to tile, put into structure
    
    ! from RedArk files:
    !
    !  force_array(:,1) = RainfSnowf   kg/m2 (1h total)	
    !  force_array(:,2) = SWdown       W/m2			
    !  force_array(:,3) = LWdown       W/m2			
    !  force_array(:,4) = Tair         degC			
    !  force_array(:,5) = Wind         m/s
    !  force_array(:,6) = Tdew         degC		

    do k=1,N_catd
       
       met_force_new(k)%Psurf   = RedArk_Psurf_in_Pa
       
       met_force_new(k)%SWdown  = force_array(b_ind(k),2)
       met_force_new(k)%LWdown  = force_array(b_ind(k),3)
       met_force_new(k)%Tair    = force_array(b_ind(k),4) + Tzero
       met_force_new(k)%Wind    = force_array(b_ind(k),5)
       
       ! get specific humidity from dew point temperature
       
       Tdew          = force_array(b_ind(k),6) + Tzero
       
       met_force_new(k)%Qair = MAPL_EQsat(Tdew,RedArk_Psurf_in_Pa)
       
       ! rain and snow:
       ! convert from RedArk 1h totals [kg/m2] (or [mm]) into precipitation 
       ! rates [kg/m2/s] 
       ! partition total precip into rain and snow according to Tair
       ! set convective precip to zero
       
       met_force_new(k)%Rainf_C = 0.
       
       if ( met_force_new(k)%Tair < Tzero ) then
          met_force_new(k)%Rainf = 0.
          met_force_new(k)%Snowf = force_array(b_ind(k),1)/dt_redark_in_seconds
       else
          met_force_new(k)%Rainf = force_array(b_ind(k),1)/dt_redark_in_seconds
          met_force_new(k)%Snowf = 0.
       end if
       
    end do
    
  end subroutine get_RedArk_ASCII
  
  ! ****************************************************************  
  
  subroutine get_RedArk_GOLD(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )
    
    ! read RedArk GOLD files and map to forcings in tile space
    ! (uses nearest neighbor interpolation)
    
    ! reichle, 9 Oct 2006
    ! reichle, 23 Feb 2009 - revised treatment of Rainf_C
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    character(*),       intent(in) :: met_path
    
    integer,              intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new
    
    real, intent(out) :: nodata_forcing
    
    ! interpolation parameters
    
    integer, parameter :: RedArk_GOLD_N_cells = 28
    
    !character(200), parameter :: RedArk_GOLD_coord_path = &
    !'/land/l_data/RedArk/coords/'
    
    character(80),  parameter :: RedArk_GOLD_coord_file = &
         'GOLD_forcing_cells_2_RedArkOSSE_tiles.dat'    
    
    ! GOLD_REDARK forcing time step in hours
    
    integer, parameter :: dt_RedArk_GOLD_in_hours = 6
    
    integer, parameter :: N_RedArk_GOLD_vars = 9
    
    real,    parameter :: nodata_RedArk_GOLD = 1.e20
    
    character(40), dimension(N_RedArk_GOLD_vars), parameter :: &
         RedArk_GOLD_name = (/ &
         'SWdown_era     ',   &   !  1 - flux
         'LWdown_era     ',   &   !  2 - flux
         'Rainf_C_gold   ',   &   !  3 - flux       ??? 
         'Rainf_gold     ',   &   !  4 - flux
         'Snowf_gold     ',   &   !  5 - flux
         'PSurf_era      ',   &   !  6 - state
         'Qair_era       ',   &   !  7 - state
         'Tair_era       ',   &   !  8 - state
         'Wind_era       ' /)     !  9 - state
    
    ! local variables
    
    real :: tol 
    
    integer, dimension(N_catd)               :: ind_gold2tile
    
    real,    dimension(RedArk_GOLD_N_cells)  :: tmp_vec

    real,    dimension(N_catd,N_RedArk_GOLD_vars) :: force_array
    
    type(date_time_type) :: date_time_tmp
    
    integer :: int_dt_RedArk_GOLD_in_seconds, i, k, this_var, tmp_tile_id
    
    character(4) :: YYYY, HHMM
    character(2) :: MM, DD
    
    character(300) :: fname
    
    character(len=*), parameter :: Iam = 'get_RedArk_GOLD'
    character(len=400) :: err_msg

    ! --------------------------------------------------------------------
    
    int_dt_RedArk_GOLD_in_seconds = 3600*dt_RedArk_GOLD_in_hours
    
    nodata_forcing = nodata_RedArk_GOLD
    
    tol = abs(nodata_forcing*nodata_tolfrac_generic)
    
    ! ----------------------------------------------    
    !
    ! load indices for nearest neighbor interpolation from GOLD RedArk vec
    ! to tile space
    
    !fname = trim(RedArk_GOLD_coord_path) // '/' // trim(RedArk_GOLD_coord_file)
    fname = trim(met_path) // '/' // trim(RedArk_GOLD_coord_file)

    open (10, file=fname, form='formatted', action='read')
    
    do i=1,N_catd
       
       read (10,*) tmp_tile_id, ind_gold2tile(i)
       
       if (tmp_tile_id /= tile_coord(i)%tile_id) then
          err_msg = 'RedArk GOLD2tile error'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

    end do
    
    ! ------------------------------------------------------
    !
    ! get forcing data
    
    do this_var = 1,N_RedArk_GOLD_vars
       
       ! time dimension 
       !
       ! First entry in GOLD_RedArk file is at 0Z for states, 
       ! with fluxes for 0Z-6Z
       ! ------>  DIFFERENT FROM GSWP2 FORMAT!!!  <------
       !
       ! At 0Z for first day of month:
       !  - for fluxes read first entry of that month
       !  - for states read first entry of that month
       ! At 6Z for first day of month:
       !  - for fluxes read second entry of that month
       !  - for states read second entry of that month
       ! and so on...
       
       select case (this_var)
          
       case (1,2,3,4,5)   ! "fluxes"
          
          date_time_tmp = date_time
          
       case (6,7,8,9)     ! "states"
          
          date_time_tmp = date_time
          
          !!call augment_date_time(-int_dt_RedArk_GOLD_in_seconds,date_time_tmp)
       case default
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error')
          
       end select
       
       ! assemble year and month strings
       
       write (YYYY,'(i4.4)') date_time_tmp%year
       write (MM,  '(i2.2)') date_time_tmp%month
       write (DD,  '(i2.2)') date_time_tmp%day
       write (HHMM,'(i4.4)') date_time_tmp%hour*100 + date_time_tmp%min
       
       ! assemble file name, open file
       
       fname = trim(met_path) // '/' // trim(RedArk_GOLD_name(this_var)) &
            // '/Y' // YYYY // '/M' // MM // '/' &
            // trim(RedArk_GOLD_name(this_var)) // '_RedArk_' // &
            YYYY // MM // DD // '_' // HHMM
       
       if(root_logit) write (logunit,*) 'opening ' // trim(fname)
       
       open(10,file=fname,form='formatted',action='read')
       
       ! read compressed data, and put into tile space
       
       do i=1,RedArk_GOLD_N_cells
          
          read (10,*) tmp_vec(i)
          
       end do
       
       close(10,status='keep')
       
       do k=1,N_catd
          
          force_array(k,this_var) = tmp_vec( ind_gold2tile(k) )
          
       end do
       
    end do
    
    ! convert variables and units of force_array to match met_force_type, 
    ! put into structure
    
    ! from RedArk GOLD files:
    !
    !  force_array(:, 1) =  SWdown        W/m2 
    !  force_array(:, 2) =  LWdown        W/m2 
    !  force_array(:, 3) =  Rainf_C       kg/m2/s   
    !  force_array(:, 4) =  Rainf         kg/m2/s
    !  force_array(:, 5) =  Snowf         kg/m2/s    
    !  force_array(:, 6) =  PSurf         Pa
    !  force_array(:, 7) =  Qair          kg/kg  
    !  force_array(:, 8) =  Tair          K  
    !  force_array(:, 9) =  Wind          m/s   

    met_force_new%SWdown  = force_array(:,1)
    met_force_new%LWdown  = force_array(:,2)
    met_force_new%Rainf_C = force_array(:,3)
    met_force_new%Rainf   = force_array(:,4)
    met_force_new%Snowf   = force_array(:,5)
    met_force_new%Psurf   = force_array(:,6)
    met_force_new%Qair    = force_array(:,7)
    met_force_new%Tair    = force_array(:,8)
    met_force_new%Wind    = force_array(:,9)
    
  end subroutine get_RedArk_GOLD

  ! ***************************************************************************
  
  subroutine get_RedArk_Princeton(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )
    
    ! read RedArk Princeton files and map to forcings in tile space
    ! (uses pre-computed nearest neighbor interpolation)
    
    ! reichle, 6 Apr 2007
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    character(*),       intent(in) :: met_path
    
    integer,              intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new
    
    real, intent(out) :: nodata_forcing
    
    ! interpolation parameters
    
    integer, parameter :: RedArk_Princeton_N_cells = 70
    
    character(80),  parameter :: RedArk_Princeton_coord_file = &
         'Princeton_forcing_cells_2_RedArkOSSE_tiles.dat'    
    
    ! Princeton_REDARK forcing time step in hours
    
    integer, parameter :: dt_RedArk_Princeton_in_hours = 3
    
    integer, parameter :: N_RedArk_Princeton_vars = 7
    
    real,    parameter :: nodata_RedArk_Princeton = -9999.     ! ?????
    
    character(40), dimension(N_RedArk_Princeton_vars), parameter :: &
         RedArk_Princeton_name = (/ &
         'dswrf',   &   !  1 - flux
         'dlwrf',   &   !  2 - flux
         'prcp ',   &   !  3 - flux
         'pres ',   &   !  4 - state
         'shum ',   &   !  5 - state
         'tas  ',   &   !  6 - state
         'wind ' /)     !  7 - state

    ! local variables
    
    real :: tol 
    
    integer, dimension(N_catd)               :: ind_princeton2tile
    
    real,    dimension(RedArk_Princeton_N_cells)  :: tmp_vec

    real,    dimension(N_catd,N_RedArk_Princeton_vars) :: force_array
    
    type(date_time_type) :: date_time_tmp
    
    integer :: int_dt_RedArk_Princeton_in_secs, i, k, this_var, tmp_tile_id
    
    character(4) :: YYYY, HHMM
    character(2) :: MM, DD
    
    character(300) :: fname

    character(len=*), parameter :: Iam = 'get_RedArk_Princeton'
    character(len=400) :: err_msg
    
    ! --------------------------------------------------------------------
    
    int_dt_RedArk_Princeton_in_secs = 3600*dt_RedArk_Princeton_in_hours
    
    nodata_forcing = nodata_RedArk_Princeton
    
    tol = abs(nodata_forcing*nodata_tolfrac_generic)

    ! ----------------------------------------------    
    !
    ! load indices for nearest neighbor interpolation from Princeton RedArk vec
    ! to tile space
    
    fname = trim(met_path) // '/' // trim(RedArk_Princeton_coord_file)

    open (10, file=fname, form='formatted', action='read')
    
    do i=1,N_catd
       
       read (10,*) tmp_tile_id, ind_princeton2tile(i)
       
       if (tmp_tile_id /= tile_coord(i)%tile_id) then
          err_msg = 'RedArk Princeton2tile error'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

    end do
    
    ! ------------------------------------------------------
    !
    ! get forcing data
    
    do this_var = 1,N_RedArk_Princeton_vars
       
       ! time dimension 
       !
       ! *** STATEMENTS BELOW STILL NEED TO BE VERIFIED!!!! - reichle, 6 Apr 2007 ****
       !
       ! First entry in Princeton_RedArk file is at 0Z for states, 
       ! with fluxes for 0Z-3Z
       ! ------>  DIFFERENT FROM GSWP2 FORMAT!!!  <------
       !
       ! At 0Z for first day of month:
       !  - for fluxes read first entry of that month
       !  - for states read first entry of that month
       ! At 3Z for first day of month:
       !  - for fluxes read second entry of that month
       !  - for states read second entry of that month
       ! and so on...
       
       select case (this_var)
          
       case (1,2,3)   ! "fluxes"
          
          date_time_tmp = date_time
          
       case (4,5,6,7)     ! "states"
          
          date_time_tmp = date_time
          
          !!call augment_date_time(-int_dt_RedArk_Princeton_in_secs,date_time_tmp)
       case default
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error')
          
       end select
       
       ! assemble year and month strings
       
       write (YYYY,'(i4.4)') date_time_tmp%year
       write (MM,  '(i2.2)') date_time_tmp%month
       write (DD,  '(i2.2)') date_time_tmp%day
       write (HHMM,'(i4.4)') date_time_tmp%hour*100 + date_time_tmp%min
       
       ! assemble file name, open file
       
       fname = trim(met_path) // '/' // trim(RedArk_Princeton_name(this_var)) &
            // '/Y' // YYYY // '/M' // MM // '/' &
            // trim(RedArk_Princeton_name(this_var)) // '_RedArk_' // &
            YYYY // MM // DD // '_' // HHMM
       
       if(root_logit) write (logunit,*) 'opening ' // trim(fname)
       
       open(10,file=fname,form='formatted',action='read')
       
       ! read compressed data, and put into tile space
       
       do i=1,RedArk_Princeton_N_cells
          
          read (10,*) tmp_vec(i)
          
       end do
       
       close(10,status='keep')
       
       do k=1,N_catd
          
          force_array(k,this_var) = tmp_vec( ind_princeton2tile(k) )
          
       end do
       
    end do
    
    ! convert variables and units of force_array to match met_force_type, 
    ! put into structure

    ! from RedArk Princeton files:    
    !
    !  force_array(:, 1) =  SWdown        W/m2 
    !  force_array(:, 2) =  LWdown        W/m2 
    !  force_array(:, 3) =  RainfSnowf    kg/m2/s
    !  force_array(:, 4) =  PSurf         Pa
    !  force_array(:, 5) =  Qair          kg/kg  
    !  force_array(:, 6) =  Tair          K  
    !  force_array(:, 7) =  Wind          m/s   
    
    met_force_new%SWdown  = force_array(:,1)
    met_force_new%LWdown  = force_array(:,2)
    met_force_new%Psurf   = force_array(:,4)
    met_force_new%Qair    = force_array(:,5)
    met_force_new%Tair    = force_array(:,6)
    met_force_new%Wind    = force_array(:,7)
    
    do k=1,N_catd
       
       ! rain and snow:
       ! partition total precip into rain and snow according to Tair
       ! set convective precip to zero
       
       met_force_new(k)%Rainf_C = 0.
       
       if ( met_force_new(k)%Tair < Tzero ) then
          met_force_new(k)%Rainf = 0.
          met_force_new(k)%Snowf = force_array(k,3)
       else
          met_force_new(k)%Rainf = force_array(k,3)
          met_force_new(k)%Snowf = 0.
       end if
       
    end do
    
  end subroutine get_RedArk_Princeton

  ! ****************************************************************  

  subroutine get_Princeton_netcdf(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )
    
    ! read Princeton_NetCDF files and extract forcing in tile space
    ! (uses nearest neighbor interpolation)
    
    ! tyamada+reichle, 19 Jul 2007
   
    implicit none
 
    type(date_time_type), intent(in) :: date_time

    character(*),       intent(in) :: met_path

    integer,              intent(in) :: N_catd

    type(tile_coord_type), dimension(:), pointer :: tile_coord ! input

    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new

    real, intent(out) :: nodata_forcing
    
    ! Princeton_netcdf grid and netcdf parameters
    
    integer, parameter :: Princeton_grid_N_lon  = 360
    integer, parameter :: Princeton_grid_N_lat  = 180
    real,    parameter :: Princeton_grid_ll_lon =   0.
    real,    parameter :: Princeton_grid_ll_lat = -90.
    real,    parameter :: Princeton_grid_dlon   = 1.
    real,    parameter :: Princeton_grid_dlat   = 1.

    ! Princeton_netcdf forcing time step in hours
    
    integer, parameter :: dt_Princeton_in_hours = 3
    integer, parameter :: nciv_data             = 5 ! (1=lon, 2=lat, 3=z, 4=time, 5=data)
    integer, parameter :: N_Princeton_vars      = 7
    real,    parameter :: nodata_Princeton      = 2.e20 
    
    character(40), dimension(N_Princeton_vars), parameter :: Princeton_name = &
         (/ &
         'dswrf', &  ! 1 - flux
         'dlwrf', &  ! 2 - flux
         'prcp ', &  ! 3 - flux
         'pres ', &  ! 4 - state
         'shum ', &  ! 5 - state
         'tas  ', &  ! 6 - state
         'wind '  &  ! 7 - state  
         /)
    
    ! local variables
    
    integer, dimension(N_catd) :: i_ind, j_ind
    
    real, dimension(Princeton_grid_N_lon, Princeton_grid_N_lat) :: tmp_grid
    
    real, dimension(N_catd, N_Princeton_vars)                   :: force_array    
    
    integer, dimension(4)      :: start, icount
    
    integer                    :: k, hours_in_year, Princeton_var, ierr, ncid
    
    real                       :: tol, this_lon, this_lat, dt_Princeton_in_seconds
    
    character(  4)             :: YYYY, HHMM
    character(  2)             :: MM, DD
    character(300)             :: fname

    character(len=*), parameter :: Iam = 'get_Princeton_netcdf'
    character(len=400) :: err_msg
    
    ! ----------------------------------------------------------------
    
    dt_Princeton_in_seconds = real(3600*dt_Princeton_in_hours)
    
    nodata_forcing = nodata_Princeton

    tol = abs(nodata_forcing*nodata_tolfrac_generic)

    ! assemble year and month strings

    write (YYYY, '(i4.4)') date_time%year
    write (MM,   '(i2.2)') date_time%month
    write (DD,   '(i2.2)') date_time%day
    write (HHMM, '(i4.4)') date_time%hour*100 + date_time%min
    
    ! set lon index
    
    start(1)  = 1
    icount(1) = 360
    
    ! set lat index
    
    start(2)  = 1
    icount(2) = 180
    
    ! set z index
    
    start(3)  = 1
    icount(3) = 1
    
    ! get time index
    
    if ( (date_time%min/=0) .or. (date_time%sec/=0) .or.   &
         (mod(date_time%hour, dt_Princeton_in_hours)/=0) ) then
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'timing ERROR!!')
       
    endif
    
    hours_in_year = (date_time%dofyr-1)*24 + date_time%hour
    
    start(4)  = hours_in_year / dt_Princeton_in_hours + 1
    icount(4) = 1
    
    ! ---------------------------------------------------------------- 
    
    do k=1,N_catd
       
       this_lon = tile_coord(k)%com_lon
       this_lat = tile_coord(k)%com_lat

       ! change lon units for compatibility with Princeton netcdf
       ! grid which starts at the Greenwich Meridian and goes
       ! eastward (lon=0:360)
       
       if (this_lon<0.)  this_lon = this_lon + 360.
       
       i_ind(k) = ceiling((this_lon - Princeton_grid_ll_lon)/Princeton_grid_dlon)
       
       j_ind(k) = ceiling((this_lat - Princeton_grid_ll_lat)/Princeton_grid_dlat)
       
       ! not sure this is quite right -- reichle, 24 Feb 2009

       i_ind(k) = mod(i_ind(k), Princeton_grid_N_lon)

       ! fixed per suggestion from Greg Walker, - reichle, 26 Aug 2015

       if(i_ind(k) < 1) i_ind(k) = i_ind(k) + Princeton_grid_N_lon 

    enddo
    
    ! ----------------------------------------------------------------
    !
    ! open input file
    
    do Princeton_var = 1,N_Princeton_vars
       
       fname = trim(met_path) // '/' // trim(Princeton_name(Princeton_var))  &
            // '_3hourly_' // YYYY // '-' // YYYY // '.nc'
       
       if(root_logit) write(logunit,*) 'opening' // trim(fname)

       ierr = NF_OPEN(fname, NF_NOWRITE, ncid)
       
       if (ierr/=0) then
          err_msg = 'error opening netcdf file'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       ierr = NF_GET_VARA_REAL(ncid, nciv_data, start, icount, tmp_grid)
       
       ierr = NF_CLOSE(ncid)
       
       do k = 1, N_catd
          
          force_array(k, Princeton_var) = tmp_grid(i_ind(k), j_ind(k))      
          
       enddo
       
    enddo
    
    ! ---------------------------------------------------------------- 
    
    met_force_new%SWdown = force_array(:, 1)
    met_force_new%LWdown = force_array(:, 2)
    met_force_new%Psurf  = force_array(:, 4)
    met_force_new%Qair   = force_array(:, 5)
    met_force_new%Tair   = force_array(:, 6)
    met_force_new%Wind   = force_array(:, 7)
    
    do k=1, N_catd
       
       met_force_new(k)%Rainf_C = 0.
       
       if ( met_force_new(k)%Tair < Tzero ) then
          met_force_new(k)%Rainf = 0.
          met_force_new(k)%Snowf = force_array(k,3)
       else
          met_force_new(k)%Rainf = force_array(k,3)
          met_force_new(k)%Snowf = 0.
       endif
       
    enddo
    
  end subroutine get_Princeton_netcdf
  
  ! ****************************************************************  

  subroutine get_conus_netcdf(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )
    
    ! read conus_NetCDF files and extract forcing in tile space
    ! (uses nearest neighbor interpolation)
    
    ! sarith+reichle, 19 Jul 2007
   
    implicit none
 
    type(date_time_type), intent(in) :: date_time

    character(*),       intent(in) :: met_path

    integer,              intent(in) :: N_catd

    type(tile_coord_type), dimension(:), pointer :: tile_coord ! input

    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new

    real, intent(out) :: nodata_forcing
    
    ! conus_netcdf grid and netcdf parameters
    
    integer, parameter :: conus_grid_N_lon  = 115
    integer, parameter :: conus_grid_N_lat  = 48
    real,    parameter :: conus_grid_ll_lon = -125.
    real,    parameter :: conus_grid_ll_lat = 25.
    real,    parameter :: conus_grid_dlon   = 0.5
    real,    parameter :: conus_grid_dlat   = 0.5

    ! conus_netcdf forcing time step in hours
    
    integer, parameter :: dt_conus_in_hours = 1
 
    integer, parameter :: N_conus_vars      = 7
    real,    parameter :: nodata_conus      = 1.e20 
    
    character(40), dimension(N_conus_vars), parameter :: conus_name = &
         (/ &
         'fsds ',  &  ! 1 - flux
         'flds ',  &  ! 2 - flux
         'precs',  &  ! 3 - flux
         'tbot ',  &  ! 4 - state
         'wind ',  &  ! 5 - state
         'psrf ',  &  ! 6 - state
         'qbot '   &  ! 7 - state  
         /)
    
    integer :: nciv_data 
    
    ! local variables
    
    integer, dimension(N_catd) :: i_ind, j_ind
    
    real, dimension(conus_grid_N_lon, conus_grid_N_lat) :: tmp_grid
    
    real, dimension(N_catd, N_conus_vars)               :: force_array    
    
    integer, dimension(3)      :: start, icount
    
    integer                    :: k, hours_in_month, conus_var, ierr, ncid
    
    real                       :: tol, this_lon, this_lat, dt_conus_in_seconds
    
    character(  4)             :: YYYY, HHMM
    character(  2)             :: MM, DD
    character(300)             :: fname
    
    character(len=*), parameter :: Iam = 'get_conus_netcdf'
    character(len=400) :: err_msg

    ! ----------------------------------------------------------------
    
    dt_conus_in_seconds = real(3600*dt_conus_in_hours)
    
    nodata_forcing = nodata_conus

    tol = abs(nodata_forcing*nodata_tolfrac_generic)

    ! assemble year and month strings

    write (YYYY, '(i4.4)') date_time%year
    write (MM,   '(i2.2)') date_time%month
    write (DD,   '(i2.2)') date_time%day
    write (HHMM, '(i4.4)') date_time%hour*100 + date_time%min
    
    ! set lon index
    
    start(1)  = 1
    icount(1) = conus_grid_N_lon
    
    ! set lat index
    
    start(2)  = 1
    icount(2) = conus_grid_N_lat
    
    ! set z index
    
    !start(3)  = 1
    !icount(3) = 1
    
    ! get time index
    
    if ( (date_time%min/=0) .or. (date_time%sec/=0) .or.   &
         (mod(date_time%hour, dt_conus_in_hours)/=0) ) then
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'timing ERROR!!')
       
    endif

    hours_in_month = (date_time%day-1)*24 + date_time%hour    
    
    start(3)  = hours_in_month / dt_conus_in_hours + 1
    icount(3) = 1
    
    ! ---------------------------------------------------------------- 
    
    do k=1,N_catd
       
       this_lon = tile_coord(k)%com_lon
       this_lat = tile_coord(k)%com_lat
       
       i_ind(k) = ceiling((this_lon - conus_grid_ll_lon)/conus_grid_dlon)
       
       j_ind(k) = ceiling((this_lat - conus_grid_ll_lat)/conus_grid_dlat)
       
    enddo
    
    ! ----------------------------------------------------------------
    !
    ! open input file
    
    fname = trim(met_path) // '/' // YYYY//'-'//MM//'.nc'
    
    if(root_logit) write (logunit,*) 'opening' // trim(fname)
    
    ierr = NF_OPEN(fname, NF_NOWRITE, ncid)
    
    if (ierr/=0) then
       err_msg = 'error opening conus_netcdf file'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    do conus_var = 1,N_conus_vars		
       
       ! nciv_data = order of fields in netcdf file:
       !
       ! 5=lat, 6=lon, 7=fsds, 8=flds, 9=precs, 10=tbot, 11=wind, 12=psrf, 13=qbot
       
       nciv_data = 6 + conus_var
       
       ierr = NF_GET_VARA_REAL(ncid, nciv_data , start, icount, tmp_grid)
       
       do k = 1, N_catd
          
          force_array(k, conus_var) = tmp_grid(i_ind(k), j_ind(k))      
          
       end do
       
    enddo
    
    ierr = NF_CLOSE(ncid)
    
    ! ---------------------------------------------------------------- 
    
    met_force_new%SWdown = force_array(:, 1)
    met_force_new%LWdown = force_array(:, 2)
    met_force_new%Psurf  = force_array(:, 6)
    met_force_new%Qair   = force_array(:, 7)
    met_force_new%Tair   = force_array(:, 4)
    met_force_new%Wind   = force_array(:, 5)
    
    do k=1, N_catd
       
       met_force_new(k)%Rainf_C = 0.
       
       if ( met_force_new(k)%Tair < Tzero ) then
          met_force_new(k)%Rainf = 0.
          met_force_new(k)%Snowf = force_array(k,3)
       else
          met_force_new(k)%Rainf = force_array(k,3)
          met_force_new(k)%Snowf = 0.
       endif
       
    enddo
    
  end subroutine get_conus_netcdf


  ! ****************************************************************  
  
  subroutine get_ERA5_LIS(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )

    ! Read ERA5 NetCDF files maintained and shared by the NASA LIS group
    ! and based on forcing data put together originally for LDAS-Monde.
    !
    ! Forcing data are interpolated into tile space using nearest-neighbor;
    ! select MET_HINTERP accordingly during run configuration.
    !
    ! Note that the LDAS-Monde data used here *differ* from the published 
    ! ERA5 data that are available through the Copernicus Climate Change 
    ! Service (C3S) Climate Data Store.
    !
    ! Both the LDAS-Monde and the public data are on 1/4-deg lat/lon grids,
    ! but the grids are offset by 1/2 grid cell in both the lat and long directions.
    ! Also, the grid of the public data has 721 latitude points.
    !
    ! Moreover, the LDAS-Monde data provide temperature and humidity for the
    ! lowest atmospheric model level (~10 m), which are not available from the
    ! C3S Climate Data Store.
    !    
    ! borescan+reichle, 14 Apr 2021
    !
    ! ----------------------------------------------------------
    
    use netcdf
    implicit none
    include 'mpif.h'

    type(date_time_type), intent(in) :: date_time
    character(*),         intent(in) :: met_path
    integer,              intent(in) :: N_catd

    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input

    type(met_force_type) , dimension(N_catd), intent(inout) :: met_force_new
    real,                                     intent(out)   :: nodata_forcing

    ! ERA5 grid and netcdf parameters
    
    integer, parameter :: era5_grid_N_lon   = 1440
    integer, parameter :: era5_grid_N_lat   =  720
    real,    parameter :: era5_grid_ll_lon  = -180.0000
    real,    parameter :: era5_grid_ll_lat  =  -90.0000 
    real,    parameter :: era5_grid_dlon    =    0.25
    real,    parameter :: era5_grid_dlat    =    0.25
    integer, parameter :: N_era5_compressed = 340819
    integer, parameter :: N_era5_vars       = 9
    real,    parameter :: nodata_era5       = 1.e20
    
    character(40), dimension(N_era5_vars), parameter :: era5_name = (/ &
         'DIR_SWdown',   &   !  (1)
         'LWdown    ',   &   !  (2)
         'Snowf     ',   &   !  (3)
         'Rainf     ',   &   !  (4)
         'Tair      ',   &   !  (5)
         'Qair      ',   &   !  (6)
         'Wind      ',   &   !  (7)
         'PSurf     ',   &   !  (8) 
         'UREF      '/)      !  (9)  at 10m   

    ! local variables

    real,    dimension(era5_grid_N_lon,era5_grid_N_lat) :: tmp_grid
    integer, dimension(N_era5_compressed)               :: land_i_era5, land_j_era5, p2g
    integer, dimension(N_catd)                          :: i_ind, j_ind
    real,    dimension(N_era5_compressed)               :: tmp_vec
    real,    dimension(N_catd,N_era5_vars)              :: force_array
    integer, dimension(2)                               :: start, icount

    integer :: k, n, hours_in_month, era5_var, ierr, ncid, era5_varid, kk
    real    :: tol,  this_lon, this_lat  

    character(4)                :: YYYY
    character(2)                :: MM
    character(300)              :: fname
    character(len=*), parameter :: Iam = 'get_ERA5_LIS'
    character(len=400)          :: err_msg

    ! ----------------------------------------------------------------------------------------
    
    nodata_forcing = nodata_era5

    tol = abs(nodata_forcing*nodata_tolfrac_generic)

    ! assemble year and month strings
    write (YYYY,'(i4.4)') date_time%year
    write (MM,  '(i2.2)') date_time%month

    ! compressed space dimension (always read global vector)
    
    start(1)  = 1
    icount(1) = N_era5_compressed

    ! time dimension (first entry in ERA5_LIS file is at 00Z)

    if ( (date_time%min/=0) .or. (date_time%sec/=0) ) then

       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'timing ERROR!!')

    end if

    hours_in_month = (date_time%day-1)*24 + date_time%hour

    start(2)  = hours_in_month + 1
    icount(2) = 1

    ! compute indices for the nearest neighbor interpolation from ERA5 grid 
    ! to tile space

    do k=1,N_catd
       
       ! ll_lon and ll_lat refer to lower left corner of grid cell
       ! (as opposed to the grid point in the center of the grid cell)
       
       this_lon = tile_coord(k)%com_lon
       this_lat = tile_coord(k)%com_lat
       
       ! add small offsets to avoid unpredictable assignment of
       ! regularly spaced tiles (such as from the EASEv2 tile space)
       ! to ERA5 grid cells along certain lat/lon values
       ! (that is, make it possible for post-processing scripts in other
       !  languages to exactly reproduce the mapping that is done here)
       ! TO DO: add if statement so the offset is only applied when
       !        the model is run in the EASE grid tile space
       
       this_lon = this_lon + 0.0001
       this_lat = this_lat + 0.0001
       
       i_ind(k) = ceiling( (this_lon - era5_grid_ll_lon)/era5_grid_dlon )
       j_ind(k) = ceiling( (this_lat - era5_grid_ll_lat)/era5_grid_dlat )

    end do
    
    ! read parameters (same for all data variables and time steps)

    ! First read the point2grid (p2g) variable from the mapping.nc file
    ! 'point2grid' is used to calculate i and j indices for the tmp_grid

    fname = trim(met_path) // '/mapping.nc'

    if(root_logit) write (logunit,*) 'get compression index vector from ' // trim(fname)

    ierr = NF90_OPEN(fname,NF90_NOWRITE,ncid)
    if (ierr/=0) then
       err_msg = 'error opening netcdf file'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! read P2G from the nc file
    ierr = NF90_INQ_VARID(ncid,'P2G',era5_varid)
    ierr = NF90_GET_VAR(ncid, era5_varid, p2g)
    
    ! calculate i, j indices on era5_grid from 1-dim p2g index

    p2g = p2g-1
        
    land_i_era5 = mod( p2g,era5_grid_N_lon ) + 1
    land_j_era5 =      p2g/era5_grid_N_lon   + 1
    
    ! close NC file
    ierr = NF90_CLOSE(ncid)

    ! read Forcing data (for the date corresponding to LIS file) 

    ! open file
    fname = trim(met_path) // '/FORCING_' // YYYY // MM // '.nc'

    if(root_logit) write (logunit,*) 'opening ' // trim(fname)

    ierr = NF90_OPEN(fname,NF90_NOWRITE,ncid)

    if (ierr/=0) then
       err_msg = 'error opening netcdf file'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! loop through forcing variables

    do era5_var = 1,N_era5_vars
       
       ierr = NF90_INQ_VARID(ncid,trim(era5_name(era5_var)),era5_varid)   ! get var ID
       ierr = NF90_GET_VAR(ncid,era5_varid,tmp_vec,start,icount)          ! read variable
       
       ! put data on global grid
       
       tmp_grid  = nodata_forcing    ! initialize 
       
       do n=1,N_era5_compressed
          tmp_grid(land_i_era5(n),land_j_era5(n)) = tmp_vec(n)
       end do
       
       ! interpolate to tile space
       do k=1,N_catd
          force_array(k,era5_var) = tmp_grid(i_ind(k),j_ind(k))
       end do
       
    end do ! era5_var

    ! close NC file
    ierr = NF90_CLOSE(ncid)

    ! All variables in ERA5_LIS files have the units needed by met_force_type. 
    !
    !  force_array(:, 1) = DIR_SWdown W/m2     ; flux  ;(ssrd)
    !  force_array(:, 2) = LWdown     W/m2     ; flux  ;(strd)
    !  force_array(:, 3) = Snowf      kg/m2/s  ; flux  ;(sf)
    !  force_array(:, 4) = Rainf      kg/m2/s  ; flux  ;(cp+lsp-sf)
    !  force_array(:, 5) = Tair       K        ; state ;(corresponds to T at the lowest model layer (~10m))
    !  force_array(:, 6) = Qair       kg/kg    ; state ;(corresponds to Q at the lowest model layer (~10m))
    !  force_array(:, 7) = Wind       m/s      ; state ;(sqrt(u2+v2))
    !  force_array(:, 8) = PSurf      Pa       ; state ;(derived from lnsp)

    met_force_new%SWdown   = force_array(:,1)
    met_force_new%LWdown   = force_array(:,2)
    met_force_new%Snowf    = force_array(:,3)        
    met_force_new%Rainf    = force_array(:,4) 
    met_force_new%Rainf_C  = 0.                     ! always set convective precip to zero
    met_force_new%Tair     = force_array(:,5)
    met_force_new%Qair     = force_array(:,6)
    met_force_new%Wind     = force_array(:,7)
    met_force_new%Psurf    = force_array(:,8)
    met_force_new%RefH     = force_array(:,9)
    
    ! do not touch %PARdrct, and %PARdffs, which are not avaibable from ERA5_LIS
    
  end subroutine get_ERA5_LIS
  
  ! ****************************************************************  
  
  subroutine get_GLDAS_2x2_5_netcdf(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )
    
    ! read GLDAS_NetCDF files and extract forcings in tile space
    ! (uses nearest neighbor interpolation)
    
    ! reichle, 30 Jun 2005
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    character(*),       intent(in) :: met_path
    
    integer,              intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new
    
    real, intent(out) :: nodata_forcing
    
    ! GLDAS grid and netcdf parameters 
    
    integer, parameter :: gldas_grid_N_lon  = 144
    integer, parameter :: gldas_grid_N_lat  =  76
    real,    parameter :: gldas_grid_ll_lon = -181.25
    real,    parameter :: gldas_grid_ll_lat = -61.
    real,    parameter :: gldas_grid_dlon   = 2.5
    real,    parameter :: gldas_grid_dlat   = 2.
    
    integer, parameter :: N_gldas_compressed = 3023
    
    ! GLDAS forcing time step in hours

    integer, parameter :: dt_gldas_in_hours = 3
    
    integer, parameter :: nciv_land_i = 3
    integer, parameter :: nciv_land_j = 4
    integer, parameter :: nciv_data   = 7
    
    integer, parameter :: N_gldas_vars = 10
    
    real,    parameter :: nodata_gldas = 1.e20
    
    character(40), dimension(N_gldas_vars), parameter :: gldas_name = (/ &
         'SWdown_gldas      ',   &   !  1
         'LWdown_gldas      ',   &   !  2
         'Snowf_gldas       ',   &   !  3
         'Rainf_gldas       ',   &   !  4
         'Tair_gldas        ',   &   !  5
         'Qair_gldas        ',   &   !  6
         'Wind_E_gldas      ',   &   !  7
         'Wind_N_gldas      ',   &   !  8
         'PSurf_gldas       ',   &   !  9
         'RainfSnowf_C_gldas' /)     ! 10  ??? 
    
    ! local variables
    
    real :: tol 
    
    real, dimension(gldas_grid_N_lon,gldas_grid_N_lat) :: tmp_grid
    
    integer, dimension(N_gldas_compressed)   :: land_i_gldas, land_j_gldas
    integer, dimension(N_catd)               :: i_ind, j_ind
    
    real,    dimension(N_gldas_compressed)   :: tmp_vec
    
    real,    dimension(N_catd,N_gldas_vars) :: force_array
    
    integer, dimension(2) :: start, icount
    
    integer :: k, n, hours_in_month, gldas_var, ierr, ncid
    
    real    :: this_lon, this_lat, dt_gldas_in_seconds
    
    character(4) :: YYYY
    character(2) :: MM
    
    character(300) :: fname
    
    character(len=*), parameter :: Iam = 'get_GLDAS_netcdf'
    character(len=400) :: err_msg

    ! --------------------------------------------------------------------
    
    dt_gldas_in_seconds = real(3600*dt_gldas_in_hours)
    
    nodata_forcing = nodata_gldas
    
    tol = abs(nodata_forcing*nodata_tolfrac_generic)
    
    ! assemble year and month strings
    
    write (YYYY,'(i4.4)') date_time%year
    write (MM,  '(i2.2)') date_time%month
    
    ! find out which data are needed
    
    ! compressed space dimension (always read global vector)
    
    start(1)  = 1
    icount(1) = N_gldas_compressed
    
    ! time dimension (first entry in GLDAS_NetCDF file is at 0Z)
    
    if ( (date_time%min/=0) .or. (date_time%sec/=0) .or.           &
         (mod(date_time%hour,dt_gldas_in_hours)/=0)        ) then
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'timing ERROR!!')
              
    end if
    
    hours_in_month = (date_time%day-1)*24 + date_time%hour
    
    start(2)  = hours_in_month / dt_gldas_in_hours + 1 
    icount(2) = 1
    
    ! ----------------------------------------------    
    !
    ! compute indices for nearest neighbor interpolation from GLDAS grid 
    ! to tile space
    !
    ! (NOTE: this should at some point be replaced with a regridding
    !  subroutine that interpolates from the
    !  native forcing grid to the GCM atmospheric grid that is used
    !  to cut catchments into tiles - then "standard" grid2tile
    !  using tile_coord%atm_i and tile_coord%atm_j applies. 
    !  reichle, 26 May 2005)
    
    do k=1,N_catd
       
       ! ll_lon and ll_lat refer to lower left corner of grid cell
       ! (as opposed to the grid point in the center of the grid cell)
       
       this_lon = tile_coord(k)%com_lon
       this_lat = tile_coord(k)%com_lat
       
       i_ind(k) = ceiling( (this_lon - gldas_grid_ll_lon)/gldas_grid_dlon )
       j_ind(k) = ceiling( (this_lat - gldas_grid_ll_lat)/gldas_grid_dlat )
       
       ! NOTE: For a "date line on center" grid and (180-dlon/2) < lon < 180 
       ! we now have i_ind=(grid%N_lon+1) 
       ! This needs to be fixed.
       
       if (i_ind(k)>gldas_grid_N_lon)  i_ind(k)=1
       
    end do
    
    ! ------------------------------------------------------
    !
    ! read compression parameters (same for all data variables and time steps)
    
    gldas_var = 1
    
    fname = trim(met_path) // trim(gldas_name(gldas_var)) // '/' // YYYY     &
         // '/' // trim(gldas_name(gldas_var)) // '.' // YYYY // MM // '.nc'
    
    if(root_logit) write (logunit,*) 'get netcdf compression params from ' // trim(fname)
    
    ierr = NF_OPEN(fname,NF_NOWRITE,ncid)
    
    if (ierr/=0) then
       err_msg = 'error opening netcdf file'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    ierr = NF_GET_VARA_INT(ncid, nciv_land_i, start, icount, land_i_gldas)
    ierr = NF_GET_VARA_INT(ncid, nciv_land_j, start, icount, land_j_gldas)
    
    ierr = NF_CLOSE(ncid)

    ! ------------------------------------------------------
    !
    ! get forcing data
    
    do gldas_var = 1,N_gldas_vars
       
       ! open file, read compressed data, and put on global grid
       
       fname = trim(met_path) // trim(gldas_name(gldas_var)) // '/' // YYYY  &
            // '/' // trim(gldas_name(gldas_var)) // '.' // YYYY // MM //    &
            '.nc'
       
       if(root_logit) write (logunit,*) 'opening ' // trim(fname)
       
       ierr = NF_OPEN(fname,NF_NOWRITE,ncid)
       
       if (ierr/=0) then
          err_msg = 'error opening netcdf file'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       ierr = NF_GET_VARA_REAL(ncid, nciv_data, start, icount, tmp_vec )
       
       ierr = NF_CLOSE(ncid)

       tmp_grid = nodata_forcing

       do n=1,N_gldas_compressed
          
          tmp_grid(land_i_gldas(n), land_j_gldas(n) ) = tmp_vec(n)
          
       end do
       
       ! interpolate to tile space
       
       ! (NOTE: This should at some point be replaced with a regridding
       !  subroutine that interpolates from the
       !  native forcing grid to the GCM atmospheric grid that is used
       !  to cut catchments into tiles - then "standard" grid2tile
       !  using tile_coord%atm_i and tile_coord%atm_j applies. 
       !  reichle, 26 May 2005)
       
       do k=1,N_catd
          
          force_array(k,gldas_var) = tmp_grid(i_ind(k), j_ind(k))
          
       end do
       
    end do
    
    ! convert variables and units of force_array to match met_force_type, 
    ! put into structure
    
    ! from GLDAS files:
    !
    !  force_array(:, 1) = SWdown        W/m2 
    !  force_array(:, 2) = LWdown        W/m2 
    !  force_array(:, 3) = Snowf         kg/m2   (3h total)
    !  force_array(:, 4) = Rainf         kg/m2   (3h total)
    !  force_array(:, 5) = Tair          K    
    !  force_array(:, 6) = Qair          kg/kg
    !  force_array(:, 7) = Wind_E        m/s  
    !  force_array(:, 8) = Wind_N        m/s  
    !  force_array(:, 9) = PSurf         Pa   
    !  force_array(:,10) = RainfSnowf_C  kg/m2   (3h total)       ???
    
    met_force_new%SWdown  = force_array(:,1)
    met_force_new%LWdown  = force_array(:,2)
    met_force_new%Tair    = force_array(:,5)
    met_force_new%Qair    = force_array(:,6)
    met_force_new%Psurf   = force_array(:,9)
    
    ! get wind speed from wind vector
    
    do k=1,N_catd
       
       if ( abs(force_array(k,7)-nodata_gldas)<tol .or.           &
            abs(force_array(k,8)-nodata_gldas)<tol      ) then
          
          met_force_new(k)%Wind = nodata_forcing
          
       else
          
          met_force_new(k)%Wind = &
               sqrt( force_array(k,7)**2 + force_array(k,8)**2 )
          
       end if
       
    end do
    
    ! rain and snow:
    ! convert from GLDAS 3h totals [kg/m2] (or [mm]) into precipitation 
    ! rates [kg/m2/s] 
    
    do k=1,N_catd
       
       ! snowfall
       
       if (abs(force_array(k,3)-nodata_gldas)<tol) then
          
          met_force_new(k)%Snowf   = nodata_forcing
          
       else
          
          met_force_new(k)%Snowf   = force_array(k,3)/dt_gldas_in_seconds
          
       end if
       
       ! rainfall
       
       if (abs(force_array(k,4)-nodata_gldas)<tol) then
          
          met_force_new(k)%Rainf   = nodata_forcing
          
       else
          
          met_force_new(k)%Rainf   = force_array(k,4)/dt_gldas_in_seconds
          
       end if
       
       ! always set convective precip to zero (until it is clear
       ! whether GLDAS "RainfSnowf_C" a.k.a. "ACPCP" does indeed contain
       ! snow and was properly corrected with CMAP precip obs).
       
       met_force_new(k)%Rainf_C = 0.
       
    end do
    
  end subroutine get_GLDAS_2x2_5_netcdf

  ! ****************************************************************  
  
  subroutine get_Viviana_OK_precip(unitnumber, date_time, met_path, met_tag, &
       N_catd, tile_coord, met_force_new, ens_id )
    
    ! read Viviana's OK precip files and extract precip in tile space
    !
    ! NOTE: This subroutine only deals with precip.  Upon input, "met_force_new" 
    !       must contain "background" forcings for all other fields.
    !
    ! NOTE: If present, the optional argument "ens_id" indicates whether (and which)
    !       perturbed precip fields (after pre-processing with SREM2D) should be read
    !
    ! Assume that there are only "good" data in Viviana's OK precip files
    ! ("nodata" value is not used)
    !
    ! (uses nearest neighbor interpolation)
    
    ! vmaggion & reichle, 17 July 2008
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    character(*),       intent(in) :: met_path
    character(*),       intent(in) :: met_tag
    
    integer,              intent(in) :: unitnumber, N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input

    ! note "inout" for met_force_new

    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new
    
    ! NOTE: If present, the optional argument "ens_id" indicates whether and which 
    !       perturbed precip fields (after pre-processing with SREM2D) should be read
    
    integer, intent(in), optional :: ens_id
    
    ! Viviana precip OK grid parameters 
    !
    ! Viviana's ASCII forcing files are stored as follows:
    ! First line of forcing file is northernmost latitude band
    ! of the OK domain going from West to East.

    integer, parameter :: viviOK_grid_N_lon  =   22
    integer, parameter :: viviOK_grid_N_lat  =   10
    real,    parameter :: viviOK_grid_ll_lon = -100.0
    real,    parameter :: viviOK_grid_ll_lat =   34.5
    real,    parameter :: viviOK_grid_dlon   =    0.25
    real,    parameter :: viviOK_grid_dlat   =    0.25
    
    ! forcing time step in hours
    
    integer, parameter :: dt_viviOK_in_hours = 3
    
    character(40), parameter :: viviOK_name = 'precip_OK'

    ! local variables
    
    real, dimension(viviOK_grid_N_lon,viviOK_grid_N_lat) :: tmp_grid
    
    real, dimension(viviOK_grid_N_lon)                   :: tmp_vec
    
    integer, dimension(N_catd) :: i_ind, j_ind
    
    integer :: k, i, j
    
    real    :: this_lon, this_lat, dt_viviOK_in_seconds, tmpreal
    
    character(4) :: YYYY, ens_id_string
    character(2) :: MM, DD, HH
    
    character(300) :: fname

    character(len=*), parameter :: Iam = 'get_Viviana_OK_precip'
    character(len=400) :: err_msg
    
    ! --------------------------------------------------------------------
    
    dt_viviOK_in_seconds = real(3600*dt_viviOK_in_hours)
    
    ! assemble year and month strings
    
    write (YYYY,'(i4.4)') date_time%year
    write (MM,  '(i2.2)') date_time%month
    write (DD,  '(i2.2)') date_time%day
    write (HH,  '(i2.2)') date_time%hour
   
    ! ----------------------------------------------    
    !
    ! compute indices for nearest neighbor interpolation from viviOK grid 
    ! to tile space
    
    do k=1,N_catd
       
       ! ll_lon and ll_lat refer to lower left corner of grid cell
       ! (as opposed to the grid point in the center of the grid cell)
       
       this_lon = tile_coord(k)%com_lon
       this_lat = tile_coord(k)%com_lat
       
       i_ind(k) = ceiling( (this_lon - viviOK_grid_ll_lon)/viviOK_grid_dlon )
       j_ind(k) = ceiling( (this_lat - viviOK_grid_ll_lat)/viviOK_grid_dlat )
       
    end do
    
    ! ------------------------------------------------------
    !
    ! get precip data
    !
    ! open file, read data
    
    if     (index(met_tag, 'Viviana_OK_nopert')/=0) then
       
       fname = trim(met_path) // '/' // trim(viviOK_name) //                   &
            '/Y' // YYYY  // '/M' // MM //                                     &
            '/precip_' // YYYY // MM // DD // '_' // HH // 'z.txt'
       
    elseif (index(met_tag, 'Viviana_OK_pert')/=0) then
       
       if (present(ens_id)) then
          
          ! read perturbed forcing field (after processing through SREM2D)
          
          write (ens_id_string,'(i4.4)') ens_id
          
          fname = trim(met_path) // '/' // trim(viviOK_name) //                   &
               '/ens' // ens_id_string // '/Y' // YYYY  // '/M' // MM //             &
               '/precip_' // YYYY // MM // DD // '_' // HH // 'z.txt'
          
       else
          
          err_msg = 'need optional input argument "ens_id"'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end if
       
    else
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown met_tag')

    end if
       

    !unitnumber = 10
    !
    !if (present(ens_id)) unitnumber = unitnumber + ens_id
    
    if(root_logit)  write (logunit,*) 'opening ', trim(fname)
    
    open(unitnumber, file=fname, form='formatted', action='read', status='old')
    
    ! First line of ASCII forcing file is northernmost latitude band
    ! of the OK domain going from West to East.
    
    do j=viviOK_grid_N_lat,1,-1
       
       read (unitnumber,*) (tmp_vec(i), i=1,viviOK_grid_N_lon)
       
       tmp_grid(:,j) = tmp_vec
       
    end do

    close(unitnumber,status='keep')

    ! interpolate to tile space
    
    do k=1,N_catd
       
       ! also convert from viviOK 3h average rain rate [mm/h] into "met_force" 
       ! precipitation rates [kg/m2/s] 
       
       tmpreal = tmp_grid(i_ind(k),j_ind(k))/3600.
       
       ! set convective rainfall to zero, separate rain and snow
       
       met_force_new(k)%Rainf_C = 0.
       
       if ( met_force_new(k)%Tair < Tzero ) then
          met_force_new(k)%Rainf = 0.
          met_force_new(k)%Snowf = tmpreal
       else
          met_force_new(k)%Rainf = tmpreal
          met_force_new(k)%Snowf = 0.
       endif
       
    enddo
    
  end subroutine get_Viviana_OK_precip

  ! ************************************************************************
  
  subroutine get_GEOSs2s(date_time, met_path, met_tag, N_catd, tile_coord, &
       met_hinterp, met_force_new, nodata_forcing, PAR_available)
    
    ! read forcing derived from GEOS S2S output and map to tile space
    ! (using nearest neighbor or bilinear interpolation)
    !
    ! forcing derived through post-processing of daily average output from S2S
    ! hindcasts/forecasts ("FCST") or the "AODAS" used for initialization
    ! (see doc/README.MetForcing_and_BCS.md)
    !
    ! implementation follows LDASsa subroutines get_Princeton_netcdf() and get_GEOS(),
    ! fzeng, 24 Jun 2019
    !
    ! jkolassa,jmpark,reichle, 10 May - 14 June 2021:
    !   modified for GEOSldas; added AODAS; met_tag encodes ID of ensemble
    !   member and initial month/day
    !   modified to compute precipitation components from total precipitation
    !   and air temperature (pre-processing of precipitation components was faulty)
    
    use netcdf
    implicit none
    include 'mpif.h'
    
    type(date_time_type), intent(in) :: date_time
    
    character(*),         intent(in) :: met_path
    character(*),         intent(in) :: met_tag
    
    integer,              intent(in) :: N_catd, met_hinterp
    
    type(tile_coord_type), dimension(:), pointer       :: tile_coord     ! input
    
    type(met_force_type),  dimension(:), intent(inout) :: met_force_new
    
    real,                 intent(out) :: nodata_forcing
    
    logical,              intent(out) :: PAR_available

    ! Hindcast grid and netcdf parameters
    
    integer, parameter :: GEOSgcm_grid_N_lon  = 576
    integer, parameter :: GEOSgcm_grid_N_lat  = 361
    real,    parameter :: GEOSgcm_grid_dlon   = 0.625
    real,    parameter :: GEOSgcm_grid_dlat   = 0.5
    real,    parameter :: GEOSgcm_grid_ll_lon = -180. - GEOSgcm_grid_dlon/2.
    real,    parameter :: GEOSgcm_grid_ll_lat =  -90. - GEOSgcm_grid_dlat/2.
    
    real,    parameter :: nodata_GEOSgcm = -9999.
    
    integer, parameter :: dt_GEOSgcm_in_hours_FCST  =  6
    integer, parameter :: dt_GEOSgcm_in_hours_AODAS =  1
    
    integer, parameter :: N_GEOSgcm_vars_FCST       = 10
    integer, parameter :: N_GEOSgcm_vars_AODAS      =  8
    
    ! variable names in "FCST" forcing files match those in MERRA-2 file specs
    
    character(40), dimension(N_GEOSgcm_vars_FCST) :: GEOSgcm_name_FCST = &
         (/             &
         'PRECTOTCORR', &  !  1 - flux,  [gauge-corr]_total_precipitation,     kg m-2 s-1
         'LWGAB      ', &  !  2 - flux,  surface_absorbed_longwave_radiation,  W m-2
         'SWGDN      ', &  !  3 - flux,  surface_incoming_shortwave_flux,      W m-2
         'PARDR      ', &  !  4 - flux,  surface_downwelling_par_beam_flux,    W m-2
         'PARDF      ', &  !  5 - flux,  surface_downwelling_par_diffuse_flux, W m-2
         'PS         ', &  !  6 - state, surface_pressure,                     Pa
         'QLML       ', &  !  7 - state, surface_specific_humidity,            kg kg-1
         'TLML       ', &  !  8 - state, surface_air_temperature,              K
         'SPEED      ', &  !  9 - state, surface_wind_speed,                   m s-1
         'HLML       '  &  ! 10 - state, surface_layer_height,                 m
         /)

    ! variable names in "AODAS" forcing files match those in S2S file specs
    
    character(40), dimension(N_GEOSgcm_vars_AODAS) :: GEOSgcm_name_AODAS = &
         (/             &
         'TPREC      ', &  !  1 - flux,  total_precipitation,                  kg m-2 s-1
         'LWS        ', &  !  2 - flux,  surface_absorbed_longwave_radiation,  W m-2
         'SWGDWN     ', &  !  3 - flux,  surface_incoming_shortwave_flux,      W m-2
         'PS         ', &  !  4 - state, surface_pressure,                     Pa
         'QA         ', &  !  5 - state, surface_specific_humidity,            kg kg-1
         'TA         ', &  !  6 - state, surface_air_temperature,              K
         'SPEED      ', &  !  7 - state, surface_wind_speed,                   m s-1
         'HLML       '  &  !  8 - state, surface_layer_height,                 m
         /)
         
    ! local variables
    
    character(40), dimension(:), allocatable  :: GEOSgcm_name
    
    integer :: N_GEOSgcm_vars, dt_GEOSgcm_in_hours, N_lon_tmp, N_lat_tmp
    
    real    :: fnbr(2,2)
    
    integer, pointer :: i1(:), i2(:), j1(:), j2(:)
    real,    pointer :: x1(:), x2(:), y1(:), y2(:)
    
    real,    dimension(:,:), allocatable  :: force_array
    
    integer, dimension(3)      :: iicount, iistart
    integer                    :: k, hours_in_month, GEOSgcm_var, nciv_data
    integer                    :: fid, rc, nv_id, status
    
    real                       :: tol, this_lon, this_lat, Tair_tmp, Totprec_tmp
    
    character(  5)             :: init_MMMDD
    character(  4)             :: YYYY, ens_num
    character(  2)             :: MM, DD
    character(300)             :: fname
    
    character( 10)             :: lat_str = 'latitude'
    character( 10)             :: lon_str = 'longitude'
    
    logical                    :: FCST  = .false.
    logical                    :: AODAS = .false.
    
    character(len=*), parameter :: Iam = 'get_GEOSs2s'
    character(len=400)          :: err_msg

    ! --------------------------------------------------------------------
    !
    ! preparations

    ! parse met_tag

    ! 123456789012345678901234
    ! GEOSs2sFCST__ensX__MMMDD
    ! GEOSs2sAODAS

    if     (index(met_tag(8:11), 'FCST')/=0) then

       FCST          = .true.
       PAR_available = .true.

       dt_GEOSgcm_in_hours = dt_GEOSgcm_in_hours_FCST
       N_GEOSgcm_vars      = N_GEOSgcm_vars_FCST

       allocate(GEOSgcm_name(N_GEOSgcm_vars))

       GEOSgcm_name = GEOSgcm_name_FCST

       ! parse ensemble and initialization month from met_tag

       ens_num    = trim(met_tag(14:17))
       init_MMMDD = trim(met_tag(20:24))

    elseif (index(met_tag(8:12), 'AODAS')/=0) then

       AODAS         = .true.
       PAR_available = .false.

       dt_GEOSgcm_in_hours = dt_GEOSgcm_in_hours_AODAS
       N_GEOSgcm_vars      = N_GEOSgcm_vars_AODAS

       allocate(GEOSgcm_name(N_GEOSgcm_vars))

       GEOSgcm_name = GEOSgcm_name_AODAS

    else

       err_msg = "unknown met_tag"
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if

    nodata_forcing = nodata_GEOSgcm

    tol = abs(nodata_forcing*nodata_tolfrac_generic)

    ! assemble year and month strings

    write (YYYY, '(i4.4)') date_time%year
    write (MM,   '(i2.2)') date_time%month
    write (DD,   '(i2.2)') date_time%day

    ! set lon index

    iistart(1) = 1
    iicount(1) = GEOSgcm_grid_N_lon

    ! set lat index
    iistart(2) = 1
    iicount(2) = GEOSgcm_grid_N_lat

    ! get time index

    if ( (date_time%min/=0) .or. (date_time%sec/=0) .or.              &
         (mod(date_time%hour, dt_GEOSgcm_in_hours)/=0) ) then
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'timing ERROR!!')
       
    endif
    
    hours_in_month = (date_time%day-1)*24 + date_time%hour
    
    iistart(3) = hours_in_month / dt_GEOSgcm_in_hours + 1
    iicount(3) = 1
    
    ! ----------------------------------------------------------------
    !
    ! open input file

    if     (FCST)  then
       fname = trim(met_path) // '/' // YYYY // '/' // init_MMMDD // '/' // ens_num // '/GEOS5.' // YYYY // MM // '.nc4'
    elseif (AODAS) then
       fname = trim(met_path) // '/S2S_hourly_' // YYYY// MM// '.nc4'
    endif

    call GEOS_openfile(FileOpenedHash,fname,fid,tile_coord,met_hinterp,rc,lat_str,lon_str)
    if (rc<0) then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error opening file')
    endif

    ! assign indices from met forcing interpolation

    i1=>local_info%i1
    i2=>local_info%i2
    j1=>local_info%j1
    j2=>local_info%j2
    x1=>local_info%x1
    x2=>local_info%x2
    y1=>local_info%y1
    y2=>local_info%y2

    ! allocate force_array

    allocate(force_array(N_catd,N_GEOSgcm_vars))
    force_array = nodata_forcing

    ! loop over variables
    
    do GEOSgcm_var = 1,N_GEOSgcm_vars
       
       ! init shared memory
       N_lon_tmp = -1
       N_lat_tmp = -1
       if (associated(ptrShForce)) then
          N_lon_tmp = size(ptrShForce,1)
          N_lat_tmp = size(ptrShForce,2)
       endif
       if(  (N_lon_tmp /= GEOSgcm_grid_N_lon) .or.          &
            (N_lat_tmp /= GEOSgcm_grid_N_lat)       ) then
          call MAPL_SyncSharedMemory(rc=status)
          VERIFY_(status)
          if (associated(ptrShForce)) then
             call MAPL_DeallocNodeArray(ptrShForce,rc=status)
             VERIFY_(status)
          endif
          call MAPL_AllocateShared(ptrShForce,(/GEOSgcm_grid_N_lon,GEOSgcm_grid_N_lat/),TransRoot= .true.,rc=status)
          VERIFY_(status)
       end if

       call MAPL_SyncSharedMemory(rc=status)
       VERIFY_(status)
       
       ! read variable from netcdf file
       if (MAPL_AmNodeRoot .or. (.not. MAPL_ShmInitialized)) then
          rc= NF90_INQ_VARID( fid, trim(GEOSgcm_name(GEOSgcm_var)), nv_id)
          _ASSERT( rc == nf90_noerr, "nf90 error")
          rc= NF90_GET_VAR( fid, nv_id, ptrShForce, start=iistart,count=iicount)
       end if
       
       call MAPL_SyncSharedMemory(rc=status)
       
       ! map variable array to force array using chosen met interpolation method
       
       select case (MET_HINTERP)
          
       case(0)  ! nearest neighbor interpolation
          
          do k = 1, N_catd
             force_array(k, GEOSgcm_var) = ptrShForce(i1(k), j1(k))
          end do
          
       case (1)  ! bilinear interpolation
          
          do k=1,N_catd
             this_lon = tile_coord(k)%com_lon
             this_lat = tile_coord(k)%com_lat
             
             fnbr(1,1) = ptrShForce(i1(k),j1(k))
             fnbr(1,2) = ptrShForce(i1(k),j2(k))
             fnbr(2,1) = ptrShForce(i2(k),j1(k))
             fnbr(2,2) = ptrShForce(i2(k),j2(k))
             
             !DEC$ FORCEINLINE
             force_array(k,GEOSgcm_var) = BilinearInterpolation(this_lon, this_lat, &
                  x1(k), x2(k), y1(k), y2(k), fnbr, nodata_forcing, tol)
          end do
          
       case default
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unsupported MET_HINTERP option')
          
       end select
    end do ! GEOSgcm_var
    
    ! ----------------------------------------------------------------

    ! convert variables and units of force_array to match met_force_type,
    ! put into structure

    ! from GEOSgcm files:
    !
    !                     FCST/Hindcast
    !
    ! force_array(:, 1) = PRECTOTCORR    kg/m2/s ([gauge-corr]_total_precipitation, rainfall+snowfall)
    ! force_array(:, 2) = LWGAB          W/m2    (surface_absorbed_longwave_radiation)
    ! force_array(:, 3) = SWGDN          W/m2    (surface_incoming_shortwave_flux)
    ! force_array(:, 4) = PARDR          W/m2    (surface_downwelling_par_beam_flux)
    ! force_array(:, 5) = PARDF          W/m2    (surface_downwelling_par_diffuse_flux)
    ! force_array(:, 6) = PS             Pa      (surface_pressure)
    ! force_array(:, 7) = QLML           kg/kg   (surface_specific_humidity)
    ! force_array(:, 8) = TLML           K       (surface_air_temperature)
    ! force_array(:, 9) = SPEED          m/s     (surface_wind_speed)
    ! force_array(:,10) = HLML           m       (surface_layer_height)

    !                     AODAS
    !
    ! force_array(:, 1) = TPREC          kg/m2/s (total_precipitation, rainfall+snowfall)
    ! force_array(:, 2) = LWS            W/m2    (surface_absorbed_longwave_radiation)
    ! force_array(:, 3) = SWGDWN         W/m2    (surface_incoming_shortwave_flux)
    ! force_array(:, 4) = PS             Pa      (surface_pressure)
    ! force_array(:, 5) = QA             kg/kg   (surface_specific_humidity)
    ! force_array(:, 6) = TA             K       (surface_air_temperature)
    ! force_array(:, 7) = SPEED          m/s     (surface_wind_speed)
    ! force_array(:, 8) = HLML           m       (surface_layer_height)
    
    if     (FCST) then
       met_force_new%LWdown    = force_array(:, 2)
       met_force_new%SWdown    = force_array(:, 3)
       met_force_new%PARdrct   = force_array(:, 4)
       met_force_new%PARdffs   = force_array(:, 5)
       met_force_new%Psurf     = force_array(:, 6)
       met_force_new%Qair      = force_array(:, 7)
       met_force_new%Tair      = force_array(:, 8)
       met_force_new%Wind      = force_array(:, 9)
       met_force_new%RefH      = force_array(:,10)
    elseif (AODAS) then
       met_force_new%LWdown    = force_array(:, 2)
       met_force_new%SWdown    = force_array(:, 3)
       met_force_new%Psurf     = force_array(:, 4)
       met_force_new%Qair      = force_array(:, 5)
       met_force_new%Tair      = force_array(:, 6)
       met_force_new%Wind      = force_array(:, 7)
       met_force_new%RefH      = force_array(:, 8)
    end if
    
    ! rainfall
    ! Rainf = convective rainfall + large-scale rainfall (total liquid precip)
    ! Snowf = convective snowfall + large-scale snowfall (total solid precip)
    
    ! revised to compute precipitation components from total precipitation
    ! and air temperature because pre-processing of precipitation components
    ! was faulty - 14 Jun 2021
    
    do k=1,N_catd
       
       ! "where" statements would probably be better than a do loop, but left for later
       ! - reichle, 14 Jun 2021
       
       Tair_tmp    = met_force_new(k)%Tair  
       
       Totprec_tmp = force_array(k,1)       ! total precip is element 1 in FCST and AODAS
       
       if ( (abs(Tair_tmp   -nodata_GEOSgcm)<tol) .or.     &
            (abs(Totprec_tmp-nodata_GEOSgcm)<tol)          &
            ) then
          
          ! If Tair or total precip is undefined, no decision about precipitation 
          ! partitioning can be made and precip forcing is set to nodata_forcing.
          
          met_force_new(k)%Rainf_C = nodata_forcing
          met_force_new(k)%Snowf   = nodata_forcing
          met_force_new(k)%Rainf   = nodata_forcing
          
       else
          
          ! Use Tair to split total precipitation into rainfall and snowfall.
          
          if     (Tair_tmp <= Tzero) then ! Tair <= 0C
             
             met_force_new(k)%Snowf   = Totprec_tmp
             met_force_new(k)%Rainf   = 0.
             
          else                            ! Tair > 0C
             
             met_force_new(k)%Snowf   = 0.
             met_force_new(k)%Rainf   = Totprec_tmp
             
          end if

          ! Assign 0 to convective rainfall because only the total rainfall is used by the
          ! Catchment model (for now).  This may have to be revised for future model versions.
          
          met_force_new(k)%Rainf_C = 0.
          
       end if
       
    end do
    
    deallocate(force_array)
    deallocate(GEOSgcm_name)
    
  end subroutine get_GEOSs2s
  
  ! *************************************************************************
  
  subroutine get_GEOS( date_time, force_dtstep, met_path, met_tag,             &
       N_catd, tile_coord, MET_HINTERP, AEROSOL_DEPOSITION,                    &
       supported_option_MET_HINTERP,                                           &
       supported_option_AEROSOL_DEPOSITION,                                    &
       met_force_new, nodata_forcing, PAR_available, MERRA_file_specs,         &
       init )
    
    ! reichle,  5 March 2008 - adapted from get_GEOSgcm_gfio to work with DAS
    !                           and MERRA file specs
    ! reichle, 21 March 2008 - overhauled for time-interpolation of Tair etc
    !                          when read from MERRA tavg files
    ! qliu+reichle, 12 Aug 2008 - for MERRA, use TLML, QLML, ULML, VLML instead
    !                             of 2m variables
    !                           - different number of variables for DAS, MERRA
    !
    ! reichle, 23 Feb 2009 - read ParDrct, ParDffs from MERRA files
    !                      - new output variable "MERRA_file_specs"
    !
    ! reichle,  5 Mar 2009 - deleted ParDrct, ParDffs after testing found no impact
    ! 
    ! reichle,  1 Dec 2009 - optionally read netcdf files with corrected precip
    !                      - parse "met_tag" for seamless integration across 
    !                         MERRA streams; assemble MERRA "met_path"
    !                      - use only "lfo" files for MERRA
    !                       
    ! reichle, 20 Dec 2011 - reinstated "PARdrct" and "PARdffs" for MERRA-Land file specs
    !
    ! reichle, 27 Feb 2012 - renamed subroutine
    !                      - revised "DAS_defs", now called "G5DAS_defs"
    !                      - parse "met_tag" to decide whether to use "MERRA_defs" 
    !                        (rather than check for presence of "diag_sfc" file)
    ! pchakrab+reichle, 
    !          13 Jan 2014 - added bilinear interpolation option
    !
    ! reichle, 27 Jul 2015 - added MERRA-2 forcing
    !
    ! jkolassa+reichle, 17 Dec 2019 - fixed "met_path", "prec_path", and "met_tag" at stream boundaries
    !                                  (for MERRA-2, always point to publicly available files)
    !                               - updated comments
    !
    !
    ! -----------------------------------
    !
    ! Read surface met forcing from MERRA, MERRA-2 or GEOS-5 DAS (i.e., FP) based on parsing of "met_tag".
    !
    ! SEE "LDASsa_default_inputs_driver.nml" for more documentation of "met_tag" 
    ! and "met_path".
    !
    ! Try reading met files in directory "met_path/" first.  
    ! If this fails, try again in "met_path/met_tag/*/Yyyyy/Mmm/"
    !
    ! Read GEOSgcm (GEOS5) hdf or nc4 files (or nc files for corrected MERRA precip)
    ! and extract forcings in tile space (uses nearest neighbor interpolation).
    !
    ! Time convention for "met_force_new" as stated in get_forcing() does
    ! NOT apply to get_GEOS(), which reads forcing states (such as
    ! Tair) at date_time (=date_time_inst) and "previous" (or "backward-looking") forcing fluxes 
    ! (such as SWdown) at date_time_bkwd, rather than "subsequent" (or "forward-looking") 
    ! forcing fluxes (except for init=.true.).
    !
    ! Example: if date_time=1z, met_force_new will contain TLML at 1z
    !          and SWGDN for average from 0z to 1z (as stored in 0:30z file)
    !
    ! When LDASsa is integrated within the coupled GEOS5 DAS, initial (time-avg) 
    ! "tavg1_2d_*_Nx" files are not available.  Use optional "init" flag to 
    ! deal with this situation.  [NOT SURE THIS STILL MAKES SENSE, reichle, 17 Dec 2019]

    use netcdf
    implicit none

    ! intent in:
    
    type(date_time_type), intent(in) :: date_time           ! date/time of 'inst' forcing
    
    integer,              intent(in) :: force_dtstep

    ! e.g.: met_path = '/land/reichle/GEOS5_land_forcings/'
    !       met_tag  = 'js4rt_b7p1'                         (GEOSgcm exp label)
    
    character(*),         intent(in) :: met_path
    character(*),         intent(in) :: met_tag
    
    integer,              intent(in) :: N_catd, MET_HINTERP, AEROSOL_DEPOSITION

    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input

    ! intent inout:

    logical, intent(inout) :: supported_option_MET_HINTERP
    logical, intent(inout) :: supported_option_AEROSOL_DEPOSITION
    
    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new

    ! intent out:
    
    real,    intent(out) :: nodata_forcing
    
    logical, intent(out) :: PAR_available              
    logical, intent(out) :: MERRA_file_specs       ! original MERRA specs, not MERRA-2

    ! optional:
    
    logical, intent(in), optional :: init

    ! -----------------------------------------
    
    ! local variables
    
    integer, parameter :: N_G5DAS_vars   = 12   ! same as for MERRA-2 (excl Aerosol vars)
    integer, parameter :: N_MERRA_vars   = 13
    integer, parameter :: N_MERRA2_vars  = 12   ! same as for G5DAS (excl Aerosol vars)
    integer, parameter :: N_Aerosol_vars = 60   ! additional aerosol forcing vars for GOSWIM (w/ MERRA-2 only for now)

    integer, parameter :: N_MERRA2plusAerosol_vars = N_MERRA2_vars + N_Aerosol_vars
    
    integer, parameter :: N_defs_cols  =  5    

    real,    parameter :: nodata_GEOSgcm = 1.e15   
    
    character(40), dimension(N_G5DAS_vars,              N_defs_cols) :: G5DAS_defs
    character(40), dimension(N_MERRA_vars,              N_defs_cols) :: MERRA_defs
    character(40), dimension(N_MERRA2plusAerosol_vars,  N_defs_cols) :: M2INT_defs
    character(40), dimension(N_MERRA2plusAerosol_vars,  N_defs_cols) :: M2COR_defs

    character(40), dimension(:,:), allocatable :: GEOSgcm_defs

    ! NOTE: met_path, prec_path, and met_tag for current ('inst') time and fwd and bkwd 'tavg' 
    !       times differ at stream boundaries -- jkolassa+reichle, 17 Dec 2019

    type(date_time_type) :: date_time_inst, date_time_fwd, date_time_bkwd, date_time_tmp
    
    character(200)       :: met_path_inst,  met_path_fwd,  met_path_bkwd,  met_path_tmp
    character(200)       :: prec_path_inst, prec_path_fwd, prec_path_bkwd, prec_path_tmp
    character( 80)       :: met_tag_inst,   met_tag_fwd,   met_tag_bkwd,   met_tag_tmp  

    character(  3)       :: met_file_ext
    character(  3)       :: precip_corr_file_ext

    integer :: N_GEOSgcm_vars, N_lon_tmp, N_lat_tmp    

    real    :: this_lon, this_lat, tmp_lon, tmp_lat

    real    :: tol 
    
    real    :: fnbr(2,2)
    
    integer, pointer :: i1(:), i2(:), j1(:), j2(:)
    real,    pointer :: x1(:), x2(:), y1(:), y2(:)

    real,    dimension(:,:),      allocatable :: force_array
    
    integer :: j, k, GEOSgcm_var, fid, km, lm, nvars, ngatts, rc, YYYYMMDD, HHMMSS

    logical :: minimize_shift, use_prec_corr, use_Predictor, tmp_init

    logical :: daily_met_files
    
    integer :: nv_id, ierr, icount(3), istart(3),lonid,latid

    character(len=*), parameter :: Iam = 'get_GEOS'
    integer :: status
    character(len=400) :: err_msg
    character(len=300) :: fname_full
    logical :: file_exists, single_time_in_file

    ! -----------------------------------------------------------------------
    !
    ! define GEOS5 file specs 
    !
    ! columns of GEOSgcm_defs are as follows:
    !  1 - short variable name in gfio file
    !  2 - averaging mode in G5DAS/MERRA file ('tavg' or 'inst')
    !  3 - file tag (eg. 'bkg.sfc' or 'diag_sfc', or 'tavg1_2d_rad_Nx')
    !  4 - file dir ('ana' or 'diag')
    !  5 - treated as S="state" or F="flux" in subroutine interpolate_to_timestep()
    
    ! G5DAS (a.k.a. FP) file specs (default, unless otherwise specified via "met_tag")
    !
    ! lfo_inst/tavg data available from 11 Jun 2013 (start of GEOS-5 ADAS version 5.11)

    G5DAS_defs( 1,:)=[character(len=40):: 'SWGDN   ','tavg','tavg1_2d_lfo_Nx','diag','F'] 
    G5DAS_defs( 2,:)=[character(len=40):: 'LWGAB   ','tavg','tavg1_2d_lfo_Nx','diag','F']
    G5DAS_defs( 3,:)=[character(len=40):: 'PARDR   ','tavg','tavg1_2d_lfo_Nx','diag','F']
    G5DAS_defs( 4,:)=[character(len=40):: 'PARDF   ','tavg','tavg1_2d_lfo_Nx','diag','F']
    G5DAS_defs( 5,:)=[character(len=40):: 'PRECCU  ','tavg','tavg1_2d_lfo_Nx','diag','F']
    G5DAS_defs( 6,:)=[character(len=40):: 'PRECLS  ','tavg','tavg1_2d_lfo_Nx','diag','F']  
    G5DAS_defs( 7,:)=[character(len=40):: 'PRECSNO ','tavg','tavg1_2d_lfo_Nx','diag','F'] 
    G5DAS_defs( 8,:)=[character(len=40):: 'PS      ','inst','inst1_2d_lfo_Nx','diag','S']  
    G5DAS_defs( 9,:)=[character(len=40):: 'HLML    ','inst','inst1_2d_lfo_Nx','diag','S']
    G5DAS_defs(10,:)=[character(len=40):: 'TLML    ','inst','inst1_2d_lfo_Nx','diag','S']    
    G5DAS_defs(11,:)=[character(len=40):: 'QLML    ','inst','inst1_2d_lfo_Nx','diag','S']    
    G5DAS_defs(12,:)=[character(len=40):: 'SPEEDLML','inst','inst1_2d_lfo_Nx','diag','S']    


    ! MERRA-2 file specs with uncorrected (AGCM) precip from the "int" Collection
    !  (ie, the precip generated by the AGCM within the MERRA-2 system)
    !
    ! NOTE: This is *NOT* the precipitation seen by the land surface in the MERRA-2 system.
    !
    ! NOTE: Use SWGDN from the "rad" Collection because SWGDN in MERRA-2 "lfo"
    !       is averaged over land tiles only, unlike all other variables,
    !       which are global, as is SWGDN in the FP "lfo" files.
    !       - reichle, 7 Dec 2015

    M2INT_defs( 1,:)=[character(len=40):: 'SWGDN   ','tavg','tavg1_2d_rad_Nx','diag','F']  ! use "rad" Collection
    M2INT_defs( 2,:)=[character(len=40):: 'LWGAB   ','tavg','tavg1_2d_lfo_Nx','diag','F']
    M2INT_defs( 3,:)=[character(len=40):: 'PARDR   ','tavg','tavg1_2d_lfo_Nx','diag','F']
    M2INT_defs( 4,:)=[character(len=40):: 'PARDF   ','tavg','tavg1_2d_lfo_Nx','diag','F']
    M2INT_defs( 5,:)=[character(len=40):: 'PRECCU  ','tavg','tavg1_2d_int_Nx','diag','F']  ! uncorrected
    M2INT_defs( 6,:)=[character(len=40):: 'PRECLS  ','tavg','tavg1_2d_int_Nx','diag','F']  ! uncorrected
    M2INT_defs( 7,:)=[character(len=40):: 'PRECSN  ','tavg','tavg1_2d_int_Nx','diag','F']  ! uncorrected
    M2INT_defs( 8,:)=[character(len=40):: 'PS      ','inst','inst1_2d_lfo_Nx','diag','S']  
    M2INT_defs( 9,:)=[character(len=40):: 'HLML    ','inst','inst1_2d_lfo_Nx','diag','S']
    M2INT_defs(10,:)=[character(len=40):: 'TLML    ','inst','inst1_2d_lfo_Nx','diag','S']    
    M2INT_defs(11,:)=[character(len=40):: 'QLML    ','inst','inst1_2d_lfo_Nx','diag','S']    
    M2INT_defs(12,:)=[character(len=40):: 'SPEEDLML','inst','inst1_2d_lfo_Nx','diag','S']    

    M2INT_defs(13,:)=[character(len=40):: 'DUDP001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(14,:)=[character(len=40):: 'DUDP002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(15,:)=[character(len=40):: 'DUDP003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(16,:)=[character(len=40):: 'DUDP004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(17,:)=[character(len=40):: 'DUDP005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(18,:)=[character(len=40):: 'DUSV001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(19,:)=[character(len=40):: 'DUSV002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(20,:)=[character(len=40):: 'DUSV003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(21,:)=[character(len=40):: 'DUSV004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(22,:)=[character(len=40):: 'DUSV005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(23,:)=[character(len=40):: 'DUWT001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(24,:)=[character(len=40):: 'DUWT002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(25,:)=[character(len=40):: 'DUWT003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(26,:)=[character(len=40):: 'DUWT004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(27,:)=[character(len=40):: 'DUWT005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(28,:)=[character(len=40):: 'DUSD001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(29,:)=[character(len=40):: 'DUSD002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(30,:)=[character(len=40):: 'DUSD003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(31,:)=[character(len=40):: 'DUSD004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(32,:)=[character(len=40):: 'DUSD005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(33,:)=[character(len=40):: 'BCDP001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(34,:)=[character(len=40):: 'BCDP002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(35,:)=[character(len=40):: 'BCSV001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(36,:)=[character(len=40):: 'BCSV002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(37,:)=[character(len=40):: 'BCWT001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(38,:)=[character(len=40):: 'BCWT002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(39,:)=[character(len=40):: 'BCSD001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(40,:)=[character(len=40):: 'BCSD002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(41,:)=[character(len=40):: 'OCDP001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(42,:)=[character(len=40):: 'OCDP002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(43,:)=[character(len=40):: 'OCSV001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(44,:)=[character(len=40):: 'OCSV002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(45,:)=[character(len=40):: 'OCWT001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(46,:)=[character(len=40):: 'OCWT002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(47,:)=[character(len=40):: 'OCSD001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(48,:)=[character(len=40):: 'OCSD002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(49,:)=[character(len=40):: 'SUDP003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(50,:)=[character(len=40):: 'SUSV003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(51,:)=[character(len=40):: 'SUWT003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(52,:)=[character(len=40):: 'SUSD003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(53,:)=[character(len=40):: 'SSDP001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(54,:)=[character(len=40):: 'SSDP002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(55,:)=[character(len=40):: 'SSDP003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(56,:)=[character(len=40):: 'SSDP004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(57,:)=[character(len=40):: 'SSDP005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(58,:)=[character(len=40):: 'SSSV001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(59,:)=[character(len=40):: 'SSSV002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(60,:)=[character(len=40):: 'SSSV003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(61,:)=[character(len=40):: 'SSSV004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(62,:)=[character(len=40):: 'SSSV005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(63,:)=[character(len=40):: 'SSWT001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(64,:)=[character(len=40):: 'SSWT002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(65,:)=[character(len=40):: 'SSWT003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(66,:)=[character(len=40):: 'SSWT004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(67,:)=[character(len=40):: 'SSWT005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(68,:)=[character(len=40):: 'SSSD001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(69,:)=[character(len=40):: 'SSSD002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(70,:)=[character(len=40):: 'SSSD003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(71,:)=[character(len=40):: 'SSSD004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2INT_defs(72,:)=[character(len=40):: 'SSSD005    ','tavg','tavg1_2d_adg_Nx','diag','F']


    ! MERRA-2 file specs with corrected precip, which could be either
    !  - native (ie, the precip seen by the land surface in the MERRA-2 system), or
    !  - corrected in post-processing using MERRA-2 (uncorrected) precip as the background
    ! The default is to use MERRA-2 native precip corrections.  If the "met_tag" includes
    ! an optional "__prec[xyz]" string, the precip corrections specified by [xyz] are used.
    !
    ! NOTE: This is  *NOT* the same as the corrected precipitation of the off-line 
    !       spin-up run used to generate the MERRA-2 land surface initial conditions 
    !       for each stream.  These precip files used for that have a MERRA background.
    !
    ! NOTE: Use SWGDN from the "rad" Collection (see comment above).

    M2COR_defs( 1,:)=[character(len=40):: 'SWGDN      ','tavg','tavg1_2d_rad_Nx','diag','F']  ! use "rad" Collection
    M2COR_defs( 2,:)=[character(len=40):: 'LWGAB      ','tavg','tavg1_2d_lfo_Nx','diag','F']
    M2COR_defs( 3,:)=[character(len=40):: 'PARDR      ','tavg','tavg1_2d_lfo_Nx','diag','F']
    M2COR_defs( 4,:)=[character(len=40):: 'PARDF      ','tavg','tavg1_2d_lfo_Nx','diag','F']
    M2COR_defs( 5,:)=[character(len=40):: 'PRECCUCORR ','tavg','tavg1_2d_lfo_Nx','diag','F']  ! MERRA-2 built-in corrections
    M2COR_defs( 6,:)=[character(len=40):: 'PRECLSCORR ','tavg','tavg1_2d_lfo_Nx','diag','F']  ! MERRA-2 built-in corrections  
    M2COR_defs( 7,:)=[character(len=40):: 'PRECSNOCORR','tavg','tavg1_2d_lfo_Nx','diag','F']  ! MERRA-2 built-in corrections 
    M2COR_defs( 8,:)=[character(len=40):: 'PS         ','inst','inst1_2d_lfo_Nx','diag','S']  
    M2COR_defs( 9,:)=[character(len=40):: 'HLML       ','inst','inst1_2d_lfo_Nx','diag','S']
    M2COR_defs(10,:)=[character(len=40):: 'TLML       ','inst','inst1_2d_lfo_Nx','diag','S']    
    M2COR_defs(11,:)=[character(len=40):: 'QLML       ','inst','inst1_2d_lfo_Nx','diag','S']    
    M2COR_defs(12,:)=[character(len=40):: 'SPEEDLML   ','inst','inst1_2d_lfo_Nx','diag','S']    

    M2COR_defs(13,:)=[character(len=40):: 'DUDP001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(14,:)=[character(len=40):: 'DUDP002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(15,:)=[character(len=40):: 'DUDP003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(16,:)=[character(len=40):: 'DUDP004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(17,:)=[character(len=40):: 'DUDP005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(18,:)=[character(len=40):: 'DUSV001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(19,:)=[character(len=40):: 'DUSV002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(20,:)=[character(len=40):: 'DUSV003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(21,:)=[character(len=40):: 'DUSV004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(22,:)=[character(len=40):: 'DUSV005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(23,:)=[character(len=40):: 'DUWT001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(24,:)=[character(len=40):: 'DUWT002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(25,:)=[character(len=40):: 'DUWT003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(26,:)=[character(len=40):: 'DUWT004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(27,:)=[character(len=40):: 'DUWT005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(28,:)=[character(len=40):: 'DUSD001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(29,:)=[character(len=40):: 'DUSD002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(30,:)=[character(len=40):: 'DUSD003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(31,:)=[character(len=40):: 'DUSD004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(32,:)=[character(len=40):: 'DUSD005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(33,:)=[character(len=40):: 'BCDP001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(34,:)=[character(len=40):: 'BCDP002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(35,:)=[character(len=40):: 'BCSV001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(36,:)=[character(len=40):: 'BCSV002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(37,:)=[character(len=40):: 'BCWT001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(38,:)=[character(len=40):: 'BCWT002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(39,:)=[character(len=40):: 'BCSD001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(40,:)=[character(len=40):: 'BCSD002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(41,:)=[character(len=40):: 'OCDP001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(42,:)=[character(len=40):: 'OCDP002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(43,:)=[character(len=40):: 'OCSV001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(44,:)=[character(len=40):: 'OCSV002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(45,:)=[character(len=40):: 'OCWT001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(46,:)=[character(len=40):: 'OCWT002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(47,:)=[character(len=40):: 'OCSD001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(48,:)=[character(len=40):: 'OCSD002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(49,:)=[character(len=40):: 'SUDP003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(50,:)=[character(len=40):: 'SUSV003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(51,:)=[character(len=40):: 'SUWT003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(52,:)=[character(len=40):: 'SUSD003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(53,:)=[character(len=40):: 'SSDP001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(54,:)=[character(len=40):: 'SSDP002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(55,:)=[character(len=40):: 'SSDP003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(56,:)=[character(len=40):: 'SSDP004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(57,:)=[character(len=40):: 'SSDP005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(58,:)=[character(len=40):: 'SSSV001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(59,:)=[character(len=40):: 'SSSV002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(60,:)=[character(len=40):: 'SSSV003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(61,:)=[character(len=40):: 'SSSV004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(62,:)=[character(len=40):: 'SSSV005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(63,:)=[character(len=40):: 'SSWT001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(64,:)=[character(len=40):: 'SSWT002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(65,:)=[character(len=40):: 'SSWT003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(66,:)=[character(len=40):: 'SSWT004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(67,:)=[character(len=40):: 'SSWT005    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(68,:)=[character(len=40):: 'SSSD001    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(69,:)=[character(len=40):: 'SSSD002    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(70,:)=[character(len=40):: 'SSSD003    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(71,:)=[character(len=40):: 'SSSD004    ','tavg','tavg1_2d_adg_Nx','diag','F']
    M2COR_defs(72,:)=[character(len=40):: 'SSSD005    ','tavg','tavg1_2d_adg_Nx','diag','F']

    
    ! MERRA file specs
    !
    ! use *only* "tavg" files b/c "bkg.sfc" files are available only every 6h
    !
    ! - replaced 'SWGNT' from MERRA 'rad' file with 'SWLAND' from 'lnd' file
    ! - changed 'ULML', 'VLML' from 'S' to 'F' 
    ! - added 'PARDR', 'PARDF'
    ! reichle, 24 Feb 2009
    ! - deleted 'PARDR', 'PARDF' after testing showed no impact (REINSTATED Dec 2011)
    ! reichle,  5 Mar 2009
    ! - use "lfo" files
    ! reichle,  1 Dec 2009
    
    !                                                                       MERRA
    !                                                                     collection
    
    MERRA_defs( 1,:)=[character(len=40):: 'SWGDN  ','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "rad"
    MERRA_defs( 2,:)=[character(len=40):: 'LWGAB  ','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "rad"
    MERRA_defs( 3,:)=[character(len=40):: 'PARDR  ','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "lnd"
    MERRA_defs( 4,:)=[character(len=40):: 'PARDF  ','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "lnd"
    MERRA_defs( 5,:)=[character(len=40):: 'PRECTOT','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "lnd"
    MERRA_defs( 6,:)=[character(len=40):: 'PRECCON','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "flx"
    MERRA_defs( 7,:)=[character(len=40):: 'PRECSNO','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "lnd"
    MERRA_defs( 8,:)=[character(len=40):: 'PS     ','tavg','tavg1_2d_lfo_Nx','diag','S']    ! "slv"
    MERRA_defs( 9,:)=[character(len=40):: 'HLML   ','tavg','tavg1_2d_lfo_Nx','diag','S']    ! "flx"
    MERRA_defs(10,:)=[character(len=40):: 'TLML   ','tavg','tavg1_2d_lfo_Nx','diag','S']    ! "flx"
    MERRA_defs(11,:)=[character(len=40):: 'QLML   ','tavg','tavg1_2d_lfo_Nx','diag','S']    ! "flx"
    MERRA_defs(12,:)=[character(len=40):: 'ULML   ','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "flx"
    MERRA_defs(13,:)=[character(len=40):: 'VLML   ','tavg','tavg1_2d_lfo_Nx','diag','F']    ! "flx"

    ! --------------------------------------------------------------------
    !
    ! preparations
    
    tmp_init = .false.
    
    if (present(init))  tmp_init = init
    
    use_prec_corr = .false.  ! use corrected precip dataset (other than native "M2COR_defs")
    
    ! use same nodata-value on input and output so that "nodata-check" can be
    ! omitted when no arithmetic is needed
    
    nodata_forcing = nodata_GEOSgcm 
    
    tol = abs(nodata_forcing*nodata_tolfrac_generic)    

    ! all GEOS forcing datasets provide PAR (so far)

    PAR_available = .true.
    
    ! input variable "date_time" is for reading instantaneous ('inst') forcing variables

    date_time_inst = date_time
    
    ! assemble date/time structures for tavg files 
    ! (e.g. 0z-1z time average is in file with 0:30z timestamp)
    
    date_time_bkwd = date_time_inst
    
    call augment_date_time( -force_dtstep/2, date_time_bkwd )

    date_time_fwd  = date_time_inst
    
    call augment_date_time( +force_dtstep/2, date_time_fwd )

    ! ---------------------------------------------------------------------------
    
    ! determine which file specs should be used (MERRA, MERRA-2, or G5DAS)
    
    ! initialize to most likely values, overwrite below as needed
    
    MERRA_file_specs     = .false.
    
    met_file_ext         = 'nc4'       
    
    daily_met_files      = .false.
    
    precip_corr_file_ext = 'nc4'
    
    if     (met_tag(4:8)=='merra') then   ! MERRA
       
       ! AEROSOL_DEPOSITION /= 0 is NOT supported

       N_GEOSgcm_vars = N_MERRA_vars
       
       allocate(GEOSgcm_defs(N_GEOSgcm_vars,N_defs_cols))

       GEOSgcm_defs(1:N_GEOSgcm_vars,:) = MERRA_defs
       
       MERRA_file_specs = .true.

       met_file_ext = 'hdf'

       precip_corr_file_ext = 'nc '
              
       call parse_MERRA_met_tag( met_path, met_tag, date_time_inst,          &
            met_path_inst, prec_path_inst, met_tag_inst, use_prec_corr )

       call parse_MERRA_met_tag( met_path, met_tag, date_time_fwd,           &
            met_path_fwd,  prec_path_fwd,  met_tag_fwd,  use_prec_corr )

       call parse_MERRA_met_tag( met_path, met_tag, date_time_bkwd,          &
            met_path_bkwd, prec_path_bkwd, met_tag_bkwd, use_prec_corr )

    elseif (met_tag(1:2)=='M2') then      ! MERRA-2
       
       select case (AEROSOL_DEPOSITION)

       case (0)  ! no aerosols (default)
          
          N_GEOSgcm_vars = N_MERRA2_vars
          
       case (1)  ! with aerosols
          
          supported_option_AEROSOL_DEPOSITION = .true.  
          
          N_GEOSgcm_vars = N_MERRA2plusAerosol_vars 

       case default
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unsupported AEROSOL_DEPOSITION option')
          
       end select
       
       allocate(GEOSgcm_defs(N_GEOSgcm_vars,N_defs_cols))
       
       if     (met_tag(1:5)=='M2INT') then
          
          GEOSgcm_defs(1:N_GEOSgcm_vars,:) = M2INT_defs(1:N_GEOSgcm_vars,:)

       elseif (met_tag(1:5)=='M2COR') then
          
          GEOSgcm_defs(1:N_GEOSgcm_vars,:) = M2COR_defs(1:N_GEOSgcm_vars,:)

       else
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown "M2[xxx]" met_tag')
 
       end if
       
       daily_met_files = .true.
       
       call parse_MERRA2_met_tag( met_path, met_tag, date_time_inst,         &
            met_path_inst, prec_path_inst, met_tag_inst, use_prec_corr )
  
       call parse_MERRA2_met_tag( met_path, met_tag, date_time_fwd,          &
            met_path_fwd,  prec_path_fwd,  met_tag_fwd,  use_prec_corr )

       call parse_MERRA2_met_tag( met_path, met_tag, date_time_bkwd,         &
            met_path_bkwd, prec_path_bkwd, met_tag_bkwd, use_prec_corr )
       
    else     ! GEOS ADAS (FP)
       
       ! AEROSOL_DEPOSITION /= 0 is NOT supported
       
       N_GEOSgcm_vars = N_G5DAS_vars 
       
       allocate(GEOSgcm_defs(N_GEOSgcm_vars,N_defs_cols))

       GEOSgcm_defs(1:N_G5DAS_vars,  :) = G5DAS_defs
       
       call parse_G5DAS_met_tag( met_path, met_tag, date_time_inst,          &
            met_path_inst, prec_path_inst, met_tag_inst, use_prec_corr,      &
            use_Predictor )

       call parse_G5DAS_met_tag( met_path, met_tag, date_time_fwd,           &
            met_path_fwd,  prec_path_fwd,  met_tag_fwd,  use_prec_corr,      &
            use_Predictor )

       call parse_G5DAS_met_tag( met_path, met_tag, date_time_bkwd,          &
            met_path_bkwd, prec_path_bkwd, met_tag_bkwd, use_prec_corr,      &
            use_Predictor )
              
       if (use_Predictor) then
          
          ! append "+-" to GCM file tag (ie, replace "Nx" with "Nx+-")
          
          do j=1,N_GEOSgcm_vars
             
             GEOSgcm_defs(j,3) = trim(GEOSgcm_defs(j,3)) // '+-'
             
          end do
          
       end if
       
    end if
    
    allocate(force_array(N_catd,N_GEOSgcm_vars))
    
    ! ---------------------------------------------------------------------------
    !
    ! get forcing data

    do GEOSgcm_var = 1,N_GEOSgcm_vars

       ! open GEOS file (G5DAS or MERRA or MERRA-2)
       ! 
       ! Initial "tavg1_2d_*_Nx" files may not be available.  In this case,
       ! use first available file.  For G5DAS file specs, only "PS" is affected 
       ! (because the fluxes from the missing initial "tavg1_2d_*_Nx" file are not needed).
       ! For MERRA file specs, air temp, humidity, wind, and pressure are affected.
       !
       ! if (init==.false.) the j loop ends for j=1 (either successfully open file
       ! or stop).
       !
       ! if (init==.true.) make second attempt (j=2) to allow for possibly 
       ! missing "diag_sfc" or "tavg" file at date_time_bkwd (and try reading 
       ! the file at date_time_fwd).
       
       do j=1,2
          
          ! determine time stamp on file and corresponding met_path, prec_path, & met_tag
          
          if      (trim(GEOSgcm_defs(GEOSgcm_var,2))=='tavg') then

             if (j==1) then
                date_time_tmp  = date_time_bkwd
                met_path_tmp   = met_path_bkwd
                prec_path_tmp  = prec_path_bkwd
                met_tag_tmp    = met_tag_bkwd
             else
                date_time_tmp  = date_time_fwd
                met_path_tmp   = met_path_fwd
                prec_path_tmp  = prec_path_fwd
                met_tag_tmp    = met_tag_fwd
             end if
             
          else if (trim(GEOSgcm_defs(GEOSgcm_var,2))=='inst' ) then
             
             date_time_tmp     = date_time_inst
             met_path_tmp      = met_path_inst
             prec_path_tmp     = prec_path_inst
             met_tag_tmp       = met_tag_inst
          else
             
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown GEOSgcm_defs(2)')
             
          end if
          
          YYYYMMDD = date_time_tmp%year*10000+date_time_tmp%month*100+date_time_tmp%day
          HHMMSS   = date_time_tmp%hour*10000+date_time_tmp%min*100  +date_time_tmp%sec
          
          ! determine forcing file name (with path) 

          if ( (use_prec_corr) .and. (GEOSgcm_defs(GEOSgcm_var,1)(1:4)=='PREC') ) then
             
             call get_GEOS_corr_prec_filename(fname_full,file_exists,date_time_tmp,        &
                  prec_path_tmp, met_tag_tmp, GEOSgcm_defs(GEOSgcm_var,:), precip_corr_file_ext )

             single_time_in_file = .true.  ! corr precip files are always hourly (incl. MERRA-2)

          else
             
             call get_GEOS_forcing_filename(fname_full,file_exists,date_time_tmp, &
                  daily_met_files, met_path_tmp, met_tag_tmp,                     &
                  GEOSgcm_defs(GEOSgcm_var,:), met_file_ext)

             single_time_in_file = .not. daily_met_files  ! MERRA-2 files are daily files

          end if
          
          if ( file_exists) then
             
             exit  ! exit j loop after successfully finding file
             
          elseif (                                                               &
               (j==1)                                       .and.                &
               (tmp_init)                                   .and.                &
               (trim(GEOSgcm_defs(GEOSgcm_var,2))=='tavg')  .and.                &
               (root_logit)                                       ) then
             
             if (.not. MERRA_file_specs)  write (logunit,'(400A)')               &
                  'NOTE: Initialization. Data from tavg file are not used '  //  &
                  'with lfo inst/tavg forcing, but dummy values must be '    //  &
                  'read from some file for backward compatibility with '     //  &
                  'MERRA forcing.'
             
             write (logunit,*) 'try again with different file...'
             
          else
             
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error finding met forcing file')
             
          end if
          
       end do  ! j=1,2

       ! open file, extract coord info, prep horizontal interpolation info (if not done already)

       call GEOS_openfile(FileOpenedHash,fname_full,fid,tile_coord,MET_HINTERP)

       !fid = ptrNode%fid 

       i1=>local_info%i1
       i2=>local_info%i2
       j1=>local_info%j1
       j2=>local_info%j2
       x1=>local_info%x1
       x2=>local_info%x2
       y1=>local_info%y1
       y2=>local_info%y2

       ! ----------------------------------------------    
       !
       ! process grid dimensions
       ! NOTE: corrected precipitation forcing from separate netcdf can be on different grid
       
       ! init shared memory
       N_lon_tmp = -1
       N_lat_tmp = -1
       if (associated(ptrShForce)) then
          N_lon_tmp = size(ptrShForce,1)
          N_lat_tmp = size(ptrShForce,2)
       endif
       if ( N_lon_tmp /= local_info%N_lon .or.            &
            N_lat_tmp /= local_info%N_lat      ) then
          call MAPL_SyncSharedMemory(rc=status)
          VERIFY_(status)
          if (associated(ptrShForce)) then
             call MAPL_DeallocNodeArray(ptrShForce,rc=status)
             VERIFY_(status)
          endif
          call MAPL_AllocateShared(ptrShForce,(/local_info%N_lon,local_info%N_lat/),TransRoot= .true.,rc=status)
          VERIFY_(status)
          call MAPL_SyncSharedMemory(rc=status)
          VERIFY_(status)
       end if
       
       ! ----------------------------------------------    
       !
       ! read global gridded field of given variable
       
       call LDAS_GetVar( fid, trim(GEOSgcm_defs(GEOSgcm_var,1)),                  &
               YYYYMMDD, HHMMSS, ptrShForce, single_time_in_file, local_info, rc)
       if (rc<0) then
           call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error reading gfio file')
       endif 

       ! interpolate to tile space
       
       select case (MET_HINTERP)

       case (0)  ! nearest-neighbor interpolation
          
          do k=1,N_catd
             
             force_array(k,GEOSgcm_var) = ptrShForce(i1(k), j1(k))

          end do

       case (1)  ! bilinear interpolation
          
          supported_option_MET_HINTERP = .true. 
          
          do k=1,N_catd
             this_lon = tile_coord(k)%com_lon
             this_lat = tile_coord(k)%com_lat
             
             fnbr(1,1) = ptrShForce(i1(k),j1(k))
             fnbr(1,2) = ptrShForce(i1(k),j2(k))
             fnbr(2,1) = ptrShForce(i2(k),j1(k))
             fnbr(2,2) = ptrShForce(i2(k),j2(k))
             
             !DEC$ FORCEINLINE
             force_array(k,GEOSgcm_var) = BilinearInterpolation(this_lon, this_lat, &
                  x1(k), x2(k), y1(k), y2(k), fnbr, nodata_forcing, tol)
          end do

       case default

          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unsupported MET_HINTERP option')
          
       end select

       ! ----------------------------------------------    
       !
       ! For non-flux "S" forcing variables that are read from MERRA "tavg" files
       ! *optionally* read forward-looking 'tavg' file (if available) and interpolate 
       ! in time.
       !
       ! Doing so minimizes the time shift between the "true" (but unavailable) MERRA
       ! instantaneous forcing values and their off-line time interpolated equivalent
       ! -- at the expense of a dampened diurnal cycle. 
       !
       ! reichle+qliu,  8 Oct 2008

       minimize_shift = .true.
       
       ! minimize_shift should *only* affect "MERRA" forcing - reichle, 27 Feb 2012

       if  ((minimize_shift)                            .and.                    &
            (trim(GEOSgcm_defs(GEOSgcm_var,2))=='tavg') .and.                    &
            (trim(GEOSgcm_defs(GEOSgcm_var,5))=='S')           ) then
          
          date_time_tmp = date_time_fwd
          
          ! open file
          
          call get_GEOS_forcing_filename( fname_full,file_exists,date_time_tmp, daily_met_files,           &
               met_path_tmp, met_tag_tmp,                                        &
               GEOSgcm_defs(GEOSgcm_var,:), met_file_ext)

          call GEOS_openfile(FileOpenedHash,fname_full,fid,tile_coord,MET_HINTERP)

          !fid = ptrNode%fid

          if (fid>0) then

             i1=>local_info%i1
             i2=>local_info%i2
             j1=>local_info%j1
             j2=>local_info%j2
             x1=>local_info%x1
             x2=>local_info%x2
             y1=>local_info%y1
             y2=>local_info%y2
             
             YYYYMMDD = date_time_tmp%year*10000+date_time_tmp%month*100+date_time_tmp%day
             HHMMSS   = date_time_tmp%hour*10000+date_time_tmp%min*100  +date_time_tmp%sec
             
             ! read global gridded field of given variable
             
             call LDAS_GetVar( fid, trim(GEOSgcm_defs(GEOSgcm_var,1)),         &
               YYYYMMDD, HHMMSS, ptrShForce, .false.,local_info ,rc)
             
             if (rc<0) then
                err_msg = 'error reading gfio file'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if

             ! interpolate to tile space and in time
       
             select case (MET_HINTERP)
                
             case (0)  ! nearest-neighbor interpolation
                
                do k=1,N_catd

                   if ( abs(force_array(k,GEOSgcm_var) -nodata_GEOSgcm)<tol .or.          & 
                        abs(ptrShForce(i1(k),j1(k))    -nodata_GEOSgcm)<tol      ) then
                      
                      force_array(k,GEOSgcm_var) = nodata_GEOSgcm

                   else                
                      
                      force_array(k,GEOSgcm_var) =                                        &
                           0.5*(force_array(k,GEOSgcm_var) + ptrShForce(i1(k),j1(k)))

                   end if

                end do

             case (1)  ! bilinear interpolation

                do k=1,N_catd
                   if ( abs(force_array(k,GEOSgcm_var)-nodata_GEOSgcm)<tol .or.          & 
                        abs(ptrShForce(i1(k),j1(k))   -nodata_GEOSgcm)<tol       ) then
                      force_array(k,GEOSgcm_var) = nodata_GEOSgcm
                   else                
                      this_lon = tile_coord(k)%com_lon
                      this_lat = tile_coord(k)%com_lat
                      fnbr(1,1) = ptrShForce(i1(k),j1(k))
                      fnbr(1,2) = ptrShForce(i1(k),j2(k))
                      fnbr(2,1) = ptrShForce(i2(k),j1(k))
                      fnbr(2,2) = ptrShForce(i2(k),j2(k))
                      
                      !DEC$ FORCEINLINE
                      force_array(k,GEOSgcm_var) = 0.5*(force_array(k,GEOSgcm_var) +     &
                           BilinearInterpolation(this_lon, this_lat, x1(k), x2(k),       &
                           y1(k), y2(k), fnbr, nodata_forcing, tol))
                   end if
                end do
                
             case default
                
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unsupported MET_HINTERP option')
                
             end select
             
          end if    ! if (fid>0)
          
       end if       ! if (minimize_shift) .and. [...]
       
    end do          ! do GEOSgcm_var = 1,N_GEOSgcm_vars

    call FileOpenedHash%free( GEOS_closefile,.false. )

    deallocate(GEOSgcm_defs)

    ! --------------------------------------------------------------------
    
    ! convert variables and units of force_array to match met_force_type, 
    ! put into structure
    
    ! from GEOSgcm files:
    !
    !                     G5DAS     
    !                     M2INT      MERRA
    !                     M2COR
    !
    ! force_array(:, 1) = SWGDN      SWGDN   W/m2    (downward shortwave)     
    ! force_array(:, 2) = LWGAB      LWGAB   W/m2    ("absorbed" longwave)
    ! force_array(:, 3) = PARDR      PARDR   W/m2    (direct  PAR)
    ! force_array(:, 4) = PARDF      PARDF   W/m2    (diffuse PAR)
    ! force_array(:, 5) = PRECCU[*]  PRECTOT kg/m2/s (*see below*)
    ! force_array(:, 6) = PRECLS[*]  PRECCON kg/m2/s (*see below*)
    ! force_array(:, 7) = PRECSN[*]  PRECSNO kg/m2/s (*see below*)
    ! force_array(:, 8) = PS         PS      Pa      (surface air pressure)      
    ! force_array(:, 9) = HLML       HLML    m       (height of lowest model level "LML")
    ! force_array(:,10) = TLML       TLML    K       (air temperature   at LML)
    ! force_array(:,11) = QLML       QLML    kg/kg   (air spec humidity at LML)
    ! force_array(:,12) = SPEEDLML   ULML    m/s     (wind speed/U-wind at LML)
    ! force_array(:,13) = n/a        VLML    m/s     (           V-wind at LML)
    !
    !  PRECTOT           kg/m2/s  (total       rain+snow) = PRECCU+PRECLS+PRECSNO
    !  PRECCON           kg/m2/s  (convective  rain+snow)
    !  PRECCU            kg/m2/s  (convective  rain)
    !  PRECLS            kg/m2/s  (large-scale rain)
    !  PRECSNO           kg/m2/s  (total       snow)
    
    met_force_new%SWdown    = force_array(:, 1)
    met_force_new%LWdown    = force_array(:, 2)
    met_force_new%PARdrct   = force_array(:, 3)
    met_force_new%PARdffs   = force_array(:, 4)
    
    met_force_new%Psurf     = force_array(:, 8)
    
    met_force_new%RefH      = force_array(:, 9)
    met_force_new%Tair      = force_array(:,10)
    met_force_new%Qair      = force_array(:,11)
    
    do k=1,N_catd

       ! get wind speed       
       
       if (MERRA_file_specs) then
          
          if ( abs(force_array(k,12)-nodata_GEOSgcm)<tol .or.           &
               abs(force_array(k,13)-nodata_GEOSgcm)<tol      ) then
             
             met_force_new(k)%Wind = nodata_forcing
             
          else
             
             met_force_new(k)%Wind = &
                  sqrt( force_array(k,12)**2 + force_array(k,13)**2 )
             
          end if
          
       else  ! G5DAS file specs
          
          met_force_new(k)%Wind = force_array(k,12)
          
       end if
       
       
       ! rainfall
       
       if ( (abs(force_array(k,5)-nodata_GEOSgcm)<tol) .or.               &
            (abs(force_array(k,6)-nodata_GEOSgcm)<tol) .or.               &
            (abs(force_array(k,7)-nodata_GEOSgcm)<tol)       )  then
          
          met_force_new(k)%Rainf   = nodata_forcing
          met_force_new(k)%Rainf_C = nodata_forcing
          met_force_new(k)%Snowf   = nodata_forcing
          
       else
          
          if (MERRA_file_specs) then
             
             if (force_array(k,5)>0) then
                
                met_force_new(k)%Snowf = force_array(k,7)

                ! total_rain = total_precip - total_snow
                
                met_force_new(k)%Rainf = force_array(k,5) - force_array(k,7)
                
                ! conv_rain = (conv_precip/total_precip) * total_rain
                
                met_force_new(k)%Rainf_C = &
                     force_array(k,6)/force_array(k,5)*met_force_new(k)%Rainf
          
             else   
                
                met_force_new(k)%Rainf   = 0.
                met_force_new(k)%Rainf_C = 0.
                met_force_new(k)%Snowf   = 0.
                
             end if
             
          else
             
             ! G5DAS file specs
             
             met_force_new(k)%Rainf   = force_array(k,5)+force_array(k,6)
             met_force_new(k)%Rainf_C = force_array(k,5)
             met_force_new(k)%Snowf   = force_array(k,7)
             
          end if
          
       end if

    end do

    if(AEROSOL_DEPOSITION /=0) then
       met_force_new%DUDP001   = force_array(:,13)
       met_force_new%DUDP002   = force_array(:,14)
       met_force_new%DUDP003   = force_array(:,15)
       met_force_new%DUDP004   = force_array(:,16)
       met_force_new%DUDP005   = force_array(:,17)
       met_force_new%DUSV001   = force_array(:,18)
       met_force_new%DUSV002   = force_array(:,19)
       met_force_new%DUSV003   = force_array(:,20)
       met_force_new%DUSV004   = force_array(:,21)
       met_force_new%DUSV005   = force_array(:,22)
       met_force_new%DUWT001   = force_array(:,23)
       met_force_new%DUWT002   = force_array(:,24)
       met_force_new%DUWT003   = force_array(:,25)
       met_force_new%DUWT004   = force_array(:,26)
       met_force_new%DUWT005   = force_array(:,27)
       met_force_new%DUSD001   = force_array(:,28)
       met_force_new%DUSD002   = force_array(:,29)
       met_force_new%DUSD003   = force_array(:,30)
       met_force_new%DUSD004   = force_array(:,31)
       met_force_new%DUSD005   = force_array(:,32)
       met_force_new%BCDP001   = force_array(:,33)
       met_force_new%BCDP002   = force_array(:,34)
       met_force_new%BCSV001   = force_array(:,35)
       met_force_new%BCSV002   = force_array(:,36)
       met_force_new%BCWT001   = force_array(:,37)
       met_force_new%BCWT002   = force_array(:,38)
       met_force_new%BCSD001   = force_array(:,39)
       met_force_new%BCSD002   = force_array(:,40)
       met_force_new%OCDP001   = force_array(:,41)
       met_force_new%OCDP002   = force_array(:,42)
       met_force_new%OCSV001   = force_array(:,43)
       met_force_new%OCSV002   = force_array(:,44)
       met_force_new%OCWT001   = force_array(:,45)
       met_force_new%OCWT002   = force_array(:,46)
       met_force_new%OCSD001   = force_array(:,47)
       met_force_new%OCSD002   = force_array(:,48)
       met_force_new%SUDP003   = force_array(:,49)
       met_force_new%SUSV003   = force_array(:,50)
       met_force_new%SUWT003   = force_array(:,51)
       met_force_new%SUSD003   = force_array(:,52)
       met_force_new%SSDP001   = force_array(:,53)
       met_force_new%SSDP002   = force_array(:,54)       
       met_force_new%SSDP003   = force_array(:,55)
       met_force_new%SSDP004   = force_array(:,56)
       met_force_new%SSDP005   = force_array(:,57)
       met_force_new%SSSV001   = force_array(:,58)
       met_force_new%SSSV002   = force_array(:,59)
       met_force_new%SSSV003   = force_array(:,60)
       met_force_new%SSSV004   = force_array(:,61)
       met_force_new%SSSV005   = force_array(:,62)
       met_force_new%SSWT001   = force_array(:,63)
       met_force_new%SSWT002   = force_array(:,64)
       met_force_new%SSWT003   = force_array(:,65)
       met_force_new%SSWT004   = force_array(:,66)
       met_force_new%SSWT005   = force_array(:,67)
       met_force_new%SSSD001   = force_array(:,68)
       met_force_new%SSSD002   = force_array(:,69)
       met_force_new%SSSD003   = force_array(:,70)
       met_force_new%SSSD004   = force_array(:,71)
       met_force_new%SSSD005   = force_array(:,72)
    endif                                      
    
    deallocate(force_array)
    
  end subroutine get_GEOS
 
! ******************************************************************
  subroutine LDAS_GetVar(fid, vname, yyyymmdd, hhmmss, &
                         ptrShForce,single_time_in_file,local_info, rc)
     use netcdf
     implicit none
     include 'mpif.h'

     integer,intent(in)           ::  fid                 ! File handle
     character(len=*), intent(in) ::  vname               ! Variable name
     integer, intent(in)          ::  yyyymmdd            ! Year-month-day, e.g., 19971003
     integer,intent(in)           ::  hhmmss              ! Hour-minute-second, e.g., 120000
     logical,intent(in)           ::  single_time_in_file ! if true, no time index is necessary
     type(local_grid),intent(in)  ::  local_info
     !OUTPUT PARAMETERS:
     real,pointer,intent(inout)   ::  ptrShForce(:,:)     ! Gridded data read for this time
     integer,intent(out)          ::  rc

     ! local
     integer                      :: begDate, begTime, seconds, minutes, incSecs
     integer                      :: iistart(3), iicount(3), timeIndex
     integer                      :: istart(4), icount(4)     ! cs grid
     real,         pointer        :: tmpShared(:,:,:,:) ! cs grid
     type(c_ptr)                  ::   c_address
     integer                      :: nv_id,imin, jmin, imax, jmax,ierr
     integer                      :: DiffDate
     integer                      :: status
     character(*), parameter      :: Iam="LDAS_getvar" 
     logical                      :: isCubed

     rc = 0
     isCubed = .false.
     if(local_info%N_lat == 6*local_info%N_lon) then
        isCubed = .true.
        istart = 1
        icount(1) = local_info%N_lon
        icount(2) = local_info%N_lon
        icount(3) = 6
        icount(4) = 1
     else
        iistart = 1
        iicount(1) = local_info%N_lon
        iicount(2) = local_info%N_lat
        iicount(3) = 1
     endif

     if (.not. single_time_in_file ) then   ! determine start index
        call GetBegDateTime ( fid, begDate, begTime, incSecs, rc )
        if (rc .NE. 0) then
           print* ,"LDAS_GetVar: could not determine begin_date/begin_time"
           return
        endif
        seconds = DiffDate (begDate, begTime, yyyymmdd, hhmmss)
        ! Make sure input time are valid (assume time is not periodic)
        if (seconds .LT. 0) then
           print *, 'LDAS_GetVar: Error code from diffdate.  Problem with date/time.'
           rc = -7
           return
        endif

        if ( MOD (seconds,60) .eq. 0 ) then
           minutes = seconds / 60
        else
           print *, 'LDAS_GetVar: Currently, times must fall on minute boundaries.'
           rc = -6
           return
        endif

        ! Determine the time index from the offset and time increment.
        if ( MOD (seconds, incSecs) .ne. 0 ) then
           print *, 'GFIO_getvar: Absolute time of ',seconds,' not ',  &
                'possible with an interval of ',incSecs
           rc = -2
           return
        else
           timeIndex = seconds/incSecs + 1
        endif
        iistart(3)=timeIndex
        istart(4) =timeIndex
     endif
     ! node root read and share
     call MAPL_SyncSharedMemory(rc=status)

     if (MAPL_AmNodeRoot .or. (.not. MAPL_ShmInitialized)) then
        rc= NF90_INQ_VARID( fid, vname, nv_id)
        _ASSERT( rc == nf90_noerr, "nf90 error")
        if (isCubed) then
          c_address = c_loc(ptrShForce(1,1))
          call c_f_pointer(c_address,tmpShared,shape=icount)
          rc= NF90_GET_VAR( fid, nv_id, tmpShared, start=istart,count=icount) 
        else
          rc= NF90_GET_VAR( fid, nv_id, ptrShForce, start=iistart,count=iicount) 
        endif
        _ASSERT( rc == nf90_noerr, "nf90 error")
     endif

     call MAPL_SyncSharedMemory(rc=status)

  end subroutine LDAS_GetVar

  ! ****************************************************************
  
  function BilinearInterpolation(x,y,x1,x2,y1,y2,fnbr,UNDEF,tol) result(fxy)

    ! pchakrab: function to compute (bilinear) interpolated
    ! value f(x,y). If we know the value at 4 points
    ! f11=f(x1,y1), f12=f(x1,y2), f21=f(x2,y1) and f22=f(x2,y2),
    ! the interpolated value, fxy, is computed as
    !
    ! f11(x2-x)(y2-y) + f21(x-x1)(y2-y) + f12(x2-x)(y-y1) + f22(x-x1)(y-y1)
    ! ---------------------------------------------------------------------
    !                       (x2-x1)(y2-y1)
    !
    ! NOTE 1: UNDEF is the nodata value (1e15)
    ! NOTE 2: If all the neighbouring f values are UNDEF, fxy = UNDEF
    !         else, the UNDEF value at a corner is replaced by an
    !         average of non-UNDEF values before fxy is computed

    real, intent(in) :: x, y, x1, x2, y1, y2, fnbr(2,2)
    real, intent(in) :: UNDEF, tol
    real :: fxy ! output

    ! local
    real :: floc(2,2), fsum, dx1, dx2, dy1, dy2
    logical :: NoData(2,2)
    integer :: i, j, numNoData

    floc = fnbr

    ! check for nodata values
    ! and compute the sum of defined f values
    NoData = .false.
    numNodata = 0
    fsum = 0.0
    do j=1,2
       do i=1,2
          if (abs(floc(i,j)-UNDEF)<tol) then
             NoData(i,j) = .true.
             numNodata = numNodata + 1
          else
             fsum = fsum + floc(i,j)
          end if
       end do
    end do
 
    ! 0<=numNodata<=4
    if (numNodata==4) then
       fxy = UNDEF
    else 
       if (numNodata>0) then
          where (NoData) floc = fsum/real(4-numNoData)
       end if
       dx1 = x-x1
       dx2 = x2-x
       dy1 = y-y1
       dy2 = y2-y
       fxy = (                                       &
            floc(1,1)*dx2*dy2 + floc(2,1)*dx1*dy2 +  &
            floc(1,2)*dx2*dy1 + floc(2,2)*dx1*dy1    &
            ) / ((x2-x1)*(y2-y1))
    end if

  end function BilinearInterpolation

  ! ****************************************************************

  subroutine parse_MERRA_met_tag( met_path_in, met_tag_in, date_time, &
       met_path_default, met_path_prec, met_tag_out, use_prec_corr )
    
    ! reichle, 1 Dec 2009
    
    ! parse MERRA "met_tag", extract MERRA stream, assemble data paths
    !
    ! Convention for driver_inputs*nml file or command line arguments:
    !
    !   met_path = "/*/merra_land/"
    !
    ! eg, on discover:  /gpfsm/dnb51/projects/p15/iau/merra_land/
    !     on discover:  /discover/nobackup/qliu/merra_land/
    !     on land01:    /merra_land/
    !
    !
    !   met_tag  = "d5_merra_[STREAM]__[GCM-TAG]{__prec[PREC]}"
    !
    ! where {__prec[PREC]} is optional and where
    !
    !       STREAM  = 'jan79', 'jan89', 'jan98', or 'cross'
    !       GCM-TAG = 'GEOSdas-2_1_4', ...
    !       PREC    = 'CMAPvS', 'GPCPv1.1', ...
    !
    ! examples:
    !
    !   STREAM  = 'jan89'            : use only Stream 2 MERRA data
    !   STREAM  = 'cross'            : integrate across more than one stream
    !   GCM-TAG = 'GEOSdas-2_1_4'    : tag of standard MERRA (may differ for replays)
    !   PREC                         : identifier for precip corrections to MERRA forcing
    !
    ! ---------------------------------------------------------------------------    

    implicit none
    
    character(*), intent(in)  :: met_path_in
    character(*), intent(in)  :: met_tag_in
    
    type(date_time_type), intent(in) :: date_time
    
    character(200), intent(out) :: met_path_default, met_path_prec
    character( 80), intent(out) :: met_tag_out
    
    logical,        intent(out) :: use_prec_corr
    
    ! local variables

    integer :: is, ie
    
    type(date_time_type) :: dt1, dt2
    
    character( 5) :: stream
    character(80) :: gcm_tag, prec_tag

    character(len=*), parameter :: Iam = 'parse_MERRA_met_tag'
    character(len=400) :: err_msg
    
    ! ----------------------------------------------------------
    
    ! define intervals that determine which MERRA stream is used
    ! in integrations that "cross" multiple streams
    
    ! Stream 1 ('jan79') -->  16 Dec 1978 - 31 Dec 1989
    ! Stream 2 ('jan89') -->   1 Jan 1990 - 31 Dec 1998
    ! Stream 3 ('jan98') -->   1 Jan 1999 - present
    
    ! dates before dt1 use Stream 1
    
    dt1%year  = 1993
    dt1%month = 1
    dt1%day   = 1
    dt1%hour  = 0
    dt1%min   = 0
    dt1%sec   = 0
    
    ! otherwise, dates before dt2 use Stream 2
    
    dt2%year  = 2001
    dt2%month = 1
    dt2%day   = 1
    dt2%hour  = 0
    dt2%min   = 0
    dt2%sec   = 0
    
    ! ----------------------------------------------------
    
    ! initialize
    
    met_tag_out = repeat(' ', len(met_tag_out))
    
    ! define which stream to use
    
    if (met_tag_in(10:14)=='cross') then
       
       if     (datetime_lt_refdatetime( date_time, dt1 )) then
          
          stream = 'jan79'
          
       elseif (datetime_lt_refdatetime( date_time, dt2 )) then
          
          stream = 'jan89'

       else
          
          stream = 'jan98'
          
       end if
       
       met_tag_out = met_tag_in(1:9) // stream
       
    else
       
       met_tag_out = met_tag_in(1:14)       
       
    end if
    
    ! -----------------------------------------------------
    !
    ! identify GCM tag and which precip corrections to use, 
    ! assemble met_path accordingly
    !
    !   met_tag = "d5_merra_[STREAM]__[GCM-TAG]{__prec[PREC]}"
    !
    ! where {__prec[PREC]} is optional

    is = index( met_tag_in, '__')
    ie = index( met_tag_in, '__', .true.)
    
    if (is/=ie) then    ! using precip corrections
       
       gcm_tag  = met_tag_in(is+2:ie-1)

       met_path_default = trim(met_path_in) // '/MERRA_land_forcing/' // &
            trim(gcm_tag) // '/'

       prec_tag = met_tag_in(ie+6:len(met_tag_in))

       met_path_prec    = trim(met_path_in) // '/precip_corr_' // trim(prec_tag) // &
            '/' //  trim(gcm_tag) // '/'
       
       use_prec_corr = .true.
       
    else                ! not using precip corrections

       gcm_tag  = met_tag_in(is+2:len(met_tag_in))
       
       met_path_default = trim(met_path_in) // '/MERRA_land_forcing/' // &
            trim(gcm_tag) // '/'
       
       prec_tag = repeat(' ', len(prec_tag))
       
       met_path_prec = met_path_default
       
       use_prec_corr = .false.

       ! check if prec_tag was accidentally appended with a single underscore
       
       if (len_trim(met_tag_in)>29) then
          
          err_msg = 'questionable met_tag_in, not enough double underscores'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
   
       end if

    end if

  end subroutine parse_MERRA_met_tag
  
  ! ****************************************************************

  subroutine parse_MERRA2_met_tag( met_path_in, met_tag_in, date_time, &
       met_path_default, met_path_prec, met_tag_out, use_prec_corr )
    
    ! reichle, 27 Jul 2015
    
    ! parse MERRA2 "met_tag", extract MERRA stream, assemble data paths
    !
    !   met_tag  = "M2[xxx]_[STREAM]{__prec[PREC]}"
    !
    ! where {__prec[PREC]} is optional and where
    !
    !       [xxx]   = 'GCM' or 'COR'
    !       STREAM  = '100', '200', '300', '400', or 'cross'
    !       PREC    = 'CMAPvS', 'GPCPv1.1', ...
    !
    ! examples:
    !
    !   STREAM  = '200'              : use only Stream 2 MERRA data
    !   STREAM  = 'cross'            : integrate across more than one stream
    !   PREC                         : identifier for corrected precip data
    !
    ! ---------------------------------------------------------------------------    

    implicit none
    
    character(*), intent(in)  :: met_path_in
    character(*), intent(in)  :: met_tag_in
    
    type(date_time_type), intent(in) :: date_time
    
    character(200), intent(out) :: met_path_default, met_path_prec
    character( 80), intent(out) :: met_tag_out
    
    logical,        intent(out) :: use_prec_corr
    
    ! local variables

    integer :: is
    
    type(date_time_type) :: dt1, dt2, dt3
    
    character(10) :: stream
    character(80) :: prec_tag
    
    character(len=*), parameter :: Iam = 'parse_MERRA2_met_tag'
    character(len=400) :: err_msg

    ! ----------------------------------------------------------
    
    ! define intervals that determine which MERRA-2 stream is used
    ! in integrations that "cross" multiple streams
    !    
    ! 1/1/1980 - 12/31/1991:  MERRA2_100 (Stream 1)
    ! 1/1/1992 - 12/31/2000:  MERRA2_200 (Stream 2)
    ! 1/1/2001 - 12/31/2010:  MERRA2_300 (Stream 3)
    ! 1/1/2011 - present:     MERRA2_400 (Stream 4)
    
    ! dates before dt1 use Stream 1
    
    dt1%year  = 1992
    dt1%month = 1
    dt1%day   = 1
    dt1%hour  = 0
    dt1%min   = 0
    dt1%sec   = 0
    
    ! otherwise, dates before dt2 use Stream 2
    
    dt2%year  = 2001
    dt2%month = 1
    dt2%day   = 1
    dt2%hour  = 0
    dt2%min   = 0
    dt2%sec   = 0

    ! otherwise, dates before dt3 use Stream 3
    
    dt3%year  = 2011
    dt3%month = 1
    dt3%day   = 1
    dt3%hour  = 0
    dt3%min   = 0
    dt3%sec   = 0

    ! ----------------------------------------------------
    
    ! initialize
    
    met_tag_out = repeat(' ', len(met_tag_out))
    
    stream      = repeat(' ', len(stream     ))
    
    ! define which stream to use
    
    if (met_tag_in(7:11)=='cross') then
       
       if     (datetime_lt_refdatetime( date_time, dt1 )) then
          
          stream = 'MERRA2_100'
          
       elseif (datetime_lt_refdatetime( date_time, dt2 )) then
          
          stream = 'MERRA2_200'

       elseif (datetime_lt_refdatetime( date_time, dt3 )) then
          
          stream = 'MERRA2_300'

       else

          stream = 'MERRA2_400'
          
       end if
       
       met_tag_out = trim(stream)
       
    else
       
       met_tag_out = 'MERRA2_' // met_tag_in(7:9)       
       
    end if
    
    met_path_default = trim(met_path_in) // '/' 
    
    ! -----------------------------------------------------
    !
    ! identify which precip corrections to use, 
    ! assemble met_path accordingly
    !
    !   met_tag = "M2[xxx]_[STREAM]{__prec[PREC]}"
    !
    ! where {__prec[PREC]} is optional

    is = index( met_tag_in, '__')
    
    if (is>0) then    ! using precip corrections
       
       prec_tag      = met_tag_in(is+6:len(met_tag_in))
       
       met_path_prec = trim(met_path_in) // '/precip_corr_' // trim(prec_tag) // '/'
       
       use_prec_corr = .true.
       
    else              ! not using precip corrections
       
       prec_tag      = repeat(' ', len(prec_tag))
       
       met_path_prec = met_path_default
       
       use_prec_corr = .false.

       ! check if prec_tag was accidentally appended with a single underscore
       
       if (len_trim(met_tag_in)>11) then
          
          err_msg = 'questionable met_tag_in, not enough double underscores'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
   
       end if
       
    end if

  end subroutine parse_MERRA2_met_tag
  
  ! ****************************************************************
  
  subroutine parse_G5DAS_met_tag( met_path_in, met_tag_in, date_time, &
       met_path_default, met_path_prec, met_tag_out, use_prec_corr,   &
       use_Predictor )
    
    ! parse G5DAS "met_tag"
    !
    ! Convention for driver_inputs*nml file or command line arguments:
    !
    !   met_tag  = "[G5DAS-NAME]{__prec[PREC]}{__Nx+-}"
    !
    ! where {__prec[PREC]} and {__Nx+-} are optional and where
    !
    !  G5DAS-NAME : name of standard G5DAS forcing (must not contain "__"!)
    !               e.g., "e5110_fp", "d591_rpit1", "d591_fpit", ...
    !               for cross-stream forcing, use "cross_d5124_RPFPIT" or "cross_FP" 
    !  PREC       : identifier for precip corrections to G5DAS forcing (eg., 'GPCPv1.1')
    !  
    ! If {__Nx+-} is present, set flag for use forcing files from the DAS/GCM Predictor 
    !  segment. (Default is to use forcing from Corrector segment.)
    !
    ! reichle,  3 Jun 2013
    ! reichle, 20 Sep 2013: restructured, added "use_Predictor"
    ! reichle, 23 Sep 2013: added "cross-stream" capability
    ! reichle,  5 Dec 2014: updated "cross-stream" dates
    ! reichle, 14 Sep 2015: revised FP "cross-stream" dates
    ! reichle, 28 Dec 2016: added 5.12.4 RPIT/FPIT streams
    ! reichle, 19 Jan 2017: added FP transition from e5131 to f516
    ! reichle, 30 Oct 2017: added FP transition from f516 to f517
    ! reichle,  9 Jul 2018: added FP transition from f517 to f521
    ! reichle, 10 Oct 2019: added FP transition from f521 to f522
    ! reichle, 17 Jan 2020: added FP transition from f522 to f525
    ! reichle,  3 Apr 2020: added FP transition from f525 to f525_p5
    !    
    ! ---------------------------------------------------------------------------    

    implicit none

    character(*),       intent(in)  :: met_path_in
    character(*),       intent(in)  :: met_tag_in

    type(date_time_type), intent(in)  :: date_time

    character(200),       intent(out) :: met_path_default, met_path_prec
    character( 80),       intent(out) :: met_tag_out
    
    logical,              intent(out) :: use_prec_corr
    logical,              intent(out) :: use_Predictor

    ! local variables

    integer                     :: is, ie, ii, N_opt_tag
    
    character(80)               :: prec_tag, stream
    character(80), dimension(2) :: tmp_tag

    type(date_time_type)        :: dt_end_d591_rpit1
    type(date_time_type)        :: dt_end_d591_rpit2
    type(date_time_type)        :: dt_end_d591_rpit3
    type(date_time_type)        :: dt_end_d591_fpit

    type(date_time_type)        :: dt_end_d5124_rpit1
    type(date_time_type)        :: dt_end_d5124_rpit2
    type(date_time_type)        :: dt_end_d5124_rpit3

    type(date_time_type)        :: dt_end_e5110_fp
    type(date_time_type)        :: dt_end_e5130_fp
    type(date_time_type)        :: dt_end_e5131_fp
    type(date_time_type)        :: dt_end_f516_fp
    type(date_time_type)        :: dt_end_f517_fp   
    type(date_time_type)        :: dt_end_f521_fp   
    type(date_time_type)        :: dt_end_f522_fp
    type(date_time_type)        :: dt_end_f525_fp
    type(date_time_type)        :: dt_end_f525_p5_fp

    character(len=*), parameter :: Iam = 'parse_G5DAS_met_tag'
    character(len=400) :: err_msg

    ! ----------------------------------------------------------
    !
    ! define transition times between RP-IT, FP-IT or FP streams 
    ! if "cross-stream" forcing is requested
    !
    !            | stream start |  stream end (as of 5 Dec 2014)
    ! ----------------------------------------
    ! d591_rpit1 |  1 Jan 2000  |   1 Jun 2007          
    ! d591_rpit2 |  1 Jun 2006  |  30 Dec 2011          
    ! d591_rpit3 |  1 Jan 2011  |   6 May 2013
    ! d591_fpit  | 30 Dec 2012  |   (present) 

    dt_end_d591_rpit1%year      = 2007
    dt_end_d591_rpit1%month     = 6
    dt_end_d591_rpit1%day       = 1
    dt_end_d591_rpit1%hour      = 0
    dt_end_d591_rpit1%min       = 0
    dt_end_d591_rpit1%sec       = 0

    dt_end_d591_rpit2%year      = 2011
    dt_end_d591_rpit2%month     = 12
    dt_end_d591_rpit2%day       = 1
    dt_end_d591_rpit2%hour      = 0
    dt_end_d591_rpit2%min       = 0
    dt_end_d591_rpit2%sec       = 0

    dt_end_d591_rpit3%year      = 2013
    dt_end_d591_rpit3%month     = 5
    dt_end_d591_rpit3%day       = 1
    dt_end_d591_rpit3%hour      = 0
    dt_end_d591_rpit3%min       = 0
    dt_end_d591_rpit3%sec       = 0
    
    dt_end_d591_fpit%year       = 9999
    dt_end_d591_fpit%month      = 1
    dt_end_d591_fpit%day        = 1
    dt_end_d591_fpit%hour       = 0
    dt_end_d591_fpit%min        = 0
    dt_end_d591_fpit%sec        = 0

    !                  | stream start |  stream end (as of 28 Dec 2016)
    ! ----------------------------------------
    ! d5124_rpit1      |  1 Jan 2000  |   1 Jan 2004          
    ! d5124_rpit2      |  1 Jan 2004  |   1 Jan 2012          
    ! d5124_rpit3      |  1 Jan 2012  |   (present)

    dt_end_d5124_rpit1%year      = 2004
    dt_end_d5124_rpit1%month     = 1
    dt_end_d5124_rpit1%day       = 1
    dt_end_d5124_rpit1%hour      = 0
    dt_end_d5124_rpit1%min       = 0
    dt_end_d5124_rpit1%sec       = 0

    dt_end_d5124_rpit2%year      = 2012
    dt_end_d5124_rpit2%month     = 1
    dt_end_d5124_rpit2%day       = 1
    dt_end_d5124_rpit2%hour      = 0
    dt_end_d5124_rpit2%min       = 0
    dt_end_d5124_rpit2%sec       = 0

    dt_end_d5124_rpit3%year      = 9999
    dt_end_d5124_rpit3%month     = 1
    dt_end_d5124_rpit3%day       = 1
    dt_end_d5124_rpit3%hour      = 0
    dt_end_d5124_rpit3%min       = 0
    dt_end_d5124_rpit3%sec       = 0

    ! ---------------------------------        
    !
    ! FP streams
    !
    ! Stream start/end in terms of availability 
    ! in the GEOS-5 archive (approximately):  
    !
    !            | stream start | stream end 
    ! ----------------------------------------
    ! e5110_fp   | 11 Jun 2013  | 19 Aug 2014 
    ! e5130_fp   |  1 Aug 2014  |  4 May 2015
    ! e5131_fp   |  1 May 2015  | 24 Jan 2017
    ! f516_fp    | 24 Jan 2017  |  1 Nov 2017
    ! f517_fp    |  1 Nov 2017  | 11 Jul 2018
    ! f521_fp    | 11 Jul 2018  | 13 Mar 2019
    ! f522_fp    | 13 Mar 2019  | 30 Jan 2020
    ! f525_fp    | 30 Jan 2020  |  7 Apr 2020
    ! f525_p5_fp |  7 Apr 2020  |   (present)
    !
    ! Official stream transition times (as defined
    ! by GMAO ops group) are:
    !
    ! FP e5110 --> e5130   : 20 Aug 2014, 6z ADAS analysis
    ! FP e5130 --> e5131   :  1 May 2015, 6z ADAS analysis
    ! FP e5131 --> f516    : 24 Jan 2017, 6z ADAS analysis
    ! FP f516  --> f517    :  1 Nov 2017, 6z ADAS analysis
    ! FP f517  --> f521    : 11 Jul 2018, 6z ADAS analysis
    ! FP f521  --> f522    : 13 Mar 2019, 6z ADAS analysis
    ! FP f522  --> f525    : 30 Jan 2020, 6z ADAS analysis
    ! FP f525  --> f525_p5 :  7 Apr 2020, 6z ADAS analysis
    !
    ! Official stream transition times refer to the definition
    ! of the official FP files with generic file names on the 
    ! NCCS data portal. 
    !
    ! Note that "lfo" files for the 6z analysis start from 4z, 
    ! that is, the exact transition time is slightly different
    ! from the LDAS perspective.
    !
    ! Define LDASsa "cross" streams dates for 0z of the 
    ! first day that contains only "new" data for 
    ! compatibility with the SMAP L4 Ops system and because of the
    ! asymmetric availability of the data on the NCCS data portal.
    !
    ! - reichle, 14 Sep 2015
    !
    ! Revised Aug 2014 cross-over date to Aug 20 at 0z because 
    ! e5110 data are not available for Aug 20
    ! - reichle+qliu, 29 Jan 2016

    dt_end_e5110_fp%year   = 2014
    dt_end_e5110_fp%month  = 8
    dt_end_e5110_fp%day    = 20
    dt_end_e5110_fp%hour   = 0
    dt_end_e5110_fp%min    = 0
    dt_end_e5110_fp%sec    = 0  

    dt_end_e5130_fp%year   = 2015
    dt_end_e5130_fp%month  = 5
    dt_end_e5130_fp%day    = 2
    dt_end_e5130_fp%hour   = 0
    dt_end_e5130_fp%min    = 0
    dt_end_e5130_fp%sec    = 0  

    dt_end_e5131_fp%year   = 2017
    dt_end_e5131_fp%month  = 1
    dt_end_e5131_fp%day    = 24
    dt_end_e5131_fp%hour   = 3
    dt_end_e5131_fp%min    = 0
    dt_end_e5131_fp%sec    = 0

    dt_end_f516_fp%year    = 2017
    dt_end_f516_fp%month   = 11
    dt_end_f516_fp%day     = 1
    dt_end_f516_fp%hour    = 3
    dt_end_f516_fp%min     = 0
    dt_end_f516_fp%sec     = 0

    dt_end_f517_fp%year    = 2018
    dt_end_f517_fp%month   = 7
    dt_end_f517_fp%day     = 11
    dt_end_f517_fp%hour    = 3
    dt_end_f517_fp%min     = 0
    dt_end_f517_fp%sec     = 0  

    dt_end_f521_fp%year    = 2019
    dt_end_f521_fp%month   = 3
    dt_end_f521_fp%day     = 13
    dt_end_f521_fp%hour    = 3
    dt_end_f521_fp%min     = 0
    dt_end_f521_fp%sec     = 0

    dt_end_f522_fp%year    = 2020
    dt_end_f522_fp%month   = 1
    dt_end_f522_fp%day     = 30
    dt_end_f522_fp%hour    = 3
    dt_end_f522_fp%min     = 0
    dt_end_f522_fp%sec     = 0  

    dt_end_f525_fp%year    = 2020
    dt_end_f525_fp%month   = 4
    dt_end_f525_fp%day     = 7
    dt_end_f525_fp%hour    = 3
    dt_end_f525_fp%min     = 0
    dt_end_f525_fp%sec     = 0  

    dt_end_f525_p5_fp%year = 9999
    dt_end_f525_p5_fp%month= 1
    dt_end_f525_p5_fp%day  = 1
    dt_end_f525_p5_fp%hour = 0
    dt_end_f525_p5_fp%min  = 0
    dt_end_f525_p5_fp%sec  = 0  

    ! ----------------------------------------------------

    ! initialize
    
    met_tag_out   = repeat(' ', len(met_tag_out))

    stream        = repeat(' ', len(stream     ))

    use_prec_corr = .false.

    use_Predictor = .false.

    ! -----------------------------------------------------
    !
    ! identify "GCM" tag and whether to use precip corrections and/or
    ! forcing from the DAS/GCM Predictor segment (default: Corrector segment)
    !
    !   met_tag  = "[G5DAS-NAME]{__prec[PREC]}{__Nx+-}"
    !
    ! where {__prec[PREC]} and {__Nx+-} are optional
    
    is = index( met_tag_in, '__')
    ie = index( met_tag_in, '__', .true.)
    
    ! determine how many optional tag segments are present
    
    if     (is==0) then    
       
       N_opt_tag   = 0
       
       met_tag_out = met_tag_in
       
    elseif (is==ie) then
       
       N_opt_tag   = 1
       
       met_tag_out = met_tag_in(1:(is-1))
       
       tmp_tag(1)  = met_tag_in((is+2):len(met_tag_in))
       
    else

       ! make sure there are at most two optional tag segments

       ! "ii" should be the index for the next "__" after "is"
       
       ii = index( met_tag_in(is+2:len(met_tag_in)), '__')
       
       if (is+1+ii==ie) then
          
          N_opt_tag  = 2
          
          met_tag_out = met_tag_in(1:(is-1))
          
          tmp_tag(1)  = met_tag_in((is+2):ie-1)
          
          tmp_tag(2)  = met_tag_in((ie+2):len(met_tag_in))
          
       else
          
          err_msg = 'invalid met_tag_in, too many double underscores'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                    
       end if

    end if


    ! resolve cross-stream if requested
    
    if     (met_tag_out(1:17)=='cross_d591_RPFPIT') then
       
       if     (datetime_lt_refdatetime( date_time, dt_end_d591_rpit1 )) then
          
          stream = 'd591_rpit1_jan00'   ! use d591 RP-IT stream 1
          
       elseif (datetime_lt_refdatetime( date_time, dt_end_d591_rpit2 )) then
          
          stream = 'd591_rpit2_jun06'   ! use d591 RP-IT stream 2

       elseif (datetime_lt_refdatetime( date_time, dt_end_d591_rpit3 )) then
          
          stream = 'd591_rpit3_jan11'   ! use d591 RP-IT stream 3
          
       else
          
          stream = 'd591_fpit'          ! use d591 FP-IT 
          
       end if

    elseif (met_tag_out(1:18)=='cross_d5124_RPFPIT') then

       if     (datetime_lt_refdatetime( date_time, dt_end_d5124_rpit1 )) then

          stream = 'd5124_rpit_jan00'  ! use d5124 RP-IT stream 1

       elseif (datetime_lt_refdatetime( date_time, dt_end_d5124_rpit2 )) then

          stream = 'd5124_rpit_jan04'  ! use d5124 RP-IT stream 2

       else

          stream = 'd5124_rpit_jan12'  ! use d5124 RP-IT stream 3 

       end if

    elseif (met_tag_out(1:8)=='cross_FP') then
       
       if     (datetime_lt_refdatetime( date_time, dt_end_e5110_fp )) then
          
          stream = 'e5110_fp'           ! use GEOS-5.11.0 output

       elseif (datetime_lt_refdatetime( date_time, dt_end_e5130_fp )) then
          
          stream = 'e5130_fp'           ! use GEOS-5.13.0 output

       elseif (datetime_le_refdatetime( date_time, dt_end_e5131_fp )) then

          ! Note "less-than-or-equal" (_le_) above

          stream = 'e5131_fp'           ! use GEOS-5.13.1 output

       elseif (datetime_le_refdatetime( date_time, dt_end_f516_fp )) then

          ! Note "less-than-or-equal" (_le_) above

          stream = 'f516_fp'            ! use GEOS-5.16.x output

       elseif (datetime_le_refdatetime( date_time, dt_end_f517_fp )) then
          
          ! Note "less-than-or-equal" (_le_) above
          
          stream = 'f517_fp'            ! use GEOS-5.17.x output

       elseif (datetime_le_refdatetime( date_time, dt_end_f521_fp )) then
          
          ! Note "less-than-or-equal" (_le_) above
          
          stream = 'f521_fp'            ! use GEOS-5.21.x output

       elseif (datetime_le_refdatetime( date_time, dt_end_f522_fp )) then
          
          ! Note "less-than-or-equal" (_le_) above
          
          stream = 'f522_fp'            ! use GEOS-5.22.x output

       elseif (datetime_le_refdatetime( date_time, dt_end_f525_fp )) then
          
          ! Note "less-than-or-equal" (_le_) above
          
          stream = 'f525_fp'            ! use GEOS-5.25.x output

       else
          
          stream = 'f525_p5_fp'         ! use GEOS-5.25_p5.x output

       end if

    else
       
       stream = met_tag_out
       
    end if
    
    met_tag_out = stream
    

    ! interpret optional tag segments
    
    do ii=1,N_opt_tag
       
       if     (tmp_tag(ii)(1:4)=='prec') then
          
          use_prec_corr = .true.
          
          prec_tag      = tmp_tag(ii)(5:len(tmp_tag(ii)))
          
       elseif (tmp_tag(ii)(1:4)=='Nx+-') then
          
          use_Predictor = .true.

       else
          
          err_msg = 'invalid met_tag_in, unknown optional tag segment'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end if
       
    end do
    
    ! get path to met forcing (except precip)
    
    met_path_default = trim(met_path_in) // '/' 
    
    ! get path to precip
    
    if (use_prec_corr) then      ! using precip corrections
       
       met_path_prec = trim(met_path_in) // '/precip_corr_' // trim(prec_tag) // '/'
       
    else                         ! *not* using precip corrections
       
       met_path_prec = met_path_default
       
    end if
    
    ! Double-check if optional tag segments were somehow missed, e.g., 
    !  because they were accidentally appended with single underscores.
    ! Assumes that "prec" and "Nx+-" never appear in GEOS-5 product names.
    ! Does NOT protect against spelling errors, e.g., "__perc" or "__Nx-+".
    ! - reichle, 24 Nov 2015
    
    if ((.not. use_prec_corr) .and. (index( met_tag_in, 'prec')>0)) then
       
       err_msg = 'questionable met_tag_in: includes "prec" but use_prec_corr=.false.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    
    end if

    if ((.not. use_Predictor) .and. (index( met_tag_in, 'Nx+-')>0)) then
       
       err_msg = 'questionable met_tag_in: includes "Nx+-" but use_Predictor=.false.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    
    end if
               
  end subroutine parse_G5DAS_met_tag
  
  ! ****************************************************************
 
  subroutine get_GEOS_forcing_filename(fname_full,file_exists, date_time, daily_file, met_path, met_tag, &
       GEOSgcm_defs, file_ext)

    ! reichle, 27 Jul 2015 - added "daily_file" option 
    !  (because MERRA-2 provides aggregated daily files that contain 24 hourly fields)

    ! reichle,  7 May 2020 - added "seamless" GEOS FP file names (so far, only "asm" assimilation files)

    implicit none
    
    character(*),                     intent(inout) :: fname_full
    logical,                          intent(out)   :: file_exists
    type(date_time_type),             intent(in)    :: date_time
    logical,                          intent(in)    :: daily_file
    character(*),                     intent(in)    :: met_path
    character(*),                     intent(in)    :: met_tag
    character( 40),     dimension(5), intent(in)    :: GEOSgcm_defs
    character(*),                     intent(in)    :: file_ext
    
    ! local variables
    
    character(300) :: fname, fname_full_tmp1, fname_full_tmp2
    character( 14) :: time_stamp
    character(  4) :: YYYY,  HHMM, day_dir
    character(  2) :: MM,    DD  

    integer        :: tmpind, tmpindend

    character(len=*), parameter :: Iam = 'get_GEOS_forcing_filename'
    
    ! assemble date/time strings
    
    write (YYYY,'(i4.4)') date_time%year
    write (MM,  '(i2.2)') date_time%month
    write (DD,  '(i2.2)') date_time%day
    write (HHMM,'(i4.4)') date_time%hour*100+date_time%min
    
    ! deal with absence of "minutes" in "bkg.sfc" file name
    ! (replace minutes with blanks, then trim, see below)
    
    if (trim(GEOSgcm_defs(3))=='bkg.sfc')  HHMM = HHMM(1:2) // '  '  

    ! assemble file name

    time_stamp = repeat(' ', len(time_stamp))

    if (daily_file) then
       
       time_stamp(1:8) = YYYY // MM // DD
       
    else
       
       time_stamp      = YYYY // MM // DD // '_' // trim(HHMM) // 'z'
       
    end if
    
    if (trim(met_tag(1:11))=='GEOS.fp.asm') then

       ! GEOS FP with generic file names, e.g., 
       !
       !    GEOS.fp.asm.inst1_2d_lfo_Nx.20200507_0000.V01.nc4
       !
       ! for now, always use product counter V01 
       !  (as of 7 May 2020, no V02 or higher was issued for GEOS FP "lfo" products 
       !   going back to Jun 2013)

       fname = trim(met_tag) // '.' // trim(GEOSgcm_defs(3)) // '.' // &
            trim(time_stamp(1:13)) // '.V01.' // trim(file_ext)
       
       ! archived files are stored in daily directories 

       day_dir = 'D' // DD // '/'

    else
       
       ! GEOS FP with experiment-specific file names and MERRA-2, e.g.,
       !
       !    f525_p5_fp.inst1_2d_lfo_Nx.20200507_0000z.nc4
       !    MERRA2_400.inst1_2d_lfo_Nx.20200507.nc4
       
       fname = trim(met_tag) // '.' // trim(GEOSgcm_defs(3)) // '.' // &
            trim(time_stamp) // '.' // trim(file_ext)

       ! archived files are stored in monthly directories
       !   (per LDAS legacy convention for GEOS FP with experiment-specific file names)
       
       day_dir = repeat(' ', len(day_dir))
       
    end if
    
    ! ----------------------------------------------
    ! 
    ! find suitable file in a couple of places

    file_exists = .false.                              ! initialize
    
          
    ! first try: year/month[/day] directory
          
    fname_full = trim(met_path) // '/' // trim(met_tag) // '/' //         &
         trim(GEOSgcm_defs(4)) // '/Y' // YYYY  // '/M' // MM // '/' //   &
         trim(day_dir) // trim(fname)
    
    inquire(file=fname_full, exist=file_exists)
    
    if (file_exists) return                            ! done

    fname_full_tmp1 = trim(fname_full)                 ! remember for error log below
          
        
    ! second try: ./met_path (coupled land-atm DAS)          
    
    fname_full = trim(met_path) // '/' // trim(fname)
    
    inquire(file=fname_full, exist=file_exists)
    
    if (file_exists) return                            ! done
      
    fname_full_tmp2 = trim(fname_full)                 ! remember for error log below


    ! last try: for GEOS FP with generic file names, try product counter '.V02.' in year/month/day dir
 
    if (trim(met_tag(1:11))=='GEOS.fp.asm') then

       fname_full = fname_full_tmp1                    ! from first try

       ! GEOS.fp.asm.inst1_2d_lfo_Nx.20200507_0000.V01.nc4
       !                                           1234567
       
       tmpindend = len_trim(fname_full)
       tmpind    = len_trim(file_ext)
       
       tmpind    = tmpindend - tmpind - 3
       
       fname_full( tmpind:tmpind+2 ) = 'V02'           ! --> *.V02.nc4
              
       inquire(file=fname_full, exist=file_exists)
    
    end if
   

    ! if no file was found, report file names that were tried
    
    if (.not. file_exists) then
       if(root_logit) then
          print '(400A)',  trim(Iam) // ': Could not find any of the following files:'
          print '(400A)',  trim(fname_full_tmp1)
          print '(400A)',  trim(fname_full_tmp2)
          print '(400A)',  trim(fname_full)
       endif
    endif    

  end subroutine get_GEOS_forcing_filename

  !**********************************************************

  subroutine GEOS_openfile(FileOpenedHash, fname_full, fid, tile_coord, m_hinterp, rc, lat_str, lon_str)

    ! open GEOS forcing file, extract coord info, prep horizontal interpolation info (if not done already)
    !
    ! ASSUMPTIONS (as of 23 Jun 2021):
    ! - forcing grid is global
    ! - if lat/lon forcing grid, pole is on center of grid cell and dateline is on center of grid cell
    ! - if cube-sphere forcing grid, tile space (tile_coord) is associated with same cube-sphere grid
    
    use netcdf
    implicit none
    include 'mpif.h'
    type(Hash_Table),                    intent(inout) :: FileOpenedHash
    character(*),                        intent(in)    :: fname_full
    integer,                             intent(out)   :: fid
    type(tile_coord_type), dimension(:), intent(in)    :: tile_coord
    integer,                             intent(in)    :: m_hinterp
    integer,               optional,     intent(out)   :: rc
    character(*),          optional,     intent(in)    :: lat_str, lon_str
    
    integer                 :: N_lat, N_lon, N_cat, N_f
    integer, allocatable    :: i1(:), i2(:), j1(:), j2(:)
    real,    allocatable    :: x1(:), x2(:), y1(:), y2(:)
    integer                 :: ierr, k
    integer                 :: latid, lonid, nfid, xdimid
    real                    :: dlon, dlat, ll_lon, ll_lat, this_lon, this_lat, tmp_lon, tmp_lat
    integer                 :: icur, jcur, inew, jnew
    real                    :: xcur, ycur, xnew, ynew
    character(len=100)      :: err_msg
    character(*), parameter :: Iam="GEOS_openfile"
    logical                 :: isCubed
    ! add mpi
    type(ESMF_VM)           :: vm
    integer                 :: comm, total_prcs, myrank
    integer                 :: status

    ! initialize hash table, if not already initialized
    
    call FileOpenedHash%init()
    
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    call ESMF_VmGet(vm, mpicommunicator=comm, rc=status)
    VERIFY_(status)

    ! get file ID from hash table
    
    call FileOpenedHash%get(fname_full,fid)
    
    if( fid == -9999 ) then ! not open yet
       ierr=nf90_open(fname_full,NF90_NOWRITE, fid) ! open file  
       if(root_logit) then
          write(logunit,'(400A)') "opening file: "//trim(fname_full)
       endif
       _ASSERT( ierr == nf90_noerr, "nf90 error")
       call FileOpenedHash%put(fname_full,fid)      ! record file ID in hash table
    endif

    ! --------------------------------------------------
    !
    ! determine grid dimension, spacing, and offset/extent
    !
    ! ONLY the grid dimensions are read from the file;
    ! the grid spacing and offset/extent are determined based on hard-wired assumptions!
    
    ! check if forcing file is on cs grid
    ierr =  nf90_inq_dimid(fid,"nf",nfid)   ! check if number-of-faces (nf) dimension exists 
    
    if (ierr == nf90_noerr) then ! it is cubed-sphere grid if face dimension is found
       
       ! cube-sphere grid: need N_f (number-of-faces), N_lon, N_lat
       
       ierr  =  nf90_inq_dimid(fid,"Xdim",xdimid)
       _ASSERT( ierr == nf90_noerr, "nf90 error")
       ierr  =  nf90_Inquire_Dimension(fid,nfid,  len=N_f)
       _ASSERT( ierr == nf90_noerr, "nf90 error")
       _ASSERT( N_f == 6, "number of (cubed-sphere) faces not equal to 6")
       ierr  =  nf90_Inquire_Dimension(fid,xdimid,len=N_lon)
       _ASSERT( ierr == nf90_noerr, "nf90 error")
       _ASSERT( N_lon == im_world_cs, "forcing on cube-sphere grid: forcing grid dimension must match native grid dimension (grid associated with tile space)")
       N_lat = N_f*N_lon
       _ASSERT( m_hinterp == 0, "forcing on cubed-sphere grid requires nearest-neighbor interpolation (m_hinterp = 0)")
       isCubed = .true.       

    else

       ! lat/lon grid: need N_lon, N_lat, dlon, dlat, ll_lon, ll_lat
       
       if (present(lat_str) .and. present(lon_str)) then
          ierr =  nf90_inq_dimid(fid,trim(lat_str),latid)
          ierr =  nf90_inq_dimid(fid,trim(lon_str),lonid)
       else
          ierr =  nf90_inq_dimid(fid,"lat",latid)
          ierr =  nf90_inq_dimid(fid,"lon",lonid)
       endif
       ierr =  nf90_Inquire_Dimension(fid,latid,len=N_lat)
       ierr =  nf90_Inquire_Dimension(fid,lonid,len=N_lon)
       isCubed = .false.       
       
       ! assume global grid w/ dateline on center and pole on center  
       
       dlon   =  360./real(N_lon)
       dlat   =  180./real(N_lat-1)
       ll_lon = -180. - dlon/2.
       ll_lat =  -90. - dlat/2.

    endif

    ! -----------------------------------------------------------------

    ! prep interpolation from forcing grid to model tiles
    
    N_cat = size(tile_coord,1)
    
    ! if dimensions are the same, no need to recalculate the local_info
    if( local_info%N_lat == N_lat .and. local_info%N_lon == N_lon ) then
       RETURN_(ESMF_SUCCESS)
    endif
    
    allocate(i1(N_cat),j1(N_cat))
    allocate(i2(N_cat),j2(N_cat),x1(N_cat),x2(N_cat),y1(N_cat),y2(N_cat))
    
    select case (m_hinterp)
       
    case (0)  ! nearest-neighbor
       
       ! compute indices for nearest neighbor interpolation from GEOSgcm grid to tile space
       
       if( isCubed ) then ! cs grid
          
          ! cube-sphere grid of forcing data must match cube-sphere grid associated with tile space
          do k=1,N_cat
             i1(k) = tile_coord(k)%i_indg
             j1(k) = tile_coord(k)%j_indg
          enddo

       else

          do k=1,N_cat
             
             ! ll_lon and ll_lat refer to lower left corner of grid cell
             ! (as opposed to the grid point in the center of the grid cell)
             
             this_lon = tile_coord(k)%com_lon
             this_lat = tile_coord(k)%com_lat
             
             i1(k) = ceiling((this_lon - ll_lon)/dlon)
             j1(k) = ceiling((this_lat - ll_lat)/dlat)
             
             ! NOTE: For a "date line on center" grid and (180-dlon/2) < lon < 180
             !  we now have i1=(grid%N_lon+1)
             ! This needs to be fixed as follows:
             
             if (i1(k)> N_lon)  i1(k)=1
             
          end do
          
       endif ! cs grid

    case (1)  ! bilinear interpolation
       
       ! compute indices of nearest neighbors needed for bilinear
       ! interpolation from GEOSgcm grid to tile space
       
       do k=1,N_cat
          
          ! ll_lon and ll_lat refer to lower left corner of grid cell
          ! (as opposed to the grid point in the center of the grid cell)
          
          ! pchakrab: For bilinear interpolation, for each tile, we need:
          !  x1, x2, y1, y2 (defining the co-ords of four neighbors) and
          !  i1, i2, j1, j2 (defining the indices of four neighbors)
          
          ! find nearest neighbor forcing grid cell ("1")

          ! com of kth tile
          this_lon = tile_coord(k)%com_lon
          this_lat = tile_coord(k)%com_lat
          icur =  ceiling((this_lon - ll_lon)/dlon)
          jcur =  ceiling((this_lat - ll_lat)/dlat)
          
          ! wrap-around
          if (icur>N_lon) icur = 1
          if (jcur>N_lat) then
             err_msg = "encountered tile near the poles"
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if

          ! determine lon and lat of forcing grid cell center
          xcur = real(icur-1)*dlon - 180.0    ! ASSUMES dateline on center!
          ycur = real(jcur-1)*dlat -  90.0    ! ASSUMES pole on center!

          ! finalize "1"
          i1(k) = icur
          j1(k) = jcur
          x1(k) = xcur    ! lon of grid cell center
          y1(k) = ycur    ! lat of grid cell center
          
          ! find forcing grid cell ("2") diagonally across from icur, jcur

          ! to this end, determine quadrant of forcing grid cell that contains com of tile (this_lon, this_lat)
          
          ! quadrant is found by determining nearest-neighbor indices (inew, jnew) when tile is shifted 
          ! by 1/2 of the grid-spacing to the northeast (upper right)
          tmp_lon = this_lon + 0.5*dlon               
          tmp_lat = this_lat + 0.5*dlat
          inew =  ceiling((tmp_lon  - ll_lon)/dlon)   ! nearest-neighbor i index of shifted location
          jnew =  ceiling((tmp_lat  - ll_lat)/dlat)   ! nearest-neighbor j index of shifted location
          if (inew==icur) inew = inew - 1             ! if shift results in same lon index, go west (left)
          if (jnew==jcur) jnew = jnew - 1             ! if shift results in same lat index, go south (down)
          
          ! determine center lon and lat of forcing grid cell ("2");
          ! must do this BEFORE wrap-around (such that xnew=x2 will be outside of [-180:180] near dateline),
          ! otherwise distance calculation would not work near dateline
          xnew = real(inew-1)*dlon - 180.0            ! ASSUMES dateline on center!
          ynew = real(jnew-1)*dlat -  90.0            ! ASSUMES pole on center!
          
          ! wrap-around
          if (inew==0) inew = N_lon
          if (inew>N_lon) inew = 1
          if ((jnew==0) .or. (jnew>N_lat)) then
             err_msg = "encountered tile near the poles"
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          i2(k) = inew
          j2(k) = jnew
          x2(k) = xnew    ! lon of grid cell center (note that this intentionally 
          y2(k) = ynew    ! lat of grid cell center
          
       end do
       
    case default
       
       err_msg = "unknown horizontal interpolation method"
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select

    local_info%N_lat = N_lat
    local_info%N_lon = N_lon
    local_info%N_cat = N_cat
    call move_alloc(i1,local_info%i1)
    call move_alloc(i2,local_info%i2)
    call move_alloc(j1,local_info%j1)
    call move_alloc(j2,local_info%j2)
    call move_alloc(x1,local_info%x1)
    call move_alloc(x2,local_info%x2)
    call move_alloc(y1,local_info%y1)
    call move_alloc(y2,local_info%y2)
    
    RETURN_(ESMF_SUCCESS)
    
  end subroutine GEOS_openfile

  ! ****************************************************************
  
  subroutine GEOS_closefile(fid)
     use netcdf
     implicit none
     integer,intent (in) :: fid
     integer :: ierr

     ierr = nf90_close(fid)
     if(ierr /= nf90_noerr) then
        print *, " error GEOS_closefile"
        stop 2
     endif
     
  endsubroutine 
  
  ! ****************************************************************

  subroutine get_GEOS_corr_prec_filename(fname_full,file_exists, date_time, met_path, met_tag, &
       GEOSgcm_defs, file_ext )
    
    implicit none
    character(*),                 intent(inout) :: fname_full
    logical,intent(out)                         :: file_exists 
    type(date_time_type),         intent(in)    :: date_time
    character(*),                 intent(in)    :: met_path
    character(*),                 intent(in)    :: met_tag
    character( 40), dimension(5), intent(in)    :: GEOSgcm_defs
    character(*),                 intent(in)    :: file_ext
        
    ! local variables

    character(100) :: fname
    character(200) :: fdir
    character(300) :: fname_full_tmp1, fname_full_tmp2
    character(  4) :: YYYY,  HHMM
    character(  2) :: MM,    DD

    integer        :: tmpind, tmpindend

    character(len=*), parameter :: Iam = 'get_GEOS_corr_prec_filename'

    ! assemble date/time strings
    
    write (YYYY,'(i4.4)') date_time%year
    write (MM,  '(i2.2)') date_time%month
    write (DD,  '(i2.2)') date_time%day
    write (HHMM,'(i4.4)') date_time%hour*100+date_time%min
    
    ! assemble file name

    if (trim(met_tag(1:11))=='GEOS.fp.asm') then
       
       ! for now, always use product counter V01 
       !  (as of 7 May 2020, no V02 or higher was issued for GEOS FP "lfo" products 
       !   going back to Jun 2013)

       fname = trim(met_tag) // '.' // trim(GEOSgcm_defs(3)) // '_corr.' //         &
            YYYY // MM // DD // '_' // trim(HHMM) // '.V01.' // trim(file_ext)
       
    else
       
       fname = trim(met_tag) // '.' // trim(GEOSgcm_defs(3)) // '_corr.' //         &
            YYYY // MM // DD // '_' // trim(HHMM) // 'z.'    // trim(file_ext)
       
    end if

    ! assemble dir name with "/Yyy" (year) dir but without "/Mmm" (month) dir

    fdir = trim(met_path) // '/' // trim(met_tag)   // '/' //                       &
         trim(GEOSgcm_defs(4)) // '/' // 'Y' // YYYY // '/'
    
    ! -----------------------------------------------------------------------
    
    file_exists = .false.                              ! initialize
    

    ! first try:  look for file in year/month dir
    !  (LDAS standard for corrected G5DAS precip)
    
    fname_full = trim(fdir) // 'M' // MM // '/' // trim(fname)
    
    inquire(file=fname_full, exist=file_exists)
    
    if (file_exists) return                            ! done
    
    fname_full_tmp1 = trim(fname_full)                 ! remember for error log below
    

    ! second try: *without* "/Mmm" (month) dir 

    ! THIS TRY IS PROBABLY OBSOLETE BUT COULD EASILY BE TWEAKED TO LOOK 
    ! IN year/month/day DIRECTORY (WHICH POSSIBLY APPLIES TO CORRECTED
    ! PRECIP FROM GMAO OPS THAT IS INPUT INTO MERRA-2)
    ! - reichle 10 May 2020

    fname_full = trim(fdir) // trim(fname)

    inquire(file=fname_full, exist=file_exists)

    if (file_exists) return                            ! done
    
    fname_full_tmp2 = trim(fname_full)                 ! remember for error log below


    ! last try: for GEOS FP with generic file names, try product counter '.V02.' in year/month dir
 
    if (trim(met_tag(1:11))=='GEOS.fp.asm') then

       fname_full = fname_full_tmp1                    ! from first try

       ! GEOS.fp.asm.inst1_2d_lfo_Nx.20200507_0000.V01.nc4
       !                                           1234567
       
       tmpindend = len_trim(fname_full)
       tmpind    = len_trim(file_ext)
       
       tmpind    = tmpindend - tmpind - 3

       fname_full( tmpind:tmpind+2 ) = 'V02'           ! --> *.V02.nc4
       
       inquire(file=fname_full, exist=file_exists)
    
    end if
    
    
    ! if no file was found, report file names that were tried
    
    if( .not. file_exists ) then
       if(root_logit) then 
          print '(400A)',  trim(Iam) // ': Could not find any of the following files:'
          print '(400A)',  trim(fname_full_tmp1)
          print '(400A)',  trim(fname_full_tmp2)
          print '(400A)',  trim(fname_full)
       endif
    endif
    
  end subroutine get_GEOS_corr_prec_filename
 
  ! ****************************************************************
  
  subroutine get_GSWP2_1x1_netcdf(date_time, met_path, N_catd, tile_coord, &
       met_force_new, nodata_forcing )
    
    ! read GSWP2_NetCDF files and extract forcings in tile space
    ! (uses nearest neighbor interpolation)
    
    ! reichle, 28 Jul 2005
    
    implicit none
    
    type(date_time_type), intent(in) :: date_time
    
    character(*),       intent(in) :: met_path
    
    integer,              intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(met_force_type), dimension(N_catd), intent(inout) :: met_force_new
    
    real, intent(out) :: nodata_forcing
    
    ! GEOS5-DAS grid and netcdf parameters 
    
    integer, parameter :: gswp2_grid_N_lon  = 360
    integer, parameter :: gswp2_grid_N_lat  = 150
    real,    parameter :: gswp2_grid_ll_lon = -180.
    real,    parameter :: gswp2_grid_ll_lat = -60.
    real,    parameter :: gswp2_grid_dlon   = 1.
    real,    parameter :: gswp2_grid_dlat   = 1.
    
    integer, parameter :: N_gswp2_compressed = 15238
    
    ! GSWP2 forcing time step in hours

    integer, parameter :: dt_gswp2_in_hours = 3
    
    integer, parameter :: nciv_land = 3
    integer, parameter :: nciv_data = 6
    
    integer, parameter :: N_gswp2_vars = 9
    
    real,    parameter :: nodata_gswp2 = 1.e20
    
    character(40), dimension(N_gswp2_vars), parameter :: gswp2_name = (/ &
         'SWdown_srb     ',   &   !  1 - flux
         'LWdown_srb     ',   &   !  2 - flux
         'Rainf_C_gswp   ',   &   !  3 - flux
         'Rainf_gswp     ',   &   !  4 - flux
         'Snowf_gswp     ',   &   !  5 - flux
         'PSurf_ecor     ',   &   !  6 - state
         'Qair_cru       ',   &   !  7 - state
         'Tair_cru       ',   &   !  8 - state
         'Wind_ncep      ' /)     !  9 - state
    
    ! local variables
    
    real :: tol 
    
    real, dimension(gswp2_grid_N_lon,gswp2_grid_N_lat) :: tmp_grid
    
    integer, dimension(N_gswp2_compressed)   :: land_gswp2
    integer, dimension(N_gswp2_compressed)   :: land_i_gswp2, land_j_gswp2
    integer, dimension(N_catd)               :: i_ind, j_ind
    
    real,    dimension(N_gswp2_compressed)   :: tmp_vec
    
    real,    dimension(N_catd,N_gswp2_vars) :: force_array
    
    integer, dimension(2) :: start, icount
    
    integer :: k, n, hours_in_month, gswp2_var, ierr, ncid

    integer :: int_dt_gswp2_in_seconds
    
    real    :: this_lon, this_lat
    
    type(date_time_type) :: date_time_tmp
    
    character(4) :: YYYY
    character(2) :: MM
    
    character(300) :: fname

    character(len=*), parameter :: Iam = 'get_GSWP2_1x1_netcdf'
    character(len=400) :: err_msg
    
    ! --------------------------------------------------------------------
    
    int_dt_gswp2_in_seconds = 3600*dt_gswp2_in_hours
    
    nodata_forcing = nodata_gswp2
    
    tol = abs(nodata_forcing*nodata_tolfrac_generic)    

    ! ----------------------------------------------    
    !
    ! compute indices for nearest neighbor interpolation from GSWP2 grid 
    ! to tile space
    !
    ! (NOTE: this should at some point be replaced with a regridding
    !  subroutine that interpolates from the
    !  native forcing grid to the GCM atmospheric grid that is used
    !  to cut catchments into tiles - then "standard" grid2tile
    !  using tile_coord%atm_i and tile_coord%atm_j applies. 
    !  reichle, 26 May 2005)
    
    do k=1,N_catd
       
       ! ll_lon and ll_lat refer to lower left corner of grid cell
       ! (as opposed to the grid point in the center of the grid cell)
       
       this_lon = tile_coord(k)%com_lon
       this_lat = tile_coord(k)%com_lat

       ! i_ind, j_ind count eastward and northward
       ! (note that lat/lon coordinates in GSWP2 netcdf files
       !  are eastward and southward)
       
       i_ind(k) = ceiling( (this_lon - gswp2_grid_ll_lon)/gswp2_grid_dlon )
       j_ind(k) = ceiling( (this_lat - gswp2_grid_ll_lat)/gswp2_grid_dlat )
              
    end do
    
    ! ------------------------------------------------------
    !
    ! space dimension is same for all variables
    
    start(1)  = 1
    icount(1) = N_gswp2_compressed
   
    ! check for possible error with time
    
    if ( (date_time%min/=0) .or. (date_time%sec/=0) .or.           &
         (mod(date_time%hour,dt_gswp2_in_hours)/=0)        ) then
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'timing ERROR!!')
       
    end if
       
    ! ------------------------------------------------------
    !
    ! get forcing data
    
    do gswp2_var = 1,N_gswp2_vars
       
       ! time dimension 
       !
       ! First entry in GSWP2_NetCDF file is at 3Z, with fluxes for 0Z-3Z
       !
       ! At 0Z for first day of month:
       !  - for fluxes read first entry of that month
       !  - for states read last entry of preceding month
       ! At 3Z for first day of month:
       !  - for fluxes read second entry of that month
       !  - for states read first entry of that month
       ! and so on...
       
       select case (gswp2_var)
          
       case (1,2,3,4,5)   ! "fluxes"
          
          date_time_tmp = date_time
          
       case (6,7,8,9)     ! "states"
          
          date_time_tmp = date_time
          
          call augment_date_time( -int_dt_gswp2_in_seconds, date_time_tmp )
          
       case default
          
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error')
          
       end select
       
       hours_in_month = (date_time_tmp%day-1)*24 + date_time_tmp%hour
       
       start(2)  = hours_in_month / dt_gswp2_in_hours + 1 
       icount(2) = 1
       
       ! assemble year and month strings
       
       write (YYYY,'(i4.4)') date_time_tmp%year
       write (MM,  '(i2.2)') date_time_tmp%month
              
       ! assemble file name, open file

       fname = trim(met_path) // trim(gswp2_name(gswp2_var)) // '/'     &
            // '/' // trim(gswp2_name(gswp2_var)) // YYYY // MM // '.nc'
       
       if(root_logit) write (logunit,*) 'opening ' // trim(fname)
       
       ierr = NF_OPEN(fname,NF_NOWRITE,ncid)
       
       if (ierr/=0) then
          err_msg = 'error opening netcdf file'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
                     
       ! ------------------------------------------------------
       !
       ! read compression parameters (same for all data variables and time steps)
       
       if (gswp2_var == 1) then
          
          if(root_logit) write (logunit,*) 'get netcdf compression params from ' // trim(fname)
          
          if (ierr/=0) then
             err_msg = 'error opening netcdf file'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if
          
          ierr = NF_GET_VARA_INT(ncid, nciv_land, start, icount, land_gswp2)

          ! land_i_gswp2, land_j_gswp2 count eastward and southward
          ! (note that lat/lon coordinates in GSWP2 netcdf files
          !  are eastward and southward)
          
          do k=1,N_gswp2_compressed
             
             land_j_gswp2(k) = (land_gswp2(k)-1)/gswp2_grid_N_lon  + 1
             land_i_gswp2(k) = &
                  land_gswp2(k) - (land_j_gswp2(k)-1)*gswp2_grid_N_lon 
             
          end do
          
       end if
       
       ! ------------------------------------------------------
       !
       ! read compressed data, and put on global grid
       
       ierr = NF_GET_VARA_REAL(ncid, nciv_data, start, icount, tmp_vec )
       
       ierr = NF_CLOSE(ncid)
       
       tmp_grid = nodata_forcing

       do n=1,N_gswp2_compressed
          
          tmp_grid(land_i_gswp2(n), land_j_gswp2(n) ) = tmp_vec(n)
          
       end do
       
       ! flip tmp_grid 
       ! (land_j_gswp2 counts southward, whereas j_ind counts northward)
       
       tmp_grid = tmp_grid(:,gswp2_grid_N_lat:1:-1)
       
       !!do k=1,gswp2_grid_N_lon
       !!   do n=1,gswp2_grid_N_lat
       !!      
       !!      write (999,*) k, n, tmp_grid(k,n)
       !!
       !!   end do
       !!end do
       !!write (logunit,*) 'debug stop here'
       !!stop
       
       ! interpolate to tile space
       
       ! (NOTE: This should at some point be replaced with a regridding
       !  subroutine that interpolates from the
       !  native forcing grid to the GCM atmospheric grid that is used
       !  to cut catchments into tiles - then "standard" grid2tile
       !  using tile_coord%atm_i and tile_coord%atm_j applies. 
       !  reichle, 26 May 2005)
       
       do k=1,N_catd
          
          force_array(k,gswp2_var) = tmp_grid(i_ind(k), j_ind(k))
          
       end do
       
    end do
    
    ! convert variables and units of force_array to match met_force_type, 
    ! put into structure
    
    ! from GSWP2 files:
    !
    !  force_array(:, 1) =  SWdown_srb    W/m2 
    !  force_array(:, 2) =  LWdown_srb    W/m2 
    !  force_array(:, 3) =  Rainf_C_gswp  kg/m2/s
    !  force_array(:, 4) =  Rainf_gswp    kg/m2/s
    !  force_array(:, 5) =  Snowf_gswp    kg/m2/s    
    !  force_array(:, 6) =  PSurf_ecor    Pa
    !  force_array(:, 7) =  Qair_cru      kg/kg  
    !  force_array(:, 8) =  Tair_cru      K  
    !  force_array(:, 9) =  Wind_ncep     m/s   

    met_force_new%SWdown  = force_array(:,1)
    met_force_new%LWdown  = force_array(:,2)
    met_force_new%Rainf_C = force_array(:,3)
    met_force_new%Rainf   = force_array(:,4)
    met_force_new%Snowf   = force_array(:,5)
    met_force_new%Psurf   = force_array(:,6)
    met_force_new%Qair    = force_array(:,7)
    met_force_new%Tair    = force_array(:,8)
    met_force_new%Wind    = force_array(:,9)
    
  end subroutine get_GSWP2_1x1_netcdf
  
  ! ****************************************************************
   
  subroutine check_forcing_nodata( N_catd, tile_coord, nodata_forcing, met_force )

    ! Check for nodata values in met_force and fill with "neighboring" data if sensible;
    !   otherwise, check_forcing_nodata() will abort.
    ! Upon successful exit, there will *not* be nodata-values in met_force except
    !   possibly RefH, which is not checked.
    !
    ! (Note: subroutine repair_forcing() checks for unphysical values.)
    !
    ! Owing to differences in land masks, some land-only forcing datasets may have
    !  only nodata-values for some GEOS land tiles.  There may also be intermittent
    !  nodata-values in the forcing dataset.    
    ! If a nodata-value is encountered, use the value from the "next" tile, where
    !  "next" is next in tile order, provided the "next" tile is within "max_distance".
    !  Abort if this does not work.
    ! The "next" tile approach is used to avoid the costly determination of the nearest
    !  tile, which would require communications across processors and long loops.
    ! For details, see helper subroutine check_forcing_nodata_2(), which also offers 
    !  a compile-time switch to create a list of all tiles that lack forcing data. 
    !
    ! reichle, 13 May 2003
    ! reichle, 13 Jun 2005
    ! reichle, 12 Feb 2021 - added "max_distance" limit, revised comments
    ! reichle, 22 Apr 2021 - PAR must now also have "good" data at this stage, added to checks
    
    implicit none
    
    integer, intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    real,    intent(in) :: nodata_forcing
    
    type(met_force_type), dimension(N_catd), intent(inout) :: met_force
    
    ! ------------------------------------------------------------
    
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%Tair   )
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%Qair   )
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%Psurf  )
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%Rainf_C)
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%Rainf  )
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%Snowf  )
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%LWdown )
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%SWdown )
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%PARdrct)
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%PARdffs)
    call check_forcing_nodata_2(N_catd,tile_coord,nodata_forcing,met_force%Wind   )
    
    ! do NOT call check_forcing_nodata_2() for "RefH" (should not have any problems)
    
  end subroutine check_forcing_nodata
  
  ! *****************************************************************
  
  subroutine check_forcing_nodata_2( N_catd, tile_coord, nodata_forcing, force_vec )
    
    ! helper subroutine for check_forcing_nodata()
    !
    ! reichle, 13 Jun 2005
    ! reichle, 12 Feb 2021 - added "max_distance" limit, revised comments
    
    implicit none
    
    integer, intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    real,    intent(in) :: nodata_forcing

    real, dimension(:), intent(inout) :: force_vec
    
    ! local variables    

    real, parameter :: max_distance = 0.5     ! [degrees]

    real :: tol, distsq_next
    
    integer :: i, i_next, N_exclude
    
    ! Set the following logical to .true. to generate an "ExcludeList"
    ! for a given forcing data set.  The list contains all tiles for
    ! which the specified forcing dataset has only nodata-values
    ! within "max_distance".  The list will end up in the file
    ! "fort.9999".  Next, run "ldas_setup" again and exclude the
    ! tiles in this list from the simulation domain (see "ExcludeList"
    ! in "exeinp" setup file.)
    
    logical, parameter :: create_ExcludeList = .false.
    
    
    character(len=*), parameter :: Iam = 'check_forcing_nodata_2'
    character(len=400) :: err_msg
    
    ! ------------------------------------------------------------
    
    N_exclude = 0
    
    tol = abs(nodata_forcing*nodata_tolfrac_generic)    

    do i=1,N_catd
       
       ! nodata-value checks 
       
       i_next = min(i+1,N_catd) 
       
       if (abs(force_vec(i)-nodata_forcing)<tol) then   ! forcing is nodata at tile i
          
          if (create_ExcludeList) then
             
             N_exclude = N_exclude + 1
             write (9999,*) tile_coord(i)%tile_id
             
          else

             ! compute square of distance to "next" tile
             
             distsq_next =                                                    &
                  (tile_coord(i)%com_lon - tile_coord(i_next)%com_lon)**2 +   &
                  (tile_coord(i)%com_lat - tile_coord(i_next)%com_lat)**2 
             
             if ( (abs(force_vec(i_next)-nodata_forcing)>tol) .and.                &
                  (distsq_next .le. max_distance**2)                 )    then

                ! "next" tile has good forcing data and is within "max_distance"
                !  --> use forcing from "next" tile for tile i, add note in log file.
                
                if (root_logit)  write (logunit,*) 'forcing has nodata-value in tile ID = ', &
                     tile_coord(i)%tile_id, '. Using forcing from nearby tile.'
                force_vec(i)=force_vec(i_next)
                
             else

                ! cannot find forcing data for tile i, abort with message
                
                write (tmpstring100,*) tile_coord(i)%tile_id
                err_msg = 'forcing has nodata-value in tile ID = ' // trim(tmpstring100) //     &
                     '. No good forcing nearby. ' //                                            &
                     'Use compile-time switch "create_ExcludeList" to create ' //               &
                     'a complete list for use in "ldas_setup".'  
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             end if
             
          end if
          
       end if
       
    end do
    
    if (create_ExcludeList) then
       if(root_logit)  write (logunit,*) '---------------------------------------------------------------'
       if(root_logit)  write (logunit,*) ' found N_exclude = ',N_exclude, ' tiles that should be in ExcludeList'
       err_msg = 'ExcludeList now in file fort.9999.  Use this info in ldas_setup.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
  end subroutine check_forcing_nodata_2
  
  ! ******************************************************************

  type(date_time_type) function shift_forcing_date( met_tag, date_time )

    ! shift date_time by years or days, useful for twin experiments
    !
    ! examples: 
    !
    !  met_tag = "RedArk_ASCII_shift-1year"  uses 1980 forcing for 1981 and so on
    !  met_tag = "RedArk_ASCII_shift+3day"   uses Jan 4 forcing for Jan 1 and so on
    !
    ! reichle, 6 Apr 2007
    
    implicit none
    
    character(80)  :: met_tag
    
    type(date_time_type) :: date_time, date_time_tmp

    integer :: is, ie, shift

    character(300) :: tmpstring300

    ! initialize
    
    date_time_tmp = date_time
    
    ! check whether met_tag asks for shift
    
    is = index( met_tag, 'shift' )
    
    if (is>0) then
       
       ! make sure this is only used for pre-specified forcing data sets
       ! for now permit use with all "RedArk" data sets (reichle, 6 Apr 2007)
       
       if (index(met_tag,'RedArk')==0) then
          
          tmpstring300 = 'shift_forcing_date(): Are you sure? ' // &
               'If so, edit source code and recompile.'                    
          
          if(root_logit) write (logunit,*) tmpstring300
          write(0,*) tmpstring300
          stop

       end if

       ie = index( met_tag, 'year')
       
       if (ie>0) then
          
          read (met_tag(is+5:ie-1),'(i1)') shift

          date_time_tmp%year = date_time%year + shift
          
          ! deal with leap year issues
          
          if ( is_leap_year(date_time_tmp%year) .or. &
               is_leap_year(date_time%year)                  ) then

             if (date_time%month==2 .and. date_time%day==29)  &
                  date_time_tmp%day = 28
             
             
             call get_dofyr_pentad( date_time_tmp )
             
          end if
          
       end if
       
       ie = index( met_tag, 'day')
       
       if (ie>0) then
          
          read (met_tag(is+5:ie-1),'(i1)') shift
          
          call augment_date_time( 86400*shift, date_time_tmp )

       end if
    
    end if
    
    shift_forcing_date = date_time_tmp
    
  end function shift_forcing_date

  ! ******************************************************************
    
end module LDAS_ForceMod 

! -------------------------------------------------------------------

#if 0

program ut_parse_G5DAS_met_tag
  
  implicit none

  integer, parameter :: N_met_tag_in=6

  ! "in"
  character(200) :: met_path_in
  character( 80) :: met_tag_in
  
  ! "out"
  character(200) :: met_path_default, met_path_prec
  character( 80) :: met_tag_out
  
  logical        :: use_prec_corr
  logical        :: use_Predictor
  

  ! other
  integer :: ii
  character( 80), dimension(N_met_tag_in) :: met_tag_in_vec
  

  met_path_in = 'mymetpathin'
  
  met_tag_in_vec = (/                              &
       'gcmexpname' ,                          &
       'gcmexpname__Nx+-' ,                    &
       'gcmexpname__precCORRPREC',             &
       'gcmexpname__precCORRPREC__Nx+-' ,            &
       'gcmexpname__Nx+-__precCORRPREC' ,            &
       'gcmexpname__qrecCORRPREC__Nx+-'    /)            
  

  do ii=1,N_met_tag_in
     
     met_tag_in = met_tag_in_vec(ii)

     write (*,*) '----------------------------------------------'
     write (*,*) 'ii               = ', ii                  
     write (*,*) 'met_tag_in       = ', trim(met_tag_in)          

     call parse_G5DAS_met_tag( met_path_in, met_tag_in, &
          met_path_default, met_path_prec, met_tag_out, use_prec_corr, use_Predictor )
          
     write (*,*) 'met_path_default = ', trim(met_path_default)    
     write (*,*) 'met_path_prec    = ', trim(met_path_prec)       
     write (*,*) 'met_tag_out      = ', trim(met_tag_out)         
     write (*,*) 'use_prec_corr    = ', use_prec_corr       
     write (*,*) 'use_Predictor    = ', use_Predictor
     
  end do
  
end program ut_parse_G5DAS_met_tag

#endif

#if 0

program test_shift_forcing_date
  
  use clsm_ensdrv_functions  
  use date_time_util
  
  type(date_time_type) :: date_time, date_time_tmp
  
  integer :: dtstep, iter, i
  
  character(80) :: met_tag

  ! -------
  
  met_tag = 'RedArk_OSSE_shift+3day'
  met_tag = 'Princeton_shift+3day'
  
  dtstep = 21600

  iter = 120
    
  date_time%year    =  1992        ! 4-digit year
  date_time%month   =     2        ! month in year
  date_time%day     =     1        ! day in month
  date_time%hour    =     3        ! hour of day
  date_time%min     =     0        ! minute of hour
  date_time%sec     =     0        ! seconds of minute
  date_time%pentad  = -9999        ! pentad of year
  date_time%dofyr   = -9999        ! day of year
  
  do i=1,iter
     
     call augment_date_time(dtstep, date_time)
     
     date_time_tmp = shift_forcing_date(met_tag, date_time)
     
     write (*,'(16i5)') date_time, date_time_tmp
     
     
  end do

end program test_shift_forcing_date

#endif

! *********** EOF **************************************************
