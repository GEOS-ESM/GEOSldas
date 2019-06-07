#include "MAPL_ErrLog.h"

module LDAS_InterpMod

  use ESMF
  use MAPL_mod
  use LDAS_DateTimeMod, only: date_time_type, datetime2_minus_datetime1
  use LDAS_DriverTypes, only: met_force_type, assignment (=)
  use LDAS_ensdrv_Globals, only: nodata_generic, nodata_tol_generic

  implicit none

  private

  public :: metforcing_tinterp

  real, parameter :: min_zth = 0.0001

contains

  subroutine metforcing_tinterp(                                                &
       lons,                                                                    &
       lats,                                                                    &
       zth,                                                                     &
       zenav,                                                                   &
       force_time_prv,                                                          &
       model_time_nxt,                                                          &
       force_dtstep,                                                            &
       mf_prv,                                                                  &
       mf_nxt,                                                                  &
       mf_ntp,                                                                  &
       AEROSOL_DEPOSITION,                                                      &
       rc                                                                       &
       )

    ! Interpolates the forcing data to current timestep.
    !
    ! model_time_nxt = date_time at the *end* of model integration time step
    !
    ! "mf"   = "met_force"
    !
    ! "mf_prv" = at prv forcing time
    ! "mf_nxt" = at nxt forcing time
    ! "mf_ntp" = at current ("interpolated") time
    !
    ! NOTE: time avg radiative fluxes for the interval between "prv"
    !       and "nxt" time must be stored in mf_prv
    !
    ! reichle, 14 May 2003
    ! reichle, 11 Sep 2007 - albedo changes
    ! reichle, 23 Feb 2009 - add ParDrct, ParDffs from MERRA forcing
    ! reichle,  6 Mar 2009 - added fractional day-of-year (fdofyr) for monthly interp
    !                      - deleted ParDrct, ParDffs after testing found no impact
    ! reichle, 20 Dec 2011 - reinstated "PARdrct", "PARdffs" for MERRA-Land file specs
    !                      - cleanup
    ! reichle, 26 Jul 2013 - revised GRN, LAI and albedo scaling parameter inputs

    implicit none

    real, intent(in) :: lats(:) ! tile lats
    real, intent(in) :: lons(:) ! tile lons
    ! fix potential inconsistency between zth and zenav owing to 300s time
    ! step used in MAPL_SunGetInsolation()
    ! in ==>inout
    real, intent(inout) :: zth(:) ! zenith angle (>=0)
    real, intent(inout) :: zenav(:) ! avg zenith angle (>=0) over forcing period
    type(date_time_type), intent(in) :: model_time_nxt
    type(date_time_type), intent(in) :: force_time_prv
    integer, intent(in) :: force_dtstep
    type(met_force_type), intent(in) :: mf_prv(:)
    type(met_force_type), intent(in) :: mf_nxt(:)
    type(met_force_type), intent(out) :: mf_ntp(:)
    integer,intent(in)             :: AEROSOL_DEPOSITION
    integer, optional, intent(out) :: rc

    real :: w
    integer :: secs_since_prv, secs_in_day
    integer :: n, nTiles
    real :: sunang_ntp, tmpreal
    character(len=*), parameter :: Iam = 'metforcing_tinterp'

    ! get secs_in_day from hh:mm:ss
    secs_in_day = model_time_nxt%hour*3600 + model_time_nxt%min*60 &
         + model_time_nxt%sec

    ! weight for forcing "states" interpolation
    ! (temperature, humidity, pressure, wind)
    secs_since_prv = datetime2_minus_datetime1( force_time_prv, model_time_nxt )

    ! use INTEGER DIVISION such that w changes from 0. to 1.
    ! halfway through the current forcing interval, that is,
    !
    !   w = 0.  if                  secs_since_prv <  force_dtstep/2
    !   w = 0.5 if                  secs_since_prv == force_dtstep/2
    !   w = 1.  if force_dtstep/2 < secs_since_prv <= force_dtstep
    !
    ! For example, using 15 min model time steps and hourly forcing,
    ! the time interpolation weights are as follows:
    !
    !   secs_since_prv:         900        1800        2700        3600
    !   w:                        0.          0.5         1.          1.
    !
    ! Note that w=0.5 for secs_since_prv==force_dtstep/2 (at the mid-point).
    if (secs_since_prv==force_dtstep/2) then
       w = 0.5
    else
       w = real( (secs_since_prv-1)/(force_dtstep/2) )
    end if

    nTiles = size(mf_prv)

    do n=1,nTiles

       ! initialize
       mf_ntp(n) = nodata_generic

       ! STATES
       !
       ! temperature, humidity, pressure and wind

       mf_ntp(n)%Tair  = (1.-w)*mf_prv(n)%Tair  + w*mf_nxt(n)%Tair
       mf_ntp(n)%Qair  = (1.-w)*mf_prv(n)%Qair  + w*mf_nxt(n)%Qair
       mf_ntp(n)%Psurf = (1.-w)*mf_prv(n)%Psurf + w*mf_nxt(n)%Psurf
       mf_ntp(n)%RefH  = (1.-w)*mf_prv(n)%RefH  + w*mf_nxt(n)%RefH

       ! Wind

       ! LDASsa CVS tags between reichle-LDASsa_m2-10 and reichle-LDASsa_m2-13_p2
       !  worked with *inst*lfo* and *tavg*lfo* G5DAS forcing (as opposed to just
       !  *tavg*lfo* from MERRA) and G5DAS Wind was read from *inst*lfo* files.
       !  But for G5DAS forcing (as for MERRA), Wind was treated as a time-average
       !  field, which implied that, e.g., the instantaneous G5DAS Wind at 0z was
       !  used to force the land at 0:20z, 0:40z, and 1z.
       ! In these tags, Wind was treated as a time-average field whenever
       !  force_dtstep<=3601, which included MERRA, G5DAS, RedArk, and CONUS
       !  forcing.
       ! To make things more consistent for G5DAS winds, the "if" statement now
       !  checks whether Wind at date_time_nxt is available (which is true G5DAS
       !  and false for MERRA).  The revised "if" statement does not change how
       !  Wind is interpolated in MERRA (because Wind is unavailable at date_time_nxt)
       !  and in GLDAS, Princeton, and GWSP (because force_dtstep>3601).
       !  The revised "if" statement does change how CONUS and RedArk Wind data
       !  are interpolated (now treated as instantaneous, previously treated as
       !  time-average).
       ! In summary, as of this change, Wind is treated as instantaneous for all
       !  forcing data *except* MERRA.
       !
       ! - reichle, 31 Jan 2014

       ! if (force_dtstep_real > 3601.) then
       if (abs(mf_nxt(n)%Wind-nodata_generic)>nodata_tol_generic) then
          ! treat Wind as instantaneous fields (all forcing data sets *except* MERRA)
          mf_ntp(n)%Wind  = (1.-w)*mf_prv(n)%Wind  + w*mf_nxt(n)%Wind
       else
          ! treat Wind as time-average fields (MERRA)
          mf_ntp(n)%Wind  = mf_prv(n)%Wind
       end if

       ! FLUXES

       ! precipitation
       mf_ntp(n)%Rainf_C = mf_prv(n)%Rainf_C
       mf_ntp(n)%Rainf   = mf_prv(n)%Rainf
       mf_ntp(n)%Snowf   = mf_prv(n)%Snowf

       ! incoming radiation
       mf_ntp(n)%LWdown = mf_prv(n)%LWdown

       ! changed min sun-angle from 0.01 to 0.0001 for consistency with CatchGridComp
       ! reichle, 23 Feb 2009
       ! changed min sun-angle back to 0.01 after testing
       ! reichle, 27 Feb 2009
       sunang_ntp = max(zth(n), min_zth)

       ! changed minimum SWdown to 0. from 0.00001 - reichle, 28 Aug 2008
       ! fix potential inconsistency between zth and zenav owing to 300s time
       ! step used in MAPL_SunGetInsolation()
       if (abs(zenav(n)) <= 0.000001) then
          zth(n) = 0.
          zenav(n) = 0.
       endif

       if (zth(n) > 0.) then

          if (zenav(n) <= 0.) then
             RETURN_(ESMF_FAILURE)
          end if

          tmpreal = zth(n)/zenav(n)
          mf_ntp(n)%SWdown = mf_prv(n)%SWdown*tmpreal

          ! mf_ntp%SWnet only used if mf_prv%SWnet is not no-d ata value;
          ! protect multiplication for any no-data-value because it could
          ! fail (floating point excess) if huge number is used as nodata value
          if(abs(mf_prv(n)%SWnet  -nodata_generic)>nodata_tol_generic)  &
               mf_ntp(n)%SWnet   = mf_prv(n)%SWnet*tmpreal
          if(abs(mf_prv(n)%PARdrct-nodata_generic)>nodata_tol_generic) then
             ! assume that PARdffs is available whenever PARdrct is
             mf_ntp(n)%PARdrct = mf_prv(n)%PARdrct*tmpreal
             mf_ntp(n)%PARdffs = mf_prv(n)%PARdffs*tmpreal
          end if

       elseif ((zth(n) <= 0.) .and. (zenav(n) <= 0.)) then

          mf_ntp(n)%SWdown  = max(0., mf_prv(n)%SWdown)
          mf_ntp(n)%SWnet   = max(0., mf_prv(n)%SWnet)    ! no-data handling done below
          mf_ntp(n)%PARdrct = max(0., mf_prv(n)%PARdrct)  ! no-data handling done below
          mf_ntp(n)%PARdffs = max(0., mf_prv(n)%PARdffs)  ! no-data handling done below

       else

          mf_ntp(n)%SWdown  = 0.
          mf_ntp(n)%SWnet   = 0.   ! no-data handling done below
          mf_ntp(n)%PARdrct = 0.   ! no-data handling done below
          mf_ntp(n)%PARdffs = 0.   ! no-data handling done below

       end if

       ! cap shortwave radiation at (cosine of) sun angle times solar constant
       ! reichle, 14 Aug 2002
       mf_ntp(n)%SWdown = min( mf_ntp(n)%SWdown, 1360.*sunang_ntp )

       ! cap SWnet at SWdown
       mf_ntp(n)%SWnet  = min( mf_ntp(n)%SWnet, mf_ntp(n)%SWdown )

       ! reinstate no-data-values
       if(abs(mf_prv(n)%SWnet-nodata_generic)<nodata_tol_generic)  &
            mf_ntp(n)%SWnet = nodata_generic
       if ( (abs(mf_prv(n)%PARdrct-nodata_generic)<nodata_tol_generic) ) then
          ! assume that PARdffs is no-data whenever PARdrct is
          mf_ntp(n)%PARdrct = nodata_generic
          mf_ntp(n)%PARdffs = nodata_generic
       end if

    end do

    if (AEROSOL_DEPOSITION /=0 .and. nTiles >=1 ) then
       mf_ntp%DUDP001     =    mf_prv%DUDP001
       mf_ntp%DUDP002     =    mf_prv%DUDP002
       mf_ntp%DUDP003     =    mf_prv%DUDP003
       mf_ntp%DUDP004     =    mf_prv%DUDP004
       mf_ntp%DUDP005     =    mf_prv%DUDP005
       mf_ntp%DUSV001     =    mf_prv%DUSV001
       mf_ntp%DUSV002     =    mf_prv%DUSV002
       mf_ntp%DUSV003     =    mf_prv%DUSV003
       mf_ntp%DUSV004     =    mf_prv%DUSV004
       mf_ntp%DUSV005     =    mf_prv%DUSV005
       mf_ntp%DUWT001     =    mf_prv%DUWT001
       mf_ntp%DUWT002     =    mf_prv%DUWT002
       mf_ntp%DUWT003     =    mf_prv%DUWT003
       mf_ntp%DUWT004     =    mf_prv%DUWT004
       mf_ntp%DUWT005     =    mf_prv%DUWT005
       mf_ntp%DUSD001     =    mf_prv%DUSD001
       mf_ntp%DUSD002     =    mf_prv%DUSD002
       mf_ntp%DUSD003     =    mf_prv%DUSD003
       mf_ntp%DUSD004     =    mf_prv%DUSD004
       mf_ntp%DUSD005     =    mf_prv%DUSD005
       mf_ntp%BCDP001     =    mf_prv%BCDP001
       mf_ntp%BCDP002     =    mf_prv%BCDP002
       mf_ntp%BCSV001     =    mf_prv%BCSV001
       mf_ntp%BCSV002     =    mf_prv%BCSV002
       mf_ntp%BCWT001     =    mf_prv%BCWT001
       mf_ntp%BCWT002     =    mf_prv%BCWT002
       mf_ntp%BCSD001     =    mf_prv%BCSD001
       mf_ntp%BCSD002     =    mf_prv%BCSD002
       mf_ntp%OCDP001     =    mf_prv%OCDP001
       mf_ntp%OCDP002     =    mf_prv%OCDP002
       mf_ntp%OCSV001     =    mf_prv%OCSV001
       mf_ntp%OCSV002     =    mf_prv%OCSV002
       mf_ntp%OCWT001     =    mf_prv%OCWT001
       mf_ntp%OCWT002     =    mf_prv%OCWT002
       mf_ntp%OCSD001     =    mf_prv%OCSD001
       mf_ntp%OCSD002     =    mf_prv%OCSD002
       mf_ntp%SUDP003     =    mf_prv%SUDP003
       mf_ntp%SUSV003     =    mf_prv%SUSV003
       mf_ntp%SUWT003     =    mf_prv%SUWT003
       mf_ntp%SUSD003     =    mf_prv%SUSD003
       mf_ntp%SSDP001     =    mf_prv%SSDP001
       mf_ntp%SSDP002     =    mf_prv%SSDP002
       mf_ntp%SSDP003     =    mf_prv%SSDP003
       mf_ntp%SSDP004     =    mf_prv%SSDP004
       mf_ntp%SSDP005     =    mf_prv%SSDP005
       mf_ntp%SSSV001     =    mf_prv%SSSV001
       mf_ntp%SSSV002     =    mf_prv%SSSV002
       mf_ntp%SSSV003     =    mf_prv%SSSV003
       mf_ntp%SSSV004     =    mf_prv%SSSV004
       mf_ntp%SSSV005     =    mf_prv%SSSV005
       mf_ntp%SSWT001     =    mf_prv%SSWT001
       mf_ntp%SSWT002     =    mf_prv%SSWT002
       mf_ntp%SSWT003     =    mf_prv%SSWT003
       mf_ntp%SSWT004     =    mf_prv%SSWT004
       mf_ntp%SSWT005     =    mf_prv%SSWT005
       mf_ntp%SSSD001     =    mf_prv%SSSD001
       mf_ntp%SSSD002     =    mf_prv%SSSD002
       mf_ntp%SSSD003     =    mf_prv%SSSD003
       mf_ntp%SSSD004     =    mf_prv%SSSD004
       mf_ntp%SSSD005     =    mf_prv%SSSD005
    endif

    RETURN_(ESMF_SUCCESS)

  end subroutine metforcing_tinterp

end module LDAS_InterpMod
