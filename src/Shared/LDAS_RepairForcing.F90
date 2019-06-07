module RepairForcingMod

  use LDAS_DriverTypes, only: met_force_type
  use LDAS_TileCoordType, only: tile_coord_type
  use LDAS_ExceptionsMod, only: ldas_abort, LDAS_GENERIC_ERROR
  use MAPL_SatVaporMod, only: MAPL_EQsat
  use LDAS_ensdrv_Globals, only: logunit, nodata_generic, nodata_tol_generic
  use MAPL_ConstantsMod, only: stefan_boltzmann=>MAPL_STFBOL
  use LDAS_ensdrv_Globals, only: master_logit
  implicit none

  private

  public :: repair_forcing

  real, parameter :: SWDN_MAX = 1360. ! W/m2
  real, parameter :: LWDN_EMISS_MIN = 0.5   ! min effective emissivity for LWdown
  real, parameter :: LWDN_EMISS_MAX = 1.0   ! max effective emissivity for LWdown

contains

  subroutine repair_forcing( N_catd, met_force, echo, tile_coord, fieldname, &
       unlimited_Qair, unlimited_LWdown )

    ! check forcing for unphysical values, reset, print optional warning
    !
    ! if optional input "fieldname" is present, only the corresponding
    ! forcing fields will be repaired
    !
    ! reichle, 30 Mar 2004
    ! reichle+qliu,  8 Oct 2008 - added optional input "unlimited_Qair"
    ! reichle,      11 Feb 2009 - added optional input "unlimited_LWdown"
    ! reichle,       2 Dec 2009 - eliminated "trim(field)==" to avoid memory
    !                              leak due to Absoft 9 bug
    ! reichle,      24 Dec 2013 - limited number of warnings to N_tile_warn_max
    !                             STILL TO DO: fix format specification for warning
    !                                          statements to avoid line breaks within
    !                                          a given statement

    implicit none

    integer, intent(in) :: N_catd

    type(met_force_type), dimension(N_catd), intent(inout) :: met_force

    logical, intent(in), optional :: echo, unlimited_Qair, unlimited_LWdown

    type(tile_coord_type), dimension(:), pointer, optional :: tile_coord ! in

    character(*), intent(in), optional :: fieldname

    ! local variables

    ! turn warnings off if warnings have been printed for N_warn_max tiles

    integer, parameter :: N_tile_warn_max = 3

    ! min/max values for allowable range of forcing fields

    real, parameter :: min_Tair           =    190. ! [K]
    real, parameter :: max_Tair           =    340. ! [K]

    real, parameter :: max_PSurf          = 115000. ! [Pa]

    ! slack parameters beyond which warnings are printed (if requested)

    ! specific humidity is often a tiny bit above saturated specific humidity
    ! (but sometimes much larger...)
    !
    ! ALWAYS limit Qair <= Qair_sat
    !
    ! *warning* ONLY for Qair <= tmpfac*Qair_sat
    !
    ! tmpfac=1.02 corresponds approximately to relative humidity above 1.02

    real, parameter :: tmpfac_warn_Qair   = 1.2    ! 1.02

    real, parameter :: tmpadd_warn_SWnet  = 1.e-2                 ! [W/m2]
    real, parameter :: tmpadd_warn_PAR    = tmpadd_warn_SWnet

    real, parameter :: tmp_warn_Prec      = 3.e-10  ! [m/s]  (1.e-10m/s ~ 3mm/year)


    logical :: warn, unlimited_Qair_tmp, unlimited_LWdown_tmp, problem_tile

    integer :: i, kk

    real    :: tmp_LWdown, min_LWdown, max_LWdown, Qair_sat, tmp_maxPar

    character(10) :: SWDN_MAX_string
    character(10) :: min_Tair_string, max_Tair_string, max_Psurf_string
    character(13) :: tmpstr13a, tmpstr13b
    character(16) :: tile_id_str, tmpstr16
    character(40) :: field

    character(len=*), parameter :: Iam = 'repair_forcing'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------

    if (present(echo)) then
       if (present(tile_coord)) then
          warn = echo
       else
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'inconsistent optional args')
       end if
    else
       warn = .false.
    end if

    if (present(unlimited_Qair)) then
       unlimited_Qair_tmp = unlimited_Qair
    else
       unlimited_Qair_tmp = .false.
    end if

    if (present(unlimited_LWdown)) then
       unlimited_LWdown_tmp = unlimited_LWdown
    else
       unlimited_LWdown_tmp = .false.
    end if

    ! --------------------------------

    if (present(fieldname)) then
       field = fieldname
    else
       field = 'all'
    end if

    ! --------------------------------

    if (warn) then

       write (SWDN_MAX_string, '(f10.1)') SWDN_MAX

       write (min_Tair_string, '(f10.1)') min_Tair
       write (max_Tair_string, '(f10.1)') max_Tair
       write (max_PSurf_string,'(f10.1)') max_PSurf

    end if

    ! --------------------------------

    kk = 0   ! counter for number of problematic tiles

    do i=1,N_catd

       problem_tile = .false.

       if (warn) write (tile_id_str,'(i16)') tile_coord(i)%tile_id  ! convert integer to string

       ! precip and snowfall must be non-negative

       ! fix Rainf first, otherwise cannot use Rainf to constrain Rainf_C
       ! reichle, 11 Aug 2010

       if (field(1:3)=='all' .or. field(1:7)=='Rainf  '  ) then

          if ((warn) .and. (met_force(i)%Rainf < -tmp_warn_Prec)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%Rainf  ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)') 'repair_forcing: Rainf < 0. in tile ID ' //  &
                  tile_id_str // ': met_force(i)%Rainf = ' // tmpstr13a

             problem_tile=.true.

          end if

          met_force(i)%Rainf  = max( 0.,       met_force(i)%Rainf)

       end if

       if (field(1:3)=='all' .or. field(1:7)=='Rainf_C') then

          if ((warn) .and. (met_force(i)%Rainf_C < -tmp_warn_Prec)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%Rainf_C  ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)') 'repair_forcing: Rainf_C < 0. in tile ID ' //&
                  tile_id_str // ': met_force(i)%Rainf_C = ' // tmpstr13a

             problem_tile=.true.

          end if

          met_force(i)%Rainf_C  = max( 0.,       met_force(i)%Rainf_C)

          ! make sure convective precip is less than total precip

          if ((warn) .and. (met_force(i)%Rainf+tmp_warn_Prec < met_force(i)%Rainf_C)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%Rainf    ! convert real to string
             write (tmpstr13b,'(e13.5)') met_force(i)%Rainf_C  ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)') 'repair_forcing: Rainf < Rainf_C in tile ID ' //  &
                  tile_id_str // ': met_force(i)%Rainf = ' // tmpstr13a //              &
                  ', met_force(i)%Rainf_C = ' // tmpstr13b

             problem_tile=.true.

          end if

          met_force(i)%Rainf_C = min( met_force(i)%Rainf, met_force(i)%Rainf_C)

       end if

       if (field(1:3)=='all' .or. field(1:7)=='Snowf  '  ) then

          if ((warn) .and. (met_force(i)%Snowf < -tmp_warn_Prec)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%Snowf    ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)') 'repair_forcing: Snowf < 0. in tile ID ' //&
                  tile_id_str // ': met_force(i)%Snowf = ' // tmpstr13a

             problem_tile=.true.

          end if

          met_force(i)%Snowf  = max( 0.,       met_force(i)%Snowf)

       end if

       ! --------------------------------

       if (field(1:3)=='all' .or. field(1:7)=='Tair   ') then

          ! temperatures must not be too low or too high

          ! NOTE: "warn" is turned on when repair_forcing is called first
          !        time after the forcing has been read from files

          if ((warn) .and. (met_force(i)%Tair < 190.)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%Tair    ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)')                                                 &
                  'repair_forcing: Tair < '//min_Tair_string//' in tile ID ' //       &
                  tile_id_str // ': met_force(i)%Tair = ' // tmpstr13a

             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'Tair too low')

          end if

          if ((warn) .and. (met_force(i)%Tair > 340.)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%Tair    ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)')                                                 &
                  'repair_forcing: Tair > '//max_Tair_string//' in tile ID ' //       &
                  tile_id_str // ': met_force(i)%Tair = ' // tmpstr13a

             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'Tair too high')

          end if

       end if

       ! --------------------------------

       if (field(1:3)=='all' .or. field(1:7)=='Psurf  ') then

          ! surface air pressure must not be too high

          if ((warn) .and. (met_force(i)%Psurf > max_Psurf)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%PSurf    ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)')                                                  &
                  'repair_forcing: Psurf > '//max_PSurf_string//' in tile ID ' //      &
                  tile_id_str // ': met_force(i)%PSurf = ' // tmpstr13a

             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'Psurf too high')

          end if

       end if

       ! --------------------------------

       ! specific humidity must not be negative

       if (field(1:3)=='all' .or. field(1:7)=='Qair   ') then

          if ((warn) .and. (met_force(i)%Qair < 0.)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%Qair    ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)') 'repair_forcing: Qair < 0. in tile ID ' //      &
                  tile_id_str // ': met_force(i)%Qair = ' // tmpstr13a

             problem_tile=.true.

          end if

          met_force(i)%Qair = max( 0., met_force(i)%Qair )

          ! specific humidity should not exceed saturated specific humidity

          if (.not. unlimited_Qair_tmp) then

             Qair_sat = MAPL_EQsat(met_force(i)%Tair, met_force(i)%Psurf )

             if ((warn) .and. (met_force(i)%Qair > tmpfac_warn_Qair*Qair_sat)) then

                write (tmpstr13a,'(e13.5)') met_force(i)%Qair    ! convert real to string
                write (tmpstr13b,'(e13.5)') met_force(i)%Qair/Qair_sat

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: Qair > Qair_sat in tile ID ' // &
                     tile_id_str // ': met_force(i)%Qair = ' // tmpstr13a //              &
                     ', met_force(i)%Qair/Qair_sat = ' // tmpstr13b

                problem_tile=.true.

             end if

             met_force(i)%Qair = min(Qair_sat, met_force(i)%Qair)

          end if

       end if

       ! --------------------------------

       ! wind must be positive
       ! (zero wind creates problem in turbulence calculations;
       !  warn only if wind is negative)

       if (field(1:3)=='all' .or. field(1:7)=='Wind   ') then

          if ((warn) .and. (met_force(i)%Wind < 0.)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%Wind    ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)') 'repair_forcing: Wind < 0. in tile ID ' //&
                  tile_id_str // ': met_force(i)%Wind = ' // tmpstr13a

             problem_tile=.true.

          end if

          met_force(i)%Wind   = max( 0.0001,   met_force(i)%Wind)

       end if

       ! --------------------------------

       if (field(1:3)=='all' .or. field(1:7)=='LWdown ') then

          if (.not. unlimited_LWdown_tmp) then

             ! make sure radiation is between min and max

             tmp_LWdown = met_force(i)%Tair*met_force(i)%Tair
             tmp_LWdown = stefan_boltzmann*tmp_LWdown*tmp_LWdown

             min_LWdown = LWDN_EMISS_MIN*tmp_LWdown
             max_LWdown = LWDN_EMISS_MAX*tmp_LWdown

             if ((warn) .and. (met_force(i)%LWdown < min_LWdown)) then

                write (tmpstr13a,'(e13.5)') met_force(i)%LWdown    ! convert real to string
                write (tmpstr13b,'(e13.5)') min_LWdown

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: LWdown < min_LWdown in tile ID ' // &
                     tile_id_str // ': met_force(i)%LWdown = ' // tmpstr13a //                &
                     ', min_LWdown = ' // tmpstr13b

                problem_tile=.true.

             end if

             met_force(i)%LWdown = max( min_LWdown, met_force(i)%LWdown)

             if ((warn) .and. (met_force(i)%LWdown > max_LWdown)) then

                write (tmpstr13a,'(e13.5)') met_force(i)%LWdown    ! convert real to string
                write (tmpstr13b,'(e13.5)') max_LWdown

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: LWdown > max_LWdown in tile ID ' // &
                     tile_id_str // ': met_force(i)%LWdown = ' // tmpstr13a //                &
                     ', max_LWdown = ' // tmpstr13b

                problem_tile=.true.

             end if

             met_force(i)%LWdown = min( max_LWdown, met_force(i)%LWdown)

          end if

       end if

       ! -----------------------------------

       if (field(1:3)=='all' .or. field(1:7)=='SWdown ') then

          if ((warn) .and. (met_force(i)%SWdown < 0. )) then

             write (tmpstr13a,'(e13.5)') met_force(i)%SWdown    ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)') 'repair_forcing: SWdown < 0. in tile ID ' //&
                  tile_id_str // ': met_force(i)%SWdown = ' // tmpstr13a

             problem_tile=.true.

          end if

          met_force(i)%SWdown  = max( 0.,       met_force(i)%SWdown)

          ! shortwave must be less than solar constant
          ! (see also subroutine interpolate_to_timestep() for interpolated
          !  shortwave forcing)

          if ((warn) .and. (met_force(i)%SWdown > SWDN_MAX)) then

             write (tmpstr13a,'(e13.5)') met_force(i)%SWdown    ! convert real to string

             if (master_logit) &
             write (logunit,'(200A)') 'repair_forcing: SWdown > ' // SWDN_MAX_string // &
                  ' in tile ID ' //  tile_id_str // ': met_force(i)%SWdown = ' // tmpstr13a

             problem_tile=.true.

          end if

          met_force(i)%SWdown  = min( SWDN_MAX, met_force(i)%SWdown)

       end if

       ! -----------------------------------

       ! SWnet is no-data-value for most forcing data sets (except MERRA, G5DAS)

       if(abs(met_force(i)%SWnet-nodata_generic)>nodata_tol_generic) then

          if (field(1:3)=='all' .or. field(1:7)=='SWnet  ') then

             if ((warn) .and. (met_force(i)%SWnet < 0. )) then

                write (tmpstr13a,'(e13.5)') met_force(i)%SWnet    ! convert real to string

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: SWnet < 0. in tile ID ' //&
                     tile_id_str // ': met_force(i)%SWnet = ' // tmpstr13a

                problem_tile=.true.

             end if

             met_force(i)%SWnet  = max( 0.,       met_force(i)%SWnet)

             ! net solar radiation must be less than solar incoming radiation

             if ( (warn) .and. &
                  (met_force(i)%SWnet > met_force(i)%SWdown+tmpadd_warn_SWnet)) then

                write (tmpstr13a,'(e13.5)') met_force(i)%SWnet    ! convert real to string
                write (tmpstr13b,'(e13.5)') met_force(i)%SWdown   ! convert real to string

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: SWnet > SWdown in tile ID ' //  &
                     tile_id_str // ': met_force(i)%SWnet = ' // tmpstr13a //             &
                     ', met_force(i)%SWdown = ' // tmpstr13b

                problem_tile=.true.

             end if

             met_force(i)%SWnet  = min( met_force(i)%SWnet,met_force(i)%SWdown)

          end if
       end if

       ! -----------------------------------

       ! PARdffs is no-data-value for most forcing data sets (except MERRA, G5DAS)

       if(abs(met_force(i)%PARdffs-nodata_generic)>nodata_tol_generic) then

          if (field(1:3)=='all' .or. field(1:7)=='PARdffs') then

             if ((warn) .and. (met_force(i)%PARdffs < 0. )) then

                write (tmpstr13a,'(e13.5)') met_force(i)%PARdffs    ! convert real to string

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: PARdffs < 0. in tile ID ' //&
                     tile_id_str // ': met_force(i)%PARdffs = ' // tmpstr13a

                problem_tile=.true.

             end if

             met_force(i)%PARdffs  = max( 0.,       met_force(i)%PARdffs)

             ! upper threshold for PARdffs = 0.6 of solar incoming radiation
             ! (after analysis if May and Dec cases from MERRA Scout)
             ! - reichle, 24 Feb 2009)

             ! updated to threshold of 0.8 after brief analysis of hourly MERRA data
             ! for Jul 2003 (courtesy of Greg Walker)
             ! - reichle, 20 Dec 2011

             tmp_maxPar = 0.8*met_force(i)%SWdown

             if ((warn) .and. (met_force(i)%PARdffs > tmp_maxPar+tmpadd_warn_PAR )) then

                write (tmpstr13a,'(e13.5)') tmp_maxPar           ! convert real to string
                write (tmpstr13b,'(e13.5)') met_force(i)%PARdffs ! convert real to string

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: PARdffs > ' // tmpstr13a //   &
                     ' in tile ID ' // tile_id_str // ': met_force(i)%PARdffs = ' //    &
                     tmpstr13b

                problem_tile=.true.

             end if

             met_force(i)%PARdffs  = min( met_force(i)%PARdffs,tmp_maxPar)

          end if
       end if

       ! -----------------------------------

       ! PARdrct is no-data-value for most forcing data sets (except MERRA, G5DAS)
       !
       ! MUST "repair" PARdffs *before* PARdrct

       if(abs(met_force(i)%PARdrct-nodata_generic)>nodata_tol_generic) then

          if (field(1:3)=='all' .or. field(1:7)=='PARdrct') then

             if ((warn) .and. (met_force(i)%PARdrct < 0. )) then

                write (tmpstr13a,'(e13.5)') met_force(i)%PARdrct ! convert real to string

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: PARdrct < 0. in tile ID ' //&
                     tile_id_str // ': met_force(i)%PARdrct = ' // tmpstr13a

                problem_tile=.true.

             end if

             met_force(i)%PARdrct  = max( 0.,       met_force(i)%PARdrct)

             ! upper threshold for PARdrct = SWdown - PARdffs

             tmp_maxPar = met_force(i)%SWdown - met_force(i)%PARdffs

             if ((warn) .and. (met_force(i)%PARdrct > tmp_maxPar+tmpadd_warn_PAR )) then

                write (tmpstr13a,'(e13.5)') tmp_maxPar           ! convert real to string
                write (tmpstr13b,'(e13.5)') met_force(i)%PARdrct ! convert real to string

                if (master_logit) &
                write (logunit,'(200A)') 'repair_forcing: PARdrct > ' // tmpstr13a //   &
                     ' in tile ID ' // tile_id_str // ': met_force(i)%PARdrct = ' //    &
                     tmpstr13b

                problem_tile=.true.

             end if

             met_force(i)%PARdrct  = min( met_force(i)%PARdrct,tmp_maxPar)

          end if
       end if

       ! ------------------------------------------------------
       !
       ! count problematic tiles

       if (problem_tile) kk=kk+1

       ! turn off warnings if number of problem tiles gets too large

       if ((warn) .and. (kk>N_tile_warn_max)) then

          warn = .false.  ! turn off warnings for the remainder of the loop through tiles

          write (tmpstr16,'(i16)') kk  ! convert integer to string

          if (master_logit) &
          write (logunit,'(200A)')                                         &
               'repair_forcing: turning OFF warnings after detecting ' //  &
               trim(tmpstr16) // ' tiles with problematic forcing'

       end if

       ! ------------------------------------------------------

    end do

  end subroutine repair_forcing

end module RepairForcingMod
