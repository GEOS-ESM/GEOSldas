
module LDAS_DateTimeMod 
  
  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: date_time_type
  public :: date_time_print
  public :: is_leap_year
  public :: days_in_month
  public :: get_dofyr_pentad
  public :: augment_date_time
  public :: datetime2_minus_datetime1
  public :: datetime_eq_refdatetime
  public :: datetime_le_refdatetime
  public :: datetime_lt_refdatetime
  public :: date_time2string
  public :: datetime_to_J2000seconds
  public :: J2000seconds_to_datetime
  
  ! ---------------------------------------------------------------------
  
  ! WARNING: do not confuse date_time_type with the f90 intrisic
  !          function "date_and_time()"
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: date_time_type              
     integer :: year               ! 4-digit year
     integer :: month              ! month in year
     integer :: day                ! day in month
     integer :: hour               ! hour of day
     integer :: min                ! minute of hour
     integer :: sec                ! seconds of minute
     integer :: pentad             ! pentad of year
     integer :: dofyr              ! day of year
  end type date_time_type

contains  
  ! **********************************************************

  function date_time_print(dt) result (dtstr)

    implicit none

    type(date_time_type), intent(in) :: dt
    character(len=19) :: dtstr ! output

    write(dtstr,298) dt%year, dt%month, dt%day, dt%hour, dt%min, dt%sec
298 format(i4.4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2)

  end function date_time_print 
  ! **********************************************************
  
  integer function days_in_month(year, month)
    
    ! return the number of days in a given month
    
    implicit none
    
    integer :: year, month
    
    integer, dimension(12), parameter :: days_in_month_leap = &
         (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    
    integer, dimension(12), parameter :: days_in_month_nonleap = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    
    if (is_leap_year(year)) then
       days_in_month = days_in_month_leap(month) 
    else
       days_in_month = days_in_month_nonleap(month) 
    end if
    
  end function days_in_month
    
  ! ------------------------------------------------------------------
  
  integer function days_in_year(year)
    
    ! return the number of days in a given year
    
    implicit none
    
    integer :: year
    
    if (is_leap_year(year)) then
       days_in_year = 366
    else
       days_in_year = 365
    end if
    
  end function days_in_year
  
  ! ------------------------------------------------------------------
  
  logical function is_leap_year(year)
    
    implicit none

    integer :: year
    
    if (mod(year,4) /= 0) then 
       is_leap_year = .false.
    else if (mod(year,400) == 0) then
       is_leap_year = .true.
    else if (mod(year,100) == 0) then 
       is_leap_year = .false.
    else
       is_leap_year = .true.
    end if
    
  end function is_leap_year

  ! ------------------------------------------------------------------
  
  integer function pentad_of_year(day_of_year, year)
    
    implicit none
    
    integer :: day_of_year, year
    
    ! determine pentad
    
    if ((is_leap_year(year)) .and. day_of_year>=59) then
       
       pentad_of_year = (day_of_year-2)/5+1
       
    else
       
       pentad_of_year = (day_of_year-1)/5+1
       
    end if
    
  end function pentad_of_year
  
  ! ------------------------------------------------------------------
    
  subroutine get_dofyr_pentad( date_time ) 
    
    ! compute dofyr and pentad for date_time
    
    implicit none
    
    type(date_time_type), intent(inout) :: date_time
    
    integer :: i    ! local
    
    ! -------------------------------------------------------------
    
    date_time%dofyr = date_time%day
    
    do i=1,date_time%month-1 
       date_time%dofyr = date_time%dofyr + days_in_month(date_time%year,i)
    end do
    
    date_time%pentad = pentad_of_year(date_time%dofyr, date_time%year)
    
  end subroutine get_dofyr_pentad
  
  ! **********************************************************
  
  subroutine augment_date_time( dtstep, date_time )
    
    ! add dtstep (in seconds) to date_time
    !
    ! hh:mm:ss from 00:00:00 to 23:59:59
    !
    ! pentad runs from 1 to 73 (pentad 12 is 6 days long in leap years)
    ! 
    ! for convenience, use "secs_in_day" - this slightly awkward code already
    ! existed at time of writing
    ! maybe one day can write directly "date_time = date_time + dtstep"
    !
    ! reichle, 15 May 2003
    ! reichle, 22 Jun 2005 - added check on dtstep
    ! reichle, 28 Jun 2005 - arbitrary non-negative dtstep
    ! reichle, 25 Jul 2005 - arbitrary dtstep (incl negative)
    !
    ! ------------------------------------------------------------
    
    implicit none
    
    integer, intent(in) :: dtstep
    
    type(date_time_type), intent(inout) :: date_time
    
    ! local
    
    integer :: last_day, secs_in_day, secs_in_day_tmp, dtstep_left, dtstep_tmp
    
    ! ------------------------------------------------------------

    dtstep_left = dtstep

    if (dtstep==0) then
       
       return

    elseif (dtstep>0) then
       
       do while (dtstep_left>0)
          
          ! increase by one day at a time
          
          dtstep_tmp  = min( dtstep_left, 86400 )
          
          dtstep_left = dtstep_left - dtstep_tmp
          
          ! find out how many days in month
          
          last_day = days_in_month( date_time%year, date_time%month) 
          
          ! compute secs_in_day from hh:mm:ss
          
          secs_in_day = date_time%hour*3600 + date_time%min*60 + date_time%sec
          
          ! augment 
          
          secs_in_day = secs_in_day + dtstep_tmp
          
          ! compute new hh:mm:ss from secs_in_day
          
          date_time%hour = mod(secs_in_day,86400)/3600
          date_time%min  = mod(secs_in_day,86400)/60 - date_time%hour*60
          date_time%sec  = mod(secs_in_day,86400) - date_time%hour*3600     &
               - date_time%min*60
          
          ! augment year/month/day and dofyr as necessary
          
          if (secs_in_day>=86400) then 
             
             if (date_time%day==last_day) then
                
                if (date_time%month==12) then
                   
                   date_time%year  = date_time%year + 1
                   date_time%month = 1
                   date_time%day   = 1
                   
                else
                   
                   date_time%month = date_time%month + 1
                   date_time%day   = 1
                   
                end if
                
             else
                
                date_time%day   = date_time%day   + 1
                
             end if
             
          end if
          
       end do
       
    else               ! negative dtstep
       
       do while (dtstep_left<0)
          
          ! decrease by one day at a time
          
          dtstep_tmp  = max( dtstep_left, -86400 )
          
          dtstep_left = dtstep_left - dtstep_tmp
          
          ! compute secs_in_day from hh:mm:ss
          
          secs_in_day = date_time%hour*3600 + date_time%min*60 + date_time%sec
          
          ! augment 
          
          secs_in_day = secs_in_day + dtstep_tmp
          
          ! compute new hh:mm:ss from secs_in_day
          
          secs_in_day_tmp = secs_in_day + 86400;
          
          date_time%hour = mod(secs_in_day_tmp,86400)/3600
          date_time%min  = mod(secs_in_day_tmp,86400)/60 - date_time%hour*60
          date_time%sec  = mod(secs_in_day_tmp,86400) - date_time%hour*3600  &
               - date_time%min*60
          
  
          ! augment year/month/day and dofyr as necessary
          
          if ( secs_in_day < 0 ) then
             
             if (date_time%day==1) then
                
                if (date_time%month==1) then
                   
                   date_time%year  = date_time%year - 1
                   date_time%month = 12
                   date_time%day   = 31
                   
                else
	  
                   date_time%month = date_time%month - 1
	  
                   ! get number of days in previous month  
                   
                   date_time%day = &
                        days_in_month(date_time%year, date_time%month)
                   
                end if
                
             else
                
                date_time%day   = date_time%day - 1
                
             end if
             
          end if
          
       end do
       
    end if
    
    ! determine dofyr and pentad
    
    call get_dofyr_pentad( date_time )

  end subroutine augment_date_time
  
  ! ***********************************************************
  
  integer function datetime2_minus_datetime1( d1, d2 )
    
    ! compute difference between date_time_2 and date_time_1 in seconds
    !
    ! this function relies on the proper setting of date_time%dofyr !!
    ! (use get_dofyr_pentad() to initialize, note that augment_date_time()
    !  also takes care of dofyr)
    !
    ! reichle, 6 Jun 2005
    !
    ! -----------------------------------------------------------------
    
    implicit none
    
    type(date_time_type), intent(in) :: d1, d2
    
    ! locals
    
    type(date_time_type) :: de, dl
    
    integer :: fac, secs, y, secs_in_year_de, secs_in_year_dl
    
    ! -------------------------------------------------------------------
    !
    ! make sure type integer is not out of range
    !
    ! [integer*4 only allows for differences of up to ~68 years]

    if ( (abs(d2%year-d1%year)+1) > ((huge(secs)/86400)/366) ) then
       
       write (*,*) 'datetime2_minus_datetime1(): integer out of range.'
       write (*,*) 'STOPPING.'
       stop
       
    end if
    
    ! -------------------------------------------------------------------
    !
    ! check whether d1 and d2 are equal
    
    if (datetime_eq_refdatetime( d2, d1)) then
       
       datetime2_minus_datetime1 = 0
       
    else
       
       ! find out which is earlier
       
       if (datetime_le_refdatetime( d2, d1 )) then
          
          de  = d2                 ! "earlier"
          dl  = d1                 ! "later"
          fac = -1
          
       else
          
          de  = d1                 ! "earlier"
          dl  = d2                 ! "later"
          fac =  1
          
       end if
       
       ! seconds from beginning of de%year (or dl%year) to de (or dl):
       
       secs_in_year_de = (de%dofyr-1)*86400 + de%hour*3600 + de%min*60 + de%sec
       
       secs_in_year_dl = (dl%dofyr-1)*86400 + dl%hour*3600 + dl%min*60 + dl%sec
       
       
       if (dl%year>de%year) then  
          
          ! count seconds from de until end of de%year
          
          secs = days_in_year(de%year)*86400 - secs_in_year_de
          
          ! count seconds in years between de%year and dl%year
          
          do y=de%year+1,dl%year-1
             
             secs = secs + days_in_year(y)*86400
             
          end do
          
          ! add seconds from beginning of dl%year to dl
          
          secs = secs + secs_in_year_dl
          
       else
          
          secs = secs_in_year_dl - secs_in_year_de
          
       end if
       
       ! attach proper sign to difference
       
       datetime2_minus_datetime1 = fac*secs
       
    end if
    
  end function datetime2_minus_datetime1
  
  ! ********************************************************************
  
  logical function datetime_eq_refdatetime( datetime, refdatetime )

    implicit none
    
    ! find out whether datetime is equal to refdatetime
    !
    ! reichle, 25 May 2005
    
    type(date_time_type), intent(in) :: datetime, refdatetime
    
    if (datetime2_gt_eq_lt_datetime1(datetime,refdatetime)==0) then
       
       datetime_eq_refdatetime = .true.
       
    else
       
       datetime_eq_refdatetime = .false.
       
    end if
    
  end function datetime_eq_refdatetime
  
  ! ********************************************************************

  logical function datetime_le_refdatetime( datetime, refdatetime )

    implicit none

    ! find out whether datetime is earlier than or equal to refdatetime
    !
    ! reichle, 20 May 2005
    
    type(date_time_type), intent(in) :: datetime, refdatetime
    
    if (datetime2_gt_eq_lt_datetime1(datetime,refdatetime)>=0) then
       
       datetime_le_refdatetime = .true.
       
    else
       
       datetime_le_refdatetime = .false.
       
    end if
    
  end function datetime_le_refdatetime
  
  ! ********************************************************************

  logical function datetime_lt_refdatetime( datetime, refdatetime )

    implicit none

    ! find out whether datetime is earlier than (and not equal to) refdatetime
    !
    ! reichle, 20 May 2005
    
    type(date_time_type), intent(in) :: datetime, refdatetime
    
    if (datetime2_gt_eq_lt_datetime1(datetime,refdatetime)>0) then
       
       datetime_lt_refdatetime = .true.
       
    else
       
       datetime_lt_refdatetime = .false.
       
    end if
    
  end function datetime_lt_refdatetime
  
  ! ********************************************************************
  
  integer function datetime2_gt_eq_lt_datetime1( d1, d2  )
    
    implicit none
    
    ! compute "date_time_2 - date_time_1" to find out which is earlier
    !
    !  date_time_1 < date_time_2 : datetime2_gt_eq_lt_datetime1 =  1
    !  date_time_1 = date_time_2 : datetime2_gt_eq_lt_datetime1 =  0
    !  date_time_1 > date_time_2 : datetime2_gt_eq_lt_datetime1 = -1
    !
    ! reichle, 20 May 2005
    !
    ! -----------------------------------------------------------------
    
    type(date_time_type), intent(in) :: d1, d2
    
    ! ------------------------------------------------------------------
    
    if (d1%year < d2%year) then
       
       datetime2_gt_eq_lt_datetime1 = 1
       
    elseif (d1%year > d2%year) then
       
       datetime2_gt_eq_lt_datetime1 = -1
       
    elseif (d1%month < d2%month) then
       
       datetime2_gt_eq_lt_datetime1 = 1
       
    elseif (d1%month > d2%month) then
       
       datetime2_gt_eq_lt_datetime1 = -1
       
    elseif (d1%day < d2%day) then
       
       datetime2_gt_eq_lt_datetime1 = 1
       
    elseif (d1%day > d2%day) then
       
       datetime2_gt_eq_lt_datetime1 = -1

    elseif (d1%hour < d2%hour) then
       
       datetime2_gt_eq_lt_datetime1 = 1
       
    elseif (d1%hour > d2%hour) then
       
       datetime2_gt_eq_lt_datetime1 = -1

    elseif (d1%min < d2%min) then
          
       datetime2_gt_eq_lt_datetime1 = 1
       
    elseif (d1%min > d2%min) then
       
       datetime2_gt_eq_lt_datetime1 = -1

    elseif (d1%sec < d2%sec) then
       
       datetime2_gt_eq_lt_datetime1 = 1
       
    elseif (d1%sec > d2%sec) then
       
       datetime2_gt_eq_lt_datetime1 = -1

    else
       
       datetime2_gt_eq_lt_datetime1 = 0

    end if
        
  end function datetime2_gt_eq_lt_datetime1

  ! ********************************************************************
  
  character(16) function date_time2string( date_time )
    
    ! Generate "YYYYMMDD_HHMMSSz" string from date_time structure 
        
    implicit none
    
    type(date_time_type) :: date_time
    
    character(4) :: char_year
    character(2) :: char_month
    character(2) :: char_day
    character(2) :: char_hour
    character(2) :: char_min
    character(2) :: char_sec
    
    write(char_year,  '(i4.4)') date_time%year
    write(char_month, '(i2.2)') date_time%month
    write(char_day,   '(i2.2)') date_time%day
    write(char_hour,  '(i2.2)') date_time%hour
    write(char_min,   '(i2.2)') date_time%min
    write(char_sec,   '(i2.2)') date_time%sec
    
    date_time2string = char_year // char_month // char_day // '_' &
         // char_hour // char_min // char_sec // 'z'
    
  end function date_time2string
  
  ! ********************************************************************
  
  real*8 function datetime_to_J2000seconds( date_time, epoch_id )

    ! reichle, 14 Jan 2014
    
    implicit none

    type(date_time_type) :: date_time

    character(4)         :: epoch_id

    integer              :: tmp_secs

    ! --------------------------------------------
        
    ! get integer seconds 
    
    tmp_secs = datetime2_minus_datetime1( J2000_epoch(epoch_id), date_time )
    
    ! convert to double precision floating point value 
    
    datetime_to_J2000seconds = real(tmp_secs,kind(0.0D0)) 

    ! correct for 0.816 milliseconds (if using "UT12" definition of J2000 Epoch)

    if (epoch_id=='UT12')  datetime_to_J2000seconds = datetime_to_J2000seconds - 0.816D0
    
  end function datetime_to_J2000seconds

  ! ********************************************************************
  
  type(date_time_type) function J2000seconds_to_datetime( J2000_seconds, epoch_id )
    
    ! reichle, 14 Jan 2014
    
    implicit none
    
    real*8               :: J2000_seconds

    type(date_time_type) :: date_time
    
    character(4)         :: epoch_id

    integer              :: tmp_secs

    ! --------------------------------------------

    date_time = J2000_epoch( epoch_id )
    
    tmp_secs = nint(J2000_seconds)
    
    if ( J2000_seconds > real(huge(tmp_secs),kind(0.0D0)) ) then
       
       write (*,*) 'J2000seconds_to_datetime(): J2000_seconds out of range.'
       write (*,*) 'STOPPING.'
       stop
       
    end if
    
    call augment_date_time( tmp_secs, date_time )
    
    J2000seconds_to_datetime = date_time
    
  end function J2000seconds_to_datetime

  ! ********************************************************************

  type(date_time_type) function J2000_epoch( epoch_id )
    
    ! reichle, 30 Jun 2015
    
    implicit none
    
    character(4) :: epoch_id

    character(len=*), parameter :: Iam = 'J2000_epoch'
 
    ! -------------------------------------
    !
    ! definition of J2000 epochs
    !
    ! "J2000 seconds" are elapsed seconds since J2000 Epoch, which is either
    !
    !   - "UT12":  11:58:55.816 on 1 Jan 2000 in Coordinated Universal Time (UTC), or
    !   - "TT12":  12:00:00.000 on 1 Jan 2000 in Terrestrial Time (TT), or
    !   - "UT00":  00:00:00.000 on 1 Jan 2000 in Coordinated Universal Time (UTC)
    !
    ! NOTE: Per SMAP L1C_TB data products specs document, SMAP time stamps use "UT12"
    !        but sample granules appear to be using "TT12".
    ! NOTE: Per Clara Draper (30 Jun 2015), the nc4 ASCAT soil moisture retrieval 
    !        product uses "UT00".
    
    type(date_time_type), parameter :: J2000_UT12 = date_time_type( &
         year    = 2000,                                            &
         month   =    1,                                            &
         day     =    1,                                            &
         hour    =   11,                                            &
         min     =   58,                                            &
         sec     =   55,                                            &  ! rounded down
         pentad  =    1,                                            &
         dofyr   =    1 )

    type(date_time_type), parameter :: J2000_TT12 = date_time_type( &
         year    = 2000,                                            &
         month   =    1,                                            &
         day     =    1,                                            &
         hour    =   12,                                            &
         min     =    0,                                            &
         sec     =    0,                                            & 
         pentad  =    1,                                            &
         dofyr   =    1 )

    type(date_time_type), parameter :: J2000_UT00 = date_time_type( &
         year    = 2000,                                            &
         month   =    1,                                            &
         day     =    1,                                            &
         hour    =    0,                                            &
         min     =    0,                                            &
         sec     =    0,                                            & 
         pentad  =    1,                                            &
         dofyr   =    1 )
    
    ! ----------------------------
    
    select case (epoch_id)
       
    case ('UT12')
       J2000_epoch = J2000_UT12
       
    case ('TT12')
       J2000_epoch = J2000_TT12
       
    case ('UT00')
       J2000_epoch = J2000_UT00

    case default
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown J2000 epoch_id')
       
    end select
    
  end function J2000_epoch
  
  ! ********************************************************************

end module LDAS_DateTimeMod 

! ******************************************************************

! driver routines for testing

#if 0

program test
  
  use date_time_util
  
  implicit none
  
  type(date_time_type) :: start_time, end_time, date_time, date_time_tmp
  
  integer :: diff

  real*8  :: J2000_seconds = -987303600.0D0

  character(4) :: J2000_epoch_id


  start_time%year    =  1992        ! 4-digit year
  start_time%month   =    11        ! month in year
  start_time%day     =     1        ! day in month
  start_time%hour    =     0        ! hour of day
  start_time%min     =     0        ! minute of hour
  start_time%sec     =     0        ! seconds of minute
  start_time%pentad  = -9999        ! pentad of year
  start_time%dofyr   = -9999        ! day of year
  
  end_time%year      =  2000        ! 4-digit year
  end_time%month     =     1        ! month in year
  end_time%day       =     1        ! day in month
  end_time%hour      =    11        ! hour of day
  end_time%min       =    58        ! minute of hour
  end_time%sec       =    55        ! seconds of minute
  end_time%pentad    = -9999        ! pentad of year
  end_time%dofyr     = -9999        ! day of year

  date_time%year     =  2014        ! 4-digit year
  date_time%month    =     1        ! month in year
  date_time%day      =    14        ! day in month
  date_time%hour     =    12        ! hour of day
  date_time%min      =    34        ! minute of hour
  date_time%sec      =    56        ! seconds of minute
  date_time%pentad   = -9999        ! pentad of year
  date_time%dofyr    = -9999        ! day of year


  write (*,*) huge(diff)

  
  call get_dofyr_pentad( start_time ) 
  call get_dofyr_pentad( end_time ) 
  call get_dofyr_pentad( date_time ) 
  
  write (*,*) start_time
  write (*,*) end_time
  
  write (*,*) datetime_lt_refdatetime(start_time,end_time)
  write (*,*) datetime2_minus_datetime1(start_time,end_time)

  start_time = end_time

  write (*,*) start_time
  write (*,*) end_time
  
  write (*,*) datetime_lt_refdatetime(start_time,end_time)
  write (*,*) datetime2_minus_datetime1(start_time,end_time)  
  
  call augment_date_time( 1, start_time )
  
  write (*,*) start_time
  write (*,*) end_time

  write (*,*) datetime_lt_refdatetime(start_time,end_time)
  write (*,*) datetime2_minus_datetime1(start_time,end_time)

  write (*,*) '----------------'
  write (*,*) start_time
  
  call augment_date_time( 86401, start_time )
  
  write (*,*) start_time
  
  call augment_date_time( 32*86401, start_time )
  
  write (*,*) start_time

  call augment_date_time( -1, start_time )
  
  write (*,*) start_time

  call augment_date_time( -3600, start_time )
  
  write (*,*) start_time

  call augment_date_time( -86402, start_time )
  
  write (*,*) start_time

  call augment_date_time( -30*86401, start_time )
  
  write (*,*) start_time
  
  ! ------------------------------

  J2000_epoch_id = 'UT12'

  write (*,*) '-------------------------------'

  write (*,*) 'J2000_epoch_id = ', J2000_epoch_id

  write (*,*) '-------------------------------'
  write (*,*) start_time
  write (*,*) datetime_to_J2000seconds( start_time, J2000_epoch_id )
  write (*,*) '-------------------------------'
  write (*,*) end_time
  write (*,*) datetime_to_J2000seconds( end_time,   J2000_epoch_id )
  write (*,*) '-------------------------------'
  write (*,*) date_time
  write (*,*) datetime_to_J2000seconds( date_time,  J2000_epoch_id )
  write (*,*) '-------------------------------'

  ! ------------------------------

  date_time_tmp = end_time

  write (*,*) end_time
  write (*,*) date_time_tmp
  write (*,*) date_time

  diff = datetime2_minus_datetime1( end_time, date_time )

  write (*,*) diff
  write (*,*) datetime_to_J2000seconds( date_time,  J2000_epoch_id )
  
  call augment_date_time( diff, date_time_tmp )
  
  write (*,*) date_time_tmp
  
  ! ----------------------------------

  write (*,*) '-------------------------------'
  write (*,*) J2000_seconds
  
  write (*,*) J2000seconds_to_datetime( J2000_seconds, J2000_epoch_id )
  write (*,*) '-------------------------------'
  

end program test

#endif

! ****** EOF *******************************************************

