
module LDAS_ensdrv_functions
 
  use LDAS_DateTimeMod,                   ONLY:   &
       date_time_type,                            &
       get_dofyr_pentad,                          &
       augment_date_time,                         &
       is_leap_year

  use LDAS_ensdrv_Globals,                ONLY:   &
       logit,                                     &
       logunit

  use LDAS_ExceptionsMod,                 ONLY:   &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR
  
  implicit none
  
  private

  public :: get_io_filename
  public :: is_in_rectangle

  character(300), private :: tmpstring300

contains  
  
  ! ********************************************************************
  
  character(300) function get_io_filename( io_path, exp_id, file_tag,   &
       date_time, dir_name, ens_id, option, file_ext, no_subdirs )
    
    ! compose file name for input/output 
    !
    ! file name = io_path/dir_name/[ensXXXX]/Yyyyy/Mmm/
    !             "exp_id"."file_tag".[ensXXXX.]YYYYMMDD_HHMMz"file_ext"
    !    
    ! example: iopath = /disk1/output/run1/GLOBAL/  (incl exp_id and domain name)
    !
    ! example  file_name = iopath/rs/ens0001/Y2005/M05/
    !                       run1.ens0001.catch_ldas_rst.20050510_0600z.bin
    !
    ! NOTE: if ens_id<0 then ens_id_string==ensXXXX="ens_avg"
    !       if ens_id is not present ensXXXX='' (empty string)
    !
    ! NOTE: option=1: return "io_path/dir_name/[ensXXXX]/*.ext"
    !       option=2: return "io_path/dir_name/[ensXXXX]/Yyyyy/Mmm/*.YYYYMM.ext"
    !       option=3: return "io_path/dir_name/[ensXXXX]/Yyyyy/Mmm/*.YYYYpPP.ext"
    !       option=4: return "io_path/dir_name/[ensXXXX]/Yyyyy/Mmm/*.YYYYMMDD.ext"
    !       option=5: return "io_path/dir_name/[ensXXXX]/Yyyyy/Mmm/*.YYYYMMDD_HHMMz.ext"
    !
    ! optional arguments:
    ! 
    !  date_time                      date and time info (must be present unless option=1) 
    !  option       default=5         controls date/time directories and string
    !  dir_name     default='cat'     what type of output (eg "rs/", "rc_out/", "ana/")
    !  ens_id       default=''        (see note above)
    !  file_ext     default='.bin'    file name extension
    !  no_subdirs   default=.false.   if .true., omit all sub-directories after "io_path"
    !
    ! reichle,  2 Sep 2008 - overhaul for new dir and file name conventions
    ! reichle, 28 Apr 2020 - added optional "no_subdirs" to facilitate writing to ./scratch

    implicit none
    
    character(*)                   :: io_path
    character(*)                   :: exp_id, file_tag       ! (eg "catch_ldas_rst")
    
    type(date_time_type), optional :: date_time
    
    integer,              optional :: ens_id       
    integer,              optional :: option
    
    character(*),         optional :: dir_name      ! default = 'cat'
    character(*),         optional :: file_ext      ! default = '.bin'

    logical,              optional :: no_subdirs    ! default = .false.  

    character(len=*), parameter    :: Iam = 'get_io_filename'    

    ! locals
    
    integer        :: tmp_option

    character(300) :: tmp_string
    character( 40) :: tmp_dir_name, tmp_file_ext, date_time_string
    character(  8) :: ens_id_string
    character(  4) :: YYYY, MMDD, HHMM, tmpstring4
    character(  2) :: PP
    logical        :: tmp_no_subdirs   
 
    ! --------------------------------------------------------
    !
    ! initialize optional arguments
 
    if (present(option)) then
       tmp_option = option
    else
       tmp_option = 5
    end if

    if (present(dir_name)) then
       tmp_dir_name = dir_name
    else
       tmp_dir_name = 'cat'
    end if

    if (present(file_ext)) then
       tmp_file_ext = file_ext
    else
       tmp_file_ext = '.bin'
    end if

    if (present(no_subdirs)) then
       tmp_no_subdirs = no_subdirs
    else
       tmp_no_subdirs = .false.
    end if
    
    ! create date/time strings
    
    if (tmp_option==1) then    
       
       date_time_string = ''
       
    else
       
       if (present(date_time)) then
          
          write (YYYY,'(i4.4)') date_time%year
          write (MMDD,'(i4.4)') date_time%month*100 + date_time%day
          write (HHMM,'(i4.4)') date_time%hour*100  + date_time%min    
          
          if (tmp_option==3) then
             
             ! determine %pentad if out of range
             ! - this might happen if only %year/%month/%day were set in "date_time"
             ! - if %pentad is within range, assume that %year/%month/%day and %pentad
             !     are consistent and do nothing
             
             if (date_time%pentad<1 .or. date_time%pentad>73)  call get_dofyr_pentad(date_time)
             
             write (PP,'(i2.2)') date_time%pentad
             
          end if
             
       else

          tmpstring300 = 'get_io_filename(): need optional argument date_time'
          
          if (logit) write (logunit,*) tmpstring300
          write(6,*) tmpstring300
          write(0,*) tmpstring300
          stop

       end if
                    
       if     (tmp_option==2) then
          
          date_time_string =  '.' // YYYY // MMDD(1:2)
          
       elseif (tmp_option==3) then
          
          date_time_string =  '.' // YYYY //  'p' // PP
          
       elseif (tmp_option==4) then
          
          date_time_string =  '.' // YYYY // MMDD
          
       elseif (tmp_option==5) then
          
          date_time_string =  '.' // YYYY // MMDD // '_' // HHMM // 'z'
          
       end if
       
    end if
    
    ! create ens ID string 
    
    if (present(ens_id)) then
       
       if (ens_id<0) then
          
          ens_id_string = '.ens_avg'
       
       else
       
          write (tmpstring4,'(i4.4)') ens_id
          
          ens_id_string = '.ens' // tmpstring4 
       
       end if
    
    else
    
       ens_id_string = ''       

    end if
    
    ! compose output path
    
    tmp_string = trim(io_path) // '/'
    
    if (.not. tmp_no_subdirs) then
       
       tmp_string = trim(tmp_string) // trim(tmp_dir_name) // '/' //  & 
            trim(ens_id_string(2:8)) // '/' 
       
       if (tmp_option>1)  &
            tmp_string = trim(tmp_string) // '/Y'// YYYY // '/M'//MMDD(1:2) // '/' 
       
    end if
       
    ! append file name to path

    get_io_filename = trim(tmp_string) // trim(exp_id) //        &
         trim(ens_id_string) // '.' // trim(file_tag) //         &
         trim(date_time_string) // trim(tmp_file_ext)
    
  end function get_io_filename

  ! ****************************************************************
    
!  character(2) function int2char2(int_in)
!    
!    ! Generates a length-2 character from an integer
!    
!    implicit none
!    
!    integer :: int_in
!    
!    write(int2char2,  '(i2.2)') int_in
!    
!  end function int2char2

  ! *****************************************************************
  
!  character(2) function int2char4(int_in)
!    
!    ! Generates a length-4 character from an integer
!    
!    implicit none
!    
!    integer :: int_in
!    
!    write(int2char4,  '(i4.4)') int_in
!    
!  end function int2char4
  
  ! ********************************************************************  
    
!  character(7) function timetag(int_day,int_hour)
!    
!    ! Generates a character time tag from an integer day and integer hour
!    
!    implicit none
!    integer :: int_day, int_hour
!    
!    character(3) :: char_day
!    character(4) :: char_hour
!    
!    write(char_day,  '(i3.3)') int_day
!    write(char_hour, '(i4.4)') int_hour
!    
!    timetag = char_day//char_hour
!    
!    return
!    
!  end function timetag
  
  ! ********************************************************************
  
!  character(6) function my_date(month,year)
!    
!    ! Generates a character date tag from an integer year and integer month
!    
!    implicit none
!    integer month,year
!    
!    character(4) :: char_year
!    character(2) :: char_month
!    
!    write(char_year,  '(i4.4)') year
!    write(char_month, '(i2.2)') month
!    
!    my_date = char_year//char_month
!    
!    return
!    
!  end function my_date
  
  ! ******************************************************************

  logical function is_in_rectangle(                                  &
       this_lon, this_lat, ll_lon, ll_lat, ur_lon, ur_lat        )
    
    ! determine whether point (this_lon, this_lat) is in rectangle defined by
    ! ll_lon, ll_lat, ur_lon, ur_lat
    
    ! - reichle, 2 Aug 2011
    
    implicit none
    
    real :: this_lon, this_lat, ll_lon, ll_lat, ur_lon, ur_lat        
    
    if ( (this_lon >= ll_lon) .and.        &
         (this_lon <= ur_lon) .and.        &
         (this_lat >= ll_lat) .and.        & 
         (this_lat <= ur_lat)       )    then
       is_in_rectangle = .true. 
    else
       is_in_rectangle = .false.
    end if
    
  end function is_in_rectangle

  ! ******************************************************************
  
end module LDAS_ensdrv_functions


! ******************************************************************

! driver routines for testing

#if 0

program test_get_io_filename
  
  use date_time_util
  use clsm_ensdrv_functions
  
  implicit none
  
  type(date_time_type) :: date_time
  
  character(300) :: io_path
  character(40)  :: io_run, file_tag, dir_name, file_ext
  
  integer :: option, ens_id

  ! ---------------------------------
  
  io_path = './output/ens_prop/N_AMER/run1/'
  io_run  = 'run1'
  file_tag = 'ldas_driver_inputs'
  dir_name = 'rc_out'
  file_ext = '.nml'

  ens_id = 12

  date_time%year    =  1992        ! 4-digit year
  date_time%month   =    11        ! month in year
  date_time%day     =     1        ! day in month
  date_time%hour    =     3        ! hour of day
  date_time%min     =     0        ! minute of hour
  date_time%sec     =     0        ! seconds of minute
  date_time%pentad  = -9999        ! pentad of year
  date_time%dofyr   = -9999        ! day of year


  call get_dofyr_pentad(date_time)

  write (*,*) get_io_filename( io_path, io_run, &
       date_time, file_tag, file_ext=file_ext )
  
  do option=1,7
     
     
     write (*,*) option, get_io_filename( io_path, io_run, &
          date_time, file_tag, option=option, &
          file_ext=file_ext)
     
  end do

end program test_get_io_filename

#endif



! ****** EOF *******************************************************

