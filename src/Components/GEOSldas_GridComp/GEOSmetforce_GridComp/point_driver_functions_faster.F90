module daily_forcing_mod
    use LDAS_DriverTypes, only: global_force_type
    type(global_force_type) :: daily_force
end module daily_forcing_mod

MODULE point_driver_functions
  
  USE LDAS_DateTimeMod, only: date_time_type

  USE LDAS_DriverTypes, only: met_force_type, assignment(=)
 
  USE netcdf

  USE daily_forcing_mod

  implicit none      
 
  ! INCLUDE 'netcdf.inc'
 
  private
  public :: write_time_varying, create_indata_files, check, HANDLE_ERR

  contains

    SUBROUTINE check(istatus)
       USE netcdf
       IMPLICIT NONE
       INTEGER, INTENT (IN) :: istatus
       IF (istatus /= nf90_noerr) THEN
         write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
       END IF
    END SUBROUTINE check 
    
    subroutine create_indata_files (N_cat,filename)
        
        use netcdf
        use daily_forcing_mod

        integer,                    intent (in)     ::  N_cat
        character (*),              intent (in)     ::  filename
        integer, dimension(8)                       ::  date_time_values
        character (22)                              ::  time_stamp
        character (40)                              ::  MyName

        integer ::  ncid, CellID, status, d2(2), vid, i, Dim2ID

        !write (*,*) 'using create_indata_files'

        call check(nf90_create (filename, NF90_CLOBBER, ncid))

        call check(nf90_def_dim (ncid, 'tile', N_cat,          CellID))
        call check(nf90_def_dim (ncid, 'time', NF90_UNLIMITED, Dim2ID))

        d2(1) = CellID
        d2(2) = Dim2ID

        call check(nf90_def_var(ncid, 'Tair'    , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'Qair'    , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'Psurf'   , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'Rainf_C' , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'Rainf'   , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'Snowf'   , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'LWdown'  , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'SWdown'  , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'PARdrct' , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'PARdffs' , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'Wind'    , NF90_FLOAT, d2, vid))    
        call check(nf90_def_var(ncid, 'RefH'    , NF90_FLOAT, d2, vid))  
        call check(nf90_def_var(ncid, 'date_int', NF90_DOUBLE, d2, vid))  


        !  Global attributes
        
        call getenv ("MYNAME"        ,MYNAME        )
        call date_and_time(VALUES=date_time_values)
         
        write (time_stamp,'(i4.4,a1,i2.2,a1,i2.2,1x,a2,1x,i2.2,a1,i2.2,a1,i2.2)')      &
             date_time_values(1),'-',date_time_values(2),'-',date_time_values(3),'at', &
             date_time_values(5),':',date_time_values(6),':',date_time_values(7)
         
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'CreatedBy', trim(MYNAME)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'Date'   ,   trim(time_stamp))) 
        call check(nf90_enddef(ncid))
        call check(nf90_close (ncid))


    end subroutine create_indata_files


    subroutine write_time_varying (local_id, filename, time, met_f,sunang,num_steps_day_75min)
        
        use netcdf
        integer,                    intent (in) ::  local_id(:), num_steps_day_75min
        type (date_time_type),      intent (in) ::  time
        character (*),              intent (in) ::  filename
        type (met_force_type),      intent (in) ::  met_f(:)
        real, dimension (:    ),    intent (in) ::  sunang
        integer ::  ncid, status, new_vid,met_fVarId, curr_varid, len_time, time_dimid, tile_dimid, len_tile, next_time
        integer :: i, stat
        character(len=99) :: char_a,char_b,char_c, dim_name
        real, allocatable                       :: var_val(:), empty_var(:), return_val(:)
        real, dimension(:,:), allocatable       :: pull_val, final_val, final_val_trans
        double precision, dimension(:,:), allocatable       :: final_val_double, final_val_trans_double

        character(len=4) :: YYYY
        character(len=2) :: MM,DD,hh,min,ss
        character(len=14) :: time_str
        double precision  :: time_double

        ! concatenate the time and convert it to an integer
        write (YYYY,'(i4.4)') time%year
        write (MM,'(i2.2)') time%month
        write (DD,'(i2.2)') time%day
        write (hh,'(i2.2)') time%hour
        write (min,'(i2.2)') time%min
        write (ss,'(i2.2)') time%sec

        time_str = YYYY//MM//DD//hh//min//ss
        read (time_str,*) time_double
        call check(nf90_open (trim(filename), nf90_write, ncid))
        call check(nf90_inq_dimid (ncid,'time',time_dimid))
        call check(nf90_inquire_dimension (ncid,time_dimid,dim_name,len_time))
        next_time = len_time+1
        call check(nf90_inq_dimid (ncid,'tile',tile_dimid))
        call check(nf90_inquire_dimension (ncid,tile_dimid,dim_name,len_tile))
        
        allocate(final_val(num_steps_day_75min,len_tile))
        allocate(final_val_trans(len_tile,num_steps_day_75min))
        allocate(final_val_double(num_steps_day_75min,len_tile))
        allocate(final_val_trans_double(len_tile,num_steps_day_75min))
        final_val = -9999
        final_val_trans = -9999
        final_val_double = -9999
        final_val_trans_double = -9999
        
        call check(nf90_inq_varid(ncid, 'Tair', curr_varid))
        final_val = daily_force%Tair
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999

        !write(*,*) 'starting Qair'
        call check(nf90_inq_varid(ncid, 'Qair', curr_varid))
        final_val = daily_force%Qair
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        
        !write(*,*) 'psurf'
        call check(nf90_inq_varid(ncid, 'Psurf', curr_varid))
        final_val = daily_force%Psurf
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'rainfc'
        call check(nf90_inq_varid(ncid, 'Rainf_C', curr_varid))
        final_val = daily_force%Rainf_C
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'rainf'
        call check(nf90_inq_varid(ncid, 'Rainf', curr_varid))
        final_val = daily_force%Rainf
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'snowf'
        call check(nf90_inq_varid(ncid, 'Snowf', curr_varid))
        final_val = daily_force%Snowf
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'lwdonw'
        call check(nf90_inq_varid(ncid, 'LWdown', curr_varid))
        final_val = daily_force%LWdown
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'swdown'
        call check(nf90_inq_varid(ncid, 'SWdown', curr_varid))
        final_val = daily_force%SWdown
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'pardrct'
        call check(nf90_inq_varid(ncid, 'PARdrct', curr_varid))
        final_val = daily_force%PARdrct
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'pardffs'
        call check(nf90_inq_varid(ncid, 'PARdffs', curr_varid))
        final_val = daily_force%PARdffs
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'wind'
        call check(nf90_inq_varid(ncid, 'Wind', curr_varid))
        final_val = daily_force%Wind
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999
        !write(*,*) 'refh'
        call check(nf90_inq_varid(ncid, 'RefH', curr_varid))
        final_val = daily_force%RefH
        final_val_trans = transpose(final_val)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        final_val = -9999
        final_val_trans = -9999

        !write(*,*) 'adding date_int'
        call check(nf90_inq_varid(ncid, 'date_int', curr_varid))
        final_val_double = daily_force%date_int
        final_val_trans_double = transpose(final_val_double)
        call check(nf90_put_var(ncid, curr_varid, final_val_trans_double, start = (/1,next_time/),count = (/len_tile,num_steps_day_75min/)))
        !write(*,*) 'final_val_trans_double(1,1)'
        !write(*,*) final_val_trans_double(1,1)
        
        
        ! ADD SOMETHING TO KEEP TRACK OF TIME IN A MEANINGFUL WAY!!!  
        !call check(nf90_inq_varid(ncid, 'date_int', curr_varid))
        !var_val = time_double
        !empty_var = time_double
        !empty_var(local_id) = var_val
        !final_val_double(1,:) = time_double
        !call check(nf90_put_var(ncid, curr_varid, final_val_double, start = (/1,next_time/),count = (/len_tile,1/)))
        
        !write(*,*) 'after adding time'

        ! close the file
        call check(nf90_close(ncid))

    end subroutine write_time_varying
       
    subroutine VID (ncid, VNAME, new_vid) 
      
      integer, intent (in)      :: ncid
      character(*), intent (in) :: VNAME
      integer, intent(out)      :: new_vid
      integer                   :: status
 
      call check(nf90_inq_varid (ncid, trim(VNAME) ,new_vid))
      IF (STATUS .NE. NF90_NOERR) &
           CALL HANDLE_ERR(STATUS, trim(VNAME))  
      
    end subroutine VID

    SUBROUTINE HANDLE_ERR(STATUS, Line)
 
      INTEGER,      INTENT (IN) :: STATUS
      CHARACTER(*), INTENT (IN) :: Line
 
      IF (STATUS .NE. NF90_NOERR) THEN
         PRINT *, trim(Line),': ',NF90_STRERROR(STATUS)
         STOP 'Stopped'
      ENDIF

    END SUBROUTINE HANDLE_ERR
      
    
  end module point_driver_functions
