MODULE point_driver_functions
  
  USE LDAS_DateTimeMod, only: date_time_type

  USE LDAS_DriverTypes, only: met_force_type, assignment(=)
 
  USE netcdf

  implicit none      
 
  ! INCLUDE 'netcdf.inc'
 
  private
  public :: write_time_varying, create_indata_files

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
        integer,                    intent (in)     ::  N_cat
        character (*),              intent (in)     ::  filename
        integer, dimension(8)                       ::  date_time_values
        character (22)                              ::  time_stamp
        character (40)                              ::  MyName

        integer ::  ncid, CellID, status, d2(2), vid, i, Dim2ID

        write (*,*) 'using create_indata_files'

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
        !call check(nf90_def_var(ncid, 'sunang'  , NF90_FLOAT, d2, vid))  


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
    

    subroutine write_time_varying (local_id, filename, time, met_f,sunang)
        
        use netcdf
        integer,                    intent (in) ::  local_id(:)
        type (date_time_type),      intent (in) ::  time
        character (*),              intent (in) ::  filename
        type (met_force_type),      intent (in) ::  met_f(:)
        real, dimension (:    ),    intent (in) ::  sunang
        integer ::  ncid, status, new_vid,met_fVarId, curr_varid, len_time, time_dimid, tile_dimid, len_tile
        integer :: a,b,c,time_con, i
        character(len=99) :: char_a,char_b,char_c, dim_name
        real, allocatable                       :: var_val(:), empty_var(:), return_val(:)
        real, dimension(:,:), allocatable       :: pull_val, final_val
        
        call check(nf90_open (trim(filename), nf90_write, ncid))

        call check(nf90_inq_dimid (ncid,'time',time_dimid))
        call check(nf90_inquire_dimension (ncid,time_dimid,dim_name,len_time))
        call check(nf90_inq_dimid (ncid,'tile',tile_dimid))
        call check(nf90_inquire_dimension (ncid,tile_dimid,dim_name,len_tile))
        
        call check(nf90_inq_varid(ncid, 'Tair', curr_varid))
        var_val = met_f%Tair
        allocate(pull_val(1,len_tile))
        !write(*,*) 'before pull'
        if (len_time > 0)then
            call check(nf90_get_var(ncid, curr_varid, pull_val, start = (/1,len_time/),count = (/len_tile,1/)))
        else
            pull_val=-9999
        endif
        write(*,*) count(pull_val == -9998)
        write(*,*) size(local_id)
        if (count(pull_val==-9999) == 0) then
            ! al = -9999
            len_time = len_time + 1
        endif 
        !write (*,*) 'after pull'
        allocate(empty_var(len_tile))
        empty_var = 0
        empty_var(local_id) = var_val
        empty_var = merge(empty_var,pull_val(1,:),empty_var>0)
        allocate(final_val(1,len_tile))
        final_val(1,:) = empty_var
        !write(*,*) 'before put'
        call check(nf90_put_var(ncid, curr_varid, final_val, start = (/1,len_time/),count = (/len_tile,1/)))
        !write (*,*) 'after put'
        
        
        
        !call check(nf90_inq_varid(ncid, 'Qair', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%Qair)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'Psurf', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%Psurf)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'Rainf_C', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%Rainf_C)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'Rainf', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%Rainf)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'Snowf', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%Snowf)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'LWdown', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%LWdown)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'SWdown', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%SWdown)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'PARdrct', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%PARdrct)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'PARdffs', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%PARdffs)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'Wind', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%Wind)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'RefH', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, met_f%RefH)) !, start = (/N_Cat,time_con/)))

        !call check(nf90_inq_varid(ncid, 'sunang', curr_varid))
        !call check(nf90_put_var(ncid, curr_varid, sunang)) !, start = (/N_Cat,time_con/)))
        
        ! ADD SOMETHING TO KEEP TRACK OF TIME IN A MEANINGFUL WAY!!! 

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
