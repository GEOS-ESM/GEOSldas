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

    subroutine write_time_varying (N_Cat, filename, time, met_f, sunang)  
  
       use netcdf
       integer,                                  intent (in) :: N_Cat
       type (date_time_type),                    intent (in) :: time
       character (*),                            intent (in) :: filename
       type (met_force_type),                    intent (in) :: met_f(:)
       real,                  dimension (N_Cat), intent (in) :: sunang
       integer                                               :: ncid, status, new_vid
      
       integer :: a,b,c,time_con
       character(len=99) :: char_a,char_b,char_c

       a = time%day
       b = time%min

       write(unit=char_a,fmt=*)a
       write(unit=char_b,fmt=*)b

       char_c = trim(adjustl(char_a))//trim(adjustl(char_b))

       read(unit=char_c,fmt=*)c 
       time_con = c

       call check(nf90_open (trim(filename),nf90_write,ncid))
       write (6,*) 'write'
       write (6,*) 'the ncid is:'
       write (6,*) ncid
       write (6,*) 'the time is:'
       write (6,*) time_con
       call VID(ncid,'time',new_vid)
       call check(nf90_put_var(ncid,new_vid, time_con       ))
       write (6,*) 'the ncid is:'
       write (6,*) ncid
       write (6,*) 'the VID for time is:'
       ! write (6,*) VID(ncid,'time',new_vid)
       call VID(ncid,'Tair',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%Tair     ))
       call VID(ncid,'Qair',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%Qair     ))
       call VID(ncid,'Psurf',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%Psurf    ))
       call VID(ncid,'Rainf_C',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%Rainf_C  ))
       call VID(ncid,'Rainf',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%Rainf    ))
       call VID(ncid,'Snowf',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%Snowf    ))
       call VID(ncid,'LWdown',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%LWdown   ))
       call VID(ncid,'SWdown',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%SWdown   ))
       call VID(ncid,'PARdrct',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%PARdrct  ))
       call VID(ncid,'PARdffs',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%PARdffs  ))
       call VID(ncid,'Wind',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%Wind     ))
       call VID(ncid,'RefH',new_vid)
       call check(nf90_put_var(ncid,new_vid, met_f%RefH     ))
       call VID(ncid,'sunang',new_vid)
       call check(nf90_put_var(ncid,new_vid, sunang         ))
  
       call check(nf90_close (ncid))
  
     end subroutine write_time_varying


    subroutine create_indata_files (N_cat, first_time, filename)

      use netcdf
      integer, intent (in)      :: N_cat
      logical, intent (in)      :: first_time
      character (*), intent (in):: filename
      integer, dimension(8)     :: date_time_values
      character (22)            :: time_stamp
      character (40)            :: MyName

      integer :: ncid, CellID, status, d2(2), vid, i,Dim2ID

      ! cat_param
      ! ---------
     

      ! time_varing
      ! -----------
         
      call check(nf90_create (filename, NF90_CLOBBER, ncid))     
      call check(nf90_def_dim (ncid,'tile',   N_cat        ,CellID))
      call check(nf90_def_dim (ncid,'time', NF90_UNLIMITED ,Dim2ID))
      d2(1)  = CellID
      d2(2)  = Dim2ID 
      ! call check(nf90_def_dim(ncid, 'time'    , NF90_INT  , Dim2ID, vid))    
      ! call check(nf90_def_dim(ncid, 'lai'     , d2, vid            ))    
      ! call check(nf90_def_var(ncid, 'grn'     , NF90_FLOAT, d2, vid))    
      call check(nf90_def_var(ncid, 'time'    , NF90_FLOAT, d2, vid))
      write (6,*) 'the vid for time is:'
      write (6,*) vid
      write (6,*) 'the ncid is:'
      write (6,*) ncid 
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
      call check(nf90_def_var(ncid, 'sunang'  , NF90_FLOAT, d2, vid))  
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
