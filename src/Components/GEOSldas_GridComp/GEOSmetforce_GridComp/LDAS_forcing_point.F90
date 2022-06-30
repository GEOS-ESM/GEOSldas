#include "MAPL_Generic.h"


module global_forcing_mod
    use LDAS_DriverTypes, only: global_force_type
    type(global_force_type) :: met_force_new
end module global_forcing_mod

module LDAS_Point_ForceMod
  
  use netcdf

  use point_driver_functions, only: check, HANDLE_ERR 
  use LDAS_DateTimeMod, only: date_time_type
  use LDAS_ExceptionsMod, only: ldas_abort, LDAS_GENERIC_ERROR

  implicit none

  private

  public :: get_forcing_point

contains

  subroutine get_forcing_point(met_path,date_time,local_id)
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! TWR NOTES:                                                        !!
    !!!! DON'T NEED TO DO ANY INTERPOLATION--THIS WILL ALREADY BE DONE   !!
    !!!!!! WHEN POINT FORCING FILES ARE WRITTEN                          !!
    !!!! 
    ! Format to write files into to make this easiest:
    !!! Tile space
    !!! Already interpolated
    ! Next steps:
    !!! Where does this get called in the model?
    !!! Change the line of code there and try to compile
    !!! Find the errors :)
    !!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    use netcdf
    use global_forcing_mod

    implicit none

    ! intent in:

    character(*),         intent(in) :: met_path !path to metforcing
    ! character(*),         intent(in) :: met_tag  !name of metforcing file
    type(date_time_type), intent(in) :: date_time! time
    integer,              intent(in) :: local_id(:)



    ! intent out:
    
    ! intent inout:
    
    !!!! THIS IS WHAT WILL BE RETURNED. THE GLOBAL VARIABLE THAT CAN BE STORED W/IN CATCHMENT BETWEEN STEPS
    !!!! I DO NOT KNOW WHAT TYPE THIS WILL BE -- UP TO SOFTWARE PEOPLE
    !!!! WILL BE 2D -- FIRST DIMENSION EVERY TIME STEP THAT YEAR, SECOND DIMENSION THE TILE NUMBER
    !!!! WILL NEED TO HAVE ONE OF THESE VARIABLES FOR EACH PARALELLIZED PIECE OF THE MODEL
    !!!!    AT LEAST THIS IS MY BEST UNDERSTANDING AT THIS POINT
    !type(met_force_type),  dimension(:), intent(out) :: met_force_new !array to write forcing data to
    
    ! local variables:

    character(4)                       :: YYYY !current year for which want data
    character(400)                     :: filename !name of file you are opening

    ! character( 10)             :: lat_str = 'latitude'
    ! character( 10)             :: lon_str = 'longitude'
    
    ! integer                    :: fid, rc, nv_id, status
    
    character(len=*), parameter :: Iam = 'read_point_forcing'
    integer                     :: ncid, ndim, nvar, natt, unlid, varid, len1, len2 ! identifiers used by nf90
    character(400)              :: varname ! name of variable, used by nf90
    double precision, allocatable, dimension(:,:)    :: values  ! values read in by nf90





    ! begin reading in the files

    ! check to make sure that you are reading in the correct file
    ! THIS MAY NEED TO BE CHANGED BASED OFF OF HOW I DECIDE TO INDEX TIME
    write(*,*) 'before abort call'
    if ( (date_time%min/=0) .or. (date_time%sec/=0) .or.       &
        (date_time%hour/=0)) then

        call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'input file incorrect!!')

    endif
    write(*,*) 'after abort call'

    !write (YYYY, '(i4.4)') date_time%year

    ! LEFT OFF AT LDAS_FORCING LINE 2829!!!!!!!!!!
    
    !filename = trim(met_path) // '/' // YYYY // 'point_forcings.nc4'
    !filename = trim(filename)
    
    ! NEED TO LOOK INTO HOW TO USE GEOS_openfile TO OPEN THIS FILE
    ! THIS IS THE NEXT STEP:
    !   DO I NEED TO USE GEOS_openfile??? I DO NOT NEED TO USE A LOT
    !       CANT USE GEOS_openfile!!!! (only for global forcing)
    !   OF THE STUFF IN JANA'S CODE BECAUSE MY VARIABLES HAVE ALREADY  
    !   BEEN INTERPOLATED IN TIME AND ARE ALREADY AT TILE SPACE.
    !   THEREFORE, WHAT IS THE EASIEST WAY TO PASS THESE TO THE REST
    !   OF THE MODEL????

    ! JUST USE THE NF90_ COMMANDS TO OPEN THE FILE AND ASSIGN THE
    ! VARIABLES TO SOMETHING THAT FORTRAN CAN STORE GLOBALLY!!!
    ! SIMPLE!!! ALL WE NEED TO DO!!!!

    !*****************************************************************************************!
    ! EVERYTHING HERE DOESN'T NEED TO BE USED. KEEPING HERE AS A REMINDER.
    ! call GEOS_openfile(FileOpenedHash,filename,fid,tile_coord,met_hinterp,rc,lat_str,lon_str)
    !
    ! if (rc<0) then
    !   call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'error opening file')
    ! endif
    !****************************************************************************************!

    ! I am going to write this assuming that we have:
        ! A netcdf file with all necessary variables
        ! each timestep is one model timestep
        ! index of the tilespace dimension is the tilespace value 
    write(*,*) 'nf90 open'
    call check(nf90_open(met_path, nf90_nowrite, ncid)) !now we have the ncid
    write(*,*) 'nf90 inquire'
    call check(nf90_inquire(ncid, ndim, nvar, natt, unlid))
    ! do varid = 1,nvar ! A good idea, but this won't work with how met_force_type works. must be done manually
    ! will set this up to work with Tair, then go from there
    write(*,*) 'inq dimension 1'
    call check(nf90_inquire_dimension(ncid,1,len=len1))
    write(*,*) 'inq dimension 2'
    call check(nf90_inquire_dimension(ncid,2,len=len2))
    write(*,*) 'allocate statements'
    allocate(met_force_new%Tair(len1,len2))
    allocate(met_force_new%Qair(len1,len2)) 
    allocate(met_force_new%Psurf(len1,len2))
    allocate(met_force_new%Rainf_C(len1,len2))
    allocate(met_force_new%Rainf(len1,len2))
    allocate(met_force_new%Snowf(len1,len2))
    allocate(met_force_new%LWdown(len1,len2))
    allocate(met_force_new%SWdown(len1,len2))  
    allocate(met_force_new%PARdrct(len1,len2))
    allocate(met_force_new%PARdffs(len1,len2))
    allocate(met_force_new%Wind(len1,len2))
    allocate(met_force_new%RefH(len1,len2))
    allocate(met_force_new%date_int(len1,len2))
    allocate(values(len1,len2))
    write(*,*) 'assigning values'
    write(*,*) 'finding varid'
    call check(nf90_inq_varid(ncid,'Tair',varid))
    write(*,*) 'getting values'
    call check(nf90_get_var(ncid, varid, values))
    write(*,*) 'putting in met_force_new'
    met_force_new%Tair = values
    write(*,*) 'done for Tair'
     
    call check(nf90_inq_varid(ncid,'Qair',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%Qair = values
    write(*,*) 'Qair' 
    call check(nf90_inq_varid(ncid,'Psurf',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%Psurf = values
    write(*,*) 'Psurf' 
    
    call check(nf90_inq_varid(ncid,'Rainf_C',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%Rainf_C = values
    write(*,*) 'Rainf_C'
    call check(nf90_inq_varid(ncid,'Rainf',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%Rainf = values
    write(*,*) 'Rainf'
    call check(nf90_inq_varid(ncid,'Snowf',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%Snowf = values
    write(*,*) 'Snowf'
    call check(nf90_inq_varid(ncid,'LWdown',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%LWdown = values
    write(*,*) 'LWdown'
    call check(nf90_inq_varid(ncid,'SWdown',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%SWdown = values
    write(*,*) 'SWdown'
    call check(nf90_inq_varid(ncid,'PARdrct',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%PARdrct = values
    write(*,*) 'PARdrct'
    call check(nf90_inq_varid(ncid,'PARdffs',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%PARdffs = values
    write(*,*) 'PARdffs'
    call check(nf90_inq_varid(ncid,'Wind',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%Wind = values
    write(*,*) 'Wind'
    call check(nf90_inq_varid(ncid,'RefH',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%RefH = values
    write(*,*) 'RefH'
    call check(nf90_inq_varid(ncid,'date_int',varid))
    call check(nf90_get_var(ncid, varid, values))
    met_force_new%date_int = values
    write(*,*) 'tim_int'
    write(*,*) values(1,1)
    write(*,*) met_force_new%date_int(1,1)
    ! Once debugged, implement this same procedure for all variables
    !   in met_force_new
    write(*,*) 'closing'
    call check(nf90_close(ncid))
    
    


    end subroutine get_forcing_point


end module LDAS_Point_ForceMod    


 






