#include "MAPL_Generic.h"

module clsm_ensdrv_out_routines
  
  ! collection of LDASsa output subroutines 
  !
  ! (originally in clsm_ensdrv_drv_routines.F90)
  !
  ! reichle, 22 Aug 2014

  use LDAS_ensdrv_globals,              ONLY:     &
       log_master_only,                           &
       logunit,                                   &
       logit
  
  use catch_types,                      ONLY:     &
       cat_param_type,                            &
       assignment (=),                            &
       operator (+),                              &
       operator (/)

  use LDAS_TileCoordType,               ONLY:     &
       tile_coord_type

  use mwRTM_types,                      ONLY:     &
       mwRTM_param_type

  use LDAS_ensdrv_mpi,                  ONLY:     &
       master_proc,                               &
       numprocs

  use LDAS_DateTimeMod,                 ONLY:    &
       date_time_type
  
  use LDAS_ensdrv_init_routines,        ONLY:     &
       clsm_ensdrv_get_command_line,              &
       add_domain_to_path

  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename
       
  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use, intrinsic :: iso_fortran_env, only : output_unit

  implicit none

  private

  public :: init_log             
  public :: GEOS_output_smapL4SMlmc   

contains
  
  ! ********************************************************************

  subroutine init_log( myid, numprocs, master_proc )
    
    ! open file for output log, write a few things

    ! changed logic so that error messages generated before init_log()
    ! has completed do not get lost
    ! - reichle, 29 Aug 2014

    implicit none
    
    integer, intent(in) :: myid, numprocs
    logical, intent(in) :: master_proc
    
    ! ------------------------------------------------------------------------
    !
    ! local variables

    type(date_time_type) :: start_time
    
    integer        :: istat
    
    character(300) :: fname
    character(200) :: work_path, io_path
    character(40)  :: exp_domain, exp_id, dir_name, file_tag, file_ext
    
    character(4)   :: myid_string
    
    character(8)   :: date_string
    character(10)  :: time_string
    
    character(len=*), parameter :: Iam = 'init_log'
    character(len=400) :: err_msg

    ! -----------------------------------------------------------------------
    
    ! interpret parameters from clsm_ensdrv_glob_param

    if (log_master_only .and. (.not. master_proc)) then
       
       logit = .false.
       
    else
       
       logit = .true.
       
    end if

    ! stop if logunit is stdout and output is requested for *all* processors
    
    if ( (.not. log_master_only) .and. (logunit==output_unit) ) then
       
       err_msg = 'logunit=output_unit (stdout) together with logging *all* procs is disabled'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    ! if requested, open a log file that is different from stdout
    ! (this would typically be needed if a separate log files 
    !  is requested for each processor)
    
    if (logit) then
       
       if (logunit/=output_unit) then
          
          ! get command line arguments
          !
          ! NOTE: If ldas_abort() is called from clsm_ensdrv_get_command_line() at this
          !       time, the error message should appear in a file named "fort.[logunit]"
          
          call clsm_ensdrv_get_command_line(                 &
               start_time=start_time,                        &
               work_path=work_path, exp_domain=exp_domain,   &
               exp_id=exp_id )         
          
          ! augment work_path (must be same as in read_driver_inputs() )
          
          io_path = add_domain_to_path( work_path, exp_domain )
          
          write (myid_string,'(i4.4)') myid
          
          dir_name = 'rc_out'
          file_tag = 'ldas_log' 
          file_ext = '.txt'
          
          if (.not. master_proc) then
             
             file_tag = trim(file_tag) // '_PE' // myid_string
             
          end if
          
          ! NOTE: If ldas_abort() is called from get_io_filename() at this time, 
          !       the error message should appear in a file named "fort.[logunit]"      
          
          fname = get_io_filename( io_path, exp_id, file_tag, date_time=start_time, &
               dir_name=dir_name, file_ext=file_ext )
          
          open (logunit, file=trim(fname), form='formatted', action='write',        &
               status='new', iostat=istat)
          
          if (istat/=0) then
             
             ! this call to ldas_abort() should create a file named "fort.[logunit]"
             
             err_msg = 'ERROR opening log file (perhaps it already exists)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             
          end if
          
          write (logunit,*)
          write (logunit,'(400A)') 'logfile: ' // trim(fname)
          
       else
          
          write (logunit,*)
          write (logunit,*) 'using stdout for log messages' 
          
       end if
       
       write (logunit,*) 
       write (logunit,*) 'Offline ensemble driver for Catchment model'
       write (logunit,*) 
       
       ! echo wall clock
       
       call date_and_time(date_string, time_string)
       
       write (logunit,*) 'started at ', date_string, ', ', time_string
       write (logunit,*)
       
       ! echo MPI environment
       
       write (logunit,*) "process ", myid, " of ", numprocs, " is alive"
       write (logunit,*)
       write (logunit,*) "process ", myid, ": master_proc=", master_proc
       write (logunit,*)
       
    end if  ! if (logit)
    
  end subroutine init_log
  
  ! ********************************************************************

  subroutine GEOS_output_smapL4SMlmc( GC, date_time, work_path, exp_id, &    
       N_catl, tile_coord_l, cat_param, mwRTM_param ) 
    
    ! write SMAP L4_SM "lmc" (land model constants) file collection
    ! as binary, tile-space output
    !
    ! requires "full" domain inputs
    !
    ! reichle, 26 Apr 2012
    ! reichle, 27 May 2014: - changed wilting point output from "clsm_wpwet" to "clsm_wp"
    ! 
    ! -------------------------------------------------------------------
    use ESMF
    USE MAPL_MOD 

    implicit none
    type(ESMF_GridComp),intent(inout) :: GC 
    type(date_time_type), intent(in) :: date_time
    
    character(*), intent(in) :: work_path
    
    character(*),  intent(in) :: exp_id

    integer, intent(in) :: N_catl

    type(tile_coord_type),  dimension(:), intent(in) :: tile_coord_l
    
    type(cat_param_type),   dimension(:), intent(in) :: cat_param
    
    type(mwRTM_param_type), dimension(:), intent(in) :: mwRTM_param

    ! ----------------------------
    
    ! local variables
    
    character(300) :: fname
    character( 40) :: file_tag, dir_name
    
    integer :: n
    
    real, dimension(N_catl) :: dztsurf, clsm_wp

    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(MAPL_LocStream) :: locstream
    type(ESMF_Grid)                 :: TILEGRID
    integer, pointer                :: mask(:)
    integer :: rc, status,unit
    character(*),parameter :: Iam="GEOS_output_smapL4SMlmc"


    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_Get(MAPL, LocStream=locstream,rc=status)
    VERIFY_(status)

    call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TileMaskGet(tilegrid,  mask, rc=status)
    VERIFY_(STATUS)

    ! ------------------------------------------------------------------
    !
    ! compute dztsurf

    dztsurf = 0.05   ! now 0.05 m everywhere due to revised CSOIL_2 in subroutine catchment()

    ! convert wilting point from wetness units to volumetric units

    clsm_wp = cat_param%wpwet * cat_param%poros

    ! -------------------

    file_tag = 'ldas_smapL4SMlmc'
    dir_name = 'rc_out'
    
    fname = get_io_filename( work_path, exp_id, file_tag, &
         date_time=date_time, dir_name=dir_name )
    
    unit = GETFILE( trim(fname), form="unformatted", RC=STATUS )
    VERIFY_(STATUS)
    
    if (logit) write (logunit,'(400A)') 'Writing SMAP L4_SM lmc file ' // trim(fname)
    
    ! --------------------
    call MAPL_VarWrite(unit, tilegrid,tile_coord_l(:)%frac_cell,  mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,tile_coord_l(:)%elev     ,  mask=mask, rc=status); VERIFY_(STATUS)
    ! for dzsf, dzrz, and dzpr change units from mm (or kg/m2) to m
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzsf/1000. ,  mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzrz/1000. ,  mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzpr/1000. ,  mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,dztsurf(:) ,               mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzgt(1) ,     mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzgt(2) ,     mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzgt(3) ,     mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzgt(4) ,     mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzgt(5) ,     mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%dzgt(6) ,     mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%poros ,       mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,clsm_wp(:) ,               mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%cdcr1 ,       mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%cdcr2 ,       mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%vegcls,     mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%soilcls ,   mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%sand ,      mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%clay ,      mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%poros ,     mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%wang_wt ,   mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%wang_wp ,   mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_hmin ,  mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_hmax ,  mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_wmin ,  mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_wmax ,  mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_Nrh ,   mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_Nrv ,   mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%rgh_polmix, mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%omega,      mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%bh,         mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%bv,         mask=mask, rc=status); VERIFY_(STATUS)
    call MAPL_VarWrite(unit, tilegrid,mwRTM_param(:)%lewt,       mask=mask, rc=status); VERIFY_(STATUS)                                                      
    call MAPL_VarWrite(unit, tilegrid,cat_param(:)%veghght,       mask=mask, rc=status); VERIFY_(STATUS)
    
    call FREE_FILE(unit, RC=STATUS); VERIFY_(STATUS)
                                          
    if (logit) write (logunit,*) 'done writing'
    if (logit) write (logunit,*) 
    
  end subroutine GEOS_output_smapL4SMlmc

  ! ********************************************************************

end module clsm_ensdrv_out_routines

! *********** EOF ******************************************************
