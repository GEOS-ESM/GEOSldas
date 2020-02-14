#include "MAPL_Generic.h"

module clsm_ensdrv_out_routines
  
  ! collection of LDASsa output subroutines 
  !
  ! (originally in clsm_ensdrv_drv_routines.F90)
  !
  ! reichle, 22 Aug 2014

  use LDAS_ensdrv_globals,           ONLY:     &
       log_master_only,                           &
       logunit,                                   &
       logit,                                     &
       nodata_generic,                            &
       nodata_tol_generic,                        &
       N_bits_shaved
  
  use catch_constants,                  ONLY:     &
       N_snow => CATCH_N_SNOW,                    &
       N_gt   => CATCH_N_GT
  
  use MAPL_ConstantsMod,                ONLY:     &
       alhe   => MAPL_ALHL,                       &
       alhs   => MAPL_ALHS,                       &
       Tzero  => MAPL_TICE

  use LDAS_DriverTypes,                     ONLY:     &
       met_force_type,                            &
       veg_param_type,                            &
       bal_diagn_type,                            &
       out_choice_type,                           &
       out_choice_time_type,                      &
       out_dtstep_type,                           &
       assignment (=),                            &
       operator (+),                              &  
       operator (/)

  use catch_types,                      ONLY:     &
       cat_param_type,                            &
       cat_progn_type,                            &
       cat_diagS_type,                            &
       cat_diagF_type,                            &
       assignment (=),                            &
       operator (+),                              &
       operator (/)

  use LDAS_TileCoordType,                 ONLY:   &
       tile_coord_type,                           &
       grid_def_type,                             &
       io_grid_def_type

  use mwRTM_types,                      ONLY:     &
       mwRTM_param_type,                          &
       io_mwRTM_param_type

  use mwRTM_routines,                   ONLY:     &
       catch2mwRTM_vars,                          &
       mwRTM_get_Tb

  use LDAS_ensdrv_mpi,                  ONLY:     &
       master_proc,                               &
       numprocs

  use LDAS_DateTimeMod,                   ONLY:     &
       date_time_type,                            &
       is_leap_year,                              &
       days_in_month,                             &
       datetime_eq_refdatetime,                   &
       augment_date_time
  
  use clsm_ensdrv_drv_routines,         ONLY:     &
       l2f_real

  use LDAS_ensdrv_init_routines,        ONLY:     &
       clsm_ensdrv_get_command_line,              &
       add_domain_to_path

  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename
       
  use LDAS_TileCoordRoutines,              ONLY:     & 
       tile2grid
  
  use LDAS_ExceptionsMod,                  ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use, intrinsic :: iso_fortran_env, only : output_unit

  implicit none

  private

  public :: init_log             
  public :: output_catparam      
  public :: output_mwRTMparam    
  public :: output_smapL4SMlmc   
  public :: GEOS_output_smapL4SMlmc   
  public :: output_calcs         
  public :: output_write         
  public :: get_ensstd_filenames 
  public :: check_output_times   
  public :: get_land_mask_ij        

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

  subroutine output_catparam( date_time, work_path, exp_id, N_catd, cat_param ) 
    
    ! writes cat_param for domain to file
    !
    ! reichle, 21 Jan 2004
    ! 
    ! -------------------------------------------------------------------
    
    implicit none

    type(date_time_type), intent(in) :: date_time

    character(200), intent(in) :: work_path
    
    character(40),  intent(in) :: exp_id

    integer, intent(in) :: N_catd
    
    type(cat_param_type), dimension(N_catd), intent(in) :: cat_param
    
    ! ----------------------------
    
    ! local variables
    
    character(300) :: fname
    character( 40) :: file_tag, dir_name

    integer :: n, k
    
    ! ------------------------------------------------------------------

    file_tag = 'ldas_catparam'
    dir_name = 'rc_out'

    fname = get_io_filename( work_path, exp_id, file_tag, date_time=date_time, &
         dir_name=dir_name )
    
    open(10, file=fname, form='unformatted', status='unknown', action='write')
    
    if (logit) write (logunit,'(400A)') 'Writing catparam file ' // trim(fname)
    
    write (10) (cat_param(n)%dpth,       n=1,N_catd)
    
    write (10) (cat_param(n)%dzsf,       n=1,N_catd)  
    write (10) (cat_param(n)%dzrz,       n=1,N_catd)  
    write (10) (cat_param(n)%dzpr,       n=1,N_catd)  
    
    do k=1,N_gt
       write (10) (cat_param(n)%dzgt(k), n=1,N_catd)
    end do

    write (10) (cat_param(n)%poros,      n=1,N_catd) 
    write (10) (cat_param(n)%cond,       n=1,N_catd)  
    write (10) (cat_param(n)%psis,       n=1,N_catd)  
    write (10) (cat_param(n)%bee,        n=1,N_catd)   
    
    write (10) (cat_param(n)%wpwet,      n=1,N_catd)  
    
    write (10) (cat_param(n)%gnu,        n=1,N_catd)    
    
    write (10) (cat_param(n)%vgwmax,     n=1,N_catd) 
     
    write (10) (cat_param(n)%vegcls,     n=1,N_catd)  
    write (10) (cat_param(n)%soilcls30,  n=1,N_catd) 
    write (10) (cat_param(n)%soilcls100, n=1,N_catd)  
    
    write (10) (cat_param(n)%bf1,        n=1,N_catd)
    write (10) (cat_param(n)%bf2,        n=1,N_catd)
    write (10) (cat_param(n)%bf3,        n=1,N_catd)
    write (10) (cat_param(n)%cdcr1,      n=1,N_catd)
    write (10) (cat_param(n)%cdcr2,      n=1,N_catd)
    write (10) (cat_param(n)%ars1,       n=1,N_catd)
    write (10) (cat_param(n)%ars2,       n=1,N_catd)
    write (10) (cat_param(n)%ars3,       n=1,N_catd)
    write (10) (cat_param(n)%ara1,       n=1,N_catd)
    write (10) (cat_param(n)%ara2,       n=1,N_catd)
    write (10) (cat_param(n)%ara3,       n=1,N_catd)
    write (10) (cat_param(n)%ara4,       n=1,N_catd)
    write (10) (cat_param(n)%arw1,       n=1,N_catd)
    write (10) (cat_param(n)%arw2,       n=1,N_catd)
    write (10) (cat_param(n)%arw3,       n=1,N_catd)
    write (10) (cat_param(n)%arw4,       n=1,N_catd)
    write (10) (cat_param(n)%tsa1,       n=1,N_catd)
    write (10) (cat_param(n)%tsa2,       n=1,N_catd)
    write (10) (cat_param(n)%tsb1,       n=1,N_catd)
    write (10) (cat_param(n)%tsb2,       n=1,N_catd)
    write (10) (cat_param(n)%atau,       n=1,N_catd)
    write (10) (cat_param(n)%btau,       n=1,N_catd)

    write (10) (cat_param(n)%gravel30,   n=1,N_catd)
    write (10) (cat_param(n)%orgC30  ,   n=1,N_catd)
    write (10) (cat_param(n)%orgC    ,   n=1,N_catd)
    write (10) (cat_param(n)%sand30  ,   n=1,N_catd)
    write (10) (cat_param(n)%clay30  ,   n=1,N_catd)
    write (10) (cat_param(n)%sand    ,   n=1,N_catd)
    write (10) (cat_param(n)%clay    ,   n=1,N_catd)
    write (10) (cat_param(n)%wpwet30 ,   n=1,N_catd)
    write (10) (cat_param(n)%poros30 ,   n=1,N_catd)

    write (10) (cat_param(n)%veghght ,   n=1,N_catd)

    close (10,status='keep')

    if (logit) write (logunit,*) 'done writing'
    if (logit) write (logunit,*) 
    
  end subroutine output_catparam

  ! ********************************************************************

  subroutine output_mwRTMparam( date_time, work_path, exp_id, N_catd, mwRTM_param ) 
    
    ! writes mwRTM_param for domain to file
    !
    ! reichle, 1 Jun 2011
    ! 
    ! -------------------------------------------------------------------
    
    implicit none

    type(date_time_type),                      intent(in   ) :: date_time

    character(200),                            intent(in   ) :: work_path
    
    character(40),                             intent(in   ) :: exp_id

    integer,                                   intent(in   ) :: N_catd
    
    type(mwRTM_param_type), dimension(N_catd), intent(inout) :: mwRTM_param
    
    ! ----------------------------
    
    ! local variables
    
    character(300) :: fname
    character( 40) :: file_tag, dir_name

    ! ------------------------------------------------------------------

    file_tag = 'ldas_mwRTMparam'
    dir_name = 'rc_out'

    fname = get_io_filename( work_path, exp_id, file_tag, date_time=date_time, &
         dir_name=dir_name )
    
    open(10, file=fname, form='unformatted', status='unknown', action='write')
    
    if (logit) write (logunit,'(400A)') 'Writing mwRTMparam file ' // trim(fname)
    
    call io_mwRTM_param_type( 'w', 10, N_catd, mwRTM_param )
    
    close (10,status='keep')
    
    if (logit) write (logunit,*) 'done writing'
    if (logit) write (logunit,*) 
    
  end subroutine output_mwRTMparam

  ! ********************************************************************

  subroutine output_smapL4SMlmc( date_time, work_path, exp_id, &    
       N_catf, tile_coord, cat_param, mwRTM_param ) 
    
    ! write SMAP L4_SM "lmc" (land model constants) file collection
    ! as binary, tile-space output
    !
    ! requires "full" domain inputs
    !
    ! reichle, 26 Apr 2012
    ! reichle, 27 May 2014: - changed wilting point output from "clsm_wpwet" to "clsm_wp"
    ! 
    ! -------------------------------------------------------------------
    
    implicit none

    type(date_time_type), intent(in) :: date_time
    
    character(200), intent(in) :: work_path
    
    character(40),  intent(in) :: exp_id

    integer, intent(in) :: N_catf

    type(tile_coord_type),  dimension(N_catf), intent(in) :: tile_coord
    
    type(cat_param_type),   dimension(N_catf), intent(in) :: cat_param
    
    type(mwRTM_param_type), dimension(N_catf), intent(in) :: mwRTM_param

    ! ----------------------------
    
    ! local variables
    
    character(300) :: fname
    character( 40) :: file_tag, dir_name
    
    integer :: n
    
    real, dimension(N_catf) :: dztsurf, clsm_wp
    
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
    
    open(10, file=fname, form='unformatted', status='unknown', action='write')
    
    if (logit) write (logunit,'(400A)') 'Writing SMAP L4_SM lmc file ' // trim(fname)
    
    ! --------------------

    write (10) (tile_coord(n)%frac_cell,   n=1,N_catf)     !  1: real
    write (10) (tile_coord(n)%elev,        n=1,N_catf)     !  2: real

    ! for dzsf, dzrz, and dzpr change units from mm (or kg/m2) to m

    write (10) (cat_param(n)%dzsf/1000.,   n=1,N_catf)     !  3: real 
    write (10) (cat_param(n)%dzrz/1000.,   n=1,N_catf)     !  4: real 
    write (10) (cat_param(n)%dzpr/1000.,   n=1,N_catf)     !  5: real 
                                                      
    write (10) (dztsurf(n),                n=1,N_catf)     !  6: real 
                                                      
    write (10) (cat_param(n)%dzgt(1),      n=1,N_catf)     !  7: real 
    write (10) (cat_param(n)%dzgt(2),      n=1,N_catf)     !  8: real 
    write (10) (cat_param(n)%dzgt(3),      n=1,N_catf)     !  9: real 
    write (10) (cat_param(n)%dzgt(4),      n=1,N_catf)     ! 10: real 
    write (10) (cat_param(n)%dzgt(5),      n=1,N_catf)     ! 11: real 
    write (10) (cat_param(n)%dzgt(6),      n=1,N_catf)     ! 12: real 
                                                      
    write (10) (cat_param(n)%poros,        n=1,N_catf)     ! 13: real 
    write (10) (clsm_wp(n),                n=1,N_catf)     ! 14: real 
                                                      
    write (10) (cat_param(n)%cdcr1,        n=1,N_catf)     ! 15: real 
    write (10) (cat_param(n)%cdcr2,        n=1,N_catf)     ! 16: real 

    write (10) (mwRTM_param(n)%vegcls,     n=1,N_catf)     ! 17: integer !!!
    write (10) (mwRTM_param(n)%soilcls,    n=1,N_catf)     ! 18: integer !!!
                                              
    write (10) (mwRTM_param(n)%sand,       n=1,N_catf)     ! 19: real 
    write (10) (mwRTM_param(n)%clay,       n=1,N_catf)     ! 20: real 
                                              
    write (10) (mwRTM_param(n)%poros,      n=1,N_catf)     ! 21: real 
                                              
    write (10) (mwRTM_param(n)%wang_wt,    n=1,N_catf)     ! 22: real 
    write (10) (mwRTM_param(n)%wang_wp,    n=1,N_catf)     ! 23: real 
                                              
    write (10) (mwRTM_param(n)%rgh_hmin,   n=1,N_catf)     ! 24: real 
    write (10) (mwRTM_param(n)%rgh_hmax,   n=1,N_catf)     ! 25: real 
    write (10) (mwRTM_param(n)%rgh_wmin,   n=1,N_catf)     ! 26: real 
    write (10) (mwRTM_param(n)%rgh_wmax,   n=1,N_catf)     ! 27: real 
    write (10) (mwRTM_param(n)%rgh_Nrh,    n=1,N_catf)     ! 28: real 
    write (10) (mwRTM_param(n)%rgh_Nrv,    n=1,N_catf)     ! 29: real 
    write (10) (mwRTM_param(n)%rgh_polmix, n=1,N_catf)     ! 30: real 
                                              
    write (10) (mwRTM_param(n)%omega,      n=1,N_catf)     ! 31: real 
                                              
    write (10) (mwRTM_param(n)%bh,         n=1,N_catf)     ! 32: real 
    write (10) (mwRTM_param(n)%bv,         n=1,N_catf)     ! 33: real 
    write (10) (mwRTM_param(n)%lewt,       n=1,N_catf)     ! 34: real 

    write (10) (cat_param(n)%veghght,      n=1,N_catf)     ! 35: real
                                                      
    close (10,status='keep')
    
    if (logit) write (logunit,*) 'done writing'
    if (logit) write (logunit,*) 
    
  end subroutine output_smapL4SMlmc

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

  subroutine output_calcs(                                       &
       option, out_collection_ID, Nt,                            &
       cat_progn,     cat_diagS,     cat_diagF,                  &
       met_force,     veg_param,     bal_diagn,                  &
       cat_progn_avg, cat_diagS_avg, cat_diagF_avg,              &
       met_force_avg, veg_param_avg, bal_diagn_avg,              &
       cat_param, mwRTM_param, tile_coord,                       &
       out_choice, date_time, work_path, exp_id, interval,       &
       model_dtstep,                                             &
       out_tile, out_grid, fname_tile, fname_grid, tile_data,    &
       out_dtstep_xhourly, ens_id )
    
    ! calculations and preparation for inst, xhourly, daily, pentad, 
    ! and monthly output
    !
    ! option = 'ini' : initialize
    ! option = 'add' : add values of current time step into "_avg" 
    ! option = 'out' : normalize and output avg
    !
    ! reichle, 29 Sep 2009 - reformulated based on old subroutine calc_tavg()
    !                        see also new subroutine write_output()
    ! reichle,  9 Dec 2011 - revised using new types "veg_param" and "bal_diagn"
    ! reichle, 28 Dec 2011 - removed field "totalb" from "cat_diagn" structure
    ! reichle, 31 Oct 2013 - split "cat_diagn" structure into "cat_diagS" and "cat_diagF"
    !
    ! ---------------------------------------------------------------------
    
    implicit none
    
    integer,      intent(in) :: out_collection_ID

    integer,      intent(in) :: Nt    ! previously named N_catd
        
    character(3), intent(in) :: option
    
    type(cat_progn_type), dimension(Nt), intent(in)    :: cat_progn
    type(cat_diagS_type), dimension(Nt), intent(in)    :: cat_diagS
    type(cat_diagF_type), dimension(Nt), intent(in)    :: cat_diagF
    type(met_force_type), dimension(Nt), intent(in)    :: met_force    
    type(veg_param_type), dimension(Nt), intent(in)    :: veg_param
    type(bal_diagn_type), dimension(Nt), intent(in)    :: bal_diagn
    
    type(cat_progn_type), dimension(Nt), intent(inout) :: cat_progn_avg
    type(cat_diagS_type), dimension(Nt), intent(inout) :: cat_diagS_avg
    type(cat_diagF_type), dimension(Nt), intent(inout) :: cat_diagF_avg
    type(met_force_type), dimension(Nt), intent(inout) :: met_force_avg
    type(veg_param_type), dimension(Nt), intent(inout) :: veg_param_avg
    type(bal_diagn_type), dimension(Nt), intent(inout) :: bal_diagn_avg

    ! optional inputs (needed when writing average to file)

    type(cat_param_type),   dimension(Nt), intent(in), optional :: cat_param
    
    type(mwRTM_param_type), dimension(Nt), intent(in), optional :: mwRTM_param

    type(tile_coord_type),  dimension(Nt), intent(in), optional :: tile_coord
    
    type(out_choice_type),                 intent(in), optional :: out_choice
    
    type(date_time_type),                  intent(in), optional :: date_time
    
    character(200),                        intent(in), optional :: work_path
    character(40),                         intent(in), optional :: exp_id
    
    ! what averaging interval is used? (need to know to construct file name)
    !
    ! inst   : interval = 'i'  -- for inst output ONLY call with option 'out'
    ! xhourly: interval = 'x'
    ! daily  : interval = 'd'
    ! pentad : interval = 'p'
    ! monthly: interval = 'm'
    
    character, intent(in), optional :: interval
    
    logical, intent(out), optional :: out_tile, out_grid
    
    character(300), intent(out), optional :: fname_tile, fname_grid

    ! changed tile_data to pointer so that compile with "-check bounds" works
    ! - reichle, 8 Feb 2013

    real, dimension(:,:), pointer, optional :: tile_data  ! intent(out)
                                                          ! dimension(Nt,N_out_fields) 
    
    integer, intent(in), optional :: model_dtstep, out_dtstep_xhourly, ens_id
        
    ! ----------------------------------------
    
    ! local variables
    
    integer              :: n, k, ens_id_tmp
    
    real                 :: n_steps, totalb, ar4, freq, inc_angle
    
    real, parameter      :: daylen = 86400.
    
    character(40)        :: file_tag
    
    logical              :: out_wetness, muststop, incl_atm_terms
    
    type(date_time_type) :: date_time_tmp

    real, dimension(Nt)  :: tmpreal, SWE, sfmc_mwRTM, tsoil_mwRTM, Tb_h, Tb_v
    
    type(cat_diagS_type) :: cat_diagS_tmp

    character(len=*), parameter :: Iam = 'output_calcs'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------

    out_wetness = .false.

    ! Removed "out_wetness" from LDASsa nml inputs and hard-wired 
    ! to "false".
    ! If needed for backward compatibility, add replicates of 
    ! Collections 1, 2, 3, or 10, assign new Collection IDs, 
    ! and set out_wetness=.true. below.
    !
    ! - reichle, 27 Aug 2014
    
    !select case (out_collection_ID)
    !   
    !   case ( "list of new_collection_IDs" )  
    !   
    !case default
    !   
    !   out_wetness = .false.
    !   
    !end select
    
    ! ------------------------------------------------------------

    if (option=='ini') then    ! initialize
       
       do n=1,Nt
          
          cat_progn_avg(n) = 0.
          cat_diagS_avg(n) = 0.
          cat_diagF_avg(n) = 0.
          met_force_avg(n) = 0.
          veg_param_avg(n) = 0.
          bal_diagn_avg(n) = 0.
          
       end do
       
    else if (option=='add') then   ! sum up for average
       
       do n=1,Nt

          cat_progn_avg(n) = cat_progn_avg(n) + cat_progn(n)

          ! ----------------
          
          ! In catchment() the snow temperature is set to TSURF if ASNOW=0.
          ! Exclude snow-free temperatures from longer-term (eg, monthly)
          ! time averages by weighting snow temperature with snow cover
          ! fraction.
          ! - reichle, 29 Feb 2012
          
          cat_diagS_tmp = cat_diagS(n)
          
          cat_diagS_tmp%tpsn(1:N_snow) = cat_diagS(n)%tpsn(1:N_snow)*cat_diagS(n)%asnow
          
          cat_diagS_avg(n) = cat_diagS_avg(n) + cat_diagS_tmp

          ! ----------------

          cat_diagF_avg(n) = cat_diagF_avg(n) + cat_diagF(n)
          
          met_force_avg(n) = met_force_avg(n) + met_force(n)
          
          veg_param_avg(n) = veg_param_avg(n) + veg_param(n)

          bal_diagn_avg(n) = bal_diagn_avg(n) + bal_diagn(n)
          
       end do
       
    else if (option=='out') then   ! finalize, output, re-initialize averages
       
       ! prepare for output 
       
       ! make sure all optional arguments are present that are required for 'out'
       
       muststop = .false.
       
       if (.not. present(out_choice))   muststop=.true.       
       if (.not. present(date_time))    muststop=.true.
       if (.not. present(work_path))    muststop=.true.
       if (.not. present(exp_id))       muststop=.true.
       if (.not. present(interval))     muststop=.true.
       if (.not. present(model_dtstep)) muststop=.true.
       
       if (muststop) then
          err_msg = 'optional input arguments missing'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       if (.not. present(ens_id)) then
          
          ens_id_tmp = -1   ! default to ensemble average 
          
       else
          
          ens_id_tmp = ens_id
          
       end if
       
       
       ! put together normalization factors and file names, also determine 
       ! whether tile or grid output is desired
       
       out_tile = .false.
       out_grid = .false.
       
       select case (interval)
          
       case ('i')     ! instantaneous
          
          n_steps = 1.

          if (out_choice%inst%tile) then
             
             out_tile = .true.
             
             file_tag = 'ldas_tile_inst_out'
             
             fname_tile = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time, ens_id=ens_id_tmp )
             
          end if
          
          if (out_choice%inst%grid) then
             
             out_grid = .true.
             
             file_tag = 'ldas_grid_inst_out'
             
             fname_grid = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time, ens_id=ens_id_tmp )
             
          end if
          
       case ('x')     ! xhourly
          
          if (present(out_dtstep_xhourly)) then
             
             n_steps = real(out_dtstep_xhourly/model_dtstep)
             
          else

             err_msg = 'need out_dtstep_xhourly'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             
          end if
             
          ! shift date/time so that mid-point of time averaging interval is
          ! used for time-tagging the output file
          
          date_time_tmp = date_time
          
          call augment_date_time( -out_dtstep_xhourly/2, date_time_tmp )
          
          if (out_choice%xhourly%tile) then
             
             out_tile = .true.
             
             file_tag = 'ldas_tile_xhourly_out'
             
             fname_tile = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time_tmp, ens_id=ens_id_tmp )
             
          end if
          
          if (out_choice%xhourly%grid) then
             
             out_grid = .true.
             
             file_tag = 'ldas_grid_xhourly_out'
             
             fname_grid = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time_tmp, ens_id=ens_id_tmp )
             
          end if
          
       case ('d')     ! daily
          
          n_steps = real(86400/model_dtstep)
          
          if (out_choice%daily%tile) then
             
             out_tile = .true.
             
             file_tag = 'ldas_tile_daily_out'  
             
             fname_tile = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time, option=4, ens_id=ens_id_tmp )
             
          end if
          
          if (out_choice%daily%grid) then
             
             out_grid = .true.
             
             file_tag = 'ldas_grid_daily_out'
             
             fname_grid = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time, option=4, ens_id=ens_id_tmp )
             
          end if

       case ('p')     ! pentad
          
          if (date_time%pentad==12 .and. is_leap_year(date_time%year)) then
             
             n_steps = real(6*86400/model_dtstep)
             
          else
             n_steps = real(5*86400/model_dtstep)
             
          end if
          
          if (out_choice%pentad%tile) then
             
             out_tile = .true.
             
             file_tag = 'ldas_tile_pentad_out'
             
             fname_tile = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time, option=3, ens_id=ens_id_tmp )
             
          end if

          if (out_choice%pentad%grid) then
             
             out_grid = .true.
             
             file_tag = 'ldas_grid_pentad_out'
             
             fname_grid = get_io_filename( work_path, exp_id, file_tag,  &
                  date_time=date_time, option=3, ens_id=ens_id_tmp )
             
          end if
          
       case ('m')     ! monthly
          
          n_steps = real(days_in_month(date_time%year,date_time%month)*86400/model_dtstep)
          
          if (out_choice%monthly%tile) then
             
             out_tile = .true.
             
             file_tag = 'ldas_tile_monthly_out'
             
             fname_tile = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time, option=2, ens_id=ens_id_tmp )
   
          end if

          if (out_choice%monthly%grid) then
             
             out_grid = .true.
             
             file_tag = 'ldas_grid_monthly_out'
             
             fname_grid = get_io_filename( work_path, exp_id, file_tag, &
                  date_time=date_time, option=2, ens_id=ens_id_tmp )
   
          end if
          
       case default
          
          err_msg = 'unknown averaging interval'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end select   ! case (interval)

       ! -----------------------
       !
       ! normalize, fix no-data-values, change units
       
       do n=1,Nt
          
          ! normalize by number of time steps
          
          cat_progn_avg(n) = cat_progn_avg(n) / n_steps
          cat_diagS_avg(n) = cat_diagS_avg(n) / n_steps
          cat_diagF_avg(n) = cat_diagF_avg(n) / n_steps
          met_force_avg(n) = met_force_avg(n) / n_steps
          veg_param_avg(n) = veg_param_avg(n) / n_steps
          bal_diagn_avg(n) = bal_diagn_avg(n) / n_steps
          
          ! In catchment() the snow temperature is set to TSURF if ASNOW=0.
          ! Exclude snow-free temperatures from longer-term (eg, monthly)
          ! time averages by weighting snow temperature with snow cover
          ! fraction.
          ! - reichle, 29 Feb 2012
          
          if (cat_diagS_avg(n)%asnow>0.) then
             
             ! normalize asnow-weighted time-average 
             ! (except for instantaneous ('i') output)
             ! -reichle+csdraper, 31 Oct 2013
             
             if (interval/='i')                                            &
                  cat_diagS_avg(n)%tpsn(1:N_snow) =                        &
                  cat_diagS_avg(n)%tpsn(1:N_snow)/cat_diagS_avg(n)%asnow
             
          else
             
             cat_diagS_avg(n)%tpsn(1:N_snow) = nodata_generic             
             
          end if
          
          ! change sub-tile canopy temperatures and spec humidities to 
          ! no-data-values when corresponding area fraction is zero
          ! reichle, 29 Feb 2012
          
          ar4 = 1. - cat_diagS_avg(n)%ar1 - cat_diagS_avg(n)%ar2
          
          if (cat_diagS_avg(n)%ar1<=0.)  cat_progn_avg(n)%tc1 = nodata_generic
          if (cat_diagS_avg(n)%ar2<=0.)  cat_progn_avg(n)%tc2 = nodata_generic
          if (                 ar4<=0.)  cat_progn_avg(n)%tc4 = nodata_generic
          
          if (cat_diagS_avg(n)%ar1<=0.)  cat_progn_avg(n)%qa1 = nodata_generic
          if (cat_diagS_avg(n)%ar2<=0.)  cat_progn_avg(n)%qa2 = nodata_generic
          if (                 ar4<=0.)  cat_progn_avg(n)%qa4 = nodata_generic
          
          
          ! change units for selected outputs 
          !
          ! reichle, 29 Feb 2012: MUST add no-data-check if changing units of 
          !                       tpsn, tc1, tc2, tc4, qa1, qa2, qa4 
          
          select case (out_collection_ID)

          case (1,2,3,10)
             
             ! convert [kg/m2/s] into [mm/day]
             
             met_force_avg(n)%Rainf_C = met_force_avg(n)%Rainf_C * daylen
             met_force_avg(n)%Rainf   = met_force_avg(n)%Rainf   * daylen
             met_force_avg(n)%Snowf   = met_force_avg(n)%Snowf   * daylen
             
             cat_diagF_avg(n)%evap   = cat_diagF_avg(n)%evap   * daylen
             
             cat_diagF_avg(n)%runoff = cat_diagF_avg(n)%runoff * daylen
             cat_diagF_avg(n)%runsrf = cat_diagF_avg(n)%runsrf * daylen
             cat_diagF_avg(n)%bflow  = cat_diagF_avg(n)%bflow  * daylen
             
             cat_diagF_avg(n)%snmelt = cat_diagF_avg(n)%snmelt * daylen
             
             bal_diagn_avg(n)%wchng  = bal_diagn_avg(n)%wchng  * daylen
             bal_diagn_avg(n)%wincr  = bal_diagn_avg(n)%wincr  * daylen
             
             ! convert [W/m2] into [mm/day] 
             
             cat_diagF_avg(n)%eint   = cat_diagF_avg(n)%eint * daylen/alhe
             cat_diagF_avg(n)%eveg   = cat_diagF_avg(n)%eveg * daylen/alhe
             cat_diagF_avg(n)%esno   = cat_diagF_avg(n)%esno * daylen/alhs
             cat_diagF_avg(n)%esoi   = cat_diagF_avg(n)%esoi * daylen/alhe
             
             ! units of qinfil should probably be changed
             
             ! why are units of energy_bal not changed?  
             !  perhaps b/c division by seconds works out just fine (Joule -> Watts )??

          case (4,5,6,7,8,9,11)
             
             ! no unit changes

          case default
             
             err_msg = 'unknown out_collection_ID'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
             
          end select
          
          ! ALWAYS convert units of tp from deg C to Kelvin
          
          do k=1,N_gt
             
             cat_diagS_avg(n)%tp(k) = cat_diagS_avg(n)%tp(k) + Tzero
             
          end do
          
       end do
       
       ! ----------------------------------------------------------
       !
       ! assemble data that are actually output
       
       select case (out_collection_ID)
          
       case (1,10)    !  1 - N_out_fields_inst = N_out_fields_tavg = 44
                      ! 10 - N_out_fields_inst = N_out_fields_tavg = 46

          !  1 - Legacy LDASsa off-line output specs ("44 variables")
          ! 10 - Legacy plus t2m, q2m 
          
          do n=1,Nt
             
             if (met_force_avg(n)%SWdown > 1e-4) then
                totalb = cat_diagF_avg(n)%swup/met_force_avg(n)%SWdown
             else
                totalb = nodata_generic
             end if
             
             tile_data(n, 1) = met_force_avg(n)%Tair                
             tile_data(n, 2) = met_force_avg(n)%Qair                
             tile_data(n, 3) = met_force_avg(n)%Psurf               
             tile_data(n, 4) = met_force_avg(n)%Rainf_C             
             tile_data(n, 5) = met_force_avg(n)%Rainf               
             tile_data(n, 6) = met_force_avg(n)%Snowf               
             tile_data(n, 7) = met_force_avg(n)%LWdown              
             tile_data(n, 8) = met_force_avg(n)%SWdown              
             tile_data(n, 9) = met_force_avg(n)%Wind                
             
             tile_data(n,10) = cat_progn_avg(n)%capac               
             tile_data(n,11) = cat_progn_avg(n)%srfexc              
             tile_data(n,12) = cat_progn_avg(n)%rzexc               
             tile_data(n,13) = cat_progn_avg(n)%catdef              
             tile_data(n,14) = sum(cat_progn_avg(n)%wesn(1:N_snow)) 
             tile_data(n,15) = sum(cat_progn_avg(n)%sndz(1:N_snow)) 
             
             tile_data(n,16) = cat_diagS_avg(n)%ar1                 
             tile_data(n,17) = cat_diagS_avg(n)%ar2                 
             tile_data(n,18) = cat_diagS_avg(n)%asnow               
             
             if (out_wetness) then
                
                tile_data(n,19) = max(min(cat_diagS_avg(n)%sfmc/cat_param(n)%poros,1.),0.)
                tile_data(n,20) = max(min(cat_diagS_avg(n)%rzmc/cat_param(n)%poros,1.),0.)
                tile_data(n,21) = max(min(cat_diagS_avg(n)%prmc/cat_param(n)%poros,1.),0.)
                
             else
                
                tile_data(n,19) = cat_diagS_avg(n)%sfmc                
                tile_data(n,20) = cat_diagS_avg(n)%rzmc                
                tile_data(n,21) = cat_diagS_avg(n)%prmc                
                
             end if
             
             tile_data(n,22) = cat_diagS_avg(n)%tsurf 
             tile_data(n,23) = cat_diagS_avg(n)%tp(1)               
             tile_data(n,24) = cat_diagS_avg(n)%tp(N_gt)        
             tile_data(n,25) = cat_diagS_avg(n)%tpsn(1)             
             tile_data(n,26) = cat_diagS_avg(n)%tpsn(N_snow)        

             tile_data(n,27) = cat_diagF_avg(n)%shflux              
             tile_data(n,28) = cat_diagF_avg(n)%lhflux              
             tile_data(n,29) = cat_diagF_avg(n)%ghflux              
             tile_data(n,30) = cat_diagF_avg(n)%evap                
             tile_data(n,31) = cat_diagF_avg(n)%eint                
             tile_data(n,32) = cat_diagF_avg(n)%eveg                
             tile_data(n,33) = cat_diagF_avg(n)%esoi                
             tile_data(n,34) = cat_diagF_avg(n)%esno                
             tile_data(n,35) = cat_diagF_avg(n)%runoff              
             tile_data(n,36) = cat_diagF_avg(n)%runsrf              
             tile_data(n,37) = cat_diagF_avg(n)%bflow               
             tile_data(n,38) = cat_diagF_avg(n)%snmelt              
             tile_data(n,39) = cat_diagF_avg(n)%lwup                
             tile_data(n,40) = cat_diagF_avg(n)%swup                
             tile_data(n,41) = cat_diagF_avg(n)%qinfil              

             tile_data(n,42) = totalb              

             tile_data(n,43) = bal_diagn_avg(n)%wincr
             tile_data(n,44) = bal_diagn_avg(n)%eincr 

             if (out_collection_ID==10) then 
                
                tile_data(n,45) = cat_diagF_avg(n)%t2m 
                tile_data(n,46) = cat_diagF_avg(n)%q2m 

             end if
                          
          end do

       case (2)    ! N_out_fields_inst = N_out_fields_tavg = 6
          
          ! Specs for SMAP Nature Run v02 (Feb 2011)   

          do n=1,Nt
             
             if (out_wetness) then
                
                tile_data(n,1) = max(min(cat_diagS_avg(n)%sfmc/cat_param(n)%poros,1.),0.)
                tile_data(n,2) = max(min(cat_diagS_avg(n)%rzmc/cat_param(n)%poros,1.),0.)
                
             else
                
                tile_data(n,1) = cat_diagS_avg(n)%sfmc                
                tile_data(n,2) = cat_diagS_avg(n)%rzmc                
                
             end if
             
             tile_data(n,3) = cat_diagS_avg(n)%tsurf               
             tile_data(n,4) = cat_diagS_avg(n)%tp(1)               
             tile_data(n,5) = sum(cat_progn_avg(n)%wesn(1:N_snow)) 
             tile_data(n,6) = met_force_avg(n)%Rainf+met_force_avg(n)%Snowf
             
          end do

       case (3)    ! N_out_fields_inst = N_out_fields_tavg = 8

          ! for L-band mwRTM calibration (before Dec 2013), SMOS DA
          
          do n=1,Nt

             if (out_wetness) then

                tile_data(n,1) = max(min(cat_diagS_avg(n)%sfmc/cat_param(n)%poros,1.),0.)
                tile_data(n,2) = max(min(cat_diagS_avg(n)%rzmc/cat_param(n)%poros,1.),0.)

             else

                tile_data(n,1) = cat_diagS_avg(n)%sfmc
                tile_data(n,2) = cat_diagS_avg(n)%rzmc

             end if

             tile_data(n,3) = cat_diagS_avg(n)%tsurf
             tile_data(n,4) = cat_diagS_avg(n)%tp(1)
             tile_data(n,5) = sum(cat_progn_avg(n)%wesn(1:N_snow))
             tile_data(n,6) = met_force_avg(n)%Rainf+met_force_avg(n)%Snowf
             tile_data(n,7) = met_force_avg(n)%Tair
             tile_data(n,8) = cat_progn_avg(n)%capac

          end do


       case (4,5)    ! N_out_fields_inst = N_out_fields_tavg = 50, 59
          
          ! 50 MERRA-Land "mld" outputs  *or*  59 = 50 MERRA-Land "mld" outputs plus 9 additional fields
          !
          ! reichle, 29 Feb 2012: updated to reflect final "mld" file specs (incl TSURF)
                       
          tile_data(1:Nt, 1) = veg_param_avg%grn                                  ! GRN 	Fraction 
          tile_data(1:Nt, 2) = veg_param_avg%lai                                  ! LAI 	m2 m-2 
                                                                                 
          tile_data(1:Nt, 3) = max(min(cat_diagS_avg%prmc/cat_param%poros,1.),0.) ! GWETPROF 	Fraction
          tile_data(1:Nt, 4) = max(min(cat_diagS_avg%rzmc/cat_param%poros,1.),0.) ! GWETROOT 	Fraction
          tile_data(1:Nt, 5) = max(min(cat_diagS_avg%sfmc/cat_param%poros,1.),0.) ! GWETTOP 	Fraction
                                                                                 
          tile_data(1:Nt, 6) = cat_diagS_avg%prmc                                 ! PRMC	m3/m3
          tile_data(1:Nt, 7) = cat_diagS_avg%rzmc                                 ! RZMC	m3/m3
          tile_data(1:Nt, 8) = cat_diagS_avg%sfmc                                 ! SFMC	m3/m3

          tile_data(1:Nt, 9) = cat_diagS_avg%tsurf                                ! TSURF 	K
          tile_data(1:Nt,10) = cat_diagS_avg%tpsn(1)                              ! TPSNOW 	K 

          tile_data(1:Nt,11) = cat_progn_avg%tc2                                  ! TUNST 	K 
          tile_data(1:Nt,12) = cat_progn_avg%tc1                                  ! TSAT 	K 
          tile_data(1:Nt,13) = cat_progn_avg%tc4                                  ! TWLT 	K 

          tile_data(1:Nt,14) = cat_diagS_avg%tp(1)                                ! TSOIL1	K 
          tile_data(1:Nt,15) = cat_diagS_avg%tp(2)                                ! TSOIL2	K 
          tile_data(1:Nt,16) = cat_diagS_avg%tp(3)                                ! TSOIL3	K 
          tile_data(1:Nt,17) = cat_diagS_avg%tp(4)                                ! TSOIL4	K 
          tile_data(1:Nt,18) = cat_diagS_avg%tp(5)                                ! TSOIL5	K 
          tile_data(1:Nt,19) = cat_diagS_avg%tp(6)                                ! TSOIL6	K 

          tile_data(1:Nt,20) = met_force_avg%Snowf                                ! PRECSNO  	kg m-2 s-1
          tile_data(1:Nt,21) = met_force_avg%Rainf + met_force_avg%Snowf          ! PRECTOT 	kg m-2 s-1

          do n=1,Nt
             
             tile_data(n,22) = sum(cat_progn_avg(n)%wesn(1:N_snow))               ! SNOMAS	kg m-2
             tile_data(n,23) = sum(cat_progn_avg(n)%sndz(1:N_snow))               ! SNODP	m 
             
          end do
          
          tile_data(1:Nt,24) = cat_diagF_avg%esoi                                 ! EVPSOIL  	W m-2
          tile_data(1:Nt,25) = cat_diagF_avg%eveg                                 ! EVPTRNS  	W m-2
          tile_data(1:Nt,26) = cat_diagF_avg%eint                                 ! EVPINTR 	W m-2
          tile_data(1:Nt,27) = cat_diagF_avg%esno                                 ! EVPSBLN 	W m-2

          tile_data(1:Nt,28) = cat_diagF_avg%runsrf                               ! RUNOFF 	kg m-2 s-1
          tile_data(1:Nt,29) = cat_diagF_avg%bflow                                ! BASEFLOW	kg m-2 s-1
          tile_data(1:Nt,30) = cat_diagF_avg%snmelt                               ! SMLAND 	kg m-2 s-1
          tile_data(1:Nt,31) = cat_diagF_avg%qinfil                               ! QINFIL	kg m-2 s-1

          ! Note: ar1+ar2+ar4=1 but need FRSAT+FRUNST+FRWLT+FRSNO=1

          tmpreal(1:Nt) = max(min((1.-cat_diagS_avg%asnow),1.),0.)  ! precompute snow-free fraction
          
          tile_data(1:Nt,32) = max(min(tmpreal*cat_diagS_avg%ar2,  1.),0.)        ! FRUNST 	Fraction
          tile_data(1:Nt,33) = max(min(tmpreal*cat_diagS_avg%ar1,  1.),0.)        ! FRSAT 	Fraction
          tile_data(1:Nt,34) = max(min(        cat_diagS_avg%asnow,1.),0.)        ! FRSNO 	Fraction
          
          tmpreal = tmpreal-tile_data(1:Nt,32)-tile_data(1:Nt,33)   ! compute wilting fraction   
          
          ! tmpreal = 1.-sum(tile_data(1:Nt,31:33),dim=2)
          ! tmpreal = tmpreal * (1.-cat_diagS_avg%ar1-cat_diagS_avg%ar2) 
          
          tile_data(1:Nt,35) = max(min(tmpreal,                    1.),0.)        ! FRWLT 	fraction
          

          tile_data(1:Nt,36) = met_force_avg%PARdffs                              ! PARDF 	W m-2
          tile_data(1:Nt,37) = met_force_avg%PARdrct                              ! PARDR 	W m-2

          tile_data(1:Nt,38) = cat_diagF_avg%shflux                               ! SHLAND 	W m-2
          tile_data(1:Nt,39) = cat_diagF_avg%lhflux                               ! LHLAND 	W m-2
          tile_data(1:Nt,40) = cat_diagF_avg%evap                                 ! EVLAND 	kg m-2 s-1

          tile_data(1:Nt,41) = met_force_avg%LWdown - cat_diagF_avg%lwup          ! LWLAND 	W m-2
          tile_data(1:Nt,42) = met_force_avg%SWdown - cat_diagF_avg%swup          ! SWLAND 	W m-2

          tile_data(1:Nt,43) = cat_diagF_avg%ghflux                               ! GHLAND 	W m-2

          tile_data(1:Nt,44) = bal_diagn_avg%wtotl                                ! TWLAND 	kg m-2
          tile_data(1:Nt,45) = bal_diagn_avg%etotl                                ! TELAND 	J m-2 
          tile_data(1:Nt,46) = bal_diagn_avg%wchng                                ! WCHANGE  	kg m-2 s-1
          tile_data(1:Nt,47) = bal_diagn_avg%echng                                ! ECHANGE 	W m-2  

          tile_data(1:Nt,48) = 0.                                                 ! SPLAND	W m-2  
          tile_data(1:Nt,49) = 0.                                                 ! SPWATR	kg m-2 s-1
          tile_data(1:Nt,50) = cat_diagF_avg%hsnacc                               ! SPSNOW	W m-2  
          
          if (out_collection_ID==5) then
             
             ! select additional outputs 

             tile_data(1:Nt,51) = met_force_avg%Tair                              ! TLML  	K
             tile_data(1:Nt,52) = met_force_avg%Qair                              ! QLML  	kg kg-1
             tile_data(1:Nt,53) = met_force_avg%LWdown                            ! LWGAB 	W m-2
             tile_data(1:Nt,54) = met_force_avg%SWdown                            ! SWGDN 	W m-2
                                                                                 
             tile_data(1:Nt,55) = cat_progn_avg%srfexc                            !       	kg m-2 
             tile_data(1:Nt,56) = cat_progn_avg%rzexc                             !       	kg m-2  
             tile_data(1:Nt,57) = cat_progn_avg%catdef                            !       	kg m-2    
                                                                                 
             tile_data(1:Nt,58) = bal_diagn_avg%wincr                             !       	kg m-2
             tile_data(1:Nt,59) = bal_diagn_avg%eincr                             !       	J m-2       
                          
          end if


       case (6)      ! N_out_fields_inst = N_out_fields_tavg = 40   (EXCL sm in pctl units!)
          
          ! output fields of SMAP L4_SM gph collection, 
          ! order follows Table 9 of SMAP L4_SM Data Products Specification Document (PSD; revised Jun 2014) 
          !
          ! NOTE: rootzone and profile soil moisture outputs in units of percentiles are appended in post-processing
          !
          ! - reichle,  9 Apr 2013
          ! - reichle, 26 May 2014 - revised output specs: replaced sm in pctl units (fields 1-3) with sm in volumetric units
             
          tile_data(1:Nt, 1) = cat_diagS_avg%sfmc                                 ! sm_surface                        m3 m-3
          tile_data(1:Nt, 2) = cat_diagS_avg%rzmc                                 ! sm_rootzone                       m3 m-3
          tile_data(1:Nt, 3) = cat_diagS_avg%prmc                                 ! sm_profile                        m3 m-3
          
          tile_data(1:Nt, 4) = max(min(cat_diagS_avg%sfmc/cat_param%poros,1.),0.) ! sm_surface_wetness                dimensionless
          tile_data(1:Nt, 5) = max(min(cat_diagS_avg%rzmc/cat_param%poros,1.),0.) ! sm_rootzone_wetness               dimensionless
          tile_data(1:Nt, 6) = max(min(cat_diagS_avg%prmc/cat_param%poros,1.),0.) ! sm_profile_wetness                dimensionless
                                                                                                                      
          tile_data(1:Nt, 7) = cat_diagS_avg%tsurf                                ! surface_temp                      K          
          tile_data(1:Nt, 8) = cat_diagS_avg%tp(1)                                ! soil_temp_layer1                  K 
          tile_data(1:Nt, 9) = cat_diagS_avg%tp(2)                                ! soil_temp_layer2                  K 
          tile_data(1:Nt,10) = cat_diagS_avg%tp(3)                                ! soil_temp_layer3                  K 
          tile_data(1:Nt,11) = cat_diagS_avg%tp(4)                                ! soil_temp_layer4                  K 
          tile_data(1:Nt,12) = cat_diagS_avg%tp(5)                                ! soil_temp_layer5                  K 
          tile_data(1:Nt,13) = cat_diagS_avg%tp(6)                                ! soil_temp_layer6                  K 
                                                                                                                      
          do n=1,Nt                                                                                               
                                                                                                                      
             tile_data(n,14) = sum(cat_progn_avg(n)%wesn(1:N_snow))               ! snow_mass                         kg m-2
             tile_data(n,15) = sum(cat_progn_avg(n)%sndz(1:N_snow))               ! snow_depth                        m 
                                                                                                                      
          end do                                                                                                      
                                                                                                                      
          tile_data(1:Nt,16) = cat_diagF_avg%evap                                 ! land_evapotranspiration_flux      kg m-2 s-1
          tile_data(1:Nt,17) = cat_diagF_avg%runsrf                               ! overland_runoff_flux              kg m-2 s-1
          tile_data(1:Nt,18) = cat_diagF_avg%bflow                                ! baseflow_flux                     kg m-2 s-1
          tile_data(1:Nt,19) = cat_diagF_avg%snmelt                               ! snow_melt_flux                    kg m-2 s-1
          tile_data(1:Nt,20) = cat_diagF_avg%qinfil                               ! soil_water_infiltration_flux      kg m-2 s-1
                                                                                                                      
                                                                                                                      
          ! Note: ar1+ar2+ar4=1 but need FRSAT+FRUNST+FRWLT+FRSNO=1                                                   
                                                                                                                      
          tmpreal(1:Nt) = max(min((1.-cat_diagS_avg%asnow),1.),0.)  ! precompute snow-free fraction                   
                                                                                                                      
          tile_data(1:Nt,21) = max(min(tmpreal*cat_diagS_avg%ar1,  1.),0.)        ! land_fraction_saturated           dimensionless
          tile_data(1:Nt,22) = max(min(tmpreal*cat_diagS_avg%ar2,  1.),0.)        ! land_fraction_unsaturated         dimensionless
                                                                                                                      
          tmpreal = tmpreal-tile_data(1:Nt,21)-tile_data(1:Nt,22)   ! compute wilting fraction                        
                                                                                                                      
          tile_data(1:Nt,23) = max(min(tmpreal,                    1.),0.)        ! land_fraction_wilting             dimensionless
                                                                                                                      
          tile_data(1:Nt,24) = max(min(        cat_diagS_avg%asnow,1.),0.)        ! land_fraction_snow_covered        dimensionless
                                                                                                                      
          tile_data(1:Nt,25) = cat_diagF_avg%shflux                               ! heat_flux_sensible                W m-2
          tile_data(1:Nt,26) = cat_diagF_avg%lhflux                               ! heat_flux_latent                  W m-2
          tile_data(1:Nt,27) = cat_diagF_avg%ghflux                               ! heat_flux_ground                  W m-2
                                                                                                                      
          tile_data(1:Nt,28) = met_force_avg%SWdown - cat_diagF_avg%swup          ! net_downward_shortwave_flux       W m-2
          tile_data(1:Nt,29) = met_force_avg%LWdown - cat_diagF_avg%lwup          ! net_downward_longwave_flux        W m-2

          tile_data(1:Nt,30) = met_force_avg%SWdown                               ! radiation_shortwave_downward_flux W m-2
          tile_data(1:Nt,31) = met_force_avg%LWdown                               ! radiation_longwave_absorbed_flux  W m-2

          tile_data(1:Nt,32) = met_force_avg%Rainf + met_force_avg%Snowf          ! precipitation_total_surface_flux  kg m-2 s-1
          tile_data(1:Nt,33) = met_force_avg%Snowf                                ! snowfall_surface_flux             kg m-2 s-1

          tile_data(1:Nt,34) = met_force_avg%Psurf                                ! surface_pressure                  Pa
          
          tile_data(1:Nt,35) = met_force_avg%RefH                                 ! height_lowatmmodlay               m
          tile_data(1:Nt,36) = met_force_avg%Tair                                 ! temp_lowatmmodlay                 K
          tile_data(1:Nt,37) = met_force_avg%Qair                                 ! specific_humidity_lowatmmodlay    kg kg-1
          tile_data(1:Nt,38) = met_force_avg%Wind                                 ! windspeed_lowatmmodlay            m s-1
          
          tile_data(1:Nt,39) = veg_param_avg%grn                                  ! vegetation_greenness_fraction     dimensionless
          tile_data(1:Nt,40) = veg_param_avg%lai                                  ! leaf_area_index                   m2 m-2 

          
       case (7,8)  ! N_out_fields_inst = 4; N_out_fields_tavg = 4, 6
          
          ! Specs for SMAP Nature Run v03 - reichle, 11 Dec 2013
          ! Modified                      - reichle,  6 Feb 2014
          ! Bug fix: units of tp(1)       - reichle, 19 Feb 2014
          ! Added tpsn(1)                 - reichle,  4 Mar 2014
      
          ! compute snow mass
          
          do n=1,Nt                                                                                               
             
             SWE(n) = sum(cat_progn_avg(n)%wesn(1:N_snow))
             
          end do

          ! generate different output for "inst" and "tavg" files
          
          select case(interval)
             
          case ('i')     ! instantaneous
             
             ! convert Catchment model variables into inputs suitable for the mwRTM 
             ! NOTE: input tp must be in degree Celsius!

             call catch2mwRTM_vars( Nt, cat_param%vegcls, cat_param%poros,    &
                  mwRTM_param%poros, cat_diagS_avg%sfmc, cat_diagS_avg%tsurf,     &
                  cat_diagS_avg%tp(1)-Tzero, sfmc_mwRTM, tsoil_mwRTM )
             
             ! calculate brightness temperatures
             ! (tau-omega model as in De Lannoy et al. 2013 [doi:10.1175/JHM-D-12-092.1]
             !  but without Pellarin atmospheric corrections)

             freq           = 1.41e9
             
             inc_angle      = 40.

             incl_atm_terms = .false.

             call mwRTM_get_Tb(Nt, freq, inc_angle, mwRTM_param, tile_coord%elev,      &
                  veg_param_avg%lai, sfmc_mwRTM, tsoil_mwRTM, SWE, met_force_avg%Tair, &
                  incl_atm_terms,                                                      &
                  Tb_h, Tb_v )
             
             ! fill tile_data
                             
             tile_data(1:Nt,1) = cat_diagS_avg%tsurf                         ! surface_temp                      K              
             tile_data(1:Nt,2) = cat_diagS_avg%tp(1)                         ! soil_temp_layer1                  K             
             tile_data(1:Nt,3) = Tb_h                                        ! tb_h [at above freq, inc_angle]   K             
             tile_data(1:Nt,4) = Tb_v                                        ! tb_v [at above freq, inc_angle]   K
             
             if (out_collection_ID==8) then
                
                tile_data(1:Nt,5) = cat_diagS_avg%tpsn(1)                       ! snow_temp_layer1                  K
                
             end if
             
          case default  ! time-average output (any averaging interval)
             
             tile_data(1:Nt,1) = cat_diagS_avg%sfmc                          ! sm_surface                        m3 m-3              
             tile_data(1:Nt,2) = cat_diagS_avg%rzmc                          ! sm_rootzone                       m3 m-3              
             tile_data(1:Nt,3) = cat_diagS_avg%prmc                          ! sm_profile                        m3 m-3              
             tile_data(1:Nt,4) = cat_diagS_avg%tp(1)                         ! soil_temp_layer1                  K             

             if (out_collection_ID==8) then
                
                tile_data(1:Nt,5) = SWE                                         ! snow_mass                         kg m-2
                tile_data(1:Nt,6) = met_force_avg%Rainf + met_force_avg%Snowf   ! precipitation_total_surface_flux  kg m-2 s-1

             end if                

          end select  
          
       case (9)    ! N_out_fields_inst = 6;  N_out_fields_tavg = 2
          
          ! for L-band mwRTM calibration (Dec 2013)
          ! Renamed from 8 to 9 - reichle,  6 Feb 2014
                    
          ! generate different output for "inst" and "tavg" files
          
          select case(interval)
             
          case ('i')     ! instantaneous
                       
             tile_data(1:Nt,1) = cat_diagS_avg%sfmc                          ! sm_surface                        m3 m-3
             tile_data(1:Nt,2) = cat_diagS_avg%rzmc                          ! sm_rootzone                       m3 m-3 
             tile_data(1:Nt,3) = cat_diagS_avg%tsurf                         ! surface_temp                      K   
             tile_data(1:Nt,4) = cat_diagS_avg%tp(1)                         ! soil_temp_layer1                  K   
             tile_data(1:Nt,5) = met_force_avg%Tair                          ! temp_lowatmmodlay                 K
             tile_data(1:Nt,6) = veg_param_avg%lai                           ! leaf_area_index                   m2 m-2 

          case default  ! time-average output (any averaging interval)
             
             ! compute snow mass
             
             do n=1,Nt                                                                                               
                
                SWE(n) = sum(cat_progn_avg(n)%wesn(1:N_snow))
                
             end do
             
             tile_data(1:Nt,1) = SWE                                         ! snow_mass                         kg m-2
             tile_data(1:Nt,2) = met_force_avg%Rainf + met_force_avg%Snowf   ! precipitation_total_surface_flux  kg m-2 s-1                     
             
          end select
          
       case (11)    ! N_out_fields_inst = N_out_fields_tavg = 5

          ! for SMOS pre-processing using Gabrielle De Lannoy's matlab routines (Dec 2014)
          ! - reichle, 30 Mar 2015
          
          tile_data(1:Nt,1) = cat_diagS_avg%tsurf                         ! surface_temp            K   
          tile_data(1:Nt,2) = cat_diagS_avg%tp(1)                         ! soil_temp_layer1        K   
          tile_data(1:Nt,3) = met_force_avg%Psurf                         ! surface_pressure        Pa
          tile_data(1:Nt,4) = cat_diagF_avg%t2m                           ! temp_2m                 K
          tile_data(1:Nt,5) = cat_diagF_avg%q2m                           ! specific_humidity_2m    kg kg-1
            
       case default
          
          err_msg = 'unknown out_collection_ID'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          
       end select     ! case (out_collection_ID)
       
       ! ----------------------------------------------------------------
       !
       ! re-initialize
              
       do n=1,Nt
          
          cat_progn_avg(n) = 0.
          cat_diagS_avg(n) = 0.
          cat_diagF_avg(n) = 0.
          met_force_avg(n) = 0.          
          veg_param_avg(n) = 0.
          bal_diagn_avg(n) = 0.
          
       end do
       
       ! ----------------------------------------------------------------
       
    else
       
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'bad option')
       
    end if       ! option 'ini', 'add', or 'out'

  end subroutine output_calcs

  ! ********************************************************************
  
  subroutine output_write( out_tile, out_grid, fname_tile, fname_grid,  &
       out_collection_ID, N_out_fields,                                 &
       N_catl, N_catf, N_land_mask, tile_coord_f, tile_grid_f,          &
       N_catl_vec, low_ind, land_mask_i, land_mask_j, tile_data_l )
    
    ! reichle, 23 Dec 2011
    
    ! revised output subroutines to accomodate LAI-weighted greenness (GRN) and
    ! for general clean-up
    
    implicit none
    
    logical,                      intent(in) :: out_tile, out_grid
    
    character(300),               intent(in) :: fname_tile, fname_grid
    
    integer,                      intent(in) :: out_collection_ID, N_out_fields
    integer,                      intent(in) :: N_catl, N_catf, N_land_mask
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord_f ! intent(in), N_catf
        
    type(grid_def_type),          intent(in) :: tile_grid_f
        
    integer, dimension(:), intent(in) :: N_catl_vec, low_ind
    
    integer, dimension(tile_grid_f%N_lon*tile_grid_f%N_lat), intent(in) :: &
         land_mask_i, land_mask_j    
    
    real,    dimension(N_catl,N_out_fields),                 intent(in) :: tile_data_l
    
    ! local variables

    integer                                                :: i
    real,   dimension(N_catf)                              :: tile_data_f
    real,   dimension(N_catf)                              :: tile_data_f_tmp
    real,   dimension(tile_grid_f%N_lon,tile_grid_f%N_lat) :: grid_data
    real,   dimension(tile_grid_f%N_lon,tile_grid_f%N_lat) :: grid_data_tmp
    
    ! ---------------------------------------------------------
    
    do i=1,N_out_fields   ! write output one field at a time
       
       ! gatherv tile data from local to full domain and map to grid
       
       call l2f_real( N_catl, N_catf, N_catl_vec, low_ind,                          & 
            tile_data_l(:,i), tile_data_f)
       
       if (master_proc)                                                             &
            call tile2grid( N_catf, tile_coord_f, tile_grid_f, tile_data_f,         &
            grid_data,                                                              &
            no_data_value=nodata_generic, no_data_tol=nodata_tol_generic)
       
       ! special case: LAI-weighted greenness for MERRA-Land and SMAP file specs
       
       if ( (i==1)                        .and.                &
            (                                                  &
            out_collection_ID==4  .or.                         &
            out_collection_ID==5  .or.                         &
            out_collection_ID==6                               &
            )                                     ) then
          
          ! i=1: GRN             
          ! i=2: LAI             

          ! gatherv LAI into tile_data_f_tmp
          
          call l2f_real( N_catl, N_catf, N_catl_vec, low_ind, tile_data_l(:,2),     &
               tile_data_f_tmp)
          
          if (master_proc)  then

             ! get gridded LAI
             
             call tile2grid( N_catf, tile_coord_f, tile_grid_f, tile_data_f_tmp,    &
                  grid_data_tmp,                                                    &
                  no_data_value=nodata_generic, no_data_tol=nodata_tol_generic)
             
             ! get LAI-weighted GRN
             
             tile_data_f_tmp = tile_data_f*tile_data_f_tmp   ! = GRN*LAI
             
             call tile2grid( N_catf, tile_coord_f, tile_grid_f, tile_data_f_tmp,    &
                  grid_data,                                                        &
                  no_data_value=nodata_generic, no_data_tol=nodata_tol_generic)
             
             ! set to no-data-values when gridded LAI is zero or no-data, 
             !  otherwise normalize, ie., compute  [GRN*LAI] / LAI
             ! [edited to avoid division by zero, -reichle+csdraper, 29 Jan 2016]

             where (                                                                &
                  grid_data_tmp < 1.e-10                                      .or.  &
                  abs(grid_data_tmp-nodata_generic)<nodata_tol_generic              &
                  )                       
                
                grid_data = nodata_generic

             elsewhere
                
                grid_data = grid_data/grid_data_tmp    ! normalize                

             end where
             
          end if
          
       end if
       
       ! write field i (ie, tile_data_f and grid_data) to file
       
       if (master_proc)                                                             &
            call write_output_field( out_tile, out_grid,                            &
            fname_tile, fname_grid, N_out_fields, i, N_catf, tile_grid_f,           &
            tile_data_f, grid_data, N_land_mask, land_mask_i, land_mask_j )
       
    end do
    
  end subroutine output_write

  ! ********************************************************************
  
  subroutine write_output_field( out_tile, out_grid, fname_tile, fname_grid,        &
       N_out_fields, out_field_num, N_catd, tile_grid_d,                            &
       tile_data, grid_data, N_land_mask, land_mask_i, land_mask_j )
    
    ! write (one field of) tile and/or gridded output for instantaneous or 
    ! time average data
    !
    ! reichle, 29 Sep 2009
    ! reichle, 23 Dec 2011 - revised for MERRA-Land file specs
    !
    ! ---------------------------------------------------------------------
    
    implicit none
    
    logical,               intent(in) :: out_tile, out_grid
    
    character(300),        intent(in) :: fname_tile, fname_grid
    
    integer,               intent(in) :: N_out_fields, out_field_num, N_catd
    
    type(grid_def_type),   intent(in) :: tile_grid_d
    
    real,                  intent(in), dimension(N_catd) :: tile_data
    
    real, intent(in), dimension(tile_grid_d%N_lon,tile_grid_d%N_lat) :: grid_data
    
    integer,               intent(in) :: N_land_mask
    
    integer, dimension(tile_grid_d%N_lon*tile_grid_d%N_lat), intent(in) :: &
         land_mask_i, land_mask_j
    
    ! local variables
    
    integer :: unitnumber_tile = 10
    integer :: unitnumber_grid = 11
    integer :: n
    
    type(grid_def_type)                                  :: tmp_grid_def
    
    real, dimension(N_catd)                              :: tile_data_shaved
    
    real, dimension(tile_grid_d%N_lon*tile_grid_d%N_lat) :: tmpvec
        
    ! --------------------------------------------------------
    !
    ! write tile output
    
    if (out_tile) then
       
       if (out_field_num==1) then
          
          open(unitnumber_tile, file=fname_tile, form='unformatted', &
               action='write', status='unknown')
          
       else
          
          inquire(file=fname_tile,number=unitnumber_tile)
          
       end if
       
       ! degrade least significant bits in return for better gzip compressoion

       call shave_bits(N_catd, tile_data, tile_data_shaved)
       
       ! write to file
       
       write(unitnumber_tile) (tile_data_shaved(n), n=1,N_catd)
       
       if (out_field_num==N_out_fields)  close(unitnumber_tile, status='keep')
       
    end if
    
    ! --------------------------------------------------------
    !
    ! write (compressed) gridded output
    
    if (out_grid) then
       
       if (out_field_num==1) then
          
          open(unitnumber_grid, file=fname_grid, form='unformatted', &
               action='write', status='unknown')
          
          ! write file header (one additional header line as of 21 May 2010!!)
          
          tmp_grid_def =  tile_grid_d
          
          call io_grid_def_type( 'w', unitnumber_grid, tmp_grid_def )
          
          write (unitnumber_grid) N_land_mask
          
          write (unitnumber_grid) (land_mask_i(n), n=1,N_land_mask)
          write (unitnumber_grid) (land_mask_j(n), n=1,N_land_mask)

       else
          
          inquire(file=fname_grid,number=unitnumber_grid)
          
       end if
       
       ! ------------------------------------------------------------
       !
       ! compress grid data to land-only vector
       
       do n=1,N_land_mask
          
          tmpvec(n) = grid_data(land_mask_i(n),land_mask_j(n))
          
       end do
       
       ! degrade least significant bits in return for better gzip compressoion
       
       call shave_bits(N_land_mask, tmpvec, tmpvec)
       
       ! write to file
       
       write (unitnumber_grid) (tmpvec(n), n=1,N_land_mask)
       
       ! close file after last field has been written
       
       if (out_field_num==N_out_fields)  close(unitnumber_grid, status='keep')
       
    end if
    
  end subroutine write_output_field
  
  ! ********************************************************************
  
  subroutine shave_bits(N_data, data, data_shaved)
    
    ! reichle, 22 Dec 2011
    
    implicit none

    integer,                    intent(in) :: N_data
    
    real,    dimension(N_data), intent(in)  :: data
    real,    dimension(N_data), intent(out) :: data_shaved
    
    ! local variables

    integer :: rc
    
    integer, external :: ShaveMantissa32

    character(len=*), parameter :: Iam = 'shave_bits'
    character(len=400) :: err_msg

    ! --------------------------------------------------------
    
    ! int ShaveMantissa32 ( a, ain, len, xbits, has_undef, undef, chunksize )

    ! int32   len;        /* size of a[] */
    ! int     xbits;      /* number of bits to excludes from mantissa */
    ! int     has_undef;  /* whether field has missing (undef) values */ 
    ! int32   chunksize;  /* find mid range over chunksize chunks     */ 
    ! float32 undef;      /* missing (undefined) value */
    ! float32 ain[];      /* input array */
    !
    ! float32 a[];    // output "shaved" array; can share storage with ain[]

    rc = 1
    
    if (N_bits_shaved==0) then
                     
       data_shaved = data
       
       rc = 0
       
    elseif (N_bits_shaved>0 .and. N_bits_shaved<=12) then
       
       rc = ShaveMantissa32( data_shaved, data, N_data,     &
            N_bits_shaved, .true., nodata_generic, N_data )
       
    else
       
       err_msg = 'shaving more than 12 bits is not recommended'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    if (rc/=0) call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'rc/=0')
    
  end subroutine shave_bits
  
  ! ********************************************************************
  
  subroutine get_ensstd_filenames( out_choice, date_time_new, work_path, exp_id, &
       interval, out_tile, out_grid, fname_tile, fname_grid, out_dtstep_xhourly )

    implicit none
    
    type(out_choice_type), intent(in)    :: out_choice
    
    type(date_time_type),  intent(in)    :: date_time_new
    
    character(200),        intent(in)    :: work_path
    character(40),         intent(in)    :: exp_id
    
    ! what averaging interval is used? (need to know to construct file name)
    !
    ! inst   : interval = 'i'  -- for inst output ONLY call with option 'out'
    ! xhourly: interval = 'x'
    ! daily  : interval = 'd'
    ! pentad : interval = 'p'
    ! monthly: interval = 'm'
    
    character,      intent(in)           :: interval
    
    logical,        intent(out)          :: out_tile, out_grid
    
    character(300), intent(out)          :: fname_tile, fname_grid

    integer,        intent(in), optional :: out_dtstep_xhourly
    
    ! local variables
    
    character(40)        :: dir_name, file_tag
    
    type(date_time_type) :: date_time_tmp
    
    character(len=*), parameter :: Iam = 'get_ensstd_filenames'
    character(len=400) :: err_msg

    ! ------------------------------------------------------

    out_tile = .false.
    out_grid = .false.
    
    dir_name = 'ana'
    
    select case (interval)
       
    case ('i')     ! instantaneous
       
       if (out_choice%inst%tile) then
          
          out_tile = .true.
          
          file_tag = 'ldas_tile_inst_ensstd'
          
          fname_tile = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_new,                              &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
       if (out_choice%inst%grid) then
          
          out_grid = .true.
          
          file_tag = 'ldas_grid_inst_ensstd'
          
          fname_grid = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_new,                              &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
    case ('x')     ! xhourly
       
       if (.not. present(out_dtstep_xhourly)) then
          err_msg = 'optional input argument missing'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       ! shift date/time so that mid-point of time averaging interval is
       ! used for time-tagging the output file
       
       date_time_tmp = date_time_new
       
       call augment_date_time( -out_dtstep_xhourly/2, date_time_tmp )
       
       if (out_choice%xhourly%tile) then
          
          out_tile = .true.
          
          file_tag = 'ldas_tile_xhourly_ensstd'
          
          fname_tile = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_tmp,                              &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
       if (out_choice%xhourly%grid) then
          
          out_grid = .true.
          
          file_tag = 'ldas_grid_xhourly_ensstd'
          
          fname_grid = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_tmp,                              &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
    case ('d')     ! daily
       
       if (out_choice%daily%tile) then
          
          out_tile = .true.
          
          file_tag = 'ldas_tile_daily_ensstd'  
          
          fname_tile = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_new, option=4,                    &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
       if (out_choice%daily%grid) then
          
          out_grid = .true.
          
          file_tag = 'ldas_grid_daily_ensstd'
          
          fname_grid = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_new, option=4,                    &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
    case ('p')     ! pentad
       
       if (out_choice%pentad%tile) then
          
          out_tile = .true.
          
          file_tag = 'ldas_tile_pentad_ensstd'
          
          fname_tile = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_new, option=3,                    &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
       if (out_choice%pentad%grid) then
          
          out_grid = .true.
          
          file_tag = 'ldas_grid_pentad_ensstd'
          
          fname_grid = get_io_filename( work_path, exp_id, file_tag,  &
               date_time=date_time_new, option=3,                     &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
    case ('m')     ! monthly
       
       if (out_choice%monthly%tile) then
          
          out_tile = .true.
          
          file_tag = 'ldas_tile_monthly_ensstd'
          
          fname_tile = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_new, option=2,                    &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
       if (out_choice%monthly%grid) then
          
          out_grid = .true.
          
          file_tag = 'ldas_grid_monthly_ensstd'
          
          fname_grid = get_io_filename( work_path, exp_id, file_tag, &
               date_time=date_time_new, option=2,                    &
               dir_name=dir_name, ens_id=-1 )
          
       end if
       
    case default
       
       err_msg = 'unknown averaging interval'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
  end subroutine get_ensstd_filenames
  
  ! ********************************************************************
  
  subroutine check_output_times( out_dtstep, date_time_new, date_time_old, end_time, &
       out_time )
    
    ! reichle, 2 Oct 2009
    
    implicit none
    
    type(out_dtstep_type),      intent(in)  :: out_dtstep
    
    type(date_time_type),       intent(in)  :: date_time_new, date_time_old, end_time
    
    type(out_choice_time_type), intent(out) :: out_time
    
    ! local variables
    
    integer :: secs_in_day

    logical :: new_day, new_pentad, new_month
    
    ! --------------------------------------------------
    
    out_time%rstrt         = .false.
    out_time%inst          = .false. 
    out_time%xhourly       = .false.           
    out_time%daily         = .false.           
    out_time%pentad        = .false.           
    out_time%monthly       = .false.           
    out_time%any_non_rstrt = .false. 
    
    secs_in_day = date_time_new%hour*3600 + date_time_new%min*60 + date_time_new%sec
    
    new_day     = (secs_in_day==0)
    
    new_pentad  = (date_time_new%pentad /= date_time_old%pentad)  
    
    new_month   = (new_day .and. date_time_new%day==1)
    
    ! check if rstrt output is needed

    ! write restarts
    !  - at beginning of month (always)
    !  - at appropriate time steps as requested
    !  - at end_time of simulation (always)
    
    if ( new_month                                                          &
         .or.                                                               &
         (out_dtstep%rstrt>0 .and. mod(secs_in_day,out_dtstep%rstrt)==0)    & 
         .or.                                                               &
         datetime_eq_refdatetime(date_time_new,end_time)                    &
         )                                                                  &
         out_time%rstrt = .true.
    
    ! inst
    
    if (out_dtstep%inst/=0) then
       
       if (mod(secs_in_day,out_dtstep%inst)   ==0)     out_time%inst = .true.
       
    end if
    
    ! xhourly

    if (out_dtstep%xhourly/=0) then
       
       if (mod(secs_in_day,out_dtstep%xhourly)==0)     out_time%xhourly = .true.
       
    end if
    
    ! daily, pentad, monthly

    if (new_day)   out_time%daily   = .true.
    
    if (date_time_new%pentad /= date_time_old%pentad)  out_time%pentad  = .true.
    
    if (new_month) out_time%monthly = .true.
    
    ! is there any non-rstrt output? 
    
    if ( out_time%inst    .or. & 
         out_time%xhourly .or. & 
         out_time%daily   .or. & 
         out_time%pentad  .or. & 
         out_time%monthly        )   out_time%any_non_rstrt = .true. 
    
  end subroutine check_output_times
  
  ! ********************************************************************
  
  subroutine get_land_mask_ij( N_catd, tile_coord, tile_grid_d, &
       N_land_mask, land_mask_i, land_mask_j )
    
    implicit none
    
    integer, intent(in) :: N_catd
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid_d
        
    integer, intent(inout) :: N_land_mask
    
    integer, dimension(tile_grid_d%N_lon*tile_grid_d%N_lat), intent(out) :: &
         land_mask_i, land_mask_j
        
    ! -----------------------------------------------------
    
    ! local variables
    
    integer :: i, j

    real, dimension(:), allocatable :: tile_data_tmp

    real, dimension(tile_grid_d%N_lon,tile_grid_d%N_lat) :: grid_data
        
    ! -------------------------------------------------------------
    
    ! map a vector full of "good" tiles to grid
    
    allocate(tile_data_tmp(N_catd))
    
    tile_data_tmp = (nodata_generic + 1000.*nodata_tol_generic)
    
    call tile2grid( N_catd, tile_coord, tile_grid_d, tile_data_tmp,               &
         grid_data, no_data_value=nodata_generic, no_data_tol=nodata_tol_generic)
    
    deallocate(tile_data_tmp)
    
    ! see which grid boxes are "good"
    
    N_land_mask = 0
    
    do j=1,tile_grid_d%N_lat
       do i=1,tile_grid_d%N_lon
          
          if (abs(grid_data(i,j)-nodata_generic)>nodata_tol_generic) then
             
             N_land_mask = N_land_mask+1
             
             land_mask_i(N_land_mask) = i
             land_mask_j(N_land_mask) = j
             
          end if
          
       end do
    end do
    
  end subroutine get_land_mask_ij

  ! ********************************************************************

end module clsm_ensdrv_out_routines

! *********** EOF ******************************************************
