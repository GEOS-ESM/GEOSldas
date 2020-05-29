module LDAS_ensdrv_Globals

  ! global parameters for LDAS ens driver
  !
  ! reichle, 25 Mar 2004
  ! reichle,  6 May 2005
  ! reichle, 29 Nov 2010 - deleted N_outselect (obsolete)
  !                      - added sfc_turb_scheme (choose Louis or Helfand Monin-Obukhov)
  ! reichle,  5 Apr 2013 - removed N_out_fields as global parameter
  ! wjiang+reichle, 
  !          21 May 2020 - added "LDAS_is_nodata" function, checks if "nodata_generic" or "MAPL_UNDEF"
  
  use, intrinsic :: iso_fortran_env, only : output_unit

  use MAPL_BaseMod,                  only : MAPL_UNDEF
  
  implicit none
  
  private

  public :: nodata_generic
  public :: nodata_tolfrac_generic
  public :: nodata_tol_generic
  public :: LDAS_is_nodata
  public :: logunit
  public :: logit
  public :: master_logit
  public :: log_master_only
  
  public :: echo_clsm_ensdrv_glob_param
  public :: write_status

  ! ----------------------------------------------------------------------
      
  ! generic no-data-value
  
  real, parameter :: nodata_generic         = -9999.
  real, parameter :: nodata_tolfrac_generic = 1.e-4
  
  real :: nodata_tol_generic     = abs(nodata_generic*nodata_tolfrac_generic)
  real :: MAPL_UNDEF_tol_generic = abs(MAPL_UNDEF    *nodata_tolfrac_generic) 

  ! ----------------------------------------------------------------
  !
  ! log file

  ! Avoid I/O buffering for the log file by setting "logunit" to stdout,
  ! then redirect stdout to the log file via driver script ("ldsetup").
  ! In case of an error during execution, this allows capturing of all log messages 
  ! until the job terminates.
  !
  ! NOTE: "logunit=stdout" is disabled if log messages are requested from *all* processors
  !       (that is, for "log_master_only=.false.") to avoid garbled output

  integer, parameter :: logunit         = output_unit ! defined in iso_fortran_env
  
  logical, parameter :: log_master_only = .true.
  
  logical            :: logit,master_logit

  
contains
  
  subroutine echo_clsm_ensdrv_glob_param()
    
    ! echo all global parameters 
    
    ! call only AFTER opening log file!!!
    
    implicit none
    
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)
    write (logunit,*) 'echo_clsm_ensdrv_glob_param():'
    write (logunit,*)
    write (logunit,*) 'nodata_generic          = ',   nodata_generic
    write (logunit,*)
    write (logunit,*) 'nodata_tolfrac_generic  = ',   nodata_tolfrac_generic
    write (logunit,*)
    write (logunit,*) 'nodata_tol_generic      = ',   nodata_tol_generic
    write (logunit,*)
    write (logunit,*) 'logunit                 = ',   logunit
    write (logunit,*)
    write (logunit,*) 'log_master_only         = ',   log_master_only
    write (logunit,*)
    write (logunit,*) 'logit                   = ',   logit
    write (logunit,*)
    
    write (logunit,*)
    write (logunit,*) 'end echo_clsm_ensdrv_glob_param()'
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)
    
  end subroutine echo_clsm_ensdrv_glob_param
  
  ! ********************************************************************
  
  subroutine write_status(lenkf_status)
    
    ! write status message (success/failure) to designated file 
    !
    ! hardwired filename, used by ADAS scripts
    
    ! Draper, reichle, 27 Feb 2012
    
    implicit none
    
    logical, intent(in)  :: lenkf_status 

    ! --------------------------------------------------
    
    open( unit=10, file='lenkf_job_completed.txt' )
    
    if (lenkf_status) then 
       
       write (10,*) 'SUCCEEDED'
       
    else

       write (10,*) 'FAILED'

    endif
    
    close(unit=10)
    
  end subroutine write_status

  ! ********************************************************************
  
  elemental   function LDAS_is_nodata(data) result(no_data)
    
    real,   intent(in) :: data
    logical            :: no_data
    
    no_data =                                                         &
         ( abs(data-nodata_generic) < nodata_tol_generic    ) .or.    &
         ( abs(data-MAPL_UNDEF)     < MAPL_UNDEF_tol_generic)     
    
  end function LDAS_is_nodata

  ! *************************************************************
  
end module LDAS_ensdrv_Globals


!======== EOF ==============================================================
