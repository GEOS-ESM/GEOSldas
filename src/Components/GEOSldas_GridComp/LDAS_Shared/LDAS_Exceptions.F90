! clsm_exceptions.F90

! pchakrab, xx July 2014

module LDAS_ExceptionsMod

  use ESMF
  use MAPL_Mod

  use, intrinsic :: iso_fortran_env, only: logunit => output_unit

  implicit none

  private

  integer, parameter, public :: LDAS_GENERIC_ERROR = 3000
  integer, parameter, public :: LDAS_FILE_NOT_FOUND = 3001
  !integer, parameter, public :: LDAS_INVALID_VALUE = 3002
  integer, parameter, public :: LDAS_GENERIC_WARNING = 6000

  ! more error/warning codes here

  public :: ldas_abort
  public :: ldas_warn

contains

  subroutine ldas_abort(err_code, calling_fn, message)
    ! input/output
    integer, intent(in) :: err_code
    character(len=*), intent(in) :: calling_fn
    character(len=*), intent(in) :: message

    ! local
    character(len=10) :: err_code_str ! largest 4B integer has 10 digits
    type(ESMF_VM) :: vm
    integer :: mpierr
    integer :: comm
    integer :: status

    ! write status (failed) file
    call write_status(.false.)

    write(err_code_str, '(i10)') err_code ! err_code from int to str
    write(logunit, *) 'LDAS ERROR (' // &
         trim(adjustl(err_code_str)) // ') from ' // &
         trim(calling_fn) // ': ' // trim(message)

#ifdef LDAS_MPI
    ! abort
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(status)
    call ESMF_VmGet(vm, mpicommunicator=comm, rc=status)
    VERIFY_(status)
    call MPI_Abort(comm, err_code, mpierr)
#else
    stop
#endif

  end subroutine ldas_abort


  subroutine ldas_warn(warn_code, calling_fn, message)
    ! input/output
    integer, intent(in) :: warn_code
    character(len=*), intent(in) :: calling_fn
    character(len=*), intent(in) :: message

    ! local
    character(len=10) :: warn_code_str ! largest 4B integer has 10 digits

    write(warn_code_str, '(i10)') warn_code
    write(logunit, *) 'LDAS WARNING (' // &
         trim(adjustl(warn_code_str)) // ') from ' // &
         trim(calling_fn) // ': ' // trim(message)

  end subroutine ldas_warn


  subroutine write_status(lenkf_status)

    ! write status message (success/failure) to designated file
    ! hardwired filename, used by ADAS scripts

    ! Draper, reichle, 27 Feb 2012

    implicit none

    logical, intent(in)  :: lenkf_status

    open( unit=10, file='lenkf_job_completed.txt' )

    if (lenkf_status) then
       write (10,*) 'SUCCEEDED'
    else
       write (10,*) 'FAILED'
    endif

    close(unit=10)

  end subroutine write_status


end module LDAS_ExceptionsMod
