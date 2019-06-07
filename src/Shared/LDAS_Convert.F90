#include "MAPL_Generic.h"

module LDAS_ConvertMod

  use ESMF
  use MAPL_mod
  use LDAS_DateTimeMod, only: date_time_type
  use LDAS_DateTimeMod, only: get_dofyr_pentad

  implicit none

  private

  public :: esmf2ldas

  interface esmf2ldas
     module procedure esmf2ldas_time
  end interface esmf2ldas

contains

  subroutine esmf2ldas_time(esmf_dt, ldas_dt, rc)
    
    type(ESMF_Time), intent(in) :: esmf_dt
    type(date_time_type), intent(out) :: ldas_dt
    integer, optional, intent(out) :: rc

    character(len=*), parameter :: Iam = 'emsf2ldas_time'
    integer :: status

    call ESMF_TimeGet(                                                          &
         esmf_dt,                                                               &
         YY=ldas_dt%year,                                                       &
         MM=ldas_dt%month,                                                      &
         DD=ldas_dt%day,                                                        &
         H=ldas_dt%hour,                                                        &
         M=ldas_dt%min,                                                         &
         S=ldas_dt%sec,                                                         &
         rc=status                                                              &
         )
    VERIFY_(status)

    call get_dofyr_pentad( ldas_dt )

    RETURN_(ESMF_SUCCESS)

  end subroutine esmf2ldas_time

end module LDAS_ConvertMod
