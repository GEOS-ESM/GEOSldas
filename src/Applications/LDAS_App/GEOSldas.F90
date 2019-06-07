#define I_AM_MAIN

#include "MAPL_Generic.h"

program LDAS_Main

  use MAPL_Mod
  use GEOS_LDASGridCompMod, only:  ROOT_SetServices => SetServices
  
  implicit none
  
  integer :: status
  character(len=9) :: Iam="LDAS_Main"
  
  call MAPL_CAP(ROOT_SetServices, FinalFile='EGRESS.ldas', rc=status)
  VERIFY_(status)
  
end program LDAS_Main
