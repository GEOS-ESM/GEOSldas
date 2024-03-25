#define I_AM_MAIN

#include "MAPL_Generic.h"

program LDAS_Main


   ! !USES:
   use MAPL
   use GEOS_LDASGridCompMod, only:  ROOT_SetServices => SetServices

   implicit none

   character(len=*), parameter :: Iam = "LDAS_Main"
   type (MAPL_Cap) :: cap
   type (MAPL_FargparseCLI) :: cli
   type (MAPL_CapOptions) :: cap_options
   integer :: status

!EOP
!----------------------------------------------------------------------
!BOC

   cli = MAPL_FargparseCLI()
   cap_options = MAPL_CapOptions(cli)
   cap_options%egress_file = 'EGRESS.ldas'

   cap = MAPL_Cap('LDAS', ROOT_SetServices, cap_options = cap_options)
   call cap%run(_RC)

   !call MAPL_CAP(ROOT_SetServices, FinalFile='EGRESS.ldas', rc=status)
   !VERIFY_(status)

end program LDAS_Main
