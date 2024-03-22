
! how to use :
! ./preprocess_ldas option arg1 arg2 arg3

program main

  use preprocess_ldas_routines,     ONLY:    &
       createf2g,                            &
       createLocalTilefile,                  &
       createLocalBC,                        &
       createLocalVegRestart,                &
       createLocalmwRTMRestart,              &
       createLocalCatchRestart,              &
       correctEase,                          &
       convert_pert_rst,                     &
       optimize_latlon
  
  implicit none
  
  character(len=20 ) :: option
  character(len=512) :: arg1
  character(len=512) :: arg2
  character(len=512) :: arg3
  character(len=512) :: arg4
  character(len=512) :: arg5
  character(len=512) :: arg6
  character(len=512) :: arg7
  character(len=512) :: arg8
  
  character(len=512) :: orig_tile
  character(len=512) :: new_tile
  character(len=512) :: domain_def_file
  character(len=512) :: catch_def_file
  character(len=512) :: out_path 
  character(len=512) :: exp_id 
  character(len=512) :: orig_catch
  character(len=512) :: new_rtm
  character(len=512) :: orig_rtm
  character(len=512) :: new_catch
  character(len=512) :: orig_BC
  character(len=512) :: new_BC
  character(len=512) :: orig_Veg
  character(len=512) :: new_veg
  character(len=512) :: orig_ease
  character(len=512) :: new_ease
  character(len=512) :: f2g_file
  character(len=12 ) :: ymdhm
  character(len=12 ) :: SURFLAY
  
  call get_command_argument(1,option)
  call get_command_argument(2,arg1)
  call get_command_argument(3,arg2)
  call get_command_argument(4,arg3)
  call get_command_argument(5,arg4)
  call get_command_argument(6,arg5)
  call get_command_argument(7,arg6)
  call get_command_argument(8,arg7)
  call get_command_argument(9,arg8)
  
  if( trim(option) == "c_f2g") then

     ! (1) generate 'f2g.txt'
     ! (2) generate tile.domain if it is local

     orig_tile       = arg1
     domain_def_file = arg2
     out_path        = arg3
     catch_def_file  = arg4
     exp_id          = arg5
     ymdhm           = trim(adjustl(arg6))
     SURFLAY         = trim(adjustl(arg7))
     f2g_file        = arg8

     call createf2g(orig_tile,domain_def_file,trim(out_path),catch_def_file,trim(exp_id),ymdhm, SURFLAY, f2g_file)
     
  else if (trim(option) == "c_localtile") then
     
     orig_tile  = arg1
     new_tile   = arg2
     f2g_file   = arg3
     call createLocalTilefile(f2g_file, orig_tile,new_tile)
     
  else if (trim(option) == "c_localbc" ) then
     
     orig_BC  = arg1
     new_BC   = arg2
     f2g_file = arg3

     call createLocalBC(f2g_file, orig_BC, new_BC)
     
  else if (trim(option) == "c_localvegrst") then

     orig_veg = arg1
     new_veg  = arg2
     f2g_file = arg3

     call  createLocalVegRestart(f2g_file, orig_veg, new_veg)      

  else if (trim(option) == "c_localmwrtmrst") then

     orig_rtm = arg1
     new_rtm  = arg2
     f2g_file = arg3

     call  createLocalmwRTMRestart(f2g_file, orig_rtm, new_rtm)      

  else if (trim(option) == "c_localcatchrst") then

     orig_catch = arg1
     new_catch  = arg2
     f2g_file   = arg3

     call createLocalCatchRestart(f2g_file, orig_catch, new_catch)

  else if (trim(option)=="correctease") then

     orig_ease = arg1 
     new_ease  = arg2

     call correctEase(orig_ease,new_ease) 

  else if (trim(option)=="c_convert_pert") then

     out_path  = arg3
     exp_id    = arg4

     call convert_pert_rst(arg1,arg2, out_path,exp_id)
     
  else if (trim(option) == "optimize") then
     
      
     call optimize_latlon(arg1,arg2, arg3, arg4)
     
  else
     
     print*, " wrong preprocess option:",option
     
  end if
  
end program main

! ====================== EOF =======================================================
