! how to use :
! ./preprocess_ldas option arg1 arg2 arg3

module preprocess_module
  
  use netcdf
  
  use MAPL
  
  use LDAS_TileCoordType,              ONLY:   &
       tile_coord_type,                        &
       grid_def_type,                          &
       io_grid_def_type,                       &
       tile_typ_land
  
  use nr_ran2_gasdev,                  ONLY:   &
       NRANDSEED
  
  use LDAS_DateTimeMod,                ONLY:   &
       date_time_type
  
  use force_and_cat_progn_pert_types,  ONLY:   &
       N_progn_pert_max,                       &
       N_force_pert_max
  
  use catch_types,                     ONLY:   &
       cat_param_type,                         &
       N_gt
  
  use preprocess_ldas_subs,            ONLY:   &
       LDAS_read_til_file
  
  use LDAS_ensdrv_functions,           ONLY:   &
       get_io_filename
  
  use LDAS_ensdrv_init_routines,       ONLY:   &
       domain_setup,                           &
       read_cat_param
  
  use LDAS_ensdrv_init_routines,       ONLY:   &
       io_domain_files
  
  use gFTL_StringVector
  
  use pFIO
  
end module preprocess_module

program main
   use preprocess_module
   implicit none
   character(len=20) :: option
   character(len=200) :: arg1
   character(len=200) :: arg2
   character(len=200) :: arg3
   character(len=200) :: arg4
   character(len=200) :: arg5
   character(len=200) :: arg6
   character(len=200) :: arg7

   character(len=200) :: orig_tile
   character(len=200) :: new_tile
   character(len=200) :: domain_def_file
   character(len=200) :: catch_def_file
   character(len=200) :: out_path 
   character(len=200) :: exp_id 
   character(len=200) :: orig_catch
   character(len=200) :: new_rtm
   character(len=200) :: orig_rtm
   character(len=200) :: new_catch
   character(len=200) :: orig_BC
   character(len=200) :: new_BC
   character(len=200) :: orig_Veg
   character(len=200) :: new_veg
   character(len=200) :: orig_ease
   character(len=200) :: new_ease
   character(len=12)  :: ymdhm
   character(len=12)  :: SURFLAY

   call get_command_argument(1,option)
   call get_command_argument(2,arg1)
   call get_command_argument(3,arg2)
   call get_command_argument(4,arg3)
   call get_command_argument(5,arg4)
   call get_command_argument(6,arg5)
   call get_command_argument(7,arg6)
   call get_command_argument(8,arg7)

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
      call createf2g(orig_tile,domain_def_file,trim(out_path),catch_def_file,trim(exp_id),ymdhm, SURFLAY)

   else if (trim(option) == "c_localtile") then

      orig_tile  = arg1
      new_tile   = arg2
      call createLocalTilefile(orig_tile,new_tile)

   else if (trim(option) == "c_localbc" ) then

      orig_BC = arg1
      new_BC =  arg2
      call createLocalBC(orig_BC, new_BC)

   else if (trim(option) == "c_localvegrst") then
      orig_veg = arg1
      new_veg  = arg2
      call  createLocalVegRestart(orig_veg, new_veg)      
   else if (trim(option) == "c_localmwrtmrst") then
      orig_rtm = arg1
      new_rtm  = arg2
      call  createLocalmwRTMRestart(orig_rtm, new_rtm)      
   else if (trim(option) == "c_localcatchrst") then
      orig_catch = arg1
      new_catch  = arg2
      call createLocalCatchRestart(orig_catch, new_catch)
   else if (trim(option)=="correctease") then
       orig_ease = arg1 
       new_ease  = arg2
       call correctEase(orig_ease,new_ease) 
   else if (trim(option)=="c_convert_pert") then
       out_path  = arg3
       exp_id    = arg4
       call convert_pert_rst(arg1,arg2, out_path,exp_id)
   else if (trim(option) == "optimize") then

       call optimize_latlon(arg1,arg2)
   else
       print*, " wrong preprocess option:",option
   end if

   contains
      function upcase(string) result(upper)
         character(len=*), intent(in) :: string
         character(len=len(string)) :: upper
         integer :: j
         do j = 1,len(string)
            if(string(j:j) >= "a" .and. string(j:j) <= "z") then
               upper(j:j) = achar(iachar(string(j:j)) - 32)
            else
               upper(j:j) = string(j:j)
            end if
         end do
      end function upcase
end program

subroutine createf2g(orig_tile,domain_def,out_path,catch_def_file,exp_id,ymdhm, SURFLAY)
   use preprocess_module
   implicit none
   character(*) :: orig_tile
   character(*) :: domain_def
   character(*) :: out_path
   character(*) :: catch_def_file
   character(*) :: exp_id
   character(*) :: ymdhm
   character(*) :: SURFLAY

   real :: minlon,maxlon,minlat,maxlat
   character(len=200):: exclude_file,include_file
   character(len=300):: bcs_path
   logical :: file_exist
   logical :: d_exist,c_exist

   integer :: n

   type(grid_def_type) :: tile_grid_g,tile_grid_d
   type(tile_coord_type), dimension(:), pointer :: tile_coord_g => null()
   type(tile_coord_type), dimension(:), pointer :: tile_coord_d => null()
   integer, dimension(:), pointer     :: f2g => null()
   integer, dimension(:), pointer     :: d2g => null()
   integer, dimension(:), pointer     :: d2f => null()
   integer :: N_catg, N_catd,n1,n2,N_catf

   type(cat_param_type),  dimension(:), allocatable :: cp
   real :: dzsf

   namelist / domain_inputs /                              &
         minlon, maxlon,minlat,maxlat,                     &
         exclude_file,include_file

   inquire(file=trim(orig_tile),exist=file_exist)
   if( .not. file_exist) stop ("original tile file not exist")

   inquire(file=trim(domain_def),exist=d_exist)
   if( .not. d_exist) then
       print*,"no domain definition file"
   endif
       
   inquire(file=trim(catch_def_file),exist=c_exist)
   if( .not. c_exist) then
      print*,"no catchment definition file:" , catch_def_file
   endif 


   if(d_exist) then
       open (10, file=trim(domain_def), delim='apostrophe', action='read', status='old')
       read (10, nml= domain_inputs)
       close(10)
   else
       minlon = -180.
       maxlon = 180.
       minlat = -90.
       maxlat = 90.
       exclude_file = ' '
       include_file = ' '
   endif

   call LDAS_read_til_file(orig_tile,catch_def_file,tile_grid_g,tile_coord_g,f2g)

   N_catg=size(tile_coord_g)

   ! include and exclude files are absolute

   call domain_setup(                                               &
       N_catg, tile_coord_g,                                        &
       tile_grid_g,                                                 &
      ' ', exclude_file, ' ', include_file,                         &
       trim(out_path), 'exp_domain ', trim(exp_id),                 &
       minlon, minlat, maxlon, maxlat,                              &
       N_catd, d2g, tile_coord_d,                                   &
       tile_grid_d )

   allocate(cp(N_catd))

   read(SURFLAY,*) dzsf
   print*, "SURFLAY: ", dzsf
   n1 = index(catch_def_file,'/clsm/')
   bcs_path(1:n1-1) = catch_def_file(1:n1-1)
   call read_cat_param( N_catg, N_catd, d2g, tile_coord_d, dzsf, bcs_path(1:n1-1), bcs_path(1:n1-1),bcs_path(1:n1-1),  &
        cp )
   call write_cat_param(cp,N_catd)

   allocate(d2f(N_catd))
   d2f = 0
   N_catf = size(f2g)
   if( N_catf /= N_catg) then
      n = 1
      do n1 = 1,N_catd
         do n2 = n, N_catf
            if (d2g(n1) == f2g(n2)) then
               d2f(n1) = n2
               n = n2+1
               exit
            endif
         enddo
      enddo
      if(any(d2f == 0)) stop " Domain includes those excluded tiles"
      print*," f2g now is d2f "
   else
     d2f = d2g
   endif
   open(40,file='f2g.txt',form='formatted',action='write')
   write(40,*)N_catf
   write(40,*)N_catd
   do n=1,N_catd
       write(40,*)d2f(n)
   enddo
   do n=1,N_catd
       write(40,*)d2g(n)
   enddo
   close(40)
   if (associated(f2g)) deallocate(f2g)
   if (associated(d2g)) deallocate(d2g)
   if (associated(d2f)) deallocate(d2f)

contains

   subroutine write_cat_param(cat_param, N_catd)
    type(cat_param_type), intent(in) :: cat_param(:)
    integer,intent(in) :: N_catd
    character(len=300):: fname
    type(date_time_type) :: start_time
    
    integer :: k,n

    read(ymdhm(1:4),*) start_time%year  ! 4-digit year
    read(ymdhm(5:6),*) start_time%month ! month in year
    read(ymdhm(7:8),*) start_time%day   ! day in month
    read(ymdhm(9:10),*) start_time%hour  ! hour of day
    read(ymdhm(11:12),*) start_time%min  
    start_time%sec = 0
    start_time%pentad  = -9999        ! pentad of year
    start_time%dofyr   = -9999   

    fname = get_io_filename(trim(out_path), trim(exp_id),'ldas_catparam', date_time=start_time, &
            dir_name='rc_out', file_ext='.bin')

    open(10, file=trim(fname), form='unformatted', status='unknown', action='write')

    print*, 'Writing catparam file : ' // trim(fname)

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

    write (10) (real(cat_param(n)%vegcls),     n=1,N_catd)
    write (10) (real(cat_param(n)%soilcls30),  n=1,N_catd)
    write (10) (real(cat_param(n)%soilcls100), n=1,N_catd)

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

   end subroutine write_cat_param

end subroutine createf2g

subroutine readsize(N_catg,N_catf)
   use preprocess_module
   implicit none
   integer,intent(out) :: N_catg
   integer,intent(out) :: N_catf

   logical :: file_exist

   inquire(file=trim('f2g.txt'),exist=file_exist)
   if(file_exist) then
      open(40,file='f2g.txt',form='formatted',action='read',status='old')
      read(40,*)N_catg
      read(40,*)N_catf
      close(40)
   else
      print*, " wrong, no f2g.txt"
   endif
end subroutine readsize

subroutine readf2g(N_catf,f2g)
   use preprocess_module
   implicit none
   integer,intent(in) :: N_catf
   integer,dimension(N_catf),intent(inout) :: f2g

   integer :: N_catg
   logical :: file_exist
   integer :: local_size,n

   inquire(file=trim('f2g.txt'),exist=file_exist)
   if(file_exist) then
      open(40,file='f2g.txt',form='formatted',action='read',status='old')
      read(40,*)N_catg
      read(40,*)local_size

      if(local_size /= N_catf) print*, "wrong f2g.txt"

      if(N_catg == N_catf) then
         close(40)
         return
      endif

      do n=1,N_catf
         read(40,*)f2g(n)
      enddo
      close(40)
    ! call MAPL_sort(this%f2g)
   else
      print*, " wrong, no f2g.txt"
   endif
   
end subroutine readf2g

subroutine createLocalTilefile(orig_tile,new_tile)
   use preprocess_module
   implicit none
   character(*) :: orig_tile
   character(*) :: new_tile
   character(len=200):: line

   logical :: file_exist

   integer, dimension(:),allocatable :: f2g 
   integer :: N_catg, N_catf,n,stat, ty
   integer :: N_tile,N_grid,g_id

   inquire(file=trim(orig_tile),exist=file_exist)
   if( .not. file_exist) stop ("original tile file not exist")

   ! Set default local tile file name
   call readsize(N_catg,N_catf)
   if(N_catg == N_catf) then
      print*, "It is global domain..."
      return
   endif
   allocate(f2g(N_catf))
   call readf2g(N_catf,f2g)   

   open(40,file=trim(orig_tile),action="read")
   open(50,file=trim(new_tile),action="write")

   ! put the head back for tile file
   do n=1,5
       read(40,'(A)') line
       if(n==1) then
          read(line,*) N_tile
       endif
       if(n==2) then
          read(line,*) N_grid
       endif
       write(50,'(A)') trim(line)
   enddo
   if (N_grid ==2) then
       do n=1,3
          read(40,'(A)') line
          write(50,'(A)') trim(line)
       enddo
   endif

   g_id = 0
   do while(.true.)
       read(40,'(A)',IOSTAT=stat) line
       if(IS_IOSTAT_END(stat)) exit
       ! just read the first four
       read(line,*) ty
       if( ty == tile_typ_land ) then
           n=index(line,'100')
       ! here g_id is the global land tiles
           g_id=g_id+1
           if(.not. any( f2g(:) == g_id)) then
       ! add 1000 to it so it will be excluded
               line(n-1:n+2)='1100'
           endif
       endif
       write(50,'(A)') trim(line)
   enddo
   close(40)
   close(50)

end subroutine createLocalTilefile

subroutine createLocalBC(orig_BC, new_BC)
   use preprocess_module
   implicit none
   character(*),intent(in) :: orig_BC
   character(*),intent(in) :: new_BC

   real,dimension(14) :: tmprealvec14
   real,allocatable ::   tmpvec(:)  
   integer :: istat, N_catg,N_catf
   integer,dimension(:),allocatable :: f2g
  
   call readsize(N_catg,N_catf)
   if(N_catg==N_catf) return
   allocate(f2g(N_catf))
   call readf2g(N_catf,f2g)  

   allocate(tmpvec(N_catg))
   open(10,file=trim(orig_BC),form='unformatted',action='read',status='old',iostat=istat)
   open(20,file=trim(new_BC),form='unformatted',action='write')

   do while(.true.)
      read(10,iostat=istat) tmprealvec14
      if(IS_IOSTAT_END(istat)) exit
      read(10) tmpvec
      write(20) tmprealvec14
      write(20) tmpvec(f2g)
   enddo
   close(10)
   close(20)
   deallocate(tmpvec)
end subroutine createLocalBC

subroutine createLocalCatchRestart(orig_catch, new_catch)
   use preprocess_module
   implicit none
   character(*),intent(in):: orig_catch
   character(*),intent(in):: new_catch
   integer,parameter :: subtile=4
   integer :: istat, filetype, rc,i, j, ndims
   real,allocatable :: tmp1(:)
   real,allocatable :: tmp2(:,:)
   type(Netcdf4_FileFormatter) :: InFmt,OutFmt
   type(FileMetadata)        :: OutCfg
   type(FileMetadata)        :: InCfg
   integer                   :: dim1,dim2
   type(StringVariableMap), pointer :: variables
   type(Variable), pointer :: var
   type(StringVariableMapIterator) :: var_iter
   type(StringVector), pointer :: var_dimensions
   character(len=:), pointer :: vname,dname
   integer ::n, N_catg,N_catf
   integer,dimension(:),allocatable :: f2g

   call readsize(N_catg,N_catf)
   if(N_catg == N_catf) return
   allocate(f2g(N_catf))
   call readf2g(N_catf,f2g) 

   allocate(tmp1(N_catg))
   allocate(tmp2(N_catg,subtile))

     ! check file type

   call MAPL_NCIOGetFileType(orig_catch, filetype,rc=rc)

   if (filetype /= 0) then

       print*, "Catchment restart is binary"

         ! binary 

      open(10,file=trim(orig_catch),form='unformatted',action='read',status='old',iostat=istat)
      open(20,file=trim(new_catch),form='unformatted',action='write')
         
      do n=1,30
         read(10) tmp1
         write(20) tmp1(f2g)
      enddo

      do n=1,2
         read(10) tmp2
         write(20) tmp2(f2g,:)
      enddo

      do n=1,20
         read(10) tmp1
         write(20) tmp1(f2g)
      enddo
         ! note : the offline restart does not have the last five variables
      do n=1,4
         read(10,iostat=istat) tmp2
         if(.not. IS_IOSTAT_END(istat)) write(20) tmp2(f2g,:)
      enddo
         ! 57 WW
      read(10,iostat=istat) tmp2
      if(.not. IS_IOSTAT_END(istat)) write(20) tmp2(f2g,:)
         
      close(10)
      close(20)
   else
         
         ! filetype = 0 : nc4 output file will also be nc4
      
      call InFmt%open(trim(orig_catch), pFIO_READ,rc=rc)
      InCfg  = InFmt%read(rc=rc)
      OutCfg = InCfg
      
      call OutCfg%modify_dimension('tile', size(f2g), rc=rc)

      call OutFmt%create(trim(new_catch),rc=rc)
      call OutFmt%write(OutCfg,rc=rc)

      variables => InCfg%get_variables()
      var_iter = variables%begin()
      do while (var_iter /= variables%end())

         vname => var_iter%key()
         var => var_iter%value()
         var_dimensions => var%get_dimensions()

         ndims = var_dimensions%size()

         if (trim(vname) =='time') then
             call var_iter%next()
             cycle
         endif

         if (ndims == 1) then
            call MAPL_VarRead (InFmt,vname,tmp1)
            call MAPL_VarWrite(OutFmt,vname,tmp1(f2g))
         else if (ndims == 2) then

            dname => var%get_ith_dimension(2)
            dim1=InCfg%get_dimension(dname)
            do j=1,dim1
               call MAPL_VarRead ( InFmt,vname,tmp1 ,offset1=j)
               call MAPL_VarWrite(OutFmt,vname,tmp1(f2g),offset1=j)
            enddo

         else if (ndims == 3) then

            dname => var%get_ith_dimension(2)
            dim1=InCfg%get_dimension(dname)
            dname => var%get_ith_dimension(3)
            dim2=InCfg%get_dimension(dname)
            do i=1,dim2
               do j=1,dim1
                 call MAPL_VarRead ( InFmt,vname,tmp1 ,offset1=j,offset2=i)
                 call MAPL_VarWrite(OutFmt,vname,tmp1(f2g) ,offset1=j,offset2=i)
              enddo
           enddo

         end if
         call var_iter%next()
      enddo
      call inFmt%close(rc=rc)
      call OutFmt%close(rc=rc)
   end if ! file type nc4
   print*, "done create local catchment restart"
end subroutine createLocalCatchRestart

subroutine createLocalmwRTMRestart(orig_mwrtm, new_mwrtm)
   use preprocess_module
   implicit none
   character(*),intent(in):: orig_mwrtm
   character(*),intent(in):: new_mwrtm
   integer,parameter :: subtile=4
   integer :: rc
   real,allocatable :: tmp1(:)
   type(Netcdf4_FileFormatter) :: InFmt,OutFmt
   type(FileMetadata)        :: OutCfg
   type(FileMetadata)        :: InCfg
   
   type(StringVariableMap), pointer :: variables
   type(StringVariableMapIterator) :: var_iter
   character(len=:), pointer :: vname
   integer :: N_catg,N_catf
   integer,dimension(:),allocatable :: f2g

   call readsize(N_catg,N_catf)
   if(N_catg == N_catf) return
   allocate(f2g(N_catf))
   call readf2g(N_catf,f2g) 

   allocate(tmp1(N_catg))
         
   ! nc4 in and out file will also be nc4
   call InFmt%open(trim(orig_mwrtm), pFIO_READ,rc=rc)
   InCfg = InFmt%read(rc=rc)
   OutCfg = InCfg

   call OutCfg%modify_dimension('tile', size(f2g), rc=rc)

   call OutFmt%create(trim(new_mwrtm),rc=rc)
   call OutFmt%write(OutCfg,rc=rc)

   variables => InCfg%get_variables()
   var_iter = variables%begin()
   do while (var_iter /= variables%end())
      vname => var_iter%key()
      call MAPL_VarRead (InFmt,vname,tmp1)
      call MAPL_VarWrite(OutFmt,vname,tmp1(f2g))
      call var_iter%next()
   enddo

   call inFmt%close(rc=rc)
   call OutFmt%close(rc=rc)

   deallocate(f2g,tmp1)
         
end subroutine createLocalmwRTMRestart

subroutine createLocalVegRestart(orig_veg, new_veg)
   use preprocess_module
   implicit none
   character(*),intent(in):: orig_veg
   character(*),intent(in):: new_veg
   integer :: istat
   real,allocatable :: rity(:)
   real,allocatable :: z2(:)
   real,allocatable :: ascatz0(:)
   real,allocatable :: tmp(:)

   integer :: N_catg,N_catf
   integer,dimension(:),allocatable :: f2g
   integer :: filetype
   type(Netcdf4_FileFormatter) :: InFmt,OutFmt
   type(FileMetadata)        :: OutCfg
   type(FileMetadata)        :: InCfg

   type(StringVariableMap), pointer :: variables
   type(StringVariableMapIterator) :: var_iter
   character(len=:), pointer :: vname
   integer :: rc

   call readsize(N_catg,N_catf)
   if(N_catg == N_catf) return
   allocate(f2g(N_catf))
   call readf2g(N_catf,f2g)  
    
   allocate(rity(N_catg))
   allocate(z2(N_catg))
   allocate(ascatz0(N_catg))

   call MAPL_NCIOGetFileType(orig_veg, filetype,rc=rc)

   if (filetype /=0) then
      open(10,file=trim(orig_veg),form='unformatted',action='read',status='old',iostat=istat)
      open(20,file=trim(new_veg),form='unformatted',action='write')
      read(10) rity 
      read(10) z2 
      read(10) ascatz0 
      write(20) rity(f2g)
      write(20) z2(f2g) 
      write(20) ascatz0(f2g) 

      close(10)
      close(20)
   else
    ! nc4 in and out file will also be nc4
      call InFmt%open(trim(orig_veg), pFIO_READ,rc=rc)
      InCfg = InFmt%read(rc=rc)
      OutCfg = InCfg

      call OutCfg%modify_dimension('tile', size(f2g), rc=rc)

      call OutFmt%create(trim(new_veg),rc=rc)
      call OutFmt%write(OutCfg,rc=rc)

      variables => InCfg%get_variables()
      var_iter = variables%begin()
      allocate(tmp(N_catg))
      do while (var_iter /= variables%end())
         vname => var_iter%key()
         call MAPL_VarRead (InFmt,vname,tmp)
         call MAPL_VarWrite(OutFmt,vname,tmp(f2g))
         call var_iter%next()
      enddo

      call inFmt%close(rc=rc)
      call OutFmt%close(rc=rc)
      deallocate(tmp)
   endif
   deallocate(f2g)

end subroutine createLocalVegRestart

subroutine correctEase(orig_ease,new_ease)
   implicit none
   character(*),intent(in) :: orig_ease
   character(*),intent(in) :: new_ease
   logical :: file_exist,is_oldEASE
   integer :: i, N_tile, N_grid
   character(len=200) :: tmpline

   inquire(file=trim(orig_ease),exist=file_exist)
   if( .not. file_exist) stop (" no ease_tile_file")

   open(55,file=trim(orig_ease),action='read')
   read(55,*) N_tile
   read(55,*) N_grid
   read(55,*)
   read(55,*)
   read(55,*)
   read(55,'(A)') tmpline
   close(55)

   is_oldEASE= .false.
   if(N_grid ==1 .and. index(tmpline,'OCEAN')/=0) is_oldEASE=.true.

   if( is_oldEASE) then
      open(55,file=trim(orig_ease),action='read')
      open(56,file=trim(new_ease),action='write')
      do i =1,5
         read(55,'(A)')tmpline
         write(56,'(A)')trim(tmpline)
      enddo
      read(55,*)
      read(55,*)
      read(55,*)
      do i=1,N_tile
         read(55,'(A)')tmpline
         write(56,'(A)')trim(tmpline)
      enddo
      close(56)
      close(55)
   end if
end subroutine correctEase

! This program optimized the domain setup by producing 
!
! NX: N_proc
! NY: 1
! IMS: i0 i1 12 .....
! run this program with:
! ./a.out tile_file N_proc

subroutine optimize_latlon(fname,arg)
   implicit none
   character(*) :: fname,arg
   integer :: N_proc
   integer :: N_tile,N_lon,N_lat,N_grid
   integer,allocatable :: landPosition(:)
   integer,allocatable :: IMS(:),JMS(:)
   integer,allocatable :: local_land(:)
   integer :: total_land
   integer :: n,typ,tmpint
   real ::  tmpreal
   integer :: avg_land,n0,local
   integer :: i,s,e,j,k,n1,n2
   logical :: file_exist
   character(len=100):: tmpLine
   character(len=100):: gridname
   real :: rate,rates(60),maxf(60)
   integer :: IMGLOB, JMGLOB
   integer :: face(6),face_land(6)

   inquire(file=trim(fname),exist=file_exist)
   if( .not. file_exist) stop ( "tile file not exist")
   read (arg,*) N_proc
  

   open (10, file=trim(fname), form='formatted', action='read')
   read (10,*) N_tile
   read (10,*) N_grid         ! some number (?)
   read (10,*) gridname       ! some string describing tile definition grid (?)
   read (10,*) N_lon
   read (10,*) n_lat

   if (index(gridname,"CF") /=0) then

      IMGLOB = N_lon
      JMGLOB = N_lat
      if(JMGLOB/6 /= IMGLOB) stop " wrong im, jm"

      allocate(landPosition(JMGLOB))
      landPosition = 0
      total_land   = 0

      if(N_grid == 2) then
         read (10,*)          ! some string describing ocean grid                   (?)
         read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
         read (10,*)
      endif
      
      do n = 1,N_tile
         read (10,*)  &
           typ,           &   !  1
           tmpreal,       &   !  2  *
           tmpreal,       &   !  3
           tmpreal,       &   !  4
           i ,            &   !  5
           j ! ,        &   !  6
           !tmpreal,       &   !  7
           !tmpint,        &   !  8
           !tmpreal,       &   !  9  *
           !tmpint,        &   ! 10
           !tmpreal,       &   ! 11
           !tmpint       ! 12  * (previously "tile_id")
         if(typ==100) then
           total_land=total_land+1
           landPosition(j) = landPosition(j)+1
         endif
         ! assume all lands are at the beginning 
         if(typ /= 100 .and. typ /= 1100 ) exit

      enddo
      close(10)

      if(mod(N_proc,6) /=0) then
        print*,"WARNING: ntasks should be adjusted to multiple of 6 for cubed-sphere grid :",N_proc
        N_proc = N_proc-mod(N_proc,6)
      endif

      print*, "total tiles", total_land

      if(sum(landPosition) /= total_land) print*, "wrong counting of land"

      do k=1,6
        n1 = (k-1)*IMGLOB+1
        n2 = k*IMGLOB
        face_land(k) = sum(landPosition(n1:n2)) 
        face(k) = nint(1.0*face_land(k)/total_land * N_proc)
        if ( face(k) == 0) face(k) = 1
      enddo

      ! now make sure sum(face) == N_proc
      k=sum(face)-N_proc

      if (k < 0) then
        do i=1, -k
           n=minloc(face,DIM=1)
           face(n) = face(n)+1
        enddo
      else
        do i = 1,k
           n=maxloc(face,DIM=1)
           face(n) = face(n)-1
        enddo
      endif

      if (sum(face) /= N_proc) stop " wrong proc face"

   ! 2) each process should have average land tiles

      ALLOCATE(JMS(N_proc))
      allocate(local_land(N_Proc))
      JMS = 0
      local_land = 0

      local  = 0
      n0     = 0
      j = 0
      do k=1,6
        n1 = (k-1)*IMGLOB+1
        n2 = k*IMGLOB

        do i=1,60
           rates(i) = (i-1)*0.1
        enddo

        maxf=rms_cs(rates)
        i=minloc(maxf,DIM=1)
        rate = rates(i)
        avg_land = ceiling(1.0*face_land(k)/face(k))
        avg_land = avg_land - nint(rate*avg_land/face(k))

        tmpint = 0
        j = j+face(k) 

        do n = n1,n2
           tmpint=tmpint+landPosition(n)
           if((local+1) == j .and. n < n2) cycle
           if((tmpint .ge. avg_land) .or. (n==n2)) then
              local = local + 1
              local_land(local)=tmpint
              JMS(local)=n-n0
              tmpint=0
              n0=n
           endif
        enddo
        local = j
      enddo

      ! adjust JMS.rc make make sure no process has 0 or 1
      j = 1
      do k = 1,6
         n1 = j
         n2 = j+face(k)-1
         do i = n1,n2
            if(JMS(i) == 0) then
               n = maxloc(JMS(n1:n2),DIM=1)
               JMS(i) = 1
               JMS(n+n1-1) = JMS(n+n1-1)-1
            endif
            if(JMS(i) == 1) then
               n = maxloc(JMS(n1:n2),DIM=1)
               JMS(i) = 2
               JMS(n+n1-1) = JMS(n+n1-1)-1
            endif
         enddo
         j=j+face(k)
      enddo   

      print*,"land_distribute: ",local_land
      print*, "JMS.rc", JMS
      if( sum(JMS) /= JMGLOB) then
        print*, sum(JMS), JMGLOB
        stop ("wrong cs-domain distribution")
      endif
      tmpint = 0
      k = 0
      do n = 1, N_proc
        tmpint= tmpint+JMS(n)
        if( tmpint == IMGLOB) then
          k=k+1
          tmpint = 0
        endif
      enddo 
      
      if( k /=6 ) stop ("one or more processes may accross the face")

      open(10,file="optimized_distribution",action='write')
      write(10,'(A)') "GEOSldas.GRIDNAME:  " // trim(gridname)
      write(10,'(A)') "GEOSldas.GRID_TYPE:  Cubed-Sphere"
      write(10,'(A)') "GEOSldas.NF:  6"
      write(10,'(A,I6)') "GEOSldas.IM_WORLD: ", IMGLOB
      write(10,'(A)') "GEOSldas.LM:   1"
      write(10,'(A,I5)') "NY: ",N_proc
      write(10,'(A)') "NX:   1"
      write(10,'(A)') "GEOSldas.JMS_FILE:    JMS.rc"
      close(10)
       
      open(10,file="JMS.rc",action='write')
      write(10,'(I5,I5)') N_proc, maxval(face)
      do n=1,N_proc
        write(10,'(I8)') JMS(n)
      enddo
      close(10)

   else

      allocate(IMS(N_Proc))
      allocate(local_land(N_Proc))
      IMS=0
      local_land = 0
   
      if(N_grid == 2) then
         read (10,*)          ! some string describing ocean grid                   (?)
         read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
         read (10,*)   
         read(10,'(A)') tmpLine
      else
         read(10,'(A)') tmpLine
         if (index(tmpLine,"OCEAN") /=0) then
            read (10,*)          ! some string describing ocean grid                   (?)
            read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
            read (10,*)   
            read(10,'(A)') tmpLine
         endif
      endif
  
      if (index(gridname,'EASE') /=0) then
         s=0
         e=N_lon-1
      else
         s=1
         e=N_lon
      endif
      allocate(landPosition(s:e))

      landPosition=0
      total_land= 0

   ! 1) read through tile file, put the land tile into the N_lon of bucket

      read (tmpLine,*)  &
        typ,           &   !  1
        tmpreal,       &   !  2  *
        tmpreal,       &   !  3
        tmpreal,       &   !  4
        i !,             &   !  5
      if(typ==100) then
         total_land=total_land+1
         landPosition(i) = landPosition(i)+1
      endif

      do n = 2,N_tile
         read (10,*)  &
           typ,           &   !  1
           tmpreal,       &   !  2  *
           tmpreal,       &   !  3
           tmpreal,       &   !  4
           i !,             &   !  5
           !tmpint,        &   !  6
           !tmpreal,       &   !  7
           !tmpint,        &   !  8
           !tmpreal,       &   !  9  *
           !tmpint,        &   ! 10
           !tmpreal,       &   ! 11
           !tmpint       ! 12  * (previously "tile_id")
         if(typ==100) then
           total_land=total_land+1
           landPosition(i) = landPosition(i)+1
         endif
         ! assume all lands are at the beginning 
         if(typ /= 100 .and. typ /=1100 ) exit
      enddo

      close(10)

      if(sum(landPosition) /= total_land) print*, "wrong counting of land"

      do n=1,60
        rates(n) = (n-1)*0.15
      enddo

      maxf=rms(rates)
      n=minloc(maxf,DIM=1)
      rate = rates(n)

   ! 2) each process should have average land tiles

      avg_land = ceiling(1.0*total_land/N_proc)
      print*,"avg_land",avg_land
      
      ! rate is used to readjust the avg_land
      ! in case that the last processors don't have any land tiles,
      ! we can increase ther rates
      
      avg_land = avg_land - nint(rate*avg_land/N_proc)
      print*,"re adjust the avg_land",avg_land
      tmpint = 0
      local = 1
      n0 = s-1
      do n=s,e
        tmpint=tmpint+landPosition(n)
        if(local == N_proc .and. n < e) cycle ! all lefteover goes to the last process
        if((tmpint .ge. avg_land) .or. (n==e)) then
          local_land(local)=tmpint
          IMS(local)=n-n0
          tmpint=0
          n0=n
          local = local + 1
        endif
      enddo
      print*,"rms rate: ", rms(rate)

      print*,"land_distribute: ",local_land

      if( sum(local_land) /= total_land) stop ("wrong distribution")
      if( sum(IMS) /= N_lon) stop ("wrong domain distribution")
 
      open(10,file="optimized_distribution",action='write')
      write(10,'(A)') "GEOSldas.GRID_TYPE:  LatLon"
      write(10,'(A)') "GEOSldas.GRIDNAME:  "//trim(gridname)
      write(10,'(A)') "GEOSldas.LM:  1"
      write(10,'(A)') "GEOSldas.POLE:  PE"
      write(10,'(A)') "GEOSldas.DATELINE:  DE"
      write(10,'(A,I6)') "GEOSldas.IM_WORLD: ",N_lon
      write(10,'(A,I6)') "GEOSldas.JM_WORLD: ",N_lat

      write(10,'(A,I5)') "NX: ",N_proc
      write(10,'(A)') "NY:   1"

      write(10,'(A)') "GEOSldas.IMS_FILE:    IMS.rc"
      close(10)
       
      open(10,file="IMS.rc",action='write')
      write(10,'(I5)') N_proc
      do n=1,N_proc
        write(10,'(I8)') IMS(n)
      enddo
      close(10)

   endif

   contains 

      elemental function rms(rates) result (f)
         real :: f
         real,intent(in) :: rates
         integer :: tmpint,local
         integer :: n0,proc,n
         integer :: avg_land
         integer,allocatable :: local_land(:)

         allocate (local_land(N_proc))
         local_land = 0
         avg_land = ceiling(1.0*total_land/N_proc)
         avg_land = avg_land -nint(rates*avg_land/N_proc)

         tmpint = 0
         local = 1
         n0 = s-1
         do n=s,e
           tmpint=tmpint+landPosition(n)
           if(local == N_proc .and. n < e) cycle ! all lefteover goes to the last process
           if((tmpint .ge. avg_land) .or. (n==e)) then
             local_land(local)=tmpint
             tmpint=0
             n0=n
             local = local + 1
           endif
         enddo
         f = 0.0
         do proc = 1, N_proc
            f =max(f,1.0*abs(local_land(proc)-avg_land))
         enddo
         deallocate(local_land)
      end function

      elemental function rms_cs(rates) result (f)
         real :: f
         real,intent(in) :: rates
         integer :: tmpint,local
         integer :: proc,n
         integer :: avg_land
         integer,allocatable :: local_land(:)
         integer :: n1,n2

         allocate (local_land(face(k)))
         local_land = 0
         avg_land = ceiling(1.0*face_land(k)/face(k))
         avg_land = avg_land -nint(rates*avg_land/face(k))
         if (avg_land <=0) then
             f = face_land(k)
             return
         endif

         tmpint = 0
         local = 1

         n1 = (k-1)*IMGLOB+1
         n2 = k*IMGLOB
         tmpint = 0
         do n = n1,n2
           tmpint=tmpint+landPosition(n)
           if(local == face(k) .and. n < n2) cycle ! all lefteover goes to the last process
           if((tmpint .ge. avg_land) .or. (n==n2)) then
              local_land(local)= tmpint
              tmpint=0
              local = local + 1
            endif
         enddo

         f = 0.0
         do proc = 1, face(k)
            ! punish for no land tiles
            f =max(f,1.0*abs(local_land(proc)-avg_land))
         enddo
         deallocate(local_land)
      end function

end subroutine optimize_latlon

subroutine convert_pert_rst(pfile_name,pfile_nc4,in_path,exp_id)
    use preprocess_module
    implicit none
    character(*),intent(in) :: pfile_name
    character(*),intent(in) :: in_path
    character(*),intent(in) :: exp_id
    character(*),intent(in) :: pfile_nc4

    integer :: N_catf,N_lon,N_lat,N_lonf,N_latf
    integer :: N_force_pert,N_progn_pert
    integer,pointer :: f2g(:)
   
    type(tile_coord_type), dimension(:), pointer :: tile_coord_f => null()

    type(grid_def_type) :: pert_grid_g
    type(grid_def_type) :: pert_grid_f 
    integer :: RC,istat
    integer,allocatable :: Pert_rseed(:)
    real,allocatable :: Force_pert_ntrmdt_f(:,:,:)
    real,allocatable :: Progn_pert_ntrmdt_f(:,:,:)
    
    call io_domain_files('r',in_path, trim(exp_id),N_catf,f2g,tile_coord_f,pert_grid_g,pert_grid_f,RC) 

    N_lon = pert_grid_g%N_lon
    N_lat = pert_grid_g%N_lat
    N_lonf= pert_grid_f%N_lon
    N_latf= pert_grid_f%N_lat

    call i_pert_ldas(RC)

    call o_pert_GEOSldas(rc)

   contains 

       subroutine i_pert_ldas(rc)
          integer,intent(inout),optional :: rc

          integer :: nrandseed_tmp
          type(grid_def_type) :: pert_grid_f_tmp 
          character(len=*), parameter :: Iam = 'io_pert_rstrt'
          integer :: k
          real, allocatable :: real_tmp(:)
 
          open(10, file=pfile_name, convert='big_endian',form='unformatted', status='old', &
          action='read', iostat=istat)

          ! one additional header line (as of 21 May 2010)!!!

          call io_grid_def_type( 'r', 10, pert_grid_f_tmp )

          read (10) nrandseed_tmp, N_force_pert, N_progn_pert

          ! check whether entries in file match passed arguments
          ! (check does *not* include *_pert_param!)

          if ( (nrandseed_tmp          /= NRANDSEED) ) then !          .or.               &
          !     (N_force_pert_tmp       /= N_force_pert)         .or.               &
          !     (N_progn_pert_tmp       /= N_progn_pert) ) then
               stop 'pert.rstrt file not compatible (1)'
          end if

          allocate(Pert_rseed(NRANDSEED))
          allocate(Force_pert_ntrmdt_f(N_lonf,N_latf, N_Force_pert)) 
          allocate(Progn_pert_ntrmdt_f(N_lonf,N_latf, N_Progn_pert))

          if ( index(pert_grid_f%gridtype,'LatLon')/=0  .or.         &
               index(pert_grid_f%gridtype,'LATLON')/=0  .or.         &
               index(pert_grid_f%gridtype,'latlon')/=0       ) then

             if ( (pert_grid_f_tmp%N_lon  /= pert_grid_f%N_lon)    .or.            &
                  (pert_grid_f_tmp%N_lat  /= pert_grid_f%N_lat)    .or.            &
                  (abs(pert_grid_f_tmp%ll_lon - pert_grid_f%ll_lon) > 1e-4) .or.   &
                  (abs(pert_grid_f_tmp%ll_lat - pert_grid_f%ll_lat) > 1e-4) .or.   &
                  (abs(pert_grid_f_tmp%dlon   - pert_grid_f%dlon)   > 1e-4) .or.   &
                  (abs(pert_grid_f_tmp%dlat   - pert_grid_f%dlat)   > 1e-4)      ) then
                stop 'pert.rstrt file not compatible (2)'
             end if

          else

             if ( index(pert_grid_f_tmp%gridtype,pert_grid_f%gridtype)==0  .or.   &
                  (pert_grid_f_tmp%N_lon  /= pert_grid_f%N_lon)            .or.   &
                  (pert_grid_f_tmp%N_lat  /= pert_grid_f%N_lat)            .or.   &
                  (pert_grid_f_tmp%i_offg /= pert_grid_f%i_offg)           .or.   &
                  (pert_grid_f_tmp%j_offg /= pert_grid_f%j_offg)                ) then
                stop 'pert.rstrt file not compatible (3)'
             end if

          end if

          ! reading
          read (10) Pert_rseed(:)
          allocate(real_tmp(N_lonf*N_latf))
          do k=1,N_force_pert
             !read (10) ((Force_pert_ntrmdt_f(i,j,k), i=1,N_lonf),j=1,N_latf)
             read (10) real_tmp(:)
             Force_pert_ntrmdt_f(:,:,k) = reshape(real_tmp,[N_lonf, N_latf])
          end do

          do k=1,N_progn_pert
             !read (10) ((Progn_pert_ntrmdt_f(i,j,k), i=1,N_lonf),j=1,N_latf)
             read (10) real_tmp(:)
             Progn_pert_ntrmdt_f(:,:,k) = reshape(real_tmp,[N_lonf, N_latf])
          end do
          
          close(10)
          deallocate(real_tmp)
          rc = 0
      end subroutine i_pert_ldas

      subroutine o_pert_GEOSldas(rc)
         integer,intent(inout) :: rc
         integer :: NCFOutID, STATUS
         integer :: seeddim,latdim, londim, Nforce,NProgn
         integer :: dims(3), seedid,forceid,prognid        
         integer :: xstart, ystart
         integer :: shuffle, deflate, deflate_level
         real    :: fill_value

         fill_value = -9999. !1.0e+15
         shuffle = 1
         deflate = 1
         deflate_level = 2

         !1) create file
         status = NF90_CREATE (trim(pfile_nc4), NF90_NOCLOBBER + NF90_HDF5, NCFOutID)

         !2) define dims
        ! status = NF_DEF_DIM(NCFOutID, 'nprogn' , N_progn_pert_max, Nprogn)
         status = NF90_DEF_DIM(NCFOutID, 'nseed' , NRANDSEED, seeddim)
         status = NF90_DEF_DIM(NCFOutID, 'lat' , N_lat, latdim)
         status = NF90_DEF_DIM(NCFOutID, 'lon' , N_lon, londim)
         status = NF90_DEF_DIM(NCFOutID, 'nforce' , N_force_pert_max, Nforce)
         status = NF90_DEF_DIM(NCFOutID, 'nprogn' , N_progn_pert_max, Nprogn)

         ! 3) define vars
         status = NF90_DEF_VAR(NCFOutID,'pert_rseed',NF90_DOUBLE,seeddim,seedid)
         status = NF90_PUT_ATT(NCFOutID, seedid, 'long_name','perturbations_rseed')
         status = NF90_PUT_ATT(NCFOutID, seedid, 'units', '1')

         dims(1)= londim
         dims(2)= latdim
         dims(3)= Nforce

         status = NF90_DEF_VAR(NCFOutID,'fpert_ntrmdt',NF90_REAL, dims, forceid)
         !status = nf90_def_var_deflate(NCFOutID, forceid, shuffle, deflate, deflate_level)
         status = NF90_PUT_ATT(NCFOutID, forceid, 'long_name', 'force_pert_intermediate')
         status = NF90_PUT_ATT(NCFOutID, forceid, 'units', '1')
         status = nf90_put_att(NCFOutID, forceid, '_FillValue', fill_value)
         dims(1)= londim
         dims(2)= latdim
         dims(3)= Nprogn

         status = NF90_DEF_VAR(NCFOutID, 'ppert_ntrmdt', NF90_REAL, dims, prognid)
         !status = nf90_def_var_deflate(NCFOutID, prognid, shuffle, deflate, deflate_level)
         status = NF90_PUT_ATT(NCFOutID, prognid, 'long_name', 'progn_pert_intermediate')
         status = NF90_PUT_ATT(NCFOutID, prognid, 'units', '1')
         status = nf90_put_att(NCFOutID, prognid, '_FillValue', fill_value)


         status = nf90_enddef(NCFOutID)
         ! 4) writing

         status= NF90_PUT_VAR(NCFOutID,seedid ,real(Pert_rseed,kind=8))

         xstart = 1 + pert_grid_f%i_offg
         ystart = 1 + pert_grid_f%j_offg

         ! will change to MAPL default 1.0e+15
         !do i = 1, N_lonf
         !   do j = 1, N_latf
         !      do k = 1, N_force_pert
         !         if (Force_pert_ntrmdt_f(i,j,k) < -9998) Force_pert_ntrmdt_f(i,j,k)=fill_value
         !      enddo
         !   enddo
         !enddo

         status= NF90_PUT_VAR(NCFOutID, forceid, Force_pert_ntrmdt_f, start=[xstart,ystart,1], &
                 count=[N_lonf,N_latf,N_force_pert])

         ! will change to MAPL default 1.0e+15
         !do i = 1, N_lonf
         !   do j = 1, N_latf
         !      do k = 1, N_progn_pert
         !         if (Progn_pert_ntrmdt_f(i,j,k) < -9998) Progn_pert_ntrmdt_f(i,j,k)=fill_value
         !      enddo
         !   enddo
         !enddo

         status= NF90_PUT_VAR(NCFOutID, prognid, Progn_pert_ntrmdt_f, start=[xstart,ystart,1], &
                 count=[N_lonf,N_latf,N_progn_pert])

         STATUS   = NF90_CLOSE (NCFOutID)
   
         deallocate(Force_pert_ntrmdt_f, Progn_pert_ntrmdt_f)

         rc = status
      end subroutine o_pert_GEOSldas
end subroutine convert_pert_rst
