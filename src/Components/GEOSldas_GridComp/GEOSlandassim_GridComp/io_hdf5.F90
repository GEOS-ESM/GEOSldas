! io_hdf5.F90

! pchakrab, 16 Jan 2014

!---------------------------
! STEPS to read 1D real data
! (For more details, please see the example at the end of this file)
!
! use io_hdf5
! type(hdf5read) :: h5r
! real, allocatable, dimension(:) :: data1D
!
! call h5r%openFile(file)
! for each dset in file
!     call h5r%queryDataset(dsetName, dsetRank, dsetSize)
!     allocate(data1D(dsetSize(1)))
!     call h5r%readDataset(data1D) [READ THE DATASET JUST QUERIED]
! end for
! call h5r%closeFile
! deallocate(data1D)
!---------------------------

module io_hdf5

  use hdf5
  use iso_fortran_env

  implicit none

  private

  integer, parameter :: UNINIT_INT = -99999
  character(len=*), parameter :: UNINIT_STR = ""

  type, public :: hdf5read
     private
     character(len=256) :: file_name = UNINIT_STR
     integer(hid_t) :: file_id = UNINIT_INT
     character(len=256) :: dset_name = UNINIT_STR
     integer(hid_t) :: dset_id = UNINIT_INT, dspace_id = UNINIT_INT, dtype_id = UNINIT_INT
     integer :: dset_rank = UNINIT_INT
     ! 7 is the max dimension of a fortran array
     integer(hsize_t) :: dset_size(7) = UNINIT_INT, dset_max_size(7) = UNINIT_INT
   contains
     ! public
     procedure, public  :: openFile
     procedure, public  :: closeFile
     procedure, public  :: queryDataset
     generic,   public  :: readDataset => readDataset1DReal, readDataset1DReal8, readDataset1DInt, readDataset1DChar24, readDataset2DReal
     ! private
     procedure, private :: readDataset1DReal
     procedure, private :: readDataset1DReal8
     procedure, private :: readDataset1DInt
     procedure, private :: readDataset1DChar24
     procedure, private :: readDataset2DReal
     procedure, private :: uninitDataset
  end type hdf5read

contains

  ! open file
  subroutine openFile(this, filename)

    ! input/output variables
    ! NEED class(hdf5read) instead of type(hdf5read)
    class (hdf5read), intent(inout) :: this
    character(len=*), intent(in) :: filename

    ! local variable
    integer :: hdf5err

    ! set obj param val
    this%file_name = filename

    ! initialize fortran interface
    call h5open_f(hdf5err)
    call checkErrCode_('h5open_f', hdf5err)

    ! open existing file
    call h5fopen_f(this%file_name, H5F_ACC_RDONLY_F, this%file_id, hdf5err)
    call checkErrCode_('h5fopen_f', hdf5err)

  end subroutine openFile

  ! close already opened file
  subroutine closeFile(this)

    ! input/output variables
    class (hdf5read), intent(inout) :: this

    ! local variable
    integer :: hdf5err

    ! ensure that dataset has been closed
    if (this%dset_name/=UNINIT_STR) stop "ERROR: Open dataset needs to be closed first. Stopping!"

    ! close file
    call h5fclose_f(this%file_id, hdf5err)
    call checkErrCode_('h5fclose_f', hdf5err)
    this%file_name = UNINIT_STR
    this%file_id = UNINIT_INT

    ! close fortran interface
    call h5close_f(hdf5err)
    call checkErrCode_('h5close_f', hdf5err)
    
  end subroutine closeFile

  ! query dataset for number of dims and its shape
  subroutine queryDataset(this, dsetName, dsetRank, dsetSize)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    character(len=*), intent(in) :: dsetName
    integer, intent(out) :: dsetRank
    integer, intent(out) :: dsetSize(7)

    ! local variable
    integer :: hdf5err

    ! ensure that file_name is set i.e. openFile
    ! must have been called prior to this routine
    if (this%file_name==UNINIT_STR) stop "ERROR: No open file available. Stopping!"

    ! set obj param val
    this%dset_name = dsetname

    ! open datset from already opened file
    call h5dopen_f(this%file_id, this%dset_name, this%dset_id, hdf5err)
    call checkErrCode_('h5dopen_f', hdf5err)

    ! get dataspace id
    call h5dget_space_f(this%dset_id, this%dspace_id, hdf5err)
    call checkErrCode_('h5dget_space_f', hdf5err)

    ! get num of dimensions
    call h5sget_simple_extent_ndims_f(this%dspace_id, this%dset_rank, hdf5err)
    call checkErrCode_('h5sget_simple_extent_ndims_f', hdf5err)
    dsetRank = this%dset_rank

    ! get size of array
    call h5sget_simple_extent_dims_f(this%dspace_id, this%dset_size, this%dset_max_size, hdf5err)
    call checkErrCode_('h5sget_simple_extent_dims_f', hdf5err)
    dsetSize = this%dset_size

  end subroutine queryDataset


  ! uninitalize dataset
  subroutine uninitDataset(this)

    ! input/output variables
    class (hdf5read), intent(inout) :: this

    ! un-initialize everything related to
    ! the dataset queried/read
    this%dset_name = UNINIT_STR
    this%dset_id = UNINIT_INT
    this%dspace_id = UNINIT_INT
    this%dset_rank = UNINIT_INT
    this%dset_size = UNINIT_INT
    this%dset_max_size = UNINIT_INT
    this%dtype_id = UNINIT_INT

  end subroutine uninitDataset


  ! read the dataset that was queried earlier
  subroutine readDataset1DChar24(this, dataChar)
    
    ! input/output variables
    class (hdf5read), intent(inout) :: this
    character(len=24), intent(out) :: dataChar(:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       call h5dread_f(this%dset_id, this%dtype_id, dataChar, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DChar24


  ! read the dataset that was queried earlier
  subroutine readDataset1DReal(this, data1D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    real, intent(out) :: data1D(:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       call h5dread_f(this%dset_id, this%dtype_id, data1D, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DReal


  ! read the dataset that was queried earlier
  subroutine readDataset1DReal8(this, data1D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    real(REAL64), intent(out) :: data1D(:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       call h5dread_f(this%dset_id, this%dtype_id, data1D, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DReal8


  subroutine readDataset1DInt(this, data1D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    integer, intent(out) :: data1D(:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       !call h5dread_f(this%dset_id, this%dtype_id, data1D, this%dset_size, hdf5err)
       call h5dread_f(this%dset_id, H5T_NATIVE_INTEGER, data1D, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DInt


  subroutine readDataset2DReal(this, data2D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    real, intent(out) :: data2D(:,:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       call h5dread_f(this%dset_id, this%dtype_id, data2D, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset2DReal

  ! check return code
  ! (not part of class hdf5read)
  subroutine checkErrCode_(routineName, hdf5errCode)

    ! input/output variables
    character(len=*), intent(in) :: routineName
    integer, intent(in) :: hdf5errCode

    if (hdf5errCode<0) then
       write(*,*) 'ERROR: ', routineName, ' returned NEGATIVE err code. Stopping!'
       stop
    end if

  end subroutine checkErrCode_

end module io_hdf5


! *****************************************************************

#ifdef TEST_IOHDF5

program test_read

  use io_hdf5

  implicit none

  character(len=*), parameter :: file_name = '/discover/nobackup/projects/gmao/smap/SMAP_L4/SMAP/L1C_TB/Y2001/M07/D20/SMAP_L1C_TB_02915_D_20010720T002132_D04003_000.h5'
  character(len=300) :: dsetName

  type(hdf5read) :: h5r
  integer :: dsetRank, dsetSize(7), i

  type myDataType
     real, pointer, dimension(:) :: tb_h_aft => null()     ! aft
     real, pointer, dimension(:) :: lon => null()
     integer, pointer, dimension(:) :: row => null()
     integer, pointer, dimension(:) :: flag => null()
     real(REAL64), pointer, dimension(:) :: tb_time => null()
     character(len=24), pointer, dimension(:) :: tb_time_utc_aft => null()
  end type MyDataType
  type(MyDataType), dimension(1) :: data
  
  print *, 'HDF5 file: ', trim(file_name)
  print *, ''

  ! open file
  call h5r%openFile(file_name)


  ! query dataset + allocate space + read data
  dsetName = '/Global_Projection/cell_tb_h_aft'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%tb_h_aft(dsetSize(1)))
  call h5r%readDataset(data(1)%tb_h_aft)
  print *, trim(dsetname),'(1:10)'
  print *, data(1)%tb_h_aft(1:10)
  print *, ''

  ! query dataset + allocate space + read data
  dsetName = '/Global_Projection/cell_lon'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%lon(dsetSize(1)))
  call h5r%readDataset(data(1)%lon)
  print *, trim(dsetname),'(1:10)'
  print *, data(1)%lon(1:10)
  print *, ''

  ! query dataset + allocate space + read data
  dsetName = '/Global_Projection/cell_row'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%row(dsetSize(1)))
  call h5r%readDataset(data(1)%row)
  print *, trim(dsetname),'(201:210)'
  print *, data(1)%row(201:210)
  print *, ''

  ! query dataset + allocate space + read data
  dsetName = '/Global_Projection/cell_tb_time_seconds_aft'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%tb_time(dsetSize(1)))
  call h5r%readDataset(data(1)%tb_time)
  print *, trim(dsetname),'(1:10)'
  print *, data(1)%tb_time(1:10)
  print *, ''

  dsetName = '/Global_Projection/cell_tb_time_utc_aft'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%tb_time_utc_aft(dsetSize(1)))
  call h5r%readDataset(data(1)%tb_time_utc_aft)
  print *, trim(dsetname),'(241:250)'
  do i=241,250
     print *, data(1)%tb_time_utc_aft(i)
  end do

  dsetName = '/Global_Projection/cell_tb_qual_flag_h_aft'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%flag(dsetSize(1)))
  call h5r%readDataset(data(1)%flag)
  print *, trim(dsetname),'(241:250)'
  do i=241,250
     print *, data(1)%flag(i)
  end do


  ! close file
  call h5r%closeFile


  ! deallocate memory
  if (associated(data(1)%tb_h_aft)) deallocate(data(1)%tb_h_aft)
  if (associated(data(1)%lon)) deallocate(data(1)%lon)
  if (associated(data(1)%row)) deallocate(data(1)%row)
  if (associated(data(1)%flag)) deallocate(data(1)%flag)
  if (associated(data(1)%tb_time)) deallocate(data(1)%tb_time)
  if (associated(data(1)%tb_time_utc_aft)) deallocate(data(1)%tb_time_utc_aft)

end program test_read

#endif

! =================== EOF ==========================================
