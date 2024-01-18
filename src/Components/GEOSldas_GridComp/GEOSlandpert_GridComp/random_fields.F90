! =========================================================================
!
! random_fields.f90
!
! random field generator in 2d:
!  generate a pair of random fields in 2d with zero mean
!
! subroutines rfg2d_fft() and sqrt_gauss_spectrum are translated from 
!  C++ code rfg2d.C written for MIT EnKF work by reichle
!  (see janus:~reichle/nasa/EnKF)
!
! covariance is specified through its spectrum, so far only Gaussian
!
! IMPORTANT: read comments for function rfg2d_fft()
!
! written for NSIPP - EnKF
! Type:   f90
! Author: Rolf Reichle
! Date:   2 Nov 2001
!
! reichle, 18 Feb 2005 - updated for use with module nr_ran2_gasdev
!                        deleted use of module select_kinds
! reichle,  8 Jun 2005 - added nr_fft should CXML fft not be available
! pchakrab, 17 Jun 2013 - redesigned to make the module object-oriented (F03)

! use intel mkl fft when available
#ifdef MKL_AVAILABLE
#include "mkl_dfti.f90"
#else
#define NR_FALLBACK
#endif

#include "MAPL_ErrLog.h"
#include "unused_dummy.H"

module random_fieldsMod
  
#ifdef MKL_AVAILABLE
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer, c_ptr, c_sizeof, C_NULL_PTR
  use mpi
  use MKL_DFTI
#else
  use nr_fft,                           ONLY:     &
       fourn
#endif
  
  use nr_ran2_gasdev,                   ONLY:     &
       NRANDSEED,                                 &
       nr_ran2_2d,                                &
       nr_gasdev

  use MAPL_ExceptionHandling

  implicit none
  
  private

  real, parameter :: TWO_PI = 2.*3.14159265
  real, parameter :: SQRT2 = sqrt(2.0)
  
  type, public :: random_fields
     private
     integer :: N_x, N_y
     integer :: N_x_fft, N_y_fft ! computed by calc_fft_grid
     real, allocatable :: field1_fft(:,:), field2_fft(:,:)
     integer :: fft_lens(2) ! length of each dim for 2D transform
#ifdef MKL_AVAILABLE
     integer :: comm
     integer :: node_comm
     integer :: win
     type (c_ptr) :: base_address 
     
     type(DFTI_DESCRIPTOR), pointer :: Desc_Handle
     type(DFTI_DESCRIPTOR), pointer :: Desc_Handle_dim1
     type(DFTI_DESCRIPTOR), pointer :: Desc_Handle_dim2
     integer, allocatable :: dim1_counts(:)
     integer, allocatable :: dim2_counts(:)
#endif

   contains
     
     !     procedure, public  :: initialize
     procedure, public  :: finalize
     procedure, public  :: rfg2d_fft
     procedure, public  :: generate_white_field
     procedure, private :: sqrt_gauss_spectrum_2d
#ifdef MKL_AVAILABLE
     procedure, private :: win_allocate
     procedure, private :: win_deallocate
#endif
  end type random_fields

  interface random_fields
     module procedure new_random_fields
  end interface random_fields
  
contains
  
  ! constructor (set parameter values), allocate memory
  function new_random_fields(Nx, Ny, Nx_fft, Ny_fft, comm, rc) result (rf)
    
    ! input/output variables [NEED class(random_fields)
    !   instead of type(random_fields)] - F2003 quirk?!?
    type(random_fields)            :: rf
    integer,           intent(in)  :: Nx, Ny, Nx_fft, Ny_fft
    integer, optional, intent(in)  :: comm
    integer, optional, intent(out) :: rc 
    
    ! local variables
    integer :: status, ierror
    integer :: rank, npes, local_dim1, local_dim2, remainder
    integer :: Stride(2)
    
    ! set obj param vals
    rf%N_x = Nx
    rf%N_y = Ny

    ! ensure N_x_fft, N_y_fft are powers of two
    rf%N_x_fft = Nx_fft
    rf%N_y_fft = Ny_fft


    ! allocate memory
    allocate(rf%field1_fft(rf%N_x_fft, rf%N_y_fft))
    allocate(rf%field2_fft(rf%N_x_fft, rf%N_y_fft))

#ifdef MKL_AVAILABLE
    if (present(comm)) then
       rf%comm = comm
       
       call MPI_Comm_split_type(rf%comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, rf%Node_Comm,ierror)
       call MPI_Comm_size(rf%Node_Comm, npes, ierror)
       call MPI_Comm_rank(rf%Node_Comm, rank, ierror)
       
       if (npes > minval([Nx_fft, Ny_fft]) ) then
          print*, " Two many processors are acquired  in a node for parallel FFT"
          print*, " The number of processors acquired in a node should be smaller than or equal to FFT grid size: ", minval([Nx_fft, Ny_fft]) 
          _FAIL('Parallel FFT failed')
       endif

       call rf%win_allocate(Nx_fft, Ny_fft, _RC)

       ! distribution of the grid for fft
       allocate(rf%dim1_counts(npes),rf%dim2_counts(npes))
       local_dim1 = Nx_fft/npes
       rf%dim1_counts = local_dim1
       remainder = mod(Nx_fft, npes)
       rf%dim1_counts(1:remainder) = local_dim1 + 1
       local_dim1 = rf%dim1_counts(rank+1) 

       local_dim2 = Ny_fft/npes
       rf%dim2_counts = local_dim2
       remainder = mod(Ny_fft, npes)
       rf%dim2_counts(1:remainder) = local_dim2 + 1
       local_dim2 = rf%dim2_counts(rank+1)


       status = DftiCreateDescriptor(rf%Desc_Handle_Dim1, DFTI_SINGLE,&
                                DFTI_COMPLEX, 1, Nx_fft )
       _VERIFY(status)
       status = DftiCreateDescriptor(rf%Desc_Handle_Dim2, DFTI_SINGLE,&
                                DFTI_COMPLEX, 1, Ny_fft )
       _VERIFY(status)

       ! perform local_dim2 one-dimensional transforms along 1st dimension
       status = DftiSetValue( rf%Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, local_dim2 )
       _VERIFY(status)
       status = DftiSetValue( rf%Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, Nx_fft )
       _VERIFY(status)
       status = DftiSetValue( rf%Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, Nx_fft )
       _VERIFY(status)
       status = DftiCommitDescriptor( rf%Desc_Handle_Dim1 )
       _VERIFY(status)
       ! status = DftiComputeForward( rf%Desc_Handle_Dim1, X )
       ! local_dim1 one-dimensional transforms along 2nd dimension
       Stride(1) = 0; Stride(2) = local_dim1
       status = DftiSetValue( rf%Desc_Handle_Dim2, DFTI_NUMBER_OF_TRANSFORMS, local_dim1)
       _VERIFY(status)
       status = DftiSetValue( rf%Desc_Handle_Dim2, DFTI_INPUT_DISTANCE, 1 )
       _VERIFY(status)
       status = DftiSetValue( rf%Desc_Handle_Dim2, DFTI_OUTPUT_DISTANCE, 1 )
       _VERIFY(status)
       status = DftiSetValue( rf%Desc_Handle_Dim2, DFTI_INPUT_STRIDES, Stride )
       _VERIFY(status)
       status = DftiSetValue( rf%Desc_Handle_Dim2, DFTI_OUTPUT_STRIDES, Stride )
       _VERIFY(status)
       status = DftiCommitDescriptor( rf%Desc_Handle_Dim2 )
       _VERIFY(status)
       !status = DftiComputeForward( rf%Desc_Handle_Dim2, X )
    else
       rf%comm = MPI_COMM_NULL
       ! allocate mem and init mkl dft
       status = DftiCreateDescriptor(rf%Desc_Handle, DFTI_SINGLE, DFTI_COMPLEX, 2, [Nx_fft, Ny_fft])
       _VERIFY(status)

       ! initialize for actual dft computation
       status = DftiCommitDescriptor(rf%Desc_Handle)
       _VERIFY(status)
    endif
#endif
    _RETURN(_SUCCESS)
  end function new_random_fields
  
  ! **************************************************************************
  
  ! destructor - deallocate memory
  subroutine finalize(this, rc)

    ! input/output variables
    class(random_fields), intent(inout) :: this
    integer, optional, intent(out) :: rc
    ! local variable
    integer :: status

    ! deallocate memory
    if(allocated(this%field1_fft)) deallocate(this%field1_fft)
    if(allocated(this%field2_fft)) deallocate(this%field2_fft)
    
#ifdef MKL_AVAILABLE
    if (this%comm == MPI_COMM_NULL) then
       status = DftiFreeDescriptor(this%Desc_Handle)
       _VERIFY(status)
    else
       
       status = DftiFreeDescriptor(this%Desc_Handle_dim1)
       _VERIFY(status)
       status = DftiFreeDescriptor(this%Desc_Handle_dim2)
       _VERIFY(status)
       
       call this%win_deallocate( _RC)
       
       deallocate(this%dim1_counts, this%dim2_counts)
    endif
#endif

  end subroutine finalize
  

  ! subroutine sqrt_gauss_spectrum_2d()
  !
  ! get SQUARE ROOT of 2d Gaussian spectrum (incl volume element)
  !
  ! 2d Gaussian spectrum:
  !   
  ! S(kx,ky) = variance
  !            *
  !            lambda_x*lambda_y/(2*pi) 
  !            * 
  !            exp( -(kx^2*lambda_x^2 + ky^2*lambda_y^2)/2 )
  ! 
  ! return: sqrt( S*dkx*dky )
  !
  ! that is return the SQUARE ROOT of the spectrum multiplied with the
  !  square root of the volume element d2k=dkx*dky of the ifft integral
  !
  ! spectrum is returned in "wrap-around" order compatible with CXML and 
  !  matlab FFT algorithms
  ! 
  ! inputs:
  !  variance : variance desired for complex field, if pair of real fields 
  !             is used each field must eventually be multiplied with sqrt(2)
  !  N_x      : number of nodes in x direction
  !  N_y      : number of nodes in y direction
  !  dkx      : wave number spacing in x direction
  !  dky      : wave number spacing in y direction
  !  lambda_x : decorrelation length in x direction 
  !  lambda_y : decorrelation length in y direction
  !
  ! modifies this%field1_fft

  subroutine sqrt_gauss_spectrum_2d(this, lx, ly, dx, dy)

    ! input/output variables
    class(random_fields), intent(inout) :: this
    real,                 intent(in)    :: lx, ly, dx, dy

    ! local variables
    real    :: dkx, dky, fac, lamx2dkx2, lamy2dky2
    real    :: lx2kx2(this%N_x_fft), ly2ky2(this%N_y_fft)
    integer :: i, j, i1, i2, rank, ierror
    real    :: var

    var = 1.0
    ! start
    dkx = (TWO_PI)/(float(this%N_x_fft)*dx)
    dky = (TWO_PI)/(float(this%N_y_fft)*dy)

    ! factor includes sqrt of volume element of ifft integral
    fac = sqrt(var*lx*ly/(TWO_PI)*dkx*dky )
    lamx2dkx2 = lx*lx*dkx*dkx
    lamy2dky2 = ly*ly*dky*dky

    ! precompute (lambda_x*k_x)^2 in "wrap-around"
    ! order suitable for CXML fft
    do i=1,(this%N_x_fft/2)
       lx2kx2(i) = lamx2dkx2*(i-1)*(i-1)
    end do
    do i=(this%N_x_fft/2+1),this%N_x_fft
       lx2kx2(i) = lamx2dkx2*(this%N_x_fft-i+1)*(this%N_x_fft-i+1)
    end do

    ! precompute (lambda_y*k_y)^2 in "wrap-around"
    ! order suitable for CXML fft
    do j=1,(this%N_y_fft/2)
       ly2ky2(j) = lamy2dky2*(j-1)*(j-1)
    end do
    do j=(this%N_y_fft/2+1),this%N_y_fft
       ly2ky2(j) = lamy2dky2*(this%N_y_fft-j+1)*(this%N_y_fft-j+1)
    end do

    ! assemble spectrum in "wrap-around" order
    i1 = 1
    i2 = this%N_x_fft
    if (this%comm /= MPI_COMM_NULL) then
       call MPI_COMM_Rank(this%node_comm, rank, ierror)
       i1 = sum(this%dim1_counts(1:rank)) + 1
       i2 = sum(this%dim1_counts(1:rank+1))
    endif
 
    do j=1,this%N_y_fft
       this%field1_fft(i1:i2,j) = fac*exp(-.25*(lx2kx2(i1:i2)+ly2ky2(j)))  
    end do

    return

  end subroutine sqrt_gauss_spectrum_2d
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! ----------------------------------------------------------------------
  !
  ! subroutine rfg2d_fft()
  !  
  ! generate a pair of 2d zero-mean random fields using FFT method
  ! (so far only Gaussian covariance implemented)
  !
  ! NOTE: implemented with index counters of type int, must have
  !
  !          N_x*N_y < maximum integer on given machine
  !
  !       on yama/alborz (at MIT) int varies from -2147483648...2147483647
  !       -> can handle up to N_x*N_y = (32768)^2
  !       if larger fields are needed, rewrite with type long int etc
  !
  ! NOTE: The fft method for generating random fields produces 
  !       fields that are periodic with the size of the domain,
  !       that is the field is (almost) the same on each boundary.
  !       This introduces unwanted correlations at lags shorter than
  !       the domain size. Therefore, only a part of the generated
  !       field is usable. As a rule of thumb, the fields should be
  !       generated on a grid that is two correlation lenghts bigger
  !       than the field on which the grid is desired. Then cut out
  !       fields of the necessary size.
  !       This procedure is included in rfg2d_fft().
  !
  ! NOTE: The variance specified as input is the theoretical variance of 
  !       the complex field that is obtained from the inverse fft of the 
  !       realization dZ.
  !       The sample variance of this *complex* field is a FIXED (non-random) 
  !       number which depends on the size of the domain, the grid spacing,
  !       and the correlation length (but not on the random seed!!).
  !       (This number is non-random because the variance is the integral
  !        of the absolute value of the spectrum. In this integral the 
  !        randomness disappears because we only choose a random phase angle.)
  !       For vanishing discretization and spectral truncation error,
  !       this number converges to the theoretical (input) value.
  !
  !       This function is set up to use the pair of two real fields, where
  !
  !                field1 = sqrt(2)*real(ifft(dZ))
  !                field2 = sqrt(2)*imag(ifft(dZ)).
  !     
  !       The factor sqrt(2) re-scales the variance of field1 and field2
  !       such that for vanishing discretization error and spectral
  !       truncation error each field converges to the specified theoretical
  !       (input) variance. 
  !       NOTE: The sum of the sample variances of the two real fields 
  !       is equal to the (FIXED) sample variance of the complex field
  !       (before re-scaling with sqrt(2)). 
  !       The individual sample variances within each pair vary from 
  !       realization to realization.
  
  subroutine rfg2d_fft(this, rseed, rfield, rfield2, lx, ly, dx, dy)
    
    ! input/output variables
    class(random_fields),                  intent(inout) :: this ! ffield*_fft is modified
    integer,                               intent(inout) :: rseed(NRANDSEED) ! nr_ran2 modifies rseed
    real,    dimension(this%N_x,this%N_y), intent(out)   :: rfield, rfield2
    real,                                  intent(in)    :: lx, ly, dx, dy

    ! local variables
    !real :: theta, ran_num ! rng
    real,allocatable,dimension(:,:) :: theta, ran_num ! rng
    integer :: N_xy_fft, k ! fft
    integer :: i, j
    integer :: N_x_fft, N_y_fft
    real :: N_xy_fft_real
#ifdef MKL_AVAILABLE
    integer :: status
    complex, allocatable :: z_inout(:)
    complex, pointer :: tmp_field(:,:)
    complex, pointer :: tmp_field_dim1(:,:)
    complex, pointer :: tmp_field_dim2(:,:)
    integer :: n1, n2, npes, rank, ldim1, ldim2, ierror
    complex, pointer  :: X(:)
    type (c_ptr) :: cptr
#else
    real, allocatable :: tmpdata(:)
#endif

    ! start
    N_x_fft = this%N_x_fft
    N_y_fft = this%N_y_fft

    ! follow Ruan & McLaughlin, 1998:
    ! compute dZ = H * exp(i*theta) * sqrt(d2k)     
    ! start with square root of spectrum (factor H*sqrt(d2k)), put into field1
    ! modify this%field1_fft
    call this%sqrt_gauss_spectrum_2d(lx, ly, dx, dy)

    ! multiply by random phase angle
   !! do j=1,N_y_fft
   !!   do i=1,N_x_fft
   !!      call nr_ran2(rseed, ran_num)
   !!       theta = (TWO_PI)*ran_num ! random phase angle
   !!       this%field2_fft(i,j) = sin(theta)*this%field1_fft(i,j)
   !!       this%field1_fft(i,j) = cos(theta)*this%field1_fft(i,j)          
   !!    end do
   !! end do

    allocate(  theta(N_x_fft, N_y_fft))
    allocate(ran_num(N_x_fft, N_y_fft))

    call nr_ran2_2d(N_x_fft, N_y_fft, rseed, ran_num)
    theta = (TWO_PI)*ran_num ! random phase angle
    n1 = 1
    n2 = N_x_fft
#ifdef MKL_AVAILABLE
    if (this%comm /= MPI_COMM_NULL) then
       call MPI_comm_rank(this%node_comm, rank, ierror)
       n1 = sum(this%dim1_counts(1:rank)) + 1
       n2 = sum(this%dim1_counts(1:rank+1)) 
    endif
#endif
    this%field2_fft(n1:n2,:) = sin(theta(n1:n2,:))*this%field1_fft(n1:n2,:)
    this%field1_fft(n1:n2,:) = cos(theta(n1:n2,:))*this%field1_fft(n1:n2,:)

    deallocate(  theta)
    deallocate(ran_num)


    ! force dZ(1,1) to zero 
    ! (zero mean random field)
    this%field1_fft(1,1) = 0.
    this%field2_fft(1,1) = 0.
    
    ! apply 2D FFT
    N_xy_fft      = N_x_fft*N_y_fft
    N_xy_fft_real = real(N_xy_fft)

#ifdef MKL_AVAILABLE
    ! use MKL FFT
    ! fill temporary 1D array
    if (this%comm == MPI_COMM_NULL) then
       allocate(z_inout(N_xy_fft))
       k = 0
       do j=1,N_y_fft
          do i=1,N_x_fft
             k=k+1
             z_inout(k) = cmplx(this%field1_fft(i,j),this%field2_fft(i,j))
          end do
       end do    

       ! compute in-place backward transform (scale=1)
       ! NOTE: MKL backward transform is the same as NR forward transform
       status = DftiComputeBackward(this%Desc_Handle, z_inout)
       if (status/= DFTI_NO_ERROR) call quit('DftiComputeBackward failed!')

       ! extract random fields from z_inout
       z_inout = z_inout/N_xy_fft_real
       k = 0
       do j=1,N_y_fft
          do i=1,N_x_fft      
             k=k+1
             this%field1_fft(i,j) = real(z_inout(k))
             this%field2_fft(i,j) = aimag(z_inout(k))
          end do
       end do

       deallocate(z_inout)
    else
       call MPI_comm_size(this%node_comm, npes, ierror)
       call c_f_pointer(this%base_address, tmp_field, shape=[N_x_fft, N_y_fft])
       ldim1  = this%dim1_counts(rank+1)

       allocate(tmp_field_dim1(ldim1, N_y_fft))
       tmp_field_dim1 = cmplx(this%field1_fft(n1:n2,:),this%field2_fft(n1:n2,:))
       cptr = c_loc(tmp_field_dim1(1,1))
       call c_f_pointer (cptr, X, [ldim1*N_y_fft])
       status = DftiComputeBackward( this%Desc_Handle_Dim2, X )
       if (status/= DFTI_NO_ERROR) call quit('DftiComputeBackward dim2 failed!')
       call MPI_Barrier(this%node_comm, ierror)
       tmp_field(n1:n2,:) = tmp_field_dim1

       call MPI_Win_fence(0, this%win, ierror)

       n1 = sum(this%dim2_counts(1:rank)) + 1
       n2 = sum(this%dim2_counts(1:rank+1)) 
       ldim2  = this%dim2_counts(rank+1)
       allocate(tmp_field_dim2(N_x_fft, ldim2))
       tmp_field_dim2 = tmp_field(:,n1:n2)
       cptr = c_loc(tmp_field_dim2(1,1))
       call c_f_pointer (cptr, X, [N_x_fft*ldim2])
       status = DftiComputeBackward( this%Desc_Handle_Dim1, X )
       if (status/= DFTI_NO_ERROR) call quit('DftiComputeBackward dim1 failed!')
       tmp_field(:,n1:n2) = tmp_field_dim2/N_xy_fft_real
       
       call MPI_Win_fence(0, this%win, ierror)

       this%field1_fft = real(tmp_field)
       this%field2_fft = aimag(tmp_field)

       deallocate(tmp_field_dim1, tmp_field_dim2)
    endif
#else  
    ! use nr_fft
    ! fill tmpdata according to Figs 12.2.2
    ! and 12.4.1 of f77 NR book
    allocate(tmpdata(2*N_xy_fft))
    k=0
    do j=1,N_y_fft
       do i=1,N_x_fft
          k=k+1
          tmpdata(k) = this%field1_fft(i,j)
          k=k+1
          tmpdata(k) = this%field2_fft(i,j)
       end do
    end do
    
    ! apply nr_fft
    call fourn(tmpdata,this%fft_lens,2,1)
    
    ! extract random fields from tmpdata
    k=0
    do j=1,N_y_fft
       do i=1,N_x_fft      
          k=k+1
          this%field1_fft(i,j) = tmpdata(k)/N_xy_fft_real
          k=k+1
          this%field2_fft(i,j) = tmpdata(k)/N_xy_fft_real
       end do
    end do
    
    deallocate(tmpdata)

#endif    

    ! multiply with factor sqrt(2) to get correct variance
    ! (see above and p. 388 Ruan and McLaughlin, 1998),
    ! also multiply with N_x_fft*N_y_fft to get correct scaling,
    ! also retain only usable part of field?_fft
    ! output variables
    rfield  = SQRT2*N_xy_fft_real*this%field1_fft(1:this%N_x,1:this%N_y)
    rfield2 = SQRT2*N_xy_fft_real*this%field2_fft(1:this%N_x,1:this%N_y)
    
  end subroutine rfg2d_fft
  
  
  
  ! generate standard-normal random field that is white in space
  !
  ! note that nr_gasdev always produces a pair of random numbers
  !
  ! do not store random numbers between subsequent calls to 
  ! the random field generator - works best if fields are large
  ! (ie. avoid using this subroutine with N_x=N_y=1)
  !
  ! revised to avoid branching (if statement w/in do loop)
  ! - pchakrab+reichle, 29 Nov 2013
  !
  subroutine generate_white_field(this, rseed, rfield)
    
    implicit none

    ! input/output variables
    class(random_fields), intent(in) :: this
    integer, intent(inout) :: rseed(NRANDSEED) ! nr_gasdev modifies rseed
    real, dimension(this%N_x,this%N_y), intent(out), target :: rfield
    
    ! local variables
    integer :: Nxy, index
    real, pointer :: ptr2rfield(:)
    logical :: NxyIsOdd
    real :: tmp_real(2)

    Nxy = this%N_x*this%N_y
    
    ptr2rfield(1:Nxy) => rfield(:,:) ! ptr rank remapping

    NxyIsOdd = .false.
    if (mod(Nxy,2)==1) NxyIsOdd = .true.

    do index=1,Nxy-1,2
       call nr_gasdev(rseed, tmp_real)
       ptr2rfield(index) = tmp_real(1)
       ptr2rfield(index+1) = tmp_real(2)
    end do
    if (NxyIsOdd) then
       call nr_gasdev(rseed, tmp_real)
       ptr2rfield(Nxy) = tmp_real(1)
    end if

  end subroutine generate_white_field

  ! a local, small stop routine
  subroutine quit(message)
    
    implicit none
    
    character(*), intent(in) :: message

    write (*,*) trim(message)
    stop
    
  end subroutine quit

  subroutine win_allocate(this, nx, ny, rc)
     class(random_fields), intent(inout) :: this
     integer, intent(in) :: nx, ny
     integer, optional, intent(out) :: rc
     complex :: dummy
     integer(kind=MPI_ADDRESS_KIND) :: windowsize
     integer :: disp_unit,status, Rank
     integer(kind=MPI_ADDRESS_KIND) :: n_bytes 


     call MPI_Comm_rank( this%node_comm, rank, status)
     n_bytes = nx*ny*c_sizeof(dummy)
     windowsize = 0_MPI_ADDRESS_KIND
     if (Rank == 0) windowsize = n_bytes
     disp_unit  = 4
     call MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, this%node_comm, &
              this%base_address, this%win, status)
     _VERIFY(status)
     if (rank /=0)  CALL MPI_Win_shared_query(this%win, 0, windowsize, disp_unit, this%base_address, status)
     call MPI_Win_fence(0, this%win, status)
     _VERIFY(status)
     call MPI_Barrier(this%node_comm, status)
     _VERIFY(status)
     _RETURN(_SUCCESS)
  end subroutine win_allocate

  subroutine win_deallocate(this, rc)
     class(random_fields), intent(inout) :: this
     integer, optional, intent(out) :: rc
     integer :: status
     call MPI_Win_fence(0, this%win, status)
     _VERIFY(status)
     call MPI_Win_free(this%win, status)
     _VERIFY(status)
     call MPI_comm_free(this%node_comm, status)
     _VERIFY(status)

  end subroutine win_deallocate
  
end module Random_fieldsMod

module StringRandom_fieldsMapMod
   use Random_fieldsMod

#include "types/key_deferredLengthString.inc"
#define _value type (random_fields)
#define _value_equal_defined

#define _map StringRandom_fieldsMap
#define _iterator StringRandom_fieldsMapIterator

#define _alt

#include "templates/map.inc"

#undef _alt
#undef _iterator
#undef _map
#undef _value
#undef _key
#undef _value_equal_defined
end module StringRandom_fieldsMapMod


#ifdef TEST_RFG2D

!program test_rfg2d
!  
!  use Random_fieldsMod
!  use nr_ran2_gasdev
!  
!  implicit none
!  
!  integer :: N_x, N_y, i, j, n_e, N_e_tot
!  real :: dx, dy, lx, ly, var
!  real, allocatable, dimension(:,:) :: field1, field2
!  
!  character(300) :: file_name
!  character(10)  :: n_e_string
!  character(100) :: output_format
!  character(10)  :: tmp_string
!
!  integer :: RSEEDCONST
!  integer, dimension(NRANDSEED) :: rseed
!  
!  character(5) :: fft_tag
!    
!  ! instance of random_fields
!  type(random_fields) :: rf
!
!  ! start
!  RSEEDCONST = -777
!  rseed(1) = RSEEDCONST
!  write (*,*) RSEEDCONST
!  call init_randseed(rseed)
!
!  N_x = 144
!  N_y = 91
!  dx = 5000.
!  dy = 5000.
!  lx = 45000.
!  ly = 45000.
!  var = 1.
!  
!  
!  allocate(field1(N_x,N_y))
!  allocate(field2(N_x,N_y))
!  
!#ifdef MKL_AVAILABLE
!  fft_tag = 'mklx.'
!#else
!  fft_tag = 'nrxx.'
!#endif
!
!  ! get N_e fields
!  N_e_tot = 10
!  do n_e=1,N_e_tot,2
!
!     rf = random_fields(N_x, N_y, Nx_fft, Ny_fft)
!     call rf%rfg2d_fft(rseed, field1, field2, lx, ly, dx, dy)
!     !call rf%generate_white_field(rseed, field1)
!     call rf%finalize
!      
!     ! write to file
!     ! field1
!     write(n_e_string,  '(i3.3)') n_e
!     file_name = 'rf.'//fft_tag// n_e_string(1:len_trim(n_e_string)) // '.dat'
!     write(tmp_string, '(i3.3)') N_y
!     output_format = '(' // tmp_string(1:len_trim(tmp_string)) // '(1x,e13.5))'
!     open (10,file=file_name,status='unknown')
!     do i=1,N_x
!        write (10,output_format(1:len_trim(output_format))) (field1(i,j), j=1,N_y)
!     end do
!     close (10,status='keep')
!
!     ! field2
!     write(n_e_string,  '(i3.3)') n_e+1
!     file_name = 'rf.' //fft_tag// n_e_string(1:len_trim(n_e_string)) // '.dat'
!     write(tmp_string, '(i3.3)') N_y
!     output_format = '(' // tmp_string(1:len_trim(tmp_string)) // '(1x,e13.5))'
!     open (10,file=file_name,status='unknown')
!     do i=1,N_x
!        write (10,output_format(1:len_trim(output_format))) (field2(i,j), j=1,N_y)
!     end do
!     close (10,status='keep')
!     
!  end do
!
!end program test_rfg2d


#endif


! ======= EOF ==================================================

