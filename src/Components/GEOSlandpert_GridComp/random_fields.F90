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

module random_fields_class

#ifdef MKL_AVAILABLE
  use MKL_DFTI
#else
  use nr_fft,                           ONLY:     &
       fourn
#endif
  
  use nr_ran2_gasdev,                   ONLY:     &
       NRANDSEED,                                 &
       nr_ran2_2d,                                   &
       nr_gasdev
  
  implicit none
  
  private

  real, parameter :: TWO_PI = 2.*3.14159265
  real, parameter :: SQRT2 = sqrt(2.0)
  
  type, public :: random_fields
     private
     integer :: N_x, N_y
     real :: var, lx, ly, dx, dy
     integer :: N_x_fft, N_y_fft ! computed by calc_fft_grid
     real, allocatable :: field1_fft(:,:), field2_fft(:,:)
     integer :: fft_lens(2) ! length of each dim for 2D transform
#ifdef MKL_AVAILABLE
     type(DFTI_DESCRIPTOR), pointer :: Desc_Handle
#endif
   contains
     procedure, public  :: initialize
     procedure, public  :: finalize
     procedure, public  :: rfg2d_fft
     procedure, public  :: generate_white_field
     procedure, private :: sqrt_gauss_spectrum_2d
     procedure, private :: calc_fft_grid
  end type random_fields
  
contains

  ! constructor (set parameter values), allocate memory
  subroutine initialize(this, Nx, Ny, var, lx, ly, dx, dy)

    ! input/output variables [NEED class(random_fields)
    ! instead of type(random_fields)] - F2003 quirk?!?
    class(random_fields), intent(inout) :: this
    integer, intent(in) :: Nx, Ny
    real, intent(in) :: var, lx, ly, dx, dy
    
    ! local variable
    integer :: mklstat

    ! set obj param vals
    this%N_x = Nx
    this%N_y = Ny
    this%var = var
    this%dx = dx
    this%dy = dy
    this%lx = lx
    this%ly = ly

    ! calculate fft grid (N_x_fft, N_y_fft)
    call this%calc_fft_grid

    ! lengths of transform in each dimension
    this%fft_lens(1) = this%N_x_fft
    this%fft_lens(2) = this%N_y_fft

    ! allocate memory
    allocate(this%field1_fft(this%N_x_fft, this%N_y_fft))
    allocate(this%field2_fft(this%N_x_fft, this%N_y_fft))

#ifdef MKL_AVAILABLE
    ! allocate mem and init mkl dft
    mklstat = DftiCreateDescriptor(this%Desc_Handle, DFTI_SINGLE, DFTI_COMPLEX, 2, this%fft_lens)
    if (mklstat/=DFTI_NO_ERROR) call quit('DftiCreateDescriptor failed!')

    ! initialize for actual dft computation
    mklstat = DftiCommitDescriptor(this%Desc_Handle)
    if (mklstat/=DFTI_NO_ERROR) call quit('DftiCommitDescriptor failed!')
#endif

  end subroutine initialize



  ! destructor - deallocate memory
  subroutine finalize(this)

    ! input/output variables
    class(random_fields), intent(inout) :: this

    ! local variable
    integer :: mklstat

    ! deallocate memory
    if(allocated(this%field1_fft)) deallocate(this%field1_fft)
    if(allocated(this%field2_fft)) deallocate(this%field2_fft)

#ifdef MKL_AVAILABLE
    mklstat = DftiFreeDescriptor(this%Desc_Handle)
    if (mklstat/=DFTI_NO_ERROR) call quit('DftiFreeDescriptor failed!')
#endif

  end subroutine finalize



  ! calculate fft grid (N_x_fft, N_y_fft) that extends
  ! beyond the desired random field by about two correlation
  ! lengths. its dimensions should be powers of 2
  subroutine calc_fft_grid(this)

    ! input/output variables
    class(random_fields), intent(inout) :: this

    ! local variables
    real, parameter :: mult_of_xcorr = 2.
    real, parameter :: mult_of_ycorr = 2.
    integer :: Nx_fft, Ny_fft

    ! add minimum required correlation lengths 
    Nx_fft = this%N_x + ceiling(mult_of_xcorr*this%lx/this%dx)
    Ny_fft = this%N_y + ceiling(mult_of_ycorr*this%ly/this%dy)

    ! ensure N_x_fft, N_y_fft are powers of two
    this%N_x_fft = 2**ceiling(log(real(Nx_fft))/log(2.))
    this%N_y_fft = 2**ceiling(log(real(Ny_fft))/log(2.))

#if TEST_RFG2D
    write (*,*)
    write (*,*) 'desired random field:'
    write (*,*) 'N_x   = ', this%N_x, ' N_y   = ', this%N_y
    write (*,*) 'dx    = ', this%dx,  ' dy    = ', this%dy
    write (*,*) 'xcorr = ', this%lx,  ' ycorr = ', this%ly
    write (*,*)
    write (*,*) 'grid used for fft: '
    write (*,*) 'N_x_fft = ', this%N_x_fft,  ' N_y_fft = ', this%N_y_fft
    write (*,*)
#endif

  end subroutine calc_fft_grid



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
  subroutine sqrt_gauss_spectrum_2d(this)

    ! input/output variables
    class(random_fields), intent(inout) :: this

    ! local variables
    real :: dkx, dky, fac, lamx2dkx2, lamy2dky2
    real :: lx2kx2(this%N_x_fft), ly2ky2(this%N_y_fft)
    integer :: i, j

    ! start
    dkx = (TWO_PI)/(float(this%N_x_fft)*this%dx)
    dky = (TWO_PI)/(float(this%N_y_fft)*this%dy)

    ! factor includes sqrt of volume element of ifft integral
    fac = sqrt(this%var*this%lx*this%ly/(TWO_PI)*dkx*dky )
    lamx2dkx2 = this%lx*this%lx*dkx*dkx
    lamy2dky2 = this%ly*this%ly*dky*dky

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
    do j=1,this%N_y_fft
       do i=1,this%N_x_fft
          this%field1_fft(i,j) = fac*exp(-.25*(lx2kx2(i)+ly2ky2(j)))  
       end do
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
  !
  subroutine rfg2d_fft(this, rseed, rfield, rfield2)

    ! input/output variables
    class(random_fields), intent(inout) :: this ! ffield*_fft is modified
    integer, intent(inout) :: rseed(NRANDSEED) ! nr_ran2 modifies rseed
    real, dimension(this%N_x,this%N_y), intent(out) :: rfield, rfield2

    ! local variables
    !real :: theta, ran_num ! rng
    real,allocatable,dimension(:,:) :: theta, ran_num ! rng
    integer :: N_xy_fft, k ! fft
    integer :: i, j
    integer :: N_x_fft, N_y_fft
    real :: N_xy_fft_real
#ifdef MKL_AVAILABLE
    integer :: mklstat
    complex, allocatable :: z_inout(:)
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
    call this%sqrt_gauss_spectrum_2d

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
    this%field2_fft = sin(theta)*this%field1_fft
    this%field1_fft = cos(theta)*this%field1_fft

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
    mklstat = DftiComputeBackward(this%Desc_Handle, z_inout)
    if (mklstat/= DFTI_NO_ERROR) call quit('DftiComputeBackward failed!')

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

  
end module random_fields_class



#ifdef TEST_RFG2D

program test_rfg2d
  
  use random_fields_class
  use nr_ran2_gasdev
  
  implicit none
  
  integer :: N_x, N_y, i, j, n_e, N_e_tot
  real :: dx, dy, lx, ly, var
  real, allocatable, dimension(:,:) :: field1, field2
  
  character(300) :: file_name
  character(10)  :: n_e_string
  character(100) :: output_format
  character(10)  :: tmp_string

  integer :: RSEEDCONST
  integer, dimension(NRANDSEED) :: rseed
  
  character(5) :: fft_tag
    
  ! instance of random_fields
  type(random_fields) :: rf

  ! start
  RSEEDCONST = -777
  rseed(1) = RSEEDCONST
  write (*,*) RSEEDCONST
  call init_randseed(rseed)

  N_x = 144
  N_y = 91
  dx = 5000.
  dy = 5000.
  lx = 45000.
  ly = 45000.
  var = 1.
  
  
  allocate(field1(N_x,N_y))
  allocate(field2(N_x,N_y))
  
#ifdef MKL_AVAILABLE
  fft_tag = 'mklx.'
#else
  fft_tag = 'nrxx.'
#endif

  ! get N_e fields
  N_e_tot = 10
  do n_e=1,N_e_tot,2

     call rf%initialize(N_x, N_y, var, lx, ly, dx, dy)
     call rf%rfg2d_fft(rseed, field1, field2)
     !call rf%generate_white_field(rseed, field1)
     call rf%finalize
      
     ! write to file
     ! field1
     write(n_e_string,  '(i3.3)') n_e
     file_name = 'rf.'//fft_tag// n_e_string(1:len_trim(n_e_string)) // '.dat'
     write(tmp_string, '(i3.3)') N_y
     output_format = '(' // tmp_string(1:len_trim(tmp_string)) // '(1x,e13.5))'
     open (10,file=file_name,status='unknown')
     do i=1,N_x
        write (10,output_format(1:len_trim(output_format))) (field1(i,j), j=1,N_y)
     end do
     close (10,status='keep')

     ! field2
     write(n_e_string,  '(i3.3)') n_e+1
     file_name = 'rf.' //fft_tag// n_e_string(1:len_trim(n_e_string)) // '.dat'
     write(tmp_string, '(i3.3)') N_y
     output_format = '(' // tmp_string(1:len_trim(tmp_string)) // '(1x,e13.5))'
     open (10,file=file_name,status='unknown')
     do i=1,N_x
        write (10,output_format(1:len_trim(output_format))) (field2(i,j), j=1,N_y)
     end do
     close (10,status='keep')
     
  end do

end program test_rfg2d


#endif


! ======= EOF ==================================================




