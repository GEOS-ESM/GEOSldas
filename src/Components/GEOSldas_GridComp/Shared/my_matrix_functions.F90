
! this file contains a collection of matrix operation subroutines 
!
! reichle,  1 May 01
! reichle, 18 Apr 06 - renamed subroutine sort()

module my_matrix_functions
  
  implicit none

  private

  public :: row_variance
  public :: row_std
  public :: adjust_mean
  !public :: matrix_std
  public :: unique_rows_3col
  
contains

  ! Changed algorithm to compute (co)variance in subroutines row_covariance() and
  ! row_variance().  Original subroutines commented out below.
  ! wjiang + reichle, 25 Nov 2020
  
  subroutine row_covariance( M, N, A, B, covar )
    ! compute covariance of each row of two M-by-N matrices A and B
    implicit none
    integer, intent(in) :: M, N
    real, intent(in), dimension(M,N) :: A
    real, intent(in), dimension(M,N) :: B
    real, intent(out), dimension(M)  :: covar
    ! locals
    integer :: i, j
    real :: mean_a, mean_b    
    ! -------------------------------------------------------------
    do i=1,M
       mean_a = sum(A(i,:))/N
       mean_b = sum(B(i,:))/N
       covar(i) = 0.
       do j=1,N
          covar(i) = covar(i) + (A(i,j)-mean_a)*(B(i,j)-mean_b)
       end do
    end do
    covar = covar/(N-1)
  end subroutine row_covariance

  ! **********************************************************************
  
  subroutine row_variance( M, N, A, var, mean )
    ! compute variance of each row of an M-by-N matrix A
    ! optionally output mean values
    implicit none
    integer, intent(in) :: M, N
    real, intent(in), dimension(M,N) :: A
    real, intent(out), dimension(M)  :: var
    real, intent(out), dimension(M), optional  :: mean
    ! locals
    integer :: i, j
    real, dimension(M) ::  mean_tmp
    ! -------------------------------------------------------------
    do i=1,M
       mean_tmp(i) = sum(A(i,:))/N
       var(i) = 0
       do j=1,N
          var(i) = var(i) + (A(i,j)-mean_tmp(i))**2
       end do
    end do
    var = var/(N-1)
    if (present(mean))  mean = mean_tmp
  end subroutine row_variance

  ! **********************************************************************

#if 0
  ! The following is a version of subroutine matrix_std() that is
  ! consistent with the revised algorithm for (co)variance computation.
  ! This subroutine is not used.
  ! Note that land_pert.F90 contains another (commented-out) version.
  ! wjiang + reichle, 25 Nov 2020
  
  subroutine matrix_std( N_row, N_col, A, std )
    ! compute std of all elements of N_row by N_col matrix A
    implicit none
    integer, intent(in) :: N_row, N_col
    real, intent(inout), dimension(N_row,N_col)  :: A
    real, intent(out) :: std
    ! ----------------------------
    ! locals
    integer ::  N
    real :: mean
    ! ------------------------------------------------------------
    N = N_row*N_col
    ! compute sample std
    mean = sum(A)/N
    std  = sqrt(sum((A-mean)*(A-mean))/(N-1))
  end subroutine matrix_std

#endif
  
!!  subroutine row_covariance( M, N, A, B, covar )
!!    
!!    ! compute covariance of each row of two M-by-N matrices A and B
!!    
!!    implicit none
!!    
!!    integer, intent(in) :: M, N
!!    
!!    real, intent(in), dimension(M,N) :: A
!!
!!    real, intent(in), dimension(M,N) :: B
!!    
!!    real, intent(out), dimension(M)  :: covar    
!!    
!!    ! locals
!!    
!!    integer :: i, j
!!    
!!    real :: x2, N_real, N_real_minus_one
!!    
!!    ! -------------------------------------------------------------
!!    
!!    N_real = real(N)
!!    
!!    N_real_minus_one = real(N-1)
!!    
!!    do i=1,M
!!       
!!       x2 = 0.0
!!       do j=1,N
!!          x2 = x2 + A(i,j)*B(i,j)
!!       end do
!!       
!!       covar(i) = ( x2 - (sum(A(i,:))*sum(B(i,:)))/N_real )/N_real_minus_one
!!       
!!    end do
!!    
!!  end subroutine row_covariance
!!  
!!  
!!  ! -------------------------------------------------------------------
!!  
!!  subroutine row_variance( M, N, A, var, mean )
!!    
!!    ! compute variance of each row of an M-by-N matrix A
!!    
!!    ! reichle, 16 Jun 2011: added optional output of mean
!!
!!    implicit none
!!    
!!    integer, intent(in) :: M, N
!!    
!!    real, intent(in), dimension(M,N) :: A
!!    
!!    real, intent(out), dimension(M)  :: var    
!!
!!    real, intent(out), dimension(M), optional  :: mean
!!    
!!    ! locals
!!    
!!    integer :: i, j
!!    
!!    real :: x2, N_real, N_real_minus_one
!!
!!    real, dimension(M) :: xm
!!    
!!    ! -------------------------------------------------------------
!!    
!!    N_real = real(N)
!!    
!!    N_real_minus_one = real(N-1)
!!    
!!    do i=1,M
!!       
!!       x2 = 0.0
!!       do j=1,N
!!          x2 = x2 + A(i,j)*A(i,j)
!!       end do
!!
!!       xm(i) = sum(A(i,:))
!!
!!       var(i) = ( x2 - (xm(i)**2)/N_real )/N_real_minus_one
!!       
!!    end do
!!    
!!    ! deal with possible round-off errors
!!    ! reichle, 24 Sep 2004
!!
!!    var = max(var,0.)
!!    
!!    if (present(mean))  mean = xm/N_real
!!    
!!  end subroutine row_variance
  
  ! **********************************************************************

  subroutine row_std( M, N, A, std, mean )
    
    ! compute standard deviation of each row of an M-by-N matrix A
    
    ! reichle,  2 May 2013: added optional output of mean

    implicit none
    
    integer, intent(in) :: M, N
    
    real, intent(in), dimension(M,N) :: A
    
    real, intent(out), dimension(M)  :: std

    real, intent(out), dimension(M), optional  :: mean

    if (present(mean)) then 

       call row_variance( M, N, A, std, mean )
       
    else
       
       call row_variance( M, N, A, std )
       
    end if

    std = sqrt(std)
    
  end subroutine row_std
    
  ! **********************************************************************
  
  subroutine row_third_moment( M, N, A, third_moment )
    
    ! compute third moment of each row of an M-by-N matrix A
    
    ! third_moment = 1/N * sum_{i=1}^N (x_i - mean(x))**3 
    
    implicit none
    
    integer, intent(in) :: M, N
    
    real, intent(in), dimension(M,N) :: A
    
    real, intent(out), dimension(M)  :: third_moment
    
    ! locals
    
    integer :: i, j
    
    real :: x3, mx, N_real
    
    ! -------------------------------------------------------------
    
    N_real = real(N)
    
    do i=1,M
       
       mx = sum(A(i,:))/N_real
       
       x3 = 0.0
       do j=1,N
          x3 = x3 + (A(i,j)-mx)**3
       end do
       
       third_moment(i) = x3/N_real
       
    end do
    
  end subroutine row_third_moment

  
  ! **********************************************************************

#if 0

  ! This subroutine was commented out because it is not currently needed but
  ! contains a call to the obsolete "nr_sort()".  If needed again, replace the 
  ! call to "nr_sort()" with call to MAPL sort routine.
  ! - reichle, 25 Aug 2014
  
  subroutine five_number_summary( M, N, A, five_numbers )
    
    ! get five number summary (median, lower and upper quartiles, min and max)
    ! of each row of an M-by-N data matrix A
    !
    ! inputs:
    !   M  : number of rows of data = number of different data types
    !   N  : number of columns of data = number of ensemble members
    !   A  : M-by-N matrix, left unchanged by this function
    !
    ! outputs:
    !   five_numbers : M-by-5 matrix containing statistical summary:
    !                  column 1: min
    !                  column 2: lower quartile
    !                  column 3: median
    !                  column 4: upper quartile
    !                  column 5: max
    !
    ! Type:   f90
    ! Author: Rolf Reichle
    ! Date:   2 May 2001
    
    implicit none
    
    integer, intent(in) :: M, N
    
    real, intent(in), dimension(M,N)  :: A
    
    real, intent(out), dimension(M,5)  :: five_numbers
    
    ! ----------------------------
    
    ! locals
    
    integer i,d
    
    real, dimension(N) :: tmpvec
    
    do i=1,M
       
       ! put i-th row of data into tmpvec
       
       tmpvec = A(i,:)
       
       ! sort tmpvec in ascending order
       
       call nr_sort( N, tmpvec )
     
       ! get min, max, median and quartiles
       
       ! min and max
       
       five_numbers(i,1) = tmpvec(1)       ! min
       
       five_numbers(i,5) = tmpvec(N)   ! max
       
       ! median
       
       if (mod(N,2) == 0) then
          five_numbers(i,3) = .5*(tmpvec(N/2)+tmpvec(N/2+1))
       else
          five_numbers(i,3) = tmpvec(N/2+1)
       end if
       
       ! quartiles 
       ! (follows Robert Johnson, "Elementary Statistics", PWS-Kent, p69, 1988)
       
       if (mod(N,4) == 0) then
          
          d = N/4
          
          five_numbers(i,2) = .5*(tmpvec(d)+tmpvec(d+1))             ! lower
          
          five_numbers(i,4) = .5*(tmpvec(N-d)+tmpvec(N-d+1)) ! upper
          
       else
          
          d = N/4+1         
          
          five_numbers(i,2) = tmpvec(d)           ! lower
          
          five_numbers(i,4) = tmpvec(N-d)     ! upper
          
       end if
       
    end do
    
  end subroutine five_number_summary

#endif

  ! ------------------------------------------------------------------
  
  subroutine adjust_mean( N_row, N_col, A, M )
    
    ! adjust N_row by N_col matrix A such that 
    ! mean over columns for each row is given by the
    ! corresponding element in vector M of length N_row
    ! 
    ! vector of mean values M is optional input, if not present 
    ! zero mean is assumed
    
    implicit none
    
    integer, intent(in) :: N_row, N_col
    
    real, intent(inout), dimension(N_row,N_col)  :: A
    
    real, intent(in), optional, dimension(N_row) :: M
    
    ! ----------------------------
    
    ! locals
    
    integer i
    
    real, dimension(N_row) :: correction
    
    ! ------------------------------------------------------------
    
    if (present(M)) then
       correction = M - sum(A,2)/real(N_col) 
    else
       correction = - sum(A,2)/real(N_col) 
    end if
    
    do i=1,N_col
       A(:,i) = A(:,i) + correction
    end do
    
  end subroutine adjust_mean
  
  ! ------------------------------------------------------------------

#if 0
  ! This subroutine is not used.
  ! Note that land_pert.F90 contains another (commented-out) version.
  ! wjiang + reichle, 25 Nov 2020

  subroutine adjust_std( N_row, N_col, A, std )
    
    ! adjust N_row by N_col matrix A such that (sample) standard deviation
    !  of all elements is exactly equal to std
    ! 
    ! std is optional input, if not present std=1 is assumed
    
    implicit none
    
    integer, intent(in) :: N_row, N_col
    
    real, intent(inout), dimension(N_row,N_col)  :: A
    
    real, intent(in), optional :: std
    
    ! ----------------------------
    
    ! locals
    
    integer :: i, j
    
    real :: correction, sample_std
    
    ! ------------------------------------------------------------
    
    ! compute sample std

    call matrix_std( N_row, N_col, A, sample_std )
    
    if (present(std)) then
       correction = std/sample_std
    else
       correction = 1./sample_std
    end if
    
    do i=1,N_row
       do j=1,N_col
          A(i,j) = correction*A(i,j)
       end do
    end do
    
  end subroutine adjust_std
#endif
  
  ! ------------------------------------------------------------------
  
!!  subroutine matrix_std( N_row, N_col, A, std )
!!    
!!    ! compute std of all elements of N_row by N_col matrix A
!!    
!!    implicit none
!!    
!!    integer, intent(in) :: N_row, N_col
!!    
!!    real, intent(inout), dimension(N_row,N_col)  :: A
!!    
!!    real, intent(out) :: std
!!    
!!    ! ----------------------------
!!    
!!    ! locals
!!    
!!    integer :: i, j
!!    
!!    real :: x2, m, N_real, N_real_minus_one
!!    
!!    ! ------------------------------------------------------------
!!    
!!    N_real = real(N_row)*real(N_col)
!!    
!!    N_real_minus_one = N_real - 1.
!!    
!!    ! compute sample std
!!    
!!    x2 = 0.0
!!    m  = 0.0
!!    
!!    do i=1,N_row
!!       do j=1,N_col
!!          m  = m  + A(i,j)
!!          x2 = x2 + A(i,j)*A(i,j)
!!       end do
!!    end do
!!    
!!    std = sqrt( ( x2 - m**2/N_real )/N_real_minus_one )
!!        
!!  end subroutine matrix_std


  ! ****************************************************************
  
  subroutine unique_rows_2col( N_rows, A, N_unique_rows, ind_A2U )
    
    ! Identify unique rows in 2-column matrix A (N_rows-by-2).  Unique rows 
    ! are returned in the first (N_unique_rows) of matrix A.  Also returned 
    ! is the index vector from the original row indices of A to the 
    ! unique rows stored in the returned A (upon return, only the first 
    ! N_unique_rows of A are meaningful).
    !
    ! uses nr_indexx.f
    !
    ! reichle, 14 Jun 2012
    ! reichle, 25 Oct 2012: fixed case N_rows=0 upon input
    !
    ! --------------------------------------------------------------------
    
    implicit none
    
    integer,                           intent(in)    :: N_rows
    
    real,    dimension(N_rows,2),      intent(inout) :: A
    
    integer,                           intent(out)   :: N_unique_rows
    
    integer, dimension(N_rows),        intent(out)   :: ind_A2U
    
    ! local variables
    
    real,    parameter         :: tol_frac = 1.e-5
    
    integer                    :: i, iunique, istart, iend, N_tmp
    
    real                       :: tol1, tol2
    
    integer, dimension(N_rows) :: indx, indx2
    
    ! --------------------------------------------------------------------

    if (N_rows==0) then

       N_unique_rows = 0

       return   ! nothing else left to do
              
    end if
        
    ! -------------------------------------
    !
    ! tolerances for check of (in)equality of real numbers
    
    tol1 = tol_frac * sum(abs(A(:,1))) / real(N_rows)
    tol2 = tol_frac * sum(abs(A(:,2))) / real(N_rows)
    
    ! -------------------------------------
    !
    ! sort A according to *first* column   
    
    call nr_indexx(N_rows, A(:,1), indx)
    
    ! apply sort
    
    A = A(indx, :)
    
    !write (*,*) 'after first column sort'
    !do j=1,N_rows
    !   write (*,*) A(j,:)
    !end do
    
    ! -------------------------------------
    !
    ! for each block of identical entries in sorted A(:,1), sort according 
    ! to *second* column
    
    istart = 1  ! start counter for block of identical entries in sorted A(:,1)
    
    do i=1,(N_rows-1)
       
       ! compare pairs of subsequent elements in first column
       
       if ( (abs(A(i,1)-A(i+1,1))>tol1) .or. (i==(N_rows-1)) ) then
          
          ! reached new block [ie, A(i,1)/=A(i+1,1)] or 
          ! reached final pair of elements 
          
          iend = i
          
          ! special treatment for final pair of elements: 
          ! if elements are identical, include final row in block to be sorted,
          ! ignore otherwise (no need to sort final row if distinct from others)
          
          if ( (abs(A(i,1)-A(i+1,1))<tol1) .and. (i==(N_rows-1)) ) then
             
             iend=i+1   ! because A(N_rows-1,1)==A(N_rows,1)
             
          end if
          
          ! within block (istart:iend) sort according to second column
          
          N_tmp = iend - istart + 1
          
          if (N_tmp>1) then   ! sort only if more than one element in block
             
             call nr_indexx(N_tmp, A(istart:iend,2), indx2(1:N_tmp))
             
             ! apply sort 
             
             A(   istart:iend,:) = A(   indx2(1:N_tmp)+istart-1,:)
             indx(istart:iend  ) = indx(indx2(1:N_tmp)+istart-1  )
             
             !write (*,*) 'after second column sort'
             !write (*,*) i, N_tmp
             !do j=1,N_rows
             !   write (*,*) A(j,:)
             !end do
             
          end if
          
          ! re-init
          
          istart = i+1
          
       end if
       
    end do
    
    ! -------------------------------------
    !  
    ! eliminate identical rows
    
    iunique = 1
    
    do i=1,(N_rows-1)
       
       ! record mapping from original row indices to unique rows
       
       ind_A2U(indx(i)) = iunique
       
       ! special treatment for last pair of elements
       
       if ( (i==(N_rows-1))               .and.            &
            (abs(A(i,1)-A(i+1,1))<tol1)   .and.            &
            (abs(A(i,2)-A(i+1,2))<tol2)          ) then
          
          ind_A2U(indx(i+1)) = iunique       ! last two rows identical
          
       else    
          
          ind_A2U(indx(i+1)) = iunique+1     ! last two rows different
          
       end if
       
       ! compare pairs of subsequent rows
       
       if ( (abs(A(i,1)-A(i+1,1))>tol1)   .or.             &
            (abs(A(i,2)-A(i+1,2))>tol2)          ) then
          
          ! found "new" row, pull forward
          
          iunique = iunique+1
          
          A(iunique,:) = A(i+1,:)
                    
       end if
       
    end do
    
    N_unique_rows = iunique
    
  end subroutine unique_rows_2col

  ! ****************************************************************
  
  subroutine unique_rows_3col( N_rows, A, N_unique_rows, ind_A2U )
    
    ! Identify unique rows in 3-column matrix A (N_rows-by-3).  Unique rows 
    ! are returned in the first (N_unique_rows) of matrix A.  Also returned 
    ! is the index vector from the original row indices of A to the 
    ! unique rows stored in the returned A (upon return, only the first 
    ! N_unique_rows of A are meaningful).
    !
    ! uses subroutine unique_rows_2col()
    !
    ! equivalent MATLAB code {Version: 7.14.0.739 (R2012a)}:
    !
    !   [ A_unique, ind_tmp, ind_A2U ] = unique( A_orig, 'rows');
    !
    ! reichle, 31 Mar 2015
    !
    ! --------------------------------------------------------------------
    
    implicit none
    
    integer,                           intent(in)    :: N_rows
    
    real,    dimension(N_rows,3),      intent(inout) :: A
    
    integer,                           intent(out)   :: N_unique_rows
    
    integer, dimension(N_rows),        intent(out)   :: ind_A2U
    
    ! local variables

    integer                         :: N_step1,       N_step2

    integer, dimension(N_rows)      :: ind_A2U_step1, ind_A2U_step2  

    real,    dimension(N_rows,2)    :: A_step1,       A_step2
    
    ! --------------------------------------------------------------------
    
    if (N_rows==0) then
       
       N_unique_rows = 0
       
       return   ! nothing else left to do
       
    end if
    
    ! Step 1:
    ! determine unique rows of submatrix A(:,1:2)
    
    A_step1       = A(:,1:2)
    
    call unique_rows_2col( N_rows, A_step1, N_step1, ind_A2U_step1 )

    ! Step 2:
    ! assemble temporary 2-column matrix with indicator of unique
    ! rows of A(:,1:2) in column 1 and A(:,3) in column 2
    
    A_step2(:,1)  = real( ind_A2U_step1 )
    A_step2(:,2)  = A(:,3)
    
    call unique_rows_2col( N_rows, A_step2, N_step2, ind_A2U_step2 )
    
    ! Step 3:
    ! finalize

    N_unique_rows = N_step2

    ind_A2U       = ind_A2U_step2

    A(:,1:2)      = A_step1(nint(A_step2(1:N_step2,1)),:)
    A(:,  3)      = A_step2(             1:N_step2    ,2)
    
  end subroutine unique_rows_3col

  ! *****************************************

end module my_matrix_functions

! ************************************************************************

! driver programs for testing

#if 0

program test_adjust_std
  
  use my_matrix_functions
  
  implicit none
  
  integer, parameter :: M=3, N=30
  
  integer :: i, j

  real, dimension(M,N) :: A
  
  real :: std
 
  call random_number( A )

  call matrix_std( M, N, A, std)
  write (*,*) std
  
  do i=1,M
     
     write (*,'(30(e13.5))') (A(i,j), j=1,N)
     
  end do
  
  call adjust_std(M,N,A,9.999)

  call matrix_std( M, N, A, std)
  write (*,*) std

  write (*,*)
  
  do i=1,M
     
     write (*,'(30(e13.5))') (A(i,j), j=1,N)
     
  end do
  
  call adjust_std(M,N,A)

  call matrix_std( M, N, A, std)
  write (*,*) std  

  write (*,*)
  
  do i=1,M
     
     write (*,'(30(e13.5))') (A(i,j), j=1,N)
     
  end do
  
end program test_adjust_std

#endif

! ---------------------------------------

#if 0

program test_adjust_mean
  
  use my_matrix_functions
  
  implicit none
  
  integer, parameter :: M=3, N=10
  
  integer :: i, j

  real, dimension(M,N) :: A
  
  real, dimension(M) :: mean_values
  
  call random_number( A )
  
  call random_number(mean_values)

  do i=1,M
     
     write (*,'(30(e13.5))') (A(i,j), j=1,N)
     
  end do

  do i=1,M
     
     write (*,'(30(e13.5))') mean_values(i)
     
  end do

  call adjust_mean(M,N,A,mean_values)
  
  write (*,*)
  
  do i=1,M
     
     write (*,'(30(e13.5))') (A(i,j), j=1,N)
     
  end do
  
  call adjust_mean(M,N,A)
  
  write (*,*)
  
  do i=1,M
     
     write (*,'(30(e13.5))') (A(i,j), j=1,N)
     
  end do
  

end program test_adjust_mean

#endif

! driver routine for testing

#if 0

program test_row_variance

  use my_matrix_functions
  
  implicit none
  
  integer, parameter :: M=3, N=10

  integer :: i, j

  real, dimension(M,N) :: A
  
  real, dimension(M) :: var
  
  call random_number( A )
  
   
  do i=1,M
     
     write (*,'(30(e12.5))') (A(i,j), j=1,N)
     
  end do
  
  call row_variance(M,N,A,var)

  write (*,*)
  
  do i=1,M
     
     write (*,'(1(e12.5))') var(i)
     
  end do
  

end program test_row_variance

#endif

! driver routine for testing

#if 0

program test_row_covariance
  
  use my_matrix_functions
  
  implicit none
  
  integer, parameter :: M=5, N=10
  
  integer :: i, j
  
  real, dimension(M,N) :: A, B
  
  real, dimension(M) :: covar
  
  call random_number( A )
  call random_number( B )

  write (*,*) 'A=['

  do i=1,M
     
     write (*,'(30(e12.5))') (A(i,j), j=1,N)
     
  end do
  
  write (*,*) ']'

  write (*,*)

  write (*,*) 'B=['

  do i=1,M
     
     write (*,'(30(e12.5))') (B(i,j), j=1,N)
     
  end do

  write (*,*) ']'

  call row_covariance(M,N,A,B,covar)
  
  write (*,*)
  
  do i=1,M
     
     write (*,'(1(e12.5))') covar(i)
     
  end do
  
  
end program test_row_covariance

#endif


#if 0

program test_five_number_summary
  
  use my_matrix_functions
  
  implicit none
  
  integer :: i,j
  
  integer, parameter :: N_data = 10, N_ens = 30
  
  real, dimension(N_data,N_ens) :: my_data
  
  real, dimension(N_data,5) :: five_numbers
  
  call random_number(my_data)
  
  do i=1,N_data
     
     write (*,'(30(e12.5))') (my_data(i,j), j=1,N_ens)
     
  end do

  call five_number_summary(N_data, N_ens, my_data, five_numbers)
  
  write (*,*) size(five_numbers)
  
  do i=1,N_data
     
     write (*,'(5(e12.5))') (five_numbers(i,j), j=1,5)
     
  end do
  

end program test_five_number_summary
  
#endif


#if 0

program test_unique_rows
  
  implicit none
  
  integer, parameter :: N_rows = 55    ! 6
  integer, parameter :: N_cols = 2
  
  real, dimension(N_rows,N_cols) :: A, A_orig
  
  integer :: N_unique_rows, i, j
  
  integer, dimension(N_rows) :: ind_A2U
  
  ! ---------------------------------------------
  
  A(1,:) = (/ 1, 2 /)
  A(2,:) = (/ 2, 2 /)
  A(3,:) = (/ 3, 2 /)
  A(4,:) = (/ 1, 2 /)
  A(5,:) = (/ 2, 2 /)
  A(6,:) = (/ 1, 2 /)

  if (N_rows==55) then
    
     A( 1,:) = (/     0.20, 32.5 /)
     A( 2,:) = (/     1.41, 37.5 /)
     A( 3,:) = (/     1.41, 42.5 /)
     A( 4,:) = (/     1.41, 47.5 /)
     A( 5,:) = (/     1.41, 52.5 /)
     A( 6,:) = (/     1.41, 57.5 /)
     A( 7,:) = (/     1.41, 32.5 /)
     A( 8,:) = (/     1.41, 37.5 /)
     A( 9,:) = (/     1.41, 42.5 /)
     A(10,:) = (/     1.41, 47.5 /)
     A(11,:) = (/     1.41, 52.5 /)
     A(12,:) = (/     1.41, 57.5 /)
     A(13,:) = (/     1.41, 32.5 /)
     A(14,:) = (/     1.41, 37.5 /)
     A(15,:) = (/     1.41, 42.5 /)
     A(16,:) = (/     1.41, 47.5 /)
     A(17,:) = (/     1.41, 52.5 /)
     A(18,:) = (/     1.41, 57.5 /)
     A(19,:) = (/     1.41, 32.5 /)
     A(20,:) = (/     1.41, 37.5 /)
     A(21,:) = (/    23.00, 42.5 /)
     A(22,:) = (/     1.41, 47.5 /)
     A(23,:) = (/     1.41, 52.5 /)
     A(24,:) = (/     1.41, 57.5 /)
     A(25,:) = (/    18.00, 32.5 /)
     A(26,:) = (/    18.00, 42.5 /)
     A(27,:) = (/    18.00, 52.5 /)
     A(28,:) = (/    18.00, 32.5 /)
     A(29,:) = (/    18.00, 42.5 /)
     A(30,:) = (/    18.00, 52.5 /)
     A(31,:) = (/     6.70, 32.5 /)
     A(32,:) = (/     6.70, 37.5 /)
     A(33,:) = (/     6.70, 42.5 /)
     A(34,:) = (/     6.70, 47.5 /)
     A(35,:) = (/     6.70, 52.5 /)
     A(36,:) = (/     6.70, 57.5 /)
     A(37,:) = (/     6.70, 32.5 /)
     A(38,:) = (/     6.70, 37.5 /)
     A(39,:) = (/     6.70, 42.5 /)
     A(40,:) = (/     6.70, 47.5 /)
     A(41,:) = (/     6.70, 52.5 /)
     A(42,:) = (/     6.70, 57.5 /)
     A(43,:) = (/     6.70, 32.5 /)
     A(44,:) = (/     6.70, 37.5 /)
     A(45,:) = (/     6.70, 42.5 /)
     A(46,:) = (/     6.70, 47.5 /)
     A(47,:) = (/     6.70, 52.5 /)
     A(48,:) = (/     6.70, 57.5 /)
     A(49,:) = (/     6.70, 32.5 /)
     A(50,:) = (/     6.70, 37.5 /)
     A(51,:) = (/     6.70, 42.5 /)
     A(52,:) = (/     6.70, 47.5 /)
     A(53,:) = (/     6.70, 52.5 /)
     A(54,:) = (/     6.70, 57.5 /)     
     A(55,:) = (/    87.00, 57.5 /)
     
     
  end if
  
  ! ------------------------------------------------
  
  A_orig = A

  do i=1,N_rows
     write (*,*) A(i,:)
  end do
  
  call unique_rows_2col( N_rows, N_cols, A, N_unique_rows, ind_A2U)
  
  write (*,*) N_unique_rows
  do i=1,N_unique_rows
     write (*,*) A(i,:)
  end do
  do i=1,N_rows
     write (*,*) A_orig(i,:), ind_A2U(i)
  end do
  
end program test_unique_rows

#endif


#if 0

program test_unique_rows_3col
  
  ! find unique rows in 3-column matrix
  !
  ! - reichle, 31 Mar 2015
  
  use my_matrix_functions, ONLY: unique_rows_3col
  
  implicit none
  
  integer, parameter :: N_rows = 55   !9  
  integer, parameter :: N_cols = 3
  
  integer            :: N_unique_rows, i

  real,    dimension(N_rows,N_cols) :: A, A_orig

  logical, dimension(N_rows,N_cols) :: final_check
  
  integer, dimension(N_rows)        :: ind_A2U
  
  ! ---------------------------------------------
  
  A(1,:) = (/ 1.5, 7.5, 4.5 /)
  A(2,:) = (/ 2.5, 7.5, 6.5 /)
  A(3,:) = (/ 3.5, 7.5, 4.5 /)
  A(4,:) = (/ 1.5, 3.5, 6.5 /)
  A(5,:) = (/ 2.5, 3.5, 6.5 /)
  A(6,:) = (/ 3.5, 7.5, 4.5 /)
  A(7,:) = (/ 1.5, 3.5, 6.5 /)
  A(8,:) = (/ 2.5, 7.5, 6.5 /)
  A(9,:) = (/ 3.5, 7.5, 6.5 /)

  if (N_rows==55) then
     
     A( 1,:) = (/     0.20, 32.5,  0 /)
     A( 2,:) = (/     1.41, 37.5, 19 /)
     A( 3,:) = (/     1.41, 42.5, 19 /)
     A( 4,:) = (/     1.41, 47.5, 19 /)
     A( 5,:) = (/     1.41, 52.5,  0 /)
     A( 6,:) = (/     1.41, 57.5,  0 /)
     A( 7,:) = (/     1.41, 32.5,  0 /)
     A( 8,:) = (/     1.41, 37.5,  0 /)
     A( 9,:) = (/     1.41, 42.5,  0 /)
     A(10,:) = (/     1.41, 47.5,  0 /)
     A(11,:) = (/     1.41, 52.5,  0 /)
     A(12,:) = (/     1.41, 57.5,  0 /)
     A(13,:) = (/     1.41, 32.5,  0 /)
     A(14,:) = (/     1.41, 37.5,  0 /)
     A(15,:) = (/     1.41, 42.5,  0 /)
     A(16,:) = (/     1.41, 47.5,  0 /)
     A(17,:) = (/     1.41, 52.5,  0 /)
     A(18,:) = (/     1.41, 57.5,  0 /)
     A(19,:) = (/     1.41, 32.5,  0 /)
     A(20,:) = (/     1.41, 37.5,  0 /)
     A(21,:) = (/    23.00, 42.5,  0 /)
     A(22,:) = (/     1.41, 47.5,  0 /)
     A(23,:) = (/     1.41, 52.5,  0 /)
     A(24,:) = (/     1.41, 57.5,  0 /)
     A(25,:) = (/    18.00, 32.5,  0 /)
     A(26,:) = (/    18.00, 42.5,  0 /)
     A(27,:) = (/    18.00, 52.5,  0 /)
     A(28,:) = (/    18.00, 32.5,  0 /)
     A(29,:) = (/    18.00, 42.5,  0 /)
     A(30,:) = (/    18.00, 52.5,  0 /)
     A(31,:) = (/     6.70, 32.5,  3 /)
     A(32,:) = (/     6.70, 37.5,  3 /)
     A(33,:) = (/     6.70, 42.5,  3 /)
     A(34,:) = (/     6.70, 47.5,  0 /)
     A(35,:) = (/     6.70, 52.5,  3 /)
     A(36,:) = (/     6.70, 57.5,  3 /)
     A(37,:) = (/     6.70, 32.5,  3 /)
     A(38,:) = (/     6.70, 37.5,  3 /)
     A(39,:) = (/     6.70, 42.5,  3 /)
     A(40,:) = (/     6.70, 47.5,  3 /)
     A(41,:) = (/     6.70, 52.5,  3 /)
     A(42,:) = (/     6.70, 57.5,  3 /)
     A(43,:) = (/     6.70, 32.5,  0 /)
     A(44,:) = (/     6.70, 37.5,  0 /)
     A(45,:) = (/     6.70, 42.5,  0 /)
     A(46,:) = (/     6.70, 47.5,  0 /)
     A(47,:) = (/     6.70, 52.5,  0 /)
     A(48,:) = (/     6.70, 57.5,  0 /)
     A(49,:) = (/     6.70, 32.5,  0 /)
     A(50,:) = (/     6.70, 37.5,  0 /)
     A(51,:) = (/     6.70, 42.5,  0 /)
     A(52,:) = (/     6.70, 47.5,  0 /)
     A(53,:) = (/     6.70, 52.5,  0 /)
     A(54,:) = (/     6.70, 57.5,  0 /)     
     A(55,:) = (/    87.00, 57.5,  0 /)
     
  end if

  ! ------------------------------------------------

  A_orig = A

  do i=1,N_rows
     write (*,*) A(i,:)
  end do
    
  call unique_rows_3col( N_rows, A, N_unique_rows, ind_A2U )
  
  write (*,*)   
  write (*,*) N_unique_rows
  do i=1,N_unique_rows
     write (*,*) A(i,:)                     ! Tb_FreqAngRTMid
  end do
  
  write (*,*) 
  do i=1,N_rows
     write (*,*) A_orig(i,:), ind_A2U(i)    ! ind_Tbspecies2TbuniqFreqAngRTMid
  end do
  
  ! final check
  
  final_check = ( A_orig == A(ind_A2U,:) )
  
  write (*,*)
  do i=1,N_rows
     write (*,*) final_check(i,:) 
  end do
  
  write (*,*) 
  if (all(final_check)) then
     
     write (*,*) 'success'
     
  else
     
     write (*,*) 'ERROR - failed final check'

  end if
    
end program test_unique_rows_3col

#endif

! ***** EOF **************************************************************
