
!********************************************************************** 
!
! nr_ran2_gasdev.f90
!
! adapted Numerical Recipes random number generator ran2() and gasdev()
!
! use ran2() instead of rand() and/or ran1() for longer period (~1e18)
!
! gasdev() uses ran2() instead of ran1()
!
! eliminate all save attributes, convert functions into subroutines,
! pass seed state vector into subroutines, enclose in module
!
! reichle, 16 Feb 2005
! reichle, 18 Feb 2005 - make nr_ran2() public, change call statement
! reichle, 24 Mar 2016 - added vectorized version nr_ran2_2d() by Matt Thompson
!
! *********************************************************************

module nr_ran2_gasdev
  
  implicit none
  
  private
  
  public :: NRANDSEED
  public :: nr_ran2_2d
  public :: nr_gasdev
  public :: init_randseed
  
  integer, parameter :: NRANDSEED = 35
  
  integer, parameter :: NTAB = NRANDSEED-3
  integer, parameter :: IM1=2147483563
  integer, parameter :: IM2=2147483399
  real,    parameter :: AM=1./IM1
  integer, parameter :: IMM1=IM1-1
  integer, parameter :: IA1=40014
  integer, parameter :: IA2=40692
  integer, parameter :: IQ1=53668
  integer, parameter :: IQ2=52774
  integer, parameter :: IR1=12211
  integer, parameter :: IR2=3791
  integer, parameter :: NDIV=1+IMM1/NTAB
  real,    parameter :: EPS=1.2e-7
  real,    parameter :: RNMX=1.-EPS
  
  ! RNMX should approximate the largest floating value that is less than 1.
  
contains
  
  !**********************************************************************
  !
  ! ran2()
  ! 
  ! Long period (>2!1e18) random number generator of L Ecuyer with
  ! Bays-Durham shuffle and added safeguards. Returns a uniform
  ! random deviate between 0.0 and 1.0 (exclusive of the endpoint
  ! values). 
  
  subroutine nr_ran2(rseed, ran_num)
    
    implicit none
    
    integer, dimension(NRANDSEED), intent(inout) :: rseed
    
    real, intent(out) :: ran_num
    
    ! local variables
    
    integer :: idum, idum2, iy
    
    integer, dimension(NTAB) :: iv
    
    integer :: j, k
    
    ! -------------------------------------------------------
    
    idum  = rseed(1)           
    idum2 = rseed(2)           
    iy    = rseed(3)           
    iv    = rseed(4:NRANDSEED) 
    
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran_num=min(AM*iy,RNMX)
    
    rseed(1)           = idum
    rseed(2)           = idum2
    rseed(3)           = iy
    rseed(4:NRANDSEED) = iv
    
    return
    
  END subroutine nr_ran2

  !**********************************************************************  
  ! 
  ! vectorized version of nr_ran2(), provided by Matt Thompson
  ! - reichle, 24 Mar 2016

  subroutine nr_ran2_2d(nx, ny, rseed, ran_num)
    
    implicit none
    
    integer, intent(in) :: nx, ny

    integer, dimension(NRANDSEED), intent(inout) :: rseed
    
    real, dimension(nx, ny), intent(out) :: ran_num
    
    ! local variables
    
    integer :: idum, idum2, iy
    
    integer, dimension(NTAB) :: iv
    
    integer :: j, k, icnt, jcnt
    
    ! -------------------------------------------------------

    do jcnt = 1, ny
       do icnt = 1, nx
          
          idum  = rseed(1)           
          idum2 = rseed(2)           
          iy    = rseed(3)           
          iv    = rseed(4:NRANDSEED) 
          
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          k=idum2/IQ2
          idum2=IA2*(idum2-k*IQ2)-k*IR2
          if (idum2.lt.0) idum2=idum2+IM2
          j=1+iy/NDIV
          iy=iv(j)-idum2
          iv(j)=idum
          if(iy.lt.1)iy=iy+IMM1
          ran_num(icnt,jcnt)=min(AM*iy,RNMX)
          
          rseed(1)           = idum
          rseed(2)           = idum2
          rseed(3)           = iy
          rseed(4:NRANDSEED) = iv

       end do
    end do
    
    return
    
  END subroutine nr_ran2_2d
  
  !*************************************************************************
  !
  ! init_randseed()
  !
  ! initialize by calling with negative integer rseed(1)
  ! and fill in the other NRANDSEED-1 integers (stored in idum2, iy, iv)
  
  subroutine init_randseed( rseed )
    
    implicit none
    
    integer, dimension(NRANDSEED), intent(inout) :: rseed
    
    ! local variables
         
    integer :: idum, idum2, iy
    
    integer, dimension(NTAB) :: iv
       
    integer :: j, k
    
    ! ------------------------------------------------------------
    
    idum = rseed(1)
    
    if (idum<=0) then
       idum=max(-idum,1)
       idum2=idum
       do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
       end do
       iy=iv(1)
    else
       write (*,*) 'init_randseed(): initialize by calling with rseed(1)<0'
       write (*,*) 'STOPPING.'
       stop       
    end if
    
    rseed(1)           = idum
    rseed(2)           = idum2
    rseed(3)           = iy
    rseed(4:NRANDSEED) = iv
    
  end subroutine init_randseed
  
  !******************************************************************
  !
  ! gasdev() adapted to use ran2()
  !
  ! Returns TWO normally distributed deviates with zero mean and unit
  ! variance, using ran2() as the source of uniform deviates.  
  !
  ! Use init_randseed() to initialize.
  
  subroutine nr_gasdev(rseed, gasdev) 
    
    implicit none
    
    integer, dimension(NRANDSEED), intent(inout) :: rseed
    
    real, dimension(2), intent(out) :: gasdev
    
    ! local variables
    
    real :: fac, rsq, v1, v2
    
    ! ---------------
    
1   call nr_ran2(rseed, v1)
    call nr_ran2(rseed, v2)
    v1=2.*v1-1.
    v2=2.*v2-1.
    rsq=v1**2+v2**2
    if(rsq.ge.1..or.rsq.eq.0.)goto 1
    fac=sqrt(-2.*log(rsq)/rsq)
    gasdev(1)=v1*fac
    gasdev(2)=v2*fac
    
    return
    
  end subroutine nr_gasdev
  
  ! ************************************************************
  
end module nr_ran2_gasdev
  
! *****************************************************************
!
! driver for testing module nr_ran2_gasdev

#if 0 

program test_my_random_numbers
  
  use nr_ran2_gasdev

  implicit none
  
  integer :: i, RSEEDCONST
  
  integer, dimension(NRANDSEED) :: rseed
  
  real, dimension(2) :: x
  
  ! --------------------------------
  
  RSEEDCONST = -777
  
  rseed(1) = RSEEDCONST
  
  write (*,*) RSEEDCONST
  
  call init_randseed(rseed)
  
  x    = .0
  
  do i=1,10
     
     !x = gasdev(RSEED)
     
     call nr_gasdev(rseed, x)
     
     write (*,*) 'rr ', x(2), rseed(1)
     write (*,*) 'rr ', x(1), rseed(1)
     
  end do
  
end program test_my_random_numbers

#endif


! ******* EOF **********************************************************
