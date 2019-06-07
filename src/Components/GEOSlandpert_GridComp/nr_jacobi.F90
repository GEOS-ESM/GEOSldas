
module nr_jacobi

  ! Jacobi transformation
  !
  ! adapted from f77 Numerical Recipes by reichle, 7 Jun 2005
  
  implicit none
  
  private
  
  public :: jacobi

contains
  
  subroutine jacobi(amat,n,d,v,nrot)
    
    ! Computes all eigenvalues and eigenvectors of a real symmetric matrix 
    ! amat, which is of size n by n. The vector d returns the 
    ! eigenvalues of amat in its first n elements. v is a matrix with the same
    ! dimensions as amat, whose columns contain, on output, 
    ! the normalized eigenvectors of amat. nrot returns the number of Jacobi 
    ! rotations that were required.
    
    ! The Jacobi method is absolutely foolproof for all real symmetric 
    ! matrices. For matrices of order greater than about 10, say, the 
    ! algorithm is slower, by a significant constant factor, than the QR 
    ! method. However, the Jacobi algorithm is much simpler than 
    ! the more efficient methods. We thus recommend it for matrices of 
    ! moderate order, where expense is not a major consideration.
    
    ! Eigenvector decomposition:
    !
    !    Amat = V*diag(D)*transpose(V)
    
    integer,                   intent(in)               :: n
    real,    dimension(n,n),   intent(in)               :: amat
    real,    dimension(n,n),   intent(out)              :: v
    real,    dimension(n),     intent(out)              :: d
    integer,                   intent(out), optional    :: nrot
    
    ! local variables
    
    integer :: i,ip,iq,j
    
    real :: c, g, h, s, sm, t, tau, theta, tresh
    
    real, dimension(n,n) :: a
    real, dimension(n)   :: b, z
    
    ! ------------------------------------------------------------
    !
    ! make sure this is not used for large matrices
    ! (comment out this block if you want to apply nr_jacobi to
    !  large matrices)
    
    if (n>100) then
       write (*,*) 'nr_jacobi(): Jacobi not efficient for large matrices'
       write (*,*) 'use of QR factorization suggested'
       write (*,*) 'STOPPING.'
       stop
    end if
       
    ! initialize temporary matrix a 
    
    a = amat
    
    ! Initialize to the identity matrix.
    
    do ip=1,n
       do iq=1,n
          v(ip,iq)=0.
       end do
       v(ip,ip)=1.
    end do
    
    ! Initialize b and d to the diagonal of a.
    ! The vector z will accumulate terms of the form ta(p,q)
    !  as in equation (11.1.14). 

    do ip=1,n
       b(ip)=a(ip,ip)
       d(ip)=b(ip)
       z(ip)=0.
    end do
    if (present(nrot)) nrot=0
    do i=1,50
       sm=0.
       do ip=1,n-1  ! Sum off-diagonal elements.
          do iq=ip+1,n
             sm=sm+abs(a(ip,iq))
          end do
       end do

       ! normal return (relies on quadratic convergence to machine underflow. 
       
       if(sm.eq.0.)  &
            return 
       
       if(i.lt.4)then
          tresh=0.2*sm/n**2  ! ...on first three sweeps
       else
          tresh=0.           ! ...thereafter
       end if
       do ip=1,n-1
          do iq=ip+1,n
             g=100.*abs(a(ip,iq))
             if((i.gt.4).and.(abs(d(ip))+ &
                  g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
                a(ip,iq)=0.
             else if (abs(a(ip,iq)).gt.tresh)then
                h=d(iq)-d(ip)
                if(abs(h)+g.eq.abs(h))then
                   t=a(ip,iq)/h
                else
                   theta=0.5*h/a(ip,iq)
                   t=1./(abs(theta)+sqrt(1.+theta**2))
                   if(theta.lt.0.)t=-t
                end if
                c=1./sqrt(1+t**2)
                s=t*c
                tau=s/(1.+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0.
                do j=1,ip-1
                   g=a(j,ip)
                   h=a(j,iq)
                   a(j,ip)=g-s*(h+g*tau)
                   a(j,iq)=h+s*(g-h*tau)
                end do
                do j=ip+1,iq-1
                   g=a(ip,j)
                   h=a(j,iq)
                   a(ip,j)=g-s*(h+g*tau)
                   a(j,iq)=h+s*(g-h*tau)
                end do
                do j=iq+1,n
                   g=a(ip,j)
                   h=a(iq,j)
                   a(ip,j)=g-s*(h+g*tau)
                   a(iq,j)=h+s*(g-h*tau)
                end do
                do j=1,n
                   g=v(j,ip)
                   h=v(j,iq)
                   v(j,ip)=g-s*(h+g*tau)
                   v(j,iq)=h+s*(g-h*tau)
                end do
                if (present(nrot)) nrot=nrot+1
             end if
          end do
       end do
       do ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
       end do
    end do
    pause 'too many iterations in jacobi'
    return
    
  end subroutine jacobi

end module nr_jacobi

! ***********************************************************

#if 0

program test
  
  use nr_jacobi
  
  implicit none
  
  integer, parameter :: N = 3

  real, dimension(N,N) :: A, V
  
  real, dimension(N)   :: D
  
  integer :: i,j,k

  character(8)  :: date_string
  character(10) :: time_string

  ! ------------------------------------------------
  
  A(1,1) = 1.
  A(1,2) = -.8
  A(1,3) =  .5
  A(2,2) = 1.
  A(2,3) = -.5
  A(3,3) = 1.
  A(2,1) = A(1,2)
  A(3,2) = A(2,3)
  A(3,1) = A(1,3)
  
  do i=1,N
     write (*,*) A(i,:)     
  end do

  write (*,*)

  call date_and_time(date_string, time_string)
    
  write (*,*) 'started looping at ', date_string, ', ', time_string

  do k=1,1000000
     
     call jacobi(A,N,D,V)

  end do

  call date_and_time(date_string, time_string)
    
  write (*,*) 'ended   looping at ', date_string, ', ', time_string

  do i=1,N
     write (*,*) A(i,:)     
  end do

  write (*,*)  

  do i=1,N
     write (*,*) V(i,:)     
  end do

  write (*,*)

  write (*,*) D
    
  write (*,*)
  
end program test
  
#endif

! ============== EOF =========================================
