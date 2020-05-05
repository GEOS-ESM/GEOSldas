
! this file contains a collection of general Ensemble Kalman filter
! subroutines and compact support subroutines
!
! reichle,      20 Apr 2001
! reichle,      18 Mar 2004 - optional arguments
! reichle,      27 Jan 2005 - eliminated use of module select_kinds
! reichle,      19 Jul 2005 - merged compact_support.f90 and enkf_general.f90
! reichle,       1 Aug 2005 - eliminated tile_coord
! reichle,      18 Oct 2005 - return increments instead of updated State
! reichle+qliu, 29 Apr 2020 - added forecast error covariance inflation

! use intel mkl lapack when available
#ifdef MKL_AVAILABLE
#include "lapack.f90"
#endif

module enkf_general
  
#ifdef MKL_AVAILABLE
  use lapack95, only: getrf, getrs
#endif   

  use enkf_types,                      ONLY:      &
       obs_type

  use LDAS_ExceptionsMod,              ONLY:      &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR
    
  implicit none
  
  private
  
  public :: enkf_increments
  
contains
  
  ! *** new version for 3D EnKF ***
  
  subroutine enkf_increments( &
       N_state, N_obs, N_ens, &
       Observations, Obs_pred, Obs_err, Obs_cov, &
       State_incr, &
       State_lon, State_lat, xcompact, ycompact, fcsterr_inflation_fac )
    
    ! perform EnKF update
    !
    ! IMPORTANT:
    ! on input, State_incr must contain State_minus(1:N_state,1:N_ens)
    ! on output, State_incr contains the increments
    !
    ! if optional inputs State_lon, State_lat, xcompact, and ycompact
    ! are present, Hadamard product is applied to HPHt and PHt
    !
    ! if optional input fcsterr_inflation_fac is present, this subroutine inflates
    ! the forecast error covariance and returns increments that must be applied to 
    ! the state vector *before* inflation
    ! (i.e., do not inflate the state vector outside of this subroutine)
    
    implicit none
    
    integer, intent(in) :: N_state, N_obs, N_ens
    
    type(obs_type), intent(in), dimension(N_obs) :: Observations 
    
    real, intent(in), dimension(N_obs,N_ens) :: Obs_pred
    real, intent(in), dimension(N_obs,N_ens) :: Obs_err
    real, intent(in), dimension(N_obs,N_obs) :: Obs_cov    
    
    real, intent(inout), dimension(N_state,N_ens) :: State_incr
    
    ! optional inputs

    real, dimension(N_state), intent(in), optional :: State_lon, State_lat 
    
    real, intent(in), optional :: xcompact               ! [deg] longitude
    real, intent(in), optional :: ycompact               ! [deg] latitude
    real, intent(in), optional :: fcsterr_inflation_fac  ! forecast error covariance inflation 

    ! -----------------------------
    
    ! locals

    character(len=*), parameter :: Iam = 'enkf_increments'
    
    integer :: n_e, i, ii, jj, kk, lapack_info

    real :: PHt_ij, dx, dy

    real :: inflation_factor
    
    real, dimension(N_state,N_ens) :: State_prime
    real, dimension(N_state)       :: State_bar
    real, dimension(N_state)       :: State_incr_tmp
    
    real, dimension(N_obs,N_ens)   :: Obs_pred_prime
    real, dimension(N_obs)         :: Obs_pred_bar
    real, dimension(N_obs)         :: rhs
    
    real, dimension(N_ens)         :: weights    
    
    real, dimension(N_obs,N_obs)   :: Repr_matrix
    
    integer,          dimension(N_obs)         :: indx

    logical :: apply_hadamard
    
    ! ------------------------------------------------------------------

    ! deal with optional argument
    
    if (present(fcsterr_inflation_fac)) then

       inflation_factor = fcsterr_inflation_fac

    else

       inflation_factor = -9999.

    end if

    ! find out whether Hadamard product should be applied
    
    apply_hadamard = (               &
         present(State_lon)    .and. &
         present(State_lat)    .and. &
         present(xcompact)     .and. &
         present(ycompact)              )
    
    ! ----------------------
    
    ! IMPORTANT: on input, State_incr contains State_minus(1:N_state,1:N_ens)
    
    ! compute ensemble mean Ybar at current update time
    
    State_bar = sum( State_incr, 2) / real(N_ens)
    
    ! finalize matrix Y_prime = Y - Ybar
    
    do n_e=1,N_ens
       
       State_prime(:,n_e) = State_incr(:,n_e) - State_bar
       
    end do

    if (inflation_factor > 0.)  State_prime = inflation_factor * State_prime        

    ! --------------------
    
    ! compute ensemble mean H*Ybar
    
    Obs_pred_bar = sum( Obs_pred, 2) / real(N_ens)
    
    ! finalize matrix Q_prime = H*Y - H*ybar
    
    do n_e=1,N_ens
       
       Obs_pred_prime(:,n_e) = Obs_pred(:,n_e) - Obs_pred_bar
       
    end do

    if (inflation_factor > 0.)  Obs_pred_prime = inflation_factor * Obs_pred_prime

    ! --------------------

    ! form repr matrix HPHt = Q_prime*(Q_prime)t/(N_e-1)
    
    Repr_matrix = &
         (matmul(Obs_pred_prime,transpose(Obs_pred_prime))) &
         /real(N_ens-1)

    ! reichle, 18 Mar 2004:
    ! maybe Hadamard product should be applied *after* adding Obs_cov
    ! to representer matrix? only matters if Obs_cov is not diagonal...
    
    if (apply_hadamard) &
         call hadamard_for_repr_matrix( N_obs, Observations, &
         xcompact, ycompact, Repr_Matrix )
    
    ! form matrix W = HPHt+ R
    
    Repr_matrix = Repr_matrix + Obs_cov

    ! maybe later: save representer matrix into file
    
    ! decompose W once (look into LAPACKs sgelss/dgelss, ask Christian)
    
#ifdef MKL_AVAILABLE
    call getrf(Repr_Matrix, indx, info=lapack_info)

    if (lapack_info .ne. 0) &
         call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'singular matrix after getrf')

#else
    call ludcmp( Repr_Matrix, N_obs, indx ) 
#endif

    ! --------------------------------------
           
    ! update each ensemble member
    
    do n_e=1,N_ens
       
       ! use random measurement error field, get Zpert = Z + v,
       ! compute right hand side rhs = Zpert - H*y^f for system equation,
       
       do i=1,N_obs
         
          rhs(i) = Observations(i)%obs + Obs_err(i,n_e) - Obs_pred(i,n_e)
          
          ! in case of inflation, correct for the fact that "Obs_pred" was not inflated above
          !  because it is intent(in)

          if (inflation_factor > 0.)  rhs(i) = rhs(i) - (1.-inflation_factor)*Obs_pred_prime(i,n_e)
          
       end do
       
       ! solve W*b = Zpert - H*y^f
       ! after back subst, rhs contains (HPHt+R)^(-1)*delta_y(n_e)       
#ifdef MKL_AVAILABLE
       call getrs(Repr_matrix, indx, rhs, trans='N', info=lapack_info)

       if (lapack_info .ne. 0) pause 'something went wrong after getrs'
#else
       call lubksb( Repr_matrix, N_obs, indx, rhs )
#endif
       
       if (.not. apply_hadamard) then
          
          ! update *without* Hadamard product
          
          ! compute w = (Q_prime)t b
          
          weights = matmul( transpose(Obs_pred_prime), rhs)
          
          ! compute new initial cond y^a = y^f + (Y-ybar)*w/(N_e-1)
          
          ! start with (Y-ybar)*w, write into Y_Vector
          
          State_incr_tmp    = matmul( State_prime, weights)
          
          State_incr(:,n_e) = State_incr_tmp
          
       else
          
          ! update *with* Hadamard product
          
          State_incr(:,n_e) = 0.           ! State_incr = analysis - forecast
          
          do ii=1,N_state
             
             do jj=1,N_obs
                
                ! compute [PHt]_ij  (normalize later)
                
                PHt_ij = 0.
                
                do kk=1,N_ens
                   
                   PHt_ij = PHt_ij + State_prime(ii,kk)*Obs_pred_prime(jj,kk)
                   
                end do
                
                ! apply Hadamard factor
                   
                dx=State_lon(ii)-Observations(jj)%lon
                dy=State_lat(ii)-Observations(jj)%lat
                
                ! multiply [PHt]_ij with Hadamard factor
                
                PHt_ij = &
                     PHt_ij * get_gaspari_cohn( dx, dy, xcompact, ycompact )
                   
                State_incr(ii,n_e) = State_incr(ii,n_e) + PHt_ij*rhs(jj)
                
             end do
             
          end do
          
       end if
       
       
       ! finish computation of increment for ensemble member n_e
       ! (normalization is NOT ensemble average, see Eq above)

       State_incr(:,n_e) = State_incr(:,n_e)/real(N_ens-1)

       ! correct for the fact that the increment will be applied to 
       ! the un-inflated state vector outside of this subroutine
       ! (note also that State_prime here has been inflated!)

       if (inflation_factor > 0.) & 
          State_incr(:,n_e) = State_incr(:,n_e) + (1.-1./inflation_factor)*State_prime(:,n_e)
       
    end do
    
  end subroutine enkf_increments
  

  ! *********************************************************************
  
  subroutine hadamard_for_repr_matrix( N_obs, Observations,             & 
       xcompact, ycompact,                                              &
       Repr_Matrix )
    
    implicit none
    
    integer, intent(in) :: N_obs
    
    type(obs_type), intent(in), dimension(N_obs) :: Observations       
    
    real, intent(in) :: xcompact       ! [deg] longitude
    real, intent(in) :: ycompact       ! [deg] latitude
    
    real, dimension(N_obs,N_obs), intent(inout) :: Repr_matrix
    
    ! locals
    
    integer :: i, j
    
    real :: tmpreal, dx, dy
    
    ! ----------------------------------
    
    do i=1,N_obs
       do j=i+1,N_obs
          
          dx = Observations(i)%lon - Observations(j)%lon
          dy = Observations(i)%lat - Observations(j)%lat
          
          tmpreal = get_gaspari_cohn( dx, dy, xcompact, ycompact ) 
          
          Repr_matrix(i,j) = tmpreal * Repr_matrix(i,j)
          Repr_matrix(j,i) = tmpreal * Repr_matrix(j,i)
          
       end do
    end do

  end subroutine hadamard_for_repr_matrix
  
  ! *********************************************************************  

  !DEC$ ATTRIBUTES FORCEINLINE :: get_gaspari_cohn
  !
  function get_gaspari_cohn(dx, dy, xcompact, ycompact) result(rslt)
    
    ! evaluate 5th-order polynomial from Gaspari & Cohn, 1999, Eq (4.10)
    !
    ! get_gaspari_cohn() uses a generalized *an*isotropic Gaspari & Cohn
    ! approach (essentially coordinate stretching, see handwritten
    ! notes for details)
    !
    ! d = separation distance relative to the distance at which all
    ! correlations vanish. In the isotropic case, Gaspari & Cohn, 1999,
    ! Eq. (4.10)
    !
    !    d = sqrt(dx**2 + dy**2) / (2*c) = |z| / (2*c)    
    !
    ! or in the anisotropic case
    !
    !    d = sqrt( (dx/xcompact)**2 + (dy/ycompact)**2 )
    !
    ! *** Use |z|/c = 2*d. All correlations vanish for d > 1. ***
    !
    ! for a given lat/lon distance (dx and dy, resp.), compute the
    !  anisotropic compact support (Gaspari & Cohn) weights
    !
    ! input distances must be in degrees latitude/longitude
    !  
    !  dx = longitude separation of two points       [deg]
    !  dy = latitude separation of two points        [deg]
    !
    !  xcompact = longitude scale of compact support [deg]
    !  ycompact = latitude scale of compact support  [deg]
    !
    ! All correlations vanish outside of an ellipse with semi-axes 
    ! xcompact and ycompact, ie Gaspari & Cohn weights vanish 
    ! for d > 1 (note the factor 2!)
    !
    ! When the anisotropic case is reduced back to the isotropic case,
    ! (ie if xcompact==ycompact) then c = xcompact/2 = ycompact/2.
    !
    ! pchakrab, rreichle: revised, 17 Sep 2013
    ! 
    ! ------------------------------------------------------------------
            
    implicit none
    
    real, intent(in) :: dx, dy, xcompact, ycompact

    real             :: rslt ! returned value

    ! local variables
    
    real             :: d, dsq

    real, parameter  :: tol = 1e-3
    
    ! ---------------------------------------------------------
    
    if ( (abs(dx)>xcompact) .or. (abs(dy)>ycompact) ) then
       
       rslt = 0.   ! nothing to do
       
    else
       
       ! compute (anisotropic) distance relative to compact support
       !
       !      d = sqrt( (dx/xcompact)**2 + (dy/ycompact)**2 )
       !
       ! NOTE: multiply d by 2 to return to Gaspari & Cohn, 1999, notation
       
       dsq = 4.0*((dx*dx)/(xcompact*xcompact) + (dy*dy)/(ycompact*ycompact))
       
       d = sqrt(dsq)

       if (d >= 2.) then

          rslt = 0.

       else if (d <= tol) then              

          rslt = 1.

       else if (d <= 1.) then

          ! y = -.25*d**5 + .5*d**4 + .625*d**3 - 5./3.*d**2 + 1.
          
          ! rslt = d*d *( d*( d*( -.25*d + .5) + .625) -5./3.) + 1.

          rslt = dsq*(dsq*(-0.25*d + 0.5) + 0.625*d - 5.0/3.0) + 1.0
          
       else
          
          ! y = d**5/12. - .5*d**4 + .625*d**3 + 5./3.*d**2 - 5.*d + 4. - 2./3./d

          rslt = d*( d*( d*( d*( d/12. - .5) + .625) + 5./3.) -5.) + 4. - 2./3./d

       end if
       
    end if
    
  end function get_gaspari_cohn
  
  ! ************************************************************

end module enkf_general

! *******************************************************************

#if 0

program test_gaspari_cohn

  ! reichle, 17 Sep 2013: updated (only works if module stuff is removed)

  implicit none
  
  real :: x, y, get_gaspari_cohn

  real :: xcompact = 5.
  real :: ycompact = 4.

  integer :: i,j,N
  
  N = 40
    
  do i=-N,N
     
     x = 7.*float(i)/float(N)
     
     do j=-N,N
        
        y = 7.*float(j)/float(N)
        
        write (999,*) x, y, get_gaspari_cohn(x,y,xcompact,ycompact)
        
     end do
  end do
      
end program test_gaspari_cohn

#endif   
  

! ********** EOF *********************************************************
