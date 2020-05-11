

module clsm_ensdrv_drv_routines
  
  ! collection of subroutines for enkf_driver written in f90
  ! reichle, 10 May 2005
  !
  ! reichle, 13 Aug 2008 - moved forcing subroutines into 
  !                         clsm_ensdrv_force_routines.F90
  ! reichle, 28 Oct 2008 - added soilcls30 and soilcls100
  !                      - optimized restart-to-exp-domain mapping in initialize_model()
  ! reichle,  5 Apr 2013 - revised treatment of output collections

  use catch_constants,                  ONLY:     &
       N_snow => CATCH_N_SNOW,                    &
       N_gt   => CATCH_N_GT

  use catch_incr,                       ONLY:     &
       check_catch_progn
  
  use catch_types,                      ONLY:     &
       cat_param_type,                            &
       cat_progn_type

  use LDAS_ensdrv_mpi,                  ONLY:     &
       mpicomm,                                   &
       mpierr,                                    &
       numprocs,                                  &
       master_proc

  implicit none

  include 'mpif.h'

  private
  
  public :: check_cat_progn
  public :: l2f_real
  public :: f2l_real
  public :: f2l_real8
  public :: f2l_logical

  character(10), private :: tmpstring10
  character(40), private :: tmpstring40

contains
  
  ! ********************************************************************
    
  subroutine check_cat_progn( N_cat, cat_param, cat_progn )
    
    ! wrapper for subroutine check_catch_progn() which has been 
    ! moved to "catch_iau.F90" in GEOScatch_GridComp - reichle, 3 Apr 2012
    
    implicit none
    
    integer, intent(in) :: N_cat
    
    type(cat_param_type), dimension(N_cat), intent(in)    :: cat_param
    
    type(cat_progn_type), dimension(N_cat), intent(inout) :: cat_progn
    
    ! ----------------------------------------------------------------
    
    ! local variables
    
    integer :: k 
    
    real, dimension(N_snow,N_cat) :: wesn
    real, dimension(N_snow,N_cat) :: htsn
    real, dimension(N_snow,N_cat) :: sndz

    real, dimension(N_gt,  N_cat) :: ghtcnt
    
    ! ----------------------------------------------------------------
    
    ! copy select cat_progn fields into 2-d arrays
    
    do k=1,N_gt
       
       GHTCNT(k,:) = cat_progn%ght(k)
       
    end do
    
    do k=1,N_snow
       
       WESN(k,:)   = cat_progn%wesn(k)
       HTSN(k,:)   = cat_progn%htsn(k)
       SNDZ(k,:)   = cat_progn%sndz(k)
       
    end do

    ! check for consistency and unphysical values
    
    call check_catch_progn( N_cat, cat_param%vegcls, cat_param%dzsf,         &
         cat_param%vgwmax,  cat_param%cdcr1, cat_param%cdcr2,                &
         cat_param%psis, cat_param%bee, cat_param%poros, cat_param%wpwet,    &
         cat_param%ars1, cat_param%ars2, cat_param%ars3,                     &
         cat_param%ara1, cat_param%ara2, cat_param%ara3, cat_param%ara4,     &
         cat_param%arw1, cat_param%arw2, cat_param%arw3, cat_param%arw4,     &
         cat_progn%tc1, cat_progn%tc2, cat_progn%tc4,                        &
         cat_progn%qa1, cat_progn%qa2, cat_progn%qa4,                        &
         cat_progn%capac, cat_progn%catdef,                                  &
         cat_progn%rzexc, cat_progn%srfexc,                                  &
         ghtcnt, wesn, htsn, sndz            )
    
    ! copy 2-d arrays back into cat_progn fields
    
    do k=1,N_gt
       
       cat_progn%ght(k)  = ghtcnt(k,:)
       
    end do
    
    do k=1,N_snow
       
       cat_progn%wesn(k) = WESN(k,:)
       cat_progn%htsn(k) = HTSN(k,:)
       cat_progn%sndz(k) = SNDZ(k,:)
       
    end do

  end subroutine check_cat_progn

  ! *************************************************************
    
  subroutine l2f_real( N_l, N_f, N_l_vec, low_ind, data_l, data_f)
    
    ! wrapper for MPI_GATHERV applied to MPI_REAL vector

    ! reichle, 23 Dec 2011
    
    implicit none
  
    integer,                      intent(in)  :: N_l, N_f
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    real,    dimension(N_l),      intent(in)  :: data_l
    real,    dimension(N_f),      intent(out) :: data_f
    
#ifdef LDAS_MPI
    call MPI_GATHERV( &
         data_l, N_l,                MPI_REAL,          &
         data_f, N_l_vec, low_ind-1, MPI_REAL,          &
         0, mpicomm, mpierr )

    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_f = data_l
#endif              
    
  end subroutine l2f_real

  ! *************************************************************
  
  subroutine l2f_real8( N_l, N_f, N_l_vec, low_ind, data_l, data_f)
    
    ! wrapper for MPI_GATHERV applied to MPI_REAL8 vector
    
    ! reichle, 31 Jan 2014
    
    implicit none
    
    integer,                      intent(in)  :: N_l, N_f
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    real*8,  dimension(N_l),      intent(in)  :: data_l
    real*8,  dimension(N_f),      intent(out) :: data_f
    
#ifdef LDAS_MPI
    call MPI_GATHERV( &
         data_l, N_l,                MPI_REAL8,         &
         data_f, N_l_vec, low_ind-1, MPI_REAL8,         &
         0, mpicomm, mpierr )
    
    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_f = data_l
#endif              
    
  end subroutine l2f_real8
  
  ! *************************************************************
  
  subroutine f2l_real( N_f, N_l, N_l_vec, low_ind, data_f, data_l)
    
    ! wrapper for MPI_SCATTERV applied to MPI_REAL vector

    ! reichle, 25 Jul 2013
    
    implicit none
  
    integer,                      intent(in)  :: N_f, N_l
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    real,    dimension(N_f),      intent(in)  :: data_f
    real,    dimension(N_l),      intent(out) :: data_l
    
#ifdef LDAS_MPI
    call MPI_SCATTERV( &
         data_f, N_l_vec, low_ind-1, MPI_REAL,          &
         data_l, N_l,                MPI_REAL,          &
         0, mpicomm, mpierr )

    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_l = data_f
#endif              
    
  end subroutine f2l_real
  
  ! *************************************************************
  
  subroutine f2l_real8( N_f, N_l, N_l_vec, low_ind, data_f, data_l)
    
    ! wrapper for MPI_SCATTERV applied to MPI_REAL8 vector

    ! reichle, 31 Jan 2014
    
    implicit none
  
    integer,                      intent(in)  :: N_f, N_l
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    real*8,  dimension(N_f),      intent(in)  :: data_f
    real*8,  dimension(N_l),      intent(out) :: data_l
    
#ifdef LDAS_MPI
    call MPI_SCATTERV( &
         data_f, N_l_vec, low_ind-1, MPI_REAL8,          &
         data_l, N_l,                MPI_REAL8,          &
         0, mpicomm, mpierr )

    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_l = data_f
#endif              
    
  end subroutine f2l_real8
  
  ! *************************************************************
  
  subroutine f2l_logical( N_f, N_l, N_l_vec, low_ind, data_f, data_l)
    
    ! wrapper for MPI_SCATTERV applied to MPI_LOGICAL vector

    ! reichle,  6 Jun 2016
    
    implicit none
  
    integer,                      intent(in)  :: N_f, N_l
    integer, dimension(numprocs), intent(in)  :: N_l_vec, low_ind  
    logical, dimension(N_f),      intent(in)  :: data_f
    logical, dimension(N_l),      intent(out) :: data_l
    
#ifdef LDAS_MPI
    call MPI_SCATTERV( &
         data_f, N_l_vec, low_ind-1, MPI_LOGICAL,        &
         data_l, N_l,                MPI_LOGICAL,        &
         0, mpicomm, mpierr )

    call MPI_BARRIER( mpicomm, mpierr )
#else                 
    data_l = data_f
#endif              
    
  end subroutine f2l_logical
  
  ! *************************************************************

end module CLSM_ensdrv_drv_routines

! *********** EOF **************************************************
