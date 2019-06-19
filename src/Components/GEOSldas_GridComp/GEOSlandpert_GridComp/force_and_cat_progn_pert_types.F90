
module force_and_cat_progn_pert_types
  
  ! types and subroutines in this module are primarily meant to 
  ! facilitate assembly of force_pert_param and progn_pert_param
  ! structures that are the parameter inputs to the generation of
  ! perturbations
  
  ! reichle,  1 Jun 2005
  ! reichle,  6 Dec 2013 - introduced "progn_pert_type" 
  !                        (no longer use "cat_progn_type" for progn perts)
  ! reichle, 21 Nov 2014 - renamed force_pert_type fields for consistency w/ met_force_type
  !                          %tmp2m --> %tair  (but note lower-case!)
  !                          %dpt2m --> %qair  (but note lower-case!)
  !                          %wnd   --> %wind  (but note lower-case!)

  use catch_constants,                  ONLY:     &
       N_gt => CATCH_N_GT

  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: N_force_pert_max
  
  public :: force_pert_real_type
  public :: force_pert_logi_type
  public :: force_pert_char_type
  
  public :: force_pert_ccor_type  

  public :: N_progn_pert_max
  
  public :: progn_pert_real_type
  public :: progn_pert_logi_type
  public :: progn_pert_char_type
  
  public :: progn_pert_ccor_type  
  
  public :: struct2vec_force_pert, struct2mat_force_pert_ccor
  public :: struct2vec_progn_pert, struct2mat_progn_pert_ccor
  
  public :: assignment (=)
  
  ! ----------------------------------------------------------------------
  !
  ! force_pert_type and progn_pert_type are used to gather input data and 
  !  assemble force_pert_param and progn_pert_param
  !
  ! ----------------------------------------------------------------------


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer :: N_force_pert_max = 6   ! must equal the number of fields below
  
  type :: force_pert_real_type
     
     real :: pcp
     real :: sw
     real :: lw
     real :: tair
     real :: qair
     real :: wind   
     
  end type force_pert_real_type
  
  ! -----------------------------------------

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  type :: force_pert_logi_type
     
     logical :: pcp
     logical :: sw
     logical :: lw
     logical :: tair
     logical :: qair
     logical :: wind   
     
  end type force_pert_logi_type

  ! -----------------------------------------

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  type :: force_pert_char_type
     
     character(40) :: pcp
     character(40) :: sw
     character(40) :: lw
     character(40) :: tair
     character(40) :: qair
     character(40) :: wind   
     
  end type force_pert_char_type
  
  ! -----------------------------------------
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
  type :: force_pert_ccor_type
     
     type(force_pert_real_type) :: pcp
     type(force_pert_real_type) :: sw   
     type(force_pert_real_type) :: lw   
     type(force_pert_real_type) :: tair
     type(force_pert_real_type) :: qair
     type(force_pert_real_type) :: wind  

  end type force_pert_ccor_type
  

  ! --------------------------------------------------------------------------------

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer :: N_progn_pert_max = 5+N_gt   ! must equal the number of fields below
  
  type :: progn_pert_real_type
     
     real                  :: catdef
     real                  :: rzexc
     real                  :: srfexc
     real                  :: snow  
     real                  :: tc
     real, dimension(N_gt) :: ght      
     
  end type progn_pert_real_type
  
  ! -----------------------------------------

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  type :: progn_pert_logi_type
     
     logical                  :: catdef 
     logical                  :: rzexc  
     logical                  :: srfexc 
     logical                  :: snow   
     logical                  :: tc     
     logical, dimension(N_gt) :: ght    
     
  end type progn_pert_logi_type

  ! -----------------------------------------

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  type :: progn_pert_char_type
     
     character(40)                  :: catdef 
     character(40)                  :: rzexc  
     character(40)                  :: srfexc 
     character(40)                  :: snow   
     character(40)                  :: tc     
     character(40), dimension(N_gt) :: ght    
     
  end type progn_pert_char_type
  
  ! -----------------------------------------
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WARNING: When modifying this derived type make sure that the corresponding
  !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
  !          any subroutines or operators defined herein
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
  type :: progn_pert_ccor_type
     
     type(progn_pert_real_type)                  :: catdef 
     type(progn_pert_real_type)                  :: rzexc  
     type(progn_pert_real_type)                  :: srfexc 
     type(progn_pert_real_type)                  :: snow   
     type(progn_pert_real_type)                  :: tc     
     type(progn_pert_real_type), dimension(N_gt) :: ght    

  end type progn_pert_ccor_type
    
  ! ----------------------------------------------------------------------
  !
  ! interface definitions

  interface struct2vec_force_pert
     module procedure struct2vec_force_pert_real
     module procedure struct2vec_force_pert_logi
     module procedure struct2vec_force_pert_char
  end interface
  
  ! -----------------------------------------
  
  interface struct2vec_progn_pert
     module procedure struct2vec_progn_pert_real
     module procedure struct2vec_progn_pert_logi
     module procedure struct2vec_progn_pert_char
  end interface

  ! -----------------------------------------

  interface assignment (=)
     module procedure scalar2force_pert_real
     module procedure scalar2force_pert_ccor
     
     module procedure scalar2progn_pert_real
     module procedure scalar2progn_pert_ccor
  end interface

  ! *********************************************************************

contains
  
  subroutine struct2vec_force_pert_real(force_pert_struct,force_pert_vec)
    
    implicit none
    
    type(force_pert_real_type), intent(in) :: force_pert_struct
    
    real, dimension(N_force_pert_max), intent(out) :: force_pert_vec
    
    force_pert_vec(1) = force_pert_struct%pcp
    force_pert_vec(2) = force_pert_struct%sw
    force_pert_vec(3) = force_pert_struct%lw
    force_pert_vec(4) = force_pert_struct%tair
    force_pert_vec(5) = force_pert_struct%qair
    force_pert_vec(6) = force_pert_struct%wind
    
  end subroutine struct2vec_force_pert_real
  
  ! -----------------------------------------
  
  subroutine struct2vec_force_pert_logi(force_pert_struct,force_pert_vec)
    
    implicit none
    
    type(force_pert_logi_type),           intent(in)  :: force_pert_struct
    
    logical, dimension(N_force_pert_max), intent(out) :: force_pert_vec
    
    force_pert_vec(1) = force_pert_struct%pcp
    force_pert_vec(2) = force_pert_struct%sw
    force_pert_vec(3) = force_pert_struct%lw
    force_pert_vec(4) = force_pert_struct%tair
    force_pert_vec(5) = force_pert_struct%qair
    force_pert_vec(6) = force_pert_struct%wind
    
  end subroutine struct2vec_force_pert_logi
  
  ! -----------------------------------------

  subroutine struct2vec_force_pert_char(force_pert_struct,force_pert_vec)
    
    implicit none
    
    type(force_pert_char_type),                 intent(in)  :: force_pert_struct
    
    character(40), dimension(N_force_pert_max), intent(out) :: force_pert_vec
    
    force_pert_vec(1) = force_pert_struct%pcp
    force_pert_vec(2) = force_pert_struct%sw
    force_pert_vec(3) = force_pert_struct%lw
    force_pert_vec(4) = force_pert_struct%tair
    force_pert_vec(5) = force_pert_struct%qair
    force_pert_vec(6) = force_pert_struct%wind  
    
  end subroutine struct2vec_force_pert_char
  

  ! -------------------------------------------------------------------------
  
  subroutine struct2mat_force_pert_ccor(force_pert_ccor_struct, &
       force_pert_ccor_mat )
    
    implicit none
    
    type(force_pert_ccor_type) :: force_pert_ccor_struct
    
    real, dimension(N_force_pert_max,N_force_pert_max), intent(out) :: &
         force_pert_ccor_mat
        
    call struct2vec_force_pert( force_pert_ccor_struct%pcp , force_pert_ccor_mat(1,:) )
    call struct2vec_force_pert( force_pert_ccor_struct%sw  , force_pert_ccor_mat(2,:) )
    call struct2vec_force_pert( force_pert_ccor_struct%lw  , force_pert_ccor_mat(3,:) )
    call struct2vec_force_pert( force_pert_ccor_struct%tair, force_pert_ccor_mat(4,:) )
    call struct2vec_force_pert( force_pert_ccor_struct%qair, force_pert_ccor_mat(5,:) )
    call struct2vec_force_pert( force_pert_ccor_struct%wind, force_pert_ccor_mat(6,:) )
    
  end subroutine struct2mat_force_pert_ccor
  
  ! -------------------------------------------------------------------------
  
  subroutine scalar2force_pert_real( force_pert_real, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(force_pert_real_type), intent(out) :: force_pert_real
    
    force_pert_real%pcp      = scalar
    force_pert_real%sw       = scalar
    force_pert_real%lw       = scalar
    force_pert_real%tair     = scalar
    force_pert_real%qair     = scalar
    force_pert_real%wind     = scalar
        
  end subroutine scalar2force_pert_real
  
  ! ------------------------------------------
  
  subroutine scalar2force_pert_ccor( force_pert_ccor, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(force_pert_ccor_type), intent(out) :: force_pert_ccor
    
    force_pert_ccor%pcp      = scalar
    force_pert_ccor%sw       = scalar
    force_pert_ccor%lw       = scalar
    force_pert_ccor%tair     = scalar
    force_pert_ccor%qair     = scalar
    force_pert_ccor%wind     = scalar
        
  end subroutine scalar2force_pert_ccor
  
  
  ! *********************************************************************

  subroutine struct2vec_progn_pert_real( progn_pert_struct, progn_pert_vec )
    
    implicit none
    
    type(progn_pert_real_type), intent(in) :: progn_pert_struct
    
    real, dimension(N_progn_pert_max), intent(out) :: progn_pert_vec

    integer :: i
    
    progn_pert_vec( 1) = progn_pert_struct%catdef 
    progn_pert_vec( 2) = progn_pert_struct%rzexc  
    progn_pert_vec( 3) = progn_pert_struct%srfexc 
    progn_pert_vec( 4) = progn_pert_struct%snow   
    progn_pert_vec( 5) = progn_pert_struct%tc     
    
    do i=1,N_gt
       progn_pert_vec(5+i) = progn_pert_struct%ght(i)
    end do
    
  end subroutine struct2vec_progn_pert_real
  
  ! -----------------------------------------

  subroutine struct2vec_progn_pert_logi( progn_pert_struct, progn_pert_vec )
    
    implicit none
    
    type(progn_pert_logi_type), intent(in) :: progn_pert_struct
    
    logical, dimension(N_progn_pert_max), intent(out) :: progn_pert_vec

    integer :: i
    
    progn_pert_vec( 1) = progn_pert_struct%catdef 
    progn_pert_vec( 2) = progn_pert_struct%rzexc  
    progn_pert_vec( 3) = progn_pert_struct%srfexc 
    progn_pert_vec( 4) = progn_pert_struct%snow   
    progn_pert_vec( 5) = progn_pert_struct%tc     
    
    do i=1,N_gt
       progn_pert_vec(5+i) = progn_pert_struct%ght(i)
    end do
    
  end subroutine struct2vec_progn_pert_logi
  
  ! -----------------------------------------

  subroutine struct2vec_progn_pert_char( progn_pert_struct, progn_pert_vec )
    
    implicit none
    
    type(progn_pert_char_type), intent(in) :: progn_pert_struct
    
    character(40), dimension(N_progn_pert_max), intent(out) :: progn_pert_vec

    integer :: i
    
    progn_pert_vec( 1) = progn_pert_struct%catdef 
    progn_pert_vec( 2) = progn_pert_struct%rzexc  
    progn_pert_vec( 3) = progn_pert_struct%srfexc 
    progn_pert_vec( 4) = progn_pert_struct%snow   
    progn_pert_vec( 5) = progn_pert_struct%tc     
    
    do i=1,N_gt
       progn_pert_vec(5+i) = progn_pert_struct%ght(i)
    end do
    
  end subroutine struct2vec_progn_pert_char
  
  ! -----------------------------------------

  subroutine struct2mat_progn_pert_ccor( progn_pert_ccor_struct, &
       progn_pert_ccor_mat )
    
    implicit none
    
    type(progn_pert_ccor_type), intent(in) :: progn_pert_ccor_struct
    
    real, dimension(N_progn_pert_max,N_progn_pert_max), intent(out) :: &
         progn_pert_ccor_mat

    integer :: i

    call struct2vec_progn_pert( progn_pert_ccor_struct%catdef , progn_pert_ccor_mat(  1,:) )
    call struct2vec_progn_pert( progn_pert_ccor_struct%rzexc  , progn_pert_ccor_mat(  2,:) )
    call struct2vec_progn_pert( progn_pert_ccor_struct%srfexc , progn_pert_ccor_mat(  3,:) )
    call struct2vec_progn_pert( progn_pert_ccor_struct%snow   , progn_pert_ccor_mat(  4,:) )
    call struct2vec_progn_pert( progn_pert_ccor_struct%tc     , progn_pert_ccor_mat(  5,:) )
    
    do i=1,N_gt
       call struct2vec_progn_pert(progn_pert_ccor_struct%ght(i),progn_pert_ccor_mat(5+i,:) )
    end do

  end subroutine struct2mat_progn_pert_ccor

  ! -----------------------------------------
     
  subroutine scalar2progn_pert_real( progn_pert_real, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(progn_pert_real_type), intent(out) :: progn_pert_real

    integer :: i

    progn_pert_real%catdef       = scalar
    progn_pert_real%rzexc        = scalar
    progn_pert_real%srfexc       = scalar
    progn_pert_real%snow         = scalar
    progn_pert_real%tc           = scalar

    do i=1,N_gt
       progn_pert_real%ght(i)    = scalar
    end do
    
  end subroutine scalar2progn_pert_real
      
  ! -----------------------------------------

  subroutine scalar2progn_pert_ccor( progn_pert_ccor, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(progn_pert_ccor_type), intent(out) :: progn_pert_ccor

    integer :: i

    progn_pert_ccor%catdef       = scalar
    progn_pert_ccor%rzexc        = scalar
    progn_pert_ccor%srfexc       = scalar
    progn_pert_ccor%snow         = scalar
    progn_pert_ccor%tc           = scalar

    do i=1,N_gt
       progn_pert_ccor%ght(i)    = scalar
    end do
    
  end subroutine scalar2progn_pert_ccor
      
end module force_and_cat_progn_pert_types

! *******************************************************************

#if 0

program test

  use force_and_cat_progn_pert_types

  implicit none
  
  type(force_pert_real_type) :: force_pert_struct_real
  
  real, dimension(:), allocatable :: force_pert_vec_real

  type(force_pert_logical_type) :: force_pert_struct_logical
  
  logical, dimension(:), allocatable :: force_pert_vec_logical
  
  type(force_pert_character_type) :: force_pert_struct_character
  
  character(40), dimension(:), allocatable :: force_pert_vec_character



  type(force_pert_ccorr_type) :: force_pert_ccorr
      
  real, dimension(:,:), allocatable :: force_pert_ccorr_mat
  

  
  !--------------------------------------------------------
  

  allocate(force_pert_vec_real(N_force_pert_max))
  

  force_pert_struct_real%pcp     =1. 
  force_pert_struct_real%sw      =2.
  force_pert_struct_real%lw      =3.
  force_pert_struct_real%tair    =4.
  force_pert_struct_real%qair    =5.
  force_pert_struct_real%wind    =6.

  call struct2vec_force_pert( force_pert_struct_real, force_pert_vec_real)
  
  write (*,*) force_pert_vec_real


  allocate(force_pert_vec_logical(N_force_pert_max))
  
  
  force_pert_struct_logical%pcp     =.true. 
  force_pert_struct_logical%sw      =.true.
  force_pert_struct_logical%lw      =.true.
  force_pert_struct_logical%tair    =.false.
  force_pert_struct_logical%qair    =.true.
  force_pert_struct_logical%wind    =.true.
  
  call struct2vec_force_pert( force_pert_struct_logical, force_pert_vec_logical)
  
  write (*,*) force_pert_vec_logical
  
  allocate(force_pert_vec_character(N_force_pert_max))
  
  
  force_pert_struct_character%pcp     ='1.asdf'
  force_pert_struct_character%sw      ='2asdf.'
  force_pert_struct_character%lw      ='3.asdf'
  force_pert_struct_character%tair    ='4.asdf'
  force_pert_struct_character%qair    ='5.asdf'
  force_pert_struct_character%wind    ='6.asdf'

  call struct2vec_force_pert( force_pert_struct_character, force_pert_vec_character)
  
  write (*,*) force_pert_vec_character

  ! -----------------------------------
  
  force_pert_ccorr%pcp%pcp  =1.
  force_pert_ccorr%pcp%sw   =2.
  force_pert_ccorr%pcp%lw   =3.
  force_pert_ccorr%pcp%tair =4.
  force_pert_ccorr%pcp%qair =5.
  force_pert_ccorr%pcp%wind =6.

  force_pert_ccorr%sw%pcp  =10.
  force_pert_ccorr%sw%sw   =20.
  force_pert_ccorr%sw%lw   =30.
  force_pert_ccorr%sw%tair =40.
  force_pert_ccorr%sw%qair =50.
  force_pert_ccorr%sw%wind =60.

  force_pert_ccorr%lw%pcp  =100.
  force_pert_ccorr%lw%sw   =200.
  force_pert_ccorr%lw%lw   =300.
  force_pert_ccorr%lw%tair =400.
  force_pert_ccorr%lw%qair =500.
  force_pert_ccorr%lw%wind =600.

  force_pert_ccorr%tair%pcp  =1000.
  force_pert_ccorr%tair%sw   =2000.
  force_pert_ccorr%tair%lw   =3000.
  force_pert_ccorr%tair%tair =4000.
  force_pert_ccorr%tair%qair =5000.
  force_pert_ccorr%tair%wind =6000.

  force_pert_ccorr%qair%pcp  =10000.
  force_pert_ccorr%qair%sw   =20000.
  force_pert_ccorr%qair%lw   =30000.
  force_pert_ccorr%qair%tair =40000.
  force_pert_ccorr%qair%qair =50000.
  force_pert_ccorr%qair%wind =60000.

  force_pert_ccorr%wind%pcp  =100000.
  force_pert_ccorr%wind%sw   =200000.
  force_pert_ccorr%wind%lw   =300000.
  force_pert_ccorr%wind%tair =400000.
  force_pert_ccorr%wind%qair =500000.
  force_pert_ccorr%wind%wind =600000.

  allocate(force_pert_ccorr_mat(N_force_pert_max,N_force_pert_max))
  
  !! force_pert_ccorr = 1.11111

  call struct2mat_force_pert_ccorr(force_pert_ccorr, &
       force_pert_ccorr_mat )


  write (*,*) 'force_pert_ccorr=', force_pert_ccorr
  
  write (*,*) 'force_pert_ccorr_mat=', force_pert_ccorr_mat

  
end program test
   

#endif

! =============== EOF ================================================
