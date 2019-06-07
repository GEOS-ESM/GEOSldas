
module mwRTM_types

  ! definition of types and associated operators for microwave radiative transfer model
  !
  ! IMPORTANT:
  ! When adding a field to any of the derived types, must also update
  ! the associated assignment and operator definitions.
  ! THERE IS NO WARNING/ERROR IF OPERATOR IS NOT DEFINED FOR ALL FIELDS!
  !
  ! reichle, 16 May 2011
  ! reichle, 21 Oct 2011 - added field "poros" to "mwRTM_param_type"
  !
  ! --------------------------------------------------------------------------
  
  use LDAS_ensdrv_globals,           ONLY:     &
       nodata_generic,                            &
       nodata_tol_generic

  use ldas_exceptionsMod,                  ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: mwRTM_param_type
  public :: io_mwRTM_param_type
  public :: mwRTM_param_nodata_check
  
  public :: assignment (=)
  
  ! ---------------------------------------------------------
  
  ! model parameters
  
  type :: mwRTM_param_type
     
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     integer :: vegcls      ! land cover/veg class (can differ from cat_param%vegcls   ! )
     integer :: soilcls     ! soil class           (can differ from cat_param%soilcls30! )
     
     real    :: sand        ! sand fraction                                 [0-1]
     real    :: clay        ! clay fraction                                 [0-1]
     
     real    :: poros       ! porosity (can differ from cat_param%poros!)   [m3/m3]

     ! parameters relating to soil mixing model
     
     real    :: wang_wt     ! Wang dielectric model: transition s.m.        [m3/m3]
     real    :: wang_wp     ! Wang dielectric model: wilting point s.m.     [m3/m3]
          
     ! parameters relating to emissivity of rough surface
     
     real    :: rgh_hmin    ! min roughness                                 [dim-less]
     real    :: rgh_hmax    ! max roughness                                 [dim-less]
     real    :: rgh_wmin    ! soil moisture for transition to hmax          [m3/m3]
     real    :: rgh_wmax    ! soil moisture for transition to hmin          [m3/m3]
     real    :: rgh_Nrh     ! h-pol exponent for inc angle parameterization [dim-less]
     real    :: rgh_Nrv     ! v-pol exponent for inc angle parameterization [dim-less]
     real    :: rgh_polmix  ! polarization mixing parameter                 [dim-less]
     
     real    :: omega       ! single scattering albedo                      [dim-less]

     ! parameters relating to vegetation opacity (tau)
     
     real    :: bh          ! veg b parameter (h-pol)  (tau = b*VWC)        [dim-less]
     real    :: bv          ! veg b parameter (v-pol)  (tau = b*VWC)        [dim-less]
     real    :: lewt        ! VWC = lewt*LAI                                [kg/m2]

     
  end type mwRTM_param_type
  
  ! ---------------------------------------------------------
  
  interface assignment (=)
     module procedure scalar2mwRTM_param
  end interface
  
contains

  ! **********************************************************************
  
  subroutine io_mwRTM_param_type( action, unitnum, N_tile, mwp, N_catg, tile_id, d2g ) 
    
    ! read/write mwRTM_param for domain from/to file
    !
    ! write: write mwRTM params for N_tile tiles in domain
    !
    ! read:  read mwRTM params for 
    
    !
    ! reichle, 1  Jun 2011
    ! reichle, 21 Oct 2011 -- added "read" for mwRTM params
    ! reichle, 22 Oct 2012 -- added check for tile_id to mwRTM param input
    ! 
    ! -------------------------------------------------------------------
    
    implicit none
    
    character,                                 intent(in)    :: action
       
    integer,                                   intent(in)    :: unitnum
    
    integer,                                   intent(in)    :: N_tile  ! =N_catd
    
    type(mwRTM_param_type), dimension(N_tile), intent(inout) :: mwp

    integer, optional,                         intent(in)    :: N_catg    

    integer, optional,      dimension(N_tile), intent(in)    :: tile_id, d2g

    ! local variables
    
    integer :: n, N_tmp

    integer, dimension(:), allocatable :: tmpint
    real,    dimension(:), allocatable :: tmpreal

    character(len=*), parameter :: Iam = 'io_mwRTM_param_type'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------------

    select case (action)
       
    case ('w','W')     ! write mwp for all tiles in domain
       
       write (unitnum) N_tile
       
       write (unitnum) (mwp(n)%vegcls    ,    n=1,N_tile)  ! integer
       write (unitnum) (mwp(n)%soilcls   ,    n=1,N_tile)  ! integer
       
       write (unitnum) (mwp(n)%sand      ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%clay      ,    n=1,N_tile)  ! real 

       write (unitnum) (mwp(n)%poros     ,    n=1,N_tile)  ! real 
       
       write (unitnum) (mwp(n)%wang_wt   ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%wang_wp   ,    n=1,N_tile)  ! real 
       
       write (unitnum) (mwp(n)%rgh_hmin  ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%rgh_hmax  ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%rgh_wmin  ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%rgh_wmax  ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%rgh_Nrh   ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%rgh_Nrv   ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%rgh_polmix,    n=1,N_tile)  ! real 
       
       write (unitnum) (mwp(n)%omega     ,    n=1,N_tile)  ! real 
       
       write (unitnum) (mwp(n)%bh        ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%bv        ,    n=1,N_tile)  ! real 
       write (unitnum) (mwp(n)%lewt      ,    n=1,N_tile)  ! real 
       
    case ('r','R')

       ! read the parameters for all global tiles (similar to read_land_parameters())
       
       ! optional inputs N_catg and d2g must be present
       
       if ( (.not. present(N_catg )) .or.      &
            (.not. present(tile_id)) .or.      &
            (.not. present(d2g    ))         ) then
          err_msg = 'missing optional inputs N_catg, tile_id, d2g'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       ! read how many tiles are in file and double-check against N_catg
       
       read (unitnum) N_tmp
       
       if (N_tmp .ne. N_catg) then
          err_msg = 'number of tiles in file .ne. N_catg'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       ! allocate tmp vectors
       
       allocate(tmpint( N_catg))
       allocate(tmpreal(N_catg))
       
       ! read tile IDs (first record)
       
       read (unitnum) tmpint;

       ! make sure tile IDs match (works only for "SiB2_V2" and newer versions)
       
       if (any(tile_id/=tmpint(d2g(1:N_tile)))) then
          err_msg = 'mismatch of tile IDs for mwRTM_parameters'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       ! read and subset from global tile space to domain
       
       read (unitnum) tmpint;   mwp(1:N_tile)%vegcls     =  tmpint( d2g(1:N_tile))
       read (unitnum) tmpint;   mwp(1:N_tile)%soilcls    =  tmpint( d2g(1:N_tile))
       
       read (unitnum) tmpreal;  mwp(1:N_tile)%sand       =  tmpreal(d2g(1:N_tile))
       read (unitnum) tmpreal;  mwp(1:N_tile)%clay       =  tmpreal(d2g(1:N_tile))
       
       read (unitnum) tmpreal;  mwp(1:N_tile)%poros      =  tmpreal(d2g(1:N_tile))
       
       read (unitnum) tmpreal;  mwp(1:N_tile)%wang_wt    =  tmpreal(d2g(1:N_tile)) 
       read (unitnum) tmpreal;  mwp(1:N_tile)%wang_wp    =  tmpreal(d2g(1:N_tile))
       
       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_hmin   =  tmpreal(d2g(1:N_tile)) 
       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_hmax   =  tmpreal(d2g(1:N_tile))
       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_wmin   =  tmpreal(d2g(1:N_tile))
       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_wmax   =  tmpreal(d2g(1:N_tile))
       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_Nrh    =  tmpreal(d2g(1:N_tile))
       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_Nrv    =  tmpreal(d2g(1:N_tile))
       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_polmix =  tmpreal(d2g(1:N_tile))    
       
       read (unitnum) tmpreal;  mwp(1:N_tile)%omega      =  tmpreal(d2g(1:N_tile))
       
       read (unitnum) tmpreal;  mwp(1:N_tile)%bh         =  tmpreal(d2g(1:N_tile))
       read (unitnum) tmpreal;  mwp(1:N_tile)%bv         =  tmpreal(d2g(1:N_tile))
       read (unitnum) tmpreal;  mwp(1:N_tile)%lewt       =  tmpreal(d2g(1:N_tile))
       
       ! clean up

       deallocate(tmpreal)
       deallocate(tmpint)
       
    case default
       
       err_msg = 'unknown action ' // action
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
  end subroutine io_mwRTM_param_type
  
  ! ************************************************************
  
  subroutine scalar2mwRTM_param( mwRTM_param, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(mwRTM_param_type), intent(out) :: mwRTM_param
    
    ! ---------------------

    mwRTM_param%vegcls      = nint(scalar)
    mwRTM_param%soilcls     = nint(scalar)

    mwRTM_param%sand        = scalar
    mwRTM_param%clay        = scalar
    mwRTM_param%poros       = scalar
    mwRTM_param%wang_wt     = scalar
    mwRTM_param%wang_wp     = scalar
    mwRTM_param%rgh_hmin    = scalar
    mwRTM_param%rgh_hmax    = scalar
    mwRTM_param%rgh_wmin    = scalar
    mwRTM_param%rgh_wmax    = scalar
    mwRTM_param%rgh_Nrh     = scalar
    mwRTM_param%rgh_Nrv     = scalar
    mwRTM_param%rgh_polmix  = scalar
    mwRTM_param%omega       = scalar
    mwRTM_param%bh          = scalar
    mwRTM_param%bv          = scalar
    mwRTM_param%lewt        = scalar
    
  end subroutine scalar2mwRTM_param

  ! ************************************************************
  
  subroutine mwRTM_param_nodata_check( mwp, is_nodata )
    
    implicit none
    
    type(mwRTM_param_type), intent(inout) :: mwp
    
    logical,                intent(  out) :: is_nodata
    
    ! local variables
    
    real :: realvegcls, realsoilcls
    
    ! ---------------------
    
    realvegcls  = real(mwp%vegcls)
    realsoilcls = real(mwp%soilcls)
    
    if ( (abs(realvegcls    -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(realsoilcls   -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%sand      -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%clay      -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%poros     -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%wang_wt   -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%wang_wp   -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%rgh_hmin  -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%rgh_hmax  -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%rgh_wmin  -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%rgh_wmax  -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%rgh_Nrh   -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%rgh_Nrv   -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%rgh_polmix-nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%omega     -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%bh        -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%bv        -nodata_generic)<nodata_tol_generic) .or.        &
         (abs(mwp%lewt      -nodata_generic)<nodata_tol_generic)      ) then
       
       mwp = nodata_generic
       
       is_nodata = .true.

    else
       
       is_nodata = .false.
       
    end if
    
  end subroutine mwRTM_param_nodata_check
  
end module mwRTM_types
  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if 0

! driver routines for testing

program test_mwRTM_types
  
  use mwRTM_types
  
  implicit none
  
  type(mwRTM_param_type) :: mwRTM_param
  
  mwRTM_param = -9999.
  
  write (*,*) mwRTM_param

  mwRTM_param%lewt = 0.5

  write (*,*) mwRTM_param

  call mwRTM_param_nodata_check(mwRTM_param)
  
  write (*,*) mwRTM_param  

end program test_mwRTM_types

#endif

! ========================== EOF ==================================

