
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
       nodata_generic,                         &
       nodata_tol_generic,                     &
       LDAS_is_nodata

  use ldas_exceptionsMod,            ONLY:     &
       ldas_abort,                             &
       LDAS_GENERIC_ERROR

  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: mwRTM_param_type
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

     ! parameters relating to vegetation opacity 
     
     real    :: bh          ! veg b parameter (h-pol)  (tau = b*VWC)        [dim-less]
     real    :: bv          ! veg b parameter (v-pol)  (tau = b*VWC)        [dim-less]
     real    :: lewt        ! VWC = lewt*LAI                                [kg/m2]
     real    :: vegopacity  ! veg opacity = tau/cos(inc_angle)              [dim-less]
     
  end type mwRTM_param_type
  
  ! ---------------------------------------------------------
  
  interface assignment (=)
     module procedure scalar2mwRTM_param
  end interface
  
contains

  ! **********************************************************************

  ! Subroutine io_mwRTM_param_type() reads and writes binary mwRTM files and is no longer used.
  !
  ! The subroutine has been replaced by:
  ! - Applications/LDAS_App/mwrtm_bin2nc4.F90 converts mwRTM files from binary to nc4
  ! - get_mwrtm_param() in GEOS_LandAssimGridComp.F90 converts the internal state
  !    variables of the Land Assim GridComp into the mwRTM structure.
  !
  ! reichle, 4 Aug 2020

!  subroutine io_mwRTM_param_type( action, unitnum, N_tile, mwp, N_catg, tile_id, d2g ) 
!    
!    ! read/write mwRTM_param for domain from/to file
!    !
!    ! write: write mwRTM params for N_tile tiles in domain
!    !
!    ! read:  read mwRTM params for 
!    
!    !
!    ! reichle, 1  Jun 2011
!    ! reichle, 21 Oct 2011 -- added "read" for mwRTM params
!    ! reichle, 22 Oct 2012 -- added check for tile_id to mwRTM param input
!    ! 
!    ! -------------------------------------------------------------------
!    
!    implicit none
!    
!    character,                                 intent(in)    :: action
!       
!    integer,                                   intent(in)    :: unitnum
!    
!    integer,                                   intent(in)    :: N_tile  ! =N_catd
!    
!    type(mwRTM_param_type), dimension(N_tile), intent(inout) :: mwp
!
!    integer, optional,                         intent(in)    :: N_catg    
!
!    integer, optional,      dimension(N_tile), intent(in)    :: tile_id, d2g
!
!    ! local variables
!    
!    integer :: n, N_tmp
!
!    integer, dimension(:), allocatable :: tmpint
!    real,    dimension(:), allocatable :: tmpreal
!
!    character(len=*), parameter :: Iam = 'io_mwRTM_param_type'
!    character(len=400) :: err_msg
!
!    ! ------------------------------------------------------------------
!
!    select case (action)
!       
!    case ('w','W')     ! write mwp for all tiles in domain
!       
!       write (unitnum) N_tile
!       
!       write (unitnum) (mwp(n)%vegcls    ,    n=1,N_tile)  ! integer
!       write (unitnum) (mwp(n)%soilcls   ,    n=1,N_tile)  ! integer
!       
!       write (unitnum) (mwp(n)%sand      ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%clay      ,    n=1,N_tile)  ! real 
!
!       write (unitnum) (mwp(n)%poros     ,    n=1,N_tile)  ! real 
!       
!       write (unitnum) (mwp(n)%wang_wt   ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%wang_wp   ,    n=1,N_tile)  ! real 
!       
!       write (unitnum) (mwp(n)%rgh_hmin  ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%rgh_hmax  ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%rgh_wmin  ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%rgh_wmax  ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%rgh_Nrh   ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%rgh_Nrv   ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%rgh_polmix,    n=1,N_tile)  ! real 
!       
!       write (unitnum) (mwp(n)%omega     ,    n=1,N_tile)  ! real 
!       
!       write (unitnum) (mwp(n)%bh        ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%bv        ,    n=1,N_tile)  ! real 
!       write (unitnum) (mwp(n)%lewt      ,    n=1,N_tile)  ! real 
!       
!    case ('r','R')
!
!       ! read the parameters for all global tiles (similar to read_land_parameters())
!       
!       ! optional inputs N_catg and d2g must be present
!       
!       if ( (.not. present(N_catg )) .or.      &
!            (.not. present(tile_id)) .or.      &
!            (.not. present(d2g    ))         ) then
!          err_msg = 'missing optional inputs N_catg, tile_id, d2g'
!          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
!       end if
!       
!       ! read how many tiles are in file and double-check against N_catg
!       
!       read (unitnum) N_tmp
!       
!       if (N_tmp .ne. N_catg) then
!          err_msg = 'number of tiles in file .ne. N_catg'
!          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
!       end if
!       
!       ! allocate tmp vectors
!       
!       allocate(tmpint( N_catg))
!       allocate(tmpreal(N_catg))
!       
!       ! read tile IDs (first record)
!       
!       read (unitnum) tmpint;
!
!       ! make sure tile IDs match (works only for "SiB2_V2" and newer versions)
!       
!       if (any(tile_id/=tmpint(d2g(1:N_tile)))) then
!          err_msg = 'mismatch of tile IDs for mwRTM_parameters'
!          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
!       end if
!       
!       ! read and subset from global tile space to domain
!       
!       read (unitnum) tmpint;   mwp(1:N_tile)%vegcls     =  tmpint( d2g(1:N_tile))
!       read (unitnum) tmpint;   mwp(1:N_tile)%soilcls    =  tmpint( d2g(1:N_tile))
!       
!       read (unitnum) tmpreal;  mwp(1:N_tile)%sand       =  tmpreal(d2g(1:N_tile))
!       read (unitnum) tmpreal;  mwp(1:N_tile)%clay       =  tmpreal(d2g(1:N_tile))
!       
!       read (unitnum) tmpreal;  mwp(1:N_tile)%poros      =  tmpreal(d2g(1:N_tile))
!       
!       read (unitnum) tmpreal;  mwp(1:N_tile)%wang_wt    =  tmpreal(d2g(1:N_tile)) 
!       read (unitnum) tmpreal;  mwp(1:N_tile)%wang_wp    =  tmpreal(d2g(1:N_tile))
!       
!       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_hmin   =  tmpreal(d2g(1:N_tile)) 
!       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_hmax   =  tmpreal(d2g(1:N_tile))
!       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_wmin   =  tmpreal(d2g(1:N_tile))
!       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_wmax   =  tmpreal(d2g(1:N_tile))
!       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_Nrh    =  tmpreal(d2g(1:N_tile))
!       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_Nrv    =  tmpreal(d2g(1:N_tile))
!       read (unitnum) tmpreal;  mwp(1:N_tile)%rgh_polmix =  tmpreal(d2g(1:N_tile))    
!       
!       read (unitnum) tmpreal;  mwp(1:N_tile)%omega      =  tmpreal(d2g(1:N_tile))
!       
!       read (unitnum) tmpreal;  mwp(1:N_tile)%bh         =  tmpreal(d2g(1:N_tile))
!       read (unitnum) tmpreal;  mwp(1:N_tile)%bv         =  tmpreal(d2g(1:N_tile))
!       read (unitnum) tmpreal;  mwp(1:N_tile)%lewt       =  tmpreal(d2g(1:N_tile))
!       
!       ! clean up
!
!       deallocate(tmpreal)
!       deallocate(tmpint)
!       
!    case default
!       
!       err_msg = 'unknown action ' // action
!       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
!       
!    end select
!    
!  end subroutine io_mwRTM_param_type
  
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
    mwRTM_param%vegopacity  = scalar
    
  end subroutine scalar2mwRTM_param

  ! ************************************************************
  
  subroutine mwRTM_param_nodata_check( mwp, mwp_nodata )

    ! check microwave radiative transfer model parameters for no-data values
    !
    ! if there is a no-data value in any required field, set "all" fields
    !   within the corresponding group of parameters to no-data 
    !
    ! vegetation attenuation parameters can come from either of two sources:
    ! - static look-up table or calibrated parameters (bh, bv, lewt), to be combined
    !   with time-varying LAI
    ! - vegetation opacity (from file); varies with time
    !
    ! preprocessing of the mwRTM restart and vegopacity files must ensure that 
    !   either (bh,bv,lewt) or (vegopacity) is no-data
    !
    ! - reichle, 13 July 2021  (revised for using vegopacity from file)
    
    implicit none
    
    type(mwRTM_param_type), intent(inout) :: mwp
    
    logical,                intent(  out) :: mwp_nodata
    
    ! local variables
    
    logical            :: veg_LUT_params_nodata, veg_params_nodata, other_params_nodata

    real               :: realvegcls, realsoilcls
    
    character(len=400) :: err_msg

    ! -----------------------------------------------------------------------------

    ! Group 1: Parameters related to vegetation attenuation
    !
    ! need either (bh, bv, lewt) or (vegopacity)
    !
    ! check if static look-up table (LUT) or calibrated parameters are available
    
    veg_LUT_params_nodata =                           &
         (                                            &
         LDAS_is_nodata( mwp%bh         ) .or.        &       
         LDAS_is_nodata( mwp%bv         ) .or.        &
         LDAS_is_nodata( mwp%lewt       )             &
         )
    
    if ( (.not. veg_LUT_params_nodata) .and. (.not. LDAS_is_nodata( mwp%vegopacity )) ) then
       
       ! inconsistent mwRTM restart and vegopacity files: 
       ! for a given tile, (bh, bv, lewt) from mwRTM restart and (vegopacity) from file must
       ! not both have good values
       
       err_msg = 'inconsistent mwRTM restart and vegopacity files: found good values for (bh,bv,lewt) *and* (vegopacity)'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    if ( veg_LUT_params_nodata ) then  
       
       ! make sure all related fields in this group are no-data if at least one is nodata
       
       mwp%bh   = nodata_generic
       mwp%bv   = nodata_generic
       mwp%lewt = nodata_generic
       
    end if
         
    ! veg_params_nodata = .true. if missing vegetation attenuation info altogether
    
    veg_params_nodata =                               &
         (                                            &
         veg_LUT_params_nodata                        &
         .and.                                        &
         LDAS_is_nodata( mwp%vegopacity )             &  
         )

    ! -----------------------------------------------------------------------------

    ! Group 2: Parameters for the rest of the tau-omega equations
    
    realvegcls  = real(mwp%vegcls)
    realsoilcls = real(mwp%soilcls)

    other_params_nodata =                             &
         (                                            &
         LDAS_is_nodata( realvegcls     ) .or.        &
         LDAS_is_nodata( realsoilcls    ) .or.        &
         LDAS_is_nodata( mwp%sand       ) .or.        &
         LDAS_is_nodata( mwp%clay       ) .or.        &
         LDAS_is_nodata( mwp%poros      ) .or.        &
         LDAS_is_nodata( mwp%wang_wt    ) .or.        &
         LDAS_is_nodata( mwp%wang_wp    ) .or.        &
         LDAS_is_nodata( mwp%rgh_hmin   ) .or.        &
         LDAS_is_nodata( mwp%rgh_hmax   ) .or.        &
         LDAS_is_nodata( mwp%rgh_wmin   ) .or.        &
         LDAS_is_nodata( mwp%rgh_wmax   ) .or.        &
         LDAS_is_nodata( mwp%rgh_Nrh    ) .or.        &
         LDAS_is_nodata( mwp%rgh_Nrv    ) .or.        &
         LDAS_is_nodata( mwp%rgh_polmix ) .or.        &
         LDAS_is_nodata( mwp%omega      )             &
         )                               

    if ( other_params_nodata ) then
       
       ! make sure all related fields in this group are no-data if at least one is nodata
       
       mwp%realvegcls     = nodata_generic
       mwp%realsoilcls    = nodata_generic
       mwp%mwp%sand       = nodata_generic
       mwp%mwp%clay       = nodata_generic
       mwp%mwp%poros      = nodata_generic
       mwp%mwp%wang_wt    = nodata_generic
       mwp%mwp%wang_wp    = nodata_generic
       mwp%mwp%rgh_hmin   = nodata_generic
       mwp%mwp%rgh_hmax   = nodata_generic
       mwp%mwp%rgh_wmin   = nodata_generic
       mwp%mwp%rgh_wmax   = nodata_generic
       mwp%mwp%rgh_Nrh    = nodata_generic
       mwp%mwp%rgh_Nrv    = nodata_generic
       mwp%mwp%rgh_polmix = nodata_generic
       mwp%mwp%omega      = nodata_generic
       
    end if
    
    ! -----------------------------------------------------------------------------
    
    ! need both groups for full tau-omega calculations:
    
    if ( veg_params_nodata .or. other_params_nodata ) then
       
       mwp_nodata = .true.
       
    else
       
       mwp_nodata = .false.
       
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

