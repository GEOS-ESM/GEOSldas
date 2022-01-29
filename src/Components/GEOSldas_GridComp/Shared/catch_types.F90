
module catch_types

  ! definition of types and associated operators for Catchment Model
  !
  ! IMPORTANT:
  ! When adding a field to any of the derived types, must also update
  ! the associated assignment and operator definitions.
  ! THERE IS NO WARNING/ERROR IF OPERATOR IS NOT DEFINED FOR ALL FIELDS!
  !
  ! reichle, 21 May 2003
  ! reichle, 25 Jan 2005 - added cat_force_type
  ! reichle, 28 Oct 2010 - added soilcls30 and soilcls100
  ! reichle,  9 Dec 2011 - removed water/energy balance terms from cat_diagn
  !                        (now done in new "bal_diagn_type" in driver_types.F90)
  ! reichle, 28 Dec 2011 - removed field totalb from cat_diagn structure
  !                        (now done via swup/SWdown)
  ! reichle, 30 Oct 2013 - removed field rzeq from cat_diagn structure
  ! reichle, 31 Oct 2013 - split "cat_diagn" structure into "cat_diagS" and "cat_diagF"
  ! reichle, 16 Nov 2015 - added vegetation height
  !
  ! --------------------------------------------------------------------------
  
  use catch_constants,                  ONLY:     &
       N_snow        => CATCH_N_SNOW,             &
       N_gt          => CATCH_N_GT     
  
  implicit none
  
  ! everything is private by default unless made public
  
  private

    
  public :: N_cat_progn,    N_cat_diagS,    N_cat_diagF
  public :: cat_progn_type, cat_diagS_type, cat_diagF_type
  public :: cat_param_type, cat_force_type
  
  public :: assignment (=), operator (/), operator (+)

  public :: catprogn2wesn, catprogn2htsn, catprogn2sndz, catprogn2ghtcnt
  ! -------------------------------------------------------------------------
  !
  ! N_cat_progn = total # states in Catchment model (including 3*N_snow states 
  ! for snow and N_gt states for ground temperature)
  ! N_snow        => CATCH_N_SNOW,             &
  ! N_gt          => CATCH_N_GT
  !integer,parameter :: N_gt   = 6
  !integer,parameter :: N_snow = 3
 
  integer, parameter :: N_cat_progn = 10 + N_gt + 3*N_snow

  ! --------------------------------------------------------------------------
  
  ! Catchment model prognostic variables
  
  type :: cat_progn_type

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     ! horizontally, the surface is divided into four fractions:
     !
     ! "1" - saturated
     ! "2" - unsaturated but not stressed
     ! "4" - stressed
     ! "S" - snow
     !
     ! ------------------------------------------------------------
     
     real :: tc1     ! surface/canopy temperature 
     real :: tc2
     real :: tc4
     
     real :: qa1     ! specific humidity in canopy air
     real :: qa2
     real :: qa4
     
     real :: capac   ! canopy interception water
     
     real :: catdef  ! catchment deficit
     real :: rzexc   ! root zone excess
     real :: srfexc  ! surface excess
     
     real, dimension(N_gt)   :: ght     ! ground heat content
     
     real, dimension(N_snow) :: wesn    ! snow water equivalent
     real, dimension(N_snow) :: htsn    ! snow heat content
     real, dimension(N_snow) :: sndz    ! snow depth
  
  end type cat_progn_type
  
  ! ---------------------------------------------------------
  
  ! Catchment model diagnostic variables

  ! Catchment model diagnostics are split into two groups:
  !
  ! cat_diagS = diagnostic "state" variables that can be computed from prognostics
  ! cat_diagF = diagnostic variables such as "fluxes" that are outputs of subroutine 
  !              catchment() but cannot be computed directly from prognostics only
  
  integer, parameter :: N_cat_diagS = 7 + N_gt + N_snow
  
  type :: cat_diagS_type
     
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     real :: ar1      ! area fraction of saturated zone
     real :: ar2      ! area fraction of unsaturated and unstressed zone

     real :: asnow    ! area fraction of snow
     
     real :: sfmc     ! surface moisture content
     real :: rzmc     ! root zone moisture content
     real :: prmc     ! profile moisture content
     
     real :: tsurf    ! mean surface temperature over entire catchment 

     real, dimension(N_gt)   :: tp     ! temperature of soil layers
     
     real, dimension(N_snow) :: tpsn   ! temperature of snow layers
     
  end type cat_diagS_type
  
  ! --------------------------

  integer, parameter :: N_cat_diagF = 22
  
  type :: cat_diagF_type
     
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     real :: shflux   ! sensible heat flux
     real :: lhflux   ! total latent heat flux
     real :: ghflux   ! ground heat flux to top soil layer
     
     real :: evap     ! total evaporation 
     real :: eint     ! interception loss
     real :: esoi     ! evaporation from bare soil
     real :: eveg     ! transpiration 
     real :: esno     ! evaporation from snow
     
     real :: runoff   ! total runoff
     real :: runsrf   ! surface runoff
     real :: bflow    ! baseflow
     
     real :: snmelt   ! snow melt
     
     real :: lwup     ! outgoing/upward longwave radiation
     real :: swup     ! outgoing/upward shortwave radiation
     
     real :: qinfil   ! infiltration

     real :: hsnacc   ! accounting term for energy related to snowfall etc.
     real :: evacc    ! accounting term for evaporation   (see catchment()) 
     real :: shacc    ! accounting term for sensible heat (see catchment()) 
     real :: lhacc    ! accounting term for latent heat   (see catchment())
     real :: eacc_0   ! accounting term for oscillations  (see catchment())  

     ! t2m and q2m depend on fluxes and cannot be computed from prognostics only

     real :: t2m      ! air temperature at 2m above the displacement height
     real :: q2m      ! specific humidity at 2m above the displacement height
          
  end type cat_diagF_type
  
  ! ---------------------------------------------------------
  
  ! Catchment model parameters
  
  type :: cat_param_type
     
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     real :: dpth  ! depth to bedrock from data file (dpth/=dzpr in general!)      [mm]
     
     ! layer thicknesses for soil moisture model (in [mm]!!!!) 
     
     real :: dzsf  ! "surface layer"                     formerly zdep1            [mm]
     real :: dzrz  ! "root zone layer"                   formerly zdep2            [mm]
     real :: dzpr  ! "profile layer" (unsaturated zone)  formerly zdep3            [mm]
     
     ! layer thicknesses for ground temperature model (in [m]!!!!)
     !
     ! dzgt SHOULD REPLACE data dz /.../ STATEMENT IN gndtp0 AND gndtmp
     !
     real, dimension(N_gt) :: dzgt            !                                    [m]      
                                                                                         
     ! soil hydraulic parameters                                                         
                                                                                         
     real :: poros   ! porosity                                                    [m3 m-3]  
     real :: cond    ! saturated hydraulic conductivity                            [m s-1]    
     real :: psis    ! Clapp-Hornberger parameter                                  [m H2O]  
     real :: bee     ! Clapp-Hornberger parameter                                  [-]      
                                                                                         
     real :: wpwet   ! wilting poing wetness                                       [-]      
                                                                                         
     real :: gnu     ! vertical decay factor for transmissivity                    [m-1]    
     
     ! constant parameters related to vegetation
     
     real :: vgwmax      ! max amount of water available to vegetation             [kg m-2]

     ! veg and soil classes
     
     integer :: vegcls         ! vegetation class                                  [-]
     integer :: soilcls30      ! soil_class_top (0- 30cm)                          [-]
     integer :: soilcls100     ! soil_class_com (0-100cm)                          [-]
          
     ! parameters specific to Catchment Model
     ! (Equation and Figure numbers refer to Ducharne et al., 2000, doi:10.1029/2000JD900328)
     
     real :: bf1       ! baseflow parameter (A in Eq 9)                            [kg m-4]
     real :: bf2       ! baseflow parameter (B in Eq 9)                            [m]
     real :: bf3       ! baseflow parameter (XBAR in Eq 8)                         [log(m)]

     real :: cdcr1     ! catdef threshold (water table at bedrock)                 [kg m-2]
     real :: cdcr2     ! catdef threshold ()                                       [kg m-2]

     ! area partitioning parameters

     real :: ars1 ! A          in Eq 12 for Asat                                   [m2 kg-1]                              
     real :: ars2 ! B          in Eq 12 for Asat                                   [m2 kg-1]   
     real :: ars3 ! C          in Eq 12 for Asat                                   [m4 kg-2]   
     real :: ara1 ! A          in Eq 14 of segment1 if skew < 0.25 else ara1=ara3  [m2 kg-1]
     real :: ara2 ! B          in Eq 14 of segment1 if skew < 0.25 else ara2=ara4  [-]       
     real :: ara3 ! A          in Eq 14 of segment1 if skew < 0.25                 [m2 kg-1]
     real :: ara4 ! B          in Eq 14 of segment1 if skew < 0.25                 [-]                                            
     real :: arw1 ! A          in Eq 12 for THETA0                                 [m2 kg-1]        
     real :: arw2 ! B          in Eq 12 for THETA0                                 [m2 kg-1] 
     real :: arw3 ! C          in Eq 12 for THETA0                                 [m4 kg-2]                         
     real :: arw4 ! Y_infinity in Eq 12 for THETA0                                 [-]                       

     ! time scale param for moisture transfer between root zone and water table 
     
     real :: tsa1 ! atau1 in Eq 16 for root zone excess > 0 (Fig 6)                [-] 
     real :: tsa2 ! atau1 in Eq 16 for root zone excess < 0 (Fig 6)                [-] 
     real :: tsb1 ! btau1 in Eq 16 for root zone excess > 0 (Fig 6)                [-] 
     real :: tsb2 ! btau1 in Eq 16 for root zone excess < 0 (Fig 6)                [-] 

     ! time scale param for moisture transfer between surface excess and root zone excess

     real :: atau ! atau2 in Eq 17                                                 [-]
     real :: btau ! btau2 in Eq 17                                                 [-]
     
     ! additional soil parameters from recent versions of "soil_param.dat"
     ! (eg. for use in calibration of the microwave radiative transfer model)
     ! - reichle,  1 Apr 2015

     real :: gravel30 ! gravel                in 0- 30cm layer                     [percent by vol]
     real :: orgC30   ! organic carbon        in 0- 30cm layer                     [percent by weight] 
     real :: orgC     ! organic carbon        in 0-100cm layer                     [percent by weight]
     real :: sand30   ! sand fraction         in 0- 30cm layer                     [percent by weight]
     real :: clay30   ! clay fraction         in 0- 30cm layer                     [percent by weight]
     real :: sand     ! sand fraction         in 0-100cm layer                     [percent by weight]
     real :: clay     ! clay fraction         in 0-100cm layer                     [percent by weight]
     real :: wpwet30  ! wilting point wetness in 0- 30cm layer                     [-]
     real :: poros30  ! porosity              in 0- 30cm layer                     [m3 m-3]

     ! static (time-invariant) vegetation parameters

     real :: veghght  ! vegetation height                                          [m]

  end type cat_param_type
  
  ! ---------------------------------------------------------
  !
  ! input forcings (or boundary conditions) and related variables 
  !
  ! horizontally, the surface is divided into four fractions:
  !
  ! "1" - saturated
  ! "2" - unsaturated but not stressed
  ! "4" - stressed
  ! "S" - snow
  
  type :: cat_force_type
     
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     real :: TRAINC        ! convective rain rate
     real :: TRAINL        ! large-scale rain rate
     real :: TSNOW         ! snowfall
     real :: UM            ! wind
     real :: ETURB1 
     real :: DEDQA1 
     real :: DEDTC1 
     real :: HSTURB1
     real :: DHSDQA1 
     real :: DHSDTC1
     real :: ETURB2 
     real :: DEDQA2 
     real :: DEDTC2 
     real :: HSTURB2
     real :: DHSDQA2 
     real :: DHSDTC2
     real :: ETURB4 
     real :: DEDQA4 
     real :: DEDTC4 
     real :: HSTURB4
     real :: DHSDQA4 
     real :: DHSDTC4
     real :: ETURBS 
     real :: DEDQAS 
     real :: DEDTCS 
     real :: HSTURBS
     real :: DHSDQAS 
     real :: DHSDTCS
     real :: TM            ! 2m temperature
     real :: QM            ! 2m humidity
     real :: ra1 
     real :: ra2 
     real :: ra4 
     real :: raS 
     real :: SUNANG        ! sun angle
     real :: PARDIR        ! direct photosynthetically active radiation
     real :: PARDIF        ! diffuse photosynthetically active radiation
     real :: SWNETF        ! net shortwave radiation (?)
     real :: SWNETS        ! net shortwave radiation (?)
     real :: HLWDWN        ! downward longwave radiation
     real :: PSUR          ! surface pressure
     real :: ZLAI          ! leaf area index
     real :: GREEN         ! greenness
     real :: Z2            
     real :: SQSCAT 
     real :: RSOIL1 
     real :: RSOIL2   
     real :: RDC  
     real :: QSAT1 
     real :: DQS1 
     real :: ALW1 
     real :: BLW1
     real :: QSAT2 
     real :: DQS2 
     real :: ALW2 
     real :: BLW2
     real :: QSAT4 
     real :: DQS4 
     real :: ALW4 
     real :: BLW4
     real :: QSATS 
     real :: DQSS 
     real :: ALWS 
     real :: BLWS                  

  end type cat_force_type
  
  ! ----------------------------------------------------------------
  
  interface assignment (=)
     module procedure scalar2cat_progn
     module procedure scalar2cat_diagS
     module procedure scalar2cat_diagF
     module procedure scalar2cat_param
     module procedure scalar2cat_force
  end interface
  
  interface operator (/)
     module procedure cat_progn_div_scalar
     module procedure cat_diagS_div_scalar
     module procedure cat_diagF_div_scalar
     module procedure cat_force_div_scalar
  end interface

  interface operator (+)
     module procedure add_cat_progn
     module procedure add_cat_diagS
     module procedure add_cat_diagF
     module procedure add_cat_force
  end interface

  ! ----------------------------------------------------------------  

contains
  
  subroutine scalar2cat_diagS( cat_diagS, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(cat_diagS_type), intent(out) :: cat_diagS
    
    integer :: i     ! local

    cat_diagS%ar1    = scalar
    cat_diagS%ar2    = scalar

    cat_diagS%asnow  = scalar
    
    cat_diagS%sfmc   = scalar
    cat_diagS%rzmc   = scalar
    cat_diagS%prmc   = scalar    

    cat_diagS%tsurf  = scalar

    do i=1,N_gt
       cat_diagS%tp(i)  = scalar
    end do

    do i=1,N_snow
       cat_diagS%tpsn(i)= scalar
    end do

  end subroutine scalar2cat_diagS
  
  ! -----------------------------------------------------------

  subroutine scalar2cat_diagF( cat_diagF, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(cat_diagF_type), intent(out) :: cat_diagF

    cat_diagF%shflux = scalar
    cat_diagF%lhflux = scalar
    cat_diagF%ghflux = scalar

    cat_diagF%evap   = scalar
    cat_diagF%eint   = scalar
    cat_diagF%esoi   = scalar
    cat_diagF%eveg   = scalar
    cat_diagF%esno   = scalar
    
    cat_diagF%runoff = scalar
    cat_diagF%runsrf = scalar
    cat_diagF%bflow  = scalar

    cat_diagF%snmelt = scalar

    cat_diagF%lwup   = scalar
    cat_diagF%swup   = scalar
    
    cat_diagF%qinfil = scalar
        
    cat_diagF%hsnacc = scalar 
    cat_diagF%evacc  = scalar 
    cat_diagF%shacc  = scalar 
    cat_diagF%lhacc  = scalar
    cat_diagF%eacc_0 = scalar

    cat_diagF%t2m    = scalar
    cat_diagF%q2m    = scalar

  end subroutine scalar2cat_diagF
  
  ! -----------------------------------------------------------

  function cat_diagS_div_scalar( cat_diagS, scalar )
    
    implicit none

    type(cat_diagS_type)             :: cat_diagS_div_scalar
    type(cat_diagS_type), intent(in) :: cat_diagS
    
    real, intent(in) :: scalar
    
    integer :: i     ! local
    
    cat_diagS_div_scalar%ar1    =     cat_diagS%ar1    / scalar
    cat_diagS_div_scalar%ar2    =     cat_diagS%ar2    / scalar

    cat_diagS_div_scalar%asnow  =     cat_diagS%asnow  / scalar

    cat_diagS_div_scalar%sfmc   =     cat_diagS%sfmc   / scalar
    cat_diagS_div_scalar%rzmc   =     cat_diagS%rzmc   / scalar
    cat_diagS_div_scalar%prmc   =     cat_diagS%prmc   / scalar
    
    cat_diagS_div_scalar%tsurf  =     cat_diagS%tsurf  / scalar
    
    do i=1,N_gt
       cat_diagS_div_scalar%tp(i)  =     cat_diagS%tp(i)  / scalar
    end do

    do i=1,N_snow
       cat_diagS_div_scalar%tpsn(i)=     cat_diagS%tpsn(i)/ scalar
    end do
 
  end function cat_diagS_div_scalar

  ! -----------------------------------------------------------

  function cat_diagF_div_scalar( cat_diagF, scalar )
    
    implicit none

    type(cat_diagF_type)             :: cat_diagF_div_scalar
    type(cat_diagF_type), intent(in) :: cat_diagF
    
    real, intent(in) :: scalar
    
    cat_diagF_div_scalar%shflux = cat_diagF%shflux / scalar
    cat_diagF_div_scalar%lhflux = cat_diagF%lhflux / scalar
    cat_diagF_div_scalar%ghflux = cat_diagF%ghflux / scalar

    cat_diagF_div_scalar%evap   = cat_diagF%evap   / scalar
    cat_diagF_div_scalar%eint   = cat_diagF%eint   / scalar
    cat_diagF_div_scalar%esoi   = cat_diagF%esoi   / scalar
    cat_diagF_div_scalar%eveg   = cat_diagF%eveg   / scalar
    cat_diagF_div_scalar%esno   = cat_diagF%esno   / scalar

    
    cat_diagF_div_scalar%runoff = cat_diagF%runoff / scalar
    cat_diagF_div_scalar%runsrf = cat_diagF%runsrf / scalar
    cat_diagF_div_scalar%bflow  = cat_diagF%bflow  / scalar

    cat_diagF_div_scalar%snmelt = cat_diagF%snmelt / scalar
    
    cat_diagF_div_scalar%lwup   = cat_diagF%lwup   / scalar
    cat_diagF_div_scalar%swup   = cat_diagF%swup   / scalar

    cat_diagF_div_scalar%qinfil = cat_diagF%qinfil / scalar
    
    cat_diagF_div_scalar%hsnacc = cat_diagF%hsnacc / scalar 
    cat_diagF_div_scalar%evacc  = cat_diagF%evacc  / scalar 
    cat_diagF_div_scalar%shacc  = cat_diagF%shacc  / scalar 
    cat_diagF_div_scalar%lhacc  = cat_diagF%lhacc  / scalar 
    cat_diagF_div_scalar%eacc_0 = cat_diagF%eacc_0 / scalar 
   
    cat_diagF_div_scalar%t2m    = cat_diagF%t2m    / scalar
    cat_diagF_div_scalar%q2m    = cat_diagF%q2m    / scalar

    
  end function cat_diagF_div_scalar

  ! -----------------------------------------------------------

  function add_cat_diagS( cat_diagS_1, cat_diagS_2 )
    
    implicit none

    type(cat_diagS_type)             :: add_cat_diagS
    type(cat_diagS_type), intent(in) :: cat_diagS_1, cat_diagS_2

    integer :: i     ! local
    
    add_cat_diagS%ar1    =     cat_diagS_1%ar1    +     cat_diagS_2%ar1    
    add_cat_diagS%ar2    =     cat_diagS_1%ar2    +     cat_diagS_2%ar2    
    
    add_cat_diagS%asnow  =     cat_diagS_1%asnow  +     cat_diagS_2%asnow  
    
    add_cat_diagS%sfmc   =     cat_diagS_1%sfmc   +     cat_diagS_2%sfmc  
    add_cat_diagS%rzmc   =     cat_diagS_1%rzmc   +     cat_diagS_2%rzmc  
    add_cat_diagS%prmc   =     cat_diagS_1%prmc   +     cat_diagS_2%prmc  

    add_cat_diagS%tsurf  =     cat_diagS_1%tsurf  +     cat_diagS_2%tsurf 

    do i=1,N_gt
       add_cat_diagS%tp(i)  =     cat_diagS_1%tp(i)  +     cat_diagS_2%tp(i)  
    end do
    
    do i=1,N_snow
       add_cat_diagS%tpsn(i)=     cat_diagS_1%tpsn(i)+     cat_diagS_2%tpsn(i)
    end do

  end function add_cat_diagS

  ! -----------------------------------------------------------

  function add_cat_diagF( cat_diagF_1, cat_diagF_2 )
    
    implicit none

    type(cat_diagF_type)             :: add_cat_diagF
    type(cat_diagF_type), intent(in) :: cat_diagF_1, cat_diagF_2

    add_cat_diagF%shflux = cat_diagF_1%shflux + cat_diagF_2%shflux 
    add_cat_diagF%lhflux = cat_diagF_1%lhflux + cat_diagF_2%lhflux
    add_cat_diagF%ghflux = cat_diagF_1%ghflux + cat_diagF_2%ghflux 

    add_cat_diagF%evap   = cat_diagF_1%evap   + cat_diagF_2%evap   
    add_cat_diagF%eint   = cat_diagF_1%eint   + cat_diagF_2%eint   
    add_cat_diagF%esoi   = cat_diagF_1%esoi   + cat_diagF_2%esoi   
    add_cat_diagF%eveg   = cat_diagF_1%eveg   + cat_diagF_2%eveg   
    add_cat_diagF%esno   = cat_diagF_1%esno   + cat_diagF_2%esno   

    add_cat_diagF%runoff = cat_diagF_1%runoff + cat_diagF_2%runoff 
    add_cat_diagF%runsrf = cat_diagF_1%runsrf + cat_diagF_2%runsrf 
    add_cat_diagF%bflow  = cat_diagF_1%bflow  + cat_diagF_2%bflow  
    
    add_cat_diagF%snmelt = cat_diagF_1%snmelt + cat_diagF_2%snmelt  

    add_cat_diagF%lwup   = cat_diagF_1%lwup   + cat_diagF_2%lwup  
    add_cat_diagF%swup   = cat_diagF_1%swup   + cat_diagF_2%swup  

    add_cat_diagF%qinfil = cat_diagF_1%qinfil + cat_diagF_2%qinfil 

    add_cat_diagF%hsnacc = cat_diagF_1%hsnacc + cat_diagF_2%hsnacc 
    add_cat_diagF%evacc  = cat_diagF_1%evacc  + cat_diagF_2%evacc 
    add_cat_diagF%shacc  = cat_diagF_1%shacc  + cat_diagF_2%shacc 
    add_cat_diagF%lhacc  = cat_diagF_1%lhacc  + cat_diagF_2%lhacc 
    add_cat_diagF%eacc_0 = cat_diagF_1%eacc_0 + cat_diagF_2%eacc_0 

    add_cat_diagF%t2m    = cat_diagF_1%t2m    + cat_diagF_2%t2m
    add_cat_diagF%q2m    = cat_diagF_1%q2m    + cat_diagF_2%q2m

    
  end function add_cat_diagF

  ! *******************************************************************
  
  subroutine scalar2cat_progn( cat_progn, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(cat_progn_type), intent(out) :: cat_progn
    
    integer :: i     ! local
    
    cat_progn%tc1    = scalar
    cat_progn%tc2    = scalar
    cat_progn%tc4    = scalar
    cat_progn%qa1    = scalar
    cat_progn%qa2    = scalar
    cat_progn%qa4    = scalar
    cat_progn%capac  = scalar
    cat_progn%catdef = scalar
    cat_progn%rzexc  = scalar
    cat_progn%srfexc = scalar
    
    do i=1,N_gt
       cat_progn%ght(i)  = scalar
    end do
    
    do i=1,N_snow
       cat_progn%wesn(i) = scalar
       cat_progn%htsn(i) = scalar
       cat_progn%sndz(i) = scalar
    end do
    
  end subroutine scalar2cat_progn
  
  ! ---------------------------------------------------
  
  function cat_progn_div_scalar( cat_progn, scalar )
    
    implicit none

    type(cat_progn_type)             :: cat_progn_div_scalar
    type(cat_progn_type), intent(in) :: cat_progn
    
    real, intent(in) :: scalar
    
    integer :: i       ! local

    cat_progn_div_scalar%tc1    =     cat_progn%tc1    / scalar
    cat_progn_div_scalar%tc2    =     cat_progn%tc2    / scalar
    cat_progn_div_scalar%tc4    =     cat_progn%tc4    / scalar
    cat_progn_div_scalar%qa1    =     cat_progn%qa1    / scalar
    cat_progn_div_scalar%qa2    =     cat_progn%qa2    / scalar
    cat_progn_div_scalar%qa4    =     cat_progn%qa4    / scalar
    cat_progn_div_scalar%capac  =     cat_progn%capac  / scalar
    cat_progn_div_scalar%catdef =     cat_progn%catdef / scalar
    cat_progn_div_scalar%rzexc  =     cat_progn%rzexc  / scalar
    cat_progn_div_scalar%srfexc =     cat_progn%srfexc / scalar
    
    do i=1,N_gt
       cat_progn_div_scalar%ght(i)  =     cat_progn%ght(i)  / scalar
    end do
    
    do i=1,N_snow
       cat_progn_div_scalar%wesn(i) =     cat_progn%wesn(i) / scalar
       cat_progn_div_scalar%htsn(i) =     cat_progn%htsn(i) / scalar
       cat_progn_div_scalar%sndz(i) =     cat_progn%sndz(i) / scalar
    end do
    
  end function cat_progn_div_scalar

  ! -----------------------------------------------------------

  function add_cat_progn( cat_progn_1, cat_progn_2 )
    
    implicit none

    type(cat_progn_type)             :: add_cat_progn
    type(cat_progn_type), intent(in) :: cat_progn_1, cat_progn_2

    integer :: i     ! local
    
    add_cat_progn%tc1    =     cat_progn_1%tc1    +     cat_progn_2%tc1
    add_cat_progn%tc2    =     cat_progn_1%tc2    +     cat_progn_2%tc2   
    add_cat_progn%tc4    =     cat_progn_1%tc4    +     cat_progn_2%tc4 
    add_cat_progn%qa1    =     cat_progn_1%qa1    +     cat_progn_2%qa1 
    add_cat_progn%qa2    =     cat_progn_1%qa2    +     cat_progn_2%qa2   
    add_cat_progn%qa4    =     cat_progn_1%qa4    +     cat_progn_2%qa4   
    add_cat_progn%capac  =     cat_progn_1%capac  +     cat_progn_2%capac   
    add_cat_progn%catdef =     cat_progn_1%catdef +     cat_progn_2%catdef 
    add_cat_progn%rzexc  =     cat_progn_1%rzexc  +     cat_progn_2%rzexc  
    add_cat_progn%srfexc =     cat_progn_1%srfexc +     cat_progn_2%srfexc 

    do i=1,N_gt
       add_cat_progn%ght(i)   = cat_progn_1%ght(i)  +   cat_progn_2%ght(i)
    end do
    
    do i=1,N_snow
       add_cat_progn%wesn(i)  = cat_progn_1%wesn(i)  +  cat_progn_2%wesn(i)  
       add_cat_progn%htsn(i)  = cat_progn_1%htsn(i)  +  cat_progn_2%htsn(i)  
       add_cat_progn%sndz(i)  = cat_progn_1%sndz(i)  +  cat_progn_2%sndz(i)  
    end do
    
    
  end function add_cat_progn
  
  ! ****************************************************
    
  subroutine scalar2cat_force( cat_force, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(cat_force_type), intent(out) :: cat_force
    
    cat_force%TRAINC       = scalar
    cat_force%TRAINL       = scalar
    cat_force%TSNOW        = scalar
    cat_force%UM           = scalar
    cat_force%ETURB1       = scalar
    cat_force%DEDQA1       = scalar
    cat_force%DEDTC1       = scalar
    cat_force%HSTURB1      = scalar
    cat_force%DHSDQA1      = scalar
    cat_force%DHSDTC1      = scalar
    cat_force%ETURB2       = scalar
    cat_force%DEDQA2       = scalar
    cat_force%DEDTC2       = scalar
    cat_force%HSTURB2      = scalar
    cat_force%DHSDQA2      = scalar
    cat_force%DHSDTC2      = scalar
    cat_force%ETURB4       = scalar
    cat_force%DEDQA4       = scalar
    cat_force%DEDTC4       = scalar
    cat_force%HSTURB4      = scalar
    cat_force%DHSDQA4      = scalar
    cat_force%DHSDTC4      = scalar
    cat_force%ETURBS       = scalar
    cat_force%DEDQAS       = scalar
    cat_force%DEDTCS       = scalar
    cat_force%HSTURBS      = scalar
    cat_force%DHSDQAS      = scalar
    cat_force%DHSDTCS      = scalar
    cat_force%TM           = scalar
    cat_force%QM           = scalar
    cat_force%ra1          = scalar
    cat_force%ra2          = scalar
    cat_force%ra4          = scalar
    cat_force%raS          = scalar
    cat_force%SUNANG       = scalar
    cat_force%PARDIR       = scalar
    cat_force%PARDIF       = scalar
    cat_force%SWNETF       = scalar
    cat_force%SWNETS       = scalar
    cat_force%HLWDWN       = scalar
    cat_force%PSUR         = scalar
    cat_force%ZLAI         = scalar
    cat_force%GREEN        = scalar
    cat_force%Z2           = scalar
    cat_force%SQSCAT       = scalar
    cat_force%RSOIL1       = scalar
    cat_force%RSOIL2       = scalar
    cat_force%RDC          = scalar
    cat_force%QSAT1        = scalar
    cat_force%DQS1         = scalar
    cat_force%ALW1         = scalar
    cat_force%BLW1         = scalar
    cat_force%QSAT2        = scalar
    cat_force%DQS2         = scalar
    cat_force%ALW2         = scalar
    cat_force%BLW2         = scalar
    cat_force%QSAT4        = scalar
    cat_force%DQS4         = scalar
    cat_force%ALW4         = scalar
    cat_force%BLW4         = scalar
    cat_force%QSATS        = scalar
    cat_force%DQSS         = scalar
    cat_force%ALWS         = scalar
    cat_force%BLWS         = scalar
    
  end subroutine scalar2cat_force
  
  ! -------------------------------------------------------------
  
  function cat_force_div_scalar( cat_force, scalar )
    
    implicit none
    
    type(cat_force_type)             :: cat_force_div_scalar
    type(cat_force_type), intent(in) :: cat_force
    
    real, intent(in) :: scalar
    
    cat_force_div_scalar%TRAINC  = cat_force%TRAINC       / scalar
    cat_force_div_scalar%TRAINL  = cat_force%TRAINL       / scalar
    cat_force_div_scalar%TSNOW   = cat_force%TSNOW        / scalar
    cat_force_div_scalar%UM      = cat_force%UM           / scalar
    cat_force_div_scalar%ETURB1  = cat_force%ETURB1       / scalar
    cat_force_div_scalar%DEDQA1  = cat_force%DEDQA1       / scalar
    cat_force_div_scalar%DEDTC1  = cat_force%DEDTC1       / scalar
    cat_force_div_scalar%HSTURB1 = cat_force%HSTURB1      / scalar
    cat_force_div_scalar%DHSDQA1 = cat_force%DHSDQA1      / scalar
    cat_force_div_scalar%DHSDTC1 = cat_force%DHSDTC1      / scalar
    cat_force_div_scalar%ETURB2  = cat_force%ETURB2       / scalar
    cat_force_div_scalar%DEDQA2  = cat_force%DEDQA2       / scalar
    cat_force_div_scalar%DEDTC2  = cat_force%DEDTC2       / scalar
    cat_force_div_scalar%HSTURB2 = cat_force%HSTURB2      / scalar
    cat_force_div_scalar%DHSDQA2 = cat_force%DHSDQA2      / scalar
    cat_force_div_scalar%DHSDTC2 = cat_force%DHSDTC2      / scalar
    cat_force_div_scalar%ETURB4  = cat_force%ETURB4       / scalar
    cat_force_div_scalar%DEDQA4  = cat_force%DEDQA4       / scalar
    cat_force_div_scalar%DEDTC4  = cat_force%DEDTC4       / scalar
    cat_force_div_scalar%HSTURB4 = cat_force%HSTURB4      / scalar
    cat_force_div_scalar%DHSDQA4 = cat_force%DHSDQA4      / scalar
    cat_force_div_scalar%DHSDTC4 = cat_force%DHSDTC4      / scalar
    cat_force_div_scalar%ETURBS  = cat_force%ETURBS       / scalar
    cat_force_div_scalar%DEDQAS  = cat_force%DEDQAS       / scalar
    cat_force_div_scalar%DEDTCS  = cat_force%DEDTCS       / scalar
    cat_force_div_scalar%HSTURBS = cat_force%HSTURBS      / scalar
    cat_force_div_scalar%DHSDQAS = cat_force%DHSDQAS      / scalar
    cat_force_div_scalar%DHSDTCS = cat_force%DHSDTCS      / scalar
    cat_force_div_scalar%TM      = cat_force%TM           / scalar
    cat_force_div_scalar%QM      = cat_force%QM           / scalar
    cat_force_div_scalar%ra1     = cat_force%ra1          / scalar
    cat_force_div_scalar%ra2     = cat_force%ra2          / scalar
    cat_force_div_scalar%ra4     = cat_force%ra4          / scalar
    cat_force_div_scalar%raS     = cat_force%raS          / scalar
    cat_force_div_scalar%SUNANG  = cat_force%SUNANG       / scalar
    cat_force_div_scalar%PARDIR  = cat_force%PARDIR       / scalar
    cat_force_div_scalar%PARDIF  = cat_force%PARDIF       / scalar
    cat_force_div_scalar%SWNETF  = cat_force%SWNETF       / scalar
    cat_force_div_scalar%SWNETS  = cat_force%SWNETS       / scalar
    cat_force_div_scalar%HLWDWN  = cat_force%HLWDWN       / scalar
    cat_force_div_scalar%PSUR    = cat_force%PSUR         / scalar
    cat_force_div_scalar%ZLAI    = cat_force%ZLAI         / scalar
    cat_force_div_scalar%GREEN   = cat_force%GREEN        / scalar
    cat_force_div_scalar%Z2      = cat_force%Z2           / scalar
    cat_force_div_scalar%SQSCAT  = cat_force%SQSCAT       / scalar
    cat_force_div_scalar%RSOIL1  = cat_force%RSOIL1       / scalar
    cat_force_div_scalar%RSOIL2  = cat_force%RSOIL2       / scalar
    cat_force_div_scalar%RDC     = cat_force%RDC          / scalar
    cat_force_div_scalar%QSAT1   = cat_force%QSAT1        / scalar
    cat_force_div_scalar%DQS1    = cat_force%DQS1         / scalar
    cat_force_div_scalar%ALW1    = cat_force%ALW1         / scalar
    cat_force_div_scalar%BLW1    = cat_force%BLW1         / scalar
    cat_force_div_scalar%QSAT2   = cat_force%QSAT2        / scalar
    cat_force_div_scalar%DQS2    = cat_force%DQS2         / scalar
    cat_force_div_scalar%ALW2    = cat_force%ALW2         / scalar
    cat_force_div_scalar%BLW2    = cat_force%BLW2         / scalar
    cat_force_div_scalar%QSAT4   = cat_force%QSAT4        / scalar
    cat_force_div_scalar%DQS4    = cat_force%DQS4         / scalar
    cat_force_div_scalar%ALW4    = cat_force%ALW4         / scalar
    cat_force_div_scalar%BLW4    = cat_force%BLW4         / scalar
    cat_force_div_scalar%QSATS   = cat_force%QSATS        / scalar
    cat_force_div_scalar%DQSS    = cat_force%DQSS         / scalar
    cat_force_div_scalar%ALWS    = cat_force%ALWS         / scalar
    cat_force_div_scalar%BLWS    = cat_force%BLWS         / scalar

  end function cat_force_div_scalar
  
   ! -----------------------------------------------------------

  function add_cat_force( cat_force_1, cat_force_2 )
    
    implicit none
    
    type(cat_force_type)             :: add_cat_force
    type(cat_force_type), intent(in) :: cat_force_1, cat_force_2
    
    add_cat_force%TRAINC  = cat_force_1%TRAINC      + cat_force_2%TRAINC   
    add_cat_force%TRAINL  = cat_force_1%TRAINL      + cat_force_2%TRAINL   
    add_cat_force%TSNOW   = cat_force_1%TSNOW       + cat_force_2%TSNOW    
    add_cat_force%UM      = cat_force_1%UM          + cat_force_2%UM       
    add_cat_force%ETURB1  = cat_force_1%ETURB1      + cat_force_2%ETURB1   
    add_cat_force%DEDQA1  = cat_force_1%DEDQA1      + cat_force_2%DEDQA1   
    add_cat_force%DEDTC1  = cat_force_1%DEDTC1      + cat_force_2%DEDTC1   
    add_cat_force%HSTURB1 = cat_force_1%HSTURB1     + cat_force_2%HSTURB1  
    add_cat_force%DHSDQA1 = cat_force_1%DHSDQA1     + cat_force_2%DHSDQA1  
    add_cat_force%DHSDTC1 = cat_force_1%DHSDTC1     + cat_force_2%DHSDTC1  
    add_cat_force%ETURB2  = cat_force_1%ETURB2      + cat_force_2%ETURB2   
    add_cat_force%DEDQA2  = cat_force_1%DEDQA2      + cat_force_2%DEDQA2   
    add_cat_force%DEDTC2  = cat_force_1%DEDTC2      + cat_force_2%DEDTC2   
    add_cat_force%HSTURB2 = cat_force_1%HSTURB2     + cat_force_2%HSTURB2  
    add_cat_force%DHSDQA2 = cat_force_1%DHSDQA2     + cat_force_2%DHSDQA2  
    add_cat_force%DHSDTC2 = cat_force_1%DHSDTC2     + cat_force_2%DHSDTC2  
    add_cat_force%ETURB4  = cat_force_1%ETURB4      + cat_force_2%ETURB4   
    add_cat_force%DEDQA4  = cat_force_1%DEDQA4      + cat_force_2%DEDQA4   
    add_cat_force%DEDTC4  = cat_force_1%DEDTC4      + cat_force_2%DEDTC4   
    add_cat_force%HSTURB4 = cat_force_1%HSTURB4     + cat_force_2%HSTURB4  
    add_cat_force%DHSDQA4 = cat_force_1%DHSDQA4     + cat_force_2%DHSDQA4  
    add_cat_force%DHSDTC4 = cat_force_1%DHSDTC4     + cat_force_2%DHSDTC4  
    add_cat_force%ETURBS  = cat_force_1%ETURBS      + cat_force_2%ETURBS   
    add_cat_force%DEDQAS  = cat_force_1%DEDQAS      + cat_force_2%DEDQAS   
    add_cat_force%DEDTCS  = cat_force_1%DEDTCS      + cat_force_2%DEDTCS   
    add_cat_force%HSTURBS = cat_force_1%HSTURBS     + cat_force_2%HSTURBS  
    add_cat_force%DHSDQAS = cat_force_1%DHSDQAS     + cat_force_2%DHSDQAS  
    add_cat_force%DHSDTCS = cat_force_1%DHSDTCS     + cat_force_2%DHSDTCS  
    add_cat_force%TM      = cat_force_1%TM          + cat_force_2%TM       
    add_cat_force%QM      = cat_force_1%QM          + cat_force_2%QM       
    add_cat_force%ra1     = cat_force_1%ra1         + cat_force_2%ra1      
    add_cat_force%ra2     = cat_force_1%ra2         + cat_force_2%ra2      
    add_cat_force%ra4     = cat_force_1%ra4         + cat_force_2%ra4      
    add_cat_force%raS     = cat_force_1%raS         + cat_force_2%raS      
    add_cat_force%SUNANG  = cat_force_1%SUNANG      + cat_force_2%SUNANG   
    add_cat_force%PARDIR  = cat_force_1%PARDIR      + cat_force_2%PARDIR   
    add_cat_force%PARDIF  = cat_force_1%PARDIF      + cat_force_2%PARDIF   
    add_cat_force%SWNETF  = cat_force_1%SWNETF      + cat_force_2%SWNETF   
    add_cat_force%SWNETS  = cat_force_1%SWNETS      + cat_force_2%SWNETS   
    add_cat_force%HLWDWN  = cat_force_1%HLWDWN      + cat_force_2%HLWDWN   
    add_cat_force%PSUR    = cat_force_1%PSUR        + cat_force_2%PSUR     
    add_cat_force%ZLAI    = cat_force_1%ZLAI        + cat_force_2%ZLAI     
    add_cat_force%GREEN   = cat_force_1%GREEN       + cat_force_2%GREEN    
    add_cat_force%Z2      = cat_force_1%Z2          + cat_force_2%Z2       
    add_cat_force%SQSCAT  = cat_force_1%SQSCAT      + cat_force_2%SQSCAT   
    add_cat_force%RSOIL1  = cat_force_1%RSOIL1      + cat_force_2%RSOIL1   
    add_cat_force%RSOIL2  = cat_force_1%RSOIL2      + cat_force_2%RSOIL2   
    add_cat_force%RDC     = cat_force_1%RDC         + cat_force_2%RDC      
    add_cat_force%QSAT1   = cat_force_1%QSAT1       + cat_force_2%QSAT1    
    add_cat_force%DQS1    = cat_force_1%DQS1        + cat_force_2%DQS1     
    add_cat_force%ALW1    = cat_force_1%ALW1        + cat_force_2%ALW1     
    add_cat_force%BLW1    = cat_force_1%BLW1        + cat_force_2%BLW1     
    add_cat_force%QSAT2   = cat_force_1%QSAT2       + cat_force_2%QSAT2    
    add_cat_force%DQS2    = cat_force_1%DQS2        + cat_force_2%DQS2     
    add_cat_force%ALW2    = cat_force_1%ALW2        + cat_force_2%ALW2     
    add_cat_force%BLW2    = cat_force_1%BLW2        + cat_force_2%BLW2     
    add_cat_force%QSAT4   = cat_force_1%QSAT4       + cat_force_2%QSAT4    
    add_cat_force%DQS4    = cat_force_1%DQS4        + cat_force_2%DQS4     
    add_cat_force%ALW4    = cat_force_1%ALW4        + cat_force_2%ALW4     
    add_cat_force%BLW4    = cat_force_1%BLW4        + cat_force_2%BLW4     
    add_cat_force%QSATS   = cat_force_1%QSATS       + cat_force_2%QSATS    
    add_cat_force%DQSS    = cat_force_1%DQSS        + cat_force_2%DQSS     
    add_cat_force%ALWS    = cat_force_1%ALWS        + cat_force_2%ALWS     
    add_cat_force%BLWS    = cat_force_1%BLWS        + cat_force_2%BLWS     

  end function add_cat_force
  
  ! ************************************************************
  
  subroutine scalar2cat_param( cat_param, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(cat_param_type), intent(out) :: cat_param
    
    integer :: i     ! local

    ! ---------------------

    cat_param%dpth       = scalar 
                         
    cat_param%dzsf       = scalar 
    cat_param%dzrz       = scalar 
    cat_param%dzpr       = scalar 
                         
    do i=1,N_gt          
       cat_param%dzgt(i) = scalar
    end do               
                         
    cat_param%poros      = scalar
    cat_param%cond       = scalar
    cat_param%psis       = scalar
    cat_param%bee        = scalar
                         
    cat_param%wpwet      = scalar
                         
    cat_param%gnu        = scalar 
        
    cat_param%vgwmax     = scalar
    
    cat_param%vegcls     = nint(scalar)
    cat_param%soilcls30  = nint(scalar)
    cat_param%soilcls100 = nint(scalar)
        
    cat_param%bf1        = scalar
    cat_param%bf2        = scalar
    cat_param%bf3        = scalar
    cat_param%cdcr1      = scalar
    cat_param%cdcr2      = scalar
    cat_param%ars1       = scalar
    cat_param%ars2       = scalar
    cat_param%ars3       = scalar
    cat_param%ara1       = scalar
    cat_param%ara2       = scalar
    cat_param%ara3       = scalar
    cat_param%ara4       = scalar
    cat_param%arw1       = scalar
    cat_param%arw2       = scalar
    cat_param%arw3       = scalar
    cat_param%arw4       = scalar
    cat_param%tsa1       = scalar
    cat_param%tsa2       = scalar
    cat_param%tsb1       = scalar
    cat_param%tsb2       = scalar
    cat_param%atau       = scalar
    cat_param%btau       = scalar
                         
    cat_param%gravel30   = scalar
    cat_param%orgC30     = scalar
    cat_param%orgC       = scalar
    cat_param%sand30     = scalar
    cat_param%clay30     = scalar
    cat_param%sand       = scalar
    cat_param%clay       = scalar
    cat_param%wpwet30    = scalar
    cat_param%poros30    = scalar

    cat_param%veghght    = scalar

  end subroutine scalar2cat_param

  ! **************************************************************************
  !
  ! utilities to convert from "cat_progn" type to regular arrays
  !
  ! NOTE: the functions catprogn2xxx can only be used within the argument list 
  !       of a subroutine when the prognostic variables are "intent(in)"
  !
  ! **************************************************************************
  
  function catprogn2wesn(N_cat, cat_progn)
    
    implicit none
    
    integer,                                intent(in) :: N_cat
    type(cat_progn_type), dimension(N_cat), intent(in) :: cat_progn
    
    real, dimension(N_snow,N_cat) :: catprogn2wesn
    
    ! local variables
    
    integer :: i
    
    ! --------------------------------

    do i=1,N_snow
       
       catprogn2wesn(i,:) = cat_progn(:)%wesn(i)
       
    end do

  end function catprogn2wesn
  
  ! ***********************************************************************

  function catprogn2htsn(N_cat, cat_progn)
    
    implicit none
    
    integer,                                intent(in) :: N_cat
    type(cat_progn_type), dimension(N_cat), intent(in) :: cat_progn
    
    real, dimension(N_snow,N_cat) :: catprogn2htsn
    
    ! local variables
    
    integer :: i
    
    ! --------------------------------
    
    do i=1,N_snow
       
       catprogn2htsn(i,:) = cat_progn(:)%htsn(i)
       
    end do

  end function catprogn2htsn

  ! ***********************************************************************

  function catprogn2sndz(N_cat, cat_progn)
    
    implicit none
    
    integer,                                intent(in) :: N_cat
    type(cat_progn_type), dimension(N_cat), intent(in) :: cat_progn
    
    real, dimension(N_snow,N_cat) :: catprogn2sndz
    
    ! local variables
    
    integer :: i
    
    ! --------------------------------

    do i=1,N_snow
       
       catprogn2sndz(i,:) = cat_progn(:)%sndz(i)
       
    end do

  end function catprogn2sndz
  
  ! ***********************************************************************

  function catprogn2ghtcnt(N_cat, cat_progn)
    
    implicit none
    
    integer,                                intent(in) :: N_cat
    type(cat_progn_type), dimension(N_cat), intent(in) :: cat_progn
    
    real, dimension(N_gt,N_cat) :: catprogn2ghtcnt
    
    ! local variables
    
    integer :: i
    
    ! --------------------------------

    do i=1,N_gt
       
       catprogn2ghtcnt(i,:) = cat_progn(:)%ght(i)
       
    end do

  end function catprogn2ghtcnt
  
  ! ***********************************************************************
    
end module catch_types
  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if 0

! driver routines for testing

program test_catch_types

  ! use module catch_types
  
  implicit none
  
  type(cat_diagS_type) :: cat_diagS_1, cat_diagS_2
  
  cat_diagS_1 = 1.
  cat_diagS_2 = 2.
  
  write (*,*) cat_diagS_1
  write (*,*) cat_diagS_2
  
  cat_diagS_2 = cat_diagS_1 + cat_diagS_2

  write (*,*) cat_diagS_1
  write (*,*) cat_diagS_2

end program test_catch_types

#endif

! ========================== EOF ==================================
