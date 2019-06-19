
module LDAS_EASE_conv
  
  ! ==========================================================================
  !
  ! easeV1_conv.F90 - FORTRAN routines for conversion of azimuthal 
  !                   equal area and equal area cylindrical grid coordinates
  ! 
  ! 30-Jan-1992 H.Maybee
  ! 20-Mar-1992 Ken Knowles  303-492-0644  knowles@kryos.colorado.edu
  ! 16-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
  !              Copied from nsmconv.f, changed resolutions from 
  !              40-20-10 km to 25-12.5 km
  ! 21-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
  !              Fixed sign of Southern latitudes in ease_inverse.
  ! 12-Sep-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
  ! 	       Changed grid cell size. Changed "c","f" to "l","h"
  ! 25-Oct-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
  ! 	       Changed row size from 587 to 586 for Mercator projection
  ! 11-May-2011 reichle: Changed "smap" to "easeV1".  
  !                      Added SSM/I and AMSR-E "M25" grid.
  !                      So far ONLY for cylindrical grids.
  !                      Converted from *.f to *.F90 module
  ! 
  ! $Log$
  ! Revision 1.1.2.3  2018/09/13 20:42:50  wjiang
  ! change M25
  !
  ! Revision 1.1.2.2  2017/09/18 15:10:25  wjiang
  ! "fix" strange compiler erro on comment line.
  !
  ! Revision 1.1.2.1  2017/01/19 19:35:58  wjiang
  ! 1)add EASE grid support
  ! 2)add ensemble average for HISTORY
  !
  ! Revision 1.2  2014-08-26 17:33:55  rreichle
  ! - clean-up of *.F90 in GEOSlana_GridComp:
  ! - make sure all modules include "private" statement
  ! - follow all "use" statements with "ONLY"
  ! - removed unused files (esat_qsat.F90, nr_sort.f)
  ! - removed unused variables
  !
  ! Revision 1.1  2011-05-11 21:58:46  rreichle
  !
  ! Adding utilities to map between EASE grids and lat/lon coordinates.
  !
  ! Revision 1.3  1994/11/01 23:40:43  brodzik
  ! Replaced all references to 'ease' with 'smap'
  ! Replaced all references to 'smap' with 'easeV1' -- reichle
  ! 
  ! ==========================================================================

  implicit none
  
  private

  public :: ease_convert
  public :: ease_inverse
  public :: easeV1_convert
  public :: easeV1_inverse
  public :: easeV2_convert
  public :: easeV2_inverse
  public :: easeV2_extent


  ! ***NEVER*** change these constants to GEOS-5 MAPL constants!!!!
    
  ! radius of the earth (km), authalic sphere based on International datum 
  
  real*8, parameter :: RE_km = 6371.228
  
  ! scale factor for standard paralles at +/-30.00 degrees
  
  real*8, parameter :: COS_PHI1 = .866025403
  
  real*8, parameter :: PI  = 3.14159265358979323846
 
  ! ==========================================================================
  !
  ! easeV2_conv.F90 - FORTRAN routines for converting grid coordinates
  !                   (latitude/longitude <--> row/column indices)
  !                   of the Equal Area Scalable Earth, version 2 (EASEv2) grid
  !
  !    ***** ONLY cylindrical ('M') projection implemented *****
  !
  ! Ported from Steven Chan's matlab code (smapease2inverse.m, 
  ! smapease2forward.m), which has been ported from NSIDC's IDL code
  ! (wgs84_convert.pro, wgs84_inverse.pro) available from  
  ! ftp://sidads.colorado.edu/pub/tools/easegrid/geolocation_tools/
  !
  ! 04-Apr-2013 - reichle
  !
  ! Official references:
  !  doi:10.3390/ijgi1010032
  !  doi:10.3390/ijgi3031154 -- correction of M25 "map_scale_m" parameters!
  !
  ! 04-Apr-2013 - reichle
  ! 11-Sep-2018 - reichle, mgirotto -- added 'M25' grid parameters 
  !
  ! ==========================================================================


  ! ***NEVER*** change these constants to GEOS-5 MAPL constants!!!!
  
  ! radius of the earth (m) and map eccentricity
  
  real*8, parameter :: map_equatorial_radius_m         = 6378137.0 
  
  real*8, parameter :: map_eccentricity                = 0.081819190843
  
  
  real*8, parameter :: e2      = map_eccentricity * map_eccentricity
  real*8, parameter :: e4      = e2 * e2
  real*8, parameter :: e6      = e2 * e4
  
  real*8, parameter :: epsilon = 1.e-6
  
  real*8, parameter :: map_reference_longitude         =   0.0  ! 'M', 'N', 'S'
  
  ! constants for 'N' and 'S' (azimuthal) projections
  
  real*8, parameter :: N_map_reference_latitude        =  90.0  
  real*8, parameter :: S_map_reference_latitude        = -90.0
  
  ! constants for 'M' (cylindrical) projection
  
  real*8, parameter :: M_map_reference_latitude        =   0.0
  real*8, parameter :: M_map_second_reference_latitude =  30.0
  
  real*8, parameter :: M_sin_phi1 = sin(M_map_second_reference_latitude*PI/180.)
  real*8, parameter :: M_cos_phi1 = cos(M_map_second_reference_latitude*PI/180.)
  
  real*8, parameter :: M_kz = M_cos_phi1/sqrt(1.0-e2*M_sin_phi1*M_sin_phi1)
  
  
contains  
  
  subroutine ease_convert (gridname, lat, lon, r, s)
    character*(*), intent(in)  :: gridname
    real,          intent(in)  :: lat, lon
    real,          intent(out) :: r, s
    character(3)  :: grid

    if (index(gridname,'M36') /=0 ) then
        grid='M36'
    else if (index(gridname,'M25') /=0 ) then
        grid='M25'
    else if (index(gridname,'M09') /=0 ) then
        grid='M09'
    else if (index(gridname,'M03') /=0 ) then
        grid='M03'
    else if (index(gridname,'M01') /=0 ) then
        grid='M01'
    endif
    
    if(index(gridname,'EASEv2') /=0) then
        call easeV2_convert(grid,lat,lon,r,s)
    else if(index(gridname,'EASE') /=0) then
        call easeV1_convert(grid,lat,lon,r,s)
    else
       print*,"wrong gridname: "//gridname
    endif 
  end subroutine

  subroutine ease_inverse (gridname, r, s, lat, lon)
    character*(*), intent(in)  :: gridname
    real,          intent(in)  :: r, s
    real,          intent(out) :: lat, lon
    character(3)  :: grid

    if (index(gridname,'M36') /=0 ) then
        grid='M36'
    else if (index(gridname,'M25') /=0 ) then
        grid='M25'
    else if (index(gridname,'M09') /=0 ) then
        grid='M09'
    else if (index(gridname,'M03') /=0 ) then
        grid='M03'
    else if (index(gridname,'M01') /=0 ) then
        grid='M01'
    endif
    
    if(index(gridname,'EASEv2') /=0) then
        call easeV2_inverse(grid,r,s,lat,lon)
    else if(index(gridname,'EASE') /=0) then
        call easeV1_inverse(grid,r,s,lat,lon)
    else
       print*,"wrong gridname: "//gridname
    endif 
  end subroutine ease_inverse

  ! *******************************************************************
  
  subroutine easeV1_convert (grid, lat, lon, r, s)
    
    ! convert geographic coordinates (spherical earth) to 
    ! azimuthal equal area or equal area cylindrical grid coordinates
    ! 
    ! status = easeV1_convert (grid, lat, lon, r, s)
    ! 
    ! input : grid - projection name '[M][xx]'
    !            where xx = approximate resolution [km]
    !               ie xx = "01", "03", "09", "36"       (SMAP)
    !               or xx = "12", "25"                   (SSM/I, AMSR-E)
    ! 	    lat, lon = geo. coords. (decimal degrees)
    ! 
    ! output: r, s - column, row coordinates
    ! 
    ! result: status = 0 indicates normal successful completion
    ! 		-1 indicates error status (point not on grid)
    ! 
    ! --------------------------------------------------------------------------
        
    character*(*), intent(in)  :: grid
    real,          intent(in)  :: lat, lon
    real,          intent(out) :: r, s

    ! local variables
    
    integer :: cols, rows
    real*8  :: Rg, phi, lam, rho, CELL_km, r0, s0
    
    ! ---------------------------------------------------------------------
    
    call easeV1_get_params( grid, CELL_km, cols, rows, r0, s0, Rg )
    
    phi = lat*PI/180.   ! convert from degree to radians
    lam = lon*PI/180.   ! convert from degree to radians
    
    if (grid(1:1).eq.'N') then
       rho = 2 * Rg * sin(PI/4. - phi/2.)
       r = r0 + rho * sin(lam)
       s = s0 + rho * cos(lam)
       
    else if (grid(1:1).eq.'S') then
       rho = 2 * Rg * cos(PI/4. - phi/2.)
       r = r0 + rho * sin(lam)
       s = s0 - rho * cos(lam)
       
    else if (grid(1:1).eq.'M') then
       r = r0 + Rg * lam * COS_PHI1
       s = s0 - Rg * sin(phi) / COS_PHI1
       
    endif
        
  end subroutine easeV1_convert
  
  ! *******************************************************************
  
  subroutine easeV1_inverse (grid, r, s, lat, lon)
    
    ! convert azimuthal equal area or equal area cylindrical 
    ! grid coordinates to geographic coordinates (spherical earth)
    ! 
    ! status = easeV1_inverse (grid, r, s, lat, lon)
    ! 
    ! input : grid - projection name '[M][xx]'
    !            where xx = approximate resolution [km]
    !               ie xx = "01", "03", "09", "36"       (SMAP)
    !               or xx = "12", "25"                   (SSM/I, AMSR-E)
    ! 	    r, s - column, row coordinates
    ! 
    ! output: lat, lon = geo. coords. (decimal degrees)
    ! 
    ! result: status = 0 indicates normal successful completion
    ! 		-1 indicates error status (point not on grid)
    ! 
    ! --------------------------------------------------------------------------

    character*(*), intent(in)  :: grid
    real,          intent(in)  :: r, s
    real,          intent(out) :: lat, lon

    ! local variables
    
    integer :: cols, rows
    real*8    :: Rg, phi, lam, rho, CELL_km, r0, s0
    real*8    :: gamma, beta, epsilon, x, y, c
    real*8    :: sinphi1, cosphi1

    ! ---------------------------------------------------------------------
    
    call easeV1_get_params( grid, CELL_km, cols, rows, r0, s0, Rg )
        
    x = r - r0
    y = -(s - s0)
    
    if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then 
       rho = sqrt(x*x + y*y)
       if (rho.eq.0.0) then
          if (grid(1:1).eq.'N') lat = 90.0 
          if (grid(1:1).eq.'S') lat = -90.0 
          lon = 0.0
       else
          if (grid(1:1).eq.'N') then
             sinphi1 = sin(PI/2.)
             cosphi1 = cos(PI/2.)
             if (y.eq.0.) then
                if (r.le.r0) lam = -PI/2.
                if (r.gt.r0) lam = PI/2.
             else
                lam = atan2(x,-y)
             endif
          else if (grid(1:1).eq.'S') then
             sinphi1 = sin(-PI/2.)
             cosphi1 = cos(-PI/2.)
             if (y.eq.0.) then
                if (r.le.r0) lam = -PI/2.
                if (r.gt.r0) lam = PI/2.
             else
                lam = atan2(x,y)
             endif
          endif
          gamma = rho/(2 * Rg)
          if (abs(gamma) .gt. 1.) return
          c = 2 * asin(gamma)
          beta = cos(c) * sinphi1 + y * sin(c) * (cosphi1/rho)
          if (abs(beta).gt.1.) return
          phi = asin(beta)
          lat = phi*180./PI   ! convert from radians to degree
          lon = lam*180./PI   ! convert from radians to degree
       endif
       
    else if (grid(1:1).eq.'M') then
       
       ! 	  allow .5 cell tolerance in arcsin function
       ! 	  so that grid coordinates which are less than .5 cells
       ! 	  above 90.00N or below 90.00S are given a lat of 90.00
       
       epsilon = 1 + 0.5/Rg
       beta = y*COS_PHI1/Rg
       if (abs(beta).gt.epsilon) return
       if (beta.le.-1.) then
          phi = -PI/2.
       else if (beta.ge.1.) then
          phi = PI/2.
       else
          phi = asin(beta)
       endif
       lam = x/COS_PHI1/Rg
       lat = phi*180./PI   ! convert from radians to degree
       lon = lam*180./PI   ! convert from radians to degree
    endif
        
  end subroutine easeV1_inverse
  
  ! *******************************************************************

  subroutine easeV1_get_params( grid, CELL_km, cols, rows, r0, s0, Rg )
    
    implicit none
    
    character*(*), intent(in)  :: grid
    real*8,        intent(out) :: CELL_km, r0, s0, Rg
    integer,       intent(out) :: cols, rows
    
    ! --------------------------------------------------------
    !
    ! r0,s0 are defined such that cells at all scales have 
    ! coincident center points
    ! 
    !c        r0 = (cols-1)/2. * scale
    !c        s0 = (rows-1)/2. * scale
    !
    ! --------------------------------------------------------
    
    if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then

       print *,'Polar projections not implemented yet'
       stop

    else if (grid(1:1).eq.'M') then

       if      (grid .eq. 'M36') then ! SMAP 36 km grid
          CELL_km = 36.00040279063    ! nominal cell size in kilometers
          cols = 963
          rows = 408
          r0 = 481.0
          s0 = 203.5
       
       else if (grid .eq. 'M25') then ! SSM/I, AMSR-E 25 km grid
          CELL_km = 25.067525         ! nominal cell size in kilometers
          cols = 1383
          rows = 586
          r0 = 691.0
          s0 = 292.5
       
       else if (grid .eq. 'M09') then ! SMAP  9 km grid
          CELL_km = 9.00010069766     ! nominal cell size in kilometers
          cols = 3852
          rows = 1632
          r0 = 1925.5
          s0 = 815.5
       
       else if (grid .eq. 'M03') then ! SMAP  3 km grid
          CELL_km = 3.00003356589     ! nominal cell size in kilometers
          cols = 11556
          rows = 4896
          r0 = 5777.5
          s0 = 2447.5
       
       else if (grid .eq. 'M01') then ! SMAP  1 km grid
          CELL_km = 1.00001118863     ! nominal cell size in kilometers
          cols = 34668
          rows = 14688
          r0 = 17333.5
          s0 = 7343.5
       
       else
       
          print *,'easeV1_convert: unknown resolution: ',grid
          stop
       
       endif

    else
       
       print *, 'easeV1_convert: unknown projection: ', grid
       stop
       
    endif
        
    Rg = RE_km/CELL_km
    
  end subroutine easeV1_get_params
  

  
  subroutine easeV2_convert (grid, lat, lon, col_ind, row_ind)
    
    ! convert geographic coordinates (spherical earth) to 
    ! azimuthal equal area or equal area cylindrical grid coordinates
    ! 
    ! *** NOTE order of calling arguments:  "lat-lon-lon-lat" ***
    !
    ! useage: call easeV2_convert (grid, lat, lon, r, s)
    ! 
    ! input : grid - projection name '[M][xx]'
    !            where xx = approximate resolution [km]
    !               ie xx = "01", "03", "09", "36"       (SMAP)
    ! 	      lat, lon = geo. coords. (decimal degrees)
    ! 
    ! output: col_ind, row_ind - column, row coordinates
    ! 
    ! --------------------------------------------------------------------------
        
    character*(*), intent(in)  :: grid
    real,          intent(in)  :: lat, lon
    real,          intent(out) :: col_ind, row_ind

    ! local variables
    
    integer :: cols, rows
    real*8  :: dlon, phi, lam, map_scale_m, r0, s0, ms, x, y, sin_phi, q
    
    ! ---------------------------------------------------------------------
    
    call easeV2_get_params( grid, map_scale_m, cols, rows, r0, s0 )

    dlon = lon
    
    if (abs(map_reference_longitude)>epsilon) then
       
       dlon = lon - map_reference_longitude
       
    end if

    if (dlon .lt. -180.0) dlon = dlon + 360.0
    if (dlon .gt.  180.0) dlon = dlon - 360.0
    
    phi =  lat*PI/180.   ! convert from degree to radians
    lam = dlon*PI/180.   ! convert from degree to radians
    
    sin_phi = sin(phi)
    
    ms      = map_eccentricity*sin_phi
    
    q = (1. - e2)*                                                     &
         (                                                             &
         (sin_phi /(1. - e2*sin_phi*sin_phi))                          &
         -                                                             &
         .5/map_eccentricity*log((1.-ms)/(1.+ms))                      &
         )
    
    ! note: "qp" only needed for 'N' and 'S' projections
    
    if      (grid(1:1).eq.'M') then
       
       x =  map_equatorial_radius_m*M_kz*lam
       
       y = (map_equatorial_radius_m*q)/(2.*M_kz)
       
    else
       
       print *,'Polar projections not implemented yet'
       stop
       
    endif
    
    row_ind = s0 - (y/map_scale_m)
    col_ind = r0 + (x/map_scale_m)
    
  end subroutine easeV2_convert
  
  ! *******************************************************************
  
  subroutine easeV2_inverse (grid, r, s, lat, lon)
    
    ! convert azimuthal equal area or equal area cylindrical 
    ! grid coordinates to geographic coordinates (spherical earth)
    ! 
    ! *** NOTE order of calling arguments:  "lon-lat-lat-lon" ***
    !
    ! useage: call easeV1_inverse (grid, r, s, lat, lon)
    ! 
    ! input : grid - projection name '[M][xx]'
    !            where xx = approximate resolution [km]
    !               ie xx = "01", "03", "09", "36"       (SMAP)
    ! 	      r, s - column, row coordinates
    ! 
    ! output: lat, lon = geo. coords. (decimal degrees)
    ! 
    ! --------------------------------------------------------------------------

    character*(*), intent(in)  :: grid
    real,          intent(in)  :: r, s
    real,          intent(out) :: lat, lon

    ! local variables
    
    integer   :: cols, rows
    real*8    :: phi, lam, map_scale_m, r0, s0, beta, x, y, qp
    
    ! ---------------------------------------------------------------------
    
    call easeV2_get_params( grid, map_scale_m, cols, rows, r0, s0 )
    
    x =  (r - r0)*map_scale_m
    y = -(s - s0)*map_scale_m
    
    qp = (1. - e2)*                                                           &
         (                                                                    &
         (1./(1.-e2))                                                         &
         -                                                                    &
         .5/map_eccentricity*log((1.-map_eccentricity)/(1.+map_eccentricity)) &
         )
    
    if      (grid(1:1).eq.'M') then
       
       beta = asin(2.*y*M_kz/(map_equatorial_radius_m*qp))
       
       lam  = x/(map_equatorial_radius_m*M_kz)
       
    else
       
       print *,'Polar projections not implemented yet'
       stop
       
    endif
    
    phi = beta                                                              &
         + ( ( e2/3.       + 31./180.*e4 + 517./ 5040.*e6 )*sin(2.*beta) )  &
         + ( (               23./360.*e4 + 251./ 3780.*e6 )*sin(4.*beta) )  &
         + ( (                             761./45360.*e6 )*sin(6.*beta) )
    
    lat = phi*180./PI                            ! convert from radians to degree
    lon = lam*180./PI + map_reference_longitude  ! convert from radians to degree
    
    if (lon .lt. -180.0) lon = lon + 360.0
    if (lon .gt.  180.0) lon = lon - 360.0
    
  end subroutine easeV2_inverse
  
  ! *******************************************************************
  
  subroutine easeV2_get_params( grid, map_scale_m, cols, rows, r0, s0 )
    
    implicit none
    
    character*(*), intent(in)  :: grid
    real*8,        intent(out) :: map_scale_m, r0, s0
    integer,       intent(out) :: cols, rows
    
    
    if (grid(1:1).eq.'M') then
       
       if      (grid .eq. 'M36') then      ! SMAP 36 km grid
          
          map_scale_m = 36032.220840584    ! nominal cell size in meters
          cols = 964
          rows = 406
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2.0


       else if (grid .eq. 'M25') then      ! 25 km grid  

          map_scale_m = 25025.2600000      ! nominal cell size in meters (see doi:10.3390/ijgi3031154)
          cols = 1388
          rows =  584
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2

       else if (grid .eq. 'M09') then      ! SMAP  9 km grid

          map_scale_m = 9008.055210146     ! nominal cell size in meters
          cols = 3856
          rows = 1624
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2.0
          
       else if (grid .eq. 'M03') then      ! SMAP  3 km grid

          map_scale_m = 3002.6850700487    ! nominal cell size in meters
          cols = 11568
          rows = 4872
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2.0
          
       else if (grid .eq. 'M01') then      ! SMAP  1 km grid

          map_scale_m = 1000.89502334956   ! nominal cell size in meters
          cols = 34704
          rows = 14616
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2.0
       
       else
          
          print *,'easeV2_convert: unknown resolution: ',grid
          stop
       
       endif

    else if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then
       
       print *,'Polar projections not implemented yet'
       stop
       
    else
       
       print *, 'easeV2_convert: unknown projection: ', grid
       stop
       
    endif
           
  end subroutine easeV2_get_params
  
  ! *******************************************************************
  
  subroutine easeV2_extent( grid, N_cols, N_rows )

    ! simple wrapper to get N_cols (N_lon) and N_rows (N_lat)
    
    implicit none
    
    character*(*), intent(in)  :: grid
    integer,       intent(out) :: N_cols, N_rows
    
    ! local variables
    
    real*8                     :: map_scale_m, r0, s0
    
    ! ------------------------------------------------
    
    call easeV2_get_params( grid, map_scale_m, N_cols, N_rows, r0, s0 )
    
  end subroutine easeV2_extent

  ! *******************************************************************

end module LDAS_EASE_conv

! =============================== EOF =================================

