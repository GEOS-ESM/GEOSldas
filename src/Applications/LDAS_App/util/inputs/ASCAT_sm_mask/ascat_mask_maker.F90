! This program produces a combined mask for use with the assimilation of ASCAT soil moisture retrievals in GEOSldas. 
! The combined mask is based on component masks from:
!
! Lindorfer, R., Wagner, W., Hahn, S., Kim, H., Vreugdenhil, M., Gruber, A., Fischer, M., & Trnka, M. (2023). 
! Global Scale Maps of Subsurface Scattering Signals Impacting ASCAT Soil Moisture Retrievals (1.0.0) [Data set]. 
! TU Wien. https://doi.org/10.48436/9a2y9-e5z14
!
! The program provides the possibility to combine different masks (default is combination of subsurface and wetland masks)
! and interpolates onto a regular grid with a (hardwired) 0.1 degree lat/lon spacing and -90/-180 degree lower left 
! corner used for quick indexing in the ASCAT observation reader QC routine.
!
! Author: AM Fox, March, 2024

program ascat_mask_maker

    use netcdf

    implicit none

    integer                            :: ncid, varid, dimid, ierr, len, N_gpi, dimids(2), N_lon, N_lat
    integer                            :: i, j, closest_index, mask_mode
    integer, dimension(:), allocatable :: cold_mask, wet_mask, veg_mask, subsurface_mask, combined_mask

    integer(kind=1), dimension(:,:), allocatable :: mask_out
    integer(kind=1)                              :: missing_value

    real, dimension(:), allocatable :: asc_lon, asc_lat
    real, dimension(:), allocatable :: lon, lat, distances
    real                            :: d_lon, d_lat, ll_lon, ll_lat

    character(200)                  :: fname_in, mask_description, fname_out

    ! --------------------------------------------------------------------------------
    !
    ! hardwired variables

    ! ASCAT soil moisture mask file from Lindorfer et al 2023 
    
    fname_in = '/discover/nobackup/projects/gmao/smap/SMAP_Nature/ASCAT_EUMETSAT/Mask/subsurface_scattering_ASCAT_ERA5_Land.nc'

    ! Specification of how to combine the masks
    ! Mask_mode = 1 (default) combines subsurface and wetland masks
    ! Mask_mode = 2 uses only the subsurface mask
    ! Mask_mode = 3 uses only the wetland mask
    ! Mask_mode = 4 combines subsurface, wetland and vegetation masks

    mask_mode = 1
    
    ! Specification of output grid and missing value 
    
    d_lon  =    0.1
    d_lat  =    0.1
    ll_lon = -180.0   ! longitude of boundary of lower-left grid cell (longitude of lower left corner of grid)
    ll_lat =  -90.0   ! latitude  of boundary of lower-left grid cell (latitude  of lower left corner of grid)
    
    missing_value = -128

    ! Specify the NetCDF file name for the output mask

    fname_out = 'ascat_combined_mask_p1.nc'
    
    ! -------------------------------
    
    ! Open the NetCDF input file
    ierr = nf90_open(fname_in, nf90_nowrite, ncid)
    if (ierr /= nf90_noerr) stop 'Error opening file: ' // trim(fname_in)
    
    ! Data in original mask file are on the 12.5 km fixed Earth grid used for ASCAT (WARP5 grid) and
    !   stored in the NetCDF file as 1-dimensional arrays of length N_gpi (over land only).

    ! Get the dimension ID
    ierr = nf90_inq_dimid(ncid, 'gpi', dimid)
    if (ierr /= nf90_noerr) stop 'Error getting dimension ID'
    
    ! Get the length of the dimension
    ierr = nf90_inquire_dimension(ncid, dimid, len = N_gpi)
    if (ierr /= nf90_noerr) stop 'Error inquiring dimension'
    
    print*, 'N_gpi = ', N_gpi

    ! Allocate the arrays
    allocate(asc_lon(        N_gpi))
    allocate(asc_lat(        N_gpi))
    allocate(cold_mask(      N_gpi))
    allocate(wet_mask(       N_gpi))
    allocate(veg_mask(       N_gpi))
    allocate(subsurface_mask(N_gpi))
    allocate(combined_mask(  N_gpi))
    allocate(distances(      N_gpi))

    ! Get the variable IDs and read the variables
    ierr = nf90_inq_varid(ncid, 'lon', varid)
    ierr = nf90_get_var(ncid, varid, asc_lon)
    ierr = nf90_inq_varid(ncid, 'lat', varid)
    ierr = nf90_get_var(ncid, varid, asc_lat)
    ierr = nf90_inq_varid(ncid, 'cold_mask', varid)
    ierr = nf90_get_var(ncid, varid, cold_mask)
    ierr = nf90_inq_varid(ncid, 'wet_mask', varid)
    ierr = nf90_get_var(ncid, varid, wet_mask)
    ierr = nf90_inq_varid(ncid, 'veg_mask', varid)
    ierr = nf90_get_var(ncid, varid, veg_mask)
    ierr = nf90_inq_varid(ncid, 'subsurface_mask', varid)
    ierr = nf90_get_var(ncid, varid, subsurface_mask)

    ! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= nf90_noerr) stop 'Error closing file'

    ! Combine the masks (1-dim arrays)
    select case (mask_mode)
    case (1)
        ! Combine wet_mask and subsurface_mask
        mask_description = 'Combined subsurface and wetland mask'
        where (wet_mask /= 0)
            combined_mask = 1
        elsewhere
            combined_mask = subsurface_mask
        end where
    case (2)
        ! Use only subsurface_mask
        mask_description = 'Used only subsurface mask'
        combined_mask = subsurface_mask
    case (3)
        ! Use only wet_mask
        mask_description = 'Used only wetland mask'
        combined_mask = wet_mask
    case (4)
        ! Combine subsurface_mask, wet_mask, and veg_mask
        mask_description = 'Combined subsurface, wetland, and vegetation mask'
        where (wet_mask /= 0 .or. veg_mask /= 0)
            combined_mask = 1
        elsewhere
            combined_mask = subsurface_mask
        end where
    end select

    ! Re-map "combined_mask" from WARP5 input grid to regular lat/lon output grid (2-dim array)

    N_lon = nint( 360.0 / d_lon ) 
    N_lat = nint( 180.0 / d_lat ) 

    allocate(lon(N_lon))
    allocate(lat(N_lat))

    lon = [((ll_lon + (d_lon / 2)) + i * d_lon, i = 0, N_lon - 1)]  ! NB using grid cell centers for nearest neighbor search
    lat = [((ll_lat + (d_lat / 2)) + i * d_lat, i = 0, N_lat - 1)]

    allocate(mask_out( N_lon, N_lat))

    do i = 1, N_lon
        print*, lon(i)
        do j = 1, N_lat
            distances     = (asc_lon - lon(i))**2 + (asc_lat - lat(j))**2
            closest_index = minloc(distances, dim = 1)
            if (distances(closest_index) > 0.14**2) then
                mask_out(i, j) = missing_value                      ! Note: ASCAT EUMETSAT reader masks everything .ne. 0
            else
                mask_out(i, j) = combined_mask(closest_index)
            end if
        end do
    end do

    ! Write out the mask to netcdf
    ierr = nf90_create(fname_out, nf90_clobber, ncid)
    if (ierr /= nf90_noerr) stop 'Error creating file: ' // fname_out

    ! Define the dimensions
    ierr = nf90_def_dim(ncid, 'lon', N_lon, dimids(1))
    ierr = nf90_def_dim(ncid, 'lat', N_lat, dimids(2))

    ! Define the global attributes
    ierr = nf90_put_att(ncid, nf90_global, 'title', 'ASCAT combined mask')
    ierr = nf90_put_att(ncid, nf90_global, 'source', 'Lindorfer et al 2023 doi:10.48436/9a2y9-e5z14')
    ierr = nf90_put_att(ncid, nf90_global, 'description', mask_description)

    ! Define the variables
    ierr = nf90_def_var(ncid, 'lat', nf90_real, dimids(2), varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'latitude')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'grid cell center latitude')
    ierr = nf90_put_att(ncid, varid, 'units', 'degrees_north')
    ierr = nf90_put_att(ncid, varid, 'axis', 'Y')

    ierr = nf90_def_var(ncid, 'lon', nf90_real, dimids(1), varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'longitude')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'grid cell center longitude')
    ierr = nf90_put_att(ncid, varid, 'units', 'degrees_east')
    ierr = nf90_put_att(ncid, varid, 'axis', 'X')

    ierr = nf90_def_var(ncid, 'mask', nf90_byte, dimids, varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'combined_mask')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'combined mask')
    ierr = nf90_put_att(ncid, varid, 'units', 'boolean')
    ierr = nf90_put_att(ncid, varid, '_FillValue', missing_value)

    ierr = nf90_def_var(ncid, 'll_lon', nf90_real, varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'longitude of lower left corner')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'longitude of lower left corner')
    ierr = nf90_put_att(ncid, varid, 'units', 'degrees_east')
    ierr = nf90_put_att(ncid, varid, 'axis', 'X')

    ierr = nf90_def_var(ncid, 'll_lat', nf90_real, varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'latitude of lower left corner')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'latitude of lower left corner')
    ierr = nf90_put_att(ncid, varid, 'units', 'degrees_north')
    ierr = nf90_put_att(ncid, varid, 'axis', 'Y')

    ierr = nf90_def_var(ncid, 'd_lon', nf90_real, varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'longitude grid spacing')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'longitude grid spacing')
    ierr = nf90_put_att(ncid, varid, 'units', 'degrees')
    ierr = nf90_put_att(ncid, varid, 'axis', 'X')

    ierr = nf90_def_var(ncid, 'd_lat', nf90_real, varid)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'latitude grid spacing')
    ierr = nf90_put_att(ncid, varid, 'units', 'degrees')
    ierr = nf90_put_att(ncid, varid, 'axis', 'Y')

    ! End define mode
    ierr = nf90_enddef(ncid)

    ! Write the variables
    ierr = nf90_inq_varid(ncid, 'lat', varid)
    ierr = nf90_put_var(ncid, varid, lat)
    ierr = nf90_inq_varid(ncid, 'lon', varid)
    ierr = nf90_put_var(ncid, varid, lon)
    ierr = nf90_inq_varid(ncid, 'mask', varid)
    ierr = nf90_put_var(ncid, varid, mask_out)
    if (ierr /= nf90_noerr) stop 'Error writing variable'
    ierr = nf90_inq_varid(ncid, 'll_lon', varid)
    ierr = nf90_put_var(ncid, varid, ll_lon)
    ierr = nf90_inq_varid(ncid, 'll_lat', varid)
    ierr = nf90_put_var(ncid, varid, ll_lat)
    ierr = nf90_inq_varid(ncid, 'd_lon', varid)
    ierr = nf90_put_var(ncid, varid, d_lon)
    ierr = nf90_inq_varid(ncid, 'd_lat', varid)
    ierr = nf90_put_var(ncid, varid, d_lat)

    ! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= nf90_noerr) stop 'Error closing file'

end program ascat_mask_maker
