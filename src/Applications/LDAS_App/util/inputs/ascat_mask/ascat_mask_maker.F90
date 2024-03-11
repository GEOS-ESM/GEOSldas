program ascat_mask_maker

    use netcdf

    implicit none

    integer :: ncid, varid, dimid, ierr, len, N_gpi, dimids(2)
    integer :: i, j, closest_index

    integer(kind=1), dimension(:,:), allocatable :: mask_out
    integer(kind=1)                              :: missing_value

    real, dimension(:), allocatable :: asc_lon, asc_lat, cold_mask, wet_mask, veg_mask, subsurface_mask, combined_mask
    real, dimension(:), allocatable :: lon, lat, distances
    real                            :: d_lon, d_lat, ll_lon, ll_lat

    ! Open the NetCDF file
    ierr = nf90_open('/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/clsm/subsurface_scattering_ASCAT_ERA5_Land.nc', nf90_nowrite, ncid)
    if (ierr /= nf90_noerr) stop 'Error opening file'

    ! Get the dimension ID
    ierr = nf90_inq_dimid(ncid, 'gpi', dimid)
    if (ierr /= nf90_noerr) stop 'Error getting dimension ID'
    
    ! Get the length of the dimension
    ierr = nf90_inquire_dimension(ncid, dimid, len = N_gpi)
    if (ierr /= nf90_noerr) stop 'Error inquiring dimension'
    
    print*, 'N_gpi = ', N_gpi

    ! Allocate the arrays
    allocate(asc_lon(N_gpi))
    allocate(asc_lat(N_gpi))
    allocate(cold_mask(N_gpi))
    allocate(wet_mask(N_gpi))
    allocate(veg_mask(N_gpi))
    allocate(subsurface_mask(N_gpi))
    allocate(combined_mask(N_gpi))

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

    ! Combine the masks
    where (wet_mask == 1)
       combined_mask = 1
    elsewhere
       combined_mask = subsurface_mask
    end where

    d_lon = 0.1
    d_lat = 0.1
    ll_lon = -180.0
    ll_lat = -90.0
    missing_value = -128

    allocate(lon(int(360.0  / d_lon)))
    allocate(lat(int(180.0  / d_lat)))

    lon = [(ll_lon + i * d_lon, i = 0, size(lon) - 1)]
    lat = [(ll_lat + i * d_lat, i = 0, size(lat) - 1)]

    allocate(mask_out(size(lon), size(lat)))

    do i = 1, size(lon)
        print*, lon(i)
        do j = 1, size(lat)
            allocate(distances(size(asc_lon)))
            distances = sqrt((asc_lon - lon(i))**2 + (asc_lat - lat(j))**2)
            closest_index = minloc(distances, dim = 1)
            if (distances(closest_index) > 0.14) then
                mask_out(i, j) = missing_value
            else
                mask_out(i, j) = combined_mask(closest_index)
            end if
            deallocate(distances)
        end do
    end do

    ! Write out the mask to netcdf
    ierr = nf90_create('ascat_combined_mask_p1_f90.nc', nf90_clobber, ncid)
    if (ierr /= nf90_noerr) stop 'Error creating file'

    ! Define the dimensions
    ierr = nf90_def_dim(ncid, 'lon', size(lon), dimids(1))
    ierr = nf90_def_dim(ncid, 'lat', size(lat), dimids(2))

    ! Define the variables
    ierr = nf90_def_var(ncid, 'lat', nf90_real, dimids(2), varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'latitude')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'latitude')
    ierr = nf90_put_att(ncid, varid, 'units', 'degrees_north')
    ierr = nf90_put_att(ncid, varid, 'axis', 'Y')

    ierr = nf90_def_var(ncid, 'lon', nf90_real, dimids(1), varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'longitude')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'longitude')
    ierr = nf90_put_att(ncid, varid, 'units', 'degrees_east')
    ierr = nf90_put_att(ncid, varid, 'axis', 'X')

    ierr = nf90_def_var(ncid, 'mask', nf90_byte, dimids, varid)
    ierr = nf90_put_att(ncid, varid, 'standard_name', 'subsurface_mask')
    ierr = nf90_put_att(ncid, varid, 'long_name', 'Mask accounting for subsurface scattering')
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