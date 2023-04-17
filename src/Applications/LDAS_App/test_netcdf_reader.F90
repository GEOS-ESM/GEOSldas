program read_netcdf
    use netcdf
    implicit none
    
    integer :: ncid, varid, dimid, status, i, j
    integer :: grid, angle, tile, pentad
    integer :: asc_flag, version, N_grid, N_angle
    real(kind=8) :: start_time(pentad), end_time(pentad)
    real(kind=8) :: o_mean(pentad, grid), o_std(pentad, grid)
    real(kind=8) :: m_mean(pentad, grid), m_std(pentad, grid)
    real(kind=8) :: n_data(pentad, grid), av_angle_bin(angle)
    real(kind=8) :: colind(grid), rowind(grid), lon(grid), lat(grid)
    
    ! Open the netCDF file
    status = nf_open("filename.nc", nf_nowrite, ncid)
    if (status /= nf90_noerr) then
      write(*,*) "Error opening netCDF file"
      stop
    endif
    
    ! Get dimension IDs
    status = nf_inq_dimid(ncid, "grid", dimid)
    status = nf_inq_dimlen(ncid, dimid, grid)
    status = nf_inq_dimid(ncid, "angle", dimid)
    status = nf_inq_dimlen(ncid, dimid, angle)
    status = nf_inq_dimid(ncid, "tile", dimid)
    status = nf_inq_dimlen(ncid, dimid, tile)
    status = nf_inq_dimid(ncid, "pentad", dimid)
    status = nf_inq_dimlen(ncid, dimid, pentad)
    
    ! Get variable IDs
    status = nf_inq_varid(ncid, "asc_flag", varid)
    status = nf_inq_varid(ncid, "version", varid)
    status = nf_inq_varid(ncid, "pentad", varid)
    status = nf_inq_varid(ncid, "start_time", varid)
    status = nf_inq_varid(ncid, "end_time", varid)
    status = nf_inq_varid(ncid, "N_grid", varid)
    status = nf_inq_varid(ncid, "N_angle", varid)
    status = nf_inq_varid(ncid, "obs_num", varid)
    status = nf_inq_varid(ncid, "av_angle_bin", varid)
    status = nf_inq_varid(ncid, "colind", varid)
    status = nf_inq_varid(ncid, "rowind", varid)
    status = nf_inq_varid(ncid, "lon", varid)
    status = nf_inq_varid(ncid, "lat", varid)
    status = nf_inq_varid(ncid, "o_mean", varid)
    status = nf_inq_varid(ncid, "o_std", varid)
    status = nf_inq_varid(ncid, "m_mean", varid)
    status = nf_inq_varid(ncid, "m_std", varid)
    status = nf_inq_varid(ncid, "n_data", varid)
    
    ! Read variables
    status = nf_get_var(ncid, nf_inq_varid(ncid, "asc_flag", varid), asc_flag)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "version", varid), version)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "pentad", varid), pentad)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "start_time", varid), start_time)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "end_time", varid), end_time)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "N_grid", varid), N_grid)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "N_angle", varid), N_angle)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "obs_num", varid), obs_num)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "av_angle_bin", varid), av_angle_bin)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "colind", varid), colind)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "rowind", varid), rowind)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "lon", varid), lon)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "lat", varid), lat)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "o_mean", varid), o_mean)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "o_std", varid), o_std)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "m_mean", varid), m_mean)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "m_std", varid), m_std)
    status = nf_get_var(ncid, nf_inq_varid(ncid, "n_data", varid), n_data)
    
    ! Close the netCDF file
    status = nf_close(ncid)
    
    ! Print some information
    write(,) "Number of grids: ", grid
    write(,) "Number of pentads: ", pentad
    write(,) "Start time for pentad 1: ", start_time(1)
    write(,) "End time for pentad 1: ", end_time(1)
    
    end program read_netcdf