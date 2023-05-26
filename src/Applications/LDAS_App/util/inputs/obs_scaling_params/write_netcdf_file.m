    function [] = write_netcdf_file(fname, colind, rowind,...
        lon_out, lat_out,... 
        av_angle_bin, data, asc_flag,...
        pentad, ...
        start_time, end_time, overwrite, N_out_fields, ...
        write_ind_latlon, data_product,...
        obsnum)  %last argument is optional

    int_precision   = 'NC_INT';      % precision of fortran tag
    float_precision = 'NC_DOUBLE';    % precision of data in input file
    
    version = 0;

    % check dimensions
    if size(data,1)~=N_out_fields
        error('ERROR: size of data incompatible with N_out_fields')
    end

    % check for presence of optional input "overwrite"
    if ~exist('overwrite','var')
        overwrite = 0;        % default: do NOT overwrite existing files
    end

    % check if file exists
    if exist(fname,'file')
        if overwrite==0
            disp(['RETURNING!!! -- NOT OVERWRITING EXISTING FILE ', fname])
            return
        else
            disp(['OVERWRITING ', fname])
        end
    else
        disp(['writing ', fname])
    end
    
    % Convert the cell arrays to matrices using cell2mat (to deal with all_pentad case)
    year_mat = cell2mat({start_time.year}');
    month_mat = cell2mat({start_time.month}');
    day_mat = cell2mat({start_time.day}');
    hour_mat = cell2mat({start_time.hour}');
    min_mat = cell2mat({start_time.min}');
    sec_mat = cell2mat({start_time.sec}');
    % Use the matrices as input to the datetime function
    d = datetime(year_mat, month_mat, day_mat, hour_mat, min_mat, sec_mat);
    % Convert to serial date number
    serialNum = datenum(d);
    % Subtract serial date number of January 1, 1950
    daysSince1950 = serialNum - datenum('January 1, 1950');
    tmp_start_time = daysSince1950;
    
    % Convert the cell arrays to matrices using cell2mat (to deal with all_pentad case)
    year_mat = cell2mat({end_time.year}');
    month_mat = cell2mat({end_time.month}');
    day_mat = cell2mat({end_time.day}');
    hour_mat = cell2mat({end_time.hour}');
    min_mat = cell2mat({end_time.min}');
    sec_mat = cell2mat({end_time.sec}');
    % Use the matrices as input to the datetime function
    d = datetime(year_mat, month_mat, day_mat, hour_mat, min_mat, sec_mat);
    % Convert to serial date number
    serialNum = datenum(d);
    % Subtract serial date number of January 1, 1950
    daysSince1950 = serialNum - datenum('January 1, 1950');
    tmp_end_time = daysSince1950;

    % determine number of grid cells ; further check dimensions
    N_grid = size(data,2);
    N_angle= 1;
    if ndims(data) ==2
        N_pentad = 1;
    else
        N_pentad = 73;
    end

    if (strcmp(write_ind_latlon,'latlon_id') && nargin == 14)
        if( size(obsnum,1) ~= N_grid )
            error('tile_id dimensions ??')
        end
        if ( size(obsnum,2) > 1)
            disp(['# subgridcells per gridcell: ',num2str(size(obsnum,2))]);
        end
    end

% create netCDF file
netcdf.setDefaultFormat('FORMAT_NETCDF4');
ncid = netcdf.create(fname, 'NETCDF4');

% define dimensions
dimid_grid = netcdf.defDim(ncid, 'grid', N_grid);
dimid_angle = netcdf.defDim(ncid, 'angle', N_angle);
dimid_tile = netcdf.defDim(ncid, 'tile', size(obsnum,2));
dimid_pentad = netcdf.defDim(ncid, 'pentad', N_pentad);

% define variables
varid_asc_flag = netcdf.defVar(ncid, 'asc_flag', int_precision, []);
varid_version = netcdf.defVar(ncid, 'version', int_precision, []);
varid_pentad = netcdf.defVar(ncid, 'pentad', int_precision, [dimid_pentad]);

varid_start_time = netcdf.defVar(ncid, 'start_time', float_precision, [dimid_pentad]);
netcdf.putAtt(ncid, varid_start_time, 'standard_name','start time');
netcdf.putAtt(ncid, varid_start_time, 'long_name','start time');
netcdf.putAtt(ncid, varid_start_time, 'axis','T');
netcdf.putAtt(ncid, varid_start_time, 'units','days since 1950-01-01 00:00:00.0 +0000');

varid_end_time = netcdf.defVar(ncid, 'end_time',float_precision, [dimid_pentad]);
netcdf.putAtt(ncid, varid_end_time, 'standard_name','end time');
netcdf.putAtt(ncid, varid_end_time, 'long_name','end time');
netcdf.putAtt(ncid, varid_end_time, 'axis','T');
netcdf.putAtt(ncid, varid_end_time, 'units','days since 1950-01-01 00:00:00.0 +0000');

varid_N_grid = netcdf.defVar(ncid, 'N_grid', int_precision, []);
varid_N_angle = netcdf.defVar(ncid, 'N_angle', int_precision, []);
varid_obsnum = netcdf.defVar(ncid, 'obs_num', int_precision, [dimid_grid dimid_tile]);
varid_av_angle_bin = netcdf.defVar(ncid, 'av_angle_bin', float_precision, [dimid_angle]);
varid_colind = netcdf.defVar(ncid, 'colind', float_precision, [dimid_grid]);
varid_rowind = netcdf.defVar(ncid, 'rowind', float_precision, [dimid_grid]);
varid_lon = netcdf.defVar(ncid, 'lon', float_precision, [dimid_grid]);
varid_lat = netcdf.defVar(ncid, 'lat', float_precision, [dimid_grid]);

% Create a new group called 'data'
groupname = 'data';
grpid = netcdf.defGrp(ncid,groupname);

varid_om = netcdf.defVar(grpid, 'o_mean', float_precision, [dimid_grid dimid_pentad]);
varid_ov = netcdf.defVar(grpid, 'o_std', float_precision, [dimid_grid dimid_pentad]);
varid_mm = netcdf.defVar(grpid, 'm_mean', float_precision, [dimid_grid dimid_pentad]);
varid_mv =  netcdf.defVar(grpid, 'm_std', float_precision, [dimid_grid dimid_pentad]);
varid_ndata =  netcdf.defVar(grpid, 'n_data', float_precision, [dimid_grid dimid_pentad]);

% end define mode
netcdf.endDef(ncid);

% write data
netcdf.putVar(ncid, varid_asc_flag, asc_flag);
netcdf.putVar(ncid, varid_version, version);
netcdf.putVar(ncid, varid_pentad, pentad);
netcdf.putVar(ncid, varid_start_time, tmp_start_time);
netcdf.putVar(ncid, varid_end_time, tmp_end_time);
netcdf.putVar(ncid, varid_N_grid, N_grid);
netcdf.putVar(ncid, varid_N_angle, N_angle);

if (~(strcmp(data_product,'scaling') && strcmp(write_ind_latlon,'latlon_id') && nargin == 14))
    netcdf.putVar(ncid, varid_obsnum, obsnum);
else
    netcdf.putVar(ncid, varid_obsnum, permute(obsnum, [2 1]));
end

netcdf.putVar(ncid, varid_av_angle_bin, squeeze(av_angle_bin(:)));

    if (N_grid >= 1)

netcdf.putVar(ncid, varid_colind, colind);
netcdf.putVar(ncid, varid_rowind, rowind);
netcdf.putVar(ncid, varid_lon, lon_out);
netcdf.putVar(ncid, varid_lat, lat_out);

if N_pentad ==1
    netcdf.putVar(grpid,varid_om,data(1,:));
    netcdf.putVar(grpid,varid_ov,data(2,:));
    netcdf.putVar(grpid,varid_mm,data(3,:));
    netcdf.putVar(grpid,varid_mv,data(4,:));
    netcdf.putVar(grpid,varid_ndata,data(5,:));
else
    netcdf.putVar(grpid,varid_om,data(1,:,:));
    netcdf.putVar(grpid,varid_ov,data(2,:,:));
    netcdf.putVar(grpid,varid_mm,data(3,:,:));
    netcdf.putVar(grpid,varid_mv,data(4,:,:));
    netcdf.putVar(grpid,varid_ndata,data(5,:,:));
end
    else
        netcdf.putVar(ncid, varid_colind, 0.0);
        netcdf.putVar(ncid, varid_rowind, 0.0);
        netcdf.putVar(ncid, varid_lon, 0.0);
        netcdf.putVar(ncid, varid_lat, 0.0)

 netcdf.putVar(grpid,varid_om,-999.0);
netcdf.putVar(grpid,varid_ov,-999.0);
netcdf.putVar(grpid,varid_mm,-999.0);
netcdf.putVar(grpid,varid_mv,-999.0);
netcdf.putVar(grpid,varid_ndata,-999.0);
    end

% close netCDF file
netcdf.close(ncid);

    end



