    function [] = write_netcdf_file(fname, colind, rowind,...
        av_angle_bin, data, asc_flag,...
        version, ...
        start_time, end_time, overwrite, N_out_fields, ...
        write_ind_latlon, data_product,...
        tile_id)  %last argument is optional

    int_precision   = 'NC_INT';      % precision of fortran tag
    float_precision = 'NC_DOUBLE';    % precision of data in input file

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

    % Calculate tmp_start_time tmp_end_time in days since January 1 1950
    d = datetime(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.min, start_time.sec);
    % Convert to serial date number
    serialNum = datenum(d);
    % Subtract serial date number of January 1, 1950
    daysSince1950 = serialNum - datenum('January 1, 1950');
    tmp_start_time = daysSince1950;
    
    d = datetime(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.min, end_time.sec);
    % Convert to serial date number
    serialNum = datenum(d);
    % Subtract serial date number of January 1, 1950
    daysSince1950 = serialNum - datenum('January 1, 1950');
    tmp_end_time = daysSince1950;

    % determine number of grid cells ; further check dimensions
    N_grid = size(data,2);
    N_angle= 1;
    N_pentad = 1;

    if (strcmp(write_ind_latlon,'latlon_id') && nargin == 14)
        if( size(tile_id,1) ~= N_grid )
            error('tile_id dimensions ??')
        end
        if ( size(tile_id,2) > 1)
            disp(['# subgridcells per gridcell: ',num2str(size(tile_id,2))]);
        end
    end

% create netCDF file
netcdf.setDefaultFormat('FORMAT_NETCDF4');
ncid = netcdf.create(fname, 'NETCDF4');

% define dimensions
dimid_grid = netcdf.defDim(ncid, 'grid', N_grid);
dimid_angle = netcdf.defDim(ncid, 'angle', N_angle);
dimid_tile = netcdf.defDim(ncid, 'tile', size(tile_id,2));
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
varid_tile_id = netcdf.defVar(ncid, 'tile_id', int_precision, [dimid_grid dimid_tile]);
varid_av_angle_bin = netcdf.defVar(ncid, 'av_angle_bin', float_precision, [dimid_angle]);
varid_colind = netcdf.defVar(ncid, 'colind', float_precision, [dimid_grid]);
varid_rowind = netcdf.defVar(ncid, 'rowind', float_precision, [dimid_grid]);

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
netcdf.putVar(ncid, varid_pentad, start_time.pentad);
netcdf.putVar(ncid, varid_start_time, tmp_start_time);
netcdf.putVar(ncid, varid_end_time, tmp_end_time);
netcdf.putVar(ncid, varid_N_grid, N_grid);
netcdf.putVar(ncid, varid_N_angle, N_angle);

if (~(strcmp(data_product,'scaling') && strcmp(write_ind_latlon,'latlon_id') && nargin == 14))
    netcdf.putVar(ncid, varid_tile_id, tile_id);
else
    netcdf.putVar(ncid, varid_tile_id, permute(tile_id, [2 1]));
end

netcdf.putVar(ncid, varid_av_angle_bin, squeeze(av_angle_bin(:)));

    if (N_grid >= 1)

netcdf.putVar(ncid, varid_colind, colind);
netcdf.putVar(ncid, varid_rowind, rowind);

 netcdf.putVar(grpid,varid_om,data(1,:));
netcdf.putVar(grpid,varid_ov,data(2,:));
netcdf.putVar(grpid,varid_mm,data(3,:));
netcdf.putVar(grpid,varid_mv,data(4,:));
netcdf.putVar(grpid,varid_ndata,data(5,:));

    else
        netcdf.putVar(ncid, varid_colind, 0.0);
netcdf.putVar(ncid, varid_rowind, 0.0);

 netcdf.putVar(grpid,varid_om,-999.0);
netcdf.putVar(grpid,varid_ov,-999.0);
netcdf.putVar(grpid,varid_mm,-999.0);
netcdf.putVar(grpid,varid_mv,-999.0);
netcdf.putVar(grpid,varid_ndata,-999.0);
    end

% close netCDF file
netcdf.close(ncid);

    end



