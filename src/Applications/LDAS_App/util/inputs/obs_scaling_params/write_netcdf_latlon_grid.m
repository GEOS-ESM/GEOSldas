function [] = write_netcdf_latlon_grid( fname, colind, rowind, ll_lons, ll_lats,                 ...
    data, pentad, start_time, end_time, overwrite, N_out_fields, ll_lon, ll_lat, d_lon, d_lat )        

    int_precision   = 'NC_INT';       % precision of fortran tag
    float_precision = 'NC_DOUBLE';    % precision of data in input file

    % Define the compression level (0-9, where 0 is no compression and 9 is maximum compression)
    compression_level = 5;
    
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
    year_mat  = cell2mat({start_time.year}');
    month_mat = cell2mat({start_time.month}');
    day_mat   = cell2mat({start_time.day}');
    hour_mat  = cell2mat({start_time.hour}');
    min_mat   = cell2mat({start_time.min}');
    sec_mat   = cell2mat({start_time.sec}');
    % Use the matrices as input to the datetime function
    d = datetime(year_mat, month_mat, day_mat, hour_mat, min_mat, sec_mat);
    % Convert to serial date number
    serialNum = datenum(d);
    % Subtract serial date number of January 1, 1950
    daysSince1950 = serialNum - datenum('January 1, 1950');
    tmp_start_time = daysSince1950;
    
    % Convert the cell arrays to matrices using cell2mat (to deal with all_pentad case)
    year_mat  = cell2mat({end_time.year}');
    month_mat = cell2mat({end_time.month}');
    day_mat   = cell2mat({end_time.day}');
    hour_mat  = cell2mat({end_time.hour}');
    min_mat   = cell2mat({end_time.min}');
    sec_mat   = cell2mat({end_time.sec}');
    % Use the matrices as input to the datetime function
    d = datetime(year_mat, month_mat, day_mat, hour_mat, min_mat, sec_mat);
    % Convert to serial date number
    serialNum = datenum(d);
    % Subtract serial date number of January 1, 1950
    daysSince1950 = serialNum - datenum('January 1, 1950');
    tmp_end_time = daysSince1950;
    
    N_lon = length(ll_lons);
    N_lat = length(ll_lats);
    
    % Have we got multiple pentads
    if ismatrix(data)
        N_pentad = 1;
    else
        N_pentad = size(data,3);
    end

    % create netCDF file
    netcdf.setDefaultFormat('FORMAT_NETCDF4');
    ncid = netcdf.create(fname, 'NETCDF4');
    
    % define dimensions
    dimid_pentad = netcdf.defDim(ncid, 'pentad', N_pentad);
    dimid_lon    = netcdf.defDim(ncid, 'lon',    N_lon);
    dimid_lat    = netcdf.defDim(ncid, 'lat',    N_lat);
    
    % define variables
    
    varid_version    = netcdf.defVar(ncid, 'version', int_precision,  []);

    varid_ll_lon     = netcdf.defVar(ncid, 'll_lon',     float_precision, []);
    netcdf.putAtt(ncid, varid_ll_lon,      'standard_name', 'longitude of lower left corner');
    netcdf.putAtt(ncid, varid_ll_lon,      'long_name',     'longitude of lower left corner');
    netcdf.putAtt(ncid, varid_ll_lon,      'units',         'degrees_east');
    netcdf.putAtt(ncid, varid_ll_lon,      'axis',          'X');
    
    varid_ll_lat     = netcdf.defVar(ncid, 'll_lat',     float_precision, []);
    netcdf.putAtt(ncid, varid_ll_lat,      'standard_name', 'latitude of lower left corner');
    netcdf.putAtt(ncid, varid_ll_lat,      'long_name',     'latitude of lower left corner');
    netcdf.putAtt(ncid, varid_ll_lat,      'units',         'degrees_north');
    netcdf.putAtt(ncid, varid_ll_lat,      'axis',          'Y');
    
    varid_d_lon      = netcdf.defVar(ncid, 'd_lon',      float_precision, []);
    netcdf.putAtt(ncid, varid_d_lon,       'standard_name', 'longitude grid spacing');
    netcdf.putAtt(ncid, varid_d_lon,       'long_name',     'longitude grid spacing');
    netcdf.putAtt(ncid, varid_d_lon,       'units',         'degrees');
    netcdf.putAtt(ncid, varid_d_lon,       'axis',          'X');
    
    varid_d_lat      = netcdf.defVar(ncid, 'd_lat',      float_precision, []);
    netcdf.putAtt(ncid, varid_d_lat,       'standard_name', 'latitude grid spacing');
    netcdf.putAtt(ncid, varid_d_lat,       'long_name',     'latitude grid spacing');
    netcdf.putAtt(ncid, varid_d_lat,       'units',         'degrees');
    netcdf.putAtt(ncid, varid_d_lat,       'axis',          'Y');
    
    varid_pentad     = netcdf.defVar(ncid, 'pentad',     int_precision,   [dimid_pentad]);
    netcdf.putAtt(ncid, varid_pentad,      'standard_name', 'pentad');
    netcdf.putAtt(ncid, varid_pentad,      'long_name',     'pentad');
    netcdf.putAtt(ncid, varid_pentad,      'units',         '1');
    netcdf.putAtt(ncid, varid_pentad,      'axis',          'T');
    
    varid_start_time = netcdf.defVar(ncid, 'start_time', float_precision, [dimid_pentad]);
    netcdf.putAtt(ncid, varid_start_time,  'standard_name', 'start time');
    netcdf.putAtt(ncid, varid_start_time,  'long_name',     'start time');
    netcdf.putAtt(ncid, varid_start_time,  'axis',          'T');
    netcdf.putAtt(ncid, varid_start_time,  'units',         'days since 1950-01-01 00:00:00.0 +0000');
    
    varid_end_time   = netcdf.defVar(ncid, 'end_time',   float_precision, [dimid_pentad]);
    netcdf.putAtt(ncid, varid_end_time,    'standard_name', 'end time');
    netcdf.putAtt(ncid, varid_end_time,    'long_name',     'end time');
    netcdf.putAtt(ncid, varid_end_time,    'axis',          'T');
    netcdf.putAtt(ncid, varid_end_time,    'units',         'days since 1950-01-01 00:00:00.0 +0000');
       
    varid_om         = netcdf.defVar(ncid, 'o_mean',     float_precision, [dimid_lat dimid_lon dimid_pentad]);
    netcdf.defVarDeflate(ncid,varid_om,true,true,compression_level);
    netcdf.putAtt(ncid, varid_om,          'standard_name', 'observation mean');
    netcdf.putAtt(ncid, varid_om,          'long_name',     'Observation mean for pentad calculated over all years for window length');
    netcdf.putAtt(ncid, varid_om,          'units',         'Degree of saturation (0-1)');
    
    varid_ov         = netcdf.defVar(ncid, 'o_std',      float_precision, [dimid_lat dimid_lon dimid_pentad]);
    netcdf.defVarDeflate(ncid,varid_ov,true,true,compression_level);
    netcdf.putAtt(ncid, varid_ov,          'standard_name', 'observation standard deviation');
    netcdf.putAtt(ncid, varid_ov,          'long_name',     'Observation standard deviation for pentad calculated over all years for window length');
    netcdf.putAtt(ncid, varid_ov,          'units',         'Degree of saturation (0-1)');
    
    varid_mm         = netcdf.defVar(ncid, 'm_mean',     float_precision, [dimid_lat dimid_lon dimid_pentad]);
    netcdf.defVarDeflate(ncid,varid_mm,true,true,compression_level);
    netcdf.putAtt(ncid, varid_mm,          'standard_name', 'model mean');
    netcdf.putAtt(ncid, varid_mm,          'long_name',     'Model mean for pentad calculated over all years for window length');
    netcdf.putAtt(ncid, varid_mm,          'units',         'Surface soil moisture (m^3 m^-3)'); 
    
    varid_mv         = netcdf.defVar(ncid, 'm_std',      float_precision, [dimid_lat dimid_lon dimid_pentad]);
    netcdf.defVarDeflate(ncid,varid_mv,true,true,compression_level);
    netcdf.putAtt(ncid, varid_mv,          'standard_name', 'model standard deviation');
    netcdf.putAtt(ncid, varid_mv,          'long_name',     'Model standard deviation for pentad calculated over all years for window length');
    netcdf.putAtt(ncid, varid_mv,          'units',         'Surface soil moisture (m^3 m^-3)');
    
    varid_mi         = netcdf.defVar(ncid, 'm_min',      float_precision, [dimid_lat dimid_lon]);
    netcdf.defVarDeflate(ncid,varid_mi,true,true,compression_level);
    netcdf.putAtt(ncid, varid_mi,          'standard_name', 'model minimum');
    netcdf.putAtt(ncid, varid_mi,          'long_name',     'Model minimum calculated over all years');
    netcdf.putAtt(ncid, varid_mi,          'units',         'Surface soil moisture (m^3 m^-3)');
    
    varid_ma         = netcdf.defVar(ncid, 'm_max',      float_precision, [dimid_lat dimid_lon]);
    netcdf.defVarDeflate(ncid,varid_ma,true,true,compression_level);
    netcdf.putAtt(ncid, varid_ma,          'standard_name', 'model maximum');
    netcdf.putAtt(ncid, varid_ma,          'long_name',     'Model maximum calculated over all years');
    netcdf.putAtt(ncid, varid_ma,          'units',         'Surface soil moisture (m^3 m^-3)');
    
    varid_ndata      = netcdf.defVar(ncid, 'n_data',     float_precision, [dimid_lat dimid_lon dimid_pentad]);
    netcdf.defVarDeflate(ncid,varid_ndata,true,true,compression_level);
    netcdf.putAtt(ncid, varid_ndata,       'standard_name', 'number of data points');
    netcdf.putAtt(ncid, varid_ndata,       'long_name',     'Number of data points for pentad calculated over all years for window length');
    netcdf.putAtt(ncid, varid_ndata,       'units',         '1');
    
    % end define mode
    netcdf.endDef(ncid);
    
    % write data
    netcdf.putVar(ncid, varid_pentad,     pentad);
    netcdf.putVar(ncid, varid_start_time, tmp_start_time);
    netcdf.putVar(ncid, varid_end_time,   tmp_end_time);
    
    netcdf.putVar(ncid, varid_ll_lon,     ll_lon);
    netcdf.putVar(ncid, varid_ll_lat,     ll_lat);
    netcdf.putVar(ncid, varid_d_lon,      d_lon);
    netcdf.putVar(ncid, varid_d_lat,      d_lat);
    
    if N_pentad ==1
    
        data_out = ones(N_out_fields,N_lat,N_lon         ) * -999.0;
    
        for n = 1:N_out_fields
            for i = 1:length(colind)
                data_out(n,rowind(i),colind(i)) = data(n,i);
            end
        end
    
        netcdf.putVar(ncid,varid_om,        data_out(1,:,:)        );
        netcdf.putVar(ncid,varid_ov,        data_out(2,:,:)        );
        netcdf.putVar(ncid,varid_mm,        data_out(3,:,:)        );
        netcdf.putVar(ncid,varid_mv,        data_out(4,:,:)        );
        netcdf.putVar(ncid,varid_mi,        data_out(6,:,:)        );
        netcdf.putVar(ncid,varid_ma,        data_out(7,:,:)        );
        netcdf.putVar(ncid,varid_ndata,     data_out(5,:,:)        );
        
    else
    
        data_out = ones(N_out_fields,N_lat,N_lon,N_pentad) * -999.0;
    
        for n = 1:N_out_fields
            for i = 1:length(colind)
                data_out(n,rowind(i),colind(i),:) = data(n,i,:);
            end
        end
    
        netcdf.putVar(ncid,varid_om,        data_out(1,:,:,:)      );
        netcdf.putVar(ncid,varid_ov,        data_out(2,:,:,:)      );
        netcdf.putVar(ncid,varid_mm,        data_out(3,:,:,:)      );
        netcdf.putVar(ncid,varid_mv,        data_out(4,:,:,:)      );
        netcdf.putVar(ncid,varid_ndata,     data_out(5,:,:,:)      );
        netcdf.putVar(ncid,varid_ma,    max(data_out(7,:,:,:),[],4));      % Max over all pentads, always only 2D
        
        min_data                          = squeeze(data_out(6, :, :,:)); 
        min_data(min_data < -9998)        = NaN;                           % Switch current missing value to NaN before calculating min
        min_data_out                      = min(min_data,[],3); 
        min_data_out(isnan(min_data_out)) = -9999.;
        
        netcdf.putVar(ncid,varid_mi,        min_data_out);                 % Min over all pentads, always only 2D
    end

    % close netCDF file
    netcdf.close(ncid);
    
end

% ================ EOF ========================================================
