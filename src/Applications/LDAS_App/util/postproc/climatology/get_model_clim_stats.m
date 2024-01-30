function [] = get_model_clim_stats( fieldname,              ...
    exp_path, exp_run, domain,                              ...
    start_year, end_year, start_HHSS,                       ...
    out_freq, w_out_freq, file_tag,                         ...
    out_wetness,                                            ...
    N_data_min, time_of_day_in_hours )

%=======================================================================
%
% Compute mean, stdv and CDF-stats of simulated fields from 
% - tile-based "nc4" files 
% - or gridded "h5" files
% for a selection of fieldname.
%
% The main purpose of this function is to aggregate the information 
% from the "cat" files so that the climatology statistics can 
% be used for the computation of percentiles.
%
%
% Need to fix those file_tag-s... at some point, we had to calculate cli-files 
% based on all different file types... ugly
%
%
% One file with statistics is generated for every pentad or month.
%
% fieldname:    (single) model field to be processed 
% start_year:   start year for each month (12 entries!)
% end_year:     end year for each month (12 entries!)
% start_HHSS:   hour, min, sec of first file (depends on output resolution)
% out_freq:     'monthly' or 'pentad' climatology files
% w_out_freq:   number of months or pentads used in the temporal smoothing
% N_data_min:   minimum number of data points to calculate a good stat
%
% 'cli'-file has a similar format as the SMOS-data files and
% observation scaling files for assimilation.
%
% ==HEADER==
% N_tiles N_fields N_stat
% fields
%
% ==DATA==
%
% %%%% tile ID --> for all tiles (1:N_tiles, not sorted by tile_id!)
% lat
% lon
%
% [for fields]
%   mean
%   std
%   min
%   max
%   N_data
%   CDF_parameter_1     OR  UL_1st_percentile
%   CDF_parameter_2         UL_2st_percentile  
%   CDF_parameter_...       UL_...  
%   CDF_parameter_N         UL_99st_percentile  
% [end fields]
%
% GDL, feb 2014
%
% -------------------------------------------------------------------
% begin user-defined inputs
% -------------------------------------------------------------------

nodata     = -9999;    
nodata_tol = 1e-4;     

ens_tag = 'ens0000';

dtstep  = 3*60*60; % hardwired 3-hourly if x-hourly output is given

% when reading h5-data

datagroup_name = '/Geophysical_Data/';

% output specs

overwrite   =  1;

percentiles = [1:99];

N_stat      = 5;  %mean, stdv, min, max, N_data
N_CDF       = length(percentiles);

N_stat_CDF  = N_stat+N_CDF; 

write_ind_latlon = 'latlon_id'; %'latlon';

% -------------------------------------------------------------------
% end user-defined inputs
% -------------------------------------------------------------------

N_fields   = 1;

fieldno    = 1;
                                  
field_tag = ['_',fieldname];

% -------------------------------------------------------------------

% determine number of entries in smoothing time window

if (std(end_year-start_year) ~= 0)
    error('same number of years should contribute to each month')
end

if strcmp(out_freq,'pentad')
    n_days = 5;  
elseif strcmp(out_freq,'monthly')
    n_days = 365/12.0; % could adjust this and work w/ min_days=28, max_days=31
end

disp(['smoothing window is ',num2str(n_days*w_out_freq),' days']);

n_time_count = round(w_out_freq * n_days * (max(end_year-start_year)+1) *...
                     ((24*60*60)./dtstep));

% auxiliary start-end_time to get 1 year climatology:
% - loop over all days in any non-leap year (365 days)
% - make sure to loop into the next year to cover all climatology pentads
% or months w/ the complete smoothing window

start_time        = start_HHSS;
start_time.year   = 2014; 
start_time.month  = 1; 
start_time.day    = 1; 

start_time        = get_dofyr_pentad(start_time); %ini correct pentad

end_time          = augment_date_time((365 + w_out_freq * n_days)*24*60*60, ...
                                      start_time);

% effective period used in climatology calculation

start_time_true       = start_time;
start_time_true.year  = min(start_year);
tmp                   = find(start_year==min(start_year));
start_time_true.month = tmp(1);
start_time_true       = get_dofyr_pentad( start_time_true );

end_time_true         = end_time;
end_time_true.year    = max(end_year);
tmp                   = find(start_year==max(start_year));
end_time_true.month   = tmp(end);
end_time_true.day     = days_in_month( end_time_true.year, end_time_true.month);
end_time_true         = get_dofyr_pentad( end_time_true );

% assemble input and output paths

inpath  = [ exp_path, '/', exp_run, '/output/', domain ];
outpath = [ inpath, '/stats/cli/'   ];

% create outpath if it doesn't exist

if ~exist(outpath,'dir')
  eval(['!mkdir -p ', outpath]);
end

% -------------------------------------------------------------		  

% load catchment coordinates

fname = [inpath, '/rc_out/', exp_run, '.ldas_tilecoord.bin'];

[ tile_coord ] = read_tilecoord( fname );

N_tile = tile_coord.N_tile;

lat_out            = tile_coord.com_lat;
lon_out            = tile_coord.com_lon;
tile_coord_tile_id = tile_coord.tile_id;

% determine if conversion of soil moisture variables to wetness is
% needed;  if yes, get porosity from cat_param file

convert_to_wetness = 0;

if contains(field_tag,'sm_') & ~contains(field_tag,'wet') & out_wetness
  
  field_tag          = [field_tag,'_wet'];
  convert_to_wetness = 1;
  
  catfname = [exp_path,'/',exp_run,'/output/',domain,'/rc_out/','/Y2015/M01/',...
      exp_run,'.ldas_catparam.','20150101_0000','z.bin'];

  cat_param = read_catparam( catfname, N_tile );
  
  poros = cat_param.poros; clear cat_param
    
end
    
% -------------------------------------------------------------

% assemble output file name

if strcmp(out_freq,'pentad')

  fname_out_base = [ outpath, '/', 'cli_',             ...
          num2str(start_time_true.year),  '_p',        ...
          num2str(start_time_true.pentad),'_',         ...
          num2str(end_time_true.year),    '_p',        ...
          num2str(end_time_true.pentad),  '_',         ...
          'W_', num2str(w_out_freq),'p_',              ...
          'Nmin_',   num2str(round(N_data_min),'%d'),  ...
          field_tag];

else strcmp(out_freq,'monthly')

  fname_out_base = [ outpath, '/', 'cli_',             ...
          num2str(start_time_true.year),  '_M',        ...
          num2str(start_time_true.month),'_',          ...
          num2str(end_time_true.year),    '_M',        ...
          num2str(end_time_true.month),  '_',          ...
          'W_', num2str(w_out_freq),'M_',              ...
	  'Nmin_',   num2str(round(N_data_min),'%d'),  ...
          field_tag];

end

if exist( 'time_of_day_in_hours', 'var')
  
  fname_out_base = [fname_out_base, '_', num2str(time_of_day_in_hours,'%2.2d'), 'z'];
  
end

% -------------------------------------------------------------

% initialize output statistics

hist_data  = zeros(N_tile,N_fields,n_time_count);

time_count = 0;

data_out   = NaN+zeros(N_fields,N_tile,N_stat_CDF);

% -------------------------------------------------------------		  

disp('climatology calculation')

time_new = start_time;

while 1
    
  if (time_new.year  == end_time.year  &...
      time_new.month == end_time.month &...
      time_new.day   == end_time.day   &...
      time_new.hour  == end_time.hour  &...
      time_new.min   == end_time.min   &...
      time_new.sec   == end_time.sec  )
    break
  end

  % augment date_time

  time_old   = time_new;
  pentad_old = time_new.pentad;
  month_old  = time_new.month;

  time_new   = augment_date_time(dtstep, time_old);
  pentad_new = time_new.pentad;
  month_new  = time_new.month;

  % check if diurnal stats are needed

  if exist('time_of_day_in_hours','var')
    tmp_hour = time_of_day_in_hours;
  else
    tmp_hour = time_old.hour;   % all hours of day will be included
  end

  if time_old.hour==tmp_hour

    minute  = time_old.min;  % floor( (seconds_in_day-hour*3600)/60 );
    seconds = time_old.sec;  % seconds_in_day-hour*3600-minute*60;

    if (seconds~=0)
      error('something is wrong! (seconds~=0)')
    end

    for year = start_year(time_old.month):end_year(time_old.month)

      time_count =  time_count+1;
      
      YYYYMMDD = [ num2str(year,           '%4.4d'),    ...
                   num2str(time_old.month, '%2.2d'),    ...
                   num2str(time_old.day,   '%2.2d')  ];    

      HHMM     = [ num2str(time_old.hour,  '%2.2d'),    ...
                   num2str(time_old.min,   '%2.2d')  ]; 

      fname = [ inpath,                                 ...
        '/cat/', ens_tag,                               ...
        '/Y',   num2str(year, '%4.4d'),                 ...
        '/M',   num2str(time_old.month,'%2.2d'),        ...
        '/', exp_run,                                   ...
        '.', file_tag, '.', YYYYMMDD, '.nc4' ];

      if ~exist(fname,'file')
        
        % try again with "_[HHMM]z" inserted into file name
        
        fname = [ inpath,                               ...
        '/cat/', ens_tag,                               ...
        '/Y',   num2str(year, '%4.4d'),                 ...
        '/M',   num2str(time_old.month,'%2.2d'),        ...
        '/', exp_run,                                   ...
        '.', file_tag, '.', YYYYMMDD, '_', HHMM, 'z.nc4' ];
      
      end

      disp(['reading ',fieldname,' from ',fname])
      data_tmp = ncread(fname, fieldname);

      if size(data_tmp,2) == 8   % hard-wired 3-hourly time step??
         tile_data_tmp = data_tmp(:,ceil(time_old.hour/3.)); clear data_tmp
      elseif size(data_tmp,2) == 1
         tile_data_tmp = data_tmp; clear data_tmp
      else
         error(['data size is incorrect from ', fname])
      end

      for s=1:N_fields

        if ~convert_to_wetness
          tile_data_tmp_1D = tile_data_tmp(:);
        else
          tile_data_tmp_1D = tile_data_tmp(:)./poros(:); 
        end
        
        good_data        = find(~(abs(tile_data_tmp_1D - nodata) < nodata_tol));

        tile_data_tmp_1D = tile_data_tmp_1D(good_data);

        %Keep a record of time series
        
        total_bin_good_ind = (time_count-1).*(N_fields*N_tile) + ...
                             (s-1)*N_tile+good_data';
        
        hist_data(total_bin_good_ind) = tile_data_tmp_1D;
            
      end

    end % loop through years

  end  % time_of_day_in_hours
 
  % check if output needs to be written
    
  if (time_count == n_time_count )
    
    % write output
   
    for s=1:N_fields

      %edges  = edge_min(s):edge_dx(s):edge_max(s);
    
      for tile=1:N_tile        

        tmp = reshape(squeeze(hist_data(tile, s, :)),1,[]);
        
        data_out(s,tile,1) = mean(      tmp,"omitnan");   % mean
        data_out(s,tile,2) = std(       tmp,"omitnan");   % stdv
        data_out(s,tile,3) = min(       tmp          );   % min
        data_out(s,tile,4) = max(       tmp          );   % max
        data_out(s,tile,5) = sum(~isnan(tmp)         );   % N_data
        
        % determine the CDF-parameters, or the edges for each
        % percentile
            
        perc = round(percentiles./100*data_out(s,tile,5));

        tmp  = sort(tmp);

        data_out(s,tile,N_stat+1:N_stat+N_CDF) = tmp(perc);
              
      end
    
    end  

    bad_ind  = find(data_out(:,:,5)<N_data_min);
    
    for ff=1:N_stat_CDF
        data_out(bad_ind,ff) = NaN;
    end 
    
    % remove NaN before writing a file
    
    data_out(isnan(data_out)) = nodata;
    data_out(isinf(data_out)) = nodata;
    
    date_time_tmp = augment_date_time( -floor(w_out_freq*n_days*(24*60*60)/2.0), time_old );

    fname_out = [fname_out_base, '_p', num2str(date_time_tmp.pentad,'%2.2d'), '.bin'];

    % check whether output file exists

    if (exist(fname_out)==2 && overwrite) 

      disp(['output file exists. overwriting', fname_out])

    elseif (exist(fname_out)==2 && ~overwrite) 

      disp(['output file exists. not overwriting. returning'])
      return

    else

      disp(['creating ', fname_out])

    end

    write_seqbin_clim_pctl_file(fname_out, lon_out, lat_out, ...
            data_out(:,:,:), fieldno(:,1), N_stat_CDF, ...  
            overwrite, ...
            write_ind_latlon, 'cli', tile_coord_tile_id)

    % re-initialize
    
    data_out   = NaN+0*data_out;  
    
    % shift time window

    shift      = n_days * (max(end_year-start_year)+1) *...
                 ((24*60*60)./dtstep);
             
    time_count                   = time_count-shift;
    hist_data(:,:,1:time_count)  = hist_data(:,:,(1+shift):n_time_count);
      
  end  % if (time_count == n_time_count )
    
end    % while 1 (time loop)

% ==================== EOF ==============================================
