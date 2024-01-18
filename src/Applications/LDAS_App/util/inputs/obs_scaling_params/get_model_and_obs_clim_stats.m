
function [] = get_model_and_obs_clim_stats( varname,               ...
    run_months, exp_path, exp_run, domain, start_year, end_year,   ...
    dt_assim, t0_assim, species, obs_param,                        ...
    hscale, inc_angle, int_Asc, w_days, Ndata_min, prefix,         ...
    convert_grid , time_of_day_in_hours )

%
% get_model_and_obs_clim_stats.m
%
% Compute mean, stdv of model and observations from tile-based 
% "innov" files for a selection of species.
%
% The main purpose of this function is to aggregate the information 
% from the "innov" files so that the climatology statistics can 
% be used in scaling of the observations before assimilation.
%
% One file with statistics is generated for every DOY (1,...,365).
% The temporal smoothing/averaging window (w_days) is given in days.
%
% Stats output file in a similar format as the SMOS-data files 
%
% ==HEADER==
% N_tiles N_angles
% angles
%
% ==DATA==
%
% %%%% tile ID --> for all tiles (1:N_tiles, not sorted by tile_id!)
% lat
% lon
%
% [for angles]
%   mean_obs_H
%   std_obs_H
%   mean_mod_H
%   std_mod_H
%   N_data_H
%
%   mean_obs_V
%   std_obs_V
%   mean_mod_V
%   std_mod_V
%   N_data_V
%
%   H_obs
%   H_mod
%   V_obs
%   V_mod
% [end angles]
%
% GDL, 10 sept 2012
%
% GDL, aug 2013: added 'convert_grid' (= EASEv2_M36, EASE_M36, ...)
%                to project the Obs (always M36 for SMOS) and Fcst to M36 
%                M09 obs are administered by tiles (0) that could be anywhere
%                around the center of the observed pixel (M36)
%                -----------
%                | X X X X |
%                | X O O X |
%                | X O O X |
%                | X X X X |
%                -----------
% GDL, jan 2014: the above issue that "any" M09 tile in the center (0)
%                could potentially administer the M36 obs is not true
%                anymore with later LDASsa-tags. 
%                => no need to pass on 'convert_grid' for tags later than
%                the summer of 2013
% reichle, qliu, 13 July 2022:
%                "convert_grid" is still needed to limit the number of tiles
%                in the scaling parameter file.  With "convert_grid" turned on, 
%                only the M09 tile to the northeast of the M36 center point 
%                is kept in the scaling parameter file, consistent with the 
%                "tmp_shift_lat" and "tmp_shift_lon" operations in the SMOS
%                and SMAP Fortran readers.  (It is not clear if this matlab 
%                function works properly if there is no M09 [land] tile 
%                immediately to the northeast of the M36 center point.  In 
%                such a case, the Fortran reader assigns the nearest M09
%                land tile as the tile that administers the obs.)
%                Presumably, scaling parameters for all M09 tiles could be kept
%                if they are stored in (compressed) nc4 format.  In this case,
%                the NaN values for the scaling parameters of 15 out of each 16
%                M09 tiles can be compressed to almost nothing.
%
% -------------------------------------------------------------------
% begin user-defined inputs
% -------------------------------------------------------------------

% obs species to be processed (see ens_upd_inputs.nml for a list)
%
% (only observation species that represent observations of the same
%  model prognostic or diagnostic can be processed together!)

nodata     = -9999;    
nodata_tol = 1e-4;     
  
% minimum number of data points to include in window statistics

% N_data_min = w_days/10.; % initial screening on minimum # points in a
%                          % window to calculate a mean or stdv
% include a final decision about "good" stats later, when merging years

% no-data-value for points that don't have good statistics

no_data_stats = -9999.;

disp('ASSUMING EASEv2 M36 observations');

if ~isempty(strfind(domain,'M36')) && isempty(strfind(obs_param(species(1)).descr, '_E'))
  tol = 1E-3;
else
  tol = 2;
end

% output specs

overwrite = 1;

Nf = 5;                                       %5 fields per polarization

N_out_fields = 2*Nf+4; %14; 

write_ind_latlon = 'latlon_id';   %'latlon';

N_angle = length(inc_angle);
N_pol   = 2;

tmp_shift_lon    = 0.01;
tmp_shift_lat    = 0.005;

store_all_M09inM36 = 0;
print_each_DOY     = 0;

% -------------------------------------------------------------------
% end user-defined inputs
% -------------------------------------------------------------------

% assemble input and output paths

%inpath  = [ exp_path, '/output/', exp_run, '/', domain ];
inpath  = [ exp_path, '/', exp_run, '/output/', domain ];

outpath = [ inpath, '/stats/z_score_clim/'   ];

% create outpath if it doesn't exist

if exist(outpath)~=2
  eval(['!mkdir -p ', outpath]);
end

% -------------------------------------------------------------

% assemble output file name

ind  = find(start_year == min(start_year));
mi_m = min(run_months(ind));
ind  = find(end_year == max(end_year));
ma_m = max(run_months(ind));

D(1) = 1;
P(1) = 1;
if mi_m > 1
  D(1) = sum(days_in_month( 2014, [1:mi_m-1]))+1;
  P(1) = ceil(D(1)/5);
end
if ma_m > 1
  D(2) = sum(days_in_month( 2014, [1:ma_m]));
else
  D(2) = 1;
end
P(2) = floor(D(2)/5);

if run_months(1) ~= run_months(end) && run_months(2) ~= run_months(end)
  disp('WARNING: incomplete pentad-windows; loop through additional months to get complete pentads');
end

fname_out_base = [ outpath, '/', prefix,                                      ...
                   num2str(min(start_year)),'_doy',num2str(D(1)),'_',         ...
                   num2str(max(end_year)),  '_doy',num2str(D(2)),             ...
                   '_hscale_', num2str(hscale,'%2.2f'), '_',                  ...
                   'W_', num2str(w_days),'d_Nmin_', num2str(Ndata_min)];
      
fname_out_base_p = [ outpath, '/', prefix,                                           ...
                     num2str(min(start_year)),'_p',num2str(P(1)),'_',                ...
                     num2str(max(end_year)),  '_p',num2str(P(2)),                    ...
                     '_hscale_', num2str(hscale,'%2.2f'), '_',                       ...
                     'W_', num2str(round(w_days/5)),'p_Nmin_', num2str(Ndata_min)];

%fname_out_base = [fname_out_base, spec_tag];

if (int_Asc == 1)
  Orbit_tag = '_A'; %'_Asc';
else
  Orbit_tag = '_D'; %'_Desc';
end

fname_out_base   = [fname_out_base,   Orbit_tag];
fname_out_base_p = [fname_out_base_p, Orbit_tag];

if exist( 'time_of_day_in_hours', 'var')
  
  fname_out_base   = [fname_out_base,   '_', num2str(time_of_day_in_hours,'%2.2d'), 'z'];
  fname_out_base_p = [fname_out_base_p, '_', num2str(time_of_day_in_hours,'%2.2d'), 'z'];
  
end

% -------------------------------------------------------------		  

% load catchment coordinates

fname = [inpath, '/rc_out/', exp_run, '.ldas_tilecoord.bin'];
fnameg= [inpath, '/rc_out/', exp_run, '.ldas_tilegrids.bin'];

[ tile_coord ] = read_tilecoord( fname  );
[ tile_grid ]  = read_tilegrids( fnameg );

N_tile = length(tile_coord.tile_id);

% -------------------------------------------------------------

% determine tiles to whose statistics the current obs will contribute to

disp('pre-computing index for regional averaging')

central_lat        = tile_coord.com_lat;
central_lon        = tile_coord.com_lon;
tile_coord_tile_id = tile_coord.tile_id;

if (exist('convert_grid'))
  
  %1) convert to M36 EASE indices
  %2) convert back to lat/lon at center of obs
  if (~isempty(strfind(convert_grid, 'M36')) && ~isempty(strfind(convert_grid, 'EASEv2')))
    gridid = 'M36';
    [central_row,central_col] = EASEv2_latlon2ind(central_lat,central_lon,gridid,1);
    [central_lat,central_lon] = EASEv2_ind2latlon(central_row,central_col,gridid);
  elseif (~isempty(strfind(convert_grid, 'M36')) && ~isempty(strfind(convert_grid, 'EASEv1')))
    error('Must provide smapeasev1_latlon2ind() and smapeasev1_ind2latlon()!')
    gridid = 'M36';
    [central_row,central_col] = smapeasev1_latlon2ind(central_lat,central_lon,gridid);
    [central_lat,central_lon] = smapeasev1_ind2latlon(central_row,central_col,gridid);
  else
    error(['Unable to convert to ',convert_grid])
  end
  
  row_col_tmp         = [central_row central_col];
  [unique_rc, ia, ic] = unique(row_col_tmp,'rows');
  
  max_Hx_c = length(find(mode(ic)==ic));
  
  %know which exact M09 tiles are actually administering the obs 
  %-------------------
  tmp_lon = central_lon(ia)+tmp_shift_lon;
  tmp_lat = central_lat(ia)+tmp_shift_lat;
  
  [N_tile_in_cell_ij, tile_num_in_cell_ij] = get_tile_num_in_cell_ij(      ...
      tile_coord, tile_grid);
  
  this_FOV = 20;
  option   = 'FOV_in_km';
  %overwrite ia with actual administering tile number            
  [ia] = get_tile_num_for_obs( tile_coord, tile_grid,                      ...
                               N_tile_in_cell_ij, tile_num_in_cell_ij,     ...
                               option, this_FOV, tmp_lat, tmp_lon);
  
   ia = ia(ia>0 & ~isnan(ia));
   
   obsnum     = NaN+zeros(length(ic),1);
   obsnum(ia) = [1:length(ia)];        
   
   N_tile_obs = length(ia);
   
   %-------------------
   
   if store_all_M09inM36
     
     %Not maintained/elaborated
     tile_coord_tile_id = zeros(N_tile_obs,max_Hx_c);
     
     disp(['centralizing obs on ',convert_grid,' grid before doing stats: max ',num2str(max_Hx_c),'tiles per obs cell'])
     
     for i=1:N_tile_obs
       
       tmp_ind = find(row_col_tmp(:,1) == unique_rc(i,1)  & row_col_tmp(:,2) == unique_rc(i,2));
       
       tile_coord_tile_id(i,1:length(tmp_ind)) = tile_coord.tile_id(tmp_ind);
       
     end
     
   else
     
     tile_coord_tile_id = tile_coord.tile_id(ia); 
     
   end
   
else
  
  N_tile_obs = N_tile;
  ia         = 1:N_tile;
  ic         = 1:N_tile;
  obsnum     = 1:N_tile;        
  
end

lon_out = tile_coord.com_lon(ia); %NaN+zeros(N_tile,1);
lat_out = tile_coord.com_lat(ia); %NaN+zeros(N_tile,1);

if hscale>0
  
  for i=1:N_tile_obs
    
    this_lat = lat_out(i);
    this_lon = lon_out(i);
    
    tmp_sq_distance =                  ...
        (central_lon - this_lon).^2 +  ...
        (central_lat - this_lat).^2;
    
    hscale_ind{i} = find( tmp_sq_distance <= hscale^2 );
  end
  
else 
  
   hscale_ind = num2cell(ia);

end
  

% initialize output statistics
% Note: Rolf suggests to have all species as one dimension, rather than 
%       N_pol and N_angle be specified here. Then subsample specifically
%       when the files are written out.

o_data   = NaN+zeros(N_pol,N_tile_obs,N_angle,w_days);
m_data   = NaN+zeros(N_pol,N_tile_obs,N_angle,w_days);
o_data2  = NaN+zeros(N_pol,N_tile_obs,N_angle,w_days);
m_data2  = NaN+zeros(N_pol,N_tile_obs,N_angle,w_days);
N_data   = NaN+zeros(N_pol,N_tile_obs,N_angle,w_days);

data_out = NaN+zeros(N_out_fields,N_tile_obs,N_angle);

% -------------------------------------------------------------		  

% make sure t0_assim is *first* analysis time in a day

t0_assim = mod( t0_assim, dt_assim );

count = 0;

for imonth = 1:length(run_months)

  month = run_months(imonth); 
  
  for day = 1:days_in_month( 2014, month) %2014 = random non-leap year
    
    if count < w_days
      count = count + 1;
    else
      count = w_days;
    end
    
    for seconds_in_day = t0_assim:dt_assim:(86400-1)
      
      hour    = floor(seconds_in_day/3600);
      
      % check if diurnal stats are needed
      
      if exist('time_of_day_in_hours','var')
        tmp_hour = time_of_day_in_hours;
      else
        tmp_hour = hour;       % all hours of day will be included
      end
      
      if hour==tmp_hour
        
        minute  = floor( (seconds_in_day-hour*3600)/60 );
        
        seconds = seconds_in_day-hour*3600-minute*60;
        
        if (seconds~=0)
               input('something is wrong! Ctrl-c now')
        end
        
        for year = start_year(imonth):end_year(imonth)
          
          YYYYMMDD = [ num2str(year,   '%4.4d'),     ...
                       num2str(month,  '%2.2d'),     ...
                       num2str(day,    '%2.2d')  ];    
          
          HHMM     = [ num2str(hour,   '%2.2d'),     ...
                       num2str(minute, '%2.2d')  ]; 
          
          % read innov files
          
          fname = [ inpath, '/ana/ens_avg/',                  ...
                    'Y', YYYYMMDD(1:4), '/',                  ...
                    'M', YYYYMMDD(5:6), '/',                  ...
                    exp_run, '.ens_avg.ldas_ObsFcstAna.',     ...
                    YYYYMMDD, '_', HHMM, 'z.bin' ];
          
          ifp = fopen( fname, 'r', 'l' );          
          
          if (ifp > 0)           % Proceed only if file exists (e.g. irregular SMOS swaths!)
            
            fclose(ifp);
            
            [ date_time,              ...
              obs_assim,              ...
              obs_species,            ...
              obs_tilenum,            ...
              obs_lon,                ...
              obs_lat,                ...
              obs_obs,                ...
              obs_obsvar,             ...
              obs_fcst,               ...
              obs_fcstvar,            ...
              obs_ana,                ...
              obs_anavar              ...
            ] =                     ...
                read_ObsFcstAna( fname );
            
            % remove tiles where obs_fcst is no-data (note: read_ObsFcstAna() returns NaN)
            
            idx = isnan(obs_fcst);
            
            obs_assim(  idx) = [];
            obs_species(idx) = [];
            obs_tilenum(idx) = [];
            obs_lon(    idx) = [];
            obs_lat(    idx) = [];
            obs_obs(    idx) = [];
            obs_obsvar( idx) = [];
            obs_fcst(   idx) = [];
            obs_fcstvar(idx) = [];
            obs_ana(    idx) = [];
            obs_anavar( idx) = [];
            
            % extract species of interest
            
            ind = [];
            
            for this_species = species
              
              ind = find( obs_species == this_species);
              
              if (~isempty(ind))
                
                obs_tilenum_i = obs_tilenum(ind); 
                obs_obs_i     = obs_obs(ind);  
                obs_fcst_i    = obs_fcst(ind);
                obs_lon_i     = obs_lon(ind);
                obs_lat_i     = obs_lat(ind);
                
                % Check if any location receives more than 1 obs (or 1 species)
                
                tmp = sort(obs_tilenum_i);
                same_tile = find(diff(tmp)==0);
                
                if (~isempty(same_tile))
                  error('multiple obs of the same species at one location? - only last one in line is used');
                end
                
                % Organize the data in a big matrix
                
                angle = obs_param(this_species == [obs_param.species]).ang;
                pol   = obs_param(this_species == [obs_param.species]).pol;
                
                % pol intrinsically gives an index
                % now find the index for the angle
                angle_i = find(angle(1) == inc_angle);
                
                % Only writes lat-lon at exact obs locations, but with
                % hscale>0, these obs are spread outside their exact
                % location. This allows to calculate stats at lan-lons 
                % where no obs are available.
                
                %lon_out(obs_tilenum_i) = obs_lon_i;
                %lat_out(obs_tilenum_i) = obs_lat_i;
                
                % obs_lat/lon are the actual M36 lat/lons, *not* the
                % administering tiles, so the lat/lons for the obs and those
                % in the tile_coord would not be identical.
                % Still, they should be in the
                % neighbourhood, so check here if that is true.
                if (any(abs(tile_coord.com_lat(obs_tilenum_i)-obs_lat_i) > tol) || ...
                    any(abs(tile_coord.com_lon(obs_tilenum_i)-obs_lon_i) > tol) )
                  error('Something wrong with tile_lat/lon')
                end
                
                % map model tiles (e.g. all M09) to observation administering
                % tiles (could be a reduced subset of all M09)
                % --------------------------------------------------------
                obs_i = obsnum(obs_tilenum_i); 
                % --------------------------------------------------------
                
                if (hscale == 0)
                  
                  % 11 May 2015: sum the obs and fcst within each day;
                  % and across years!
                  % some obs can be found at multiple hours within a day
                  % e.g. at the poles.
                  % **sum(...,"omitnan") of NaNs** results in zero, this need to be
                  % taken care of
                  o_data( pol(1),obs_i,angle_i,count) = sum([o_data( pol(1),obs_i,angle_i,count); obs_obs_i'          ], "omitnan");
                  m_data( pol(1),obs_i,angle_i,count) = sum([m_data( pol(1),obs_i,angle_i,count); obs_fcst_i'         ], "omitnan");
                  
                  % X^2
                  o_data2(pol(1),obs_i,angle_i,count) = sum([o_data2(pol(1),obs_i,angle_i,count); obs_obs_i'.^2       ], "omitnan");
                  m_data2(pol(1),obs_i,angle_i,count) = sum([m_data2(pol(1),obs_i,angle_i,count); obs_fcst_i'.^2      ], "omitnan");
                  
                  % Sum of obs or model elements at each location
                  N_data(pol(1), obs_i,angle_i,count) = sum([N_data( pol(1),obs_i,angle_i,count); ~isnan([obs_obs_i])'], "omitnan");
                  
                else
                  
                  for i_ind = 1:length(obs_obs_i)
                    
                    % introduce a spatial effect of each observation on 
                    % neighbouring statistics (through hscale)
                    s_eff = unique(hscale_ind{obs_i(i_ind)});
                    %hscale_ind =[obs space] % 
                    
                    % Sum of X
                    o_data(pol(1),s_eff,angle_i,count) = ...
                        sum([o_data( pol(1),s_eff,angle_i,count); repmat(        obs_obs_i( i_ind),   1,length(s_eff))], "omitnan");
                    m_data(pol(1),s_eff,angle_i,count) = ...
                        sum([m_data( pol(1),s_eff,angle_i,count); repmat(        obs_fcst_i(i_ind),   1,length(s_eff))], "omitnan");
                    
                    % Sum of X^2
                    o_data2(pol(1),s_eff,angle_i,count) = ...
                        sum([o_data2(pol(1),s_eff,angle_i,count); repmat(        obs_obs_i( i_ind).^2,1,length(s_eff))], "omitnan");
                    m_data2(pol(1),s_eff,angle_i,count) = ...
                        sum([m_data2(pol(1),s_eff,angle_i,count); repmat(        obs_fcst_i(i_ind).^2,1,length(s_eff))], "omitnan");
                    
                    % Sum of obs or model elements at each location
                    N_data(pol(1),s_eff,angle_i,count) = ...
                        sum([N_data( pol(1),s_eff,angle_i,count); repmat(~isnan([obs_obs_i(i_ind)]),  1,length(s_eff))], "omitnan");
                    
                  end
                  
                end  % (hscale == 0)
                
              end  % ~isempty(ind)
              
            end  % species
            
          end  % if file present
          
        end  % loop over multiple years
        
      end  % hour == tmp_hour (time_of_day_in_hours)
      
    end  % seconds_in_day
    
    %count = count+1;
    
    if count >= w_days %wait initially until enough data is built up
      
      end_time.year  = 2014;
      end_time.month = month;
      end_time.day   = day;  
      end_time.hour  = hour;  
      end_time.min   = minute;
      end_time.sec   = seconds;
    
      start_time     = augment_date_time( -floor(w_days*(24*60*60)), end_time );
      
      % At the end of each day, collect the obs and fcst of the last
      % w_day period, and write out a statistics-file at [w_day - floor(w_day/2)]
      
      o_data(abs(o_data - nodata) <= nodata_tol) = NaN;
      m_data(abs(o_data - nodata) <= nodata_tol) = NaN;
      
      % data_out = zeros(N_out_fields,1:N_tiles,N_angle);
      
      for pol=[0 1]
        
        pp = pol*Nf;
        
        N_hscale_window         = sum(N_data(1+pol,:,:,1:w_days),                           4,"omitnan");
        
        if w_days == 95
          N_hscale_inner_window = sum(N_data(1+pol,:,:,((w_days+1)/2-15):((w_days+1)/2+15)),4,"omitnan");
        end
        
        % OBSERVATIONS
        %----------------
        % o_data is a sum over neighbouring obs above; 
        % here then take a sum over the time steps in the window
        data_out(1+pp,:,:) = sum(  o_data(   1+pol,:,:,1:w_days),4,"omitnan");
        
        % then make the average, by dividing over the sum of the number of
        % timesteps and influencing obs at each location
        data_out(1+pp,:,:) = data_out(       1+pp, :,:)./N_hscale_window;   
        
        %stdv_H = sqrt(E[X^2] - E[X]^2)
        data_out(2+pp,:,:) = sum(  o_data2(  1+pol,:,:,1:w_days),4,"omitnan");
        data_out(2+pp,:,:) =       data_out( 2+pp, :,:)./N_hscale_window;
        data_out(2+pp,:,:) = sqrt( data_out( 2+pp, :,:) - data_out(1+pp,:,:).^2);
        
        % MODEL
        %----------------
        data_out(3+pp,:,:) = sum(  m_data(   1+pol,:,:,1:w_days),4,"omitnan");
        data_out(3+pp,:,:) =       data_out( 3+pp, :,:)./N_hscale_window;            
        
        data_out(4+pp,:,:) = sum(  m_data2(  1+pol,:,:,1:w_days),4,"omitnan");
        data_out(4+pp,:,:) =       data_out( 4+pp, :,:)./N_hscale_window;
        data_out(4+pp,:,:) = sqrt( data_out( 4+pp, :,:) - data_out(3+pp,:,:).^2);
        
        data_out(5+pp,:,:) = N_hscale_window;
        
        % Toss out stats that are based on too little data
        
        data_out(  [1:5]+pp,N_hscale_window       < Ndata_min      ) = NaN;
        
        if w_days == 95
          data_out([1:5]+pp,N_hscale_inner_window < (Ndata_min/2.5)) = NaN;
        end
        
      end
      
      % Get the actual obs/model at the center point (for debugging only!!)
      
      data_out(11,:,:) = o_data(1,:,:,w_days-floor(w_days/2.0))./N_data(1,:,:,w_days-floor(w_days/2.0));
      data_out(12,:,:) = m_data(1,:,:,w_days-floor(w_days/2.0))./N_data(1,:,:,w_days-floor(w_days/2.0));
      data_out(13,:,:) = o_data(2,:,:,w_days-floor(w_days/2.0))./N_data(2,:,:,w_days-floor(w_days/2.0));
      data_out(14,:,:) = m_data(2,:,:,w_days-floor(w_days/2.0))./N_data(2,:,:,w_days-floor(w_days/2.0));
      
      % Get rid of NaN before writing a file
    
      data_out(isnan(data_out)) = nodata;
      %lon_out(isnan(lon_out))   = nodata;
      %lat_out(isnan(lat_out))   = nodata;
      
      % write output file 
      
      date_time = end_time;
      date_time = augment_date_time( -floor(w_days*(24*60*60)/2.0), date_time );
      
      % always 365 files
      
      DOY      = date_time.dofyr;
      
      if(is_leap_year(date_time.year) && DOY>=59)
        
        DOY = DOY-1;
        
        error('This code should never hit a leap year');
        
      end
    
    
      fname_out = [fname_out_base, '_DOY', num2str(DOY,'%3.3d'), '.bin'];
      
      % check whether output file exists
      
      if (exist(fname_out)==2 && overwrite) 
        
        disp(['output file exists. overwriting', fname_out])
        
      elseif (exist(fname_out)==2 && ~overwrite) 
        
        disp(['output file exists. not overwriting. returning'])
        disp(['writing ', fname_out])
        return
        
      else
        
        disp(['creating ', fname_out])
        
      end
      
      % write output for each DOY, sorted by all tiles
      
      if print_each_DOY
        
        write_seqbin_file(fname_out, lon_out, lat_out,                ...
                          inc_angle, data_out(:,:,:), int_Asc, 0,     ...  % instead of writing the version#, write Ndata_min=0
                          start_time, end_time, overwrite,            ...
                          N_out_fields, write_ind_latlon, 'scaling',  ...
                          tile_coord_tile_id)
      else
        
        % if DOY is at middle of pentad, then copy the DOY to a pentad file
        % DOY = pentad*5 - 2; ==> pentad = (DOY + 2)/5;
        
        pentad = (DOY + 2)/5;
        
        if mod((DOY + 2),5) == 0
          
          write_seqbin_file(fname_out, lon_out, lat_out,                ...
                            inc_angle, data_out(:,:,:), int_Asc, 0,     ...  
                            start_time, end_time, overwrite,            ...
                            N_out_fields, write_ind_latlon, 'scaling',  ...
                            tile_coord_tile_id)
          
          fname_out_p = [fname_out_base_p, '_p', num2str(pentad,'%2.2d'), '.bin'];
          
          copyfile(fname_out,fname_out_p);
          
        end
        
      end
      
      %clear idx_keep lon_out_write lat_out_write data_out_write tile_coord_tile_id_write
      
      % shift the window by one day and make room for the next day at the end      
      
      o_data( :,:,:,1:w_days-1) = o_data( :,:,:,2:w_days);
      m_data( :,:,:,1:w_days-1) = m_data( :,:,:,2:w_days);
      o_data2(:,:,:,1:w_days-1) = o_data2(:,:,:,2:w_days);
      m_data2(:,:,:,1:w_days-1) = m_data2(:,:,:,2:w_days);
      N_data( :,:,:,1:w_days-1) = N_data( :,:,:,2:w_days);
      
      o_data( :,:,:,w_days)  = NaN;
      m_data( :,:,:,w_days)  = NaN;
      o_data2(:,:,:,w_days)  = NaN;
      m_data2(:,:,:,w_days)  = NaN;
      N_data( :,:,:,w_days)  = NaN;
      
      data_out = NaN+0.0.*data_out;      
      
    end
    
  end        % day
end        % month


% ==================== EOF ==============================================
