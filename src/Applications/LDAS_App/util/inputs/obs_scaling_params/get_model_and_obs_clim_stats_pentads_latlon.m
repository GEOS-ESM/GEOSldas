
function [] = get_model_and_obs_clim_stats_pentads_latlon( species_names,              ...
    run_months, exp_path, exp_run, domain, start_year, end_year,   ...
    dt_assim, t0_assim, species, combine_species_stats, obs_param, ...
    hscale, w_days, Ndata_min, prefix, ...
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
% We calcualte the bias correction factors and write on a regular 0.25 degree 
% lat/lon grid as there is no regular grid for ASCAT observations

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

disp('ASSUMING ACAT observations/undefined observation grid');

% output specs

overwrite    =  1;

if combine_species_stats
    N_species = 1;
else
    N_species = length(species);
end

Nf = 5;                  %5 fields per species
N_pentads = 73;

write_ind_latlon = 'latlon_id'; %'latlon';

% tmp_shift_lon    = 0.01;
% tmp_shift_lat    = 0.005;

% More user switches that could be moved
store_all_025_latlon = 0;
print_each_DOY = 0;
print_each_pentad = 1;
print_all_pentads = 1;

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

fname_out_base = [ outpath, '/', prefix,   ...
          num2str(min(start_year)),'_doy',num2str(D(1)),'_',...
          num2str(max(end_year)),  '_doy',num2str(D(2)),...
	      '_W_', num2str(w_days),'d_Nmin_', num2str(Ndata_min)];
      
fname_out_base_p = [ outpath, '/', prefix,   ...
          num2str(min(start_year)),'_p',num2str(P(1)),'_',...
          num2str(max(end_year)),  '_p',num2str(P(2)),...
	      '_W_', num2str(round(w_days/5)),'p_Nmin_', num2str(Ndata_min)];

%==============================================================
% Some clunky code to maintain backwards compatibility with adding orbit
% tag

% Initialize counters for cells ending in "_A" and cells ending in "_D"
a_count = 0;
d_count = 0;

% Loop through each cell in the array
for i = 1:numel(species_names)
    % Check if the text in the cell ends with either "_A" or "_D"
    if endsWith(species_names{i}, '_A')
        % If it ends with "_A", increment the "_A" counter
        a_count = a_count + 1;
    elseif endsWith(species_names{i}, '_D')
        % If it ends with "_D", increment the "_D" counter
        d_count = d_count + 1;
    end
    
    if startsWith(species_names{i}, 'SMAP')
        inc_angle = [40.0];
    else
        inc_angle = [-999.9];
    end
    
end

% Determine the output based on the values of the "_A" and "_D" counters
if a_count == numel(species_names)
    % Both cells end in "_A"
    disp('All species are "_A"');
    Orbit_tag = '_A'; 
     int_Asc = 1;
elseif d_count == numel(species_names)
    % Both cells end in "_D"
    disp('All species are "_D"');
    Orbit_tag = '_D'; 
    int_Asc = 2;
elseif a_count > 0 && d_count > 0
    % There is a mix of "_A" and "_D"
    disp('Species have a mix of "_A" and "_D"');
    Orbit_tag = '_AD'; 
    int_Asc = 3;
else
    % Neither cell ends in "_A" or "_D"
    disp('Neither cell ends in "_A" or "_D"');
    Orbit_tag = '_NoOrbits';
    int_Asc = 4;
end

fname_out_base   = [fname_out_base, Orbit_tag];
fname_out_base_p = [fname_out_base_p, Orbit_tag];
   
if exist( 'time_of_day_in_hours', 'var')
  
  fname_out_base   = [fname_out_base, '_', num2str(time_of_day_in_hours,'%2.2d'), 'z'];
  fname_out_base_p = [fname_out_base_p, '_', num2str(time_of_day_in_hours,'%2.2d'), 'z'];
  
end

%======================================================

% -------------------------------------------------------------
% Define our 1/4 degree lat/lon grid
n_lon = 1440;
n_lat = 720;
ll_lon = -180;
ll_lat = -90;
d_lon = 0.25;
d_lat = 0.25;
ll_lons = linspace(-180, 179.75, 1440);
ll_lats  = linspace(-90, 89.75, 720);
i_lon = (1:1440);
j_lat =  (1:720);

grid_idx = zeros(1036800,5);
grid_idx(:,1) = (1:1036800);
cnt = 0;
for i = 1:n_lon
    for j = 1:n_lat
        cnt = cnt+1;
        grid_idx(cnt,2) = i_lon(i);
        grid_idx(cnt,3) = j_lat(j);
        grid_idx(cnt,4) = ll_lons(i);
        grid_idx(cnt,5) = ll_lats(j);        
    end
end
        
% -------------------------------------------------------------

obsnum = grid_idx(:,1); 
i_out = grid_idx(:,2);
j_out = grid_idx(:,3);
lon_out    = grid_idx(:,4);
lat_out    = grid_idx(:,5);
N_gridcells = length(grid_idx);

% Not sure about this as aren't centering lat/lons
if hscale>0

    for i=1:N_gridcells

        this_lat = lat_out(i);
        this_lon = lon_out(i);

        tmp_sq_distance =              ...
        (lon_out - this_lon).^2 +  ...
        (lat_out - this_lat).^2;

        hscale_ind{i} = find( tmp_sq_distance <= hscale^2 );
    end

else 

    hscale_ind = num2cell(obsnum);

end
  

% initialize output statistics
% Note: Rolf suggests to have all species as one dimension, rather than 
%       N_pol and N_angle be specified here. Then subsample specifically
%       when the files are written out.

o_data     = NaN+zeros(N_species, N_gridcells, w_days);
m_data     = NaN+zeros(N_species, N_gridcells, w_days);
o_data2    = NaN+zeros(N_species, N_gridcells, w_days);
m_data2    = NaN+zeros(N_species, N_gridcells, w_days);
N_data     = NaN+zeros(N_species, N_gridcells, w_days);

data_out   = NaN+zeros(N_species, Nf, N_gridcells, N_pentads);
data2D = NaN+zeros(Nf, N_gridcells);

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
                      exp_run, '.ens_avg.ldas_ObsFcstAna.',        ...
                      YYYYMMDD, '_', HHMM, 'z.bin' ];

           ifp = fopen( fname, 'r', 'l' );          

           if (ifp > 0)           %Proceed only if file exists (e.g. irregular SMOS swaths!)

           fclose(ifp);

           [date_time,           ...
            obs_assim,         ...
            obs_species,      ...
            obs_tilenum,      ...
            obs_lon,              ...
            obs_lat,               ...
            obs_obs,             ...
            obs_obsvar,        ...
            obs_fcst,             ...
            obs_fcstvar,       ...
            obs_ana,             ...
            obs_anavar        ...
                 ] =                    ...
                read_ObsFcstAna( fname );
            
          % remove tiles when there is no obs_fcst (obs_fcst == 0 in innov output when
          % missing)
          
          idx = find(obs_fcst == 0);
          obs_assim(idx) = [];
          obs_species(idx) =  [];
          obs_tilenum(idx) =[];
          obs_lon(idx) =[];
          obs_lat(idx) = [];
          obs_obs(idx) = [];
          obs_obsvar(idx) = [];
          obs_fcst(idx) = [];
          obs_fcstvar(idx) = [];
          obs_ana(idx) = [];
          obs_anavar(idx) = [];

          % extract species of interest
          
          ind = [];

          %for this_species = species
          for scnt = 1:N_species
              
              if combine_species_stats % We are  combining stats
                  ind = find(ismember(obs_species, species));
                  this_species = species(scnt); % Only first species in list. But only used in determining angle and pol for species, which shouldn't vary between species being combined
              else
                  ind = find( obs_species == this_species);
              end

             if (~isempty(ind))

                obs_tilenum_i = obs_tilenum(ind); 
                obs_obs_i     = obs_obs(ind);  
                obs_fcst_i    = obs_fcst(ind);
                obs_lon_i     = obs_lon(ind);
                obs_lat_i     = obs_lat(ind);

                % Check if any location receives more than 1 obs (or 1 species)

                tmp = sort(obs_tilenum_i);
                same_tile = find(diff(tmp)==0, 1);

                if (~isempty(same_tile) && combine_species_stats==0)
                    error('multiple obs of the same species at one location? - only last one in line is used');
                end

                % Organize the data in a big matrix

                pol   = obs_param(this_species == [obs_param.species]).pol;
              
                % Only writes lat-lon at exact obs locations, but with
                % hscale>0, these obs are spread outside their exact
                % location. This allows to calculate stats at lan-lons 
                % where no obs are available.

                %lon_out(obs_tilenum_i) = obs_lon_i;
                %lat_out(obs_tilenum_i) = obs_lat_i;

                %obs_lat/lon are the actual M36 lat/lons, *not* the
                %administering tiles, so the lat/lons for the obs and those
                %in the tile_coord would not be identical.
                %Still, they should be in the
                %neighbourhood, so check here if that is true.
%                 if (any(abs(tile_coord.com_lat(obs_tilenum_i)-obs_lat_i) > tol) || ...
%                     any(abs(tile_coord.com_lon(obs_tilenum_i)-obs_lon_i) > tol) )
%                     error('Something wrong with tile_lat/lon')
%                 end

                %map model tiles (e.g. all M09) to observation administering
                %tiles (could be a reduced subset of all M09)
                %--------------------------------------------------------
%                 obs_i = obsnum(obs_tilenum_i); 
                %--------------------------------------------------------
                
                % Put obs lat/lon on our grid and figure out obsnum/grid
                % index
                
                i_idx = floor((obs_lon_i - ll_lon)/d_lon) + 1;
                j_idx = floor((obs_lat_i - ll_lat)/d_lat) + 1;
                [~, obs_idx] = ismember([i_idx, j_idx], [i_out, j_out], 'rows');
                
                obs_i = obsnum(obs_idx); 
                        
                if (hscale == 0)

                        %11 May 2015: sum the obs and fcst within each day;
                        %and across years!
                        %some obs can be found at multiple hours within a day
                        %e.g. at the poles.
                        %**nansum of NaN's** result in zero, this need to be
                        %taken care of
                        o_data(scnt,obs_i,count) = nansum([o_data(scnt,obs_i,count); obs_obs_i' ]);
                        m_data(scnt,obs_i,count) = nansum([m_data(scnt,obs_i,count); obs_fcst_i']);
                            
                        %X^2
                        o_data2(scnt,obs_i,count) = nansum([o_data2(scnt,obs_i,count); obs_obs_i'.^2 ]);
                        m_data2(scnt,obs_i,count) = nansum([m_data2(scnt,obs_i,count); obs_fcst_i'.^2]);
                            
                        %Sum of obs or model elements at each location
                        N_data(scnt,obs_i,count) = nansum([N_data(scnt,obs_i,count); ~isnan([obs_obs_i])']);

                else

                    for i_ind = 1:length(obs_obs_i)

                        %introduce a spatial effect of each observation on 
                        %neighbouring statistics (through hscale)
                        s_eff = unique(hscale_ind{obs_i(i_ind)});
                                      %hscale_ind =[obs space] % 

                        %Sum of X
                        o_data(scnt,obs_i,count) = ...
                            nansum([o_data(scnt,obs_i,count); repmat(obs_obs_i(i_ind),1,length(s_eff))]);
                        m_data(scnt,obs_i,count) = ...
                            nansum([m_data(pol(1),s_eff,angle_i,count); repmat(obs_fcst_i(i_ind),1,length(s_eff))]);

                        %Sum of X^2
                        o_data2(scnt,obs_i,count) = ...
                            nansum([o_data2(scnt,obs_i,count); repmat(obs_obs_i(i_ind).^2,1,length(s_eff))]);
                        m_data2(scnt,obs_i,count) = ...
                            nansum([m_data2(scnt,obs_i,count); repmat(obs_fcst_i(i_ind).^2,1,length(s_eff))]);

                        %Sum of obs or model elements at each location
                        N_data(scnt,obs_i,count) = ...
                            nansum([N_data(scnt,obs_i,count);  repmat(~isnan([obs_obs_i(i_ind)]),1,length(s_eff)) ]);

                    end

                end

            end

          end

          end  % if file present

      end % loop over multiple years

  end  % time_of_day_in_hours

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
  
  for i = 1:N_species
      
      N_hscale_window = nansum(N_data(i,:,1:w_days),3);

      if w_days == 95
         N_hscale_inner_window = nansum(N_data(i,:,((w_days+1)/2-15):((w_days+1)/2+15)),3);
      end
      
      % OBSERVATIONS
      %----------------
      %o_data is a sum over neighbouring obs above; 
      %here then take a sum over the time steps in the window
      data2D(1,:) = nansum(o_data(i,:,1:w_days),3);

      %then make the average, by dividing over the sum of the number of
      %timesteps and influencing obs at each location
      data2D(1,:) = data2D(1,:)./N_hscale_window;   

      %stdv_H = sqrt(E[X^2] - E[X]^2)
      data2D(2,:) = nansum(o_data2(i,:,1:w_days),3);
      data2D(2,:) = data2D(2,:)./N_hscale_window;
      data2D(2,:) = sqrt( data2D(2,:) - data2D(1,:).^2);

      % MODEL
      %----------------
      data2D(3,:) = nansum(m_data(i,:,1:w_days),3);
      data2D(3,:) = data2D(3,:)./N_hscale_window;            

      data2D(4,:) = nansum(m_data2(i,:,1:w_days),3);
      data2D(4,:) = data2D(4,:)./N_hscale_window;
      data2D(4,:) = sqrt( data2D(4,:) - data2D(3,:).^2);

      data2D(5,:) = N_hscale_window;

      % Toss out stats that are based on too little data

      data2D([1:Nf],N_hscale_window<Ndata_min) = NaN;
        

    % Get rid of NaN before writing a file

  data_out(isnan(data_out)) = nodata;
  
   
 startidx = strfind(fname_out_base, 'z_score_clim//');
 endidx = startidx + length('z_score_clim//');
  
 if combine_species_stats % We are  combining stats
     fname_out_base_s = [fname_out_base(1:startidx-1) 'z_score_clim/combined_', fname_out_base(endidx:end)];
 else
     fname_out_base_s = [fname_out_base(1:startidx-1) 'z_score_clim/', char(species_names(i)),'_', fname_out_base(endidx:end)];
 end
  
 
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

  
  fname_out = [fname_out_base_s, '_DOY', num2str(DOY,'%3.3d'), '.nc4'];

  % compress data before writing in file. 
  
  %idx_keep = find(any(abs(data_out -nodata) > nodata_tol,1));
  %lon_out_write = lon_out(idx_keep);
  %lat_out_write = lat_out(idx_keep);
  %data_out_write = data_out(:,idx_keep);
  %tile_coord_tile_id_write = tile_coord_tile_id(idx_keep);
  
  
  % write output for each DOY, sorted by all tile
  if print_each_DOY
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
      write_netcdf_file_2D_grid(fname_out, i_out, j_out, lon_out, lat_out, ...
                inc_angle, data2D, int_Asc, pentad, ...  %instead of writing the version#, write pentad
                start_time, end_time, overwrite, ...
                Nf, write_ind_latlon, 'scaling',...
                obsnum)
  else
      % if DOY is at middle of pentad, then copy the DOY to a pentad file
      % DOY = pentad*5 - 2; ==> pentad = (DOY + 2)/5;
      pentad = (DOY + 2)/5;
      if mod((DOY + 2),5) == 0
          data_out(i,:,:,pentad) = data2D;
          start_time_p(pentad) = start_time;
          end_time_p(pentad) = end_time;
          if print_each_pentad
              fname_out = [fname_out_base_p, '_p', num2str(pentad,'%2.2d'), '.nc4'];
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
              write_netcdf_file_2D_grid(fname_out, i_out, j_out, lon_out, lat_out, ...
                  inc_angle, data2D, int_Asc, pentad, ...  
                  start_time, end_time, overwrite, ...
                  Nf, write_ind_latlon, 'scaling',...
                  obsnum)
          end
      end
  
  end
  
  %clear idx_keep lon_out_write lat_out_write data_out_write tile_coord_tile_id_write
  
  % shift the window by one day and make room for the next day at the end      

  o_data(:,:,1:w_days-1)  = o_data(:,:,2:w_days);
  m_data(:,:,1:w_days-1)  = m_data(:,:,2:w_days);
  o_data2(:,:,1:w_days-1) = o_data2(:,:,2:w_days);
  m_data2(:,:,1:w_days-1) = m_data2(:,:,2:w_days);
  N_data(:,:,1:w_days-1)  = N_data(:,:,2:w_days);

  o_data(:,:,w_days)  = NaN;
  m_data(:,:,w_days)  = NaN;
  o_data2(:,:,w_days) = NaN;
  m_data2(:,:,w_days) = NaN;
  N_data(:,:,w_days)  = NaN;

  data2D = NaN+0.0.*data2D;      
  
  end
  
  end
  
end        % day
end        % month

if print_all_pentads
    for i = 1:N_species
        data_o = squeeze(data_out(i,:,:,:));
        
         if combine_species_stats % We are  combining stats
             fname_out = [fname_out_base(1:startidx-1) 'z_score_clim/combined_all_pentads_', fname_out_base(endidx:end),'.nc4'];
         else
             fname_out = [fname_out_base(1:startidx-1) 'z_score_clim/', char(species_names(i)),'_all_pentads_', fname_out_base(endidx:end),'.nc4'];
         end

        write_netcdf_file_2D_grid(fname_out, i_out, j_out, lon_out, lat_out, ...
                  inc_angle, data_o, int_Asc, [1:73], ...  
                  start_time_p, end_time_p, overwrite, ...
                  Nf, write_ind_latlon, 'scaling',...
                  obsnum)
    end
end


% ==================== EOF ==============================================
