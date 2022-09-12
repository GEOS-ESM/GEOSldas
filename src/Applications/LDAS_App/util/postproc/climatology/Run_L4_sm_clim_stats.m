% Compute statistics of root-zone and profile soil moisture, which are needed to compute
% L4SM soil moisture percentile output.

clear
addpath ../../shared/matlab
%-------------------------------------------------------------------------

% experiment information 
exp_path = {'/home/qliu/smap/SMAP_Nature/SMAP_Nature_v10/'};
exp_run  = {'SMAP_Nature_v10.0'}; 

domain   = 'SMAP_EASEv2_M09_GLOBAL'; 

file_tag = {'tavg3_1d_lnr_Nt'};

% climatological period start and end year 
start_year(1:12) = [repmat(2001,12,1)]; %start year for each month,the short start year was due to CPCU data inconsistance 
end_year(1:12)   = [repmat(2021,12,1)]; %end year for each month

% only use 1 variable if running into memory issues 
field_names = {'sm_rootzone','sm_profile'};   

% convert soil moisture variables in "field_names" (if any) into
% wetness units
% (needed for L4SM because L4 ops script "prcntl.py" expects clim in 
%  wetness units, but clim is computed from Nature Run, which only 
%  outputs volumetric soil moisture) 
out_wetness = 1;

% linked to the output resolution
start_HHSS.hour = 1;
start_HHSS.min  = 30;
start_HHSS.sec  = 0;

out_freq   = 'pentad';   % pentad or monthly for now

% now define the smoothing window based on the number of years of the clim period
if end_year(1) - start_year(1) > 9
   w_out_freq = 5;  % smoothing window (number of pentads or months)
else
   w_out_freq = 11;
end

% minimum number of data requirement is based on number of years and window size
N_data_min = 2 * w_out_freq *(end_year(1)-start_year(1)+1); %per out_freq  

% ------------------

for n=1:length(exp_run)

  for ff = 1:length(field_names)

    if exist('time_of_day_in_hours','var')
    
      for j=1:length(time_of_day_in_hours)

        get_model_clim_stats( field_names{ff},                                         ...
                              exp_path{n}, exp_run{n}, domain, start_year, end_year,   ...
                              out_freq, w_out_freq, file_tag{n},                       ...
                              out_wetness,                                             ...
                              N_data_min, time_of_day_in_hours(j) )

      end
      
    else
      
      get_model_clim_stats( field_names{ff},                                           ...
                            exp_path{n}, exp_run{n}, domain,                           ...
                            start_year, end_year, start_HHSS,                          ...
                            out_freq, w_out_freq, file_tag{n},                         ...
                            out_wetness,                                               ...
                            N_data_min )
      
    end
    
  end 
end
  

% ============= EOF ====================================================


