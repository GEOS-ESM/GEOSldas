clear

% -------------------------------------------------------------------
%                     Begin user-defined inputs
% -------------------------------------------------------------------

% addpath('../../shared/matlab/');
addpath('/discover/nobackup/amfox/current_GEOSldas/GEOSldas/src/Components/GEOSldas_GridComp/GEOSldas_App/util/shared/matlab')

% Define the Open Loop experiment path, run name, domain, and output prefix

exp_path   = '/discover/nobackup/amfox/Experiments/OLv7_M36_ascat';
exp_run    = {'OLv7_M36_ascat'};
domain     = 'SMAP_EASEv2_M36_GLOBAL';
prefix_out = 'M36_zscore_stats_';

% Define the Open Loop experiment start and end dates

start_month =    4;
start_year  = 2015;
end_month   =    3;
end_year    = 2021;

% Define the species names

species_names = {'ASCAT_META_SM','ASCAT_METB_SM','ASCAT_METC_SM'};

% Define whether to combine species

combine_species_stats = 1;   % 1 to combine all species into single set of statistics

% Define the grid resolution (degrees)

grid_resolution = 0.25;

% Define moving window size over which statistics are calculated,
% and minimum number of data points required to calculate statistics

w_days    = 75;
Ndata_min =  5;

% Define the assimilation time step and initial time

dt_assim = 3*60*60;
t0_assim = 0;

% Define print intervals

print_each_DOY    = 1;
print_each_pentad = 0;
print_all_pentads = 1;

% Define output directory (takes form "domain"/stats/"out_dir")
out_dir = 'z_score_clim_quarter_degree';

% Define the months to run over, 1:12, plus a number of months required to complete the window
run_months = [1:12 1:ceil(w_days/30)];

% -------------------------------------------------------------------
%                     End user-defined inputs
% -------------------------------------------------------------------

% Calculate the earliest and latest years for each month in the experiment
earliest_year = zeros(length(run_months),1);
latest_year   = zeros(length(run_months),1);

cnt = 0;
for month = run_months
    % Initialize the earliest and latest year variables
    cnt = cnt + 1;

    % Check if the current year/month combination is earlier than the earliest
    if datenum(start_year, month, 1) < datenum(start_year, start_month, 1)
        earliest_year(cnt) = start_year+1;
    else
        earliest_year(cnt) = start_year;
    end

    % Check if the current year/month combination is later than the latest
    if datenum(end_year, month, 1) > datenum(end_year, end_month,1)
        latest_year(cnt) = end_year-1;
    else
        latest_year(cnt) = end_year;
    end
end

% assume "ldas_obsparam" file is available at 0z on first day of start_month/start_year

YYYY = num2str( start_year,  '%4.4d' );
MM   = num2str( start_month, '%2.2d' );

obs_param_fname = [exp_path, '/', exp_run{1}, '/output/', domain, '/rc_out/Y', YYYY,        ...
                   '/M', MM, '/',exp_run{1}, '.ldas_obsparam.', YYYY, MM, '01_0000z.txt'];

[N_obs_param, obs_param ] = read_obsparam(obs_param_fname);

species =[];

for i = 1:length(species_names)
    add_species = obs_param(strcmp(species_names(i),{obs_param.descr})).species;
    species = union(species,add_species);
end

if combine_species_stats
    disp('Calculating stats by combining multiple species');
end

% Calculate the climatology statistics

get_model_and_obs_clim_stats_latlon_grid( species_names, run_months, exp_path, exp_run{1}, domain, earliest_year, ...
                                          latest_year, dt_assim, t0_assim, species, combine_species_stats,        ...
                                          grid_resolution, w_days, Ndata_min, prefix_out, print_each_DOY,         ...
                                          print_each_pentad, print_all_pentads, out_dir );

% ================= EOF =========================================================================
