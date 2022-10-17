
function [ana_time, obs_descr, read_obs_subr, N_obs_read, obs_fname] = read_obslog(fname)

% read LDASsa obs log files
%
% reichle,  2 Jan 2014
% reichle, 11 Feb 2021 - corrected for post-launch file format
%
% -------------------------------------------------------------

% read file

disp(['reading from ', fname])

fid = fopen( fname );

data = textscan(fid, '%s%s%s%d%s',              ...
                'Delimiter', ',',             ...
                'HeaderLines', 3   );

fclose(fid);

disp('done reading file')

% make sure that final line contained 'EOF', remove from data

if strcmp(data{1}{end},'EOF')
  
  data{1}(end) = [];
  data{2}(end) = [];
  data{3}(end) = [];
  data{4}(end) = [];
  data{5}(end) = [];

else
  
  error('read_obslog(): ERROR reading data')
  
end

% extract data into output variables


ana_time      = data{1};
obs_descr     = data{2};
read_obs_subr = data{3};
N_obs_read    = data{4};
obs_fname     = data{5};
  

% ================== EOF =====================================
