
function [obs_descr, read_obs_subr, N_obs_read, obs_fname] = read_obslog(fname)

% read LDASsa obs log files
%
% reichle, 2 Jan 2014
%
% -------------------------------------------------------------

% read file

disp(['reading from ', fname])

fid = fopen( fname );

data = textscan(fid, '%s%s%d%s',              ...
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

else
  
  error('read_obslog(): ERROR reading data')
  
end

% extract data into output variables

obs_descr     = data{1};
read_obs_subr = data{2};
N_obs_read    = data{3};
obs_fname     = data{4};
  

% ================== EOF =====================================
