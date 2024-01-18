
function [ aup, units ] = read_smapL4SMaup( fname, N_tile, isLDASsa );

% reichle, 26 Apr 2013
% reichle,  5 Feb 2014 - added tb_[h/v]_obs_time_sec
% reichle, 21 Mar 2015 - changed units of soil moisture output from wetness
%                         [dimensionless] to volumetric soil moisture [m3/m3]
% reichle, 28 Jul 2022 - cleaned up LDASsa/GEOSldas switch for commit into GEOSldas repo

% NOTE: For large files this reader is inefficient (slow execution,
% excessive memory demand)  due to the use of a matlab structure
% array. If better performance is needed, convert to reading 
% data into a regular matrix (as opposed to a structure array).

if ~exist('isLDASsa','var')  is_LDASsa = 0; end  % default is GEOSldas output

N_param     = 31;           % number of records

dbl_records = [1 2];        % double precision records
  
int_records = [3 4 5 6];    % integer records


% ----------------------------------------------------------------

int_precision   = 'int32';      % precision of fortran tag and integer data in input file
float_precision = 'float32';    % precision of real data in input file
dbl_precision   = 'float64';    % precision of real*8 data in input file

disp(['read_smapL4SMaup.m: reading from ', fname])
  
if isLDASsa ~= 0
  machfmt = 'b'; % big-endian, LDASsa
else
  machfmt = 'l'; % little-endian, GEOSldas
end

ifp = fopen( fname, 'r', machfmt );

tmp_data = NaN*ones(N_param,N_tile);

for i=1:N_param
  
  fortran_tag = fread( ifp, 1, int_precision );

  if any(i==dbl_records)

    N_bytes = 8;

  else

    N_bytes = 4;
  
  end

  if (N_bytes*N_tile ~= fortran_tag)

    error('read_smapL4SMaup.m: inconsistent N_tile')

  end
    
  if     any(i==int_records)
    tmp         = fread( ifp, [1 N_tile], int_precision );
  elseif any(i==dbl_records)
    tmp         = fread( ifp, [1 N_tile], dbl_precision );
  else
    tmp         = fread( ifp, [1 N_tile], float_precision );
  end

  fortran_tag = fread( ifp, 1, int_precision );
  
  tmp_data(i,:) = tmp;

end

fclose(ifp);

% ---------------------------------------------------------

disp(['read_smapL4SMaup.m: assembling structure array'])

aup.tb_h_obs_time_sec                   = tmp_data( 1,:)';  units{ 1} = '[s]';
aup.tb_v_obs_time_sec                   = tmp_data( 2,:)';  units{ 2} = '[s]';
aup.tb_h_resolution_flag                = tmp_data( 3,:)';  units{ 3} = '[dimensionless]';
aup.tb_v_resolution_flag                = tmp_data( 4,:)';  units{ 4} = '[dimensionless]';
aup.tb_h_orbit_flag                     = tmp_data( 5,:)';  units{ 5} = '[dimensionless]';
aup.tb_v_orbit_flag                     = tmp_data( 6,:)';  units{ 6} = '[dimensionless]';
aup.tb_h_obs                            = tmp_data( 7,:)';  units{ 7} = '[K]'; 
aup.tb_v_obs                            = tmp_data( 8,:)';  units{ 8} = '[K]'; 

aup.tb_h_obs_assim                      = tmp_data( 9,:)';  units{ 9} = '[K]';   
aup.tb_v_obs_assim                      = tmp_data(10,:)';  units{10} = '[K]';   
aup.tb_h_obs_errstd                     = tmp_data(11,:)';  units{11} = '[K]'; 
aup.tb_v_obs_errstd                     = tmp_data(12,:)';  units{12} = '[K]'; 

aup.tb_h_forecast                       = tmp_data(13,:)';  units{13} = '[K]'; 
aup.tb_v_forecast                       = tmp_data(14,:)';  units{14} = '[K]'; 
aup.tb_h_forecast_ensstd                = tmp_data(15,:)';  units{15} = '[K]'; 
aup.tb_v_forecast_ensstd                = tmp_data(16,:)';  units{16} = '[K]'; 

aup.sm_surface_forecast                 = tmp_data(17,:)';  units{17} = '[m3 m-3]'; 
aup.sm_rootzone_forecast                = tmp_data(18,:)';  units{18} = '[m3 m-3]';       
aup.sm_profile_forecast                 = tmp_data(19,:)';  units{19} = '[m3 m-3]'; 
aup.surface_temp_forecast               = tmp_data(20,:)';  units{20} = '[K]'; 
aup.soil_temp_layer1_forecast           = tmp_data(21,:)';  units{21} = '[K]'; 

aup.sm_surface_analysis                 = tmp_data(22,:)';  units{22} = '[m3 m-3]';      
aup.sm_rootzone_analysis                = tmp_data(23,:)';  units{23} = '[m3 m-3]'; 
aup.sm_profile_analysis                 = tmp_data(24,:)';  units{24} = '[m3 m-3]'; 
aup.surface_temp_analysis               = tmp_data(25,:)';  units{25} = '[K]';       
aup.soil_temp_layer1_analysis           = tmp_data(26,:)';  units{26} = '[K]';  

aup.sm_surface_analysis_ensstd          = tmp_data(27,:)';  units{27} = '[m3 m-3]';       
aup.sm_rootzone_analysis_ensstd         = tmp_data(28,:)';  units{28} = '[m3 m-3]'; 
aup.sm_profile_analysis_ensstd          = tmp_data(29,:)';  units{29} = '[m3 m-3]'; 
aup.surface_temp_analysis_ensstd        = tmp_data(30,:)';  units{30} = '[K]'; 
aup.soil_temp_layer1_analysis_ensstd    = tmp_data(31,:)';  units{31} = '[K]'; 
  
% =========== EOF ===========================================

































