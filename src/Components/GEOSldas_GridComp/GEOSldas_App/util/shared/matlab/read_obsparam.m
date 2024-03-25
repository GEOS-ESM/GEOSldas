function [ N_obs_param, obs_param ] = read_obs_param( fname )
 
% Get observation parameters
% Format as in module enkf_types, subroutine write_obs_param
%
% Gabrielle De Lannoy  - 26 Oct 2011
%
%  1 Dec 2011 - reichle: minor modifications and check-in to CVS
%
%  8 Jun 2017 - reichle: added "flistpath" and "flistname" 
%
% ------------------------------------------------------------------     

fid = fopen(fname);

disp(['Reading ',fname]);

N_obs_param = fscanf(fid, '%d ', 1);

for i=1:N_obs_param

    obs_param(i).descr           = fscanf(fid, '%s ', 1);
    obs_param(i).species         = fscanf(fid, '%f ', 1);    
    obs_param(i).orbit           = fscanf(fid, '%f ', 1);  %1=A, 2=D
    obs_param(i).pol             = fscanf(fid, '%f ', 1);  %1=H, 2=V
    obs_param(i).N_ang           = fscanf(fid, '%f ', 1);
    
    obs_param(i).ang             = fscanf(fid, '%f ', obs_param(i).N_ang);

    obs_param(i).freq            = fscanf(fid, '%f ', 1);
    obs_param(i).FOV             = fscanf(fid, '%f ', 1);
    obs_param(i).FOV_units       = fscanf(fid, '%s ', 1);
    obs_param(i).assim           = fscanf(fid, '%s ', 1);
    obs_param(i).scale           = fscanf(fid, '%s ', 1);
    obs_param(i).getinnov        = fscanf(fid, '%s ', 1);
    obs_param(i).RTM_ID          = fscanf(fid, '%f ', 1);
    obs_param(i).bias_Npar       = fscanf(fid, '%f ', 1);
    obs_param(i).bias_trel       = fscanf(fid, '%f ', 1);
    obs_param(i).bias_tcut       = fscanf(fid, '%f ', 1);
    obs_param(i).nodata          = fscanf(fid, '%f ', 1);    
    obs_param(i).varname         = fscanf(fid, '%s ', 1);
    obs_param(i).units           = fscanf(fid, '%s ', 1);
    obs_param(i).path            = fscanf(fid, '%s ', 1);
    obs_param(i).name            = fscanf(fid, '%s ', 1);
    obs_param(i).maskpath        = fscanf(fid, '%s ', 1);
    obs_param(i).maskname        = fscanf(fid, '%s ', 1);
    obs_param(i).scalepath       = fscanf(fid, '%s ', 1);
    obs_param(i).scalename       = fscanf(fid, '%s ', 1);                     
    obs_param(i).flistpath       = fscanf(fid, '%s ', 1);
    obs_param(i).flistname       = fscanf(fid, '%s ', 1);                     
    obs_param(i).errstd          = fscanf(fid, '%f ', 1);
    obs_param(i).std_normal_max  = fscanf(fid, '%f ', 1);
    obs_param(i).zeromean        = fscanf(fid, '%s ', 1);
    obs_param(i).coarsen_pert    = fscanf(fid, '%s ', 1);
    obs_param(i).xcorr           = fscanf(fid, '%f ', 1);
    obs_param(i).ycorr           = fscanf(fid, '%f ', 1);
    obs_param(i).adapt           = fscanf(fid, '%f ', 1);
        
    % remove leading and trailing quotes from strings
    
    obs_param(i).descr           = obs_param(i).descr(    2:end-1);
    obs_param(i).FOV_units       = obs_param(i).FOV_units(2:end-1);
    obs_param(i).varname         = obs_param(i).varname(  2:end-1);
    obs_param(i).units           = obs_param(i).units(    2:end-1);
    obs_param(i).path            = obs_param(i).path(     2:end-1);
    obs_param(i).name            = obs_param(i).name(     2:end-1);
    obs_param(i).maskpath        = obs_param(i).maskpath( 2:end-1);
    obs_param(i).maskname        = obs_param(i).maskname( 2:end-1);
    obs_param(i).scalepath       = obs_param(i).scalepath(2:end-1);
    obs_param(i).scalename       = obs_param(i).scalename(2:end-1);
    obs_param(i).flistpath       = obs_param(i).flistpath(2:end-1);
    obs_param(i).flistname       = obs_param(i).flistname(2:end-1);
    
end

fclose(fid);

disp(['Done reading obs_param for ',num2str(N_obs_param),' species']);

% =========================== EOF ====================================
