function [date_time,              ...
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
 	 ] =                      ...
    read_ObsFcstAna( fname, isLDASsa )

%
% read_ObsFcstAna.m can be used to read "ObsFcstAna" files that 
%  contain Observations and observation-space model forecasts and 
%  analysis data
%
% data format:
%  see f90 subroutine output_ObsFcstAna() in module clsm_ensupd_enkf_update
%
% reichle,  4 Oct 2011
%
% ------------------------------------------------------------------

int_precision     = 'int32';      % precision of fortran tag
float_precision   = 'float32';    % precision of data in input file
logical_precision = 'int32';      % precision of data in input file

% initialize outputs in case file does not exist or is empty

nodata = -9999;

date_time   = struct('year',   nodata, ...
		     'month',  nodata, ...
		     'day',    nodata, ...
		     'hour',   nodata, ...
		     'min',    nodata, ...
		     'sec',    nodata, ...
		     'dofyr',  nodata, ...
		     'pentad', nodata );

obs_assim               = [];
obs_species             = [];
obs_tilenum             = [];
obs_lon                 = [];
obs_lat                 = [];
obs_obs                 = [];
obs_obsvar              = [];
obs_fcst                = [];
obs_fcstvar             = [];
obs_ana                 = [];
obs_anavar              = [];
      
if exist('isLDASsa','var') && isLDASsa == 1
  machfmt = 'b'; % big-endian, LDASsa
else
  machfmt = 'l'; % little-endian, GEOSldas
end

% read file if it exists

if exist(fname)==2
  
  disp(['reading from ', fname  ])
  
  ifp = fopen( fname, 'r', machfmt );       
  
  % read N_obs and time stamp entry
  
  fortran_tag  = fread( ifp, 1, int_precision );
  N_obs        = fread( ifp, 1, int_precision );
  year         = fread( ifp, 1, int_precision );
  month        = fread( ifp, 1, int_precision );
  day          = fread( ifp, 1, int_precision );
  hour         = fread( ifp, 1, int_precision );
  minute       = fread( ifp, 1, int_precision );
  second       = fread( ifp, 1, int_precision );
  dofyr        = fread( ifp, 1, int_precision );
  pentad       = fread( ifp, 1, int_precision );
  fortran_tag  = fread( ifp, 1, int_precision );
  
  date_time.year   = year;
  date_time.month  = month;
  date_time.day    = day;
  date_time.hour   = hour;
  date_time.min    = minute;
  date_time.sec    = second;
  date_time.dofyr  = dofyr;
  date_time.pentad = pentad;
  
  % read observation assim flag
  
  fortran_tag  = fread( ifp, 1, int_precision );
  tmp_data     = fread( ifp, [N_obs 1], logical_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	
  
  obs_assim = zeros( N_obs, 1);
  obs_assim( tmp_data~= 0 ) = 1;

  % read species information
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_species  = fread( ifp, [N_obs 1], int_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	
    
  % read tile number information
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_tilenum  = fread( ifp, [N_obs 1], int_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	

  % read longitude
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_lon      = fread( ifp, [N_obs 1], float_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	

  % read latitude
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_lat      = fread( ifp, [N_obs 1], float_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	
  

  % read observation value
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_obs      = fread( ifp, [N_obs 1], float_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	
    
  % read observation variance
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_obsvar   = fread( ifp, [N_obs 1], float_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	


  % read observation-space model forecast value
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_fcst     = fread( ifp, [N_obs 1], float_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	
    
  % read observation-space model forecast variance
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_fcstvar  = fread( ifp, [N_obs 1], float_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	


  % read observation-space analysis value
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_ana      = fread( ifp, [N_obs 1], float_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	
    
  % read observation-space analysis variance
  
  fortran_tag  = fread( ifp, 1, int_precision );
  obs_anavar   = fread( ifp, [N_obs 1], float_precision );
  fortran_tag  = fread( ifp, 1, int_precision );	
  
  
  % no-data check 
  %  - single ensemble member integrations yield obs_obsvar==nodata)  
  %  - in some cases obs_fcst (a.k.a. Obs_pred) is no-data-value, 
  %     eg. SMOS Tb when snow is present)

  obs_obsvar(  obs_obsvar    == nodata ) = NaN;	  

  obs_fcst(    obs_fcst      == nodata ) = NaN;	  
  obs_fcstvar( obs_fcstvar   == nodata ) = NaN;	  

  obs_ana(     obs_ana       == nodata ) = NaN;	  
  obs_anavar(  obs_anavar    == nodata ) = NaN;	  

  
  % close file  
  
  fclose(ifp);
  
else            % if exist(fname)==2
  
  disp(['file does not exist: ', fname])
  
end 


% ======= EOF ==========================================================
