
function [ cat_param, cat_param_units ] = read_catparam( fname, N_tile, isLDASsa );

% reichle,  2 Jun 2006
% reichle, 16 Jul 2010 - added vegcls lookup table
% reichle, 28 Oct 2010 - added soilcls*
%                      - changed cat_param structure from "vector of
%                        structures" to "structure of vectors"
% reichle,  1 Apr 2015 - added new soil parameter fields (file_format==3)
%                      - added cat_param_units
% reichle, 28 Jul 2022 - cleaned up LDASsa/GEOSldas switch for commit into GEOSldas repo

% NOTE: For large files this reader is inefficient (slow execution,
% excessive memory demand)  due to the use of a matlab structure
% array. If better performance is required, convert to reading 
% data into a regular matrix (as opposed to a structure array).
%
%
% parameter "vegcls" is land cover type:
%
%  vegcls = 1:  broadleaf evergreen trees
%  vegcls = 2:  broadleaf deciduous trees
%  vegcls = 3:  needleleaf trees
%  vegcls = 4:  grassland
%  vegcls = 5:  broadleaf shrubs
%  vegcls = 6:  dwarf trees
%  vegcls = 7:  bare soil
%  vegcls = 8:  desert soil
%
%
% parameters "soilcls30" and "soilcls100" are 0-30 cm and 0-100 cm soil class:
%
% Two similar but different look-up tables were used for the MERRA and 
% Fortuna versions of the Catchment model.
%
% The first look-up table is from http://www.iges.org/gswp2/.  [Note that
% the GSWP-2 documentation refers to these values as "Cosby" (even though
% they are different from the Cosby et al 1984 values, see below.]
% The GSWP-2 values were used by Sarith in the Richards' equation
% pre-processing steps and are used in subroutine catchment() (via Sarith's
% Catchment model parameter files):
%
%   Soil   Soil             B     Porosity Wilting      Psis   *Surface* 
%   Class  Type                   (v/v)    Point(v/v)   (m)    Ks (m/s)
%
%    1  :  Sand             3.30  0.373	  0.089	    -0.05	0.0285
%    2  :  Loamy Sand       3.80  0.386   0.132     -0.07	0.0204
%    3  :  Sandy Loam       4.34  0.419   0.205     -0.16	0.0097
%    4  :  Silt Loam        5.25  0.476   0.355     -0.65	0.0027
%    5  :  Silt                  	
%    6  :  Loam             5.96  0.437   0.339     -0.24	0.0054
%    7  :  Sandy Clay Loam  7.32  0.412   0.379     -0.12	0.0073
%    8  :  Silty Clay Loam  8.41  0.478   0.521     -0.63	0.0017
%    9  :  Clay Loam        8.34  0.447   0.472     -0.28	0.0032
%   10  :  Sandy Clay       9.70  0.415   0.480     -0.12	0.0050
%   11  :  Silty Clay      10.78  0.478   0.598     -0.58	0.0012
%   12  :  Clay            12.93  0.450   0.613     -0.27	0.0015
%
% The second look-up table is from Cosby et al. WRR 1984.  These values
% were used by Randy in some of the Catchment model paramter pre-processing 
% steps (other than the Richards' equation solver):
%
%   Soil   Soil             B     Porosity Wilting      Psis   *Surface* 
%   Class  Type                   (v/v)    Point(v/v)   (m)    Ks (m/s)
%   
%    1  :  Sand             2.79  0.339    0.0218     -0.0692  0.0542931
%    2  :  Loamy Sand       4.26  0.421    0.0599     -0.0363  0.0163963
%    3  :  Sandy Loam       4.74  0.434    0.1002     -0.1413  0.00609179
%    4  :  Silt Loam        5.33  0.476    0.1772     -0.7586  0.00327148
%    5  :  Silt             5.33  0.476    0.1772     -0.7586  0.00327148
%    6  :  Loam             5.25  0.439    0.1393     -0.3548  0.00393319
%    7  :  Sandy Clay Loam  6.77  0.404    0.1438     -0.1349  0.00518495
%    8  :  Silty Clay Loam  8.72  0.464    0.2477     -0.6166  0.00236998
%    9  :  Clay Loam        8.17  0.465    0.2144     -0.2630  0.00284934
%   10  :  Sandy Clay      10.73  0.406    0.2053     -0.0977  0.00840901
%   11  :  Silty Clay      10.39  0.468    0.2597     -0.3236  0.00156583
%   12  :  Clay            11.55  0.468    0.2240     -0.0389  0.00113434
%
%
% In both cases, the *surface* Ks that is needed as input to subroutine 
% catchment() was extrapolated from the Ks provided by GSWP-2 or Cosby et al
% using the vertical conductivity decay factor "gnu" (which might be 
% inconsistent with "gnu" used elsewhere).
%
%
% Starting in late 2014, revised soil parameters can be used in LDASsa.
% For details see De Lannoy et al., 2014, doi:10.1002/2014MS000330.
%
% ------------------------------------------------------------------

if ~exist('isLDASsa','var')  isLDASsa = 0; end  % default is GEOSldas output

% for backward compatibility, back out number of parameters in file
% from file size:

% file size = N_param * (N_tile + 2) * bytes_per_datapoint

tmps = dir(fname);

if isLDASsa ~= 0
  machfmt = 'b'; % big-endian, LDASsa
else
  machfmt = 'l'; % little-endian, GEOSldas
end

N_param = tmps.bytes/((N_tile+2)*4);

if     N_param==40
  
  file_format = 1;

  if isLDASsa ~= 0
     int_columns = 18;             % vegcls
  else
     int_columns = [];             % GEOSldas files contain only real*4 numbers 
  end

elseif N_param==42 | N_param==51 | N_param==52
  
  file_format = 2;
 
  if isLDASsa ~= 0 
     int_columns = [ 18 19 20 ];   % vegcls, soilcls30, soilcls100
  else
     int_columns = [];             % GEOSldas files contain only real*4 numbers
  end
  
else
  
  error('read_catparam.m: something wrong with file size or format')
  
end

disp(['read_catparam.m: expecting ', num2str(N_param), ...
      ' parameters in file with file_format ', num2str(file_format)])

% ----------------------------------------------------------------

int_precision   = 'int32';      % precision of fortran tag
float_precision = 'float32';    % precision of data in input file

disp(['read_catparam.m: reading from ', fname])

ifp = fopen( fname, 'r', machfmt );

for i=1:N_param
  
  fortran_tag = fread( ifp, 1, int_precision );
  if any(i==int_columns)
    tmp         = fread( ifp, [1 N_tile], int_precision );
  else
    tmp         = fread( ifp, [1 N_tile], float_precision );
  end
  fortran_tag = fread( ifp, 1, int_precision );
  
  tmp_data(i,:) = tmp;
  
end

disp(['read_catparam.m: assembling structure array'])

switch file_format
  
 case {1}
  
  cat_param.dpth       = tmp_data( 1,:)';    cat_param_units.dpth       = '[mm]';     
                                                                                      
  cat_param.dzsf       = tmp_data( 2,:)';    cat_param_units.dzsf       = '[mm]';     
  cat_param.dzrz       = tmp_data( 3,:)';    cat_param_units.dzrz       = '[mm]';     
  cat_param.dzpr       = tmp_data( 4,:)';    cat_param_units.dzpr       = '[mm]';     
                                                                                      
  cat_param.dzgt(:,1)  = tmp_data( 5,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,2)  = tmp_data( 6,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,3)  = tmp_data( 7,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,4)  = tmp_data( 8,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,5)  = tmp_data( 9,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,6)  = tmp_data(10,:)';    cat_param_units.dzgt       = '[m]';      
                                                                                 	 
  cat_param.poros      = tmp_data(11,:)';    cat_param_units.poros      = '[m3 m-3]'; 
  cat_param.cond       = tmp_data(12,:)';    cat_param_units.cond       = '[m s-1]';  
  cat_param.psis       = tmp_data(13,:)';    cat_param_units.psis       = '[m H2O]';  
  cat_param.bee        = tmp_data(14,:)';    cat_param_units.bee        = '[-]';      
                                                                                 	 
  cat_param.wpwet      = tmp_data(15,:)';    cat_param_units.wpwet      = '[-]';      
                                                                                 	 
  cat_param.gnu        = tmp_data(16,:)';    cat_param_units.gnu        = '[m-1]';    
                                                                                      
  cat_param.vgwmax     = tmp_data(17,:)';    cat_param_units.vgwmax     = '[kg m-2]'; 
                                                                                      
  cat_param.vegcls     = tmp_data(18,:)';    cat_param_units.vegcls     = '[-]';      
                                             
  cat_param.bf1        = tmp_data(19,:)';    cat_param_units.bf1        = '[kg m-4]'; 
  cat_param.bf2        = tmp_data(20,:)';    cat_param_units.bf2        = '[m]';                                               
  cat_param.bf3        = tmp_data(21,:)';    cat_param_units.bf3        = '[log(m)]'; 
  cat_param.cdcr1      = tmp_data(22,:)';    cat_param_units.cdcr1      = '[kg m-2]'; 
  cat_param.cdcr2      = tmp_data(23,:)';    cat_param_units.cdcr2      = '[kg m-2]'; 
  cat_param.ars1       = tmp_data(24,:)';    cat_param_units.ars1       = '[m2 kg-1]';
  cat_param.ars2       = tmp_data(25,:)';    cat_param_units.ars2       = '[m2 kg-1]';
  cat_param.ars3       = tmp_data(26,:)';    cat_param_units.ars3       = '[m4 kg-2]';
  cat_param.ara1       = tmp_data(27,:)';    cat_param_units.ara1       = '[m2 kg-1]';
  cat_param.ara2       = tmp_data(28,:)';    cat_param_units.ara2       = '[-]';      
  cat_param.ara3       = tmp_data(29,:)';    cat_param_units.ara3       = '[m2 kg-1]';
  cat_param.ara4       = tmp_data(30,:)';    cat_param_units.ara4       = '[-]';      
  cat_param.arw1       = tmp_data(31,:)';    cat_param_units.arw1       = '[m2 kg-1]';
  cat_param.arw2       = tmp_data(32,:)';    cat_param_units.arw2       = '[m2 kg-1]';
  cat_param.arw3       = tmp_data(33,:)';    cat_param_units.arw3       = '[m4 kg-2]';
  cat_param.arw4       = tmp_data(34,:)';    cat_param_units.arw4       = '[-]';      
  cat_param.tsa1       = tmp_data(35,:)';    cat_param_units.tsa1       = '[-]';      
  cat_param.tsa2       = tmp_data(36,:)';    cat_param_units.tsa2       = '[-]';      
  cat_param.tsb1       = tmp_data(37,:)';    cat_param_units.tsb1       = '[-]';      
  cat_param.tsb2       = tmp_data(38,:)';    cat_param_units.tsb2       = '[-]';      
  cat_param.atau       = tmp_data(39,:)';    cat_param_units.atau       = '[-]';      
  cat_param.btau       = tmp_data(40,:)';    cat_param_units.btau       = '[-]';      
                                               
 case {2}

  cat_param.dpth       = tmp_data( 1,:)';    cat_param_units.dpth       = '[mm]';         
                                                                                      
  cat_param.dzsf       = tmp_data( 2,:)';    cat_param_units.dzsf       = '[mm]';     
  cat_param.dzrz       = tmp_data( 3,:)';    cat_param_units.dzrz       = '[mm]';     
  cat_param.dzpr       = tmp_data( 4,:)';    cat_param_units.dzpr       = '[mm]';     
                                                                                      
  cat_param.dzgt(:,1)  = tmp_data( 5,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,2)  = tmp_data( 6,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,3)  = tmp_data( 7,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,4)  = tmp_data( 8,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,5)  = tmp_data( 9,:)';    cat_param_units.dzgt       = '[m]';      
  cat_param.dzgt(:,6)  = tmp_data(10,:)';    cat_param_units.dzgt       = '[m]';      
                                                                                 	 
  cat_param.poros      = tmp_data(11,:)';    cat_param_units.poros      = '[m3 m-3]'; 
  cat_param.cond       = tmp_data(12,:)';    cat_param_units.cond       = '[m s-1]';  
  cat_param.psis       = tmp_data(13,:)';    cat_param_units.psis       = '[m H2O]';  
  cat_param.bee        = tmp_data(14,:)';    cat_param_units.bee        = '[-]';      
                                                                                 	 
  cat_param.wpwet      = tmp_data(15,:)';    cat_param_units.wpwet      = '[-]';      
                                                                                 	 
  cat_param.gnu        = tmp_data(16,:)';    cat_param_units.gnu        = '[m-1]';    
                                                                                      
  cat_param.vgwmax     = tmp_data(17,:)';    cat_param_units.vgwmax     = '[kg m-2]'; 
                                                                                      
  cat_param.vegcls     = tmp_data(18,:)';    cat_param_units.vegcls     = '[-]';      
  cat_param.soilcls30  = tmp_data(19,:)';    cat_param_units.soilcls30  = '[-]';      
  cat_param.soilcls100 = tmp_data(20,:)';    cat_param_units.soilcls100 = '[-]';      
                                                                                      
  cat_param.bf1        = tmp_data(21,:)';    cat_param_units.bf1        = '[kg m-4]'; 
  cat_param.bf2        = tmp_data(22,:)';    cat_param_units.bf2        = '[m]';      
  cat_param.bf3        = tmp_data(23,:)';    cat_param_units.bf3        = '[log(m)]'; 
  cat_param.cdcr1      = tmp_data(24,:)';    cat_param_units.cdcr1      = '[kg m-2]'; 
  cat_param.cdcr2      = tmp_data(25,:)';    cat_param_units.cdcr2      = '[kg m-2]'; 
  cat_param.ars1       = tmp_data(26,:)';    cat_param_units.ars1       = '[m2 kg-1]';
  cat_param.ars2       = tmp_data(27,:)';    cat_param_units.ars2       = '[m2 kg-1]';
  cat_param.ars3       = tmp_data(28,:)';    cat_param_units.ars3       = '[m4 kg-2]';
  cat_param.ara1       = tmp_data(29,:)';    cat_param_units.ara1       = '[m2 kg-1]';
  cat_param.ara2       = tmp_data(30,:)';    cat_param_units.ara2       = '[-]';      
  cat_param.ara3       = tmp_data(31,:)';    cat_param_units.ara3       = '[m2 kg-1]';
  cat_param.ara4       = tmp_data(32,:)';    cat_param_units.ara4       = '[-]';      
  cat_param.arw1       = tmp_data(33,:)';    cat_param_units.arw1       = '[m2 kg-1]';
  cat_param.arw2       = tmp_data(34,:)';    cat_param_units.arw2       = '[m2 kg-1]';
  cat_param.arw3       = tmp_data(35,:)';    cat_param_units.arw3       = '[m4 kg-2]';
  cat_param.arw4       = tmp_data(36,:)';    cat_param_units.arw4       = '[-]';      
  cat_param.tsa1       = tmp_data(37,:)';    cat_param_units.tsa1       = '[-]';      
  cat_param.tsa2       = tmp_data(38,:)';    cat_param_units.tsa2       = '[-]';      
  cat_param.tsb1       = tmp_data(39,:)';    cat_param_units.tsb1       = '[-]';      
  cat_param.tsb2       = tmp_data(40,:)';    cat_param_units.tsb2       = '[-]';      
  cat_param.atau       = tmp_data(41,:)';    cat_param_units.atau       = '[-]';      
  cat_param.btau       = tmp_data(42,:)';    cat_param_units.btau       = '[-]';      
                                                                                      
  if N_param==51 | N_param==52
                                                                                      
    cat_param.gravel30 = tmp_data(43,:)';    cat_param_units.gravel30 = '[%vol]';   
    cat_param.orgC30   = tmp_data(44,:)';    cat_param_units.orgC30   = '[%weight]';
    cat_param.orgC     = tmp_data(45,:)';    cat_param_units.orgC     = '[%weight]';
    cat_param.sand30   = tmp_data(46,:)';    cat_param_units.sand30   = '[%weight]';
    cat_param.clay30   = tmp_data(47,:)';    cat_param_units.clay30   = '[%weight]';
    cat_param.sand     = tmp_data(48,:)';    cat_param_units.sand     = '[%weight]';
    cat_param.clay     = tmp_data(49,:)';    cat_param_units.clay     = '[%weight]';
    cat_param.wpwet30  = tmp_data(50,:)';    cat_param_units.wpwet30  = '[-]';      
    cat_param.poros30  = tmp_data(51,:)';    cat_param_units.poros30  = '[m3 m-3]'; 
    
  end

  if N_param==52
                                                                                      
    cat_param.veghght  = tmp_data(52,:)';    cat_param_units.veghght  = '[m]';   

  end
  
 otherwise
  
  error('read_catparam.m: something wrong with file size or format')  
  
end

fclose(ifp);


% =========== EOF ===========================================

