
function [ lmc, units ] = read_smapL4SMlmc( fname, N_tile, isLDASsa);

% reichle, 26 Apr 2013

% NOTE: For large files this reader is inefficient (slow execution,
% excessive memory demand)  due to the use of a matlab structure
% array. If you need better performance, convert to reading 
% data into a regular matrix (as opposed to a structure array).
%
% GDL,     17 Feb 2014: - added units
%                       - revised fieldnames for consistency with L4_SM Product Specs Doc
% reichle, 27 May 2014: - changed wilting point output from "clsm_wpwet" to "clsm_wp"
% reichle, 17 Nov 2015: - added "veghght" output
% reichle, 28 Jul 2022 - cleaned up LDASsa/GEOSldas switch for commit into GEOSldas repo

% ----------------------------------------------------------------

if ~exist('isLDASsa','var')  isLDASsa = 0; end  % default is GEOSldas output

% for backward compatibility, back out number of parameters in file
% from file size:

% file size = N_param * (N_tile + 2) * bytes_per_datapoint

tmps = dir(fname);

N_param = tmps.bytes/((N_tile+2)*4);

if     N_param==34 | N_param==35
  
  int_records = [17 18];

else
  
  error('read_smapL4SMlmc.m: something wrong with file size or format')
  
end

disp(['read_smapL4SMlmc.m: expecting ', num2str(N_param), ' parameters in file'])

% ----------------------------------------------------------------

int_precision   = 'int32';      % precision of fortran tag
float_precision = 'float32';    % precision of data in input file

disp(['read_smapL4SMlmc.m: reading from ', fname])

if isLDASsa ~= 0
  machfmt = 'b'; % big-endian, LDASsa
else
  machfmt = 'l'; % little-endian, GEOSldas
end

ifp = fopen( fname, 'r', machfmt );

for i=1:N_param
  
  fortran_tag = fread( ifp, 1, int_precision );

  if (4*N_tile ~= fortran_tag)

    error('read_smapL4SMlmc.m: inconsistent N_tile')

  end

  if any(i==int_records)
    tmp         = fread( ifp, [1 N_tile], int_precision );
  else
    tmp         = fread( ifp, [1 N_tile], float_precision );
  end

  fortran_tag = fread( ifp, 1, int_precision );
  
  tmp_data(i,:) = tmp;
  
end

fclose(ifp);

% ---------------------------------------------------------

disp(['read_smapL4SMlmc.m: assembling structure array'])

lmc.cell_land_fraction  = tmp_data( 1,:)';  units{ 1} = '[dimensionless]';
lmc.cell_elevation      = tmp_data( 2,:)';  units{ 2} = '[m]';

lmc.clsm_dzsf           = tmp_data( 3,:)';  units{ 3} = '[m]';
lmc.clsm_dzrz      	= tmp_data( 4,:)';  units{ 4} = '[m]';
lmc.clsm_dzpr      	= tmp_data( 5,:)';  units{ 5} = '[m]';
			                  
lmc.clsm_dztsurf   	= tmp_data( 6,:)';  units{ 6} = '[m]';
			                  
lmc.clsm_dzgt1   	= tmp_data( 7,:)';  units{ 7} = '[m]';
lmc.clsm_dzgt2   	= tmp_data( 8,:)';  units{ 8} = '[m]';
lmc.clsm_dzgt3   	= tmp_data( 9,:)';  units{ 9} = '[m]';
lmc.clsm_dzgt4   	= tmp_data(10,:)';  units{10} = '[m]';
lmc.clsm_dzgt5   	= tmp_data(11,:)';  units{11} = '[m]';
lmc.clsm_dzgt6   	= tmp_data(12,:)';  units{12} = '[m]';
			                  
lmc.clsm_poros    	= tmp_data(13,:)';  units{13} = '[m3 m-3]';
lmc.clsm_wp     	= tmp_data(14,:)';  units{14} = '[m3 m-3]';
			                  
lmc.clsm_cdcr1   	= tmp_data(15,:)';  units{15} = '[kg m-2]';
lmc.clsm_cdcr2   	= tmp_data(16,:)';  units{16} = '[kg m-2]';


lmc.mwrtm_vegcls        = tmp_data(17,:)';  units{17} = '[dimensionless]';
lmc.mwrtm_soilcls      	= tmp_data(18,:)';  units{18} = '[dimensionless]';
			                  
lmc.mwrtm_sand         	= tmp_data(19,:)';  units{19} = '[dimensionless]';
lmc.mwrtm_clay          = tmp_data(20,:)';  units{20} = '[dimensionless]';
lmc.mwrtm_poros         = tmp_data(21,:)';  units{21} = '[m3 m-3]';
			                  
lmc.mwrtm_wangwt    	= tmp_data(22,:)';  units{22} = '[m3 m-3]';
lmc.mwrtm_wangwp    	= tmp_data(23,:)';  units{23} = '[m3 m-3]';
			                  
lmc.mwrtm_rghhmin   	= tmp_data(24,:)';  units{24} = '[dimensionless]';
lmc.mwrtm_rghhmax   	= tmp_data(25,:)';  units{25} = '[dimensionless]';
lmc.mwrtm_rghwmin   	= tmp_data(26,:)';  units{26} = '[m3 m-3]';
lmc.mwrtm_rghwmax   	= tmp_data(27,:)';  units{27} = '[m3 m-3]';
lmc.mwrtm_rghnrh     	= tmp_data(28,:)';  units{28} = '[dimensionless]';
lmc.mwrtm_rghnrv      	= tmp_data(29,:)';  units{29} = '[dimensionless]';
lmc.mwrtm_rghpolmix   	= tmp_data(30,:)';  units{30} = '[dimensionless]';
			                  
lmc.mwrtm_omega        	= tmp_data(31,:)';  units{31} = '[dimensionless]';
			                  
lmc.mwrtm_bh           	= tmp_data(32,:)';  units{32} = '[dimensionless]';
lmc.mwrtm_bv            = tmp_data(33,:)';  units{33} = '[dimensionless]';
lmc.mwrtm_lewt        	= tmp_data(34,:)';  units{34} = '[kg m-2]';

if N_param==35

  lmc.clsm_veghght   	= tmp_data(35,:)';  units{35} = '[m]';

end

% =========== EOF ===========================================

