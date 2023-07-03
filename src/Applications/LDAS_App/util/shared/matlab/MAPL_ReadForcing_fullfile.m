function [tile_data,N_tile,start_time,end_time] = ...
                        MAPL_ReadForcing_fullfile(fname,nodata,nodata_tolfrac)

% Matlab version of MAPL_ReadForcing().  So far only reads complete file!
%
% Reads binary climatology or time series files (e.g., lai.dat, vegopacity.dat). 

% Q. Liu,   18 Jul 2022 
% rreichle, 29 Jul 2022
                    
% -------------------------------------------------------------------------
%
% check whether no-data variables are available on input 

if ~exist('nodata',         'var'),  nodata         = 1.e15;  end   % default: MAPL_UNDEF
if ~exist('nodata_tolfrac', 'var'),  nodata_tolfrac = 1.e-4;  end

nodata_tol = abs( nodata*nodata_tolfrac ); 

% -------------------------------------------------------------------------

int_precision   = 'int32';      % precision of fortran tag
float_precision = 'float32';    % precision of data in input file

disp(['reading from ', fname])

ifp = fopen( fname, 'r', 'l' );

% time loop (continue until reach end-of-file) 

nn=0;

while 1

   nn = nn +1;
   
   % read header

   fortran_tag = fread( ifp,  1, int_precision );
   header      = fread( ifp, 14, float_precision );
   fortran_tag = fread( ifp,  1, int_precision );
   
   % safe way to detect end-of-file
   if isempty(header)
       break
   end

   start_time(nn).year  = header( 1);
   start_time(nn).month = header( 2);
   start_time(nn).day   = header( 3);
   start_time(nn).hour  = header( 4);
   start_time(nn).min   = header( 5);
   start_time(nn).sec   = header( 6);

   end_time(  nn).year  = header( 7);
   end_time(  nn).month = header( 8);
   end_time(  nn).day   = header( 9);       
   end_time(  nn).hour  = header(10);
   end_time(  nn).min   = header(11);
   end_time(  nn).sec   = header(12);       

   N_tile = header(13);

   % read science data

   fortran_tag = fread( ifp,  1,         int_precision );
   tmp_data    = fread( ifp, [1 N_tile], float_precision );
   fortran_tag = fread( ifp,  1,         int_precision );

   tile_data(nn,:) = tmp_data;

end

% replace nodata values with NaN

tile_data( abs( single(tile_data) - single(nodata) ) < nodata_tol ) = NaN;

fclose(ifp);

% ========================= EOF =============================================
