
function [ tile_grid_g, tile_grid_d ] = read_tilegrids( fname, isLDASsa )

% read tile grid definitions for "global" and "domain" grids 
% from *_tilegrids.[ext] file written by LDASsa
%
% reichle, 8 July 2010
% reichle, 7 Jan  2014 - added capability to read binary "tilegrids" files
%                        file extension: ".txt" --> ASCII file
%                                        ".bin" --> binary file
% jperket,  4 Dec 2017 - added flag for LDASsa, big-endian format
% reichle, 28 Jul 2022 - cleaned up LDASsa/GEOSldas switch for commit into GEOSldas repo
%
% -------------------------------------------------------------

if ~exist('isLDASsa','var')  isLDASsa = 0; end  % default is GEOSldas output

int_precision   = 'int32';      % precision of fortran tag
float_precision = 'float32';    % precision of data in input file

if isLDASsa ~= 0
  machfmt = 'b'; % big-endian, LDASsa
else
  machfmt = 'l'; % little-endian, GEOSldas
end

% determine file name extension 

file_ext = deblank(fname);

file_ext = file_ext(end-3:end);

% ---------------------------------
%
% read file

if strcmp(file_ext,'.txt')
  
  % read ASCII file
  
  disp(['reading from ', fname])
  
  % open *_tilegrids.txt file
  
  ifp = fopen( fname, 'rt' );
  
  % read contents into cell array
  
  A=textscan(ifp,'%s');
  
  fclose(ifp);

  disp('done reading file')
  
  % --------------------------------------------------
  
  % re-assemble into long string and evaluate

  C=[];
  
  for i=1:length(A{1})
    
    C=[C, A{1}{i}];
    
  end
  
  eval(C)  % now have structures "tile_grid_g" and "tile_grid_d" 
  
elseif strcmp(file_ext,'.bin')

  % read binary file
  
  disp(['reading from ', fname])
  
  % open *_tilegrids.txt file
  
  ifp = fopen( fname, 'r', machfmt);
  
  % read contents 
  
  % first record: "global" grid (tile_grid_g)
  
  fortran_tag           = fread( ifp,  1, int_precision );
  tile_grid_g.gridtype  = fread( ifp, 40, 'uint8=>char' ); 
  tile_grid_g.ind_base  = fread( ifp,  1, int_precision );
  tile_grid_g.i_dir     = fread( ifp,  1, int_precision );    
  tile_grid_g.j_dir     = fread( ifp,  1, int_precision );    
  tile_grid_g.N_lon     = fread( ifp,  1, int_precision );    
  tile_grid_g.N_lat     = fread( ifp,  1, int_precision );     
  tile_grid_g.i_offg    = fread( ifp,  1, int_precision );   
  tile_grid_g.j_offg    = fread( ifp,  1, int_precision );  
  tile_grid_g.ll_lon    = fread( ifp,  1, float_precision );   
  tile_grid_g.ll_lat    = fread( ifp,  1, float_precision );  
  tile_grid_g.ur_lon    = fread( ifp,  1, float_precision );   
  tile_grid_g.ur_lat    = fread( ifp,  1, float_precision );  
  tile_grid_g.dlon      = fread( ifp,  1, float_precision );     
  tile_grid_g.dlat      = fread( ifp,  1, float_precision );         
  fortran_tag           = fread( ifp,  1, int_precision );

  tile_grid_g.gridtype  = deblank(tile_grid_g.gridtype');
  
  % second record: "domain" grid (tile_grid_d)
  
  fortran_tag           = fread( ifp,  1, int_precision );
  tile_grid_d.gridtype  = fread( ifp, 40, 'uint8=>char' ); 
  tile_grid_d.ind_base  = fread( ifp,  1, int_precision );
  tile_grid_d.i_dir     = fread( ifp,  1, int_precision );    
  tile_grid_d.j_dir     = fread( ifp,  1, int_precision );   
  tile_grid_d.N_lon     = fread( ifp,  1, int_precision );    
  tile_grid_d.N_lat     = fread( ifp,  1, int_precision );   
  tile_grid_d.i_offg    = fread( ifp,  1, int_precision );   
  tile_grid_d.j_offg    = fread( ifp,  1, int_precision );  
  tile_grid_d.ll_lon    = fread( ifp,  1, float_precision );   
  tile_grid_d.ll_lat    = fread( ifp,  1, float_precision );  
  tile_grid_d.ur_lon    = fread( ifp,  1, float_precision );   
  tile_grid_d.ur_lat    = fread( ifp,  1, float_precision );  
  tile_grid_d.dlon      = fread( ifp,  1, float_precision );     
  tile_grid_d.dlat      = fread( ifp,  1, float_precision );         
  fortran_tag           = fread( ifp,  1, int_precision );

  tile_grid_d.gridtype  = deblank(tile_grid_d.gridtype');
  
  % close file
  
  fclose(ifp);

  disp('done reading file')
  
else

  error('read_tilegrids.m: ERROR - unknown file extension')
  
end

% ================== EOF =====================================

