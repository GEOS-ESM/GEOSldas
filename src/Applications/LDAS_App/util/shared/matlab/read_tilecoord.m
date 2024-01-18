
function [tile_coord ] = read_tilecoord( fname, bin2txt, isLDASsa )

% read tile coordinates from *_tilecoord.[ext] file written by LDASsa
%
% reichle, 29 Jun 2005
% GDL,     22 Jun 2010 - changed i_atm/j_atm/frac_atm to i_indg/j_indg/frac_cell
% reichle, 31 May 2011 - accomodate new field "elev" (elevation)
% reichle,  7 Jan 2014 - added capability to read binary "tilecoord" files
%                         and to convert a binary file into a txt file
%                        file extension: ".txt" --> ASCII file
%                                        ".bin" --> binary file
%                        ASCII option maintains backward compatibility
%
% jperket,  1 Dec 2017 - added flag for LDASsa, big-endian format
% reichle, 28 Jul 2022 - cleaned up LDASsa/GEOSldas switch for commit into GEOSldas repo
%
% -------------------------------------------------------------

if ~exist('isLDASsa','var')  isLDASsa = 0; end  % default is GEOSldas output

int_precision   = 'int32';      % precision of fortran tag
float_precision = 'float32';    % precision of data in input file

% deal with "optional" bin2txt argument

if ~exist('bin2txt','var') 

  bin2txt = 0;

end

if isLDASsa ~= 0
  machfmt = 'b'; % big-endian, LDASsa
else
  machfmt = 'l'; % little-endian, GEOSldas
end

% determine file name extension 

file_ext = deblank(fname);

file_ext = file_ext(end-3:end);

if strcmp(file_ext,'.txt')

  is_binary = 0; 

  if bin2txt
   
    error('read_tilecoord.m: ERROR -- bin2txt conversion ', ...
          'requires input file name for bin file');

  end

elseif strcmp(file_ext,'.bin')

  is_binary = 1;

else

  error('read_tilecoord.m: ERROR - unknown file extension')

end


% ---------------------------------
%
% read file

disp(['reading from ', fname])

if is_binary

  % open *_tilecoord.bin file

  ifp = fopen( fname, 'r', machfmt);

  fortran_tag          = fread( ifp,  1, int_precision );
  tile_coord.N_tile    = fread( ifp,  1, int_precision );
  fortran_tag          = fread( ifp,  1, int_precision );

  Nt = tile_coord.N_tile;
  
  fortran_tag          = fread( ifp,  1, int_precision );
  tile_coord.tile_id   = fread( ifp, Nt, int_precision );
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.typ       = fread( ifp, Nt, int_precision );
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.pfaf      = fread( ifp, Nt, int_precision );      
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.com_lon   = fread( ifp, Nt, float_precision );   
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.com_lat   = fread( ifp, Nt, float_precision );   
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.min_lon   = fread( ifp, Nt, float_precision );   
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.max_lon   = fread( ifp, Nt, float_precision );   
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.min_lat   = fread( ifp, Nt, float_precision );    
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.max_lat   = fread( ifp, Nt, float_precision );   
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.i_indg    = fread( ifp, Nt, int_precision );     
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.j_indg    = fread( ifp, Nt, int_precision );     
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.frac_cell = fread( ifp, Nt, float_precision );  
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.frac_pfaf = fread( ifp, Nt, float_precision ); 
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.area      = fread( ifp, Nt, float_precision );
  fortran_tag          = fread( ifp,  1, int_precision );
  
  fortran_tag          = fread( ifp,  1, int_precision );  
  tile_coord.elev      = fread( ifp, Nt, float_precision );
  fortran_tag          = fread( ifp,  1, int_precision );


  % if requested, convert to ASCII (txt) file
  
  if bin2txt 

    fname_out = deblank(fname);

    fname_out = [fname_out(1:end-4), '.txt'];

    % open *_tilecoord.txt file

    disp(['writing to ', fname])

    ofp = fopen( fname_out, 'wt' );

    fprintf( ofp, '%i\n', tile_coord.N_tile  );

    for ii=1:tile_coord.N_tile

      fprintf( ofp,['%8i%10i%9i',                           ...
                    '%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f', ... 
		    '%5i%5i%13.6f%13.6f%13.4f%13.4f\n'],    ... 
              [tile_coord.tile_id(ii),     ...
               tile_coord.typ(ii),         ...
               tile_coord.pfaf(ii),        ...      
               tile_coord.com_lon(ii),     ...   
               tile_coord.com_lat(ii),     ...   
               tile_coord.min_lon(ii),     ...   
               tile_coord.max_lon(ii),     ...   
               tile_coord.min_lat(ii),     ...    
               tile_coord.max_lat(ii),     ...   
               tile_coord.i_indg(ii),      ...     
               tile_coord.j_indg(ii),      ...     
               tile_coord.frac_cell(ii),   ...  
               tile_coord.frac_pfaf(ii),   ... 
               tile_coord.area(ii),        ...      
	       tile_coord.elev(ii)       ]);

    end

    fclose(ofp);

    disp('done writing file')
    
  end

else

  % open *_tilecoord.txt file

  ifp = fopen( fname, 'rt' );

  tmpdata = fscanf( ifp, '%f' );

  % --------------------------------------------------

  % process data

  tile_coord.N_tile        = tmpdata(1);

  N_cols                   = (length(tmpdata)-1)/tile_coord.N_tile; 

  tmpdata = reshape(tmpdata(2:end), [N_cols, tile_coord.N_tile])';

  tile_coord.tile_id       = tmpdata(:, 1);
  tile_coord.typ           = tmpdata(:, 2);
  tile_coord.pfaf          = tmpdata(:, 3);      
  tile_coord.com_lon       = tmpdata(:, 4);   
  tile_coord.com_lat       = tmpdata(:, 5);   
  tile_coord.min_lon       = tmpdata(:, 6);   
  tile_coord.max_lon       = tmpdata(:, 7);   
  tile_coord.min_lat       = tmpdata(:, 8);   
  tile_coord.max_lat       = tmpdata(:, 9);   
  tile_coord.i_indg        = tmpdata(:,10);     
  tile_coord.j_indg        = tmpdata(:,11);     
  tile_coord.frac_cell     = tmpdata(:,12);  
  tile_coord.frac_pfaf     = tmpdata(:,13); 
  tile_coord.area          = tmpdata(:,14);      

  if N_cols==15

     tile_coord.elev       = tmpdata(:,15);		
   
  end

end

% close file

fclose(ifp);

disp('done reading file')

% =========== EOF ========================================


