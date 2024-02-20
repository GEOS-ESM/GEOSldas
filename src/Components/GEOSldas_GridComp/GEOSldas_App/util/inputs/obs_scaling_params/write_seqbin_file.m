function [] = write_seqbin_file(fname, colind, rowind,...
    av_angle_bin, data, asc_flag, ...
    version, ...
    start_time, end_time, overwrite, N_out_fields, ...
    write_ind_latlon, data_product,...
    tile_id)  %last argument is optional

% write "fortran sequential" tile tavg files (identical to LDASsa output)
%
% optional input:
%
%   overwrite = 0 -- do NOT overwrite existing files, print warning
%                    message, return
%   overwrite = 1 -- overwrite existing files, print warning message
%
% De Lannoy,  4 Oct 2010
% De Lannoy, 26 Sep 2012: added optional argument of tile_id
%           used to write scaling files, with ''latlon_id''.
% De Lannoy, 25 Oct 2012: added the processor version number,
%           inserted after the Asc_flag.
% ------------------------------------------------------------------

%N_out_fields     % 1 - Col-index, 0-based;
                  % 2 - Row-index, 0-based; 
                  %OR (for nearest neighbout)
                  % 1 - Lon;
                  % 2 - Lat; 
                  
%N_out_fields     % 1 - Tbh; 
                  % 2 - Tbv; 
  
                  % 3 - heterogeneity index Tbh
                  % 4 - heterogeneity index Tbv
                  
                  % 5 - # SMOS pixels in EASE grid pixel Tbh
                  % 6 - # SMOS pixels in EASE grid pixel Tbv
     
                  % 7 - RA Tbh
                  % 8 - RA Tbv
 
                  %=> repeated for T3 and T4 (9-16)
                 
                  %OR FOR SMUDP2:
                  
                  % 1 - SM
                  % 2 - ST
                  % 3 - opacity      
                  % 4 - Tbh; 
                  % 5 - Tbv; 
                  
                  % 6 - SM RSTD
                  % 7 - ST RSTD
                  % 8 - opac RSTD
                  
                  % 9 - stdv in SM (grid cell averaging)
                  % 10 - # small SMOS pixels inside 1 EASE grid cell;
                         
                  % 11 - omega; scattering albedo
                  % 12 - diff_albedos (om_H-om_V) 
                  % 13 - max_roughness
                  % 14 - RSTD omega
                  % 15 - RSTD diff_omega
                  % 16 - RSTD max_roughness
                  
int_precision   = 'int32';      % precision of fortran tag
float_precision = 'float32';    % precision of data in input file

% check dimensions

if size(data,1)~=N_out_fields

  error('ERROR: size of data incompatible with N_out_fields')

end

% check for presence of optional input "overwrite"

if ~exist('overwrite','var')

  overwrite = 0;        % default: do NOT overwrite existing files

end

% check if file exists

if exist(fname,'file')

  if overwrite==0

    disp(['RETURNING!!! -- NOT OVERWRITING EXISTING FILE ', fname])

    return

  else

    disp(['OVERWRITING ', fname])

  end

else
  disp(['writing ', fname])

end

% open file

ifp = fopen( fname, 'w', 'b' );

% determine number of grid cells ; further check dimensions

N_grid = size(data,2);
N_angle= 1;

if (length(size(data)) == 3)
  N_angle = size(data,3);
  data_org = data;
  if (N_angle ~= length(av_angle_bin))
    disp(['ERROR in N_angle'])
    return
  end
end

if (strcmp(write_ind_latlon,'latlon_id') && nargin == 14)
    
    if( size(tile_id,1) ~= N_grid )
        error('tile_id dimensions ??')
    end
    if ( size(tile_id,2) > 1)
        disp(['# subgridcells per gridcell: ',num2str(size(tile_id,2))]);
    end
    
end


% write all records

fortran_tag = 2*4;   % length of each record in bytes

count = fwrite( ifp, fortran_tag,    int_precision );
count = fwrite( ifp, [asc_flag version], int_precision );
count = fwrite( ifp, fortran_tag,    int_precision );

fortran_tag = 5*4;   % length of each record in bytes

count = fwrite( ifp, fortran_tag,    int_precision );
count = fwrite( ifp, [start_time.year, start_time.month, ...
    start_time.day, start_time.hour, start_time.min], int_precision );
count = fwrite( ifp, fortran_tag,    int_precision );

count = fwrite( ifp, fortran_tag,    int_precision );
count = fwrite( ifp, [end_time.year, end_time.month, ...
    end_time.day, end_time.hour, end_time.min], int_precision );
count = fwrite( ifp, fortran_tag,    int_precision );


if (~(strcmp(data_product,'scaling') && strcmp(write_ind_latlon,'latlon_id') && nargin == 14))
    fortran_tag = 2*4;   % length of each record in bytes
    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, [N_grid N_angle], int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );
else
    fortran_tag = 3*4;   % length of each record in bytes
    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, [N_grid N_angle size(tile_id,2)], int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );
end

if (N_grid >= 1)

  fortran_tag = N_angle*4;  %angles for which the output fields will be repeated below

  count = fwrite( ifp, fortran_tag,    int_precision );
  count = fwrite( ifp, squeeze(av_angle_bin(:)), float_precision );
  count = fwrite( ifp, fortran_tag,    int_precision );

  fortran_tag = N_grid*4;   % length of each record in bytes

  if (strcmp(write_ind_latlon,'ind') )
        
    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, round(colind(:)), int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision ); 
      
    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, round(rowind(:)), int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision ); 

  elseif (strcmp(write_ind_latlon,'latlon') )        
        
    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, colind(:),     float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, rowind(:),     float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

  elseif (strcmp(write_ind_latlon,'latlon_id') && nargin == 14)      

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, colind(:),     float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, rowind(:),     float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );      
             
    for i=1:size(tile_id,2)
    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, round(tile_id(:,i)),  int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );      
    end
    
  else      

    error('output-arguments do not line up')
      
  end

  fortran_tag = N_grid*4; 

  for i=1:N_out_fields
      
      for j=1:N_angle
          
          if (N_angle > 1)
              data = squeeze(data_org(:,:,j));
          end
          
          count = fwrite( ifp, fortran_tag,    int_precision );
          count = fwrite( ifp, data(i,:), float_precision );
          count = fwrite( ifp, fortran_tag,    int_precision );
          
          
      end
      
  end

else

  fortran_tag = N_angle*4;

  count = fwrite( ifp, fortran_tag,    int_precision );
  count = fwrite( ifp, squeeze(av_angle_bin(:)), float_precision );
  count = fwrite( ifp, fortran_tag,    int_precision );

  fortran_tag = 4;   % length of each record in bytes
  
  if (strcmp(write_ind_latlon,'ind') )

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0, int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0, int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

  elseif (strcmp(write_ind_latlon,'latlon') )

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0.0, float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0.0, float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

  elseif (strcmp(write_ind_latlon,'latlon_id') && nargin == 14)

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0.0, float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0.0, float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0, int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

  else

    error('output-arguments do not line up')

  end

  for i=1:N_out_fields
      
      for j=1:N_angle
          
          count = fwrite( ifp, fortran_tag,    int_precision );
          count = fwrite( ifp, -999.0, float_precision );
          count = fwrite( ifp, fortran_tag,    int_precision );
          
      end
      
  end
  
end

fclose(ifp);


