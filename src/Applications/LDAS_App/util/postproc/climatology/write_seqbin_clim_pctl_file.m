function [] = write_seqbin_clim_pctl_file(fname, colind, rowind,...
    data, fieldno, N_stat,...
    overwrite, ...
    write_ind_latlon, file_type, tile_id)  %last argument is optional

% write "fortran binary sequential" tile files with climatology info 
% or percentile output
%
% optional input:
%
%   overwrite = 0 -- do NOT overwrite existing files, print warning
%                    message, return
%   overwrite = 1 -- overwrite existing files, print warning message
%
% De Lannoy, 27 Feb 2014: adopted from write_seqbin_file.m
% ------------------------------------------------------------------
                  
nodata_val      = -9999.0;    
                  
int_precision   = 'int32';      % precision of fortran tag
float_precision = 'float32';    % precision of data in input file

% check dimensions

if size(data,1)~=length(fieldno)

  error('ERROR: size of data incompatible with N_fields')

end

N_grid  = size(data,2);

N_field = length(fieldno);

if (length(size(data)) == 3)

  N_stat_tmp = size(data,3);

  data_org = data;
  
  if (N_stat_tmp ~= N_stat)
    disp(['ERROR in N_stat ',num2str(N_stat_tmp),' vs ',num2str(N_stat)])
    return
  end
  
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


if (strcmp(write_ind_latlon,'latlon_id') && exist('tile_id','var'))
    
  if ( size(tile_id,1) ~= N_grid )
    error('tile_id dimensions ??')
  end
  if ( size(tile_id,2) > 1)
    disp(['# subgridcells per gridcell: ',num2str(size(tile_id,2))]);
  end
    
end


% write all records

if (strcmp(file_type,'cli'))

  % dimensions

  fortran_tag = 2*4;   % length of each record in bytes

  count = fwrite( ifp, fortran_tag,      int_precision );
  count = fwrite( ifp, [N_grid N_field], int_precision );
  count = fwrite( ifp, fortran_tag,      int_precision );

end

if (N_grid >= 1)

  if (strcmp(file_type,'cli') )
    
    fortran_tag = N_grid*4;   % length of each record in bytes

    if (strcmp(write_ind_latlon,'ind') )

      count = fwrite( ifp, fortran_tag,         int_precision );
      count = fwrite( ifp, round(colind(:)),    int_precision );
      count = fwrite( ifp, fortran_tag,         int_precision ); 

      count = fwrite( ifp, fortran_tag,         int_precision );
      count = fwrite( ifp, round(rowind(:)),    int_precision );
      count = fwrite( ifp, fortran_tag,         int_precision ); 

    elseif (strcmp(write_ind_latlon,'latlon') )        

      count = fwrite( ifp, fortran_tag,         int_precision );
      count = fwrite( ifp, colind(:),           float_precision );
      count = fwrite( ifp, fortran_tag,         int_precision );

      count = fwrite( ifp, fortran_tag,         int_precision );
      count = fwrite( ifp, rowind(:),           float_precision );
      count = fwrite( ifp, fortran_tag,         int_precision );

    elseif (strcmp(write_ind_latlon,'latlon_id') )      

      count = fwrite( ifp, fortran_tag,         int_precision );
      count = fwrite( ifp, colind(:),           float_precision );
      count = fwrite( ifp, fortran_tag,         int_precision );

      count = fwrite( ifp, fortran_tag,         int_precision );
      count = fwrite( ifp, rowind(:),           float_precision );
      count = fwrite( ifp, fortran_tag,         int_precision );      

      for i=1:size(tile_id,2)
      count = fwrite( ifp, fortran_tag,         int_precision );
      count = fwrite( ifp, round(tile_id(:,i)), int_precision );
      count = fwrite( ifp, fortran_tag,         int_precision );      
      end

    else      

      error('output-arguments do not line up')

    end

  end
  
  fortran_tag = N_grid*4; 

  for j=1:N_stat

    for i=1:N_field
    
      if (N_stat > 1)
        data = squeeze(data_org(:,:,j));
      end
    
      if ( j == 5 && strcmp(file_type,'cli'))
        count = fwrite( ifp, fortran_tag,      int_precision );
        count = fwrite( ifp, round(data(i,:)), int_precision );
        count = fwrite( ifp, fortran_tag,      int_precision );
      else            
        count = fwrite( ifp, fortran_tag,      int_precision );
        count = fwrite( ifp, data(i,:),        float_precision );
        count = fwrite( ifp, fortran_tag,      int_precision );
      end

    end

  end

else

  if (strcmp(file_type,'cli'))
      
  fortran_tag = 4;   % length of each record in bytes
  
  if (strcmp(write_ind_latlon,'ind') )

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0,              int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0,              int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

  elseif (strcmp(write_ind_latlon,'latlon') )

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0.0,            float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0.0,            float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

  elseif (strcmp(write_ind_latlon,'latlon_id'))

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0.0,            float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0.0,            float_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

    count = fwrite( ifp, fortran_tag,    int_precision );
    count = fwrite( ifp, 0,              int_precision );
    count = fwrite( ifp, fortran_tag,    int_precision );

  else

    error('output-arguments do not line up')

  end

  end
  
  for j=1:N_stat

    for i=1:N_field 

      if ( j == 5 && strcmp(file_type,'cli'))
        count = fwrite( ifp, fortran_tag,    int_precision );
        count = fwrite( ifp, 0,              int_precision );
        count = fwrite( ifp, fortran_tag,    int_precision );
      else           
        count = fwrite( ifp, fortran_tag,    int_precision );
        count = fwrite( ifp, nodata_val,    float_precision );
        count = fwrite( ifp, fortran_tag,    int_precision );
      end

    end

  end
  
end

fclose(ifp);

%=========================EOF====================================
