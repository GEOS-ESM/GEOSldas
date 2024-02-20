function [data, tile_id, col_ind, row_ind, fieldno] = ...
    read_seqbin_clim_pctl_file(fname, N_field, N_stat, read_ind_latlon ,file_type, ifp)

% Gabrielle De Lannoy, 27 Feb 2014
%
% ------------------------------------------------------------------

int_precision   = 'int32';      % precision of fortran tag
float_precision = 'float32';    % precision of data in input file

disp(['reading from ', fname])

if ~exist('ifp','var')
  ifp = fopen( fname, 'r', 'b' );
end

data    = NaN;
tile_id = NaN;
col_ind = NaN;
row_ind = NaN;
fieldno = NaN;

if (strcmp(file_type,'cli'))

N_f = N_field;

fortran_tag  = fread( ifp, 1,     int_precision );
tmp_data_int = fread( ifp, [1 2], int_precision );
fortran_tag  = fread( ifp, 1,     int_precision );

N_grid       = tmp_data_int(1);
N_field      = tmp_data_int(2); 

if N_field~=N_f
  error('N_fields in clim file?')
end

disp(['N_grid and N_field:',num2str(N_grid),' and ' ,num2str(N_field)]);

end

% read all records

if (N_grid > 1)

  if (strcmp(file_type,'cli'))

    if (strcmp(read_ind_latlon,'ind'))

        fortran_tag = fread( ifp, 1,          int_precision );
        col_ind     = fread( ifp, [1 N_grid], int_precision );
        fortran_tag = fread( ifp, 1,          int_precision );

        fortran_tag = fread( ifp, 1,          int_precision );
        row_ind     = fread( ifp, [1 N_grid], int_precision );
        fortran_tag = fread( ifp, 1,          int_precision );

    elseif (strcmp(read_ind_latlon,'latlon'))

        fortran_tag = fread( ifp, 1,          int_precision );
        col_ind     = fread( ifp, [1 N_grid], float_precision );
        fortran_tag = fread( ifp, 1,          int_precision );

        fortran_tag = fread( ifp, 1,          int_precision );
        row_ind     = fread( ifp, [1 N_grid], float_precision );
        fortran_tag = fread( ifp, 1,          int_precision );

    elseif (strcmp(read_ind_latlon,'latlon_id') )

        fortran_tag = fread( ifp, 1,          int_precision );
        col_ind     = fread( ifp, [1 N_grid], float_precision );
        fortran_tag = fread( ifp, 1,          int_precision );

        fortran_tag = fread( ifp, 1,          int_precision );
        row_ind     = fread( ifp, [1 N_grid], float_precision );
        fortran_tag = fread( ifp, 1,          int_precision );

        fortran_tag = fread( ifp, 1,          int_precision );
        tile_id     = fread( ifp, [1 N_grid], int_precision );  
        fortran_tag = fread( ifp, 1,          int_precision );       

    else

        error('not sure how the file looks like, based on the combination of input-specs')

    end

  end

  data = NaN*ones(N_field,N_grid,N_stat);
  
  for j=1:N_stat

     for i=1:N_field 

       if (j == 5  && strcmp(file_type,'cli')) 
         
         fortran_tag = fread( ifp, 1,          int_precision );
         tmp_data    = fread( ifp, [1 N_grid], int_precision );
         fortran_tag = fread( ifp, 1,          int_precision ); 
       
       else        

         fortran_tag = fread( ifp, 1,          int_precision );
         tmp_data    = fread( ifp, [1 N_grid], float_precision );
         fortran_tag = fread( ifp, 1,          int_precision );
       
       end
       
       data(i,1:N_grid,j) = tmp_data(1:N_grid);
       
    end

  end

  data = squeeze(data);

else
    
  data      = NaN*ones(N_field,1);
  col_ind   = NaN;
  row_ind   = NaN;
  fieldno   = NaN;
    
end

fclose(ifp);
  
% ======= EOF ==========================================================
