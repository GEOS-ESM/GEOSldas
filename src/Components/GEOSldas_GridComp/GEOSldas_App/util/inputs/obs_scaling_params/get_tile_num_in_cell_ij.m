%=========================================================================

function [N_tile_in_cell_ij, tile_num_in_cell_ij] = ...
                get_tile_num_in_cell_ij(tile_coord, tile_grid)
            
    % Initialize

    max_tile_in_cell = 10; %for EASE grids, this is really just 1

    tile_num_in_cell_ij = NaN+zeros(tile_grid.N_lon,tile_grid.N_lat,max_tile_in_cell);

    % adjust for 0-based indexing (eg., EASE grids)

    off_i = tile_grid.i_offg + (tile_grid.ind_base - 1);
    off_j = tile_grid.j_offg + (tile_grid.ind_base - 1);

    % (re-)initialize

    N_tile_in_cell_ij = zeros(tile_grid.N_lon,tile_grid.N_lat);

    for n=1:tile_coord.N_tile

       i = tile_coord.i_indg(n) - off_i;
       j = tile_coord.j_indg(n) - off_j;

       N_tile_in_cell_ij(i,j) = N_tile_in_cell_ij(i,j) + 1;

       k = N_tile_in_cell_ij(i,j);

       tile_num_in_cell_ij(i,j,k) = n;

    end 

    max_N = max(max(N_tile_in_cell_ij));

    disp(['Maximum number of tiles in tile def grid cell = ', num2str(max_N)]);

    tile_num_in_cell_ij = tile_num_in_cell_ij(:,:,1:max_N);

end

%=========================================================================