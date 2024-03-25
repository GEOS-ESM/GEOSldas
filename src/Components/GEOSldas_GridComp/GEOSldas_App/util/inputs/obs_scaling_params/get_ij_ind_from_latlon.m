%=========================================================================

function [i_ind,j_ind] = get_ij_ind_from_latlon( tile_grid, lat, lon)

    if (strcmp(tile_grid.gridtype,'EASEv2_M36'))
        %row, col
        [j_indg,i_indg] = ...
            EASEv2_latlon2ind(lat,lon,'M36',1);
    elseif (strcmp(tile_grid.gridtype,'EASEv2_M09'))
        %row, col
        [j_indg,i_indg] = ...
            EASEv2_latlon2ind(lat,lon,'M09',1);
    else
        error('not ready for this grid');
    end

    % convert to index into array defined by tile_grid_d

    i_ind = i_indg - tile_grid.i_offg - (tile_grid.ind_base - 1);
    j_ind = j_indg - tile_grid.j_offg - (tile_grid.ind_base - 1);

end

%=========================================================================
