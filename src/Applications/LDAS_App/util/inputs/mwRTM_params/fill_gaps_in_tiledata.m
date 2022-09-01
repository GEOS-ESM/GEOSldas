function [tile_data_filled] = fill_gaps_in_tiledata(tile_coord, tile_grid_g, tile_data, N_cells, iscube )

% Fill missing values in tile-space data with the mean value (excl. NaNs) of surrounding grid cells.
%
% N_cells is number of grid cells averaged in each linear dimension; e.g., N_cells=5 averages across
%  a 5-by-5 neighborhood  (-2, -1, 0, 1, 2 in each direction).
%
% Uses tile2grid() --> Should work for EASE[_v2] and lat/lon tile spaces.  
%                      Probably needs work for cube-sphere tile space!!! 
%
% iscube:  required input to alert user that fn is not ready for data in cube-sphere tile space
%
% tile_data[_filled] = N_fields-by-N_tile
%
% Q. Liu,  19 Jul 2022 
% reichle, 29 Jul 2022 - minor clean-up and generalization
%
% -----------------------------------------------------------------------------------------------

if ~exist('iscube'), error('Must specify if data is in cube-sphere tile space.'), end

if iscube, error('Function not ready for data in cube-sphere tile space.'), end

tc = tile_coord;
tg = tile_grid_g;

%if     strcmp(EASEv2_grid,'M09')
%    N_cells = 5;
%elseif strcmp(EASEv2_grid,'M36')
%    N_cells = 3;
%else
%    error('input grid invalid, use only M09 or M36')
%end

if size(tile_data,2) ~= tc.N_tile
    error('N_tile incorrect, input data size should be [N_fields,N_tile].')
end

grid_data        = tile2grid(tile_data, tc, tg);

N_f              = size(tile_data,1);
N_lon            = size(grid_data,1);
N_lat            = size(grid_data,2);

tile_data_filled = tile_data;

Nshift           = floor(N_cells/2);

for ff = 1:N_f

    grid  = grid_data(:,:,ff);
    d_sum = zeros(size(grid));
    N_sum = zeros(size(grid));

    for xshift  = -Nshift:Nshift
        for yshift = -Nshift:Nshift

            shift = circshift(grid,[xshift, yshift]);

            if xshift < 0
                shift(end+xshift+1:end,:) = NaN;
            elseif xshift > 0
                shift(1:xshift,        :) = NaN;
            end
            if yshift < 0
                shift(:,end+yshift+1:end) = NaN;
            elseif yshift > 0
                shift(:,1:yshift        ) = NaN;
            end

            d_sum(~isnan(shift)) = d_sum(~isnan(shift)) + shift(~isnan(shift));
            N_sum(~isnan(shift)) = N_sum(~isnan(shift)) + 1;

            clear shift
        end
    end

    coarse                = d_sum ./ N_sum;
    coarse(isinf(coarse)) = NaN; 
    
    idx_nan = find(isnan(tile_data(ff,:)));

    for i = 1:length(idx_nan)

        tile_data_filled(ff,idx_nan(i))                                ...
            =                                                          ...
            coarse(tc.i_indg(idx_nan(i))+1, tc.j_indg(idx_nan(i))+1);
    end
    
   clear grid d_sum smooth N_sum 

end

% ===================== EOF =====================================
