%=========================================================================

function [tile_num] = get_tile_num_for_obs(tile_coord, tile_grid,...
                N_tile_in_cell_ij, tile_num_in_cell_ij,    ...
                option, this_FOV, lat, lon)
    
    N_dat = length(lat);
    
% find one tile for each obs that "administers" the obs 

%  get "max_dist" in deg lat/lon from field-of-view (FOV) 
% 
%  "max_dist" = Maximum distance allowed between obs lat/lon and tile com_lat/com_lon
%  when searching for a tile to which the obs will be assigned.
%  
%  NOTE: Subroutine get_tile_num_from_latlon() computes distances in Minkowski norm.

    if ~isempty(strfind(option,'FOV_in_deg'))

       max_dist_y          = this_FOV;
       max_dist_x(1:N_dat) = this_FOV;

    elseif ~isempty(strfind(option,'FOV_in_km'))

       % convert from [km] (FOV) to [deg] (max_dist_*)

       [max_dist_x, max_dist_y] = dist_km2deg( this_FOV, lat);

    else

       error('unknown FOV_option')

    end 

    if (max_dist_y<0. || any(max_dist_x<0.))  
        error('encountered negative max_dist');
    end

    tile_num = zeros(N_dat,1);       
    
    for n=1:N_dat
        
       % make sure lat/lon is *inside* tile_grid (to within "max_dist"),
       % otherwise do nothing

       if ( tile_grid.ll_lat           <= (lat(n)+max_dist_y   )  &&...
            tile_grid.ll_lon           <= (lon(n)+max_dist_x(n))  &&...
            (lat(n)-max_dist_y   ) <= tile_grid.ur_lat            &&...
            (lon(n)-max_dist_x(n)) <= tile_grid.ur_lon    ) 

          % min_dist = distance betw lat/lon in question and center-of-mass of
          %            matching tile 

          min_dist_x = 1.e10;    % initialize 
          min_dist_y = 1.e10;    % initialize 

          % determine grid cell that contains lat/lon 

          [i_ind,j_ind] = get_ij_ind_from_latlon( tile_grid, lat(n), lon(n));

          % make sure that i/j_ind is still within bounds 
          % (works in conjunction with if statement above re. ll/ur_lat/lon)

          i_ind = min( max(i_ind, 1), tile_grid.N_lon );
          j_ind = min( max(j_ind, 1), tile_grid.N_lat );

          % map from i_ind, j_ind to tile_num

          if   ( ~isempty(strfind(tile_grid.gridtype, 'EASE_M'))   ||    ...
                 ~isempty(strfind(tile_grid.gridtype, 'EASE-M'))   ||    ...
                 ~isempty(strfind(tile_grid.gridtype, 'EASEv2-M')) ||    ...
                 ~isempty(strfind(tile_grid.gridtype, 'EASEv2_M'))     )

             % ASSUMPTION: tiles match EASE or EASEv2 grid cells exactly
             %             (unless "outside" the domain, eg. water surface)

             if     (N_tile_in_cell_ij(i_ind,j_ind)==1) 

                tile_num(n)=tile_num_in_cell_ij(i_ind,j_ind,1);

                min_dist_x = abs(lon(n) - tile_coord.com_lon(tile_num(n)));
                min_dist_y = abs(lat(n) - tile_coord.com_lat(tile_num(n)));

             elseif (N_tile_in_cell_ij(i_ind,j_ind)==0) 

                % Do nothing.  If given EASE or EASEv2 grid cell is not land, 
                % tile_num will not change from its initialized value.

             else

                error( 'something wrong for EASE grid');
                
             end
                 
          else
             error('not ready');
          end
         
          if (tile_num(n)>0) 

             outside_bbox  = ( ...
                  lon(n) < tile_coord.min_lon(tile_num(n))   || ...
                  lon(n) > tile_coord.max_lon(tile_num(n))   || ...
                  lat(n) < tile_coord.min_lat(tile_num(n))   || ...
                  lat(n) > tile_coord.max_lat(tile_num(n)) );

             too_far_away = ( ...
                  min_dist_x > max_dist_x(n)   || ...
                  min_dist_y > max_dist_y );

             % keep tile_num unless obs is outside the bounding box *and* too far away

             if (outside_bbox && too_far_away)  
                 tile_num(n) = NaN; 
             end

          end 

        end
            
    end
    
end

%=========================================================================            
