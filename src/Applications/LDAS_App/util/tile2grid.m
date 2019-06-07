function [ grid_data ]  =                                            ...  
    tile2grid( tile_data, tile_coord, tile_grid, nodata, nodata_tol )

% Mapping from tile to grid is based on fields "i_indg" and
% "j_indg" of tilecoord structure which are in reference to 
% the *global* grid that underlies the tile definitions.
% Therefore, the input variable "tile_grid" must refer to *global*
% grid.

% reichle, 26 Jan 2006
% reichle, 25 Jul 2006 - expanded for RedArk_OSSE 
% GDL,     22 Jun 2010 - adapted for latest LDAS-tag
% reichle,  8 Jul 2010 - use "tile_grid" as input, not bkwd-compatible!
%
% -----------------------------------------------------------------

% check whether no-data variables are available on input 

if ~exist('nodata'),            nodata     = -9999;    end
if ~exist('nodata_tol'),        nodata_tol = 1e-4;     end
  
% -----------------------------------------------------

N_fields = size(tile_data,1);

% minimal check for consistency between tile_data and tile_coord 

if (size(tile_data,2)~=tile_coord.N_tile)
  
  input('tile2grid.m: Something wrong with N_tile, ctrl-c now!')
  
end

% ------------------------------------------------------

% initialize

grid_data = zeros( tile_grid.N_lon, tile_grid.N_lat, N_fields );

for k=1:N_fields
  
  wgrid = zeros( tile_grid.N_lon, tile_grid.N_lat);
  
  % loop through tile space
  
  for n=1:tile_coord.N_tile
    
    i = tile_coord.i_indg(n) - (tile_grid.i_offg - (1-tile_grid.ind_base));
    j = tile_coord.j_indg(n) - (tile_grid.j_offg - (1-tile_grid.ind_base));
    
    w = tile_coord.frac_cell(n);
    
    if (abs(tile_data(k,n)-nodata)>nodata_tol) 
      
      grid_data(i,j,k) = grid_data(i,j,k) + w*tile_data(k,n);
      
      wgrid(i,j) = wgrid(i,j) + w;
      
    end
    
  end 
  
  % normalize and set no-data-value
  
  for i=1:tile_grid.N_lon
    
    for j=1:tile_grid.N_lat
      
      if (wgrid(i,j)>0.)
	
        grid_data(i,j,k) = grid_data(i,j,k)/wgrid(i,j);
        
      else
	
        grid_data(i,j,k) = nodata;
        
      end 
      
    end
  end
  
end


grid_data(find(grid_data==nodata)) = NaN;

% ================ EOF =========================================
