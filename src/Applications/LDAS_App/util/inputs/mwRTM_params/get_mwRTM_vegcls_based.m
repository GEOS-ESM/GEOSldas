function [mwRTMparam] = get_mwRTM_vegcls_based( tile_coord, dominant_M36vegcls,resolution)

% Helper function to get RTM parameter values based on the vegcls lookup table 
% before writing in the mwRTM_params.nc4 file. The function is called in 
% Write_mwRTM_nc4_file.m

% Q. Liu 20 Jul 2022


% Vegetation-based RTM-parameters
lookup_option = 'Lit4'; 

%----------------------------------------------------------------------------

% Vegetation classes
% read in M36 & M09  vegcls for both M36 and M09 mwRTM
vegcls_fname  = '/discover/nobackup/qliu/gdelanno_RTM/SMAP_aux/EASEV2/dominantIGBP36km.406x964.uint8';

ifp     = fopen( vegcls_fname, 'r', 'b' );
veg_clsM36 = fread( ifp, [406 964] ,'uint8');
fclose(ifp);

vegcls_fname  = '/discover/nobackup/qliu/gdelanno_RTM/SMAP_aux/EASEV2/dominantIGBP09km.1624x3856.uint8';

ifp     = fopen( vegcls_fname, 'r', 'b' );
veg_clsM09 = fread( ifp, [1624 3856] ,'uint8');
fclose(ifp);

%  list of parameters can be based on vegcls lookup table
fn_mwRTMparam  = { ...
                'vegcls',   ...  
                'rgh_Nrh', 'rgh_Nrv', ...
                'rgh_polmix','omega',   'bh',      'bv',      'lewt'};

for i=1:length(fn_mwRTMparam)
          
    mwRTMparam.(fn_mwRTMparam{i})  = [NaN+zeros(length(tile_coord.N_tile),1)];
    
end

[ veg_lookup, tmp ] = get_mwRTM_lookup(lookup_option);
clear tmp

%====================================================================

%order of tiles depend on LDASsa-output; tile_coord.txt
%output from this code should *always* be sorted as [1:N_tile], that is -
%even if re-ordered output from an LDASsa-run is used as input, the
%resulting mwRTMparam.xxxx(:) is always monotonically increasingly sorted

[row_M09,col_M09] = EASEv2_latlon2ind(tile_coord.com_lat,tile_coord.com_lon,'M09',1);
row_M09 = row_M09 +1;
col_M09 = col_M09 +1;

[row_M36,col_M36] = EASEv2_latlon2ind(tile_coord.com_lat,tile_coord.com_lon,'M36',1);
row_M36 = row_M36 + 1;
col_M36 = col_M36 + 1;

for tile = 1 : tile_coord.N_tile   

        %since m10_p3, the tile-order has changed, so tile<>tile_id !
        tileid = tile_coord.tile_id(tile);
        
        if contains(resolution,'_M09') && ~dominant_M36vegcls
            
            jj = row_M09(tile);
            ii = col_M09(tile);
            i_vegcls = int32(veg_clsM09(jj,ii));
            
        else
            
            % First, rely on M36 database
            jj = row_M36(tile);
            ii = col_M36(tile);
            i_vegcls = int32(veg_clsM36(jj,ii));
            
            % Next, fill missing vegcls with the dominant vegcls of the
            % M09 grids
            
            if (i_vegcls<=0 || i_vegcls>200 || isnan(i_vegcls))
                
                if contains(resolution,'M09')
                    tmp_ind  = find(row_M36 == jj & col_M36 == ii);
                    for t = 1:length(tmp_ind)
                        i_vegcls(t) = int32(veg_clsM09(row_M09(tmp_ind(t)),col_M09(tmp_ind(t))));
                    end
                else
                    t=0;
                    for j1 = ((jj-1)*4+1):(jj*4)
                        for i1 = ((ii-1)*4+1):(ii*4)
                            t = t+1;
                            i_vegcls(t) = int32(veg_clsM09(j1,i1));
                        end
                    end
                end
                
                ind_veg  = find(~(i_vegcls<=0 | i_vegcls>200 | isnan(i_vegcls)));
                if ~isempty(ind_veg)
                    i_vegcls = i_vegcls(ind_veg);
                end
                i_vegcls = mode(i_vegcls);
                
            end
        end
        %====for valid vegetation classes in IGBP:====
        if (i_vegcls > 0 )
            
            mwRTMparam.vegcls(tileid,1)     = i_vegcls;  %IGBP
            
            mwRTMparam.rgh_Nrh(tileid,1)    = veg_lookup.rgh_Nrh(i_vegcls);
            mwRTMparam.rgh_Nrv(tileid,1)    = veg_lookup.rgh_Nrv(i_vegcls);
            
            mwRTMparam.rgh_polmix(tileid,1) = 0;
            mwRTMparam.omega(tileid,1)      = veg_lookup.omega(i_vegcls);
            mwRTMparam.bh(tileid,1)         = veg_lookup.bh(i_vegcls);
            mwRTMparam.bv(tileid,1)         = veg_lookup.bv(i_vegcls);
            mwRTMparam.lewt(tileid,1)       = veg_lookup.lewt(i_vegcls);
            
        else
            
            for i = 1:length(fn_mwRTMparam)
                mwRTMparam.(fn_mwRTMparam{i})(tileid,1) = NaN;
            end
            
        end
        
end

%===============================EOF======================================
