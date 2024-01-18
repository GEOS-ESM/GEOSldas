
function out_mwRTM  = ...
             get_L2_RTM_constants_tile_data(tile_coord,L2_file_tag,...
                                    L2_version,start_time, end_time)

% function to compute the 2 RTM constant variables: Albedo and Roughness(H) based on the
% preprocessed data. Constants are taken as the long term temporal mean
% with maximum coverage across Asc/Desc passes.

% Q. Liu 18 Jul 2022

L2_Ascdes_all = {'_A_','_D_'};

out_Para = {'Albedo','Roughness'};

%L2_file_tag = 'L2_SM_P_E';
%L2_version = 'R18290';

if strcmp(L2_file_tag(end-1:end),'_E')
    resolution = 'M09';
    out_Nlon = 3856;
    out_Nlat = 1624;
else
    resolution = 'M36';
    out_Nlon = 964;
    out_Nlat = 406;
end

% provide a GEOSldas with proper tile information
%if strcmp(resolution,'M36')
%    L4_path = '/discover/nobackup/projects/gmao/smap/SMAP_Nature/SMAP_Nature_v9.x/';
%    L4_version = 'SMAP_Nature_v9.1_M36';
%    based_on_h5 = 0;
%    out_Nlon = 3856/4;
%    out_Nlat = 1624/4;
%    ftilecoord = [L4_path,L4_version,'/output/SMAP_EASEv2_M36_GLOBAL/rc_out/', ...
%    L4_version,'.ldas_tilecoord.bin'];
%else
%    L4_path = '/css/smapl4/public/L4_Products/L4_SM/';
%    L4_version = 'Vv6030';
%    based_on_h5 = 1; 
%    out_Nlon = 3856;
%    out_Nlat = 1624;
%    ftilecoord = [L4_path,L4_version,'/rc_out/SPL4SM_', L4_version,'.ldas_tilecoord.bin'];

%end

%  daily mat file path
mat_path = '/discover/nobackup/qliu/matlab/SMAP/L2L4/VOD/QC_frozen_RFI/';

int_precision   = 'int32';
float_precision = 'float32';

% read tile info for binary output
tc = tile_coord;
%tc =  read_tilecoord( ftilecoord);

% 2 parameters
for iPara = 1:length(out_Para)
    
    % 2 overpasses
    for iAD = 1:2
        
        L2_Ascdes = L2_Ascdes_all{iAD};
        L2_qc_yes = 1;
        
        dtstep = 10800;
        
        % time period for computing climatology
        %start_time.year  = 2015;
        %start_time.month = 4;
        %start_time.day   = 1;
        
%         start_time.hour  = 1;
%         start_time.min   = 30;
%         start_time.sec   = 0;
%         
%         end_time.year  = 2022;
%         end_time.month = 4;
%         end_time.day   = 1;
%         end_time.hour  = start_time.hour;
%         end_time.min   = start_time.min;
%         end_time.sec   = start_time.sec;
%         
%         start_time = get_dofyr_pentad(start_time);
%         end_time   = get_dofyr_pentad(end_time);

        % -----------------------------------------------------------------------
        % read time series of SMAP L2 fields
        if end_time.month ==1
            time_tag = [num2str(start_time.year,'%4.4d'),num2str(start_time.month,'%2.2d'), ...
                '_',num2str(end_time.year-1,'%4.4d'),'12'];
        else
            time_tag = [num2str(start_time.year,'%4.4d'),num2str(start_time.month,'%2.2d'), ...
                '_',num2str(end_time.year,'%4.4d'),num2str(end_time.month-1,'%2.2d')];
        end
        
        % file to  save long term mean tile data
        fname_clim = [mat_path,'/',out_Para{iPara},'_clim_L2_',L2_version,L2_Ascdes,resolution,'tile_',time_tag,'.bin'];
        
        L2_para_clim_sum = zeros(out_Nlon,out_Nlat);
        N_L2_clim_sum = zeros(out_Nlon,out_Nlat);
        
        % only do time loop if no previously saved climatology file is found
        if ~exist(fname_clim,'file')
            
            start_time.day = 1;
            start_time.hour = 0;
            start_time.min  = 0;
            start_time.sec  = 0;
            
            end_time.day = 1; 
            date_time = start_time;
            
            while 1
                
                if (date_time.year ==end_time.year   && ...
                        date_time.month==end_time.month  && ...
                        date_time.day  ==end_time.day    )
                    break
                end
                
                dt_tag = [num2str(date_time.year,'%4.4d'), ...
                    num2str(date_time.month,'%2.2d'), ...
                    num2str(date_time.day,'%2.2d')];
                
                mat_fname = [mat_path,'L2DCA_RTM_',L2_file_tag,'_',L2_version,L2_Ascdes, dt_tag,'.mat'];
                
                if exist(mat_fname,'file')
                    
                    disp(['loading ', mat_fname])
                    if contains(fname_clim,'Albedo_')
                        load(mat_fname,'L2_omg')
                        L2_para = L2_omg; clear L2_omg
                    elseif contains(fname_clim,'Roughness_')
                        load(mat_fname,'L2_h')
                        L2_para = L2_h; clear L2_h
                    else
                        error(['unknown clim fname ',fname_clim])
                    end
                                        
                    % compute climatology
                    tmp = L2_para;  
                    tmp(isnan(tmp)) = 0;
                    L2_para_clim_sum(:,:) = squeeze(L2_para_clim_sum(:,:)) + tmp; clear tmp
                    
                    tmp = ~isnan(L2_para);
                    N_L2_clim_sum(:,:)     = squeeze(N_L2_clim_sum(:,:)) + tmp; clear tmp
                    
                    clear L2_para
                    
                else
                    
                    error('daily mat file not found')
                    
                end
                
                date_time = augment_date_time(86400, date_time);
                
            end
            
            L2_para_clim = L2_para_clim_sum ./N_L2_clim_sum;
            L2_para_clim(N_L2_clim_sum < 1) = NaN; 
            
            % regrid to til grid
            
            L2_para_tile = NaN + ones(tc.N_tile,1);
            for k = 1:tc.N_tile
                L2_para_tile(k,1) = L2_para_clim(tc.i_indg(k)+1, tc.j_indg(k)+1);
            end
            
            ifp = fopen(fname_clim,'w','l');
            
            
                tile_data = L2_para_tile; % should be non-negative?
                
                tile_data(isnan(tile_data)) = 1.e15;
                
                fwrite( ifp, tc.N_tile*4, int_precision );% fortran_tag
                fwrite( ifp, tile_data(:), float_precision );
                fwrite( ifp, tc.N_tile*4, int_precision );% fortran_tag
                
                clear  tile_data
           
            fclose(ifp)
        end
        
    end
    
    % combine Asc and Desc data for maximum spatial coverage
    data_clim_tile = NaN + ones(tc.N_tile,2);
    for iAD = 1:2
        
        L2_Ascdes = L2_Ascdes_all{iAD};
        
        fname = [mat_path,'/',out_Para{iPara},'_clim_L2_',L2_version,L2_Ascdes,resolution,'tile_',time_tag,'.bin'];
         
        disp(['read ',fname])
        ifp = fopen(fname,'r','l');
        
        fortran_tag = fread( ifp, 1, int_precision );
        tmp = fread( ifp, tc.N_tile, float_precision );
        fortran_tag = fread( ifp, 1, int_precision );
        
        data_clim_tile(:,iAD) = tmp;
        
        fclose(ifp);
    end
    
    % fllValue =  1.e15
    data_clim_tile(data_clim_tile > 10.) = NaN;
    data_clim_tile(data_clim_tile < 0.) = 0.; %  set small negative values to 0
    
    % averaging  A, D values
    tile_data = mean(data_clim_tile,2,"omitnan");
    
    eval(['out_mwRTM.',out_Para{iPara},'=tile_data;'])
end
