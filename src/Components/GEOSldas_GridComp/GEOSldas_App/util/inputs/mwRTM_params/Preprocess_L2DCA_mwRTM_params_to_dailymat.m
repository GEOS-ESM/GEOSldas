% Script to read SMAP L2 files and extract RTM variables (albedo, vegopacity, roughness)
% to store in global EASEv2 grid daily composite (Asc and Desc separately) mat files for future use.

% Q. Liu 18 Jul 2022

clear

% add path to matlab functions in src/Components/GEOSldas_GridComp/GEOSldas_App/util/shared/matlab/
addpath('../../shared/matlab/');

L2_Ascdes = {'_A_','_D_'};
L2_qc_yes = 1;

% L2 file information. Use 'L2_SM_P' for M36 resolution and 'L2_SM_P_E' for M09
SMAP_product = 'L2_SM_P'; %'L2_SM_P';
L2_path = ['/discover/nobackup/projects/gmao/smap/SMAP_L4/SMAP/OPS/',SMAP_product,'/'] ;

L2_version = 'R18290'; %'R17000';

out_path = '/discover/nobackup/qliu/matlab/SMAP/L2L4/VOD/QC_frozen_RFI/';

L2_dtstep = 10800;

start_time.year  = 2021;
start_time.month = 4;
start_time.day   = 1;

end_time.year  = 2022;
end_time.month = 4;
end_time.day   = 1;

M09_Nlon = 3856;
M09_Nlat = 1624;

M36_Nlon = M09_Nlon/4;
M36_Nlat = M09_Nlat/4;

if strcmp(SMAP_product(end-1:end),'_E')
    out_Nlon = M09_Nlon;
    out_Nlat = M09_Nlat;
else
    out_Nlon = M36_Nlon;
    out_Nlat = M36_Nlat;
end

% -----------------------------------------------------------------------
% Read L2 data
for iorb = 1:length(L2_Ascdes)

    date_time = start_time;
    date_time.hour  = 0;
    date_time.min   = 0;
    date_time.sec   = 0;

    fname_L2_pre = [];

    while 1

        if (date_time.year ==end_time.year   && ...
                date_time.month==end_time.month  && ...
                date_time.day  ==end_time.day    )
            break
        end

        outfile_tag = [num2str(date_time.year,'%4.4d'), ...
            num2str(date_time.month,'%2.2d'), ...
            num2str(date_time.day,'%2.2d')];

        mat_fname = [out_path,'L2DCA_RTM_',SMAP_product,'_',L2_version,L2_Ascdes{iorb}, outfile_tag,'.mat'];

        if exist(mat_fname,'file')

            disp(['mat file exist ',mat_fname])

        else

            L2_tau = NaN + ones(out_Nlon,out_Nlat);
            L2_omg = NaN + ones(out_Nlon,out_Nlat);
            L2_h   = NaN + ones(out_Nlon,out_Nlat);

            L2_data_path = [L2_path, '/Y', num2str(date_time.year,'%4.4d'), ...
                '/M', num2str(date_time.month, '%2.2d'), ...
                '/D', num2str(date_time.day, '%2.2d')];

            % get list of all files in subdirectory of given date
            L2_files_all = dir([L2_data_path, '/SMAP_',SMAP_product,'*',L2_Ascdes{iorb},'*_',L2_version,'_*.h5']);

            % check if multiple versions of the same half  orbit file
            % only keep the data with the highest version id

            L2_files = {};
            kk = 0; % counter  of the final file list
            % Remove duplicate files from the list when multiple versions exist
            for ff = 1:length(L2_files_all)
                % check if v002 or higher exist
                if str2num(L2_files_all(ff).name(end-4:end-3)) > 1
                    L2_files{kk} = L2_files_all(ff).name;
                    % if v002 or higher matches previous file in final flist, replace
                    % previous file in list. Only increase the final file counter when there is no duplicates
                    if ~strcmp(L2_files{kk}(1:end-5), L2_files_all(ff).name(1:end-5))
                        kk = kk + 1;
                    end
                else
                    kk = kk + 1;
                    L2_files{kk} = L2_files_all(ff).name;
                end
            end

            clear L2_fiels_all L2_fname
            if ~isempty(L2_files)

                for ifile = 1:length(L2_files)
                    L2_fname{ifile} = [L2_data_path,'/', L2_files{ifile}];
                end

                if ~isempty(fname_L2_pre)
                    L2_fname((ifile+1):(ifile+length(fname_L2_pre))) = fname_L2_pre;
                end

                fname_L2_pre = [];

                ii = 0;
                for ifile = 1: length(L2_fname)

                    fname = L2_fname{ifile};

                    disp(fname)

                    L2_row = h5read(fname,'/Soil_Moisture_Retrieval_Data/EASE_row_index'); %zero-based
                    L2_row = L2_row + 1;
                    L2_col = h5read(fname,'/Soil_Moisture_Retrieval_Data/EASE_column_index');
                    L2_col = L2_col + 1;

                    L2_utc_seconds = h5read(fname,'/Soil_Moisture_Retrieval_Data/tb_time_seconds');

                    L2_vod = h5read(fname,'/Soil_Moisture_Retrieval_Data/vegetation_opacity_option3');
                    fill_value = h5readatt(fname,'/Soil_Moisture_Retrieval_Data/vegetation_opacity_option3','_FillValue');
                    L2_vod(L2_vod == fill_value) = NaN;

                    L2_alb = h5read(fname,'/Soil_Moisture_Retrieval_Data/albedo_option3');
                    fill_value = h5readatt(fname,'/Soil_Moisture_Retrieval_Data/albedo_option3','_FillValue');
                    L2_alb(L2_alb == fill_value) = NaN;

                    L2_rough = h5read(fname,'/Soil_Moisture_Retrieval_Data/roughness_coefficient_option3');
                    fill_value = h5readatt(fname,'/Soil_Moisture_Retrieval_Data/roughness_coefficient_option3','_FillValue');
                    L2_rough(L2_rough == fill_value) = NaN;

                    % quality flag
                    L2_qf = h5read(fname,'/Soil_Moisture_Retrieval_Data/retrieval_qual_flag_option3');

                    % surface status land = 0, nonland = 1
                    L2_ss = h5read(fname,'/Soil_Moisture_Retrieval_Data/grid_surface_status');

                    % surface flag
                    L2_sf = h5read(fname,'/Soil_Moisture_Retrieval_Data/surface_flag');

                    L2_rfi_h = h5read(fname,'/Soil_Moisture_Retrieval_Data/tb_qual_flag_h');
                    L2_rfi_v = h5read(fname,'/Soil_Moisture_Retrieval_Data/tb_qual_flag_v');

                    % exclude points according to quality flag
                    if L2_qc_yes

                        % QC based on retrieval quality flag
                        %L2_rt = bitget(L2_qf, 1); % only use retrieval_recommended
                        L2_rt = bitget(L2_qf, 3); % use retrieval_succeeded

                        % QC b ased on surface flag
                        L2_frozen_model = bitget(L2_sf, 9); % model frozen ground
                        L2_snow = bitget(L2_sf,6); % snow and ice
                        L2_pice = bitget(L2_sf,7); % permanent snow and ice
                        L2_rfi_h_qf = bitget(L2_rfi_h,1); % quality flag RFI H
                        L2_rfi_v_qf = bitget(L2_rfi_v,1); % quality flag RFI V
                        idx = find(L2_rt == 0 & L2_frozen_model ==0 & ...
                            L2_snow == 0 & L2_pice ==0 & L2_ss == 0 & ...
                            L2_rfi_h_qf == 0 & L2_rfi_v_qf == 0);

                        % only keep data/coord that pass QC
                        L2_vod = L2_vod(idx);
                        L2_alb = L2_alb(idx);
                        L2_rough = L2_rough(idx);
                        L2_row = L2_row(idx);
                        L2_col = L2_col(idx);
                        L2_utc_seconds = L2_utc_seconds(idx);
                        clear idx
                    end

                    if ~isempty(L2_vod)

                        % round date_time to nearest 3 hourly UTC
                        utc_t2k = round(double(L2_utc_seconds)/L2_dtstep)*L2_dtstep;

                        [yr, doy, mm, dd, hr, mn] = J2000_to_DateTime( utc_t2k );

                        % use points for current UTC day only
                        idx = find(yr == date_time.year & mm == date_time.month & ...
                            dd == date_time.day);

                        % points across UTC days will be saved in next daily file
                        if length(idx) < length(L2_utc_seconds) && ifile <= length(L2_files)
                            disp('L2 across UTC days')
                            ii = ii + 1;
                            fname_L2_pre{ii} = fname;
                        end

                        L2_vod = L2_vod(idx);
                        L2_alb = L2_alb(idx);
                        L2_rough = L2_rough(idx);
                        L2_row = L2_row(idx);
                        L2_col = L2_col(idx);
                        L2_utc_seconds = L2_utc_seconds(idx);
                        hr = hr(idx);
                        clear idx

                        % Map L2 to 2d grid
                        for k = 1:length(L2_vod)
                            this_col = L2_col(k);
                            this_row = L2_row(k);

                            L2_tau(this_col, this_row) = L2_vod(k);
                            L2_omg(this_col, this_row) = L2_alb(k);
                            L2_h(this_col, this_row) = L2_rough(k);
                        end
                    end
                end

            else

                disp(['no L2 data found in', L2_data_path])
                pause
            end

            save(mat_fname,'L2_tau','L2_omg','L2_h')
            clear L2_tau L2_omg L2_h

        end

        date_time = augment_date_time(86400, date_time);
    end
end

% ------------------------EOF-----------------------------------------
