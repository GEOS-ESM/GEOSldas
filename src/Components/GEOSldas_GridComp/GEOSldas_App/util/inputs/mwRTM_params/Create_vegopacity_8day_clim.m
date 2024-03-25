% script to create 8-day climatology of vegetation opacity for L-band microwave
%   radiative transfer model (mwRTM)
%
% requires pre-processing of SMAP L2 data into daily *.mat files using
%   Preprocess_L2DCA_mwRTM_params_to_dailymat.m
%
% output files written in MAPL_ReadForcing format
%
% qliu + rreichle, 29 Jul 2022
%
% -------------------------------------------------------------------------------------

clear

% add path to matlab functions in src/Components/GEOSldas_GridComp/GEOSldas_App/util/shared/matlab/
addpath('../../shared/matlab/');

L2_Ascdes_all  = {'_A_','_D_'};

out_Para       = 'VOD';

L2_file_tag    = 'L2_SM_P';
L2_version     = 'R18290';

if strcmp(L2_file_tag(end-1:end),'_E')
    resolution = 'M09';
else
    resolution = 'M36';
end

out_path = '/discover/nobackup/qliu/matlab/SMAP/L2L4/VOD/QC_frozen_RFI/';

fill_small_gaps =  1;

% provide a GEOSldas simulation with matching tile information
if strcmp(resolution,'M36')
    L4_path = '/discover/nobackup/projects/gmao/smap/SMAP_Nature/SMAP_Nature_v9.x/';
    L4_version = 'SMAP_Nature_v9.1_M36';
    based_on_h5 = 0;
    out_Nlon = 3856/4;
    out_Nlat = 1624/4;
else
    L4_path = '/css/smapl4/public/L4_Products/L4_SM/';
    L4_version = 'Vv6030';
    based_on_h5 = 1;
    out_Nlon = 3856;
    out_Nlat = 1624;
end

ftilecoord = [L4_path,L4_version,'/output/SMAP_EASEv2_',resolution,'_GLOBAL/rc_out/', ...
    L4_version,'.ldas_tilecoord.bin'];
ftilegrids = [L4_path,L4_version,'/output/SMAP_EASEv2_',resolution,'_GLOBAL/rc_out/', ...
    L4_version,'.ldas_tilegrids.bin'];

if ~exist(ftilecoord,'file')
    ftilecoord = [L4_path,L4_version,'/rc_out/SPL4SM_', L4_version,'.ldas_tilecoord.bin'];
    ftilegrids = [L4_path,L4_version,'/rc_out/SPL4SM_', L4_version,'.ldas_tilegrids.bin'];
end
% read tile info for binary output
tc =  read_tilecoord( ftilecoord);
tg =  read_tilegrids( ftilegrids);

int_precision   = 'int32';
float_precision = 'float32';

for iAD = 1:2

    L2_Ascdes = L2_Ascdes_all{iAD};
    L2_qc_yes = 1;

    dtstep = 10800;

    % time period for computing climatology
    start_time.year  = 2015;
    start_time.month = 4;
    start_time.day   = 1;

    start_time.hour  = 1;
    start_time.min   = 30;
    start_time.sec   = 0;

    end_time.year    = 2022;
    end_time.month   = 4;
    end_time.day     = 1;
    end_time.hour    = start_time.hour;
    end_time.min     = start_time.min;
    end_time.sec     = start_time.sec;

    start_time       = get_dofyr_pentad(start_time);
    end_time         = get_dofyr_pentad(end_time);

    % lookup table of month and day of first day in 8-day average (non-leap year)

    clim_8d_m1 = [    1  1  1  1  2  2  2  2 ...
        3  3  3  3  4  4  4 ...
        5  5  5  5  6  6  6  6 ...
        7  7  7  7  8  8  8  8 ...
        9  9  9  9  10  10 10 ...
        11  11  11  11  12  12  12  12];

    clim_8d_d1 = [    1     9    17    25     2    10    18    26  ...
        6    14    22    30     7    15    23   ...
        1     9    17    25     2    10    18    26 ...
        4    12    20    28     5    13    21    29 ...
        6    14    22    30     8    16    24  ...
        1     9    17    25     3    11    19    27    ];

    clim_8d_m2 = [clim_8d_m1(2:46) 1];
    clim_8d_d2 = [clim_8d_d1(2:46) 1];

    % -----------------------------------------------------------------------
    % read from preprocessed daily mat file
    if end_time.month ==1
        time_tag = [num2str(start_time.year,'%4.4d'),num2str(start_time.month,'%2.2d'), ...
            '_',num2str(end_time.year-1,'%4.4d'),'12'];
    else
        time_tag = [num2str(start_time.year,'%4.4d'),num2str(start_time.month,'%2.2d'), ...
            '_',num2str(end_time.year,'%4.4d'),num2str(end_time.month-1,'%2.2d')];
    end

    fname_clim = [out_path,'/',out_Para,'_clim_L2_',L2_version,L2_Ascdes,'8d_',resolution,'tile_',time_tag,'_w24d.bin'];

    L2_tau_clim_sum = zeros(46,out_Nlon,out_Nlat);
    N_L2_clim_sum = zeros(46,out_Nlon,out_Nlat);

    if exist(fname_clim,'file')

        disp(['found preprocessed climatology file ',fname_clim])

    else

        date_time = start_time;

        while 1

            if (date_time.year ==end_time.year   && ...
                    date_time.month==end_time.month  && ...
                    date_time.day  ==end_time.day    )
                break
            end

            outfile_tag = [num2str(date_time.year,'%4.4d'), ...
                num2str(date_time.month,'%2.2d'), ...
                num2str(date_time.day,'%2.2d')];

            mat_fname = [out_path,'L2DCA_RTM_',L2_file_tag,'_',L2_version,L2_Ascdes, outfile_tag,'.mat'];

            if exist(mat_fname,'file')

                disp(['loading ', mat_fname])
                if contains(fname_clim,'VOD_')
                    load(mat_fname,'L2_tau')
                elseif contains(fname_clim,'Albedo_')
                    load(mat_fname,'L2_omg')
                    L2_tau = L2_omg; clear L2_omg
                elseif contains(fname_clim,'Roughness_')
                    load(mat_fname,'L2_h')
                    L2_tau = L2_h; clear L2_h
                else
                    error(['unknown clim fname ',fname_clim])
                end

                d_idx = find(date_time.month == clim_8d_m1 & date_time.day >= clim_8d_d1);
                if ~isempty(d_idx)
                    d_idx = d_idx(end);
                else
                    d_idx = find(date_time.month == clim_8d_m2 & date_time.day < clim_8d_d2);
                    d_idx = d_idx(1);
                end

                % compute climatology
                tmp = max(L2_tau,0);  % make sure no negative values in  tau
                tmp(isnan(tmp)) = 0;

                L2_tau_clim_sum(d_idx,:,:) = squeeze(L2_tau_clim_sum(d_idx,:,:)) + tmp; clear tmp
                tmp = ~isnan(L2_tau);
                N_L2_clim_sum(d_idx,:,:)     = squeeze(N_L2_clim_sum(d_idx,:,:)) + tmp; clear tmp

                clear L2_tau

            else

                error('daily mat file not found')

            end

            date_time = augment_date_time(86400, date_time);

        end

        L2_tau_clim = L2_tau_clim_sum ./N_L2_clim_sum;


        % regrid to til grid

        L2_tau_tile = NaN + ones(46, tc.N_tile);
        for k = 1:tc.N_tile
            L2_tau_tile(:,k) = L2_tau_clim(:,tc.i_indg(k)+1, tc.j_indg(k)+1);
        end

        ifp = fopen(fname_clim,'w','l');

        for n =  1: 48
            if n == 1
                y1 = 0;
                y2 = 1;

                nidx = 46;
                nidx_pre = 45;
                nidx_nxt = 1;
            elseif n == 2
                y1 = 1;
                y2 = 1;
                nidx = n-1;
                nidx_pre = 46;
                nidx_nxt = n;
            elseif n == 47
                y1 = 1;
                y2 = 2;

                nidx = n-1;
                nidx_pre = n-2;
                nidx_nxt = 1;
            elseif n == 48
                y1 = 2;
                y2 = 2;
                nidx = 1;
                nidx_pre = 46;
                nidx_nxt = 2;
            else
                y1 = 1;
                y2 = 1;
                nidx = n-1;
                nidx_pre = n-2;
                nidx_nxt = n;
            end
            m1 = clim_8d_m1(nidx);
            m2 = clim_8d_m2(nidx);
            d1 = clim_8d_d1(nidx);
            d2 = clim_8d_d2(nidx);

            header = [y1 m1 d1 0 0 0 y2 m2 d2 0 0 0 tc.N_tile 1];

            tile_data = mean(L2_tau_tile([nidx_pre nidx nidx_nxt],:),1,"omitnan");
            tile_data(isnan(tile_data)) = 1.e15;

            fwrite( ifp, 14*4, int_precision ); % fortran_tag
            fwrite( ifp, header, float_precision );
            fwrite( ifp, 14*4, int_precision ); % fortran_tag

            fwrite( ifp, tc.N_tile*4, int_precision );% fortran_tag
            fwrite( ifp, tile_data(:), float_precision );
            fwrite( ifp, tc.N_tile*4, int_precision );% fortran_tag

            clear header tile_data
        end
        fclose(ifp)
    end

end

% ========================
% The  final vegopacity.bin file contains data averaged (mean) of Asc
% and Des tau climatology. VOD is climatology,the other 2 parameters are
% time constant with maximum spatial coverage
data_clim_tile = NaN + ones(48, tc.N_tile,2);
for iAD = 1:2

    L2_Ascdes = L2_Ascdes_all{iAD};

    fname = [out_path,'/',out_Para,'_clim_L2_',L2_version,L2_Ascdes,'8d_',resolution,'tile_',time_tag,'_w24d.bin'];

    disp(['read ',fname])
    ifp = fopen(fname,'r','l');

    for n = 1:48

        fortran_tag = fread( ifp, 1, int_precision );
        tmp = fread( ifp, 14, float_precision );
        header(n,:) = tmp;
        fortran_tag = fread( ifp, 1, int_precision );

        fortran_tag = fread( ifp, 1, int_precision );
        tmp = fread( ifp, tc.N_tile, float_precision );
        fortran_tag = fread( ifp, 1, int_precision );

        data_clim_tile(n,:,iAD) = tmp;

    end

    fclose(ifp);
end

data_clim_tile(data_clim_tile > 10.) = NaN;
data_clim_tile(data_clim_tile <  0.) = 0.;    %  set small negative values to 0

% averaging  A, D values
tile_data = mean(data_clim_tile,3,"omitnan");
fname_out = strrep(fname, L2_Ascdes,'_AD_');
if fill_small_gaps

    if     strcmp(resolution,'M09')
        N_cells = 5;
        iscube  = 0;
    elseif strcmp(resolution,'M36')
        N_cells = 3;
        iscube  = 0;
    else
        error('invalid resolution, use ''M09'' or ''M36'' only')
    end

    tmpstr  = num2str(N_cells);

    fname_out = [fname_out(1:end-4),'_',tmpstr,'gx',tmpstr,'gfilled_test.bin'];

    tile_data = fill_gaps_in_tiledata(tc, tg, tile_data, N_cells, iscube );

end

tile_data(isnan(tile_data)) = 1.e15;   % fillValue =  1.e15

disp(['write ',fname_out])

ifp = fopen(fname_out,'w','l');
for n = 1:48

    % write header

    fwrite( ifp, 14*4,           int_precision );
    fwrite( ifp, header(n,:),    float_precision );
    fwrite( ifp, 14*4,           int_precision );

    % write science data

    fwrite( ifp, tc.N_tile*4,    int_precision );
    fwrite( ifp, tile_data(n,:), float_precision );
    fwrite( ifp, tc.N_tile*4,    int_precision );
end
fclose(ifp);

% ============================ EOF ======================================
