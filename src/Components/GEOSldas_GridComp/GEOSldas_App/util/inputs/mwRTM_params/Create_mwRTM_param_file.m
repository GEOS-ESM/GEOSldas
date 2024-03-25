% ---------------------------------------------------------------------------
% This script is to generate the mwRTM_param.nc4 file for the GEOSldas mwRTM.
% All constant mwRTM parameters are in this file. They come from 3 different
% sources: cat_param, vegcls lookup table and preprocessed L2DCA daily mat files.
% Therefore, need to run Preprocess_L2DCA_mwRTM_into_dailymat.m before running
% this script.

% qliu + rreichle, 29 Jul 2022

% ----------------------------------------------------------------------------

clear

% add path to matlab functions in src/Components/GEOSldas_GridComp/GEOSldas_App/util/shared/matlab/
addpath('../../shared/matlab/');

% option to fill small gaps based on neighboring grids, 1 is recommended.
fill_small_gaps = 1;

% fill value in output file
fillValue = single(1.e15);

% resolution of output parameters: only works with "M09" or "M36".
EASEv2_grid = 'M09';

% older version mwRTM_param.nc4 file for parameter names and attributes
fname_in = ['/home/qliu/smap/SMAP_Nature/bcs/RTM_params/RTMParam_SMAP_L4SM_v004/SMAP_EASEv2_',EASEv2_grid,'/mwRTM_param.nc4'];

% target mwRTM_param file name
fname_out = ['/home/qliu/smap/SMAP_Nature/bcs/RTM_params/RTMParam_L2_omega_H_tmp/SMAP_EASEv2_',EASEv2_grid,'/mwRTM_param_L2_omega_H_fillValue_bhbvlewt.nc4'];

% get inputs for fill_gaps_in_tiledata() below
if fill_small_gaps

    if     strcmp(EASEv2_grid,'M09')
        N_cells = 5;
        iscube  = 0;
    elseif strcmp(EASEv2_grid,'M36')
        N_cells = 3;
        iscube  = 0;
    else
        error('invalid resolution, use ''M09'' or ''M36'' only')
    end

    tmpstr  = num2str(N_cells);

    fname_out = strrep(fname_out,'.nc4','_',tmpstr,'gx',tmpstr,'gfilled.nc4');
end

% Do not overwrite if file exists
if exist(fname_out,'file')

   disp(['file exist ',fname_out])
   return

end

% GEOSldas experiment for tilecoord and cat_params
exp_path = '/home/qliu/smap/SMAP_Nature/SMAP_Nature_v10/';

if strcmp(EASEv2_grid,'M09')
   exp_run = 'SMAP_Nature_v10.0';
   domain = 'SMAP_EASEv2_M09_GLOBAL';
else
   exp_run =  'SMAP_Nature_v10.0_M36';
   domain = 'SMAP_EASEv2_M36_GLOBAL';
end

fname_tc = [exp_path,exp_run,'/output/',domain,'/rc_out/',exp_run,'.ldas_tilecoord.bin'];
fname_tg = [exp_path,exp_run,'/output/',domain,'/rc_out/',exp_run,'.ldas_tilegrids.bin'];
fname_catparam = [exp_path,exp_run,'/output/',domain,'/rc_out/Y2015/M04/',exp_run,'.ldas_catparam.20150401_0000z.bin'];

% If use L4 products in /css/smapl4/
if ~exist(fname_tc, 'file')
    fname_tc = [exp_path,exp_run,'/rc_out/',exp_run,'.ldas_tilecoord.bin'];
    fname_tg = [exp_path,exp_run,'/rc_out/',exp_run,'.ldas_tilegrids.bin'];
    fname_catparam = [exp_path,exp_run,'/rc_out/Y2015/M04/',exp_run,'.ldas_catparam.20150401_0000z.bin'];
end

tc = read_tilecoord(fname_tc);

% double check for tile order, may not work if exp_run uses older bcs
% version
if max(abs(transpose([1:tc.N_tile])-tc.tile_id)) > 0
   error('tile order is not strictly tile_id ascending, need to modify script to reorder')
   return
end

tg = read_tilegrids(fname_tg);
cat_param = read_catparam(fname_catparam, tc.N_tile);

% L2RTM parameter source information
L2_version = 'R18290';
L2_start_time.year =  2015; L2_start_time.month = 4; L2_start_time.day = 1;
L2_end_time.year   =  2022; L2_end_time.month   = 4; L2_end_time.day   = 1;

if strcmp(EASEv2_grid,'M36')
    L2_file_tag = 'L2_SM_P';
else
    L2_file_tag = 'L2_SM_P_E';
end

% L2DCA based parameters
L2_param = get_L2_RTM_constants_tile_data(tc,L2_file_tag,...
                                    L2_version,L2_start_time, L2_end_time);
omega = L2_param.Albedo;
hparam = L2_param.Roughness;
clear  L2_param

% cat_param based parameters
mwRTMparam.soilcls = int32(cat_param.soilcls30);
mwRTMparam.sand   = cat_param.sand30/100.;
mwRTMparam.clay   = cat_param.clay30/100.;

% there are 0 values in poros30 for unknown reasons, set 0 to next minimum
% value (0.3741)
mwRTMparam.poros  = max(0.3741, cat_param.poros30);
mwRTMparam.wang_wp  = cat_param.wpwet30 .* cat_param.poros30;
mwRTMparam.wang_wt  = 0.49*mwRTMparam.wang_wp + 0.165;
mwRTMparam.rgh_wmin  = mwRTMparam.wang_wt;
mwRTMparam.rgh_wmax  = mwRTMparam.poros;

% Initialize the input and output file interface

netcdf.setDefaultFormat('FORMAT_NETCDF4');

fin_id = netcdf.open(fname_in, 'NOWRITE');
fout_id = netcdf.create(fname_out, 'NETCDF4');

if fout_id < 0, error(['Creating ' fname_out 'failed']); end

finfo = ncinfo(fname_in); netcdf.close(fin_id)

% Define Dimension (tile)
Dim_id = netcdf.defDim(fout_id,'tile',tc.N_tile);

nvar_in_file = length(finfo.Variables);

for i=1: nvar_in_file

    data_name = finfo.Variables(i).Name;
    data_type = finfo.Variables(i).Datatype;

    data_type = 'float';

    data_size = finfo.Variables(i).Size;

    varid(i) = netcdf.defVar(fout_id, data_name, data_type, Dim_id );

    netcdf.defVarFill(fout_id, varid(i), false, fillValue);

    n_attr = length(finfo.Variables(i).Attributes);

    for iv = 1:n_attr
        att_name  = finfo.Variables(i).Attributes(iv).Name;
        att_value = finfo.Variables(i).Attributes(iv).Value;

        netcdf.putAtt(fout_id, varid(i), att_name, att_value);
    end

    netcdf.endDef(fout_id);

    startVAR = repmat([0], 1, length(data_size));
    countVAR = data_size;

    % get parameter values from their respective sources
    % Total of 18 mwRTM parameters:
    % 8 from cat_param: SOILCLS, SOIL, CLAY, POROS, WANGWT, WANTWP, RGHWMIN,  RGHWMAX
    % 4 from vegcls lookup table  : VEGCLS, RGHNRH, RGHNRV,POLMIX
    % 3 from L2RTM: RGHHMIN, RGHHMAX, OMEGA
    % 3 are set to fillValue: BH,BV, OMEGA

    if strcmp(data_name,'MWRTM_OMEGA')
        if fill_small_gaps
            omega_filled = fill_gaps_in_tiledata(tc, tg, transpose(omega), N_cells, iscube );
            data = omega_filled;
        else
            data = omega;
        end
    elseif contains(data_name, 'MWRTM_RGHHM')
        if fill_small_gaps
            hparam_filled = = fill_gaps_in_tiledata(tc, tg, transpose(hparam), N_cells, iscube );
            data = hparam_filled;
        else
            data = hparam;
        end
    elseif strcmp(data_name,'MWRTM_SOILCLS')
        data  = mwRTMparam.soilcls;
    elseif strcmp(data_name,'MWRTM_SAND')
        data  = mwRTMparam.sand;
    elseif strcmp(data_name,'MWRTM_CLAY')
        data  = mwRTMparam.clay;
    elseif strcmp(data_name,'MWRTM_POROS')
        data  = mwRTMparam.poros;
    elseif strcmp(data_name,'MWRTM_WANGWT')
        data  = mwRTMparam.wang_wt;
    elseif strcmp(data_name,'MWRTM_WANGWP')
        data  = mwRTMparam.wang_wp;
    elseif strcmp(data_name,'MWRTM_RGHWMIN')
        data  = mwRTMparam.rgh_wmin;
    elseif strcmp(data_name,'MWRTM_RGHWMAX')
        data  = mwRTMparam.rgh_wmax;
    elseif strcmp(data_name,'MWRTM_LEWT') || contains(data_name, 'MWRTM_B')
        data = fillValue .* ones(data_size,1);
    else
        if strcmp(EASEv2_grid,'M09')
            dominant_M36vegcls = 1;
        else
            dominant_M36vegcls = 0;
        end
        tmp_rtm = get_mwRTM_vegcls_based( tc, dominant_M36vegcls,['EASEv2_',EASEv2_grid]);
        if strcmp(data_name,'MWRTM_RGHNRH')
            data = tmp_rtm.rgh_Nrh;
        elseif strcmp(data_name,'MWRTM_RGHNRV')
            data = tmp_rtm.rgh_Nrv;
        elseif strcmp(data_name,'MWRTM_VEGCLS')
            data = tmp_rtm.vegcls;
        elseif strcmp(data_name,'MWRTM_RGHPOLMIX')
            data = tmp_rtm.rgh_polmix;
        end
    end

    % earlier fillValue was -9999. replace with new fillValue (1.e15)
    data(abs(data-(-9999.)) < abs(-9999.*1e-4)) = NaN;

    data(isnan(data)) = fillValue;

    netcdf.putVar(fout_id, varid(i), startVAR, countVAR, data); clear data

    netcdf.reDef(fout_id);

end

netcdf.endDef(fout_id);

netcdf.close(fout_id);

disp(['done writing ',fname_out])

% --------------------------EOF--------------------------------------------
