function [] = write_smapL4SMqa( gph_aup_lmc_fnames, tilecoord_fname, tilegrids_fname )

% THE FOLLOWING PATH SHOULD BE ADDED IN MATLAB SCRIPT THAT CALLS THIS FUNCTION
% add path to matlab functions in src/Components/GEOSldas_GridComp/GEOSldas_App/util/shared/matlab/
addpath('../shared/matlab/');

% Generate *.qa files from SMAP L4_SM "gph" or "aup" granules.
%
% Inputs:
%
%  gph_and_aup_fnames : cell array of gph and/or aup granule (file) names
%  tilecoord_fname    : "tilecoord" file name
%  tilegrids_fname    : "tilegrids" file name
%
% The "tilecoord" and "tilegrids" information must work with all
% input granules that are to be converted.
%
% For each "gph" or "aup" granule, the *.qa file is written into
% the directory that holds the granule.
%
% Currently operates on binary LDASsa aup files. (Could be changed
% but would require h5 readers for "gph" and "aup" granules.)
%
% TBD: Insert official "h5" granule name into *.qa files.
%
% reichle,   12 Feb 2014
% de lannoy, 17 Feb 2014: added lmc
% reichle,   18 Feb 2014: minor edits and clean-up
%
% ------------------------------------------------------------------------------------

% ######################### SAMPLE DRIVER SCRIPT #####################################
%
% exppath   = '/hydro/gdelanno/SMAP_Delivered/SMAP_L4_SM_D00500/';
%
% expid     = 'SMAP_D00500_L4_SM_synth_e001';
%
% expdom    = 'SMAP_EASEv2_M09_GLOBAL';
%
% yyyymm    = '200107';
%
% gphpath   = [exppath, '/', expid, '/', expdom, '/cat/ens_avg/', ...
%              '/Y', yyyymm(1:4), '/M', yyyymm(5:6), '/'];
%
% auppath   = [exppath, '/', expid, '/', expdom, '/ana/ens_avg/', ...
%              '/Y', yyyymm(1:4), '/M', yyyymm(5:6), '/'];
%
% lmcpath   = [exppath, '/', expid, '/', expdom, '/rc_out/', ...
%              '/Y', yyyymm(1:4), '/M', yyyymm(5:6), '/'];
%
% tcpath    = [exppath, '/', expid, '/', expdom, '/rc_out/'];
%
% tc_fname  = [tcpath, expid, '.ldas_tilecoord.bin'];
% tg_fname  = [tcpath, expid, '.ldas_tilegrids.bin'];
%
% gphnames  = ...
%     {[expid,'.ens_avg.ldas_tile_xhourly_out.20010725_0130z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_xhourly_out.20010725_0430z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_xhourly_out.20010725_0730z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_xhourly_out.20010725_1030z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_xhourly_out.20010725_1330z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_xhourly_out.20010725_1630z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_xhourly_out.20010725_1930z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_xhourly_out.20010725_2230z.bin']};
%
% aupnames  = ...
%     {[expid,'.ens_avg.ldas_tile_inst_smapL4SMaup.20010725_0300z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_inst_smapL4SMaup.20010725_0600z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_inst_smapL4SMaup.20010725_0900z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_inst_smapL4SMaup.20010725_1200z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_inst_smapL4SMaup.20010725_1500z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_inst_smapL4SMaup.20010725_1800z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_inst_smapL4SMaup.20010725_2100z.bin'], ...
%      [expid,'.ens_avg.ldas_tile_inst_smapL4SMaup.20010726_0000z.bin']};
%
% % there should only be one "lmc" file:
%
% lmcname = [expid,'.ldas_smapL4SMlmc.20010725_0000z.bin'];
%
% % concatenate all file names into one cell array
%
% for ii=1:length(gphnames)
%
%   fnames{ii} = [gphpath, gphnames{ii}];
%
% end
%
% ii_off = length(gphnames);
%
% for ii=1:length(aupnames)
%
%   fnames{ii+ii_off} = [auppath, aupnames{ii}];
%
% end
%
% ii_off = ii_off + length(aupnames);
%
% fnames{ii_off+1} = [lmcpath, lmcname];
%
% write_smapL4SMqa( fnames, tc_fname, tg_fname );
%
% ##################### END SAMPLE DRIVER SCRIPT #####################################

% ------------------------------------------------------------------------------------
%
% make sure that list of input file names is a cell array

if ~iscell(gph_aup_lmc_fnames)

  error('write_smapL4SMqa.m: input list of file names must be a cell array')

end

% ----------------------------------------------

% read tile coordinate information

tile_coord = read_tilecoord(tilecoord_fname);

% Make sure LDASsa output was for EASEv2_M09 tile space

[ tile_grid_g, tile_grid_d ] = read_tilegrids(tilegrids_fname);

if isempty(strfind(tile_grid_d.gridtype, 'EASEv2_M09'))
  error('Expecting aup file in EASEv2_M09 tile space');
end

% ----------------------------------------------

% process each file in list of input files

for ii=1:length(gph_aup_lmc_fnames)

  this_fname = gph_aup_lmc_fnames{ii};

  disp(['processing: ', this_fname])

  % parse input string name to decide whether a "gph" or "aup" granule is to be converted

  if     any( findstr( this_fname , 'ldas_tile_xhourly_out'      ))

    get_gph_qa( this_fname, tile_coord );

  elseif any( findstr( this_fname , 'ldas_tile_inst_smapL4SMaup' ))

    get_aup_qa( this_fname, tile_coord );

  elseif any( findstr( this_fname , 'ldas_smapL4SMlmc' ))

    get_lmc_qa( this_fname, tile_coord );

  else

    error('write_smapL4SMqa.m: something wrong with input file name')

  end

  disp(['---------------------------------------------------------'])

end


% *********************************************************************************************
% *********************************************************************************************
% *********************************************************************************************


function [] = get_gph_qa( gph_fname, tile_coord )

%====================================================================
%
% Matlab function to produce .qa-files for SMAP L4_SM *gph* output.
%
% 30jan14: Gabrielle De Lannoy - initial draft
%  1feb14: Gabrielle De Lannoy - use land fraction to calculate stats
%                              - text edits, formatting
%  7feb14: Gabrielle De Lannoy - edits
%
% - currently operates on binary LDASsa gph files
% - TBD: file name change from "bin" to official "h5" granule name
%
%====================================================================
%
% [QA] -- SMAP_L4_SM_PSD p.29:
%
% "...
%  The QA file contains statistical information that will enable users
%  to better assess the quality of the associated granule.
%  QA products bear exactly the same name as the products [(.h5)]
%  that they represent. The only difference in names is the extension.
%  The extension for all QA products is *.qa.
%  ..."
%
% [gph] -- SMAP_L4_SM_PSD p.4:
%  "...
% The first Collection is a series of 3-hourly time average geophysical ("gph")
% land surface fields that are output by the L4_SM algorithm. This
% Collection will be of primary interest to most users.
%  ..."
%
%====================================================================

check_on       = 1; %1 = LDASsa sanity checks, write warnings
                    %2 = LDASsa sanity checks, stop if check fails

nodata_val     = -9999;
nodata_tol     = 1E-4;

tablefields    = {'Fieldname','Units',...
                  'Mean','Std-dev','Min','Max','N'};

Nstat          = length(tablefields)-2;

str_l          = 51;
unt_l          = 16;
num_l          = 13;
num_s_l        = 9;

str_f          = num2str(str_l);
unt_f          = num2str(unt_l);
num_f          = num2str(num_l);
num_s_f        = num2str(num_s_l);

delim          = ',';

tableformat    = [ '%-',str_f,'s',delim,'%-',unt_f,'s'];
tableformat_sc = [ '%-',str_f,'s',delim,'%-',unt_f,'s'];

for f=1:Nstat
    if f~=Nstat
        tableformat    = [tableformat,    delim,'%',num_f,'.4f'];
        tableformat_sc = [tableformat_sc, delim,'%',num_f,'.4e'];
    else
        tableformat    = [tableformat,    delim,'%',num_s_f,'d\n'];
        tableformat_sc = [tableformat_sc, delim,'%',num_s_f,'d\n'];
    end
end

out_collection_ID = 6;

%N_out_fields = 40;  % for raw LDASsa output (EXCL. sm in pctl units)
N_out_fields = 42;  % for post-processed LDASsa output (INCL. sm in pctl units)

%====================================================================
%
% Grid information

N_gridcells_M09 = length(tile_coord.com_lat);

weights         = tile_coord.frac_cell;

%====================================================================

% expect LDASsa *ldas_tile_xhourly_out*.bin file (not "h5" output)

% read LDASsa binary gph file and check fieldnames

[fn, units] = get_data_tag(out_collection_ID, N_out_fields);

tile_data   = read_tile_data(gph_fname, tile_coord.N_tile, N_out_fields);

if tile_coord.N_tile ~= size(tile_data,2)
  error('Number of land tiles does not match the length of simulated vectors');
end

out_fname = [gph_fname(1:end-4), '.qa' ];

% assemble placeholder h5 file name

%%ind = findstr(gph_fname,'/');

%%h5_fname = [gph_fname(ind(end)+1:end-4), '.h5'];

h5_fname = '<L4_SM_GPH_GRANULE_NAME.h5>';

%============================================================
% OUTPUT
%============================================================

ofp = fopen( out_fname, 'w' );

disp(['writing ',out_fname]);

% header information: 4 lines

fprintf(ofp, ['%s\n'],...
        ['Quality Assessment for SMAP L4_SM Granule ', h5_fname]);
fprintf(ofp, ['%s%8d\n\n'],...
        ['Number of L4_SM EASEv2  9 km land grid cells = '], N_gridcells_M09);

% comma-delimited table:
% - 4 header lines (observation space)
% - X variable lines
% - footnotes

%=================================================================
% Model space
%=================================================================

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','='  );
end
fprintf(ofp, '\n%s\n',...
        'Geophysical variables');
for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','='  );
end
fprintf(ofp, '\n');

for f = 1:Nstat+2
  if f==1             fprintf(ofp, ['%-',str_f,'s',delim,''], [tablefields{f},''     ]); end
  if f==2             fprintf(ofp, ['%-',unt_f,'s',delim,''], [tablefields{f},' (*1)']); end
  if f==3 || f==4     fprintf(ofp, ['%', num_f,'s',delim,''], [tablefields{f},' (*2)']); end
  if f>4 && f<Nstat+2 fprintf(ofp, ['%', num_f,'s',delim,''], tablefields{f}); end
  if f==Nstat+2       fprintf(ofp, ['%', num_s_f,'s\n'],      [tablefields{f},' (*3)']); end
end

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','-'  );
end
fprintf(ofp, '\n');

%-----Raw gph-fields----------------------------------------------

for f = 1:length(fn)

  if strcmp(units{f},'[kg m-2 s-1]') || strcmp(units{f},'[kg kg-1]')

    fprintf(ofp, tableformat_sc, ...
            fn{f}, units{f}, ...
            stats_array(tile_data(f,:)',weights, Nstat));

  else

    fprintf(ofp, tableformat, ...
            fn{f}, units{f}, ...
            stats_array(tile_data(f,:)',weights, Nstat));

  end

end

%-----Footnotes---------------------------------------------------

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','-'  );
end

fprintf(ofp, '\n%s\n\n', ...
        'See SMAP L4_SM Data Products Specification Document for additional information.');

fprintf(ofp, '%s\n', ...
        '(*1) Units are valid for all statistics except N [dimensionless].');

fprintf(ofp, '%s\n',...
        '(*2) Mean and std-dev statistics are weighted by the land fraction of each grid cell.');

fprintf(ofp, '%s\n', ...
        '(*3) N is the number of 9 km EASEv2 grid cells that contribute to the statistics.');

fclose(ofp);



% *********************************************************************************************
% *********************************************************************************************
% *********************************************************************************************

function [] = get_aup_qa( aup_fname, tile_coord )

%====================================================================
%
% Matlab function to produce .qa-files for SMAP L4_SM *aup* output.
%
% 16jan14: Gabrielle De Lannoy - Initial draft
% 28jan14: Rolf Reichle - removed time loop etc, converted to "function"
% 30jan14: Gabrielle De Lannoy - moved headers
%  1feb14: Gabrielle De Lannoy - use land fraction to calculate stats
%                              - text edits, formatting
%  4feb14: Gabrielle De Lannoy - new aup structure
%  7feb14: Gabrielle De Lannoy - edits
%
% - currently operates on binary LDASsa aup files
% - TBD: file name change from "bin" to official "h5" granule name
%
%====================================================================
%
% [QA] -- SMAP_L4_SM_PSD p.29:
%
% "...
%  The QA file contains statistical information that will enable users
%  to better assess the quality of the associated granule.
%  QA products bear exactly the same name as the products [(.h5)]
%  that they represent. The only difference in names is the extension.
%  The extension for all QA products is *.qa.
%  ..."
%
% [aup] -- SMAP_L4_SM_PSD p.4:
%  "...
%  The second Collection provides diagnostics from the land surface analysis
%  updates ("aup"). This Collection consists of a series of 3-hourly
%  instantaneous (or snapshot) files that contain the assimilated SMAP
%  observations, the corresponding land model predictions and analysis
%  estimates, and additional data assimilation diagnostics.
%  ..."
%
%====================================================================

check_on       = 1; %1 = LDASsa sanity checks, write warnings
                    %2 = LDASsa sanity checks, stop if check fails

nodata_val     = -9999;
nodata_tol     = 1E-4;

res_flag_tag   = {'36km', '09km'};
orb_flag_tag   = {'AD', 'A', 'D'};
obs_pol        = {'h' , 'v'};

tablefields    = {'Fieldname','Units',...
                  'Mean','Std-dev','Min','Max','N'};

Nstat          = length(tablefields)-2;

str_l          = 51;
unt_l          = 16;
num_l          = 13;
num_s_l        = 9;

str_f          = num2str(str_l);
unt_f          = num2str(unt_l);
num_f          = num2str(num_l);
num_s_f        = num2str(num_s_l);

delim          = ',';

tableformat    = [ '%-',str_f,'s',delim,'%-',unt_f,'s'];
tableformat_sc = [ '%-',str_f,'s',delim,'%-',unt_f,'s'];

for f=1:Nstat
    if f~=Nstat
        tableformat    = [tableformat,    delim,'%',num_f,'.4f'];
        tableformat_sc = [tableformat_sc, delim,'%',num_f,'.4e'];
    else
        tableformat    = [tableformat,    delim,'%',num_s_f,'d\n'];
        tableformat_sc = [tableformat_sc, delim,'%',num_s_f,'d\n'];
    end
end

AmF_threshold(17:19) = 1E-5; % soil moisture
AmF_threshold(20:21) = 1E-3; % temperature

%====================================================================
%
% Grid information

% Map M09 tiles to M36 grid cells to verify the M36 obs

[M36_row,M36_col]   = ...
    EASEv2_latlon2ind(tile_coord.com_lat,tile_coord.com_lon,'M36',1);

M36_row_col         = [M36_row M36_col];
[unique_rc, ia, ic] = unique(M36_row_col,'rows');

N_gridcells_M36     = length(unique_rc);
N_gridcells_M09     = length(tile_coord.com_lat);

weights             = tile_coord.frac_cell;

%====================================================================

% expect LDASsa *aup*.bin file (not "h5" output)

% read LDASsa binary aup file and check fieldnames

[aup, units] = read_smapL4SMaup( aup_fname, tile_coord.N_tile );

fn   = fieldnames(aup);
for f = 1:length(units)
  if strcmp(units{f},'[-]')
    units{f} = '[dimensionless]';
  end
end

if tile_coord.N_tile ~= length(aup.sm_surface_wetness_forecast)
  error('Number of land tiles does not match the length of simulated vectors');
end

out_fname = [aup_fname(1:end-4), '.qa' ];

% assemble placeholder h5 file name

%%ind = findstr(aup_fname,'/');

%%h5_fname = [aup_fname(ind(end)+1:end-4), '.h5'];

h5_fname = '<L4_SM_AUP_GRANULE_NAME.h5>';

%============================================================
% OUTPUT
%============================================================

ofp = fopen( out_fname, 'w' );

disp(['writing ',out_fname]);

% header information: 4 lines

fprintf(ofp, ['%s\n'],...
        ['Quality Assessment for SMAP L4_SM Granule ', h5_fname]);
fprintf(ofp, ['%s%8d\n%s%8d\n\n'],...
        ['Number of L4_SM EASEv2  9 km land grid cells = '], N_gridcells_M09,...
        ['Number of L4_SM EASEv2 36 km land grid cells = '], N_gridcells_M36);

% comma-delimited table:
% - 4 header lines (observation space)
% - X variable lines
% - footnotes
% - 4 header lines (model space)
% - X variable lines
% - footnotes

%=================================================================
% Observation space
%=================================================================

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','='  );
end
fprintf(ofp, '\n%s\n',...
        'Brightness temperatures ("EnKF observation space")');
for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','='  );
end
fprintf(ofp, '\n');

for f = 1:Nstat+2
  if f==1             fprintf(ofp, ['%-',str_f,'s',delim,''], [tablefields{f},' (*1)']); end
  if f==2             fprintf(ofp, ['%-',unt_f,'s',delim,''], [tablefields{f},' (*2)']); end
  if f==3 || f==4     fprintf(ofp, ['%', num_f,'s',delim,''], [tablefields{f},' (*3)']); end
  if f>4 && f<Nstat+2 fprintf(ofp, ['%', num_f,'s',delim,''], tablefields{f}); end
  if f==Nstat+2       fprintf(ofp, ['%', num_s_f,'s\n'],      [tablefields{f},' (*4)']); end
end

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','-'  );
end
fprintf(ofp, '\n');

%-----Raw aup-fields----------------------------------------------

for f = 7:16

  var = getfield(aup, fn{f});
  var(abs(var-nodata_val)<nodata_tol) = NaN;
  aup.(fn{f}) = var; %replace original -9999 with NaN for future calcs

  % Tb obs: - split up by original resolution (L1C=1, L2AP=2)
  %         - split up by ascending/descending/both orbit directions

  if (isempty(strfind(fn{f},'tb_')))

    fprintf(ofp, tableformat, ...
            fn{f}, units{f}, stats_array(var, weights, Nstat));
    disp(['WARNING: Expecting tb-fields in the first block of variables']);

  else

    if (~isempty(strfind(fn{f},'tb_h')))
      res_flag = aup.tb_h_resolution_flag;
      orb_flag = aup.tb_h_orbit_flag;
    elseif (~isempty(strfind(fn{f},'tb_v')))
      res_flag = aup.tb_v_resolution_flag;
      orb_flag = aup.tb_v_orbit_flag;
    else
      error('unknown resolution_flag')
    end

    for r=1:length(res_flag_tag)
    for d=1:length(orb_flag_tag)

      % Reduce M36 observations, partitioned to M09 in the aup-file,
      % to the actual M36 resolution.
      % Check if the partitioned M09 obs are identical (stdv ~ 0 K)
      % inside their corresponding M36 grid cell.

      %M36
      if strcmp(res_flag_tag{r},'36km')

        if (strcmp(orb_flag_tag{d},'AD'))

          subset     = (res_flag==r & ~isnan(var));
          out_orb_flag_tag = '';

        else

          subset     = (res_flag==r & orb_flag==d-1 & ~isnan(var));
          out_orb_flag_tag = ['_',orb_flag_tag{d}];

        end

        M36_row_col_obs = M36_row_col(subset,:);

        tmp_var = map_M09toM36(    var(subset),  M36_row_col_obs, 0, check_on );
        tmp_w   = map_M09toM36( weights(subset), M36_row_col_obs, 1, 0 );

        if isempty(tmp_w)

            tmp_var = NaN;
            tmp_w   = NaN;

        end

      %M09
      else

        if (strcmp(orb_flag_tag{d},'D'))

          subset  = (res_flag==r & orb_flag==d-1);
          out_orb_flag_tag = ['_',orb_flag_tag{d}];

          tmp_var = var(subset);
          tmp_w   = weights(subset);

          if isempty(tmp_w)

            tmp_var = NaN;
            tmp_w   = NaN;

          end

        else

          tmp_var = nodata_val;
          tmp_w   = -1;

        end

      end

      if tmp_w ~= -1

          fprintf(ofp, tableformat, ...
              [fn{f},'_',res_flag_tag{r},out_orb_flag_tag], ...
              units{f}, stats_array(tmp_var, tmp_w, Nstat));

      end

    end
    end

  end

end

%-----Obs-minus-Forecast (innovations)----------------------------

for p=1:length(obs_pol);

  cmd = ['var      = aup.tb_',obs_pol{p},'_obs_assim - aup.tb_',obs_pol{p},'_forecast;'];
  eval(cmd);
  cmd = ['res_flag = aup.tb_',obs_pol{p},'_resolution_flag;'];
  eval(cmd);
  cmd = ['orb_flag = aup.tb_',obs_pol{p},'_orbit_flag;'];
  eval(cmd);

  for r=1:length(res_flag_tag)
  for d=1:length(orb_flag_tag)

      if strcmp(res_flag_tag{r},'36km')

        if (strcmp(orb_flag_tag{d},'AD'))
          subset     = (res_flag==r & ~isnan(var));
          out_orb_flag_tag = '';
        else
          subset     = (res_flag==r & orb_flag==d-1 & ~isnan(var));
          out_orb_flag_tag = ['_',orb_flag_tag{d}];
        end

        M36_row_col_obs = M36_row_col(subset,:);

        tmp_var = map_M09toM36(    var(subset),  M36_row_col_obs, 0, check_on );
        tmp_w   = map_M09toM36( weights(subset), M36_row_col_obs, 1, 0 );

        if isempty(tmp_w)

            tmp_var = NaN;
            tmp_w   = NaN;

        end

      else

        if (strcmp(orb_flag_tag{d},'D'))

          subset  = (res_flag==r & orb_flag==d-1);

          out_orb_flag_tag = ['_',orb_flag_tag{d}];

          tmp_var = var(subset);
          tmp_w   = weights(subset);

          if isempty(tmp_w)

            tmp_var = NaN;
            tmp_w   = NaN;

          end

        else

          tmp_var = nodata_val;
          tmp_w   = -1;

        end

      end

      if tmp_w~=-1

        fprintf(ofp, tableformat, ...
            ['tb_',obs_pol{p},'_obs_assim_minus_forecast_',res_flag_tag{r},out_orb_flag_tag], ...
            '[K]',stats_array(tmp_var, tmp_w, Nstat));

      end

  end
  end

end

%-----Normalized Obs-minus-Forecast (normalized innovations)------

for p=1:length(obs_pol);

  cmd = ['var      = (aup.tb_',obs_pol{p},'_obs_assim - aup.tb_',obs_pol{p},'_forecast)', ...
         './sqrt(aup.tb_',obs_pol{p},'_obs_errstd.^2 + aup.tb_',obs_pol{p},'_forecast_ensstd.^2);'];
  eval(cmd);
  cmd = ['res_flag = aup.tb_',obs_pol{p},'_resolution_flag;'];
  eval(cmd);
  cmd = ['orb_flag = aup.tb_',obs_pol{p},'_orbit_flag;'];
  eval(cmd);

  for r=1:length(res_flag_tag)
  for d=1:length(orb_flag_tag)

      if strcmp(res_flag_tag{r},'36km')

        if (strcmp(orb_flag_tag{d},'AD'))
          subset     = (res_flag==r & ~isnan(var));
          out_orb_flag_tag = '';
        else
          subset     = (res_flag==r & orb_flag==d-1 & ~isnan(var));
          out_orb_flag_tag = ['_',orb_flag_tag{d}];
        end

        M36_row_col_obs = M36_row_col(subset,:);

        tmp_var = map_M09toM36(    var(subset),  M36_row_col_obs, 0, check_on );
        tmp_w   = map_M09toM36( weights(subset), M36_row_col_obs, 1, 0 );

        if isempty(tmp_w)

            tmp_var = NaN;
            tmp_w   = NaN;

        end

      else

        if (strcmp(orb_flag_tag{d},'D'))
          subset  = (res_flag==r & orb_flag==d-1);
          out_orb_flag_tag = ['_',orb_flag_tag{d}];

          tmp_var = var(subset);
          tmp_w   = weights(subset);

          if isempty(tmp_w)

            tmp_var = NaN;
            tmp_w   = NaN;

          end

        else
          tmp_var = nodata_val;
          tmp_w   = -1;
        end

      end

      if tmp_w~=-1

        fprintf(ofp, tableformat, ...
            ['tb_',obs_pol{p},'_norm_obs_assim_minus_forecast_',res_flag_tag{r},out_orb_flag_tag],...
            '[K K-1]', stats_array(tmp_var, tmp_w, Nstat));

      end

  end
  end

end

%-----Footnotes---------------------------------------------------

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','-'  );
end

fprintf(ofp, '\n%s\n\n', ...
        'See SMAP L4_SM Data Products Specification Document for additional information.');

fprintf(ofp, '%s\n     %s\n     %s\n     %s\n     %s\n     %s\n     %s\n     %s\n     %s\n     %s\n     %s\n     %s\n', ...
        ['(*1) Fieldnames that contain "_h" and "_v" are for H-polarization and',...
        ' V-polarization brightness temperature (Tb) data, respectively.'],...
        ['"tb_h_obs" and "tb_v_obs" are observed Tbs obtained from SMAP L1C_TB and L2_SM_AP files',...
        ' after quality control based on'],...
        ['information from the same files.'],...
        ['"tb_h_obs_assim" and "tb_v_obs_assim" are observed Tbs that were assimilated',...
        ' in the L4_SM system after land model-based'],...
        ['quality control and climatological adjustment (scaling).'],...
        ['Fieldnames that contain "_36km" or "_09km" provide statistics for Tbs from',...
        ' SMAP L1C_TB files or L2_SM_AP files, respectively, '],...
        ['and corresponding model Tbs.'],...
        ['Fieldnames that contain "_A" or "_D" provide statistics that are masked to Tbs from ascending',...
        ' or descending orbits, respectively.'],...
        ['Some observations at very high latitudes may have resulted',...
        ' from averaging over both ascending and descending orbits.'],...
        ['Only descending orbits are available',...
        ' for 9 km Tb observations from SMAP L2_SM_AP files.'],...
        ['Fieldnames that contain "norm_obs_assim_minus_forecast" provide statistics for',...
        ' normalized Tb innovations, defined as'],...
        ['(tb_obs_assim - tb_forecast)/sqrt(tb_obs_errstd^2 + tb_forecast_ensstd^2).']);

fprintf(ofp, '%s\n',...
        '(*2) Units are valid for all statistics except N [dimensionless].');

fprintf(ofp, '%s\n',...
        '(*3) Mean and std-dev statistics are weighted by the land fraction in each grid cell.');

fprintf(ofp, '%s\n     %s\n\n',...
        ['(*4) N is the number of EASEv2 grid cells that contribute to the statistics.'],...
        ['For fieldnames containing "_36km" and "_09km", N is the number of contributing',...
        ' 36 km and 9 km EASEv2 grid cells, respectively.']);

%=================================================================
% Model space
%=================================================================

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','=');
end
fprintf(ofp, '\n%s\n',...
        'Geophysical variables ("EnKF state space")');
for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','=');
end
fprintf(ofp, '\n');

for f = 1:Nstat+2
  if f==1             fprintf(ofp, ['%-',str_f,'s',delim,''], [tablefields{f},' (*5)']); end
  if f==2             fprintf(ofp, ['%-',unt_f,'s',delim,''], [tablefields{f},' (*6)']); end
  if f==3 || f==4     fprintf(ofp, ['%', num_f,'s',delim,''], [tablefields{f},' (*7)']); end
  if f>4 && f<Nstat+2 fprintf(ofp, ['%', num_f,'s',delim,''], tablefields{f}); end
  if f==Nstat+2       fprintf(ofp, ['%', num_s_f,'s\n'],     [tablefields{f},' (*8)']); end
end

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','-');
end
fprintf(ofp, '\n');

%-----Get joint increment mask------------------------------------

for f = 17:21

  cmd = ['var = (aup.',fn{f+5},'-aup.',fn{f},');'];
  eval(cmd);

  tmp_subset = abs(var) > AmF_threshold(f);

  if f==17
    subset_exclzero = tmp_subset;
  else
    subset_exclzero = (subset_exclzero | tmp_subset);
  end

end

%-----Raw aup-fields----------------------------------------------

for f = 17:length(fn)

  var = getfield(aup, fn{f});
  var(abs(var-nodata_val)<nodata_tol) = NaN;
  aup.(fn{f}) = var; %replace original -9999 with NaN for future calcs

  if (isempty(strfind(fn{f},'obs')) && isempty(strfind(fn{f},'tb')))

    fprintf(ofp, tableformat, ...
            fn{f}, units{f}, stats_array(var, weights, Nstat));

    % only diagnose where non-zero increments are found

    var    = var(subset_exclzero);
    tmp_w  = weights(subset_exclzero);

    fprintf(ofp, tableformat, ...
            [fn{f},'_masked'], units{f}, stats_array(var, tmp_w, Nstat));

  else

    error('No variables in observation space expected');

  end

end


%-----Increments (analysis-forecast)------------------------------

for f = 17:21

  cmd = ['var = (aup.',fn{f+5},'-aup.',fn{f},');'];
  eval(cmd);

  if ~strcmp(fn{f+5}(1:end-length('analysis')),fn{f}(1:end-length('forecast')))
    error('AmF: incorrect fieldnames');
  end

  fprintf(ofp, tableformat_sc, ...
          ['analysis_minus_forecast_',fn{f+5}(1:end-length('analysis')-1)], ...
          units{f}, stats_array(var, weights, Nstat));

  % only diagnose non-zero increments

  var   = var(subset_exclzero);
  tmp_w = weights(subset_exclzero);

  fprintf(ofp, tableformat_sc, ...
          ['analysis_minus_forecast_',fn{f+5}(1:end-length('analysis')-1),'_masked'],...
          units{f}, stats_array(var, tmp_w, Nstat));

end

%-----Footnotes---------------------------------------------------

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','-'  );
end

fprintf(ofp, '\n%s\n\n', ...
        'See SMAP L4_SM Data Products Specification Document for additional information.');

fprintf(ofp, '%s\n     %s\n     %s\n', ...
        ['(*5) For fieldnames ending in "_masked", statistics are masked to areas',...
        ' with non-zero increments, where zero is defined to be within'],...
        ['+/-10e-5 [dimensionless] for soil moisture variables and +/-10e-3 [K] for temperature variables.'],...
        ['Fieldnames starting with "analysis_minus_forecast" provide statistics for analysis increments.']);

fprintf(ofp, '%s\n', ...
        '(*6) Units are valid for all statistics except N [dimensionless].');

fprintf(ofp, '%s\n', ...
        '(*7) Mean and std-dev statistics are weighted by the land fraction of each grid cell.');

fprintf(ofp, '%s\n', ...
        '(*8) N is the number of 9 km EASEv2 grid cells that contribute to the statistics.');

fclose(ofp);


% *********************************************************************************************
% *********************************************************************************************
% *********************************************************************************************

function [] = get_lmc_qa( lmc_fname, tile_coord )

%====================================================================
%
% Matlab function to produce .qa-files for SMAP L4_SM *lmv* output.
%
% 17feb14: Gabrielle De Lannoy - Initial draft
%
% - currently operates on binary LDASsa lmc files
% - TBD: file name change from "bin" to official "h5" granule name
%
%====================================================================
%
% [QA] -- SMAP_L4_SM_PSD p.29:
%
% "...
%  The QA file contains statistical information that will enable users
%  to better assess the quality of the associated granule.
%  QA products bear exactly the same name as the products [(.h5)]
%  that they represent. The only difference in names is the extension.
%  The extension for all QA products is *.qa.
%  ..."
%
% [aup] -- SMAP_L4_SM_PSD p.4:
%  "...
%  The third Collection provides static (time-invariant) land surface model
%  constants ("lmc") that will be needed by some users for further
%  interpretation of the geophysical land surface fields. This Collection
%  consists of only one granule (file) per L4_SM data product release (or data
%  product version).
%  ..."
%
%====================================================================

nodata_val     = -9999;
nodata_tol     = 1E-4;

tablefields    = {'Fieldname','Units',...
                  'Mean','Std-dev','Min','Max','N'};

Nstat          = length(tablefields)-2;

str_l          = 51;
unt_l          = 16;
num_l          = 13;
num_s_l        = 9;

str_f          = num2str(str_l);
unt_f          = num2str(unt_l);
num_f          = num2str(num_l);
num_s_f        = num2str(num_s_l);

delim          = ',';

tableformat    = [ '%-',str_f,'s',delim,'%-',unt_f,'s'];
tableformat_sc = [ '%-',str_f,'s',delim,'%-',unt_f,'s'];

for f=1:Nstat
    if f~=Nstat
        tableformat    = [tableformat,    delim,'%',num_f,'.4f'];
        tableformat_sc = [tableformat_sc, delim,'%',num_f,'.4e'];
    else
        tableformat    = [tableformat,    delim,'%',num_s_f,'d\n'];
        tableformat_sc = [tableformat_sc, delim,'%',num_s_f,'d\n'];
    end
end

%====================================================================
%
% Grid information

N_gridcells_M09 = length(tile_coord.com_lat);

weights         = tile_coord.frac_cell;

%====================================================================

% expect LDASsa *ldas_smapL4SMlmc*.bin file (not "h5" output)

% read LDASsa binary gph file and check fieldnames

[lmc, units] = read_smapL4SMlmc ( lmc_fname , tile_coord.N_tile);

fn   = fieldnames(lmc);

if tile_coord.N_tile ~= length(lmc.(fn{1}))
  error('Number of land tiles does not match the length of parameter vectors');
end

out_fname = [lmc_fname(1:end-4), '.qa' ];

% assemble placeholder h5 file name

%%ind = findstr(gph_fname,'/');

%%h5_fname = [gph_fname(ind(end)+1:end-4), '.h5'];

h5_fname = '<L4_SM_LMC_GRANULE_NAME.h5>';

%============================================================
% OUTPUT
%============================================================

ofp = fopen( out_fname, 'w' );

disp(['writing ',out_fname]);

% header information: 4 lines

fprintf(ofp, ['%s\n'],...
        ['Quality Assessment for SMAP L4_SM Granule ', h5_fname]);
fprintf(ofp, ['%s%8d\n\n'],...
        ['Number of L4_SM EASEv2  9 km land grid cells = '], N_gridcells_M09);

% comma-delimited table:
% - 4 header lines (observation space)
% - X variable lines
% - footnotes

%=================================================================
% Model space (parameters)
%=================================================================

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','='  );
end
fprintf(ofp, '\n%s\n',...
        'Land model constants');
for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','='  );
end
fprintf(ofp, '\n');

for f = 1:Nstat+2
  if f==1             fprintf(ofp, ['%-',str_f,'s',delim,''], [tablefields{f},''     ]); end
  if f==2             fprintf(ofp, ['%-',unt_f,'s',delim,''], [tablefields{f},' (*1)']); end
  if f==3 || f==4     fprintf(ofp, ['%', num_f,'s',delim,''], [tablefields{f},' (*2)']); end
  if f>4 && f<Nstat+2 fprintf(ofp, ['%', num_f,'s',delim,''], tablefields{f}); end
  if f==Nstat+2       fprintf(ofp, ['%', num_s_f,'s\n'],      [tablefields{f},' (*3)']); end
end

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','-'  );
end
fprintf(ofp, '\n');

%-----Raw lmc-fields----------------------------------------------

for f = 1:length(fn)

   var = getfield(lmc, fn{f});

   if ~(strcmp(fn{f},'mwRTM_vegcls') || strcmp(fn{f},'mwRTM_soilcls'))

    fprintf(ofp, tableformat, ...
            fn{f}, units{f}, ...
            stats_array(var,weights, Nstat));

   else

     tmp_stats      = stats_array(var,weights, Nstat);
     tmp_stats(1:2) = nodata_val;

     fprintf(ofp, tableformat, ...
            fn{f}, units{f}, ...
            tmp_stats);

   end

end

%-----Footnotes---------------------------------------------------

for f = 1:str_l+1+unt_l+Nstat+(Nstat-1)*num_l+num_s_l
  fprintf(ofp, '%s','-'  );
end

fprintf(ofp, '\n%s\n\n', ...
        'See SMAP L4_SM Data Products Specification Document for additional information.');

fprintf(ofp, '%s\n', ...
        '(*1) Units are valid for all statistics except N [dimensionless].');

fprintf(ofp, '%s\n',...
        '(*2) Mean and std-dev statistics are weighted by the land fraction of each grid cell.');

fprintf(ofp, '%s\n', ...
        '(*3) N is the number of 9 km EASEv2 grid cells that contribute to the statistics.');

fclose(ofp);


% *********************************************************************************************
% *********************************************************************************************
% *********************************************************************************************

function [array] = map_M09toM36( array, M36_row_col, return_mean, check_on )

% Reduce M36 observations, partitioned to M09 in the aup-file,
% to the actual M36 resolution.
%
% Check if the partitioned M09 obs are identical (stdv ~ 0 K)
% inside their corresponding M36 grid cell.
%
% Input: array:       1-dimensional
%        M36_row_col: colum 1= M36 row-ind, columns = M36 col-ind
%        return_mean: provide average of all M09 values inside a M36 grid cell
%        check_on:    either stop or send a warning if all
%                     elements inside one M36 grid cell are not identical
%
%
% GDL, 24 Jan 2014
% GDL,  1 Feb 2014: Option to return the average M09 value across the M36 pixel,
%                   rather than one of the values.
%                   This is needed to find the land_fraction across a M36
%                   gridcell.
%
% ---------------------------------------------------------------------

if length(array) ~= size(M36_row_col,1)
  error('input arrays need to have the same length');
end

stdv_Tbobs_tol      = 1E-2; %stdv of M09 obs within M36 grid cells

[unique_rc, ia, ic] = unique(M36_row_col,'rows');

if (check_on > 0 || return_mean)

    if return_mean
      array_out           = NaN*zeros(length(ia),1);
    end

    N_bad = 0;
    for i=1:size(unique_rc,1)

        if check_on > 0
        if (std(array( ...
                 M36_row_col(:,1) == unique_rc(i,1) & ...
                 M36_row_col(:,2) == unique_rc(i,2) ))...
            > stdv_Tbobs_tol)
            N_bad = N_bad+1;
        end
        end

        if return_mean
           array_out(i) = mean(array( ...
                 M36_row_col(:,1) == unique_rc(i,1) & ...
                 M36_row_col(:,2) == unique_rc(i,2) ));
        end

    end

    if N_bad > 0
        if check_on == 2
            error([num2str(N_bad),' M36 grid cells out of ',...
               num2str(length(unique_rc)),' contain M09 obs with std-dev<>0'])
        else
            disp(['WARNING: ',num2str(N_bad),' M36 grid cells out of ',...
               num2str(length(unique_rc)),' contain M09 obs with std-dev<>0'])
        end
    end

end

if ~return_mean
   array = array(ia);
else
   array = array_out;
end

% *********************************************************************************************

function [ stats_out ] = stats_array( array_in, weights, Nstat )

% calculate elementary summary statistics for array_in
%
% input:  array_in = numerical array (1-dimensional)
%         weights  = weight for each element in array_in (1-dimensional)
%
% output: array_out = [mean stdv min max N_data frac]
%
% GDL, 16 Jan 2014
% GDL, 30 Jan 2014: remove NaN prior to the calculation of stats
% GDL,  1 Feb 2014: - weighted statistics for mean and stdv
%                   - add 6th statistic (frac):
%                     the actual used fraction of N_data (based on weights)
% GDL,  7 Feb 2014: added Nstat as input, to limit the returned output.
%                   should be revised in the future to select specific output,
%                   rather than to limit the output fields
% ---------------------------------------------------------------------

nodata_val   = -9999;
nodata_tol   = 1E-5;

stats_out    = nodata_val+zeros(6,1);
stats_out(5) = 0 ;
stats_out(6) = 0 ;

if length(weights) ~= length(array_in)
  error('input vectors need to have the same dimensions')
end

nodata_ind  = (isnan(array_in) | ...
               abs(array_in - nodata_val)<nodata_tol | ...
               isnan(weights) | ...
               abs(weights - nodata_val)<nodata_tol);

array_in   = array_in(~nodata_ind);
weights    = weights(~nodata_ind);

if ~isempty(array_in)

  %actual fraction of data points used in stats
  stats_out(6) = mean(weights);

  weights    = weights./sum(weights);

  %-------------

  %weighted mean
  stats_out(1) = sum(weights.* array_in);

  %weighted stdv, without statistical bias correction
  stats_out(2) = sqrt(sum(weights.*(array_in-stats_out(1)).^2));

  %min and max
  stats_out(3) = min(array_in);
  stats_out(4) = max(array_in);

  %total number of data points
  stats_out(5) = length(array_in);

end

%only return the first Nstat statistics
stats_out = stats_out(1:Nstat);

% ========= EOF =========================================================
