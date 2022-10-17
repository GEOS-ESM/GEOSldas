%=========================================================================
% script to convert SMAPL4 SM climatology statistic binary output to nc4
% file. 
%
% The tile_ids are in the same order as in ldas_tilecoord.txt 
%
% Q. Liu       - Jun 20, 2016
% G. De Lannoy - Apr 21, 2015
%=========================================================================

clear all

write_GEOSldas_tileorder = 0;

time_stamp = '2001_p1_2021_p73';
specs_tag  = '_W_5p_Nmin_210';

bin_fpath  = '/home/qliu/smap/SMAP_Nature/SMAP_Nature_v10/SMAP_Nature_v10.0/output/SMAP_EASEv2_M09_GLOBAL/stats/cli/';

file_out   = ['L4SM_NRv10.0_cli_',time_stamp,specs_tag,'_sm_wetness_EASEv2_M09.nc4'];

if write_GEOSldas_tileorder
   file_out = ['GEOSldas_',file_out];
end

fname_out  = [bin_fpath, file_out]; 

% read an arbitrary file to get header info
fname = [bin_fpath, 'cli_',time_stamp,specs_tag,'_sm_rootzone_wet_p06.bin'];

N_stat = 99 + 5;   % 99 percentile values + 5 stats (mean, stdv,min, max, N_data) 

[data_tmp, tile_id, lon, lat] = ...
    read_seqbin_clim_pctl_file(fname, 1, N_stat, 'latlon_id','cli');

N_tile = length(tile_id);

if write_GEOSldas_tileorder

    [tile_id_out, tile_idx] = sort(tile_id,'ascend');
    lon_out = lon(tile_idx);
    lat_out = lat(tile_idx);
else
    tile_id_out = tile_id;
    lon_out = lon;
    lat_out = lat;
end

clear data_tmp

netcdf.setDefaultFormat('FORMAT_NETCDF4');
fout_id = netcdf.create(fname_out, 'NETCDF4');

if fout_id < 0, error(['Creating ' fname_out 'failed']); end

% Setup global attributes

NC_GLOBAL = netcdf.getConstant('GLOBAL');

netcdf.putAtt(fout_id, NC_GLOBAL, 'Title', ['SMAP L4 SM pentad clim. ',time_stamp,' statistics on EASEv2 M09']);

netcdf.putAtt(fout_id, NC_GLOBAL, 'Filename', file_out);
netcdf.putAtt(fout_id, NC_GLOBAL, 'Institution', 'NASA GMAO');
netcdf.putAtt(fout_id, NC_GLOBAL, 'History', ['File written by matlab-r2021a on ',datestr(now)]);
netcdf.putAtt(fout_id, NC_GLOBAL, 'Contact', 'NASA/GMAO Rolf Reichle');
netcdf.putAtt(fout_id, NC_GLOBAL, 'Comments', 'NETCDF-4');

% Define dimensions:

dimid1 = netcdf.defDim(fout_id,'tile',length(tile_id_out));
dimid2 = netcdf.defDim(fout_id,'percentile_wetness' , 99);
dimid3 = netcdf.defDim(fout_id,'pentad', netcdf.getConstant('UNLIMITED'));

% Define global variables:

varid = netcdf.defVar(fout_id,'tile_id','int',dimid1);
netcdf.putAtt(fout_id,varid,'standard_name','tile_id');
netcdf.putAtt(fout_id,varid,'long_name','tile_id');
netcdf.putAtt(fout_id,varid,'units','-');
netcdf.putVar(fout_id,varid,tile_id_out);

varid = netcdf.defVar(fout_id,'lon','double',dimid1);
netcdf.putAtt(fout_id,varid,'standard_name','longitude');
netcdf.putAtt(fout_id,varid,'long_name','longitude');
netcdf.putAtt(fout_id,varid,'units','degrees_east');
netcdf.putVar(fout_id,varid,lon_out);

varid = netcdf.defVar(fout_id,'lat','double',dimid1);
netcdf.putAtt(fout_id,varid,'standard_name','latitude');
netcdf.putAtt(fout_id,varid,'long_name','latitude');
netcdf.putAtt(fout_id,varid,'units','degrees_north');
netcdf.putVar(fout_id,varid,lat_out);

% Synchronize global
netcdf.sync(fout_id)

% Define groups and variables in each group

vars = {'mean', 'stdv', 'min', 'max', 'N_data', 'percentile_UL'};

da_group_id(1) = netcdf.defGrp(fout_id, 'rootzone_wetness_cli_stat');

da_group_id(2) = netcdf.defGrp(fout_id, 'profile_wetness_cli_stat');

fillValue = single(1.e15);

DeflateLevel = 5;

% Put data into the data group:

% Insert data:

for i=1:length(da_group_id) %loop through groups
    
    start_time = 0;
    
    for k = 1:73 % loop through pentad
        
        clear data
        
        if i==1
            
            fname = [bin_fpath, 'cli_',time_stamp,specs_tag,'_sm_rootzone_wet_p', ...
                num2str(k, '%2.2d'), '.bin'];
            
        else
            
            fname = [bin_fpath, 'cli_',time_stamp,specs_tag,'_sm_profile_wet_p', ...
                num2str(k, '%2.2d'), '.bin'];
            
        end
        
        [data_tmp, tile_id_old, lon_old, lat_old] = ...
            read_seqbin_clim_pctl_file(fname, 1, 104, 'latlon_id','cli');

        clear tile_id_old lon_old lat_old

        if write_GEOSldas_tileorder
           data = data_tmp(tile_idx,:);
        else
           data = data_tmp;
        end

        clear data_tmp 
        
        for iv = [1 3 4 6:size(data,2)]
            data(data(:,iv)>1, iv) = 1.;
        end
               
        
        data(data<-9998) = fillValue;
        data(isnan(data)) = fillValue;
        
        
        for iv = 1:4 % loop through variables
            
            if k == 1 

                varid(iv) = netcdf.defVar(da_group_id(i),vars{iv},'float',[dimid1,dimid3]);
                netcdf.putAtt(da_group_id(i),varid(iv),'name', vars{iv});
                netcdf.putAtt(da_group_id(i),varid(iv),'units','wetness');
                netcdf.defVarFill(da_group_id(i),varid(iv),false,fillValue); 
                netcdf.defVarDeflate(da_group_id(i),varid(iv),true,true, DeflateLevel);

            end
            
            netcdf.putVar(da_group_id(i),varid(iv),[0,start_time], [N_tile,1], squeeze(data(:,iv)));
        
        end
        
        if k == 1

            varid(5) = netcdf.defVar(da_group_id(i),vars{5},'int',[dimid1,dimid3]);
            netcdf.putAtt(da_group_id(i),varid(5),'name', vars{5});
            netcdf.putAtt(da_group_id(i),varid(5),'units','-');
            netcdf.defVarDeflate(da_group_id(i),varid(5),true,true, DeflateLevel);

        end

        netcdf.putVar(da_group_id(i),varid(5),[0,start_time], [N_tile,1],squeeze(data(:,5)));
        
        if k==1

            varid(6) = netcdf.defVar(da_group_id(i),vars{6},'float',[dimid2,dimid1,dimid3]);
            netcdf.putAtt(da_group_id(i),varid(6),'name', vars{6});
            netcdf.putAtt(da_group_id(i),varid(6),'units','wetness');
            netcdf.defVarFill(da_group_id(i),varid(6),false,fillValue);
            netcdf.defVarDeflate(da_group_id(i),varid(6),true,true, DeflateLevel);

        end
        
        netcdf.putVar(da_group_id(i),varid(6),[0,0,start_time], [99,N_tile, 1],squeeze(data(:,6:104)'));
        
        start_time = start_time +1;
    end
    
    netcdf.sync(da_group_id(i))
    
    
end

% Synchronize:
netcdf.close(fout_id)

%====================================EOF===================================



