import os
import sys
import subprocess
import shutil
import glob
import ruamel.yaml
import shlex
from remap_utils import *

def remap_config_ldas(config, RESTART_str, RESTART_PATH, RESTART_ID):
   yyyymmddhh = config['input']['shared']['yyyymmddhh']
   MODEL      = config['input']['surface']['catch_model']
   out_dir    = config['output']['shared']['out_dir']
   # MERRA2
   if RESTART_str == "M" :
     config['input']['shared'] = merra2_expid(config['input']['shared'])
     config['input']['shared']['rst_dir'] = out_dir+ '/merra2_tmp_'+ yyyymmddhh
     config['input']['surface']['wemin'] = 26
     config['input']['shared']['bcs_dir'] = '/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Ganymed-4_0/Ganymed-4_0_MERRA-2/CF0180x6C_DE1440xPE0720/'
     
   if RESTART_str == "G" :
     # WY note: it is a bad idea to overload restart_path and restart_id
     config['input']['surface']['catch_tilefile'] = os.path.realpath(RESTART_ID+'scratch/tile.data')
     config['input']['shared']['rst_dir']  = os.path.dirname(RESTART_PATH) + '/'
     config['input']['surface']['wemin'] = 13
     if 'NLv' in config['input']['surface']['catch_tilefile'] : config['input']['surface']['wemin'] = 26
     
   if RESTART_str == "F" :
     date_16 = 20170124 
     date_17 = 20171101 
     date_21 = 20180711 
     date_22 = 20190313
     date_25 = 20200130
     date_27 = 20210225
     date_29 = 20220228
     expdata = int(yyyymmddhh[0:8])
     if (expdate < date_16):
       print( "WARNING : FP restarts before $date_16 are not availale.")
       print( "          Please select RESTART: M and use MERRA-2, instead.")
       sys.exit(1)

     config['input']['shared']['bcs_dir'] = '/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus/Icarus_Ostia/CF0720x6C_CF0720x6C/'
     config['input']['surface']['wemin'] = 26
     config['input']['shared']['rst_dir'] = out_dir+'/InData'+ '/'
     suffix = '_21z.tar'

     if ((date_16 <= expdate) and (expdate < date_17)):
        fpver = 'GEOS-5.16/GEOSadas-5_16/'
        fplab = 'f516_fp'
        config['input']['shared']['bcs_dir'] = '/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Ganymed-4_0/Ganymed-4_0_Ostia/CF0720x6C_DE2880xPE1440/'
        suffix = '_21z.bin'

     if ((date_17 <= expdate) and (expdate < date_21)):
        fpver = 'GEOS-5.17/GEOSadas-5_17/'
        fplab = 'f517_fp'

     if ((date_21 <= expdate) and (expdate < date_22)):
        fpver = 'GEOS-5.21/GEOSadas-5_21/'
        fplab = 'f521_fp'

     if ((date_22 <= expdate) and (expdate < date_25)):
        fpver = 'GEOS-5.22/GEOSadas-5_22/'
        fplab = 'f522_fp'

     if (date_25 <= expdate and expdate < date_27) :
        fpver = 'GEOS-5.25/GEOSadas-5_25/'
        fplab = 'f525land_fpp'
        config['input']['surface']['wemin'] = 13

     if (date_27 <= expdate and expdate < date_29) :
        fpver = 'GEOS-5.27/GEOSadas-5_27/'
        fplab = 'f5271_fpp'
        config['input']['surface']['wemin'] = 13
     
     if (date_29 <= expdate ) :
        fpver = 'GEOS-5.29/GEOSadas-5_29/'
        fplab = 'f5293_fpp'
        config['input']['surface']['wemin'] = 13


     rstfile = '/archive/u/dao_ops/'+fpver+'/'+fplab+'/rs/Y'+YYYY+'/M'+MM+'/'+fplab+'.'+MODEL+'_internal_rst.'+YYYYMMDD + suffix
     fname = os.path.basename(rstfile)
     dest = out_dir+'/InData'+ '/'+fname
     print("Copy file "+ rstfile +" to " + out_dir+'/InData'+ '/')
     shutil.copy(rstfile, dest)

     if suffix == "_21z.tar" :
        new_rst = out_dir+'/InData'+ '/'+ fplab+'.'+MODEL+'_internal_rst.'+YYYYMMDD+'_21z.nc4'
        subprocess.call(['tar', '-xvf', dest, new_rst])

   return config
