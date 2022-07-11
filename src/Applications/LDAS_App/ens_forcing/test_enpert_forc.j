#!/usr/bin/csh
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --job-name=regd
#SBATCH --output=/discover/nobackup/user/test/test.o%j.txt
#SBATCH --export=NONE
#SBATCH --qos=debug

########################################################
## user:   user
## GEOSldas: GEOSldas ROOTdir 
## testdir:  user testdir
## GEOSadas: GEOSadas ROOTdir 
########################################################

setenv  GEOSBIN /discover/nobackup/user/GEOSldas/install/bin  
source $GEOSBIN/g5_modules
module load python/GEOSpyD/Ana2019.03_py3.7
if ( -e /etc/os-release ) then
  module load nco/4.8.1
else
  module load other/nco-4.6.8-gcc-5.3-sp3
endif

setenv RUNDIR /discover/nobackup/user/testdir/ 
setenv ADAS_EXPDIR /discover/nobackup/user/GEOSadas 
setenv NENS 24  
setenv GRID PE180x1080-CF

$GEOSBIN/enpert_forc.csh 

