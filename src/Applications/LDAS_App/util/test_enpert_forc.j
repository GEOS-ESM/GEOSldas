#!/usr/bin/csh
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --job-name=regd
#SBATCH --output=/discover/nobackup/qzhang/tmp/testregrid/regridUt.o%j.txt
#SBATCH --export=NONE
#SBATCH --qos=debug

setenv  GEOSBIN /discover/nobackup/projects/c612/user/qzhang/git_check/ensforc/GEOSldas/install/bin  
source $GEOSBIN/g5_modules
module load python/GEOSpyD/Ana2019.03_py3.7
if ( -e /etc/os-release ) then
  module load nco/4.8.1
else
  module load other/nco-4.6.8-gcc-5.3-sp3
endif

setenv RUNDIR /discover/nobackup/qzhang/tmp/testregrid  
setenv ADAS_EXPDIR /discover/nobackup/projects/c612/user/qzhang/LA5294tag 
setenv NENS 24  
setenv GRID PE180x1080-CF

$RUNDIR/enpert_forc.csh 

