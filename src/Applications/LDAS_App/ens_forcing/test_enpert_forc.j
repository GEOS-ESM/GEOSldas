#!/usr/bin/csh

# sample script for testing "enpert_forc.csh"
#
# must edit paths before submitting this script
#
# -----------------------------------------------------------

#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --job-name=test_enpert_forc
#SBATCH --output=/discover/nobackup/[USER]/test_enpert_forc.o%j.txt
#SBATCH --export=NONE
#SBATCH --qos=debug

# set environment variables needed by enpert_forc.csh

setenv GEOSBIN     /discover/nobackup/[USER]/GEOSldas/install/bin  
setenv ADAS_EXPDIR /discover/nobackup/[USER]/[ADAS_EXPID]/
setenv NENS        24  
setenv GRID        PE180x1080-CF

# load modules

source $GEOSBIN/g5_modules

# python should come with ESMA_env g5_modules
# module load python/GEOSpyD/Ana2019.03_py3.7

if ( -e /etc/os-release ) then
  module load nco/4.8.1
else
  module load other/nco-4.6.8-gcc-5.3-sp3
endif

# execute test

$GEOSBIN/enpert_forc.csh 

# ====================== EOF ===================================
