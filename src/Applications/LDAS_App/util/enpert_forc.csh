#!/usr/bin/csh -f 

## env : ADAS_EXPDIR,  GRID,  NENS, RUNDIR

set force_cntr =  "${ADAS_EXPDIR}/recycle/holdforc"
set force_orig =  "${ADAS_EXPDIR}/atmens"   
set force_rgd = "${ADAS_EXPDIR}/atmens/rgdlfo" 
set outgrid =  "${GRID}"
mkdir $force_rgd

$RUNDIR/regrid_forc.csh  $force_orig $force_rgd $outgrid    
rm -rf  $force_orig/tmp* 
python  $RUNDIR/ensemble_forc.py $force_rgd $force_cntr $NENS 
 
