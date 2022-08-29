#!/usr/bin/csh -f 

# Script to create re-centered ensemble land forcing files (*inst1_2d_lfo*, *tavg1_2d_lfo*),
# typically used in LADAS:
#
#   1. Regrid coarse-resolution lfo files from the ADAS atmospheric ensemble to the higher 
#      resolution of the single-member central (or deterministic) simulation.
#
#   2. Compute perturbations from regridded lfo files.
#
#   3. Apply perturbations to lfo files from central simulation. 
#
# Requires environment variables: 
#   ADAS_EXPDIR - path the ADAS experiment directory (ensemble of lfo files)
#   GRID        - target grid of forcing files (typically that of deterministic ADAS simulation)
#   NENS        - number of (LDAS) ensemble members (<= number of ADAS ensemble members)
#   GEOSBIN     - GEOSldas "ROOT" directory (typically GEOSldas/install/bin) 
#
# Input data sets:
#   1. lfo files from coarse-resolution ADAS ensemble
#   2. lfo files from higher-resolution deterministic ADAS simulation
#
# Output data set:
#   1. ensemble of lfo files at resolution of deterministic ADAS simulation
#
# ------------------------------------------------------------------------------------

set force_cntr = "${ADAS_EXPDIR}/recycle/holdforc"
set force_orig = "${ADAS_EXPDIR}/atmens"   
set force_rgd  = "${ADAS_EXPDIR}/atmens/rgdlfo" 
set outgrid    = "${GRID}"

mkdir   $force_rgd

$GEOSBIN/regrid_forc.csh  $force_orig $force_rgd $outgrid    

rm -rf  $force_orig/tmp* 

python  $GEOSBIN/ensemble_forc.py $force_rgd $force_cntr $NENS 

cd $force_rgd
@ inens = 0
while ($inens < $NENS)
      @ inens ++
      if ($inens <10) then
          set ENSDIR = `echo mem00${inens}`
      else if($inens<100) then
          set ENSDIR=`echo mem0${inens}`
      endif
      cd ${ENSDIR}
      /bin/rm -rf  *lfo*nc4
      $GEOSBIN/stripname nc4.Cpert nc4
      cd $force_rgd
end

