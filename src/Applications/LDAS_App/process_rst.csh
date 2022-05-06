#!/bin/csh -f

setenv  SPONSORID    $1
setenv  EXPID        $2
setenv  out_dir      $3
setenv  out_bcsdir   $4
setenv  out_tilefile $5
setenv  model        $6
setenv  HAVE_RESTART $7
setenv  YYYYMMDDHH   $8 
setenv  in_rstfile   $9
setenv  in_tilefile  $10
setenv  surflay      $11
setenv  in_wemin     $12
setenv  out_wemin    $13

set PWD = `pwd`
setenv INSTDIR `echo $PWD | rev | cut -d'/' -f2- | rev`
set Bin = $PWD

set YYYYMMDD = `echo $YYYYMMDDHH | cut -c1-8`

set YYYY = `echo $YYYYMMDD | cut -c1-4`
set   MM = `echo $YYYYMMDD | cut -c5-6`
set   DD = `echo $YYYYMMDD | cut -c7-8`

mkdir -p $out_dir/InData/
cd $out_dir

if ($HAVE_RESTART == M) then
  set ARCDIR = /archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/d5124_m2_jan79/rs/Y ; set mlable = jan79
  if ($YYYY > 1991) set ARCDIR = /archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/d5124_m2_jan91/rs/Y ; set mlable = jan91
  if ($YYYY > 2000) set ARCDIR = /archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/d5124_m2_jan00/rs/Y ; set mlable = jan00
  if ($YYYY > 2010) set ARCDIR = /archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/d5124_m2_jan10/rs/Y ; set mlable = jan10	
  set rstfile = ${ARCDIR}${YYYY}/M${MM}/d5124_m2_${mlable}.${MODEL}_internal_rst.${YYYYMMDD}_21z.bin
  dmget $rstfile
  set INTILFILE = /gpfsm/dnb02/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_MERRA-2/CF0180x6C_DE1440xPE0720/CF0180x6C_DE1440xPE0720-Pfafstetter.til
  set  WEMIN_IN = 26
  /bin/cp -p $rstfile $out_dir/InData/M2Restart
  setenv in_wemin $WEMIN_IN
  setenv in_tilefile $INTILEFILE
  setenv in_rstfile $out_dir/InData/M2Restart
endif

if ($HAVE_RESTART == G) then
  set rstfile = `echo $RESTART_PATH | rev | cut -c 2- | rev`
  set INTILFILE = `readlink $RESTART_ID/scratch/tile.data`
  set  WEMIN_IN = 13
  if ( `echo $INTILFILE | grep -n 'NLv3'` == '' ) set  WEMIN_IN = 26
  /bin/cp -p $rstfile $out_dir/InData/M2Restart
  setenv in_wemin $WEMIN_IN
  setenv in_tilefile $INTILEFILE
  setenv in_rstfile $out_dir/InData/M2Restart
endif

if ($HAVE_RESTART == F) then
   set date_16 = `date -d"2017-1-24" +%Y%m%d`
   set date_17 = `date -d"2017-11-1" +%Y%m%d`
   set date_21 = `date -d"2018-7-11" +%Y%m%d`
   set date_22 = `date -d"2019-3-13" +%Y%m%d`
   set date_25 = `date -d"2020-1-30" +%Y%m%d`
   set expdate = `date -d"$YYYY-$MM-$DD" +%Y%m%d`

   if ($expdate < $date_16) then
      echo "WARNING : FP restarts before $date_16 are not availale."
      echo "          Please select RESTART: M and use MERRA-2, instead."
      exit
   endif

   if (($expdate >= $date_16) && ($expdate < $date_17)) then
      set fpver = GEOS-5.16/GEOSadas-5_16/
      set fplab = f516_fp
      set INTILFILE = /discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_Ostia/CF0720x6C_DE2880xPE1440/CF0720x6C_DE2880xPE1440-Pfafstetter.til
      set  WEMIN_IN = 26
      set   rstfile = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.${MODEL}_internal_rst.${YYYYMMDD}_21z.bin
      dmget $rstfile
      /bin/cp -p $rstfile $out_dir/InData/M2Restart
   endif
   if (($expdate >= $date_17) && ($expdate < $date_21)) then
      set fpver = GEOS-5.17/GEOSadas-5_17/
      set fplab = f517_fp
      set INTILFILE = /discover/nobackup/ltakacs/bcs/Icarus/Icarus_Ostia/CF0720x6C_CF0720x6C/CF0720x6C_CF0720x6C-Pfafstetter.til
      set  WEMIN_IN = 26
      set TARFILE = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.rst.${YYYYMMDD}_21z.tar
      dmget $TARFILE
      set   rstfile = ${fplab}.${MODEL}_internal_rst.${YYYYMMDD}_21z.nc4
      tar -xvf $TARFILE $rstfile && /bin/mv $rstfile $out_dir/InData/M2Restart	 	    
   endif
   if (($expdate >= $date_21) && ($expdate < $date_22)) then
      set fpver = GEOS-5.21/GEOSadas-5_21/
      set fplab = f521_fp
      set INTILFILE = /discover/nobackup/ltakacs/bcs/Icarus/Icarus_Ostia/CF0720x6C_CF0720x6C/CF0720x6C_CF0720x6C-Pfafstetter.til
      set  WEMIN_IN = 26
      set TARFILE = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.rst.${YYYYMMDD}_21z.tar
      dmget $TARFILE
      set   rstfile = ${fplab}.${MODEL}_internal_rst.${YYYYMMDD}_21z.nc4
      tar -xvf $TARFILE $rstfile && /bin/mv $rstfile $out_dir/InData/M2Restart	 	    
   endif
   if (($expdate >= $date_22) && ($expdate < $date_25)) then
      set fpver = GEOS-5.22/GEOSadas-5_22/
      set fplab = f522_fpp
      set INTILFILE = /discover/nobackup/ltakacs/bcs/Icarus/Icarus_Ostia/CF0720x6C_CF0720x6C/CF0720x6C_CF0720x6C-Pfafstetter.til
      set  WEMIN_IN = 26
      set TARFILE = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.rst.${YYYYMMDD}_21z.tar
      dmget $TARFILE
      set   rstfile = ${fplab}.${MODEL}_internal_rst.${YYYYMMDD}_21z.nc4
      tar -xvf $TARFILE $rstfile && /bin/mv $rstfile $out_dir/InData/M2Restart	 	    
  endif
  if ($expdate >= $date_25) then
     set fpver = GEOS-5.25/GEOSadas-5_25/
     set fplab = f525land_fpp
     set INTILFILE = /discover/nobackup/ltakacs/bcs/Icarus-NLv3/Icarus-NLv3_Ostia/CF0720x6C_CF0720x6C/CF0720x6C_CF0720x6C-Pfafstetter.til
     set  WEMIN_IN = 13
     set TARFILE = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.rst.${YYYYMMDD}_21z.tar
     dmget $TARFILE
     set   rstfile = ${fplab}.${MODEL}_internal_rst.${YYYYMMDD}_21z.nc4
     tar -xvf $TARFILE $rstfile && /bin/mv $rstfile $out_dir/InData/M2Restart
  endif
  setenv in_wemin $WEMIN_IN
  setenv in_tilefile $INTILEFILE
  setenv in_rstfile $out_dir/InData/M2Restart
endif
 
echo ' '

cat << _EOI3_ > mkLDASsa.j
#!/bin/csh -fx
 
#SBATCH --account=${SPONSORID}
#SBATCH --time=1:00:00
#SBATCH --ntasks=56
#SBATCH --job-name=mkLDAS
#SBATCH --constraint=hasw
#SBATCH --qos=debug
#SBATCH --output=mkLDAS.o
#SBATCH --error=mkLDAS.e
 
source $INSTDIR/bin/g5_modules
setenv OMPI_MCA_shmem_mmap_enable_nfs_warning 0
limit stacksize unlimited
set echo

${Bin}/esma_mpirun -np 56 ${Bin}/mk_catchANDcnRestarts.x -model ${model} -time ${YYYYMMDDHH} -in_tilefile ${in_tilefile}  -out_bcs ${out_bcsdir} -out_tilefile ${out_tilefile} -out_dir ${out_dir} -surflay ${surflay} -in_wemin ${in_wemin} -out_wemin ${out_wemin} -in_rst ${in_rstfile} -out_rst ${model}_internal_rst.$YYYYMMDD  

echo DONE > done_rst_file

_EOI3_

sbatch mkLDASsa.j
cd $PWD
