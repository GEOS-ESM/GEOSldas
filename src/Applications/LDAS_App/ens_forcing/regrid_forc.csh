#!/usr/bin/csh -f 

set echo 
setenv MYNAME regrid_forc.csh

if ( $#argv < 3 ) then
   echo " "
   echo " SYNOPSIS "
   echo " "
   echo "  $MYNAME  force_in force_rgd outgrid"
   echo " "
   echo " where"
   echo "   force_in  -  path to  orginal forcing, e.g., $ADAS_EXPDIR/atmens/"
   echo "   force_rgd  -  path to regrid forcing, e.g., $ADAS_EXPDIR/atmens/rgdlfo"
   echo "   outgrid  - output grid , e.g., PE180x1080-CF ") 
   echo "   NENS  ensemble number , e.g., 24 should be env var ")
   exit(0)
endif 

set  ogdir = $1 
set  rgdir = $2
set  outgrid = $3 

set NENSin = $NENS 
echo "check NENS, $NENSin "

mkdir $ogdir/tmpd

@ inens = 0
while ($inens < $NENSin)
      @ inens ++
      if ($inens <10) then
          set ENSDIR = `echo mem00${inens}`
      else if($inens<100) then
          set ENSDIR=`echo mem0${inens}`
      endif  
      /bin/ln -s $ogdir/ensdiag/${ENSDIR} $ogdir/tmpd/${ENSDIR}
      mkdir $rgdir/${ENSDIR} 
end

cd $ogdir
/bin/ls -1 tmpd/mem*/*inst*lfo*nc4| awk 'NR == 1 {printf $0} NR > 1 {printf ", %s",$0} END {printf "\n"}' > tmpf
cat tmpf | sed 's/, /,/g' > tmpff
set ininst1 = (`cat tmpff  `)
set outinst1 = (`cat tmpff | sed 's/tmpd/rgdlfo/g' `)

/bin/ls -1 tmpd/mem*/*tavg*lfo*nc4| awk 'NR == 1 {printf $0} NR > 1 {printf ", %s",$0} END {printf "\n"}' > tmpf
cat tmpf | sed 's/, /,/g' > tmpff
set intavg1 = (`cat tmpff `)
set outtavg1 = (`cat tmpff | sed 's/tmpd/rgdlfo/g' `)

set mympi = "mpirun"

   $mympi -np 24 $GEOSBIN/Regrid_Util.x  -i  $ininst1 \
                        -o $outinst1 \
                        -nx 4 -ny 6  \
                        -ogrid $outgrid
  
   $mympi -np 24 $GEOSBIN/Regrid_Util.x  -i  $intavg1 \
                        -o $outtavg1 \
                        -nx 4 -ny 6  \
                        -ogrid $outgrid
 
