#!/bin/csh -f

setenv  SPONSORID    $1
setenv  EXPID        $2
setenv  EXPDIR       $3
setenv  BCSDIR       $4
setenv  TILFILE      $5
setenv  LSM_CHOICE   $6
setenv  HAVE_RESTART $7
setenv  YYYYMMDD     $8 
setenv  RESTART_ID   $9
setenv  RESTART_DOMAIN $10
setenv  RESTART_PATH $11
setenv  NUMENS       $12
setenv  RUN_IRRIG    $13
setenv  SURFLAY      $14
setenv  WEMIN_IN     $15
setenv  WEMIN_OUT    $16
setenv  RESTART_short ${RESTART_PATH}/${RESTART_ID}/output/${RESTART_DOMAIN}/
setenv  PARAM_FILE `ls $RESTART_short/rc_out/*/*/*ldas_catparam* | head -1`

set PWD=`pwd`
setenv INSTDIR `echo $PWD | rev | cut -d'/' -f2- | rev`

if($LSM_CHOICE == 1) then 
   set MODEL=catch
endif

if($LSM_CHOICE == 2) then
   set MODEL=catchcn
endif

switch ($HAVE_RESTART) 

case [0] :

   echo 
   echo "########################################################"
   echo 
   echo "          NOTE:   IMPORTANT IMPORTANT IMPROTANT          "
   echo "   							  "
   echo "           In the absence of a restart file, 		  "
   echo "      a global catch(cn)_internal_rst file is being	  "
   echo "       produced by regridding SMAP_Nature_v05 restarts	  "
   echo "      from SMAP 9km grid to the BCS grid. The start date "
   echo "              will be forced January 1st.		  "
   echo "							  "
   echo "               PLEASE SPIN UP BEFORE USE		  "
   echo
   echo "########################################################" 
   echo

## No restarts : regrid from archived SMAP M09 restarts

    mkdir -p $EXPDIR/$EXPID/mk_restarts/OutData1/
    mkdir -p $EXPDIR/$EXPID/mk_restarts/OutData2/
    ln -s $BCSDIR/$TILFILE $EXPDIR/$EXPID/mk_restarts/OutData1/OutTileFile
    ln -s $BCSDIR/$TILFILE $EXPDIR/$EXPID/mk_restarts/OutData2/OutTileFile
    ln -s $BCSDIR/clsm $EXPDIR/$EXPID/mk_restarts/OutData2/clsm
    ln -s $INSTDIR/bin $EXPDIR/$EXPID/mk_restarts/

    cd $EXPDIR/$EXPID/mk_restarts/

    cat << _EOI_ > mkLDASsa.j
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
if ( -e /etc/os-release ) then
  module load nco/4.8.1
else
  module load other/nco-4.6.8-gcc-5.3-sp3 
endif
#setenv OMPI_MCA_shmem_mmap_enable_nfs_warning 0
#setenv MKL_CBWR SSE4_2 # ensure zero-diff across archs
#setenv MV2_ON_DEMAND_THRESHOLD 8192 # MVAPICH2
setenv LAIFILE `find ${BCSDIR}/lai_clim*`
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BASEDIR}/Linux/lib

limit stacksize unlimited
 
$INSTDIR/bin/esma_mpirun -np 56 bin/mk_GEOSldasRestarts.x -a ${SPONSORID} -b ${BCSDIR} -t ${TILFILE} -m ${MODEL} -s ${SURFLAY} -j Y

sleep 3

if($LSM_CHOICE == 1) then
   /bin/cp OutData1/catch_internal_rst OutData2/catch_internal_rst
else
   /bin/cp OutData1/catchcn_internal_rst OutData2/catchcn_internal_rst
endif

$INSTDIR/bin/esma_mpirun -np 56 bin/mk_GEOSldasRestarts.x -a ${SPONSORID} -b ${BCSDIR} -t ${TILFILE} -m ${MODEL} -s ${SURFLAY} -j Y

_EOI_

    if($LSM_CHOICE == 1) sed -i '$ a\bin/Scale_Catch OutData1/catch_internal_rst OutData2/catch_internal_rst catch_internal_rst $SURFLAY $WEMIN_IN $WEMIN_OUT \'  mkLDASsa.j
    if($LSM_CHOICE == 2) sed -i '$ a\bin/Scale_CatchCN OutData1/catchcn_internal_rst OutData2/catchcn_internal_rst catchcn_internal_rst $SURFLAY $WEMIN_IN $WEMIN_OUT \'  mkLDASsa.j

    sed -i '$ a\ \'  mkLDASsa.j
    sed -i '$ a\## Done creating catch*_internal_rst file \'  mkLDASsa.j
    sed -i '$ a\ \'  mkLDASsa.j

    cat << _EOI2_ > mkLDASsa.j2
sleep 2

if($LSM_CHOICE == 1) then
   if (-f irrigation_internal_rst && $RUN_IRRIG == 1) then 
      ncks -4  -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF irrigation_internal_rst -A catch_internal_rst
   endif
   ln -s  catch_internal_rst catch_internal_rst.$YYYYMMDD
else
   if (-f irrigation_internal_rst && $RUN_IRRIG == 1) then 
      ncks -4 -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF irrigation_internal_rst -A catchcn_internal_rst 
   endif
   ln -s  catchcn_internal_rst catchcn_internal_rst.$YYYYMMDD
endif

echo DONE > done_rst_file
_EOI2_

    cat mkLDASsa.j2 >>mkLDASsa.j
    rm mkLDASsa.j2
    sbatch mkLDASsa.j
    cd $PWD
    breaksw

case [1]:

    set coordfile=${RESTART_short}/rc_out/${RESTART_ID}.ldas_tilecoord.bin
    if (-e $coordfile ) then
        set ENDI=`file -b $coordfile | rev | cut -d' ' -f2 | rev`
    else
        echo "WARNING :  no file $coordfile. DEFAULT GEOSldas restart"
        set ENDI = LSB
    endif


    ## restart is from old LDAS which produce big endian binary
    if($ENDI == MSB) then
        echo ' '
        mkdir -p $EXPDIR/$EXPID/mk_restarts/OutData2/
        ln -s $BCSDIR/$TILFILE $EXPDIR/$EXPID/mk_restarts/OutData2/OutTileFile
        ln -s $BCSDIR/clsm $EXPDIR/$EXPID/mk_restarts/OutData2/clsm
        ln -s $INSTDIR/bin $EXPDIR/$EXPID/mk_restarts/

        cd $EXPDIR/$EXPID/mk_restarts/
        echo '#\!/bin/csh -f ' > this.file
        echo 'source $INSTDIR/bin/g5_modules' >> this.file
        echo 'if ( -e /etc/os-release ) then' >> this.file
        echo '    module load nco/4.8.1'      >> this.file
        echo 'else'                           >> this.file
        echo '    module load other/nco-4.6.8-gcc-5.3-sp3 ' >> this.file
        echo 'endif' >> this.file
        #echo 'setenv OMPI_MCA_shmem_mmap_enable_nfs_warning 0' >> this.file
        echo 'setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BASEDIR}/Linux/lib' >> this.file

        set mpi_mpmd = "${INSTDIR}/bin/esma_mpirun -np 1 bin/mk_GEOSldasRestarts.x -b ${BCSDIR} -d ${YYYYMMDD} -e ${RESTART_ID} -k 0000 -l ${RESTART_short} -m ${MODEL} -s ${SURFLAY} -r Y -t ${TILFILE}"
        set j = 1
        while ($j < $NUMENS)
           set ENS = `printf '%04d' $j`
           set mpi_mpmd = "${mpi_mpmd} : -np 1 bin/mk_GEOSldasRestarts.x -b ${BCSDIR} -d ${YYYYMMDD} -e ${RESTART_ID} -k ${ENS} -l ${RESTART_short} -m ${MODEL} -s ${SURFLAY} -r Y -t ${TILFILE}"
           @ j++
        end
        echo $mpi_mpmd >> this.file

        set j = 0
        while ($j < $NUMENS)
           set ENS = `printf '%04d' $j`
           echo  "ncks -4  -O -h -x -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF ${MODEL}${ENS}_internal_rst.${YYYYMMDD} ${MODEL}${ENS}_internal_rst.${YYYYMMDD} &"  >> this.file
           @ j++
        end
        echo 'wait' >> this.file

        chmod +x this.file
        ./this.file
        rm -f this.file
        echo DONE > done_rst_file
        cd $PWD
        sleep 1
     else
    ## restart is from GEOSldas little endian binary
       if( ! -e $EXPDIR/$EXPID/mk_restarts ) mkdir -p $EXPDIR/$EXPID/mk_restarts/
       echo DONE > $EXPDIR/$EXPID/mk_restarts/done_rst_file
       sleep 1
     endif
    
     breaksw

case [2]:

    echo ' '
    mkdir -p $EXPDIR/$EXPID/mk_restarts/OutData1/
    mkdir -p $EXPDIR/$EXPID/mk_restarts/OutData2/
    ln -s $BCSDIR/$TILFILE $EXPDIR/$EXPID/mk_restarts/OutData1/OutTileFile
    ln -s $BCSDIR/$TILFILE $EXPDIR/$EXPID/mk_restarts/OutData2/OutTileFile
    ln -s $BCSDIR/clsm $EXPDIR/$EXPID/mk_restarts/OutData2/clsm
    ln -s $INSTDIR/bin $EXPDIR/$EXPID/mk_restarts/
    
    cd $EXPDIR/$EXPID/mk_restarts/

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
if ( -e /etc/os-release ) then
  module load nco/4.8.1
else
  module load other/nco-4.6.8-gcc-5.3-sp3 
endif
setenv LAIFILE `find ${BCSDIR}/lai_clim*`
limit stacksize unlimited
 
$INSTDIR/bin/esma_mpirun -np 56 bin/mk_GEOSldasRestarts.x -b ${BCSDIR} -d ${YYYYMMDD} -e ${RESTART_ID} -l ${RESTART_short} -t ${TILFILE} -m ${MODEL} -s $SURFLAY -j Y -r R -p ${PARAM_FILE}
sleep 3

_EOI3_

    if($LSM_CHOICE == 1) sed -i '$ a\bin/Scale_Catch OutData1/catch_internal_rst OutData2/catch_internal_rst catch_internal_rst $SURFLAY $WEMIN_IN $WEMIN_OUT \'  mkLDASsa.j
    if($LSM_CHOICE == 2) sed -i '$ a\bin/Scale_CatchCN OutData1/catchcn_internal_rst OutData2/catchcn_internal_rst catchcn_internal_rst $SURFLAY $WEMIN_IN $WEMIN_OUT \'  mkLDASsa.j

    sed -i '$ a\ \'  mkLDASsa.j
    sed -i '$ a\## Done creating catch*_internal_rst file \'  mkLDASsa.j
    sed -i '$ a\ \'  mkLDASsa.j

    cat << _EOI4_ > mkLDASsa.j2
sleep 2

if($LSM_CHOICE == 1) then
   if (-f irrigation_internal_rst && $RUN_IRRIG == 1) then 
      ncks -4  -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF irrigation_internal_rst -A catch_internal_rst
   endif
   ln -s  catch_internal_rst catch_internal_rst.$YYYYMMDD
else
   if (-f irrigation_internal_rst && $RUN_IRRIG == 1) then 
      ncks -4 -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF irrigation_internal_rst -A catchcn_internal_rst 
   endif
   ln -s  catchcn_internal_rst catchcn_internal_rst.$YYYYMMDD
endif

echo DONE > done_rst_file
_EOI4_

    cat mkLDASsa.j2 >>mkLDASsa.j
    rm mkLDASsa.j2
    sbatch mkLDASsa.j
    cd $PWD
    breaksw

case [FGM]:

    set YYYY = `echo $YYYYMMDD | cut -c1-4`
    set   MM = `echo $YYYYMMDD | cut -c5-6`
    set   DD = `echo $YYYYMMDD | cut -c7-8`

    mkdir -p $EXPDIR/$EXPID/mk_restarts/InData/
    mkdir -p $EXPDIR/$EXPID/mk_restarts/OutData.1/
    mkdir -p $EXPDIR/$EXPID/mk_restarts/OutData.2/

    if ($HAVE_RESTART == M) then
	set ARCDIR = /archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/d5124_m2_jan79/rs/Y ; set mlable = jan79
	if ($YYYY > 1991) set ARCDIR = /archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/d5124_m2_jan91/rs/Y ; set mlable = jan91
	if ($YYYY > 2000) set ARCDIR = /archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/d5124_m2_jan00/rs/Y ; set mlable = jan00
	if ($YYYY > 2010) set ARCDIR = /archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/d5124_m2_jan10/rs/Y ; set mlable = jan10	
	set rstfile = ${ARCDIR}${YYYY}/M${MM}/d5124_m2_${mlable}.catch_internal_rst.${YYYYMMDD}_21z.bin
	dmget $rstfile
	set INTILFILE = /gpfsm/dnb02/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_MERRA-2/CF0180x6C_DE1440xPE0720/CF0180x6C_DE1440xPE0720-Pfafstetter.til
	set  WEMIN_IN = 26
	/bin/cp -p $rstfile $EXPDIR/$EXPID/mk_restarts/InData/M2Restart
    endif

    if ($HAVE_RESTART == G) then
	set rstfile = `echo $RESTART_PATH | rev | cut -c 2- | rev`
	set INTILFILE = `readlink $RESTART_ID/scratch/tile.data`
	if ( `echo $INTILFILE | grep -n 'NLv3'` == '' ) set  WEMIN_IN = 26
	/bin/cp -p $rstfile $EXPDIR/$EXPID/mk_restarts/InData/M2Restart
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
	    set   rstfile = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.catch_internal_rst.${YYYYMMDD}_21z.bin
	    dmget $rstfile
	    /bin/cp -p $rstfile $EXPDIR/$EXPID/mk_restarts/InData/M2Restart
	endif
	if (($expdate >= $date_17) && ($expdate < $date_21)) then
	    set fpver = GEOS-5.17/GEOSadas-5_17/
	    set fplab = f517_fp
	    set INTILFILE = /discover/nobackup/ltakacs/bcs/Icarus/Icarus_Ostia/CF0720x6C_CF0720x6C/CF0720x6C_CF0720x6C-Pfafstetter.til
  	    set  WEMIN_IN = 26
	    set TARFILE = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.rst.${YYYYMMDD}_21z.tar
	    dmget $TARFILE
	    set   rstfile = ${fplab}.catch_internal_rst.${YYYYMMDD}_21z.nc4
	    tar -xvf $TARFILE $rstfile && /bin/mv $rstfile $EXPDIR/$EXPID/mk_restarts/InData/M2Restart	 	    
	endif
	if (($expdate >= $date_21) && ($expdate < $date_22)) then
	    set fpver = GEOS-5.21/GEOSadas-5_21/
	    set fplab = f521_fp
	    set INTILFILE = /discover/nobackup/ltakacs/bcs/Icarus/Icarus_Ostia/CF0720x6C_CF0720x6C/CF0720x6C_CF0720x6C-Pfafstetter.til
  	    set  WEMIN_IN = 26
	    set TARFILE = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.rst.${YYYYMMDD}_21z.tar
	    dmget $TARFILE
	    set   rstfile = ${fplab}.catch_internal_rst.${YYYYMMDD}_21z.nc4
	    tar -xvf $TARFILE $rstfile && /bin/mv $rstfile $EXPDIR/$EXPID/mk_restarts/InData/M2Restart	 	    
	endif
	if (($expdate >= $date_22) && ($expdate < $date_25)) then
	    set fpver = GEOS-5.22/GEOSadas-5_22/
	    set fplab = f522_fpp
	    set INTILFILE = /discover/nobackup/ltakacs/bcs/Icarus/Icarus_Ostia/CF0720x6C_CF0720x6C/CF0720x6C_CF0720x6C-Pfafstetter.til
  	    set  WEMIN_IN = 26
	    set TARFILE = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.rst.${YYYYMMDD}_21z.tar
	    dmget $TARFILE
	    set   rstfile = ${fplab}.catch_internal_rst.${YYYYMMDD}_21z.nc4
	    tar -xvf $TARFILE $rstfile && /bin/mv $rstfile $EXPDIR/$EXPID/mk_restarts/InData/M2Restart	 	    
	endif
	if ($expdate >= $date_25) then
	    set fpver = GEOS-5.25/GEOSadas-5_25/
	    set fplab = f525land_fpp
	    set INTILFILE = /discover/nobackup/ltakacs/bcs/Icarus-NLv3/Icarus-NLv3_Ostia/CF0720x6C_CF0720x6C/CF0720x6C_CF0720x6C-Pfafstetter.til
  	    set  WEMIN_IN = 13
	    set TARFILE = /archive/u/dao_ops/$fpver/${fplab}/rs/Y${YYYY}/M${MM}/${fplab}.rst.${YYYYMMDD}_21z.tar
	    dmget $TARFILE
	    set   rstfile = ${fplab}.catch_internal_rst.${YYYYMMDD}_21z.nc4
	    tar -xvf $TARFILE $rstfile && /bin/mv $rstfile $EXPDIR/$EXPID/mk_restarts/InData/M2Restart
	endif

    endif
    
    /bin/ln -s $INTILFILE $EXPDIR/$EXPID/mk_restarts/InData/InTilFile
    /bin/ln -s $BCSDIR/$TILFILE $EXPDIR/$EXPID/mk_restarts/OutData.1/OutTilFile
    /bin/ln -s $BCSDIR/$TILFILE $EXPDIR/$EXPID/mk_restarts/OutData.2/OutTilFile
    /bin/ln -s $BCSDIR/clsm $EXPDIR/$EXPID/mk_restarts/OutData.2/clsm
    /bin/ln -s $INSTDIR/bin $EXPDIR/$EXPID/mk_restarts/
    cd $EXPDIR/$EXPID/mk_restarts/

    cat << _EOI5_ > mkLDASsa.j
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
#setenv MKL_CBWR SSE4_2 # ensure zero-diff across archs
#setenv MV2_ON_DEMAND_THRESHOLD 8192 # MVAPICH2
setenv LAIFILE `find ${BCSDIR}/lai_clim*`
if ( -e /etc/os-release ) then
  module load nco/4.8.1
else
  module load other/nco-4.6.8-gcc-5.3-sp3 
endif
limit stacksize unlimited
 
/bin/ln -s OutData.1 OutData
if($LSM_CHOICE == 1) then 
$INSTDIR/bin/esma_mpirun -np 56 bin/mk_CatchRestarts OutData/OutTilFile InData/InTilFile InData/M2Restart $SURFLAY 4 
else
$INSTDIR/bin/esma_mpirun -np 56 bin/mk_CatchCNRestarts OutData/OutTilFile InData/InTilFile InData/M2Restart $SURFLAY $YYYYMMDD 4 
endif
/bin/rm OutData

/bin/ln -s OutData.2 OutData
if($LSM_CHOICE == 1) then 
$INSTDIR/bin/esma_mpirun -np 1 bin/mk_CatchRestarts OutData/OutTilFile OutData.1/OutTilFile OutData.1/M2Restart $SURFLAY 4 
else
$INSTDIR/bin/esma_mpirun -np 1 bin/mk_CatchCNRestarts OutData/OutTilFile OutData.1/OutTilFile OutData.1/M2Restart $SURFLAY $YYYYMMDD 4 
endif
/bin/rm OutData

if($LSM_CHOICE == 1) then 
bin/Scale_Catch OutData.1/M2Restart OutData.2/M2Restart catch_internal_rst $SURFLAY 26 $WEMIN_OUT
else
bin/Scale_CatchCN OutData.1/M2Restart OutData.2/M2Restart catchcn_internal_rst $SURFLAY 26 $WEMIN_OUT
endif

if (-f irrigation_internal_rst && $RUN_IRRIG == 1) then 
    ncks -4  -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF irrigation_internal_rst -A catch_internal_rst
endif
/bin/ln -s  catch_internal_rst catch_internal_rst.$YYYYMMDD

echo DONE > done_rst_file  
_EOI5_

    sbatch mkLDASsa.j
    cd $PWD
    breaksw

default :
    echo $HAVE_RESTART is not implemented
endsw
