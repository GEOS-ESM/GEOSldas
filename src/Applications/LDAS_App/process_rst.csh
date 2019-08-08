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
setenv  RESTART_short ${RESTART_PATH}/${RESTART_ID}/output/${RESTART_DOMAIN}/

set PWD=`pwd`

set THISDIR = `echo $PWD | rev | cut -d'/' -f1 | rev`
if ( $THISDIR == LDAS_App ) setenv ESMADIR `echo $PWD | rev | cut -d'/' -f4- | rev`
if ( $THISDIR == bin ) setenv ESMADIR `echo $PWD | rev | cut -d'/' -f3- | rev`
#setenv ESMADIR2 /gpfsm/dnb02/smahanam/LDAS/GEOSldas/IU_V24C05_GEOSldas_GOSWIMplus/GEOSagcm/


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
    ln -s $ESMADIR/install/bin $EXPDIR/$EXPID/mk_restarts/

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
 
source $ESMADIR/install/bin/g5_modules
setenv OMPI_MCA_shmem_mmap_enable_nfs_warning 0
#setenv MKL_CBWR SSE4_2 # ensure zero-diff across archs
#setenv MV2_ON_DEMAND_THRESHOLD 8192 # MVAPICH2
setenv LAIFILE `find ${BCSDIR}/lai_clim*`
setenv PATH $PATH\:/usr/local/other/SLES11.3/nco/4.6.8/gcc-5.3-sp3/bin/
limit stacksize unlimited
 
#mpirun -map-by core --mca btl ^vader -np 56 bin/mk_LDASsaRestarts -a ${SPONSORID} -b ${BCSDIR} -t ${TILFILE} -m ${MODEL} -s 50 -j Y
$ESMADIR/install/bin/esma_mpirun -np 56 bin/mk_LDASsaRestarts -a ${SPONSORID} -b ${BCSDIR} -t ${TILFILE} -m ${MODEL} -s 50 -j Y
sleep 3

if($LSM_CHOICE == 1) then
   /bin/cp OutData1/catch_internal_rst OutData2/catch_internal_rst
else
   /bin/cp OutData1/catchcn_internal_rst OutData2/catchcn_internal_rst
endif

#mpirun -map-by core --mca btl ^vader -np 56 bin/mk_LDASsaRestarts -a ${SPONSORID} -b ${BCSDIR} -t ${TILFILE} -m ${MODEL} -s 50 -j Y
$ESMADIR/install/bin/esma_mpirun -np 56 bin/mk_LDASsaRestarts -a ${SPONSORID} -b ${BCSDIR} -t ${TILFILE} -m ${MODEL} -s 50 -j Y

_EOI_

    if($LSM_CHOICE == 1) sed -i '$ a\bin/Scale_Catch OutData1/catch_internal_rst OutData2/catch_internal_rst catch_internal_rst 50 \'  mkLDASsa.j
    if($LSM_CHOICE == 2) sed -i '$ a\bin/Scale_CatchCN OutData1/catchcn_internal_rst OutData2/catchcn_internal_rst catchcn_internal_rst 50 \'  mkLDASsa.j
    #if($LSM_CHOICE == 2) sed -i '$ a\bin/Scale_Catch OutData1/catchcn_internal_clm45 OutData2/catchcn_internal_clm45 catchcn_internal_clm45 50 \'  mkLDASsa.j

    sed -i '$ a\ \'  mkLDASsa.j
    sed -i '$ a\## Done creating catch*_internal_rst file \'  mkLDASsa.j
    sed -i '$ a\ \'  mkLDASsa.j

    cat << _EOI2_ > mkLDASsa.j2
sleep 2

if($LSM_CHOICE == 1) then
   if (-f irrigation_internal_rst && $RUN_IRRIG == 1) then 
      ncks -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF irrigation_internal_rst -A catch_internal_rst
   endif
   ln -s  catch_internal_rst catch_internal_rst.$YYYYMMDD
else
   if (-f irrigation_internal_rst && $RUN_IRRIG == 1) then 
      ncks -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF irrigation_internal_rst -A catchcn_internal_rst 
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
        echo 'restart is from old LDAS'
        mkdir -p $EXPDIR/$EXPID/mk_restarts/OutData2/
        ln -s $BCSDIR/$TILFILE $EXPDIR/$EXPID/mk_restarts/OutData2/OutTileFile
        ln -s $BCSDIR/clsm $EXPDIR/$EXPID/mk_restarts/OutData2/clsm
        ln -s $ESMADIR/install/bin $EXPDIR/$EXPID/mk_restarts/

        cd $EXPDIR/$EXPID/mk_restarts/
        rm -f this.file
        echo '#\!/bin/csh -f ' > this.file
        echo '#SBATCH --account='${SPONSORID} >> this.file
        echo '#SBATCH --time=1:00:00' >> this.file
        echo '#SBATCH --ntasks=1' >> this.file
        echo '#SBATCH --job-name=LDASrst2GEOSrst' >> this.file
        echo '#SBATCH --constraint=hasw' >> this.file
        echo '#SBATCH --qos=debug' >> this.file
        echo '#SBATCH --output=LDASrst2GEOSrst.o' >> this.file
        echo '#SBATCH --error=LDASrst2GEOSrst.e' >> this.file
        echo 'source $ESMADIR/install/bin/g5_modules' >> this.file
        echo 'setenv OMPI_MCA_shmem_mmap_enable_nfs_warning 0' >> this.file
        echo 'setenv PATH $PATH\:/usr/local/other/SLES11.3/nco/4.6.8/gcc-5.3-sp3/bin/' >> this.file

        set j = 0
        while ($j < $NUMENS)
           set ENS = `printf '%04d' $j`
           echo $ESMADIR/install/bin/esma_mpirun -np 1 bin/mk_LDASsaRestarts -b ${BCSDIR} -d ${YYYYMMDD} -e ${RESTART_ID} -k ${ENS} -l ${RESTART_short} -m ${MODEL} -s 50 -r Y -t ${TILFILE}  >> this.file
           echo  ncks -O -h -x -v IRRIGFRAC,PADDYFRAC,LAIMIN,LAIMAX,CLMPT,CLMST,CLMPF,CLMSF ${MODEL}${ENS}_internal_rst.${YYYYMMDD} ${MODEL}${ENS}_internal_rst.${YYYYMMDD} >> this.file 
           @ j++
        end
        
        echo 'echo DONE > done_rst_file' >> this.file
        chmod +x this.file
        sbatch this.file
       # ./this.file
       # echo DONE > done_rst_file
        cd $PWD
        sleep 1
     else
    ## restart is from GEOSldas little endian binary
       if( ! -e $EXPDIR/$EXPID/mk_restarts ) mkdir -p $EXPDIR/$EXPID/mk_restarts/
       echo DONE > $EXPDIR/$EXPID/mk_restarts/done_rst_file
       sleep 1
     endif
    
     breaksw
default :
    echo $HAVE_RESTART is not implemented
endsw
