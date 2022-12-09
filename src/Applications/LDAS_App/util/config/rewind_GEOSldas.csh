#!/bin/csh

# rewind existing GEOSldas run to specified date/time 

setenv MYNAME rewind_GEOSldas.csh

if ( $#argv < 4 ) then
   echo " "
   echo " NAME "
   echo " "
   echo "   $MYNAME - rewind existing GEOSldas run to restart time of nymd nhms"
   echo " "
   echo " SYNOPSIS "
   echo " "
   echo "   $MYNAME nymd nhms expid exppath "
   echo " "
   echo " where "
   echo "   nymd    - restart date, as YYYYMMDD "
   echo "   time    - restart time, as HHMMSS "
   echo "   expid   - experiment name, e.g., ldas4coup "
   echo "   exppath - run directory path, e.g., /discover/nobackup/[user]/ "
   echo " "
   echo " DESCRIPTION "
   echo " "
   echo "   This procedure rewinds and resets the GEOSldas experiment "
   echo "   specified by expid and exppath to the restart date/time " 
   echo "   specified by nymd and nhms. "
   echo " "
   echo " Example of valid command line: "
   echo "   $MYNAME 20170829 210000 ldas4coup /discover/nobackup/qzhang "
   exit(0)
endif

set nymd  = $1
set nhms  = $2 

echo " ymd = $nymd  " 
echo " hms = $nhms  "

set yin = `echo $nymd | cut -c1-4`
set min = `echo $nymd | cut -c5-6`
set hin = `echo $nhms | cut -c1-4`

set date = ${nymd}_${hin} 

set expid   = $3
set rundir  = $4
echo " expid  = $expid  "

cd ${rundir}/${expid}/run
set nmem  = `grep NUM_LDAS_ENSEMBLE:  LDAS.rc | cut -d':' -f2`
cd ${rundir}/${expid} 
set grid  = `ls output`

## rewind links to restart files  

@ NENS = $nmem 

set rsout = ${rundir}/${expid}/output/${grid}/rs 
set rstin = ${rundir}/${expid}/input/restart 
cd $rstin 

/bin/rm -rf catch*_internal_rst
/bin/rm -rf landpert*_internal_rst
/bin/rm -rf landassim_obspertrseed*_rst

@ inens = 1

while ($inens <= $NENS)

    if     ($inens <   10) then
        set ENSDIR = `echo          ens000${inens}`
        set catin  = `echo        catch000${inens}`
        set pertin = `echo     landpert000${inens}`
        set seedin = `echo obspertrseed000${inens}`
    else if($inens <  100) then
        set ENSDIR = `echo           ens00${inens}`
        set catin  = `echo         catch00${inens}`
        set pertin = `echo      landpert00${inens}`
        set seedin = `echo  obspertrseed00${inens}`
    else if($inens < 1000) then
        set ENSDIR = `echo            ens0${inens}`
        set catin  = `echo          catch0${inens}`
        set pertin = `echo       landpert0${inens}`
        set seedin = `echo   obspertrseed0${inens}`
    else
        set ENSDIR = `echo             ens${inens}`
        set catin  = `echo           catch${inens}`
        set pertin = `echo        landpert${inens}`
        set seedin = `echo    obspertrseed${inens}`
    endif

    /bin/ln -s ${rsout}/${ENSDIR}/Y${yin}/M${min}/${expid}.catch_internal_rst.${date} ${catin}_internal_rst

    if (-e ${rsout}/${ENSDIR}/Y${yin}/M${min}/${expid}.landpert_internal_rst.${date}.gz ) then  
      gunzip ${rsout}/${ENSDIR}/Y${yin}/M${min}/${expid}.landpert_internal_rst.${date}.gz 
    endif 

    /bin/ln -s ${rsout}/${ENSDIR}/Y${yin}/M${min}/${expid}.landpert_internal_rst.${date} ${pertin}_internal_rst

    /bin/ln -s ${rsout}/${ENSDIR}/Y${yin}/M${min}/${expid}.landassim_obspertrseed_rst.${date} landassim_${seedin}_rst

    @ inens ++
end

## -- remove records in rc_out 
cd  ${rundir}/${expid}/output/${grid}/rc_out/Y${yin}/M${min} 
/bin/rm -rf *smapL4*${date}z.bin 
/bin/rm -rf *.${date}z.txt
/bin/rm -rf *.${date}z.nml

## -- reset cap_restart 
cd  ${rundir}/${expid}/run 
/bin/rm -rf cap_restart
echo $nymd ${nhms} > cap_restart 

## EOF ####################################################################
