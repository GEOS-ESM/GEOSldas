# run from /path/to/GEOSldas/ to install the just-built directory
exeinp=GEOSldas_CN45_exeinp_default.txt
batinp=GEOSldas_CN45_batinp_default.txt
exp_dir_gen=/discover/nobackup/trobinet/exps/
exp_to_rm=/discover/nobackup/trobinet/exps/GEOSldas_CN45_quick
exp_name=GEOSldas_CN45_quick

# use variables to set up
yes n | ./parallel_build.csh
cd ./install/bin
rm -rf $exp_to_rm
source g5_modules.sh
./ldas_setup setup $exp_dir_gen ./$exeinp ./$batinp
cp ../../src/Applications/LDAS_App/HISTORY.rc /discover/nobackup/trobinet/exps/$exp_name/run
cd /discover/nobackup/trobinet/exps/$exp_name/run
#cp /discover/nobackup/trobinet/GEOSldas_EFPH_develop/GEOSldas/src/Applications/LDAS_App/HISTORY.rc .
sed -i "/EXPID:/c EXPID:  $exp_name" HISTORY.rc
#./lenkf.j
sbatch lenkf.j
cd ../scratch
touch GEOSldas_log_txt
tail -f GEOSldas_log_txt
