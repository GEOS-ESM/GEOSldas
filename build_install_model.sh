# run from /path/to/GEOSldas/ to install the just-built directory
exeinp=GEOSldas_CN45_exeinp_EFPH.txt
batinp=GEOSldas_CN45_batinp_EFPH.txt
exp_dir_gen=/discover/nobackup/trobinet/exps/
exp_to_rm=/discover/nobackup/trobinet/exps/GEOSldas_CN45_EFPH
exp_name=GEOSldas_CN45_EFPH

# use variables to set up
./parallel_build_install.csh
# I ELIMINATED CLEAN OPTION HERE
# TO GO BACK AND ALLOW FOR CLEANING AGAIN, YOU MUST:
# GO TO ./@env/build.csh
# Go to line 500-ish
# Comment out "set do_clean = n", uncomment the similar line just above it
cd ./install/bin
rm -rf $exp_to_rm
source g5_modules.sh
./ldas_setup setup $exp_dir_gen ./$exeinp ./$batinp
cd /discover/nobackup/trobinet/exps/$exp_name/run
cp /discover/nobackup/trobinet/GEOSldas_EFPH_develop/GEOSldas/src/Applications/LDAS_App/HISTORY.rc .
sbatch lenkf.j
cd ../
cd ./scratch
touch GEOSldas_log_txt
tail -f GEOSldas_log_txt
