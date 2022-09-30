# run from /path/to/GEOSldas/ to install the just-built directory
exeinp=GEOSldas_CN45_exeinp_PSO.txt
batinp=GEOSldas_CN45_batinp_PSO.txt
exp_dir_gen=/discover/nobackup/trobinet/exps/
exp_to_rm=/discover/nobackup/trobinet/exps/GEOSldas_CN45_PSO_v3
exp_name=GEOSldas_CN45_PSO_v3

# use variables to set up
./parallel_build.csh
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
sbatch lenkf.j
