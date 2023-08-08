#!/bin/bash
# run from /path/to/GEOSldas/ to install the just-built directory

# to be set by the user
model_dir_gen="/discover/nobackup/trobinet/GEOSldas_pso_g1_ksat/GEOSldas"
exp_dir_gen="/discover/nobackup/trobinet/exps"
exp_name="GEOSldas_CN45_med_default"
include_name="/discover/nobackup/trobinet/pso/step_1_choose_tiles/outputs/include"
num_members=1
# directly based off of user input (shouldn't need to edit)
exeinp="${model_dir_gen}/GEOSldas_CN45_exeinp_default.txt"
batinp="${model_dir_gen}/GEOSldas_CN45_batinp_default.txt"
exp_to_rm="${exp_dir_gen}/${exp_name}"
exp_dir="${exp_dir_gen}/${exp_name}"
hist_rc_loc="${model_dir_gen}/src/Applications/LDAS_App/HISTORY.rc"
exp_to_rm="${exp_dir_gen}/${exp_name}"
include_to_place="${model_dir_gen}/src/Applications/LDAS_App/include"

# use variables to set up
num_mems_idx=$(($num_members-1))
yes n | ./parallel_build.csh
read -p "Are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
  exit 1
fi
rm -rf $exp_to_rm
mkdir $exp_dir
cp $include_name $include_to_place
source "${model_dir_gen}/install/bin/g5_modules.sh"
for (( c=0; c<=$num_mems_idx; c++ ))
do
  echo $c
  sed -i "17s/.*/EXP_ID: $c/" $exeinp
  #sed -i "13s/.*/EXP_ID: $c/" ./src/Applications/LDAS_App/HISTORY.rc
  #cd "${model_dir_gen}/install/bin"
  cd "${model_dir_gen}/install/bin"
  echo "ldas setup command:"
  echo ""${model_dir_gen}/install/bin/ldas_setup" setup $exp_dir $exeinp $batinp"
  "${model_dir_gen}/install/bin/ldas_setup" setup $exp_dir $exeinp $batinp
  cd ../../
  #cp ../../src/Applications/LDAS_App/HISTORY.rc $exp_dir_gen/$exp_name/run
  #cd $exp_dir/$c/run
  #sed -i "/EXPID:/c EXPID:  $exp_name" HISTORY.rc
  #./lenkf.j
  echo "$c" > "${exp_dir}/${c}/run/ens_num.txt"
  #cp $hist_rc_loc .
  sbatch "${exp_dir}/${c}/run/lenkf.j"
  #cd "${exp_dir}/${c}/scratch"
  touch "${exp_dir}/${c}/scratch/GEOSldas_log_txt"
  #cd $cur_dir
done
tail -f "${exp_dir}/0/scratch/GEOSldas_log_txt"
