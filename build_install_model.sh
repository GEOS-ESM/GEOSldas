#!/bin/bash
# run from /path/to/GEOSldas/ to install the just-built directory

# to be set by the user
# what is the directory where the model information is stored
model_dir_gen="/shared/models/GEOSldas_pso_g1_ksat/GEOSldas"
# what is the directory where we should place the experiment
exp_dir_gen="/lustre/catchment/exps"
# what is the name of this experiment
exp_name="GEOSldas_CN45_pso_g1_et_v2"
# where is the list of tiles to run?
include_name="/shared/pso/step_1_choose_tiles/outputs/include"
# how many particles are we running for the pso?
# if want just a normal run, make this 1
num_members=10

# directly based off of user input (shouldn't need to edit)
# what is the exeinp file
exeinp="${model_dir_gen}/GEOSldas_CN45_exeinp_default.txt"
# what is the batinp file
batinp="${model_dir_gen}/GEOSldas_CN45_batinp_default.txt"
# what is the full directory where we will palce the experiment
exp_dir="${exp_dir_gen}/${exp_name}"
# where is the history.rc file?
hist_rc_loc="${model_dir_gen}/src/Applications/LDAS_App/HISTORY.rc"
# what experiment do we need to remove to set this up
exp_to_rm="${exp_dir_gen}/${exp_name}"
# we will copy the include file to this model. where are we copying to?
include_to_place="${model_dir_gen}/src/Applications/LDAS_App/include"

# need to change ulimit first for compute nodes
ulimit -s unlimited
# need to get the idx for number of members to run, since we will start
# counting at zeros
num_mems_idx=$(($num_members-1))
# ask the user if they want to rebuild
read -rep $'Re-build model? <y/n>\n' rebuild
if [ $rebuild = y ]; then
    # if they want to rebuild, call manual_build
    ./manual_build.sh
    # make sure they want to continue; will allow them to stop if there is an
    # error in manual_build
    read -rep $'Are you sure you want to continue? <y/n>\n' sure
    if [ $sure = n ]; then
      exit 1
    fi
fi
# remove the experiment that is needed to build new directories
rm -rf $exp_to_rm
# create the new experiment dir
mkdir $exp_dir
# copy the include file
cp $include_name $include_to_place
# loop over every pso member to set up that directory
for (( c=0; c<=$num_mems_idx; c++ ))
do
  echo $c
  # change the exp id to be this particle number
  sed -i "17s/.*/EXP_ID: $c/" $exeinp
  # do this also for HISTORY.rc
  sed -i "13s/.*/EXPID:  $c/" "${model_dir_gen}/src/Applications/LDAS_App/HISTORY.rc"
  # go to install/bin
  cd "${model_dir_gen}/install/bin"
  # run ldas setup
  echo "ldas setup command:"
  echo ""${model_dir_gen}/install/bin/ldas_setup" setup $exp_dir $exeinp $batinp"
  "${model_dir_gen}/install/bin/ldas_setup" setup $exp_dir $exeinp $batinp
  # go back to top directory
  cd ../../
  echo "$c" > "${exp_dir}/${c}/run/ens_num.txt"
  # add the correct partition to lenkf.j
  sed -i '18 i #SBATCH --partition=catch-m6a4xl-spot' "${exp_dir}/${c}/run/lenkf.j"
  touch "${exp_dir}/${c}/scratch/GEOSldas_log_txt"
done
#tail -f "${exp_dir}/0/scratch/GEOSldas_log_txt"
# NOTE: have to sbatch everything manually on AWS; doesn't work from compute
# node
