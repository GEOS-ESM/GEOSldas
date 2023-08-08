#!/bin/bash
# user must update num_parts
num_parts=10
# restart value is passed to script
# if restart is 0, ./lenf.j will be run
# if restart is 1, sbatch lenkf.j will be run
restart=$1
num_parts_idx=$(($num_parts-1))
echo $num_parts
echo $num_parts_idx
if [ $restart -eq 1 ]
then
  for (( c=0; c<=$num_parts_idx; c++ ))
  do
    sbatch lenkf.j
  done
