#!/bin/bash
#SBATCH --job-name=start_jobs
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=./start_jobs_out.txt
#SBATCH --partition=catch-m6a4xl-demand


# things to edit
total_particles=30
particles_per_batch=30
sleep_hours=0
sleep_min=0
####

# initialize the counter
n=0
# initialize the batch counter
this_batch_max=0
# how many we want to start in the first batch
this_batch_max=$(( $this_batch_max + $particles_per_batch ))


# while we haven't yet started all particles
while [ "$n" -lt "$total_particles" ];
do
    # while we are still starting particles in this batch
    while [ "$n" -lt "$this_batch_max" ];
    do
        # start the particles
        cd ./$n/run
        sbatch lenkf.j
        echo "started $n"
        cd ../../
        # update the counter
        n=$(( $n + 1 ))
        sleep .2
    done
    # update what particle we stop at for our next batch
    this_batch_max=$(( $this_batch_max + $particles_per_batch))
    # if we have started all of our particles exit so that we don't
    # wait for nothing
    if [ "$this_batch_max" -gt "$total_particles" ]; then
        exit 0
    fi
    # otherwise sleep the specified time and restart
    sleep "$sleep_hours"h "$sleep_min"m 0s
done
