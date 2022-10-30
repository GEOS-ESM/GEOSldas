# script that will run the PSO algorithm

#activate conda environment

# import libraries
from os.path import exists
import pandas as pd
import numpy as np
import os
import glob
from netCDF4 import Dataset
import sys
from compute_objective_test import compute_objective



# set important parameters
convergence_threshold = 15 #num iterations with no change to global_best to converge
nparticles = 20 #fits architecture of NCCS?
nparameters = 5


# do PSO calculations
#PSO parameters
w = 1 #SEE WHAT WAS USED PREVIOUSLY; inertia weight that slows down the velocity of particle at previous step
c1 = 1 #SEE WHAT WAS USED PREVIOUSLY; weight put on the individual particle's best position (higher=more exploration)
c2 = 1 #SEE WHAT WAS USED PREVIOUSLY; weight put on the global best position (higher=faster convergence)\
r1 = np.random.random() #should this random number be between 0 and 1 or should I scale?
r2 = np.random.random() #should this random number be between 0 and 1 or should I scale?

num_iterations = 0
num_without_restart = 0
iterations_without_change = 0
positions = np.random.rand(nparticles,nparameters)*5
velocities = np.random.rand(nparticles,nparameters)*5
global_best_positions_last = positions
global_best_objective_last = float('inf')
best_position_objective_last = np.zeros((nparticles))
best_position_objective_last[:] = float('inf')
best_position_last = positions


while iterations_without_change < convergence_threshold:
    objective_result = compute_objective(nparticles,nparameters,positions) #objective function to be filled in
    #best global position
    #update global best position
    objective_result_min_idx = np.argmin(objective_result)
    objective_result_min = objective_result[objective_result_min_idx]
    if objective_result_min < global_best_objective_last:
        global_best_positions = positions[objective_result_min_idx,:]
        global_best_objective = objective_result_min

    #best local position
    #update local best position
    #local_objective_result = np.random.random(nparticles) #objective_function(best_position)
    local_improved_idx = np.where(objective_result < best_position_objective_last)
    best_position = np.copy(best_position_last)
    best_position_objective = np.copy(best_position_objective_last)
    best_position[local_improved_idx] = positions[local_improved_idx]
    best_position_objective[local_improved_idx] = objective_result[local_improved_idx]

    #update velocity
    velocity_next = (w*velocities + c1*r2*(best_position - positions) +
                     c2*r2*(np.tile(global_best_positions,(nparticles,1)) -
                            positions)
                    )

    #update positions
    positions_next = positions + velocity_next

    #save the updated arrays to their respective .csv files to be called in the next iteration
    if np.all(global_best_positions_last == global_best_positions):
        iterations_without_change += 1
    else:
        iterations_without_change = 0

    print('iterations without change')
    print(iterations_without_change)
    print('iterations')
    print(num_iterations)

    #check if convergence has occured
    if iterations_without_change == convergence_threshold:
        is_converged = 1
    else:
        is_converged = 0

    #write to PSO_tracker.csv to keep track of variables
    num_iterations += 1
    num_without_restart += 1
    # more will go here as I figure out what PSO values we will need in the next iteration....

    global_best_positions_last = global_best_positions
    global_best_objective_last = global_best_objective
    best_position_objective_last = best_position_objective
    best_position_last = best_position

    positions = positions_next
    velocities = velocity_next

print('global_best_position')
print(global_best_positions)
print('global_best_objective')
print(global_best_objective)
