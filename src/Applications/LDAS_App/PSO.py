# script that will run the PSO algorithm

# import libraries
from os.path import exists
import pandas as pd
import numpy as np
import os

# set important parameters
convergence_threshold = 10 #num iterations with no change to global_best to converge
max_iterations = 1000
restart_every = 500
nparticles = 28 #fits architecture of NCCS?
nparameters = 5
#parameter 1: g1 for forests
#parameter 2: g1 for croplands
#parameter 3: g1 for grasslands
#parameter 4: g1 for savannas
#parameter 5: g1 for shrublands

# find out where we are in # of iterations
if exists('iteration_track.csv') == True:
    iteration_track = pd.read_csv('iteration_track.csv')
    num_iterations = iteration_track['iterations']
    num_without_restart = iteration_track['num_without_restart']
    iterations_without_change = iteration_track['iterations_without_change']
else:
    num_iterations = 0
    num_without_restart = 0
    iterations_without_change = 0
    iteration_track = pd.DataFrame({})

# do PSO calculations
#PSO parameters
w = 0.5 #SEE WHAT WAS USED PREVIOUSLY; inertia weight that slows down the velocity of particle at previous step
c1 = 1 #SEE WHAT WAS USED PREVIOUSLY; weight put on the individual particle's best position (higher=more exploration)
c2 = 1 #SEE WHAT WAS USED PREVIOUSLY; weight put on the global best position (higher=faster convergence)\
r1 = np.random.random() #should this random number be between 0 and 1 or should I scale?
r2 = np.random.random() #should this random number be between 0 and 1 or should I scale?


#check if we need to initialize
if exists('position_vals.csv') == True:
    positions = np.genfromtxt('position_vals.csv',delimiter=',')
    velocities = np.genfromtxt('velocity_vals.csv',delimiter=',')
    best_position = np.genfromtxt('best_position.csv',delimiter=',')
    global_best_positions = np.genfromtxt('global_best_positions.csv',delimiter=',')
    global_best_objective = np.genfromtxt('global_best_objective.csv',delimiter=',')
else:
    #position needs to be updated to span all expected values of g1. Currently, everything is between 0 and 1 (random).
    #velocity needs to be updated to ?????
    positions = np.random.random((nparticles,nparameters))
    velocities = np.random.random((nparticles,nparameters))
    best_position = positions
    objective_result = objective_function(positions) #need to define this function
    global_best_idx = np.argmin(objective_result) #make sure that this axis is correct
    assert len(global_best_idx) == nparameters
    global_best = np.tile(positions(global_best_idx,:))

#update velocity
velocity_next = w*velocities + c1*r2*(best_position - positions) + c2*r2*(global_best - positions)

#update positions
positions_next = positions + velocity_next

#evaluate the objective function here (still need to work on how to do that)
#will result in an array that is objective_result = np_array(nparticles,1)
objective_result = objective_function(positions_next) #need to write this function

#update global best position
global_objective_result = objective_function(global_best)
objective_result_min_idx = np.argmin(objective_result)
objective_result_min = positions_next(objective_result_min_idx)
global_objective_result_min = global_objective_result[0,:]
global_best_next = np.copy(global_best)
for param in len(objective_result_min):
    if objective_result_min[param] < global_objective_result_min[param]:
        global_best_next[param] = positions_next[objective_result_min_idx[param]]

#update local best position
local_objective_result = objective_function(best_position)
local_improved_idx = np.where(objective_result < local_objective_result)
best_position_next = np.copy(best_position)
best_position_next[local_improved_idx] = positions_next[local_improved_idx]

#save the updated arrays to their respective .csv files to be called in the next iteration
np.savetxt('position_vals.csv',positions_next,delimiter=',')
np.savetxt('velocity_vals.csv',velocity_next,delimiter=',')
np.savetxt('best_position.csv',best_position_next,delimiter=',')
np.savetxt('global_best.csv',global_best_next,delimiter=',')

#check if global best position vector has changed
if global_best_next == global_best:
    iterations_without_change += 1
else:
    iterations_without_change == 0

#check if convergence has occured
if iterations_without_change == convergence_threshold:
    is_converged = 1
else
    is_converged = 0

#write to PSO_tracker.csv to keep track of variables
num_iterations += 1
num_without_restart += 1
# more will go here as I figure out what PSO values we will need in the next iteration....

iteration_track['iterations'] = num_iterations
iteration_track['num_without_restart'] = num_without_restart
iteration_track['iterations_without_change'] = iterations_without_change
iteration_track.to_csv('iteration_track.csv')

# determine value to return to lenkf.j
if (is_converged == 1):
    convergence = 0
    #print(convergence)
elif (is_converged == 0 and (num_iterations <= max_iterations)
      and (num_without_restart < restart_every)):
    convergence = 1
    #print(convergence)
elif (is_converged == 0 and (num_iterations <= max_iterations)
      and (num_without_restart >= restart_every)):
    convergence = 2
    iteration_track['num_without_restart'] = 0
    iteration_track.to_csv('iteration_track.csv')
    #print(convergence)
elif (is_converged == 0 and (num_iterations > max_iterations):
    convergence = 3
    #print(convergence)


#assign the ensemble members their updated parameter values convergence != 3:
    # METHOD TO ASSIGN
    # do this later
    hello = 1


# return the convergence value
print(convergence)
