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
from compute_objective_v2 import compute_objective


exp_dir = sys.argv[1]

# set important parameters
convergence_threshold = 10 #num iterations with no change to global_best to converge
max_iterations = 1000
restart_every = 500
nparticles = 2 #fits architecture of NCCS?
nparameters = 5
#parameter 1: j for forests
#parameter 2: j for croplands
#parameter 3: j for grasslands
#parameter 4: j for savannas
#parameter 5: j for shrublands

# do PSO calculations
#PSO parameters
w = 0.5 #SEE WHAT WAS USED PREVIOUSLY; inertia weight that slows down the velocity of particle at previous step
c1 = 1 #SEE WHAT WAS USED PREVIOUSLY; weight put on the individual particle's best position (higher=more exploration)
c2 = 1 #SEE WHAT WAS USED PREVIOUSLY; weight put on the global best position (higher=faster convergence)\
r1 = np.random.random() #should this random number be between 0 and 1 or should I scale?
r2 = np.random.random() #should this random number be between 0 and 1 or should I scale?

#check if we need to initialize
#check_exists = exp_dir+'run/position_vals.csv'
if exists(exp_dir+'/run/position_vals.csv') == False:
    #position needs to be updated to span all expected values of g1. Currently, everything is between 0 and 1 (random).
    #velocity needs to be updated to ?????
    positions = np.random.random((nparticles,nparameters))*2 #all random values between 0 and 2
    velocities = np.random.random((nparticles,nparameters))

    #save positions and velocities and iteration data here
    num_iterations = 0
    num_without_restart = 0
    iterations_without_change = 0
    iteration_track = pd.DataFrame({})

    iteration_track['iterations'] = [num_iterations]
    iteration_track['num_without_restart'] = [num_without_restart]
    iteration_track['iterations_without_change'] = [iterations_without_change]
    save_to = exp_dir+'run/iteration_track.csv'
    #print(exp_dir)
    #print(save_to)
    iteration_track.to_csv(exp_dir+'/run/iteration_track.csv')

    np.savetxt(exp_dir+'/run/position_vals.csv',positions,delimiter=',')
    np.savetxt(exp_dir+'/run/velocity_vals.csv',velocities,delimiter=',')

    #return convergence code to submit job
    convergence = 1
    with open(exp_dir+'/run/output.txt','a') as f:
        f.write('num_iterations')
        f.write('\n')
        f.write(str(num_iterations))
        f.write('\n')
        f.write('positions')
        f.write('\n')
        f.write(str(positions))
        f.write('\n')
        f.write('velocities')
        f.write('\n')
        f.write(str(velocities))
        f.write('\n')
    print(convergence)

elif ((exists(exp_dir+'/run/position_vals.csv') == True) and (exists(exp_dir+'/run/best_position_objective.csv') == False)):
    iteration_track = pd.read_csv(exp_dir+'/run/iteration_track.csv')
    num_iterations = iteration_track['iterations'].loc[0]
    num_without_restart = iteration_track['num_without_restart'].loc[0]
    iterations_without_change = iteration_track['iterations_without_change'].loc[0]

    positions = np.genfromtxt(exp_dir+'/run/position_vals.csv',delimiter=',')
    velocities = np.genfromtxt(exp_dir+'/run/velocity_vals.csv',delimiter=',')
    best_position = np.copy(positions)
    objective_result = compute_objective(exp_dir,nparticles) #returns a np array from the different particles
    best_position_objective = np.copy(objective_result)
    #print(objective_result)
    global_best_idx = np.argmin(objective_result)
    #print(global_best_idx)
    global_best_positions = positions[global_best_idx,:]
    #print(global_best_positions)
    global_best_objective = objective_result[global_best_idx]
    #print(global_best_objective)

    with open(exp_dir+'/run/output.txt','a') as f:
        f.write(str(objective_result))
        f.write('\n')
        f.write(str(global_best_idx))
        f.write('\n')
        f.write(str(global_best_objective))
        f.write('\n')


    #update velocity
    velocity_next = (w*velocities + c1*r2*(best_position - positions) + 
                     c2*r2*(np.tile(global_best_positions,(nparticles,1)) - 
                            positions))

    #update positions
    positions_next = positions + velocity_next

    #save all data here
    #save the updated arrays to their respective .csv files to be called in the next iteration
    global_best_objective = [global_best_objective]
    np.savetxt(exp_dir+'/run/position_vals.csv',positions_next,delimiter=',')
    np.savetxt(exp_dir+'/run/velocity_vals.csv',velocity_next,delimiter=',')
    np.savetxt(exp_dir+'/run/best_position.csv',best_position,delimiter=',')
    np.savetxt(exp_dir+'/run/best_position_objective.csv',best_position_objective,delimiter=',')
    np.savetxt(exp_dir+'/run/global_best_positions.csv',global_best_positions,delimiter=',')
    np.savetxt(exp_dir+'/run/global_best_objective.csv',global_best_objective,delimiter=',')

    iteration_track['iterations'] = [num_iterations]
    iteration_track['num_without_restart'] = [num_without_restart]
    iteration_track['iterations_without_change'] = [iterations_without_change]
    iteration_track.to_csv(exp_dir+'/run/iteration_track.csv')

    #return convergence code to submit job
    convergence = 1
    with open(exp_dir+'/run/output.txt','a') as f:
        f.write('num_iterations')
        f.write('\n')
        f.write(str(num_iterations))
        f.write('\n')
        f.write('positions_next')
        f.write('\n')
        f.write(str(positions_next))
        f.write('\n')
        f.write('velocity_next')
        f.write('\n')
        f.write(str(velocity_next))
        f.write('\n')
        f.write('best_position')
        f.write('\n')
        f.write(str(best_position))
        f.write('\n')
        f.write('best_position_objective')
        f.write('\n')
        f.write(str(best_position_objective))
        f.write('\n')
        f.write('global_best_positions')
        f.write('\n')
        f.write(str(global_best_positions))
        f.write('\n')
        f.write('global_best_objective')
        f.write('\n')
        f.write(str(global_best_objective))
    print(convergence)

else:
    #run as layed out previuosly
    iteration_track = pd.read_csv(exp_dir+'/run/iteration_track.csv')
    num_iterations = iteration_track['iterations'].loc[0]
    num_without_restart = iteration_track['num_without_restart'].loc[0]
    iterations_without_change = iteration_track['iterations_without_change'].loc[0]

    positions = np.genfromtxt(exp_dir+'/run/position_vals.csv',delimiter=',')
    velocities = np.genfromtxt(exp_dir+'/run/velocity_vals.csv',delimiter=',')
    best_position_last = np.genfromtxt(exp_dir+'/run/best_position.csv',delimiter=',')
    best_position_objective_last = np.genfromtxt(exp_dir+'/run/best_position_objective.csv',delimiter=',')
    global_best_positions_last = np.genfromtxt(exp_dir+'/run/global_best_positions.csv',delimiter=',')
    global_best_objective_last = np.genfromtxt(exp_dir+'/run/global_best_objective.csv',delimiter=',')


    #objective function
    objective_result = compute_objective(exp_dir,nparticles) #objective function to be filled in

    with open(exp_dir+'/run/output.txt','a') as f:
        f.write('global_best_objective_last')
        f.write(str(global_best_objective_last))
        f.write('\n')
        f.write(str(np.shape(global_best_objective_last)))
        f.write(str(type(global_best_objective_last)))
        f.write('\n')
        f.write(str(objective_result))

    #best global position
    #update global best position
    objective_result_min_idx = np.argmin(objective_result)
    objective_result_min = objective_result[objective_result_min_idx]
    
    global_objective_result_min_last = global_best_objective_last

    global_best_positions = np.copy(global_best_positions_last)
    global_best_objective = np.copy(global_best_objective_last)
    #print(objective_result_min)
    #print(global_objective_result_min_last)
    if objective_result_min < global_objective_result_min_last:
        global_best_positions = positions[objective_result_min_idx,:]
        global_best_objective = objective_result_min
        #print('global best updated')
        #print('old best positions global')
        #print(global_best_positions_last)
        #print('new best positions global')
        #print(global_best_positions)
        #print('old best objective global')
        #print(global_best_objective_last)
        #print('new best objective global')
        #print(global_best_objective)


    #best local position
    #update local best position
    #local_objective_result = np.random.random(nparticles) #objective_function(best_position)
    local_improved_idx = np.where(objective_result < best_position_objective_last)
    best_position = np.copy(best_position_last)
    best_position_objective = np.copy(best_position_objective_last)
    best_position[local_improved_idx] = positions[local_improved_idx]
    best_position_objective[local_improved_idx] = objective_result[local_improved_idx]
    #print('old best positions local')
    #print(best_position_last)
    #print('new best positions local')
    #print(best_position)
    #print('old best objective local')
    #print(best_position_objective_last)
    #print('new best objective local')
    #print(best_position_objective)

    #update velocity
    velocity_next = (w*velocities + c1*r2*(best_position - positions) +
                     c2*r2*(np.tile(global_best_positions,(nparticles,1)) -
                            positions))

    #update positions
    positions_next = positions + velocity_next

    #save the updated arrays to their respective .csv files to be called in the next iteration
    global_best_objective = [global_best_objective]
    np.savetxt(exp_dir+'/run/position_vals.csv',positions_next,delimiter=',')
    np.savetxt(exp_dir+'/run/velocity_vals.csv',velocity_next,delimiter=',')
    np.savetxt(exp_dir+'/run/best_position.csv',best_position,delimiter=',')
    np.savetxt(exp_dir+'/run/best_position_objective.csv',best_position_objective,delimiter=',')
    np.savetxt(exp_dir+'/run/global_best_positions.csv',global_best_positions,delimiter=',')
    np.savetxt(exp_dir+'/run/global_best_objective.csv',global_best_objective,delimiter=',')

    #check if global best position vector has changed
    #print(global_best_positions_last)
    #print(global_best_positions)
    if np.all(global_best_positions_last == global_best_positions):
        iterations_without_change += 1
    else:
        iterations_without_change = 0

    #print('iterations without change')
    #print(iterations_without_change)
    #print('convergence threshold')
    #print(convergence_threshold)
    #print('num iterations')
    #print(num_iterations)
    #print('num no restart')
    #print(num_without_restart)

    #check if convergence has occured
    if iterations_without_change == convergence_threshold:
        is_converged = 1
    else:
        is_converged = 0

    #write to PSO_tracker.csv to keep track of variables
    num_iterations += 1
    num_without_restart += 1
    # more will go here as I figure out what PSO values we will need in the next iteration....

    iteration_track['iterations'] = [num_iterations]
    iteration_track['num_without_restart'] = [num_without_restart]
    iteration_track['iterations_without_change'] = [iterations_without_change]
    iteration_track.to_csv(exp_dir+'/run/iteration_track.csv')

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
        iteration_track.to_csv(exp_dir+'/run/iteration_track.csv')
        #print(convergence)
    elif (is_converged == 0 and (num_iterations > max_iterations)):
        convergence = 3
        #print(convergence)


    with open(exp_dir+'/run/output.txt','a') as f:
        f.write('num_iterations')
        f.write('\n')
        f.write(str(num_iterations))
        f.write('\n')
        f.write('positions_next')
        f.write('\n')
        f.write(str(positions_next))
        f.write('\n')
        f.write('velocity_next')
        f.write('\n')
        f.write(str(velocity_next))
        f.write('\n')
        f.write('best_position')
        f.write('\n')
        f.write(str(best_position))
        f.write('\n')
        f.write('best_position_objective')
        f.write('\n')
        f.write(str(best_position_objective))
        f.write('\n')
        f.write('global_best_positions')
        f.write('\n')
        f.write(str(global_best_positions))
        f.write('\n')
        f.write('global_best_objective')
        f.write('\n')
        f.write(str(global_best_objective))

    print(convergence)
