import numpy as np
import sys
import os
import random
import pandas as pd
import pickle
import netCDF4 as nc

class pso:
    def __init__(self,params):
        '''
        pass the needed parameters into init and save to self
        '''
        self.params = params
    def get_convergence(self):
        '''
        function that drives the PSO class. runs all of the following functions
        in this class.
        '''
        first_iteration = self.check_iteration_and_initialize()
        # if this was the first iteration we are done here
        if first_iteration:
            # return a dummy convergence code to keep the script running
            convergence = 1
        # otherwise continue on, calculate objective, etc.
        else:
            self.run_pso()
            sys.exit()
        # return the convergence code that was generated
        return convergence
    def check_iteration_and_initialize(self):
        # extract the relevant information from self
        exp_dir = self.params['exp_dir']
        num_particles = self.params['num_particles']
        num_params = self.params['num_particles']
        param_range = self.params['param_range']
        load_previous_pso = self.params['load_previous_pso']
        # check to see if a position_vals.csv file exists
        exists = os.path.exists(
            os.path.join(exp_dir,'../','positions.csv')
        )
        # if it does not exist, then we need to initialize the files
        if exists == False:
            # if we're not going to load from a previous PSO experiment, let's
            # initialize everything
            if load_previous_pso == False:
                # let's initialize particles as randome within the range
                # and let's intialize velocity as a random number not to exceed
                # 1/10
                # of the range of the particle
                positions = np.zeros((num_particles,num_params))
                velocities = np.zeros((num_particles,num_params))
                # let's do this for each particle
                for p in range(num_particles):
                    # let's do this for each parameter in this particle
                    for r,ran in enumerate(param_range):
                        positions[p,r] = random.uniform(ran[0],ran[1])
                        full_range = ran[1] - ran[0]
                        acceptable_range = full_range/10
                        velocities[p,r] = random.uniform(
                            -acceptable_range,acceptable_range
                        )
                # we need local best positions, which at this point are each
                # individual's best positions
                local_best_positions = positions
                # we also need global best positions, which we don't know yet, so
                # we will set to infinity
                global_best_positions = np.zeros(num_params)
                global_best_positions[:] = float('inf')
                # we need to assign the local and global best objectives as
                # infinity since we don't know these yet
                local_best_obj = np.zeros(num_particles)
                local_best_obj[:] = float('inf')
                global_best_obj = [float('inf')]
                # initialize the iteration counters
                num_iterations = 1
                num_without_restart = 1
                iterations_without_change = 0
                # create the df that holds all of the iteration information
                iteration_track = pd.DataFrame({})
                iteration_track['iterations'] = [num_iterations]
                iteration_track['num_without_restart'] = [num_without_restart]
                iteration_track['iterations_without_change'] = [
                    iterations_without_change
                ]

                all_info = {
                    'iteration_{}'.format(num_iterations):{
                        'positions':positions,
                        'velocities':velocities,
                        'local_best_positions':local_best_positions,
                        'global_best_positions':global_best_positions,
                        'local_best_obj':local_best_obj,
                        'global_best_obj':global_best_obj,
                        'iteration_track':iteration_track
                    }
                }
                print(positions)
                print(velocities)
                print(iteration_track)
            # if we are loading from another experiment, load those pso values
            else:
                with open(os.path.join(pso_restart_dir,'all_info.pkl'),'rb') as f:
                    all_info = pickle.load(f)

                iteration_track = pd.read_csv(
                    os.path.join(
                        pso_restart_dir,'iteration_track.csv'
                    )
                )
                positions = np.genfromtxt(
                    os.path.join(
                        pso_restart_dir,'positions.csv'
                    )
                    ,delimiter=','
                )
                velocities = np.genfromtxt(
                    os.path.join(
                        pso_restart_dir,'velocities.csv'
                    )
                    ,delimiter=','
                )
                best_positions = np.genfromtxt(
                    os.path.join(
                        pso_restart_dir,'local_best_positions.csv'
                    )
                    ,delimiter=','
                )
                best_positions_objective = np.genfromtxt(
                    os.path.join(
                        pso_restart_dir,'local_best_obj.csv'
                    )
                    ,delimiter=','
                )
                global_best_positions = np.genfromtxt(
                    os.path.join(
                        pso_restart_dir,'global_best_positions.csv'
                    )
                    ,delimiter=','
                )
                global_best_objective = np.genfromtxt(
                    os.path.join(
                        pso_restart_dir,'global_best_obj.csv'
                    )
                    ,delimiter=','
                )
            # now let's save these so that they can be used to run the model
            with open(os.path.join(exp_dir,'../','all_info.pkl'),'wb') as f:
                pickle.dump(all_info,f)

            iteration_track.to_csv(
                os.path.join(
                    exp_dir,'../','iteration_track.csv'
                )
            )
            np.savetxt(
                os.path.join(
                    exp_dir,'../','positions.csv'
                )
                ,positions,delimiter=','
            )
            np.savetxt(
                os.path.join(
                    exp_dir,'../','velocities.csv'
                )
                ,velocities,delimiter=','
            )
            np.savetxt(
                os.path.join(
                    exp_dir,'../','local_best_positions.csv'
                )
                ,local_best_positions,delimiter=','
            )
            np.savetxt(
                os.path.join(
                    exp_dir,'../','local_best_obj.csv'
                )
                ,local_best_obj,delimiter=','
            )
            np.savetxt(
                os.path.join(
                    exp_dir,'../','global_best_positions.csv'
                )
                ,global_best_positions,delimiter=','
            )
            np.savetxt(
                os.path.join(
                    exp_dir,'../','global_best_obj.csv'
                )
                ,global_best_obj,delimiter=','
            )
            # confirm that this was the first iteration
            first_iteration = True
        else:
            # go further because this was not the first iteration
            first_iteration = False
        # return information on the first iteration so that the model knows
        # whether to continue running or if it has done it's job
        return first_iteration
    def run_pso(self):
        # extract the relevant information from self
        exp_dir = self.params['exp_dir']
        spinup_iterations = self.params['spinup_iter']
        # load the information that we need
        # first define where the pso information has been saved
        save_location = os.path.join(
            exp_dir,'../'
        )
        with open(os.path.join(save_location,'all_info.pkl'),'rb') as f:
            all_info = pickle.load(f)

        iteration_track = pd.read_csv(
            os.path.join(
                save_location,'iteration_track.csv'
            )
        )
        positions = np.genfromtxt(
            os.path.join(
                save_location,'positions.csv'
            )
            ,delimiter=','
        )
        velocities = np.genfromtxt(
            os.path.join(
                save_location,'velocities.csv'
            )
            ,delimiter=','
        )
        best_positions = np.genfromtxt(
            os.path.join(
                save_location,'local_best_positions.csv'
            )
            ,delimiter=','
        )
        best_positions_objective = np.genfromtxt(
            os.path.join(
                save_location,'local_best_obj.csv'
            )
            ,delimiter=','
        )
        global_best_positions = np.genfromtxt(
            os.path.join(
                save_location,'global_best_positions.csv'
            )
            ,delimiter=','
        )
        global_best_objective = np.genfromtxt(
            os.path.join(
                save_location,'global_best_obj.csv'
            )
            ,delimiter=','
        )
        # get all of our interation information
        num_iterations = iteration_track['iterations'].loc[0]
        num_without_restart = iteration_track['num_without_restart'].loc[0]
        iterations_without_change = (
            iteration_track['iterations_without_change'].loc[0]
        )
        # check that we have completed the requisite number of spinup
        # iterations
        # logical variable that saves whether we are in a spinup iteration
        in_spinup = False
        if (num_iterations < spinup_iterations):
            in_spinup = True
            num_iterations += 1
            num_without_restart += 1
            iterations_without_change = 1
            is_converged = 0
            # update iteration tracker 
            iteration_track['iterations'] = [num_iterations]
            iteration_track['num_without_restart'] = [num_without_restart]
            iteration_track['iterations_without_change'] = [
                iterations_without_change
            ]
        # if not in spinup, change particle positions according to the pso
        obj_out = self.compute_objective()
        print('done')
        sys.exit()



        ## copy the output files to a new location so that they can 
        ## be analyzed if desired. otherwise they will be overwritten in the
        ## next iteration of the PSO
        ## only do this if in spinup iterations or if this was a new best result
        #if iterations_without_change == 0 or in_spinup == True:
        #    for par in range(nparticles):
        #        # name of the experiment
        #        exp_name = os.path.split(os.path.split(exp_dir)[0])[1]
        #        # path where the catchment output files are currently stored
        #        src_path = os.path.join(
        #            exp_dir,'../',str(global_best_part),
        #            'output/SMAP_EASEv2_M36/cat/ens0000/'
        #        )
        #        # destination to copy the files tol
        #        dst_path = os.path.join(
        #            '/discover/nobackup/projects/medComplex/pso_outputs',
        #            exp_name,'num_'+str(num_iterations),
        #            str(global_best_part)
        #        )
        #        if os.path.exists(dst_path) and os.path.isdir(dst_path):
        #            shutil.rmtree(dst_path)
        #        destination = shutil.copytree(src_path, dst_path)


        ## determine value to return to lenkf.j
        #if (is_converged == 1):
        #    convergence = 0
        #elif (is_converged == 0 and (num_iterations <= max_iterations)
        #      and (num_without_restart < restart_every)):
        #    convergence = 1
        #elif (is_converged == 0 and (num_iterations <= max_iterations)
        #      and (num_without_restart >= restart_every)):
        #    convergence = 2
        #    iteration_track['num_without_restart'] = 0
        #    iteration_track.to_csv(exp_dir+'/run/iteration_track.csv')
        #elif (is_converged == 0 and (num_iterations > max_iterations)):
        #    convergence = 3

    def compute_objective(self):
        # get information needed from self
        exp_dir = self.params['exp_dir']
        num_particles = self.params['num_particles']
        truth_fname = self.params['truth_fname']
        # for the times, need to convert them to how they are saved in the
        # truth datasets
        start = self.params['start'].strftime('%Y-%m-%d')
        end = self.params['end'].strftime('%Y-%m-%d')
        # set the timedelta. what is the timestep at whcihc Catchment-CN saves
        # model output?
        delta = datetime.timedelta(days=1)
        # intialize the array to hold objective output for each particle
        obj_out = np.zeros(num_particles)
        # get the truth information
        truth_data = pd.read_csv(truth_fname)
        # set the index and isolate only the times that we want to look over
        truth_data = truth_data.set_index('time')
        truth_data = truth_data.loc[start:end]
        # save as a numpy array
        truth_np = np.array(truth_data)
        # get the number of steps and number of pixels from the size of the
        # truth dataset
        num_steps,num_pixels = np.shape(truth_np)
        print(truth_data)
        print(num_steps)
        print(num_pixels)
        # get model output and calculate rmse for each particle
        for part in range(num_particles):
            # initiatilze array to hold et
            this_part_et = np.zeros((num_steps,num_pixels))
            # set the base dir where the model output is stored
            base_dir = os.path.join(
                exp_dir,'../',str(part),'output/SMAP_EASEv2_M36/cat/ens0000'
            )
            # get the output from each model simulated day
            for s in range(num_steps):
                # get the strings needed to navigate to this day
                curr_str = curr.strftime("%Y%m%d")
                curr_yr = curr.strftime("%Y")
                curr_mon = curr.strftime("%m")
                # and put together the name for the specific file
                this_file = os.path.join(
                    root,
                    'Y'+curr_yr,
                    'M'+curr_mon,
                    str(ens) + '.tavg24_1d_lnd_Nt.'+curr_str+'_1200z.nc4'
                )
                # get the data
                # and assign to the correct array
                data = nc.Dataset(this_file,mode='r')
                this_et = np.array(data['LHLAND'][:])
                # save to the larger array
                this_part_et[s,:] = this_et
                # update the time for the next step
                curr += delta
            # we only want to compute the objective function where the model
            # thinks there is positive ET, so we will get the idx here
            neg_idx = np.where(this_part_et < 0)
            # calculate the objective value for this particle
            this_obj = np.sqrt(
                (
                    (
                        this_part_et[neg_idx] - truth_np[neg_idx]
                    )**2
                ).mean()
            )
    def save_as_csv(self,data,fname):
        pass
    def test_pso(self):
        pass
    def plot_particle_movement(self):
        pass
