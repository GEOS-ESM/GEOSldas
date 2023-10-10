import numpy as np
import sys
import os
import random
import pandas as pd
import pickle
import netCDF4 as nc
import datetime
import copy
import shutil
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import pickle as pkl

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
            convergence = self.run_pso()
        # return the convergence code that was generated
        return convergence
    def check_iteration_and_initialize(self):
        # extract the relevant information from self
        exp_dir = self.params['exp_dir']
        num_particles = self.params['num_particles']
        num_params = self.params['num_params']
        param_range = self.params['param_range']
        load_previous_pso = self.params['load_previous_pso']
        save_path = self.params['save_path']
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

                all_info = {}
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
            # delete a previous save path with this name if it exists
            if os.path.exists(save_path) and os.path.isdir(save_path):
                shutil.rmtree(save_path)
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
        conv_thresh = self.params['conv_thresh']
        save_path = self.params['save_path']
        max_iterations = self.params['max_iter']
        restart_every = self.params['restart_every']
        w = self.params['w']
        c1 = self.params['c1']
        c2 = self.params['c2']
        num_particles = self.params['num_particles']
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
        local_best_positions = np.genfromtxt(
            os.path.join(
                save_location,'local_best_positions.csv'
            )
            ,delimiter=','
        )
        local_best_obj = np.genfromtxt(
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
        global_best_obj = np.genfromtxt(
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
        # set the all_info dictionary key for this iteration
        all_info['iteration_{}'.format(num_iterations)] = {}
        # check that we have completed the requisite number of spinup
        # iterations
        # logical variable that saves whether we are in a spinup iteration
        in_spinup = True
        if (num_iterations > spinup_iterations):
            in_spinup = False
            all_obj_out = self.compute_objective(num_iterations)
            obj_out = all_obj_out[0]
            obj_out_norm = all_obj_out[1]
            et_obj_out = all_obj_out[2]
            et_obj_out_norm = all_obj_out[3]
            strm_obj_out = all_obj_out[4]
            strm_obj_out_norm = all_obj_out[5]
            #obj_out = np.arange(10)
            # find the minimum objective value for checking if global best has been
            # improved
            obj_out_min_idx = np.argmin(obj_out_norm)
            obj_out_min = obj_out_norm[obj_out_min_idx]
            # update global best if it is new
            if obj_out_min < global_best_obj:
                # mark that we have a new global best
                global_best_changed = True
                # update the global best values
                global_best_positions = positions[obj_out_min_idx,:]
                global_best_obj = obj_out_min
                global_best_part = obj_out_min_idx
            else:
                # mark that we don't have a new global best
                global_best_changed = False
            # now find which particles have found new local best positions
            local_improved_idx = np.where(obj_out_norm < local_best_obj)
            # update the local best positions for these locations
            local_best_positions[local_improved_idx] = positions[
                local_improved_idx
            ]
            local_best_obj[local_improved_idx] = obj_out_norm[local_improved_idx]
            # check if global best position vector has changed
            if global_best_changed:
                iterations_without_change = 0
            else:
                iterations_without_change += 1
            # let's save the best particle from this experiment
            src_path = os.path.join(
                exp_dir,'../',str(obj_out_min_idx),
                'output/SMAP_EASEv2_M36/cat/ens0000/'
            )
            # destination to copy the files tol
            dst_path = os.path.join(
                save_path,
                'num_'+str(num_iterations),
                str(obj_out_min_idx)
            )
            destination = shutil.copytree(src_path, dst_path)
            # check if convergence has occured
            if iterations_without_change >= conv_thresh:
                is_converged = 1
            else:
                is_converged = 0
            # now let's update velocities and positions for the next iteration
            # first get our two random numbers; will be between 0.75 and 1.25
            r1 = np.random.random()/2 + 0.75
            r2 = np.random.random()/2 + 0.75
            # calculate the velocity for the next step
            velocities_next = (
                w*velocities +
                c1*r1*(local_best_positions - positions) +
                (
                    c2*r2*(
                        np.tile(global_best_positions,(num_particles,1)) -
                        positions
                    )
                )
            )
            # calculate the positions for the next step
            positions_next = positions + velocities_next
        # if we haven't completled the requisite number of spinup iterations,
        # positions and velocities will not change
        else:
            velocities_next = copy.deepcopy(velocities)
            positions_next = copy.deepcopy(positions)
        # determine value to return to lenkf.j
        if (is_converged == 1):
            convergence = 0
        elif (is_converged == 0 and (num_iterations <= max_iterations)
              and (num_without_restart < restart_every)):
            convergence = 1
        elif (is_converged == 0 and (num_iterations <= max_iterations)
              and (num_without_restart >= restart_every)):
            convergence = 2
            iteration_track['num_without_restart'] = 0
            iteration_track.to_csv(exp_dir+'/run/iteration_track.csv')
        elif (is_converged == 0 and (num_iterations > max_iterations)):
            convergence = 3
        # add the information that we are interested in to all_info for
        # safekeeping
        # let's track positions and velocities
        all_info['iteration_{}'.format(num_iterations)]['positions'] = (
            positions
        )
        all_info['iteration_{}'.format(num_iterations)]['velocities'] = (
            velocities
        )
        # let's track whether we are in spinup or not
        all_info['iteration_{}'.format(num_iterations)]['in_spinup'] = (
            in_spinup
        )
        # let's track the global best positions, objectives, and particles for
        # this iteration
        all_info['iteration_{}'.format(num_iterations)]['global_best_positions'] = (
            global_best_positions
        )
        all_info['iteration_{}'.format(num_iterations)]['global_best_objective'] = (
            global_best_obj
        )
        # track these local best positions and objectives
        all_info['iteration_{}'.format(num_iterations)]['local_best_positions'] = (
            local_best_positions
        )
        all_info['iteration_{}'.format(num_iterations)]['local_best_obj'] = (
            local_best_obj
        )
        # save all of this iteration information for safekeeping
        all_info['iteration_{}'.format(num_iterations)]['iterations'] = (
            num_iterations
        )
        all_info['iteration_{}'.format(num_iterations)]['iterations_without_restart'] = (
            num_without_restart
        )
        all_info['iteration_{}'.format(num_iterations)]['iterations_without_change'] = (
            iterations_without_change
        )
        # save all of the objective information for safekeeping
        all_info['iteration_{}'.format(num_iterations)]['obj_out'] = (
            obj_out
        )
        all_info['iteration_{}'.format(num_iterations)]['obj_out_norm'] = (
            obj_out_norm
        )
        all_info['iteration_{}'.format(num_iterations)]['et_obj_out'] = (
            et_obj_out
        )
        all_info['iteration_{}'.format(num_iterations)]['et_obj_out_norm'] = (
            et_obj_out_norm
        )
        all_info['iteration_{}'.format(num_iterations)]['strm_obj_out'] = (
            strm_obj_out
        )
        all_info['iteration_{}'.format(num_iterations)]['strm_obj_out_norm'] = (
            strm_obj_out_norm
        )
        # update the iterations information as we are now moving to the
        # next step
        num_iterations += 1
        num_without_restart += 1
        # put these in iteration_track
        iteration_track['iterations'] = [num_iterations]
        iteration_track['num_without_restart'] = [num_without_restart]
        iteration_track['iterations_without_change'] = [iterations_without_change]
        # save everything
        # some things need to be formatted before being saved
        global_best_obj = [global_best_obj]
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
            ,positions_next,delimiter=','
        )
        np.savetxt(
            os.path.join(
                exp_dir,'../','velocities.csv'
            )
            ,velocities_next,delimiter=','
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
        # return the convergence code
        return convergence
    def compute_objective(self,num_iterations):
        # get information needed from self
        exp_dir = self.params['exp_dir']
        num_particles = self.params['num_particles']
        et_truth_fname = self.params['et_truth_fname']
        stream_truth_fname = self.params['stream_truth_fname']
        et_constraint = self.params['et_constraint']
        streamflow_constraint = self.params['streamflow_constraint']
        spinup_iterations = self.params['spinup_iter']
        exp_dir = self.params['exp_dir']
        inter_info_fname = self.params['inter_info_fname']
        all_pix_fname = self.params['all_pix_fname']
        # for the times, need to convert them to how they are saved in the
        # truth datasets
        start = self.params['start'].strftime('%Y-%m-%d')
        end = self.params['end'].strftime('%Y-%m-%d')
        start_stream = int(self.params['start'].strftime('%Y%m'))
        end_stream = int(self.params['end'].strftime('%Y%m'))
        # set the timedelta. what is the timestep at whcihc Catchment-CN saves
        # model output?
        delta = datetime.timedelta(days=1)
        # intialize the array to hold objective output for each particle
        obj_out = np.zeros(num_particles)
        obj_out_norm = np.zeros(num_particles)
        et_obj_out = np.zeros(num_particles)
        et_obj_out_norm = np.zeros(num_particles)
        strm_obj_out = np.zeros(num_particles)
        strm_obj_out_norm = np.zeros(num_particles)
        # the et constraint info will give us good information about the way
        # that the model is structured (total number of days, total number of
        # tiles running) so we are going to load that truth dataset whether it
        # is a contraint or not
        # get the et truth information
        et_truth = pd.read_csv(et_truth_fname)
        # set the index and isolate only the times that we want to
        # calculate truth for
        et_truth = et_truth.set_index('time')
        et_truth = et_truth.loc[start:end]
        # make this a numpy array
        et_truth_np = np.array(et_truth)
        # get number of steps and pixels from ET truth dataset
        model_num_steps,model_num_pixels = np.shape(et_truth_np)
        # were also going to need the list of all of the pixels that we are
        # running
        tiles = pd.read_csv(all_pix_fname,header=None)
        # turn this into an np array
        tiles = np.array(tiles).astype(int)
        # make a nice np array
        tiles = tiles.T
        all_tiles = tiles[0]

        # get the truth information
        stream_truth = pd.read_csv(stream_truth_fname)
        # set the index and isolate only times that we want
        stream_truth = stream_truth.set_index('time')
        start_fmt = int(self.params['start'].strftime('%Y%m'))
        end_fmt = int(self.params['end'].strftime('%Y%m'))
        stream_truth = stream_truth.loc[start_fmt:end_fmt]
        # make this a numpy array
        stream_truth_np = np.array(stream_truth)
        # get number of catchments and number of steps from catchment truth
        # dataset
        stream_num_steps,stream_num_hucs = np.shape(stream_truth_np)
        # get model output and calculate rmse for each particle
        for part in range(num_particles):
            if et_constraint:
                # initilize the et model output array
                this_part_et = np.zeros((model_num_steps,model_num_pixels))
            if streamflow_constraint:
                this_part_stream_day = np.zeros(
                    (model_num_steps,model_num_pixels)
                )
                this_part_stream_mon = np.zeros(
                    (stream_num_steps,stream_num_hucs)
                )
            # set the base dir where the model output is stored
            base_dir = os.path.join(
                exp_dir,'../',str(part),'output/SMAP_EASEv2_M36/cat/ens0000'
            )
            # set the start date for calculating the objective function
            curr = copy.deepcopy(self.params['start'])
            # get the output from each model simulated day
            for s in range(model_num_steps):
                # get the strings needed to navigate to this day
                curr_str = curr.strftime("%Y%m%d")
                curr_yr = curr.strftime("%Y")
                curr_mon = curr.strftime("%m")
                # and put together the name for the specific file
                this_file = os.path.join(
                    base_dir,
                    'Y'+curr_yr,
                    'M'+curr_mon,
                    str(part) + '.tavg24_1d_lnd_Nt.'+curr_str+'_1200z.nc4'
                )
                # get the data
                # and assign to the correct array
                data = nc.Dataset(this_file,mode='r')
                # if we are using an ET constraint, get that data
                if et_constraint:
                    this_et = np.array(data['LHLAND'][:])
                    # save to the larger array
                    this_part_et[s,:] = this_et
                if streamflow_constraint:
                    this_runoff = np.array(data['RUNOFF'])
                    this_baseflow = np.array(data['BASEFLOW'])
                    this_strm = this_runoff + this_baseflow
                    # lets convert units here
                    # need to go from kg/m2/s to mm/day
                    this_strm = this_strm*86400 # now in mm/day
                    this_part_stream_day[s,:] = this_strm
                    # load up the intersection data that we are going to need for the
                    # weighted averaging
                    with open(inter_info_fname,'rb') as f:
                        intersection_info = pkl.load(f)
                curr += datetime.timedelta(days=1)
            if streamflow_constraint:
                # now let's average the model streamflow to monthly resolution
                # first let's get the watersheds that we're working with here
                watersheds_str = list(intersection_info.keys())
                # get the list of watersheds as integers
                watersheds = [int(k) for k in intersection_info.keys()]
                # the number of watersheds
                num_watersheds = len(watersheds)
                # let's calcualte for each watershed
                for w,wat in enumerate(watersheds):
                    # get tiles and corresponding percentages for each of these
                    # tiles
                    # get the tiles in this watershed
                    this_tiles = intersection_info[watersheds_str[w]][0]
                    # get the percent of each of these tiles in this 
                    # watershed
                    this_perc = intersection_info[watersheds_str[w]][1]
                    # index for the tiles of interest
                    this_tiles_idx = np.zeros(0,dtype=np.int64)
                    for t,ti in enumerate(this_tiles):
                        this_ti_idx = int(np.where(all_tiles == ti)[0][0])
                        this_tiles_idx = np.append(this_tiles_idx,this_ti_idx)
                    # set the date information that we need to get the files
                    curr = copy.deepcopy(self.params['start'])
                    last_month = curr.strftime('%m')
                    days_elapsed = 0
                    months_elapsed = 0
                    last_days_elapsed = 0
                    # loop over all days until we hit the endo f the period
                    while curr <= (self.params['end'] + datetime.timedelta(days=1)):
                        # get this month
                        this_month = curr.strftime('%m')
                        # if it's a new month
                        if this_month != last_month:
                            # copy the part of the record that is just this month
                            this_mon_strm = copy.deepcopy(
                                this_part_stream_day[
                                    last_days_elapsed:days_elapsed,
                                    this_tiles_idx
                                ]
                            )
                            # take the weighted average across all tiles
                            this_mon_strm_avg = np.average(
                                this_mon_strm,
                                axis=1,
                                weights=this_perc
                            )
                            # take the sum across all days since each day is
                            # already in mm/month
                            this_mon_strm_avg = np.sum(
                                this_mon_strm
                            )
                            # assign this to the correct month and watershed
                            this_part_stream_mon[months_elapsed,w] = (
                                this_mon_strm_avg
                            )
                            # this is where we stopped last time
                            last_days_elapsed = copy.deepcopy(days_elapsed)
                            last_month = curr.strftime('%m')
                            # another month has elapsed
                            months_elapsed += 1
                        # another day has elapsed
                        days_elapsed += 1
                        # update the current date
                        curr += datetime.timedelta(days=1)
            if et_constraint:
                usable_data_idx = np.where(
                    (this_part_et > 0) &
                    (np.isnan(et_truth_np) == False)
                )
                this_et_obj = np.sqrt(
                    (
                        (
                            this_part_et[usable_data_idx] -
                            et_truth_np[usable_data_idx]
                        )**2
                    ).mean()
                )
                # now we need to normalize this by average et in fluxcom
                # dataset
                avg_fluxcom_et  = np.mean(et_truth_np[usable_data_idx])
                this_et_obj_norm = this_et_obj/avg_fluxcom_et
            if streamflow_constraint:
                this_strm_obj = np.sqrt(
                    (
                        (
                            this_part_stream_mon -
                            stream_truth_np
                        )**2
                    ).mean()
                )
                # and now we need to normalize by average streamflow for
                # waterwatch
                avg_waterwatch_strm = np.mean(stream_truth_np)
                this_strm_obj_norm = this_strm_obj/avg_waterwatch_strm
            # now let's calculate the final objective value for this particle
            # based off of what we are using as constraints.
            # if only using ET or only use streamflow as constraints, then the
            # particle obj will just be hte obj from those calculations
            # if using both ET and streamflow, the obj will be normalized by
            # the obj for that constraint as recorded on the first iteration
            # for the first particle to
            # allow both streamflow and et to contribute equally to the PSO
            if et_constraint and not streamflow_constraint:
                obj_out[part] = this_et_obj
                obj_out_norm[part] = this_et_obj_norm
                et_obj_out[part] = this_et_obj
                et_obj_out_norm[part] = this_et_obj_norm
            if streamflow_constraint and not et_constraint:
                obj_out[part] = this_strm_obj
                obj_out_norm[part] = this_strm_obj_norm
                strm_obj_out[part] = this_strm_obj
                strm_obj_out_norm[part] = this_strm_obj_norm
            if et_constraint and streamflow_constraint:
                # get where we will save or load the normalize values from
                #et_norm_fname = os.path.join(
                #    exp_dir,'../','et_norm.csv'
                #)
                #stream_norm_fname = os.path.join(
                #    exp_dir,'../','stream_norm.csv'
                #)
                #if (num_iterations - spinup_iterations) == 1 and part == 0:
                #    # get the values we will use to normalize
                #    et_norm = this_et_obj
                #    stream_norm = this_strm_obj
                #    # save these so they can be loaded by future particles and
                #    # iterations
                #    # for et norm
                #    with open(et_norm_fname,'w') as f:
                #        f.write('{}'.format(et_norm))
                #    # for streamflow norm
                #    with open(stream_norm_fname,'w') as f:
                #        f.write('{}'.format(stream_norm))
                #else:
                #    # load up the normalized values
                #    with open(et_norm_fname,'r') as f:
                #        et_norm = f.read()
                #    et_norm = float(et_norm)
                #    with open(stream_norm_fname,'r') as f:
                #        stream_norm = f.read()
                #    stream_norm = float(stream_norm)
                ## normalize with these values
                #et_obj_norm = this_et_obj/et_norm
                #stream_obj_norm = this_strm_obj/stream_norm
                # add them together so that we get the final obj
                et_weight = self.params['obj_weights'][0]
                strm_weight = self.params['obj_weights'][1]
                this_tot_obj = this_et_obj + this_strm_obj
                this_tot_obj_norm = (
                    et_weight*this_et_obj_norm +
                    strm_weight*this_strm_obj_norm
                )
                obj_out[part] = this_tot_obj
                obj_out_norm[part] = this_tot_obj_norm
                et_obj_out[part] = this_et_obj
                et_obj_out_norm[part] = this_et_obj_norm
                strm_obj_out[part] = this_strm_obj
                strm_obj_out_norm[part] = this_strm_obj_norm
        # return all objective values
        out_vals = [
            obj_out,obj_out_norm,et_obj_out,et_obj_out_norm,
            strm_obj_out,strm_obj_out_norm
        ]
        return out_vals
    def test_pso(self):
        '''
        To test the PSO, simply run an optimization that only runs the
        model for a couple of days.
        To debug the PSO, simply run the first iteration of such an
        optimization. Then this first output can be used to debug runtime
        errors. Errors in logic can then be addressed using the testing method
        listed above.
        '''
        pass
    def plot_pso_progress(self):
        '''
        What plots would be helpful?
        For each parameter and after each iteration, plot the location of
        the particles with each iteration.
        For each iteration, plot the movement of that particles' objective
        output.
        For each iteration, plot the movement of the global best objective
        value.
        For each iteration, plot each particles' local best objective value.

        '''
        
