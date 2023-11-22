from pso import pso
import sys
import datetime

def main():
    # the experiment directory is passed in from the command line
    exp_dir = sys.argv[1]
    # location to save outputs of new global best model runs
    save_path = '/shared/pso_outputs/fake'
    # should Fluxcom ET be used as a constraint in objective function?
    et_constraint = True
    # should USGS WaterWatch runoff be used as a constraint in objective
    # function?
    streamflow_constraint = True
    # the location of the ET truth file for the PSO
    et_truth_fname = (
        '/shared/pso/step_3_process_fluxcom/outputs/' +
        'le_truth_fluxcom_rs_CRUNCEP_ensemble_watts_per_m2_2001-01-01_' +
        '2015-12-31_selected_tiles.csv'
    )
    # the location of the streamflow truth file for the PSO
    stream_truth_fname = (
        '/shared/pso/step_3.1_process_streamflow/outputs/' +
        'streamflow_data_mm_per_month_200101_201512.csv'
    )
    # location of the intersection info .pkl file that relates catchment tiles
    # to huc watershed
    inter_info_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersection_info.pkl'
    )
    # location where list of all pixels being run is stored
    all_pix_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    # date range over which to calculate the objective function
    # this is inclusive
    start = datetime.date(2002,1,1)
    end = datetime.date(2002,12,31)
    # number of iterations with no change to global best before convergence
    conv_thresh = 1
    # number of iterations at which to declare does not converge
    max_iter = 1
    # how often do we submit a new job?
    restart_every = 1
    # how many particles are we running?
    num_particles = 1
    # how many parameters for each particle?:
    num_params = 10
    # what are the names of the different paramaeters?
        #parameter 1: aj in g1 EF for forests
        #parameter 2: aj in g1 EF for croplands
        #parameter 3: aj in g1 EF for grasslands
        #parameter 4: aj in g1 EF for savannas
        #parameter 5: aj in g1 EF for shrublands
        #parameter 6: alpha in Ksat EF
        #parameter 7: beta in Ksat EF
        #parameter 8: constant_1 in Ksat EF
        #parameter 9: constant_2 in Ksat EF
        #parameter 10: sand_exp in Ksat EF
    parameter_names = [
        'a1_forests',
        'a1_croplands',
        'a1_grasslands',
        'a1_savannas',
        'a1_shrublands',
        'a0_forests',
        'a0_croplands',
        'a0_grasslands',
        'a0_savannas',
        'a0_shrublands'
    ]
    # what are the starting ranges for the different paramters?
    param_range = [
        [0,15], # parameter 1
        [0,15], # parameter 2
        [0,15], # parameter 3
        [0,15], # parameter 4
        [0,15], # parameter 5
        [-.25,.25],   # parameter 6
        [-.25,.25],   # parameter 7
        [-.25,.25],   # parameter 8
        [-.25,.25], # parameter 9
        [-.25,.25] # parameter 10
    ]
    # what are the objective function weights? This is multiplied by the
    # normalized objective output for [0]:et and [1]:streamflow. These two
    # weights should sum to one
    obj_weights = [0.5,0.5]
    # PSO hyper-parameters
    # intertia; weights the velocity at the previous step compared to 
    # velocity produced by local and global best positions
    w = 1
    # individual weight; weights the velocity change based off of individual 
    # particle's best position
    c1 = 0.5
    # global weight; weights the velocity chagne based off of global 
    # particles's best position
    c2 = 0.2
    # how many spinup iterations should we do before changing initial
    # parameters?
    spinup_iter = 0
    # where should we write PSO-specific text updates?
    text_out = 'pso_out.txt'
    # are we starting from a previous PSO run?
    load_previous_pso = False
    # if we are starting from a previous PSO run, what is the base directory of
    # this?
    previous_pso_dir = ''

    # compile all these parameters into one dictionary to pass to pso script
    params = {
        'exp_dir':exp_dir,
        'conv_thresh':conv_thresh,
        'max_iter':max_iter,
        'restart_every':restart_every,
        'num_particles':num_particles,
        'num_params':num_params,
        'param_range':param_range,
        'w':w,
        'c1':c1,
        'c2':c2,
        'spinup_iter':spinup_iter,
        'text_out':text_out,
        'load_previous_pso':load_previous_pso,
        'previous_pso_dir':previous_pso_dir,
        'et_truth_fname':et_truth_fname,
        'stream_truth_fname':stream_truth_fname,
        'start':start,
        'end':end,
        'save_path':save_path,
        'et_constraint':et_constraint,
        'streamflow_constraint':streamflow_constraint,
        'param_names':parameter_names,
        'inter_info_fname':inter_info_fname,
        'all_pix_fname':all_pix_fname,
        'obj_weights':obj_weights
    }

    # do the pso calculations
    this_pso = pso(params)
    convergence_code = this_pso.get_convergence()
    # return the convergence code to the bash script
    print(convergence_code)

if __name__ == '__main__':
    main()
