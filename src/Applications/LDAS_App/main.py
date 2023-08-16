from pso import pso
import sys
import datetime

def main():
    # the experiment directory is passed in from the command line
    exp_dir = sys.argv[1]
    # the location of the file for the PSO
    truth_fname = (
        '/shared/pso/step_3_process_fluxcom/outputs/' +
        'le_truth_fluxcom_rs_meteo_ensemble_watts_per_m2_2001-01-01_' +
        '2009-12-31_selected_tiles.csv'
    )
    # date range over which to calculate the objective function
    start = datetime.date(2001,1,1)
    end = datetime.date(2001,1,3)
    # number of iterations with no change to global best before convergence
    conv_thresh = 5
    # number of iterations at which to declare does not converge
    max_iter = 1000
    # how often do we submit a new job?
    restart_every = 1
    # how many particles are we running?
    num_particles = 10
    # how many parameters for each particle?:
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
    num_params = 10
    # what are the proper ranges for the different paramters?
    param_range = [
        [0.1,4], # parameter 1
        [0.1,4], # parameter 2
        [0.1,4], # parameter 3
        [0.1,4], # parameter 4
        [0.1,4], # parameter 5
        [2,6],   # parameter 6
        [3,7],   # parameter 7
        [2,5],   # parameter 8
        [0.5,3], # parameter 9
        [.05,.3] # parameter 10
    ]
    # PSO hyper-parameters
    # intertia; weights the velocity at the previous step compared to 
    # velocity produced by local and global best positions
    w = 1
    # individual weight; weights the velocity change based off of individual 
    # particle's best position
    c1 = 0.7
    # global weight; weights the velocity chagne based off of global 
    # particles's best position
    c2 = 1.3
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
        'truth_fname':truth_fname,
        'start':start,
        'end':end
    }

    # do the pso calculations
    this_pso = pso(params)
    convergence_code = this_pso.get_convergence()
    # return the convergence code to the bash script
    print(convergence_code)

if __name__ == '__main__':
    main()
