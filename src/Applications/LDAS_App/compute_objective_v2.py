import numpy as np
import pandas as pd
import glob
from netCDF4 import Dataset

def compute_objective(exp_dir,num_ens):
    # THIS WILL COMPUTE THE OBJECTIVE FUNCTION WHEN:
        ## the output is daily data
        ## the observations are daily data
    objective_output = np.zeros(num_ens)
    #flux_et = pd.read_csv(exp_dir+'/output/flux_tower_et.csv')
    flux_et = pd.read_csv(
        '../../../misc_and_testing/create_env_data/fluxcluster_et_final_19800101_19800102_TESTING.csv'
    )
    flux_et = flux_et.drop('time',axis=1)
    flux_et = np.array(flux_et)
    flux_et = flux_et[:,1:]
    #print('flux_et')
    #print(flux_et)

    # get index where we have flux ET measurements;
    # only perform the objective calculation at these locations
    flux_data_idx = np.where(flux_et != -9999)

    # get shaes for later and initialize
    num_times, num_pixels = np.shape(flux_et)
    model_et = np.zeros(num_times)
    for ens in range(num_ens):
        model_et_all = np.zeros((num_times,num_pixels))
        time = 0
        root = exp_dir+'/output/SMAP_EASEv2_M36/cat/ens'+f"{ens:04}"
        for filename in glob.iglob(root+'/**/*.nc4',recursive=True):
            if (('lnd' in filename) and ('monthly' not in filename)):
                # get the model data necesarry
                data = Dataset(filename,mode='r')
                model_et_date = np.array(data['EVPTRNS'][:]) #will already be in ascending pixel order by default                
                model_et_date = model_et_date
                #print('np.shape(model_et_date)')
                #print(np.shape(model_et_date))
                #print('model_et_date')
                #print(model_et_date)
                #print(filename)
                model_et_all[time,:] = model_et_date
                time += 1
                #with open(exp_dir+'/run/output.txt','a') as f:
                #    f.write('str(model_et_date)')
                #    f.write('\n')
                #    f.write(str(model_et_date))
                #    f.write('\n')
                #    f.write('model_et_all')
                #    f.write('\n')
                #    f.write(str(model_et_all))
                #    f.write('\n')
                #    f.write('filename')
                #    f.write('\n')
                #    f.write(filename)
                #    f.write('\n')
        #get the flux data necesarry
        #print(np.shape(flux_et))
        #print(np.shape(model_et_all))
        #print(flux_et)
        #print(model_et_all)
        #print(model_et_all[0,:])
        curr_obj_val = np.sqrt(((np.sum(flux_et[flux_data_idx] - model_et_all[flux_data_idx]))**2)/(np.shape(flux_data_idx)[1]))
        #print(curr_obj_val)
        #assign to array to be returned
        objective_output[ens] = curr_obj_val
        
        with open(exp_dir+'/run/output.txt','a') as f:
            f.write('model_et_all')
            f.write('\n')
            f.write(str(model_et_all))
            f.write('\n')
        #print('curr_obj_val')
        #print(curr_obj_val)

    #with open(exp_dir+'/run/output.txt','a') as f:
    #    f.write('np.shape(flux_et)')
    #    f.write('\n')
    #    f.write(str(np.shape(flux_et)))
    #    f.write('\n')
    #    f.write('np.shape(model_et_all)')
    #    f.write('\n')
    #    f.write(str(np.shape(model_et_all)))
    #    f.write('\n')
    #    f.write('str(flux_et)')
    #    f.write(str(flux_et))
    #    f.write('str(model_et_all)')
    #    f.write(str(model_et_all))

    #print('objecitve_output')
    #print(objective_output)
    return objective_output
