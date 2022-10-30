import numpy as np
import pandas as pd
import glob
from netCDF4 import Dataset
import sys

def compute_objective(num_parts,num_params,positions):
    objective_output = np.zeros(num_parts)
    truth_vals = np.array((1,2,3,-9999,5))
    for i in range(num_parts):
        # get index where we have flux ET measurements;
        # only perform the objective calculation at these locations
        active_pos = np.where(truth_vals != -9999)
        active_pos=0
        this_pos = positions[i,:]
        curr_obj_val = np.sqrt(((this_pos[active_pos] - truth_vals[active_pos]) ** 2).mean())
        #assign to array to be returned
        objective_output[i] = curr_obj_val
        
    return objective_output
