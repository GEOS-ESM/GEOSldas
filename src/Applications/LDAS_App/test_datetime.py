import pandas as pd
import numpy as np
from datetime import datetime

flux_et = pd.read_csv(
    '/discover/nobackup/trobinet/misc_and_testing/create_env_data/le_pixels_weighted_2005-01-01_2005-12-31.csv'
)
#print(flux_et)
times = np.array(flux_et['time'])
flux_et = np.array(flux_et)
print(times)

grow_months = [5,6,7,8,9,10]
time_idx = np.zeros(0,dtype=int)
for i,time in enumerate(times):
    time_datetime = datetime.strptime(time, '%Y-%m-%d')
    if time_datetime.month in grow_months:
        print(i)
        print(time_datetime.month)
        time_idx = np.append(time_idx,i)
print(time_idx)
flux_et_time = flux_et[time_idx,:]
print(flux_et_time)
