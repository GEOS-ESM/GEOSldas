import numpy as np
def fake_mod(num_times,num_pixels,positions,ens):
    curr_pred = np.zeros((num_times,num_pixels))
    curr_pred[0,:] = positions[ens,0]**2
    curr_pred[1,:] = positions[ens,1]/5
    correct_ans = np.zeros((num_times,num_pixels))
    correct_ans[0,:] = 9
    correct_ans[1,:] = 2
    return (curr_pred,correct_ans)
