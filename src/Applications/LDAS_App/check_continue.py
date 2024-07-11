import numpy as np
import sys
import os

class check:
    def __init__(self,exp_dir):
        self.exp_dir = exp_dir
        #print('initialized')
    def check_run_done(self,num_mems):
        #print('inside check done')
        num_mems_idx = num_mems - 1
        total_ens = 0
        for i in range(num_mems):
            this_fname = os.path.join(self.exp_dir,'../','finished_'+str(i)+'.txt')
            if os.path.isfile(this_fname) == True:
                total_ens += 1
            with open(os.path.join(self.exp_dir,'run/','check_runs.txt'),'a') as f:
                f.write('this_fname')
                f.write('\n')
                f.write(str(this_fname))
                f.write('\n')
                f.write('total_ens')
                f.write('\n')
                f.write(str(total_ens))
                f.write('\n')
        if (total_ens == num_mems):
            with open(os.path.join(self.exp_dir,'run/','check_runs.txt'),'a') as f:
                f.write('got here in convergence')
                f.write('\n')
            print(1)
        else:
            print(0)
    def check_pso_done(self):
        #print('in pso done')
        pso_finished_fname = os.path.join(self.exp_dir,'../','pso_finished.txt')
        #print(pso_finished_fname)
        if os.path.isfile(pso_finished_fname) == True:
            pso = 1
        else:
            pso = 0
        with open(os.path.join(self.exp_dir,'run','check_pso.txt'),'a') as f:
            f.write('pso_finished_fname')
            f.write('\n')
            f.write(str(pso_finished_fname))
            f.write('\n')
            f.write('pso')
            f.write('\n')
            f.write(str(pso))
            f.write('\n')
        if pso == 1:
            print(1)
        else:
            print(0)

def main(check_code,exp_dir,num_mems):
    #print('in main')
    ch = check(exp_dir)
    #print(check_code)
    #print(type(check_code))
    if (check_code == 0):
        #print('about to run check')
        ch.check_run_done(num_mems)
        #print('ran check')
    elif (check_code == 1):
        ch.check_pso_done()

if __name__ == '__main__':
     check_code = int(sys.argv[1])
     exp_dir  = sys.argv[2]
     num_mems = 30
     #print(check_code)
     #print(exp_dir)
     main(check_code,exp_dir,num_mems)
