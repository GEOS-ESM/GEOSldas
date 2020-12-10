#!/usr/bin/env python
#
# module load python/GEOSpyD/Ana2019.03_py3.7
# module  load nco/4.8.1
#
import os
import glob
import subprocess as sp

def averaging_forcing(in_path, out_path, nens):
   """ The ensemble number will be appended to in_path
       it starts from 001
       out_path will be created if it does not exists """ 
   if not os.path.exists(out_path):
      os.makedirs(out_path) 
   files_list=[]
   for i in range(1,nens+1):
      sfx = f'{i:03}'
      folder  = in_path+sfx
      fs      = sorted(glob.glob(folder+'/*.nc4'))
      files_list.append(fs)
   for fs in zip(*files_list):
      k = 1
      # verify filenames are the same.
      for f in fs:
         n = f.rindex('/')
         if (k==1):
            f0 = f[n+1:]
            k  = k+1
         f1 = f[n+1:]
         assert f0 == f1, "averaging different files. Each folder should have same files"
      fnames = " ".join(fs)
      cmd = 'ncea ' + fnames + ' ' + out_path +'/' + f0         
      sp.call(cmd, shell=True)
   
if __name__ == '__main__' :
   # These there arguments should be changed accordingly
   in_path  = '/discover/nobackup/qzhang/forcing/stage/ensdiag/mem'
   nens     = 32
   out_path = './mem_avg'

   averaging_forcing(in_path, out_path, nens)
