#!/usr/bin/env python
#
# module load python/GEOSpyD/Ana2019.03_py3.7
# module load nco/4.8.1
#
# Script for creating ensemble-perturbed  land forcing (lfo) files.
# Usage:  ensemble_forc.py  [in_path]  [cntr_path]  [nens]
#
# where 
# 
#   in_path  : path to ensemble lfo files
#   cntr_path : path to deterministic lfo files
#   nens     : number of ensemble members
#===========================================

import sys
import os
import glob
import subprocess as sp
import numpy as np
from netCDF4 import Dataset


def averaging_forcing(in_path, avg_path, nens):
   """ The ensemble number will be appended to in_path starting from 001.
       The out_path will be created if it does not exist. """ 
   print ( " average ensemble " ) 
   if not os.path.exists(avg_path):
      os.makedirs(avg_path) 
   files_list=[]
   for i in range(1,nens+1):
      sfx = '%03d'%(i)
      folder  = in_path+'/mem'+sfx
      fs      = sorted(glob.glob(folder+'/*lfo*.nc4'))
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
      cmd = 'ncea ' + fnames + ' ' + avg_path +'/' + f0         
      sp.call(cmd, shell=True)

def mean_diff(ensm_path, cntr_path): 
   print ( " prepare central - ensmean for recentering ") 
   files_list=[]
   folder  = ensm_path
   fs      = sorted(glob.glob(folder+'/*lfo*.nc4'))
   files_list.append(fs)
   for fs in zip(*files_list):
      k = 1
      for f in fs:
         n = f.rindex('/')
         if (k==1):
            f0 = f[n+1:]
            k  = k+1
         f1 = f[n+1:]
         assert f0 == f1, "central and ensmean should have same filenames"
         cmd = 'ncbo --op_typ=sbt ' + cntr_path + '/' + f1 + ' ' +  f  + ' ' + f +'.dif'
         #print ( cmd )
         sp.call(cmd, shell=True)

def recenter_forc(in_path, avg_path, cntr_path, nens):
   print ( "recenter ens forcing to central forcing " ) 
   Varlist = ["PRECCU", "PRECLS", "PRECSNO", "SWGDN" ]
   print ( "Varlist for multipl-recentering: ", Varlist) 

   mean_diff(avg_path, cntr_path)
   
   files_list1=[]
   for i in range(1,nens+1):
      sfx = '%03d'%(i)
      folder  = in_path+'/mem'+sfx
      fs1      = sorted(glob.glob(folder+'/*lfo*.nc4')) 
      files_list1.append(fs1)
   for fs in zip(*files_list1):
      k = 1
      for f in fs:
         n = f.rindex('/')
         if (k==1):
            f0 = f[n+1:]
            k  = k+1
         f1 = f[n+1:]
         assert f0 == f1, "memeber and ensmean should have same filenames"
         cmd = 'ncbo --op_typ=add ' + f + ' ' + avg_path +'/' + f1 + '.dif' + ' ' + f +'.Cpert'
         sp.call(cmd, shell=True) 

   files_list=[]
   for i in range(1,nens+1):
      sfx = '%03d'%(i)
      folder  = in_path+'/mem'+sfx
      fs      = sorted(glob.glob(folder+'/*tavg*lfo*.nc4'))
      files_list.append(fs)
   for fs in zip(*files_list):
      k = 1
      for f in fs:
         n = f.rindex('/')
         if (k==1):
            f0 = f[n+1:]
            m_file = avg_path+'/'+f0 
            Var_m = readlfo(m_file, Varlist)
            c_file = cntr_path+'/'+f0
            Var_c = readlfo(c_file, Varlist)
            k  = k+1
         f1 = f[n+1:]
         fout = f + '.Cpert' 
         Var_e = readlfo(f,  Varlist)
         Var_ux = fact_multpl(Var_m,Var_e,Var_c)
         upd_vars(fout, Var_ux, Varlist)

def readlfo(ncfile, Varlist):
    wkfile = Dataset(ncfile, "r")
    nydim = len(wkfile.dimensions["Ydim"])
    nxdim = len(wkfile.dimensions["Xdim"])
    nfdim = len(wkfile.dimensions["nf"])
    dim3 = [nfdim, nydim, nxdim]
    nvar = len(Varlist)
    vardata = np.zeros((nvar, dim3[0], dim3[1], dim3[2]), dtype=float)
    for rdvar in iter(range(nvar)):
            vardata[rdvar] = wkfile.variables[Varlist[rdvar]][0]

    wkfile.close()
    return vardata

def upd_vars(ncfile, Dataupd, Varupd):

    wkfile = Dataset(ncfile,'r+')
    nvar = len(Varupd)
    for rdvar in iter(range(nvar)) :
        wkfile.variables[Varupd[rdvar]][0] = Dataupd[rdvar]

    wkfile.close()

def fact_multpl(Var_dn,Var_up,Var_cn): 

    dim = np.shape(Var_dn) 
    ndim = len(dim)  
    dim3 = [ dim[1], dim[2],dim[3] ] 
    ndim3 = len(dim3)

    #sum precip [0, 1, 2 from Varlist(4)  ) 
    Var_3pup = np.zeros(dim3, dtype=float)
    Var_3pdn = np.zeros(dim3, dtype=float)
    Var_3pup[:,:,:] = Var_up[0,:,:,:] + Var_up[1,:,:,:] + Var_up[2,:,:,:]
    Var_3pdn[:,:,:] = Var_dn[0,:,:,:] + Var_dn[1,:,:,:] + Var_dn[2,:,:,:]
    nvec = 1 
    for id in iter(range(ndim3)):
        nvec = dim3[id]*nvec 

    vecdn = np.reshape(Var_3pdn,nvec)
    vecup = np.reshape(Var_3pup,nvec)
    vecnew = np.zeros(nvec,dtype=float)
    
    ind = np.argwhere(vecdn>1.0e-10)
    vecnew [ind] =  vecup[ind] / vecdn[ind]
     
    Var_uu = np.zeros((dim3), dtype=float)
    Var_uu = np.reshape(vecnew,(dim3)) 
    # apply the same factor to 3 precip 
    Var_upx = np.zeros((dim), dtype=float)
    Var_upx[0,:,:,:] = Var_uu[:,:,:] * Var_cn[0,:,:,:]
    Var_upx[1,:,:,:] = Var_uu[:,:,:] * Var_cn[1,:,:,:]
    Var_upx[2,:,:,:] = Var_uu[:,:,:] * Var_cn[2,:,:,:]
    # deal with swgdn 
    Var_3pup[:,:,:] = Var_up[3,:,:,:]
    Var_3pdn[:,:,:] = Var_dn[3,:,:,:]
    vecdn = np.reshape(Var_3pdn,nvec)
    vecup = np.reshape(Var_3pup,nvec)
    vecnew = np.zeros(nvec,dtype=float)
    ind = np.argwhere(vecdn!=0)
    vecnew [ind] =  vecup[ind] / vecdn[ind]
    Var_uu = np.zeros((dim3), dtype=float)
    Var_uu = np.reshape(vecnew,(dim3))
    Var_upx[3,:,:,:] = Var_uu[:,:,:] * Var_cn[3,:,:,:] 
    
    return  Var_upx

   
if __name__ == '__main__' :
 
   in_path  =     sys.argv[1]
   cntr_path =    sys.argv[2]
   nens     = int(sys.argv[3])
   mean_path = in_path + '/lfomean' 
   averaging_forcing(in_path,mean_path,  nens) 
   recenter_forc(in_path, mean_path, cntr_path, nens)
