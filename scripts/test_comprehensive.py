#!/usr/bin/env python
# coding: utf-8

import sys as sys
import os

n_proc = 8
os.environ["OMP_NUM_THREADS"] = str(n_proc) # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = str(n_proc) # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = str(n_proc) # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = str(n_proc) # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = str(n_proc) # export NUMEXPR_NUM_THREADS=6

sys.path.append('..')
import SHT, fastSHT, time

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

import importlib
importlib.reload(SHT)

def test_t2alm(nside, lmax, nsim, niter=0, seed=23333, compare=True):
    np.random.seed(seed)
    maps = np.asfortranarray(np.random.rand(npix, nsim))

    sht = SHT.SHT(nside, lmax, nsim, niter)

    alms = sht.t2alm_old(maps)

    if(compare == False):
        return
    
    alms_hp = sht.convert_alm_healpy(alms)
    alms_hp = (alms_hp[0,:,:] + 1j * alms_hp[1,:,:])
    
    start = time.time()
    alms2_hp = np.array([hp.map2alm(maps[:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    
    cl = np.array([hp.alm2cl(alms_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i,:]) for i in range(nsim)])
    
    max_err = (np.abs(cl2 - cl) / cl2.mean()).max()

    print('Max relative cl-TT error for the t2alm test is: ' + str(max_err))
    return max_err

# In[5]:


def test_qu2eb(nside, lmax, nsim, niter=0, seed=23333, compare=True):
    np.random.seed(seed)
    T = np.asfortranarray(np.random.rand(npix, nsim))
    Q = np.asfortranarray(np.random.rand(npix, nsim))
    U = np.asfortranarray(np.random.rand(npix, nsim))
    
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)
    
    start = time.time()
    almEs, almBs = sht.qu2eb(Q, U)

    if(compare == False):
        return

    almEs_hp = sht.convert_alm_healpy(almEs)
    almEs_hp = (almEs_hp[0,:,:] + 1j * almEs_hp[1,:,:])
    almBs_hp = sht.convert_alm_healpy(almBs)
    almBs_hp = (almBs_hp[0,:,:] + 1j * almBs_hp[1,:,:])
    
    maps = np.asfortranarray( [T, Q, U ] )
    
    start = time.time()
    alms2_hp = np.array([hp.map2alm(maps[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    
    cl = np.array([hp.alm2cl(almEs_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i, 1, :]) for i in range(nsim)])
    
    max_errE = (np.abs(cl2 - cl) / cl.mean()).max()
    print('Max relative cl-EE error in the qu2eb test is: ' + str(max_errE))
    
    cl = np.array([hp.alm2cl(almBs_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i, 2, :]) for i in range(nsim)])
    
    max_errB = (np.abs(cl2 - cl) / cl.mean()).max()
    print('Max relative cl-BB error in the qu2eb test is: ' + str(max_errB))
    return max_errE, max_errB


nside_list = [32, 128, 1024]
nsim = 8
nn = len(nside_list)
max_err = 0

nside_rec, lmax_rec, errT_rec, errE_rec, errB_rec = [], [], [], [], []
for i in range(nn):
    nside = nside_list[i]
    npix = 12 * nside ** 2
    lmax_list = [1*nside-1, 1*nside, 1*nside+1, 2*nside-1, 2*nside, 2*nside+1, 
                 3*nside-1, 3*nside, 3*nside+1, 4*nside-1, 4*nside]
    mm = len(lmax_list)
    for j in range(mm):
        lmax = lmax_list[j]
        print(" ")
        print("Testing nside=%i, lmax=%i" %(nside, lmax))
        errT = test_t2alm(nside, lmax, nsim, niter=1, compare=True)
        errE, errB = test_qu2eb(nside, lmax, nsim, niter=1, compare=True)
        nside_rec.append(nside)
        lmax_rec.append(lmax)
        errT_rec.append(errT)
        errE_rec.append(errE)
        errB_rec.append(errB)

print("%16s %16s %16s %16s %16s" %("Nside", "Lmax", "Err_T", "Err_E", "Err_B"))
for i in range(len(nside_rec)):
    print("%16i %16i %16.5e %16.5e %16.5e" %(nside_rec[i], lmax_rec[i], errT_rec[i], errE_rec[i], errB_rec[i]))

buff = np.array([errT_rec, errE_rec, errB_rec])
print("*********************************************************")
print("Summary: the maximum relative error in all tests is: %16.6e" %(np.amax(buff)))
print("*********************************************************")

        


