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

def test_t2alm(seed=23333):
    np.random.seed(seed)
    T = np.asfortranarray(np.random.rand(npix, nsim))

    sht = SHT.SHT(nside, lmax, nsim, niter)

    alms = sht.t2alm_old(T)

    alms1_hp = sht.convert_alm_healpy(alms)
    alms1_hp = (alms1_hp[0,:,:] + 1j * alms1_hp[1,:,:])
    
    alms2_hp = np.array([hp.map2alm(T[:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    
    cl1 = np.array([hp.alm2cl(alms1_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i,:]) for i in range(nsim)])
    
    max_errTT = (np.abs(cl2 - cl1) / cl2.mean()).max()

    maps1 = sht.alm2t(alms)
    maps2 = np.array([hp.alm2map(alms1_hp[:,i], nside, pol=False) for i in range(nsim)])
    max_errT = (np.abs(maps1-maps2.transpose())/np.std(maps2)).max()

    return max_errT, max_errTT


def test_qu2eb(seed=23333):
    np.random.seed(seed)
    T = np.asfortranarray(np.random.rand(npix, nsim))
    Q = np.asfortranarray(np.random.rand(npix, nsim))
    U = np.asfortranarray(np.random.rand(npix, nsim))
    
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)
    
    start = time.time()
    almEs, almBs = sht.qu2eb(Q, U)

    almEs_hp = sht.convert_alm_healpy(almEs)
    almEs_hp = (almEs_hp[0,:,:] + 1j * almEs_hp[1,:,:])
    almBs_hp = sht.convert_alm_healpy(almBs)
    almBs_hp = (almBs_hp[0,:,:] + 1j * almBs_hp[1,:,:])
    
    maps = np.asfortranarray( [T, Q, U ] )
    
    start = time.time()
    alms2_hp = np.array([hp.map2alm(maps[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    
    cl = np.array([hp.alm2cl(almEs_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i, 1, :]) for i in range(nsim)])
    max_errEE = (np.abs(cl2 - cl) / cl.mean()).max()
    
    cl = np.array([hp.alm2cl(almBs_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i, 2, :]) for i in range(nsim)])
    max_errBB = (np.abs(cl2 - cl) / cl.mean()).max()

    q1, u1 = sht.eb2qu(almEs, almBs)
    maps2 = np.array([hp.alm2map(alms2_hp[i,:,:], nside, pol=True) for i in range(nsim)])
    q2 = maps2[:,1,:].transpose()
    u2 = maps2[:,2,:].transpose()

    max_errQ = (np.abs(q2 - q1)/q1.std()).max()
    max_errU = (np.abs(u2 - u1)/u1.std()).max()

    return max_errQ, max_errU, max_errEE, max_errBB


nside_list = [32, 128, 1024]
nsim = 8
niter = 1
max_err = 0

nside_rec, lmax_rec, errT_rec, errTT_rec, errQ_rec, errU_rec, errEE_rec, errBB_rec = [], [], [], [], [], [], [], []
for i in range(len(nside_list)):
    nside = nside_list[i]
    npix = 12 * nside ** 2
    lmax_list = [1*nside-1, 1*nside, 1*nside+1, 2*nside-1, 2*nside, 2*nside+1, 
                 3*nside-1, 3*nside, 3*nside+1, 4*nside-1, 4*nside]
    mm = len(lmax_list)
    for j in range(mm):
        lmax = lmax_list[j]
        print(" ")
        print("Testing nside=%i, lmax=%i" %(nside, lmax))
        errT, errTT = test_t2alm()
        errQ, errU, errEE, errBB = test_qu2eb()
        print('Max relative map-T error for the alm2t test is: ' + str(errT ))
        print('Max relative map-Q error for the alm2t test is: ' + str(errQ ))
        print('Max relative map-U error for the alm2t test is: ' + str(errU ))
        print('Max relative cl-TT error for the t2alm test is: ' + str(errTT))
        print('Max relative cl-EE error for the t2alm test is: ' + str(errEE))
        print('Max relative cl-BB error for the t2alm test is: ' + str(errBB))
        nside_rec.append(nside)
        lmax_rec.append(lmax)
        errT_rec.append(errT)
        errTT_rec.append(errTT)
        errQ_rec.append(errQ)
        errU_rec.append(errU)
        errEE_rec.append(errEE)
        errBB_rec.append(errBB)

print("%16s %16s %16s %16s %16s %16s %16s %16s" %("Nside", "Lmax", "Err_T", "Err_Q", "Err_U", 'Err_TT', "Err_EE", "Err_BB"))
for i in range(len(nside_rec)):
    print("%16i %16i %16.5e %16.5e %16.5e %16.5e %16.5e %16.5e" 
        %(nside_rec[i], lmax_rec[i], errT_rec[i], errQ_rec[i], errU_rec[i], errTT_rec[i], errEE_rec[i], errBB_rec[i]))

buff = np.array([errT_rec, errQ_rec, errU_rec, errTT_rec, errEE_rec, errBB_rec])
print("*********************************************************")
print("Summary: the maximum relative error in all tests is: %16.6e" %(np.amax(buff)))
print("*********************************************************")

        


