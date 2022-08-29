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

    alms = sht.t2alm(T)

    alms_hp1 = sht.convert_alm_healpy(alms)
    
    alms_hp2 = alms_hp1*0
    for i in range(nsim):
        alms_hp2[:,i] = hp.map2alm(T[:,i], lmax=lmax, iter=niter)
    
    cl1 = np.array([hp.alm2cl(alms_hp1[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms_hp2[:,i]) for i in range(nsim)])
    
    max_errTT = (np.abs(cl2 - cl1) / cl2.mean()).max()

    maps1 = sht.alm2t(alms)
    maps2 = maps1*0
    for i in range(nsim):
        maps2[:,i] = hp.alm2map(alms_hp1[:,i], nside, pol=False)
    max_errT = (np.abs(maps1-maps2)/np.std(maps2)).max()

    return max_errT, max_errTT


def test_qu2eb(seed=23333):
    np.random.seed(seed)
    T = np.asfortranarray(np.random.rand(npix, nsim))
    Q = np.asfortranarray(np.random.rand(npix, nsim))
    U = np.asfortranarray(np.random.rand(npix, nsim))
    
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)
    
    almEs_0, almBs_0 = sht.qu2e(Q, U), sht.qu2b(Q, U)
    almEs_1, almBs_1 = sht.qu2eb(Q, U)

    almEs_hp0 = sht.convert_alm_healpy(almEs_0)
    almBs_hp0 = sht.convert_alm_healpy(almBs_0)
    almEs_hp1 = sht.convert_alm_healpy(almEs_1)
    almBs_hp1 = sht.convert_alm_healpy(almBs_1)
    
    maps = np.asfortranarray( [T, Q, U] )
    
    # note: qu2e and qu2b has no iteration
    alms_hp2 = np.array([hp.map2alm(maps[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    alms_hp3 = np.array([hp.map2alm(maps[:,:,i], lmax=lmax, iter=0) for i in range(nsim)])
    almEs_hp2 = alms_hp2[:,1,:].transpose()
    almBs_hp2 = alms_hp2[:,2,:].transpose()
    almEs_hp3 = alms_hp3[:,1,:].transpose()
    almBs_hp3 = alms_hp3[:,2,:].transpose()
    
    cl0 = np.array([hp.alm2cl(almEs_hp0[:,i]) for i in range(nsim)])
    cl1 = np.array([hp.alm2cl(almEs_hp1[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(almEs_hp2[:,i]) for i in range(nsim)])
    cl3 = np.array([hp.alm2cl(almEs_hp3[:,i]) for i in range(nsim)])
    err0 = (np.abs(cl3 - cl0) / cl2.mean()).max()
    err1 = (np.abs(cl2 - cl1) / cl2.mean()).max()
    max_errEE = max(err0, err1)
    
    
    cl0 = np.array([hp.alm2cl(almBs_hp0[:,i]) for i in range(nsim)])
    cl1 = np.array([hp.alm2cl(almBs_hp1[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(almBs_hp2[:,i]) for i in range(nsim)])
    cl3 = np.array([hp.alm2cl(almBs_hp3[:,i]) for i in range(nsim)])
    err0 = (np.abs(cl3 - cl0) / cl2.mean()).max()
    err1 = (np.abs(cl2 - cl1) / cl2.mean()).max()
    max_errBB = max(err0, err1)

    Q1, U1 = sht.eb2qu(almEs_1, almBs_1)
    maps2 = np.array([hp.alm2map(alms_hp2[i,:,:], nside, pol=True) for i in range(nsim)])
    Q2 = maps2[:,1,:].transpose()
    U2 = maps2[:,2,:].transpose()

    max_errQ = (np.abs(Q2 - Q1)/Q1.std()).max()
    max_errU = (np.abs(U2 - U1)/U1.std()).max()

    almEs_tmp = sht.qu2e(Q, U)
    almBs_tmp = sht.qu2b(Q, U)

    return max_errQ, max_errU, max_errEE, max_errBB


nside_list = [32, 128, 1024]
nsim = 4
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

        


