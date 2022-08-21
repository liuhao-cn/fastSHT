#!/usr/bin/env python
# coding: utf-8

# In[1]:
import sys as sys
import os

# default parameters
nside = 64
nsim = 2000
n_proc = 8
niter = 3
compare = False

# the command line input will overwrite the defaults
if len(sys.argv)>1:
    nside = int(sys.argv[1])
if len(sys.argv)>2:
    nsim = int(sys.argv[2])
if len(sys.argv)>3:
    n_proc = int(sys.argv[3])
if len(sys.argv)>4:
    niter = int(sys.argv[4])
if len(sys.argv)>5:
    compare = sys.argv[5].lower() == "true"

print(" ")
print("Working with the following parameters:")
print("Nside = %i, Nsim = %i, n_proc = %i, Niter = %i, comparison with HEALPix = %s" 
    %(nside, nsim, n_proc, niter, compare))
print("Note: this test will only show the time cost. For accuracies use test_comprehensive.py")
print(" ")

os.environ["OMP_NUM_THREADS"] = str(n_proc) # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = str(n_proc) # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = str(n_proc) # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = str(n_proc) # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = str(n_proc) # export NUMEXPR_NUM_THREADS=6

sys.path.append('..')
import SHT
import fastSHT


import numpy as np

import healpy as hp
import matplotlib.pyplot as plt
import time

import numba
from numba import cuda

import importlib
importlib.reload(SHT)


nrep = 3
lmax = 3*nside - 1
npix = 12 * nside ** 2

def test_t2alm(seed=23333):
    print('Testing t2alm...')

    np.random.seed(seed)
    maps = np.asfortranarray(np.random.rand(npix, nsim))

    start = time.time()
    sht = SHT.SHT(nside, lmax, nsim, niter)
    #alms = numba.cuda.pinned_array((nsim, lmax+1, lmax+1), dtype=np.double, strides=None, order='F')
    alms = np.empty((nsim, lmax+1, lmax+1), dtype=np.double,  order='F')
    end = time.time() - start
    print('Time cost for memory initialization is ' + str(end))

    start = time.time()
    for i in range(nrep):
        sht.t2alm_old(maps, alms)
    end = time.time() - start
    print('Calculation time cost for fastSHT is ' + str(end / nrep))

    #print(time.sleep(10))
    if(compare==False):
        return
    
    alms_hp = sht.convert_alm_healpy(alms); del alms
    alms_hp = (alms_hp[0,:,:] + 1j * alms_hp[1,:,:])
    alm_shape = alms_hp.shape[0]
    del alms_hp

    start = time.time()
    #alms2_hp = np.array([hp.map2alm(maps[:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    overhead = 0
    roverhead = 0
    for j in range(nrep):
        for i in range(nsim):
            ostart = time.time()
            maps[:,i].astype(np.float64, order='C', copy=True)
            rstart = time.time()
            np.empty((alm_shape),order='C', dtype=np.complex128)
            end = time.time()
            overhead += end - ostart
            roverhead += end - rstart
            hp.map2alm(maps[:,i], lmax=lmax, iter=niter)

    tot_time = time.time() - start
    print('Time cost for memory initialization is ' + str(roverhead/nrep))
    print('Calculation time cost for healpy is ' + str((tot_time - overhead)/nrep ))


def make_mask(nside, upper_lat = 5):
    npix = nside**2 * 12
    mask = np.ones(npix)
    
    for i in range(npix):
        theta, phi = hp.pix2ang(nside, i, lonlat=True)
        if phi < 5 and phi > -5:
            mask[i] = 0
    return mask

# In[7]:

def fix_EB(Q, U, mask, nside, lmax, niter=0, seed=23333):
    
    T = np.asfortranarray(np.random.rand(npix, nsim))

    vid = (np.arange(len(mask))[mask == 1])
    nv = len(vid)

    maps_in = np.array( [T, Q, U ])
    maps_in[:,mask==0,:] = hp.UNSEEN
    
    start = time.time()
    
    alms2_hp = np.array([hp.map2alm(maps_in[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    
    BO = np.array([hp.alm2map(alms2_hp[i,2,:], nside=nside, lmax=lmax) for i in range(nsim)])
    
    alms2_hp[:,2,:] = 0
    alms2_hp[:,0,:] = 0
    
    maps = np.array([hp.alm2map(alms2_hp[i,:,:], nside=nside, lmax=lmax) for i in range(nsim)])
    maps[:,:,mask==0] = 0
    alms2_hp = np.array([hp.map2alm(maps[i,:,:], lmax=lmax, iter=niter) for i in range(nsim)])
    alms2_hp = alms2_hp[:,2,:]
    BT = np.array([hp.alm2map(alms2_hp[i,:], nside=nside, lmax=lmax) for i in range(nsim)])
    BC = np.zeros(BT.shape)
    for i in range(nsim):
        x = BT[i, vid]
        y = BO[i, vid]
        coe = np.polyfit(x, y, 1)
        
        # mx = np.sum(x) / nv
        # my = np.sum(y) / nv
        # cxx = np.sum( (x - mx)**2 )
        # cxy = np.sum( (y - my) * (x - mx) )
        # coe = [cxy / cxx, my  - mx * (cxy / cxx) ]
        
        BC[i, vid] = BO[i, vid] - BT[i, vid] * coe[0] - coe[1]
        
    # print('Time cost for Healpy is ' + str(time.time() - start))
    return BC

def test_fix_EB(seed=23333):
    npix = 12*nside**2
    print('Testing fix_EB...')
    np.random.seed(seed)
    Q = np.asfortranarray(np.random.rand(npix, nsim))
    U = np.asfortranarray(np.random.rand(npix, nsim))
    
    mask = make_mask(nside)

    start = time.time()
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)

    Bmap = sht.fix_eb(Q, U, mask)
    print('Time cost for fastSHT is ' + str(time.time() - start))
    if(compare == False):
        return
    
    start = time.time()
    BC = fix_EB(Q, U, mask, nside, lmax, niter, seed)
    print('Time cost for healpy is ' + str(time.time() - start))

    return (Bmap, BC)

def test_qu2eb(seed=23333):

    print('Testing qu2eb...')
    np.random.seed(seed)
    
    T = np.asfortranarray(np.random.rand(npix, nsim))
    Q = np.asfortranarray(np.random.rand(npix, nsim))
    U = np.asfortranarray(np.random.rand(npix, nsim))
    
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)
    
    start = time.time()
    almEs, almBs = sht.qu2eb(Q, U)
    print('Time cost for fastSHT is ' + str(time.time() - start))

    if(compare == False):
        return

    almEs_hp = sht.convert_alm_healpy(almEs)
    almEs_hp = (almEs_hp[0,:,:] + 1j * almEs_hp[1,:,:])
    almBs_hp = sht.convert_alm_healpy(almBs)
    almBs_hp = (almBs_hp[0,:,:] + 1j * almBs_hp[1,:,:])
    
    maps = np.asfortranarray( [T, Q, U ])
    
    start = time.time()
    alms2_hp = np.array([hp.map2alm(maps[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    print('Time cost for healpy is ' + str(time.time() - start))
    

# test_t2alm()
# test_qu2eb()
test_fix_EB()
