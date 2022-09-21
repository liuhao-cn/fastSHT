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
test_type = 't2alm'
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
    test_type = sys.argv[5].lower()
if len(sys.argv)>6:
    compare = sys.argv[6].lower() == "true"

if test_type=='t2alm':
    tstr = 't2alm'
elif test_type=='qu2eb':
    tstr = 'qu2eb'
else:
    tstr = 'fix-eb'


print(" ")
print("Working with the following parameters:")
print("Nside = %i, Nsim = %i, n_proc = %i, Niter = %i, comparison with HEALPix = %s, test type=%s" 
    %(nside, nsim, n_proc, niter, compare, tstr))
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
        sht.t2alm(maps, alms)
    end = time.time() - start
    print('Calculation time cost for fastSHT is ' + str(end / nrep))

    #print(time.sleep(10))
    if(compare==False):
        return
    
    alms_hp = sht.convert_alm_healpy(alms); del alms
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
    T = np.zeros((npix,))

    vid = (np.arange(len(mask))[mask == 1])
    BC = np.zeros((nsim,npix))

    start = time.time()
    for i in range(nsim):
        maps_in = np.array( [T, Q[:,i], U[:,i]])
        maps_in[:,mask==0] = hp.UNSEEN

        alms2_hp = hp.map2alm(maps_in, lmax=lmax, iter=niter)
        BO = hp.alm2map(alms2_hp[2,:], nside=nside, lmax=lmax)

        alms2_hp[2,:] = 0
        alms2_hp[0,:] = 0
        maps = hp.alm2map(alms2_hp, nside=nside, lmax=lmax)
        maps[:,mask==0] = 0
        alms2_hp = hp.map2alm(maps, lmax=lmax, iter=niter)
        BT = hp.alm2map(alms2_hp[2,:], nside=nside, lmax=lmax)

        x = BT[vid]
        y = BO[vid]
        coe = np.polyfit(x, y, 1)

        BC[i, vid] = y - x * coe[0] - coe[1]
    print('Time cost for Healpy is ' + str(time.time() - start))

    return BC

def test_fix_EB(seed=23333):
    npix = 12*nside**2
    print('Testing fix_EB...')
    np.random.seed(seed)
    Q = np.asfortranarray(np.random.rand(npix, nsim))
    U = np.asfortranarray(np.random.rand(npix, nsim))
    
    mask = make_mask(nside)

    start = time.time()
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True, all_buff=True)

    Bmap = sht.fix_eb(Q, U, mask)
    print('Time cost for fastSHT is ' + str(time.time() - start))
    if(compare == False):
        return
    
    BC = fix_EB(Q, U, mask, nside, lmax, niter, seed)

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
    almBs_hp = sht.convert_alm_healpy(almBs)
    
    maps = np.asfortranarray( [T, Q, U ])
    
    start = time.time()
    alms2_hp = np.array([hp.map2alm(maps[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    print('Time cost for healpy is ' + str(time.time() - start))
    

if test_type=='t2alm':
    test_t2alm()
elif test_type=='qu2eb':
    test_qu2eb()
else:
    test_fix_EB()
