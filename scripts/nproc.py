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
compare = True


# the command line input will overwrite the defaults
if len(sys.argv)>1:
    nside = int(sys.argv[1])
if len(sys.argv)>2:
    nsim = int(sys.argv[2])
if len(sys.argv)>3:
    n_proc = int(sys.argv[3])
if len(sys.argv)>4:
    niter = int(sys.argv[4])


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


nrep = 1
lmax = 3*nside - 1
npix = 12 * nside ** 2

def test_t2alm(seed=23333):
    # print('Testing t2alm...')

    np.random.seed(seed)
    maps = np.asfortranarray(np.random.rand(npix, nsim))

    start = time.time()
    sht = SHT.SHT(nside, lmax, nsim, niter)
    #alms = numba.cuda.pinned_array((nsim, lmax+1, lmax+1), dtype=np.double, strides=None, order='F')
    alms = np.empty((nsim, lmax+1, lmax+1), dtype=np.double,  order='F')
    end = time.time() - start
    # print('Time cost for memory initialization is ' + str(end))

    start = time.time()
    sht.t2alm(maps, alms)
    end1 = time.time() - start
    # print('Calculation time cost for fastSHT is ' + str(end / nrep))

    #print(time.sleep(10))
    if(compare==False):
        return
    
    alms_hp = sht.convert_alm_healpy(alms); del alms
    alm_shape = alms_hp.shape[0]
    del alms_hp

    start = time.time()
    for i in range(nsim):
        hp.map2alm(maps[:,i], lmax=lmax, iter=niter)
    end2 = time.time() - start
    
    return end1, end2

end1, end2 = test_t2alm()

print("t2alm: Nside = %i, Nsim = %i, n_proc = %2i, fastSHT-CPU = %6.2f, healpy = %6.2f" 
    %(nside, nsim, n_proc, end1, end2))